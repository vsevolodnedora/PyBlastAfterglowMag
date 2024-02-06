import os
import re
import h5py
import copy
import numpy as np
from scipy import ndimage, interpolate, spatial
from glob import glob

def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx

class ProcessRawSkymap():
    """
    Processes raw skymaps from PyBlastAfterglowas
    - Compute flux centroid positions from raw data as weighted averages of intensity
    - Combine raw skymaps from multiple shells by
        i) interpolating them onto the uniform grid and summing up (optional: applying Gaussian smoothing kernel after)
        ii) making a 2D histogram of intensity to mimic the photon collection of a CCD (how the observer would see it)  (optional: applying Gaussian smoothing kernel after)
    - Computing X and Y integrated distributions from interpolated skymap (integrating intensity along the other axis)
    - Computing image size at flux centroid positions as "full-width-half-maximum" using integrated distributions
    - saving the image properties for all times and frequencies in a single small-size file
    """

    def __init__(self, conf : dict, verbose=False):
        self.verb = verbose
        # Settings for the filter to be convolved with the skymap (smoothing) (if type==None, nothing done)
        self.conf = conf
        self.conf_intp_filter = self.conf["intp_filter"]#{ "type":"gaussian", "sigma":2, "mode":'reflect' }
        self.conf_hist_filter = self.conf["hist_filter"]#{ "type":"gaussian", "sigma":2, "mode":'reflect' }


    @staticmethod
    def _apply_filter(image : np.ndarray, conf : dict) -> np.ndarray:
        """ convolves the image with a given filter """
        type=conf["type"]
        if type is None:
            return image
        elif type=="uniform":
            size=conf["size"]
            image = ndimage.uniform_filter(image, size=size )
        elif type == "gaussian":
            sigma=conf["sigma"]
            mode = conf["mode"]
            if image.ndim == 1:
                image = ndimage.filters.gaussian_filter(image, sigma=[sigma], mode=mode)
            elif image.ndim == 2:
                image = ndimage.filters.gaussian_filter(image, sigma=[sigma,sigma], mode=mode)
        else:
            raise KeyError("apply_filter type is not recognize: {}".format(type))
        return image

    @staticmethod
    def _get_skymap_fwhm(grid, latAvDist, cc, fac=0.5):
        """ Assuming a bell-shaped function, find its X1 and X2 where f = fac*f_at_cc """
        assert len(grid) == len(latAvDist)
        val = latAvDist[find_nearest_index(grid, cc)]
        val =  np.interp(cc, grid, latAvDist)# latAvDist[find_nearest_index(grid, cc)]
        max_val = np.max(latAvDist)
        x1 = grid[np.argmin(latAvDist < val * fac)]
        # x1 = np.interp(val * fac, latAvDist, grid) #  grid[np.argmin(latAvDist < val * fac)]
        # plt.plot(grid, latAvDist)
        # plt.axhline(y=val * fac)
        # plt.show()
        x2 = grid[::-1][np.argmin(latAvDist[::-1] < val * fac)]

        if (x2 < x1):
            raise ValueError(f"x2={x2} < x1={x1}")

        return (x1, x2)

    @staticmethod
    def _interp(xxs : np.ndarray, yys : np.ndarray, fluxes : np.ndarray, x_grid : np.ndarray, y_grid : np.ndarray,
                method='linear'):
        """ using unstractured D-D interpolation """

        assert xxs.ndim == 1, f" input shapes should be 1. Given {xxs.ndim}"
        assert yys.ndim == 1, f" input shapes should be 1. Given {xxs.ndim}"
        assert fluxes.ndim == 1, f" input shapes should be 1. Given {xxs.ndim}"

        # assert len(xxs) > 3, f"Size of xxs={len(xxs)}"
        # assert len(yys) > 3, f"Size of yys={len(xxs)}"
        # assert len(fluxes) > 3, f"Size of fluxes={len(xxs)}"

        try:
            res = interpolate.griddata(np.vstack((xxs, yys)).T, fluxes,
                                   np.array((x_grid, y_grid)).T, method=method, fill_value=0.)
        except:
            raise ValueError("Failed interpolation")

        return res

    @staticmethod
    def _compute_position_of_the_flux_centroid(all_x_arrs, all_y_arrs, all_z_arrs):
        """ flux centroid calculation """
        xcs = np.average(np.concatenate(all_x_arrs), weights=np.concatenate(all_z_arrs))
        xcs_m = xcs
        ycs = np.average(np.concatenate(all_y_arrs), weights=np.concatenate(all_z_arrs))
        ycs_m = ycs
        return(xcs_m, ycs_m) # in milli-arc-seconds


    def _get_data_extend(self, xs, ys, datas):
        xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
        ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
        i_min, i_max = [], []
        i_shells = []
        ncells = []

        nshells = len(xs)

        for ish in range(nshells):
            xrs_i = np.array(xs[ish])
            yrs_i = np.array(ys[ish])
            int_i = np.array(datas[ish])

            # skip empty shells
            if ((np.sum(int_i) == 0)):
                continue

            # if (len(xrs_i) % 2 > 0):
            #     raise ValueError(f"expected to get an even number for ncells. Got:{len(xrs_i)}")

            ncells.append( int( len(xrs_i) / 2 ) )
            i_shells.append(ish)
            #
            # if (type=="principle"):
            #     xrs_i = xrs_i
            #     yrs_i = yrs_i[:ncells[-1]]
            #     int_i = int_i[:ncells[-1]]
            # elif(type=="counter"):
            #     xrs_i = xrs_i
            #     yrs_i = yrs_i[ncells[-1]:]
            #     int_i = int_i[ncells[-1]:]
            # elif(type=="both"):
            #     pass
            # else:
            #     raise KeyError(f"type={type} is not recognized")

            xmin_neg.append( xrs_i[xrs_i < 0].min() if len(xrs_i[xrs_i < 0]) > 0 else 0 )
            xmin_pos.append( xrs_i[xrs_i > 0].min() if len(xrs_i[xrs_i > 0]) > 0 else 0 )
            xmax.append( xrs_i.max() )
            xmin.append( xrs_i.min() )
            ymin_neg.append( yrs_i[yrs_i < 0].min() if len(yrs_i[yrs_i < 0]) > 0 else 0 )
            ymin_pos.append( yrs_i[yrs_i > 0].min() if len(yrs_i[yrs_i > 0]) > 0 else 0 )
            ymax.append( yrs_i.max() )
            ymin.append( yrs_i.min() )
            i_min.append( int_i.min() )
            i_max.append( int_i.max() )

        xmin_neg = np.array(xmin_neg)
        xmin_pos = np.array(xmin_pos)
        xmax = np.array(xmax)
        xmin = np.array(xmin)
        ymin_neg = np.array(ymin_neg)
        ymin_pos = np.array(ymin_pos)
        ymax = np.array(ymax)
        ymin = np.array(ymin)
        i_min = np.array(i_min)
        i_max = np.array(i_max)

        if self.verb:
            for i in range(len(i_shells)):
                print("\t PRINCIPLE JET")
                print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
                                                                                      xmax[i]))
                print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
                                                                                      ymax[i]))
            print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
            print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
            print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
        #
        # if (xmin==ymin==0):
        #     raise ValueError()
        # if (ymax==xmax==0):
        #     raise ValueError()
        return(xmin, xmax, ymin, ymax, i_shells, ncells)


    def _get_skymap(self, in_f : h5py.File, return_sph_coords : bool = False, remove_zeros : bool = True) \
            -> tuple[list[np.ndarray],list[np.ndarray],list[np.ndarray]]:
        """
        Read 'raw_skymap.h5' file and do the following:
            1. Extract X, Y, I data for each shell where sum(I)>0
            2. Assuming that ncells = len(X_i)/2, divide data into principle and counter jet data
                for each shell (doubling the number of shells in effect) WARNING this copies data
            3. Remove I=0 entries from arrays to simplify future analysis
        :param in_f:
        :param return_sph_coords:
        :param remove_zeros:
        :return: list of Xs, Ys, Intensity for each 2*nshells,
                 with [:nshells] for pricnciple and [nshells:] for counter jet
        """
        nshells = int(in_f.attrs["nshells"])

        # locate shells with data
        i_shells = []
        for ish in range(nshells):
            sumi = np.sum(np.array(in_f[f"shell={ish}"]["intensity"]))
            if (sumi == 0):
                if self.verb:
                    print(f"Skipping empty sum(intensity)=0 shell ish={ish}")
                continue
            i_shells.append(ish)

        # extract row data
        ncells = []
        all_xrs, all_yrs, all_zz = [], [], []
        all_theta, all_phi, all_r = [], [], []
        # loop over non-empy shells and collect data
        for ii, ish in enumerate(i_shells):
            # cartesian coordiantes on the projected plane
            xrs_i = np.array(in_f[f"shell={ish}"]["xrs"]) # mas
            yrs_i = np.array(in_f[f"shell={ish}"]["yrs"]) # mas
            int_i = np.array(in_f[f"shell={ish}"]["intensity"]) # mJy / mas^2

            all_xrs.append(xrs_i)
            all_yrs.append(yrs_i)
            all_zz.append(int_i)

            if (return_sph_coords):

                ctheta_i = np.array(in_f[f"shell={ish}"]["ctheta"])
                cphi_i = np.array(in_f[f"shell={ish}"]["cphi"])
                r_i = np.array(in_f[f"shell={ish}"]["r"])
                mu_i = np.array(in_f[f"shell={ish}"]["mu"])

                all_theta.append(ctheta_i)
                all_phi.append(cphi_i)
                all_r.append(r_i)

            ncells.append( int(len(xrs_i) / 2) )
            if (len(xrs_i) % 2 > 0):
                raise ValueError(f"len(xrs) is expected to be even (2*ncells). Got={len(xrs_i)}")

        pass
        # collect result into lists, removing zeross if needed
        all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj, maxs = [], [], [], []
        all_ctheta_pjcj, all_cphi_pjcj, all_r_pjcj = [], [], []
        for ii, ish in enumerate(i_shells):
            _xrs_pj = all_xrs[ii][:ncells[ii]]
            _yrs_pj = all_yrs[ii][:ncells[ii]]
            _zz_pj  = all_zz[ii][:ncells[ii]]
            maxs.append(np.max(_zz_pj))

            all_xrs_pjcj.append(_xrs_pj[ _zz_pj > 0 ] if (remove_zeros) else _xrs_pj )
            all_yrs_pjcj.append(_yrs_pj[ _zz_pj > 0 ] if (remove_zeros) else _yrs_pj )
            all_zz_pjcj.append(  _zz_pj[ _zz_pj > 0 ] if (remove_zeros) else _zz_pj )

            # if (len(all_xrs_pjcj[-1]) < 3):
            #     raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
            #                      f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_pj  = all_theta[ii][:ncells[ii]]
                _cphi_pj    = all_phi[ii][:ncells[ii]]
                _r_pj       = all_r[ii][:ncells[ii]]

                all_ctheta_pjcj.append(_ctheta_pj[ _zz_pj > 0 ] if (remove_zeros) else _ctheta_pj )
                all_cphi_pjcj.append(  _cphi_pj[   _zz_pj > 0 ] if (remove_zeros) else _cphi_pj )
                all_r_pjcj.append(     _r_pj[      _zz_pj > 0 ] if (remove_zeros) else _r_pj )

        if self.verb: print(f"Principle only max_intensity = {maxs}")

        # process counter jet
        for ii, ish in enumerate(i_shells):
            _xrs_cj = all_xrs[ii][ncells[ii]:]
            _yrs_cj = all_yrs[ii][ncells[ii]:]
            _zz_cj  = all_zz[ii][ncells[ii]:]
            maxs.append(np.max(_zz_cj))

            all_xrs_pjcj.append(_xrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _xrs_cj )
            all_yrs_pjcj.append(_yrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _yrs_cj )
            all_zz_pjcj.append(  _zz_cj[ _zz_cj > 0 ] if (remove_zeros) else _zz_cj )

            if (len(all_xrs_pjcj[-1]) < 3):
                raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
                                 f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_cj  = all_theta[ii][ncells[ii]:]
                _cphi_cj    = all_phi[ii][ncells[ii]:]
                _r_cj       = all_r[ii][ncells[ii]:]

                all_ctheta_pjcj.append(_ctheta_cj[ _zz_cj > 0 ] if (remove_zeros) else _ctheta_cj )
                all_cphi_pjcj.append(  _cphi_cj[   _zz_cj > 0 ] if (remove_zeros) else _cphi_cj )
                all_r_pjcj.append(     _r_cj[      _zz_cj > 0 ] if (remove_zeros) else _r_cj )

        if self.verb: print(f"Principle & counter only max_intensity = {maxs}")

        if (self.verb):
            print("Given Data From SkyMap file:")
            for ii, ish in enumerate(i_shells):
                if len(all_xrs_pjcj[ii]) > 0:
                    print(f"\tSHELL={ish} Principle: extend "
                          f"X=[{np.min(all_xrs_pjcj[ii])}, {np.max(all_xrs_pjcj[ii])}] "
                          f"Y=[{np.min(all_yrs_pjcj[ii])}, {np.max(all_yrs_pjcj[ii])}] "
                          f"Z=[{np.min(all_zz_pjcj[ii])}, {np.max(all_zz_pjcj[ii])}]")
                    nshells_ = int(len(all_xrs_pjcj) / 2)
                    # if (nshells+ii > len(all_xrs_pjcj)-1):
                    #     raise ValueError()

                    print(f"\tSHELL={ish} Counter: extend "
                          f"X=[{np.min(all_xrs_pjcj[nshells_+ii])}, {np.max(all_xrs_pjcj[nshells_+ii])}] "
                          f"Y=[{np.min(all_yrs_pjcj[nshells_+ii])}, {np.max(all_yrs_pjcj[nshells_+ii])}] "
                          f"Z=[{np.min(all_zz_pjcj[nshells_+ii])}, {np.max(all_zz_pjcj[nshells_+ii])}]")
                else:
                    print(f"\tSHELL={ish} EMPTY SHELL")


        if (return_sph_coords):
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj, all_ctheta_pjcj, all_cphi_pjcj, all_r_pjcj)
        else:
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj)

    def _lateral_distributions(self, grid_x, grid_y, image, collapse_axis='y', method="integ") -> np.ndarray:
        """ compute image extend along X and Y directions """
        if collapse_axis == 'x':
            points = image.shape[1]
            latAvDist = np.array([image[:, ii].mean() for ii in range(points)])
            # latAvDist2 = np.array([np.sum(image[:-1, ii]*np.diff(grid_y[:,ii]))/np.sum(np.diff(grid_y[:,ii])) for ii in range(points)])
            # latAvDist2 = np.array([np.sum(image[:-1, ii]*np.diff(np.abs(grid_y[:,ii])))/np.sum(np.diff(np.abs(grid_y[:,ii]))) for ii in range(points)])

            # latAvDist2 = np.zeros(points)
            # for ii in range(points):
            #     diff = grid_x[:, ii]
            #     dx = np.diff(diff)[0]
            #     arr = image[:, ii]
            #     latAvDist2[ii] = np.trapz(arr, diff, dx=dx)
            latAvDist2 = np.array( [ np.trapz(image[:, ii], grid_x, dx=np.diff(grid_x)[0] ) for ii in range(points) ] )

            # latAvDist2 = np.array([np.trapz(image[:, ii], grid_y[:,ii]) for ii in range(points)]) # , dx=np.abs(np.diff(grid_y[:,ii])[0])
            latMaxDist = np.array([image[:, ii].max() for ii in range(points)])

        elif collapse_axis == 'y':
            points = image.shape[0]
            latAvDist = np.array([image[ii, :].mean() for ii in range(points)])
            # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(grid_x[ii,:]))/np.sum(np.diff(grid_x[ii,:])) for ii in range(points)])
            # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(np.abs(grid_x[ii,:])))/np.sum(np.diff(np.abs(grid_x[ii,:]))) for ii in range(points)])
            latAvDist2 = np.array( [ np.trapz(image[ii, :], grid_y, dx=np.diff(grid_y)[0] ) for ii in range(points) ] )

            latMaxDist = np.array([image[ii, :].max() for ii in range(points)])

        else:
            raise KeyError()


        if method == "integ":
            return latAvDist2
        elif method == "ave":
            return latAvDist
        elif method == "max":
            return latMaxDist
        else:
            raise KeyError(f"method {method} is not recognized")

    def _total_hist_skymaps(self, xs, ys, datas, edges_x, edges_y) -> np.ndarray:
        # init the total intensity array for histogram and for interpolation
        zz_hist = np.zeros((len(edges_x) - 1, len(edges_y) - 1))
        for i, (x_i, y_i, data_i) in enumerate(zip(xs, ys, datas)):
            if self.verb: print(f"Histogram processing: shell={i} [{i}/{len(xs)}]")
            # binn the shell data onto a uniform grid
            if len(data_i) > 0:
                i_zz, _ = np.histogramdd(tuple([x_i, y_i]), bins=tuple([edges_x, edges_y]), weights=data_i)
                zz_hist += i_zz
            else:
                if self.verb:
                    print(f"Skipping shell i={i} in computing histogram. It is empty.")
        return zz_hist

    def _total_intep_skymap(self, xs, ys, datas, centers_X, centers_Y, method='linear') -> np.ndarray:
        # init the total intensity array for histogram and for interpolation
        zz_int = np.zeros_like(centers_X).T

        for i, (x_i, y_i, data_i) in enumerate(zip(xs, ys, datas)):
            if self.verb: print(f"Re-interpolating: shell={i} [{i}/{len(xs)}]")
            # interpolate the shell data onto uniform grid
            if len(data_i) > 0:
                i_zz = ProcessRawSkymap._interp(x_i, y_i, data_i, centers_X, centers_Y, method)
                # i_zz *= (len(x_i)*len(y_i)) / (len(centers_x)*len(centers_y))
                zz_int += i_zz
            if self.verb:
                print(f"Skipping shell i={i} in computing interpolation. It is empty.")
        return zz_int



    def _combine_images_from_shells(self, xs, ys, datas, edges_x, edges_y):

        # prepare grid for interpolation
        centers_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        centers_y = 0.5 * (edges_y[1:] + edges_y[:-1])

        nshells = int(len(xs) / 2)
        if (len(xs) % 2 > 0):
            raise ValueError(f"expected len(xs) to be even 2*nshells got={len(xs)}")

        # prep grid for the image

        if (self.verb):
            print("Given Data From SkyMap file:")
            for ish in range(nshells):
                if len(datas[ish]) > 0:
                    print(f"\tSHELL={ish} Principle: extend X=[{np.min(xs[ish])}, {np.max(xs[ish])}] "
                          f"Y=[{np.min(ys[ish])}, {np.max(ys[ish])}] "
                          f"Z=[{np.min(datas[ish])}, {np.max(datas[ish])}]")
                    print(f"\tSHELL={ish} Counter: extend X=[{np.min(xs[nshells+ish])}, {np.max(xs[nshells+ish])}] "
                          f"Y=[{np.min(ys[nshells+ish])}, {np.max(ys[nshells+ish])}] "
                          f"Z=[{np.min(datas[nshells+ish])}, {np.max(datas[nshells+ish])}]")
                else:
                    print(f"\tSHELL={ish} EMPTY SHELL")

        centers_X, centers_Y = np.meshgrid(centers_x, centers_y)
        # centers_X, centers_Y = centers_X.flatten(), centers_Y.flatten()

        # process principle jet
        zz_pj_intp = self._total_intep_skymap(xs[:nshells], ys[:nshells], datas[:nshells], centers_X, centers_Y)
        zz_pj_hist = self._total_hist_skymaps(xs[:nshells], ys[:nshells], datas[:nshells], edges_x, edges_y)

        # process counter jet
        zz_cj_intp = self._total_intep_skymap(xs[nshells:], ys[nshells:], datas[nshells:], centers_X, centers_Y)
        zz_cj_hist = self._total_hist_skymaps(xs[nshells:], ys[nshells:], datas[nshells:], edges_x, edges_y)


        if (self.verb):
            print("SkyMap After Interpolation")
            print(f"\tPrinciple: extend X=[{np.min(centers_x)}, {np.max(centers_x)}] "
                  f"Y=[{np.min(centers_y)}, {np.max(centers_y)}]")
            print(f"\tZ_hist=[{np.min(zz_pj_hist)}, {np.max(zz_pj_hist)}]"
                  f"Z_intp=[{np.min(zz_pj_intp)}, {np.max(zz_pj_intp)}]" )
            print(f"\tZ_hist=[{np.min(zz_cj_hist)}, {np.max(zz_cj_hist)}]"
                  f"Z_intp=[{np.min(zz_cj_intp)}, {np.max(zz_cj_intp)}]" )

        # sum the result
        if (zz_pj_intp.shape != zz_cj_intp.shape):
            raise ValueError("Shape mismatch")
        if (zz_pj_hist.shape != zz_cj_hist.shape):
            raise ValueError("Shape mismatch")
        if (zz_pj_intp.shape != zz_pj_hist.shape):
            raise ValueError("Shape mismatch")
        zz_intp = zz_pj_intp + zz_cj_intp
        zz_hist = zz_pj_hist + zz_cj_hist

        # get hist bin centers
        xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
        yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])

        return (xx_grid, yy_grid, zz_intp, zz_hist)

    def _process_file(self, in_f : h5py.File, group : h5py.Group,  edges_x=None, edges_y=None):
        # collate data from all shells from a given file
        xs, ys, datas = self._get_skymap(in_f)
        assert len(xs) == len(ys), f"shells must have the same xs{len(xs)} ys={len(ys)}"
        assert len(xs) == len(datas), f"shells must have the same xs{len(xs)} datas={len(datas)}"

        # compute sky map central mass using weighted average
        xc, yc = self._compute_position_of_the_flux_centroid(xs, ys, datas)

        # get the data extend
        xmin, xmax, ymin, ymax, i_shells, ncells = ( self._get_data_extend( xs, ys, datas ) )

        # prepare grid from the histogram (create new grid if needed)
        extend_grid = self.conf["extend_grid"]
        if edges_x is None:
            edges_x = np.linspace(xmin.min()*extend_grid, xmax.max()*extend_grid, num=self.conf["nx"])
        if edges_y is None:
            edges_y = np.linspace(ymin.min()*extend_grid, ymax.max()*extend_grid, num=self.conf["ny"])

        # combine images 1) after interpolation 2) after histogram binning
        grid_x, grid_y, image_intp, image_hist = \
            (self._combine_images_from_shells(xs, ys, datas, edges_x, edges_y))

        # apply smoothing kernal to the images
        image_intp = self._apply_filter(image_intp, self.conf_intp_filter)
        image_hist = self._apply_filter(image_hist, self.conf_hist_filter)

        # compute image X and Y extend
        dist_x = self._lateral_distributions(
            grid_x, grid_y, image_intp, collapse_axis='y', method=self.conf["lat_dist_method"])
        dist_y = self._lateral_distributions(
            grid_x, grid_y, image_intp, collapse_axis='x', method=self.conf["lat_dist_method"])

        # compute skymap size ax Xc and Yc from its X and Y extend (using integrated distributions)
        x1, x2 = self._get_skymap_fwhm(grid_x, dist_x, xc, fac=self.conf["fwhm_fac"])
        y1, y2 = self._get_skymap_fwhm(grid_y, dist_y, yc, fac=self.conf["fwhm_fac"])

        # save the result in the file

        group.create_dataset("image_hist", data=image_hist)
        group.create_dataset("image_intp", data=image_intp)
        group.create_dataset("grid_x", data=grid_x)
        group.create_dataset("grid_y", data=grid_y)
        group.create_dataset("dist_x", data=dist_x)
        group.create_dataset("dist_y", data=dist_y)
        # add atributes to the file
        group.attrs.create("xc", data=xc)
        group.attrs.create("yc", data=yc)
        group.attrs.create("x1", data=x1)
        group.attrs.create("x2", data=x2)
        group.attrs.create("y1", data=y1)
        group.attrs.create("y2", data=y2)
        # copy attributes from the source file (model parameters)
        for key in in_f.attrs.keys():
            group.attrs.create(key, data=in_f.attrs[key])

        # in case there was a removal of a shell (empty)
        group.attrs["nshells"] = len(datas) / 2

    def process_singles(self, infpaths="raw_skymap_*.h5", outfpath="skymap.h5",remove_input=False,
                        edges_x = None, edges_y = None):
        """
        # Kursbuch auf Seite 92 und 93
            # for each row file (for each time and frequency) 390 ubung 4 Ein frage, Was meininen Sie? In ihrem ort wolk haben?
            # Mit welhes afolk haben? Were Hatter Wurde benutzen?
            # Text for ubung 4
            # Halb zeite
            # Kursbuch auf Seite 92 und 93
            # Mit welcher Geschäftsidee könnte man als Existenzgründer in Potsdam oder Berlin Erfolg haben?
            # Erörtere die Chancen und Risiken!

        :param infpaths:
        :param outfpath:
        :param remove_input:
        :param edges_x:
        :param edges_y:
        :return:
        """

        files = glob(infpaths)
        if (len(files)==0):
            raise FileNotFoundError(f" Not raw skympas found. Searching: {infpaths}")
        files = sorted(files, key=lambda key : int(re.findall(r'\d+', key)[-2]))
        # check if output file already exists
        if (os.path.isfile(outfpath)):
            if self.verb: print(f"File already exists {outfpath}")

        out_f = h5py.File(outfpath,"w")
        times = []
        freqs = []
        for fl in files:
            print(f"Processing {fl}")
            in_f = h5py.File(fl, "r")

            time = float(in_f.attrs["time"])
            freq = float(in_f.attrs["freq"])
            key = "time={:.4e} freq={:.4e}".format(time, freq)
            group = out_f.create_group(key)
            self._process_file(in_f=in_f, group=group, edges_x=edges_x, edges_y=edges_y)
            times.append(time)
            freqs.append(freq)

            in_f.close()
        out_f.create_dataset("times", data=np.array(times))
        out_f.create_dataset("freqs", data=np.array(freqs))
        out_f.close()

        # remove input files (save space)
        if remove_input:
            for fl in files:
                if self.verb:
                    print(f"Deleting file... {fl}")
                os.remove(fl)
