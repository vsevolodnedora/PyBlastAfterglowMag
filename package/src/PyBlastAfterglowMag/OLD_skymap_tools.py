import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy

from .utils import cgs, get_beta, find_nearest_index

def compute_position_of_the_flux_centroid(all_x_arrs, all_y_arrs, all_z_arrs, d_l):
    rad2mas = 2.062648062470964e+08
    xcs = np.average(np.concatenate(all_x_arrs), weights=np.concatenate(all_z_arrs))
    xcs_m = xcs# * rad2mas / d_l
    ycs = np.average(np.concatenate(all_y_arrs), weights=np.concatenate(all_z_arrs))
    ycs_m = ycs# * rad2mas / d_l
    return(xcs_m, ycs_m) # in milli-arc-seconds

def interp(xxs, yys, fluxes, x_grid, y_grid, method='linear'):
    return interpolate.griddata(np.vstack((xxs, yys)).T, fluxes, np.array((x_grid, y_grid)).T, method=method, fill_value=0.)

def lateral_distributions(grid_x, grid_y, image, collapse_axis='y'):
    # if collapse_axis == "y":
    #     return (grid_x[:, 0], latAvDist, latMaxDist)
    # else:
    #     return (grid_y[0, :], latAvDist, latMaxDist)
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
        latAvDist2 = np.array( [ np.trapz(image[:, ii], grid_x[:, ii], dx=np.diff(grid_x[:, ii])[0] ) for ii in range(points) ] )

        # latAvDist2 = np.array([np.trapz(image[:, ii], grid_y[:,ii]) for ii in range(points)]) # , dx=np.abs(np.diff(grid_y[:,ii])[0])
        latMaxDist = np.array([image[:, ii].max() for ii in range(points)])

    elif collapse_axis == 'y':
        points = image.shape[0]
        latAvDist = np.array([image[ii, :].mean() for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(grid_x[ii,:]))/np.sum(np.diff(grid_x[ii,:])) for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(np.abs(grid_x[ii,:])))/np.sum(np.diff(np.abs(grid_x[ii,:]))) for ii in range(points)])
        latAvDist2 = np.array( [ np.trapz(image[ii, :], grid_y[ii, :], dx=np.diff(grid_y[ii, :])[0] ) for ii in range(points) ] )

        latMaxDist = np.array([image[ii, :].max() for ii in range(points)])

    else:
        raise KeyError()

    return (latAvDist, latAvDist2, latMaxDist)

def image_slice(fluxes, xxs, yys, position, nn=50, axis='y', inter='linear', fac=1):
    nn = complex(0, nn)
    if axis == 'y':
        grid_x, grid_y = np.mgrid[position:position:1j, yys.min():yys.max():nn]
    else:
        grid_x, grid_y = np.mgrid[xxs.min():xxs.max():nn, position:position:1j]

    slice = interpolate.gdd(np.array([xxs, yys]).T * fac, fluxes,
                            (grid_x * fac, grid_y * fac), method=inter, fill_value=0)

    return grid_x, grid_y, slice








def TOREMOVE_combine_images_old(xs, ys, datas, verbose=False, hist_or_int="int", shells=False, nx=200, ny=100, extend=2,
                       retrun_edges=False, edges_x=None, edges_y=None):

    if shells:
        assert len(xs) == len(ys)
        assert len(xs) == len(datas)

        nshells = len(xs)


        xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
        ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
        i_min, i_max = [], []
        i_shells = []
        for ish in range(nshells):
            xrs_i = np.array(xs[ish])
            yrs_i = np.array(ys[ish])
            int_i = np.array(datas[ish])

            # skip empty shells
            if ((np.sum(int_i) == 0)):
                continue

            if (len(xrs_i) % 2 > 0):
                raise ValueError(f"expected to get an even number for ncells. Got:{ len(xrs_i) }")
            ncells = int( len(xrs_i) / 2 )

            if (verbose):
                print(f"Shell = {ish}/{nshells} ncells={ncells}")
                print(f"Principle: extend X=[{np.min(xrs_i[:ncells])}, {np.max(xrs_i[:ncells])}] "
                      f"Y=[{np.min(yrs_i[:ncells])}, {np.max(yrs_i[:ncells])}]"
                      f"Z=[{np.min(int_i[:ncells])}, {np.max(int_i[:ncells])}]")
                print(f"Counter: extend X=[{np.min(xrs_i[ncells:])}, {np.max(xrs_i[ncells:])}] "
                      f"Y=[{np.min(yrs_i[ncells:])}, {np.max(yrs_i[ncells:])}]"
                      f"Z=[{np.min(int_i[ncells:])}, {np.max(int_i[ncells:])}]")

            i_shells.append(ish)
            xmin_neg.append( xrs_i[xrs_i < 0].min() if len(xrs_i[xrs_i < 0]) > 0 else 0 )
            xmin_pos.append( xrs_i[xrs_i > 0].min() if len(xrs_i[xrs_i > 0]) > 0 else 0 )
            xmax.append( xrs_i.max() )
            xmin.append( xrs_i.min() )
            ymin_neg.append( yrs_i[yrs_i < 0].min() if len(yrs_i[yrs_i < 0]) > 0 else 0 )
            ymin_pos.append( yrs_i[yrs_i > 0].min() if len(yrs_i[yrs_i > 0]) > 0 else 0 )
            ymax.append(yrs_i.max())
            ymin.append(yrs_i.min())
            i_min.append(int_i.min())
            i_max.append(int_i.max())
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
        if verbose:
            for i in range(len(i_shells)):
                print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
                                                                                      xmax[i]))
                print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
                                                                                      ymax[i]))
            print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
            print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
            print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
        x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
                                 np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
        y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
                                 np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))

        if edges_x is None: edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx)
        if edges_y is None: edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny)

        if verbose:
            print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                                                                              x_grid[x_grid > 0].min(),
                                                                              x_grid.max()))
            print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                                                                              y_grid[y_grid > 0].min(),
                                                                              y_grid.max()))
        xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
        if hist_or_int == "int":
            zz = np.zeros_like((xx_grid))
        else:
            zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))
        # interpolate onto the grid that covers all images
        all_xrs, all_yrs, all_zz = [], [], []
        for ii, ish in enumerate(i_shells):  # range(len(i_shells))
            if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
            xrs_i = np.array(xs[ish])
            yrs_i = np.array(ys[ish])
            int_i = np.array(datas[ish])  # * dfile.attrs["d_L"] ** 2

            if hist_or_int == "int":
                i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
            else:
                i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)

            zz += i_zz
        xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
        yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])

        if retrun_edges:
            return (xx_grid, yy_grid, zz, edges_x, edges_y)
        else:
            return (xx_grid, yy_grid, zz)
    else:

        xs = np.concatenate(xs)
        ys = np.concatenate(ys)
        datas = np.concatenate(datas)

        if hist_or_int == "int":
            nx = np.complex(0, nx)
            ny = np.complex(0, ny)
            grid_x, grid_y = np.mgrid[xs.min()*extend:xs.max()*extend:nx, ys.min()*extend:ys.max()*extend:ny]
            i_zz = interp(xs, ys, datas, grid_x, grid_y, 'linear')
            i_zz = i_zz
        elif hist_or_int == "both":

            nx = np.complex(0, nx)
            ny = np.complex(0, ny)
            grid_x, grid_y = np.mgrid[xs.min():xs.max():nx, ys.min():ys.max():ny]
            i_zz = interp(xs, ys, datas, grid_x, grid_y, 'linear')

            nx = np.complex(0, nx + 1)
            ny = np.complex(0, ny + 1)
            edges_x = np.mgrid[xs.min():xs.max():nx]
            edges_y = np.mgrid[ys.min():ys.max():ny]
            # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
            #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
            i_zz, _ = np.histogramdd(tuple([grid_x.flatten(), grid_y.flatten()]), bins=tuple([edges_x, edges_y]),
                                     weights=i_zz.T.flatten())
            grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
            grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
            # print(i_zz.shape)
        else:

            # nx = 200
            # ny = 100
            nx = np.complex(0, nx + 1)
            ny = np.complex(0, ny + 1)
            edges_x = np.mgrid[xs.min()*extend:xs.max()*extend:nx]
            edges_y = np.mgrid[ys.min()*extend:ys.max()*extend:ny]
            # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
            #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]

            xs = xs[datas > 0]
            ys = ys[datas > 0]
            datas = datas[datas>0]

            i_zz, _ = np.histogramdd(tuple([xs, ys]), bins=tuple([edges_x, edges_y]), weights=datas,normed='density')
            grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
            grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
            print(i_zz.shape)

        # i_zz = ndimage.uniform_filter(i_zz, )
        # i_zz = ndimage.filters.gaussian_filter(i_zz, [10,10], mode='reflect')

        return (grid_x, grid_y, i_zz)
    #
    # assert len(xs) == len(ys)
    # assert len(xs) == len(datas)
    # nshells = len(xs)
    # nx = 1000
    # ny = 1000
    # xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
    # ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
    # i_min, i_max = [], []
    # i_shells = []
    # for ish in range(nshells):
    #     xrs_i = np.array(xs[ish])
    #     yrs_i = np.array(ys[ish])
    #     int_i = np.array(datas[ish])
    #     if (np.sum(int_i) == 0):
    #         continue
    #     i_shells.append(ish)
    #     xmin_neg.append(xrs_i[xrs_i < 0].min())
    #     xmin_pos.append(xrs_i[xrs_i > 0].min())
    #     xmax.append(xrs_i.max())
    #     xmin.append(xrs_i.min())
    #     ymin_neg.append(yrs_i[yrs_i < 0].min())
    #     ymin_pos.append(yrs_i[yrs_i > 0].min())
    #     ymax.append(yrs_i.max())
    #     ymin.append(yrs_i.min())
    #     i_min.append(int_i.min())
    #     i_max.append(int_i.max())
    # xmin_neg = np.array(xmin_neg)
    # xmin_pos = np.array(xmin_pos)
    # xmax = np.array(xmax)
    # ymin_neg = np.array(ymin_neg)
    # ymin_pos = np.array(ymin_pos)
    # ymax = np.array(ymax)
    # i_min = np.array(i_min)
    # i_max = np.array(i_max)
    # if verbose:
    #     for i in range(len(i_shells)):
    #         print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
    #                                                                               xmax[i]))
    #         print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
    #                                                                               ymax[i]))
    #     print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
    #     print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
    #     print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
    # x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
    #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
    # y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
    #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
    #
    # edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
    # edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
    #
    # if verbose:
    #     print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
    #                                                                       x_grid[x_grid > 0].min(),
    #                                                                       x_grid.max()))
    #     print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
    #                                                                       y_grid[y_grid > 0].min(),
    #                                                                       y_grid.max()))
    # xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    # if hist_or_int == "int":
    #     zz = np.zeros_like((xx_grid))
    # else:
    #     zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))
    # # interpolate onto the grid that covers all images
    # all_xrs, all_yrs, all_zz = [], [], []
    # for ii, ish in enumerate(i_shells):  # range(len(i_shells))
    #     if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
    #     xrs_i = np.array(xs[ish])
    #     yrs_i = np.array(ys[ish])
    #     int_i = np.array(datas[ish])  # * dfile.attrs["d_L"] ** 2
    #
    #     if hist_or_int == "int":
    #         # nx = np.complex(0, nx)
    #         # ny = np.complex(0, ny)
    #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #         # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #     else:
    #         # nx = 100
    #         # ny = 100
    #         # nx = np.complex(0, nx + 1)
    #         # ny = np.complex(0, ny + 1)
    #         # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #         # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #         i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #         # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #         # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #         # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #
    #     zz += i_zz
    #     all_xrs.append(xrs_i)
    #     all_yrs.append(yrs_i)
    #     all_zz.append(int_i)
    #     #
    #     # if hist_or_int == "int":
    #     #     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #     #     zz += i_zz
    #     #     all_xrs.append(xrs_i)
    #     #     all_yrs.append(yrs_i)
    #     #     all_zz.append(int_i)
    #     # else:
    #     #     nx = 2000
    #     #     ny = 1000
    #     #     nx = np.complex(0, nx + 1)
    #     #     ny = np.complex(0, ny + 1)
    #     #     edges_x = np.mgrid[x_grid.min():x_grid.max():nx]
    #     #     edges_y = np.mgrid[y_grid.min():yrs_i.max():ny]
    #     #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #     #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #     #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #     #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #     #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #     #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    # # if verbose:
    #     # print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
    #     #       .format(zz.min(), zz.max(), np.sum(zz), float(
    #     #     np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))

    return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)

def TOREMOVE_combine_images_old_new(xs, ys, datas, verbose=False, hist_or_int="int", shells=False, nx=200, ny=100, extend=2,
                           retrun_edges=False, edges_x=None, edges_y=None):

    if (not shells):
        raise KeyError(" shells should be true (delete this keyword) ")

    assert len(xs) == len(ys), f"shells must have the same xs{len(xs)} ys={len(ys)}"
    assert len(xs) == len(datas), f"shells must have the same xs{len(xs)} datas={len(datas)}"

    # principle jet
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

        if (len(xrs_i) % 2 > 0):
            raise ValueError(f"expected to get an even number for ncells. Got:{len(xrs_i)}")

        ncells.append( int( len(xrs_i) / 2 ) )
        i_shells.append(ish)

        # xrs_i = xrs_i
        # yrs_i = yrs_i[:ncells]
        # int_i = int_i[:ncells]

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

    if verbose:
        for i in range(len(i_shells)):
            print("\t PRINCIPLE JET")
            print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
                                                                                  xmax[i]))
            print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
                                                                                  ymax[i]))
        print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
        print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
        print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))

    # create grid that covers the image
    x_grid = np.concatenate((
        -1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
        np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx))
    )
    y_grid = np.concatenate((
        -1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
        np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny))
    )

    # create edges got hist
    if edges_x is None: edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx)
    if edges_y is None: edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny)
    if verbose:
        print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                                                                          x_grid[x_grid > 0].min(),
                                                                          x_grid.max()))
        print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                                                                          y_grid[y_grid > 0].min(),
                                                                          y_grid.max()))
    xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

    centers_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    centers_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    centers_x_, centers_y_ = np.meshgrid(centers_x, centers_y)
    # edges_x_, edges_y_ = np.meshgrid(edges_x, edges_y)


    # if hist_or_int == "int": zz = np.zeros_like((xx_grid))
    # else:
    zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))

    # interpolate onto the grid that covers all images
    for ii, ish in enumerate(i_shells):
        if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))

        # process principle jet

        ncells_pj = ncells[ii]
        xrs_i = np.array(xs[ish][:ncells_pj])
        yrs_i = np.array(ys[ish][:ncells_pj])
        int_i = np.array(datas[ish][:ncells_pj])  # * dfile.attrs["d_L"] ** 2
        if (np.sum(int_i)>0):
            if (hist_or_int == "int"):
                i_zz_pj = interp(xrs_i, yrs_i, int_i, centers_x_, centers_y_, 'linear')
                # i_zz_pj = np.reshape(i_zz_pj, newshape=zz.shape)
            elif (hist_or_int == "hist"):
                i_zz_pj, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
            elif (hist_or_int == "both"):
                k = 2
                edges_x_int = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx*k)
                edges_y_int = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny*k)
                edges_x_int, edges_y_int = np.meshgrid(edges_x_int, edges_y_int)
                i_zz_pj = interp(xrs_i, yrs_i, int_i, edges_x_int, edges_y_int, 'linear')

                edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx)
                edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny)
                # edges_x_, edges_y_ = np.meshgrid(edges_x_int, edges_y_int)
                i_zz_pj, _ = np.histogramdd(tuple([edges_x_int.flatten(), edges_y_int.flatten()]), bins=tuple([edges_x, edges_y]), weights=i_zz_pj.T.flatten())
            else:
                raise KeyError(f"hist_or_int={hist_or_int} is not recognized")
            zz += i_zz_pj

        # process counter jet

        xrs_i = np.array(xs[ish][ncells_pj:])
        yrs_i = np.array(ys[ish][ncells_pj:])
        int_i = np.array(datas[ish][ncells_pj:])  # * dfile.attrs["d_L"] ** 2
        if (np.sum(int_i)>0):
            if (hist_or_int == "int"):
                i_zz_cj = interp(xrs_i, yrs_i, int_i, centers_x_, centers_y_, 'linear')
            elif (hist_or_int == "hist"):
                i_zz_cj, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
            elif (hist_or_int == "both"):
                k = 2
                edges_x_int = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx*k)
                edges_y_int = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny*k)
                edges_x_int, edges_y_int = np.meshgrid(edges_x_int, edges_y_int)
                i_zz_cj = interp(xrs_i, yrs_i, int_i, edges_x_int, edges_y_int, 'linear')

                edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx)
                edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny)
                # edges_x_, edges_y_ = np.meshgrid(edges_x_int, edges_y_int)
                i_zz_cj, _ = np.histogramdd(tuple([edges_x_int.flatten(), edges_y_int.flatten()]), bins=tuple([edges_x, edges_y]), weights=i_zz_cj.T.flatten())
            else:
                raise KeyError(f"hist_or_int={hist_or_int} is not recognized")
            zz += i_zz_cj

    xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
    yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])
    if retrun_edges:
        return (xx_grid, yy_grid, zz, edges_x, edges_y)
    else:
        return (xx_grid, yy_grid, zz)

def _make_grid(xmin_pos, xmax, nx, ymin_pos, ymax, ny, log_grid=True):
    if (log_grid):
        x_grid = np.concatenate((
            -1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
            np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx))
        )
        y_grid = np.concatenate((
            -1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
            np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny))
        )
    else:
        x_grid = np.concatenate((
            -1. * np.logspace(xmin_pos.min(), xmax.max(), nx)[::-1],
            np.logspace(xmin_pos.min(), xmax.max(), nx))
        )
        y_grid = np.concatenate((
            -1. * np.logspace(ymin_pos.min(), ymax.max(), ny)[::-1],
            np.logspace(ymin_pos.min(), ymax.max(), ny))
        )
    return (x_grid, y_grid)
def _get_data_extend(xs, ys, datas, type="principle", verbose=False):
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

    if verbose:
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

def _reinterpolate_jet(xs, ys, datas, verbose=False, nx=200, ny=100, type="principle",
                       extend=2., pre_int=True, k = 2, edges_x=None, edges_y=None):

    # get data extend
    xmin, xmax, ymin, ymax, i_shells, ncells = _get_data_extend(xs, ys, datas, type=type, verbose=verbose)

    # prep grid for pre-interpolation (depth controlled by k)
    extend_ = 1.
    edges_x_int = np.linspace(xmin.min()*extend_, xmax.max()*extend_, num=nx*k)
    edges_y_int = np.linspace(ymin.min()*extend_, ymax.max()*extend_, num=ny*k)
    edges_x_int, edges_y_int = np.meshgrid(edges_x_int, edges_y_int)

    # prep the grid for histogram binning
    k = 1
    if edges_x is None: edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx*k)
    if edges_y is None: edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny*k)

    # init results
    zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))

    for ii, ish in enumerate(i_shells):
        if verbose: print(f"Pocessing: shell={ish} [{ii}/{len(i_shells)}] (Pre_int={pre_int} k={k})")
        ncells_pj = ncells[ii]
        # get the data
        if (type=="principle"):
            xrs_i = np.array(xs[ish][:ncells_pj])
            yrs_i = np.array(ys[ish][:ncells_pj])
            int_i = np.array(datas[ish][:ncells_pj])
        elif (type=="counter"):
            xrs_i = np.array(xs[ish][ncells_pj:])
            yrs_i = np.array(ys[ish][ncells_pj:])
            int_i = np.array(datas[ish][ncells_pj:])
        elif (type=="both"):
            xrs_i = np.array(xs[ish])
            yrs_i = np.array(ys[ish])
            int_i = np.array(datas[ish])
        else:
            raise KeyError(f"type={type} is not recognized")

        if (pre_int):
            # pre-inteprolate the data

            i_zz = interp(xrs_i, yrs_i, int_i, edges_x_int, edges_y_int, 'linear') * (len(xrs_i)*len(yrs_i)) / (len(edges_x_int[:,0])*len(edges_y_int[0,:]))

            # bin the data onto a histogram
            i_zz, _ = np.histogramdd(tuple([edges_x_int.flatten(), edges_y_int.flatten()]),
                                     bins=tuple([edges_x, edges_y]), weights=i_zz.T.flatten())
        else:
            # bin the data onto a histogram
            i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)

        # add results
        zz += i_zz
    return (edges_x, edges_y, zz)


def combine_images(xs, ys, datas, verbose=False, pre_int=True, k=2, nx=200, ny=100, extend=2,
                   retrun_edges=False, edges_x=None, edges_y=None):

    assert len(xs) == len(ys), f"shells must have the same xs{len(xs)} ys={len(ys)}"
    assert len(xs) == len(datas), f"shells must have the same xs{len(xs)} datas={len(datas)}"

    # get total extend of the data
    xmin, xmax, ymin, ymax, i_shells, ncells = _get_data_extend(
        xs, ys, datas, type="both", verbose=False)

    nshells = int(len(xs) / 2)
    if (len(xs) % 2 > 0):
        raise ValueError(f"expected len(xs) to be even 2*nshells got={len(xs)}")

    # prep grid for the image
    k_ = 1
    if edges_x is None: edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx*k_)
    if edges_y is None: edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny*k_)

    if (verbose):
        print("Given Data From SkyMap file:")
        for ish in range(nshells):
            print(f"\tSHELL={ish} Principle: extend X=[{np.min(xs[ish])}, {np.max(xs[ish])}] "
                  f"Y=[{np.min(ys[ish])}, {np.max(ys[ish])}] "
                  f"Z=[{np.min(datas[ish])}, {np.max(datas[ish])}]")
            print(f"\tSHELL={ish} Counter: extend X=[{np.min(xs[nshells+ish])}, {np.max(xs[nshells+ish])}] "
                  f"Y=[{np.min(ys[nshells+ish])}, {np.max(ys[nshells+ish])}] "
                  f"Z=[{np.min(datas[nshells+ish])}, {np.max(datas[nshells+ish])}]")

    # process principle jet
    edges_x_pj, edges_y_pj, zz_pj = _reinterpolate_jet(
        xs[:nshells], ys[:nshells], datas[:nshells], verbose=verbose,
        nx=nx, ny=ny, type="both", extend=extend, pre_int=pre_int, k=k,
        edges_x=edges_x, edges_y=edges_y
    )

    # process counter jet
    edges_x_cj, edges_y_cj, zz_cj = _reinterpolate_jet(
        xs[nshells:], ys[nshells:], datas[nshells:], verbose=verbose,
        nx=nx, ny=ny, type="both", extend=extend, pre_int=pre_int, k=k,
        edges_x=edges_x, edges_y=edges_y
    )

    if (verbose):
        print("SkyMap After Interpolation")
        print(f"\tPrinciple: extend X=[{np.min(edges_x_pj)}, {np.max(edges_x_pj)}] "
              f"Y=[{np.min(edges_y_pj)}, {np.max(edges_y_pj)}]"
              f"Z=[{np.min(zz_pj)}, {np.max(zz_pj)}]")
        print(f"\tCounter: extend X=[{np.min(edges_x_cj)}, {np.max(edges_x_cj)}] "
              f"Y=[{np.min(edges_y_cj)}, {np.max(edges_y_cj)}]"
              f"Z=[{np.min(zz_cj)}, {np.max(zz_cj)}]")

    # sum the result
    if (zz_pj.shape != zz_cj.shape):
        raise ValueError("Shape mismatch")
    zz = zz_pj + zz_cj

    # get hist bin centers
    xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
    yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])
    if retrun_edges:
        return (xx_grid, yy_grid, zz, edges_x, edges_y)
    else:
        return (xx_grid, yy_grid, zz)




# def get_ej_skymap(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False, int_or_hist="int"):
#
#
#     nx = 200
#     ny = 100
#     min_mu = 1e-4 # limit for mu that if too small blow up the intensity for some reason...
#     min_mu_frac = 0.85
#     # limit_data = True
#     times = self.get_ej_skymap_times()
#     freqs = self.get_ej_skymap_freqs()
#     dfile = self.get_ej_skymap_obj()
#     nshells = int(dfile.attrs["nshells"])
#     d_l = float(self.get_ej_skymap_obj().attrs["d_L"])
#     if ((not time is None)and(not time in times)):
#         raise ValueError("time={} is not in the list for skypams={}".format(time, times))
#     if((not freq is None)and(not freq in freqs)):
#         raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
#     if ((not ishell is None)and(ishell > nshells-1)):
#         raise ValueError("shell={} is beying skymap nshells={}".format(ishell, nshells))
#     if ((not time is None) and (not freq is None)):
#         print(dfile.keys())
#         ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
#         if (not ishell is None):
#             r_i = np.array(ddfile["r"][ishell])
#             mu_i = np.array(ddfile["mu"][ishell ])
#             xrs_i = np.array(ddfile["xrs"][ishell]) * cgs.rad2mas / d_l # m -> mas
#             yrs_i = np.array(ddfile["yrs"][ishell]) * cgs.rad2mas / d_l # m -> mas
#             int_i = np.array(ddfile["intensity"][ishell]) * (d_l ** 2 / cgs.rad2mas ** 2) #  -> mJy / mas^2
#             gam_i = np.array(ddfile["gamma"][ishell])
#             B_i = np.array(ddfile["B"][ishell])
#             tb_i = np.array(ddfile["tburst"][ishell])
#             theta_i = np.array(ddfile["theta"][ishell])
#             phi_i = np.array(ddfile["phi"][ishell])
#
#             if remove_mu:
#                 int_i *= np.abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
#
#             # plt.figure(figsize=(4.6, 3.2))
#             # plt.semilogy(mu_i, int_i / int_i.max(), '.')
#             # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
#             # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
#             # plt.xlabel(
#             #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
#             # plt.ylabel(r"$I/I_{\rm max}$")
#             # plt.tight_layout()
#             # plt.show()
#
#             # nx = np.complex(0, nx)
#             # ny = np.complex(0, ny)
#             # grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx, yrs_i.min():yrs_i.max():ny]
#             # i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
#             # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
#
#             if int_or_hist == "int":
#                 nx = np.complex(0, nx)
#                 ny = np.complex(0, ny)
#                 grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
#                                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                 i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
#                 return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
#             else:
#                 nx = 100
#                 ny = 100
#                 nx = np.complex(0, nx + 1)
#                 ny = np.complex(0, ny + 1)
#                 edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
#                 edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                 # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
#                 #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                 i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
#                 grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
#                 grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
#                 return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
#                 # print(i_zz.shape)
#
#         else:
#             # assess what is the combined grid extend (for all images)
#             xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
#             ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
#             i_min, i_max = [], []
#             i_shells = []
#             for ish in range(nshells):
#                 xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
#                 yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
#                 int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
#                 if (np.sum(int_i) == 0):
#                     continue
#                 i_shells.append(ish)
#                 xmin_neg.append(xrs_i[xrs_i < 0].min())
#                 xmin_pos.append(xrs_i[xrs_i > 0].min())
#                 xmax.append(xrs_i.max())
#                 xmin.append(xrs_i.min())
#                 ymin_neg.append(yrs_i[yrs_i < 0].min())
#                 ymin_pos.append(yrs_i[yrs_i > 0].min())
#                 ymax.append(yrs_i.max())
#                 ymin.append(yrs_i.min())
#                 i_min.append(int_i.min())
#                 i_max.append(int_i.max())
#             xmin_neg = np.array(xmin_neg)
#             xmin_pos = np.array(xmin_pos)
#             xmax = np.array(xmax)
#             xmin = np.array(xmin)
#             ymin_neg = np.array(ymin_neg)
#             ymin_pos = np.array(ymin_pos)
#             ymax = np.array(ymax)
#             ymin = np.array(ymin)
#             i_min = np.array(i_min)
#             i_max = np.array(i_max)
#
#             if verbose:
#                 for i in range(len(i_shells)):
#                     print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i], xmax[i]))
#                     print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],  ymax[i]))
#                 print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
#                 print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
#                 print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
#
#             edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
#             edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
#
#
#             x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
#                                      np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
#             y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
#                                      np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
#
#             # edges_x = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), 101)[::-1],
#             #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
#             # edges_y = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), 101)[::-1],
#             #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
#             # plt.plot(edges_x, edges_y, marker='.', ls='none')
#             # plt.show()
#
#             if verbose:
#                 print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
#                                                                                 x_grid[x_grid > 0].min(),
#                                                                                 x_grid.max()))
#                 print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
#                                                                                 y_grid[y_grid > 0].min(),
#                                                                                 y_grid.max()))
#             xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
#             if int_or_hist == "int":
#                 zz = np.zeros_like((xx_grid))
#             else:
#                 zz = np.zeros((len(edges_x)-1,len(edges_y)-1))
#             # interpolate onto the grid that covers all images
#             all_xrs, all_yrs, all_zz = [], [], []
#             for ii, ish in enumerate(i_shells):  # range(len(i_shells))
#                 if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
#                 r_i = np.array(ddfile["r"][ish])
#                 mu_i = np.array(ddfile["mu"][ish])
#                 xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
#                 yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
#                 int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2) # * dfile.attrs["d_L"] ** 2
#                 gam_i = np.array(ddfile["gamma"][ish])
#                 B_i = np.array(ddfile["B"][ish])
#                 tb_i = np.array(ddfile["tburst"][ish])
#                 theta_i = np.array(ddfile["theta"][ish])
#                 theta_ij = np.array(ddfile["theta_j"][ish])
#                 theta0_i = np.array(ddfile["theta0"][ish])
#                 phi_i = np.array(ddfile["phi"][ish])
#
#
#                 # plt.plot(mu_i, theta0_i, marker='.', ls='none')
#                 # plt.show()
#                 # int_i[mu_i < 5e-2] = 0.
#                 # int_i *= abs(mu_i)
#
#                 if remove_mu:
#                     # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
#                     # int_i[idx1] =
#                     # print("len(gam[idx1])={}".format(gam_i[idx1]))
#                     # idx = int_i > int_i.max() * min_mu_frac
#                     # idx = idx1
#                     # if verbose:
#                     #     if (len(mu_i[idx1]) != len(mu_i[idx])):
#                     #         print('\t', mu_i[idx1])
#                     #         print('\t', mu_i[idx])
#                     #         # exit(1)
#                     # if verbose:
#                     #     if len(idx) < 2:
#                     #         print("No excess = {}".format(ii))
#                     #     print("len(gam[idx])={}".format(gam_i[idx]))
#                     #     print("Gamma={}".format(gam_i[idx]))
#                     #     print("beta={}".format(get_beta(gam_i[idx])))
#                     #     print("mu_i={}".format(mu_i[idx]))
#                     #     print("r_i={}".format(r_i[idx]))
#                     #     print("B={}".format(B_i[idx]))
#                     #     print("tb_i={}".format(tb_i[idx]))
#                     # print(np.abs(mu_i).min())
#                     # print(theta_i[np.abs(mu_i) < 1e-5])
#                     # print(phi_i[np.abs(mu_i) < 1e-5])
#                     # print(np.abs(mu_i).min())
#
#                     # int_i[mu_i < 1e-3] = 0.
#                     int_i *= abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
#                     # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
#                     # plt.figure(figsize=(4.6,3.2))
#                     # plt.semilogy(mu_i,int_i/int_i.max(), '.')
#                     # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
#                     # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
#                     # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
#                     # plt.ylabel(r"$I/I_{\rm max}$")
#                     # plt.tight_layout()
#                     # plt.show()
#                     # int_i[idx] = 0.
#                 # int_i /= abs(mu_i)[abs(mu_i)>0.01]
#                 # else:
#                 #     int_i2 = int_i
#
#                 if int_or_hist == "int":
#                     # nx = np.complex(0, nx)
#                     # ny = np.complex(0, ny)
#                     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
#                     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
#                     #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
#                 else:
#                     # nx = 100
#                     # ny = 100
#                     # nx = np.complex(0, nx + 1)
#                     # ny = np.complex(0, ny + 1)
#                     # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
#                     # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
#                     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
#                     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
#                     # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
#                     # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
#                     #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
#
#                 # i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
#                 zz += i_zz
#                 all_xrs.append(xrs_i)
#                 all_yrs.append(yrs_i)
#                 all_zz.append(int_i)
#             if verbose:
#                 print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
#                       .format(zz.min(), zz.max(), np.sum(zz), float(np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
#
#             all_xrs = np.concatenate(all_xrs)
#             all_yrs = np.concatenate(all_yrs)
#             all_zz = np.concatenate(all_zz)
#             if int_or_hist == "int":
#
#                 return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
#             else:
#                 grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
#                 grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
#                 return (grid_x, grid_y, zz, all_xrs, all_yrs, all_zz)
#     else:
#         raise NotImplementedError("Not implemented")
def get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y", fac=1.0, nx=1000, ny=1000, extend=2):

    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)
    all_fluxes = np.concatenate(all_fluxes)

    # nx, ny = len(all_x), len(all_y)
    nx = np.complex(0, nx + 1)
    ny = np.complex(0, ny + 1)
    # edges_x = np.mgrid[all_x.min() * 1.0:all_x.max() * 1.0:nx]
    # edges_y = np.mgrid[all_y.min() * 1.0:all_y.max() * 1.0:ny]
    #
    # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])

    # grid_X, grid_Y = np.meshgrid(grid_x,grid_y,indexing="ij")

    # image = interp(all_x, all_y, all_fluxes, grid_X, grid_Y, method='linear')
    fac = 1.0
    grid_x, grid_y = np.mgrid[all_x.min()*extend:all_x.max()*extend:nx,
                     all_y.min()*extend:all_y.max()*extend:ny]

    # nx = np.complex(0, nx + 1)
    # ny = np.complex(0, ny + 1)
    # edges_x = np.mgrid[all_x.min() * extend:all_y.max() * extend:nx]
    # edges_y = np.mgrid[all_x.min() * extend:all_y.max() * extend:ny]

    # if (np.sum(all_x)==0. or np.sum(all_y)==0):
    if (np.sum(all_fluxes)==0.):
        raise ValueError(" x_arr or y_arr arrays in the image all FULL 0. Cannot re-interpolate!")
    if (all_x.min()==0 and all_x.max()==0):
        raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")
    if (all_y.min()==0 and all_y.max()==0):
        raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")
    try:
        image = interpolate.griddata(np.array([all_x, all_y]).T * fac, all_fluxes,
                                     (grid_x * fac, grid_y * fac), method='linear', fill_value=0)
    except:
        raise ValueError("failed interpolation.")
    latAvDist, latAvDist2, latMaxDist = lateral_distributions(grid_x, grid_y, image, collapse_axis=collapse_axis)
    if collapse_axis == "y":
        return (grid_x[:, 0], latAvDist, latAvDist2, latMaxDist)
    else:
        return (grid_y[0, :], latAvDist, latAvDist2, latMaxDist)

    # int_x, int_y, int_zz, all_x, all_y, all_fluxes \
    #     = pb.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True, int_or_hist="hist")
    if axis == "x":
        # nx = 200
        # ny = 100
        nx = np.complex(0, nx + 1)
        ny = np.complex(0, ny + 1)
        edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
        edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
        # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
        #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
        grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
        # tmp = []
        # for i in range(len(grid_x)):
        #     mask = i_zz[i, :] >  1e-3*i_zz[i, :].max()
        #     tmp.append(np.sum( (i_zz[i, :] * np.abs(np.diff(edges_y)))[mask] )  / np.sum(np.diff(edges_y)[mask]) )
        #     # tmp.append(np.average())
        # y_sum_ii = np.array(tmp)

        y_sum_ii = np.sum(i_zz * np.diff(edges_y), axis=1) / np.sum(np.diff(edges_y))
        max_ii = np.max(y_sum_ii)
        x1 = grid_x[np.argmin(y_sum_ii < max_ii * fac)]
        x2 = grid_x[::-1][np.argmin(y_sum_ii[::-1] < max_ii * fac)]
        # assert x2 >= x1
        return (x1, x2, grid_x, y_sum_ii)
    elif axis == "z":
        # nx = 200
        # ny = 100
        nx = np.complex(0, nx + 1)
        ny = np.complex(0, ny + 1)
        edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
        edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
        # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
        #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
        grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
        tmp = []
        for i in range(len(grid_y)):
            tmp.append(np.sum(i_zz[:, i] * np.diff(edges_x)) / np.sum(np.diff(edges_x)))
        tmp = np.array(tmp)

        x_sum_ii = np.sum(i_zz, axis=0)
        max_ii = np.max(x_sum_ii)
        y1 = grid_y[np.argmin(x_sum_ii < max_ii * fac)]
        y2 = grid_y[::-1][np.argmin(x_sum_ii[::-1] < max_ii * fac)]
        assert y2 > y1
        return (y1, y2, grid_y, x_sum_ii)
    else:
        raise KeyError()

#
#
# def _get_lc_for_freq(self, freq, times, freqs):
#     _times, _fluxes = [], []
#     for i in range(len(times)):
#         if (times[i] == freq):
#             _times.append(times[i])
#             _fluxes.append(freqs[i])
#     _times = np.array(_times)
#     _fluxes = np.array(_fluxes)
#     if (len(_times) < 1):
#         raise ValueError("Freq={} not found in data. Avaialable={}".format(freq, set(freqs)))
#     return (times, freqs)#(_times, _fluxes / 1e3 / 1e23)  # [s, mJy -> CGS]
# def _check_if_loaded_ej_lc(self):
#     if (self.ej_lc is None):
#         self.ej_lc = h5py.File(self.res_dir + self.ejecta_prefix + self.def_name_lc)
#
# def get_ej_lc_2d_arr(self, freq=None, ishell=None, ilayer=None):
#     dfile = self.get_ej_lc_obj()
#     v_n = list( dfile.keys() )[0]
#     times_ejecta = np.array(dfile["times"])
#     freqs_ejecta = np.array(dfile["freqs"])
#     fluxes_ejecta = np.array(dfile[v_n]["totalflux"])
#     if (freq is None and ishell is None and ilayer is None):
#         return fluxes_ejecta / 1e3 / 1e23
#     elif (ishell is None and ilayer is None):
#         return self._get_lc_for_freq(freq, times_ejecta, freqs_ejecta)
#
#     layer = "shell={} layer={}".format(ishell, ilayer)
#     if (not layer in list(dfile.keys())):
#         raise NameError("Layer {} (aka '{}') is not in the ejecta comov.spec. file.\n Available: {}"
#                         .format(ilayer, layer, dfile.keys()))
#     if (not v_n in dfile[layer].keys()):
#         raise NameError("v_n {} is not in the ejecta comov.spec. dfile[{}].keys() \n Avaialble: {}"
#                         .format(v_n, layer, dfile[layer].keys()))
#     if (freq is None):
#         return times_ejecta, fluxes_ejecta / 1e3 / 1e23  # [s, mJy -> CGS]
#     else:
#         _times, _fluxes = [], []
#         for i in range(len(times_ejecta)):
#             if (freqs_ejecta[i] == freq):
#                 _times.append(times_ejecta[i])
#                 _fluxes.append(fluxes_ejecta[i])
#         _times = np.array(_times)
#         _fluxes = np.array(_fluxes)
#         if (len(_times) < 1):
#             raise ValueError("Freq={} not found in data. Avaialable={}".format(freq, set(freqs_ejecta)))
#         return (_times, _fluxes / 1e3 / 1e23)  # [s, mJy -> CGS]
#
#     return np.array(dfile[layer][v_n])
def get_skymap_fwhm(grid, latAvDist, cc, fac=0.5):
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
    return (x1, x2)

def smooth_interpolated_skymap_with_gaussian_kernel(i_zz, type="uniform", sigma=3):
    if type=="uniform":
        i_zz = ndimage.uniform_filter(i_zz, size=sigma )
    elif type == "gaussian":
        if i_zz.ndim == 1:
            i_zz = ndimage.filters.gaussian_filter(i_zz, sigma=[sigma], mode='reflect')
        elif i_zz.ndim == 2:
            i_zz = ndimage.filters.gaussian_filter(i_zz, sigma=[sigma,sigma], mode='reflect')
    else:
        raise KeyError("smooth type is not recognize: {}".format(type))
    return i_zz





def get_skymap_spec_idx(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False):

    freqs = self.get_skymap_freqs()
    assert len(freqs) > 1
    idx = find_nearest_index(freqs, freq)
    dfile = self.get_skymap_obj()
    nshells = int(dfile.attrs["nshells"])

    all_xrs, all_yrs, all_zz = \
        self.get_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

    all_xrs_m1, all_yrs_m1, all_zz_m1 = \
        self.get_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

    assert len(all_zz) == len(all_zz_m1)

    all_values = []
    for i in range(len(all_zz_m1)):
        ffreq_i = np.full(all_zz[i].shape, freqs[idx])
        ffreq_im1 = np.full(all_zz_m1[i].shape, freqs[idx - 1])
        num = np.log10(all_zz[i]) - np.log10(all_zz_m1[i])
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        all_values.append(values)

    return (all_xrs, all_yrs, all_values)


