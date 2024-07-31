import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy
import shutil
import subprocess

from .utils import cgs, find_nearest_index
from .parfile_tools import read_parfile
from .OLD_skymap_tools import compute_position_of_the_flux_centroid
from .skymap_process import ProcessRawSkymap
from .id_analytic import JetStruct
from .parfile_tools import create_parfile

''' 
pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .
'''


class Base:
    def __init__(self,workingdir : str,readparfileforpaths : bool,parfile : str,verbose : bool):
        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir = workingdir

        self.fpath_id = None
        self.fpath_dyn = None
        self.fpath_spec = None
        self.fpath_light_curve = None
        self.fpath_sky_map = None

        self.id_dfile = None
        self.dyn_dfile = None
        self.spec_dfile = None
        self.lc_dfile = None
        self.skymap_dfile = None

        self.verb = verbose

    def read_grb_part_parfile(self,parfile="parfile.par"):
        grb_pars, grb_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                          sep1="# ---------------------- GRB afterglow ----------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_ejecta_id" in grb_opts.keys(): self.fpath_id = self.res_dir + grb_opts["fname_ejecta_id"]
        if "fname_dyn" in grb_opts.keys(): self.fpath_dyn = self.res_dir + grb_opts["fname_dyn"]
        if "fname_spectrum" in grb_opts.keys(): self.fpath_spec = self.res_dir + grb_opts["fname_spectrum"]
        if "fname_light_curve" in grb_opts.keys(): self.fpath_light_curve = self.res_dir + grb_opts["fname_light_curve"]
        if "fname_sky_map" in grb_opts.keys(): self.fpath_sky_map = self.res_dir + grb_opts["fname_sky_map"]
        return (grb_pars,grb_opts)
    def read_kn_part_parfile(self,parfile="parfile.par"):
        kn_pars, kn_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                        sep1="# ----------------------- kN afterglow ----------------------",
                                        sep2="# --------------------------- END ---------------------------")
        if "fname_ejecta_id" in kn_opts.keys(): self.fpath_id = self.res_dir + kn_opts["fname_ejecta_id"]
        if "fname_dyn" in kn_opts.keys(): self.fpath_dyn = self.res_dir + kn_opts["fname_dyn"]
        if "fname_spec" in kn_opts.keys(): self.fpath_spec = self.res_dir + kn_opts["fname_spec"]
        if "fname_light_curve" in kn_opts.keys(): self.fpath_light_curve = self.res_dir + kn_opts["fname_light_curve"]
        if "fname_sky_map" in kn_opts.keys(): self.fpath_sky_map = self.res_dir + kn_opts["fname_sky_map"]
        return (kn_pars,kn_opts)
    def read_pwn_part_parfile(self,parfile="parfile.par"):
        kn_pars, kn_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                        sep1="# --------------------------- PWN ---------------------------",
                                        sep2="# --------------------------- END ---------------------------")
        if "fname_ejecta_id" in kn_opts.keys(): self.fpath_id = self.res_dir + kn_opts["fname_ejecta_id"]
        if "fname_dyn" in kn_opts.keys(): self.fpath_dyn = self.res_dir + kn_opts["fname_dyn"]
        if "fname_spec" in kn_opts.keys(): self.fpath_spec = self.res_dir + kn_opts["fname_spec"]
        if "fname_light_curve" in kn_opts.keys(): self.fpath_light_curve = self.res_dir + kn_opts["fname_light_curve"]
        if "fname_sky_map" in kn_opts.keys(): self.fpath_sky_map = self.res_dir + kn_opts["fname_sky_map"]
        return (kn_pars,kn_opts)
    def clear(self):
        # self.overwrite = True
        if (not self.dyn_dfile is None):
            self.dyn_dfile.close()
            self.dyn_dfile = None
        if (not self.spec_dfile is None):
            self.spec_dfile.close()
            self.spec_dfile = None
        if (not self.lc_dfile is None):
            self.lc_dfile.close()
            self.lc_dfile = None
        if (not self.skymap_dfile is None):
            self.skymap_dfile.close()
            self.skymap_dfile = None

class Skymap:
    def __init__(self, dfile : h5py.Group):
        self.flux = dfile.attrs["flux"]
        self.xc = dfile.attrs["xc"]
        self.yc = dfile.attrs["yc"]
        self.x1 = dfile.attrs["x1"]
        self.x2 = dfile.attrs["x2"]
        self.y1 = dfile.attrs["y1"]
        self.y2 = dfile.attrs["y2"]
        self.grid_x = np.array(dfile["grid_x"])
        self.grid_y = np.array(dfile["grid_y"])
        self.dist_x = np.array(dfile["dist_x"])
        self.dist_y = np.array(dfile["dist_y"])
        self.im_intp = np.array(dfile["image_intp"])
        self.im_hist = np.array(dfile["image_hist"])
        self.time = dfile.attrs["time"]
        self.freq = dfile.attrs["freq"]


class Ejecta(Base):
    def __init__(self,workingdir : str, readparfileforpaths : bool, parfile : str, type : str, verbose : bool):
        super().__init__(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile,verbose=verbose)

        if not os.path.isdir(workingdir):
            raise IOError("Working directory not found {}".format(workingdir))
        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir = workingdir
        self.prefix = type+"_"

        if readparfileforpaths:
            self.reload_parfile(type=type)
        else:
            self.fpath_id = self.res_dir + self.prefix + "ejecta_id.h5"
            self.fpath_dyn = self.res_dir + self.prefix + "dynamics_layers.h5"
            self.fpath_spec = self.res_dir + self.prefix + "spectra.h5"
            self.fpath_light_curve = self.res_dir + self.prefix + "lightcurves_layers.h5"
            self.fpath_sky_map = self.res_dir + self.prefix + "skymap.h5"

    def reload_parfile(self, type : str) -> None:
        if (type=="kn"):
            self.pars, self.opts = self.read_kn_part_parfile( self.parfile )
        elif(type=="grb"):
            self.pars, self.opts = self.read_grb_part_parfile( self.parfile )
        elif(type=="pwn"):
            self.pars, self.opts = self.read_pwn_part_parfile( self.parfile )
        else:
            raise KeyError("not implemented")

            # ejecta id file
    def _ckeck_if_loaded_id_obj(self):
        if (self.fpath_id is None):
            raise IOError("self.fpath_id is not set")
        if (self.id_dfile is None):
            self.id_dfile = h5py.File(self.fpath_id, "r")
    def _ckeck_if_loaded_dyn_obj(self):
        if (self.fpath_dyn is None):
            raise IOError("self.fpath_kn_dyn is not set")
        if (self.dyn_dfile is None):
            self.dyn_dfile = h5py.File(self.fpath_dyn)

    # ejecta spectrum
    def _check_if_loaded_spec(self):
        if (self.fpath_spec is None):
            raise IOError("self.fpath_spec is not set")
        if (self.spec_dfile is None):
            self.spec_dfile = h5py.File(self.fpath_spec)
    # ejecta lightcurves
    def _check_if_loaded_lc(self):
        if (self.fpath_light_curve is None):
            raise IOError("self.fpath_light_curve is not set")
        if (self.lc_dfile is None):
            self.lc_dfile = h5py.File(self.fpath_light_curve)
    # ejecta skymaps
    def _check_if_loaded_skymap(self):
        if (self.fpath_sky_map is None):
            raise IOError("self.fpath_kn_sky_map is not set")
        if (self.skymap_dfile is None):
            try:
                self.skymap_dfile = h5py.File(self.fpath_sky_map)
            except OSError:
                raise OSError("failed to load the file: \n {}".format(self.fpath_sky_map))

    def _get_1d_or_2d_array(self, dfile : h5py.File, ishell : int, ilayer : int, v_n : str) -> np.ndarray:
        nlayers = int(dfile.attrs["nlayers"])
        nshells = int(dfile.attrs["nshells"])
        if ((ishell is None) and (ilayer is None)):
            raise ValueError("both ishell and ilayer cannot be None for this data")
        elif ((not ishell is None) and (ilayer is None)):
            arr = []
            for il in range(nlayers):
                layer = f"shell={ishell} layer={il}"
                # arr.append(self.get_dyn_arr(v_n, ishell=ishell, ilayer=il))
                arr = np.array(dfile[layer][v_n])
            arr = np.reshape(np.array(arr), newshape=(nlayers, len(arr[0])))
            if (np.sum(arr) == 0):
                print(f"Warning np.sum(arr)=0 for ishell={ishell} ilayer={ilayer}")
            return arr
        elif ((ishell is None) and (not ilayer is None)):
            arr = []
            for ish in range(nshells):
                layer = f"shell={ish} layer={ilayer}"
                # arr.append(self.get_dyn_arr(v_n, ishell=ish, ilayer=ilayer))
                arr = np.array(dfile[layer][v_n])
            arr = np.reshape(np.array(arr), newshape=(nshells, len(arr[0])))
            if (np.sum(arr) == 0):
                print(f"Warning np.sum(arr)=0 for ishell={ishell} ilayer={ilayer}")
            return arr
        elif ((not ishell is None) and (not ilayer is None)):
            # group = dfile["shell={} layer={}".format(ishell, ilayer)]
            layer = f"shell={ishell} layer={ilayer}"
            if (not layer in list(dfile.keys())):
                raise NameError(f"Layer {ilayer} (key '{layer}') is not in the ejecta file nlayer={nlayers} nshells={nshells}")
            if not (v_n in dfile[layer].keys()):
                raise KeyError(f"key {v_n} is not in the list: {dfile[layer].keys()}")
            arr = np.array(dfile[layer][v_n])
            if (np.sum(arr) == 0):
                print(f"Warning np.sum(arr)=0 for ishell={ishell} ilayer={ilayer}")
            return arr
        else:
            raise NameError()

    # --------- Initial Data -----------

    def get_id_obj(self):
        self._ckeck_if_loaded_id_obj()
        return self.id_dfile

    def get_id_attr(self, v_n : str) -> str or np.float64:
        dfile = self.get_id_obj()
        if (not v_n in dfile.attrs.keys()):
            raise KeyError(f"key={v_n} is not in the list={dfile.attrs.keys()}")
        return dfile.attrs[v_n]

    # --------- Dynamics -------------

    def get_dyn_obj(self):
        self._ckeck_if_loaded_dyn_obj()
        return self.dyn_dfile

    def get_dyn_arr(self, v_n : str, ishell : int or None = None, ilayer : int or None = None) -> np.ndarray:
        self._ckeck_if_loaded_dyn_obj()
        return self._get_1d_or_2d_array(self.get_dyn_obj(),v_n=v_n,ishell=ishell,ilayer=ilayer)

    # --------- Light Curves & Spectra ---------

    def get_lc_obj(self, spec : bool = False) -> h5py.File:
        if spec:
            self._check_if_loaded_spec()
            return self.spec_dfile
        else:
            self._check_if_loaded_lc()
            return self.lc_dfile

    def get_grid(self, key:str, unique : bool = True, spec : bool = False) -> np.ndarray:
        dfile = self.get_lc_obj(spec=spec)
        if not key in dfile.keys():
            raise KeyError(f"Key={key} is not found in dfile.keys():\n"
                           f"{dfile.keys()}")
        arr = np.array(dfile[key])
        if (not unique):
            return arr
        arr_u = np.unique(arr)
        if len(arr_u) == 0:
            raise ValueError(f"no unique values found in data={key} spec={spec} \n {arr}")
        return np.array(arr_u)

    def get_gams(self,unique : bool = True) -> np.ndarray:
        return self.get_grid(key="gams",unique=unique,spec=True)
    #     dfile = self.get_lc_obj(spec=True)
    #     arr = np.array(dfile["gams"])
    #     if (not unique):
    #         return arr
    #     arr_u = np.unique(arr)
    #     if len(arr_u) == 0:
    #         raise ValueError("no unique gams found in data \n {}".format(arr))
    #     return np.array(arr_u)# np.array(dfile["freqs"])
    #
    def get_lc_times(self, unique : bool = True, spec : bool = False) -> np.ndarray:
        return self.get_grid(key="times",unique=unique, spec=spec)
    #     dfile = self.get_lc_obj(spec=spec)
    #     arr = np.array(dfile["times"])
    #     if (not unique):
    #         return arr
    #     arr_u = np.unique(arr)
    #     if len(arr_u) == 0:
    #         raise ValueError("no unique times found in array \n {}".format(arr))
    #     return arr_u
    #
    def get_lc_freqs(self,unique : bool = True, spec : bool = False) -> np.ndarray:
        return self.get_grid(key="freqs",unique=unique, spec=spec)
    #     dfile = self.get_lc_obj(spec=spec)
    #     arr = np.array(dfile["freqs"])
    #     if (not unique):
    #         return arr
    #     arr_u = np.unique(arr)
    #     if len(arr_u) == 0:
    #         raise ValueError("no unique freqs found in light curve \n {}".format(arr))
    #     return np.array(arr_u)# np.array(dfile["freqs"])

    def _get_closest_grid_val(self, val, key:str, uvals:np.ndarray)->float:
        if (not val in uvals):
            if (val > uvals.max()):
                raise ValueError(f"requested {key}={val} > dfile unique {key}.max()={uvals.max()}")
            if (val < uvals.min()):
                raise ValueError(f"requested {key}={val} < dfile unique {key}.min()={uvals.min()}")
            _val = uvals[find_nearest_index(uvals, val)]
            if self.verb:
                print(f"Warning: {val}={val} is not in unique {val} {uvals} Using {val}={_val}")
        else:
            _val = uvals[int(np.where(uvals==val)[0])]
        return _val

    def OLD_get_lc_totalflux(self, freq : float or None = None, time : float or None = None, spec : bool = False) -> np.ndarray:
        dfile = self.get_lc_obj(spec=spec)
        utimes = self.get_lc_times(spec=spec,unique=True)
        ufreqs = self.get_lc_freqs(spec=spec,unique=True)
        fluxes =  np.array(dfile["total_power"]) if spec else np.array(dfile["total_fluxes"])
        # return np.zeros(0,)
        # def _get_time() -> float:
        #     if (not time in utimes):
        #         if (time > utimes.max()):
        #             raise ValueError(f"requested time={time} > dfile times.max()={utimes.max()}")
        #         if (time < utimes.min()):
        #             raise ValueError(f"requested time={time} < dfile times.min()={utimes.min()}")
        #         _time = utimes[find_nearest_index(utimes, time)]
        #         if self.verb:
        #             print(f"Warning: time={time} is not in {utimes} Using time={_time}")
        #     else:
        #         _time = utimes[int(np.where(utimes==time)[0])]
        #     return _time
        #
        # def _get_freq() -> float:
        #     if (not freq in ufreqs):
        #         if (freq > ufreqs.max()):
        #             raise ValueError(f"requested freq={freq} > dfile freqs.max()={ufreqs.max()}")
        #         if (freq < ufreqs.min()):
        #             raise ValueError(f"requested freq={freq} < dfile freqs.min()={ufreqs.min()}")
        #         _freq = ufreqs[find_nearest_index(ufreqs, freq)]
        #         if self.verb: print(f"Warning: freq={freq} is not in {ufreqs} Using freq={_freq}")
        #     else:
        #         _freq = ufreqs[int(np.where(ufreqs==freq)[0])]
        #     return _freq

        if ((freq is None) and (not time is None)):
            # light curve mode
            _time = self._get_closest_grid_val(val=time, key="time", uvals=utimes)
            arr = fluxes[np.where(self.get_lc_times(spec=spec,unique=False) == _time)]
            return arr
        elif ((not freq is None) and (time is None)):
            # spectrum mode
            _freq = self._get_closest_grid_val(val=freq, key="freq", uvals=ufreqs)
            arr = fluxes[np.where(self.get_lc_freqs(spec=spec,unique=False) == _freq)]
            return arr
        elif ((not time is None) and (not freq is None)):
            _time = self._get_closest_grid_val(val=time, key="time", uvals=utimes)
            _freq = self._get_closest_grid_val(val=freq, key="freq", uvals=ufreqs)
            # flux at a time and freq mocde
            arr = fluxes[np.where(((self.get_lc_freqs(spec=spec,unique=False) == _freq).astype(int) *
                                   (self.get_lc_times(spec=spec,unique=False) == _time).astype(int)).astype(bool))]
            return arr
        else:
            arr = np.vstack(( [fluxes[np.where(self.get_lc_freqs(spec=spec,unique=False)==_freq)] for _freq in ufreqs] ))
            return arr


        # if (not time is None):
        #     if (not time in utimes):
        #         _time = utimes[find_nearest_index(utimes, time)]
        #         if self.verb: print(f"Warning: time={time} is not in {utimes} Using time={_time}")
        #     else:
        #         _time = utimes[int(np.where(utimes==time)[0])]
        #
        #
        # if (freq is None):
        #     arr = np.vstack(( [fluxes[np.where(self.get_lc_freqs(spec=spec,unique=False)==_freq)] for _freq in ufreqs] ))
        #     if (not time is None):
        #         raise NotImplementedError("method is not implemented")
        #     return arr
        # else :
        #     if (freq > ufreqs.max()):
        #         raise ValueError(f"requested freq={freq} > dfile freqs.max()={ufreqs.max()}")
        #     if (freq < ufreqs.min()):
        #         raise ValueError(f"requested freq={freq} < dfile freqs.min()={ufreqs.min()}")
        #     if (not freq in ufreqs):
        #         _freq = ufreqs[find_nearest_index(ufreqs, freq)]
        #         if self.verb: print(f"Warning: freq={freq} is not in {ufreqs} Using freq={_freq}")
        #     else:
        #         _freq = ufreqs[int(np.where(ufreqs==freq)[0])]
        #
        #     if (not time is None):
        #         _i = self.get_lc_times(spec=spec,unique=False) == _time
        #         _j = self.get_lc_freqs(spec=spec,unique=False) == _freq
        #         arr = fluxes[np.where(((self.get_lc_freqs(spec=spec,unique=False) == _freq).astype(int) *
        #                                (self.get_lc_times(spec=spec,unique=False) == _time).astype(int)).astype(bool))]
        #     else:
        #         arr = fluxes[np.where(self.get_lc_freqs(spec=spec,unique=False) == _freq)]
        #     return arr


        # if (freq is None):
        #     return fluxes
        # if (not freq in freqs):
        #     raise ValueError("freq={} not found in freqs={}".format(freq, freqs))
        # arr = fluxes[freqs==freq]
        # if (len(arr)==0):
        #     raise ValueError("no fluxes found for freq={}".format(freq))
        # return arr

        #
        # nlayers = int(dfile.attrs["nlayers"])
        # nshells = int(dfile.attrs["nshells"])
        # times = self.get_ej_lc_times()
        # freqs = self.get_ej_lc_freqs()
        # try:
        #     key = str("totalflux at freq={:.4e}".format(freq)).replace('.', ',')
        #     if not key in dfile.keys():
        #         raise NameError("Not found: {} among keys:{}".format(key, [key for key in dfile.keys() if key.__contains__("totalflux")]))
        # except NameError:
        #     key = str("totalflux at freq={:.4e}".format(freq))
        #     if not key in dfile.keys():
        #         raise NameError("Not found for ej. lightcurve: {} among keys:{} \n ejecta_prefix:{}"
        #                         .format(key, [key for key in dfile.keys() if key.__contains__("totalflux")], self.ejecta_prefix))
        # except:
        #     raise NameError()
        # return np.array(dfile[key])

    def _get_lc_shell_layer(self, key:str, ishell:int, ilayer:int, spec=False)->np.ndarray:

        dfile = self.get_lc_obj(spec=spec)
        grp = dfile["shell={} layer={}".format(ishell, ilayer)]
        if not key in grp.keys():
            raise KeyError(f"key={key} is not in the dfile: {grp.keys()}")
        arr = np.array(grp[key])
        return arr

    def _get_lc_for_mask(self, key:str, mask:np.ndarray or None,
                         x_arr:np.ndarray,
                         ishell:int or None, ilayer:int or None, spec=False, sum_shells_layers=False):

        dfile = self.get_lc_obj(spec=spec)
        if not ("nlayers" in dfile.attrs.keys()):
            raise KeyError(f"key= nlayers is not in dfile.attrs.keys()=[{dfile.attrs.keys()}]")
        nlayers = int(dfile.attrs["nlayers"])
        if not ("nshells" in dfile.attrs.keys()):
            raise KeyError(f"key= nshells is not in dfile.attrs.keys()=[{dfile.attrs.keys()}]")
        nshells = int(dfile.attrs["nshells"])
        # compute the array of a required shape
        if not ishell is None and not ilayer is None:
            res = self._get_lc_shell_layer(key=key,ishell=ishell, ilayer=ilayer, spec=spec)
            if not mask is None: res = res[mask]
            return res
        elif ishell is None and not ilayer is None:
            res = [self._get_lc_shell_layer(key=key,ishell=i, ilayer=ilayer, spec=spec) for i in range(nshells)]
            if not mask is None: res = [res_i[mask] for res_i in res]
            res = np.reshape(np.array(res),newshape=(nshells,len(x_arr)))
            if sum_shells_layers: res = np.sum(res,axis=(0))
            return res
        elif not ishell is None and ilayer is None:
            res = [self._get_lc_shell_layer(key=key,ishell=ishell, ilayer=i, spec=spec) for i in range(nlayers)]
            if not mask is None: res = [res_i[mask] for res_i in res]
            res = np.reshape(np.array(res),newshape=(nlayers, len(x_arr)))
            if sum_shells_layers: res = np.sum(res,axis=(0))
            return res
        else:

            res = [[self._get_lc_shell_layer(key=key,ishell=i, ilayer=j, spec=spec)
                   for i in range(nshells)] for j in range(nlayers)]
            if not mask is None: res = [[res_i[mask] for res_i in res_j] for res_j in res]
            if sum_shells_layers:
                res_ = np.zeros_like(res[0][0])
                for ishell in range(nshells):
                    for ilayer in range(nlayers):
                        res_ += res[ilayer][ishell]
            else:
                res_= np.reshape(np.array(res),newshape=(nshells,nlayers,len(x_arr)))
            return np.array(res_)

    def get_lc(self,
               key:str="fluxdens", xkey:str="freqs", key_time:str= "times",
               freq:float or None = None, time:float or None=None,
               ishell:int or None=None, ilayer:int or None=None, sum_shells_layers:bool=True,
               spec:bool=False):
        """

        :param key: str. Options(spec=True): n_ele, synch, ssa, ssc. Options (spec=False): fluxdens
        :param xkey: str. Options(spec=True): gams, freqs. Options (spec=False): freqs
        :param key_time: str. Options(spec=True): times_gams, times_freqs. Options (spec=False): times
        :param freq: float or None (frequency) [hz]
        :param time: float or None (time) [s]
        :param ishell: int or None : number of the shell (velocity structure)
        :param ilayer: int or None : number of the layer (angualr structure)
        :param sum_shells_layers: bool (if ishell or ilayer = None, sum the spectra/lcs over them)
        :param spec: bool (use spectral dfile or lc dfile)
        :return: np.ndarray
        """
        dfile = self.get_lc_obj(spec=spec)
        if not xkey in dfile.keys():
            raise KeyError(f"xkey={xkey} is not recognized. Avaialble: {dfile.keys()}")

        utimes = self.get_grid(key=key_time, unique=True, spec=spec)
        ufreqs = self.get_grid(key=xkey, unique=True, spec=spec)

        # light curve
        if (time is None) and (not freq is None):
            # get freq that is in the h5file
            _freq = self._get_closest_grid_val(val=freq, key=xkey, uvals=ufreqs)
            # get mask for this freq
            mask = np.array(self.get_grid(key=xkey,spec=spec,unique=False) == _freq,dtype=bool)
            # get light curve for this mask and for this shell and layer
            res = self._get_lc_for_mask(key=key,mask=mask, x_arr=utimes, ishell=ishell, ilayer=ilayer, spec=spec,
                                        sum_shells_layers=sum_shells_layers)
            return res

        # spectrum
        if (not time is None) and (freq is None):
            # get time that is in the h5file
            _time = self._get_closest_grid_val(val=time, key=key_time, uvals=utimes)
            # get mask for this fre,q
            mask = np.array(self.get_grid(key=key_time, spec=spec, unique=False) == _time, dtype=bool)
            # get spectrum for this mask and for this shell and layer
            res = self._get_lc_for_mask(key=key, mask=mask, x_arr=ufreqs, ishell=ishell, ilayer=ilayer, spec=spec,
                                        sum_shells_layers=sum_shells_layers)
            return res

        if (time is None and freq is None) and (not ishell is None) and (not ilayer is None):
            if xkey=="gams":
                res = self._get_lc_shell_layer(key=key, ishell=ishell, ilayer=ilayer, spec=spec)
                res = np.reshape(res, newshape=(len(utimes),len(ufreqs))).T
            else:
                res = self._get_lc_shell_layer(key=key,ishell=ishell, ilayer=ilayer, spec=spec)
                res = np.reshape(res, newshape=(len(utimes),len(ufreqs))).T
                if sum_shells_layers: res = np.sum(res,axis=(0,1))
            return res

        if (time is None and freq is None) and (ishell is None) and (ilayer is None):
            res = self._get_lc_for_mask(key=key, mask=None, x_arr=ufreqs, ishell=ishell, ilayer=ilayer, spec=spec,
                                        sum_shells_layers=sum_shells_layers)
            res = np.reshape(res, (len(utimes),len(ufreqs))).T
            return res

        # total spectrum from all shells/layers
        # if (time is None and freq is None) and (ishell is None) and (ilayer is None):
        #     res = self._get_lc_for_mask()


        raise ValueError("Either time of freq must be specified to compute lc/spec")

    def OLD_get_lc(self, freq=None, time=None, ishell=None, ilayer=None, spec=False) -> np.ndarray:
        dfile = self.get_lc_obj(spec=spec)
        if not ("nlayers" in dfile.attrs.keys()):
            raise KeyError(f"key= nlayers is not in dfile.attrs.keys()=[{dfile.attrs.keys()}]")
        nlayers = int(dfile.attrs["nlayers"])
        if not ("nshells" in dfile.attrs.keys()):
            raise KeyError(f"key= nshells is not in dfile.attrs.keys()=[{dfile.attrs.keys()}]")
        nshells = int(dfile.attrs["nshells"])
        utimes = self.get_lc_times(spec=spec,unique=True)
        ufreqs = self.get_lc_freqs(spec=spec,unique=True)


        tidx = None
        if (not time is None):
            if (time > utimes.max()):
                raise ValueError(f"requested time={time} > dfile times.max()={utimes.max()}")
            if (time < utimes.min()):
                raise ValueError(f"requested time={time} < dfile times.min()={utimes.min()}")
            if (not time in utimes):
                _time = utimes[find_nearest_index(utimes, time)]
                if self.verb: print(f"Warning: time={time} is not in {utimes} Using time={_time}")
            else:
                _time = utimes[int(np.where(utimes==time)[0])]
        else:
            _time = time
        if (not freq is None):
            if (freq > ufreqs.max()):
                raise ValueError(f"requested freq={freq} > dfile freqs.max()={ufreqs.max()}")
            if (freq < ufreqs.min()):
                raise ValueError(f"requested freq={freq} < dfile freqs.min()={ufreqs.min()}")
            if (not freq in ufreqs):
                _freq = ufreqs[find_nearest_index(ufreqs, freq)]
                if self.verb: print(f"Warning: freq={freq} is not in {ufreqs} Using freq={_freq}")
            else:
                _freq = ufreqs[int(np.where(ufreqs==freq)[0])]
        else:
            _freq = freq
        # if (freq is None and time is None):
        #     raise KeyError("Only one, time or freq can be none (for 2D output")

        times = self.get_lc_times(spec=spec,unique=False)
        freqs = self.get_lc_freqs(spec=spec,unique=False)



        if (freq is None):
            # spectum
            if ((ishell is None) and (ilayer is None)):
                fluxes2d = []
                for ifreq in ufreqs:
                    fluxes2d.append(self.get_lc_totalflux(freq=ifreq,time=_time,spec=spec))  # [freq,time]
                fluxes2d = np.reshape(np.array(fluxes2d), (len(freqs), len(times)))
                if (not time is None): return fluxes2d[tidx,:]
                else: return fluxes2d
            elif ((ishell is None) and (not ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                for ish in range(nshells):
                    arr = np.array(dfile["shell={} layer={}".format(ish, ilayer)])
                    arr = np.reshape(arr, (len(times),len(freqs)))
                    fluxes2d += arr  # [freq,time]
                if (not time is None): return fluxes2d[tidx,:]
                else: return fluxes2d
            elif ((not ishell is None) and (ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                for il in range(nlayers):
                    arr = np.array(dfile["shell={} layer={}".format(ishell, il)])
                    arr = np.reshape(arr, (len(times),len(freqs)))
                    fluxes2d += arr  # [freq,time]
                if (not time is None): return fluxes2d[tidx,:]
                else: return fluxes2d
            elif ((not ishell is None) and (not ilayer is None)):
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                arr = np.array(dfile["shell={} layer={}".format(ishell, ilayer)])
                if (time is None):
                    data = np.vstack((
                        [self.get_lc(freq=_freq,ishell=ishell,ilayer=ilayer,spec=spec)[1:-1]
                         for _freq in self.get_lc_freqs(spec=spec,unique=True)]
                    ))
                    return data

                arr = arr[np.where(times==_time)]
                return arr
                # arr = np.reshape(arr, (len(utimes),len(ufreqs)))
                # fluxes2d = arr  # [freq,time]
                # if (not time is None): return arr[tidx,:]
                # else: return fluxes2d
            else:
                raise NameError()
        else:
            # light curves
            if (not _freq in ufreqs):
                raise ValueError("freq:{} is not in ej_lc Given:{}".format(_freq, ufreqs))
            # ifreq = find_nearest_index(self.get_ej_lc_freqs(), freq)
            if ((ishell is None) and (ilayer is None)):
                return self.get_lc_totalflux(freq=_freq,time=_time,spec=spec)
            elif ((ishell is None) and (not ilayer is None)):
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                if (ilayer > nlayers - 1):
                    raise ValueError("Layer={} is > nlayers={}".format(ilayer, nlayers - 1))
                # fluxes1d = np.zeros_like(times)
                # for ish in range(nshells):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ish, ilayer)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for ish in range(nshells):
                    arr = np.array(dfile["shell={} layer={}".format(ish, ilayer)])
                    arr = arr[freqs==_freq]
                    fluxes2d.append(arr)  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nshells, len(times)))
            elif ((not ishell is None) and (ilayer is None)):
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                if (ishell > nshells - 1):
                    raise ValueError("Shell={} is > nshell={}".format(ishell, nshells - 1))
                # fluxes1d = np.zeros_like(times)
                # for il in range(nlayers):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ishell, il)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for il in range(nlayers):
                    arr = np.array(dfile["shell={} layer={}".format(ishell, il)])
                    arr = arr[freqs==_freq]
                    fluxes2d.append(arr)  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nlayers, len(times)))
            elif ((not ishell is None) and (not ilayer is None)):
                if self.verb: print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                arr = np.array(dfile["shell={} layer={}".format(ishell, ilayer)])
                arr = arr[np.where(freqs==_freq)]
                fluxes1d = arr  # [freq,time]
                if (len(fluxes1d)!=len(utimes)):
                    raise ValueError("size mismatch")
                return fluxes1d
            else:
                raise NameError()


    def _alpha(self, freq, freqm1, ishell, ilayer, spec=False) -> np.ndarray:
        values_i = self.get_lc(freq=freq, ishell=ishell, ilayer=ilayer,spec=spec)
        values_im1 = self.get_lc(freq=freqm1, ishell=ishell, ilayer=ilayer,spec=spec)
        ffreq_i = np.full(values_i.shape, freq)
        ffreq_im1 = np.full(values_im1.shape, freqm1)
        num = np.log10(values_i) - np.log10(values_im1)
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        return values

    def get_lc_spec_idx(self, freq, ishell=None, ilayer=None, freq1=None, freq2=None, spec=False):
        freqs = self.get_lc_freqs(spec=spec)
        times = self.get_lc_times(spec=spec)
        if ((freq is None) and (ishell is None) and (ilayer is None)):
            if (not freq1 is None) or (not freq2 is None):
                raise KeyError("For 2d spectral index freq1 and freq2 should be None (for now)")
            arr = []
            for ifreq in range(1, len(freqs)):
                arr.append( self._alpha(freq=freqs[ifreq], freqm1=freqs[ifreq - 1], ishell=ishell, ilayer=ilayer, spec=spec) )
            return np.reshape(np.array(arr), newshape=(len(freqs) - 1, len(times)))
        else:
            idx = find_nearest_index(freqs, freq)
            _freq1 = freqs[idx] if freq1 is None else freq1
            _freq2 = freqs[idx - 1] if freq2 is None else freq2
            return self._alpha(freq=_freq1, freqm1=_freq2, ishell=ishell, ilayer=ilayer, spec=spec)

    def get_lc_temp_idx(self, freq, ishell, ilayer, spec=False):
        freqs = self.get_lc_freqs(spec=spec)
        t = self.get_lc_times(spec=spec)
        Lnu = self.get_lc(freq=freq, ishell=ishell, ilayer=ilayer,spec=spec)
        val = np.zeros_like(Lnu)
        val[1:-1] = np.log10(Lnu[2:] / Lnu[:-2]) / np.log10(t[2:] / t[:-2])
        val[0] = np.log10(Lnu[1] / Lnu[0]) / np.log10(t[1] / t[0])
        val[-1] = np.log10(Lnu[-1] / Lnu[-2]) / np.log10(t[-1] / t[-2])
        return val

    # ------------ spectrum -------------
    #
    # def get_spec_obj(self):
    #     return self.get_lc_obj(spec=True)
    #     # self._check_if_loaded_spec()
    #     # return self.spec_dfile
    #
    # def get_spec_times)
    #
    # def get_spec_2d_arr(self, v_n="em", ishell=0, ilayer=0):
    #     dfile = self.get_spec_obj()
    #     layer = "shell={} layer={}".format(ishell, ilayer)
    #     if (not layer in list(dfile.keys())):
    #         raise NameError("Layer {} (aka '{}') is not in the ejecta comov.spec. file.\n Available: {}"
    #                         .format(ilayer, layer, dfile.keys()))
    #     if (not v_n in dfile[layer].keys()):
    #         raise NameError("v_n {} is not in the ejecta comov.spec. dfile[{}].keys() \n Avaialble: {}"
    #                         .format(v_n, layer, dfile[layer].keys()))
    #     return np.array(dfile[layer][v_n])
    #
    # def get_spec_times(self):
    #     dfile = self.get_spec_obj()
    #     return np.array(dfile["times"])
    #
    # def get_spec_freqs(self):
    #     dfile = self.get_spec_obj()
    #     return np.array(dfile["freqs"])
    #
    # def get_spec_1d_arr(self, freq=None, time=None, v_n="em", ishell=0, ilayer=0):
    #     arr = self.get_spec_2d_arr(v_n=v_n, ishell=ishell, ilayer=ilayer)
    #     freqs = self.get_spec_freqs()
    #     times = self.get_spec_times()
    #     if (not freq is None):
    #         if (freq < freqs.min() or freq > freqs.max()):
    #             raise ValueError("Freq={} (jet comov cpec) is not in the limit of freqs[{}, {}]"
    #                              .format(freq, freqs.min(), freqs.max()))
    #         if (not freq in freqs):
    #             idx = find_nearest_index(freqs, freq)
    #         else:
    #             idx = int(np.where(freqs == freq)[0])
    #         arr_ = arr[idx, :]
    #         assert len(arr_) == len(times)
    #         return arr_
    #     if (not times is None):
    #         if (time < times.min() or time > times.max()):
    #             raise ValueError("Time={} (jet comov cpec) is not in the limit of times[{}, {}]"
    #                              .format(time, times.min(), times.max()))
    #         if (not time in times):
    #             idx = find_nearest_index(time, times)
    #         else:
    #             idx = int(np.where(time == times)[0][0])
    #         arr_ = arr[:, idx]
    #         assert len(arr_) == len(freqs)
    #         return arr_
    #
    # def get_spec_2d_arr_layers(self, freq=None, time=None, v_n="em", ishell=None, ilayer=None):
    #     obj_dyn = self.get_dyn_obj()
    #     arr = []
    #     x_arr = []
    #     # x_arr = self.get_ej_dyn_arr("ctheta",ishell=ishell,ilayer=0)
    #     if (not ishell is None and ilayer is None):
    #         y_arr = self.get_dyn_arr("R", ishell=ishell, ilayer=0)
    #         # freqs = self.get_ej_spec_freqs()
    #         # times = self.get_ej_spec_times()
    #         nlayers = int(obj_dyn.attrs["nlayers"])
    #         for il in range(nlayers):
    #             x_arr.append(self.get_dyn_arr("ctheta", ishell=ishell, ilayer=il)[0])
    #             arr.append(self.get_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ishell, ilayer=il))
    #         x_arr = np.array(x_arr)
    #         arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
    #     else:
    #         y_arr = self.get_dyn_arr("R", ishell=0, ilayer=ilayer)
    #         nshells = int(obj_dyn.attrs["nshells"])
    #         for ish in range(nshells):
    #             x_arr.append(self.get_dyn_arr("Gamma", ishell=ish, ilayer=ilayer)[0])
    #             arr.append(self.get_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ish, ilayer=ilayer))
    #         x_arr = np.array(x_arr)
    #         arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
    #     return arr
    pass
    # ---------- Sky Maps --------------

    def get_skymap_obj(self) -> h5py.File:
        self._check_if_loaded_skymap()
        return self.skymap_dfile

    def get_skymap_times(self, unique:bool=True) -> np.ndarray:
        self._check_if_loaded_skymap()
        # print(self.ej_skymap.keys())
        res = np.array(self.skymap_dfile["times"])
        if unique: return np.unique(res)
        else: return res

    def get_skymap_freqs(self,unique:bool=True) -> np.ndarray:
        self._check_if_loaded_skymap()
        res = np.array(self.skymap_dfile["freqs"])
        if unique: return np.unique(res)
        else: return res

    def check_skymap_time(self, time : float) -> float:
        times = self.get_skymap_times()
        if not time in times:
            if (time < self.get_skymap_times().min()):
                raise ValueError(f"time {time} < times.min()={times.min()}")
            if (time > self.get_skymap_times().max()):
                raise ValueError(f"time {time} > times.max()={times.max()}")
                # self.get_skymap_totfluxes(freq=freq, shell=None, time=None)
            time = times[find_nearest_index(times, time)]
        return time

    def check_skymap_freq(self, freq : float) -> float:
        freqs = self.get_skymap_freqs()
        if not freq in freqs:
            if (freq < self.get_skymap_times().min()):
                raise ValueError(f"freq {freq} < freqs.min()={freqs.min()}")
            if (freq > self.get_skymap_times().max()):
                raise ValueError(f"freq {freq} > freqs.max()={freqs.max()}")
                # self.get_skymap_totfluxes(freq=freq, shell=None, time=None)
            freq = freqs[find_nearest_index(freqs, freq)]
        return freq

    def get_skymap_attr(self, v_n : str, freq : float, time : float) -> np.number:

        self.check_skymap_time(time=time)
        self.check_skymap_freq(freq=freq)

        dfile = self.get_skymap_obj()
        ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]

        if (not v_n in ddfile.attrs.keys()):
            raise KeyError(f"key '{v_n}' is not found in \n {ddfile.attrs.keys()}")

        return ddfile.attrs[v_n]

    def get_skymap(self, time : float, freq : float) -> Skymap:

        time = self.check_skymap_time(time=time)
        freq = self.check_skymap_freq(freq=freq)
        # if (not type in ["hist","intp"]):
        #     raise KeyError(f"skymap type '{type}' is not recognized.")
        dfile = self.get_skymap_obj()
        key = "time={:.4e} freq={:.4e}".format(time, freq)
        if not key in dfile.keys():
            raise KeyError(f"key={key} not in keys for skymap: {dfile.keys()}")
        return Skymap(dfile[key])
        # dfile = self.get_skymap_obj()
        # ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
        #
        # skymap = np.array(ddfile["image_"+type], dtype=np.float64)
        # grid_x = np.array(ddfile["grid_x"], dtype=np.float64)
        # grid_y = np.array(ddfile["grid_y"], dtype=np.float64)
        #
        # return (grid_x, grid_y, skymap)


class Magnetar:
    def __init__(self,workingdir:str,readparfileforpaths:bool,parfile:str,verbose:bool):
        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir = workingdir
        self.fpath_mag = None
        self.verb = verbose
        if readparfileforpaths:
            self.mag_pars, self.mag_opts = self.read_magnetar_part_parfile( self.parfile )
        else:
            self.fpath_mag = self.res_dir + "magnetar.h5"

        self.mag_dfile = None

    def reload_parfile(self):
        self.mag_pars, self.mag_opts = self.read_magnetar_part_parfile( self.parfile )

    def read_magnetar_part_parfile(self, parfile="parfile.par"):
        mag_pars, mag_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                          sep1="# ------------------------ Magnetar -------------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_mag" in mag_opts.keys(): self.fpath_mag = self.res_dir + mag_opts["fname_mag"]
        return (mag_pars, mag_opts)

    def _check_if_loaded_mag_obj(self):
        if (self.fpath_mag is None):
            raise IOError("self.fpath_mag is not set")
        if (self.mag_dfile is None):
            self.mag_dfile = h5py.File(self.fpath_mag)
    # magnetar
    def get_mag_obj(self):
        self._check_if_loaded_mag_obj()
        return self.mag_dfile

    def clear(self):
        # self.overwrite = True
        if (not self.mag_dfile is None):
            self.mag_dfile.close()
            self.mag_dfile = None


class PyBlastAfterglow:
    '''
        Process output_uniform_grb files: load, extract for a specific way
    '''
    def __init__(self, workingdir : str, readparfileforpaths : bool = True,
                 parfile : str = "parfile.par", verbose = False):
        # super().__init__(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile)

        self.KN = Ejecta(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile,type="kn",verbose=verbose)
        self.GRB = Ejecta(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile,type="grb",verbose=verbose)
        self.PWN = Ejecta(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile,type="pwn",verbose=verbose)
        self.MAG = Magnetar(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile,verbose=verbose)

        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir = workingdir

        if readparfileforpaths:
            self.main_pars, self.main_opts = self.read_main_part_parfile( self.parfile )



    def read_main_part_parfile(self, parfile : str = "parfile.par") -> tuple[dict,dict]:
        main_pars, main_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                            sep1="# -------------------------- main ---------------------------",
                                            sep2="# --------------------------- END ---------------------------")
        return (main_pars,main_opts)

    def reload_parfile(self) -> None:
        self.main_pars, self.main_opts = self.read_main_part_parfile()
        self.KN.reload_parfile(type="kn")
        self.GRB.reload_parfile(type="grb")
        self.MAG.reload_parfile()

    def run(self, path_to_cpp_executable : str, loglevel : str = "info") -> None:
        # this mess is because I did not figure out how $PATH thing works...
        # curdir = os.getcwd()
        # curdir = os.path.dirname(os.path.abspath(__file__))
        # pbadir = curdir.split("PyBlastAfterglowMag")[0]
        # path_to_cpp_executable = pbadir+"PyBlastAfterglowMag"+"/src/pba.out"
        # print(os.getcwd())
        # os.chdir("../../../src/")
        # path_to_executable = "pba.out"
        if not os.path.isfile(path_to_cpp_executable):
            raise IOError("executable is not found: {}".format(path_to_cpp_executable))
        # subprocess.call(path_to_executable, input="")
        # print("{} {} {} {}".format(path_to_cpp_executable, self.workingdir, self.parfile, self.loglevel))
        # subprocess.run(path_to_cpp_executable, input=self.workingdir)
        subprocess.check_call([path_to_cpp_executable, self.workingdir, self.parfile, loglevel])

    def NEW_run(self, P: dict, run: bool = True, process_skymaps: bool = True,
            path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info"):
        """

        :param P: dict with parameters
        :param run:
        :param process_skymaps:
        :param path_to_cpp_executable:
        :param loglevel:
        :return:

        Example P = dict(
            main=dict(n_ism = 1e-2), # main parameters (ISM, grid, etc)
            grb = dict(
                type = "a" # "a" or "pw" for adaptive or piece-wise EATS integrator
                grb=dict(save_dynamics='yes',do_mphys_in_situ='no',do_lc = "no",ebl_tbl_fpath='none',
                    struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
                    skymap_conf=dict(nx=64,ny32,extend_grid=2,fwhm_fac=.5,lat_dist_method="integ",
                        intp_filter=dict(type=None,sigma=2,mode="reflect"), # of "gaussian"
                        hist_filter=dict(type=None,sigma=2,mode="reflect")
                    )
                )
            )
        )

        """
        """
        :
        
                conf = {"nx": 64, "ny": 32, "extend_grid": 2, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                    "intp_filter": {"type": None, "sigma": 2, "mode": 'reflect'},  # "gaussian"
                    "hist_filter": {"type": None, "sigma": 2, "mode": 'reflect'}}
        :param working_dir:
        :param struct:
        :param P:
        :param type:
        :param run:
        :return:
        """
        # clean he temporary direcotry
        if run and os.path.isdir(self.workingdir):
            shutil.rmtree(self.workingdir)
        if not os.path.isdir(self.workingdir):
            os.mkdir(self.workingdir)

        # generate initial data for blast waves
        struct = copy.deepcopy(P["grb"]["struct"])
        del P["grb"]["struct"]
        if (struct["type"] in ["tophat","guassian"]):
            pba_id = JetStruct(
                n_layers_pw=80 if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
                n_layers_a=(1 if struct["struct"] == "tophat" else
                            (20 if not "n_layers_a" in struct.keys() else struct["n_layers_a"])))
        else:
            raise KeyError("Not implemented")

        # save piece-wise EATS ID
        id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=self.workingdir + "id_pw.h5")

        # save adaptive EATS ID
        id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=self.workingdir + "id_a.h5")

        # create new parfile
        P = copy.deepcopy(P)
        type = P["grb"]["type"]
        del P["grb"]["type"]
        P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
        P["grb"]["method_eats"] = "piece-wise" if type == "pw" else "adaptive"
        if (struct["struct"]=="tophat"): P["grb"]["nsublayers"] = 35 # for skymap resolution
        grb_skymap_config = copy.deepcopy(P["grb"]["skymap_conf"])
        del P["grb"]["skymap_conf"]
        create_parfile(working_dir=self.workingdir, P=P)

        # run the code with given parfile
        if not os.path.isfile(path_to_cpp_executable):
            raise IOError("executable is not found: {}".format(path_to_cpp_executable))
        subprocess.check_call([path_to_cpp_executable, self.workingdir, self.parfile, loglevel])

        # process skymap
        if (process_skymaps and self.GRB.opts["do_skymap"] == "yes"):
            prep = ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
            prep.process_singles(infpaths=self.workingdir + "raw_skymap_*.h5",
                                 outfpath=self.GRB.fpath_sky_map,
                                 remove_input=False)


    def clear(self):
        self.KN.clear()
        self.GRB.clear()
        self.MAG.clear()




''' parallel runs TOBE REMOVED '''

class REMOVE_ME:

    def OLD_get_skymap_totfluxes(self, freq, shell=None, time=None):
        self._check_if_loaded_skymap()
        if time is None:
            if (shell is None):
                return np.array(self.skymap_dfile["totalflux at freq={:.4e}".format(freq)])
            else:
                return np.array(self.skymap_dfile["totalflux at freq={:.4e} shell={}".format(freq, shell)])
        else:
            if (shell is None):
                arr = np.array(self.skymap_dfile["totalflux at freq={:.4e}".format(freq)])
                if (len(arr)==1):
                    return arr[0]
                if (time < self.get_skymap_times().min()):
                    raise ValueError(f"time {time} < get_skymap_times().min()={self.get_skymap_times().min()}")
                if (time > self.get_skymap_times().max()):
                    raise ValueError(f"time {time} > get_skymap_times().max()={self.get_skymap_times().max()}")
                # self.get_skymap_totfluxes(freq=freq, shell=None, time=None)
                val = arr[find_nearest_index(self.get_skymap_times(), time)]
                return val
            else:
                raise KeyError("Not finished...")


    def OLD_get_skymap_cm(self, all_xrs, all_yrs, all_zz):
        dfile = self.get_skymap_obj()
        xc_m, yc_m = compute_position_of_the_flux_centroid(all_xrs, all_yrs, all_zz, float(dfile.attrs["d_l"]))
        return (xc_m, yc_m)



    def OLD_get_combained_skymaps_adjusted_to_other(self, time, freq, other_pb_instance, nx=100, ny=50):

        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, retrun_edges=True)
        all_x_m1, all_y_m1, all_fluxes_m1 \
            = other_pb_instance.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        _, _, int_zz_m1 = other_pb_instance.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                                           hist_or_int="hist", shells=True, nx=nx,
                                                           ny=ny, edges_x=edges_x, edges_y=edges_y)

        return (int_x, int_y, int_zz, int_zz_m1)

    def OLD_get_combined_spectral_map(self, time, freq, nx=100, ny=50, extend=2):
        freqs = self.get_skymap_freqs()
        assert len(freqs) > 2
        idx = find_nearest_index(freqs, freq)
        assert idx != len(freqs)-1
        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, extend=extend, retrun_edges=True)

        all_x_m1, all_y_m1, all_fluxes_m1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        all_x_p1, all_y_p1, all_fluxes_p1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        # ffreq_i = np.full(int_zz.shape, freqs[idx])
        ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        ffreq = np.full(int_zz.shape, freqs[idx])
        # num = np.log10(int_zz) - np.log10(int_zz_m1)
        num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        # denum = np.log10(ffreq) - np.log10(ffreq_im1)
        denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        # int_zz = num / denum
        int_zz = 1.0 * num / denum / 2.0

        return (int_x, int_y, int_zz)


    def OLD_TOREMOVE_get_skymap_old(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False, renormalize=True, normtype="pw"):

        # nx = 200
        # ny = 100
        # min_mu = 1e-4 # limit for mu that if too small blow up the intensity for some reason...
        # min_mu_frac = 0.85
        # # limit_data = True
        times = self.get_skymap_times()
        freqs = self.get_skymap_freqs()
        dfile = self.get_skymap_obj()
        nshells = int(dfile.attrs["nshells"])
        d_l = float(self.get_skymap_obj().attrs["d_l"])
        if ((not time is None) and (not time in times)):
            raise ValueError(
                "time={} day is not in the list for skypams={} days".format(time / cgs.day, times / cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
        if ((not ishell is None) and (ishell > nshells - 1)):
            raise ValueError("shell={} is beying skymap nshells={}".format(ishell, nshells))
        if ((not time is None) and (not freq is None)):
            # print(dfile.keys())
            ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
            if (not (ishell is None)):
                # r_i = np.array(ddfile["r"][ishell])
                mu_i = np.array(ddfile["mu"][ishell])
                xrs_i = np.array(ddfile["xrs"][ishell])# * cgs.rad2mas / d_l  # m -> mas
                yrs_i = np.array(ddfile["yrs"][ishell])# * cgs.rad2mas / d_l  # m -> mas
                int_i = np.array(ddfile["intensity"][ishell])# * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                # gam_i = np.array(ddfile["gamma"][ishell])
                # B_i = np.array(ddfile["B"][ishell])
                # tb_i = np.array(ddfile["tburst"][ishell])
                # theta_i = np.array(ddfile["theta"][ishell])
                # phi_i = np.array(ddfile["phi"][ishell])

                # mu_i = mu_i[int_i > 0]
                # xrs_i = xrs_i[int_i > 0] # m -> mas
                # yrs_i = yrs_i[int_i > 0]  # m -> mas
                # int_i = int_i[int_i > 0]  # -> mJy / mas^2

                ncells = int(len(xrs_i)/2)


                if (np.sum(int_i)==0.):
                    raise ValueError(" x_arr or y_arr arrays in the image all FULL 0. Cannot re-interpolate!")
                if (xrs_i.min()==0 and xrs_i.max()==0):
                    raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")
                if (yrs_i.min()==0 and yrs_i.max()==0):
                    raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")

                if remove_mu:
                    if self.verb: print("Removing 'mu' from ejecta skymap")
                    int_i *= np.abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!

                if renormalize:
                    if self.verb: print("Renormalizing ejecta skymap (shell my shell separately)")
                    fnus = self.get_skymap_totfluxes(freq=freq, shell=ishell)
                    fnu = fnus[find_nearest_index(self.get_skymap_times(), time)]
                    all_fluxes_arr = np.array(int_i)
                    delta_x = np.array(xrs_i).max() - np.array(xrs_i).min()
                    delta_y = np.array(yrs_i).max() - np.array(yrs_i).min()
                    dfnu = fnu / (delta_x * delta_y)
                    if verbose:
                        print("Ejecta shell {}".format(ishell))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(xrs_i).min(), np.array(xrs_i).max()))
                        print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(yrs_i).min(), np.array(yrs_i).max()))
                        print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                    int_i = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                    fnus_tot = int_i + fnus

                return (xrs_i, yrs_i, int_i)

                # plt.figure(figsize=(4.6, 3.2))
                # plt.semilogy(mu_i, int_i / int_i.max(), '.')
                # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
                # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                # plt.xlabel(
                #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                # plt.ylabel(r"$I/I_{\rm max}$")
                # plt.tight_layout()
                # plt.show()

                # nx = np.complex(0, nx)
                # ny = np.complex(0, ny)
                # grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx, yrs_i.min():yrs_i.max():ny]
                # i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

                # if int_or_hist == "int":
                #     nx = np.complex(0, nx)
                #     ny = np.complex(0, ny)
                #     grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #                      yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                # else:
                #     nx = 100
                #     ny = 100
                #     nx = np.complex(0, nx + 1)
                #     ny = np.complex(0, ny + 1)
                #     edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                #     edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #     # print(i_zz.shape)
            else:
                i_shells = []
                # total_fluxes = []
                for ish in range(nshells):
                    xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                    # mu_i = np.array(ddfile["mu"][ish])

                    sumi = np.sum(int_i)


                    if (sumi == 0):
                        continue

                    nnonzero = np.count_nonzero(sumi)

                    if (xrs_i.min()==0 and xrs_i.max()==0 and nnonzero != 1):
                        raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")
                    if (yrs_i.min()==0 and yrs_i.max()==0 and nnonzero != 1):
                        raise ValueError(" x_arr arrays in the image all FULL 0. Cannot re-interpolate!")


                    i_shells.append(ish)
                    # if remove_mu:
                    #     print("Removing 'mu' from ejecta skymap")
                    #     int_i *= np.abs( mu_i )
                all_xrs, all_yrs, all_zz = [], [], []
                for ii, ish in enumerate(i_shells):
                    mu_i = np.array(ddfile["mu"][ish])
                    xrs_i = np.array(ddfile["xrs"][ish]) #* cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) #* cgs.rad2mas / d_l  # m -> mas
                    # rs_i = np.array(ddfile["r"][ish])
                    # cthetas_i = np.array(ddfile["ctheta"][ish])
                    # cphis_i = np.array(ddfile["cphi"][ish])
                    int_i = np.array(ddfile["intensity"][ish]) #* ( d_l ** 2 / cgs.rad2mas ** 2 )  # * dfile.attrs["d_L"] ** 2

                    # mu_i = mu_i[int_i > 0]
                    # xrs_i = xrs_i[int_i > 0] # m -> mas
                    # yrs_i = yrs_i[int_i > 0]  # m -> mas
                    # int_i = int_i[int_i > 0]  # -> mJy / mas^2

                    if remove_mu:
                        # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                        # int_i[idx1] =
                        # print("len(gam[idx1])={}".format(gam_i[idx1]))
                        # idx = int_i > int_i.max() * min_mu_frac
                        # idx = idx1
                        # if verbose:
                        #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                        #         print('\t', mu_i[idx1])
                        #         print('\t', mu_i[idx])
                        #         # exit(1)
                        # if verbose:
                        #     if len(idx) < 2:
                        #         print("No excess = {}".format(ii))
                        #     print("len(gam[idx])={}".format(gam_i[idx]))
                        #     print("Gamma={}".format(gam_i[idx]))
                        #     print("beta={}".format(get_beta(gam_i[idx])))
                        #     print("mu_i={}".format(mu_i[idx]))
                        #     print("r_i={}".format(r_i[idx]))
                        #     print("B={}".format(B_i[idx]))
                        #     print("tb_i={}".format(tb_i[idx]))
                        # print(np.abs(mu_i).min())
                        # print(theta_i[np.abs(mu_i) < 1e-5])
                        # print(phi_i[np.abs(mu_i) < 1e-5])
                        # print(np.abs(mu_i).min())
                        # plt.figure(figsize=(4.6, 3.2))
                        # plt.semilogy(mu_i, int_i / int_i.max(), '.')
                        # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
                        # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                        # plt.xlabel(
                        #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                        # plt.ylabel(r"$I/I_{\rm max}$")
                        # plt.tight_layout()
                        # plt.show()

                        # int_i[mu_i < 1e-3] = 0.
                        int_i *= abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                        # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                        # import matplotlib.pyplot as plt
                        # plt.figure(figsize=(4.6,3.2))
                        # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                        # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                        # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                        # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                        # plt.ylabel(r"$I/I_{\rm max}$")
                        # plt.tight_layout()
                        # plt.show()
                        # int_i[idx] = 0.
                    all_xrs.append(xrs_i)
                    all_yrs.append(yrs_i)
                    all_zz.append(int_i)

                    import matplotlib.pyplot as plt
                    from matplotlib.colors import Normalize, LogNorm
                    from matplotlib import cm, rc, rcParams
                    import matplotlib.colors as colors
                    from mpl_toolkits.axes_grid1 import make_axes_locatable
                    from mpl_toolkits.axes_grid1 import ImageGrid
                    from mpl_toolkits.axisartist.grid_finder import MaxNLocator
                    from matplotlib.colors import BoundaryNorm
                    from matplotlib.ticker import MaxNLocator
                    # plt.close()
                    # layers = int(self.get_dyn_obj().attrs["nlayers"])-1
                    # plt.semilogx(self.get_dyn_arr("tt", ilayer=0) / cgs.day, self.get_dyn_arr("ctheta", ilayer=0),
                    #              label="nl={}".format(0))
                    # plt.semilogx(self.get_dyn_arr("tt", ilayer=layers - 1) / cgs.day,
                    #              self.get_dyn_arr("ctheta", ilayer=layers - 1), label="nl={}".format(layers - 1))
                    # plt.semilogx(self.get_dyn_arr("tt", ilayer=layers - int(layers / 2)) / cgs.day,
                    #              self.get_dyn_arr("ctheta", ilayer=layers - int(layers / 2)),
                    #              label="nl={}".format(layers - int(layers / 2)))
                    # plt.xlabel("Time [days]")
                    # plt.grid()
                    # plt.legend()
                    # plt.show()

                    # plt.close()
                    # fig = plt.figure()
                    # cmap = cm.get_cmap('inferno')
                    # my_norm = LogNorm(int_i.max() * 1e-1, int_i.max())
                    # ax = fig.add_subplot(projection='3d')
                    # # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
                    # # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(rs_i).flatten(), c=cmap(my_norm(int_i.flatten())))
                    # # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
                    # # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cthetas_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
                    # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),np.log10(int_i.flatten()),  c=cmap(my_norm(int_i.flatten())))
                    # # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cphis_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
                    # ax.set_xlabel('X Label')
                    # ax.set_ylabel('Y Label')
                    # ax.set_zlabel('I Label')
                    # ax.set_zlim(np.log10(int_i.flatten()).max()*1.-2,np.log10(int_i.flatten()).max())
                    # plt.show()

                # Problem: changing nlayers changes the image/image, Fnu per pixel; Solution:
                if (renormalize and normtype=="pw"):
                    print("Renormalizing ejecta skymap (shell by shell separately)")
                    fnus_tot = np.zeros_like(self.get_skymap_times())
                    for i, i_ish in enumerate(i_shells):
                        fnus = self.get_skymap_totfluxes(freq=freq, shell=i_ish)
                        fnu = fnus[find_nearest_index(self.get_skymap_times(), time)]
                        all_fluxes_arr = np.array(all_zz[i])
                        delta_x = np.array(all_xrs[i]).max() - np.array(all_xrs[i]).min()
                        delta_y = np.array(all_yrs[i]).max() - np.array(all_yrs[i]).min()
                        dfnu = fnu / (delta_x * delta_y)
                        if (~np.isfinite(dfnu)):
                            for i, val in enumerate(np.array(all_yrs[i])):
                                if (~np.isfinite(val)):
                                    print(f"i={i} val={val}")
                            raise ValueError(f"dfnu={dfnu} ishell={i_ish}")
                        if verbose:
                            print("SHELL {}".format(i_ish))
                            print("\tfnu = {:.2e} ".format(fnu))
                            print("\tdfnu = {:.2e} ".format(dfnu))
                            print("\tall_fluxes_arr.max() = {:.2e} ".format(all_fluxes_arr.max()))
                            print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(all_xrs[i]).min(), np.array(all_xrs[i]).max()))
                            print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(all_yrs[i]).min(), np.array(all_yrs[i]).max()))
                            print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                            print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        all_zz[i] = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                        fnus_tot = fnus_tot + fnus
                elif (renormalize and normtype=="a"):
                    print("Renormalizing ejecta skymap (shell by shell separately)")
                    fnus_tot = np.zeros_like(self.get_skymap_times())
                    for i, i_ish in enumerate(i_shells):
                        fnus = self.get_skymap_totfluxes(freq=freq, shell=i_ish)
                        fnu = fnus[find_nearest_index(self.get_skymap_times(), time)]
                        all_fluxes_arr_old = np.array(all_zz[i])
                        delta_x = np.array(all_xrs[i]).max() - np.array(all_xrs[i]).min()
                        delta_y = np.array(all_yrs[i]).max() - np.array(all_yrs[i]).min()
                        # grad_x = np.gradient(all_xrs[i])
                        # diff_x = np.append(np.diff(all_xrs[i]),0)
                        # diff_y = np.append(np.diff(all_yrs[i]),0)
                        all_fluxes_arr = np.array(all_zz[i])# * np.abs(diff_x) * np.abs(diff_y)
                        dfnu = fnu / (delta_x * delta_y)
                        if verbose:
                            print("SHELL {}".format(i_ish))
                            print("\tfnu = {:.2e} ".format(fnu))
                            print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(all_xrs).min(), np.array(all_xrs).max()))
                            print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(all_yrs).min(), np.array(all_yrs).max()))
                            print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                            print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        all_zz[i] = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                        fnus_tot = fnus_tot + fnus
                elif (renormalize):
                    raise KeyError(f"norm type {normtype} is not recognized")
                return (all_xrs, all_yrs, all_zz)

                # assess what is the combined grid extend (for all images)
                # xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
                # ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
                # i_min, i_max = [], []
                # i_shells = []
                # for ish in range(nshells):
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
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
                # xmin = np.array(xmin)
                # ymin_neg = np.array(ymin_neg)
                # ymin_pos = np.array(ymin_pos)
                # ymax = np.array(ymax)
                # ymin = np.array(ymin)
                # i_min = np.array(i_min)
                # i_max = np.array(i_max)
                #
                # if verbose:
                #     for i in range(len(i_shells)):
                #         print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i], xmax[i]))
                #         print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],  ymax[i]))
                #     print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
                #     print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
                #     print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
                #
                # edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
                # edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
                #
                #
                # x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
                #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
                #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                #
                # # edges_x = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # # edges_y = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                # # plt.plot(edges_x, edges_y, marker='.', ls='none')
                # # plt.show()
                #
                # if verbose:
                #     print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                #                                                                     x_grid[x_grid > 0].min(),
                #                                                                     x_grid.max()))
                #     print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                #                                                                     y_grid[y_grid > 0].min(),
                #                                                                     y_grid.max()))
                # xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
                # if int_or_hist == "int":
                #     zz = np.zeros_like((xx_grid))
                # else:
                #     zz = np.zeros((len(edges_x)-1,len(edges_y)-1))
                # # interpolate onto the grid that covers all images
                # all_xrs, all_yrs, all_zz = [], [], []
                # for ii, ish in enumerate(i_shells):  # range(len(i_shells))
                #     if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
                #     r_i = np.array(ddfile["r"][ish])
                #     mu_i = np.array(ddfile["mu"][ish])
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2) # * dfile.attrs["d_L"] ** 2
                #     gam_i = np.array(ddfile["gamma"][ish])
                #     B_i = np.array(ddfile["B"][ish])
                #     tb_i = np.array(ddfile["tburst"][ish])
                #     theta_i = np.array(ddfile["theta"][ish])
                #     theta_ij = np.array(ddfile["theta_j"][ish])
                #     theta0_i = np.array(ddfile["theta0"][ish])
                #     phi_i = np.array(ddfile["phi"][ish])
                #
                #
                #     # plt.plot(mu_i, theta0_i, marker='.', ls='none')
                #     # plt.show()
                #     # int_i[mu_i < 5e-2] = 0.
                #     # int_i *= abs(mu_i)
                #
                #     if remove_mu:
                #         # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                #         # int_i[idx1] =
                #         # print("len(gam[idx1])={}".format(gam_i[idx1]))
                #         # idx = int_i > int_i.max() * min_mu_frac
                #         # idx = idx1
                #         # if verbose:
                #         #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                #         #         print('\t', mu_i[idx1])
                #         #         print('\t', mu_i[idx])
                #         #         # exit(1)
                #         # if verbose:
                #         #     if len(idx) < 2:
                #         #         print("No excess = {}".format(ii))
                #         #     print("len(gam[idx])={}".format(gam_i[idx]))
                #         #     print("Gamma={}".format(gam_i[idx]))
                #         #     print("beta={}".format(get_beta(gam_i[idx])))
                #         #     print("mu_i={}".format(mu_i[idx]))
                #         #     print("r_i={}".format(r_i[idx]))
                #         #     print("B={}".format(B_i[idx]))
                #         #     print("tb_i={}".format(tb_i[idx]))
                #         # print(np.abs(mu_i).min())
                #         # print(theta_i[np.abs(mu_i) < 1e-5])
                #         # print(phi_i[np.abs(mu_i) < 1e-5])
                #         # print(np.abs(mu_i).min())
                #
                #         # int_i[mu_i < 1e-3] = 0.
                #         int_i *= abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                #         # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                #         # plt.figure(figsize=(4.6,3.2))
                #         # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                #         # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                #         # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                #         # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                #         # plt.ylabel(r"$I/I_{\rm max}$")
                #         # plt.tight_layout()
                #         # plt.show()
                #         # int_i[idx] = 0.
                #     # int_i /= abs(mu_i)[abs(mu_i)>0.01]
                #     # else:
                #     #     int_i2 = int_i
                #
                #     if int_or_hist == "int":
                #         # nx = np.complex(0, nx)
                #         # ny = np.complex(0, ny)
                #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
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
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #
                #     # i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #     zz += i_zz
                #     all_xrs.append(xrs_i)
                #     all_yrs.append(yrs_i)
                #     all_zz.append(int_i)
                # if verbose:
                #     print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
                #           .format(zz.min(), zz.max(), np.sum(zz), float(np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
                #
                # all_xrs = np.concatenate(all_xrs)
                # all_yrs = np.concatenate(all_yrs)
                # all_zz = np.concatenate(all_zz)
                # if int_or_hist == "int":
                #
                #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
                # else:
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, zz, all_xrs, all_yrs, all_zz)
        else:
            raise NotImplementedError("Not implemented")

    def OLD_TOREMOVE_get_skymap(self, time, freq, verbose=False, remove_zeros=True, return_sph_coords=False):
        times = self.get_skymap_times()
        freqs = self.get_skymap_freqs()
        dfile = self.get_skymap_obj()
        nshells = int(dfile.attrs["nshells"])
        # d_l = float(self.get_skymap_obj().attrs["d_l"])
        if ((not time is None) and (not time in times)):
            raise ValueError(
                "time={} day is not in the list for skypams={} days".format(time / cgs.day, times / cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))

        ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
        i_shells = []

        # loop over shells and get shells with sum(int)>0
        for ish in range(nshells):
            sumi = np.sum(np.array(ddfile["intensity"][ish]))
            if (sumi == 0):
                if verbose:
                    print(f"Skipping empty sum(intensity)=0 shell ish={ish}")
                continue
            i_shells.append(ish)

        ncells = []
        all_xrs, all_yrs, all_zz = [], [], []
        all_theta, all_phi, all_r = [], [], []
        # loop over non-empy shells and collect data
        for ii, ish in enumerate(i_shells):
            mu_i = np.array(ddfile["mu"][ish])
            # cartesian coordiantes on the projected plane
            xrs_i = np.array(ddfile["xrs"][ish])# * cgs.rad2mas / d_l  # m -> mas
            yrs_i = np.array(ddfile["yrs"][ish])# * cgs.rad2mas / d_l  # m -> mas
            int_i = np.array(ddfile["intensity"][ish])# * ( d_l ** 2 / cgs.rad2mas ** 2 )  # * dfile.attrs["d_L"] ** 2
            # spherical coordinates
            ctheta_i = np.array(ddfile["ctheta"][ish])
            cphi_i = np.array(ddfile["cphi"][ish])
            r_i = np.array(ddfile["r"][ish])

            all_xrs.append(xrs_i)
            all_yrs.append(yrs_i)
            all_zz.append(int_i)

            if (return_sph_coords):
                all_theta.append(ctheta_i)
                all_phi.append(cphi_i)
                all_r.append(r_i)

            ncells.append( int(len(xrs_i) / 2) )
            if (len(xrs_i) % 2 > 0):
                raise ValueError(f"len(xrs) is expected to be even (2*ncells). Got={len(xrs_i)}")

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

            if (len(all_xrs_pjcj[-1]) == 0):
                raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
                                 f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_pj  = all_theta[ii][:ncells[ii]]
                _cphi_pj    = all_phi[ii][:ncells[ii]]
                _r_pj       = all_r[ii][:ncells[ii]]

                all_ctheta_pjcj.append(_ctheta_pj[ _zz_pj > 0 ] if (remove_zeros) else _ctheta_pj )
                all_cphi_pjcj.append(  _cphi_pj[   _zz_pj > 0 ] if (remove_zeros) else _cphi_pj )
                all_r_pjcj.append(     _r_pj[      _zz_pj > 0 ] if (remove_zeros) else _r_pj )

        print(f"Principle only maxs = {maxs}")

        # process counter jet
        for ii, ish in enumerate(i_shells):
            _xrs_cj = all_xrs[ii][ncells[ii]:]
            _yrs_cj = all_yrs[ii][ncells[ii]:]
            _zz_cj  = all_zz[ii][ncells[ii]:]
            maxs.append(np.max(_zz_cj))

            all_xrs_pjcj.append(_xrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _xrs_cj )
            all_yrs_pjcj.append(_yrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _yrs_cj )
            all_zz_pjcj.append(  _zz_cj[ _zz_cj > 0 ] if (remove_zeros) else _zz_cj )

            if (len(all_xrs_pjcj[-1]) == 0):
                raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
                                 f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_cj  = all_theta[ii][ncells[ii]:]
                _cphi_cj    = all_phi[ii][ncells[ii]:]
                _r_cj       = all_r[ii][ncells[ii]:]

                all_ctheta_pjcj.append(_ctheta_cj[ _zz_cj > 0 ] if (remove_zeros) else _ctheta_cj )
                all_cphi_pjcj.append(  _cphi_cj[   _zz_cj > 0 ] if (remove_zeros) else _cphi_cj )
                all_r_pjcj.append(     _r_cj[      _zz_cj > 0 ] if (remove_zeros) else _r_cj )

        print(f"Principle & counter only maxs = {maxs}")

        if (return_sph_coords):
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj, all_ctheta_pjcj, all_cphi_pjcj, all_r_pjcj)
        else:
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj)


    def OLD_get_skymap_OLD(self, time : float, freq : float, verbose=False, remove_zeros=True, return_sph_coords=False):
        times = self.get_skymap_times()
        freqs = self.get_skymap_freqs()
        dfile = self.get_skymap_obj()
        nshells = int(dfile.attrs["nshells"])
        if (time is None):
            raise NotImplementedError("time=None or freq=None is not longer supported")

        # ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
        i_shells = []

        # loop over shells and get shells with sum(int)>0
        for ish in range(nshells):
            ddfile = dfile["shell={} time={:.4e} freq={:.4e}".format(ish, time, freq)]
            sumi = np.sum(np.array(ddfile["intensity"]))
            if (sumi == 0):
                if verbose:
                    print(f"Skipping empty sum(intensity)=0 shell ish={ish}")
                continue
            i_shells.append(ish)

        ncells = []
        all_xrs, all_yrs, all_zz = [], [], []
        all_theta, all_phi, all_r = [], [], []
        # loop over non-empy shells and collect data
        for ii, ish in enumerate(i_shells):
            ddfile = dfile["shell={} time={:.4e} freq={:.4e}".format(ish, time, freq)]
            mu_i = np.array(ddfile["mu"])
            # cartesian coordiantes on the projected plane
            xrs_i = np.array(ddfile["xrs"])# * cgs.rad2mas / d_l  # m -> mas
            yrs_i = np.array(ddfile["yrs"])# * cgs.rad2mas / d_l  # m -> mas
            int_i = np.array(ddfile["intensity"])# * ( d_l ** 2 / cgs.rad2mas ** 2 )  # * dfile.attrs["d_L"] ** 2
            # spherical coordinates
            ctheta_i = np.array(ddfile["ctheta"])
            cphi_i = np.array(ddfile["cphi"])
            r_i = np.array(ddfile["r"])

            all_xrs.append(xrs_i)
            all_yrs.append(yrs_i)
            all_zz.append(int_i)

            if (return_sph_coords):
                all_theta.append(ctheta_i)
                all_phi.append(cphi_i)
                all_r.append(r_i)

            ncells.append( int(len(xrs_i) / 2) )
            if (len(xrs_i) % 2 > 0):
                raise ValueError(f"len(xrs) is expected to be even (2*ncells). Got={len(xrs_i)}")

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

            if (len(all_xrs_pjcj[-1]) == 0):
                raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
                                 f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_pj  = all_theta[ii][:ncells[ii]]
                _cphi_pj    = all_phi[ii][:ncells[ii]]
                _r_pj       = all_r[ii][:ncells[ii]]

                all_ctheta_pjcj.append(_ctheta_pj[ _zz_pj > 0 ] if (remove_zeros) else _ctheta_pj )
                all_cphi_pjcj.append(  _cphi_pj[   _zz_pj > 0 ] if (remove_zeros) else _cphi_pj )
                all_r_pjcj.append(     _r_pj[      _zz_pj > 0 ] if (remove_zeros) else _r_pj )

        print(f"Principle only maxs = {maxs}")

        # process counter jet
        for ii, ish in enumerate(i_shells):
            _xrs_cj = all_xrs[ii][ncells[ii]:]
            _yrs_cj = all_yrs[ii][ncells[ii]:]
            _zz_cj  = all_zz[ii][ncells[ii]:]
            maxs.append(np.max(_zz_cj))

            all_xrs_pjcj.append(_xrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _xrs_cj )
            all_yrs_pjcj.append(_yrs_cj[ _zz_cj > 0 ] if (remove_zeros) else _yrs_cj )
            all_zz_pjcj.append(  _zz_cj[ _zz_cj > 0 ] if (remove_zeros) else _zz_cj )

            if (len(all_xrs_pjcj[-1]) == 0):
                raise ValueError(f"Empty shell {ish} ncells[ii]={ncells[ii]} "
                                 f"len(all_xrs[ii][:ncells[ii]]={all_xrs[ii][:ncells[ii]]}); after 'remove_zeros' {len(all_zz_pjcj[-1])}")

            if (return_sph_coords):
                _ctheta_cj  = all_theta[ii][ncells[ii]:]
                _cphi_cj    = all_phi[ii][ncells[ii]:]
                _r_cj       = all_r[ii][ncells[ii]:]

                all_ctheta_pjcj.append(_ctheta_cj[ _zz_cj > 0 ] if (remove_zeros) else _ctheta_cj )
                all_cphi_pjcj.append(  _cphi_cj[   _zz_cj > 0 ] if (remove_zeros) else _cphi_cj )
                all_r_pjcj.append(     _r_cj[      _zz_cj > 0 ] if (remove_zeros) else _r_cj )

        print(f"Principle & counter only maxs = {maxs}")

        if (return_sph_coords):
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj, all_ctheta_pjcj, all_cphi_pjcj, all_r_pjcj)
        else:
            return (all_xrs_pjcj, all_yrs_pjcj, all_zz_pjcj)

class PBA_BASE_OLD:
    '''
        Pars initial data and parameters into C++ code interface
        Load result files
    '''
    def __init__(self, workingdir, readparfileforpaths=True,parfile="pafile.par"):

        if not os.path.isdir(workingdir):
            raise IOError("Working directory not found {}".format(workingdir))
        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir_mag = workingdir
        self.res_dir_kn = workingdir
        self.res_dir_grb = workingdir
        self.grb_prefix = "grb_"
        self.kn_prefix = "kn_"
        # self.loglevel = loglevel
        # ------------------ MAIN
        if readparfileforpaths:
            self.main_pars, self.main_opts = self.read_main_part_parfile( self.parfile)
        # ------------------ MAGNETAR
        self.fpath_mag = None
        if readparfileforpaths:
            self.mag_pars, self.mag_opts = self.read_magnetar_part_parfile( self.parfile)
        else:
            self.fpath_mag = self.res_dir_mag + "magnetar.h5"
        # ------------------ GRB
        self.fpath_grb_dyn = None
        self.fpath_grb_spec = None
        self.fpath_grb_light_curve = None
        self.fpath_grb_sky_map = None
        if readparfileforpaths:
            self.grb_pars, self.grb_opts = self.read_grb_part_parfile( self.parfile)
        else:
            self.fpath_grb_dyn = self.res_dir_kn + self.grb_prefix + "dynamics_layers.h5"
            self.fpath_grb_spec = self.res_dir_kn + self.grb_prefix + "spectra.h5"
            self.fpath_grb_light_curve = self.res_dir_kn + self.grb_prefix + "lightcurves_layers.h5"
            self.fpath_grb_sky_map = self.res_dir_kn + self.grb_prefix + "skymap.h5"
        # ------------------ EJECTA
        self.fpath_kn_dyn = None
        self.fpath_kn_spec = None
        self.fpath_kn_light_curve = None
        self.fpath_kn_sky_map = None
        if readparfileforpaths:
            self.kn_pars, self.kn_opts = self.read_kn_part_parfile( self.parfile)
        else:
            self.fpath_kn_dyn = self.res_dir_kn + self.kn_prefix + "dynamics_layers.h5"
            self.fpath_kn_spec = self.res_dir_kn + self.kn_prefix + "spectra_layers.h5"
            self.fpath_kn_light_curve = self.res_dir_kn + self.kn_prefix + "lightcurves_layers.h5"
            self.fpath_kn_sky_map = self.res_dir_kn + self.kn_prefix + "skymap.h5"
        # ------------------- MAGNETAR
        # TODO
        # -----------------
        self.mag_dfile = None
        self.grb_dyn_dfile = None
        self.grb_spec_dfile = None
        self.grb_lc_dfile = None
        self.grb_skymap_dfile = None
        self.kn_dyn_dfile = None
        self.kn_spec_dfile = None
        self.kn_lc_dfile = None
        self.kn_skymap_dfile = None
        # ----------------

    def clear(self):
        # self.overwrite = True
        if (not self.mag_dfile is None):
            self.mag_dfile.close()
            self.mag_dfile = None
        if (not self.grb_dyn_dfile is None):
            self.grb_dyn_dfile.close()
            self.grb_dyn_dfile = None
        if (not self.grb_spec_dfile is None):
            self.grb_spec_dfile.close()
            self.grb_spec_dfile = None
        if (not self.grb_lc_dfile is None):
            self.grb_lc_dfile.close()
            self.grb_lc_dfile = None
        if (not self.grb_skymap_dfile is None):
            self.grb_skymap_dfile.close()
            self.grb_skymap_dfile = None
        if (not self.kn_dyn_dfile is None):
            self.kn_dyn_dfile.close()
            self.kn_dyn_dfile = None
        if (not self.kn_spec_dfile is None):
            self.kn_spec_dfile.close()
            self.kn_spec_dfile = None
        if (not self.kn_lc_dfile is None):
            self.kn_lc_dfile.close()
            self.kn_lc_dfile = None
        if (not self.kn_skymap_dfile is None):
            self.kn_skymap_dfile.close()
            self.kn_skymap_dfile = None
    def read_main_part_parfile(self, parfile="parfile.par"):
        main_pars, main_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                            sep1="# -------------------------- main ---------------------------",
                                            sep2="# --------------------------- END ---------------------------")
        return (main_pars,main_opts)
    def read_magnetar_part_parfile(self, parfile="parfile.par"):
        mag_pars, mag_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                          sep1="# ------------------------ Magnetar -------------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_mag" in mag_opts.keys(): self.fpath_mag = self.res_dir_mag + mag_opts["fname_mag"]
        return (mag_pars, mag_opts)
    def read_grb_part_parfile(self,parfile="parfile.par"):
        grb_pars, grb_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                          sep1="# ---------------------- GRB afterglow ----------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_dyn" in grb_opts.keys(): self.fpath_grb_dyn = self.res_dir_grb + grb_opts["fname_dyn"]
        if "fname_spec" in grb_opts.keys(): self.fpath_grb_spec = self.res_dir_grb + grb_opts["fname_spec"]
        if "fname_light_curve" in grb_opts.keys(): self.fpath_grb_light_curve = self.res_dir_grb + grb_opts["fname_light_curve"]
        if "fname_sky_map" in grb_opts.keys(): self.fpath_grb_sky_map = self.res_dir_grb + grb_opts["fname_sky_map"]
        return (grb_pars,grb_opts)
    def read_kn_part_parfile(self,parfile="parfile.par"):
        kn_pars, kn_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                        sep1="# ----------------------- kN afterglow ----------------------",
                                        sep2="# --------------------------- END ---------------------------")
        if "fname_dyn" in kn_opts.keys(): self.fpath_kn_dyn = self.res_dir_kn + kn_opts["fname_dyn"]
        if "fname_spec" in kn_opts.keys(): self.fpath_kn_spec = self.res_dir_kn + kn_opts["fname_spec"]
        if "fname_light_curve" in kn_opts.keys(): self.fpath_kn_light_curve = self.res_dir_kn + kn_opts["fname_light_curve"]
        if "fname_sky_map" in kn_opts.keys(): self.fpath_kn_sky_map = self.res_dir_kn + kn_opts["fname_sky_map"]
        return (kn_pars,kn_opts)
    # def modify_parfile(self, base_parfile,
    #                    new_main_pars:dict, new_main_opts:dict,
    #                    new_grb_pars:dict, new_grb_opts:dict,
    #                    new_kn_pars:dict, new_kn_opts:dict):
    #     if(len(new_main_pars.keys())+len(new_main_opts.keys())>0):
    #         modify_parfile(newpars=new_main_pars,newopts=new_main_opts,workingdir=self.workingdir,
    #                        comment="#",fname=de,
    #                        newfname="parfile2.par",
    #                        sep1="# -------------------------- main ---------------------------",
    #                        sep2="# --------------------------- END ---------------------------")
    #         copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #         copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #         os.remove(self.workingdir+"parfile2.par")
    #
    #
    #
    #
    #     if (reload_parfile): self.main_pars, self.main_opts = self.read_main_part_parfile()


    # def modify_main_part_parfile(self, newpars : dict, newopts : dict, parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     # if save_old:
    #     #     copyfile(self.workingdir+parfile,self.workingdir+"old_{}".format(parfile))
    #     if (keep_old and parfile == newparfile):
    #         raise NameError("Cannot keep old parfile if the new name is the same as old")
    #     copyfile(self.workingdir+parfile,self.workingdir+"tmp_{}".format(parfile))
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",
    #                    fname="tmp_{}".format(parfile),
    #                    newfname="tmp_mod_{}".format(newparfile),
    #                    sep1="# -------------------------- main ---------------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     if not keep_old:
    #         os.remove(self.workingdir+parfile)
    #     copyfile(self.workingdir+"tmp_mod_{}".format(newparfile),self.workingdir+newparfile)
    #     os.remove(self.workingdir+"tmp_{}".format(parfile))
    #     os.remove(self.workingdir+"tmp_mod_{}".format(newparfile))
    #     # if not save_old:
    #     #     os.remove(self.workingdir+"parfile2.par")
    #     # copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     # copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     # os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.main_pars, self.main_opts = self.read_main_part_parfile()
    # def modify_grb_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     if (keep_old and parfile == newparfile):
    #         raise NameError("Cannot keep old parfile if the new name is the same as old")
    #     copyfile(self.workingdir+parfile,self.workingdir+"tmp_{}".format(parfile))
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",
    #                    fname="tmp_{}".format(parfile),
    #                    newfname="tmp_mod_{}".format(newparfile),
    #                    sep1="# ---------------------- GRB afterglow ----------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     if not keep_old:
    #         os.remove(self.workingdir+parfile)
    #     copyfile(self.workingdir+"tmp_mod_{}".format(newparfile),self.workingdir+newparfile)
    #     os.remove(self.workingdir+"tmp_{}".format(parfile))
    #     os.remove(self.workingdir+"tmp_mod_{}".format(newparfile))
    #     # copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     # copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     # os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
    # def modify_mag_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",fname=parfile,
    #                    newfname="parfile2.par",
    #                    sep1="# ------------------------ Magnetar -------------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
    # def modify_kn_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",fname=parfile,
    #                    newfname="parfile2.par",
    #                    sep1="# ----------------------- kN afterglow ----------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.kn_pars, self.kn_opts = self.read_kn_part_parfile()
    def reload_parfile(self):
        self.main_pars, self.main_opts = self.read_main_part_parfile()
        self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
        self.kn_pars, self.kn_opts = self.read_kn_part_parfile()
    def run(self, loglevel="info"):
        # this mess is because I did not figure out how $PATH thing works...
        curdir = os.getcwd()
        pbadir = curdir.split("PyBlastAfterglowMag")[0]
        path_to_cpp_executable = pbadir+"PyBlastAfterglowMag"+"/src/pba.out"
        # print(os.getcwd())
        # os.chdir("../../../src/")
        # path_to_executable = "pba.out"
        if not os.path.isfile(path_to_cpp_executable):
            raise IOError("pba.out executable is not found: {}".format(path_to_cpp_executable))
        # subprocess.call(path_to_executable, input="")
        # print("{} {} {} {}".format(path_to_cpp_executable, self.workingdir, self.parfile, self.loglevel))
        # subprocess.run(path_to_cpp_executable, input=self.workingdir)
        subprocess.check_call([path_to_cpp_executable, self.workingdir, self.parfile, loglevel])
    ''' --------- magnetar ---------- '''
    def _check_if_loaded_mag_obj(self):
        if (self.fpath_mag is None):
            raise IOError("self.fpath_mag is not set")
        if (self.mag_dfile is None):
            self.mag_dfile = h5py.File(self.fpath_mag)
    ''' --------- jet (only layers) -------------- '''
    # jet dynamics
    def _ckeck_if_loaded_j_dyn_obj(self):
        if (self.fpath_grb_dyn is None):
            raise IOError("self.fpath_grb_dyn is not set")
        if (self.grb_dyn_dfile is None):
            self.grb_dyn_dfile = h5py.File(self.fpath_grb_dyn)
    # jet spectrum
    def _check_if_loaded_j_spec(self):
        if (self.fpath_grb_spec is None):
            raise IOError("self.fpath_grb_spec is not set")
        if (self.grb_spec_dfile is None):
            self.grb_spec_dfile = h5py.File(self.fpath_grb_spec)
    # jet lightcurves
    def _check_if_loaded_jet_lc(self):
        if (self.fpath_grb_light_curve is None):
            raise IOError("self.fpath_grb_light_curve is not set")
        if (self.grb_lc_dfile is None):
            self.grb_lc_dfile = h5py.File(self.fpath_grb_light_curve)
    # jet skymaps
    def _check_if_loaded_jet_skymap(self):
        if (self.fpath_grb_sky_map is None):
            raise IOError("self.fpath_grb_sky_map is not set")
        if (self.grb_skymap_dfile is None):
            self.grb_skymap_dfile = h5py.File(self.fpath_grb_sky_map)

    ''' -------- ejecta (shells and layers) -------- '''

    # ejecta dynamics
    def _ckeck_if_loaded_ej_dyn_obj(self):
        if (self.fpath_kn_dyn is None):
            raise IOError("self.fpath_kn_dyn is not set")
        if (self.kn_dyn_dfile is None):
            self.kn_dyn_dfile = h5py.File(self.fpath_kn_dyn)
    # ejecta spectrum
    def _check_if_loaded_ej_spec(self):
        if (self.fpath_kn_spec is None):
            raise IOError("self.fpath_kn_spec is not set")
        if (self.kn_spec_dfile is None):
            self.kn_spec_dfile = h5py.File(self.fpath_kn_spec)
    # ejecta lightcurves
    def _check_if_loaded_ej_lc(self):
        if (self.fpath_kn_light_curve is None):
            raise IOError("self.fpath_kn_light_curve is not set")
        if (self.kn_lc_dfile is None):
            self.kn_lc_dfile = h5py.File(self.fpath_kn_light_curve)
    # ejecta skymaps
    def _check_if_loaded_ej_skymap(self):
        if (self.fpath_kn_sky_map is None):
            raise IOError("self.fpath_kn_sky_map is not set")
        if (self.kn_skymap_dfile is None):
            try:
                self.kn_skymap_dfile = h5py.File(self.fpath_kn_sky_map)
            except OSError:
                raise OSError("failed to load the file: \n {}".format(self.fpath_kn_sky_map))
class BPA_METHODS_OLD(PBA_BASE_OLD):
    '''
        Process output_uniform_grb files: load, extract for a specific way
    '''
    def __init__(self, workingdir, readparfileforpaths=True, parfile="parfile.par"):
        super().__init__(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile)

    # magnetar
    def get_mag_obj(self):
        self._check_if_loaded_mag_obj()
        return self.mag_dfile


    # jet dynamics

    def get_jet_dyn_obj(self, ilayer=None):
        self._ckeck_if_loaded_j_dyn_obj()
        if (ilayer is None):
            return self.grb_dyn_dfile
        else:
            return self.grb_dyn_dfile["layer={}".format(ilayer)]

    def get_jet_dyn_arr(self, v_n, ilayer=0):
        self._ckeck_if_loaded_j_dyn_obj()
        dfile = self.grb_dyn_dfile
        layer = "layer={}".format(ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the jet dyn. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the jet dyn. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    # jet spectrum

    def get_jet_spec_obj(self):
        self._check_if_loaded_j_spec()
        return self.grb_spec_dfile

    def get_jet_spec_2d_arr(self, v_n="em", ilayer=0):
        dfile = self.get_jet_spec_obj()
        layer = "layer={}".format(ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the jet dyn. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the jet dyn. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    def get_jet_spec_times(self):
        dfile = self.get_jet_spec_obj()
        return np.array(dfile["times"])

    def get_jet_spec_freqs(self):
        dfile = self.get_jet_spec_obj()
        return np.array(dfile["freqs"])

    def get_jet_spec_1d_arr(self, freq=None, time=None, v_n="em", ilayer=0):
        arr = self.get_jet_spec_2d_arr(v_n=v_n, ilayer=ilayer)
        freqs = self.get_jet_spec_freqs()
        times = self.get_jet_spec_times()
        if (not freq is None):
            if (freq < freqs.min() or freq > freqs.max()):
                raise ValueError("Freq={} (jet comov cpec) is not in the limit of freqs[{}, {}]"
                                 .format(freq, freqs.min(), freqs.max()))
            if (not freq in freqs):
                idx = find_nearest_index(freqs, freq)
            else:
                idx = int(np.where(freqs == freq)[0])
            arr_ = arr[idx, :]
            assert len(arr_) == len(times)
            return arr_
        if (not times is None):
            if (time < times.min() or time > times.max()):
                raise ValueError("Time={} (jet comov cpec) is not in the limit of times[{}, {}]"
                                 .format(time, times.min(), times.max()))
            if (not time in times):
                idx = find_nearest_index(time, times)
            else:
                idx = int(np.where(time == times)[0][0])
            arr_ = arr[:, idx]
            assert len(arr_) == len(freqs)
            return arr_

    # jet lightcurves

    def get_jet_lc_obj(self):
        self._check_if_loaded_jet_lc()
        return self.grb_lc_dfile

    def get_jet_lc_times(self):
        dfile = self.get_jet_lc_obj()
        arr = np.array(dfile["times"])
        arr_u = np.unique(arr)
        if len(arr_u) == 0:
            raise ValueError("no unique times found in grb light curve")
        return arr_u

    def get_jet_lc_freqs(self):
        dfile = self.get_jet_lc_obj()
        arr = np.array(dfile["freqs"])
        arr_u = np.unique(arr)
        if len(arr_u) == 0:
            raise ValueError("no unique freqs found in grb light curve")
        return np.array(dfile["freqs"])

    def get_jet_lc_totalflux(self, freq=None):
        dfile = self.get_jet_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        times = self.get_jet_lc_times()
        freqs = self.get_jet_lc_freqs()
        fluxes = np.array(dfile["total_fluxes"])
        if freq is None:
            return fluxes
        if (not freq in freqs):
            raise ValueError("freq={} not found in freqs={}".format(freq, freqs))
        arr = fluxes[freqs==freq]
        if (len(arr)==0):
            raise ValueError("no fluxes found for freq={}".format(freq))
        return arr

    def get_jet_lc(self, freq=None, ilayer=None):
        dfile = self.get_jet_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        times = self.get_jet_lc_times()
        freqs = self.get_jet_lc_freqs()

        if (freq is None):
            # spectum
            if (ilayer is None):
                # fluxes2d = np.zeros((len(freqs), len(times)))
                # for ish in range(nshells):
                #     for il in range(nlayers):
                #         fluxes2d += np.array(dfile["shell={} layer={}".format(ish, il)][v_n])# [freq,time]
                # return fluxes2d
                fluxes2d = []
                for ifreq in freqs:
                    fluxes2d.append(self.get_jet_lc_totalflux(freq=ifreq))  # [freq,time]
                fluxes2d = np.reshape(np.array(fluxes2d), (len(freqs), len(times)))
                return fluxes2d
            elif (not ilayer is None):
                print("WARNING UNTESTED NEW OUTPUT FROM LIGHT CURVE FOR LAYER")
                arr = np.array(dfile["layer={}".format(ilayer)])
                arr = np.reshape(arr,newshape=(len(times),len(freqs)))
                return arr
            else:
                raise NameError()
        else:
            # light curves
            if (not freq in self.get_jet_lc_freqs()):
                raise ValueError("freq:{} is not in ej_lc Given:{}".format(freq, self.get_jet_lc_freqs()))
            # ifreq = find_nearest_index(self.get_jet_lc_freqs(), freq)
            if (ilayer is None):
                return self.get_jet_lc_totalflux(freq=freq)
            elif (not ilayer is None):
                print("WARNING UNTESTED NEW OUTPUT FROM LIGHT CURVE FOR LAYER")
                arr = np.array(dfile["layer={}".format(ilayer)])
                fluxes1d = arr[freqs==freq]
                # fluxes1d = np.array(dfile["layer={}".format(ilayer)][ifreq])  # [freq,time]
                return fluxes1d
            else:
                raise NameError()

    def _alpha_jet(self, freq, freqm1, ilayer, v_n):
        values_i = self.get_jet_lc(freq=freq, ilayer=ilayer, v_n=v_n)
        values_im1 = self.get_jet_lc(freq=freqm1, ilayer=ilayer, v_n=v_n)
        ffreq_i = np.full(values_i.shape, freq)
        ffreq_im1 = np.full(values_im1.shape, freqm1)
        num = np.log10(values_i) - np.log10(values_im1)
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        return values

    def get_jet_lc_spec_idx(self, freq, ishell=None, ilayer=None):
        freqs = self.get_jet_lc_freqs()
        times = self.get_jet_lc_times()
        if ((freq is None) and (ilayer is None)):
            arr = []
            for ifreq in range(1, len(freqs)):
                arr.append(self._alpha_jet(freq=freqs[ifreq], freqm1=freqs[ifreq - 1], ilayer=ilayer, v_n=None))
            return np.reshape(np.array(arr), newshape=(len(freqs) - 1, len(times)))
        else:
            idx = find_nearest_index(freqs, freq)
            return self._alpha(freq=freqs[idx], freqm1=freqs[idx - 1], ishell=ishell, ilayer=ilayer, v_n=None)

    def get_jet_lc_temp_idx(self, freq, ishell, ilayer):
        freqs = self.get_jet_lc_freqs()
        t = self.get_jet_lc_times()
        Lnu = self.get_jet_lc(freq=freq, ilayer=ilayer, v_n=None)
        val = np.zeros_like(Lnu)
        val[1:-1] = np.log10(Lnu[2:] / Lnu[:-2]) / np.log10(t[2:] / t[:-2])
        val[0] = np.log10(Lnu[1] / Lnu[0]) / np.log10(t[1] / t[0])
        val[-1] = np.log10(Lnu[-1] / Lnu[-2]) / np.log10(t[-1] / t[-2])
        return val

    # jet skymaps

    def get_jet_skymap_obj(self):
        self._check_if_loaded_jet_skymap()
        return self.grb_skymap_dfile

    def get_jet_skymap_times(self):
        self._check_if_loaded_jet_skymap()
        return np.array(self.grb_skymap_dfile["times"])

    def get_jet_skymap_freqs(self):
        self._check_if_loaded_jet_skymap()
        return np.array(self.grb_skymap_dfile["freqs"])

    def get_jet_skymap_cm(self, all_xrs, all_yrs, all_zz):
        dfile = self.get_jet_skymap_obj()
        xc_m, yc_m = compute_position_of_the_flux_centroid(all_xrs, all_yrs, all_zz, float(dfile.attrs["d_l"]))
        return (xc_m, yc_m)

    # def get_jet_skymap(self, time=None, freq=None, verbose=False, remove_mu=False, int_or_hist="int"):
    #     ish = 0
    #     nx = 2000
    #     ny = 1000
    #     min_mu = 1e-5  # limit for mu that if too small blow up the intensity for some reason...
    #     min_mu_frac = 0.85
    #     # limit_data = False
    #     times = self.get_jet_skymap_times()
    #     freqs = self.get_jet_skymap_freqs()
    #     dfile = self.get_jet_skymap_obj()
    #     # nshells = int(dfile.attrs["nshells"])
    #     if ((not time is None) and (not time in times)):
    #         raise ValueError("time={} is not in the list for skypams={}".format(time, times))
    #     if ((not freq is None) and (not freq in freqs)):
    #         raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
    #     if ((not time is None) and (not freq is None)):
    #         print(dfile.keys())
    #         key = "time={:.4e} freq={:.4e}".format(time, freq)
    #         try:
    #             ddfile = dfile[key]
    #         except KeyError:
    #             try:
    #                 ddfile = dfile[key.replace(".", ",")]
    #             except KeyError:
    #                     raise KeyError("Failed to find ky for jet skymap : {} Avaialble:\n{}".format(key, dfile.keys()))
    #         r_i = np.array(ddfile["r"][ish])
    #         mu_i = np.array(ddfile["mu"][ish])
    #         xrs_i = np.array(ddfile["xrs"][ish])
    #         yrs_i = np.array(ddfile["yrs"][ish])
    #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
    #         gam_i = np.array(ddfile["gamma"][ish])
    #         B_i = np.array(ddfile["B"][ish])
    #         tb_i = np.array(ddfile["tburst"][ish])
    #         theta_i = np.array(ddfile["theta"][ish])
    #         theta_ij = np.array(ddfile["theta_j"][ish])
    #         theta0_i = np.array(ddfile["theta0"][ish])
    #         phi_i = np.array(ddfile["phi"][ish])
    #
    #         # remove nans
    #         # int_i = _int_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #         # yrs_i = _yrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #         # xrs_i = _xrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #
    #
    #         # print(np.arccos(mu_i[np.abs(mu_i) < 1e-4]))
    #         print("Min={} max={}".format(np.min(np.abs(mu_i)), np.max(np.abs(mu_i))))
    #
    #         # plt.plot(mu_i, theta0_i, marker='.', ls='none')
    #         # plt.semilogy(mu_i, int_i, marker='.', ls='none')
    #         # plt.show()
    #
    #         if verbose:
    #
    #
    #             _mu = np.arccos(mu_i)
    #             x = np.exp(1 + abs(mu_i).max() / -abs(mu_i))
    #             idx1 = abs(mu_i) < min_mu
    #             int_i2 = np.copy(int_i)  * abs(mu_i) #/ (1. - 2.*np.sin(_mu/2)) #/ mu_i# (1.-mu_i*get_beta(gam_i))/(mu_i-get_beta(gam_i))#np.abs(mu_i) / mu_i**2
    #             #fac = np.exp(1 / (1 / abs(mu_i))-1)
    #             # int_i2[abs(mu_i) < 1e-3] *= 1e-3
    #             #int_i2[abs(mu_i) < 1e-1] *= abs(mu_i)[abs(mu_i) < 1e-1]
    #             # int_i2[abs(mu_i) < 1e-5] *= 1e-5
    #             # int_i2[abs(mu_i) < 1e-6] *= 1e-6
    #             #idx1 = phi_i == 0.
    #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
    #             # int_i[idx1] =
    #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
    #             # idx = int_i > int_i.max() * min_mu_frac
    #             idx = idx1
    #             if verbose:
    #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
    #                     print('\t', mu_i[idx1])
    #                     print('\t', mu_i[idx])
    #                     # exit(1)
    #             if verbose:
    #                 # if len(idx) < 2:
    #                 #     print("No excess = {}".format(ii))
    #                 print("len(gam[idx])={}".format(gam_i[idx]))
    #                 print("Gamma={}".format(gam_i[idx]))
    #                 print("beta={}".format(get_beta(gam_i[idx])))
    #                 print("mu_i={}".format(mu_i[idx]))
    #                 print("r_i={}".format(r_i[idx]))
    #                 print("B={}".format(B_i[idx]))
    #                 print("tb_i={}".format(tb_i[idx]))
    #             print(np.abs(mu_i).min())
    #             print(theta_i[np.abs(mu_i) < 1e-5])
    #             print(phi_i[np.abs(mu_i) < 1e-5])
    #
    #             fig, axes = plt.subplots(figsize=(4.6+2, 3.2), ncols=4, nrows=1)
    #             ax = axes[0]
    #             ax.plot(mu_i, int_i / int_i.max(), '.')
    #             ax.plot(mu_i, int_i2 / int_i.max(), 'o', zorder=-1)
    #             if len(mu_i[idx]) > 0: ax.plot(mu_i[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[1]
    #             ax.plot(np.sin(phi_i), int_i / int_i.max(), '.')
    #             ax.plot(np.sin(phi_i), int_i2 / int_i.max(), 'o', zorder=-1)
    #             if len(mu_i[idx]) > 0: ax.plot(np.sin(phi_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\phi$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[2]
    #             ax.plot(np.sin(theta_i), int_i / int_i.max(), '.')
    #             ax.plot(np.sin(theta_i), int_i2 / int_i2.max(), 'o')
    #             if len(mu_i[idx]) > 0: ax.plot(np.sin(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\sin\theta$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[3]
    #             ax.plot(np.cos(theta_i), int_i / int_i.max(), '.')
    #             ax.plot(np.cos(theta_i), int_i2 / int_i.max(), 'o')
    #             if len(mu_i[idx]) > 0: ax.plot(np.cos(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\cos\theta$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             plt.tight_layout()
    #             plt.show()
    #
    #         # plt.plot(theta0_i, mu_i, marker='.', ls='none')
    #         # plt.show()
    #
    #
    #         if remove_mu:
    #             int_i *= np.abs(mu_i)
    #
    #         # int_i[~np.isfinite(int_i)] = 0.0
    #         # int_i[idx] = 0.
    #         # int_i = int_i2
    #
    #         if int_or_hist == "int":
    #             nx = np.complex(0, nx)
    #             ny = np.complex(0, ny)
    #             grid_x, grid_y = np.mgrid[xrs_i.min()*1.2:xrs_i.max()*1.2:nx, yrs_i.min()*1.2:yrs_i.max()*1.2:ny]
    #             i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #         elif int_or_hist == "both":
    #
    #             nx = 800
    #             ny = 400
    #
    #             nx = np.complex(0, nx)
    #             ny = np.complex(0, ny)
    #             grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx,
    #                              yrs_i.min():yrs_i.max():ny]
    #             i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #
    #
    #             nx = np.complex(0, nx + 1)
    #             ny = np.complex(0, ny + 1)
    #             edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #             edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #             #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             i_zz, _ = np.histogramdd(tuple([grid_x.flatten(), grid_y.flatten()]), bins=tuple([edges_x, edges_y]), weights=i_zz.T.flatten())
    #             grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #             grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #             print(i_zz.shape)
    #         else:
    #
    #
    #             nx = 200
    #             ny = 100
    #             nx = np.complex(0, nx+1)
    #             ny = np.complex(0, ny+1)
    #             edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #             edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #             #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #             grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #             grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #             print(i_zz.shape)
    #
    #
    #         # i_zz = ndimage.uniform_filter(i_zz, )
    #         # i_zz = ndimage.filters.gaussian_filter(i_zz, [10,10], mode='reflect')
    #
    #         return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #
    #         #
    #         # else:
    #         #     # assess what is the combined grid extend (for all images)
    #         #     xmin_neg, xmin_pos, xmax = [], [], []
    #         #     ymin_neg, ymin_pos, ymax = [], [], []
    #         #     i_min, i_max = [], []
    #         #     i_shells = []
    #         #     for ish in range(nshells):
    #         #         xrs_i = np.array(ddfile["xrs"][ish])
    #         #         yrs_i = np.array(ddfile["yrs"][ish])
    #         #         int_i = np.array(ddfile["intensity"][ish])
    #         #         if (np.sum(int_i) == 0):
    #         #             continue
    #         #         i_shells.append(ish)
    #         #         xmin_neg.append(xrs_i.min())
    #         #         xmin_pos.append(xrs_i[xrs_i > 0].min())
    #         #         xmax.append(xrs_i.max())
    #         #         ymin_neg.append(yrs_i.min())
    #         #         ymin_pos.append(yrs_i[yrs_i > 0].min())
    #         #         ymax.append(yrs_i.max())
    #         #         i_min.append(int_i.min())
    #         #         i_max.append(int_i.max())
    #         #     xmin_neg = np.array(xmin_neg)
    #         #     xmin_pos = np.array(xmin_pos)
    #         #     xmax = np.array(xmax)
    #         #     ymin_neg = np.array(ymin_neg)
    #         #     ymin_pos = np.array(ymin_pos)
    #         #     ymax = np.array(ymax)
    #         #     i_min = np.array(i_min)
    #         #     i_max = np.array(i_max)
    #         #     if verbose:
    #         #         for i in range(len(i_shells)):
    #         #             print(
    #         #                 "\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
    #         #                                                                                 xmax[i]))
    #         #             print(
    #         #                 "\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
    #         #                                                                                 ymax[i]))
    #         #         print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
    #         #         print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
    #         #         print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
    #         #     x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
    #         #                              np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
    #         #     y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
    #         #                              np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
    #         #     if verbose:
    #         #         print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
    #         #                                                                           x_grid[x_grid > 0].min(),
    #         #                                                                           x_grid.max()))
    #         #         print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
    #         #                                                                           y_grid[y_grid > 0].min(),
    #         #                                                                           y_grid.max()))
    #         #     xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    #         #     zz = np.zeros_like((xx_grid))
    #         #     # interpolate onto the grid that covers all images
    #         #     all_xrs, all_yrs, all_zz = [], [], []
    #         #     for ii, ish in enumerate(i_shells):  # range(len(i_shells))
    #         #         if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
    #         #         r_i = np.array(ddfile["r"][ish])
    #         #         mu_i = np.array(ddfile["mu"][ish])
    #         #         xrs_i = np.array(ddfile["xrs"][ish])
    #         #         yrs_i = np.array(ddfile["yrs"][ish])
    #         #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
    #         #         gam_i = np.array(ddfile["gamma"][ish])
    #         #         B_i = np.array(ddfile["B"][ish])
    #         #         tb_i = np.array(ddfile["tburst"][ish])
    #         #         theta_i = np.array(ddfile["theta"][ish])
    #         #         if limit_data:
    #         #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
    #         #             # int_i[idx1] =
    #         #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
    #         #             idx = int_i > int_i.max() * min_mu_frac
    #         #             idx = idx1
    #         #             if verbose:
    #         #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
    #         #                     print('\t', mu_i[idx1])
    #         #                     print('\t', mu_i[idx])
    #         #                     # exit(1)
    #         #             if verbose:
    #         #                 if len(idx) < 2:
    #         #                     print("No excess = {}".format(ii))
    #         #                 print("len(gam[idx])={}".format(gam_i[idx]))
    #         #                 print("Gamma={}".format(gam_i[idx]))
    #         #                 print("beta={}".format(get_beta(gam_i[idx])))
    #         #                 print("mu_i={}".format(mu_i[idx]))
    #         #                 print("r_i={}".format(r_i[idx]))
    #         #                 print("B={}".format(B_i[idx]))
    #         #                 print("tb_i={}".format(tb_i[idx]))
    #         #             # print(gam_i[idx], Beta(gam_i[idx]))
    #         #
    #         #             # int_i[mu_i < 1e-3] = 0.
    #         #             # int_i2 = int_i * abs(mu_i)# * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
    #         #             # plt.figure(figsize=(4.6,3.2))
    #         #             # plt.semilogy(mu_i,int_i/int_i.max(), '.')
    #         #             # plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
    #         #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #         #             # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #         #             # plt.ylabel(r"$I/I_{\rm max}$")
    #         #             # plt.tight_layout()
    #         #             # plt.show()
    #         #             int_i[idx] = 0.
    #         #         # int_i /= abs(mu_i)[abs(mu_i)>0.01]
    #         #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #         #         zz += i_zz
    #         #         all_xrs.append(xrs_i)
    #         #         all_yrs.append(yrs_i)
    #         #         all_zz.append(int_i)
    #         #     if verbose:
    #         #         print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
    #         #               .format(zz.min(), zz.max(), np.sum(zz), float(
    #         #             np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
    #         #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
    #     else:
    #         raise NotImplementedError("Not implemented")
    def get_jet_skymap(self, time=None, freq=None, verbose=False, remove_mu=False, plots=False, renormalize=True):
        ish = 0
        # nx = 1000
        # ny = 500
        min_mu = 1e-5  # limit for mu that if too small blow up the intensity for some reason...
        # min_mu_frac = 0.85
        # limit_data = False
        d_l = float(self.get_jet_skymap_obj().attrs["d_l"])
        # theta_obs = float(self.get_jet_skymap_obj().attrs["thetaObs"])
        times = self.get_jet_skymap_times()
        freqs = self.get_jet_skymap_freqs()
        dfile = self.get_jet_skymap_obj()
        # nshells = int(dfile.attrs["nshells"])
        if ((not time is None) and (not time in times)):
            raise ValueError("time={} ({} days) is not in the list for skypams={} ({} days)"
                             .format(time, time/cgs.day, times, times/cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
        if ((not time is None) and (not freq is None)):
            # print(dfile.keys())
            key = "time={:.4e} freq={:.4e}".format(time, freq)
            try:
                ddfile = dfile[key]
            except KeyError:
                try:
                    ddfile = dfile[key.replace(".", ",")]
                except KeyError:
                    raise KeyError("Failed to find ky for jet skymap : {} Avaialble:\n{}".format(key, dfile.keys()))
            # r_i = np.array(ddfile["r"][ish])
            mu_i = np.array(ddfile["mu"][ish])
            xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
            yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
            int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
            # gam_i = np.array(ddfile["gamma"][ish])
            # B_i = np.array(ddfile["B"][ish])
            # tb_i = np.array(ddfile["tburst"][ish])
            # theta_i = np.array(ddfile["theta"][ish])
            # theta_ij = np.array(ddfile["theta_j"][ish])
            # theta0_i = np.array(ddfile["theta0"][ish])
            # phi_i = np.array(ddfile["phi"][ish])
            # tt_i = np.array(ddfile["tt"][ish])

            # xrs =  -1. * np.cos(theta_obs) * np.sin(theta_i) * np.sin(phi_i) + np.sin(theta_obs) * np.cos(theta_i)
            # yrs = np.sin(theta_i) * np.cos(phi_i)

            # xrs_i = xrs * r_i * cgs.rad2mas / d_l
            # yrs_i = yrs * r_i * cgs.rad2mas / d_l

            # xrs_i[np.isnan(int_i)] = 0.
            # yrs_i[np.isnan(int_i)] = 0.
            # int_i[np.isnan(int_i)] = 0.
            # int_i[phi_i < np.pi ] = 0.
            # yrs_i[phi_i < np.pi] = 0.
            # xrs_i[phi_i < np.pi] = 0.

            if plots:
                layers = int(self.get_jet_dyn_obj().attrs["nlayers"])
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=0) / cgs.day, self.get_jet_dyn_arr("ctheta", ilayer=0),
                             label="nl={}".format(0))
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=layers - 1) / cgs.day,
                             self.get_jet_dyn_arr("ctheta", ilayer=layers - 1), label="nl={}".format(layers - 1))
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=layers - int(layers / 2)) / cgs.day,
                             self.get_jet_dyn_arr("ctheta", ilayer=layers - int(layers / 2)),
                             label="nl={}".format(layers - int(layers / 2)))
                plt.xlabel("Time [days]")
                plt.grid()
                plt.legend()
                plt.show()

                fig = plt.figure()
                cmap = cm.get_cmap('viridis')
                my_norm = LogNorm(int_i.max() * 1e-2, int_i.max())
                ax = fig.add_subplot(projection='3d')
                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
                ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(r_i).flatten(), c=cmap(my_norm(int_i.flatten())))
                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
                ax.set_xlabel('X Label')
                ax.set_ylabel('Y Label')
                ax.set_zlabel('I Label')
                plt.show()

            # remove nans
            # int_i = _int_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
            # yrs_i = _yrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
            # xrs_i = _xrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]

            # print(np.arccos(mu_i[np.abs(mu_i) < 1e-4]))
            # print("Jet skymap min={} max={}".format(np.min(np.abs(mu_i)), np.max(np.abs(mu_i))))

            # plt.plot(mu_i, theta0_i, marker='.', ls='none')
            # plt.semilogy(mu_i, int_i, marker='.', ls='none')
            # plt.show()

            if verbose:

                _mu = np.arccos(mu_i)
                x = np.exp(1 + abs(mu_i).max() / -abs(mu_i))
                idx1 = abs(mu_i) < min_mu
                int_i2 = np.copy(int_i) * abs(
                    mu_i)  # / (1. - 2.*np.sin(_mu/2)) #/ mu_i# (1.-mu_i*get_beta(gam_i))/(mu_i-get_beta(gam_i))#np.abs(mu_i) / mu_i**2
                # fac = np.exp(1 / (1 / abs(mu_i))-1)
                # int_i2[abs(mu_i) < 1e-3] *= 1e-3
                # int_i2[abs(mu_i) < 1e-1] *= abs(mu_i)[abs(mu_i) < 1e-1]
                # int_i2[abs(mu_i) < 1e-5] *= 1e-5
                # int_i2[abs(mu_i) < 1e-6] *= 1e-6
                # idx1 = phi_i == 0.
                idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
                # int_i[idx1] =
                # print("len(gam[idx1])={}".format(gam_i[idx1]))
                # idx = int_i > int_i.max() * min_mu_frac
                idx = idx1
                if verbose:
                    if (len(mu_i[idx1]) != len(mu_i[idx])):
                        print('\t', mu_i[idx1])
                        print('\t', mu_i[idx])
                        # exit(1)
                if verbose:
                    # if len(idx) < 2:
                    #     print("No excess = {}".format(ii))
                    print("len(gam[idx])={}".format(gam_i[idx]))
                    print("Gamma={}".format(gam_i[idx]))
                    print("beta={}".format(get_beta(gam_i[idx])))
                    print("mu_i={}".format(mu_i[idx]))
                    print("r_i={}".format(r_i[idx]))
                    print("B={}".format(B_i[idx]))
                    print("tb_i={}".format(tb_i[idx]))
                print(np.abs(mu_i).min())
                print(theta_i[np.abs(mu_i) < 1e-5])
                print(phi_i[np.abs(mu_i) < 1e-5])

                fig, axes = plt.subplots(figsize=(4.6 + 2, 3.2), ncols=4, nrows=1)
                ax = axes[0]
                ax.plot(mu_i, int_i / int_i.max(), '.')
                ax.plot(mu_i, int_i2 / int_i.max(), 'o', zorder=-1)
                if len(mu_i[idx]) > 0: ax.plot(mu_i[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(
                    r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[1]
                ax.plot(np.sin(phi_i), int_i / int_i.max(), '.')
                ax.plot(np.sin(phi_i), int_i2 / int_i.max(), 'o', zorder=-1)
                if len(mu_i[idx]) > 0: ax.plot(np.sin(phi_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\phi$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[2]
                ax.plot(np.sin(theta_i), int_i / int_i.max(), '.')
                ax.plot(np.sin(theta_i), int_i2 / int_i2.max(), 'o')
                if len(mu_i[idx]) > 0: ax.plot(np.sin(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\sin\theta$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[3]
                ax.plot(np.cos(theta_i), int_i / int_i.max(), '.')
                ax.plot(np.cos(theta_i), int_i2 / int_i.max(), 'o')
                if len(mu_i[idx]) > 0: ax.plot(np.cos(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\cos\theta$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                plt.tight_layout()
                plt.show()

            # plt.plot(theta0_i, mu_i, marker='.', ls='none')
            # plt.show()

            if remove_mu:
                print("Removing 'mu' from jet skymap")
                int_i *= np.abs(mu_i)

            # int_i[~np.isfinite(int_i)] = 0.0
            # int_i[idx] = 0.
            # int_i = int_i2

            if renormalize:
                print("Renormalizing jet skymap to total flux")
                # Problem: 'all_fluxes_jet' decreased for increasing nlayers; Solution: Below
                all_fluxes_jet = np.array(int_i)
                fnus = self.get_jet_skymap_totfluxes(freq=freq)
                fnu = fnus[find_nearest_index(self.get_jet_skymap_times(), time)]
                delta_x = xrs_i.max() - xrs_i.min()
                delta_y = yrs_i.max() - yrs_i.min()
                dfnu = fnu / (delta_x * delta_y)
                if verbose:
                    print("Jet")
                    print("\tfnu = {:.2e} ".format(fnu))
                    print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(xrs_i).min(), np.array(xrs_i).max()))
                    print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(yrs_i).min(), np.array(yrs_i).max()))
                    print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                    print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                int_i = (all_fluxes_jet / all_fluxes_jet.max()) * dfnu

            return (xrs_i, yrs_i, int_i)

            #
            # else:
            #     # assess what is the combined grid extend (for all images)
            #     xmin_neg, xmin_pos, xmax = [], [], []
            #     ymin_neg, ymin_pos, ymax = [], [], []
            #     i_min, i_max = [], []
            #     i_shells = []
            #     for ish in range(nshells):
            #         xrs_i = np.array(ddfile["xrs"][ish])
            #         yrs_i = np.array(ddfile["yrs"][ish])
            #         int_i = np.array(ddfile["intensity"][ish])
            #         if (np.sum(int_i) == 0):
            #             continue
            #         i_shells.append(ish)
            #         xmin_neg.append(xrs_i.min())
            #         xmin_pos.append(xrs_i[xrs_i > 0].min())
            #         xmax.append(xrs_i.max())
            #         ymin_neg.append(yrs_i.min())
            #         ymin_pos.append(yrs_i[yrs_i > 0].min())
            #         ymax.append(yrs_i.max())
            #         i_min.append(int_i.min())
            #         i_max.append(int_i.max())
            #     xmin_neg = np.array(xmin_neg)
            #     xmin_pos = np.array(xmin_pos)
            #     xmax = np.array(xmax)
            #     ymin_neg = np.array(ymin_neg)
            #     ymin_pos = np.array(ymin_pos)
            #     ymax = np.array(ymax)
            #     i_min = np.array(i_min)
            #     i_max = np.array(i_max)
            #     if verbose:
            #         for i in range(len(i_shells)):
            #             print(
            #                 "\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
            #                                                                                 xmax[i]))
            #             print(
            #                 "\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
            #                                                                                 ymax[i]))
            #         print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
            #         print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
            #         print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
            #     x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
            #                              np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
            #     y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
            #                              np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
            #     if verbose:
            #         print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
            #                                                                           x_grid[x_grid > 0].min(),
            #                                                                           x_grid.max()))
            #         print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
            #                                                                           y_grid[y_grid > 0].min(),
            #                                                                           y_grid.max()))
            #     xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
            #     zz = np.zeros_like((xx_grid))
            #     # interpolate onto the grid that covers all images
            #     all_xrs, all_yrs, all_zz = [], [], []
            #     for ii, ish in enumerate(i_shells):  # range(len(i_shells))
            #         if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
            #         r_i = np.array(ddfile["r"][ish])
            #         mu_i = np.array(ddfile["mu"][ish])
            #         xrs_i = np.array(ddfile["xrs"][ish])
            #         yrs_i = np.array(ddfile["yrs"][ish])
            #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
            #         gam_i = np.array(ddfile["gamma"][ish])
            #         B_i = np.array(ddfile["B"][ish])
            #         tb_i = np.array(ddfile["tburst"][ish])
            #         theta_i = np.array(ddfile["theta"][ish])
            #         if limit_data:
            #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
            #             # int_i[idx1] =
            #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
            #             idx = int_i > int_i.max() * min_mu_frac
            #             idx = idx1
            #             if verbose:
            #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
            #                     print('\t', mu_i[idx1])
            #                     print('\t', mu_i[idx])
            #                     # exit(1)
            #             if verbose:
            #                 if len(idx) < 2:
            #                     print("No excess = {}".format(ii))
            #                 print("len(gam[idx])={}".format(gam_i[idx]))
            #                 print("Gamma={}".format(gam_i[idx]))
            #                 print("beta={}".format(get_beta(gam_i[idx])))
            #                 print("mu_i={}".format(mu_i[idx]))
            #                 print("r_i={}".format(r_i[idx]))
            #                 print("B={}".format(B_i[idx]))
            #                 print("tb_i={}".format(tb_i[idx]))
            #             # print(gam_i[idx], Beta(gam_i[idx]))
            #
            #             # int_i[mu_i < 1e-3] = 0.
            #             # int_i2 = int_i * abs(mu_i)# * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
            #             # plt.figure(figsize=(4.6,3.2))
            #             # plt.semilogy(mu_i,int_i/int_i.max(), '.')
            #             # plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
            #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
            #             # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
            #             # plt.ylabel(r"$I/I_{\rm max}$")
            #             # plt.tight_layout()
            #             # plt.show()
            #             int_i[idx] = 0.
            #         # int_i /= abs(mu_i)[abs(mu_i)>0.01]
            #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
            #         zz += i_zz
            #         all_xrs.append(xrs_i)
            #         all_yrs.append(yrs_i)
            #         all_zz.append(int_i)
            #     if verbose:
            #         print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
            #               .format(zz.min(), zz.max(), np.sum(zz), float(
            #             np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
            #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
        else:
            raise NotImplementedError("Not implemented")

    def get_jet_skymap_totfluxes(self, freq, time=None):

        self._check_if_loaded_jet_skymap()
        if time is None:
            return np.array(self.grb_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
        else:
            arr = np.array(self.grb_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
            times = self.get_jet_skymap_times()
            if len(times) > 1:
                assert self.get_jet_skymap_times().min() < time < self.get_jet_skymap_times().max()
                val = arr[find_nearest_index(self.get_jet_skymap_totfluxes(freq=freq, time=None), time)]
                return val
            else:
                return arr[0]

    ''' -------- ejecta (shells and layers) -------- '''

    # ejecta dynamics

    def get_ej_dyn_obj(self):
        self._ckeck_if_loaded_ej_dyn_obj()
        return self.kn_dyn_dfile

    def get_ej_dyn_arr(self, v_n, ishell=None, ilayer=None):
        self._ckeck_if_loaded_ej_dyn_obj()
        obj_dyn = self.get_ej_dyn_obj()
        nlayers = int(obj_dyn.attrs["nlayers"])
        nshells = int(obj_dyn.attrs["nshells"])
        dfile = self.kn_dyn_dfile
        if ((ishell is None) and (ilayer is None)):
            raise ValueError("both ishell and ilayer cannot be None for dynamical data")
        elif ((not ishell is None) and (ilayer is None)):
            arr = []
            for il in range(nlayers):
                arr.append(self.get_ej_dyn_arr(v_n, ishell=ishell, ilayer=il))
            arr = np.reshape(np.array(arr), newshape=(nlayers, len(arr[0])))
            return arr
        elif ((ishell is None) and (not ilayer is None)):
            arr = []
            for ish in range(nshells):
                arr.append(self.get_ej_dyn_arr(v_n, ishell=ish, ilayer=ilayer))
            arr = np.reshape(np.array(arr), newshape=(nshells, len(arr[0])))
            return arr
        elif ((not ishell is None) and (not ilayer is None)):
            layer = "shell={} layer={}".format(ishell, ilayer)
            if (not layer in list(dfile.keys())):
                raise NameError("Layer {} (aka '{}') is not in the ejecta dyn. file.\n Available: {}"
                                .format(ilayer, layer, dfile.keys()))
            if (not v_n in dfile[layer].keys()):
                raise NameError("v_n {} is not in the ejecta dyn. dfile[{}].keys() \n Avaialble: {}"
                                .format(v_n, layer, dfile[layer].keys()))
            return np.array(dfile[layer][v_n])
        else:
            raise NameError()

    def get_ej_dyn_1d_arr_layers(self, v_n="em", idx=0, ishell=None, ilayer=None):
        obj_dyn = self.get_ej_dyn_obj()
        nlayers = int(obj_dyn.attrs["nlayers"])
        nshells = int(obj_dyn.attrs["nshells"])
        x_arr = []
        if (not ishell is None and ilayer is None):
            for il in range(nlayers):
                x_arr.append(self.get_ej_dyn_arr(v_n, ishell=ishell, ilayer=il)[idx])
        else:
            for ish in range(nshells):
                x_arr.append(self.get_ej_dyn_arr(v_n, ishell=ish, ilayer=ilayer)[idx])
        return np.array(x_arr)

    # ejecta spectrum

    def get_ej_spec_obj(self):
        self._check_if_loaded_ej_spec()
        return self.kn_spec_dfile

    def get_ej_spec_2d_arr(self, v_n="em", ishell=0, ilayer=0):
        dfile = self.get_ej_spec_obj()
        layer = "shell={} layer={}".format(ishell, ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the ejecta comov.spec. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the ejecta comov.spec. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    def get_ej_spec_times(self):
        dfile = self.get_ej_spec_obj()
        return np.array(dfile["times"])

    def get_ej_spec_freqs(self):
        dfile = self.get_ej_spec_obj()
        return np.array(dfile["freqs"])

    def get_ej_spec_1d_arr(self, freq=None, time=None, v_n="em", ishell=0, ilayer=0):
        arr = self.get_ej_spec_2d_arr(v_n=v_n, ishell=ishell, ilayer=ilayer)
        freqs = self.get_ej_spec_freqs()
        times = self.get_ej_spec_times()
        if (not freq is None):
            if (freq < freqs.min() or freq > freqs.max()):
                raise ValueError("Freq={} (jet comov cpec) is not in the limit of freqs[{}, {}]"
                                 .format(freq, freqs.min(), freqs.max()))
            if (not freq in freqs):
                idx = find_nearest_index(freqs, freq)
            else:
                idx = int(np.where(freqs == freq)[0])
            arr_ = arr[idx, :]
            assert len(arr_) == len(times)
            return arr_
        if (not times is None):
            if (time < times.min() or time > times.max()):
                raise ValueError("Time={} (jet comov cpec) is not in the limit of times[{}, {}]"
                                 .format(time, times.min(), times.max()))
            if (not time in times):
                idx = find_nearest_index(time, times)
            else:
                idx = int(np.where(time == times)[0][0])
            arr_ = arr[:, idx]
            assert len(arr_) == len(freqs)
            return arr_

    def get_ej_spec_2d_arr_layers(self, freq=None, time=None, v_n="em", ishell=None, ilayer=None):
        obj_dyn = self.get_ej_dyn_obj()
        arr = []
        x_arr = []
        # x_arr = self.get_ej_dyn_arr("ctheta",ishell=ishell,ilayer=0)
        if (not ishell is None and ilayer is None):
            y_arr = self.get_ej_dyn_arr("R", ishell=ishell, ilayer=0)
            # freqs = self.get_ej_spec_freqs()
            # times = self.get_ej_spec_times()
            nlayers = int(obj_dyn.attrs["nlayers"])
            for il in range(nlayers):
                x_arr.append(self.get_ej_dyn_arr("ctheta", ishell=ishell, ilayer=il)[0])
                arr.append(self.get_ej_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ishell, ilayer=il))
            x_arr = np.array(x_arr)
            arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
        else:
            y_arr = self.get_ej_dyn_arr("R", ishell=0, ilayer=ilayer)
            nshells = int(obj_dyn.attrs["nshells"])
            for ish in range(nshells):
                x_arr.append(self.get_ej_dyn_arr("Gamma", ishell=ish, ilayer=ilayer)[0])
                arr.append(self.get_ej_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ish, ilayer=ilayer))
            x_arr = np.array(x_arr)
            arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
        return arr

    # ejecta lightcurves

    def get_ej_lc_obj(self):
        self._check_if_loaded_ej_lc()
        return self.kn_lc_dfile

    def get_ej_lc_times(self):
        dfile = self.get_ej_lc_obj()
        arr = np.array(dfile["times"])
        arr_u = np.unique(arr)
        if len(arr_u) == 0:
            raise ValueError("no unique times found in kn light curve")
        return arr_u

    def get_ej_lc_freqs(self):
        dfile = self.get_ej_lc_obj()
        arr = np.array(dfile["freqs"])
        arr_u = np.unique(arr)
        if len(arr_u) == 0:
            raise ValueError("no unique freqs found in kn light curve")
        return np.array(dfile["freqs"])

    def get_ej_lc_totalflux(self, freq=None):
        dfile = self.get_ej_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        nshells = int(dfile.attrs["nshells"])
        times = self.get_ej_lc_times()
        freqs = self.get_ej_lc_freqs()
        fluxes = np.array(dfile["total_fluxes"])
        if (freq is None):
            return fluxes
        if (not freq in freqs):
            raise ValueError("freq={} not found in freqs={}".format(freq, freqs))
        arr = fluxes[freqs==freq]
        if (len(arr)==0):
            raise ValueError("no fluxes found for freq={}".format(freq))
        return arr

        #
        # nlayers = int(dfile.attrs["nlayers"])
        # nshells = int(dfile.attrs["nshells"])
        # times = self.get_ej_lc_times()
        # freqs = self.get_ej_lc_freqs()
        # try:
        #     key = str("totalflux at freq={:.4e}".format(freq)).replace('.', ',')
        #     if not key in dfile.keys():
        #         raise NameError("Not found: {} among keys:{}".format(key, [key for key in dfile.keys() if key.__contains__("totalflux")]))
        # except NameError:
        #     key = str("totalflux at freq={:.4e}".format(freq))
        #     if not key in dfile.keys():
        #         raise NameError("Not found for ej. lightcurve: {} among keys:{} \n ejecta_prefix:{}"
        #                         .format(key, [key for key in dfile.keys() if key.__contains__("totalflux")], self.ejecta_prefix))
        # except:
        #     raise NameError()
        # return np.array(dfile[key])

    def get_ej_lc(self, freq=None, ishell=None, ilayer=None):
        dfile = self.get_ej_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        nshells = int(dfile.attrs["nshells"])
        times = self.get_ej_lc_times()
        freqs = self.get_ej_lc_freqs()

        if (freq is None):
            # spectum
            if ((ishell is None) and (ilayer is None)):
                fluxes2d = []
                for ifreq in freqs:
                    fluxes2d.append(self.get_ej_lc_totalflux(freq=ifreq))  # [freq,time]
                fluxes2d = np.reshape(np.array(fluxes2d), (len(freqs), len(times)))
                return fluxes2d
            elif ((ishell is None) and (not ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                for ish in range(nshells):
                    arr = np.array(dfile["shell={} layer={}".format(ish, ilayer)])
                    arr = np.reshape(arr, (len(times),len(freqs)))
                    fluxes2d += arr  # [freq,time]
                return fluxes2d
            elif ((not ishell is None) and (ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                for il in range(nlayers):
                    arr = np.array(dfile["shell={} layer={}".format(ishell, il)])
                    arr = np.reshape(arr, (len(times),len(freqs)))
                    fluxes2d += arr  # [freq,time]
                return fluxes2d
            elif ((not ishell is None) and (not ilayer is None)):
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                arr = np.array(dfile["shell={} layer={}".format(ishell, ilayer)])
                arr = np.reshape(arr, (len(times),len(freqs)))
                fluxes2d = arr  # [freq,time]
                return fluxes2d
            else:
                raise NameError()
        else:
            # light curves
            if (not freq in self.get_ej_lc_freqs()):
                raise ValueError("freq:{} is not in ej_lc Given:{}".format(freq, self.get_ej_lc_freqs()))
            # ifreq = find_nearest_index(self.get_ej_lc_freqs(), freq)
            if ((ishell is None) and (ilayer is None)):
                return self.get_ej_lc_totalflux(freq=freq)
            elif ((ishell is None) and (not ilayer is None)):
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                if (ilayer > nlayers - 1):
                    raise ValueError("Layer={} is > nlayers={}".format(ilayer, nlayers - 1))
                # fluxes1d = np.zeros_like(times)
                # for ish in range(nshells):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ish, ilayer)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for ish in range(nshells):
                    arr = np.array(dfile["shell={} layer={}".format(ish, ilayer)])
                    arr = arr[freqs==freq]
                    fluxes2d.append(arr)  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nshells, len(times)))
            elif ((not ishell is None) and (ilayer is None)):
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                if (ishell > nshells - 1):
                    raise ValueError("Shell={} is > nshell={}".format(ishell, nshells - 1))
                # fluxes1d = np.zeros_like(times)
                # for il in range(nlayers):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ishell, il)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for il in range(nlayers):
                    arr = np.array(dfile["shell={} layer={}".format(ishell, il)])
                    arr = arr[freqs==freq]
                    fluxes2d.append(arr)  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nlayers, len(times)))
            elif ((not ishell is None) and (not ilayer is None)):
                print("UNTESTED PART OF CDOE AFTER NEW LIGHT CURVE OUTPUT")
                arr = np.array(dfile["shell={} layer={}".format(ishell, ilayer)])
                arr = arr[freqs==freq]
                fluxes1d = arr  # [freq,time]
                return fluxes1d
            else:
                raise NameError()

    def _alpha(self, freq, freqm1, ishell, ilayer, v_n):
        values_i = self.get_ej_lc(freq=freq, ishell=ishell, ilayer=ilayer, v_n=v_n)
        values_im1 = self.get_ej_lc(freq=freqm1, ishell=ishell, ilayer=ilayer, v_n=v_n)
        ffreq_i = np.full(values_i.shape, freq)
        ffreq_im1 = np.full(values_im1.shape, freqm1)
        num = np.log10(values_i) - np.log10(values_im1)
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        return values

    def get_ej_lc_spec_idx(self, freq, ishell=None, ilayer=None, freq1=None, freq2=None):
        freqs = self.get_ej_lc_freqs()
        times = self.get_ej_lc_times()
        if ((freq is None) and (ishell is None) and (ilayer is None)):
            if (not freq1 is None) or (not freq2 is None):
                raise KeyError("For 2d spectral index freq1 and freq2 should be None (for now)")
            arr = []
            for ifreq in range(1, len(freqs)):
                arr.append( self._alpha(freq=freqs[ifreq], freqm1=freqs[ifreq - 1], ishell=ishell, ilayer=ilayer, v_n=None) )
            return np.reshape(np.array(arr), newshape=(len(freqs) - 1, len(times)))
        else:
            idx = find_nearest_index(freqs, freq)
            _freq1 = freqs[idx] if freq1 is None else freq1
            _freq2 = freqs[idx - 1] if freq2 is None else freq2
            return self._alpha(freq=_freq1, freqm1=_freq2, ishell=ishell, ilayer=ilayer, v_n=None)

    def get_ej_lc_temp_idx(self, freq, ishell, ilayer):
        freqs = self.get_ej_lc_freqs()
        t = self.get_ej_lc_times()
        Lnu = self.get_ej_lc(freq=freq, ishell=ishell, ilayer=ilayer, v_n=None)
        val = np.zeros_like(Lnu)
        val[1:-1] = np.log10(Lnu[2:] / Lnu[:-2]) / np.log10(t[2:] / t[:-2])
        val[0] = np.log10(Lnu[1] / Lnu[0]) / np.log10(t[1] / t[0])
        val[-1] = np.log10(Lnu[-1] / Lnu[-2]) / np.log10(t[-1] / t[-2])
        return val

    # ejecta skymaps

    def get_ej_skymap_obj(self):
        self._check_if_loaded_ej_skymap()
        return self.kn_skymap_dfile

    def get_ej_skymap_times(self):
        self._check_if_loaded_ej_skymap()
        # print(self.ej_skymap.keys())
        return np.array(self.kn_skymap_dfile["times"])

    def get_ej_skymap_freqs(self):
        self._check_if_loaded_ej_skymap()
        return np.array(self.kn_skymap_dfile["freqs"])

    def get_ej_skymap_totfluxes(self, freq, shell=None, time=None):
        self._check_if_loaded_ej_skymap()
        if time is None:
            if (shell is None):
                return np.array(self.kn_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
            else:
                return np.array(self.kn_skymap_dfile["totalflux at freq={:.4e} shell={}".format(freq, shell)])
        else:
            if (shell is None):
                arr = np.array(self.kn_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
                if (len(arr)==1):
                    return arr[0]
                assert self.get_ej_skymap_times().min() < time < self.get_ej_skymap_times().max()
                val = arr[find_nearest_index(self.get_ej_skymap_totfluxes(freq=freq, shell=None, time=None), time)]
                return val
            else:
                raise KeyError("Not finished...")

    def get_ej_skymap_cm(self, all_xrs, all_yrs, all_zz):
        dfile = self.get_ej_skymap_obj()
        # _x = np.concatenate(all_xrs)
        # _y = np.concatenate(all_yrs)
        # _z = np.concatenate(all_zz)
        xc_m, yc_m = compute_position_of_the_flux_centroid(all_xrs, all_yrs, all_zz, float(dfile.attrs["d_l"]))
        return (xc_m, yc_m)

    @staticmethod
    def combine_images(xs, ys, datas, verbose=False, hist_or_int="int", shells=False, nx=200, ny=100, extend=2,
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
                if (np.sum(int_i) == 0):
                    continue
                i_shells.append(ish)
                xmin_neg.append(xrs_i[xrs_i < 0].min())
                xmin_pos.append(xrs_i[xrs_i > 0].min())
                xmax.append(xrs_i.max())
                xmin.append(xrs_i.min())
                ymin_neg.append(yrs_i[yrs_i < 0].min())
                ymin_pos.append(yrs_i[yrs_i > 0].min())
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
                    # nx = np.complex(0, nx)
                    # ny = np.complex(0, ny)
                    # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                    #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                    # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                else:
                    # nx = 100
                    # ny = 100
                    # nx = np.complex(0, nx + 1)
                    # ny = np.complex(0, ny + 1)
                    # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                    # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                    #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                    # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                    # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                    # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

                zz += i_zz
                # all_xrs.append(xrs_i)
                # all_yrs.append(yrs_i)
                # all_zz.append(int_i)
                #
                # if hist_or_int == "int":
                #     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #     zz += i_zz
                #     all_xrs.append(xrs_i)
                #     all_yrs.append(yrs_i)
                #     all_zz.append(int_i)
                # else:
                #     nx = 2000
                #     ny = 1000
                #     nx = np.complex(0, nx + 1)
                #     ny = np.complex(0, ny + 1)
                #     edges_x = np.mgrid[x_grid.min():x_grid.max():nx]
                #     edges_y = np.mgrid[y_grid.min():yrs_i.max():ny]
                #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
            xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
            yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])
            #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

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
                i_zz, _ = np.histogramdd(tuple([xs, ys]), bins=tuple([edges_x, edges_y]), weights=datas)
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

    def get_combained_ej_skymaps_adjusted_to_other(self, time, freq, other_pb_instance, nx=100, ny=50):

        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, retrun_edges=True)
        all_x_m1, all_y_m1, all_fluxes_m1 \
            = other_pb_instance.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        _, _, int_zz_m1 = other_pb_instance.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                                           hist_or_int="hist", shells=True, nx=nx,
                                                           ny=ny, edges_x=edges_x, edges_y=edges_y)

        return (int_x, int_y, int_zz, int_zz_m1)

    def get_combined_ej_spectral_map(self, time, freq, nx=100, ny=50, extend=2):
        freqs = self.get_ej_skymap_freqs()
        assert len(freqs) > 2
        idx = find_nearest_index(freqs, freq)
        assert idx != len(freqs)-1
        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, extend=extend, retrun_edges=True)

        all_x_m1, all_y_m1, all_fluxes_m1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        all_x_p1, all_y_p1, all_fluxes_p1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        # ffreq_i = np.full(int_zz.shape, freqs[idx])
        ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        ffreq = np.full(int_zz.shape, freqs[idx])
        # num = np.log10(int_zz) - np.log10(int_zz_m1)
        num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        # denum = np.log10(ffreq) - np.log10(ffreq_im1)
        denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        # int_zz = num / denum
        int_zz = 1.0 * num / denum / 2.0

        return (int_x, int_y, int_zz)

    def get_combined_jet_spectral_map(self, time, freq, nx=100, ny=50, extend=2):

        # TODO NOT TESTED

        freqs = self.get_jet_skymap_freqs()
        assert len(freqs) > 2
        idx = find_nearest_index(freqs, freq)
        assert idx != len(freqs)-1
        all_x, all_y, all_fluxes \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=False, nx=nx,
                                                                     ny=ny, extend=extend, retrun_edges=True)

        all_x_m1, all_y_m1, all_fluxes_m1 \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                              hist_or_int="hist", shells=False, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        all_x_p1, all_y_p1, all_fluxes_p1 \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
                                              hist_or_int="hist", shells=False, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        # ffreq_i = np.full(int_zz.shape, freqs[idx])
        ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        int_zz = 1. * num / denum / 2.

        return (int_x, int_y, int_zz)

    def get_combined_kn_grb_spectral_map(self, time, freq, nx=100, ny=50, extend=2):

        # TODO NOT TESTED

        raise ValueError(" I AM NOT DONE!!! ")

        # freqs = self.get_ej_skymap_freqs()
        # assert len(freqs) > 2
        # idx = find_nearest_index(freqs, freq)
        # assert idx != len(freqs)-1
        # all_x, all_y, all_fluxes \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        # int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
        #                                                            hist_or_int="hist", shells=False, nx=nx,
        #                                                            ny=ny, extend=extend, retrun_edges=True)
        #
        # all_x_m1, all_y_m1, all_fluxes_m1 \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        # _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
        #                                     hist_or_int="hist", shells=False, nx=nx,
        #                                     ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)
        #
        # all_x_p1, all_y_p1, all_fluxes_p1 \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        # _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
        #                                     hist_or_int="hist", shells=False, nx=nx,
        #                                     ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)
        #
        # # ffreq_i = np.full(int_zz.shape, freqs[idx])
        # ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        # ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        # num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        # denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        # int_zz = 1. * num / denum / 2.
        #
        # return (int_x, int_y, int_zz)

    def get_ej_skymap(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False, renormalize=True):

        # nx = 200
        # ny = 100
        # min_mu = 1e-4 # limit for mu that if too small blow up the intensity for some reason...
        # min_mu_frac = 0.85
        # # limit_data = True
        times = self.get_ej_skymap_times()
        freqs = self.get_ej_skymap_freqs()
        dfile = self.get_ej_skymap_obj()
        nshells = int(dfile.attrs["nshells"])
        d_l = float(self.get_ej_skymap_obj().attrs["d_l"])
        if ((not time is None) and (not time in times)):
            raise ValueError(
                "time={} day is not in the list for skypams={} days".format(time / cgs.day, times / cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
        if ((not ishell is None) and (ishell > nshells - 1)):
            raise ValueError("shell={} is beying skymap nshells={}".format(ishell, nshells))
        if ((not time is None) and (not freq is None)):
            # print(dfile.keys())
            ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
            if (not ishell is None):
                r_i = np.array(ddfile["r"][ishell])
                mu_i = np.array(ddfile["mu"][ishell])
                xrs_i = np.array(ddfile["xrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
                yrs_i = np.array(ddfile["yrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
                int_i = np.array(ddfile["intensity"][ishell]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                gam_i = np.array(ddfile["gamma"][ishell])
                B_i = np.array(ddfile["B"][ishell])
                tb_i = np.array(ddfile["tburst"][ishell])
                theta_i = np.array(ddfile["theta"][ishell])
                phi_i = np.array(ddfile["phi"][ishell])

                if remove_mu:
                    print("Removing 'mu' from ejecta skymap")
                    int_i *= np.abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!

                if renormalize:
                    print("Renormalizing ejecta skymap (shell my shell separately)")
                    fnus = self.get_ej_skymap_totfluxes(freq=freq, shell=ishell)
                    fnu = fnus[find_nearest_index(self.get_ej_skymap_times(), time)]
                    all_fluxes_arr = np.array(int_i)
                    delta_x = np.array(xrs_i).max() - np.array(xrs_i).min()
                    delta_y = np.array(yrs_i).max() - np.array(yrs_i).min()
                    dfnu = fnu / (delta_x * delta_y)
                    if verbose:
                        print("Ejecta shell {}".format(ishell))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(xrs_i).min(), np.array(xrs_i).max()))
                        print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(yrs_i).min(), np.array(yrs_i).max()))
                        print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                    int_i = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                    fnus_tot = int_i + fnus

                return (xrs_i, yrs_i, int_i)

                # plt.figure(figsize=(4.6, 3.2))
                # plt.semilogy(mu_i, int_i / int_i.max(), '.')
                # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
                # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                # plt.xlabel(
                #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                # plt.ylabel(r"$I/I_{\rm max}$")
                # plt.tight_layout()
                # plt.show()

                # nx = np.complex(0, nx)
                # ny = np.complex(0, ny)
                # grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx, yrs_i.min():yrs_i.max():ny]
                # i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

                # if int_or_hist == "int":
                #     nx = np.complex(0, nx)
                #     ny = np.complex(0, ny)
                #     grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #                      yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                # else:
                #     nx = 100
                #     ny = 100
                #     nx = np.complex(0, nx + 1)
                #     ny = np.complex(0, ny + 1)
                #     edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                #     edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #     # print(i_zz.shape)
            else:
                i_shells = []
                # total_fluxes = []
                for ish in range(nshells):
                    xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                    if (np.sum(int_i) == 0):
                        continue
                    i_shells.append(ish)
                if remove_mu:
                    print("Removing 'mu' from ejecta skymap")
                all_xrs, all_yrs, all_zz = [], [], []
                for ii, ish in enumerate(i_shells):
                    mu_i = np.array(ddfile["mu"][ish])
                    xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    int_i = np.array(ddfile["intensity"][ish]) * ( d_l ** 2 / cgs.rad2mas ** 2 )  # * dfile.attrs["d_L"] ** 2
                    if remove_mu:
                        # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                        # int_i[idx1] =
                        # print("len(gam[idx1])={}".format(gam_i[idx1]))
                        # idx = int_i > int_i.max() * min_mu_frac
                        # idx = idx1
                        # if verbose:
                        #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                        #         print('\t', mu_i[idx1])
                        #         print('\t', mu_i[idx])
                        #         # exit(1)
                        # if verbose:
                        #     if len(idx) < 2:
                        #         print("No excess = {}".format(ii))
                        #     print("len(gam[idx])={}".format(gam_i[idx]))
                        #     print("Gamma={}".format(gam_i[idx]))
                        #     print("beta={}".format(get_beta(gam_i[idx])))
                        #     print("mu_i={}".format(mu_i[idx]))
                        #     print("r_i={}".format(r_i[idx]))
                        #     print("B={}".format(B_i[idx]))
                        #     print("tb_i={}".format(tb_i[idx]))
                        # print(np.abs(mu_i).min())
                        # print(theta_i[np.abs(mu_i) < 1e-5])
                        # print(phi_i[np.abs(mu_i) < 1e-5])
                        # print(np.abs(mu_i).min())

                        # int_i[mu_i < 1e-3] = 0.
                        int_i *= abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                        # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                        # plt.figure(figsize=(4.6,3.2))
                        # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                        # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                        # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                        # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                        # plt.ylabel(r"$I/I_{\rm max}$")
                        # plt.tight_layout()
                        # plt.show()
                        # int_i[idx] = 0.
                    all_xrs.append(xrs_i)
                    all_yrs.append(yrs_i)
                    all_zz.append(int_i)
                # Problem: changing nlayers changes the image/image, Fnu per pixel; Solution:
                if renormalize:
                    print("Renormalizing ejecta skymap (shell by shell separately)")
                    fnus_tot = np.zeros_like(self.get_ej_skymap_times())
                    for i, i_ish in enumerate(i_shells):
                        fnus = self.get_ej_skymap_totfluxes(freq=freq, shell=i_ish)
                        fnu = fnus[find_nearest_index(self.get_ej_skymap_times(), time)]
                        all_fluxes_arr = np.array(all_zz[i])
                        delta_x = np.array(all_xrs[i]).max() - np.array(all_xrs[i]).min()
                        delta_y = np.array(all_yrs[i]).max() - np.array(all_yrs[i]).min()
                        dfnu = fnu / (delta_x * delta_y)
                        if verbose:
                            print("SHELL {}".format(i_ish))
                            print("\tfnu = {:.2e} ".format(fnu))
                            print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(all_xrs).min(), np.array(all_xrs).max()))
                            print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(all_yrs).min(), np.array(all_yrs).max()))
                            print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                            print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        all_zz[i] = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                        fnus_tot = fnus_tot + fnus

                return (all_xrs, all_yrs, all_zz)

                # assess what is the combined grid extend (for all images)
                # xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
                # ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
                # i_min, i_max = [], []
                # i_shells = []
                # for ish in range(nshells):
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
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
                # xmin = np.array(xmin)
                # ymin_neg = np.array(ymin_neg)
                # ymin_pos = np.array(ymin_pos)
                # ymax = np.array(ymax)
                # ymin = np.array(ymin)
                # i_min = np.array(i_min)
                # i_max = np.array(i_max)
                #
                # if verbose:
                #     for i in range(len(i_shells)):
                #         print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i], xmax[i]))
                #         print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],  ymax[i]))
                #     print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
                #     print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
                #     print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
                #
                # edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
                # edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
                #
                #
                # x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
                #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
                #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                #
                # # edges_x = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # # edges_y = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                # # plt.plot(edges_x, edges_y, marker='.', ls='none')
                # # plt.show()
                #
                # if verbose:
                #     print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                #                                                                     x_grid[x_grid > 0].min(),
                #                                                                     x_grid.max()))
                #     print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                #                                                                     y_grid[y_grid > 0].min(),
                #                                                                     y_grid.max()))
                # xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
                # if int_or_hist == "int":
                #     zz = np.zeros_like((xx_grid))
                # else:
                #     zz = np.zeros((len(edges_x)-1,len(edges_y)-1))
                # # interpolate onto the grid that covers all images
                # all_xrs, all_yrs, all_zz = [], [], []
                # for ii, ish in enumerate(i_shells):  # range(len(i_shells))
                #     if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
                #     r_i = np.array(ddfile["r"][ish])
                #     mu_i = np.array(ddfile["mu"][ish])
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2) # * dfile.attrs["d_L"] ** 2
                #     gam_i = np.array(ddfile["gamma"][ish])
                #     B_i = np.array(ddfile["B"][ish])
                #     tb_i = np.array(ddfile["tburst"][ish])
                #     theta_i = np.array(ddfile["theta"][ish])
                #     theta_ij = np.array(ddfile["theta_j"][ish])
                #     theta0_i = np.array(ddfile["theta0"][ish])
                #     phi_i = np.array(ddfile["phi"][ish])
                #
                #
                #     # plt.plot(mu_i, theta0_i, marker='.', ls='none')
                #     # plt.show()
                #     # int_i[mu_i < 5e-2] = 0.
                #     # int_i *= abs(mu_i)
                #
                #     if remove_mu:
                #         # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                #         # int_i[idx1] =
                #         # print("len(gam[idx1])={}".format(gam_i[idx1]))
                #         # idx = int_i > int_i.max() * min_mu_frac
                #         # idx = idx1
                #         # if verbose:
                #         #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                #         #         print('\t', mu_i[idx1])
                #         #         print('\t', mu_i[idx])
                #         #         # exit(1)
                #         # if verbose:
                #         #     if len(idx) < 2:
                #         #         print("No excess = {}".format(ii))
                #         #     print("len(gam[idx])={}".format(gam_i[idx]))
                #         #     print("Gamma={}".format(gam_i[idx]))
                #         #     print("beta={}".format(get_beta(gam_i[idx])))
                #         #     print("mu_i={}".format(mu_i[idx]))
                #         #     print("r_i={}".format(r_i[idx]))
                #         #     print("B={}".format(B_i[idx]))
                #         #     print("tb_i={}".format(tb_i[idx]))
                #         # print(np.abs(mu_i).min())
                #         # print(theta_i[np.abs(mu_i) < 1e-5])
                #         # print(phi_i[np.abs(mu_i) < 1e-5])
                #         # print(np.abs(mu_i).min())
                #
                #         # int_i[mu_i < 1e-3] = 0.
                #         int_i *= abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                #         # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                #         # plt.figure(figsize=(4.6,3.2))
                #         # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                #         # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                #         # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                #         # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                #         # plt.ylabel(r"$I/I_{\rm max}$")
                #         # plt.tight_layout()
                #         # plt.show()
                #         # int_i[idx] = 0.
                #     # int_i /= abs(mu_i)[abs(mu_i)>0.01]
                #     # else:
                #     #     int_i2 = int_i
                #
                #     if int_or_hist == "int":
                #         # nx = np.complex(0, nx)
                #         # ny = np.complex(0, ny)
                #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
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
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #
                #     # i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #     zz += i_zz
                #     all_xrs.append(xrs_i)
                #     all_yrs.append(yrs_i)
                #     all_zz.append(int_i)
                # if verbose:
                #     print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
                #           .format(zz.min(), zz.max(), np.sum(zz), float(np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
                #
                # all_xrs = np.concatenate(all_xrs)
                # all_yrs = np.concatenate(all_yrs)
                # all_zz = np.concatenate(all_zz)
                # if int_or_hist == "int":
                #
                #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
                # else:
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, zz, all_xrs, all_yrs, all_zz)
        else:
            raise NotImplementedError("Not implemented")



    def get_ej_skymap_spec_idx(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False):

        freqs = self.get_ej_skymap_freqs()
        assert len(freqs) > 1
        idx = find_nearest_index(freqs, freq)
        dfile = self.get_ej_skymap_obj()
        nshells = int(dfile.attrs["nshells"])

        all_xrs, all_yrs, all_zz = \
            self.get_ej_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

        all_xrs_m1, all_yrs_m1, all_zz_m1 = \
            self.get_ej_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

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
    @staticmethod
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

        image = interpolate.griddata(np.array([all_x, all_y]).T * fac, all_fluxes,
                                     (grid_x * fac, grid_y * fac), method='linear', fill_value=0)
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
    @staticmethod
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

    @staticmethod
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