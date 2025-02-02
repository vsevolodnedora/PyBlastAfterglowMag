'''

    This script reads output files from postprocessing Kenta Kiuchi code that
    Input files are: * ejecta.h5 *
    These files two-dimensional histograms of the ejecta velocity distribution.
    The script removes one half of the sphere from the angular distribution,
    and parts of it where mass is 0 or velocity is negative, e.g., it clears the
    data so PyBlastAfterglow does not need to deal with unphysical initial data.

    Usage
    python3 ./id_maker_from_thc_outflow.py -i path/to/hist_or_corr.h5 -o path/to/output.h5 -m corr_or_hist -l 30 --factor 2.

'''
import matplotlib.pyplot as plt
import numpy as np
import h5py
import pandas as pd
from scipy import interpolate
import copy
import glob
import re
import os
import sys
import argparse
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from .id_tools import (reinterpolate_hist, reinterpolate_hist2, compute_ek_corr)
from .utils import (cgs, find_nearest_index,
                    GammaFromMom,GammaFromBeta,BetaFromMom,BetaFromGamma,MomFromBeta,MomFromGamma)

class ProcessRawFiles:
    """
        Processes files from Kenta Kiuchi BNS merger simulations
    """
    def __init__(self, files : list[str], verbose:bool,mode:str="mass"):

        self.files = files

        if (mode == "mass"):
            self.key_v_asymptotic = "v_asymptotic"
            self.key_theta = "theta"
            self.key_Mejecta = "Mejecta"
            self.expected_keys = ['Mejecta', 'T_ejecta', 'Ye_ejecta', 'entropy_ejecta',
                                  'internal_energy_ejecta', 'pressure_ejecta', 'rho_ejecta', 'time']
            self.key_to_key = {
                "T_ejecta":"temp",
                "Ye_ejecta":"ye",
                "entropy_ejecta":"entr",
                "internal_energy_ejecta":"eps",
                "pressure_ejecta":"press",
                "rho_ejecta":"rho",
                "Mejecta":"mass"
            }

        elif mode == "mdot":
            self.key_v_asymptotic = "R_ext"
            self.key_theta = "theta"
            self.key_Mejecta = "Mdot"
            self.expected_keys = ['Mdot', 'Mdot2', 'Mdot3', 'Mdot4', 'Mdot5', 'time']
            self.key_to_key = {
                "Mdot":"mdot_total", # Mdot : Total
                "Mdot2":"mdot_slow", # Mdot2: v^r 0.6c
                "Mdot3":"mdot_slow", # Mdot3: \beta \Gamma < 0.1
                "Mdot4":"mdot_mid", # Mdot4: 0.1 < \beta \Gamma < 1
                "Mdot5":"mdot_fast", # Mdot5: 1 < \beta Gamma
            }

        else:
            raise KeyError(f"Only 'mass' or 'mdot' modes are supported. Given = {mode}")

        self.verb=verbose
        if len(files) == 0:
            raise ValueError(f"Files list is empty")

    def _get_data(self) -> tuple[np.ndarray,np.ndarray,tuple[np.ndarray,list[dict]]]:
        pars_list = []
        texts = []
        v_inf = np.zeros(0,)
        thetas = np.zeros(0,)

        _total_mass = 0.

        for i, fl in enumerate(self.files):
            if self.verb: print("\t Processing File {}".format(fl))
            dfile = h5py.File(fl, "r")
            if self.verb:
                print("\t Keys in the dfile: {}".format(dfile.keys()))

            if (len(dfile.keys()) == 0):
                print(f"\t\t SKIPPING FILE (no keys in it) {fl}")
                continue

            if self.verb:
                print("\t {}              = {} [{}, {}]".format(self.key_theta,
                                                                np.array(dfile[self.key_theta]).shape,
                                                                np.array(dfile[self.key_theta])[0],
                                                                np.array(dfile[self.key_theta])[-1]))
                print("\t {}       = {} [{}, {}]".format(self.key_v_asymptotic,
                                                         np.array(dfile[self.key_v_asymptotic]).shape,
                                                         np.array(dfile[self.key_v_asymptotic])[0],
                                                         np.array(dfile[self.key_v_asymptotic])[-1]))
                print("\t dfile['data1'].keys= {}".format(dfile["data1"].keys()))

            # sort the keys in the file (data groups for different extraction times)
            sort_by = lambda k: int(re.findall(r'\d+', k)[0]) if k.__contains__("data") else -1
            keys = sorted(dfile.keys(), key=sort_by)
            tkeys = [key for key in keys if key.__contains__("data")]

            # check what data groups are not empty
            times = []
            failed_to_get_time = []
            for key in keys:
                if key.__contains__("data"):
                    try:
                        times.append(np.array(dfile[key]["time"], dtype=np.float64)[0])
                    except KeyError:
                        failed_to_get_time.append(key)
                        if self.verb:
                            print(f"Count not extract time from key={key} dfile[key]={dfile[key]}")
            if self.verb:
                print(f"Failed to extract time for {len(failed_to_get_time)}/{len(keys)}")

            # get unique times from the dataframe (some might be too close)
            unique_times = list(set(np.around(times, 0)))
            if len(unique_times) == 1:
                idxs = [0]
            else:
                idxs = [find_nearest_index(times, u_time) for u_time in unique_times]

            # process data from each group (each timestep)
            for idx in idxs:

                # load 2D histogram axis data (velocity and angle)
                v_inf = np.array(dfile[self.key_v_asymptotic], dtype=np.float64)
                if not (np.all(v_inf[:-1] <=v_inf[1:])):
                    sort_index = np.argsort(v_inf)
                    v_inf = np.sort(v_inf)
                else:
                    sort_index = np.argsort(v_inf)
                thetas = np.array(dfile[self.key_theta], dtype=np.float64)
                # if (self.key_v_asymptotic == "R_ext"): v_inf *= 0.4816 # ul = 0.4816 # km code -> km


                # load the histogram weights (mass)
                mass = np.array(dfile[tkeys[idx]][self.key_Mejecta], dtype=np.float64)

                if self.verb:
                    print("Processing: time={} key={} {}={}"
                          .format(times[idx], tkeys[idx], self.key_Mejecta, np.sum(mass)))

                res = {}
                for key, new_key in self.key_to_key.items():
                    if key in list(dfile[tkeys[idx]].keys()):
                        arr = np.array(dfile[tkeys[idx]][key], dtype=np.float64)
                        # apply units
                        if self.verb: print(f"\tFound '{key}' sahpe={arr.shape} min={arr.min()} max={arr.max()} sum={arr.sum()}")
                        if key == "T_ejecta": arr *= 11604525006.17 # MeV -> Kelvin
                        if key == "rho_ejecta": arr *= 5.807e18
                        # if key in ["Mdot","Mdot2","Mdot3","Mdot4","Mdot5"]: arr *= (0.326 / 1.607e-6) # (um / ut) um = 0.326 ut = 1.607e-6 -> Msun/s

                        # apply sorting
                        res[new_key] = arr[:,sort_index]

                # ek = compute_ek_corr(v_inf, res[mass]).T

                pars_list.append(copy.deepcopy(res))
                texts.append(times[idx])

        if self.verb:
            print(" N of Total iterations : {}".format(len(pars_list)))
        if self.verb:
            print("Total times            : {}".format(np.array(texts,dtype=float)))

        # sort the list with respect to time
        sorted_pars_list = []
        texts = np.array(texts)
        sorted_texts = np.array(np.sort(texts))
        for _time in sorted_texts:
            for i in range(len(pars_list)):
                if (_time == np.array(texts)[i]):
                    sorted_pars_list.append(pars_list[i])
        assert len(pars_list) == len(sorted_pars_list)

        summed_vals = {}
        for key in sorted_pars_list[0].keys():
            summed_vals[key] = 0.
        for i, _t in enumerate(sorted_texts):
            for key in sorted_pars_list[i].keys():
                summed_vals[key] += np.sum(sorted_pars_list[i][key])

        if self.verb:
            print(" -------------------------------------- ")
            print(f"Total interations: {len(sorted_pars_list)} [{sorted_texts[0]:.1f}-{sorted_texts[-1]:.1f}]")
            for key in sorted_pars_list[0].keys():
                print(f"key={key} sum()={summed_vals[key]}")
            print(" -------------------------------------- ")

        return (v_inf, thetas, (sorted_texts, sorted_pars_list))

    def process_save(self, outfnmae : str):

        vinf, thetas, (times, datas) = self._get_data()

        # save data as a single file
        try:
            with h5py.File(outfnmae,"w") as f:
                f.create_dataset(self.key_v_asymptotic,data=vinf)
                f.create_dataset(self.key_theta,data=thetas)
                f.create_dataset("text",data=times)
                for t, d in zip(times, datas):
                    group = f.create_group(name="time={:.4f}".format(t))
                    for key, arr in d.items():
                        group.create_dataset(name=key, data=arr)
        except BlockingIOError:
            raise BlockingIOError(f"Cannot open file to write error {outfnmae}")

class EjectaData:
    def __init__(self, fpath : str, verbose : bool, mode : str = "mass"):
        self.verb = verbose
        self.fpath = fpath
        self.dfile = None
        self.texts = np.zeros(0,)
        self.vinf = np.zeros(0,)
        self.theta = np.zeros(0,)
        self.v_ns = ["temp","ye","entr","eps","press","rho","mass"]

        self._load(mode)

    def _load(self, mode : str):
        if (not os.path.isfile(self.fpath)):
            raise FileNotFoundError(f"Collated ejecta file is not found. {self.fpath}")
        if self.dfile is None:
            if (mode == "mass"):
                self.dfile = h5py.File(self.fpath, "r")
                self.texts = np.array(self.dfile["text"])
                self.vinf = np.array(self.dfile["v_asymptotic"])
                self.theta = np.array(self.dfile["theta"])
            elif (mode == "mdot"):
                self.dfile = h5py.File(self.fpath, "r")
                self.texts = np.array(self.dfile["text"])
                self.rext = np.array(self.dfile["R_ext"])
                self.theta = np.array(self.dfile["theta"])
            else:
                raise KeyError(f"mode={mode} is not supported")
        if (self.theta.max() > 3. * np.pi / 4.):
            if self.verb:
                print(f"WARNING: theta grid extends to {self.theta.max()} for  {self.fpath} ")
                print("Performing averaging between two hemispheres")
            self.idx = len(self.theta)//2
        else:
            self.idx = -1
        if self.idx > 0 :
            self.theta = self.theta[:self.idx+1]

        #     if self.verb:
        #         print(f"WARNING: theta grid extends to {self.theta.max()} for  {self.fpath} ")
        #         print("Performing averaging between two hemispheres")
        #     fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(4.6,3.2),layout='constrained',sharey='all')
        #     cmap = plt.get_cmap('viridis')
        #     mass = self.dfile['mass']
        #     norm = LogNorm(vmin=mass.max()*1e-7, vmax=ek.max())
        #     im = ax.pcolormesh(ctheta*180./np.pi,mom, ek, cmap=cmap, norm=norm)
        #     # ax0.set_yscale('log')
        #     # fig.colorbar(im, ax=ax)
        #     ax.set_yscale('log')

            # raise ValueError("theta grid is too extensive")

    def get_theta(self):
        return self.theta

    def get_vinf(self):
        return self.vinf

    def get_rext(self):
        return self.rext

    def get(self, v_n : str, text : float):
        if (text not in self.texts):
            if (text < self.texts.min()):
                raise ValueError(f"text={text} < texts.min()={self.texts.min()}")
            if (text > self.texts.max()):
                raise ValueError(f"text={text} > texts.max()={self.texts.max()}")
            text = self.texts[find_nearest_index(self.texts, text)]
        ddfile = self.dfile["time={:.4f}".format(text)]
        if not v_n in ddfile.keys():
            raise KeyError(f"key={v_n} is not in the dfile.keys():\n {ddfile.keys()}")
        ddfile = self.dfile["time={:.4f}".format(text)]
        arr = np.array(ddfile[v_n])[:self.idx,:]
        if self.idx > 0 :
            return arr[:self.idx,:]
        else:
            return arr
        # check if the data is on one hemisphere or not
        # if (self.idx > -1):
        #     if self.verb:
        #         print(f"WARNING: theta grid extends to {self.theta.max()} for  {self.fpath} ")
        #         print("Performing averaging between two hemispheres")
        #     assert len(arr) == len(self.theta)
        #     assert len(arr[0]) == len(self.vinf)
        #
        #     fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(4.6,3.2),layout='constrained',sharey='all')
        #     cmap = plt.get_cmap('viridis')
        #     norm = LogNorm(vmin=arr.max()*1e-7, vmax=arr.max())
        #     im = ax.pcolormesh(self.theta,self.vinf, arr.T, cmap=cmap, norm=norm)
        #     # ax0.set_yscale('log')
        #     # fig.colorbar(im, ax=ax)
        #     # ax.set_yscale('log')
        #     plt.show()
        # return arr

    def total_mass(self) -> np.ndarray:
        mass = np.array([np.sum(self.get(v_n="mass",text=text)) for text in self.texts])
        return mass

    def get_vinf_mask(self, crit : str) -> np.ndarray:
        if (crit is None):
            return np.ones_like(self.vinf).astype(bool)
        elif (crit == "fast"):
            mask = MomFromBeta(self.vinf) > 1# self.vinf * GammaFromBeta(self.vinf) > 1
        elif (crit == "mid"):
            mask = (MomFromBeta(self.vinf)<=1 & MomFromBeta(self.vinf)>0.1) # (self.vinf * GammaFromBeta(self.vinf) <= 1) & \
                      #(self.vinf * GammaFromBeta(self.vinf) > 0.1)
        elif (crit == "slow"):
            mask = MomFromBeta(self.vinf)<=0.1#(self.vinf * GammaFromBeta(self.vinf) <= 0.1)
        else:
            raise KeyError()
        return mask

    def total_mass_vs_text(self,crit:str="fast") -> np.ndarray:
        mask = self.get_vinf_mask(crit=crit)
        mass = np.array([np.sum(self.get(v_n="mass",text=text)[:, mask]) for text in self.texts])
        return mass


    def getText(self) -> np.ndarray:
        return np.array(self.dfile["text"])

class Data():

    def __init__(self,
                 fpath_rhomax:str,
                 fpath_mdot:str):
        # file names
        self.fpath_rhomax = fpath_rhomax
        self.fpath_mdot = fpath_mdot

        # dataframes
        self.df_rho = pd.DataFrame()
        self.df_mdot = pd.DataFrame()

        self._load_rhomax_df()
        self._load_mdot_df()

    def _load_rhomax_df(self) -> None:
        if os.path.isfile(self.fpath_rhomax):
            # units
            rho_nuc = 2.7e14
            ut = 1.607e-6
            urho = 5.807e18
            # load data
            t_rho_max, rho_max = np.loadtxt(self.fpath_rhomax, unpack=True, usecols=(0,1))
            # make dataframe
            data = {
                "time" : t_rho_max * ut * 1e3,
                "rho_max": rho_max * urho / rho_nuc
            }
            self.df_rho = pd.DataFrame.from_dict(data=data)
        else:
            print("File not found: {}".format(self.fpath_rhomax))

    def _read_mdot_table2(self) -> np.ndarray:
        table = []
        try:
            table = np.loadtxt(self.fpath_mdot, unpack=True).T
        except:
            print("failed to read table. Trying line by line")
            with open(self.fpath_mdot, "r") as file:
                lines = file.readlines()
                for il, line in enumerate(lines):
                    vals = line.split()
                    nums = []
                    for val in vals:
                        try:
                            nums.append(np.float64(val))
                        except:
                            nums.append(np.nan)
                    table.append(nums)
                table = np.reshape(np.array(table),newshape=(il, len(nums)))
        return (table)
    def get_rext(self) -> np.ndarray:
        ul = 0.4816 # km
        rext = np.array([50, 100, 300, 1000, 3000, 4000, 10000, 1229.3])
        return np.array(rext*ul, dtype=np.int64)
    def _load_mdot_df(self) -> None:
        if os.path.isfile(self.fpath_mdot):
            mdot_table = self._read_mdot_table2()

            # convert the units from what Kenta was using to CGS + Msun
            ul = 0.4816 # km
            ut = 1.607e-6
            um = 0.326
            umdot = um / ut
            mdot_table[:,0] *= (ut * 1e3) # -> s -> ms
            mdot_table[:,1:41] *= umdot # to Msun/s
            rext = np.array([50, 100, 300, 1000, 3000, 4000, 10000, 1229.3]) # 8 radii  in code unit.
            rext *= ul # code -> km

            # create dict for data
            df = {
                "time":mdot_table[:,0],
            }
            i = 1
            for r in rext:
                df[f"mdot_r{int(r)}"] = mdot_table[:, i + len(rext) * 0] # [:,1:9] 2-9 Mdot_ext (8 extraction radii in code unit, all contribution)
                df[f"mdot_vr_r{int(r)}"] = mdot_table[:, i + len(rext) * 1] # [:,9:17] 10-17 Mdot_ext2 (8 extraction radii in code unit, vr > 0.6c)
                df[f"mdot_slow_r{int(r)}"] = mdot_table[:, i + len(rext) * 2] # 18-25 Mdot_ext3 (8 extraction radii in code unit, \gamma_beta <= 0.1)
                df[f"mdot_mid_r{int(r)}"] = mdot_table[:, i + len(rext) * 3] # 26-33 Mdot_ext4 (8 extraction radii in code unit, 0.1 <= \gamma_beta < 1))
                df[f"mdot_fast_r{int(r)}"] = mdot_table[:, i + len(rext) * 4] # 34-41 Mdot_ext5 (8 extraction radii in code unit, \gamma_beta >= 1)
                i+=1
            i+=len(rext)*4
            # velocity
            df["vave_slow"] = mdot_table[:, i] # 42 mass_averaged_velocity (\gamma_beta <= 0.1)
            df["vave_mid"]  = mdot_table[:, i+1] # 43 mass_averaged_velocity (0.1 <= \gamma_beta < 1)
            df["vave_fast"] = mdot_table[:, i+2] # 44 mass_averaged_velocity (\gamma_beta >= 1)

            # create dataframe
            self.df_mdot = pd.DataFrame.from_dict(data=df,orient="columns")
            self.df_mdot.set_index(keys="time")
        else:
            print(f"File not found {self.fpath_mdot}")
    def ret_time(self, r0 : float, vave_key : str) -> float:
        """
            r0 : rm
            speed : c
        """
        vinf = self.df_mdot[vave_key].values
        if np.sum(vinf)==0:
            raise ValueError("All vinf=0")
        t_ret_0 = (r0/1000./vinf/cgs.c) # s
        return self.df_mdot["time"].values - t_ret_0

class EjStruct(EjectaData):
    def __init__(self, fpath : str, verbose : bool):
        super().__init__(fpath, verbose)

    def get_r0(self, vinf, theta, mass, rho, t0 : float, method : str) -> np.ndarray:
        r = np.zeros_like(mass)
        if (method == "from_rho"):
            if (t0 < 0):
                raise ValueError(f" Not set t0={t0} must be > 0 in seconds for: r_base = t0 * cgs.c * vinf[0]")
            for ith in range(len(theta)):
                r_base = t0 * cgs.c * vinf[0]
                r[0,ith] = r_base
                for ir in range(1, len(vinf)):
                    if (mass[ir,ith]==0. or rho[ir,ith]==0.):
                        print(f"Error ir={ir} ith={ith} mass={mass[ir,ith]} rho={rho[ir,ith]}. Setting mass to 0.")
                        mass[ir,ith]=0.
                        continue
                    r_i = (3./4.) * (1./np.pi) * len(theta) * mass[ir,ith] / rho[ir,ith] + r[ir - 1,ith] ** 3
                    r_i = r_i**(1./3.)
                    r[ir,ith] = r_i
                    if ((r_i <= r[ir-1,ith]) or (~np.isfinite(r[ir,ith]))):
                        raise ValueError()
        elif (method == "from_beta"):
            t = t0
            for ith in range(len(theta)):
                for ir in range(len(vinf)):
                    r[ir,ith] = vinf[ir] * cgs.c * t # TODO THis is theta independent!
        else:
            raise KeyError(f"method={method} is not recognized")
        return r

    def get_2D_id(self, text : float, method_r0 : str, t0 = None,
                  force_spherical:bool=False,
                  new_theta_len:None or int = None, new_vinf_len:None or int = None) -> dict:
        res = {}
        masses = self.get(v_n="mass", text=text)[:, :] # [itheta ivinf]



        # fill with average values for each angle if angular structure is not needed
        if (force_spherical):
            sums = np.zeros(len(masses[0,:]))
            for ivinf in range(len(masses[0,:])):
                masses[:,ivinf] = np.sum(masses[:,ivinf])/len(masses[:,ivinf])
                sums[ivinf] = np.sum(masses[:,ivinf])

            # print(repr(self.vinf))
            # print(repr(sums))

        v_inf, thetas, masses = reinterpolate_hist2(self.vinf, self.theta, masses,
                                                    new_theta_len=new_theta_len,
                                                    new_vinf_len=new_vinf_len,
                                                    mass_conserving=True)
        masses = masses.T # -> [i_vinf, i_theta]

        thetas = 0.5 * (thetas[1:] + thetas[:-1])
        v_inf  = 0.5 * (v_inf[1:] + v_inf[:-1])
        mask = ((v_inf > 0) & (v_inf < 1) & (np.sum(masses, axis=1) > 0))
        for v_n in self.v_ns:
            if (v_n != "mass"):
                data = self.get(v_n=v_n, text=text)[:, :]

                # fill with average values for each angle if angular structure is not needed
                if (force_spherical):
                    for ivinf in range(len(data[0,:])):
                        data[:,ivinf] = np.sum(data[:,ivinf])/len(data[:,ivinf])

                _, _, data = reinterpolate_hist2(self.vinf, self.theta, data,
                                                            new_theta_len=new_theta_len,
                                                            new_vinf_len=new_vinf_len,
                                                            mass_conserving=True)
                res[v_n] = data[:, mask].T

        res["mass"] = masses[mask,:] * cgs.solar_m
        res["r"] = self.get_r0(vinf=v_inf[mask], theta=thetas, mass=masses[mask,:], rho=res["rho"],
                               t0=t0, method=method_r0)
        res["ek"] = np.column_stack([masses[mask, i] * cgs.solar_m * v_inf[mask]**2 * cgs.c * cgs.c for i in range(len(thetas))])
        mom = MomFromBeta(v_inf[mask])
        thetas,mom = np.meshgrid(thetas,mom)
        if not (mom.shape == res["ek"].shape):
            raise ValueError("{} != {}".format(mom.shape,res["ek"].shape))

        def plot_final_profile(ctheta, mom, data2d, cmap='jet',vmin=None,vmax=None):
            from matplotlib.colors import LogNorm
            if vmin is None: vmin = data2d[(data2d > 0) & (np.isfinite(data2d))].min()
            if vmax is None: vmax = data2d[(data2d > 0) & (np.isfinite(data2d))].max()
            norm = LogNorm(vmin, vmax)
            fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(4.6,3.2))
            ctheta *= (180. / np.pi)
            im = ax.pcolor(mom, ctheta, data2d, cmap=cmap, norm=norm, shading='auto')
            # ax0.axhline(y=1, linestyle='--', color='gray')

            # adjust the bottom subplot
            ax.set_ylim(0, 90)
            ax.set_xlim(1e-2, 4)
            ax.set_xscale('log')
            ax.set_yscale('linear')
            ax.set_xlabel(r"$\Gamma\beta$", fontsize=12)
            ax.set_ylabel(r"Polar angle", fontsize=12)
            plt.show()

        # plot_final_profile(ctheta=thetas,mom=mom,data2d=res["ek"])
        # plot_final_profile(ctheta=thetas,mom=mom,data2d=res["mass"])

        res["mom"] = mom
        res["theta"] = thetas
        res["ctheta"] = thetas



        # for v_n in res.keys():
        #     res[v_n] = np.copy(res[v_n].T)
        # print(res["mom"][:,0])
        # print(res["mom"][0,:])
        # print(res["ctheta"][0,:])
        return res

    @staticmethod
    def plot_init_profile(mom : np.ndarray, ctheta : np.ndarray, mass : np.ndarray,
                          xmin=0,xmax=90,ymin=1e-2,ymax=6,vmin=1e-12,vmax=1e-6,
                          norm_mode="log", cmap = plt.get_cmap('RdYlBu_r'),
                          xscale="linear",yscale="linear",
                          subplot_mode="sum",
                          title=None, figpath=None):

        fontsize=12
        mom = mom.copy()
        ctheta = ctheta.copy()
        mass = mass.copy()

        # moms = data["mom"]
        # ctheta = data["ctheta"] * 180 / cgs.pi
        # eks = data["ek"]

        ctheta *= 180 / cgs.pi

        fig = plt.figure(figsize=(4.5 + 1, 3.6 + 3))
        # fig.suptitle(r"BLh* $(1.259+1.482)M_{\odot}$")

        ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.59 - 0.12])
        ax1 = fig.add_axes([0.16, 0.61, 0.81 - 0.15, 0.91 - 0.61])
        cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.59 - 0.12])

        #                   x1    y1    delta x       delta y
        # ax1 = fig.add_axes([0.16, 0.45, 0.81 - 0.15, 0.59 - 0.12])
        # ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.91 - 0.61])
        # cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.91 - 0.61])

        # top panel
        # for t in tasks:
        #     hist_vinf, hist_mass = t["data"].load_vinf_hist()
        #     hist_mom = hist_vinf * get_Gamma(hist_vinf)
        #     hist_eks = np.cumsum(np.array(0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m)[::-1])[::-1]

        ax1.plot(ctheta, np.sum(mass,axis=0), color='black', ls='-', drawstyle='steps')

        if (not title is None): ax1.set_title(title)

        ###  mooley
        # mool_mom = np.linspace(0.5 * get_Gamma(0.5), 0.8 * get_Gamma(0.8), 20)
        # mool_ek = 5e50 * (mool_mom / 0.4) ** -5
        # _l, = ax1.plot(mool_mom, mool_ek, color="gray", ls="-")
        # ax1.text(0.30, 0.90, "Mooley+17", color='black', transform=ax1.transAxes, fontsize=12)

        # tasks = [1, 1, 1]
        # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="-", label=r"$q=1.00$")
        # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="--", label=r"$q=1.43$")
        # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls=":", label=r"$q=1.82$")
        # han, lab = ax1.get_legend_handles_labels()
        # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
        #                          **{"fancybox": False, "loc": 'lower left',
        #                            # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #                            "shadow": "False", "ncol": 1, "fontsize": 11,
        #                            "framealpha": 0., "borderaxespad": 0., "frameon": False}))
        #
        # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
        #                           **{"fancybox": False, "loc": 'upper right',
        #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        # "shadow": "False", "ncol": 1, "fontsize": 11,
        # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
        #
        # ax1.add_artist(ax1.legend(han[len(han) - len(tasks):], lab[len(lab) - len(tasks):],
        #                           **{"fancybox": False, "loc": 'center right',
        #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        # "shadow": "False", "ncol": 1, "fontsize": 11,
        # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
        #
        # fit *= np.sum(mdens) / np.sum(fit)
        # fit_x, fit_y = _fit(mom, ekdens)
        # ax1.plot(fit_x, fit_y, 'k--', label=r'$\sin^2(\theta)$')
        # ax1.legend(**{"fancybox": False, "loc": 'upper right',
        #                # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #                "shadow": "False", "ncol": 1, "fontsize": 11,
        #                "framealpha": 0., "borderaxespad": 0., "frameon": False})

        if norm_mode=="log":
            ax1.set_yscale("log")
        else:
            ax1.set_yscale("linear")
        # ax1.set_ylim(ymin, ymax)
        # ax1.yaxis.tick_right()
        # ax1.yaxis.tick_left()
        ax1.set_ylabel(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
        ax1.get_yaxis().set_label_coords(-0.15, 0.5)

        ax1.set_xlim(xmin, xmax)
        ax1.xaxis.set_ticklabels([])

        ax1.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12,
                        direction='in',
                        bottom=True, top=True, left=True, right=True)
        ax1.minorticks_on()

        ax11 = ax1.twinx()

        mask = mom > 1

        if len(mom[mask])>0:
            # axx = np.sum(betas[mask, np.newaxis] * eks[mask, :],axis=0)
            # ayy = np.sum(eks[mask, :],axis=0)
            # ax11.plot(ctheta, axx/ayy, color="gray", ls="-")
            # if (subplot_mode=="sum"): _yarr = np.sum(ctheta*eks[mask, :],axis=0)
            # elif (subplot_mode=="ave"): _yarr = axx/ayy
            # else:raise KeyError("subplot_mode is not recognized")
            ax11.plot(ctheta, np.sum(mass[mask, :],axis=0), color="gray", ls="-")
        if norm_mode=="log":
            ax11.set_yscale("log")
        else:
            ax11.set_yscale("linear")
        ax11.set_ylabel(r"$M_{\rm ej}(\Gamma\beta>1)$ [M$_{\odot}$]", fontsize=fontsize, color="gray")
        ax11.minorticks_on()
        ax11.tick_params(axis='both', which='both', labelleft=False,
                         labelright=True, tick1On=False, tick2On=True,
                         labelsize=12,
                         direction='in',
                         bottom=True, top=True, left=True, right=True)


        # bottom panel
        # mask = vinf > 0.6
        # eks2 = np.zeros_like(eks)
        # for i in range(len(vinf)):
        #     if vinf[i] > 0.6:
        #         eks2[:, i] = eks[:, i]

        # import h5py
        # dset = h5py.File("/home/vsevolod/Desktop/Hajela/BLh_q100.h5", 'w')
        # dset.create_dataset(name="beta", data=vinf)
        # dset.create_dataset(name="theta", data=theta * 180 / np.pi)
        # dset.create_dataset(name="Ek(>Gamma_beta)", data=eks)
        # dset.close()

        if (norm_mode=="log"):
            norm = LogNorm(mass[(mass > 0) & (np.isfinite(mass))].min(), mass[(mass > 0) & (np.isfinite(mass))].max())
        elif (norm_mode=="linear"):
            norm = Normalize(mass[(mass > 0) & (np.isfinite(mass))].min(), mass[(mass > 0) & (np.isfinite(mass))].max())
        elif (norm_mode=="levels"):
            levels = MaxNLocator(nbins=15).tick_values(mass.min(), mass.max())
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        else:
            raise KeyError(" norm_mode is not recognized ")


        im = ax0.pcolor(ctheta, mom, mass, cmap=cmap, norm=norm, shading='auto')
        # cbar = fig.colorbar(im, ax=ax, extend='both')
        # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

        ax0.axhline(y=1, linestyle='--', color='gray')

        ax0.set_ylim(ymin, ymax)
        ax0.set_xlim(xmin, xmax)
        ax0.set_xscale(xscale)
        ax0.set_yscale(yscale)

        ax0.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
        ax0.set_xlabel(r"Polar angle", fontsize=fontsize)
        ax0.get_yaxis().set_label_coords(-0.15, 0.5)
        # ax0.text(0.75, 0.88, lbl, color='white', transform=ax0.transAxes, fontsize=fontsize)
        ax0.minorticks_on()
        ax0.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12,
                        direction='in',
                        bottom=True, top=True, left=True, right=True)

        # ax0.text(0.05, 0.05, models.print_fancy_label(name), color='white', transform=ax0.transAxes)

        # ekdens = np.sum(eks, axis=0)
        # ax1.step(0.5 * (theta[1:] + theta[:-1]), mdens, where='mid', color='red')
        # hist_vinf, hist_mass = o_data.load_vinf_hist()
        # hist_mom = hist_vinf * get_Gamma(hist_vinf)
        # hist_eks = 0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m
        # ax1.step(mom, ekdens, where='mid', color='red')
        # ax1.step(hist_mom, hist_eks, where='mid', color='black')

        cbar = plt.colorbar(im, cax=cax, norm=norm)
        cbar.set_label(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
        cbar.ax.minorticks_off()
        # plt.savefig(PAPERPATH + "kinetic_energy_struct_models.pdf")
        # plt.savefig(FIGPATH + "kinetic_energy_struct_models.png")
        if not figpath is None: plt.savefig(figpath+".png", dpi=256)
        if not figpath is None: plt.savefig(figpath+".pdf")
        # plt.savefig(sys.argv[0].replace(".py", "_") + name.lower() + ".pdf")
        plt.show()
        plt.close()







    #
    #
    #
    #
    # # eks = np.log10(eks)
    # fig, ax = plt.subplots(figsize=(4.6, 2.8), ncols=1, nrows=2, sharex="all")
    #
    # # fig.add_axes([0.6, .1, .35, .3])
    #
    # ax = ax[1]
    #
    # # ax = plt.subplot(111, polar=True)
    # levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
    # cmap = plt.get_cmap('RdYlBu_r')
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    # norm = LogNorm(eks[(eks>0)&(np.isfinite(eks))].min(), eks[(eks>0)&(np.isfinite(eks))].max())
    #
    # im = ax.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")
    #
    # ax.set_xlabel(r"Polar angle, $\theta$ [deg]")
    # ax.set_ylabel(r"$\Gamma\beta$")
    # ax.set_yscale("log")
    # ax.set_ylim(1e-1, 4)
    # if (not title is None):
    #     ax.set_title(title)
    #
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.minorticks_on()
    #
    # # ax.set_yscale("log")
    #
    # # cbar.set_title("x")
    # # ax.set_title('pcolormesh with levels')
    # # print(ax.get_rmin(), ax.get_rmax())
    # # ax.set_rmax(10)
    # # ax.set_rmin(20)
    # # print(ax.get_rmin(), ax.get_rmax())
    #
    # # max_theta = 90
    # # ax.set_thetamax(max_theta)
    # # ax.set_rlim(10,20)
    #
    # # ax.set_rlim(Rs.min(), Rs.max())
    # # ax.set_rscale("log")
    # # ax.set_rscale('log')
    # #
    # # ticklabels = ax.get_yticklabels()
    # # labels = range(80, 0, -10)
    # # for i in range(0, len(labels)):
    # #     ticklabels[i] = str(labels[i])
    # # ax.set_yticklabels(ticklabels)
    # plt.tight_layout()
    # if (save_figs): plt.savefig(FIGPATH + figname + ".png", dpi=256)
    # # if (save_figs): plt.savefig(PAPERPATH + figname + ".pdf")
    # plt.show()








class OLD_DATA():
    def get_rhomax(self):
        if not (os.path.isfile(self.fpath_rhomax)):
            print(f"File not found {self.fpath_rhomax}")
            return (np.zeros(0,), np.zeros(0,))

        rho_nuc = 2.7e14
        ut = 1.607e-6

        t_rho_max, rho_max = np.loadtxt(self.fpath_rhomax, unpack=True, usecols=(0,1))
        return (t_rho_max * ut * 1e3, rho_max / rho_nuc) # [s,rho/rho_nuc]

    def _read_mdot_table(self):
        names = None
        table = []
        try:
            table = np.loadtxt(self.fpath_mdot, unpack=True).T
            with open(self.fpath_mdot, "r") as file:
                names = file.readline()
        except:
            print("failed to read table. Trying line by line")
            names = None
            table = []
            n = 0
            with open(self.fpath_mdot, "r") as file:
                names = file.readline()
                lines = file.readlines()
                for line in lines:
                    vals = line.split()
                    # if len(vals) != len(names):
                    #     print(line)
                    #     raise ValueError(f"len(vals)={len(vals)} != len(names)={len(names)}")
                    nums = []
                    for val in vals:
                        try:
                            nums.append(np.float64(val))
                        except:
                            nums.append(np.nan)
                    table.append(nums)
                    n+=1
            table = np.reshape(np.array(table),newshape=(n,len(names)))
        print(table.shape)
        return (names, table)
    def get_mdot(self, r_ext=1000):
        if not (os.path.isfile(self.fpath_mdot)):
            print(f"File not found {self.fpath_mdot}")
            return (np.zeros(0,), np.zeros(0,))

        # convert the units from what Kenta was using to CGS + Msun
        ul = 0.4816 # km
        ut = 1.607e-6
        um = 0.326
        umdot = um / ut

        # table = np.loadtxt(self.fpath_mdot, unpack=True).T
        names, table = self._read_mdot_table()

        names = names.split()[1:]
        iv_n = lambda v_n : int(list(names).index(v_n)) if v_n in names else print("{} not found".format(v_n))
        get_mdot = lambda r_ext : table[:, iv_n("Mdot(rext={})".format(r_ext))] * umdot# Msun
        time = table[:, 0] * ut * 1e3 # -> s -> ms
        return (time, get_mdot(r_ext))

    def _read_mdot_table2(self):
        table = []
        try:
            table = np.loadtxt(self.fpath_mdot, unpack=True).T
        except:
            print("failed to read table. Trying line by line")
            with open(self.fpath_mdot, "r") as file:
                lines = file.readlines()
                for il, line in enumerate(lines):
                    vals = line.split()
                    nums = []
                    for val in vals:
                        try:
                            nums.append(np.float64(val))
                        except:
                            nums.append(np.nan)
                    table.append(nums)
                table = np.reshape(np.array(table),newshape=(il, len(nums)))
        return (table)

    def _read_mdot_table(self):
        if not (os.path.isfile(self.fpath_mdot)):
            print(f"File not found {self.fpath_mdot}")
            return (np.zeros(0,), np.zeros(0,))

        mdot_table = self._read_mdot_table2()

        # convert the units from what Kenta was using to CGS + Msun
        ul = 0.4816 # km
        ut = 1.607e-6
        um = 0.326
        umdot = um / ut
        mdot_table[:,0] *= (ut * 1e3) # -> s -> ms
        mdot_table[:,1:41] *= umdot # to Msun/s
        rext = np.array([50, 100, 300, 1000, 3000, 4000, 10000, 1229.3]) # 8 radii  in code unit.
        rext *= (ul * 1000) # code -> meters

        # create dict for data
        df = {
            "time":mdot_table[:,0],
        }
        i = 1
        for r in rext:
            df["mdot_r{:.1f}".format(r)] = mdot_table[:, i + len(rext) * 0] # [:,1:9] 2-9 Mdot_ext (8 extraction radii in code unit, all contribution)
            df["mdot_vr_r{:.1f}".format(r)] = mdot_table[:, i + len(rext) * 1] # [:,9:17] 10-17Mdot_ext2 (8 extraction radii in code unit, vr > 0.6c)
            df["mdot_slow_r{:.1f}".format(r)] = mdot_table[:, i + len(rext) * 2] # 18-25 Mdot_ext3 (8 extraction radii in code unit, \gamma_beta <= 0.1)
            df["mdot_mid_r{:.1f}".format(r)] = mdot_table[:, i + len(rext) * 3] # 26-33 Mdot_ext4 (8 extraction radii in code unit, 0.1 <= \gamma_beta < 1))
            df["mdot_fast_r{:.1f}".format(r)] = mdot_table[:, i + len(rext) * 3] # 34-41 Mdot_ext5 (8 extraction radii in code unit, \gamma_beta >= 1)
            i+=1
        # velocity
        df["vave_slow"] = mdot_table[:, 41] # 42 mass_averaged_velocity (\gamma_beta <= 0.1)
        df["vave_mid"]  = mdot_table[:, 42] # 43 mass_averaged_velocity (0.1 <= \gamma_beta < 1)
        df["vave_fast"] = mdot_table[:, 43] # 44 mass_averaged_velocity (\gamma_beta >= 1)

        # create dataframe
        self.mdot_df = pd.DataFrame.from_dict(data=df,orient="columns")
        self.mdot_df.set_index(keys="time")

    def get_mdot_at_r(self, r : np.float64, mode : str):
        if (not r in self.rext_mdot):
            idx = find_nearest_index(r,self.rext_mdot)
        else:
            idx = int(np.where(self.rext_mdot == r)[0])
        if mode == "vr>0.6":
            return self.table[:,1:9]    # 2-9 Mdot_ext (8 extraction radii in code unit, all contribution)

    def get_mdot2(self):
        if not (os.path.isfile(self.fpath_mdot)):
            print(f"File not found {self.fpath_mdot}")
            return (np.zeros(0,), np.zeros(0,))
        names, table = self._read_mdot_table2()



        print(table.shape)
        rext = np.array([50, 100, 300, 1000, 3000, 4000, 10000, 1229.3]) #  in code unit.
        time = table[:,0]           # 1 t (code unit)
        Mdot_ext  = table[:,1:9]    # 2-9 Mdot_ext (8 extraction radii in code unit, all contribution)
        Mdot_ext2 = table[:,9:17]   # 10-17Mdot_ext2 (8 extraction radii in code unit, vr > 0.6c)
        Mdot_ext3 = table[:,17:25]  # 18-25 Mdot_ext3 (8 extraction radii in code unit, \gamma_beta <= 0.1)
        Mdot_ext4 = table[:,25:33]  # 26-33 Mdot_ext4 (8 extraction radii in code unit, 0.1 <= \gamma_beta < 1))
        Mdot_ext5 = table[:,33:41]  # 34-41 Mdot_ext5 (8 extraction radii in code unit, \gamma_beta >= 1)
        vave      = table[:,41]     # 42 mass_averaged_velocity (\gamma_beta <= 0.1)
        vave1     = table[:,42]     # 43 mass_averaged_velocity (0.1 <= \gamma_beta < 1)
        vave3     = table[:,43]     # 44 mass_averaged_velocity (\gamma_beta >= 1)

        # convert the units from what Kenta was using to CGS + Msun
        ul = 0.4816 # km
        ut = 1.607e-6
        um = 0.326
        umdot = um / ut

        time *= (ut * 1e3) # -> s -> ms
        for arr in [Mdot_ext,Mdot_ext2,Mdot_ext3,Mdot_ext4,Mdot_ext5]:
            arr *= umdot # to Msun/s

        rext *= (ul * 1000) # code -> meters

        return (names, table)

def OLD_get_ej_data_for_text(files : list[str],
                         req_times=np.array([25]),
                         new_theta_len=None,
                         new_vinf_len=None,
                         verbose = True):
    """
        Load Kenta's data for various extraction times and get it
    """

    # sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    # files = sorted(glob.glob(datadir + "ejecta*", recursive=True), key=sort_by)
    if verbose: print("Ejecta files found: {}".format(files))
    vals = {"masses": [],
            "texts": [],
            "fast_masses": [],
            "v_ave": [],
            "v_fast_ave": [],
            "theta_rms":[],
            "fast_theta_rms":[]
            }
    pars_list = []

    #pb = PBA_TST2(do_ej_ele=True, do_ej_lc=True)
    for file in files:
        if verbose: print("\t Processing File {}".format(file))
        dfile = h5py.File(file, "r")
        if verbose: print("\t Keys in the dfile: {}".format(dfile.keys()))
        if verbose:print("\t theta              = {} [{}, {}]".format(np.array(dfile["theta"]).shape, np.array(dfile["theta"])[0],
                                                  np.array(dfile["theta"])[-1]))
        if verbose:print("\t v_asymptotic       = {} [{}, {}]".format(np.array(dfile["v_asymptotic"]).shape,
                                                  np.array(dfile["v_asymptotic"])[0],
                                                  np.array(dfile["v_asymptotic"])[-1]))
        if verbose:print("\t dfile['data1'].keys= {}".format(dfile["data1"].keys()))

        sort_by = lambda k: int(re.findall(r'\d+', k)[0]) if k.__contains__("data") else -1
        keys = sorted(dfile.keys(), key=sort_by)
        tkeys = [key for key in keys if key.__contains__("data")]
        times = []
        failed_to_get_time = []
        for key in keys:
            if key.__contains__("data"):
                try:
                    times.append(np.array(dfile[key]["time"], dtype=np.float64)[0])
                except KeyError:
                    failed_to_get_time.append(key)
                    print(f"Count not extract time from key={key} dfile[key]={dfile[key]}")
        print(f"Failed to extract time for {len(failed_to_get_time)}/{len(keys)}")
        # times = np.array(
        #     [np.array(dfile[key]["time"], dtype=np.float64)[0] for key in keys if key.__contains__("data")]
        # )
        # print("All times    = {}".format(times))
        unique_times = list(set(np.around(times, 0)))
        if len(unique_times) == 1:
            idxs = [0]
        else:
            idxs = [find_nearest_index(times, u_time) for u_time in unique_times]
        # print("Selected times: {}".format(times[idxs]))
        for idx in idxs:

            # extract arrays from h5 file (note that data may be in float32 not float64)
            v_inf = np.array(dfile["v_asymptotic"], dtype=np.float64)
            thetas = np.array(dfile["theta"], dtype=np.float64)
            mass = np.array(dfile[tkeys[idx]]["Mejecta"], dtype=np.float64)

            if verbose: print("Processing: time={} key={} tot_mass={}".format(times[idx], tkeys[idx], np.sum(mass)))

            # compute additional averaged quantities
            tot_mass = np.sum(mass)
            mass_fast = mass[:, v_inf * GammaFromBeta(v_inf) > 1]
            fast_ej_mass = np.sum(mass[:, v_inf * GammaFromBeta(v_inf) > 1])
            v_ave = np.sum(np.sum(mass, axis=0) * v_inf) / np.sum(mass)
            v_ave_fast = np.sum(np.sum(mass_fast, axis=0) * v_inf[v_inf * GammaFromBeta(v_inf) > 1]) / np.sum(mass_fast)
            theta_rms = (180. / np.pi) * np.sqrt(np.sum(np.sum(mass, axis=1) * thetas ** 2) / np.sum(mass))
            theta_rms_fast = (180. / np.pi) * np.sqrt(np.sum(np.sum(mass_fast, axis=1) * thetas ** 2) / np.sum(mass_fast))

            # collect data for each timestep, separate them by millisecond (integer)
            if (not int(times[idx]) in np.array(vals["texts"],dtype=np.int)):


                vals["texts"].append(times[idx])
                vals["masses"].append(tot_mass)
                vals["v_ave"].append(v_ave)
                vals["v_fast_ave"].append(v_ave_fast)
                vals["fast_masses"].append(fast_ej_mass)
                vals["theta_rms"].append(theta_rms)
                vals["fast_theta_rms"].append(theta_rms_fast)

                # ----------------------------------------

                # check if outdir exists
                # pars = copy.deepcopy(run_pars)
                # if (res_dir is None) and ("res_dir" in pars.keys()) and (not pars["res_dir"] is None):
                #     outdir = datadir + '/' + main_dir
                #     if not (os.path.isdir(outdir)):
                #         os.mkdir(outdir)
                #     outdir += res_dir
                #     if not (os.path.isdir(outdir)):
                #         os.mkdir(outdir)
                #
                #     pars["res_dir"] = outdir

                ''' create dict{} for runnin the model '''
                pars = {}
                pars["text"] = times[idx]

                # pars = {
                #     "n0": np.power(10,-3.51), "p": 2.05, "eps_e": 0.1, "eps_b": 1e-2, "eps_t":1.,
                #     "theta_obs": 30. * np.pi / 180, "timegrid": np.logspace(1., 6., 150) * cgs.day,
                #     "freq": [3e9], "z": 0.0099, "d_l": 41.3e6 * cgs.pc,
                #     "time": times[idx],
                #     "method_Up": "useEint2", "method_dgdr": "our",
                #     "method_shock_vel":"shockVel", "method_synchrotron":"Marg21",
                #     "method_lf_min":"useTheta","method_nonreldist":"useGm",
                #     "emissivity":"em", "absorption":"abs", "use_ssa":"yes",
                #     "t_arr": np.logspace(1., 6., 150) * cgs.day, "freq_arr": np.array([1e9, 3e9, 2.41e+17])
                #     "ejecta_prefix": pb.ejecta_prefix + "eps_t1_lfcut_tex{}_".format(int(times[idx]))
                # }

                # pars["method_Up"] = "useEint2"  # default "useGamma"
                # pars["method_dgdr"] = "our"  # default "peer"
                # pars["method_shock_vel"] = "shockVel"
                # pars["method_synchrotron"] = "Marg21"
                # pars["method_lf_min"] = "useTheta"
                # pars["method_nonreldist"] = "useGm"
                # pars["eps_t"] = 1.
                # pars["emissivity"] = emissivity
                # pars["absorption"] = "abs"
                # pars["use_ssa"] = "yes"
                # pars["t_arr"] = np.logspace(1., 6., 150) * cgs.day
                # pars["freq_arr"] = freqs  # np.geomspace(1e7, 1e15, 60)
                # # pars["p"] = 2.05
                # pars["n0"] = pb.test_gauss_jet["n0"]
                #pars["ejecta_prefix"] += "eps_t1_lfcut_tex{}_".format(int(times[idx]))
                #print(pars["ejecta_prefix"])

                # thetas, masses = reinterpolate_hist(thetas, masses[:-1, :], new_theta_len=new_theta_len)

                temp_ej = None
                if verbose: print(dfile[tkeys[idx]].keys())
                if ("T_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'T_ejecta'")
                    temp_ej_ = np.array(dfile[tkeys[idx]]["T_ejecta"], dtype=np.float64)
                    temp_ej_ *= 11604525006.17 # MeV -> Kelvin
                    v_inf_, thetas_, temp_ej = reinterpolate_hist2(v_inf, thetas, temp_ej_[:, :],
                                                                   new_theta_len=new_theta_len,
                                                                   new_vinf_len=new_vinf_len,
                                                                   mass_conserving=False)

                eps_ej = None
                if verbose: print(dfile[tkeys[idx]].keys())
                if ("internal_energy_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'internal_energy_ejecta'")
                    eps_ej_ = np.array(dfile[tkeys[idx]]["internal_energy_ejecta"], dtype=np.float64)
                    # temp_ej_ *= 11604525006.17 # MeV -> Kelvin
                    v_inf_, thetas_, eps_ej = reinterpolate_hist2(v_inf, thetas, eps_ej_[:, :],
                                                                   new_theta_len=new_theta_len,
                                                                   new_vinf_len=new_vinf_len,
                                                                   mass_conserving=False)

                press_ej = None
                if verbose: print(dfile[tkeys[idx]].keys())
                if ("pressure_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'pressure_ejecta'")
                    press_ej_ = np.array(dfile[tkeys[idx]]["pressure_ejecta"], dtype=np.float64)
                    # press_ej *= 11604525006.17 # MeV -> Kelvin
                    v_inf_, thetas_, press_ej = reinterpolate_hist2(v_inf, thetas, press_ej_[:, :],
                                                                   new_theta_len=new_theta_len,
                                                                   new_vinf_len=new_vinf_len,
                                                                   mass_conserving=False)
                entr_ej = None
                if verbose: print(dfile[tkeys[idx]].keys())
                if ("entropy_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'entropy_ejecta'")
                    entr_ej_ = np.array(dfile[tkeys[idx]]["entropy_ejecta"], dtype=np.float64)
                    # press_ej *= 11604525006.17 # MeV -> Kelvin
                    v_inf_, thetas_, entr_ej = reinterpolate_hist2(v_inf, thetas, entr_ej_[:, :],
                                                                    new_theta_len=new_theta_len,
                                                                    new_vinf_len=new_vinf_len,
                                                                    mass_conserving=False)

                ye_ej = None
                if ("Ye_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'Ye_ejecta'")
                    ye_ej_ = np.array(dfile[tkeys[idx]]["Ye_ejecta"], dtype=np.float64)
                    v_inf_, thetas_, ye_ej = reinterpolate_hist2(v_inf, thetas, ye_ej_[:, :],
                                                                 new_theta_len=new_theta_len,
                                                                 new_vinf_len=new_vinf_len,
                                                                 mass_conserving=False)

                rho_ej = None
                if ("rho_ejecta" in list(dfile[tkeys[idx]].keys())):
                    if verbose: print("Found 'rho_ejecta'")
                    rho_ej_ = np.array(dfile[tkeys[idx]]["rho_ejecta"], dtype=np.float64)
                    rho_ej_ *= 5.807e18 # Code units -> CGS
                    v_inf_, thetas_, rho_ej = reinterpolate_hist2(v_inf, thetas, rho_ej_[:, :],
                                                                  new_theta_len=new_theta_len,
                                                                  new_vinf_len=new_vinf_len,
                                                                  mass_conserving=False)


                masses = np.array(dfile[tkeys[idx]]["Mejecta"], dtype=np.float64)  # / cgs.solar_m
                v_inf, thetas, masses = reinterpolate_hist2(v_inf, thetas, masses[:, :],
                                                            new_theta_len=new_theta_len,
                                                            new_vinf_len=new_vinf_len,
                                                            mass_conserving=True)


                thetas = 0.5 * (thetas[1:] + thetas[:-1])
                v_inf  = 0.5 * (v_inf[1:] + v_inf[:-1])
                # v_inf = 0.5 * (v_inf[1:])
                # masses = masses[::-1, :]
                # masses = masses[thetas > 0, :]
                # thetas = thetas[thetas > 0]
                if verbose: print("Total mass:{}".format(np.sum(masses)))
                _betas = []
                _masses = []
                _temps = []
                _yes = []
                _rhos = []
                _press = []
                _eps = []
                _entr = []
                for iv in range(len(v_inf)):
                    if ((np.sum(masses[:, iv]) > 0.) & (v_inf[iv] > 0.) & (v_inf[iv] < 1.)):
                        _betas.append(v_inf[iv])
                        _masses.append(masses[:, iv])
                        if(not temp_ej is None): _temps.append(temp_ej[:, iv])
                        if(not ye_ej is None):_yes.append(ye_ej[:, iv])
                        if(not rho_ej is None):_rhos.append(rho_ej[:, iv])
                        if(not press_ej is None):_press.append(press_ej[:, iv])
                        if(not eps_ej is None):_eps.append(eps_ej[:, iv])
                        if(not entr_ej is None):_entr.append(entr_ej[:, iv])

                _betas = np.array(_betas)
                _masses = np.reshape(np.array(_masses), newshape=(len(_betas), len(thetas)))
                # ----------------------
                if(not temp_ej is None):
                    _temps = np.reshape(np.array(_temps), newshape=(len(_betas), len(thetas)))
                else:
                    _temps = np.zeros_like(_masses)
                if(not press_ej is None):
                    _press = np.reshape(np.array(_press), newshape=(len(_betas), len(thetas)))
                else:
                    _press = np.zeros_like(_masses)
                if(not ye_ej is None):
                    _yes = np.reshape(np.array(_yes), newshape=(len(_betas), len(thetas)))
                else:
                    _yes = np.zeros_like(_masses)
                if(not rho_ej is None):
                    _rhos = np.reshape(np.array(_rhos), newshape=(len(_betas), len(thetas)))
                else:
                    _rhos = np.zeros_like(_masses)
                if(not eps_ej is None):
                    _eps = np.reshape(np.array(_eps), newshape=(len(_betas), len(thetas)))
                else:
                    _eps = np.zeros_like(_masses)
                if(not entr_ej is None):
                    _entr = np.reshape(np.array(_entr), newshape=(len(_betas), len(thetas)))
                else:
                    _entr = np.zeros_like(_masses)
                # EjectaEk.compute_ek_corr(_betas, _masses)
                # print(_betas)
                # print(thetas)
                # print(_masses)
                pars["thetas"], pars["betas"], pars["masses"], pars["temp"], pars["ye"], pars["rho"], pars["press"], pars["eps"], pars["entr"] = \
                    thetas, _betas, _masses.T, _temps.T, _yes.T, _rhos.T, _press.T, _eps.T, _entr.T

                # if do_dist_plots:
                #     plot_init_profile(pars["thetas"], pars["betas"], pars["masses"].T,
                #                       figname=FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
                #                       title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(times[idx])))
                # print(np.diff(np.cos(pars["thetas"])))

                pars_list.append(copy.deepcopy(pars))

    if verbose: print("Total iterations : {}".format(len(pars_list)))
    if verbose: print("Total times      : {}".format(np.array(vals["texts"],dtype=int)))

    sorted_texts = np.sort(np.array(vals["texts"]))
    sorted_pars_list, sorted_vals = [], []
    sorted_vals = {"masses": [],
                   "texts": [],
                   "fast_masses": [],
                   "v_ave": [],
                   "v_fast_ave": [],
                   "theta_rms": [],
                   "fast_theta_rms": []
                   }
    for _time in sorted_texts:
        for i in range(len(pars_list)):
            if (_time == np.array(vals["texts"])[i]):
                for v_n in vals.keys():
                    sorted_vals[v_n].append(vals[v_n][i])
                sorted_pars_list.append(pars_list[i])
    assert len(pars_list) == len(sorted_pars_list)

    # -------------------------------------------------------------------

    # _distribute_and_run(pars_list[::2], n_cpu)

    if (not req_times is None):
        idxes = np.array(sorted(list(set([find_nearest_index(np.array(sorted_vals["texts"]), t) for t in req_times]))),
                         dtype=np.int64)
        selected_par_list = [sorted_pars_list[idx] for idx in idxes]
        selected_vals = {}
        for key in sorted_vals.keys():
            selected_vals[key] = [sorted_vals[key][idx] for idx in idxes]

        # print(np.array(sorted_vals["texts"])[idxes])
        if verbose: print("Selected iterations : {}".format(len(selected_par_list)))
        if verbose: print("Selected times      : {}".format(np.array([par["text"] for par in selected_par_list])))
        return (selected_par_list, selected_vals)
    else:
        return (sorted_pars_list, sorted_vals)


def OLD_prepare_kn_ej_id_2d(files : list[str],
                        outfpaths : list[str],
                        dist="pw",
                        req_times=np.array([25]),
                        new_theta_len=None,
                        new_vinf_len=None,
                        r0type="fromrho", t0=.1, r0frac=0.5,
                        verbose = True,
                        ):

    if (len(outfpaths) != len(req_times)):
        raise ValueError(" number of output files should be equal to requested times")

    if (dist != "pw"):
        raise NotImplementedError(" ID for other EATS methods are not available")

    selected_par_list, sorted_vals = \
        OLD_get_ej_data_for_text(files=files, req_times=req_times,
                             new_theta_len=new_theta_len, new_vinf_len=new_vinf_len, verbose = verbose)

    for pars, outfpath in zip(selected_par_list, outfpaths):
        theta_corr2 = pars["thetas"]
        vinf_corr2 = pars["betas"]
        ek_corr2 = compute_ek_corr(pars["betas"], pars["masses"]).T
        if verbose: print(theta_corr2.shape, vinf_corr2.shape, ek_corr2.shape)

        mass_corr = pars["masses"].T * cgs.solar_m
        temp_corr = pars["temp"].T
        ye_corr = pars["ye"].T
        rho_corr = pars["rho"].T
        entr_corr = pars["entr"].T
        eps_corr = pars["eps"].T
        press_corr = pars["press"].T

        if (r0type == "fromrho"):
            if (t0 < 0):
                raise ValueError(f" Not set t0={t0} must be > 0 in seconds for: r_base = t0 * cgs.c * vinf_corr2[0]")
            r = np.zeros_like(mass_corr)
            for ith in range(len(theta_corr2)):
                r_base = t0 * cgs.c * vinf_corr2[0]
                r[0,ith] = r_base
                for ir in range(1,len(vinf_corr2)):
                    if (mass_corr[ir,ith]==0. or rho_corr[ir,ith]==0.):
                        print(f"Error ir={ir} ith={ith} mass={mass_corr[ir,ith]} rho={rho_corr[ir,ith]}. Setting mass to 0.")
                        mass_corr[ir,ith]=0.
                        continue
                    r_i = (3./4.) * (1./np.pi) * len(theta_corr2) * mass_corr[ir,ith] / rho_corr[ir,ith] + r[ir-1,ith]**3
                    r_i = r_i**(1./3.)
                    r[ir,ith] = r_i
                    if ((r_i <= r[ir-1,ith]) or (~np.isfinite(r[ir,ith]))):
                        raise ValueError()
        elif (r0type == "frombeta"):
            r = np.zeros_like(mass_corr)
            t = t0
            for ith in range(len(theta_corr2)):
                for ir in range(len(vinf_corr2)):
                    r[ir,ith] = vinf_corr2[ir] * cgs.c * t # TODO THis is theta independent!
        else:
            raise KeyError(f"r0type={r0type} is not recognized")


        # fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(5,9),sharex='all')
        # im=axes[0].pcolormesh(theta_corr2* 180 / np.pi,vinf_corr2,mass_corr,norm=LogNorm(mass_corr[mass_corr>0].max()*1e-3,mass_corr[mass_corr>0].max()))
        # axes[0].set_title("Mass")
        # axes[0].set_ylabel("beta")
        # # axes[0].set_xticklabels(["{:.1f}".format(val) for val in theta_corr2 * 180 / np.pi])
        # # axes[0].set_yticklabels(["{:.2f}".format(val) for val in vinf_corr2])
        # # axes[0].set_yscale('log')
        # plt.colorbar(im)
        # im=axes[1].pcolormesh(theta_corr2* 180 / np.pi,vinf_corr2,rho_corr,norm=LogNorm(rho_corr[rho_corr>0].max()*1e-3,rho_corr[rho_corr>0].max()))
        # axes[1].set_title("rho")
        # axes[1].set_ylabel("beta")
        # plt.colorbar(im)
        # im=axes[2].pcolormesh(theta_corr2* 180 / np.pi,vinf_corr2,r,norm=LogNorm(r[r>0].max()*1e-3,r[r>0].max()))
        # axes[2].set_title("r")
        # axes[2].set_ylabel("beta")
        # axes[2].set_xlabel("theta")
        # plt.colorbar(im)
        # plt.show()


        ctheta_corr3 = np.zeros_like(ek_corr2)
        theta_corr3 = np.zeros_like(ek_corr2)
        for imom in range(len(vinf_corr2)):
            ctheta_corr3[imom,:] = theta_corr2
            theta_corr3[imom,:] = theta_corr2
        mom_corr3 = np.zeros_like(ek_corr2)
        for ith in range(len(theta_corr2)):
            mom_corr3[:,ith]=np.array( vinf_corr2*GammaFromBeta(vinf_corr2))

        # if (r0type == "fromrho"):
        #     r = np.zeros_like(mass_corr)
        #     for ith in range(len(ctheta_corr3[0,:])):
        #         idx = 0
        #         k = r0frac#0.5
        #         r[idx,ith] = (k*(3/4./np.pi)*rho_corr[idx,ith]*mass_corr[idx,ith])**(1./3.)
        #
        #         if (r[idx,ith] == 0):
        #             raise ValueError()
        #         for ir in range(1,len(mom_corr3[:,0]),1):
        #             _val = (3./4./np.pi)*rho_corr[ir,ith]*mass_corr[ir,ith]
        #             if (_val < 0):
        #                 raise ValueError(f"val={_val}")
        #             _rm1 = r[ir-1,ith]**3
        #             if(mass_corr[ir,ith]>0):
        #                 r[ir,ith] = (_rm1 + _val)**(1./3.)
        #             if ((r[ir-1,ith]>r[ir,ith])and(r[ir,ith]>0)):
        #                 raise ValueError(f"ir={ir} r[ir-1,ith]={r[ir-1,ith]} r[ir,ith]={r[ir,ith]}")
        # elif (r0type == "frombeta"):
        #     r = np.zeros_like(mass_corr)
        #     t = t0
        #     for ith in range(len(ctheta_corr3[0,:])):
        #         for ir in range(0,len(mom_corr3[:,0]),1):
        #             r[ir,ith] =  BetFromMom(mom_corr3[ir,ith])*cgs.c * t
        # else:
        #     raise KeyError(f"r0type={r0type} is not recognized")

        # check that radii are ordered

        # for i in range(len(ctheta_corr3[0,:])):
        #     for j in range(len(r[:,0])-1):
        #         if ((r[j+1,i] > 0) and (not (r[j+1,i] > r[j,i]))):
        #             print (f"i={i} j={j} and r[j+1,i]={r[j+1,i]} and r[j,i]={r[j,i]}")
        #             print(r)
        #             exit(1)
                # assert r[j+1,i] > r[j,i]

        # self.o_pba.setEjectaStructNumeric(theta_corr2, vinf_corr2, ek_corr2, fac, True, self.pars_kn["eats_method"])

        if verbose: print(len(theta_corr2), theta_corr2)
        dfile = h5py.File(outfpath, "w")
        dfile.create_dataset("r",data=r)
        dfile.create_dataset("theta",data=theta_corr3)
        dfile.create_dataset("ctheta",data=theta_corr3)
        dfile.create_dataset("mom",data=mom_corr3)
        dfile.create_dataset("ek",data=ek_corr2)
        dfile.create_dataset("mass",data=mass_corr)
        dfile.create_dataset("ye",data=ye_corr)
        dfile.create_dataset("rho",data=rho_corr)
        dfile.create_dataset("temp",data=temp_corr)
        dfile.create_dataset("press",data=press_corr)
        dfile.create_dataset("eps",data=eps_corr)
        dfile.create_dataset("entr",data=entr_corr)

        dfile.close()
        if verbose: print("file saved: {}".format(outfpath))

def OLD_load_init_data(fpath):
    dfile = h5py.File(fpath, "r")
    r_corr2 = np.array(dfile["r"],dtype=np.float64)
    theta_corr2 = np.array(dfile["theta"],dtype=np.float64)
    ctheta_corr2 = np.array(dfile["ctheta"],dtype=np.float64)
    mom_corr2 = np.array(dfile["mom"],dtype=np.float64)
    ek_corr2 = np.array(dfile["ek"],dtype=np.float64)
    mass_corr2 = np.array(dfile["mass"],dtype=np.float64)
    ye_corr2 = np.array(dfile["ye"],dtype=np.float64)
    rho_corr2 = np.array(dfile["rho"],dtype=np.float64)
    temp_corr2 = np.array(dfile["temp"],dtype=np.float64)
    press_corr2 = np.array(dfile["press"],dtype=np.float64)
    eps_corr2 = np.array(dfile["eps"],dtype=np.float64)
    entr_corr2 = np.array(dfile["entr"],dtype=np.float64)
    dfile.close()
    return (r_corr2, mom_corr2, theta_corr2, ctheta_corr2, ek_corr2, mass_corr2, ye_corr2, rho_corr2, temp_corr2, press_corr2, eps_corr2, entr_corr2)

def OLD_plot_init_profile(ctheta, betas, eks,
                      xmin=0,xmax=90,ymin=1e-2,ymax=6,vmin=1e-12,vmax=1e-6,
                      norm_mode="log", cmap = plt.get_cmap('RdYlBu_r'),
                      xscale="linear",yscale="linear",
                      subplot_mode="sum",
                      title=None, figpath=None):
    fontsize=12
    gammas = get_Gamma(betas)
    moms = gammas * betas
    ctheta = ctheta * 180 / cgs.pi

    fig = plt.figure(figsize=(4.5 + 1, 3.6 + 3))
    # fig.suptitle(r"BLh* $(1.259+1.482)M_{\odot}$")

    ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.59 - 0.12])
    ax1 = fig.add_axes([0.16, 0.61, 0.81 - 0.15, 0.91 - 0.61])
    cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.59 - 0.12])

    #                   x1    y1    delta x       delta y
    # ax1 = fig.add_axes([0.16, 0.45, 0.81 - 0.15, 0.59 - 0.12])
    # ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.91 - 0.61])
    # cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.91 - 0.61])

    # top panel
    # for t in tasks:
    #     hist_vinf, hist_mass = t["data"].load_vinf_hist()
    #     hist_mom = hist_vinf * get_Gamma(hist_vinf)
    #     hist_eks = np.cumsum(np.array(0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m)[::-1])[::-1]

    ax1.plot(ctheta, np.sum(eks,axis=0), color='black', ls='-', drawstyle='steps')

    if (not title is None): ax1.set_title(title)

    ###  mooley
    # mool_mom = np.linspace(0.5 * get_Gamma(0.5), 0.8 * get_Gamma(0.8), 20)
    # mool_ek = 5e50 * (mool_mom / 0.4) ** -5
    # _l, = ax1.plot(mool_mom, mool_ek, color="gray", ls="-")
    # ax1.text(0.30, 0.90, "Mooley+17", color='black', transform=ax1.transAxes, fontsize=12)

    # tasks = [1, 1, 1]
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="-", label=r"$q=1.00$")
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="--", label=r"$q=1.43$")
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls=":", label=r"$q=1.82$")
    # han, lab = ax1.get_legend_handles_labels()
    # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
    #                          **{"fancybox": False, "loc": 'lower left',
    #                            # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                            "shadow": "False", "ncol": 1, "fontsize": 11,
    #                            "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
    #                           **{"fancybox": False, "loc": 'upper right',
    #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # "shadow": "False", "ncol": 1, "fontsize": 11,
    # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # ax1.add_artist(ax1.legend(han[len(han) - len(tasks):], lab[len(lab) - len(tasks):],
    #                           **{"fancybox": False, "loc": 'center right',
    #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # "shadow": "False", "ncol": 1, "fontsize": 11,
    # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # fit *= np.sum(mdens) / np.sum(fit)
    # fit_x, fit_y = _fit(mom, ekdens)
    # ax1.plot(fit_x, fit_y, 'k--', label=r'$\sin^2(\theta)$')
    # ax1.legend(**{"fancybox": False, "loc": 'upper right',
    #                # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                "shadow": "False", "ncol": 1, "fontsize": 11,
    #                "framealpha": 0., "borderaxespad": 0., "frameon": False})

    if norm_mode=="log":
        ax1.set_yscale("log")
    else:
        ax1.set_yscale("linear")
    # ax1.set_ylim(ymin, ymax)
    # ax1.yaxis.tick_right()
    # ax1.yaxis.tick_left()
    ax1.set_ylabel(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    ax1.get_yaxis().set_label_coords(-0.15, 0.5)

    ax1.set_xlim(xmin, xmax)
    ax1.xaxis.set_ticklabels([])

    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.minorticks_on()

    ax11 = ax1.twinx()
    mask = betas*get_Gamma(betas)>1
    if len(betas[mask])>0:
        axx = np.sum(betas[mask, np.newaxis] * eks[mask, :],axis=0)
        ayy = np.sum(eks[mask, :],axis=0)
        # ax11.plot(ctheta, axx/ayy, color="gray", ls="-")
        if (subplot_mode=="sum"): _yarr = np.sum(ctheta*eks[mask, :],axis=0)
        elif (subplot_mode=="ave"): _yarr = axx/ayy
        else:raise KeyError("subplot_mode is not recognized")
        ax11.plot(ctheta, np.sum(eks[mask, :],axis=0), color="gray", ls="-")
    if norm_mode=="log":
        ax11.set_yscale("log")
    else:
        ax11.set_yscale("linear")
    ax11.set_ylabel(r"$M_{\rm ej}(\Gamma\beta>1)$ [M$_{\odot}$]", fontsize=fontsize, color="gray")
    ax11.minorticks_on()
    ax11.tick_params(axis='both', which='both', labelleft=False,
                     labelright=True, tick1On=False, tick2On=True,
                     labelsize=12,
                     direction='in',
                     bottom=True, top=True, left=True, right=True)


    # bottom panel
    # mask = vinf > 0.6
    # eks2 = np.zeros_like(eks)
    # for i in range(len(vinf)):
    #     if vinf[i] > 0.6:
    #         eks2[:, i] = eks[:, i]

    # import h5py
    # dset = h5py.File("/home/vsevolod/Desktop/Hajela/BLh_q100.h5", 'w')
    # dset.create_dataset(name="beta", data=vinf)
    # dset.create_dataset(name="theta", data=theta * 180 / np.pi)
    # dset.create_dataset(name="Ek(>Gamma_beta)", data=eks)
    # dset.close()

    if (norm_mode=="log"):
        norm = LogNorm(eks[(eks > 0) & (np.isfinite(eks))].min(), eks[(eks > 0) & (np.isfinite(eks))].max())
    elif (norm_mode=="linear"):
        norm = Normalize(eks[(eks > 0) & (np.isfinite(eks))].min(), eks[(eks > 0) & (np.isfinite(eks))].max())
    elif (norm_mode=="levels"):
        levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    else:
        raise KeyError(" norm_mode is not recognized ")


    im = ax0.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

    ax0.axhline(y=1, linestyle='--', color='gray')

    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)
    ax0.set_xscale(xscale)
    ax0.set_yscale(yscale)

    ax0.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
    ax0.set_xlabel(r"Polar angle", fontsize=fontsize)
    ax0.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, lbl, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    # ax0.text(0.05, 0.05, models.print_fancy_label(name), color='white', transform=ax0.transAxes)

    # ekdens = np.sum(eks, axis=0)
    # ax1.step(0.5 * (theta[1:] + theta[:-1]), mdens, where='mid', color='red')
    # hist_vinf, hist_mass = o_data.load_vinf_hist()
    # hist_mom = hist_vinf * get_Gamma(hist_vinf)
    # hist_eks = 0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m
    # ax1.step(mom, ekdens, where='mid', color='red')
    # ax1.step(hist_mom, hist_eks, where='mid', color='black')

    cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar.set_label(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.ax.minorticks_off()
    # plt.savefig(PAPERPATH + "kinetic_energy_struct_models.pdf")
    # plt.savefig(FIGPATH + "kinetic_energy_struct_models.png")
    if not figpath is None: plt.savefig(figpath+".png", dpi=256)
    if not figpath is None: plt.savefig(figpath+".pdf")
    # plt.savefig(sys.argv[0].replace(".py", "_") + name.lower() + ".pdf")
    plt.show()
    plt.close()







    #
    #
    #
    #
    # # eks = np.log10(eks)
    # fig, ax = plt.subplots(figsize=(4.6, 2.8), ncols=1, nrows=2, sharex="all")
    #
    # # fig.add_axes([0.6, .1, .35, .3])
    #
    # ax = ax[1]
    #
    # # ax = plt.subplot(111, polar=True)
    # levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
    # cmap = plt.get_cmap('RdYlBu_r')
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    # norm = LogNorm(eks[(eks>0)&(np.isfinite(eks))].min(), eks[(eks>0)&(np.isfinite(eks))].max())
    #
    # im = ax.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")
    #
    # ax.set_xlabel(r"Polar angle, $\theta$ [deg]")
    # ax.set_ylabel(r"$\Gamma\beta$")
    # ax.set_yscale("log")
    # ax.set_ylim(1e-1, 4)
    # if (not title is None):
    #     ax.set_title(title)
    #
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.minorticks_on()
    #
    # # ax.set_yscale("log")
    #
    # # cbar.set_title("x")
    # # ax.set_title('pcolormesh with levels')
    # # print(ax.get_rmin(), ax.get_rmax())
    # # ax.set_rmax(10)
    # # ax.set_rmin(20)
    # # print(ax.get_rmin(), ax.get_rmax())
    #
    # # max_theta = 90
    # # ax.set_thetamax(max_theta)
    # # ax.set_rlim(10,20)
    #
    # # ax.set_rlim(Rs.min(), Rs.max())
    # # ax.set_rscale("log")
    # # ax.set_rscale('log')
    # #
    # # ticklabels = ax.get_yticklabels()
    # # labels = range(80, 0, -10)
    # # for i in range(0, len(labels)):
    # #     ticklabels[i] = str(labels[i])
    # # ax.set_yticklabels(ticklabels)
    # plt.tight_layout()
    # if (save_figs): plt.savefig(FIGPATH + figname + ".png", dpi=256)
    # # if (save_figs): plt.savefig(PAPERPATH + figname + ".pdf")
    # plt.show()


def plot2(vals : dict, figpath = None):
    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(4.6,2+2.8), sharex="all")
    ax = axes[0]

    ax.plot(np.array(vals["texts"]), np.array(vals["v_ave"])*get_Gamma(vals["v_ave"]), 'x', color='blue')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=False,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_ylabel(r"Mass-averaged $\langle \Gamma\beta \rangle$", color="blue")
    # ax.set_xlabel(r"$t_{\rm ext}$ [ms]")

    ax1 = ax.twinx()
    ax1.plot(np.array(vals["texts"]), np.array(vals["v_fast_ave"])*get_Gamma(vals["v_fast_ave"]), 'o', color="red")
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=False, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax1.set_ylabel(r"Mass-averaged $\langle \Gamma\beta(\Gamma\beta>1) \rangle$",color="red")
    ax1.set_xlim(vals["texts"][0], vals["texts"][-1])
    # plt.tight_layout()
    # plt.savefig(FIGPATH+figppath+"average_velocity_evolution"+"png",dpi=256)
    # plt.show()

    # -------------------------------------------------------------------
    # fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4.6, 2.8))
    ax = axes[1]
    ax.plot(np.array(vals["texts"]), np.array(vals["theta_rms"]), 'x', color='blue', label="Total ejecta")
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=False,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_ylabel(r"RMS half-openning angle $\langle \theta_{\rm RMS} \rangle$", color="red")

    ax.plot(np.array(vals["texts"]), np.array(vals["fast_theta_rms"]), 'd', color="blue", label=r"Fast tail, $\Gamma\beta>1$")
    ax.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax.legend(**{"fancybox": False, "loc": 'lower right',
                 # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                 "shadow": "False", "ncol": 1, "fontsize": 12 - 2, "columnspacing": 0.4,
                 "framealpha": 0., "borderaxespad": 0., "frameon": False})
    ax1 = ax.twinx()
    ax1.plot(np.array(vals["texts"]), np.array(vals["fast_masses"]), 'x', color="red")
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=False, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.set_ylabel(r"$\langle M_{\rm ej}(\Gamma\beta>1)$", color="red")
    # ax1.set_yscale("log")
    ax.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax.set_ylabel(r"$\theta_{\rm RMS} = \sqrt{ \Sigma(m_i\theta_i)/\Sigma(m_i) }$", color="blue")
    ax.set_xlim(vals["texts"][0], vals["texts"][-1])
    plt.tight_layout()
    if not figpath is None: plt.savefig(figpath + ".png", dpi=256)
    if not figpath is None: plt.savefig(figpath + ".pdf")
    plt.show()


