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
from scipy import interpolate
import copy
import glob
import re
import os
import sys
import argparse

from .id_maker_tools import (reinterpolate_hist, compute_ek_corr)
from .utils import (cgs, get_Gamma, get_Beta, find_nearest_index)


def get_ej_data_for_text(datadir : str,
                         req_times=np.array([25]),
                         new_theta_len=None,
                         verbose = True):
    """
        Load Kenta's data for various extraction times and get it
    """

    sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    files = sorted(glob.glob(datadir + "ejecta*", recursive=True), key=sort_by)
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
        times = np.array(
            [np.array(dfile[key]["time"], dtype=np.float64)[0] for key in keys if key.__contains__("data")]
        )
        # print("All times    = {}".format(times))
        unique_times = set(np.around(times, 0))
        idxs = [find_nearest_index(times, u_time) for u_time in unique_times]
        print("Selected times: {}".format(times[idxs]))
        for idx in idxs:

            # extract arrays from h5 file (note that data may be in float32 not float64)
            v_inf = np.array(dfile["v_asymptotic"], dtype=np.float64)
            thetas = np.array(dfile["theta"], dtype=np.float64)
            mass = np.array(dfile[tkeys[idx]]["Mejecta"], dtype=np.float64)

            if verbose: print("Processing: time={} key={} tot_mass={}".format(times[idx], tkeys[idx], np.sum(mass)))

            # compute additional averaged quantities
            tot_mass = np.sum(mass)
            mass_fast = mass[:, v_inf * get_Gamma(v_inf) > 1]
            fast_ej_mass = np.sum(mass[:, v_inf * get_Gamma(v_inf) > 1])
            v_ave = np.sum(np.sum(mass, axis=0) * v_inf) / np.sum(mass)
            v_ave_fast = np.sum(np.sum(mass_fast, axis=0) * v_inf[v_inf * get_Gamma(v_inf) > 1]) / np.sum(mass_fast)
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

                masses = np.array(dfile[tkeys[idx]]["Mejecta"], dtype=np.float64)  # / cgs.solar_m
                thetas, masses = reinterpolate_hist(thetas, masses[:-1, :], new_theta_len=new_theta_len)
                thetas = 0.5 * (thetas[1:] + thetas[:-1])
                # v_inf = 0.5 * (v_inf[1:])
                # masses = masses[::-1, :]
                # masses = masses[thetas > 0, :]
                # thetas = thetas[thetas > 0]
                print("Total mass:{}".format(np.sum(masses)))
                _betas = []
                _masses = []
                for iv in range(len(v_inf)):
                    if ((np.sum(masses[:, iv]) > 0.) & (v_inf[iv] > 0.) & (v_inf[iv] < 1.)):
                        _betas.append(v_inf[iv])
                        _masses.append(masses[:, iv])
                _betas = np.array(_betas)
                _masses = np.reshape(np.array(_masses), newshape=(len(_betas), len(thetas)))
                # EjectaEk.compute_ek_corr(_betas, _masses)
                print(_betas)
                print(thetas)
                print(_masses)
                pars["thetas"], pars["betas"], pars["masses"] = thetas, _betas, _masses.T

                # if do_dist_plots:
                #     plot_init_profile(pars["thetas"], pars["betas"], pars["masses"].T,
                #                       figname=FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
                #                       title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(times[idx])))
                # print(np.diff(np.cos(pars["thetas"])))

                pars_list.append(copy.deepcopy(pars))

    print("Total iterations : {}".format(len(pars_list)))
    print("Total times      : {}".format(np.array(vals["texts"],dtype=int)))

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

        print(np.array(sorted_vals["texts"])[idxes])
        print("Selected iterations : {}".format(len(selected_par_list)))
        print("Selected times      : {}".format(np.array([par["text"] for par in selected_par_list])))
        return (selected_par_list, selected_vals)
    else:
        return (sorted_pars_list, sorted_vals)

def prepare_kn_ej_id_2d(datadir : str,
                        outfpaths : list[str],
                        dist="pw",
                        req_times=np.array([25]),
                        new_theta_len=None,
                        verbose = True,
                        ):

    if (len(outfpaths) != len(req_times)):
        raise ValueError(" number of output files should be equal to requested times")

    if (dist != "pw"):
        raise NotImplementedError(" ID for other EATS methods are not available")

    selected_par_list, sorted_vals = \
        get_ej_data_for_text(datadir=datadir, req_times=req_times, new_theta_len=new_theta_len, verbose = verbose)

    for pars, outfpath in zip(selected_par_list, outfpaths):
        theta_corr2 = pars["thetas"]
        vinf_corr2 = pars["betas"]
        ek_corr2 = compute_ek_corr(pars["betas"], pars["masses"]).T
        if verbose: print(theta_corr2.shape, vinf_corr2.shape, ek_corr2.shape)

        # self.o_pba.setEjectaStructNumeric(theta_corr2, vinf_corr2, ek_corr2, fac, True, self.pars_kn["eats_method"])

        if verbose: print(len(theta_corr2), theta_corr2)
        dfile = h5py.File(outfpath, "w")
        dfile.create_dataset("theta",data=theta_corr2)
        dfile.create_dataset("vel_inf",data=vinf_corr2)
        dfile.create_dataset("ek",data=ek_corr2)
        dfile.close()
        if verbose: print("file saved: {}".format(outfpath))
