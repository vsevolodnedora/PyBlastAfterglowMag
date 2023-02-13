from multiprocessing import Pool
import contextlib
from numpy import *
from numpy.linalg import solve, det
from numpy.random import uniform
from scipy.interpolate import interp1d
import emcee
import os
import copy
import hashlib
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from shutil import copyfile

import numpy as np


try:
    from PyBlastAfterglowMag.interface import (BPA_METHODS, modify_parfile_par_opt)
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import (BPA_METHODS, modify_parfile_par_opt)
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import latex_float, cgs
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

class GRB170817A(object):
    """ observational data collected by Gavin Lamb """

    def __init__(self):

        # load data
        times, freqs, data, errs = self._load_data()

        # unique freqs in the data
        ufreqs = np.unique(freqs)

        # concatenate the data for unique frequencies
        # data_ord = np.array([])
        # err_ord = np.array([])
        #
        # for frequency in freqs:
        #     # print(frequency)
        #     data_ord = np.concatenate([data_ord, data[freq == frequency]])
        #     err_ord = np.concatenate([err_ord, err[freq == frequency]])
        # assert np.shape(err_ord) == np.shape(data_ord)

        print("--- Observations for GRB170817 ---")
        print("Total amount of freq: {}".format(len(freqs)))
        print("Total amount of time: {}".format(len(times)))
        print("Total amount of data: {}".format(len(data)))
        print("Total amount of errs: {}".format(len(errs)))
        print("Unique frequencies   ({})".format(len(ufreqs)))
        print(ufreqs)
        print("Data per frequency:")
        for ifreq, freq in enumerate(ufreqs):
            print("freq={:.2e} N={:d}".format(freq, len(data[freq==freqs])))

        self.times = times
        self.data = data
        self.errs = errs
        self.freqs = freqs
        self.ufreqs = ufreqs

    def __call__(self, freq):
        return self.get(freq)

    def _load_data(self):
        """
            Chandra21 :: https://ui.adsabs.harvard.edu/abs/2020GCN.29055....1H/abstract
            Hajela:2020
        :return:
        """

        # Updated optical and X-ray data from Fong ea 2019, and Hajela ea 2019
        time = np.array([9.2, 14.9, 16.4, 17.4, 18.3, 18.7, 19.4, 21.4, 22.4,
                         23.4, 24.2, 31.3, 35.3, 39.2, 46.3, 53.3, 54.3, 57.2,
                         65.9, 66.6, 67.2, 72.2, 75.5, 75.5, 77.6, 79.2, 80.1,
                         92.4, 93.1, 93.1, 93.2, 97.1, 107., 107., 109., 109.,
                         111., 112., 115., 115., 115., 125., 125., 126., 133.,
                         137., 149., 150., 152., 158., 161., 163., 163., 163.,
                         163., 165., 167., 170., 172., 183., 197., 197., 207.,
                         209., 216., 217., 217., 217., 217., 218., 218., 222.,
                         229., 252., 257., 259., 261., 267., 267., 273., 273.,
                         289., 289., 294., 297., 298., 320., 324., 328., 357.,
                         359., 362., 380., 489., 545., 580., 581., 741., 767.,
                         938.,
                         1211 # Chandra21
                         ])  # time in days

        # flux in mJy
        data = np.array(
            [5.66e-04, 6.58e-04, 1.87e+01, 1.51e+01, 1.45e+01, 1.54e+01, 1.59e+01, 1.36e+01, 2.25e+01, 2.00e+01,
             2.56e+01,
             3.40e+01, 4.40e+01, 2.28e+01, 4.40e+01, 3.20e+01, 4.80e+01, 6.10e+01,
             1.48e+02, 9.80e+01, 4.26e+01, 5.80e+01, 3.59e+01, 3.96e+01, 7.70e+01, 4.50e+01, 4.17e+01, 3.17e+01,
             9.80e+01,
             7.00e+01, 2.60e+01, 1.99e+02, 1.27e+02, 5.32e+01, 2.96e-03, 1.09e-01,
             1.11e-01, 6.29e+01, 9.62e+01, 5.12e+01, 4.12e+01, 5.82e+01, 1.28e+02, 2.21e+02, 3.37e-03, 8.40e-02,
             6.06e+01,
             9.00e+01, 1.84e+02, 3.03e-03, 2.66e-03, 9.73e+01, 6.73e+01, 4.74e+01,
             3.96e+01, 9.10e-02, 5.79e+01, 1.13e-01, 8.50e-02, 2.11e+02, 7.59e+01, 8.93e+01, 4.20e+01, 8.20e-02,
             3.63e+01,
             6.05e+01, 4.17e+01, 3.26e+01, 2.47e+01, 6.47e+01, 6.30e-02, 3.97e+01,
             4.80e+01, 7.13e+01, 4.32e+01, 1.55e-03, 6.26e+01, 2.50e+01, 4.03e+01, 3.48e+01, 2.72e+01, 3.63e+01,
             2.70e+01,
             3.12e+01, 4.40e-02, 2.34e+01, 2.31e+01, 4.72e+01, 3.40e-02, 9.70e-04,
             1.55e+01, 2.70e-02, 3.79e+01, 1.48e+01, 5.90e+00, 1.80e+01, 3.54e-04, 2.68e-04, 4.90e+00, 1.95e-04,
             3.46e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
             ])

        # frequency of obs.
        freq = np.array([2.41e+17, 2.41e+17, 3.00e+09, 3.00e+09, 3.00e+09, 7.25e+09,
                         6.20e+09, 6.20e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         1.50e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         6.70e+08, 1.30e+09, 6.00e+09, 4.50e+09, 7.35e+09, 7.35e+09,
                         1.40e+09, 4.50e+09, 6.00e+09, 7.25e+09, 1.50e+09, 3.00e+09,
                         1.50e+10, 6.70e+08, 1.30e+09, 1.30e+09, 2.41e+17, 3.80e+14,
                         5.06e+14, 6.00e+09, 3.00e+09, 1.00e+10, 1.50e+10, 7.25e+09,
                         1.30e+09, 6.70e+08, 2.41e+17, 5.06e+14, 7.25e+09, 5.10e+09,
                         1.30e+09, 2.41e+17, 2.41e+17, 3.00e+09, 6.00e+09, 1.00e+10,
                         1.50e+10, 5.06e+14, 7.25e+09, 3.80e+14, 5.06e+14, 6.50e+08,
                         3.00e+09, 1.30e+09, 5.00e+09, 5.06e+14, 1.00e+10, 3.00e+09,
                         6.00e+09, 1.00e+10, 1.50e+10, 3.00e+09, 5.06e+14, 7.25e+09,
                         4.50e+09, 1.30e+09, 3.00e+09, 2.41e+17, 1.30e+09, 7.25e+09,
                         3.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09,
                         5.06e+14, 7.25e+09, 7.25e+09, 1.30e+09, 5.06e+14, 2.41e+17,
                         7.25e+09, 5.06e+14, 1.30e+09, 3.00e+09, 6.00e+09, 7.25e+09,
                         2.41e+17, 2.41e+17, 3.00e+09, 2.41e+17,
                         2.41e+17 # Chandra21
                         ])

        # error on flux
        err = np.array([1.70e-04, 1.30e-04, 6.30e+00, 3.90e+00, 3.70e+00, 4.80e+00,
                        5.50e+00, 2.90e+00, 3.40e+00, 3.10e+00, 2.90e+00, 3.60e+00,
                        1.00e+01, 2.60e+00, 4.00e+00, 4.00e+00, 6.00e+00, 9.00e+00,
                        2.20e+01, 2.00e+01, 4.10e+00, 5.00e+00, 4.30e+00, 7.00e+00,
                        1.90e+01, 7.00e+00, 4.70e+00, 4.30e+00, 1.40e+01, 5.70e+00,
                        4.40e+00, 1.60e+01, 1.80e+01, 4.50e+00, 2.60e-04, 1.70e-02,
                        1.90e-02, 3.20e+00, 8.00e+00, 3.40e+00, 1.90e+00, 5.00e+00,
                        2.10e+01, 1.90e+01, 4.00e-04, 1.80e-02, 4.30e+00, 3.00e+01,
                        1.90e+01, 2.60e-04, 2.70e-04, 1.13e+01, 4.10e+00, 3.60e+00,
                        2.00e+00, 1.60e-02, 6.90e+00, 1.90e-02, 1.70e-02, 3.40e+01,
                        5.20e+00, 1.39e+01, 1.20e+01, 2.00e-02, 3.60e+00, 7.50e+00,
                        7.50e+00, 4.00e+00, 3.10e+00, 2.70e+00, 1.80e-02, 7.20e+00,
                        6.00e+00, 6.70e+00, 5.80e+00, 1.90e-04, 7.00e+00, 4.10e+00,
                        2.70e+00, 4.90e+00, 2.10e+00, 3.90e+00, 2.80e+00, 3.60e+00,
                        1.40e-02, 4.20e+00, 4.00e+00, 1.28e+01, 1.10e-02, 1.90e-04,
                        5.00e+00, 7.00e-03, 1.18e+01, 2.90e+00, 1.90e+00, 4.20e+00,
                        9.00e-05, 9.00e-05, 1.80e+00, 7.00e-05,
                        1.3e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
                        ])

        assert np.shape(time) == np.shape(freq)
        assert np.shape(data) == np.shape(err)

        # data *= 1.e-29 # muJy -> ergs
        # err *= 1.e-29 # muJy -> ergs

        return(time, freq, data, err) # [days, Hz, ergs, ergs]

    def get(self, freq):
        if not freq in self.ufreqs:
            raise NameError("freq: {} is not in data. Available are: {}".format(freq, self.ufreqs))
        mask = freq == self.freqs
        if len(mask) < 1:
            raise ValueError("Failed to get data for freq:{}".format(freq))
        return (self.times[mask], self.data[mask], self.errs[mask])

    def get_chandra(self):
        # Time since Merger (days) Fluxd (mJy) Err (mJy)
        data = np.array(
            [[2.19000, 1.40000e-07,  0.00000,     0.00000],
             [9.19679,  2.10000e-07,  8.72000e-08, -9.02000e-08],
             [15.3900,  6.44000e-07,  1.72000e-07, -1.20000e-07],
             [108.386,  2.21000e-06,  2.14000e-07, -2.02000e-07],
             [157.755,  2.41000e-06,  2.62000e-07, -2.11000e-07],
             [259.665,  1.07000e-06,  1.64000e-07, -1.55000e-07],
             [358.609,  9.12000e-07,  2.10000e-07, -1.73000e-07],
             [581.818,  2.14000e-07,  1.09000e-07, -7.72000e-08],
             [741.478,  1.26000e-07,  6.64000e-08, -5.39000e-08],
             [939.310,  1.54977e-07,  8.46022e-08, -6.26098e-08],
             [1234.11,  2.16000e-07,  5.45000e-08, -7.95000e-08]
             ])
        # return [day] [erg] [erg] [erg]
        errs = np.abs(data[:, 2:][:,::-1]).T
        # errs[:, 0] += data[:, 1]
        # errs[:, 1] += data[:, 1]
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, errs / 1e3 / 1e23, uplims)##,  data[:, 3] * 1e3 * 1e23)

    def get_vla_3ggz(self):

        data = np.array(
            [
                # Time since Merger (days) Fluxd (mJy) Err (mJy)
                [3.34000,     0.032000002,       0.0000000],
                [16.4200,     0.018700000,    0.0063000000],
                [17.3900,     0.015100000,    0.0038999999],
                [18.3300,     0.014500000,    0.0037000000],
                [22.3600,     0.022500001,    0.0034000000],
                [24.2600,     0.025599999,    0.0029000000],
                [31.3200,     0.034000002,    0.0035999999],
                [46.2600,     0.044000000,    0.0040000002],
                [54.2700,     0.048000000,    0.0060000001],
                [57.2200,     0.061000001,    0.0089999996],
                [93.1300,     0.070000000,    0.0057000001],
                ##[196.790     0.081793795,    0.0081793800],
                [196.790,     0.073213301,    0.0066756521],
                [115.200,     0.075690837,     0.019037671],
                [115.200,      0.10307770,     0.011835645],
                ##[163.070,     0.096567894,     0.020509181],
                [163.070,     0.098128255,     0.018721838],
                [216.910,     0.068999998,     0.015000000],
                [220.000,     0.064700000,    0.0027000001],
                [256.760,     0.055000000,     0.012300000],
                [267.000,     0.040300000,    0.0027000001],
                [272.670,     0.043900002,     0.010500000],
                [288.610,     0.046399999,     0.011400000],
                [294.000,     0.031199999,    0.0035999999],
                [489.000,     0.014800000,    0.0029000000],
                [767.000,    0.0049000001,    0.0018000000],
                # [1221.50,    0.0054000001,      0.00000000],
                [1243.000,   2.86*1e-3,          0.99*1e-3], # https://arxiv.org/pdf/2103.04821.pdf
                [4.6*cgs.year/cgs.day,   4.5*1e-3,          1.1*1e-3] # https://arxiv.org/pdf/2205.14788.pdf
            ]
        )
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        # data[:, 2] = 0.00200000
        uplims[-1] = True
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, data[:, 2] / 1e3 / 1e23, uplims)

def compute_likelihood(P, v_ns, scales, obs_data):
    pars = {}
    for i, v_n in enumerate( v_ns ):
        pars[v_n] = P[i]

    times = np.array(obs_data["times"])
    freqs = np.array(obs_data["freqs"])

    if pars["theta_c"] > pars["theta_w"]:
        # return (np.array([-np.inf for i in range(len(times))]), -np.inf)
        return -np.inf

    if "theta_c" in pars.keys():
        pars["theta_c"] =pars["theta_c"] * np.pi / 180.
    if "theta_w" in pars.keys():
        pars["theta_w"] = pars["theta_w"] * np.pi / 180.
    if "theta_obs" in pars.keys():
        pars["theta_obs"] = pars["theta_obs"] * np.pi / 180.

    for v_n, scale in zip(v_ns, scales):
        if scale == "log_uniform":
            pars[v_n] = 10.**pars[v_n]
    # pars["Eiso_c"] = 10**pars["Eiso_c"]
    # pars["eps_e"] = 10**pars["eps_e"]
    # pars["eps_b"] = 10**pars["eps_b"]
    # pars["n_ism"] = 10**pars["n_ism"]

    if (pars["Gamma0c"]) < 1.:
        raise ValueError()

    # assert P["Eiso_c"] > 1e45
    # assert P["Gamma0c"] > 10.
    # assert P["n_ism"] > 1e-5
    # assert P["Eiso_c"] > 1e45
    # assert P["Eiso_c"] > 1e45

    workdir = os.getcwd()+'/'+"working/"
    newparfilenameroot = "EisoC[Eiso_c]_Gamma0c[Gamma0c]_thetaC[theta_c]_thetaW[theta_w]_" \
                         "theta[theta_obs]_nism[n_ism]_p[p]_epse[eps_e]_epsb[eps_b]_"
    # we need to make a parfile for a given parameter set. Reuse the old code to when we iterate over values
    iter_pars_dict = copy.deepcopy(pars)
    for key in iter_pars_dict.keys():
        iter_pars_dict[key] = [iter_pars_dict[key]]
    iter_pars = list(iter_pars_dict.keys())
    new_pars = set_parlists_for_pars(iter_pars_keys=iter_pars,
                                     iter_pars=iter_pars_dict,
                                     fname=newparfilenameroot)
    if len(new_pars) > 1:
        raise ValueError("Size of the pars is > 1. Should not happen")
    new_pars = new_pars[0]

    hash = hashlib.sha1(str.encode("gauss_"+new_pars["name"]+"parfile.par")).hexdigest()

    # set parameters for the model
    main_pars = copy.deepcopy(new_pars)
    grb_pars = copy.deepcopy(new_pars)
    grb_opts = {"fname_light_curve": "tophat_"+hash+"lc.h5"}
    # parfile_name = "tophat_"+new_pars["name"]+"parfile.par"
    parfile_name = "tophat_"+hash+".par"

    # separate parameters between main and grb ones
    main_keys = ["n_ism","theta_obs"]
    grb_keys = ["Eiso_c", "Gamma0c", "theta_c", "theta_w", "p", "eps_e", "eps_b"]
    for key in grb_keys:
        main_pars.pop(key)
    for key in main_keys:
        grb_pars.pop(key)

    # 'raw_array' that one time corresponds to one frequency (arrys or equal length)]
    main_opts = {}
    main_opts["lc_freqs"] = "array " + "".join(["{:.2e} ".format(val) for val in freqs])
    main_opts["lc_times"] = "array " + "".join(["{:.2e} ".format(val) for val in times])

    if os.path.isfile(parfile_name):
        raise FileExistsError("parfile already exists: {}".format(parfile_name))

    # create the parfile for this run
    modify_parfile_par_opt(part="main", newpars=main_pars, newopts=main_opts,
                           workingdir=workdir, parfile="parfile.par", newparfile=parfile_name, keep_old=True)
    modify_parfile_par_opt(part="grb", newpars=grb_pars, newopts=grb_opts,
                           workingdir=workdir, parfile=parfile_name, newparfile=parfile_name, keep_old=False)

    # run the code. The code may fail for some reason, so save the parfile that led to failure
    try:
        pba = BPA_METHODS(workingdir=workdir, readparfileforpaths=True, parfile=parfile_name)
        pba.run(loglevel="err")
        dfile = pba.get_jet_lc_obj()
        _times = np.array(dfile["times"]) # [s]
        _freqs = np.array(dfile["freqs"]) # [Hz]
        _fluxes =  np.array(dfile["fluxes"]) # [mJy]
        dfile.close()
        pba.clear()
    except:
        print("PBA failed. Saving parfile for future analysis: {}".format(workdir+"failed_"+parfile_name))
        copyfile(workdir+parfile_name, workdir+"failed_"+parfile_name)
        return -np.inf

    if len(times)!=len(_times):
        raise ValueError("model licht curve has {} times while input data has {} times".format(len(_times),len(times)))
    if len(freqs)!=len(_freqs):
        raise ValueError("model licht curve has {} freqs while input data has {} freqs".format(len(_times),len(times)))

    # os.remove(workdir+"tmp*")
    os.remove(workdir+parfile_name)
    os.remove(workdir+grb_opts["fname_light_curve"])

    # TODO make loading skymap and Xc calcs.

    # data_err = _fluxes
    # delxmu = data_err - pars["obs_fluxes"]
    # nparams =
    # covMatrix = np.diag(data_err*data_err)
    # loglikelihood = -0.5*np.dot(delxmu, np.solve(covMatrix,delxmu)) -0.5*nparams*np.log(2.*pi) - 0.5*np.log(np.det(covMatrix))

    # times = np.array(pars["times"])
    # freqs = np.array(pars["freqs"])
    obs_fluxes = obs_data["fluxes"]
    obs_fluxes_errs = obs_data["flux_errs"]
    loglikelihood = -0.5*(np.sum(((obs_fluxes-_fluxes)/obs_fluxes_errs)**2.))#+np.log(err**2)
    # fig, axes = plt.subplots(ncols=1, nrows=2)
    # axes[0].loglog( obs_data["times"], ((obs_fluxes-_fluxes)/obs_fluxes_errs)**2., marker='.', color='black', label=f'{loglikelihood}')
    # axes[1].loglog( obs_data["times"], _fluxes, marker='.', color='black', label=f'{loglikelihood}')
    # axes[1].loglog( obs_data["times"], obs_fluxes, marker='.', color='gray', label=f'{loglikelihood}')
    # axes[0].legend()
    # axes[1].set_ylim(1e-7,1e-3)
    # # plt.show()
    # plt.savefig(workdir+hash+".png",dpi=256)
    # plt.close(fig)
    return loglikelihood


    # xM = 0.
    # return (_fluxes, xM)

class MyFit():
    def __init__(self, workingdir, settings, obs_data):

        self.nwalkers = 50

        self.ranges =   [par["range"]   for par in settings]
        self.v_ns =     [par["name"]    for par in settings]
        self.scales =   [par["prior"]   for par in settings]
        self.labels =   [par["label"]   for par in settings]
        self.nparameters = len(self.v_ns)
        self.input_pars = settings
        self.obs_data = obs_data

        p0 = np.zeros([self.nwalkers,self.nparameters],np.float64)
        # Initial parameter guess
        print ("Initial parameter guesses")
        seed = np.random.get_state()[1][0]
        for i, set_v_n in enumerate(settings):
            p0[:,i] = set_v_n["guess"] + np.random.randn(self.nwalkers) * set_v_n["mult"]
            # print(set_v_n["name"], " ", set_v_n["guess"], " ", set_v_n["mult"], p0[:,i])
            # continue
        # print("-"*20)
        # print(repr(p0))
        # print("-"*20)
        # p0[:, 0] = 51 + np.random.randn(self.nwalkers)*0.1  # sample E1
        # p0[:, 1] = 100 + np.random.randn(nwalkers)*10  # sample G1
        # p0[:, 2] = 0.12 + np.random.randn(nwalkers)*0.01  # sample thc1
        # p0[:, 3] = 0.91 + np.random.randn(nwalkers)*0.01  # sample cos(incl1)
        # p0[:, 4] = -2.1 + np.random.randn(nwalkers)*0.1  # sample EB1
        # p0[:, 5] = -2.1 + np.random.randn(nwalkers)*0.1  # sample Ee1
        # p0[:, 6] = -3.1 + np.random.randn(nwalkers)*0.1  # sample n1
        # p0[:, 7] = 2.16 + np.random.randn(nwalkers)*0.01  # sample p

        # Multiprocess
        fpath = workingdir+"backend-refreshed-AG_all.h5" #THIS IS A SPECIAL FILE THAT MEANS YOU CAN (I CAN SHOW YOU) RESTART A RUN WHERE YOU LEFT OFF
        print(f"Creating backend file {fpath}")

        backend = emcee.backends.HDFBackend(fpath)
        print("Initial size: {0}".format(backend.iteration))

        # backend.reset(self.nwalkers, self.nparameters)

        burnsteps = 1000
        nsteps = 10000

        sampler = None
        with Pool(processes=12) as pool:
        # pool = Pool()
            sampler = emcee.EnsembleSampler(self.nwalkers, self.nparameters, self.lnprob,
                                            pool=pool, backend=backend, runtime_sortingfn=self.sort_on_runtime)
            ""
            # Burn in run
            print ("Burning in")
            pos, prob, state = sampler.run_mcmc(p0, burnsteps, progress=True)
            sampler.reset()
            ""
            # Production run
            print ("Production run")
            sampler.run_mcmc(pos, nsteps, progress=True)
            # Result['Chain'] = finalpos.chain
            # Result['LnProbability'] = finalprob # self._Sampler.lnprobability
            # Result['AcceptanceFraction'] =finalstate.acceptance_fraction
            # # TheChain = Result['Chain']
            # # LnProbability = Result['LnProbability']
            # samples = (sampler.get_chain(discard=50, flat=True))

        # Write full MCMC to file
        fchain = workingdir+"chain.csv"
        print("Writing full chain into file: {}".format(fchain))
        with open(fchain, 'w') as f:
            f.write(f"{self.nparameters}, {self.nwalkers}, {nsteps}\n")
            for j in range(nsteps):
                for i in range(self.nwalkers):
                    for k in range(self.nparameters):
                        f.write(f"{sampler.chain[i, j, k]:.6f}, ")
                    f.write(f"{sampler.lnprobability[i, j]:.6f}\n")

        # Write each individual parameter to it's own file
        for k in range(self.nparameters):
            print("Writing chain for: {}".format(self.input_pars[k]["name"]))
            fname = workingdir+"chain_{}.csv".format(self.input_pars[k]["name"])
            with open(fname, 'w') as f:
                for j in range(nsteps):
                    for i in range(self.nwalkers):
                        if i == (self.nwalkers-1):
                            f.write(f"{sampler.chain[i, j, k]:.6f}\n")
                        else:
                            f.write(f"{sampler.chain[i,j,k]:.6f}, ")

        # Write probability to it's own file
        fname = workingdir+"lnp.csv"
        print(f"Writing probability file: {fname}")
        with open(fname, 'w') as f:
            for j in range(nsteps):
                for i in range(self.nwalkers):
                    if i == (self.nwalkers-1):
                        f.write(f"{sampler.lnprobability[i, j]:.6f}\n")
                    else:
                        f.write(f"{sampler.lnprobability[i,j]:.6f}, ")

        # Acceptance fraction and autocorrelation time
        try:
            tau = sampler.get_autocorr_time()
            print(
                f"Mean acceptance fraction: {np.mean(sampler.acceptance_fraction)}\n" +
                f"Average auto-correlation time: {np.mean(tau):.3f}"
            )
            fname = workingdir+"info.json"
            print(f"Writing file: {fname}")
            info = {
                "Npars": self.nparameters,
                "Nwalk": self.nwalkers,
                "Nstep": int(nsteps),
                "seed": int(seed),
                "MeanAcceptanceFraction": np.mean(sampler.acceptance_fraction),
                "AverageAuto-correlationTime": np.mean(tau)
            }
            # info["acceptance_fraction"] = np.mean(sampler.acceptance_fraction)
            # info["tau"] = tau
            with open(fname, "w") as f:
                json.dump(info, f)
        except:
            print("Failed to compute tau or save the file")
        try:
            fname = workingdir+"trace.pdf"
            print(f"Plotting trace plot: {fname}")
            self.create_trace_plot(sampler, self.nparameters, nsteps, self.nwalkers, fname, self.labels)
        except:
            print("Failed to plot the result")
        print("Finished")
        # print("Saving {}")
        # np.savetxt(workingdir+'samples.txt', samples)
        # cPickle.dump(samples, open("samples.pkl", "w"))

    def initialize_walkers(self):
        initializations = zeros([self.nwalkers, self.nparameters])
        for ii in range(self.nparameters):
            initializations[:,ii] = uniform(self.ranges[ii][0], self.ranges[ii][1], self.nwalkers)
        return initializations

    def lnprior(self, P):
        for i, set in enumerate(self.input_pars):
            if (P[i] < self.ranges[i][0]) or (P[i] > self.ranges[i][1]):
                # print(f"P[i]={P[i]} < min={self.ranges[i][0]} OR p[i]={P[i]} > max={self.ranges[i][1]}")
                return -np.inf
        return 0.

    def lnprob(self, P):
        lp = self.lnprior(P)
        if np.isinf(lp):
            return -np.inf
        loglikelihood = compute_likelihood(P, self.v_ns, self.scales, self.obs_data)
        return lp + loglikelihood

    @staticmethod
    def sort_on_runtime(pos):
        """
        Function to sort chain runtimes at execution.
        """
        p = np.atleast_2d(pos)
        idx = np.argsort(p[:, 0])[::-1]

        return (p[idx], idx)

    @staticmethod
    def create_trace_plot(sampler, Npars, Nstep, Nwalk, fplot, names):
        fig, axes = plt.subplots(Npars+1, 1, sharex=True, figsize=(6,8))

        for i in range(Nwalk):
            axes[0].plot(range(Nstep), sampler.get_log_prob()[:, i], c='gray', alpha=0.4)
        axes[0].yaxis.set_major_locator(MaxNLocator(4, prune='lower'))
        axes[0].tick_params(axis='both', which='major', labelsize=10)
        axes[0].set_ylabel('$\ln (p)$', fontsize=12)

        for i in range(Npars):
            for j in range(Nwalk):
                axes[i + 1].plot(range(Nstep), sampler.get_chain()[:, j, i], c='gray', alpha=0.4)
            axes[i + 1].yaxis.set_major_locator(MaxNLocator(4, prune='lower'))
            axes[i + 1].tick_params(axis='both', which='major', labelsize=10)
            axes[i + 1].set_ylabel(names[i], fontsize=12)

        axes[-1].set_xlabel('Model Number', fontsize=12)
        fig.tight_layout(h_pad=0.1)
        fig.savefig('{0}'.format(fplot))
        plt.clf()

def main():
    # get observational data
    workingdir = os.getcwd()+"/working/"
    gw = GRB170817A()
    data_all = {"times":gw.times*cgs.day, # [s]
                "freqs":gw.freqs,         # [Hz]
                "fluxes":gw.data,  # [mJy]
                "flux_errs":gw.errs # [mJy]
                }

    settings = [
        {"name": "Eiso_c",      "range":[50., 53.],       "guess":52,     "mult":0.1,   "prior":"log_uniform",      "label":r"$E_{\rm iso}$"},
        {"name": "Gamma0c",     "range":[100.,1000.],   "guess":100.,   "mult":0.1,   "prior":"linear_uniform",   "label":r"$\Gamma_0$"},
        {"name": "theta_c",     "range":[3., 30.],      "guess":5.,     "mult":0.1,   "prior":"linear_uniform",   "label":r"$\theta_c$"},
        {"name": "theta_w",     "range":[5., 30.],      "guess":16.,    "mult":0.1,   "prior":"linear_uniform",   "label":r"$\theta_w$"},
        {"name": "n_ism",       "range":[-4., 0.],       "guess":-3.2,   "mult":0.1,   "prior":"log_uniform",      "label":r"$n_{\rm ISM}$"},
        {"name": "theta_obs",   "range":[0., 30.],      "guess":21.,    "mult":0.1,   "prior":"linear_uniform",   "label":r"$\theta_{\rm obs}$"},
        {"name": "eps_e",       "range":[-4.,-0.5],      "guess":-2,     "mult":0.1,   "prior":"log_uniform",      "label":r"$\epsilon_{e}$"},
        {"name": "eps_b",       "range":[-4.,-0.5],      "guess":-3,     "mult":0.1,   "prior":"log_uniform",      "label":r"$\epsilon_{b$"},
        {"name": "p",           "range":[2.05, 2.30],   "guess":2.16,   "mult":0.1,   "prior":"linear_uniform",   "label":r"$p$"}
    ]
    fit = MyFit( workingdir, settings, data_all )

if __name__ == '__main__':
    main()