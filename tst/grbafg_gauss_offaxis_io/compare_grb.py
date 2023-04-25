# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

# from PyBlastAfterglowMag import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import PyBlastAfterglow
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
    from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
        from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

afterglowpy = True

try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"


def main():
    workdir = os.getcwd()+'/'
    pba = PyBlastAfterglow(workingdir=os.curdir+'/',readparfileforpaths=True)

    prepare_grb_ej_id_1d({"struct":"gaussian",
        "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
        "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10}, type="pw",
        outfpath=workdir+"gauss_grb_id.h5")

    pba.run()
    # print(pba.fpath_kn_light_curve)
    # print(pba.get_jet_lc_totalflux(freq=1e9))

    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"./afgpy_grb170817.txt",unpack=True)
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

    tts, fluxes = np.loadtxt(curdir+"./jelib_grb170817.txt",unpack=True)
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')

    ax.plot(pba.GRB.get_lc_times(), pba.GRB.get_lc_totalflux(freq=3.e9), ls='-', color="black", label="PBA")
    print(pba.GRB.get_lc_obj().keys())
    print(pba.GRB.get_lc_obj().attrs.keys())
    for il in range(int(pba.GRB.get_lc_obj().attrs["nlayers"])):
        ax.plot(pba.GRB.get_lc_times(), pba.GRB.get_lc(freq=3.e9, ishell=0, ilayer=il), ls='-', color="gray")

    print(pba.GRB.get_lc_totalflux(freq=3.e9))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()

if __name__ == '__main__':

    main()



def tst_tophat_old():
    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    fname = "test_data_grb/lcs_tophat.h5"

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    # colnames, data = load_data(curdir + fname)
    # nx, ny = data.shape
    # t_arr = data[:, 0]
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    ax.plot(our_times, our_fluxes, ls='-', color='black', label='PyBlastAfterglow [pw]')
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()

def compareTopHatLightCurves():

    ''' Run 'compareTopHatLightCurves()' first '''
    fname = "test_data_grb/lcs_tophat.h5"
    if (afterglowpy) : Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
                            'specType':    0,                  # Basic Synchrotron Spectrum
                            'counterjet':  0,
                            'spread':      0,
                            'thetaObs':    0.0,   # Viewing angle in radians
                            'E0':          1.0e52, # Isotropic-equivalent energy in erg
                            'g0':          1000,
                            'thetaCore':   0.2,    # Half-opening angle in radians
                            'thetaWing':   0.2,
                            'n0':          1e-3,    # circumburst density in cm^{-3}
                            'p':           2.2,    # electron energy distribution index
                            'epsilon_e':   0.1,    # epsilon_e
                            'epsilon_B':   0.01,   # epsilon_B
                            'xi_N':        1.0,    # Fraction of electrons accelerated
                            'd_L':         3.09e26, # Luminosity distance in cm
                            'z':           0.0099}   # redshift

    t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9
    if (afterglowpy) : Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


    if (afterglowpy) :
        ax.plot(t, Fnu, color='gray', label='afterglowpy')

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    ax.plot(our_times, our_fluxes, color="black", ls="-", label="PyBlastAfterglow [PW]")
    #
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    #         break

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Top-hat [1GGz]")
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()

    # -------------------------------------------------

    plt.show()
# compareTopHatLightCurves()

def compareTopHatLightCurvesMethods():

    ''' Run 'compareTopHatLightCurves()' first '''
    fnames = {
        "test_data_grb/lcs_tophat_a_obs.h5":"blue",
        "test_data_grb/lcs_tophat_pw_obs.h5":"cyan",
        "test_data_grb/lcs_tophat_a_comov.h5":"green",
        "test_data_grb/lcs_tophat_pw_comov.h5":"lime"
    }
    if (afterglowpy) : Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
                            'specType':    0,                  # Basic Synchrotron Spectrum
                            'counterjet':  0,
                            'spread':      0,
                            'thetaObs':    0.0,   # Viewing angle in radians
                            'E0':          1.0e52, # Isotropic-equivalent energy in erg
                            'g0':          1000,
                            'thetaCore':   0.2,    # Half-opening angle in radians
                            'thetaWing':   0.2,
                            'n0':          1e-3,    # circumburst density in cm^{-3}
                            'p':           2.2,    # electron energy distribution index
                            'epsilon_e':   0.1,    # epsilon_e
                            'epsilon_B':   0.01,   # epsilon_B
                            'xi_N':        1.0,    # Fraction of electrons accelerated
                            'd_L':         3.09e26, # Luminosity distance in cm
                            'z':           0.0099}   # redshift

    t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9
    if (afterglowpy) : Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


    if (afterglowpy) :
        ax.plot(t, Fnu, color='gray', label='afterglowpy')

    for fname in fnames.keys():
        dfile = h5py.File(curdir+fname, "r")
        for key in dfile.keys():
            print(key)
        our_times = np.array( dfile["times"] )
        our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

        ax.plot(our_times, our_fluxes, color=fnames[fname], ls="-", label=fname.split("/")[-1])
    #
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    #         break

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Top-hat [1GGz]")
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()

    # -------------------------------------------------

    plt.show()
compareTopHatLightCurvesMethods()

def compareGaussianOffAxisGflatSStructLightCurves():
    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    fname = "test_data_grb/lcs_FernadOffAxisGflatSSstruct_methods.h5"

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    # colnames, data = load_data(curdir + fname)
    # nx, ny = data.shape
    # t_arr = data[:, 0]
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    ax.plot(our_times, our_fluxes, ls='-', color='black', label='PyBlastAfterglow [pw]')
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()
# compareGaussianOffAxisGflatSStructLightCurves()

def compareGaussianOffAxisGflatSStructLightCurvesMethods():

    fnames = {
        "test_data_grb/lcs_FernandOffAxisGauss_a_obs.h5":"blue",
        "test_data_grb/lcs_FernandOffAxisGauss_pw_obs.h5":"cyan",
        "test_data_grb/lcs_FernandOffAxisGauss_a_comov.h5":"green",
        "test_data_grb/lcs_FernandOffAxisGauss_pw_comov.h5":"lime"
    }

    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    # fname = "test_data_grb/lcs_FernadOffAxisGflatSSstruct_methods.h5"

    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')

    for fname in fnames.keys():
        dfile = h5py.File(curdir+fname, "r")
        for key in dfile.keys():
            print(key)
        our_times = np.array( dfile["times"] )
        our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

        # colnames, data = load_data(curdir + fname)
        # nx, ny = data.shape
        # t_arr = data[:, 0]
        # icolor1, icolor2 = 0, 0
        # for icol in range(1,len(colnames)):
        #     name = colnames[icol]
        #     row = data[:,icol]
        #     if name.__contains__("[A]"):
        #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
        #         icolor1 += 1
        #     elif name.__contains__("[P]"):
        #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
        #         icolor2 += 1

        ax.plot(our_times, our_fluxes, ls='-', color=fnames[fname], label=fname.split("/")[-1])
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()
compareGaussianOffAxisGflatSStructLightCurvesMethods()
