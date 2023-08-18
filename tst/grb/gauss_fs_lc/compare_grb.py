# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.cm import ScalarMappable
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



try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"


def main():
    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    # ----------------------- piece-wise --------------------
    workdir = os.getcwd()+'/'
    prepare_grb_ej_id_1d({"struct":"gaussian",
                          "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
                          "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10}, type="pw",
                         outfpath=workdir+"gauss_grb_id.h5")
    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                           newpars={},
                           newopts={"method_eats":"piece-wise", "method_spread":"AA"},
                           parfile="parfile.par", newparfile="parfile.par", keep_old=False)

    pba_pw = PyBlastAfterglow(workingdir=os.curdir+'/',readparfileforpaths=True, parfile="parfile.par")
    pba_pw.run()

    ax.plot(pba_pw.GRB.get_lc_times(), pba_pw.GRB.get_lc_totalflux(freq=3.e9), ls='-', color="green", label="PBA [PW]")
    # print(pba_pw.GRB.get_lc_obj().keys())
    # print(pba_pw.GRB.get_lc_obj().attrs.keys())
    cmap = cm.get_cmap('Blues')
    norm = Normalize(vmin=0,vmax=int(pba_pw.GRB.get_lc_obj().attrs["nlayers"]))
    for il in range(int(pba_pw.GRB.get_lc_obj().attrs["nlayers"])):
        ax.plot(pba_pw.GRB.get_lc_times(), pba_pw.GRB.get_lc(freq=3.e9, ishell=0, ilayer=il), ls='-', color=cmap(norm(il)), alpha=0.7)


    # -------------------- ADAPTIVE ------------------
    # workdir = os.getcwd()+'/'
    prepare_grb_ej_id_1d({"struct":"gaussian",
                          "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
                          "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10}, type="a",
                         outfpath=workdir+"gauss_grb_id.h5")
    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                           newpars={},
                           newopts={"method_eats":"adaptive", "method_spread":"AFGPY"},
                           parfile="parfile.par", newparfile="parfile.par", keep_old=False)

    pba_a = PyBlastAfterglow(workingdir=os.curdir+'/',readparfileforpaths=True, parfile="parfile.par")
    pba_a.run()

    ax.plot(pba_a.GRB.get_lc_times(), pba_a.GRB.get_lc_totalflux(freq=3.e9), ls='-', color="red", label="PBA [A]")
    # print(pba_pw.GRB.get_lc_obj().keys())
    # print(pba_pw.GRB.get_lc_obj().attrs.keys())
    cmap = cm.get_cmap('Reds')
    norm = Normalize(vmin=0,vmax=int(pba_a.GRB.get_lc_obj().attrs["nlayers"]))
    for il in range(int(pba_a.GRB.get_lc_obj().attrs["nlayers"])):
        ax.plot(pba_a.GRB.get_lc_times(), pba_a.GRB.get_lc(freq=3.e9, ishell=0, ilayer=il), ls='-', color=cmap(norm(il)), alpha=0.7)


    # print(pba_pw.fpath_kn_light_curve)
    # print(pba_pw.get_jet_lc_totalflux(freq=1e9))

    # --------- Other Models -----------
    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"./afgpy_grb170817.txt",unpack=True)
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

    tts, fluxes = np.loadtxt(curdir+"./jelib_grb170817.txt",unpack=True)
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')


    # print(pba_pw.GRB.get_lc_totalflux(freq=3.e9))
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
# compareTopHatLightCurvesMethods()

def compareTopHatLightCurvesMethods():

    ''' Run 'compareTopHatLightCurves()' first '''
    fnames = {
        "test_data_grb/lcs_tophat_a_obs.h5":"blue",
        "test_data_grb/lcs_tophat_pw_obs.h5":"cyan",
        "test_data_grb/lcs_tophat_a_comov.h5":"green",
        "test_data_grb/lcs_tophat_pw_comov.h5":"lime"
    }
    Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
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
    Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


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
