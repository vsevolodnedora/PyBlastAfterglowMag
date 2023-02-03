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
from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
from package.src.PyBlastAfterglowMag.interface import cgs
from package.src.PyBlastAfterglowMag.utils import latex_float

afterglowpy = True

try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

def tst_against_afgpy(withSpread = False,
                      savefig = "compare_uniform_afgpy.png",
                      load_data = True):

    pba = PBA(workingdir=os.getcwd()+"/", readparfileforpaths=True)

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes

    # pba_0 = PBA(os.getcwd()+"/", readparfileforpaths=True)
    # pba_016 = PBA(os.getcwd()+"/", readparfileforpaths=True)

    lls, lbls = [], []
    for (i_thetaobs, i_freq, i_color) in [
        # (thetaObs, freqobs, "blue"),
        (0.16, 1e9, "orange"),
        (0, 1e18, "green"),
        (0.16, 1.e18, "red"),
        (0, 1e9, "gray"),
    ]:

        pba.modify_main_part_parfile(newpars={"theta_obs":i_thetaobs},newopts={})
        pba.modify_grb_part_parfile(newpars={},newopts={"method_synchrotron":"Joh06",
                                                        "fname_light_curve":"tophat_{}_joh06.h5".format(
                                                            str(i_thetaobs).replace(".",""))})
        pba.reload_parfile()

        pba.run()

        ax.plot(pba.get_jet_lc_times() / cgs.day,
                pba.get_jet_lc_totalflux(freq=i_freq), color=i_color, ls='-',
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        pba.clear()
        pba.modify_main_part_parfile(newpars={"theta_obs":i_thetaobs},newopts={})
        pba.modify_grb_part_parfile(newpars={},newopts={"method_synchrotron":"WSPN99",
                                                        "fname_light_curve":"tophat_{}_WSPN99.h5".format(
                                                            str(i_thetaobs).replace(".",""))})
        pba.reload_parfile()

        pba.run()

        ax.plot(pba.get_jet_lc_times() / cgs.day,
                pba.get_jet_lc_totalflux(freq=i_freq), color=i_color, ls=':',
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        pba.clear()

        if load_data:
            if withSpread:
                fname = "afterglowpy_theta{:d}_lognu{:d}_spread.txt".format(int(i_thetaobs * 180 / np.pi),
                                                                            int(np.log10(i_freq)))
            else:
                fname = "afterglowpy_theta{:d}_lognu{:d}.txt".format(int(i_thetaobs * 180 / np.pi),
                                                                     int(np.log10(i_freq)))
            _t, _ref_F_afgpy, _ref_F = np.loadtxt(os.getcwd() + '/' + fname, unpack=True)
            _ll, = ax.plot(_t / cgs.day, _ref_F_afgpy, color=i_color, ls='--')
            lls.append(_ll)
            lbls.append(r"$\nu=$" + r"${}$ Hz ".format(latex_float(i_freq))
                        + r"$\theta_{\rm obs}=$" + r"{:.1f} deg".format(i_thetaobs * 180 / np.pi))

        # break

    l11, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='-', label=r"\texttt{PBA} with J06")
    l12, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls=':', label=r"\texttt{PBA} with WSPN99")
    l13, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"\texttt{afterglowpy}")

    legend1 = plt.legend([l11, l12, l13],
                         # [r"\& J\'{o}hannesson+06", r"\& WSPN+99", r"\texttt{afterglowpy}"],
                         [r"\texttt{PBA} with J06", r"\texttt{PBA} with WSPN99", r"\texttt{afterglowpy}"],
                         loc="center", bbox_to_anchor=(0.78, 0.56), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)
    legend2 = plt.legend(lls, lbls,
                         loc="center", bbox_to_anchor=(0.4, 0.16), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")

    # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   # labelsize=plotdic["fontsize"],
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [mJy]")
    ax.set_xlim(1e-1, 1e3)
    ax.set_ylim(1e-9, 1e2)
    # ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(PAPERPATH + save)
    if savefig: plt.savefig(os.getcwd() + '/' + savefig, dpi=256)

    plt.show()

if __name__ == '__main__':
    tst_against_afgpy()
    exit(0)

















def main():

    # afterglowpy
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

    # pyblastafterglow
    pba = PBA(workingdir=os.curdir+'/',readparfileforpaths=True)
    print(pba.main_pars)
    pba.modify_main_part_parfile(newpars=dict({"theta_obs":0.16}),newopts={})
    pba.modify_grb_part_parfile(newpars={},newopts={"fname_light_curve":"grb_lc_theta016.h5"})
    pba.run()
    exit(1)
    print(pba.main_pars)
    print(pba.get_jet_lc_totalflux(freq=1e9))

    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"./afgpy_grb170817.txt",unpack=True)
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

    tts, fluxes = np.loadtxt(curdir+"./jelib_grb170817.txt",unpack=True)
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')

    ax.plot(pba.get_jet_lc_times(), pba.get_jet_lc_totalflux(freq=3.e9), ls='-', color="black", label="PBA")

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
