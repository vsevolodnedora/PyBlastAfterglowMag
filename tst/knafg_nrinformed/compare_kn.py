# import PyBlastAfterglowMag
import shutil

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

from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
from package.src.PyBlastAfterglowMag.interface import cgs, modify_parfile_par_opt
from package.src.PyBlastAfterglowMag.id_maker_from_thc_ourflow import prepare_kn_ej_id_1d, prepare_kn_ej_id_2d

def main():
    # prepare ejecta ID from NR output
    workdir = os.getcwd()+'/'
    shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes
    pba = PBA(workingdir=os.getcwd()+'/',readparfileforpaths=True)

    tasks = [
        {"corr_fpath":"corr_vel_inf_theta_DD2_M135135_M0.h5",
         "hist_fpath":"hist_vel_inf_DD2_M135135_M0.dat",
         "outfpath":os.getcwd()+'/'+"ejecta_id_DD2_M135135_M0.h5",
         "new_main_pars":{"n_ism":1e-1,"d_l":100e6 * cgs.pc,"z":0.001},
         "new_kn_pars":{"p":2.5,"eps_e":1.e-1,"eps_b":1e-1,"theta_obs":0.0},
         "fname_light_curve":"lc_DD2_M135135_M0.h5",
         "fname_ejecta_id":"ejecta_id_DD2_M135135_M0.h5",
         "fpath_kenta_res":os.getcwd()+'/'+"DD2_n01_100Mpc.dat",
         "color":"blue","label":"DD2 M135135 M0"
         },
        # {"hist_fpath":"corr_vel_inf_theta_BHBlp_M135135_M0.h5",
        #  "outfpath":os.getcwd()+'/'+"ejecta_id_BHBlp_M135135_M0.h5",
        #  "new_main_pars":{"n_ism":1e-1,"d_l":100e6 * cgs.pc,"z":0.001},
        #  "new_kn_pars":{"p":2.5,"eps_e":1.e-1,"eps_b":1e-1,"theta_obs":90.0},
        #  "fname_light_curve":"lc_BHBlp_M135135_M0.h5",
        #  "fname_ejecta_id":"ejecta_id_BHBlp_M135135_M0.h5",
        #  "fpath_kenta_res":os.getcwd()+'/'+"BHBlp_n01_100Mpc.dat",
        #  "color":"green","label":"BHBlp M135135 M0"
        #  }
        # {"hist_fpath":"hist_vel_inf_SFHo_M135135_M0.dat",
        #  "outfpath":os.getcwd()+'/'+"ejecta_id_SFHo_M135135_M0.h5",
        #  "new_main_pars":{"n_ism":5e-3,"d_l":41.6e6 * cgs.pc,"z":0.0099},
        #  "new_kn_pars":{"p":2.16,"eps_e":1.e-1,"eps_b":1e-2,"theta_obs":0.0},
        #  "fname_light_curve":"lc_SFHo_M135135_M0.h5",
        #  "fname_ejecta_id":"ejecta_id_SFHo_M135135_M0.h5",
        #  "fpath_kenta_res":os.getcwd()+'/'+"SFHo_n0005.dat",
        #  "color":"red","label":"SFHo M135135 M0"
        #  }
    ]

    for task in tasks :

        shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")

        kenta_t, kenta_f = np.loadtxt(task["fpath_kenta_res"], unpack=True, usecols=(0, 2))
        ax.plot(kenta_t[1:-1], kenta_f[1:-1], **{"color": task["color"], "ls": ":", "lw": 0.8, "label": task["label"]})

        # prepare_kn_ej_id_1d(nlayers=30,hist_fpath=task["hist_fpath"],outfpath=task["outfpath"])
        prepare_kn_ej_id_2d(nlayers=30,corr_fpath=task["corr_fpath"],outfpath=task["outfpath"], dist="pw")
        # prepare_kn_ej_id_1d(nlayers=30,hist_fpath=task["hist_fpath"],outfpath=task["outfpath"])


        modify_parfile_par_opt(workingdir=workdir, part="main",
                               newpars=task["new_main_pars"],newopts={},
                               parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        modify_parfile_par_opt(workingdir=workdir, part="kn", newpars=task["new_kn_pars"],
                                   newopts={"method_eos":"Nava13","method_GammaSh":"useJKwithGammaRel",
                                            "method_Delta":"useJoh06","method_Rsh":"useGammaSh",
                                            "method_dmdr":"usingdthdr","method_Up":"useGamma",
                                            "method_dgdr":"peer","method_shock_vel":"sameAsBW",
                                            "method_synchrotron":"WSPN99","method_lf_min":"useU_e",
                                            "method_ne":"useNe", # usenprime
                                            "method_nonreldist":"none","emissivity":"em","absorption":"abs",
                                            "fname_light_curve":task["fname_light_curve"],
                                            "fname_ejecta_id":task["fname_ejecta_id"]},
                               parfile="parfile.par",newparfile="parfile.par",keep_old=False
                              )
        pba.reload_parfile()
        pba.run()
        ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
                **{"color": task["color"], "ls": "--", "lw": 0.8, "label": task["label"]})
        pba.clear()

        shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")

        modify_parfile_par_opt(workingdir=workdir, part="main",newpars=task["new_main_pars"],newopts={},
                               parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        modify_parfile_par_opt(workingdir=workdir, part="kn", newpars=task["new_kn_pars"],
                                   newopts={"method_GammaSh":"useJKwithGammaRel","method_Delta":"useJoh06",
                                            "method_Rsh":"useGammaSh","method_dmdr":"usingdthdr",
                                            "method_Up":"useEint2","method_dgdr":"our",
                                            "method_shock_vel":"shockVel","method_synchrotron":"Joh06",
                                            "method_lf_min":"useU_e","emissivity":"em","absorption":"abs",
                                            "method_ne":"usenprime",
                                            "fname_light_curve":task["fname_light_curve"],
                                            "fname_ejecta_id":task["fname_ejecta_id"]},
                               parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        pba.reload_parfile()
        pba.run()
        ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
                **{"color": task["color"], "ls": "-", "lw": 0.8, "label": task["label"]})
        pba.clear()


    pass

    # ---------------------------------------------
    l11, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='-')
    l12, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--')
    l13, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls=':')

    legend1 = plt.legend([l11, l12, l13],
                         # [r"\& Our Eq. \& J\'{o}hannesson+06", r"\& Peer+12 \& WSPN+99", r"Hotokezaka+15 code"],
                         [r"\texttt{PBA} with J06", r"\texttt{PBA} with P12 \& WSPN99", r"H15"],
                         loc="upper right", fancybox=False, shadow=False, ncol=1, framealpha=0, fontsize=10,
                         borderaxespad=0., frameon=False)
    # legend2 = plt.legend(lls, lbls,
    #                      loc="lower center", fancybox=False, shadow=False, ncol=1, framealpha=0, fontsize=10,
    #                      borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    # ax.add_artist(legend2)
    # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")

    # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=11,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
    ax.set_xlim(1e1, 1e4)
    ax.set_ylim(1e-1, 1e3)
    ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(figdir + save, dpi=256)
    # if save_figs: plt.savefig(PAPERPATH + "figname" + ".pdf")
    # if save_figs: plt.savefig(FIGPATH + figname + ".png", dpi=256)

    plt.show()
    #
    #
    #
    # # ------------------- | DD2 M135135 M0 | ---------------------
    #
    # prepare_kn_ej_id_1d(nlayers=30,
    #                     hist_fpath="hist_vel_inf_DD2_M135135_M0.dat",
    #                     outfpath=os.getcwd()+'/'+"ejecta_id_DD2_M135135_M0.h5")
    # pba.modify_main_part_parfile(newpars={"n_ism":1e-1,"d_l":100e6 * cgs.pc,"z":0.001},newopts={})
    # pba.modify_kn_part_parfile(
    #     newpars={"p":2.5,"eps_e":1.e-1,"eps_b":1e-1,"theta_obs":0.0},
    #     newopts={"method_Up":"useGamma","method_dgdr":"peer","method_shock_vel":"sameAsBW",
    #              "method_synchrotron":"WSPN99","method_lf_min":"useU_e",
    #              "method_nonreldist":"none","emissivity":"em_pl",
    #              "fname_ejecta_id":"ejecta_id_DD2_M135135_M0.h5",
    #              "fname_light_curve":"lc_DD2_M135135_M0.h5"}
    # )
    # pba.run()
    #
    #
    # ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
    #          **{"color": "blue", "ls": "--", "lw": 0.8, "label": "DD2 M135135 LK"})
    # pba.clear()
    #
    # # ------------------- | BHBlp M135135 M0 | ---------------------
    #
    # prepare_kn_ej_id_1d(nlayers=30,
    #                     hist_fpath="hist_vel_inf_BHBlp_M135135_M0.dat",
    #                     outfpath=os.getcwd()+'/'+"ejecta_id_BHBlp_M135135_M0.h5")
    # pba.modify_main_part_parfile(newpars={"n_ism":1e-1,"d_l":100e6 * cgs.pc,"z":0.001},newopts={})
    # pba.modify_kn_part_parfile(
    #     newpars={"p":2.5,"eps_e":1.e-1,"eps_b":1e-1,"theta_obs":0.0},
    #     newopts={"method_Up":"useGamma","method_dgdr":"peer","method_shock_vel":"sameAsBW",
    #              "method_synchrotron":"WSPN99","method_lf_min":"useU_e",
    #              "method_nonreldist":"none","emissivity":"em_pl",
    #              "fname_ejecta_id":"ejecta_id_BHBlp_M135135_M0.h5",
    #              "fname_light_curve":"lc_BHBlp_M135135_LK.h5"}
    # )
    # pba.run()
    #
    # kenta_t, kenta_f = np.loadtxt(os.getcwd()+'/'+"BHBlp_n01_100Mpc.dat", unpack=True, usecols=(0, 2))
    # ax.plot(kenta_t[1:-1], kenta_f[1:-1], **{"color": "green", "ls": ":", "lw": 0.8, "label": "BHBlp M135135 LK"})
    #
    # ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
    #         **{"color": "green", "ls": "--", "lw": 0.8, "label": "BHBlp M135135 LK"})
    # pba.clear()
    #
    # # ------------------- | SFHo M135135 M0 | ---------------------
    #
    # kenta_t, kenta_f = np.loadtxt(os.getcwd()+'/'+"SFHo_n0005.dat", unpack=True, usecols=(0, 2))
    # ax.plot(kenta_t[1:-1], kenta_f[1:-1], **{"color": "red", "ls": ":", "lw": 0.8, "label": "SFHo M135135 M0"})
    #
    # # prepare_kn_ej_id_2d(nlayers=30,
    # #                     corr_fpath="corr_vel_inf_theta_SFHo_M135135_M0.h5",
    # #                     outfpath=os.getcwd()+'/'+"ejecta_id_SFHo_M135135_M0.h5")
    # prepare_kn_ej_id_1d(nlayers=30,
    #                     hist_fpath="hist_vel_inf_SFHo_M135135_M0.dat",
    #                     outfpath=os.getcwd()+'/'+"ejecta_id_SFHo_M135135_M0.h5")
    # pba.modify_main_part_parfile(newpars={"n_ism":5.e-3, "d_l":41.6e6 * cgs.pc,"z":0.0099},newopts={})
    # pba.modify_kn_part_parfile(
    #     newpars={"p":2.16,"eps_e":1.e-1,"eps_b":1.e-2,"theta_obs":0.0},
    #     newopts={"method_Up":"useGamma","method_dgdr":"peer","method_shock_vel":"sameAsBW",
    #              "method_synchrotron":"WSPN99","method_lf_min":"useU_e",
    #              "method_nonreldist":"none","emissivity":"em_pl",
    #              "fname_ejecta_id":"ejecta_id_SFHo_M135135_M0.h5",
    #              "fname_light_curve":"lc_SFHo_M135135_M0.h5"}
    # )
    # pba.run()
    #
    # ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
    #     **{"color": "red", "ls": "--", "lw": 0.8, "label": "SFHo M135135 M0"})
    # pba.clear()
    # # ---------------------------------------------
    # pba.modify_main_part_parfile(newpars={"n_ism":5.e-3, "d_l":41.6e6 * cgs.pc,"z":0.0099},newopts={})
    # pba.modify_kn_part_parfile(
    #     newpars={"p":2.16,"eps_e":1.e-1,"eps_b":1.e-2,"theta_obs":0.0},
    #     newopts={"method_Up":"useEint2","method_dgdr":"our",
    #              "method_shock_vel":"shockVel","method_synchrotron":"Joh06",
    #              "method_lf_min":"useU_e",
    #              "method_nonreldist":"none","emissivity":"em_pl",
    #              "fname_ejecta_id":"ejecta_id_SFHo_M135135_M0.h5",
    #              "fname_light_curve":"lc_SFHo_M135135_M0.h5",
    #              "method_GammaSh":"useJKwithGammaRel",
    #              "method_Delta":"useJoh06","method_Rsh":"useGammaSh",
    #              "method_dmdr":"usingdthdr"}
    # )
    # pba.run()
    #
    # ax.plot(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9) * 1e3,
    #         **{"color": "red", "ls": "-", "lw": 0.8, "label": "SFHo M135135 M0"})
    # pba.clear()
    #
    # # ---------------------------------------------
    # l11, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='-')
    # l12, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--')
    # l13, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls=':')
    #
    # legend1 = plt.legend([l11, l12, l13],
    #                      # [r"\& Our Eq. \& J\'{o}hannesson+06", r"\& Peer+12 \& WSPN+99", r"Hotokezaka+15 code"],
    #                      [r"\texttt{PBA} with J06", r"\texttt{PBA} with P12 \& WSPN99", r"H15"],
    #                      loc="upper right", fancybox=False, shadow=False, ncol=1, framealpha=0, fontsize=10,
    #                      borderaxespad=0., frameon=False)
    # # legend2 = plt.legend(lls, lbls,
    # #                      loc="lower center", fancybox=False, shadow=False, ncol=1, framealpha=0, fontsize=10,
    # #                      borderaxespad=0., frameon=False)
    #
    # ax.add_artist(legend1)
    # # ax.add_artist(legend2)
    # # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")
    #
    # # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=11,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.minorticks_on()
    # ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    # ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
    # ax.set_xlim(1e1, 1e4)
    # ax.set_ylim(1e-1, 1e3)
    # # ax.legend()
    # plt.tight_layout()
    # # if save_figs: plt.savefig(figdir + save, dpi=256)
    # # if save_figs: plt.savefig(PAPERPATH + "figname" + ".pdf")
    # # if save_figs: plt.savefig(FIGPATH + figname + ".png", dpi=256)
    #
    # plt.show()



if __name__ == '__main__':
    main()