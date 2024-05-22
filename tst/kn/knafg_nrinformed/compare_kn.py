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

# try:
#     import package.src.PyBlastAfterglowMag as PBA
# except ImportError:
#     try:
#         import PyBlastAfterglowMag as PBA
#     except:
#         raise ImportError("Cannot import PyBlastAfterglowMag")
import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs,d2d

def main():
    # prepare ejecta ID from NR output
    workdir = os.getcwd()+'/'
    shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes

    tasks = [
        {"corr_fpath":"corr_vel_inf_theta_DD2_M135135_M0.h5",
         "hist_fpath":"hist_vel_inf_DD2_M135135_M0.dat",
         "outfpath":os.getcwd()+'/'+"ejecta_id_DD2_M135135_M0.h5",
         "new_pars":{"n_ism":1e-1,"d_l":100e6 * PBA.utils.cgs.pc,"z":0.001,"theta_obs":0.},
         "new_kn_pars":{"p_fs":2.5,"eps_e_fs":1.e-1,"eps_b_fs":1e-1},
         "fname_light_curve":"lc_DD2_M135135_M0.h5",
         "fname_ejecta_id":"ejecta_id_DD2_M135135_M0.h5",
         "fpath_kenta_res":os.getcwd()+'/'+"DD2_n01_100Mpc.dat",
         "color":"blue","label":"DD2 M135135 M0",
         "workingdir":os.getcwd()+'/'+'workingdir_DD2_M135135_M0/'
         },
        {"corr_fpath":"corr_vel_inf_theta_DD2_M135135_M0.h5",
         "hist_fpath":"hist_vel_inf_DD2_M135135_M0.dat",
         "outfpath":os.getcwd()+'/'+"ejecta_id_DD2_M135135_M0.h5",
         "new_pars":{"n_ism":0.00031,"d_l":1.27e+26,"z":0.0099,"theta_obs":0.785398},
         "new_kn_pars":{"p_fs":2.05,"eps_e_fs":0.1,"eps_b_fs":0.01},
         "fname_light_curve":"lc_DD2_M135135_M0.h5",
         "fname_ejecta_id":"ejecta_id_DD2_M135135_M0.h5",
         "fpath_kenta_res":os.getcwd()+'/'+"DD2_n01_100Mpc.dat",
         "color":"blue","label":"DD2 M135135 M0",
         "workingdir":os.getcwd()+'/'+'workingdir_DD2_M135135_M0/'
         },
        # {"hist_fpath":"corr_vel_inf_theta_BHBlp_M135135_M0.h5",
        #  "outfpath":os.getcwd()+'/'+"ejecta_id_BHBlp_M135135_M0.h5",
        #  "new_main_pars":{"n_ism":1e-1,"d_l":100e6 * cgs.pc,"z":0.001,"theta_obs":90.0},
        #  "new_kn_pars":{"p_fs":2.5,"eps_e_fs":1.e-1,"eps_b_fs":1e-1},
        #  "fname_light_curve":"lc_BHBlp_M135135_M0.h5",
        #  "fname_ejecta_id":"ejecta_id_BHBlp_M135135_M0.h5",
        #  "fpath_kenta_res":os.getcwd()+'/'+"BHBlp_n01_100Mpc.dat",
        #  "color":"green","label":"BHBlp M135135 M0",
        #  "workingdir":os.getcwd()+'/'+'workingdir_BHBlp_M135135_M0/'
        #  },
        # {"hist_fpath":"hist_vel_inf_SFHo_M135135_M0.dat",
        #  "outfpath":os.getcwd()+'/'+"ejecta_id_SFHo_M135135_M0.h5",
        #  "new_main_pars":{"n_ism":5e-3,"d_l":41.6e6 * cgs.pc,"z":0.0099,"theta_obs":0.0},
        #  "new_kn_pars":{"p_fs":2.16,"eps_e_fs":1.e-1,"eps_b_fs":1e-2},
        #  "fname_light_curve":"lc_SFHo_M135135_M0.h5",
        #  "fname_ejecta_id":"ejecta_id_SFHo_M135135_M0.h5",
        #  "fpath_kenta_res":os.getcwd()+'/'+"SFHo_n0005.dat",
        #  "color":"red","label":"SFHo M135135 M0",
        #  "workingdir":os.getcwd()+'/'+'workingdir_SFHo_M135135_M0/'
        #  }
    ]

    for task in tasks :

        kenta_t, kenta_f = np.loadtxt(task["fpath_kenta_res"], unpack=True, usecols=(0, 2))
        ax.plot(kenta_t[1:-1], kenta_f[1:-1], **{"color": task["color"], "ls": ":", "lw": 0.8, "label": task["label"]})

        pba = PBA.wrappers.run_kn(
            working_dir=task["workingdir"],
            struct=dict(struct="numeric",n_layers_pw=None,corr_fpath_david=task["corr_fpath"]),
            P=dict(main=task["new_pars"],
                   kn=task["new_kn_pars"] | dict(method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                                                 use_1d_id="no",do_skymap="no", do_lc="yes",
                                                 ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                                                 method_ne_fs="usenprime",method_comp_mode="observFlux")),
            run=True
        )

        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
                **{"color": task["color"], "ls": "--", "lw": 0.8, "label": task["label"]})
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
    plt.savefig(os.getcwd()+"/figure.png",dpi=256)
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

def main_skymap():
    pba = PBA.wrappers.run_kn(
        working_dir=os.getcwd()+'/tmp1/',
        struct=dict(struct="numeric",n_layers_pw=None,
                    corr_fpath_david=os.getcwd()+"/"+"corr_vel_inf_theta_DD2_M135135_M0.h5"),
        P=dict(main=dict(n_ism=0.00031,d_l=1.27e+26,z=0.0099,theta_obs=np.pi/2.,#0.785398,
                         skymap_freqs='array 1.e9 2.4e9 5.e9',
                         skymap_times='array logspace 3e4 3e8 32'),
               kn=dict(p_fs=2.05,eps_e_fs=0.1,eps_b_fs=0.01)
                  |
                  dict(method_synchrotron_fs="Joh06",method_ele_fs="analytic",use_1d_id="no",do_skymap="yes", do_lc="no",
                       ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                       method_ne_fs="usenprime",method_comp_mode="observFlux",skymap_remove_mu="yes",
                       skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
                                        intp_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect'),  # "gaussian"
                                        hist_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect')))),
        run=True
    )
    config = {
        "gridspec": {
            "width_ratios": (4, 2), "height_ratios": (2, 4),
            "left": 0.14, "right": 0.95, "bottom": 0.1, "top": 0.96, "wspace": 0.05, "hspace": 0.05
        },
        "figname": "out", "paperpath": os.getcwd()+'/', "figfpath": os.getcwd()+'/', "save_pdf": False, "save_figs": True,
        "show_figs": True,
        "grid": False,
        "figsize": (4.8, 4.8),
        "type":"hist",
        "cm": {"color": 'yellow', "marker": "o"},
        "ysize": {"capsize": 2, "color": "yellow", "lw": 0.5},
        "xsize": {"capsize": 2, "color": "yellow", "lw": 0.5},
        "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                       "norm": ("linear", "0.1max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan},
        "xlim": (-1., 1.0), "ylim": (-1.0, 1.0),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "black",
        "histy_backgound_color": "black",
        "plot_grids":True,
        "histx_lim":(1e-4,1e-2),
        "histy_lim":(1e-4,1e-2)
    }
    time= 120*cgs.day; freq= 1e9
    PBA.skymap_plotting_tools.full_plot_skymap_with_hists(
        skymap=pba.KN.get_skymap(time=time, freq=freq),
        conf=config
    )

def main_OLD():
    # prepare ejecta ID from NR output
    workdir = os.getcwd()+'/'
    shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes
    pba = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True)

    tasks = [
        {"corr_fpath":"corr_vel_inf_theta_DD2_M135135_M0.h5",
         "hist_fpath":"hist_vel_inf_DD2_M135135_M0.dat",
         "outfpath":os.getcwd()+'/'+"ejecta_id_DD2_M135135_M0.h5",
         "new_main_pars":{"n_ism":1e-1,"d_l":100e6 * PBA.utils.cgs.pc,"z":0.001},
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
        PBA.id_david.prepare_kn_ej_id_2d(nlayers=None,corr_fpath=task["corr_fpath"],outfpath=task["outfpath"], dist="pw")
        # prepare_kn_ej_id_1d(nlayers=30,hist_fpath=task["hist_fpath"],outfpath=task["outfpath"])


        PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="main",
                                                 newpars=task["new_main_pars"],newopts={},
                                                 parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="kn", newpars=task["new_kn_pars"],
                                                 newopts={"method_eos":"Nava13","method_GammaSh":"useGammaShock",#"useJKwithGammaRel",
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
        pba.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out", loglevel="info")
        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
                **{"color": task["color"], "ls": "--", "lw": 0.8, "label": task["label"]})
        pba.clear()

        shutil.copyfile(os.getcwd()+"/default_parfile.par",os.getcwd()+"/parfile.par")

        PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="main",newpars=task["new_main_pars"],newopts={},
                                                 parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="kn", newpars=task["new_kn_pars"],
                                                 newopts={"method_GammaSh":"useGammaShock",#"useJKwithGammaRel",
                                                          "method_Delta":"useJoh06",
                                                          "method_Rsh":"useGammaSh","method_dmdr":"usingdthdr",
                                                          "method_Up":"useEint2","method_dgdr":"our",
                                                          "method_shock_vel":"shockVel","method_synchrotron":"Joh06",
                                                          "method_lf_min":"useU_e","emissivity":"em","absorption":"abs",
                                                          "method_ne":"usenprime",
                                                          "fname_light_curve":task["fname_light_curve"],
                                                          "fname_ejecta_id":task["fname_ejecta_id"]},
                                                 parfile="parfile.par",newparfile="parfile.par",keep_old=False)
        pba.reload_parfile()
        pba.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out", loglevel="info")
        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc_totalflux(freq=3.e9) * 1e3,
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
    plt.savefig(os.getcwd()+"/figure.png",dpi=256)
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
    main_skymap()
    # main()