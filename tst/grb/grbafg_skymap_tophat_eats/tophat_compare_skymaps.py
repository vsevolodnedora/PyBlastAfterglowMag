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

# try:
#     from PyBlastAfterglowMag.interface import modify_parfile_par_opt
#     from PyBlastAfterglowMag.interface import PyBlastAfterglow
#     from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#     from PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
#                                            BetFromMom, GamFromMom, MomFromGam)
#     from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     from PyBlastAfterglowMag.skymap_tools import \
#         (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
#                                                            BetFromMom, GamFromMom, MomFromGam)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#         from package.src.PyBlastAfterglowMag.skymap_tools import \
#             (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt, combine_images_old
from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
                                                   BetFromMom, GamFromMom, MomFromGam)
from package.src.PyBlastAfterglowMag.id_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
from package.src.PyBlastAfterglowMag.skymap_plotting_tools import \
    (plot_skymaps, plot_skymap_properties_evolution, plot_one_skymap_with_dists, precompute_skymaps,
     combine_images, get_skymap_lat_dist, plot_skymap_with_hists, get_skymap_fwhm)




def compare_skymaps(resolutions=((80,100,120),
                                 # (21,41,81,101,121),
                                 # (21,81,121,161,201),
                                 (381,401,521),
                                 ('red','orange','yellow', 'cyan', 'lime'))):

    fig,axes = plt.subplots(ncols=2,nrows=len(resolutions[0])+1,sharex='all',sharey='row',figsize=(5,10))
    times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
    time_ = 1600
    freq = 1e9
    tmp={"hist_nx": 81, "hist_ny": 81, "spec": False,
         "smooth": {},  # {"type": "gaussian", "sigma": 10},
         "cm": {"color": 'cyan', "marker": "+"},
         "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
         "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
         "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                        "norm": ("linear", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": 0}
         }
    settings={
        "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "black",
        "histy_backgound_color": "black",
        "plot_grids": True,
        "histx_lim":(1e-2, 1e1),
        "histy_lim":(1e-2, 1e1)
    }

    theta_w = 0.05# np.pi/2.
    nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
    nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0
    nres = np.max([nres_pw,nres_a])
    tot_fluxes = {}
    for i in range(nres):
        nlayer_a = resolutions[1][i] if nres_a > 0 else 0
        nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0
        color=resolutions[2][i] if nres_pw > 0 else resolutions[2]
        tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color

        # ------ Piece Wise -------
        if (nres_pw > 0):
            prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                                  "nlayers_pw": nlayer_pw, "nlayers_a": 1, "struct":"tophat"},type='pw',outfpath="tophat_grb_id_pw.h5")

            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                                   newpars={},
                                   newopts={"fname_ejecta_id":"tophat_grb_id_pw.h5","method_eats":"piece-wise",
                                            "fname_sky_map":"skymap_tophat_pw.h5", "fname_light_curve":"lc_tophat_pw.h5"},
                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
            pba_pw = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
            pba_pw.run(loglevel='info')

            # all_x_jet, all_y_jet, all_fluxes_jet \
            #     = pba_pw.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True, remove_zeros=True)
            all_x_jet, all_y_jet, all_fluxes_jet \
                = pba_pw.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True)
            # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
            #                                             pre_int=True, k=5, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=1.1)
            int_x_j, int_y_j, int_zz_j = combine_images_old(all_x_jet, all_y_jet, all_fluxes_jet, shells=True, verbose=True,
                                                           hist_or_int="hist", nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                   collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                   collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_m_j, yc_m_j = pba_pw.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
            if (nres_a>0 and nres_pw>0):
                ax_main = axes[i+1][0]; ax_histx=axes[0,0]; ax_histy=None
            else:
                ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
            im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
                                        grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)
            # ax_cbar = axes[i+1][0]
            # divider = make_axes_locatable(ax_cbar)
            # cax = divider.append_axes('left', size='99%', pad=0.9)
            # plt.delaxes(ax_cbar)
            # cbar = plt.colorbar(im, cax=cax,
            #                     # format='%.0e',ticks=ticks
            #                     orientation='vertical',
            #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
            #                     )

        # ------ Adaptive -------
        if (nres_a > 0):
            prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                                  "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                                   newpars={"nsublayers":nlayer_a},
                                   newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive",
                                            "fname_sky_map":"skymap_tophat_a.h5", "fname_light_curve":"lc_tophat_a.h5"},
                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
            pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
            pba_a.run(loglevel='info')

            all_x_jet, all_y_jet, all_fluxes_jet \
                = pba_a.GRB.get_skymap(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=True, renormalize=False, normtype='pw', remove_zeros=True)
            int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet, verbose=True,
                                                        pre_int=True, k=2, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                              collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                              collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
            if (nres_a>0 and nres_pw>0):
                ax_main = axes[i+1][1]; ax_histx=axes[0,1]; ax_histy=None
            else:
                ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
            im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
                                        grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)

        tot_fluxes[f"nlayer_a={nlayer_a} nlayer_pw={nlayer_pw}"] = f" PW={pba_pw.GRB.get_skymap_totfluxes(freq=freq,time=time_* cgs.day)}" + \
                                                                       f" lc:({pba_pw.GRB.get_lc(freq=freq,time=time_* cgs.day)})" + \
                                                                       f" A={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=time_* cgs.day)}" + \
                                                                       f" lc:({pba_a.GRB.get_lc(freq=freq,time=time_* cgs.day)})"

        # --- COLOBAR
            # ax_cbar = axes[1,i+1]
            # divider = make_axes_locatable(ax_cbar)
            # cax = divider.append_axes('right', size='99%', pad=0.9)
            # plt.delaxes(ax_cbar)
            # cbar = plt.colorbar(im, cax=cax,
            #                     # format='%.0e',ticks=ticks
            #                     orientation='vertical',
            #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
            #                     )
    for key, val in tot_fluxes.items():
        print(f"{key} {val}")
    # pass

    if hasattr(axes[0],'__len__'):
        axes_ = axes[0,:]
        axes[0,0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)  # , color='blue')
    else:
        axes_ = axes
        axes[0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
    for i_, ax_histx in enumerate(axes_):
        ax_histx.set_yscale("linear")
        ax_histx.set_xscale("linear")
        ax_histx.minorticks_on()
        ax_histx.set_yscale("log")
        # # ax_histx.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
        if (i_ == 0):
            ax_histx.tick_params(axis='both', which='both', labelleft=True,
                                labelright=False, tick1On=True, tick2On=True,
                                labelsize=12, color="white",
                                direction='in',
                                bottom=True, top=True, left=True, right=True)
        else:
            ax_histx.tick_params(axis='both', which='both', labelleft=False,
                                labelright=False, tick1On=True, tick2On=True,
                                labelsize=12, color="white",
                                direction='in',
                                bottom=True, top=True, left=True, right=True)
        ax_histx.set_facecolor(settings["histx_backgound_color"])

        if settings["plot_grids"]:
            ax_histx.grid()
        # ax_histx.set_ylim(*settings["histx_lim"])
    # ax_histy.set_yscale("linear")
    # ax_histy.set_xscale("linear")
    # ax_histy.minorticks_on()
    # ax_histy.set_xscale("log")
    # ax_histy.set_xlim(ax_histx.get_ylim())
    # # ax_histy.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    # ax_histy.set_xlabel(r"$I_{\nu;\rm m}(z)$", fontsize=12)  # , color='blue')
    # ax_histy.tick_params(axis='both', which='both', labelleft=False,
    #                      labelright=False, tick1On=True, tick2On=True,
    #                      labelsize=12, color="white",
    #                      direction='in',
    #                      bottom=True, top=True, left=True, right=False)
    # ax_histy.set_facecolor(settings["histy_backgound_color"])


    if hasattr(axes[0],'__len__'):
        for i_, axes_main in enumerate(zip(axes[1:,0],axes[1:,1])):
            for j_, ax_main in enumerate(axes_main):
                ax_main.set_yscale("linear")
                ax_main.set_xscale("linear")
                ax_main.minorticks_on()
                ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
                ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
                if (j_ == 0):
                    ax_main.tick_params(axis='both', which='both', labelleft=True,
                                        labelright=False, tick1On=True, tick2On=True,
                                        labelsize=12, color="white",
                                        direction='in',
                                        bottom=True, top=True, left=True, right=True)
                    ax_main.set_ylabel("$Z$ [mas]")
                else:
                    ax_main.tick_params(axis='both', which='both', labelleft=False,
                                        labelright=False, tick1On=True, tick2On=True,
                                        labelsize=12, color="white",
                                        direction='in',
                                        bottom=True, top=True, left=True, right=True)


                if settings["plot_grids"]:
                    ax_main.grid()

    plt.tight_layout()
    plt.show()
def compare_skymaps_theta_im_max(nlayer_a = 321, theta_maxs=((0.2, 0.4, .9, 1.2,  1.57),
                                 ('red','orange','yellow', 'cyan', 'lime'))):

    fig,axes = plt.subplots(ncols=1,nrows=len(theta_maxs[0])+1,sharex='all',figsize=(5,10))
    times = [1.,10.,40.,100.,200.]
    freq = 1e9
    tmp={"hist_nx": 71, "hist_ny": 71, "spec": False,
         "smooth": {},  # {"type": "gaussian", "sigma": 10},
         "cm": {"color": 'cyan', "marker": "+"},
         "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
         "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
         "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                        "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
         }
    settings={
        "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "black",
        "histy_backgound_color": "black",
        "plot_grids": True,
        "histx_lim":(1e-2, 1e1),
        "histy_lim":(1e-2, 1e1)
    }

    theta_w = 0.2# np.pi/2.
    tot_fluxes = {}
    for i in range(len(theta_maxs[0])):
        theta_max = theta_maxs[0][i]
        # nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0


        nlayer_pw = 80

        color=theta_maxs[1][i]
        tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color

        prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                              "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars={"nsublayers":nlayer_a, "im_max_theta":theta_max},
                               newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive"},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        pba_a.run(loglevel='info')

        all_x_jet, all_y_jet, all_fluxes_jet \
            = pba_a.GRB.get_skymap(time=times[3] * cgs.day, freq=freq, verbose=False, remove_mu=True, renormalize=True, normtype='pw')
        int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
                                                    hist_or_int="hist", shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
        grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                               collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
        grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                               collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
        xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
        ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
        im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
                                    grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)

        tot_fluxes[f"nlayer_a={nlayer_a} theta_max={theta_max}"] = f" A={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=times[3]* cgs.day)}"
        # --- COLOBAR
        # ax_cbar = axes[1,i+1]
        # divider = make_axes_locatable(ax_cbar)
        # cax = divider.append_axes('right', size='99%', pad=0.9)
        # plt.delaxes(ax_cbar)
        # cbar = plt.colorbar(im, cax=cax,
        #                     # format='%.0e',ticks=ticks
        #                     orientation='vertical',
        #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
        #                     )
    for key, val in tot_fluxes.items():
        print(f"{key} {val}")
    # pass

    axes[0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
    axes[0].set_yscale("linear")
    axes[0].set_xscale("linear")
    axes[0].minorticks_on()
    axes[0].set_yscale("log")
    axes[0].tick_params(axis='both', which='both', labelleft=True,
                         labelright=False, tick1On=True, tick2On=True,
                         labelsize=12, color="white",
                         direction='in',
                         bottom=True, top=True, left=True, right=True)
    axes[0].set_facecolor(settings["histx_backgound_color"])
    for j_, ax_main in enumerate(axes[1:]):
        ax_main.set_yscale("linear")
        ax_main.set_xscale("linear")
        ax_main.minorticks_on()
        ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
        ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
        ax_main.tick_params(axis='both', which='both', labelleft=True,
                            labelright=False, tick1On=True, tick2On=True,
                            labelsize=12, color="white",
                            direction='in',
                            bottom=True, top=True, left=True, right=True)
        ax_main.set_ylabel("$Z$ [mas]")


        if settings["plot_grids"]:
            ax_main.grid()

    plt.tight_layout()
    plt.show()
def compare_skymaps_3d(resolutions=((20,40,80,100,120),
                                 # (21,41,81,101,121),
                                 (21,81,121,161,201),
                                 # (121,241,381,401,521),
                                 ('red','orange','yellow', 'cyan', 'lime')),
                                log=True):

    theta_w = 0.2# np.pi/2.
    times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
    time = 1600
    freq = 1e9

    nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
    nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0

    fig = plt.figure(figsize=plt.figaspect(0.5))

    nres = np.max([nres_pw,nres_a])
    for i in range(nres):
        nlayer_a = resolutions[1][i] if nres_a > 0 else 0
        nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0


        prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                              "nlayers_pw": nlayer_pw, "nlayers_a": 1, "struct":"tophat"},type='pw',outfpath="tophat_grb_id_pw.h5")

        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars={},
                               newopts={"fname_ejecta_id":"tophat_grb_id_pw.h5","method_eats":"piece-wise",
                                        "fname_sky_map":"skymap_tophat_pw.h5", "fname_light_curve":"lc_tophat_pw.h5"},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        pba_pw = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        pba_pw.run(loglevel='info')

        all_x_jet_pw, all_y_jet_pw, all_fluxes_jet_pw \
            = pba_pw.GRB.get_skymap_old(time=time * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True)

        ax_a = fig.add_subplot(1, 2, 1, projection='3d')

        for (x_jet_pw, y_jet_pw, fluxes_jet_pw) in zip(all_x_jet_pw, all_y_jet_pw, all_fluxes_jet_pw):
            # fluxes_jet_pw+=1
            # fluxes_jet_pw[fluxes_jet_pw<1e-4] = 0.
            cmap = cm.get_cmap('inferno')
            if log:
                my_norm = LogNorm(fluxes_jet_pw.max() * 1e-2, fluxes_jet_pw.max())
                c = cmap(my_norm(fluxes_jet_pw.flatten()))
                fluxes_jet_pw = np.log10(fluxes_jet_pw)
            else:
                my_norm = Normalize(fluxes_jet_pw.max() * 1e-2, fluxes_jet_pw.max())
                c = cmap(my_norm(fluxes_jet_pw.flatten()))
                fluxes_jet_pw = fluxes_jet_pw
            ax_a.scatter(x_jet_pw.flatten(), y_jet_pw.flatten(), fluxes_jet_pw.flatten(),  c=c)
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cphis_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
            ax_a.set_xlabel('X Label')
            ax_a.set_ylabel('Y Label')
            ax_a.set_zlabel('I Label')
            # ax_a.set_zlim(np.log10(fluxes_jet_pw.flatten()).max()*1e-2,np.log10(fluxes_jet_pw.flatten()).max())
        ax_a.set_title(f"Tot.Flux={pba_pw.GRB.get_skymap_totfluxes(freq=freq,time=time * cgs.day)}")
        # ---------------------- #

        prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                              "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars={"nsublayers":nlayer_a},
                               newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive",
                                        "fname_sky_map":"skymap_tophat_a.h5", "fname_light_curve":"lc_tophat_a.h5"},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        pba_a.run(loglevel='info')

        all_x_jet_a, all_y_jet_a, all_fluxes_jet_a \
            = pba_a.GRB.get_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True, renormalize=False, normtype='pw', remove_zeros=True)



        ax_b = fig.add_subplot(1, 2, 2, projection='3d')
        xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a)

        my_norm_a = None
        for (x_jet_a, y_jet_a, fluxes_jet_a) in zip(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a):
            # fluxes_jet_a += 1
            # fluxes_jet_a[fluxes_jet_a<1e-4] = 0.
            cmap = cm.get_cmap('inferno')
            if log:
                if (my_norm_a is None): my_norm_a = LogNorm(fluxes_jet_a.max() * 1e-2, fluxes_jet_a.max())
                c = cmap(my_norm_a(fluxes_jet_a.flatten()))
                fluxes_jet_a = np.log10(fluxes_jet_a)
            else:
                if (my_norm_a is None): my_norm_a = Normalize(fluxes_jet_a.max() * 1e-2, fluxes_jet_a.max())
                c = cmap(my_norm_a(fluxes_jet_a.flatten()))
                fluxes_jet_a = fluxes_jet_a
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(rs_i).flatten(), c=cmap(my_norm(int_i.flatten())))
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cthetas_i.flatten(),  c=cmap(my_norm(int_i.flatten())))

            ax_b.scatter(x_jet_a.flatten(), y_jet_a.flatten(), fluxes_jet_a.flatten(),  c=c)
            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cphis_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
            ax_b.set_xlabel('X Label')
            ax_b.set_ylabel('Y Label')
            ax_b.set_zlabel('I Label')
            ax_b.scatter(np.full_like(fluxes_jet_a,xc_m_j), np.full_like(fluxes_jet_a,yc_m_j), fluxes_jet_a.flatten())

        # ax_b.set_zlim(np.log10(fluxes_jet_a.flatten()).max()*1e-2,np.log10(fluxes_jet_a.flatten()).max())
        ax_b.set_title(f"Tot.Flux={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=time * cgs.day)}")


    print(f"Piecewise: I[{fluxes_jet_pw.min(),fluxes_jet_pw.max()}]")
    print(f"Adaptive: I[{fluxes_jet_a.min(),fluxes_jet_a.max()}]")

    # plt.tight_layout()
    plt.show()
def compare_skymaps_3d_theta_im_max(theta_maxs=((0.2, 0.4, .9, 1.2,  1.57),
                                                 (121,321,521,621,721),
                                                 ('red','orange','yellow', 'cyan', 'lime')),
                                    extend=2, nx = 41, ny = 31):

    theta_w = 0.2# np.pi/2.
    times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
    time_ = 1600
    freq = 1e9

    # nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
    # nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0

    fig = plt.figure(figsize=(15,10))#plt.figaspect(0.5),fig)
    res = {}
    # nres = np.max([nres_pw,nres_a])

    my_norm = None
    my_norm2 = None
    for i in range(len(theta_maxs[0])):
        # nlayer_a = resolutions[1][i] if nres_a > 0 else 0
        # nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0

        # ---------------------- #
        nlayer_a = theta_maxs[1][i]
        prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                              "nlayers_pw": 60, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars={"nsublayers":nlayer_a, "im_max_theta":theta_maxs[0][i]},
                               newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive"},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        pba_a.run(loglevel='info')

        # x_jet_a, y_jet_a, fluxes_jet_a \
        #     = pba_a.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=True, remove_mu=False,
        #                                renormalize=False, normtype='a')
        x_jet_a, y_jet_a, fluxes_jet_a \
            = pba_a.GRB.get_skymap(time=time_ * cgs.day, freq=freq, verbose=True, remove_mu=False,
                                   renormalize=False, normtype='a', remove_zeros=True)

        xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(x_jet_a, y_jet_a, fluxes_jet_a)
        print(f"xc={xc_m_j} yc={yc_m_j}")

        # ax_b = fig.add_subplot(len(theta_maxs[0]), 2, i*2+1, projection='3d')
        ax_b = fig.add_subplot(len(theta_maxs[0]), 2, i*2+1)

        # for (x_jet_a, y_jet_a, fluxes_jet_a) in zip(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a):
        # fluxes_jet_a += 1
        # fluxes_jet_a[fluxes_jet_a<1e-4] = 0.
        cmap = cm.get_cmap('viridis')
        # my_norm = Normalize(fluxes_jet_a.max() * 1e-1, fluxes_jet_a.max())
        if my_norm is None: my_norm = Normalize(fluxes_jet_a[0].max() * 1e-1, fluxes_jet_a[0].max())
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
        # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(rs_i).flatten(), c=cmap(my_norm(int_i.flatten())))
        # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
        # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cthetas_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
        # fluxes_jet_=fluxes_jet_a;#np.log10(fluxes_jet_a)
        c = ax_b.scatter(x_jet_a[0].flatten(), y_jet_a[0].flatten(), c=cmap(my_norm(fluxes_jet_a[0].flatten())))
        c = ax_b.scatter(x_jet_a[1].flatten(), y_jet_a[1].flatten(), c=cmap(my_norm(fluxes_jet_a[1].flatten())))
        # cbar = plt.colorbar(ScalarMappable(cmap=c.get_cmap(), norm=c.norm), ax=ax_b)

        int_x_j, int_y_j, int_zz_j = combine_images(x_jet_a, y_jet_a, fluxes_jet_a,
                                                    pre_int=True, k=3, nx=nx, ny=ny, extend=2)
        #
        # int_x_j, int_y_j, int_zz_j = combine_images_old(x_jet_a, y_jet_a, fluxes_jet_a,
        #                                                 shells=True, hist_or_int='hist', nx=nx, ny=ny, extend=2)
        grid_X, grid_Y = np.meshgrid(int_x_j, int_y_j, indexing='ij')
        max = 5.5
        if my_norm2 is None: my_norm2 = LogNorm(int_zz_j.max() * 1e-1, int_zz_j.max())
        # ax_x = fig.add_subplot(len(theta_maxs[0]), 2, (i+1)*2, projection='3d')
        ax_x = fig.add_subplot(len(theta_maxs[0]), 2, (i+1)*2)
        if int_x_j.ndim == 1:
            int_x_j, int_y_j = np.meshgrid(int_x_j, int_y_j, indexing='ij')

        c = ax_x.scatter(int_x_j.flatten(), int_y_j.flatten(), c=cmap(my_norm2(int_zz_j.flatten())), marker='s')
            # ax_x.pcolormesh(grid_X, grid_Y, i_zz,  cmap=cmap,norm=my_norm)
            # cbar = plt.colorbar(ScalarMappable(cmap=c.get_cmap(), norm=c.norm), ax=ax_x)

        # xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(x_jet_a, y_jet_a, fluxes_jet_a)
        ax_b.plot([xc_m_j], [yc_m_j], color="lime", marker='o', ms=10)
        ax_b.axhline(y=0, color='gray')

        ax_b.set_title(f"Tot.Flux={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=time_ * cgs.day)}")
        res[f"thmax={theta_maxs[0][i]}"] = f"xmin={x_jet_a[0].min():.2e} xmax={x_jet_a[0].min():.2e} " \
                                           f"| ymin={y_jet_a[0].min():.2e} ymax={y_jet_a[0].min():.2e} " \
                                            f"| Imin={fluxes_jet_a[0].min():.2e} Imax={fluxes_jet_a[0].max():.2e}" \
                                            f"| IHISTmin={int_zz_j.min():.2e} IHISTmax={int_zz_j.max():.2e}" \

    # print(f"Piecewise: I[{fluxes_jet_pw.min(),fluxes_jet_pw.max()}]")
    print(f"Adaptive: I[{fluxes_jet_a[0].min(),fluxes_jet_a[0].max()}]")

    for key, val in res.items():
        print(f"{key} {val}")

    # plt.tight_layout()
    plt.show()

def compare_skymap_evolution(resolutions=((21,81,121,161,201),
                                          ('red','orange','yellow', 'cyan', 'lime')),
                                          method_spread="AFGPY"):
    times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
    nlayer_pw = 80
    theta_w = 0.2
    freq = 1.e9
    tmp = {"hist_nx":91, "hist_ny":91}

    res = {
        "pw":{"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),"xsize":np.zeros_like(times), "ysize":np.zeros_like(times),"fnu":np.zeros_like(times)}
    }
    for ir, sublayers in enumerate(resolutions[0]):
        res[f"a{sublayers}"] = {"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),"xsize":np.zeros_like(times), "ysize":np.zeros_like(times),"fnu":np.zeros_like(times)}


    prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                          "nlayers_pw": nlayer_pw, "nlayers_a": 1, "struct":"tophat"},type='pw',outfpath="tophat_grb_id_pw.h5")

    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                           parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                           newpars={},
                           newopts={"fname_ejecta_id":"tophat_grb_id_pw.h5","method_eats":"piece-wise", "method_spread":method_spread},
                           parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    pba_pw = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
    pba_pw.run(loglevel='info')

    # Piece wise
    for it, t in enumerate(times):
        # Plot peice-wise

        all_x_jet, all_y_jet, all_fluxes_jet \
            = pba_pw.GRB.get_skymap(time=t * cgs.day, freq=freq, verbose=False, remove_mu=False, renormalize=True)
        # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
        #                                             hist_or_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)

        xc_m_j, yc_m_j = pba_pw.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)

        grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                               collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
        y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
        grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                               collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
        x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)

        res["pw"]["xc"][it] = xc_m_j
        res["pw"]["yc"][it] = yc_m_j
        res["pw"]["xsize"][it] = x2-x1
        res["pw"]["ysize"][it] = y2-y1
        res["pw"]["fnu"][it] = pba_pw.GRB.get_skymap_totfluxes(freq=freq,time=t*cgs.day)


    # adaptive
    for ir, nsublayers in enumerate(resolutions[0]):

        prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
                              "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars={"nsublayers": nsublayers},
                               newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive", "method_spread":method_spread},
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        pba_a.run(loglevel='info')

        for it, t in enumerate(times):

            all_x_jet, all_y_jet, all_fluxes_jet \
                = pba_a.GRB.get_skymap(time=t * cgs.day, freq=freq, verbose=False, remove_mu=False, renormalize=False, normtype='a')
            # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
            #                                             pre_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)

            xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)

            grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                   collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)

            grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                   collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)

            res[f"a{nsublayers}"]["xc"][it] = xc_m_j
            res[f"a{nsublayers}"]["yc"][it] = yc_m_j
            res[f"a{nsublayers}"]["xsize"][it] = x2-x1
            res[f"a{nsublayers}"]["ysize"][it] = y2-y1
            res[f"a{nsublayers}"]["fnu"][it] = pba_a.GRB.get_skymap_totfluxes(freq=freq,time=t*cgs.day)

    # plotting
    fig, axes=plt.subplots(ncols=1,nrows=5,sharex='all',figsize=(10,7))
    axes[0].plot(res["pw"]["t"], res["pw"]["xc"], color='black', ls='-', lw=2, alpha=1, marker='x')
    axes[1].plot(res["pw"]["t"], res["pw"]["yc"], color='black', ls='-', lw=2, alpha=1, marker='x')
    axes[2].plot(res["pw"]["t"], res["pw"]["xsize"], color='black', ls='-', lw=2, alpha=1, marker='x')
    axes[3].plot(res["pw"]["t"], res["pw"]["ysize"], color='black', ls='-', lw=2, alpha=1, marker='x')
    axes[4].plot(res["pw"]["t"], res["pw"]["fnu"], color='black', ls='-', lw=2, alpha=1, marker='x')
    for nsublayers, color in zip(resolutions[0],resolutions[1]):
        axes[0].plot(res[f"a{nsublayers}"]["t"], res[f"a{nsublayers}"]["xc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
        axes[1].plot(res[f"a{nsublayers}"]["t"], res[f"a{nsublayers}"]["yc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
        axes[2].plot(res[f"a{nsublayers}"]["t"], res[f"a{nsublayers}"]["xsize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
        axes[3].plot(res[f"a{nsublayers}"]["t"], res[f"a{nsublayers}"]["ysize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
        axes[4].plot(res[f"a{nsublayers}"]["t"], res[f"a{nsublayers}"]["fnu"], color=color, ls='-', lw=2, alpha=.7, marker='o')

    for ax,lbl in zip(axes,["xc", "yc", "xsize", "ysize"]):
        ax.set_ylabel(lbl)

    axes[0].set_ylim(-2,2)
    axes[1].set_ylim(-2,2)
    axes[4].set_yscale('log')
    axes[len(axes)-1].set_xlabel("time [d]")

    plt.show()

def main():
    # task_kn_skymap_with_dist_one_time(method_eats="piece-wise",type="pw")
    # task_kn_skymap_with_dist_one_time(method_eats="adaptive",type="a")
    compare_skymaps(resolutions=((80,100,120,140),
                                 # (21,41,81,101,121),
                                 # (81,),
                                 # (220,),
                                 (101,141,181,221),
                                 # (121,),
                                 ('red','orange','yellow', 'cyan', 'lime')))
    # compare_skymaps_3d(resolutions=((80,),
    #                              (701,),
    #                              # (21,81,121,161,201),
    #                              # (121,241,381,401,521),
    #                              ('red','orange','yellow', 'cyan', 'lime')))
    # compare_skymaps(resolutions=((),
    #                              # (21,41,81,101,121),
    #                              (121,),
    #                              # (121,241,381,401,521),
    #                              ('orange')))
    compare_skymaps_theta_im_max(nlayer_a = 21, theta_maxs=((0.2, 0.4, .9, 1.2,  1.57),
                                 ('red','orange','yellow', 'cyan', 'lime')))

    # compare_skymaps_3d(resolutions=((21,),
    #                              (31,),
    #                              # (21,81,121,161,201),
    #                              # (121,241,381,401,521),
    #                              ('red','orange','yellow', 'cyan', 'lime')))
    # compare_skymaps_3d_theta_im_max(theta_maxs=((1.57,),
    #                                             (81,),
    #                                             ('red', 'blue', 'orange')))


    compare_skymap_evolution()
if __name__ == '__main__':
    main()