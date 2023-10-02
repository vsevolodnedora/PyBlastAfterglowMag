# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import gridspec
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

# import PyBlastAfterglowMag as PBA
import package.src.PyBlastAfterglowMag as PBA
# plt.style.use('seaborn-v0_8')

# plt.style.use('dark_background')

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
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
#                                                            BetFromMom, GamFromMom, MomFromGam)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#         from package.src.PyBlastAfterglowMag.skymap_plotting_tools import \
#             (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

day = 86400



def task_kn_skymap_with_dist_one_time_OLD():

    workdir = os.getcwd()+'/'

    PBA.id_maker_analytic.prepare_grb_ej_id_1d({"struct":"gaussian",
                          "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
                          "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10}, type="pw",
                          outfpath=workdir+"gauss_grb_id.h5")

    pba = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True,parfile="parfile.par")

    pba.reload_parfile()

    pba.run()

    time = 240
    workdir = os.getcwd() + '/'
    sim = "BLh_M13641364_M0_LK_SR"
    # # res_dir_kn = "/media/vsevolod/Data/postprocessed5/{}/outflow_1/geo//afterglow/tst_spec_maps/".format(sim)
    # res_dir_kn = workdir#root + "BLh_M13641364_M0_LK_SR" + "/outflow_1/" + MASK + "/" + main_dir + res_dir
    # res_dir_grb = workdir#root + "BLh_M13641364_M0_LK_SR" + "/outflow_1/" + MASK + "/" + main_dir + res_dir
    # prefix_pl = "skymap_kn_example"#"ejecta_nlayers100_theta450_nism000031_p205_epse01_epsb001_epst1_dl00413"
    # prefix_th = "skymap_kn_example"#"ejecta_nlayers100_theta450_nism000031_p205_epse01_epsb001_epst1_dl00413"
    # jet_prefix = "skymap_grb_gauss"#"jet_G_nism000031_thetaobs450_epse007079_epsb000525_p216_epst00_nlayers100_dl00413_"
    # # prefix_root = prefix_pl.replace("theta{}_", "")
    task_to_plot = {
        "workingdir" : workdir,
        "parfile1":  "parfile.par",
        "parfile2":  "parfile.par",
        "time": time, "freq": 1e9, "title":"BLh $q=1.00$"
    }

    fname = workdir+ "precomputed_grb_skymap"#"int_"+sim + "_with_hists_time{}_".format(time) + prefix_th
    settings = {
        "gridspec": {
            "width_ratios": (4, 2), "height_ratios": (2, 4),
            "left": 0.14, "right": 0.95, "bottom": 0.1, "top": 0.96, "wspace": 0.05, "hspace": 0.05
        },
        # "figname": fname,
        "paperpath": os.getcwd()+'/', "figfpath": fname, "save_pdf": True, "save_figs": True,
        "show_figs": True,
        "grid": False,
        "figsize": (4.8, 4.8),
        # "workingdir" : workdir,
        "rerun":False, "precompute": True,
        # "precompute_fpath": "/media/vsevolod/Data/postprocessed5/tmp/" + sim + "_single_time{}".format(time) + ".h5",
        "precompute_fpath": fname + ".h5",#res_dir_kn+"int_skymaps_kn_" + prefix_root + ".h5",
        # "precompute_fpath": OUTDIR + fname + ".h5",
        "kn_skymap": {
            # "hist_nx": 1024, "hist_ny": 512, "spec": False,
            # "smooth": {},  # {"type": "gaussian", "sigma": 10},
            # # "smooth": {"type": "gaussian", "sigma": 10},#"smooth": {"type": "uniform", "sigma": 5},
            # "cm": {"color": 'yellow', "marker": "o"},
            # "ysize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            # # "pcolormesh": {"cmap": 'viridis', "norm": ("log", "0.2max", "1max"), "facecolor": 0.0, "alpha": 1.0},
            # "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
            #                "norm": ("log", "0.2max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
            # {}
        },
        "grb_skymap": {
            "hist_nx": 256, "hist_ny": 128, "spec": False,
            "smooth": {},  # {"type": "gaussian", "sigma": 10},
            "cm": {"color": 'cyan', "marker": "+"},
            "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                           "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
        },
        "kn_grb_skymap": {
            # "hist_nx": 1024, "hist_ny": 512, "spec": False,
            # "smooth": {},  # {"type": "gaussian", "sigma": 10},  # "smooth": {"type": "uniform", "sigma": 5},
            # "cm": {"color": 'lime', "marker": "o"},
            # "ysize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # # "pcolormesh": {"cmap": 'viridis', "norm": ("log", "0.2max", "1max"), "facecolor": 0.0, "alpha": 1.0},
            # "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
            #                "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
        },
        "kn_w_skymap": {
            # "hist_nx": 125, "hist_ny": 75,
            # "cm": {"color": 'red', "marker": "x"},
            # "ysize": {"capsize": 2, "color": "orange", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "orange", "lw": 0.5},
            # "pcolormesh":{"cmap":'viridis', "norm":("log","0.2max","1max"), "facecolor":0.0, "alpha":1.0}
        },
        "kn_skymap_ratio": {
            # "pcolormesh": {"cmap": 'viridis', "norm": ("twoslope", 0, 1, 4), "facecolor": 0.0, "alpha": 1.0}  # {}
        },
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
    PBA.skymap_plotting_tools.plot_one_skymap_with_dists(task_to_plot=task_to_plot, settings=settings)

def plot_stacked_heatmaps(heatmaps, transparency=0.5):
    """
    Plot a set of heatmaps as stacked 2D plots in 3D.

    :param heatmaps: List of 2D numpy arrays representing heatmaps.
    :param transparency: Transparency value of each heatmap (0: fully transparent, 1: fully opaque).
    :return: None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    num_maps = len(heatmaps)
    cmap = plt.cm.viridis
    min_val = min(hm.min() for hm in heatmaps)
    max_val = max(hm.max() for hm in heatmaps)
    norm = plt.Normalize(min_val, max_val)

    for idx, heatmap in enumerate(heatmaps):
        x = np.arange(heatmap.shape[1])
        y = np.arange(heatmap.shape[0])
        x, y = np.meshgrid(x, y)
        z = np.ones_like(x) * idx

        color = cmap(norm(heatmap))
        color[..., -1] = transparency  # set alpha (transparency)

        surf = ax.plot_surface(x, y, z, facecolors=color, rstride=1, cstride=1, shade=False)

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    fig.colorbar(mappable, ax=ax)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Layer')

    plt.show()

def plot_skymaps_3d(skymaps : list[PBA.interface.Skymap]):

    # g_min = np.min([np.min(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_gmax = np.max([np.max(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_min = g_gmax * 1e-3
    print(f"g_min={g_min} g_gmax={g_gmax}")

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
                     gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)

    # ax = fig.add_subplot(211,projection='3d')
    #
    #
    # for skymap in skymaps:
    #     ax.plot(skymap.grid_x, np.log10(skymap.dist_x), zs=np.log10(skymap.time), zdir='x', color="white")
    # # ax.set_zscale('log')
    #
    #
    # ax.set_ylabel(r"$X$ [mas]")
    # ax.set_zlabel(r"$Z$ [mas]")
    # ax.set_xlabel(r"time [s]")
    # ax.set_facecolor("black")
    # ax.grid(False)
    # ax.w_xaxis.pane.fill = False
    # ax.w_yaxis.pane.fill = False
    # ax.w_zaxis.pane.fill = False
    # ax.set_box_aspect(aspect = (4,2,2))
    # ax.view_init(elev=5, azim=150, roll=0)

    # ---------------

    # Create the 4th color-rendered dimension
    cmap = cm.get_cmap('Greys')
    cmap.set_under("white")
    # cmap.set_over("red")
    scam = plt.cm.ScalarMappable(
        norm=cm.colors.LogNorm(g_min, g_gmax),
        cmap=cmap # see https://matplotlib.org/examples/color/colormaps_reference.html
    )


    # ax = fig.add_subplot(projection='3d')
    zmin = 0.
    xmin, xmax, ymin,ymax = 0,0,0,0
    for skymap in skymaps:
        X, Y = np.meshgrid(skymap.grid_x, skymap.grid_y)
        if (zmin > np.min(Y)): zmin = np.min(Y)
        if (xmin > np.min(X)): xmin = np.min(X)
        if (ymin > np.min(Y)): ymin = np.min(Y)
        if (xmax < np.max(X)): xmax = np.max(X)
        if (ymax < np.max(Y)): ymax = np.max(Y)

        Z = np.full_like(X,fill_value=np.log10(skymap.time / PBA.utils.cgs.day))
        G = skymap.im_hist.T
        # X = X[G > g_min]
        # Y = Y[G > g_min]
        Z[G < g_min] = np.nan
        G[G < g_min] = g_min
        scam.set_array([])
        facecolors = scam.to_rgba(G)
        # ax.scatter(Z, X, Y, color=facecolors)
        ax.plot_surface(
            Z, X, Y,
            facecolors  = facecolors,
            antialiased = True,
            rstride=1, cstride=1, alpha=.6, shade=False,

        )
    for skymap in skymaps:
        ax.plot(skymap.xc, zmin, zs=np.log10(skymap.time/PBA.utils.cgs.day), zdir='x', marker="o", ms=1, color="black")
        ax.plot([skymap.x1,skymap.x2], [zmin,zmin], zs=np.log10(skymap.time/PBA.utils.cgs.day),
                zdir='x', ls="--", lw=.6, color="black")

        zmin_ = -.5
        ax.plot(np.log10(skymap.time/PBA.utils.cgs.day), zmin_, zs=skymap.yc, zdir='z', marker="o", ms=1, color="black")
        ax.plot([np.log10(skymap.time/PBA.utils.cgs.day),np.log10(skymap.time/PBA.utils.cgs.day)],
                [zmin_,zmin_], zs=[skymap.y1,skymap.y2],
                zdir='z', ls="--", lw=.6, color="black")

        # ax.plot(skymap.grid_x, skymap.dist_x, zs=np.log10(skymap.time), zdir='x', label='curve in (x, y)')
        # ax.plot(skymap.grid_y, skymap.dist_y, zs=np.log10(skymap.time), zdir='z', label='curve in (x, y)')
    times = [np.log10(skymap.time/PBA.utils.cgs.day) for skymap in skymaps]
    # xmin = times.min(),
    # xmax = times.max()
    # n = len(times)
    # for t
    ax.set_ylim(0,11)
    ax.set_ylabel(r"$X$ [mas]",fontsize=12,labelpad=10)
    ax.set_zlabel(r"$Z$ [mas]",fontsize=12,labelpad=5)
    ax.set_xlabel(r"$\log(t_{\rm obs})$ [day]",fontsize=12,labelpad=18)
    ax.minorticks_on()
    # ax.set_facecolor("black")
    ax.grid(False)
    ax.w_xaxis.pane.fill = False
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False
    # ax.w_zaxis.line.set_visible(False)
    ax.set_box_aspect(aspect = (6,2,2))
    ax.view_init(elev=10, azim=150, roll=0)
    ax.tick_params(direction='in', length=10, width=2, colors='black',
                   grid_color='gray', grid_alpha=0.1,which="both", axis="both",labelsize=12)
    # x_scale=4
    # y_scale=1
    # z_scale=1
    #
    # scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    # scale=scale*(1.0/scale.max())
    # scale[3,3]=1.0
    # def short_proj():
    #     return np.dot(Axes3D.get_proj(ax), scale)
    # fig.set_facecolor('black')
    # fig.subplots_adjust(left=-.2,top=.2)  # plot outside the normal area
    # fig.tight_layout()

    # ax.get_proj=short_proj
    # ax.set_title(r"Gaussian jet with $\texttt{PW}$ method",loc='center',color="black")
    # plt.title(r"Gaussian jet with $\texttt{PW}$ method",loc='center',color="black")
    # ax.text(x=.8, y=.1, z=.6, s=r'Gaussian jet with $\texttt{PW}$ method', horizontalalignment='center',
    #      verticalalignment='center', transform=ax.transAxes)
    # ax.annotate(r'Gaussian jet with $\texttt{PW}$ method', xy=(2, 1), xytext=(-200, 200), textcoords='offset points', ha='left', bbox=dict(boxstyle='circle', fc='green', alpha=0.7),
    #          arrowprops=dict(arrowstyle='->'))
    fig.suptitle(r"Gaussian jet with $\texttt{PW}$ method", y=0.90, fontsize=16)

    plt.savefig(os.getcwd()+'/'+"3d_example.png",dpi=256)
    plt.show()


def task_kn_skymap_with_dist_one_time():

    time=200*day; freq=1e9

    workdir = os.getcwd()+'/'

    # prepare initial data
    struct = {"struct":"gaussian",
              "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
              "theta_c": 0.085, "theta_w": 0.2618}
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=100, n_layers_a=10)
    id_dict = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, outfpath=workdir+"gauss_grb_id.h5")

    # initialize interface
    pba = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True,parfile="parfile.par")
    pba.reload_parfile()

    # run
    pba.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out", loglevel="info")

    # process skymap
    conf = {"nx":64, "ny":32, "extend_grid":1.1, "fwhm_fac":0.5, "lat_dist_method":"integ",
            "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
            "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }}
    prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=True)
    prep.process_singles(infpaths=workdir+"raw_skymap_*.h5", outfpath=pba.GRB.fpath_sky_map, remove_input=True)

    # plot skymap
    config = {
            "gridspec": {
                "width_ratios": (4, 2), "height_ratios": (2, 4),
                "left": 0.14, "right": 0.95, "bottom": 0.1, "top": 0.96, "wspace": 0.05, "hspace": 0.05
            },
            "show_figs": True, "save_figs":False, "save_pdf":False, "figfpath":workdir,
            "grid": False,
            "figsize": (4.8, 4.8),
            "type" : "hist",
            "cm": {"color": 'cyan', "marker": "+"},
            "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                           "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan},
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
    PBA.skymap_plotting_tools.full_plot_skymap_with_hists(skymap=pba.GRB.get_skymap(time=time, freq=freq), conf=config)

    skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()[:-4]]
    # times = pba.GRB.get_skymap_times()
    # skymaps = [pba.GRB.get_skymap(times[0], freq=freq),
    #            pba.GRB.get_skymap(times[1], freq=freq),
    #            pba.GRB.get_skymap(times[2], freq=freq),
    #            pba.GRB.get_skymap(times[3], freq=freq),
    #            pba.GRB.get_skymap(times[4], freq=freq),
    #            pba.GRB.get_skymap(times[5], freq=freq),
    #            pba.GRB.get_skymap(times[6], freq=freq),
    #            ]

    plot_skymaps_3d(skymaps)


def main():
    task_kn_skymap_with_dist_one_time()

if __name__ == '__main__':
    main()