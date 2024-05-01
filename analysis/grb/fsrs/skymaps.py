import copy

from more_itertools import numeric_range

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs

import os, shutil, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, SymLogNorm, TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator, AutoLocator
from matplotlib import ticker
import matplotlib.ticker as plticker
import numpy as np
from scipy import special
from scipy import integrate
from matplotlib import cm
from package.src.PyBlastAfterglowMag import Ejecta

working_dir = os.getcwd() + '/tmp1/'
fig_dir = os.getcwd() + '/figs/'


def d2d(default: dict, new: dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_


def run(working_dir: str, struct: dict, P: dict, type: str = "a", run: bool = True,
        process_skymaps: bool = True) -> PBA.PyBlastAfterglow:
    """
            conf = {"nx": 64, "ny": 32, "extend_grid": 2, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                "intp_filter": {"type": None, "sigma": 2, "mode": 'reflect'},  # "gaussian"
                "hist_filter": {"type": None, "sigma": 2, "mode": 'reflect'}}
    :param working_dir:
    :param struct:
    :param P:
    :param type:
    :param run:
    :return:
    """
    # clean he temporary direcotry
    if run and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    # generate initial data for blast waves
    struct = copy.deepcopy(struct)
    pba_id = PBA.id_analytic.JetStruct(
        n_layers_pw=80 if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
        n_layers_a=(1 if struct["struct"] == "tophat" else
                    (20 if not "n_layers_a" in struct.keys() else struct["n_layers_a"])))

    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_pw.h5")

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_a.h5")

    # create new parfile
    P = copy.deepcopy(P)
    P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
    P["grb"]["method_eats"] = "piece-wise" if type == "pw" else "adaptive"
    if (struct["struct"]=="tophat"): P["grb"]["nsublayers"] = 35 # for skymap resolution
    grb_skymap_config = copy.deepcopy(P["grb"]["skymap_conf"])
    del P["grb"]["skymap_conf"]
    PBA.parfile_tools.create_parfile(working_dir=working_dir, P=P)

    # instantiate PyBlastAfterglow
    pba = PBA.interface.PyBlastAfterglow(workingdir=working_dir)

    # run the code with given parfile
    if run:
        pba.run(
            path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
            loglevel="info"
        )

    # process skymap
    if (process_skymaps and pba.GRB.opts["do_skymap"] == "yes"):
        prep = PBA.skymap_process.ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=False)

    return pba


def plot_skymaps_3d(ax: plt.axes, skymaps: list[PBA.Skymap]):
    # g_min = np.min([np.min(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_gmax = np.max([np.max(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_min = g_gmax * 1e-5
    print(f"g_min={g_min} g_gmax={g_gmax}")

    # fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
    #                        gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)

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
    cmap = cm.get_cmap('jet')
    cmap.set_under("white")
    # cmap.set_over("red")
    scam = plt.cm.ScalarMappable(
        norm=cm.colors.LogNorm(g_min, g_gmax),
        cmap=cmap  # see https://matplotlib.org/examples/color/colormaps_reference.html
    )

    # ax = fig.add_subplot(projection='3d')
    zmin = 0.
    xmin, xmax, ymin, ymax = 0, 0, 0, 0
    for skymap in skymaps:
        X, Y = np.meshgrid(skymap.grid_x, skymap.grid_y)
        if (zmin > np.min(Y)): zmin = np.min(Y)
        if (xmin > np.min(X)): xmin = np.min(X)
        if (ymin > np.min(Y)): ymin = np.min(Y)
        if (xmax < np.max(X)): xmax = np.max(X)
        if (ymax < np.max(Y)): ymax = np.max(Y)

        Z = np.full_like(X, fill_value=np.log10(skymap.time / PBA.utils.cgs.day))
        G = skymap.im_hist.T
        # X = X[G > g_min]
        # Y = Y[G > g_min]
        Z[G < g_min] = np.nan
        G[G < g_min] = g_min
        Y[X < 0] = 0.
        # Z[X < 0] = 0.
        G[X < 0] = 0.
        X[X < 0] = 0.

        scam.set_array([])
        facecolors = scam.to_rgba(G)
        # ax.scatter(Z, X, Y, color=facecolors)
        ax.plot_surface(
            Z, X, Y,
            facecolors=facecolors,
            antialiased=True,  # True
            rstride=1, cstride=1, alpha=.01, shade=False,

        )
    for skymap in skymaps:
        ax.plot(skymap.xc, zmin, zs=np.log10(skymap.time / cgs.day), zdir='x', marker="o", ms=1, color="black")
        ax.plot([skymap.x1, skymap.x2], [zmin, zmin], zs=np.log10(skymap.time / cgs.day),
                zdir='x', ls="--", lw=.6, color="black")

        zmin_ = -.5
        ax.plot(np.log10(skymap.time / cgs.day), zmin_, zs=skymap.yc, zdir='z', marker="o", ms=1, color="black")
        ax.plot([np.log10(skymap.time / cgs.day), np.log10(skymap.time / cgs.day)],
                [zmin_, zmin_], zs=[skymap.y1, skymap.y2],
                zdir='z', ls="--", lw=.6, color="black")

        # ax.plot(skymap.grid_x, skymap.dist_x, zs=np.log10(skymap.time), zdir='x', label='curve in (x, y)')
        # ax.plot(skymap.grid_y, skymap.dist_y, zs=np.log10(skymap.time), zdir='z', label='curve in (x, y)')
    times = [np.log10(skymap.time / cgs.day) for skymap in skymaps]
    # xmin = times.min(),
    # xmax = times.max()
    # n = len(times)
    # for t
    ax.set_ylim(0, 4)
    ax.set_ylabel(r"$X$ [mas]", fontsize=12, labelpad=10)
    ax.set_zlabel(r"$Z$ [mas]", fontsize=12, labelpad=5)
    ax.set_xlabel(r"$\log(t_{\rm obs})$ [day]", fontsize=12, labelpad=18)
    ax.minorticks_on()
    # ax.set_facecolor("black")
    ax.grid(False)
    ax.w_xaxis.pane.fill = False
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False
    # ax.w_zaxis.line.set_visible(False)
    ax.set_box_aspect(aspect=(6, 2, 2))
    ax.view_init(elev=10, azim=150, roll=0)
    ax.tick_params(direction='in', length=10, width=2, colors='black',
                   grid_color='gray', grid_alpha=0.1, which="both", axis="both", labelsize=12)
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


def plot_3d_skmap_stack(do_run: bool, plot: bool, struct: dict, P: dict, name="tophat_skymap"):
    pba = run(
        working_dir=os.getcwd() + f"/working_dirs/{name}/",
        struct=struct, P=d2d(default=P, new=dict(
            main=dict(theta_obs=np.pi / 4.),
            grb=dict(  #method_ssc_fs='numeric',
                #use_ssa_fs='yes'
            ))),
        type="a", run=do_run
    )

    freq = pba.GRB.get_skymap_freqs()[0]
    skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()[3::][:-6]]
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(5.2, 4),
                           gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1))  # figsize=(4.2,4)figsize=(5,4)

    plot_skymaps_3d(ax=ax, skymaps=skymaps)
    plt.show()

    freq = pars["obs_freq"]
    pars_a = copy.deepcopy(pars)
    pars_a["mom0_frac_when_start_spread"] = 0.95
    opts_a["do_skymap"] = "yes"
    opts_a["fname_sky_map"] = f"skymap_a_for_3d_plot.h5"
    pba_a = self.run_a(struct=struct, pars=pars_a, opts={}, opts_grb=opts_a)

    skymaps = [pba_a.GRB.get_skymap(time, freq=freq) for time in pba_a.GRB.get_skymap_times()]
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize=(5.2, 4),
                           gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1))  # figsize=(4.2,4)figsize=(5,4)

    plot_skymaps_3d(ax=ax, skymaps=skymaps)

    fig.suptitle(r"Gaussian jet with $\texttt{A}$ method", y=0.90, fontsize=16)
    if not figpath is None:
        _figpath = figpath + "_A"
        print("Saving:\n {}".format(_figpath))
        if save_pdf: plt.savefig(_figpath + ".pdf")
        plt.savefig(_figpath + ".png", dpi=256)
    if show_fig: plt.show()


def plot_skymap(do_run: bool, process_skymaps: bool, task: dict, struct: dict, P: dict, name="tophat_skymap_hists"):
    workingdir = f"/working_dirs/{name}/"
    pba = run(
        working_dir=os.getcwd() + workingdir,
        struct=struct, P=d2d(default=P, new=dict(
            main=dict(theta_obs=np.pi / 4.),
            grb=dict(do_rs="yes",bw_type="fsrs",
                #method_ssc_fs='numeric',
                #use_ssa_fs='yes'
                # nsublayers=35
            ))),
        type='a', run=do_run, process_skymaps=process_skymaps
    )
    # freq = pba.GRB.get_skymap_freqs()[0]
    # skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
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
        "ysize": {"capsize": 2, "color": "black", "lw": 0.5},
        "xsize": {"capsize": 2, "color": "black", "lw": 0.5},
        "pcolormesh": {"cmap": 'jet', "set_under": 'white', "set_over": None, "set_rasterized": True,
                       "norm": ("log", "0.001max", "1max"), "facecolor": None, "alpha": 1.0, "isnan": np.nan},
        "xlim": (-1., 1.0), "ylim": (-1.0, 1.0),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "white",
        "histy_backgound_color": "white",
        "plot_grids":True,
        "histx_lim":(1e-4,1e-2),
        "histy_lim":(1e-4,1e-2)
    }
    time= 500.*cgs.day
    freq= 1.e9
    PBA.skymap_plotting_tools.full_plot_skymap_with_hists(skymap=pba.GRB.get_skymap(time=time, freq=freq), conf=config)

def plot_skymap_properties(do_run: bool, process_skymaps: bool, task: dict, struct: dict, P: dict, name="tophat_skymap_hists"):
    fig, axes = plt.subplots(ncols=1, nrows=2, sharex='all', figsize=(5, 3), layout='constrained')

# ------------ RESOLUTION / METHODS ---------------
def plot_skymaps_comparison_tophat(do_run: bool, process_skymaps: bool, struct: dict, P: dict):
    freq= 1.e9
    name="tophat_skymap_props"
    # task1 = dict(plot=dict(ls='-',color='black',marker='.',label='Default'),
    #              workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa_ssc/",
    #              pars=dict(main=dict(theta_obs=np.pi / 4.),
    #                        grb=dict(do_rs="yes",bw_type="fsrs",do_rs_radiation='no',method_ssc_fs='numeric', use_ssa_fs='yes')))
    # task2 = dict(plot=dict(ls='-',color='blue',marker='.',label=r'FS \& RS'),
    #              workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc_fsrs/",
    #              pars=dict(main=dict(theta_obs=np.pi / 4.),
    #                        grb=dict(do_rs="yes",bw_type="fsrs",method_ssc_fs='numeric', use_ssa_fs='yes')))
    task1 = dict(plot=dict(ls='-',color='black',marker='.',label='Default'),
                 workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa_ssc_fs/",
                 pars=dict(main=dict(),
                           grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes',method_ele_fs='mix')))
    task2 = dict(plot=dict(ls='-',color='blue',marker='.',label=r'FS \& RS'),
                 workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc_fs/",
                 pars=dict(main=dict(),
                           grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes')))

    pba1 = run(
        working_dir=task1["workingdir"], struct=struct, P=d2d(default=P, new=task1["pars"]),
        type='a', run=do_run, process_skymaps=False
    )

    pba2 = run(
        working_dir=task2["workingdir"], struct=struct, P=d2d(default=P, new=task2["pars"]),
        type='a', run=do_run, process_skymaps=False
    )

    if (process_skymaps):
        prep = PBA.skymap_process.ProcessRawSkymap(conf=copy.deepcopy(P["grb"]["skymap_conf"]), verbose=False)
        prep.process_multiples(infpaths=(pba1.workingdir+"raw_skymap_*.h5",
                                         pba2.workingdir+"raw_skymap_*.h5"),
                             outfpaths =(pba1.GRB.fpath_sky_map,
                                         pba2.GRB.fpath_sky_map))




    # PBA.skymap_plotting_tools.full_plot_skymap_with_hists(skymap=pba.GRB.get_skymap(time=time, freq=freq), conf=config)
    # plt.imshow(im)


    fig, axes = plt.subplots(ncols=4, nrows=2, sharex='col', sharey='row', figsize=(9, 4), layout='constrained')
    # norm = SymLogNorm(linthresh=task.get("vmin", 1e-6),
    #                   vmin=task.get("vmin", 1e-6),
    #                   vmax=task.get("vmax", 1e6), base=10)
    # norm = LogNorm(vmin=task.get("vmin", skymap2.im_hist[skymap2.im_hist>0].max()),
    #                vmax=task.get("vmax", skymap2.im_hist[skymap2.im_hist>0].max()*1e-6))
    # norm = LogNorm(vmin=task.get("vmin", 0.5),
    #                vmax=task.get("vmax", 1.5))
    # cmap = plt.get_cmap(task.get("cmap", 'jet'))
    # cmap.set_under(task.get("set_under", 'white'), alpha=0.2)
    # cmap.set_over(task.get("set_over", 'white'), alpha=0.2)


    times = [25,35,50,75]
    # im = im[~np.isfinite(im)] = np.nan
    for i, (ax_hist,ax_main) in enumerate(axes.T):
        time= times[i]*cgs.day

        skymap1 = pba1.GRB.get_skymap(time=time, freq=freq)
        skymap2 = pba2.GRB.get_skymap(time=time, freq=freq)



        ax_hist.plot(skymap2.grid_x, skymap1.dist_x * 1e3, color='blue', lw=1., ls='-', drawstyle='steps', label='FS')
        ax_hist.plot(skymap2.grid_x, (skymap2.dist_x-skymap1.dist_x) * 1e3, color='green', lw=1., ls='-', drawstyle='steps',label='RS')
        ax_hist.axvline(x=skymap2.xc,color='black',linewidth=0.7,linestyle='--',label='$x_c$')
        ax_hist.set_yscale('log')

        # ax_hist.set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)

        # im = im[~np.isfinite(im)] = np.nan
        # norm = LogNorm(vmin=task.get("vmin", 1e-4),  vmax=task.get("vmax", 1e2))
        norm = LogNorm(vmin=skymap2.im_hist[skymap2.im_hist>0].max()*1e-4,
                       vmax=skymap2.im_hist[skymap2.im_hist>0].max())
        cmap = plt.get_cmap('jet')
        _c = ax_main.pcolormesh(skymap1.grid_x, skymap2.grid_y, skymap2.im_hist.T, cmap=cmap, norm=norm)
        _c.set_rasterized(True)
        ax_main.plot([skymap1.xc,],[skymap1.yc,],marker='o',color='cyan',fillstyle='none',label='FS $x_c$',linestyle='none')
        ax_main.plot([skymap2.xc,],[skymap2.yc,],marker='s',color='cyan',fillstyle='none',label='FS \& RS $x_c$',linestyle='none')
        # ax_main.errorbar([skymap2.x1, skymap2.x2],
        #                  [skymap2.yc, skymap2.yc],
        #                  yerr=[skymap2.grid_y.max() / 10, skymap2.grid_y.max() / 10], **dict(capsize=2, color="white",lw=0.5))
        # ax_main.errorbar([skymap2.xc, skymap2.xc],
        #                  [skymap2.y1, skymap2.y2],
        #                  xerr=[skymap2.grid_x.max() / 30, skymap2.grid_x.max() / 30], **dict(capsize=2, color="white",lw=0.5))
        # val1 = np.array(((skymap2.im_hist-skymap1.im_hist) > skymap1.im_hist), dtype=np.float64)
        # val1[val1 == 0] = np.nan
        # # val1[val1 == 1] = np.nan
        # cs = ax_main.contourf(skymap1.grid_x, skymap2.grid_y, val1.T, hatches=['..'], alpha=0.0, edgecolor="r", facecolor="blue",
        #                  color="green")  #,zorder=0

        # ax_main.set_xlabel("$X$ [mas]", fontsize=12)
        # ax_main.set_ylabel("$Z$ [mas]", fontsize=12)
        ax_hist.grid(ls=':')
        ax_main.axvline(x=0,color='black',linestyle='dotted',linewidth=0.5)
        ax_main.axhline(y=0,color='black',linestyle='dotted',linewidth=0.5)

        ax_main.minorticks_on()
        ax_hist.minorticks_on()
        ax_hist.tick_params(axis='both', which='both', direction='in', labelsize=11,tick1On=True, tick2On=True,)
        ax_main.tick_params(axis='both', which='both', direction='in', labelsize=11,tick1On=True, tick2On=True,)

        ax_main.set_xlim(skymap2.grid_x.min(),skymap2.grid_x.max()*.75)
        ax_main.set_ylim(skymap2.grid_y.min()*.75,skymap2.grid_y.max()*.75)

        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=1.)
        ax_main.text(0.4, 0.1, "{:.0f} days".format(round(skymap1.time/cgs.day,0)), fontsize=12, bbox=bbox,
                transform=ax_main.transAxes, horizontalalignment='right')

    cbar = fig.colorbar(_c, ax=ax_main, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.ax.set_title(r"$I/I_{\rm max}$", size=12)

    axes[0][0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
    axes[1][0].set_ylabel("$Z$ [mas]", fontsize=12)
    for ax_main in axes[1][:]:
        ax_main.set_xlabel("$X$ [mas]", fontsize=12)
    axes[0][0].legend(fancybox=False, loc='upper left', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=1, fontsize=12, labelcolor='black',
              framealpha=0.0, borderaxespad=0.)
    axes[1][0].legend(fancybox=False, loc='upper right', columnspacing=0.8,
                      # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      shadow=False, ncol=1, fontsize=12, labelcolor='black',
                      framealpha=0.0, borderaxespad=0.)

    # Create the new axis for marginal X and Y labels
    # ax = fig.add_subplot(111, frameon=False)

    # Disable ticks. using ax.tick_params() works as well
    # ax.set_xticks([])
    # ax.set_yticks([])

    # Set X and Y label. Add labelpad so that the text does not overlap the ticks
    # ax.set_xlabel("The X label", labelpad=40, fontsize=12)
    # ax.set_ylabel("The Y label", labelpad=40, fontsize=12)


    plt.savefig(os.getcwd() + '/figs/' + name + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + name + '.pdf')
    if plot: plt.show()
    plt.close(fig)


    name="tophat_skymap_props2"


    times = [25,35,50,75]
    fig, axes = plt.subplots(ncols=1, nrows=3, sharex='all', figsize=(5, 5), layout='constrained')

    skymaps1 = [pba1.GRB.get_skymap(time, freq=freq) for time in pba1.GRB.get_skymap_times()]
    skymaps2 = [pba2.GRB.get_skymap(time, freq=freq) for time in pba2.GRB.get_skymap_times()]

    axes[0].plot([skymap.time / cgs.day for skymap in skymaps1],
                 [skymap2.xc-skymap1.xc for (skymap1,skymap2) in zip(skymaps1, skymaps2)],
                 color='gray', fillstyle='none', ls='-', label=r'FS')# marker='o',
    axes[1].plot([skymap.time / cgs.day for skymap in skymaps1],
                 [abs(skymap2.x2 - skymap2.x1)-abs(skymap1.x2 - skymap1.x1) for (skymap1,skymap2) in zip(skymaps1,skymaps2)],
                 color='gray', fillstyle='none', ls='-', label=r'FS \& RS')
    axes[2].plot([skymap.time / cgs.day for skymap in skymaps1],
                 [skymap1.flux for skymap1,skymap2 in zip(skymaps1,skymaps2)],
                 color='gray',  fillstyle='none', ls='-', label=r'FS')
    axes[2].plot([skymap.time / cgs.day for skymap in skymaps1],
                 [skymap2.flux for skymap1,skymap2 in zip(skymaps1,skymaps2)],
                 color='black',fillstyle='none', ls='-', label=r'FS \& RS')
    axes[2].set_ylim(1e-4,10)
    # axes[0].plot([skymap.time / cgs.day for skymap in skymaps1],
    #              [skymap.xc for skymap in skymaps1],
    #              color='gray', marker='o', fillstyle='none', ls='none', label=r'FS')
    # axes[1].plot([skymap.time / cgs.day for skymap in skymaps1],
    #              [abs(skymap.x2 - skymap.x1) for skymap in skymaps1],
    #              color='gray', marker='o', fillstyle='none', ls='none', label=r'FS \& RS')
    # axes[2].plot([skymap.time / cgs.day for skymap in skymaps1],
    #              [skymap.flux for skymap in skymaps1],
    #              color='gray', marker='o', fillstyle='none', ls='none', label=r'FS \& RS')

    # axes[0].plot([skymap.time / cgs.day for skymap in skymaps2],
    #              [skymap.xc for skymap in skymaps2],
    #              color='black', marker='s', fillstyle='none', ls='none', label=r'FS')
    # axes[1].plot([skymap.time / cgs.day for skymap in skymaps2],
    #              [abs(skymap.x2 - skymap.x1) for skymap in skymaps2],
    #              color='black', marker='s', fillstyle='none', ls='none', label=r'FS \& RS')
    # axes[2].plot([skymap.time / cgs.day for skymap in skymaps2],
    #              [skymap.flux for skymap in skymaps2],
    #              color='black', marker='s', fillstyle='none', ls='none', label=r'FS \& RS')
    # workingdir = f"/working_dirs/{name+'_num_ssa_ssc_pp'}/"
    # pba = run(
    #     working_dir=os.getcwd() + workingdir,
    #     struct=struct, P=d2d(default=P, new=dict(
    #         main=dict(theta_obs=np.pi / 4.),
    #         grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
    #     type='a', run=do_run, process_skymaps=process_skymaps
    # )
    # freq = pba.GRB.get_skymap_freqs()[0]
    # skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
    # axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps], color='black',
    #              marker='.', ls='-')
    # axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
    #              [abs(skymap.x2 - skymap.x1) for skymap in skymaps],
    #              color='black', marker='.', ls='-')


    # _l1, = axes[0].plot([0, 0], [1, 1], color='black', marker='.', ls='-')
    # _l2, = axes[0].plot([0, 0], [1, 1], color='blue', marker='.', ls='-')
    # _l3, = axes[0].plot([0, 0], [1, 1], color='green', marker='.', ls='-')
    # legend1 = axes[0].legend([_l1, _l2, _l3], ['Default', "No SSC", "No SSA"],
    #                          fancybox=False, loc='upper center', columnspacing=0.8,
    #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #                          framealpha=0.0, borderaxespad=0.)

    # _l1, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='-')
    # _l2, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='--')
    # # _l3, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls=':')
    # legend2 = axes[0].legend([_l1, _l2], ['SR', "HR"],
    #                          fancybox=False, loc='upper left', columnspacing=0.8,
    #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #                          framealpha=0.0, borderaxespad=0.)

    # axes[0].add_artist(legend1)
    # axes[0].add_artist(legend2)

    for (i, ax) in enumerate(axes):
        ax.set_xscale('log')
        ax.grid(ls=':')
        # ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        # ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
        #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #           shadow=False, ncol=4, fontsize=12, labelcolor='black',
        #           framealpha=0.8, borderaxespad=0.)
        # ax.set_rasterized(True)
    axes[0].set_ylabel(r"$x_c$ [mas]", fontsize=12)
    axes[1].set_ylabel(r"$\Delta_x$ [mas]", fontsize=12)
    axes[2].set_ylabel(r"$F_{\nu}$ [mJy]", fontsize=12)
    axes[2].set_yscale("log")

    axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]", fontsize=12)
    axes[-1].legend(fancybox=False, loc='upper left', columnspacing=0.8,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol=4, fontsize=12, labelcolor='black',
                   framealpha=0.0, borderaxespad=0.)
    # ax.set_title(task["title"], fontsize=12)
    # ax.set_xlim(*task['xlim'])
    axes[0].set_title("Top-hat jet",fontsize=14)
    plt.savefig(os.getcwd() + '/figs/' + name + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + name + '.pdf')
    if plot: plt.show()
    plt.close(fig)

    #
    #
    #
    # name="tophat_skymap_props"
    # show=True
    # tasks = [
    #     dict(plot=dict(ls='-',color='black',marker='.',label='Default'),
    #          workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa_ssc/",
    #          pars=dict(main=dict(theta_obs=np.pi / 4.),
    #                    grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
    #     dict(plot=dict(ls='-',color='blue',marker='.',label=r'FS \& RS'),
    #          workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc_fsrs/",
    #          pars=dict(main=dict(theta_obs=np.pi / 4.),
    #                    grb=dict(do_rs="yes",bw_type="fsrs",method_ssc_fs='numeric', use_ssa_fs='yes'))),
    #
    #     # dict(plot=dict(ls='-',color='blue',marker='.',label='Mix'),
    #     #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc/",
    #     #      pars=dict(main=dict(theta_obs=np.pi / 4.),
    #     #                grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes',method_ele_fs='mix'))),
    #
    #     # dict(plot=dict(ls='-',color='blue',marker='.',label='No SSC'),
    #     #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa/",
    #     #      pars=dict(main=dict(theta_obs=np.pi / 4.),
    #     #                grb=dict(method_ssc_fs='none', use_ssa_fs='yes'))),
    #     # dict(plot=dict(ls='-',color='green',marker='.',label='No SSA'),
    #     #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssc/",
    #     #      pars=dict(main=dict(theta_obs=np.pi / 4.),
    #     #                grb=dict(method_ssc_fs='numeric', use_ssa_fs='no'))),
    # ]
    #
    # fig, axes = plt.subplots(ncols=1, nrows=2, sharex='all', figsize=(5, 3), layout='constrained')
    # for task in tasks:
    #     pba = run(
    #         working_dir=task["workingdir"], struct=struct, P=d2d(default=P, new=task["pars"]),
    #         type='a', run=do_run, process_skymaps=process_skymaps
    #     )
    #     freq = pba.GRB.get_skymap_freqs()[0]
    #     skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
    #     axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps],
    #                  color=task["plot"]["color"],
    #                  marker=task["plot"]["marker"],
    #                  ls=task["plot"]["ls"],
    #                  label=task["plot"]["label"])
    #     axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
    #                  [abs(skymap.x2 - skymap.x1) for skymap in skymaps],
    #                  color=task["plot"]["color"],
    #                  marker=task["plot"]["marker"],
    #                  ls=task["plot"]["ls"])
    #
    # # workingdir = f"/working_dirs/{name+'_num_ssa_ssc_pp'}/"
    # # pba = run(
    # #     working_dir=os.getcwd() + workingdir,
    # #     struct=struct, P=d2d(default=P, new=dict(
    # #         main=dict(theta_obs=np.pi / 4.),
    # #         grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
    # #     type='a', run=do_run, process_skymaps=process_skymaps
    # # )
    # # freq = pba.GRB.get_skymap_freqs()[0]
    # # skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
    # # axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps], color='black',
    # #              marker='.', ls='-')
    # # axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
    # #              [abs(skymap.x2 - skymap.x1) for skymap in skymaps],
    # #              color='black', marker='.', ls='-')
    #
    #
    # # _l1, = axes[0].plot([0, 0], [1, 1], color='black', marker='.', ls='-')
    # # _l2, = axes[0].plot([0, 0], [1, 1], color='blue', marker='.', ls='-')
    # # _l3, = axes[0].plot([0, 0], [1, 1], color='green', marker='.', ls='-')
    # # legend1 = axes[0].legend([_l1, _l2, _l3], ['Default', "No SSC", "No SSA"],
    # #                          fancybox=False, loc='upper center', columnspacing=0.8,
    # #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    # #                          framealpha=0.0, borderaxespad=0.)
    #
    # # _l1, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='-')
    # # _l2, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='--')
    # # # _l3, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls=':')
    # # legend2 = axes[0].legend([_l1, _l2], ['SR', "HR"],
    # #                          fancybox=False, loc='upper left', columnspacing=0.8,
    # #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    # #                          framealpha=0.0, borderaxespad=0.)
    #
    # # axes[0].add_artist(legend1)
    # # axes[0].add_artist(legend2)
    #
    # for (i, ax) in enumerate(axes):
    #     ax.set_xscale('log')
    #     ax.grid(ls=':')
    #     # ax.set_yscale('log')
    #     ax.minorticks_on()
    #     ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    #     # ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
    #     #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #     #           shadow=False, ncol=4, fontsize=12, labelcolor='black',
    #     #           framealpha=0.8, borderaxespad=0.)
    #     # ax.set_rasterized(True)
    # axes[0].set_ylabel(r"$x_c$ [mas]", fontsize=12)
    # axes[1].set_ylabel(r"$\Delta_x$ [mas]", fontsize=12)
    # axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]", fontsize=12)
    # axes[0].legend(fancybox=False, loc='upper left', columnspacing=0.8,
    #                 # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                 shadow=False, ncol=4, fontsize=12, labelcolor='black',
    #                 framealpha=0.0, borderaxespad=0.)
    # # ax.set_title(task["title"], fontsize=12)
    # # ax.set_xlim(*task['xlim'])
    # axes[0].set_title("Top-hat jet",fontsize=14)
    # plt.savefig(os.getcwd() + '/figs/' + name + '.png', dpi=256)
    # plt.savefig(os.getcwd() + '/figs/' + name + '.pdf')
    # if show: plt.show()
    # if plot: plt.show()
    # plt.close(fig)

def plot_skymaps_comparison_tophat_evol(do_run: bool, process_skymaps: bool, task: dict, struct: dict, P: dict):

    name="tophat_skymap_props"
    show=True
    tasks = [
        dict(plot=dict(ls='-',color='orange',marker='.',label='Default'),
             workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa_ssc/",
             pars=dict(main=dict(theta_obs=np.pi / 4.),
                       grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
        dict(plot=dict(ls='-',color='blue',marker='.',label=r'FS \& RS'),
             workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc_fsrs/",
             pars=dict(main=dict(theta_obs=np.pi / 4.),
                       grb=dict(do_rs="yes",bw_type="fsrs",method_ssc_fs='numeric', use_ssa_fs='yes'))),

        # dict(plot=dict(ls='-',color='blue',marker='.',label='Mix'),
        #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_mix_ssa_ssc/",
        #      pars=dict(main=dict(theta_obs=np.pi / 4.),
        #                grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes',method_ele_fs='mix'))),

        # dict(plot=dict(ls='-',color='blue',marker='.',label='No SSC'),
        #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssa/",
        #      pars=dict(main=dict(theta_obs=np.pi / 4.),
        #                grb=dict(method_ssc_fs='none', use_ssa_fs='yes'))),
        # dict(plot=dict(ls='-',color='green',marker='.',label='No SSA'),
        #      workingdir=os.getcwd()+f"/working_dirs/{name}"+"_num_ssc/",
        #      pars=dict(main=dict(theta_obs=np.pi / 4.),
        #                grb=dict(method_ssc_fs='numeric', use_ssa_fs='no'))),
    ]
    times = [25,35,50,75]
    fig, axes = plt.subplots(ncols=1, nrows=2, sharex='all', figsize=(5, 3), layout='constrained')
    for task in tasks:
        pba = run(
            working_dir=task["workingdir"], struct=struct, P=d2d(default=P, new=task["pars"]),
            type='a', run=do_run, process_skymaps=process_skymaps
        )
        freq = pba.GRB.get_skymap_freqs()[0]
        skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
        axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps],
                     color=task["plot"]["color"],
                     marker=task["plot"]["marker"],
                     ls=task["plot"]["ls"],
                     label=task["plot"]["label"])
        axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
                     [abs(skymap.x2 - skymap.x1) for skymap in skymaps],
                     color=task["plot"]["color"],
                     marker=task["plot"]["marker"],
                     ls=task["plot"]["ls"])

    # workingdir = f"/working_dirs/{name+'_num_ssa_ssc_pp'}/"
    # pba = run(
    #     working_dir=os.getcwd() + workingdir,
    #     struct=struct, P=d2d(default=P, new=dict(
    #         main=dict(theta_obs=np.pi / 4.),
    #         grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
    #     type='a', run=do_run, process_skymaps=process_skymaps
    # )
    # freq = pba.GRB.get_skymap_freqs()[0]
    # skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
    # axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps], color='black',
    #              marker='.', ls='-')
    # axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
    #              [abs(skymap.x2 - skymap.x1) for skymap in skymaps],
    #              color='black', marker='.', ls='-')


    # _l1, = axes[0].plot([0, 0], [1, 1], color='black', marker='.', ls='-')
    # _l2, = axes[0].plot([0, 0], [1, 1], color='blue', marker='.', ls='-')
    # _l3, = axes[0].plot([0, 0], [1, 1], color='green', marker='.', ls='-')
    # legend1 = axes[0].legend([_l1, _l2, _l3], ['Default', "No SSC", "No SSA"],
    #                          fancybox=False, loc='upper center', columnspacing=0.8,
    #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #                          framealpha=0.0, borderaxespad=0.)

    # _l1, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='-')
    # _l2, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='--')
    # # _l3, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls=':')
    # legend2 = axes[0].legend([_l1, _l2], ['SR', "HR"],
    #                          fancybox=False, loc='upper left', columnspacing=0.8,
    #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                          shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #                          framealpha=0.0, borderaxespad=0.)

    # axes[0].add_artist(legend1)
    # axes[0].add_artist(legend2)

    for (i, ax) in enumerate(axes):
        ax.set_xscale('log')
        ax.grid(ls=':')
        # ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        # ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
        #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #           shadow=False, ncol=4, fontsize=12, labelcolor='black',
        #           framealpha=0.8, borderaxespad=0.)
        # ax.set_rasterized(True)
    axes[0].set_ylabel(r"$x_c$ [mas]", fontsize=12)
    axes[1].set_ylabel(r"$\Delta_x$ [mas]", fontsize=12)
    axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]", fontsize=12)
    axes[0].legend(fancybox=False, loc='upper left', columnspacing=0.8,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol=4, fontsize=12, labelcolor='black',
                   framealpha=0.0, borderaxespad=0.)
    # ax.set_title(task["title"], fontsize=12)
    # ax.set_xlim(*task['xlim'])
    axes[0].set_title("Top-hat jet",fontsize=14)
    plt.savefig(os.getcwd() + '/figs/' + name + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + name + '.pdf')
    if show: plt.show()
    if plot: plt.show()
    plt.close(fig)
def plot_skymaps_comparison_gaussian(do_run: bool, process_skymaps: bool, task: dict, struct: dict, P: dict, name="guass_skymap"):
    fig, axes = plt.subplots(ncols=1, nrows=2, sharex='all', figsize=(5, 3), layout='constrained')

    for (color, eats_type) in zip(['blue', 'green'], ["a", "pw"]):
        for (ls, (res_a, res_pw)) in zip(['-', '--'], [(21, 80), (51, 160)]):
            struct["n_layers_pw"] = res_pw
            struct["n_layers_a"] = res_a
            workingdir = f"/working_dirs/{name + '_' + eats_type + '_' + (str(res_a) if eats_type == 'a' else str(res_pw))}/"
            pba = run(
                working_dir=os.getcwd() + workingdir,
                struct=struct, P=d2d(default=P, new=dict(
                    main=dict(theta_obs=np.pi / 4.),
                    grb=dict(
                        #method_ssc_fs='numeric',
                        #use_ssa_fs='yes'
                        # nsublayers=res_a
                    ))),
                type=eats_type, run=do_run, process_skymaps=process_skymaps
            )
            freq = pba.GRB.get_skymap_freqs()[0]
            skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]
            axes[0].plot([skymap.time / cgs.day for skymap in skymaps], [skymap.xc for skymap in skymaps], color=color,
                         marker='.', ls=ls)
            # ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.yc for skymap in skymaps],marker='.',ls='-')
            axes[1].plot([skymap.time / cgs.day for skymap in skymaps],
                         [abs(skymap.x2 - skymap.x1) for skymap in skymaps], color=color, marker='.', ls=ls)
            # ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.y2-skymap.y1 for skymap in skymaps],marker='.',ls='none')

    _l1, = axes[0].plot([0, 0], [1, 1], color='blue', marker='.', ls='-')
    _l2, = axes[0].plot([0, 0], [1, 1], color='green', marker='.', ls='-')
    legend1 = axes[0].legend([_l1, _l2], ['A', "PW"],
                             fancybox=False, loc='upper center', columnspacing=0.8,
                             # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                             shadow=False, ncol=1, fontsize=12, labelcolor='black',
                             framealpha=0.0, borderaxespad=0.)

    _l1, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='-')
    _l2, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls='--')
    # _l3, = axes[0].plot([0, 0], [1, 1], color='gray', marker='.', ls=':')
    legend2 = axes[0].legend([_l1, _l2], ['SR', "HR"],
                             fancybox=False, loc='upper left', columnspacing=0.8,
                             # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                             shadow=False, ncol=1, fontsize=12, labelcolor='black',
                             framealpha=0.0, borderaxespad=0.)

    axes[0].add_artist(legend1)
    # axes[0].add_artist(legend2)

    for (i, ax) in enumerate(axes):
        ax.set_xscale('log')
        ax.grid(ls=':')
        # ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        # ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
        #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #           shadow=False, ncol=4, fontsize=12, labelcolor='black',
        #           framealpha=0.8, borderaxespad=0.)
        # ax.set_rasterized(True)
    axes[0].set_ylabel(r"$x_c$ [mas]", fontsize=12)
    axes[1].set_ylabel(r"$\Delta_x$ [mas]", fontsize=12)
    axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]", fontsize=12)
    axes[0].set_title("Gaussian jet",fontsize=14)

    # axes[0].legend(fancybox=False, loc='upper right', columnspacing=0.8,
        #                 # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #                 shadow=False, ncol=4, fontsize=12, labelcolor='black',
        #                 framealpha=0.8, borderaxespad=0.)
        # ax.set_title(task["title"], fontsize=12)
        # ax.set_xlim(*task['xlim'])
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)

if __name__ == '__main__':
    do_run = False
    process_skymaps = False
    plot = True
    struct = dict(struct="tophat", Eiso_c=1.e53, Gamma0c=400., M0c=-1., theta_c=0.1, theta_w=0.1)
    # struct = dict(struct="gaussian",Eiso_c=1.e53, Gamma0c=400., M0c= -1.,theta_c= 0.1, theta_w= 0.3)
    P = dict(
        main=dict(n_ism=1., tb0=3e3, ntb=1000, rtol=5e-7, theta_obs=np.pi / 4.,
                  # lc_freqs='array logspace 1e8 1e29 96',
                  # lc_times='array logspace 3e3 1e10 128',
                  skymap_freqs='array 1.e9 2.4e9 5.e9',
                  # skymap_times='array 86400 864000 3.456e6 8.64e6 1.728e7 2.592e+7 3.456e+7 5.184e+7 6.912e+7 8.64e+7 1.1232e+8 1.3824e+8 1.728e+8 2.16e+8 2.592e+8 3.456e+8 4.32e+8 6.048e+8 8.64e+8',
                  #                   0.5   1     2      4      8      16       32       48       64      96       128      172      256      364      512      742      900      1024
                  # skymap_times='array 43200 86400 172800 345600 691200 1.382e+6 2.765e+6 4.147e+6 5.53e+6 8.294e+6 1.106e+7 1.486e+7 2.212e+7 3.145e+7 4.424e+7 6.411e+7 7.776e+7 8.8474e+7',
                  skymap_times='array logspace 3e4 3e8 32',
                  ),
        grb=dict(save_dynamics='no', save_spec='no', do_lc='no', do_skymap='yes',
                 # method_nonrel_dist_fs='none',
                 # method_nonrel_dist_rs='none',
                 eps_e_fs=0.1, eps_b_fs=0.001, p_fs=2.2,
                 gamma_max_fs=4e7, method_gamma_max_fs="useConst",
                 gamma_max_rs=4e7, method_gamma_max_rs="useConst",
                 max_substeps_fs=1000, max_substeps_rs=1000,
                 # method_synchrotron_fs="Joh06",
                 # method_synchrotron_rs="Joh06",
                 # ngam_fs=1001,gam1_rs=1,gam2_rs=1e4,ngam_rs=1001,
                 # eps_b_fs = 1e-7,
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 # method_synchrotron_fs="Bessel",
                 # method_synchrotron_rs="Bessel",
                 # method_ele_fs='analytic', method_ne_fs="usenprime",
                 # method_ele_rs="analytic", method_ne_rs="usenprime",
                 # num_ele_use_adi_loss_fs='no',
                 # num_ele_use_adi_loss_rs='no',
                 gam1_fs=1., gam2_fs=1e8, ngam_fs=401,
                 gam1_rs=1., gam2_rs=1e8, ngam_rs=401,
                 freq1_fs=1e6, freq2_fs=1e32, nfreq_fs=401,
                 freq1_rs=1e6, freq2_rs=1e32, nfreq_rs=401,
                 # ebl_tbl_fpath="none"
                 skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
                                  intp_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect'),  # "gaussian"
                                  hist_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect'))
                 )
    )

    # plot_skymap(do_run=do_run, process_skymaps=process_skymaps,
    #             task=dict(figname="skymap_props_compare", show=True),
    #             struct=struct, P=P)
    plot_skymaps_comparison_tophat(do_run=do_run, process_skymaps=process_skymaps,
                                   struct=struct, P=P)
    # plot_skymaps_comparison_tophat_evol(do_run=do_run, process_skymaps=process_skymaps,
    #                                     task=dict(show=True),
    #                                     struct=struct, P=P)

    # ------------------------------------------------------------------

    P["grb"]["method_synchrotron_fs"] = "Joh06"
    P["grb"]["method_synchrotron_rs"] = "Joh06"
    P["grb"]["method_ele_fs"] = "analytic"
    P["grb"]["method_ele_rs"] = "analytic"
    P["grb"]["method_ne_fs"] = "usenprime"
    P["grb"]["method_ne_rs"] = "usenprime"

# plot_3d_skmap_stack(do_run=do_run, plot=plot, struct=struct, P=P)
    # plot_skymaps_comparison_tophat(do_run=do_run, process_skymaps=process_skymaps,
    #              task=dict(figname="skymap_props_compare", show=True),
    #              struct=struct, P=P)
    #
    # struct = dict(struct="gaussian", Eiso_c=1.e53, Gamma0c=400., M0c=-1., theta_c=0.1, theta_w=0.3)
    # plot_skymaps_comparison_gaussian(do_run=do_run, process_skymaps=process_skymaps,
    #                                task=dict(figname="gauss_skymap_props_compare", show=True),
    #                                struct=struct, P=P)
