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


def run(working_dir: str, struct: dict, P: dict, type: str = "a", run: bool = True) -> PBA.PyBlastAfterglow:
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
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80,
                                       n_layers_a=1 if struct["struct"] == "tophat" else 20)

    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_pw.h5")

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_a.h5")

    # create new parfile
    P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
    grb_skymap_config = copy.deepcopy( P["grb"]["skymap_conf"] )
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
    if (pba.GRB.opts["do_skymap"] == "yes"):
        prep = PBA.skymap_process.ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=False)

    return pba

def plot_skymaps_3d(ax : plt.axes, skymaps : list[PBA.Skymap]):

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
            facecolors  = facecolors,
            antialiased = True, # True
            rstride=1, cstride=1, alpha=.01, shade=False,

        )
    for skymap in skymaps:
        ax.plot(skymap.xc, zmin, zs=np.log10(skymap.time/cgs.day), zdir='x', marker="o", ms=1, color="black")
        ax.plot([skymap.x1,skymap.x2], [zmin,zmin], zs=np.log10(skymap.time/cgs.day),
                zdir='x', ls="--", lw=.6, color="black")

        zmin_ = -.5
        ax.plot(np.log10(skymap.time/cgs.day), zmin_, zs=skymap.yc, zdir='z', marker="o", ms=1, color="black")
        ax.plot([np.log10(skymap.time/cgs.day),np.log10(skymap.time/cgs.day)],
                [zmin_,zmin_], zs=[skymap.y1,skymap.y2],
                zdir='z', ls="--", lw=.6, color="black")

        # ax.plot(skymap.grid_x, skymap.dist_x, zs=np.log10(skymap.time), zdir='x', label='curve in (x, y)')
        # ax.plot(skymap.grid_y, skymap.dist_y, zs=np.log10(skymap.time), zdir='z', label='curve in (x, y)')
    times = [np.log10(skymap.time/cgs.day) for skymap in skymaps]
    # xmin = times.min(),
    # xmax = times.max()
    # n = len(times)
    # for t
    ax.set_ylim(0,4)
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

def plot_3d_skmap_stack(do_run: bool, plot: bool, struct: dict, P: dict, name="tophat_skymap"):

    pba = run(
        working_dir=os.getcwd() + f"/working_dirs/{name}/",
        struct=struct, P=d2d(default=P, new=dict(
            main=dict(theta_obs=np.pi/4.),
            grb=dict(#method_ssc_fs='numeric',
                     #use_ssa_fs='yes'
            ))),
        type="a", run=do_run
    )

    freq = pba.GRB.get_skymap_freqs()[0]
    skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()[3::][:-6]]
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
                           gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)

    plot_skymaps_3d(ax=ax,skymaps=skymaps)
    plt.show()



    freq = pars["obs_freq"]
    pars_a = copy.deepcopy(pars)
    pars_a["mom0_frac_when_start_spread"] = 0.95
    opts_a["do_skymap"] = "yes"
    opts_a["fname_sky_map"]=f"skymap_a_for_3d_plot.h5"
    pba_a = self.run_a(struct=struct, pars=pars_a, opts={}, opts_grb=opts_a)

    skymaps = [pba_a.GRB.get_skymap(time, freq=freq) for time in pba_a.GRB.get_skymap_times()]
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
                           gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)

    plot_skymaps_3d(ax=ax,skymaps=skymaps)

    fig.suptitle(r"Gaussian jet with $\texttt{A}$ method", y=0.90, fontsize=16)
    if not figpath is None:
        _figpath = figpath+"_A"
        print("Saving:\n {}".format(_figpath))
        if save_pdf: plt.savefig(_figpath+".pdf")
        plt.savefig(_figpath+".png", dpi=256)
    if show_fig: plt.show()

def plot_skymaps(do_run: bool, task: dict, struct: dict, P: dict, name="tophat_skymap"):

    pba = run(
        working_dir=os.getcwd() + f"/working_dirs/{name}/",
        struct=struct, P=d2d(default=P, new=dict(
            main=dict(theta_obs=np.pi/4.),
            grb=dict(#method_ssc_fs='numeric',
                #use_ssa_fs='yes'
            ))),
        type="a", run=do_run
    )

    freq = pba.GRB.get_skymap_freqs()[0]
    skymaps = [pba.GRB.get_skymap(time, freq=freq) for time in pba.GRB.get_skymap_times()]

    fig,ax = plt.subplots(ncols=1,nrows=1)


    ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.xc for skymap in skymaps],marker='.',ls='-')
    # ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.yc for skymap in skymaps],marker='.',ls='-')
    ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.x2-skymap.x1 for skymap in skymaps],marker='.',ls='-')
    # ax.plot([skymap.time/cgs.day for skymap in skymaps],[skymap.y2-skymap.y1 for skymap in skymaps],marker='.',ls='none')
    print(pba.GRB.get_skymap_times())

    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=4, fontsize=12, labelcolor='black',
              framealpha=0.8, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_ylabel(r"$x_c$", fontsize=12)
    ax.set_xlabel(r"$t_{\rm obs}$ [s]", fontsize=12)
    # ax.set_title(task["title"], fontsize=12)
    # ax.set_xlim(*task['xlim'])
    # plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    # plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    # if task["show"]: plt.show()
    # if plot: plt.show()
    # plt.close(fig)

    plt.show()


if __name__ == '__main__':
    do_run = False
    plot = True
    struct = dict(struct="tophat", Eiso_c=1.e53, Gamma0c=400., M0c=-1., theta_c=0.1, theta_w=0.1)
    # struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    P = dict(
        main=dict(n_ism=1., tb0=3e3, ntb=1000, rtol=5e-7, theta_obs=0,
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
                 method_synchrotron_fs="Joh06",
                 method_synchrotron_rs="Joh06",
                 # ngam_fs=1001,gam1_rs=1,gam2_rs=1e4,ngam_rs=1001,
                 # eps_b_fs = 1e-7,
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 # method_synchrotron_fs="Bessel",
                 # method_synchrotron_rs="Bessel",
                 method_ele_fs='analytic',
                 method_ele_rs="analytic",
                 # num_ele_use_adi_loss_fs='no',
                 # num_ele_use_adi_loss_rs='no',
                 gam1_fs=1., gam2_fs=1e8, ngam_fs=401,
                 gam1_rs=1., gam2_rs=1e8, ngam_rs=401,
                 freq1_fs=1e6, freq2_fs=1e32, nfreq_fs=401,
                 freq1_rs=1e6, freq2_rs=1e32, nfreq_rs=401,
                 # ebl_tbl_fpath="none"
                 skymap_conf=dict(nx=64,ny=32,extend_grid=2,fwhm_fac=0.5,lat_dist_method="integ",
                                  intp_filter=dict(type=None, sigma=2, mode='reflect'),  # "gaussian"
                                  hist_filter=dict(type=None,sigma=2, mode='reflect'))
                 )
    )

    # plot_3d_skmap_stack(do_run=do_run, plot=plot, struct=struct, P=P)
    plot_skymaps(do_run=do_run, task=dict(), struct=struct, P=P)