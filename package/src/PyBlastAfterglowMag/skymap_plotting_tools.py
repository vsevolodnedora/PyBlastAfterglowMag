import copy
import os.path
import numpy as np
import h5py
import hashlib
from scipy import ndimage, interpolate
from glob import glob

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec


rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')
rcParams['font.size'] = 8
# plt.style.use('fivethirtyeight')

from .interface import Skymap


def plot_pcolomesh(ax, task_dict, int_x, int_y, int_zz, outfname=None):
    norm_min = task_dict["norm"][1]
    norm_max = task_dict["norm"][2]

    if not hasattr(norm_min, "__contains__"):
        pass
    else:
        if norm_min.__contains__("min"):
            norm_min = float(norm_min.replace("min", "")) * int_zz[np.isfinite(int_zz)].min()
        elif norm_min.__contains__("max"):
            norm_min = float(norm_min.replace("max", "")) * int_zz[np.isfinite(int_zz)].max()
        else:
            raise KeyError("norm_min does not contain 'min' or 'max' ")

    if not hasattr(norm_max, "__contains__"):
        pass
    else:
        if norm_max.__contains__("max"):
            norm_max = float(norm_max.replace("max", "")) * int_zz[np.isfinite(int_zz)].max()
        elif norm_max.__contains__("min"):
            norm_max = float(norm_max.replace("min", "")) * int_zz[np.isfinite(int_zz)].min()
        else:
            raise KeyError("norm_max does not contain 'min' or 'max' ")

    if task_dict["norm"][0] == "linear":
        norm = Normalize(norm_min, norm_max)
    elif task_dict["norm"][0] == "log":
        norm = LogNorm(norm_min, norm_max)
    elif task_dict["norm"][0] == "twoslope":
        norm = colors.TwoSlopeNorm(vmin=task_dict["norm"][1], vcenter=task_dict["norm"][2], vmax=task_dict["norm"][3])
    else:
        raise KeyError("norm is wrong")

    # levels = MaxNLocator(nbins=40).tick_values(int_ration.min(), int_ration.max())
    # # levels = MaxNLocator(nbins=40).tick_values(-5, 1)
    # cmap = plt.get_cmap(tmp["pcolormesh"]["cmap"])
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    cmap = cm.get_cmap(task_dict["cmap"])
    if not task_dict["set_under"] is None: cmap.set_under(task_dict["set_under"])
    if not task_dict["set_over"] is None: cmap.set_over(task_dict["set_over"])
    im = ax.pcolormesh(int_x, int_y, int_zz, norm=norm, cmap=cmap, alpha=float(task_dict["alpha"]))
    # im = ax_main.contourf(int_x, int_y, int_zz, norm=norm, cmap=cmap, alpha=1.0)
    # im.set_rasteraized(True)
    if not (task_dict["facecolor"] is None):
        ax.set_facecolor(cmap(float(task_dict["facecolor"])))

    ''' --------- SAVE FILE --------- '''
    # if not outfname is None:
    #     arr = copy.deepcopy(int_zz)
    #     arr[~np.isfinite(arr)] = 0.0
    #     print("max={} min={}".format(arr.max(), arr.min()))
    #     arr = arr / arr.max()
    #     # arr = arr.T
    #     out = np.zeros((len(int_x) * len(int_y), 3))
    #     i = 0
    #     for ix in range(len(int_x)):
    #         for iy in range(len(int_y)):
    #             out[i][0] = int_x[ix]
    #             out[i][1] = int_y[iy]
    #             out[i][2] = arr[iy, ix]
    #             i = i+1
    #     print(out.shape)
    #     np.savetxt(outfname + ".txt", out, fmt=("%.4f", "%.4f", "%.4e"), header="X[mas] Z[mas] I[normalized]",
    #                footer="# X has {} points Z has {} points".format(len(int_x), len(int_y)))
    #     # arr = int_zz/int_zz.max()
    #
    #     dfile = h5py.File(outfname + ".h5", "w")
    #     dfile.create_dataset(name="X[mas]", data=int_x)
    #     dfile.create_dataset(name="Z[mas]", data=int_y)
    #     dfile.create_dataset(name="intensity[normalized]", data=arr)
    #     dfile.close()
    #     exit(1)

    return im


def plot_skymap_with_hists(skymap, tmp, ax_main, ax_histx, ax_histy):
    # assert x2 > x1
    # ax.axvline(x=x1, color='gray', linestyle='dotted')
    # ax.axvline(x=x2, color='gray', linestyle='dotted')
    # ax.axvline(x=xcs_m, color='gray', linestyle='solid')
    if ("cm" in tmp.keys()) and len(tmp["cm"].keys()) > 0:
        # ---
        if (not (ax_main is None)): ax_main.plot(skymap.xc, skymap.yc, **tmp["cm"])
        if (not (ax_histx is None)): ax_histx.axvline(x=skymap.xc, color=tmp["cm"]["color"], linestyle='dashed')
        if (not (ax_histy is None)): ax_histy.axhline(y=skymap.yc, color=tmp["cm"]["color"], linestyle='dashed')
    if ("ysize" in tmp.keys()) and len(tmp["ysize"].keys()) > 0:
        if (not (ax_main is None)): ax_main.errorbar([skymap.xc, skymap.xc],
                                                     [skymap.y1, skymap.y2],
                                                     xerr=[skymap.grid_x.max() / 10, skymap.grid_x.max() / 10], **tmp["ysize"])
        if (not (ax_histy is None)): ax_histy.plot(skymap.dist_y * 1e3, skymap.grid_y, lw=1., ls='-',
                                                   drawstyle='steps', color=tmp["ysize"]["color"])
        if (not (ax_histy is None)): ax_histy.axhline(y=skymap.y2, color=tmp["ysize"]["color"], linestyle='dotted')
        if (not (ax_histy is None)): ax_histy.axhline(y=skymap.y1, color=tmp["ysize"]["color"], linestyle='dotted')
    if ("xsize" in tmp.keys()) and len(tmp["xsize"].keys()) > 0:

        if (not (ax_main is None)): ax_main.errorbar([skymap.x1, skymap.x2],
                                                     [skymap.yc, skymap.yc],
                                                     yerr=[skymap.grid_y.max() / 10, skymap.grid_y.max() / 10], **tmp["xsize"])
        if (not (ax_histx is None)): ax_histx.plot(skymap.grid_x, skymap.dist_x * 1e3, lw=1., ls='-', drawstyle='steps',
                                                   color=tmp["xsize"]["color"])  # , color='blue')
        if (not (ax_histx is None)): ax_histx.axvline(x=skymap.x2, color=tmp["xsize"]["color"], linestyle='dotted')
        if (not (ax_histx is None)): ax_histx.axvline(x=skymap.x1, color=tmp["xsize"]["color"], linestyle='dotted')
    # --------------------
    # if len(tmp["pcolormesh"].keys()) > 0:
    #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
    #     int_zz = int_zz.T
    #     int_zz[~np.isfinite(int_zz)] = 0.
    #     if (len(tmp["smooth"].keys()) > 0):
    #         int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"], sigma=tmp["smooth"]["sigma"])
    #     plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
    if ("pcolormesh" in tmp.keys()) and "pcolormesh" in tmp.keys():
        # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        int_zz = skymap.im_hist if tmp["type"] == "hist" else skymap.im_intp

        int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])

        im = plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"],
                                    int_x=skymap.grid_x, int_y=skymap.grid_y, int_zz=int_zz.T, outfname=None)
        if tmp["pcolormesh"]["set_rasterized"]: ax_main.set_rasterized(im)
        return im


def full_plot_skymap_with_hists(skymap : Skymap, conf : dict) -> None:
    settings = copy.deepcopy(conf)
    fig = plt.figure(figsize=settings["figsize"])
    gs = fig.add_gridspec(2, 2, **settings["gridspec"])
    # ---
    ax_main = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)
    ax_cbar = fig.add_subplot(gs[0, 1], sharey=ax_main)

    # plot the colormesh with histograms
    im = plot_skymap_with_hists(skymap, conf, ax_main, ax_histx, ax_histy)

    # make the plot pretty
    ax_histx.set_yscale("linear")
    ax_histx.set_xscale("linear")
    ax_histx.minorticks_on()
    ax_histx.set_yscale("log")
    # ax_histx.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    ax_histx.set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)  # , color='blue')
    ax_histx.tick_params(axis='both', which='both', labelleft=True, labelbottom=False,
                         labelright=False, tick1On=True, tick2On=True,
                         labelsize=12, color="white",
                         direction='in',
                         bottom=False, top=True, left=True, right=True)
    ax_histx.set_facecolor(settings["histx_backgound_color"])
    ax_histy.set_yscale("linear")
    ax_histy.set_xscale("linear")
    ax_histy.minorticks_on()
    ax_histy.set_xscale("log")
    ax_histy.set_xlim(ax_histx.get_ylim())
    # ax_histy.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    ax_histy.set_xlabel(r"$I_{\nu;\rm m}(z)$", fontsize=12)  # , color='blue')
    ax_histy.tick_params(axis='both', which='both', labelleft=False,
                         labelright=False, tick1On=True, tick2On=True,
                         labelsize=12, color="white",
                         direction='in',
                         bottom=True, top=True, left=True, right=False)
    ax_histy.set_facecolor(settings["histy_backgound_color"])
    ax_main.set_yscale("linear")
    ax_main.set_xscale("linear")
    ax_main.minorticks_on()
    ax_main.set_xlabel(settings["xlabel"], fontsize=12)
    ax_main.set_ylabel(settings["ylabel"], fontsize=12)
    # ax_main.set_xlim(int_x.min() * 1.1, int_x.max() * 1.1)
    # ax_main.set_ylim(int_y.min() * 1.1, int_y.max() * 1.1)
    if "xlim" in settings.keys() and len(settings["xlim"]) == 2: ax_main.set_xlim(settings["xlim"])
    else: ax_main.set_xlim(skymap.grid_x.min() * 1.1, skymap.grid_x.max() * 1.1)
    if "ylim" in settings.keys() and len(settings["ylim"]) == 2: ax_main.set_ylim(settings["ylim"])
    else: ax_main.set_ylim(skymap.grid_y.min() * 1.1, skymap.grid_y.max() * 1.1)
    ax_main.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12,
                        direction='in',
                        bottom=True, top=True, left=True, right=True)
    ax_main.minorticks_on()
    ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
    ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
    ax_main.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, color="white",
                        direction='in',
                        bottom=True, top=True, left=True, right=True)

    # if "text" in settings.keys() and settings["text"] == "full":
    #     tot_flux = float(self.ej.get_skymap_totfluxes(freq) * 1e3)
    #     ax_main.text(1.05, 0.6, r"$\\I_{\rm max}$="
    #                  + r"${}$".format(latex_float(int_zz.max() * 1e3)) + r"\,[$\mu$Jy/mas$^2$]\\"
    #                  + r"$F_{\nu}=$" + r"${}$".format(latex_float(tot_flux)) + r"\,[$\mu$Jy]\\"
    #                  + r" $\nu=$" + r"${}$".format(latex_float(freq)) + r"\,[Hz] \\"
    #                  + r"$\theta_{\rm obs}=" + r"{:.1f}".format(theta_obs) + r"$ [deg]\\"
    #                  + r"$t_{\rm obs}=" + r"{:.1f}".format(time) + "$ [day]"
    #                  , color='black', transform=ax.transAxes, fontsize=fontsize - 2)

    if settings["plot_grids"]:
        ax_main.grid()
        ax_histx.grid()
        ax_histy.grid()
    ax_histx.set_ylim(*settings["histx_lim"])
    ax_histy.set_xlim(*settings["histy_lim"])

    # ax.tick_params(  # axis='both',
    #     which='both',  # labelleft=True,
    #     # labelright=False,
    #     tick1On=True, tick2On=True,
    #     labelsize=12,
    #     direction='in',
    #     color="black",
    #     #               bottom=True, top=True, left=True, right=True
    # )
    # ax.set_xlim(settings["xlim"])
    # ax.set_ylim(settings["ylim"])
    # if settings["grid"]:
    #     ax.set_xlim(min(settings["xlim"]) * 1.1, max(settings["xlim"]) * 1.1)
    #     ax.set_ylim(min(settings["ylim"]) * 1.1, max(settings["ylim"]) * 1.1)
    #
    # # ax.tick_params(axis='both', which='both', labelleft=True,
    # #                labelright=False, tick1On=True, tick2On=True,
    # #                labelsize=11, color="white",
    # #                direction='in',
    # #                bottom=True, top=True, left=True, right=True)
    # # ax.set_title(r"$t_{\rm obs}=$" + r"{:d} day".format(time))
    # # ax.set_title(line["label"])
    #
    # ax.axhline(y=0, linestyle=':', linewidth=0.4, color="black")
    # ax.axvline(x=0, linestyle=':', linewidth=0.4, color="black")
    #
    # if (len(settings["title"].keys()) > 0):
    #     if settings["title"]["title"] == "time_fluxratio":
    #         times = pb.get_ej_skymap_times()
    #         fluxes = pb.get_ej_skymap_totfluxes(freq=freq)
    #         fluxes_w = pb_w.get_ej_skymap_totfluxes(freq=freq)
    #         fluxes_ration = fluxes_w / fluxes
    #         idx = find_nearest_index(times, time * cgs.day)
    #         # ax.set_title("t={:.0f} d. F/F={}".format(time, fluxes_ration[idx]))
    #         ax.set_title("t={:.0f} d. ".format(time) +
    #                      r"$F_{\nu}^{\rm w}/F_{\nu}^{\rm w/o}=$" +
    #                      "{:.2f}".format(fluxes_ration[idx]))
    #     elif settings["title"]["title"] == "time":
    #         ax.set_title("t={:.0f} d. ".format(time))
    #
    # # ------------------------------------------------------------------------------------------------------------------
    #
    #
    # ''' ------------------------------------- '''
    #
    # fig = plt.figure(figsize=(5, 5))
    #
    # # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
    # # the size of the marginal axes and the main axes in both directions.
    # # Also adjust the subplot parameters for a square plot.
    # gs = fig.add_gridspec(2, 2, width_ratios=(4, 2), height_ratios=(2, 4),
    #                       left=0.13, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)
    #
    # ax_main = fig.add_subplot(gs[1, 0])
    # ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
    # ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)
    #
    # ''' -------------------------------------- '''
    #
    # # nx = 200
    # # ny = 100
    # # nx = np.complex(0, nx + 1)
    # # ny = np.complex(0, ny + 1)
    # # edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
    # # edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
    # # # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    # # #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    # # i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
    # # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    # # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    # # tmp = []
    # # for i in range(len(grid_x)):
    # #     tmp.append(np.sum(i_zz[i,:] * np.diff(edges_y))/np.sum(np.diff(edges_y)))
    # # tmp = np.array(tmp)
    # #
    # # y_sum_ii = np.sum(i_zz * 1e3, axis=1)
    # # max_ii = np.max(y_sum_ii)
    # # x1 = grid_x[np.argmin(y_sum_ii < max_ii * 0.5)]
    # # x2 = grid_x[::-1][np.argmin(y_sum_ii[::-1] < max_ii * 0.5)]
    #
    # grid_x, i_zz_x, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y")
    # xc, yc = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
    # x2, x1 = pb.get_skymap_fwhm(grid_x, i_zz_x, xc)
    #
    # ax = ax_histx
    # ax.plot(grid_x, i_zz_x * 1e3, lw=1., ls='-', drawstyle='steps', color="black")  # , color='blue')
    # # ax.plot(grid_x, np.sum(i_zz * 1e3, axis=1) / np.sum(np.diff(edges_y)), lw=1., ls='--', drawstyle='steps')#, color='blue')
    # # ax.plot(int_x, int_zz.max(axis=1), lw=1., ls='-')
    # # ax.plot(int_x, np.sum(int_zz * int_y, axis=1)/np.sum(int_y,axis=1), lw=1., ls='--')
    # tot_flux = float(pb.get_ej_skymap_totfluxes(freq) * 1e3)
    # ax.text(1.05, 0.6, r"$\\I_{\rm max}$="
    #         + r"${}$".format(latex_float(int_zz.max() * 1e3)) + r"\,[$\mu$Jy/mas$^2$]\\"
    #         + r"$F_{\nu}=$" + r"${}$".format(latex_float(tot_flux)) + r"\,[$\mu$Jy]\\"
    #         + r" $\nu=$" + r"${}$".format(latex_float(freq)) + r"\,[Hz] \\"
    #         + r"$\theta_{\rm obs}=" + r"{:.1f}".format(theta_obs) + r"$ [deg]\\"
    #         + r"$t_{\rm obs}=" + r"{:.1f}".format(time) + "$ [day]"
    #         , color='black', transform=ax.transAxes, fontsize=fontsize - 2)
    # ax.set_yscale("log")
    # ax.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    # ax.tick_params(axis='both', which='both', labelleft=True, labelbottom=False,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize,
    #                direction='in',
    #                bottom=False, top=True, left=True, right=True)
    # # ax1=ax.twinx()
    # # ax1.plot(grid_x, tmp, lw=1., ls='-', drawstyle='steps', color='red')
    # # ax1.set_yscale("log")
    # # ax1.set_ylabel(r"$\sum_{y} I \Delta y / \sum_{y} \Delta y$ [$\mu$Jy/mas$^2$]", fontsize=12, color='red')
    # ax.axvline(x=x1, color='gray', linestyle='dotted')
    # ax.axvline(x=x2, color='gray', linestyle='dotted')
    # ax.axvline(x=xcs_m, color='gray', linestyle='solid')
    #
    # ''' --------------------------------------- '''
    #
    # # nx = 200
    # # ny = 100
    # # nx = np.complex(0, nx + 1)
    # # ny = np.complex(0, ny + 1)
    # # edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
    # # edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
    # # # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    # # #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    # # i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
    # # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    # # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    # # tmp = []
    # # for i in range(len(grid_y)):
    # #     tmp.append(np.sum(i_zz[:, i] * np.diff(edges_x)) / np.sum(np.diff(edges_x)))
    # # tmp = np.array(tmp)
    # #
    # # x_sum_ii = np.sum(i_zz * 1e3, axis=0)
    # # max_ii = np.max(x_sum_ii)
    # # y1 = grid_y[np.argmin(x_sum_ii < max_ii * 0.5)]
    # # y2 = grid_y[::-1][np.argmin(x_sum_ii[::-1] < max_ii * 0.5)]
    # # y1, y2, grid_y, i_zz_y = pb.get_ej_skymap_fwhm(all_x, all_y, all_fluxes,axis="z")
    # grid_y, i_zz_y, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="x")
    # xc, yc = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
    # y2, y1 = pb.get_skymap_fwhm(grid_y, i_zz_y, yc)
    #
    # ax = ax_histy
    # ax.plot(i_zz_y * 1e3, grid_y, lw=1., ls='-', drawstyle='steps', color='black')
    # # ax.plot(int_x, int_zz.max(axis=1), lw=1., ls='-')
    # # ax.plot(int_x, np.sum(int_zz * int_y, axis=1)/np.sum(int_y,axis=1), lw=1., ls='--')
    # ax.set_xscale("log")
    # ax.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    # ax.tick_params(axis='both', which='both', labelleft=False,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=False)
    # # ax1 = ax.twiny()
    # # ax1.plot(tmp, grid_y, lw=1., ls='-', drawstyle='steps', color='red')
    # # ax1.set_xscale("log")
    # # ax1.set_xlabel(r"$\sum_{y} I \Delta y / \sum_{y} \Delta y$ [$\mu$Jy/mas$^2$]", fontsize=12, color='red')
    # ax.axhline(y=y1, color='gray', linestyle='dotted')
    # ax.axhline(y=y2, color='gray', linestyle='dotted')
    # ax.axhline(y=ycs_m, color='gray', linestyle='solid')
    #
    # # ax.yaxis.set_ticklabels([])
    #
    # ''' --------------------------------------- '''
    #
    # # my_norm = LogNorm(1e-12, 1e-8)
    # # my_norm = LogNorm(1e-2, 1e0)
    # # my_norm = Normalize(int_zz.min(),int_zz.max())
    # # my_norm = LogNorm(int_zz.max()*1e-2,int_zz.max())
    # cmap = cm.get_cmap('viridis')
    #
    # ax = ax_main
    #
    # # if not (title is None): ax.set_title(title)
    # im = _plot_image(ax, int_x, int_y, int_zz, LogNorm(int_zz.max() * 1e-2, int_zz.max()),
    #                  levels=np.geomspace(int_zz.max() * 1e-2, int_zz.max(), 50))
    #
    # # center of mass
    # ax.plot(xcs_m, ycs_m, "x", color='white')
    #
    # # _i = np.argmax(x_latAvDist)
    # # xm = (int_x)[np.argmax(x_latAvDist), np.argmax(y_latAvDist)]
    # # ym = (int_y)[np.argmax(x_latAvDist), np.argmax(y_latAvDist)]
    # # ax.plot(_x, fwhm, ls='none', marker='.')
    # # ax.scatter([xm], [ym], s=80, marker='o', color="white", facecolors='none')
    # # _fwhm = fwhm
    #
    # # im.cmap.set_over("black")
    # # im.cmap.set_over("gray")
    # #
    # ax.set_yscale("linear")
    # ax.set_xscale("linear")
    # #
    # ax.set_xlabel("X [mas]", fontsize=12)
    # ax.set_ylabel("Z [mas]", fontsize=12)
    # #
    # # ax.set_xlim(grid_x.min(), grid_x.max())
    # # ax.set_ylim(grid_y.min(), grid_y.max())
    # #
    # ax.set_xlim(int_x.min() * 1.1, int_x.max() * 1.1)
    # ax.set_ylim(int_y.min() * 1.1, int_y.max() * 1.1)
    # #
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.minorticks_on()
    # ax.axhline(y=0, linestyle=':', linewidth=0.4)
    # ax.axvline(x=0, linestyle=':', linewidth=0.4)
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize, color="white",
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # --------- CBAR ----------------
    divider = make_axes_locatable(ax_cbar)
    cax = divider.append_axes('left', size='99%', pad=0.9)
    plt.delaxes(ax_cbar)
    cbar = plt.colorbar(im, cax=cax,
                        # format='%.0e',ticks=ticks
                        orientation='vertical',
                        # label=r"$I_{\nu}$ [mJy/mas$^2$]"
                        )
    cbar.set_label(label=r"$I_{\nu}$ [mJy/mas$^2$]", fontsize=12)
    cax.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    # print("plotted: \n")
    # plt.tight_layout()
    figname = settings["figfpath"]  # make_prefix_for_pars(pars)
    if settings["save_figs"]:
        print(r"Saving:\n{}".format(figname + ".png"))
        plt.savefig(".png", dpi=256)
    if settings["save_pdf"]:
        print(r"Saving:\n{}".format(figname + ".pdf"))
        plt.savefig(figname + ".pdf")
    if settings["show_figs"]: plt.show()

