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
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
import os
from matplotlib import cm
import os

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import BPA_METHODS
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
    from PyBlastAfterglowMag.id_maker_from_kenta_bns import prepare_kn_ej_id_2d
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import BPA_METHODS
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
        from package.src.PyBlastAfterglowMag.id_maker_from_kenta_bns import prepare_kn_ej_id_2d
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")



def plot2(vals : dict, figpath = None):
    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(4.6,2+2.8), sharex="all")
    ax = axes[0]

    ax.plot(np.array(vals["texts"]), np.array(vals["v_ave"])*get_Gamma(vals["v_ave"]), 'x', color='blue')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=False,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_ylabel(r"Mass-averaged $\langle \Gamma\beta \rangle$", color="blue")
    # ax.set_xlabel(r"$t_{\rm ext}$ [ms]")

    ax1 = ax.twinx()
    ax1.plot(np.array(vals["texts"]), np.array(vals["v_fast_ave"])*get_Gamma(vals["v_fast_ave"]), 'o', color="red")
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=False, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax1.set_ylabel(r"Mass-averaged $\langle \Gamma\beta(\Gamma\beta>1) \rangle$",color="red")
    ax1.set_xlim(vals["texts"][0], vals["texts"][-1])
    # plt.tight_layout()
    # plt.savefig(FIGPATH+figppath+"average_velocity_evolution"+"png",dpi=256)
    # plt.show()

    # -------------------------------------------------------------------
    # fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4.6, 2.8))
    ax = axes[1]
    ax.plot(np.array(vals["texts"]), np.array(vals["theta_rms"]), 'x', color='blue', label="Total ejecta")
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=False,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_ylabel(r"RMS half-openning angle $\langle \theta_{\rm RMS} \rangle$", color="red")

    ax.plot(np.array(vals["texts"]), np.array(vals["fast_theta_rms"]), 'd', color="blue", label=r"Fast tail, $\Gamma\beta>1$")
    ax.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax.legend(**{"fancybox": False, "loc": 'lower right',
                 # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                 "shadow": "False", "ncol": 1, "fontsize": 12 - 2, "columnspacing": 0.4,
                 "framealpha": 0., "borderaxespad": 0., "frameon": False})
    ax1 = ax.twinx()
    ax1.plot(np.array(vals["texts"]), np.array(vals["fast_masses"]), 'x', color="red")
    ax1.minorticks_on()
    ax1.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=False, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.set_ylabel(r"$\langle M_{\rm ej}(\Gamma\beta>1)$", color="red")
    # ax1.set_yscale("log")
    ax.set_xlabel(r"$t_{\rm ext}$ [ms]")
    ax.set_ylabel(r"$\theta_{\rm RMS} = \sqrt{ \Sigma(m_i\theta_i)/\Sigma(m_i) }$", color="blue")
    ax.set_xlim(vals["texts"][0], vals["texts"][-1])
    plt.tight_layout()
    if not figpath is None: plt.savefig(figpath + ".png", dpi=256)
    if not figpath is None: plt.savefig(figpath + ".pdf")
    plt.show()

def plot_init_profile(ctheta, betas, eks,
                      xmin=0,xmax=90,ymin=1e-2,ymax=6,vmin=1e-12,vmax=1e-6,
                      title=None, figpath=None):
    fontsize=12
    gammas = get_Gamma(betas)
    moms = gammas * betas
    ctheta = ctheta * 180 / cgs.pi

    fig = plt.figure(figsize=(4.5 + 1, 3.6 + 3))
    # fig.suptitle(r"BLh* $(1.259+1.482)M_{\odot}$")

    ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.59 - 0.12])
    ax1 = fig.add_axes([0.16, 0.61, 0.81 - 0.15, 0.91 - 0.61])
    cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.59 - 0.12])

    #                   x1    y1    delta x       delta y
    # ax1 = fig.add_axes([0.16, 0.45, 0.81 - 0.15, 0.59 - 0.12])
    # ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.91 - 0.61])
    # cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.91 - 0.61])

    # top panel
    # for t in tasks:
    #     hist_vinf, hist_mass = t["data"].load_vinf_hist()
    #     hist_mom = hist_vinf * get_Gamma(hist_vinf)
    #     hist_eks = np.cumsum(np.array(0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m)[::-1])[::-1]

    ax1.plot(ctheta, np.sum(eks,axis=0), color='black', ls='-', drawstyle='steps')

    if (not title is None): ax1.set_title(title)

    ###  mooley
    # mool_mom = np.linspace(0.5 * get_Gamma(0.5), 0.8 * get_Gamma(0.8), 20)
    # mool_ek = 5e50 * (mool_mom / 0.4) ** -5
    # _l, = ax1.plot(mool_mom, mool_ek, color="gray", ls="-")
    # ax1.text(0.30, 0.90, "Mooley+17", color='black', transform=ax1.transAxes, fontsize=12)

    # tasks = [1, 1, 1]
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="-", label=r"$q=1.00$")
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls="--", label=r"$q=1.43$")
    # ax1.plot([-1., -1., ], [-1., -1., ], color="gray", ls=":", label=r"$q=1.82$")
    # han, lab = ax1.get_legend_handles_labels()
    # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
    #                          **{"fancybox": False, "loc": 'lower left',
    #                            # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                            "shadow": "False", "ncol": 1, "fontsize": 11,
    #                            "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # ax1.add_artist(ax1.legend(han[:-1 * len(tasks)], lab[:-1 * len(tasks)],
    #                           **{"fancybox": False, "loc": 'upper right',
    #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # "shadow": "False", "ncol": 1, "fontsize": 11,
    # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # ax1.add_artist(ax1.legend(han[len(han) - len(tasks):], lab[len(lab) - len(tasks):],
    #                           **{"fancybox": False, "loc": 'center right',
    #                              "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    # "shadow": "False", "ncol": 1, "fontsize": 11,
    # "framealpha": 0., "borderaxespad": 0., "frameon": False}))
    #
    # fit *= np.sum(mdens) / np.sum(fit)
    # fit_x, fit_y = _fit(mom, ekdens)
    # ax1.plot(fit_x, fit_y, 'k--', label=r'$\sin^2(\theta)$')
    # ax1.legend(**{"fancybox": False, "loc": 'upper right',
    #                # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                "shadow": "False", "ncol": 1, "fontsize": 11,
    #                "framealpha": 0., "borderaxespad": 0., "frameon": False})

    ax1.set_yscale("log")
    # ax1.set_ylim(ymin, ymax)
    # ax1.yaxis.tick_right()
    # ax1.yaxis.tick_left()
    ax1.set_ylabel(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    ax1.get_yaxis().set_label_coords(-0.15, 0.5)

    ax1.set_xlim(xmin, xmax)
    ax1.xaxis.set_ticklabels([])

    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.minorticks_on()

    ax11 = ax1.twinx()
    mask = betas*get_Gamma(betas)>1
    if len(betas[mask])>0:
        axx = np.sum(betas[mask, np.newaxis] * eks[mask, :],axis=0)
        ayy = np.sum(eks[mask, :],axis=0)
        # ax11.plot(ctheta, axx/ayy, color="gray", ls="-")
        ax11.plot(ctheta, np.sum(eks[mask, :],axis=0), color="gray", ls="-")
    ax11.set_yscale("log")
    ax11.set_ylabel(r"$M_{\rm ej}(\Gamma\beta>1)$ [M$_{\odot}$]", fontsize=fontsize, color="gray")
    ax11.minorticks_on()
    ax11.tick_params(axis='both', which='both', labelleft=False,
                     labelright=True, tick1On=False, tick2On=True,
                     labelsize=12,
                     direction='in',
                     bottom=True, top=True, left=True, right=True)


    # bottom panel
    # mask = vinf > 0.6
    # eks2 = np.zeros_like(eks)
    # for i in range(len(vinf)):
    #     if vinf[i] > 0.6:
    #         eks2[:, i] = eks[:, i]

    # import h5py
    # dset = h5py.File("/home/vsevolod/Desktop/Hajela/BLh_q100.h5", 'w')
    # dset.create_dataset(name="beta", data=vinf)
    # dset.create_dataset(name="theta", data=theta * 180 / np.pi)
    # dset.create_dataset(name="Ek(>Gamma_beta)", data=eks)
    # dset.close()

    levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
    cmap = plt.get_cmap('RdYlBu_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = LogNorm(eks[(eks > 0) & (np.isfinite(eks))].min(), eks[(eks > 0) & (np.isfinite(eks))].max())

    im = ax0.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

    ax0.axhline(y=1, linestyle='--', color='gray')

    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)

    ax0.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
    ax0.set_xlabel(r"Polar angle", fontsize=fontsize)
    ax0.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, lbl, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    # ax0.text(0.05, 0.05, models.print_fancy_label(name), color='white', transform=ax0.transAxes)

    # ekdens = np.sum(eks, axis=0)
    # ax1.step(0.5 * (theta[1:] + theta[:-1]), mdens, where='mid', color='red')
    # hist_vinf, hist_mass = o_data.load_vinf_hist()
    # hist_mom = hist_vinf * get_Gamma(hist_vinf)
    # hist_eks = 0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m
    # ax1.step(mom, ekdens, where='mid', color='red')
    # ax1.step(hist_mom, hist_eks, where='mid', color='black')

    cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar.set_label(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.ax.minorticks_off()
    # plt.savefig(PAPERPATH + "kinetic_energy_struct_models.pdf")
    # plt.savefig(FIGPATH + "kinetic_energy_struct_models.png")
    if not figpath is None: plt.savefig(figpath+".png", dpi=256)
    if not figpath is None: plt.savefig(figpath+".pdf")
    # plt.savefig(sys.argv[0].replace(".py", "_") + name.lower() + ".pdf")
    plt.show()
    plt.close()







    #
    #
    #
    #
    # # eks = np.log10(eks)
    # fig, ax = plt.subplots(figsize=(4.6, 2.8), ncols=1, nrows=2, sharex="all")
    #
    # # fig.add_axes([0.6, .1, .35, .3])
    #
    # ax = ax[1]
    #
    # # ax = plt.subplot(111, polar=True)
    # levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
    # cmap = plt.get_cmap('RdYlBu_r')
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    # norm = LogNorm(eks[(eks>0)&(np.isfinite(eks))].min(), eks[(eks>0)&(np.isfinite(eks))].max())
    #
    # im = ax.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")
    #
    # ax.set_xlabel(r"Polar angle, $\theta$ [deg]")
    # ax.set_ylabel(r"$\Gamma\beta$")
    # ax.set_yscale("log")
    # ax.set_ylim(1e-1, 4)
    # if (not title is None):
    #     ax.set_title(title)
    #
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.minorticks_on()
    #
    # # ax.set_yscale("log")
    #
    # # cbar.set_title("x")
    # # ax.set_title('pcolormesh with levels')
    # # print(ax.get_rmin(), ax.get_rmax())
    # # ax.set_rmax(10)
    # # ax.set_rmin(20)
    # # print(ax.get_rmin(), ax.get_rmax())
    #
    # # max_theta = 90
    # # ax.set_thetamax(max_theta)
    # # ax.set_rlim(10,20)
    #
    # # ax.set_rlim(Rs.min(), Rs.max())
    # # ax.set_rscale("log")
    # # ax.set_rscale('log')
    # #
    # # ticklabels = ax.get_yticklabels()
    # # labels = range(80, 0, -10)
    # # for i in range(0, len(labels)):
    # #     ticklabels[i] = str(labels[i])
    # # ax.set_yticklabels(ticklabels)
    # plt.tight_layout()
    # if (save_figs): plt.savefig(FIGPATH + figname + ".png", dpi=256)
    # # if (save_figs): plt.savefig(PAPERPATH + figname + ".pdf")
    # plt.show()

def main():
    workdir = os.getcwd()+'/'
    ###  locate and extract the data from the original ejecta profiles from Kenta
    # path_to_original_data = "/media/vsevolod/data/KentaData/SFHo_13_14_150m_11/" #
    # if not os.path.isdir(path_to_original_data):
    #     raise FileExistsError("Input ejecta data is not found")

    # lead the original data and prepare the ID file for PyBlastAfterglow
    text = 25.
    label = f"corr_id_SFHo_13_14_150m_11_text{int(text)}"
    # prepare_kn_ej_id_2d(datadir=path_to_original_data,
    #                     outfpaths=[workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5"],
    #                     req_times=np.array([text]),
    #                     new_theta_len=None,
    #                     verbose=True,
    #                     dist="pw")

    # data_par_list, data_par_list_all = get_ej_data_for_text(datadir=path_to_original_data,
    #                                                        req_times=np.array([text]),
    #                                                        new_theta_len=None,
    #                                                        verbose=True)
    # data_par_list = data_par_list[0]
    # plot_init_profile(data_par_list["thetas"], data_par_list["betas"], data_par_list["masses"].T,
    #                   figpath=None,#FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
    #                   title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(text)))
    modify_parfile_par_opt(workingdir=workdir, part="kn", newpars={},
                           newopts={"fname_light_curve":f"{label}.lc",
                                    "fname_ejecta_id":f"{label}.h5"},
                           parfile="parfile.par",newparfile="parfile.par",keep_old=False
                           )
    pba = BPA_METHODS(workingdir=os.getcwd()+'/',readparfileforpaths=True)
    pba.reload_parfile()
    pba.run(loglevel="info")

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes
    ax.loglog(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9),
            **{"color":"black", "ls": "--", "lw": 0.8, "label": label})

    pba.clear()
    plt.show()


if __name__ == '__main__':
    main()