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
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc, rcParams
from matplotlib.patches import Rectangle
import matplotlib.colors as colors
from matplotlib import cm
import os
from matplotlib import cm
import os
import re
from glob import glob

# try:
#     from PyBlastAfterglowMag.interface import modify_parfile_par_opt
#     from PyBlastAfterglowMag.interface import PyBlastAfterglow
#     from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#     from PyBlastAfterglowMag.utils import \
#         (latex_float, cgs, get_beta, get_Gamma,get_Beta,BetFromMom,GamFromMom,MomFromGam)
#     from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     from PyBlastAfterglowMag.id_maker_from_kenta_bns import \
#         prepare_kn_ej_id_2d, load_init_data
#     from PyBlastAfterglowMag.skymap_tools import \
#         (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import \
#             (latex_float, cgs, get_beta, get_Gamma,get_Beta,BetFromMom,GamFromMom,MomFromGam)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#         from package.src.PyBlastAfterglowMag.id_maker_from_kenta_bns import \
#             prepare_kn_ej_id_2d, load_init_data
#         from package.src.PyBlastAfterglowMag.skymap_tools import \
#             (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
from package.src.PyBlastAfterglowMag.utils import \
    (latex_float, cgs, get_beta, get_Gamma,get_Beta,BetFromMom,GamFromMom,MomFromGam)
from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
from package.src.PyBlastAfterglowMag.id_maker_from_kenta_bns import \
    prepare_kn_ej_id_2d, load_init_data
from package.src.PyBlastAfterglowMag.skymap_tools import \
    (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)

def plot_init_profile(ctheta, betas, eks,
                      xmin=0,xmax=90,ymin=1e-2,ymax=6,vmin=1e-12,vmax=1e-6,
                      norm_mode="log", cmap = plt.get_cmap('RdYlBu_r'),
                      xscale="linear",yscale="linear",
                      subplot_mode="sum",
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

    if norm_mode=="log":
        ax1.set_yscale("log")
    else:
        ax1.set_yscale("linear")
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
        if (subplot_mode=="sum"): _yarr = np.sum(ctheta*eks[mask, :],axis=0)
        elif (subplot_mode=="ave"): _yarr = axx/ayy
        else:raise KeyError("subplot_mode is not recognized")
        ax11.plot(ctheta, np.sum(eks[mask, :],axis=0), color="gray", ls="-")
    if norm_mode=="log":
        ax11.set_yscale("log")
    else:
        ax11.set_yscale("linear")
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

    if (norm_mode=="log"):
        norm = LogNorm(eks[(eks > 0) & (np.isfinite(eks))].min(), eks[(eks > 0) & (np.isfinite(eks))].max())
    elif (norm_mode=="linear"):
        norm = Normalize(eks[(eks > 0) & (np.isfinite(eks))].min(), eks[(eks > 0) & (np.isfinite(eks))].max())
    elif (norm_mode=="levels"):
        levels = MaxNLocator(nbins=15).tick_values(eks.min(), eks.max())
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    else:
        raise KeyError(" norm_mode is not recognized ")


    im = ax0.pcolor(ctheta, moms, eks, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

    ax0.axhline(y=1, linestyle='--', color='gray')

    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)
    ax0.set_xscale(xscale)
    ax0.set_yscale(yscale)

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

def plot_skymaps_():
    workdir = os.getcwd()+'/'
    pba = PyBlastAfterglow(workingdir=workdir,readparfileforpaths=True,parfile="parfile.par")

    plt.loglog(pba.PWN.get_skymap_times(),pba.PWN.get_skymap_totfluxes(freq=pba.PWN.get_skymap_freqs()[0]))
    plt.show()

    print(f"times={pba.PWN.get_skymap_times()}")
    print(f"freqs={pba.PWN.get_skymap_freqs()}")
    dfile = pba.PWN.get_skymap_obj()#(time=pba.PWN.get_skymap_times()[0],freq=pba.PWN.get_skymap_freqs()[0])
    ddfile = dfile["time={:.4e} freq={:.4e}".format(pba.PWN.get_skymap_times()[0], pba.PWN.get_skymap_freqs()[0])]
    ishell=0
    d_l = float(dfile.attrs["d_l"])
    r_i = np.array(ddfile["r"][ishell])
    mu_i = np.array(ddfile["mu"][ishell])
    ctheta = np.array(ddfile["ctheta"][ishell])
    phi = np.array(ddfile["phi"][ishell])
    tau_compton_i = np.array(ddfile["tau_compton"][ishell])
    tau_bh_i = np.array(ddfile["tau_bh"][ishell])
    tau_bf_i = np.array(ddfile["tau_bf"][ishell])
    xrs_i = np.array(ddfile["xrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
    yrs_i = np.array(ddfile["yrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
    int_i = np.array(ddfile["intensity"][ishell]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2

    x1 = r_i * np.sin(phi) * np.cos(ctheta)
    y1 = r_i * np.sin(phi) * np.sin(ctheta)
    z1 = r_i * np.cos(phi)

    fig = plt.figure()
    cmap = cm.get_cmap('viridis')
    my_norm = LogNorm(int_i.max() * 1e-2, int_i.max())
    ax = fig.add_subplot(projection='3d')
    # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
    # ax.scatter(np.log10(x1).flatten(), np.log10(y1).flatten(), np.log10(z1).flatten(), c=cmap(my_norm(int_i.flatten())))
    ax.scatter(ctheta.flatten(), phi.flatten(), np.log10(r_i).flatten(), c=cmap(my_norm(int_i.flatten())))
    # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('I Label')

    ax = fig.add_subplot()
    ax.loglog(pba.PWN.get_skymap_times(),pba.PWN.get_skymap_totfluxes(freq=pba.PWN.get_skymap_freqs()[0]))
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.show()





    fig = plt.figure()
    cmap = cm.get_cmap('viridis')
    my_norm = LogNorm(int_i.max() * 1e-2, int_i.max())
    ax = fig.add_subplot(projection='3d')
    # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
    ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(r_i).flatten(), c=cmap(my_norm(int_i.flatten())))
    # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('I Label')
    plt.show()

def plot_spectrum():
    workdir = os.getcwd()+'/'
    pba = PyBlastAfterglow(workingdir=workdir,readparfileforpaths=True,parfile="parfile.par")

    # data = np.vstack((
    #     [pba.PWN.get_lc(freq=freq,ishell= None,ilayer=None,spec=True)[1:-1] for freq in pba.PWN.get_lc_freqs(spec=True)]
    # ))
    data = pba.PWN.get_lc(ishell=0,ilayer=0,spec=True);#/np.reshape(pba.PWN.get_lc_freqs(spec=True),(-1,1))
    x_arr = pba.PWN.get_lc_times(spec=True)[1:-1]
    y_arr = pba.PWN.get_lc_freqs(spec=True)
    min_val = data.max()/1e10#//data[data>0].min()
    max_val = data.max()
    levels = np.logspace(np.log10(min_val), np.log10(max_val),  20)
    fig, axes = plt.subplots(nrows=1,ncols=1)
    ax = axes
    # cs = ax.contourf(pba.PWN.get_lc_times(spec=True)[1:-1],
    #                   pba.PWN.get_lc_freqs(spec=True),
    #                   data,
    #                   norm=LogNorm(vmin=min_val, vmax=max_val),
    #                   locator = ticker.LogLocator(),
    #                   antialiased=True, alpha=1., linewidths=0.4)
    cf = ax.contourf(x_arr/cgs.day, y_arr, data, levels=levels,
                      cmap='viridis', norm=LogNorm(vmin=min_val, vmax=max_val),
                      #locator=ticker.LogLocator(),
                      antialiased=True, alpha=1., linewidths=0.4)  # ,vmin=1e10,vmax=1e30)
    ax2_divider = make_axes_locatable(ax)
    cax2 = ax2_divider.append_axes("top", size="7%", pad="2%")

    ticks = levels[::4]
    cbar = fig.colorbar(cf, #location='top',
                        # ax=ax1,
                        cax=cax2,#pad=0.1
                        orientation="horizontal",
                        # format='%.1e',
                        ticks=ticks
                        )
    cbar.ax.set_xticklabels([r"$10^{" + "{:.0f}".format(i) + "}$" for i in np.log10(ticks)])  # add the labels
    cax2.xaxis.set_ticks_position("top")
    cax2.set_title(r"Intensity $I_{\nu'}'$ [erg $s^{-1}$ cm$^{-2}$ Hz$^{-1}$]", fontsize=12)
    cax2.tick_params(direction='in',labelsize=12)
    ax.set_yscale("log")
    ax.set_xlim(x_arr[0]/cgs.day, x_arr[-1]/cgs.day)
    ax.set_ylim(1e9, 1e22)
    ax.set_xscale("log")
    ax.set_ylabel(r"Frequency $\nu'$ [Hz]", fontsize=12)
    ax.set_xlabel(r"Time, days", fontsize=12)
    ax.get_yaxis().set_label_coords(-0.15, 0.5)

    # ax.xaxis.set_ticklabels([])
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=False,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    plt.savefig(os.getcwd()+"/old_pwn_spec2_nossa.png")
    plt.show()
    # ax = plt.colorbar(cs)

    # plt.title('matplotlib.pyplot.contourf() Example')
    # plt.show()

    print(f"times {pba.PWN.get_lc_times(spec=True).shape}")
    print(f"freqs {pba.PWN.get_lc_freqs(spec=True).shape}")
    print(data.shape)

    print( pba.PWN.get_lc(freq=3e9,ishell=0,ilayer=0,spec=True) )

    # print(pba.PWN.get_lc_times(unique=True,spec=True))
    # print(pba.PWN.get_lc_freqs(unique=True,spec=True))

    print(pba.PWN.get_lc(freq=3.e9,ishell=0,ilayer=0,spec=True)/pba.PWN.get_lc(freq=2.4e18,ishell=0,ilayer=0,spec=True))

    # plt.loglog(pba.PWN.get_lc_times(spec=True,unique=True)[1:-1],pba.PWN.get_lc(freq=3.e9,ishell=0,ilayer=0,spec=True)[1:-1])#/pba.PWN.get_lc_obj(spec=True).attrs["d_l"]**2)
    # plt.loglog(pba.PWN.get_lc_times(spec=True,unique=True)[1:-1],pba.PWN.get_lc(freq=2.4e18,ishell=0,ilayer=0,spec=True)[1:-1])#/pba.PWN.get_lc_obj(spec=True).attrs["d_l"]**2)
    # plt.loglog(pba.PWN.get_lc_times(spec=True,unique=True)[1:-1],pba.PWN.get_lc(freq=2.4e22,ishell=0,ilayer=0,spec=True)[1:-1])#/pba.PWN.get_lc_obj(spec=True).attrs["d_l"]**2)
    plt.loglog(pba.PWN.get_lc_freqs(spec=True,unique=True),pba.PWN.get_lc(time=5e2,ishell=0,ilayer=0,spec=True))


    # plt.savefig(os.getcwd()+"/old_pwn_spec.png")
    plt.show()

def plot_lc():
    workdir = os.getcwd()+'/'
    pba = PyBlastAfterglow(workingdir=workdir,readparfileforpaths=True,parfile="parfile.par")
    freqs = pba.PWN.get_lc_freqs()
    times = pba.PWN.get_lc_times()
    cmap = cm.Reds_r
    nlayers = int(pba.PWN.get_lc_obj().attrs["nlayers"])
    nshells = 1
    mynorm = Normalize(vmin=0,vmax=nlayers*nshells)#norm(len(ishells)*len(ilayers))
    print(pba.PWN.get_lc(freq=freqs[0],ishell=0,ilayer=0))
    fig, axes = plt.subplots(ncols=1,nrows=3)
    for il, layer in enumerate(range(nlayers)):
        color=cmap(mynorm(int(il)))
        axes[0].plot(times,pba.PWN.get_lc(freq=freqs[0],ishell=0,ilayer=il),color=color,ls='-')
        axes[0].plot(times,pba.PWN.get_lc(freq=freqs[1],ishell=0,ilayer=il),color=color,ls='--')
        # axes[1].plot(pba.PWN.get_dyn_arr("tburst",ishell=0,ilayer=il),
        #              pba.PWN.get_dyn_arr("Rw",ishell=0,ilayer=il),color=color,ls='-')
        axes[1].plot(pba.PWN.get_dyn_arr("tburst",ishell=0,ilayer=il),
                     pba.PWN.get_dyn_arr("Gamma",ishell=0,ilayer=il),color=color,ls='-')
        axes[2].plot(pba.KN.get_dyn_arr("tburst",ishell=0,ilayer=il),
                     pba.KN.get_dyn_arr("Gamma",ishell=0,ilayer=il),color=color,ls='-')

    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xscale("log")

    plt.legend()
    plt.show()

def run():
    workdir = os.getcwd()+'/'
    path_to_original_data = "/media/vsevolod/data/KentaData/SFHo_13_14_150m_11/" #
    dfile = h5py.File(workdir+"kenta_ejecta_13.h5")
    print(dfile.keys())

    text = 70.
    label = f"corr_id_SFHo_13_14_150m_11_text{int(text)}"
    sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    files = sorted(glob(workdir + "kenta_ejecta*", recursive=True), key=sort_by)

    prepare_kn_ej_id_2d(files=files,
                        outfpaths=[workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5"],
                        req_times=np.array([text]),
                        new_theta_len=20,
                        new_vinf_len=None,
                        verbose=True,
                        r0type="fromrho", r0frac=0.5, t0=-1,
                        dist="pw")

    r_, mom_, theta_, ctheta_, ek_, mass_, ye_, s_, rho_, temp_ \
        = load_init_data(workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5")

    plot_init_profile(ctheta_[0,:], mom_[:,0], rho_,
                      figpath=None,#FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
                      norm_mode="log",
                      subplot_mode="ave",
                      title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(text)))

    modify_parfile_par_opt(workingdir=workdir, part="kn", newpars={},
                           newopts={"fname_light_curve":f"{label}.lc",
                                    "fname_ejecta_id":f"{label}.h5"},
                           parfile="parfile.par",newparfile="parfile.par",keep_old=False )

    pba = PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True)
    pba.reload_parfile()
    pba.run(loglevel="info")


def main():
    # os.getcwd()
    # run()
    plot_lc()
    # plot_spectrum()
    # plot_skymaps_()

if __name__ == '__main__':
    main()