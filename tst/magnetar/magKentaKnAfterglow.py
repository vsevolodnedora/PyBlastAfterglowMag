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

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import PyBlastAfterglow
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import \
        (latex_float, cgs, get_beta, get_Gamma,get_Beta,BetFromMom,GamFromMom,MomFromGam)
    from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
    from PyBlastAfterglowMag.id_maker_from_kenta_bns import \
        prepare_kn_ej_id_2d, load_init_data
    from PyBlastAfterglowMag.skymap_tools import \
        (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
except ImportError:
    try:
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
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

# import matplotlib.pyplot as plt
# N = 1000
# y = np.zeros(N)
# plt.semilogx(np.geomspace(1, 1000, N, endpoint=True), y + 1, 'o')
# plt.semilogx(np.geomspace(1, 1000, N, endpoint=False), y + 2, 'o')
# plt.semilogx(np.logspace(1, 1000, N, endpoint=False), y + 3, 'o')
# plt.axis([0.5, 2000, 0, 4])
# plt.grid(True, color='0.7', linestyle='-', which='both', axis='both')
# plt.show()



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



def main_old():
    workdir = os.getcwd()+'/'
    ### locate and extract the data from the original ejecta profiles from Kenta
    path_to_original_data = "/media/vsevolod/data/KentaData/SFHo_13_14_150m_11/" #
    # if not os.path.isdir(path_to_original_data):
    #     raise FileExistsError("Input ejecta data is not found")

    dfile = h5py.File(workdir+"kenta_ejecta_13.h5")
    print(dfile.keys())

    # print(np.gradient([0.1,0.3], edge_order=1))

    # ## lead the original data and prepare the ID file for PyBlastAfterglow
    # text = 25.
    # label = f"corr_id_SFHo_13_14_150m_11_text{int(text)}"
    # sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    # files = sorted(glob("/media/vsevolod/data/KentaData/SFHo_13_14_150m_11/" + "ejecta*", recursive=True), key=sort_by)

    text = 70.
    label = f"corr_id_SFHo_13_14_150m_11_text{int(text)}"
    sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    files = sorted(glob(workdir + "kenta_ejecta*", recursive=True), key=sort_by)

    prepare_kn_ej_id_2d(files=files,
                        outfpaths=[workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5"],
                        req_times=np.array([text]),
                        new_theta_len=None,
                        new_vinf_len=None,
                        verbose=True,
                        r0type="fromrho",
                        dist="pw")



    #
    mom_, theta_, ctheta_, ek_, mass_, ye_, s_, rho_, temp_\
        = load_init_data(workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5")
    r = np.zeros_like(mass_)
    for ith in range(len(ctheta_[0,:])):
        idx = 0
        k = 0.5
        r[idx,ith] = (k*(3/4./np.pi)*rho_[idx,ith]*mass_[idx,ith])**(1./3.)

        if (r[idx,ith] == 0):
            raise ValueError()
        for ir in range(1,len(mom_[:,0]),1):
            _val = (3./4./np.pi)*rho_[ir,ith]*mass_[ir,ith]
            _rm1 = r[ir-1,ith]**3
            if(mass_[ir,ith]>0):
                r[ir,ith] = (_rm1 + _val)**(1./3.)

    # t = 70. * 1e-3
    # for ith in range(len(ctheta_[0,:])):
    #     for ir in range(0,len(mom_[:,0]),1):
    #         r[ir,ith] =  BetFromMom(mom_[ir,ith])*cgs.c * t


    # r = np.zeros_like(mass_)
    # for ith in range(len(ctheta_)):
    #     idx = np.argmax(np.where(mass_[:,ith] > 0))
    #     if idx != len(mom_)-1: print(idx, mass_[idx,ith], mass_[idx+1,ith])
    #     if (len(mass_[:,ith]>0.)==1):
    #         raise ValueError()
    #
    #     k = 100.
    #     r[idx,ith] = (k*(3/4./np.pi)*rho_[idx,ith]*mass_[idx,ith])**(1./3.)
    #     for i in range(idx-1, 0, -1):
    #         _rp1 = r[i+1,ith]**3
    #         _val = (3./4./np.pi)*rho_[i,ith]*mass_[i,ith]
    #         if (_rp1 < _val):
    #             raise ValueError()
    #         r[i,ith] = (_rp1 - _val)**(1./3.)
    #         y = 1
    #     print(r[:,ith])
    #     x = 1

        # if (idx == 1):
        #     print(mass_[:,ith])

    # exit(1)
    # k = 1./2.
    # rn = (k*(3/4/np.pi)*rho_[-1,:]*mass_[-1,:])**(1./3.)
    # r = np.zeros_like(mass_)
    # r[-1,:] = rn
    # for i in range(len(mom_)-2, 0, -1):
    #     print(i)
    #     r[i,:] = (r[i+1,:]**3 - (3./4./np.pi)*rho_[i,:]*mass_[i,:])**(1./3.)




    plot_init_profile(ctheta_[0,:], mom_[:,0], r,
                      figpath=None,#FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
                      norm_mode="log",
                      subplot_mode="ave",
                      title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(text)))

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

    # dfile = h5py.File(label+".h5",'a')
    # dfile = h5py.File(label+".h5",'a')
    # print(np.array(dfile["mom"]))
    # dfile.create_dataset(name="mom", data=np.array(dfile["vel_inf"])*get_Gamma(np.array(dfile["vel_inf"])))
    # dfile.create_dataset(name="mom", data=np.array(dfile["vel_inf"])*get_Gamma(np.array(dfile["vel_inf"])))
    # dfile.create_dataset(name="s", data=np.full_like(np.array(dfile["ek"]),fill_value=10.))
    # dfile.close()


    pba = PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True)
    pba.reload_parfile()
    pba.run(loglevel="info")



    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    # ax = axes
    # ax.loglog(pba.get_ej_lc_times() / cgs.day, pba.get_ej_lc_totalflux(freq=3.e9),
    #         **{"color":"black", "ls": "--", "lw": 0.8, "label": label})
    #
    # pba.clear()
    # plt.show()

def main():

    # x_arr = np.array([ 6.62534e-19, 9.57652e-19, 1.38423e-18, 2.00082e-18, 2.89206e-18, 4.1803e-18, 6.04238e-18, 8.73389e-18, 1.26243e-17, 1.82477e-17, 2.63759e-17, 3.81248e-17, 5.51071e-17, 7.96541e-17, 1.15135e-16, 1.66421e-16, 2.40551e-16, 3.47703e-16, 5.02583e-16, 7.26454e-16, 1.05005e-15, 1.51778e-15, 2.19386e-15, 3.17109e-15, 4.58361e-15, 6.62534e-15, 9.57652e-15, 1.38423e-14, 2.00082e-14, 2.89206e-14, 4.1803e-14, 6.04238e-14, 8.73389e-14, 1.26243e-13, 1.82477e-13, 2.63759e-13, 3.81248e-13, 5.51071e-13, 7.96541e-13, 1.15135e-12, 1.66421e-12, 2.40551e-12, 3.47703e-12, 5.02583e-12, 7.26454e-12, 1.05005e-11, 1.51778e-11, 2.19386e-11, 3.17109e-11, 4.58361e-11, 6.62534e-11, 9.57652e-11, 1.38423e-10, 2.00082e-10, 2.89206e-10, 4.1803e-10, 6.04238e-10, 8.73389e-10, 1.26243e-09, 1.82477e-09, 2.63759e-09, 3.81248e-09, 5.51071e-09, 7.96541e-09, 1.15135e-08, 1.66421e-08, 2.40551e-08, 3.47703e-08, 5.02583e-08, 7.26454e-08, 1.05005e-07, 1.51778e-07, 2.19386e-07, 3.17109e-07, 4.58361e-07, 6.62534e-07, 9.57652e-07, 1.38423e-06, 2.00082e-06, 2.89206e-06, 4.1803e-06, 6.04238e-06, 8.73389e-06, 1.26243e-05, 1.82477e-05, 2.63759e-05, 3.81248e-05, 5.51071e-05, 7.96541e-05, 0.000115135, 0.000166421, 0.000240551, 0.000347703, 0.000502583, 0.000726454, 0.00105005, 0.00151778, 0.00219386, 0.00317109, 0.00458361])
    # y_arr = np.array([ 9.89229e+10, 8.22805e+10, 6.8438e+10, 5.69242e+10, 4.73475e+10, 3.93819e+10, 3.27565e+10, 2.72456e+10, 2.26619e+10, 1.88494e+10, 1.56782e+10, 1.30406e+10, 1.08467e+10, 9.02188e+09, 7.50407e+09, 6.24162e+09, 5.19155e+09, 4.31814e+09, 3.59167e+09, 2.98743e+09, 2.48483e+09, 2.06679e+09, 1.71908e+09, 1.42987e+09, 1.18932e+09, 9.89229e+08, 8.22805e+08, 6.8438e+08, 5.69242e+08, 4.73475e+08, 3.93819e+08, 3.27565e+08, 2.72456e+08, 2.26619e+08, 1.88494e+08, 1.56782e+08, 1.30406e+08, 1.08467e+08, 9.02188e+07, 7.50407e+07, 6.24162e+07, 5.19155e+07, 4.31814e+07, 3.59167e+07, 2.98743e+07, 2.48483e+07, 2.06679e+07, 1.71908e+07, 1.42987e+07, 1.18932e+07, 9.89229e+06, 8.22805e+06, 6.8438e+06, 5.69242e+06, 4.73475e+06, 3.93819e+06, 3.27565e+06, 2.72456e+06, 2.26619e+06, 1.88494e+06, 1.56782e+06, 1.30406e+06, 1.08467e+06, 902188, 750407, 624162, 519155, 431814, 359167, 298743, 248483, 206679, 171908, 142987, 118932, 98922.9, 82280.5, 68438, 56924.2, 47347.5, 39381.9, 32756.5, 27245.6, 22661.9, 18849.4, 15678.2, 13040.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    #
    # plt.loglog(x_arr, y_arr, ls='--', label='afg')
    # plt.show()


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
                        new_theta_len=4,
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


if __name__ == '__main__':
    main()