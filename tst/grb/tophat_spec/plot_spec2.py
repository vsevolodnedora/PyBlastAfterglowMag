# import PyBlastAfterglowMag
import numpy as np
import h5py
from pyrsistent import pbag
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm, SymLogNorm, TwoSlopeNorm, CenteredNorm
from matplotlib import cm
import os
import copy
import gc
# try:
#     from PyBlastAfterglowMag.interface import modify_parfile_par_opt
#     from PyBlastAfterglowMag.interface import PyBlastAfterglow
#     from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#     from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
#     from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

# try:
#     import package.src.PyBlastAfterglowMag as PBA
# except ImportError:
#     try:
#         import PyBlastAfterglowMag as PBA
#     except:
#         raise ImportError("Cannot import PyBlastAfterglowMag")
import package.src.PyBlastAfterglowMag as PBA

curdir = os.getcwd()+'/'

def run_2_pba(i_thetaobs = 0., run=True, rs_dyn="yes", rs_rad="yes"):
    th = str(i_thetaobs).replace(".","")

    # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
    #                                          parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
    #                                          newpars={"gam1":1,"gam2":1e8,"ngam":250},
    #                                          newopts={"method_synchrotron_fs":"Dermer09",
    #                                                   "method_synchrotron_rs":"Dermer09",
    #                                                   "method_comp_mode":"comovSpec",
    #                                                   "method_eats":"adaptive",
    #                                                   "method_ne_fs":"useNe",
    #                                                   "method_ne_rs":"useNe",
    #                                                   "method_gamma_c_fs" : "useTcomov",
    #                                                   "method_gamma_c_rs" : "useTcomov",
    #                                                   "method_gamma_min_fs":"useNumericGamma",
    #                                                   "method_gamma_min_rs":"useNumericGamma",
    #                                                   "method_nonrel_dist_fs":"use_gamma_min",
    #                                                   "method_nonrel_dist_rs":"use_gamma_min",
    #                                                   "limit_lf_min_to1_fs":"yes",
    #                                                   "limit_lf_min_to1_rs":"yes",
    #                                                   "method_ele_fs":"analytic",
    #                                                   "method_ele_rs":"analytic",
    #                                                   "do_rs" : rs_dyn,
    #                                                   "do_rs_radiation": rs_rad,
    #                                                   "bw_type":"fsrs" if rs_dyn=="yes" else "fs",
    #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
    #                                                   "fname_dyn":f"dyn_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5",
    #                                                   "fname_light_curve":f"tophat_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5",
    #                                                   "fname_spectrum":f"spectrum_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5"},
    #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    # pba_an = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
    # if run: pba_an.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")

    # ----------------------------------------

    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_synchrotron_rs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "method_ne_rs":"useNe",
                                                      "method_gamma_c_fs" : "useTcomov",
                                                      "method_gamma_c_rs" : "useTcomov",
                                                      "method_gamma_min_fs":"useNumericGamma",
                                                      "method_gamma_min_rs":"useNumericGamma",
                                                      "method_nonrel_dist_fs":"use_gamma_min",
                                                      "method_nonrel_dist_rs":"use_gamma_min",
                                                      "limit_lf_min_to1_fs":"yes",
                                                      "limit_lf_min_to1_rs":"yes",
                                                      "method_ele_fs":"mix",
                                                      "method_ele_rs":"mix",
                                                      "do_rs" : rs_dyn,
                                                      "do_rs_radiation": rs_rad,
                                                      "bw_type":"fsrs" if rs_dyn=="yes" else "fs",
                                                      "fname_ejecta_id":"tophat_grb_id_a.h5",
                                                      "fname_dyn":f"dyn_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5",
                                                      "fname_light_curve":f"tophat_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5",
                                                      "fname_spectrum":f"spectrum_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_an.h5"},
                                             parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    pba_mix = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
    if run: pba_mix.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")

    # ----------------------------------------

    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_synchrotron_rs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "method_gamma_c_fs" : "useTcomov",
                                                      "method_gamma_c_rs" : "useTcomov",
                                                      "method_gamma_min_fs":"useNumericGamma",
                                                      "method_gamma_min_rs":"useNumericGamma",
                                                      "method_nonrel_dist_fs":"use_gamma_min",
                                                      "method_nonrel_dist_rs":"use_gamma_min",
                                                      "limit_lf_min_to1_fs":"yes",
                                                      "limit_lf_min_to1_rs":"yes",
                                                      "method_ele_fs":"numeric",
                                                      "method_ne_rs":"useNe",
                                                      "method_ele_rs":"numeric",
                                                      "do_rs" : rs_dyn,
                                                      "do_rs_radiation": rs_rad,
                                                      "bw_type":"fsrs" if rs_dyn=="yes" else "fs",
                                                      "fname_ejecta_id":"tophat_grb_id_a.h5",
                                                      "fname_dyn":f"dyn_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_num.h5",
                                                      "fname_light_curve":f"tophat_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_num.h5",
                                                      "fname_spectrum":f"spectrum_{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_num.h5"},
                                             parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    pba_num = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
    if run: pba_num.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
    # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
    #         pba_a2.GRB.get_lc_totalflux(freq=i_freq), color=i_color, ls='-', lw=1,
    #         label=r" $\nu$={:.1e} Ne numeric".format(i_freq))
    # pba_a2.clear()
    return (pba_mix, pba_num)

def plot_two_spectra(ej:PBA.Ejecta, ej_an:PBA.Ejecta, fs_or_rs :str = "fs", ele_syn_ssc:str = "ele",
                     xkey="times_gams", ykey="gams",
                     task={"ys":(1e9,1e18,1e22), "colors_freqs":("blue","green","red"), "ylim":(1,1e9),
                           "xs":(1.e3,1.e5,1.e7), "colors_freqs":("blue","green","red"), "xlim":(2e3, 1e8),
                           "figname":"spec_fs_ele.png"}):
    mp = 1.6726e-24
    qe = 4.803204e-10
    me = 9.1094e-28
    c  = 2.99792458e10

    spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs,xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,spec=True)
    spec_an=ej_an.get_lc(key=ele_syn_ssc+'_'+fs_or_rs,xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,spec=True)


    m2 = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3",ishell=0,ilayer=0)
    m2_an = ej_an.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3",ishell=0,ilayer=0)
    dm = np.diff(m2)
    dm_an = np.diff(m2_an)

    xs = ej.get_grid(key=xkey,spec=True)
    ys = ej.get_grid(key=ykey,spec=True)

    xs_an = ej_an.get_grid(key=xkey,spec=True)
    ys_an = ej_an.get_grid(key=ykey,spec=True)

    n_ele = np.trapz(y=spec, x=ys, axis=1)#[:,np.newaxis]
    n_ele_an = np.trapz(y=spec_an, x=ys_an, axis=1)#[:,np.newaxis]

    fig,ax = plt.subplots(ncols=1,nrows=1)

    ax.plot(xs[:-1],dm/mp/n_ele[:-1])
    ax.plot(xs[:-1],dm_an/mp/n_ele_an[:-1])
    # ax.plot(xs,ej.get_dyn_arr(v_n="gamma_c",ishell=0,ilayer=0))
    # ax.plot(xs,ej.get_dyn_arr(v_n="M2",ishell=0,ilayer=0))
    # ax.plot(xs[:-1],ej.get_dyn_arr(v_n="M2",ishell=0,ilayer=0))
    # ax.plot(xs_an,ej_an.get_dyn_arr(v_n="gamma_c",ishell=0,ilayer=0))
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.show()


    if ele_syn_ssc == "n_ele":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        spec *= np.power(ys, 2)[np.newaxis,:]

        spec_an = spec_an / np.trapz(y=spec_an, x=ys_an, axis=1)[:,np.newaxis]
        spec_an *= np.power(ys_an, 2)[np.newaxis,:]
    else:
        spec *= ys[np.newaxis,:]
        spec_an *= ys_an[np.newaxis,:]

    spec[~np.isfinite(spec)] = 1e-100
    spec[spec <= 0] = 1e-100
    spec_an[~np.isfinite(spec_an)] = 1e-100
    spec_an[spec_an <= 0] = 1e-100

    # spec = spec.T
    # spec_an = spec_an.T

    fig, axes = plt.subplots(ncols=2, nrows=4, figsize=(8, 6), #sharey="row",sharex="col",
                             gridspec_kw=dict(width_ratios=[0.5,1.0],hspace=-0.1,wspace=0.01),#dict(height_ratios=[0.5,1.2]),
                             layout='constrained')

    # axes[0,0].axis("off")

    # -----------| Time Evolution Of Spectra At A Given Freq/Gamma |----------
    # X : inverted Flux/Number Y : Freq/Gams
    ax = axes[0,1]
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])

    ax = ax.twinx()
    colors= task["colors_ys"]
    indexes_ys = [PBA.utils.find_nearest_index(ys,freq) for freq in task["ys"]]
    lines = []
    for idx, color in zip(indexes_ys, colors):
        vals = spec[:, idx]
        lns, = ax.plot(xs, vals, color=color, linewidth=1.0, linestyle="-")
        lines.append(lns)
        vals_an = spec_an[:, idx]
        ax.plot(xs_an, vals_an, color=color, linewidth=1.0, linestyle="--")

    legend1 = ax.legend(lines, [task["ylabel"]+" = {:.2e}".format(ye) for ye in task["ys"]], loc=1)
    lin_, = ax.plot([0,0],[1,1],color='gray',ls='-')
    lin2_, = ax.plot([0,0],[1,1],color='gray',ls='--')
    legend2 = ax.legend([lin_,lin2_], ["numeric","analytic"], loc=4)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.set_xscale('log')
    ax.set_yscale('log')

    if not task["ylim"]:
        ax.set_ylim(np.max(spec[:, indexes_ys])*1e-5,
                    np.max(spec[:, indexes_ys])*10)
    else:
        ax.set_ylim(ys[0], ys[-1])
    ax.set_xlim(xs[0],xs[-1])
    ax.set_ylabel(task["zlabel"])
    ax.minorticks_on()
    ax.set_xticks([])

    # -------------------| Spectral Subplot |------------------------------

    ax = axes[1,0]
    colors= task["colors_xs"]
    indexes_xs = [PBA.utils.find_nearest_index(xs,time) for time in task["xs"]]
    for idx, color in zip(indexes_xs, colors):
        vals = spec[idx, :]
        ax.plot(vals, ys, color=color, linewidth=1.0, linestyle="-")

        vals_an = spec_an[idx, :]
        ax.plot(vals_an, ys_an, color=color, linewidth=1.0, linestyle="--")

    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim(freqs[0], freqs[-1])
    ax.set_xlabel(task["zlabel"])#(r'$\nu j_{\nu}$ [erg/s]')
    ax.set_ylabel(task["ylabel"])
    ax.legend(loc='upper right', prop={'size': 9})  # , bbox_to_anchor=(1.1, 1.1)
    if not task["xlim"]:
        ax.set_xlim(np.max(spec[indexes_xs, :])*10,
                    np.max(spec[indexes_xs, :])*1e-7)
    else:
        ax.set_xlim(*task["xlim"])
    ax.set_ylim(ys[0],ys[-1])
    ax.minorticks_on()

    # --------------------------------------------------------------------------------
    # X : time Y : Gams/Freqs
    ax = axes[1,1]
    for idx, color in zip(indexes_xs, task["colors_xs"]):
        ax.plot([xs[idx], xs[idx]], [ys[0],ys[15]],color=color, linewidth=2, linestyle="-")
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    # ax = axes[2, 0]  # [1, 1]  # row,col
    # spec_syn = spec_syn * freqs[np.newaxis, :]
    norm = LogNorm(vmin=spec.max() * 1e-6, vmax=spec.max() * 10)
    _c = ax.pcolormesh(xs, ys, spec.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(task["ylabel"])
    ax.minorticks_on()
    ax.set_xlim(xs[0], xs[-1])

    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)

    if ele_syn_ssc == "n_ele":
        ax.plot(xs,gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs,gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs,gamma_max, color='black', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min * gamma_min * gamToNuFactor
        nu_c = XpS * gamma_c * gamma_c * gamToNuFactor
        nu_max = XpS * gamma_max * gamma_max * gamToNuFactor
        ax.plot(xs,nu_min, color='black', linewidth=1, linestyle=":", label=r"$\nu_{\rm m}$")
        ax.plot(xs,nu_c, color='black', linewidth=1, linestyle="--", label=r"$\nu_{\rm c}$")
        ax.plot(xs,nu_max, color='black', linewidth=2, linestyle="-.", label=r"$\nu_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])

    # else:
    #     ax2 = ax.twinx()
    #     ax2.set_yscale("log")
    #     ax2.set_zorder(0)
    #     ax2.grid(False)
    #     ax2.tick_params(axis = 'y')
    #     ax2.plot(xs,gamma_min, color='orange', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
    #     ax2.plot(xs,gamma_c, color='orange', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
    #     ax2.set_ylim(ys_an[0],ys_an[-1])
    # ax.set_xticks([])
    # ax.set_yticks([])

    # if fs_or_rs == "rs":
    #     Gamma_rs = ej.get_dyn_arr(v_n="GammaRsh",ishell=0,ilayer=0)
    #     idx = np.argmax(np.where(Gamma_rs != 0.))
    #     ax.plot([xs[idx], xs[idx]],[ys[0],ys[10]],color="gray", linewidth=2, linestyle="-")
    #     gamma_min_rs = ej.get_dyn_arr(v_n="gamma_min_rs",ishell=0,ilayer=0)
    #     gamma_c_rs = ej.get_dyn_arr(v_n="gamma_c_rs",ishell=0,ilayer=0)
    #     # ax2 = ax.twinx()
    #     ax2.plot(xs, gamma_min_rs, color='magenta', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m;\, rs}$")
    #     ax2.plot(xs, gamma_c_rs, color='magenta', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c;\, rs}$")
    #     # ax2.set_ylim(1e-2, 1e9)

    # ---------------------------------------------------------------

    # axes[2,0].axis("off")

    # ---------------------------------------------------------------
    # X : time Y : Gams/Freqs
    ax = axes[2,1]
    for idx, color in zip(indexes_xs, task["colors_xs"]):
        ax.plot([xs_an[idx], xs_an[idx]], [ys_an[0],ys_an[15]],color=color, linewidth=2, linestyle="-")
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs_an[0], xs_an[15]], [ys_an[idx],ys_an[idx]],color=color, linewidth=2, linestyle="-")
    # ax = axes[2, 0]  # [1, 1]  # row,col
    # spec_syn = spec_syn * freqs[np.newaxis, :]
    # norm = LogNorm(vmin=spec_an.max() * 1e-6, vmax=spec_an.max() * 10)
    _c = ax.pcolormesh(xs_an, ys_an, spec_an.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(task["xlabel"])
    ax.set_ylabel(task["ylabel"])
    # ax.set_xlabel(r'$\nu$ [Hz]')
    ax.set_xlim(xs_an[0], xs_an[-1])
    ax.minorticks_on()
    gamma_min_an = ej_an.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c_an = ej_an.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max_an = ej_an.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B_an = ej_an.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)

    if ele_syn_ssc == "n_ele":
        ax.plot(xs_an,gamma_min_an, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs_an,gamma_c_an, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs_an,gamma_max_an, color='black', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B_an) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min_an * gamma_min_an * gamToNuFactor
        nu_c = XpS * gamma_c_an * gamma_c_an * gamToNuFactor
        nu_max = XpS * gamma_max_an * gamma_max_an * gamToNuFactor
        ax.plot(xs,nu_min, color='black', linewidth=1, linestyle=":", label=r"$\nu_{\rm m}$")
        ax.plot(xs,nu_c, color='black', linewidth=1, linestyle="--", label=r"$\nu_{\rm c}$")
        ax.plot(xs,nu_max, color='black', linewidth=2, linestyle="-.", label=r"$\nu_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])

        # ax2 = ax.twinx()
        # ax2.set_yscale("log")
        # ax2.set_zorder(0)
        # ax2.grid(False)
        # ax2.tick_params(axis = 'y')
        # ax2.plot(xs_an,gamma_min_an_an, color='orange', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        # ax2.plot(xs_an,gamma_c_an_an, color='orange', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        # ax2.set_ylim(ys_an[0],ys_an[-1])
    # ax.set_xticks([])
    # ax.set_yticks([])

    # if fs_or_rs == "rs":
    #     Gamma_rs = ej_an.get_dyn_arr(v_n="GammaRsh",ishell=0,ilayer=0)
    #     idx = np.argmax(np.where(Gamma_rs != 0.))
    #     ax.plot([xs[idx], xs[idx]],[ys[0],ys[10]],
    #             color="gray", linewidth=2, linestyle="-")
    #     gamma_min_rs = ej_an.get_dyn_arr(v_n="gamma_min_rs",ishell=0,ilayer=0)
    #     gamma_c_rs = ej_an.get_dyn_arr(v_n="gamma_c_rs",ishell=0,ilayer=0)
    #     # ax2 = ax.twinx()
    #     ax2.plot(xs_an, gamma_min_rs, color='magenta', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m;\, rs}$")
    #     ax2.plot(xs_an, gamma_c_rs, color='magenta', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c;\, rs}$")
    #     # ax2.set_ylim(1e-2, 1e9)
    cbar = fig.colorbar(_c, ax=[axes[1,1],axes[2,1]], shrink=0.95,label=task["zlabel"],pad=-0.05)

    ax = axes[3,1]
    ratio = spec_an.T / spec.T
    # norm = LogNorm(vmin=ratio.min(), vmax=ratio.max())
    norm=SymLogNorm(linthresh=1e-2, vmin=1e-2, vmax=1.e2, base=10) # linthresh=0.03, linscale=0.03,
    # norm=CenteredNorm(halfrange=1.e2, vcenter=1.) # linthresh=0.03, linscale=0.03,
    # norm=TwoSlopeNorm(vmin=1e-3, vmax=1.e3, vcenter=1.) # linthresh=0.03, linscale=0.03,
    _c = ax.pcolormesh(xs_an, ys_an, ratio, cmap='RdBu_r', norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95,label="Analytc/Numeric",pad=-0.05)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(task["xlabel"])
    ax.set_ylabel(task["ylabel"])
    ax.set_xlim(xs_an[0], xs_an[-1])
    ax.set_ylim(ys_an[0],ys_an[-1])
    # fig.colorbar(_c, ax=ax,shrink=0.9,pad=.01,label=task["zlabel"])
    # plt.subplots_adjust(wspace=-0.2, hspace=-0.2)
    # plt.tight_layout()
    plt.savefig(os.getcwd()+'/'+task["figname"])
    plt.show()

def plot_compare_two_spectra_new(ej:PBA.Ejecta, ej_an:PBA.Ejecta, name1:str= "Numeric", name2="Analytic", is_spec=True,
                                 fs_or_rs :str = "fs", ele_syn_ssc:str = "ele",
                                 xkey="times_gams", ykey="gams",
                                 task=dict()):
    mp = 1.6726e-24
    qe = 4.803204e-10
    me = 9.1094e-28
    c  = 2.99792458e10


    fig, axes = plt.subplots(ncols=1, nrows=5, figsize=(7, 8), sharex="col", #sharey="row",sharex="col",
                             gridspec_kw=dict(height_ratios=[1.0,1.,1.,1.,0.4],hspace=-0.1,wspace=0.01),#dict(height_ratios=[0.5,1.2]),
                             layout='constrained')


    spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                   xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                   spec=is_spec,sum_shells_layers=False)
    xs = ej.get_grid(key=xkey,spec=is_spec)
    ys = ej.get_grid(key=ykey,spec=is_spec)
    if ele_syn_ssc == "n_ele":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        spec *= np.power(ys, 2)[np.newaxis,:]
    else:
        spec *= ys[np.newaxis,:]
    spec[~np.isfinite(spec)] = 1e-100
    spec[spec <= 0] = 1e-100


    spec_an=ej_an.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                         xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                         spec=is_spec,sum_shells_layers=False)
    xs_an = ej_an.get_grid(key=xkey,spec=is_spec)
    ys_an = ej_an.get_grid(key=ykey,spec=is_spec)
    if ele_syn_ssc == "n_ele":
        spec_an = spec_an / np.trapz(y=spec_an, x=ys_an, axis=1)[:,np.newaxis]
        spec_an *= np.power(ys_an, 2)[np.newaxis,:]
    else:
        spec_an *= ys_an[np.newaxis,:]
    spec_an[~np.isfinite(spec_an)] = 1e-100
    spec_an[spec_an <= 0] = 1e-100

    # ------------------------------------

    ax = axes[0] #axes[1,1]
    ax.set_title(task["title"],fontsize=12)

    colors= task["colors_ys"]
    indexes_ys = [PBA.utils.find_nearest_index(ys,freq) for freq in task["ys"]]
    lines = []
    for idx, color in zip(indexes_ys, colors):
        vals = spec[:, idx]
        lns, = ax.plot(xs, vals, color=color, linewidth=1.0, linestyle="-")
        lines.append(lns)
        vals_an = spec_an[:, idx]
        ax.plot(xs_an, vals_an, color=color, linewidth=1.0, linestyle="--")

    legend1 = ax.legend(lines, [task["ylabel"]+" = {:.2e}".format(ye) for ye in task["ys"]], loc=1)
    lin_, = ax.plot([0,0],[1,1],color='gray',ls='-')
    lin2_, = ax.plot([0,0],[1,1],color='gray',ls='--')
    legend2 = ax.legend([lin_,lin2_], [name1,name2], loc=4)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    if not task["ylim"]:
        ax.set_ylim(np.max(spec[:, indexes_ys])*1e-5,
                    np.max(spec[:, indexes_ys])*10)
    else:
        ax.set_ylim(ys[0], ys[-1])
    ax.set_xlim(xs[0],xs[-1])
    ax.set_ylabel(task["zlabel"],fontsize=12)
    ax.minorticks_on()


    # ---------------------------------------------------------------

    ax = axes[1] #axes[1,1]
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.97, 0.05, name1, fontsize=12, bbox=props,
            transform=ax.transAxes, horizontalalignment='right')

    norm = LogNorm(vmin=spec.max() * 1e-6, vmax=spec.max() * 10)
    _c = ax.pcolormesh(xs, ys, spec.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(task["ylabel"],fontsize=12)
    ax.minorticks_on()
    ax.set_xlim(xs[0], xs[-1])

    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
    Gamma = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs=="fs" else "Gamma",ishell=0,ilayer=0)
    z = float(ej.get_lc_obj(spec=is_spec).attrs["z"])

    if ele_syn_ssc == "n_ele":
        ax.plot(xs,gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs,gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs,gamma_max, color='black', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys[0],ys[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min * gamma_min * gamToNuFactor
        nu_c = XpS * gamma_c * gamma_c * gamToNuFactor
        nu_max = XpS * gamma_max * gamma_max * gamToNuFactor

        # observer frame
        if is_spec:
            ax.plot(xs,nu_min, color='black', linewidth=1, linestyle=":", label=r"$\nu'_{\rm m}$" if is_spec else r"$\nu_{\rm m}$")
            ax.plot(xs,nu_c, color='black', linewidth=1, linestyle="--", label=r"$\nu'_{\rm c}$"if is_spec else r"$\nu_{\rm c}$")
            ax.plot(xs,nu_max, color='black', linewidth=2, linestyle="-.", label=r"$\nu'_{\rm M}$"if is_spec else r"$\nu_{\rm M}$")
            ax.set_ylim(ys[0],ys[-1])
        else:
            pass
            # nu_min *= (1 + z) / Gamma
            # nu_c *= (1 + z) / Gamma
            # nu_max *= (1 + z) / Gamma


    # ---------------------------------------------------------------



    ax = axes[2]#axes[2,1]
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.97, 0.05, name2, fontsize=12, bbox=props,
            transform=ax.transAxes, horizontalalignment='right')

    _c = ax.pcolormesh(xs_an, ys_an, spec_an.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(task["ylabel"],fontsize=12)

    ax.set_xlim(xs_an[0], xs_an[-1])
    ax.minorticks_on()
    gamma_min_an = ej_an.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c_an = ej_an.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max_an = ej_an.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B_an = ej_an.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
    Gamma_an = ej_an.get_dyn_arr(v_n="Gamma" if fs_or_rs=="fs" else "Gamma",ishell=0,ilayer=0)

    if ele_syn_ssc == "n_ele":
        ax.plot(xs_an,gamma_min_an, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs_an,gamma_c_an, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs_an,gamma_max_an, color='black', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B_an) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min_an * gamma_min_an * gamToNuFactor
        nu_c = XpS * gamma_c_an * gamma_c_an * gamToNuFactor
        nu_max = XpS * gamma_max_an * gamma_max_an * gamToNuFactor

        # observer frame
        if is_spec:
            ax.plot(xs,nu_min, color='black', linewidth=1, linestyle=":", label=r"$\nu'_{\rm m}$" if is_spec else r"$\nu_{\rm m}$")
            ax.plot(xs,nu_c, color='black', linewidth=1, linestyle="--", label=r"$\nu'_{\rm c}$" if is_spec else r"$\nu_{\rm c}$")
            ax.plot(xs,nu_max, color='black', linewidth=2, linestyle="-.", label=r"$\nu'_{\rm M}$" if is_spec else r"$\nu_{\rm M}$")
            ax.set_ylim(ys_an[0],ys_an[-1])
        else:
            pass
            # nu_min *= (1 + z) / Gamma
            # nu_c *= (1 + z) / Gamma
            # nu_max *= (1 + z) / Gamma


    cbar = fig.colorbar(_c, ax=[axes[1],axes[2]], shrink=0.95,label=task["zlabel"],pad=0.05)

    # ---------------------------------------------------------------

    ax = axes[3]
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.97, 0.05, f"{name2}/{name1}", fontsize=12, bbox=props,
                 transform=ax.transAxes, horizontalalignment='right')
    ratio = spec_an.T / spec.T
    # norm = LogNorm(vmin=ratio.min(), vmax=ratio.max())
    norm=SymLogNorm(linthresh=1e-2, vmin=1e-2, vmax=1.e2, base=10) # linthresh=0.03, linscale=0.03,
    # norm=CenteredNorm(halfrange=1.e2, vcenter=1.) # linthresh=0.03, linscale=0.03,
    # norm=TwoSlopeNorm(vmin=1e-3, vmax=1.e3, vcenter=1.) # linthresh=0.03, linscale=0.03,
    _c = ax.pcolormesh(xs_an, ys_an, ratio, cmap='RdBu_r', norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95,label=f"{name2}/{name1}",pad=0.05)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlabel(task["xlabel"])
    ax.set_ylabel(task["ylabel"],fontsize=12)
    ax.set_xlim(xs_an[0], xs_an[-1])
    ax.set_ylim(ys_an[0],ys_an[-1])
    # fig.colorbar(_c, ax=ax,shrink=0.9,pad=.01,label=task["zlabel"])
    # plt.subplots_adjust(wspace=-0.2, hspace=-0.2)
    # plt.tight_layout()


    ax = axes[4]
    n_ele = np.trapz(y=spec, x=ys, axis=1)#[:,np.newaxis]
    n_ele_an = np.trapz(y=spec_an, x=ys_an, axis=1)#[:,np.newaxis]
    ax.plot(xs, n_ele/n_ele_an, color="black",ls='-',lw=1.5, label=f"integrated {name2}/{name1}")
    # ax.plot(xs, n_ele_an, color="black",ls='-',lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(task["xlabel"],fontsize=12)
    # ax.set_ylabel(label="Analytic/Numeric")
    ax.set_xlim(xs_an[0], xs_an[-1])
    ax.set_ylim(1e-2,1e2)
    ax.grid(linestyle=':')
    ax.legend()


    plt.savefig(os.getcwd()+'/'+task["figname"])
    if task["show"]: plt.show()

# def plot_electron_conservation(ej:PBA.Ejecta, ej_an:PBA.Ejecta):
#     fig,axes = plt.subplots(ncols=1,nrows=1)
#
#     n_ele = np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
def plot_two_observer_spectra(ej:PBA.Ejecta, ej_an:PBA.Ejecta, fs_or_rs :str = "fs", task={}):

    spec=ej.get_lc(key="fluxdens"+'_'+fs_or_rs,xkey="freqs",ykey="times",freq=None,time=None,ishell=0,ilayer=0,spec=False)
    spec_an=ej_an.get_lc(key="fluxdens"+'_'+fs_or_rs,xkey="freqs",ykey="times",freq=None,time=None,ishell=0,ilayer=0,spec=False)


    m2 = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3",ishell=0,ilayer=0)
    m2_an = ej_an.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3",ishell=0,ilayer=0)
    dm = np.diff(m2)
    dm_an = np.diff(m2_an)

    xs = ej.get_grid(key=xkey,spec=False)
    ys = ej.get_grid(key=ykey,spec=False)

    xs_an = ej_an.get_grid(key=xkey,spec=False)
    ys_an = ej_an.get_grid(key=ykey,spec=False)

    fig = plt.figure(figsize=(6, 6))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    # Draw the scatter plot and marginals.

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # colormesh
    ax.scatter(x, y)


def plot_comoving_spectrum():

    # run

    pba_fs_mix, pba_fs_num = run_2_pba(run=False, rs_dyn="no", rs_rad="no")
    pba_fs_mix_fsrs_dyn, pba_fs_num_fsrs_dyn = run_2_pba(run=False, rs_dyn="yes", rs_rad="no")
    pba_fsrs_mix_fsrs, pba_fsrs_num_fsrs = run_2_pba(run=False, rs_dyn="yes", rs_rad="yes")

    ''' RESOLUTION PLOT || COMOVING SPECTRA '''

    plot_compare_two_spectra_new(
        ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
        "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
        "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
        "figname":"abstract_ele_fs__fs_fs__num_vs_ana.png","show":True, "title":"FS; FS dynamics; numeric vs analytic"}
    )
    plot_compare_two_spectra_new(
        ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$j_{\rm nu}'$ [cgs]",
              "figname":"abstract_synch_fs__fs_fs__num_vs_ana.png","show":True, "title":"FS; FS dynamics; numeric vs analytic"}
    )

    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #           "figname":"abstract_ele_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$j_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_synch_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )

    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs.GRB, ej_an=pba_fsrs_num_fsrs.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #           "figname":"abstract_ele_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs.GRB, ej_an=pba_fsrs_num_fsrs.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$j_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_synch_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )





    plot_compare_two_spectra_new(
        ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]",
              "figname":"abstract_ele_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":"FS; with/without RS dynamics"}
    )
    plot_compare_two_spectra_new(
        ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]",
              "figname":"res_fig1.png","show":True, "title": r"FS; with/without RS \& RS dynamics"}
    )



    plot_compare_two_spectra_new(
        ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="mix", name2="num", is_spec=False,
        fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
        "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
        "zlabel":r"$\nu f_{\nu}$ ", "figname":"res_fig1.png",    "show":True,
        "title":"FS; Numeric Vs Analytic; FS dynamics"}
    )
    plot_compare_two_spectra_new(
        ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="mix", name2="num", is_spec=False,
        fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu f_{\nu}$ ", "figname":"res_fig1.png",    "show":True,
              "title":r"FS; Numeric Vs Analytic; FS \& RS dynamics"}
    )




    plot_compare_two_spectra_new(ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="mix", name2="num", is_spec=False,
                                 fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
                                 task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu f_{\nu}$ ", "figname":"spec_fs_ele.png","show":True,
                                       "title":"FS; Numeric Vs Analytic; FS dynamics"})
    plot_compare_two_spectra_new(ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="mix", name2="num", is_spec=False,
                                 fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
                                 task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu f_{\nu}$ ", "figname":"spec_fs_ele.png","show":True})
    plot_compare_two_spectra_new(ej=pba_fsrs_mix_fsrs.GRB, ej_an=pba_fsrs_num_fsrs.GRB, name1="mix", name2="num", is_spec=False,
                                 fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
                                 task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu f_{\nu}$ ", "figname":"spec_fs_ele.png",
                                       "title":"FS; Numeric Vs Analytic; FS \& RS dynamics", "show":True})

    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fs_num.GRB, name1="w RS", name2="w/o RS", is_spec=False,
                                 fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
                                 task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu f_{\nu}$ ", "figname":"spec_fs_ele.png","show":True})
    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fsrs_num_fsrs.GRB, name1="RSd", name2="RS", is_spec=False,
                                 fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs",
                                 task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu f_{\nu}$ ", "figname":"spec_fs_ele.png","show":True})

    # compare FS emission when RS dynamics is ON/OFF (ANALYTIC)
    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fs_num.GRB, name1="w RS", name2="w/o RS", is_spec=True,
                                 fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
                                 task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
                                 # "xs":(7.e3,1.e5,1.e7), "colors_xs":("cyan","green","red"),
                                 "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                 "zlabel":r"$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]",
                                 "figname":"fig1.png","show":True, "title":"FS; with/without RS dynamics"})
    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fs_num.GRB, name1="w RS", name2="w/o RS", is_spec=True,
                                 fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
                                 task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("cyan","green","red"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]",
                                       "figname":"fig2.png","show":True,
                                       "title":"RS; with/without RS dynamics"})
    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fs_num.GRB, name1="w RS", name2="w/o RS", is_spec=True,
                                 fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
                                 task={"ys":(1.e9,1e18,1e22), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
                                 # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                 "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                 "zlabel":r"$\nu j'_{\nu}$ ", "figname":"fig3.png","show":True,
                                       "title":"FS; with/without RS dynamics"})
    plot_compare_two_spectra_new(ej=pba_fs_num_fsrs_dyn.GRB, ej_an=pba_fs_num.GRB, name1="w RS", name2="w/o RS", is_spec=True,
                                 fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
                                 task={"ys":(1.e9,1e18,1e22), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
                                       # "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"),
                                       "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                                       "zlabel":r"$\nu j'_{\nu}$ ", "figname":"fig3.png","show":True,
                                       "title":"FS; with/without RS dynamics"})

    # compare FS emission when RS dynamics is ON/OFF (NUMERIC)
    plot_compare_two_spectra_new(ej=pba_fs_num.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, fs_or_rs ="fs", ele_syn_ssc ="n_ele", xkey="times_gams", ykey="gams",
                                 task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}$",
                               "xs":(7.e3,1.e5,1.e7), "colors_xs":("cyan","green","red"), "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                               "zlabel":r'$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]',
                               "figname":"spec_fs_ele.png"})


    plot_compare_two_spectra_new(ej=pba_fs_num.GRB, ej_an=pba_fs_mix.GRB, fs_or_rs ="fs", ele_syn_ssc ="n_ele", xkey="times_gams", ykey="gams",
                                 task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}$",
                           "xs":(7.e3,1.e5,1.e7), "colors_xs":("cyan","green","red"), "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
                           "zlabel":r'$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]',
                           "figname":"spec_fs_ele.png"})
    # plot_two_spectra_new(ej=pba_fsrs_num.GRB, ej_an=pba_fsrs_mix.GRB, fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams",
    #                      task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}$",
    #                            "xs":(7.e3,1.e5,1.e7), "colors_xs":("cyan","green","red"), "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #                            "zlabel":r'$N_{e;\, \rm tot}^{-1} [ \gamma_{e}^2 dN_{e}/d\gamma_{e}$ ]',
    #                            "figname":"spec_fs_ele.png"})
    # plot_two_spectra_new(ej=pba_fs_num.GRB, ej_an=pba_fs_mix.GRB, fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs",
    #                  task={"ys":(1.e9,1e18,1e22), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #                        "xs":(7.e3,1.e5,1.e7), "colors_xs":("lime","orange","magenta"), "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #                        "zlabel":r"$\nu j_{\nu}$ ",
    #                        "figname":"spec_fs_ele.png"})
    # plot_two_observer_spectra(ej=pba_num.GRB, ej_an=pba_mix.GRB)

if __name__ == '__main__':
    plot_comoving_spectrum()