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
import shutil
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

def run_2_pba(n_ism=1e-2,i_thetaobs = 0., run=True, rs_dyn="yes", rs_rad="yes", ssa="yes", adi="yes"):
    working_dir = os.getcwd()+'/'+'working_dir/'
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    shutil.copyfile(os.getcwd()+'/'+"parfile_def.par",
                    working_dir+"parfile_def.par")

    struct = {"struct":"tophat",
              "Eiso_c":1.e52, "Gamma0c": 350., "M0c": -1.,
              "theta_c": np.pi / 10, "theta_w": np.pi / 10}
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80, n_layers_a=1)

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir+"tophat_grb_id_a.h5")


    th = str(i_thetaobs).replace(".","")
    nsim = str(np.log10(n_ism)).replace(".","")

    key = f"logn{nsim}_th{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_ssa-{ssa}_an"
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="main",newpars={"theta_obs":i_thetaobs,"n_ism":n_ism},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_synchrotron_rs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "method_ne_rs":"useNe",
                                                      "use_ssa_fs":ssa,
                                                      "use_ssa_rs":ssa,
                                                      "method_gamma_c_fs" : "useTcomov",
                                                      "method_gamma_c_rs" : "useTcomov",
                                                      "method_gamma_min_fs":"useNumericGamma",
                                                      "method_gamma_min_rs":"useNumericGamma",
                                                      "method_nonrel_dist_fs":"use_gamma_min",
                                                      "method_nonrel_dist_rs":"use_gamma_min",
                                                      "limit_lf_min_to1_fs":"yes",
                                                      "limit_lf_min_to1_rs":"yes",
                                                      "method_ele_fs":"analytic",
                                                      "method_ele_rs":"analytic",
                                                      "do_rs" : rs_dyn,
                                                      "do_rs_radiation": rs_rad,
                                                      "bw_type":"fsrs" if rs_dyn=="yes" else "fs",
                                                      "fname_ejecta_id":"tophat_grb_id_a.h5",
                                                      "fname_dyn":f"dyn_{key}.h5",
                                                      "fname_light_curve":f"lc_{key}.h5",
                                                      "fname_spectrum":f"spectrum_{key}.h5"},
                                             parfile="parfile.par", newparfile=f"parfile_{key}.par", keep_old=False)
    pba_an = PBA.interface.PyBlastAfterglow(workingdir=working_dir, parfile=f"parfile_{key}.par")
    if run: pba_an.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")

    # ----------------------------------------
    key = f"logn{nsim}_th{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_ssa-{ssa}_mix"
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="main",newpars={"theta_obs":i_thetaobs,"n_ism":n_ism},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_synchrotron_rs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "method_ne_rs":"useNe",
                                                      "use_ssa_fs":ssa,
                                                      "use_ssa_rs":ssa,
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
                                                      "fname_dyn":f"dyn_{key}.h5",
                                                      "fname_light_curve":f"lc_{key}.h5",
                                                      "fname_spectrum":f"spectrum_{key}.h5"},
                                             parfile="parfile.par", newparfile=f"parfile_{key}.par", keep_old=False)
    pba_mix = PBA.interface.PyBlastAfterglow(workingdir=working_dir, parfile=f"parfile_{key}.par")
    if run: pba_mix.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")

    # ----------------------------------------
    key = f"logn{nsim}_th{th}_rsdyn-{rs_dyn}_rsrad-{rs_rad}_ssa-{ssa}_adi-{adi}_num"
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="main",newpars={"theta_obs":i_thetaobs,"n_ism":n_ism},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=working_dir, part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_synchrotron_rs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "use_ssa_fs":ssa,
                                                      "use_ssa_rs":ssa,
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
                                                      "num_ele_use_adi_loss_fs":adi,
                                                      "num_ele_use_adi_loss_rs":adi,
                                                      "do_rs" : rs_dyn,
                                                      "do_rs_radiation": rs_rad,
                                                      "bw_type":"fsrs" if rs_dyn=="yes" else "fs",
                                                      "fname_ejecta_id":"tophat_grb_id_a.h5",
                                                      "fname_dyn":f"dyn_{key}.h5",
                                                      "fname_light_curve":f"lc_{key}.h5",
                                                      "fname_spectrum":f"spectrum_{key}.h5"},
                                             parfile="parfile.par", newparfile=f"parfile_{key}.par", keep_old=False)
    pba_num = PBA.interface.PyBlastAfterglow(workingdir=working_dir, parfile=f"parfile_{key}.par")
    if run: pba_num.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
    # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
    #         pba_a2.GRB.get_lc_totalflux(freq=i_freq), color=i_color, ls='-', lw=1,
    #         label=r" $\nu$={:.1e} Ne numeric".format(i_freq))
    # pba_a2.clear()
    return (pba_an, pba_mix, pba_num)

def plot_two_spectra(ej:PBA.Ejecta, ej_an:PBA.Ejecta, fs_or_rs :str = "fs", ele_syn_ssc:str = "ele",
                     xkey="times_gams", ykey="gams",
                     task={"ys":(1e9,1e18,1e22), "colors_freqs":("blue","green","red"), "ylim":(1,1e9),
                           "xs":(1.e3,1.e5,1.e7), "colors_freqs":("blue","green","red"), "xlim":(2e3, 1e8),
                           "figname":"spec_fs_ele.png"}):
    mp = 1.6726e-24
    qe = 4.803204e-10
    me = 9.1094e-28
    c  = 2.99792458e10

    spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs, xkey=xkey, key_time=ykey, freq=None, time=None, ishell=0, ilayer=0, spec=True)
    spec_an=ej_an.get_lc(key=ele_syn_ssc+'_'+fs_or_rs, xkey=xkey, key_time=ykey, freq=None, time=None, ishell=0, ilayer=0, spec=True)


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
def _normalize_spec(spec:np.ndarray, xs:np.ndarray, ys:np.ndarray, norm_method:None or str, mask_negatives=True):
    if not norm_method:
        return spec
    if norm_method == "/integ *y^2":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        spec *= np.power(ys, 2)[np.newaxis,:]
    elif norm_method == "*y":
        spec *= ys[np.newaxis,:]
    elif norm_method == "*y^2":
        spec *= ys[np.newaxis,:] ** 2
    elif norm_method == "/integ":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
    else:
        raise KeyError("Norm method is not recognized")
    if mask_negatives:
        spec[~np.isfinite(spec)] = 1e-100
        spec[spec <= 0] = 1e-100
    return spec
def plot_compare_two_spectra_new(ej:PBA.Ejecta, ej_an:PBA.Ejecta,
                                 ej_sustract:PBA.Ejecta or None = None,
                                 ej_an_sustract:PBA.Ejecta or None = None,
                                 name1:str= "Numeric", name2="Analytic", is_spec=True,
                                 fs_or_rs :str = "fs", ele_syn_ssc:str = "ele", norm_method:str or None="integ_gam2",
                                 xkey="times_gams", ykey="gams",
                                 task=dict()):
    mp = 1.6726e-24
    qe = 4.803204e-10
    me = 9.1094e-28
    c  = 2.99792458e10


    # plt.loglog(ej.get_dyn_arr(v_n='B_rs',ishell=0,ilayer=0))
    # plt.loglog(ej.get_dyn_arr(v_n='B',ishell=0,ilayer=0))
    # plt.show()


    fig, axes = plt.subplots(ncols=1, nrows=5, figsize=(7, 8), sharex="col", #sharey="row",sharex="col",
                             gridspec_kw=dict(height_ratios=[1.0,1.,1.,1.,0.4],hspace=-0.1,wspace=0.01),#dict(height_ratios=[0.5,1.2]),
                             layout='constrained')


    spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                   xkey=xkey, key_time=ykey, freq=None, time=None, ishell=0, ilayer=0,
                   spec=is_spec, sum_shells_layers=False)
    xs = ej.get_grid(key=xkey,spec=is_spec)
    ys = ej.get_grid(key=ykey,spec=is_spec)
    if fs_or_rs == "rs" and not ej_sustract is None:
        spec_sub=ej_sustract.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                       xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                       spec=is_spec,sum_shells_layers=False)
        xs_ = ej_sustract.get_grid(key=xkey,spec=is_spec)
        ys_ = ej_sustract.get_grid(key=ykey,spec=is_spec)
        if not np.array_equal(xs,xs_):
            raise ValueError("grid mismatch")
        if not np.array_equal(ys,ys_):
            raise ValueError("grid mismatch")
        spec -= spec_sub
    spec = _normalize_spec(spec=spec,xs=xs,ys=ys,norm_method=norm_method,mask_negatives=True)

    spec_an=ej_an.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                         xkey=xkey, key_time=ykey, freq=None, time=None, ishell=0, ilayer=0,
                         spec=is_spec, sum_shells_layers=False)
    xs_an = ej_an.get_grid(key=xkey,spec=is_spec)
    ys_an = ej_an.get_grid(key=ykey,spec=is_spec)
    if fs_or_rs == "rs" and not ej_sustract is None:
        spec_sub_an=ej_an_sustract.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                                    xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                                    spec=is_spec,sum_shells_layers=False)
        xs_an_ = ej_an_sustract.get_grid(key=xkey,spec=is_spec)
        ys_an_ = ej_an_sustract.get_grid(key=ykey,spec=is_spec)
        if not np.array_equal(xs_an,xs_an_):
            raise ValueError("grid mismatch")
        if not np.array_equal(ys_an,ys_an_):
            raise ValueError("grid mismatch")
        spec_an -= spec_sub_an
    spec_an = _normalize_spec(spec=spec_an,xs=xs_an,ys=ys_an,norm_method=norm_method,mask_negatives=True)

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
    # ax_ = ax.twinx()
    # ax_.plot(xs_an,ej_an.get_dyn_arr(v_n="B_rs",ilayer=0,ishell=0),color='gray',ls='-',label="B")
    # # ax_.plot(xs_an,ej_an.get_dyn_arr(v_n="GammaRsh",ilayer=0,ishell=0),color='gray',ls='--',label="GammaRsh")
    # ax_.set_xscale("log"); ax_.set_yscale('log'); ax_.legend(), ax_.set_ylim(1e-2,1e3)

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

def plot_(ej_1:PBA.Ejecta, ej_2:PBA.Ejecta,
          fs_or_rs_1 :str = "fs", fs_or_rs_2 :str = "fs"):

    fig, ax = plt.subplots(ncols=1,nrows=1)

    # r = ej.get_dyn_arr(v_n="R",ishell=0,ilayer=0)
    # dr = ej.get_dyn_arr(v_n="thickness",ishell=0,ilayer=0)
    sigmaT = 6.6524e-25
    me = 9.1094e-28
    c  = 2.99792458e10
    gamma_max = ej_1.get_dyn_arr(v_n="gamma_max" if fs_or_rs_1 == "fs" else "gamma_max_rs", ishell=0, ilayer=0)
    B = ej_1.get_dyn_arr(v_n="B" if fs_or_rs_1 == "fs" else "B_rs", ishell=0, ilayer=0)
    tcool_syn = sigmaT * gamma_max * B * B / (6. * np.pi * me * c)

    ax.plot(ej_1.get_dyn_arr(v_n="tburst", ishell=0, ilayer=0), tcool_syn, ls='-',color='blue', label="tcool_syn-"+fs_or_rs_1)

    m = ej_1.get_dyn_arr(v_n="M2" if fs_or_rs_1 == "fs" else "M3", ishell=0, ilayer=0)
    rho = ej_1.get_dyn_arr(v_n="rho2" if fs_or_rs_1 == "fs" else "rho3", ishell=0, ilayer=0)
    t = ej_1.get_dyn_arr(v_n="tcomov", ishell=0, ilayer=0)
    vol = m[:-1]/rho[:-1]
    vol_p1 = m[1:]/rho[1:]
    dt = np.diff(t)
    dlnVdt = 1. / dt * (1. - vol / vol_p1)
    tcool_adi = (gamma_max[1:]*gamma_max[1:]-1.)/(3.*gamma_max[1:]*gamma_max[1:])*dlnVdt

    ax.plot(ej_1.get_dyn_arr(v_n="tburst", ishell=0, ilayer=0)[1:], tcool_adi, ls='-',color='red', label="tcool_adi-"+fs_or_rs_1)



    gamma_max = ej_2.get_dyn_arr(v_n="gamma_max" if fs_or_rs_2 == "fs" else "gamma_max_rs", ishell=0, ilayer=0)
    B = ej_2.get_dyn_arr(v_n="B" if fs_or_rs_2 == "fs" else "B_rs", ishell=0, ilayer=0)
    tcool_syn = sigmaT * gamma_max * B * B / (6. * np.pi * me * c)

    ax.plot(ej_2.get_dyn_arr(v_n="tburst", ishell=0, ilayer=0), tcool_syn, ls='--',color='blue', label="tcool_syn-"+fs_or_rs_2)

    m = ej_2.get_dyn_arr(v_n="M2" if fs_or_rs_2 == "fs" else "M3", ishell=0, ilayer=0)
    rho = ej_2.get_dyn_arr(v_n="rho2" if fs_or_rs_2 == "fs" else "rho3", ishell=0, ilayer=0)
    t = ej_2.get_dyn_arr(v_n="tcomov", ishell=0, ilayer=0)
    vol = m[:-1]/rho[:-1]
    vol_p1 = m[1:]/rho[1:]
    dt = np.diff(t)
    dlnVdt = 1. / dt * (1. - vol / vol_p1)
    tcool_adi = (gamma_max[1:]*gamma_max[1:]-1.)/(3.*gamma_max[1:]*gamma_max[1:])*dlnVdt

    ax.plot(ej_2.get_dyn_arr(v_n="tburst", ishell=0, ilayer=0)[1:], tcool_adi, ls='--',color='red', label="tcool_adi-"+fs_or_rs_2)

    plt.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("tburst",fontsize=12)
    plt.show()

def plot_comoving_spectrum():

    # run
    run = False
    # n_ism = 1.e-3
    n_ism = 10.
    pba_fs_an, pba_fs_mix, pba_fs_num = run_2_pba(n_ism=n_ism,run=run, rs_dyn="no", rs_rad="no", ssa="no", adi="yes")
    _, _, pba_fs_noadi_num = run_2_pba(n_ism=n_ism,run=run, rs_dyn="no", rs_rad="no", ssa="no", adi="no")
    pba_fs_ssa_an, pba_fs_ssa_mix, pba_fs_ssa_num = run_2_pba(n_ism=n_ism,run=run, rs_dyn="no", rs_rad="no", ssa="yes", adi="yes")
    pba_fs_an_fsrs_dyn, pba_fs_mix_fsrs_dyn, pba_fs_num_fsrs_dyn = run_2_pba(n_ism=n_ism,run=run, rs_dyn="yes", rs_rad="no",ssa="no", adi="yes")
    pba_fs_an_fsrs_ssa_dyn, pba_fs_mix_fsrs_ssa_dyn, pba_fs_num_fsrs_ssa_dyn = run_2_pba(n_ism=n_ism,run=run, rs_dyn="yes", rs_rad="no",ssa="yes", adi="yes")
    pba_fsrs_an_fsrs_dyn, pba_fsrs_mix_fsrs_dyn, pba_fsrs_num_fsrs_dyn = run_2_pba(n_ism=n_ism,run=run, rs_dyn="yes", rs_rad="yes", ssa="no", adi="yes")
    _, _, pba_fsrs_num_fsrs_noadi_dyn = run_2_pba(n_ism=n_ism,run=run, rs_dyn="yes", rs_rad="yes", ssa="no", adi="no")
    pba_fsrs_an_fsrs_ssa_dyn, pba_fsrs_mix_fsrs_ssa_dyn, pba_fsrs_num_fsrs_ssa_dyn = run_2_pba(n_ism=n_ism,run=run, rs_dyn="yes", rs_rad="yes", ssa="yes", adi="yes")

    if run: exit(0)




    ''' EFFECT OF ADDIABATIC COOLING ON ELECTRON/SYNCHROTRON SPECTRA '''
    # plot_(ej_1=pba_fs_num.GRB, ej_2=pba_fsrs_num_fsrs_dyn.GRB, fs_or_rs_1="fs", fs_or_rs_2="rs")
    ''' FS '''
    plot_compare_two_spectra_new(
        ej=pba_fs_noadi_num.GRB, ej_an=pba_fs_num.GRB, name1="no Adi", name2="w Adi", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y^2",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
        "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
        "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
        "figname":"abstract_ele_fs__fs_fs__num_adi_vs_noadi.png","show":True, "title":"FS; FS dynamics; numeric; adiabatic cooling"}
    )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_noadi_num.GRB, ej_an=pba_fs_num.GRB, name1="no Adi", name2="w Adi", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$j_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_synch_fs__fs_fs__num__adi_vs_noadi.png","show":True, "title":"FS; FS dynamics; numeric; adiabatic cooling"}
    # )
    ''' RS '''
    # plot_(ej_1=pba_fsrs_num_fsrs_dyn.GRB, ej_2=pba_fsrs_num_fsrs_dyn.GRB, fs_or_rs_1="fs", fs_or_rs_2="rs")
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_num_fsrs_noadi_dyn.GRB, ej_an=pba_fsrs_num_fsrs_dyn.GRB, name1="no Adi", name2="w Adi", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y^2",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #           "figname":"abstract_ele_rs__fsrs_fsrs__num__adi_vs_noadi.png","show":True, "title":r"RS; FS \& RS dynamics; numeric; adiabatic cooling"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_num_fsrs_noadi_dyn.GRB, ej_an=pba_fsrs_num_fsrs_dyn.GRB,  name1="no Adi", name2="w Adi", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$\nu' j_{nu'}'$ [cgs]",
    #           "figname":"abstract_synch_rs__fsrs_fsrs__num__adi_vs_noadi.png","show":True, "title":r"RS; FS \& RS dynamics; numeric; adiabatic cooling"}
    # )


    ''' ANALYTIC VS NUMERIC PLOT || COMOVING SPECTRA '''
    ''' FS '''
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y^2",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #     "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #     "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #     "figname":"abstract_ele_fs__fs_fs__num_vs_ana.png","show":True, "title":"FS; FS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$j_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_synch_fs__fs_fs__num_vs_ana.png","show":True, "title":"FS; FS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "ssa", xkey="times_freqs", ykey="freqs", norm_method="*y^2",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$\nu^{'2} a_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_ssa_fs__fs_fs__num_vs_ana.png","show":True, "title":"FS; FS dynamics; numeric vs analytic"}
    # )

    ''' FS with RS dynamics '''
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y^2",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #           "figname":"abstract_ele_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$j_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_synch_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_mix_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "fs", ele_syn_ssc = "ssa", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$\nu^{'2} a_{\rm nu}'$ [cgs]",
    #           "figname":"abstract_ssa_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )

    ''' RS '''
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs_dyn.GRB, ej_an=pba_fsrs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y^2",
    #     task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\gamma_{e}'$",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$N_{e;\, \rm tot}^{'-1} [ \gamma_{e}^{'2} dN'_{e}/d\gamma'_{e}$ ]",
    #           "figname":"abstract_ele_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs_dyn.GRB, ej_an=pba_fsrs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$\nu' j_{nu'}'$ [cgs]",
    #           "figname":"abstract_synch_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs_ssa_dyn.GRB, ej_an=pba_fsrs_num_fsrs_ssa_dyn.GRB, name1="analytic", name2="numeric", is_spec=True,
    #     fs_or_rs = "rs", ele_syn_ssc = "ssa", xkey="times_freqs", ykey="freqs", norm_method="*y^2",
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
    #           "zlabel":r"$\nu^{'2} a_{nu'}'$ [cgs]",
    #           "figname":"abstract_ssa_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )



    ''' ANALYTIC VS NUMERIC PLOT || OBSERVER SPECTRA '''

    # plot_compare_two_spectra_new(
    #     ej=pba_fs_an.GRB, ej_an=pba_fs_num.GRB, name1="analytic", name2="numeric", is_spec=False,
    #     fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs", norm_method=None,
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm obs}$ [s]",
    #           "zlabel":r"$F_{\rm nu}$ [mJy]",
    #           "figname":"abstract_fluxdens_fs__fs_fs__num_vs_ana.png","show":True, "title":r"FS; FS dynamics; numeric vs analytic"}
    # )
    # #
    # plot_compare_two_spectra_new(
    #     ej=pba_fs_an_fsrs_dyn.GRB, ej_an=pba_fs_num_fsrs_dyn.GRB, name1="analytic", name2="numeric", is_spec=False,
    #     fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs", norm_method=None,
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm obs}$ [s]",
    #           "zlabel":r"$F_{\rm nu}$ [mJy]",
    #           "figname":"abstract_fluxdens_fs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"FS; FS \& RS dynamics; numeric vs analytic"}
    # )
    #
    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_mix_fsrs_dyn.GRB, ej_an=pba_fsrs_num_fsrs_dyn.GRB,
    #     ej_sustract=pba_fs_mix_fsrs_dyn.GRB, ej_an_sustract=pba_fs_num_fsrs_dyn.GRB,
    #     name1="analytic", name2="numeric", is_spec=False,
    #     fs_or_rs = "rs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs", norm_method=None,
    #     task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm obs}$ [s]",
    #           "zlabel":r"$F_{\rm nu}$ [mJy]",
    #           "figname":"abstract_fluxdens_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )


    ''' EFFECT OF SELF-ABSORPTION ON OBSERVED SPECTRUM '''

    # plot_compare_two_spectra_new(
    #     ej=pba_fs_num.GRB, ej_an=pba_fs_ssa_num.GRB, name1="w/o SSA", name2="SSA", is_spec=False,
    #     fs_or_rs = "fs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs", norm_method=None,
    #     task={"ys":(1.e9,1e10,1e11), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm obs}$ [s]",
    #           "zlabel":r"$F_{\rm nu}$ [mJy]",
    #           "figname":"abstract_ssa_fluxdens_fs__fs_fs__num_vs_ana.png","show":True, "title":r"FS; FS dynamics; numeric vs analytic"}
    # )

    # plot_compare_two_spectra_new(
    #     ej=pba_fsrs_num_fsrs_dyn.GRB, ej_an=pba_fsrs_num_fsrs_ssa_dyn.GRB,
    #     ej_sustract=pba_fs_num_fsrs_dyn.GRB, ej_an_sustract=pba_fs_num_fsrs_ssa_dyn.GRB,
    #     name1="w/o SSA", name2="SSA", is_spec=False,
    #     fs_or_rs = "rs", ele_syn_ssc = "fluxdens", xkey="times", ykey="freqs", norm_method=None,
    #     task={"ys":(1.e9,1e10,1e11), "colors_ys":("cyan","green","red"), "ylim":(), "ylabel":r"$\nu$ [Hz]",
    #           "xlim":(), "xlabel":r"$t_{\rm obs}$ [s]",
    #           "zlabel":r"$F_{\rm nu}$ [mJy]",
    #           "figname":"abstract_ssa_fluxdens_rs__fsrs_fsrs__num_vs_ana.png","show":True, "title":r"RS; FS \& RS dynamics; numeric vs analytic"}
    # )


if __name__ == '__main__':
    plot_comoving_spectrum()