# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import PyBlastAfterglow
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
    from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
        from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")


try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"


def tst_dynamics_fsrs(withSpread = False,
                      savefig = "compare_uniform_afgpy.png"):
    # prepare ID
    prepare_grb_ej_id_1d({"Eisso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
                          "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
                         type="pw",outfpath="tophat_grb_id.h5")
    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
                           newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")


    fig, axes = plt.subplots(ncols=1,nrows=1)
    ax = axes

    cmap = cm.viridis
    ishells=(1,)
    ilayers=(0,10,22),
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))

    ax.plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
            pba.GRB.get_dyn_arr(v_n="Gamma",ishell=0,ilayer=0),color="black",ls="-")
    plt.show()
def plot_ejecta_layers(ishells=(0,), ilayers=(0,25,49),
                       v_n_x = "R", v_n_ys = ("rho", "mom"), colors_by="layers",legend=False,
                       figname="dyn_layers_fsrs.png"):

    # prepare ID
    workdir = os.getcwd()+"/"
    # prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
    #                       "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
    #                      type="pw",outfpath="tophat_grb_id.h5")
    prepare_grb_ej_id_1d({"Eiso_c":2*1.e53, "Gamma0c": 1000., "M0c": -1.,"theta_c": np.pi/8., "theta_w": np.pi/2.,
                          "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},
                         type="pw",outfpath="tophat_grb_id.h5")

    # run fs-only model
    modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                    "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_tophad.h5"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba_fs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
    pba_fs.run(loglevel="info")

    # run fsrs model
    modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_light_curve":"lc_grb_tophad.h5"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba_fsrs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
    pba_fsrs.run(loglevel="info")

    # print(pba_fs.GRB.get_dyn_arr(v_n="M3",ishell=0,ilayer=0))

    # plot

    layers = []
    for i in ishells:
        for j in ilayers:
            layers.append("shell={} layer={}".format(i,j))

    # v_ns = ["Gamma"]

    # dfile = h5py.File(curdir+"magnetar_driven_ej.h5", "r")
    # print(dfile.keys())
    # print(dfile["layer=0"]["M2"])

    fid, axes = plt.subplots(ncols=1, nrows=len(v_n_ys), figsize=(6,6),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    for iv_n, v_n in enumerate(v_n_ys):
        ax = axes[iv_n] if len(v_n_ys) > 1 else axes
        i = 0
        for il, layer in enumerate(layers):
            x_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n_x,ishell=ishells[0],ilayer=il)#  np.array(dfile[layer][v_n_x])
            y_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n,ishell=ishells[0],ilayer=il)#np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
            i=i+1
        # --------------------------------
        i = 0
        for il, layer in enumerate(layers):
            x_arr = pba_fsrs.GRB.get_dyn_arr(v_n=v_n_x,ishell=ishells[0],ilayer=il)#  np.array(dfile[layer][v_n_x])
            y_arr = pba_fsrs.GRB.get_dyn_arr(v_n=v_n,ishell=ishells[0],ilayer=il)#np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            ax.plot(x_arr, y_arr, ls='--', color=color, label=layer)
            i=i+1

        ax.set_xlabel(v_n_x)
        if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]")
        ax.set_ylabel(v_n)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(5e-1,4e3)
    if legend: plt.legend()
    plt.savefig(workdir+figname, dpi=256)
    plt.show()



def tst_against_afgpy(withSpread = False,
                      savefig = "compare_uniform_afgpy.png",
                      load_data = True):

    # pba = PBA(workingdir=os.getcwd()+"/", readparfileforpaths=True)
# modify_parfile_par_opt
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax = axes

    # pba_0 = PBA(os.getcwd()+"/", readparfileforpaths=True)
    # pba_016 = PBA(os.getcwd()+"/", readparfileforpaths=True)

    prepare_grb_ej_id_1d({"Eisso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
                          "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
                         type="pw",outfpath="tophat_grb_id.h5")

    lls, lbls = [], []
    for (i_thetaobs, i_freq, i_color) in [
        # (thetaObs, freqobs, "blue"),
        (0.16, 1e9, "orange"),
        (0, 1e18, "green"),
        (0.16, 1.e18, "red"),
        (0, 1e9, "gray"),
    ]:

        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
                               parfile="default_parfile.par", newparfile="parfile2.par", keep_old=True)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
                               newopts={"method_synchrotron":"Joh06", "fname_light_curve":"tophat_{}_joh06.h5"
                               .format( str(i_thetaobs).replace(".",""))},
                               parfile="parfile2.par", newparfile="parfile.par",keep_old=False)
        pba = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
        # pba.reload_parfile()


        pba.run()

        ax.plot(pba.GRB.get_lc_times() / cgs.day,
                pba.GRB.get_lc_totalflux(freq=i_freq), color=i_color, ls='-',
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        pba.clear()
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
                               parfile="default_parfile.par", newparfile="parfile2.par", keep_old=True)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",newpars={},
                               newopts={"method_synchrotron":"WSPN99", "fname_light_curve":"tophat_{}_WSPN99.h5"
                               .format( str(i_thetaobs).replace(".",""))},
                               parfile="parfile2.par", newparfile="parfile.par",keep_old=False)
        pba.reload_parfile()

        pba.run()

        ax.plot(pba.GRB.get_lc_times() / cgs.day,
                pba.GRB.get_lc_totalflux(freq=i_freq), color=i_color, ls=':',
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        pba.clear()

        if load_data:
            if withSpread:
                fname = "afterglowpy_theta{:d}_lognu{:d}_spread.txt".format(int(i_thetaobs * 180 / np.pi),
                                                                            int(np.log10(i_freq)))
            else:
                fname = "afterglowpy_theta{:d}_lognu{:d}.txt".format(int(i_thetaobs * 180 / np.pi),
                                                                     int(np.log10(i_freq)))
            _t, _ref_F_afgpy, _ref_F = np.loadtxt(os.getcwd() + '/' + fname, unpack=True)
            _ll, = ax.plot(_t / cgs.day, _ref_F_afgpy, color=i_color, ls='--')
            lls.append(_ll)
            lbls.append(r"$\nu=$" + r"${}$ Hz ".format(latex_float(i_freq))
                        + r"$\theta_{\rm obs}=$" + r"{:.1f} deg".format(i_thetaobs * 180 / np.pi))

        # break

    l11, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='-', label=r"\texttt{PBA} with J06")
    l12, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls=':', label=r"\texttt{PBA} with WSPN99")
    l13, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"\texttt{afterglowpy}")

    legend1 = plt.legend([l11, l12, l13],
                         # [r"\& J\'{o}hannesson+06", r"\& WSPN+99", r"\texttt{afterglowpy}"],
                         [r"\texttt{PBA} with J06", r"\texttt{PBA} with WSPN99", r"\texttt{afterglowpy}"],
                         loc="center", bbox_to_anchor=(0.78, 0.56), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)
    legend2 = plt.legend(lls, lbls,
                         loc="center", bbox_to_anchor=(0.4, 0.16), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")

    # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   # labelsize=plotdic["fontsize"],
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [mJy]")
    ax.set_xlim(1e-1, 1e3)
    ax.set_ylim(1e-9, 1e2)
    # ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(PAPERPATH + save)
    if savefig: plt.savefig(os.getcwd() + '/' + savefig, dpi=256)

    plt.show()

if __name__ == '__main__':
    plot_ejecta_layers(ishells=(0,), ilayers=(0,),
                       v_n_x = "R", v_n_ys = ("B", "B_rs", "gamma_min", "gamma_min_rs","gamma_c", "gamma_c_rs"), colors_by="layers",legend=False,
                       figname="dyn_layers_fs.png")
    # tst_against_afgpy()
    exit(0)