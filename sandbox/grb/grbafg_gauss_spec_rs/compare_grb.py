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
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
        from package.src.PyBlastAfterglowMag.id_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")


try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"


pars = {"Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.4,
        "nlayers_pw": 50, "nlayers_a": 0, "struct":"gaussian"}


def plot_ejecta_layers(ishells=(0,), ilayers=(0,25,49), nlayers_a=20,
                       v_n_x = "R", v_n_ys = ("rho", "mom"), colors_by="layers",legend=False,
                       figname="dyn_layers_fsrs.png", run_fs_only=False):

    # prepare ID
    workdir = os.getcwd()+"/"
    # prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
    #                       "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
    #                      type="pw",outfpath="tophat_grb_id.h5")
    theta_h = np.pi/2.
    one_min_cos = 2. * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h)
    ang_size_layer = 2.0 * np.pi * one_min_cos / (4.0 * np.pi)
    pars["nlayers_a"] = nlayers_a
    prepare_grb_ej_id_1d(pars, type="a",outfpath="gauss_grb_id.h5")

    ### run fs-only model
    if(run_fs_only):
        modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                               newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no", "spec_var_out":"em",
                                        "fname_ejecta_id":"gauss_grb_id.h5",
                                        "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_tophad.h5", "method_eats": "adaptive"},
                               parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
        pba_fs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
        pba_fs.run(loglevel="info")

    # run fsrs model
    modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes", "spec_var_out":"em_rs",
                                    "fname_ejecta_id":"gauss_grb_id.h5",
                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_spectrum":"lc_grb_tophad.h5", "method_eats": "adaptive"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba_fsrs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
    pba_fsrs.run(loglevel="info")

    # ref = RefData(workdir)


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

    fid, axes = plt.subplots(ncols=1, nrows=len(v_n_ys), figsize=(4.6+2,4.2+2),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap_fs = cm.Blues
    cmap_fsrs = cm.Reds
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    for iv_n, v_n in enumerate(v_n_ys):
        ax = axes[iv_n] if len(v_n_ys) > 1 else axes
        i = 0
        if(run_fs_only):
            for il, layer in enumerate(layers):
                x_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n_x,ishell=ishells[0],ilayer=ilayers[il])#  np.array(dfile[layer][v_n_x])
                y_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n,ishell=ishells[0],ilayer=ilayers[il])#np.array(dfile[layer][v_n])
                y_arr = y_arr[x_arr > 0]
                x_arr = x_arr[x_arr > 0]
                # if (v_n == "R"):
                # y_arr = y_arr/y_arr.max()
                if (colors_by=="layers"): color=cmap_fs(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
                else: color=cmap_fs(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
                if (v_n_x == "tburst"): x_arr /=cgs.day;
                ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
                i=i+1
        # --------------------------------
        i = 0
        for il, layer in enumerate(layers):
            x_arr = pba_fsrs.GRB.get_dyn_arr(v_n=v_n_x,ishell=ishells[0],ilayer=ilayers[il])#  np.array(dfile[layer][v_n_x])
            y_arr = pba_fsrs.GRB.get_dyn_arr(v_n=v_n,ishell=ishells[0],ilayer=ilayers[il])#np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap_fsrs(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap_fsrs(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            ax.plot(x_arr, y_arr, ls='--', color=color, label=layer)
            i=i+1

            # --- plot ref data
            # x_arr = ref.get(il, v_n_x)
            # if (v_n_x == "tburst"): x_arr /=cgs.day;
            # ax.plot(x_arr, ref.get(il, v_n), ls=':', color=color, lw=2.)

        if (v_n=="B"):
            ax.set_ylim(1e-3,10)
        elif(v_n=="B_rs"):
            ax.set_ylim(1e-2,1e4)
        elif(v_n=="gamma_min"):
            ax.set_ylim(30,1e6)
        elif(v_n=="gamma_min_rs"):
            ax.set_ylim(1e-5,5e3)
        elif(v_n=="gamma_c"):
            ax.set_ylim(73e3,1e9)
        elif(v_n=="gamma_c_rs"):
            ax.set_ylim(1,1e9)

        ax.set_ylabel(v_n,fontsize=12)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.set_xlim(1e-2,1e7)
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(5e-1,4e3)
    if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]",fontsize=12)
    if legend: plt.legend()
    plt.tight_layout()
    plt.savefig(workdir+figname, dpi=256)
    plt.show()

def plot_ejecta_layers_spec(freq=1e18,ishells=(0,), ilayers=(0,25,49),nlayers_a=20,colors_by="layers",legend=False,
                       figname="dyn_layers_fsrs.png", run_fs_only=True,type="a",method_eats="adaptive"):

    # prepare ID
    workdir = os.getcwd()+"/"
    # prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
    #                       "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
    #                      type="pw",outfpath="tophat_grb_id.h5")
    theta_h = np.pi/2.
    one_min_cos = 2. * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h)
    ang_size_layer = 2.0 * np.pi * one_min_cos / (4.0 * np.pi)

    pars["nlayers_a"] = nlayers_a
    prepare_grb_ej_id_1d(pars, type=type,outfpath="gauss_grb_id.h5")

    ### run fs-only model
    if(run_fs_only):
        modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                               newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                        "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_tophad.h5",
                                        "fname_ejecta_id":"gauss_grb_id.h5",
                                        "fname_spectrum_layers":"spec_fs_layers.h5", "method_eats": method_eats,
                                        "spec_var_out":"em"},
                               parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
        pba_fs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
        pba_fs.run(loglevel="info")

    # run fsrs model
    modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_spectrum":"lc_grb_tophad.h5",
                                    "fname_ejecta_id":"gauss_grb_id.h5",
                                    "fname_spectrum_layers":"spec_fsrs_layers.h5", "method_eats": method_eats,
                                    "spec_var_out":"em_rs"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba_fsrs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
    pba_fsrs.run(loglevel="info")

    # ref = RefData(workdir)


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

    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(4.6,4.2),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap_fs = cm.Blues
    cmap_fsrs = cm.Reds
    # mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))
    mynorm = Normalize(vmin=0,vmax=ilayers[-1])#norm(len(ishells)*len(ilayers))

    ax = axes[0]
    ax.set_ylim(5e-31,2e-8)
    ax.set_ylabel("Emissivity", fontsize=12)
    i = 0
    if(run_fs_only):
        for il, layer in enumerate(layers):
            x_arr = pba_fs.GRB.get_lc_times(spec=True)#  np.array(dfile[layer][v_n_x])
            y_arr = pba_fs.GRB.get_lc(freq=freq,ishell=ishells[0],ilayer=ilayers[il],spec=True)#np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            # color=cmap(mynorm(int(i)))
            color=cmap_fs(mynorm(int(ilayers[il])))
            x_arr /=cgs.day
            ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
            i=i+1
    # fig.colorbar(cm.ScalarMappable(norm=mynorm, cmap=cmap),ax=axes,label="Layers")
    # --------------------------------
    i = 0
    for il, layer in enumerate(layers):
        x_arr = pba_fsrs.GRB.get_lc_times(spec=True)#  np.array(dfile[layer][v_n_x])
        y_arr = pba_fsrs.GRB.get_lc(freq=freq,ishell=ishells[0],ilayer=ilayers[il],spec=True)#np.array(dfile[layer][v_n])
        y_arr = y_arr[x_arr > 0]
        x_arr = x_arr[x_arr > 0]
        # if (v_n == "R"):
        # y_arr = y_arr/y_arr.max()
        # color=cmap(mynorm(int(i)))
        color=cmap_fsrs(mynorm(int(ilayers[il])))
        x_arr /=cgs.day
        ax.plot(x_arr, y_arr, ls='--', color=color, label=layer)
        i=i+1

        # --- plot ref data
        # x_arr = ref.get(il, v_n_x)
        # if (v_n_x == "tburst"): x_arr /=cgs.day;
        # ax.plot(x_arr, ref.get(il, v_n), ls=':', color=color, lw=2.)

    ax = axes[1]
    ax.set_ylabel("Total Emissivity", fontsize=12)
    ax.set_ylim(1e-25,1e-6)
    if(run_fs_only):
        ax.plot(pba_fs.GRB.get_lc_times(spec=True)/cgs.day,
                pba_fs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='-', color='blue', label='FS')
    ax.plot(pba_fsrs.GRB.get_lc_times(spec=True)/cgs.day,
            pba_fsrs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='--', color='red', label='RS')
    ax.legend()

    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()

        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_xlim(1e-2,1e5)
    ax.set_xlabel("tburst" + " [day]",fontsize=12)
    plt.tight_layout()
    if legend: plt.legend()
    # plt.title("Comoving Spectrum Gaussian Jet")
    plt.savefig(workdir+figname, dpi=256)
    plt.show()

def plot_tst_total_spec_resolution(freq=1e9, nlayers=(10,20,40,80,120),legend=False,
                                   figname="dyn_layers_fsrs.png", include_fs_only=True,type="pw",method_eats="piece-wise"):
    workdir = os.getcwd()+"/"

    fid, ax = plt.subplots(ncols=1, nrows=1, figsize=(4.6+2,3.6+2),sharex="all")
    cmap_rs = cm.Reds
    cmap_fs = cm.Blues
    mynorm = Normalize(vmin=0,vmax=nlayers[-1])#norm(len(ishells)*len(ilayers))
    for i, i_nlayers in enumerate(nlayers):
        pars["nlayers_a"] = i_nlayers
        prepare_grb_ej_id_1d(pars, type=type,outfpath="gauss_grb_id.h5")

        ### run fs-only model
        if(include_fs_only):
            modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                                   newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                            "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"spec_grb_gauss.h5",
                                            "fname_spectrum_layers":"spec_fs_layers.h5", "method_eats": method_eats,
                                            "fname_ejecta_id":"gauss_grb_id.h5",
                                            "spec_var_out":"em"},
                                   parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
            pba_fs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
            pba_fs.run(loglevel="info")

        # run fsrs model
        modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                               newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                        "fname_dyn":"dyn_bw_fsrs.h5","fname_spectrum":"spec_grb_gauss.h5",
                                        "fname_spectrum_layers":"spec_fsrs_layers.h5", "method_eats": method_eats,
                                        "fname_ejecta_id":"gauss_grb_id.h5",
                                        "spec_var_out": "em_rs"},
                               parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
        pba_fsrs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
        pba_fsrs.run(loglevel="info")

        color_fs=cmap_fs(mynorm(int(i_nlayers)))
        color_rs=cmap_rs(mynorm(int(i_nlayers)))
        if(include_fs_only):
            label="FS" if i==len(nlayers)-1 else None
            ax.plot(pba_fs.GRB.get_lc_times(spec=True)/cgs.day,
                    pba_fs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='-', color=color_fs, label=label)
        label="RS" if i==len(nlayers)-1 else None
        ax.plot(pba_fsrs.GRB.get_lc_times(spec=True)/cgs.day,
                pba_fsrs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='--', color=color_rs, label=label, lw=2)
    ax.legend()


    ax.set_ylabel("Total Emissivity", fontsize=12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # axes[iv_n].legend()
    ax.set_xlim(1e-2,1e7)
    ax.grid()
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("tburst" + " [day]",fontsize=12)
    if legend: plt.legend()
    plt.colorbar(cm.ScalarMappable(norm=mynorm, cmap=plt.cm.Reds),label="Layers")
    plt.tight_layout()
    plt.title("Comoving Spectrum Gaussian Jet")
    plt.tight_layout()
    plt.savefig(workdir+figname, dpi=256)
    plt.show()

if __name__ == '__main__':
    # plot_ejecta_layers(ishells=(0,), ilayers=(0,10,20,40,49), nlayers_a=50,
    #                    v_n_x = "tburst", v_n_ys = ("B", "B_rs", "gamma_min", "gamma_min_rs","gamma_c", "gamma_c_rs"), colors_by="layers",legend=False,run_fs_only=True,
    #                    figname="dyn_layers_fs.png")
    plot_ejecta_layers_spec(freq=1e9, ishells=(0,), ilayers=(0,10,20,30,40,49),nlayers_a=50,type="a",method_eats="adaptive")
    #
    # plot_tst_total_spec_resolution(freq=1e9, nlayers=(40,80,120,160,200), type="pw",method_eats="piece-wise")
    plot_tst_total_spec_resolution(freq=1e9, nlayers=(10,30,50,70),type="a",method_eats="adaptive", include_fs_only=True)
    exit(0)