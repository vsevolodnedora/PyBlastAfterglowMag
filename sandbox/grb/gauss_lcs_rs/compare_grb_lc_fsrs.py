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

import package.src.PyBlastAfterglowMag as PBA

try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

pars = {"Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.4,
        "nlayers_pw": 0, "nlayers_a": 0, "struct":"gaussian"}

# pars =  {"struct":"gaussian",
#  "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
#  "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10}

class RefData():
    def __init__(self,workdir:str,fname:str):
        self.workdir = workdir
        self.keys = ["tburst",
                     "tcomov",
                     "Gamma",
                     "Eint2",
                     "Eint3",
                     "theta",
                     "Erad2",
                     "Erad3",
                     "Esh2",
                     "Esh3",
                     "Ead2",
                     "Ead3",
                     "M2",
                     "M3",
                     "deltaR4"]
        self.keys = [
            "R", "rho", "dlnrhodr", "Gamma", "Eint2", "U_e", "theta", "ctheta", "Erad2", "Esh2", "Ead2", "M2",
            "tcomov", "tburst", "tt", "delta", "W",
            "Gamma43", "Eint3", "U_e3", "Erad3", "Esh3", "Ead3", "M3", "rho4", "deltaR4", "W3"
        ]
        self.refdata = None
        self.fname = fname
    def idx(self, key : str) -> int:
        return self.keys.index(key)
    def load(self) -> None:
        # self.refdata = np.loadtxt(self.workdir+"reference_fsrs.txt")
        self.refdata = h5py.File(self.workdir+self.fname,'r')
    def get(self,freq:float,theta_obs:float) -> [np.ndarray, np.ndarray]:
        if (self.refdata is None): self.load()
        group = self.refdata["{:.2f}deg {:.2e}Hz".format(freq, theta_obs)]
        times = np.array(group["time"])
        fluxdens = np.array(group["fluxdens"])
        return (times, fluxdens)
    def close(self) -> None:
        self.refdata.close()

structure = {"Eiso_c": 1.e53, "Gamma0c": 1000., "M0c": -1., "theta_c": 0.1, "theta_w": 0.3,
             "nlayers_pw": 50, "nlayers_a": 10, "struct": "gaussian"}

def plot_ejecta_layers(ishells=(0,), ilayers=(0,25,49),
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

    pba_id = PBA.id_maker_analytic.JetStruct(n_layers_pw=50,n_layers_a=10)
    pba_id_dict, id_pars = pba_id.get_1D_id(type="adaptive",pars=structure)
    pba_id.save_1d_id(id_dict=pba_id_dict, id_pars=id_pars, outfpath=workdir+"tophat_grb_id.h5")


        ### run fs-only model
    if(run_fs_only):
        PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                               newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                        "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_tophad.h5", "method_eats": "adaptive"},
                               parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
        pba_fs = PBA.interface.PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
        pba_fs.run(loglevel="info")

    # run fsrs model
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_spectrum":"lc_grb_tophad.h5", "method_eats": "adaptive"},
                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
    pba_fsrs = PBA.interface.PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
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

    fid, axes = plt.subplots(ncols=1, nrows=len(v_n_ys), figsize=(4.6,4.2),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    for iv_n, v_n in enumerate(v_n_ys):
        ax = axes[iv_n] if len(v_n_ys) > 1 else axes
        i = 0
        if(run_fs_only):
            for il, layer in enumerate(layers):
                x_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n_x,ishell=ishells[0],ilayer=il)#  np.array(dfile[layer][v_n_x])
                y_arr = pba_fs.GRB.get_dyn_arr(v_n=v_n,ishell=ishells[0],ilayer=il)#np.array(dfile[layer][v_n])
                y_arr = y_arr[x_arr > 0]
                x_arr = x_arr[x_arr > 0]
                # if (v_n == "R"):
                # y_arr = y_arr/y_arr.max()
                if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
                else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
                if (v_n_x == "tburst"): x_arr /=PBA.utils.cgs.day;
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
            if (v_n_x == "tburst"): x_arr /=PBA.utils.cgs.day;
            ax.plot(x_arr, y_arr, ls='--', color=color, label=layer)
            i=i+1

            # --- plot ref data
            # x_arr = ref.get(il, v_n_x)
            # if (v_n_x == "tburst"): x_arr /=cgs.day;
            # ax.plot(x_arr, ref.get(il, v_n), ls=':', color=color, lw=2.)



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

def plot_ejecta_layers_spec(freq=1e18,ishells=(0,), ilayers=(0,25,49),colors_by="layers",legend=False,
                       figname="dyn_layers_fsrs.png", run_fs_only=True):

    # prepare ID
    workdir = os.getcwd()+"/"
    # prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
    #                       "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},
    #                      type="pw",outfpath="tophat_grb_id.h5")
    theta_h = np.pi/2.
    one_min_cos = 2. * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h)
    ang_size_layer = 2.0 * np.pi * one_min_cos / (4.0 * np.pi)
    prepare_grb_ej_id_1d({"Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.3,
                          "nlayers_pw": 100, "nlayers_a": 40, "struct":"gaussian"},
                           type="pw",outfpath="tophat_grb_id.h5")

    ### run fs-only model
    if(run_fs_only):
        modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                               newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                        "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_tophad.h5",
                                        "fname_spectrum_layers":"spec_fs_layers.h5", "method_eats": "adaptive"},
                               parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
        pba_fs = PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
        pba_fs.run(loglevel="info")

    # run fsrs model
    modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_spectrum":"lc_grb_tophad.h5",
                                    "fname_spectrum_layers":"spec_fsrs_layers.h5", "method_eats": "piece-wise"},
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

    fid, axes = plt.subplots(ncols=1, nrows=2, figsize=(4.6,4.2),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    ax = axes[0]
    i = 0
    if(run_fs_only):
        for il, layer in enumerate(layers):
            x_arr = pba_fs.GRB.get_lc_times(spec=True)#  np.array(dfile[layer][v_n_x])
            y_arr = pba_fs.GRB.get_lc(freq=freq,ishell=ishells[0],ilayer=il,spec=True)#np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            x_arr /=cgs.day
            ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
            i=i+1
    # --------------------------------
    i = 0
    for il, layer in enumerate(layers):
        x_arr = pba_fsrs.GRB.get_lc_times(spec=True)#  np.array(dfile[layer][v_n_x])
        y_arr = pba_fsrs.GRB.get_lc(freq=freq,ishell=ishells[0],ilayer=il,spec=True)#np.array(dfile[layer][v_n])
        y_arr = y_arr[x_arr > 0]
        x_arr = x_arr[x_arr > 0]
        # if (v_n == "R"):
        # y_arr = y_arr/y_arr.max()
        if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
        else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
        x_arr /=cgs.day
        ax.plot(x_arr, y_arr, ls='--', color=color, label=layer)
        i=i+1

        # --- plot ref data
        # x_arr = ref.get(il, v_n_x)
        # if (v_n_x == "tburst"): x_arr /=cgs.day;
        # ax.plot(x_arr, ref.get(il, v_n), ls=':', color=color, lw=2.)

    ax = axes[1]
    if(run_fs_only):
        ax.plot(pba_fs.GRB.get_lc_times(spec=True)/cgs.day,
                pba_fs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='-', color='gray', label='FS')
    ax.plot(pba_fsrs.GRB.get_lc_times(spec=True)/cgs.day,
            pba_fsrs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=True), ls='--', color='black', label='FSRS')

    for ax in axes:
        ax.set_ylabel("Emissivity", fontsize=12)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.set_xlim(1e-2,1e7)
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")


    ax.set_xlabel("tburst" + " [day]",fontsize=12)
    if legend: plt.legend()
    plt.tight_layout()
    plt.savefig(workdir+figname, dpi=256)
    plt.show()

def plot_tst_total_spec_resolution(freq=1e9, nlayers=(10,20,40,80,120),legend=False,
                                   figname="dyn_layers_fsrs.png", run_fs_only=True,run_fsrs_only=False,
                                   plot_ref=True,plot_afgpy=True,
                                   method_eats="piece-wise"):
    workdir = os.getcwd()+"/"

    fid, ax = plt.subplots(ncols=1, nrows=1, figsize=(4.6,4.2),sharex="all")
    cmap_fs = cm.Blues_r
    cmap_fsrs = cm.Reds_r
    mynorm = Normalize(vmin=0,vmax=len(nlayers))#norm(len(ishells)*len(ilayers))

    ref = RefData(workdir=workdir,fname="reference_lc_spread.h5")

    methods_spread={"methods":["AFGPY"],#,"Adi","AFGPY","AA"],
                    "ls":["-"],#,"--",":"]
                    }
    for (method_spread, ls) in zip(methods_spread["methods"], methods_spread["ls"]):
        for (theta, lw, color) in zip([0, 0.9, 1.5], [1, 1, 1], ["blue", "green", "red"]):
            for i, i_nlayers in enumerate(nlayers):
                pars["nlayers_a"] = i_nlayers
                pars["nlayers_pw"] = i_nlayers

                pba_id = PBA.id_maker_analytic.JetStruct(n_layers_pw=i_nlayers,n_layers_a=i_nlayers)
                pba_id_dict, id_pars = pba_id.get_1D_id(type=method_eats,pars=structure)
                pba_id.save_1d_id(id_dict=pba_id_dict, id_pars=id_pars, outfpath=workdir+"gauss_grb_id.h5")


                PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="main", newpars={"theta_obs":theta}, newopts={},
                                       parfile="default_parfile.par", newparfile="default_parfile.par",keep_old=False)
                ### run fs-only model
                if(run_fs_only):
                    PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                                           newopts={"rhs_type":"grb_fs", "outfpath":"grb_fs.h5", "do_rs":"no",
                                                    "fname_ejecta_id":"gauss_grb_id.h5", "method_spread":method_spread,
                                                    "fname_dyn":"dyn_bw_fs.h5","fname_light_curve":"lc_grb_fs.h5",
                                                    "fname_light_curve_layers":"lc_grb_fs_layers.h5", "method_eats": method_eats,
                                                    # "method_comp_mode": "observFlux", "do_spec":"no"
                                                    },
                                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
                    pba_fs = PBA.interface.PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
                    pba_fs.run(loglevel="info")

                # run fsrs model
                if (run_fsrs_only):
                    PBA.parfile_tools.modify_parfile_par_opt(workingdir=workdir, part="grb", newpars={},
                                           newopts={"rhs_type":"grb_fsrs", "outfpath":"grb_fsrs.h5", "do_rs":"yes",
                                                    "fname_ejecta_id":"gauss_grb_id.h5", "method_spread":method_spread,
                                                    "fname_dyn":"dyn_bw_fsrs.h5","fname_light_curve":"lc_grb_fsrs.h5",
                                                    "fname_light_curve_layers":"lc_grb_fsrs_layers.h5", "method_eats": method_eats
                                                    },
                                           parfile="default_parfile.par", newparfile="parfile.par",keep_old=True)
                    pba_fsrs = PBA.interface.PyBlastAfterglow(workingdir=workdir, readparfileforpaths=True, parfile="parfile.par")
                    pba_fsrs.run(loglevel="info")

                color_fs = cmap_fs(mynorm(int(i))) if color is None else color
                color_fsrs = cmap_fsrs(mynorm(int(i))) if color is None else color

                if(run_fs_only):
                    ax.plot(pba_fs.GRB.get_lc_times(spec=False)/PBA.utils.cgs.day,
                            # pba_fs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=False),
                    pba_fs.GRB.get_lc_totalflux(freq=freq,time=None,spec=False),
                    ls=ls, color=color_fs, lw=lw, label='FS')

                if (run_fsrs_only):
                    ax.plot(pba_fsrs.GRB.get_lc_times(spec=False)/PBA.utils.cgs.day,
                            # pba_fsrs.GRB.get_lc(freq=freq,ishell=None,ilayer=None,spec=False),
                            pba_fsrs.GRB.get_lc_totalflux(freq=freq,spec=False),
                            ls=ls, color=color_fsrs, lw=lw, label='FSRS')

                if (plot_ref):
                    times, fluxes = ref.get(freq=freq, theta_obs=theta)
                    ax.plot(times/PBA.utils.cgs.day, fluxes, color='gray', lw=0.8, ls=':')
                if (plot_afgpy):
                    # -------- Afterglopy --------------
                    Z = {'jetType':     grb.jet.TopHat if structure["struct"] == "tophat" else grb.jet.Gaussian,     # Top-Hat jet
                         'specType':    0,                  # Basic Synchrotron Spectrum
                         'counterjet':  1,
                         'spread':      7,
                         'thetaObs':    theta,   # Viewing angle in radians
                         'E0':          structure["Eiso_c"], # Isotropic-equivalent energy in erg
                         'g0':          structure["Gamma0c"],
                         'thetaCore':   structure["theta_c"],    # Half-opening angle in radians
                         'thetaWing':   structure["theta_w"],
                         'n0':          pba_fs.main_pars["n_ism"],    # circumburst density in cm^{-3}
                         'p':           pba_fs.GRB.pars["p"],    # electron energy distribution index
                         'epsilon_e':   pba_fs.GRB.pars["eps_e"],    # epsilon_e
                         'epsilon_B':   pba_fs.GRB.pars["eps_b"],   # epsilon_B
                         'xi_N':        1.0,    # Fraction of electrons accelerated
                         'd_L':         pba_fs.main_pars["d_l"], # Luminosity distance in cm
                         'z':           pba_fs.main_pars["z"]}   # redshift

                    t = np.geomspace(1.0 * 86400.0, 1.0e5 * 86400.0, 100)
                    nu = np.empty(t.shape)
                    nu[:] = freq
                    Fnu = grb.fluxDensity(t, nu, **Z)
                    ax.plot(t/PBA.utils.cgs.day, Fnu, color='gray', lw=0.8, ls='-')


    ax.set_ylabel("Flux Density [mJy]", fontsize=12)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # axes[iv_n].legend()
    ax.set_xlim(1e-2,1e7)
    ax.grid()
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("time" + " [day]",fontsize=12)
    if legend: plt.legend()
    plt.tight_layout()
    plt.savefig(workdir+figname, dpi=256)
    plt.show()

if __name__ == '__main__':

    plot_tst_total_spec_resolution(freq=1e9, nlayers=(20,),method_eats="adaptive")#30,50,70


    # plot_ejecta_layers(ishells=(0,), ilayers=(0,),
    #                    v_n_x = "tburst", v_n_ys = ("B", "B_rs", "gamma_min", "gamma_min_rs","gamma_c", "gamma_c_rs"), colors_by="layers",legend=False,
    #                    figname="dyn_layers_fs.png")
    # plot_ejecta_layers_spec(freq=1e9, ishells=(0,), ilayers=(0,2,4,6,9))
    #
    # plot_tst_total_spec_resolution(freq=1e9, nlayers=(80,120,160,200), type="pw",method_eats="piece-wise")
    exit(0)