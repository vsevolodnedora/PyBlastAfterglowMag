# import PyBlastAfterglowMag
import numpy as np
import h5py
import shutil
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.cm import ScalarMappable
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

# from PyBlastAfterglowMag import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag

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
parfiledir = os.getcwd().replace("structured","")
paperfigdir = "/home/vsevolod/Work/GIT/overleaf/grb_afg/figs/"


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
    def get(self,il,key) -> np.ndarray:
        if (self.refdata is None):
            self.load()
        if (key in self.keys):
            return np.array(np.array(self.refdata[f"layer={il}"][key]))
        elif (key=="mom"):
            return np.array(np.array(self.refdata[f"layer={il}"]["Gamma"])) * \
                get_beta(np.array(np.array(self.refdata[f"layer={il}"]["Gamma"])))
        else:
            raise KeyError(f"key={key} is not recognized")
    def close(self) -> None:
        self.refdata.close()

    def nlayers(self):
        if (self.refdata is None):
            self.load()
        return int(len(self.refdata.keys()))

class RefDataLC():
    def __init__(self,workdir:str,fname:str):
        self.workdir = workdir
        self.keys = ["times", "fluxes"]
        self.refdata = None
        self.fname = fname
    def idx(self, key : str) -> int:
        return self.keys.index(key)
    def load(self) -> None:
        # self.refdata = np.loadtxt(self.workdir+"reference_fsrs.txt")
        self.refdata = h5py.File(self.workdir+self.fname,'r')
    def get(self,il) -> np.ndarray:
        if (il == -1):
            arr = np.zeros_like(np.array(self.refdata[f"layer={0}"]["fluxes"]))
            for _il in range(self.nlayers()):
                arr += np.array(self.refdata[f"layer={_il}"]["fluxes"])
            return arr
        if (self.refdata is None):
            self.load()
        return np.array(self.refdata[f"layer={il}"]["fluxes"])

    def close(self) -> None:
        self.refdata.close()

    def times(self):
        if (self.refdata is None):
            self.load()
        return np.array(np.array(self.refdata[f"layer={0}"]["times"]))

    def nlayers(self):
        if (self.refdata is None):
            self.load()
        return int(len(self.refdata.keys()))

class TestBases():

    def run_pw(self, pars : dict, struct : dict, opts : dict, opts_grb : dict) -> PyBlastAfterglow:
        workdir = os.getcwd()+'/'
        # copy the main parfile into workdir
        shutil.copy(parfiledir+"parfile_def.par", curdir+"parfile_def.par")
        prepare_grb_ej_id_1d(struct, type="pw", outfpath=workdir+"gauss_grb_id.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",
                               newpars=pars,
                               newopts=opts,
                               parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars=pars,
                               newopts=opts_grb,
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)

        pba_pw = PyBlastAfterglow(workingdir=os.curdir+'/',readparfileforpaths=True, parfile="parfile.par")
        pba_pw.run()
        # remove the default parfile from the working dir
        os.remove(curdir+"parfile_def.par")
        return pba_pw

    def run_a(self, pars : dict, struct : dict, opts : dict, opts_grb : dict) -> PyBlastAfterglow:
        workdir = os.getcwd()+'/'
        # copy the main parfile into workdir
        shutil.copy(parfiledir+"parfile_def.par", curdir+"parfile_def.par")
        prepare_grb_ej_id_1d(struct, type="a", outfpath=workdir+"gauss_grb_id.h5")
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",
                               newpars=pars,
                               newopts=opts,
                               parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
        modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                               newpars=pars,
                               newopts=opts_grb,
                               parfile="parfile.par", newparfile="parfile.par", keep_old=False)

        pba_a = PyBlastAfterglow(workingdir=os.curdir+'/',readparfileforpaths=True, parfile="parfile.par")
        pba_a.run()
        # remove the default parfile from the working dir
        os.remove(curdir+"parfile_def.par")
        return pba_a

    def plot_lcs(self, ax, pars : dict, pba : PyBlastAfterglow, layers = (), plot={}, plot_layer={}):
        ls = plot["ls"] if "ls" in plot.keys() else '-'
        color = plot["color"] if "color" in plot.keys() else 'red'
        alpha = plot["alpha"] if "alpha" in plot.keys() else 1.0
        label = plot["label"] if "label" in plot.keys() else None
        ax.plot(pba.GRB.get_lc_times(), pba.GRB.get_lc_totalflux(freq=pars["obs_freq"]),
                ls=ls, color=color, label=label, alpha=alpha)
        cmap = plot_layer["cmap"] if "cmap" in plot_layer.keys() else 'viridis'
        cmap = cm.get_cmap(cmap)
        vmin = plot_layer["vmin"] if "vmin" in plot_layer.keys() else 0
        vmax = plot_layer["vmax"] if "vmax" in plot_layer.keys() else int(pba.GRB.get_lc_obj().attrs["nlayers"])
        norm = Normalize(vmin=vmin, vmax=vmax)
        for il in range(int(pba.GRB.get_lc_obj().attrs["nlayers"])):
            if (((len(layers) > 0) and (il in layers)) or (len(layers)==0)):
                ls = plot_layer["ls"] if "ls" in plot_layer.keys() else '-'
                color = plot_layer["color"] if "color" in plot_layer.keys() else 'red'
                color = cmap(norm(il)) if "cmap" in plot_layer.keys() else color
                alpha = plot_layer["alpha"] if "alpha" in plot_layer.keys() else 0.9
                label = plot_layer["label"] if "label" in plot_layer.keys() else None
                ax.plot(pba.GRB.get_lc_times(), pba.GRB.get_lc(freq=pars["obs_freq"], ishell=0, ilayer=il),
                        ls=ls, color=color, alpha=alpha, label=label)

    def plot_dyn(self, ax, pba : PyBlastAfterglow, v_n_x : str, v_n_y : str, layers=(), plot_layer = {}):
        cmap = plot_layer["cmap"] if "cmap" in plot_layer.keys() else 'viridis'
        cmap = cm.get_cmap(cmap)
        llayers = pba.GRB.get_dyn_obj().attrs["nlayers"]
        vmin = plot_layer["vmin"] if "vmin" in plot_layer.keys() else 0
        vmax = plot_layer["vmax"] if "vmax" in plot_layer.keys() else int(pba.GRB.get_dyn_obj().attrs["nlayers"])
        norm = Normalize(vmin=vmin, vmax=vmax)
        for il in range(int(pba.GRB.get_dyn_obj().attrs["nlayers"])):
            if (((len(layers) > 0) and (il in layers)) or (len(layers)==0)):
                ls = plot_layer["ls"] if "ls" in plot_layer.keys() else '-'
                color = plot_layer["color"] if "color" in plot_layer.keys() else 'red'
                color = cmap(norm(il)) if "cmap" in plot_layer.keys() else color
                alpha = plot_layer["alpha"] if "alpha" in plot_layer.keys() else 0.9
                label = plot_layer["label"] if "label" in plot_layer.keys() else None
                ax.plot(pba.GRB.get_dyn_arr(v_n=v_n_x,ishell=0,ilayer=il),
                        pba.GRB.get_dyn_arr(v_n=v_n_y,ishell=0,ilayer=il), ls=ls, color=color, alpha=alpha, label=label)

class TestCasesFS(TestBases):
    name_grb170817a_like_event = "Gauss 170817-like off-axis"
    struct_grb170817a_like_event = {
        "struct":"gaussian",
        "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
        "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 10
    }
    pars_grb170817a_like_event = {
        "obs_freq":3e9,
        "n_ism":0.00031,    "eps_e":0.0708,
        "d_l":1.27e+26,     "eps_b": 0.0052,
        "z": 0.0099,        "p":2.16,
        "theta_obs": 0.3752
    }

    name_1 = "Gauss off-axis"
    struct_1 = {
        "struct":"gaussian",
        "Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,
        "theta_c": 0.1, "theta_w": 0.4, "nlayers_pw": 150, "nlayers_a": 52
    }
    pars_1 = {
        "obs_freq":3e9,
        "n_ism": 1e-2,      "eps_e": 0.1,
        "d_l": 3.09e26,     "eps_b": 0.01,
        "z": 0.028,         "p": 2.2,
        "theta_obs": 0.9
    }
    opts_a_1 = {"method_eats":"adaptive", "method_spread":"AFGPY"}
    opts_pw_1 = {"method_eats":"piece-wise", "method_spread":"AA"}

    def __init__(self):
        pass

    def plot_lcs_ref(self, ax, ref : RefDataLC, nlayers : int, layers = (),):
        cmap = cm.get_cmap('Greys')
        norm = Normalize(vmin=-50,vmax=60)
        nlayers = ref.nlayers()
        nlayers_a = int(nlayers)
        if (int(nlayers)!=int(nlayers_a)):
            raise ValueError(f"For ref deta expected nlayers={nlayers_a} got nlayer={nlayers}")
        ax.plot(ref.times(), ref.get(-1), ls='-', color='black', lw=0.9, zorder=-1, label=r"\texttt{afterglowpy}")
        for il in range(int(nlayers)):
            if (((len(layers) > 0) and (il in layers)) or (len(layers)==0)):
                # --- plot ref data
                x_arr = ref.times()
                ax.plot(x_arr, ref.get(il), ls='-.', color=cmap(norm(il)), lw=.8, zorder=-1)

    def plot_dyn_ref(self, ax, ref : RefData, v_n_x : str, v_n_y : str, nlayers : int, layers=()):
        cmap = cm.get_cmap('Greys')
        norm = Normalize(vmin=-50,vmax=60)
        nlayers = ref.nlayers()
        nlayers_a = int(nlayers)
        if (int(nlayers)!=int(nlayers_a)):
            raise ValueError(f"For ref deta expected nlayers={nlayers_a} got nlayer={nlayers}")
        for il in range(int(nlayers)):
            if (((len(layers) > 0) and (il in layers)) or (len(layers)==0)):
                # --- plot ref data
                x_arr = ref.get(il, v_n_x)
                ax.plot(x_arr, ref.get(il, v_n_y), ls=':', color=cmap(norm(il)), lw=1., zorder=-1)

    def plot_170817_like(self, pars : dict, struct : dict, title : str):

        # run the code for given pars
        pba_pw = self.run_pw(struct=struct, pars=pars, opts={}, opts_grb=self.opts_pw_1)
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)

        fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

        # self.plot_lcs(ax, pars=pars, pba_pw=pba_pw, pba_a=pba_a, ref = None)
        self.plot_lcs(ax=ax, pars=pars, pba=pba_pw, layers = (),
                      plot={"ls":'-', "color":"green", "label":"PBA [PW]"},
                      plot_layer={"ls":'-', "cmap":"inferno", "alpha":.5})
        self.plot_lcs(ax=ax, pars=pars, pba=pba_a, layers = (),
                      plot={"ls":'-', "color":"red", "label":"PBA [A]"},
                      plot_layer={"ls":'-', "cmap":"viridis", "alpha":.5})

        # --------- Reference Models -----------
        tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"./afgpy_grb170817.txt",unpack=True)
        ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

        tts, fluxes = np.loadtxt(curdir+"./jelib_grb170817.txt",unpack=True)
        ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')

        # -------- Afterglopy --------------
        Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
             'specType':    0,                  # Basic Synchrotron Spectrum
             'counterjet':  1,
             'spread':      7,
             'thetaObs':    pars["theta_obs"],   # Viewing angle in radians
             'E0':          struct["Eiso_c"], # Isotropic-equivalent energy in erg
             'g0':          struct["Gamma0c"],
             'thetaCore':   struct["theta_c"],    # Half-opening angle in radians
             'thetaWing':   struct["theta_w"],
             'n0':          pars["n_ism"],    # circumburst density in cm^{-3}
             'p':           pars["p"],    # electron energy distribution index
             'epsilon_e':   pars["eps_e"],    # epsilon_e
             'epsilon_B':   pars["eps_b"],   # epsilon_B
             'xi_N':        1.0,    # Fraction of electrons accelerated
             'd_L':         pars["d_l"], # Luminosity distance in cm
             'z':           pars["z"]}   # redshift

        t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
        nu = np.empty(t.shape)
        nu[:] = 3.0e9
        Fnu = grb.fluxDensity(t, nu, **Z)

        # plot
        ax.plot(t, Fnu, ls='-', color='gray', label='afterglopy')

        ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("Flux density [mJy]", fontsize=12)
        ax.set_title(title)
        ax.set_xlim(1e5,1e8)
        ax.set_ylim(1e-4,1)
        ax.grid()
        plt.show()

    def plot_generic(self, pars : dict, struct : dict, title : str):
        # run the code for given pars
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)

        ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_layer.h5")
        ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_layer.h5")

        fig, axes = plt.subplots(figsize=(9,9.5), ncols=1, nrows=3)

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = (),
                      plot={"ls":'-', "color":"red", "label":"PBA [A]"},
                      plot_layer={"ls":'-', "cmap":"viridis", "alpha":.5})
        self.plot_lcs_ref(ax=axes[0], ref=ref_lc, nlayers=pba_a.GRB.get_lc_obj().attrs["nlayers"])

        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=(), plot_layer={"ls":'-', "cmap":"viridis", "alpha":.5})
        self.plot_dyn_ref(axes[1], ref=ref, v_n_x="tburst", v_n_y="mom", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"])

        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=(), plot_layer={"ls":'-', "cmap":"viridis", "alpha":.5})
        self.plot_dyn_ref(axes[2], ref=ref, v_n_x="tburst", v_n_y="theta", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"])


        # -------- Afterglopy --------------
        Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
             'specType':    0,                  # Basic Synchrotron Spectrum
             'counterjet':  1,
             'spread':      7,
             'thetaObs':    pars["theta_obs"],   # Viewing angle in radians
             'E0':          struct["Eiso_c"], # Isotropic-equivalent energy in erg
             'g0':          struct["Gamma0c"],
             'thetaCore':   struct["theta_c"],    # Half-opening angle in radians
             'thetaWing':   struct["theta_w"],
             'n0':          pars["n_ism"],    # circumburst density in cm^{-3}
             'p':           pars["p"],    # electron energy distribution index
             'epsilon_e':   pars["eps_e"],    # epsilon_e
             'epsilon_B':   pars["eps_b"],   # epsilon_B
             'xi_N':        1.0,    # Fraction of electrons accelerated
             'd_L':         pars["d_l"], # Luminosity distance in cm
             'z':           pars["z"]}   # redshift

        t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
        nu = np.empty(t.shape)
        nu[:] = pars["obs_freq"]
        Fnu = grb.fluxDensity(t, nu, **Z)

        # plot
        ax = axes[0]
        ax.plot(t, Fnu, ls='-', color='gray', label='afterglopy')
        ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("Flux density [mJy]", fontsize=12)
        ax.set_title(title)
        ax.set_xlim(1e5,1e8)
        ax.set_ylim(1e-4,1e0)
        ax.grid()

        ax = axes[1]
        ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")

        ax = axes[2]
        ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("linear")

        plt.show()

    def paper_plot_compare_spreading(self, pars : dict, struct : dict, title : str):
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_layer.h5")

        ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_layer_GamInf.h5")
        ref = RefData(workdir=curdir, fname="reference_afgpy_dyn_GamInf.h5")

        fig, axes = plt.subplots(figsize=(6,6.5), ncols=1, nrows=3)

        pars["Gamma0_frac_when_start_spread"] = .1
        self.opts_a_1["method_limit_spread"] = "Gamma0Frac"

        layers=(0,45)

        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)
        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"red", "label":r"R19"},
                      plot_layer={"ls":'-.', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        self.opts_a_1["method_spread"] = "AA"
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"blue", "label":"GP12"},
                      plot_layer={"ls":'-.', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        self.opts_a_1["method_spread"] = "Adi"
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"green", "label":"HDL99"},
                      plot_layer={"ls":'-.', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'--', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'--', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------
        self.plot_lcs_ref(ax=axes[0], ref=ref_lc, nlayers=pba_a.GRB.get_lc_obj().attrs["nlayers"], layers=layers)
        self.plot_dyn_ref(axes[1], ref=ref, v_n_x="tburst", v_n_y="mom", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"], layers=layers)
        self.plot_dyn_ref(axes[2], ref=ref, v_n_x="tburst", v_n_y="theta", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"], layers=layers)

        # -------- Afterglopy --------------
        # Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
        #      'specType':    0,                  # Basic Synchrotron Spectrum
        #      'counterjet':  1,
        #      'spread':      7,
        #      'thetaObs':    pars["theta_obs"],   # Viewing angle in radians
        #      'E0':          struct["Eiso_c"], # Isotropic-equivalent energy in erg
        #      'g0':          struct["Gamma0c"],
        #      'thetaCore':   struct["theta_c"],    # Half-opening angle in radians
        #      'thetaWing':   struct["theta_w"],
        #      'n0':          pars["n_ism"],    # circumburst density in cm^{-3}
        #      'p':           pars["p"],    # electron energy distribution index
        #      'epsilon_e':   pars["eps_e"],    # epsilon_e
        #      'epsilon_B':   pars["eps_b"],   # epsilon_B
        #      'xi_N':        1.0,    # Fraction of electrons accelerated
        #      'd_L':         pars["d_l"], # Luminosity distance in cm
        #      'z':           pars["z"]}   # redshift
        #
        # t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
        # nu = np.empty(t.shape)
        # nu[:] = pars["obs_freq"]
        # Fnu = grb.fluxDensity(t, nu, **Z)
        #
        # axes[0].plot(t, Fnu, ls='-', color='gray', label='afterglopy')


        # plot
        ax = axes[0]
        # ax.grid()
        ax.legend(fancybox=True, loc='upper left',
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=1, fontsize=12,
                  framealpha=0., borderaxespad=0.)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("Flux density [mJy]", fontsize=12)
        # ax.set_title(title)
        ax.set_xlim(1e6,1e8)
        ax.set_ylim(1e-4,1e0)
        # ax.grid()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12,
                        direction='in',
                        bottom=True, top=True, left=True, right=True)
        # ax.set_facecolor("pink")

        ax = axes[1]
        # ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(r"$\Gamma\beta$", fontsize=12)
        ax.set_xlim(1e6,1e11)
        ax.set_ylim(1e-3,1e5)
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        # ax.set_facecolor("pink")

        ax = axes[2]
        # ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("linear")
        ax.set_ylabel(r"$\omega$ [rad]", fontsize=12)
        ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
        ax.set_xlim(1e6,1e11)
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        # ax.set_facecolor("pink")

        plt.tight_layout()
        print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        plt.show()

class TestCasesRS(TestBases):
    name_1 = "Gauss off-axis"
    struct_1 = {
        "struct":"gaussian",
        "Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,
        "theta_c": 0.1, "theta_w": 0.4, "nlayers_pw": 150, "nlayers_a": 52
    }
    pars_1 = {
        "obs_freq":3e9,
        "n_ism": 1e-2,      "eps_e": 0.1,   "eps_e_rs": 0.2,
        "d_l": 3.09e26,     "eps_b": 0.01,  "eps_b_rs": 0.02,
        "z": 0.028,         "p": 2.2,       "p_rs": 2.4,
        "theta_obs": 0.9,
    }
    opts_1 = {"rtol": 1e-6, "ntb":10000, "iout": 10}
    opts_a_1 = {"method_eats":"adaptive", "method_spread":"AFGPY", "rhs_type":"grb_fs", "do_rs": "no"}

    def __init__(self):
        super().__init__()

    def paper_plot_compare_fsrs(self, pars : dict, struct : dict, title : str):
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_layer.h5")

        fig, axes = plt.subplots(figsize=(6,6.5), ncols=1, nrows=3)

        layers=(0,45)

        pba_a = self.run_a(struct=struct, pars=pars, opts=self.opts_1, opts_grb=self.opts_a_1)
        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"red", "label":r"R19"},
                      plot_layer={"ls":'-.', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        self.opts_a_1["rhs_type"] = "grb_fsrs"
        self.opts_a_1["do_rs"] = "yes"
        pba_a = self.run_a(struct=struct, pars=pars, opts=self.opts_1, opts_grb=self.opts_a_1)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"blue", "label":"GP12"},
                      plot_layer={"ls":'-.', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})


        # -------- Afterglopy --------------
        Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
             'specType':    0,                  # Basic Synchrotron Spectrum
             'counterjet':  1,
             'spread':      7,
             'thetaObs':    pars["theta_obs"],   # Viewing angle in radians
             'E0':          struct["Eiso_c"], # Isotropic-equivalent energy in erg
             'g0':          struct["Gamma0c"],
             'thetaCore':   struct["theta_c"],    # Half-opening angle in radians
             'thetaWing':   struct["theta_w"],
             'n0':          pars["n_ism"],    # circumburst density in cm^{-3}
             'p':           pars["p"],    # electron energy distribution index
             'epsilon_e':   pars["eps_e"],    # epsilon_e
             'epsilon_B':   pars["eps_b"],   # epsilon_B
             'xi_N':        1.0,    # Fraction of electrons accelerated
             'd_L':         pars["d_l"], # Luminosity distance in cm
             'z':           pars["z"]}   # redshift

        t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
        nu = np.empty(t.shape)
        nu[:] = pars["obs_freq"]
        Fnu = grb.fluxDensity(t, nu, **Z)

        axes[0].plot(t, Fnu, ls='-', color='gray', label='afterglopy')


        # plot
        ax = axes[0]
        # ax.grid()
        ax.legend(fancybox=True, loc='upper left',
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=1, fontsize=12,
                  framealpha=0., borderaxespad=0.)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("Flux density [mJy]", fontsize=12)
        # ax.set_title(title)
        ax.set_xlim(1e6,1e8)
        ax.set_ylim(1e-4,1e0)
        # ax.grid()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.set_facecolor("pink")

        ax = axes[1]
        # ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(r"$\Gamma\beta$", fontsize=12)
        ax.set_xlim(1e6,1e11)
        ax.set_ylim(1e-3,1e5)
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        # ax.set_facecolor("pink")

        ax = axes[2]
        # ax.grid()
        ax.legend()
        ax.set_xscale("log")
        ax.set_yscale("linear")
        ax.set_ylabel(r"$\omega$ [rad]", fontsize=12)
        ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
        ax.set_xlim(1e6,1e11)
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        # ax.set_facecolor("pink")

        plt.tight_layout()
        print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        plt.show()

def main():
    grb = TestCasesFS()
    # grb.plot_170817_like(struct=grb.struct_grb170817a_like_event, pars=grb.pars_grb170817a_like_event, title=grb.name_grb170817a_like_event)
    # grb.plot_generic(struct=grb.struct_1, pars=grb.pars_1, title=grb.name_1)
    # grb.plot_generic(struct=grb.struct_1, pars=grb.pars_1, title=grb.name_1)
    # grb.paper_plot_compare_spreading(struct=grb.struct_1, pars=grb.pars_1, title=grb.name_1)

    grbrs = TestCasesRS()
    grbrs.paper_plot_compare_fsrs(struct=grbrs.struct_1, pars=grbrs.pars_1, title=grbrs.name_1)

if __name__ == '__main__':
    main()




def tst_tophat_old():
    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    fname = "test_data_grb/lcs_tophat.h5"

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    # colnames, data = load_data(curdir + fname)
    # nx, ny = data.shape
    # t_arr = data[:, 0]
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    ax.plot(our_times, our_fluxes, ls='-', color='black', label='PyBlastAfterglow [pw]')
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()

def compareTopHatLightCurves():

    ''' Run 'compareTopHatLightCurves()' first '''
    fname = "test_data_grb/lcs_tophat.h5"
    if (afterglowpy) : Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
                            'specType':    0,                  # Basic Synchrotron Spectrum
                            'counterjet':  0,
                            'spread':      0,
                            'thetaObs':    0.0,   # Viewing angle in radians
                            'E0':          1.0e52, # Isotropic-equivalent energy in erg
                            'g0':          1000,
                            'thetaCore':   0.2,    # Half-opening angle in radians
                            'thetaWing':   0.2,
                            'n0':          1e-3,    # circumburst density in cm^{-3}
                            'p':           2.2,    # electron energy distribution index
                            'epsilon_e':   0.1,    # epsilon_e
                            'epsilon_B':   0.01,   # epsilon_B
                            'xi_N':        1.0,    # Fraction of electrons accelerated
                            'd_L':         3.09e26, # Luminosity distance in cm
                            'z':           0.0099}   # redshift

    t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9
    if (afterglowpy) : Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


    if (afterglowpy) :
        ax.plot(t, Fnu, color='gray', label='afterglowpy')

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    ax.plot(our_times, our_fluxes, color="black", ls="-", label="PyBlastAfterglow [PW]")
    #
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    #         break

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Top-hat [1GGz]")
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()

    # -------------------------------------------------

    plt.show()
# compareTopHatLightCurves()

def compareTopHatLightCurvesMethods():

    ''' Run 'compareTopHatLightCurves()' first '''
    fnames = {
        "test_data_grb/lcs_tophat_a_obs.h5":"blue",
        "test_data_grb/lcs_tophat_pw_obs.h5":"cyan",
        "test_data_grb/lcs_tophat_a_comov.h5":"green",
        "test_data_grb/lcs_tophat_pw_comov.h5":"lime"
    }
    if (afterglowpy) : Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
                            'specType':    0,                  # Basic Synchrotron Spectrum
                            'counterjet':  0,
                            'spread':      0,
                            'thetaObs':    0.0,   # Viewing angle in radians
                            'E0':          1.0e52, # Isotropic-equivalent energy in erg
                            'g0':          1000,
                            'thetaCore':   0.2,    # Half-opening angle in radians
                            'thetaWing':   0.2,
                            'n0':          1e-3,    # circumburst density in cm^{-3}
                            'p':           2.2,    # electron energy distribution index
                            'epsilon_e':   0.1,    # epsilon_e
                            'epsilon_B':   0.01,   # epsilon_B
                            'xi_N':        1.0,    # Fraction of electrons accelerated
                            'd_L':         3.09e26, # Luminosity distance in cm
                            'z':           0.0099}   # redshift

    t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9
    if (afterglowpy) : Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


    if (afterglowpy) :
        ax.plot(t, Fnu, color='gray', label='afterglowpy')

    for fname in fnames.keys():
        dfile = h5py.File(curdir+fname, "r")
        for key in dfile.keys():
            print(key)
        our_times = np.array( dfile["times"] )
        our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

        ax.plot(our_times, our_fluxes, color=fnames[fname], ls="-", label=fname.split("/")[-1])
    #
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    #         break

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Top-hat [1GGz]")
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()

    # -------------------------------------------------

    plt.show()
# compareTopHatLightCurvesMethods()

def compareTopHatLightCurvesMethods():

    ''' Run 'compareTopHatLightCurves()' first '''
    fnames = {
        "test_data_grb/lcs_tophat_a_obs.h5":"blue",
        "test_data_grb/lcs_tophat_pw_obs.h5":"cyan",
        "test_data_grb/lcs_tophat_a_comov.h5":"green",
        "test_data_grb/lcs_tophat_pw_comov.h5":"lime"
    }
    Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
         'specType':    0,                  # Basic Synchrotron Spectrum
         'counterjet':  0,
         'spread':      0,
         'thetaObs':    0.0,   # Viewing angle in radians
         'E0':          1.0e52, # Isotropic-equivalent energy in erg
         'g0':          1000,
         'thetaCore':   0.2,    # Half-opening angle in radians
         'thetaWing':   0.2,
         'n0':          1e-3,    # circumburst density in cm^{-3}
         'p':           2.2,    # electron energy distribution index
         'epsilon_e':   0.1,    # epsilon_e
         'epsilon_B':   0.01,   # epsilon_B
         'xi_N':        1.0,    # Fraction of electrons accelerated
         'd_L':         3.09e26, # Luminosity distance in cm
         'z':           0.0099}   # redshift

    t = np.geomspace(1.0 * 86400.0, 1.0e3 * 86400.0, 100)
    nu = np.empty(t.shape)
    nu[:] = 3.0e9
    Fnu = grb.fluxDensity(t, nu, **Z)

    fig, ax = plt.subplots(figsize=(6, 4.5), ncols=1, nrows=1)


    ax.plot(t, Fnu, color='gray', label='afterglowpy')

    for fname in fnames.keys():
        dfile = h5py.File(curdir+fname, "r")
        for key in dfile.keys():
            print(key)
        our_times = np.array( dfile["times"] )
        our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

        ax.plot(our_times, our_fluxes, color=fnames[fname], ls="-", label=fname.split("/")[-1])
    #
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    #         break

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Top-hat [1GGz]")
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()

    # -------------------------------------------------

    plt.show()
compareTopHatLightCurvesMethods()




def compareGaussianOffAxisGflatSStructLightCurves():
    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    fname = "test_data_grb/lcs_FernadOffAxisGflatSSstruct_methods.h5"

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)

    dfile = h5py.File(curdir+fname, "r")
    for key in dfile.keys():
        print(key)
    our_times = np.array( dfile["times"] )
    our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

    # colnames, data = load_data(curdir + fname)
    # nx, ny = data.shape
    # t_arr = data[:, 0]
    # icolor1, icolor2 = 0, 0
    # for icol in range(1,len(colnames)):
    #     name = colnames[icol]
    #     row = data[:,icol]
    #     if name.__contains__("[A]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
    #         icolor1 += 1
    #     elif name.__contains__("[P]"):
    #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
    #         icolor2 += 1
    ax.plot(our_times, our_fluxes, ls='-', color='black', label='PyBlastAfterglow [pw]')
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()
# compareGaussianOffAxisGflatSStructLightCurves()

def compareGaussianOffAxisGflatSStructLightCurvesMethods():

    fnames = {
        "test_data_grb/lcs_FernandOffAxisGauss_a_obs.h5":"blue",
        "test_data_grb/lcs_FernandOffAxisGauss_pw_obs.h5":"cyan",
        "test_data_grb/lcs_FernandOffAxisGauss_a_comov.h5":"green",
        "test_data_grb/lcs_FernandOffAxisGauss_pw_comov.h5":"lime"
    }

    fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

    # fname = "test_data_grb/lcs_FernadOffAxisGflatSSstruct_methods.h5"

    tts_afgpy, ffs_afgpy = np.loadtxt(curdir+"test_data_grb/afgpy_grb170817.txt",unpack=True)
    ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

    tts, fluxes = np.loadtxt(curdir+"test_data_grb/jelib_grb170817.txt",unpack=True)
    ax.plot(tts, fluxes*1e26, ls='--', color='gray', label='Joelib')

    for fname in fnames.keys():
        dfile = h5py.File(curdir+fname, "r")
        for key in dfile.keys():
            print(key)
        our_times = np.array( dfile["times"] )
        our_fluxes = np.array( dfile["totalflux at freq=3.0000e+09"] )

        # colnames, data = load_data(curdir + fname)
        # nx, ny = data.shape
        # t_arr = data[:, 0]
        # icolor1, icolor2 = 0, 0
        # for icol in range(1,len(colnames)):
        #     name = colnames[icol]
        #     row = data[:,icol]
        #     if name.__contains__("[A]"):
        #         ax.loglog(t_arr, row, color=get_color(icolor1), ls='-', label=name)
        #         icolor1 += 1
        #     elif name.__contains__("[P]"):
        #         ax.loglog(t_arr, row, color=get_color(icolor2), ls='--', label=name)
        #         icolor2 += 1

        ax.plot(our_times, our_fluxes, ls='-', color=fnames[fname], label=fname.split("/")[-1])
    # r"jet $E_0$={:.1e} $\Gamma_0$={:.0f} $\theta$={:.2f} $\epsilon_e$={:.2f} $\epsilon_b$={:.2f} $p$={:.2f} $\nu$={:.2e} Hz".format(
    # Z["E0"],Z["g0"],Z["thetaWing"],Z["epsilon_e"],Z["epsilon_B"],Z["p"],nu[0]))
    ax.grid()
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("time [s]", fontsize=12)
    ax.set_ylabel("Flux density [mJy]", fontsize=12)
    ax.set_title(r"Gaussian jet Off-axis Gflat [3GGz]")
    ax.grid()
    plt.show()
compareGaussianOffAxisGflatSStructLightCurvesMethods()
