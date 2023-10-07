import copy
import os.path
import numpy as np
import h5py
import shutil
import hashlib
from scipy import ndimage, interpolate
from glob import glob
from sklearn.metrics import mean_squared_error

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


from .interface import PyBlastAfterglow, Skymap, Ejecta
from .utils import cgs, get_beta, get_Beta, get_Gamma
from .id_analytic import JetStruct
from .parfile_tools import read_parfile, modify_parfile, modify_parfile_par_opt
from .skymap_process import ProcessRawSkymap
from .skymap_plotting_tools import plot_skymap_with_hists, plot_pcolomesh, full_plot_skymap_with_hists

try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

def plot_skymaps_3d(ax : plt.axes, skymaps : list[Skymap]):

    # g_min = np.min([np.min(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_gmax = np.max([np.max(skymap.im_hist[skymap.im_hist > 0 & np.isfinite(skymap.im_hist)]) for skymap in skymaps])
    g_min = g_gmax * 1e-4
    print(f"g_min={g_min} g_gmax={g_gmax}")

    # fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
    #                        gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)

    # ax = fig.add_subplot(211,projection='3d')
    #
    #
    # for skymap in skymaps:
    #     ax.plot(skymap.grid_x, np.log10(skymap.dist_x), zs=np.log10(skymap.time), zdir='x', color="white")
    # # ax.set_zscale('log')
    #
    #
    # ax.set_ylabel(r"$X$ [mas]")
    # ax.set_zlabel(r"$Z$ [mas]")
    # ax.set_xlabel(r"time [s]")
    # ax.set_facecolor("black")
    # ax.grid(False)
    # ax.w_xaxis.pane.fill = False
    # ax.w_yaxis.pane.fill = False
    # ax.w_zaxis.pane.fill = False
    # ax.set_box_aspect(aspect = (4,2,2))
    # ax.view_init(elev=5, azim=150, roll=0)

    # ---------------

    # Create the 4th color-rendered dimension
    cmap = cm.get_cmap('Greys')
    cmap.set_under("white")
    # cmap.set_over("red")
    scam = plt.cm.ScalarMappable(
        norm=cm.colors.LogNorm(g_min, g_gmax),
        cmap=cmap # see https://matplotlib.org/examples/color/colormaps_reference.html
    )


    # ax = fig.add_subplot(projection='3d')
    zmin = 0.
    xmin, xmax, ymin,ymax = 0,0,0,0
    for skymap in skymaps:
        X, Y = np.meshgrid(skymap.grid_x, skymap.grid_y)
        if (zmin > np.min(Y)): zmin = np.min(Y)
        if (xmin > np.min(X)): xmin = np.min(X)
        if (ymin > np.min(Y)): ymin = np.min(Y)
        if (xmax < np.max(X)): xmax = np.max(X)
        if (ymax < np.max(Y)): ymax = np.max(Y)

        Z = np.full_like(X,fill_value=np.log10(skymap.time / cgs.day))
        G = skymap.im_hist.T
        # X = X[G > g_min]
        # Y = Y[G > g_min]
        Z[G < g_min] = np.nan
        G[G < g_min] = g_min
        scam.set_array([])
        facecolors = scam.to_rgba(G)
        # ax.scatter(Z, X, Y, color=facecolors)
        ax.plot_surface(
            Z, X, Y,
            facecolors  = facecolors,
            antialiased = True,
            rstride=1, cstride=1, alpha=.6, shade=False,

        )
    for skymap in skymaps:
        ax.plot(skymap.xc, zmin, zs=np.log10(skymap.time/cgs.day), zdir='x', marker="o", ms=1, color="black")
        ax.plot([skymap.x1,skymap.x2], [zmin,zmin], zs=np.log10(skymap.time/cgs.day),
                zdir='x', ls="--", lw=.6, color="black")

        zmin_ = -.5
        ax.plot(np.log10(skymap.time/cgs.day), zmin_, zs=skymap.yc, zdir='z', marker="o", ms=1, color="black")
        ax.plot([np.log10(skymap.time/cgs.day),np.log10(skymap.time/cgs.day)],
                [zmin_,zmin_], zs=[skymap.y1,skymap.y2],
                zdir='z', ls="--", lw=.6, color="black")

        # ax.plot(skymap.grid_x, skymap.dist_x, zs=np.log10(skymap.time), zdir='x', label='curve in (x, y)')
        # ax.plot(skymap.grid_y, skymap.dist_y, zs=np.log10(skymap.time), zdir='z', label='curve in (x, y)')
    times = [np.log10(skymap.time/cgs.day) for skymap in skymaps]
    # xmin = times.min(),
    # xmax = times.max()
    # n = len(times)
    # for t
    ax.set_ylim(0,11)
    ax.set_ylabel(r"$X$ [mas]",fontsize=12,labelpad=10)
    ax.set_zlabel(r"$Z$ [mas]",fontsize=12,labelpad=5)
    ax.set_xlabel(r"$\log(t_{\rm obs})$ [day]",fontsize=12,labelpad=18)
    ax.minorticks_on()
    # ax.set_facecolor("black")
    ax.grid(False)
    ax.w_xaxis.pane.fill = False
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False
    # ax.w_zaxis.line.set_visible(False)
    ax.set_box_aspect(aspect = (6,2,2))
    ax.view_init(elev=10, azim=150, roll=0)
    ax.tick_params(direction='in', length=10, width=2, colors='black',
                   grid_color='gray', grid_alpha=0.1,which="both", axis="both",labelsize=12)
    # x_scale=4
    # y_scale=1
    # z_scale=1
    #
    # scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    # scale=scale*(1.0/scale.max())
    # scale[3,3]=1.0
    # def short_proj():
    #     return np.dot(Axes3D.get_proj(ax), scale)
    # fig.set_facecolor('black')
    # fig.subplots_adjust(left=-.2,top=.2)  # plot outside the normal area
    # fig.tight_layout()

    # ax.get_proj=short_proj
    # ax.set_title(r"Gaussian jet with $\texttt{PW}$ method",loc='center',color="black")
    # plt.title(r"Gaussian jet with $\texttt{PW}$ method",loc='center',color="black")
    # ax.text(x=.8, y=.1, z=.6, s=r'Gaussian jet with $\texttt{PW}$ method', horizontalalignment='center',
    #      verticalalignment='center', transform=ax.transAxes)
    # ax.annotate(r'Gaussian jet with $\texttt{PW}$ method', xy=(2, 1), xytext=(-200, 200), textcoords='offset points', ha='left', bbox=dict(boxstyle='circle', fc='green', alpha=0.7),
    #          arrowprops=dict(arrowstyle='->'))

class RefData():
    ''' Load .h5 files with data dynamics data from afterglowpy  '''
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
    ''' load light curve .h5 files from afterglowpy '''
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


class Base():
    ''' change parfile and run PyBlastAfterglow for piece-wise [pw] or adaptive [a] eats itegration  '''

    conf = {"nx":128, "ny":64, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
            "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
            "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }}

    def __init__(self, default_parfile_fpath, workingdir):
        if not os.path.isdir(workingdir):
            raise IOError(f"Workingdir dir does not exists: {workingdir}")
        self.workingdir = workingdir
        if not os.path.isfile(default_parfile_fpath):
            raise IOError(f"default_parfile_fpath not exists: {default_parfile_fpath}")
        self.default_parfile_fpath = default_parfile_fpath

    def run_pw(self, pars : dict, struct : dict, opts : dict, opts_grb : dict) -> PyBlastAfterglow:
        # workdir = os.getcwd()+'/'
        # copy the main parfile into workdir
        parfilefname = self.default_parfile_fpath.split("/")[-1]
        shutil.copy(self.default_parfile_fpath, self.workingdir+parfilefname)
        # prepare initial data
        pba_id = JetStruct(n_layers_pw=struct["nlayers_pw"], n_layers_a=struct["nlayers_a"])
        id_dict, id_pars = pba_id.get_1D_id(pars=struct,type="piece-wise")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=self.workingdir+"gauss_grb_id.h5")

        # modify parfile
        modify_parfile_par_opt(workingdir=self.workingdir, part="main",
                                                 newpars=pars,
                                                 newopts=opts,
                                                 parfile=parfilefname, newparfile="parfile.par", keep_old=True)
        modify_parfile_par_opt(workingdir=self.workingdir, part="grb",
                                                 newpars=pars,
                                                 newopts=opts_grb,
                                                 parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        if not os.path.isfile(self.workingdir+"parfile.par"):
            raise FileNotFoundError(f"Parfile not found in {self.workingdir}"+"parfile.par")
        # init model interface
        pba_pw = PyBlastAfterglow(workingdir=self.workingdir,readparfileforpaths=True, parfile="parfile.par")
        # run
        pba_pw.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
                   loglevel="info")
        # postprocess skymaps
        if (pba_pw.GRB.opts["do_skymap"]=="yes"):
            prep = ProcessRawSkymap(conf=self.conf, verbose=False)
            prep.process_singles(infpaths=self.workingdir+"raw_skymap_*.h5",
                                 outfpath=pba_pw.GRB.fpath_sky_map,
                                 remove_input=False)

        # remove the default parfile from the working dir
        os.remove(self.workingdir+parfilefname)
        return pba_pw

    def run_a(self, pars : dict, struct : dict, opts : dict, opts_grb : dict) -> PyBlastAfterglow:
        # workdir = os.getcwd()+'/'
        parfilefname = self.default_parfile_fpath.split("/")[-1]
        # copy the main parfile into workdir
        shutil.copy(self.default_parfile_fpath, self.workingdir+parfilefname)
        # prepare initial data
        pba_id = JetStruct(n_layers_pw=struct["nlayers_pw"], n_layers_a=struct["nlayers_a"])
        id_dict, id_pars = pba_id.get_1D_id(pars=struct,type="adaptive")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=self.workingdir+"gauss_grb_id.h5")

        # modify parfile
        modify_parfile_par_opt(workingdir=self.workingdir, part="main",
                                                 newpars=pars,
                                                 newopts=opts,
                                                 parfile=parfilefname, newparfile="parfile.par", keep_old=True)
        modify_parfile_par_opt(workingdir=self.workingdir, part="grb",
                                                 newpars=pars,
                                                 newopts=opts_grb,
                                                 parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        if not os.path.isfile(self.workingdir+"parfile.par"):
            raise FileNotFoundError(f"Parfile not found in {self.workingdir}"+"parfile.par")
        pba_a = PyBlastAfterglow(workingdir=self.workingdir,readparfileforpaths=True, parfile="parfile.par")
        pba_a.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
                  loglevel="info")
        # postprocess skymaps
        if (pba_a.GRB.opts["do_skymap"]=="yes"):
            prep = ProcessRawSkymap(conf=self.conf, verbose=True)
            prep.process_singles(infpaths=self.workingdir+"raw_skymap_*.h5",
                                 outfpath=pba_a.GRB.fpath_sky_map,
                                 remove_input=False)
        # remove the default parfile from the working dir
        os.remove(self.workingdir+parfilefname)
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


class CasesFS(Base):

    def __init__(self, default_parfile_fpath, workingdir):
        super().__init__(default_parfile_fpath, workingdir)
        # if not (os.path.isdir(figsdir)):
        #     raise NotADirectoryError(f"Dire for figs is not found {figsdir}")
        # self.figsdir = figsdir

    def plot_lcs_ref(self, ax, ref : RefDataLC, nlayers : int, layers = (),):
        cmap = cm.get_cmap('Greys_r')
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

    def plot_generic(self, pars : dict, opts_a : dict, struct : dict, title : str,
                     ref_dyn_fname : str, ref_lc_fname : str,
                     figpath : str or None, save_pdf = True, show_fig = False):
        # run the code for given pars
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)

        ref = RefData(workdir=self.workingdir, fname=ref_dyn_fname)#"reference_afgpy_dyn.h5")
        ref_lc = RefDataLC(workdir=self.workingdir, fname=ref_lc_fname)#"reference_lc_0deg_layer.h5")

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
        # Z = {'jetType':     grb.jet.TopHat if struct["struct"] == "tophat" else grb.jet.Gaussian,     # Top-Hat jet
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

        # plot
        ax = axes[0]
        # ax.plot(t, Fnu, ls='-', color='gray', label='afterglopy')

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

        if not figpath is None:
            print("Saving:\n {}".format(figpath))
            if save_pdf: plt.savefig(figpath+".pdf")
            plt.savefig(figpath+".png", dpi=256)
        if show_fig: plt.show()

    def paper_plot_compare_spreading(self, pars : dict, opts_a:dict, struct : dict, title : str,
                                     ref_dyn_fname : str or None, ref_lc_fname : str or None,
                                     figpath : str or None, save_pdf = True, show_fig = True, layers=()):
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_0deg_layer.h5")

        if ref_lc_fname is None:
            if (pars["theta_obs"]==0): _fname="reference_lc_0deg_layer.h5"
            elif(pars["theta_obs"]==0.785): _fname="reference_lc_0785deg_layer.h5"
            elif(pars["theta_obs"]==1.57): _fname="reference_lc_157deg_layer.h5"
            else:raise KeyError("no ref.data for theta_obs={}".format(pars["theta_obs"]))
            ref_lc = RefDataLC(workdir=self.workingdir, fname=_fname)
        else:
            ref_lc = RefDataLC(workdir=self.workingdir, fname=ref_lc_fname)

        if ref_dyn_fname is None:
            ref = RefData(workdir=self.workingdir, fname="reference_afgpy_dyn.h5")
        else:
            ref = RefData(workdir=self.workingdir, fname=ref_dyn_fname)


        fig, axes = plt.subplots(figsize=(6,6.5), ncols=1, nrows=3)

        pars["Gamma0_frac_when_start_spread"] = .1
        opts_a["method_limit_spread"] = "Gamma0Frac"

        # layers=()


        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"red", "label":r"R19"},
                      plot_layer={"ls":'-.', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        opts_a["method_spread"] = "AA"
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"blue", "label":"GP12"},
                      plot_layer={"ls":'-.', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        opts_a["method_spread"] = "Adi"
        pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"green", "label":"HDL99"},
                      plot_layer={"ls":'-.', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'--', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'--', "cmap":"Greens", "alpha":.9, "vmin":-50, "vmax":60})

        # # ---------------------------------------------
        # opts_a["method_spread"] = "None"
        # pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
        # # ---------------------------------------------
        #
        # self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
        #               plot={"ls":'-', "color":"grey", "label":"HDL99"},
        #               plot_layer={"ls":'-.', "cmap":"Greys", "alpha":.9, "vmin":-50, "vmax":60})
        # self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
        #               plot_layer={"ls":'--', "cmap":"Greys", "alpha":.9, "vmin":-50, "vmax":60})
        # self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
        #               plot_layer={"ls":'--', "cmap":"Greys", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------
        self.plot_lcs_ref(ax=axes[0], ref=ref_lc, nlayers=pba_a.GRB.get_lc_obj().attrs["nlayers"], layers=layers)
        self.plot_dyn_ref(axes[1], ref=ref, v_n_x="tburst", v_n_y="mom", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"], layers=layers)
        self.plot_dyn_ref(axes[2], ref=ref, v_n_x="tburst", v_n_y="theta", nlayers=pba_a.GRB.get_dyn_obj().attrs["nlayers"], layers=layers)

        # -------- Afterglopy --------------
        # Z = {'jetType':     grb.jet.TopHat if struct["struct"] == "tophat" else grb.jet.Gaussian,     # Top-Hat jet
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
        # t = np.geomspace(1.0 * 864.0, 1.0e3 * 86400.0, 100)
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
        ax.set_xlim(1e4,1e8)
        ax.set_ylim(1e-1,2e2)
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


        if not figpath is None:
            print("Saving:\n {}".format(figpath))
            if save_pdf: plt.savefig(figpath+".pdf")
            plt.savefig(figpath+".png", dpi=256)
        if show_fig: plt.show()

    def plot_170817_like(self, pars : dict, opts_a:dict, opts_pw:dict, struct : dict, title : str,
                         figpath : str, save_pdf = True, show_fig = False):
        # run the code for given pars


        fig, ax = plt.subplots(figsize=(9,2.5), ncols=1, nrows=1)

        # --- piece-wise
        pars_pw = copy.deepcopy(pars)
        pars_pw["mom0_frac_when_start_spread"] = 0.1
        opts_pw["fname_light_curve"]=f"lc_pw_res_170817.h5"
        opts_pw["fname_light_curve_layers"]=f"lc_dense_pw_res_170817.h5"
        pba_pw = self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=opts_pw)
        self.plot_lcs(ax=ax, pars=pars_pw, pba=pba_pw, layers = (),
                      plot={"ls":'-', "color":"green", "label":"PBA [PW]"},
                      plot_layer={"ls":'-', "cmap":"inferno", "alpha":.5})

        # --- adaptive
        pars_a = copy.deepcopy(pars)
        pars_a["mom0_frac_when_start_spread"] = 0.95
        opts_a["fname_light_curve"]=f"lc_a_res_170817.h5"
        opts_a["fname_light_curve_layers"]=f"lc_dense_a_res_170817.h5"
        pba_a = self.run_a(struct=struct, pars=pars_a, opts={}, opts_grb=opts_a)
        self.plot_lcs(ax=ax, pars=pars_a, pba=pba_a, layers = (),
                      plot={"ls":'-', "color":"red", "label":"PBA [A]"},
                      plot_layer={"ls":'-', "cmap":"viridis", "alpha":.5})

        # --------- Reference Models -----------
        tts_afgpy, ffs_afgpy = np.loadtxt(self.workingdir+"./afgpy_grb170817.txt",unpack=True)
        ax.plot(tts_afgpy, ffs_afgpy, ls=':', color='gray', label='afterglopy')

        tts, fluxes = np.loadtxt(self.workingdir+"./jelib_grb170817.txt",unpack=True)
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

        if not figpath is None:
            print("Saving:\n {}".format(figpath))
            if save_pdf: plt.savefig(figpath+".pdf")
            plt.savefig(figpath+".png", dpi=256)
        if show_fig: plt.show()

    def plot_3d_stacked_skymaps(self, pars : dict, opts_a:dict, opts_pw:dict, struct : dict, title : str,
                                figpath : str, save_pdf = True, show_fig = False):
        # --- piece-wise
        freq = pars["obs_freq"]
        pars_pw = copy.deepcopy(pars)
        pars_pw["mom0_frac_when_start_spread"] = 0.1
        opts_pw["do_skymap"] = "yes"
        opts_pw["fname_sky_map"]=f"skymap_pw_for_3d_plot.h5"
        pba_pw = self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=opts_pw)

        skymaps = [pba_pw.GRB.get_skymap(time, freq=freq) for time in pba_pw.GRB.get_skymap_times()[:-4]]
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
                               gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)
        plot_skymaps_3d(ax=ax, skymaps=skymaps)
        fig.suptitle(r"Gaussian jet with $\texttt{PW}$ method", y=0.90, fontsize=16)
        if not figpath is None:
            _figpath = figpath+"_PW"
            print("Saving:\n {}".format(_figpath))
            if save_pdf: plt.savefig(_figpath+".pdf")
            plt.savefig(_figpath+".png", dpi=256)
        if show_fig: plt.show()

        # --- adaptive
        freq = pars["obs_freq"]
        pars_a = copy.deepcopy(pars)
        pars_a["mom0_frac_when_start_spread"] = 0.95
        opts_a["do_skymap"] = "yes"
        opts_a["fname_sky_map"]=f"skymap_a_for_3d_plot.h5"
        pba_a = self.run_a(struct=struct, pars=pars_a, opts={}, opts_grb=opts_a)

        skymaps = [pba_a.GRB.get_skymap(time, freq=freq) for time in pba_a.GRB.get_skymap_times()[:-4]]
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'},figsize=(5.2,4),
                               gridspec_kw=dict(top=1.2, bottom=-.2, left=-.2, right=1)) # figsize=(4.2,4)figsize=(5,4)
        plot_skymaps_3d(ax=ax,skymaps=skymaps)
        fig.suptitle(r"Gaussian jet with $\texttt{A}$ method", y=0.90, fontsize=16)
        if not figpath is None:
            _figpath = figpath+"_A"
            print("Saving:\n {}".format(_figpath))
            if save_pdf: plt.savefig(_figpath+".pdf")
            plt.savefig(_figpath+".png", dpi=256)
        if show_fig: plt.show()


    # --- SKYMAPS ----

    # def compare_skymaps_theta_im_max(self, nlayer_a = 321, theta_maxs=((0.2, 0.4, .9, 1.2,  1.57),
    #                                                                    ('red','orange','yellow', 'cyan', 'lime'))):
    #
    #     fig,axes = plt.subplots(ncols=1,nrows=len(theta_maxs[0])+1,sharex='all',figsize=(5,10))
    #     times = [1.,10.,40.,100.,200.]
    #     freq = 1e9
    #     tmp={"hist_nx": 71, "hist_ny": 71, "spec": False,
    #          "smooth": {},  # {"type": "gaussian", "sigma": 10},
    #          "cm": {"color": 'cyan', "marker": "+"},
    #          "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
    #          "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
    #          "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
    #                         "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
    #          }
    #     settings={
    #         "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
    #         "title": {"title": "time_fluxratio"},  # "time_fluxratio"
    #         "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
    #         "xlabel": "x [mas]", "ylabel": "z [mas]",
    #         "histx_backgound_color": "black",
    #         "histy_backgound_color": "black",
    #         "plot_grids": True,
    #         "histx_lim":(1e-2, 1e1),
    #         "histy_lim":(1e-2, 1e1)
    #     }
    #
    #     theta_w = 0.2# np.pi/2.
    #     tot_fluxes = {}
    #     for i in range(len(theta_maxs[0])):
    #         theta_max = theta_maxs[0][i]
    #         # nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0
    #
    #         nlayer_pw = 80
    #
    #         color=theta_maxs[1][i]
    #         tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color
    #
    #         prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
    #                               "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
    #         modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
    #                                parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    #         modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
    #                                newpars={"nsublayers":nlayer_a, "im_max_theta":theta_max},
    #                                newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive"},
    #                                parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    #         pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
    #         pba_a.run(loglevel='info')
    #
    #         all_x_jet, all_y_jet, all_fluxes_jet \
    #             = pba_a.GRB.get_skymap(time=times[3] * cgs.day, freq=freq, verbose=False, remove_mu=True, renormalize=True, normtype='pw')
    #         int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
    #                                                     hist_or_int="hist", shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
    #         grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
    #                                                                collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
    #         grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
    #                                                                collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
    #         xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
    #         ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
    #         im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
    #                                      grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)
    #
    #         tot_fluxes[f"nlayer_a={nlayer_a} theta_max={theta_max}"] = f" A={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=times[3]* cgs.day)}"
    #         # --- COLOBAR
    #         # ax_cbar = axes[1,i+1]
    #         # divider = make_axes_locatable(ax_cbar)
    #         # cax = divider.append_axes('right', size='99%', pad=0.9)
    #         # plt.delaxes(ax_cbar)
    #         # cbar = plt.colorbar(im, cax=cax,
    #         #                     # format='%.0e',ticks=ticks
    #         #                     orientation='vertical',
    #         #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
    #         #                     )
    #     for key, val in tot_fluxes.items():
    #         print(f"{key} {val}")
    #     # pass
    #
    #     axes[0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
    #     axes[0].set_yscale("linear")
    #     axes[0].set_xscale("linear")
    #     axes[0].minorticks_on()
    #     axes[0].set_yscale("log")
    #     axes[0].tick_params(axis='both', which='both', labelleft=True,
    #                         labelright=False, tick1On=True, tick2On=True,
    #                         labelsize=12, color="white",
    #                         direction='in',
    #                         bottom=True, top=True, left=True, right=True)
    #     axes[0].set_facecolor(settings["histx_backgound_color"])
    #     for j_, ax_main in enumerate(axes[1:]):
    #         ax_main.set_yscale("linear")
    #         ax_main.set_xscale("linear")
    #         ax_main.minorticks_on()
    #         ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
    #         ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
    #         ax_main.tick_params(axis='both', which='both', labelleft=True,
    #                             labelright=False, tick1On=True, tick2On=True,
    #                             labelsize=12, color="white",
    #                             direction='in',
    #                             bottom=True, top=True, left=True, right=True)
    #         ax_main.set_ylabel("$Z$ [mas]")
    #
    #
    #         if settings["plot_grids"]:
    #             ax_main.grid()
    #
    #     plt.tight_layout()
    #     plt.show()

 #    def compare_skymaps_3d(self, resolutions=((20,40,80,100,120),
 #                                              # (21,41,81,101,121),
 #                                              (21,81,121,161,201),
 #                                              # (121,241,381,401,521),
 #                                              ('red','orange','yellow', 'cyan', 'lime')), log=True):
 #
 #        theta_w = 0.2# np.pi/2.
 #        times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
 #        time = 1600
 #        freq = 1e9
 #
 #        nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
 #        nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0
 #
 #        fig = plt.figure(figsize=plt.figaspect(0.5))
 #
 #        nres = np.max([nres_pw,nres_a])
 #        for i in range(nres):
 #            nlayer_a = resolutions[1][i] if nres_a > 0 else 0
 #            nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0
 #
 #
 #            prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
 #                                  "nlayers_pw": nlayer_pw, "nlayers_a": 1, "struct":"tophat"},type='pw',outfpath="tophat_grb_id_pw.h5")
 #
 #            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
 #                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
 #            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
 #                                   newpars={},
 #                                   newopts={"fname_ejecta_id":"tophat_grb_id_pw.h5","method_eats":"piece-wise",
 #                                            "fname_sky_map":"skymap_tophat_pw.h5", "fname_light_curve":"lc_tophat_pw.h5"},
 #                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
 #            pba_pw = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
 #            pba_pw.run(loglevel='info')
 #
 #            all_x_jet_pw, all_y_jet_pw, all_fluxes_jet_pw \
 #                = pba_pw.GRB.get_skymap_old(time=time * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True)
 #
 #            ax_a = fig.add_subplot(1, 2, 1, projection='3d')
 #
 #            for (x_jet_pw, y_jet_pw, fluxes_jet_pw) in zip(all_x_jet_pw, all_y_jet_pw, all_fluxes_jet_pw):
 #                # fluxes_jet_pw+=1
 #                # fluxes_jet_pw[fluxes_jet_pw<1e-4] = 0.
 #                cmap = cm.get_cmap('inferno')
 #                if log:
 #                    my_norm = LogNorm(fluxes_jet_pw.max() * 1e-2, fluxes_jet_pw.max())
 #                    c = cmap(my_norm(fluxes_jet_pw.flatten()))
 #                    fluxes_jet_pw = np.log10(fluxes_jet_pw)
 #                else:
 #                    my_norm = Normalize(fluxes_jet_pw.max() * 1e-2, fluxes_jet_pw.max())
 #                    c = cmap(my_norm(fluxes_jet_pw.flatten()))
 #                    fluxes_jet_pw = fluxes_jet_pw
 #                ax_a.scatter(x_jet_pw.flatten(), y_jet_pw.flatten(), fluxes_jet_pw.flatten(),  c=c)
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cphis_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #                ax_a.set_xlabel('X Label')
 #                ax_a.set_ylabel('Y Label')
 #                ax_a.set_zlabel('I Label')
 #                # ax_a.set_zlim(np.log10(fluxes_jet_pw.flatten()).max()*1e-2,np.log10(fluxes_jet_pw.flatten()).max())
 #            ax_a.set_title(f"Tot.Flux={pba_pw.GRB.get_skymap_totfluxes(freq=freq,time=time * cgs.day)}")
 #            # ---------------------- #
 #
 #            prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": theta_w, "theta_w": theta_w,
 #                                  "nlayers_pw": 1, "nlayers_a": 1, "struct":"tophat"},type='a',outfpath="tophat_grb_id_a.h5")
 #            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
 #                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
 #            modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
 #                                   newpars={"nsublayers":nlayer_a},
 #                                   newopts={"fname_ejecta_id":"tophat_grb_id_a.h5","method_eats":"adaptive",
 #                                            "fname_sky_map":"skymap_tophat_a.h5", "fname_light_curve":"lc_tophat_a.h5"},
 #                                   parfile="parfile.par", newparfile="parfile.par", keep_old=False)
 #            pba_a = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")
 #            pba_a.run(loglevel='info')
 #
 #            all_x_jet_a, all_y_jet_a, all_fluxes_jet_a \
 #                = pba_a.GRB.get_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True, renormalize=False, normtype='pw', remove_zeros=True)
 #
 #
 #
 #            ax_b = fig.add_subplot(1, 2, 2, projection='3d')
 #            xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a)
 #
 #            my_norm_a = None
 #            for (x_jet_a, y_jet_a, fluxes_jet_a) in zip(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a):
 #                # fluxes_jet_a += 1
 #                # fluxes_jet_a[fluxes_jet_a<1e-4] = 0.
 #                cmap = cm.get_cmap('inferno')
 #                if log:
 #                    if (my_norm_a is None): my_norm_a = LogNorm(fluxes_jet_a.max() * 1e-2, fluxes_jet_a.max())
 #                    c = cmap(my_norm_a(fluxes_jet_a.flatten()))
 #                    fluxes_jet_a = np.log10(fluxes_jet_a)
 #                else:
 #                    if (my_norm_a is None): my_norm_a = Normalize(fluxes_jet_a.max() * 1e-2, fluxes_jet_a.max())
 #                    c = cmap(my_norm_a(fluxes_jet_a.flatten()))
 #                    fluxes_jet_a = fluxes_jet_a
 #                # ax = fig.add_subplot(projection='3d')
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(rs_i).flatten(), c=cmap(my_norm(int_i.flatten())))
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cthetas_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #
 #                ax_b.scatter(x_jet_a.flatten(), y_jet_a.flatten(), fluxes_jet_a.flatten(),  c=c)
 #                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cphis_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #                ax_b.set_xlabel('X Label')
 #                ax_b.set_ylabel('Y Label')
 #                ax_b.set_zlabel('I Label')
 #                ax_b.scatter(np.full_like(fluxes_jet_a,xc_m_j), np.full_like(fluxes_jet_a,yc_m_j), fluxes_jet_a.flatten())
 #
 #            # ax_b.set_zlim(np.log10(fluxes_jet_a.flatten()).max()*1e-2,np.log10(fluxes_jet_a.flatten()).max())
 #            ax_b.set_title(f"Tot.Flux={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=time * cgs.day)}")
 #
 #
 #        print(f"Piecewise: I[{fluxes_jet_pw.min(),fluxes_jet_pw.max()}]")
 #        print(f"Adaptive: I[{fluxes_jet_a.min(),fluxes_jet_a.max()}]")
 #
 #        # plt.tight_layout()
 #        plt.show()
 #
 #    def compare_skymaps_3d_theta_im_max(self, pars : dict, opts_a : dict, struct : dict, title : str,
 #                                        theta_maxs=(
 #                                                (1.57,), #0.4, .9, 1.2,  1.57),
 #                                                # (21,121,121,121,121),
 #                                                ((1,22),),#121,121,121,121),
 #                                                ('red','orange','yellow', 'cyan', 'lime')),
 #                                        extend=2, nx = 81, ny = 81):
 #
 #        opts_a["do_skymap"] = "yes"
 #
 #        theta_w = 0.2# np.pi/2.
 #        times = np.array([1.,10.,40.,100.,200.,400.,800.,1600.])
 #        time_ = 7000#1600
 #        freq = pars["obs_freq"]
 #        # nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
 #        # nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0
 #
 #        fig = plt.figure(figsize=(15,10))#plt.figaspect(0.5),fig)
 #        res = {}
 #        # nres = np.max([nres_pw,nres_a])
 #
 #        my_norm = None
 #        my_norm2 = None
 #        for i in range(len(theta_maxs[0])):
 #            struct["nlayers_a"] = theta_maxs[1][i][0]
 #            pars["nsublayers"] = theta_maxs[1][i][1]
 #            pars["im_max_theta"] = theta_maxs[0][i]
 #
 #            pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
 #            x_jet_a, y_jet_a, fluxes_jet_a, ctheta_jet, cphi_jet, r_jet \
 #                = pba_a.GRB.get_skymap(time=time_ * cgs.day, freq=freq, verbose=True, remove_mu=True,
 #                                       renormalize=False, normtype='a', remove_zeros=True, return_sph_coords=True)
 #
 #            # pba_a = self.run_pw(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
 #            # x_jet_a, y_jet_a, fluxes_jet_a \
 #            #     = pba_a.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=True, remove_mu=False,
 #            #                                renormalize=True, normtype='pw')
 #
 #
 #            # _x_jet_a = np.empty(0,)
 #            # _y_jet_a = np.empty(0,)
 #            # _fluxes_jet_a = np.empty(0,)
 #            # for _x, _y, _val, _cth, _cph, _r in zip(x_jet_a[0], y_jet_a[0], fluxes_jet_a[0], ctheta_jet[0], cphi_jet[0], r_jet[0]):
 #            #     print(f"\tx={_x} y={_y} z={_val}  ctheta={_cth} _phi={_cph} _r={_r}")
 #            #     _x_jet_a = np.append(_x_jet_a,  _x)
 #            #     _y_jet_a = np.append(_x_jet_a,  _y)
 #            #     _fluxes_jet_a = np.append(_x_jet_a,  _val)
 #            #
 #            #
 #            #
 #            # _x_jet_a =x_jet_a[0]#[y_jet_a[0] != 0]
 #            # _fluxes_jet_a =fluxes_jet_a[0]#[y_jet_a[0] != 0]
 #            # _y_jet_a =y_jet_a[0]#[y_jet_a[0] != 0]
 #            #
 #            # print( f">0 = {len( _y_jet_a[_y_jet_a>0] )} <0 {len( _y_jet_a[_y_jet_a<0] )}" )
 #            # print( f">0 = {len( np.unique(_y_jet_a[_y_jet_a>0]) )} <0 {len( np.unique(_y_jet_a[_y_jet_a<0]) )}" )
 #            #
 #            # N = 10
 #            # xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm([_x_jet_a], [_y_jet_a], [_fluxes_jet_a])
 #            xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(x_jet_a, y_jet_a, fluxes_jet_a)
 #            print(f"xc={xc_m_j} yc={yc_m_j}")
 #
 #
 #            # # for (x_jet_a, y_jet_a, fluxes_jet_a) in zip(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a):
 #            # # fluxes_jet_a += 1
 #            # # fluxes_jet_a[fluxes_jet_a<1e-4] = 0.
 #            # cmap = cm.get_cmap('viridis')
 #            # # my_norm = Normalize(fluxes_jet_a.max() * 1e-1, fluxes_jet_a.max())
 #            # if my_norm is None: my_norm = Normalize(fluxes_jet_a[0].max() * 1e-1, fluxes_jet_a[0].max())
 #            # ##  ax_c = fig.add_subplot(len(theta_maxs[0]), 3, i*3+1, projection='3d')
 #            # ax_x = fig.add_subplot(len(theta_maxs[0]), 3, 1)
 #            # ax_y = fig.add_subplot(len(theta_maxs[0]), 3, 2)
 #            # _nshells = int(len(x_jet_a)/2)
 #            # for ish in range(_nshells):
 #            #     c = ax_x.scatter(x_jet_a[ish].flatten(), fluxes_jet_a[ish].flatten(), c=cmap(my_norm(fluxes_jet_a[ish].flatten())), alpha=0.3)
 #            #     c = ax_x.scatter(x_jet_a[_nshells+ish].flatten(), fluxes_jet_a[_nshells+ish].flatten(), c=cmap(my_norm(fluxes_jet_a[_nshells+ish].flatten())), alpha=0.3)
 #            #
 #            #     c = ax_y.scatter(y_jet_a[ish].flatten(), fluxes_jet_a[ish].flatten(), c=cmap(my_norm(fluxes_jet_a[ish].flatten())), alpha=0.3)
 #            #     c = ax_y.scatter(y_jet_a[_nshells+ish].flatten(), fluxes_jet_a[_nshells+ish].flatten(), c=cmap(my_norm(fluxes_jet_a[_nshells+ish].flatten())), alpha=0.3)
 #            # ax_x.axvline(x=xc_m_j)
 #            # ax_x.grid()
 #            # ax_y.axvline(x=yc_m_j)
 #            # ax_y.grid()
 #
 #
 #            # for (x_jet_a, y_jet_a, fluxes_jet_a) in zip(all_x_jet_a, all_y_jet_a, all_fluxes_jet_a):
 #            # fluxes_jet_a += 1
 #            # fluxes_jet_a[fluxes_jet_a<1e-4] = 0.
 #
 #            cmap_pj = cm.get_cmap('viridis')
 #            cmap_cj = cm.get_cmap('inferno')
 #            # my_norm = Normalize(fluxes_jet_a.max() * 1e-1, fluxes_jet_a.max())
 #            if my_norm is None: my_norm = Normalize(fluxes_jet_a[0].max() * 1e-1, fluxes_jet_a[0].max())
 #            # ax_c = fig.add_subplot(len(theta_maxs[0]), 3, i*3+1, projection='3d')
 #            ax_c = fig.add_subplot(len(theta_maxs[0]), 3, 1, projection='3d')
 #            _nshells = int(len(x_jet_a)/2)
 #            for ish in range(_nshells):
 #                if my_norm is None: my_norm = Normalize(fluxes_jet_a[ish].max() * 1e-1, fluxes_jet_a[ish].max())
 #                c = ax_c.scatter(x_jet_a[ish].flatten(), y_jet_a[ish].flatten(), fluxes_jet_a[ish].flatten(),
 #                                 c=cmap_pj(my_norm(fluxes_jet_a[ish].flatten())), alpha=.5)
 #                c = ax_c.scatter(x_jet_a[_nshells+ish].flatten(), y_jet_a[_nshells+ish].flatten(), fluxes_jet_a[_nshells+ish].flatten(),
 #                                 c=cmap_cj(my_norm(fluxes_jet_a[_nshells+ish].flatten())), alpha=.5)
 #
 #
 #
 #
 #            # #
 #            # ax_b = fig.add_subplot(len(theta_maxs[0]), 2, i*2+1, projection='3d')
 #            cmap_pj = cm.get_cmap('viridis')
 #            cmap_cj = cm.get_cmap('inferno')
 #            # if my_norm is None: my_norm = Normalize(fluxes_jet_a[0].max() * 1e-1, fluxes_jet_a[0].max())
 #            ax_b = fig.add_subplot(len(theta_maxs[0]), 3, i*2+2)
 #            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
 #            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(rs_i).flatten(), c=cmap(my_norm(int_i.flatten())))
 #            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #            # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),cthetas_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
 #            # fluxes_jet_=fluxes_jet_a;#np.log10(fluxes_jet_a)
 #            _nshells = int(len(x_jet_a)/2)
 #            for ish in range(_nshells):
 #                if my_norm is None: my_norm = Normalize(fluxes_jet_a[ish].max() * 1e-1, fluxes_jet_a[ish].max())
 #                c = ax_b.scatter(x_jet_a[ish].flatten(), y_jet_a[ish].flatten(), c=cmap_pj(my_norm(fluxes_jet_a[ish].flatten())), alpha=.5)
 #                c = ax_b.scatter(x_jet_a[_nshells+ish].flatten(), y_jet_a[_nshells+ish].flatten(), c=cmap_cj(my_norm(fluxes_jet_a[_nshells+ish].flatten())), alpha=.5)
 #            # cbar = plt.colorbar(ScalarMappable(cmap=c.get_cmap(), norm=c.norm), ax=ax_b)
 #
 #            int_x_j, int_y_j, int_zz_j = combine_images(x_jet_a, y_jet_a, fluxes_jet_a, pre_int=True, k=10, nx=nx, ny=ny, extend=2, verbose=True)
 #
 #            # int_x_j, int_y_j, int_zz_j = combine_images_old(x_jet_a, y_jet_a, fluxes_jet_a,
 #            #                                                 shells=True, hist_or_int='hist', nx=nx, ny=ny, extend=2)
 #
 #            grid_X, grid_Y = np.meshgrid(int_x_j, int_y_j, indexing='ij')
 #            max = 5.5
 #            if my_norm2 is None: my_norm2 = LogNorm(int_zz_j.max() * 1e-1, int_zz_j.max())
 #            # ax_x = fig.add_subplot(len(theta_maxs[0]), 2, (i+1)*2, projection='3d')
 #            ax_x = fig.add_subplot(len(theta_maxs[0]), 3, (i+1)*3)
 #            if int_x_j.ndim == 1:
 #                int_x_j, int_y_j = np.meshgrid(int_x_j, int_y_j, indexing='ij')
 #
 #            cmap = cm.get_cmap('inferno')
 #            c = ax_x.scatter(int_x_j.flatten(), int_y_j.flatten(), c=cmap(my_norm2(int_zz_j.flatten())), marker='s', alpha=.5)
 #            # ax_x.pcolormesh(grid_X, grid_Y, i_zz,  cmap=cmap,norm=my_norm)
 #            # cbar = plt.colorbar(ScalarMappable(cmap=c.get_cmap(), norm=c.norm), ax=ax_x)
 #
 #            # xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(x_jet_a, y_jet_a, fluxes_jet_a)
 #            ax_b.plot([xc_m_j], [yc_m_j], color="lime", marker='o', ms=10)
 #            ax_b.axhline(y=0, color='gray')
 #            ax_b.set_title(f"Tot.Flux={pba_a.GRB.get_skymap_totfluxes(freq=freq,time=time_ * cgs.day)}")
 #
 #            res[f"thmax={theta_maxs[0][i]}"] = f"xmin={x_jet_a[0].min():.2e} xmax={x_jet_a[0].min():.2e} " \
 #                                               f"| ymin={y_jet_a[0].min():.2e} ymax={y_jet_a[0].min():.2e} " \
 #                                               f"| Imin={fluxes_jet_a[0].min():.2e} Imax={fluxes_jet_a[0].max():.2e}" \
 #                                               f"| IHISTmin={int_zz_j.min():.2e} IHISTmax={int_zz_j.max():.2e}" \
 # \
 # \
 # \
 # \
 #                # print(f"Piecewise: I[{fluxes_jet_pw.min(),fluxes_jet_pw.max()}]")
 #            print(f"Adaptive: I[{fluxes_jet_a[0].min(),fluxes_jet_a[0].max()}]")
 #
 #        for key, val in res.items():
 #            print(f"{key} {val}")
 #
 #        # plt.tight_layout()
 #        plt.show()
 #    def compare_skymaps(self, pars : dict, opts_a:dict, opts_pw:dict, struct : dict, title : str,
 #                        resolutions : tuple[tuple],
 #                        figpath : str or None,
 #                        times : np.ndarray, show_fig : bool, save_pdf : bool):
 #
 #        pba_pw = []
 #        pba_a = []
 #        nres_pw = len(resolutions[0]) if hasattr(resolutions[0],'__len__') else 0
 #        nres_a = len(resolutions[1]) if hasattr(resolutions[1],'__len__') else 0
 #        nres = np.max([nres_pw,nres_a])
 #        for i in range(nres):
 #            nlayer_a = resolutions[1][i] if nres_a > 0 else 0
 #            nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0
 #            # ------ Piece Wise -------
 #            if (nres_pw > 0):
 #                struct["nlayers_pw"] = nlayer_pw
 #                opts_pw["do_skymap"] = "yes"
 #                pars_pw = copy.deepcopy(pars)
 #                pars_pw["mom0_frac_when_start_spread"] = 0.1
 #                opts_pw["fname_sky_map"]=f"skymap_pw_res{nlayer_pw}.h5"
 #                pba_pw.append( self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=opts_pw) )
 #            # ------ Adaptive -------
 #            if (nres_a > 0):
 #                struct["nlayers_a"] = nlayer_a[0]
 #                pars["ntheta"] = nlayer_a[1]
 #                pars["nphi"] = nlayer_a[2]
 #                pars["mom0_frac_when_start_spread"] = 0.95
 #                opts_a["do_skymap"] = "yes"
 #                opts_a["fname_sky_map"]=f"skymap_a_res{nlayer_a[0]}-{nlayer_a[1]}-{nlayer_a[2]}.h5"
 #                pba_a.append( self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a) )
 #        print("\tComputed. Plotting...")
 #
 #        freq = pars["obs_freq"]
 #        tmp={"hist_nx": 81, "hist_ny": 81, "spec": False,
 #             "smooth": {},  # {"type": "gaussian", "sigma": 10},
 #             "cm": {"color": 'cyan', "marker": "+"},
 #             "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
 #             "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
 #             "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
 #                            "norm": ("linear", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": 0}
 #             }
 #        settings={
 #            "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
 #            "title": {"title": "time_fluxratio"},  # "time_fluxratio"
 #            "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
 #            "xlabel": "x [mas]", "ylabel": "z [mas]",
 #            "histx_backgound_color": "black",
 #            "histy_backgound_color": "black",
 #            "plot_grids": True,
 #            "histx_lim":(1e-2, 1e1),
 #            "histy_lim":(1e-2, 1e1)
 #        }
 #
 #
 #        for time_ in times:
 #            fig,axes = plt.subplots(ncols=2,nrows=len(resolutions[0])+1,sharex='all',sharey='row',figsize=(5,10))
 #            tot_fluxes = {}
 #            for i in range(nres):
 #                nlayer_a = resolutions[1][i] if nres_a > 0 else 0
 #                nlayer_pw = resolutions[0][i] if nres_pw > 0 else 0
 #                color=resolutions[2][i] if nres_pw > 0 else resolutions[2]
 #                tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color
 #
 #                # ------ Piece Wise -------
 #                if (nres_pw > 0):
 #                    # struct["nlayers_pw"] = nlayer_pw
 #                    # self.opts_pw_1["do_skymap"] = "yes"
 #                    # pars_pw = copy.deepcopy(pars)
 #                    # self.opts_pw_1["fname_sky_map"]=f"skymap_pw_res{nlayer_pw}.h5"
 #                    # pba_pw = self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=self.opts_pw_1)
 #
 #                    # all_x_jet, all_y_jet, all_fluxes_jet \
 #                    #     = pba_pw.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True, remove_zeros=True)
 #                    all_x_jet, all_y_jet, all_fluxes_jet \
 #                        = pba_pw[i].GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=True,renormalize=True)
 #                    int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                pre_int=False, k=5, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=1.1)
 #                    # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet, shells=True, verbose=True,
 #                    #                                                 hist_or_int="hist", nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
 #                    grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                           collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                    grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                           collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                    xc_m_j, yc_m_j = pba_pw[i].GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
 #                    if (nres_a>0 and nres_pw>0):
 #                        ax_main = axes[i+1][0]; ax_histx=axes[0,0]; ax_histy=None
 #                    else:
 #                        ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
 #                    im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
 #                                                grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)
 #                    ax_main.annotate(str(int(nlayer_pw)), color="white", xy=(0.95, 0),xycoords='axes fraction',
 #                                     fontsize=14, horizontalalignment='right', verticalalignment='bottom')
 #                    # ax_cbar = axes[i+1][0]
 #                    # divider = make_axes_locatable(ax_cbar)
 #                    # cax = divider.append_axes('left', size='99%', pad=0.9)
 #                    # plt.delaxes(ax_cbar)
 #                    # cbar = plt.colorbar(im, cax=cax,
 #                    #                     # format='%.0e',ticks=ticks
 #                    #                     orientation='vertical',
 #                    #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
 #                    #                     )
 #
 #                # ------ Adaptive -------
 #                if (nres_a > 0):
 #                    # pars["nsublayers"] = nlayer_a
 #                    # self.opts_a_1["do_skymap"] = "yes"
 #                    # self.opts_a_1["fname_sky_map"]=f"skymap_a_res{nlayer_a}.h5"
 #                    # pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=self.opts_a_1)
 #
 #                    all_x_jet, all_y_jet, all_fluxes_jet \
 #                        = pba_a[i].GRB.get_skymap(time=time_ * cgs.day, freq=freq, verbose=True, remove_mu=False,
 #                                                  renormalize=True, normtype='a', remove_zeros=False)
 #                    int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet, verbose=True,
 #                                                                pre_int=True, k=10, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
 #                    grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                           collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                    grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                           collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                    xc_m_j, yc_m_j = pba_a[i].GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
 #                    if (nres_a>0 and nres_pw>0):
 #                        ax_main = axes[i+1][1]; ax_histx=axes[0,1]; ax_histy=None
 #                    else:
 #                        ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
 #                    im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
 #                                                grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)
 #                    ax_main.annotate(f"{int(nlayer_a[0])}-{int(nlayer_a[1])}-{int(nlayer_a[2])}", color="white", xy=(0.95, 0),xycoords='axes fraction',
 #                                     fontsize=14, horizontalalignment='right', verticalalignment='bottom')
 #
 #                tot_fluxes[f"nlayer_a={nlayer_a[0]}-{nlayer_a[1]}-{nlayer_a[2]} nlayer_pw={nlayer_pw}"] = (
 #                        f" PW={pba_pw[-1].GRB.get_skymap_totfluxes(freq=freq,time=time_* cgs.day)}" + \
 #                        f" lc:({pba_pw[-1].GRB.get_lc(freq=freq,time=time_* cgs.day)})" + \
 #                        f" A={pba_a[-1].GRB.get_skymap_totfluxes(freq=freq,time=time_* cgs.day)}" + \
 #                        f" lc:({pba_a[-1].GRB.get_lc(freq=freq,time=time_* cgs.day)})")
 #
 #                if (i==nres-1):
 #                    ax_main.set_xlabel("")
 #
 #                # --- COLOBAR
 #                # ax_cbar = axes[1,i+1]
 #                # divider = make_axes_locatable(ax_cbar)
 #                # cax = divider.append_axes('right', size='99%', pad=0.9)
 #                # plt.delaxes(ax_cbar)
 #                # cbar = plt.colorbar(im, cax=cax,
 #                #                     # format='%.0e',ticks=ticks
 #                #                     orientation='vertical',
 #                #                     # label=r"$I_{\nu}$ [mJy/mas$^2$]"
 #                #                     )
 #            for key, val in tot_fluxes.items():
 #                print(f"{key} {val}")
 #
 #            if hasattr(axes[0],'__len__'):
 #                axes_ = axes[0,:]
 #                axes[0,0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)  # , color='blue')
 #            else:
 #                axes_ = axes
 #                axes[0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
 #            for i_, ax_histx in enumerate(axes_):
 #                ax_histx.set_yscale("linear")
 #                ax_histx.set_xscale("linear")
 #                ax_histx.minorticks_on()
 #                ax_histx.set_yscale("log")
 #                # # ax_histx.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
 #                if (i_ == 0):
 #                    ax_histx.tick_params(axis='both', which='both', labelleft=True,
 #                                         labelright=False, tick1On=True, tick2On=True,
 #                                         labelsize=12, color="white",
 #                                         direction='in',
 #                                         bottom=True, top=True, left=True, right=True)
 #                else:
 #                    ax_histx.tick_params(axis='both', which='both', labelleft=False,
 #                                         labelright=False, tick1On=True, tick2On=True,
 #                                         labelsize=12, color="white",
 #                                         direction='in',
 #                                         bottom=True, top=True, left=True, right=True)
 #                ax_histx.set_facecolor(settings["histx_backgound_color"])
 #
 #                if settings["plot_grids"]:
 #                    ax_histx.grid()
 #                # ax_histx.set_ylim(*settings["histx_lim"])
 #            # ax_histy.set_yscale("linear")
 #            # ax_histy.set_xscale("linear")
 #            # ax_histy.minorticks_on()
 #            # ax_histy.set_xscale("log")
 #            # ax_histy.set_xlim(ax_histx.get_ylim())
 #            # # ax_histy.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
 #            # ax_histy.set_xlabel(r"$I_{\nu;\rm m}(z)$", fontsize=12)  # , color='blue')
 #            # ax_histy.tick_params(axis='both', which='both', labelleft=False,
 #            #                      labelright=False, tick1On=True, tick2On=True,
 #            #                      labelsize=12, color="white",
 #            #                      direction='in',
 #            #                      bottom=True, top=True, left=True, right=False)
 #            # ax_histy.set_facecolor(settings["histy_backgound_color"])
 #
 #
 #            if hasattr(axes[0],'__len__'):
 #                for i_, axes_main in enumerate(zip(axes[1:,0],axes[1:,1])):
 #                    for j_, ax_main in enumerate(axes_main):
 #                        ax_main.set_yscale("linear")
 #                        ax_main.set_xscale("linear")
 #                        ax_main.minorticks_on()
 #                        ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
 #                        ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
 #                        if (j_ == 0):
 #                            ax_main.tick_params(axis='both', which='both', labelleft=True,
 #                                                labelright=False, tick1On=True, tick2On=True,
 #                                                labelsize=12, color="white",
 #                                                direction='in',
 #                                                bottom=True, top=True, left=True, right=True)
 #                            ax_main.set_ylabel("$Z$ [mas]",fontsize=12)
 #                        else:
 #                            ax_main.tick_params(axis='both', which='both', labelleft=False,
 #                                                labelright=False, tick1On=True, tick2On=True,
 #                                                labelsize=12, color="white",
 #                                                direction='in',
 #                                                bottom=True, top=True, left=True, right=True)
 #                        if settings["plot_grids"]:
 #                            ax_main.grid(color='gray', linestyle='-', linewidth=0.6)
 #                axes[-1,0].set_xlabel("$X$ [mas]",fontsize=12)
 #                axes[-1,1].set_xlabel("$X$ [mas]",fontsize=12)
 #            # plt.tight_layout()
 #            # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
 #            plt.subplots_adjust(left = 0.2, top = 0.95, right = 0.95, bottom = 0.05, hspace = 0.1, wspace = 0.1)
 #            fig.suptitle(f'Time={int(time_)} d.'+r" $\theta_{\rm obs}$="+str(int(pars["theta_obs"]*180/np.pi))+" deg.", fontsize=14)
 #
 #            if not figpath is None:
 #                figpath_=figpath+ f"_t{int(time_)}d"
 #                print("Saving:\n {}".format(figpath_))
 #                if save_pdf: plt.savefig(figpath_+".pdf")
 #                plt.savefig(figpath_+".png", dpi=256)
 #            if show_fig: plt.show()
 #
 #    def compare_skymap_evolution(self, pars : dict, opts_a:dict, opts_pw:dict, struct : dict, title : str,
 #                                 resolutions_pw : tuple[tuple],
 #                                 resolutions_a : tuple[tuple],
 #                                 figpath : str,
 #                                 plots : tuple,
 #                                 show_fig : bool, save_pdf : bool):
 #        if len(plots)==0:
 #            raise ValueError("nothing to plot")
 #        freq = pars["obs_freq"]
 #        tmp = {"hist_nx":91, "hist_ny":91}
 #
 #        times = np.copy(pars["skymap_times"])
 #
 #        res = {}
 #        for ir, layers in enumerate(resolutions_pw[0]):
 #            res[f"pw{layers}"] = {"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),
 #                                  "xsize":np.zeros_like(times), "ysize":np.zeros_like(times),
 #                                  "fnu":np.zeros_like(times)}
 #
 #        for ir, sublayers in enumerate(resolutions_a[0]):
 #            res[f"a{sublayers[0]}-{sublayers[1]}-{sublayers[2]}"] \
 #                = {"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),
 #                   "xsize":np.zeros_like(times), "ysize":np.zeros_like(times),
 #                   "fnu":np.zeros_like(times)}
 #
 #        # Piece wise
 #        for ir, layers in enumerate(resolutions_pw[0]):
 #
 #            struct["nlayers_pw"] = layers
 #            opts_pw["do_skymap"] = "yes"
 #            pars_pw = copy.deepcopy(pars)
 #            pars_pw["mom0_frac_when_start_spread"] = 0.1
 #            pba_pw = self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=opts_pw)
 #
 #            for it, t in enumerate(times):
 #                # Plot peice-wise
 #
 #                # all_x_jet, all_y_jet, all_fluxes_jet = pba_pw.GRB.get_skymap(time=t * cgs.day, freq=freq, verbose=False, remove_mu=False, renormalize=True)
 #                all_x_jet, all_y_jet, all_fluxes_jet = pba_pw.GRB.get_skymap_old(
 #                    time=t * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=False)
 #
 #                # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
 #                #                                             hist_or_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
 #
 #                xc_m_j, yc_m_j = pba_pw.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
 #
 #                grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                       collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
 #                grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                       collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
 #
 #                res[f"pw{layers}"]["xc"][it] = xc_m_j
 #                res[f"pw{layers}"]["yc"][it] = yc_m_j
 #                res[f"pw{layers}"]["xsize"][it] = x2-x1
 #                res[f"pw{layers}"]["ysize"][it] = y2-y1
 #                res[f"pw{layers}"]["fnu"][it] = pba_pw.GRB.get_skymap_totfluxes(freq=freq,time=t*cgs.day)
 #
 #
 #        # adaptive
 #        for ir, nll in enumerate(resolutions_a[0]):
 #
 #            pars["ntheta"] = nll[1]
 #            pars["nphi"] = nll[2]
 #            struct["nlayers_a"] = nll[0]
 #            opts_a["do_skymap"] = "yes"
 #            pba_a = self.run_a(struct=struct, pars=pars, opts={}, opts_grb=opts_a)
 #
 #            for it, t in enumerate(times):
 #
 #                all_x_jet, all_y_jet, all_fluxes_jet \
 #                    = pba_a.GRB.get_skymap(
 #                    time=t * cgs.day, freq=freq, verbose=False, remove_mu=False, renormalize=False, normtype='a')
 #                # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
 #                #                                             pre_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
 #
 #                xc_m_j, yc_m_j = pba_a.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
 #
 #                grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                       collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
 #
 #                grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
 #                                                                       collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
 #                x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
 #
 #                res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xc"][it] = xc_m_j
 #                res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["yc"][it] = yc_m_j
 #                res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xsize"][it] = x2-x1
 #                res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["ysize"][it] = y2-y1
 #                res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["fnu"][it] = pba_a.GRB.get_skymap_totfluxes(freq=freq,time=t*cgs.day)
 #
 #        # plotting
 #        fig, axes=plt.subplots(ncols=1,nrows=len(plots),sharex='all',figsize=(6,7))
 #        for nlayers, color, lbl in zip(resolutions_pw[0],resolutions_pw[1],resolutions_pw[2]):
 #            p = 0
 #            if ("xc" in plots):
 #                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["xc"], color=color, ls='-', lw=2, alpha=1, marker='x')
 #                p=p+1
 #            if ("yc" in plots):
 #                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["yc"], color=color, ls='-', lw=2, alpha=1, marker='x')
 #                p=p+1
 #            if ("xsize" in plots):
 #                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["xsize"], color=color, ls='-', lw=2, alpha=1, marker='x')
 #                p=p+1
 #            if ("ysize" in plots):
 #                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["ysize"], color=color, ls='-', lw=2, alpha=1, marker='x')
 #                p=p+1
 #            if ("fluxes" in plots):
 #                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["fnu"], color=color, ls='-', lw=2, alpha=1, marker='x', label=f"PW{nlayers}")
 #                #p=p+1
 #        for nll, color, lbl in zip(resolutions_a[0],resolutions_a[1],resolutions_a[2]):
 #            p = 0
 #            if ("xc" in plots):
 #                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
 #                p=p+1
 #            if ("yc" in plots):
 #                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["yc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
 #                p=p+1
 #            if ("xsize" in plots):
 #                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xsize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
 #                p=p+1
 #            if ("ysize" in plots):
 #                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["ysize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
 #                p=p+1
 #            if ("fluxes" in plots):
 #                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["fnu"], color=color, ls='-', lw=2, alpha=.7, marker='o',
 #                             label=f"A{nll[0]}-{nll[1]}-{nll[2]}")
 #                #p=p+1
 #        labels = {"xc":r"$X_{c}$ [mas]", "yc":r"$Z_{c}$ [mas]", "xsize":r"FWHM$_{x}$", "ysize":r"FWHM$_{z}$","fluxes":r"$F_{\nu}$ [mJy]"}
 #        for ax, key in zip(axes,plots):
 #            if (key in plots):
 #                ax.set_ylabel(labels[key])
 #            if (key == "fluxes"):
 #                axes[p].set_yscale('log')
 #                ax.legend(fancybox=True, loc='lower right',
 #                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
 #                          shadow=False, ncol=2, fontsize=12,
 #                          framealpha=0., borderaxespad=0.)
 #
 #            ax.minorticks_on()
 #            ax.tick_params(axis='both', which='both', labelleft=True,
 #                           labelright=False, tick1On=True, tick2On=True,
 #                           labelsize=12,
 #                           direction='in',
 #                           bottom=True, top=True, left=True, right=True)
 #
 #            # ax.set_xscale("log")
 #
 #
 #        # axes[0].set_ylim(-2,2)
 #        # axes[1].set_ylim(-2,2)
 #        axes[len(axes)-1].set_xlabel("time [d]")
 #
 #        plt.tight_layout()
 #
 #        if not figpath is None:
 #            print("Saving:\n {}".format(figpath))
 #            if save_pdf: plt.savefig(figpath+".pdf")
 #            plt.savefig(figpath+".png", dpi=256)
 #        if show_fig: plt.show()
    def compare_skymap_evolution_plot_each_timestep(
            self, pars : dict, opts_a:dict, opts_pw:dict, struct : dict, title : str,
            resolutions_pw : tuple[tuple], resolutions_a : tuple[tuple],
            figfpath_skymap : str, figfpath_evolve : str, fpath_out : str,
            plots : tuple, show_fig : bool, save_pdf : bool, save_res : bool):
        if len(plots)==0:
            raise ValueError("nothing to plot")
        freq = pars["obs_freq"]

        tmp={
            "type":"hist",
            "cm": {"color": 'cyan', "marker": "+"},
            "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                           "norm": ("linear", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": 0}
        }
        settings={
            "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
            "title": {"title": "time_fluxratio"},  # "time_fluxratio"
            "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
            "xlabel": "x [mas]", "ylabel": "z [mas]",
            "histx_backgound_color": "black",
            "histy_backgound_color": "black",
            "plot_grids": True,
            "histx_lim":(1e-2, 1e1),
            "histy_lim":(1e-2, 1e1)
        }

        times = np.copy(pars["skymap_times"])

        res = {}
        for ir, layers in enumerate(resolutions_pw[0]):
            res[f"pw{layers}"] = {"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),
                                  "xsize":np.zeros_like(times), "ysize":np.zeros_like(times),
                                  "fnu":np.zeros_like(times)
                                  # "int_x_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                                  # "int_y_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                                  # "int_zz_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                                  }

        for ir, sublayers in enumerate(resolutions_a[0]):
            res[f"a{sublayers[0]}-{sublayers[1]}-{sublayers[2]}"] \
                = {"t":times,"xc":np.zeros_like(times),"yc":np.zeros_like(times),
                   "xsize":np.zeros_like(times), "ysize":np.zeros_like(times),
                   "fnu":np.zeros_like(times)
                   # "int_x_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                   # "int_y_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                   # "int_zz_j":np.zeros((len(times),tmp["hist_nx"]-1,tmp["hist_ny"]-1)),#[np.empty((tmp["hist_nx"],tmp["hist_ny"])) for i in range(len(times))],
                   }

        pba_pw = []
        pba_a = []

        # --- RUN ---
        # Piece wise
        for ir, layers in enumerate(resolutions_pw[0]):
            struct["nlayers_pw"] = layers
            _pars = copy.deepcopy(pars)
            _opts = copy.deepcopy(opts_pw)
            _opts["do_skymap"] = "yes"
            _pars["mom0_frac_when_start_spread"] = 0.1
            _opts["fname_sky_map"]=f"skymap_pw_res{layers}.h5"
            pba_pw.append(
                self.run_pw(struct=struct, pars=_pars, opts={}, opts_grb=_opts)
            )
        # adaptive
        for ir, nll in enumerate(resolutions_a[0]):
            struct["nlayers_a"] = nll[0]
            _pars = copy.deepcopy(pars)
            _opts_a = copy.deepcopy(opts_a)
            _pars["ntheta"] = nll[1]
            _pars["nphi"] = nll[2]
            _pars["mom0_frac_when_start_spread"] = 0.95
            _opts_a["do_skymap"] = "yes"
            _opts_a["fname_sky_map"]=f"skymap_a_res{nll[0]}-{nll[1]}-{nll[2]}.h5"
            pba_a.append(
                self.run_a(struct=struct, pars=_pars, opts={}, opts_grb=_opts_a)
            )

        # --- EXTRACT SKYAP DATA & PLOT SKYMAPS  ---
        nres_pw = len(resolutions_a[0]) if hasattr(resolutions_a[0],'__len__') else 0
        nres_a = len(resolutions_a[1]) if hasattr(resolutions_a[1],'__len__') else 0
        nres = np.max([nres_pw,nres_a])
        for it, t in enumerate(times):
            fig,axes = plt.subplots(ncols=2,nrows=len(resolutions_a[0])+1,sharex='all',sharey='row',figsize=(5,10))
            tot_fluxes = {}
            for i, layers in enumerate(resolutions_pw[0]):

                color=resolutions_pw[1][i] if nres_pw > 0 else resolutions_pw[1]
                tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color
                nlayer_pw = resolutions_pw[0][i] if nres_pw > 0 else 0

                # ------ Piece Wise -------
                if (nres_pw > 0):
                    # # struct["nlayers_pw"] = nlayer_pw
                    # # self.opts_pw_1["do_skymap"] = "yes"
                    # # pars_pw = copy.deepcopy(pars)
                    # # self.opts_pw_1["fname_sky_map"]=f"skymap_pw_res{nlayer_pw}.h5"
                    # # pba_pw = self.run_pw(struct=struct, pars=pars_pw, opts={}, opts_grb=self.opts_pw_1)
                    #
                    # # all_x_jet, all_y_jet, all_fluxes_jet \
                    # #     = pba_pw.GRB.get_skymap_old(time=time_ * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=True, remove_zeros=True)
                    # # all_x_jet, all_y_jet, all_fluxes_jet \
                    # #     = pba_pw[i].GRB.get_skymap_old(time=t * cgs.day, freq=freq, verbose=False, remove_mu=False,renormalize=False)
                    # all_x_jet, all_y_jet, all_fluxes_jet \
                    #     = pba_pw[i].GRB.get_skymap(time=t * cgs.day, freq=freq, verbose=False)
                    # # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
                    # #                                             pre_int=True, k=5, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=1.1)
                    #
                    # xc_m_j, yc_m_j = pba_pw[i].GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
                    # # xc_m_j_ = pba_pw[i].GRB.get_skymap_property(v_n="xc", freq=freq, time=t*cgs.day)
                    # # yc_m_j_ = pba_pw[i].GRB.get_skymap_property(v_n="yc", freq=freq, time=t*cgs.day)
                    #
                    # # int_x_j, int_y_j, int_zz_j = combine_images_old(all_x_jet, all_y_jet, all_fluxes_jet, shells=True, verbose=True,
                    # #                                                 hist_or_int="hist", nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
                    # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet, verbose=True,
                    #                                                 pre_int=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
                    # grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                    #                                                        collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
                    # y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
                    #
                    # grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                    #                                                        collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
                    # x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
                    #
                    if (nres_a>0 and nres_pw>0):
                        ax_main = axes[i+1][0]; ax_histx=axes[0,0]; ax_histy=None
                    else:
                        ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
                    # im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
                    #                              grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)

                    skymap = pba_pw[i].GRB.get_skymap(time=t*cgs.day, freq=freq)
                    im = plot_skymap_with_hists(
                        skymap=skymap, tmp=tmp, ax_main=ax_main, ax_histx=ax_histx, ax_histy=ax_histy
                    )
                    ax_main.annotate(str(int(nlayer_pw)), color="white", xy=(0.95, 0),xycoords='axes fraction',
                                     fontsize=14, horizontalalignment='right', verticalalignment='bottom')

                    # save data for evol. plot
                    res[f"pw{layers}"]["xc"][it] = skymap.xc
                    res[f"pw{layers}"]["yc"][it] = skymap.yc
                    res[f"pw{layers}"]["xsize"][it] = skymap.x2-skymap.x1
                    res[f"pw{layers}"]["ysize"][it] = skymap.y2-skymap.y1
                    res[f"pw{layers}"]["fnu"][it] = skymap.flux
                    # res[f"pw{layers}"]["int_x_j"][it,:,:] = np.copy(int_x_j)
                    # res[f"pw{layers}"]["int_y_j"][it,:,:] = np.copy(int_y_j)
                    # res[f"pw{layers}"]["int_zz_j"][it,:,:] = np.copy(int_zz_j)

            for i, nll in enumerate(resolutions_a[0]):
                color=resolutions_a[1][i] if nres_pw > 0 else resolutions_a[1]
                tmp["cm"]["color"]=color; tmp["ysize"]["color"]=color; tmp["xsize"]["color"]=color
                # ------ Adaptive -------
                if (nres_a > 0):
                    # all_x_jet, all_y_jet, all_fluxes_jet \
                    #     = pba_a[i].GRB.get_skymap(time=t * cgs.day, freq=freq, verbose=True, remove_zeros=True)
                    #
                    # xc_m_j, yc_m_j = pba_a[i].GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
                    # # xc_m_j_ = pba_a[i].GRB.get_skymap_property(v_n="xc", freq=freq, time=t*cgs.day)
                    # # yc_m_j_ = pba_a[i].GRB.get_skymap_property(v_n="yc", freq=freq, time=t*cgs.day)
                    #
                    # int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet, verbose=True,
                    #                                             pre_int=False, k=10, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
                    # # int_x_j, int_y_j, int_zz_j = combine_images_old(all_x_jet, all_y_jet, all_fluxes_jet, shells=True, verbose=True,
                    # #                                                 hist_or_int="hist", nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
                    # grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                    #                                                        collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
                    # y1, y2 = get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
                    #
                    # grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                    #                                                        collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
                    # x1, x2 = get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
                    #
                    if (nres_a>0 and nres_pw>0):
                        ax_main = axes[i+1][1]; ax_histx=axes[0,1]; ax_histy=None
                    else:
                        ax_main = axes[i+1]; ax_histx=axes[0]; ax_histy=None
                    # im = plot_skymap_with_hists(ax_main, ax_histx, ax_histy, _, tmp,
                    #                              grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)


                    skymap = pba_a[i].GRB.get_skymap(time=t*cgs.day, freq=freq)
                    im = plot_skymap_with_hists(
                        skymap=skymap, tmp=tmp, ax_main=ax_main, ax_histx=ax_histx, ax_histy=ax_histy
                    )
                    ax_main.annotate(f"{int(nll[0])}-{int(nll[1])}-{int(nll[2])}",
                                     color="white", xy=(0.95, 0),xycoords='axes fraction',
                                     fontsize=14, horizontalalignment='right', verticalalignment='bottom')

                    res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xc"][it] = skymap.xc
                    res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["yc"][it] = skymap.yc
                    res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xsize"][it] = skymap.x2-skymap.x1
                    res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["ysize"][it] = skymap.y2-skymap.y1
                    res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["fnu"][it] = skymap.flux
                    # res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["int_x_j"][it,:,:] = np.copy(int_x_j)
                    # res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["int_y_j"][it,:,:] =  np.copy(int_y_j)
                    # res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["int_zz_j"][it,:,:] =  np.copy(int_zz_j)
                    # second column should not have xlabels in top plot

                    ax_main.set_xlim(skymap.grid_x.min()*2, skymap.grid_x.max()*2)

                    if (i==nres-1):
                        ax_main.set_xlabel("")


            # collect also total fluxes for each image
            # if (len(resolutions_a[0]) == len(resolutions_pw[0])):
            #     for i in range(len(resolutions_a[0])):
            #         tot_fluxes[f"nlayer_a={resolutions_a[i][0]}-{resolutions_a[i][1]}-{resolutions_a[i][2]} nlayer_pw={resolutions_pw[i][0]}"] = (
            #         f" PW={pba_pw[i].GRB.get_skymap_totfluxes(freq=freq,time=t* cgs.day)}" + \
            #             f" lc:({pba_pw[i].GRB.get_lc(freq=freq,time=t * cgs.day)})" + \
            #             f"   A={pba_a[i].GRB.get_skymap_totfluxes(freq=freq,time=t* cgs.day)}" + \
            #             f" lc:({pba_a[i].GRB.get_lc(freq=freq,time=t * cgs.day)})")

            print("------------------------------------------------")
            print(f"time={t}")
            for key, val in tot_fluxes.items():
                print(f"\t {key} {val}")

            if hasattr(axes[0],'__len__'):
                axes_ = axes[0,:]
                axes[0,0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)  # , color='blue')
            else:
                axes_ = axes
                axes[0].set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)
            for i_, ax_histx in enumerate(axes_):
                ax_histx.set_yscale("linear")
                ax_histx.set_xscale("linear")
                ax_histx.minorticks_on()
                ax_histx.set_yscale("log")
                # # ax_histx.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
                if (i_ == 0):
                    ax_histx.tick_params(axis='both', which='both', labelleft=True,
                                         labelright=False, tick1On=True, tick2On=True,
                                         labelsize=12, color="white",
                                         direction='in',
                                         bottom=True, top=True, left=True, right=True)
                else:
                    ax_histx.tick_params(axis='both', which='both', labelleft=False,
                                         labelright=False, tick1On=True, tick2On=True,
                                         labelsize=12, color="white",
                                         direction='in',
                                         bottom=True, top=True, left=True, right=True)
                ax_histx.set_facecolor(settings["histx_backgound_color"])

                if settings["plot_grids"]:
                    ax_histx.grid()
                # ax_histx.set_ylim(*settings["histx_lim"])

            if hasattr(axes[0],'__len__'):
                for i_, axes_main in enumerate(zip(axes[1:,0],axes[1:,1])):
                    for j_, ax_main in enumerate(axes_main):
                        ax_main.set_yscale("linear")
                        ax_main.set_xscale("linear")
                        ax_main.minorticks_on()
                        ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
                        ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
                        if (j_ == 0):
                            ax_main.tick_params(axis='both', which='both', labelleft=True,
                                                labelright=False, tick1On=True, tick2On=True,
                                                labelsize=12, color="white",
                                                direction='in',
                                                bottom=True, top=True, left=True, right=True)
                            ax_main.set_ylabel("$Z$ [mas]",fontsize=12)
                        else:
                            ax_main.tick_params(axis='both', which='both', labelleft=False,
                                                labelright=False, tick1On=True, tick2On=True,
                                                labelsize=12, color="white",
                                                direction='in',
                                                bottom=True, top=True, left=True, right=True)
                        if settings["plot_grids"]:
                            ax_main.grid(color='gray', linestyle='-', linewidth=0.6)

                        ax_main.set_xlim(skymap.grid_x.min()*2, skymap.grid_x.max()*2)
                        ax_main.set_ylim(skymap.grid_y.min()*2, skymap.grid_y.max()*2)

                axes[-1,0].set_xlabel("$X$ [mas]",fontsize=12)
                axes[-1,1].set_xlabel("$X$ [mas]",fontsize=12)


            plt.subplots_adjust(left = 0.2, top = 0.95, right = 0.95, bottom = 0.05, hspace = 0.1, wspace = 0.1)
            fig.suptitle(f'Time={int(t)} d.'+r" $\theta_{\rm obs}$="+str(int(pars["theta_obs"]*180/np.pi))+" deg.", fontsize=14)

            if not figfpath_skymap is None:
                figpath_=figfpath_skymap+ f"_t{int(t)}d"
                print("Saving:\n {}".format(figpath_))
                if save_pdf: plt.savefig(figpath_+".pdf")
                plt.savefig(figpath_+".png", dpi=256)
            if show_fig: plt.show()


        # --- PLOT EVOLUTION SKYMAP PROPERTIES ---

        fig, axes=plt.subplots(ncols=1,nrows=len(plots),sharex='all',figsize=(6,7))
        for nlayers, color in zip(resolutions_pw[0],resolutions_pw[1]):
            # color = "black"
            p = 0
            if ("xc" in plots):
                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["xc"], color=color, ls=':', lw=2, alpha=1, marker='x')
                p=p+1
            if ("yc" in plots):
                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["yc"], color=color, ls=':', lw=2, alpha=1, marker='x')
                p=p+1
            if ("xsize" in plots):
                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["xsize"], color=color, ls=':', lw=2, alpha=1, marker='x')
                p=p+1
            if ("ysize" in plots):
                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["ysize"], color=color, ls=':', lw=2, alpha=1, marker='x')
                p=p+1
            if ("fluxes" in plots):
                axes[p].plot(res[f"pw{nlayers}"]["t"], res[f"pw{nlayers}"]["fnu"], color=color, ls=':', lw=2, alpha=1, marker='x', label=f"PW{nlayers}")
                #p=p+1
        for nll, color in zip(resolutions_a[0],resolutions_a[1]):
            p = 0
            if ("xc" in plots):
                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
                p=p+1
            if ("yc" in plots):
                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["yc"], color=color, ls='-', lw=2, alpha=.7, marker='o')
                p=p+1
            if ("xsize" in plots):
                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["xsize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
                p=p+1
            if ("ysize" in plots):
                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["ysize"], color=color, ls='-', lw=2, alpha=.7, marker='o')
                p=p+1
            if ("fluxes" in plots):
                axes[p].plot(res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["t"], res[f"a{nll[0]}-{nll[1]}-{nll[2]}"]["fnu"], color=color, ls='-', lw=2, alpha=.7, marker='o',
                             label=f"A{nll[0]}-{nll[1]}-{nll[2]}")
                #p=p+1
        labels = {"xc":r"$X_{c}$ [mas]", "yc":r"$Z_{c}$ [mas]", "xsize":r"FWHM$_{x}$", "ysize":r"FWHM$_{z}$","fluxes":r"$F_{\nu}$ [mJy]"}
        for ax, key in zip(axes,plots):
            if (key in plots):
                ax.set_ylabel(labels[key])
            if (key == "fluxes"):
                axes[p].set_yscale('log')
                ax.legend(fancybox=True, loc='lower right',
                          # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                          shadow=False, ncol=2, fontsize=12,
                          framealpha=0., borderaxespad=0.)

            ax.minorticks_on()
            ax.tick_params(axis='both', which='both', labelleft=True,
                           labelright=False, tick1On=True, tick2On=True,
                           labelsize=12,
                           direction='in',
                           bottom=True, top=True, left=True, right=True)


        # axes[0].set_ylim(-2,2)
        # axes[1].set_ylim(-2,2)
        axes[len(axes)-1].set_xlabel("time [d]")

        plt.tight_layout()

        if not figfpath_evolve is None:
            print("Saving:\n {}".format(figfpath_evolve))
            if save_pdf: plt.savefig(figfpath_evolve+".pdf")
            plt.savefig(figfpath_evolve+".png", dpi=256)
        if show_fig: plt.show()
        if save_res:
            fname = self.workingdir+"out/"+fpath_out+".h5"
            print(f"saving out_file: {fname}")
            dfile = h5py.File(name=fname, mode='w')
            for ires, resolution in enumerate(res.keys()):
                group = dfile.create_group(name=resolution)
                for ik, key in enumerate(res[resolution]):
                    group.create_dataset(name=key, data=res[resolution][key])
            dfile.close()


class CasesFSRS(Base):
    # name_1 = "Tophat off-axis"
    # struct_1 = {
    #     "struct":"tophat",
    #     "Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,
    #     "theta_c": 0.1, "theta_w": 0.1, "nlayers_pw": 150, "nlayers_a": 1
    # }
    # pars_1 = {
    #     "obs_freq":3e9,
    #     "n_ism": 1e-2,      "eps_e": 0.1,   "eps_e_rs": 0.2,
    #     "d_l": 3.09e26,     "eps_b": 0.01,  "eps_b_rs": 0.02,
    #     "z": 0.028,         "p": 2.2,       "p_rs": 2.4,
    #     "theta_obs": 0.2,
    # }
    # opts_1 = {"rtol": 1e-6, "ntb":10000, "iout": 10}
    # opts_a_1 = {"method_eats":"adaptive", "method_spread":"AFGPY", "rhs_type":"bw_fs", "do_rs": "no"}

    def __init__(self, default_parfile_fpath, workingdir):
        super().__init__(default_parfile_fpath, workingdir)

    def paper_plot_compare_fsrs(self, pars : dict, struct : dict, layers : tuple, figfpath : str, save_pdf : bool, show_fig : bool):
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_0deg_layer.h5")

        fig, axes = plt.subplots(figsize=(6,6.5), ncols=1, nrows=3)

        opts_ = copy.deepcopy(self.opts_a_1)
        opts_["bw_type"] = "fs"
        opts_["do_rs"] = "no"
        opts_["fname_dyn"] = f"dyn_a_fs.h5"
        opts_["fname_sky_map"] = f"skymap_a_fs.h5"
        opts_["fname_light_curve"] = f"lc_a_fs.h5"
        opts_["fname_light_curve_layers"] = f"lc_a_fs_layers.h5"
        pba_a = self.run_a(struct=struct, pars=pars, opts=self.opts_1, opts_grb=opts_)
        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"red", "label":r"R19"},
                      plot_layer={"ls":'-.', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":'-', "cmap":"Reds", "alpha":.9, "vmin":-50, "vmax":60})

        # ---------------------------------------------
        opts_ = copy.deepcopy(self.opts_a_1)
        opts_["bw_type"] = "fsrs"
        opts_["do_rs"] = "yes"
        opts_["fname_dyn"] = f"dyn_a_fsrs.h5"
        opts_["fname_sky_map"] = f"skymap_a_fsrs.h5"
        opts_["fname_light_curve"] = f"lc_a_fsrs.h5"
        opts_["fname_light_curve_layers"] = f"lc_a_fsrs_layers.h5"
        pba_a = self.run_a(struct=struct, pars=pars, opts=self.opts_1, opts_grb=opts_)
        # ---------------------------------------------

        self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a, layers = layers,
                      plot={"ls":'-', "color":"blue", "label":"GP12"},
                      plot_layer={"ls":'-.', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[1], pba=pba_a, v_n_x="tburst", v_n_y="mom", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})
        self.plot_dyn(axes[2], pba=pba_a, v_n_x="tburst", v_n_y="theta", layers=layers,
                      plot_layer={"ls":':', "cmap":"Blues", "alpha":.9, "vmin":-50, "vmax":60})


        # -------- Afterglopy --------------
        Z = {'jetType':     grb.jet.TopHat,     # Top-Hat jet
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
        ax.set_xlim(1e5,1e8)
        ax.set_ylim(1e-2,2e2)
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
        # print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        # plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        # plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        # plt.show()
        if not figfpath is None:
            print("Saving:\n {}".format(figfpath))
            if save_pdf: plt.savefig(figfpath+".pdf")
            plt.savefig(figfpath+".png", dpi=256)
        if show_fig: plt.show()

    def paper_plot_resolution_rs(self, pars : dict, opts_a : dict, struct : dict, layers : tuple, resolutions_a : tuple[tuple],
                                 figfpath : str, save_pdf : bool, show_fig : bool):
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_0deg_layer.h5")

        fig, ax = plt.subplots(figsize=(6,6.5), ncols=1, nrows=1)

        # ---------------------------------------------
        pba_a = []
        for ir, nll in enumerate(resolutions_a[0]):
            struct["nlayers_a"] = nll[0]
            _pars = copy.deepcopy(pars)
            _opts_a = copy.deepcopy(opts_a)
            _pars["ntheta"] = nll[1]
            _pars["nphi"] = nll[2]
            _pars["mom0_frac_when_start_spread"] = 0.95
            _opts_a["do_rs"] = "yes"
            _opts_a["bw_type"] = "fsrs"
            _opts_a["fname_dyn"]=f"dyn_fsrs_res{nll[0]}-{nll[1]}-{nll[2]}.h5"
            _opts_a["fname_light_curve"]=f"lc_fsrs_res{nll[0]}-{nll[1]}-{nll[2]}.h5"
            _opts_a["fname_light_curve_layers"]=f"lc_dense_fsrs{nll[0]}-{nll[1]}-{nll[2]}.h5"
            pba_a.append(
                self.run_a(struct=struct, pars=_pars, opts={}, opts_grb=_opts_a)
            )
            color = resolutions_a[1][ir]
            self.plot_lcs(ax=ax, pars=pars, pba=pba_a[ir], layers = layers,
                          plot={"ls":'-', "color":color, "label":f"{int(nll[0])}-{int(nll[1])}-{int(nll[2])}"},
                          plot_layer={"ls":'-.', "color":color, "alpha":.9, "vmin":-50, "vmax":60})


        ax.legend(fancybox=True, loc='upper left',
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=1, fontsize=12,
                  framealpha=0., borderaxespad=0.)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("Flux density [mJy]", fontsize=12)
        # ax.set_title(title)
        ax.set_xlim(1e5,1e8)
        ax.set_ylim(1e-2,2e2)
        # ax.grid()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        # ax.set_facecolor("pink")
        ax.set_facecolor("gray")
        plt.tight_layout()

        # print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        # plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        # plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        # plt.show()
        if not figfpath is None:
            print("Saving:\n {}".format(figfpath))
            if save_pdf: plt.savefig(figfpath+".pdf")
            plt.savefig(figfpath+".png", dpi=256)
        if show_fig: plt.show()

    def paper_plot_resolution2_rs(self, pars : dict, opts_a : dict, struct : dict, resolutions_a : tuple[tuple],
                                  figfpath : str, save_pdf : bool, show_fig : bool) -> dict:
        # run the code for given pars

        # ref = RefData(workdir=curdir, fname="reference_afgpy_dyn.h5")
        # ref_lc = RefDataLC(workdir=curdir, fname="reference_lc_0deg_layer.h5")

        fig, axes = plt.subplots(figsize=(5,3*2), ncols=1, nrows=2, sharex="all")

        # ---------------------------------------------
        pba_a = []
        lightcurves = {}
        for ir, nll in enumerate(resolutions_a[0]):
            struct["nlayers_a"] = nll[0]
            _pars = copy.deepcopy(pars)
            _opts_a = copy.deepcopy(opts_a)
            _pars["ntheta"] = nll[1]
            _pars["nphi"] = nll[2]
            _pars["mom0_frac_when_start_spread"] = 0.95
            _opts_a["do_rs"] = "yes"
            _opts_a["bw_type"] = "fsrs"
            _opts_a["fname_dyn"]=f"dyn_fsrs_res{nll[0]}-{nll[1]}-{nll[2]}.h5"
            _opts_a["fname_light_curve"]=f"lc_fsrs_res{nll[0]}-{nll[1]}-{nll[2]}.h5"
            _opts_a["fname_light_curve_layers"]=f"lc_dense_fsrs{nll[0]}-{nll[1]}-{nll[2]}.h5"
            pba_a.append(
                self.run_a(struct=struct, pars=_pars, opts={}, opts_grb=_opts_a)
            )
            lightcurves[f"il={nll[0]}"] = pba_a[-1].GRB.get_lc_totalflux(freq=pars["obs_freq"], time=None, spec=False)
        lightcurves["time"] = pba_a[-1].GRB.get_lc_times(unique=True,spec=False)

        times = pba_a[-1].GRB.get_lc_times(unique=True,spec=False)
        fluxes = pba_a[-1].GRB.get_lc_totalflux(freq=pars["obs_freq"], time=None, spec=False)

        ax = axes[0]
        for ir, nll in enumerate(resolutions_a[0]):
            color = resolutions_a[1][ir]
            fluxes_i = pba_a[ir].GRB.get_lc_totalflux(freq=pars["obs_freq"], time=None, spec=False)
            ax.plot(times, fluxes_i, color=color, ls='-', alpha=.9,label=f"N={nll[0]}")

        ax.legend(fancybox=True, loc='upper left',
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=1, fontsize=12,
                  framealpha=0., borderaxespad=0.)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel(r"$F_{\nu}$ [mJy]", fontsize=12)
        # ax.set_title(title)
        ax.set_xlim(times.min(),times.max())
        ax.set_ylim(fluxes.max()*1e-4, fluxes.max()*1.1)
        ax.grid()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()


        ax = axes[1]
        # ratios = []
        # ress = []
        for ir, nll in enumerate(resolutions_a[0][:-1]):
            color = resolutions_a[1][ir]
            fluxes_i = pba_a[ir].GRB.get_lc_totalflux(freq=pars["obs_freq"], time=None, spec=False)
            ratio = (fluxes - fluxes_i) / fluxes
            rms = mean_squared_error(fluxes, fluxes_i, squared=False)
            # ratios.append(ratio)
            ax.plot(times, ratio, color=color, ls='-', alpha=.9,label=f"N={nll[0]} RMS={rms:.2e}")
            # ress.append(nll[0])


        ax.legend(fancybox=True, loc='upper left',
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=1, fontsize=12,
                  framealpha=0., borderaxespad=0.)
        ax.set_xscale("log")
        ax.set_yscale("linear")
        ax.set_xlabel("time [s]", fontsize=12)
        ax.set_ylabel("$\Delta F$", fontsize=12)
        # ax.set_title(title)
        ax.set_xlim(times.min(),times.max())
        ax.set_ylim(-1., 1.)
        ax.grid()
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        plt.tight_layout()

        # print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        # plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        # plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        # plt.show()
        if not figfpath is None:
            print("Saving:\n {}".format(figfpath))
            if save_pdf: plt.savefig(figfpath+".pdf")
            plt.savefig(figfpath+".png", dpi=256)
        if show_fig: plt.show()

        return lightcurves


    def compare_grbs(self, pars : dict, opts_a : dict, struct : dict, layers:tuple,
                     setups : tuple[dict], figfpath : str, save_pdf : bool, show_fig : bool):
        """
        ({"n_ism":1e-4,"color":"blue","label":r"$n_{\rm ISM}=$"+"$10^{-4}$ cm$^{-3}$"},
         {"n_ism":1e0,"color":"red","label":r"$n_{\rm ISM}=$"+"$1$ cm$^{-3}$"})
        :param pars:
        :param opts_a:
        :param struct:
        :param vary:
        :param figfpath:
        :param save_pdf:
        :param show_fig:
        :return:
        """
        fig, axes = plt.subplots(figsize=(6,6.5), ncols=1, nrows=3)
        pba_a = []
        for isetup, setup in enumerate(setups):
            _pars = copy.deepcopy(pars)
            _opts_a = copy.deepcopy(opts_a)
            _pars["ntb"] = 3000
            _pars["mom0_frac_when_start_spread"] = 0.95
            _opts_a["do_rs"] = "yes"
            _opts_a["bw_type"] = "fsrs"
            _opts_a["fname_dyn"]=f"dyn_fsrs_setup{isetup}.h5"
            _opts_a["fname_light_curve"]=f"lc_fsrs_setup{isetup}.h5"
            _opts_a["fname_light_curve_layers"]=f"lc_dense_fsrs_setup{isetup}.h5"
            for __par in setup.keys():
                if __par in struct.keys(): struct[__par] = setup[__par]
                if __par in pars.keys(): pars[__par] = setup[__par]
            pba_a.append(
                self.run_a(struct=struct, pars=_pars, opts={}, opts_grb=_opts_a)
            )
            color = setup["color"]
            label = setup["label"]
            cmap = setup["cmap"]
            self.plot_lcs(ax=axes[0], pars=pars, pba=pba_a[isetup], layers = layers,
                          plot={"ls":'-', "color":color, "label":label},
                          plot_layer={"ls":'-.', "color":color, "alpha":.9, "vmin":-50, "vmax":60})
            self.plot_dyn(axes[1], pba=pba_a[isetup], v_n_x="tburst", v_n_y="GammaFsh", layers=layers,
                          plot_layer={"ls":'-', "cmap":cmap, "alpha":.9, "vmin":-50, "vmax":60})
            self.plot_dyn(axes[2], pba=pba_a[isetup], v_n_x="tburst", v_n_y="GammaRsh", layers=layers,
                          plot_layer={"ls":'-', "cmap":cmap, "alpha":.9, "vmin":-50, "vmax":60})

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
        ax.set_xlim(1e5,1e8)
        ax.set_ylim(1e-5,1e-1)
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
        ax.set_ylabel(r"$\Gamma_{\rm fsh}$", fontsize=12)
        ax.set_xlim(1e6,1e11)
        ax.set_ylim(1e0,1e3)
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
        ax.set_yscale("log")
        ax.set_ylabel(r"$\Gamma_{\rm rsh}$", fontsize=12)
        ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
        ax.set_xlim(1.e6,1e11)
        ax.set_ylim(1.e0,100.)
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.minorticks_on()
        # ax.set_facecolor("pink")

        plt.tight_layout()
        # print("Saving:\n {}".format(paperfigdir+"abstract_spread_lcs_dyn.pdf"))
        # plt.savefig(paperfigdir+"abstract_spread_lcs_dyn.pdf")
        # plt.savefig(paperfigdir+"abstract_gaus_spread_lcs_dyn.png", dpi=256)
        # plt.show()
        if not figfpath is None:
            print("Saving:\n {}".format(figfpath))
            if save_pdf: plt.savefig(figfpath+".pdf")
            plt.savefig(figfpath+".png", dpi=256)
        if show_fig: plt.show()

# -----------------------------

def runset_for_skymap(tsk, default_parfile_fpath : str, workingdir : str, figdir : str):

    # parfiledir = os.getcwd().replace("tophat","")
    # parfiledir = os.getcwd().replace("structured","")

    # ------------------- SKYMAPS WITHOUT SPREADING ---------------

    grb = CasesFS(default_parfile_fpath=default_parfile_fpath, workingdir=workingdir)

    tsk.pars["theta_obs"] = 0.#0.785
    tsk.opts_a["method_spread"] = "None"
    tsk.opts_pw["method_spread"] = "None"
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_nospread_0deg_res",
    #                     # times=np.array([5000,]),
    #                     show_fig=False,save_fig=True)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_nospread_0deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_nospread_0deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_nospread_0deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    tsk.pars["theta_obs"] = 0.785
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_nospread_45deg_res",
    #                     show_fig=False)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_nospread_45deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_nospread_45deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_nospread_45deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    tsk.pars["theta_obs"] = 1.57
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_nospread_90deg_res",
    #                     show_fig=False)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_nospread_90deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_nospread_90deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_nospread_90deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    # ----------------------------------- | WITH SPREADING | ----------------------------------

    grb = CasesFS(default_parfile_fpath=default_parfile_fpath, workingdir=workingdir)
    tsk.opts_a["method_spread"] = "AFGPY"
    tsk.opts_pw["method_spread"] = "AA"
    tsk.pars["theta_obs"] = 0.
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_spread_0deg_res",
    #                     show_fig=False,save_fig=True)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_spread_0deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_spread_0deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_spread_0deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    # grb = TestCasesFS(parfiledir)
    tsk.pars["theta_obs"] = 0.785
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_spread_45deg_res",
    #                     show_fig=False,save_fig=True)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_spread_45deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_spread_45deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_spread_45deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    # grb = TestCasesFS(parfiledir)
    tsk.pars["theta_obs"] = 1.57
    # grb.compare_skymaps(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
    #                     resolutions=tsk.resolutions_compare_plot,
    #                     figname=f"abstract_skymap_{tsk.figname}_spread_90deg_res",
    #                     show_fig=False,save_fig=True)
    grb.compare_skymap_evolution_plot_each_timestep(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.name_1,
                                                    resolutions_pw=tsk.resolutions_pw, resolutions_a=tsk.resolutions_a,
                                                    figfpath_skymap=figdir+f"abstract_skymap_{tsk.figname}_spread_90deg_res",
                                                    figfpath_evolve=figdir+f"abstract_skymap_{tsk.figname}_spread_90deg_propres",
                                                    fpath_out=workingdir+f"abstract_skymap_{tsk.figname}_spread_90deg_propres",
                                                    plots=("xc","xsize","ysize","fluxes"),
                                                    show_fig=False,save_pdf=True,save_res=False)

    # -------------------------------------------------------------------------------
    print(f"DONE for {tsk.figname}")