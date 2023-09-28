"""
    Uses PyBlastAfterglowMag code
    From iside `\package\' run:
    pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .
    to install the postprocessing/run utilities

"""
import copy

import numpy as np
import h5py
import shutil
from glob import glob
from multiprocessing import Pool
import plotly.graph_objects as go
import plotly.io as pio
# png_renderer = pio.renderers["png"]
# png_renderer.width = 500
# png_renderer.height = 500
#
# pio.renderers.default = "png"

from scipy.integrate import ode
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')

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
from matplotlib import cm
import os
import re

# import PyBlastAfterglowMag as PBA
import package.src.PyBlastAfterglowMag as PBA
from settings import *

def parallel_run_parameter_exploration(iter_pars_dict : dict, prefix : str, ending : str,
                                       root_workingdir : str, source_parfiledir : str, parfilename : str):

    # check if the output dir does exist
    if not (os.path.isdir(root_workingdir)):
        raise FileNotFoundError(f"root_workingdir is not found: {root_workingdir}")
    if not (os.path.isdir(source_parfiledir)):
        raise FileNotFoundError(f"parfiledir is not found: {source_parfiledir}")
    if not (os.path.isfile(source_parfiledir + parfilename)):
        raise FileNotFoundError(f"parfile in parfiledir is not found: {source_parfiledir + parfilename}")

    # generate directory names for the runs (concacenate parameters to iterate over)
    workdir_template = "".join([f"{key}[{key}]_" for key in iter_pars_dict.keys()])

    # add user-specified prefix if needed
    workdir_template = prefix + workdir_template

    # create lists of parameters by doing permutations on the given 'iter_pars_dict'
    new_pars = PBA.parfile_tools.set_parlists_for_pars(
        iter_pars_keys=list(iter_pars_dict.keys()), iter_pars=iter_pars_dict, fname=workdir_template)

    # duplicate the list for 'main' 'kn' and 'grb' parameters
    new_pars_main = copy.deepcopy(new_pars)
    new_pars_grb = copy.deepcopy(new_pars)
    new_pars_kn = copy.deepcopy(new_pars)

    # main loop:
    working_dirs = []
    for i in range(len(new_pars_main)):

        # create directory for each run
        workingdir = root_workingdir + new_pars_main["name"] + ending + "/"
        working_dirs.append(workingdir)
        if not (os.path.isdir(workingdir)):
            os.mkdir(workingdir)
        else:
            print(f"directory already existis {workingdir}")

        # check if parfile is already in the working dir
        if not (os.path.isfile(workingdir+parfilename)):
            print(f"Copyting {source_parfiledir + parfilename} -> {workingdir + parfilename}")
            shutil.copyfile(source_parfiledir + parfilename, workingdir + parfilename)
        else:
            print(f"parfile already exists {workingdir+parfilename}")

        # modify workingdir parfile in accordance with the current iteration
        PBA.parfile_tools.modify_parfile_par_opt(part="main", newpars=new_pars[i], newopts={},
                                                 workingdir=workingdir, parfile=parfilename,
                                                 newparfile="parfile.par", keep_old=True)
        PBA.parfile_tools.modify_parfile_par_opt(part="kn", newpars=new_pars[i], newopts={},
                                                 workingdir=workingdir, parfile="parfile.par",
                                                 newparfile="parfile.par", keep_old=False)
    return (new_pars, working_dirs)

class Ejecta():
    def __init__(self, properties : dict):
        self.sim = properties

    def process_raw_ejecta_files(self, infiles : str = "ejecta_*.h5", fname_output : str = "ej_collated.h5"):
        # dir = "/media/vsevolod/T7/work/KentaData/"
        # simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
        name = self.sim["name"]
        datadir = self.sim["datadir"]
        if (not os.path.isdir(datadir)):
            raise FileNotFoundError(f"Datadir for {name} is not found Looking for {datadir}")
        files = glob(datadir + infiles)
        if (len(files) == 0):
            raise FileNotFoundError("Files not found")
        id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True)
        collated_ej_data_fpath = datadir + fname_output
        id.process_save(collated_ej_data_fpath)


class ParallelRuns():
    def __init__(self,
                 afg_data_dir : str = "/media/vsevolod/T7/work/afterglowKentaProject/",
                 dir_for_run : str= "run1/") :
        # check if output dir where to store simulation results exists
        if (not os.path.isdir(afg_data_dir)):
            raise FileNotFoundError(f"directory for simulation output not found: {afg_data_dir}")

        # create the directory for this set of runs
        self.outdir = afg_data_dir + dir_for_run
        if (not os.path.isdir(self.outdir)):
            os.mkdir(self.outdir)
        else:
            print("Directory for this set of runs already exists: {}".format(self.outdir))

        # copy the parfile from HERE to mother dir for runs
        # shutil.copyfile(os.getcwd() + "/" + "parfile.par", outdir + "parfile.par")
        self.ejecta_id_fname = None
        self.working_dirs = []
        self.new_pars = []

    def setup_working_dirs(self, iter_pars_dict : dict, dirname_prefix : str = "kn_", dirname_ending : str = "afg") -> None:
        """
        iter_pars_dict = {
            "n_ism": [1.0, 0.1, 0.01, 0.001],
            "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
            "p": [2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
            "eps_e": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
            "eps_b": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
        }
        :return:
        """
        if (iter_pars_dict.keys() == 0):
            print("Empty iter_pars_dict is given")

        # generate directory names for the runs (concacenate parameters to iterate over)
        workdir_template = "".join(["{}[{}]_".format(key.replace("_",""),key) for key in iter_pars_dict.keys()])
        # add user-specified prefix if needed
        workdir_template = dirname_prefix + workdir_template
        # create lists of parameters by doing permutations on the given 'iter_pars_dict'
        self.new_pars = PBA.parfile_tools.set_parlists_for_pars(
            iter_pars_keys=list(iter_pars_dict.keys()), iter_pars=iter_pars_dict, fname=workdir_template)

        # create directories for each run and copy there the parfile from outdir
        ending = "afg"
        parfilename = "parfile.par"
        if len(self.new_pars) and len(self.new_pars[0].keys()) > 0:
            self.working_dirs = [self.outdir + new_par["name"] + dirname_ending + "/" for new_par in self.new_pars]
        else:
            self.working_dirs = [self.outdir + dirname_ending + "/"]
        print(f"NOTE {len(self.working_dirs)} directories will be created!")
        print(f"\t Example {self.working_dirs[0]}")

    def setup_parfiles(self, parfilename : str = "parfile.par", fname_ejecta_id : str = "ejecta_id.h5") -> None:

        self.ejecta_id_fname = fname_ejecta_id

        if (~os.path.isfile(self.outdir + parfilename)):
            raise FileNotFoundError(f"Source parfile not found {self.outdir + parfilename}")

        # Change Parfile in each directory with new parameters
        for workingdir, pars in zip(self.working_dirs, self.new_pars):

            # check if dirs for runs exist
            if not (os.path.isdir(workingdir)):
                print(f"\t Creating... {workingdir}")
                os.mkdir(workingdir)
            else:
                print(f"\tdirectory already existis {workingdir}")

            # check if parfile is already in the working dir
            if not (os.path.isfile(workingdir+parfilename)):
                print(f"\tCopyting {self.outdir + parfilename} -> {workingdir + parfilename}")
                shutil.copyfile(self.outdir + parfilename, workingdir + parfilename)
            else:
                print(f"\tparfile already exists {workingdir+parfilename}")

            # modify workingdirs parfile in accordance with the current iteration
            if len(pars.keys() > 0):
                print(f"\tModifying {workingdir + parfilename}")
                PBA.parfile_tools.modify_parfile_par_opt(part="main", newpars=pars,
                                                         newopts={},
                                                         workingdir=workingdir, parfile=parfilename,
                                                         newparfile="parfile.par", keep_old=True)
                PBA.parfile_tools.modify_parfile_par_opt(part="kn", newpars=pars,
                                                         newopts={"fname_ejecta_id": fname_ejecta_id},
                                                         workingdir=workingdir, parfile="parfile.par",
                                                         newparfile="parfile.par", keep_old=False)
            else:
                print(f"\tNo pars given, not parfiles are modified")

    def setup_id(self, text : float = 25,
                 fpath_collated_ejecta : str = "/media/vsevolod/T7/work/KentaData/" + "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/" + "ej_collated.h5") -> None:

        if (~os.path.isfile(fpath_collated_ejecta)):
            raise FileNotFoundError(f"Collated ejecta file not found: {fpath_collated_ejecta}")
        if (len(self.working_dirs)==0):
            raise RuntimeError(f"No directories exist for running: {len(self.working_dirs)}")

        id = PBA.id_kenta.EjStruct(fpath=fpath_collated_ejecta, verbose=True)
        id_dict = id.get_2D_id(text=text, method_r0="from_beta", t0=1e3, new_theta_len=None, new_vinf_len=None)
        id.plot_init_profile(mom=id_dict["mom"][:, 0], ctheta=id_dict["ctheta"][0, :], mass=id_dict["ek"])
        for working_dir_i in self.working_dirs:
            self.ejecta_id_fpath = working_dir_i + self.ejecta_id_fname
            with h5py.File(self.ejecta_id_fpath, "w") as dfile:
                dfile.attrs.create("text",data=text)
                for key, data in id_dict.items():
                    dfile.create_dataset(name=key, data=np.copy(data))

    def launch_runs(self, n_cpu : int, skymap_postprocess_conf : dict, loglevel : str = "info") -> None:
        """

        skymap_postprocess_conf = {
            "nx":48, "ny":32, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
            "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
            "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }
        }

        :param n_cpu:
        :param skymap_postprocess_conf:
        :return:
        """

        pba_parallel = PBA.parallel_runs.ParallelRunDispatcher(
            path_to_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
            loglevel=loglevel,
            working_dirs=self.working_dirs,
            parfile_name="parfile.par",
            skymap_postprocess_conf=skymap_postprocess_conf)
        if (n_cpu == 1):
            for i, pars in enumerate(self.working_dirs):
                pba_parallel(i)
        else:
            if (n_cpu is None):
                ncpus = os.cpu_count()
            else:
                ncpus = int(n_cpu)
            try:
                pool = Pool(ncpus)  # on 8 processors
                pool.map(pba_parallel, range(len(self.working_dirs)))  # distribute tasks for each clkass
                # pool.map(pb, *[(pars, kn_pars, grb_pars) for pars, kn_pars, grb_pars in zip(pars_list, grb_pars_list, kn_pars_list)])  # distribute tasks for each clkass
            finally:  # To make sure processes are closed in the end, even if errors happen
                pool.close()
                pool.join()

def run():

    nr_data_dir = "/media/vsevolod/T7/work/KentaData/"
    nr_sim_name = "SFHoTim276_13_14_0025_150mstg_B0_HLLC"

    afg_data_dir = "/media/vsevolod/T7/work/afterglowKentaProject/"

    dir_for_run = "run1/"


    # check if inout dir with data exists
    if (not os.path.isdir(nr_data_dir)):
        raise FileNotFoundError(f"nr_data_dir not found: {nr_data_dir}")
    # check if inout dir with data exists
    if (not os.path.isdir(nr_data_dir+nr_sim_name)):
        raise FileNotFoundError(f"nr_data_dir for simulation data not found: {nr_data_dir+nr_sim_name}")
    # check if output dir where to store simulation results exists
    if (not os.path.isdir(afg_data_dir)):
        raise FileNotFoundError(f"directory for simulation output not found: {afg_data_dir}")

    # create the directory for this set of runs
    outdir = afg_data_dir + dir_for_run
    if (not os.path.isdir(outdir)):
        os.mkdir(outdir)
    else:
        print("Directory for this set of runs already exists: {}".format(outdir))

    # copy the parfile from HERE to mother dir for runs
    shutil.copyfile(os.getcwd() + "/" + "parfile.par", outdir + "parfile.par")

    iter_pars_dict = {
        "n_ism": [1.0, 0.1, 0.01, 0.001],
        "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
        "p": [2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
        "eps_e": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        "eps_b": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
    }
    iter_pars_dict = {}

    # generate directory names for the runs (concacenate parameters to iterate over)
    workdir_template = "".join(["{}[{}]_".format(key.replace("_",""),key) for key in iter_pars_dict.keys()])
    # add user-specified prefix if needed
    prefix = "kn_"
    workdir_template = prefix + workdir_template
    # create lists of parameters by doing permutations on the given 'iter_pars_dict'
    new_pars = PBA.parfile_tools.set_parlists_for_pars(
        iter_pars_keys=list(iter_pars_dict.keys()), iter_pars=iter_pars_dict, fname=workdir_template)

    # create directories for each run and copy there the parfile from outdir
    ending = "afg"
    parfilename = "parfile.par"
    if len(new_pars) and len(new_pars[0].keys()) > 0:
        working_dirs = [outdir + new_par["name"] + ending + "/" for new_par in new_pars]
    else:
        working_dirs = [outdir + ending + "/"]
    print(f"NOTE {len(working_dirs)} directories will be created!")
    print(f"\t Example {working_dirs[0]}")

    # Change Parfile in each directory with new parameters
    for workingdir, pars in zip(working_dirs, new_pars):
        if not (os.path.isdir(workingdir)):
            print(f"\t Creating... {workingdir}")
            os.mkdir(workingdir)
        else:
            print(f"\tdirectory already existis {workingdir}")
        # check if parfile is already in the working dir
        if not (os.path.isfile(workingdir+parfilename)):
            print(f"\tCopyting {outdir + parfilename} -> {workingdir + parfilename}")
            shutil.copyfile(outdir + parfilename, workingdir + parfilename)
        else:
            print(f"\tparfile already exists {workingdir+parfilename}")

        # modify workingdir parfile in accordance with the current iteration
        if len(pars.keys()):
            print(f"\tModifying {workingdir + parfilename}")
            PBA.parfile_tools.modify_parfile_par_opt(part="main", newpars=pars, newopts={},
                                                     workingdir=workingdir, parfile=parfilename,
                                                     newparfile="parfile.par", keep_old=True)
            PBA.parfile_tools.modify_parfile_par_opt(part="kn", newpars=pars, newopts={},
                                                     workingdir=workingdir, parfile="parfile.par",
                                                     newparfile="parfile.par", keep_old=False)


    # process and collate initial data (NEEDED only ONCE)
    dir = "/media/vsevolod/T7/work/KentaData/"
    simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
    files = glob(dir + simname + "ejecta_*.h5")
    if (len(files) == 0):
        raise FileNotFoundError("Files not found")
    id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True)
    collated_ej_data_fpath = dir + simname + "ej_collated.h5"
    id.process_save(collated_ej_data_fpath)

    # extract data from collated file for a given extraction time and save it into working dir.
    text = 25
    id = PBA.id_kenta.EjStruct(fpath=collated_ej_data_fpath, verbose=True)
    id_dict = id.get_2D_id(text=text, method_r0="from_beta", t0=1e3, new_theta_len=None, new_vinf_len=None)
    id.plot_init_profile(mom=id_dict["mom"][:, 0], ctheta=id_dict["ctheta"][0, :], mass=id_dict["ek"])
    for working_dir_i in working_dirs:
        outfnmae = working_dir_i + "ejecta_id.h5"
        with h5py.File(outfnmae, "w") as dfile:
            dfile.attrs.create("text",data=text)
            for key, data in id_dict.items():
                dfile.create_dataset(name=key, data=np.copy(data))

    # save data
    # mom = id_dict["mom"][:,0].flatten()
    # # ctheta = id_dict["ctheta"].flatten()
    # ek = np.sum(id_dict["ek"],axis=1).flatten()
    # mom = mom[ek>0 & np.isfinite(ek)]
    # # ctheta = ctheta[ek>0 & np.isfinite(ek)]
    # ek = ek[ek>0 & np.isfinite(ek)]
    # out = np.column_stack((mom, ek))
    # np.savetxt(os.getcwd()+'/'+"data.csv", out, delimiter=",")
    # exit(1)


    # run PyBlastAfterglow in parallel in several directories at the same time
    n_cpu = 1
    skymap_postprocess_conf = {
        "nx":48, "ny":32, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
        "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
        "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }}
    pba_parallel = PBA.parallel_runs.ParallelRunDispatcher(
        path_to_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out", loglevel="info",
        working_dirs=working_dirs, parfile_name="parfile.par", skymap_postprocess_conf=skymap_postprocess_conf)
    if (n_cpu == 1):
        for i, pars in enumerate(working_dirs):
            pba_parallel(i)
    else:
        if (n_cpu is None):
            ncpus = os.cpu_count()
        else:
            ncpus = int(n_cpu)
        try:
            pool = Pool(ncpus)  # on 8 processors
            pool.map(pba_parallel, range(len(working_dirs)))  # distribute tasks for each clkass
            # pool.map(pb, *[(pars, kn_pars, grb_pars) for pars, kn_pars, grb_pars in zip(pars_list, grb_pars_list, kn_pars_list)])  # distribute tasks for each clkass
        finally:  # To make sure processes are closed in the end, even if errors happen
            pool.close()
            pool.join()


def plot_dyn(dfile, ishells=(1,10,30), ilayers=(0, 10, 22), v_n_x="R", v_n_ys=("rho", "mom"),
             colors_by="layers", legend=False):
    layers = []
    for i in ishells:
        for j in ilayers:
            layers.append("shell={} layer={}".format(i, j))
    # print(f"nlayers={dfile.attrs.keys()}")
    # print(f"nlayers={dfile.keys()}")
    fid, axes = plt.subplots(ncols=1, nrows=len(v_n_ys), figsize=(6, 6), sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis
    mynorm = Normalize(vmin=0, vmax=len(ishells) * len(ilayers))  # norm(len(ishells)*len(ilayers))

    print(dfile.keys())

    for iv_n, v_n in enumerate(v_n_ys):
        i = 0
        ax = axes[iv_n] if len(v_n_ys) > 1 else axes
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            y_arr = np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]

            if (colors_by == "layers"):
                color = cmap(mynorm(int(i)))  # color=cmap(norm(int(layer.split("layer=")[-1])))
            else:
                color = cmap(mynorm(int(i)))  # color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))

            if (v_n_x == "tburst"): x_arr /= PBA.utils.cgs.day;
            ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
            i = i + 1
        ax.set_xlabel(v_n_x)
        if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]")
        ax.set_ylabel(v_n)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.grid(visible=True)
        ax.minorticks_on()
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(5e-1,4e3)
    if legend:
        plt.legend()
    plt.show()
def plot():
    afg_data_dir = "/media/vsevolod/T7/work/afterglowKentaProject/"
    dir_for_run = "run1/"
    outdir = afg_data_dir + dir_for_run
    ending = "afg"
    parfilename = "parfile.par"

    pba = PBA.interface.PyBlastAfterglow(workingdir=outdir+ending+'/',readparfileforpaths=True,parfile=parfilename)

    plot_dyn(pba.KN.get_dyn_obj(),ishells=(0,10,20,30),ilayers=(1,),v_n_x="R",v_n_ys=("mom","Eint2"))


    times = pba.KN.get_lc_times(unique=True,spec=False)
    freqs = pba.KN.get_lc_freqs(unique=True,spec=False)
    fluxes = pba.KN.get_lc_totalflux(freq=freqs[0],time=None,spec=False)

    fig, ax = plt.subplots(ncols=1,nrows=1)
    ax.loglog(times,fluxes)
    ax.grid(visible=True)
    plt.title('ax.grid(True)', family='monospace', fontsize='small')
    plt.show()


if __name__ == '__main__':
    # run()
    plot()