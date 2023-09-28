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
from multiprocessing import Pool

from settings import DATADIR


try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

from settings import SIMULATIONS
from paths import *

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

    def setup_1_working_dirs(self, iter_pars_dict : dict, dirname_prefix : str = "kn_", dirname_ending : str = "afg") -> None:
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

    def setup_2_parfiles(self, parfilename : str = "parfile.par", fname_ejecta_id : str = "ejecta_id.h5") -> None:

        self.ejecta_id_fname = fname_ejecta_id

        if (not os.path.isfile(self.outdir + parfilename)):
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
            if len(pars.keys()) > 0:
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

    def setup_3_id(self, text : float = 25,
                 fpath_collated_ejecta : str = "/media/vsevolod/T7/work/KentaData/" + "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/" + "ej_collated.h5") -> None:

        if (not os.path.isfile(fpath_collated_ejecta)):
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
            path_to_executable=EXECUTABLE,
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

def example_run(sim, text : int):
    pr = ParallelRuns(
        afg_data_dir=AFGRUNDIR,
        dir_for_run=sim["name"]+"/")
    iter_pars_dict = {
        "n_ism": [1.0, 0.1, 0.01, 0.001]
        # "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
        # "p": [2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
        # "eps_e": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        # "eps_b": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
    }
    pr.setup_1_working_dirs(iter_pars_dict=iter_pars_dict,dirname_prefix="kn_",dirname_ending=f"text{text}")

    pr.setup_2_parfiles(parfilename="parfile_def.par")

    pr.setup_3_id(text=text, fpath_collated_ejecta=sim["datadir"] + "ej_collated.h5")

    skymap_postprocess_conf = {
        "nx":256, "ny":128, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
        "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
        "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }
    }
    pr.launch_runs(n_cpu=4, skymap_postprocess_conf=skymap_postprocess_conf, loglevel="err")

    print("RUNS FINISHED")

def main():
    times = [10., 20., 40., 60., 80., 100., 120., 140., 160., 180., 200., 220., 240., 260., 280.,
             300., 320., 360., 400.,
             460., 500., 560., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500.,
             1600., 1800., 2000., 2300.,
             2500., 2800., 3000., 3600., 4000., 4600., 5000., 6000., 7000., 8000., 9000., 10000.,
             12000., 15000.,
             20000.]
    # for t in times:
    #     print(t*PBA.utils.cgs.day,end=" ")
    # print("\n")
    example_run( SIMULATIONS[1], text=25 )


if __name__ == '__main__':
    main()