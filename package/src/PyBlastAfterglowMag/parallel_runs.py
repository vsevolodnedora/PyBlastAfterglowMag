import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy
import subprocess
from shutil import copyfile
import shutil
from multiprocessing import Pool

from .utils import cgs, find_nearest_index

from .interface import PyBlastAfterglow
from .skymap_process import ProcessRawSkymap
from .parfile_tools import modify_parfile_par_opt, set_parlists_for_pars
from .id_analytic import JetStruct
from .id_kenta import EjStruct

class ParallelRunDispatcher:
    def __init__(self,path_to_executable:str,parfile_name:str,loglevel:str,
                 working_dirs : list[str], skymap_postprocess_conf : dict):
        assert len(working_dirs) > 0, "no simulation direrctories given"
        for sim_dir in working_dirs:
            if not os.path.isdir(sim_dir):
                raise FileNotFoundError(f"Simulation dir not found: {sim_dir}")
            if not os.path.isfile(sim_dir+parfile_name):
                raise FileNotFoundError(f"parfile not found: {sim_dir+parfile_name}")
        if not os.path.isfile(path_to_executable):
            raise FileNotFoundError(f"path_to_executable is invalid: {path_to_executable}")
        self.sim_dirs = working_dirs
        self.parfile_name = parfile_name
        self.skymap_postprocess_conf = skymap_postprocess_conf
        self.path_to_executable = path_to_executable
        self.loglevel = loglevel
        self.skymap_postprocess = len(skymap_postprocess_conf.keys()) > 0

    def __call__(self, idx):
        ''' run PyBlastAfterglow with parfile in a given working directory '''
        if (idx > (len(self.sim_dirs) - 1)):
            raise ValueError("simdir given {} while index requiested {}".format(len(self.sim_dirs), idx))
        sim_dir = self.sim_dirs[idx]
        if not (os.path.isfile(sim_dir+self.parfile_name)):
            raise FileNotFoundError("parfile not found {}".format(sim_dir+self.parfile_name))

        pba = PyBlastAfterglow(workingdir=sim_dir, readparfileforpaths=True, parfile=self.parfile_name)

        pba.run(path_to_cpp_executable=self.path_to_executable,loglevel=self.loglevel)

        if (self.skymap_postprocess):
            outfpath = None
            if ((not pba.KN.fpath_sky_map is None) and (pba.GRB.fpath_sky_map is None)):
                outfpath = pba.KN.fpath_sky_map
            elif ((pba.KN.fpath_sky_map is None) and (not pba.GRB.fpath_sky_map is None)):
                outfpath = pba.GRB.fpath_sky_map
            else:
                raise KeyError("Expected either pba.KN.fpath_sky_map != None or pba.KN.fpath_sky_map != None")
            prep = ProcessRawSkymap(conf=self.skymap_postprocess_conf, verbose=False)
            prep.process_singles(infpaths=sim_dir+"raw_skymap_*.h5", outfpath=outfpath, remove_input=True)

        pba.clear()

class ParallelRuns():
    def __init__(self,
                 afg_data_dir : str = "/media/vsevolod/T7/work/afterglowKentaProject/", # where to find main parfile
                 dir_for_run : str= "run1/" # where to setup all the runs (one per its own directory
                 ):
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
        # create lists of parameters by doing permutations on the given 'iter_pars_dict'
        self.new_pars = set_parlists_for_pars(
            iter_pars_keys=list(iter_pars_dict.keys()),
            iter_pars=iter_pars_dict)

        # create directories for each run and copy there the parfile from outdir
        if len(self.new_pars) and len(self.new_pars[0].keys()) > 0:
            self.working_dirs = [self.outdir + dirname_prefix + '_' + new_par["name"] + dirname_ending + "/"
                                 for new_par in self.new_pars]
        else:
            self.working_dirs = [self.outdir + dirname_prefix + dirname_ending + "/"]
        print(f"NOTE {len(self.working_dirs)} directories will be created!")
        print(f"\t Example {self.working_dirs[0]}")

    def setup_2_parfiles(self, parfilename : str = "parfile.par", fname_ejecta_id : str = "ejecta_id.h5", type="kn") -> None:

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
                modify_parfile_par_opt(
                    part="main", newpars=pars, newopts={}, workingdir=workingdir,
                    parfile=parfilename,newparfile="parfile.par", keep_old=True)
                modify_parfile_par_opt(
                    part=type, newpars=pars, newopts={"fname_ejecta_id": fname_ejecta_id},
                    workingdir=workingdir, parfile="parfile.par", newparfile="parfile.par", keep_old=False)
            else:
                print(f"\tNo pars given, not parfiles are modified")

    def setup_3_id_knej(self, id_pars : dict, fpath_collated_ejecta : str, show_id : bool) -> None:
        """

        :param id_pars: {"text":25, "method_r0":"from_beta", "t0":1e3, "new_theta_len":None, "new_vinf_len":None}
        :param fpath_collated_ejecta: "/media/vsevolod/T7/work/KentaData/"
                                                      + "SFHoTim276_12_15_0025_150mstg_B0_HLLC"
                                                      + "/" + "ej_collated.h5"
        :param show_id: False
        :return:
        """
        if (not os.path.isfile(fpath_collated_ejecta)):
            raise FileNotFoundError(f"Collated ejecta file not found: {fpath_collated_ejecta}")
        if (len(self.working_dirs)==0):
            raise RuntimeError(f"No directories exist for running: {len(self.working_dirs)}")

        id = EjStruct(fpath=fpath_collated_ejecta, verbose=True)
        id_dict = id.get_2D_id(**id_pars)
        if show_id:
            id.plot_init_profile(mom=id_dict["mom"][:, 0], ctheta=id_dict["ctheta"][0, :], mass=id_dict["ek"])
        for working_dir_i in self.working_dirs:
            self.ejecta_id_fpath = working_dir_i + self.ejecta_id_fname
            with h5py.File(self.ejecta_id_fpath, "w") as dfile:
                dfile.attrs.create("text",data=id_pars["text"])
                for key, data in id_dict.items():
                    dfile.create_dataset(name=key, data=np.copy(data))

    def setup_3_id_grbej(self, type_eats : str = "adaptive", overwrite : bool = True) -> None:

        if (len(self.working_dirs)==0):
            raise RuntimeError(f"No directories exist for running: {len(self.working_dirs)}")

        # Change Parfile in each directory with new parameters
        for workingdir, pars in zip(self.working_dirs, self.new_pars):
            n_layers_pw = pars["nlayers_pw"]
            nlayers_a = pars["nlayers_a"]
            pba_id = JetStruct(n_layers_pw=n_layers_pw, n_layers_a=nlayers_a)
            id_dict, id_pars = pba_id.get_1D_id(pars=pars,type=type_eats)
            ejecta_id_fpath = workingdir + self.ejecta_id_fname
            if (os.path.isfile(ejecta_id_fpath) and not overwrite):
                print(f"Skipping creating ID (file already exists) {ejecta_id_fpath}")
            else:
                print(f"Creating ID {ejecta_id_fpath}")
                pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=ejecta_id_fpath)


    def launch_runs(self, n_cpu : int,
                    skymap_postprocess_conf : dict,
                    path_to_executable : str,
                    loglevel : str = "info") -> None:
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

        pba_parallel = ParallelRunDispatcher(
            path_to_executable=path_to_executable,
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

def distribute_and_parallel_run(path_to_executable:str, working_dirs:list[str], parfile_name:str, n_cpu:int, pba_loglevel:int):
    ''' run multiple instances of PyBlastAfterglow in different directories '''
    pba_parallel = ParallelRunDispatcher(
        path_to_executable=path_to_executable,working_dirs=working_dirs,
        parfile_name=parfile_name,loglevel=pba_loglevel,
        skymap_postprocess_conf=None
    )
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
    print(f"Finished 'distribute_and_parallel_run()' for N={len(working_dirs)} working directories")

