import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy
import subprocess
from shutil import copyfile
from multiprocessing import Pool

from .utils import cgs, get_beta, find_nearest_index

from .interface import PyBlastAfterglow
from .skymap_process import ProcessRawSkymap

class ParallelRunDispatcher:
    def __init__(self, working_dirs : list[str], parfile_name:str, skymap_postprocess_conf : dict):
        assert len(working_dirs) > 0, "no simulation direrctories given"
        for sim_dir in working_dirs:
            if not os.path.isdir(sim_dir):
                raise FileNotFoundError(f"Simulation dir not found: {sim_dir}")
            if not os.path.isfile(sim_dir+parfile_name):
                raise FileNotFoundError(f"parfile not found: {sim_dir+parfile_name}")
        self.sim_dirs = working_dirs
        self.parfile_name = parfile_name
        self.skymap_postprocess_conf = skymap_postprocess_conf

    def __call__(self, idx):
        ''' run PyBlastAfterglow with parfile in a given working directory '''
        if (idx > (len(self.sim_dirs) - 1)):
            raise ValueError("simdir given {} while index requiested {}".format(len(self.sim_dirs), idx))
        sim_dir = self.sim_dirs[idx]
        if not (os.path.isfile(sim_dir+self.parfile_name)):
            raise FileNotFoundError("parfile not found {}".format(sim_dir+self.parfile_name))

        pba = PyBlastAfterglow(workingdir=sim_dir, readparfileforpaths=True, parfile=self.parfile_name)

        pba.run()

        if (pba.KN.opts["do_skymap"]=="yes"):
            prep = ProcessRawSkymap(conf=self.skymap_postprocess_conf, verbose=False)
            prep.process_singles(infpaths=sim_dir+"raw_skymap_*.h5",
                                 outfpath=pba.KN.fpath_sky_map, remove_input=True)

        pba.clear()

def distribute_and_parallel_run(working_dirs:list[str], parfile_name:str, n_cpu:int):
    ''' run multiple instances of PyBlastAfterglow in different directories '''
    pba_parallel = ParallelRunDispatcher(working_dirs=working_dirs, parfile_name=parfile_name)
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

