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

def example_run(sim, text : int):
    pr = PBA.parallel_runs.ParallelRuns(
        afg_data_dir=AFGRUNDIR,
        dir_for_run=sim["name"]+"/")

    iter_pars_dict = {
        "n_ism": [1.0, 0.1, 0.01, 0.001],
        "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
        "p": [2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
        "eps_e": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        "eps_b": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
        "eps_t": [1., 0.5, 0.1, 0.01]
    }
    pr.setup_1_working_dirs(iter_pars_dict=iter_pars_dict,dirname_prefix="kn_",dirname_ending=f"text{text}")

    pr.setup_2_parfiles(parfilename="parfile_def.par")

    id_pars = {"text":25, "method_r0":"from_beta", "t0":1e3, "new_theta_len":40, "new_vinf_len":None}
    pr.setup_3_id_knej(id_pars=id_pars, fpath_collated_ejecta=sim["datadir"] + "ej_collated.h5", show_id=False)

    skymap_postprocess_conf = {
        "nx":256, "ny":128, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
        "intp_filter":{ "type":None, "sigma":2, "mode":'reflect' }, # "gaussian"
        "hist_filter":{ "type":None, "sigma":2, "mode":'reflect' }
    }

    pr.launch_runs(n_cpu=1,
                   path_to_executable=EXECUTABLE,
                   skymap_postprocess_conf={},
                   loglevel="info")

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
    example_run( SIMULATIONS["SFHo_q125_res150"], text=25 )


if __name__ == '__main__':
    main()