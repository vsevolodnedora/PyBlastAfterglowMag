import numpy as np
from sympy.polys.benchmarks.bench_solvers import k2

np.set_printoptions(precision=2)

import h5py
import pandas as pd
from glob import glob
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import json

from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import differential_evolution

import pysr
import sympy
from pysr import PySRRegressor
from sklearn.model_selection import train_test_split

# settings
#plt.style.use("fivethirtyeight")

# Supress runtime warning
import warnings
warnings.filterwarnings("ignore")

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

DATA_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/kenta_data/"

# load the metadata
with open(DATA_PATH+"metadata.json") as json_file:
    json_data = json.load(json_file)
    print(json_data)
SIMS = pd.DataFrame.from_dict(json_data).T
SIMS.set_index("name")
# select only new simulations
df = SIMS[SIMS["given_time"] == "new"]

get_ej_data = lambda name : DATA_PATH+name+'/'+"ej_collated.h5"

class ProcessRaw():
    def __init__(self, simumlation : dict):
        self.sim = simumlation
        self.dfile = None

    def process_raw_ejecta_files(self, infiles : str = "ejecta_*.h5", fname_output : str = "ej_collated.h5",
                                 mode:str="mass", overwrite:bool=False) -> None:
        # dir = "/media/vsevolod/T7/work/KentaData/"
        # simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
        name = self.sim["name"]
        datadir = self.sim["datadir"]
        collated_ej_data_fpath = datadir + fname_output
        # skip if file exists and no overwrite option set
        if (os.path.isfile(collated_ej_data_fpath)and(not overwrite)):
            return None
        if (not os.path.isdir(datadir)):
            raise FileNotFoundError(f"Datadir for {name} is not found Looking for {datadir}")
        files = glob(datadir + infiles)
        if (len(files) == 0):
            raise FileNotFoundError(f"Files {infiles} not found in {datadir}")
        id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True, mode=mode)
        id.process_save(collated_ej_data_fpath)

    def getProcessed(self, fname_output : str = "ej_collated.h5") -> h5py.File:
        if (self.dfile is None):
            self.dfile = h5py.File(self.sim["datadir"] + fname_output,"r")
        return self.dfile

    def getText(self) -> np.ndarray:
        return np.array(self.dfile["text"])

    def getData(self, v_n : str, text : int or None):
        if not ("time={}".format(text) in self.dfile.keys()):
            print(self.dfile.keys())
            raise KeyError("time={} is not found in dfile.keys() See above")
        if text == None:
            data = np.stack([self.dfile["time={}".format(time)][v_n] for time in self.getText()],axis=2)
        else:
            np.array(self.dfile["time={}".format(text)][v_n])

def process(sim : dict) -> None:
    ej = ProcessRaw(simumlation = sim)
    ej.process_raw_ejecta_files()
    print("Data collation is successful")


def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx

def task_process():
    pr = ProcessRaw(simumlation=SIMS.T["DD2_135_135_res150_floor"])
    pr.process_raw_ejecta_files(infiles= "ejecta_*.h5", fname_output= "ej_collated.h5", mode="mass",overwrite=False)
    #pr.process_raw_ejecta_files(infiles="Mdot_ejecta_*.h5", fname_output="mdot_ej_collated.h5", mode="mdot",overwrite=False)
    # Collate data for simulations (USE ONLY ONCE)\
    for sim_key, sim_dic in SIMS.items():
        pr = ProcessRaw(simumlation=sim_dic)
        #pr.process_raw_ejecta_files(infiles= "ejecta_*.h5", fname_output= "ej_collated.h5", mode="mass",overwrite=False)
        #pr.process_raw_ejecta_files(infiles="Mdot_ejecta_*.h5", fname_output="mdot_ej_collated.h5", mode="mdot",overwrite=False)

if __name__ == '__main__':
    task_process()