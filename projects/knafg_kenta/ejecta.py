"""
    Uses PyBlastAfterglowMag code
    From iside `\package\' run:
    pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .
    to install the postprocessing/run utilities

"""
import copy

import numpy as np
import h5py
import pandas as pd
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
from glob import glob
import os
import re

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")


from settings import *

class ProcessRaw():
    def __init__(self, simumlation : dict):
        self.sim = simumlation
        self.dfile = None

    def process_raw_ejecta_files(self, infiles : str = "ejecta_*.h5", fname_output : str = "ej_collated.h5", mode:str="mass"):
        # dir = "/media/vsevolod/T7/work/KentaData/"
        # simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
        name = self.sim["name"]
        datadir = self.sim["datadir"]
        if (not os.path.isdir(datadir)):
            raise FileNotFoundError(f"Datadir for {name} is not found Looking for {datadir}")
        files = glob(datadir + infiles)
        if (len(files) == 0):
            raise FileNotFoundError(f"Files {infiles} not found in {datadir}")
        id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True, mode=mode)
        collated_ej_data_fpath = datadir + fname_output
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




def main():

    for sim_key, sim_dic in SIMULATIONS.items():
        pr = ProcessRaw(simumlation=sim_dic)
        pr.process_raw_ejecta_files(infiles= "ejecta_*.h5", fname_output= "ej_collated.h5", mode="mass")

    exit(1);

    ej = ProcessRaw(simumlation = SIMULATIONS["BHBLp_q1_res150"])
    ej.process_raw_ejecta_files(infiles="Mdot_ejecta_*.h5",
                                fname_output="mdot_ej_collated.h5",
                                mode="mdot")
    print("Data collation is successful")
    exit(0)

    sim_dic = SIMULATIONS["BHBLp_q1_res150"]

    ej_data = PBA.id_kenta.EjectaData(sim_dic["datadir"]+"ej_collated.h5",verbose=True)



    data = PBA.id_kenta.Data(fpath_rhomax=sim_dic["datadir"]+sim_dic["rhomax"],
                             fpath_mdot=sim_dic["datadir"]+sim_dic["mdot_extract"])
    df = pd.merge(data.df_mdot,data.df_rho,on="time")
    print(df.shape)

    print(data.df_mdot.keys())
    # print(data.df_mdot.describe())
    print(data.df_mdot.isnull().sum(axis=0))

    rext = data.get_rext()
    fig, ax = plt.subplots(ncols=1,nrows=1)

    for i, r in enumerate(rext):
        try:
            ax.plot(data.ret_time(r0=r,vave_key="vave_fast"), data.df_mdot[f"mdot_fast_r{r}"].values)
        except:
            print(f"failed for r={r}")
        # df.plot("time",f"mdot_fast_r{r}",ax=ax,color=color_pal[i])
    ax2 = ax.twinx()
    data.df_rho.plot("time",f"rho_max",ax=ax2,color="black")

    ax.set_yscale("log")
    ax2.set_yscale("log")
    # ax.set_ylim(1e-5,1e1)
    plt.show()

    sim_dic = SIMULATIONS["BHBLp_q1_res150"]
    rhofpath = sim_dic["datadir"]+sim_dic["rhomax"]
    mdotfpath = sim_dic["datadir"]+sim_dic["mdot_extract"]
    ej_data = PBA.id_kenta.EjectaData(fpath=sim_dic["datadir"]+"ej_collated.h5",verbose=True)
    data = PBA.id_kenta.Data(fpath_rhomax=rhofpath,fpath_mdot=mdotfpath)
    print(data.df_mdot.keys())

    fig, ax = plt.subplots(ncols=1,nrows=1)
    ax.plot(ej_data.getText(), ej_data.total_mass(), color="gray",label="Total")
    ax.plot(ej_data.getText(), ej_data.total_mass_fasttail(),color="black",label=r"$\Gamma\beta>1$")
    ax2 = ax.twinx()
    data.df_rho.plot("time","rho_max",ax=ax2,color="green",label=r"$\rho_{\rm max}$")
    data.df_mdot.plot("time","mdot_r481",ax=ax2,color="green",label=r"$\dot{M}_{\rm ej}$")

    # ax2.plot(*data.get_rhomax(),color="green",label=r"$\rho_{\rm max}$")
    # ax2.plot(*data.get_mdot2(),color="green",label=r"$\dot{M}_{\rm ej}$")
    ax2.set_yscale("log")
    plt.legend()
    ax.set_yscale("log")
    plt.show()

    # for sim in SIMULATIONS:
        # process_data( sim )

    df = pd.DataFrame(SIMULATIONS).T
    # df["n"] = {"SFHo_q1_res150":1,"SFHo_q111_res150":1,"SFHo_q116_res150":1,"SFHo_q116_res200":1,"SFHo_q125_res150":1}

    df["tmin"],df["tmax"] = {}, {}
    for sim, sim_dic in df.iterrows():
        data = ProcessRaw(sim_dic)
        dfile = data.getProcessed()
        # print(dfile.keys())
        # print(dfile.attrs.keys())
        df["tmin"][sim] = np.array(dfile["text"]).min()
        df["tmax"][sim] = np.array(dfile["text"]).max()

    print(df[["tmin","tmax"]])

if __name__ == '__main__':
    main()

def save_data_for_gilad():
    workdir = os.getcwd()+'/'
    # path_to_original_data = "/media/vsevolod/data/KentaData/SFHo_13_14_150m_11/" #
    # dfile = h5py.File(workdir+"kenta_ejecta_13.h5")
    # print(dfile.keys())

    text = 70.
    label = f"corr_id_SFHo_13_14_150m_11_text{int(text)}"
    sort_by = lambda k: int(re.findall(r'\d+', str(k.split("/")[-1]))[0])
    files = sorted(glob(DATADIR + "ejecta*.h5", recursive=True), key=sort_by)
    if (len(files) == 0):
        raise FileNotFoundError(f"not found in: {DATADIR}")

    PBA.id_maker_from_kenta_bns.prepare_kn_ej_id_2d(files=files,
                                                    outfpaths=[workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5"],
                                                    req_times=np.array([text]),
                                                    new_theta_len=None,
                                                    new_vinf_len=None,
                                                    verbose=True,
                                                    r0type="frombeta",#"fromrho",
                                                    r0frac=0.5,
                                                    t0=-1,
                                                    dist="pw")

    r_, mom_, theta_, ctheta_, ek_, mass_, ye_, rho_, temp_, press_, eps_, entr_ \
        = PBA.id_maker_from_kenta_bns.load_init_data(workdir+f"corr_id_SFHo_13_14_150m_11_text{int(text)}.h5")
    print(repr(eps_))
    # PBA.id_maker_from_kenta_bns.plot_init_profile(ctheta_[0,:], mom_[:,0], press_,
    #                   figpath=None,#FIGPATH+"ang_mass_dist" + r"_text{}".format(int(times[idx])),
    #                   norm_mode="log",
    #                   subplot_mode="ave",
    #                   title=r"$t_{\rm ext}="+r"{}$ [ms]".format(int(text)))

    idx=0
    np.savetxt(os.getcwd()+f"/ejecta_layer{idx}.dat",X=np.column_stack((mom_[:,idx],mass_[:,idx],rho_[:,idx],eps_[:,idx],press_[:,idx])),
               header="GammaBeta Mass rho eps pres")
    idx = len(mom_[0,:])-1
    np.savetxt(os.getcwd()+f"/ejecta_layer{idx}.dat",X=np.column_stack((mom_[:,idx],mass_[:,idx],rho_[:,idx],eps_[:,idx],press_[:,idx])),
               header="GammaBeta Mass rho eps pres")



