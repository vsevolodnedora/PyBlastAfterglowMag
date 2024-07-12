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

def get_text_max_mass(ej:PBA.id_kenta.EjectaData, crit:str='fast')->float:
    return float( ej.getText()[np.argmax(ej.total_mass_vs_text(crit=crit))] )
def get_text_max_mom(ej:PBA.id_kenta.EjectaData, crit:str='fast')->float:
    vals = []
    for it, text in enumerate(ej.getText()):
        ej_mass = ej.get(v_n="mass",text=text)
        mask = ej.get_vinf_mask(crit=crit)
        ej_vinf = ej.get_vinf()

        mass_ = np.sum(ej_mass[:,mask],axis=0)
        mom = PBA.MomFromBeta( ej_vinf[mask] ) #e j_vinf[mask]*PBA.utils.get_Gamma(ej_vinf[mask])
        idx_ = np.argmax(mom[mass_ > 0.]) if len(mom[mass_ > 0.]) > 0 else 0

        vals.append(mom[idx_])
    return float ( ej.getText()[np.argmax(vals)] )

# compute mass-averaged quantites for EACH text
def get_vave_theta_rms(sim_dic:dict, crit:str) -> pd.DataFrame:
    vals = {"vave":[], "mom_max":[], "ye_ave":[], "theta_rms":[],
            "entr_ave":[],"eps_ave":[],"temp_ave":[],
            "time":[], "mass":[], "ek":[]}

    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)

    for it, text in enumerate(ej.getText()):

        ej_mass = ej.get(v_n="mass",text=text)

        mask = ej.get_vinf_mask(crit=crit)

        vals["mass"].append( np.sum(ej_mass[:,mask]) )

        ej_vinf = ej.get_vinf()
        vals["vave"].append( np.sum(ej_mass[:,mask]*ej_vinf[mask])/np.sum(ej_mass[:,mask]) )

        mass_ = np.sum(ej_mass[:,mask],axis=0)
        mom = ej_vinf[mask]*PBA.utils.get_Gamma(ej_vinf[mask])
        idx_ = np.argmax(mom[mass_ > 0.]) if len(mom[mass_ > 0.]) > 0 else 0
        vals["mom_max"].append(mom[idx_])

        solar_m = 1.989e+33
        c = 2.9979e10
        vals["ek"].append( np.sum( ej_mass[:,mask] * solar_m * (ej_vinf[np.newaxis, mask] * c) ** 2 ))


        for v_n in ["ye","entr","eps","temp"]:
            ave = ej.get(v_n=v_n,text=text)
            vals[f"{v_n}_ave"].append( np.sum(ej_mass[:,mask]*ave[:,mask])/np.sum(ej_mass[:,mask]) )

        thetas = ej.get_theta()
        vals["theta_rms"].append( (180. / np.pi) * np.sqrt(np.sum(np.sum(ej_mass[:,mask], axis=1) * thetas ** 2) / np.sum(ej_mass[:,mask])) )
    # vals = {key:np.array(vals) for (key,val) in vals.items()}
    # print(vals["mom_max"])

    tmerg = sim_dic["tmerg"]
    times = ej.getText()
    vals["time"] = times - tmerg

    return pd.DataFrame.from_dict(vals)

def plot_all_sim_ejecta_massave_vel(crit=None, yscale="linear", figsize=(5,10), xlim:tuple=(-2,40),ylim:tuple=(0,1e-2),
                                    do_vave=False,do_mom_max=False, do_theta=False, do_ye=False,
                                    figname:str= "mej_vej_theta.png"):
    nrows = sum([1,do_ye,do_vave,do_mom_max,do_theta])
    fig, ax = plt.subplots(ncols=1,nrows=nrows,sharex="all",figsize=figsize,layout='constrained')
    if (nrows == 1):
        ax = [ax]
    #ax2 = ax.twinx()
    for idx, sim_dic in enumerate(df.iterrows()):
        sim_dic = sim_dic[1]

        df_ave = get_vave_theta_rms(sim_dic=sim_dic, crit=crit)
        time = df_ave["time"]

        ax[0].plot(time[time>0], df_ave["mass"][time>0], color=sim_dic["color"], label=sim_dic["label"], ls=sim_dic["ls"], lw=.8)

        # ---

        if do_vave: ax[1].plot(time[time>0], PBA.MomFromBeta(df_ave["vave"][time>0]), color=sim_dic["color"], ls=sim_dic["ls"], lw=.8)
        if do_mom_max:
            N = 3
            moving_average = np.convolve(df_ave["mom_max"][time>0], np.ones(N)/N, mode='valid')
            moving_average_t = np.convolve(time[time>0], np.ones(N)/N, mode='valid')
            ax[2].plot(moving_average_t, moving_average, color=sim_dic["color"], ls=sim_dic["ls"], lw=.8)
        if do_theta: ax[3].plot(time[time>0], df_ave["theta_rms"][time>0], color=sim_dic["color"], ls=sim_dic["ls"],lw=.8)
        if do_ye: ax[4].plot(time[time>0], df_ave["ye_ave"][time>0], color=sim_dic["color"], ls=sim_dic["ls"], lw=.8)
        # ax[4].plot(df_ave["time"], df_ave["entr_ave"], color=sim_dic["color"], ls=sim_dic["ls"],lw=.8)
        # ax[5].plot(df_ave["time"], df_ave["eps_ave"], color=sim_dic["color"], ls=sim_dic["ls"], lw=.8)



    for _ax in ax:
        _ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        _ax.minorticks_on()
        _ax.grid(ls=':',lw=0.8)
    # ax[0].set_ylim(*ylim)
    ax[0].set_yscale(yscale)
    ax[0].legend(fancybox=False,loc= 'upper center',columnspacing=0.4,
                 #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                 shadow=False, ncol= 2, fontsize= 12,
                 framealpha=0., borderaxespad= 0., frameon=False)
    #ax[5].set_yscale("log")
    # ax[-1].set_xlabel("Extraction time [ms]",fontsize=14)
    ax[-1].set_xlabel(r"$t_{\rm ext} - t_{\rm merg}$ [ms]", fontsize=12)
    ax[0].set_ylabel(
        r"$M_{\rm ft}$ [M$_{\odot}$]" if crit=="fast" else r"$M_{\rm ej}$ [M$_{\odot}$]",
        fontsize=12)
    if do_vave: ax[1].set_ylabel(
        r"$\langle \Gamma_{\rm ft}\beta_{\rm ft} \rangle$ [c]" if crit=="fast" else r"$\langle \Gamma\beta \rangle$ [c]",
        fontsize=12)
    if do_mom_max: ax[2].set_ylabel(
        r"$(\Gamma_{\rm ft}\beta_{\rm ft})_{\rm max}$" if crit=="fast" else r"$(\Gamma\beta)_{\rm max}$",
        fontsize=12)
    if do_theta: ax[3].set_ylabel(
        r"$\langle \theta_{\rm RMS; ft} \rangle$ [deg]" if crit=="fast" else r"$\langle \theta_{\rm RMS} \rangle$ [deg]",
        fontsize=12)
    if do_ye: ax[4].set_ylabel(
        r"$\langle Y_{e;\,\rm ft} \rangle$" if crit=="fast" else r"$\langle Y_{e;\,\rm ej} \rangle$"
        ,fontsize=12)
    ax[-1].set_xlim(*xlim)
    ax[0].set_ylim(*ylim)
    # ax[4].set_ylabel(r"$\langle s \rangle$ [$k_b/$baryon",fontsize=14)
    # ax[5].set_ylabel(r"$\langle \epsilon \rangle$",fontsize=14)
    # ax[0].set_title(title,fontsize=12)
    # ax.grid(which="both",axis="both")
    plt.tight_layout()
    plt.savefig(os.getcwd() +'/figs/' + figname, dpi=512)
    plt.savefig(os.getcwd() +'/figs/' + figname.replace(".png", ".pdf"))
    plt.show()

def print_table(crit:str="fast",maximize:str="mass") -> pd.DataFrame:
    df_res = { name : {} for name, sim_dic in df.iterrows()}
    for name, sim_dic in df.iterrows():

        # sim_dic = sim_dic[1]
        df_res[name]["name"] = sim_dic["name"]
        df_res[name]["label"] = sim_dic["label"]

        df_ave = get_vave_theta_rms(sim_dic=sim_dic, crit=crit)
        if maximize == 'mass':
            # Select text where ejecta mass is maximum
            idx = df_ave["mass"].argmax()
        elif maximize=="mom":
            idx = df_ave["mom_max"].argmax()
        else:
            raise KeyError('maximize must be mass or mom')
        # print(f"{sim_dic['name']} i={idx}")

        df_res[name]["text"] = df_ave["time"][idx] #- sim_dic["tmerg"]
        df_res[name]["mass"] = df_ave["mass"][idx]
        df_res[name]["vave"] = df_ave["vave"][idx]
        df_res[name]["mom_max"] = df_ave["mom_max"][idx]
        df_res[name]["ye_ave"] = df_ave["ye_ave"][idx]
        df_res[name]["eps_ave"] = df_ave["eps_ave"][idx]
        df_res[name]["temp_ave"] = df_ave["temp_ave"][idx]
        df_res[name]["theta_rms"] = df_ave["theta_rms"][idx]
        df_res[name]["ek"] = df_ave["ek"][idx]
    df_res = pd.DataFrame.from_dict(df_res).T
    return df_res

if __name__ == '__main__':

    print("FAST TAIL PROPERTIES AT EXTRACTION TIME WHEN MASS IS MAX ")
    df_res = print_table(crit="fast", maximize='mom')
    df_res.to_csv(os.getcwd()+'/'+"ejecta_fasttail_vals_at_massmax.csv",index=True)

