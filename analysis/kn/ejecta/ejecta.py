import numpy as np;
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

get_Gamma = lambda beta: np.float64(np.sqrt(1. / (1. - np.float64(beta) ** 2.)))

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


def piecewise_linear(x, x0, x1, y0, k1, k2, k3):
    condlist = [x < x0,
                (x >= x0) & (x < x1),
                x >= x1]
    funclist = [lambda x: k1 * x + y0 - k1 * x0,
                lambda x: k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                lambda x: k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                ]
    return np.piecewise(x, condlist, funclist)

def piecewise_power(x, x0, x1, y0, k1, k2, k3):
    condlist = [x < x0,
                (x >= x0) & (x < x1),
                x >= x1]
    funclist = [lambda x: y0*(x/x0)**k1,# k1 * x + y0 - k1 * x0,
                lambda x: y0*(x/x0)**k2,#k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                lambda x: y0*(x/x1)**k3 * (x1/x0)**k2 #k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                ]
    return np.piecewise(x, condlist, funclist)
class Fit1D:
    def __init__(self):
        pass

    def compute_chi2(self, x, y, y_pred, n):
        x = 10**x
        y_values = 10**y
        y_fit = 10**y_pred

        ss_res = np.sum((y_values - y_fit) ** 2)
        reduced_chi_squared = ss_res / (len(y_values) - n)
        return np.log10( reduced_chi_squared )

    def compute_r2(self,x,y,y_pred,n):
        x=10**x
        y_values=10**y
        y_fit =10**y_pred

        ss_res = np.sum((y_values - y_fit) ** 2)
        ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        adjusted_r_squared = 1 - (1 - r_squared) * (len(y_values) - 1) / (len(y_values) - n - 1)
        return adjusted_r_squared


    def fit(self, x_values:np.ndarray, y_values:np.ndarray):
        # Find the index of the maximum y-value
        max_index = np.argmax(y_values)
        # Corresponding x-value for the maximum point
        x_max = x_values[max_index]
        # Finding the closest index to x=0
        zero_index = np.argmin(np.abs(x_values))
        # The x-value closest to zero
        x_zero = x_values[zero_index]
        # Initial guesses for the parameters
        initial_guesses = [x_max, x_zero, max(y_values), 2, -4, -8]

        # Fit the model to the data
        popt, _ = curve_fit(piecewise_linear, x_values, y_values, p0=initial_guesses)
        print(popt)
        # Extract the optimized parameters
        x0_opt, x1_opt, y0_opt, k1_opt, k2_opt, k3_opt = popt

        # Generate x-values for plotting the fitted function
        x_fit = x_values#np.linspace(np.min(x_values), np.max(x_values), 1000)
        y_fit = piecewise_linear(x_fit, *popt)
        #
        # ss_res = np.sum((y_values - y_fit) ** 2)
        # reduced_chi_squared = ss_res / (len(y_values) - len(popt))
        reduced_chi_squared = self.compute_chi2(x_values,y_values,y_fit,len(popt))

        # ss_res = np.sum((y_values - y_fit) ** 2)
        # ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
        # r_squared = 1 - (ss_res / ss_tot)
        # adjusted_r_squared = 1 - (1 - r_squared) * (len(y_values) - 1) / (len(y_values) - len(popt) - 1)
        adjusted_r_squared = self.compute_r2(x_values,y_values,y_fit,len(popt))

        print(f"Chi2={reduced_chi_squared} R2={adjusted_r_squared}")

        # Plot the original data and the fitted curve
        # plt.figure(figsize=(10, 6))
        # plt.scatter(x_values, y_values, color='blue', label='Data Points')
        # plt.plot(x_fit, y_fit, 'r-', label='Fitted Piece-wise Linear Function')
        # plt.title('Fit of Piece-wise Linear Function to Data')
        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.axvline(x=x0_opt, color='green', linestyle='--', label='Breakpoint at x0')
        # plt.axvline(x=x1_opt, color='purple', linestyle='--', label='Breakpoint at x_zero')
        # plt.legend()
        # plt.grid(True)
        # plt.show()

        return x_fit, y_fit, dict(chi2=reduced_chi_squared, r2=adjusted_r_squared,
                                  x0=10**x0_opt,x1=10**x1_opt,y0=10**y0_opt,k1=k1_opt,k2=k2_opt,k3=k3_opt
                                  # pars=popt
                                  )

def fit_data(x, y,name:str):
    np.savetxt(os.getcwd()+'/'+name+"log_mom_log_ek.txt",X=np.column_stack((x,y)),fmt="%.3f")

    o_fit = Fit1D()
    return o_fit.fit(x, y)


# def plot_all_sim_ejecta_mass_evol(crit=None,yscale="linear",ylim=(0,0.04),title="Volume integrated ejecta mass",
#                              figname="figname"):
#     fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(7,4))
#     #ax2 = ax.twinx()
#     for idx, sim_dic in enumerate(df.iterrows()):
#         sim_dic = sim_dic[1]
#         ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)
#         if (crit is None):
#             mass = ej.total_mass()
#         else:
#             mass = ej.total_mass_vs_text(crit=crit)
#         tmerg = sim_dic["tmerg"]
#         if sim_dic["given_time"] == "new":
#             ax.plot(ej.getText()-tmerg, mass,
#                     color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"])
#         else:
#             ax.plot(ej.getText()-tmerg, mass, color=sim_dic["color"], ls=sim_dic["ls"],
#                     label=sim_dic["label"])
#         # data = PBA.id_kenta.Data(fpath_rhomax=sim_dic["datadir"]+sim_dic["rhomax"],
#         #                          fpath_mdot=sim_dic["datadir"]+sim_dic["mdot_extract"])
#         # ax.plot(ej.getText(), ej.total_mass_fasttail(),color=color_pal[sim_dic["idx"]],
#         #         label=r"$\Gamma\beta>1$",ls="--")
#         # ax2.plot(*data.get_rhomax(),color=color_pal[sim_dic["idx"]],label=r"$\rho_{\rm max}$",ls=":")
#         # ax2.plot(*data.get_mdot(),color=color_pal[sim_dic["idx"]],label=r"$\dot{M}$",ls="-.")
#     ax.tick_params(labelsize=12)
#     ax.set_ylim(*ylim)
#     ax.set_yscale(yscale)
#     ax.legend(fontsize=12,ncol=2)
#     # ax.set_yscale("log")
#     ax.set_xlabel(r"$t_{\rm ext} - t_{\rm merg}$ [ms]",fontsize=14)
#     ax.set_ylabel(r"Ejecta mass $[M_{\odot}]$",fontsize=14)
#     ax.set_title(title,fontsize=14)
#     # ax.grid(which="both",axis="both")
#     plt.tight_layout()
#     plt.savefig(os.getcwd()+f'/figs/{figname}.png',dpi=256)
#     plt.savefig(os.getcwd()+f'/figs/{figname}.pdf')
#     plt.show()
def get_text_max_mass(ej:PBA.id_kenta.EjectaData, crit:str='fast')->float:
    return float( ej.getText()[np.argmax(ej.total_mass_vs_text(crit=crit))] )
def get_text_max_mom(ej:PBA.id_kenta.EjectaData, crit:str='fast')->float:
    vals = []
    for it, text in enumerate(ej.getText()):
        ej_mass = ej.get(v_n="mass",text=text)
        mask = ej.get_vinf_mask(crit=crit)
        ej_vinf = ej.get_vinf()

        mass_ = np.sum(ej_mass[:,mask],axis=0)
        mom = ej_vinf[mask]*PBA.utils.get_Gamma(ej_vinf[mask])
        idx_ = np.argmax(mom[mass_ > 0.]) if len(mom[mass_ > 0.]) > 0 else 0

        vals.append(mom[idx_])
    return float ( ej.getText()[np.argmax(vals)] )
def get_2d_ek(ej:PBA.id_kenta.EjectaData, text:float)->tuple[np.ndarray,np.ndarray,np.ndarray]:
    thetas = ej.get_theta()
    vinf = ej.get_vinf()
    vinf = vinf[:-1]
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom = get_Gamma(vinf_c) * vinf_c
    ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    ek = ej.get(v_n="mass",text=text) #* PBA.utils.cgs.solar_m
    ek = 0.5 * (ek[:,1:] + ek[:,1:]) [:,:-1]
    ek = 0.5 * (ek[1:,:] + ek[1:,:])
    ek *= PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2
    ek = ek.T
    return (mom, ctheta, ek)
def plot_all_sim_ejecta_mass(yscale="linear",xscale="linear",ylim=(0,0.04),xlim=(1e-3,4),title="Volume integrated ejecta mass",
                             figname="figname"):
    do_cumulative = True
    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val

    fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(4.6,2*3.2),layout='constrained',sharex='col',#sharex='col',sharey='row',
                             gridspec_kw={'height_ratios': [2,2,1]})

    df_fit = pd.read_csv(os.getcwd()+'/'+'piecewise_line_fits.csv')

    infos = dict()
    idx = 0

    x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
    axins = axes[1].inset_axes(
        [0.05, 0.1, 0.5, 0.5],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

    # sub region of the original image
    x1, x2, y1, y2 = 1e-1, 5e-1, 1e48, 4e49
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.grid(ls=':')
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.set_xlabel('')
    axins.xaxis.set_visible(False)
    axins.set_xticks([])
    axins.minorticks_on()
    axins.tick_params(axis='both', which='both', labelleft=False,
                      labelright=True, tick1On=True, tick2On=True,
                      labelsize=12,
                      direction='in',
                      bottom=True, top=True, left=True, right=True)
    axes[1].indicate_inset_zoom(axins, edgecolor="black")

    # for idx, sim_dic in enumerate(df.iterrows()):
    for (sim_dic, fit_dic) in zip(df.iterrows(), df_fit.iterrows()):
        sim_dic = sim_dic[1]
        fit_dic = fit_dic[1]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)
        mom, _, ek = get_2d_ek(ej=ej,
                               # text=get_text_max_mass(ej=ej,crit='fast')
                               text=get_text_max_mom(ej=ej,crit='fast')
                               )
        ek = np.sum(ek,axis=1)
        # ej_mass = ej.get(v_n="mass",text=get_text_max_mass(ej=ej,crit='fast'))
        # ek = np.cumsum(np.sum(ek,axis=0)[::-1] * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)[::-1]/np.sum(ek * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)
        # ek = np.cumsum(ek[::-1])[::-1]#/np.sum(ek)


        l_mom, l_ek = np.log10(mom), np.log10(ek)
        mask = np.isfinite(l_ek)
        l_mom = l_mom[mask]
        l_ek = l_ek[mask]

        N = 5
        l_ek = np.convolve(l_ek, np.ones(N)/N, mode='valid')
        l_mom = np.convolve(l_mom, np.ones(N)/N, mode='valid')


        axes[0].plot(10**l_mom, get_cumulative(10**l_ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')
        axes[1].plot(10**l_mom, 10**l_ek, color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')

        l_mom_pred, l_ek_pred, info = fit_data(l_mom, l_ek,name=sim_dic["name"])
        info["label"] = sim_dic["label"]
        infos[sim_dic["name"]] = info
        axes[0].plot(10**l_mom_pred, get_cumulative(10**l_ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        axes[1].plot(10**l_mom_pred, 10**l_ek_pred, color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        axes[2].plot(10**l_mom_pred, (10**l_ek-10**l_ek_pred) / 10**l_ek, color=sim_dic["color"], ls=sim_dic["ls"], lw=0.7)

        axins.plot(10**l_mom, 10**l_ek, color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)
        axins.plot(10**l_mom_pred, 10**l_ek_pred, color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')

        '''
        y_pred = piecewise_power(10**l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
                                 k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        # y_pred = 10**piecewise_linear(l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
        #                          k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        axes[0].plot(10**l_mom, get_cumulative(y_pred), color='gray', ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        '''
        # axes[0].plot(mom, get_cumulative(ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')

        # mom_pred, ek_pred, info = fit_data(mom, ek,name=sim_dic["name"])
        # infos[sim_dic["name"]] = info
        # axes[0].plot(mom_pred, get_cumulative(ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        # axes[1].plot(mom_pred, (ek-ek_pred) / ek, color=sim_dic["color"], ls=sim_dic["ls"], lw=0.7)


    axes[0].plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-', label='Simulation',lw=1.2)#, lw=0.7, drawstyle='steps')
    axes[0].plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-',label='Fit',lw=.6)#, lw=0.7, drawstyle='steps')


    info = pd.DataFrame.from_dict(infos).T
    info.to_csv(os.getcwd()+'/piecewise_line_fits.csv')
    info["y0"] = ["${}$".format(PBA.utils.latex_float(y0)) for y0 in info["y0"]]
    del info["chi2"]
    info = info[['label', 'x0', 'x1', 'y0', 'k1', 'k2', 'k3', 'r2']]
    print(info.to_latex(float_format="%.2f",index=False))

    for ax in axes:
        ax.tick_params(labelsize=12)
        # ax.set_ylim(*ylim)

        ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        ax.minorticks_on()
        ax.grid(ls=':',lw=0.8)
        ax.set_xscale('log')

    #for idx in range(len(df)):
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[0].set_ylim(*ylim)
    axes[1].set_ylim(*ylim)
    axes[2].set_ylim(-2.5,1.1)
    axes[-1].set_xlim(*xlim)

    # ax.set_yscale(yscale)
    # ax.set_xscale(xscale)

    n = 2
    ax = axes[0]
    han, lab = ax.get_legend_handles_labels()
    ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
                            **dict(fancybox=False,loc= 'lower left',columnspacing=0.4,
                            #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                            shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
                            **dict(fancybox=False,loc= 'upper right',columnspacing=0.4,
                            #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                            shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    # ax.set_yscale("log")
    axes[-1].set_xlabel(r"$\Gamma\beta$",fontsize=12)
    axes[0].set_ylabel(r"$E_{\rm k} (> \Gamma \beta)$ [erg]",fontsize=12)
    axes[1].set_ylabel(r"$E_{\rm k}$ [erg]",fontsize=12)
    # ax.set_title(title,fontsize=12)
    axes[-1].set_ylabel(r"$\Delta E_{\rm k}$ [erg]",fontsize=12)

    plt.tight_layout()
    plt.savefig(os.getcwd()+f'/figs/{figname}.png',dpi=256)
    plt.savefig(os.getcwd()+f'/figs/{figname}.pdf')
    plt.show()

def plot_sim_ekecta_mass(sim_dic:dict, text=20):
    ej = PBA.id_kenta.EjStruct(get_ej_data(sim_dic["name"]),verbose=True)
    data = ej.get_2D_id(text=text, method_r0="from_beta",t0=1e3)
    ctheta = data["ctheta"][0,:]
    mom = data["mom"][:,0]
    ek = ej.get_2D_id(text=text, method_r0="from_beta",t0=1e3)["mass"] / PBA.utils.cgs.solar_m
    ek *= PBA.cgs.solar_m * PBA.utils.BetFromMom(mom)[:,np.newaxis]**2 * PBA.cgs.c**2

    thetas = ej.get_theta()
    vinf = ej.get_vinf()
    vinf = vinf[:-1]
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom = get_Gamma(vinf_c) * vinf_c
    ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    ek = ej.get(v_n="mass",text=text) #* PBA.utils.cgs.solar_m
    ek = 0.5 * (ek[:,1:] + ek[:,1:]) [:,:-1]
    ek = 0.5 * (ek[1:,:] + ek[1:,:])
    ek *= PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2
    ek = ek.T




    print(f"ctheta:{ctheta.shape} mom:{mom.shape} ek:{ek.shape}")

    low = 30. * np.pi / 180.
    mid = 60. * np.pi / 180.

    ek_tot = np.sum(ek, axis=1)
    ek_equatorial = np.sum(ek[:,(ctheta >= mid)], axis=1)
    ek_mid = np.sum(ek[:,(ctheta > low) & (ctheta <= mid)], axis=1)
    ek_polar = np.sum(ek[:,(ctheta <= low)], axis=1)

    default_pysr_params = dict(
        procs=4,
        populations=8,
        population_size=50,
        ncycles_per_iteration=500,
        early_stop_condition=(
            "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            # Stop early if we find a good and simple equation
        ),
        model_selection="best",
        weight_randomize=0.1, # ^ Randomize the tree much more frequently
    )
    model = PySRRegressor(
        niterations=100,
        binary_operators=["+", "*", "-", "/"],
        unary_operators=["cos", "exp", "sin", "square", "cube", "cos2(x)=cos(x)^2", "sin2(x)=sin(x)^2"],
        extra_sympy_mappings={"cos2": lambda x: sympy.cos(x)**2,
                              "sin2": lambda x: sympy.sin(x)**2},
        **default_pysr_params,
    )

    log_ek = np.log10(ek_tot)
    mom_ = mom[np.isfinite(log_ek)]
    log_ek = log_ek[np.isfinite(log_ek)]
    log_ek = log_ek[mom_>0.1]
    mom_ = mom_[mom_>0.1]
    model.fit(mom_.reshape(-1, 1),log_ek.reshape(-1, 1))#(X, y)
    print("Best:", model.sympy())
    print(model)
    # for i in range(4):
    #     print(f"Eq{i}: ", model.sympy(i))
    ypredict = model.predict(mom_.reshape(-1, 1))
    print("Default selection MSE:", np.power(ypredict - log_ek, 2).mean())
    ypredict = 10**ypredict.flatten()

    # exit(1)

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all',layout='constrained')

    ax = axes[0]
    do_cumulative = True
    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    ax.plot(mom,get_cumulative(ek_tot),color='black',label=r'Total')
    ax.plot(mom_,get_cumulative(ypredict),color='gray')
    ax.plot(mom,get_cumulative(ek_equatorial),label=r'$\theta > 60$')
    ax.plot(mom,get_cumulative(ek_mid),label=r'$30 < \theta < 60$')
    ax.plot(mom,get_cumulative(ek_polar),label=r'$\theta < 30$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2,4)
    # ax.set_ylim(1e21,1e31)
    ax.legend(fancybox=False,loc= 'lower left',columnspacing=0.4,
              #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 1, fontsize= 12,
              framealpha=0., borderaxespad= 0., frameon=False)
    ax.set_ylabel(r"$E_k$")

    ax = axes[1]
    ax.pcolormesh(mom, ctheta*180/np.pi, ek.T, cmap=plt.get_cmap('jet'),
                  norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    ax.set_xscale('log')
    ax.set_xlabel(r"$\Gamma\beta$")
    ax.set_ylabel(r"Polar angle [deg]")

    plt.show()

def plot_sim_ekecta_mass__(sim_dic:dict):
    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)

    ej = PBA.id_kenta.EjStruct(get_ej_data(sim_dic["name"]),verbose=True)
    data = ej.get_2D_id(text=ej.getText()[0], method_r0="from_beta",t0=1e3)

    ctheta = data["ctheta"][0,:]
    mom = data["mom"][:,0]
    ek = np.column_stack([np.cumsum(ej.get_2D_id(text=text, method_r0="from_beta",t0=1e3),axis=0)] for text in ej.getText())

    cmap = plt.get_cmap('jet')

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')

    # ax = axes[0]
    im = axes[0].pcolormesh(ej.getText()-sim_dic["tmerg"], ctheta, ek.T, cmap=plt.get_cmap('jet'),
                            norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=axes[0])
    plt.show()

    texts = ej.getText()
    vinf = ej.get_vinf()
    thetas = ej.get_theta()
    print(ej.get_vinf().shape, ej.get_theta().shape,ej.get(v_n="mass",text=ej.getText()[0]).shape)

    # mass_for_vin/f = lambda vinf, text : ej.get(v_n="mass",text=text)[:,find_nearest_index(ej.get_vinf(),vinf)]

    # mass_vs_vinf = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=0)
    # mass_vs_theta = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=1)

    # res = []
    # for text in texts:
    #     res.append(mass_vs_vinf(text))
    # res = np.reshape(np.array(res),newshape=(len(texts),(len(vinf))))
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom_c = get_Gamma(vinf_c) * vinf_c
    mom_c = mom_c[:-1]
    theta_c = 0.5 * (thetas[1:] + thetas[:-1])



    def get_mass_vs_vinf(angle:str):
        if angle == "equatorial":
            theta_l = 60. * np.pi / 180
            theta_h = 90. * np.pi / 180
            np.column_stack([
                np.sum(ej.get(v_n="mass",text=text)[thetas ], axis=0) for text in texts
            ]).T

    mass_vs_vinf = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=0) for text in texts]).T
    mass_vs_theta = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=1) for text in texts]).T

    mass_vs_vinf_c = 0.5 * (mass_vs_vinf[:,1:] + mass_vs_vinf[:,1:]) [:,:-1]
    mass_vs_theta_c = 0.5 * (mass_vs_theta[:,1:] + mass_vs_theta[:,1:])


    cmap = plt.get_cmap('jet')

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')

    # ax = axes[0]
    im = axes[0].pcolormesh(texts-sim_dic["tmerg"], mom_c, mass_vs_vinf_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_vinf_c.max()*1e-7, vmax=mass_vs_vinf_c.max()))

    im = axes[1].pcolormesh(texts-sim_dic["tmerg"], theta_c, mass_vs_theta_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_theta_c.max()*1e-7, vmax=mass_vs_theta_c.max()))
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])
    plt.show()


def plot_sim_ekecta_mass_old(sim_dic:dict):
    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)
    texts = ej.getText()
    vinf = ej.get_vinf()
    thetas = ej.get_theta()
    print(ej.get_vinf().shape, ej.get_theta().shape,ej.get(v_n="mass",text=ej.getText()[0]).shape)
    vinf_centers = 0.5 * (vinf[1:]+vinf[:-1])
    vinf_centers = vinf_centers[:-1]
    moms = vinf_centers*PBA.utils.get_Gamma(vinf_centers)
    ej_data_vs_mom,ej_data_vs_theta = [],[]
    for text in texts:
        ej_mass = ej.get(v_n="mass",text=text)
        ej_ek = ej_mass * PBA.utils.cgs.solar_m * (vinf[np.newaxis, :] * PBA.utils.cgs.c) ** 2
        ej_ek_centers = 0.5*(ej_ek[:,1:] + ej_ek[:,:-1])[:,:-1]
        ej_cumulative_ek_vs_mom = np.cumsum(np.sum(ej_ek_centers,axis=0)[::-1])[::-1]/np.sum(ej_ek_centers)
        ej_data_vs_mom.append(ej_cumulative_ek_vs_mom)

        ej_ek_centers = 0.5*(ej_ek[1:,:] + ej_ek[:-1,:])
        ej_data_vs_theta.append(np.sum(ej_ek_centers,axis=0)/np.sum(ej_ek_centers))

    ej_data_vs_mom = np.reshape(np.array(ej_data_vs_mom),newshape=(len(texts),len(ej_data_vs_mom[0])))
    ej_data_vs_theta = np.reshape(np.array(ej_data_vs_theta),newshape=(len(texts),len(ej_data_vs_theta[0])))
    ej_data_vs_mom[~np.isfinite(ej_data_vs_mom)] = 0.
    ej_data_vs_theta[~np.isfinite(ej_data_vs_theta)] = 0.

    levels = MaxNLocator(nbins=15).tick_values(ej_data_vs_mom.min(), ej_data_vs_mom.max())
    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = LogNorm(vmin=ej_data_vs_mom.max()*1e-7, vmax=ej_data_vs_mom.max())

    fig, axes = plt.subplots(nrows=2,ncols=1)

    ax = axes[0]
    im = ax.pcolormesh(texts-sim_dic["tmerg"], moms, ej_data_vs_mom.T, cmap=cmap, norm=norm)
    im = ax.pcolormesh(texts-sim_dic["tmerg"], thetas, ej_data_vs_theta.T, cmap=cmap, norm=norm)
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=ax)
    ax.set_title('pcolormesh with levels')
    plt.show()

def plot_rho_max(df:pd.DataFrame):

    fig,ax = plt.subplots(ncols=1,nrows=1)

    for name, sim_dic in df.iterrows():
        print(sim_dic["name"], sim_dic["rhomax"], sim_dic["mdot_extract"])
        if (sim_dic["rhomax"]!=-1):
            data = PBA.id_kenta.Data(fpath_rhomax=DATA_PATH+sim_dic["name"]+'/'+sim_dic["rhomax"],
                                     fpath_mdot=DATA_PATH+sim_dic["name"]+'/'+sim_dic["mdot_extract"])
            # df = data.df_mdot
            tmerg = sim_dic["tmerg"]
            # times = df["time"]
            # times = times - tmerg

            ax.plot(data.df_rho["time"]-tmerg, data.df_rho[f"rho_max"], color=sim_dic['color'],lw=1,alpha=.8,label=sim_dic['label'],ls=sim_dic['ls'])
    ax.grid(color='gray',ls='--',alpha=.5)
    ax.legend()
    ax.set_yscale("log")
    plt.show()



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

        if do_vave: ax[1].plot(time[time>0], df_ave["vave"][time>0]*get_Gamma(df_ave["vave"][time>0]), color=sim_dic["color"], ls=sim_dic["ls"], lw=.8)
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


# Definition of a retarded time for ejecta mass flux
def get_tret(r : float, vinf_ave0 : float =.7) -> float:
    vinf_ave0 = .7
    return (r*1000.*100./vinf_ave0/PBA.utils.cgs.c) * 1e3 # ms
def get_tret_arr(r : float, vinf_ave0 : np.ndarray) -> np.ndarray:
    # vinf_ave0 = .7
    return (r*1000.*100./vinf_ave0/PBA.utils.cgs.c) * 1e3 # ms
def plot_rho_mdot_rext(ylim=(1e-7, 1e-2), ylim2=(1, 6), xlim=(-2, 2.), crit="fast") -> None:

    # color_pal_ = sns.color_palette(n_colors=10)
    ncols = len([(name, dic) for name, dic in df.iterrows() if dic["rhomax"]!=-1])
    fig, axes = plt.subplots(ncols=ncols,
                             nrows=1, figsize=(12,2.5), sharey='all',layout='constrained')
    i = 0
    for name, sim_dic in df.iterrows():

        if sim_dic["rhomax"] == -1 or sim_dic["mdot_extract"] == -1:
            continue

        data = PBA.id_kenta.Data(fpath_rhomax=DATA_PATH+'/'+sim_dic["name"]+'/'+sim_dic["rhomax"],
                                 fpath_mdot=DATA_PATH+'/'+sim_dic["name"]+'/'+sim_dic["mdot_extract"])
        rext = sorted(data.get_rext())
        rext = rext[:-2][::-1]
        norm = Normalize(vmin=0,vmax=len(rext)-2)
        colors=['red','orange','green','blue']

        time = data.df_mdot["time"]
        #color_pal = sns.color_palette(n_colors=len(rext)+1)

        tmerg = sim_dic["tmerg"]

        ax = axes[i]
        for j, r in enumerate(rext):
            tret = get_tret(r=r) # ms
            tret = get_tret_arr(r=r, vinf_ave0=data.df_mdot[f"vave_{crit}"]) # ms
            # ax.plot(time - tret- tmerg, data.df_mdot[f"mdot_fast_r{r}"],
            #             color=color_pal_[j], label=r"$R_{\rm ext}$"+f"={r} [km]",lw=1,alpha=.9)
            vals = data.df_mdot[f"mdot_{crit}_r{r}"]
            if np.sum(vals) > 1e-10:
                ax.plot(time - tret - tmerg, vals, #color=color_pal_[j],
                        color=colors[j],#plt.cm.get_cmap('seismic')(norm(j)),
                        label=r"$R_{\rm ext}$"+f"={r}",lw=0.7,alpha=.8)

        # for i, r in enumerate(rext):
        # ax.plot(data.df_mdot["time"], data.df_mdot[f"mdot_fast_r{r}"],color=color_pal[i])
        # ax.plot(data.ret_time(r0=r,vave_key="vave_fast"), data.df_mdot[f"mdot_fast_r{r}"],color=color_pal[i])
        # df.plot("time",f"mdot_fast_r{r}",ax=ax,color=color_pal[i])


        ax2 = ax.twinx()
        ax2.plot(data.df_rho["time"]-tmerg, data.df_rho[f"rho_max"], color='black',lw=1.0,alpha=.8)

        # df.plot("time",f"rho_max",ax=ax2,color="black")

        # ax.set_xlabel("time [ms]", fontsize=12)
        #ax.set_ylabel("Mdot [Msun/s]", fontsize=12)
        #ax2.set_ylabel("rho max / rho nuc", fontsize=12)
        ax.set_yscale("log")
        # ax2.set_yscale("log")
        # ax.grid(color='black')
        # ax2.grid(color='gray',ls=':',alpha=.5)
        ax2.tick_params(labelsize=12,which='both',direction='in',labelright=True,tick2On=True,right=True,tick1On=True)
        # ax.legend(loc='lower right', fontsize=12)
        ax.tick_params(labelsize=12,which='both',direction='in')
        if len(ylim) > 0 : ax.set_ylim(*ylim)
        if len(ylim2) > 0 : ax2.set_ylim(*ylim2)
        if len(xlim) > 0 : ax.set_xlim(*xlim)

        ax.set_title(sim_dic['label'], fontsize=12)#(r"$\dot{M}$ vs $\rho_{\rm max}$" + f" for {sim_dic['label']}", fontsize=14)

        if i!=ncols-1:
            ax2.yaxis.set_ticklabels([])
        if i == 0:
            ax.legend(fancybox=False,loc= 'upper left',columnspacing=0.2,
                      #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      shadow=False, ncol= 1, fontsize= 11,
                      framealpha=0., borderaxespad= 0., frameon=False)
        i+=1

    # axes[-1].legend(loc='lower right', fontsize=12)
    axes[0].set_ylabel(r"$\dot{M}_{\rm ft}$ [$M_{\odot}$/s]", fontsize=12)
    for ax in axes:
        ax.grid(ls=':',lw=0.7)
        ax.set_xlabel(r"$t-t_{\rm merg}$ [ms]", fontsize=12)
        ax.minorticks_on()
        ax2.minorticks_on()

    ax2.set_ylabel(r"$\rho_{\rm max}/\rho_{\rm nuc}$", color="black",fontsize=12)

    # fig.text(0.5, 0.001, "time [ms]", fontsize=14, ha='center')
    # # plt.tight_layout()
    # # fig.text(0.5, 0.001, "time [ms]", fontsize=14, ha='center')
    # plt.savefig(os.getcwd()+'/figs/'+"mdot_rhomax.png",dpi=256)
    # plt.savefig(os.getcwd()+'/figs/s'+"mdot_rext_rhomax.pdf")
    # plt.show()

    # data.df_mdot.plot("time","mdot_fast_r48")
    # data.df_mdot.plot("time","mdot_fast_r144")
    # data.df_mdot.plot("time","mdot_fast_r482")
    # data.df_mdot.plot("time","mdot_fast_r144")
    plt.savefig(os.getcwd()+'/figs/'+'mdot_rhomax.pdf')
    plt.show()


if __name__ == '__main__':
    # process new simulation data -> collated.h5
    # task_process()

    ''' ejecta properties at a given extraction time '''
    # plot_all_sim_ejecta_mass_evol(crit=None,yscale="linear",ylim=(0., 0.01),figname="ejecta_mass_evol",
    #                          title="Volume integrated ejecta mass")
    plot_all_sim_ejecta_mass(figname="cumulative_ejecta_mass",
                             # ylim=(1e43,1e51),xlim=(1e-2,4),
                             ylim=(1e43,1e51),xlim=(1e-2,4),
                             yscale='log',xscale='linear',
                             title="Cumulative normalized ejecta mass distribution")
    # plot_sim_ekecta_mass(sim_dic=df.loc["SFHo_13_14_res150"])
    # plot_sim_ekecta_mass(sim_dic=df.loc["DD2_135_135_res150"])

    # plot_rho_max(df=df)

    ''' Ej. prop. as a function of extraction time '''
    # plot_all_sim_ejecta_massave_vel(yscale='linear', figsize=(4.6,3.2), xlim=(-2,35),
    #                                 do_vave=False, do_theta=False, do_ye=False,
    #                                 figname="mej_evol.png")
    # plot_all_sim_ejecta_massave_vel(crit='fast',yscale='log', figsize=(4.6,3*3.2), xlim=(-2,35), ylim=(1e-8,1e-4),
    #                                 do_vave=True, do_mom_max=True, do_theta=True, do_ye=True,
    #                                 figname="mej_evol_fast_tail.png")

    ''' Mass flux & rho_max as a function of time '''
    # plot_rho_mdot_rext()