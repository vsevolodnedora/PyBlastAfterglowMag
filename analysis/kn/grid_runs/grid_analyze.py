import copy
import shutil,json,os,h5py

import joblib
import pandas as pd
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import shap
import seaborn as sns
from matplotlib.cm import ScalarMappable
import plotly.express as px
from multiprocessing import Pool

from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.interpolate import interp1d
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm,SymLogNorm
from matplotlib import cm

# try:
#     import package.src.PyBlastAfterglowMag as PBA
# except ImportError:
#     try:
#         import PyBlastAfterglowMag as PBA
#     except:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs

# --------------------------------------------------------------------------------
DATA_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/kenta_data/"
get_ej_data = lambda name : DATA_PATH+name+'/'+"ej_collated.h5"
get_runs_data = lambda name : os.getcwd()+'/runs/'+name+'/collated.parquet'

# load the metadata
with open(DATA_PATH+"metadata.json") as json_file:
    json_data = json.load(json_file)
    print(json_data)
SIMS = pd.DataFrame.from_dict(json_data).T
SIMS.set_index("name")
# select only new simulations
df = SIMS[SIMS["given_time"] == "new"]

# -------------------------------------------------------------------------------

EJ_TEXT_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/ejecta/output/"
df_text = pd.read_csv(EJ_TEXT_PATH+"ejecta_fasttail_vals_at_massmax.csv",index_col=0)
df_ej = pd.read_csv(EJ_TEXT_PATH+"ejecta_vals_at_tend.csv",index_col=0)

# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------

df.drop("DD2_135_135_res150_floor",axis=0,inplace=True) # No grid runs was performed for it

# --------------------------------------------------------------------------------

texts = {"BHBLpTim326_135_135_45km_150mstg_B0_HLLC":30,
         "DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS":22,
         "SFHoTim276_12_15_0025_150mstg_B15_HLLD_CT_GS_onFugaku":28,
         "SFHoTim276_13_14_0025_150mstg_B0_HLLC":26,
         "SFHoTim276_135_135_45km_150mstg_B0_FUKA":32}
features={
    "theta_obs":r"$\theta_{\rm obs}$",
    "n_ism":r"$n_{\rm ISM}$",
    "eps_e":r"$\epsilon_{\rm e}$",
    "eps_b":r"$\epsilon_{\rm b}$",
    "eps_t":r"$\epsilon_{\rm T}$",
    "freq":r"$\nu_{\rm obs}$",
    "time":r"$t_{\rm obs}$",
    "p":r"$p$",
    "text":r"$t_{\rm ext}$",
    "ek":r"$E_{\rm ek}$",
    "y0":r"$\mathcal{E}_0$",
    "x0":r"$\mathcal{M}_0$",
    "x1":r"$\mathcal{M}_1$",
    "k1":"$k_1$",
    "k2":"$k_2$",
    "k3":"$k_3$"
}

def convert_csv_to_paraquet(remove_csv:bool=True):
    for sim, sim_dict in df.iterrows():
        if not os.path.isfile(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.csv'):
            raise FileNotFoundError("Not found {}".format(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.csv'))

        df_ = pd.read_csv(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.csv',index_col=0)
        # print(df_)
        df_.to_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.parquet')
        # df__ = pd.read_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.parquet')
        # print(df__)
        # break
        if remove_csv:
            print("Deleting {}".format(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.csv'))
            os.remove(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.csv')



def create_peak_data(features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p"),
                     text_:float or None=None):


    for sim, sim_dict in df.iterrows():
        df_ = pd.read_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated.parquet')
        text = texts[sim_dict['name']]

        # remove unnecessary features
        if not text is None:
            print(f'Extraction times: {df_["text"].unique()}')
            df_ = df.loc[df["text"] == text]
            df_.drop(["text"], axis=1, inplace=True)
            features.remove("text")

        # groupby and save
        idx = df_.groupby(list(features))["flux"].idxmax()
        df_peak = df_.loc[idx]
        df_peak.to_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated_peaks.parquet')
        print(f"{sim} {df_.shape} -> {df_peak.shape} Done.")

def create_one_peak_data_dataframe(features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p")):
    vals = {feature:[] for feature in features}
    for sim, sim_dict in df.iterrows():
        df_ = pd.read_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated_peaks.parquet')
        for feature in features:
            vals_i = df_[feature].unique().tolist()
            if not vals[feature]: vals[feature] = vals_i
            else:
                for val in vals[feature]:
                    if not val in vals_i:
                        raise ValueError(f"Feautre={feature} expected vals={vals[feature]} gor val={val}")
    print("All clean")
def print_df(keys = ("eps_e","eps_b","eps_t","freq","p","theta_obs","n_ism"),drop_text:bool=True):
    df_runs = pd.read_parquet(get_runs_data(sim_dict['name']))
    if drop_text: df_runs.drop('text',axis=1,inplace=True)
    for key in keys:
        print(f"{key} {pd.Series(df_runs[key]).unique()}")
    # print( pd.Series(df_runs["eps_e"]).unique() )
    print(len(df_runs))
    return df_runs

def analyze_perm_importance(mode:str or None="study_flux_xgboost",
                            title="Permutation Importances (test set) for XGBoost model",
                            combined=True):


    df_imp = pd.DataFrame(columns = df.index)
    # df_imp.columns = df.keys()
    fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(4.6,3.2))
    for sim, sim_dict in df.iterrows():
        text = texts[sim_dict['name']]
        if not combined:
            fpath = os.getcwd()+'/'+'runs/'+sim_dict['name']+'/'+f'{mode}_text{text}/' + "perm_importances.csv"
        else:
            fpath = os.getcwd()+'/'+'runs/'+f'combined_{mode}/' + "perm_importances.csv"
        importances = pd.read_csv(fpath)
        importances = importances.transpose()
        # for key, val_dict in importances.iterrows():
            # print(float(list(val_dict)[0]))
            # ax.plot(key, float(list(val_dict)[0]),color='black')
        print(importances[0])
        ax.plot(
            [features[idx] for idx in importances.index],
            importances[0].values,color=sim_dict['color'],marker=sim_dict['marker'],
            ls='none',fillstyle='none',ms=12,label=sim_dict['label']
        )
        # print(importances.iloc[0])
        # df_imp[sim] = importances.iloc[0]
        # df_imp = df_imp.merge(right=importances.iloc[0])
        # print(importances)
        # ax = importances.plot.box(ax=ax,vert=True, whis=10)
        # break
    df_imp = df_imp.transpose()

    ax.set_title(title,fontsize= 12)
    ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
    # ax.tick_params(axis='x', which='minor', bottom=False)
    # ax.grid()
    # ax.minorticks_on()
    # ax.axvline(x=0, color="k", linestyle="--")
    ax.set_ylabel("Decrease in accuracy score",fontsize= 12)
    ax.legend(**dict(fancybox=False,loc= 'upper left',columnspacing=0.4,
                      #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)
               )
    # ax.figure.tight_layout()
    if not combined:fpath = os.getcwd()+'/figs/'+"perm_importances_"+mode+".pdf"
    else:fpath = os.getcwd()+'/figs/'+"perm_importances_combined_"+mode+".pdf"
    plt.savefig(fpath)
    plt.show()

def analyze_perm_importance_all(mode:str or None="study_flux_xgboost",
                            title="Permutation Importances (test set) for XGBoost model",
                            combined=True):


    df_imp = pd.DataFrame(columns = df.index)
    # df_imp.columns = df.keys()
    fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(4.6,3.2))

    fpath = os.getcwd()+'/'+'runs/'+f'combined_{mode}/' + "perm_importances.csv"
    importances = pd.read_csv(fpath)
    importances = importances.transpose()
    # for key, val_dict in importances.iterrows():
    # print(float(list(val_dict)[0]))
    # ax.plot(key, float(list(val_dict)[0]),color='black')
    print(importances[0])
    ax.plot(
        [features[idx] for idx in importances.index],
        importances[0].values,color='black',marker='s',
        ls='none',fillstyle='none',ms=12
    )
        # print(importances.iloc[0])
        # df_imp[sim] = importances.iloc[0]
        # df_imp = df_imp.merge(right=importances.iloc[0])
        # print(importances)
        # ax = importances.plot.box(ax=ax,vert=True, whis=10)
        # break
    # df_imp = df_imp.transpose()

    ax.set_title(title,fontsize= 12)
    ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
    # ax.tick_params(axis='x', which='minor', bottom=False)
    # ax.grid()
    # ax.minorticks_on()
    # ax.axvline(x=0, color="k", linestyle="--")
    ax.set_ylabel("Decrease in accuracy score",fontsize= 12)
    ax.legend(**dict(fancybox=False,loc= 'upper left',columnspacing=0.4,
                     #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                     shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)
              )
    # ax.figure.tight_layout()
    if not combined:fpath = os.getcwd()+'/figs/'+"perm_importances_"+mode+".pdf"
    else:fpath = os.getcwd()+'/figs/'+"perm_importances_combined_"+mode+".pdf"
    plt.savefig(fpath)
    plt.show()


def plot_lcs_sim(sim_:str,mode="study_flux_xgboost"):
    from train_rf_optuna import Data, plot_lcs
    tasks = [
        {"pars":{"eps_e":0.001,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":1e0}, "all_freqs":True},
        {"pars":{"eps_e":0.01,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":0.01}, "all_freqs":True},
        {"pars":{"eps_e":0.1,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":0.001}, "all_freqs":True},
    ]
    for t in tasks:
        title = [f"{features[v]}={val} " for v, val in t["pars"].items()]
        title = "".join(title)
        t["title"]=title
    for sim, sim_dict in df.iterrows():
        if not sim == sim_:
            continue
        text = texts[sim_dict['name']]
        fpath = os.getcwd()+'/'+'runs/'+sim_dict['name']+'/'+f'{mode}_text{text}/'
        data = Data(working_dir=fpath)
        model = joblib.load(fpath+"final_model.pkl")

        plot_lcs(tasks,model,data,os.getcwd()+'/figs/'+sim+"_"+mode+"_lcs",ylim=(1e-4,3e3))


def plot_shap_value(mode:str,color_bar_label:str,run:bool):
    from train_xgb_peak import Data

    # compute
    for sim, sim_dict in df.iterrows():
        text = texts[sim_dict['name']]

        fpath = os.getcwd()+'/'+'runs/'+sim_dict['name']+'/'+f'{mode}_text{text}/'

        data = Data(working_dir=fpath)

        # get data normalized according to the config files (from prep_data())
        X_norm, y_norm = data.get_normalized_train_data()

        model = joblib.load(fpath+"final_model.pkl")

        # explainer = shap.Explainer(model.predict,feature_names=data.features)
        if run:
            print("Creating Explainer")
            explainer = shap.Explainer(
                model.predict, X_norm, feature_names=[features[idx] for idx in data.features]
            )
            print(f"Computing Shap Values {sim} {mode}")
            shap_values = explainer(X_norm)
            joblib.dump(shap_values,os.getcwd()+'/output/'+f"{sim}_{mode}_shap_values.pkl")
        else:
            shap_values = joblib.load(os.getcwd()+'/output/'+f"{sim}_{mode}_shap_values.pkl")
            print("Plotting")
            print(shap_values)
            fig,ax = plt.subplots(ncols=1,nrows=1,layout='constrained',figsize=(5,3))
            plt.grid()
            shap.plots.beeswarm(shap_values,color_bar_label=color_bar_label,show=False,log_scale=False,alpha=0.6,plot_size=(5,3),color='jet')
            # shap.summary_plot(shap_values, X_norm)
            # shap.summary_plot(shap_values[0], X_norm)
            plt.savefig(os.getcwd()+'/figs/'+sim+"_"+mode+"_shap_values"+".pdf")
            plt.show()

            # shap.summary_plot(shap_values, plot_type='violin',color_bar_label=color_bar_label)
            # plt.show()
def plot_shap_value_all(mode:str,color_bar_label:str,run:bool):
    from train_xgb_peak import Data

    # compute

    fpath = os.getcwd()+'/'+'runs/'+f'combined_{mode}/'

    data = Data(working_dir=fpath)

    # get data normalized according to the config files (from prep_data())
    X_norm, y_norm = data.get_normalized_train_data()

    model = joblib.load(fpath+"final_model.pkl")

    # explainer = shap.Explainer(model.predict,feature_names=data.features)
    if run:
        print("Creating Explainer")
        explainer = shap.Explainer(
            model.predict, X_norm, feature_names=[features[idx] for idx in data.features]
        )
        print(f"Computing Shap Values {'combined'} {mode}")
        shap_values = explainer(X_norm)
        joblib.dump(shap_values,os.getcwd()+'/output/'+f"{'combined'}_{mode}_shap_values.pkl")
    else:
        shap_values = joblib.load(os.getcwd()+'/output/'+f"{'combined'}_{mode}_shap_values.pkl")
        print("Plotting")
        print(shap_values)
        fig,ax = plt.subplots(ncols=1,nrows=1,layout='constrained',figsize=(5,4))
        plt.grid()
        shap.plots.beeswarm(shap_values,color_bar_label=color_bar_label,
                            show=False,log_scale=False,alpha=0.6,plot_size=(5,3),color='jet',
                            max_display=13)
        # shap.summary_plot(shap_values, X_norm)
        # shap.summary_plot(shap_values[0], X_norm)
        plt.savefig(os.getcwd()+'/figs/'+'combined'+"_"+mode+"_shap_values"+".pdf")
        plt.savefig(os.getcwd()+'/figs/'+'combined'+"_"+mode+"_shap_values"+".png",dpi=256)
        plt.show()

        # shap.summary_plot(shap_values, plot_type='violin',color_bar_label=color_bar_label)
        # plt.show()


class GRB170817A(object):

    def __init__(self):

        # load data
        times, freqs, data, errs = self._load_data()

        # unique freqs in the data
        ufreqs = np.unique(freqs)

        # concatenate the data for unique frequencies
        # data_ord = np.array([])
        # err_ord = np.array([])
        #
        # for frequency in freqs:
        #     # print(frequency)
        #     data_ord = np.concatenate([data_ord, data[freq == frequency]])
        #     err_ord = np.concatenate([err_ord, err[freq == frequency]])
        # assert np.shape(err_ord) == np.shape(data_ord)

        # print("--- Observations for GRB170817 ---")
        # print("Total amount of time: {}".format(len(times)))
        # print("Total amount of data: {}".format(len(data)))
        # print("Total amount of errs: {}".format(len(errs)))
        # print("Unique frequencies   ({})".format(len(ufreqs)))
        # print(ufreqs)
        # print("Data per frequency:")
        # for ifreq, freq in enumerate(ufreqs):
        #     print("freq={:.2e} N={:d}".format(freq, len(data[freq==freqs])))

        self.times = times
        self.data = data
        self.errs = errs
        self.freqs = freqs
        self.ufreqs = ufreqs

    def __call__(self, freq):
        return self.get(freq)

    def _load_data(self):
        """
            Chandra21 :: https://ui.adsabs.harvard.edu/abs/2020GCN.29055....1H/abstract
            Hajela:2020
        :return:
        """

        # Updated optical and X-ray data from Fong ea 2019, and Hajela ea 2019
        time = np.array([9.2, 14.9, 16.4, 17.4, 18.3, 18.7,
                         19.4, 21.4, 22.4, 23.4, 24.2, 31.3,
                         35.3, 39.2, 46.3, 53.3, 54.3, 57.2,
                         65.9, 66.6, 67.2, 72.2, 75.5, 75.5,
                         77.6, 79.2, 80.1, 92.4, 93.1, 93.1,
                         93.2, 97.1, 107., 107., 109., 109.,
                         111., 112., 115., 115., 115., 125.,
                         125., 126., 133., 137., 149., 150.,
                         152., 158., 161., 163., 163., 163.,
                         163., 165., 167., 170., 172., 183.,
                         197., 197., 207., 209., 216., 217.,
                         217., 217., 217., 218., 218., 222.,
                         229., 252., 257., 259., 261., 267.,
                         267., 273., 273., 289., 289., 294.,
                         297., 298., 320., 324., 328., 357.,
                         359., 362., 380., 489., 545., 580.,
                         581., 741., 767., 938., 1211 # Chandra21
                         ])  # time in days

        # flux in mJy
        data = np.array(
            [5.66e-04, 6.58e-04, 1.87e+01, 1.51e+01, 1.45e+01, 1.54e+01,
             1.59e+01, 1.36e+01, 2.25e+01, 2.00e+01, 2.56e+01, 3.40e+01,
             4.40e+01, 2.28e+01, 4.40e+01, 3.20e+01, 4.80e+01, 6.10e+01,
             1.48e+02, 9.80e+01, 4.26e+01, 5.80e+01, 3.59e+01, 3.96e+01,
             7.70e+01, 4.50e+01, 4.17e+01, 3.17e+01, 9.80e+01, 7.00e+01,
             2.60e+01, 1.99e+02, 1.27e+02, 5.32e+01, 2.96e-03, 1.09e-01,
             1.11e-01, 6.29e+01, 9.62e+01, 5.12e+01, 4.12e+01, 5.82e+01,
             1.28e+02, 2.21e+02, 3.37e-03, 8.40e-02, 6.06e+01, 9.00e+01,
             1.84e+02, 3.03e-03, 2.66e-03, 9.73e+01, 6.73e+01, 4.74e+01,
             3.96e+01, 9.10e-02, 5.79e+01, 1.13e-01, 8.50e-02, 2.11e+02,
             7.59e+01, 8.93e+01, 4.20e+01, 8.20e-02, 3.63e+01, 6.05e+01,
             4.17e+01, 3.26e+01, 2.47e+01, 6.47e+01, 6.30e-02, 3.97e+01,
             4.80e+01, 7.13e+01, 4.32e+01, 1.55e-03, 6.26e+01, 2.50e+01,
             4.03e+01, 3.48e+01, 2.72e+01, 3.63e+01, 2.70e+01, 3.12e+01,
             4.40e-02, 2.34e+01, 2.31e+01, 4.72e+01, 3.40e-02, 9.70e-04,
             1.55e+01, 2.70e-02, 3.79e+01, 1.48e+01, 5.90e+00, 1.80e+01,
             3.54e-04, 2.68e-04, 4.90e+00, 1.95e-04, 3.46e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
             ])

        # frequency of obs.
        freq = np.array([2.41e+17, 2.41e+17, 3.00e+09, 3.00e+09, 3.00e+09, 7.25e+09,
                         6.20e+09, 6.20e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         1.50e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         6.70e+08, 1.30e+09, 6.00e+09, 4.50e+09, 7.35e+09, 7.35e+09,
                         1.40e+09, 4.50e+09, 6.00e+09, 7.25e+09, 1.50e+09, 3.00e+09,
                         1.50e+10, 6.70e+08, 1.30e+09, 1.30e+09, 2.41e+17, 3.80e+14,
                         5.06e+14, 6.00e+09, 3.00e+09, 1.00e+10, 1.50e+10, 7.25e+09,
                         1.30e+09, 6.70e+08, 2.41e+17, 5.06e+14, 7.25e+09, 5.10e+09,
                         1.30e+09, 2.41e+17, 2.41e+17, 3.00e+09, 6.00e+09, 1.00e+10,
                         1.50e+10, 5.06e+14, 7.25e+09, 3.80e+14, 5.06e+14, 6.50e+08,
                         3.00e+09, 1.30e+09, 5.00e+09, 5.06e+14, 1.00e+10, 3.00e+09,
                         6.00e+09, 1.00e+10, 1.50e+10, 3.00e+09, 5.06e+14, 7.25e+09,
                         4.50e+09, 1.30e+09, 3.00e+09, 2.41e+17, 1.30e+09, 7.25e+09,
                         3.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09,
                         5.06e+14, 7.25e+09, 7.25e+09, 1.30e+09, 5.06e+14, 2.41e+17,
                         7.25e+09, 5.06e+14, 1.30e+09, 3.00e+09, 6.00e+09, 7.25e+09,
                         2.41e+17, 2.41e+17, 3.00e+09, 2.41e+17, 2.41e+17 # Chandra21
                         ])

        # error on flux
        err = np.array([1.70e-04, 1.30e-04, 6.30e+00, 3.90e+00, 3.70e+00, 4.80e+00,
                        5.50e+00, 2.90e+00, 3.40e+00, 3.10e+00, 2.90e+00, 3.60e+00,
                        1.00e+01, 2.60e+00, 4.00e+00, 4.00e+00, 6.00e+00, 9.00e+00,
                        2.20e+01, 2.00e+01, 4.10e+00, 5.00e+00, 4.30e+00, 7.00e+00,
                        1.90e+01, 7.00e+00, 4.70e+00, 4.30e+00, 1.40e+01, 5.70e+00,
                        4.40e+00, 1.60e+01, 1.80e+01, 4.50e+00, 2.60e-04, 1.70e-02,
                        1.90e-02, 3.20e+00, 8.00e+00, 3.40e+00, 1.90e+00, 5.00e+00,
                        2.10e+01, 1.90e+01, 4.00e-04, 1.80e-02, 4.30e+00, 3.00e+01,
                        1.90e+01, 2.60e-04, 2.70e-04, 1.13e+01, 4.10e+00, 3.60e+00,
                        2.00e+00, 1.60e-02, 6.90e+00, 1.90e-02, 1.70e-02, 3.40e+01,
                        5.20e+00, 1.39e+01, 1.20e+01, 2.00e-02, 3.60e+00, 7.50e+00,
                        7.50e+00, 4.00e+00, 3.10e+00, 2.70e+00, 1.80e-02, 7.20e+00,
                        6.00e+00, 6.70e+00, 5.80e+00, 1.90e-04, 7.00e+00, 4.10e+00,
                        2.70e+00, 4.90e+00, 2.10e+00, 3.90e+00, 2.80e+00, 3.60e+00,
                        1.40e-02, 4.20e+00, 4.00e+00, 1.28e+01, 1.10e-02, 1.90e-04,
                        5.00e+00, 7.00e-03, 1.18e+01, 2.90e+00, 1.90e+00, 4.20e+00,
                        9.00e-05, 9.00e-05, 1.80e+00, 7.00e-05, 1.3e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
                        ])

        assert np.shape(time) == np.shape(freq)
        assert np.shape(data) == np.shape(err)

        data *= 1.e-29 # muJy -> ergs
        err *= 1.e-29 # muJy -> ergs

        return(time, freq, data, err) # [days, Hz, ergs, ergs]

    def get(self, freq):
        if not freq in self.ufreqs:
            raise NameError("freq: {} is not in data. Available are: {}".format(freq, self.ufreqs))
        mask = freq == self.freqs
        if len(mask) < 1:
            raise ValueError("Failed to get data for freq:{}".format(freq))
        return (self.times[mask], self.data[mask], self.errs[mask], np.zeros_like(self.times[mask]))

    def get_chandra(self):
        # Time since Merger (days) Fluxd (mJy) Err (mJy)
        data = np.array(
            [[2.19000, 1.40000e-07,  0.00000,     0.00000],
             [9.19679,  2.10000e-07,  8.72000e-08, -9.02000e-08],
             [15.3900,  6.44000e-07,  1.72000e-07, -1.20000e-07],
             [108.386,  2.21000e-06,  2.14000e-07, -2.02000e-07],
             [157.755,  2.41000e-06,  2.62000e-07, -2.11000e-07],
             [259.665,  1.07000e-06,  1.64000e-07, -1.55000e-07],
             [358.609,  9.12000e-07,  2.10000e-07, -1.73000e-07],
             [581.818,  2.14000e-07,  1.09000e-07, -7.72000e-08],
             [741.478,  1.26000e-07,  6.64000e-08, -5.39000e-08],
             [939.310,  1.54977e-07,  8.46022e-08, -6.26098e-08],
             [1234.11,  2.16000e-07,  5.45000e-08, -7.95000e-08]
             ])
        # return [day] [erg] [erg] [erg]
        errs = np.abs(data[:, 2:][:,::-1]).T
        # errs[:, 0] += data[:, 1]
        # errs[:, 1] += data[:, 1]
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, errs / 1e3 / 1e23, uplims)##,  data[:, 3] * 1e3 * 1e23)

    def get_vla_3ggz(self):

        data = np.array(
            [
                # Time since Merger (days) Fluxd (mJy) Err (mJy)
                [3.34000,     0.032000002,       0.0000000],
                [16.4200,     0.018700000,    0.0063000000],
                [17.3900,     0.015100000,    0.0038999999],
                [18.3300,     0.014500000,    0.0037000000],
                [22.3600,     0.022500001,    0.0034000000],
                [24.2600,     0.025599999,    0.0029000000],
                [31.3200,     0.034000002,    0.0035999999],
                [46.2600,     0.044000000,    0.0040000002],
                [54.2700,     0.048000000,    0.0060000001],
                [57.2200,     0.061000001,    0.0089999996],
                [93.1300,     0.070000000,    0.0057000001],
                ##[196.790     0.081793795,    0.0081793800],
                [196.790,     0.073213301,    0.0066756521],
                [115.200,     0.075690837,     0.019037671],
                [115.200,      0.10307770,     0.011835645],
                ##[163.070,     0.096567894,     0.020509181],
                [163.070,     0.098128255,     0.018721838],
                [216.910,     0.068999998,     0.015000000],
                [220.000,     0.064700000,    0.0027000001],
                [256.760,     0.055000000,     0.012300000],
                [267.000,     0.040300000,    0.0027000001],
                [272.670,     0.043900002,     0.010500000],
                [288.610,     0.046399999,     0.011400000],
                [294.000,     0.031199999,    0.0035999999],
                [489.000,     0.014800000,    0.0029000000],
                [767.000,    0.0049000001,    0.0018000000],
                # [1221.50,    0.0054000001,      0.00000000],
                [1243.000,   2.86*1e-3,          0.99*1e-3], # https://arxiv.org/pdf/2103.04821.pdf
                [4.6*cgs.year/cgs.day,   4.5*1e-3,          1.1*1e-3] # https://arxiv.org/pdf/2205.14788.pdf
            ]
        )
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        # data[:, 2] = 0.00200000
        uplims[-1] = True
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, data[:, 2] / 1e3 / 1e23, uplims)
def plot_grb170817a_lc(ax, freq, color):
    o_obs = GRB170817A()
    if (freq == 3e9):
        x_obs, y_obs, yerr_obs, uplims = o_obs.get_vla_3ggz()
        lbl = r"GRB170817A VLA 3\,GHz"
    elif (freq == 2.41e+17):
        x_obs, y_obs, yerr_obs, uplims = o_obs.get_chandra()
        lbl = r"GRB170817A Chandra 1\,keV"
    elif freq in o_obs.freqs:
        x_obs, y_obs, yerr_obs, uplims = o_obs.get(freq)
        lbl = "GRB170817A"
    else:
        return
    # plotdic = {}
    plotdic = {
        "x_vals": x_obs, "y_vals": y_obs * 1e23 * 1e3, "yerr_vals": yerr_obs * 1e23 * 1e3,
        "plot": {'ecolor': 'gray', 'mec': color,  "uplims":uplims,
                 'marker': 'o', 'markersize': 8,
                 "capsize": 2, "linestyle": "None", 'mfc': "white",
                 "label": lbl, "zorder": 100},
        "legend": {"fancybox": False, "loc": 'lower left',
                   # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   "shadow": "False", "ncol": 1, "fontsize": 12,
                   "framealpha": 0., "borderaxespad": 0., "frameon": False}
    }
    tmp = copy.deepcopy(plotdic)
    obs_times, obs_fluxes, obs_errs = tmp["x_vals"], tmp["y_vals"], tmp["yerr_vals"]
    tmptmp = copy.deepcopy(plotdic["plot"])
    if "label" in tmptmp.keys():
        lbl = tmptmp["label"]
        del tmptmp["label"]
    _l = ax.errorbar(obs_times, obs_fluxes, yerr=obs_errs, **tmptmp)
    # tmptmp["marker"] = "v"
    # for ulim in uplims:
    #     ax.plot(ulim)
    if "label" in plotdic["plot"].keys():
        obs_legend = copy.deepcopy(tmp["legend"])
        obs_legend["loc"] = "lower right"
        leg2 = ax.legend([_l], [lbl], **obs_legend)
        ax.add_artist(leg2)

def compare_with_grb170817(sim_:str,mode="study_flux_xgboost"):
    from train_rf_optuna import Data, plot_lcs
    o_obs = GRB170817A()
    for sim, sim_dict in df.iterrows():
        if not sim == sim_:
            continue
        text = texts[sim_dict['name']]
        fpath = os.getcwd()+'/'+'runs/'+sim_dict['name']+'/'+f'{mode}_text{text}/'
        data = Data(working_dir=fpath)
        times = data.times #


        t_obs, y_obs, _, _ = o_obs.get_vla_3ggz()
        y_obs *= 1e23 * 1e6 # erj -> muJy


if __name__ == '__main__':
    # create_peak_data()
    # create_one_peak_data_dataframe()

    # analyze_perm_importance(mode="xgb_model_peak_time",title=r"Permutation Importances for $t_{\rm peak}$",combined=False)
    # analyze_perm_importance(mode="xgb_model_peak_flux",title=r"Permutation Importances for $F_{\rm peak}$",combined=False)
    # analyze_perm_importance_all(mode="time",title=r"Permutation Importances for $t_{\rm peak}$",combined=True)
    # analyze_perm_importance_all(mode="flux",title=r"Permutation Importances for $F_{\rm peak}$",combined=True)
    # plot_lcs_sim(sim_="SFHo_13_14_res150")
    # plot_shap_value(mode="xgb_model_peak_time",
    #                 color_bar_label=r'$t_{\rm peak}$',run=False)
    # plot_shap_value(mode="xgb_model_peak_flux",
    #                 color_bar_label=r'$F_{\rm peak}$',run=False)
    plot_shap_value_all(mode="time",
                    color_bar_label=r'$t_{\rm peak}$',run=False)
    plot_shap_value_all(mode="flux",
                    color_bar_label=r'$F_{\rm peak}$',run=False)
    exit(1)

    x = np.linspace(0, 10, 100)

    for sim, sim_dict in df.iterrows():
        print_df(drop_text=True if sim != 'SFHo_13_14_res150' else False,
                 keys = ("eps_e","eps_b","eps_t","freq","p","theta_obs","n_ism") if sim != 'SFHo_13_14_res150'
                 else ("eps_e","eps_b","eps_t","freq","p","theta_obs","n_ism", "text"))
    exit(1)

        # df_runs = pd.read_parquet(get_runs_data(sim_dict['name']))
        # df_runs.drop('text',axis=1,inplace=True)
        # keys = ["eps_e","eps_b","eps_t","freq","p","theta_obs","n_ism"]
        # for key in keys:
        #     print(f"{key} {pd.Series(df_runs[key]).unique()}")
        # # print( pd.Series(df_runs["eps_e"]).unique() )
        # print(len(df_runs))
        # break
        #
        # # fig, axes = plt.subplots(figsize=(10, 10), sharex='col', sharey='row', ncols=len(df.columns), nrows=len(df.columns))
        # # print(df_runs.columns)
        #
        # # collect peak times and peak fluxes
        # keys = list(df_runs.keys())
        # keys.remove("time")
        # keys.remove("flux")
        # idx = df_runs.groupby(keys)["flux"].idxmax()
        # print(idx)
        #
        # max_fluxes = df_runs.loc[idx]
        # print(max_fluxes.shape)
        # print(max_fluxes.head())
        #
        # # for key in ["eps_e","eps_b","eps_t","freq","flux","time"]:
        # #     max_fluxes[key] = np.log10(max_fluxes[key])
        #
        # keys_ = ["eps_e","eps_b","eps_t","freq"]
        # key_ = "time"
        # max_fluxes.drop("flux",axis=1,inplace=True)
        # # max_fluxes = max_fluxes.assign(time = lambda x: np.log10(x['time']))
        # # sns.pairplot(max_fluxes, hue='time',corner=True)
        # # fig = px.parallel_coordinates(max_fluxes,
        # #                               color="time",
        # #                               # labels={col: col.replace('_', ' ') for col in df.columns},
        # #                               color_continuous_scale=px.colors.diverging.Tealrose,
        # #                               title="Parallel Coordinates Plot of Features and Target")
        # # for i,v_n_1 in enumerate(keys_):
        # #     for j,v_n_2 in enumerate(keys_[::-1]):
        # #         if i<j:
        # #             axes[i, j].axis('off')
        # #         else:
        # #             sns.pairplot(df_runs, hue='species')
        # #             axes[i, j].plot(x, np.sin((i+j) *x))
        # # pivot_table = max_fluxes.pivot_table(values='time', index='eps_b', columns='n_ism', aggfunc='mean')
        # # print(pivot_table)
        # #
        # # heatmap = sns.heatmap(pivot_table, annot=False, cmap='plasma', fmt=".1f")
        # #
        # # plt.show()
        #
        # features = ["freq","eps_e"]
        # fig, axes = plt.subplots(nrows=len(features), ncols=len(features), figsize=(12, 12),
        #                          layout='constrained',sharex='col',sharey='row')  # Adjust size as needed
        # for i,v_n_1 in enumerate(features):
        #     for j,v_n_2 in enumerate(features):
        #         # if (i==j):
        #         #     continue
        #         if i<=j:
        #             # axes[i, j].axis('off')
        #             continue
        #         else:
        #             print(v_n_1,v_n_2)
        #             pivot_table = max_fluxes.pivot_table(values='time', index=v_n_1, columns=v_n_2, aggfunc='mean')
        #             table = df_runs
        #             heatmap = sns.heatmap(pivot_table, annot=True, cmap='plasma', fmt=".2e",ax=axes[i, j],cbar=False)
        #             axes[i, j].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{float(x):.2e}'))
        #             axes[i, j].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{float(x):.2e}'))
        #             # axes[i, j].plot(x, np.sin((i+j) *x))
        #
        # # fig.show()
        # plt.show()
        # break