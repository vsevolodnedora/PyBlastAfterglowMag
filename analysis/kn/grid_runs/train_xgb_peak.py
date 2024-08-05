import copy
import shutil,json,os,h5py
import pandas as pd
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
from matplotlib.cm import ScalarMappable
import plotly.express as px
from multiprocessing import Pool

from scipy.interpolate import interp1d
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm,SymLogNorm
from matplotlib import cm

from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.inspection import permutation_importance
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score, KFold, ShuffleSplit
import xgboost as xgb
from optuna.samplers import TPESampler
from sklearn.model_selection import RepeatedKFold
import joblib
import gc

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

# -------------------------------------------------------------------------------

df.drop("DD2_135_135_res150_floor",axis=0,inplace=True) # No grid runs was performed for it

# -------------------------------------------------------------------------------

class Data():
    def __init__(self, working_dir,load_time_freq:bool):
        self.working_dir = working_dir
        # load training data information
        with open(self.working_dir+"train_data_info.json", 'r') as infile:
            self.info = json.load(infile)
        self.features = self.info['features']

        # create xscaler
        if (self.info["x_scaler"].__contains__(".pkl")):
            scaler = joblib.load(self.working_dir+'x_scaler.pkl')
            self.transform_x = lambda X: scaler.transform(X)
            self.inverse_transform_x = lambda X_norm: scaler.inverse_transform(X_norm)
        elif (self.info["x_scaler"] == "none"):
            self.transform_x = lambda X: X
            self.inverse_transform_x = lambda X_norm: X_norm
        else: raise NotImplementedError("Not implemented yet")

        # read yscaler
        if (self.info["y_scaler"] == "log10"):
            self.transform_y = lambda y: np.log10(y)
            self.inverse_transform_y = lambda y_norm: np.power(10., y_norm)
        else: raise NotImplementedError("Not implemented yet")

        if load_time_freq:
            # load train data meta
            with h5py.File(self.working_dir+"train_data_meta.h5", "r") as f:
                self.times = np.array(f["times"])
                self.freqs = np.array(f["freqs"])
            print(f"times={self.times.shape} freq={self.freqs.shape}")

    def _load_all_data(self):
        with h5py.File(self.working_dir+"X_train.h5", "r") as f:
            X = np.array(f["X"])
        with h5py.File(self.working_dir+"y_train.h5", "r") as f:
            y = np.array(f["y"])
        return (X, y)

    def get_normalized_train_data(self):
        X, y = self._load_all_data()
        X_norm = self.transform_x(X)
        y_norm = self.transform_y(y)
        return (X_norm, y_norm)


def prep_data(working_dir, df:pd.DataFrame, features_names:list, target:str):

    # extract X and y as arrays
    X = df.copy()
    y = np.array( X.loc[:,target] )
    X = np.array( X.loc[:,features_names] )
    # print(f"X shape: {X.shape}, y shape={y.shape}")

    # save metadata of the original file
    with h5py.File(working_dir+"train_data_meta.h5", "w") as f:
        # f.create_dataset(name="times", data=np.array( df["time"].unique()))
        f.create_dataset(name="freqs", data=np.array( df["freq"].unique()))
        f.create_dataset(name="features", data=features_names)

        # save metadata of the original file
    with h5py.File(working_dir+"X_train.h5", "w") as f:
        f.create_dataset(name="X", data=X)
    with h5py.File(working_dir+"y_train.h5", "w") as f:
        f.create_dataset(name="y", data=y)

    # normalize X
    scaler_x = preprocessing.MinMaxScaler()
    scaler_x.fit(X)
    # norm_X = scaler_x.transform(X)
    fname_x = working_dir+'x_scaler.pkl'
    joblib.dump(scaler_x, fname_x)
    # scaler_x = joblib.load(working_dir+'x_scaler.pkl')

    # normalize y
    # y = np.log10(y)

    # save info for future
    with open(working_dir+"train_data_info.json", 'w') as outfile:
        json.dump(
            {
                "target": target,
                "x_scaler": "none",
                "y_scaler": "log10",
                "features": features_names
            },
            outfile)

def process_simulation_data(sim, sim_dict,
                            features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p"),
                            text:float or None = None,
                            mode:str="flux"):
    # load full paraquet file
    fpath = os.getcwd()+'/runs/'+sim_dict['name']+'/collated.parquet'
    df_data = pd.read_parquet(fpath)
    print(f"File loaded: {fpath}")

    # remove unnecessary features
    if (not ("text" in features)) and (not (text is None)):
        print(f'Extraction times: {df_data["text"].unique()}')
        texts = list(df_data["text"].unique())
        if not text in texts:
            raise ValueError(f"Req. text={text} is not in dataframe_texts={texts}")
        df_peak = df_data.loc[df_data["text"] == text]
        df_peak.drop(["text"], axis=1, inplace=True)
    elif ("text" in features) and (text is None):
        print(f'Using text as a feature: {df_data["text"].unique()}')
    else:
        raise KeyError("Mode is not recognzied")

    # get peak data
    idx = df_data.groupby(list(features))["flux"].idxmax()
    df_peak = df_data.loc[idx]
    df_peak.to_parquet(os.getcwd()+'/runs/'+sim_dict['name']+'/collated_peaks.parquet')
    print(f"Peak data extracted {df_data.shape} -> {df_peak.shape} ")

    if mode == "time" :
        print("Removing flux data from dataframe to predict time of the peak only")
        df_peak.drop(["flux"], axis=1, inplace=True)
    elif mode == "flux":
        print("Removing time data from dataframe to predict flux of the peak only")
        df_peak.drop(["time"], axis=1, inplace=True)

    # save data (save X.h5 and Y.h5 files)
    working_dir = fpath.replace('collated.parquet',f'xgb_model_peak_{mode}_text{text}/')
    if os.path.isdir(working_dir): shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir): os.mkdir(working_dir)
    prep_data(working_dir,df=df_peak,features_names=features,target=mode)

    # load X.h5 data and Y.h5
    data = Data(working_dir=working_dir)

    # get data normalized according to the config files (from prep_data())
    X_norm, y_norm = data.get_normalized_train_data()

    # define model
    n_repeats,n_splis=1,4
    params = dict(n_estimators=50,max_depth=10,learning_rate=0.2,n_jobs=6)
    model = xgb.XGBRegressor(**params)

    # perform corss-validation
    rkf = RepeatedKFold(n_splits=n_splis, n_repeats=n_repeats, random_state=42)
    y_pred = np.zeros_like(y_norm)
    for train_index, test_index in rkf.split(X_norm):
        X_A, X_B = X_norm[train_index, :], X_norm[test_index, :]
        y_A, y_B = y_norm[train_index], y_norm[test_index]
        model.fit(
            X_A,
            y_A,
            eval_set=[(X_B, y_B)],
            verbose=0,
        )
        y_pred[test_index] += model.predict(X_B)
    y_pred /= n_repeats
    rmse = np.sqrt(mean_squared_error(y_norm, y_pred))

    print(f"{n_splis}-fold cross-validatiaon finished with RMSE={rmse}")

    # train final model using complete dataset
    model = xgb.XGBRegressor(**params)
    model.fit(X_norm, y_norm)
    rmse = np.sqrt(mean_squared_error(y_norm, model.predict(X_norm)))
    with open(working_dir+'final_model_rmse.txt', 'w') as f:
        f.write(f'rmse = {rmse:.5e}')

    # saving final model
    fpath = working_dir + "final_model.pkl"
    joblib.dump(model, fpath) # compress=0

    # performing permutation importance study
    print("Performing permutation importance analysis: ")
    perm_importance = permutation_importance(model, X_norm, y_norm)
    sorted_importances_idx = perm_importance.importances_mean.argsort()
    importances = pd.DataFrame(
        perm_importance.importances[sorted_importances_idx].T,
        columns=[data.features[idx] for idx in sorted_importances_idx],
    )
    fpath = working_dir + "perm_importances.csv"
    importances.to_csv(fpath, index=False)

    # plot permutation importance
    ax = importances.plot.box(vert=False, whis=10)
    ax.set_title("Permutation Importances (test set)")
    ax.axvline(x=0, color="k", linestyle="--")
    ax.set_xlabel("Decrease in accuracy score")
    ax.figure.tight_layout()
    plt.savefig(working_dir+"perm_importances.png",dpi=256)

    del X_norm
    del y_norm
    gc.collect()
    print(f"Analysis for {sim} is completed.")

def main():
    texts = {"BHBLpTim326_135_135_45km_150mstg_B0_HLLC":30,
             "DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS":22,
             "SFHoTim276_12_15_0025_150mstg_B15_HLLD_CT_GS_onFugaku":28,
             "SFHoTim276_13_14_0025_150mstg_B0_HLLC":26,
             "SFHoTim276_135_135_45km_150mstg_B0_FUKA":32}

    for sim,sim_dict in df.iterrows():
        process_simulation_data(sim=sim,sim_dict=sim_dict,
                                features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p"),
                                text=texts[sim_dict['name']],
                                mode='flux')
        process_simulation_data(sim=sim,sim_dict=sim_dict,
                                features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p"),
                                text=texts[sim_dict['name']],
                                mode='time')
        # break
    # exit(1)
    # special analysis of the SFHo q=1.08 data
    sim = "SFHo_13_14_res150"
    process_simulation_data(sim=sim,sim_dict=df.loc[sim],
                            features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p","text"),
                            text=None,
                            mode='flux')
    process_simulation_data(sim=sim,sim_dict=df.loc[sim],
                            features=("eps_e","eps_b","eps_t","n_ism","theta_obs","freq","p","text"),
                            text=None,
                            mode='time')

if __name__ == '__main__':
    main()