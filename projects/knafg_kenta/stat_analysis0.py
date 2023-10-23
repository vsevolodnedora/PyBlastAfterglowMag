import numpy as np
import h5py
import pandas as pd
from glob import glob
import os
import re
import seaborn as sns
from glob import glob
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.ticker as mticker
# settings
import PyBlastAfterglowMag as PBA
from settings import *
from paths import *
from tqdm import tqdm

import os
import warnings
from pathlib import Path

from IPython.display import display
from pandas.api.types import CategoricalDtype

from category_encoders import MEstimateEncoder
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_selection import mutual_info_regression
from sklearn.model_selection import KFold, cross_val_score
from xgboost import XGBRegressor
from sklearn.inspection import permutation_importance
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import clone
from rfpimp import plot_corr_heatmap, feature_corr_matrix
from scipy.stats import spearmanr


color_pal = sns.color_palette()
plt.style.use("fivethirtyeight")

class DataSimulationRuns():
    def __init__(self, runs_dir : str, verbose : bool = False) -> None:
        self.df = pd.DataFrame()
        self.verb = verbose
        self.simafgpath = runs_dir
        self.fname_lc = "lc_kn.h5"

    def collate_kn_lcs(self, frame : dict, outfname : str):
        '''
            1. Locate fi

            process a simulation data into a dataframe (only Kilonova and only Light Curves)
        '''
        if (not os.path.isdir(self.simafgpath)):
            raise FileNotFoundError(f"Folder not fould: {self.simafgpath}")

        # collect light curve files
        lc_files = glob(self.simafgpath + "/*/" + self.fname_lc)
        if self.verb: print(f"Files found: {len(lc_files)}")
        if self.verb: print(f"Example: {lc_files[0]}")

        # load 0 file
        file0 = h5py.File(lc_files[0],mode="r")
        if self.verb: print(file0.keys())
        if self.verb: print(file0.attrs.keys())
        if self.verb: print(f"Times = {len(np.array(file0['times']))}")
        if self.verb: print(f"Freqs = {len(np.array(file0['freqs']))}")

        for file in tqdm(lc_files):
            workingdir = file.replace(self.fname_lc,"")
            #print(f"Processing {i}/{len(files)}: {workingdir}")
            if not os.path.isdir(workingdir):
                raise FileNotFoundError(f"workingdir not found: {workingdir}")
            pba = PBA.interface.PyBlastAfterglow(workingdir=workingdir,verbose=False)
            times = pba.KN.get_lc_times(spec=False)
            freqs = pba.KN.get_lc_freqs(spec=False,unique=True)
            dfile = pba.KN.get_lc_obj()
            for freq in freqs:
                for time in times:
                    flux = pba.KN.get_lc_totalflux(freq=freq,time=time,spec=False)
                    for key in frame.keys():
                        if (not key in ["time","freq","flux"]):
                            frame[key].append(np.float64(dfile.attrs[key]))
                    frame["time"].append(np.float64(time))
                    frame["freq"].append(np.float64(freq))
                    frame["flux"].append(np.float64(flux))

        # convert dictionary to dataframe
        if self.verb: print(f"Finished processing {len(lc_files)} files")
        self.df = pd.DataFrame.from_dict(frame)
        if (not outfname is None):
            if self.verb: print(f"Saving data into: {self.simafgpath+outfname}")
            self.df.to_csv(self.simafgpath+outfname)

    # load processed data
    def get_df(self, outfname : str = "collated.csv", index_col : int = 0) -> pd.DataFrame:
        if len(self.df.keys()) == 0:
            if self.verb:
                print(f"Loading data: {self.simafgpath+outfname}")
            self.df = pd.read_csv(self.simafgpath+outfname, index_col=index_col)
        return self.df

class StatAnalysis():

    def __init__(self, df : pd.DataFrame, features : list[str]):
        self.df = df
        self.features = features

    def train_rf(self, key : str, rf : RandomForestRegressor) -> RandomForestRegressor:

        X = self.df.copy()
        y = X.pop(key)
        X = X.loc[:,self.features]

        rf.fit(X, y)
        print(f"RF train accuracy: {rf.score(X, y):.3f}")

        return rf

    def make_permutation_importance_plot(self, key : str, ax : plt.axes or None,
                                         rf_fitted : RandomForestRegressor,
                                         n_repeats=20, random_state=42, n_jobs=2) \
            -> tuple[plt.axes,pd.DataFrame] :

        if ax is None:
            f, ax = plt.subplots(1,1,figsize=(15,5))

        X = self.df.copy()
        y = X.pop(key)
        X = X.loc[:,self.features]

        result_train = permutation_importance(
            rf_fitted, X, y,
            n_repeats=n_repeats, random_state=random_state, n_jobs=n_jobs
        )

        sorted_importances_idx_train = result_train.importances_mean.argsort()

        importances_train = pd.DataFrame(
            result_train.importances[sorted_importances_idx_train].T,
            columns=X.columns[sorted_importances_idx_train],
        )

        importances_train.plot.box(vert=False, whis=10, ax = ax)
        ax.set_title("Permutation Importances (train set)")
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Decrease in accuracy score")
        if not ax is None:
            ax.figure.tight_layout()
            plt.show()

        return (ax, importances_train)

    def make_dropcol_importances(self, key : str,
                                 rf : RandomForestRegressor) -> pd.DataFrame:

        X_init = self.df.copy()
        y_init = X_init.pop(key)
        X_init = X_init.loc[:,self.features]

        rf_ = clone(rf)
        rf_.random_state = 42
        rf_.fit(X_init, y_init)

        #use out of bag error as performance measurement
        baseline = rf_.oob_score_
        imp = []
        for col in X_init.columns:
            X = X_init.drop(col, axis=1)
            rf_ = clone(rf)
            rf_.random_state = 42
            rf_.fit(X, y_init)
            o = rf_.oob_score_
            imp.append(baseline - o)
        imp = np.array(imp)
        I = pd.DataFrame(
            data={'Feature':X_init.columns,
                  'Importance':imp})
        I = I.set_index('Feature')
        I = I.sort_values('Importance', ascending=True)
        return I

def analysis():
    df = pd.read_csv("/media/vsevolod/T7/work/afterglowKentaProject/"+
                     "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + '/' +
                     "collated.csv", index_col=0)

    log_features = ["eps_e","eps_b","n_ism","freq"]
    df[log_features] = df[log_features].apply(np.log10)

    X = df.loc[:,['eps_e', 'eps_b', 'n_ism', 'theta_obs', 'freq']]
    print(X.shape)
    corr, pvals = spearmanr(X)
    print(repr(corr))

    print(feature_corr_matrix(X))

    viz = plot_corr_heatmap(df=X, figsize=(7,5), precision=10)
    viz.view()


    print(df.keys())
    ann = StatAnalysis(df = df, features=['eps_e', 'eps_b', 'n_ism', 'theta_obs', 'freq'])

    rf = RandomForestRegressor(
        n_estimators=200,
        n_jobs=-1,
        min_samples_leaf = 20,
        oob_score=True,
        random_state = 42)

    imp = ann.make_dropcol_importances(key = "flux", rf = rf)
    imp.plot(kind = 'barh')
    plt.show()

    rf = ann.train_rf(key="flux", rf=rf)

    ann.make_permutation_importance_plot(key="flux",
                                         ax=None,
                                         rf_fitted=rf)

def main():
    analysis()

if __name__ == '__main__':
    main()