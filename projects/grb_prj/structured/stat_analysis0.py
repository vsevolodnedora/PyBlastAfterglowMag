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

class ProcessRunsData:

    def __init__(self, workingdirs : list, key_to_key_map : dict):
        self.workingdirs = workingdirs
        if len(self.workingdirs) == 0:
            raise ValueError("No workingdirs are given")
        self.key_to_key_map = key_to_key_map
        if len(key_to_key_map.keys()) == 0:
            raise ValueError("Empty 'from_to_keys' dictionary is given")

        self.out = {}
        for in_key, out_key in self.key_to_key_map.items():
            self.out[out_key] = []
    def _process_workdir(self, i : int, do_id : bool, do_lc : bool, do_skymap : bool, verbose : bool = False):
        workingdir = self.workingdirs[i]
        pba = PBA.interface.PyBlastAfterglow(workingdir=workingdir,verbose=verbose)

        if do_skymap:
            raise NotImplementedError("not implemented yet")

        if do_lc:
            times = pba.GRB.get_lc_times(spec=False)
            freqs = pba.GRB.get_lc_freqs(spec=False, unique=True)
            lc_attrs = pba.GRB.get_lc_obj().attrs
            for freq in freqs:
                for time in times:

                    for in_key, out_key in self.key_to_key_map.items():
                        if (in_key in lc_attrs.keys()):
                            self.out[out_key].append(np.float64(lc_attrs[in_key]))

                    if do_id:
                        id_attrs = pba.GRB.get_id_obj().attrs
                        for in_key, out_key in self.key_to_key_map.items():
                            if ((in_key in id_attrs.keys()) and (not in_key in lc_attrs.keys())):
                                self.out[out_key].append( id_attrs[in_key] )


                    flux = pba.GRB.get_lc_totalflux(freq=freq,time=time,spec=False)
                    if "time" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["time"]].append( np.float64(time) )
                    if "freq" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["freq"]].append( np.float64(freq) )
                    if "flux" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["flux"]].append( np.float64(flux) )

            for _, key in self.key_to_key_map.items():
                if (len(self.out[key])) != len(times) * len(freqs) * (i + 1):
                    raise ValueError("Size mismatch")


    def process(self, do_id : bool, do_lc : bool, do_skymap : bool, verbose : bool = False):
        for i in tqdm(range(len(self.workingdirs))):
            if (not os.path.isdir(self.workingdirs[i])):
                raise IOError(f"Directory is not found: {self.workingdirs[i]}")
            self._process_workdir(i=i,do_id=do_id, do_lc=do_lc, do_skymap=do_skymap,verbose=verbose)
        print(f"Processing {len(self.workingdirs)} complete")

    def get_df(self) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(self.out)
        return df

    def save_df_as_csv(self):
        df = self.get_df()
        df.to_csv()





def process_id_files(workdirs : list[str],
                     id_fname : str,
                     frame : dict) -> dict:

    for workingdir in tqdm(workdirs):
        #print(f"Processing {i}/{len(files)}: {workingdir}")
        if not os.path.isdir(workingdir):
            raise FileNotFoundError(f"workingdir not found: {workingdir}")
        with h5py.File(workingdir+id_fname, mode="r") as f:
            for key in frame.keys():
                if key in f.attrs.keys():
                    frame[key].appen(np.float64(f.attrs[key]))

    print(f"Finished processing {len(workdirs)} files")
    return frame

def process_lc_files(workdirs : list[str],
                     frame : dict) -> dict:

    for workingdir in tqdm(workdirs):
        #print(f"Processing {i}/{len(files)}: {workingdir}")
        if not os.path.isdir(workingdir):
            raise FileNotFoundError(f"workingdir not found: {workingdir}")
        pba = PBA.interface.PyBlastAfterglow(workingdir=workingdir,verbose=False)
        times = pba.GRB.get_lc_times(spec=False)
        freqs = pba.GRB.get_lc_freqs(spec=False,unique=True)
        dfile = pba.GRB.get_lc_obj()
        for freq in freqs:
            for time in times:
                flux = pba.GRB.get_lc_totalflux(freq=freq,time=time,spec=False)
                for key in frame.keys():
                    if (not key in ["time","freq","flux"]):
                        frame[key].append(np.float64(dfile.attrs[key]))
                frame["time"].append(np.float64(time))
                frame["freq"].append(np.float64(freq))
                frame["flux"].append(np.float64(flux))
    print(f"Finished processing {len(workdirs)} files")
    return frame

class DataSimulationRuns():
    def __init__(self, runs_dir : str) -> None:
        self.df = None
        self.simafgpath = runs_dir

    def process_lc_files(self, lc_files : list[str], lc_fname : str, frame : dict):
        for file in tqdm(lc_files):
            workingdir = file.replace(lc_fname,"")
            #print(f"Processing {i}/{len(files)}: {workingdir}")
            if not os.path.isdir(workingdir):
                raise FileNotFoundError(f"workingdir not found: {workingdir}")
            pba = PBA.interface.PyBlastAfterglow(workingdir=workingdir,verbose=False)
            times = pba.GRB.get_lc_times(spec=False)
            freqs = pba.GRB.get_lc_freqs(spec=False,unique=True)
            dfile = pba.GRB.get_lc_obj()
            for freq in freqs:
                for time in times:
                    flux = pba.GRB.get_lc_totalflux(freq=freq,time=time,spec=False)
                    for key in frame.keys():
                        if (not key in ["time","freq","flux"]):
                            frame[key].append(np.float64(dfile.attrs[key]))
                    frame["time"].append(np.float64(time))
                    frame["freq"].append(np.float64(freq))
                    frame["flux"].append(np.float64(flux))

                    # convert dictionary to dataframe
        print(f"Finished processing {len(lc_files)} files")
        return frame

    # process a simulation data into a dataframe
    def process_runs(self, frame : dict, outfname : str, lc_fname : str):
        '''
            Process all files for a given simuslation
        '''
        if (not os.path.isdir(self.simafgpath)):
            raise FileNotFoundError(f"Folder not fould: {self.simafgpath}")

        # collect light curve files
        lc_files = glob(self.simafgpath + "/*/" + lc_fname)
        print(f"Files found: {len(lc_files)}")
        print(f"Example: {lc_files[0]}")

        # load 0 file
        file0 = h5py.File(lc_files[0],mode="r")
        print(file0.keys())
        print(file0.attrs.keys())
        print(f"Times = {len(np.array(file0['times']))}")
        print(f"Freqs = {len(np.array(file0['freqs']))}")

        for file in tqdm(lc_files):
            workingdir = file.replace(lc_fname,"")
            #print(f"Processing {i}/{len(files)}: {workingdir}")
            if not os.path.isdir(workingdir):
                raise FileNotFoundError(f"workingdir not found: {workingdir}")
            pba = PBA.interface.PyBlastAfterglow(workingdir=workingdir,verbose=False)
            times = pba.GRB.get_lc_times(spec=False)
            freqs = pba.GRB.get_lc_freqs(spec=False,unique=True)
            dfile = pba.GRB.get_lc_obj()
            for freq in freqs:
                for time in times:
                    flux = pba.GRB.get_lc_totalflux(freq=freq,time=time,spec=False)
                    for key in frame.keys():
                        if (not key in ["time","freq","flux"]):
                            frame[key].append(np.float64(dfile.attrs[key]))
                    frame["time"].append(np.float64(time))
                    frame["freq"].append(np.float64(freq))
                    frame["flux"].append(np.float64(flux))

                    # convert dictionary to dataframe
        print(f"Finished processing {len(lc_files)} files")
        self.df = pd.DataFrame.from_dict(frame)
        if (not outfname is None):
            print(f"Saving data into: {self.simafgpath+outfname}")
            self.df.to_csv(self.simafgpath+outfname)

    # load processed data
    def load_processed_data(self, outfname="collated.csv") -> pd.DataFrame:
        if self.df is None:
            self.df = pd.read_csv(self.simafgpath+outfname, index_col=0)
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

    runs_dir = "/media/vsevolod/T7/work/prj_grb_afterglow/struct_fsrs_resolution/"
    lc_fname = "out_lc.h5"
    id_fname = "ejecta_id.h5"

    lc_files = glob(runs_dir + "*/" + lc_fname)
    id_files = glob(runs_dir + "*/" + id_fname)
    if len(lc_files) != len(id_files):
        raise ValueError("File number mismatch")
    workingdirs = [lc_file.replace(lc_fname,"") for lc_file in lc_files]
    print(f"Number of working dirs: {len(workingdirs)}")


    print(f"Files found: {len(lc_files)}")
    print(f"Example: {lc_files[0]}")
    file0 = h5py.File(lc_files[0],mode="r")
    print(f"LC Keys: {file0.keys()}")
    print(f"LC Attrs: {file0.attrs.keys()}")
    file0.close()

    file0 = h5py.File(id_files[0],mode="r")
    print(f"ID Keys: {file0.keys()}")
    print(f"ID Attrs: {file0.attrs.keys()}")
    file0.close()

    prd = ProcessRunsData(workingdirs=workingdirs, key_to_key_map={
        "Eiso_c"    : "Eiso_c",
        "Gamma0c"   : "Gamma0c",
        "theta_c"   : "theta_c",
        "theta_w"   : "theta_w",
        "nlayers_a" : "nlayers_a",
        "n_ism"     : "n_ism",
        "time"      : "time",
        "flux"      : "flux"
    })

    prd.process(do_id=True, do_lc=True, do_skymap=False, verbose=False)

    ds = DataSimulationRuns(runs_dir="/media/vsevolod/T7/work/prj_grb_afterglow/struct_fsrs_resolution/")
    ds.process_runs(frame = { "eps_e":[], "eps_b":[], "n_ism":[],
                                           "theta_obs":[], "freq":[], "time":[], "flux":[] },
                                  outfname="collated.csv")


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