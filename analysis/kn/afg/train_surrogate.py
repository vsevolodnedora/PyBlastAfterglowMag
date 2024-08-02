import pandas as pd
import os.path
import os
import h5py
import copy
import numpy as np
import tqdm
import itertools
import joblib
import gc
import optuna
from optuna.trial import TrialState
from optuna.integration import XGBoostPruningCallback
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LogNorm, Normalize
import json

from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

from sklearn.inspection import permutation_importance
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score, KFold, ShuffleSplit
import xgboost as xgb
from optuna.samplers import TPESampler
from sklearn.model_selection import RepeatedKFold
from xgboost import XGBRegressor
from sklearn.ensemble import RandomForestRegressor

def prep_data(working_dir, df:pd.DataFrame, features_names:list):

    # extract X and y as arrays
    X = df.copy()
    y = np.array( X.loc[:,"flux"] )
    X = np.array( X.loc[:,features_names] )
    print(f"X shape: {X.shape}, y shape={y.shape}")

    # save metadata of the original file
    with h5py.File(working_dir+"train_data_meta.h5", "w") as f:
        f.create_dataset(name="times", data=np.array( df["time"].unique()))
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
                "target": "flux",
                "x_scaler": "none",
                "y_scaler": "log10",
                "features": features_names
            },
            outfile)

class Data():
    def __init__(self, working_dir):
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

    def get_x_for_pars(self, pars:dict):
        pars["time"] = 0.
        pars = [pars[key] for key in self.features] # to preserve the order for inference
        _pars = np.vstack(([np.array(pars) for _ in range(len(self.times))]))
        _pars[:,-1] = self.times
        # print(_pars)
        return self.transform_x( _pars )

    def get_data_lc(self, pars):
        X, y = self._load_all_data()
        mask = np.ones_like(y, dtype=bool)
        x_ = self.get_x_for_pars(pars)
        # run for all features except time
        for i in range(len(self.features)-1):
            i_mask = X[:,i] == x_[0,i]
            mask = mask & i_mask
        # assert np.sum(mask) == len(self.times)
        lc = y[mask]
        return lc



class RFObjective():
    base_pars={
        "random_state":23,
        "oob_score":False,
    }
    def __init__(self, X, y, n_jobs):
        self.X = X
        self.y = y
        self.n_jobs = n_jobs


    def model(self, params:dict, n_jobs:int):
        params["n_jobs"] = n_jobs
        return RandomForestRegressor(**params)
    def train_eval_model(self, params:dict, n_splits=10, n_repeats=1):
        model = self.model(params, int(self.n_jobs / 4) if self.n_jobs >= 4 else self.n_jobs)
        # K-Fold Cross-Validation
        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        # Perform cross-validation and return the average score
        scores = cross_val_score(model, self.X, self.y,
                                 cv=kfold,
                                 scoring='neg_mean_squared_error',
                                 n_jobs=4)
        return (-1.0 * np.mean(scores))
        #
        # rkf = RepeatedKFold(
        #     n_splits=n_splits, n_repeats=n_repeats, random_state=42
        # )
        # X_values = self.X
        # y_values = self.y
        # y_pred = np.zeros_like(y_values)
        # for train_index, test_index in rkf.split(X_values):
        #     X_A, X_B = X_values[train_index, :], X_values[test_index, :]
        #     y_A, y_B = y_values[train_index], y_values[test_index]
        #     model.fit(
        #         X_A,
        #         y_A
        #         # eval_set=[(X_B, y_B)],
        #         # verbose=0,
        #     )
        #     y_pred[test_index] += model.predict(X_B)
        # y_pred /= n_repeats
        # pruning_callback = XGBoostPruningCallback(trial, 'valid-aft-nloglik')

        # return np.sqrt(mean_squared_error(self.y, y_pred))

    def __call__(self, trial:optuna.trial.Trial):
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 5, 100),
            "max_depth": trial.suggest_int("max_depth", 2, 20),
            'max_features': trial.suggest_int("max_features", 2, 20),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 15),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 2, 15),
        }
        params.update(copy.deepcopy(self.base_pars))

        # pruning_callback = XGBoostPruningCallback(trial, "validation_0-rmse")

        # params["callbacks"] = [pruning_callback]

        rmse = self.train_eval_model(params=params)
        print(f"\t trial={trial.number} model=rf njobs={self.n_jobs} mean(rmse)={rmse}")
        return rmse

class XGBoostObjetive():
    """
        Class for optuna_study training RF model
    """
    base_pars = {
        "verbosity": 0,  # 0 (silent) - 3 (debug)
        "objective": "reg:squarederror",
        "seed": 42,
        "eval_metric":"rmse"
        # "early_stopping_rounds":50,
    }
    def __init__(self, X, y, n_jobs):
        self.X = X
        self.y = y
        self.n_jobs = n_jobs

    def model(self,params:dict, n_jobs:int):
        params["n_jobs"] = n_jobs
        return xgb.XGBRegressor(**params)

    def train_eval_model(self, params:dict, n_splits=10, n_repeats=1):
        model=self.model(params, self.n_jobs)
        rkf = RepeatedKFold(
            n_splits=n_splits, n_repeats=n_repeats, random_state=42
        )
        X_values = self.X
        y_values = self.y
        y_pred = np.zeros_like(y_values)
        for train_index, test_index in rkf.split(X_values):
            X_A, X_B = X_values[train_index, :], X_values[test_index, :]
            y_A, y_B = y_values[train_index], y_values[test_index]
            model.fit(
                X_A,
                y_A,
                eval_set=[(X_B, y_B)],
                verbose=0,
            )
            y_pred[test_index] += model.predict(X_B)
        y_pred /= n_repeats
        # pruning_callback = XGBoostPruningCallback(trial, 'valid-aft-nloglik')

        return np.sqrt(mean_squared_error(self.y, y_pred))

    def __call__(self, trial:optuna.trial.Trial):
        params = {
            "n_estimators": trial.suggest_int("n_estimators", 10, 700),
            "max_depth": trial.suggest_int("max_depth", 4, 20),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.99),
            # "colsample_bytree": trial.suggest_float("colsample_bytree", 0.3, 0.7),
            "subsample": trial.suggest_float("subsample", 0.01, 1),
            "alpha": trial.suggest_float("alpha", 0., 100.),
            "lambda": trial.suggest_float("lambda", 0., 100.),
            # "gamma": trial.suggest_float("gamma", 1e-7, 1.0,log=True),
            "min_child_weight": trial.suggest_float("min_child_weight", 0, 1000)
        }
        params.update(copy.deepcopy(self.base_pars))

        pruning_callback = XGBoostPruningCallback(trial, "validation_0-rmse")

        params["callbacks"] = [pruning_callback]

        rmse = self.train_eval_model(params=params)
        print(f"\t trial={trial.number} model=xboost njobs={self.n_jobs} rmse={rmse}")
        return rmse
        # # Define hyperparameter search space
        # base_params = {
        #     'verbosity': 0,
        #     'objective': 'reg:squarederror',
        #     'eval_metric': 'rmse',
        #     'tree_method': 'hist' # Faster histogram optimized approximate greedy algorithm.
        # }  # Hyperparameters common to all trials
        # params = {
        #     'learning_rate': trial.suggest_loguniform('learning_rate', 0.01, 1.0),
        #     'max_depth': trial.suggest_int('max_depth', 6, 10), # Extremely prone to overfitting!
        #     'n_estimators': trial.suggest_int('n_estimators', 400, 4000, 400), # Extremely prone to overfitting!
        #     'eta': trial.suggest_float('eta', 0.007, 0.013), # Most important parameter.
        #     'subsample': trial.suggest_discrete_uniform('subsample', 0.2, 0.9, 0.1),
        #     'colsample_bytree': trial.suggest_discrete_uniform('colsample_bytree', 0.2, 0.9, 0.1),
        #     'colsample_bylevel': tarial.suggest_discrete_uniform('colsample_bylevel', 0.2, 0.9, 0.1),
        #     'min_child_weight': trial.suggest_loguniform('min_child_weight', 1e-4, 1e4), # I've had trouble with LB score until tuning this.
        #     'reg_lambda': trial.suggest_loguniform('reg_lambda', 1e-4, 1e4), # L2 regularization
        #     'reg_alpha': trial.suggest_loguniform('reg_alpha', 1e-4, 1e4), # L1 regularization
        #     'gamma': trial.suggest_loguniform('gamma', 1e-4, 1e4),
        # }
        # params.update(base_params)
        # pruning_callback = optuna.integration.XGBoostPruningCallback(
        #     trial, f'valid-{base_params["eval_metric"]}'
        # )
        #
        # bst = xgb.train(params, self.dtrain, num_boost_round=10000,
        #                 evals=[(self.dtrain, 'train'),
        #                        (self.dvalid, 'valid')],
        #                 early_stopping_rounds=50,
        #                 verbose_eval=False,
        #                 callbacks=[pruning_callback])
        # if bst.best_iteration >= 25:
        #     return bst.best_score
        # else:
        #     return np.inf  # Reject models with < 25 trees

def save_study_results(working_dir:str, study:optuna.study.Study):

    # Find number of pruned and completed trials
    pruned_trials = study.get_trials(deepcopy=False, states=[TrialState.PRUNED])
    complete_trials = study.get_trials(deepcopy=False, states=[TrialState.COMPLETE])

    with open(working_dir + "summary.txt", 'w') as file:
        # Display the study statistics
        file.write("\nStudy statistics: \n")
        file.write(f"  Number of finished trials: {len(study.trials)}\n")
        file.write(f"  Number of pruned trials: {len(pruned_trials)}\n")
        file.write(f"  Number of complete trials: {len(complete_trials)}\n")

        trial = study.best_trial
        file.write("\nBest trial:\n")
        file.write(f" Value: {trial.value}\n")
        file.write(f" Numer: {trial.number}\n")
        file.write("  Params: \n")
        for key, value in trial.params.items():
            file.write("    {}: {}\n".format(key, value))

        # Find the most important hyperparameters
        most_important_parameters = optuna.importance.get_param_importances(study, target=None)
        # Display the most important hyperparameters
        file.write('\nMost important hyperparameters:\n')
        for key, value in most_important_parameters.items():
            file.write('  {}:{}{:.2f}%\n'.format(key, (15-len(key))*' ', value*100))

    # Save results to csv file
    df = study.trials_dataframe().drop(['datetime_start',
                                        'datetime_complete',
                                        'duration'], axis=1)  # Exclude columns
    df = df.loc[df['state'] == 'COMPLETE']        # Keep only results that did not prune
    df = df.drop('state', axis=1)                 # Exclude state column
    df = df.sort_values('value')                  # Sort based on accuracy
    df.to_csv(working_dir + 'optuna_results.csv', index=False)  # Save to csv file
    # Display results in a dataframe
    print("\nOverall Results (ordered by loss):\n {}".format(df))


def plot_lcs(tasks, model, data:Data, figfpath):
    times = data.times
    freqs = data.freqs
    norm = LogNorm(vmin=np.min(freqs),
                   vmax=np.max(freqs))
    cmap_name = "coolwarm_r" # 'coolwarm_r'
    cmap = plt.get_cmap(cmap_name)

    fig, axes = plt.subplots(2, len(tasks), figsize=(15,5), sharex='all',
                             gridspec_kw={'height_ratios': [3, 1]})

    for i, task in enumerate(tasks):
        # pars = [req_pars[feat] for feat in features_names]
        l2s = []
        for freq in freqs:
            req_pars = copy.deepcopy(task["pars"])
            req_pars["freq"] = freq

            # Create a boolean mask based on the dictionary
            # mask = pd.Series(True, index=df.index)
            # for col, value in req_pars.items():
            #     i_mask = df[col] == value
            #     mask = mask & i_mask
            #
            # lc = np.log10( np.array(df["flux"][mask]) ).flatten()
            # x_ = data.get_x_for_pars(req_pars)
            # X, y = data._load_all_data()w

            lc = np.log10( data.get_data_lc(req_pars) )

            if not len(lc) == len(times):
                raise ValueError(f"Size mismatch for lc={lc.shape}, times={times.shape}")

            l11, = axes[0, i].plot(times/86400, lc,ls='-', color='gray', alpha=0.5, label=f"Original") # for second legend
            l21, = axes[0, i].plot(times/86400, lc,ls='-', color=cmap(norm(freq)), alpha=0.5, label=f"{freq/1e9:.1f} GHz")
            l2s.append(l21)

            # # get light curve from model
            # lc_nn = np.log10( inference(rf, req_pars, times) )
            lc_nn = model.predict( data.get_x_for_pars(req_pars) )
            lc_nn = data.inverse_transform_y( lc_nn )
            lc_nn = np.log10(lc_nn)
            # print("hi")

            l12, = axes[0, i].plot(times/86400, lc_nn,ls='--', color='gray', alpha=0.5, label=f"cVAE") # for second legend
            l22, = axes[0, i].plot(times/86400, lc_nn,ls='--', color=cmap(norm(freq)), alpha=0.5)


            # plot difference
            axes[1, i].plot(times/86400, lc-lc_nn, ls='-', color=cmap(norm(freq)), alpha=0.5)
            # del mask
        axes[0, i].set_xscale("log")
        axes[1, i].set_xlabel(r"times [day]")

    axes[0,0].set_ylabel(r'$\log_{10}(F_{\nu})$')
    axes[1,0].set_ylabel(r'$\Delta\log_{10}(F_{\nu})$')

    first_legend = axes[0,0].legend(handles=l2s, loc='upper right')

    axes[0,0].add_artist(first_legend)
    axes[0,0].legend(handles=[l11,l12], loc='lower right')

    plt.tight_layout()
    plt.savefig(figfpath,dpi=256)
    # plt.show()

def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx
def plot_violin(data:Data, model, figfpath):

    req_times = np.array([0.1, 1., 10., 100., 1000., 1e4]) * 86400.

    X, y = data._load_all_data()
    y = data.transform_y(y)

    times = data.times
    lcs = np.reshape(y,
                     newshape=(len(y)//len(times),
                               len(times)))
    yhat = model.predict(X)
    lcs_nn = np.reshape(yhat,
                        newshape=(len(yhat)//len(times),
                                  len(times)))
    allfreqs = X[:,data.features.index("freq")]
    allfreqs = np.reshape(allfreqs,
                          newshape=(len(allfreqs)//len(times),
                                    len(times)))


    cmap_name = "coolwarm_r" # 'coolwarm_r'
    cmap = plt.get_cmap(cmap_name)

    # freqs = np.array(df["freq"].unique())
    # times = np.array(df["time"].unique())
    norm = LogNorm(vmin=np.min(data.freqs),
                   vmax=np.max(data.freqs))

    # log_lcs = np.log10(lcs)
    # log_nn_lcs = np.log10(lcs_nn)
    delta = lcs - lcs_nn
    print(delta.shape)

    fig, ax = plt.subplots(2, 3, figsize=(12, 5), sharex="all", sharey="all")
    ax = ax.flatten()
    for ifreq, freq in enumerate(data.freqs):
        i_mask1 = (allfreqs == freq).astype(int)#X[]#(np.array(df["freq"]) == freq).astype(bool)

        _delta = delta * i_mask1

        time_indeces = [find_nearest_index(times, t) for t in req_times]
        _delta = _delta[:, time_indeces]


        color = cmap(norm(data.freqs[0]))

        if np.sum(_delta) == 0:
            raise ValueError(f"np.sum(delta) == 0 delta={_delta.shape}")
        # print(_delta.shape)
        violin = ax[ifreq].violinplot(_delta, positions=range(len(req_times)),
                                      showextrema=False, showmedians=True)


        for pc in violin['bodies']:
            pc.set_facecolor(color)
        violin['cmedians'].set_color(color)
        for it, t in enumerate(req_times):
            ax[ifreq].vlines(it, np.quantile(_delta[:,it], 0.025), np.quantile(_delta[:,it], 0.975),
                             color=color, linestyle='-', alpha=.8)

        # ax[ifreq].hlines([-1,0,1], 0.1, 6.5, colors='gray', linestyles=['dashed', 'dotted', 'dashed'], alpha=0.5)


        ax[ifreq].set_xticks(np.arange(0, len(req_times)))
        # print(ax[ifreq].get_xticklabels(), ax[ifreq])
        _str = lambda t : '{:.1f}'.format(t/86400.) if t/86400. < 1 else '{:.0f}'.format(t/86400.)
        ax[ifreq].set_xticklabels([_str(t) for t in req_times])

        ax[ifreq].annotate(f"{freq/1e9:.1f} GHz", xy=(1, 1),xycoords='axes fraction',
                           fontsize=12, horizontalalignment='right', verticalalignment='bottom')

        ax[ifreq].set_ylim(-0.5,0.5)


    # Create the new axis for marginal X and Y labels
    ax = fig.add_subplot(111, frameon=False)

    # Disable ticks. using ax.tick_params() works as well
    ax.set_xticks([])
    ax.set_yticks([])

    # Set X and Y label. Add labelpad so that the text does not overlap the ticks
    ax.set_xlabel(r"Time [days]", labelpad=20, fontsize=12)
    ax.set_ylabel(r"$\Delta \log_{10}(F_{\nu})$", labelpad=40, fontsize=12)

    plt.tight_layout()
    plt.savefig(figfpath,dpi=256)
    # plt.show()

def main_run(fpath_to_collated:str=os.getcwd()+'/'+"collated.csv",
             path_to_save_run:str=os.getcwd()+'/',
             target_key:str="flux",
             model_name:str="rf",
             text:int or None=32,
             sample:int or None=None,
             n_jobs:int=32,
             study_name="example-study"):
    """
    fpath_to_collated = "/media/vsevolod/T7/work/prj_kn_afterglow/SFHoTim276_135_135_45km_150mstg_B0_FUKA/collated.csv"

    :return:
    """

    # Load collated csv file
    assert os.path.isfile(fpath_to_collated), "Collated file not found"
    df = pd.read_csv(fpath_to_collated, index_col=0)
    print(f"File loaded: {fpath_to_collated} {print(df.info(memory_usage='deep'))}")

    # select a sample of data for testing
    # if not sample is None:
    #     df = df.sample(n=sample)
    df = df.head(n=sample)
    # get feature names
    features_names = [col for col in list(df.columns) if col not in [target_key]] # here time is included

    # remove unnecessary features
    if not text is None:
        print(f'Extraction times: {df["text"].unique()}')
        df = df.loc[df["text"] == text]
        df.drop(["text"], axis=1, inplace=True)
        features_names.remove("text")

    # create directory for the current run
    optuna_study_dir = path_to_save_run + f'study_{target_key}_{model_name}_text{text}/'
    if not os.path.isdir(optuna_study_dir): os.mkdir(optuna_study_dir)

    # save data for study (save X.h5 and Y.h5 files)
    prep_data(working_dir=optuna_study_dir, df=df, features_names=features_names)

    # load X.h5 data and Y.h5
    data = Data(working_dir=optuna_study_dir)

    # get data normalized according to the config files (from prep_data())
    X_norm, y_norm = data.get_normalized_train_data()

    # define objective function (__call__()) to be optimized
    if model_name == "xgboost":
        objective = XGBoostObjetive(X=X_norm, y=y_norm, n_jobs=n_jobs)
    elif model_name == "rf":
        objective = RFObjective(X=X_norm, y=y_norm, n_jobs=n_jobs)
    else:
        raise KeyError(f"Model name is not recognized: {model_name}")


    # instantiate the study
    fpath = optuna_study_dir + "study.pkl"
    if not os.path.isfile(fpath):
        print(f"Performing study... {fpath}")
        sampler = TPESampler(seed=42)
        study = optuna.create_study(
            study_name=study_name,
            direction='minimize',
            # pruner=optuna.pruners.MedianPruner(n_warmup_steps=5),
            sampler=sampler
        )
        # study = optuna.create_study(sampler=optuna.samplers.GridSampler(objective.search_space))


        # run all trials
        study.optimize(objective,
                       n_trials=200,
                       callbacks=[lambda study, trial: gc.collect()])
        # save the whole study
        joblib.dump(study, fpath)
    else:
        print(f"Study already exists. Loading... {fpath}")
        study = joblib.load(fpath)

    # record results
    save_study_results(working_dir=optuna_study_dir, study=study)

    params = {}
    params.update(objective.base_pars)
    params.update(study.best_trial.params)



    # Re-run training with the best hyperparameter combination
    fpath = optuna_study_dir + "final_model.pkl"
    if not os.path.isfile(fpath):
        print('Training model on the entire dataset = {}'.format(params))
        model = objective.model(params, n_jobs=n_jobs)
        model.fit(X_norm, y_norm)
        rmse = np.sqrt(mean_squared_error(y_norm, model.predict(X_norm)))
        with open(optuna_study_dir+'final_model_rmse.txt', 'w') as f:
            f.write(f'rmse = {rmse:.5e}')
        # model, rmse = objective.train_eval_model(params=params)
        print('Saving the best trial... rmse = {}'.format(rmse))
        # model.save_model('best_model.json')
        joblib.dump(model, fpath) # compress=0
    else:
        print(f"Best trial model already exists. Loading... {fpath}")
        model = joblib.load(fpath)


    # model = XGBRegressor()/RandomForestRegressor()
    # model.load_model(working_dir+'best_model.json')
    # data = Data(working_dir=working_dir)

    # plot model performance
    if sample is None:
        plot_violin(data, model, optuna_study_dir+"violin.png")
        tasks = [
            {"pars":{"eps_e":0.001,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":1e0}, "all_freqs":True},
            {"pars":{"eps_e":0.01,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":0.01}, "all_freqs":True},
            {"pars":{"eps_e":0.1,"eps_b":0.01,"eps_t":1.,"p":2.2,"theta_obs":0.,"n_ism":0.001}, "all_freqs":True},
        ]
        plot_lcs(tasks, model, data, optuna_study_dir+"lcs.png")



    # compute permutation importance
    fpath = optuna_study_dir + "perm_importances.csv"
    if not os.path.isfile(fpath):
        print("Performing permutation importance analysis: ")
        perm_importance = permutation_importance(model, X_norm, y_norm)
        # sorted_idx = perm_importance.importances_mean.argsort()

        sorted_importances_idx = perm_importance.importances_mean.argsort()
        importances = pd.DataFrame(
            perm_importance.importances[sorted_importances_idx].T,
            columns=[data.features[idx] for idx in sorted_importances_idx],
        )
        importances.to_csv(fpath, index=False)
    else:
        print(f"Permutation importance is found. Loading... {fpath}")
        importances = pd.read_csv(fpath)

    # plot permutation importance
    ax = importances.plot.box(vert=False, whis=10)
    ax.set_title("Permutation Importances (test set)")
    ax.axvline(x=0, color="k", linestyle="--")
    ax.set_xlabel("Decrease in accuracy score")
    ax.figure.tight_layout()
    plt.savefig(optuna_study_dir+"perm_importances.png",dpi=256)

    data = None
    del X_norm
    del y_norm
    gc.collect()

if __name__ == '__main__':
    tasks = [
        {"name":"BHBLpTim326_135_135_45km_150mstg_B0_HLLC", "text":[30]},
        {"name":"DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS", "text":[22,18]},
        {"name":"SFHoTim276_12_15_0025_150mstg_B15_HLLD_CT_GS_onFugaku","text":[28]},
        {"name":"SFHoTim276_13_14_0025_150mstg_B0_HLLC","text":[22, 28, 30, 26, 20, 24]},
        {"name":"SFHoTim276_135_135_45km_150mstg_B0_FUKA","text":[32]}
    ]
    # conda activate PyBlastAfterglow
    failed = []
    for task in tasks:
        sim = task["name"]
        texts = task["text"]
        for text in texts:
            for method in ["rf","xgboost"]:
                root_path = os.getcwd() + '/' + sim + '/'
                main_run(
                    fpath_to_collated=root_path+"collated.csv",
                    path_to_save_run=root_path,
                    target_key="flux",
                    # model_name="xgboost",
                    model_name=method,
                    text=text,
                    sample=None,#=150*10
                    n_jobs=32,
                    study_name=sim + '_' + f"study_{'flux'}_{method}_text{text}",
                )

    print(f"Completed. Failed: {len(failed)}")
    print(failed)

    #
    #
    # root_path = os.getcwd() + '/' + sim + '/'
    # main_run(
    #     fpath_to_collated=root_path+"collated.csv",
    #     path_to_save_run=root_path,
    #     target_key="flux",
    #     # model_name="xgboost",
    #     model_name="rf",
    #     text=32,
    #     sample=None,#=150*10
    #     n_jobs=32
    # )

