import shutil,json,os,h5py
import pandas as pd
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
from multiprocessing import Pool

from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm

# try:
#     import package.src.PyBlastAfterglowMag as PBA
# except ImportError:
#     try:
#         import PyBlastAfterglowMag as PBA
#     except:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

import package.src.PyBlastAfterglowMag as PBA

# --------------------------------------------------------------------------------
DATA_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/kenta_data/"
get_ej_data = lambda name : DATA_PATH+name+'/'+"ej_collated.h5"

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
# def piecewise_power(x, x0, x1, y0, k1, k2, k3):
#     condlist = [x < x0,
#                 (x >= x0) & (x < x1),
#                 x >= x1]
#     funclist = [lambda x: y0*(x/x0)**k1,# k1 * x + y0 - k1 * x0,
#                 lambda x: y0*(x/x0)**k2,#k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
#                 # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
#                 lambda x: y0*(x/x1)**k3 * (x1/x0)**k2 #k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
#                 ]
#     return np.piecewise(x, condlist, funclist)
def piecewise_power(x, x0, x1, x2, y0, k1, k2, k3):
    condlist = [x < x0,
                (x >= x0) & (x < x1),
                (x >= x1) & (x < x2),
                x >= x2]
    funclist = [lambda x: y0*(x/x0)**k1, # k1 * x + y0 - k1 * x0,
                lambda x: y0,
                lambda x: y0*(x/x1)**k2, #k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                lambda x: y0*(x/x2)**k3 * (x2/x1)**k2 #k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                ]
    return np.piecewise(x, condlist, funclist)


def runs(task:dict):
    pba = PBA.wrappers.run_kn(
        working_dir=task["working_dir"],
        struct=task["struct"],
        P=task["P"],
        run=task['run']
    )
    pba.clear()

def compare_nr_and_fit(run:bool,run_fit:bool,sim_:str or None,xlim:tuple or None,ylim0:tuple or None,ylim1:tuple or None,
                       figname:str, P:dict,fit_func:str,suffix:str):

    # do_cumulative = True
    # log_type = 2
    # log_type_y = 2

    # get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    # do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    # un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))
    #
    # do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    # un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))


    df_fit = pd.read_csv(EJ_TEXT_PATH+f"piecewise_line_{fit_func}.csv",index_col=0)

    '''
    tasks = []
    n_shells = {}
    fig, axes = plt.subplots(ncols=1,nrows=2,sharex='all')
    for sim, sim_dict in df.iterrows():
        if sim_ and not sim == sim_:
            continue
        fit_dict = df_fit.loc[sim]
        text_dict = df_text.loc[sim]
        data = PBA.id_kenta.EjStruct(fpath=get_ej_data(sim_dict['name']),verbose=True)
        vinf = data.get_vinf()
        masses = data.get(v_n="mass",text=float(sim_dict["tmerg"])+float(text_dict["text"]))
        sums = np.zeros(len(masses[0,:]))
        sums_ek = np.zeros(len(masses[0,:]))
        for ivinf in range(len(masses[0,:])):
            masses[:,ivinf] = np.sum(masses[:,ivinf])/len(masses[:,ivinf])
            sums[ivinf] = np.sum(masses[:,ivinf])
            sums_ek[ivinf] = np.sum(masses[:,ivinf]) * (vinf[ivinf]**2*PBA.cgs.c**2*PBA.utils.cgs.solar_m)
        # ek = sums * PBA.cgs.c**2 * vinf *
        # vinf = 0.5 * (vinf[1:]+vinf[:-1])
        vinf_ = np.linspace(0.01,1, 102)
        mom=PBA.MomFromBeta(vinf_)[:-1]
        # plt.plot(range(len(vinf)),vinf,marker='.')
        # plt.plot(range(len(vinf_)),vinf_,marker='x')
        # plt.show()
        n_shells[sim] = len(vinf)
        # mom = np.linspace(*mom_lim, len(vinf))
        coeffs = fit_dict[["x0", "x1", "x2", "y0", "k1", "k2", "k3"]]
        # for i in [0,1,2]: coeffs[i] = un_log(coeffs[i])


        ek = piecewise_power(mom, *coeffs)
        mass = ek / (PBA.BetaFromMom(mom)*PBA.utils.cgs.c)**2 / PBA.utils.cgs.solar_m

        # mass = piecewise_power(mom, *coeffs) / PBA.cgs.solar_m

        axes[0].plot(PBA.utils.MomFromBeta(vinf[vinf<1.]),sums[vinf<1.],color=sim_dict['color'],ls=sim_dict['ls'],label=sim_dict['label'],lw=2)
        axes[0].plot(mom,mass,color=sim_dict['color'],ls=sim_dict['ls'],lw=1.)

        # eks =

        axes[1].plot(PBA.utils.MomFromBeta(vinf[vinf<1.]),sums_ek[vinf<1.],color=sim_dict['color'],ls=sim_dict['ls'],label=sim_dict['label'],lw=2)
        axes[1].plot(mom,ek,color=sim_dict['color'],ls=sim_dict['ls'],lw=1.)

    # for ax in axes:
    #     ax.set_xscale('log')
        # ax.set_yscale('log')
    plt.legend()
    plt.show()
    '''
    tasks = []
    for sim, sim_dict in df.iterrows():
        if sim_ and not sim == sim_:
            continue

        fit_dict = df_fit.loc[sim]
        text_dict = df_text.loc[sim]

        # using NR data
        for force_spherical, mode in zip([False,True],["asym","sph"]):
            tasks.append(
                dict(working_dir=os.getcwd()+'/working_dirs/'+f'{sim}_{mode}_{suffix}/',
                     struct=dict(struct="numeric",
                                 n_layers_pw=30,
                                 corr_fpath_kenta=get_ej_data(sim_dict['name']),
                                 text=float(sim_dict["tmerg"])+float(text_dict["text"]),
                                 t0=1e3,
                                 force_spherical=force_spherical
                                 ),
                     P=P,
                     run=run,
                     label=f"{sim} {mode}"
                     )
            )

        # # using fitting function
        # coeffs = fit_dict[["x0", "x1", "x2", "y0", "k1", "k2", "k3"]]
        # # mom = np.linspace(*mom_lim, n_shells[sim])
        # vinf_ = np.linspace(0.01,1, n_shells[sim])
        # mom=PBA.MomFromBeta(vinf_)[:-1]
        # ek = piecewise_power(mom, *coeffs)
        # mass = ek / (PBA.cgs.solar_m * PBA.cgs.c**2 * PBA.BetaFromMom(mom)**2)#piecewise_power(mom, *coeffs) / PBA.cgs.solar_m

        _, _, l_mom, l_ek = np.loadtxt(EJ_TEXT_PATH+sim_dict['name']+f"_log_mom_log_ek_sph_and_fit_{fit_func}.txt",unpack=True)
        mom,ek=10**l_mom,10**l_ek
        mass = ek / (PBA.cgs.solar_m * PBA.cgs.c**2 * PBA.BetaFromMom(mom)**2)#piecewise_power(mom, *coeffs) / PBA.cgs.solar_m
        tasks.append(
            dict(
                working_dir=os.getcwd()+'/working_dirs/'+f'{sim}_{fit_func}_{suffix}/',
                struct=dict(struct="numeric",
                            dist='pw',
                            n_layers_pw=30,
                            hist_1d_data=dict(mom=mom, ek=ek, mass=mass),
                            t0=1e3,
                            ),
                P=P,
                run=run_fit,
                label=f"{sim} {'fit'}"
            )
        )

    # run all simulations asynch
    if run or run_fit:
        # equal load run
        nmax = 14
        if len(tasks) > nmax:
            nn = int(len(tasks) / nmax) + 1
            ncpu = len(tasks) // nn
        else:
            ncpu = len(tasks)
        with Pool(ncpu) as pool:
            result = pool.map(runs, tasks)

    # plot the result
    fig, axes = plt.subplots(ncols=1,nrows=2,figsize=(4.6,1.5*3.2),layout='constrained',sharex='col',#sharex='col',sharey='row',
                             gridspec_kw={'height_ratios': [2,1]})

    # pba_dict = dict()
    # for sim, sim_dict in df.iterrows():
    #     pba_dict[sim] = [task for task in tasks if task['label'].__contains__(sim)]

    fit_perform = {}
    for sim, sim_dict in df.iterrows():
        if sim_ and sim != sim_:
            continue

        task_asym = [task for task in tasks if task["label"]==f"{sim} {'asym'}"][0]
        pba_asym = PBA.wrappers.run_kn(working_dir=task_asym["working_dir"],struct=task_asym["struct"],P=task_asym["P"],run=False)

        task_sph = [task for task in tasks if task["label"].__contains__(f"{sim} {'sph'}")][0]
        pba_sph = PBA.wrappers.run_kn(working_dir=task_sph["working_dir"],struct=task_sph["struct"],P=task_sph["P"],run=False)

        task_fit = [task for task in tasks if task["label"].__contains__(f"{sim} {'fit'}")][0]
        pba_fit = PBA.wrappers.run_kn(working_dir=task_fit["working_dir"],struct=task_fit["struct"],P=task_fit["P"],run=False)

        axes[0].plot(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,   pba_asym.KN.get_lc(freq=3.e9) * 1e3 ,
                     color=sim_dict["color"],ls=sim_dict['ls'], lw=1.2, label=sim_dict['label'])
        # axes[0].plot(pba_sph.KN.get_lc_times() / PBA.utils.cgs.day,    pba_sph.KN.get_lc(freq=3.e9) * 1e3 ,
        #              color=sim_dict["color"],ls=sim_dict['ls'], lw=1.2)
        axes[0].plot(pba_fit.KN.get_lc_times() / PBA.utils.cgs.day,    pba_fit.KN.get_lc(freq=3.e9) * 1e3 ,
                     color=sim_dict["color"], ls=sim_dict['ls'], lw=0.6)
        print(pba_asym.KN.get_lc(freq=3.e9) * 1e3 )
        axes[0].fill_between(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,
                             pba_asym.KN.get_lc(freq=3.e9) * 1e3 ,
                             pba_sph.KN.get_lc(freq=3.e9) * 1e3,
                             color=sim_dict["color"],ls=sim_dict['ls'], alpha=0.5)

        axes[1].plot(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,
                     (np.log10(pba_sph.KN.get_lc(freq=3.e9)) - np.log10(pba_fit.KN.get_lc(freq=3.e9))),
                     color=sim_dict["color"],ls=sim_dict['ls'], lw=1.)

        fit_perform[sim] = dict(
            # Calculate R-squared
            r_squared = r2_score(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
                                 np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            # Calculate MSE and RMSE
            mse = mean_squared_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
                                     np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            rmse = np.sqrt(mean_squared_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
                                              np.log10(pba_fit.KN.get_lc(freq=3.e9)))),
            # Calculate MAE
            mae = mean_absolute_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
                                      np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            # sum of squared residuals
            sse = np.sum((np.log10(pba_sph.KN.get_lc(freq=3.e9))-np.log10(pba_fit.KN.get_lc(freq=3.e9)))**2)
        )

    df_fit_perform = pd.DataFrame.from_dict(fit_perform).T
    df_fit_perform.to_csv(os.getcwd()+f'/nr_{fit_func}_lcs_{suffix}.csv',index=True)
    print(df_fit_perform)

    axes[0].set_yscale("log")
    for ax in axes:
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.minorticks_on()
        ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        ax.minorticks_on()

    n = 2
    ax = axes[0]
    ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-', label='Simulation',lw=1.2)#, lw=0.7, drawstyle='steps')
    ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-',label='Fit',lw=0.6)#, lw=0.7, drawstyle='steps')
    han, lab = ax.get_legend_handles_labels()
    ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
                            **dict(fancybox=False,loc= 'lower center',columnspacing=0.4,
                                   #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                                   shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
                            **dict(fancybox=False,loc= 'upper left',columnspacing=0.4,
                                   #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                                   shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))

    if xlim: ax.set_xlim(*xlim)
    if ylim0: axes[0].set_ylim(*ylim0)
    if ylim1: axes[1].set_ylim(*ylim1)
    axes[0].set_ylabel(r"$F_{\nu}$ [$\mu$Jy]",fontsize=12)
    axes[1].set_ylabel(r"$\log_{10}(F_{\nu\,;\rm sph}) - \log_{10}(F_{\nu\,;\rm fit})$",fontsize=12)
    axes[1].axhline(y=0,color='gray',linestyle='-')
    axes[1].grid(color='gray',linestyle=':')
    axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]",fontsize=12)

    figname = os.getcwd()+f'/figs/{figname}_{suffix}_{fit_func}'
    plt.savefig(figname+'.png',dpi=256)
    plt.savefig(figname+'.pdf')
    plt.show()

    plt.show()
def plot_nr_and_fits(run:bool,run_fit:bool,sim_:str or None,xlim:tuple or None,ylim0:tuple or None,ylim1:tuple or None,
                       figname:str, P:dict,suffix:str,show_legends:bool):

    # do_cumulative = True
    # log_type = 2
    # log_type_y = 2

    # get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    # do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    # un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))
    #
    # do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    # un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))

    # fit_funcs = [
    #     dict(name="3segFit", label="3SegFit", lw=0.6),
    #     dict(name="4segFit", label="4SegFit", lw=0.3)
    # ]
    fit_funcs = ["3segFit","4segFit"]




    # df_fit- = pd.read_csv(EJ_TEXT_PATH+f"piecewise_line_{fit_func}.csv",index_col=0)

    '''
    tasks = []
    n_shells = {}
    fig, axes = plt.subplots(ncols=1,nrows=2,sharex='all')
    for sim, sim_dict in df.iterrows():
        if sim_ and not sim == sim_:
            continue
        fit_dict = df_fit.loc[sim]
        text_dict = df_text.loc[sim]
        data = PBA.id_kenta.EjStruct(fpath=get_ej_data(sim_dict['name']),verbose=True)
        vinf = data.get_vinf()
        masses = data.get(v_n="mass",text=float(sim_dict["tmerg"])+float(text_dict["text"]))
        sums = np.zeros(len(masses[0,:]))
        sums_ek = np.zeros(len(masses[0,:]))
        for ivinf in range(len(masses[0,:])):
            masses[:,ivinf] = np.sum(masses[:,ivinf])/len(masses[:,ivinf])
            sums[ivinf] = np.sum(masses[:,ivinf])
            sums_ek[ivinf] = np.sum(masses[:,ivinf]) * (vinf[ivinf]**2*PBA.cgs.c**2*PBA.utils.cgs.solar_m)
        # ek = sums * PBA.cgs.c**2 * vinf *
        # vinf = 0.5 * (vinf[1:]+vinf[:-1])
        vinf_ = np.linspace(0.01,1, 102)
        mom=PBA.MomFromBeta(vinf_)[:-1]
        # plt.plot(range(len(vinf)),vinf,marker='.')
        # plt.plot(range(len(vinf_)),vinf_,marker='x')
        # plt.show()
        n_shells[sim] = len(vinf)
        # mom = np.linspace(*mom_lim, len(vinf))
        coeffs = fit_dict[["x0", "x1", "x2", "y0", "k1", "k2", "k3"]]
        # for i in [0,1,2]: coeffs[i] = un_log(coeffs[i])


        ek = piecewise_power(mom, *coeffs)
        mass = ek / (PBA.BetaFromMom(mom)*PBA.utils.cgs.c)**2 / PBA.utils.cgs.solar_m

        # mass = piecewise_power(mom, *coeffs) / PBA.cgs.solar_m

        axes[0].plot(PBA.utils.MomFromBeta(vinf[vinf<1.]),sums[vinf<1.],color=sim_dict['color'],ls=sim_dict['ls'],label=sim_dict['label'],lw=2)
        axes[0].plot(mom,mass,color=sim_dict['color'],ls=sim_dict['ls'],lw=1.)

        # eks =

        axes[1].plot(PBA.utils.MomFromBeta(vinf[vinf<1.]),sums_ek[vinf<1.],color=sim_dict['color'],ls=sim_dict['ls'],label=sim_dict['label'],lw=2)
        axes[1].plot(mom,ek,color=sim_dict['color'],ls=sim_dict['ls'],lw=1.)

    # for ax in axes:
    #     ax.set_xscale('log')
        # ax.set_yscale('log')
    plt.legend()
    plt.show()
    '''
    tasks = []
    for sim, sim_dict in df.iterrows():
        if sim_ and not sim == sim_:
            continue

        # fit_dict = df_fit.loc[sim]
        text_dict = df_text.loc[sim]

        # using NR data
        for force_spherical, mode in zip([False,True],["asym","sph"]):
            tasks.append(
                dict(working_dir=os.getcwd()+'/working_dirs/'+f'{sim}_{mode}_{suffix}/',
                     struct=dict(struct="numeric", n_layers_pw=30,
                                 corr_fpath_kenta=get_ej_data(sim_dict['name']),
                                 text=float(sim_dict["tmerg"])+float(text_dict["text"]),
                                 t0=1e3,
                                 force_spherical=force_spherical
                                 ),
                     P=P,
                     run=run,
                     label=f"{sim} {mode}"
                     )
            )

        # # using fitting function
        # coeffs = fit_dict[["x0", "x1", "x2", "y0", "k1", "k2", "k3"]]
        # # mom = np.linspace(*mom_lim, n_shells[sim])
        # vinf_ = np.linspace(0.01,1, n_shells[sim])
        # mom=PBA.MomFromBeta(vinf_)[:-1]
        # ek = piecewise_power(mom, *coeffs)
        # mass = ek / (PBA.cgs.solar_m * PBA.cgs.c**2 * PBA.BetaFromMom(mom)**2)#piecewise_power(mom, *coeffs) / PBA.cgs.solar_m

        for fit_func in fit_funcs:
            _, _, l_mom, l_ek = np.loadtxt(EJ_TEXT_PATH+sim_dict['name']+f"_log_mom_log_ek_sph_and_fit_{fit_func}.txt",
                                           unpack=True)
            mom,ek=10**l_mom,10**l_ek
            mass = ek / (PBA.cgs.solar_m * PBA.cgs.c**2 * PBA.BetaFromMom(mom)**2)#piecewise_power(mom, *coeffs) / PBA.cgs.solar_m
            tasks.append(
                dict(
                    working_dir=os.getcwd()+'/working_dirs/'+f'{sim}_{fit_func}_{suffix}/',
                    struct=dict(struct="numeric",
                                dist='pw',
                                n_layers_pw=30,
                                hist_1d_data=dict(mom=mom, ek=ek, mass=mass),
                                t0=1e3,
                                ),
                    P=P,
                    run=run_fit,
                    label=f"{sim} {fit_func}"
                )
            )

    # run all simulations asynch
    if run or run_fit:
        # equal load run
        nmax = 14
        if len(tasks) > nmax:
            nn = int(len(tasks) / nmax) + 1
            ncpu = len(tasks) // nn
        else:
            ncpu = len(tasks)
        with Pool(ncpu) as pool:
            result = pool.map(runs, tasks)

    # plot the result
    fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(4.6,2*3.2),layout='constrained',sharex='col',#sharex='col',sharey='row',
                             gridspec_kw={'height_ratios': [2,1,1]})

    # pba_dict = dict()
    # for sim, sim_dict in df.iterrows():
    #     pba_dict[sim] = [task for task in tasks if task['label'].__contains__(sim)]

    lc = lambda pba : pba.KN.get_lc(freq=3.e9) * 1e3
    loglc = lambda pba : np.log10( pba.KN.get_lc(freq=3.e9) * 1e3 )
    t = lambda pba : pba_asym.KN.get_lc_times() / PBA.utils.cgs.day

    fit_perform = {}
    for sim, sim_dict in df.iterrows():
        if sim_ and sim != sim_:
            continue

        task_asym = [task for task in tasks if task["label"]==f"{sim} {'asym'}"][0]
        pba_asym = PBA.wrappers.run_kn(working_dir=task_asym["working_dir"],struct=task_asym["struct"],P=task_asym["P"],run=False)

        task_sph = [task for task in tasks if task["label"].__contains__(f"{sim} {'sph'}")][0]
        pba_sph = PBA.wrappers.run_kn(working_dir=task_sph["working_dir"],struct=task_sph["struct"],P=task_sph["P"],run=False)

        axes[0].plot(t(pba_asym), lc(pba_asym), color=sim_dict["color"],ls=sim_dict['ls'], lw=1.2, label=sim_dict['label'])
        axes[0].fill_between(t(pba_asym), lc(pba_asym), lc(pba_sph), color=sim_dict["color"],ls=sim_dict['ls'], alpha=0.5)

        for (fit_func,lw,ax_fit) in zip(fit_funcs,[0.7,0.3],[axes[1],axes[2]]):
            task_fit = [task for task in tasks if task["label"].__contains__(f"{sim} {fit_func}")][0]
            pba_fit = PBA.wrappers.run_kn(working_dir=task_fit["working_dir"],struct=task_fit["struct"],P=task_fit["P"],run=False)
            axes[0].plot(t(pba_fit), lc(pba_fit), color=sim_dict["color"], ls=sim_dict['ls'], lw=lw)

            ax_fit.plot(t(pba_asym), (loglc(pba_sph) - loglc(pba_fit)),  color=sim_dict["color"], ls=sim_dict['ls'], lw=1)



        # task_fit = [task for task in tasks if task["label"].__contains__(f"{sim} {'fit'}")][0]
        # pba_fit = PBA.wrappers.run_kn(working_dir=task_fit["working_dir"],struct=task_fit["struct"],P=task_fit["P"],run=False)
        #
        # task_fit = [task for task in tasks if task["label"].__contains__(f"{sim} {'fit'}")][0]
        # pba_fit = PBA.wrappers.run_kn(working_dir=task_fit["working_dir"],struct=task_fit["struct"],P=task_fit["P"],run=False)
        #
        #
        # axes[0].plot(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,   pba_asym.KN.get_lc(freq=3.e9) * 1e3 ,
        #              color=sim_dict["color"],ls=sim_dict['ls'], lw=1.2, label=sim_dict['label'])
        # # axes[0].plot(pba_sph.KN.get_lc_times() / PBA.utils.cgs.day,    pba_sph.KN.get_lc(freq=3.e9) * 1e3 ,
        # #              color=sim_dict["color"],ls=sim_dict['ls'], lw=1.2)
        # axes[0].plot(pba_fit.KN.get_lc_times() / PBA.utils.cgs.day,    pba_fit.KN.get_lc(freq=3.e9) * 1e3 ,
        #              color=sim_dict["color"], ls=sim_dict['ls'], lw=0.6)
        # print(pba_asym.KN.get_lc(freq=3.e9) * 1e3 )
        # axes[0].fill_between(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,
        #                      pba_asym.KN.get_lc(freq=3.e9) * 1e3 ,
        #                      pba_sph.KN.get_lc(freq=3.e9) * 1e3,
        #                      color=sim_dict["color"],ls=sim_dict['ls'], alpha=0.5)
        #
        # axes[1].plot(pba_asym.KN.get_lc_times() / PBA.utils.cgs.day,
        #              (np.log10(pba_sph.KN.get_lc(freq=3.e9)) - np.log10(pba_fit.KN.get_lc(freq=3.e9))),
        #              color=sim_dict["color"],ls=sim_dict['ls'], lw=1.)

            # fit_perform[sim] = dict(
            #     # Calculate R-squared
            #     r_squared = r2_score(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
            #                          np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            #     # Calculate MSE and RMSE
            #     mse = mean_squared_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
            #                              np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            #     rmse = np.sqrt(mean_squared_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
            #                                       np.log10(pba_fit.KN.get_lc(freq=3.e9)))),
            #     # Calculate MAE
            #     mae = mean_absolute_error(np.log10(pba_sph.KN.get_lc(freq=3.e9)),
            #                               np.log10(pba_fit.KN.get_lc(freq=3.e9))),
            #     # sum of squared residuals
            #     sse = np.sum((np.log10(pba_sph.KN.get_lc(freq=3.e9))-np.log10(pba_fit.KN.get_lc(freq=3.e9)))**2)
            # )

    # df_fit_perform = pd.DataFrame.from_dict(fit_perform).T
    # df_fit_perform.to_csv(os.getcwd()+f'/nr_{fit_func}_lcs_{suffix}.csv',index=True)
    # print(df_fit_perform)

    axes[0].set_yscale("log")
    for ax in axes:
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.minorticks_on()
        ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        ax.minorticks_on()

    n = 3
    ax = axes[0]
    ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-', label='Simulation',lw=1.2)#, lw=0.7, drawstyle='steps')
    ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-',label='3SegFit',lw=0.4)#, lw=0.7, drawstyle='steps')
    ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-',label='4SegFit',lw=0.3)#, lw=0.7, drawstyle='steps')
    if show_legends:
        han, lab = ax.get_legend_handles_labels()
        ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
                                **dict(fancybox=False,loc= 'lower center',columnspacing=0.4,
                                       #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                                       shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
        ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
                                **dict(fancybox=False,loc= 'upper left',columnspacing=0.4,
                                       #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                                       shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))

    if xlim: ax.set_xlim(*xlim)
    if ylim0: axes[0].set_ylim(*ylim0)
    if ylim1: axes[1].set_ylim(*ylim1)
    if ylim1: axes[2].set_ylim(*ylim1)
    axes[0].set_ylabel(r"$F_{\nu}$ [$\mu$Jy]",fontsize=12)
    axes[1].set_ylabel(r"$\Delta\log_{10}F_{\nu;\,\rm sph,\,3SegFit}$",fontsize=12)#(r"$\log_{10}(F_{\nu\,;\rm sph}) - \log_{10}(F_{\nu\,;\rm fit})$"
    axes[2].set_ylabel(r"$\Delta\log_{10}F_{\nu;\,\rm sph,\,4SegFit}$",fontsize=12)#(r"$\log_{10}(F_{\nu\,;\rm sph}) - \log_{10}(F_{\nu\,;\rm fit})$",fontsize=12)
    axes[1].axhline(y=0,color='gray',linestyle='-')
    axes[1].grid(color='gray',linestyle=':')
    axes[2].grid(color='gray',linestyle=':')
    axes[-1].set_xlabel(r"$t_{\rm obs}$ [day]",fontsize=12)

    figname = os.getcwd()+f'/figs/{figname}_{suffix}_{"both_fits"}'
    plt.savefig(figname+'.png',dpi=256)
    plt.savefig(figname+'.pdf')
    plt.show()

    plt.show()

def main(sim_dict:pd.Series, coeffs_dict:pd.Series, ej_data_dict:pd.Series):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))

    data = PBA.id_kenta.EjStruct(fpath=get_ej_data(sim_dict['name']),verbose=True)
    vinf = data.get_vinf()
    masses = data.get(v_n="mass",text=float(sim_dict["tmerg"])+float(ej_data_dict["text"]))
    sums = np.zeros(len(masses[0,:]))
    for ivinf in range(len(masses[0,:])):
        masses[:,ivinf] = np.sum(masses[:,ivinf])/len(masses[:,ivinf])
        sums[ivinf] = np.sum(masses[:,ivinf])
    # ek = sums * PBA.cgs.c**2 * vinf *

    mom = np.logspace(-2.1, np.log10(4,), 100)
    coeffs = coeffs_dict[["x0", "x1", "y0", "k1", "k2", "k3"]]
    ek = piecewise_power(mom, *coeffs)
    mass = ek / (PBA.BetaFromMom(mom)*PBA.utils.cgs.c)**2 / PBA.utils.cgs.solar_m

    # plt.loglog(PBA.utils.MomFromBeta(vinf[vinf<1.]),sums[vinf<1.])
    # plt.loglog(mom,mass)
    # plt.show()
    P = dict(main=dict(n_ism=1e-1,d_l=100e6 * PBA.utils.cgs.pc, z=0.001,tb0 = 1e4, tb1 = 1e12,
                         theta_obs=np.pi/4.,integrator="DOP853",
                         lc_freqs = "array 3e9 1e18",
                         lc_times = "array logspace 1e5 1e9 100"),
               kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,
                       method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                       use_1d_id="no",do_skymap="no", do_lc="yes",
                       ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                       method_ne_fs="usenprime",method_comp_mode="observFlux",
                       method_gamma_min_fs='useU_e'
                       ))
    tasks = []
    for force_spherical, working_dir in zip([True,False],["sph","asym"]):
        tasks.append(
            dict(working_dir=os.getcwd()+'/'+f'working_dir_{working_dir}/',
                 struct=dict(struct="numeric",
                             n_layers_pw=30,
                             corr_fpath_kenta=get_ej_data(sim_dict['name']),
                             text=float(sim_dict["tmerg"])+float(ej_data_dict["text"]),
                             t0=1e3,
                             force_spherical=force_spherical
                             ),
                 P=P
            )
        )
    tasks.append(
        dict(
            working_dir=os.getcwd()+'/'+f'working_dir_{"fit"}/',
            struct=dict(struct="numeric",
                        dist='pw',
                        n_layers_pw=30,
                        hist_1d_data=dict(mom=mom, ek=ek),
                        t0=1e3,
                        ),
            P=P
        )
    )

    ncpu = min(len(tasks),14)
    with Pool(ncpu) as pool:
        result = pool.map(runs, tasks)


    for task, ls, label in zip(tasks,['-','--',':'],['NR',"Sph","fit"]):
        pba = PBA.wrappers.run_kn(working_dir=task["working_dir"],struct=task["struct"],P=task["P"],run=False)
        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
                color=sim_dict["color"],ls=ls,lw=1,label=label)
        nshells =  int( pba.KN.get_lc_obj().attrs["nshells"] )
        for ishell in np.arange(0,nshells,step=10)[:nshells]:
            ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9,ishell=ishell) * 1e3,
                    color='gray',ls=ls,lw=0.7)
        pba.clear()

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
    ax.set_xlim(1e-1, 1e5)
    ax.set_ylim(1e-1, 1e3)
    ax.legend()

    plt.show()

    for force_spherical, working_dir, ls in zip([True,False],["sph","asym"],["--",'-']):
        pba = PBA.wrappers.run_kn(
            working_dir=os.getcwd()+'/'+f'working_dir_{working_dir}/',
            struct=dict(struct="numeric",
                        n_layers_pw=30,
                        corr_fpath_kenta=get_ej_data(sim_dict['name']),
                        text=float(sim_dict["tmerg"])+float(ej_data_dict["text"]),
                        t0=1e3,
                        force_spherical=force_spherical
                        ),
            P=dict(main=dict(n_ism=1e-1,d_l=100e6 * PBA.utils.cgs.pc, z=0.001,tb0 = 1e4, tb1 = 1e12,
                             theta_obs=np.pi/4.,integrator="DOP853",
                             lc_freqs = "array 3e9 1e18",
                             lc_times = "array logspace 1e5 1e9 100"),
                   kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,
                           method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                           use_1d_id="no",do_skymap="no", do_lc="yes",
                           ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                           method_ne_fs="usenprime",method_comp_mode="observFlux",
                           method_gamma_min_fs='useU_e'
                           )),
            run=True)

        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
                color=sim_dict["color"],ls=ls,lw=1,label=working_dir)
        nshells =  int( pba.KN.get_lc_obj().attrs["nshells"] )
        for ishell in np.arange(0,nshells,step=10)[:nshells]:
            ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9,ishell=ishell) * 1e3,
                    color='gray',ls=ls,lw=0.7)
        pba.clear()

    mom = np.logspace(-2, np.log10(4,), nshells)
    coeffs = coeffs_dict[["x0", "x1", "y0", "k1", "k2", "k3"]]
    ek = piecewise_power(mom,*coeffs)
    # plt.loglog(mom,ek)
    # plt.show()

    pba = PBA.wrappers.run_kn(
        working_dir=os.getcwd()+'/'+f'working_dir_{"fit"}/',
        struct=dict(struct="numeric",
                    dist='pw',
                    n_layers_pw=30,
                    hist_1d_data=dict(mom=mom, ek=ek),
                    t0=1e3,
                    ),
        P=dict(main=dict(n_ism=1e-1,d_l=100e6 * PBA.utils.cgs.pc, z=0.001,tb0 = 1e4, tb1 = 1e12,
                         theta_obs=np.pi/4.,integrator="DOP853",
                         lc_freqs = "array 3e9 1e18",
                         lc_times = "array logspace 1e5 1e9 100"),
               kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,
                       method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                       use_1d_id="no",do_skymap="no", do_lc="yes",
                       ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                       method_ne_fs="usenprime",method_comp_mode="observFlux",
                       method_gamma_min_fs='useU_e'
                       )),
        run=True)

    ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
            color=sim_dict["color"],ls=':',lw=1,label='fit')
    for ishell in np.arange(0,nshells,step=10)[:nshells]:
        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9,ishell=ishell) * 1e3,
                color='gray',ls=':',lw=0.7)
    pba.clear()


    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
    ax.set_xlim(1e-1, 1e5)
    ax.set_ylim(1e-1, 1e3)
    ax.legend()

    plt.show()

def load_tables_print_table():
    # df_fit = pd.read_csv(os.getcwd()+f'/output/'+'nr_3segFit_lcs_joh06.csv',index_col=0)
    # df_fit = pd.read_csv(os.getcwd()+f'/output/'+'nr_4segFit_lcs_joh06.csv',index_col=0)
    # df_fit = pd.read_csv(os.getcwd()+f'/output/'+'nr_3segFit_lcs_marg21.csv',index_col=0)
    df_fit = pd.read_csv(os.getcwd()+f'/output/'+'nr_4segFit_lcs_marg21.csv',index_col=0)
    print(df_fit)
    # df_fit.rename(columns=dict(sse="3seg_joh06"),inplace=True)
    # df_fit=df_fit.filter(["3seg_joh06"])
    # df_fit.merge(pd.read_csv(os.getcwd()+'/'+'nr_4segFit_lcs_joh06.csv',index_col=0),right_index=True,left_index=True)
    # df_fit_sse = df_fit[["sse"]]
    # print(df_fit_sse)
    # print(df_fit)

# def analyze_2d_ejecta_data(sim_:str):
#     for sim, sim_dict in df.iterrows():
#         if sim_ and not sim == sim_:
#             continue


if __name__ == '__main__':
    # main(df.loc["SFHo_135_135_res150_new"],
    #      df_fit.loc["SFHo_135_135_res150_new"],
    #      df_text.loc["SFHo_135_135_res150_new"])

    P = dict(main=dict(n_ism=1e-1,d_l=100e6 * PBA.utils.cgs.pc, z=0.001, tb0 = 1e4, tb1 = 1e12,
                       theta_obs=np.pi/4.,integrator="DOP853",
                       lc_freqs = "array 3e9 1e18",
                       lc_times = "array logspace 1e5 1e9 100"),
             kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,eps_t_fs=1.,
                     method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                     use_1d_id="no",do_skymap="no", do_lc="yes",
                     ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                     method_ne_fs="usenprime",method_comp_mode="observFlux",
                     method_gamma_min_fs='useU_e',method_gamma_c_fs="useTe"
                     ))
    # compare_nr_and_fit(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.2,0.2),
    #                    fit_func = "3segFit",figname='radio_lcs_nr_fit',suffix="joh06")
    # compare_nr_and_fit(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.2,0.2),
    #                    fit_func = "4segFit",figname='radio_lcs_nr_fit',suffix="joh06")


    P = dict(main=dict(n_ism=0.1,d_l=100e6 * PBA.utils.cgs.pc, z=0.001, tb0 = 1e4, tb1 = 1e14,
                       theta_obs=np.pi/4.,integrator="DOP853",
                       lc_freqs = "array 3e9 1e18",
                       lc_times = "array logspace 1e5 1e11 100"),
             kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,eps_t_fs=1.,
                     method_synchrotron_fs="Marg21",method_ele_fs="analytic",
                     use_1d_id="no",do_skymap="no", do_lc="yes",
                     ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",#use_Sironi",
                     method_ne_fs="usenprime",method_comp_mode="observFlux",
                     method_gamma_min_fs='useTheta',method_B_fs='useMAG21'
                     ,method_gamma_c_fs="useTe"
                     ))
    # compare_nr_and_fit(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.3,0.3),
    #                    fit_func = "3segFit",figname='radio_lcs_nr_fit',suffix="marg21")
    # compare_nr_and_fit(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.3,0.3),
    #                    fit_func = "4segFit",figname='radio_lcs_nr_fit',suffix="marg21")

    # load_tables_print_table()

    # plot_nr_and_fits(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.2,0.2),
    #                  figname='radio_lcs_nr_fit',suffix="joh06",show_legends=True)
    # plot_nr_and_fits(P=P,run=False,run_fit=False,sim_=None,xlim=(1e0,1e4),ylim0=(7e-1,1e3),ylim1=(-0.2,0.2),
    #                  figname='radio_lcs_nr_fit',suffix="marg21",show_legends=False)

    # analyze_2d_ejecta_data()