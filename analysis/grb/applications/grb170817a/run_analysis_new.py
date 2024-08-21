# import PyBlastAfterglowMag
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import time as timer

from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
# from skimage.metrics import mean_squared_error

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs
from sklearn.metrics import root_mean_squared_error

# temprary dir to save PyBlastAfterglow output into
working_dir = os.getcwd() + '/working_dirs/'
fig_dir = os.getcwd() + '/figs/'

def grb170817_vla_hubble_skymap_data():
    """
    Make a dataset of observed flux centroid positions for 8, 75, 206 and 230 days after merger.
    Here both, VLA and HST, HSA data is given.
    Radio data is expected to be given at 4.5e9 Hz; while optical data at 8 days is presumably at 4.6e14 Hz
    Resulting dataset has the following structure

    # time[days] freq[Hz] xc[mas] yc[mas] xc_err_left xc_err_right yc_err_low yc_err_up
    # ...        ...      ...     ...     ...         ...          ...        ...
    :return:
    """
    # extracted from https://www.nature.com/articles/s41586-022-05145-7 Fig.1
    data = np.array([
        [6.92063492063492,	1.68247422680412],
        [6.91609977324263,	1.28041237113402],
        [6.91609977324263,	2.07835051546392],
        [7.08390022675737,	1.68247422680412],
        [6.75283446712018,	1.68247422680412],
        [5.83673469387755,	1.52783505154639],
        [5.83673469387755,	1.29278350515464],
        [5.84126984126984,	1.77525773195876],
        [6.04988662131519,	1.5340206185567],
        [5.63265306122449,	1.54020618556701],
        [3.99546485260771,	1.50309278350515],
        [4.0,	1.10103092783505],
        [4.0,	1.90515463917526],
        [4.11337868480726,	1.50309278350515],
        [3.87755102040816,	1.50309278350515],
        [1.37868480725624,	1.58350515463918],
        [1.38321995464853,	1.40412371134021],
        [1.38321995464853,	1.76907216494845],
        [1.70068027210884,	1.58969072164948],
        [1.06575963718821,	1.58969072164948],
    ])
    times = np.array([8., 75., 206., 230.])
    freqs = np.array([4.6e14, 4.5e9, 4.5e9, 4.5e9])

    xs = data[0::5,0]
    left = data[3::5,0]
    right = data[4::5,0]

    ys = data[0::5,1]
    low = data[1::5,1]
    up = data[2::5,1]

    # print('x', xs)
    # print('left', left)
    # print('right', right)
    #
    # print('ys',ys)
    # print('lower',low)
    # print('upper',up)
    # print('---------------------')

    # noramize
    low-=ys[2]; up-=ys[2]; ys-=ys[2]# set so that data at 75 days is at [0,0]
    left-=xs[2]; right-=xs[2];xs-=xs[2] # set so that data at 75 days is at [0,0]

    # print('x', xs)
    # print('left', left)
    # print('right', right)
    #
    # print('ys',ys)
    # print('lower',low)
    # print('upper',up)

    left = left - xs
    right = xs - right

    low = ys - low
    up = up - ys

    res = np.column_stack(
        (times, freqs, xs, ys, left, right, low, up)
    )
    return (res)

def grb170817_radio_xray_wang23_data(freq:float,mult:float,data_dir:str):
    if (freq == 5.06e14 or freq == "optical") and mult == 100:
        return np.loadtxt(data_dir+"grb_170817_Optical_x100.txt")
    elif (freq==2.41e17 or freq=="1KeV" or freq=="xray") and mult == 800:
        return np.loadtxt(data_dir+"grb_170817_1KeV_x800.txt")
    elif (freq==3e9 or freq=="3GHz") and mult == 6:
        return np.loadtxt(data_dir+"grb_170817_3GHz_x6.txt")
    elif (freq==6e9 or freq=="6GHz") and mult == 1:
        return np.loadtxt(data_dir+"grb_170817_6GHz.txt")
    else:
        raise FileNotFoundError(f"No data for freq={freq} and mult={mult}")

def compare_lcs_and_skymap_props(pba:PBA.PyBlastAfterglow,plot:dict,skymap_freq = 4.5e9):

    # ---------
    # ---------

    freqs = ()



    freqs = (3e9, 6e9, 5.06e+14, 2.41e+17)
    colors=("red","orange","green","blue")
    mults = (6, 1, 100, 800)
    lbls = (r"3\,GHz ($\times6$)", r"6\,GHz", r"optical ($\times100$)", r"xray ($\times800$)")

    fig, axes = plt.subplots(nrows=4,ncols=1,figsize=(5,7),layout='constrained')

    ''' -------------- light curves ------------- '''

    i=0
    errors = dict()
    errors_lc = [[],[]]
    for freq, color, mult, lbl in zip(freqs, colors, mults, lbls):
        lc = pba.GRB.get_lc(freq=freq)
        times = pba.GRB.get_lc_times()
        data = grb170817_radio_xray_wang23_data(freq=freq,mult=mult, data_dir=os.getcwd()+'/observations/')
        mask = times/cgs.day > data[:,0].min()
        axes[i].plot(times[mask]/cgs.day, lc[mask]*mult,color=color,ls='-')#, label="${}$".format(PBA.utils.latex_float(freq)))
        axes[i].errorbar(data[:,0], data[:,1], yerr=data[:,2:].T,
                         ecolor= 'gray', mec= color,  #uplims=uplims,
                         marker= 'o', markersize= 8,
                         capsize= 2, linestyle= "None", mfc= "white",
                         label= lbl, zorder= 100)
        lc_intep = interp1d(times/cgs.day, lc,bounds_error=False,fill_value=0)(data[:,0])
        errors[f"lc {lbl}"] = root_mean_squared_error(np.log10(data[:,1]),np.log10(lc_intep)) # log-residuals
        errors_lc[0].append( np.log10(data[:,1]))
        errors_lc[1].append( np.log10(lc_intep) )
    errors[f"lc"] = root_mean_squared_error(
        np.concatenate(errors_lc[0]).flatten(), np.concatenate(errors_lc[1]).flatten()
    )

    axes[i].legend(fancybox=False, loc='upper left', columnspacing=0.8,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol=1, fontsize=12, labelcolor='black',
                   framealpha=0.0, borderaxespad=0.)
    axes[i].set_yscale('log')
    # axes[i].set_xlabel("time [day]")
    # axes[i].grid(ls=':')
    axes[i].set_ylabel(r"$F_{\nu}$ [mJy]",fontsize=12)
    # ax.legend(fancybox=False, loc='lower right', columnspacing=0.8,
    #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #           shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #           framealpha=0.0, borderaxespad=0.)
    # axes[i].set_xlim(1,2e3)
    axes[i].set_ylim(2e-5,2e0)
    # axes[i].set_xlim(5,1e3)
    # plt.show()

    ''' ----------- sky maps xc ------------- '''
    # data from https://arxiv.org/abs/2210.06568 (fig 2)
    skymap_times = (74.5719, 205.4533, 229.2961)
    skymap_xcs = (2.4054, 4.089, 5.0562)
    skymap_xc_errs = np.array([2.4054-2.0239, 2.7802-2.4054, 4.089-3.6593, 4.4974-4.089, 5.0562-4.6610, 5.4719-5.0562])

    times = pba.GRB.get_skymap_times()
    skymaps = [pba.GRB.get_skymap(time=t, freq=skymap_freq) for t in times]
    xcs = [skymap.xc for skymap in skymaps]
    # delta_x = [abs(skymap.x2-skymap.x1) for skymap in skymaps]
    # fluxes = [skymap.flux for skymap in skymaps]

    xcs_interp = interp1d(times/cgs.day, np.array(xcs), fill_value=0)(np.array(skymap_times))
    errors["xc"] = root_mean_squared_error(np.array(skymap_xcs), np.array(xcs_interp))

    # fig,axes = plt.subplots(ncols=1,nrows=2,sharex='all',layout='constrained',figsize=(6,4))
    i+=1
    axes[i].plot(times/cgs.day, xcs, color='black', ls='-', marker='.')
    axes[i].set_ylabel(r"$x_c$ [mas]",fontsize=12)
    axes[i].set_ylim(0, 7)
    # axes[i].grid(ls=':')
    # axes[i].set_xlim(5,1e3)
    # axes[i].set_xscale('log')


    axes[i].errorbar(skymap_times,
                     skymap_xcs,
                     # xerr=np.array([2, 3, 2, 3, 2, 3]).reshape(-1,3),
                     yerr=skymap_xc_errs.reshape(-1,3),
                     **{'ecolor': 'gray', 'mec': 'gray', 'marker': 'o', 'markersize': 8, "capsize": 2,
                        "linestyle": "None", 'mfc': "white", "label": r"GRB\,170817A", "zorder": 100})
    if plot["show_subluminal"]:
        sublum_dxcs = np.zeros_like(times)
        for j in range(len(times)):
            beta_app_sub = 1.
            sublum_dxcs[j] = beta_app_sub / cgs.rad2mas * times[j]
        sublum_xcs = np.cumsum(sublum_dxcs)#[:-1] + sublum_dxcs[1:]
        axes[i].plot(times[1:]/cgs.day, sublum_xcs[1:], color='gray', ls=':')
        axes[i].fill_between(times[1:]/cgs.day, sublum_xcs[1:], where=sublum_xcs[1:]>=0, interpolate=True, color='gray',
                             label=r"subluminal")
    axes[i].legend(**{"fancybox": False, "loc": 'upper left', "shadow": "False", "ncol": 1, "fontsize": 12,
                      "framealpha": 0., "borderaxespad": 0., "frameon": False})

    ''' ----------- sky maps beta_app ------------- '''

    # data from
    obs_times_beta_app = [(75-8)/2, (206-8)/2, (230-8)/2]
    obs_beta_app = [7.6, 4.7, 5.2]
    obs_beta_app_time_err = np.array([(75-8)/2-8, 75-(75-8)/2, (206-8)/2-8, 206-(206-8)/2, (230-8)/2-8, 230-(230-8)/2])
    obs_beta_app_beta_err = np.array([1.3,1.3, 0.6,0.6, 0.5,0.5])
    beta_app = np.diff(xcs)*cgs.rad2mas / times[1:] # / cgs.c

    # average for smoothness
    N = 3
    beta_app_averaged = np.convolve(beta_app, np.ones(N)/N, mode='valid')
    time_for_beta_app_averaged =  np.convolve(times[1:]/cgs.day, np.ones(N)/N, mode='valid')

    beta_app_averaged_interp = interp1d(time_for_beta_app_averaged, beta_app_averaged, fill_value=0)(np.array(obs_times_beta_app))
    errors["beta_app"] = root_mean_squared_error(np.array(obs_beta_app), beta_app_averaged_interp)

    i+=1
    axes[i].plot(time_for_beta_app_averaged, beta_app_averaged, color='black', ls='-', marker='.')
    axes[i].set_ylabel(r"$\beta_{\rm app}$ [c]",fontsize=12)
    axes[i].set_ylim(-1,8)
    # axes[i].grid(ls=':')
    # axes[i].set_xlim(1,1e4)
    # axes[i].set_xscale('log')
    #4.1±0.5 # Mooley et al 2018
    axes[i].errorbar(obs_times_beta_app,obs_beta_app,
                     xerr=obs_beta_app_time_err.reshape(-1,3),
                     yerr=obs_beta_app_beta_err.reshape(-1,3),
                     **{'ecolor': 'gray', 'mec': 'gray', 'marker': 'o', 'markersize': 8, "capsize": 2,
                        "linestyle": "None", 'mfc': "white", "label": "GRB\, 170817A", "zorder": 100})

    axes[i].legend(**{"fancybox": False, "loc": 'lower left', "shadow": "False", "ncol": 1, "fontsize": 12,
                      "framealpha": 0., "borderaxespad": 0., "frameon": False})
    for ax in axes:
        ax.set_xlim(1,2e3)
        ax.grid(ls=':')
        ax.set_xscale('log')
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in',
                       bottom=True, top=True, left=True, right=True)
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([])
    axes[2].set_xlabel(r"$t_{\rm obs}$ [day]",fontsize=12)
    # plt.show()


    # -----------



    # fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(5,2),layout='constrained')

    ''' ----------- sky maps ------------- '''

    obs = grb170817_vla_hubble_skymap_data()
    obs_times = obs[:,0]
    obs_freqs = obs[:,1]
    obs_xcs = obs[:,2]

    ref_skymap = pba.GRB.get_skymap(time=obs_times[1]*cgs.day, freq=obs_freqs[1])
    offset = ref_skymap.xc

    ax =axes[3]
    for t,nu,obs_xc in zip(obs_times[::-1],obs_freqs[::-1],obs_xcs):
        if not nu in pba.GRB.get_skymap_freqs(unique=True):
            print(f"ERROR freq={nu} is not in the pba.skymaps={pba.GRB.get_skymap_freqs(unique=True)}")
            continue
        skymap=pba.GRB.get_skymap(time=t*cgs.day, freq=nu)
        int_zz = skymap.im_hist; int_x=skymap.grid_x-offset; int_y=skymap.grid_y
        xc = skymap.xc-offset
        yc = skymap.yc
        # obs_xc = obs_xc-obs_xcs[1]
        for i,x in enumerate(int_x):
            for j,y in enumerate(int_y):
                if (abs(x - xc) > 1. or abs(y - yc) > 1.):
                    int_zz[i,j] = np.nan

        norm = Normalize(int_zz[np.isfinite(int_zz)].max()*1e-2, int_zz[np.isfinite(int_zz)].max())
        cmap = cm.get_cmap('Greys')
        cmap.set_under('white')
        cmap.set_over('black')
        im = ax.pcolormesh(int_x, int_y, int_zz.T, norm=norm, cmap=cmap, alpha=0.7)
        # ax.set_facecolor(cmap(float(0.)))
        ax.plot(xc, yc, **dict(marker='s',color='black',markersize=12,fillstyle='none'))
        errors[f"xc {t}"] = xc - obs_xc
        # ax.errorbar([skymap.xc, skymap.xc], [skymap.y1, skymap.y2],
        #             xerr=[0.1, 0.1], **dict(color='black',ls="-"))
        # ax.errorbar([skymap.x1, skymap.x2], [skymap.yc, skymap.yc],
        #             yerr=[0.1, 0.1], **dict(color='black',ls="-"))
    for t,x in zip(obs_times,[0.85,0.5,0.3,0.15]):
        bbox = dict(boxstyle='round', fc='white', ec='black', alpha=1.)
        ax.text(x, 0.80, f"{int(t)}", fontsize=12, bbox=bbox,
                transform=ax.transAxes, horizontalalignment='right')

    # observational data of grb170817a (data from https://arxiv.org/abs/2210.06568, fig 1)
    obs = grb170817_vla_hubble_skymap_data()
    ax.errorbar(obs[:,2],obs[:,3],
                xerr=np.column_stack((obs[:,4],obs[:,5])).T,
                yerr=np.column_stack((obs[:,6],obs[:,7])).T,
                **{'ecolor': 'gray', 'mec': 'gray', 'marker': 'o', 'markersize': 8, "capsize": 2,
                   "linestyle": "None", 'mfc': "white", "label": "HST, HSA and gVLBI", "zorder": 100})



    # im.set_rasteraized(True)
    print(errors)
    with open(plot['outpath']+'summary.txt', 'w') as f:
        f.write('Error metrics for the run (log-RMSE, skymap errors, beta_errs, skymap position errors)\n')
        for key, val in errors.items():
            f.write(f"{key} = {val:.2f} \n")

    ax.set_yscale("linear")
    ax.set_xscale("linear")
    if "xlim" in plot.keys() and len(plot["xlim"]) == 2: ax.set_xlim(plot["xlim"])
    else: ax.set_xlim(skymap.grid_x.min() * 1.1, skymap.grid_x.max() * 1.1)
    if "ylim" in plot.keys() and len(plot["ylim"]) == 2: ax.set_ylim(plot["ylim"])
    else: ax.set_ylim(skymap.grid_y.min() * 1.1, skymap.grid_y.max() * 1.1)
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=12, direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.minorticks_on()
    ax.axhline(y=0, linestyle=':', linewidth=0.4)
    ax.axvline(x=0, linestyle=':', linewidth=0.4)
    ax.grid(ls=':')
    ax.set_xlabel("$X$ [mas]",fontsize=12)
    ax.set_ylabel("$Z$ [mas]",fontsize=12)
    ax.set_aspect(aspect=1/2)
    ax.legend(fancybox=False, loc='lower right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=1, fontsize=12, labelcolor='black',
              framealpha=0.0, borderaxespad=0.)
    plt.subplots_adjust(wspace=0.01)
    if "figname" in plot.keys():
        plt.savefig(plot['outpath']+plot["figname"]+'.png',dpi=256)
    plt.savefig(plot['outpath']+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()

if __name__ == '__main__':
    # https://arxiv.org/pdf/1808.00469
    # https://www.nature.com/articles/s41586-022-05145-7 (https://arxiv.org/abs/2210.06568) (figure 1,2 data)
    struct = dict(struct="gaussian",Eiso_c=np.power(10, 54.1), Gamma0c= 300., M0c= -1.,
                  theta_c= np.deg2rad(3.50), theta_w= np.deg2rad(25.))
    skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
                     intp_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'),  # "gaussian"
                     hist_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'))


    working_dir = working_dir + 'tmp/'
    start_time = timer.perf_counter()
    pba = PBA.wrappers.run_grb(
        working_dir=working_dir,
        P=dict(main=dict(lc_freqs='array logspace 1e8 1e29 96',
                         lc_times='array logspace 3e3 1e10 128',
                         tb0=3e2,tb1=1e14,ntb=1000,
                         skymap_freqs="array 4.5e9 4.6e14", # Radio of VLA and Optical of HST
                         skymap_times="array logspace 1e3 1e8 100",
                         n_ism = np.power(10, -1.60),theta_obs=np.deg2rad(20.8),d_l = 1.27e+26, z = 0.0099),
               grb=dict(structure=struct,eats_type='a',do_rs='no',bw_type='fs',do_lc='yes',do_skymap='yes',
                        eps_e_fs = np.power(10,-3.42), eps_b_fs =np.power(10,-4.02), p_fs = 2.10,
                        eps_e_rs = np.power(10,-3.42), eps_b_rs =np.power(10,-4.02), p_rs = 2.10,
                        method_ssc_fs='none',method_pp_fs='none',use_ssa_fs='no',
                        method_ssc_rs='none',method_pp_rs='none',use_ssa_rs='no',
                        ebl_tbl_fpath='none',
                        skymap_conf=skymap_conf)),
        run=False,
        process_skymaps=False,
        loglevel="err",#"err"
    )
    finish_time = timer.perf_counter()
    print(f"PyBlastAfterglow finished in {finish_time-start_time:.2f} seconds. ({(finish_time-start_time)/60.:.2f} min.) ")

    pba = PBA.PyBlastAfterglow(workingdir=working_dir)
    compare_lcs_and_skymap_props(pba=pba,
                                 plot = dict(xlim=(4, -4), ylim=(-2, 2), show=True, figname="lcs_170817a_gauss",
                                             outpath=working_dir,
                                             show_subluminal=False))

    # struct = dict(struct="gaussian",Eiso_c=np.power(10, 54.2), Gamma0c= 300., M0c= -1.,
    #               theta_c= np.deg2rad(4.0), theta_w= np.deg2rad(40.))
    # skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
    #                  intp_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'),  # "gaussian"
    #                  hist_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'))
    # compare_lcs_and_skymap_props(
    #     # struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
    #     pp = dict(main=dict(lc_freqs='array logspace 1e8 1e29 96',
    #                         lc_times='array logspace 3e3 1e10 128',
    #                         tb0=3e1,tb1=1e15,ntb=2000,iout=2,
    #                         skymap_freqs="array 4.5e9 4.6e14", # Radio of VLA and Optical of HST
    #                         skymap_times="array logspace 1e3 1e8 100",
    #                         n_ism = np.power(10, -1.95),theta_obs=np.deg2rad(20.),d_l = 1.27e+26, z = 0.0099),
    #               grb=dict(structure=struct,eats_type='a',do_rs='no',bw_type='fs',do_lc='yes',do_skymap='yes',
    #                        eps_e_fs = np.power(10,-3.42), eps_b_fs =np.power(10,-4.02), p_fs = 2.10,
    #                        eps_e_rs = np.power(10,-3.42), eps_b_rs =np.power(10,-4.02), p_rs = 2.10,
    #                        method_ssc_fs='none',method_pp_fs='none',use_ssa_fs='no',
    #                        method_ssc_rs='none',method_pp_rs='none',use_ssa_rs='no',
    #                        ebl_tbl_fpath='none',do_mphys_in_situ='no',do_mphys_in_ppr='yes',
    #                        skymap_conf=skymap_conf)),
    #     plot = dict(xlim=(4, -4), ylim=(-2, 2), show=True, figname="lcs_170817a_gauss"),
    #     working_dir=working_dir+"tmp/",
    #     run=True,process_skymaps=True
    # )