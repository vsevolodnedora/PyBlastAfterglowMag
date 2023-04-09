# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import PyBlastAfterglow
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
                                                       BetFromMom, GamFromMom, MomFromGam)
    from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
    from PyBlastAfterglowMag.id_maker_from_thc_ourflow import prepare_kn_ej_id_2d
    from PyBlastAfterglowMag.skymap_tools import \
        (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
                                                           BetFromMom, GamFromMom, MomFromGam)
        from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
        from package.src.PyBlastAfterglowMag.id_maker_from_thc_ourflow import prepare_kn_ej_id_2d
        from package.src.PyBlastAfterglowMag.skymap_tools import \
            (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

def id_dyn_test():
    workdir = os.getcwd()+'/'

    dfile = h5py.File(workdir+"ejecta_2d_id_BLh_M13641364_M0_LK_SR.h5")
    dfile_ref = h5py.File(workdir+"ejecta_2d_id_BLh_M13641364_M0_LK_SR_ref.h5")

    # fig, ax = plt.subplots(ncols=1, nrows=1)
    # ax.plot(np.array(dfile["ctheta"]),BetFromMom(np.array(dfile["mom"])), 'x')
    # ax.plot(np.array(dfile_ref["theta"]),np.array(dfile_ref["vel_inf"]), '.')
    # ax.plot(np.arange(len(dfile["ctheta"])),np.array(dfile["ctheta"]),BetFromMom(np.array(dfile["mom"])), 'x')
    # ax.plot(np.arange(len(dfile["ctheta"])),np.array(dfile["ctheta"]),BetFromMom(np.array(dfile["mom"])), 'x')
    # ax.plot(np.array(dfile_ref["theta"]),np.array(dfile_ref["vel_inf"]), '.')
    # plt.show()

    # fig, ax = plt.subplots(ncols=1, nrows=1)
    # ax.plot(np.arange(len(dfile["ctheta"])),np.array(dfile["ctheta"]),marker='x',ls='none',color='red',label="our")
    # ax.plot(np.arange(len(dfile_ref["theta"])),np.array(dfile_ref["theta"]),marker='.',ls='none',color='red',label="ref")
    #
    # ax.plot(np.arange(len(dfile["mom"])),BetFromMom(np.array(dfile["mom"])),marker='x',ls='none',color='blue')
    # ax.plot(np.arange(len(dfile_ref["vel_inf"])),np.array(dfile_ref["vel_inf"]),marker='.',ls='none',color='blue')
    # plt.legend()
    # plt.show()

    fig, ax = plt.subplots(ncols=1, nrows=1)
    dfile = h5py.File(workdir+"dyn_kn.h5")
    fname_ref = h5py.File(workdir+"dyn_kn_ref.h5")
    nlayers = int(dfile.attrs["nlayers"])
    nshells = 20#int(dfile.attrs["nshells"])
    for il in range(nlayers-1):
        for ish in [nshells]:#range(nshells):
            key = f"shell={ish} layer={il}"
            # print(dfile[key].attrs["Gamma0"])
            ax.plot([dfile[key].attrs["Gamma0"]], [dfile[key].attrs["ctheta0"]],marker='x',label="our")
            ax.plot([fname_ref[key].attrs["Gamma0"]], [fname_ref[key].attrs["ctheta0"]],marker='.',label="ref")


            # ax.plot(np.array(dfile[key]["R"]),
            #         np.array(dfile[key]["Eint2"]))
            # ax.plot(np.array(fname_ref[key]["R"]),
            #         np.array(fname_ref[key]["Eint2"]), ls=':')
    # ax.axhline(y=np.pi/2.)

    # ax.set_xscale("log")
    # ax.set_yscale("log")
    plt.legend()

    plt.show()

def skymap_test():
    workdir = os.getcwd()+'/'
    dfile = h5py.File(workdir+"skymap_kn.h5")
    dfile_ref = h5py.File(workdir+"skymap_kn_ref.h5")
    print(dfile.keys())
    print(dfile_ref.keys())
    # exit(1)


def main():
    # skymap_test()
    # id_dyn_test()

    workdir = os.getcwd()+'/'
    # Pre-process the intput data
    fpath_2d_table = "corr_vel_inf_theta_BLh_M13641364_M0_LK_SR.h5"
    outfpath_2d_id = "ejecta_2d_id_BLh_M13641364_M0_LK_SR.h5"
    prepare_kn_ej_id_2d(nlayers=100,corr_fpath=workdir+fpath_2d_table,outfpath=workdir+outfpath_2d_id, dist="pw")

    # fpath_1d_table = "hist_vel_inf_BLh_M13641364_M0_LK_SR.dat"
    # outfpath_1d_id = "ejecta_1d_id_BLh_M13641364_M0_LK_SR.h5"
    # prepare_kn_ej_id_1d(nlayers=30, hist_fpath=workdir+fpath_1d_table, outfpath=workdir+outfpath_1d_id)

    pba = PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True,parfile="parfile.par")

    modify_parfile_par_opt(workingdir=workdir,part="kn",newpars={}, newopts={ "fname_ejecta_id":outfpath_2d_id},
                           parfile="parfile.par",newparfile="parfile.par",keep_old=False)

    pba.reload_parfile()
    pba.run()



    # plt.loglog(pba.KN.get_dyn_arr("R",ishell=29,ilayer=1),
    #            pba.KN.get_dyn_arr("ctheta",ishell=29,ilayer=1))
    # plt.show()

    workdir = os.getcwd() + '/'
    task_to_plot = {
        "workingdir" : workdir,
        "parfile1":  "parfile.par",
        "parfile2":  "parfile.par",
        "time": 120, "freq": 1e9, "title":"BLh $q=1.00$"
    }

    fname = workdir+ "precomputed_kn_skymap"#"int_"+sim + "_with_hists_time{}_".format(time) + prefix_th
    settings = {
        "gridspec": {
            "width_ratios": (4, 2), "height_ratios": (2, 4),
            "left": 0.14, "right": 0.95, "bottom": 0.1, "top": 0.96, "wspace": 0.05, "hspace": 0.05
        },
        "figname": fname,
        "paperpath": os.getcwd()+'/', "figfpath": fname, "save_pdf": True, "save_figs": True,
        "show_figs": True,
        "grid": False,
        "figsize": (4.8, 4.8),
        "rerun":False, "precompute": True,
        # "precompute_fpath": "/media/vsevolod/Data/postprocessed5/tmp/" + sim + "_single_time{}".format(time) + ".h5",
        "precompute_fpath": fname + ".h5",
        "kn_skymap": {
            "hist_nx": 256, "hist_ny": 128, "spec": False,
            "smooth": {"type": "gaussian", "sigma": 3},#"smooth": {"type": "uniform", "sigma": 5},
            "cm": {"color": 'yellow', "marker": "o"},
            "ysize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            "xsize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            # "pcolormesh": {"cmap": 'viridis', "norm": ("log", "0.2max", "1max"), "facecolor": 0.0, "alpha": 1.0},
            "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                           "norm": ("linear", "0.1max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
            # {}
        },
        "grb_skymap": {
            # "hist_nx": 125, "hist_ny": 75,
            # "cm": {"color": 'cyan', "marker": "+"},
            # "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            # "pcolormesh":{"cmap":'viridis', "norm":("log","0.2max","1max"), "facecolor":0.0, "alpha":1.0}
        },
        "kn_grb_skymap": {
            # "hist_nx": 125, "hist_ny": 75,
            # "cm": {"color": 'lime', "marker": "d"},
            # "ysize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # "pcolormesh":{"cmap":'viridis', "norm":("log","0.2max","1max"), "facecolor":0.0, "alpha":1.0}
        },
        "kn_w_skymap": {
            # "hist_nx": 125, "hist_ny": 75,
            # "cm": {"color": 'red', "marker": "x"},
            # "ysize": {"capsize": 2, "color": "orange", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "orange", "lw": 0.5},
            # "pcolormesh":{"cmap":'viridis', "norm":("log","0.2max","1max"), "facecolor":0.0, "alpha":1.0}
        },
        "kn_skymap_ratio": {
            # "pcolormesh": {"cmap": 'viridis', "norm": ("twoslope", 0, 1, 4), "facecolor": 0.0, "alpha": 1.0}  # {}
        },
        "xlim": (-1., 1.0), "ylim": (-1.0, 1.0),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "black",
        "histy_backgound_color": "black",
        "plot_grids":True,
        "histx_lim":(1e-4,1e-2),
        "histy_lim":(1e-4,1e-2)
    }
    plot_one_skymap_with_dists(task_to_plot=task_to_plot, settings=settings)

if __name__ == '__main__':
    main()