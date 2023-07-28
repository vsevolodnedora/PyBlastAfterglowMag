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

# try:
#     from PyBlastAfterglowMag.interface import modify_parfile_par_opt
#     from PyBlastAfterglowMag.interface import PyBlastAfterglow
#     from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#     from PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
#                                            BetFromMom, GamFromMom, MomFromGam)
#     from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     from PyBlastAfterglowMag.skymap_tools import \
#         (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
#                                                            BetFromMom, GamFromMom, MomFromGam)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#         from package.src.PyBlastAfterglowMag.skymap_tools import \
#             (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma,
                                                   BetFromMom, GamFromMom, MomFromGam)
from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
from package.src.PyBlastAfterglowMag.skymap_tools import \
    (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)
def task_kn_skymap_with_dist_one_time():

    workdir = os.getcwd()+'/'
    prepare_grb_ej_id_1d({"Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1,
                          "nlayers_pw": 50, "nlayers_a": 1, "struct":"tophat"},type="a",outfpath="tophat_grb_id.h5")
    modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={},newopts={},
                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    pba = PyBlastAfterglow(workingdir=os.getcwd()+"/", readparfileforpaths=True, parfile="parfile.par")

    # prepare_grb_ej_id_1d({"struct":"tophat",
    #                       "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
    #                       "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":350,"nlayers_a": 10}, type="a",
    #                       outfpath=workdir+"gauss_grb_id.h5")

    # pba = PyBlastAfterglow(workingdir=os.getcwd()+'/',readparfileforpaths=True,parfile="parfile.par")

    pba.reload_parfile()

    pba.run()



    times = [1.,10.,40.,100.,200.]
    workdir = os.getcwd() + '/'
    sim = "BLh_M13641364_M0_LK_SR"
    # # res_dir_kn = "/media/vsevolod/Data/postprocessed5/{}/outflow_1/geo//afterglow/tst_spec_maps/".format(sim)
    # res_dir_kn = workdir#root + "BLh_M13641364_M0_LK_SR" + "/outflow_1/" + MASK + "/" + main_dir + res_dir
    # res_dir_grb = workdir#root + "BLh_M13641364_M0_LK_SR" + "/outflow_1/" + MASK + "/" + main_dir + res_dir
    # prefix_pl = "skymap_kn_example"#"ejecta_nlayers100_theta450_nism000031_p205_epse01_epsb001_epst1_dl00413"
    # prefix_th = "skymap_kn_example"#"ejecta_nlayers100_theta450_nism000031_p205_epse01_epsb001_epst1_dl00413"
    # jet_prefix = "skymap_grb_gauss"#"jet_G_nism000031_thetaobs450_epse007079_epsb000525_p216_epst00_nlayers100_dl00413_"
    # # prefix_root = prefix_pl.replace("theta{}_", "")
    task_to_plot = {
        "workingdir" : workdir,
        "parfile1":  "parfile.par",
        "parfile2":  "parfile.par",
        "time": times[3], "freq": 1e9, "title":"BLh $q=1.00$"
    }

    fname = workdir+ "precomputed_grb_skymap"#"int_"+sim + "_with_hists_time{}_".format(time) + prefix_th
    settings = {
        "gridspec": {
            "width_ratios": (4, 2), "height_ratios": (2, 4),
            "left": 0.14, "right": 0.95, "bottom": 0.1, "top": 0.96, "wspace": 0.05, "hspace": 0.05
        },
        # "figname": fname,
        "paperpath": os.getcwd()+'/', "figfpath": fname, "save_pdf": True, "save_figs": True,
        "show_figs": True,
        "grid": False,
        "figsize": (4.8, 4.8),
        # "workingdir" : workdir,
        "rerun":False, "precompute": True,
        # "precompute_fpath": "/media/vsevolod/Data/postprocessed5/tmp/" + sim + "_single_time{}".format(time) + ".h5",
        "precompute_fpath": fname + ".h5",#res_dir_kn+"int_skymaps_kn_" + prefix_root + ".h5",
        # "precompute_fpath": OUTDIR + fname + ".h5",
        "kn_skymap": {
            # "hist_nx": 1024, "hist_ny": 512, "spec": False,
            # "smooth": {},  # {"type": "gaussian", "sigma": 10},
            # # "smooth": {"type": "gaussian", "sigma": 10},#"smooth": {"type": "uniform", "sigma": 5},
            # "cm": {"color": 'yellow', "marker": "o"},
            # "ysize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "yellow", "lw": 0.5},
            # # "pcolormesh": {"cmap": 'viridis', "norm": ("log", "0.2max", "1max"), "facecolor": 0.0, "alpha": 1.0},
            # "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
            #                "norm": ("log", "0.2max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
            # {}
        },
        "grb_skymap": {
            "hist_nx": 128, "hist_ny": 128, "spec": False,
            "smooth": {},  # {"type": "gaussian", "sigma": 10},
            "cm": {"color": 'cyan', "marker": "+"},
            "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
            "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                           "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
        },
        "kn_grb_skymap": {
            # "hist_nx": 1024, "hist_ny": 512, "spec": False,
            # "smooth": {},  # {"type": "gaussian", "sigma": 10},  # "smooth": {"type": "uniform", "sigma": 5},
            # "cm": {"color": 'lime', "marker": "o"},
            # "ysize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # "xsize": {"capsize": 2, "color": "lime", "lw": 0.5},
            # # "pcolormesh": {"cmap": 'viridis', "norm": ("log", "0.2max", "1max"), "facecolor": 0.0, "alpha": 1.0},
            # "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
            #                "norm": ("log", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": np.nan}
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
        "xlim": (-0.5, 5.0), "ylim": (-2.5, 2.5),
        "title": {"title": "time_fluxratio"},  # "time_fluxratio"
        "cbar_title": r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$',
        "xlabel": "x [mas]", "ylabel": "z [mas]",
        "histx_backgound_color": "black",
        "histy_backgound_color": "black",
        "plot_grids": True,
        "histx_lim":(1e-2, 1e1),
        "histy_lim":(1e-2, 1e1)
    }
    plot_one_skymap_with_dists(task_to_plot=task_to_plot, settings=settings)

def main():
    task_kn_skymap_with_dist_one_time()

if __name__ == '__main__':
    main()