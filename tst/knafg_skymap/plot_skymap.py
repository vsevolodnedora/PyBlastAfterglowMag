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

# from PyBlastAfterglowMag import BPA_METHODS as PBA
from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
from package.src.PyBlastAfterglowMag.interface import cgs
from package.src.PyBlastAfterglowMag.utils import latex_float
from package.src.PyBlastAfterglowMag.id_maker_from_thc_ourflow import \
    (prepare_kn_ej_id_1d,prepare_kn_ej_id_2d)
from package.src.PyBlastAfterglowMag.skymap_tools import \
    (plot_skymaps,plot_skymap_properties_evolution,plot_one_skymap_with_dists,precompute_skymaps)

def main():
    workdir = os.getcwd()+'/'
    # Pre-process the intput data
    fpath_2d_table = "corr_vel_inf_theta_BLh_M13641364_M0_LK_SR.h5"
    outfpath_2d_id = "ejecta_2d_id_BLh_M13641364_M0_LK_SR.h5"
    prepare_kn_ej_id_2d(nlayers=100,corr_fpath=workdir+fpath_2d_table,outfpath=workdir+outfpath_2d_id, dist="pw")

    fpath_1d_table = "hist_vel_inf_BLh_M13641364_M0_LK_SR.dat"
    outfpath_1d_id = "ejecta_1d_id_BLh_M13641364_M0_LK_SR.h5"
    prepare_kn_ej_id_1d(nlayers=30, hist_fpath=workdir+fpath_1d_table, outfpath=workdir+outfpath_1d_id)

    pba = PBA(workingdir=os.getcwd()+'/',readparfileforpaths=True,parfile="parfile.par")

    pba.modify_kn_part_parfile(newpars={}, newopts={ "fname_ejecta_id":outfpath_2d_id} )
    pba.reload_parfile()
    pba.run()

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
        "rerun":False, "precompute": False,
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