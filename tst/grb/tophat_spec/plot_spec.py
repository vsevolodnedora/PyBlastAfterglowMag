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
#     from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
#     from PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")

# try:
#     import package.src.PyBlastAfterglowMag as PBA
# except ImportError:
#     try:
#         import PyBlastAfterglowMag as PBA
#     except:
#         raise ImportError("Cannot import PyBlastAfterglowMag")
import package.src.PyBlastAfterglowMag as PBA
curdir = os.getcwd()+'/'


def plot_spec():

    # prepare initial data (piecewise and adaptive)
    struct = {"struct":"tophat", "Eiso_c":1.e52, "Gamma0c": 150., "M0c": -1.,"theta_c": 0.1, "theta_w": 0.1}
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80, n_layers_a=1)
    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=curdir+"tophat_grb_id_pw.h5")
    
    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=curdir+"tophat_grb_id_a.h5")

    i_thetaobs = 0.

    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
                                             parfile="parfile_def.par", newparfile="parfile.par", keep_old=True)
    PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
                                             newpars={"gam1":1,"gam2":1e8,"ngam":250},
                                             newopts={"method_synchrotron_fs":"Dermer09",
                                                      "method_comp_mode":"comovSpec",
                                                      "method_eats":"adaptive",
                                                      "method_ne_fs":"useNe",
                                                      "method_ele_fs":"numeric",
                                                      "fname_ejecta_id":"tophat_grb_id_a.h5",
                                                      "fname_spec":"tophat_spec_{}_a.h5",
                                                      "fname_light_curve":"tophat_{}_a.h5"
                                             .format( str(i_thetaobs).replace(".",""))},
                                             parfile="parfile.par", newparfile="parfile.par", keep_old=False)
    pba_a2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
    # pba_a2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")

    spec = pba_a2.GRB.get_lc(key="n_ele_fs",xkey="gams",ykey="times_gams",freq=None,time=1e10,ishell=0,ilayer=0,spec=True)

    gams = pba_a2.GRB.get_gams(unique=True)
    ts = pba_a2.GRB.get_grid(key="times_gams",spec=True)

    plt.loglog(gams, spec)
    plt.show()

    spec[~np.isfinite(spec)] = 1e-100
    norm = LogNorm(vmin=spec.max() * 1e-5, vmax=spec.max())

    fig, ax = plt.subplots(ncols=1,nrows=1)
    _c = ax.pcolormesh(gams, ts, spec, cmap='viridis', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylabel(r'comoving time [s]', fontsize=12)
    ax.set_xlabel(r'$\gamma$', fontsize=12)
    ax.set_xlim(.5, gams[-1])
    cbar = fig.colorbar(_c, ax=ax, shrink=0.9, pad=.01, label=r'$N_{\rm tot}^{-1} [ \gamma_e^2 dn/d\gamma$ ]')
    cbar.ax.tick_params(labelsize=12)
    ax.tick_params(direction="in", which="both")
    plt.tight_layout()
    plt.show()

    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5.6, 4.2))
    # ax = axes
if __name__ == '__main__':
    plot_spec()