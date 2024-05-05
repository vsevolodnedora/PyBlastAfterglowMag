import copy

from more_itertools import numeric_range

try:
    import package.src.PyBlastAfterglowMag as PBA
    from package.src.PyBlastAfterglowMag.utils import cgs
except ImportError:
    import PyBlastAfterglowMag as PBA
    from PyBlastAfterglowMag.utils import cgs


import os, shutil, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, SymLogNorm, TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator, AutoLocator
from matplotlib import ticker
import matplotlib.ticker as plticker
import numpy as np
import time
from scipy import special
from scipy import integrate
from itertools import product
from multiprocessing import Pool
from tqdm import tqdm

working_dir = os.getcwd() + '/tmp1/'
fig_dir = os.getcwd() + '/figs/'
path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out"

def d2d(default: dict, new: dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_

def run(P: dict, run: bool = True, process_skymaps: bool = True, overwrite=False) -> PBA.PyBlastAfterglow:
    """
            conf = {"nx": 64, "ny": 32, "extend_grid": 2, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                "intp_filter": {"type": None, "sigma": 2, "mode": 'reflect'},  # "gaussian"
                "hist_filter": {"type": None, "sigma": 2, "mode": 'reflect'}}
    :param working_dir:
    :param struct:
    :param P:
    :param type:
    :param run:
    :return:
    """
    # clean he temporary direcotry
    working_dir = copy.deepcopy(P["working_dir"])
    del P["working_dir"]
    if overwrite and (run and os.path.isdir(working_dir)):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    # generate initial data for blast waves
    struct = copy.deepcopy(P["struct"])
    del P["struct"]
    pba_id = PBA.id_analytic.JetStruct(
        n_layers_pw=80 if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
        n_layers_a=(1 if struct["struct"] == "tophat" else
                    (20 if not "n_layers_a" in struct.keys() else struct["n_layers_a"])))

    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_pw.h5")

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_a.h5")

    # create new parfile
    P = copy.deepcopy(P)
    P["grb"]["fname_ejecta_id"] = "id_a.h5" if struct["type"] == "a" else "id_pw.h5"
    P["grb"]["method_eats"] = "piece-wise" if struct["type"] == "pw" else "adaptive"

    if (struct["struct"]=="tophat"): P["grb"]["nsublayers"] = 35 # for skymap resolution
    grb_skymap_config = copy.deepcopy(P["grb"]["skymap_conf"])
    del P["grb"]["skymap_conf"]
    PBA.parfile_tools.create_parfile(working_dir=working_dir, P=P)

    # instantiate PyBlastAfterglow
    pba = PBA.interface.PyBlastAfterglow(workingdir=working_dir)

    # check if result exists
    is_exists = True
    if (pba.GRB.opts["do_lc"]=="yes"): is_exists = is_exists & os.path.isfile(pba.GRB.fpath_light_curve)
    if (pba.GRB.opts["do_skymap"]=="yes"): is_exists = is_exists & os.path.isfile(pba.GRB.fpath_sky_map)
    if is_exists and not overwrite:
        print(f"Result already exists: {working_dir}")

    # run the code with given parfile nism10_Eisoc550_Gamma0c6000_thetaw100_pfs22_epsefs01_epsbfs00001_
    if ((is_exists and overwrite and run) or (not is_exists and run)):
        pba.run(
            path_to_cpp_executable=path_to_cpp_executable,
            loglevel="err"
        )

    if (pba.GRB.opts["do_lc"]=="yes" and not os.path.isfile(pba.GRB.fpath_light_curve)):
        raise FileNotFoundError("failed to locate spectrum")
    if (pba.GRB.opts["do_skymap"]=="yes" and not os.path.isfile(pba.GRB.fpath_sky_map)):
        raise FileNotFoundError("failed to locate sky maps")

    # process skymap
    if (process_skymaps and pba.GRB.opts["do_skymap"] == "yes"):
        prep = PBA.skymap_process.ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=False)

    return pba

def plot_spectrum(ej:PBA.Ejecta, task_i:dict):

    figname = ej.workingdir + task_i["figname"] + '.png'
    if (os.path.isfile(figname)):
        return

    # get spectrum
    spec = ej.get_lc(key="fluxdens",
                     xkey="freqs",
                     key_time="times",
                     freq=None,
                     time=None,
                     ishell=None,
                     ilayer=None,
                     spec=False,
                     sum_shells_layers=True)
    xs = ej.get_grid(key="freqs", spec=False)
    times = ej.get_grid(key="times", spec=False)
    # get nu*Fnu
    spec *= xs[:, np.newaxis]


    # plot
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 3), layout='constrained', sharex='all', sharey='all')
    norm = LogNorm(vmin=task_i.get("vmin", spec.max() * 1e-6),
                   vmax=task_i.get("vmax", spec.max() * 2))
    cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
    cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
    cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)
    spec = np.ma.masked_where(spec < norm.vmin, spec)
    spec = np.ma.masked_where(spec > norm.vmax, spec)
    _c = ax.contourf(times, xs, spec, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
    # _c = ax.pcolormesh(times, xs, spec, cmap=cmap, norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"], size=12)
    ax.set_ylabel(task_i["ylabel"], fontsize=12)
    if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
    if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=4, fontsize=12, labelcolor='black',
              framealpha=0.4, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_xlabel(task_i["xlabel"], fontsize=12)
    plt.savefig(figname, dpi=256)
    # plt.savefig(os.getcwd() + '/figs/' + task_i["figname"] + '.pdf')
    if task_i["show"]: plt.show()
    plt.close(fig)

class DispatcherRunPyBlastAfterglow:
    def __init__(self,default_pars:dict, list_pars:list[dict]):
        self.default = default_pars
        self.pars = list_pars
    def __call__(self, idx):
        start_time = time.perf_counter()
        # merge default parameters and the required parameters
        pars = copy.deepcopy( self.pars[idx] )
        P = copy.deepcopy( self.default )
        P["working_dir"] = P["working_dir"] + pars["name"] + '/'
        for par,val in pars.items():
            # check struct
            if par in P["struct"].keys():   P["struct"][par] = val
            if par in P["main"].keys():     P["main"][par] = val
            if par in P["grb"].keys():      P["grb"][par] = val
        # run simulation with given parameters & plot final spectrum
        pba = run(P=P,run=True,process_skymaps=True)
        plot_spectrum(pba.GRB,
                      task_i=dict(ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' F'_{\rm total}$ [cgs]",
                                  xlabel=r"$t_{\rm obs}$ [s]", figname="spectrum",show=False))
        finish_time = time.perf_counter()
        with open(pba.GRB.workingdir+'runtime.txt', 'w') as f:
            f.write("{:.3f}".format(finish_time-start_time))
        print(f"Run #{idx} finished in {finish_time-start_time:.2f} seconds. ")
        pba.clear()

def get_str_val(v_n, val):
    # if v_n == "theta_obs":
    #     val = "{:.0f}".format(val * 180.0 / np.pi) # rad -> deg
    # elif v_n == "d_l":
    #     val = "{:.1f}".format(val / cgs.pc / 1.e9)
    # else:
    #     val = str(val)
    #
    # return val.replace(".", "")
    if ((v_n == "theta_obs") or (v_n == "theta_c") or (v_n == "theta_w")):
        val = "{:.1f}".format(val / np.pi * 180.) # int(val / np.pi * 180.)
    elif ((v_n == "Eiso_c") or ((v_n == "Eiso_c"))):
        val = np.log10(val)
    elif (v_n == "d_l"):
        val = val / 1.e9 / cgs.pc
    else:
        pass
    if len(str(val)) > 7:
        val = "{:.5f}".format(val)
    val = str(val).replace(".", "")
    return val

def grid_run_tophat():
    ncpus = 4
    # struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    P = dict(
        working_dir=os.getcwd() + f"/working_dirs/",
        struct=dict(type="a",struct="tophat", Eiso_c=1.e53, Gamma0c=400., M0c=-1., theta_c=0.1, theta_w=0.1,
                    n_layers_pw=80, n_layers_a=1),
        main=dict(n_ism=1., tb0=1., tb1=1e8, ntb=1000, rtol=5e-7, theta_obs=0, z=0.4245, d_l=2.3e9*cgs.pc,
                  lc_freqs='array logspace 1e9 1e29 96',
                  lc_times='array logspace 10 1e6 128'),
        grb=dict(save_dynamics='yes', save_spec='no', do_lc='yes',
                 # method_nonrel_dist_fs='none',
                 # method_nonrel_dist_rs='none',s
                 eps_e_fs=0.1, eps_b_fs=0.001, p_fs=2.2,
                 gamma_max_fs=4e7, method_gamma_max_fs="useConst",
                 gamma_max_rs=4e7, method_gamma_max_rs="useConst",
                 max_substeps_fs=1000, max_substeps_rs=1000,
                 # ngam_fs=1001,gam1_rs=1,gam2_rs=1e4,ngam_rs=1001,
                 # eps_b_fs = 1e-7,
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 # method_synchrotron_fs="Bessel",
                 # method_synchrotron_rs="Bessel",
                 # method_ele_fs='mix',
                 # method_ele_rs="mix",
                 method_ssc_fs='numeric',
                 method_ssc_rs='numeric',
                 use_ssa_fs='no',
                 use_ssa_rs='no',
                 # num_ele_use_adi_loss_fs='no',
                 # num_ele_use_adi_loss_rs='no',
                 gam1_fs=1., gam2_fs=1e8, ngam_fs=401,
                 gam1_rs=1., gam2_rs=1e8, ngam_rs=401,
                 freq1_fs=1e6, freq2_fs=1e33, nfreq_fs=401,
                 freq1_rs=1e6, freq2_rs=1e33, nfreq_rs=401,
                 # ebl_tbl_fpath="none"
                 skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
                                  intp_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect'),  # "gaussian"
                                  hist_filter=dict(type='gaussian', size=2, sigma=1.5, mode='reflect'))
                 )
    )
    if not os.path.isdir(P["working_dir"]):
        os.mkdir(P["working_dir"])

    # grb190114c
    # Distance is from H.E.S.S collaboration paper
    # n_ism = 1 Eiso = 3e53
    # Angle: expected 13 https://www.aanda.org/articles/aa/pdf/2022/03/aa41788-21.pdf
    pars = {
        "n_ism":    np.array([1.0, 0.5, 0.1, 0.05, 0.01, 0.001]),
        # "theta_obs":np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
        "Eiso_c":   np.array([1e51, 1.e52, 1.e53, 1e54, 1e55]),
        "Gamma0c":  np.array([100., 300., 600., 1000.]),
        # "theta_c":  np.array([5., 10., 15., 20.]) * np.pi / 180.,
        "theta_w":  np.array([5., 10., 15., 20.]) * np.pi / 180.,
        "p_fs":     np.array([2.2, 2.4, 2.6, 2.8]),  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
        "eps_e_fs": np.array([1e-1, 1e-2, 1e-3]),  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        "eps_b_fs": np.array([1e-2, 1e-3, 1e-4, 1e-5]),  # [0.001, 0.005, 0.01, 0.05, 0.1],
    }
    ranges = [pars[key] for key in pars.keys()]
    pars_cobinations = []
    all_combinations = product(*ranges)
    for combination in all_combinations:
        if (("theta_c" in pars) and ("theta_w" in pars)):
            # skip if jet core exceeds jet winds (artifact of permutations)
            if combination[list(pars.keys()).index("theta_c")] > combination[list(pars.keys()).index("theta_w")]:
                continue
        # create a dict with {par:value} for each parameter and value in current permitation
        pars_cobinations.append({par:val for par,val in zip(list(pars.keys()), combination)})
        # create a str containing the par_value for each par and value (used later to label the run)
        pars_cobinations[-1]["name"] = "".join([par.replace("_","")+get_str_val(par,val)+'_'
                                               for par,val in zip(pars.keys(),combination)])
    print(f"Total runs to be performed {len(pars_cobinations)}")

    # initialize the run class
    run_grb = DispatcherRunPyBlastAfterglow(P, pars_cobinations)

    # perform simulations in parallel
    start_time = time.perf_counter()
    with Pool(ncpus) as pool:
        pool.map(run_grb, range(len(pars_cobinations)))
    finish_time = time.perf_counter()

    print(f"Program finished in {finish_time-start_time:.2f} seconds. ")


if __name__ == '__main__':
    grid_run_tophat()