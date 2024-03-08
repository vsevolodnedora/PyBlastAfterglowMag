# import PyBlastAfterglowMag
import copy
import shutil

import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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
afterglowpy = True
try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

def run(working_dir:str, struct:dict, P:dict, type:str="a") -> PBA.PyBlastAfterglow:
    # clean he temporary direcotry
    if os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)

    # generate initial data for blast waves
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80,
                                       n_layers_a=1 if struct["struct"]=="tophat" else 20)

    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir+"id_pw.h5")

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir+"id_a.h5")

    # create new parfile
    P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
    PBA.parfile_tools.create_parfile(working_dir=working_dir, P=P)

    # instantiate PyBlastAfterglow
    pba = PBA.interface.PyBlastAfterglow(workingdir=working_dir)

    # run the code with given parfile
    pba.run(
        path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
        loglevel="info"
    )

    # process skymap
    if (pba.GRB.opts["do_skymap"]=="yes"):
        conf = {"nx":128, "ny":64, "extend_grid":1, "fwhm_fac":0.5, "lat_dist_method":"integ",
                "intp_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }, # "gaussian"
                "hist_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }}
        prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
        prep.process_singles(infpaths=working_dir+"raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=True)

    return pba

def mrg(dict2:dict, dict1:dict):
    return {k: v | dict2[k] for k, v in dict1.items() if k in dict2.keys()}

def run_afgpy(freq:float, struct:dict, pba:PBA.PyBlastAfterglow):
    # afterglowpy
    Z = {'jetType':     grb.jet.TopHat if struct["struct"] == "tophat" else grb.jet.Gaussian,     # Top-Hat jet
         'specType':    0,                  # Basic Synchrotron Spectrum
         'counterjet':  1,
         'spread':      7 if pba.GRB.opts["method_spread"] != "None" else -1,
         'thetaObs':    pba.main_pars["theta_obs"],   # Viewing angle in radians
         'E0':          struct["Eiso_c"], # Isotropic-equivalent energy in erg
         'g0':          struct["Gamma0c"],
         'thetaCore':   struct["theta_c"],    # Half-opening angle in radians
         'thetaWing':   struct["theta_w"],
         'n0':          pba.main_pars["n_ism"],    # circumburst density in cm^{-3}
         'p':           pba.GRB.pars["p_fs"],    # electron energy distribution index
         'epsilon_e':   pba.GRB.pars["eps_e_fs"],    # epsilon_e
         'epsilon_B':   pba.GRB.pars["eps_b_fs"],   # epsilon_B
         'xi_N':        1.0,    # Fraction of electrons accelerated
         'd_L':         pba.main_pars["d_l"], # Luminosity distance in cm
         'z':           pba.main_pars["z"]}   # redshift
    tday = np.logspace(-2, 3, 100)
    tsecond = tday * 3600 * 24
    nu = np.empty(tsecond.shape)
    nu[:] = freq
    Fnu = grb.fluxDensity(tsecond, nu, **Z)
    # _ll, = ax.plot(t / PBA.utils.cgs.day, Fnu, color=i_color, ls='-')
    # lls.append(_ll)
    # lbls.append(r"$\nu=$" + r"${}$ Hz ".format(PBA.utils.latex_float(i_freq))
    #             + r"$\theta_{\rm obs}=$" + r"{:.1f} deg".format(i_thetaobs * 180 / np.pi))
    return (tsecond, Fnu)

def run_jetsim(freq:float, struct:dict, pba:PBA.PyBlastAfterglow):
    # pba_ = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")

    import jetsimpy
    P = dict(
        Eiso = struct["Eiso_c"],        # (Jet) Isotropic equivalent energy
        lf = struct["Gamma0c"],           # (Jet) Lorentz factor
        theta_c = struct["theta_c"],      # (Jet) half opening angle
        n0 = pba.main_pars["n_ism"],             # (ISM) constant number density
        k = 0.0,            # (ISM) wind power index
        A = 0,              # (ISM) wind amplitude
        eps_e = pba.GRB.pars["eps_e_fs"],        # (Radiation) epsilon_e
        eps_b = pba.GRB.pars["eps_b_fs"],       # (Radiation) epsilon_b
        p = pba.GRB.pars["p_fs"],           # (Radiation) electron power index
        theta_v = pba.main_pars["theta_obs"],      # (Radiation) viewing angle
        d = pba.main_pars["d_l"]/PBA.utils.cgs.pc/1e6,         # (radiation) distance (Mpc)
        z = pba.main_pars["z"],            # (radiation) redshift
        xi_N = 1.0,       # (radiation customized) total fraction of accelerated electrons
        b = 0,              # (radiation) magnetic field anisotropy
    )
    # tabulated energy/LF structure
    if struct["struct"] == "tophat":
        theta = np.linspace(0, np.pi, 1000)
        Eiso = np.full_like(theta,P["Eiso"]) # P["Eiso"] * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2)
        lf = np.full_like(theta, P["lf"]) #(P["lf"] - 1) * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2) + 1
        Eiso[theta > struct["theta_w"]] = 0.
    elif struct["struct"] == "gaussian":
        theta = np.linspace(0, np.pi, 1000)
        Eiso = P["Eiso"] * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2)
        lf = (P["lf"] - 1) * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2) + 1
        # Eiso[theta > struct["theta_w"]] = 0.
    else:
        raise KeyError("structure is not recognized")

    if pba.GRB.opts["do_lc"] == 'yes':
        jet1 = jetsimpy.Afterglow(
            theta,           # array of theta
            Eiso,            # array of isotropic equivalent energy
            lf,              # array of initial lorentz factor
            P["A"],          # scale of wind density
            P["n0"],         # constant number density
            spread=pba.GRB.opts["method_spread"] != "None",    # (default = True) with/without spreading effect
            coast=True,      # (default = True) with/without coasting. If this is "False", the initial lorentz factor data will be omitted.
        )
        # define the observing time and frequency
        tday = np.logspace(-2, 3, 100)
        tsecond = tday * 3600 * 24
        # nu = 3e9
        # calculate the afterglow flux density (unit: mJy)
        flux1 = jet1.FluxDensity(
            tsecond,           # [second] observing time span
            freq,                # [Hz]     observing frequency
            copy.deepcopy(P),                 # parameter dictionary for radiation
            rtol=1e-2,         # (default=1e-2) integration error tolerance
            model="sync",      # default radiation model
        )
        # _t, _ref_F_afgpy, _ref_F = np.loadtxt(os.getcwd() + '/' + fname, unpack=True)
        # _ll, = ax.plot(tsecond / PBA.utils.cgs.day, flux1, color=i_color, ls='-')
        # lls.append(_ll)
        # lbls.append(r"$\nu=$" + r"${}$ Hz ".format(PBA.utils.latex_float(i_freq))
        #             + r"$\theta_{\rm obs}=$" + r"{:.1f} deg".format(i_thetaobs * 180 / np.pi))
        return (tsecond, flux1)
    if pba.GRB.opts["do_skymap"] == 'yes':
        theta_c_min = struct["theta_c"]
        arcsinhcells = np.linspace(0, np.arcsinh(np.pi / theta_c_min), 200)
        cells = np.sinh(arcsinhcells) * theta_c_min
        cells[-1] = np.pi
        afterglow = jetsimpy.Afterglow(
            theta,
            Eiso,           # (required) energy distribution
            lf,               # (required) LF distribution
            0.0,
            P["n0"],
            spread=True,    # (optional) lateral expansion effect
            tmin=10,
            tmax=3.2e10,         # (optional) PDE maximum time
            grid=cells,
            coast=True
        )
        # solve image
        image = afterglow.solve_image(
            float(pba.main_opts["skymap_times"].split()[-1]),
            float(pba.main_opts["skymap_freqs"].split()[-1]),
            inum=100, pnum=20, para=P
        )
        # setup image scales & criticle points
        width = image.half_width
        centroid = image.offset
        yscale = image.yscale
        xscale = image.xscale
        # intensity and polarization data
        I = np.array(image.intensity_image).T
        # orientation
        orientation = 2
        # intensity scale
        Inorm = I / I.max()
        Inorm = np.rot90(Inorm, orientation) # rotate


        # im = ax.imshow(Inorm, interpolation='gaussian', cmap="inferno", origin='lower', extent=[-width,     width, -width, width])
        return image

def plot_jetsim_skymap(ax, image):
    # setup image scales & criticle points
    width = image.half_width
    centroid = image.offset
    yscale = image.yscale
    xscale = image.xscale

    # orientation
    orientation = 2

    # intensity and polarization data
    I = np.array(image.intensity_image).T

    # intensity scale
    Inorm = I / I.max()
    Inorm = np.rot90(Inorm, orientation) # rotate

    # plot intensity
    im = ax.imshow(Inorm, interpolation='gaussian', cmap="inferno", origin='lower', extent=[-width, width, -width, width])
    im.set_clim(vmin=0.0, vmax=1)

    # offset and scale
    xc = centroid
    yc = 0.0
    sigmax = xscale
    sigmay = yscale
    ax.scatter(xc, yc, marker="+", color="white", label="centroid")
    ax.add_patch(Ellipse((xc, yc), sigmax * 2, sigmay * 2, edgecolor="white", fill=False,   linestyle="--", linewidth=1))


def compare_lcs(struct:dict, pp:dict, plot:dict):

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5.6, 4.2))
    ax = axes

    working_dir = os.getcwd() + '/tmp1/'
    fig_dir = os.getcwd() + '/figs/'
    # prepare initial data (piecewise and adaptive)

    lls, lbls = [], []
    # for (i_thetaobs, i_freq, i_ls) in [
    #     # (thetaObs, freqobs, "blue"),
    #     (0.16, 1e9, "-"),
    #     (0, 1e18, "--"),
    #     (0.16, 1.e18, "-."),
    #     (0, 1e9, ":")
    # ]:

    for iter in plot["iters"]:
        i_freq,i_thetaobs,i_ls = iter["freq"], iter["theta_obs"], iter["ls"]
        # default : Analytic
        pba = run(working_dir=working_dir,struct=struct, type="a", P=mrg(pp,{
            "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
            "grb":{"method_ele_fs":"analytic",
                   "method_synchrotron_fs":"WSPN99"}}))
        ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
                pba.GRB.get_lc(freq=i_freq), color='pink', ls=i_ls, lw=1.,
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        # default : Semi-Analytic
        pba = run(working_dir=working_dir,struct=struct, type="a", P=mrg(pp,{
            "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
            "grb":{"method_ele_fs":"mix"}}))
        ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
                pba.GRB.get_lc(freq=i_freq), color='blue', ls=i_ls, lw=1.,
                label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        # default : Numeric
        # pba = run(working_dir=working_dir,struct=struct, type="a", P=mrg(pp,{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{}}))
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color='green', ls=i_ls, lw=1.,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        ''' -------------------- REFERENCES --------------- '''
        t_, f_ = run_afgpy(i_freq, struct, pba)
        _ll, = ax.plot(t_ / PBA.utils.cgs.day, f_, color='black', ls=i_ls, lw=.8)

        t_, f_ = run_jetsim(i_freq, struct, pba)
        _ll, = ax.plot(t_ / PBA.utils.cgs.day, f_, color='gray', ls=i_ls, lw=.8)

        lls.append(_ll)
        lbls.append(r"$\nu=$" + r"${}$ Hz ".format(PBA.utils.latex_float(i_freq))
                    + r"$\theta_{\rm obs}=$" + r"{:.1f} deg".format(i_thetaobs * 180 / np.pi))

        ''' -------------------- PIECE - WISE ------------------- '''

        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_eats":"piece-wise","method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_eats":"piece-wise","method_ne_fs":"usenprime",
        #            "method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_eats":"piece-wise","method_ne_fs":"useNe",
        #            "method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"comovSpec","method_eats":"piece-wise","method_ne_fs":"usenprime",
        #            "method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # ''' -------------------- ADAPTIVE ------------------- '''
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_ne_fs":"usenprime","method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"observFlux","method_ne_fs":"useNe","method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        #
        # pba = run(working_dir=working_dir,struct=struct, type="pw", P=P|{
        #     "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
        #     "grb":{"method_comp_mode":"comovSpec","method_ne_fs":"usenprime","method_ele_fs":"analytic"}})
        # ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.5,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))


        #
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
        #                        parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
        #                        newopts={"method_synchrotron_fs":"Joh06",
        #                                 "method_comp_mode":"observFlux",
        #                                 "method_eats":"piece-wise",
        #                                 "method_ne_fs":"useNe",
        #                                 "method_ele_fs":"analytic",
        #                                 "fname_ejecta_id":"tophat_grb_id_pw.h5",
        #                                 "fname_light_curve":"tophat_{}_pw.h5".format( str(i_thetaobs).replace(".",""))},
        #                        parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        #
        #
        # pba_pw = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_pw.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_pw.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_pw.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_pw.clear()
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
        #                                          newopts={"method_synchrotron_fs":"Joh06",
        #                                                   "method_comp_mode":"observFlux",
        #                                                   "method_eats":"piece-wise",
        #                                                   "method_ne_fs":"usenprime",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_pw.h5",
        #                                                   "fname_light_curve":"tophat_{}_pw.h5".format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_pw = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # # pba.reload_parfile()
        # pba_pw.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_pw.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_pw.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.0,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_pw.clear()
        # # # -------------------------------------------------------
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
        #                                          newopts={"method_synchrotron_fs":"WSPN99",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"piece-wise",
        #                                                   "method_ne_fs":"useNe",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_pw.h5",
        #                                                   "fname_light_curve":"tophat_{}_pw.h5".format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_pw2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_pw2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_pw2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_pw2.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=1.5)
        # pba_pw2.clear()
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
        #                                          newopts={"method_synchrotron_fs":"Joh06",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"piece-wise",
        #                                                   "method_ne_fs":"usenprime",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_pw.h5",
        #                                                   "fname_light_curve":"tophat_{}_pw.h5".format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_pw2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # # pba.reload_parfile()
        # pba_pw2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_pw2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_pw2.GRB.get_lc(freq=i_freq), color=i_color, ls=':', lw=2.0)
        # pba_pw2.clear()

        # --------------------------------------------------------

        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main", newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb", newpars={},
        #                                          newopts={"method_synchrotron_fs":"Dermer09",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"piece-wise",
        #                                                   "method_ne_fs":"useNe",
        #                                                   "method_ele_fs":"numeric",
        #                                                   "fname_ejecta_id":"tophat_grb_id_pw.h5",
        #                                                   "fname_light_curve":"tophat_{}_pw.h5".format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_pw2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # # pba.reload_parfile()
        # pba_pw2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_pw2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_pw2.GRB.get_lc_totalflux(freq=i_freq), color=i_color, ls='-.', lw=1.0)
        # pba_pw2.clear()

        # ========================================================
        #
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                        parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",newpars={},
        #                        newopts={"method_synchrotron_fs":"Joh06",
        #                                 "method_comp_mode":"observFlux",
        #                                 "method_eats":"adaptive",
        #                                 "method_ne_fs":"useNe",
        #                                 "method_ele_fs":"analytic",
        #                                 "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                 "fname_light_curve":"tophat_{}_a.h5"
        #                        .format( str(i_thetaobs).replace(".",""))},
        #                        parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a.GRB.get_lc(freq=i_freq), color=i_color, ls='--', lw=0.5, # densely dashed
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a.clear()
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",newpars={},
        #                                          newopts={"method_synchrotron_fs":"WSPN99",
        #                                                   "method_comp_mode":"observFlux",
        #                                                   "method_eats":"adaptive",
        #                                                   "method_ne_fs":"usenprime",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                                   "fname_light_curve":"tophat_{}_a.h5"
        #                                          .format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a.GRB.get_lc(freq=i_freq), color=i_color, ls='--', lw=1.0,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a.clear()
        # # # -------------------------------------------------------
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",newpars={},
        #                                          newopts={"method_synchrotron_fs":"WSPN99",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"adaptive",
        #                                                   "method_ne_fs":"useNe",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                                   "fname_light_curve":"tophat_{}_a.h5"
        #                                          .format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a2.GRB.get_lc(freq=i_freq), color=i_color, ls='--', lw=1.5, # densely dashed
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a2.clear()
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",newpars={},
        #                                          newopts={"method_synchrotron_fs":"Joh06",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"adaptive",
        #                                                   "method_ne_fs":"usenprime",
        #                                                   "method_ele_fs":"analytic",
        #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                                   "fname_light_curve":"tophat_{}_a.h5"
        #                                          .format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a2.GRB.get_lc(freq=i_freq), color=i_color, ls='--', lw=2,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a2.clear()


        # -------------------------------------------------------
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
        #                                          newpars={"gam1":1,"gam2":1e8,"ngam":250},
        #                                          newopts={"method_synchrotron_fs":"Dermer09",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"adaptive",
        #                                                   "method_ne_fs":"useNe",
        #                                                   "method_ele_fs":"mix",
        #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                                   "fname_light_curve":"tophat_{}_a.h5"
        #                                          .format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a2.GRB.get_lc(freq=i_freq), color=i_color, ls='-.', lw=1,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a2.clear()
        #
        # # -------------------------------------------------------
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="main",newpars={"theta_obs":i_thetaobs},newopts={},
        #                                          parfile="default_parfile.par", newparfile="parfile.par", keep_old=True)
        # PBA.parfile_tools.modify_parfile_par_opt(workingdir=os.getcwd()+"/", part="grb",
        #                                          newpars={"gam1":1,"gam2":1e8,"ngam":250},
        #                                          newopts={"method_synchrotron_fs":"Dermer09",
        #                                                   "method_comp_mode":"comovSpec",
        #                                                   "method_eats":"adaptive",
        #                                                   "method_ne_fs":"useNe",
        #                                                   "method_ele_fs":"numeric",
        #                                                   "fname_ejecta_id":"tophat_grb_id_a.h5",
        #                                                   "fname_light_curve":"tophat_{}_a.h5"
        #                                          .format( str(i_thetaobs).replace(".",""))},
        #                                          parfile="parfile.par", newparfile="parfile.par", keep_old=False)
        # pba_a2 = PBA.interface.PyBlastAfterglow(workingdir=os.getcwd()+"/", parfile="parfile.par")
        # pba_a2.run(path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",loglevel="info")
        # ax.plot(pba_a2.GRB.get_lc_times() / PBA.utils.cgs.day,
        #         pba_a2.GRB.get_lc(freq=i_freq), color=i_color, ls='-.', lw=2,
        #         label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))
        # pba_a2.clear()

        # break


    l11, = ax.plot([-1., -1.], [-2., -2.], color='gray', ls='-', lw=1.)
    l12, = ax.plot([-1., -1.], [-2., -2.], color='black', ls='-', lw=1.)
    l13, = ax.plot([-1., -1.], [-2., -2.], color='green', ls='-',  lw=1.)
    l14, = ax.plot([-1., -1.], [-2., -2.], color='blue', ls='-',  lw=1.)
    # l13, = ax.plot([-1., -1.], [-2., -2.], color='pink', ls='-',  lw=0.5, label=r"\texttt{afterglowpy}")

    legend1 = plt.legend([l11, l12, l14],
                         # [r"\& J\'{o}hannesson+06", r"\& WSPN+99", r"\texttt{afterglowpy}"],
                         [r"\texttt{jetsim}", r"\texttt{afterglopy}", r"\texttt{PyBlastAfterglow}*"],
                         loc="center", bbox_to_anchor=(0.78, 0.56), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)
    legend2 = plt.legend(lls, lbls,
                         loc="center", bbox_to_anchor=(0.4, 0.16), fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")

    # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   # labelsize=plotdic["fontsize"],
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [mJy]")
    ax.set_xlim(1e-1, 1e3)
    ax.set_ylim(1e-9, 1e2)
    # ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(PAPERPATH + save)
    if not plot["figname"] is None: plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)

    plt.show()


def compare_skymaps(struct:dict, pp:dict, plot:dict):
    working_dir = os.getcwd() + '/tmp1/'
    pba = run(working_dir=working_dir,struct=struct, type="a", P=mrg(pp,{
        "main":{}, "grb":{}}))
    skymap = pba.GRB.get_skymap(
        time=float(pp["main"]["skymap_times"].split()[-1]),
        freq=float(pp["main"]["skymap_freqs"].split()[-1])
    )
    tmp={
        "type":"hist",
        "cm": {"color": 'cyan', "marker": "+"},
        "ysize": {"capsize": 2, "color": "cyan", "lw": 0.5},
        "xsize": {"capsize": 2, "color": "cyan", "lw": 0.5},
        "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                       "norm": ("linear", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": 0}
    }
    fig,axes = plt.subplots(ncols=2,nrows=1,figsize=(7,4), sharex='all',sharey='all')
                            # gridspec_kw={'width_ratios': [1, 1],
                            #              'height_ratios': [1, ]})

    # --- reference ---
    image = run_jetsim(freq=float(pp["main"]["skymap_freqs"].split()[-1]), struct=struct, pba=pba)


    # ax_main = axes[1,0]
    # ax_histx = axes[0,0]
    # ax_histy = axes[1,1]
    im = PBA.skymap_plotting_tools.plot_skymap_with_hists(
        skymap=skymap, tmp=tmp, ax_main=axes[0], ax_histx=None, ax_histy=None
    )
    axes[0].scatter(0, 0, marker="*", color="white", label="origin")
    axes[0].set_aspect('equal')

    plot_jetsim_skymap(axes[1], image)
    axes[1].scatter(0, 0, marker="*", color="white", label="origin")
    axes[1].set_aspect('equal')
    for k in range(len(axes)):
        axes[k].set_xlabel(r"$X$ [mas]", fontsize=12)
        if k == 0: axes[k].set_ylabel(r"$Z$ [mas]", fontsize=12)
        axes[k].tick_params(axis='both', which='both', labelleft=True,
                            labelright=False, tick1On=True, tick2On=True,
                            labelsize=12,
                            direction='in',
                            bottom=True, top=True, left=True, right=True)
        axes[k].minorticks_on()
        axes[k].axhline(y=0, linestyle=':', linewidth=0.4)
        axes[k].axvline(x=0, linestyle=':', linewidth=0.4)
        axes[k].set_facecolor('black')
        # axes[k].set_xlim(-1,3)
        # axes[k].set_ylim(-1,3)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # compare_jetsim_afgpy(
    #     struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
    #     pp = dict(main=dict(n_ism = 1e-2),
    #               grb=dict()),
    #     plot = dict(iters=[
    #         dict(theta_obs=0.16,freq=1.e9,ls='-'),
    #         dict(theta_obs=0.0,freq=1.e18,ls='--'),
    #         dict(theta_obs=0.16,freq=1.e18,ls='-.'),
    #         dict(theta_obs=0.0,freq=1.e9,ls=':')
    #     ], figname="fig1")
    # )
    # compare_jetsim_afgpy(
    #     struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
    #     pp = dict(main=dict(n_ism = 0.00031,d_l = 1.27e+26, z = 0.0099),
    #               grb=dict(eps_e_fs = 0.0708, eps_b_fs = 0.0052, p_fs = 2.16)),
    #     plot = dict(iters=[
    #         dict(theta_obs=0.3752,freq=1.e9,ls='-'),
    #         dict(theta_obs=0.0,freq=1.e18,ls='--'),
    #         dict(theta_obs=0.3752,freq=1.e18,ls='-.'),
    #         dict(theta_obs=0.0,freq=1.e9,ls=':')
    #         ], figname="fig2")
    # )
    compare_skymaps(
        struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
        pp = dict(main=dict(n_ism = 0.00031,d_l = 1.27e+26, z = 0.0099, theta_obs = 0.3752,
                            skymap_freqs=f"array 1.e9", skymap_times=f"array {230.*86400.}"),
                  grb=dict(eps_e_fs = 0.0708, eps_b_fs = 0.0052, p_fs = 2.16,
                           do_skymap="yes",do_lc="no",
                           method_ele_fs="analytic",  method_synchrotron_fs="WSPN99")),
        plot = dict()
    )


