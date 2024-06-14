# import PyBlastAfterglowMag
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os

import package.src.PyBlastAfterglowMag as PBA

try:
    import afterglowpy as grb
except:
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

try:
    import jetsimpy
except:
    print("Error! could not import jetsim")

# temprary dir to save PyBlastAfterglow output into
working_dir = os.getcwd() + '/working_dirs/'
fig_dir = os.getcwd() + '/figs/'

# def run(working_dir:str, struct:dict, P:dict, type:str="a") -> PBA.PyBlastAfterglow:
#     # clean he temporary direcotry
#     if os.path.isdir(working_dir):
#         shutil.rmtree(working_dir)
#     os.mkdir(working_dir)
#
#     # generate initial data for blast waves
#     pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80,
#                                        n_layers_a=1 if struct["struct"]=="tophat" else 20)
#
#     # save piece-wise EATS ID
#     id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
#     pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir+"id_pw.h5")
#
#     # save adaptive EATS ID
#     id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
#     pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir+"id_a.h5")
#
#     # create new parfile
#     P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
#     PBA.parfile_tools.create_parfile(working_dir=working_dir, P=P)
#
#     # instantiate PyBlastAfterglow
#     pba = PBA.interface.PyBlastAfterglow(workingdir=working_dir)
#
#     # run the code with given parfile
#     pba.run(
#         path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
#         loglevel="err"
#     )
#
#     # process skymap
#     if (pba.GRB.opts["do_skymap"]=="yes"):
#         conf = {"nx":128, "ny":64, "extend_grid":1.1, "fwhm_fac":0.5, "lat_dist_method":"integ",
#                 "intp_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }, # "gaussian"
#                 "hist_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }}
#         prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
#         prep.process_singles(infpaths=working_dir+"raw_skymap_*.h5",
#                              outfpath=pba.GRB.fpath_sky_map,
#                              remove_input=True)
#
#     return pba

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

    P = dict(
        Eiso = struct["Eiso_c"],        # (Jet) Isotropic equivalent energy
        lf = struct["Gamma0c"],           # (Jet) Lorentz factor
        # theta_c = struct["theta_c"],      # (Jet) half opening angle
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
        theta = np.linspace(0, struct["theta_w"]*1.01, 500)
        Eiso = np.full_like(theta,P["Eiso"]) # P["Eiso"] * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2)
        lf = np.full_like(theta, P["lf"]) #(P["lf"] - 1) * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2) + 1
        lf[theta > struct["theta_c"]] = 0.
        Eiso[theta > struct["theta_c"]] = 0.
    elif struct["struct"] == "gaussian":
        theta = np.linspace(0, np.pi, 1000)
        Eiso = P["Eiso"] * np.exp(- 0.5 * (theta / struct["theta_c"]) ** 2)
        lf = (P["lf"] - 1) * np.exp(- 0.5 * (theta / struct["theta_c"]) ** 2) + 1
        # Eiso[theta > struct["theta_w"]] = 0.
    else:
        raise KeyError("structure is not recognized")

    if pba.GRB.opts["save_dynamics"] == "yes":
        jet1 = jetsimpy.Afterglow(
            theta,           # array of theta
            Eiso,            # array of isotropic equivalent energy
            lf,              # array of initial lorentz factor
            # 0.,
            P["A"],          # scale of wind density
            P["n0"],         # constant number density
            tmin=pba.main_pars["tb0"],
            tmax=pba.main_pars["tb1"],
            spread=pba.GRB.opts["method_spread"] != "None",    # (default = True) with/without spreading effect
            coast=True,      # (default = True) with/without coasting. If this is "False", the initial lorentz factor data will be omitted.
        )
        # t_arr = np.logspace(np.log10(pba.main_pars["tb0"]),
        #                     np.log10(pba.main_pars["tb1"]),
        #                     pba.main_pars["ntb"])
        # mom = jet1.beta_gamma(t_arr,)
        return jet1

    if pba.GRB.opts["do_lc"] == 'yes':
        jet1 = jetsimpy.Afterglow(
            theta,           # array of theta
            Eiso,            # array of isotropic equivalent energy
            lf,              # array of initial lorentz factor
            # 0.,
            P["A"],          # scale of wind density
            P["n0"],         # constant number density
            tmin=pba.main_pars["tb0"],
            tmax=pba.main_pars["tb1"],
            spread=pba.GRB.opts["method_spread"] != "None",    # (default = True) with/without spreading effect
            coast=False      # (default = True) with/without coasting. If this is "False", the initial lorentz factor data will be omitted.
            # tail=False
        )
        # define the observing time and frequency
        tday = np.logspace(-1, 3, 100)
        tsecond = tday * 3600 * 24
        # nu = 3e9
        # calculate the afterglow flux density (unit: mJy)
        flux1 = jet1.FluxDensity(
            tsecond,           # [second] observing time span
            freq,                # [Hz]     observing frequency
            copy.deepcopy(P),                 # parameter dictionary for radiation
            rtol=1e-4,         # (default=1e-2) integration error tolerance
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
    im.set_clim(vmin=0.01, vmax=1)

    # offset and scale
    xc = centroid
    yc = 0.0
    sigmax = xscale
    sigmay = yscale
    ax.scatter(xc, yc, marker="+", color="white", label="centroid")
    ax.add_patch(Ellipse((xc, yc), sigmax * 2, sigmay * 2, edgecolor="white", fill=False,   linestyle="--", linewidth=1))



def compare_lcs(pp:dict, plot:dict,working_dir:str, run:bool=True):

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5.2, 3.4))
    ax = axes

    # prepare initial data (piecewise and adaptive)

    lls, lbls = [], []
    # for (i_thetaobs, i_freq, i_ls) in [
    #     # (thetaObs, freqobs, "blue"),
    #     (0.16, 1e9, "-"),
    #     (0, 1e18, "--"),
    #     (0.16, 1.e18, "-."),
    #     (0, 1e9, ":")
    # ]:

    for i, iter in enumerate(plot["iters"]):
        i_freq,i_thetaobs,i_ls = iter["freq"], iter["theta_obs"], iter["ls"]

        # default : Analytic
        if plot["plot_analytic"]:
            pba = PBA.wrappers.run_grb(working_dir=working_dir+f'analytic_{i}/',run=run, P=mrg(pp, {
                "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
                "grb":{"method_ele_fs":"analytic",
                       "method_synchrotron_fs":"WSPN99"}}))
            ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
                    pba.GRB.get_lc(freq=i_freq), color='red', ls=i_ls, lw=1.,
                    label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        # default : Semi-Analytic
        if plot["plot_semi_analytic"]:
            pba = PBA.wrappers.run_grb(working_dir=working_dir+f'mix_{i}/',run=run, P=mrg(pp, {
                "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
                "grb":{"method_ele_fs":"mix"}})) # "method_synchrotron_fs":"GSL"
            ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
                    pba.GRB.get_lc(freq=i_freq), color='green', ls=i_ls, lw=1.,
                    label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        # default : Numeric
        if plot["plot_numeric"]:
            pba = PBA.wrappers.run_grb(working_dir=working_dir+f'num_{i}/', run=run, P=mrg(pp, {
                "main":{"theta_obs":i_thetaobs,"lc_freqs":f"array {i_freq}"},
                "grb":{}}))
            ax.plot(pba.GRB.get_lc_times() / PBA.utils.cgs.day,
                    pba.GRB.get_lc(freq=i_freq), color='blue', ls=i_ls, lw=1.,
                    label=r"$\theta_{obs}=$" + "{:.2f}".format(i_thetaobs) + r" $\nu$={:.1e}".format(i_freq))

        ''' -------------------- REFERENCES --------------- '''

        t_, f_ = run_afgpy(i_freq, pp["grb"]["structure"], pba)
        _ll, = ax.plot(t_ / PBA.utils.cgs.day, f_, color='black', ls=i_ls, lw=.8)

        t_, f_ = run_jetsim(i_freq,  pp["grb"]["structure"], pba)
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
    l13, = ax.plot([-1., -1.], [-2., -2.], color='red', ls='-',  lw=1.)
    l14, = ax.plot([-1., -1.], [-2., -2.], color='green', ls='-',  lw=1.)
    l15, = ax.plot([-1., -1.], [-2., -2.], color='blue', ls='-',  lw=1.)
    lines, labels = [l11, l12], [r"\texttt{jetsimpy}", r"\texttt{afterglopy}"]
    if plot["plot_analytic"]:
        lines.append(l13)
        labels.append(r"\texttt{PyBlastAfterglow}*")
    if plot["plot_semi_analytic"]:
        lines.append(l14)
        labels.append(r"\texttt{PyBlastAfterglow}")
    if plot["plot_numeric"]:
        lines.append(l15)
        labels.append(r"\texttt{PyBlastAfterglow}")
    legend1 = plt.legend(lines, labels,
                         loc="center", bbox_to_anchor=plot["bbox_to_anchor_1"], fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)
    legend2 = plt.legend(lls, lbls,
                         loc="center", bbox_to_anchor=plot["bbox_to_anchor_2"], fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    ax.add_artist(legend2)
    # ax.plot([-1., -1.], [-2., -2.], color='gray', ls='--', label=r"afterglowpy")

    # ax.set_title(r"Comparison with afterglowpy (e.g. Fig.2)")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=13,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]",fontsize=13)
    ax.set_ylabel(r"$F_{\nu}$ [mJy]",fontsize=13)
    ax.set_xlim(plot["xlim"])
    ax.set_ylim(plot["ylim"])
    # ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(PAPERPATH + save)
    if "figname" in plot.keys():
        plt.savefig(fig_dir+plot["figname"]+'.png',dpi=256)
        plt.savefig(fig_dir+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()

def compare_dyn(pp:dict, plot:dict,working_dir:str):
    pba = PBA.wrappers.run_grb(working_dir=working_dir, P=copy.deepcopy(pp), loglevel="err")
    jetsimpy = run_jetsim(freq=-1.,struct= pp["grb"]["structure"],pba=pba)
    fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(5.,3.))

    def plot_momentum(ax, color, ilayer):
        ax.plot(
            pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=ilayer),
            pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=ilayer),
            color=color,ls='-',lw=1.
        )
        times = np.logspace(np.log10(jetsimpy.tmin),np.log10(jetsimpy.tmax),200)
        thetas = np.array(
            [pba.GRB.get_dyn_arr(v_n="theta",ishell=0,ilayer=ilayer)[
                 PBA.utils.find_nearest_index(
                     pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=ilayer), t_
                 )
             ]
             for t_ in times]
        )
        # mom = jetsimpy.beta_gamma(
        #     t=times,
        #     theta=[
        #         # thetas[0]
        #         # pba.GRB.get_dyn_obj()[f"shell=0 layer={ilayer}"].attrs["theta_c_l"]
        #         pba.GRB.get_dyn_obj()[f"shell=0 layer={ilayer}"].attrs["theta_c"]
        #         +
        #     ]
        # )
        theta_c = pba.GRB.get_dyn_obj()[f"shell=0 layer={ilayer}"].attrs["theta_c"]
        mom = np.array([jetsimpy.beta_gamma(
            t_,
            theta_c + .5 * theta_ - .5 * thetas[0]).flatten() for (t_,theta_) in zip(times,thetas)])

        ax.plot(
            times,
            mom,
            color=color,ls='--',lw=0.7
        )

    def plot_(ax):
        if  pp["grb"]["structure"]["struct"] == 'tophat':
            plot_momentum(ax,color=plot["colors"],ilayer=plot["layers"])
        elif  pp["grb"]["structure"]["struct"] == 'gaussian':
            pba = PBA.wrappers.run_grb(working_dir=working_dir, P=copy.deepcopy(pp))
            jetsimpy = run_jetsim(freq=-1.,struct=pp["grb"]["structure"],pba=pba)
            for ilayer, color in zip(plot["layers"],plot["colors"]):
                plot_momentum(ax,color=color,ilayer=ilayer)

    plot_(ax)

    # inset axes....
    if plot["include_zoom_in"]:
        x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
        axins = ax.inset_axes(
            [0.05, 0.05, 0.45, 0.45],
            xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
        plot_(axins)
        # sub region of the original image
        x1, x2, y1, y2 = 5e7, 2e9, 1e-1, 5.
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xscale('log')
        axins.set_yscale('log')
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        axins.minorticks_on()
        axins.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        ax.indicate_inset_zoom(axins, edgecolor="black")

    ax.set_xlim(1e3, 1e10)
    ax.set_ylim(1e-2, 5e2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r"$\Gamma\beta$", fontsize=12)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax.minorticks_on()

    l11, = ax.plot([-1., -1.], [-2., -2.], color='black', ls='--', lw=1.)
    l12, = ax.plot([-1., -1.], [-2., -2.], color='black', ls='-', lw=1.)

    legend1 = plt.legend([l11, l12],
                         # [r"\& J\'{o}hannesson+06", r"\& WSPN+99", r"\texttt{afterglowpy}"],
                         [r"\texttt{jetsimpy}", r"\texttt{PyBlastAfterglow}"],
                         loc="upper right", # bbox_to_anchor=(0.78, 0.56),
                         fancybox=False, shadow=False, ncol=1,
                         fontsize=10, framealpha=0, borderaxespad=0., frameon=False)
    # legend2 = plt.legend(lls, lbls,
    #                      loc="center", bbox_to_anchor=(0.4, 0.16), fancybox=False, shadow=False, ncol=1,
    #                      fontsize=10, framealpha=0, borderaxespad=0., frameon=False)

    ax.add_artist(legend1)
    # ax.add_artist(legend2)
    plt.tight_layout()
    if "figname" in plot.keys():
        plt.savefig(fig_dir+plot["figname"]+'.png',dpi=256)
        plt.savefig(fig_dir+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()

def compare_skymaps(pp:dict, plot:dict):
    pba = PBA.wrappers.run_grb(working_dir=working_dir, P=mrg(pp, {
        "main":{}, "grb":{}}))
    skymap = pba.GRB.get_skymap(
        time=float(pp["main"]["skymap_times"].split()[-1]),
        freq=float(pp["main"]["skymap_freqs"].split()[-1])
    )
    tmp={
        "type":"hist",
        "cm": {"color": 'cyan', "marker": "+"},
        "ysize": {"capsize": 1, "color": "white", "lw": 0.3},
        "xsize": {"capsize": 1, "color": "white", "lw": 0.3},
        "pcolormesh": {"cmap": 'inferno', "set_under": None, "set_over": None, "set_rasterized": True,
                       "norm": ("linear", "0.01max", "1max"), "facecolor": 0, "alpha": 1.0, "isnan": 0}
    }
    fig,axes = plt.subplots(ncols=2,nrows=1,figsize=(7,4), sharex='all',sharey='all')
                            # gridspec_kw={'width_ratios': [1, 1],
                            #              'height_ratios': [1, ]})

    # --- reference ---
    image = run_jetsim(freq=float(pp["main"]["skymap_freqs"].split()[-1]), struct=pp["grb"]["structure"], pba=pba)


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
        axes[k].tick_params(axis='both', which='both', labelbottom=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12, color="white",
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
        axes[k].minorticks_on()
        # axes[k].axhline(y=0, linestyle=':', linewidth=0.4)
        # axes[k].axvline(x=0, linestyle=':', linewidth=0.4)
        axes[k].set_facecolor('black')
        # axes[k].set_xlim(-1,3)
        # axes[k].set_ylim(-1,3)
    axes[0].text(.1, .1, r"\texttt{PyBlastAfterglow}*",
                 color='white',
                 transform=axes[0].transAxes, fontsize=12)
    axes[1].text(.1, .1, r"\texttt{jetsimpy}",
                 color='white',
                 transform=axes[1].transAxes, fontsize=12)
    plt.tight_layout()
    if "figname" in plot.keys():
        plt.savefig(fig_dir+plot["figname"]+'.png',dpi=256)
        plt.savefig(fig_dir+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()


if __name__ == '__main__':
    show = True
    ''' tophat jet '''
    # struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    # compare_dyn(
    #     pp = dict(main=dict(n_ism = 1e-2),
    #               grb=dict(structure=struct,eats_type='a',save_dynamics='yes',do_mphys_in_situ="no",do_lc = "no")),
    #     plot=dict(
    #         layers=0,
    #         colors='blue',
    #         include_zoom_in=False,
    #         show=show,
    #         figname="dyn_tophat"
    #     ),
    #     working_dir=working_dir+"tmp_tophat_dyn_",
    # )
    struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    compare_lcs(
        pp = dict(main=dict(n_ism = 1e-2,ntb=3000),
                  grb=dict(structure=struct,eats_type='a',ebl_tbl_fpath='none')),
        plot = dict(plot_analytic=True,plot_semi_analytic=True,plot_numeric=False,iters=[
            dict(theta_obs=0.16,freq=1.e9,ls='-'),
            dict(theta_obs=0.0,freq=1.e18,ls='--'),
            dict(theta_obs=0.16,freq=1.e18,ls='-.'),
            dict(theta_obs=0.0,freq=1.e9,ls=':')
        ], xlim=(1e-1, 1e3), ylim=(1e-9, 5e2),
                    bbox_to_anchor_2=(0.35, 0.16), bbox_to_anchor_1=(0.78, 0.52), show=show, figname="lcs_tophat"),
        working_dir=working_dir+"tmp_tophat_",
        run=False
    )

    ''' gaussian jet'''
    struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618)
    compare_dyn(
        pp = dict(main=dict(n_ism = 0.00031,d_l = 1.27e+26, z = 0.0099),
                            grb=dict(structure=struct,eats_type='a',eps_e_fs = 0.0708, eps_b_fs = 0.0052, p_fs = 2.16,
                                     save_dynamics='yes',do_mphys_in_situ='no',do_lc = "no",ebl_tbl_fpath='none')),
        plot=dict(
            layers=[0,8,16,19],
            colors=['blue','green','orange','red'],
            include_zoom_in=True,
            show=show,
            figname="dyn_gauss"
        ),
        working_dir=working_dir+"tmp_gauss_dyn_",
    )
    compare_lcs(
        # struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
        pp = dict(main=dict(n_ism = 0.00031,d_l = 1.27e+26, z = 0.0099),
                  grb=dict(structure=struct,eats_type='a',eps_e_fs = 0.0708, eps_b_fs = 0.0052, p_fs = 2.16,ebl_tbl_fpath='none')),
        plot = dict(plot_analytic=True,plot_semi_analytic=True,plot_numeric=False,iters=[
            dict(theta_obs=0.3752,freq=1.e9,ls='-'),
            dict(theta_obs=0.0,freq=1.e18,ls='--'),
            dict(theta_obs=0.3752,freq=1.e18,ls='-.'),
            dict(theta_obs=0.0,freq=1.e9,ls=':')
            ], xlim=(1e-1, 1e3), ylim=(1e-11, 1e2),
                    bbox_to_anchor_2=(0.65, 0.16), bbox_to_anchor_1=(0.78, 0.65), show=show, figname="lcs_gauss"),
        working_dir=working_dir+"tmp_gauss_",
        run=False
    )
    struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618)
    compare_skymaps(
        pp = dict(main=dict(n_ism = 0.00031,d_l = 1.27e+26, z = 0.0099, theta_obs = 0.3752,
                            skymap_freqs=f"array 1.e9", skymap_times=f"array {75.*86400.}"),
                  grb=dict(structure=struct,eats_type='a',eps_e_fs = 0.0708, eps_b_fs = 0.0052, p_fs = 2.16,
                           do_skymap="yes",do_lc="no",
                           method_ele_fs="analytic",  method_synchrotron_fs="WSPN99",
                           ebl_tbl_fpath='none',
                           )),
        plot = dict( show=show, figname="skymaps_gauss" )
    )


