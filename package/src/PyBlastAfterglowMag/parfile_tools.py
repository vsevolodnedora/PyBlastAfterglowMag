import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy
import subprocess
from shutil import copyfile
from multiprocessing import Pool
import itertools

from .utils import cgs, get_beta, find_nearest_index



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
        val = val / 1.e6 / cgs.pc
    else:
        pass
    if len(str(val)) > 7:
        val = "{:.5f}".format(val)
    val = str(val).replace(".", "")
    return val

def _apply_pars(pars, keys : list, vals : list, prefix_key : str, prefix : str):
    _pars = copy.deepcopy(pars)
    for (key, val) in zip(keys, vals):
        _pars[key] = val
    _pars[prefix_key] = prefix
    return _pars

def set_parlists_for_pars(iter_pars_keys:list[str], iter_pars : dict):
    """ Generate a list of dictionaries containing all possible permutations
        of parameters and returns the list of dictionaries

        `iter_pars_keys` : list[str] list of parameter names that should be iterated over
        :iter_pars: dict that must contain all iter_pars_keys with lists as values e.g.,
            iter_pars_keys : {
            par1: [1,2,3,4]
            par2: [4,6,7,8]
        }
        """
    ranges = [iter_pars[key] for key in iter_pars_keys]
    result = []
    # Generate all possible combinations of values for the 5 parameters
    all_combinations = itertools.product(*ranges)
    for combination in all_combinations:
        # create a dict with {par:value} for each parameter and value in current permitation
        result.append({par:val for par,val in zip(iter_pars_keys,combination)})
        # create a str containing the par_value for each par and value (used later to label the simulation)
        result[-1]["name"] = "".join([par.replace("_","")+get_str_val(par,val)+'_'
                                      for par,val in zip(iter_pars_keys,combination)])
    return result

    # print(f"Total combinations {len(all_combinations)}")

def OLD_set_parlists_for_pars(iter_pars_keys, iter_pars : dict, fname : str)->list[dict]:
    pars = {}
    prefix_key = "name"
    # fname = pars["ejecta_prefix"]
    # if "[text]" in fname:
    #     fname = fname.replace("[text]", "{:.0f}".format(pars["text"]))
    # if "[nlayers]" in fname:
    #     fname = fname.replace("[nlayers]", "{:.0f}".format(int(pars["nlayers"])))
    # if "[d_l]" in fname:
    #     fname = fname.replace("[d_l]", get_str_val("d_l", pars["d_l"]))
    # iterate over parameters that needed to be varied
    if len(iter_pars_keys) > 0:
        n_iter_pars = len(iter_pars_keys)
        # if n_iter_pars != 6: raise ValueError("n_iter_pars = {} is not supported".format(n_iter_pars))
        # print("iterating over {} pars \n{}".format(n_iter_pars, iter_pars_keys))

        ikey = 0
        k0s = iter_pars_keys[ikey]
        par_list = []
        for v0s in iter_pars[k0s]:
            _k0s = "[" + k0s + "]"
            if not _k0s in fname: raise KeyError("_k0s={} is missing in fname={}".format(_k0s,fname))
            fname0 = fname.replace(_k0s, get_str_val(k0s,v0s))  # if _k0s in fname else exit(1)
            if len(iter_pars_keys) == 1:
                par_list.append(_apply_pars(pars,
                                            keys=[k0s],
                                            vals=[v0s],
                                            prefix_key=prefix_key, prefix=fname0))
                continue
                # return par_list

            k1s = iter_pars_keys[ikey+1]
            for v1s in iter_pars[k1s]:
                _k1s = "[" + k1s + "]"
                if not _k1s in fname: raise KeyError("k1s = {} is missing from fname = {}".format(_k1s, fname))
                fname1 = fname0.replace(_k1s, get_str_val(k1s,v1s))  # if _k0s in fname else exit(1)
                if len(iter_pars_keys) == 2:
                    par_list.append(_apply_pars(pars,
                                                keys=[k0s, k1s],
                                                vals=[v0s, v1s],
                                                prefix_key=prefix_key, prefix=fname1))
                    continue
                    # return par_list

                k2s = iter_pars_keys[ikey+2]
                for v2s in iter_pars[k2s]:
                    _k2s = "[" + k2s + "]"
                    if not _k2s in fname: raise KeyError("k2s = {} is missing from fname = {}".format(_k1s, fname))
                    fname2 = fname1.replace(_k2s, get_str_val(k2s,v2s))  # if _k0s in fname else exit(1)
                    if len(iter_pars_keys) == 3:
                        par_list.append(_apply_pars(pars,
                                                    keys=[k0s, k1s, k2s],
                                                    vals=[v0s, v1s, v2s],
                                                    prefix_key=prefix_key, prefix=fname2))
                        continue

                    k3s = iter_pars_keys[ikey+3]
                    for v3s in iter_pars[k3s]:
                        _k3s = "[" + k3s + "]"
                        if not _k3s in fname: raise KeyError("k3s = {} is missing from fname = {}".format(_k1s, fname))
                        fname3 = fname2.replace(_k3s, get_str_val(k3s,v3s))  # if _k0s in fname else exit(1)
                        if len(iter_pars_keys) == 4:
                            par_list.append(_apply_pars(pars,
                                                        keys=[k0s, k1s, k2s, k3s],
                                                        vals=[v0s, v1s, v2s, v3s],
                                                        prefix_key=prefix_key, prefix=fname3))
                            continue

                        k4s = iter_pars_keys[ikey+4]
                        for v4s in iter_pars[k4s]:
                            _k4s = "[" + k4s + "]"
                            if not _k4s in fname: raise KeyError( "k4s = {} is missing from fname = {}".format(_k1s, fname))
                            fname4 = fname3.replace(_k4s, get_str_val(k4s,v4s))  # if _k0s in fname else exit(1)
                            if len(iter_pars_keys) == 5:
                                par_list.append(_apply_pars(pars,
                                                            keys=[k0s, k1s, k2s, k3s, k4s],
                                                            vals=[v0s, v1s, v2s, v3s, v4s],
                                                            prefix_key=prefix_key, prefix=fname4))
                                continue

                            k5s = iter_pars_keys[ikey+5]
                            for v5s in iter_pars[k5s]:
                                _k5s = "[" + k5s + "]"
                                if not _k5s in fname: raise KeyError( "k5s = {} is missing from fname = {}".format(_k1s, fname))
                                fname5 = fname4.replace(_k5s, get_str_val(k5s,v5s))  # if _k0s in fname else exit(1)
                                if len(iter_pars_keys) == 6:
                                    par_list.append(_apply_pars(pars,
                                                                keys=[k0s, k1s, k2s, k3s, k4s, k5s],
                                                                vals=[v0s, v1s, v2s, v3s, v4s, v5s],
                                                                prefix_key=prefix_key, prefix=fname5))
                                    continue

                                k6s = iter_pars_keys[ikey + 6]
                                for v6s in iter_pars[k6s]:
                                    _k6s = "[" + k6s + "]"
                                    if not _k5s in fname: raise KeyError(
                                        "k6s = {} is missing from fname = {}".format(_k6s, fname))
                                    fname6 = fname5.replace(_k6s, get_str_val(k6s, v6s))  # if _k0s in fname else exit(1)
                                    if len(iter_pars_keys) == 7:
                                        par_list.append(_apply_pars(pars,
                                                                    keys=[k0s, k1s, k2s, k3s, k4s, k5s, k6s],
                                                                    vals=[v0s, v1s, v2s, v3s, v4s, v5s, v6s],
                                                                    prefix_key=prefix_key, prefix=fname6))
                                        continue

                                    k7s = iter_pars_keys[ikey + 7]
                                    for v7s in iter_pars[k7s]:
                                        _k7s = "[" + k7s + "]"
                                        if not _k7s in fname: raise KeyError(
                                            "k7s = {} is missing from fname = {}".format(_k7s, fname))
                                        fname7 = fname6.replace(_k7s, get_str_val(k7s, v7s))  # if _k0s in fname else exit(1)
                                        if len(iter_pars_keys) == 8:
                                            par_list.append(_apply_pars(pars,
                                                                        keys=[k0s, k1s, k2s, k3s, k4s, k5s, k6s, k7s],
                                                                        vals=[v0s, v1s, v2s, v3s, v4s, v5s, v6s, v7s],
                                                                        prefix_key=prefix_key, prefix=fname7))
                                            continue

                                        k8s = iter_pars_keys[ikey + 8]
                                        for v8s in iter_pars[k8s]:
                                            _k8s = "[" + k8s + "]"
                                            if not _k8s in fname: raise KeyError(
                                                "k8s = {} is missing from fname = {}".format(_k8s, fname))
                                            fname8 = fname7.replace(_k8s, get_str_val(k8s, v8s))  # if _k0s in fname else exit(1)
                                            if len(iter_pars_keys) == 9:
                                                par_list.append(_apply_pars(pars,
                                                                            keys=[k0s, k1s, k2s, k3s, k4s, k5s, k6s, k7s, k8s],
                                                                            vals=[v0s, v1s, v2s, v3s, v4s, v5s, v6s, v7s, v8s],
                                                                            prefix_key=prefix_key, prefix=fname8))
                                                continue

                                            k9s = iter_pars_keys[ikey + 9]
                                            for v9s in iter_pars[k9s]:
                                                _k9s = "[" + k9s + "]"
                                                if not _k9s in fname: raise KeyError(
                                                    "k9s = {} is missing from fname = {}".format(_k9s, fname))
                                                fname9 = fname8.replace(_k9s, get_str_val(k9s, v9s))  # if _k0s in fname else exit(1)
                                                if len(iter_pars_keys) == 10:
                                                    par_list.append(_apply_pars(pars,
                                                                                keys=[k0s, k1s, k2s, k3s, k4s, k5s, k6s, k7s, k8s, k9s],
                                                                                vals=[v0s, v1s, v2s, v3s, v4s, v5s, v6s, v7s, v8s, v9s],
                                                                                prefix_key=prefix_key, prefix=fname9))
                                                    continue

    else:
        par_list = [pars]
    return (par_list)


def read_parfile(workingdir, fname="parfile.par", comment="#", sep1="", sep2=""):
    par_sep = "* Parameters"
    opt_sep = "* Settings"
    par_var_sep = " = "
    lines = []
    pars = {}
    opts = {}
    reading_pars = False
    reading_opts = False
    with open(workingdir+fname) as file:
        lines = [line.rstrip() for line in file]
    skip = True
    for iline, line in enumerate(lines):
        if (line == sep1):
            skip = False
        if (line == sep2):
            skip = True
        if (not skip):
            if (len(line) == 1 or line == '' or line[0] == comment):
                continue
            if (line == par_sep):
                reading_pars = True
                reading_opts = False
                continue
            if (line == opt_sep):
                reading_pars = False
                reading_opts = True
                continue
            if (reading_pars):
                if (line.__contains__("#")):
                    line = line.split("#")[0]
                if (not line.__contains__(par_var_sep)):
                    raise KeyError("line does not contain par-var separator {} : {}".format(par_var_sep, line))
                name=line.split(par_var_sep)[0]
                val=line.split(par_var_sep)[-1]
                val = np.double(val)
                pars[name] = val
            if (reading_opts):
                if (line.__contains__("#")):
                    line = line.split("#")[0]
                if (not line.__contains__(par_var_sep)):
                    raise KeyError("line does not contain par-var separator {} : {}".format(par_var_sep, line))
                name=line.split(par_var_sep)[0]
                val=line.split(par_var_sep)[-1]
                opts[name] = val
            if (iline==len(lines)-1):
                return ({}, {})
    # if (len(pars.keys())==0):
    #     raise ValueError("Empty pars.dict something went wrong.")
    # if (len(opts.keys())==0):
    #     raise ValueError("empty opts.dict: something went wrong.")
    return (pars, opts)


def modify_parfile(newpars : dict, newopts : dict, workingdir, comment="#",
                   fname="parfile.par", newfname="parfile2.par", sep1="", sep2="",verbose=False):
    par_sep = "* Parameters"
    opt_sep = "* Settings"
    par_var_sep = " = "
    new_lines = []
    reading_pars = False
    reading_opts = False
    with open(workingdir+fname) as file:
        lines = [line.rstrip() for line in file]
    skip = True
    for iline, line in enumerate(lines):
        new_lines.append(line)
        if (line == sep1):
            skip = False
        if (line == sep2):
            skip = True
        if (not skip):
            if (len(line) == 1 or line == '' or line[0] == comment):
                continue
            if (line == par_sep):
                reading_pars = True
                reading_opts = False
                continue
            if (line == opt_sep):
                reading_pars = False
                reading_opts = True
                continue
            if (reading_pars):
                for key in newpars.keys():
                    if ((line.split(par_var_sep)[0]) == key):
                        _comment = ""
                        if (line.__contains__("#")):
                            _comment = " # "+ line.split(" # ")[-1]
                        new_lines[iline] = line.replace(line.split(par_var_sep)[-1],str(newpars[key])) + _comment
                        i = 1
            if (reading_opts):
                for key in newopts.keys():
                    if ((line.split(par_var_sep)[0]) == key):
                        _comment = ""
                        if (line.__contains__("#")):
                            _comment = " # "+ line.split(" # ")[-1]
                        new_lines[iline] = line.replace(line.split(par_var_sep)[-1],str(newopts[key])) + _comment
    # for line in new_lines :
    #     line += "\n"
    with open(workingdir+newfname, 'w') as f:
        for line in new_lines:
            f.write(f"{line}\n")
    if verbose:
        print("saved {}".format(workingdir+newfname))


def modify_parfile_par_opt(workingdir : str, part : str, newpars : dict, newopts : dict,
                           parfile="parfile.par",newparfile="parfile2.par",keep_old=True, verbose=False):
    if (keep_old and parfile == newparfile):
        raise NameError("Cannot keep old parfile if the new name is the same as old")
    if not (os.path.isfile(workingdir+parfile)):
        raise FileNotFoundError("parfile {} not found".format(workingdir+parfile))
    copyfile(workingdir+parfile,workingdir+"tmp_{}".format(newparfile))
    # -----------------------------------------------------------------------
    if (part=="main"):
        sep1="# -------------------------- main ---------------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="grb"):
        sep1="# ---------------------- GRB afterglow ----------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="kn"):
        sep1="# ----------------------- kN afterglow ----------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="magnetar"):
        sep1="# ------------------------ Magnetar -------------------------"
        sep2="# --------------------------- END ---------------------------"
    else:
        raise NameError("no part: {} in parfile".format(part))
    # -------------------------------------------------------------------------
    modify_parfile(newpars=newpars,newopts=newopts,workingdir=workingdir,comment="#",
                   fname="tmp_{}".format(newparfile),
                   newfname="tmp_mod_{}".format(newparfile),
                   sep1=sep1,
                   sep2=sep2,
                   verbose=verbose)
    if not keep_old:
        os.remove(workingdir+parfile)
    copyfile(workingdir+"tmp_mod_{}".format(newparfile),workingdir+newparfile)
    os.remove(workingdir+"tmp_{}".format(newparfile))
    os.remove(workingdir+"tmp_mod_{}".format(newparfile))

class Defaults:

    parfile_main_part = dict(
        pars=dict(
            tb0 = 1.e3,     # [numeric] start of the time grid (burster frame) [s]
            tb1 = 1.e12,    # [numeric] end of the time grid (burster frame) [s]
            ntb = 1000,     # [numeric] number of grid points for ODE solver
            iout = 1,       # [numeric] keep (store) solution at 'iout'th iteration
            rtol = 1e-8,    # [numeric] relative tolerance for ODE solver
            nmax = 1000,    # [numeric] maximum number of iteration for adaptive ODE solver
            A0 = -1,        # [ISM] wind environment constant, keep =0 for uniform ISM
            s = -1,         # [ISM] wind environment slope, keep =0 for uniform ISM
            r_ej = -1,      # [ISM] radius at which first break in density profile [cm]
            r_ism = -1,     # [ISM] radius at which second break in density profile [cm]
            n_ism = 1.,    # [ISM] ism number density if it is constnat [cm^-3]
            d_l = 3.09e26,  # [source] luminocity distance to the source
            z = 0.028,      # [source] redshift of the source
            theta_obs = 0   # [source] observer angle with respect to the polar axis
        ),
        opts=dict(
            do_average_solution = "no", # [numeric] store only averaged ODE result for 'n_store_substeps' (see below)
            lc_freqs = "array 1e9 1e18", # [obs] frequencies to compute light curve [Hz]
            lc_times = "array logspace 3.e3 1.e10 100", # [obs] times to compute light curve [s]
            lc_use_freq_to_time = "no", # [obs] assume freq-to-time 1-1 relation (lc=len(times)*len(freqs) otherwise)
            skymap_freqs = "array 1e9 1e18", # [obs] frequencies to compute skymap [Hz]
            skymap_times = "array logspace 8.e5 2e9 50", # [obs] times to compute skymaps [s]
            integrator = "DOP853E", # [numeric] ODE solver (DOP853E uses adaptive step-size)
        )
    )
    parfile_grb_part = dict(
        pars=dict(
            # --- Forward Shock Microphysics ---
            eps_e_fs = 0.1,
            eps_b_fs = 0.01,
            eps_t_fs = 0.,
            p_fs = 2.2,
            epsilon_e_rad = -1,    # [FS] fraction of the Esh2 removed due to radiation per timestep (dynamics)
            gamma_max_fs = 1e7,    # [numeric] Used only if 'method_gamma_max_fs=useConst'; must be < gam2
            max_substeps_fs = 2000,# [numeric] Number of cooling substeps in electron evolution (between evol.steps)
            gam1_fs = 1.,      # [numeric] lower lim for comoving electron spectrum
            gam2_fs = 1.e8,    # [numeric] upper lim for comoving electron spectrum
            ngam_fs = 451,     # [numeric] size of the electron grid points for Chang-Cooper scheme
            freq1_fs = 1.e5,   # [numeric] lower lim for comoving synchrotron spectrum
            freq2_fs = 1.e22,  # [numeric] uppers lim for comoving synchrotron spectrum
            nfreq_fs = 401,    # [numeric] size of the freq. grid points for Chang-Cooper scheme
            # --- Reverse shock Microphsyics ---
            eps_e_rs = 0.1,
            eps_b_rs = 0.01,
            eps_t_rs = 0.,
            p_rs = 2.2,
            epsilon_e_rad_rs = -1, # [RS] fraction of the Esh2 removed due to radiation per timestep (dynamics)
            gamma_max_rs = 1e7,    # [numeric] Used only if 'method_gamma_max_fs=useConst'; must be < gam2
            max_substeps_rs = 2000,# [numeric] Number of cooling substeps in electron evolution (between evol.steps)
            gam1_rs = 1.,      # [numeric] lower lim for comoving electron spectrum
            gam2_rs = 1.e7,    # [numeric] upper lim for comoving electron spectrum
            ngam_rs = 451,     # [numeric] size of the electron grid points for Chang-Cooper scheme
            freq1_rs = 1.e5,   # [numeric] lower lim for comoving synchrotron spectrum
            freq2_rs = 1.e22,  # [numeric] uppers lim for comoving synchrotron spectrum
            nfreq_rs = 401,    # [numeric] size of the freq. grid points for Chang-Cooper scheme
            # -------------------
            n_store_substeps = 10,  # use n steps of ODE solver to average over and store (used if iout >> 1)
            tprompt = 1.e3,         # [RS] duration of the ejection (for RS initial width Delta=tprompt*c)
            a = 1,                  # [spread] if method_spread="AA", controls dtheta/dR slope
            rs_shutOff_criterion_rho = 1e-50, # [RS] criterion for rho4 when to shut down the reverse shock
            min_Gamma0_for_rs=0.,   # [RS] If initial Gamma0 of a BW (layer) < this value, use 'fs' RHS not 'fsrs'
            mom0_frac_when_start_spread = 0.9, # [spread] frac, when to allow spread, \Gamma\beta < frac * Gamma\beta_0
            rs_Gamma0_frac_no_exceed = .92, # [RS] if Gamma > frac*Gamma0; set dGammadR = 0 (prevent error acceleration)
            save_dyn_every_it = 10, # [numeric] if to save dynamics, save every it'th iteration,
            rtol_phi = 1e-6,        # [eats] relative tolerance for adaptive quadrature for EATS integration
            rtol_theta = 1e-6,      # [eats] relative tolerance for adaptive quadrature for EATS integration
            nmax_phi = 1000,        # [eats] maximum number of adaptive quadrature iterations to find solution
            nmax_theta = 1000,      # [eats] maximum number of adaptive quadrature iterations to find solution
            theta_max = np.pi/2.,   # [eats] maximum extend of each hemispheres
            beta_min_fs = 1e-5,     # [numeric] if < betaShock; do not compute any microphsyics (FS)
            beta_min_rs = 1e-5,     # [numeric] if < betaShock; do not compute any microphsyics (RS)
            # --- Skymap adaptive calculation; resize untill 'min_sublayers' each of which has 'min_non_zero_cells'
            nsublayers = 10,        # [numeric] initial division of a theta-layer into sublayers (may give I=0 cells -> 'redo')
            frac_to_increase = 1.5, # [numeric] frac. to increase nsublayers id skymap is not resolved (see below)
            max_restarts = 10,      # [numeric] max number of increasing 'nsublayers' to resolve the jet
            min_sublayers = 3,      # [numeric] criterion, min number sublayers each of which has 'min_non_zero_cells'
            min_non_zero_cells = 3, # [numeric] criterion, min number of phi-cells with intensity > 0 for jet to be resolved
            im_max_theta = 1.5708   # [numeric] max value in intensity calculations for skymap
        ),
        opts=dict(
            run_bws = "yes",        # [task] evolve the blastwaves (if no, expected that load_dynamics=yes)
            save_dynamics = "no",   # [task] save blastwaves evolution history
            load_dynamics = "no",   # [task] load the blastwave dynamics from a file
            do_mphys_in_situ= "yes",# [task] compute electrons/syc/ssc comov.spec after each it. of BW evolution
            do_mphys_in_ppr = "no", # [task] compute electrons/syc/ssc comov.spec after full BW evolution has finisheds
            # do_ele = "yes",         # [task] compute comoving electron spectrum (numerically or analytically)
            # do_spec = "no",         # [task] compute comoving synchrotron/SSC spectrum (numerically or analytically)
            save_spec = "no",       # [task] save comoving ele/synchrotron/SSC spectrum (numerically or analytically)
            do_lc = "yes",          # [task] compute & save light curves
            do_skymap = "no",       # [task] compute & save raw skymaps
            save_raw_skymap="yes",  # [task] currently only raw, unstructured images can be saved. (UNFINISHED)
            skymap_remove_mu = "no",# [task] remove 'mu' from skymap calculation
            counter_jet = "yes",    # [numeric] do include counter jet as well in LCs and Sky Maps

            do_rs = "no",               # [RS] include RS into consideration (main switch)
            do_rs_radiation="yes",      # [RS] if RS is included, compute also the radiation from RS (adds to total LC)
            bw_type = "fs",             # [numeric] type pf the blastwave RHS to use, e.g. fs - forward shock only
            init_deltaR4="no",          # [numeric] set deltaR4[0] = c*beta0*tprompt for ODE solver
            exponential_rho4="yes",     # [numeric] use exponential rho4 decay as exp(-Delta4/Delta0)
            method_collision = "none",  # [numeric] include blastwave collision (UNFINISHED)
            method_eats = "adaptive",   # [numeric] main switch for blastwave discretezation (piece-wise or adaptive)
            method_quad = "CADRE",      # [numeric] EATS quadrature method (CADRE is adaptive)
            method_comp_mode = "comovSpec", # [numeric] interpolated comoving spectra, or compute in-situe
            allow_termination = "no",   # [numerc] continue if one of the blastwaves fails (ODE solver fails)

            do_thermrad_loss = "no",    # [numeric] include thermal radiation from ejecta (UNFINISHED)
            do_eninj_inside_rhs = "no", # [numeric] magnetar-driven ejecta; (UNFINISHED)

            use_1d_id = "yes",          # [I/O] type of the initail data, if 'yes' expects 1D arrays with E,Gamma...
            fname_ejecta_id = "id.h5",  # [I/O] file name (in working_dir) with initial data
            load_r0 = "no",             # [I/O] use R0 from the file instead of computing it as R0=beta0 * tb0 * c
            fname_dyn = "dyn.h5",       # [I/O] file name (in working_dir) to save dynamics
            fname_spectrum = "spec.h5", # [I/O] file name (in working_dir) to save the comoving ele/radition spectrum
            fname_light_curve = "lc.h5",# [I/O] file name (in working_dir) to save light curve
            fname_sky_map = "skymap.h5",# [I/O] file name (in working_dir) to save raw light curve

            ebl_tbl_fpath = "../../../data/EBL/Franceschini18/table.h5", # if not "none" use F_nu = F_nu*exp(-tau(z,nu))

            do_nucinj = "no", # [numeric] include r-process heating in ejecta (UNFINISHED)

            method_spread = "our",              # [spread] method for lateral spreading
            method_limit_spread = "Mom0Frac",   # [numeric] how to limit spreading of the blastwave
            method_dgdr = "our",                # [numeric] choice of equation for BW dynamical evolution dGamma/dR
            method_eos = "Nava13",              # [numeric] choice of EOS for the blast wave
            method_dmdr = "usingdthdr",         # [numeric] choice of equation for accreted mass dm/dr

            use_dens_prof_behind_jet_for_ejecta = "no", # [numeric] include jet in ejecta mode (UNFINISHED)

            # --- Forward Shock ---
            # method_radius_fs = "useGammaR",       # [numeric] how to ge radius for forward shock (use GammaShock or not)
            use_adiabLoss = "yes",              # [numeric] include blast wave adiabatic lossess (FS)
            method_Gamma_fs = "useGammaShock",  # [numeric] compute GammaShock via EOS or assume = to Gamma
            method_Up_fs = "useEint2",          # [numeric] compute internal energy from Eint2 or Gamma
            method_thickness_fs = "useJoh06",   # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
            method_vel_fs = "sameAsBW",         # [numeric] "shockVel" in EATS, compute abberation using GammaShock or Gamma
            method_ele_fs = "numeric",          # [numeric] assume analytical electron profile or evolve
            num_ele_use_adi_loss_fs="yes",      # [numeric] include adiabatic cooling term into kinetic eq. for ele. evolution
            method_ne_fs = "useNe",             # [numeric] compute emissivities using Ne or nprime
            method_nonrel_dist_fs ="use_Sironi",# [numeric] include Deep Newtonian regime for electron dist.
            method_gamma_min_fs = "useNumeric", # [numeric] how to compute gamma_min
            method_gamma_c_fs = "useTcomov",    # [numeric] how to compute gamma_c
            method_gamma_max_fs = "useB",       # [numeric] how to compute gamma_max
            method_B_fs = "useU_b" ,            # [numeric] how to compute magnetic field
            method_synchrotron_fs = "CSYN",     # [numeric] method for the synchrotron radiation
            use_ssa_fs = "no",                  # [numeric] include SSA
            method_ssc_fs = "none",             # [numeric] method for SSC
            method_pp_fs = "none",              # [numeric] method for pair-production
            method_tau_fs = "smooth",           # [numeric] method for optical depth calculation in shock

            # --- Reverse Shock ---
            # method_radius_rs = "sameAsR",     # [numeric] how to ge radius for reverse shock (use GammaShock or not)
            use_adiabLoss_rs = "yes",           # [numeric] include blast wave adiabatic lossess (RS)
            method_Gamma_rs = "useGammaShock",  # [numeric] compute GammaShock via EOS or assume = to Gamma
            method_Up_rs = "useEint2",          # [numeric] compute internal energy from Eint2 or Gamma
            method_thickness_rs = "useJoh06",   # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
            method_vel_rs = "sameAsBW",         # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
            method_ele_rs = "numeric",          # [numeric] assume analytical electron profile or evolve
            num_ele_use_adi_loss_rs="yes",      # [numeric] include adiabatic cooling term into kinetic eq. for ele. evolution
            method_ne_rs = "useNe",             # [numeric] compute emissivities using Ne or nprime
            method_nonrel_dist_rs ="use_Sironi",# [numeric] include Deep Newtonian regime for electron dist.
            method_gamma_min_rs = "useNumeric", # [numeric] how to compute gamma_min
            method_gamma_c_rs = "useTcomov",    # [numeric] how to compute gamma_c
            method_gamma_max_rs = "useB",       # [numeric] how to compute gamma_max
            method_B_rs = "useU_b",             # [numeric] how to compute magnetic field
            method_synchrotron_rs = "CSYN",     # [numeric] method for the synchrotron radiation
            use_ssa_rs = "no",                  # [numeric] include SSA
            method_ssc_rs = "none",             # [numeric] method for SSC
            method_pp_rs = "none",              # [numeric] method for pair-production
            method_tau_rs = "smooth",           # [numeric] method for optical depth calculation in shock
        )
    )

    parfile_kn_part = dict(
        pars = dict(),
        opts = dict()
    )

    parfile_mag_part = dict(
        pars = dict(),
        opts = dict()
    )

def _create_parfile_part(lines:list,part:str,sep1:str,sep2:str,default:dict,new:dict):
    sep_pars = "* Parameters"
    sep_opts = "* Settings"
    # update main parameters if they are in the 'P'
    if part in new.keys():
        # update if needed
        for k, v in new[part].items():
            # check parameters
            if k in default["pars"].keys():
                default["pars"][k] = v
            elif k in default["opts"].keys():
                default["opts"][k] = v
            else:
                raise KeyError("key = {} is not in the {} parameters: {} \n or  options: {} \n"
                               .format(k,part,default["pars"].keys(), default["opts"].keys()))
    # Creat Main Section of the Parfile
    lines.append("\n"); lines.append(sep1)
    lines.append("\n"); lines.append(sep_pars); lines.append("\n")
    for k, v in default["pars"].items():
        lines.append(f"{k} = {v}")
    lines.append("\n"); lines.append(sep_opts); lines.append("\n")
    for k, v in default["opts"].items():
        lines.append(f"{k} = {v}")
    lines.append("\n"); lines.append(sep2); lines.append("\n")
    return lines
def create_parfile(working_dir : str, P : dict):
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    lines = []
    lines = _create_parfile_part(lines=lines,part="main",
                         sep1="# -------------------------- main ---------------------------",
                         sep2="# --------------------------- END ---------------------------",
                         default=copy.deepcopy(Defaults.parfile_main_part),
                         new=P)

    if ("grb" in P.keys()):
        lines = _create_parfile_part(
            lines=lines,part="grb",
            sep1="# ---------------------- GRB afterglow ----------------------",
            sep2="# --------------------------- END ---------------------------",
            default=copy.deepcopy(Defaults.parfile_grb_part),
            new=P
        )

    if ("kn" in P.keys()):
        lines = _create_parfile_part(
            lines=lines,part="kn",
            sep1="# ----------------------- kN afterglow ----------------------",
            sep2="# --------------------------- END ---------------------------",
            default=copy.deepcopy(Defaults.parfile_kn_part),
            new=P
        )

    if ("magnetar" in P.keys()):
        _create_parfile_part(
            lines=lines,part="magnetar",
            sep1="# ------------------------ Magnetar -------------------------",
            sep2="# --------------------------- END ---------------------------",
            default=copy.deepcopy(Defaults.parfile_mag_part),
            new=P
        )

    with open(working_dir+"parfile.par", 'w') as f:
        for line in lines:
            f.write(f"{line}\n")
