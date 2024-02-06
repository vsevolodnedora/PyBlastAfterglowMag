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
