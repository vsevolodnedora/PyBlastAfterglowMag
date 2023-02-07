#! python3
import numpy as np
import os
import copy

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import distribute_and_run
    from PyBlastAfterglowMag.utils import latex_float, cgs
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import distribute_and_run
        from package.src.PyBlastAfterglowMag.utils import latex_float, cgs
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

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
def _apply_pars(pars, keys : list, vals : list, prefix_key : str, prefix : str):
    _pars = copy.deepcopy(pars)
    for (key, val) in zip(keys, vals):
        _pars[key] = val
    _pars[prefix_key] = prefix
    return _pars
def set_parlists_for_pars(iter_pars_keys, iter_pars : dict, fname : str):
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
        print("iterating over {} pars \n{}".format(n_iter_pars, iter_pars_keys))

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

                                                # _pars = copy.deepcopy(pars)
                                                # _pars[k0s] = v0s
                                                # _pars[k1s] = v1s
                                                # _pars[k2s] = v2s
                                                # _pars[k3s] = v3s
                                                # _pars[k4s] = v4s
                                                # _pars[k5s] = v5s
                                                # _pars[k6s] = v6s
                                                # _pars[k7s] = v7s
                                                # _pars[k8s] = v8s
                                                # _pars[k9s] = v9s
                                                # _pars[prefix_key] = fname9
                                                # par_list.append(_pars)
                                                # print(fname5)
    else:
        par_list = [pars]
    return (par_list)
def remove_theta_w_less_theta_c(parlist : list):
    n1 = len(parlist)
    for i, par in enumerate(parlist):
        if float(par["theta_w"]) < float(par["theta_c"]):
            parlist.pop(i)
    n2 = len(parlist)
    print("Removed {} settings for theta_w < theta_c".format(n1-n2))
    return parlist

def run_list_grbafg_parallel(n_cpu=12):
    working_dir = os.getcwd()+'/'
    iter_pars = ["Eiso_c", "Gamma0c", "theta_c", "theta_w", "p", "eps_e", "eps_b","n_ism","theta_obs"]
    iter_pars_dict = {
        "n_ism": [1.0, 0.1, 0.01, 0.001],
        "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
        "Eiso_c": [1.e50, 1.e51, 1.e52, 1.e53],
        "Gamma0c":[100., 300., 600., 1000.],
        "theta_c": np.array([5., 10., 15., 20.]) * np.pi / 180.,
        "theta_w": np.array([5., 10., 15., 20.]) * np.pi / 180.,
        "p": [2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
        "eps_e": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        "eps_b": [0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
    }
    # iter_pars_dict = {
    #     "n_ism": [1.0],#[1.0, 0.1, 0.01, 0.001],
    #     "theta_obs": np.array([0., 15., 45.0, 60., 75., 90.]) * np.pi / 180.0,  # [75*np.pi/180]#
    #     "Eiso_c": [1e51],#[1e50, 1e51, 1e52, 1e53],
    #     "Gamma0c":[100.],#[100., 300, 600, 1000],
    #     "theta_c": np.array([5.]) * np.pi / 180.,#np.array([5., 10., 15., 20.]) * np.pi / 180.,
    #     "theta_w": np.array([5.]) * np.pi / 180.,#np.array([5., 10., 15., 20.]) * np.pi / 180.,
    #     "p": [2.6],#[2.2, 2.4, 2.6, 2.8],  # [2.05, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0],
    #     "eps_e": [0.1],#[0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
    #     "eps_b": [0.001]#[0.5, 0.1, 0.01, 0.001],  # [0.001, 0.005, 0.01, 0.05, 0.1],
    # }

    # a name of the parfile and light curve will be related to the actial parameters vaied
    newparfilenameroot = "EisoC[Eiso_c]_Gamma0c[Gamma0c]_thetaC[theta_c]_thetaW[theta_w]_" \
                         "theta[theta_obs]_nism[n_ism]_p[p]_epse[eps_e]_epsb[eps_b]_"
    # separate dict is needed for pars and opts; so we make two copies
    new_main_parts = set_parlists_for_pars(iter_pars_keys=iter_pars,
                                         iter_pars=iter_pars_dict, fname=newparfilenameroot)
    # new_main_parts = remove_theta_w_less_theta_c(new_main_parts)
    new_grb_parts = set_parlists_for_pars(iter_pars_keys=iter_pars,
                                         iter_pars=iter_pars_dict, fname=newparfilenameroot)
    # new_grb_parts = remove_theta_w_less_theta_c(new_grb_parts)
    # the light curve name is actually an option, so should be in a separate dict...
    new_grb_opts = [{"fname_light_curve":""} for i in range(len(new_grb_parts))]

    # now we need to remove grb pars from main pars dict and vise versa
    main_keys = ["n_ism","theta_obs"]
    grb_keys = ["Eiso_c", "Gamma0c", "theta_c", "theta_w", "p", "eps_e", "eps_b"]
    # and collect the parfile names to pass them later ito parallel executor/drive
    parfile_fname_list = []
    for i in range(len(new_main_parts)):
        for key in grb_keys:
            new_main_parts[i].pop(key)
        for key in main_keys:
            new_grb_parts[i].pop(key)
        # NOTE this is where the light curve name is made. Tophat is assumed.
        new_grb_opts[i]["fname_light_curve"] = "tophat_"+new_main_parts[i]["name"]+"lc.h5"
        parfile_fname_list.append("tophat_"+new_main_parts[i]["name"]+"parfile.par")
        # if parfile already exists, assume that the calculation has already been done, so skip it
        if (os.path.isfile(new_grb_opts[i]["fname_light_curve"])):
            del new_main_parts[i]
            del new_grb_parts[i]
            del new_grb_opts[i]
            del parfile_fname_list[i]
            print("Skipping {}".format(new_grb_opts[i]["fname_light_curve"]))
        # jet core cannot be wider than wings. Remove incorrect settings
        if (float(new_grb_opts[i]["theta_w"]) < float(new_grb_opts[i]["theta_c"])):
            del new_main_parts[i]
            del new_grb_parts[i]
            del new_grb_opts[i]
            del parfile_fname_list[i]
            print("Skipping incorrect {}".format(new_grb_opts[i]["fname_light_curve"]))
        # if required parfile name is not found, make one by modifying the 'default' parfile.par
        if not(os.path.isfile(parfile_fname_list[i])):
            print("Creating {}".format(parfile_fname_list[i]))
            modify_parfile_par_opt(part="main", newpars=new_main_parts[i], newopts={},
                                   workingdir=working_dir, parfile="parfile.par", newparfile=parfile_fname_list[i],keep_old=True)
            modify_parfile_par_opt(part="grb", newpars=new_grb_parts[i], newopts=new_grb_opts[i],
                                   workingdir=working_dir, parfile=parfile_fname_list[i], newparfile=parfile_fname_list[i],keep_old=False)
    # run PyBlastAfterglow for all parfiles needed
    print("Parfiles created. Starting runs...")
    if len(parfile_fname_list) > 1:
        distribute_and_run(working_dir=working_dir, list_parfiles=parfile_fname_list, n_cpu=n_cpu)
    print("Runs finished successfully")


    # for parfile in parfile_fname_list:
    #     pba = PBA(workingdir=working_dir,readparfileforpaths=True,parfile=parfile)
    #     plt.loglog(pba.get_jet_lc_times(),pba.get_jet_lc_totalflux(1e9))
    # plt.show()

def main():
    run_list_grbafg_parallel(n_cpu=12)

if __name__ == '__main__':
    main()