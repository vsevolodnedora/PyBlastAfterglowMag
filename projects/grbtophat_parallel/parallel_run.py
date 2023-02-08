#! python3
import numpy as np
import os
import copy

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import latex_float, cgs
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")

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
    del_idx = np.array([], dtype=int)
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
            del_idx = np.append(del_idx, i)
        # jet core cannot be wider than wings. Remove incorrect settings
        if (float(new_grb_parts[i]["theta_w"]) < float(new_grb_parts[i]["theta_c"])):
            del_idx = np.append(del_idx, i)

    # check if there is anything to run
    if (len(del_idx) == len(parfile_fname_list)):
        print("No parfiles left to run. Exiting...")
        exit(0)

    print("Removing incorrect/existing parfiles N={}/{}".format(len(del_idx),len(new_main_parts)))
    # print("Creating parfiles for remaining settings {}".format(len(new_main_parts)))

    clean_parfile_list = []
    for i in range(len(new_main_parts)):
        if (int(i) in del_idx):
            continue

        clean_parfile_list.append(parfile_fname_list[i])
        if not(os.path.isfile(parfile_fname_list[i])):
            print("Creating {}".format(parfile_fname_list[i]))
            modify_parfile_par_opt(part="main", newpars=new_main_parts[i], newopts={},
                                   workingdir=working_dir, parfile="parfile.par", newparfile=parfile_fname_list[i],keep_old=True)
            modify_parfile_par_opt(part="grb", newpars=new_grb_parts[i], newopts=new_grb_opts[i],
                                   workingdir=working_dir, parfile=parfile_fname_list[i], newparfile=parfile_fname_list[i],keep_old=False)

    # run PyBlastAfterglow for all parfiles needed
    print("Parfiles collected N={}/{} Starting runs...".format(len(clean_parfile_list),len(parfile_fname_list)))
    if len(parfile_fname_list) > 1:
        distribute_and_run(working_dir=working_dir, list_parfiles=clean_parfile_list, n_cpu=n_cpu)
    print("Runs finished successfully")


    # for parfile in parfile_fname_list:
    #     pba = PBA(workingdir=working_dir,readparfileforpaths=True,parfile=parfile)
    #     plt.loglog(pba.get_jet_lc_times(),pba.get_jet_lc_totalflux(1e9))
    # plt.show()

def main():
    run_list_grbafg_parallel(n_cpu=12)

if __name__ == '__main__':
    main()