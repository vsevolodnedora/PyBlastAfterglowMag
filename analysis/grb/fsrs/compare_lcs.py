import package.src.PyBlastAfterglowMag as PBA
import os,shutil,matplotlib.pyplot as plt
import numpy as np
working_dir = os.getcwd()+'/tmp1/'
fig_dir = os.getcwd()+'/figs/'

def mrg(dict2:dict, dict1:dict):
    return {k: v | dict2[k] for k, v in dict1.items() if k in dict2.keys()}

def run(working_dir:str, struct:dict, P:dict, type:str="a", run:bool=True) -> PBA.PyBlastAfterglow:
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
    if run:
        pba.run(
            path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
            loglevel="info"
        )

    # process skymap
    if (run and pba.GRB.opts["do_skymap"]=="yes"):
        conf = {"nx":128, "ny":64, "extend_grid":1.1, "fwhm_fac":0.5, "lat_dist_method":"integ",
                "intp_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }, # "gaussian"
                "hist_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }}
        prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
        prep.process_singles(infpaths=working_dir+"raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=True)

    return pba

def plot(struct:dict,pp:dict,plot:dict):


    fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(4.6,3.2))

    pba = run(working_dir=working_dir,struct=struct,P=mrg(dict(
        main=dict(),
        grb=dict(do_rs_radiation='no',do_rs='no',bw_type='fs')),pp),type="a",run=True)
    ax.plot(
        pba.GRB.get_lc_times(unique=True,spec=False),
        pba.GRB.get_lc(key="fluxdens", xkey="freqs", key_time="times", freq=1e9, time=None,
                       ishell=0, ilayer=0, sum_shells_layers=True, spec=False),
        color='blue',label='FS'

    )

    pba = run(working_dir=working_dir,struct=struct,P=mrg(dict(
        main=dict(),
        grb=dict(do_rs_radiation='yes',do_rs='yes',bw_type='fsrs')),pp),type="a",run=True)
    ax.plot(
        pba.GRB.get_lc_times(unique=True,spec=False),
        pba.GRB.get_lc(key="fluxdens", xkey="freqs", key_time="times", freq=1e9, time=None,
                       ishell=0, ilayer=0, sum_shells_layers=True, spec=False),
        color='green',label='FS \& RS'
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"$F_{\rm nu}$")
    ax.set_xlabel(r"$t_{\rm obs}$")
    plt.tight_layout()
    plt.legend()
    plt.show()

def plot_burster(struct:dict,pp:dict,plot:dict):

    v_ns = (
        dict(v_ns=("rho2","rho3"),ylabel=r'$\rho$'),
        dict(v_ns=("U_p", "U_p3"),ylabel=r'$U_{\rm p}$',ylim=(1e-1,1e7)),
        dict(v_ns=("Gamma", "Gamma43"),ylabel='$\Gamma$',ylim=(1e-1,1e3)),
        # ("GammaFsh", "GammaRsh"),
        dict(v_ns=("gamma_min","gamma_min_rs"),ylabel=r"$\gamma_{\rm m}$",ylim=(1e-1,2e4)),
        dict(v_ns=("accel_frac","accel_frac_rs"),ylabel=r"$\xi_{\rm DN}$",ylim=(1e-3,2)),
        dict(v_ns=("gamma_c","gamma_c_rs"),ylabel=r"$\gamma_{\rm c}$",ylim=(1.,1e7)),
    )
    xlim=(3e3,1e9)

    fig,axes = plt.subplots(ncols=1,nrows=len(v_ns),figsize=(5,8),sharex='all')

    pba = run(working_dir=working_dir,struct=struct,P=mrg(dict(
        main=dict(),
        grb=dict(do_rs_radiation='no',do_rs='no',bw_type='fs')),pp),type="a",run=True)
    for i, d in enumerate(v_ns):
        ax = axes[i]
        ax.plot(
            pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=0),
            pba.GRB.get_dyn_arr(v_n=d["v_ns"][0],ishell=0,ilayer=0),
            color='blue',label='FS', ls='-'
        )

    pba = run(working_dir=working_dir,struct=struct,P=mrg(dict(
        main=dict(),
        grb=dict(do_rs_radiation='yes',do_rs='yes',bw_type='fsrs')),pp),type="a",run=True)
    for i, d in enumerate(v_ns):
        ax = axes[i]
        ax.plot(
            pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=0),
            pba.GRB.get_dyn_arr(v_n=d["v_ns"][0],ishell=0,ilayer=0),
            color='blue',label='FS \& RS', ls=':', lw=1.5
        )
        ax.plot(
            pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=0),
            pba.GRB.get_dyn_arr(v_n=d["v_ns"][1],ishell=0,ilayer=0),
            color='green',label='RS', ls='-'
        )

    for i, d in enumerate(v_ns):
        ax = axes[i]
        ax.set_ylabel(d['ylabel'])
        ax.set_xscale("log")
        ax.set_yscale("log")
        if 'ylim' in d.keys():
            ax.set_ylim(*d['ylim'])
        ax.legend()
    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]")
    axes[-1].set_xlim(*xlim)
    # ax.set_ylabel(r"$F_{\rm nu}$")
    # ax.set_xlabel(r"$t_{\rm obs}$")
    plt.tight_layout()
    # plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_burster(struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
         pp = dict(main=dict(n_ism=1, tb0=3e3, ntb=1000,rtol=1e-7,theta_obs=0),
                   grb=dict(
                       # method_ele_fs='analytic',method_ele_rs='analytic',method_synchrotron_fs='WSPN99',
                       save_dynamics='yes'
                   #method_spread='None'
                   )),
         plot=dict())