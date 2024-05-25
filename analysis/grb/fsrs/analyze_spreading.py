import package.src.PyBlastAfterglowMag as PBA
import os,shutil,copy,matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm
from matplotlib import cm
import numpy as np
working_dir = os.getcwd()+'/tmp1/'
fig_dir = os.getcwd()+'/figs/'

# def run(working_dir:str, struct:dict, P:dict, type:str="a", run:bool=True) -> PBA.PyBlastAfterglow:
#     # clean he temporary direcotry
#     if run and os.path.isdir(working_dir):
#         shutil.rmtree(working_dir)
#     if not os.path.isdir(working_dir):
#         os.mkdir(working_dir)
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
#     if run:
#         pba.run(
#             path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
#             loglevel="info"
#         )
#
#     # process skymap
#     if (run and pba.GRB.opts["do_skymap"]=="yes"):
#         conf = {"nx":128, "ny":64, "extend_grid":1.1, "fwhm_fac":0.5, "lat_dist_method":"integ",
#                 "intp_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }, # "gaussian"
#                 "hist_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }}
#         prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
#         prep.process_singles(infpaths=working_dir+"raw_skymap_*.h5",
#                              outfpath=pba.GRB.fpath_sky_map,
#                              remove_input=True)
#
#     return pba
def d2d(default:dict,new:dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_


def plot_dyn(ax, ej:PBA.Ejecta, v_n_x : str, v_n_y : str, layers=(), plot_layer = {}):
    cmap = plot_layer["cmap"] if "cmap" in plot_layer.keys() else 'viridis'
    cmap = cm.get_cmap(cmap)
    llayers = ej.get_dyn_obj().attrs["nlayers"]
    vmin = plot_layer["vmin"] if "vmin" in plot_layer.keys() else 0
    vmax = plot_layer["vmax"] if "vmax" in plot_layer.keys() else int(ej.get_dyn_obj().attrs["nlayers"])
    norm = Normalize(vmin=vmin, vmax=vmax)
    for il in range(int(ej.get_dyn_obj().attrs["nlayers"])):
        if (((len(layers) > 0) and (il in layers)) or (len(layers)==0)):
            ls = plot_layer["ls"] if "ls" in plot_layer.keys() else '-'
            lw = plot_layer["lw"] if "lw" in plot_layer.keys() else 1.
            color = plot_layer["color"] if "color" in plot_layer.keys() else 'red'
            color = cmap(norm(il)) if "cmap" in plot_layer.keys() else color
            alpha = plot_layer["alpha"] if "alpha" in plot_layer.keys() else 0.9
            label = plot_layer["label"] if "label" in plot_layer.keys() else f"layer={il}"
            ax.plot(ej.get_dyn_arr(v_n=v_n_x,ishell=0,ilayer=il),
                    ej.get_dyn_arr(v_n=v_n_y,ishell=0,ilayer=il), lw=lw, ls=ls, color=color, alpha=alpha, label=label)

if __name__ == '__main__':
    figname = "gauss_spread_comparison"
    run = False
    # struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618)
    struct = dict(struct="gaussian",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.3)
    P=dict(main=dict(n_ism = 1., tb0=1e4, ntb=4000,rtol=1e-6,
                     lc_freqs = "array 1e9 1e18"),
           grb=dict(save_dynamics='yes',do_rs='no',bw_type='fs',do_mphys_in_situ = "no",do_lc = "no",
                    do_rs_radiation="no",
                    # method_spread='None'
                    # exponential_rho4='no'
                    )
           )
    pba_default = PBA.wrappers.run_grb(
        working_dir=os.getcwd()+"/working_dirs/dyn_spread_default/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict())),
        type="a",run=run
    )
    # pba_gp12 = PBA.wrappers.run_grb(
    #     working_dir=os.getcwd()+"/working_dirs/dyn_spread_gp12/",
    #     struct=struct, P = d2d(default=P, new=dict(grb=dict(method_spread='AA',a=1))),
    #     type="a",run=run
    # )
    pba_hu00 = PBA.wrappers.run_grb(
        working_dir=os.getcwd()+"/working_dirs/dyn_spread_hu00/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(method_spread='Adi'))),
        type="a",run=run
    )
    pba_afgpy = PBA.wrappers.run_grb(
        working_dir=os.getcwd()+"/working_dirs/dyn_spread_afgpy/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(method_spread='AFGPY'))),
        type="a",run=run
    )


    fig,ax = plt.subplots(ncols=1,nrows=1,layout='constrained',figsize=(5,2.5))
    plot_dyn(ax=ax,ej=pba_default.GRB,v_n_x="R",v_n_y="theta",layers=(0,10,15,20),
             plot_layer=dict(cmap="jet",ls='-'))
    # plot_dyn(ax=ax,ej=pba_gp12.GRB,v_n_x="R",v_n_y="theta",layers=(0,6,12,19),
    #          plot_layer=dict(cmap="jet",ls='--'))
    plot_dyn(ax=ax,ej=pba_hu00.GRB,v_n_x="R",v_n_y="theta",layers=(0,10,15,20),
             plot_layer=dict(cmap="jet",ls=':',lw=1.5,label=None))
    plot_dyn(ax=ax,ej=pba_afgpy.GRB,v_n_x="R",v_n_y="theta",layers=(0,10,15,20),
             plot_layer=dict(cmap="jet",ls='-.',lw=2.,label=None))
    ax.set_xscale("log")
    ax.set_xscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=False, loc='upper left', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=1, fontsize=12, labelcolor='black',
              framealpha=0.0, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_ylabel("$\omega$ [rad]", fontsize=12)
    # ax.set_title(task["title"], fontsize=12)
    ax.set_xlim(1e17,4e19)
    ax.grid(ls=':')
    ax.set_xlabel("$R$ [cm]", fontsize=12)
    plt.savefig(os.getcwd() + '/figs/' + figname + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + figname + '.pdf')
    plt.show()
    plt.close(fig)


