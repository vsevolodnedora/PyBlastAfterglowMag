import copy

import package.src.PyBlastAfterglowMag as PBA
import os,shutil,matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize,TwoSlopeNorm,SymLogNorm
import numpy as np
working_dir = os.getcwd()+'/tmp1/'
fig_dir = os.getcwd()+'/figs/'


def d2d(default:dict,new:dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_

def run(working_dir:str, struct:dict, P:dict, type:str="a", run:bool=True) -> PBA.PyBlastAfterglow:
    # clean he temporary direcotry
    if run and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
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

def _normalize_spec(spec:np.ndarray, xs:np.ndarray, ys:np.ndarray, norm_method:None or str, mask_negatives=True):
    if not norm_method:
        return spec
    if norm_method == "/integ *y^2":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        spec *= np.power(ys, 2)[np.newaxis,:]
    elif norm_method == "/integ *y":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        spec *= np.power(ys, 1)[np.newaxis,:]
    elif norm_method == "*y":
        spec *= ys[np.newaxis,:]
    elif norm_method == "*y^2":
        spec *= ys[np.newaxis,:] ** 2
    elif norm_method == "/integ":
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
    else:
        raise KeyError("Norm method is not recognized")
    if mask_negatives:
        spec[~np.isfinite(spec)] = 1e-100
        spec[spec <= 0] = 1e-100
    return spec

def plot_compare_two_spectra_new(ej:PBA.Ejecta, ej_an:PBA.Ejecta,
                                 ej_sustract:PBA.Ejecta or None = None,
                                 ej_an_sustract:PBA.Ejecta or None = None,
                                 name1:str= "Numeric", name2="Analytic", is_spec=True,
                                 fs_or_rs :str = "fs", ele_syn_ssc:str = "ele", norm_method:str or None="integ_gam2",
                                 xkey="times_gams", ykey="gams",
                                 task=dict()):
    mp = 1.6726e-24
    qe = 4.803204e-10
    me = 9.1094e-28
    c  = 2.99792458e10


    # plt.loglog(ej.get_dyn_arr(v_n='B_rs',ishell=0,ilayer=0))
    # plt.loglog(ej.get_dyn_arr(v_n='B',ishell=0,ilayer=0))
    # plt.show()



    fig, axes = plt.subplots(**task['subplots'])
    iplot = 0

    # -------------------------------------------

    spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                   xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                   spec=is_spec,sum_shells_layers=False)
    xs = ej.get_grid(key=xkey,spec=is_spec)
    ys = ej.get_grid(key=ykey,spec=is_spec)
    if fs_or_rs == "rs" and not ej_sustract is None:
        spec_sub=ej_sustract.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                                    xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                                    spec=is_spec,sum_shells_layers=False)
        xs_ = ej_sustract.get_grid(key=xkey,spec=is_spec)
        ys_ = ej_sustract.get_grid(key=ykey,spec=is_spec)
        if not np.array_equal(xs,xs_):
            raise ValueError("grid mismatch")
        if not np.array_equal(ys,ys_):
            raise ValueError("grid mismatch")
        spec -= spec_sub
    spec = _normalize_spec(spec=spec,xs=xs,ys=ys,norm_method=norm_method,mask_negatives=True)

    spec_an=ej_an.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                         xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                         spec=is_spec,sum_shells_layers=False)
    xs_an = ej_an.get_grid(key=xkey,spec=is_spec)
    ys_an = ej_an.get_grid(key=ykey,spec=is_spec)
    if fs_or_rs == "rs" and not ej_sustract is None:
        spec_sub_an=ej_an_sustract.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                                          xkey=xkey,ykey=ykey,freq=None,time=None,ishell=0,ilayer=0,
                                          spec=is_spec,sum_shells_layers=False)
        xs_an_ = ej_an_sustract.get_grid(key=xkey,spec=is_spec)
        ys_an_ = ej_an_sustract.get_grid(key=ykey,spec=is_spec)
        if not np.array_equal(xs_an,xs_an_):
            raise ValueError("grid mismatch")
        if not np.array_equal(ys_an,ys_an_):
            raise ValueError("grid mismatch")
        spec_an -= spec_sub_an
    spec_an = _normalize_spec(spec=spec_an,xs=xs_an,ys=ys_an,norm_method=norm_method,mask_negatives=True)

    # --------------------------------------------

    ax = axes[iplot] #axes[1,1]
    ax.set_title(task["title"],fontsize=12)

    colors= task["colors_ys"]
    indexes_ys = [PBA.utils.find_nearest_index(ys,freq) for freq in task["ys"]]
    lines = []
    for idx, color in zip(indexes_ys, colors):
        vals = spec[:, idx]
        lns, = ax.plot(xs, vals, color=color, linewidth=1.0, linestyle="-")
        lines.append(lns)
        vals_an = spec_an[:, idx]
        ax.plot(xs_an, vals_an, color=color, linewidth=1.0, linestyle="--")

    legend1 = ax.legend(lines, [task["ylabel"]+f"=${PBA.utils.latex_float(ye,format='{0:.0e}')}$" for ye in task["ys"]],
                        loc='upper right',fontsize=11,fancybox=False,shadow=False,framealpha=0., borderaxespad=0.,
                        ncols=3,columnspacing=0.8)
    lin_, = ax.plot([0,0],[1,1],color='gray',ls='-')
    lin2_, = ax.plot([0,0],[1,1],color='gray',ls='--')
    legend2 = ax.legend([lin_,lin2_], [name1,name2],
                        loc='lower right',fontsize=11,fancybox=False,shadow=False,framealpha=0., borderaxespad=0.,
                        ncols=2,columnspacing=0.8)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    if not task["ylim"]:
        ax.set_ylim(np.max(spec[:, indexes_ys])*1e-7,
                    np.max(spec[:, indexes_ys])*10)
    else:
        ax.set_ylim(ys[0], ys[-1])
    ax.set_xlim(xs[0],xs[-1]) if not task["xlim"] else ax.set_xlim(*task["xlim"])
    ax.set_ylabel(task["zlabel"],fontsize=12)
    ax.minorticks_on()
    iplot += 1

    # ---------------------------------------------------------------

    ax = axes[iplot] #axes[1,1]
    # ax_ = ax.twinx()
    # ax_.plot(xs_an,ej_an.get_dyn_arr(v_n="B_rs",ilayer=0,ishell=0),color='gray',ls='-',label="B")
    # # ax_.plot(xs_an,ej_an.get_dyn_arr(v_n="GammaRsh",ilayer=0,ishell=0),color='gray',ls='--',label="GammaRsh")
    # ax_.set_xscale("log"); ax_.set_yscale('log'); ax_.legend(), ax_.set_ylim(1e-2,1e3)

    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.97, 0.07, name1, fontsize=12, bbox=props,
            transform=ax.transAxes, horizontalalignment='right')

    norm = LogNorm(vmin=spec.max() * 1e-6, vmax=spec.max() * 10)
    _c = ax.pcolormesh(xs, ys, spec.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(task["ylabel"],fontsize=12)
    ax.minorticks_on()
    ax.set_xlim(xs[0],xs[-1]) if not task["xlim"] else ax.set_xlim(*task["xlim"])
    ax.set_rasterized(True)

    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
    Gamma = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs=="fs" else "Gamma",ishell=0,ilayer=0)
    z = float(ej.get_lc_obj(spec=is_spec).attrs["z"])

    if ele_syn_ssc == "n_ele":
        ax.plot(xs,gamma_min, color='red', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs,gamma_c, color='red', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs,gamma_max, color='red', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys[0],ys[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min * gamma_min * gamToNuFactor
        nu_c = XpS * gamma_c * gamma_c * gamToNuFactor
        nu_max = XpS * gamma_max * gamma_max * gamToNuFactor

        # observer frame
        if is_spec:
            ax.plot(xs,nu_min, color='red', linewidth=1, linestyle=":", label=r"$\nu'_{\rm m}$" if is_spec else r"$\nu_{\rm m}$")
            ax.plot(xs,nu_c, color='red', linewidth=1, linestyle="--", label=r"$\nu'_{\rm c}$"if is_spec else r"$\nu_{\rm c}$")
            ax.plot(xs,nu_max, color='red', linewidth=2, linestyle="-.", label=r"$\nu'_{\rm M}$"if is_spec else r"$\nu_{\rm M}$")
            ax.set_ylim(ys[0],ys[-1])
        else:
            pass
            # nu_min *= (1 + z) / Gamma
            # nu_c *= (1 + z) / Gamma
            # nu_max *= (1 + z) / Gamma

    ax.legend(fancybox=True, loc='lower left',columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=True, ncol= 3,
              fontsize=12, labelcolor='white',
              framealpha=0., borderaxespad=0.)

    iplot+=1

    # ---------------------------------------------------------------



    ax = axes[iplot]#axes[2,1]
    for idx, color in zip(indexes_ys, task["colors_ys"]):
        ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.97, 0.07, name2, fontsize=12, bbox=props,
            transform=ax.transAxes, horizontalalignment='right')

    _c = ax.pcolormesh(xs_an, ys_an, spec_an.T, cmap='jet', norm=norm)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(task["ylabel"],fontsize=12)

    ax.set_xlim(xs[0],xs[-1]) if not task["xlim"] else ax.set_xlim(*task["xlim"])
    ax.minorticks_on()
    ax.set_rasterized(True)

    gamma_min_an = ej_an.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
    gamma_c_an = ej_an.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
    gamma_max_an = ej_an.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
    B_an = ej_an.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
    Gamma_an = ej_an.get_dyn_arr(v_n="Gamma" if fs_or_rs=="fs" else "Gamma",ishell=0,ilayer=0)

    if ele_syn_ssc == "n_ele":
        ax.plot(xs_an,gamma_min_an, color='red', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        ax.plot(xs_an,gamma_c_an, color='red', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        ax.plot(xs_an,gamma_max_an, color='red', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")
        ax.set_ylim(ys_an[0],ys_an[-1])
    else:
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B_an) / (me * c)
        p = ej.pars[f"p_{fs_or_rs}"]
        XpS = 0.06 + 0.28 * p
        nu_min = XpS * gamma_min_an * gamma_min_an * gamToNuFactor
        nu_c = XpS * gamma_c_an * gamma_c_an * gamToNuFactor
        nu_max = XpS * gamma_max_an * gamma_max_an * gamToNuFactor

        # observer frame
        if is_spec:
            ax.plot(xs,nu_min, color='red', linewidth=1, linestyle=":", label=r"$\nu'_{\rm m}$" if is_spec else r"$\nu_{\rm m}$")
            ax.plot(xs,nu_c, color='red', linewidth=1, linestyle="--", label=r"$\nu'_{\rm c}$" if is_spec else r"$\nu_{\rm c}$")
            ax.plot(xs,nu_max, color='red', linewidth=2, linestyle="-.", label=r"$\nu'_{\rm M}$" if is_spec else r"$\nu_{\rm M}$")
            ax.set_ylim(ys_an[0],ys_an[-1])
        else:
            pass
            # nu_min *= (1 + z) / Gamma
            # nu_c *= (1 + z) / Gamma
            # nu_max *= (1 + z) / Gamma


    cbar = fig.colorbar(_c, ax=[axes[1],axes[2]], shrink=0.95,label=task["zlabel"],pad=0.05)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(task["zlabel"],size=12)

    iplot += 1

    # ---------------------------------------------------------------
    if task['plot_ratio']:
        ax = axes[iplot]
        for idx, color in zip(indexes_ys, task["colors_ys"]):
            ax.plot([xs[0], xs[15]], [ys[idx],ys[idx]],color=color, linewidth=2, linestyle="-")
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.97, 0.07, f"{name2}/{name1}", fontsize=12, bbox=props,
                transform=ax.transAxes, horizontalalignment='right')
        ratio = spec_an.T / spec.T
        # norm = LogNorm(vmin=ratio.min(), vmax=ratio.max())
        norm=SymLogNorm(linthresh=1e-2, vmin=1e-2, vmax=1.e2, base=10) # linthresh=0.03, linscale=0.03,
        # norm=CenteredNorm(halfrange=1.e2, vcenter=1.) # linthresh=0.03, linscale=0.03,
        # norm=TwoSlopeNorm(vmin=1e-3, vmax=1.e3, vcenter=1.) # linthresh=0.03, linscale=0.03,
        _c = ax.pcolormesh(xs_an, ys_an, ratio, cmap='RdBu_r', norm=norm)
        cbar = fig.colorbar(_c, ax=ax, shrink=0.95,label=f"{name2}/{name1}",pad=0.05)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(label=f"{name2}/{name1}",size=12)#,weight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlabel(task["xlabel"])
        ax.set_ylabel(task["ylabel"],fontsize=12)
        ax.set_xlim(xs_an[0],xs_an[-1]) if not task["xlim"] else ax.set_xlim(*task["xlim"])
        ax.set_ylim(ys_an[0],ys_an[-1])
        ax.set_rasterized(True)
        # fig.colorbar(_c, ax=ax,shrink=0.9,pad=.01,label=task["zlabel"])
        # plt.subplots_adjust(wspace=-0.2, hspace=-0.2)
        # plt.tight_layout()
        iplot+=1
    else:
        ax.set_xlabel(task["xlabel"],fontsize=12)
    # ---------------------------------------------------------------

    if task['plot_total_ratio']:
        ax = axes[iplot]
        n_ele = np.trapz(y=spec, x=ys, axis=1)#[:,np.newaxis]
        n_ele_an = np.trapz(y=spec_an, x=ys_an, axis=1)#[:,np.newaxis]
        ax.plot(xs, n_ele/n_ele_an, color="black",ls='-',lw=1.5, label=f"integrated {name2}/{name1}")
        # ax.plot(xs, n_ele_an, color="black",ls='-',lw=1.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(task["xlabel"],fontsize=12)
        # ax.set_ylabel(label="Analytic/Numeric")
        ax.set_xlim(xs_an[0],xs_an[-1]) if not task["xlim"] else ax.set_xlim(*task["xlim"])
        ax.set_ylim(1e-2,1e2)
        ax.grid(linestyle=':')

        ax.legend(fancybox=True, loc='best',columnspacing=1.2,
                       # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                       shadow=False, ncol= 1,
                       fontsize=12,
                       framealpha=0., borderaxespad=0.)
    else:
        ax.set_xlabel(task["xlabel"],fontsize=12)
    # -------------------------------------------------------------

    for ax in axes:
        ax.tick_params(axis='both',direction='in',labelsize=12)

    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.png',dpi=256)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.pdf')
    if task["show"]: plt.show()

def plot_for_fs(dyn_fs__rad_fs__mix,dyn_fs__rad_fs__num,dyn_fs__rad_fs__num__noadi,dyn_fsrs__rad_fs__num__ssc):

    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fs__rad_fs__num__noadi.GRB, name1="Adi", name2="noAdi", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\gamma_{e}$",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$(\gamma_{e} dN_{e}/d\gamma_{e}) / N_{e;\, \rm tot}$",
              "figname":"abstract_ele_fs__fs_fs__num__adi_vs_noadi","show":True,
              "title":"FS; FS dynamics; adi VS noAdi",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fs__rad_fs__num__noadi.GRB, name1="Adi", name2="noAdi", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu' j_{\rm nu'}'$ [cgs]",
              "figname":"abstract_synch_fs__fs_fs__num__adi_vs_noadi","show":True,
              "title":"FS; FS dynamics; adi VS noAdi",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )

    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fs__rad_fs__mix.GRB, name1="Numeric", name2="Mixed", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\gamma_{e}$",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$(\gamma_{e} dN_{e}/d\gamma_{e}) / N_{e;\, \rm tot}$",
              "figname":"abstract_ele_fs__fs_fs__num_vs_mix","show":True,
              "title":"FS; FS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fs__rad_fs__mix.GRB, name1="Numeric", name2="Mixed", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
        task={"ys":(1.e9,1e12,1e19), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu' j_{\rm nu'}'$ [cgs]",
              "figname":"abstract_synch_fs__fs_fs__num_vs_mix","show":True,
              "title":"FS; FS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )

    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fsrs__rad_fs__num__ssc.GRB, name1="noSSC", name2="SSC", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y",
        task={"ys":(1.e7,5e7,8e7), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\gamma_{e}$",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$(\gamma_{e} dN_{e}/d\gamma_{e}) / N_{e;\, \rm tot}$",
              "figname":"abstract_ele_fs__fs_fs__noSSC_vs_SSC","show":True,
              "title":"FS; FS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_fs__num.GRB, ej_an=dyn_fsrs__rad_fs__num__ssc.GRB, name1="noSSC", name2="SSC", is_spec=True,
        fs_or_rs = "fs", ele_syn_ssc = "total", xkey="times_freqs", ykey="freqs", norm_method="*y",
        task={"ys":(1e12,1e19,1e23), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu' j_{\rm nu'}'$ [cgs]",
              "figname":"abstract_ssc_fs__fs_fs__SSC_vs_SSC","show":True,
              "title":"FS; FS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
def plot_for_rs(dyn_fs__rad_rs__mix,dyn_fs__rad_rs__num,dyn_fs__rad_rs__num__noadi,dyn_fsrs__rad_rs__num__ssc):

    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_rs__num.GRB, ej_an=dyn_fs__rad_rs__num__noadi.GRB, name1="Adi", name2="noAdi", is_spec=True,
        fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\gamma_{e}$",
              "xlim":(3e3,1e7), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$(\gamma_{e} dN_{e}/d\gamma_{e}) / N_{e;\, \rm tot}$",
              "figname":"abstract_ele_fsrs__rs_rs__num__adi_vs_noadi","show":True,
              "title":r"RS; FS \& RS dynamics; adi VS noAdi",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_rs__num.GRB, ej_an=dyn_fs__rad_rs__num__noadi.GRB, name1="Adi", name2="noAdi", is_spec=True,
        fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(3e3,1e7), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu' j_{\rm nu'}'$ [cgs]",
              "figname":"abstract_synch_fsrs__rs_rs__num__adi_vs_noadi","show":True,
              "title":r"RS; FS \& RS dynamics; adi VS noAdi",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )

    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_rs__num.GRB, ej_an=dyn_fs__rad_rs__mix.GRB, name1="Numeric", name2="Mixed", is_spec=True,
        fs_or_rs = "rs", ele_syn_ssc = "n_ele", xkey="times_gams", ykey="gams", norm_method="/integ *y",
        task={"ys":(1.e1,1e3,1e6), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\gamma_{e}$",
              "xlim":(3e3,1e7), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$(\gamma_{e} dN_{e}/d\gamma_{e}) / N_{e;\, \rm tot}$",
              "figname":"abstract_ele_fsrs__rs_rs__num_vs_mix","show":True,
              "title":r"RS; FS \& RS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )
    plot_compare_two_spectra_new(
        ej=dyn_fs__rad_rs__num.GRB, ej_an=dyn_fs__rad_rs__mix.GRB, name1="Numeric", name2="Mixed", is_spec=True,
        fs_or_rs = "rs", ele_syn_ssc = "synch", xkey="times_freqs", ykey="freqs", norm_method="*y",
        task={"ys":(1.e9,1e12,1e18), "colors_ys":("cyan","green","orange"), "ylim":(), "ylabel":r"$\nu'$ [Hz]",
              "xlim":(3e3,1e7), "xlabel":r"$t_{\rm burst}$ [s]",
              "zlabel":r"$\nu' j_{\rm nu'}'$ [cgs]",
              "figname":"abstract_synch_fsrs__rs_rs__num_vs_mix","show":True,
              "title":r"RS; FS \& RS dynamics; numeric VS mixed",
              "subplots":dict(
                  ncols=1, nrows=5-1, figsize=(6,7),
                  sharex="col", #sharey="row",sharex="col",
                  gridspec_kw=dict(height_ratios=[1.0,1.,1., 1.]),# 1.,0.4]),#dict(height_ratios=[0.5,1.2]),
                  layout='constrained'),
              "plot_ratio":True, "plot_total_ratio":False
              }
    )

if __name__ == '__main__':
    do_run = True
    struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    P = dict(
        main=dict(n_ism=1., tb0=3e3, ntb=1000,rtol=1e-7,theta_obs=0,
                  lc_freqs='array logspace 1e8 1e27 64',
                  lc_times='array logspace 3e3 1e10 128'),
        grb=dict(save_dynamics='yes',do_spec='yes',save_spec='yes',do_lc='yes',
                 method_nonrel_dist_fs='none',
                 method_nonrel_dist_rs='none',
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 freq1=1e6,freq2=1e30,nfreq=400)
    )

    # --- fs -- fs ---
    dyn_fs__rad_fs__num = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num/",
        struct=struct, P = P, type="a",run=do_run
    )
    dyn_fs__rad_fs__num__noadi = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num__noadi/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(num_ele_use_adi_loss_fs='no'))),
        type="a",run=do_run
    )
    dyn_fs__rad_fs__num__ssc = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num__ssc/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(method_ssc_fs='numeric'))),
        type="a",run=do_run
    )
    dyn_fs__rad_fs__mix = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__mix/",
        struct=struct,
        P = d2d(default=P, new=dict(grb=dict(method_ele_fs='mix'))),
        type="a",run=do_run
    )
    plot_for_fs(dyn_fs__rad_fs__mix,dyn_fs__rad_fs__num,dyn_fs__rad_fs__num__noadi,dyn_fs__rad_fs__num__ssc)

    # --- fsrs -- rs ---
    dyn_fsrs__rad_rs__num = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__num/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs'))),
        type="a",run=do_run
    )
    dyn_fsrs__rad_rs__num__noadi = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__num__noadi/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            num_ele_use_adi_loss_fs='no',
                                                            num_ele_use_adi_loss_rs='no'))),
        type="a",run=do_run
    )
    dyn_fsrs__rad_rs__num__ssc = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__num__ssc/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            method_ssc_fs='numeric',
                                                            method_ssc_rs='numeric'))),
        type="a",run=do_run
    )
    dyn_fsrs__rad_rs__mix = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__mix/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            method_ele_fs='mix',
                                                            method_ele_rs='mix'))),
        type="a",run=do_run
    )
    plot_for_rs(dyn_fsrs__rad_rs__mix,dyn_fsrs__rad_rs__num,dyn_fsrs__rad_rs__num__noadi,dyn_fsrs__rad_rs__num__ssc)
