import copy

import package.src.PyBlastAfterglowMag as PBA
import os,shutil,matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize,TwoSlopeNorm,SymLogNorm,TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator,AutoLocator
from matplotlib import ticker
import matplotlib.ticker as plticker
import numpy as np
from scipy import special

from package.src.PyBlastAfterglowMag import Ejecta

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
    elif norm_method == "*y /integ":
        spec *= np.power(ys, 1)[np.newaxis,:]
        spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
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


class PlotSpectra:
    def __init__(self):
        pass

    @staticmethod
    def _get_spectrum(ej:PBA.Ejecta,v_n:str,fs_or_rs:str,norm_method:str,
                      ishell:int=0,ilayer:int=0,sum_shells_layers:bool=False):

        if v_n in ["n_ele","gam_dot_syn","gam_dot_adi","gam_dot_ssc","yg"]:
            xkey = "times_gams"
            ykey = "gams"
        else:
            xkey = "times_freqs"
            ykey = "freqs"

        if v_n == "fluxdense": is_spec = False
        else: is_spec = True

        ele_syn_ssc = v_n
        if v_n == "syn_tau": ele_syn_ssc = "syn_a"

        if v_n == "total_f":
            spec_syn=ej.get_lc(key="syn_j"+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                           xkey=xkey,ykey=ykey,freq=None,time=None,
                           ishell=ishell,ilayer=ilayer,
                           spec=is_spec,sum_shells_layers=sum_shells_layers)
            spec_ssc=ej.get_lc(key="ssc_j"+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                               xkey=xkey,ykey=ykey,freq=None,time=None,
                               ishell=ishell,ilayer=ilayer,
                               spec=is_spec,sum_shells_layers=sum_shells_layers)
            spec = spec_syn + spec_ssc
        elif v_n == "total_i":
            spec_syn=ej.get_lc(key="syn_i"+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                               xkey=xkey,ykey=ykey,freq=None,time=None,
                               ishell=ishell,ilayer=ilayer,
                               spec=is_spec,sum_shells_layers=sum_shells_layers)
            spec_ssc=ej.get_lc(key="ssc_i"+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                               xkey=xkey,ykey=ykey,freq=None,time=None,
                               ishell=ishell,ilayer=ilayer,
                               spec=is_spec,sum_shells_layers=sum_shells_layers)
            spec = spec_syn + spec_ssc
        else:
            spec=ej.get_lc(key=ele_syn_ssc+'_'+fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                           xkey=xkey,ykey=ykey,freq=None,time=None,
                           ishell=ishell,ilayer=ilayer,
                           spec=is_spec,sum_shells_layers=sum_shells_layers)
        # if v_n == "Yg":
        #     spec = PlotSpectra._get_spectrum(ej=ej,v_n="total_f",fs_or_rs=fs_or_rs,norm_method=norm_method,
        #                                      ishell=ishell,ilayer=ilayer,sum_shells_layers=sum_shells_layers)

        if v_n == "syn_tau":
            dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs=="fs" else "thichness_rs",ishell=ishell,ilayer=ilayer)
            Gamma = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs=="fs" else "Gamma43",ishell=ishell,ilayer=ilayer)
            dr_comov = dr * Gamma # thickness as R / Gamma in comoving frame (dr is in observer frame so R/Gamma^2
            spec = spec * dr_comov[:, np.newaxis]
        if (np.sum(spec) == 0):
            raise ValueError("(np.sum(spec) == 0)")
        xs = ej.get_grid(key=xkey,spec=is_spec)
        ys = ej.get_grid(key=ykey,spec=is_spec)

        spec = _normalize_spec(spec=spec,xs=xs,ys=ys,norm_method=norm_method,mask_negatives=True)

        return (xs, ys, spec)
    @staticmethod
    def _gam_to_nu(gam:np.ndarray,B:np.ndarray,p:float,syn_or_ssc:str):
        """ Using Johanesson+2006 for XpS and Sari & Esin for SSC """
        qe = 4.803204e-10
        me = 9.1094e-28
        c  = 2.99792458e10
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
        XpS = 0.06 + 0.28 * p
        nu = XpS * gam**2 * gamToNuFactor

        if syn_or_ssc.__contains__("ssc"): return 4.* gam ** 2 * nu * np.sqrt(2)/3.
        return nu
        # nu_min = XpS * gamma_min * gamma_min * gamToNuFactor
        # nu_c = XpS * gamma_c * gamma_c * gamToNuFactor
        # nu_max = XpS * gamma_max * gamma_max * gamToNuFactor

    @staticmethod
    def _plot_ele_spectrum(ax, fig, xs:np.ndarray, ys:np.ndarray, spec:np.ndarray,
                           ej:PBA.Ejecta, fs_or_rs:str, task_i:dict):
        # ys=PBA.utils.get_beta(ys)*ys
        # task_i = task[v_n]
        if task_i['norm'] == "LogNorm":
            norm = LogNorm(vmin=spec.max() * 1e-5, vmax=spec.max())
        elif task_i['norm'] == "SymLogNorm":
            norm=SymLogNorm(linthresh=task_i.get("vmin", 1e-4), vmin=task_i.get("vmin", 1e-4),
                            vmax=task_i.get("vmax", 1e4), base=10)
            # norm=TwoSlopeNorm(vmin=0.1,vcenter=1,vmax=10)
        else:
            raise KeyError(f"norm {task_i['norm']} is not recognized")
        cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
        cmap.set_under(task_i.get("set_under",'white'),alpha=0.2)
        cmap.set_over(task_i.get("set_over",'white'),alpha=0.2)
        spec = np.ma.masked_where(spec < spec.max() * 1e-5, spec)
        if task_i["mode"]=="contour":
            _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)
        cbar = fig.colorbar(_c, ax=ax, shrink=0.95,pad=0.01,extend='both')# orientation = 'horizontal')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label(task_i["zlabel"],size=12)
        # cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(n=6))
        # cbar.ax.yaxis.set_major_locator(AutoLocator())
        # cbar.ax.set_yticks(np.logspace(np.log10(spec.max() * 1e-6),np.log10(spec.max() * 10),5))
        # cbar.ax.set_yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e1,1e2,1e3])
        ax.set_ylabel(task_i["ylabel"],fontsize=12)
        if task_i["plot_gm"]:
            gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
            ax.plot(xs,gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        if task_i["plot_gc"]:
            gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
            ax.plot(xs,gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        if task_i["plot_gM"]:
            gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
            ax.plot(xs,gamma_max, color='black', linewidth=1, linestyle="-.", label=r"$\gamma_{\rm M}$")
        if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
        if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])
        return ax

    @staticmethod
    def _plot_rad_spectrum(ax, fig, xs:np.ndarray, ys:np.ndarray, spec:np.ndarray,
                           ej:PBA.Ejecta, fs_or_rs:str, v_n:str, task_i:dict):
        # task_i = task[v_n]
        if task_i['norm'] == "LogNorm":
            norm = LogNorm(vmin=task_i.get("vmin", spec.max() * 1e-5), vmax=task_i.get("vmax", spec.max() * 2))
        elif task_i['norm'] == "SymLogNorm":
            norm=SymLogNorm(linthresh=task_i.get("vmin", 1e-4),
                            vmin=task_i.get("vmin", 1e-4),
                            vmax=task_i.get("vmax", 1e4), base=10)
        else:
            raise KeyError(f"norm {task_i['norm']} is not recognized")

        cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
        cmap.set_under(task_i.get("set_under",'white'),alpha=0.2)
        cmap.set_over(task_i.get("set_over",'white'),alpha=0.2)

        if not v_n.__contains__("tau") and task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec > norm.vmax, spec)
            _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)

        if task_i["mode"] == "contour":
            _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)

        cbar = fig.colorbar(_c, ax=ax, shrink=0.95,pad=0.01,extend='both')# orientation = 'horizontal')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label(task_i["zlabel"],size=12)
        # tick_locator = ticker.MaxNLocator(nbins=15)
        # cbar.locator = tick_locator
        # cbar.set_norm(norm)
        # cbar.update_normal(_c)
        # cbar.update_ticks()
        ax.set_ylabel(task_i["ylabel"],fontsize=12)
        if task_i["plot_num"]:
            gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            ax.plot(xs,PlotSpectra._gam_to_nu(gam=gamma_min,B=B,p=p,syn_or_ssc=v_n),
                    color='black', linewidth=1, linestyle=":",
                    label=r"$\nu'_{\rm m}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm m;\, ssc}$")
        if task_i["plot_nuc"]:
            gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            ax.plot(xs,PlotSpectra._gam_to_nu(gam=gamma_c,B=B,p=p,syn_or_ssc=v_n),
                    color='black', linewidth=1, linestyle="--",
                    label=r"$\nu'_{\rm c}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm c;\, ssc}$")
        if task_i["plot_nuM"]:
            gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=0,ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            ax.plot(xs,PlotSpectra._gam_to_nu(gam=gamma_max,B=B,p=p,syn_or_ssc=v_n),
                    color='black', linewidth=1, linestyle="-.",
                    label=r"$\nu'_{\rm M}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm M;\, ssc}$")
        if task_i["plot_nua"]:
            Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs=="fs" else "GammaRsh",ishell=0,ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            nu_a = 3.41 * 1e9 / Gamma_sh * (p+2)/(3.*p+2)*(p-1)**(8/5)/(p-2) # From Warren:2018lyx

            theta = ej.get_dyn_arr(v_n="theta",ishell=0,ilayer=0)
            mask = np.array(Gamma_sh < 0.95 * Gamma_sh[0], dtype=int) * np.array(theta < 1.1 * theta[0], dtype=int)
            mask = np.array(mask,dtype=bool)
            if np.sum(mask) > 0:
                ax.plot(xs[mask],nu_a[mask], color='black', linewidth=1, linestyle=":", label=r"$\nu'_{\rm a}$")
        if task_i["plot_tau1"]:
            Gamma = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs=="fs" else "GammaRsh",ishell=0,ilayer=0)
            xs, ys, tau = PlotSpectra._get_spectrum(ej=ej,v_n="syn_tau",fs_or_rs=fs_or_rs,norm_method=None)
            ax.plot(
                xs,
                [ys[PBA.utils.find_nearest_index(tau[i, :], 1.)] for i in range(len(Gamma))], # [xs,ys]
                color='black', lw=0.8, ls='-',
                label=r"$\tau_{\nu'}'=1$"
            )
        if task_i["plot_max"]:
            ax.plot(
                xs,
                [ys[idx] for idx in np.argmax(spec,axis=1)], # [xs,ys]
                color='gray', lw=0.9, ls='-',
                label=r"$\nu'_{\rm p}$"
            )

        if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
        if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])

        return ax


def plot_spectra_evolution(ej:PBA.Ejecta, fs_or_rs:str, tasks:dict,title:str or None,
                           figsize:tuple, figname:str,show:bool):
    n_tasks = len(tasks.keys())
    fig, axes = plt.subplots(ncols=1,nrows=n_tasks,sharex='all',layout='constrained',figsize=figsize)
    if not hasattr(axes,'__len__'):
        axes = [axes]
    i_plot = 0

    pl = PlotSpectra()

    for v_n, task in tasks.items():
        if v_n in ["n_ele","gam_dot_syn","gam_dot_adi","gam_dot_ssc","yg"]:
            xs, ys, spec = pl._get_spectrum(ej=ej,v_n=v_n,fs_or_rs=fs_or_rs,norm_method=task["norm_method"])
            pl._plot_ele_spectrum(ax=axes[i_plot], fig=fig, xs=xs,ys=ys,spec=spec, ej=ej, fs_or_rs=fs_or_rs, task_i=task)
        else:
            xs, ys, spec = pl._get_spectrum(ej=ej,v_n=v_n,fs_or_rs=fs_or_rs,norm_method=task["norm_method"])
            pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec, fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n=v_n, task_i=task)
        i_plot += 1

    #
    #
    # # plot ele spectrum
    #
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="n_ele",fs_or_rs=fs_or_rs,norm_method=task["n_ele"]["norm_method"])
    # pl._plot_ele_spectrum(ax=axes[i_plot], fig=fig, xs=xs,ys=ys,spec=spec, ej=ej, fs_or_rs=fs_or_rs, task_i=task['n_ele'])
    # i_plot += 1
    #
    # # plot synch spectrum
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="syn_j",fs_or_rs=fs_or_rs,norm_method=task["syn_j"]["norm_method"])
    # pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec, fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n="syn_j", task=task)
    # i_plot += 1
    #
    # # plot ssc spectrum
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="ssc_j",fs_or_rs=fs_or_rs,norm_method=task["ssc_j"]["norm_method"])
    # pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec,fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n="ssc_j", task=task)
    # i_plot += 1
    #
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="syn_f",fs_or_rs=fs_or_rs,norm_method=task["syn_f"]["norm_method"])
    # pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec,fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n="syn_f", task=task)
    # i_plot += 1
    #
    # # plot ssa spectrum
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="syn_a",fs_or_rs=fs_or_rs,norm_method=task["syn_a"]["norm_method"])
    # pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec,fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n="syn_a", task=task)
    # i_plot += 1
    #
    # # plot ssa spectrum
    # xs, ys, spec = pl._get_spectrum(ej=ej,v_n="syn_tau",fs_or_rs=fs_or_rs,norm_method=task["syn_tau"]["norm_method"])
    # pl._plot_rad_spectrum(ax=axes[i_plot], xs=xs,ys=ys,spec=spec,fig=fig, ej=ej, fs_or_rs=fs_or_rs, v_n="syn_tau", task=task)
    # i_plot += 1

    for ax in axes:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='both',direction='in',labelsize=11)
        ax.legend(fancybox=True, loc='upper right',columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol= 4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        # start, end = ax.get_xlim()
        # ax.xaxis.set_ticks(np.logspace(np.log10(start), np.log10(end), 9))
        # loc = plticker.MultipleLocator(base=100.0) # this locator puts ticks at regular intervals
        # ax.xaxis.set_major_locator(loc)
        # ax.locator_params(axis='y',nbins=10)
    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    axes[0].set_title(title, fontsize=12)
    plt.savefig(os.getcwd()+'/figs/'+figname+'.png',dpi=256)
    plt.savefig(os.getcwd()+'/figs/'+figname+'.pdf')
    if show: plt.show()
    if plot: plt.show()
    plt.close(fig)
def plot_spectra_evolution_ratio(ej1:PBA.Ejecta, ej2:PBA.Ejecta, v_n:str, fs_or_rs:str, task:dict):

    pl = PlotSpectra()

    fig, ax = plt.subplots(ncols=1,nrows=1,sharex='all',layout='constrained',figsize=(5.4,3.0))
    task_i = task[v_n]
    xs, ys, spec1 = pl._get_spectrum(ej=ej1,v_n=v_n,fs_or_rs=fs_or_rs,norm_method=task_i["norm_method"])
    xs, ys, spec2 = pl._get_spectrum(ej=ej2,v_n=v_n,fs_or_rs=fs_or_rs,norm_method=task_i["norm_method"])

    spec = spec1 / spec2
    if v_n == "n_ele":
        pl._plot_ele_spectrum(ax=ax, fig=fig, xs=xs,ys=ys,spec=spec, ej=ej1, fs_or_rs=fs_or_rs, task_i=task[v_n])
    else:
        pl._plot_rad_spectrum(ax=ax, fig=fig, xs=xs,ys=ys,spec=spec, ej=ej1, fs_or_rs=fs_or_rs, v_n=v_n, task_i=task[v_n])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both',which='both',direction='in',labelsize=11)
    ax.legend(fancybox=True, loc='upper right',columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 4, fontsize=12, labelcolor='black',
              framealpha=0.4, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.set_title(task["title"], fontsize=12)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.png',dpi=256)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)

def plot_emission_region_prop(ej:PBA.Ejecta,fs_or_rs:str,xlim:tuple,task:dict,ishell=0,ilayer=0):
    fig,axes = plt.subplots(ncols=1,nrows=2, figsize=(5,4.5),layout='constrained',sharex="all")
    #     double T = source.dr / CGS::c; // escape time (See Huang+2022)
    # double fac = T / source.vol;
    x = ej.get_dyn_arr(v_n="tburst",ishell=ishell,ilayer=ilayer)
    r = ej.get_dyn_arr(v_n="R",ishell=ishell,ilayer=ilayer)
    Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs=="fs" else "GammaRsh",ishell=ishell,ilayer=ilayer)
    dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs=="fs" else "thichness_rs",ishell=ishell,ilayer=ilayer)
    mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs=="fs" else "M3",ishell=ishell,ilayer=ilayer)
    dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs=="fs" else "rho3",ishell=ishell,ilayer=ilayer)
    vol = mass / dens
    tau_tomps = PBA.utils.cgs.sigmaT * (mass/PBA.utils.cgs.mp / vol) * dr * Gamma_sh
    tmp = dr*Gamma_sh*PBA.utils.cgs.c / vol
    # ax.plot(
    #     x, PBA.utils.cgs.c/r**2
    # )


    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=ishell,ilayer=ilayer)
    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=ishell,ilayer=ilayer)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs=="fs" else "B_rs",ishell=ishell,ilayer=ilayer)
    delta_t_syn = PBA.utils.cgs.sigmaT * gamma_max * gamma_max * B * B / (6. * np.pi * PBA.utils.cgs.me * PBA.utils.cgs.c)
    tcomov = ej.get_dyn_arr(v_n="tcomov",ishell=ishell,ilayer=ilayer)
    delta_t_adi = ((gamma_max*gamma_max-1.)/(3.*gamma_max))[:-1] * \
                  ((tcomov[1:]-tcomov[:-1])**-1 * (1. - vol[:-1] / vol[1:]))




    sp = PlotSpectra()


    xs, ys, spec = sp._get_spectrum(ej=ej,v_n="n_ele",fs_or_rs=fs_or_rs,
                                    ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    # Energy density in relativistic particles https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf
    tot_ele = np.trapz(spec*ys[np.newaxis,:]/vol[:,np.newaxis]*PBA.utils.cgs.me*PBA.utils.cgs.c**2,x=ys,axis=1)

    xs, ys, spec = sp._get_spectrum(ej=ej,v_n="syn_f",fs_or_rs=fs_or_rs,
                                    ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    tot_syn = np.trapz(PBA.utils.cgs.h*spec*ys[np.newaxis,:],x=ys,axis=1) # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    sp = PlotSpectra()
    xs, ys, spec = sp._get_spectrum(ej=ej,v_n="ssc_f",fs_or_rs=fs_or_rs,
                                    ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    tot_ssc = np.trapz(PBA.utils.cgs.h*spec*ys[np.newaxis,:],x=ys,axis=1) # PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    ax = axes[0]
    ax.plot(xs, tot_ele,color='red',label=r"$u'_{\rm ele}$")
    ax.set_ylabel(r"Energy density of electrons",color='red',fontsize=12)
    # ax.legend()

    ax1 = ax.twinx()
    # ax1.plot(xs, tot_ssc/tot_syn,color='blue',label=r'$u_{\rm syn}$')
    ax1.plot(xs, tot_syn,color='blue',label=r"$u'_{\rm syn}$")
    ax1.plot(xs, tot_ssc,color='green',label=r"$u'_{\rm ssc}$")
    # ax1.plot(xs, tot_ele,color='red',label=r'$u_{\rm ele}$')
    ax1.set_ylabel(r"Energy density of radiation",color='black',fontsize=12)
    # ax1.legend()
    # ax1.plot(x, gamma_min)


    ax = axes[1]
    ax.plot(x, tmp, color='black', label=r"$\Delta t' / V'$")
    ax.set_ylabel(r"Escape time / V' ",color='black',fontsize=12)
    # ax.legend()

    ax2 = ax.twinx()
    ax2.plot(xs, tau_tomps,color='magenta',label=r'$\tau_{\rm comp}$')
    ax2.set_ylabel(r"Compton optical depth",color='magenta',fontsize=12)
    # ax2.legend()


    # ax = axes[2]
    # # ax.plot(x, delta_t_syn, color='blue', label=r"$\dot{\gamma}'_{\rm syn}$")
    # # ax.plot(x[:-1], delta_t_adi, color='red', label=r"$\dot{\gamma}'_{\rm adi}$")
    #
    # delta_t_syn = PBA.utils.cgs.sigmaT * gamma_min * gamma_min * B * B / (6. * np.pi * PBA.utils.cgs.me * PBA.utils.cgs.c)
    # delta_t_adi = ((gamma_min*gamma_min-1.)/(3.*gamma_min))[:-1] * \
    #               ((tcomov[1:]-tcomov[:-1])**-1 * (1. - vol[:-1] / vol[1:]))
    # ax.plot(x, delta_t_syn, color='blue', ls='--', label=r"$\dot{\gamma}'_{\rm syn}$")
    # ax.plot(x[:-1], delta_t_adi, color='red',  ls='--',label=r"$\dot{\gamma}'_{\rm adi}$")

    for ax,loc in zip([axes[0],ax1, axes[1], ax2],
                      ["lower left","upper right","center left","center right"]):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='both',direction='in',labelsize=11)
        ax.legend(fancybox=True, loc=loc,columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol= 4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)

    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]",fontsize=12)
    # ax.set_ylabel(r"$\tau_{\rm e}$")
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.png',dpi=256)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)

def plot_cooling_terms(ej:PBA.Ejecta,fs_or_rs:str,xlim:tuple,task:dict,ishell=0,ilayer=0):
    task_i=task
    sp = PlotSpectra()
    xs, ys, gam_dot_syn = sp._get_spectrum(
        ej=ej,v_n="gam_dot_syn",fs_or_rs=fs_or_rs,
        ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    xs, ys, gam_dot_adi = sp._get_spectrum(
        ej=ej,v_n="gam_dot_adi",fs_or_rs=fs_or_rs,
        ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    xs, ys, gam_dot_ssc = sp._get_spectrum(
        ej=ej,v_n="gam_dot_ssc",fs_or_rs=fs_or_rs,
        ishell=ishell,ilayer=ilayer,sum_shells_layers=False,norm_method=None)
    gam_tot = gam_dot_syn + gam_dot_adi + gam_dot_ssc

    fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(5,3),layout='constrained')
    norm = LogNorm(vmin=gam_tot.max()*1e-20, vmax=gam_tot.max() * 1)
    cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
    cmap.set_under(task_i.get("set_under",'white'),alpha=0.2)
    cmap.set_over(task_i.get("set_over",'white'),alpha=0.2)
    _c = ax.pcolormesh(xs, ys, gam_tot.T, cmap=cmap, norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95,pad=0.01,extend='both')# orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"],size=12)

    # use contourf() with proper hatch pattern and alpha value
    val1 = np.array((gam_dot_syn>gam_dot_adi) * (gam_dot_syn>gam_dot_ssc), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(xs, ys, val1.T  , hatches=['//'],   alpha=0.0,color = "green") #,zorder=0

    val1 = np.array((gam_dot_adi>gam_dot_syn) * (gam_dot_adi>gam_dot_ssc), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(xs, ys, val1.T  , hatches=["\\\\"],  alpha=0.0) #,zorder=0

    val1 = np.array((gam_dot_ssc>gam_dot_syn) * (gam_dot_ssc>gam_dot_adi), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(xs, ys, val1.T  , hatches=['++'],  alpha=0.0, edgecolor = "r",facecolor="blue", color = "green") #,zorder=0


    if task_i["plot_gm"]:
        gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
        ax.plot(xs,gamma_min, color='blue', linewidth=2, linestyle=":", label=r"$\gamma_{\rm m}$")
    if task_i["plot_gc"]:
        gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
        ax.plot(xs,gamma_c, color='blue', linewidth=2, linestyle="--", label=r"$\gamma_{\rm c}$")
    if task_i["plot_gM"]:
        gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
        ax.plot(xs,gamma_max, color='blue', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both',which='both',direction='in',labelsize=11)
    ax.legend(fancybox=True, loc='upper right',columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 4, fontsize=12, labelcolor='black',
              framealpha=0.8, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_ylabel(r"$\gamma_{e}$", fontsize=12)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.set_title(task["title"], fontsize=12)
    ax.set_xlim(*xlim)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.png',dpi=256)
    plt.savefig(os.getcwd()+'/figs/'+task["figname"]+'.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def tasks_fs(do_run:bool, plot:bool, struct:dict, P:dict):
    dyn_fs__rad_fs__num__ssa__ssc = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num__ssa__ssc/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(method_ssc_fs='numeric',use_ssa_fs='yes'))),
        type="a",run=do_run
    )

    # plot_cooling_terms(dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs="fs",xlim=(3e3,1e9),
    #                    task=dict(show=True,figname="ele_cool_terms",zlabel=r"$\dot{\gamma}$",title="Cooling terms",
    #                              set_under='blue',set_over='red',
    #                              plot_gm=True,plot_gc=True,plot_gM=True))
    # #
    # plot_emission_region_prop(dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs="fs",xlim=(3e3,1e9),
    #                           task=dict(show=True,figname="rad_ele_energy_dens_fs"))

    plot_spectra_evolution(
        ej=dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs='fs',
        title="FS comoving spectra evolution", figname="spec_dyn_fs__rad_fs__num__ssa__ssc", show=plot,
        figsize=(5,4),
        tasks=dict(
            # n_ele=dict(norm_method='/integ *y',ylabel=r"$\gamma_{e}$",
            #            zlabel=r"$(\gamma_{e} N_{e}) / N_{e;\, \rm tot}$",
            #            norm="LogNorm",mode="contour", xlim=(3e3,1e9),
            #            plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False),
            # yg=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$Y_g$ [cgs]",
            #         norm="LogNorm",
            #         plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False),
            # syn_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm syn}$ [cgs]",
            #            norm="LogNorm",ylim=(1e12,1e24), xlim=(3e3,1e9),vmin=1e-9,vmax=1e-3,mode="contour",
            #            plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True),
            # ssc_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm ssc}$ [cgs]",
            #          norm="LogNorm", ylim=(1e18,1e30), xlim=(3e3,1e9),vmin=1e-9,vmax=1e-3,mode="contour",
            #          plot_num=True,plot_nuc=True,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=True),
            # syn_a=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\alpha'_{\rm syn}$ [cgs]",
            #          norm="LogNorm",mode="contour",
            #          plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
            #          ylim=(1e6,1e12), xlim=(3e3,1e9)),
            syn_tau=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\tau'_{\rm ssa}$",
                     norm="SymLogNorm",mode=None,
                     plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=True,plot_tau1=True,plot_max=False,
                     ylim=(1e6,1e14), xlim=(3e3,1e9)),
            total_i=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' I'_{\rm total}$ [cgs]",
                     norm="LogNorm",mode="contour",
                     plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
                     ylim=(1e7,1e28), xlim=(3e3,1e9), vmin=1e-15,vmax=1e-3),# vmin=1e-17,vmax=1e-3
        ))


def tasks_rs(do_run:bool, plot:bool, struct:dict, P:dict):
    dyn_fs__rad_rs__num__ssa__ssc = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_rs__num__ssa__ssc/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            use_ssa_fs='yes',
                                                            use_ssa_rs='yes',
                                                            method_ssc_fs='numeric',
                                                            method_ssc_rs='numeric'))),
        type="a",run=do_run
    )

    plot_cooling_terms(dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs="rs",xlim=(3e3,1e9),
                       task=dict(show=True,figname="ele_cool_terms_rs",zlabel=r"$\dot{\gamma}$",title="Cooling terms",
                                 set_under='blue',set_over='red',
                                 plot_gm=True,plot_gc=True,plot_gM=True))
    #
    plot_emission_region_prop(dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs="rs",xlim=(3e3,1e9),
                              task=dict(show=True,figname="rad_ele_energy_dens_rs"))

    # plot_spectra_evolution(
    #     ej=dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs='rs',
    #     title="RS comoving spectra evolution", figname="spec_dyn_fs__rad_rs__num__ssa__ssc", show=plot,
    #     tasks=dict(
    #         n_ele=dict(norm_method='/integ *y',ylabel=r"$\gamma_{e}$",
    #                    zlabel=r"$(\gamma_{e} N_{e}) / N_{e;\, \rm tot}$",
    #                    norm="LogNorm",
    #                    plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False,
    #                    ylim=(1,1e5)),
    #         synch=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm syn}$ [cgs]",
    #                    norm="LogNorm",
    #                    plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True),
    #         ssc=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm ssc}$ [cgs]",
    #                  norm="LogNorm",
    #                  plot_num=True,plot_nuc=True,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=True),
    #         ssa=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\alpha'_{\rm syn}$ [cgs]",
    #                  norm="LogNorm",
    #                  plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
    #                  ylim=(1e6,1e14)),
    #         tau=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\tau'_{\rm ssa}$",
    #                  norm="SymLogNorm",
    #                  plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=True,plot_tau1=True,plot_max=False,
    #                  ylim=(1e6,1e14), xlim=(3e3,1e7))
    #     ))
    plot_spectra_evolution(
        ej=dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs='rs',
        title="FS comoving spectra evolution", figname="spec_dyn_fsrs__rad_rs__num__ssa__ssc", show=plot,
        tasks=dict(
            n_ele=dict(norm_method='/integ *y',ylabel=r"$\gamma_{e}$",
                       zlabel=r"$(\gamma_{e} N_{e}) / N_{e;\, \rm tot}$",
                       norm="LogNorm",mode=None,#"contour",
                       plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False,
                       ylim=(1,1e4)),
            # yg=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$Y_g$ [cgs]",
            #         norm="LogNorm",
            #         plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False),
            syn_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm syn}$ [cgs]",
                       norm="LogNorm",ylim=(1e7,1e24),mode=None,#,vmin=1e-9,vmax=1e-3 "contour",
                       plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True),
            ssc_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm ssc}$ [cgs]",
                       norm="LogNorm", ylim=(1e16,1e28),mode=None,#,vmin=1e-9,vmax=1e-3 "contour",
                       plot_num=True,plot_nuc=True,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=True),
            # syn_a=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\alpha'_{\rm syn}$ [cgs]",
            #          norm="LogNorm",mode="contour",
            #          plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
            #          ylim=(1e6,1e12)),
            syn_tau=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\tau'_{\rm ssa}$",
                         norm="SymLogNorm",mode=None,
                         plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=True,plot_tau1=True,plot_max=False,
                         ylim=(1e7,1e14), xlim=(3e3,1e7)),
            total_i=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' I'_{\rm total}$ [cgs]",
                         norm="LogNorm",mode=None,#"contour",
                         plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
                         ylim=(1e7,1e28), xlim=(3e3,1e7)),# ,vmin=1e-15,vmax=1e-3 vmin=1e-17,vmax=1e-3
        ))


def tasks_fs_comparison(do_run:bool, plot:bool, struct:dict, P:dict):


    # --- fs -- fs ---
    dyn_fs__rad_fs__num = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num/",
        struct=struct, P = P,
        type="a",run=do_run
    )
    dyn_fs__rad_fs__ana__ssa = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__ana_ssa/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(method_ele_fs='analytic',use_ssa_fs='yes',method_synchrotron_fs="Dermer09"))),#dyn_fs__rad_rs__num__ssa
        type="a",run=do_run
    )
    dyn_fs__rad_fs__num__ssa = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__num_ssa/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(use_ssa_fs='yes'))),#method_synchrotron_fs="Dermer09"
        type="a",run=do_run
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

    for (model1, model2, name, label) in [
        (dyn_fs__rad_fs__num.GRB, dyn_fs__rad_fs__mix.GRB,"num_mix","Numeric/Mixed"),
        (dyn_fs__rad_fs__num.GRB, dyn_fs__rad_fs__num__noadi.GRB,"adi_noadi","Adi/noAdi"),
        (dyn_fs__rad_fs__num__ssa.GRB, dyn_fs__rad_fs__num.GRB,"ssa_nossa","SSA/noSSA"),
        (dyn_fs__rad_fs__num__ssc.GRB, dyn_fs__rad_fs__num.GRB,"ssc_nossc","SSC/noSSC"),
    ]:
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="n_ele", fs_or_rs="fs", task=dict(
                title="FS comoving electron spectra evolution ratio",
                figname=f"spec_fs_ele_ratio_{name}", show=plot,
                n_ele=dict(norm_method=None,ylabel=r"$\gamma_{e}$",
                           zlabel=label, vmin=1e-2,vmax=1e2, cmap="RdBu_r",xlim=(3e3,8e6),set_under='blue',set_over='red',mode=None,
                           norm="SymLogNorm",plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False)))
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="synch", fs_or_rs="fs", task=dict(
                title="FS comoving synchrotron spectra evolution ratio",
                figname=f"spec_fs_synch_ratio_{name}", show=plot,
                synch=dict(norm_method=None,ylabel=r"$\nu'$",
                           zlabel=label, vmin=1e-2,vmax=1e2, cmap="RdBu_r", xlim=(3e3,8e6),
                           norm="SymLogNorm",set_under='blue',set_over='red',mode=None,
                           plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True)))
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="int", fs_or_rs="fs", task=dict(
                title="FS comoving intensity spectra evolution ratio",
                figname=f"spec_fs_int_ratio_{name}", show=plot,
                int=dict(norm_method=None,ylabel=r"$\nu'$",
                         zlabel=label, vmin=1e-2,vmax=1e2, cmap="RdBu_r", xlim=(3e3,8e6), #ylim=(1e6,1e14),
                         norm="SymLogNorm",set_under='blue',set_over='red',mode=None,
                         plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=False)))

def tasks_rs_comparison(do_run:bool, plot:bool, struct:dict, P:dict):

    # --- fsrs -- rs ---
    dyn_fsrs__rad_rs__num = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__num/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs'))),
        type="a",run=do_run
    )
    dyn_fsrs__rad_rs__ana__ssa = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__ana__ssa/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            method_synchrotron_fs="Dermer09",
                                                            method_synchrotron_rs="Dermer09",
                                                            use_ssa_fs='yes',
                                                            use_ssa_rs='yes'))),
        type="a",run=do_run
    )
    dyn_fsrs__rad_rs__num__ssa = run(
        working_dir=os.getcwd()+"/working_dirs/dyn_fsrs__rad_rs__num__ssa/",
        struct=struct, P = d2d(default=P, new=dict(grb=dict(do_rs='yes',bw_type='fsrs',
                                                            use_ssa_fs='yes',
                                                            use_ssa_rs='yes'))),
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

    for (model1, model2, name, label) in [
        (dyn_fsrs__rad_rs__num.GRB, dyn_fsrs__rad_rs__mix.GRB,"num_mix","Numeric/Mixed"),
        (dyn_fsrs__rad_rs__num.GRB, dyn_fsrs__rad_rs__num__noadi.GRB,"adi_noadi","Adi/noAdi"),
        (dyn_fsrs__rad_rs__num__ssa.GRB, dyn_fsrs__rad_rs__num.GRB,"ssa_nossa","SSA/noSSA"),
        (dyn_fsrs__rad_rs__num__ssc.GRB, dyn_fsrs__rad_rs__num.GRB,"ssc_nossc","SSC/noSSC"),
    ]:
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="n_ele", fs_or_rs="rs", task=dict(
                title="RS comoving electron spectra evolution ratio",
                figname=f"spec_rs_ele_ratio_{name}", show=plot,
                n_ele=dict(norm_method=None,ylabel=r"$\gamma_{e}$",
                           zlabel=label, vmin=1e-1,vmax=1e1, cmap="RdBu_r",xlim=(3e3,8e6),set_under='blue',set_over='red',
                           norm="SymLogNorm",plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False)))
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="synch", fs_or_rs="rs", task=dict(
                title="RS comoving synchrotron spectra evolution ratio",
                figname=f"spec_rs_synch_ratio_{name}", show=plot,
                synch=dict(norm_method=None,ylabel=r"$\nu'$",
                           zlabel=label, vmin=1e-1,vmax=1e1, cmap="RdBu_r", xlim=(3e3,8e6),set_under='blue',set_over='red',
                           norm="SymLogNorm",
                           plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True)))
        plot_spectra_evolution_ratio(
            ej1=model1,ej2=model2,
            v_n="int", fs_or_rs="rs", task=dict(
                title="RS comoving intensity spectra evolution ratio",
                figname=f"spec_rs_int_ratio_{name}", show=plot,
                int=dict(norm_method=None,ylabel=r"$\nu'$",
                         zlabel=label, vmin=1e-1,vmax=1e1, cmap="RdBu_r", xlim=(3e3,8e6), ylim=(1e6,1e14),
                         norm="SymLogNorm", set_under='blue',set_over='red',
                         plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True)))



if __name__ == '__main__':
    do_run = False; plot = True
    struct = dict(struct="tophat",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    # struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    P = dict(
        main=dict(n_ism=1., tb0=3e3, ntb=1000,rtol=1e-7,theta_obs=0,
                  lc_freqs='array logspace 1e8 1e27 64',
                  lc_times='array logspace 3e3 1e10 128'),
        grb=dict(save_dynamics='yes',do_spec='yes',save_spec='yes',do_lc='yes',
                 method_nonrel_dist_fs='none',
                 method_nonrel_dist_rs='none',
                 eps_e_fs=0.1,eps_b_fs=1e-3,#gamma_max_fs=4e7,
                 # eps_b_fs = 1e-7,
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 freq1=1e6,freq2=1e30,nfreq=400)
    )

    # --- fs -- fs ---
    tasks_fs(do_run=do_run, plot=plot, struct=struct, P=P)
    # tasks_fs_comparison(do_run=do_run, plot=plot, struct=struct, P=P)

    # --- fsrs -- rs ---
    # tasks_rs(do_run=do_run, plot=plot, struct=struct, P=P)
    # tasks_rs_comparison(do_run=do_run, plot=plot,struct=struct, P=P)
