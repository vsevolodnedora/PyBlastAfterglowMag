import copy

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs,d2d
import os, shutil, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, SymLogNorm, TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator, AutoLocator
from matplotlib import ticker
import matplotlib.ticker as plticker
import numpy as np
from scipy import special
from scipy import integrate

from .interface import PyBlastAfterglow, Ejecta

def _normalize_spec(spec: np.ndarray, xs: np.ndarray, norm_method: None or str, mask_negatives=True):

    if (spec.ndim == 1):
        spec = spec[:, np.newaxis]

    if not norm_method:
        return spec
    if norm_method == "/integ *y^2":
        spec = spec / np.trapz(y=spec, x=xs, axis=0)[:, np.newaxis]
        spec *= np.power(xs, 2)[:, np.newaxis]
    elif norm_method == "/integ *y":
        integ = integrate.simps(y=spec, x=xs, axis=0)
        spec = spec / integ[np.newaxis, :]
        spec = spec * np.power(xs, 1)[:, np.newaxis]
    elif norm_method == "*y /integ":
        spec = spec * np.power(xs, 1)[np.newaxis, :]
        integ = integrate.simps(y=spec, x=xs, axis=0)
        spec = spec / integ[np.newaxis,:]
    elif norm_method == "*y":
        spec *= xs[:, np.newaxis]
    elif norm_method == "*y^2":
        spec *= xs[:, np.newaxis] ** 2
    elif norm_method == "/integ":
        # spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        integ = integrate.simps(y=spec, x=xs, axis=0)
        spec = spec / integ[np.newaxis,:]
    else:
        raise KeyError("Norm method is not recognized")
    if mask_negatives:
        spec[~np.isfinite(spec)] = 1e-100
        spec[spec <= 0] = 1e-100
    return spec

class PlotSpectra:
    qe = 4.803204e-10
    me = 9.1094e-28
    c = 2.99792458e10

    def __init__(self):
        pass

    @staticmethod
    def _get_spectrum(ej: PBA.Ejecta, v_n: str, fs_or_rs: str, norm_method: str,
                      ishell: int = 0, ilayer: int = 0, time: float or None =None, sum_shells_layers: bool = False):

        if v_n in ["n_ele", "gam_dot_syn", "gam_dot_adi", "gam_dot_ssc", "yg"]:
            xkey = "gams"
            key_time = "times_gams"
        else:
            xkey = "freqs"
            key_time = "times_freqs"

        if v_n == "fluxdens":
            is_spec = False
            xkey = "freqs"
            key_time = "times"
        else:
            is_spec = True

        ele_syn_ssc = v_n
        if v_n == "syn_tau": ele_syn_ssc = "syn_a"

        spec = ej.get_lc(key=ele_syn_ssc + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
                         xkey=xkey, key_time=key_time, freq=None, time=time,
                         ishell=ishell, ilayer=ilayer,
                         spec=is_spec, sum_shells_layers=sum_shells_layers)

        # if v_n == "total_f":
        #     spec_syn = ej.get_lc(key="syn_j" + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
        #                          xkey=xkey, key_time=key_time, freq=None, time=None,
        #                          ishell=ishell, ilayer=ilayer,
        #                          spec=is_spec, sum_shells_layers=sum_shells_layers)
        #     spec_ssc = ej.get_lc(key="ssc_j" + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
        #                          xkey=xkey, key_time=key_time, freq=None, time=None,
        #                          ishell=ishell, ilayer=ilayer,
        #                          spec=is_spec, sum_shells_layers=sum_shells_layers)
        #     spec = spec_syn + spec_ssc
        # elif v_n == "total_i":
        #     spec_syn = ej.get_lc(key="syn_i" + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
        #                          xkey=xkey, key_time=key_time, freq=None, time=None,
        #                          ishell=ishell, ilayer=ilayer,
        #                          spec=is_spec, sum_shells_layers=sum_shells_layers)
        #     spec_ssc = ej.get_lc(key="ssc_i" + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
        #                          xkey=xkey, key_time=key_time, freq=None, time=None,
        #                          ishell=ishell, ilayer=ilayer,
        #                          spec=is_spec, sum_shells_layers=sum_shells_layers)
        #     spec = spec_syn + spec_ssc
        # else:
        #     spec = ej.get_lc(key=ele_syn_ssc + '_' + fs_or_rs if ele_syn_ssc != "fluxdens" else "fluxdens",
        #                      xkey=xkey, key_time=key_time, freq=None, time=None,
        #                      ishell=ishell, ilayer=ilayer,
        #                      spec=is_spec, sum_shells_layers=sum_shells_layers)
        # if v_n == "Yg":
        #     spec = PlotSpectra._get_spectrum(ej=ej,v_n="total_f",fs_or_rs=fs_or_rs,norm_method=norm_method,
        #                                      ishell=ishell,ilayer=ilayer,sum_shells_layers=sum_shells_layers)

        if v_n == "syn_tau":
            dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thichness_rs", ishell=ishell, ilayer=ilayer)
            Gamma = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs == "fs" else "Gamma43", ishell=ishell, ilayer=ilayer)
            dr_comov = dr * Gamma  # thickness as R / Gamma in comoving frame (dr is in observer frame so R/Gamma^2
            spec = spec * dr_comov[np.newaxis, :]
        if (np.sum(spec) == 0):
            raise ValueError(f"v_n={v_n} fs_rs={fs_or_rs} ishell={ishell} ilayer={ilayer} (np.sum(spec) == 0)")
        xs = ej.get_grid(key=xkey, spec=is_spec)
        times = ej.get_grid(key=key_time, spec=is_spec)
        if not time is None:
            times = [times[PBA.utils.find_nearest_index(times, time)]]

        spec = _normalize_spec(spec=spec, xs=xs, norm_method=norm_method, mask_negatives=True)

        if not time is None:
            return xs, np.ndarray.flatten(spec)

        return (xs, times, spec)

    # @staticmethod
    # def _gam_to_nu(gam:np.ndarray,B:np.ndarray,p:float,syn_or_ssc:str,regime:str):
    #     """ Using Johanesson+2006 for XpS and Sari & Esin for SSC """
    #     gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
    #     if regime=="fast":
    #         XpS = 0.06 + 0.28 * p
    #     else:
    #         XpS = 0.455 + 0.08 * p
    #     nu = XpS * gam**2 * gamToNuFactor
    #
    #     if syn_or_ssc.__contains__("ssc"): return 4.* gam ** 2 * nu * np.sqrt(2)/3.
    #     return nu
    #     # nu_min = XpS * gamma_min * gamma_min * gamToNuFactor
    #     # nu_c = XpS * gamma_c * gamma_c * gamToNuFactor
    #     # nu_max = XpS * gamma_max * gamma_max * gamToNuFactor
    @staticmethod
    def _nu_m(gm: np.ndarray, gc: np.ndarray, B: np.ndarray, p: float, ssc: bool):
        qe = 4.803204e-10
        me = 9.1094e-28
        c = 2.99792458e10
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)
        XpS = np.full_like(gm, fill_value=1.)  # fast
        # XpS[gm<gc] = 0.06 + 0.28 * p # Fast
        # XpS[gm>gc] = 0.455 + 0.08 * p # Slow
        nu_m = XpS * gm ** 2 * gamToNuFactor
        if ssc:
            return 4. * gm ** 2 * nu_m * np.sqrt(2) / 3.  # Sari & Esin + 02
        else:
            return nu_m

    @staticmethod
    def _nu_c(gm: np.ndarray, gc: np.ndarray, B: np.ndarray, p: float, ssc: bool):
        qe = 4.803204e-10
        me = 9.1094e-28
        c = 2.99792458e10
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)  # qe * B / (2.*np.pi*me*c) # Ludovica
        XpS = np.full_like(gm, fill_value=1)  # fast
        # XpS[gm>gc] = 0.06 + 0.28 * p # Fast
        # XpS[gm<gc] = 0.455 + 0.08 * p # Slow
        nu_m = XpS * gc ** 2 * gamToNuFactor
        if ssc:
            return 4. * gc ** 2 * nu_m * np.sqrt(2) / 3.  # Sari & Esin + 02
        else:
            return nu_m

    @staticmethod
    def _nu_M(gM: np.ndarray, B: np.ndarray, p: float, ssc: bool):
        qe = 4.803204e-10
        me = 9.1094e-28
        c = 2.99792458e10
        gamToNuFactor = (3.0 / (4.0 * np.pi)) * (qe * B) / (me * c)  # qe * B / (2.*np.pi*me*c) # Ludovica
        nu_m = gM ** 2 * gamToNuFactor
        if ssc:
            return 4. * gM ** 2 * nu_m * np.sqrt(2) / 3.  # Sari & Esin + 02
        else:
            return nu_m

    @staticmethod
    def _pmax(nu_m:float,nu_c:float,B:float,p:float):
        c,pi,mp,me,kB = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28 , 1.38065e-16
        qe = 4.803204e-10
        phipF = 0.54 + 0.08*p
        phipS = 1.89 - 0.935*p + 0.17*p**2
        PmaxF = phipF * 2.234*qe**3*B/me/c**2
        PmaxS = phipS * 11.17*(p-1)*qe**3*B/(3*p-1)/me/c**2
        return (PmaxS if nu_m < nu_c else PmaxF)

    @staticmethod
    def _plot_ele_spectrum(ax, fig, times: np.ndarray, ys: np.ndarray, spec: np.ndarray,
                           ej: PBA.Ejecta, fs_or_rs: str, task_i: dict):
        # ys=PBA.utils.get_beta(ys)*ys
        # task_i = task[v_n]
        if task_i['norm'] == "LogNorm":
            norm = LogNorm(vmin=task_i.get("vmin", spec.max() * 1e-5),
                           vmax=task_i.get("vmax", spec.max()))
        elif task_i['norm'] == "SymLogNorm":
            norm = SymLogNorm(linthresh=task_i.get("vmin", 1e-4), vmin=task_i.get("vmin", 1e-4),
                              vmax=task_i.get("vmax", 1e4), base=10)
            # norm=TwoSlopeNorm(vmin=0.1,vcenter=1,vmax=10)
        else:
            raise KeyError(f"norm {task_i['norm']} is not recognized")

        if task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec > norm.vmax, spec)

        cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
        cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
        cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)
        # spec = np.ma.masked_where(spec < spec.max() * 1e-5, spec)
        if task_i["mode"] == "contour":
            _c = ax.contourf(times, ys, spec, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(times, ys, spec, cmap=cmap, norm=norm)
        cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label(task_i["zlabel"], size=12)
        # cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(n=6))
        # cbar.ax.yaxis.set_major_locator(AutoLocator())
        # cbar.ax.set_yticks(np.logspace(np.log10(spec.max() * 1e-6),np.log10(spec.max() * 10),5))
        # cbar.ax.set_yticks([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e1,1e2,1e3])
        ax.set_ylabel(task_i["ylabel"], fontsize=12)
        if task_i["plot_gm"]:
            gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=0, ilayer=0)
            ax.plot(times, gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        if task_i["plot_gc"]:
            gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=0, ilayer=0)
            ax.plot(times, gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        if task_i["plot_gM"]:
            gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=0, ilayer=0)
            ax.plot(times, gamma_max, color='black', linewidth=1, linestyle="-.", label=r"$\gamma_{\rm M}$")
        if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
        if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])
        return ax

    @staticmethod
    def nuprime_to_nu(nuprime: np.ndarray, ej: PBA.Ejecta):
        z = 0
        Gamma = ej.get_dyn_arr(v_n="Gamma", ishell=0, ilayer=0)
        val = nuprime * Gamma / (1 + z)
        return val

    @staticmethod
    def _plot_rad_spectrum(ax, fig, times: np.ndarray, ys: np.ndarray, spec: np.ndarray,
                           ej: PBA.Ejecta, fs_or_rs: str, v_n: str, task_i: dict):

        if task_i["plot_num"]:
            gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=0, ilayer=0)
            gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=0, ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=0, ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            val = PlotSpectra._nu_m(gm=gamma_min, gc=gamma_c, B=B, p=p, ssc=v_n.__contains__("ssc"))
            if v_n == "fluxdens":
                val = PlotSpectra.nuprime_to_nu(nuprime=val, ej=ej)
            ax.plot(times, val, color='black', linewidth=1, linestyle=":",
                    label=r"$\nu'_{\rm m}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm m;\, ssc}$")
        if task_i["plot_nuc"]:
            gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=0, ilayer=0)
            gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=0, ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=0, ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            val = PlotSpectra._nu_c(gm=gamma_min, gc=gamma_c, B=B, p=p, ssc=v_n.__contains__("ssc"))
            if v_n == "fluxdens":
                val = PlotSpectra.nuprime_to_nu(nuprime=val, ej=ej)
            ax.plot(times, val,
                    color='black', linewidth=1, linestyle="--",
                    label=r"$\nu'_{\rm c}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm c;\, ssc}$")
        if task_i["plot_nuM"]:
            gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=0, ilayer=0)
            B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=0, ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            val = PlotSpectra._nu_M(gM=gamma_max, B=B, p=p, ssc=v_n.__contains__("ssc"))
            if v_n == "fluxdens":
                val = PlotSpectra.nuprime_to_nu(nuprime=val, ej=ej)
            ax.plot(times, val,
                    color='black', linewidth=1, linestyle="-.",
                    label=r"$\nu'_{\rm M}$" if not v_n.__contains__("ssc") else r"$\nu'_{\rm M;\, ssc}$")
        if task_i["plot_nua"]:
            Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=0, ilayer=0)
            p = ej.pars[f"p_{fs_or_rs}"]
            nu_a = 3.41 * 1e9 / Gamma_sh * (p + 2) / (3. * p + 2) * (p - 1) ** (8 / 5) / (p - 2)  # From Warren:2018lyx

            theta = ej.get_dyn_arr(v_n="theta", ishell=0, ilayer=0)
            mask = np.array(Gamma_sh < 0.95 * Gamma_sh[0], dtype=int) * np.array(theta < 1.1 * theta[0], dtype=int)
            mask = np.array(mask, dtype=bool)
            if v_n == "fluxdens":
                nu_a = PlotSpectra.nuprime_to_nu(nuprime=nu_a[mask], ej=ej)
            else:
                nu_a = nu_a[mask]
            if np.sum(mask) > 0:
                ax.plot(times[mask], nu_a, color='black', linewidth=1, linestyle=":", label=r"$\nu'_{\rm a}$")
        if task_i["plot_tau1"]:
            Gamma = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=0, ilayer=0)
            ys, times, tau = PlotSpectra._get_spectrum(ej=ej, v_n="syn_tau", fs_or_rs=fs_or_rs, norm_method=None)
            freqs_tau_1_num = np.array([ys[PBA.utils.find_nearest_index(tau[:, i], 1.)] for i in range(len(Gamma))])
            ax.plot(
                times,
                freqs_tau_1_num,  # [xs,ys]
                color='black', lw=0.8, ls='-',
                label=r"$\tau_{\nu'}'=1$"
            )
        if task_i["plot_max"]:
            ax.plot(
                times,
                [ys[idx] for idx in np.argmax(spec, axis=0)],  # [xs,ys]
                color='gray', lw=0.9, ls='-',
                label=r"$\nu'_{\rm p}$"
            )

        if (task_i["plot_nua"] and task_i["plot_tau1"]):
            print(nu_a / freqs_tau_1_num[mask])

        # task_i = task[v_n]
        if task_i['norm'] == "LogNorm":
            norm = LogNorm(vmin=task_i.get("vmin", spec.max() * 1e-5), vmax=task_i.get("vmax", spec.max() * 2))
        elif task_i['norm'] == "SymLogNorm":
            norm = SymLogNorm(linthresh=task_i.get("vmin", 1e-4),
                              vmin=task_i.get("vmin", 1e-4),
                              vmax=task_i.get("vmax", 1e4), base=10)
        else:
            raise KeyError(f"norm {task_i['norm']} is not recognized")

        cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
        # _x = task_i.get("set_under", 'white')
        cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
        cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)

        if not v_n.__contains__("tau") and task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec > norm.vmax, spec)
            # if np.sum(np.asarray(spec,dtype='int')) :
            #     raise ValueError("nothing to plot")
            # _c = ax.pcolormesh(xs, ys, spec, cmap=cmap, norm=norm)

        if task_i["mode"] == "contour":
            _c = ax.contourf(times, ys, spec, cmap=cmap, locator=ticker.LogLocator(), norm=norm,extend='min')
            # _c.cmap.set_under(task_i.get("set_under", 'white'))
            # _c.changed()
        else:
            _c = ax.pcolormesh(times[::3], ys[::3], spec[::3,::3], cmap=cmap, norm=norm)
            # ax.set_rasterized(_c)

        cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label(task_i["zlabel"], size=12)
        # tick_locator = ticker.MaxNLocator(nbins=15)
        # cbar.locator = tick_locator
        # cbar.set_norm(norm)
        # cbar.update_normal(_c)
        # cbar.update_ticks()
        ax.set_ylabel(task_i["ylabel"], fontsize=12)

        if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
        if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])

        return ax


def plot_spectra_evolution(ej: PBA.Ejecta, fs_or_rs: str, tasks: dict, title: str or None, legend: dict or None,
                           figsize: tuple, figname: str or None, show: bool):
    n_tasks = len(tasks.keys())
    fig, axes = plt.subplots(ncols=1, nrows=n_tasks, sharex='all', layout='constrained', figsize=figsize)
    if not hasattr(axes, '__len__'):
        axes = [axes]
    i_plot = 0

    pl = PlotSpectra()

    for v_n, task in tasks.items():
        if v_n in ["n_ele", "gam_dot_syn", "gam_dot_adi", "gam_dot_ssc", "yg"]:
            ys, times, spec = pl._get_spectrum(ej=ej, v_n=v_n, fs_or_rs=fs_or_rs, norm_method=task["norm_method"])
            pl._plot_ele_spectrum(ax=axes[i_plot], fig=fig, times=times, ys=ys, spec=spec, ej=ej, fs_or_rs=fs_or_rs,
                                  task_i=task)
        else:
            ys, times, spec = pl._get_spectrum(ej=ej, v_n=v_n, fs_or_rs=fs_or_rs, norm_method=task["norm_method"])
            pl._plot_rad_spectrum(ax=axes[i_plot], times=times, ys=ys, spec=spec, fig=fig, ej=ej, fs_or_rs=fs_or_rs,
                                  v_n=v_n, task_i=task)
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
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        if legend: ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
                             # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                             shadow=False, ncol=4, fontsize=12, labelcolor='black',
                             framealpha=0.8, borderaxespad=0.)
        # ax.set_rasterized(True)
        # start, end = ax.get_xlim()
        # ax.xaxis.set_ticks(np.logspace(np.log10(start), np.log10(end), 9))
        # loc = plticker.MultipleLocator(base=100.0) # this locator puts ticks at regular intervals
        # ax.xaxis.set_major_locator(loc)
        # ax.locator_params(axis='y',nbins=10)
    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    axes[0].set_title(title, fontsize=12)
    if not figname is None:
        plt.savefig(os.getcwd() + '/figs/' + figname + '.png', dpi=256)
        plt.savefig(os.getcwd() + '/figs/' + figname + '.pdf')
    if show: plt.show()

    plt.close(fig)

