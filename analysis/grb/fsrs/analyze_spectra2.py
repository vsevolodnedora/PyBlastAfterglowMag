import copy

from more_itertools import numeric_range

import package.src.PyBlastAfterglowMag as PBA
import os, shutil, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, SymLogNorm, TwoSlopeNorm
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator, AutoLocator
from matplotlib import ticker
import matplotlib.ticker as plticker
import numpy as np
from scipy import special
from scipy import integrate

from package.src.PyBlastAfterglowMag import Ejecta

working_dir = os.getcwd() + '/tmp1/'
fig_dir = os.getcwd() + '/figs/'


"""
Papers:
THE SIGNATURE OF A WIND REVERSE SHOCK IN GAMMA-RAY BURST AFTERGLOWS (Asaf Pe’er1 and Ralph A. M. J. Wijers1)
https://iopscience.iop.org/article/10.1086/500969/pdf

GAMMA: a new method for modeling relativistic hydrodynamics and
non-thermal emission on a moving mesh (Eliot H. Ayache,1★ Hendrik J. van Eerten,1† Rupert W. Eardley 1)
https://arxiv.org/pdf/2104.09397.pdf

Probing particle acceleration at trans-relativistic shocks with off-axis gamma-ray burst afterglows 
(Kazuya Takahashi, Kunihito Ioka, Yutaka Ohira, Hendrik J van Eerten)
https://academic.oup.com/mnras/article/517/4/5541/6769826 (October 2022)

Gamma-Ray Burst Afterglows In The Multi-Messenger Era: Numerical Models and Closure Relations
Geoffrey Ryan,1 , ∗ Hendrik van Eerten,2 Luigi Piro,3 and Eleonora Troja4, 5
https://arxiv.org/pdf/1909.11691.pdf

A Semi-analytic Formulation for Relativistic Blast Waves with a Long-lived Reverse Shock
(Z. Lucas Uhm1,2)
https://arxiv.org/pdf/1003.1115.pdf (May 2011)

REFRESHED SHOCKS AND AFTERGLOW LONGEVITY IN GAMMA-RAY BURSTS
https://iopscience.iop.org/article/10.1086/311244/pdf (1998)
M. J. Rees & P. Me´sza´ros 

The origin of very-high-energy gamma-rays from GRB 221009A:
implications for reverse shock proton synchrotron emission
B. Theodore Zhang (张兵),1★ Kohta Murase,2,3,4,1 Kunihito Ioka,1 and Bing Zhang (张冰) 5,6
https://arxiv.org/pdf/2311.13671.pdf

GAMMA-RAY BURST REVERSE SHOCK EMISSION IN EARLY RADIO AFTERGLOWS
Lekshmi Resmi1 and Bing Zhang2
https://web.archive.org/web/20170508183725id_/http://www.physics.unlv.edu/~bzhang/rz16.pdf (2016)a

On the existence of a reverse shock in magnetized GRB ejecta
(D. Giannios1, P. Mimica2, and M. A. Aloy2
https://arxiv.org/pdf/0711.1980.pdf (2007)

SYNCHROTRON SELF-INVERSE COMPTON RADIATION FROM REVERSE SHOCK ON GRB 120326A
https://iopscience.iop.org/article/10.1088/0004-637X/789/2/146/pdf
Yuji Urata1, Kuiyun Huang2,3, Satoko Takahashi2,4,5, Myungshin Im6, Kazutaka Yamaoka7,8,
Makoto Tashiro9, Jae-Woo Kim6, Minsung Jang6, and Soojong Pak 
"""


def d2d(default: dict, new: dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_


def run(working_dir: str, struct: dict, P: dict, type: str = "a", run: bool = True) -> PBA.PyBlastAfterglow:
    # clean he temporary direcotry
    if run and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    # generate initial data for blast waves
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=80,
                                       n_layers_a=1 if struct["struct"] == "tophat" else 20)

    # save piece-wise EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_pw.h5")

    # save adaptive EATS ID
    id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
    pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_a.h5")

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
    if (run and pba.GRB.opts["do_skymap"] == "yes"):
        conf = {"nx": 128, "ny": 64, "extend_grid": 1.1, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                "intp_filter": {"type": 'gaussian', "sigma": 2, "mode": 'reflect'},  # "gaussian"
                "hist_filter": {"type": 'gaussian', "sigma": 2, "mode": 'reflect'}}
        prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=True)

    return pba


def _normalize_spec(spec: np.ndarray, xs: np.ndarray, times: np.ndarray, norm_method: None or str, mask_negatives=True):
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
        spec = spec / integ[:, np.newaxis]
    elif norm_method == "*y":
        spec *= xs[:, np.newaxis]
    elif norm_method == "*y^2":
        spec *= xs[:, np.newaxis] ** 2
    elif norm_method == "/integ":
        # spec = spec / np.trapz(y=spec, x=ys, axis=1)[:,np.newaxis]
        integ = integrate.simps(y=spec, x=xs, axis=0)
        spec = spec / integ[:, np.newaxis]
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
                      ishell: int = 0, ilayer: int = 0, sum_shells_layers: bool = False):

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
                         xkey=xkey, key_time=key_time, freq=None, time=None,
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

        spec = _normalize_spec(spec=spec, xs=xs, times=times, norm_method=norm_method, mask_negatives=True)

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
        cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
        cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)

        if not v_n.__contains__("tau") and task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec > norm.vmax, spec)
            # if np.sum(np.asarray(spec,dtype='int')) :
            #     raise ValueError("nothing to plot")
            # _c = ax.pcolormesh(xs, ys, spec, cmap=cmap, norm=norm)

        if task_i["mode"] == "contour":
            _c = ax.contourf(times, ys, spec, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(times, ys, spec, cmap=cmap, norm=norm)

        cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
        cbar.ax.tick_params(labelsize=11)
        cbar.set_label(task_i["zlabel"], size=12)
        # tick_locator = ticker.MaxNLocator(nbins=15)
        # cbar.locator = tick_locator
        # cbar.set_norm(norm)
        # cbar.update_normal(_c)
        # cbar.update_ticks()
        ax.set_ylabel(task_i["ylabel"], fontsize=12)
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
            ax.plot(
                times,
                [ys[PBA.utils.find_nearest_index(tau[:, i], 1.)] for i in range(len(Gamma))],  # [xs,ys]
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

        if "ylim" in task_i.keys(): ax.set_ylim(*task_i["ylim"])
        if "xlim" in task_i.keys(): ax.set_xlim(*task_i["xlim"])

        return ax


def plot_spectra_evolution(ej: PBA.Ejecta, fs_or_rs: str, tasks: dict, title: str or None,
                           figsize: tuple, figname: str, show: bool):
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
        ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        # start, end = ax.get_xlim()
        # ax.xaxis.set_ticks(np.logspace(np.log10(start), np.log10(end), 9))
        # loc = plticker.MultipleLocator(base=100.0) # this locator puts ticks at regular intervals
        # ax.xaxis.set_major_locator(loc)
        # ax.locator_params(axis='y',nbins=10)
    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    axes[0].set_title(title, fontsize=12)
    plt.savefig(os.getcwd() + '/figs/' + figname + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + figname + '.pdf')
    if show: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_spectra_evolution_ratio(ej1: PBA.Ejecta, ej2: PBA.Ejecta, v_n: str, fs_or_rs: str, task: dict):
    pl = PlotSpectra()

    fig, ax = plt.subplots(ncols=1, nrows=1, sharex='all', layout='constrained', figsize=(5.4, 3.0))
    task_i = task[v_n]
    ys, times, spec1 = pl._get_spectrum(ej=ej1, v_n=v_n, fs_or_rs=fs_or_rs, norm_method=task_i["norm_method"])
    ys, times, spec2 = pl._get_spectrum(ej=ej2, v_n=v_n, fs_or_rs=fs_or_rs, norm_method=task_i["norm_method"])

    spec = spec1 / spec2
    if v_n == "n_ele":
        pl._plot_ele_spectrum(ax=ax, fig=fig, times=times, ys=ys, spec=spec, ej=ej1, fs_or_rs=fs_or_rs,
                              task_i=task[v_n])
        # if task_i["plot_gm"]:
        #     gamma_min = ej1.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        # if task_i["plot_gc"]:
        #     gamma_c = ej1.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        # if task_i["plot_gM"]:
        #     gamma_max = ej1.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_max, color='black', linewidth=1, linestyle="-.", label=r"$\gamma_{\rm M}$")
        # if task_i["plot_gm"]:
        #     gamma_min = ej2.get_dyn_arr(v_n="gamma_min" if fs_or_rs=="fs" else "gamma_min_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_min, color='black', linewidth=1, linestyle=":", label=r"$\gamma_{\rm m}$")
        # if task_i["plot_gc"]:
        #     gamma_c = ej2.get_dyn_arr(v_n="gamma_c" if fs_or_rs=="fs" else "gamma_c_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_c, color='black', linewidth=1, linestyle="--", label=r"$\gamma_{\rm c}$")
        # if task_i["plot_gM"]:
        #     gamma_max = ej2.get_dyn_arr(v_n="gamma_max" if fs_or_rs=="fs" else "gamma_max_rs",ishell=0,ilayer=0)
        #     ax.plot(xs,gamma_max, color='black', linewidth=1, linestyle="-.", label=r"$\gamma_{\rm M}$")
    else:
        pl._plot_rad_spectrum(ax=ax, fig=fig, times=times, ys=ys, spec=spec, ej=ej1, fs_or_rs=fs_or_rs, v_n=v_n,
                              task_i=task[v_n])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=4, fontsize=12, labelcolor='black',
              framealpha=0.4, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.set_title(task["title"], fontsize=12)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_emission_region_prop(ej: PBA.Ejecta, fs_or_rs: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(5, 4.5), layout='constrained', sharex="all")
    #     double T = source.dr / CGS::c; // escape time (See Huang+2022)
    # double fac = T / source.vol;
    x = ej.get_dyn_arr(v_n="tburst", ishell=ishell, ilayer=ilayer)
    r = ej.get_dyn_arr(v_n="R", ishell=ishell, ilayer=ilayer)
    Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=ishell, ilayer=ilayer)
    dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thichness_rs", ishell=ishell, ilayer=ilayer)
    mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
    dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
    vol = mass / dens
    tau_tomps = PBA.utils.cgs.sigmaT * (mass / PBA.utils.cgs.mp / vol) * dr * Gamma_sh
    tmp = dr * Gamma_sh * PBA.utils.cgs.c / vol
    # ax.plot(
    #     x, PBA.utils.cgs.c/r**2
    # )

    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell, ilayer=ilayer)
    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell, ilayer=ilayer)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
    delta_t_syn = PBA.utils.cgs.sigmaT * gamma_max * gamma_max * B * B / (
                6. * np.pi * PBA.utils.cgs.me * PBA.utils.cgs.c)
    tcomov = ej.get_dyn_arr(v_n="tcomov", ishell=ishell, ilayer=ilayer)
    delta_t_adi = ((gamma_max * gamma_max - 1.) / (3. * gamma_max))[:-1] * \
                  ((tcomov[1:] - tcomov[:-1]) ** -1 * (1. - vol[:-1] / vol[1:]))

    sp = PlotSpectra()

    gams, times, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    # Energy density in relativistic particles https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf
    tot_ele = np.trapz(spec * gams[:, np.newaxis] / vol[np.newaxis, :] * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2,
                       x=gams, axis=0)

    gams, times, spec = sp._get_spectrum(ej=ej, v_n="syn_f", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_syn = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                       axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    sp = PlotSpectra()
    gams, times, spec = sp._get_spectrum(ej=ej, v_n="ssc_f", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_ssc = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                       axis=0)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    ax = axes[0]
    ax.plot(times, tot_ele, color='red', label=r"$u'_{\rm ele}$")
    ax.set_ylabel(r"Energy density of electrons", color='red', fontsize=12)
    # ax.legend()

    ax1 = ax.twinx()
    # ax1.plot(xs, tot_ssc/tot_syn,color='blue',label=r'$u_{\rm syn}$')
    ax1.plot(times, tot_syn, color='blue', label=r"$u'_{\rm syn}$")
    ax1.plot(times, tot_ssc, color='green', label=r"$u'_{\rm ssc}$")
    # ax1.plot(xs, tot_ele,color='red',label=r'$u_{\rm ele}$')
    ax1.set_ylabel(r"Energy density of radiation", color='black', fontsize=12)
    # ax1.legend()
    # ax1.plot(x, gamma_min)

    ax = axes[1]
    ax.plot(x, tmp, color='black', label=r"$\Delta t' / V'$")
    ax.set_ylabel(r"Escape time / V' ", color='black', fontsize=12)
    # ax.legend()

    ax2 = ax.twinx()
    ax2.plot(times, tau_tomps, color='magenta', label=r'$\tau_{\rm comp}$')
    ax2.set_ylabel(r"Compton optical depth", color='magenta', fontsize=12)
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

    for ax, loc in zip([axes[0], ax1, axes[1], ax2],
                       ["lower left", "upper right", "center left", "center right"]):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc=loc, columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)

    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    # ax.set_ylabel(r"$\tau_{\rm e}$")
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_emission_region_prop_energies(ej: PBA.Ejecta, fs_or_rs: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(5, 4.5), layout='constrained', sharex="all")

    mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
    dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
    vol = mass / dens

    Gamma = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs == "fs" else "Gamma43", ishell=ishell, ilayer=ilayer)
    Esh = ej.get_dyn_arr(v_n="Esh2" if fs_or_rs == "fs" else "Esh3", ishell=ishell, ilayer=ilayer)
    dEsh = np.diff(Esh)
    dEsh = np.insert(dEsh, 0, 0)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
    R = ej.get_dyn_arr(v_n="R" if fs_or_rs == "fs" else "R", ishell=ishell, ilayer=ilayer)
    deltaR = np.diff(R)
    deltaR = np.insert(deltaR, 0, 0)
    beta = PBA.utils.get_beta(Gamma)
    one_over_beta = 1. / beta
    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell,
                               ilayer=ilayer)
    gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=ishell, ilayer=ilayer)
    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell,
                               ilayer=ilayer)
    m = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
    tcomov = ej.get_dyn_arr(v_n="tcomov" if fs_or_rs == "fs" else "tcomov", ishell=ishell, ilayer=ilayer)
    dtcomov = np.diff(tcomov)
    dtcomov = np.insert(dtcomov, 0, 0)
    dm = np.diff(m)
    dm = np.insert(dm, 0, 0) / PBA.utils.cgs.mp
    p = float(ej.pars["p_fs" if fs_or_rs == "fs" else "p_rs"])
    eps_e = float(ej.pars["eps_e_fs" if fs_or_rs == "fs" else "eps_e_rs"])
    epsilon = np.zeros_like(B)
    dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thickness_rs", ishell=ishell, ilayer=ilayer)
    GammaSh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=ishell, ilayer=ilayer)

    drcomov = dr * GammaSh
    dt = drcomov / PBA.utils.cgs.c
    def get_rad_en(v_n):

        sp = PlotSpectra()

        gams, times, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                             ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        tot_ele = np.trapz(spec, x=gams, axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

        freqs, times, spec = sp._get_spectrum(ej=ej, v_n="syn_j", fs_or_rs=fs_or_rs,
                                             ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        # tot_syn = np.trapz(PBA.utils.cgs.h * freqs[:, np.newaxis] * spec, x=freqs, axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]
        tot_syn = np.trapz(spec, x=freqs, axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

        sp = PlotSpectra()
        freqs, times, spec = sp._get_spectrum(ej=ej, v_n="ssc_j", fs_or_rs=fs_or_rs,
                                             ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        # tot_ssc = np.trapz(PBA.utils.cgs.h * freqs[:, np.newaxis] * spec, x=freqs,axis=0)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]
        tot_ssc = np.trapz(spec, x=freqs, axis=0)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]
        # dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)

        tot = tot_syn + tot_ssc

        tot *= tot_ele
        # tot_syn *= dtcomov
        tot *= dt

        return tot*2.# * tot_ele / dens * PBA.utils.cgs.mp)



    def get_cont_rad_loss():
        Gammas = ej.get_dyn_arr(v_n="Gamma" if fs_or_rs == "fs" else "Gamma43", ishell=ishell, ilayer=ilayer)
        Bs = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
        Rs = ej.get_dyn_arr(v_n="R" if fs_or_rs == "fs" else "R", ishell=ishell, ilayer=ilayer)
        betas = PBA.utils.get_beta(Gammas)
        gamma_mins = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell,
                                    ilayer=ilayer)
        gamma_cs = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=ishell, ilayer=ilayer)
        gamma_maxs = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell,
                                    ilayer=ilayer)
        ms = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
        dmdr = np.diff(ms) / np.diff(Rs)
        n_eles = dmdr / PBA.utils.cgs.mp
        p = float(ej.pars["p_fs" if fs_or_rs == "fs" else "p_rs"])

        dErad2dR_cont = np.zeros_like(Gammas)
        i = -1
        for (B, gamma_min, gamma_c, gamma_max, Gamma, beta, n_ele) in \
                zip(Bs, gamma_mins, gamma_cs, gamma_maxs, Gammas, betas, n_eles):
            i += 1
            # if (gamma_c > gamma_min):  ### Slow cooling
            #     N0_total = n_ele * ((gamma_c ** (1 - p) - gamma_min ** (1 - p)) / (1 - p) - gamma_c * (gamma_max ** (-p) - gamma_c ** (-p)) / p)
            #     # N0_total = M2/mp*(p-1)*gamma_min**(p-1)
            #     # N0_total = Eint2 / c**2 /  me * (p-2) * gamma_min**(p-2) * ModVar.epsilone
            #     if gamma_c > gamma_max:
            #         dErad2dR_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B ** 2 / 6 / np.pi * (gamma_max ** (3 - p) - gamma_min ** (3 - p)) / (3 - p) / beta / PBA.utils.cgs.c ** 2 / Gamma
            #     else:
            #         dErad2dR_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B ** 2 / 6 / np.pi * ((gamma_c ** (3 - p) - gamma_min ** (3 - p)) / (3 - p) + (gamma_c * gamma_max ** (2 - p) - gamma_c ** (3 - p)) / (2 - p)) / beta / PBA.utils.cgs.c ** 2 / Gamma
            # else:  ### Fast cooling
            #
            #     N0_total = n_ele * p / (gamma_min ** (-p) - gamma_max ** (-p))
            #     # N0_total = M2 / mp * (gamma_min**(1-p)*(gamma_c_w**-1 - gamma_min**-1) + gamma_min**(-p)/p)
            #     # N0_total = M2 / mp * gamma_min**(-p)/p
            #     # dErad2dR += N0_total * sigmaT*B**2/6/pi * (gamma_min**(1-p)*(gamma_min-gamma_c_w) +
            #     # (gamma_max**(2-p)-gamma_min**(2-p))/(2-p)) / beta/c**2/Gamma
            #
            #     dErad2dR_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B ** 2 / 6 / np.pi * (gamma_max ** (2 - p) - gamma_min ** (2 - p)) / (2 - p) / beta / PBA.utils.cgs.c ** 2 / Gamma
            if gamma_c > gamma_min:  # Slow cooling
                N0_total = n_ele * ((gamma_c ** (1 - p) - gamma_min ** (1 - p)) / (1 - p) - gamma_c * (
                            gamma_max ** (-p) - gamma_c ** (-p)) / p)
                dErad2dR_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B ** 2 / 6 / PBA.utils.cgs.pi * (
                            (gamma_c ** (3 - p) - gamma_min ** (3 - p)) / (3 - p) + (
                                gamma_c * gamma_max ** (2 - p) - gamma_c ** (3 - p)) / (
                                        2 - p)) / beta / PBA.utils.cgs.c ** 2 / Gamma
            else:  # Fast cooling
                N0_total = n_ele * p / (gamma_min ** (-p) - gamma_max ** (-p))
                dErad2dR_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B ** 2 / 6 / PBA.utils.cgs.pi * (
                            gamma_max ** (2 - p) - gamma_min ** (2 - p)) / (2 - p) / beta / PBA.utils.cgs.c ** 2 / Gamma

        return dErad2dR_cont

    def i_eps():

        cooltime_integ = False
        remix_radloss = False
        continuous_radloss = True


        # gamma_c_w = gamma_c_w[i]
        #
        # V2 = M2[i] / (4 * Gamma[i] * rho[i])
        # B = np.sqrt(8 * np.pi * pars.eB * Eint2[i] / V2)
        Erad_cont = np.zeros_like(Gamma)
        for i in range(len(Gamma)-1):
            mue = PBA.utils.cgs.me / PBA.utils.cgs.mp

            # dE2sh_el = (Gamma[i] - 1) * dm[i] * PBA.utils.cgs.c ** 2 * eps_e
            dE2sh_el = (Esh[i+1]-Esh[i]) * eps_e

            # dn = dm[i] / PBA.utils.cgs.mp
            dn = (m[i+1]-m[i]) / PBA.utils.cgs.mp

            # gamma_max = (6 * np.pi * PBA.utils.cgs.qe / PBA.utils.cgs.sigmaT / B[i]) ** .5
            # gamma_min = (p - 2) / (p - 1) * (eps_e / mue * (Gamma[i] - 1) + 1)
            # if gamma_c[i] > 0:
            #     gmin_func = lambda gmin: self.gmin_fzero(gmin, self.gamma_c_w[i], gamma_max,
            #                                              p, self.pars.epsilone, self.Gamma[i] - 1, self.pars.mue)
            #     try:
            #         gamma_min = optimize.newton(gmin_func, gamma_min)
            #     except:
            #         res = optimize.fsolve(gmin_func, gamma_min)
            #         if len(res) > 1:
            #             raise ValueError('lenght of res > 1 !!!!!!!!!!!!!!')
            #             # raw_input(res)
            #         else:
            #             gamma_min = res[0]
            if gamma_c[i] > gamma_min[i]:
                dn_inj = dn / \
                         ((gamma_c[i] ** (-p - 1) - gamma_min[i] ** (-p + 1)) / (1 - p) -
                          (gamma_c[i] * gamma_max[i] ** (-p) - gamma_c[i] ** (-p + 1)) / p)
            else:
                dn_inj = dn / \
                         ((gamma_max[i] ** (1 - p) - gamma_min[i] ** (1 - p)) /
                          (1 - p) + gamma_min[i] ** (2 - p) / gamma_c[i] - gamma_min[i] ** (1 - p))
            gamma_c_mean = np.sqrt(gamma_min[i] * gamma_max[i])
            if gamma_c_mean > gamma_max[i]:
                gamma_c_mean = gamma_max[i]
            dErad2_ = 0.

            if cooltime_integ:
                if gamma_min[i] < gamma_max[i] and dE2sh_el > 0.:
                    dErad2_ += dn_inj * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2 * (
                            (gamma_c_mean ** (2 - p) - gamma_max[i] ** (2 - p)) / (3 - p) / (p - 2) -
                            (gamma_min[i] ** (3 - p) * (gamma_c_mean ** -1 - gamma_max[i] ** -1)) / (3 - p) +
                            (gamma_c_mean ** (2 - p) - gamma_max[i] ** (2 - p)) / (2 - p) ** 2 -
                            gamma_max[i] ** (2 - p) * np.log(gamma_c_mean / gamma_max[i]) / (2 - p))

            if remix_radloss and gamma_c[i] > 0:
                if gamma_min[i] > gamma_max[i]:
                    dErad2_ += dE2sh_el
                else:
                    if gamma_c[i] < gamma_max[i]:  ### No losses when gamma_c_w > gamma_max
                        dErad2_ += dn_inj * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2 * \
                                   (gamma_min[i] ** (2 - p) / (p - 2) - (
                                           gamma_c[i] * gamma_max[i] ** (1 - p) - gamma_c[i] ** (2 - p)
                                   )
                                    / (1 - p) - (gamma_c[i] ** (2 - p) - gamma_min[i] ** (2 - p)) / (2 - p))
                    else:
                        raise ValueError("val")


            if continuous_radloss and gamma_c[i] > 0:
                if gamma_c[i] > gamma_min[i]:  ### Slow cooling
                    N0_total = m[i] / PBA.utils.cgs.mp * (p - 1) * gamma_min[i] ** (p - 1)
                    # N0_total = Eint2 / c**2 /  me * (p-2) * gamma_min**(p-2) * epsilone
                    if gamma_c[i] > gamma_max[i]:
                        Erad_cont[i] + N0_total * PBA.utils.cgs.sigmaT * B[i] ** 2 * PBA.utils.cgs.c / 6 / np.pi * \
                                   (gamma_max[i] ** (3 - p) - gamma_min[i] ** (3 - p)) / \
                                   (3 - p) #* dtcomov[i]  # deltaR[i] * one_over_beta[i] / PBA.utils.cgs.c / Gamma[i]
                    else:
                        Erad_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B[i] ** 2 * PBA.utils.cgs.c / 6 / np.pi * (
                                (gamma_c[i] ** (3 - p) - gamma_min[i] ** (3 - p)) / (3 - p) +
                                (gamma_c[i] * gamma_max[i] ** (2 - p) - gamma_c[i] ** (3 - p)) /
                                (2 - p)) #* dtcomov[i]  # deltaR[i] * one_over_beta[i] / PBA.utils.cgs.c / Gamma[i]
                else:  ### Fast cooling
                    N0_total = m[i] / PBA.utils.cgs.mp * p / (gamma_min[i] ** (-p) - gamma_max[i] ** (-p))
                    Erad_cont[i] = N0_total * PBA.utils.cgs.sigmaT * B[i] ** 2 / 6 / np.pi * \
                               (gamma_max[i] ** (2 - p) - gamma_min[i] ** (2 - p)) / \
                               (2 - p) #* dtcomov[i]  # deltaR[i] * one_over_beta[i] / PBA.utils.cgs.c / Gamma[i]
                if (i > 0):
                    dErad2_ += Erad_cont[i]
                    # raise ValueError("[Fast cooling] not implemented")
            if gamma_c_mean > gamma_max[i]:
                print('gamma_c_mean larger than gamma_max!')
            if dErad2_ < 0:
                raise ValueError("dErad2 < 0:")
            epsilon[i] = dErad2_ / dE2sh_el
        #print("\teps_rad[{}]: {}".format(i, self.epsilon[i]))
        # dErad2[i] = dErad2_
        # assert epsilon_rad[i] > -1e-5 and epsilon_rad[i] < 1.01


        return epsilon

        #self.dErad2[i] = self.epsilon[i] * self.pars.epsilone * (self.Gamma[i] - 1) * self.dM2[i] * cgs.c ** 2

    # energies


    dt = dr * GammaSh / PBA.utils.cgs.c

    x = ej.get_dyn_arr(v_n="tburst", ishell=ishell, ilayer=ilayer)
    axes[0].plot(x, ej.get_dyn_arr(v_n="Eint2" if fs_or_rs == "fs" else "Eint3", ishell=ishell, ilayer=ilayer) / vol,
                 color='black', ls='-')
    axes[0].plot(x, ej.get_dyn_arr(v_n="Esh2" if fs_or_rs == "fs" else "Esh3", ishell=ishell, ilayer=ilayer) / vol,
                 color='black', ls='--')
    tmp = ej.get_dyn_obj()["shell=0 layer=0"].attrs["M0"] * PBA.utils.cgs.c ** 2
    axes[0].plot(x,
                 -ej.get_dyn_arr(v_n="Ead2" if fs_or_rs == "fs" else "Ead3", ishell=ishell, ilayer=ilayer) / tmp / vol,
                 color='gray', ls='--')
    axes[0].plot(x, get_rad_en('total'), color='black', ls=':')
    # axes[0].plot(x, get_cont_rad_loss() , color='black',ls='-.')

    axes[1].plot(x, i_eps(), color='black', ls='--')
    eps_e = float(ej.pars["eps_e_fs" if fs_or_rs == "fs" else "eps_e_rs"])
    ue = ej.get_dyn_arr(v_n="Esh2" if fs_or_rs == "fs" else "Esh3", ishell=ishell, ilayer=ilayer) * eps_e
    u_rad = get_rad_en('total')
    axes[1].plot(x, u_rad / ue, color='black', ls=':')
    axes[1].set_ylim(1e-3, 2)
    axes[1].axhline(y=1., color='gray',linestyle='dotted')

    for i, ax in enumerate(axes):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc='best', columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)
    plt.show()

    def get_rad_en(v_n):
        sp = PlotSpectra()
        gams, times, spec = sp._get_spectrum(ej=ej, v_n="syn_f", fs_or_rs=fs_or_rs,
                                             ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        tot_syn = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                           axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

        sp = PlotSpectra()
        gams, times, spec = sp._get_spectrum(ej=ej, v_n="ssc_f", fs_or_rs=fs_or_rs,
                                             ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        tot_ssc = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                           axis=0)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    #     double T = source.dr / CGS::c; // escape time (See Huang+2022)
    # double fac = T / source.vol;
    x = ej.get_dyn_arr(v_n="tburst", ishell=ishell, ilayer=ilayer)
    r = ej.get_dyn_arr(v_n="R", ishell=ishell, ilayer=ilayer)
    Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=ishell, ilayer=ilayer)
    dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thichness_rs", ishell=ishell, ilayer=ilayer)
    mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
    dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
    vol = mass / dens
    tau_tomps = PBA.utils.cgs.sigmaT * (mass / PBA.utils.cgs.mp / vol) * dr * Gamma_sh
    tmp = dr * Gamma_sh * PBA.utils.cgs.c / vol
    # ax.plot(
    #     x, PBA.utils.cgs.c/r**2
    # )

    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell, ilayer=ilayer)
    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell, ilayer=ilayer)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
    delta_t_syn = PBA.utils.cgs.sigmaT * gamma_max * gamma_max * B * B / (
                6. * np.pi * PBA.utils.cgs.me * PBA.utils.cgs.c)
    tcomov = ej.get_dyn_arr(v_n="tcomov", ishell=ishell, ilayer=ilayer)
    delta_t_adi = ((gamma_max * gamma_max - 1.) / (3. * gamma_max))[:-1] * \
                  ((tcomov[1:] - tcomov[:-1]) ** -1 * (1. - vol[:-1] / vol[1:]))

    sp = PlotSpectra()

    gams, times, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    # Energy density in relativistic particles https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf
    tot_ele = np.trapz(spec * gams[:, np.newaxis] / vol[np.newaxis, :] * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2,
                       x=gams, axis=0)

    gams, times, spec = sp._get_spectrum(ej=ej, v_n="syn_f", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_syn = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                       axis=0)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    sp = PlotSpectra()
    gams, times, spec = sp._get_spectrum(ej=ej, v_n="ssc_f", fs_or_rs=fs_or_rs,
                                         ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_ssc = np.trapz(PBA.utils.cgs.h * spec * gams[:, np.newaxis], x=gams,
                       axis=0)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    ax = axes[0]
    ax.plot(times, tot_ele, color='red', label=r"$u'_{\rm ele}$")
    ax.set_ylabel(r"Energy density of electrons", color='red', fontsize=12)
    # ax.legend()

    ax1 = ax.twinx()
    # ax1.plot(xs, tot_ssc/tot_syn,color='blue',label=r'$u_{\rm syn}$')
    ax1.plot(times, tot_syn, color='blue', label=r"$u'_{\rm syn}$")
    ax1.plot(times, tot_ssc, color='green', label=r"$u'_{\rm ssc}$")
    # ax1.plot(xs, tot_ele,color='red',label=r'$u_{\rm ele}$')
    ax1.set_ylabel(r"Energy density of radiation", color='black', fontsize=12)
    # ax1.legend()
    # ax1.plot(x, gamma_min)

    ax = axes[1]
    ax.plot(x, tmp, color='black', label=r"$\Delta t' / V'$")
    ax.set_ylabel(r"Escape time / V' ", color='black', fontsize=12)
    # ax.legend()

    ax2 = ax.twinx()
    ax2.plot(times, tau_tomps, color='magenta', label=r'$\tau_{\rm comp}$')
    ax2.set_ylabel(r"Compton optical depth", color='magenta', fontsize=12)
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

    for ax, loc in zip([axes[0], ax1, axes[1], ax2],
                       ["lower left", "upper right", "center left", "center right"]):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc=loc, columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)

    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    # ax.set_ylabel(r"$\tau_{\rm e}$")
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_emission_region_prop_3(ej: PBA.Ejecta, fs_or_rs: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    sp = PlotSpectra()
    fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(5, 4.5), layout='constrained', sharex="all")

    # ej = pba.GRB

    def get_total_ele_n(ishell, ilayer, fs_or_rs):
        mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
        xs, ys, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        total_ele = integrate.simps(y=spec, x=xs, axis=0)
        ksi_dn = ej.get_dyn_arr(v_n="accel_frac" if fs_or_rs == "fs" else "accel_frac_rs", ishell=ishell, ilayer=ilayer)
        return total_ele / (mass / PBA.utils.cgs.mp)

    def get_total_ele_en(ishell, ilayer, fs_or_rs):
        mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
        dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
        vol = mass / dens
        xs, ys, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
        # Energy density in relativistic particles https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf
        tot_ele = np.trapz(spec * xs[:, np.newaxis] * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2, x=xs, axis=0)

        ue = ej.get_dyn_arr(v_n="U_p" if fs_or_rs == "fs" else "U_p3", ishell=ishell, ilayer=ilayer)
        eps_e = ej.pars["eps_e_fs"] if fs_or_rs == "fs" else ej.pars["eps_e_rs"]
        ksi_dn = ej.get_dyn_arr(v_n="accel_frac" if fs_or_rs == "fs" else "accel_frac_rs", ishell=ishell, ilayer=ilayer)
        return tot_ele / (eps_e * ue * vol / ksi_dn)

    x = ej.get_dyn_arr(v_n="tburst", ishell=ishell, ilayer=ilayer)

    # plot ele number
    ax = axes[0]
    for i, il in enumerate(task["layers"]):
        ax.plot(x, get_total_ele_n(0, il, "fs"), color=task["colors"][i], ls='-')
        ax.plot(x, get_total_ele_n(0, il, "rs"), color=task["colors"][i], ls='--')

    # plot ele energy
    ax = axes[1]
    for i, il in enumerate(task["layers"]):
        ax.plot(x, get_total_ele_en(0, il, "fs"), color=task["colors"][i], ls='-')
        ax.plot(x, get_total_ele_en(0, il, "rs"), color=task["colors"][i], ls='--')

    # plot xi
    ax = axes[2]
    for i, il in enumerate(task["layers"]):
        ax.plot(x, ej.get_dyn_arr(v_n="accel_frac", ishell=0, ilayer=il), color=task["colors"][i], ls='-')
        ax.plot(x, ej.get_dyn_arr(v_n="accel_frac_rs", ishell=0, ilayer=il), color=task["colors"][i], ls='--')

    for i, ax in enumerate(axes):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc="best", columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)
        ax.grid(ls=':')
        if (f"ylim{i}" in task.keys()): ax.set_ylim(*task[f"ylim{i}"])

    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    # ax.set_ylabel(r"$\tau_{\rm e}$")
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)

    #     double T = source.dr / CGS::c; // escape time (See Huang+2022)
    # double fac = T / source.vol;
    x = ej.get_dyn_arr(v_n="tburst", ishell=ishell, ilayer=ilayer)
    r = ej.get_dyn_arr(v_n="R", ishell=ishell, ilayer=ilayer)
    Gamma_sh = ej.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=ishell, ilayer=ilayer)
    dr = ej.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thichness_rs", ishell=ishell, ilayer=ilayer)
    mass = ej.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
    dens = ej.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
    vol = mass / dens
    tau_tomps = PBA.utils.cgs.sigmaT * (mass / PBA.utils.cgs.mp / vol) * dr * Gamma_sh
    tmp = dr * Gamma_sh * PBA.utils.cgs.c / vol
    # ax.plot(
    #     x, PBA.utils.cgs.c/r**2
    # )

    gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell, ilayer=ilayer)
    gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell, ilayer=ilayer)
    B = ej.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
    delta_t_syn = PBA.utils.cgs.sigmaT * gamma_max * gamma_max * B * B / (
                6. * np.pi * PBA.utils.cgs.me * PBA.utils.cgs.c)
    tcomov = ej.get_dyn_arr(v_n="tcomov", ishell=ishell, ilayer=ilayer)
    delta_t_adi = ((gamma_max * gamma_max - 1.) / (3. * gamma_max))[:-1] * \
                  ((tcomov[1:] - tcomov[:-1]) ** -1 * (1. - vol[:-1] / vol[1:]))

    sp = PlotSpectra()

    xs, ys, spec = sp._get_spectrum(ej=ej, v_n="n_ele", fs_or_rs=fs_or_rs,
                                    ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    # Energy density in relativistic particles https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf
    tot_ele = np.trapz(spec * ys[np.newaxis, :] / vol[:, np.newaxis] * PBA.utils.cgs.me * PBA.utils.cgs.c ** 2, x=ys,
                       axis=1)

    xs, ys, spec = sp._get_spectrum(ej=ej, v_n="syn_f", fs_or_rs=fs_or_rs,
                                    ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_syn = np.trapz(PBA.utils.cgs.h * spec * ys[np.newaxis, :], x=ys,
                       axis=1)  # * PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    sp = PlotSpectra()
    xs, ys, spec = sp._get_spectrum(ej=ej, v_n="ssc_f", fs_or_rs=fs_or_rs,
                                    ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    tot_ssc = np.trapz(PBA.utils.cgs.h * spec * ys[np.newaxis, :], x=ys,
                       axis=1)  # PBA.utils.cgs.h*spec*ys[np.newaxis,:]

    ax = axes[0]
    ax.plot(xs, tot_ele, color='red', label=r"$u'_{\rm ele}$")
    ax.set_ylabel(r"Energy density of electrons", color='red', fontsize=12)
    # ax.legend()

    ax1 = ax.twinx()
    # ax1.plot(xs, tot_ssc/tot_syn,color='blue',label=r'$u_{\rm syn}$')
    ax1.plot(xs, tot_syn, color='blue', label=r"$u'_{\rm syn}$")
    ax1.plot(xs, tot_ssc, color='green', label=r"$u'_{\rm ssc}$")
    # ax1.plot(xs, tot_ele,color='red',label=r'$u_{\rm ele}$')
    ax1.set_ylabel(r"Energy density of radiation", color='black', fontsize=12)
    # ax1.legend()
    # ax1.plot(x, gamma_min)

    ax = axes[1]
    ax.plot(x, tmp, color='black', label=r"$\Delta t' / V'$")
    ax.set_ylabel(r"Escape time / V' ", color='black', fontsize=12)
    # ax.legend()

    ax2 = ax.twinx()
    ax2.plot(xs, tau_tomps, color='magenta', label=r'$\tau_{\rm comp}$')
    ax2.set_ylabel(r"Compton optical depth", color='magenta', fontsize=12)
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

    for ax, loc in zip([axes[0], ax1, axes[1], ax2],
                       ["lower left", "upper right", "center left", "center right"]):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc=loc, columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.4, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_xlim(*xlim)

    axes[-1].set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    # ax.set_ylabel(r"$\tau_{\rm e}$")
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_total_observed_spectrum_OLD(do_run: bool, norm_method: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    ej_fs = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num__ssa__ssc__fluxdens/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          do_rs_radiation="no",
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )
    ej_fsrs = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_fsrs__num__ssa__ssc__fluxdens/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )

    task_i = task
    sp = PlotSpectra()
    xs, ys, spec_fs = sp._get_spectrum(
        ej=ej_fs.GRB, v_n="fluxdens", fs_or_rs="rs",
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=norm_method)
    xs, ys, spec_fsrs = sp._get_spectrum(
        ej=ej_fsrs.GRB, v_n="fluxdens", fs_or_rs="rs",
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=norm_method)

    tmp = spec_fsrs / spec_fs

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 3), layout='constrained')
    norm = LogNorm(vmin=task_i.get("vmin", spec_fsrs.max() * 1e-20),
                   vmax=task_i.get("vmax", spec_fsrs.max() * 1))
    cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
    cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
    cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)

    spec = spec_fs;  #spec_fsrs - spec_fs

    if task_i["mode"] == "contour":
        spec = np.ma.masked_where(spec < norm.vmin, spec)
        spec = np.ma.masked_where(spec_fsrs > norm.vmax, spec)
        # _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)

    if task_i["mode"] == "contour":
        # _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
    else:
        _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"], size=12)

    # use contourf() with proper hatch pattern and alpha value
    val1 = np.array((spec_fsrs > 1.1 * spec_fs), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(xs, ys, val1.T, hatches=['..'], alpha=0.0, color="green")  #,zorder=0
    #
    # val1 = np.array((gam_dot_adi>gam_dot_syn) * (gam_dot_adi>gam_dot_ssc), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=["\\\\"],  alpha=0.0) #,zorder=0
    #
    # val1 = np.array((gam_dot_ssc>gam_dot_syn) * (gam_dot_ssc>gam_dot_adi), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=['++'],  alpha=0.0, edgecolor = "r",facecolor="blue", color = "green") #,zorder=0

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=4, fontsize=12, labelcolor='black',
              framealpha=0.8, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_ylabel(r"$\gamma_{e}$", fontsize=12)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.set_title(task["title"], fontsize=12)
    ax.set_xlim(*xlim)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_total_observed_spectrum(do_run: bool, norm_method: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    ej_fs = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num__ssa__ssc__fluxdens/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          do_rs_radiation="no",
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )
    ej_fsrs = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_fsrs__num__ssa__ssc__fluxdens/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )

    def plot_one(ax, xs, ys, spec, norm, cmap, text):
        if task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec_fsrs > norm.vmax, spec)
            # _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)

        if task_i["mode"] == "contour":
            # _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
            _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        ax.text(0.95, 0.07, text, fontsize=12, bbox=bbox,
                transform=ax.transAxes, horizontalalignment='right')
        return _c

    task_i = task
    sp = PlotSpectra()
    xs, ys, spec_fs = sp._get_spectrum(
        ej=ej_fs.GRB, v_n="fluxdens", fs_or_rs="rs",
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=norm_method)
    xs, ys, spec_fsrs = sp._get_spectrum(
        ej=ej_fsrs.GRB, v_n="fluxdens", fs_or_rs="rs",
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=norm_method)

    spec_rs = spec_fsrs - spec_fs

    fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(5, 6), layout='constrained', sharex='all', sharey='all')
    norm = LogNorm(vmin=task_i.get("vmin", spec_fsrs.max() * 1e-20),
                   vmax=task_i.get("vmax", spec_fsrs.max() * 1))
    cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
    cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
    cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)

    _c = plot_one(axes[0], xs, ys, spec_fs, norm, cmap, r"FS")
    _c = plot_one(axes[1], xs, ys, spec_rs, norm, cmap, r"RS")
    _c = plot_one(axes[2], xs, ys, spec_fsrs, norm, cmap, r"FS \& RS")

    cbar = fig.colorbar(_c, ax=axes, shrink=0.95, pad=0.01, aspect=30, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"], size=12)

    # use contourf() with proper hatch pattern and alpha value
    val1 = np.array((spec_rs > spec_fs), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = axes[-1].contourf(xs, ys, val1.T, hatches=['..'], alpha=0.0, color="green")  #,zorder=0
    #
    # val1 = np.array((gam_dot_adi>gam_dot_syn) * (gam_dot_adi>gam_dot_ssc), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=["\\\\"],  alpha=0.0) #,zorder=0
    #
    # val1 = np.array((gam_dot_ssc>gam_dot_syn) * (gam_dot_ssc>gam_dot_adi), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=['++'],  alpha=0.0, edgecolor = "r",facecolor="blue", color = "green") #,zorder=0
    for ax in axes:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.8, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_ylabel(task_i["ylabel"], fontsize=12)
        # ax.set_title(task["title"], fontsize=12)
        ax.set_xlim(*xlim)
    axes[-1].set_xlabel(task_i["xlabel"], fontsize=12)

    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def plot_cooling_terms(ej: PBA.Ejecta, fs_or_rs: str, xlim: tuple, task: dict, ishell=0, ilayer=0):
    task_i = task
    sp = PlotSpectra()
    gams, times, gam_dot_syn = sp._get_spectrum(
        ej=ej, v_n="gam_dot_syn", fs_or_rs=fs_or_rs,
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    gams, times, gam_dot_adi = sp._get_spectrum(
        ej=ej, v_n="gam_dot_adi", fs_or_rs=fs_or_rs,
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    gams, times, gam_dot_ssc = sp._get_spectrum(
        ej=ej, v_n="gam_dot_ssc", fs_or_rs=fs_or_rs,
        ishell=ishell, ilayer=ilayer, sum_shells_layers=False, norm_method=None)
    gam_tot = gam_dot_syn + gam_dot_adi + gam_dot_ssc

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 3), layout='constrained')
    norm = LogNorm(vmin=gam_tot.max() * 1e-20, vmax=gam_tot.max() * 1)
    cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
    cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
    cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)
    _c = ax.pcolormesh(times, gams, gam_tot, cmap=cmap, norm=norm)
    cbar = fig.colorbar(_c, ax=ax, shrink=0.95, pad=0.01, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"], size=12)

    # use contourf() with proper hatch pattern and alpha value
    val1 = np.array((gam_dot_syn > gam_dot_adi) * (gam_dot_syn > gam_dot_ssc), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(times, gams, val1, hatches=['//'], alpha=0.0, color="green")  #,zorder=0

    val1 = np.array((gam_dot_adi > gam_dot_syn) * (gam_dot_adi > gam_dot_ssc), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(times, gams, val1, hatches=["\\\\"], alpha=0.0)  #,zorder=0

    val1 = np.array((gam_dot_ssc > gam_dot_syn) * (gam_dot_ssc > gam_dot_adi), dtype=np.float64)
    val1[val1 < 1] = np.nan
    cs = ax.contourf(times, gams, val1, hatches=['++'], alpha=0.0, edgecolor="r", facecolor="blue",
                     color="green")  #,zorder=0

    if task_i["plot_gm"]:
        gamma_min = ej.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=0, ilayer=0)
        ax.plot(times, gamma_min, color='blue', linewidth=2, linestyle=":", label=r"$\gamma_{\rm m}$")
    if task_i["plot_gc"]:
        gamma_c = ej.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=0, ilayer=0)
        ax.plot(times, gamma_c, color='blue', linewidth=2, linestyle="--", label=r"$\gamma_{\rm c}$")
    if task_i["plot_gM"]:
        gamma_max = ej.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=0, ilayer=0)
        ax.plot(times, gamma_max, color='blue', linewidth=2, linestyle="-.", label=r"$\gamma_{\rm M}$")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=4, fontsize=12, labelcolor='black',
              framealpha=0.8, borderaxespad=0.)
    ax.set_rasterized(True)
    ax.set_ylabel(r"$\gamma_{e}$", fontsize=12)
    ax.set_xlabel(r"$t_{\rm burst}$ [s]", fontsize=12)
    ax.set_title(task["title"], fontsize=12)
    ax.set_xlim(*xlim)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def tasks_fs(do_run: bool, plot: bool, struct: dict, P: dict,
             name: str = "dyn_fs__rad_fs__num__ssa__ssc") -> PBA.PyBlastAfterglow:
    dyn_fs__rad_fs__num__ssa__ssc = run(
        working_dir=os.getcwd() + f"/working_dirs/{name}/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(method_ssc_fs='numeric', use_ssa_fs='yes'))),
        type="a", run=do_run
    )
    # raise ValueError("Computed!!! :) ")

    # plot_cooling_terms(dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs="fs",xlim=(3e3,1e9),
    #                    task=dict(show=True,figname="ele_cool_terms",zlabel=r"$\dot{\gamma}$",title="Cooling terms",
    #                              set_under='blue',set_over='red',
    #                              plot_gm=True,plot_gc=True,plot_gM=True))

    # plot_emission_region_prop(dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs="fs",xlim=(3e3,1e9),
    #                           task=dict(show=True,figname="rad_ele_energy_dens_fs"))
    plot_emission_region_prop_energies(dyn_fs__rad_fs__num__ssa__ssc.GRB, fs_or_rs="fs", xlim=(3e3, 1e9),
                                       task=dict(show=True, figname="rad_ele_energy_dens_fs"))

    # plot_spectra_evolution(
    #     ej=dyn_fs__rad_fs__num__ssa__ssc.GRB, fs_or_rs='fs',
    #     title="FS comoving spectra evolution", figname="spec_dyn_fs__rad_fs__num__ssa__ssc", show=plot,
    #     figsize=(5, 8),
    #     tasks=dict(
    #         n_ele=dict(norm_method='/integ *y',  #'/integ *y',
    #                    ylabel=r"$\gamma_{e}$",
    #                    zlabel=r"$(\gamma_{e} N_{e}) / N_{e;\, \rm tot}$",
    #                    norm="LogNorm", mode="contour", xlim=(3e3, 1e9), vmin=1e-6, vmax=1e0,
    #                    plot_gm=True, plot_gc=True, plot_gM=True, plot_max=False),
    #         # yg=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$Y_g$ [cgs]",
    #         #         norm="LogNorm",
    #         #         plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False),
    #         syn_j=dict(norm_method="*y", ylabel=r"$\nu'$ [Hz]", zlabel=r"$\nu' j'_{\rm syn}$ [cgs]",
    #                    norm="LogNorm", ylim=(1e12, 1e28), xlim=(3e3, 1e9), vmin=1e-9, vmax=1e-3, mode="contour",
    #                    plot_num=True, plot_nuc=True, plot_nuM=True, plot_nua=False, plot_tau1=False, plot_max=True),
    #         ssc_j=dict(norm_method="*y", ylabel=r"$\nu'$ [Hz]", zlabel=r"$\nu' j'_{\rm ssc}$ [cgs]",
    #                    norm="LogNorm", ylim=(1e18, 1e30), xlim=(3e3, 1e9), vmin=1e-9, vmax=1e-3, mode="contour",
    #                    plot_num=True, plot_nuc=True, plot_nuM=False, plot_nua=False, plot_tau1=False, plot_max=True),
    #         # syn_a=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\alpha'_{\rm syn}$ [cgs]",
    #         #          norm="LogNorm",mode="contour",
    #         #          plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
    #         #          ylim=(1e6,1e12), xlim=(3e3,1e9)),
    #         syn_tau=dict(norm_method=None, ylabel=r"$\nu'$ [Hz]", zlabel=r"$\tau'_{\rm ssa}$",
    #                      norm="SymLogNorm", mode=None,
    #                      plot_num=False, plot_nuc=False, plot_nuM=False, plot_nua=True, plot_tau1=True, plot_max=False,
    #                      ylim=(1e6, 1e14), xlim=(3e3, 1e9)),
    #         total_j=dict(norm_method="*y", ylabel=r"$\nu'$ [Hz]", zlabel=r"$\nu' I'_{\rm total}$ [cgs]",
    #                      norm="LogNorm", mode="contour",
    #                      plot_num=False, plot_nuc=False, plot_nuM=False, plot_nua=False, plot_tau1=False,
    #                      plot_max=False,
    #                      ylim=(1e12, 1e28), xlim=(3e3, 1e9), vmin=1e-9, vmax=1e-3),
    #                      # vmin=1e-15, vmax=1e-3),  # vmin=1e-17,vmax=1e-3
    #     ))
    plot_spectra_evolution(
        ej=dyn_fs__rad_fs__num__ssa__ssc.GRB,fs_or_rs='fs',
        title="FS comoving spectra evolution", figname="spec_dyn_fs__rad_fs__num__ssa__ssc", show=plot,
        figsize=(6,4),
        tasks=dict(
            fluxdens=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' F'_{\rm total}$ [cgs]",
                         norm="LogNorm",mode="contour",
                         plot_num=False,plot_nuc=False,plot_nuM=False ,plot_nua=False,plot_tau1=False,plot_max=False,
                         ylim=(1e7,1e28), xlim=(3e3,1e9), vmin=1e10,vmax=1e21),# vmin=1e-17,vmax=1e-3
        ))
    return dyn_fs__rad_fs__num__ssa__ssc


def tasks_rs(do_run: bool, plot: bool, struct: dict, P: dict) -> PBA.PyBlastAfterglow:
    dyn_fs__rad_rs__num__ssa__ssc = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_rs__num__ssa__ssc/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )

    # plot_cooling_terms(dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs="rs",xlim=(3e3,3e7),
    #                    task=dict(show=True,figname="ele_cool_terms_rs",zlabel=r"$\dot{\gamma}$",title="Cooling terms",
    #                              set_under='blue',set_over='red',
    #                              plot_gm=True,plot_gc=True,plot_gM=True))

    # plot_emission_region_prop(dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs="rs",xlim=(3e3,3e7),
    #                           task=dict(show=True,figname="rad_ele_energy_dens_rs"))
    #
    # plot_spectra_evolution(
    #     ej=dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs='rs',
    #     title="FS comoving spectra evolution", figname="spec_dyn_fsrs__rad_rs__num__ssa__ssc", show=plot,
    #     figsize=(5,6),
    #     tasks=dict(
    #         n_ele=dict(norm_method="/integ *y",#'*y /integ',
    #                    ylabel=r"$\gamma_{e}$",
    #                    zlabel=r"$(\gamma_{e} N_{e}) / N_{e;\, \rm tot}$",
    #                    norm="LogNorm",mode="contour",#"contour",
    #                    plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False,
    #                    vmin=1e-6,vmax=1e0,ylim=(1,1e4)),
    #         # yg=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$Y_g$ [cgs]",
    #         #         norm="LogNorm",
    #         #         plot_gm=True,plot_gc=True,plot_gM=True,plot_max=False),
    #         syn_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm syn}$ [cgs]",
    #                    norm="LogNorm",ylim=(1e7,1e25),vmin=1e-12,vmax=1e-6,mode="contour",
    #                    plot_num=True,plot_nuc=True,plot_nuM=True,plot_nua=False,plot_tau1=False,plot_max=True),
    #         # ssc_j=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' j'_{\rm ssc}$ [cgs]",
    #         #            norm="LogNorm", ylim=(1e16,1e28),mode=None,#,vmin=1e-9,vmax=1e-3 "contour",
    #         #            plot_num=True,plot_nuc=True,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=True),
    #         # syn_a=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\alpha'_{\rm syn}$ [cgs]",
    #         #          norm="LogNorm",mode="contour",
    #         #          plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
    #         #          ylim=(1e6,1e12)),
    #         syn_tau=dict(norm_method=None,ylabel= r"$\nu'$ [Hz]",zlabel=r"$\tau'_{\rm ssa}$",
    #                      norm="SymLogNorm",mode=None,
    #                      plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=True,plot_tau1=True,plot_max=False,
    #                      ylim=(1e7,1e14), xlim=(3e3,3e7)),
    #         total_i=dict(norm_method="*y",ylabel= r"$\nu'$ [Hz]",zlabel=r"$\nu' I'_{\rm total}$ [cgs]",
    #                      norm="LogNorm",vmin=1e-12,vmax=1e-6,mode="contour",
    #                      plot_num=False,plot_nuc=False,plot_nuM=False,plot_nua=False,plot_tau1=False,plot_max=False,
    #                      ylim=(1e7,1e25), xlim=(3e3,3e7)),# ,vmin=1e-15,vmax=1e-3 vmin=1e-17,vmax=1e-3
    #     ))
    plot_spectra_evolution(
        ej=dyn_fs__rad_rs__num__ssa__ssc.GRB, fs_or_rs='rs',
        title="FS comoving spectra evolution", figname="spec_dyn_fsrs__fluxdens_rs__num__ssa__ssc", show=plot,
        figsize=(6, 4),
        tasks=dict(fluxdens=dict(norm_method="*y", ylabel=r"$\nu'$ [Hz]", zlabel=r"$\nu' F'_{\rm total}$ [cgs]",
                                 norm="LogNorm", vmin=1e10, vmax=1e21, mode="contour",
                                 plot_num=False, plot_nuc=False, plot_nuM=False, plot_nua=False, plot_tau1=False,
                                 plot_max=False,
                                 ylim=(1e7, 1e28), xlim=(3e3, 1e9)),  # ,vmin=1e-15,vmax=1e-3 vmin=1e-17,vmax=1e-3
                   ))

    plot_emission_region_prop_3(ej=dyn_fs__rad_rs__num__ssa__ssc.GRB, fs_or_rs="fs",
                                xlim=(3e3, 1e9), task=dict(layers=(0,),
                                                           colors=("black",),
                                                           ylim0=(0.1, 10), ylim1=(0.1, 10), ylim2=(0.1, 10),
                                                           figname="tophat_bw_shock_props", show=True
                                                           ))

    return dyn_fs__rad_rs__num__ssa__ssc


def plot_total_observed_spectrum_2(ej_fs: PBA.Ejecta, ej_fsrs: PBA.Ejecta, norm_method: str, xlim: tuple, task: dict):
    task_i = task
    sp = PlotSpectra()
    fig, axes = plt.subplots(ncols=1, nrows=len(task_i["layers"]) + 1,
                             figsize=(5, 9), layout='constrained', sharex='all', sharey='all')

    def plot_one(ax, text, ilayer, ishell, sum_shells_layers, fs_rs_tot):

        ys, times, spec_fs = sp._get_spectrum(
            ej=ej_fs, v_n=task_i["v_n"], fs_or_rs="rs",
            ishell=ishell, ilayer=ilayer, sum_shells_layers=sum_shells_layers, norm_method=norm_method)
        ys, times, spec_fsrs = sp._get_spectrum(
            ej=ej_fsrs, v_n=task_i["v_n"], fs_or_rs="rs",
            ishell=ishell, ilayer=ilayer, sum_shells_layers=sum_shells_layers, norm_method=norm_method)

        spec_rs = spec_fsrs - spec_fs
        if fs_rs_tot == "total":
            spec = spec_fsrs
        elif fs_rs_tot == "fs":
            spec = spec_fs
        elif fs_rs_tot == "rs":
            spec = spec_rs

        norm = LogNorm(vmin=task_i.get("vmin", spec_fsrs.max() * 1e-20),
                       vmax=task_i.get("vmax", spec_fsrs.max() * 1))
        cmap = plt.get_cmap(task_i.get("cmap", 'jet'))
        cmap.set_under(task_i.get("set_under", 'white'), alpha=0.2)
        cmap.set_over(task_i.get("set_over", 'white'), alpha=0.2)

        if task_i["mode"] == "contour":
            spec = np.ma.masked_where(spec < norm.vmin, spec)
            spec = np.ma.masked_where(spec_fsrs > norm.vmax, spec)
            # _c = ax.pcolormesh(xs, ys, spec.T, cmap=cmap, norm=norm)

        if task_i["mode"] == "contour":
            # _c = ax.contourf(xs, ys, spec.T, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
            _c = ax.contourf(times, ys, spec, cmap=cmap, locator=ticker.LogLocator(), norm=norm)
        else:
            _c = ax.pcolormesh(times, ys, spec, cmap=cmap, norm=norm)
        bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
        ax.text(0.95, 0.07, text, fontsize=12, bbox=bbox,
                transform=ax.transAxes, horizontalalignment='right')

        # use contourf() with proper hatch pattern and alpha value
        val1 = np.array((spec_rs > spec_fs), dtype=np.float64)
        val1[val1 < 1] = np.nan
        cs = ax.contourf(times, ys, val1, hatches=['..'], alpha=0.0, color="green")  #,zorder=0

        return _c

    for i, ilayer in enumerate(task_i["layers"][::-1]):
        _c = plot_one(axes[i], f"il={ilayer}", ilayer, 0, False, task_i["fs_rs_total"])
    _c = plot_one(axes[-1], f"total", None, None, True, task_i["fs_rs_total"])

    cbar = fig.colorbar(_c, ax=axes, shrink=0.95, pad=0.01, aspect=40, extend='both')  # orientation = 'horizontal')
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label(task_i["zlabel"], size=12)

    #
    # val1 = np.array((gam_dot_adi>gam_dot_syn) * (gam_dot_adi>gam_dot_ssc), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=["\\\\"],  alpha=0.0) #,zorder=0
    #
    # val1 = np.array((gam_dot_ssc>gam_dot_syn) * (gam_dot_ssc>gam_dot_adi), dtype=np.float64)
    # val1[val1 < 1] = np.nan
    # cs = ax.contourf(xs, ys, val1.T  , hatches=['++'],  alpha=0.0, edgecolor = "r",facecolor="blue", color = "green") #,zorder=0
    for ax in axes:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
        ax.legend(fancybox=True, loc='upper right', columnspacing=0.8,
                  # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  shadow=False, ncol=4, fontsize=12, labelcolor='black',
                  framealpha=0.8, borderaxespad=0.)
        ax.set_rasterized(True)
        ax.set_ylabel(task_i["ylabel"], fontsize=12)
        # ax.set_title(task["title"], fontsize=12)
        ax.set_xlim(*xlim)
    axes[-1].set_xlabel(task_i["xlabel"], fontsize=12)

    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + task["figname"] + '.pdf')
    if task["show"]: plt.show()
    if plot: plt.show()
    plt.close(fig)


def tasks_gauss(do_run: bool, plot: bool, struct: dict, P: dict) -> PBA.PyBlastAfterglow:
    dyn_fs__rad_rs__num__ssa__ssc = run(
        working_dir=os.getcwd() + "/working_dirs/gauss_dyn_fs__rad_rs__num__ssa__ssc/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )
    #
    # plot_emission_region_prop_3(ej=dyn_fs__rad_rs__num__ssa__ssc.GRB,fs_or_rs="fs",
    #                             xlim=(3e3,1e9),task=dict(layers=(0,5,10,15,19),
    #                                                      colors=("red","orange","yellow","green","blue"),
    #                                                      ylim0=(0.1,10),ylim1=(0.1,10),ylim2=(0.1,10),
    #                                                      figname="tophat_bw_shock_props",show=plot
    #                                                      ))

    dyn_fs__rad_fs__num__ssa__ssc = run(
        working_dir=os.getcwd() + "/working_dirs/gauss_dyn_fs__rad_fs__num__ssa__ssc/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes',
                                                          do_rs_radiation="no",
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )

    plot_total_observed_spectrum_2(
        ej_fs=dyn_fs__rad_fs__num__ssa__ssc.GRB,
        ej_fsrs=dyn_fs__rad_rs__num__ssa__ssc.GRB, norm_method="*y",
        xlim=(3e3, 3e7), task=dict(show=True, figname="gauss_fluxdens_plot",
                                   v_n="fluxdens", fs_rs_total="total",
                                   norm="LogNorm", vmin=1e10, vmax=1e19, mode="contour",
                                   ylabel=r"$\nu$ [Hz]", xlabel=r"$t_{\rm obs}$ [s]",
                                   zlabel=r"$\nu F_{\rm total}$ [mJy]", title="Gauss jet; FS \& RS",
                                   set_under='blue', set_over='red', layers=(0, 5, 10, 19),
                                   plot_gm=False, plot_gc=False, plot_gM=False))


def tasks_fs_comparison(do_run: bool, plot: bool, struct: dict, P: dict):
    # --- fs -- fs ---
    dyn_fs__rad_fs__num = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_fs__num/",
        struct=struct, P=P,
        type="a", run=do_run
    )
    # dyn_fs__rad_fs__ana__ssa = run(
    #     working_dir=os.getcwd()+"/working_dirs/dyn_fs__rad_fs__ana_ssa/",
    #     struct=struct, P = d2d(default=P, new=dict(grb=dict(method_ele_fs='analytic',use_ssa_fs='yes',method_synchrotron_fs="Dermer09"))),#dyn_fs__rad_rs__num__ssa
    #     type="a",run=do_run
    # )
    dyn_fs__rad_fs__num__ssa = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_fs__num_ssa/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(use_ssa_fs='yes'))),  #method_synchrotron_fs="Dermer09"
        type="a", run=do_run
    )
    dyn_fs__rad_fs__num__noadi = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_fs__num__noadi/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(num_ele_use_adi_loss_fs='no'))),
        type="a", run=do_run
    )
    dyn_fs__rad_fs__num__ssc = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_fs__num__ssc/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(method_ssc_fs='numeric'))),
        type="a", run=do_run
    )
    dyn_fs__rad_fs__mix = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fs__rad_fs__mix/",
        struct=struct,
        P=d2d(default=P, new=dict(grb=dict(method_ele_fs='mix'))),
        type="a", run=do_run
    )

    for (model1, model2, name, label) in [
        (dyn_fs__rad_fs__num.GRB, dyn_fs__rad_fs__mix.GRB, "num_mix", "Numeric/Mixed"),
        (dyn_fs__rad_fs__num.GRB, dyn_fs__rad_fs__num__noadi.GRB, "adi_noadi", "Adi/noAdi"),
        (dyn_fs__rad_fs__num__ssa.GRB, dyn_fs__rad_fs__num.GRB, "ssa_nossa", "SSA/noSSA"),
        (dyn_fs__rad_fs__num__ssc.GRB, dyn_fs__rad_fs__num.GRB, "ssc_nossc", "SSC/noSSC"),
    ]:
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="n_ele", fs_or_rs="fs", task=dict(
                title="FS comoving electron spectra evolution ratio",
                figname=f"spec_fs_ele_ratio_{name}", show=plot,
                n_ele=dict(norm_method=None, ylabel=r"$\gamma_{e}$",
                           zlabel=label, vmin=1e-2, vmax=1e2, cmap="RdBu_r", xlim=(3e3, 8e6), set_under='blue',
                           set_over='red', mode=None,
                           norm="SymLogNorm", plot_gm=True, plot_gc=True, plot_gM=True, plot_max=False)))
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="syn_j", fs_or_rs="fs", task=dict(
                title="FS comoving synchrotron spectra evolution ratio",
                figname=f"spec_fs_synch_ratio_{name}", show=plot,
                syn_j=dict(norm_method=None, ylabel=r"$\nu'$",
                           zlabel=label, vmin=1e-2, vmax=1e2, cmap="RdBu_r", xlim=(3e3, 8e6),
                           norm="SymLogNorm", set_under='blue', set_over='red', mode=None,
                           plot_num=True, plot_nuc=True, plot_nuM=True, plot_nua=False, plot_tau1=False,
                           plot_max=True)))
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="syn_i", fs_or_rs="fs", task=dict(
                title="FS comoving intensity spectra evolution ratio",
                figname=f"spec_fs_int_ratio_{name}", show=plot,
                syn_i=dict(norm_method=None, ylabel=r"$\nu'$",
                           zlabel=label, vmin=1e-2, vmax=1e2, cmap="RdBu_r", xlim=(3e3, 8e6),  #ylim=(1e6,1e14),
                           norm="SymLogNorm", set_under='blue', set_over='red', mode=None,
                           plot_num=True, plot_nuc=True, plot_nuM=True, plot_nua=False, plot_tau1=False,
                           plot_max=False)))


def tasks_rs_comparison(do_run: bool, plot: bool, struct: dict, P: dict):
    # --- fsrs -- rs ---
    dyn_fsrs__rad_rs__num = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs'))),
        type="a", run=do_run
    )
    dyn_fsrs__rad_rs__ana__ssa = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__ana__ssa/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          method_synchrotron_fs="Dermer09",
                                                          method_synchrotron_rs="Dermer09",
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes'))),
        type="a", run=do_run
    )
    dyn_fsrs__rad_rs__num__ssa = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num__ssa/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          use_ssa_fs='yes',
                                                          use_ssa_rs='yes'))),
        type="a", run=do_run
    )
    dyn_fsrs__rad_rs__num__noadi = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num__noadi/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          num_ele_use_adi_loss_fs='no',
                                                          num_ele_use_adi_loss_rs='no'))),
        type="a", run=do_run
    )
    dyn_fsrs__rad_rs__num__ssc = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__num__ssc/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          method_ssc_fs='numeric',
                                                          method_ssc_rs='numeric'))),
        type="a", run=do_run
    )
    dyn_fsrs__rad_rs__mix = run(
        working_dir=os.getcwd() + "/working_dirs/dyn_fsrs__rad_rs__mix/",
        struct=struct, P=d2d(default=P, new=dict(grb=dict(do_rs='yes', bw_type='fsrs',
                                                          method_ele_fs='mix',
                                                          method_ele_rs='mix'))),
        type="a", run=do_run
    )

    for (model1, model2, name, label) in [
        (dyn_fsrs__rad_rs__num.GRB, dyn_fsrs__rad_rs__mix.GRB, "num_mix", "Numeric/Mixed"),
        (dyn_fsrs__rad_rs__num.GRB, dyn_fsrs__rad_rs__num__noadi.GRB, "adi_noadi", "Adi/noAdi"),
        (dyn_fsrs__rad_rs__num__ssa.GRB, dyn_fsrs__rad_rs__num.GRB, "ssa_nossa", "SSA/noSSA"),
        (dyn_fsrs__rad_rs__num__ssc.GRB, dyn_fsrs__rad_rs__num.GRB, "ssc_nossc", "SSC/noSSC"),
    ]:
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="n_ele", fs_or_rs="rs", task=dict(
                title="RS comoving electron spectra evolution ratio",
                figname=f"spec_rs_ele_ratio_{name}", show=plot,
                n_ele=dict(norm_method=None, ylabel=r"$\gamma_{e}$",
                           zlabel=label, vmin=1e-1, vmax=1e1, cmap="RdBu_r", xlim=(3e3, 8e6), set_under='blue',
                           set_over='red',
                           norm="SymLogNorm", plot_gm=True, plot_gc=True, plot_gM=True, plot_max=False)))
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="synch", fs_or_rs="rs", task=dict(
                title="RS comoving synchrotron spectra evolution ratio",
                figname=f"spec_rs_synch_ratio_{name}", show=plot,
                synch=dict(norm_method=None, ylabel=r"$\nu'$",
                           zlabel=label, vmin=1e-1, vmax=1e1, cmap="RdBu_r", xlim=(3e3, 8e6), set_under='blue',
                           set_over='red',
                           norm="SymLogNorm",
                           plot_num=True, plot_nuc=True, plot_nuM=True, plot_nua=False, plot_tau1=False,
                           plot_max=True)))
        plot_spectra_evolution_ratio(
            ej1=model1, ej2=model2,
            v_n="int", fs_or_rs="rs", task=dict(
                title="RS comoving intensity spectra evolution ratio",
                figname=f"spec_rs_int_ratio_{name}", show=plot,
                int=dict(norm_method=None, ylabel=r"$\nu'$",
                         zlabel=label, vmin=1e-1, vmax=1e1, cmap="RdBu_r", xlim=(3e3, 8e6), ylim=(1e6, 1e14),
                         norm="SymLogNorm", set_under='blue', set_over='red',
                         plot_num=True, plot_nuc=True, plot_nuM=True, plot_nua=False, plot_tau1=False, plot_max=True)))


if __name__ == '__main__':
    do_run = False
    plot = True
    struct = dict(struct="tophat", Eiso_c=1.e53, Gamma0c=400., M0c=-1., theta_c=0.1, theta_w=0.1)
    # struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1)
    P = dict(
        main=dict(n_ism=1., tb0=3e3, ntb=1000, rtol=5e-7, theta_obs=0,
                  lc_freqs='array logspace 1e8 1e29 96',
                  lc_times='array logspace 3e3 1e10 128'),
        grb=dict(save_dynamics='yes', save_spec='yes', do_lc='yes',
                 # method_nonrel_dist_fs='none',
                 # method_nonrel_dist_rs='none',
                 eps_e_fs=0.1, eps_b_fs=0.1, p_fs=2.1,
                 gamma_max_fs=4e7, method_gamma_max_fs="useConst",
                 gamma_max_rs=4e7, method_gamma_max_rs="useConst",
                 max_substeps_fs=1000, max_substeps_rs=1000,
                 # ngam_fs=1001,gam1_rs=1,gam2_rs=1e4,ngam_rs=1001,
                 # eps_b_fs = 1e-7,
                 # method_gamma_max_fs='useConst',method_gamma_max_rs='useConst',
                 # method_synchrotron_fs="Bessel",
                 # method_synchrotron_rs="Bessel",
                 # method_ele_fs='mix',
                 # method_ele_rs="mix",
                 # num_ele_use_adi_loss_fs='no',
                 # num_ele_use_adi_loss_rs='no',
                 gam1_fs=1., gam2_fs=1e8, ngam_fs=401,
                 gam1_rs=1., gam2_rs=1e8, ngam_rs=401,
                 freq1_fs=1e6, freq2_fs=1e32, nfreq_fs=401,
                 freq1_rs=1e6, freq2_rs=1e32, nfreq_rs=401,
                 # ebl_tbl_fpath="none"
                 )
    )

    # --- fs -- fs ---
    pba_fs = tasks_fs(do_run=do_run, plot=plot, struct=struct, P=P)
    # tasks_fs_comparison(do_run=do_run, plot=plot, struct=struct, P=P)

    # --- fsrs -- rs ---
    # pba_fsrs = tasks_rs(do_run=do_run, plot=plot, struct=struct, P=P)
    # tasks_rs_comparison(do_run=do_run, plot=plot, struct=struct, P=P)

    # plot_total_observed_spectrum(
    #     do_run=do_run,norm_method="*y",
    #     xlim=(3e3,3e7), task=dict(show=True,figname="tophat_fluxdens_plot",
    #                               norm="LogNorm",vmin=1e10,vmax=1e21,mode="contour",
    #                               ylabel=r"$\nu$ [Hz]", xlabel=r"$t_{\rm obs}$ [s]",
    #                               zlabel=r"$\nu F_{\rm total}$ [mJy]",title="Tophat jet; FS \& RS",
    #                               set_under='blue',set_over='red',
    #                               plot_gm=False,plot_gc=False,plot_gM=False))

    # struct_gauss = dict(struct="gaussian",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.3)
    # pba_gauss = tasks_gauss(do_run=do_run, plot=plot, struct=struct_gauss, P=P)
