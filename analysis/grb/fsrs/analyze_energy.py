import package.src.PyBlastAfterglowMag as PBA
import os,shutil,matplotlib.pyplot as plt
import numpy as np
import copy
working_dir = os.getcwd()+'/tmp1/'
fig_dir = os.getcwd()+'/figs/'

class PlotSpectra:
    qe = 4.803204e-10
    me = 9.1094e-28
    c = 2.99792458e10

    def __init__(self):
        pass



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
    def nuprime_to_nu(nuprime: np.ndarray, ej: PBA.Ejecta):
        z = 0
        Gamma = ej.get_dyn_arr(v_n="Gamma", ishell=0, ilayer=0)
        val = nuprime * Gamma / (1 + z)
        return val



def d2d(default: dict, new: dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_

def gamma_adi(Gamma, beta):
    """ Adiabatic index of the fluid From Nava 2013 paper """
    return (4. + 1. / Gamma) / 3.

def GammaEff(Gamma, gammaAdi):
    return (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

def run(working_dir:str, struct:dict, P:dict, type:str="a", do_run:bool=True) -> PBA.PyBlastAfterglow:
    # clean he temporary direcotry
    if os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    os.mkdir(working_dir)

    # generate initial data for blast waves
    pba_id = PBA.id_analytic.JetStruct(n_layers_pw=81,
                                       n_layers_a=1 if struct["struct"]=="tophat" else 21)

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
    if do_run:
        pba.run(
            path_to_cpp_executable="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
            loglevel="info"
        )
        # process skymap
        if (pba.GRB.opts["do_skymap"]=="yes"):
            conf = {"nx":128, "ny":64, "extend_grid":1.1, "fwhm_fac":0.5, "lat_dist_method":"integ",
                    "intp_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }, # "gaussian"
                    "hist_filter":{ "type":'gaussian', "sigma":2, "mode":'reflect' }}
            prep = PBA.skymap_process.ProcessRawSkymap(conf=conf, verbose=False)
            prep.process_singles(infpaths=working_dir+"raw_skymap_*.h5",
                                 outfpath=pba.GRB.fpath_sky_map,
                                 remove_input=True)

    return pba

def plot_fs_energy(struct:dict,pp:dict,plot:dict):
    pba = run(working_dir=working_dir,struct=struct,P=pp,type='a')
    #
    fig, axes = plt.subplots(figsize=(5.5,5.), ncols=1, nrows=2, sharex="all")
    if not hasattr(axes,'__len__'): axes = [axes]


    E0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["E0"])
    M0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["M0"])
    Rd = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["Rd"])
    Gamma0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["Gamma0"])
    theta_b0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["theta_b0"])
    print(f"E0 = {E0}, M0 = {M0}, Rd = {Rd} n_ism={pba.main_pars['n_ism']}, "
          f"Gamma0 = {Gamma0}, theta0 = {theta_b0}")

    R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0)
    Rsh = pba.GRB.get_dyn_arr(v_n="Rsh",ishell=0,ilayer=0)
    theta = pba.GRB.get_dyn_arr(v_n="theta",ishell=0,ilayer=0)
    Gamma = pba.GRB.get_dyn_arr(v_n="Gamma",ishell=0,ilayer=0)
    mom = pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0)
    thickness = pba.GRB.get_dyn_arr(v_n="thickness",ishell=0,ilayer=0)
    Gamma_fs = pba.GRB.get_dyn_arr(v_n="GammaFsh",ishell=0,ilayer=0)
    beta = pba.GRB.get_dyn_arr(v_n="beta",ishell=0,ilayer=0)
    Eint2 = pba.GRB.get_dyn_arr(v_n="Eint2",ishell=0,ilayer=0)
    M2 = pba.GRB.get_dyn_arr(v_n="M2",ishell=0,ilayer=0)
    tburst = pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=0)
    gamAdi = gamma_adi(Gamma, beta)
    GammaEff2 = GammaEff(Gamma, gamAdi)

    # energy
    Ekin4 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * np.full_like(R,M0)
    Ekin2 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * M2
    Eint2 *= GammaEff2

    # plot
    axes[0].plot(R, Eint2/E0, label=r"$E_{\rm int; 2}$", color='blue', ls='--')
    axes[0].plot(R, Ekin2/E0, label=r"$E_{\rm kin; 2}$", color='blue', ls=':')

    if (pba.GRB.opts["do_rs"] == "yes"):
        Gamma43 = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=0)
        Eint3 = pba.GRB.get_dyn_arr(v_n="Eint3",ishell=0,ilayer=0)
        M3 = pba.GRB.get_dyn_arr(v_n="M3",ishell=0,ilayer=0)
        Gamma_rs = pba.GRB.get_dyn_arr(v_n="GammaRsh",ishell=0,ilayer=0)
        thickness_rs = pba.GRB.get_dyn_arr(v_n="thichness_rs",ishell=0,ilayer=0)
        delraR4 = pba.GRB.get_dyn_arr(v_n="deltaR4",ishell=0,ilayer=0)
        Rrs = pba.GRB.get_dyn_arr(v_n="Rsh_rs",ishell=0,ilayer=0)

        gamAdi3 = gamma_adi(Gamma43, PBA.utils.get_beta(Gamma43))
        GammaEff3 = GammaEff(Gamma43, gamAdi3)

        GammaEff3[~np.isfinite(GammaEff3)] = 0.
        Eint3 *= GammaEff3

        Ekin4 = (Gamma0 - 1.0) * (M0 - M3) * PBA.utils.cgs.c ** 2
        Ekin3 = (Gamma - 1.0) * M3 * PBA.utils.cgs.c ** 2 # TODO Is it Gamma43 or Gamma here????

        # axes[0].plot(tburst, Etot3/E0, label=r"$E_{\rm tot; 3}$", color='green')
        axes[0].plot(R, Eint3/E0, label=r"$E_{\rm int; 3}$", color='green', ls="--")
        axes[0].plot(R[Gamma43>0], (Ekin3/E0)[Gamma43>0], label=r"$E_{\rm kin; 3}$", color='green', ls=":")
        axes[0].plot(R, Ekin4/E0, label=r"$E_{\rm kin; 4}$", color='green', ls="-.")

    if (pba.GRB.opts["do_rs"] == "yes"):
        axes[0].plot(R, (Ekin2 + Eint2 + Ekin4 + Ekin3 + Eint3)/E0, ls='-', label=r"$E_{\rm tot}$", color='black')
    else:
        axes[0].plot(R, Ekin4/E0, label=r"$E_{\rm kin; 3,4}$", color='blue', ls='-.')
        axes[0].plot(R, (Ekin2 + Eint2 + Ekin4)/E0, ls='-', label=r"$E_{\rm tot}$", color='black')
    # axes[0].set_ylim(*plot["ylim"])
    axes[0].set_ylabel(r"Energy/$E_0$", fontsize=12)
    # axes[0].axhline(y=1, color='black', linestyle='-', linewidth=0.5)

    axes[0].text(0.97, 0.54, plot["text"], fontsize=12, # bbox=bbox,
                 transform=axes[0].transAxes, horizontalalignment='right')

    # axes[1].axvline(x=tburst[PBA.utils.find_nearest_index(R,Rd)],
    #                 color='gray',linestyle=':',label='$t_{\rm dec}$')
    # tmp = 2/((M2+M0)**2*(Gamma[0]+1)/(M0**2*(Gamma0-1))-1)+1
    def get_Rdec2(E, nn, Gamma):
        rdec = (3. / (4. * np.pi) * 1. / (PBA.utils.cgs.c ** 2. * PBA.utils.cgs.mp) * E / (nn * Gamma ** 2.)) ** (1. / 3.)
        return rdec
    if plot["rdec"]:
        rdec = get_Rdec2(struct["Eiso_c"],nn=pba.main_pars["n_ism"],Gamma=struct["Gamma0c"])
        # tdec = tburst[PBA.utils.find_nearest_index(R,rdec)],
        axes[1].axvline(x=rdec,color='gray',linestyle='dotted',label=r"$R_{\rm dec}$")

    def get_bm79(E0, A0, s, R):
        Gamma = np.sqrt((17. - 4.*s) * E0 / (16 * np.pi * A0 * (PBA.utils.cgs.c ** 2) * R ** (3 - s))) # FIXED
        beta = (1 - 1. / Gamma ** 2) ** (1 / 2.)
        # return R[Gamma * beta > 1.], Gamma[Gamma * beta > 1.], beta[Gamma * beta > 1.]

        return Gamma*beta
    if plot["bm"]:
        Gamma_bm = get_bm79(E0=struct["Eiso_c"], A0=pba.main_pars["n_ism"]*PBA.utils.cgs.mp, s=0, R=R)
        axes[1].plot(R[mom>1.],Gamma_bm[mom>1.],color='gray',ls='-.',label='BM76')

    def get_sedovteylor(Rdec, beta0, R):
        # Rdec = float(((3 - s) * M0 / (4 * np.pi * A0 * Gamma0)) ** (1 / (3. - s)))

        beta = np.zeros(len(R))
        beta.fill(beta0)
        beta[R > Rdec] = beta0 * (R[R > Rdec] / Rdec) ** (-3 / 2)
        Gamma_st = 1 / (1 - beta ** 2) ** (1 / 2)

        # return R[R > Rdec] + Rdec, Gamma_st[R > Rdec], beta[R > Rdec]
        return Gamma_st*beta
    # beta_st = get_sedovteylor(5e16,#R[PBA.utils.find_nearest_index(Gamma,1.1)],
    #                           beta0=PBA.get_beta(Gamma[0]), R=R)
    # axes[1].plot(R,beta_st,color='gray',ls='--',label='ST')

    if plot["theta_spread_0"]:
        itheta = np.argmax(theta!=theta[0])
        axes[1].axvline(x=R[itheta],color='black',linestyle='dashdot',label=r"$\omega\neq\omega_0$")
    if plot["theta_spread_1"]:
        itheta = np.argmax(np.where(theta<0.99*np.pi/2.))
        axes[1].axvline(x=R[itheta],color='black',linestyle='dotted',label=r"$\omega=\pi/2$")

    # axes[1].plot(tburst, Gamma, color='blue', linestyle='-', label=r"$\Gamma$")
    # axes[1].plot(tburst, Gamma_fs, color='blue', linestyle='--', label=r"$\Gamma_{\rm fs}$")
    # axes[1].plot(tburst, thickness, color='blue', linestyle=':', label=r"$\Delta_{\rm fs}$")
    # axes[1].plot(tburst, R, color='gray', linestyle='-', label=r"$R$")
    # axes[1].plot(tburst, Rsh, color='blue', linestyle='--', label=r"$R_{fs}$", lw=2)
    axes[1].plot(R, mom, color='blue', linestyle='-', label=r"$\Gamma\beta$")

    if (pba.GRB.opts["do_rs"] == "yes"):
        # axes[1].plot(tburst, Gamma_rs, color='green', linestyle='--', label=r"$\Gamma_{\rm rs}$")
        # axes[1].plot(tburst, thickness_rs, color='green', linestyle='-.', label=r"$\Delta_{\rm rs}$")
        # axes[1].plot(tburst, delraR4, color='green', linestyle=':', label=r"$\Delta_{4}$")
        # axes[1].plot(tburst, Rrs, color='green', linestyle='--', label=r"$R_{\rm rs}$")
        axes[1].plot(R, Gamma43*PBA.utils.get_beta(Gamma43), color='green', linestyle='-', label=r"$\Gamma_{43}\beta_{43}$")


    axes[1].set_ylabel(r"$\Gamma\beta$",fontsize=12)

    for i, ax in enumerate(axes):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        ax.set_ylim(*plot[f"ylim{i+1}"])
        ax.set_xlim(*plot["xlim"])
        ax.grid(ls=':')
    axes[0].legend(fancybox=True, loc='best', columnspacing=1.2,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 1,
              fontsize=12,
              framealpha=0., borderaxespad=0.)
    axes[1].legend(fancybox=True, loc='best',columnspacing=1.2,
                    # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                    shadow=False, ncol= 3 if pba.GRB.opts["do_rs"]=="yes" else 1,
                    fontsize=12,
                    framealpha=0., borderaxespad=0.)

    # ax2 =axes[1].twinx()
    # # ax2.plot(R,pba.GRB.get_dyn_arr(v_n="rho4",ishell=0,ilayer=0))
    # ax2.plot(R,pba.GRB.get_dyn_arr(v_n="rho4",ishell=0,ilayer=0))
    # ax2.set_yscale('log')

    axes[-1].set_xlabel(r"$R$ [cm]", fontsize=12)
    plt.tight_layout()
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()

def plot_fs_energy2(struct:dict,pp:dict,plot:dict):
    pba = run(working_dir=working_dir,struct=struct,P=pp,type='a')
    #
    fig, axes = plt.subplots(figsize=(5.5,6.), ncols=1, nrows=3, sharex="all",layout='constrained',
                             gridspec_kw=dict(height_ratios=[1,2,2]))
    if not hasattr(axes,'__len__'): axes = [axes]

    def get_total_energy(ishell,ilayer):
        key = f"shell={ishell} layer={ilayer}"
        Gamma = pba.GRB.get_dyn_arr(v_n="Gamma",ishell=ishell,ilayer=ilayer)
        beta = pba.GRB.get_dyn_arr(v_n="beta",ishell=ishell,ilayer=ilayer)
        M0 = float(pba.GRB.get_dyn_obj()[key].attrs["M0"])
        E0 = float(pba.GRB.get_dyn_obj()[key].attrs["E0"])
        Gamma0 = float(pba.GRB.get_dyn_obj()[key].attrs["Gamma0"])
        Eint2 = pba.GRB.get_dyn_arr(v_n="Eint2",ishell=ishell,ilayer=ilayer)
        M2 = pba.GRB.get_dyn_arr(v_n="M2",ishell=ishell,ilayer=ilayer)
        gamAdi = gamma_adi(Gamma, beta)
        GammaEff2 = GammaEff(Gamma, gamAdi)

        Ekin4 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * np.full_like(Gamma,M0)
        Ekin2 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * M2
        Eint2 *= GammaEff2

        Etot = Ekin2 + Eint2 + Ekin4

        if (pba.GRB.opts["do_rs"] == "yes"):
            M3 = pba.GRB.get_dyn_arr(v_n="M3",ishell=ishell,ilayer=ilayer)
            Eint3 = pba.GRB.get_dyn_arr(v_n="Eint3",ishell=0,ilayer=0)
            Gamma43 = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=ishell,ilayer=ilayer)

            gamAdi3 = gamma_adi(Gamma43, PBA.utils.get_beta(Gamma43))
            GammaEff3 = GammaEff(Gamma43, gamAdi3)

            GammaEff3[~np.isfinite(GammaEff3)] = 0.
            Eint3 *= GammaEff3

            Ekin4 = (Gamma0 - 1.0) * (M0 - M3) * PBA.utils.cgs.c ** 2
            # Ekin4 = (Gamma0 - 1.0) * (M0 - M3) * PBA.utils.cgs.c ** 2
            Ekin3 = (Gamma - 1.0) * M3 * PBA.utils.cgs.c ** 2 # TODO is it Gamma or Gamma43 here????

            Etot = Ekin2 + Eint2 + Ekin4 + Ekin3 + Eint3

        return Etot/E0

    # plot total energy
    ax = axes[0]
    ax.set_yscale("log")
    for i, ilayer in enumerate(plot["layers"]):
        R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=ilayer)
        ax.plot(R,get_total_energy(0,ilayer),color=plot["colors"][i],ls='-')
    ax.set_ylabel(r"$E_{\rm tot}/E_{0}$", fontsize=12)

    # plot momenta
    ax = axes[1]
    ax.set_yscale("log")
    for i, ilayer in enumerate(plot["layers"]):
        R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=ilayer)
        ax.plot(R,pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=ilayer),
                color=plot["colors"][i],ls='-')
        if (pba.GRB.opts["do_rs"] == "yes"):
            ax.plot(R,pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=ilayer) * \
                    PBA.utils.get_beta(pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=ilayer)),
                    color=plot["colors"][i],ls='--')
    ax.set_ylabel(r"Momentum", fontsize=12)
    ax.plot([0,0],[1,1],color='gray',ls='-',label=r"$\Gamma\beta$")
    ax.plot([0,0],[1,1],color='gray',ls='--',label=r"$\Gamma_{43}\beta_{43}$")
    ax.legend(fancybox=True, loc='upper right', columnspacing=1.2,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol= 2,
                   fontsize=12,
                   framealpha=0., borderaxespad=0.)
    # # plot momenta
    # ax = axes[2]
    # ax.set_yscale("log")
    # for i, ilayer in enumerate(plot["layers"]):
    #     R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=ilayer)
    #     ax.plot(R,pba.GRB.get_dyn_arr(v_n="M2",ishell=0,ilayer=ilayer),
    #             color=plot["colors"][i],ls='-')
    #     if (pba.GRB.opts["do_rs"] == "yes"):
    #         ax.plot(R,pba.GRB.get_dyn_arr(v_n="M3",ishell=0,ilayer=ilayer),
    #                 color=plot["colors"][i],ls='--')
    # ax.set_ylabel(r"Momentum", fontsize=12)
    # ax.plot([0,0],[1,1],color='gray',ls='-',label=r"$\Gamma\beta$")
    # ax.plot([0,0],[1,1],color='gray',ls='--',label=r"$\Gamma_{43}\beta_{43}$")
    # ax.legend(fancybox=True, loc='upper right', columnspacing=1.2,
    #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #           shadow=False, ncol= 2,
    #           fontsize=12,
    #           framealpha=0., borderaxespad=0.)
    # plot spreading
    ax = axes[2]
    ax.set_yscale("linear")
    for i, ilayer in enumerate(plot["layers"]):
        R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=ilayer)
        print(pba.GRB.get_dyn_obj()[f"shell={0} layer={ilayer}"].attrs.keys())
        theta_c_l = float(pba.GRB.get_dyn_obj()[f"shell={0} layer={ilayer}"].attrs["theta_c_l"])
        ax.plot(R,theta_c_l+pba.GRB.get_dyn_arr(v_n="theta",ishell=0,ilayer=ilayer),
                color=plot["colors"][i],ls='-',label=f"layer {ilayer}")
    ax.set_ylabel(r"$\omega$ [rad]", fontsize=12)
    ax.legend(fancybox=True, loc='upper left', columnspacing=1.2,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 3,
              fontsize=12,
              framealpha=0., borderaxespad=0.)

    for i, ax in enumerate(axes):
        ax.set_xscale("log")
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        if (f"ylim{i+1}" in plot.keys()): ax.set_ylim(*plot[f"ylim{i+1}"])
        ax.set_xlim(*plot["xlim"])
        ax.grid(ls=':')

    axes[-1].set_xlabel(r"$R$ [cm]", fontsize=12)
    plt.tight_layout()
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()
    plt.close(fig)




    E0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["E0"])
    M0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["M0"])
    Rd = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["Rd"])
    Gamma0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["Gamma0"])
    theta_b0 = float(pba.GRB.get_dyn_obj()["shell=0 layer=0"].attrs["theta_b0"])
    # print(f"E0 = {E0}, M0 = {M0}, Rd = {Rd} n_ism={pba.main_pars['n_ism']}, "
    #       f"Gamma0 = {Gamma0}, theta0 = {theta_b0}")

    R = pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0)
    Rsh = pba.GRB.get_dyn_arr(v_n="Rsh",ishell=0,ilayer=0)
    theta = pba.GRB.get_dyn_arr(v_n="theta",ishell=0,ilayer=0)
    Gamma = pba.GRB.get_dyn_arr(v_n="Gamma",ishell=0,ilayer=0)
    mom = pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0)
    thickness = pba.GRB.get_dyn_arr(v_n="thickness",ishell=0,ilayer=0)
    Gamma_fs = pba.GRB.get_dyn_arr(v_n="GammaFsh",ishell=0,ilayer=0)
    beta = pba.GRB.get_dyn_arr(v_n="beta",ishell=0,ilayer=0)
    Eint2 = pba.GRB.get_dyn_arr(v_n="Eint2",ishell=0,ilayer=0)
    M2 = pba.GRB.get_dyn_arr(v_n="M2",ishell=0,ilayer=0)
    tburst = pba.GRB.get_dyn_arr(v_n="tburst",ishell=0,ilayer=0)
    gamAdi = gamma_adi(Gamma, beta)
    GammaEff2 = GammaEff(Gamma, gamAdi)

    # energy
    Ekin4 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * np.full_like(R,M0)
    Ekin2 = (Gamma - 1) * PBA.utils.cgs.c ** 2 * M2
    Eint2 *= GammaEff2

    # plot
    axes[0].plot(R, Eint2/E0, label=r"$E_{\rm int; 2}$", color='blue', ls='--')
    axes[0].plot(R, Ekin2/E0, label=r"$E_{\rm kin; 2}$", color='blue', ls=':')

    if (pba.GRB.opts["do_rs"] == "yes"):
        Gamma43 = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=0)
        Eint3 = pba.GRB.get_dyn_arr(v_n="Eint3",ishell=0,ilayer=0)
        M3 = pba.GRB.get_dyn_arr(v_n="M3",ishell=0,ilayer=0)
        Gamma_rs = pba.GRB.get_dyn_arr(v_n="GammaRsh",ishell=0,ilayer=0)
        thickness_rs = pba.GRB.get_dyn_arr(v_n="thichness_rs",ishell=0,ilayer=0)
        delraR4 = pba.GRB.get_dyn_arr(v_n="deltaR4",ishell=0,ilayer=0)
        Rrs = pba.GRB.get_dyn_arr(v_n="Rsh_rs",ishell=0,ilayer=0)

        gamAdi3 = gamma_adi(Gamma43, Gamma43)
        GammaEff3 = GammaEff(Gamma43, gamAdi3)

        GammaEff3[~np.isfinite(GammaEff3)] = 0.
        Eint3 *= GammaEff3

        Ekin4 = (Gamma0 - 1.0) * (M0 - M3) * PBA.utils.cgs.c ** 2
        Ekin3 = (Gamma43 - 1.0) * M3 * PBA.utils.cgs.c ** 2

        # axes[0].plot(tburst, Etot3/E0, label=r"$E_{\rm tot; 3}$", color='green')
        axes[0].plot(R, Eint3/E0, label=r"$E_{\rm int; 3}$", color='green', ls="--")
        axes[0].plot(R, Ekin3/E0, label=r"$E_{\rm kin; 3}$", color='green', ls=":")
        axes[0].plot(R, Ekin4/E0, label=r"$E_{\rm kin; 4}$", color='green', ls="-.")

    if (pba.GRB.opts["do_rs"] == "yes"):
        axes[0].plot(R, (Ekin2 + Eint2 + Ekin4 + Ekin3 + Eint3)/E0, ls='-', label=r"$E_{\rm tot}$", color='black')
    else:
        axes[0].plot(R, Ekin4/E0, label=r"$E_{\rm kin; 3,4}$", color='blue', ls='-.')
        axes[0].plot(R, (Ekin2 + Eint2 + Ekin4)/E0, ls='-', label=r"$E_{\rm tot}$", color='black')
    # axes[0].set_ylim(*plot["ylim"])
    axes[0].set_ylabel(r"Energy/$E_0$", fontsize=12)
    # axes[0].axhline(y=1, color='black', linestyle='-', linewidth=0.5)

    axes[0].text(0.97, 0.54, plot["text"], fontsize=12, # bbox=bbox,
                 transform=axes[0].transAxes, horizontalalignment='right')

    # axes[1].axvline(x=tburst[PBA.utils.find_nearest_index(R,Rd)],
    #                 color='gray',linestyle=':',label='$t_{\rm dec}$')
    # tmp = 2/((M2+M0)**2*(Gamma[0]+1)/(M0**2*(Gamma0-1))-1)+1
    def get_Rdec2(E, nn, Gamma):
        rdec = (3. / (4. * np.pi) * 1. / (PBA.utils.cgs.c ** 2. * PBA.utils.cgs.mp) * E / (nn * Gamma ** 2.)) ** (1. / 3.)
        return rdec
    if plot["rdec"]:
        rdec = get_Rdec2(struct["Eiso_c"],nn=pba.main_pars["n_ism"],Gamma=struct["Gamma0c"])
        # tdec = tburst[PBA.utils.find_nearest_index(R,rdec)],
        axes[1].axvline(x=rdec,color='gray',linestyle='dotted',label=r"$R_{\rm dec}$")

    def get_bm79(E0, A0, s, R):
        Gamma = np.sqrt((17. - 4.*s) * E0 / (16 * np.pi * A0 * (PBA.utils.cgs.c ** 2) * R ** (3 - s))) # FIXED
        beta = (1 - 1. / Gamma ** 2) ** (1 / 2.)
        # return R[Gamma * beta > 1.], Gamma[Gamma * beta > 1.], beta[Gamma * beta > 1.]

        return Gamma*beta
    if plot["bm"]:
        Gamma_bm = get_bm79(E0=struct["Eiso_c"], A0=pba.main_pars["n_ism"]*PBA.utils.cgs.mp, s=0, R=R)
        axes[1].plot(R[mom>1.],Gamma_bm[mom>1.],color='gray',ls='-.',label='BM76')

    def get_sedovteylor(Rdec, beta0, R):
        # Rdec = float(((3 - s) * M0 / (4 * np.pi * A0 * Gamma0)) ** (1 / (3. - s)))

        beta = np.zeros(len(R))
        beta.fill(beta0)
        beta[R > Rdec] = beta0 * (R[R > Rdec] / Rdec) ** (-3 / 2)
        Gamma_st = 1 / (1 - beta ** 2) ** (1 / 2)

        # return R[R > Rdec] + Rdec, Gamma_st[R > Rdec], beta[R > Rdec]
        return Gamma_st*beta
    # beta_st = get_sedovteylor(5e16,#R[PBA.utils.find_nearest_index(Gamma,1.1)],
    #                           beta0=PBA.get_beta(Gamma[0]), R=R)
    # axes[1].plot(R,beta_st,color='gray',ls='--',label='ST')

    if plot["theta_spread_0"]:
        itheta = np.argmax(theta!=theta[0])
        axes[1].axvline(x=R[itheta],color='black',linestyle='dashdot',label=r"$\omega\neq\omega_0$")
    if plot["theta_spread_1"]:
        itheta = np.argmax(np.where(theta<0.99*np.pi/2.))
        axes[1].axvline(x=R[itheta],color='black',linestyle='dotted',label=r"$\omega=\pi/2$")

    # axes[1].plot(tburst, Gamma, color='blue', linestyle='-', label=r"$\Gamma$")
    # axes[1].plot(tburst, Gamma_fs, color='blue', linestyle='--', label=r"$\Gamma_{\rm fs}$")
    # axes[1].plot(tburst, thickness, color='blue', linestyle=':', label=r"$\Delta_{\rm fs}$")
    # axes[1].plot(tburst, R, color='gray', linestyle='-', label=r"$R$")
    # axes[1].plot(tburst, Rsh, color='blue', linestyle='--', label=r"$R_{fs}$", lw=2)
    axes[1].plot(R, mom, color='blue', linestyle='-', label=r"$\Gamma\beta$")

    if (pba.GRB.opts["do_rs"] == "yes"):
        # axes[1].plot(tburst, Gamma_rs, color='green', linestyle='--', label=r"$\Gamma_{\rm rs}$")
        # axes[1].plot(tburst, thickness_rs, color='green', linestyle='-.', label=r"$\Delta_{\rm rs}$")
        # axes[1].plot(tburst, delraR4, color='green', linestyle=':', label=r"$\Delta_{4}$")
        # axes[1].plot(tburst, Rrs, color='green', linestyle='--', label=r"$R_{\rm rs}$")
        axes[1].plot(R, Gamma43*PBA.utils.get_beta(Gamma43), color='green', linestyle='-', label=r"$\Gamma_{43}\beta_{43}$")


    axes[1].set_ylabel(r"$\Gamma\beta$",fontsize=12)

    for i, ax in enumerate(axes):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        ax.set_ylim(*plot[f"ylim{i+1}"])
        ax.set_xlim(*plot["xlim"])
        ax.grid(ls=':')
    axes[0].legend(fancybox=True, loc='best', columnspacing=1.2,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol= 1,
                   fontsize=12,
                   framealpha=0., borderaxespad=0.)
    axes[1].legend(fancybox=True, loc='best',columnspacing=1.2,
                   # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol= 3 if pba.GRB.opts["do_rs"]=="yes" else 1,
                   fontsize=12,
                   framealpha=0., borderaxespad=0.)

    # ax2 =axes[1].twinx()
    # # ax2.plot(R,pba.GRB.get_dyn_arr(v_n="rho4",ishell=0,ilayer=0))
    # ax2.plot(R,pba.GRB.get_dyn_arr(v_n="rho4",ishell=0,ilayer=0))
    # ax2.set_yscale('log')

    axes[-1].set_xlabel(r"$R$ [cm]", fontsize=12)
    plt.tight_layout()
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()

def plot_fs_energy_rad(struct:dict,pp:dict,plot:dict,fs_or_rs="fs",ishell=0,ilayer=0):

    pba = run(working_dir=working_dir,struct=struct,
              P=d2d(default=pp, new=dict(grb=dict(epsilon_e_rad=-1,epsilon_e_rad_rs=-1))),
              type="a")


    p = float(pba.GRB.pars["p_fs" if fs_or_rs == "fs" else "p_rs"])
    eps_e = float(pba.GRB.pars["eps_e_fs" if fs_or_rs == "fs" else "eps_e_rs"])

    def i_eps_an():
        mass = pba.GRB.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
        dens = pba.GRB.get_dyn_arr(v_n="rho2" if fs_or_rs == "fs" else "rho3", ishell=ishell, ilayer=ilayer)
        vol = mass / dens

        Gamma = pba.GRB.get_dyn_arr(v_n="Gamma" if fs_or_rs == "fs" else "Gamma43", ishell=ishell, ilayer=ilayer)
        U_p = pba.GRB.get_dyn_arr(v_n="U_p" if fs_or_rs == "fs" else "U_p3", ishell=ishell, ilayer=ilayer)
        Esh = pba.GRB.get_dyn_arr(v_n="Esh2" if fs_or_rs == "fs" else "Esh3", ishell=ishell, ilayer=ilayer)
        dEsh = np.diff(Esh)
        dEsh = np.insert(dEsh, 0, 0)
        B = pba.GRB.get_dyn_arr(v_n="B" if fs_or_rs == "fs" else "B_rs", ishell=ishell, ilayer=ilayer)
        R = pba.GRB.get_dyn_arr(v_n="R" if fs_or_rs == "fs" else "R", ishell=ishell, ilayer=ilayer)
        tburst = pba.GRB.get_dyn_arr(v_n="tburst" if fs_or_rs == "fs" else "tburst", ishell=ishell, ilayer=ilayer)
        deltaR = np.diff(R)
        deltaR = np.insert(deltaR, 0, 0)
        beta = PBA.utils.get_beta(Gamma)
        one_over_beta = 1. / beta
        gamma_min = pba.GRB.get_dyn_arr(v_n="gamma_min" if fs_or_rs == "fs" else "gamma_min_rs", ishell=ishell,
                                   ilayer=ilayer)
        gamma_c = pba.GRB.get_dyn_arr(v_n="gamma_c" if fs_or_rs == "fs" else "gamma_c_rs", ishell=ishell, ilayer=ilayer)
        gamma_max = pba.GRB.get_dyn_arr(v_n="gamma_max" if fs_or_rs == "fs" else "gamma_max_rs", ishell=ishell,
                                   ilayer=ilayer)
        m = pba.GRB.get_dyn_arr(v_n="M2" if fs_or_rs == "fs" else "M3", ishell=ishell, ilayer=ilayer)
        tcomov = pba.GRB.get_dyn_arr(v_n="tcomov" if fs_or_rs == "fs" else "tcomov", ishell=ishell, ilayer=ilayer)
        dtcomov = np.diff(tcomov)
        dtcomov = np.insert(dtcomov, 0, 0)
        dm = np.diff(m)
        dm = np.insert(dm, 0, 0) / PBA.utils.cgs.mp
        p = float(pba.GRB.pars["p_fs" if fs_or_rs == "fs" else "p_rs"])
        eps_e = float(pba.GRB.pars["eps_e_fs" if fs_or_rs == "fs" else "eps_e_rs"])
        epsilon = np.zeros_like(B)
        dr = pba.GRB.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thickness_rs", ishell=ishell, ilayer=ilayer)
        GammaSh = pba.GRB.get_dyn_arr(v_n="GammaFsh" if fs_or_rs == "fs" else "GammaRsh", ishell=ishell, ilayer=ilayer)

        # nu_arr = np.logspace(6,30,601)
        nu_arr=[1e6]
        P_total = np.zeros_like(Gamma)
        eps = np.zeros_like(Gamma)

        drcomov = dr * GammaSh
        dt = drcomov / PBA.utils.cgs.c

        for i in range(len(Gamma)-1):
            nu_m = PlotSpectra._nu_m(gm=gamma_min[i],gc=gamma_c[i],B=B[i],p=p,ssc=False)
            nu_c = PlotSpectra._nu_c(gm=gamma_min[i],gc=gamma_c[i],B=B[i],p=p,ssc=False)
            nu_M = PlotSpectra._nu_M(gM=gamma_max[i],B=B[i],p=p,ssc=False)
            pmax = PlotSpectra._pmax(nu_m=nu_m,nu_c=nu_c,B=B[i],p=p)

            # analytic
            P_total[i] = 0.
            if (nu_m <= nu_c):
                term0 = (3./4.)*np.power(1/nu_m,1/3)*(np.power(nu_m,4./3.)-np.power(nu_arr[0],4./3.))
                term1 = -1/(3-p)*(2*np.power(1./nu_m, 1/2-p/2))*(np.power(nu_m,3/2-p/2)-np.power(nu_c,3/2-p/2))
                term2 = np.power(nu_c / nu_m, -(p-1)/2) * 2/(p-2)*(nu_c*np.power(nu_M,p/2) - nu_M*np.power(nu_c,p/2))*np.power(nu_c*nu_M/nu_c,-p/2)
                P_total[i] = term0 + term1 + term2
            else:
                term0 = (3./4.)*np.power(1/nu_c,1/3)*(np.power(nu_c,4./3.)-np.power(nu_arr[0],4./3.))
                term1 = -2/np.sqrt(1/nu_c)*(np.sqrt(nu_c)-np.sqrt(nu_m))
                term2 = 1/np.sqrt(nu_m / nu_c) * 2/(p-2)*(nu_m*np.power(nu_M,p/2) - nu_M*np.power(nu_m,p/2))*np.power(nu_m*nu_M/nu_m,-p/2)
                P_total[i] = term0 + term1 + term2

            # P_total[i] = 0
            # for j in range(len(nu_arr)-1):
            #     P_total[i]+=spec_wspn(nu_arr[j],nu_m,nu_c,nu_M)*(nu_arr[j+1]-nu_arr[j])
            P_total[i] *= pmax
            dm = (m[i+1]-m[i])/(R[i+1]-R[i])
            # dm = (m[i+1]-m[i])
            # P_total[i] *= (dm / PBA.utils.cgs.mp)
            P_total[i] *= (dm / PBA.utils.cgs.mp)

        for i in range(len(Gamma)-1):
            # de = (Esh[i+1]-Esh[i])/(R[i+1]-R[i])
            de = (Esh[i+1]-Esh[i])/(R[i+1]-R[i])
            eps[i] = (P_total[i] * dt[i]) / ((de*eps_e))
            # eps[i] = (P_total[i]) / ((de*eps_e))
        return eps

    fig, axes = plt.subplots(figsize=(5,2.), ncols=1, nrows=1, sharex="all",layout='constrained',
                             # gridspec_kw=dict(height_ratios=[1,2,2])
                             )
    if not hasattr(axes,'__len__'): axes = [axes]
    # Erad = pba.GRB.get_dyn_arr(v_n="Erad2",ishell=0,ilayer=0)
    # Esh = pba.GRB.get_dyn_arr(v_n="Esh2",ishell=0,ilayer=0)
    # axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),Erad/(eps_e*Esh),color='black',ls='-')
    # axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),i_eps_an(),color='black',ls=':')

    # axes[1].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="gamma_min",ishell=0,ilayer=0),color='black',ls='-',label=r'$\gamma_{\rm m}$')
    # axes[1].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="gamma_c",ishell=0,ilayer=0),color='black',ls='--',label=r'$\gamma_{\rm c}$')
    # axes[1].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="gamma_max",ishell=0,ilayer=0),color='black',ls='-.',label=r'$\gamma_{\rm M}$')

    mom_smi =  pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0)
    if (pba.GRB.opts["do_rs"] == "yes"):
        g = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=0)
        b = PBA.utils.get_beta(g)
        mom_smi_rs = g*b

    # axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0),color='black',ls='--',label=r'$Semi-radiative$')

    pba = run(working_dir=working_dir,struct=struct,
              P=d2d(default=pp, new=dict(grb=dict(epsilon_e_rad=0,epsilon_e_rad_rs=0))),
              type="a")

    mom_adi =  pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0)
    if (pba.GRB.opts["do_rs"] == "yes"):
        g = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=0)
        b = PBA.utils.get_beta(g)
        mom_adi_rs = g*b
    # axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0),color='black',ls='-',label=r'$Adiabatic$')

    pba = run(working_dir=working_dir,struct=struct,
              P=d2d(default=pp, new=dict(grb=dict(epsilon_e_rad=1,epsilon_e_rad_rs=1))),
              type="a")

    mom_rad = pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0)
    if (pba.GRB.opts["do_rs"] == "yes"):
        g = pba.GRB.get_dyn_arr(v_n="Gamma43",ishell=0,ilayer=0)
        b = PBA.utils.get_beta(g)
        mom_rad_rs = g*b

    # axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
    #              pba.GRB.get_dyn_arr(v_n="mom",ishell=0,ilayer=0),color='black',ls='-.',label=r'$Radiative$')

    axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
                 mom_adi/mom_smi,color='black',ls='--',label=r'$\frac{\Gamma\beta|_{\rm adi}}{\Gamma\beta|_{\rm semi}}$')
    axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
                 mom_adi/mom_rad,color='black',ls='-.',label=r'$\frac{\Gamma\beta|_{\rm adi}}{\Gamma\beta|_{\rm rad}}$')
    if (pba.GRB.opts["do_rs"] == "yes"):
        axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
                     mom_adi_rs/mom_smi_rs,color='gray',ls='--',label=r'$\frac{\Gamma_{43}\beta_{43}|_{\rm adi}}{\Gamma_{43}\beta_{43}|_{\rm semi}}$')
        axes[0].plot(pba.GRB.get_dyn_arr(v_n="R",ishell=0,ilayer=0),
                     mom_adi_rs/mom_rad_rs,color='gray',ls='-.',label=r'$\frac{\Gamma_{43}\beta_{43}|_{\rm adi}}{\Gamma_{43}\beta_{43}|_{\rm rad}}$')

    axes[0].set_xlabel('$R$ [cm]',fontsize=12)
    axes[0].legend(fancybox=True, loc='upper left', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 2,
              fontsize=16,
              framealpha=0., borderaxespad=0.)
    axes[0].set_ylabel("Momentum ratio",fontsize=12)
    for i, ax in enumerate(axes):
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        if (f"ylim{i+1}" in plot.keys()): ax.set_ylim(*plot[f"ylim{i+1}"])
        if ("xlim" in plot.keys()): ax.set_xlim(*plot["xlim"])
        ax.grid(ls=':')
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()

def plot_id(struct:dict,pp:dict,plot:dict):
    pba = run(working_dir=working_dir,struct=struct,P=pp,type='a',do_run=False)
    fig,axes = plt.subplots(ncols=1,nrows=2,sharex='all',layout='constrained')

    # plot energy
    print( pba.GRB.get_id_obj().attrs.keys() )
    print( pba.GRB.get_id_obj().keys() )
    print( 'theta', np.array( pba.GRB.get_id_obj()['theta'] ) )

    theta_c = float( pba.GRB.get_id_obj().attrs["theta_c"] )
    theta_w = float( pba.GRB.get_id_obj().attrs["theta_w"] )

    theta = np.array( pba.GRB.get_id_obj()['theta'] )
    ctheta = np.array( pba.GRB.get_id_obj()['ctheta'] )
    ek = np.array( pba.GRB.get_id_obj()['ek'] )
    mom = np.array( pba.GRB.get_id_obj()['mom'] )

    ax = axes[0]
    ax.set_ylabel(r"$E_{0}$ [ergs]",fontsize=12)
    for i in range(len(theta)):
        if not i in plot['layers']:
            ax.plot(ctheta[i], ek[i], marker='_',markersize=20, color='black',ls='none')#, ls='-',drawstyle='steps-mid')
    for i, il in enumerate(plot['layers']):
        ax.plot(ctheta[il], ek[il], marker='_',markersize=20,color=plot['colors'][i],ls='none')#, ls='-',drawstyle='steps-mid')
    for i in range(len(theta)):
        ax.axvline(x=theta[i],color='gray',linestyle='dotted',linewidth=0.5)

    ax = axes[1]
    ax.set_ylabel(r"$\Gamma_{0} \beta_{0}$",fontsize=12)
    for i in range(len(theta)):
        if not i in plot['layers']:
            ax.plot(ctheta[i], mom[i], marker='_',markersize=20, color='black',ls='none')#, ls='-',drawstyle='steps-mid')
    for i, il in enumerate(plot['layers']):
        ax.plot(ctheta[il], mom[il], marker='_',markersize=20,color=plot['colors'][i],ls='none')#, ls='-',drawstyle='steps-mid')
    for i in range(len(theta)):
        ax.axvline(x=theta[i],color='gray',linestyle='dotted',linewidth=0.5)


    for i, ax in enumerate(axes):
        ax.set_yscale('log')
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12,
                       direction='in',
                       bottom=True, top=True, left=True, right=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        ax.minorticks_on()
        if (f"ylim{i+1}" in plot.keys()): ax.set_ylim(*plot[f"ylim{i+1}"])

        ax.axvline(x=theta_c,color='black',linestyle='dashed',linewidth=1.5, label=r'$\theta_{\rm c}$ [rad]')
        # ax.axvline(x=theta_w,color='gray',linestyle='solid',linewidth=1.5, label=r'$\theta_{\rm w}$ [rad]')
        ax.set_xlim(0,theta_w)
    axes[0].legend(fancybox=True, loc='lower left', columnspacing=1.2,
                       # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                       shadow=False, ncol= 1,
                       fontsize=12,
                       framealpha=0., borderaxespad=0.)
    axes[-1].set_xlabel("$\omega_{0}$ [rad]",fontsize=12)
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()




if __name__ == '__main__':
    ''' -------- TOPHAT ---------- '''
    # plot_fs_energy(
    #     struct = dict(struct="tophat",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
    #     pp = dict(main=dict(n_ism = 1, tb0=3e3),
    #               grb=dict(save_dynamics='yes',do_mphys_in_situ="no",do_lc = "no",# method_spread='None'
    #                        )),
    #     plot=dict(figname = "tophat_fs_energy", text="FS",
    #               xlim=(1e14,1e19), ylim1=(1e-4,2), ylim2=(1e-3,1e3), rdec=True, bm=True,
    #               theta_spread_0=True, theta_spread_1=True)
    # )
    # plot_fs_energy(
    #     struct = dict(struct="tophat",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
    #     pp = dict(main=dict(n_ism = 1., tb0=3e3, ntb=1000,rtol=1e-7,
    #                         lc_freqs = "array 1e9 1e18"),
    #               grb=dict(save_dynamics='yes',do_rs='yes',bw_type='fsrs',do_mphys_in_situ="no",do_lc = "no",do_rs_radiation="no",
    #                        # method_spread='AFGPY'
    #                        # exponential_rho4='no'
    #                        )),
    #     plot=dict(figname = "tophat_fsrs_energy", text="FS \& RS",
    #               xlim=(1e14,1e19), ylim1=(1e-4,2), ylim2=(1e-3,1e3), rdec=False, bm=True,method_ele_fs='mix',
    #               theta_spread_0=True, theta_spread_1=True)
    # )
    ''' ---------- RAD.LOSSES -------- '''
    # plot_fs_energy_rad(
    #     struct = dict(struct="tophat",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
    #     pp = dict(main=dict(n_ism = 1., tb0=3e3, ntb=3000,rtol=1e-7,
    #                         lc_freqs = "array 1e9 1e18"),
    #               grb=dict(save_dynamics='yes',#do_rs='yes',bw_type='fsrs',
    #                        do_mphys_in_situ="yes",do_lc = "no",do_rs_radiation="no",
    #                        method_gamma_min_fs='useU_e',
    #                        method_gamma_min_rs='useU_e',
    #                        method_ele_fs='analytic',method_synchrotron_fs='Joh06',
    #                        method_ele_rs='analytic',method_synchrotron_rs='Joh06',
    #                        eps_e_fs=0.1, eps_b_fs=0.001, p_fs=2.2,
    #                        eps_e_rs=0.1, eps_b_rs=0.001, p_rs=2.2,
    #                        gamma_max_fs=4e7, method_gamma_max_fs="useConst",
    #                        gamma_max_rs=4e7, method_gamma_max_rs="useConst",
    #                        ebl_tbl_fpath="none",method_spread='None'
    #                        )),
    #     plot=dict(figname = "tophat_fs_rad_momentum_ratio", text="FS \& RS",
    #               xlim=(1e14,1e19),
    #               # ylim1=(1e-3,2), ylim2=(1e-1,1e9),
    #               ylim1=(0.9,2),
    #               rdec=False, bm=True,method_ele_fs='mix',
    #               theta_spread_0=True, theta_spread_1=True)
    # )

    ''' ---------- GAUSSIAN --------- '''
    plot_id(
        struct = dict(struct="gaussian",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.3),
        pp = dict(main=dict(),
                  grb=dict()),
        plot=dict(figname = "gaussian_structure_adaptive", text="Structure",
                  xlim=(1e14,1e19), ylim1=(1e50,3e53), ylim2=(1e0,1e3), ylim3=(0.,np.pi/2.),#ylim4=(1e-10,np.pi/2.),
                  layers=(0,5,10,15,20), colors=("red","orange","yellow","green","blue"))
    )

    plot_fs_energy2(
        struct = dict(struct="gaussian",Eiso_c=1.e53, Gamma0c= 400., M0c= -1.,theta_c= 0.1, theta_w= 0.3),
        pp = dict(main=dict(n_ism = 1., tb0=3e3, ntb=2000,rtol=1e-7,
                            lc_freqs = "array 1e9 1e18"),
                  grb=dict(save_dynamics='yes',do_rs='yes',bw_type='fsrs',do_mphys_in_situ="no",do_lc = "no",do_rs_radiation="no",
                           # method_spread='AFGPY'
                           # exponential_rho4='no'
                           )),
        plot=dict(figname = "gaussian_fsrs_energy", text="FS \& RS",
                  xlim=(1e14,1e19), ylim1=(1e-1,10), ylim2=(1e-3,1e3), ylim3=(0.,np.pi/2.),#ylim4=(1e-10,np.pi/2.),
                  rdec=False, bm=True,method_ele_fs='mix',
                  layers=(0,5,10,15,20), colors=("red","orange","yellow","green","blue"),
                  theta_spread_0=True, theta_spread_1=True)
    )
