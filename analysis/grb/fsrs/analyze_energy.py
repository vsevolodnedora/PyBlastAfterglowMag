import package.src.PyBlastAfterglowMag as PBA
import os,shutil,matplotlib.pyplot as plt
import numpy as np
working_dir = os.getcwd()+'/tmp1/'
fig_dir = os.getcwd()+'/figs/'
def gamma_adi(Gamma, beta):
    """ Adiabatic index of the fluid From Nava 2013 paper """
    return (4. + 1. / Gamma) / 3.

def GammaEff(Gamma, gammaAdi):
    return (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

def run(working_dir:str, struct:dict, P:dict, type:str="a") -> PBA.PyBlastAfterglow:
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
        Gamma_bm = get_bm79(E0=struct["Eiso_c"]/2, A0=pba.main_pars["n_ism"]*PBA.utils.cgs.mp, s=0, R=R)
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
    axes[-1].set_xlabel(r"$R$ [cm]", fontsize=12)
    plt.tight_layout()
    plt.savefig(fig_dir+plot["figname"]+".pdf")
    plt.savefig(fig_dir+plot["figname"]+".png", dpi=256)
    plt.show()

if __name__ == '__main__':
    # plot_fs_energy(
    #     struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
    #     pp = dict(main=dict(n_ism = 10),
    #               grb=dict(save_dynamics='yes',do_ele = "no",do_lc = "no",# method_spread='None'
    #                        )),
    #     plot=dict(figname = "tophat_fs_energy", text="FS",
    #               xlim=(1e14,1e19), ylim1=(1e-4,2), ylim2=(1e-3,1e3), rdec=True, bm=True,
    #               theta_spread_0=True, theta_spread_1=True)
    # )
    plot_fs_energy(
        struct = dict(struct="tophat",Eiso_c=1.e52, Gamma0c= 350., M0c= -1.,theta_c= 0.1, theta_w= 0.1),
        pp = dict(main=dict(n_ism = 100., tb0=1e4, ntb=1000,rtol=1e-7,
                            lc_freqs = "array 1e9 1e18"),
                  grb=dict(save_dynamics='yes',do_rs='yes',bw_type='fsrs',do_ele = "no",do_lc = "no",do_rs_radiation="no",
                           #method_spread='None'
                           )),
        plot=dict(figname = "tophat_fsrs_energy", text="FS \& RS",
                  xlim=(1e14,1e19), ylim1=(1e-4,2), ylim2=(1e-3,1e3), rdec=False, bm=True,method_ele_fs='mix',
                  theta_spread_0=True, theta_spread_1=True)
    )