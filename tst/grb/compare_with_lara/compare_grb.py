# import PyBlastAfterglowMag
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt
from click.core import F
from matplotlib.patches import Ellipse
import os

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs

try:
    import afterglowpy as grb
except:
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

try:
    import jetsimpy
except:
    print("Error! could not import jetsim")

# temprary dir to save PyBlastAfterglow output into
working_dir = os.getcwd() + '/working_dirs/'
fig_dir = os.getcwd() + '/figs/'

def compare_spec(pp:dict, plot:dict,working_dir:str, run:bool=True,fs_or_rs:str="fs",time=1e4):

    pba = PBA.wrappers.run_grb(working_dir=working_dir+f'num_spec/', run=run, P=pp)
    freqs = pba.GRB.get_lc_freqs()
    spec = pba.GRB.get_lc(time=time)
    spec = freqs*spec*1e-3*1e-23

    # dr = pba.GRB.get_dyn_arr(v_n="thickness" if fs_or_rs == "fs" else "thichness_rs", ishell=0, ilayer=0)
    # Gamma = pba.GRB.get_dyn_arr(v_n="Gamma" if fs_or_rs == "fs" else "Gamma43", ishell=0, ilayer=0)
    # dr_comov = dr * Gamma  # thickness as R / Gamma in comoving frame (dr is in observer frame so R/Gamma^2
    # spec = pba.GRB.get_lc(key="syn_a" + '_' + fs_or_rs, xkey="freqs", key_time="times_freqs", freq=None, time=None,
    #                  ishell=0, ilayer=0,
    #                  spec=True, sum_shells_layers=False)
    # times = pba.GRB.get_grid(key="times_freqs", spec=True)
    # ys = pba.GRB.get_grid(key="freqs", spec=True)
    # tau = spec * dr_comov[np.newaxis, :]
    #
    # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.2, 3.4))
    # ax.plot(
    #     times,
    #     [ys[PBA.utils.find_nearest_index(tau[:, i], 1.)] for i in range(len(Gamma))],  # [xs,ys]
    #     color='black', lw=0.8, ls='-',
    #     label=r"$\tau_{\nu'}'=1$"
    # )
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # plt.show()
    #
    #
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.2, 3.4))

    # extracted form the paper
    ref_syn = np.array([
        [103226951.858901,	3.8496755945462E-22],
        [1662224407.9925,	9.20783876864956E-19],
        [195604711.04249,	2.54107756086457E-21],
        [402876943.692736,	2.03349789698026E-20],
        [763411364.353909,	1.29154966501488E-19],
        [4044795338.20542,	4.2974799486962E-18],
        [10698190939.1319,	1.71935041046893E-17],
        [32513892501.3145,	7.14893483124243E-17],
        [98816071951.023,	2.86017018886645E-16],
        [308784272767.175,	1.10107273723278E-15],
        [1346822443349.93,	4.07862464897348E-15],
        [7756310250274.17,	9.51688066627356E-15],
        [44668359215096.2,	1.76245113527495E-14],
        [287488640837590,	3.26392030451741E-14],
        [1610262027560936,	5.59640581945151E-14],
        [9.53477271083525E+015,	8.54869143640456E-14],
        [6.6701857930592E+016,	1.25650343388445E-13],
        [4.53833501334105E+017,	1.46578169930112E-13],
        [3.00321349571879E+018,	1.46578169930112E-13],
        [2.28360344579806E+019,	1.25650343388445E-13],
        [1.24402064072824E+020,	1.03641015746956E-13],
        [6.41056859642021E+020,	8.22570703582356E-14],
        [2.36667122805298E+021,	6.04452262024225E-14],
        [1.06135841841288E+022,	3.52526974415465E-14],
        [3.41003330518472E+022,	1.57013817704677E-14],
        [8.29784880405789E+022,	4.75794431400939E-15],
        [1.6622244079925E+023,	8.09104802511268E-16],
        [2.74113884733337E+023,	9.36108448607985E-17],
        [4.04479533820542E+023,	8.27108957121876E-18],
        [5.49104714405709E+023,	5.3701702284098E-19],
        [7.25011386916185E+023,	2.56213583397132E-20],
        [8.80699168290405E+023,	1.17622348016418E-21],
        [1.04049831036579E+024,	1.12248155037532E-22]
    ])
    ref_spn_syn = np.array([
        [80384518.3319958,	1.25996483093298E-22],
        [2592943797.40467,	5.41467359699691E-18],
        [553481932057.804,	5.99484250318942E-15],
        [1.26393324723418E+017,	4.83713060388845E-13],
        [3.13726875419839E+022,	7.32814453215123E-14]
    ])
    ref_ssc = np.array([
        [34783983015347.1,	1.12248155037532E-22],
        [34783983015347.1,	1.12248155037532E-22],
        [118144577379989,	1.0082871429944E-21],
        [118144577379989,	1.0082871429944E-21],
        [625967411718142,	1.55295510542585E-20],
        [625967411718142,	1.55295510542585E-20],
        [2655448793140943,	1.29154966501488E-19],
        [1.00797155358033E+016,	7.03190312955176E-19],
        [3.82612029820155E+016,	3.54471921861865E-18],
        [1.623098317592E+017,	1.65439036765474E-17],
        [9.88160719510226E+017,	8.02454745290512E-17],
        [6.18556407953358E+018,	2.75210799468352E-16],
        [4.32719928378322E+019,	7.49120999259122E-16],
        [2.78501530521364E+020,	1.88792902380624E-15],
        [2.05964846933173E+021,	4.4052092978157E-15],
        [1.32560396096054E+022,	9.15731639608505E-15],
        [7.42488369523096E+022,	1.69586267226692E-14],
        [3.72125691921704E+023,	2.39843700066182E-14],
        [2.32938556206052E+024,	3.14060373016904E-14],
        [1.77123368142049E+025,	3.26392030451741E-14],
        [9.64901004723305E+025,	2.30781988776765E-14],
        [4.83596065281448E+026,	1.11019748445346E-14],
        [1.64254486698933E+027,	5.34070471018872E-15],
        [5.27732442172512E+027,	2.20237504972001E-15],
        [1.47558936758722E+028,	7.49120999259122E-16],
        [3.21289514081771E+028,	1.87241207717822E-16],
        [5.60113270169382E+028,	3.43905930797929E-17]
    ])
    ref_se_ssc = np.array([
        [25622452084929.8,	1.08007224909123E-22],
        [719277626775213,	8.13568762530566E-20],
        [1.29955032877223E+017,	7.14893483124243E-17],
        [3.50612653996912E+022,	4.44171589718697E-14],
        [8.00660264807753E+027,	7.38887391377279E-13],
        [6.08811375626488E+028,	6.84109226951232E-13]
    ])
    # ax.plot(ref_syn[:,0], ref_syn[:,1],color='blue',ls='--',label='Miceli+22 (syn.)')
    # ax.plot(ref_spn_syn[:,0], ref_spn_syn[:,1],color='blue',ls=':',label='SPN98')
    # ax.plot(ref_ssc[:,0], ref_ssc[:,1],color='red',ls='--',label='Miceli+22 (ssc.)')
    # ax.plot(ref_se_ssc[:,0], ref_se_ssc[:,1],color='red',ls=':',label='SE01')

    # Ludovica's data: Syn_IC_spectrum_dens_N_e_ER_DAVIDE_fig5.txt 6,7 for Numerical 2,3 Analyical
    dfile = np.loadtxt(os.getcwd()+'/'+"Syn_IC_spectrum_dens_N_e_ER_DAVIDE_fig5.txt")


    def plot_(ax):
        ax.plot(freqs, spec,color='black',ls='-',lw=2,label='PyBlastAfterglow')
        ax.plot(dfile[:,0], dfile[:,0]*dfile[:,1],color='blue',ls='--',label='SPN98')
        ax.plot(dfile[:,0], dfile[:,0]*dfile[:,2],color='blue',ls=':',label='SE01')
        # ax.plot(dfile[:,0], dfile[:,0]*dfile[:,3],color='red',ls='--',label='Ref. Syn.')
        # ax.plot(dfile[:,0], dfile[:,0]*dfile[:,4],color='red',ls=':',label='Ref. SSC')
        ax.plot(dfile[:,0], dfile[:,0]*dfile[:,5],color='red',ls='-',label='Menegazzi et al.')
    plot_(ax)

    # inset axes....
    if plot["include_zoom_in"]:
        x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
        axins = ax.inset_axes(
            [0.05, 0.05, 0.45, 0.45],
            xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
        plot_(axins)
        # sub region of the original image
        x1, x2, y1, y2 = 1e20, 1e27, 1e-15, 1e-12
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xscale('log')
        axins.set_yscale('log')
        axins.set_xticklabels('')
        axins.set_yticklabels('')
        axins.minorticks_on()
        axins.tick_params(axis='both', which='both', labelleft=True,
                          labelright=False, tick1On=True, tick2On=True,
                          labelsize=12,
                          direction='in',
                          bottom=True, top=True, left=True, right=True)
        ax.indicate_inset_zoom(axins, edgecolor="black")


    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=1.)
    ax.text(0.95, 0.9, r"${}$ s.".format(PBA.utils.latex_float(time)), fontsize=12, bbox=bbox,
            transform=ax.transAxes, horizontalalignment='right')
    ax.legend(fancybox=False, loc='upper left', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=2, fontsize=12, labelcolor='black',
              framealpha=0.0, borderaxespad=0.)
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=13,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.grid(ls=':')
    ax.set_xlabel(r"$\nu$ [Hz]",fontsize=13)
    ax.set_ylabel(r"$\nu F_{\nu}$ [erg cm$^{-2}$ s$^{-1}$]",fontsize=13)
    ax.set_xlim(plot["xlim"])
    ax.set_ylim(plot["ylim"])
    # ax.legend()
    plt.tight_layout()
    # if save_figs: plt.savefig(PAPERPATH + save)
    if "figname" in plot.keys():
        plt.savefig(fig_dir+plot["figname"]+'.png',dpi=256)
        plt.savefig(fig_dir+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()


if __name__ == '__main__':
    show = True
    ''' VHE spectrum '''
    struct = dict(struct="tophat",Eiso_c=2*1.e52, Gamma0c= 400., M0c= -1.,theta_c= np.pi/2., theta_w= np.pi/2)
    compare_spec(
        pp = dict(main=dict(n_ism = 1.,tb0=1.,tb1=1e10,ntb=1000, d_l=6791.e6 * cgs.pc, z=1, theta_obs=0.,
                            lc_freqs = "array logspace 1e8 1e29 96", lc_times = "array logspace 10 1e6 128"),
                  grb=dict(structure=struct,eats_type='a',save_spec='yes',save_dynamics='yes',
                           eps_e_fs=0.05, eps_b_fs=5e-4,p_fs=2.3, freq1_fs=1e5, freq2_fs=1e33,
                           ebl_tbl_fpath='none',
                           method_spread='None',
                           method_ssc_fs='numeric',method_pp_fs='none',use_ssa_fs='yes',
                           # method_synchrotron_fs='WSPN99',
                           # method_ele_fs='analytic'
                           )),
        plot = dict(plot_analytic=True,plot_semi_analytic=True,plot_numeric=False,
                    #             iters=[
                    #     dict(theta_obs=0.16,freq=1.e9,ls='-'),
                    #     dict(theta_obs=0.0,freq=1.e18,ls='--'),
                    #     dict(theta_obs=0.16,freq=1.e18,ls='-.'),
                    #     dict(theta_obs=0.0,freq=1.e9,ls=':')
                    # ],
                    xlim=(1e7, 1e29), ylim=(1e-22, 1e-8),
                    include_zoom_in=False,
                    show=show, figname="spec_tophat_vhe"),
        working_dir=working_dir+"tmp_tophat_",
        run=False
    )