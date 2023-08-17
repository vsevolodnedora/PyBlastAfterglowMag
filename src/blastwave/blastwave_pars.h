//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLASTWAVE_PARS_H
#define SRC_BLASTWAVE_PARS_H

//#include "../synchrotron_an.h"
#include "../eats.h"

/* ------------- EQUATIONS ----------------- */

enum RHS_TYPES { iGRG_FS, iGRG_FSRS, iEJ, iEJ_PWN };


enum METHODS_Up { iuseEint2, iuseGamma }; // energy density behind the shock
enum METHOD_Delta { iuseJoh06, iuseVE12, iNoDelta }; // thickness of the shock
enum METHOD_GammaSh { iuseJK, isameAsGamma, iuseGammaRel, iuseJKwithGammaRel };
enum METHOD_RSh { isameAsR, iuseGammaSh };
enum METHOD_dmdr{ iusingA, iusingdthdR, iNodmdr };
enum METHOD_dgdr{ iour, ipeer };

enum METHODS_RAD { icomovspec, iobservflux };
enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };
enum METHOD_NE{ iusenprime, iuseNe };

enum METHOD_SINGLE_BW_DELTA {iconst, ifrac_last, ilr, ibm, ist, ibm_st };
enum METHOD_THICK_FOR_RHO { iFromRadius };



struct Pars{

    Pars(VecVector & data, VecVector & data_tmp, unsigned loglevel) : m_data(data), m_data_tmp(data) {
        p_log = std::make_unique<logger>(std::cout,std::cerr,loglevel,"PWNPars");
    }
    std::unique_ptr<SynchrotronAnalytic> p_syna = nullptr;
    std::unique_ptr<SynchrotronAnalytic> p_syna_rs = nullptr;

    std::unique_ptr<logger> p_log;
    Vector m_freq_arr{}; Vector m_synch_em{}; Vector m_synch_abs{}; Vector m_synch_em_rs{}; Vector m_synch_abs_rs{};
    VecVector & m_data;
    VecVector & m_data_tmp;

    // set a reference to the data container
    // *************************************** //

    // initial conditions (settings)
//    bool is_init = false;
    RHS_TYPES m_rhs{};
    METHODS_Up m_method_up{};
    METHOD_Delta m_method_Delta{};
    METHOD_GammaSh m_method_gamma_sh{};
    METHOD_RSh m_method_r_sh{};
    METHOD_dmdr m_method_dmdr{};
//    bool only_last_shell_dmdr= false;
    METHOD_dgdr m_method_dgdr{};
    EjectaID2::STUCT_TYPE m_method_eats{};

    METHOD_SINGLE_BW_DELTA method_single_bw_delta{};
    METHOD_THICK_FOR_RHO method_thick_for_rho{};
    /// blast wave initial conditions
    double M0 = -1.;
    double R0 = -1.;
    double tb0 = -1.;
    double Gamma0 = -1.;
    double beta0 = -1.;
    double mom0 = -1.;
    double E0 = -1.;
    double Eint0 = 0.;
    double Ye0 = -1.;
    double s0 = -1.;
    double theta_a = -1.;
    double theta_b0 = -1.;
    double ncells = -1.;
    double ctheta0 = -1.;
//        double theta_h0 = -1;
    double theta_c_l = -1.;
    double theta_c_h = -1.;
    double theta_w = -1.;
    double theta_max=-1.;
    /// deceleration radius
    double Rd = -1;

    double eps_rad = 0.;
    double dEinjdt = 0.;
    double facPSRdep = 0.;
    double dEnuc = 0.;
    double dElum = 0.;
    double kappa = 0.;
    double delta = 0.;
    double vol = 0.;
    double _last_frac=0.;// last time when Nshell>1; delta/delta
    double _last_frac_vol=0.;// last time when Nshell>1; vol/vol
    // ---
    bool adiabLoss = true;
    // ---
//    size_t comp_ix = 0;
    size_t nr = -1;
    size_t ilayer = 0;
    size_t ishell = 0;
    size_t nlayers = 0;
    size_t nshells = 0;
    size_t ii_eq = 0;
    size_t i_end_r = 0.; // index of R when evolution has stopped (for EATS)
    unsigned loglevel = 0;

    // --- PARS FOR INTERACTION WITH OTHER MODELS
    bool entry = false; // has this BW entered the 'void' left by another
    double entry_time = -1; // time when this BW entered the 'void' left by another
    double entry_r = -1; // same for radius
    // --
    double prev_x=-1; size_t prev_idx_x=0; size_t n_substeps=10; // for Derivatives, logging
    // ---
    size_t ijl = 123456;
    size_t prev_ijl = 123456;
    double first_entry_r = -1;
    double r_dist = -1;
    bool is_within = false;
    bool is_within0 = false;
    bool is_using_st_prof = false;

    bool allow_termination = false;
    bool end_evolution = false;
    bool end_spreading = false;
    bool spread_only_after_rd = false;
    double min_beta_terminate = 1.e-8;
    // ---
    double fraction_of_Gamma0_when_spread = -1.;

    /// main switch
    bool use_dens_prof_behind_jet_for_ejecta = false;
    bool use_dens_prof_inside_ejecta = false;
    int which_jet_layer_to_use = 0;

    /// exp decay && floor behaviour
    bool use_exp_rho_decay_as_floor = true;
    bool use_flat_dens_floor = false;
    double steepnes_of_exp_decay = 1.;
//    double dens_floor_frac = 1e-20;

    /// ST profile
    bool use_st_dens_profile = true;
    double Gamma_when_st_starts = 2.;

    /// BM profile
    bool use_bm_dens_profile = false;
    double fraction_of_Gamma0_when_bm_for_bm = 1.98; // if > 1 -- BM also in free coasting

    long n_fialed_electrons = 0;
    size_t i0_failed_elecctrons = 0;

    /// ----------------------------
    double d_l=-1.;
    double theta_obs=-1.;
    double z=-1.;
    METHODS_SHOCK_VEL method_shock_vel{};
    METHOD_NE m_method_ne{};
    METHODS_RAD m_method_rad{};
    /// ---------------------------
    int opacitymode = -1;
    double Z_eff=-1.; // electron mean molecular weight
    double mu_e=-1.; // effective nuclear weight
    double albd_fac = -1.;
    double mom_when_bm_to_st = 0.;
    double fraction_last_delta = 0.;
    double _b0_delta = -1., _b1_delta = -1.;
    double _b0_vol = -1., _b1_vol = -1.;
    Vector m_mu{};
    Vector ttobs{};

    // ------- PWN --------------------
    double eps_e_w{};
    double eps_mag_w{};
    double epsth_w{};
    double gamma_b_w{};
    double curr_ldip=-1.;
    double curr_lacc=-1.;
    double curr_b_pwn=-1.;


    // --------- RS ---------------
    double rs_Gamma0_frac_no_exceed = 0.98;
    bool do_rs = false;
    bool shutOff = false; /// include Reverse shock now?
    bool adiabLoss_rs = false;
    double tprompt = 0.;
    double epsilon_rad_rs = 0.;
    double rs_shutOff_criterion_rho = 0.;
    double prev_dGammadR=0;
    double prev_dM3dR=0;
    double prev_dRdt = 0;
    double min_Gamma0_for_rs = 0.;
};

/// Each BlastWave collects the following data
namespace BW{
    enum Q {
        // -- dynamics ---
        iR, iRsh, irho, idrhodr, iGammaCBM, iGammaREL, idGammaCBMdr, idGammaRELdGamma, iPcbm, idPCBMdrho, iMCBM, iCSCBM,
        iGamma, ibeta, imom, iEint2, iU_p, itheta, ictheta, iErad2, iEsh2, iEad2, iM2,
        iEint3, iErad3, iEad3, iEsh3, iM3, ideltaR4, iGamma43, iadi3, irho4, irho3, ithichness_rs, iU_p3,
        itcomov, itburst, itt, ithickness, iadi, irho2, iGammaFsh,
        // --- electrons  ---
        igm, igc, igM, iB, iTheta, iz_cool, ix, inprime, iacc_frac,
        igm_rs, igc_rs, igM_rs, iB3, iTheta_rs, iz_cool_rs, ix_rs, inprime_rs, iacc_frac_rs,
        // --- observables ---
        imu,
        // ---
        ijl, ir_dist,
        // --- energy injection
        ipsrFrac, iLmag, iLnuc,
        // --- properties of the ejecta element
        iEJr, iEJrho, iEJbeta, iEJdelta, iEJvol, iEJdtau, iEJtaucum, iEJtaucum0, iEJeth, iEJtemp, iEJlum, iEJtdiff,
        // -- PWN
        i_Rw, i_Wmom, i_Wenb, i_Wepwn, i_Wtt, i_Wb, i_WGamma, i_Wdr
    };
    const std::vector<std::string> m_vnames{
            // --- dynamics ---
            "R", "Rsh", "rho", "drhodr", "GammaRho", "GammaRel", "dGammaRhodr", "dGammaReldGamma", "PCBM", "dPCBMdrho", "MCBM", "CSCBM",
            "Gamma", "beta", "mom", "Eint2", "U_p", "theta", "ctheta", "Erad2", "Esh2", "Ead2", "M2",
            "Eint3", "Erad3", "iEad3", "Esh3", "M3", "deltaR4", "Gamma43", "adi3", "rho4", "rho3", "thichness_rs", "U_p3",
            "tcomov", "tburst", "tt", "thickness", "adi", "rho2", "GammaFsh",
            // --- electrons
            "gamma_min", "gamma_c", "gamma_max", "B", "ThetaTherm", "z_cool", "x", "nprime", "accel_frac",
            "gamma_min_rs", "gamma_c_rs", "gamma_max_rs", "B_rs", "ThetaTherm_rs", "z_cool_rs", "x_rs", "nprime_rs", "accel_frac_rs",
            // --- observables
            "mu",
            // ---
            "ijl", "r_dist",
            // --- energy injection
            "psrFrac", "Lmag", "Lnuc",
            // --- properties of the ejecta as a part of the cumShell
            "EJr", "EJrho", "EJbeta", "EJdelta", "EJvol", "EJdtau", "EJtaucum", "EJtaucum0", "EJeth", "EJtemp", "EJlum", "EJtdiff",
            // --- FOR PWN ---
            "WR", "Wmom", "Wenb", "Wepwn", "Wtt", "Wb", "WGamma", "Wdr"
    };
    static constexpr size_t NVALS = 87; // number of variables to save
//    static constexpr size_t DERIVATIVES_NVALS = 14; // number of variables to save for derivatives at EACH step

    /// ---
    enum QSH { iRs, iGammas, iDeltas, iVols, irhos, iEints, iPs};
    static const size_t NVALS_SH = 7;

}

/// Each blastwave evolution computes the following data
namespace SOL{
//    enum QS { iR, iRsh, itt, itcomov, imom, iEint2, itheta, iErad2, iEsh2, iEad2, iM2,
//             iRw, iWmom, iWenb, iWepwn, iWtt};
//    std::vector<std::string> vars { "R", "Rsh", "tt", "tcomov", "mom",
//                                    "Eint2", "theta", "Erad2", "Esh2", "Ead2", "M2",
//                                    "Rw", "momW", "enbW", "EpwnW", "ttW"};

    enum QS { iR, iRsh, itt, itcomov, iGamma, iEint2, iEint3, itheta,
            iErad2, iErad3, iEsh2, iEsh3, iEad2, iEad3, iM2, iM3, ideltaR4,
        iRw, iWmom, iWenb, iWepwn, iWtt};
    std::vector<std::string> vars { "R", "Rsh", "tt", "tcomov", "Gamma",
                                    "Eint2", "Eint3", "theta",
                                    "Erad2", "Erad3", "Esh2", "Esh3", "Ead2", "Ead3", "M2", "M3", "deltaR4",
                                    "Rw", "momW", "enbW", "EpwnW", "ttW"};

    const static int neq = 22;
}

#endif //SRC_BLASTWAVE_PARS_H
