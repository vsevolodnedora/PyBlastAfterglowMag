//
// Created by vsevolod on 05/04/23.
//

#ifndef SRC_BLASTWAVE_H
#define SRC_BLASTWAVE_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "blastwave_components.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"
#include "composition.h"

#include "eats.h"

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

struct Pars{

    Pars(VecVector & data,unsigned loglevel) : m_data(data){
        p_log = std::make_unique<logger>(std::cout,std::cerr,loglevel,"PWNPars");
    }
    std::unique_ptr<SynchrotronAnalytic> p_syna = nullptr;
    std::unique_ptr<logger> p_log;
    Vector m_freq_arr{}; Vector m_synch_em{}; Vector m_synch_abs{};
    VecVector & m_data;

    // set a reference to the data container
    // *************************************** //

    // initial conditions (settings)
//    bool is_init = false;
    METHODS_Up m_method_up{};
    METHOD_Delta m_method_Delta{};
    METHOD_GammaSh m_method_gamma_sh{};
    METHOD_RSh m_method_r_sh{};
    METHOD_dmdr m_method_dmdr{};
//    bool only_last_shell_dmdr= false;
    METHOD_dgdr m_method_dgdr{};
    EjectaID2::STUCT_TYPE m_method_eats{};

    METHOD_SINGLE_BW_DELTA method_single_bw_delta{};
    /// blast wave initial conditions
    double M0 = -1.;
    double R0 = -1.;
    double tb0 = -1.;
    double Gamma0 = -1.;
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
    size_t comp_ix = 0;
    size_t nr = -1;
    size_t ilayer = 0;
    size_t ishell = 0;
    size_t nlayers = 0;
    size_t ii_eq = 0;
    size_t i_end_r = 0.; // index of R when evolution has stopped (for EATS)
    unsigned loglevel = 0;

    // --- PARS FOR INTERACTION WITH OTHER MODELS
    bool entry = false; // has this BW entered the 'void' left by another
    double entry_time = -1; // time when this BW entered the 'void' left by another
    double entry_r = -1; // same for radius
    // --
    double prev_x=-1;
    // ---
    size_t ijl = 123456;
    size_t prev_ijl = 123456;
    double first_entry_r = -1;
    double r_dist = -1;
    bool is_within = false;
    bool is_within0 = false;
    bool is_using_st_prof = false;

    bool end_evolution = false;
    bool end_spreading = false;
    double min_beta_terminate = 1.e-8;
    // ---
    double fraction_of_Gamma0_when_spread = -1.;

    /// main switch
    bool use_dens_prof_behind_jet_for_ejecta = false;
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
};

/// Each BlastWave collects the following data
namespace BW{
    enum Q {
        // -- dynamics ---
        iR, iRsh, irho, idrhodr, iGammaCBM, iGammaREL, idGammaCBMdr, idGammaRELdGamma, iPcbm, idPCBMdrho, iMCBM, iCSCBM,
        iGamma, ibeta, imom, iEint2, iU_p, itheta, ictheta, iErad2, iEsh2, iEad2, iM2,
        itcomov, itburst, itt, ithickness, iadi, irho2, iGammaFsh,
        // --- electrons  ---
        igm, igc, igM, iB, iTheta, iz_cool, ix, inprime, iacc_frac,
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
    std::vector<std::string> m_vnames{
            // --- dynamics ---
            "R", "Rsh", "rho", "drhodr", "GammaRho", "GammaRel", "dGammaRhodr", "dGammaReldGamma", "PCBM", "dPCBMdrho", "MCBM", "CSCBM",
            "Gamma", "beta", "mom", "Eint2", "U_p", "theta", "ctheta", "Erad2", "Esh2", "Ead2", "M2",
            "tcomov", "tburst", "tt", "thickness", "adi", "rho2", "GammaFsh",
            // --- electrons
            "gamma_min", "gamma_c", "gamma_max", "B", "ThetaTherm", "z_cool", "x", "nprime", "accel_frac",
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
    static constexpr size_t NVALS = 66; // number of variables to save
}

/// Each blastwave evolution computes the following data
namespace SOL{
    enum QS { iR, iRsh, itt, itcomov, imom, iEint2, itheta, iErad2, iEsh2, iEad2, iM2,
              iRw, iWmom, iWenb, iWepwn, iWtt};
    std::vector<std::string> vars { "R", "Rsh", "tt", "tcomov", "mom",
                                    "Eint2", "theta", "Erad2", "Esh2", "Ead2", "M2",
                                    "Rw", "momW", "enbW", "EpwnW", "ttW"};
    const static int neq = 16;
}

/// Main Blastwave class
class BlastWave{
    std::unique_ptr<logger> p_log;
protected:
    Vector m_tb_arr;
    VecVector m_data{}; // container for the solution of the evolution
    Pars * p_pars = nullptr;
    std::unique_ptr<LatSpread> p_spread = nullptr;
    std::unique_ptr<EOSadi> p_eos = nullptr;
    std::unique_ptr<RhoISM> p_dens = nullptr;
    std::unique_ptr<SedovTaylor> p_sedov = nullptr;
    std::unique_ptr<BlandfordMcKee2> p_bm = nullptr;
    std::unique_ptr<NuclearAtomic> p_nuc = nullptr;
//    std::unique_ptr<BlastWaveRadiation> p_rad = nullptr;
//    std::unique_ptr<PWNPars> p_pars = nullptr;
//    std::unique_ptr<EatsPars> p_eats_pars = nullptr;
    std::unique_ptr<EATS> p_eats_fs = nullptr;
    std::unique_ptr<EATS> p_eats_opt_depth = nullptr;
    std::unique_ptr<LinearRegression> p_lr_delta = nullptr;
    std::unique_ptr<LinearRegression> p_lr_vol = nullptr;
//    int m_loglevel;
//    size_t ish = 0;
//    size_t il = 0;
    static constexpr int iters=1000;
    Vector frac_psr_dep_{iters, 0.};
public:
    bool is_initialized = false;

    enum CASES { i_INSIDE_BEHIND, i_OUTSIDE_BEHIND, i_INSIDE_ABOVE, i_OUTSIDE_ABOVE, i_AT_ZERO_INSIDE };
    BlastWave(Vector & tb_arr, size_t ishell, size_t ilayer, int loglevel ) : m_tb_arr(tb_arr){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW");
        /// First: resize the container
        if (m_data.empty()){
            m_data.resize( BW::NVALS );
        }
        /// Check if contener will be filled by evolving or loading
        if (m_tb_arr.empty()){
            // REMOVING LOGGER
            (*p_log)(LOG_WARN,AT) << " Time grid was not initialized\n";
//            std::cerr << AT  << "\n";
//            exit(1);
        }
        /// if no evolution required; do not allocate memory for each variable
        if (m_data[BW::Q::itburst].size() < 1) {
            for (auto & arr : m_data) {
                arr.resize( tb_arr.size(), 0.0);
            }
        }
        // ---------------------- Methods
//        p_pars = std::make_unique<PWNPars>(); //
        p_pars = new Pars(m_data, loglevel); //
        p_lr_delta = std::make_unique<LinearRegression>(m_data[BW::Q::iR],m_data[BW::Q::iEJdelta]);
        p_lr_vol = std::make_unique<LinearRegression>(m_data[BW::Q::iR],m_data[BW::Q::iEJvol]);
        p_spread = std::make_unique<LatSpread>();
        p_eos = std::make_unique<EOSadi>();
        p_dens = std::make_unique<RhoISM>(loglevel);
        p_sedov = std::make_unique<SedovTaylor>();
        p_bm = std::make_unique< BlandfordMcKee2>();
        p_nuc = std::make_unique<NuclearAtomic>(loglevel);
//        p_rad = std::make_unique<BlastWaveRadiation>(m_data, ishell, ilayer, loglevel);
        p_pars->p_syna = std::make_unique<SynchrotronAnalytic>(loglevel);
        /// EATS integrator for forward shock
//        p_eats_pars = new EatsPars(m_data); /// for static methods (data link)
        p_eats_fs = std::make_unique<EATS>(m_data[BW::Q::itburst],
                                           m_data[BW::Q::itt], m_data[BW::Q::iR], m_data[BW::Q::itheta],
                                           m_data[BW::Q::iGamma],m_data[BW::Q::ibeta],
                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
//        p_eats_opt_depth = std::make_unique<EATS>(m_data[BW::Q::itburst],
//                                           m_data[BW::Q::itt], m_data[BW::Q::iR], m_data[BW::Q::itheta],
//                                           m_data[BW::Q::iGamma],m_data[BW::Q::ibeta],
//                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
//                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
        p_eats_fs->setFluxFunc(fluxDensPW);
        p_eats_fs->setFluxFuncA(fluxDensA);
//        p_eats_fs->setFuncOptDepth(optDepthPW);
        /// ----------------------
        p_pars->loglevel = loglevel;
        p_pars->nr = m_tb_arr.size();
        p_pars->ilayer = ilayer;
        p_pars->ishell = ishell;
//        ish = ishell;
//        il = ilayer;
        is_initialized = true;
    }
    ~BlastWave(){delete p_pars;}
    void setParams(std::unique_ptr<EjectaID2> & id, StrDbMap & pars, StrStrMap & opts, size_t ilayer, size_t ii_eq){

        double nism, A0, s, r_ej, r_ism,  a, theta_max, epsilon_e_rad;

        // set parameters for ISM density
        nism = getDoublePar("n_ism", pars, AT, p_log, -1, true);//pars.at("nism");
        A0 = getDoublePar("A0", pars, AT,p_log,-1,false);//pars.at("A0");
        s = getDoublePar("s", pars, AT,p_log,-1,false);//pars.at("s");
        r_ej = getDoublePar("r_ej", pars, AT,p_log,-1,false);//pars.at("r_ej");
        r_ism = getDoublePar("r_ism", pars, AT,p_log,-1,false);//pars.at("r_ism");
        p_dens->setPars(nism, A0, s, r_ej, r_ism, true);

        // spreading
        a = getDoublePar("a", pars, AT,p_log,-1,false);//pars.at("a");
        theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);//pars.at("theta_max");

        // radiative losses
        epsilon_e_rad = getDoublePar("epsilon_e_rad", pars, AT,p_log,0.,false);// pars.at("epsilon_e_rad");

        // interaction parameters
        p_pars->which_jet_layer_to_use =
                (int)getDoublePar("which_jet_layer_to_use",pars,AT,p_log,1.e5,false);
        p_pars->steepnes_of_exp_decay =
                (double)getDoublePar("steepnes_of_exp_decay",pars,AT,p_log,1.,false);
        p_pars->Gamma_when_st_starts =
                (double)getDoublePar("Gamma_when_st_starts",pars,AT,p_log,2.,false);
        p_pars->fraction_of_Gamma0_when_bm_for_bm =
                (double)getDoublePar("fraction_of_Gamma0_when_bm_for_bm",pars, AT,p_log,1.98,false);
        p_pars->fraction_of_Gamma0_when_spread =
                (double)getDoublePar("fraction_of_Gamma0_when_spread",pars, AT,p_log,.75,false);

        // set options
        std::string opt;
        // mass accretion from ISM
        opt = "method_dmdr";
        METHOD_dmdr method_dmdr;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_dmdr = METHOD_dmdr::iusingdthdR;
        }
        else{
            if(opts.at(opt) == "usingA")
                method_dmdr = METHOD_dmdr::iusingA;
            else if(opts.at(opt) == "usingdthdr")
                method_dmdr = METHOD_dmdr::iusingdthdR;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      << " given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " usingA " <<" usingdthdr " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_dmdr = method_dmdr;
//        p_pars->only_last_shell_dmdr
//            = getBoolOpt("only_last_shell_dmdr", opts, AT,p_log,false, false);


        // evolution eq
        opt = "method_dgdr";
        METHOD_dgdr method_dgdr;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_dgdr = METHOD_dgdr::iour;
        }
        else{
            if(opts.at(opt) == "our")
                method_dgdr = METHOD_dgdr::iour;
            else if(opts.at(opt) == "peer")
                method_dgdr = METHOD_dgdr::ipeer;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      << " given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " our " <<" peer " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_dgdr = method_dgdr;

        // set parameters for lateral expanding
        opt = "method_spread";
        LatSpread::METHODS method_spread;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_spread = LatSpread::iNULL;
        }
        else{
            if(opts.at(opt) == "None")
                method_spread = LatSpread::iNULL;
            else if(opts.at(opt) == "AFGPY")
                method_spread = LatSpread::iAFGPY;
            else if(opts.at(opt) == "Adi")
                method_spread = LatSpread::iAdi;
            else if(opts.at(opt) == "AA")
                method_spread = LatSpread::iAA;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " None " <<" AFGPY " << " Adi " << " AA " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_spread->setPars(a,theta_max,id->theta_core,
                          id->theta_wing, method_spread);


        // set parameters for EOS
        opt = "method_eos";
        EOSadi::METHODS method_eos;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_eos = EOSadi::iNava13;
        }
        else{
            if(opts.at(opt) == "Nava13")
                method_eos = EOSadi::iNava13;
            else if(opts.at(opt) == "Peer12")
                method_eos = EOSadi::iPeer12;
            else{
                (*p_log)(LOG_WARN,AT)<< " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized \n"
                                     << " Possible options: "
                                     << " Nava13 " << " Peer12 " << "\n"
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_eos->setPars(method_eos);

        // set Nuclear Atomic pars
        p_nuc->setPars(pars, opts);

        /// method for shock radius
        opt = "method_Rsh";
        METHOD_RSh m_method_r_sh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_r_sh = METHOD_RSh::isameAsR;
        }
        else{
            if(opts.at(opt) == "sameAsR")
                m_method_r_sh = METHOD_RSh::isameAsR;
            else if(opts.at(opt) == "useGammaSh")
                m_method_r_sh = METHOD_RSh::iuseGammaSh;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " sameAsR " << " useGammaSh " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_r_sh = m_method_r_sh;

        /// method for shock velocity
        opt = "method_GammaSh";
        METHOD_GammaSh m_method_gamma_sh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_gamma_sh = METHOD_GammaSh::isameAsGamma;
        }
        else{
            if(opts.at(opt) == "useGamma")
                m_method_gamma_sh = METHOD_GammaSh::isameAsGamma;
            else if(opts.at(opt) == "useGammaRel")
                m_method_gamma_sh = METHOD_GammaSh::iuseGammaRel;
            else if(opts.at(opt) == "useJK")
                m_method_gamma_sh = METHOD_GammaSh::iuseJK;
            else if(opts.at(opt) == "useJKwithGammaRel")
                m_method_gamma_sh = METHOD_GammaSh::iuseJKwithGammaRel;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " useGamma " << " useGammaRel " << " useJK "<< " useJKwithGammaRel " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_gamma_sh = m_method_gamma_sh;

        /// method for energy density behind shock
        opt = "method_Up";
        METHODS_Up m_method_up;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_up = METHODS_Up::iuseEint2;
        }
        else{
            if(opts.at(opt) == "useEint2")
                m_method_up = METHODS_Up::iuseEint2;
            else if(opts.at(opt) == "useGamma")
                m_method_up = METHODS_Up::iuseGamma;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " useEint2 " << " useGamma " << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_up = m_method_up;

        /// method for shock thickness
        opt = "method_Delta";
        METHOD_Delta m_method_delta;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_delta = METHOD_Delta::iuseJoh06;
        }
        else{
            if(opts.at(opt) == "useJoh06")
                m_method_delta = METHOD_Delta::iuseJoh06;
            else if(opts.at(opt) == "useVE12")
                m_method_delta = METHOD_Delta::iuseVE12;
            else if(opts.at(opt) == "None")
                m_method_delta = METHOD_Delta::iNoDelta;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized \n"
                                      << " Possible options: "
                                      << " useJoh06 " << " useVE12 " << "None" << "\n"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_Delta = m_method_delta;

        /// EATS method guverns the BW discretization (piece-wise versus adaptive)
        p_pars->m_method_eats = id->method_eats;
        p_pars->nlayers = id->nlayers;

        /// set boolean pars
        p_pars->use_dens_prof_behind_jet_for_ejecta =
                getBoolOpt("use_dens_prof_behind_jet_for_ejecta", opts, AT,p_log,false, false);

        p_pars->use_exp_rho_decay_as_floor =
                getBoolOpt("use_exp_rho_decay_as_floor", opts, AT,p_log,false, false);

        p_pars->use_flat_dens_floor =
                getBoolOpt("use_flat_dens_floor", opts, AT,p_log,false, false);

        p_pars->use_st_dens_profile =
                getBoolOpt("use_st_dens_profile", opts, AT,p_log,false, false);

        p_pars->use_bm_dens_profile =
                getBoolOpt("use_bm_dens_profile", opts, AT,p_log,false, false);

        p_pars->adiabLoss =
                getBoolOpt("use_adiabLoss", opts, AT,p_log,true, false);

        /// set sedov-taylor profile (for jet to be seen by ejecta as it moves behind)
        if (p_pars->use_st_dens_profile) {
            p_sedov->setPars(1.5, 3, 0.); // TODO this should not be here and evaluated for EVERY bw...
            p_sedov->evaluate();
        }

        /// set parameters for computing observed emission

        /// -----------  set initials and constants for the blast wave ----------------------
        size_t & ish = p_pars->ishell;
        size_t & il = p_pars->ilayer;

        p_pars->E0        = id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_pw[ilayer];
        p_pars->Ye0       = id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_pw[ilayer];
        p_pars->s0        = id->get(ish,il,EjectaID2::Q::is);//latStruct.dist_s_pw[ilayer];
        p_pars->M0        = id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_pw[ilayer];
        p_pars->R0        = id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_pw[ilayer];
        p_pars->mom0      = id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_pw[ilayer];
        p_pars->tb0       = m_tb_arr.empty() ? 0 : m_tb_arr[0];
        p_pars->theta_a   = 0.;
        p_pars->theta_b0  = ((id->method_eats) == EjectaID2::ipiecewise)
                          ? id->theta_wing : id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.m_theta_w; // theta_b0
        p_pars->ctheta0   = id->get(ish,il,EjectaID2::Q::ictheta); // TODO 0.5 * (latStruct.thetas_c_l[ilayer] + latStruct.thetas_c_h[ilayer]);
        p_pars->theta_w   = id->theta_wing;//latStruct.m_theta_w; //
        p_pars->theta_max = theta_max;
        p_pars->ncells    = ((id->method_eats) == EjectaID2::ipiecewise)
                          ? (double)id->ncells : 1.;//(double) latStruct.ncells;
        p_pars->eps_rad   = epsilon_e_rad;

        p_spread->m_theta_b0 = p_pars->theta_b0;
        p_pars->prev_x = p_pars->tb0;

        p_pars->ii_eq  = ii_eq;
#if 0
        switch (id->method_eats) {

            case EjectaID2::ipiecewise:
#if 0
                if ( latStruct.dist_E0_pw.empty() || latStruct.dist_M0_pw.empty() || latStruct.dist_Mom0_pw.empty()
                     || latStruct.dist_Ye_pw.empty() || latStruct.theta_pw.empty()){
                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
                    exit(1);
                }
                (*p_log)(LOG_INFO,AT) << " Init. [pw] "
                                      << " E0="<<latStruct.dist_E0_pw[ilayer]
                                      << " Ye="<<latStruct.dist_Ye_pw[ilayer]
                                      << " s="<<latStruct.dist_s_pw[ilayer]
                                      << " M0="<<latStruct.dist_M0_pw[ilayer]
                                      << " M0m0="<<latStruct.dist_Mom0_pw[ilayer]
                                      << " beta0="<<EQS::BetFromMom(latStruct.dist_Mom0_pw[ilayer])
                                      << " tb0="<<m_tb_arr[0]
                                      << " thetab0="<<latStruct.m_theta_w
                                      << " theta0="<<latStruct.theta_pw[ilayer]
                                      << " theta1="<<latStruct.theta_pw[ilayer]
                                      << " theta_w="<<latStruct.m_theta_w
                                      << " ii_eq="<<ii_eq
                                      << " ncells="<<latStruct.ncells
                                      << "\n";

                if (latStruct.m_theta_w > theta_max){
                    (*p_log)(LOG_ERR,AT) << " theta_b0="<<latStruct.m_theta_w<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
#endif

                p_pars->E0        = id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_pw[ilayer];
                p_pars->Ye0       = id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_pw[ilayer];
                p_pars->s0        = id->get(ish,il,EjectaID2::Q::is);//latStruct.dist_s_pw[ilayer];
                p_pars->M0        = id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_pw[ilayer];
                p_pars->R0        = id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_pw[ilayer];
                p_pars->mom0      = id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_pw[ilayer];
//                p_pars->Eint0     = id->get(ish,il,EjectaID2::Q::ieint);
                p_pars->tb0       = m_tb_arr.empty() ? 0 : m_tb_arr[0];
                p_pars->theta_a   = 0.; // theta_a
                p_pars->theta_b0  = id->theta_wing;//latStruct.m_theta_w; // theta_b0
                p_pars->ctheta0   = id->get(ish,il,EjectaID2::Q::ictheta); //TODO !! 0.5 * (latStruct.theta_pw[ilayer] + latStruct.theta_pw[ilayer]);
//        p_pars->theta_h0= theta_c_h;
//                p_pars->theta_c_l = 0.;//id->get(ish,il,EjectaID2::Q::itheta_c_l);//latStruct.theta_pw[ilayer];//theta_c_l;
//                p_pars->theta_c_h = 0.,//id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.theta_pw[ilayer];
                p_pars->theta_w   = id->theta_wing;//latStruct.m_theta_w; //
                p_pars->theta_max = theta_max;
                p_pars->ncells    = (double)id->ncells;//(double) latStruct.ncells;
                p_pars->eps_rad   = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->prev_x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;

                /// initialize non-thermal radiation module
//                p_rad->setEatsPars(pars,opts,id->m_nlayers,
//                                   id->get(ish,il,EjectaID2::Q::ictheta),
//                                   0.,0.,id->theta_wing,
//                                   getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                break;

            case EjectaID2::iadaptive:
#if 0
                if ( latStruct.dist_E0_a.empty() || latStruct.dist_M0_a.empty() || latStruct.dist_Mom0_a.empty()
                     || latStruct.dist_Ye_a.empty() || latStruct.thetas_c_h.empty()){
                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
                    exit(1);
                }
                (*p_log)(LOG_INFO,AT)<<"Init. [a] "
                                     << " E0="<<latStruct.dist_E0_a[ilayer]
                                     << " Ye="<<latStruct.dist_Ye_a[ilayer]
                                     << " s="<<latStruct.dist_s_a[ilayer]
                                     << " M0="<<latStruct.dist_M0_a[ilayer]
                                     << " G0="<<latStruct.dist_Mom0_a[ilayer]
                                     << " beta0="<<EQS::BetFromMom(latStruct.dist_Mom0_a[ilayer])
                                     << " tb0="<<m_tb_arr[0]
                                     << " thetab0="<<latStruct.thetas_c_h[ilayer]
                                     << " theta0="<<latStruct.thetas_c_l[ilayer]
                                     << " theta1="<<latStruct.thetas_c_h[ilayer]
                                     << " theta_w="<<latStruct.m_theta_w
                                     << " ii_eq="<<ii_eq
                                     << " ncells="<<1.
                                     << "\n";
                double fac = 2 * std::sin(0.5 * latStruct.thetas_c_h[ilayer]) * std::sin(0.5 * latStruct.thetas_c_h[ilayer]);
                if (latStruct.thetas_c_h[ilayer] > theta_max){
                    (*p_log)(LOG_ERR,AT) << " theta_b0="<<latStruct.thetas_c_h[ilayer]<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
#endif
                p_pars->E0      = id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_a[ilayer];
                p_pars->Ye0     = id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_a[ilayer];
                p_pars->s0      = id->get(ish,il,EjectaID2::Q::is);//latStruct.dist_s_a[ilayer];
                p_pars->M0      = id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_a[ilayer];
                p_pars->R0      = id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_a[ilayer];
                p_pars->mom0    = id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_a[ilayer];
                p_pars->tb0     = m_tb_arr.empty() ? 0 : m_tb_arr[0];
                p_pars->theta_a = 0.;
                p_pars->theta_b0= id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.thetas_c_h[ilayer];
                p_pars->ctheta0 = id->get(ish,il,EjectaID2::Q::ictheta); // TODO 0.5 * (latStruct.thetas_c_l[ilayer] + latStruct.thetas_c_h[ilayer]);
//        p_pars->theta_h0= theta_c_h;
//                p_pars->theta_c_l = id->get(ish,il,EjectaID2::Q::itheta_c_l);//latStruct.thetas_c_l[ilayer];
//                p_pars->theta_c_h = id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.thetas_c_h[ilayer];
                p_pars->theta_w = 0.; //
                p_pars->theta_max = theta_max;
                p_pars->ncells  = 1.;
                p_pars->eps_rad = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->prev_x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;

                /// initialize non-thermal radiation module
//                p_rad->setEatsPars(pars,opts,id->m_nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
//                                   id->get(ish,il,EjectaID2::Q::itheta_c_l),
//                                   id->get(ish,il,EjectaID2::Q::itheta_c_h),0.,
//                                   getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                // double E0, double M0, double Gamma0, double tb0, double theta_a, double theta_b0,
                // double theta_c_l, double theta_c_h, double theta_w, double theta_max, double epsilon_e_rad,
                // size_t ii_eq,
                // double ncells
//                bw_obj.setMagPars(latStruct.dist_E0_a[ilayer], = double E0,
//                               latStruct.dist_M0_a[ilayer], = double M0,
//                               latStruct.dist_G0_a[ilayer], = double Gamma0,
//                               t_grid[0],                   = double tb0,
//                               0.,                          = double theta_a,
//                               latStruct.thetas_c_h[ilayer],= double theta_b0,
//                               latStruct.thetas_c_l[ilayer],= double theta_c_l,
//                               latStruct.thetas_c_h[ilayer],= double theta_c_h,
//                               0.,                          = double theta_w,
//                               theta_max,                   = double theta_max,
////                             latStruct.cthetas0[ilayer],  = double epsilon_e_rad,
//                               epsilon_e_rad,               = size_t ii_eq,
//                               ii_eq,
//                               1.);
                break;
        }
#endif

        /// ----------------------- set options ------------------------------
        opt = "method_ne";
        METHOD_NE methodNe;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodNe = METHOD_NE::iuseNe;
        }
        else{
            if(opts.at(opt) == "useNe")
                methodNe = METHOD_NE::iuseNe;
            else if(opts.at(opt) == "usenprime")
                methodNe = METHOD_NE::iusenprime;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useNe " << " usenprime " << "\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_ne = methodNe;

        opt = "method_shock_vel";
        METHODS_SHOCK_VEL methodsShockVel;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
        }
        else{
            if(opts.at(opt) == "sameAsBW")
                methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
            else if(opts.at(opt) == "shockVel")
                methodsShockVel = METHODS_SHOCK_VEL::ishockVel;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << " Possible options: "
                                     << " sameAsBW " << " shockVel " << "\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->method_shock_vel = methodsShockVel;

        opt = "method_comp_mode";
        METHODS_RAD methodCompMode;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodCompMode = METHODS_RAD::iobservflux;
        }
        else{
            if(opts.at(opt) == "observFlux")
                methodCompMode = METHODS_RAD::iobservflux;
            else if(opts.at(opt) == "comovSpec")
                methodCompMode = METHODS_RAD::icomovspec;
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " observFlux " << " comovSpec " << "\n";
                exit(1);
            }
        }
        p_pars->m_method_rad = methodCompMode;

        opt = "method_single_bw_delta";
        METHOD_SINGLE_BW_DELTA method_single_bw_delta;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_single_bw_delta = METHOD_SINGLE_BW_DELTA::iconst;
        }
        else{
            if(opts.at(opt) == "frac_last") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ifrac_last;
//                p_pars->fraction_last_delta =
//                        (double)getDoublePar("fraction_last_delta",pars, AT,p_log,.75,true);
            }
            else if(opts.at(opt) == "const")
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::iconst;
            else if(opts.at(opt) == "lr")
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ilr;
            else if(opts.at(opt) == "bm") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ibm;
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
            }
            else if(opts.at(opt) == "st") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ist;
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
            }
            else if(opts.at(opt) == "bm_st") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ibm_st; // TODO add an option that smoothly connects BM and ST profiles
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
                p_pars->mom_when_bm_to_st =
                        (double)getDoublePar("mom_when_bm_to_st",pars, AT,p_log,.75,true);
            }
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " frac_last " << " bm " << " st "<< " bm "<< " bm_st "<< " bm "<< "\n";
                exit(1);
            }
        }
        p_pars->method_single_bw_delta = method_single_bw_delta;

        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, false);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e14, false);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,100, false);//pars.at("nfreq");

        p_pars->m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                                  << " freqs="<<p_pars->m_freq_arr.size() << " by radii=" << p_pars->nr << " Spec. grid="
                                  << p_pars->m_freq_arr.size() * p_pars->nr << "\n";
            p_pars->m_synch_em.resize( p_pars->m_freq_arr.size() * p_pars->nr );
            p_pars->m_synch_abs.resize( p_pars->m_freq_arr.size() *p_pars->nr );
        }

        p_pars->theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        p_pars->d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        p_pars->z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
        p_pars->theta_c_l = id->get(ish,il,EjectaID2::Q::itheta_c_l);
        p_pars->theta_c_h = id->get(ish,il,EjectaID2::Q::itheta_c_h);
        p_pars->theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);

        p_eats_fs->setEatsPars(
                pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                id->get(ish,il,EjectaID2::Q::itheta_c_l),
                id->get(ish,il,EjectaID2::Q::itheta_c_h),
                id->theta_wing,
                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));
        p_pars->p_syna->setPars(pars, opts);

        /// ---------------------------------------
//        int p_pars->opacitymode=0; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
        if (p_pars->Ye0 <= 0.2)
            p_pars->opacitymode = 2;
        else if ((p_pars->Ye0 > 0.2) or (p_pars->Ye0 < 0.3))
            p_pars->opacitymode = 1;
        else if (p_pars->Ye0 >= 0.3)
            p_pars->opacitymode = 0;
        else if (p_pars->Ye0 > 0.5)
            p_pars->opacitymode = 3;
        else{
            (*p_log)(LOG_ERR,AT) << " error \n";
            exit(1);
        }

        if(p_pars->opacitymode==0){
            p_pars->Z_eff = 24.21; /* effective nuclear weight */
            p_pars->mu_e = 2.148;
        }
        else if(p_pars->opacitymode==1){
            p_pars->Z_eff = 26.74;
            p_pars->mu_e = 2.2353;
        }
        else if(p_pars->opacitymode==2){
            p_pars->Z_eff = 53.90;
            p_pars->mu_e= 2.4622;
        }
        else if(p_pars->opacitymode==3){
            p_pars->Z_eff = 7.0;
            p_pars->mu_e = 2.0;
        }

        /// ionization fraction in ejecta
        p_pars->albd_fac = getDoublePar("albd_fac", pars, AT, p_log,0.5, false);//pars.at("freq1");
        p_pars->eps_e_w = getDoublePar("eps_e_w",pars,AT,p_log,-1,true); // electron acceleration efficiency
        p_pars->eps_mag_w = getDoublePar("eps_mag_w",pars,AT,p_log,-1,true); // magnetic field efficiency
        p_pars->epsth_w = getDoublePar("eps_th_w",pars,AT,p_log,-1,true); // initial absorption fraction
        p_pars->gamma_b_w = getDoublePar("gamma_b_w",pars,AT,p_log,-1,true); //break Lorentz factor of electron injection spectrum

        /// PWN wind initial radius
//        p_pars->pwn_Rw0 = getDoublePar("Rw0", pars, AT, p_log, -1, false); // PWN radius at t=0; [km]



//        std::cout << " ["<<bw_obj.getPars()->ishell<<", "<<bw_obj.getPars()->ilayer<<"] "
//                  <<" G0="<<bw_obj.getPars()->Gamma0
//                  <<" E0="<<bw_obj.getPars()->E0
//                  <<" M0="<<bw_obj.getPars()->M0
//                  <<" ctheta0="<<bw_obj.getPars()->ctheta0
//                  <<"\n";
    }
    // ------------------------------------------------------
    Pars *& getPars(){ return p_pars; }
    std::unique_ptr<EOSadi> & getEos(){ return p_eos; }
    std::unique_ptr<LatSpread> & getSpread(){ return p_spread; }
    std::unique_ptr<RhoISM> & getDensIsm(){ return p_dens; }
    std::unique_ptr<SedovTaylor> & getSedov(){ return p_sedov; }
    std::unique_ptr<BlandfordMcKee2> & getBM(){ return p_bm; }
    std::unique_ptr<EATS> & getFsEATS(){ return p_eats_fs; }
    std::unique_ptr<LinearRegression> & getLRforDelta(){ return p_lr_delta; }
    std::unique_ptr<LinearRegression> & getLRforVol(){ return p_lr_vol; }
    // --------------------------------------------------------
    void updateNucAtomic( double * sol, const double t ){
        p_nuc->update(
                p_pars->Ye0,
                p_pars->M0 * (double)EjectaID2::CellsInLayer(p_pars->ilayer), // m_iso ( for eps_th )
                EQS::BetFromMom( sol[ p_pars->ii_eq + SOL::QS::imom ] ),
                t,
                sol[ p_pars->ii_eq + SOL::QS::iR ],
                p_pars->s0
        );
        p_pars->kappa = p_nuc->getPars()->kappa;
        p_pars->dEnuc = p_pars->M0 * p_nuc->getPars()->eps_nuc_thermalized;
    }
    void updateCurrentBpwn( const double * sol ){
        double r_w = sol[p_pars->ii_eq + SOL::iRw];
        double u_b_pwn = 3.0*sol[p_pars->ii_eq + SOL::iWepwn]/4.0/M_PI/r_w/r_w/r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
        double b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
        p_pars->curr_b_pwn = b_pwn;
    }
    void updateEnergyInjection( double ldip, double lacc ){
        if (ldip < 0 or !std::isfinite(ldip)){
            (*p_log)(LOG_ERR,AT) << " ldip < 0 or nan; ldip="<<ldip<<"\n";
        }
        if (lacc < 0 or !std::isfinite(lacc)){
            (*p_log)(LOG_ERR,AT) << " lacc < 0 or nan; lacc="<<lacc<<"\n";
        }
        // ----------------------
        p_pars->curr_ldip = ldip;
        p_pars->curr_lacc = lacc;
    }

    /// set initial condition 'ic_arr' for this blast wave using data from PWNPars{} struct
    void setInitConditions( double * ic_arr, size_t i ) {

        /// if layer does not have energy / mass -- do not evolve it
        if ((p_pars->M0 == 0.) && (p_pars->E0 == 0.)){
//            std::cerr << AT << "\n "
//                      << "[ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
//                      <<" M0=0 and E0=0 -> Ignoring this layer.\n";
            p_pars->end_evolution = true;
            for (size_t v = 0; v < SOL::neq; ++v){
                ic_arr[i+v] = 0.;
            }
            return;
        }
        p_pars->Gamma0 = EQS::GamFromMom(p_pars->mom0);
        // ****************************************
        double beta0 = EQS::BetFromMom(p_pars->mom0);
//        p_pars->R0    = p_pars->tb0 * beta0 * CGS::c;

        p_dens->evaluateRhoDrhoDrDefault(p_pars->R0, p_pars->ctheta0);
        if (p_pars->is_within0){ //p_pars->j_i0!=123456789
            // behind a jet BW
            p_dens->m_rho_ = p_dens->m_rho_floor_val * p_dens->m_rho_def;
            p_dens->m_drhodr_ = 0.;
        }
        else{
            // outside jet BW
            p_dens->m_rho_ = p_dens->m_rho_def;
            p_dens->m_drhodr_=p_dens->m_drhodr_def;
        }

        double m_M20 = (2.0 / 3.0) * CGS::pi
                       * (std::cos(p_pars->theta_a) - std::cos(p_pars->theta_b0))
                       * p_dens->m_rho_
                       * std::pow(p_pars->R0, 3)
                       / p_pars->ncells;
        double adi0 = p_eos->getGammaAdi(p_pars->Gamma0,EQS::Beta(p_pars->Gamma0));
        double GammaSh0 = EQS::GammaSh(p_pars->Gamma0,adi0);
        // ****************************************
        if ((p_pars->mom0 <= 0.) || (!std::isfinite(p_pars->mom0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " Mom0 < 0 (Mom0=" <<p_pars->Gamma0 << ") "
                                   << "Mom0="<<p_pars->mom0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                                   << " \n";
            //std::cout << "[ Error ] " << "Gamma0 < 0 (Gamma0=" <<Gamma0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->M0 <= 0.) || (!std::isfinite(p_pars->M0))){
//            std::cout << "[ WARNING ] " << "M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            // REMOVING LOGGER
            (*p_log)(LOG_WARN, AT)  << " M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            p_pars->M0 = p_pars->E0 / (p_pars->Gamma0 * CGS::c * CGS::c);
        }
        if ((p_pars->R0 <= 1.) || (!std::isfinite(p_pars->R0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " R0 <= 0 (R0=" <<p_pars->R0 << ") " << "G0="<<p_pars->Gamma0
                                   << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<") \n"
                                   << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "R0 <= 0 (R0=" <<R0 << ")\n";
            exit(1);
        }
        if ((p_dens->m_rho_) <= 0.|| (!std::isfinite(p_dens->m_rho_))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " rho0 < 0 (rho0=" <<p_dens->m_rho_ << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                                   << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                                   << " \n";
//            std::cerr << AT  << "\n";

            //std::cout << "[ Error ] " << "rho0 < 0 (rho0=" <<rho0 << ")\n";
            exit(1);
        }
        if ((p_pars->E0 <= 0.) || (!std::isfinite(p_pars->E0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << "E0 <= 0 (E0=" <<p_pars->E0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                                   << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                                   << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "E0 < 0 (E0=" <<E0 << ")\n";
            exit(1);
        }
        if ((p_pars->Gamma0 < 1.) || (!std::isfinite(p_pars->Gamma0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " Gamma0 < 1 (Gamma0=" <<p_pars->Gamma0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                                   << " \n";
            //std::cout << "[ Error ] " << "Gamma0 < 0 (Gamma0=" <<Gamma0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->theta_b0 < 0.) || (!std::isfinite(p_pars->theta_b0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " theta_b0 < 0 (theta_b0=" <<p_pars->theta_b0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                                   << " \n";
            //std::cout << "[ Error ] " << "theta_b0 < 0 (theta_b0=" <<theta_b0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->ncells < 1 )|| (!std::isfinite(p_pars->ncells))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " ncells < 1 (ncells=" <<p_pars->ncells << ")\n"
                                   << " \n";
            //std::cout << "[ Error ] " << "ncells < 1 (ncells=" <<ncells << ")\n";
            exit(1);
        }
        if (cos(p_pars->theta_a) <= cos(p_pars->theta_b0)){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " cos(theta_a) < cos(theta_b0) (theta_a="<<p_pars->theta_a
                                   << ", theta_b0="<<p_pars->theta_b0<<")\n"
                                   << " \n";
            //std::cout << "[ Error ] " <<" cos(theta_a) < cos(theta_b0) (theta_a="<<theta_a<<", theta_b0="<<theta_b0<<")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        bool use_spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        p_pars->Rd = std::pow(3./(4.*CGS::pi) * 1./(CGS::c*CGS::c*CGS::mp) *
                              p_pars->E0/( p_dens->m_rho_ / CGS::mp * p_pars->Gamma0*p_pars->Gamma0), 1./3.);

        double x = p_pars->Rd / p_pars->R0;
        // ***************************************
        // -------------- DYNAMICS ---------------
        ic_arr[i + SOL::QS::iR]      = m_tb_arr[0] * beta0 * CGS::c;
        ic_arr[i + SOL::QS::iRsh]    = m_tb_arr[0] * EQS::Beta(GammaSh0) * CGS::c;
        ic_arr[i + SOL::QS::itt]     = EQS::init_elapsed_time(p_pars->R0, p_pars->mom0, use_spread);
        ic_arr[i + SOL::QS::itcomov] = p_pars->R0 / (beta0 * p_pars->Gamma0 * CGS::c);
        ic_arr[i + SOL::QS::iEint2]  = (p_pars->Gamma0 - 1. ) * m_M20 / p_pars->M0;  //TODO Isnt it just E0 / m_M0 ???? As M0 = E0 * cgs.c ** -2 / Gamma0
        ic_arr[i + SOL::QS::iEint2]  += p_pars->Eint0 / p_pars->M0 / CGS::c / CGS::c; // add initial internal energy
        ic_arr[i + SOL::QS::imom]    = p_pars->Gamma0 * EQS::Beta(p_pars->Gamma0);//std::log( p_pars->Gamma0 );
//        ic_arr[i + QS::iGamma]  = p_pars->Gamma0;//std::log( p_pars->Gamma0 );
        ic_arr[i + SOL::QS::itheta]  = p_pars->theta_b0;
        ic_arr[i + SOL::QS::iErad2]  = 0.0;
        ic_arr[i + SOL::QS::iEsh2]   = 0.0;
        ic_arr[i + SOL::QS::iEad2]   = 0.0;
        ic_arr[i + SOL::QS::iM2]     = m_M20 / p_pars->M0;
        // ------------ PWN -------------------
        ic_arr[i + SOL::QS::iRw] = ic_arr[i + SOL::QS::iR]; // assume that PWN wind is at R0 of ej.
        if (p_pars->curr_ldip < 0 || p_pars->curr_lacc < 1){
            (*p_log)(LOG_ERR, AT)  << "p_pars->curr_ldip = " << p_pars->curr_ldip
                                   << "p_pars->curr_lacc = "<< p_pars->curr_lacc
                                   << "\n";
            exit(1);
        }
        ic_arr[i + SOL::QS::iWenb] = p_pars->eps_e_w * p_pars->curr_ldip
                                   + p_pars->epsth_w * p_pars->curr_lacc;
        ic_arr[i + SOL::QS::iWenb] = p_pars->eps_mag_w * p_pars->curr_ldip;
        ic_arr[i + SOL::QS::iWtt] = ic_arr[i + SOL::QS::itt]; // we assume that also time is the same
        // ***************************************
        for (size_t v = 0; v < SOL::neq; ++v){
            if (!std::isfinite(ic_arr[i + v])){
                (*p_log)(LOG_ERR, AT)  << " NAN in initial data for shell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                                       << " v_n="<<SOL::vars[v]<<" val="<<ic_arr[i + v]<<" Exiting...\n";
//                std::cerr << AT  << "\n";
                exit(1);
            }
        }
//        p_spread->m_theta_b0 = p_pars->theta_b0;
//        p_pars->x = p_pars->tb0;
        // ***************************************
    }
    /// check the evolution result
    void checkEvolution(){
        /// limit the evaluation to the latest 'R' that is not 0 (before termination)
        if (p_pars->end_evolution)
            return;
        size_t nr = m_data[BW::Q::iR].size();
        size_t i_end_r = nr;
        for(size_t ir = 0; ir < nr; ++ir){
            if (m_data[BW::Q::iR][ir] == 0.) {
                i_end_r = ir;
                break;
            }
        }
//        if (i_end_r == 0){
//            (*p_log)(LOG_ERR,AT) << "[il="<< p_pars->ilayer << ", ish="<< p_pars->ishell
//                                            << "] Blastwave was not evolved: i_end_r = " << i_end_r << "\n";
//            exit(1);
//        }
        p_pars->i_end_r = i_end_r;
    }
    /// add the current solution 'sol' to the 'm_data' which is Vector of Arrays (for all variables)
    void insertSolution( const double * sol, size_t it, size_t i ) {
        if (p_pars->end_evolution)
            return;
//        double mom = sol[i+QS::imom];
//        double gam = sol[i+QS::iGamma];
        p_pars->comp_ix = it; // latest computed iteration
        m_data[BW::Q::itburst][it]   = m_tb_arr[it];
        m_data[BW::Q::iR][it]        = sol[i+SOL::QS::iR]; // TODO you do not need 'i' -- this is p_pars->ii_eq
        m_data[BW::Q::iRsh][it]      = sol[i+SOL::QS::iRsh];
        m_data[BW::Q::itt][it]       = sol[i+SOL::QS::itt];
        m_data[BW::Q::imom][it]      = sol[i+SOL::QS::imom];
        m_data[BW::Q::iGamma][it]    = EQS::GamFromMom(sol[i+SOL::QS::imom]);//sol[i+QS::iGamma];
        m_data[BW::Q::itheta][it]    = sol[i+SOL::QS::itheta];
        m_data[BW::Q::iM2][it]       = sol[i+SOL::QS::iM2];
        m_data[BW::Q::itcomov][it]   = sol[i+SOL::QS::itcomov];
        m_data[BW::Q::iEad2][it]     = sol[i+SOL::QS::iEad2];
        m_data[BW::Q::iEint2][it]    = sol[i+SOL::QS::iEint2];
        m_data[BW::Q::iEsh2][it]     = sol[i+SOL::QS::iEsh2];
        m_data[BW::Q::iErad2][it]    = sol[i+SOL::QS::iErad2];
        /// --- PWN ---
        m_data[BW::Q::i_Wtt][it]     = sol[i+SOL::QS::iWtt];
        m_data[BW::Q::i_Wmom][it]    = sol[i+SOL::QS::iWmom];
        m_data[BW::Q::i_Wepwn][it]   = sol[i+SOL::QS::iWepwn];
        m_data[BW::Q::i_Wenb][it]    = sol[i+SOL::QS::iWenb];
        m_data[BW::Q::i_WGamma][it]  = EQS::GamFromMom(sol[i+SOL::QS::iWmom]);//sol[i+QS::iGamma];
        if (sol[i+SOL::QS::iR] < 1. || m_data[BW::Q::iGamma][it] < 1. || sol[i + SOL::QS::imom] < 0) {
            (*p_log)(LOG_ERR, AT)  << "Wrong value at i=" << it << " tb=" << sol[i + SOL::QS::iR]
                                   << " iR="      << sol[i + SOL::QS::iR]
                                   << " iRsh="    << sol[i + SOL::QS::iRsh]
                                   //                       << " imom="  << m_data[Q::iGamma][it]
                                   << " iMom="    << sol[i + SOL::QS::imom]
                                   << " itheta="  << sol[i + SOL::QS::itheta]
                                   << " iM2="     << sol[i + SOL::QS::iM2]
                                   << " itcomov=" << sol[i + SOL::QS::itcomov]
                                   << " iEad2="   << sol[i + SOL::QS::iEad2]
                                   << " iEint2="  << sol[i + SOL::QS::iEint2]
                                   << " iEsh2="   << sol[i + SOL::QS::iEsh2]
                                   << " iErad2="  << sol[i + SOL::QS::iErad2]
                                   << " "
                                   << " Exiting...\n";
            std::cerr << AT  << "\n";
            exit(1);
        }
//        / toto
//        updateCurrentBpwn(sol);
    }
    /// Mass and energy are evolved in units of M0 and M0c2 respectively
    void applyUnits( double * sol, size_t i ) {
        sol[i + SOL::QS::iM2]    *= p_pars->M0;
        sol[i + SOL::QS::iEint2] *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEad2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEsh2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEad2]  *= (p_pars->M0 * CGS::c * CGS::c);
    }
    /// check if to terminate the evolution
    bool isToTerminate( double * sol, size_t i ) {
        double mom = sol[i + SOL::QS::imom];
//        double igamma = sol[i+QS::iGamma];//EQS::GamFromMom(mom);
        double igamma = EQS::GamFromMom(mom);
        double ibeta = EQS::Beta(igamma);//EQS::BetFromMom(mom);
        double iEint2 = sol[i + SOL::QS::iEint2];
        bool do_terminate = false;
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
//            double beta = EQS::Beta(igamma);
            if ((iEint2 <= 0.)||(!std::isfinite(iEint2))||(mom < 0)||(ibeta < 0)){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: Eint2 < 0 or NAN ("<<iEint2<<")\n";
            }
            if ((std::abs(ibeta) < p_pars->min_beta_terminate)||(!std::isfinite(ibeta))){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: beta < beta_min or NAN ("<<ibeta<<") with beta_min= "
                        << p_pars->min_beta_terminate<<"\n";
            }
            if (do_terminate) {
                (*p_log)(LOG_WARN, AT) << " TERMINATION Previous iterations: \n"
                                       << " E0=" << string_format("%.9e", p_pars->E0)
                                       << " Gamma0=" << string_format("%.9e", p_pars->Gamma0)
                                       << " M0=" << string_format("%.9e", p_pars->M0)
                                       << " \n";
                (*p_log)(LOG_WARN, AT) << " i=Current" << " Eint2=" << string_format("%.9e", sol[i + SOL::QS::iEint2])
                                       << " Gamma=" << string_format("%.9e", igamma)
                                       << " beta=" << string_format("%.9f", ibeta)
                                       << " M2=" << string_format("%.9e", sol[i + SOL::QS::iM2])
                                       << " rho=" << string_format("%.9e", p_dens->m_rho_)
                                       << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix << " Eint2="
                                       << string_format("%.9e", getVal(BW::Q::iEint2, p_pars->comp_ix))
                                       << " Gamma=" << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix))
                                       << " beta=" << string_format("%.9f", getVal(BW::Q::ibeta, p_pars->comp_ix))
                                       << " M2=" << string_format("%.9e", getVal(BW::Q::iM2, p_pars->comp_ix))
                                       << " rho=" << string_format("%.9e", getVal(BW::Q::irho, p_pars->comp_ix))
                                       << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix - 1 << " Eint2="
                                       << string_format("%.9e", getVal(BW::Q::iEint2, p_pars->comp_ix - 1))
                                       << " Gamma=" << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix - 1))
                                       << " beta=" << string_format("%.9f", getVal(BW::Q::ibeta, p_pars->comp_ix - 1))
                                       << " M2=" << string_format("%.9e", getVal(BW::Q::iM2, p_pars->comp_ix - 1))
                                       << " rho=" << string_format("%.9e", getVal(BW::Q::irho, p_pars->comp_ix - 1))
                                       << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix - 2 << " Eint2="
                                       << string_format("%.9e", getVal(BW::Q::iEint2, p_pars->comp_ix - 2))
                                       << " Gamma=" << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix - 2))
                                       << " beta=" << string_format("%.9f", getVal(BW::Q::ibeta, p_pars->comp_ix - 2))
                                       << " M2=" << string_format("%.9e", getVal(BW::Q::iM2, p_pars->comp_ix - 2))
                                       << " rho=" << string_format("%.9e", getVal(BW::Q::irho, p_pars->comp_ix - 2))
                                       << " \n";
            }
        }
        return do_terminate;
    }
    /// check if to terminate the evolution
    bool isToStopLateralExpansion( double * sol, size_t i ) {
        double itheta = sol[i + SOL::QS::itheta];
        double iEint2 = sol[i + SOL::QS::iEint2];
        bool do_terminate = false;
        if (!p_pars->end_spreading) {
            if (itheta > 0.99 * p_pars->theta_max){
//                std::cerr << AT << " \n"
//                          << " Stopping Jet Lateral Expansion  ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer
//                          << " with Gamma0=" << p_pars->Gamma0 << " [iteration=" << p_pars->comp_ix << "] REASON :: theta > 0.99 theta_max ("<<string_format("%.3e", p_pars->theta_max)<<") \n";
//                std::cerr << " Last Gammas: " << "\n"
//                          << " Gamma[" << p_pars->comp_ix << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix)) << "\n"
//                          << " theta[" << p_pars->comp_ix << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix))
//                          << " Gamma[" << p_pars->comp_ix - 1 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 1)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 1 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 1))
//                          << " Gamma[" << p_pars->comp_ix - 2 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 2)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 2 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 2))
//                          << " Gamma[" << p_pars->comp_ix - 3 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 3)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 3 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 3))
//                          << "\n";
                do_terminate = true;
            }
        }
        return do_terminate;
    }
    /// check if the solution 'makes sense'
    bool isSolutionOk( double * sol, size_t i ) {

        bool no_issues = true;
        double mom = sol[i + SOL::QS::imom];
        double mom0 = p_pars->Gamma0 * EQS::Beta(p_pars->Gamma0);
//        double igamma = sol[i+QS::iGamma];
        double igamma = EQS::GamFromMom(mom);
//        double ibeta = EQS::Beta(igamma);//EQS::BetFromMom(mom);
        double ibeta = EQS::BetFromMom(mom);
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
            double beta   = EQS::Beta(igamma);

            if ((igamma < 1.) || (mom < 0)) {
                // REMOVING LOGGER

                //            std::cout << "theta:" << getData(Q::iR)<<"\n";
                (*p_log)(LOG_ERR,AT)  << AT << " \n"
                                      << " Gamma < 1. ishell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                                      << " with Gamma0="<<p_pars->Gamma0 << " [iteration="<<p_pars->comp_ix<<"] \n";
                (*p_log)(LOG_ERR,AT) << " Last Gammas: " << "\n";
                (*p_log)(LOG_ERR,AT)<< " Gamma[" << p_pars->comp_ix << "]="
                                    << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix << "]="
                                    << string_format("%.9e", getVal(BW::Q::ibeta, p_pars->comp_ix))
                                    << " Gamma[" << p_pars->comp_ix - 1 << "]="
                                    << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix - 1)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 1 << "]="
                                    << string_format("%.9e", getVal(BW::Q::ibeta, p_pars->comp_ix - 1))
                                    << " Gamma[" << p_pars->comp_ix - 2 << "]="
                                    << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix - 2)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 2 << "]="
                                    << string_format("%.9e", getVal(BW::Q::ibeta, p_pars->comp_ix - 2))
                                    << " Gamma[" << p_pars->comp_ix - 3 << "]="
                                    << string_format("%.9e", getVal(BW::Q::iGamma, p_pars->comp_ix - 3)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 3 << "]="
                                    << string_format("%.9e", getVal(BW::Q::ibeta, p_pars->comp_ix - 3))
                                    << "\n";
//                std::cerr << "Gamma:" << getData(Q::iGamma) << "\n";
//                std::cerr << "R:" << getData(Q::iR) << "\n";
                (*p_log)(LOG_ERR,AT)  << " Gamma cannot be less than 1. or too large. "
                                         "Found: Gamma(" << igamma << "). Gamma0="<<p_pars->Gamma0<<" \n";
                no_issues = false;
            }
//            if (igamma > p_pars->Gamma0 * 2.) {
//                // REMOVING LOGGER
//                (*p_log)(LOG_ERR,AT)  << " too large value of Gamma(" << igamma
//                           << ") > 2*Gamma0" << p_pars->Gamma0 * 2 << ")\n";
//                no_issues = false;
//            }
            if (sol[i + SOL::QS::itheta] < p_pars->theta_b0 * 0.99) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical decrease in theta(" << sol[i + SOL::QS::itheta]
                                     << ") < theta_b0(" << p_pars->theta_b0 << ")\n";
                no_issues = false;
            }
            if (sol[i + SOL::QS::itt] <= 0.) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical value of observer time (on axis) itt("
                                     << sol[i + SOL::QS::itt] << ") < 0 \n";
                no_issues = false;
            }
        }
        return no_issues;
    }
    /// right hand side of the blast wave evolution equation (ODE)
    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ) {
        // ****************************************
        double R      = Y[i+ SOL::QS::iR];
        double Rsh    = Y[i+ SOL::QS::iRsh];
//        double tcomov = Y[i+Q::itcomov];
        double mom    = Y[i+ SOL::QS::imom];
//        double Gamma  = Y[i+QS::iGamma];
        double Eint2  = Y[i+ SOL::QS::iEint2];
        double theta  = Y[i+ SOL::QS::itheta];
        double M2     = Y[i+ SOL::QS::iM2];
        // ****************************************
        if (mom < 0) {
            (*p_log)(LOG_ERR, AT) << " negative momentum\n";
            exit(1);
        }
        double Gamma = EQS::GamFromMom(mom);
//        if (Gamma <= 1.)
//            Gamma = EQS::GamFromMom(mom);
//            Gamma = 1.0000000001;
//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }
        double beta   = EQS::Beta(Gamma);

        // ****************************************
//        auto *_pars = (struct RHS_pars *) rhs_pars;
        double ctheta_ = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);//p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
        p_dens->evaluateRhoDrhoDrDefault(R, ctheta_);
        double rho = p_dens->m_rho_def / p_pars->M0;
        double drhodr = p_dens->m_drhodr_def / p_pars->M0;

//        dlnrho1dR /= p_pars->M0;
//        double _rhoi = p_pars->rho/p_pars->M0;
        // ****************************************
        double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...

        double gammaAdi  = p_eos->getGammaAdi(Gamma, beta);
        double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
        double dRshdt = EQS::Beta(GammaSh) * CGS::c;

        double dthetadr = 0.0;
        if ( (theta < 2*p_pars->theta_max)
             && (R > p_pars->Rd)
             && (Gamma < std::max(2., p_pars->Gamma0*p_pars->fraction_of_Gamma0_when_spread))
             && (!p_pars->end_spreading) ) { // &&(Gamma < p_pars->Gamma0*.95)&&
            dthetadr = p_spread->getDthetaDr(Gamma, R, gammaAdi, theta);
//            if (dthetadr > 0.){
//                int x = 1;
//            }
        }


        double dM2dR = 0;
        switch (p_pars->m_method_dmdr) {

            case iusingA:
                dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
                break;
            case iusingdthdR:
                dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
                break;
        }
        if (dM2dR < 0.){
            (*p_log)(LOG_ERR, AT) << " dM/dR < 0\n";
            exit(1);
        }

        // --- dGammadR ---
        double dGammadR = 0.,GammaEff=0.,dGammaEffdGamma=0., num=0., denum=0.;
        switch (p_pars->m_method_dgdr) {

            case iour:
                GammaEff = get_GammaEff(Gamma, gammaAdi); // gammaAdi) #(gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
                dGammaEffdGamma = get_dGammaEffdGamma(Gamma, gammaAdi); // 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.
                num = ((Gamma - 1.) * (GammaEff + 1.) * dM2dR - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR / M2 - drhodr / rho));
                denum = ((1. + M2) + Eint2 * dGammaEffdGamma + GammaEff * (gammaAdi - 1.) * Eint2 / Gamma);
                dGammadR = -num/denum;
                break;
            case ipeer:
                dGammadR = EQS::dgdr(1, Gamma, beta, M2, gammaAdi, dM2dR);
                break;
        }
//        double dmomdR = dGammadR * (Gamma / std::sqrt(Gamma * Gamma - 1.0));
        double dmomdR = dGammadR / beta;


//        if (dGammadR > 0){
//            if ( Gamma > 0.95 * p_pars->Gamma0 ){
//                dGammadR = 0.;
//            }
//        }
        // -- Energies --
        double dEsh2dR  = (Gamma - 1.0) * dM2dR; // Shocked energy;
        double dlnV2dR  = dM2dR / M2 - drhodr / rho - dGammadR / Gamma;
        double dEad2dR  = 0.0;
        if ( p_pars->adiabLoss )
            dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;
        // -- Radiative losses
        double dErad2dR = p_pars->eps_rad * dEsh2dR;
        // -- Energy equation
        double dEint2dR = dEsh2dR + dEad2dR - dErad2dR; // / (m_pars.M0 * c ** 2)
        double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
#if 0
        double dttdr = 0.;
        if (p_spread->m_method != LatSpread::METHODS::iNULL)
            dttdr = 1. / (CGS::c * beta) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
        else
            dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
#endif
        bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        double dttdr = EQS::evalElapsedTime(R,mom,dthetadr,spread);
        // ****************************************

        if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives. Exiting..." << "\n";
            exit(1);
        }
        // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
        if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0.
            || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
            (*p_log)(LOG_ERR,AT) << " nan in derivatives. "
                                 << " dRdt="<<dRdt<<"\n"
                                 << " dM2dR="<<dM2dR<<"\n"
                                 << " dthetadr=" << dthetadr << "\n"
                                 << " dGammadR=" << dGammadR << "\n"
                                 << " dthetadr=" << dRdt << "\n"
                                 << " dEsh2dR=" << dEsh2dR << "\n"
                                 << " dlnV2dR=" << dlnV2dR << "\n"
                                 << " dEad2dR=" << dEad2dR << "\n"
                                 << " dErad2dR=" << dErad2dR << "\n"
                                 << " dEint2dR=" << dEint2dR << "\n"
                                 << " dttdr=" << dttdr << "\n"
                                 << "  \n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        if (std::abs(dRdt * dGammadR) > p_pars->Gamma0){
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives. "
                                  << "i="<< i << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"] " << "\n"
                                  << " dGamma/dr > Gamma0 "
                                  <<  "dG/dr="<< dRdt * dGammadR
                                  << " Gamma0=" <<p_pars->Gamma0 << "\n"
                                  << " E0=" <<p_pars->E0 << "\n"
                                  << " M0=" <<p_pars->M0 << "\n"
                                  << " R=" <<R << " Rsh=" <<Rsh<< " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
                                  << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
            (*p_log)(LOG_ERR,AT)  << " maybe timestep is too large \n Exiting...";
//            std::cerr << AT << "\n";
            exit(1);
        }
        p_pars->prev_x = x; // update
        out_Y[i + SOL::QS::iR] = dRdt;//1.0 / beta / CGS::c;
        out_Y[i + SOL::QS::iRsh] = dRshdt;//1.0 / beta / CGS::c;
        out_Y[i + SOL::QS::itt] = dRdt * dttdr;
        out_Y[i + SOL::QS::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i + QS::iGamma] = 0.;//dRdt * dGammadR;//dRdt * dGammadR/Gamma;
        out_Y[i + SOL::QS::imom] = dRdt * dmomdR;//dRdt * dGammadR/Gamma;
        out_Y[i + SOL::QS::iEint2] = dRdt * dEint2dR;
        out_Y[i + SOL::QS::itheta] = dRdt * dthetadr;
        out_Y[i + SOL::QS::iErad2] = dRdt * dErad2dR;
        out_Y[i + SOL::QS::iEsh2] = dRdt * dEsh2dR;
        out_Y[i + SOL::QS::iEad2] = dRdt * dEad2dR;
        out_Y[i + SOL::QS::iM2] = dRdt * dM2dR;
//        if (beta>1e-5)  p_pars->end_evolution = true;
        // ****************************************
//        if (mom+out_Y[i + QS::imom])
    }
    /// RHS with ODEs for dGammadR modified for pre-accelerated ISM
    void evaluateRhsDens( double * out_Y, size_t i, double x, double const * Y ) {
//        double Gamma = Y[i+QS::iGamma];//EQS::GamFromMom(Y[i+QS::imom]);
        double mom = Y[i+ SOL::QS::imom];
        // ****************************************
        double R      = Y[i+ SOL::QS::iR];
        double Rsh    = Y[i+ SOL::QS::iRsh];
//        double tcomov = Y[i+Q::itcomov];
//        double Gamma  = std::exp(Y[i+QS::ilnGamma]);
        double Eint2  = Y[i+ SOL::QS::iEint2];
        double theta  = Y[i+ SOL::QS::itheta];
        double M2     = Y[i+ SOL::QS::iM2];
        // *****************************************
        double Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
//        double beta = EQS::Beta(Y[i+QS::iGamma]);//EQS::BetFromMom(Y[i+QS::imom]);
        double beta = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
        if (mom < 0){
            (*p_log)(LOG_ERR,AT) << "Error\n";
            mom = 1e-5;
            Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
            beta = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
        }
        // ****************************************
        if (!std::isfinite(R) || !std::isfinite(Gamma) || M2 < 0.
            || !std::isfinite(M2) || !std::isfinite(Eint2) || Eint2<0) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " << "\n"
                                  << " R="<<R<<"\n"
                                  << " M2="<<M2<<"\n"
                                  << " Gamma=" << Gamma << "\n"
                                  << " Mom=" << mom << "\n"
                                  << " Eint2=" << Eint2
                                  << " \n";
            exit(1);
        }
//        if (Gamma <= 1.) { // TODO to be removed
//            Gamma = 1.0001;
//            (*p_log)(LOG_ERR, AT) << " Gamma < 1 in RHS for kN Ejecta\n";
//        }
//        if (Gamma > p_pars->Gamma0){ // TODO to be removed
//            Gamma = p_pars->Gamma0;
//            (*p_log)(LOG_ERR, AT) << " Gamma > Gamma0 in RHS for kN Ejecta\n";
//        }
//        if (mom < 0.){
//            (*p_log)(LOG_ERR, AT) << " mom < 0 = "<< mom << " in kN RHS dynamics\n";
//            exit(1);
//        }



//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }
//        double beta   = EQS::Beta(Gamma);


        // ****************************************
        // Get ISM density and its velocity
//        double rho, m_drhodr;
        double ctheta_ = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);// = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//        p_dens->getDrhoDr( rho, m_drhodr, R,ctheta );
//        rho /= p_pars->M0;
//        m_drhodr /= p_pars->M0;
//        p_dens->evaluate(R, ctheta);
        if (p_dens->m_GammaRel < 0){
            (*p_log)(LOG_ERR,AT) << AT << " GammaRel is not set or incorrect \n";
            exit(1);
        }
        double rho = p_dens->m_rho_ / p_pars->M0;
        double drhodr = p_dens->m_drhodr_ / p_pars->M0;
        double GammaRho = p_dens->m_GammaRho;
        double GammaRel = p_dens->m_GammaRel; // Gamma
        double dGammaRelDGamma = p_dens->m_dGammaReldGamma; // 1.
        double dGammaRhodR = p_dens->m_dGammaRhodR; // 0.
        double cs_cbm = p_dens->m_CS_CBM;


        if (GammaRho < 1.) {
            (*p_log)(LOG_ERR,AT)  << "GammaRho=" << GammaRho << "\n" << " Exiting...";
//            std::cerr << AT<< "\n";
            exit(1);
        }
        if (GammaRel < 1.) {
            (*p_log)(LOG_ERR,AT)  << "GammaRel=" << GammaRel << "\n"; GammaRel = Gamma;
//            std::cerr << AT<< "\n";
            exit(1);
        }
        if (!std::isfinite(dGammaRelDGamma)) { (*p_log)(LOG_ERR,AT) << "dGammaRelDGamma=" << dGammaRelDGamma << "\n"; exit(1); }
        if (!std::isfinite(dGammaRhodR)) { (*p_log)(LOG_ERR,AT)  << "dlnGammaCBMdR=" << dGammaRhodR << "\n"; exit(1); }
        // ****************************************
        double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...
        double Gamma_ = GammaRel;
        double beta_ = EQS::Beta(Gamma_);
        double gammaAdi  = p_eos->getGammaAdi(Gamma_, beta_);//p_eos->getGammaAdi(Gamma, beta);

        double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
        double dRshdt = EQS::Beta(GammaSh) * CGS::c;

        double dthetadr = 0.0;
//        if (theta < (p_pars->theta_max * 0.999999) ) {
//            dthetadr = p_spread->getDthetaDr(Gamma, R, gammaAdi, theta);
//        }
//        double dM2dR     = EQS::dmdr(Gamma, R, p_pars->theta_a,
//                                     theta, rho, p_spread->m_aa) / p_pars->ncells;

        double dM2dR = 0.;

        switch (p_pars->m_method_dmdr) {

            case iusingA:
                /// Old equation from Gavin P. Lamb work
                dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
                break;
            case iusingdthdR:
                /// Equation motivated by Gavin P. Lamb paper
                dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
                break;
            case iNodmdr:
                /// if no mass needs to be swept up
                break;
        }
//        dM2dR*=0;
        if (dM2dR < 0.){
            (*p_log)(LOG_ERR,AT) << " dMdR < 0 in RHS dyn for kN ejecta\n";
            exit(1);
        }

        // --- Energy injection --- ||
        double xi_inj = 1.; // TODO put in parfile
        double dEindt = p_pars->dEinjdt;
        if (!std::isfinite( dEindt) || dEindt < 0){
            (*p_log)(LOG_ERR,AT) << " dEindt = "<<dEindt << "\n";
            exit(1);
        }

        double dEinjdt = dEindt / (p_pars->M0 * CGS::c * CGS::c) / p_pars->ncells;
        double dEinjdR = dEinjdt / dRdt;
        double theta_ej = 0.; // assume that ejecta is alinged with magnetar emission?..
        double Doppler = Gamma / (1. - beta * std::cos(theta_ej));
        double dEinjdR_dop = dEinjdR * Doppler;
        double dEingdR_abs = dEinjdR;// * ( 1. - std::exp(-1.*p_pars->dtau) ) * std::exp(-1.*p_pars->tau_to0);
        double dEingdR_abs_dop = dEingdR_abs / Doppler / Doppler;

        double dEnuc = p_pars->dEnuc;//->getPars()->eps_nuc_thermalized;
        double dEnuc_norm = dEnuc / (p_pars->M0 * CGS::c * CGS::c);
        double dEnucdR = dEnuc_norm / dRdt;

        double dElum = p_pars->dElum;
        double dElum_norm = dElum / (p_pars->M0 * CGS::c * CGS::c);
        double dElumdR = dElum_norm / dRdt;
        if ((std::abs(dElumdR) > std::abs(dEnucdR*1.5)) and (dEnuc > 0.))
            dElumdR = dEnucdR;

        // --- dGammadR ---
        double dGammadR = 0., GammaEff=0.,dGammaEffdGamma=0.,num1=0.,num2=0.,num3=0.,denum1=0.,denum2=0.,denom3=0.;
//        double _tmp = (1. - GammaEff / Gamma * xi_inj) * dEingdR_abs;
        switch (p_pars->m_method_dgdr) {

            case iour:
                GammaEff = get_GammaEff(Gamma_, gammaAdi); // TODO check if GammaRel should be used for this!!!
                dGammaEffdGamma = get_dGammaEffdGamma(Gamma_, gammaAdi);
                num1 = (Gamma - GammaRho + GammaEff * (GammaRel - 1.)) * dM2dR;
                num2 = - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR/M2 - drhodr/rho - dGammaRhodR / GammaRho); // - 3.*Eint2/R
                num3 = - (1. + GammaEff / Gamma * xi_inj) * dEingdR_abs;
                denum1 = (1.+M2);
                denum2 = Eint2 * dGammaEffdGamma;
                denom3 = GammaEff * (gammaAdi - 1.) * Eint2 * dGammaRelDGamma / GammaRel;
                dGammadR = -1. * (num1 + num2 + num3) / (denum1+denum2+denom3);
                break;
            case ipeer:
                dGammadR = EQS::dgdr(1, Gamma, beta, M2, gammaAdi, dM2dR);
                break;
        }



//        if ((dGammadR > 0) && (Gamma > p_pars->Gamma0*0.99)){
//            std::cerr << AT << " \n"
//                      << "dGammadR>0 ("<<dGammadR<<") and Gamma="<<Gamma<<" > Gamma0="<<p_pars->Gamma0<<"\n";
//            std::cerr << " p_dens->m_GammaRho="<<p_dens->m_GammaRho<<"\n"
//                      << " p_dens->m_GammaRel="<<p_dens->m_GammaRel<<"\n"
//                      << " p_dens->m_dGammaReldGamma="<<p_dens->m_dGammaReldGamma<<"\n"
//                      << " p_dens->m_dGammaRhodR="<<p_dens->m_dGammaRhodR<<"\n";
//            dGammadR = 0.;
//        }
        //        double dmomdR = dGammadR * (Gamma / std::sqrt(Gamma * Gamma - 1.0));
        double dmomdR = dGammadR / beta;


//        drhodr = 0.; // TODO why was this set?! This overrides the derivative
        // -- Energies --

//        double dlnV2dR  = dM2dR / M2 - m_drhodr - dGammadR / Gamma;
        double dlnV2dR  = dM2dR / M2 - drhodr / rho - (1./GammaRel)*dGammaRelDGamma*dGammadR + dGammaRhodR / GammaRho;// + 3./R;
        double dlnVdR = (3./R) - (1/Gamma * dGammadR);
        double dEad2dR  = 0.0;
        if ( p_pars->adiabLoss )
            dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;
//        double dx = dRdt * (x - p_pars->prev_x);

//        double mom_ = mom + dmomdR * dRdt * (x - p_pars->x);
//        if (mom_ < 0){
//            int x = 1;
//        }

        // --- shock energy
        double dEsh2dR = (GammaRel - 1.0) * dM2dR;
//        if (cs_cbm > 100){
//            int y = 1;
//        }
        if ((cs_cbm > EQS::Beta(GammaRel)&&(p_pars->comp_ix > p_pars->nr*1e-2))){ // TODO Subsonic flow -- no shock
            dEsh2dR *= 1e-10;
        }
//        double dEsh2dR = 0.;

//        double dEsh2dR  = (GammaRel - 1.0) * dM2dR; // Shocked energy;
        // -- Radiative losses
        double dErad2dR = p_pars->eps_rad * dEsh2dR;
        // -- Energy equation
        double dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEnucdR - dElumdR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)
//        double _x = dEint2dR/Eint2;
//        double Eint_ = Eint2 + dEint2dR * dRdt * (x-p_pars->prev_x);
//        if (Eint_ < 0){
//            dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)
//            double Eint__ = Eint2 + dEint2dR * dRdt * (x-p_pars->prev_x);
//            if (Eint__ < 0){
//                int z = 1;
//            }
//        }

        double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
//        double dttdr;
//        if (p_spread->m_method != LatSpread::METHODS::iNULL)
//            dttdr = 1. / (CGS::c * beta) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
//        else
//            dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
//        dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        double dttdr = EQS::evalElapsedTime(R,mom,dthetadr,spread);
        // ****************************************
//        if (!std::isfinite(Gamma) ||
//        !std::isfinite(beta) //||
////        !std::isfinite(Eint2) ||
////        !std::isfinite(theta) ||
////        !std::isfinite(M2)
////        || Eint2 < 0.
////        || M2 < 0.
//          ) {
//            std::cerr << AT << " wrong data at the beginning of RHS t_b=" << x << "\n";
//            std::cerr << " Gamma="<<Gamma<<"\n"
//                      << " theta="<<theta<<"\n"
//                      << " M2="<<M2<<"\n"
//                      << " Eint2="<<Eint2<<"\n"
//                      << " p_dens->m_GammaRho="<<p_dens->m_GammaRho<<"\n"
//                      << " p_dens->Gamma_rel="<<p_dens->m_GammaRel<<"\n"
//                      << " p_dens->m_dGammaRhodR="<<p_dens->m_dGammaRhodR<<"\n"
//                      << " p_dens->m_dGammaReldGamma="<<p_dens->m_dGammaReldGamma<<"\n"
//                      << " p_dens->m_CS_CBM="<<p_dens->m_CS_CBM<<"\n"
//                      << " p_dens->m_P_cbm="<<p_dens->m_P_cbm<<"\n";
//            exit(1);
//        }
        if (mom < 0.){
            (*p_log)(LOG_ERR, AT) << " mom < 0 = "<< mom << " in kN RHS dynamics\n";
//            exit(1);
        }
        if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0.
            || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " //<< "\n"
                                  << " dRdt="<<dRdt//<<"\n"
                                  << " dM2dR="<<dM2dR//<<"\n"
                                  << " dthetadr=" << dthetadr// << "\n"
                                  << " dGammadR=" << dGammadR //<< "\n"
                                  << " dthetadr=" << dRdt //<< "\n"
                                  << " dEsh2dR=" << dEsh2dR //<< "\n"
                                  << " dlnV2dR=" << dlnV2dR //<< "\n"
                                  << " dEad2dR=" << dEad2dR //<< "\n"
                                  << " dErad2dR=" << dErad2dR //<< "\n"
                                  << " dEint2dR=" << dEint2dR //<< "\n"
                                  << " dttdr=" << dttdr// << "\n"
                                  << " dGammaRelDGamma=" << dGammaRelDGamma //<< "\n"
                                  << " dGammaRhodR=" << dGammaRhodR //<< "\n"
                                  << " drhodr=" << drhodr //<< "\n"
                                  << " Gamma="<<Gamma//<< "\n"
                                  << " beta="<<beta//<< "\n"
                                  << " theta"<<theta//<< "\n"
                                  << " M2="<<M2//<< "\n"
                                  << " Eint2="<<Eint2//<< "\n"
                                  << " GammaRel"<<GammaRel//<< "\n"
                                  << " GammaRho"<<GammaRho//<< "\n"
                                  << " \n";
            exit(1);
        }
        // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
//        if (dGammadR > 0){
//            std::cout << " Gamma="<<Gamma<< " Gamma0="<<p_pars->Gamma0<<" Gamma/Gamma0="<<Gamma/p_pars->Gamma0<<"\n";
//        }
//        p_pars->x = x; // update
//        if (Gamma > 2.*p_pars->Gamma0){
//            int x = 1;
//        }
//        if ( std::abs(dRdt * dGammadR ) > 2. * p_pars->Gamma0 )+{
//            dGammadR = 0.;
//        }
//        if ( std::abs(dRdt * dGammadR ) > p_pars->Gamma0 ){
//            std::cerr << AT
//                      << "i="<< i << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"] " << "\n"
//                      << " dGamma/dr > Gamma0 "
//                      <<  "dG/dr="<< dRdt * dGammadR
//                      << " Gamma0=" <<p_pars->Gamma0 << "\n"
//                      << " theta_b0=" <<p_pars->theta_b0 << "\n"
//                      << " GammaRho=" <<GammaRho << "\n"
//                      << " GammaRel=" <<GammaRel << "\n"
//                      << " dGammaRelDGamma=" <<dGammaRelDGamma << "\n"
//                      << " dGammaRhodR=" <<dGammaRhodR << "\n"
//                      << " R=" <<R << " Rsh=" <<Rsh<< " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
//                      << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
//            std::cerr << " maybe timestep is too large \n";
////            exit(1);
//            dGammadR = 0.;
//            out_Y[i+QS::iR]      = dRdt;
//            out_Y[i+QS::iRsh]    = dRshdt;
//            out_Y[i+QS::itt]     = dRdt * dttdr;
//            out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+QS::iGamma]  = 0 * dGammadR;
//            out_Y[i+QS::iEint2]  = 0 * dEint2dR;
//            out_Y[i+QS::itheta]  = 0 * dthetadr;
//            out_Y[i+QS::iErad2]  = 0 * dErad2dR;
//            out_Y[i+QS::iEsh2]   = 0 * dEsh2dR;
//            out_Y[i+QS::iEad2]   = 0 * dEad2dR;
//            out_Y[i+QS::iM2]     = 0 * dM2dR;
////            Gamma = p_pars->Gamma0;
//        }
//        else{
//            out_Y[i+QS::iR]      = dRdt;//1.0 / beta / CGS::c;
//            out_Y[i+QS::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
//            out_Y[i+QS::itt]     = dRdt * dttdr;
//            out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+QS::iGamma]  = dRdt * dGammadR;
//            out_Y[i+QS::iEint2]  = dRdt * dEint2dR;
//            out_Y[i+QS::itheta]  = dRdt * dthetadr;
//            out_Y[i+QS::iErad2]  = dRdt * dErad2dR;
//            out_Y[i+QS::iEsh2]   = dRdt * dEsh2dR;
//            out_Y[i+QS::iEad2]   = dRdt * dEad2dR;
//            out_Y[i+QS::iM2]     = dRdt * dM2dR;
//        }
        p_pars->prev_x = x;
        out_Y[i+ SOL::QS::iR]      = dRdt;//1.0 / beta / CGS::c;
        out_Y[i+ SOL::QS::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
        out_Y[i+ SOL::QS::itt]     = dRdt * dttdr;
        out_Y[i+ SOL::QS::itcomov] = dRdt * dtcomov_dR;
        out_Y[i+ SOL::QS::imom]    = dRdt * dmomdR; //dRdt * dGammadR / Gamma;
//        out_Y[i+QS::iGamma]  = 0.;//dRdt * dGammadR; //dRdt * dGammadR / Gamma;
        out_Y[i+ SOL::QS::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+QS::iEinj]   = dRdt * dEingdR_abs;
        out_Y[i+ SOL::QS::itheta]  = dRdt * dthetadr;
        out_Y[i+ SOL::QS::iErad2]  = dRdt * dErad2dR;
        out_Y[i+ SOL::QS::iEsh2]   = dRdt * dEsh2dR;
        out_Y[i+ SOL::QS::iEad2]   = dRdt * dEad2dR;
        out_Y[i+ SOL::QS::iM2]     = dRdt * dM2dR;
//        if (Gamma)
//        out_Y[i+QS::iR]      = dRdt;//1.0 / beta / CGS::c;
//        out_Y[i+QS::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
//        out_Y[i+QS::itt]     = dRdt * dttdr;
//        out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i+QS::iGamma]  = dRdt * dGammadR;
//        out_Y[i+QS::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+QS::itheta]  = dRdt * dthetadr;
//        out_Y[i+QS::iErad2]  = dRdt * dErad2dR;
//        out_Y[i+QS::iEsh2]   = dRdt * dEsh2dR;
//        out_Y[i+QS::iEad2]   = dRdt * dEad2dR;
//        out_Y[i+QS::iM2]     = dRdt * dM2dR;
        // ****************************************
//        if (beta<1e-5)  p_pars->end_evolution = true;
    }
    /// RHS with ODEs for dGammadR modified for pre-accelerated ISM

    void evaluateRhsDensPWN( double * out_Y, size_t i, double x, double const * Y ) {
        // *************| Extract Values |***************************
        double mom    = Y[i+ SOL::QS::imom];
        double R      = Y[i+ SOL::QS::iR];
        double Rsh    = Y[i+ SOL::QS::iRsh];
        double tcomov = Y[i+ SOL::QS::itcomov];
        double Eint2  = Y[i+ SOL::QS::iEint2];
        double theta  = Y[i+ SOL::QS::itheta];
        double M2     = Y[i+ SOL::QS::iM2];
        double Gamma  = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
        double beta   = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
        /// ---- PWN ---
        double r_w   = Y[i+ SOL::iRw];
        double e_nb  = Y[i+ SOL::iWenb];
        double e_pwn = Y[i+ SOL::iWepwn];


        // *************| Check Values |****************************
        if (mom < 0){
            (*p_log)(LOG_ERR,AT) << "Error mom < 0! Resetting to default \n";
            mom = 1e-5;
            Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
            beta = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
        }
        if ((!std::isfinite(R)) || (!std::isfinite(Gamma)) ||
            (M2 < 0.) || (!std::isfinite(M2)) || (!std::isfinite(Eint2)) || (Eint2<0)) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " << "\n"
                                  << " R="<<R<<"\n"
                                  << " M2="<<M2<<"\n"
                                  << " Gamma=" << Gamma << "\n"
                                  << " Mom=" << mom << "\n"
                                  << " Eint2=" << Eint2
                                  << " \n";
            exit(1);
        }

        // ****| Compute the Downstream properties |*********************
        double bw_thickness = -1; // TODO interpolate it from BM-ST self-similar solution
        double bw_rho = -1;
        double bw_temp = -1;
        double bw_tau = -1;
        double bw_ye = -1;

        // ****| Check the downstream conditions |******************
        if ((!std::isfinite(bw_thickness)) || (!std::isfinite(bw_rho)) || (bw_thickness < 0.)
            || (!std::isfinite(bw_temp)) || (bw_tau < 0) || ((bw_ye < 0))) {
            (*p_log)(LOG_ERR,AT) << " bad value downstream " << "\n"
                                            << " bw_thickness="<<bw_thickness<<"\n"
                                            << " bw_rho="<<bw_rho<<"\n"
                                            << " bw_temp=" << bw_temp << "\n"
                                            << " bw_tau=" << bw_tau << "\n"
                                            << " \n";
            exit(1);
        }

        // ****| Update the PWN state |***************************
        double ldip = p_pars->curr_ldip;
        double lacc = p_pars->curr_lacc;
        if (r_w > R){ r_w = R; } // PWN should be bounded by the blastwave
        double v_w = std::sqrt( (7./6.) * e_nb / (4. * CGS::pi * r_w*r_w*r_w * bw_rho) );// (see Eq. 28 in Kashiyama+16)
        if (v_w > beta*CGS::c){ v_w = beta*CGS::c; } // velocity should also be bounded by BW TODO see if not
        double mom_wind = EQS::MomFromBeta(v_w/CGS::c);
        double dEnbdt = 0; // evaluateShycnhrotronSpectrum nebula energy \int(Lem * min(1, tau_T^ej * V_ej / c))dt
        if (bw_tau * (r_w / R) > CGS::c / v_w){ // Eq.[28] in Eq. 28 in Kashiyama+16
            dEnbdt = (p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
        }
        else{
            dEnbdt = (bw_tau * (r_w / R) * (v_w/CGS::c) * p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
        }
        double tdyn = R / beta / CGS::c; // dyn. timescale
        double dEpwndt = p_pars->eps_mag_w * ldip - e_pwn / tdyn; // adiabatic loss inclided
        double dttdr_w = EQS::evalElapsedTime(r_w, mom_wind, 0., false);
        if (!std::isfinite(dttdr_w)||(dttdr_w<=0)){
            (*p_log)(LOG_ERR,AT)<<"dttdr_w="<<dttdr_w<<"\n";
            exit(1);
        }
        // using the downstream properties, compute the fraction of PWN energy thermlized in BW
        double fac_psr_dep_tmp = getFacPWNdep(bw_rho,bw_thickness,bw_temp,bw_ye); // double rho_ej, double delta_ej, double T_ej, double Ye
        double dEindt = fac_psr_dep_tmp * (p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
        if (!std::isfinite(dEindt)||(dEindt<=0)){
            (*p_log)(LOG_ERR,AT)<<"dEindt="<<dEindt<<"\n";
            exit(1);
        }

        // ****| Get the Upstream properties |******************
        double rho = -1 / p_pars->M0; // TODO interpolate from ejecta structure as it moves freely through ISM
        double drhodr = -1 / p_pars->M0;
        double GammaRho = -1;
        double GammaRel = -1; // Gamma
        double dGammaRelDGamma = -1; // 1.
        double dGammaRhodR = -1; // 0.
        double cs_cbm = -1;
        if ((rho < 0)||(drhodr < 0)||(GammaRho<0)||(GammaRel<0)||(dGammaRelDGamma<0)||(dGammaRhodR<0)){
            (*p_log)(LOG_ERR,AT) << "Wrong value upstream"
                << " rho="<<rho<< " drhodr="<<drhodr<< " GammaRho="<<GammaRho<< " GammaRel="<<GammaRel
                << " dGammaRelDGamma="<<dGammaRelDGamma<< " dGammaRhodR="<<dGammaRhodR<<"\n";
            exit(1);
        }

        // ***| Check the upstream conditions |******************
        if (GammaRho < 1.) {
            (*p_log)(LOG_ERR,AT)  << "GammaRho=" << GammaRho << "\n" << " Exiting...";
//            std::cerr << AT<< "\n";
            exit(1);
        }
        if (GammaRel < 1.) {
            (*p_log)(LOG_ERR,AT)  << "GammaRel=" << GammaRel << "\n"; GammaRel = Gamma;
//            std::cerr << AT<< "\n";
            exit(1);
        }
        if (!std::isfinite(dGammaRelDGamma)) { (*p_log)(LOG_ERR,AT) << "dGammaRelDGamma=" << dGammaRelDGamma << "\n"; exit(1); }
        if (!std::isfinite(dGammaRhodR)) { (*p_log)(LOG_ERR,AT)  << "dlnGammaCBMdR=" << dGammaRhodR << "\n"; exit(1); }

        // ***| Update BW state |********************************
        double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...
        double Gamma_ = GammaRel;
        double beta_ = EQS::Beta(Gamma_);
        double gammaAdi  = p_eos->getGammaAdi(Gamma_, beta_);//p_eos->getGammaAdi(Gamma, beta);
        double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
        double dRshdt = EQS::Beta(GammaSh) * CGS::c;
        double dthetadr = 0.0; // NO lateral expansion for this BW

        double dM2dR = 0.;
        switch (p_pars->m_method_dmdr) {
            case iusingA:
                /// Old equation from Gavin P. Lamb work
                dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
                break;
            case iusingdthdR:
                /// Equation motivated by Gavin P. Lamb paper
                dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
                break;
            case iNodmdr:
                /// if no mass needs to be swept up
                break;
        }
        if (dM2dR < 0.){
            (*p_log)(LOG_ERR,AT) << " dMdR < 0 in RHS dyn for kN ejecta\n";
            exit(1);
        }

        // PWN energy injection
        double xi_inj = 1.; // Additional free parameters (used in dynamics eq.)
        double dEinjdt = dEindt / (p_pars->M0 * CGS::c * CGS::c) / p_pars->ncells;
        double dEinjdR = dEinjdt / dRdt;
        double theta_ej = 0.; // assume that ejecta is alinged with magnetar emission?..
        double Doppler = Gamma / (1. - beta * std::cos(theta_ej));
        double dEinjdR_dop = dEinjdR * Doppler;
        double dEingdR_abs = dEinjdR;// * ( 1. - std::exp(-1.*p_pars->dtau) ) * std::exp(-1.*p_pars->tau_to0);
        double dEingdR_abs_dop = dEingdR_abs / Doppler / Doppler; // in comoving frame

        // Nuclear heating by r-process
        double dEnuc = p_pars->dEnuc;//->getPars()->eps_nuc_thermalized;
        double dEnuc_norm = dEnuc / (p_pars->M0 * CGS::c * CGS::c);
        double dEnucdR = dEnuc_norm / dRdt;

        // Energy loss to thermal radiation
        double dElum = p_pars->dElum;
        double dElum_norm = dElum / (p_pars->M0 * CGS::c * CGS::c);
        double dElumdR = dElum_norm / dRdt;
        if ((std::abs(dElumdR) > std::abs(dEnucdR*1.5)) and (dEnuc > 0.)) // TODO These are artificial limiters
            dElumdR = dEnucdR;

        // Dynamics equation
        double GammaEff = get_GammaEff(Gamma_, gammaAdi); // TODO check if GammaRel should be used for this!!!
        double dGammaEffdGamma = get_dGammaEffdGamma(Gamma_, gammaAdi);
        double num1 = (Gamma - GammaRho + GammaEff * (GammaRel - 1.)) * dM2dR;
        double num2 = - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR/M2 - drhodr/rho - dGammaRhodR / GammaRho); // - 3.*Eint2/R
        double num3 = - (1. + GammaEff / Gamma * xi_inj) * dEingdR_abs; // injection term
        double denum1 = (1. + M2);
        double denum2 = Eint2 * dGammaEffdGamma;
        double denom3 = GammaEff * (gammaAdi - 1.) * Eint2 * dGammaRelDGamma / GammaRel;
        double dGammadR = -1. * (num1 + num2 + num3) / (denum1+denum2+denom3);
        double dmomdR = dGammadR / beta;

        // Adiabatic energy loss
        double dlnV2dR  = dM2dR / M2 - drhodr / rho - (1. / GammaRel) * dGammaRelDGamma * dGammadR + dGammaRhodR / GammaRho;// + 3./R;
        double dlnVdR = (3./R) - (1/Gamma * dGammadR);
        double dEad2dR  = 0.0;
        if ( p_pars->adiabLoss )
            dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;

        // Energy generation at the shock
        double dEsh2dR = (GammaRel - 1.0) * dM2dR;
        if ((cs_cbm > EQS::Beta(GammaRel)&&(p_pars->comp_ix > p_pars->nr*1e-2))){ // TODO Subsonic flow -- no shock
            dEsh2dR *= 1e-10;
        }

        // -- Radiative losses
        double dErad2dR = p_pars->eps_rad * dEsh2dR;

        // -- Energy equation
        double dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEnucdR - dElumdR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)

        double dtcomov_dR = 1.0 / beta / Gamma / CGS::c; // comoving time

        // -- time in a observer frame
        bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        double dttdr = EQS::evalElapsedTime(R,mom,dthetadr,spread);

        // ***| final checks |*************************

        if (mom < 0.){
            (*p_log)(LOG_ERR, AT) << " mom < 0 = "<< mom << " in kN RHS dynamics\n";
//            exit(1);
        }
        if ((!std::isfinite(dRdt)) || (!std::isfinite(dGammadR)) || (dM2dR < 0.)
            || (!std::isfinite(dlnV2dR)) || (!std::isfinite(dthetadr))) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " //<< "\n"
                                  << " dRdt="<<dRdt//<<"\n"
                                  << " dM2dR="<<dM2dR//<<"\n"
                                  << " dthetadr=" << dthetadr// << "\n"
                                  << " dGammadR=" << dGammadR //<< "\n"
                                  << " dthetadr=" << dRdt //<< "\n"
                                  << " dEsh2dR=" << dEsh2dR //<< "\n"
                                  << " dlnV2dR=" << dlnV2dR //<< "\n"
                                  << " dEad2dR=" << dEad2dR //<< "\n"
                                  << " dErad2dR=" << dErad2dR //<< "\n"
                                  << " dEint2dR=" << dEint2dR //<< "\n"
                                  << " dttdr=" << dttdr// << "\n"
                                  << " dGammaRelDGamma=" << dGammaRelDGamma //<< "\n"
                                  << " dGammaRhodR=" << dGammaRhodR //<< "\n"
                                  << " drhodr=" << drhodr //<< "\n"
                                  << " Gamma="<<Gamma//<< "\n"
                                  << " beta="<<beta//<< "\n"
                                  << " theta"<<theta//<< "\n"
                                  << " M2="<<M2//<< "\n"
                                  << " Eint2="<<Eint2//<< "\n"
                                  << " GammaRel"<<GammaRel//<< "\n"
                                  << " GammaRho"<<GammaRho//<< "\n"
                                  << " \n";
            exit(1);
        }

        // ***| Save the state |***********************
        p_pars->prev_x = x;
        out_Y[i+ SOL::QS::iR]      = dRdt;//1.0 / beta / CGS::c;
        out_Y[i+ SOL::QS::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
        out_Y[i+ SOL::QS::itt]     = dRdt * dttdr;
        out_Y[i+ SOL::QS::itcomov] = dRdt * dtcomov_dR;
        out_Y[i+ SOL::QS::imom]    = dRdt * dmomdR; //dRdt * dGammadR / Gamma;
        out_Y[i+ SOL::QS::iEint2]  = dRdt * dEint2dR;
        out_Y[i+ SOL::QS::itheta]  = dRdt * dthetadr;
        out_Y[i+ SOL::QS::iErad2]  = dRdt * dErad2dR;
        out_Y[i+ SOL::QS::iEsh2]   = dRdt * dEsh2dR;
        out_Y[i+ SOL::QS::iEad2]   = dRdt * dEad2dR;
        out_Y[i+ SOL::QS::iM2]     = dRdt * dM2dR;
        // --- PWN ---
        double drdtw = v_w;
        out_Y[i+ SOL::QS::iRw]     = drdtw;
        out_Y[i+ SOL::QS::iWtt]    = dttdr_w*drdtw;
        out_Y[i+ SOL::QS::iWenb]   = dEnbdt;
        out_Y[i+ SOL::QS::iWepwn]  = dEpwndt
    }
    /// --- Evaluate the density (and its velocity) at point ej_R left by blast wave at a point j_R
    void evalDensAndItsVelocityBehindBlastWave(double & rho, double & GammaCMB, double & p_cbm,
                                               double rho_def, double rho_prev,
                                               double ej_R, double j_R, double j_rho, double j_rho2,
                                               double j_Gamma, double j_Gamma0, double j_P2){

        GammaCMB = 1.;
        /// exponential decay from the point of entry
        double rho0 = getVal(BW::Q::irho, 0);
        /// if inside from the it=0 -- keep floor dens untill BM and ST profiles
        if ((rho0==(p_dens->m_rho_floor_val*rho_def))
            &&(rho_prev==(p_dens->m_rho_floor_val*rho_def))){
            rho = rho_prev;
        }
        else if (p_pars->use_exp_rho_decay_as_floor) {
            double coeff = p_pars->steepnes_of_exp_decay;
            double scale_rho = exp(1. - coeff * ej_R / getPars()->first_entry_r);
            double scale_dlnrhodr = -1. * coeff * scale_rho / getPars()->first_entry_r;
            if (scale_rho > 1.) {
                (*p_log)(LOG_ERR,AT)  << " scale_rho > 1 for ej_R=" << ej_R
                                      << "  entry_r=" << getPars()->first_entry_r << "\n"
                                      << " Exiting...";
//                std::cerr << AT  << "\n";
                exit(1);
            }
            rho = rho_def * scale_rho;
//            if (rho_prev < rho){
//                std::cerr <<AT<<"\n";
//                exit(1);
//            }
//            dlnrhodr = rho_def * scale_dlnrhodr;// analytc. derivative
        }
        else if ((p_pars->use_flat_dens_floor)){
            rho = rho_def * p_dens->m_rho_floor_val;
//            dlnrhodr = 0.;
        }

        /// check
        if (!std::isfinite(rho) || rho < 0. || rho > 1. || rho > rho_def ){
            (*p_log)(LOG_ERR,AT) << " wrong floor rho="<<rho<<"\n" << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        /// Blandford-Mckee tail from the jet blast wave (j_Gamma > GammaST and j_Gamma < GammaBM)
        double rho_bm = -1;
        double gam_cbm_bm = -1;
        if((ej_R > j_R * 0.01)&&(p_pars->use_bm_dens_profile)
           &&(j_Gamma > p_pars->Gamma_when_st_starts)
           &&(j_Gamma < (j_Gamma0 * p_pars->fraction_of_Gamma0_when_bm_for_bm))) {
            rho_bm = getBM()->rho_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2) * CGS::mp;
            gam_cbm_bm = getBM()->gamma_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2);
//            p_cbm_s = getSedov()->pressure_profile_int(ej_R, j_R, EQS::Beta(j_P2));
            (*p_log)(LOG_ERR,AT)  << " pressure is not implemented\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
            /// apply only of above the floor
            if (rho_bm > rho) { rho = rho_bm; GammaCMB = gam_cbm_bm; }
        }

        /// Sedov-Taylor tail from the jet blast wave (j_Gamma < GammaST)
        double rho_s = -1;
        double gam_cbm_s = -1;
        double p_cbm_s = -1.;
        if ((ej_R > j_R * 0.01) && // TODO this is needed to activate ST when kB BW is close to jet BW
            (p_pars->use_st_dens_profile)
            &&(j_Gamma <= p_pars->Gamma_when_st_starts)) {
            rho_s = getSedov()->rho_profile_int(ej_R, j_R, j_rho2);
            gam_cbm_s = getSedov()->Gamma_profile_int(ej_R, j_R, EQS::Beta(j_Gamma));
            p_cbm_s = getSedov()->pressure_profile_int(ej_R, j_R, j_P2);
            // apply only of above the floor
//            if ((rho_bm > rho_s)&&(rho_bm > rho)){
//                rho = rho_bm; GammaCMB = gam_cbm_bm;
//            }
//            else
            if (rho_s > rho) { // # THis is needed for some reason :)
                p_pars->is_using_st_prof = true;
                rho = rho_s; GammaCMB = gam_cbm_s; p_cbm = p_cbm_s;
            }


//            if (rho_prev < 1e-10 * rho){
//                std::cerr << AT << " rho_prev="<<rho_prev<<" rho="<<rho<<" rho0="<<rho0<<" density gradient >10 orders of magnitude\n";
////                exit(1);
//            }
        }



        if ((GammaCMB < 1.) || (!std::isfinite(GammaCMB))){
            (*p_log)(LOG_ERR,AT)  << " bad GammaCMB=" << GammaCMB << "\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        if ((((std::abs(rho_prev/rho) < 1e-4)||((std::abs(rho_prev/rho) > 1e4))))){
            //std::cerr << AT << " rho_prev="<<rho_prev<<" rho="<<rho<<" rho0="<<rho0<<" density gradient >4 orders of magnitude\n";
//            exit(1);
        }
    }

    // ---------------------------------------------------------
    size_t ntb() const { return m_tb_arr.size(); }
    Vector & getTbGrid() {return m_tb_arr;}
    Vector getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
        Vector tmp{};
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            tmp.push_back(m_tb_arr[it]);
        }
//        Vector tmp2 (tmp.data(), tmp.size());
        return std::move(tmp);
    }
    inline Vector & operator[](unsigned ll){ return this->m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return this->m_data[ivn][ir]; }
    inline double ctheta(double theta){
        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        if (theta > p_pars->theta_max ){
//            std::cerr << AT << " theta="<<theta<<" > theta_max=" << p_pars->theta_max << "\n";
//        }
//        if (std::fabs( theta - p_pars->theta_b0) > 1e-2){
//            std::cerr << AT << " theta="<<theta<<" < theta_b0=" << p_pars->theta_b0 <<"\n";
//            exit(1);
//        }
//        double ctheta = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w); // TODO WROOOONG

        double ctheta = 0.;
        if (p_pars->ilayer > 0) {
            //
            double fac0 = (double)p_pars->ilayer/(double)p_pars->nlayers;
            double fac1 = (double)(p_pars->ilayer+1)/(double)p_pars->nlayers;
//            std::cout << std::asin(CGS::pi*3/4.) << "\n";
            if (!std::isfinite(std::sin(theta))){
                (*p_log)(LOG_ERR,AT) << " sin(theta= "<<theta<<") is not finite... Exiting..." << "\n";
                exit(1);
            }

            double x2 = fac1*std::sin(theta / 2.);
            double xx2 = 2.*std::asin(x2);
            double x1 = fac0*std::sin(theta / 2.);
            double xx1 = 2.*std::asin(x1);

            ctheta = 0.5 * (xx1 + xx2);
            if (!std::isfinite(ctheta)){
                (*p_log)(LOG_ERR,AT) << "ctheta is not finite. ctheta="<<ctheta<<" Exiting..." << "\n";
                exit(1);
            }
        }
        return ctheta;
    }

    inline VecVector & getData(){ return m_data; }
    inline Vector & getData(BW::Q var){ return m_data[ var ]; }
    inline double & getVal(BW::Q var, int ix){
        auto ixx = (size_t)ix;
        if (ix == -1) { ixx = m_data[0].size()-1; }
        return m_data[var][ix];
    }
    inline double & getLastVal(BW::Q var){ return m_data[var][p_pars->comp_ix]; }
    void addOtherVars(size_t it){

        m_data[BW::Q::ictheta][it] = ctheta(m_data[BW::Q::itheta][it]);//p_pars->ctheta0 + 0.5 * (2. * m_data[Q::itheta][it] - 2. * p_pars->theta_w);
//        p_dens->evaluateRhoDrhoDr(m_data[Q::iR][it], m_data[Q::ictheta][it]);
        double rho_prev = m_data[BW::Q::irho][it-1];
        double rho = m_data[BW::Q::irho][it];
        if ((rho < 0)||(!std::isfinite(rho))){
            (*p_log)(LOG_ERR,AT)<<" negative density!\n";
            exit(1);
        }

        /// related to the jet BW density profile
        m_data[BW::Q::irho][it] = p_dens->m_rho_;
        m_data[BW::Q::idrhodr][it] = p_dens->m_drhodr_;
        m_data[BW::Q::iGammaCBM][it] = p_dens->m_GammaRho;
        m_data[BW::Q::iGammaREL][it] = p_dens->m_GammaRel;
        m_data[BW::Q::idGammaCBMdr][it] = p_dens->m_dGammaRhodR;
        m_data[BW::Q::idGammaRELdGamma][it] = p_dens->m_dGammaReldGamma;
        m_data[BW::Q::idPCBMdrho][it] = p_dens->m_dPCBMdrho;
        m_data[BW::Q::iPcbm][it] = p_dens->m_P_cbm;
        m_data[BW::Q::iMCBM][it] = p_dens->m_M_cbm;
        m_data[BW::Q::iCSCBM][it] = p_dens->m_CS_CBM;
        m_data[BW::Q::ijl][it]       = (double)p_pars->ijl;
        m_data[BW::Q::ir_dist][it]   = p_pars->r_dist;
        if (m_data[BW::Q::iGammaREL][it] < 1.){
            m_data[BW::Q::iGammaREL][it] = m_data[BW::Q::iGamma][it];
//            std::cerr << AT << "\n GammaRel="<<m_data[Q::iGammaREL][it]<<"; Exiting...\n";
//            exit(1);
        }

        // other parameters
        m_data[BW::Q::ibeta][it]     = EQS::Beta(m_data[BW::Q::iGamma][it]);

        if (p_pars->end_evolution)
            return;

        if ((it>1)&&(rho_prev < 1e-10 * m_data[BW::Q::irho][it])){
            (*p_log)(LOG_ERR,AT) << " it="<<it<<" density gradient >10 orders of magnitude\n";
//            exit(1);
        }

        m_data[BW::Q::iadi][it]      = p_eos->getGammaAdi(m_data[BW::Q::iGamma][it], // TODO ! is it adi or adi21 (using GammaRel)??
                                                      m_data[BW::Q::ibeta][it]);
        m_data[BW::Q::irho2][it]     = EQS::rho2t(m_data[BW::Q::iGamma][it], // TODO should there be a gammaRel?? with adi43??..
                                              m_data[BW::Q::iadi][it],
                                              m_data[BW::Q::irho][it]);
        /// shock front velocity
        switch (p_pars->m_method_gamma_sh) {

            case iuseJK:
                m_data[BW::Q::iGammaFsh][it] = EQS::GammaSh(m_data[BW::Q::iGamma][it],m_data[BW::Q::iadi][it]);
                break;
            case isameAsGamma:
                m_data[BW::Q::iGammaFsh][it] = m_data[BW::Q::iGamma][it];
                break;
            case iuseGammaRel:
                m_data[BW::Q::iGammaFsh][it] = m_data[BW::Q::iGammaREL][it];
                break;
            case iuseJKwithGammaRel:
                m_data[BW::Q::iGammaFsh][it] = EQS::GammaSh(m_data[BW::Q::iGammaREL][it],m_data[BW::Q::iadi][it]);
                break;
        }
        /// shock front radius
        switch (p_pars->m_method_r_sh) {

            case isameAsR:
                m_data[BW::Q::iRsh][it] = m_data[BW::Q::iR][it]; // overrude
                break;
            case iuseGammaSh:
                break;
        }
        /// shock thickness
        switch (p_pars->m_method_Delta) {

            case iuseJoh06:
                m_data[BW::Q::ithickness][it]= EQS::shock_delta_joh06(m_data[BW::Q::iR][it], m_data[BW::Q::iM2][it],
                                                                      m_data[BW::Q::itheta][it], m_data[BW::Q::iGamma][it],
                                                                      m_data[BW::Q::irho2][it], p_pars->ncells);
                break;
            case iuseVE12:
                m_data[BW::Q::ithickness][it]=EQS::shock_delta(m_data[BW::Q::iR][it],m_data[BW::Q::iGamma][it]);
                break;
            case iNoDelta:
                m_data[BW::Q::ithickness][it] = 1.;
                break;
        }
        /// shock downstream energy density
        switch(p_pars->m_method_up){
            case iuseEint2:
                m_data[BW::Q::iU_p][it] = EQS::get_U_p(m_data[BW::Q::irho2][it],
                                                       m_data[BW::Q::iM2][it],
                                                       m_data[BW::Q::iEint2][it]);
                break;
            case iuseGamma:
                m_data[BW::Q::iU_p][it]= EQS::get_U_p(m_data[BW::Q::irho2][it], m_data[BW::Q::iGammaFsh][it]);
                break;
        }

//        m_data[Q::imom][it] = m_data[Q::imom][it];//m_data[Q::iGamma][it] * m_data[Q::ibeta][it];


        /// evaluateShycnhrotronSpectrum time in the observer frame at 'i' step // TODO put this in the RHS -- DOne
//        bool use_spread = p_spread->m_method != LatSpread::iNULL;
//        if (m_data[Q::itt][0]==0.)
//            m_data[Q::itt][0] = EQS::init_elapsed_time(
//                    m_data[Q::iR][0], m_data[Q::iGamma][0], use_spread);
//        m_data[Q::itt][it] = m_data[Q::itt][0] + EQS::integrate_elapsed_time(
//                it, m_data[Q::iR], m_data[Q::iGamma], m_data[Q::itheta], use_spread);

        if ( ( m_data[BW::Q::iadi][it] < 1.) || (m_data[BW::Q::iR][it] < 1.) || (m_data[BW::Q::irho2][it] < 0.) ||
             ( m_data[BW::Q::iU_p][it] < 0) || (m_data[BW::Q::iGammaCBM][it] < 1.) ||  (m_data[BW::Q::iGammaFsh][it] < 1.) ) {
            std::cerr << AT << " \n Wrong value at i=" << it << " tb=" << m_data[BW::Q::itburst][it] << "\n"
                      << " ----------------------------------------------------------- \n"
                      << " end_evolution=" << p_pars->end_evolution << "\n"
                      << " which_jet_layer_to_use=" << p_pars->which_jet_layer_to_use << "\n"
                      << " first_entry_r=" << p_pars->first_entry_r << "\n"
                      << " use_flat_dens_floor=" << p_pars->use_flat_dens_floor<< "\n"
                      << " use_exp_rho_decay_as_floor=" << p_pars->use_exp_rho_decay_as_floor<< "\n"
                      << " use_st_dens_profile=" << p_pars->use_st_dens_profile<< "\n"
                      << " use_bm_dens_profile=" << p_pars->use_bm_dens_profile<< "\n"
                      << " is_within0=  " << p_pars->is_within0<< "\n"
                      << " is_within=   " << p_pars->is_within<< "\n"
                      << " ----------------------------------------------------------- \n"
                      << " ishell="<<p_pars->ishell << " ilayer=" << p_pars->ilayer<<" ii_eq="<<p_pars->ii_eq
                      << " Mom0="<<p_pars->mom0 <<" E0="<<p_pars->E0<<" theta0="<<p_pars->theta_b0
                      << " theta_max="<<p_pars->theta_max <<" ncells="<<p_pars->ncells<< "\n"
                      << " ----------------------------------------------------------- \n"
                      << " iR=          " << m_data[BW::Q::iR][it] << "\n"
                      << " itt=         " << m_data[BW::Q::itt][it] << "\n"
                      << " imom=        " << m_data[BW::Q::imom][it] << "\n"
                      << " iGamma=      " << m_data[BW::Q::iGamma][it] << "\n"
                      << " iGammaFsh=   " << m_data[BW::Q::iGammaFsh][it] << "\n"
                      << " iEint2=      " << m_data[BW::Q::iEint2][it] << "\n"
                      << " iEad2=       " << m_data[BW::Q::iEad2][it] << "\n"
                      << " iEsh2=       " << m_data[BW::Q::iEsh2][it] << "\n"
                      << " ictheta=     " << m_data[BW::Q::ictheta][it] << "\n"
                      << " irho=        " << m_data[BW::Q::irho][it] << "\n"
                      << " iEinj=       " << m_data[BW::Q::iEint2][it] << "\n"
                      //                      << " idlnrho1_dr="<< m_data[Q::idlnrho1_dr][it] << "\n"
                      << " idrhodr=     " << m_data[BW::Q::idrhodr][it] << "\n"
                      << " iGammaCBM=   " << m_data[BW::Q::iGammaCBM][it] << "\n"
                      << " iGammaREL=   " << m_data[BW::Q::iGammaREL][it] << "\n"
                      << " idGammaCBMdr=" << m_data[BW::Q::idGammaCBMdr][it] << "\n"
                      << " idGammaRELdGamma="<< m_data[BW::Q::idGammaRELdGamma][it] << "\n"
                      << " ibeta=       "      << m_data[BW::Q::ibeta][it] << "\n"
                      << " iadi=        "       << m_data[BW::Q::iadi][it] << "\n"
                      << " irho2=       "      << m_data[BW::Q::irho2][it] << "\n"
                      << " ithickness=  " << m_data[BW::Q::ithickness][it] << "\n"
                      << " iU_p=        "       << m_data[BW::Q::iU_p][it] << "\n"
                      << " imom=        "       << m_data[BW::Q::imom][it] << "\n";
            exit(1);
        }


        /// -------- energty injections
        m_data[BW::Q::ipsrFrac][it] = p_pars->facPSRdep;
        m_data[BW::Q::iLmag][it] = p_pars->dEinjdt;
        m_data[BW::Q::iLnuc][it] = p_pars->dEnuc;
        ///
        double r_w = m_data[BW::Q::i_Rw][it];
        double u_b_pwn = 3.0*m_data[BW::Q::i_Wepwn][it]/4.0/M_PI/r_w/r_w/r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
        double b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
        m_data[BW::Q::i_Wb][it] = b_pwn;
        m_data[BW::Q::i_Wdr][it] = it > 0 ? m_data[BW::Q::i_Rw][it] - m_data[BW::Q::i_Rw][it-1] : 0.;
    }

    /// Density profiles application based on the jet BW position and velocity
    void evalDensAndItsVelocityBehindBlastWave_Case1(
            double j_R, double j_Gamma, double j_ctheta, double j_rho, double j_rho2, double j_Gamma0, double j_P2,
            double ej_R, double ej_Gamma, double ej_theta ){

        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);

        if (ej_Gamma <= 1.) { ej_Gamma = 1. + 5e-5; } // TODO REMOVE

        size_t prev_ix = p_pars->comp_ix;

        // --- evaluate what density profile the ejecta blast wave would experience
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta); // set default values for density
        double rho_prev = getVal(BW::Q::irho, prev_ix-1);
        double rho_def = p_dens->m_rho_def;
        double drhodr_def = p_dens->m_drhodr_def;
        if (rho_def<(p_dens->m_rho_floor_val*rho_def)){
            (*p_log)(LOG_ERR,AT)  << " rho_def is not initialized\n Exiting...\n";
//            std::cerr<<AT<<"\n";
            exit(1);
        }
//        if ((rho_prev/p_dens->m_rho_) >= p_dens->m_rho_floor_val){
//            std::cerr << AT << " " << p_pars->is_within0 << " and "<< p_pars->is_within <<"\n";
//            std::cerr<<AT<<" rho_def is not initialized\n";
//            exit(1);
//        }
        double rho=rho_def, GammaCMB=1, P=0, GammaREL=ej_Gamma, drhodr=drhodr_def, dGammaRELdGamma=1, dGammaRhodR=0, dPdrho=0., cscbm=1., mcbm=0.;
        if (rho_prev<(p_dens->m_rho_floor_val*rho_def)){
            (*p_log)(LOG_ERR,AT) <<" rho_prev is not initialized\n Exiting...\n";
//            std::cerr<<AT<<"\n";
            exit(1);
        }
        if (rho_prev == 0){
            (*p_log)(LOG_ERR,AT)  << " rho[i] = 0" << "\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        /// eval dens seen by ejecta BW
        evalDensAndItsVelocityBehindBlastWave( rho, GammaCMB, P, rho_def, rho_prev, ej_R,  j_R,
                                               j_rho,  j_rho2,  j_Gamma, j_Gamma0, j_P2);

        /// Compute the density derivative (with radius)
        double r_prev = getVal(BW::Q::iR, prev_ix - 1);
//        double rho_prev = getVal(Q::irho, evaled_ix - 1);
        if (ej_R == r_prev){ (*p_log)(LOG_ERR,AT) << AT << " R = R_i-1 \n"; exit(1); }
        m_data[BW::Q::irho][prev_ix+1] = rho;
        m_data[BW::Q::iR][prev_ix+1] = ej_R;
        m_data[BW::Q::iPcbm][prev_ix + 1] = P;
        double dr = m_data[BW::Q::iR][prev_ix] - m_data[BW::Q::iR][prev_ix-1];
//            double dr1 = m_data[Q::iR][0] - m_data[Q::iR][1];
        drhodr = dydx(m_data[BW::Q::iR], m_data[BW::Q::irho], m_data[BW::Q::idrhodr],
                      dr, prev_ix+1, ej_R, true);
//        double drhodr2 = (rho - rho_prev) / (ej_R - r_prev);
        if (rho < p_dens->m_rho_floor_val*rho_def){ rho = p_dens->m_rho_floor_val*rho_def; drhodr = 0.; }
        if (!std::isfinite(rho) || !std::isfinite(drhodr)) {
            (*p_log)(LOG_ERR,AT)  <<" bad rho="<<rho<<" or m_drhodr="<<drhodr<<" \n Exiting...\n";
//            std::cerr << AT <<" \n";
            exit(1);
        }

        /// Compute the CBM density derivative (with radius)
        double GammaCBM_prev = getVal(BW::Q::iGammaCBM, prev_ix - 1);
        if (GammaCMB == GammaCBM_prev ){
            dGammaRhodR = 0;
        }
        else{
            m_data[BW::Q::iGammaCBM][prev_ix+1] = GammaCMB;
            dGammaRhodR = dydx(m_data[BW::Q::iR], m_data[BW::Q::iGammaCBM], m_data[BW::Q::idGammaCBMdr],
                               dr, prev_ix+1, ej_R, false);

//            double dGammaRhodR1 = (GammaCMB - GammaCBM_prev) / (ej_R - r_prev);
//            if (dGammaRhodR != dGammaRhodR1){
//                exit(1);
//            }
        }
        if (!std::isfinite(dGammaRhodR)){
            (*p_log)(LOG_ERR,AT)  << " Nan dGammaRhodR . setting to 0\n Exiting...";
//            std::cerr << AT << "\n";
            dGammaRhodR = 0.;
        }

        /// evaluateShycnhrotronSpectrum the relative velocity of the ISM (with respect to ejecta) and its derivative
//            double betaREL = EQS::BetaRel(EQS::Beta(ej_Gamma),EQS::Beta(GammaCMB));
//            GammaREL = EQS::Gamma(betaREL);
        GammaREL = EQS::GammaRel(ej_Gamma,GammaCMB);
        if (GammaREL < 1.) GammaREL = 1.;
        if ((GammaREL < 1.) || (!std::isfinite(GammaREL)) ){
            (*p_log)(LOG_ERR,AT)  << " bad GammaREL=" << GammaREL << " (too big for ejecta or nan)" << "\n"
                                  << " j_R=" << j_R << "\n"
                                  << " ej_R=" << ej_R << "\n"
                                  << " ej_theta=" << ej_theta << "\n"
                                  << " j_Gamma=" << j_Gamma << "\n"
                                  << " ej_Gamma=" << ej_Gamma << "\n"
                                  << " Gamma0=" << p_pars->Gamma0 << "\n"
                                  << " E0=" << p_pars->E0 << "\n"
                                  << " M0=" << p_pars->M0 << "\n"
                                  << " rho=" << rho << "\n"
                                  << " rho_prev=" << rho_prev << "\n"
                                  << " GammaCMB=" << GammaCMB << "\n"
                                  << " GammaCBM_prev=" << GammaCBM_prev << "\n"
                                  << " Exiting...";
//            std::cerr << "| Gamma evol: \n";
//            std::cerr << m_data[Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        double adi = p_eos->getGammaAdi(GammaCMB,EQS::Beta(GammaCMB));//5/3.;
        cscbm = sqrt(adi * P / rho) / CGS::c; // sound speed
        mcbm = sqrt(rho * EQS::Beta(GammaREL) * EQS::Beta(GammaREL) * CGS::c * CGS::c / adi / P); // mach number

        double ej_Gamma_prev = getVal(BW::Q::iGamma, prev_ix - 1);
        double GammaCMB_prev = getVal(BW::Q::iGammaCBM, prev_ix - 1);
        double GammaREL_prev = EQS::GammaRel(ej_Gamma_prev, GammaCMB_prev);
        if ((ej_Gamma == ej_Gamma_prev) || (GammaREL_prev < 1.) || (GammaCMB_prev == 1.)){
            dGammaRELdGamma = 1.;
        }
        else{
            m_data[BW::Q::iGammaREL][prev_ix+1] = GammaREL;
            double dGammaRel = m_data[BW::Q::iGammaREL][prev_ix] - m_data[BW::Q::iGammaREL][prev_ix-1];
            dGammaRELdGamma = dydx(m_data[BW::Q::iGammaREL], m_data[BW::Q::iGamma], m_data[BW::Q::idGammaRELdGamma],
                                   dr, prev_ix+1, GammaREL, false);
            dPdrho = dydx(m_data[BW::Q::iPcbm], m_data[BW::Q::irho], m_data[BW::Q::idPCBMdrho],
                          dr, prev_ix+1, P, false);
            double tmp = sqrt(P / rho) / CGS::c;
            int x = 1;
//            dGammaRELdGamma1 = (GammaREL - GammaREL_prev) / (ej_Gamma - ej_Gamma_prev);
        }
        if ((!std::isfinite(dGammaRELdGamma)) || (dGammaRELdGamma > 1.e2)){
//                std::cerr << AT << " Nan dGammaRELdGamma . setting to 1\n";
            dGammaRELdGamma = 1.;
        }

        // Finally, apply the computed density profile to be used in RHS
        if( (GammaCMB >= ej_Gamma) ){
//            std::cerr << AT << " GammaCMB="<<GammaCMB<<" > ej_Gamma="<<ej_Gamma<<"\n";
            GammaCMB = 1.; GammaREL = ej_Gamma; dGammaRhodR = 0; dGammaRELdGamma = 1.; //rho = 1e-70; m_drhodr = 0.;
        }

        p_dens->m_rho_ = rho;
        p_dens->m_drhodr_ = drhodr;
        p_dens->m_GammaRho = GammaCMB;
        p_dens->m_GammaRel = GammaREL;
        p_dens->m_dGammaRhodR = dGammaRhodR;
        p_dens->m_dGammaReldGamma = dGammaRELdGamma;
        p_dens->m_dPCBMdrho = dPdrho;
        p_dens->m_P_cbm = P;
        p_dens->m_M_cbm = mcbm; // 4.47 Zhang:2018
        p_dens->m_CS_CBM = cscbm; // 4.43 Zhang:2018

    }
    void set_no_ism(double ej_R, double ej_ctheta, double ej_Gamma_rel){
//        if (ej_Gamma_rel>10){ std::cerr << AT << "\n"; exit(1); }
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta);
        p_dens->m_rho_ = p_dens->m_rho_floor_val * p_dens->m_rho_def;
        p_dens->m_drhodr_ = 0.;

        p_dens->m_GammaRho = 1;
        p_dens->m_GammaRel = ej_Gamma_rel;
        p_dens->m_dGammaRhodR = 0;
        p_dens->m_dGammaReldGamma = 1;
        p_dens->m_dPCBMdrho = 0.;
        p_dens->m_P_cbm = 0.;
        p_dens->m_M_cbm = 0.;
        p_dens->m_CS_CBM = 0;
    }
    void set_standard_ism(double ej_R, double ej_ctheta, double ej_Gamma_rel){
//        if (ej_Gamma_rel>10){ std::cerr << AT << "\n"; exit(1); }
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta);
        p_dens->m_rho_ = p_dens->m_rho_def;
        p_dens->m_drhodr_ = p_dens->m_drhodr_def;
        p_dens->m_GammaRho = 1;
        p_dens->m_GammaRel = ej_Gamma_rel;
        p_dens->m_dGammaRhodR = 0;
        p_dens->m_dGammaReldGamma = 1;
        p_dens->m_dPCBMdrho = 0.;
        p_dens->m_P_cbm = 0.;
        p_dens->m_M_cbm = 0.;
        p_dens->m_CS_CBM = 0;
    }
    void prepareDensProfileFromJet(double * out_Y, size_t i, double x, double const * Y,
                                   void * _others, size_t evaled_ix){
//        auto * p_pars = (struct PWNPars *) params; // removing EATS_pars for simplicity
        auto * p_others = (std::vector<std::unique_ptr<BlastWave>> *) _others;
        auto & others = * p_others;
        // --- current ejecta bw
//        double ej_Gamma  = Y[p_pars->ii_eq + DynRadBlastWave::QS::iGamma];
        double ej_mom = Y[p_pars->ii_eq + SOL::QS::imom];
        if (ej_mom < 0.){
            (*p_log)(LOG_ERR,AT)<<" ej_mom < 0\n";
            exit(1);
        }
        double ej_Gamma  = EQS::GamFromMom(ej_mom);
//        double ej_Gamma = std::exp(ej_lnGamma);
        if (ej_Gamma < 1.) ej_Gamma = 1.0000000099999;
//        if (ej_Gamma > 10.) {
//            std::cerr << AT << " ej_Gamma="<<ej_Gamma<<"\n";
//            exit(1);
//        }
        double ej_R      = Y[p_pars->ii_eq + SOL::QS::iR];
        double theta_b0  = p_pars->theta_b0;
        double ej_theta  = Y[p_pars->ii_eq + SOL::QS::itheta];
        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);
        // -- loop over jets

        int i_ej_l = p_pars->which_jet_layer_to_use;//others.size()-1;
        bool is_within = false;
//        for (size_t ij = 0; ij < others.size(); ij++){
//            double j_R      = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            double j_Gamma  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//            double j_theta  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::itheta];
//            double j_ctheta = others[ij]->ctheta(j_theta);
//            double j_theta0 = others[ij]->theta0(j_theta); //
//            double j_theta1 = others[ij]->theta1(j_theta);
//            if ((ej_ctheta > j_theta0) && (ej_ctheta <= j_theta1)){
////            if ((ej_ctheta < j_ctheta)){
//                i_ej_l = ij; is_within = true;
////                if(i_ej_l!=p_pars->j_i0){
////                    std::cerr << AT<< " \nj_theta0"<<p_pars->j_theta0<<" j_theta1="<<p_pars->j_theta1<<" prev_layer="<<p_pars->j_i0<<"\n";
////                    std::cerr << "ctheta="<<j_theta0<<" ctheta0="<< p_pars->ctheta0 <<"\n";
////                    exit(1);
////                }
//            }
//        }
//        if (!is_within){
//            double j_R      = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            double j_Gamma  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//            double j_theta  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::itheta];
//            double j_ctheta = others[0]->ctheta(j_theta);
//            double j_theta0 = others[0]->theta0(j_theta); //
//            double j_theta1 = others[0]->theta1(j_theta);
//            if (ej_ctheta < j_theta0){
//                i_ej_l = 0.;
//                is_within = true;
//            }
//        }

//        i_ej_l = 0;//others.size()-1;
        double j_mom = Y[others[i_ej_l]->getPars()->ii_eq + SOL::QS::imom];
        double j_Gamma = EQS::GamFromMom(j_mom);
//        double j_Gamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//        double ibeta = EQS::BetFromMom(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::imom]);
        double j_theta = Y[others[i_ej_l]->getPars()->ii_eq + SOL::QS::itheta];
//        double j_lnGamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::ilnGamma];
//        double j_Gamma = std::exp(j_lnGamma);
        double j_Gamma0 = others[i_ej_l]->getPars()->Gamma0;
//        double j_ctheta = others[i_ej_l]->ctheta(j_theta);
        double j_ctheta = EjectaID2::ctheta(j_theta,
                                            others[i_ej_l]->getPars()->ilayer,
                                            others[i_ej_l]->getPars()->nlayers);//others[i_ej_l]->ctheta(j_theta);
        if (ej_ctheta < j_ctheta){ is_within = true; }

        auto & other = others[i_ej_l];
        double j_R = Y[other->getPars()->ii_eq + SOL::QS::iR];

        if ((ej_R == p_pars->R0)||(evaled_ix == 0)){
//            if (evaled_ix > 1){
//                std::cerr<<AT<<" \n";exit(1);
//            }
            if (is_within){
                // inside -- NO ISM
                set_no_ism(ej_R, ej_ctheta, ej_Gamma);
//                p_pars->is_within0 = is_within;
            }
            else {
                // outside -- standard ISM
                set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
//                p_pars->is_within0 = is_within;
            }
//            if ((p_pars->ilayer==1)&&(p_dens->m_rho_!=(p_dens->m_rho_floor_val*p_dens->m_rho_def))){
//                exit(1);
//            }
        }
        else if ((ej_R < j_R)&&(ej_R > p_pars->R0)){
            if (is_within && (evaled_ix > 0)){
                // inside and behind -- complex ISM
                if ((p_pars->first_entry_r < 0) || (p_pars->first_entry_r > ej_R))
                    p_pars->first_entry_r = ej_R;
                // if entering a new jet layer
                if (i_ej_l != p_pars->ijl) {
                    p_pars->prev_ijl = p_pars->ijl; // -1 0 1 2 3 4
                    p_pars->ijl = i_ej_l; // 0 2 3 4 5
                }

                size_t other_i = others[i_ej_l]->getPars()->ii_eq;//other->getPars()->ii_eq;
                double other_Gamma = EQS::GamFromMom(Y[other_i + SOL::QS::imom]);//std::exp( Y[other_i + DynRadBlastWave::QS::ilnGamma] );
//                double other_Gamma = Y[other_i + DynRadBlastWave::QS::iGamma];//std::exp( Y[other_i + DynRadBlastWave::QS::ilnGamma] );
                // --- density experienced by the jet blast wave and the density in the downstream (rho2)
                other->getDensIsm()->evaluateRhoDrhoDrDefault(j_R, INFINITY);
                double j_rho = other->getDensIsm()->m_rho_def;
                double j_drhodr = other->getDensIsm()->m_drhodr_def;
                double j_adi = other->getEos()->getGammaAdi(other_Gamma, EQS::Beta(other_Gamma));
                double j_rho2 = EQS::rho2t(j_Gamma, j_adi, j_rho);
                double j_V2 = Y[other_i + SOL::QS::iM2] / j_rho2; // Units -> c^2 for energy
                double j_P2 = (j_adi - 1.) * Y[other_i + SOL::QS::iEint2] / j_V2 * CGS::c * CGS::c;//# * CGS::c * CGS::c; // Units -> c^2 for energy
                double cs = sqrt(j_adi * j_P2 / j_rho2) / CGS::c;
                double cs2 = sqrt(j_adi*(j_adi-1)*(j_Gamma-1)/(1+j_adi*(j_Gamma-1)));
                // ---
//                size_t i = p_pars->ii_eq;
//                double ej_Gamma = Y[i + DynRadBlastWave::QS::iGamma];
//                double ej_R = Y[i + DynRadBlastWave::QS::iR];
//                double theta_b0 = p_pars->theta_b0;
//                double ej_theta = Y[i + DynRadBlastWave::QS::itheta];
//                double ej_ctheta = ctheta(ej_theta);
                if (ej_Gamma <= 1.) { ej_Gamma = 1.0000000000999; } // --------------------------- [ kostil ]

                evalDensAndItsVelocityBehindBlastWave_Case1(j_R, j_Gamma, j_ctheta, j_rho, j_rho2, j_Gamma0, j_P2, ej_R, ej_Gamma, ej_theta);
//                evalDensAndItsVelocityBehindBlastWave_Case1( Y, other );
            }
            else{
                // outside -- standard ISM
                set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
            }
        }
        else if (ej_R >= j_R){
            // outside -- standard ISM
            set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
        }
        else{
            (*p_log)(LOG_ERR,AT)  << " should not be entered. Exiting..."<<"\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

//        if ((p_pars->is_within)&&(!is_within)){
//            exit(1);
//        }
//        double rho =  p_dens->m_rho_;
//        p_dens->evaluateRhoDrhoDr(ej_R, ej_ctheta); // set default values for density
//        double rho_def = p_dens->m_rho;
//        double drhodr_def = p_dens->m_drhodr;
        double rho_prev    = getVal(BW::Q::irho, evaled_ix-1);
        double drhodr_prev = getVal(BW::Q::idrhodr, evaled_ix-1);
        double prev_ej_R   = getVal(BW::Q::iR, evaled_ix-1);
//        if ((rho<1e-30)&&(evaled_ix>1)&&(rho_prev<rho)){
//            std::cerr<<AT<<" rho_prev is not initialized\n";
//            exit(1);
//        }
        p_pars->is_within = is_within;
        p_pars->r_dist = j_R - ej_R;

//        if (p_pars-> < pfsolvePars->dens_floor_frac*rho_def)
        /// evaluate the RHS

//        if ((p_pars->prev_ijl!=i_ej_l)&&(evaled_ix>20)){
//            std::cerr << " ishell="<<p_pars->ishell
//                      << " ilayer="<<p_pars->ilayer
//                      << " Gamma0="<<p_pars->Gamma0
//                      << " Gamma="<<ej_Gamma
//                      << " rho_prev"<<rho_prev
//                      << " drhodr_prev"<<drhodr_prev
//                      << " rho="<<p_dens->m_rho
//                      << " drhodr="<<p_dens->m_drhodr
//                      << " ej_ctheta="<<ej_ctheta
//                      << " ej_ctheta_prev="<<getVal(Q::ictheta, evaled_ix-1)
//                      << " j_theta0="<<others[ijl]->theta0(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::itheta])
//                      << " j_theta1="<<others[ijl]->theta1(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::itheta])
//                      << "\n";
//        }

//        if (((rho_prev/p_dens->m_rho_) >= p_dens->m_rho_floor_val)){
//            exit(1);
//        }
//        if ((p_pars->ilayer==1)&&(p_dens->m_rho_!=(p_dens->m_rho_floor_val*p_dens->m_rho_def))){
//            exit(1);
//        }
//        if (p_pars->is_within0 != p_pars->is_within){
//            exit(1);
//        }
    }
//
    void evaluateRhsDensModel2(double * out_Y, size_t i, double x, double const * Y, void * others, size_t evaled_ix) {

        /// do not evaluate RHS if the evolution was terminated
        if (p_pars->end_evolution) {
            return;
        }

        /// evaluateShycnhrotronSpectrum density profile in front of the kN BW
        if (p_pars->use_dens_prof_behind_jet_for_ejecta) {
            prepareDensProfileFromJet(out_Y, i, x, Y, others, evaled_ix);
        }
        else{
            double ej_Gamma  = EQS::GamFromMom( Y[i + SOL::QS::imom] );
//            double ej_Gamma  = Y[i + DynRadBlastWave::QS::iGamma];
            if (ej_Gamma < 1.) {
                (*p_log)(LOG_ERR,AT) << "Gamma < 1\n";
                ej_Gamma = 1. + 1e-5;
            }
            double ej_R      = Y[i + SOL::QS::iR];
            double theta_b0  = p_pars->theta_b0;
            double ej_theta  = Y[i + SOL::QS::itheta];
            double ej_ctheta = EjectaID2::ctheta(ej_theta,p_pars->ilayer,p_pars->nlayers);//ctheta(ej_theta);
            set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
        }

        /// evaluate actual RHS
        evaluateRhsDens(out_Y, i, x, Y);
//        evaluateRhs(out_Y, i, x, Y);

    }



    /// ----------------- BLAST WAVE RADIATION ------------
    void addComputeForwardShockMicrophysics(size_t it){

        if (it > m_data[BW::Q::itburst].size()){
            (*p_log)(LOG_ERR,AT)<< " it="<<it<<" out of range=m_data[BW::Q::itburst].size()="
                                <<m_data[BW::Q::itburst].size()<<" mtburst.size()="<<p_pars->nr<<"\n";
            exit(1);
        }
//        auto & p_eats = p_pars; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;
        if ((m_data[BW::Q::ibeta][it] < p_syna->getPars()->beta_min)){
            if (p_pars->i_end_r > it-1)
                p_pars->i_end_r = it-1;
            return;
        }
        if ((m_data[BW::Q::iCSCBM][it] >= m_data[BW::Q::ibeta][it])){
            m_data[BW::Q::iB][it] = 0.;
            if (p_pars->i_end_r > it-1)
                p_pars->i_end_r = it-1;
            return;
        }
        // no ISM case -> no electron acceleration
        if ((m_data[BW::Q::iU_p][it] <= 0) || (m_data[BW::Q::irho2][it] <= 0)){
            if (p_pars->i0_failed_elecctrons == 0)
                p_pars->i0_failed_elecctrons = it;
            p_pars->n_fialed_electrons += 1;
//            std::cerr << AT << " at it="<<it<<" U_e="<<m_data[Q::iU_e][it]<<" and rho2="<<m_data[Q::irho2][it]<<" skipping electron calc.\n";
            return;
        }
//        if(m_data[Q::iGammaREL][it]<=1.){
//            std::cerr << AT << " at it="<<it<<" iGammaREL="<<m_data[Q::iGammaREL][it]<<" skipping electron calc.\n";
//            exit(1);
//        }
        double Gamma_; // TODO should be GammaSh
        if (m_data[BW::Q::iGammaREL][it] > 0.) Gamma_ = m_data[BW::Q::iGammaREL][it];
        else Gamma_ = m_data[BW::Q::iGamma][it];

        /// for observer frame evaluation, the rad. needs TT, while for comov. spectrum it needs tcomov
        // TODO Check if time is TT or tburst here
//        Array * tmp_time;
//        if (p_syn->getPars()->method_comp_mode==SynchrotronAnalytic::METHODS_RAD::icomovspec)
//            *tmp_time = m_data[Q::itt][it];
//        else
//            *tmp_time = m_data[Q::itburst][it];
        p_syna->evaluateElectronDistribution(m_data[BW::Q::iU_p][it],
                                             m_data[BW::Q::iGamma][it],
//                               m_data[Q::iGamma][it],
                                             m_data[BW::Q::iGammaFsh][it],
                                             m_data[BW::Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               m_data[Q::itburst][it], // emission time (TT)
                                             m_data[BW::Q::irho2][it] / CGS::mp);
        // adding
        m_data[BW::Q::iB][it]      = p_syna->getPars()->B;
        m_data[BW::Q::igm][it]     = p_syna->getPars()->gamma_min;
        m_data[BW::Q::igM][it]     = p_syna->getPars()->gamma_max;
        m_data[BW::Q::igc][it]     = p_syna->getPars()->gamma_c;
        m_data[BW::Q::iTheta][it]  = p_syna->getPars()->Theta;
//        m_data[Q::ix][it]      = p_syn->getPars()->x; // IT is not evaluated/ It depends in freq
        m_data[BW::Q::iz_cool][it] = p_syna->getPars()->z_cool;
        m_data[BW::Q::inprime][it] = p_syna->getPars()->n_prime;
        m_data[BW::Q::iacc_frac][it] = p_syna->getPars()->accel_frac;
//        if ( m_data[Q::iGamma][it] < 5 ){
//            exit(1);
//        }
        if ((!std::isfinite( m_data[BW::Q::igm][it] )) || (m_data[BW::Q::iB][it] < 0.)) {
            (*p_log)(LOG_ERR,AT) << " Wrong value at i=" << it << " tb=" << m_data[BW::Q::itburst][it]
                                 << " iR=" << m_data[BW::Q::iR][it] << "\n"
                                 << " iGamma=" << m_data[BW::Q::iGamma][it] << "\n"
                                 << " ibeta=" << m_data[BW::Q::ibeta][it] << "\n"
                                 << " iM2=" << m_data[BW::Q::iM2][it] << "\n"
                                 << " iEint2=" << m_data[BW::Q::iEint2][it] << "\n"
                                 << " iU_p=" << m_data[BW::Q::iU_p][it] << "\n"
                                 << " irho=" << m_data[BW::Q::irho][it] << "\n"
                                 << " irho2=" << m_data[BW::Q::irho2][it] << "\n"
                                 << " iB=" << m_data[BW::Q::iB][it] << "\n"
                                 << " igm=" << m_data[BW::Q::igm][it] << "\n"
                                 << " igM=" << m_data[BW::Q::igM][it] << "\n"
                                 << " igc=" << m_data[BW::Q::igc][it] << "\n"
                                 << " iTheta=" << m_data[BW::Q::iTheta][it] << "\n"
                                 //                      << " ix="     << m_data[Q::ix][it] << "\n"
                                 << " iz_cool=" << m_data[BW::Q::iz_cool][it] << "\n"
                                 << " inprime=" << m_data[BW::Q::inprime][it]
                                 << "\n";
            exit(1);
        }
    }
    void computeForwardShockElectronAnalyticVars(){
        for (size_t it = 0; it < p_pars->nr; it++){
//                std::cout << nr<<"\n";
            addComputeForwardShockMicrophysics(it);
        }
        /// check if electron spectrum failed for any reason
        if ((m_data[BW::Q::iR][0]>0)&&(p_pars->n_fialed_electrons == p_pars->nr)&&(!p_pars->end_evolution)){
            (*p_log)(LOG_ERR,AT)
                << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                << " Electron calculation failed for all iterations. Exiting...\n";
            exit(1);
        }
        else if ((m_data[BW::Q::iR][0]>0)&&(p_pars->n_fialed_electrons > 0)&&(p_pars->n_fialed_electrons < p_pars->nr)){
            (*p_log)(LOG_ERR,AT)
                << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                <<" Electron calculation failed for n=" << p_pars->n_fialed_electrons
                << " iterations starting from it=" << p_pars->i0_failed_elecctrons<<"\n";
        }
    }
    void computeForwardShockSynchrotronAnalyticSpectrum(){
        if (p_pars->m_method_rad==METHODS_RAD::icomovspec) {
            (*p_log)(LOG_INFO,AT) << " computing analytic comoving spectrum\n";
            for (size_t it = 0; it < p_pars->nr; ++it) {
                //        auto & p_eats = p_pars; // removing EATS_pars for simplicity
                /// exit if the obs. radiation method of choice does not need comoving spectrum
                if (p_pars->m_freq_arr.size() < 1){
                    (*p_log)(LOG_ERR,AT) << " array for comoving spectrum is not initialized \n Exiting...\n";
                    exit(1);
                }
                if (p_pars->m_synch_em.size() < 1){
                    (*p_log)(LOG_ERR,AT)<< " array for comoving spectrum frequencies is not initialized \n Exiting...\n";
                    exit(1);
                }
                size_t nfreq = p_pars->m_freq_arr.size();
//        (*p_log)(LOG_INFO) << " computing comoving intensity spectum for "
//                              << p_syna->getFreqArr().size() << " grid";
                /// -- check if there are any data first
                double beta_;
                if (m_data[BW::Q::iGammaREL][it] > 0.) beta_ = EQS::Beta( m_data[BW::Q::iGammaREL][it] );
                else beta_ = m_data[BW::Q::ibeta][it];

                auto & p_syna = p_pars->p_syna;
                /// if BW is not evolved to this 'it' or velocity is smaller than minimum
                if ((m_data[BW::Q::iR][it] < 1) || (beta_ < p_syna->getPars()->beta_min))
                    return;

                if (m_data[BW::Q::igm][it] == 0){
                    (*p_log)(LOG_ERR,AT)<< " in evolved blast wave, found gm = 0" << "\n";
                    exit(1);
                }

                /// evaluateShycnhrotronSpectrum emissivity and absorption for each frequency
                for (size_t ifreq = 0; ifreq < p_pars->m_freq_arr.size(); ++ifreq){
                    /// evaluateShycnhrotronSpectrum all types of emissivities and absoprtions
                    p_syna->evaluateShycnhrotronSpectrum(
                            m_data[BW::Q::irho2][it] / CGS::mp,//m_data[Q::iM2][it] / CGS::mp,
                            m_data[BW::Q::iM2][it] / CGS::mp,//m_data[Q::iM2][it] / CGS::mp,
                            m_data[BW::Q::iacc_frac][it],
                            m_data[BW::Q::iB][it],
                            m_data[BW::Q::igm][it],
                            m_data[BW::Q::igM][it],
                            m_data[BW::Q::igc][it],
                            m_data[BW::Q::iTheta][it],
                            m_data[BW::Q::iz_cool][it],
//                                       m_data[Q::irho2][it] / CGS::mp,
                            p_pars->m_freq_arr[ifreq]
                    );
                    /// add evaluated data to the storage
//            double thick_tau = EQS::shock_delta(m_data[Q::iRsh][it],m_data[Q::iGammaFsh][it]);
//            p_syna->addIntensity(thick_tau, 0., 1.);
                    p_pars->m_synch_em[ifreq + nfreq * it] = p_syna->getPars()->em;
                    p_pars->m_synch_abs[ifreq + nfreq * it] = p_syna->getPars()->abs;
                }
            }
        }
    }

#if 0
    static void optDepthPW(double & tau_Compton, double & tau_BH, double & tau_bf, double & r, double & ctheta,
                           size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                           Vector & ttobs, void * params){

        auto * p_pars = (struct PWNPars *) params;
        auto & m_data = p_pars->m_data;
        if (p_pars->i_end_r==0)
            return;

        /// interpolate ejecta properties for a given time
        double delta_ej = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iEJdelta]);
        double rho_ej = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iEJrho]);
        double Gamma = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iGamma]);
        double beta = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::ibeta]);
        ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BW::Q::ictheta]);

        /// account for relativistic motion of the shell
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
        double nu_erg = nu_obs*4.1356655385381E-15*CGS::EV_TO_ERG;

        /// evaluate optical depths
        double e_gamma = nu_erg; //! TODO CHECK!;
        double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/p_pars->mu_e/CGS::M_PRO;
        tau_Compton = rho_ej*delta_ej*Kcomp;
        /// optical depth of BH pair production
//        double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/CGS::M_PRO/r_ej/r_ej;
        double KBH = (1.0+p_pars->Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/p_pars->mu_e/CGS::M_PRO;
        tau_BH = rho_ej*delta_ej*KBH;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
//        double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma, Z_eff, opacitymode)/r_ej/r_ej;
        double Kbf = (1.0-p_pars->albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, p_pars->Z_eff, p_pars->opacitymode);
        tau_bf = rho_ej*delta_ej*Kbf;

    }
#endif

//    static void fluxDensPW(double & flux_dens, double r, double & ctheta, double theta, double phi,
//                    size_t ia, size_t ib, double ta, double tb, double mu, double t_obs, double nu_obs, void * params){

    static void fluxDensPW(double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
                           double r, double & ctheta, double theta, double phi,
                           size_t ia, size_t ib, double ta, double tb, double mu,
                           double t_obs, double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params;
        auto & m_data = p_pars->m_data;
        if (p_pars->i_end_r==0)
            return;

        if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
            Interp2d int_em(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_em);
            Interp2d int_abs(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;

            double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
            double beta = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ibeta]);
            // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor
            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z) * nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_pars->m_freq_arr, p_pars->m_freq_arr.size());
            size_t ib_nu = ia_nu + 1;
            /// interpolate the emissivity and absorption coefficines
//                double em_prime = int_em.Interpolate(nuprime, r, mth);
            double em_prime = int_em.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
//                double abs_prime = int_abs.Interpolate(nuprime, r, mth);
            double abs_prime = int_abs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
            /// convert to the laboratory frame
            double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
            double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

            /// evaluateShycnhrotronSpectrum optical depth (for this shock radius and thickness are needed)
            double GammaShock = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGammaFsh]);
            double dr = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness]);
            double dr_tau = EQS::shock_delta(r,
                                             GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

            double beta_shock;
            switch (p_pars->method_shock_vel) {

                case isameAsBW:
                    beta_shock = EQS::Beta(Gamma);
                    break;
                case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = RadiationBase::optical_depth(abs_lab, dr_tau, mu, beta_shock);
            double intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                               p_pars->p_syna->getPars()->method_tau);
            flux_dens = (intensity * r * r * dr) * (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
//        flux += flux_dens;
            /// save the result in image
            ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
            //  double ctheta = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::ictheta] ) ).Interpolate(r, mth );
//        image(Image::iintens, i) =
//                flux_dens / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
//        image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::imu, i) = mu_arr[i];
        }
        else{
            double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
            double GammaSh = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGammaFsh]);
            double rho2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::irho2]);
            double m2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iM2]);
            double frac = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iacc_frac]);
            double thick = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness]);
            theta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::itheta]);
            ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
            double B = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iB]);
            double gm = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igm]);
            double gM = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igM]);
            double gc = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igc]);
            double Theta = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iTheta]);
            double z_cool = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iz_cool]);
            double tburst = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itburst]);
            double tt = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itt]);
            double cs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iCSCBM]);

            if ((!std::isfinite(gm))||(!std::isfinite(B))||(!std::isfinite(m2))){

                (*p_pars->p_log)(LOG_ERR,AT)
                    <<"[ish="<<p_pars->ishell<<", "<<"il="<<p_pars->ilayer<<"] "
                    <<" nanss {"
                    <<" ia="<<ia<<" ib="<<ib<<" ta="<<ta<<" tb="<<tb<<" r="<<r<<" mu="<<mu
                        <<" nu_obs="<<nu_obs<<" t_obs="<<t_obs
                        <<" phi="<<phi<<" theta="<<theta<<" ctheta="<<ctheta<<" flux_dens="<<flux_dens
                        << " | " << " Gamma="<<Gamma<<" rho2="<<rho2<<" m2="<<m2<<" B="<<B
                    <<"\n";
                    ;
                exit(1);
            }

            double dFnu = 0.;
            if ((m_data[BW::Q::iB][ia] == 0.) || (m_data[BW::Q::iB][ib] == 0.)){
                dFnu = 0.;
            }
            else{
#if 0
                if ((Gamma < 1. || !std::isfinite(Gamma))
                        || (gm < 0.) || (!std::isfinite(gm))
                        || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
                        || (theta <= 0.)
                        || (rho2 < 0.) || (!std::isfinite(rho2))
                        || (thick <= 0.) || (GammaSh <= 1.)) {
                        (*p_log)(LOG_ERR,AT)<< " Error in interpolation to EATS surface: \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " GammaSh = " << GammaSh << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n";
//                    exit(1);
                    }
                    if ((B != 0.) && (!std::isfinite(rho2))) {
                        (*p_log)(LOG_ERR,AT)<< " B!=0 and rho2 is NAN \n"
                                            << " Error in data \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n"
                                            << " Exiting...\n";
                        exit(1);
                    }
#endif
                double thick_tau = EQS::shock_delta(r,GammaSh); // TODO this is added becasue in Johanneson Eq. I use ncells
                flux_dens = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                           Theta, z_cool, t_obs, mu,
                                                           r, thick,  thick_tau, nu_obs, p_pars);
            }
            flux_dens *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
        }
    }

    static void fluxDensA(double & flux_dens, double & r, double & ctheta, double theta, double phi,
                          size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;//->getAnSynch();
        auto & m_data = p_pars->m_data;
        auto & tburst = m_data[BW::Q::itburst];
        auto & r_arr = m_data[BW::Q::iR];

        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            Interp2d int_em(p_pars->m_freq_arr, r_arr, p_pars->m_synch_em);
            Interp2d int_abs(p_pars->m_freq_arr, r_arr, p_pars->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r = interpSegLog(ia, ib, t_e, tburst, r_arr);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                             << " Current R grid us ["
                                             << r_arr[0] << ", "
                                             << r_arr[tburst.size() - 1] << "] "
                                             << "and tburst arr ["
                                             << tburst[0] << ", " << tburst[p_pars->nr - 1]
                                             << "] while the requried obs_time=" << t_obs
                                             << "\n";
//                std::cerr << AT << "\n";
                return;
            }
            double Gamma = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::iGamma]);
            double beta = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::ibeta]);
            // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor
            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_pars->m_freq_arr, p_pars->m_freq_arr.size());
            size_t ib_nu = ia_nu + 1;
            /// interpolate the emissivity and absorption coefficines
//                double em_prime = int_em.Interpolate(nuprime, r, mth);
            double em_prime = int_em.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
//                double abs_prime = int_abs.Interpolate(nuprime, r, mth);
            double abs_prime = int_abs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
            /// convert to the laboratory frame
            double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
            double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

            /// evaluateShycnhrotronSpectrum optical depth (for this shock radius and thickness are needed)
            double GammaShock = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::iGammaFsh]);
            double dr = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::ithickness]);
            double dr_tau = EQS::shock_delta(r, GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

            double beta_shock;
            switch (p_pars->method_shock_vel) {

                case isameAsBW:
                    beta_shock = EQS::Beta(Gamma);
                    break;
                case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
            double intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                               p_syna->getPars()->method_tau);
            flux_dens = (intensity * r * r * dr); //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
//            dFnu+=flux_dens;
            /// save the result in image
            ctheta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::ictheta]);
            double theta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::itheta]);
        }
        else{
            double R = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iR]);
            if (!std::isfinite(R)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R is NAN in integrand for radiation" << "\n";
                // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << m_data[BW::Q::iR] << "\n";
//            std::cout << " Gamma= " << m_data[BW::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
//                return 0.;
                exit(1);
            }

            double rho = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::irho]);
            double Gamma = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::iGamma]);
            double GammaSh = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                          m_data[BW::Q::iGammaFsh]);
            double beta = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::ibeta]);
            double U_p = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM2));
            double theta = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::itheta]);
            double rho2 = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::irho2]);
            double m2 = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iM2]);
            double frac = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                       m_data[BW::Q::iacc_frac]);
            double thick = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::ithickness]);
            double gm = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igm]);
            double gM = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igM]);
            double gc = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igc]);
            double B = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iB]);
            double Theta = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::iTheta]);
            double z_cool = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                         m_data[BW::Q::iz_cool]);

            if (rho < 0. || Gamma < 1. || !std::isfinite(Gamma)
                || U_p < 0. || theta <= 0. || rho2 < 0. || thick <= 0.) {
                std::cerr << " wrong value in interpolation to EATS surface  \n"
                                             << " R = " << R << "\n"
                                             << " rho = " << rho << "\n"
                                             << " Gamma = " << Gamma << "\n"
                                             << " U_p = " << U_p << "\n"
                                             << " theta = " << theta << "\n"
                                             << " rho2 = " << rho2 << "\n"
                                             << " thick = " << thick << "\n"
                                             << " t_e = " << t_e << "\n"
                                             << " Exiting...\n";
                exit(1);
            }

            flux_dens = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2,
                                                  frac, B, gm, gM, gc, Theta, z_cool,
                                                  t_e, mu, R, thick, thick, nu_obs, params);
//            flux_dens*=(p_pars->d_l*p_pars->d_l*2.);
#if 0
            /* -- Reverse shock --- */
            double dFnu_rs = 0.0;
            if (p_pars->synch_rs){
                std::cout << AT << " Warning! Reverse shock is not finished!" << "\n";
                double rho4   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::irho4));
                double U_e3   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iU_e3));
                double rho3   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::irho2));
                double thick3 = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::ithickness));
    //            double M3     = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM3));
                double gamma43= interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iGamma43));
    //            double gammaAdi_rs = p_pars->p_eos->getGammaAdi(gamma43, EQS::Beta(gamma43));
    //            double nprime3 = 4.0 * Gamma * (rho4 / CGS::mppme);
    //            double nprime3 = p_pars->eq_rho2(Gamma, rho4 / CGS::mp, gammaAdi_rs); // TODO check if for RS here is Gamma!
                double nprime3 = rho3 / CGS::mp;
                /// evaluateShycnhrotronSpectrum the 'thickness' of the shock (emitting region)
    //            double dr_rs = thick3;

    //          // TODO check if for the shock velocity gamma43 has to be used!!! I think it is gam43! See gammaAdi calc.
                dFnu_rs = Observables::shock_synchrotron_flux_density(
                        Gamma, gamma43, rho3, U_e3, t_e, mu, R, thick3, thick3, params );
            }
            dFnu += dFnu_rs;
#endif
            if (flux_dens == 0 || !std::isfinite(flux_dens)) {
                // REMOVING LOGGER
                std::cerr << " flux density is zero ( dFnu = 0 )" << "\n";
            }

//            p_pars->o_gam = Gamma;
//            p_pars->o_r = R;
//            p_pars->o_mu = mu;
//            p_pars->o_flux = dFnu;
//            p_pars->o_theta_j = theta;
        }
    }

    static double shock_synchrotron_flux_density(
            double Gamma, double GammaShock, double m2, double rho2, double acc_frac, double B,
            double gm, double gM, double gc, double Theta, double z_cool,
            double t_e, double mu, double R, double dr, double dr_tau,
            double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;//->getAnSynch();

        /// relativistic effects on the emitting region
        double beta = EQS::Beta(Gamma);
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        double beta_shock;
        switch (p_pars->method_shock_vel) {

            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                break;
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
        dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
        dr_tau /= ashock;

        /// properties of the emitting electrons
        double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
        double Ne = m2 / CGS::mp; // numer of protons/electrons
        double nprime = rho2 / CGS::mp; // number density of protons/electrons
        double em_prime,em_lab,abs_prime,abs_lab,intensity,flux_dens,dtau;


#if 1
        switch (p_pars->m_method_ne) {
            // default (presumably more correct version)
            case iusenprime:
                /// this is informed by van Earten et al. and afterglowpy
                p_syna->evaluateShycnhrotronSpectrum(nprime, Ne, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime);
                em_prime = p_syna->getPars()->em;
                em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                abs_prime = p_syna->getPars()->abs;
                abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_syna->getPars()->method_tau);
//                flux_dens = (intensity * R * R * dr) * (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);

                flux_dens = (intensity * R * R * dr) ;
                break;
            case iuseNe: // TODO THIS SHOULD NOT BE USED (tau is incorrectly estimated)
                /// This is informed by the G. Lamb and Fernandez et al.
                p_syna->evaluateShycnhrotronSpectrum(Ne, Ne, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime);
                em_prime = p_syna->getPars()->em;
                em_lab = em_prime / (delta_D * delta_D);
                em_lab /= delta_D; // TODO this should be from 'dr'...
                abs_lab = abs_prime * delta_D; // TODO with Ne this might not work as we do not use 'dr' of the shock...
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau,mu,beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_syna->getPars()->method_tau);
//                flux_dens = intensity * (1.0 + p_eats->z) / (p_eats->d_l * p_eats->d_l) / 10;
                flux_dens = intensity / 5.; // TODO why no '2'?? // why this need /10 to fit the upper result?
                break;
        }
#endif
        return flux_dens;
    }

    bool evalEATSindexes(size_t &ia, size_t &ib, double t_obs, double z, //size_t m_i_end_r,
                         double ctheta, double cphi, double theta_obs,
                         double (*obs_angle_func)( const double &, const double &, const double & )){

//        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & m_mu =  p_pars->m_mu;
        auto & ttobs =  p_pars->ttobs;
        size_t m_i_end_r= p_pars->i_end_r;
        if (m_mu.size() < 1) {
            m_mu.resize(p_pars->i_end_r);
            ttobs.resize(p_pars->i_end_r);
        }
        /// evaluate mu array
        for (size_t i = 0; i < m_i_end_r; i++)
            m_mu[i] = ( m_data[BW::Q::itburst][i] - t_obs / (1.0 + z) ) / m_data[BW::Q::iR][i] * CGS::c;
        double mu = obs_angle_func(ctheta, cphi, theta_obs);
        for (size_t i_ = 0; i_ < m_i_end_r; i_++) {
            ttobs[i_] = m_data[BW::Q::itt][i_] + m_data[BW::Q::iR][i_] / CGS::c * (1.0 - mu);
        }
        if (t_obs < ttobs[0]) {
            std::cerr << AT << "t_obs=" << t_obs << " < ttobs_arr[0]=" << ttobs[0] << ". Extend ttobs.\n";
            return false;
        }
        if ((t_obs > ttobs[m_i_end_r - 1])) {
            if (m_i_end_r == m_mu.size()-1)
                    std::cerr <<AT << " t_obs="<<t_obs<<" > ttobs_arr[i_end="<<m_i_end_r - 1<<"]="
                              <<ttobs[m_i_end_r - 1]<<". Extend ttobs.\n";
//
//            << " time grid ends too early. "
//                                  << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                  << " while requested obs.time=" << t_obs
//                                  << " extend the grid to later time or request tobs at earlier times\n";
//                    std::cout << ttobs << std::endl;
//            exit(1);
                return false;
        }
        /// locate closest evolution points to the requested obs. time
        ia = findIndex(t_obs, ttobs, ttobs.size());
        if (ia >= m_i_end_r - 1)
            return false; // ??
        ib = ia + 1;
        return true;
        }

    void evalOpticalDepths(double & tau_comp, double & tau_BH, double & tau_bf,
                           size_t ia, size_t ib, double t_obs, double freqprime){

        double rej = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iR]);
        double rho_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iEJrho]);
        double delta_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iEJdelta]);

        /// TODO this freq. should be doppler shifted
        double e_gamma = freqprime * 6.62606957030463E-27;; // Hz -> erg //* 6.62606957030463E-27;//* 4.1356655385381E-15 * CGS::EV_TO_ERG;
        double mu_e = p_pars->mu_e;
        double Z_eff = p_pars->Z_eff;
        int opacitymode = p_pars->opacitymode;
        double albd_fac = p_pars->albd_fac;

        /// --- optical depth due to compton scattering
        double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
        double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
        tau_comp=tau_comp_;
        /// optical depth of BH pair production
        double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
        double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
        tau_BH=tau_BH_;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
        double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
        double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
        tau_bf=tau_bf_;
    }

    /// Fraction of PWN emission that thermalises in ejecta
    double facPSRdep(const double rho_ej, const double delta_ej,
                     const double T_ej, const int opacitymode){
        if (p_pars->curr_b_pwn < 0. || !std::isfinite(p_pars->curr_b_pwn)){
            (*p_log)(LOG_ERR,AT)<<" b_pwn is not set\n";
            exit(1);
        }

        double e_gamma_min;
        double e_gamma_max_tmp;
        double e_gamma_gamma_ani_tmp;
        double e_gamma_syn_b_tmp;

        e_gamma_min = 1.0*EV_TO_ERG; // eV
        e_gamma_max_tmp = PWNradiationMurase::e_gamma_max(p_pars->curr_b_pwn);
        e_gamma_gamma_ani_tmp = PWNradiationMurase::e_gamma_gamma_ani(T_ej);
        e_gamma_syn_b_tmp = PWNradiationMurase::e_gamma_syn_b(
                p_pars->curr_b_pwn,p_pars->gamma_b_w);

        if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
            e_gamma_max_tmp = e_gamma_gamma_ani_tmp;

//        const int i_max = 1000;
//        double e_tmp = e_gamma_min;
//        double del_ln_e = 0.0;

        const double albd_fac = p_pars->albd_fac;
//        int opacitymode = 0;

//        tmp(rho_ej,delta_ej,T_ej,p_pars->albd_fac,opacitymode,e_gamma_max_tmp,e_gamma_syn_b_tmp, e_gamma_min);
//        std::cout << " 11 \n";
        double frac_psr_dep_tmp = 0.0; int i = 0;
//        std::cout << " frac_psr_dep_.size(0" << frac_psr_dep_.size()<<"\n";
//        for (size_t i = 0; i <= p_pars->iterations; i++)
//            frac_psr_dep_tmp += frac_psr_dep_[i];
        int nthreads = 6;
        ///
        if (e_gamma_max_tmp > e_gamma_syn_b_tmp){
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(iters+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=iters;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(
                        e_tmp,p_pars->curr_b_pwn,p_pars->gamma_b_w,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==iters)
                    frac_psr_dep_tmp_tmp = (1.0/2.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i] = frac_psr_dep_tmp_tmp;
            }
        }
        else {
            e_gamma_max_tmp = e_gamma_syn_b_tmp;
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(iters+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=iters;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(
                        e_tmp,p_pars->curr_b_pwn,p_pars->gamma_b_w,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==iters)
                    frac_psr_dep_tmp_tmp = (1.0/3.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i]= frac_psr_dep_tmp_tmp;
            }
        }

        for (size_t i = 0; i <= iters; ++i)
            frac_psr_dep_tmp += frac_psr_dep_[i];

        if (frac_psr_dep_tmp > 1.)
            frac_psr_dep_tmp = 1.;
        if (!std::isfinite(frac_psr_dep_tmp)||frac_psr_dep_tmp < 0){
            (*p_log)(LOG_ERR,AT) << "frac_psr_dep_tmp="<<frac_psr_dep_tmp<<"\n";
            exit(1);
        }
        return frac_psr_dep_tmp;
    }
    double getFacPWNdep( double rho_ej, double delta_ej, double T_ej, double Ye){
        if(!std::isfinite(rho_ej) || rho_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: rho_ej="<<rho_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(delta_ej) || delta_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: delta_ej="<<delta_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(T_ej) || T_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: T_ej="<<T_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(Ye) || Ye < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: Ye="<<Ye<<"\n"; exit(1);
        }

        int opacitymode=0; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
        if (Ye <= 0.2)
            opacitymode = 2; // r-process heavy
        else if ((Ye > 0.2) && (Ye <= 0.3))
            opacitymode = 1; // r-process light
        else if ((Ye > 0.3) && (Ye <= 0.5))
            opacitymode = 0;//0=iron-rich
        else if (Ye > 0.5)
            opacitymode = 3; // CO
        else{
            (*p_log)(LOG_ERR,AT) << " error \n";
            exit(1);
        }
        return facPSRdep(rho_ej, delta_ej, T_ej, opacitymode);
    }

};



/// how to collide two blastwaves
class BlastWaveCollision{
    struct FsolvePars{
        EOSadi * p_eos = nullptr;
        FsolvePars(EOSadi * eos_pointer) : p_eos(eos_pointer){ }
        double m_gam1{}, m_beta1{}, m_mass1{}, m_eint1{}, m_adi1{};
        double m_gam2{}, m_beta2{}, m_mass2{}, m_eint2{}, m_adi2{};
        int getTheMostMassive() const{ return m_mass1 > m_mass2 ? 1 : 2; }
        int getTheFastest() const{ return m_gam1 > m_gam2 ? 1 : 2; }
        int getTheMostEnergetic() const{ return m_eint1 > m_eint2 ? 1 : 2; }
        /// chose which shell to delete: the slowest/less energetic/less massive
        int choseShell(){
            int idx_massive = getTheMostMassive();
            int idx_energetic = getTheMostEnergetic();
            int idx_fastest = getTheFastest();
            if (idx_massive == idx_energetic && idx_massive == idx_fastest){ return idx_massive;}
            if (idx_massive != idx_energetic && idx_massive == idx_fastest){ return idx_massive;}
            if (idx_massive == idx_energetic && idx_massive != idx_fastest){ return idx_fastest;}
            if (idx_massive != idx_energetic && idx_massive != idx_fastest){ return idx_fastest;}
            std::cerr<<AT<<" ERROR\n";
            exit(1);
        }
        void set(double gam1, double beta1, double mass1, double eint1, double adi1,
                 double gam2, double beta2, double mass2, double eint2, double adi2){
            m_gam1=gam1; m_beta1=beta1; m_mass1=mass1; m_eint1=eint1; m_adi1=adi1;
            m_gam2=gam2; m_beta2=beta2; m_mass2=mass2; m_eint2=eint2; m_adi2=adi2;
        }
    };
    FsolvePars * p_colsolve;
    std::unique_ptr<logger> p_log;
    EOSadi * p_eos = nullptr;
public:
    BlastWaveCollision(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BlastWaveCollision");
        p_eos = new EOSadi();
        p_colsolve = new FsolvePars(p_eos);
    }

    ~BlastWaveCollision(){ delete p_colsolve; delete p_eos; }

    void collideBlastWaves(std::unique_ptr<BlastWave> & bw1,
                           std::unique_ptr<BlastWave> & bw2,
                           double * Y, double rcoll, size_t il){
        if (p_colsolve->p_eos == nullptr){
            (*p_log)(LOG_ERR,AT) << " eos pointer is not set\n;";
            exit(1);
        }
        if ((bw1->getPars()->end_evolution) || (bw2->getPars()->end_evolution)){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active\n";
            exit(1);
        }
        /// get relevant parameters
        double gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + SOL::QS::imom] );
        double beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + SOL::QS::imom] );
        double m0_1 = bw1->getPars()->M0;
        double m2_1 = Y[bw1->getPars()->ii_eq + SOL::QS::iM2];
        double m2_1_ = m2_1 * m0_1;
        double eint2_1 = Y[bw1->getPars()->ii_eq + SOL::QS::iEint2];
        double eint2_1_ = eint2_1 * bw1->getPars()->M0 * CGS::c * CGS::c;
        double adi1 = bw1->getEos()->getGammaAdi(gam1, beta1);
        /// --------------------
        double gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + SOL::QS::imom] );
        double beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + SOL::QS::imom] );
        double m0_2 = bw2->getPars()->M0;
        double m2_2 = Y[bw2->getPars()->ii_eq + SOL::QS::iM2];
        double m2_2_ = m2_2 * m0_2;
        double eint2_2 = Y[bw2->getPars()->ii_eq + SOL::QS::iEint2];
        double eint2_2_ = eint2_1 * bw2->getPars()->M0 * CGS::c * CGS::c;
        double adi2 = bw2->getEos()->getGammaAdi(gam2, beta2);
        /// -------------------
        p_colsolve->set(gam1, beta1, m2_1_ + m0_1, eint2_1_, adi1,
                        gam2, beta2, m2_2_ + m0_2, eint2_2_, adi2);
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll=" << il << "] Colliding shells: "
                              << " Masses=(" << p_colsolve->m_mass1 << ", " << p_colsolve->m_mass2 << ")"
                              << " Gammas=(" << p_colsolve->m_gam1 << ", " << p_colsolve->m_gam2 << ")"
                              << " betas=(" << p_colsolve->m_beta1 << ", " << p_colsolve->m_beta2 << ")"
                              << " Eint2=(" << p_colsolve->m_eint1 << ", " << p_colsolve->m_eint2 << ")"
                              << "\n";
        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1 * p_colsolve->m_mass1 + p_colsolve->m_gam2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = p_colsolve->m_eint1 + p_colsolve->m_eint2;
        //(p_colsolve->m_eint1 * p_colsolve->m_mass1 + p_colsolve->m_eint2 * p_colsolve->m_mass2)
        /// (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        /// solve the collision
        int result = solve(i_gM, i_eM, i_mM);
        double m0_c = m0_2 + m0_1;
        double eint2_c = i_eM / (m0_c * CGS::c * CGS::c);
        double m2_c = (1 * (m2_1 + 1) / (m0_2/m0_1 + 1)) + (1 * (m2_2 + 1) / (1 + m0_1/m0_2));
        if (m2_c==1)
            m2_c = m2_1 + m2_2;
        double mom_c = i_gM * EQS::Beta(i_gM);
        /// update the shell composition (mass averaged)
        double ye_c = (bw1->getPars()->Ye0 * bw1->getPars()->M0 + bw2->getPars()->Ye0 * bw2->getPars()->M0)
                      / (bw1->getPars()->M0 + bw2->getPars()->M0);
        if ((ye_c > bw1->getPars()->Ye0) && (ye_c > bw2->getPars()->Ye0)){
            (*p_log)(LOG_ERR,AT)<< "after collision ye="<<ye_c<<" while BWs had: 1 : "
                << bw1->getPars()->Ye0<< " 2 : " << bw2->getPars()->Ye0 << "\n";
            exit(1);
        }
        /// chose which shell to update/delete
        int ish = p_colsolve->choseShell();
        double eint_before,eint_after,m2_before,m2_after,m0_before,m0_after;
        if (ish == 1){
            bw2->getPars()->end_evolution = true;
            bw1->getPars()->M0 = m0_c;
            bw1->getPars()->Ye0 = ye_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::imom] = mom_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::iEint2] = eint2_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::iM2] = m2_c;
        }
        else {
            bw1->getPars()->end_evolution = true;
            bw2->getPars()->M0 = m0_c;
            bw2->getPars()->Ye0 = ye_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::imom] = mom_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::iEint2] = eint2_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::iM2] = m2_c;
        }
        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<il<<"] Outcome for"
                              << " Eint2/M0c^2: ["<<eint2_1<<", "<<eint2_2<<"] -> "<<eint2_c
                              << " M2/M0: ["<<m2_1<<", "<<m2_2<<"] -> "<<m2_c
                              << " M0: ["<<m0_1<<", "<<m0_2<<"] -> "<<m0_c
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }


#if 0
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        double m2_1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2]
                <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        double m0_1 = bw1->getPars()->M0;
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        p_colsolve->m_mass1 = m0_1 + m2_1;

        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        double m2_2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2]
                  <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi2 = bw2->getEos()->getGammaAdi(p_colsolve->m_gam2, p_colsolve->m_beta2);
        /// apply units for the second shell (mass and energy)
        double m0_2 = bw2->getPars()->M0;
        p_colsolve->m_eint2 = p_colsolve->m_eint2 * (bw2->getPars()->M0 * CGS::c * CGS::c);
        p_colsolve->m_mass2 = m2_2 + m0_2;

        /// log the data
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll=" << il << "] Colliding shells: "
//                              << "idx1="<<idx1<<", idx2="<<idx2
                              << " Masses=(" << p_colsolve->m_mass1 << ", " << p_colsolve->m_mass2 << ")"
                              << " Gammas=(" << p_colsolve->m_gam1 << ", " << p_colsolve->m_gam2 << ")"
                              << " betas=(" << p_colsolve->m_beta1 << ", " << p_colsolve->m_beta2 << ")"
                              << " Eint2=(" << p_colsolve->m_eint1 << ", " << p_colsolve->m_eint2 << ")"
                              << "\n";

        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1 * p_colsolve->m_mass1 + p_colsolve->m_gam2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = (p_colsolve->m_eint1 * p_colsolve->m_mass1 + p_colsolve->m_eint2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);

        /// solve the collision
        int result = solve(i_gM, i_eM);
//        i_mM =

        /// chose which shell to update/delete
        int ish = p_colsolve->choseShell();
        size_t iieq; size_t iieq_other;
        double eint_before=0., eint_after=0.;
        double m2_before=0., m2_after=0.;
        double m0_before=0., m0_after=0.;
        if (ish == 1){
            iieq = bw1->getPars()->ii_eq;
            iieq_other = bw2->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            m0_before = bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            //
            bw1->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw1->getPars()->M0;//Y[iieq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::QS::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::QS::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itt] = 0;
            eint_after = i_eM / ( bw1->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::QS::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::QS::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw1->getPars()->M0;
        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            m0_before = bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::QS::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::QS::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itt] = 0;
            eint_after = i_eM / ( bw2->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::QS::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::QS::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw2->getPars()->M0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
//        double _mom = i_gM * EQS::Beta(i_gM);
//        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
//        Y[iieq + DynRadBlastWave::QS::imom] = _mom;
//        Y[iieq + DynRadBlastWave::QS::iEint2] = _eint2; // use ODE units

        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<il<<"] Outcome for"
//                              << " idx1="<<idx1<<", idx2="<<idx2 << " collision:"
                              << " Eint2/M0c^2: "<<eint_before<<" -> "<<eint_after
                              << " M2/M0: "<<m2_before<<" -> "<<m2_after
                              << " M0: "<<m0_before<<" -> "<<m0_after
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }
#endif
    }
private:
    static void func( int n, double x[], double fx[], int &iflag, void * pars ){
        auto pp = (struct FsolvePars *) pars;
        double gM = x[0];
        double eM = x[1];
//        double em = x[1];
        double gAdiM = pp->p_eos->getGammaAdi(gM, EQS::Beta(gM));
        /// total mass conservation
        double mM = pp->m_mass1 + pp->m_mass2;
        /// total energy consercation (Ek + Eint)
        double fx1 = (pp->m_gam1 * pp->m_mass1 + EQS::get_GammaEff(pp->m_gam1, pp->m_adi1) * pp->m_eint1)
                     + (pp->m_gam2 * pp->m_mass2 + EQS::get_GammaEff(pp->m_gam2, pp->m_adi2) * pp->m_eint2)
                     - (gM * mM + EQS::get_GammaEff(gM,gAdiM) * eM );
        /// total momentum conservation
        double fx2 = sqrt(pp->m_gam1 * pp->m_gam1 - 1) * (pp->m_mass1 + pp->m_adi1 * pp->m_eint1 )
                     + sqrt(pp->m_gam2 * pp->m_gam2 - 1) * (pp->m_mass2 + pp->m_adi2 * pp->m_eint2 )
                     - sqrt(gM * gM - 1) * (mM + gAdiM * eM );
//        double fx3 = em
//                     - pp->m_mass1 + pp->m_mass2;
        if (!std::isfinite(fx1) || !std::isfinite(fx2))
            iflag = -1;

        if (!std::isfinite(sqrt(pp->m_gam1 * pp->m_gam1 - 1)))
            iflag = -1;

        fx[0] = fx1;//std::log10(fx1);// * fx1; TODO THis is wrong. It should not be sqred!
        fx[1] = fx2;//std::log10(fx2);// * fx2;
//        fx[2] = fx3;
    };
    int solve(double & iGamma, double & iEint, double &im){
        double *fx;
        int iflag;
        int info;
        int lwa;
        int n = 2;
        double tol = 1e-8;
        double *wa;
        double *x;

        lwa = ( n * ( 3 * n + 13 ) ) / 2;
        fx = new double[n];
        wa = new double[lwa];
        x = new double[n];

//            (*p_log)(LOG_INFO,AT) << "\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve_test2():\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve() solves a nonlinear system of 2 equations.\n";

        x[0] = iGamma;
        x[1] = iEint;
//        x[2] = im;
        iflag = 0;
        func ( n, x, fx, iflag, p_colsolve );
//            r8vec2_print ( n, x, fx, "  Initial X, F(X)" );
        info = fsolve ( func, n, x, p_colsolve, fx, tol, wa, lwa );

        if (info<0){
            (*p_log)(LOG_ERR,AT)<< "\tFsolve failed (try setting 'tol' lower). "
                                   "Using initial guess values. New shell has "
                                <<"Gamma="<<iGamma<<" beta="<<EQS::Beta(iGamma)<<" Eint="<<iEint<<" im="<<im<<"\n";
//                i_gM = x[0];
//                i_eM = x[1];
        }
        else{
            (*p_log)(LOG_INFO,AT)<< "Fsolve successful. New shell has "
                                 <<"Gamma="<<x[0]
                                 <<" beta="<<EQS::Beta(x[0])
                                 <<" Eint="<<x[1]<<"\n";
            iGamma = x[0];
            iEint = x[1];
//            im = x[2];
        }
//            std::cout << "\n";
//            std::cout << "  Returned value of INFO = " << info << "\n";
//            r8vec2_print ( n, x, fx, "  Final X, FX" );
        delete [] fx;
        delete [] wa;
        delete [] x;
        return info;
    };
};


#endif //SRC_BLASTWAVE_H
