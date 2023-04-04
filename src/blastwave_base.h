//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_BASE_H
#define SRC_BLASTWAVE_BASE_H

#include "pch.h"
#include "utils.h"
#include "base_equations.h"
#include "interpolators.h"
#include "ode_solvers.h"
#include "quadratures.h"
#include "rootfinders.h"
#include "observer.h"
#include "synchrotron_an.h"
#include "composition.h"
#include "initial_data.h"

enum METHODS_Up { iuseEint2, iuseGamma }; // energy density behind the shock
enum METHOD_Delta { iuseJoh06, iuseVE12, iNoDelta }; // thickness of the shock
enum METHOD_GammaSh { iuseJK, isameAsGamma, iuseGammaRel, iuseJKwithGammaRel };
enum METHOD_RSh { isameAsR, iuseGammaSh };
enum METHOD_dmdr{ iusingA, iusingdthdR };
enum METHOD_dgdr{ iour, ipeer };

enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };
enum METHOD_NE{ iusenprime, iuseNe };
enum METHODS_RAD { icomovspec, iobservflux };

class ShockMicrophysics{
    std::unique_ptr<logger> p_log;
    std::unique_ptr<SynchrotronAnalytic> p_syn_an;
public:
    ShockMicrophysics(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "ShockMicrophysics");
        p_syn_an = std::make_unique<SynchrotronAnalytic>(loglevel);
    }
    std::unique_ptr<SynchrotronAnalytic> & getAnSynch(){ return p_syn_an; }
};

// -----------------| PARAMETERS |-------------------
struct Pars{
    // set a reference to the data container
    explicit Pars(VecArray & m_data,
                  std::unique_ptr<ShockMicrophysics> & p_fs,
                  std::unique_ptr<ShockMicrophysics> & p_rs,
                  std::unique_ptr<logger> & p_log)
                  : m_data(m_data), p_fs(p_fs), p_rs(p_rs), p_log(p_log) {

    }
    // *************************************** //
    std::unique_ptr<ShockMicrophysics> & p_fs;
    std::unique_ptr<ShockMicrophysics> & p_rs;
    std::unique_ptr<logger> & p_log;
    VecArray & m_data; // data container ( for use in static EATS interators)
    // *************************************** //

    // initial conditions (settings)
    bool is_init = false;
    METHODS_Up m_method_up{};
    METHOD_Delta m_method_Delta{};
    METHOD_GammaSh m_method_gamma_sh{};
    METHOD_RSh m_method_r_sh{};
    METHOD_dmdr m_method_dmdr{};
    METHOD_dgdr m_method_dgdr{};
    LatStruct::METHOD_eats m_method_eats{};
    /// For Magnetar
//    double thickness = -1;
//    double volume = -1;
//    double Rho0 = -1;
//    double kappa0 = -1;
//    double dtau = -1;
//    double tau_to0 = -1;
//    bool is_above_tau1 = false;
    /// blast wave initial conditions
    double M0 = -1.;
    double R0 = -1.;
    double tb0 = -1.;
    double Gamma0 = -1.;
    double mom0 = -1.;
    double E0 = -1.;
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
    // ---
//        double rho0 = -1.;
//        double drhodr0 = -1.;
//        double rho = 0;
//        double drhodr = 0;
    double eps_rad = 0.;
    double dEinjdt = 0.;
    double dEnuc = 0.;
    double dElum = 0.;
    double kappa = 0.;
    // ---
    bool adiabLoss = true;
    // ---
    size_t comp_ix = 0;
    size_t nr = -1;
    size_t ilayer = 0;
    size_t nlayers = 0;
    size_t ishell = 0;
    size_t ii_eq = 0;

    /// --- for work with magnetar
    size_t i_position_in_layer = 0;

    // --- PARS FOR INTERACTION WITH OTHER MODELS
    bool entry = false; // has this BW entered the 'void' left by another
    double entry_time = -1; // time when this BW entered the 'void' left by another
    double entry_r = -1; // same for radius
    // --
    double x=-1;
    // ---
    size_t ijl = 123456;
    size_t prev_ijl = 123456;
//        bool is_fist_entry = false;
    double first_entry_r = -1;
//        double r2 = -2;
    double r_dist = -1;
    bool is_within = false;
    bool is_within0 = false;
    bool is_using_st_prof = false;
//        double Gamma_rel=-1;
//        double Gamma_cbm=-1;
//        double dGamma_rel_dG=-1;
//        double dlnGammaCBMdR=-1;
    bool end_evolution = false;
    bool end_spreading = false;
    double min_beta_terminate = 1.e-8;
    // ---
//        size_t j_i0=123456789;
    double j_theta0=-1;
    double j_theta1=-1;

    double fraction_of_Gamma0_when_spread = -1.;

    /// main switch
//        bool run_jet_bws = false;
//        bool run_ejecta_bws = false;
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

    size_t i0_failed_elecctrons = 0;
    long n_fialed_electrons = 0;
    size_t i0_failed_image = 0;
    long n_failed_images = 0;

    /// -------------------------------

    // set
    bool use_t_e = false;
    double z{}, d_l{}, nu_obs{}, t_obs{}, theta_obs{};
    double (* obsangle)(const double &, const double &, const double &) = nullptr;
    double (* im_xxs)( const double &, const double &, const double & ) = nullptr;
    double (* im_yys)( const double &, const double &, const double & ) = nullptr;
    // --- for adaptive quad. method
    int nmax_phi=-1, nmax_theta=-1;
    double rtol_theta=-1., rtol_phi=-1.;
    double atol_theta=-1., atol_phi=-1.;
    int error = 0;
    double current_phi_hi=-1.;
    double current_phi_low=-1.;
    double current_theta_cone_hi=-1;
    double current_theta_cone_low=-1.;
    double cos_theta_obs=-1.;
    double sin_theta_obs=-1.;
    double phi = -1.;
    double cos_phi = -1.;
    double theta=-1., o_phi=-1, o_theta=-1., o_gam=-1, o_mu=-1, o_r=-1, o_flux=-1, o_theta_j=-1; // for map
    double freq1=-1.0, freq2=-1.0;
    size_t nfreq=0;
    METHODS_QUADRATURES method_quad{};
    METHODS_SHOCK_VEL method_shock_vel{};
    METHOD_NE m_method_ne{};
    METHODS_RAD m_method_rad{};
//        METHOD_TAU method_tau{};
    bool counter_jet = true;
    int spread_method = 1; // TODO this is the only option!
    long nevals = 0; // counter for how many times the integrand was evaluated
    // ---
//        VecArray & m_data;
//        size_t nr =-1;
    Array m_freq_arr{};
    Array m_synch_em{};
    Array m_synch_abs{};
//    std::unique_ptr<SynchrotronAnalytic> & p_syn;
//    std::unique_ptr<logger> p_log;

};
// --------------------------------------------------

/// blast wave system base (class)
class BlastWaveBase{
    std::unique_ptr<logger> p_log;
protected:
    Vector m_tb_arr;
    VecArray m_data{}; // container for the solution of the evolution
    LatSpread * p_spread = nullptr;
    EOSadi * p_eos = nullptr;
    RhoISM * p_dens = nullptr;
    SedovTaylor * p_sedov = nullptr;
    BlandfordMcKee2 * p_bm = nullptr;
    std::unique_ptr<NuclearAtomic> p_nuc = nullptr;
    std::unique_ptr<ShockMicrophysics> p_fs{};
    std::unique_ptr<ShockMicrophysics> p_rs{};
    int m_loglevel;
    Pars * p_pars;
public:
    explicit BlastWaveBase(Vector & tb_arr, size_t ishell, size_t ilayer, int loglevel )
            : m_tb_arr(tb_arr), m_loglevel(loglevel) {
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BlastWaveBase");
        // allocate the space for the the entire solution of a given blast wave
        if (m_tb_arr.size() < 1){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR,AT) << " Time grid is not initialized\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if (m_data.empty()){
            m_data.resize( NVALS );
        }
        if (m_data[Q::itburst].size() < 1) {
            for (auto & arr : m_data) {
                arr.resize( m_tb_arr.size(), 0.0);
            }
        }
        // ---------------------- Methods
//        p_syn = std::make_unique<SynchrotronAnalytic>(loglevel);// SynchrotronAnalytic(loglevel);
        p_fs = std::make_unique<ShockMicrophysics>(loglevel);
        p_rs = std::make_unique<ShockMicrophysics>(loglevel);
//        p_pars->p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "pars");
        p_pars = new Pars(m_data, p_fs, p_rs, p_log); // TODO replace 'new' with std::unique_ptr<>
        p_spread = new LatSpread();
        p_eos = new EOSadi();
        p_dens = new RhoISM(loglevel);
        p_sedov = new SedovTaylor();
        p_bm = new BlandfordMcKee2();
        p_nuc = std::make_unique<NuclearAtomic>(loglevel);
        // ----------------------
        p_pars->nr = m_tb_arr.size();
        p_pars->ilayer = ilayer;
        p_pars->ishell = ishell;
        // ----------------------
    }
    ~BlastWaveBase(){ delete p_pars; delete p_spread; delete p_eos; delete p_dens; delete p_sedov; delete p_bm; }
    // ------------------------------------------------------
    Pars *& getPars(){ return p_pars; }
    EOSadi *& getEos(){ return p_eos; }
    LatSpread *& getSpread(){ return p_spread; }
    RhoISM *& getDensIsm(){ return p_dens; }
    SedovTaylor *& getSedov(){ return p_sedov; }
    BlandfordMcKee2 *& getBM(){ return p_bm; }
    // --------------------------------------------------------
    virtual void updateNucAtomic( double * sol, double t) = 0;
    virtual void setInitConditions( double * arr, size_t i ) = 0;
    virtual bool isSolutionOk( double * sol, size_t i ) = 0;
    virtual bool isToTerminate( double * sol, size_t i ) = 0;
    virtual bool isToStopLateralExpansion( double * sol, size_t i ) = 0;
    virtual void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ) = 0;
    virtual void evaluateRhsDens( double * out_Y, size_t i, double x, double const * Y ) = 0;
    virtual void insertSolution( const double * sol, size_t it, size_t i ) = 0;
    virtual size_t getNeq() = 0;
    virtual void applyUnits( double * sol, size_t i ) = 0;
    // ---------------------------------------------------------
    enum Q {
        // --- properties of the ejecta element
        iRho, idR, iTau, iTauOut, iEinj,
        // -- dynamics ---
        iR, iRsh, irho, idrhodr, iGammaCBM, iGammaREL, idGammaCBMdr, idGammaRELdGamma, iPcbm, idPCBMdrho, iMCBM, iCSCBM,
        iGamma, ibeta, imom, iEint2, iU_p, itheta, ictheta, iErad2, iEsh2, iEad2, iM2,
        itcomov, itburst, itt, ithickness, iadi, irho2, iGammaFsh,
        // --- electrons  ---
        igm, igc, igM, iB, iTheta, iz_cool, ix, inprime, iacc_frac,
        // --- observables ---
        imu,
        // ---
        ijl, ir_dist
    };
    std::vector<std::string> m_vnames{
            // --- properties of the ejecta element
            "Rho", "dR", "Tau", "TauOut", "Einj",
            // --- dynamics ---
            "R", "Rsh", "rho", "drhodr", "GammaRho", "GammaRel", "dGammaRhodr", "dGammaReldGamma", "PCBM", "dPCBMdrho", "MCBM", "CSCBM",
            "Gamma", "beta", "mom", "Eint2", "U_p", "theta", "ctheta", "Erad2", "Esh2", "Ead2", "M2",
            "tcomov", "tburst", "tt", "thickness", "adi", "rho2", "GammaFsh",
            // --- electrons
            "gamma_min", "gamma_c", "gamma_max", "B", "ThetaTherm", "z_cool", "x", "nprime", "accel_frac",
            // --- observables
            "mu",
            // ---
            "ijl", "r_dist"
    };
    static constexpr size_t NVALS = 47; // number of variables to save
    // ---------------------------------------------------------
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
    inline Array & operator[](unsigned ll){ return this->m_data[ll]; }
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
//    inline double etot(dobule Gamma){
//        double etot = (dyn2t2.get(f"Gamma_ej_{ii}") - 1)* cgs.c ** 2 * \
//               (dyn2t2.get(f"M2_ej_{ii}") + 1e48 / (cgs.c ** 2 * dyn2t2.get(f"Gamma_ej_{ii}")[0])) + GammaEff * dyn2t2.get(f"Eint2_ej_{ii}")
//    }
    inline VecArray & getData(){ return m_data; }
    inline Array & getData(Q var){ return m_data[var]; }
    inline double & getVal(Q var, int ix){
        auto ixx = (size_t)ix;
        if (ix == -1) { ixx = m_data[0].size()-1; }
        return m_data[var][ix];
    }
    inline double & getLastVal(Q var){ return m_data[var][p_pars->comp_ix]; }
    void addOtherVars(size_t it){

        m_data[Q::ictheta][it] = ctheta(m_data[Q::itheta][it]);//p_pars->ctheta0 + 0.5 * (2. * m_data[Q::itheta][it] - 2. * p_pars->theta_w);
//        p_dens->evaluateRhoDrhoDr(m_data[Q::iR][it], m_data[Q::ictheta][it]);
        double rho_prev = m_data[Q::irho][it-1];
        double rho = m_data[Q::irho][it];


        /// related to the jet BW density profile
        m_data[Q::irho][it] = p_dens->m_rho_;
        m_data[Q::idrhodr][it] = p_dens->m_drhodr_;
        m_data[Q::iGammaCBM][it] = p_dens->m_GammaRho;
        m_data[Q::iGammaREL][it] = p_dens->m_GammaRel;
        m_data[Q::idGammaCBMdr][it] = p_dens->m_dGammaRhodR;
        m_data[Q::idGammaRELdGamma][it] = p_dens->m_dGammaReldGamma;
        m_data[Q::idPCBMdrho][it] = p_dens->m_dPCBMdrho;
        m_data[Q::iPcbm][it] = p_dens->m_P_cbm;
        m_data[Q::iMCBM][it] = p_dens->m_M_cbm;
        m_data[Q::iCSCBM][it] = p_dens->m_CS_CBM;
        m_data[Q::ijl][it]       = (double)p_pars->ijl;
        m_data[Q::ir_dist][it]   = p_pars->r_dist;
        if (m_data[Q::iGammaREL][it] < 1.){
            m_data[Q::iGammaREL][it] = m_data[Q::iGamma][it];
//            std::cerr << AT << "\n GammaRel="<<m_data[Q::iGammaREL][it]<<"; Exiting...\n";
//            exit(1);
        }

        // other parameters
        m_data[Q::ibeta][it]     = EQS::Beta(m_data[Q::iGamma][it]);

        if (p_pars->end_evolution)
            return;

        if ((it>1)&&(rho_prev < 1e-10 * m_data[Q::irho][it])){
            (*p_log)(LOG_ERR,AT) << " it="<<it<<" density gradient >10 orders of magnitude\n";
//            exit(1);
        }

        m_data[Q::iadi][it]      = p_eos->getGammaAdi(m_data[Q::iGamma][it], // TODO ! is it adi or adi21 (using GammaRel)??
                                                      m_data[Q::ibeta][it]);
        m_data[Q::irho2][it]     = EQS::rho2t(m_data[Q::iGamma][it], // TODO should there be a gammaRel?? with adi43??..
                                              m_data[Q::iadi][it],
                                              m_data[Q::irho][it]);
        /// shock front velocity
        switch (p_pars->m_method_gamma_sh) {

            case iuseJK:
                m_data[Q::iGammaFsh][it] = EQS::GammaSh(m_data[Q::iGamma][it],m_data[Q::iadi][it]);
                break;
            case isameAsGamma:
                m_data[Q::iGammaFsh][it] = m_data[Q::iGamma][it];
                break;
            case iuseGammaRel:
                m_data[Q::iGammaFsh][it] = m_data[Q::iGammaREL][it];
                break;
            case iuseJKwithGammaRel:
                m_data[Q::iGammaFsh][it] = EQS::GammaSh(m_data[Q::iGammaREL][it],m_data[Q::iadi][it]);
                break;
        }
        /// shock front radius
        switch (p_pars->m_method_r_sh) {

            case isameAsR:
                m_data[Q::iRsh][it] = m_data[Q::iR][it]; // overrude
                break;
            case iuseGammaSh:
                break;
        }
        /// shock thickness
        switch (p_pars->m_method_Delta) {

            case iuseJoh06:
                m_data[Q::ithickness][it]= EQS::shock_delta_joh06(m_data[Q::iR][it], m_data[Q::iM2][it],
                                                                  m_data[Q::itheta][it], m_data[Q::iGamma][it],
                                                                  m_data[Q::irho2][it], p_pars->ncells);
                break;
            case iuseVE12:
                m_data[Q::ithickness][it]=EQS::shock_delta(m_data[Q::iR][it],m_data[Q::iGamma][it]);
                break;
            case iNoDelta:
                m_data[Q::ithickness][it] = 1.;
                break;
        }
        /// shock downstream energy density
        switch(p_pars->m_method_up){
            case iuseEint2:
                m_data[Q::iU_p][it] = EQS::get_U_p(m_data[Q::irho2][it],  m_data[Q::iM2][it], m_data[Q::iEint2][it]);
                break;
            case iuseGamma:
                m_data[Q::iU_p][it]= EQS::get_U_p(m_data[Q::irho2][it], m_data[Q::iGammaFsh][it]);
                break;
        }

//        m_data[Q::imom][it] = m_data[Q::imom][it];//m_data[Q::iGamma][it] * m_data[Q::ibeta][it];


        /// compute time in the observer frame at 'i' step // TODO put this in the RHS -- DOne
//        bool use_spread = p_spread->m_method != LatSpread::iNULL;
//        if (m_data[Q::itt][0]==0.)
//            m_data[Q::itt][0] = EQS::init_elapsed_time(
//                    m_data[Q::iR][0], m_data[Q::iGamma][0], use_spread);
//        m_data[Q::itt][it] = m_data[Q::itt][0] + EQS::integrate_elapsed_time(
//                it, m_data[Q::iR], m_data[Q::iGamma], m_data[Q::itheta], use_spread);

        if ( ( m_data[Q::iadi][it] < 1.) || (m_data[Q::iR][it] < 1.) || (m_data[Q::irho2][it] < 0.) ||
             ( m_data[Q::iU_p][it] < 0) || (m_data[Q::iGammaCBM][it] < 1.) ||  (m_data[Q::iGammaFsh][it] < 1.) ) {
            std::cerr << AT << " \n Wrong value at i=" << it << " tb=" << m_data[Q::itburst][it] << "\n"
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
                      << " iR=          " << m_data[Q::iR][it] << "\n"
                      << " itt=         " << m_data[Q::itt][it] << "\n"
                      << " imom=        " << m_data[Q::imom][it] << "\n"
                      << " iGamma=      " << m_data[Q::iGamma][it] << "\n"
                      << " iGammaFsh=   " << m_data[Q::iGammaFsh][it] << "\n"
                      << " iEint2=      " << m_data[Q::iEint2][it] << "\n"
                      << " iEad2=       " << m_data[Q::iEad2][it] << "\n"
                      << " iEsh2=       " << m_data[Q::iEsh2][it] << "\n"
                      << " ictheta=     " << m_data[Q::ictheta][it] << "\n"
                      << " irho=        " << m_data[Q::irho][it] << "\n"
                      << " iEinj=       " << m_data[Q::iEint2][it] << "\n"
                      //                      << " idlnrho1_dr="<< m_data[Q::idlnrho1_dr][it] << "\n"
                      << " idrhodr=     " << m_data[Q::idrhodr][it] << "\n"
                      << " iGammaCBM=   " << m_data[Q::iGammaCBM][it] << "\n"
                      << " iGammaREL=   " << m_data[Q::iGammaREL][it] << "\n"
                      << " idGammaCBMdr=" << m_data[Q::idGammaCBMdr][it] << "\n"
                      << " idGammaRELdGamma="<< m_data[Q::idGammaRELdGamma][it] << "\n"
                      << " ibeta=       "      << m_data[Q::ibeta][it] << "\n"
                      << " iadi=        "       << m_data[Q::iadi][it] << "\n"
                      << " irho2=       "      << m_data[Q::irho2][it] << "\n"
                      << " ithickness=  " << m_data[Q::ithickness][it] << "\n"
                      << " iU_p=        "       << m_data[Q::iU_p][it] << "\n"
                      << " imom=        "       << m_data[Q::imom][it] << "\n";
            exit(1);
        }
    }
//    virtual void evaluateCollision(double * out_Y, size_t i, double x, double const * Y,
//                                   std::unique_ptr<BlastWaveBase> & other, size_t other_i ) = 0;

    virtual void evaluateRhsDensModel2(double * out_Y, size_t i, double x, double const * Y,
                                       void * others, size_t prev_ix) = 0; // std::vector<std::unique_ptr<BlastWaveBase>>
    void setAllParametersForOneLayer(LatStruct & latStruct, //RadBlastWave & bw_obj,
                                     StrDbMap & pars, StrStrMap & opts,
                                     size_t ilayer, size_t ii_eq){



//        if (!std::isfinite(latStruct.dist_E0_pw[ilayer])||
//            !std::isfinite(latStruct.dist_M0_pw[ilayer])||
//            !std::isfinite(latStruct.dist_Mom0_pw[ilayer])){
//            std::cerr << AT << "\n E0="<<latStruct.dist_E0_pw[ilayer]
//                      << " M0=" << latStruct.dist_M0_pw[ilayer]
//                      << " G0=" << latStruct.dist_Mom0_pw[ilayer]
//                      << " wrong value. exiting...\n";
////            std::cerr << AT << " error." << "\n";
//            exit(1);
//        }

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


        // check if unknown option is given
//        check_for_unexpected_opt(opts, listBwOpts(), "listBwOpts()");

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
        p_spread->setPars(a,theta_max,latStruct.m_theta_c,
                                    latStruct.m_theta_w, method_spread);


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
        p_pars->m_method_eats = latStruct.m_method_eats;
        p_pars->nlayers = latStruct.nlayers;

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

        // set initials and costants for the blast wave
        switch (latStruct.m_method_eats) {

            case LatStruct::i_pw:
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
                // double E0, double M0, double Gamma0, double tb0, double theta_a, double theta_b0,
                // double theta_c_l, double theta_c_h, double theta_w, double theta_max, double epsilon_e_rad,
                // size_t ii_eq,
                // double ncells
//                setMagPars(latStruct.dist_E0_pw[ilayer],  = double E0,
//                        latStruct.dist_M0_pw[ilayer], = double M0,
//                        latStruct.dist_G0_pw[ilayer], = double Gamma0,
//                        m_mag_time[0],                  = double tb0
//                        0.,                           = double theta_a
//                        latStruct.m_theta_w,          = double theta_b0,
//                        latStruct.theta_pw[ilayer],   = double theta_c_l,
//                        latStruct.theta_pw[ilayer],   = double theta_c_h,
//                        latStruct.m_theta_w,          = double theta_w,
//                        theta_max,                    = double theta_max,
//                        epsilon_e_rad,                = double epsilon_e_rad,
//                        ii_eq,                        = size_t ii_eq,
//                        (double) latStruct.ncells);   = double ncells
                if (latStruct.m_theta_w > theta_max){
                    (*p_log)(LOG_ERR,AT) << " theta_b0="<<latStruct.m_theta_w<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
//
                p_pars->E0        = latStruct.dist_E0_pw[ilayer];
                p_pars->Ye0       = latStruct.dist_Ye_pw[ilayer];
                p_pars->s0        = latStruct.dist_s_pw[ilayer];
                p_pars->M0        = latStruct.dist_M0_pw[ilayer];
                p_pars->mom0      = latStruct.dist_Mom0_pw[ilayer];
                p_pars->tb0       = m_tb_arr[0];
                p_pars->theta_a   = 0.; // theta_a
                p_pars->theta_b0  = latStruct.m_theta_w; // theta_b0
                p_pars->ctheta0   = 0.5 * (latStruct.theta_pw[ilayer] + latStruct.theta_pw[ilayer]);
//        p_pars->theta_h0= theta_c_h;
                p_pars->theta_c_l = latStruct.theta_pw[ilayer];//theta_c_l;
                p_pars->theta_c_h = latStruct.theta_pw[ilayer];
                p_pars->theta_w   = latStruct.m_theta_w; //
                p_pars->theta_max = theta_max;
                p_pars->ncells    = (double) latStruct.ncells;
                p_pars->eps_rad   = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;


                break;
            case LatStruct::i_adap:
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

                p_pars->E0      = latStruct.dist_E0_a[ilayer];
                p_pars->Ye0     = latStruct.dist_Ye_a[ilayer];
                p_pars->s0      = latStruct.dist_s_a[ilayer];
                p_pars->M0      = latStruct.dist_M0_a[ilayer];
                p_pars->mom0    = latStruct.dist_Mom0_a[ilayer];
                p_pars->tb0     = m_tb_arr[0];
                p_pars->theta_a = 0.;
                p_pars->theta_b0= latStruct.thetas_c_h[ilayer];
                p_pars->ctheta0 = 0.5 * (latStruct.thetas_c_l[ilayer] + latStruct.thetas_c_h[ilayer]);
//        p_pars->theta_h0= theta_c_h;
                p_pars->theta_c_l = latStruct.thetas_c_l[ilayer];
                p_pars->theta_c_h = latStruct.thetas_c_h[ilayer];
                p_pars->theta_w = 0.; //
                p_pars->theta_max = theta_max;
                p_pars->ncells  = 1.;
                p_pars->eps_rad = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;
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

//        std::cout << " ["<<bw_obj.getPars()->ishell<<", "<<bw_obj.getPars()->ilayer<<"] "
//                  <<" G0="<<bw_obj.getPars()->Gamma0
//                  <<" E0="<<bw_obj.getPars()->E0
//                  <<" M0="<<bw_obj.getPars()->M0
//                  <<" ctheta0="<<bw_obj.getPars()->ctheta0
//                  <<"\n";
    }

};

//    void setMagPars(double E0, double M0, double Gamma0, double tb0, double theta_a, double theta_b0,
//                 double theta_c_l, double theta_c_h, double theta_w, double theta_max, double epsilon_e_rad,
//                 size_t ii_eq,
//                 double ncells){
//        if (theta_b0 > theta_max){
//            std::cerr << " theta_b0="<<theta_b0<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
////
//        p_pars->E0      = E0;
//        p_pars->M0      = M0;
//        p_pars->Gamma0  = Gamma0;
//        p_pars->tb0     = tb0;
//        p_pars->theta_a = theta_a;
//        p_pars->theta_b0= theta_b0;
//        p_pars->ctheta0 = 0.5 * (theta_c_l + theta_c_h);
////        p_pars->theta_h0= theta_c_h;
//        p_pars->theta_c_l = theta_c_l;
//        p_pars->theta_c_h = theta_c_h;
//        p_pars->theta_w = theta_w; //
//        p_pars->theta_max = theta_max;
//        p_pars->ncells  = ncells;
//        p_pars->eps_rad = epsilon_e_rad;
//
//        p_spread->m_theta_b0 = p_pars->theta_b0;
//        p_pars->x = p_pars->tb0;
//
//        p_pars->ii_eq  = ii_eq;
//    }

//    inline double theta0(double theta){
//        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        return p_pars->theta_c_l + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//    }
//    inline double theta1(double theta){
//        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        return p_pars->theta_c_h + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//    }

//    virtual void evaluateRhsDensModel2(double * out_Y, size_t i, double x, double const * Y,
//                                       std::vector<std::unique_ptr<RadBlastWave>> & others, size_t prev_ix) = 0;
//    class ISMbehindJet{
//        std::vector<std::unique_ptr<RadBlastWave>> & j_bws;
//        std::unique_ptr<RadBlastWave> & ej_bw;
////        std::unique_ptr<logger> p_log;
//        Array j_cthetas{};
//        Array j_rs{};
//        VecArray m_rho{};
//    public:
//        explicit ISMbehindJet(std::vector<std::unique_ptr<RadBlastWave>> & j_bws,
//                              std::unique_ptr<RadBlastWave> & ej_bw, unsigned loglevel) : j_bws(j_bws), ej_bw(ej_bw) {
////            p_log = std::make_unique<logger>(std::cout, loglevel, "ISMbehindJet");
//            j_cthetas.resize( j_bws.size() );
//            j_rs.resize(j_bws.size() );
//            m_rho.resize( j_bws.size() );
//            for (size_t i = 0; i < j_bws.size(); ++i){
//                m_rho[i].resize( j_bws[i]->getTbGrid().size() );
//            }
//        }
//        void evaluateDensProf(){
//            std::cout  <<" \n" << "evaluating density profile up to it=" << ej_bw->getPars()->comp_ix << "\n";
//            for (size_t i = 0; i < j_bws.size(); ++i){
//                j_cthetas[i] = j_bws[i]->ctheta(j_bws[i]->getLastVal(Q::itheta));
//                j_rs[i] = j_bws[i]->getLastVal(Q::iR);
//                for (size_t j = 0; j < j_bws[i]->getPars()->comp_ix; ++j){
//
//                }
//            }
//        }
//        void setCurrentJet( std::unique_ptr<RadBlastWave> & ej ){
//            for (size_t i = 0; i < j_bws.size(); ++i){
//                j_cthetas[i] = j_bws[i]->ctheta(j_bws[i]->getLastVal(Q::itheta));
//                j_rs[i] = j_bws[i]->getLastVal(Q::iR);
//                for (size_t j = 0; j < ej->getTbGrid().size(); ++j){
//                    //
//
//
//                    //
////                m_rho[i][j] = i_rho;
//                }
//            }
//
//        }
//        void evalDens(RhoISM *& p_dens, double ctheta, double R){
//            for (size_t i = 0; i < j_bws.size(); ++i){
//                if (j_rs[i] <= R){
//                    std::cerr  <<" j_rs[i] <= R \n Exiting..." << "\n";
//                    std::cerr << AT << "\n";
//                    exit(1);
//                }
//
//            }
//        }
//    };
//    RadBlastWave::ISMbehindJet * p_dens_jet = nullptr;
//    void setDensJet( std::vector<std::unique_ptr<RadBlastWave>> & j_bws,
//                     std::unique_ptr<RadBlastWave> & ej_bw ){
//        p_dens_jet = new RadBlastWave::ISMbehindJet(j_bws, ej_bw, m_loglevel);
//    }
//    ~RadBlastWave(){
//        delete p_dens_jet;
//        delete p_eats; // removing EATS_pars for simplicity
//    }


#endif //SRC_BLASTWAVE_BASE_H
