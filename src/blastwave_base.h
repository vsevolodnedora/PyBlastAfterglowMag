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

/// blast wave system base (class)
class BlastWaveBase{
//    std::unique_ptr<logger> p_log;
public:
    enum METHODS_Up { iuseEint2, iuseGamma }; // energy density behind the shock
    enum METHOD_Delta { iuseJoh06, iuseVE12, iNoDelta }; // thickness of the shock
    enum METHOD_GammaSh { iuseJK, isameAsGamma, iuseGammaRel, iuseJKwithGammaRel };
    enum METHOD_RSh { isameAsR, iuseGammaSh };
    enum METHOD_dmdr{ iusingA, iusingdthdR };
    enum METHOD_dgdr{ iour, ipeer };

    enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };
    enum METHOD_NE{ iusenprime, iuseNe };
    enum METHODS_RAD { icomovspec, iobservflux };
protected:
    struct Pars{
        // set a reference to the data container
        explicit Pars(VecArray & m_data) : m_data(m_data) { }
        VecArray & m_data; // data container ( for use in static EATS interators)
        // initial conditions (settings)
        bool is_init = false;
        METHODS_Up m_method_up{};
        METHOD_Delta m_method_Delta{};
        METHOD_GammaSh m_method_gamma_sh{};
        METHOD_RSh m_method_r_sh{};
        METHOD_dmdr m_method_dmdr{};
        METHOD_dgdr m_method_dgdr{};
        LatStruct::METHOD_eats m_method_eats{};
        /// blast wave initial conditions
        double M0 = -1.;
        double R0 = -1.;
        double tb0 = -1.;
        double Gamma0 = -1.;
        double E0 = -1.;
        double theta_a = -1.;
        double theta_b0 = -1.;
        double ncells = -1.;
        double ctheta0 = -1.;
//        double theta_h0 = -1;
        double theta_c_l = -1.;
        double theta_c_h = -1.;
        double theta_w = -1.;
        double theta_max=-1;
        /// deceleration radius
        double Rd = -1;
        // ---
//        double rho0 = -1.;
//        double drhodr0 = -1.;
//        double rho = 0;
//        double drhodr = 0;
        double eps_rad = 0;
        // ---
        bool adiabLoss = true;
        // ---
        size_t comp_ix = 0;
        size_t nr = -1;
        size_t ilayer = 0;
        size_t nlayers = 0;
        size_t ishell = 0;
        size_t ii_eq = 0;
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
        double min_beta_terminate = 1.e-5;
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
        std::unique_ptr<SynchrotronAnalytic> p_syn{};
        std::unique_ptr<logger> p_log;

    };
    Pars * p_pars;
public:
    explicit BlastWaveBase(Array & tb_arr, size_t ishell, size_t ilayer, int loglevel )
            : m_tb_arr(tb_arr), m_loglevel(loglevel) {
//        p_log = std::make_unique<logger>(std::cout, loglevel, "BlastWaveBase");
        // allocate the space for the the entire solution of a given blast wave
        if (m_tb_arr.size() < 1){
            // REMOVING LOGGER
            std::cerr << " Time grid is not initialized\n";
            std::cerr << AT  << "\n";
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
        p_pars = new Pars(m_data);
        p_pars->p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "pars");
        p_spread = new LatSpread();
        p_eos = new EOSadi();
        p_dens = new RhoISM();
        p_sedov = new SedovTaylor();
        p_bm = new BlandfordMcKee2();
        // ----------------------
        p_pars->nr = m_tb_arr.size();
        p_pars->ilayer = ilayer;
        p_pars->ishell = ishell;
        // ----------------------
    }
    ~BlastWaveBase(){ delete p_pars; delete p_spread; delete p_eos; delete p_dens; delete p_sedov; delete p_bm; }
    Pars *& getPars(){ return p_pars; }
    EOSadi *& getEos(){ return p_eos; }
    LatSpread *& getSpread(){ return p_spread; }
    RhoISM *& getDensIsm(){ return p_dens; }
    SedovTaylor *& getSedov(){ return p_sedov; }
    BlandfordMcKee2 *& getBM(){ return p_bm; }
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
    static constexpr size_t NVALS = 43; // number of variables to save
    // ---------------------------------------------------------
    void setPars(double E0, double M0, double Gamma0, double tb0, double theta_a, double theta_b0,
                 double theta_c_l, double theta_c_h, double theta_w, double theta_max, double epsilon_e_rad,
                 size_t ii_eq,
                 double ncells){
        if (theta_b0 > theta_max){
            std::cerr << " theta_b0="<<theta_b0<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
            std::cerr << AT << "\n";
            exit(1);
        }

        p_pars->E0      = E0;
        p_pars->M0      = M0;
        p_pars->Gamma0  = Gamma0;
        p_pars->tb0     = tb0;
        p_pars->theta_a = theta_a;
        p_pars->theta_b0= theta_b0;
        p_pars->ctheta0 = 0.5 * (theta_c_l + theta_c_h);
//        p_pars->theta_h0= theta_c_h;
        p_pars->theta_c_l = theta_c_l;
        p_pars->theta_c_h = theta_c_h;
        p_pars->theta_w = theta_w; //
        p_pars->theta_max = theta_max;
        p_pars->ncells  = ncells;
        p_pars->eps_rad = epsilon_e_rad;

        p_spread->m_theta_b0 = p_pars->theta_b0;
        p_pars->x = p_pars->tb0;

        p_pars->ii_eq  = ii_eq;
    }
    Array & getTbGrid() {return m_tb_arr;}
    Array getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
        Vector tmp{};
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            tmp.push_back(m_tb_arr[it]);
        }
        Array tmp2 (tmp.data(), tmp.size());
        return std::move(tmp2);
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
                std::cerr << " sin(theta= "<<theta<<") is not finite... Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }

            double x2 = fac1*std::sin(theta / 2.);
            double xx2 = 2.*std::asin(x2);
            double x1 = fac0*std::sin(theta / 2.);
            double xx1 = 2.*std::asin(x1);

            ctheta = 0.5 * (xx1 + xx2);
            if (!std::isfinite(ctheta)){
                std::cerr << "ctheta is not finite. ctheta="<<ctheta<<" Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        return ctheta;
    }
//    inline double theta0(double theta){
//        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        return p_pars->theta_c_l + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//    }
//    inline double theta1(double theta){
//        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        return p_pars->theta_c_h + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//    }
    inline Array & getArr(Q var){ return m_data[var]; }
    inline double & getVal(Q var, size_t ix){ return m_data[var][ix]; }
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

        if ((it>1)&&(rho_prev < 1e-10 * rho)){
            std::cerr << AT << " it="<<it<<" density gradient >10 orders of magnitude\n";
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

        m_data[Q::imom][it] = m_data[Q::iGamma][it] * m_data[Q::ibeta][it];


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
                      << " Gamma0="<<p_pars->Gamma0 <<" E0="<<p_pars->E0<<" theta0="<<p_pars->theta_b0
                      << " theta_max="<<p_pars->theta_max <<" ncells="<<p_pars->ncells<< "\n"
                      << " ----------------------------------------------------------- \n"
                      << " iR=          "         << m_data[Q::iR][it] << "\n"
                      << " itt=         "        << m_data[Q::itt][it] << "\n"
                      << " iGamma=      " << m_data[Q::iGamma][it] << "\n"
                      << " iGammaFsh=   " << m_data[Q::iGammaFsh][it] << "\n"
                      << " iEint2=      " << m_data[Q::iEint2][it] << "\n"
                      << " iEad2=       " << m_data[Q::iEad2][it] << "\n"
                      << " iEsh2=       " << m_data[Q::iEsh2][it] << "\n"
                      << " ictheta=     " << m_data[Q::ictheta][it] << "\n"
                      << " irho=        " << m_data[Q::irho][it] << "\n"
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
protected:
    Array m_tb_arr{};
    VecArray m_data{}; // container for the solution of the evolution
    LatSpread * p_spread = nullptr;
    EOSadi * p_eos = nullptr;
    RhoISM * p_dens = nullptr;
    SedovTaylor * p_sedov = nullptr;
    BlandfordMcKee2 * p_bm = nullptr;
    int m_loglevel;
};

#endif //SRC_BLASTWAVE_BASE_H
