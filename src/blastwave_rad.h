//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_RAD_H
#define SRC_BLASTWAVE_RAD_H

#include "pch.h"
#include "utils.h"
#include "base_equations.h"
#include "interpolators.h"
#include "ode_solvers.h"
#include "quadratures.h"
#include "rootfinders.h"
#include "observer.h"
#include "synchrotron_an.h"

#include "blastwave_base.h"

static double check_t_e( double t_e, double mu, double t_obs, Array & mu_arr, int N ) {
    if(mu > mu_arr[N - 1])
    {
        printf("mu >> 1? this should not have happened\n");
        printf("   t_obs=%.6lg t_e=%.6lg mu=%.6lg mu_arr[-1]=%.6lg\n", t_obs, t_e, mu, mu_arr[N - 1]);
        return -1;
    }
    if(mu_arr[0] >= mu) // happens only if t_e very small
    {
        printf("very small mu: mu=%.3lg, mu[0]=%.3lg\n", mu, mu_arr[0]);
        //return t_obs / (1.0 - mu); // so return small t_e limit
        return -1;
    }
    return t_e;
}

/// Blast wave emission from shock
class RadBlastWave : public BlastWaveBase{
    std::unique_ptr<logger> p_log;
public:
    enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };
    enum METHOD_NE{ iusenprime, iuseNe };
    enum METHODS_RAD { icomovspec, iobservflux };

private:
    struct EatsPars{
        explicit EatsPars(VecArray & m_data) : m_data(m_data) { }
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
        VecArray & m_data;
        size_t nr=-1;
        Array m_freq_arr{};
        Array m_synch_em{};
        Array m_synch_abs{};
        std::unique_ptr<SynchrotronAnalytic> p_syn;
//        std::unique_ptr<logger> p_log;
        // ---
    };
    EatsPars * p_eats;

    /// for internal use inside the adaptive quadrature
    void parsPars(double t_obs, double nu_obs,
                  double theta_cone_low, double theta_cone_hi, double phi_low, double phi_hi,
                  double (*obs_angle)( const double &, const double &, const double & )){

        auto & m_data = p_eats->m_data;
        p_eats->nu_obs = nu_obs;
        p_eats->t_obs = t_obs;
        // -- settings
        p_eats->error = 0;
        p_eats->current_phi_hi = phi_hi;
        p_eats->current_phi_low = phi_low;
        p_eats->current_theta_cone_hi = theta_cone_hi;
        p_eats->current_theta_cone_low = theta_cone_low;
        p_eats->cos_theta_obs = cos(p_eats->theta_obs);
        p_eats->sin_theta_obs = sin(p_eats->theta_obs);
        p_eats->obsangle = obs_angle;
        // ---
        p_eats->nr = p_pars->nr;
        for (size_t i = 0; i < p_pars->nr; i++)
            m_data[Q::imu][i] = ( m_data[Q::itburst][i] - t_obs / (1.0 + p_eats->z) ) / m_data[Q::iR][i] * CGS::c;

        if((m_data[Q::imu][p_pars->nr-1] < 1. ) || (m_data[Q::imu][0] > 1.)){
//            std::cout << m_data[Q::itburst][i] << "\n";
//            std::cout << m_data[Q::imu] << "\n";
            std::cerr  << " warning" << " mu[0]=" << m_data[Q::imu][0] << " > 1 (expected <1) or \n";
            std::cerr  << " warning" << " mu[-1]=" << m_data[Q::imu][p_pars->nr-1] << " < 1 (expected >1) \n";
//            exit(1);
        }


        //        std::cerr << AT
//                  << "Observables for Gamma0=[" << (p_pars->p_dyn)->getData()[(*p_pars->p_dyn).Q::iGamma][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iGamma][p_pars->nr - 1] << "] "
//                  << " R0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][p_pars->nr - 1] << "] "
//                  << " time0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][p_pars->nr - 1] << "] "
//                  << " mu0=[" << p_pars->mu_arr[0] << ", " << p_pars->mu_arr[p_pars->nr-1]<<"] "
//                  <<"\n";

    }
    void check_pars(){
        if ((p_eats->nu_obs <= 0.) || (p_eats->z < 0) || (p_eats->z > 10) || (p_eats->d_l < 1) || (p_eats->t_obs < 1)){
            std::cerr << " error in input parameters"
                      << " nu_obs="<<p_eats->nu_obs
                      << " z="<<p_eats->z
                      << " d_l="<<p_eats->d_l
                      << " t_obs="<<p_eats->t_obs
                      << " theta_obs="<<p_eats->theta_obs
                      << " theta_cone_low="<<p_eats->current_theta_cone_low
                      << " theta_cone_hi="<<p_eats->current_theta_cone_hi
                      << " phi_low="<<p_eats->current_phi_low
                      << " phi_hi="<<p_eats->current_phi_hi
                      << "\n Exiting... \n";
            std::cerr << AT << "\n";
            exit(1);
        }
    }

public:
    RadBlastWave(Array & tb_arr, size_t ishell, size_t ilayer, int loglevel )
            : BlastWaveBase(tb_arr, ishell, ilayer, loglevel )
    {
//        p_log = std::make_unique<logger>(std::cout, loglevel, "RadBlastWave");
        p_eats = new EatsPars(m_data);
        p_eats->nr = tb_arr.size();
        p_eats->p_syn = std::make_unique<SynchrotronAnalytic>( loglevel );
//        p_eats->p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EatsPars");
//        m_cspec.resize( cspec_vnames.size() );
//        for(auto & arr : m_cspec)
//            arr.resize( arr.resize(getTbGrid()) )
//        std::unique_ptr<logger>(p_log);
        p_log = std::make_unique<logger>(std::cout, std::cerr, CurrLogLevel, "RadBlastWave");

    }

    static auto listPars(){
        return std::vector<std::string>{
                "theta_obs", "d_l", "z","freq1","freq2","nfreq",
                "nmax_phi","nmax_theta","rtol_theta","rtol_phi","atol_theta","atol_phi" };
    }
    static auto listOpts(){
        return std::vector<std::string>{ "method_quad", "method_shock_vel" "counter_jet", "method_comp_mode" };
    }

    /// set parameters used inside the integrator
    void setEatsPars(StrDbMap & pars, StrStrMap & opts){

        // set parameters
        p_eats->nmax_phi = (int)getDoublePar("nmax_phi", pars, AT, p_log,1000, false);//pars.at("nmax_phi");
        p_eats->nmax_theta = (int)getDoublePar("nmax_theta", pars, AT, p_log,1000, false);//pars.at("nmax_theta");
        p_eats->rtol_theta = getDoublePar("rtol_theta", pars, AT, p_log,1e-2, false);//pars.at("rtol_theta");
        p_eats->rtol_phi = getDoublePar("rtol_phi", pars, AT, p_log,1e-2, false);//pars.at("rtol_phi");
        p_eats->atol_theta = getDoublePar("atol_theta", pars, AT, p_log,1e-2, false);//pars.at("atol_theta");
        p_eats->atol_phi = getDoublePar("atol_phi", pars, AT, p_log,1e-2, false);//pars.at("atol_phi");
        p_eats->theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        p_eats->d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        p_eats->z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");

        p_eats->freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, false);//pars.at("freq1");
        p_eats->freq2 = getDoublePar("freq2", pars, AT, p_log,1.e14, false);//pars.at("freq2");
        p_eats->nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,100, false);//pars.at("nfreq");

        // set options
        std::string opt;

        opt = "method_ne";
        METHOD_NE methodNe;
        if ( opts.find(opt) == opts.end() ) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            methodNe = METHOD_NE::iuseNe;
        }
        else{
            if(opts.at(opt) == "useNe")
                methodNe = METHOD_NE::iuseNe;
            else if(opts.at(opt) == "usenprime")
                methodNe = METHOD_NE::iusenprime;
            else{
                std::cerr << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized \n";
                std::cerr << "Possible options: "
                          << " useNe " << " usenprime " << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_eats->m_method_ne = methodNe;

        opt = "method_quad";
        METHODS_QUADRATURES methodsQuadratures;
        if ( opts.find(opt) == opts.end() ) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
        }
        else{
            if(opts.at(opt) == "CADRE")
                methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
            else if(opts.at(opt) == "TRAP_FIXED")
                methodsQuadratures = METHODS_QUADRATURES::INT_TRAP_FIXED;
            else{
                std::cerr << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized \n";
                std::cerr << "Possible options: "
                          << " CADRE " << " TRAP_FIXED " << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_eats->method_quad = methodsQuadratures;

        opt = "method_shock_vel";
        METHODS_SHOCK_VEL methodsShockVel;
        if ( opts.find(opt) == opts.end() ) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
        }
        else{
            if(opts.at(opt) == "sameAsBW")
                methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
            else if(opts.at(opt) == "shockVel")
                methodsShockVel = METHODS_SHOCK_VEL::ishockVel;
            else{
                std::cerr << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized \n";
                std::cerr << "Possible options: "
                          << " sameAsBW " << " shockVel " << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_eats->method_shock_vel = methodsShockVel;

        p_eats->counter_jet = getBoolOpt("counter_jet", opts, AT, p_log,true);

        opt = "method_comp_mode";
        METHODS_RAD methodCompMode;
        if ( opts.find(opt) == opts.end() ) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            methodCompMode = METHODS_RAD::iobservflux;
        }
        else{
            if(opts.at(opt) == "observFlux")
                methodCompMode = METHODS_RAD::iobservflux;
            else if(opts.at(opt) == "comovSpec")
                methodCompMode = METHODS_RAD::icomovspec;
            else{
                std::cerr << AT << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized \n";
                std::cerr << "Possible options: "
                          << " observFlux " << " comovSpec " << "\n";
                exit(1);
            }
        }
        p_eats->m_method_rad = methodCompMode;


        p_eats->m_freq_arr = TOOLS::MakeLogspace(log10(p_eats->freq1), log10(p_eats->freq2),(int)p_eats->nfreq);
        if (p_eats->m_method_rad == METHODS_RAD::icomovspec){
            std::cout << " allocating comoving spectrum array (fs) "
                      << " freqs="<<p_eats->m_freq_arr.size() << " by radii=" << p_pars->nr << " Spec. grid="
                      << p_eats->m_freq_arr.size() * p_pars->nr << "\n";
            p_eats->m_synch_em.resize( p_eats->m_freq_arr.size() * p_pars->nr );
            p_eats->m_synch_abs.resize( p_eats->m_freq_arr.size() * p_pars->nr );
        }


    }
    void updateObsPars(StrDbMap & pars){
        p_eats->theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        p_eats->d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        p_eats->z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
//        check_for_unexpected_par(pars, {"theta_obs","d_l","z"});
    }

    auto *& getEatsPars() { return p_eats; }
    auto & getSynchAnPtr(){ return p_eats->p_syn; }

    static double shock_synchrotron_flux_density(double Gamma, double GammaShock, double m2, double rho2, double acc_frac, double B,
                                                 double gm, double gM, double gc, double Theta, double z_cool,
                                                 double t_e, double mu, double R, double dr, double dr_tau, void * params){
        auto * p_eats = (struct EatsPars *) params;
        auto & m_data = p_eats->m_data;
        double beta = EQS::Beta(Gamma);
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor

        double beta_shock;
        switch (p_eats->method_shock_vel) {

            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> compute shock velocity
                beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                break;
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
        dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
        dr_tau /= ashock;

        double nuprime = (1.0 + p_eats->z ) * p_eats->nu_obs * delta_D;
        double Ne = m2 / CGS::mp; // numer of protons/electrons
        double nprime = rho2 / CGS::mp; // number density of protons/electrons
        double em_prime,em_lab,abs_prime,abs_lab,intensity,flux_dens,dtau;

        switch (p_eats->m_method_ne) {
            // default (presumably more correct version)
            case iusenprime:
                /// this is informed by van Earten et al. and afterglowpy
                p_eats->p_syn->compute( nprime, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime );
                em_prime = p_eats->p_syn->get_em();
                em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                abs_prime = p_eats->p_syn->get_abs();
                abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_eats->p_syn->getPars()->method_tau);
//                flux_dens = (intensity * R * R * dr) * (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);
                flux_dens = (intensity * R * R * dr) ;
                break;
            case iuseNe:
                /// This is informed by the G. Lamb and Fernandez et al.
                p_eats->p_syn->compute( Ne, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime );
                em_prime = p_eats->p_syn->get_em();
                em_lab = em_prime / (delta_D * delta_D);
                em_lab /= delta_D; // TODO this should be from 'dr'...
                abs_lab = abs_prime * delta_D; // TODO with Ne this might not work as we do not use 'dr' of the shock...
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau,mu,beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_eats->p_syn->getPars()->method_tau);
//                flux_dens = intensity * (1.0 + p_eats->z) / (p_eats->d_l * p_eats->d_l) / 10;
                flux_dens = intensity / 5.; // TODO why no '2'?? // why this need /10 to fit the upper result?
                break;
        }
        return flux_dens;
    }



    /// evaluate flux density using adaptive integrator
    double evalFluxDensA(double t_obs, double nu_obs, double atol) {
        double fluxdens = 0.;

        parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h, 0., M_PI, obsAngle);
        check_pars();
        double Fcoeff = CGS::cgs2mJy / (4.0 * M_PI * p_eats->d_l * p_eats->d_l); // result will be in mJy
        p_eats->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
        p_eats->atol_phi = atol;//  / (2.0 * Fcoeff);
        fluxdens += 2. * Fcoeff * integrate_theta_phi(p_eats); // 2. because Integ_0^pi (not 2pi)
        if (p_eats->counter_jet){
            parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h, 0., M_PI, obsAngleCJ);
            check_pars();
            p_eats->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
            p_eats->atol_phi = atol;//  / (2.0 * Fcoeff);
            fluxdens += 2. * Fcoeff * integrate_theta_phi(p_eats);
        }
        return fluxdens;
    }
    /// check if during the quadrature integration there was an error
    static int check_error(void *params) {
        auto *fp = (struct EatsPars *) params;
        return fp->error;
//        return 0;
    }
    /// find angle at which currently jet ends (after lateral expansion)
    static double find_jet_edge(double phi, double theta_obs, //double cos_th_obs, double sin_th_obs,
                                double th_con_hi, Array & a_mu, Array & a_thj, int N,
                                double (*obs_angle)( const double &, const double &, const double & )) {
        /*
         *
         * Return the (outer) edge of the jet section for a particular obs.
         *
         * phi: double
         *      phi-coordinate of jet along which to search
         * cos_theta_obs: double
         *      cos(theta_obs) cosine of observer angle
         * sin_theta_obs: double
         *      sin(theta_obs) sine of observer angle
         * th_con_hi: double
         *
         * a_mu: double array
         *      time-ordered array of mu values for this observation.
         *      mu = c * (t_em - t_obs) / R(t_em)
         *         = cos(th_obs)*cos(theta) + sin(theta_obs)*sin(theta)*cos(phi)
         * a_thj: double array
         *      time ordered array of jet-edge values. [from dynamic evolution]
         */

//        double cos_phi = cos(phi);
//        double mu = cos(th_con_hi) * cos_th_obs + sin(th_con_hi) * sin_th_obs * cos(phi); // cos(th_obs)*cos(theta) + sin(theta_obs)*sin(theta)*cos(phi)
        double mu = obs_angle( th_con_hi, phi, theta_obs );

        // in the mu_array -- find the index of the current 'mu'
        size_t ia = findIndex(mu, a_mu, N);
        // [ia,ib] contains the boundary of the shell, based on mu

        // *if the jet has not yet spread* = if theta[ia] from dynamic evolution < th_con_hi
        if(a_thj[ia] <= th_con_hi && th_con_hi <= a_thj[ia + 1])
            return th_con_hi;

        // *if the jet started to spread* start with the assumptions

        double theta_a, theta_b;
        // if current angle is larger then th_con_hi
        if(th_con_hi < a_thj[ia]) {
            //The jet is spreading
            theta_a = th_con_hi;   // set boundary
            theta_b = 0.5 * M_PI; // set boundary
        }
        else // if th_con_hi >= a_thj[ia]
        {
            //Guessed too far out!
            theta_a = 0.0;    // set boundary
            theta_b = th_con_hi; // set boundary
        }

        int i = 0;
        // similar to search sorted alborithim :: iteratively approach the 'theta' (theta_a) that correspods to mu
        while(theta_b - theta_a > 1.0e-5 && i < 100)
        {
            double theta_c_i = 0.5 * (theta_a + theta_b);          // middle of the two boundaries of the edge [theta_a ... edge ... ... theta_b]
//            mu = cos(theta_c_i) * cos_th_obs + sin(theta_c_i) * sin_th_obs * cos(phi);  // evaluate 'mu' at the current iteration of the search
            mu = obs_angle(theta_c_i, phi, theta_obs);

            ia = findIndex(mu, a_mu, N);     //
            if(theta_c_i < a_thj[ia])
                theta_a = theta_c_i;
            else
                theta_b = theta_c_i;
            i++;
        }

        //printf("iter: %d, th0=%.6lf, th=%.6lf\n", i, th_con_hi, theta_a);

        return theta_a;
    }
    /// function to be integrated for theta and phi
    static double integrand( double i_cos_theta, double i_phi, void* params ){

        auto * p_eats = (struct EatsPars *) params;
        auto & m_data = p_eats->m_data;
        auto & tburst = m_data[BlastWaveBase::Q::itburst];


        double a_theta = arccos(i_cos_theta);
        double mu = p_eats->obsangle(a_theta, i_phi, p_eats->theta_obs);

        p_eats->theta = a_theta;
        p_eats->o_phi = i_phi;
        p_eats->o_theta = a_theta;

        double dFnu = 0.;

        /// Observed flux density evaluation (interpolate comoving spectrum)
        if (p_eats->m_method_rad == METHODS_RAD::icomovspec){
            Interp2d int_em(p_eats->m_freq_arr, m_data[BlastWaveBase::Q::iR], p_eats->m_synch_em);
            Interp2d int_abs(p_eats->m_freq_arr, m_data[BlastWaveBase::Q::iR], p_eats->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            ///----
            auto & mu_arr = m_data[BlastWaveBase::Q::imu];
            size_t ia = findIndex(mu, m_data[BlastWaveBase::Q::imu], p_eats->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, m_data[BlastWaveBase::Q::imu], m_data[BlastWaveBase::Q::itburst]);
            t_e = check_t_e(t_e, mu, p_eats->t_obs, m_data[BlastWaveBase::Q::imu], (int) p_eats->nr);
            if (t_e < 0.0) {
                // REMOVING LOGGER
                std::cerr << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::iR]);
            //  double r = ( Interp1d(ttobs, m_data[BlastWaveBase::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                std::cerr << " R <= 0. Extend R grid (increasing R0, R1). "
                          << " Current R grid us ["
                          << m_data[BlastWaveBase::Q::iR][0] << ", "
                          << m_data[BlastWaveBase::Q::iR][tburst.size() - 1] << "] "
                          << "and tburst arr ["
                          << tburst[0] << ", " << tburst[p_eats->nr - 1]
                          << "] while the requried obs_time=" << p_eats->t_obs
                          << "\n";
                std::cerr << AT << "\n";
                return 0.;
            }
            double Gamma = interpSegLog(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::iGamma]);
            double beta = interpSegLog(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::ibeta]);
            // double GammaSh = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// compute Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor
            /// compute the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_eats->z ) * p_eats->nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_eats->m_freq_arr, p_eats->m_freq_arr.size());
            size_t ib_nu = ia_nu + 1;
            /// interpolate the emissivity and absorption coefficines
//                double em_prime = int_em.Interpolate(nuprime, r, mth);
            double em_prime = int_em.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
//                double abs_prime = int_abs.Interpolate(nuprime, r, mth);
            double abs_prime = int_abs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
            /// convert to the laboratory frame
            double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
            double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

            /// compute optical depth (for this shock radius and thickness are needed)
            double GammaShock = interpSegLog(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::iGammaFsh]);
            double dr = interpSegLog(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::ithickness]);
            double dr_tau = EQS::shock_delta(r, GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

            double beta_shock;
            switch (p_eats->method_shock_vel) {

                case isameAsBW:
                    beta_shock = EQS::Beta(Gamma);
                    break;
                case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> compute shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
            double intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                               p_eats->p_syn->getPars()->method_tau);
            double flux_dens = (intensity * r * r * dr); //* (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);
            dFnu+=flux_dens;
            /// save the result in image
            double ctheta = interpSegLin(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::ictheta]);
            double theta = interpSegLin(ia, ib, t_e, tburst, m_data[BlastWaveBase::Q::itheta]);

            /// save current position
            p_eats->o_gam = Gamma;
            p_eats->o_r = r;
            p_eats->o_mu = mu;
            p_eats->o_flux = dFnu;
            p_eats->o_theta_j = theta;

        }
            /// Observed flux density evaluation (compute directly)
        else {

            size_t ia = findIndex(mu, m_data[BlastWaveBase::Q::imu], p_eats->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, m_data[BlastWaveBase::Q::imu], m_data[BlastWaveBase::Q::itburst]);
            t_e = check_t_e(t_e, mu, p_eats->t_obs, m_data[BlastWaveBase::Q::imu], (int) p_eats->nr);
            if (t_e < 0.0) {
                // REMOVING LOGGER
                std::cerr << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }

//        double R = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->r_arr, p_eats->nr);
            double R = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iR]);
            if (!std::isfinite(R)) {
                std::cerr << " R is NAN in integrand for radiation" << "\n";
                // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << m_data[BlastWaveBase::Q::iR] << "\n";
//            std::cout << " Gamma= " << m_data[BlastWaveBase::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
                return 0.;
            }

            double rho = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::irho]);
            double Gamma = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                        m_data[BlastWaveBase::Q::iGamma]);
            double GammaSh = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                          m_data[BlastWaveBase::Q::iGammaFsh]);
            double beta = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::ibeta]);
            double U_p = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::iM2));
            double theta = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                        m_data[BlastWaveBase::Q::itheta]);
            double rho2 = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::irho2]);
            double m2 = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iM2]);
            double frac = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                       m_data[BlastWaveBase::Q::iacc_frac]);
            double thick = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                        m_data[BlastWaveBase::Q::ithickness]);
            double gm = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igm]);
            double gM = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igM]);
            double gc = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igc]);
            double B = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iB]);
            double Theta = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                        m_data[BlastWaveBase::Q::iTheta]);
            double z_cool = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst],
                                         m_data[BlastWaveBase::Q::iz_cool]);

            if (rho < 0. || Gamma < 1. || U_p < 0. || theta <= 0. || rho2 < 0. || thick <= 0.) {
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
                std::cout << AT << "\n";
                exit(1);
            }

            dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2,
                                                  frac, B, gm, gM, gc, Theta, z_cool,
                                                  t_e, mu, R, thick, thick, params);
//            dFnu*=(p_eats->d_l*p_eats->d_l*2.);
#if 0
            /* -- Reverse shock --- */
            double dFnu_rs = 0.0;
            if (p_eats->synch_rs){
                std::cout << AT << " Warning! Reverse shock is not finished!" << "\n";
                double rho4   = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::irho4));
                double U_e3   = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::iU_e3));
                double rho3   = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::irho2));
                double thick3 = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::ithickness));
    //            double M3     = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::iM3));
                double gamma43= interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::iGamma43));
    //            double gammaAdi_rs = p_eats->p_eos->getGammaAdi(gamma43, EQS::Beta(gamma43));
    //            double nprime3 = 4.0 * Gamma * (rho4 / CGS::mppme);
    //            double nprime3 = p_eats->eq_rho2(Gamma, rho4 / CGS::mp, gammaAdi_rs); // TODO check if for RS here is Gamma!
                double nprime3 = rho3 / CGS::mp;
                /// compute the 'thickness' of the shock (emitting region)
    //            double dr_rs = thick3;

    //          // TODO check if for the shock velocity gamma43 has to be used!!! I think it is gam43! See gammaAdi calc.
                dFnu_rs = Observables::shock_synchrotron_flux_density(
                        Gamma, gamma43, rho3, U_e3, t_e, mu, R, thick3, thick3, params );
            }
            dFnu += dFnu_rs;
#endif
            if (dFnu == 0 || !std::isfinite(dFnu)) {
                // REMOVING LOGGER
                std::cerr << " flux density is zero ( dFnu = 0 )" << "\n";
            }

            p_eats->o_gam = Gamma;
            p_eats->o_r = R;
            p_eats->o_mu = mu;
            p_eats->o_flux = dFnu;
            p_eats->o_theta_j = theta;

        }

        return dFnu;

    }
    static double costheta_integrand( double aomct, void* params ){
        auto * p_eats = (struct EatsPars *) params;
        p_eats->nevals = p_eats->nevals + 1;
        double act = 1 - aomct; // one minus cos theta 0 -> 'doppler_d' cos theta
        double integ = integrand(act, p_eats->phi, params );
        return integ;

    }
    /// integral of the (costheta_integrand)dtheta
    static double phi_integrand( double a_phi, void* params ){
        double result;
        auto * p_eats = (struct EatsPars *) params;
        auto & m_data = p_eats->m_data;
//        std::cout <<  p_eats->mu_arr[0]  << ' ' << p_eats->mu_arr[100] << "\n";
        p_eats->phi = a_phi;
        p_eats->cos_phi = cos(a_phi);
        // INITIAL GUESSES FOR BOUNDARIES :: TO BE REFINED (based on spreading model)
        double theta_1 = p_eats->current_theta_cone_hi;  // shel upper theta boudnary
        double theta_0 = p_eats->current_theta_cone_low; // shel lower theta boundary
        double Dtheta = theta_1 - theta_0;             // m_size of the theta shell
        /// get the boundaries of the jet
        double th_0, th_1;
        if (p_eats->spread_method == 1){
            th_1 = find_jet_edge(a_phi, p_eats->theta_obs,// p_eats->cos_theta_obs, p_eats->sin_theta_obs,
                                 theta_1, m_data[BlastWaveBase::Q::imu], m_data[BlastWaveBase::Q::itheta],
                                 (int)p_eats->nr,
                                 p_eats->obsangle);
            double frac = theta_0 / theta_1;
            th_0 = frac * th_1;
            theta_0 = th_0;
            theta_1 = th_1;
            // assure that the inner and outer boundaries are within then 0.5Pi steradian
            if(theta_0 > 0.5*M_PI) theta_0 = 0.5*M_PI;
            if(theta_1 > 0.5*M_PI) theta_1 = 0.5*M_PI;
        } else {
            std::cerr  << " method not implemented. Exiting..." << "\n";
            std :: cout << AT << "\n";
            exit(1);
        }
        /// For a given phi, integrate over 1 - cos(theta) (to avoid sin(theta->0) errors)
        double sin_half_theta_0 = sin(0.5 * theta_0);
        double sin_half_theta_1 = sin(0.5 * theta_1);
        double one_minus_cos_theta_0 = 2.0 * sin_half_theta_0 * sin_half_theta_0;
        double one_minus_cos_theta_1 = 2.0 * sin_half_theta_1 * sin_half_theta_1;

        int Neval;
        double err;
        switch (p_eats->method_quad) {

            case INT_TRAP_FIXED:
                result = trap(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              params, check_error);
                break;
            case INT_TRAP_ADAPT:
                result = trap_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                    check_error, nullptr, nullptr);
                break;
            case INT_SIMP_FIXED:
                result = simp(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              params, check_error);
                break;
            case INT_SIMP_ADAPT:
                result = simp_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                    check_error, nullptr, nullptr);
                break;
            case INT_ROMB_ADAPT:
                Neval = 0;
                err = 0;
                result = romb(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              p_eats->atol_theta, p_eats->rtol_theta, params,
                              &Neval, &err, 0, check_error, nullptr, nullptr);
                break;
            case INT_TRAP_NL:
                result = trapNL_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                      check_error, nullptr, nullptr);
                break;
            case INT_HYBRID:
                result = hybrid_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                      check_error, nullptr, nullptr);
                break;
            case INT_CADRE:
                result = cadre_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                     p_eats->nmax_theta, p_eats->atol_theta,
                                     p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                     check_error, nullptr, nullptr);
                break;
            case INT_GK49_ADAPT:
                result = gk49_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                    check_error);
                break;
            case INT_GK715_ADAPT:
                result = gk715_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                     p_eats->nmax_theta, p_eats->atol_theta,
                                     p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                     check_error);
                break;
            case INT_GK1021_ADAPT:
                result = gk1021_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                      check_error);
                break;
        }

        if(result != result || result < 0.0){
            // REMOVING LOGGER
            std::cerr << " phi_integrand failed. Nan or negative. Result=" << result
                      << " for t_obs=" << p_eats->t_obs
                      << " theta_lo=" << theta_0
                      << " theta_hi=" << theta_1
                      << " phi=" << p_eats->phi
                      << "\n"
                      << " Exiting...";
            std::cerr << AT << "\n";
            return 0.;
        }

        return result;

    }
    /// integral of the (phi_integrand)dphi
    static double integrate_theta_phi( void* params ){
        auto * p_eats = (struct EatsPars *) params;
        double atol = p_eats->atol_theta;
        // check if the parameters are set TODO wrap it into "run in safe mode"
        if(p_eats->nmax_theta < 0 || p_eats->nmax_phi < 0 || p_eats->rtol_phi < 0
           || p_eats->rtol_theta < 0 || p_eats->atol_theta < 0 || p_eats->atol_phi < 0 ){
            std::cerr << " one of the tolerance parameters for adaptive quad. is not set (< 0). Given:\n"
                      << " nmax_theta=" << p_eats->nmax_theta
                      << " nmax_phi=" << p_eats->nmax_phi
                      << " rtol_phi=" << p_eats->rtol_phi
                      << " rtol_theta=" << p_eats->rtol_theta
                      << " atol_theta=" << p_eats->atol_theta
                      << " atol_phi=" << p_eats->atol_phi << "\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        double phi_0 = p_eats->current_phi_low;//0.0;
        double phi_1 = p_eats->current_phi_hi;//M_PI;

        double result;
        double phi_a;
        switch (p_eats->method_quad) {

            case INT_TRAP_FIXED:
                result = trap(&phi_integrand, phi_0, phi_1,
                              p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_TRAP_ADAPT:
                result = trap_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_SIMP_FIXED:
                result = simp(&phi_integrand, phi_0, phi_1,
                              p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_SIMP_ADAPT:
                result = simp_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_ROMB_ADAPT:
                phi_a = phi_0 + 0.5*(phi_1-phi_0);
                result = romb(&phi_integrand, phi_0, phi_a,
                              p_eats->nmax_phi, atol,
                              p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                              check_error, nullptr, nullptr);
                result += romb(&phi_integrand, phi_a, phi_1,
                               p_eats->nmax_phi,
                               (atol + p_eats->rtol_phi * result),
                               p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                               check_error, nullptr, nullptr);
                break;
            case INT_TRAP_NL:
                result = trapNL_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_HYBRID:
                result = hybrid_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      0, check_error, nullptr, nullptr);
                break;
            case INT_CADRE:
                result = cadre_adapt(&phi_integrand, phi_0, phi_1,
                                     p_eats->nmax_phi, atol,
                                     p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                     0, check_error, nullptr, nullptr);
                break;
            case INT_GK49_ADAPT:
                result = gk49_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    0, check_error);
                break;
            case INT_GK715_ADAPT:
                result = gk715_adapt(&phi_integrand, phi_0, phi_1,
                                     p_eats->nmax_phi, atol,
                                     p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                     0, check_error);
                break;
            case INT_GK1021_ADAPT:
                result = gk1021_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      0, check_error);
                break;
        }
        return result;
    }



    /// evaluate intensity/flux density distribution using piece-wise summation
    void evalImagePW(Image & image, double obs_time, double obs_freq){
        /// get a pair of images
        Image im_pj;
        Image im_cj;
//        auto tmp = computeImagePW( obs_time, obs_freq );
        computeImagePW(im_pj, im_cj, obs_time, obs_freq );
        /// create empty image with the size two (stacked images)
        Array phi_grid = LatStruct::getCphiGridPW( p_pars->ilayer );
//        Image image( 2 * phi_grid.size(), 0. );
        image.resize( 2 * phi_grid.size(), 0. );
        /// combine the two images (use full 'ncells' array, and fill only cells that correspond to this layer)
        for (size_t icell = 0; icell < phi_grid.size(); ++icell) {
            for (size_t ivn = 0; ivn < image.m_names.size(); ++ivn)
//                image(ivn, icell) = std::get<0>(tmp)(ivn,icell);
                image(ivn, icell) = im_pj(ivn,icell);

            for (size_t ivn = 0; ivn < image.m_names.size(); ++ivn)
//                image(ivn,phi_grid.size()+icell) = std::get<1>(tmp)(ivn,icell);
                image(ivn,phi_grid.size()+icell) = im_cj(ivn,icell);
        }
//        std::cout << image.gerArr(Image::imu) << "\n";
//        image.m_f_tot = std::get<0>(tmp).m_f_tot + std::get<1>(tmp).m_f_tot;
        image.m_f_tot = im_pj.m_f_tot + im_cj.m_f_tot;
//        return std::move( image );
    }
    /// compute the observed flux density distrib 'image' (for a given projection) for given time, freq, angle, distance, red shift
    void evalImageFromPW(Image & image, double t_obs, double nu_obs,
                         double (*obs_angle)( const double &, const double &, const double & ),
                         double (*im_xxs)( const double &, const double &, const double & ),
                         double (*im_yys)( const double &, const double &, const double & )){

//        Image image((size_t) p_pars->ncells);
        image.resize((size_t) p_pars->ncells);

        if ((m_data[BlastWaveBase::Q::iR][0] == 0.) && (p_pars->Gamma0 > 0)){
            std::cerr << " [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                      << " R[0]=0. Seems not evolved -> returning empty image." << "\n";
//            exit(1);
//            return std::move(image);
            return;
        }


        parsPars(t_obs, nu_obs, 0., 0., 0., 0., obs_angle);
        check_pars();

        double flux = 0.;
        auto & m_data = p_eats->m_data;
        if (p_pars->end_evolution){
            std::cerr
                    << "\n [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                    << " Evolution was terminated at ix="<<p_pars->comp_ix<<"\n"
                    << " error might occure here... [TODO] Check if limited calcs to this ix works..\n";
        }
        Array ttobs( 1e200, m_data[BlastWaveBase::Q::iR].size()  );
        Array cphis = LatStruct::getCphiGridPW(p_pars->ilayer);

        /// limit the evaluation to the latest 'R' that is not 0 (before termination)
        size_t nr = m_data[BlastWaveBase::Q::iR].size();
        size_t i_end_r = nr;
        for(size_t ir = 0; ir < nr; ++ir){
            if (m_data[BlastWaveBase::Q::iR][ir] == 0.) {
                i_end_r = ir;
                break;
            }
        }
        if (i_end_r == 0){
            std::cerr << " i_end_r = " << i_end_r << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        /// main loops (if to use pre-computed comoving emissivities or compute them here)
        if (p_eats->m_method_rad == METHODS_RAD::icomovspec) {
            /// init interpolator // TODO major speedup -- do index search for interpolation ONCE and use for both
            Interp2d int_em(p_eats->m_freq_arr, m_data[BlastWaveBase::Q::iR], p_eats->m_synch_em);
            Interp2d int_abs(p_eats->m_freq_arr, m_data[BlastWaveBase::Q::iR], p_eats->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            for (size_t i = 0; i < cphis.size(); i++){
                double phi_cell = cphis[i];
                double ctheta_cell = p_pars->ctheta0;//m_data[BlastWaveBase::Q::ictheta][0]; //cthetas[0];
                double mu = obs_angle(ctheta_cell, phi_cell, p_eats->theta_obs);
                for (size_t i_ = 0; i_ < i_end_r; i_++) {
                    ttobs[i_] = m_data[BlastWaveBase::Q::itt][i_] + m_data[BlastWaveBase::Q::iR][i_] / CGS::c * (1.0 - mu);
                }
                /// check if req. obs time is outside of the evolved times (throw error)
                if (t_obs < ttobs[0]) {
                    std::cerr << " time grid starts too late "
                              << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                              << " extend the grid to earlier time or request tobs at later times\n"
                              << " Exiting...\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }
                if ((i_end_r == nr) && (t_obs > ttobs[i_end_r - 1])) {
                    std::cerr << " time grid ends too early. "
                              << " t_grid[i_end_r-1]=" << ttobs[i_end_r - 1] << " while requested obs.time=" << t_obs << "\n"
                              << " extend the grid to later time or request tobs at earlier times\n"
                              << " Exiting...\n";
                    std::cout << ttobs << std::endl;
                    std::cerr << AT << "\n";
                    exit(1);
                }
                else if ((i_end_r < nr) && (t_obs > ttobs[i_end_r - 1])) {
                    std::cerr << " time grid was shorten to i=" << i_end_r << " from nr=" << nr
                              << " and now ends at t_grid[i_end_r-1]=" << ttobs[i_end_r - 1]
                              << " while t_obs=" << t_obs << "\n";
                    continue;
                }
                /// locate closest evolution points to the requested obs. time
                size_t ia = findIndex(t_obs, ttobs, ttobs.size());
                if (ia >= i_end_r - 1) continue; // ??
                size_t ib = ia + 1;
                /// interpolate the exact radial position of the blast that corresponds to the req. obs time
                double r = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iR]);
                //  double r = ( Interp1d(ttobs, m_data[BlastWaveBase::Q::iR] ) ).Interpolate(t_obs, mth );
                if ((r <= 0.0) || !std::isfinite(r)) {
                    std::cerr << " R <= 0. Extend R grid (increasing R0, R1). "
                              << " Current R grid us ["
                              << m_data[BlastWaveBase::Q::iR][0] << ", "
                              << m_data[BlastWaveBase::Q::iR][m_tb_arr.size() - 1] << "] "
                              << "and tobs arr ["
                              << ttobs[0] << ", " << ttobs[p_pars->nr - 1]
                              << "] while the requried obs_time=" << p_eats->t_obs
                              << "\n";
                    std::cerr << AT << "\n";
                    break;
                }
                double Gamma = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iGamma]);
                double beta = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::ibeta]);
                // double GammaSh = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGammaFsh] ) ).Interpolate(r, mth );
                /// compute Doppler factor
                double a = 1.0 - beta * mu; // beaming factor
                double delta_D = Gamma * a; // doppler factor
                /// compute the comoving obs. frequency from given one in obs. frame
                double nuprime = (1.0 + p_eats->z ) * p_eats->nu_obs * delta_D;
                size_t ia_nu = findIndex(nuprime, p_eats->m_freq_arr, p_eats->m_freq_arr.size());
                size_t ib_nu = ia_nu + 1;
                /// interpolate the emissivity and absorption coefficines
//                double em_prime = int_em.Interpolate(nuprime, r, mth);
                double em_prime = int_em.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
//                double abs_prime = int_abs.Interpolate(nuprime, r, mth);
                double abs_prime = int_abs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                /// convert to the laboratory frame
                double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

                /// compute optical depth (for this shock radius and thickness are needed)
                double GammaShock = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iGammaFsh]);
                double dr = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::ithickness]);
                double dr_tau = EQS::shock_delta(r, GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

                double beta_shock;
                switch (p_eats->method_shock_vel) {

                    case isameAsBW:
                        beta_shock = EQS::Beta(Gamma);
                        break;
                    case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> compute shock velocity
                        beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                        break;
                }
                double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
                dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
                dr_tau /= ashock;
                double dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
                double intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                                   p_eats->p_syn->getPars()->method_tau);
                double flux_dens = (intensity * r * r * dr) * (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);
                flux+=flux_dens;
                /// save the result in image
                double ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::ictheta]);
                //  double ctheta = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::ictheta] ) ).Interpolate(r, mth );
                image(Image::iintens, i) =
                        flux_dens / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
                image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_eats->theta_obs);
                image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_eats->theta_obs);
                image(Image::imu, i) = mu;
            }
            image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
        }
//        else{
//            std::cerr << AT << "\n";
//            exit(1);
//        }
#if 1
        else {
            /// always False
            if(p_eats->use_t_e) {
                /// using t_e & mu
                for (size_t i = 0; i < cphis.size(); i++) {
                    double phi_cell = cphis[i];
                    double ctheta_cell = m_data[BlastWaveBase::Q::ictheta][0]; //cthetas[0];
                    double mu = obs_angle(ctheta_cell, phi_cell, p_eats->theta_obs);

                    size_t ia = findIndex(mu, m_data[BlastWaveBase::Q::imu], p_eats->nr);
                    size_t ib = ia + 1;
                    /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
                    double t_e = interpSegLin(ia, ib, mu, m_data[BlastWaveBase::Q::imu], m_data[BlastWaveBase::Q::itburst]);
                    t_e = check_t_e(t_e, mu, p_eats->t_obs, m_data[BlastWaveBase::Q::imu], (int) p_eats->nr);
                    if (t_e < 0.0) {
                        // REMOVING LOGGER
                        std::cerr << AT << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                        exit(1);
                    }

//        double R = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->r_arr, p_eats->nr);
                    double R = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iR]);
                    if (R == 0.0 || R != R) {
                        // REMOVING LOGGER
                        std::cerr << " nan in data after interpolation: R = " << R
                                  << " R[0:5]: " << m_data[BlastWaveBase::Q::iR][0] << ","
                                  << m_data[BlastWaveBase::Q::iR][1] << ","
                                  << m_data[BlastWaveBase::Q::iR][2] << ","
                                  << m_data[BlastWaveBase::Q::iR][3] << ","
                                  << m_data[BlastWaveBase::Q::iR][4] << ";\t"
                                  << " Gamma[0:5]: "
                                  << m_data[BlastWaveBase::Q::iGamma][0] << ","
                                  << m_data[BlastWaveBase::Q::iGamma][1] << ","
                                  << m_data[BlastWaveBase::Q::iGamma][2] << ","
                                  << m_data[BlastWaveBase::Q::iGamma][3] << ","
                                  << m_data[BlastWaveBase::Q::iGamma][4] << "\n"
                                  <<" Exiting...";
                        std::cerr << AT << "\nR = " << R << "\n";
                        std::cout << "R = " << m_data[BlastWaveBase::Q::iR] << "\n";
                        std::cout << "Gamma= " << m_data[BlastWaveBase::Q::iGamma] << "\n";
                        exit(1);
                    }
                    double r = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iR]);
                    double rs = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iRsh]);
                    double rho = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::irho]);
                    double Gamma = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iGamma]);
                    double GammaSh = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iGammaFsh]);
                    double beta = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::ibeta]);
                    double U_p = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_eats->t_arr_burst, p_eats->dyn(BWDyn::iM2));
                    double theta = interpSegLin(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::itheta]);
                    double ctheta = interpSegLin(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::ictheta]);
                    double rho2 = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::irho2]);
                    double m2 = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iM2]);
                    double frac = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iacc_frac]);
                    double thick = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::ithickness]);
                    double gm = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igm]);
                    double gM = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igM]);
                    double gc = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::igc]);
                    double B = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iB]);
                    double Theta = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iTheta]);
                    double z_cool = interpSegLog(ia, ib, t_e, m_data[BlastWaveBase::Q::itburst], m_data[BlastWaveBase::Q::iz_cool]);

                    double dFnu = 0.;
                    if ((m_data[BlastWaveBase::Q::iB][ia] == 0.) || (m_data[BlastWaveBase::Q::iB][ib] == 0.)){
                        dFnu = 0.;
                    }
                    else{
                        if ((Gamma < 1. || !std::isfinite(Gamma))
                            || (gm < 0.) || (!std::isfinite(gm))
                            || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
                            || (theta <= 0.)
                            || (rho2 < 0.) || (!std::isfinite(rho2))
                            || (thick <= 0.) || (GammaSh <= 1.)) {

                            std::cerr << " wrong value in interpolation to EATS surface "
                                      << " Error in data \n"
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
                            std::cerr << " B != 0 and rho2 is NAN "
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
                            std::cerr << AT << "\n";
                            exit(1);
                        }

                        double thick_tau = EQS::shock_delta(r,GammaSh); // TODO this is added becasue in Johanneson Eq. I use ncells
                        dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                              Theta, z_cool, t_obs, mu,
                                                              r, thick,  thick_tau, p_eats);
                        dFnu *= (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);
                    }

                    flux += dFnu;

                    double mu_ = obs_angle(ctheta, phi_cell, p_eats->theta_obs);
//            if (mu < 1e-2){
//                std::cerr << AT << "\n";
//                exit(1);
//            }
                    image(Image::iintens, i) = dFnu / (r * r * std::abs(mu_)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
//                image(Image::itheta, i) = ctheta;
//                image(Image::itheta_j, i) = theta;
//                image(Image::iphi, i) = phi_cell;
                    image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_eats->theta_obs);
                    image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_eats->theta_obs);
//                image(Image::ir, i) = r;
//                image(Image::igam, i) = Gamma;
                    image(Image::imu, i) = mu_;
//                image(Image::igm, i) = gm;
//                image(Image::iB, i) = B;
//                image(Image::igc, i) = gc;
//                image(Image::itb, i) = t_e;
//                image(Image::itt, i) = 0.;

                }
                image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
//        std::cout << "Total flux " << flux << "\n";
//            return std::move(image);
            }
            else {
                /// using _tobs
                Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
                for (size_t i = 0; i < cphis.size(); i++){
                    double phi_cell = cphis[i];
                    double ctheta_cell = p_pars->ctheta0;//m_data[BlastWaveBase::Q::ictheta][0]; //cthetas[0];
                    double mu = obs_angle(ctheta_cell, phi_cell, p_eats->theta_obs);
//            ttobs = m_data[BlastWaveBase::Q::itt] + m_data[BlastWaveBase::Q::iR] / CGS::c * ( 1.0 - mu );
                    for (size_t i_ = 0; i_ < i_end_r; i_++) {
                        ttobs[i_] = m_data[BlastWaveBase::Q::itt][i_] + m_data[BlastWaveBase::Q::iR][i_] / CGS::c * (1.0 - mu);
                    }
//                print_xy_as_numpy(m_data[BlastWaveBase::Q::iR],m_data[BlastWaveBase::Q::iM2],m_data[BlastWaveBase::Q::itt].size(), 10);
                    if (t_obs < ttobs[0]) {
                        std::cerr << " time grid starts too late "
                                  << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                  << " extend the grid to earlier time or request tobs at later times\n"
                                  << " Exiting...\n";
                        std::cerr << AT << "\n";
                        exit(1);
                    }
                    if ((i_end_r == nr) && (t_obs > ttobs[i_end_r - 1])) {
                        std::cerr << " time grid ends too early. "
                                  << " t_grid[i_end_r-1]=" << ttobs[i_end_r - 1] << " while requested obs.time=" << t_obs << "\n"
                                  << " extend the grid to later time or request tobs at earlier times\n"
                                  << " Exiting...\n";
                        std::cout << ttobs << std::endl;
                        std::cerr << AT << "\n";
                        exit(1);
                    } else if ((i_end_r < nr) && (t_obs > ttobs[i_end_r - 1])) {
                        std::cerr << " time grid was shorten to i=" << i_end_r << " from nr=" << nr
                                  << " and now ends at t_grid[i_end_r-1]=" << ttobs[i_end_r - 1]
                                  << " while t_obs=" << t_obs << "\n";
                        continue;
                    }

                    size_t ia = findIndex(t_obs, ttobs, ttobs.size());
                    if (ia >= i_end_r - 1)
                        continue;
                    size_t ib = ia + 1;


//            double r = ( Interp1d(ttobs, m_data[BlastWaveBase::Q::iR] ) ).Interpolate(t_obs, mth );
//            if( r <= 0.0 || r != r ){
//                std::cerr << " R <= 0. Extend R grid (increasing R0, R1). "
//                          << " Current R grid us ["
//                          << m_data[BlastWaveBase::Q::iR][0] << ", "
//                          << m_data[BlastWaveBase::Q::iR][m_tb_arr.size()-1] << "] "
//                          << "and tobs arr ["
//                          << ttobs[0] << ", "<< ttobs[p_pars->nr-1]
//                          << "] while the requried obs_time=" << p_eats->t_obs
//                          << "\n";
//                break;
//            }
//            double Gamma = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGamma] ) ).Interpolate(r, mth );
//            double GammaSh = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGammaFsh] ) ).Interpolate(r, mth );
//            double rho2 = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::irho2] ) ).Interpolate(r, mth );
//            double thick = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::ithickness] ) ).Interpolate(r, mth );
//            double theta = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itheta] ) ).Interpolate(r, mth );
//            double ctheta = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::ictheta] ) ).Interpolate(r, mth );
//            double B = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iB] ) ).Interpolate(r, mth );
//            double gm = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igm] ) ).Interpolate(r, mth );
//            double gM = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igM] ) ).Interpolate(r, mth );
//            double gc = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igc] ) ).Interpolate(r, mth );
//            double Theta = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iTheta] ) ).Interpolate(r, mth );
//            double z_cool = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iz_cool] ) ).Interpolate(r, mth );
//            double tb = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itburst] ) ).Interpolate(r, mth );
//            double tt = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itt] ) ).Interpolate(r, mth );
//            double m2 = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iM2] ) ).Interpolate(r, mth );
//            double frac = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iacc_frac] ) ).Interpolate(r, mth );
//            double Gamma = ( Interp1d(m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGamma] ) ).Interpolate(r, mth );


//            double r     = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iR] );
                    double r = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iR]);
                    if ((r <= 0.0) || !std::isfinite(r)) {
                        std::cerr << " R <= 0. Extend R grid (increasing R0, R1). "
                                  << " Current R grid us ["
                                  << m_data[BlastWaveBase::Q::iR][0] << ", "
                                  << m_data[BlastWaveBase::Q::iR][m_tb_arr.size() - 1] << "] "
                                  << "and tobs arr ["
                                  << ttobs[0] << ", " << ttobs[p_pars->nr - 1]
                                  << "] while the requried obs_time=" << p_eats->t_obs
                                  << "\n";
                        break;
                    }

//                double Gamma = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGamma]);
//                double GammaSh = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iGammaFsh]);
//                double rho2 = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::irho2]);
//                double m2 = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iM2]);
//                double frac = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iacc_frac]);
//                double thick = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::ithickness]);
//                double theta = interpSegLin(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itheta]);
//                double ctheta = interpSegLin(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::ictheta]);
//                double B = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iB]);
//                double gm = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igm]);
//                double gM = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igM]);
//                double gc = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::igc]);
//                double Theta = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iTheta]);
//                double z_cool = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iz_cool]);
//                double tb = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itburst]);
//                double tt = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::itt]);
//                double cs = interpSegLog(ia, ib, r, m_data[BlastWaveBase::Q::iR], m_data[BlastWaveBase::Q::iCSCBM]);


                    double Gamma = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iGamma]);
                    double GammaSh = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iGammaFsh]);
                    double rho2 = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::irho2]);
                    double m2 = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iM2]);
                    double frac = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iacc_frac]);
                    double thick = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::ithickness]);
                    double theta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::itheta]);
                    double ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::ictheta]);
                    double B = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iB]);
                    double gm = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::igm]);
                    double gM = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::igM]);
                    double gc = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::igc]);
                    double Theta = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iTheta]);
                    double z_cool = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iz_cool]);
                    double tb = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::itburst]);
                    double tt = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::itt]);
                    double cs = interpSegLog(ia, ib, t_obs, ttobs, m_data[BlastWaveBase::Q::iCSCBM]);

                    double dFnu = 0.;
                    if ((m_data[BlastWaveBase::Q::iB][ia] == 0.) || (m_data[BlastWaveBase::Q::iB][ib] == 0.)){
                        dFnu = 0.;
                    }
                    else{
                        if ((Gamma < 1. || !std::isfinite(Gamma))
                            || (gm < 0.) || (!std::isfinite(gm))
                            || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
                            || (theta <= 0.)
                            || (rho2 < 0.) || (!std::isfinite(rho2))
                            || (thick <= 0.) || (GammaSh <= 1.)) {
                            std::cerr<< " Error in interpolation to EATS surface: \n"
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
                            std::cerr<< " B!=0 and rho2 is NAN \n"
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
                            std::cerr << AT << " \n";
                            exit(1);
                        }

                        double thick_tau = EQS::shock_delta(r,GammaSh); // TODO this is added becasue in Johanneson Eq. I use ncells
                        dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                              Theta, z_cool, t_obs, mu,
                                                              r, thick,  thick_tau, p_eats);
                    }
                    dFnu *= (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);
//                if (((Gamma < 1.) || (!std::isfinite(Gamma)))
//                    || (gm < 0.) || (!std::isfinite(gm))
//                    || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
//                    || (theta <= 0.)
//                    || (rho2 < 0.) || (!std::isfinite(rho2))
//                    || (thick <= 0.) || (GammaSh <= 1.)) {
//                    std::cerr << AT << " Error in data \n"
//                              << " R = " << r << "\n"
//                              << " Gamma = " << Gamma << "\n"
//                              << " GammaSh = " << GammaSh << "\n"
//                              << " gm = " << gm << "\n"
//                              << " gM = " << gM << "\n"
//                              << " gc = " << gc << "\n"
//                              << " B = " << B << "\n"
//                              << " Theta = " << Theta << "\n"
//                              << " z_cool = " << z_cool << "\n"
//                              << " theta = " << theta << "\n"
//                              << " rho2 = " << rho2 << "\n"
//                              << " thick = " << thick << "\n"
//                              << " t_obs = " << t_obs << "\n";
//                    exit(1);
//                }
//                if ((B != 0.) && (!std::isfinite(rho2))) {
//                    std::cerr << AT << " Error in data \n"
//                              << " R = " << r << "\n"
//                              << " Gamma = " << Gamma << "\n"
//                              << " gm = " << gm << "\n"
//                              << " gM = " << gM << "\n"
//                              << " gc = " << gc << "\n"
//                              << " B = " << B << "\n"
//                              << " Theta = " << Theta << "\n"
//                              << " z_cool = " << z_cool << "\n"
//                              << " theta = " << theta << "\n"
//                              << " rho2 = " << rho2 << "\n"
//                              << " thick = " << thick << "\n"
//                              << " t_obs = " << t_obs << "\n";
//                    exit(1);
//                }
//
//                // compute total shock thickness for SSA
//                double thick_tau = EQS::shock_delta(r,
//                                                    GammaSh);//r/12./Gamma/Gamma;//m_dyn.pShock()->ShockThickness(r,Gamma,theta,M2,rho2,1.);
//                // if there was no shock
//                double dFnu = 0.;
//                if (B != 0.)
//                    dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc, Theta, z_cool, t_obs,
//                                                            mu, r, thick, thick_tau, p_eats);
////                    dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, B, gm, gM, gc, Theta, z_cool, t_obs,
////                                                          mu, r, thick, thick_tau, p_eats);
////                    Gamma, Gamma, rho2, U_e, t_obs, mu, r, thick, thick_tau, p_eats);

                    flux += dFnu;

//                double mu_ = obs_angle(ctheta, phi_cell, p_eats->theta_obs);
//                if (!std::isfinite(mu_)){
//                    std::cerr << AT << "\n";
//                    exit(1);
//                }
//            if (mu < 1e-2){
//                std::cerr << AT << "\n";
//                exit(1);
//            }
                    image(Image::iintens, i) =
                            dFnu / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
//                image(Image::itheta, i) = ctheta;
//                image(Image::itheta0, i) = ctheta_cell;
//                image(Image::itheta_j, i) = theta;
//                image(Image::iphi, i) = phi_cell;
                    image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_eats->theta_obs);
                    image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_eats->theta_obs);
//                image(Image::ir, i) = r;
//                image(Image::igam, i) = Gamma;
                    image(Image::imu, i) = mu;
//                image(Image::igm, i) = gm;
//                image(Image::iB, i) = B;
//                image(Image::igc, i) = gc;
//                image(Image::itb, i) = tb;
//                image(Image::itt, i) = tt;

                }
                image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
//        std::cout << "Total flux " << flux << "\n";
//            return std::move(image);
            }
        }
#endif
    }
    /// get the observed flux density distrib 'image' for 2 projections for given time, freq, angle, distance, red shift
    void computeImagePW(Image & im_pj, Image & im_cj, double obs_time, double obs_freq){
        size_t cells_in_layer = LatStruct::CellsInLayer(p_pars->ilayer);
        /// compute image for primary jet and counter jet
//        Image im_pj = evalImageFromPW(obs_time, obs_freq, obsAngle, imageXXs, imageYYs);
        evalImageFromPW(im_pj, obs_time, obs_freq, obsAngle, imageXXs, imageYYs);
//        Image im_pj = Image( cells_in_layer, 0);
//        Image im_cj;
//        Image im_pj = Image( cells_in_layer, 0);
        if (p_eats->counter_jet) // p_eats->counter_jet
            evalImageFromPW(im_cj, obs_time, obs_freq, obsAngleCJ, imageXXsCJ, imageYYsCJ);
//            im_cj = evalImageFromPW(im_cj, obs_time, obs_freq, obsAngleCJ, imageXXsCJ, imageYYsCJ);
//        else
//            im_cj = Image( cells_in_layer, 0);
//        im_cj = Image( cells_in_layer, 0);
//        return std::move( std::tuple<Image,Image> (im_pj, im_cj) );
    }


    /// precompute electron distribution properties for each timestep
    void addElectronAnalyticVars(size_t it){
        if ((m_data[Q::ibeta][it] < p_eats->p_syn->getPars()->beta_min)){
            return;
        }
        if ((m_data[Q::iCSCBM][it] >= m_data[Q::ibeta][it])){
            m_data[Q::iB][it] = 0.;
            return;
        }
        // no ISM case -> no electron acceleration
        if ((m_data[Q::iU_p][it] <= 0) || (m_data[Q::irho2][it] <= 0)){
            if (p_pars->i0_failed_elecctrons == 0)
                p_pars->i0_failed_elecctrons = it;
            p_pars->i0_failed_elecctrons += 1;
//            std::cerr << AT << " at it="<<it<<" U_e="<<m_data[Q::iU_e][it]<<" and rho2="<<m_data[Q::irho2][it]<<" skipping electron calc.\n";
            return;
        }
//        if(m_data[Q::iGammaREL][it]<=1.){
//            std::cerr << AT << " at it="<<it<<" iGammaREL="<<m_data[Q::iGammaREL][it]<<" skipping electron calc.\n";
//            exit(1);
//        }
        double Gamma_; // TODO should be GammaSh
        if (m_data[Q::iGammaREL][it] > 0.) Gamma_ = m_data[Q::iGammaREL][it];
        else Gamma_ = m_data[Q::iGamma][it];

        /// for observer frame evaluation, the rad. needs TT, while for comov. spectrum it needs tcomov
        // TODO Check if time is TT or tburst here
//        Array * tmp_time;
//        if (p_eats->p_syn->getPars()->method_comp_mode==SynchrotronAnalytic::METHODS_RAD::icomovspec)
//            *tmp_time = m_data[Q::itt][it];
//        else
//            *tmp_time = m_data[Q::itburst][it];
        p_eats->p_syn->precompute(m_data[Q::iU_p][it],
                                  m_data[Q::iGamma][it],
//                               m_data[Q::iGamma][it],
                                  m_data[Q::iGammaFsh][it],
                                  m_data[Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               m_data[Q::itburst][it], // emission time (TT)
                                  m_data[Q::irho2][it] / CGS::mp);
        // adding
        m_data[Q::iB][it]      = p_eats->p_syn->getPars()->B;
        m_data[Q::igm][it]     = p_eats->p_syn->getPars()->gamma_min;
        m_data[Q::igM][it]     = p_eats->p_syn->getPars()->gamma_max;
        m_data[Q::igc][it]     = p_eats->p_syn->getPars()->gamma_c;
        m_data[Q::iTheta][it]  = p_eats->p_syn->getPars()->Theta;
//        m_data[Q::ix][it]      = p_eats->p_syn->getPars()->x; // IT is not evaluated/ It depends in freq
        m_data[Q::iz_cool][it] = p_eats->p_syn->getPars()->z_cool;
        m_data[Q::inprime][it] = p_eats->p_syn->getPars()->n_prime;
        m_data[Q::iacc_frac][it] = p_eats->p_syn->getPars()->accel_frac;
//        if ( m_data[Q::iGamma][it] < 5 ){
//            exit(1);
//        }
        if ((!std::isfinite( m_data[Q::igm][it] )) || (m_data[Q::iB][it] < 0.)) {
            std::cerr << " Wrong value at i=" << it << " tb=" << m_data[Q::itburst][it]
                      << " iR="     << m_data[Q::iR][it] << "\n"
                      << " iGamma=" << m_data[Q::iGamma][it] << "\n"
                      << " ibeta=" << m_data[Q::ibeta][it] << "\n"
                      << " iM2="    << m_data[Q::iM2][it] << "\n"
                      << " iEint2=" << m_data[Q::iEint2][it] << "\n"
                      << " iU_p="   << m_data[Q::iU_p][it] << "\n"
                      << " irho="   << m_data[Q::irho][it] << "\n"
                      << " irho2="  << m_data[Q::irho2][it] << "\n"
                      << " iB="     << m_data[Q::iB][it] << "\n"
                      << " igm="    << m_data[Q::igm][it] << "\n"
                      << " igM="    << m_data[Q::igM][it] << "\n"
                      << " igc="    << m_data[Q::igc][it] << "\n"
                      << " iTheta=" << m_data[Q::iTheta][it] << "\n"
                      //                      << " ix="     << m_data[Q::ix][it] << "\n"
                      << " iz_cool="<< m_data[Q::iz_cool][it] << "\n"
                      << " inprime="<< m_data[Q::inprime][it]
                      << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
    }
    void computeElectronAnalyticVars(){
        for (size_t it = 0; it < m_tb_arr.size(); ++it){
            addElectronAnalyticVars(it);
        }

        if (p_pars->n_fialed_electrons == m_tb_arr.size()-1){
            std::cerr << "Electron calculation failed for all iterations. Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        else if (p_pars->n_fialed_electrons != 0){
            std::cerr <<" Electron calculation failed for n=" << p_pars->n_fialed_electrons
                      << " iterations starting from it=" << p_pars->i0_failed_elecctrons<<"\n";
//            exit(1);
            std::cerr << AT << "\n";
        }
    }

    /// compute comoving emissivity and absorption for each ifreq and each timestep
    void computeComovingEmissivityAndAbsorption(size_t it){

        /// exit if the obs. radiation method of choice does not need comoving spectrum
        if (p_eats->m_freq_arr.size() < 1){
            std::cerr << " array for comoving spectrum is not initialized \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if (p_eats->m_synch_em.size() < 1){
            std::cerr<< " array for comoving spectrum frequencies is not initialized \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
//        (*p_log)(LOG_INFO) << " computing comoving intensity spectum for "
//                              << p_eats->p_syn->getFreqArr().size() << " grid";
        /// -- check if there are any data first
        double beta_;
        if (m_data[Q::iGammaREL][it] > 0.) beta_ = EQS::Beta( m_data[Q::iGammaREL][it] );
        else beta_ = m_data[Q::ibeta][it];

        /// if BW is not evolved to this 'it' or velocity is smaller than minimum
        if ((m_data[Q::iR][it] < 1)||(beta_ < p_eats->p_syn->getPars()->beta_min))
            return;

        if (m_data[Q::igm][it] == 0){
            std::cerr<< " in evolved blast wave, found gm = 0" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        /// compute emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < p_eats->m_freq_arr.size(); ++ifreq){
            /// compute all types of emissivities and absoprtions
            p_eats->p_syn->compute(m_data[Q::irho2][it] / CGS::mp,//m_data[Q::iM2][it] / CGS::mp,
                                   m_data[Q::iacc_frac][it],
                                   m_data[Q::iB][it],
                                   m_data[Q::igm][it],
                                   m_data[Q::igM][it],
                                   m_data[Q::igc][it],
                                   m_data[Q::iTheta][it],
                                   m_data[Q::iz_cool][it],
//                                       m_data[Q::irho2][it] / CGS::mp,
                                   p_eats->m_freq_arr[ifreq]
            );
            /// add evaluated data to the storage
//            double thick_tau = EQS::shock_delta(m_data[Q::iRsh][it],m_data[Q::iGammaFsh][it]);
//            p_eats->p_syn->addIntensity(thick_tau, 0., 1.);
            p_eats->m_synch_em[ifreq + p_eats->nfreq * it] = p_eats->p_syn->get_em();
            p_eats->m_synch_abs[ifreq + p_eats->nfreq * it] = p_eats->p_syn->get_abs();
        }
    }
    void computeSynchrotronAnalyticSpectrum(){
        if (p_eats->m_method_rad==METHODS_RAD::icomovspec) {
            std::cout << " computing analytic comoving spectrum\n";
            for (size_t it = 0; it < m_tb_arr.size(); ++it) {
                computeComovingEmissivityAndAbsorption(it);
            }
        }
    }

    /// comoving spectrum for external use
    // TODO remove and use previous (above) functions for it!!
    auto evalComovingSynchrotron( Vector & freq_arr, size_t every_it ){

        Array t_arr{};
        if (every_it == 0){
            std::cerr << " comov spectrum at every_it="<<every_it<<" cannot be evaluated.\n Exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        else if (every_it == 1){
//            t_arr = getTbGrid();
//            out_em.resize(m_tb_arr.size());
//            out_abs.resize(m_tb_arr.size());
            t_arr = m_tb_arr;
            for(auto & arr : p_eats->p_syn->m_names_)
                arr.resize( m_tb_arr.size() );
        }
        else{
            t_arr = getTbGrid(every_it);
            for(auto & arr : p_eats->p_syn->m_names_)
                arr.resize( m_tb_arr.size() );
//            out_em.resize(t_arr.size());
//            out_abs.resize(t_arr.size());
        }

        /// allocate momory for data
        std::vector< // types (em/abs)
                std::vector< // freqs
                        std::vector<double>>> // times
        data ( p_eats->p_syn->m_names_.size() );
        for (auto & arr : data) {
            arr.resize(freq_arr.size());
            for (auto & arrr : arr){
                arrr.resize( t_arr.size() );
            }
        }

        /// compute spectra and fill storage
        size_t iit = 0;
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            /// -- check if there are any data first
            double beta_;
            if (m_data[Q::iGammaREL][it] > 0.) beta_ = EQS::Beta( m_data[Q::iGammaREL][it] );
            else beta_ = m_data[Q::ibeta][it];
            if ((m_data[Q::iR][it] < 1)||(beta_ < p_eats->p_syn->getPars()->beta_min))
                continue;
            if (m_data[Q::igm][it] == 0){
                std::cerr << " in evolved blast wave, found gm = 0" << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }

            /// if so
            for (size_t ifreq = 0; ifreq < freq_arr.size(); ++ifreq){
                /// compute all types of emissivities and absoprtions
                p_eats->p_syn->compute(m_data[Q::irho2][it]/CGS::mp,//m_data[Q::iM2][it] / CGS::mp,
                                       m_data[Q::iacc_frac][it],
                                       m_data[Q::iB][it],
                                       m_data[Q::igm][it],
                                       m_data[Q::igM][it],
                                       m_data[Q::igc][it],
                                       m_data[Q::iTheta][it],
                                       m_data[Q::iz_cool][it],
//                                       m_data[Q::irho2][it] / CGS::mp,
                                       freq_arr[ifreq]
                );
                /// add evaluated data to the storage
                double thick_tau = EQS::shock_delta(m_data[Q::iRsh][it],m_data[Q::iGammaFsh][it]);
//                p_eats->p_syn->addIntensity(m_data[Q::ithickness][it], 0., 1.);
                p_eats->p_syn->addIntensity(thick_tau, 0., 1.);
                double val;
                for (size_t i_type = 0; i_type < p_eats->p_syn->m_names_.size(); ++i_type) {
                    val = p_eats->p_syn->getData()[i_type];
                    data[i_type][ifreq][iit] = val;
                }
//                data[p_eats->p_syn->i_int][ifreq][iit] =
//                        (data[p_eats->p_syn->i_em][ifreq][iit] / data[p_eats->p_syn->i_abs][ifreq][iit])
//                        * (1. - )
//                intensity = computeIntensity(em_lab, abs_lab, dr_tau, mu, beta_shock, p_eats->method_tau);


            }
            iit++;
        }
        return std::move(data);
#if 0
        for (size_t ifreq = 0; ifreq < freq_arr.size(); ++ifreq){
            for (size_t it = 0; it < t_arr.size(); ++it){
                // -- check if
                double beta_;
                if (m_data[Q::iGammaREL][it] > 0.) beta_ = EQS::Beta( m_data[Q::iGammaREL][it] );
                else beta_ = m_data[Q::ibeta][it];
                if ((m_data[Q::iR][it] < 1)||(beta_ < p_eats->p_syn->getPars()->beta_min))
                    continue;
            }
        }

        for (size_t i_type = 0; i_type < p_eats->p_syn->cspec_vnames.size(); ++i_type) {
            for ()
        }

        VecVector data_em{}; // emissivity of absorption in the comoving frame
        VecVector data_abs{}; // emissivity of absorption in the comoving frame
        for (auto & freq : freq_arr) {
//            std::cout << "inu=" << freq << "/" << freq_arr.size() << " [comoving spectrum]\n";
            auto data = evalComovingSynchrotron(freq, every_it);
            data_em.emplace_back(std::get<0>(data));
            data_abs.emplace_back(std::get<1>(data));
//            data.emplace_back(evalComovingSynchrotron(freq, itype, every_it));
        }
        return std::move(std::pair(data_em,data_abs));
#endif
    }

private:




    /// compute spectra and fill storage


/* --- ADAPTIVE ---  */
private:



public:
    virtual void evaluateCollision(double * out_Y, size_t i, double x, double const * Y,
                                   std::unique_ptr<RadBlastWave> & other, size_t other_i ) = 0;
    virtual void evaluateRhsDensModel2(double * out_Y, size_t i, double x, double const * Y,
                                       std::vector<std::unique_ptr<RadBlastWave>> & others, size_t prev_ix) = 0;
    class ISMbehindJet{
        std::vector<std::unique_ptr<RadBlastWave>> & j_bws;
        std::unique_ptr<RadBlastWave> & ej_bw;
//        std::unique_ptr<logger> p_log;
        Array j_cthetas{};
        Array j_rs{};
        VecArray m_rho{};
    public:
        explicit ISMbehindJet(std::vector<std::unique_ptr<RadBlastWave>> & j_bws,
                              std::unique_ptr<RadBlastWave> & ej_bw, unsigned loglevel) : j_bws(j_bws), ej_bw(ej_bw) {
//            p_log = std::make_unique<logger>(std::cout, loglevel, "ISMbehindJet");
            j_cthetas.resize( j_bws.size() );
            j_rs.resize(j_bws.size() );
            m_rho.resize( j_bws.size() );
            for (size_t i = 0; i < j_bws.size(); ++i){
                m_rho[i].resize( j_bws[i]->getTbGrid().size() );
            }
        }
        void evaluateDensProf(){
            std::cout  <<" \n" << "evaluating density profile up to it=" << ej_bw->getPars()->comp_ix << "\n";
            for (size_t i = 0; i < j_bws.size(); ++i){
                j_cthetas[i] = j_bws[i]->ctheta(j_bws[i]->getLastVal(Q::itheta));
                j_rs[i] = j_bws[i]->getLastVal(Q::iR);
                for (size_t j = 0; j < j_bws[i]->getPars()->comp_ix; ++j){

                }
            }
        }
        void setCurrentJet( std::unique_ptr<RadBlastWave> & ej ){
            for (size_t i = 0; i < j_bws.size(); ++i){
                j_cthetas[i] = j_bws[i]->ctheta(j_bws[i]->getLastVal(Q::itheta));
                j_rs[i] = j_bws[i]->getLastVal(Q::iR);
                for (size_t j = 0; j < ej->getTbGrid().size(); ++j){
                    //


                    //
//                m_rho[i][j] = i_rho;
                }
            }

        }
        void evalDens(RhoISM *& p_dens, double ctheta, double R){
            for (size_t i = 0; i < j_bws.size(); ++i){
                if (j_rs[i] <= R){
                    std::cerr  <<" j_rs[i] <= R \n Exiting..." << "\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }

            }
        }
    };
    RadBlastWave::ISMbehindJet * p_dens_jet = nullptr;
    void setDensJet( std::vector<std::unique_ptr<RadBlastWave>> & j_bws,
                     std::unique_ptr<RadBlastWave> & ej_bw ){
        p_dens_jet = new RadBlastWave::ISMbehindJet(j_bws, ej_bw, m_loglevel);
    }
    ~RadBlastWave(){ delete p_dens_jet; delete p_eats; }

#if 0
    //    enum COMOV_SPECTRUM_TYPE { iemissivity, iabsorption, iem_th, iem_pl, iabs_th, iabs_pl };
    void evalComovingSynchrotron( double freq, size_t every_it ){
        Array t_arr{};
        Vector out_em{};
        Vector out_abs{};
        if (every_it == 0){
            std::cerr << AT << " comov spectrum at every_it="<<every_it<<" cannot be evaluated.\n";
            exit(1);
        }
        else if (every_it == 1){
//            t_arr = getTbGrid();
//            out_em.resize(m_tb_arr.size());
//            out_abs.resize(m_tb_arr.size());
            for(auto & arr : p_eats->p_syn->m_data)
                arr.resize( m_tb_arr.size() );
        }
        else{
            t_arr = getTbGrid(every_it);
            for(auto & arr : p_eats->p_syn->m_cspec)
                arr.resize( m_tb_arr.size() );
//            out_em.resize(t_arr.size());
//            out_abs.resize(t_arr.size());
        }

        size_t ii = 0;
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){

//        for (size_t it = 0; it < m_tb_arr.size(); it ++){
            /// skip not computed part (after the termination of evolution)
            double beta_;
            if (m_data[Q::iGammaREL][it] > 0.) beta_ = EQS::Beta( m_data[Q::iGammaREL][it] );
            else beta_ = m_data[Q::ibeta][it];
            if ((m_data[Q::iR][it] < 1)||(beta_ < p_eats->p_syn->getPars()->beta_min))
                continue;

            if (m_data[Q::igm][it] == 0){
                std::cerr << AT << " gm = 0" << "\n";
                exit(1);
            }
//            std::cout << "\t it="<<it<<"/"<<m_tb_arr.size()<<" for freq="<<freq<<" [comoving spectrum]\n";
            // evaluate the spectrum
            p_eats->p_syn->compute(m_data[Q::iB][it],
                                   m_data[Q::igm][it],
                                   m_data[Q::igM][it],
                                   m_data[Q::igc][it],
                                   m_data[Q::iTheta][it],
                                   m_data[Q::iz_cool][it],
                                   m_data[Q::irho2][it] / CGS::mp,
                                   freq);
            /// insert computed emissivities/absorptions into resized storage
            for (size_t i = 0; i < p_eats->p_syn->cspec_vnames.size(); i++)
                p_eats->p_syn->m_cspec[i][ii] = p_eats->p_syn->m_data_cspec[i];

//            out_em[ii] = p_eats->p_syn->getPars()->emissivity_prime;
//            out_abs[ii] = p_eats->p_syn->getPars()->absorption_prime;

//            switch (type) {
//                case iemissivity:
//                    out[ii] = p_eats->p_syn->getPars()->emissivity_prime; // total
//                    break;
//                case iabsorption:
//                    out[ii] = p_eats->p_syn->getPars()->absorption_prime;
//                    break;
//                case iem_th:
//                    out[ii] = p_eats->p_syn->getPars()->em_th; // fpr Margalit model
//                    break;
//                case iem_pl:
//                    out[ii] = p_eats->p_syn->getPars()->em_pl;
//                    break;
//                case iabs_th:
//                    out[ii] = p_eats->p_syn->getPars()->abs_th;
//                    break;
//                case iabs_pl:
//                    out[ii] = p_eats->p_syn->getPars()->abs_pl;
//                    break;
//            }
            ii++;
        }
//        return std::move( std::pair(out_em,out_abs) );
    }
#endif

};

#endif //SRC_BLASTWAVE_RAD_H
