//
// Created by vsevolod on 21/04/23.
//

#ifndef SRC_EATS_H
#define SRC_EATS_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "blastwave_components.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"
//#include "blastwave.h"

enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };
enum METHOD_NE{ iusenprime, iuseNe };
enum METHODS_RAD { icomovspec, iobservflux };

static double check_emission_time( double t_e, double mu, double t_obs, Vector & mu_arr, int N ) {
    if(mu > mu_arr[N - 1]) {

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

/// methods to evaluateShycnhrotronSpectrum radiation from a Blastwave
class EATS{
    struct Pars{
        std::unique_ptr<logger> p_log = nullptr;
        std::unique_ptr<SynchrotronAnalytic> p_syna = nullptr;
        Pars(Vector & tburst, Vector & tt, Vector & r, Vector & theta, Vector & m_gam, Vector & m_bet,
             Vector & freq_arr, Vector & synch_em, Vector & synch_abs,
             size_t i_end_r, size_t ish, size_t il, int loglevel, void * params)
            : m_tburst(tburst), m_tt(tt),  m_r(r), m_theta(theta), m_gam(m_gam), m_bet(m_bet),
              m_freq_arr(freq_arr), m_synch_em(synch_em), m_synch_abs(synch_abs), m_params(params) {

            m_mu.resize(m_tburst.size());
//            m_gam.resize(m_mom.size());
//            m_bet.resize(m_mom.size());
//            for (size_t i = 0; i < m_mom.size(); i++){
//                m_gam[i] = GamFromMom(m_mom[i]);
//                m_bet[i] = BetFromMom(m_mom[i]);
//            }

            ishell = ish;
            ilayer= il;
            m_i_end_r = i_end_r;
//            nr = m_data[BW::Q::itburst].size();
            p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EATS_pars");
            p_syna = std::make_unique<SynchrotronAnalytic>(loglevel);
        }
        void * m_params;
        Vector & m_tburst; Vector & m_tt;
        Vector & m_r; Vector & m_theta; Vector & m_gam; Vector & m_bet;
        Vector & m_freq_arr; Vector & m_synch_em; Vector & m_synch_abs;
        Vector m_mu{};
//        VecVector & m_data;
//        Vector & tb_arr;
        double ctheta0 = -1.;
//        inline Vector & operator[](BW::Q ivn){ return this->m_data[ivn]; }
//        inline double & operator()(BW::Q ivn, size_t ir){ return this->m_data[ivn][ir]; }
        Vector getTbGrid(size_t every_it) {
            if ((every_it == 1)||(every_it==0)) return  m_tburst;
            Vector tmp{};
            for (size_t it = 0; it < nr; it = it + every_it){
                tmp.push_back(m_tburst[it]);
            }
            return std::move(tmp);
        }
        /// -----------------------------------------------------------------------------
        void parsPars(double t_obs_, double nu_obs_,
                      double theta_cone_low, double theta_cone_hi, double phi_low, double phi_hi,
                      double (*obs_angle)( const double &, const double &, const double & )){
//        auto & m_data = p_eats->m_data;
            nu_obs = nu_obs_;
            t_obs = t_obs_;
            // -- settings
            error = 0;
            current_phi_hi = phi_hi;
            current_phi_low = phi_low;
            current_theta_cone_hi = theta_cone_hi;
            current_theta_cone_low = theta_cone_low;
            cos_theta_obs = std::cos(theta_obs);
            sin_theta_obs = std::sin(theta_obs);
            obsangle = obs_angle;
            // ---
//        p_eats->nr = p_pars->nr; // removing EATS_pars for simplicity
            for (size_t i = 0; i < nr; i++)
                m_mu[i] = ( m_tburst[i] - t_obs / (1.0 + z) ) / m_r[i] * CGS::c;

            if(m_mu[nr - 1] < 1. ){
                std::cout << m_tburst << "\n";
                std::cout << m_r << "\n";
                (*p_log)(LOG_WARN,AT) << " mu[-1]=" <<m_mu[nr - 1] << " < 1 (expected >1) "
                                      << " tobs=" << t_obs
                                      << " tbutst[-1]=" << m_tburst[nr - 1]
                                      << " R[nr-1]=" << m_r[nr - 1]
                                      << "\n";
            }
            if(m_mu[0] > 1. ){
                (*p_log)(LOG_WARN,AT) << " mu[0]=" << m_mu[0] << " > 1 (expected <1) "
                                      << " tobs=" << t_obs
                                      << " tbutst[0]=" << m_tburst[0]
                                      << " R[0]=" << m_r[0]
                                      << "\n";
            }


            //        std::cerr << AT
//                  << "Observables for Gamma0=[" << (p_pars->p_dyn)->getData()[(*p_pars->p_dyn).Q::iGamma][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iGamma][p_pars->nr - 1] << "] "
//                  << " R0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][p_pars->nr - 1] << "] "
//                  << " time0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][p_pars->nr - 1] << "] "
//                  << " mu0=[" << p_pars->mu_arr[0] << ", " << p_pars->mu_arr[p_pars->nr-1]<<"] "
//                  <<"\n";

        }
        void check_pars() const {
            if ((nu_obs <= 0.) || (z < 0) || (z > 10) || (d_l < 1) || (t_obs < 1) || (theta_c_h < 0) || (theta_c_l < 0)){
                (*p_log)(LOG_ERR,AT) << " error in input parameters"
                                     << " nu_obs="<<nu_obs
                                     << " z="<<z
                                     << " d_l="<<d_l
                                     << " t_obs="<<t_obs
                                     << " theta_obs="<<theta_obs
                                     << " theta_cone_low="<<current_theta_cone_low
                                     << " theta_cone_hi="<<current_theta_cone_hi
                                     << " phi_low="<<current_phi_low
                                     << " phi_hi="<<current_phi_hi
                                     << " theta_c_h="<<theta_c_h
                                     << " theta_c_l="<<theta_c_l
                                     << "\n Exiting... \n";
                exit(1);
            }
        }
        void setEatsPars(StrDbMap & pars, StrStrMap & opts, size_t n_layers, double ctheta_0,
                         double theta_c_l_, double theta_c_h_, double theta_w_, double theta_max_){
            nlayers=n_layers;
            ctheta0=ctheta_0;
            theta_c_l=theta_c_l_;
            theta_c_h=theta_c_h_;
            theta_w=theta_w_;
            theta_max=theta_max_;

            // set parameters
            nmax_phi = (int)getDoublePar("nmax_phi", pars, AT, p_log,1000, false);//pars.at("nmax_phi");
            nmax_theta = (int)getDoublePar("nmax_theta", pars, AT, p_log,1000, false);//pars.at("nmax_theta");
            rtol_theta = getDoublePar("rtol_theta", pars, AT, p_log,1e-2, false);//pars.at("rtol_theta");
            rtol_phi = getDoublePar("rtol_phi", pars, AT, p_log,1e-2, false);//pars.at("rtol_phi");
            atol_theta = getDoublePar("atol_theta", pars, AT, p_log,1e-2, false);//pars.at("atol_theta");
            atol_phi = getDoublePar("atol_phi", pars, AT, p_log,1e-2, false);//pars.at("atol_phi");
            theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
            d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
            z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");

            freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, false);//pars.at("freq1");
            freq2 = getDoublePar("freq2", pars, AT, p_log,1.e14, false);//pars.at("freq2");
            nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,100, false);//pars.at("nfreq");

            // set options
            std::string opt;

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
            m_method_ne = methodNe;

            opt = "method_quad";
            METHODS_QUADRATURES methodsQuadratures;
            if ( opts.find(opt) == opts.end() ) {
                (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
                methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
            }
            else{
                if(opts.at(opt) == "CADRE")
                    methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
                else if(opts.at(opt) == "TRAP_FIXED")
                    methodsQuadratures = METHODS_QUADRATURES::INT_TRAP_FIXED;
                else{
                    (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                         <<" given: " << opts.at(opt)
                                         << " is not recognized. "
                                         << "Possible options: "
                                         << " CADRE " << " TRAP_FIXED " << "\n";
                    exit(1);
                }
            }
            method_quad = methodsQuadratures;

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
            method_shock_vel = methodsShockVel;

            counter_jet = getBoolOpt("counter_jet", opts, AT, p_log,true, false);

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
            m_method_rad = methodCompMode;

            m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);
            if (m_method_rad == METHODS_RAD::icomovspec){
                (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                                      << " freqs="<<m_freq_arr.size() << " by radii=" << nr << " Spec. grid="
                                      << m_freq_arr.size() * nr << "\n";
                m_synch_em.resize( m_freq_arr.size() * nr );
                m_synch_abs.resize( m_freq_arr.size() * nr );
            }


            /// set synchrotron parameters
            p_syna->setPars(pars, opts);

        }
        void updateObsPars(StrDbMap & pars){
            theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
            d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
            z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
//        check_for_unexpected_par(pars, {"theta_obs","d_l","z"});
        }

        /// -------------------------------------------------------------------------------
        size_t ishell = 0, ilayer = 0;  size_t nlayers = 0; size_t m_i_end_r = 0;
        double theta_c_l = -1.;
        double theta_c_h = -1.;
        double theta_w = -1.;
        double theta_max=-1.;
        bool use_t_e = false;
        double z{}, d_l{}, nu_obs{}, t_obs{}, theta_obs{};
        double (* obsangle)(const double &, const double &, const double &) = nullptr;
        double (* im_xxs)( const double &, const double &, const double & ) = nullptr;
        double (* im_yys)( const double &, const double &, const double & ) = nullptr;
        // ---
        void (* fluxFunc)(
                double & flux_dens, double & r, double & ctheta,
                size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                Vector ttobs, void * params
                ) = nullptr;
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
        size_t nr = 0;
//        Vector m_freq_arr{};
//        Vector m_synch_em{};
//        Vector m_synch_abs{};
        size_t i0_failed_elecctrons = 0;
        long n_fialed_electrons = 0;
    };
    std::unique_ptr<logger> p_log;
    Pars * p_pars{};
private:
    /// for internal use inside the adaptive quadrature
    void parsPars(double t_obs, double nu_obs,
                  double theta_cone_low, double theta_cone_hi, double phi_low, double phi_hi,
                  double (*obs_angle)( const double &, const double &, const double & )){
        p_pars->parsPars(t_obs,nu_obs,theta_cone_low,theta_cone_hi,phi_low,phi_hi,obs_angle);
    }
    void check_pars(){ p_pars->check_pars(); }

public:
    /// ----------------------------------------------------------------------------------------------
    EATS(Vector & tburst, Vector & tt, Vector & r, Vector & theta,Vector & m_gam, Vector & m_bet,
                       Vector & freq_arr, Vector & synch_em, Vector & synch_abs,
                       size_t & i_end_r, size_t ish, size_t il, int loglevel, void * params){
        p_pars = new Pars(tburst,tt,r,theta,m_gam,m_bet,freq_arr,synch_em,synch_abs,
                          i_end_r,ish,il,loglevel,params);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EATS");
    }
    ~EATS(){
        delete p_pars;
    }
    /// ----------------------------------------------------------------------------------------------
    void setEatsPars(StrDbMap & pars, StrStrMap & opts, size_t nlayers, double ctheta0,
                     double theta_c_l, double theta_c_h, double theta_w, double theta_max){
        p_pars->setEatsPars(pars,opts,nlayers,ctheta0,
                            theta_c_l,theta_c_h,theta_w,theta_max);
    }
    void setFluxFunc(void (* fluxFunc)( double & flux_dens, double & r, double & ctheta,
                                        size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                                        Vector ttobs, void * params )){
        p_pars->fluxFunc = fluxFunc;
    }
    void setFluxFuncA(void (* fluxFuncA)(  )){
        p_pars->fluxFuncA = fluxFuncA;
    }
    void updateObsPars(StrDbMap & pars){p_pars->updateObsPars(pars);}
    /// evaluate flux density using adaptive integrator
    double evalFluxDensA(double t_obs, double nu_obs, double atol) {
        double fluxdens = 0.;
        parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                 0., M_PI, obsAngle);
        check_pars();
        double Fcoeff = CGS::cgs2mJy / (4.0 * M_PI * p_pars->d_l * p_pars->d_l); // result will be in mJy
        p_pars->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
        p_pars->atol_phi = atol;//  / (2.0 * Fcoeff);
        fluxdens += 2. * Fcoeff * integrate_theta_phi(p_pars); // 2. because Integ_0^pi (not 2pi)
        if (p_pars->counter_jet){
            parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                     0., M_PI, obsAngleCJ);
            check_pars();
            p_pars->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
            p_pars->atol_phi = atol;//  / (2.0 * Fcoeff);
            fluxdens += 2. * Fcoeff * integrate_theta_phi(p_pars);
        }
        return fluxdens;
    }

    /// evaluate intensity/flux density distribution using piece-wise summation
    void evalImagePW(Image & image, Image & im_pj, Image & im_cj, double obs_time, double obs_freq){
        Vector phi_grid = EjectaID2::getCphiGridPW( p_pars->ilayer );
        computeImagePW(im_pj, im_cj, obs_time, obs_freq );
        /// combine the two images (use full 'ncells' array, and fill only cells that correspond to this layer)
        for (size_t icell = 0; icell < phi_grid.size(); ++icell) {
            for (size_t ivn = 0; ivn < image.m_names.size(); ++ivn)
                image(ivn, icell) = im_pj(ivn,icell);
            for (size_t ivn = 0; ivn < image.m_names.size(); ++ivn)
                image(ivn,phi_grid.size()+icell) = im_cj(ivn,icell);
        }
        image.m_f_tot = (im_pj.m_f_tot + im_cj.m_f_tot);
//        std::cout<<image.m_f_tot<<"\n";

    }


    void evalImageFromPW(Image & image, double t_obs, double nu_obs,
                         double (*obs_angle)( const double &, const double &, const double & ),
                         double (*im_xxs)( const double &, const double &, const double & ),
                         double (*im_yys)( const double &, const double &, const double & )){

//        auto & p_syna = p_pars->p_rs->getAnSynch();
//        auto & m_data = p_pars->m_data;
        if ((p_pars->m_r[0] == 0.) && (p_pars->m_gam[0] == 0.)){
            (*p_log)(LOG_WARN,AT) << " [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                                  << " R[0]=0. Seems not evolved -> returning empty image." << "\n";
            return;
        }

        parsPars(t_obs, nu_obs, 0., 0., 0., 0., obs_angle);
        check_pars();

        double flux = 0.;
//        if (p_pars->end_evolution){
//            (*p_log)(LOG_WARN,AT)
//                    << "[ish=" << p_pars->ishell << ", il="<<p_pars->ilayer << "] "
//                    << " Evolution was terminated at ix="<<p_pars->comp_ix<<" "
//                    << " Error might occure here... [TODO] Check if limited calcs to this ix works..\n";
//        }
        Vector ttobs( p_pars->m_r.size(), std::numeric_limits<double>::max() );
        Vector cphis = EjectaID2::getCphiGridPW(p_pars->ilayer);

//        /// limit the evaluation to the latest 'R' that is not 0 (before termination)
//        size_t nr = p_pars->m_r.size();
//        size_t i_end_r = nr;
//        for(size_t ir = 0; ir < nr; ++ir){
//            if (p_pars->m_r[ir] == 0.) {
//                i_end_r = ir;
//                break;
//            }
//        }
//        if (i_end_r == 0){
//            (*p_log)(LOG_ERR,AT) << " i_end_r = " << i_end_r << "\n";
//            exit(1);
//        }

//        Vector mu_arr (cphis.size(),0.0);
//        for (size_t k = 0; k < cphis.size(); k++) {
//            double _phi_cell = cphis[k];
//            double _ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
//            mu_arr[k] = obs_angle(_ctheta_cell, _phi_cell, p_pars->theta_obs);
//        }

        /// TODO make it so after some 'phi' the loop stops -- speed up the calculation
        for (size_t i = 0; i < cphis.size(); i++) {
            double phi_cell = cphis[i];
            double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
            double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
//                double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (t_obs < ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
            if ((p_pars->m_i_end_r == p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_ERR, AT) << " time grid ends too early. "
                                      << " t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                      << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to later time or request tobs at earlier times\n"
                                      << " Exiting...\n";
//                    std::cout << ttobs << std::endl;
                exit(1);
            }
            else if ((p_pars->m_i_end_r < p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->nr
                                       << " and now ends at t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << t_obs << "\n";
                continue;
            }
            /// locate closest evolution points to the requested obs. time
            size_t ia = findIndex(t_obs, ttobs, ttobs.size());
            if (ia >= p_pars->m_i_end_r - 1) continue; // ??
            size_t ib = ia + 1;
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_obs, ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->nr - 1] << "] "
                                      << "and tobs arr ["
                                      << ttobs[0] << ", " << ttobs[p_pars->nr - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                break;
            }
            /// ----------------------------------------------------
            double flux_dens; double ctheta;
            p_pars->fluxFunc(flux_dens, r, ctheta,
                             ia, ib, mu, t_obs, nu_obs,
                             ttobs, p_pars->m_params);
            /// ----------------------------------------------------
            flux += flux_dens;
            image(Image::iintens, i) =
                    flux_dens / (r * r * std::abs(mu)) *
                    CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
            image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
            image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
            image(Image::imu, i) = mu;
        }
        image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy


#if 0
        /// main loops (if to use pre-computed comoving emissivities or evaluateShycnhrotronSpectrum them here)
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
            /// init interpolator // TODO major speedup -- do index search for interpolation ONCE and use for both
//            Interp2d int_em(p_pars->m_freq_arr, p_pars->m_r, p_pars->m_synch_em);
//            Interp2d int_abs(p_pars->m_freq_arr, p_pars->m_r, p_pars->m_synch_abs);
//            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            for (size_t i = 0; i < cphis.size(); i++){
                double phi_cell = cphis[i];
                double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
                double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
//                double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
                for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                    ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu);
                }
                /// check if req. obs time is outside of the evolved times (throw error)
                if (t_obs < ttobs[0]) {
                    (*p_log)(LOG_ERR,AT) << " time grid starts too late "
                                         << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                         << " extend the grid to earlier time or request tobs at later times\n"
                                         << " Exiting...\n";
                    exit(1);
                }
                if ((p_pars->m_i_end_r == p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                    (*p_log)(LOG_ERR,AT) << " time grid ends too early. "
                                         << " t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                         << " while requested obs.time=" << t_obs << "\n"
                                         << " extend the grid to later time or request tobs at earlier times\n"
                                         << " Exiting...\n";
//                    std::cout << ttobs << std::endl;
                    exit(1);
                }
                else if ((p_pars->m_i_end_r < p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                    (*p_log)(LOG_WARN,AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                          << " from nr=" << p_pars->nr
                                          << " and now ends at t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                          << " while t_obs=" << t_obs << "\n";
                    continue;
                }
                /// locate closest evolution points to the requested obs. time
                size_t ia = findIndex(t_obs, ttobs, ttobs.size());
                if (ia >= p_pars->m_i_end_r - 1) continue; // ??
                size_t ib = ia + 1;
                /// interpolate the exact radial position of the blast that corresponds to the req. obs time
                double r = interpSegLog(ia, ib, t_obs, ttobs, p_pars->m_r);
                //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
                if ((r <= 0.0) || !std::isfinite(r)) {
                    (*p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                         << " Current R grid us ["
                                         << p_pars->m_r[0] << ", "
                                         << p_pars->m_r[p_pars->nr - 1] << "] "
                                         << "and tobs arr ["
                                         << ttobs[0] << ", " << ttobs[p_pars->nr - 1]
                                         << "] while the requried obs_time=" << p_pars->t_obs
                                         << "\n";
                    break;
                }



                double Gamma = interpSegLog(ia, ib, t_obs, ttobs, p_pars->m_gam);
                double beta = interpSegLog(ia, ib, t_obs, ttobs, p_pars->m_bet);
                // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
                /// evaluateShycnhrotronSpectrum Doppler factor
                double a = 1.0 - beta * mu; // beaming factor
                double delta_D = Gamma * a; // doppler factor
                /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
                double nuprime = (1.0 + p_pars->z ) * p_pars->nu_obs * delta_D;
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
                double GammaShock = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iGammaFsh]);
                double dr = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::ithickness]);
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
                                                                   p_pars->p_syna->getPars()->method_tau);
                double flux_dens = (intensity * r * r * dr) * (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                flux += flux_dens;
                /// save the result in image
                double ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BW::Q::ictheta]);
                //  double ctheta = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::ictheta] ) ).Interpolate(r, mth );
                image(Image::iintens, i) =
                        flux_dens / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
                image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
                image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
                image(Image::imu, i) = mu;
            }
            image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
        }
        else {
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            for (size_t i = 0; i < cphis.size(); i++){
                double phi_cell = cphis[i];
                double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
//                double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
                for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                    ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu_arr[i]);
                }
                if (t_obs < ttobs[0]) {
                    (*p_log)(LOG_ERR,AT) << " time grid starts too late "
                                         << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                         << " extend the grid to earlier time or request tobs at later times\n"
                                         << " Exiting...\n";
//                        (*p_log)(LOG_ERR,AT) << AT << "\n";
                    exit(1);
                }
                if ((p_pars->m_i_end_r == p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                    (*p_log)(LOG_ERR,AT) << " time grid ends too early. "
                                         << " t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                         << " while requested obs.time=" << t_obs << "\n"
                                         << " extend the grid to later time or request tobs at earlier times\n"
                                         << " Exiting...\n";
//                        std::cout << ttobs << std::endl;
//                        std::cerr << AT << "\n";
                    exit(1);
                }
                else if ((p_pars->m_i_end_r <p_pars->nr) && (t_obs > ttobs[p_pars->m_i_end_r - 1])) {
                    (*p_log)(LOG_WARN,AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                          << " from nr=" << p_pars->nr
                                          << " and now ends at t_grid[i_end_r-1]=" << ttobs[p_pars->m_i_end_r - 1]
                                          << " while t_obs=" << t_obs << "\n";
                    continue;
                }

                size_t ia = findIndex(t_obs, ttobs, ttobs.size());
                if (ia >= p_pars->m_i_end_r - 1)
                    continue;
                size_t ib = ia + 1 ;

                double r = interpSegLog(ia, ib, t_obs, ttobs, p_pars->m_r);
                if ((r <= 0.0) || !std::isfinite(r)) {
                    (*p_log)(LOG_WARN,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                          << " Current R grid us ["
                                          << p_pars->m_r[0] << ", "
                                          << p_pars->m_r[p_pars->nr - 1] << "] "
                                          << "and tobs arr ["
                                          << ttobs[0] << ", " << ttobs[p_pars->nr - 1]
                                          << "] while the requried obs_time=" << p_pars->t_obs
                                          << "\n";
                    break;
                }




                double Gamma = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iGamma]);
                double GammaSh = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iGammaFsh]);
                double rho2 = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::irho2]);
                double m2 = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iM2]);
                double frac = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iacc_frac]);
                double thick = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::ithickness]);
                double theta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BW::Q::itheta]);
                double ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BW::Q::ictheta]);
                double B = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iB]);
                double gm = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::igm]);
                double gM = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::igM]);
                double gc = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::igc]);
                double Theta = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iTheta]);
                double z_cool = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iz_cool]);
                double tb = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::itburst]);
                double tt = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::itt]);
                double cs = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iCSCBM]);

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
                    dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                          Theta, z_cool, t_obs, mu_arr[i],
                                                          r, thick,  thick_tau, p_pars);
                }
                dFnu *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);

                flux += dFnu;

                image(Image::iintens, i) =
                        dFnu / (r * r * std::abs(mu_arr[i])) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
                image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
                image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
                image(Image::imu, i) = mu_arr[i];
//                image(Image::ir, i) = r;
//                image(Image::igam, i) = Gamma;
//                image(Image::itheta, i) = ctheta_cell;
//                image(Image::itheta0, i) = theta;
//                image(Image::itt, i) = tt;
//                image(Image::iB, i) = B;
//                image(Image::itb, i) = tb;
//                image(Image::igc, i) = gc;
//                image(Image::igm, i) = gm;
            }
            image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
        }
#endif
    }
    /// get the observed flux density distrib 'image' for 2 projections for given time, freq, angle, distance, red shift
    void computeImagePW(Image & im_pj, Image & im_cj, double obs_time, double obs_freq){
        size_t cells_in_layer = EjectaID2::CellsInLayer(p_pars->ilayer);
        /// evaluateShycnhrotronSpectrum image for primary jet and counter jet
        evalImageFromPW(im_pj, obs_time, obs_freq, obsAngle, imageXXs, imageYYs);
        if (p_pars->counter_jet) // p_eats->counter_jet
            evalImageFromPW(im_cj, obs_time, obs_freq, obsAngleCJ, imageXXsCJ, imageYYsCJ);
//        auto x = im_pj.m_f_tot+im_pj.m_f_tot;
//        std::cout<<x<<"\n";
    }


#if 0
    auto evalForwardShockComovingSynchrotron(Vector & freq_arr, size_t every_it){
        return p_pars->evalForwardShockComovingSynchrotron(freq_arr,every_it);
    }
#endif

    /// evaluateShycnhrotronSpectrum light curve using Adapitve or Piece-Wise EATS method
    void evalForwardShockLightCurve(EjectaID2::STUCT_TYPE m_method_eats,
                                    Image & image, Image & im_pj, Image & im_cj,
                                    Vector & light_curve, Vector & times, Vector & freqs ){
//        Vector light_curve (times.size(), 0.0);
//        auto & m_data = p_pars->m_data;
        if ((p_pars->m_r[0] == 0.0) && (p_pars->m_r[p_pars->nr - 1] == 0.0)){
            (*p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell="<<p_pars->ishell<<", ilayer="<<p_pars->ilayer<<"]\n";
            std::fill(light_curve.begin(), light_curve.end(),0.0);
            return ; //std::move(light_curve);
        }

        double rtol = 1e-6;
//        Image image; Image im_pj; Image im_cj;
        for (size_t it = 0; it < times.size(); it++) {
//            (*p_log)(LOG_INFO,AT)<<" LC processing it="<<it<<"/"<<times.size()<<"\n";
            if (m_method_eats == EjectaID2::STUCT_TYPE::ipiecewise) {
//                    auto image = model->evalImagePW(obs_times[it], obs_freqs[it]);
                computeImagePW(im_pj, im_cj, times[it], freqs[it] );
                image.m_f_tot = (im_pj.m_f_tot + im_cj.m_f_tot);
//                evalImagePW(image, im_pj, im_cj, times[it], freqs[it]);
                light_curve[it] += image.m_f_tot;
            }
            else{
                double atol = light_curve[it] * rtol / (double)p_pars->nlayers;
                light_curve[it] += evalFluxDensA(times[it], freqs[it], atol);
            }
        }
//        return std::move(light_curve);
    }

private:

#if 0
    static double shock_synchrotron_flux_density(double Gamma, double GammaShock, double m2, double rho2, double acc_frac, double B,
                                                 double gm, double gM, double gc, double Theta, double z_cool,
                                                 double t_e, double mu, double R, double dr, double dr_tau, void * params){
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
        double nuprime = (1.0 + p_pars->z ) * p_pars->nu_obs * delta_D;
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
#endif



    /// check if during the quadrature integration there was an error
    static int check_error(void *params) {
        auto *fp = (struct Pars *) params; // removing EATS_pars for simplicity
        return fp->error;
//        return 0;
    }
    /// find angle at which currently jet ends (after lateral expansion)
    static double find_jet_edge(double phi, double theta_obs, //double cos_th_obs, double sin_th_obs,
                                double th_con_hi, Vector & a_mu, Vector & a_thj, int N,
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

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;//->getAnSynch();
//        auto * p_log = p_ params;
//        auto & m_data = p_pars->m_data;
        auto & tburst = p_pars->m_tburst;//m_data[BW::Q::itburst];
        auto & r_arr = p_pars->m_r;//m_data[BW::Q::itburst];
        auto & mu_arr = p_pars->m_mu;//m_data[BW::Q::itburst];

        if (r_arr[0] == 0.0 && r_arr[p_pars->nr - 1] == 0.0){
            (*p_pars->p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell=" << p_pars->ishell << ", ilayer=" << p_pars->ilayer << "]\n";
            return 0.0;
        }

        double a_theta = arccos(i_cos_theta);
        double mu = p_pars->obsangle(a_theta, i_phi, p_pars->theta_obs);

        p_pars->theta = a_theta;
        p_pars->o_phi = i_phi;
        p_pars->o_theta = a_theta;

        double dFnu = 0.;
        size_t ia = findIndex(mu, mu_arr, p_pars->nr);
        size_t ib = ia + 1;
        /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
        double t_e = interpSegLin(ia, ib, mu, mu_arr, tburst);
        t_e = check_emission_time(t_e, mu, p_pars->t_obs, mu_arr, (int) p_pars->nr);
        if (t_e < 0.0) {
            // REMOVING LOGGER
            (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
            return 0.;
        }

        double flux_dens = 0;
        p_pars->fluxFuncA(flux_dens);

#if 0
        /// Observed flux density evaluation (interpolate comoving spectrum)
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            Interp2d int_em(p_pars->m_freq_arr, r, p_pars->m_synch_em);
            Interp2d int_abs(p_pars->m_freq_arr, r, p_pars->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            ///----
//            auto & mu_arr = m_data[BW::Q::imu];
            size_t ia = findIndex(mu, mu_arr, p_pars->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, mu_arr, tburst);
            t_e = check_emission_time(t_e, mu, p_pars->t_obs, mu_arr, (int) p_pars->nr);
            if (t_e < 0.0) {
                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }


            p_pars->fluxFuncA();


            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_e, tburst, r_arr);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                             << " Current R grid us ["
                                             << r_arr[0] << ", "
                                             << r_arr[tburst.size() - 1] << "] "
                                             << "and tburst arr ["
                                             << tburst[0] << ", " << tburst[p_pars->nr - 1]
                                             << "] while the requried obs_time=" << p_pars->t_obs
                                             << "\n";
//                std::cerr << AT << "\n";
                return 0.;
            }
            double Gamma = interpSegLog(ia, ib, t_e, tburst, p_pars->m_gam);
            double beta = interpSegLog(ia, ib, t_e, tburst, p_pars->m_bet);
            // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor

            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z ) * p_pars->nu_obs * delta_D;
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
            double flux_dens = (intensity * r * r * dr); //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
            dFnu+=flux_dens;
            /// save the result in image
            double ctheta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::ictheta]);
            double theta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::itheta]);

            /// save current position
            p_pars->o_gam = Gamma;
            p_pars->o_r = r;
            p_pars->o_mu = mu;
            p_pars->o_flux = dFnu;
            p_pars->o_theta_j = theta;

        }
            /// Observed flux density evaluation (evaluateShycnhrotronSpectrum directly)
        else {

            size_t ia = findIndex(mu, p_pars->m_mu, p_pars->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, p_pars->m_mu, p_pars->m_tburst);
            t_e = check_emission_time(t_e, mu, p_pars->t_obs, p_pars->m_mu, (int) p_pars->nr);
            if (t_e < 0.0) {
//                auto & _r = m_data[BW::Q::iR];
//                auto & _mu = m_data[BW::Q::imu];
//                auto & _tt = m_data[BW::Q::itt];
//                auto & _tburst = m_data[BW::Q::itburst];
//                auto & _beta = m_data[BW::Q::ibeta];
//                (*p_pars->p_log)(LOG_ERR,AT) << "R      " << _r << "\n";
//                (*p_log)(LOG_ERR,AT) << "Mu     " << _mu << "\n";
//                (*p_log)(LOG_ERR,AT) << "tt     " << _tt << "\n";
//                (*p_log)(LOG_ERR,AT) << "tburst " << _tburst << "\n";
//                (*p_log)(LOG_ERR,AT) << "beta   " << _beta << "\n";


                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }

            double flux_dens = 0;
            p_pars->fluxFuncA(flux_dens);


//        double R = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->r_arr, p_pars->nr);
            double R = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iR]);
            if (!std::isfinite(R)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R is NAN in integrand for radiation" << "\n";
                // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << m_data[BW::Q::iR] << "\n";
//            std::cout << " Gamma= " << m_data[BW::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
                return 0.;
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

            if (rho < 0. || Gamma < 1. || U_p < 0. || theta <= 0. || rho2 < 0. || thick <= 0.) {
                (*p_pars->p_log)(LOG_ERR,AT) << " wrong value in interpolation to EATS surface  \n"
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

            dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2,
                                                  frac, B, gm, gM, gc, Theta, z_cool,
                                                  t_e, mu, R, thick, thick, params);
//            dFnu*=(p_pars->d_l*p_pars->d_l*2.);
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
            if (dFnu == 0 || !std::isfinite(dFnu)) {
                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " flux density is zero ( dFnu = 0 )" << "\n";
            }

            p_pars->o_gam = Gamma;
            p_pars->o_r = R;
            p_pars->o_mu = mu;
            p_pars->o_flux = dFnu;
            p_pars->o_theta_j = theta;

        }
#endif
        return dFnu;

    }
    static double costheta_integrand( double aomct, void* params ){
        auto * p_eats = (struct Pars *) params; // removing EATS_pars for simplicity
        p_eats->nevals = p_eats->nevals + 1;
        double act = 1 - aomct; // one minus cos theta 0 -> 'doppler_d' cos theta
        double integ = integrand(act, p_eats->phi, params );
        return integ;

    }
    /// integral of the (costheta_integrand)dtheta
    static double phi_integrand( double a_phi, void* params ){
        double result;
        auto * p_eats = (struct Pars *) params; // removing EATS_pars for simplicity
//        auto & m_data = p_eats->m_data;
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
                                 theta_1, p_eats->m_mu, p_eats->m_theta,//m_data[BW::Q::itheta],
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
            (*p_eats->p_log)(LOG_ERR,AT)  << " method not implemented. Exiting..." << "\n";
//            std :: cout << AT << "\n";
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
            (*p_eats->p_log)(LOG_ERR,AT) << " phi_integrand failed. Nan or negative. Result=" << result
                                         << " for t_obs=" << p_eats->t_obs
                                         << " theta_lo=" << theta_0
                                         << " theta_hi=" << theta_1
                                         << " phi=" << p_eats->phi
                                         << "\n"
                                         << " Exiting...";
//            std::cerr << AT << "\n";
            return 0.;
        }

        return result;

    }
    /// integral of the (phi_integrand)dphi
    static double integrate_theta_phi( void* params ){
        auto * p_eats = (struct Pars *) params; // removing EATS_pars for simplicity
        double atol = p_eats->atol_theta;
        // check if the parameters are set TODO wrap it into "run in safe mode"
        if(p_eats->nmax_theta < 0 || p_eats->nmax_phi < 0 || p_eats->rtol_phi < 0
           || p_eats->rtol_theta < 0 || p_eats->atol_theta < 0 || p_eats->atol_phi < 0 ){
            (*p_eats->p_log)(LOG_ERR,AT) << " one of the tolerance parameters for adaptive quad. is not set (< 0). Given:\n"
                                         << " nmax_theta=" << p_eats->nmax_theta
                                         << " nmax_phi=" << p_eats->nmax_phi
                                         << " rtol_phi=" << p_eats->rtol_phi
                                         << " rtol_theta=" << p_eats->rtol_theta
                                         << " atol_theta=" << p_eats->atol_theta
                                         << " atol_phi=" << p_eats->atol_phi << "\n"
                                         << " Exiting...\n";
//            std::cerr << AT << "\n";
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
};

#endif //SRC_EATS_H
