//
// Created by vsevolod on 21/04/23.
//

#ifndef SRC_EATS_H
#define SRC_EATS_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "blastwave/blastwave.h"
#include "blastwave/blastwave_components.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"
//#include "blastwave.h"


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
        Pars(Vector & tburst, Vector & tt, Vector & r, Vector & theta, Vector & m_gam, Vector & m_bet,
             Vector & freq_arr, Vector & synch_em, Vector & synch_abs,
             size_t & i_end_r, size_t ish, size_t il, int loglevel, void * params)
            : m_tburst(tburst), m_tt(tt),  m_r(r), m_theta(theta), m_gam(m_gam), m_bet(m_bet),
              m_freq_arr(freq_arr), m_synch_em(synch_em), m_synch_abs(synch_abs), m_i_end_r(i_end_r),
              m_params(params) {
            ishell = ish;
            ilayer= il;
            m_i_end_r = i_end_r;
            p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EATS_pars");
        }
        void * m_params;
        Vector & m_tburst; Vector & m_tt;
        Vector & m_r; Vector & m_theta; Vector & m_gam; Vector & m_bet;
        Vector & m_freq_arr; Vector & m_synch_em; Vector & m_synch_abs;
        Vector m_mu{};
//        VecVector & m_data;
//        Vector & tb_arr;
        double ctheta0 = -1.;
        size_t n_sublayers_image_adaptive = 10;
//        inline Vector & operator[](BW::Q ivn){ return this->m_data[ivn]; }
//        inline double & operator()(BW::Q ivn, size_t ir){ return this->m_data[ivn][ir]; }
//        Vector getTbGrid(size_t every_it) {
//            if ((every_it == 1)||(every_it==0)) return  m_tburst;
//            Vector tmp{};
//            for (size_t it = 0; it < nr; it = it + every_it){
//                tmp.push_back(m_tburst[it]);
//            }
//            return std::move(tmp);
//        }
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
            if (m_tburst.empty()){
                std::cerr << AT<< " isEmpty array\n";
                exit(1);
            }
//            nr = m_tburst.size();

            m_mu.resize(m_tburst.size());
//        p_eats->nr = p_pars->nr; // removing EATS_pars for simplicity
            for (size_t i = 0; i < m_i_end_r; i++) {
                m_mu[i] = (m_tburst[i] - t_obs / (1.0 + z)) / m_r[i] * CGS::c;
//                if (m_mu[i] == 0.){
//                    printf("error");
//                    exit(1);
//                }
            }
//            std::cout <<AT << ' '<< m_mu[0]<< ", "<<m_mu[1]<<" ... "<<m_mu[N-2]<<", "<<m_mu[N-1]<<"\n";

            if(m_mu[m_i_end_r - 1] < 1. ){
//                std::cout << " tburst = " << m_tburst << "\n";
//                std::cout << " r      = "<< m_r << "\n";
//                std::cout << " mu     = "<<m_mu<<"\n";
                if (m_i_end_r == m_mu.size()-1)
                    (*p_log)(LOG_WARN,AT) << " mu[-1]=" <<m_mu[m_i_end_r - 1] << " < 1 (expected >1) "
                                          << " tobs=" << t_obs
                                          << " tbutst[-1]=" << m_tburst[m_i_end_r - 1]
                                          << " R[nr-1]=" << m_r[m_i_end_r - 1]
                                          << "\n";
                int x = 1;
            }
            if(m_mu[0] > 1. ){
                (*p_log)(LOG_WARN,AT) << " mu[0]=" << m_mu[0] << " > 1 (expected <1) "
                                      << " tobs=" << t_obs
                                      << " tbutst[0]=" << m_tburst[0]
                                      << " R[0]=" << m_r[0]
                                      << "\n";
            }
            if(m_mu[0]==m_mu[m_mu.size()-1] and m_mu[0]==0){
                (*p_log)(LOG_ERR,AT) << " m_mu[0]=m_mu[-1]=0 for [il="
                    <<ilayer<<"] m_i_end_r="<<m_i_end_r<<"\n";
            }
            if (ttobs.empty())
                ttobs.resize( m_r.size(), std::numeric_limits<double>::max() );

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
        void setEatsPars(StrDbMap & pars, StrStrMap & opts, size_t n_layers, size_t ncells_, double ctheta_0,
                         double theta_c_l_, double theta_c_h_, double theta_w_, double theta_max_){
            nlayers=n_layers;
            ctheta0=ctheta_0;
            theta_c_l=theta_c_l_;
            theta_c_h=theta_c_h_;
            theta_w=theta_w_;
            theta_max=theta_max_;
            ncells=ncells_;

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

            counter_jet = getBoolOpt("counter_jet", opts, AT, p_log,true, false);

            // set options
            std::string opt;

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

            /// set synchrotron parameters
//            p_syna->setPars(pars, opts);

        }
        void updateObsPars(StrDbMap & pars){
            theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
            d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
            z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
//        check_for_unexpected_par(pars, {"theta_obs","d_l","z"});
        }

        /// -------------------------------------------------------------------------------
        size_t ishell = 0, ilayer = 0;  size_t nlayers = 0; size_t & m_i_end_r; size_t ncells;
        double theta_c_l = -1.;
        double theta_c_h = -1.;
        double theta_w = -1.;
        double theta_max=-1.;
        bool use_t_e = false;
        double z{}, d_l{}, nu_obs{}, t_obs{}, theta_obs{};
        double (* obsangle)(const double &, const double &, const double &) = nullptr;
        double (* im_xxs)( const double &, const double &, const double & ) = nullptr;
        double (* im_yys)( const double &, const double &, const double & ) = nullptr;
        Vector ttobs{}; // for PW eats only
        // -----------------------------
        void (* fluxFunc)(
                double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
                double r, double & ctheta, double theta, double phi,
                size_t ia, size_t ib, double ta, double tb, double mu, double t_obs, double nu_obs, void * params
                ) = nullptr;
//        void (* funcOptDepth)(
//                double & frac, double ctheta, double r,
//                double phi, double theta,
//                double phi_obs, double theta_obs, double r_obs,
//                double mu, double time, double freq, void * params
//                ) = nullptr;
        void (* fluxFuncA)(
                double & flux_dens, double & r, double & ctheta, double theta, double phi,
                size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs, void * params
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
        double theta=-1.;
        double o_phi=-1, o_theta=-1., o_gam=-1, o_mu=-1, o_r=-1, o_flux=-1, o_theta_j=-1; // for map
        double freq1=-1.0, freq2=-1.0;
        size_t nfreq=0;
        METHODS_QUADRATURES method_quad{};

//        METHOD_TAU method_tau{};
        bool counter_jet = true;
        int spread_method = 1; // TODO this is the only option!
        long nevals = 0; // counter for how many times the integrand was evaluated
//        size_t nr = 0;
        size_t n_tburst = 0;
//        Vector m_freq_arr{};
//        Vector m_synch_em{};
//        Vector m_synch_abs{};
        size_t i0_failed_elecctrons = 0;
        long n_fialed_electrons = 0;
        /// ---------------------------------------------------
//        double flux_dens=0., r=0., ctheta=0.,mu=0.,;
    };
    std::unique_ptr<logger> p_log;
    Pars * p_pars{};

    Vector cphis{};//= EjectaID2::getCphiGridPW( p_pars->ilayer );
//    Vector mu{};
//    std::vector<size_t> ia{};
//    std::vector<size_t> ib{};
//    Vector ta{};
//    Vector tb{};
//    Vector r{};
//    Vector fluxes{};

public:
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
    void setEatsPars(StrDbMap & pars, StrStrMap & opts, size_t nlayers, size_t ncells, double ctheta0,
                     double theta_c_l, double theta_c_h, double theta_w, double theta_max){
        p_pars->setEatsPars(pars,opts,nlayers,ncells,ctheta0,
                            theta_c_l,theta_c_h,theta_w,theta_max);
    }
    void setFluxFunc(void (* fluxFunc)( double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
                                        double r, double & ctheta, double theta, double phi,
                                        size_t ia, size_t ib, double ta, double tb, double mu, double t_obs, double nu_obs,
                                        void * params )){
        p_pars->fluxFunc = fluxFunc;
    }

    void setFluxFuncA(void (* fluxFuncA)(double & flux_dens, double & r, double & ctheta, double theta, double phi,
                                         size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs,
                                         void * params )){
        p_pars->fluxFuncA = fluxFuncA;
    }
    /// ----------------------------------------------------------------------------------------------
    /// evaluate flux density using adaptive integrator
    double evalFluxDensA(double t_obs, double nu_obs, double atol) {
        double fluxdens = 0.;
        parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                 0., M_PI, obsAngle);
        check_pars();
        double Fcoeff = CGS::cgs2mJy / (4.0 * M_PI * p_pars->d_l * p_pars->d_l); // result will be in mJy
        p_pars->atol_theta = atol/(2*Fcoeff*M_PI);// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
        p_pars->atol_phi = atol/(2*Fcoeff*M_PI);//  / (2.0 * Fcoeff);
        fluxdens += integrate_theta_phi(p_pars); // 2. because Integ_0^pi (not 2pi)
        if (p_pars->counter_jet){
            parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                     0., M_PI, obsAngleCJ);
            check_pars();
            p_pars->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
            p_pars->atol_phi = atol;//  / (2.0 * Fcoeff);
            fluxdens += integrate_theta_phi(p_pars);
        }
        return fluxdens;
    }

    /// evaluate intensity/flux density distribution using piece-wise summation
    void evalImagePW(Image & image, Image & im_pj, Image & im_cj, double obs_time, double obs_freq){
        Vector phi_grid = EjectaID2::getCphiGridPW( p_pars->ilayer );
        computeImagePW(im_pj, im_cj, obs_time, obs_freq );
        /// combine the two images (use full 'ncells' array, and fill only cells that correspond to this layer)
        for (size_t icell = 0; icell < phi_grid.size(); ++icell) {
            for (size_t ivn = 0; ivn < image.m_n_vn; ++ivn)
                image(ivn, icell) = im_pj(ivn,icell);
            for (size_t ivn = 0; ivn < image.m_n_vn; ++ivn)
                image(ivn,phi_grid.size()+icell) = im_cj(ivn,icell);
        }
        image.m_f_tot = (im_pj.m_f_tot + im_cj.m_f_tot);
//        std::cout<<image.m_f_tot<<"\n";
    }
    /// evaluate intensity/flux density distribution using piece-wise summation
    double evalImagePW_new(std::vector<VecVector> & out, double obs_time, double obs_freq, size_t offset){
        double flux_pj=0., flux_cj=0.;
        flux_pj = evalImageFromPW_new(out, obs_time, obs_freq, offset,
                                      obsAngle, imageXXs, imageYYs);
        if (p_pars->counter_jet) // p_eats->counter_jet
            flux_cj = evalImageFromPW_new(out, obs_time, obs_freq, p_pars->ncells + offset,
                                          obsAngleCJ, imageXXsCJ, imageYYsCJ);
        return flux_pj + flux_cj;
    }

    /// evaluate intensity/flux density distribution using adaptive summation
    void evalImageA(Image & image, Image & im_pj, Image & im_cj, double theta_l, double theta_h, size_t cil,
                    double obs_time, double obs_freq, double atol){
        computeImageA(im_pj, im_cj, theta_l, theta_h, cil, obs_time, obs_freq, atol);
        if ((im_cj.m_size!=im_pj.m_size)||(im_pj.m_size*2!=image.m_size)||((im_pj.m_size==0)||(im_cj.m_size==0))){
            (*p_log)(LOG_ERR,AT)<<" size mismatch. "
                << " im_pj.size()="<<im_pj.m_size
                << " im_cj.size()="<<im_cj.m_size
                << " image.size()="<<image.m_size<< "\n";
            exit(1);
        }
        for (size_t icell = 0; icell < cil; ++icell) {
            for (size_t ivn = 0; ivn < image.m_n_vn; ++ivn) {
                image(ivn, icell) = im_pj(ivn, icell);
//                image(IMG::Q::icj, icell) = 0;
            }
            for (size_t ivn = 0; ivn < image.m_n_vn; ++ivn) {
                image(ivn, cil + icell) = im_cj(ivn, icell);
//                image(IMG::Q::icj, icell) = 1;
            }
        }
        image.m_f_tot = im_pj.m_f_tot + im_cj.m_f_tot;//evalFluxDensA(obs_time, obs_freq, atol);
//        image.m_f_tot = evalFluxDensA(obs_time, obs_freq, atol);
    }
    /// evaluate intensity/flux density distribution using adaptive summation
    void evalImageA_new(std::vector<VecVector> & out, double obs_time, double obs_freq, size_t il, size_t nth, size_t nphi){
        /// evaluateShycnhrotronSpectrum image for primary jet and counter jet
        evalImageFromA_new(out, obs_time, obs_freq, il, nth, nphi, 0, obsAngle, imageXXs, imageYYs);
        if (p_pars->counter_jet) // p_eats->counter_jet
            evalImageFromA_new(out, obs_time, obs_freq, il, nth, nphi, nth*nphi, obsAngleCJ, imageXXsCJ, imageYYsCJ);
    }

    /// ----------------------------------------------------------------------------------------------



    bool evalEATSindexes(size_t & ia, size_t & ib,
                           double t_obs, double obs_angle, double ctheta_cell, double phi_cell,
                           double (*obs_angle_func)( const double &, const double &, const double & )){
        check_pars();

//        if (p_pars->theta_obs < 0.){
//            (*p_log)(LOG_ERR,AT)<<" theta_obs="<<p_pars->theta_obs<<"\n";
//            exit(1);
//        }
        if (p_pars->m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT)<<" p_pars->m_i_end_r="<<p_pars->m_i_end_r<<"\n";
            return false;
        }

//        double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
        double mu = obs_angle_func(ctheta_cell, phi_cell, obs_angle);
        for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
            p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu);
        }
        /// check if req. obs time is outside of the evolved times (throw error)
        if (t_obs < p_pars->ttobs[0]) {
            (*p_log)(LOG_WARN, AT) << "t_obs="<<t_obs<<" < ttobs_arr[0]="<<p_pars->ttobs[0]<<". Extend ttobs.\n";

//            " time grid starts too late "
//                                  << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
//                                  << " extend the grid to earlier time or request tobs at later times\n"
//                                  << " Exiting...\n";
            exit(1);
        }
        if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
            if (p_pars->m_i_end_r == p_pars->m_mu.size()-1)
                (*p_log)(LOG_WARN, AT) << "t_obs="<<t_obs<<" > ttobs_arr[i_end="<<p_pars->m_i_end_r - 1<<"]="
                    <<p_pars->ttobs[p_pars->m_i_end_r - 1]<<". Extend ttobs.\n";
//
//            << " time grid ends too early. "
//                                  << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                  << " while requested obs.time=" << t_obs
//                                  << " extend the grid to later time or request tobs at earlier times\n";
//                    std::cout << ttobs << std::endl;
//            exit(1);
            return false;
        }
//        else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
//            (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
//                                   << " from nr=" << p_pars->m_i_end_r
//                                   << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                   << " while t_obs=" << t_obs << "\n";
//            return;
//        }
        /// locate closest evolution points to the requested obs. time
        ia = findIndex(t_obs, p_pars->ttobs, p_pars->ttobs.size());
        if (ia >= p_pars->m_i_end_r - 1)
            return false; // ??
        ib = ia + 1;

        return true;
        ///
//        return r;
    }
    Vector & getTobs(){ return p_pars->ttobs; }

    void evalImageFromPW(Image & image, double t_obs, double nu_obs,
                         double (*obs_angle)( const double &, const double &, const double & ),
                         double (*im_xxs)( const double &, const double &, const double & ),
                         double (*im_yys)( const double &, const double &, const double & )){

        if ((p_pars->m_r[0] == 0.) && (p_pars->m_gam[0] == 0.)){
            (*p_log)(LOG_WARN,AT) << " [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                                  << " R[0]=0. Seems not evolved -> returning isEmpty image." << "\n";
            return;
        }
        if (p_pars->m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_i_end_r = 0\n";
            exit(1);
        }
        if (p_pars->m_tt[0] == 0 and p_pars->m_tt[p_pars->m_i_end_r] == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_tt = 0\n";
            exit(1);
        }

        parsPars(t_obs, nu_obs, 0., 0., 0., 0., obs_angle);
        check_pars();
        size_t cil = EjectaID2::CellsInLayer(p_pars->ilayer);
        image.m_f_tot = 0.;
#if 1
//        Vector ttobs( p_pars->m_r.size(), std::numeric_limits<double>::max() );
        cphis.resize(p_pars->m_i_end_r, std::numeric_limits<double>::max());//= EjectaID2::getCphiGridPW( p_pars->ilayer );
        EjectaID2::getCphiGridPW( cphis, p_pars->ilayer );
        Vector mu (p_pars->m_i_end_r, std::numeric_limits<double>::max());
        std::vector<size_t> ia(p_pars->m_i_end_r, 0);
        std::vector<size_t> ib(p_pars->m_i_end_r, 0);
        Vector ta(p_pars->m_i_end_r, 0);
        Vector tb(p_pars->m_i_end_r, 0);
        Vector r(p_pars->m_i_end_r, 0);
        Vector fluxes(p_pars->m_i_end_r, 0);
#endif
#if 0
        if (cphis.size()!=p_pars->m_i_end_r){
            cphis.resize(p_pars->m_i_end_r,std::numeric_limits<double>::max());
            mu.resize(p_pars->m_i_end_r,std::numeric_limits<double>::max());
            ia.resize(p_pars->m_i_end_r, 0);
            ib.resize(p_pars->m_i_end_r, 0);
            ta.resize(p_pars->m_i_end_r, 0);
            tb.resize(p_pars->m_i_end_r, 0);
            r.resize(p_pars->m_i_end_r, 0);
            fluxes.resize(p_pars->m_i_end_r, 0);
        }
        else{
            std::fill(cphis.begin(), cphis.end(),std::numeric_limits<double>::max());
            std::fill(mu.begin(), mu.end(),std::numeric_limits<double>::max());
            std::fill(r.begin(), r.end(),0);
            std::fill(fluxes.begin(), fluxes.end(),0.);
            std::fill(ia.begin(), ia.end(),0.);
            std::fill(ib.begin(), ib.end(),0.);
            std::fill(ta.begin(), ta.end(),0.);
            std::fill(tb.begin(), tb.end(),0.);
        }

        std::fill(cphis.begin(), cphis.end(),std::numeric_limits<double>::max());
        std::fill(mu.begin(), mu.end(),std::numeric_limits<double>::max());
        std::fill(r.begin(), r.end(),0);
        std::fill(fluxes.begin(), fluxes.end(),0.);
        std::fill(ia.begin(), ia.end(),0.);
        std::fill(ib.begin(), ib.end(),0.);
        std::fill(ta.begin(), ta.end(),0.);
        std::fill(tb.begin(), tb.end(),0.);
#endif


        /// Serial Loop... Check if times are accessible
        for (size_t i = 0; i < cil; i++) {
            double phi_cell = cphis[i];
            double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
            mu[i] = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu[i]);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (t_obs < p_pars->ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
            else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->m_i_end_r
                                       << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << t_obs << "\n";
                mu[i] = std::numeric_limits<double>::max();
                continue;
            }
//            ia[i] = findIndex(t_obs, p_pars->ttobs, p_pars->ttobs.size());
            int guess = i > 0 ? ia[i-1] : p_pars->ttobs.size()/2.;
            ia[i] = findClosestIndex(t_obs, p_pars->ttobs, guess);
            if (ia[i] >= p_pars->m_i_end_r - 1) {
                mu[i] = std::numeric_limits<double>::max();
                continue; // ??
            }
            ib[i] = ia[i] + 1;
            ta[i] = p_pars->ttobs[ia[i]];
            tb[i] = p_pars->ttobs[ib[i]];
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r[i] = interpSegLog(ia[i], ib[i], t_obs, p_pars->ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r[i] <= 0.0) || (!std::isfinite(r[i]))) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->m_i_end_r - 1] << "] "
                                      << "and tobs arr ["
                                      << p_pars->ttobs[0] << ", " << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                mu[i] = std::numeric_limits<double>::max();
                break;
            }
        }

        /// Parallel loop
        auto * _params = p_pars->m_params;
        double phi_cell=0., ctheta_cell=0., flux_dens=0., ctheta=0.;
        double tau_comp=-1., tau_bf=-1., tau_BH=-1.;
//#pragma omp parallel for private(phi_cell,ctheta_cell,flux_dens,ctheta) shared(t_obs,nu_obs,_params) num_threads( 6 )
        for (size_t i = 0; i < cil; i++) {
            if (mu[i] == std::numeric_limits<double>::max())
                continue;

            phi_cell = cphis[i];
            ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
#if 0
            double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (t_obs < p_pars->ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
//            if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
//                (*p_log)(LOG_ERR, AT) << " time grid ends too early. "
//                                      << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                      << " while requested obs.time=" << t_obs << "\n"
//                                      << " extend the grid to later time or request tobs at earlier times\n"
//                                      << " Exiting...\n";
////                    std::cout << ttobs << std::endl;
//                exit(1);
//            }
            else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->m_i_end_r
                                       << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << t_obs << "\n";
                continue;
            }
            /// locate closest evolution points to the requested obs. time
            size_t ia = findIndex(t_obs, p_pars->ttobs, p_pars->ttobs.size());
            if (ia >= p_pars->m_i_end_r - 1)
                continue; // ??
            size_t ib = ia + 1;
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_obs, p_pars->ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || (!std::isfinite(r))) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->m_i_end_r - 1] << "] "
                                      << "and tobs arr ["
                                      << p_pars->ttobs[0] << ", " << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                exit(1);
            }
#endif
            /// --------------------------------------------------------------------------
            p_pars->fluxFunc(flux_dens, tau_comp, tau_BH, tau_bf, r[i], ctheta, ctheta_cell, phi_cell,
                             ia[i], ib[i], ta[i], tb[i], mu[i], t_obs, nu_obs, _params);
            /// --------------------------------------------------------------------------
            fluxes[i] = flux_dens;
            image(IMG::Q::iintens, i) = flux_dens / (r[i] * r[i] * std::abs(mu[i])) * CGS::cgs2mJy;
            image(IMG::Q::ixr, i) = r[i] * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
            image(IMG::Q::iyr, i) = r[i] * im_yys(ctheta, phi_cell, p_pars->theta_obs);
            image(IMG::Q::ir, i) = r[i];
            image(IMG::Q::ictheta, i) = ctheta_cell;
            image(IMG::Q::icphi, i) = phi_cell;
            image(IMG::Q::imu, i) = mu[i];
//            image(IMG::Q::itau_comp, i) = tau_comp;
//            image(IMG::Q::itau_bh, i) = tau_BH;
//            image(IMG::Q::itau_bf, i) = tau_bf;
        }

        /// eval total flux
        double flux = std::accumulate(fluxes.begin(), fluxes.end(), decltype(fluxes)::value_type(0));
        image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
    }
    double evalImageFromPW_new(std::vector<VecVector> & out, double obs_time, double obs_freq, size_t offset,
                         double (*obs_angle)( const double &, const double &, const double & ),
                         double (*im_xxs)( const double &, const double &, const double & ),
                         double (*im_yys)( const double &, const double &, const double & )){

        /// out is [i_vn][ish][itheta_iphi]
        bool save_im = true;
        if (out.empty())
            save_im = false;

        if ((p_pars->m_r[0] == 0.) && (p_pars->m_gam[0] == 0.)){
            (*p_log)(LOG_WARN,AT) << " [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                                  << " R[0]=0. Seems not evolved -> returning isEmpty image." << "\n";
            return 0;
        }
        if (p_pars->m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_i_end_r = 0\n";
            exit(1);
        }
        if (p_pars->m_tt[0] == 0 and p_pars->m_tt[p_pars->m_i_end_r] == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_tt = 0\n";
            exit(1);
        }

        parsPars(obs_time, obs_freq, 0., 0., 0., 0., obs_angle);
        check_pars();
        size_t cil = EjectaID2::CellsInLayer(p_pars->ilayer);


        /// allocate memory for quantities that help
        Vector mu (p_pars->m_i_end_r, std::numeric_limits<double>::max());
        std::vector<size_t> ia(p_pars->m_i_end_r, 0);
        std::vector<size_t> ib(p_pars->m_i_end_r, 0);
        Vector ta(p_pars->m_i_end_r, 0);
        Vector tb(p_pars->m_i_end_r, 0);
        Vector r(p_pars->m_i_end_r, 0);
        Vector fluxes(p_pars->m_i_end_r, 0);
        /// Find the region for EATS interpoaltion
        for (size_t i = 0; i < cil; i++) {
            double cphi = (double)i * 2.0 * M_PI / (double)cil;
            double ctheta_cell = p_pars->ctheta0;
            mu[i] = obs_angle(ctheta_cell, cphi, p_pars->theta_obs);
            /// fill-in the time on EATS plane
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu[i]);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (obs_time < p_pars->ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << obs_time << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
            else if ((obs_time > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->m_i_end_r
                                       << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << obs_time << "\n";
                mu[i] = std::numeric_limits<double>::max();
                continue;
            }
            ///
            int guess = i > 0 ? (int)ia[i-1] : (int)(p_pars->ttobs.size()/2);
            ia[i] = findClosestIndex(obs_time, p_pars->ttobs, guess);
            if (ia[i] >= p_pars->m_i_end_r - 1) {
                mu[i] = std::numeric_limits<double>::max();
                continue; // ??
            }
            ib[i] = ia[i] + 1;
            ta[i] = p_pars->ttobs[ia[i]];
            tb[i] = p_pars->ttobs[ib[i]];
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r[i] = interpSegLog(ia[i], ib[i], obs_time, p_pars->ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r[i] <= 0.0) || (!std::isfinite(r[i]))) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->m_i_end_r - 1] << "] "
                                      << "and tobs arr ["
                                      << p_pars->ttobs[0] << ", " << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                mu[i] = std::numeric_limits<double>::max();
                break;
            }
        }

        auto * _params = p_pars->m_params;
        double phi_cell=0., ctheta_cell=0., flux_dens=0., ctheta=0.;
        double tau_comp=-1., tau_bf=-1., tau_BH=-1.;
        double tot_flux = 0.;
        /// Perform EATS integration
        for (size_t i = 0; i < cil; i++) {
            if (mu[i] == std::numeric_limits<double>::max())
                continue;

            phi_cell = (double) i * 2.0 * M_PI / (double) cil;
            ctheta_cell = p_pars->ctheta0;

            p_pars->fluxFunc(flux_dens, tau_comp, tau_BH, tau_bf, r[i], ctheta, ctheta_cell, phi_cell,
                             ia[i], ib[i], ta[i], tb[i], mu[i], obs_time, obs_freq, _params);

            if (!std::isfinite(flux_dens) || flux_dens < 0){
                (*p_log)(LOG_ERR,AT) << " flux_dens="<<flux_dens<<"\n";
                exit(1);
            }
            /// save data in container
            if (save_im) {
                out[IMG::Q::iintens][p_pars->ishell][offset + i] = flux_dens / (r[i] * r[i] * std::abs(mu[i])) * CGS::cgs2mJy;
                out[IMG::Q::ixr][p_pars->ishell][offset + i] = r[i] * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
                out[IMG::Q::iyr][p_pars->ishell][offset + i] = r[i] * im_yys(ctheta, phi_cell, p_pars->theta_obs);
                out[IMG::Q::ir][p_pars->ishell][offset + i] = r[i];
                out[IMG::Q::ictheta][p_pars->ishell][offset + i] = ctheta_cell;
                out[IMG::Q::icphi][p_pars->ishell][offset + i] = phi_cell;
                out[IMG::Q::imu][p_pars->ishell][offset + i] = mu[i];
            }
            tot_flux += flux_dens;
        }
        return (tot_flux * CGS::cgs2mJy);
    }
#if 1

    void evalImageFromA(Image & image, double t_obs, double nu_obs, double atol,
                        double theta_c_l, double theta_c_h, size_t cil,
                        double (*obs_angle)( const double &, const double &, const double & ),
                        double (*im_xxs)( const double &, const double &, const double & ),
                        double (*im_yys)( const double &, const double &, const double & )) {
        if (image.m_size==0){
            (*p_log)(LOG_ERR,AT) << " image is empty\n";
            exit(1);
        }
        image.m_f_tot = 0.;
        if (cil<1){
            (*p_log)(LOG_ERR,AT) << " cil < 1\n";
            exit(1);
        }
        parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                 0., 2. * M_PI, obs_angle);
        check_pars();
        double summed_intensity = 0;
        double Fcoeff = cgs2mJy / (4. * M_PI * p_pars->d_l * p_pars->d_l);
        size_t ii = 0;
        for (ii = 0; ii < cil; ii++) {

            double cphi = image.gerArr(IMG::Q::icphi)[ii];
            double ctheta = image.gerArr(IMG::Q::ictheta)[ii];

            if (ctheta < p_pars->theta_c_l)
                continue;
            // if jet is spreading, compute the upper boundary of the jet
            double th_b = find_jet_edge(cphi, p_pars->theta_obs, //p_pars->cos_theta_obs, p_pars->sin_theta_obs,
                                        p_pars->theta_c_h, p_pars->m_mu, p_pars->m_theta,//m_data[BW::Q::itheta],
                                        (int)p_pars->m_i_end_r, p_pars->obsangle);
            if (th_b >= CGS::pi/2.)
                continue;
            double th_a = ( p_pars->theta_c_l / p_pars->theta_c_h ) * th_b; // ???
            if (ctheta < th_a || ctheta > th_b)
                continue;
            // compute intensity
            double r = 0., mu = 0., gam = 0., ctheta_bw = 0.;
            double intensity = integrand(cos(ctheta ), cphi,
                                         r, mu, gam, ctheta_bw, p_pars);

//            if (ii>0){
//                double dtheta = image.gerArr(IMG::Q::ictheta)[ii]-image.gerArr(IMG::Q::ictheta)[ii-1];
////                double dphi = image.gerArr(IMG::Q::icphi)[ii]-image.gerArr(IMG::Q::icphi)[ii-1];
//                intensity*=dtheta;
//            }
            if ((r == 0)||(r!=r)){
                (*p_log)(LOG_ERR,AT) << " r=0"<<"\n";
                exit(1);
            }

            summed_intensity += intensity;

            double x = r * im_xxs(ctheta, cphi, p_pars->theta_obs);
            double y = r * im_yys(ctheta, cphi, p_pars->theta_obs);

            image(IMG::Q::iintens,ii) = intensity * Fcoeff / ( r * r * std::abs(mu) );//* phi_edges.size();
            image(IMG::Q::ir,ii) = r;
            image(IMG::Q::imu,ii) = mu;
            image(IMG::Q::ixr,ii) = x;//-1 * sin(obs_theta) * ( sin(_cthetas[ii]) * sin(_cphis[ii])) + cos(obs_theta)*cos(_cthetas[ii]);
            image(IMG::Q::iyr,ii) = y;//
            image(IMG::Q::ictheta,ii) = ctheta;//
            image(IMG::Q::icphi,ii) = cphi;//
        }

        summed_intensity *= Fcoeff;
        image.m_f_tot = summed_intensity;
    }
    void evalImageFromA_new(std::vector<VecVector> & out, double obs_time, double obs_freq, size_t il, size_t nth, size_t nphi, size_t ii_ofset,
                        double (*obs_angle)( const double &, const double &, const double & ),
                        double (*im_xxs)( const double &, const double &, const double & ),
                        double (*im_yys)( const double &, const double &, const double & )) {
        /// settings for intensity search
        double phi0 = 0.;
        double phi1 = 2. * M_PI;
        double theta0 = p_pars->theta_c_l;
        double theta1 = M_PI / 2.;

        if (out.size()==0){
            (*p_log)(LOG_ERR,AT) << " image is empty\n";
            exit(1);
        }
        if (il<0){
            (*p_log)(LOG_ERR,AT) << " il < 0\n";
            exit(1);
        }
        parsPars(obs_time, obs_freq, p_pars->theta_c_l, p_pars->theta_c_h,
                 0., 2. * M_PI, obs_angle);
        check_pars();

        size_t ii = 0;
        double summed_intensity = 0;
        double Fcoeff = cgs2mJy / (4. * M_PI * p_pars->d_l * p_pars->d_l);

        for (size_t iphi = 0; iphi < nphi; iphi++){
            double th_b_min = std::numeric_limits<double>::max();
            double th_b_max = std::numeric_limits<double>::min();
            double cphi = phi0 + (double)iphi * (phi1 - phi0) / (double)nphi;
            /// check what the extend of theta
            // if jet is spreading, compute the upper boundary of the jet
            double th_b = find_jet_edge(cphi, p_pars->theta_obs, //p_pars->cos_theta_obs, p_pars->sin_theta_obs,
                                        p_pars->theta_c_h, p_pars->m_mu, p_pars->m_theta,//m_data[BW::Q::itheta],
                                        (int)p_pars->m_i_end_r, p_pars->obsangle);
//            /// check the theta extend for this phi
//            if (th_b > th_b_max)
//                th_b_max = th_b;
//            if (th_b < th_b_min)
//                th_b_min = th_b;
//            /// TODO check if we need to check this
//            if (th_b_min < theta0)
//                th_b_min = theta0;
//
//            if ((th_b_min <= th_b_max) || (th_b_max == 0.)){
//                (*p_log)(LOG_ERR,AT)<<" no theta area to consider\n";
//                exit(1);
//            }

            double th_a = ( p_pars->theta_c_l / p_pars->theta_c_h ) * th_b; // ???
//            std::cout << AT<< " iphi="<<iphi<<" phi="<<cphi<<" th_a="<<th_a<<" th_b="<<th_b<<"\n";

            for (size_t ith = 0; ith < nth; ith++){
                double cth = th_a + (double)ith * (th_b - th_a) / (double)nth;
                // compute intensity
                double r = 0., mu = 0., gam = 0., ctheta_bw = 0.;
                double intensity = integrand(cos(cth ), cphi,
                                             r, mu, gam, ctheta_bw, p_pars);
                if ((r == 0)||(r!=r)){
                    (*p_log)(LOG_ERR,AT) << " r=0"<<"\n";
                    exit(1);
                }
                summed_intensity += intensity;
                double x = r * im_xxs(cth, cphi, p_pars->theta_obs);
                double y = r * im_yys(cth, cphi, p_pars->theta_obs);
                // ---
                out[IMG::Q::ictheta][il][ii_ofset+ii] = cth;
                out[IMG::Q::icphi][il][ii_ofset+ii] = cphi;
                out[IMG::Q::iintens][il][ii_ofset+ii] = intensity;
                out[IMG::Q::ir][il][ii_ofset+ii] = r;
                out[IMG::Q::ixr][il][ii_ofset+ii] = x;
                out[IMG::Q::iyr][il][ii_ofset+ii] = y;
                out[IMG::Q::imu][il][ii_ofset+ii] = mu;
                ii ++;
            }
        }
        summed_intensity *= Fcoeff;
    }
#endif

#if 0
    void evalImageFromA(Image & image, double t_obs, double nu_obs, double atol,
                        double theta_c_l, double theta_c_h, size_t isublayer,
                        double (*obs_angle)( const double &, const double &, const double & ),
                        double (*im_xxs)( const double &, const double &, const double & ),
                        double (*im_yys)( const double &, const double &, const double & )) {
        if ((p_pars->m_r[0] == 0.0) && (p_pars->m_r[p_pars->m_i_end_r - 1] == 0.0)) {
            (*p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell=" << p_pars->ishell << ", ilayer=" << p_pars->ilayer
                    << "]\n";
            return; //std::move(light_curve);
        }
        parsPars(t_obs, nu_obs, p_pars->theta_c_l, p_pars->theta_c_h,
                 0, 2.*CGS::pi, obs_angle);
        check_pars();

        /// update for integration
        p_pars->current_theta_cone_hi = theta_c_h;
        p_pars->current_theta_cone_low = theta_c_l;
        double dtheta = theta_c_h - theta_c_l;

        /// total cells across all sublayers
//        size_t ncells = 0;
//        for (size_t i = 0; i < _theta_c_l.size(); i++)
//            ncells += _nphis[i];
//        Vector fluxes(ncells, 0);
//        if (image.m_size != ncells) {
//            (*p_log)(LOG_ERR, AT) << " image.m_size=" << image.m_size << " while ncells=" << ncells << "\n";
//            exit(1);
//        }
        ///
        size_t cil = EjectaID2::CellsInLayer(isublayer);
        Vector phi_edges = EjectaID2::getCphiGridPW_(cil); // +1 to get boundaries
        Vector fluxes(phi_edges.size(),0.);
        size_t ii = 0;
//        double phi_low = 0.;
//        double phi_high = 2. * M_PI;
//        double dphi = (phi_high - phi_low) / (double) nphi;
        double Fcoeff = CGS::cgs2mJy / (4.0 * M_PI * p_pars->d_l * p_pars->d_l); // result will be in mJy

        for (size_t iphi = 1; iphi <= phi_edges.size(); iphi++) {
            double i_phi_0 = phi_edges[iphi-1];
            double i_phi_1 = phi_edges[iphi];
            double cphi = 0.5 * (i_phi_1 - i_phi_0);
            double dphi = i_phi_1-i_phi_0;
            p_pars->current_phi_hi = i_phi_1; // Update for integrator
            p_pars->current_phi_low = i_phi_0;
            // -------------------------------------
            p_pars->atol_theta = atol;// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
            p_pars->atol_phi = atol;//  / (2.0 * Fcoeff);
            fluxes[ii] = integrate_theta_phi(p_pars); // 2. because Integ_0^pi (not 2pi)
            /// to get the radius to the (ctheta0,cphi)
            double ctheta = theta_c_l + 0.5 * (theta_c_h - theta_c_l);


            double r = 0., mu = 0., gam = 0., ctheta_bw = 0.;
            double _ = integrand(std::cos(ctheta),
                                 i_phi_0 + cphi, r, mu, gam, ctheta_bw, p_pars);
//            fluxdens += 2. * Fcoeff * integrate_theta_phi(p_pars);
//            fluxes[ii] *= (dtheta * dphi);

            double x = r * im_xxs(ctheta, i_phi_0+cphi, p_pars->theta_obs);
            double y = r * im_yys(ctheta, i_phi_0+cphi, p_pars->theta_obs);
            image(IMG::Q::iintens, ii) = fluxes[ii];// * Fcoeff / (r * r * std::abs(mu));
            image(IMG::Q::ixr, ii) = x;
            image(IMG::Q::iyr, ii) = y;
            image(IMG::Q::ir, ii) = r;
            image(IMG::Q::ictheta, ii) = ctheta;
            image(IMG::Q::icphi, ii) = i_phi_0 + cphi;
            image(IMG::Q::imu, ii) = mu;
            if (std::abs(mu) < 1e-5){
                int x = 1;
            }
            ii++;
        }
        double flux = std::accumulate(fluxes.begin(), fluxes.end(), decltype(fluxes)::value_type(0));
        image.m_f_tot = flux;// * Fcoeff; /// flux in mJy
    }
#endif
#if 0
    void evalImageFromPW_old(Image & image, double t_obs, double nu_obs,
                         double (*obs_angle)( const double &, const double &, const double & ),
                         double (*im_xxs)( const double &, const double &, const double & ),
                         double (*im_yys)( const double &, const double &, const double & )){

        if ((p_pars->m_r[0] == 0.) && (p_pars->m_gam[0] == 0.)){
            (*p_log)(LOG_WARN,AT) << " [ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
                                  << " R[0]=0. Seems not evolved -> returning empty image." << "\n";
            return;
        }
        if (p_pars->m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_i_end_r = 0\n";
            exit(1);
        }
        if (p_pars->m_tt[0] == 0 and p_pars->m_tt[p_pars->m_i_end_r] == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_tt = 0\n";
            exit(1);
        }

        parsPars(t_obs, nu_obs, 0., 0., 0., 0., obs_angle);
        check_pars();

//        Vector ttobs( p_pars->m_r.size(), std::numeric_limits<double>::max() );
        Vector cphis = EjectaID2::getCphiGridPW( p_pars->ilayer );
        Vector mu (p_pars->m_i_end_r, std::numeric_limits<double>::max());
        std::vector<size_t> ia(p_pars->m_i_end_r, 0);
        std::vector<size_t> ib(p_pars->m_i_end_r, 0);
        Vector ta(p_pars->m_i_end_r, 0);
        Vector tb(p_pars->m_i_end_r, 0);
        Vector r(p_pars->m_i_end_r, 0);

        /// Serial Loop... Check if times are accessible
        for (size_t i = 0; i < cphis.size(); i++) {
            double phi_cell = cphis[i];
            double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
            mu[i] = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu[i]);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (t_obs < p_pars->ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
            else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->m_i_end_r
                                       << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << t_obs << "\n";
                mu[i] = std::numeric_limits<double>::max();
                continue;
            }
            ia[i] = findIndex(t_obs, p_pars->ttobs, p_pars->ttobs.size());
            if (ia[i] >= p_pars->m_i_end_r - 1) {
                mu[i] = std::numeric_limits<double>::max();
                continue; // ??
            }
            ib[i] = ia[i] + 1;
            ta[i] = p_pars->ttobs[i];
            tb[i] = p_pars->ttobs[i];
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r[i] = interpSegLog(ia[i], ib[i], t_obs, p_pars->ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r[i] <= 0.0) || (!std::isfinite(r[i]))) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->m_i_end_r - 1] << "] "
                                      << "and tobs arr ["
                                      << p_pars->ttobs[0] << ", " << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                mu[i] = std::numeric_limits<double>::max();
                break;
            }
        }

        /// Parallel loop
        double flux = 0.;
//#pragma omp parallel for num_threads( 6 )
        for (size_t i = 0; i < cphis.size(); i++) {

            double phi_cell = cphis[i];
            double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
#if 1
            double mu = obs_angle(ctheta_cell, phi_cell, p_pars->theta_obs);
            for (size_t i_ = 0; i_ < p_pars->m_i_end_r; i_++) {
                p_pars->ttobs[i_] = p_pars->m_tt[i_] + p_pars->m_r[i_] / CGS::c * (1.0 - mu);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (t_obs < p_pars->ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
//            if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
//                (*p_log)(LOG_ERR, AT) << " time grid ends too early. "
//                                      << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                      << " while requested obs.time=" << t_obs << "\n"
//                                      << " extend the grid to later time or request tobs at earlier times\n"
//                                      << " Exiting...\n";
////                    std::cout << ttobs << std::endl;
//                exit(1);
//            }
            else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
                (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
                                       << " from nr=" << p_pars->m_i_end_r
                                       << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                       << " while t_obs=" << t_obs << "\n";
                continue;
            }
            /// locate closest evolution points to the requested obs. time
            size_t ia = findIndex(t_obs, p_pars->ttobs, p_pars->ttobs.size());
            if (ia >= p_pars->m_i_end_r - 1)
                continue; // ??
            size_t ib = ia + 1;
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_obs, p_pars->ttobs, p_pars->m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || (!std::isfinite(r))) {
                (*p_log)(LOG_ERR, AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                      << " Current R grid us ["
                                      << p_pars->m_r[0] << ", "
                                      << p_pars->m_r[p_pars->m_i_end_r - 1] << "] "
                                      << "and tobs arr ["
                                      << p_pars->ttobs[0] << ", " << p_pars->ttobs[p_pars->m_i_end_r - 1]
                                      << "] while the requried obs_time=" << p_pars->t_obs
                                      << "\n";
                exit(1);
            }
#endif
            /// ----------------------------------------------------
//            if (image.m_n_vn==IMG::m_names.size()) {
            double flux_dens;
            double ctheta;
            p_pars->fluxFunc(flux_dens, r, ctheta, ctheta_cell, phi_cell,
                             ia, ib, mu, t_obs, nu_obs, p_pars->ttobs, p_pars->m_params);
            /// ----------------------------------------------------
            flux += flux_dens;
            image(IMG::Q::iintens, i) =
                    flux_dens / (r * r * std::abs(mu)) *
                    CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
            image(IMG::Q::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
            image(IMG::Q::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
            image(IMG::Q::imu, i) = mu;
//            }
//            if (image.m_n_vn==IMG_TAU::m_names.size()) {
//                double ctheta;
//                double tau_Compton,tau_BH,tau_bf;
//                p_pars->funcOptDepth(double & frac, double ctheta, double r,
//                        double phi, double theta,
//                        double phi_obs, double theta_obs, double r_obs,
//                        double mu, double time, double freq, p_pars->m_params);
//                /// ----------------------------------------------------
//                image(IMG_TAU::Q::itau_comp, i) = tau_Compton;
//                image(IMG_TAU::Q::itau_bf, i) = tau_bf;
//                image(IMG_TAU::Q::itau_bh, i) = tau_BH;
//                image(IMG_TAU::Q::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
//                image(IMG_TAU::Q::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
//                image(IMG_TAU::Q::imu, i) = mu;
//            }

        }
        image.m_f_tot = flux * CGS::cgs2mJy; /// flux in mJy
    }
#endif
    /// get the observed flux density distrib 'image' for 2 projections for given time, freq, angle, distance, red shift
    void computeImagePW(Image & im_pj, Image & im_cj, double obs_time, double obs_freq){
        /// evaluateShycnhrotronSpectrum image for primary jet and counter jet
        evalImageFromPW(im_pj, obs_time, obs_freq, obsAngle, imageXXs, imageYYs);
        if (p_pars->counter_jet) // p_eats->counter_jet
            evalImageFromPW(im_cj, obs_time, obs_freq, obsAngleCJ, imageXXsCJ, imageYYsCJ);
    }

    /// get the observed flux density distrib 'image' for 2 projections for given time, freq, angle, distance, red shift
    void computeImageA(Image & im_pj, Image & im_cj, double theta_l , double theta_h, size_t cil,
                       double obs_time, double obs_freq, double atol){
        /// evaluateShycnhrotronSpectrum image for primary jet and counter jet
        evalImageFromA(im_pj, obs_time, obs_freq, atol, theta_l, theta_h, cil,
                       obsAngle, imageXXs, imageYYs);
        if (p_pars->counter_jet) // p_eats->counter_jet
            evalImageFromA(im_cj, obs_time, obs_freq, atol, theta_l, theta_h, cil,
                           obsAngleCJ, imageXXsCJ, imageYYsCJ);
    }


    /// evaluateShycnhrotronSpectrum light curve using Adapitve or Piece-Wise EATS method
    void evalLC(EjectaID2::STUCT_TYPE m_method_eats,
                Image & image, Image & im_pj, Image & im_cj,
                Vector & light_curve, Vector & times, Vector & freqs ){
//        Vector light_curve (times.size(), 0.0);
//        auto & m_data = p_pars->m_data;
        if ((p_pars->m_r[0] == 0.0) && (p_pars->m_r[p_pars->m_i_end_r - 1] == 0.0)){
            (*p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell="<<p_pars->ishell<<", ilayer="<<p_pars->ilayer<<"]\n";
            std::fill(light_curve.begin(), light_curve.end(),0.0);
            return ; //std::move(light_curve);
        }

        double rtol = 1e-10;
//        Image image; Image im_pj; Image im_cj;
//#pragma omp parallel for shared(m_method_eats,image,im_pj,im_cj,light_curve,times,freqs) num_threads( 6 )
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
    void evalLC_new(Vector & out, EjectaID2::STUCT_TYPE m_method_eats, Vector & times, Vector & freqs ){
        double rtol = 1e-10;
        std::vector<VecVector> empty{};
        for (size_t it = 0; it < times.size(); it++) {
            if (m_method_eats == EjectaID2::STUCT_TYPE::ipiecewise)
                out[it] = evalImagePW_new(empty,times[it],freqs[it],0);

            else{
                double atol = out[it] * rtol / (double)p_pars->nlayers;
                out[it] += evalFluxDensA(times[it], freqs[it], atol);
            }
        }
    }
private:

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
//        std::cout <<AT << ' '<< a_mu[0]<< ", "<<a_mu[1]<<" ... "<<a_mu[N-2]<<", "<<a_mu[N-1]<<"\n";

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
        while(theta_b - theta_a > 1.0e-5 && i < 100) {
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
    static double integrand( double i_cos_theta, double i_phi, double & r, double & mu, double & gam, double & ctheta,
                             void* params ){

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
//        auto & p_syna = p_pars->p_syna;//->getAnSynch();
//        auto * p_log = p_ params;
//        auto & m_data = p_pars->m_data;
        auto & tburst = p_pars->m_tburst;//m_data[BW::Q::itburst];
        auto & r_arr = p_pars->m_r;//m_data[BW::Q::itburst];
        auto & mu_arr = p_pars->m_mu;//m_data[BW::Q::itburst];

        if (r_arr[0] == 0.0 && r_arr[p_pars->m_i_end_r - 1] == 0.0){
            (*p_pars->p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell=" << p_pars->ishell << ", ilayer=" << p_pars->ilayer << "]\n";
            return 0.0;
        }

        double a_theta = arccos(i_cos_theta);
        mu = p_pars->obsangle(a_theta, i_phi, p_pars->theta_obs);

        p_pars->theta = a_theta;
        p_pars->o_phi = i_phi;
        p_pars->o_theta = a_theta;
        p_pars->o_mu = mu;

//        double dFnu = 0.;
        size_t ia = findIndex(mu, mu_arr, p_pars->m_i_end_r);
        size_t ib = ia + 1;
        /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
        double t_e = interpSegLin(ia, ib, mu, mu_arr, tburst);
        t_e = check_emission_time(t_e, mu, p_pars->t_obs, mu_arr, (int) p_pars->m_i_end_r);
        if (t_e < 0.0||!std::isfinite(t_e)) {
            // REMOVING LOGGER
            (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
            return 0.;
        }
        /// -----------------------------------------
        double flux_dens = 0; //r = 0, ctheta=0.;
        p_pars->fluxFuncA(flux_dens, r, ctheta, a_theta, i_phi,
                          ia, ib, mu, t_e, p_pars->t_obs, p_pars->nu_obs, p_pars->m_params);
//        std::cout<<"fluxdens="<<flux_dens<<"\n";
        /// ----------------------------------------
        gam = interpSegLog(ia, ib, t_e, tburst, p_pars->m_gam);

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
#if 1
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



        }
#endif
        p_pars->o_gam = gam;
        p_pars->o_r = r;
        p_pars->o_flux = flux_dens;
        p_pars->o_theta_j = ctheta;
        return flux_dens;
    }
    static double costheta_integrand( double aomct, void* params ){
        auto * p_eats = (struct Pars *) params; // removing EATS_pars for simplicity
        p_eats->nevals = p_eats->nevals + 1;
        double act = 1 - aomct; // one minus cos theta 0 -> 'doppler_d' cos theta
        double r = 0., mu=0., gam=0., ctheta=0.;
        double integ = integrand(act, p_eats->phi, r, mu, gam, ctheta, params );
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
                                 p_eats->m_i_end_r,
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
        double sin_half_theta_0 = std::sin(0.5 * theta_0);
        double sin_half_theta_1 = std::sin(0.5 * theta_1);
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
//        printf("   a_phi: %.6lf (%.6le)\n", a_phi, result);
        return result;

    }
    /// integral of the (phi_integrand)dphi
    static double integrate_theta_phi( void* params ){
        auto * p_eats = (struct Pars *) params; // removing EATS_pars for simplicity
        double atol = p_eats->atol_theta;
        // check if the parameters are set
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
        double Fcoeff = CGS::cgs2mJy / (4*M_PI * p_eats->d_l* p_eats->d_l);

        double result;
        double phi_a;
        switch (p_eats->method_quad) {

            case INT_TRAP_FIXED:
                result = 2 * Fcoeff * trap(&phi_integrand, phi_0, phi_1,
                              p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_TRAP_ADAPT:
                result = 2 * Fcoeff * trap_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_SIMP_FIXED:
                result = 2 * Fcoeff * simp(&phi_integrand, phi_0, phi_1,
                              p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_SIMP_ADAPT:
                result = 2 * Fcoeff * simp_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_ROMB_ADAPT:
                phi_a = phi_0 + 0.5*(phi_1-phi_0);
                result = 2 * Fcoeff * romb(&phi_integrand, phi_0, phi_a,
                              p_eats->nmax_phi, atol,
                              p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                              check_error, nullptr, nullptr);
                result += 2 * Fcoeff * romb(&phi_integrand, phi_a, phi_1,
                               p_eats->nmax_phi,
                               (atol + p_eats->rtol_phi * result),
                               p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                               check_error, nullptr, nullptr);
                break;
            case INT_TRAP_NL:
                result = 2 * Fcoeff * trapNL_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_HYBRID:
                result = 2 * Fcoeff * hybrid_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      0, check_error, nullptr, nullptr);
                break;
            case INT_CADRE:
                result = 2 * Fcoeff * cadre_adapt(&phi_integrand, phi_0, phi_1,
                                     p_eats->nmax_phi, atol,
                                     p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                     0, check_error, nullptr, nullptr);
                break;
            case INT_GK49_ADAPT:
                result = 2 * Fcoeff * gk49_adapt(&phi_integrand, phi_0, phi_1,
                                    p_eats->nmax_phi, atol,
                                    p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                    0, check_error);
                break;
            case INT_GK715_ADAPT:
                result = 2 * Fcoeff * gk715_adapt(&phi_integrand, phi_0, phi_1,
                                     p_eats->nmax_phi, atol,
                                     p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                     0, check_error);
                break;
            case INT_GK1021_ADAPT:
                result = 2 * Fcoeff * gk1021_adapt(&phi_integrand, phi_0, phi_1,
                                      p_eats->nmax_phi, atol,
                                      p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                      0, check_error);
                break;
        }
        return result;
    }
};

#endif //SRC_EATS_H
