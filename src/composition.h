//
// Created by vsevolod on 29/03/23.
//

#ifndef SRC_COMPOSITION_H
#define SRC_COMPOSITION_H

#include "utilitites/utils.h"
#include "utilitites/pch.h"
#include "utilitites/logger.h"
#include "utilitites/interpolators.h"
#include "utilitites/H5Easy.h"

struct LippunerHeating{
    std::unique_ptr<logger> p_log;
    explicit LippunerHeating(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LippunerHeating");
    }
    std::string fpath{};
    double eps0{};
    /// ------------------------------
    Vector ye{};// = np.array(dfile["Ye"])
    Vector s{};// = np.array(dfile["s"])
    Vector tau{};// = np.array(dfile["tau"])
    /// -------------------------------
    Vector alpha{};// = np.array(dfile["alpha"])
    Vector A{};// = np.array(dfile["A"])
    Vector B1{};// = np.array(dfile["B1"])
    Vector beta1{};// = np.array(dfile["beta1"])
    Vector B2{};// = np.array(dfile["B2"])
    Vector beta2{};// = np.array(dfile["beta2"])
    Vector B3{};// = np.array(dfile["B3"])
    Vector beta3{};// = np.array(dfile["beta3"])
    /// -------------------------------
    bool is_loaded = false;
    /// -------------------------------
    void setPars(StrDbMap pars, StrStrMap opts){
        std::string workdir = getStrOpt("workingdir",opts,AT,p_log,"",true);
        std::string fname = getStrOpt("fname_nuc",opts,AT,p_log,"",true);
        eps0 = getDoublePar("eps_nuc0",pars,AT,p_log,0,true);
        fpath = workdir + fname;
    }
    /// Load Lippuner+15 table in h5 form
    void load(){
        if (is_loaded)
            return;
        if (!std::experimental::filesystem::exists(fpath))
            throw std::runtime_error("File not found. " + fpath);
        LoadH5 ldata;
        ldata.setFileName(fpath);

        ldata.setVarName("Ye");     ye = ldata.getDataVDouble();
        ldata.setVarName("s");      s = ldata.getDataVDouble();
        ldata.setVarName("tau");    tau = ldata.getDataVDouble();
        ldata.setVarName("alpha");  alpha = ldata.getDataVDouble();
        ldata.setVarName("A");      A = ldata.getDataVDouble();
        ldata.setVarName("B1");     B1 = ldata.getDataVDouble();
        ldata.setVarName("beta1");  beta1 = ldata.getDataVDouble();
        ldata.setVarName("B2");     B2 = ldata.getDataVDouble();
        ldata.setVarName("beta2");  beta2 = ldata.getDataVDouble();
        ldata.setVarName("B3");     B3 = ldata.getDataVDouble();
        ldata.setVarName("beta3");  beta3 = ldata.getDataVDouble();
        is_loaded = true;
        (*p_log)(LOG_INFO, AT) << "Lippuner heating rate table is loaded\n";
    }
    /// Interpolate Lippuner+15 table using triliniar interpolation
    double evalHeatingRate(double ye_val, double s_val, double tau_val, double time){
        if (!is_loaded) load();
        auto FA_ = TriliniarInterpolation(ye, s, tau, A);
        double A_ = FA_.interp(ye_val, s_val, tau_val);
        /// -----------------------------------------------------------
        auto Falpha_ = TriliniarInterpolation(ye, s, tau, alpha);
        double alpha_ = Falpha_.interp(ye_val, s_val, tau_val);
        /// -----------------------------------------------------------
        auto FB1_ = TriliniarInterpolation(ye, s, tau, B1);
        double B1_ = FB1_.interp(ye_val, s_val, tau_val);
        auto Fbeta1_ = TriliniarInterpolation(ye, s, tau, beta1);
        double beta1_ = Fbeta1_.interp(ye_val, s_val, tau_val);
        /// -----------------------------------------------------------
        auto FB2_ = TriliniarInterpolation(ye, s, tau, B2);
        double B2_ = FB2_.interp(ye_val, s_val, tau_val);
        auto Fbeta2_ = TriliniarInterpolation(ye, s, tau, beta2);
        double beta2_ = Fbeta2_.interp(ye_val, s_val, tau_val);
        /// -----------------------------------------------------------
        auto FB3_ = TriliniarInterpolation(ye, s, tau, B3);
        double B3_ = FB3_.interp(ye_val, s_val, tau_val);
        auto Fbeta3_ = TriliniarInterpolation(ye, s, tau, beta3);
        double beta3_ = Fbeta3_.interp(ye_val, s_val, tau_val);
        /// -----------------------------------------------------------
        // Eq.~4 in https://iopscience.iop.org/article/10.1088/0004-637X/815/2/82/pdf
        double Q = A_ * pow(time, -alpha_) \
                      + B1_ * std::exp(-time / beta1_) \
                      + B2_ * std::exp(-time / beta2_) \
                      + B3_ * std::exp(-time / beta3_);
        double res = (eps0 / 2.e18)*Q; /// convention
        return res;
    }
};

struct BarnsThermalization{

    std::unique_ptr<logger> p_log;
    /// From Barns et al (Values for mej and vej)
    Array a; Array b; Array d;
    Array x = {std::log10(1.e-3),std::log10(5e-3),std::log10(1e-2),std::log10(5e-2)};
    Array y = {0.1,0.2,0.3};
    /// From Barns et al for mej/vej^2
    Array x_barnes = {0.011, 0.025, 0.0556, 0.1, 0.111, 0.125, 0.25, 0.5, 0.5556, 1., 1.25, 5.};
    Array a_barnes = {8.16, 4.52, 3.20, 2.01, 2.19, 1.90, 1.31, 0.81, 0.95, 0.56, 0.55, 0.27};
    Array b_barnes = {1.19, 0.62, 0.45, 0.28, 0.31, 0.28, 0.21, 0.19, 0.15, 0.17, 0.13, 0.10};
    Array d_barnes = {1.52, 1.39, 1.39, 1.12, 1.32, 1.21, 1.13, 0.86, 1.13, 0.74, 0.90, 0.60};
    ///
//    METHOD_THERMALIZATION methodThermalization{};
    double thef{};

    explicit BarnsThermalization(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BarnsThermalization");

        VecVector a_ = {{2.01, 0.81, 0.56, 0.27},
                        {4.52, 1.90, 1.31, 0.55},
                        {8.16, 3.20, 2.19, 0.95}};
        VecVector b_ = {{0.28, 0.19, 0.17, 0.10},
                        {0.62, 0.28, 0.21, 0.13},
                        {1.19, 0.45, 0.31, 0.15}};
        VecVector d_ = {{1.12, 0.86, 0.74, 0.60},
                        {1.39, 1.21, 1.13, 0.90},
                        {1.52, 1.39, 1.32, 1.13}};

        a.resize(a_.size() * a_[0].size());
        size_t k = 0;
        for (size_t i = 0; i < a_.size(); ++i){
            for (size_t j = 0; j < a_.size(); ++j){
                a[k] = a_[i][j]; k++;
            }
        }
        k = 0;
        b.resize(b_.size() * b_[0].size());
        for (size_t i = 0; i < b_.size(); ++i){
            for (size_t j = 0; j < b_.size(); ++j){
                b[k] = b_[i][j]; k++;
            }
        }
        k = 0;
        d.resize(d_.size() * d_[0].size());
        for (size_t i = 0; i < d_.size(); ++i){
            for (size_t j = 0; j < d_.size(); ++j){
                d[k] = d_[i][j]; k++;
            }
        }
    }

    void setPars(StrDbMap pars, StrStrMap opts){ }

    /// Thermalization efficiency from Barns et al
    double evalThermalization(double m_iso, double beta, double time){
        Interp2d fa(x, y, a);
        Interp2d fb(x, y, b);
        Interp2d fd(x, y, d);

        double xnew = std::log10(m_iso); // mass [Msun] 4*CGS::pi/ncells
        double ynew = beta; //velocity [c]


        size_t ia = findIndex(xnew, x, x.size());
        if (xnew < x[0])
            xnew = x[0];
        if (xnew > x[x.size()-1])
            xnew = x[x.size()-1];
        if (ynew < y[0])
            ynew = y[0];
        if (ynew > y[y.size()-1])
            ynew = y[y.size()-1];
        if (ia >= x.size() - 1){
            std::cerr << " failed 2d interpolation\n"; exit(1);
        }
        size_t ib = ia + 1;
        size_t ia_ = findIndex(xnew, y, y.size());
        if (ia_ >= y.size() - 1){
            std::cerr << " failed 2d interpolation\n"; exit(1);
        }
        size_t ib_ = ia_ + 1;
        double coeff0 = fa.InterpolateBilinear(xnew, ynew, ia_, ib_, ia, ib);
        double coeff1 = fb.InterpolateBilinear(xnew, ynew, ia_, ib_, ia, ib);
        double coeff2 = fd.InterpolateBilinear(xnew, ynew, ia_, ib_, ia, ib);

//        double coeff0 = fa.InterpolateBilinear(xnew, ynew);
//        double coeff1 = fb.InterpolateBilinear(xnew, ynew);
//        double coeff2 = fd.InterpolateBilinear(xnew, ynew);


//        double coeff0 = fa.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);
//        double coeff1 = fb.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);
//        double coeff2 = fd.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);

        double time_days = time / CGS::day;
        double tmp = 2. * coeff1 * std::pow(time_days, coeff2);
        double val = 0.36 * (std::exp(-coeff0 * time_days) + std::log( 1. + tmp) / tmp );
        if (!std::isfinite(val)){
            (*p_log)(LOG_ERR,AT)<<" nan in thermalization\n";
            exit(1);
        }
        if (val < 0 || val > 1){
            (*p_log)(LOG_ERR, AT) << " check the thermalization efficiency = "<<val << "\n";
            exit(1);
        }
        return val;
    }
    double evalThermalization1D(double mass, double beta, double time){
        Interp1d fa_1d (x_barnes, a_barnes);
        Interp1d fb_1d (x_barnes, b_barnes);
        Interp1d fd_1d (x_barnes, d_barnes);

        double mass_ = mass / CGS::solar_m;
        double vel_ = beta;
        double xnew = mass_ / vel_ / vel_; // x=m/v^2

        double coeff0 = fa_1d.Interpolate(xnew, Interp1d::METHODS::iLagrangeLinear);
        double coeff1 = fb_1d.Interpolate(xnew, Interp1d::METHODS::iLagrangeLinear);
        double coeff2 = fd_1d.Interpolate(xnew, Interp1d::METHODS::iLagrangeLinear);

        double time_days = time * CGS::day;
        double tmp = 2.*coeff1 * std::pow(time_days, coeff2);
        double val = 0.36 * (std::exp(-coeff0*time_days) + std::log(1. + tmp)/tmp );
        if (val < 0 || val > 1){
            (*p_log)(LOG_ERR, AT) << " check the thermalization efficiency = "<<val << "\n";
            exit(1);
        }
        return val;
    }

};

double smoothclamp(double x, double x1, double x2, double y1, double y2) {
    /// y1 + (y2-y1)*(lambda t: np.where(t < 0 , 0, np.where( t <= 1 , 3.*t**2-2.*t**3, 1 ) ) )( np.log10(x/x1)/np.log10(x2/x1) )
    double t = std::log10(x / x1) / std::log10(x2 / x1);
    double value = 0;
    if (t < 0) {
        value = 0;
    } else if (t <= 1) {
        value = 3 * std::pow(t, 2.0) - 2 * std::pow(t, 3.0);
    } else {
        value = 1;
    }
    return y1 + (y2 - y1) * value;
}

/// Model suggested by Korobkin + 2012
struct KorobkinHeating{

    double alpha{}; double t0eps{}; double eps0{}; double sigma0;
    double cnst_a_eps_nuc{}; double cnst_b_eps_nuc{}; double cnst_t_eps_nuc{};
    double kappa{};
    std::unique_ptr<logger> p_log;
    KorobkinHeating(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "KorobkinHeating");
    }
    void setPars(StrDbMap pars, StrStrMap opts){
        /// para
        alpha = getDoublePar("alpha", pars, AT, p_log, 0.0, true);
        eps0 = getDoublePar("eps_nuc0", pars, AT, p_log, 0.0, true);
        t0eps = getDoublePar("t0eps", pars, AT, p_log, 0.0, true);
        sigma0 = getDoublePar("sigma", pars, AT, p_log, 0.0, true);
        cnst_a_eps_nuc = getDoublePar("a_eps_nuc", pars, AT, p_log, 0.0, true);
        cnst_b_eps_nuc = getDoublePar("b_eps_nuc", pars, AT, p_log, 0.0, true);
        cnst_t_eps_nuc = getDoublePar("t_eps_nuc", pars, AT, p_log, 0.0, true);
    }

    /// Nuclear Heating due to r-process
    double yevar_eps_nuc(double a_eps_nuc, double b_eps_nuc, double t_eps_nuc, double time){
        double time_ = time * CGS::day;
        double tmp = std::min(std::max(4. * time_ - 4., -20.), 20.);   /// t_eps_nuc seem,s missing!
        double val = a_eps_nuc + b_eps_nuc/(1. + std::exp(tmp));
        return val;
    }

    double calc_eps_nuc(double kappa, double time, double eps0,
                        double a_eps_nuc, double b_eps_nuc, double t_eps_nuc) {
        double weight = smoothclamp(kappa, 1., 10., 1., 0.);
        double val = eps0 * ((1. - weight) + weight * yevar_eps_nuc(a_eps_nuc, b_eps_nuc, t_eps_nuc, time));
        return val;
    }

    double evalHeatingRate(double mass, double beta, double time, double kappa){
        /// Perego et al 2017 ApJL
        double eps_nuc = calc_eps_nuc(kappa, time, eps0, cnst_a_eps_nuc, cnst_b_eps_nuc,
                                      cnst_t_eps_nuc);
        double oneoverpi = 0.31830988618;
        double res = eps_nuc * std::pow(0.5 - oneoverpi * std::atan((time-t0eps)/sigma0), alpha);
        return res;
    }
};

struct TanakaOpacity{
    Array Ye{0.01,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50};
    Array kappa{30.1,30.0,29.9,22.30,5.60,5.36,3.30,0.96,0.1};
    TanakaOpacity(){}
    void setPars(StrDbMap & pars, StrStrMap & opts){}
    double evaluateOpacity(double ye){
        if (ye < Ye[0])
            ye = Ye[0];
        if (ye > Ye[Ye.size()-1])
            ye = Ye[Ye.size()-1];
        auto interp = Interp1d(Ye, kappa);
        double kap = interp.Interpolate(ye, Interp1d::METHODS::iLinear);
        if (kap < 0 or !std::isfinite(kap)){
            std::cerr << AT <<" kap = "<<kap<<"\n";
            exit(1);
        }
        return kap;
    }
};

class NuclearAtomic{
    enum METHOD_HEATING { iPBR, iLR, iNone };
    enum METHOD_THERMALIZATION { iBKWM, iBKWM_1d, iConstTherm, iNoneTh };
    METHOD_HEATING methodHeating{};
    METHOD_THERMALIZATION methodThermalization{};
    struct Pars{ double eps_th0=0.; double kappa=0.; double eps_nuc_thermalized=0.; double eps_th=0.; };
    bool do_nucinj= false;
    /// -----------------------------------------
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    std::unique_ptr<LippunerHeating> p_lip;
    std::unique_ptr<BarnsThermalization> p_barns;
    std::unique_ptr<KorobkinHeating> p_korob;
    std::unique_ptr<TanakaOpacity> p_tanaka;
    int m_loglevel{};
public:
    NuclearAtomic(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel,"NuclearAtomic");
        p_pars = std::make_unique<Pars>();
        m_loglevel=loglevel;
    }

    void setPars(StrDbMap pars, StrStrMap opts){
        do_nucinj = getBoolOpt("do_nucinj", opts, AT, p_log, false, true);
        if (!do_nucinj)
            return;

        std::string opt = "method_heating";
        METHOD_HEATING m_methodHeating;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_methodHeating = METHOD_HEATING::iNone;
        }
        else{
            if(opts.at(opt) == "none")
                m_methodHeating = METHOD_HEATING::iNone;
            else if(opts.at(opt) == "PBR")
                m_methodHeating = METHOD_HEATING::iPBR;
            else if(opts.at(opt) == "LR")
                m_methodHeating = METHOD_HEATING::iLR;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << " Possible options: "
                                      << " none " << " PBR " << " LR "
                                      << " \n";
                exit(1);
            }
        }
        methodHeating = m_methodHeating;

        opt = "method_thermalization";
        METHOD_THERMALIZATION m_methodThermalization;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_methodThermalization = METHOD_THERMALIZATION::iConstTherm;
        }
        else{
            if(opts.at(opt) == "none")
                m_methodThermalization = METHOD_THERMALIZATION::iNoneTh;
            else if(opts.at(opt) == "const")
                m_methodThermalization = METHOD_THERMALIZATION::iConstTherm;
            else if(opts.at(opt) == "BKWM")
                m_methodThermalization = METHOD_THERMALIZATION::iBKWM;
            else if(opts.at(opt) == "BKWM_1d")
                m_methodThermalization = METHOD_THERMALIZATION::iBKWM_1d;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << " Possible options: "
                                      << " const " << " BKWM " << " BKWM_1d "
                                      << " \n";
                exit(1);
            }
        }
        methodThermalization = m_methodThermalization;


//        p_pars->eps0 = getDoublePar("eps_nuc0", pars, AT, p_log, 0.0, true);
        /// set method for thermalization calculations
        switch (methodThermalization) {
            case iBKWM:
                p_barns = std::make_unique<BarnsThermalization>(m_loglevel);
                p_barns->setPars(pars, opts);
                break;
            case iBKWM_1d:
                p_barns = std::make_unique<BarnsThermalization>(m_loglevel);
                p_barns->setPars(pars, opts);
                break;
            case iConstTherm:
                p_pars->eps_th0 = getDoublePar("eps_thermalization", pars, AT, p_log, 0.0, true);
                break;
            case iNoneTh:
                p_pars->eps_th0 = 0.;
                break;
        }
        /// set method for heating rate calculations
        switch (methodHeating) {
            case iPBR:
                p_korob = std::make_unique<KorobkinHeating>(m_loglevel);
                p_korob->setPars(pars, opts);
                break;
            case iLR:
                p_lip = std::make_unique<LippunerHeating>(m_loglevel);
                p_lip->setPars(pars, opts);
                break;
            case iNone:
                break;
        }
        /// set opacity calculator
        p_tanaka = std::make_unique<TanakaOpacity>();
    }

    void update(const double ye, const double m_iso, const double beta,
                const double time, const double radius, const double s){
        /// update opacity
        if (!do_nucinj)
            return;
        double kappa = p_tanaka->evaluateOpacity(ye);
        if (kappa < 0 or !std::isfinite(kappa)){
            (*p_log)(LOG_ERR,AT)<<" kappa = "<<kappa<<"\n";
            exit(1);
        }
        p_pars->kappa = kappa;

        /// update thermalization efficiency
        double eps_th = 0.;
        switch (methodThermalization) {
            case iBKWM:
                eps_th = p_barns->evalThermalization(m_iso/CGS::solar_m, beta, time);
                break;
            case iBKWM_1d:
                eps_th = p_barns->evalThermalization1D(m_iso/CGS::solar_m, beta, time);
                break;
            case iConstTherm:
                eps_th = p_pars->eps_th0;
                break;
        }
        /// update the nuclear heating rate
        double eps_nuc = 0.; double tau = 0.;
        switch (methodHeating) {
            case iPBR:
                eps_nuc = p_korob->evalHeatingRate(m_iso, beta, time, kappa);
                break;
            case iLR:
                tau = 0.5 * 2.718 * (radius / beta / CGS::c); // expansion timescale
                eps_nuc = p_lip->evalHeatingRate(ye, s, tau, time);
                break;
            case iNone:
                break;
        }
        /// evaluate the actual heating reate
        double eps_nuc_thermalized = eps_nuc * eps_th;
        if (!std::isfinite(eps_nuc_thermalized)){
            (*p_log)(LOG_ERR,AT)<<" nan in eps_nuc_thermalized\n";
            exit(1);
        }
        /// ---------------------------------------
//        std::cout << "eps_nuc="<<eps_nuc<<"\n";
        p_pars->eps_nuc_thermalized = eps_nuc_thermalized;
        p_pars->eps_th = eps_th;
    }

    std::unique_ptr<Pars> & getPars(){return p_pars; };

};

#endif //SRC_COMPOSITION_H
