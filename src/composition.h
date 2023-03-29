//
// Created by vsevolod on 29/03/23.
//

#ifndef SRC_COMPOSITION_H
#define SRC_COMPOSITION_H

#include "utils.h"
#include "pch.h"
#include "logger.h"
#include "interpolators.h"

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

class Composition{
    enum METHOD_HEATING { iPBR, iLR, iNone };
    enum METHOD_THERMALIZATION { iBKWM, iBKWM_1d, iConstTherm };
    struct Pars{
        double tcomov0{};
        METHOD_HEATING methodHeating{};
        METHOD_THERMALIZATION methodThermalization{};
        double thef{};
        double alpha{}; double t0{}; double eps0{}; double sigma0{};
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
        Pars(){
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

    };
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
public:
    Composition(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel,"Composition");
        p_pars = std::make_unique<Pars>();
    }

    void setPars(StrDbMap pars, StrStrMap opts){
        p_pars->methodHeating = METHOD_HEATING::iPBR;
        p_pars->methodThermalization = METHOD_THERMALIZATION::iBKWM;
        p_pars->thef = 0.1;
        p_pars->alpha = 0.;
        p_pars->t0 = 0.;
        p_pars->sigma0 = 0.;
        p_pars->eps0 = 0.;
    }
    void setEjectaPars(const double tcomov0){
        p_pars->tcomov0 = tcomov0;
    }

    /// Thermalization efficiency from Barns et al
    double evalThermalization(double mass, double beta, double time){
        if (p_pars->methodThermalization == METHOD_THERMALIZATION::iBKWM){
            Interp2d fa(p_pars->x, p_pars->y, p_pars->a);
            Interp2d fb(p_pars->x, p_pars->y, p_pars->b);
            Interp2d fd(p_pars->x, p_pars->y, p_pars->d);

            double xnew=std::log10(mass) / CGS::solar_m; // mass [Msun] 4*CGS::pi/ncells
            double ynew=beta; //velocity [c]
            double coeff0 = fa.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);
            double coeff1 = fb.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);
            double coeff2 = fd.Interpolate(xnew, ynew, Interp2d::METHODS::iLagrangeLinear);

            double time_days = time * CGS::day;
            double tmp = 2.*coeff1 * std::pow(time_days, coeff2);
            double val = 0.36 * (std::exp(-coeff0 * time_days) + std::log( 1. + tmp) / tmp );
            if (val < 0 || val > 1){
                (*p_log)(LOG_ERR, AT) << " check the thermalization efficiency = "<<val << "\n";
                exit(1);
            }
            return val;
        }
        else if (p_pars->methodThermalization == METHOD_THERMALIZATION::iBKWM_1d){
            Interp1d fa_1d (p_pars->x_barnes, p_pars->a_barnes);
            Interp1d fb_1d (p_pars->x_barnes, p_pars->b_barnes);
            Interp1d fd_1d (p_pars->x_barnes, p_pars->d_barnes);

            double mass_ = mass / CGS::solar_m;
            double vel_ = beta;
            double xnew = mass / vel_;

            double coeff0 = fa_1d.Interpolate(xnew, Interp2d::METHODS::iLagrangeLinear);
            double coeff1 = fb_1d.Interpolate(xnew, Interp2d::METHODS::iLagrangeLinear);
            double coeff2 = fd_1d.Interpolate(xnew, Interp2d::METHODS::iLagrangeLinear);

            double time_days = time * CGS::day;
            double tmp = 2.*coeff1 * std::pow(time_days, coeff2);
            double val = 0.36 * (std::exp(-coeff0*time_days) + std::log(1. + tmp)/tmp );
            if (val < 0 || val > 1){
                (*p_log)(LOG_ERR, AT) << " check the thermalization efficiency = "<<val << "\n";
                exit(1);
            }
            return val;
        }
        else{
            return p_pars->thef;
        }
    }

    /// Nuclear Heating due to r-process
    double yevar_eps_nuc(double a_eps_nuc, double b_eps_nuc, double t_eps_nuc, double time){
        double time_ = time * CGS::day;
        double tmp = std::min(std::max(4.*time_ - 4., -20.), 20.);   /// t_eps_nuc seem,s missing!
        double val = a_eps_nuc + b_eps_nuc/(1.+std::exp(tmp));
    }

    double calc_eps_nuc(double kappa, double time, double eps0,
                        double a_eps_nuc, double b_eps_nuc, double t_eps_nuc) {
        double weight = smoothclamp(kappa, 1., 10., 1., 0.);
        return eps0 * ((1. - weight) + weight * yevar_eps_nuc(a_eps_nuc, b_eps_nuc, t_eps_nuc, time))
    }

    double evalNucHeat(double mass, double beta, double time){

        double eps_th = evalThermalization(mass, beta, time);

        if (p_pars->methodHeating == METHOD_HEATING::iPBR){
            /// Perego et al 2017 ApJL
        }
        else if (p_pars->methodHeating == METHOD_HEATING::iLR){
            /// Lippuner & Roberts 2016 ApJ
        }
        else if (p_pars->methodHeating == METHOD_HEATING::iNone){
            /// None for no Ye dependence
            double eps_nuc = p_pars->eps0;
            double mass_ = mass / CGS::solar_m;
            double vel_ = beta;
            double oneoverpi = 0.31830988618;
            ///  Korobkin heating rates
            double val = eps_nuc
                       * std::pow(0.5 - oneoverpi * std::atan((time-p_pars->t0)/p_pars->sigma0), p_pars->alpha)
                       * (2.*eps_th);

        }
        else{
            (*p_log)(LOG_ERR, AT) << " wrong option for heating rate calc." << "\n";
            exit(1);
        }


//        double ejecta_co_t0 = 1.3; // eq. 15 Sun (2017)
//        double ejecta_co_sigma = 0.11;
//        double val = std::pow(0.5 - 1./CGS::pi*std::atan((tcomov-ejecta_co_t0)/ejecta_co_sigma), 1.3);
//        val *= 4e49*mej/1.e-2/CGS::solar_m;
//        return val;
    }

    void getHeatingRate()

};

#endif //SRC_COMPOSITION_H
