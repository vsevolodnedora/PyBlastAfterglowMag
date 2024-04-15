//
// Created by vsevolod on 1/15/24.
//

#ifndef SRC_KERNELS_H
#define SRC_KERNELS_H

# define M_SQRT3    1.7320508075688772	/* sqrt(3) */

#include "../utilitites/utils.h"
#include "numeric_model.h"

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_spline.h>

class SSCKernel{

    std::vector<VecVector> kernel; // [i_nu_ssc, i_gam, i_nu_syn]
    VecVector kernel_integral;     // [i_gam, i_nu_syn]
//    State & ssc; State & ele; State & syn;
    double (*m_func) (double,double,double) = nullptr;

public:

//    SSCKernel(State & ele, State & syn, State & ssc, double (*func)(double,double,double))
//            :ele(ele),syn(syn),ssc(ssc),m_func(func){}

    SSCKernel() = default;

    void allocate(State & ele, State & syn, State & ssc, double (*m_func_) (double,double,double)){
        //        std::vector<VecVector> kernel{};
        kernel.resize(ssc.numbins);
        /// allocate memory
        for (auto & vecvec : kernel) {
            vecvec.resize(ele.numbins);
            for (auto & vec : vecvec)
                vec.resize(syn.numbins, 0.);
        }
        m_func = m_func_;
    }

    void evalSSCkernel(State & ele, State & syn, State & ssc){
        /// compute kernel
//    auto * eq = ssc_kernel_nava;
        for (size_t i = 0; i < ssc.numbins; i++){
            for (size_t j = 0; j < ele.numbins; j++){
                for (size_t k = 0; k < syn.numbins; k++){
                    kernel[i][j][k] = m_func(ssc.e[i], ele.e[j], syn.e[k]);//[i_nu_ssc, i_gam, i_nu_syn]
                }
            }
        }
//        return std::move( kernel );
    }

    /*
     * Pre-compute the inner integral in SSC cooling over the grid of SSC frequencies, as
     * :math: \int_{\nu_0}^{\nu_1} \nu d\nu F(q,p)
     * See Eq.15 in Huang22
     */
    void evalSSCgridIntergral(State & ele, State & syn, State & ssc){

//        double h = 6.6261e-27; //  erg*sec

        /// allocate memory
        kernel_integral.resize(ele.numbins);
        for (auto & vec : kernel_integral)
            vec.resize(syn.numbins, 0.);

        /// compute integral over SSC. frequencies
        for (size_t j = 0; j < ele.numbins; j++){
            for (size_t k = 0; k < syn.numbins; k++){
                double res = 0;
                for (size_t i = 0; i < ssc.numbins-1; i++) {
                    res += kernel[i][j][k] * ssc.e[i+1] * ssc.de[i+1];
//                           * (ssc.e[i+1] / 8.093440820813486e-21) // freq
//                           * (ssc.de[i+1] / 8.093440820813486e-21) // dfreq
//                           * CGS::h;
                }
                kernel_integral[j][k] = res; // [i_gamma, i_energy_syn]
            }
        }
//        int x = 1;
    }

    static double sscNava(double nu, double gam, double nu_tilde){

        double en = nu * CGS::h_by_mec2; // Eq.34
        double en_tilde = nu_tilde * CGS::h_by_mec2;

        double K = 0.; // Eq.33
        if (en_tilde / (4. * gam*gam) < en and en <= en_tilde)
            K = (en / en_tilde) - (1. / (4. * gam*gam));
        else if (en_tilde < en && en < (4. * gam*gam * en_tilde / (1. + 4. * gam * en_tilde))){
            double q = en / (4. * gam * en_tilde * (gam - en));
            K = 2. * q * std::log(q) \
              + (1. + 2. * q) * (1. - q) \
              + .5 * (1. - q) * std::pow(4. * gam * en_tilde * q, 2.) / (1. + 4. * gam * en_tilde * q);
        }
        else
            K = 0;

        if (not std::isfinite(K))
            throw std::runtime_error("nan in ssc");

        return K;
    }
    /**
     * See Huang+20, equations in the paragraph after eq.15
     * @param en_ssc
     * @param gam
     * @param en_syn
     * @return
     */
    static double sscHUANG(double freq_ssc, double gam, double freq_syn){
//        double freq_ssc = en_ssc / 8.093440820813486e-21;
//        double freq_syn = en_syn / 8.093440820813486e-21;
//        double mec2 = 8.187e-7; // ergs
//        double h = 6.6261e-27; // erg*sec
        double g = gam * h * freq_syn / mec2;
        double w = h * freq_ssc / gam / mec2;
        double q = w / (4.*g*(1.-w));
        double freq_ssc_max = gam * CGS::mec2 * 4. * g / (4. * g + 1.) / CGS::h;
        double f = 0.;
        if ((freq_ssc > freq_syn) && (freq_ssc < freq_ssc_max))
            f = 2. * q * std::log(q) \
            + (1. + 2. * q) * (1. - q) \
            + (1. / 2.) * (1. - q) * std::pow(4. * q * g, 2.) / (1. + 4. * g * q);
        return f;
    }
    /**
     * @return  [i_nu_ssc, i_gam, i_nu_syn]
     */
    const auto & getKernel(){ return kernel; }

    const auto & getKernelInteg(){ return kernel_integral; }
};

/**
 * evaluate the cheb poly for the given value of x
 * @param coeff
 * @param order
 * @param a
 * @param b
 * @param x
 * @return
 */
double cheb_eval(const double * coeff, int order, double a, double b, double x){
    double d = 0.0;
    double dd = 0.0;
    double y = (2.0 * x - a - b) / (b - a);
    double y2 = 2.0 * y;

    double temp = 0.;
    for (int j = order; j>0; j--){
        temp = d;
        d = y2 * d - dd + coeff[j];
        dd = temp;
    }
    temp = d;
    d = y * d - dd + 0.5 * coeff[0];
    return d;
}

class SynKernel{
    VecVector kernel{}; // [i_freq, i_gam]
    double (*m_func) (double,double,double) = nullptr;
//    State & ele; State & syn;
public:
    SynKernel() = default;

    void allocate(size_t ngam, size_t nfreq){
        /// allocate memory for the kernel (kernel is not static, but size is)
        kernel.resize(nfreq);
        for (auto & arr : kernel)
            arr.resize(ngam, 0.);
    }
    void setFunc(double (*m_func_) (double,double,double)){ m_func = m_func_; }

    static inline double nu_crit(double B, double gam){
        double nu_crit = (3.*CGS::qe*B*gam*gam) / (4.*M_PI*CGS::me*CGS::c);
        return nu_crit;
    }

    static inline double nu_larmor(double B, double gam){
        double nu_larmor = (3.*CGS::qe*B)/(2.*pi*CGS::me*CGS::c);
        return nu_larmor;
    }

    /**
     * this implements the synchrotron kernal from GSL without relying on GSL
     * @param B
     * @param en
     * @param gam
     * @return
     */
    static double synchGSL(double B, double nu, double gam){

        double synchrotron1_data[13] = {
                30.364682982501076273,
                17.079395277408394574,
                4.560132133545072889,
                0.549281246730419979,
                0.372976075069301172e-01,
                0.161362430201041242e-02,
                0.481916772120371e-04,
                0.10512425288938e-05,
                0.174638504670e-07,
                0.22815486544e-09,
                0.240443082e-11,
                0.2086588e-13,
                0.15167e-15,
        };
        int synchrotron1_cs_order = 12;
        double synchrotron1_cs_a = -1.0;
        double synchrotron1_cs_b = 1.0;

        double synchrotron1a_data[23] = {
                2.1329305161355000985,
                0.741352864954200240e-01,
                0.86968099909964198e-02,
                0.11703826248775692e-02,
                0.1645105798619192e-03,
                0.240201021420640e-04,
                0.35827756389389e-05,
                0.5447747626984e-06,
                0.838802856196e-07,
                0.13069882684e-07,
                0.2053099071e-08,
                0.325187537e-09,
                0.517914041e-10,
                0.83002988e-11,
                0.13352728e-11,
                0.2159150e-12,
                0.349967e-13,
                0.56994e-14,
                0.9291e-15,
                0.152e-15,
                0.249e-16,
                0.41e-17,
                0.7e-18,
        };

        int synchrotron1a_cs_order = 22;
        double synchrotron1a_cs_a = -1.0;
        double synchrotron1a_cs_b = 1.0;

        double synchrotron2_data[12] = {
                0.4490721623532660844,
                0.898353677994187218e-01,
                0.81044573772151290e-02,
                0.4261716991089162e-03,
                0.147609631270746e-04,
                0.3628633615300e-06,
                0.66634807498e-08,
                0.949077166e-10,
                0.1079125e-11,
                0.10022e-13,
                0.77e-16,
                0.5e-18,
        };

        int synchrotron2_cs_order = 11;
        double synchrotron2_cs_a = -1.0;
        double synchrotron2_cs_b = 1.0;


        double ec = 1.5 * B / 4.14e13; // Bcritical
        double en = nu * CGS::h_by_mec2;
        double x = en / (ec * gam * gam);

//        double M_SQRT2 = np.sqrt(2.0);
//        double M_SQRT3 = std::sqrt(3.0);
//        double M_PI = np.pi;
//        double GSL_SQRT_DBL_EPSILON = 1.4901161193847656e-08;
//        double GSL_LOG_DBL_MIN = -7.0839641853226408e02;

        double c0 = M_PI / M_SQRT3;

        double c01 = 0.2257913526447274323630976;
        double cond1 = 2 * M_SQRT2 * 1.4901161193847656e-08;//GSL_SQRT_DBL_EPSILON;
        double cond3 = -8.0 * (-7.0839641853226408e02) / 7.0 ; // GSL_LOG_DBL_MIN

        double result;
        if (x < cond1){
            double z = std::pow(x, 1.0 / 3.0);
            double cf = 1 - 8.43812762813205e-01 * z * z;
            result = 2.14952824153447863671 * z * cf;
        }
        else if (x <= 4.0){
            double px = std::pow(x, 1.0 / 3.0);
            double px11 = std::pow(px, 11.);
            double t = x * x / 8.0 - 1.0;
            double result_c1 = cheb_eval(
                    synchrotron1_data,
                    synchrotron1_cs_order,
                    synchrotron1_cs_a,
                    synchrotron1_cs_b,
                    t
            );

            double result_c2 = cheb_eval(
                    synchrotron2_data,
                    synchrotron2_cs_order,
                    synchrotron2_cs_a,
                    synchrotron2_cs_b,
                    t
            );

            result = px * result_c1 - px11 * result_c2 - c0 * x;
        }
        else if (x < cond3){
            double t = (12.0 - x) / (x + 4.0);
            double result_c1 = cheb_eval(
                    synchrotron1a_data,
                    synchrotron1a_cs_order,
                    synchrotron1a_cs_a,
                    synchrotron1a_cs_b,
                    t
            );
            result = std::sqrt(x) * result_c1 * std::exp(c01 - x);
        }
        else
            result = 0.0;

        result *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / CGS::mec2; // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)
        return result;
    };
#if 0
    static double cycSyn(double b, double nu, double gamma){
        const double arg[47] = {
                0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.03, 0.05, 0.07, 0.1,
                0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5,
                8., 8.5, 9., 9.5, 10., 12., 14., 16., 18., 20., 25., 30., 40., 50.};
        const double var[47]	= {-1.002e+0, -9.031e-1, -7.696e-1, -6.716e-1, -5.686e-1, -4.461e-1, -3.516e-1,
                             -2.125e-1, -1.537e-1, -1.192e-1, -8.725e-2, -4.383e-2, -3.716e-2, -4.528e-2, -5.948e-2, -7.988e-2,
                             -1.035e-1, -1.296e-1, -1.586e-1, -1.838e-1, -3.507e-1, -5.214e-1, -6.990e-1, -8.861e-1, -1.073e+0,
                             -1.267e+0, -1.470e+0, -1.670e+0, -1.870e+0, -2.073e+0, -2.279e+0, -2.483e+0, -2.686e+0, -2.893e+0,
                             -3.097e+0, -3.303e+0, -3.510e+0, -3.717e+0, -4.550e+0, -5.388e+0, -6.230e+0, -7.075e+0, -7.921e+0,
                             -1.005e+1, -1.218e+1, -1.646e+1, -2.076e+1};
        double nu_c, x, emisfunc, nu_larmor, psquared, ngamma;
//        double nu = en * CGS::ergToFreq;

        auto acc_syn =  gsl_interp_accel_alloc();
        auto syn = gsl_spline_alloc(gsl_interp_cspline,47);
        gsl_spline_init(syn,arg,var,47);

        //this is in the synchrotron regime
        double charg=4.8e-10;
//        double cee =GSL_CONST_CGSM_SPEED_OF_LIGHT;
//        double emgm =GSL_CONST_CGSM_MASS_ELECTRON;
//        double emerg=(GSL_CONST_CGSM_MASS_ELECTRON*pow(GSL_CONST_CGSM_SPEED_OF_LIGHT,2.));
        if (gamma > 2.) {
            emisfunc = synchGSL(b, nu, gamma);
//            nu_c = (3.*CGS::qe*b*pow(gamma,2.))/(4.*pi*CGS::me*CGS::c);
//            x = nu/nu_c;
//            if(x <= 1.e-4){
//                emisfunc = 4.*pi*pow(x/2.,(1./3.))/(sqrt(3.)*2.68);
//            }
//            else if(x > 50.){
//                emisfunc = sqrt(pi*x/2.)*exp(-x);
//            }
//            else{
//                emisfunc = pow(10.,gsl_spline_eval(syn, x, acc_syn));
//            }
//            emisfunc *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * b / CGS::mec2;
        } else { //cyclotron regime
            nu_larmor = (CGS::qe*b)/(2.*pi*CGS::me*CGS::c);
            x = nu/nu_larmor;
            psquared = pow(gamma,2.)-1.;
            emisfunc = (2.*psquared)/(1.+3.*psquared)*exp((2.*(1.-x))/(1.+3.*psquared));
            emisfunc *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * b / CGS::mec2;
        }

        if (!std::isfinite( emisfunc) ){
            int z = 1;
        }

        return emisfunc;
    }
#endif
    static double synchBessel(double B, double nu, double gam){
//        double t = (4.*M_PI*CGS::me*CGS::c*nu) / (3.*CGS::qe*B*gam*gam);
//        double nu_s = (3.*CGS::qe*B*gam*gam) / (4.*M_PI*CGS::me*CGS::c);
        double t = nu / nu_crit(B,gam);
//        double p_norm = (2*std::sqrt(3))*CGS::qe*CGS::qe*CGS::qe*B/(CGS::me*CGS::c*CGS::c);
        double res_;
        try {
            double bessel_43 = std::cyl_bessel_k(4. / 3., t);
            double bessel_13 = std::cyl_bessel_k(1. / 3., t);
            res_ = 2*t*t * (bessel_43 * bessel_13 - 3. / 5 * t * (bessel_43*bessel_43 - bessel_13*bessel_13));
        } catch (const std::runtime_error & ){
            res_ = 0;
        }
        // factor 2 is missing...
        res_ *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / (CGS::me*CGS::c*CGS::c); // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)
        return res_;
    }

    static double cycSynch(double B, double nu, double gam){
        double emisfunc;
        if (gam > 2.) {
//            emisfunc = synchGSL(B, nu, gam);
            /// ratio of the frequency to the critical synchrotron frequency from Eq. 7.34 in DermerMenon2009, argument of R(x),
//            double nu_s = (3.*CGS::qe*B*gam*gam) / (4.*M_PI*CGS::me*CGS::c); // scale synchrotron frequency
            double x = nu / nu_crit(B,gam);
            /// Aharonian2010 Eq. D7
            double term_1_num = 1.808 * std::pow(x, 1. / 3.);
            double term_1_denom = std::sqrt(1 + 3.4 * std::pow(x, 2. / 3.));
            double term_2_num = 1. + 2.21 * std::pow(x, 2. / 3.) + 0.347 * std::pow(x, 4. / 3.);
            double term_2_denom = 1. + 1.353 * std::pow(x, 2. / 3.) + 0.217 * std::pow(x, 4. / 3.);
            emisfunc = term_1_num / term_1_denom * term_2_num / term_2_denom * std::exp(-x);
            emisfunc *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / CGS::mec2;
        } else { //cyclotron regime
            /// Petrosian (1981) ; Chisellini 1988 Eq.8 ; Lucchini Eq.16
//            double nu_larmor = (CGS::qe*B)/(2.*pi*CGS::me*CGS::c);
            double x = nu/ nu_larmor(B,gam);
            double psquared = gam*gam - 1.;
            emisfunc = (2.*psquared)/(1.+3.*psquared)*exp((2.*(1.-x))/(1.+3.*psquared)); // Eq.8 in Chiselini+98
            emisfunc *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / CGS::mec2;
        }
        return emisfunc;
    }

    /**
     * ratio of the frequency to the critical synchrotron frequency from
     * Eq. 7.34 in [DermerMenon2009]_, argument of R(x),
     * note B has to be in cgs Gauss units
     *
     * Eq. 7.45 in [Dermer2009]_, angle-averaged integrand of the radiated power, the
     * approximation of this function, given in Eq. D7 of [Aharonian2010]_, is used.
     * @param B
     * @param en
     * @param gam
     * @return
     */
    static double synchDermer(double B, double nu, double gam){

        // /nu_c = (3.*CGS::qe*b*pow(gamma,2.))/(4.*pi*CGS::me*CGS::c);
//        double nu = en * CGS::ergToFreq;
//        double en = nu * CGS::h_by_mec2;

//        double x = 4. * M_PI * en * (CGS::me*CGS::me) * (CGS::c3) / (3. * CGS::qe * CGS::h * B * gam*gam);

//        double x = 4. * M_PI * nu * CGS::me * CGS::c / (3. * CGS::qe * B * gam*gam);
        // double t = (2.*M_PI*CGS::me*CGS::c*nu) / (3.*CGS::qe*B*gam*gam);
//        double nu_s = (3.*CGS::qe*B*gam*gam) / (4.*M_PI*CGS::me*CGS::c);
//        double nu_s = (CGS::qe*B*gam*gam) / (4.*M_PI*CGS::me*CGS::c);
        double x = nu / nu_crit(B,gam);

        double term_1_num = 1.808 * std::pow(x, 1. / 3.);
        double term_1_denom = std::sqrt(1 + 3.4 * std::pow(x, 2. / 3.));
        double term_2_num = 1. + 2.21 * std::pow(x, 2. / 3.) + 0.347 * std::pow(x, 4. / 3.);
        double term_2_denom = 1. + 1.353 * std::pow(x, 2. / 3.) + 0.217 * std::pow(x, 4. / 3.);
        double res = term_1_num / term_1_denom * term_2_num / term_2_denom * std::exp(-x);

        res *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / (CGS::me*CGS::c*CGS::c); // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)

        return res;
#if 0
        double res_gsl = synchGSL(B,nu,gam);
        double res_CySyn = cycSynch(B,nu,gam);
        double res_bess = synchBessel(B,nu,gam);

        std::cout << "g="<< gam << " B="<< B << " nu="<<nu << " | x=" << x
                  << " der="<<res<<" gsl="<<res_gsl<<" cycl="<<res_CySyn<<" bess="<<res_bess<<"\n";
        if (gam < 2 and res_CySyn > 0){
            int z = 1;
        }
        return res_CySyn;

//
////        return res;
//
//        double t = (2.*M_PI*CGS::me*CGS::c*nu) / (3.*CGS::qe*B*gam*gam);
//
////        std::cout << "g="<< gam << " B="<< B << " nu="<<nu << " | x=" << x << " | t=" <<t ;
//
//
//        double p_norm = (2*std::sqrt(3))*CGS::qe*CGS::qe*CGS::qe*B/(CGS::me*CGS::c*CGS::c);
//        double res_;
//        try {
//            double bessel_43 = std::cyl_bessel_k(4. / 3., t);
//            double bessel_13 = std::cyl_bessel_k(1. / 3., t);
//            res_ = 2*t*t * (bessel_43 * bessel_13 - 3. / 5 * t * (bessel_43*bessel_43 - bessel_13*bessel_13));
//        } catch (const std::runtime_error & ){
//            res = 0;
//        }
//        res_ *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / (CGS::me*CGS::c*CGS::c); // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)
//
//        return res_;
//        if (!std::isfinite(res_)){
//            int x_ = 1;
//        }
//        double res___ = cycSyn(B, nu, gam);
//        double res__ = synchGSL(B,nu,gam);
////        res___ *= std::sqrt(3.) * std::pow(CGS::qe, 3.) * B / CGS::mec2; // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)
//
////        if (gam < 2 and res > 0){
////            int zz = 0;
////        }
//
//        std::cout<< " | res="  << res
//                  << " | /bes=" << res/res_<<" | /cs=" << res/res___ << " | /gsl=" << res/ res__ << "\n";
//
//        return res;
#endif
    }

    /**
     * Fill the kernel array with values for this value of magnetic field
     * @param B
     */
    void evalSynKernel(Vector & gams, Vector & freqs, double B){
        for (size_t k = 0; k < freqs.size(); k++){
            for (size_t j = 0; j < gams.size(); j++){
                kernel[k][j] = m_func(B, freqs[k], gams[j]);
            }
        }
    }

    /*
     * // kernel[i_freq, i_gam]
     */
    VecVector & getKernel(){return kernel;}

};

#endif //SRC_KERNELS_H


