//
// Created by vsevolod on 1/15/24.
//

#ifndef SRC_KERNELS_H
#define SRC_KERNELS_H

#include "../utilitites/utils.h"

class SSCKernel{

    std::vector<VecVector> kernel; // [i_nu_ssc, i_gam, i_nu_syn]
    VecVector kernel_integral;     // [i_gam, i_nu_syn]
    State & ssc; State & ele; State & syn;
    double (*m_func) (double,double,double);

public:

    SSCKernel(State & ele, State & syn, State & ssc, double (*func)(double,double,double))
            :ele(ele),syn(syn),ssc(ssc),m_func(func){}

    void evalSSCkernel(){
//        std::vector<VecVector> kernel{};
        kernel.resize(ssc.numbins);
        /// allocate memory
        for (auto & vecvec : kernel) {
            vecvec.resize(ele.numbins);
            for (auto & vec : vecvec)
                vec.resize(syn.numbins, 0.);
        }
        /// compute kernel
//    auto * eq = ssc_kernel_nava;
        for (size_t i = 0; i < ssc.numbins; i++){
            for (size_t j = 0; j < ele.numbins; j++){
                for (size_t k = 0; k < syn.numbins; k++){
                    kernel[i][j][k] = m_func(ssc.e[i], ele.e[j], syn.e[k]);
                }
            }
        }
//        return std::move( kernel );
    }

    void evalSSCgridIntergral(){

        double h = 6.6261e-27; //  erg*sec

        /// allocate memory
        kernel_integral.resize(ele.numbins);
        for (auto & vec : kernel_integral)
            vec.resize(syn.numbins, 0.);

        /// compute integral over SSC. frequencies
        for (size_t j = 0; j < ele.numbins; j++){
            for (size_t k = 0; k < syn.numbins; k++){
                double res = 0;
                for (size_t i = 0; i < ssc.numbins-1; i++) {
                    res += kernel[i][j][k]
                           * (ssc.e[i+1] / 8.093440820813486e-21) // freq
                           * (ssc.de[i+1] / 8.093440820813486e-21) // dfreq
                           * h;
                }
                kernel_integral[j][k] = res; // [i_gamma, i_energy_syn]
            }
        }
    }

    static double sscNava(double en, double gam, double en_tilde){
        double K = 0.;
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
    static double sscHUANG(double en_ssc, double gam, double en_syn){
        double freq_ssc = en_ssc / 8.093440820813486e-21;
        double freq_syn = en_syn / 8.093440820813486e-21;
        double mec2 = 8.187e-7; // ergs
        double h = 6.6261e-27; // erg*sec
        double g = gam * h * freq_syn / mec2;
        double w = h * freq_ssc / gam / mec2;
        double q = w / (4.*g*(1.-w));
        double freq_ssc_max = gam * mec2 * 4. * g / (4. * g + 1.) / h;
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
    double (*m_func) (double,double,double);
    State & ele; State & syn;
public:
    SynKernel(State & ele, State & syn, double (*func) (double,double,double)) : ele(ele), syn(syn), m_func(func){
        /// allocate memory for the kernel (kernel is not static, but size is)
        kernel.resize(syn.numbins);
        for (auto & arr : kernel)
            arr.resize(ele.numbins);
    }

    /**
     * this implements the synchrotron kernal from GSL without relying on GSL
     * @param B
     * @param en
     * @param gam
     * @return
     */
    static double synchGSL(double B, double en, double gam){

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
        double x = en / (ec * gam * gam);

//        double M_SQRT2 = np.sqrt(2.0);
//        double M_SQRT3 = std::sqrt(3.0);
//        double M_PI = np.pi;
        double GSL_SQRT_DBL_EPSILON = 1.4901161193847656e-08;
        double GSL_LOG_DBL_MIN = -7.0839641853226408e02;

        double c0 = M_PI / M_SQRT3;

        double c01 = 0.2257913526447274323630976;
        double cond1 = 2 * M_SQRT2 * GSL_SQRT_DBL_EPSILON;
        double cond3 = -8.0 * GSL_LOG_DBL_MIN / 7.0 ;

        if (x < cond1){
            double z = std::pow(x, 1.0 / 3.0);
            double cf = 1 - 8.43812762813205e-01 * z * z;
            return 2.14952824153447863671 * z * cf;
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

            return px * result_c1 - px11 * result_c2 - c0 * x;
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
            return std::sqrt(x) * result_c1 * std::exp(c01 - x);
        }
        else
            return 0.0;
    };

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
    static double synchDermer(double B, double en, double gam){
//        double me = 9.1094e-28;// g
//        double c = 2.99792e10; // cm/s
//        double qe = 4.8032e-10; // statcoulonb
//        double h = 6.6261e-27; // erg*sec

        double x = 4. * M_PI * en * (CGS::me*CGS::me) * (CGS::c3)
                   / (3. * CGS::qe * B * CGS::h * gam*gam);

        double term_1_num = 1.808 * std::pow(x, 1. / 3.);
        double term_1_denom = std::sqrt(1 + 3.4 * std::pow(x, 2. / 3.));
        double term_2_num = 1. + 2.21 * std::pow(x, 2. / 3.) + 0.347 * std::pow(x, 4. / 3.);
        double term_2_denom = 1. + 1.353 * std::pow(x, 2. / 3.) + 0.217 * std::pow(x, 4. / 3.);
        double res = term_1_num / term_1_denom * term_2_num / term_2_denom * std::exp(-x);
        return res;
    }

    /**
     * Fill the kernel array with values for this value of magnetic field
     * @param B
     */
    void evalSynKernel(double B){
        for (size_t k = 0; k < syn.numbins; k++){
            for (size_t j = 0; j < ele.numbins; j++){
                kernel[k][j] = m_func(B, syn.e[k], ele.e[j]);
            }
        }
    }

    /*
     * // kernel[i_freq, i_gam]
     */
    VecVector & getKernel(){return kernel;}

};

#endif //SRC_KERNELS_H
