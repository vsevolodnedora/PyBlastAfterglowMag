//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_SYNCHROTRON_AN_H
#define SRC_SYNCHROTRON_AN_H

#include "pch.h"
#include "logger.h"
#include "utils.h"
#include "interpolators.h"
#include "quadratures.h"
#include "base_equations.h"

class RadiationBase{
public:
    enum METHOD_TAU { iAPPROX, iTHICK, iSMOOTH, iSHARP };


    static inline double optical_depth(const double abs_lab, const double dr,
                                       const double mu, const double beta_shock){
        if (abs_lab <= 0.)
            return 0.;
        double dtau;
        if(mu == beta_shock)
            dtau = 1.0e100; // HUGE VAL, just in case
        else
            dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // basic dtau = abs_lab * dr;
        return dtau;
    }

    /**
    * Compute the I = emission/absorption * ( 1 - tau^{-absorption * thickness})
    * NOTE: adapted from afterglopy library
    * @param em_lab : emission
    * @param abs_lab : absorption
    * @param dr : thickness
    * @param mu : cos(angle to line of sight)
    * @param beta_shock : shock velocity
    * @return
    */
    static inline double computeIntensity(const double em_lab, double dtau, METHOD_TAU m_tau){

        if ( dtau == 0 )
            return em_lab;

        double intensity = em_lab;

        // (Signed) Optical depth through this shell.
        // if negative, face is oriented away from observer.
//        double dtau;
//        if(mu == beta_shock)
//            dtau = 1.0e100; // HUGE VAL, just in case
//        else
//            dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // basic dtau = abs_lab * dr;

        // Now that we know the optical depth, we apply it in a way
        // according to the given specType
        if(m_tau == iTHICK) {
            // Special case: use the optically thick limit *everywhere*
            if(dtau <= 0.0)
                intensity = 0.0;
            else
                intensity /= dtau;
        }
        else if (m_tau == iAPPROX){
            double u = 1. / 2. + exp(-dtau) / dtau - (1 - exp(-dtau)) / (dtau * dtau);
            if (dtau < 1e-3) {
                intensity *= 1.;
            }
            else {
                intensity *= (3. * u / dtau);
            }
//            return np.where(tau < 1e-3, 1., 3 * u / tau)
        }
        else if(m_tau == iSMOOTH) {
            // Apply self-absorption "properly"
            //
            // correction factor to emissivity from absorption
            // ( 1 - e^(-tau) ) / tau  (on front face)
            //
            // back face has extra factor ~ e^-beta_shock/(mu-beta_shock)
            //
            // for now ignoring shadowing by the front face.
            double abs_fac;
            if(dtau == 0.0)
                abs_fac = 1.0;
            else if(dtau > 0.0)
                abs_fac = -expm1(-dtau) / dtau;
            else {
                abs_fac = expm1(dtau) / dtau; //* exp(
                //abs * DR * beta_shock*mu / (mu - beta_shock));
            }

            intensity *= abs_fac;
        }
        else if(m_tau == iSHARP){
            // Apply self-absorption "simply".
            //
            // Compute flux in optically thick limit,
            // use in final result if less than optically thin calculation.
            //
            // e.g. use tau->infty limit if tau > 1.0

            // "Forward" face
            if(dtau > 1.0)
                intensity /= dtau;
                // "Back" face --> assume shadowed by front
            else if(dtau < -1.0)
                intensity = 0.0;
        }
        return intensity;
    }
};

/// From https://github.com/bmargalit/thermal-synchrotron
struct Margalit21 {

//        static constexpr size_t n = 7;
//        static inline Array temps(){
//            return { 5e8, 1e9, 2e9, 4e9, 8e9, 1.6e10, 3.2e10 };
//        }
//        static inline Array alphas(){
//            return { 0.0431, 1.121, 1.180, 1.045, 0.9774, 0.9768, 0.9788 };
//        }
//        static inline Array betas (){
//            return {10.44, -10.65, -4.008, -0.1897, 1.160, 1.095, 1.021};
//        }
//        static inline Array gammas(){
//            return {16.61, 9.169, 1.559, 0.0595, 0.2641, 0.8332, 1.031};
//        }
//        static inline Array Thetas(){
//            return temps() * CGS::kB / (CGS::me * CGS::c * CGS::c);
//        }


    /* Utility function defined in eq. (1), MQ21
    This is an approximate fitting function from Gammie & Popham (1998)

    Parameters
     __________
    Theta : float
            The dimensionless electron temperature

            Returns
     _______
            a : float
            value of a(Theta)
    */
    static inline double
    a_fun(const double &Theta) {
        double val = (6.0 + 15.0 * Theta) / (4.0 + 5.0 * Theta);
        return val;
    }

    /*
     * alculate electron temperature

        This calculates the post-shock dimensionless electron temperature following
        eqs. (2,3) of MQ21

        Parameters
        __________
        beta : float
            Shock velocity / c
        mu : float, optional
            Mean molecular weight, default is 0.62
        mu_e : float, optional
            Mean molecular weight per elecron, default is 1.18
        epsilon_T : float, optional
            Electron thermalization efficiency (<=1), default is 1e0

        Returns
        _______
        Theta : float
            Dimensionless electron temperature = kb*Te/me*c**2
     */
    static inline double
    Theta_fun(const double beta, const double mu, const double mu_e, const double epsilon_T) {

        ///  eq. (3) of MQ21:
        double Theta0 = epsilon_T * (9.0 * mu * CGS::mp / (32.0 * mu_e * CGS::me)) * beta * beta;
        /// eq. (2) of MQ21:
        double val = (5.0 * Theta0 - 6.0 + sqrt(25.0 * Theta0 * Theta0 + 180.0 * Theta0 + 36.0)) / 30.0;
        return val;
    }

    /*
    Calculate minimal Lorentz factor of non-thermal electrons; eq. (6), MQ21

    Parameters
            __________
    Theta : float
            Dimensionless electron temperature

    Returns
            _______
    gamma_m : float
            minimum Lorentz factor of power-law electrons
    */
    static inline double
    gamma_m_fun(const double Theta) {
        return 1. + a_fun(Theta) * Theta; // 1e0 +
    }

    /*
    Utility function f(Theta); eq. (5) of MQ21

    Parameters
            __________
    Theta : float
            Dimensionless electron temperature

    Returns
            _______
    f : float
            Correction term to the thermal electron distribution that is relevant
            only in the non-relativistic regime (Theta ~< 1)
    */
    static inline double
    f_fun(const double Theta) {
        //  Computes the irregular modified cylindrical Bessel function
        // (also known as modified Bessel function of the second kind) of ν and x.
        double tmp = 0;
        try {
            tmp = 2.0 * Theta * Theta / std::cyl_bessel_k(2, 1.0 / Theta);
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        return tmp;
    }

    /*
    Utility function g(Theta); eq. (8) of MQ21

    Parameters
            __________
    Theta : float
            Dimensionless electron temperature
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0

    Returns
            _______
    g : float
            Correction term to the power-law electron distribution that is relevant
    only in the non-relativistic regime (Theta ~< 1)
    */
    static inline double
    g_fun(const double Theta, double gamma_m, const double p) {
//            const double gamma_m = gamma_m_fun(Theta);
        const double val = ((p - 1.0) * gamma_m / ((p - 1.0) * gamma_m - p + 2.0))
                           * std::pow(gamma_m / (3.0 * Theta), (p - 1.0))
        ;
//                    * pow(1./3./Theta, (p - 1.0));
        return val;
    }

    /*
    Function I'(x) derived by Mahadevan et al. (1996)

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)

    Returns
            _______
    I : float
            Spectral energy distribution function (eq. 13, MQ21)
    */
    static inline double
    I_of_x(const double x, const double Theta) {

        /// Mahadevan+95 "Harmony in Electrons"
        const size_t n = 7;
        Array temps = { 5e8, 1e9, 2e9, 4e9, 8e9, 1.6e10, 3.2e10 };
        Array alphas = { 0.0431, 1.121, 1.180, 1.045, 0.9774, 0.9768, 0.9788 };
        Array betas = {10.44, -10.65, -4.008, -0.1897, 1.160, 1.095, 1.021};
        Array gammas = {16.61, 9.169, 1.559, 0.0595, 0.2641, 0.8332, 1.031};

        Array Thetas = temps * CGS::kB / (CGS::me * CGS::c * CGS::c);

        double alpha, beta, gamma;
        if (Theta < Thetas[0]){
            alpha = alphas[0]; beta=betas[0]; gamma=gammas[0];
        }
        else if (Theta > Thetas[n-1]){
            alpha=alphas[n-1]; beta=betas[n-1]; gamma=gammas[n-1];
        }
        else{
            size_t ia = findIndex(Theta, temps, n);
            size_t ib = ia + 1;
            alpha = interpSegLin(ia, ib, Theta, Thetas, alphas);
            beta = interpSegLin(ia, ib, Theta, Thetas, betas);
            gamma = interpSegLin(ia, ib, Theta, Thetas, gammas);
        }


        const double i = 4.0505 * alpha * std::pow(x, -1.0 / 6.0)
                         * (1.0 + 0.40 * beta * std::pow(x, -0.25)
                            + 0.5316 * gamma * std::pow(x, -0.5))
                         * exp(-1.8899 * std::pow(x, 1.0 / 3.0));
        return i;


//            const double i_ = 4.0505 * pow(x, -1.0 / 6.0) * (1.0 + 0.40 * pow(x, -0.25) + 0.5316 * pow(x, -0.5))
//                             * exp(-1.8899 * pow(x, 1.0 / 3.0));
//            return i_;
    }

    /*
    Prefactor to power-law synchrotron emissivity (eq. 15, MQ21)

    Parameters
            __________
    p_xi : float
            Slope of power-law electron distribution

            Returns
    _______
            Cj : float
            Synchrotron constant
    */
    static inline double
    C_j(const double p) {
        double tmp = 0;
        try {
            tmp = (tgamma((p + 5.0) / 4.0) / tgamma((p + 7.0) / 4.0)) *
                  tgamma((3.0 * p + 19.0) / 12.0) * tgamma((3.0 * p - 1.0) / 12.0) *
                  ((p - 2.0) / (p + 1.0)) * std::pow(3.0, (2.0 * p - 1.0) / 2.0)
                  * std::pow(2.0, -(7.0 - p) / 2.0) * std::pow(M_PI, -0.5);
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        return tmp;
    }

    /*
    Prefactor to power-law synchrotron absorption coefficient (eq. 17, MQ21)

    Parameters
            __________
    p_xi : float
            Slope of power-law electron distribution

            Returns
    _______
            Calpha : float
            Synchrotron constant
    */
    static inline double
    C_alpha(const double p) {
        double tmp = 0;
        try {
            tmp = (tgamma((p + 6.0) / 4.0) / tgamma((p + 8.0) / 4.0)) *
                  tgamma((3.0 * p + 2.0) / 12.0) * tgamma((3.0 * p + 22.0) / 12.0) * (p - 2.0)
                  * std::pow(3.0, (2.0 * p - 5.0) / 2.0) * std::pow(2.0, p / 2.0) * std::pow(M_PI, 3.0 / 2.0);
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        return tmp;
    }

    /*
    Low-frequency correction to power-law emissivity

    This function returns a multiplicative frequency-dependent correction term
            that modifies the high-frequency power-law emissivity in the regime
            nu <~ nu_m (where nu_m is the synchrotron frequency corresponding to the
    minimum Lorentz factor of power-law electrons, gamma_m). We adopt an
            approximate expression that smoothly interpolates between the exact high-
    and low-frequency results.

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    Theta : float
            Dimensionless electron temperature
    p_xi : float
            Slope of power-law electron distribution

            Returns
    _______
            val : float
            Correction term
    */
    static inline double
    low_freq_jpl_correction(double gamma_m, const double x, const double Theta, const double p) {
//            double gamma_m = gamma_m_fun(Theta);
        /// synchrotron constant in x<<x_m limit

        double Cj_low = 0;
        try {
            Cj_low = -std::pow(M_PI, 1.5) * (p - 2.0) /
                     (std::pow(2.0, 1.0 / 3.0) * std::pow(3.0, 1.0 / 6.0) * (3.0 * p - 1.0) * tgamma(1.0 / 3.0) *
                      tgamma(-1.0 / 3.0) * tgamma(11.0 / 6.0));
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        /// multiplicative correction term
        double corr = (Cj_low / C_j(p)) * std::pow(gamma_m / (3.0 * Theta), -(3.0 * p - 1.0) / 3.0) *
                      std::pow(x, (3.0 * p - 1.0) / 6.0);
        /// approximate interpolation with a "smoothing parameter" = s
        double s = 3.0 / p;
        double val = std::pow(1e0 + std::pow(corr, -s), -1.0 / s);
        return val;
    }

    /*
    Low-frequency correction to power-law absorption coefficient

    This function returns a multiplicative frequency-dependent correction term
            that modifies the high-frequency power-law absorption coeff in the regime
    nu <~ nu_m (where nu_m is the synchrotron frequency corresponding to the
    minimum Lorentz factor of power-law electrons, gamma_m). We adopt an
            approximate expression that smoothly interpolates between the exact high-
    and low-frequency results.

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    Theta : float
            Dimensionless electron temperature
    p_xi : float
            Slope of power-law electron distribution

            Returns
    _______
            val : float
            Correction term
    */
    static inline double
    low_freq_apl_correction(const double x, const double Theta, const double p) {

        double gamma_m = gamma_m_fun(Theta);
        /// synchrotron constant in x<<x_m limit
        double Calpha_low = 0;
        try {
            Calpha_low = -std::pow(2.0, 8.0 / 3.0) * std::pow(M_PI, 7.0 / 2.0) * (p + 2.0) * (p - 2.0) /
                         (std::pow(3.0, 19.0 / 6.0) * (3.0 * p + 2) * tgamma(1.0 / 3.0) * tgamma(-1.0 / 3.0) *
                          tgamma(11.0 / 6.0));
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        /// multiplicative correction term
        double corr = (Calpha_low / C_alpha(p)) * std::pow(gamma_m / (3.0 * Theta), -(3.0 * p + 2.0) / 3.0) *
                      std::pow(x, (3.0 * p + 2.0) / 6.0);
        /// approximate interpolation with a "smoothing parameter" = s
        double s = 3.0 / p;
        double val = std::pow(1e0 + std::pow(corr, -s), -1.0 / s);
        return val;
    }

    /*Synchrotron emissivity of power-law electrons (eqs. 14,19; MQ21)

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    n : float
            Electron number density in the emitting region (in cm^{-3})
    B : float
            Magnetic field strength (in G)
    Theta : float
            Dimensionless electron temperature
    delta : float, optional
            Fraction of energy carried by power-law electrons, default is 1e-1
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0
    z_cool : float, optional
            Normalized cooling Lorentz factor = gamma_cool/Theta (eq. 18, MQ21),
    default is np.inf (negligible cooling)

    Returns
            _______
    val : float
            Synchrotron emissivity
    */
    static inline double
    jnu_pl(const double x, const double n, const double B, const double Theta, const double gamma_m,
           const double delta, const double p, const double z_cool) {


//            double gamma_m = gamma_m_fun(Theta);
        double val = C_j(p) * (CGS::qe * CGS::qe * CGS::qe / (CGS::me * CGS::c * CGS::c))
                     * delta * n * B
                     * g_fun(Theta, gamma_m, p)
                     * std::pow(x, -(p - 1.0) / 2.0);
        /// correct emission at low-frequencies x < x_m:
        val *= low_freq_jpl_correction(gamma_m, x, Theta, p);
        /// fast-cooling correction:
        double z0 = std::pow(x, 0.5);
        val *= std::min(1e0, 1./(z0 / z_cool));
        return val;
    }

    /*
    Synchrotron absorption coeff of power-law electrons (eqs. 16,19; MQ21)

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    n : float
            Electron number density in the emitting region (in cm^{-3})
    B : float
            Magnetic field strength (in G)
    Theta : float
            Dimensionless electron temperature
    delta : float, optional
            Fraction of energy carried by power-law electrons, default is 1e-1
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0
    z_cool : float, optional
            Normalized cooling Lorentz factor = gamma_cool/Theta (eq. 18, MQ21),
    default is np.inf (negligible cooling)

    Returns
            _______
    val : float
            Synchrotron absorption coefficient
    */
    static inline double
    alphanu_pl(const double x, const double n, const double B,
               const double Theta, const double gamma_m, const double delta, const double p,
               const double z_cool) {

//            double gamma_m = gamma_m_fun(Theta);
        double val = C_alpha(p)
                     * CGS::qe
                     * (delta * n / (std::pow(Theta, 5) * B))
                     * g_fun(Theta, gamma_m, p)
                     * std::pow(x, -(p + 4.0) / 2.0);
        /// correct emission at low-frequencies x < x_m:
        val *= low_freq_apl_correction(x, Theta, p);
        /// fast-cooling correction:
        double z0 = std::pow(x, 0.5);
        val *= std::min(1e0, 1./(z0 / z_cool));
        return val;
    }

    /*
    Synchrotron emissivity of thermal electrons (eqs. 10,20; MQ21)

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    n : float
            Electron number density in the emitting region (in cm^{-3})
    B : float
            Magnetic field strength (in G)
    Theta : float
            Dimensionless electron temperature
    z_cool : float, optional
            Normalized cooling Lorentz factor = gamma_cool/Theta (eq. 18, MQ21),
    default is np.inf (negligible cooling)

    Returns
            _______
    val : float
            Synchrotron emissivity
    */
    static inline double
    jnu_th(const double x, const double n, const double B, const double Theta, const double z_cool){

        double val = (sqrt(3.0) / (8.0 * M_PI)) * (CGS::qe*CGS::qe*CGS::qe / (CGS::me * CGS::c * CGS::c))
                     * f_fun(Theta) * n * B * x * I_of_x(x, Theta);
        /// fast-cooling correction:
        double z0 = std::pow(2.0 * x, 1.0 / 3.0);
        val *= std::min( 1e0, 1./(z0/z_cool));
        return val;
    }
    /*
    Synchrotron absorption coeff of thermal electrons (eqs. 12,20; MQ21)

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    n : float
            Electron number density in the emitting region (in cm^{-3})
    B : float
            Magnetic field strength (in G)
    Theta : float
            Dimensionless electron temperature
    z_cool : float, optional
            Normalized cooling Lorentz factor = gamma_cool/Theta (eq. 18, MQ21),
    default is np.inf (negligible cooling)

    Returns
            _______
    val : float
            Synchrotron absorption coefficient
    */
    static double alphanu_th(const double x, const double n, const double B, const double Theta, const double z_cool) {
        const double i_x = I_of_x(x, Theta);
        if (i_x == 0.) // TODO check why this gives 0 and add a better constraint on which Ejecta produce termal emission
            return 0.;
        double val = (M_PI * std::pow(3.0, -3.0 / 2.0)) * CGS::qe * (n / (std::pow(Theta,5) * B)) * f_fun(Theta) * (1./x) * i_x;
        if (!std::isfinite(val)){
            std::cerr << AT <<" \n nan in alphanu_th() " << " x=" << x << " B="<<B<<" Theta="<<Theta<<" z_cool="<<z_cool << "\n";
//                exit(1);
            return 0.;
        }
        /// fast-cooling correction:
        double z0 = std::pow(2.0 * x, 1.0 / 3.0);
        val *= std::min(1e0, 1./(z0 / z_cool));
        return val;
    }
    /*
    Total (thermal+non-thermal) synchrotron optical depth

    Parameters
            __________
    x : float
            Dimensionless frequency = nu/nu_Theta (see eq. 11, MQ21)
    n : float
            Electron number density in the emitting region (in cm^{-3})
    R : float
            Characteristic size of the emitting region (in cm)
    B : float
            Magnetic field strength (in G)
    Theta : float
            Dimensionless electron temperature
    delta : float, optional
            Fraction of energy carried by power-law electrons, default is 1e-1
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0
    z_cool : float, optional
            Normalized cooling Lorentz factor = gamma_cool/Theta (eq. 18, MQ21),
    default is np.inf (negligible cooling)

    Returns
            _______
    val : float
            Synchrotron absorption coefficient
    */
    static inline double
    tau_nu(const double x, const double n, const double R, const double B, const double Theta, const double gamma_m,
           const double delta, const double p, const double z_cool) {

        double val = R * (alphanu_th(x, n, B, Theta, z_cool) +
                          alphanu_pl(x, n, B, Theta, gamma_m, delta, p, z_cool));
        return val;
    }
    /*
    Characteristic "thermal" synchrotron frequency (eq. 11, MQ21)

    Parameters
            __________
    Theta : float
            Dimensionless electron temperature
    B : float
            Magnetic field strength (in G)

    Returns
            _______
    val : float
            Synchrotron frequency
    */
    static inline double
    nu_Theta(const double Theta, const double B) {

        double val = 3.0 * Theta * Theta * CGS::qe * B / (4.0 * M_PI * CGS::me * CGS::c);
        return val;
    }
    /*
    Synchrotron specific luminosity

    This function calculates the emergent synchrotron luminosity at frequency
    nu (including synchrotron self-absorption and synchrotron cooling effects)
    as a function of the physical parameters within the emitting region. This
            form is applicable also outside the scope of shock-powered transients.

    Parameters
            __________
    nu : float
    Frequency (in Hz)
    ne : float
            Electron number density in the emitting region (in cm^{-3})
    R : float
            Characteristic size of the emitting region (in cm)
    Theta : float
            Dimensionless electron temperature
    B : float
            Magnetic field strength (in G)
    t : float
            Dynamical time = R/vv ~ source age (in s)
    delta : float, optional
            Fraction of energy carried by power-law electrons, default is 1e-1
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0
    mu : float, optional
            Mean molecular weight, default is 0.62
    mu_e : float, optional
            Mean molecular weight per elecron, default is 1.18

    Returns
            _______
    val : float
            Specific luminosity (in erg/s/Hz)
    */
    static double
    Lnu_fun(const double nu, const double ne, const double R, const double Theta, const double gamma_m,
            const double B, const double t, const double delta, const double p,
            const double mu, const double mu_e) {


        /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21):
        double z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * B * B * t)) / Theta;
        /// normalized frequency:
        double x = nu / nu_Theta(Theta, B);
        /// calculate total emissivity & optical depth:
        double j = jnu_th(x, ne, B, Theta, z_cool) + jnu_pl(x, ne, B, Theta, gamma_m, delta,  p, z_cool);
        /// calculate optical depth
//        double tau = tau_nu(x, ne, R, B, Theta, delta, p_xi, z_cool);
        double tau = R * (alphanu_th(x, ne, B, Theta, z_cool) + alphanu_pl(x, ne, B, Theta, gamma_m, delta, p, z_cool));
        /// calculate specific luminosity Lnu (eq. 21, MQ21):
        double val = 4.0 * M_PI * M_PI * R*R*R * j * (1e0 - exp(-tau)) / tau;
        /// prevent roundoff errors at tau<<1:
        if (tau < 1e-10)
            val = 4.0 * M_PI * M_PI * R*R*R * j;
        return val;
    }
    /*
    Synchrotron specific luminosity as a function of shock parameters

    This function calculates the emergent synchrotron luminosity at frequency
    nu (including synchrotron self-absorption and synchrotron cooling effects)
    within the context of a shock-powered model: the post-shock magnetic field,
    electron temperature, and electron number-density are related to the shock
    velocity and upstream (ambient medium) density using epsilon_B & epsilon_T.

    Parameters
            __________
    nu : float
    Frequency (in Hz)
    n_ism : float
            *Upstream* number density (in cm^{-3})
    R : float
            Characteristic size of the emitting region (in cm)
    vv : float
            Shock velocity (in g/s)
    epsilon_T : float, optional
            Electron thermalization efficiency, default is 1e0
    epsilon_B : float, optional
            Magnetic field amplification efficiency, default is 1e-1
    delta : float, optional
            Fraction of energy carried by power-law electrons, default is 1e-1
    p_xi : float, optional
            Slope of power-law electron distribution, default is 3.0
    mu : float, optional
            Mean molecular weight, default is 0.62
    mu_e : float, optional
            Mean molecular weight per elecron, default is 1.18

    Returns
            _______
    val : float
            Specific luminosity (in erg/s/Hz)
    */
    static double
    Lnu_shock(const double nu, const double n_ism, const double R, const double v, const double epsilon_T,
              const double epsilon_B, const double delta, const double p,
              const double mu, const double mu_e) {

        /// downstream (post-shock) electron number density:
        double ne = 4.0 * mu_e * n_ism;
        /// downstream electron temperature:
        double Theta = Theta_fun(v / c, mu, mu_e, epsilon_T);
        double gamma_m = gamma_m_fun(Theta);
        /// shock-amplified magnetic field (eq. 9; MQ21):
        double B = sqrt(9.0 * M_PI * epsilon_B * n_ism * mu * CGS::mp) * v;
        /// mean dynamical time:
        double t = R / v;
        /// calculate luminosity:
        double val = Lnu_fun(nu, ne, R, Theta, gamma_m, B, t, delta, p, mu, mu_e);
        return val;
    }
};

/// from Lamb+2818 model (from Fernandez GitHub repo)
double interpolated_xi(const double &p){

    Array p_xi = {
            1.44449807, 1.48206417, 1.48206417, 1.49279734, 1.50353051,
            1.51963027, 1.54109661, 1.55361864, 1.58402929, 1.61944876,
            1.60012905, 1.64842832, 1.66989466, 1.70746076, 1.71282735,
            1.73429369, 1.76327325, 1.79869271, 1.82552564, 1.84699198,
            1.88992467, 1.92212418, 1.95432369, 1.97579004, 2.01872272,
            2.05628882, 2.0992215 , 2.14215419, 2.19045346, 2.23875273,
            2.29778517, 2.36003756, 2.42121664, 2.48561566, 2.55001469,
            2.60904713, 2.67881274, 2.74213845, 2.80761079, 2.8709365 ,
            2.93748216, 3.00402782, 3.06198695, 3.13711915, 3.19937154,
            3.2659172 , 3.3254863 , 3.39471525, 3.45428435, 3.52780657,
            3.60186545, 3.66626448, 3.73066351, 3.78969595, 3.85838824,
            3.93352044, 4.00435937, 4.06875839, 4.13852401, 4.20184971,
            4.27268864, 4.33708767, 4.40685328, 4.47661889, 4.54101792,
            4.61078353, 4.67625588, 4.74494817, 4.80398061, 4.87911281,
            4.93277866, 5.00254428, 5.09377623, 5.19037477, 5.28697331,
            5.38357185, 5.48017039, 5.57676893, 5.67336747, 5.76996601,
            5.86656455, 5.95242991, 5.9577965
    };

    Array Xi = {
            0.99441363, 0.99511912, 0.97241493, 0.95724356, 0.94207218,
            0.92309915, 0.90708165, 0.89361622, 0.8763981 , 0.84773936,
            0.86391198, 0.8286159 , 0.81232812, 0.7972488 , 0.78083371,
            0.7650865 , 0.74646758, 0.73108898, 0.71613501, 0.70336097,
            0.68853431, 0.67473396, 0.66093361, 0.6474388 , 0.63513483,
            0.622398  , 0.60928316, 0.59724948, 0.5847657 , 0.5697232 ,
            0.55760059, 0.54334114, 0.5303826 , 0.51719725, 0.50509306,
            0.49560125, 0.48484889, 0.47591959, 0.46802836, 0.46036041,
            0.4506277 , 0.44403032, 0.43897477, 0.43224107, 0.4270633 ,
            0.4210065 , 0.41605377, 0.40991608, 0.40481919, 0.39895571,
            0.39424369, 0.39060851, 0.38715353, 0.38344588, 0.37884159,
            0.3742702 , 0.37030754, 0.36766342, 0.36473138, 0.36236108,
            0.35942552, 0.35705168, 0.35493051, 0.35118761, 0.34881378,
            0.34588174, 0.34328815, 0.34030559, 0.33821967, 0.33617096,
            0.33293143, 0.33032374, 0.328078  , 0.32577859, 0.3231188 ,
            0.32027882, 0.31815961, 0.31568001, 0.31320041, 0.31072081,
            0.30824122, 0.68419991, 0.30534679
    };

    return ( Interp1d(p_xi, Xi) ).Interpolate(p, Interp1d::iLagrangeLinear);
    //return LinInterp(m_Xp.p_xi, m_Xp.Xi, p_xi,  false);
}
double interpolated_phi(const double &p){

    Array p_phi = {
            1.00621662, 1.04516293, 1.13098906, 1.21681519, 1.25720396,
            1.35312728, 1.46419639, 1.58608392, 1.62070287, 1.69643181,
            1.81759811, 1.85798688, 1.9589588 , 2.09022229, 2.23158298,
            2.27702034, 2.38520454, 2.42342963, 2.53449874, 2.5849847 ,
            2.68595662, 2.73644258, 2.83034646, 2.88285186, 2.99896957,
            3.04440693, 3.16052464, 3.205962  , 3.31198252, 3.36751707,
            3.47353759, 3.51897495, 3.63509266, 3.68053002, 3.79664773,
            3.8420851 , 3.94810561, 3.99859157, 4.10966068, 4.15509805,
            4.27121575, 4.31665312, 4.42267363, 4.48022763, 4.58639239,
            4.63976326, 4.73871574, 4.79122114, 4.89219306, 4.94873733,
            5.02682228, 5.0890883 , 5.14462286, 5.25064337, 5.28598354,
            5.39200406, 5.44249002, 5.51821896, 5.57375351, 5.60909369,
            5.7151142 , 5.7726682 , 5.87666927, 5.93725242
    };
    Array phi = {
            0.41141367, 0.41748522, 0.43521089, 0.45395245, 0.46103901,
            0.48129824, 0.50341068, 0.52512949, 0.53192933, 0.5467828 ,
            0.56618001, 0.5703882 , 0.58454897, 0.60117452, 0.61418181,
            0.61787897, 0.62681073, 0.62896427, 0.63753144, 0.64207209,
            0.6459046 , 0.6483288 , 0.65228415, 0.65489336, 0.66095157,
            0.6610931 , 0.66491634, 0.66708966, 0.6690566 , 0.67080045,
            0.67386795, 0.67502537, 0.67698614, 0.67797425, 0.68010433,
            0.68092312, 0.68289006, 0.684637  , 0.68482304, 0.68581115,
            0.68675602, 0.68672823, 0.68903381, 0.69069177, 0.69280375,
            0.69221962, 0.69358135, 0.69375242, 0.69399543, 0.69406244,
            0.6961142 , 0.69472159, 0.69638078, 0.69648525, 0.69629432,
            0.69639879, 0.69670654, 0.69784544, 0.69730352, 0.69812849,
            0.69806364, 0.69823162, 0.69796483, 0.69792778
    };

    return ( Interp1d(p_phi, phi) ).Interpolate(p, Interp1d::iLagrangeLinear);
//        return LinInterp(m_Phip.p_xi, m_Phip.phi, p_xi,  false);
}

struct Dermer09{
    /// for Dermer+09 model
    static double brokenPowerLaw(const double gamma_e, const double gm, const double gc, const double gM,
                                 const double p1, const double p2){
        double spectrum = 0.;
        double index;
        if (gm <= gamma_e && gamma_e <= gM){
            if (gamma_e <= gc){
                index = p1; // slow cooling
            } else {
                index = p2; // fast cooling
            }
            spectrum = std::pow(gamma_e / gc, -1.0 * index);
        }
        return spectrum;
    }
    /*
     * (analytical) integrand for the synchrotron_old self-absorption:
     *  :math:`\gamma'^2 \frac{d}{d \gamma'} \left(\frac{n_e(\gamma)}{\gamma'^2}\right)`
     */
    static double brokenPowerLawSSA(const double gamma_e, const double gm, const double gc, const double gM,
                                    const double p1, const double p2){
        double spectrum = 0.;
        double index;
        if (gm <= gamma_e && gamma_e <= gM){
            if (gamma_e <= gc){
                index = p1; // slow cooling
            } else {
                index = p2; // fast cooling
            }
            spectrum = -(index + 2) / gamma_e * std::pow(gamma_e / gc, -1.0 * index);
        }
        return spectrum;
    }
    /**
    * observed peak frequency for monoenergetic electrons Eq. 7.19 in [DermerMenon2009]_
    *
    * @param B
    * @param gamma
    * @return
    */
    static double nu_synch_peak(const double &B, const double &gamma){
        return (CGS::qe * B / (2 * CGS::pi * CGS::me * CGS::c)) * (gamma * gamma);
    }
    /**
    * ratio of the frequency to the critical synchrotron frequency from
    * Eq. 7.34 in [DermerMenon2009]_, argument of R(x),
    *
    * @param B
    * @param epsilon
    * @param gamma
    * @return
    */
    static double calc_x(const double &B, const double &epsilon, const double &gamma){
        double x = 4.0 * CGS::pi * epsilon * (CGS::me * CGS::me) * (CGS::c * CGS::c * CGS::c) /
                   (3.0 * CGS::qe * B * CGS::h * (gamma * gamma));
        return x;
    }
    /**
     * Eq. 7.45 in [Dermer2009]_, angle-averaged integrand of the radiated power, the
     * approximation of this function, given in Eq. D7 of [Aharonian2010]_, is used.
     *
     * @param x = ratio of the frequency to the critical synchrotron frequency
     * @return
     */
    static double R(const double &x){

        double term_1_num = 1.808 * std::pow(x, 1.0 / 3.0);
        double term_1_denom = sqrt(1 + 3.4 * std::pow(x, 2.0 / 3.0));
        double term_2_num = 1.0 + 2.21 * std::pow(x, 2.0 / 3.0) + 0.347 * std::pow(x, 4.0 / 3.0);
        double term_2_denom = 1.0 + 1.353 * std::pow(x, 2.0 / 3.0) + 0.217 * std::pow(x, 4.0 / 3.0);

        return term_1_num / term_1_denom * term_2_num / term_2_denom * exp(-x);
    }
    /**
     * angle-averaged synchrotron power for a single electron,
     * to be folded with the electron distribution
     *
     * @param B = magnetic field strength
     * @param epsilon = nuprim * h / (m_e * c^2) where nuprim is the freq. in comoving frame
     * @param gamma_arr = electron energy distribution in terms of lorentz factors
     */
    static double single_electron_synch_power(const double &B, const double &epsilon, const double &gamma){
        double result = 0.;
        double prefactor = sqrt(3.0) * (CGS::qe * CGS::qe * CGS::qe) * B / CGS::h;
        double x = calc_x(B, epsilon, gamma);
        result = prefactor * R(x);
        return result;
    }
    double pprime(const double nuprim, const double p, double gm, const double gc, const double B, double nprim){

        double gM = 1e10;//Electrons::gamma_max(B);
        double epsilon = nuprim * CGS::h / CGS::mec2;

        // electron distribution, broken power law with two regimes depending on the order of gamma_i
        auto integrand_ele = [&](double gam){
            double p1 = gm < gc ? p : 2.0;
            double p2 = p + 1.0;
            double gmax = gM;
            double gmin = gm < gc ? gm : gc;
            double gb = gm < gc ? gc : gm;
            return brokenPowerLaw(gam, gmin, gb, gmax, p1, p2);
        };

        // compute electron distribution normalisation
        double k_e = nprim / Simpson38(gm, gM, 200, integrand_ele); // TODO replace with adative integrals

        // convolve the electron distribution with emission spectrum and itegrate
        auto integrand = [&](double gam){
            double power_e = single_electron_synch_power( B, epsilon, gam );
            double n_e = k_e * integrand_ele(gam);
            return n_e * power_e;
        };
        double power = Simpson38(gm, gM, 200, integrand); // TODO replace with adative integrals
        power *= (CGS::h / CGS::mec2);
//        if (power != power || power < 1e-40){
//            exit(1);
//        }
//        std::cout << "gm=" << gm << " gc=" << gc << " gM=" << gM << " k_e=" << k_e << " power=" << power << " nprim=" << nprim << "\n";
        return power;
    }
};

struct Bretta{
private:
    double B;
    struct Pars{
        int error = 0;
    };
public:

    /// check if during the quadrature integration there was an error
    static int check_error(void *params) {
        auto *fp = (struct Pars *) params; // removing EATS_pars for simplicity
        return fp->error;
//        return 0;
    }

    static double integrant_bessel(double x, void * pars){
        double tmp = 0;
        static double o = 5./3.;
        try {
            tmp = std::cyl_bessel_k(o, x);
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << "\n";
        }
        return tmp;
    }

    /// integrate bessel 5/3 function via adaptive quadrature
    static double bessel(const double xs, void * pars){
        int nmax = 200;
        double rtol = 1e-8;//integrant_bessel(xs, nullptr)/1e7;
        double atol = 0.;//integrant_bessel(xs, nullptr)/1e7;
        const double r = \
                cadre_adapt(&integrant_bessel, xs, 1.e5, // TODO this is ad-hock limit
                            nmax, atol,
                            rtol, pars, nullptr, nullptr,
                            0, check_error, nullptr, nullptr);
        return r * xs;

    }
    static double singleElecSynchroPower(double nu, double gamma, double nuL){
        // From gamma, B -> nu_c; then x_0 = nu / nu_s
        auto p_pars = new Pars;
        const double x0 = (nu/(3./2.*gamma*gamma*(nuL)));
        const double y0 = bessel(x0, p_pars);
        const double bes = y0;
        delete p_pars;
        return bes;
    }
    static double syncEmissivityKernPL(const double gamma, const double nuprim,
                                       const double p, const double nuL, const double gc,
                                       const Vector & gams){
        const double gam_min = gams[0];
        const double gam_max = gams[gams.size()-1];
        /// broken power law
        double ne;
        if ((gam_min <= gamma) && (gamma < gc)){
            ne = std::pow(gamma, -p);
        }
        else if (gamma >= gc) {
            ne = std::pow(gamma, -p-1.);
        }
        else{
            std::cerr << " should not be entered: "
                         " gam_min="<<gam_min<<" gam_c="<<gc<<" gam_max="<<gam_max<<" gam="<<gamma<<"\n";
            exit(1);
        }
        const double k1 = ne * singleElecSynchroPower(nuprim, gamma, nuL);
        return k1;

    }
    static double nu_L(const double b) {
        const double k = (CGS::qe) / (2. * CGS::pi * (CGS::me) * (CGS::c));
        const double nuL = k * b;
        return nuL;
    }
//    void evalSynchEmissivity(double & em, double & atten, const Vector & gammas,
//                             double nuprime, double B, double p, double gc, double n0){
//        const double nuL = Bretta::nu_L(B);
//        const double n1 = 2.*CGS::pi*std::sqrt(3.)*(CGS::qe)*(CGS::qe)*(nuL)/(CGS::c); // Eq, ?? Dermer09?
//        const double k0 = (p+2)/(8*CGS::pi*(CGS::me));
//        double r=0., al=0.;
//        for (size_t igam; igam < gammas.size(); igam++){
//            const double asso = singleElecSynchroPower(nuprime, gammas[igam], nuL)
//                              * std::pow(gammas[igam], -(p+1)); // absorption
//            const double js = syncEmissivityKernPL(gammas[igam], nuprime, p, nuL,gc,gammas); //
//            r += asso; // NOTE! replacing simpson integral with the direct sum!
//            al += js;  // NOTE! replacing simpson integral with the direct sum!
//        }
//        double I=0, att=-1;
//        if (r > 1.e-90){
//            I = n1*r*n0;
//            att = n1*al*k0*n0/std::pow(nuprime,2.);
//        }
//        em = I;
//        atten = att;
//    }
};

/// evaluate comoving emissivity and absorption for synchrotron mech.
class SynchrotronAnalytic : public RadiationBase{
/// methods
    enum METHODS { iWSPN99, iJOH06, iDER06, iMARG21, iBerrettaSynch, iBerettaSynchSSC };
    enum METHODS_LFMIN { igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };
    enum METHODS_SSA { iSSAoff, iSSAon };
    enum METHOD_NONRELDIST{ inone, iuseGm };
    enum QQ { i_em_pl, i_em_th, i_abs_pl, i_abs_th, i_em, i_abs, i_tau, i_int };
    /// parameter container for the class
    struct Pars{
        // --- in
        double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
        double mu=-1, mu_e=-1;
        bool lim_gm_to_1= true;
        double beta_min = -1;
        // --- methods
        METHODS m_sychMethod{};
        METHODS_LFMIN m_methodsLfmin{};
        METHOD_NONRELDIST m_method_nonreldist{};
        RadiationBase::METHOD_TAU method_tau{};
//        METHODS_RAD method_comp_mode{};

        METHODS_SSA m_methods_ssa{};
        QQ m_marg21opt_em = i_em_pl;
        QQ m_marg21opt_abs = i_abs_pl;
        // --- out
        double B=-1, gamma_min=-1, gamma_max=-1, gamma_c=-1, n_prime=-1, accel_frac=-1.;
        double Theta=-1, z_cool=-1, x=-1;
        static std::vector<std::string> listPars(){
            return { "eps_e", "eps_b", "eps_t", "p", "ksi_n", "mu", "mu_e", "beta_min"};
        }
        static std::vector<std::string> listOpts(){
            return { "method_synchrotron", "method_lf_min","method_nonreldist",
                     "method_tau", "use_ssa", "emissivity", "absorption" };
        }

//        double freq1 = -1;
//        double freq2 = -1;
//        size_t nfreq = 0;
    };
    Pars * p_pars;
    std::unique_ptr<logger> p_log;
    Array m_data{};
    Vector m_gamma_arr;
    Vector m_freq_arr_syn;
    Vector m_tmp_arr1{};
    Vector m_tmp_arr2{};
public:
    ///
    SynchrotronAnalytic( int loglevel ){

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronAnalytic");
        p_pars = new Pars();
        p_pars->lim_gm_to_1 = false;
        m_data.resize( m_names_.size(), -1. );
    }
    ~SynchrotronAnalytic(){ delete p_pars; }
    Pars *& getPars(){ return p_pars; }
    /// set model parameters
    static auto listPars(){ return Pars::listPars(); }
    static auto listOpts(){ return Pars::listOpts(); }
    void setPars(StrDbMap & pars, StrStrMap & opts){
        // set parameters
        p_pars->ksi_n = getDoublePar("ksi_n",pars,AT,p_log,1.,false);//pars.at("ksi_n");
        p_pars->eps_e = getDoublePar("eps_e",pars,AT,p_log,-1,true);//pars.at("eps_e");
        p_pars->eps_b = getDoublePar("eps_b",pars,AT,p_log,-1,true);//pars.at("eps_b");
        p_pars->eps_t = getDoublePar("eps_t",pars,AT,p_log,0.,true);//pars.at("eps_t");
        p_pars->p = getDoublePar("p",pars,AT,p_log,-1,true);//pars.at("p");
        p_pars->mu = getDoublePar("mu",pars,AT,p_log,0.62,false);//pars.at("mu");
        p_pars->mu_e = getDoublePar("mu_e",pars,AT,p_log,1.18,false);//pars.at("mu_e");
        p_pars->beta_min = getDoublePar("beta_min",pars,AT,p_log,1.e-5,false);//pars.at("beta_min");


        // set options
        std::string opt;
        opt = "method_synchrotron";
        METHODS val_synch;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_synch = SynchrotronAnalytic::METHODS::iJOH06;
        }
        else{
            if(opts.at(opt) == "Joh06")
                val_synch = SynchrotronAnalytic::METHODS::iJOH06;
            else if(opts.at(opt) == "WSPN99")
                val_synch = SynchrotronAnalytic::METHODS::iWSPN99;
            else if(opts.at(opt) == "Marg21")
                val_synch = SynchrotronAnalytic::METHODS::iMARG21;
            else if(opts.at(opt) == "Dermer09")
                val_synch = SynchrotronAnalytic::METHODS::iDER06;
//            else if(opts.at(opt) == "Bretta")
//                val_synch = SynchrotronAnalytic::METHODS::iBerrettaSynch;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " Joh06 " << " WSPN99 " << " Marg21 " << " Dermer09 " << "\n";
                exit(1);
            }
        }
        p_pars->m_sychMethod = val_synch;

        opt = "method_nonreldist";
        SynchrotronAnalytic::METHOD_NONRELDIST val_monreldist;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for " << opt << " is not set. Using default value.\n";
            val_monreldist = SynchrotronAnalytic::METHOD_NONRELDIST::inone;
        }
        else{
            if(opts.at(opt) == "none")
                val_monreldist = SynchrotronAnalytic::METHOD_NONRELDIST::inone;
            else if(opts.at(opt) == "useGm")
                val_monreldist = SynchrotronAnalytic::METHOD_NONRELDIST::iuseGm;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " none " << " useGm " << "\n";
                exit(1);
            }
        }
        p_pars->m_method_nonreldist = val_monreldist;

        opt = "method_lf_min";
        SynchrotronAnalytic::METHODS_LFMIN val_lfmin;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmUprime;
        }
        else{
            if(opts.at(opt) == "useU_e")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmUprime;
            else if(opts.at(opt) == "useBeta")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmNakarPiran;
            else if(opts.at(opt) == "useGamma")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmJoh06;
            else if(opts.at(opt) == "useTheta")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmMAG21;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " useU_e " << " useTheta "<< " useBeta " << " useGamma " << "\n";
                exit(1);
            }
        }
        p_pars->m_methodsLfmin = val_lfmin;

        bool tmp = getBoolOpt("use_ssa", opts, AT,p_log, false, true);
        if (tmp) p_pars->m_methods_ssa = SynchrotronAnalytic::METHODS_SSA::iSSAon;
        else p_pars->m_methods_ssa = SynchrotronAnalytic::METHODS_SSA::iSSAoff;

        opt = "emissivity";
        SynchrotronAnalytic::QQ val_em;
        if ( opts.find(opt) == opts.end() ) {
            val_em = SynchrotronAnalytic::QQ::i_em;
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
        }
        else{
            if(opts.at(opt) == "em_pl")
                val_em = SynchrotronAnalytic::QQ::i_em_pl;
            else if(opts.at(opt) == "em_th")
                val_em = SynchrotronAnalytic::QQ::i_em_th;
            else if(opts.at(opt) == "em")
                val_em = SynchrotronAnalytic::QQ::i_em;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " em_pl " << " em_th " << " em " << "\n";
                exit(1);
            }
        }
        p_pars->m_marg21opt_em = val_em;

        opt = "absorption";
        SynchrotronAnalytic::QQ val_abs;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_abs = SynchrotronAnalytic::QQ::i_abs;
        }
        else{
            if(opts.at(opt) == "abs_pl")
                val_abs = SynchrotronAnalytic::QQ::i_abs_pl;
            else if(opts.at(opt) == "abs_th")
                val_abs = SynchrotronAnalytic::QQ::i_abs_th;
            else if(opts.at(opt) == "abs")
                val_abs = SynchrotronAnalytic::QQ::i_abs;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " abs_pl " << " abs_th " << " abs " << "\n";
                exit(1);
            }
        }
        p_pars->m_marg21opt_abs = val_abs;
        // ---

        opt = "method_tau";
        METHOD_TAU methodTau;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodTau = METHOD_TAU::iSMOOTH;
        }
        else{
            if(opts.at(opt) == "thick")
                methodTau = METHOD_TAU::iTHICK;
            else if(opts.at(opt) == "sharp")
                methodTau = METHOD_TAU::iSHARP;
            else if(opts.at(opt) == "smooth")
                methodTau = METHOD_TAU::iSMOOTH;
            else if(opts.at(opt) == "approx")
                methodTau = METHOD_TAU::iAPPROX;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << " Possible options: "
                          << " approx " << " thick " << " sharp " << " smooth " << "\n";
                exit(1);
            }
        }
        p_pars->method_tau = methodTau;

        // ---
//        (*p_log)(LOG_INFO,AT) << " allocating memory for gamma_arr\n";

//        double m_gam1 = 1.;
//        double m_gam2 = 1e6;
//        int ngam = (int)getDoublePar("ngamma",pars,AT,p_log,-1.,true);
//        m_gamma_arr = TOOLS::MakeLogspaceVec(std::log10(m_gam1), std::log10(m_gam2), ngam, 10);
//        m_tmp_arr1.resize(m_gamma_arr.size(), 0.);
//        m_tmp_arr2.resize(m_gamma_arr.size(), 0.);


    }
    /// evaluate frequency independent quantities (critical LFs, Bfield, etc)
    void precompute(const double epime, const double Gamma, const double Gamma_shock, const double t_e, const double n_prime) {
        double beta = EQS::Beta(Gamma_shock);
        /// if velocity is too small shock may not be possible
        if (beta < p_pars->beta_min) {
            return;
        }

        double p = p_pars->p;
        double eps_b = p_pars->eps_b;
        double eps_e = p_pars->eps_e;
        double eps_t = p_pars->eps_t;
        double ksi_n = p_pars->ksi_n;
        // check
        if (eps_e <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<eps_e<<")\n";
            exit(1);
        }
        if (eps_b <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_b is not set (eps_e="<<eps_b<<")\n";
            exit(1);
        }
        if (eps_t < 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<eps_t<<") (even through it is only needed in Marg21 model)\n";
            exit(1);
        }
        if (p <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<p<<")\n";
            exit(1);
        }
        if (ksi_n <= 0){
            (*p_log)(LOG_ERR,AT)<< " ksi_n is not set (ksi_n="<<ksi_n<<")\n";
            exit(1);
        }
        if (n_prime <= 0){
            (*p_log)(LOG_ERR,AT) << " n_prime is not set (n_prime="<<n_prime<<")\n";
            exit(1);
        }

        // break if the vlocity is too low TODO update it with proper sound speed check

        // fraction of electrons that are accelerated
        double nprime = n_prime * ksi_n;

        /// fraction of the shock energy in mag. fields
        double U_b_prime = epime * eps_b;
        double B = sqrt(8. * CGS::pi * U_b_prime);

        /// fraction of the shock energy density in electrons
        double U_e_prime = epime * eps_e;

        /// get electron distribution boundaries
        double gM = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / B); // Kumar+14
        double gm = 0, accel_frac = 1.;
        switch (p_pars->m_methodsLfmin) {
            case igmUprime:
                gm = (p - 2.) / (p - 1.) * U_e_prime / (nprime * CGS::me * CGS::c * CGS::c); // Eq. A18 vanEarten+10 (and Sironi+13)
                break;
            case igmNakarPiran:
                gm = (p - 2.) / (p - 1.) * (CGS::mp / CGS::me) * eps_e * EQS::Beta2(Gamma_shock); // Eq.(before 4.) in Nakar&Piran 1102.1020
                break;
            case igmJoh06:
                gm = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (Gamma_shock - 1.) + 1.); // Eq. A3 in J+06
                break;
            case igmMAG21:
                double Theta = Margalit21::Theta_fun(beta, p_pars->mu, p_pars->mu_e, eps_t);
                gm = Margalit21::gamma_m_fun(Theta);

        }
        if ((p_pars->lim_gm_to_1) && (gm < 1.)) {
            gm = 1.; // Sironi et al 2013 suggestion to limi gm=1. for Deep Newtinoan regime # TODO to be removed. I did no understand it
        }

        switch (p_pars->m_method_nonreldist) {

            case inone:
                break;
            case iuseGm:
                accel_frac = (std::pow(gM, 2-p) - std::pow(gm-1, 2-p)) / (std::pow(gM,2-p) - 1.)
                             * (std::pow(gM,1-p)-1.) / (std::pow(gM,1-p) - std::pow(gm-1, 1-p));
                if (accel_frac > 1.) accel_frac = 1.;
                break;
        }

        double gc = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * t_e * B * B) / Gamma; // Eq. A19 in vanEarten+10
        /// in Margalit+21 model thermal electrons are present. Different treatment
        double Theta=-1., z_cool=-1.;
        if (p_pars->m_sychMethod == METHODS::iMARG21) {
            /* Margalit+21 arXiv:2111.00012 */
            /// downstream electron temperature:
            Theta = Margalit21::Theta_fun(beta, p_pars->mu, p_pars->mu_e, eps_t);
            /// minimum LFs for PL electrons is a average of the thermal ones
            gm = Margalit21::gamma_m_fun(Theta); // unused in the actual model. For analysis only.
            /// shock-amplified magnetic field (eq. 9; MQ21):
            B = sqrt(9.0 * M_PI * eps_b * nprime * p_pars->mu * CGS::mp) * (beta * CGS::c);
            /// mean dynamical time:
            double t = t_e;
            /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21):
            z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * B * B * t)) / Theta;
            //
            if ((Theta <= 0)||(!std::isfinite(z_cool))){
                (*p_log)(LOG_ERR,AT) << AT << " Theta = " << Theta << " z_cool="<<z_cool<<"\n";
                exit(1);
            }
        }

        /// check
        if (!std::isfinite(gm)){
            (*p_log)(LOG_ERR,AT) << " error gm nan \n";
            exit(1);
        }

        /// save the result
        p_pars->Theta = Theta;
        p_pars->z_cool = z_cool;
        p_pars->B = B;
        p_pars->gamma_min = gm;
        p_pars->gamma_max = gM;
        p_pars->gamma_c = gc;
        p_pars->n_prime = nprime;
        p_pars->accel_frac = accel_frac;

    }
    /// evaluate the comoving emissivity and absorption (frequency dependent)
    void compute(const double ndens_e, const double n_e, const double acc_frac,
                 const double B, double gm, const double gM, const double gc,
                 const double Theta, const double z_cool, const double nuprime ){

        // TODO WARNING I did replace n_prime with ne is absorption, but this might not be correct!!!

        double p = p_pars->p;
        double eps_e = p_pars->eps_e;
        double eps_t = p_pars->eps_t;
        double delta = eps_e / eps_t; // defined in Margalit+21 arXiv:2111.00012
//        double n_prime = nprime * p_pars->ksi_n;
        double ne = ndens_e * p_pars->ksi_n;

        /// compute synchrotron emissivity 'j' corrected for SSA, i.e., j * (1-e^-tau)/tau
        double x, em_th=0., em_pl=0., abs_th=0., abs_pl=0.; // for Margalit model
        if (p_pars->m_sychMethod == METHODS::iJOH06){

            double gamToNuFactor = (3.0 / (4.0 * CGS::pi)) * (CGS::qe * B) / (CGS::me * CGS::c);
            double XpF = 0.455 + 0.08 * p;
            double XpS = 0.06 + 0.28 * p;
            double phipF = 1.89 - 0.935 * p + 0.17 * (p * p);
            double phipS = 0.54 + 0.08 * p;
            double kappa1 = 2.37 - 0.3 * p;
            double kappa2 = 14.7 - 8.68 * p + 1.4 * (p * p);
            double kappa3 = 6.94 - 3.844 * p + 0.62 * (p * p);
            double kappa4 = 3.5 - 0.2 * p;
            double kappa13 = -kappa1 / 3.0;
            double kappa12 = kappa1 / 2.0;
            double kappa11 = -1.0 / kappa1;
            double kappa2p = kappa2 * (p - 1.0) / 2.0;
            double kappa12inv = -1.0 / kappa2;
            double kappa33 = -kappa3 / 3.;
            double kappa3p = kappa3 * (p - 1.0) / 2.0;
            double kappa13inv = -1.0 / kappa3;
            double kappa42 = kappa4 / 2.0;
            double kappa14 = -1.0 / kappa4;
            double nu_m, nu_c, emissivity, scaling, abs=-1., abs_scaling=1.;
            if (gm < gc){
                // slow cooling
                nu_m = XpS * gm * gm * gamToNuFactor;
                nu_c = XpS * gc * gc * gamToNuFactor;
//                double _phip = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * phipS;
                emissivity = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p)// phipS
                             * CGS::qe * CGS::qe * CGS::qe * ndens_e * B / (CGS::me * CGS::c * CGS::c);
                scaling = std::pow(std::pow(nuprime / nu_m, kappa33) + std::pow(nuprime / nu_m, kappa3p), kappa13inv)
                          * std::pow(1. + std::pow(nuprime / nu_c, kappa42), kappa14);
//                if (nuprime < nu_m) emissivity = 0.;
                /// -- SSA
                if (p_pars->m_methods_ssa!=iSSAoff) {
                    double _alpha = 7.8 * phipS * std::pow(XpS, -(4 + p) / 2.) * (p + 2) * (p - 1)
                                    * CGS::qe / CGS::mp / (p + 2 / 3.);
                    abs = _alpha * ne * CGS::mp * std::pow(gm, -5) / B;
                    if (nuprime <= nu_m)
                        abs_scaling = std::pow(nuprime / nu_m, -5 / 3.);
                    else if ((nu_m < nuprime) and (nuprime <= nu_c))
                        abs_scaling = std::pow(nuprime / nu_m, -(p + 4) / 2);
                    else if (nu_c < nuprime)
                        abs_scaling = std::pow(nu_c / nu_m, -(p + 4) / 2) * std::pow(nuprime / nu_c, -(p + 5) / 2);
                    else {
                        (*p_log)(LOG_ERR,AT) << "Error! in SSA\n";
                        exit(1);
                    }
                }
            }
            else {
                // fast cooling
                nu_m = XpF * gm * gm * gamToNuFactor;
                nu_c = XpF * gc * gc * gamToNuFactor;
                double _phip = 2.234 * phipF;
                emissivity = _phip * CGS::qe * CGS::qe * CGS::qe * ndens_e * B / (CGS::me * CGS::c * CGS::c);
                scaling = std::pow(std::pow(nuprime / nu_c, kappa13) + std::pow(nuprime / nu_c, kappa12), kappa11)
                          * std::pow(1. + std::pow(nuprime / nu_m, kappa2p), kappa12inv);
                /// --- SSA
                if (p_pars->m_methods_ssa!=iSSAoff) {
                    double _alpha = 11.7 * phipF * std::pow(XpF, -3) * CGS::qe / CGS::mp;
                    abs = _alpha * (ne * CGS::mp) * std::pow(gc, -5) / B;
                    if (nuprime <= nu_c)
                        abs_scaling = std::pow(nuprime / nu_c, -5 / 3.);
                    else if ((nu_c < nuprime) and (nuprime <= nu_m))
                        abs_scaling = std::pow(nuprime / nu_c, -3);
                    else if (nu_m < nuprime)
                        abs_scaling = std::pow(nu_m / nu_c, -3) * std::pow(nuprime / nu_m, -(p + 5) / 2);
                    else {
                        (*p_log)(LOG_ERR,AT) << "Error! in SSA\n";
                        exit(1);
                    }
                }
            }
            em_pl = emissivity * scaling;
            abs_pl = abs * abs_scaling;

//            std::cerr << AT << "\n";

        }
        else if (p_pars->m_sychMethod == METHODS::iWSPN99){



            double nu_m, nu_c, emissivity, scaling, abs=-1., abs_scaling=1.;
//            double PhiP = interpolated_phi(p);
//            emissivity = (interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi)
//                    * n_prime * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
            emissivity = (interpolated_phi(p) * sqrt(3.0) / 4.0 * CGS::pi)
                         * ndens_e * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
            double Xp = interpolated_xi(p);
            nu_m = 3.0 / ( 4.0 * CGS::pi ) * Xp * gm * gm * CGS::qe * B / ( CGS::me * CGS::c );
            nu_c = 0.286 * 3. * gc * gc * CGS::qe * B / ( 4.0 * CGS::pi * CGS::me * CGS::c );
            if (nu_m <= nu_c){//  # slow cooling
                if (nuprime < nu_m) {
                    scaling = std::pow(nuprime / nu_m, 1.0 / 3.0);
                }
                else if (nuprime >= nu_m && nuprime < nu_c) {
                    scaling = std::pow(nuprime / nu_m, -1.0 * (p - 1.0) / 2.0);
                }
                else  { // if (nuprime >= nu_c)
                    scaling = std::pow(nu_c / nu_m, -1.0 * (p - 1.0) / 2.0) * std::pow(nuprime / nu_c, -1.0 * p / 2.0);
                }
            }
            else {//  # fast cooling
                if (nuprime < nu_c){
                    scaling = std::pow(nuprime / nu_c, 1.0 / 3.0);
                }
                else if (nuprime >= nu_c && nuprime < nu_m) {
                    scaling = std::pow(nuprime / nu_c, -1.0 / 2.0);
                }
                else { // if (nuprime >= nu_m)
                    scaling = std::pow(nu_m / nu_c, -1.0 / 2.0) * std::pow(nuprime / nu_m, -p / 2.0);
                }
            }
            /// from vanEarten+2010
            if (p_pars->m_methods_ssa!=iSSAoff) {
                abs = sqrt(3) * std::pow(CGS::qe, 3) * (p - 1) * (p + 2) * ne * B
                      / (16 * M_PI * CGS::me * CGS::me * CGS::c * CGS::c * gm * nuprime * nuprime);
                if (nuprime < nu_m) // slow cooling
                    abs_scaling = std::pow(nuprime / nu_m, 1.0 / 3.0);
                else
                    abs_scaling = std::pow(nuprime / nu_m, -0.5 * p);
            }
            em_pl = emissivity * scaling;
            abs_pl = abs * abs_scaling;
//            std::cout << 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p) << " " << interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi << "\n";
//            exit(1);
        }
        else if (p_pars->m_sychMethod == METHODS::iDER06){
            double emissivity, abs=-1.;
            double epsilon = nuprime * CGS::h / CGS::mec2;
            // electron distribution, broken power law with two regimes depending on the order of gamma_i
            auto integrand_ele = [&](double gam){
                double p1 = gm < gc ? p : 2.0;
                double p2 = p + 1.0;
                double gmax = gM;
                double gmin = gm < gc ? gm : gc;
                double gb = gm < gc ? gc : gm;
                return Dermer09::brokenPowerLaw(gam, gmin, gb, gmax, p1, p2);
            };
            // compute electron distribution normalisation
            double k_e = ne / Simpson38(gm, gM, 200, integrand_ele); // TODO replace with adative integrals
            // convolve the electron distribution with emission spectrum and itegrate
            auto integrand = [&](double gam){
                double power_e = Dermer09::single_electron_synch_power( B, epsilon, gam );
                double n_e = k_e * integrand_ele(gam);
                return n_e * power_e;
            };
            double power = Simpson38(gm, gM, 200, integrand); // TODO replace with adative integrals
            emissivity = power * (CGS::h / CGS::mec2);
            /* --- SSA ---
             * Computes the syncrotron self-absorption opacity for a general set
             * of model parameters, see
             * :func:`~agnpy:sycnhrotron.Synchrotron.evaluate_sed_flux`
             * for parameters defintion.
             * Eq. before 7.122 in [DermerMenon2009]_.
             */
            if (p_pars->m_methods_ssa!=iSSAoff) {
                auto integrand_ele_ssa = [&](double gam) {
                    double p1 = gm < gc ? p : 2.0;
                    double p2 = p + 1.0;
                    double gmax = gM;
                    double gmin = gm < gc ? gm : gc;
                    double gb = gm < gc ? gc : gm;
                    return Dermer09::brokenPowerLawSSA(gam, gmin, gb, gmax, p1, p2);
                };
                // compute electron distribution normalisation
                double k_e_ssa = ne / Simpson38(gm, gM, 200, integrand_ele_ssa); // TODO replace with adative integrals
                // convolve the electron distribution with emission spectrum and itegrate
                auto integrand_ssa = [&](double gam) {
                    double power_e = Dermer09::single_electron_synch_power(B, epsilon, gam);
                    double n_e = k_e * integrand_ele_ssa(gam);
                    return n_e * power_e;
                };
                abs = Simpson38(gm, gM, 200, integrand); // TODO replace with adative integrals
                double coeff = -1 / (8 * CGS::pi * CGS::me * std::pow(epsilon, 2)) * std::pow(CGS::lambda_c / CGS::c, 3);
                abs *= coeff;
            }
            em_pl = emissivity;
            abs_pl = abs;
        }
        else if (p_pars->m_sychMethod == METHODS::iMARG21){
//            double emissivity, abs=-1.;
            /* Margalit+21 arXiv:2111.00012 */
            /// downstream (post-shock) electron number density:
            double ne_ = p_pars->mu_e * ne;//4.0 * p_pars->mu_e * p_pars->n_ism;
            /// normalized frequency:
            x = nuprime / Margalit21::nu_Theta(Theta, B);
            if (!std::isfinite(x)){
                (*p_log)(LOG_ERR,AT) << " x = "<< x << "\n";
                exit(1);
            }
            /// calculate total emissivity & optical depth:
            em_th = Margalit21::jnu_th(x, ne_, B, Theta, z_cool);
            em_pl = Margalit21::jnu_pl(x, ne_*acc_frac, B, Theta, gm, delta, p, z_cool); //TODO added 'accel_frac'
            if ((!std::isfinite(em_th))||(em_th < 0.)) em_th = 0.;
            if ((!std::isfinite(em_pl))||(em_pl < 0.)) em_pl = 0.;
//            emissivity = em_pl + em_th;
            if (p_pars->m_methods_ssa!=iSSAoff) {
                abs_th = Margalit21::alphanu_th(x, ne_, B, Theta, z_cool);
                abs_pl = Margalit21::alphanu_pl(x, ne_*acc_frac, B, Theta, gm, delta, p, z_cool);
                if (!std::isfinite(abs_th)) abs_th = 0.;
                if (!std::isfinite(abs_pl)) abs_pl = 0.;
//                abs = abs_th + abs_pl;
            }
            p_pars->x = x;
        }
        else if (p_pars->m_sychMethod == METHODS::iBerrettaSynch){
            /// fill gamma array
//            TOOLS::MakeLogspaceVec(m_gamma_arr, std::log10(gm), std::log10(gM));
            m_gamma_arr = TOOLS::MakeLogspaceVec( std::log10(gm), std::log10(gM), m_gamma_arr.size());
            std::cout << m_gamma_arr.size()<<" "<<" 1="<<m_gamma_arr[0]<<" -1="<<m_gamma_arr[m_gamma_arr.size()-1]<<"\n";
            double n0 = ne;
            const double nuL = Bretta::nu_L(B);
            const double n1 = 2.*CGS::pi*std::sqrt(3.)*(CGS::qe)*(CGS::qe)*(nuL)/(CGS::c); // Eq, ?? Dermer09?
            const double k0 = (p+2)/(8*CGS::pi*(CGS::me));
            double r=0., al=0.;
            for (size_t igam=0; igam < m_gamma_arr.size(); igam++){
                const double asso = Bretta::singleElecSynchroPower(nuprime, m_gamma_arr[igam], nuL)
                                    * std::pow(m_gamma_arr[igam], -(p+1)); // absorption
                const double js = Bretta::syncEmissivityKernPL(
                        m_gamma_arr[igam], nuprime, p, nuL, gc,m_gamma_arr); //
                r += js; // NOTE! replacing simpson integral with the direct sum!
                al += asso;  // NOTE! replacing simpson integral with the direct sum!
            }
            double I=0, att=-1;
            if (r > 1.e-90){
                I = n1*r*n0;
                att = n1*al*k0*n0/std::pow(nuprime,2.);
            }
//            em = I;
//            atten = att;

            em_pl = I;
            abs_pl = att;
        }
        else{
            (*p_log)(LOG_ERR,AT) << " synchrotoron method = " << p_pars->m_sychMethod << " is not recognized\n";
            exit(1);
        }

        m_data[i_em_pl] = em_pl; // 1e-26
        m_data[i_em_th] = em_th;
        m_data[i_em] = em_pl+em_th;

        if (p_pars->m_methods_ssa!=iSSAoff) {
            m_data[i_abs_th] = abs_th;
            m_data[i_abs_pl] = abs_pl; // 1e-17
            m_data[i_abs] = abs_th + abs_pl;
        }

        if (( em_pl < 0.) || (!std::isfinite( abs_pl )) ){
            (*p_log)(LOG_ERR,AT) << " em_pl_prime < 0 or nan ("<< em_pl<<") or \n";
            (*p_log)(LOG_ERR,AT) << " abs_pl_prime < 0 or nan ("<< abs_pl<<")\n";
            (*p_log)(LOG_ERR,AT) << " Error in data \n"
                      << " eps_e = " << eps_e << "\n"
                      << " eps_t = " << eps_t << "\n"
                      << " ne = " << ne << "\n"
                      << " gm = " << gm << "\n"
                      << " gM = " << gM << "\n"
                      << " gc = " << gc << "\n"
                      << " B = " << B << "\n"
                      << " Theta = " << Theta << "\n"
                      << " z_cool = " << z_cool << "\n"
                      << " nuprime = " << nuprime << "\n";
            exit(1);

            m_data[i_em_pl] = 0.;
            m_data[i_em_th] = 0.;
            m_data[i_em] = 0.;

            if (p_pars->m_methods_ssa!=iSSAoff) {
                m_data[i_abs_th] = -1;
                m_data[i_abs_pl] = -1.;
                m_data[i_abs] = -1.;
            }
        }




    }
    inline double get_em(){
        QQ opt = p_pars->m_marg21opt_em;
        double res = 0.;
        if (opt == QQ::i_em){
            return m_data[QQ::i_em];
        }
        if (opt == QQ::i_em_pl){
            return m_data[QQ::i_em_pl];
        }
        if (opt == QQ::i_em_th){
            return m_data[QQ::i_em_th];
        }
        (*p_log)(LOG_ERR,AT) << "should not be entered.\n";
        exit(1);
    }
    inline double get_abs(){
        QQ opt = p_pars->m_marg21opt_abs;
        double res = 0.;
        if (opt == QQ::i_abs){
            return m_data[QQ::i_abs];
        }
        if (opt == QQ::i_abs_pl){
            return m_data[QQ::i_abs_pl];
        }
        if (opt == QQ::i_abs_th){
            return m_data[QQ::i_abs_th];
        }
        (*p_log)(LOG_ERR,AT) << "should not be entered.\n";
        exit(1);
    }
    void addOptDepth(double abs, double dr, double beta_sh, double mu){
        m_data[QQ::i_tau] = RadiationBase::optical_depth(abs,dr,mu,beta_sh);
    }
    void addIntensity(double dr, double beta_sh, double mu){
//        std::cout << m_data[i_abs] << "\n";
        m_data[QQ::i_tau] = RadiationBase::optical_depth(m_data[i_abs], dr, mu, beta_sh);
        m_data[QQ::i_int] = RadiationBase::computeIntensity(m_data[i_em], m_data[QQ::i_tau], p_pars->method_tau);
    }
    inline Array & getData(){return m_data;}
    inline double getIntensity(){ return m_data[QQ::i_int]; }
//    inline Array & getSpectrum(){ return m_spectrum; }
//    inline Array & getFreqArr(){ return m_freq_grid; }
    std::vector<std::string> m_names_ { "em_pl", "em_th", "abs_pl", "abs_th", "em", "abs", "tau", "int" };
};


//    void allocateSpaceForComovingSpectrum(size_t nr){
//        /// allocate memory for comoving spectra to be evaluated
//        m_freq_grid = TOOLS::MakeLogspace(log10(p_pars->freq1),
//                                          log10(p_pars->freq2),(int)p_pars->nfreq);
//        if (p_pars->method_comp_mode == SynchrotronAnalytic::METHODS_RAD::icomovspec){
//            std::cerr << " allocating comoving spectrum array (fs) "
//                               << " freqs="<<m_freq_grid.size() << " by radii=" << nr << " Spec. grid="
//                               << m_freq_grid.size() * nr << "\n";
//            m_spectrum.resize( m_freq_grid.size() * nr );
//        }
//    }

//class SynchrotronNumeric{
//    std::unique_ptr<logger> p_log;
//    Vector farr_sync;
//    Vector freq_ssc;
//    struct Pars{
//        double f1_sync, f2_sync; size_t nf_sync;
//        double f1_ssc, f2_ssc; size_t nf_ssc;
//
//    };
//    std::unique_ptr<Pars> p_pars;
//public:
//    SynchrotronNumeric(int loglevel){
//        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronNumeric");
//        p_pars = std::make_unique<Pars>();
//    }
//    void setPars(StrDbMap pars, StrStrMap sets){
////        p_pars->f1_ssc =
//    }
//};


//class SynchrotronAnalyticComoving{
//    enum METHODS { iWSPN99, iJOH06, iDER06, iMARG21 };
//    enum METHODS_LFMIN { igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };
//    enum METHODS_SSA { iSSAoff, iSSAon };
//    enum METHOD_TAU { iTHICK, iSMOOTH, iSHARP };
//    enum QQ { i_em_pl, i_em_th, i_abs_pl, i_abs_th, i_em, i_abs };
//    struct Pars{
//        // --- in
//        double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
//        double mu=-1, mu_e=-1;
//        bool lim_gm_to_1= true;
//        double beta_min = -1;
//        // --- methods
//        METHODS m_sychMethod{};
//        METHODS_LFMIN m_methodsLfmin{};
//
//        METHODS_SSA m_methods_ssa{};
//        QQ m_marg21opt_em = i_em_pl;
//        QQ m_marg21opt_abs = i_abs_pl;
//    };
//    Pars * p_pars = nullptr;
//
//    static constexpr size_t n_vars = 5;
//    enum QS { igm, igc };
//    std::vector<std::string> m_vars{ "gm", "gc", "num", "nuc", "pmax" };
//    VecArray m_data{};// ( n_vars );
//public:
//    SynchrotronAnalyticComoving(size_t nt){
//        p_pars = new Pars();
//        allocateSpace(nt);
//    }
//    ~SynchrotronAnalyticComoving(){ delete p_pars; }
//    void allocateSpace(size_t nt){
//        m_data.resize(n_vars);
//        for (auto & arr : m_data){
//            arr.resize( nt, 0. );
//        }
//    }
//
//    void setPars(std::unordered_map<std::string,double> & pars, std::unordered_map<std::string,std::string> & opts){
//        // set parameters
//
//    }
//};

#endif //SRC_SYNCHROTRON_AN_H
