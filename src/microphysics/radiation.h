//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_RADIATION_H
#define SRC_RADIATION_H

#include "../utilitites/pch.h"
#include "../utilitites/logger.h"
#include "../utilitites/utils.h"
#include "../utilitites/interpolators.h"
#include "../utilitites/quadratures.h"
#include "../utilitites/rootfinders.h"

#include "../blastwave/blastwave_components.h"

#include "numeric_model.h"

enum METHODS_SHOCK_ELE { iShockEleAnalyt, iShockEleNum };

enum METHOD_TAU { iAPPROX, iTHICK, iSMOOTH, iSHARP };

enum METHODS_SYNCH { iWSPN99, iJOH06, iDER06, iMARG21, iNumeric, iGSL };

enum METHOD_SSC { inoSSC, iNumSSC };

enum METHOD_ELE { iEleAnalytic, iEleNumeric };

enum METHODS_LFMIN { igmNumGamma, igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };

enum METHODS_B { iBuseUb, iBasMAG21, iBuseGammaSh };

enum METHOD_LFMAX { iConst, iuseB };

enum METHODS_SSA { iSSAoff, iSSAon };

enum METHOD_NONRELDIST{ inone, iuseGm };


/// Evaluate optical depth
inline double optical_depth(const double abs_lab, const double dr,
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

/// Compute the I = emission/absorption * ( 1 - tau^{-absorption * thickness})
double computeIntensity(const double em_lab, double dtau, METHOD_TAU m_tau){

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

/// From https://github.com/bmargalit/thermal-synchrotron
struct Margalit21 {

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
    a_fun(const double Theta) {
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
    static double
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

        /// TODO add assymptotic expansion for small Theta
        if (Theta < 3.e-7)
            return 0.;

        double tmp = 0;
        try {
            tmp = 2.0 * Theta * Theta / std::cyl_bessel_k(2, 1.0 / Theta);
        } catch (const std::runtime_error& error) {
            std::cout << AT << " " << error.what() << " Theta="<<Theta<< "\n";
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
                           * std::pow(gamma_m / (3.0 * Theta), (p - 1.0));
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
    static double
    I_of_x(const double x, const double Theta) {

        /// Mahadevan+95 "Harmony in Electrons"
        const size_t n = 7;
        const Vector temps = { 5e8, 1e9, 2e9, 4e9, 8e9, 1.6e10, 3.2e10 };
        const Vector alphas = { 0.0431, 1.121, 1.180, 1.045, 0.9774, 0.9768, 0.9788 };
        const Vector betas = {10.44, -10.65, -4.008, -0.1897, 1.160, 1.095, 1.021};
        const Vector gammas = {16.61, 9.169, 1.559, 0.0595, 0.2641, 0.8332, 1.031};

        Vector Thetas (temps.size());
        for (size_t i = 0; i < temps.size(); i++)
            Thetas[i] = temps[i] * CGS::kB / (CGS::me * CGS::c * CGS::c);

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
                         * std::exp(-1.8899 * std::pow(x, 1.0 / 3.0));
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
    static double
    C_j(const double p) {
        double tmp = 0;
        try {
            tmp = (tgamma((p + 5.0) / 4.0) / tgamma((p + 7.0) / 4.0)) *
                  tgamma((3.0 * p + 19.0) / 12.0) * tgamma((3.0 * p - 1.0) / 12.0) *
                  ((p - 2.0) / (p + 1.0)) * std::pow(3.0, (2.0 * p - 1.0) / 2.0)
                  * std::pow(2.0, -(7.0 - p) / 2.0) * std::pow(M_PI, -0.5);
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << " p = "<<p<<"\n";
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
    static double
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
    static double
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
    static double
    low_freq_apl_correction(const double x, const double Theta, const double p) {

        double gamma_m = gamma_m_fun(Theta);
        /// synchrotron constant in x<<x_m limit
        double Calpha_low = 0;
        try {
            Calpha_low = -std::pow(2.0, 8.0 / 3.0) * std::pow(M_PI, 7.0 / 2.0) * (p + 2.0) * (p - 2.0) /
                         (std::pow(3.0, 19.0 / 6.0) * (3.0 * p + 2) * tgamma(1.0 / 3.0) * tgamma(-1.0 / 3.0) *
                          tgamma(11.0 / 6.0));
        } catch (const std::runtime_error& error) {
            std::cout << AT << error.what() << " x="<<x<<" Theta="<<Theta<<" p="<<p<< "\n";
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
    static double
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
    static double
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
    static double
    jnu_th(const double x, const double n, const double B, const double Theta, const double z_cool){

        double val = (sqrt(3.0) / (8.0 * M_PI))
                    * (CGS::qe*CGS::qe*CGS::qe / (CGS::me * CGS::c * CGS::c))
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
double interpolated_xi(const double p){

    const Vector p_xi = {
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

    const Vector Xi = {
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

    return ( Interp1d( p_xi, Xi) ).Interpolate(p, Interp1d::iLagrangeLinear);
    //return LinInterp(m_Xp.p_xi, m_Xp.Xi, p_xi,  false);
}
double interpolated_phi(const double p){

    const Vector p_phi = {
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
    const Vector phi = {
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
        if ((gm <= gamma_e) && (gamma_e <= gM)){
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
        double x = 4.0 * CGS::pi * epsilon * (CGS::me * CGS::me) * (CGS::c * CGS::c * CGS::c) / (3.0 * CGS::qe * B * CGS::h * (gamma * gamma));
        return x;
    }
    /**
     * Eq. 7.45 in [Dermer2009]_, angle-averaged integrand of the radiated power, the
     * approximation of this function, given in Eq. D7 of [Aharonian2010]_, is used.
     *
     * @param x = ratio of the frequency to the critical synchrotron frequency
     * @return
     */
    static double R(const double x){

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
    static double single_electron_synch_power(const double B, const double epsilon, const double gamma){
        double result = 0.;
        double prefactor = sqrt(3.0) * (CGS::qe * CGS::qe * CGS::qe) * B / CGS::h;
        double x = calc_x(B, epsilon, gamma);
        result = prefactor * R(x);
        return result;
    }

    double pprime(double nuprim, double p, double gm, double gc, double gM, double B, double nprim){

//        double gM = 1e10;//Electrons::gamma_max(B);
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

        // computeSynchrotronEmissivityAbsorptionAnalytic electron distribution normalisation
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


class ElectronAndRadiaionBase{
public:
    Vector out_spectrum{}; Vector out_specturm_ssa{};
    Vector m_freq_arr{};
    /// ------------------------------------
    double B=-1, gamma_min=-1, gamma_max=-1, gamma_c=-1;
    double n_prime=-1, eprime=-1.,Gamma=-1,GammaSh=-1, beta=-1.,t_e=-1.,n_protons=-1.;
    double accel_frac=-1.;
    double Theta=-1, z_cool=-1, x=-1;
    double em=-1., abs=-1.;
    METHOD_TAU method_tau{}; double beta_min = -1;
    METHODS_SHOCK_ELE m_eleMethod{};
    /// --------------------------------------
    void setBasePars( StrDbMap & pars, StrStrMap & opts ){

        // set parameters
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";

        ksi_n = getDoublePar("ksi_n" + fs_or_rs, pars, AT, p_log, 1., false);//pars.at("ksi_n");
        eps_e = getDoublePar("eps_e" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_e");
        eps_b = getDoublePar("eps_b" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_b");
        eps_t = getDoublePar("eps_t" + fs_or_rs, pars, AT, p_log, 0., true);//pars.at("eps_t");
        p = getDoublePar("p" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("p");
        mu = getDoublePar("mu" + fs_or_rs, pars, AT, p_log, 0.62, false);//pars.at("mu");
        mu_e = getDoublePar("mu_e" + fs_or_rs, pars, AT, p_log, 1.18, false);//pars.at("mu_e");
        beta_min = getDoublePar("beta_min" + fs_or_rs, pars, AT, p_log, 1.e-5, false);//pars.at("beta_min");
        gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars, AT, p_log, 1.e7, false);//pars.at("beta_min");

        lim_gm_to_1 = getBoolOpt("limit_lf_min_to1" + fs_or_rs, opts, AT, p_log, false, false);//pars.at("beta_min");

        // set options

        std::string opt;
        opt = "method_shock_ele" + fs_or_rs;
        METHODS_SHOCK_ELE val_ele;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_ele = METHODS_SHOCK_ELE::iShockEleAnalyt;
        }
        else{
            if(opts.at(opt) == "analytic")
                val_ele = METHODS_SHOCK_ELE::iShockEleAnalyt;
            else if(opts.at(opt) == "numeric")
                val_ele = METHODS_SHOCK_ELE::iShockEleNum;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " analytic " << " numeric " << "\n";
                exit(1);
            }
        }
        m_eleMethod = val_ele;

        opt = "method_synchrotron" + fs_or_rs;
        METHODS_SYNCH val_synch;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_synch = METHODS_SYNCH::iJOH06;
        }
        else{
            if(opts.at(opt) == "Joh06")
                val_synch = METHODS_SYNCH::iJOH06;
            else if(opts.at(opt) == "WSPN99")
                val_synch = METHODS_SYNCH::iWSPN99;
            else if(opts.at(opt) == "Marg21")
                val_synch = METHODS_SYNCH::iMARG21;
            else if(opts.at(opt) == "Dermer09")
                val_synch = METHODS_SYNCH::iDER06;
            else if(opts.at(opt) == "GSL") {
                if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT) << " options synchrotron GSL and analytic electron evol. are incompatible\n";
                    exit(1);
                }
                val_synch = METHODS_SYNCH::iGSL;

            }
//            else if(opts.at(opt) == "Bretta")
//                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iBerrettaSynch;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " Joh06 " << " WSPN99 " << " Marg21 " << " Dermer09 " << " GSL " << "\n";
                exit(1);
            }
        }
        m_sychMethod = val_synch;

        opt = "method_ssc" + fs_or_rs;
        METHOD_SSC val_ssc;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_ssc = METHOD_SSC::inoSSC;
        }
        else{
            if(opts.at(opt) == "none")
                val_ssc = METHOD_SSC::inoSSC;
            else if(opts.at(opt) == "numeric") {
                if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT) << " SSC is not supported for analytic electron dsitrib. model\n";
                    exit(1);
                }
                val_ssc = METHOD_SSC::iNumSSC;
            }
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " none " << " numeric " << "\n";
                exit(1);
            }
        }
        m_methods_ssc = val_ssc;

        opt = "method_nonreldist" + fs_or_rs;
        METHOD_NONRELDIST val_monreldist;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for " << opt << " is not set. Using default value.\n";
            val_monreldist = METHOD_NONRELDIST::inone;
        }
        else{
            if(opts.at(opt) == "none")
                val_monreldist = METHOD_NONRELDIST::inone;
            else if(opts.at(opt) == "useGm")
                val_monreldist = METHOD_NONRELDIST::iuseGm;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " none " << " useGm " << "\n";
                exit(1);
            }
        }
        m_method_nonreldist = val_monreldist;

        opt = "method_lf_min" + fs_or_rs;
        METHODS_LFMIN val_lfmin;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_lfmin = METHODS_LFMIN::igmUprime;
        }
        else{
            if(opts.at(opt) == "useU_e")
                val_lfmin = METHODS_LFMIN::igmUprime;
            else if(opts.at(opt) == "useNumericGamma")
                val_lfmin = METHODS_LFMIN::igmNumGamma;
            else if(opts.at(opt) == "useBeta")
                val_lfmin = METHODS_LFMIN::igmNakarPiran;
            else if(opts.at(opt) == "useGamma")
                val_lfmin = METHODS_LFMIN::igmJoh06;
            else if(opts.at(opt) == "useTheta")
                val_lfmin = METHODS_LFMIN::igmMAG21;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useU_e " << " useTheta "<< " useBeta " << " useGamma " << "\n";
                exit(1);
            }
        }
        m_methodsLfmin = val_lfmin;

        opt = "method_Bsh" + fs_or_rs;
        METHODS_B methodsB;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsB = METHODS_B::iBuseUb;
        }
        else{
            if(opts.at(opt) == "useUb")
                methodsB = METHODS_B::iBuseUb;
            else if(opts.at(opt) == "useMAG21")
                methodsB = METHODS_B::iBasMAG21;
            else if(opts.at(opt) == "useGammaSh")
                methodsB = METHODS_B::iBuseGammaSh;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useUb " << " useMAG21 " << " useGammaSh " << "\n";
                exit(1);
            }
        }
        m_methodsB = methodsB;

        opt = "method_lf_max" + fs_or_rs;
        METHOD_LFMAX methodLfmax;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodLfmax = METHOD_LFMAX::iuseB;
        }
        else{
            if(opts.at(opt) == "useB")
                methodLfmax = METHOD_LFMAX::iuseB;
            else if(opts.at(opt) == "useConst") {
                gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars,
                                         AT, p_log, 1.e7, true);//pars.at("beta_min");
                methodLfmax = METHOD_LFMAX::iConst;
            }
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useB " << " useConst " << "\n";
                exit(1);
            }
        }
        m_methodsLfmax = methodLfmax;

        bool tmp = getBoolOpt("use_ssa" + fs_or_rs, opts, AT, p_log, false, true);
        if (tmp) m_methods_ssa = METHODS_SSA::iSSAon;
        else m_methods_ssa = METHODS_SSA::iSSAoff;

#if 0
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
#endif
        opt = "method_tau" + fs_or_rs;
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
        method_tau = methodTau;
    }
    void allocateStorageForOutputSpectra(size_t nr, double freq1, double freq2,  size_t nfreq){
        m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);
        /// allocate space for the comoving spectrum
        out_spectrum.resize(nfreq * nr );
        out_specturm_ssa.resize(nfreq * nr );

    }
    /// store current shock properties
    void updateSockProperties(double e_prime, double Gamma_, double Gamma_shock, double t_e_,
                              double n_prime_, double n_protons_){
        /// Store current parameters

        GammaSh = Gamma_shock; // for bisect solver
        eprime = e_prime; // for bisect solver
        Gamma = Gamma_; // for bisect solver
        beta = EQS::Beta(Gamma_shock); // for bisect solver
        t_e = t_e_; // for bisect solver
        n_prime = ksi_n * n_prime_; //  for bisect solver
        n_protons = ksi_n * n_protons_;
    }
    /// In case the electrons were computed elsewhere E.g., if interpolated for EATS plane
    void setShockElectronParameters(double n_prime_, double acc_frac,
                                    double B_, double gm, double gM, double gc,
                                    double Theta_, double z_cool_){

        if (!std::isfinite(gm) || !std::isfinite(gc) || !std::isfinite(n_prime_)) {
            (*p_log)(LOG_ERR, AT) << " nans is synchrotron spectum\n";
            exit(1);
        }

        gamma_min = gm;
        gamma_max = gM;
        n_prime = n_prime_;
        B = B_;
        accel_frac = acc_frac;
        gamma_c = gc;
        z_cool = z_cool_;
        Theta = Theta_;
    }
    /// evaluate frequency independent quantities (critical LFs, Bfield, etc)
    void evaluateElectronDistributionAnalytic() {

        if (not checkParams())
            return;

        computeMagneticField();

        computeGammaMax();

        computeGammaMin();

        computeNonRelativisticFlattening();

        /// compute cooling lorentz factor
        gamma_c = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * t_e * B * B) / Gamma; // Eq. A19 in vanEarten+10
    }
    /// ---
    [[nodiscard]] ElectronAndRadiaionBase const * getThis() const { return this; }
protected:
    int m_loglevel = -1;
    std::unique_ptr<logger> p_log = nullptr;
    bool is_rs = false;
    /// --------------------------------------
    double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
    double mu=-1, mu_e=-1;
    bool lim_gm_to_1= true;
    /// --------------------------------------
    METHODS_SYNCH m_sychMethod{};
    METHODS_LFMIN m_methodsLfmin{};
    METHODS_B m_methodsB{};
    METHOD_LFMAX m_methodsLfmax{};
    METHOD_NONRELDIST m_method_nonreldist{};
    METHODS_SSA m_methods_ssa{};
    METHOD_SSC m_methods_ssc{};
    /// --------------------------------------
    ElectronAndRadiaionBase(int loglevel, bool _is_rs){
        m_loglevel = loglevel; is_rs = _is_rs;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronAnalytic");
    }
    /// --------------------------------------
    bool checkParams(){

        /// if velocity is too small shock may not be possible TODO update it with proper sound speed check
        if (beta < beta_min)
            return false;


        /// in Margalit+21 model thermal electrons are present. Different treatment :: Assert:
        if (m_sychMethod == METHODS_SYNCH::iMARG21) {
            /* Margalit+21 arXiv:2111.00012 */
            if (m_methodsLfmin != igmMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsLfmin=gmMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
            if (m_methodsB != iBasMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsB=asMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
        }

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
        if ((eprime < 0) || !std::isfinite(eprime)){
            (*p_log)(LOG_ERR,AT) << " eprime is not set (eprime=" << eprime << ")\n";
            exit(1);
        }

        return true;
    }
    /// --------------------------------------
    static double gammaMinFunc(const double &x, void * pars){
        auto * pp = (struct ElectronAndRadiaionBase *) pars;
        return (pp->p - 1) / (pp->p - 2)
            * (std::pow(pp->gamma_max, -pp->p + 2) - std::pow(x, -pp->p + 2))
            / (std::pow(pp->gamma_max, -pp->p + 1) - std::pow(x, -pp->p + 1))
            - pp->eps_e * mp / me * (pp->GammaSh - 1);
    }
    /// --------------------------------------
    void computeMagneticField(){
        /// fraction of the shock energy density in electrons
        double U_b_prime;
        switch (m_methodsB) {
            case iBuseUb:
                U_b_prime = eprime * eps_b;
                B = sqrt(8. * CGS::pi * U_b_prime);
                break;
            case iBuseGammaSh:
                /// See Nava 2013 paper or others; classical equation
                B = std::sqrt( 32. * M_PI * eps_b * (GammaSh - 1.)
                  * (GammaSh + 3. / 4) * n_prime * CGS::mp * CGS::c2);
                break;
            case iBasMAG21:
                B = sqrt(9.0 * M_PI * eps_b * n_prime * mu * CGS::mp)
                  * (beta * CGS::c);
                break;
        }
    }
    void computeGammaMax(){
        /// get electron distribution boundaries (has ot be set BEFORE gamma_min)
        switch (m_methodsLfmax) {
            case iConst:
                break;
            case iuseB:
                gamma_max = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / B); // Kumar+14
                break;
        }
    }
    void computeGammaMin(){
        /// compute injection LF (lower bound of source injection gammaMinFunc)
        int status = 0; double U_e_prime = -1;
//        void * pars = (struct ElectronAndRadiaionBase *) this; // this "this" i do not understand myself...
        switch (m_methodsLfmin) {
            case igmUprime:
                U_e_prime = eprime * eps_e;
                gamma_min = (p - 2.) / (p - 1.) * U_e_prime / (n_prime * CGS::me * CGS::c * CGS::c); // Eq. A18 vanEarten+10 (and Sironi+13)
                break;
            case igmNakarPiran:
                gamma_min = (p - 2.) / (p - 1.) * (CGS::mp / CGS::me) * eps_e * EQS::Beta2(GammaSh); // Eq.(before 4.) in Nakar&Piran 1102.1020
                break;
            case igmJoh06:
                gamma_min = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (GammaSh - 1.) + 1.); // Eq. A3 in J+06
                break;
            case igmNumGamma:
                /// solve gamma_min fun numerically; use fixed limits and number of iterations
                gamma_min = Bisect(ElectronAndRadiaionBase::gammaMinFunc,
                                   1, 1e8, 0, .001, 100, this, status);
                /// If numerical solution failed, use simple analytical solution
                if (status < 0)
                    gamma_min = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (GammaSh - 1.) + 1.);
                break;
            case igmMAG21:
                /// downstream electron temperature:
                Theta = Margalit21::Theta_fun(beta, mu, mu_e, eps_t);
                gamma_min = Margalit21::gamma_m_fun(Theta);
                break;
        }

        /// check
        if (!std::isfinite(gamma_min)){
            (*p_log)(LOG_ERR,AT) << " error gm nan \n";
            exit(1);
        }

        /// limit the min lf to 1
        if ((lim_gm_to_1) && (gamma_min < 1.))
            gamma_min = 1.; // Sironi et al 2013 suggestion to limi gm=1. for Deep Newtinoan regime # TODO to be removed. I did no understand it

        /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21): NOTE we use mean dynamical time:
        z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * B * B * t_e)) / Theta;

        if (!std::isfinite(z_cool)){
            (*p_log)(LOG_ERR,AT) << AT << " Theta = " << Theta << " z_cool="<<z_cool<<"\n";
            exit(1);
        }
    }
    /// for semi-neutonian regime, where gm ->
    void computeNonRelativisticFlattening(){
        switch (m_method_nonreldist) {
            case inone:
                break;
            case iuseGm:
                if (lim_gm_to_1) {
                    (*p_log)(LOG_ERR, AT) << " 'method_nonreldist' is incopatible with 'lim_gm_to_1' \n";
                    exit(1);
                }
                accel_frac = (std::pow(gamma_max, 2-p) - std::pow(gamma_min-1, 2-p))
                           / (std::pow(gamma_max,2-p) - 1.)
                           * (std::pow(gamma_max,1-p)-1.) / (std::pow(gamma_max,1-p)
                           - std::pow(gamma_min-1, 1-p));
                if (accel_frac > 1.)
                    accel_frac = 1.;
                break;
        }
    }

};


class ElectronAndRadiation : public ElectronAndRadiaionBase{

    /// --------------------------------------
    Source source{};
    State ele{}; // std::unique_ptr<State> ele = nullptr;
    State syn{}; //std::unique_ptr<State> syn = nullptr;
    State ssc{}; //std::unique_ptr<State> ssc = nullptr;
    SynKernel syn_kernel{};//std::unique_ptr<SSCKernel> ssc_kernel = nullptr;
    SSCKernel ssc_kernel{};//std::unique_ptr<SynKernel> syn_kernel = nullptr;
    double vol=-1.;
    double vol_p1=-1.;
    double dm=-1.;
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchJOH06(double nuprime){
//        double gamma_min=p_pars->gamma_min, gamma_c=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
//        double p = p_pars->p;
//        double em=0., abs=0.;
        // prefactors
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
        double nu_m, nu_c, emissivity, scaling, abs_scaling=1.;

        double ne = n_prime * ksi_n;

        if (gamma_min < gamma_c){
            // slow cooling
            nu_m = XpS * gamma_min * gamma_min * gamToNuFactor;
            nu_c = XpS * gamma_c * gamma_c * gamToNuFactor;
//                double _phip = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * phipS;
            emissivity = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p)// phipS
                         * CGS::qe * CGS::qe * CGS::qe * n_prime * B / (CGS::me * CGS::c * CGS::c);
            scaling = std::pow(std::pow(nuprime / nu_m, kappa33) + std::pow(nuprime / nu_m, kappa3p), kappa13inv)
                      * std::pow(1. + std::pow(nuprime / nu_c, kappa42), kappa14);
//                if (nuprime < nu_m) emissivity = 0.;
            /// -- SSA
            if (m_methods_ssa!=iSSAoff) {
                double _alpha = 7.8 * phipS * std::pow(XpS, -(4 + p) / 2.) * (p + 2) * (p - 1)
                                * CGS::qe / CGS::mp / (p + 2 / 3.);
                abs = _alpha * ne * CGS::mp * std::pow(gamma_min, -5) / B;
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
            nu_m = XpF * gamma_min * gamma_min * gamToNuFactor;
            nu_c = XpF * gamma_c * gamma_c * gamToNuFactor;
            double _phip = 2.234 * phipF;
            emissivity = _phip * CGS::qe * CGS::qe * CGS::qe * n_prime * B / (CGS::me * CGS::c * CGS::c);
            scaling = std::pow(std::pow(nuprime / nu_c, kappa13) + std::pow(nuprime / nu_c, kappa12), kappa11)
                      * std::pow(1. + std::pow(nuprime / nu_m, kappa2p), kappa12inv);
            /// --- SSA
            if (m_methods_ssa!=iSSAoff) {
                double _alpha = 11.7 * phipF * std::pow(XpF, -3) * CGS::qe / CGS::mp;
                abs = _alpha * (ne * CGS::mp) * std::pow(gamma_c, -5) / B;
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

        em = emissivity * scaling;
        abs = abs * abs_scaling;
    };
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchWSPN99(double nuprime){
//        double gamma_min=p_pars->gamma_min, gamma_c=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
//        double p = p_pars->p;
        double nu_m, nu_c, emissivity, scaling, abs_scaling=1.;
//        double em=0.,abs=0.;
//            double PhiP = interpolated_phi(p);
//            emissivity = (interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi)
//                    * n_prime * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
        double ne = n_prime * ksi_n;
        emissivity = (interpolated_phi(p) * sqrt(3.0) / 4.0 * CGS::pi) * ne * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
        double Xp = interpolated_xi(p);
        nu_m = 3.0 / ( 4.0 * CGS::pi ) * Xp * gamma_min * gamma_min * CGS::qe * B / (CGS::me * CGS::c );
        nu_c = 0.286 * 3. * gamma_c * gamma_c * CGS::qe * B / (4.0 * CGS::pi * CGS::me * CGS::c );
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
        if (m_methods_ssa!=iSSAoff) {
            abs = sqrt(3) * std::pow(CGS::qe, 3) * (p - 1) * (p + 2) * ne * B
                  / (16 * M_PI * CGS::me * CGS::me * CGS::c * CGS::c * gamma_min * nuprime * nuprime);
            if (nuprime < nu_m) // slow cooling
                abs_scaling = std::pow(nuprime / nu_m, 1.0 / 3.0);
            else
                abs_scaling = std::pow(nuprime / nu_m, -0.5 * p);
        }
        em = emissivity * scaling;
        abs = abs * abs_scaling;
//            std::cout << 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p) << " " << interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi << "\n";
//            exit(1);
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchDER06(double nuprime){

        Vector gammmas = TOOLS::MakeLogspaceVec(std::log10(gamma_min),
                                                std::log10(gamma_max),
                                                200);

        double ne = n_prime * ksi_n;
        double epsilon = nuprime * CGS::h / CGS::mec2;
        auto integrand_ele = [&](const double gam){
            double p1 = gamma_min < gamma_c ? p : 2.0;
            double p2 = p + 1.0;
            double gmax = gamma_max;
            double gmin = gamma_min < gamma_c ? gamma_min : gamma_c;
            double gb = gamma_min < gamma_c ? gamma_c : gamma_min;
            return Dermer09::brokenPowerLaw(gam, gmin, gb, gmax, p1, p2);
        };
        Vector k_e_s (gammmas.size(), 0.);
        double k_e = 0., power_e = 0., power = 0.;
        for (size_t i = 0; i < gammmas.size()-1; i++) {
            k_e_s[i] = integrand_ele(gammmas[i]);
            k_e += k_e_s[i] * (gammmas[i + 1] - gammmas[i]);
        }
        k_e = ne / k_e;

        for (size_t i = 0; i < gammmas.size()-1; i++) {
            power_e = Dermer09::single_electron_synch_power( B, epsilon, gammmas[i] );
            power += k_e_s[i] * power_e * (gammmas[i + 1] - gammmas[i]);
        }
        power *= k_e;

        em = power * (CGS::h / CGS::mec2);

        if (m_methods_ssa!=iSSAoff) {
            auto integrand_ele_ssa = [&](double gam) {
                double p1 = gamma_min < gamma_c ? p : 2.0;
                double p2 = p + 1.0;
                double gmax = gamma_max;
                double gmin = gamma_min < gamma_c ? gamma_min : gamma_c;
                double gb = gamma_min < gamma_c ? gamma_c : gamma_min;
                return Dermer09::brokenPowerLawSSA(gam, gmin, gb, gmax, p1, p2);
            };

            Vector k_e_s_ssa (gammmas.size(), 0.);
            double k_e_ssa = 0., absorption = 0.;
            for (size_t i = 0; i < gammmas.size()-1; i++) {
                k_e_s_ssa[i] = integrand_ele_ssa(gammmas[i]);
                k_e_ssa += k_e_s_ssa[i] * (gammmas[i + 1] - gammmas[i]);
            }
            k_e_ssa = ne / k_e_ssa;

            for (size_t i = 0; i < gammmas.size()-1; i++) {
                power_e = Dermer09::single_electron_synch_power(B, epsilon, gammmas[i]);
                absorption += power_e * k_e_s_ssa[i] * (gammmas[i + 1] - gammmas[i]);
            }
            absorption *= k_e; // TODO check k_e_ssa
            double coeff = -1 / (8 * CGS::pi * CGS::me * std::pow(epsilon, 2)) * std::pow(CGS::lambda_c / CGS::c, 3);
            absorption *= coeff;

            abs = absorption;
        }
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchMARG21(double nuprime){
//        double gamma_min=p_pars->gamma_min, gamma_c=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
//        double Theta=p_pars->Theta, z_cool = p_pars->z_cool, acc_frac = p_pars->accel_frac;
//        double p = p_pars->p;
        double x, em_th=0., em_pl=0., abs_th=0., abs_pl=0.; // for Margalit model
        /// Margalit+21 arXiv:2111.00012
        ///
        double ne = n_prime * ksi_n;
        double ne_ = mu_e * ne;//4.0 * p_pars->mu_e * p_pars->n_ism;
        double delta = eps_e / eps_t; // defined in Margalit+21 arXiv:2111.00012
        /// normalized frequency:
        x = nuprime / Margalit21::nu_Theta(Theta, B);
        if (!std::isfinite(x)){
            (*p_log)(LOG_ERR,AT) << " x = "<< x << "\n";
            exit(1);
        }
        /// calculate total emissivity & optical depth:
        em_th = Margalit21::jnu_th(x, ne_, B, Theta, z_cool);
        em_pl = Margalit21::jnu_pl(x, ne_ * accel_frac, B, Theta, gamma_min, delta, p, z_cool); //TODO added 'accel_frac'
        if ((!std::isfinite(em_th))||(em_th < 0.)) em_th = 0.;
        if ((!std::isfinite(em_pl))||(em_pl < 0.)) em_pl = 0.;
        if ((em_pl==0.) && (em_th==0.)){
            (*p_log)(LOG_ERR,AT) << " em_pl=em_th=0 for"
                                 << " n_prime=" << n_prime << " accel_frac" << accel_frac
                                 << " B=" << B << " gamma_min=" << gamma_min << " Theta=" << Theta
                                 << " z_cool=" << z_cool << " nuprime=" << nuprime << " x=" << x
                                 << "\n";
            exit(1);
        }
//            emissivity = em_pl + em_th;
        if (m_methods_ssa!=iSSAoff) {
            abs_th = Margalit21::alphanu_th(x, ne_, B, Theta, z_cool);
            abs_pl = Margalit21::alphanu_pl(x, ne_ * accel_frac, B, Theta, gamma_min, delta, p, z_cool);
            if (!std::isfinite(abs_th)) abs_th = 0.;
            if (!std::isfinite(abs_pl)) abs_pl = 0.;
//                abs = abs_th + abs_pl;
        }

        em=em_pl+em_th;//m_data[i_em] = em_pl + em_th;
        abs=abs_pl+abs_th;//m_data[i_abs] = abs_th + abs_pl;
    }

    /// Analytical Synchrotron Sectrum; BPL;
    void checkEmssivityAbsorption(){
        if (( em < 0.) || (!std::isfinite( em )) ){
            (*p_log)(LOG_ERR,AT) << " em_pl_prime < 0 or nan ("<< em<<") or \n";
            (*p_log)(LOG_ERR,AT) << " abs_pl_prime < 0 or nan ("<< abs<<")\n";
//            (*p_log)(LOG_ERR,AT) << " Error in data \n"
//                                 << " eps_e = " << p_pars->eps_e << "\n"
//                                 << " eps_t = " << p_pars->eps_t << "\n"
//                                 << " ne = " << p_pars->ne << "\n"
//                                 << " gm = " << p_pars->gm << "\n"
//                                 << " gM = " << p_pars->gM << "\n"
//                                 << " gc = " << p_pars->gc << "\n"
//                                 << " B = " << p_pars->B << "\n"
//                                 << " Theta = " << p_pars->Theta << "\n"
//                                 << " z_cool = " << p_pars->z_cool << "\n"
//                                 << " nuprime = " << p_pars->nuprime << "\n";
            exit(1);
        }
    }

public: // ---------------- ANALYTIC -------------------------- //

    ElectronAndRadiation(int loglevel, bool _is_rs) : ElectronAndRadiaionBase(loglevel, _is_rs){}

    void setPars( StrDbMap & pars, StrStrMap & opts, size_t nr ){
        setBasePars(pars, opts);
        // ---
        switch (m_eleMethod) {
            case iShockEleAnalyt:
                allocateStorageForAnalyticSpectra(pars, nr);
                break;
            case iShockEleNum:
                allocateStorageForNumericSpectra(pars, nr);
                break;
        }
    }

    void allocateStorageForAnalyticSpectra(StrDbMap & pars, size_t nr){
        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e28, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,200, true);//pars.at("nfreq");

        (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                              << " freqs="<<nfreq << " by radii=" << nr << " Spec. grid="
                              << nfreq * nr << "\n";

        /// allocate space for the comoving spectrum
        out_spectrum.resize(m_freq_arr.size() * nr );
        out_specturm_ssa.resize(m_freq_arr.size() * nr );

        allocateStorageForOutputSpectra(nr, freq1, freq2, nfreq);
    }

    /// evaluate the comoving emissivity and absorption (frequency dependent)
    void computeSynchrotronEmissivityAbsorptionAnalytic( double nuprime ) {

        // TODO WARNING I did replace n_prime with ne is absorption, but this might not be correct!!!

        if (m_sychMethod == METHODS_SYNCH::iJOH06)
            computeAnalyticSynchJOH06(nuprime);
        else if (m_sychMethod == METHODS_SYNCH::iWSPN99)
            computeAnalyticSynchWSPN99(nuprime);
        else if (m_sychMethod == METHODS_SYNCH::iDER06)
            computeAnalyticSynchDER06(nuprime);
        else if (m_sychMethod == METHODS_SYNCH::iMARG21)
            computeAnalyticSynchMARG21(nuprime);
        else{
            (*p_log)(LOG_ERR,AT)<<" analytic synchrotron method is not supported \n";
            exit(1);
        }

        checkEmssivityAbsorption();
    }

    /// compute spectrum for all freqs and add it to the container
    void computeSynchrotronSpectrumAnalytic(size_t it){
        size_t nfreq = m_freq_arr.size();
        /// computeSynchrotronEmissivityAbsorptionAnalytic emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
            computeSynchrotronEmissivityAbsorptionAnalytic(m_freq_arr[ifreq]);
            out_spectrum[ifreq + nfreq * it] = em;
            out_specturm_ssa[ifreq + nfreq * it] = abs;
        }
    }

public: // -------------------- NUMERIC -------------------------------- //

    void allocateStorageForNumericSpectra(StrDbMap & pars, size_t nr){

        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e28, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,200, true);//pars.at("nfreq");

        double gam1 = getDoublePar("gam1", pars, AT, p_log,1, true);//pars.at("freq1");
        double gam2 = getDoublePar("gam2", pars, AT, p_log,1.e8, true);//pars.at("freq2");
        size_t ngams = (size_t)getDoublePar("ngam", pars, AT, p_log,250, true);//pars.at("nfreq");

        ele.allocate(gam1, gam2, ngams);
        ele.build_grid_chang_cooper(); // build the grid for the implicit solver

        syn.allocate(freq1*CGS::freqToErg, freq2*CGS::freqToErg, nfreq);

        if ((m_sychMethod != METHODS_SYNCH::iGSL) and (m_sychMethod != METHODS_SYNCH::iDER06)){
            (*p_log)(LOG_ERR,AT) << " only GSL or DER06 synchrotron options are "
                                    "avaialble when evolving electrons numerically. \n";
            exit(1);
        }

        /// allocate memory for synchrotron kernel
        if (m_sychMethod == METHODS_SYNCH::iGSL)
            syn_kernel.allocate(ele, syn, SynKernel::synchGSL);
        else
            syn_kernel.allocate(ele, syn, SynKernel::synchDermer);

        /// allocate memory for SSC structure and kernel
        if (m_methods_ssc == METHOD_SSC::iNumSSC) {
            ssc.allocate(freq1 * CGS::freqToErg, freq2 * CGS::freqToErg, nfreq);
            ssc_kernel.allocate(ele, syn, ssc, SSCKernel::sscNava);
            ssc_kernel.evalSSCkernel(ele, syn, ssc);
            ssc_kernel.evalSSCgridIntergral(ele, syn, ssc);
        }

        allocateStorageForOutputSpectra(nr, freq1, freq2, nfreq);

    }

    void evaluateElectronDistributionNumeric(double dt, double d_m, double r, double dr, double rp1, double drp1){

        /// check parameters
        if (ele.e.size() < 1){
            (*p_log)(LOG_ERR,AT) << " electon grid is not initialized\n";
            exit(1);
        }
        if (B <= 0 || gamma_min <= 0 || gamma_max <= 0){
            (*p_log)(LOG_ERR,AT) << "wrong value, cannot evolve electrons: "
                                    "B="<<B<<" gm="<<gamma_min<<" gM="<<gamma_max<<"\n";
            exit(1);
        }

        /// for adiabatic cooling of electron distribution
        vol = r * r * dr;
        vol_p1 = rp1 * rp1 * drp1;
        dm = d_m;
        double dlnVdt = 1. / dt * (1. - vol / vol_p1);
        /// number of injected electrons
        double N = dm / CGS::mp / dt;
        /// compute substepping to account for electron cooling timescales
        double dlog_gam = ((ele.e[1]-ele.e[0])/ele.e[0]);
        double delta_t_syn = CGS::sigmaT * gamma_max * B * B / (6. * M_PI * CGS::me * CGS::c);
        double delta_t_adi = (1.-gamma_max*gamma_max)/(3.*gamma_max*gamma_max)*dlnVdt;
        double delta_t = dlog_gam / (delta_t_syn + delta_t_adi);

//        int n_substeps = (int) (dt / delta_t);

        /// if cooling is too slow, we still need to evolve distribution
        int max_substeps = 1000;
        int n_substeps;
        if (delta_t <= dt){
            delta_t = dt;
            n_substeps = 1;
        }
        else if (delta_t > dt / (double)max_substeps){
            delta_t = dt / (double) max_substeps;
            n_substeps = (int) (dt / delta_t);
        }
        else{
            n_substeps = (int) (dt / delta_t);
        }
//
//        n_substeps = std::max(n_substeps, 1);
//        n_substeps = std::min(n_substeps, max_substeps);
//
//
//        /// if cooling is too fast we need to limit the maximum number of iteration
//        if (delta_t > dt)
//            delta_t = dt;
//        if (delta_t > dt / (double)max_substeps) {
//            delta_t = dt / (double) max_substeps;
//            n_substeps = max_substeps;
//        }


        delta_t = std::max(delta_t, dt / (double)max_substeps);

        (*p_log)(LOG_INFO, AT) << "\t nsubsteps="<<n_substeps<<" max="<<max_substeps<<"\n";

        /// update source properties
        source.B = B;
        source.dlnVdt = dlnVdt;
        source.r = r;
        source.dr = dr;
        source.N = N;

        ChangCooper model = ChangCooper(source, ele, syn, ssc, syn_kernel, ssc_kernel);

        /// Assume source/escape do not change during substepping
        model.setSourceFunction(gamma_min, gamma_max, -p, N);
        model.setEscapeFunction(gamma_min, gamma_max);

        /// Init photon field
        Vector photon_field (syn.numbins, 0);

        for (size_t i = 0; i < (size_t)n_substeps; i++){

            /// Set gamma_dot terms (uses syn.n[] for photon density for SSC from previous solution)
            model.setHeatCoolTerms(2);

            /// Assuming sources change during substepping
//            model.setSourceFunction(m_data[Q::igm][it], m_data[Q::igM][it], index, N);
//            model.setEscapeFunction(m_data[Q::igm][it], m_data[Q::igM][it]);

            /// setup grid of differences for escape term
            model.computeDeltaJ();

            /// setup tridiagonal solver arrays
            model.setupVectors(delta_t);

            /// solve system
            model.solve(delta_t);

            /// Update photon field during electron evolution
//            model.update_radiation(); // Too long, if inside

        }

        /// Update photon field during electron evolution
        model.update_radiation(); // saves photon density in syn.n[] -> used in ssc cooling next iteration

    }

    void computeSynchrotronSpectrumNumeric(size_t it){

//        m_sychMethod = METHODS_SYNCH::iDER06;

        size_t nfreq = m_freq_arr.size();
        /// store emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {

//            computeSynchrotronEmissivityAbsorptionAnalytic(m_freq_arr[ifreq]);

            out_spectrum[ifreq + nfreq * it] = syn.j[ifreq]/n_protons*n_prime;
            out_specturm_ssa[ifreq + nfreq * it] = syn.a[ifreq]/n_protons*n_prime;

//            std::cout<<ifreq<<" em="<<em<<" j="<<out_spectrum[ifreq + nfreq * it]<<"\n";
//            std::cout<<ifreq<<" em="<<abs<<" a="<<out_specturm_ssa[ifreq + nfreq * it]<<"\n";
//            int x = 1;

        }

        /// implicitely assume that SSC and Syn grids are the same. TODO generalize
        if (m_methods_ssc != METHOD_SSC::inoSSC)
            for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
                out_spectrum[ifreq + nfreq * it] += ssc.j[ifreq];
            }
    }
};

#if 0
class RadiationNumeric : public ElectronAndRadiation {
public:
    /// --------------------------------------
    Source source{};
    State ele{}; // std::unique_ptr<State> ele = nullptr;
    State syn{}; //std::unique_ptr<State> syn = nullptr;
    State ssc{}; //std::unique_ptr<State> ssc = nullptr;
    SynKernel syn_kernel{};//std::unique_ptr<SSCKernel> ssc_kernel = nullptr;
    SSCKernel ssc_kernel{};//std::unique_ptr<SynKernel> syn_kernel = nullptr;


    RadiationNumeric(int loglevel, bool _is_rs) : ElectronAndRadiaionBase(loglevel, _is_rs) {}

    void setPars(StrDbMap & pars, StrStrMap & opts){
        m_eleMethod = METHODS_SHOCK_ELE::iShockEleNum;
        setBasePars(pars, opts);
    }

    void allocate_storage(StrDbMap & pars, size_t nr){

        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e28, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,200, true);//pars.at("nfreq");

        double gam1 = getDoublePar("gam1", pars, AT, p_log,1, true);//pars.at("freq1");
        double gam2 = getDoublePar("gam2", pars, AT, p_log,1.e8, true);//pars.at("freq2");
        size_t ngams = (size_t)getDoublePar("ngam", pars, AT, p_log,250, true);//pars.at("nfreq");

        ele.allocate(gam1, gam2, ngams);
        ele.build_grid_chang_cooper(); // build the grid for the implicit solver

        syn.allocate(freq1*CGS::freqToErg, freq2*CGS::freqToErg, nfreq);

        if ((m_sychMethod != METHODS_SYNCH::iGSL) and (m_sychMethod != METHODS_SYNCH::iDER06)){
            (*p_log)(LOG_ERR,AT) << " only GSL or DER06 synchrotron options are "
                                    "avaialble when evolving electrons numerically. \n";
            exit(1);
        }

        /// allocate memory for synchrotron kernel
        if (m_sychMethod == METHODS_SYNCH::iGSL)
            syn_kernel.allocate(ele, syn, SynKernel::synchGSL);
        else
            syn_kernel.allocate(ele, syn, SynKernel::synchDermer);

        /// allocate memory for SSC structure and kernel
        if (m_methods_ssc == METHOD_SSC::iNumSSC) {
            ssc.allocate(freq1 * CGS::freqToErg, freq2 * CGS::freqToErg, nfreq);
            ssc_kernel.allocate(ele, syn, ssc, SSCKernel::sscNava);
            ssc_kernel.evalSSCkernel(ele, syn, ssc);
            ssc_kernel.evalSSCgridIntergral(ele, syn, ssc);
        }

        allocate_output_storage(nr, freq1, freq2, nfreq);

    }

    void evaluateElectronDistributionNumeric(double dt, double dm, double r, double dr, double rp1, double drp1){

        /// check parameters
        if (ele.e.size() < 1){
            (*p_log)(LOG_ERR,AT) << " electon grid is not initialized\n";
            exit(1);
        }
        if (B <= 0 || gamma_min <= 0 || gamma_max <= 0){
            (*p_log)(LOG_ERR,AT) << "wrong value, cannot evolve electrons: "
                                    "B="<<B<<" gm="<<gamma_min<<" gM="<<gamma_max<<"\n";
            exit(1);
        }

        /// for adiabatic cooling of electron distribution
        double dlnVdt = 1. / dt * (1. - (r * r * dr) / (rp1 * rp1 * drp1));
        /// number of injected electrons
        double N = dm / CGS::mp / dt;
        /// compute substepping to account for electron cooling timescales
        double dlog_gam = ((ele.e[1]-ele.e[0])/ele.e[0]);
        double delta_t_syn = CGS::sigmaT * gamma_max * B * B / (6. * M_PI * CGS::me * CGS::c);
        double delta_t_adi = (1-gamma_max*gamma_max)/(3.*gamma_max*gamma_max)*dlnVdt;
        double delta_t = dlog_gam / (delta_t_syn + delta_t_adi);

        double n_substeps = dt / delta_t;
        /// if cooling is too slow, we still need to evolve distribution
        n_substeps = std::max(int(n_substeps), 1);

        /// if cooling is too fast we need to limit the maximum number of iteration
        if (delta_t > dt) delta_t = dt;
        delta_t = std::max(delta_t, dt / 1000);

        /// update source properties
        source.B = B;
        source.dlnVdt = dlnVdt;
        source.r = r;
        source.dr = dr;
        source.N = N;

        ChangCooper model = ChangCooper(source, ele, syn, ssc, syn_kernel, ssc_kernel);

        /// Assume source/escape do not change during substepping
        model.setSourceFunction(gamma_min, gamma_max, p, N);
        model.setEscapeFunction(gamma_min, gamma_max);

        /// Init photon field
        Vector photon_field (syn.numbins, 0);

        for (size_t i = 0; i < (size_t)n_substeps; i++){

            /// Set gamma_dot terms (uses syn.n[] for photon density for SSC from previous solution)
            model.setHeatCoolTerms(2);

            /// Assuming sources change during substepping
//            model.setSourceFunction(m_data[Q::igm][it], m_data[Q::igM][it], index, N);
//            model.setEscapeFunction(m_data[Q::igm][it], m_data[Q::igM][it]);

            /// setup grid of differences for escape term
            model.computeDeltaJ();

            /// setup tridiagonal solver arrays
            model.setupVectors(delta_t);

            /// solve system
            model.solve(delta_t);

            /// Update photon field during electron evolution
//            model.update_radiation(); // Too long, if inside

        }

        /// Update photon field during electron evolution
        model.update_radiation(); // saves photon density in syn.n[] -> used in ssc cooling next iteration

    }

    void computeSynchrotronSpectrumNumeric(size_t it){
        size_t nfreq = m_freq_arr.size();
        /// store emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
            out_spectrum[ifreq + nfreq * it] = syn.j[ifreq];
            out_specturm_ssa[ifreq + nfreq * it] = syn.a[ifreq];
        }

        /// implicitely assume that SSC and Syn grids are the same. TODO generalize
        if (m_methods_ssc != METHOD_SSC::inoSSC)
            for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
                out_spectrum[ifreq + nfreq * it] += ssc.j[ifreq];
            }
    }
};
#endif


#if 0
/// evaluate comoving emissivity and absorption for synchrotron mech.
class SynchrotronAnalytic : private ElectronAndRadiaionBase{

    int m_loglevel = -1;

    /// methods


    /// parameter container for the class
    struct Pars{
        // --- in
        double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
        double mu=-1, mu_e=-1;
        bool lim_gm_to_1= true;
        double beta_min = -1;
        // --- methods
        METHODS_SYNCH m_sychMethod{};
        METHODS_LFMIN m_methodsLfmin{};
        METHODS_B m_methodsB{};
        METHOD_LFMAX m_methodsLfmax{};
        METHOD_NONRELDIST m_method_nonreldist{};
        METHOD_TAU method_tau{};
//        METHODS_RAD method_comp_mode{};

        METHODS_SSA m_methods_ssa{};

        // --- out
        double B=-1, gamma_min=-1, gamma_max=-1, gamma_c=-1;
        double n_prime=-1, eprime=-1.,Gamma=-1,GammaSh=-1, beta=-1.,t_e=-1.;
        double accel_frac=-1.;
        double Theta=-1, z_cool=-1, x=-1;

//        double em=0.,em_th=0.,em_pl=0.;
//        double abs=0.,abs_th=0.,abs_pl=0.;

        double em=-1., abs=-1.;

    };
//    std::unique_ptr<Pars> p_pars = nullptr;
    Pars * p_pars = nullptr;
    std::unique_ptr<logger> p_log = nullptr;

    /// -------------------------------------------------------
    bool is_rs = false;
    /// -------------------------------------------------------

public:

    Vector m_freq_arr{}; Vector out_spectrum{}; Vector out_specturm_ssa{};

    SynchrotronAnalytic( int loglevel, bool _is_rs ){
        m_loglevel = loglevel;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronAnalytic");
        p_pars = new Pars();
//        p_pars = std::make_unique<Pars>();
        p_pars->lim_gm_to_1 = false;
//        m_data.resizeEachImage( Rad::m_names_.size(), -1. );
        is_rs = _is_rs;
    }
    ~SynchrotronAnalytic(){ delete p_pars; }

//    std::unique_ptr<Pars> & getPars(){ return p_pars; }
    Pars *& getPars(){ return p_pars; }

    /// set model parameters
    void setPars(StrDbMap & pars, StrStrMap & opts){

        // set parameters
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";

        p_pars->ksi_n = getDoublePar("ksi_n" + fs_or_rs, pars, AT, p_log, 1., false);//pars.at("ksi_n");
        p_pars->eps_e = getDoublePar("eps_e" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_e");
        p_pars->eps_b = getDoublePar("eps_b" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_b");
        p_pars->eps_t = getDoublePar("eps_t" + fs_or_rs, pars, AT, p_log, 0., true);//pars.at("eps_t");
        p_pars->p = getDoublePar("p" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("p");
        p_pars->mu = getDoublePar("mu" + fs_or_rs, pars, AT, p_log, 0.62, false);//pars.at("mu");
        p_pars->mu_e = getDoublePar("mu_e" + fs_or_rs, pars, AT, p_log, 1.18, false);//pars.at("mu_e");
        p_pars->beta_min = getDoublePar("beta_min" + fs_or_rs, pars, AT, p_log, 1.e-5, false);//pars.at("beta_min");
        p_pars->gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars, AT, p_log, 1.e7, false);//pars.at("beta_min");

        // set options
        std::string opt;
        opt = "method_synchrotron" + fs_or_rs;
        METHODS_SYNCH val_synch;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_synch = SynchrotronAnalytic::METHODS_SYNCH::iJOH06;
        }
        else{
            if(opts.at(opt) == "Joh06")
                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iJOH06;
            else if(opts.at(opt) == "WSPN99")
                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iWSPN99;
            else if(opts.at(opt) == "Marg21")
                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iMARG21;
            else if(opts.at(opt) == "Dermer09")
                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iDER06;
            else if(opts.at(opt) == "New")
                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iNumeric;
//            else if(opts.at(opt) == "Bretta")
//                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iBerrettaSynch;
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

        opt = "method_nonreldist" + fs_or_rs;
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

        opt = "method_lf_min" + fs_or_rs;
        SynchrotronAnalytic::METHODS_LFMIN val_lfmin;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmUprime;
        }
        else{
            if(opts.at(opt) == "useU_e")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmUprime;
            else if(opts.at(opt) == "useNumericGamma")
                val_lfmin = SynchrotronAnalytic::METHODS_LFMIN::igmNumGamma;
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

        opt = "method_Bsh" + fs_or_rs;
        SynchrotronAnalytic::METHODS_B methodsB;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsB = SynchrotronAnalytic::METHODS_B::iuseUb;
        }
        else{
            if(opts.at(opt) == "useUb")
                methodsB = SynchrotronAnalytic::METHODS_B::iuseUb;
            else if(opts.at(opt) == "useMAG21")
                methodsB = SynchrotronAnalytic::METHODS_B::iasMAG21;
            else if(opts.at(opt) == "useGammaSh")
                methodsB = SynchrotronAnalytic::METHODS_B::iuseGammaSh;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useUb " << " useMAG21 " << " useGammaSh " << "\n";
                exit(1);
            }
        }
        p_pars->m_methodsB = methodsB;

        opt = "method_lf_max" + fs_or_rs;
        SynchrotronAnalytic::METHOD_LFMAX methodLfmax;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodLfmax = SynchrotronAnalytic::METHOD_LFMAX::iuseB;
        }
        else{
            if(opts.at(opt) == "useB")
                methodLfmax = SynchrotronAnalytic::METHOD_LFMAX::iuseB;
            else if(opts.at(opt) == "useConst") {
                p_pars->gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars,
                                                 AT, p_log, 1.e7, true);//pars.at("beta_min");
                methodLfmax = SynchrotronAnalytic::METHOD_LFMAX::iConst;
            }
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useB " << " useConst " << "\n";
                exit(1);
            }
        }
        p_pars->m_methodsLfmax = methodLfmax;

        bool tmp = getBoolOpt("use_ssa" + fs_or_rs, opts, AT, p_log, false, true);
        if (tmp) p_pars->m_methods_ssa = SynchrotronAnalytic::METHODS_SSA::iSSAon;
        else p_pars->m_methods_ssa = SynchrotronAnalytic::METHODS_SSA::iSSAoff;

#if 0
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
#endif
        opt = "method_tau" + fs_or_rs;
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
        (*p_log)(LOG_INFO,AT) << " allocating memory for gamma_arr\n";

    }

    void allocateSpectrumStorage(size_t nr, StrDbMap & pars, StrStrMap & opts){

        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e14, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,100, true);//pars.at("nfreq");

        /// allocate storage feoe

        ///

        (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                              << " freqs="<<nfreq << " by radii=" << nr << " Spec. grid="
                              << nfreq * nr << "\n";

        /// allocate space for the comoving spectrum
        m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);

        out_spectrum.resize(m_freq_arr.size() * nr );
        out_specturm_ssa.resize(m_freq_arr.size() * nr );

    }

    bool checkParams(){

        /// if velocity is too small shock may not be possible TODO update it with proper sound speed check
        if (p_pars->beta < p_pars->beta_min)
            return false;


        /// in Margalit+21 model thermal electrons are present. Different treatment :: Assert:
        if (p_pars->m_sychMethod == METHODS_SYNCH::iMARG21) {
            /* Margalit+21 arXiv:2111.00012 */
            if (p_pars->m_methodsLfmin != igmMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsLfmin=gmMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
            if (p_pars->m_methodsB != iasMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsB=asMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
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
        if (p_pars->n_prime <= 0){
            (*p_log)(LOG_ERR,AT) << " n_prime is not set (n_prime="<<p_pars->n_prime<<")\n";
            exit(1);
        }
        if ((p_pars->eprime < 0) || !std::isfinite(p_pars->eprime)){
            (*p_log)(LOG_ERR,AT) << " eprime is not set (eprime=" << p_pars->eprime << ")\n";
            exit(1);
        }

        return true;
    }


    static double gammaMinFunc(const double &x, void * pars){
        auto * pp = (struct Pars *) pars;
        return (pp->p - 1) / (pp->p - 2) * (std::pow(pp->gamma_max, -pp->p + 2) - std::pow(x, -pp->p + 2))
               / (std::pow(pp->gamma_max, -pp->p + 1) - std::pow(x, -pp->p + 1)) - pp->eps_e * mp / me * (pp->GammaSh - 1);
    }


    void computeMagneticField(){
        /// fraction of the shock energy density in electrons
        double U_b_prime;
        switch (p_pars->m_methodsB) {
            case iuseUb:
                U_b_prime = p_pars->eprime * p_pars->eps_b;
                p_pars->B = sqrt(8. * CGS::pi * U_b_prime);
                break;
            case iuseGammaSh:
                /// See Nava 2013 paper or others; classical equation
                p_pars->B = std::sqrt( 32. * M_PI * p_pars->eps_b * (p_pars->GammaSh - 1.)
                               * (p_pars->GammaSh + 3. / 4) * p_pars->n_prime * CGS::mp * CGS::c2);
                break;
            case iasMAG21:
                p_pars->B = sqrt(9.0 * M_PI * p_pars->eps_b * p_pars->n_prime * p_pars->mu * CGS::mp)
                                * (p_pars->beta * CGS::c);
                break;
        }
    }

    void computeGammaMax(){
        /// get electron distribution boundaries (has ot be set BEFORE gamma_min)
        switch (p_pars->m_methodsLfmax) {
            case iConst:
                break;
            case iuseB:
                p_pars->gamma_max = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / p_pars->B); // Kumar+14
                break;
        }
    }

    void computeGammaMin(){
        /// compute injection LF (lower bound of source injection gammaMinFunc)
        double gm = 0, Theta=-1;
        int status = 0; double U_e_prime = -1;

        double Gamma_shock = p_pars->GammaSh;
        double eps_e = p_pars->eps_e;
        double p = p_pars->p;

        switch (p_pars->m_methodsLfmin) {
            case igmUprime:
                U_e_prime = p_pars->eprime * eps_e;
                gm = (p - 2.) / (p - 1.) * U_e_prime / (p_pars->n_prime * CGS::me * CGS::c * CGS::c); // Eq. A18 vanEarten+10 (and Sironi+13)
                break;
            case igmNakarPiran:
                gm = (p - 2.) / (p - 1.) * (CGS::mp / CGS::me) * eps_e * EQS::Beta2(Gamma_shock); // Eq.(before 4.) in Nakar&Piran 1102.1020
                break;
            case igmJoh06:
                gm = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (Gamma_shock - 1.) + 1.); // Eq. A3 in J+06
                break;
            case igmNumGamma:
                /// solve gamma_min fun numerically; use fixed limits and number of iterations
                gm = Bisect(gammaMinFunc, 1, 1e8, 0, .001, 100, p_pars, status);
                /// If numerical solution failed, use simple analytical solution
                if (status < 0)
                    gm = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (Gamma_shock - 1.) + 1.);
                break;
            case igmMAG21:
                /// downstream electron temperature:
                Theta = Margalit21::Theta_fun(p_pars->beta, p_pars->mu, p_pars->mu_e, p_pars->eps_t);
                gm = Margalit21::gamma_m_fun(Theta);
                break;
        }

        /// check
        if (!std::isfinite(gm)){
            (*p_log)(LOG_ERR,AT) << " error gm nan \n";
            exit(1);
        }

        /// limit the min lf to 1
        if ((p_pars->lim_gm_to_1) && (gm < 1.))
            gm = 1.; // Sironi et al 2013 suggestion to limi gm=1. for Deep Newtinoan regime # TODO to be removed. I did no understand it

        p_pars->gamma_min = gm;
        p_pars->Theta = Theta;

        /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21): NOTE we use mean dynamical time:
        p_pars->z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * p_pars->B * p_pars->B * p_pars->t_e)) / p_pars->Theta;

        if (!std::isfinite(p_pars->z_cool)){
            (*p_log)(LOG_ERR,AT) << AT << " Theta = " << p_pars->Theta << " z_cool="<<p_pars->z_cool<<"\n";
            exit(1);
        }
    }



    /// store current shock properties
    void updateSockProperties(double eprime, double Gamma, double Gamma_shock, double t_e, double n_prime){
        /// Store current parameters
        p_pars->GammaSh = Gamma_shock; // for bisect solver
        p_pars->eprime = eprime; // for bisect solver
        p_pars->Gamma = Gamma; // for bisect solver
        p_pars->beta = EQS::Beta(Gamma_shock); // for bisect solver
        p_pars->t_e = t_e; // for bisect solver
        p_pars->n_prime = p_pars->ksi_n * n_prime; //  for bisect solver
    }

    /// evaluate frequency independent quantities (critical LFs, Bfield, etc)
    void evaluateElectronDistribution() {

        if (not checkParams())
            return;

        computeMagneticField();

        computeGammaMax();

        computeGammaMin();

        computeNonRelativisticFlattening();

        /// compute cooling lorentz factor
        p_pars->gamma_c = 6. * CGS::pi * CGS::me * CGS::c
                        / (CGS::sigmaT * p_pars->t_e * p_pars->B * p_pars->B)
                        / p_pars->Gamma; // Eq. A19 in vanEarten+10

    }


    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchJOH06(double nuprime){
        double gm=p_pars->gamma_min, gc=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
        double p = p_pars->p;
        double em=0., abs=0.;
        // prefactors
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
        double nu_m, nu_c, emissivity, scaling, abs_scaling=1.;

        double ne = n_prime * p_pars->ksi_n;

        if (gm < gc){
            // slow cooling
            nu_m = XpS * gm * gm * gamToNuFactor;
            nu_c = XpS * gc * gc * gamToNuFactor;
//                double _phip = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * phipS;
            emissivity = 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p)// phipS
                         * CGS::qe * CGS::qe * CGS::qe * n_prime * B / (CGS::me * CGS::c * CGS::c);
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
            emissivity = _phip * CGS::qe * CGS::qe * CGS::qe * n_prime * B / (CGS::me * CGS::c * CGS::c);
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

        p_pars->em = emissivity * scaling;
        p_pars->abs = abs * abs_scaling;
    };
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchWSPN99(double nuprime){
        double gm=p_pars->gamma_min, gc=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
        double p = p_pars->p;
        double nu_m, nu_c, emissivity, scaling, abs_scaling=1.;
        double em=0.,abs=0.;
//            double PhiP = interpolated_phi(p);
//            emissivity = (interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi)
//                    * n_prime * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
        double ne = n_prime * p_pars->ksi_n;
        emissivity = (interpolated_phi(p) * sqrt(3.0) / 4.0 * CGS::pi) * ne * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
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
        p_pars->em = emissivity * scaling;
        p_pars->abs = abs * abs_scaling;
//            std::cout << 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p) << " " << interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi << "\n";
//            exit(1);
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchDER06(double nuprime){
        double gm=p_pars->gamma_min, gc=p_pars->gamma_c, gM=p_pars->gamma_max, B=p_pars->B, n_prime=p_pars->n_prime;
        double ne = n_prime * p_pars->ksi_n;
        double p = p_pars->p;
        double em=0.,abs=0.;
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
        // computeShycnhrotronEmissivityAbsorption electron distribution normalisation
        double k_e = ne / Simpson38(gm, gM, 200, integrand_ele); // TODO replace with adative integrals
        // convolve the electron distribution with emission spectrum and itegrate
        auto integrand = [&](double gam){
            double power_e = Dermer09::single_electron_synch_power( B, epsilon, gam );
            double n_e = k_e * integrand_ele(gam);
            return n_e * power_e;
        };
        double power = Simpson38(gm, gM, 200, integrand); // TODO replace with adative integrals
        em = power * (CGS::h / CGS::mec2);
        /* --- SSA ---
         * Computes the syncrotron self-absorption opacity for a general set
         * of model parameters, see
         * :gammaMinFunc:`~agnpy:sycnhrotron.Synchrotron.evaluate_sed_flux`
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
            // computeShycnhrotronEmissivityAbsorption electron distribution normalisation
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
        p_pars->em = em;
        p_pars->abs = abs;
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchMARG21(double nuprime){
        double gm=p_pars->gamma_min, gc=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
        double Theta=p_pars->Theta, z_cool = p_pars->z_cool, acc_frac = p_pars->accel_frac;
        double p = p_pars->p;
        double x, em_th=0., em_pl=0., abs_th=0., abs_pl=0.; // for Margalit model
        /// Margalit+21 arXiv:2111.00012
        ///
        double ne = n_prime * p_pars->ksi_n;
        double ne_ = p_pars->mu_e * ne;//4.0 * p_pars->mu_e * p_pars->n_ism;
        double delta = p_pars->eps_e / p_pars->eps_t; // defined in Margalit+21 arXiv:2111.00012
        /// normalized frequency:
        x = nuprime / Margalit21::nu_Theta(Theta, B);
        if (!std::isfinite(x)){
            (*p_log)(LOG_ERR,AT) << " x = "<< x << "\n";
            exit(1);
        }
        /// calculate total emissivity & optical depth:
        em_th = Margalit21::jnu_th(x, ne_, B, Theta, z_cool);
        em_pl = Margalit21::jnu_pl(x, ne_ * acc_frac, B, Theta, gm, delta, p, z_cool); //TODO added 'accel_frac'
        if ((!std::isfinite(em_th))||(em_th < 0.)) em_th = 0.;
        if ((!std::isfinite(em_pl))||(em_pl < 0.)) em_pl = 0.;
        if ((em_pl==0.) && (em_th==0.)){
            (*p_log)(LOG_ERR,AT) << " em_pl=em_th=0 for"
                                 << " n_prime=" << n_prime << " acc_frac" << acc_frac
                                 << " B="<<B<< " gm="<<gm<< " Theta="<<Theta
                                 << " z_cool="<<z_cool<< " nuprime="<<nuprime<< " x="<<x
                                 <<"\n";
            exit(1);
        }
//            emissivity = em_pl + em_th;
        if (p_pars->m_methods_ssa!=iSSAoff) {
            abs_th = Margalit21::alphanu_th(x, ne_, B, Theta, z_cool);
            abs_pl = Margalit21::alphanu_pl(x, ne_ * acc_frac, B, Theta, gm, delta, p, z_cool);
            if (!std::isfinite(abs_th)) abs_th = 0.;
            if (!std::isfinite(abs_pl)) abs_pl = 0.;
//                abs = abs_th + abs_pl;
        }

        p_pars->em=em_pl+em_th;//m_data[i_em] = em_pl + em_th;
        p_pars->abs=abs_pl+abs_th;//m_data[i_abs] = abs_th + abs_pl;
        p_pars->x = x;
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void checkEmssivityAbsorption(){
        if (( p_pars->em < 0.) || (!std::isfinite( p_pars->em )) ){
            (*p_log)(LOG_ERR,AT) << " em_pl_prime < 0 or nan ("<< p_pars->em<<") or \n";
            (*p_log)(LOG_ERR,AT) << " abs_pl_prime < 0 or nan ("<< p_pars->abs<<")\n";
//            (*p_log)(LOG_ERR,AT) << " Error in data \n"
//                                 << " eps_e = " << p_pars->eps_e << "\n"
//                                 << " eps_t = " << p_pars->eps_t << "\n"
//                                 << " ne = " << p_pars->ne << "\n"
//                                 << " gm = " << p_pars->gm << "\n"
//                                 << " gM = " << p_pars->gM << "\n"
//                                 << " gc = " << p_pars->gc << "\n"
//                                 << " B = " << p_pars->B << "\n"
//                                 << " Theta = " << p_pars->Theta << "\n"
//                                 << " z_cool = " << p_pars->z_cool << "\n"
//                                 << " nuprime = " << p_pars->nuprime << "\n";
            exit(1);
        }
    }



    /// In case the electrons were computed elsewhere E.g., if interpolated for EATS plane
    void setShockElectronParameters(double n_prime, double acc_frac,
                                    double B, double gm, double gM, double gc,
                                    double Theta, double z_cool){

        if (!std::isfinite(gm) || !std::isfinite(gc) || !std::isfinite(n_prime)) {
            (*p_log)(LOG_ERR, AT) << " nans is synchrotron spectum\n";
            exit(1);
        }

        p_pars->gamma_min = gm;
        p_pars->gamma_max = gM;
        p_pars->n_prime = n_prime;
        p_pars->B = B;
        p_pars->accel_frac = acc_frac;
        p_pars->gamma_c = gc;
        p_pars->z_cool = z_cool;
        p_pars->Theta = Theta;
    }

    /// evaluate the comoving emissivity and absorption (frequency dependent)
    void computeShycnhrotronEmissivityAbsorption(double nuprime ) {

        // TODO WARNING I did replace n_prime with ne is absorption, but this might not be correct!!!

        if (p_pars->m_sychMethod == METHODS_SYNCH::iJOH06)
            computeAnalyticSynchJOH06(nuprime);
        else if (p_pars->m_sychMethod == METHODS_SYNCH::iWSPN99)
            computeAnalyticSynchWSPN99(nuprime);
        else if (p_pars->m_sychMethod == METHODS_SYNCH::iDER06)
            computeAnalyticSynchDER06(nuprime);
        else if (p_pars->m_sychMethod == METHODS_SYNCH::iMARG21)
            computeAnalyticSynchMARG21(nuprime);
        else{
            (*p_log)(LOG_ERR,AT)<<" analytic synchrotron method is not supported \n";
            exit(1);
        }

        checkEmssivityAbsorption();
    }

    /// compute spectrum for all freqs and add it to the container
    void computeSynchrotronSpectrum(size_t it){
        size_t nfreq = m_freq_arr.size();
        /// computeShycnhrotronEmissivityAbsorption emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
            computeSynchrotronEmissivityAbsorptionAnalytic( m_freq_arr[ifreq] );
            out_spectrum[ifreq + nfreq * it] = p_pars->em;
            out_specturm_ssa[ifreq + nfreq * it] = p_pars->abs;
        }
    }

};
#endif






// ---------------
enum METHOD_PWN_SPEC { inumBPL, ianaBPL };
static void initialize_e_dis(double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt,
                             double *dgam_dt, double *tad, double *tsyn, double gam_max, int Nbin_e){
    double dln_gam = log(gam_max)/(double)(Nbin_e-1);
    int i;
    for (i=0;i<Nbin_e;i++){
        gam[i] = exp(dln_gam*i);
        dgam[i] = gam[i]*(exp(dln_gam)-1.);
        dN_dgam[i] = 0.;
        dN_dgam_dt[i] = 0.;
        dgam_dt[i] = 0.;
        tad[i] = 0.;
        tsyn[i] = 0.;
    }
}
static void initialize_ph_dis(Vector & freqs, double *gam_ph, double *P_nu_syn, double *alpha_nu_syn){
    int i;
//    double del_ln_gam_ph = (log(gam_ph_max)-log(gam_ph_min))/(double)(Nbin_ph-1);
    for (i=0;i<freqs.size();i++){
        gam_ph[i] = freqs[i] / (CGS::MeC2 / CGS::H); // Hz->erg //gam_ph_min*exp(del_ln_gam_ph*(double)i);
        P_nu_syn[i] = 0.;
        alpha_nu_syn[i] = 0.;
    }
}
static double syn_func_fit(double x)
{
    /* analytical fitting of synchrotron function F(x) */
    /* see http://arxiv.org/pdf/1301.6908.pdf */

    double F1 = M_PI*pow(2.0,5.0/3.0)/sqrt(3.0)/GAMMA13*pow(x,1.0/3.0);
    double F2 = sqrt(M_PI/2.0)*exp(-x)*pow(x,1.0/2.0);

    double a1_1 = -0.97947838884478688;
    double a1_2 = -0.83333239129525072;
    double a1_3 = 0.1554179602681624;
    double H_1 = a1_1*pow(x,1.0)+a1_2*pow(x,1.0/2.0)+a1_3*pow(x,1.0/3.0);
    double delta_1 = exp(H_1);

    double a2_1 = -0.0469247165562628882;
    double a2_2 = -0.70055018056462881;
    double a2_3 = 0.0103876297841949544;

    double H_2 = a2_1*pow(x,1.0)+a2_2*pow(x,1.0/2.0)+a2_3*pow(x,1.0/3.0);
    double delta_2 = 1.0-exp(H_2);

    return F1*delta_1+F2*delta_2;
}
static void calc_syn_spec(const double B, const double dr, const double vol,
                          const double *gam, const double *dgam, const double *dN_dgam,
                          const double gam_max, const int Nbin_e, const double *gam_ph,
                          double *P_nu_syn, double *alpha_nu_syn, const int Nbin_ph){

//    double nu=0.,x=0.,sin_alpha=2./3.,tau_sa=0.;
//    double integ=0.,integ_alpha=0.;
//#pragma omp parallel for default(shared) reduction (P_nu_syn,alpha_nu_syn) firstprivate(nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol)// private(nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol) shared(std::cout,std::cerr,B,r,dr,gam,dgam,dN_dgam,gam_max,Nbin_e,Nbin_ph,gam_ph,P_nu_syn,alpha_nu_syn,CGS::c,CGS::MeC2,CGS::H,CGS::ELEC,CGS::me) num_threads( 6 ) // private(i,nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol) shared(B,r,dr,gam,dgam,dN_dgam,gam_max,Nbin_e,gam_ph,P_nu_syn,alpha_nu_syn,Nbin_ph) num_threads( 6 )
#pragma omp parallel for num_threads( 6 )
    for (size_t k=0;k<Nbin_ph;k++) {
        double integ = 0.0;
        double sin_alpha=2./3.;
        double integ_alpha = 0.0;
        double nu = gam_ph[k]*CGS::MeC2/CGS::H;
//#pragma omp parallel for firstprivate(integ,integ_alpha,nu,x) num_threads( 10 )
        for (size_t i=0;i<Nbin_e;i++) {
            double x = (2.0*M_PI*nu)/(3.0*CGS::ELEC*gam[i]*gam[i]*B/2.0/CGS::me/CGS::c*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */

            //            std::cout << x << ", ";
            if ((i==0) || (i==Nbin_e-1)) {
                integ += 0.5*dN_dgam[i]*dgam[i]*syn_func_fit(x);
                integ_alpha += -0.5*sin_alpha*std::pow(gam[i],2.0)
                               * (-dN_dgam[i]/std::pow(gam[i],2.0))
                               / dgam[i]*syn_func_fit(x)
                               * dgam[i] / CGS::MeC2;
            }
            else {
                integ += dN_dgam[i]*dgam[i]*syn_func_fit(x);
                integ_alpha += -sin_alpha*std::pow(gam[i],2.0)*
                               (dN_dgam[i+1]/std::pow(gam[i+1],2.0) - dN_dgam[i]/std::pow(gam[i],2.0))
                               / dgam[i] * syn_func_fit(x) * dgam[i] / CGS::MeC2;
            }
//                if ((!std::isfinite(P_nu_syn[k]))||(P_nu_syn[k]<0)) {
//                    std::cerr << AT << " nan in pwn spec P_nu_syn[i]=" << P_nu_syn[k] << "\n";
//                    std::cout << " x=" << x << " integ=" << integ << " nu=" << nu << " dN_dgam[0]=" << dN_dgam[0]
//                              << " dgam[0]=" << dgam[0] << " gam[0]=" << gam[0] << " B=" << B
//                              << " vol=" << vol << " tau_sa=" << tau_sa << " alpha_nu_sym[0]=" << alpha_nu_syn[0] << "\n";
//                    exit(1);
//                }
        }

        P_nu_syn[k] = sqrt(3.0)*pow(CGS::ELEC,3.0)*B*sin_alpha/CGS::MeC2*integ; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
        alpha_nu_syn[k] = CGS::c*CGS::c/8.0/M_PI/nu/nu*sqrt(3.0)*pow(CGS::ELEC,3.0)*B*integ_alpha/vol; /* Eq. (6.52) of Rybicki & Lightman */

        double tau_sa = alpha_nu_syn[k] * dr;
        if (tau_sa > 1.0e-6){
            P_nu_syn[k] = (1.0-exp(-tau_sa)) * P_nu_syn[k] / tau_sa;
        }

        integ = 0.0;
        integ_alpha = 0.0;


    }
}
static double Ntot(double *dgam, double *dN_dgam, int Nbin_e) {
    int i;
    double tmp=0.;
    for (i=0;i<Nbin_e-1;i++){
        tmp += dN_dgam[i]*dgam[i];
    }
    tmp += 0.5*dN_dgam[Nbin_e-1]*dgam[Nbin_e-1];
    return tmp;
}
static double dgam_dt_ad(double gam, double t) {
    return gam/t;
}
static double dgam_dt_syn(double gam, double B) {
    // electron synchrotron energy loss rate (see e.g., Eq. 7.13 of Dermer & Menon)
    // double sin2phi = 2.0/3.0; /* averaging pitch angle */
    // double beta_par = 1.0; /* assuming that particles are relativistic */

    return 4.0/3.0*CGS::c*CGS::SIGMA_T*(B*B/8.0/M_PI)*gam*gam/CGS::MeC2;
}
static void cooling(double t, double B, double *dgam_dt, double *gam, double *tad, double *tsyn, int Nbin_e){
    int i;
    for (i=0;i<Nbin_e;i++) {
        dgam_dt[i] = dgam_dt_ad(gam[i],t)+dgam_dt_syn(gam[i],B);
        tad[i] = gam[i]/dgam_dt_ad(gam[i],t);
        tsyn[i]= gam[i]/dgam_dt_syn(gam[i],B);
    }
}
static double dN_dgam_dt_inj(double gam, double Lpsr, double epse, double gam_b, double gam_max, double p1, double p2){
    double fac_bol = 1./(2. - p1) * (1. - pow(gam_b, -2. + p1))
                   + 1./(2. - p2) * (pow(gam_max/gam_b, 2. - p2) - 1.);
    if (gam < gam_b){
        return epse * Lpsr/(MeC2*gam_b*gam_b) / fac_bol * pow(gam/gam_b,-p1);
    }
    else {
        return epse * Lpsr/(MeC2*gam_b*gam_b) / fac_bol * pow(gam/gam_b,-p2);
    }
}
static void time_evolution_e(double dt, double *gam, double *dgam, double *dN_dgam,
                             double *dN_dgam_dt, double *dgam_dt, int Nbin_e)
{
    int i;
    double dN_dgam_old[Nbin_e];
    double N_cool = 0.0;
    double dN_dgam_cool = 0.0;

    for (i=0; i<Nbin_e; i++){
        dN_dgam_old[i] = dN_dgam[i];
    }

    dN_dgam[Nbin_e-1] = (dN_dgam_old[Nbin_e-1])
                      / (1.0+dt/dgam[Nbin_e-1]*dgam_dt[Nbin_e-1])
                      + dN_dgam_dt[Nbin_e-1]*dt;
    for(i = Nbin_e-2; i>0; i--){
        dN_dgam[i] = (dN_dgam_old[i]+dN_dgam[i+1]*dt/dgam[i]*dgam_dt[i+1]) / (1.0+dt/dgam[i]*dgam_dt[i])
                   + dN_dgam_dt[i] * dt;
    }
    dN_dgam[0] = dN_dgam_old[0]
               + dN_dgam_dt[1] * dt / dgam[0] * dgam_dt[1]
               + dN_dgam_dt[0] * dt;
}
static void injection(double *gam, double *dgam, double *dN_dgam_dt, double *dgam_dt, double Lpsr, double dt,
                      double *N_inj_tot, double epse, double gam_b, double gam_max,
                      double p1, double p2, int Nbin_e){
    int i,j;
    double tmp=0.0,total_frac_depo=0.0,frac_resi=0.0,frac_depo=0.0;
    double dN_dgam_dt_eff[Nbin_e];

    /// standard integration
    for (i=0;i<Nbin_e;i++){
        dN_dgam_dt[i] = dN_dgam_dt_inj(gam[i],Lpsr,epse,gam_b,gam_max,p1,p2);
        dN_dgam_dt_eff[i] = 0.0;
        tmp += dN_dgam_dt[i] * dt * dgam[i];
    }
    tmp -= 0.5*(dN_dgam_dt[0] * dt * dgam[0] + dN_dgam_dt[Nbin_e-1 ]* dt * dgam[Nbin_e-1]);
    *N_inj_tot += tmp;

    /// ???
    for (i=Nbin_e-1; i>0; i--){
        total_frac_depo = 0.0;
        frac_resi = 1.0;
        frac_depo = 0.0;
        for (j=i; j>0; j--){
            frac_depo = frac_resi/(1.0 + dt/dgam[i] * dgam_dt[i]);
            dN_dgam_dt_eff[j] += frac_depo * dN_dgam_dt[i] * (dgam[i]/dgam[j]); // eff?

            total_frac_depo += frac_depo;
            if(total_frac_depo > 1.0)
                break;

            frac_resi = (1.-total_frac_depo);
        }
        dN_dgam_dt_eff[0] += frac_resi*dN_dgam_dt[i]*(dgam[i]/dgam[0]);
    }
    dN_dgam_dt_eff[0] += dN_dgam_dt[0];

    for (i=0;i<Nbin_e;i++){
        dN_dgam_dt[i] = dN_dgam_dt_eff[i];
    }

}
// ------------------
class MagnetarSynchrotron{
    struct Pars{
        double gam_b{};
        double p1{};
        double p2{};
        double eps_e{};
        double gam_max{};
//        int Nbin_e{};
//        int Nbin_ph{};
    };
    std::unique_ptr<Pars> p_pars = nullptr;
    std::unique_ptr<logger> p_log = nullptr;
//    VecVector & m_data;
//    Vector spectrum{}; Vector spec_gams{}; Vector spec_freqs{};
//    Vector emissivity{}; Vector absorption{};
public:
    MagnetarSynchrotron(unsigned loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "MagnetarSynchrotron");
        p_pars = std::make_unique<Pars>();
    }
    std::unique_ptr<Pars> & getPars(){ return p_pars; }
//    Vector & getSpectrum(){return spectrum;}
//    Vector & getEmissivity(){return emissivity;}
//    Vector & getAbsorption(){return absorption;}

    void setPars(StrDbMap & pars, StrStrMap & opts){
        p_pars->gam_b = getDoublePar("gam_b",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->p1 = getDoublePar("p1",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->p2 = getDoublePar("p2",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->eps_e = getDoublePar("eps_e_spec",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->gam_max = getDoublePar("gam_max",pars,AT,p_log,1.e6, false);//pars.at("ksi_n");
        /// ---
//        spec_gams = makeVecFromString(getStrOpt("spec_gams",opts,AT,p_log,"",true),p_log);
//        spec_freqs = makeVecFromString(getStrOpt("spec_freqs",opts,AT,p_log,"",true),p_log);
//        p_pars->Nbin_e = (int)spec_gams.size();
//        p_pars->Nbin_ph = (int)spec_freqs.size();

    }

    void computeComovingSpectrum(Vector & spectrum, Vector & emissivity, Vector & absorption,
                                 Vector & times, Vector & freqs, Vector & gams, Vector & Bnb,
                                 Vector & vol, Vector & drnb, Vector & Lpsr){

        if ((spectrum.size() < 1)||(emissivity.empty() < 1)||(absorption.empty() < 1)){
            spectrum.resize(times.size()*freqs.size(), 0.);
            emissivity.resize(times.size()*freqs.size(), 0.);
            absorption.resize(times.size()*freqs.size(), 0.);
        }
        /// initialize storage
        int Nbin_e = (int)gams.size();
        double gam[Nbin_e],dN_dgam[Nbin_e],dgam[Nbin_e],dN_dgam_dt[Nbin_e],dgam_dt[Nbin_e],tad[Nbin_e],tsyn[Nbin_e];
        double N_inj_tot=0.;
        initialize_e_dis(gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,tad,tsyn,p_pars->gam_max,Nbin_e);
        int Nbin_ph = (int)freqs.size();
        double gam_ph[Nbin_ph],P_nu_syn[Nbin_ph],alpha_nu_syn[Nbin_ph];
        initialize_ph_dis(freqs,gam_ph,P_nu_syn,alpha_nu_syn);
        /// ----------------
        size_t ii = 0;
        for (size_t it = 0; it < times.size(); it++) {
//            double Bnb = 0., rnb = 0., dr = 0., Lpsr = 0., dt = 0.;
            double t = times[it];
            double dt = (times[it] - times[it - 1]);
            if (it > 1)
                calc_syn_spec(Bnb[it],drnb[it],vol[it],gam,dgam,dN_dgam,p_pars->gam_max,
                              Nbin_e,gam_ph,P_nu_syn,alpha_nu_syn,Nbin_ph);
            double number_conservation = Ntot(dgam,dN_dgam,Nbin_e)/N_inj_tot;
            if (number_conservation > 1.2 or number_conservation < 0.9){
                std::cerr << AT << " conserv="<<number_conservation<<"\n";
            }
            cooling(t,Bnb[it],dgam_dt,gam,tad,tsyn,Nbin_e);
            injection(gam, dgam, dN_dgam_dt, dgam_dt, Lpsr[it], dt, &N_inj_tot,
                      p_pars->eps_e, p_pars->gam_b, p_pars->gam_max, p_pars->p1, p_pars->p2, Nbin_e);
            time_evolution_e(dt,gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,Nbin_e);
            /// --------------------------
            if ((!std::isfinite(P_nu_syn[0]))||(!std::isfinite(P_nu_syn[Nbin_ph-1])) || (P_nu_syn[0] < 0)||(P_nu_syn[Nbin_ph-1] < 0)){
                (*p_log)(LOG_ERR,AT) << " nan in pwn spec P_nu_syn[0]=" << P_nu_syn[0] << "\n";
                std::cout << "Bnb="<<Bnb<<" gam[0]="<<gam[0]<<" dgam[0]="<<dgam[0]<<" dN_dgam[0]="<<dN_dgam[0]<<" gam_max="<<p_pars->gam_max
                          <<" Nbin_e="<<Nbin_e<<" Lpsr="<<Lpsr<<" n_ing_tot="<<N_inj_tot<<" number_conserv="<<number_conservation<<"\n";
                exit(1);
            }
//                spec[ii] = gam_ph[inu] * CGS::MeC2 / CGS::H * P_nu_syn[inu]; // [erg/s/Hz]
            for (size_t inu = 0; inu < freqs.size();  inu++) {
                emissivity[inu + freqs.size() * it] = P_nu_syn[inu]; // TODO ? IS IT TRUE?
                absorption[inu + freqs.size() * it] = alpha_nu_syn[inu];
                spectrum[inu + freqs.size() * it] = gam_ph[inu] * CGS::MeC2 / CGS::H * P_nu_syn[inu]; // [erg/s/Hz]

                double abs_fac=0.;
                double tau_sa = absorption[inu + freqs.size() * it] * drnb[it];
                if(tau_sa == 0.0)
                    abs_fac = 1.0;
                else if(tau_sa > 0.0)
                    abs_fac = -expm1(-tau_sa) / tau_sa;
                else {
                    abs_fac = expm1(tau_sa) / tau_sa; //* exp(
                    //abs * DR * beta_shock*mu / (mu - beta_shock));
                }
                spectrum[inu + freqs.size() * it] = tau_sa > 1e-6 ? emissivity[inu + freqs.size() * it] * abs_fac : emissivity[inu + freqs.size() * it];

                ii++;
            }
        }
    }

};


//    void allocateSpaceForComovingSpectrum(size_t nr){
//        /// allocate memory for comoving spectra to be evaluated
//        m_freq_grid = TOOLS::MakeLogspace(log10(p_pars->freq1),
//                                          log10(p_pars->freq2),(int)p_pars->nfreq);
//        if (p_pars->method_comp_mode == SynchrotronAnalytic::METHODS_RAD::icomovspec){
//            std::cerr << " allocating comoving spectrum array (fs) "
//                               << " freqs="<<m_freq_grid.size() << " by radii=" << nr << " Spec. grid="
//                               << m_freq_grid.size() * nr << "\n";
//            m_spectrum.resizeEachImage( m_freq_grid.size() * nr );
//        }
//    }

//class SynchrotronNumeric{
//    std::unique_ptr<logger> p_log;
//    Vector farr_sync;
//    Vector freq_ssc;
//    struct PWNPars{
//        double f1_sync, f2_sync; size_t nf_sync;
//        double f1_ssc, f2_ssc; size_t nf_ssc;
//
//    };
//    std::unique_ptr<PWNPars> p_pars;
//public:
//    SynchrotronNumeric(int loglevel){
//        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronNumeric");
//        p_pars = std::make_unique<PWNPars>();
//    }
//    void setMagPars(StrDbMap pars, StrStrMap sets){
////        p_pars->f1_ssc =
//    }
//};


//class SynchrotronAnalyticComoving{
//    enum METHODS_SYNCH { iWSPN99, iJOH06, iDER06, iMARG21 };
//    enum METHODS_LFMIN { igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };
//    enum METHODS_SSA { iSSAoff, iSSAon };
//    enum METHOD_TAU { iTHICK, iSMOOTH, iSHARP };
//    enum QQ { i_em_pl, i_em_th, i_abs_pl, i_abs_th, i_em, i_abs };
//    struct PWNPars{
//        // --- in
//        double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
//        double mu=-1, mu_e=-1;
//        bool lim_gm_to_1= true;
//        double beta_min = -1;
//        // --- methods
//        METHODS_SYNCH m_sychMethod{};
//        METHODS_LFMIN m_methodsLfmin{};
//
//        METHODS_SSA m_methods_ssa{};
//        QQ m_marg21opt_em = i_em_pl;
//        QQ m_marg21opt_abs = i_abs_pl;
//    };
//    PWNPars * p_pars = nullptr;
//
//    static constexpr size_t n_vars = 5;
//    enum QS { igm, igc };
//    std::vector<std::string> m_vars{ "gm", "gc", "num", "nuc", "pmax" };
//    VecArray m_data{};// ( n_vars );
//public:
//    SynchrotronAnalyticComoving(size_t nt){
//        p_pars = new PWNPars();
//        allocateSpace(nt);
//    }
//    ~SynchrotronAnalyticComoving(){ delete p_pars; }
//    void allocateSpace(size_t nt){
//        m_data.resizeEachImage(n_vars);
//        for (auto & arr : m_data){
//            arr.resizeEachImage( nt, 0. );
//        }
//    }
//
//    void setMagPars(std::unordered_map<std::string,double> & pars, std::unordered_map<std::string,std::string> & opts){
//        // set parameters
//
//    }
//};

#endif //SRC_RADIATION_H
