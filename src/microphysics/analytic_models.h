//
// Created by vsevolod on 1/20/24.
//

#ifndef SRC_ANALYTIC_MODELS_H
#define SRC_ANALYTIC_MODELS_H

#include "../utilitites/utils.h"
#include "../utilitites/interpolators.h"
#include "../utilitites/quadratures.h"

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
        /// calculate total_rad emissivity & optical depth:
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


class SynchrotronAnalytic{
    double gamma_min,gamma_c,gamma_max,p,B;
    bool m_methods_ssa;
    std::unique_ptr<logger> & p_log;
public:
    SynchrotronAnalytic(double gamma_min, double gamma_c, double gamma_max, double B,
                        double p,bool do_ssa,std::unique_ptr<logger> & p_log)
            :gamma_min(gamma_min),gamma_c(gamma_c),gamma_max(gamma_max),B(B),p(p),m_methods_ssa(do_ssa),
             p_log(p_log){
    }

    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchJOH06(double & em, double & abs, double nuprime, double n_prime){
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
            if (m_methods_ssa) {
                double _alpha = 7.8 * phipS * std::pow(XpS, -(4 + p) / 2.) * (p + 2) * (p - 1)
                                * CGS::qe / CGS::mp / (p + 2 / 3.);
                abs = _alpha * n_prime * CGS::mp * std::pow(gamma_min, -5) / B;
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
            if (m_methods_ssa) {
                double _alpha = 11.7 * phipF * std::pow(XpF, -3) * CGS::qe / CGS::mp;
                abs = _alpha * (n_prime * CGS::mp) * std::pow(gamma_c, -5) / B;
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
    void computeAnalyticSynchWSPN99(double & em, double & abs, double nuprime, double n_prime){
//        double gamma_min=p_pars->gamma_min, gamma_c=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
//        double p = p_pars->p;
        double scaling=0., abs_scaling=1.;
//        double em=0.,abs=0.;
//            double PhiP = interpolated_phi(p);
//            emissivity = (interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi)
//                    * n_prime * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);

        double emissivity = (interpolated_phi(p) * sqrt(3.0) / 4.0 * CGS::pi)
                          * n_prime * CGS::qe * CGS::qe * CGS::qe * B / (CGS::me * CGS::c * CGS::c);
        double Xp = interpolated_xi(p);
        double nu_m = 3.0 / ( 4.0 * CGS::pi ) * Xp * gamma_min * gamma_min * CGS::qe * B / (CGS::me * CGS::c );
        double nu_c = 0.286 * 3. * gamma_c * gamma_c * CGS::qe * B / (4.0 * CGS::pi * CGS::me * CGS::c );
        if (nu_m <= nu_c){//  # slow cooling
            if (nuprime < nu_m) {
                scaling = std::pow(nuprime / nu_m, 1.0 / 3.0);
            }
            else if (nuprime >= nu_m && nuprime < nu_c) {
                scaling = std::pow(nuprime / nu_m, -1.0 * (p - 1.0) / 2.0);
            }
            else  { // if (nuprime >= nu_c)
                scaling = std::pow(nu_c / nu_m, -1.0 * (p - 1.0) / 2.0)
                        * std::pow(nuprime / nu_c, -1.0 * p / 2.0);
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
        em = emissivity * scaling;


        /// from vanEarten+2010
        if (m_methods_ssa) {
            double absorption = sqrt(3) * std::pow(CGS::qe, 3) * (p - 1) * (p + 2) * n_prime * B
                  / (16 * M_PI * CGS::me * CGS::me * CGS::c * CGS::c * gamma_min * nuprime * nuprime);
            if (nuprime < nu_m) // slow cooling
                abs_scaling = std::pow(nuprime / nu_m, 1.0 / 3.0);
            else
                abs_scaling = std::pow(nuprime / nu_m, -0.5 * p);
            abs = absorption * abs_scaling;

        }
//            std::cout << 11.17 * (p - 1.0) / (3.0 * p - 1.0) * (0.54 + 0.08 * p) << " " << interpolated_phi(p) * sqrt(3.0) /  4.0 * CGS::pi << "\n";
//            exit(1);
    }
    /// Analytical Synchrotron Sectrum; BPL;
    void computeAnalyticSynchDER06(double & em, double & abs, double nuprime, double n_prime){

        int numbins = 200;
        double gammmas[numbins];

        double step = std::exp((std::log(gamma_max) / (double)200));
        for (size_t i = 0; i < numbins; i++)
            gammmas[i] = std::pow(step,(double)i);

        double epsilon = nuprime * CGS::h / CGS::mec2;
        auto integrand_ele = [&](const double gam){
            double p1 = gamma_min < gamma_c ? p : 2.0;
            double p2 = p + 1.0;
            double gmax = gamma_max;
            double gmin = gamma_min < gamma_c ? gamma_min : gamma_c;
            double gb = gamma_min < gamma_c ? gamma_c : gamma_min;

            return Dermer09::brokenPowerLaw(gam, gmin, gb, gmax, p1, p2);
        };

        double k_e_s[numbins];
        double k_e = 0., power_e = 0., power = 0.;
        for (size_t i = 0; i < numbins-1; i++) {
            k_e_s[i] = integrand_ele(gammmas[i]);
            k_e += k_e_s[i] * (gammmas[i + 1] - gammmas[i]);
        }
        k_e = n_prime / k_e;

        for (size_t i = 0; i < numbins-1; i++) {
            power_e = Dermer09::single_electron_synch_power( B, epsilon, gammmas[i] );
            power += k_e_s[i] * (gammmas[i + 1] - gammmas[i]) * power_e;
        }
        power *= k_e;


        em = power * (CGS::h / CGS::mec2);

        if (m_methods_ssa) {
            auto integrand_ele_ssa = [&](double gam) {
                double p1 = gamma_min < gamma_c ? p : 2.0;
                double p2 = p + 1.0;
                double gmax = gamma_max;
                double gmin = gamma_min < gamma_c ? gamma_min : gamma_c;
                double gb = gamma_min < gamma_c ? gamma_c : gamma_min;
                return Dermer09::brokenPowerLawSSA(gam, gmin, gb, gmax, p1, p2);
            };

            double k_e_s_ssa[numbins];

            double k_e_ssa = 0., absorption = 0.;
            for (size_t i = 0; i < numbins-1; i++) {
                k_e_s_ssa[i] = integrand_ele_ssa(gammmas[i]);
                k_e_ssa += k_e_s_ssa[i] * (gammmas[i + 1] - gammmas[i]);
            }
            k_e_ssa = n_prime / k_e_ssa;

            for (size_t i = 0; i < numbins-1; i++) {
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
    void computeAnalyticSynchMARG21(double & em, double & abs, double nuprime, double ne_,
                                    double delta, double Theta, double z_cool, double accel_frac){
//        double gamma_min=p_pars->gamma_min, gamma_c=p_pars->gamma_c, B=p_pars->B, n_prime=p_pars->n_prime;
//        double Theta=p_pars->Theta, z_cool = p_pars->z_cool, acc_frac = p_pars->accel_frac;
//        double p = p_pars->p;
        double x, em_th=0., em_pl=0., abs_th=0., abs_pl=0.; // for Margalit model
        /// Margalit+21 arXiv:2111.00012
        ///
//        double ne = n_prime * ksi_n;
//        double ne_ = mu_e * n_prime;//4.0 * p_pars->mu_e * p_pars->n_ism;
//        double delta = eps_e / eps_t; // defined in Margalit+21 arXiv:2111.00012
        /// normalized frequency:
        x = nuprime / Margalit21::nu_Theta(Theta, B);
        if (!std::isfinite(x)){
            (*p_log)(LOG_ERR,AT) << " x = "<< x << "\n";
            exit(1);
        }
        /// calculate total_rad emissivity & optical depth:
        em_th = Margalit21::jnu_th(x, ne_, B, Theta, z_cool);
        em_pl = Margalit21::jnu_pl(x, ne_ * accel_frac, B, Theta, gamma_min, delta, p, z_cool); //TODO added 'accel_frac'
        if ((!std::isfinite(em_th))||(em_th < 0.)) em_th = 0.;
        if ((!std::isfinite(em_pl))||(em_pl < 0.)) em_pl = 0.;
        if ((em_pl==0.) && (em_th==0.)){
            (*p_log)(LOG_ERR,AT) << " em_pl=em_th=0 for"
                                 << " n_prime=" << ne_ << " accel_frac" << accel_frac
                                 << " B=" << B << " gamma_min=" << gamma_min << " Theta=" << Theta
                                 << " z_cool=" << z_cool << " nuprime=" << nuprime << " x=" << x
                                 << "\n";
            exit(1);
        }
//            emissivity = em_pl + em_th;
        em=em_pl+em_th;//m_data[i_em] = em_pl + em_th;

        if (m_methods_ssa) {
            abs_th = Margalit21::alphanu_th(x, ne_, B, Theta, z_cool);
            abs_pl = Margalit21::alphanu_pl(x, ne_ * accel_frac, B, Theta, gamma_min, delta, p, z_cool);
            if (!std::isfinite(abs_th)) abs_th = 0.;
            if (!std::isfinite(abs_pl)) abs_pl = 0.;
//                abs = abs_th + abs_pl;
            abs=abs_pl+abs_th;//m_data[i_abs] = abs_th + abs_pl;

        }
    }
};

#endif //SRC_ANALYTIC_MODELS_H
