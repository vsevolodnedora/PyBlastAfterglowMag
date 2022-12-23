//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BASE_EQUATIONS_H
#define SRC_BASE_EQUATIONS_H

#include "pch.h"
#include "utils.h"
#include "interpolators.h"
#include "quadratures.h"

/* ------------- EQUATIONS ----------------- */

static size_t findIndex( const double & x, const Array & arr, size_t N ) {
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N-2;

    unsigned int i = ((unsigned int) N) >> 1;
    unsigned int a = 0;
    unsigned int b = N-1;

    // https://stackoverflow.com/questions/4192440/is-there-any-difference-between-1u-and-1-in-c/4192469
    // untill the m_size of b-a > 1 continue shrinking the array, approaching the 'x'
    while (b-a > 1u) // 1U is an unsigned value with the single bit 0 set
    {
        i = (b+a) >> 1; // ???
        if (arr[i] > x)
            b = i;
        else
            a = i;
    }

    return (int)a;
}
static inline double interpSegLin( size_t & a, size_t & b, const double & x, Array & X, Array & Y) {
    // take two indexis, 'a' 'b' and two arrays 'X' and 'Y' and interpolate
    // between them 'Y' for at 'x" of "X'
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];
    return ya + (yb-ya) * (x-xa)/(xb-xa);
}
static inline double interpSegLog( size_t & a, size_t & b, double x, Array & X, Array & Y) {
//        std::cout << a << ' ' << b << ' '<< x << ' '<< X << ' '<< Y << ' '<< N << "\n";
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];

    return ya * std::pow(yb/ya, log(x/xa)/log(xb/xa));
}

inline namespace EQS{
    /*
 * Compute velocity in [c] from Lorentz Factor 'iGamma'
 */
    double Beta(double const &Gamma){
//        return sqrt(1.0 - pow(Gamma, -2));
        return (1. / Gamma) * sqrt( (Gamma - 1.) * (Gamma + 1.) );
    }
    double Beta2(double const &Gamma){
        return (Gamma-1.) * (Gamma+1.) / ( Gamma * Gamma );
    }
    Array Beta(Array &Gamma){
//        return sqrt(1.0 - (1.0 / (Gamma * Gamma)) );
        return (1. / Gamma) * sqrt( (Gamma - 1.) * (Gamma + 1.) );
    }

    double Gamma(const double & beta){
//        return sqrt(1.0 / (1.0 - (beta * beta)));
        return sqrt(1. + (beta * beta / (1. - beta * beta)));
    }
    Vector Gamma(Vector & beta){
//        return sqrt(1.0 / (1.0 - (beta * beta)));
        Vector res ( beta.size() );
        for(size_t i = 0; i < beta.size(); i++)
            res[i] = Gamma(beta[i]);
        return std::move( res );
    }

    double GammaRel(const double Gamma1, const double Gamma2){
//        return Gamma1 * Gamma2 * (1. - EQS::Beta(Gamma1) * EQS::Beta(Gamma2));
        return Gamma1 * Gamma2 - sqrt(Gamma1*Gamma1-1.) * sqrt(Gamma2*Gamma2-1.);
    }
    double BetaRel(const double beta1, const double beta2){
        return (1. - beta1 * beta2) * sqrt(1./(1.-beta1*beta1)) * sqrt(1./(1.-beta2*beta2));
    }

    /*
     * Compute fluid internal energy from adiabatic index (See Peer+2012)
     */
    double InternalEnergy(double const &Gamma, double const &rho, double(*AdiabIndexFunc)(double const &)){
        double ada = AdiabIndexFunc(Gamma);
        return (ada * Gamma + 1.) / (ada - 1) * (Gamma - 1.) * rho * CGS::c * CGS::c;
    }
    double InternalEnergy(double const &Gamma, double const &rho, double gammaAdi){
        double ada = gammaAdi;
        return (ada * Gamma + 1.) / (ada - 1) * (Gamma - 1.) * rho * CGS::c * CGS::c;
    }

    /*
     * compute the swept-up mass by the expanding blast wave
     * Accounts for lateral spreading if 'aa' < 0
     * https://arxiv.org/pdf/1203.5797.pdf
     * https://academic.oup.com/mnras/article/421/1/570/990386 (for introduction of 'a' parameter)
     *
     * basic form was dmdr = 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / m_pars["m_scale"]
     *
     */
    double dmdr(double &Gamma, double &RR, double &thetaE, double &theta, double &rho, double aa){

        // First term: change in swept-up mass due to the change in solid angle
        // double t1 = 0. if (aa < 0) else (1. / 3.) * np.sin(itheta) / (iGamma ** (1. + aa) * itheta ** (aa));
        double t1 = (aa < 0) ? 0.0 : sin(theta) / (3.0 * std::pow(Gamma, 1.0 + aa) * std::pow(theta, aa));
        // Second term: change in swept-up mass due to radial expansion
        double t2 = (cos(thetaE) - cos(theta));
        return 2.0 * CGS::pi * rho * (t1 + t2) * RR * RR;
    }
    double dmdr(double &Gamma, double &RR, double &dthetadR, double &theta, double &rho){
        return 2.*CGS::pi*rho*((1.- cos(theta)) + (1./3.)* sin(theta)*RR*dthetadR)*RR*RR;
    }

    /*
     * Evolution equation for the lorentz factor (see Peer+2012)
     */
    double dgdr(double const &M0, double const &Gamma, double const &beta,
                double const &mm, double const &gamma_adi, const double dmdr){
        double numerator = -(gamma_adi * (Gamma * Gamma - 1.0) - (gamma_adi - 1.0) * Gamma * beta * beta);
        double denominator = M0 + mm * (2.0 * gamma_adi * Gamma - (gamma_adi - 1) * (1. + std::pow(Gamma, -2)));
        return (numerator / denominator) * dmdr;
    }

    /*
     * Compute mass of the shall assuming relativistic, radiation dominated shell
     * M = E * c^-2 / iGamma
     */
    double EtoM(const double &E, const double &Gamma){
        return E * std::pow(CGS::c, -2) / Gamma;
    }




    /*
     * Radius at which blast wave starts to decelerate
     * Assumes constant density ISM density
     * Source: #TODO find the source
     */
    double Rdec2(const double &E, const double &nn, const double &Gamma){
        return std::pow(3. / (4. * CGS::pi) * 1. / (CGS::c * CGS::c * CGS::mp) * E / (nn * Gamma * Gamma), 1.0 / 3.0);
    }



    /*
     * Compute time in the lab frame
     */
    double init_elapsed_time(const double &r0, const double &Gamma0, const bool &use_spread){
        double beta0 = Beta(Gamma0);
        double result;
        if (use_spread){
            double dThetadr0 = 0.0;
            result= r0 / (CGS::c * beta0) * (sqrt(1. + r0 * r0 * dThetadr0 * dThetadr0)) - r0 / CGS::c;
        } else {
            result= r0 / (CGS::c * Gamma0 * Gamma0 * beta0 * (1.0 + beta0));
        }
        if (result!=result){
            std::cerr<<AT<<" Failed to init TT"
                     <<" i="<<0
                     <<" Rs[0]"<<r0
                     <<" Gamma[0]"<<Gamma0
                     <<" use_spread="<<use_spread
                     <<" result[0]="<<result
                     <<" \n";
            exit(1);
        }
        return result;
    }

    /*
     * Integrate the elapsed time in the lab frame
     */
    double integrate_elapsed_time(unsigned int i,
                                  std::valarray<double> &Rss,
                                  std::valarray<double> &Gammas,
                                  std::valarray<double> &thetas,
                                  const bool &use_spread){

        // integrate only to 'i'

        size_t n = i+1;
        double result = 0.0;
        Array Rs(n);
        for (size_t j = 0; j < n; j++) Rs[j] = Rss[j];
        Array integrand( n );
        Array dthetadr( i+1 );
        if (use_spread){
            // considering lateral spreading
            dthetadr[0] = 0.0;
            for (unsigned int j = 1; j < n; j++) {
                double beta = Beta(Gammas[j]);
                dthetadr[j] = (thetas[j] - thetas[j - 1]) / (Rs[j] - Rs[j - 1]);
                integrand[j] = 1.0 / (CGS::c * beta) * sqrt(1. + Rs[j] * Rs[j] * dthetadr[j] * dthetadr[j]) - 1. / (CGS::c);
            }
            auto func = [&](double x){ return ( Interp1d(Rs,integrand)).Interpolate(x, Interp1d::iLagrangeLinear); };
            result = Simpson38(Rs[0], Rs[i], 100, func);
        } else {
            // No spreading
            for (unsigned int j = 1; j < n; j++) {
                double beta = Beta(Gammas[j]);
                integrand[j] = 1.0 / (CGS::c * Gammas[j] * Gammas[j] * beta * (1.0 + beta));
            }
            auto func = [&](double x){ return ( Interp1d(Rs,integrand)).Interpolate(x, Interp1d::iLagrangeLinear); };
            result = Simpson38(Rs[0], Rs[i], 100, func);
        }
        if (result!=result){
            std::cerr<< AT
                     <<" Failed to compute TT"
                     <<" i="<<i
                     <<" Rs[0]"<<Rs[0]
                     <<" Rs[i]"<<Rs[i]
                     <<" Gamma[0]"<<Gammas[0]
                     <<" Gamma[i]"<<Gammas[i]
                     <<" thetas[0]"<<thetas[0]
                     <<" thetas[i]"<<thetas[i]
                     <<" use_spread="<<use_spread
                     <<" integrand[0]="<<integrand[0]
                     <<" integrand[i]="<<integrand[i]
                     <<" \n";
            exit(1);
        }
        return result;
    }

    /*
     * Compute the thickness of the shock created by the blast wave
     * See Johanesson+2006
     *
     */
    double shock_thickness(const double &mass,
                           const double &rhoprime,
                           const double &theta_a,
                           const double &theta_b,
                           const double &Gamma,
                           const double &R,
                           const double &ncells){
        double one_min_costheta = (cos(theta_a) - cos(theta_b) / ncells); // TODO :: THIS equation might be incorrect
        // rhoprime = 4. * rho * iGamma // for relativistic shock
        return mass / (2.0 * CGS::pi * one_min_costheta * rhoprime * Gamma * R * R);
    }
//    double shock_thickness(const double m2or3, const double theta_a,
//                           const double theta, const double Gamma,
//                           const double rho_or_rho4, const double R, size_t ncells){
//        double rhoprime = 4.0 * rho_or_rho4 * Gamma;
//        double one_min_costheta = (cos(theta_a) - cos(theta) / ncells); // TODO :: THIS equation might be incorrect
//        return m2or3 / (2 * CGS::pi * one_min_costheta * rhoprime * Gamma * R * R)
//    }

    /// simple argument for shock thickness (delta R) See e.g., vanEerten+2010
    double shock_delta(const double & R, const double & Gamma) {
        return R / (12. * Gamma * Gamma);
    }

    /// shock thickness from the considerations of non-uniform ISM (see Johannesson+06)
    double shock_delta_joh06(const double &R, const double &M2, const double &theta,
                             const double &Gamma, const double &rho2, const double &ncells) {
        double nprime = rho2 / CGS::mp;
        return M2 / (2. * CGS::pi * CGS::mp * R * R * (1. - cos(theta) / ncells) * Gamma * nprime);
    }

    /// Using fixed compression ration post/pre-shocked material (4.)
    double rho2(const double & Gamma, const double & rho) {
        return 4. * Gamma * rho;
    }
    /// Using 'proper' compression ratio that depends on the fluid adiabatic index 'adi'
    double rho2t(const double & Gamma, const double & adi, const double & rho){
        return (adi * Gamma + 1.) / (adi - 1.) * rho;
    }

    /// From Rankie-Hugonoid Jump Conditions (See Zhang Book)
    double GammaSh(const double & Gamma, const double & adi){
        return (Gamma+1.)*(adi*(Gamma-1.)+1.) / ( adi * (2.-adi)*(Gamma-1.) + 2. );
    }

    // *************************** NAVA et al 2013 *******************

    inline double get_GammaEff(const double &Gamma, const double &gammaAdi){
        return (gammaAdi * Gamma * Gamma - gammaAdi + 1.0) / Gamma;
    }

    inline double get_dGammaEffdGamma(const double &Gamma, const double &gammaAdi){
        return (gammaAdi * Gamma * Gamma + gammaAdi - 1.0) / (Gamma * Gamma);
    }

    double dGammadR_fs(const double &Gamma,
                       const double &gammaAdi,
                       const double &dlnrho1dR,
                       const double &M2,
                       const double &dM2dR,
                       const double &Eint2){

        double GammaEff = get_GammaEff(Gamma, gammaAdi); // gammaAdi) #(gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
        double dGammaEffdGamma = get_dGammaEffdGamma(Gamma, gammaAdi); // 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.

        double f_2 = GammaEff * (gammaAdi - 1.) * Eint2 / Gamma;
        double h_2 = GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR / M2 - dlnrho1dR);

        double dGammadR = -((Gamma - 1.) * (GammaEff + 1.) * dM2dR - h_2) / ((1. + M2) + Eint2 * dGammaEffdGamma + f_2);

        return dGammadR;
    }


    double get_U_p(const double &rhoprim, //const double &Gamma, // GammaSh
                   const double &M2, const double &Eint2){
//        double rhoprim = 4. * rho * Gamma ;     // comoving density
        double V2 = M2 / rhoprim ;              // comoving volume
        double U_p = Eint2 / V2 ;               // comoving energy density (electrons)
        // U_b = eps_b * U_e                    // comoving energy density (MF)
        return U_p;
    }
    double get_U_p(const double & rhoprim, const double & Gamma){
        return (Gamma - 1.) * rhoprim * CGS::c * CGS::c;
    }

    // *************************** NAVA et al 2013 (Reverse shock) *******************

    /*
     * Adiabatic index in the 3rd region (See Nava+2013)
     */
    inline double get_gamma43_minus_one(const double &Gamma,
                                        const double &Gamma0,
                                        const double &beta,
                                        const double &beta0,
                                        const double &beta_switch=0.999995){
        /*
         * Idea:
         *
         * gamma43_minus_one = Gamma * Gamma0 * (1 - beta * beta0)
         * However, this is numerically... difficult, as 1-beta*beta0 is subjected to truncation error a lot
         * Using expansions we write
         *
         * gamma43_minus_one = Gamma*Gamma0 - np.sqrt(Gamma0**2 - 1) * np.sqrt(Gamma**2 - 1) # -- wolfram. Creats jump
         * if beta < 0.999:
         * gamma43_minus_one = Gamma * Gamma0 * (1 - beta * beta0)
         *
         */
        double gamma43_minus_one;
        if (Gamma > beta_switch * Gamma0){
            gamma43_minus_one = 1.0 - 1.0;
        } else {
            gamma43_minus_one = Gamma * Gamma0 * \
                            (1.0 / (Gamma0 * Gamma0) + 1.0 / (Gamma * Gamma) - 1.0 / (Gamma * Gamma) / (Gamma0 * Gamma0)) / \
                            (1 + beta * beta0) - 1.0;
        }
        return gamma43_minus_one;
    }

    inline double get_dgamma43dGamma(const double &Gamma0, const double &Gamma){
        /*
         * Idea:
         * Gamma0 - Gamma * (np.sqrt(Gamma0**2-1)/np.sqrt(Gamma**2-1))
         * But numerically unstable -- truncation error
         */
        return Gamma0 - (Gamma0*Gamma0 * Gamma - Gamma) /
                        sqrt(Gamma*Gamma * Gamma0*Gamma0 - Gamma0*Gamma0 - Gamma*Gamma + 1.0);

    }

    inline double get_dGammaEff3dGamma(const double &Gamma,
                                       const double &gammaAdi3,
                                       const double &dgamma43dGamma,
                                       const double &gamma43){
        return gammaAdi3 * (1.0 + 1.0 / (Gamma * Gamma)) - dgamma43dGamma / 3.0 / gamma43*gamma43 *
                                                           (Gamma - 1.0 / Gamma) - 1.0/(Gamma * Gamma);

    }

    /*
     * See Nava+2013 for the details (appendix on the reverse shock implementation)
     */
    double get_dGammadR_fs_rs(const double &Gamma, const double &Gamma0, const double &gammaAdi,
                              const double &dlnrho1dR, const double &M2, const double &dM2dR,
                              const double &dlnrho4dR, const double &M3, const double &dM3dR,
                              const double &Eint2, const double &Eint3, const double &gammaAdi3,
                              const double &gamma43_m1){

        // Using asymptotic approximation of beta
        // 0.5 / Gamma0 + 0.5 / Gamma ** 2 * (0.5 / Gamma0 - Gamma0 - 3. * Gamma0 / 8. / Gamma ** 2)
        long double dgamma43dGamma = get_dgamma43dGamma(Gamma0, Gamma);

        // (gammaAdi3 * Gamma ** 2 - gammaAdi3 + 1) / Gamma # (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
        long double GammaEff3 = get_GammaEff(Gamma, gammaAdi3);

        // gammaAdi3 * (1. + Gamma ** -2) - dgamma43dGamma / 3. / (gamma43_m1 + 1.) ** 2 * (Gamma - 1. / Gamma) - Gamma ** -2
        long double dGammaEff3dGamma = get_dGammaEff3dGamma(Gamma, gammaAdi3, dgamma43dGamma, gamma43_m1 + 1.0);

        // 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.
        double dGammaEffdGamma = get_dGammaEffdGamma(Gamma, gammaAdi);
        double GammaEff = get_GammaEff(Gamma, gammaAdi);  // (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

        double f_2 = GammaEff * (gammaAdi - 1.0) * Eint2 / Gamma;
        double h_2 = GammaEff * (gammaAdi - 1.0) * Eint2 * (dM2dR / M2 - dlnrho1dR);

        double fh_factor3 = GammaEff3 * (gammaAdi3 - 1.0) * Eint3;
        double f_3 = fh_factor3 / (gamma43_m1 + 1.0) * dgamma43dGamma;
        double h_3 = 0.0;
        if (Eint3 != 0)
            h_3 = fh_factor3 * (dM3dR / M3 - dlnrho4dR);

        double dGammadR = -((Gamma - 1.0) * (GammaEff + 1.0) * dM2dR + (Gamma - Gamma0 + GammaEff3 * gamma43_m1) *
                                                                       dM3dR - h_2 - h_3) / ((M2 + M3) + Eint2 * dGammaEffdGamma + Eint3 * dGammaEff3dGamma + f_2 + f_3);

        return dGammadR;

    }
    /*
     * Density in the unshocked region of the ejecta (region 4 in Nava+2013)
     */
    double rho4( const double &R, const double &deltaR4, const double &tprompt, const double &beta0,
                 const double &M0, const double &theta_a, const double &theta_b0){
        double alpha_of = tprompt * beta0 * CGS::c;
        double rho4_fac_1 = M0 / (2.0 * alpha_of * CGS::pi * (cos(theta_a) - cos(theta_b0)));
        double rho4_fac = rho4_fac_1 / (R * R);
        double rho4 = rho4_fac * exp(-deltaR4 / alpha_of); // assuming exponential cutoff
        return rho4 ;
    }

    double get_U_e_rs(const double &Gamma, const double &M3, const double &Eint3, const double &rho4){
        double rho3 = 4.0 * Gamma * rho4; // TODO is it correct iwth Gamma and not Gamma43?
        double V3 = M3 / rho3;            // comoving volume
        double U_e = Eint3 / V3;          // comoving energy density
        return U_e;
    }


}

/// Blandford-McKee self-simialr solution
class BlandfordMcKee2{
    struct Pars{
        double Etot = -1;
        double k = 0;
        double Rscale = -1;
        double gamAdi = 4/3.;
        double eta     = 1.e-5;        //          eta = p/(rho*c^2)
        double n_ext = 1.e0;      // external medium number density

        double n0      = 1.;           // cm-3:    CBM number density
        double rho0    = n0*CGS::mp;       // g.cm-3:  comoving CBM mass density
        double rhoNorm = rho0;     // n0*mp_;       // g.cm-3:  comoving CBM mass density
        double lNorm = CGS::c;                     // distance normalised to c
        double vNorm = CGS::c;                     // velocity normalised to c
        double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

        double rho_ = 0;
        double Gamma_ = 0;
    };
    Pars * p_pars = nullptr;
public:
    BlandfordMcKee2(){
        p_pars = new Pars;
    }
    ~BlandfordMcKee2(){ delete p_pars; }
    void calcBM(double r, double t, double & rho, double & u, double & p){

        auto & _p = * p_pars;

        // returns normalised primitive variables for BM solution at position r (denormalised).

        // repeating from above (these are local)
        double lfacShock2 = _p.Etot * (17.- 4.*_p.k)/(8.*CGS::pi*_p.n_ext*CGS::mp*std::pow(t,3)*std::pow(CGS::c,5));
        // Blandford&McKee(1976) eq. 69
        double RShock = CGS::c * t * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2));
        // Blandford&McKee(1976) eq. 27
        double rhoa = _p.n_ext * CGS::mp * std::pow(RShock/_p.Rscale, -_p.k);

        double pa = rhoa * _p.eta * CGS::c * CGS::c;
        double ha = CGS::c * CGS::c + pa * _p.gamAdi / (_p.gamAdi - 1.) / rhoa;
        double pf = 2./3. * lfacShock2 * ha * rhoa;
        double lfacf = sqrt(fmax(1.,.5 * lfacShock2));
        double Df = 2. * lfacShock2 * rhoa;

        // Blandford&McKee(1976) eq. 8-10
        double chi = (1. + 2. * (4. - _p.k) * lfacShock2) * ( 1. - r / ( CGS::c * t) );
        p = pf*std::pow(chi, -(17.-4.*_p.k)/(12.-3.*_p.k)) / _p.pNorm;
        double lfac = sqrt(lfacf * lfacf/chi + 1);  // (+1) to ensure lfac>1
        double D = Df * std::pow(chi, -(7.-2.*_p.k)/(4.-_p.k));
        // Blandford&McKee(1976) eq. 28-30 / 65-67
        rho = D / lfac / _p.rhoNorm;
        u = CGS::c * sqrt(1.-1./(lfac*lfac))*lfac / _p.vNorm;
    }

    void eval(double R, double Rshock, double Gamma, double n_ism, double rho2){
        auto & _p = * p_pars;
        _p.n0 = n_ism ;
        _p.n_ext = n_ism ;

        double lfacShock2 = Gamma * Gamma;
//        double RShock = CGS::c * t * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2));
        double t = Rshock / ( CGS::c * (1. - 1./(1. + 2. * (4. - _p.k) * lfacShock2)) );
        double rhoa = _p.n_ext * CGS::mp;// * pow(RShock/_p.Rscale, -_p.k);
        double Df = 2. * lfacShock2 * rhoa;
        double lfacf = sqrt(fmax(1.,.5 * lfacShock2));
        double chi = (1. + 2. * (4. - _p.k) * lfacShock2) * ( 1. - R / ( CGS::c * t) );
        double lfac = sqrt(lfacf * lfacf / chi + 1);  // (+1) to ensure lfac>1
        double D = Df * std::pow(chi, -(7. - 2. * _p.k)/(4. - _p.k));
        _p.rho_ = D / lfac / _p.rhoNorm;
        _p.Gamma_ = lfac;//CGS::c * sqrt(1.-1./(lfac*lfac)) / p_pars->vNorm;
        if(_p.rho_ > rho2 / CGS::mp){
            std::cerr << " rho behind shock(" << _p.rho_ << ") > rho2(" << rho2 << ") [BM]\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        if(_p.Gamma_ > Gamma){
            std::cerr << " Gamma behind shock(" << _p.Gamma_ << ") > Gamma(" << Gamma << ") [BM]\n";
            std::cerr << AT << "\n";
            exit(1);
        }

    }

    double rho_downstream(double R, double Rshock, double Gamma, double n_ism, double rho2){
        eval(R,Rshock,Gamma,n_ism,rho2);
        return p_pars->rho_;
    }

    double gamma_downstream(double R, double Rshock, double Gamma, double n_ism, double rho2){
        eval(R,Rshock,Gamma,n_ism,rho2);
        return p_pars->Gamma_;
    }

};

/// Sedov-Taylor self-similar solution
class SedovTaylor{
    struct Pars{
        double gamma = -1;
        size_t nu = 0;
        double w = -1;
        // ------------
        double w1{}, w2{}, w3{};
        double b0{}, b1{}, b2{}, b3{}, b4{}, b5{}, b6{}, b7{}, b8{};
        double c0{}, c1{}, c2{}, c3{}, c4{}, c5{}, c6{}, c7{}, c8{};
        Array f{};//, 1e5);
        Array eta{};
        Array d{};
        Array p{};
        Array vv{};
    };
    Pars * p_pars{};
public:
    SedovTaylor() {
        p_pars = new Pars();
        p_pars->f.resize(1e5);
        p_pars->eta.resize( p_pars->f.size() );
    }
    ~SedovTaylor() { delete p_pars; }
    void setPars(double gamma, size_t nu, double w){
        p_pars->gamma = gamma;
        p_pars->nu = nu;
        p_pars->w = w;
    }
    void evaluate(){

//        std::cout << AT << " \n Computing Sedov-Taylor profile for gamma="
//                  << p_pars->gamma<<" and grid of " << p_pars->eta.size() << " \n";

        auto nu = (double)p_pars->nu;
        double gamma = p_pars->gamma;
        double w = p_pars->w;
        // Constants for the parametric equations:
        p_pars->w1 = (3 * nu - 2 + gamma * (2 - nu)) / (gamma + 1.);
        p_pars->w2 = (2. * (gamma - 1) + nu) / gamma;
        p_pars->w3 = nu * (2. - gamma);

        p_pars->b0 = 1. / (nu * gamma - nu + 2);
        p_pars->b2 = (gamma - 1.) / (gamma * (p_pars->w2 - w));
        p_pars->b3 = (nu - w) / (float(gamma) * (p_pars->w2 - w));
        p_pars->b5 = (2. * nu - w * (gamma + 1)) / (p_pars->w3 - w);
        p_pars->b6 = 2. / (nu + 2 - w);
        p_pars->b1 = p_pars->b2 + (gamma + 1.) * p_pars->b0 - p_pars->b6;
        p_pars->b4 = p_pars->b1 * (nu - w) * (nu + 2. - w) / (p_pars->w3 - w);
        p_pars->b7 = w * p_pars->b6;
        p_pars->b8 = nu * p_pars->b6;

        // simple interpolation of correct function (only for nu=1,2,3)
        p_pars->c0 = 2 * (nu - 1) * CGS::pi + (nu - 2) * (nu - 3);
        p_pars->c5 = 2. / (gamma - 1);
        p_pars->c6 = (gamma + 1) / 2.;
        p_pars->c1 = p_pars->c5 * gamma;
        p_pars->c2 = p_pars->c6 / gamma;
        p_pars->c3 = (nu * gamma - nu + 2.) / ((p_pars->w1 - w) * p_pars->c6);
        p_pars->c4 = (nu + 2. - w) * p_pars->b0 * p_pars->c6;

        // Characterize the solution
        double f_min =  p_pars->w1 > w ? p_pars->c2 : p_pars->c6;

        p_pars->f = TOOLS::MakeLogspace(log10(f_min), 0., (int)p_pars->f.size()); // log10(f_min)
//        print_x_as_numpy(p_pars->f, 100, "f");


//                p_pars->c2 ? p_pars->w1 > w : p_pars->c6
        bool use_uniform_log_grid = false;
        if (use_uniform_log_grid){
            p_pars->f = TOOLS::MakeLogspace(log10(f_min), 0., (int)p_pars->f.size());
        }
        else {
            double tmp_fmin = log10(f_min);
            double tmp_fmax = tmp_fmin;
            std::vector<size_t> segment_lengths = {10000, 1000, 100, 100};
            VecArray segments{};
            size_t tot_size = 0;
            for (size_t i_segment = 0; i_segment < segment_lengths.size() - 1; i_segment++) {
                double _p = (double) i_segment - (double) segment_lengths.size() - 1;
                double step = std::pow(10, (double) _p);
                if (tmp_fmax == 1) tmp_fmax = 0;
                //            std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<tmp_fmin + step<<"\n";
                segments.emplace_back(TOOLS::MakeLogspace(tmp_fmin, tmp_fmin + step, (int) segment_lengths[i_segment]));
                //            std::cout << segments[i_segment] << "\n";
                tmp_fmin = tmp_fmin + step;
                tot_size = tot_size + segment_lengths[i_segment];
            }
            segments.emplace_back(TOOLS::MakeLogspace(tmp_fmin, 0, (int) segment_lengths[segment_lengths.size() - 1]));
//            std::cout << segments[segments.size() - 1] << "\n";
            tot_size = tot_size + segment_lengths[segment_lengths.size() - 1];
            //        std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<0<<"\n";
            Array tmp(tot_size);
            size_t ii = 0;
            for (size_t i_segment = 0; i_segment < segments.size(); i_segment++) {
                for (size_t ip = 0; ip < segment_lengths[i_segment]; ip++) {
                    tmp[ii] = segments[i_segment][ip];
                    ii++;
                }
            }
//            std::cout << tmp << '\n';
            p_pars->f = tmp; // log10(f_min)
        }


        // Sort the etas for our interpolation function
        p_pars->eta = parametrized_eta(p_pars->f);
        std::valarray<size_t> idx = sort_indexes(p_pars->eta);
        p_pars->f = sort_by_indexes(p_pars->f, idx);
        p_pars->eta = sort_by_indexes(p_pars->eta, idx);

        p_pars->d = parametrized_d(p_pars->f);
        p_pars->p = parametrized_p(p_pars->f);
        p_pars->vv = parametrized_v(p_pars->f);

//        print_x_as_numpy(p_pars->eta,100,"eta");
//        print_x_as_numpy(p_pars->d,100,"d");

        if (p_pars->eta[0] > 0.){
            std::cerr << " in Sedov Taylow, the extend to eta is too short. Addd more 0.00000..." << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        // Finally Calculate the normalization of R_s:
//        auto tmp = std::pow(p_pars->vv, 2);
        Array integral_ = std::pow(p_pars->eta, nu - 1) * (p_pars->d * std::pow(p_pars->vv, 2.) + p_pars->p);
        Array integral( integral_.size() - 1 );
        Array deta (integral_.size() - 1);
        double integ = 0;
        for (size_t i = 0; i < integral.size(); i++){
            integral[i] = 0.5 * (integral_[i+1] + integral_[i]);
            deta[i] = (p_pars->eta[i+1] - p_pars->eta[i]);
            integ += integral[i] * deta[i];
        }
        double alpha = integ * (8 * p_pars->c0) / ((std::pow(gamma, 2) - 1.) * std::pow(nu + 2. - w, 2));
//        p_pars->_c = pow(1. / alpha, 1. / (nu + 2 - w))


    }

    Array parametrized_eta(Array & var) {
        return std::pow(var, -p_pars->b6) * std::pow((p_pars->c1 * (var - p_pars->c2)), p_pars->b2)
               * std::pow( (p_pars->c3 * (p_pars->c4 - var)), -p_pars->b1);
    }

    Array parametrized_d(Array & var) {
        return std::pow(var, -p_pars->b7)
               * ( std::pow(p_pars->c1 * (var - p_pars->c2), p_pars->b3 - p_pars->w * p_pars->b2) )
               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + p_pars->w * p_pars->b1) )
               * std::pow(p_pars->c5 * (p_pars->c6 - var),  -p_pars->b5 );
    }

    Array parametrized_p(Array & var) {
        return std::pow(var, p_pars->b8)
               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + (p_pars->w - 2) * p_pars->b1))
               * ( std::pow(p_pars->c5 * (p_pars->c6 - var), 1 - p_pars->b5));
    }

    Array parametrized_v(Array & var) {
        return parametrized_eta(var) * var;
    }

    double rho_profile_int(double r, double r_shock, double rho2) {
        // density at radius r
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
        if (ia == 0 || ib == p_pars->eta.size()-1){
//            std::cerr << "Eta[0,1,2,3,4,5]="<<p_pars->eta[0]<<", "<<p_pars->eta[1]<<", "<<p_pars->eta[2]<<", "<<p_pars->eta[3]<<", "<<p_pars->eta[4]<<"\n";
//            std::cerr << "d[0,1,2,3,4,5]="<<p_pars->d[0]<<", "<<p_pars->d[1]<<", "<<p_pars->d[2]<<", "<<p_pars->d[3]<<", "<<p_pars->d[4]<<"\n";
//            std::cerr << AT << " ia=0 or ib=n-1 for eta="<<eta<<"\n";
//            exit(1);
        }
        double dint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->d);
        if (dint > 1){
            std::cerr << " Sedov profile: rho > rho2 \n";
            std::cerr << AT << "\n";
        }
        if (!std::isfinite(dint)){
            std::cerr << " NAN in ST interpolation" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        return rho2 * dint;
    }

    double Gamma_profile_int(double r, double r_shock, double vshock){
        // velocity at radius r,
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
        double vint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->vv);
        if (vshock * (2. / (p_pars->gamma + 1.)) > 1){
            std::cerr << " Sedov profile: vv > vshock \n";
            std::cerr << AT "\n";
        }
        return EQS::Gamma( vint * vshock * (2. / (p_pars->gamma + 1.)) );
    }

    double pressure_profile_int(double r, double r_shock, double p_shock){ // p' = -(gamAdi -1)Eint2/V2
        double eta = r / r_shock;
        size_t ia = findIndex(eta, p_pars->eta, p_pars->eta.size());
        size_t ib = ia + 1;
//        std::cout<<p_pars->p<<"\n";
        double pint = interpSegLin(ia, ib, eta, p_pars->eta, p_pars->p);
        return p_shock * pint;
    }
};

/// ISM density profile that a blast wave passes through
class RhoISM{
public:
    RhoISM() = default;
    void setPars(const double nn, const double A0, const double s,
                 const double R_EJ, const double R_ISM, bool set_static_ism){
        m_nn = nn;
        m_A0 = A0;
        m_s = s;
        m_R_EJ = R_EJ;
        m_R_ISM = R_ISM;

        if (set_static_ism) {
            m_GammaRho = 1;
            m_GammaRel = -1; // Gamma
            m_dGammaReldGamma = 1; // 1.
            m_dGammaRhodR = 0.; // 0.
            m_dPCBMdrho = 0.;
            m_P_cbm = 0.;
        }
        else{
            std::cerr <<  " only static ISM can be set at the present\n";
            std::cerr << AT << "\n";
            exit(1);
        }
    }

    /// only static ISM
    void evaluateRhoDrhoDrDefault(const double R, const double theta ){
        m_R = R;
        m_theta = theta;
        if (!is_overriden){  }

        /// uniform ISM density
        if (m_nn > 0.0)  { m_rho_def = m_nn * CGS::mppme; m_drhodr_def = 0.0;}
            /// step function as a ISM density
        else {
            // Inner part of the CBM :: inside the merger ejecta (high density, constant)
            if (m_R < m_R_EJ) { m_rho_def = m_A0 * std::pow(m_R_EJ, -m_s) * CGS::mppme; m_drhodr_def = 0.0; }

            // middle part :: zone of decreasing density
            if ((m_R_EJ <= m_R) && (m_R < m_R_ISM)) { m_rho_def = m_A0 * std::pow(m_R, -m_s) * CGS::mppme, m_drhodr_def = m_rho_def * (-m_s / m_R); }

            // outer part, constant low density
            if (m_R_ISM <= m_R) { m_rho_def = m_A0 * std::pow(m_R_ISM, -m_s) * CGS::mppme, m_drhodr_def = 0.0; }
        }
        if(m_rho_def < 0){
            std::cerr<<"evaluateRhoDrhoDrDefault failed m_rho_def < 0\n";
            std::cerr<<AT<<"\n";
            exit(1);
        }
    }

//    void evaluateRhoDrhoDr2(double );

public:
    double m_rho_def = -1.;
    double m_drhodr_def = -1.;
    double m_rho_ = -1.;
    double m_drhodr_ = -1.;
    double m_GammaRho = -1.;
    double m_GammaRel = -1.;
    double m_dGammaRhodR = -1.;
    double m_dGammaReldGamma = -1.;
    double m_dPCBMdrho = -1.;
    double m_P_cbm = -1;
    double m_M_cbm = -1;
    double m_CS_CBM = 0.;

    double m_rho_floor_val = 1e-5;
//    double m_rho0 = -1;
//    double m_GammaRho0 = -1;
private:
    // ------------------------------------------
    double m_nn = -1, m_A0 = -1, m_s = -1, m_R_EJ = -1, m_R_ISM = -1;
    bool is_overriden = false;
    double m_R = -1, m_theta = -1;

//    void evaluateRhoDrhoDr2(std::vector<std::unique_ptr<RadBlastWave>> &j_bws);
};

/// Description of the fluid EOS (trans-relativistic ideal gas)
class EOSadi{
public:
    enum METHODS { iPeer12, iNava13 };
    EOSadi() = default;
    void setPars( METHODS method ){ m_method = method; }
    double getGammaAdi(const double Gamma, const double beta){
        double gammaAdi;
        switch (m_method) {

            case iPeer12:
                gammaAdi= AdiabaticIndexPeer(Gamma, beta);
                break;
            case iNava13:
                gammaAdi= AdiabaticIndexNava(Gamma);
                break;
        }
        return gammaAdi;
    }
private:
    /*
    * Adiabatic index of an ideal gas (informed by NR simulations, see Peer+2012)
    */
    static double AdiabaticIndexPeer(const double Gamma, const double beta){
        double mom = Gamma * beta;
        double theta = mom / 3.0 * (mom + 1.07 * mom * mom) / (1.0 + mom + 1.07 * mom * mom); // momentrum
        double zz = theta / (0.24 + theta);
        return (5.0 - 1.21937 * zz + 0.18203 * std::pow(zz, 2) - 0.96583 * std::pow(zz, 3) + 2.32513 * std::pow(zz, 4) -
                2.39332 * std::pow(zz, 5) + 1.07136 * std::pow(zz, 6)) / 3.0;
    }
    /*
     * Simple formula for adiabatic index of an ideal fluid (see Nava+2013)
     */
    static double AdiabaticIndexNava(const double Gamma){
        return (4.0 + 1.0 / Gamma) / 3.0;
    }
private:
    METHODS m_method{};
};

/// prescriptions for a blast wave lateral spreading
class LatSpread{
public:
    LatSpread() = default;
    enum METHODS { iNULL, iAdi, iAA, iAFGPY };
    void setPars(const double aa, const double theta_max, const double thetaC, const double thetaW, METHODS method ){
        m_aa = aa; m_theta_max=theta_max; m_thetaC = thetaC; m_thetaW=thetaW; m_method=method;
    }
    double getDthetaDr(const double &Gamma, const double &R, const double &gammaAdi, const double &theta){
        double dthetadr;
        switch (m_method) {
            case iNULL:
                dthetadr = 0.;
                break;
            case iAdi:
                dthetadr = dthetadr_Adi(Gamma, R, gammaAdi, theta);
                break;
            case iAA:
                dthetadr = dthetadr_AA(Gamma, R, gammaAdi, theta);
                break;
            case iAFGPY:
                dthetadr = dthetadr_afterglopy(Gamma, R, gammaAdi, theta);
                break;
        }
        return dthetadr;
    }

private:
    /**
 * Eq. based on the sound spead in the blast (Eqs. in doi:10.1111/j.1365-2966.2004.08165.x )
 * @param Gamma
 * @param R
 * @param gammaAdi
 * @param theta
 * @return
 */
    [[nodiscard]] double dthetadr_Adi(const double Gamma, const double R, const double gammaAdi, const double theta) const {
        //    vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
        double vperp = sqrt(gammaAdi * (gammaAdi - 1.0) * (Gamma - 1.0) / (1.0 + gammaAdi * (Gamma - 1.0))) * CGS::c;
        //        return vperp / R / Gamma / beta / cgs.c if (theta < thetamax) & useSpread else 0.
        double beta = EQS::Beta(Gamma);
        return (theta < m_theta_max) ? (vperp / R / Gamma / beta / CGS::c) : 0.0;
    }

    /**
     * Method of Granot & Piran 2012 with 'a' parameter
     * @param Gamma
     * @param R
     * @param gammaAdi
     * @param theta
     * @return
     */
    [[nodiscard]] double dthetadr_AA(const double Gamma, const double R, const double gammaAdi, const double theta) const{
        //1. / (R * Gamma ** (1. + aa) * theta ** (aa)) if (theta < thetamax) & useSpread else 0.
        double Q0 = 2.;
        double g = Gamma;
        double beta = EQS::Beta(Gamma);
        double u = Gamma * beta;
        double thC = m_thetaC;
        double thW = m_thetaW;
        double dthetadr = 0;
//        if ((theta < 0.5 * M_PI) )//&& 2. * EQS::Beta(Gamma) * Gamma * thetaC < 1.)
//            dthetadr = 1 / ( R * pow( Gamma, 1. + m_aa) * pow(theta, m_aa));
        //* pow(std::tan(theta),-m_aa) * EQS::Beta(Gamma); // This did not seem to be working
        return 1 / ( R * std::pow( Gamma, 1. + m_aa) * std::pow(theta, m_aa));
    }

    /**
     * This is option 7 from 'afterglowpy' package
     * There the eq. is given as dtheta/dt
     * we thus return dtheta/dr = dtheta/dt * (dR/dt)^-1, both of
     * which are taken from the 'afterglopy' code.
     */
    [[nodiscard]] double dthetadr_afterglopy(const double Gamma, const double R, const double gammaAdi, const double theta) const{
        // TODO this does not work at all...
        if (m_thetaC < 0){
            std::cerr << " thetaC is not set!" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        double Q0 = 2.0;
        double Q = sqrt(2.0)*3.0;
        double beta = EQS::Beta(Gamma);
        double g = Gamma;
        double u = Gamma * beta;
        double thC = m_thetaC;
        double thW = m_thetaW;
//        thC = thW;
        double th = theta;
        double th0 = m_theta_b0;
        double bes = 4*u*g/(4*u*u+3); // shock velocity
        bes = sqrt(gammaAdi * (gammaAdi - 1.0) * (Gamma - 1.0) / (1.0 + gammaAdi * (Gamma - 1.0)));
        double dRshdt = bes * CGS::c;

        double dThdt = 0.;
        if(Q0 * u * thC < 1) {
            double bew = 0.5*sqrt((2*u*u+3)/(4*u*u+3))*bes/g; // is it beta_{\pe}
            double tmp = (1. - beta * dRshdt / CGS::c) * bes * bes - std::pow( dRshdt / CGS::c  - beta, 2 );
            if (tmp > 0)
                bew = Gamma * sqrt( tmp );
            else
                return 0.;
            double fac = u*thC*Q < 1.0
                         ? 1.0
                         : Q*(1-Q0*u*thC) / (Q-Q0);
            if (th0 < thC)
                fac *= tan(0.5*th0)/tan(0.5*thC); //th0/thC;
            dThdt = fac * bew * CGS::c / R;
        }
        if (dThdt < 0.) {
            std::cerr << " dtheta/dr < 0\n ";
            std::cerr << "\n ";
            exit(1);
        }
        return dThdt * (1. / beta / CGS::c );
    };

public:
    double m_aa=-1.;
    METHODS m_method{};
    double m_theta_b0 = -1.;
private:
    double m_theta_max=-1., m_thetaC=-1, m_thetaW=-1;
};

#endif //SRC_BASE_EQUATIONS_H
