//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_COMPONENTS_H
#define SRC_BLASTWAVE_COMPONENTS_H

//#include "../utilitites/pch.h"
//#include "../utilitites/utils.h"
//#include "../utilitites/interpolators.h"
#include "../utilitites/quadratures.h"

inline namespace EQS{
    inline double MomFromGam(const double gam){
        return std::sqrt(gam * gam - 1.);
    }
    inline double GamFromMom(const double mom){
        return std::sqrt(1.0+mom*mom);
    }
    inline double BetFromMom(const double mom){
        return mom / EQS::GamFromMom(mom);
    }
    inline double MomFromBeta(const double beta){
        return beta * beta / std::sqrt(1.0 - beta * beta);
    }
    inline double BetaFromGamma(const double Gamma){
        /// return sqrt(1.0 - pow(Gamma, -2)); || std::sqrt(1.-std::pow(Gamma*Gamma,-2.))
        /// return (1. / Gamma) * sqrt( (Gamma - 1.) * (Gamma + 1.) );
        return BetFromMom(MomFromGam(Gamma));
    }
    inline double GammaFromBeta(const double beta){
        /// return sqrt(1. + (betaSh * betaSh / (1. - betaSh * betaSh)));
        /// return sqrt(1.0 / (1.0 - (betaSh * betaSh)));
        return GamFromMom(MomFromBeta(beta));
    }


    double GammaRel(const double Gamma1, const double Gamma2){
//        return Gamma1 * Gamma2 * (1. - EQS::Beta(Gamma1) * EQS::Beta(Gamma2));
        return Gamma1 * Gamma2 - sqrt(Gamma1 * Gamma1 - 1.) * sqrt(Gamma2 * Gamma2 - 1.);
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
     * computeSynchrotronEmissivityAbsorptionAnalytic the swept-up mass by the expanding blast wave
     * Accounts for lateral spreading if 'aa' < 0
     * https://arxiv.org/pdf/1203.5797.pdf
     * https://academic.oup.com/mnras/article/421/1/570/990386 (for introduction of 'a' parameter)
     *
     * basic form was dmdr = 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / m_pars["m_scale"]
     *
     */
    double dmdr(const double Gamma, const double RR, const double thetaE, const double theta, const double rho, const double aa){

        // First term: change in swept-up mass due to the change in solid angle
        // double t1 = 0. if (aa < 0) else (1. / 3.) * np.sin(itheta) / (iGamma ** (1. + aa) * itheta ** (aa));
        double t1 = (aa < 0) ? 0.0 : sin(theta) / (3.0 * std::pow(Gamma, 1.0 + aa) * std::pow(theta, aa));
        // Second term: change in swept-up mass due to radial expansion
        double t2 = (cos(thetaE) - cos(theta));
        return 2.0 * CGS::pi * rho * (t1 + t2) * RR * RR;
    }
    double dmdr(const double Gamma, const double RR, const double dthetadR, const double theta, const double rho){
        return 2.*CGS::pi*rho*((1.- cos(theta)) + (1./3.)* sin(theta)*RR*dthetadR)*RR*RR;
    }

    /*
     * Evolution equation for the lorentz factor (see Peer+2012)
     */
    double dgdr(double const &M0, double const &Gamma, long double const &beta,
                double const &mm, double const &gamma_adi, const double dmdr){
        long double numerator = -(gamma_adi * (Gamma * Gamma - 1.0) - (gamma_adi - 1.0) * Gamma * beta * beta);
        long double denominator = M0 + mm * (2.0 * gamma_adi * Gamma - (gamma_adi - 1) * (1. + std::pow(Gamma, -2)));
        return (double)(numerator / denominator) * dmdr;
    }

    /*
     * Compute mass of the shall assuming relativistic, microphysics dominated shell
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
     * Compute comoving time
     */
    double initTComov(const double R0, const double beta0, const double Gamma0){
        return R0 / (beta0 * Gamma0 * CGS::c);
    }

    /*
     * Compute time in the lab frame
     */
    double init_elapsed_time(const double &r0, const double &mom0, const bool &use_spread){
//        double beta0 = Beta(Gamma0);
        double result;
        double xx = 1. + mom0*mom0*mom0/(std::sqrt(1.+mom0*mom0));//Gamma0 * Gamma0 * beta0;
        double beta0 = BetFromMom(mom0);
        if (use_spread){
            double dThetadr0 = 0.0;
            result= r0 / (CGS::c * beta0) * (sqrt(1. + r0 * r0 * dThetadr0 * dThetadr0)) - r0 / CGS::c;
        } else {
            result= r0 / (CGS::c * xx * (1.0 + beta0));
        }
        if (result!=result){
            std::cerr<<AT<<" Failed to init TT"
                     <<" i="<<0
                     <<" Rs[0]"<<r0
                     <<" mom0[0]"<<mom0
                     <<" use_spread="<<use_spread
                     <<" result[0]="<<result
                     <<" \n";
            exit(1);
        }
        return result;
    }

    double evalElapsedTime(double R, double mom, double dthetadr, bool spread){
        double xx = 1. + mom*mom*mom/(std::sqrt(1.+mom*mom));
        double beta = BetFromMom(mom);
        double dttdr = 0.;
        if (spread)
            dttdr = 1. / (CGS::c * beta) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
        else
            dttdr = 1. / (CGS::c * xx * (1. + beta));
//        double dttdr = 1. / (CGS::c * xx * (1. + betaSh));
        return dttdr;
    }

    /*
     * Integrate the elapsed time in the lab frame
     */
    double integrate_elapsed_time(unsigned int i,
                                  Vector &Rss,
                                  Vector &Gammas,
                                  Vector &thetas,
                                  const bool &use_spread){

        // integrate only to 'i'

        size_t n = i+1;
        double result = 0.0;
        Vector Rs(n);
        for (size_t j = 0; j < n; j++) Rs[j] = Rss[j];
        Vector integrand( n );
        Vector dthetadr( i+1 );
        if (use_spread){
            // considering lateral spreading
            dthetadr[0] = 0.0;
            for (unsigned int j = 1; j < n; j++) {
                double beta = BetaFromGamma(Gammas[j]);
                dthetadr[j] = (thetas[j] - thetas[j - 1]) / (Rs[j] - Rs[j - 1]);
                integrand[j] = 1.0 / (CGS::c * beta) * sqrt(1. + Rs[j] * Rs[j] * dthetadr[j] * dthetadr[j]) - 1. / (CGS::c);
            }
            auto func = [&](double x){ return ( Interp1d(Rs,integrand)).Interpolate(x, Interp1d::iLagrangeLinear); };
            result = Simpson38(Rs[0], Rs[i], 100, func);
        } else {
            // No spreading
            for (unsigned int j = 1; j < n; j++) {
                double beta = BetaFromGamma(Gammas[j]);
                integrand[j] = 1.0 / (CGS::c * Gammas[j] * Gammas[j] * beta * (1.0 + beta));
            }
            auto func = [&](double x){ return ( Interp1d(Rs,integrand)).Interpolate(x, Interp1d::iLagrangeLinear); };
            result = Simpson38(Rs[0], Rs[i], 100, func);
        }
        if (result!=result){
            std::cerr<< AT
                     <<" Failed to computeSynchrotronEmissivityAbsorptionAnalytic TT"
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
    double shock_delta_lu(const double & R, const double & Gamma) {
        return R / Gamma;
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
        return (Gamma+1.0)*(adi*(Gamma-1.0)+1.0) / ( adi * (2.0-adi)*(Gamma-1.0) + 2.0);
    }

    double dRsh_dt (const double Gamma, const long double beta, const double GammaSh){
        long double betaSh = EQS::BetaFromGamma(GammaSh);
        long double dRsh_dt_old = betaSh * CGS::c;
//        long double dRsh_dt_new = betaSh / (Gamma * GammaSh * (1. - betaSh * betaSh)) * CGS::c;// / Gamma;
        return (double)dRsh_dt_old;
//        std::cout << "dr= " << dRsh_dt_old << " !/ " << dRsh_dt_new << "\n";
    }

    // *************************** NAVA et al 2013 *******************

    inline double get_GammaEff(const double Gamma, const double gammaAdi){
        return (gammaAdi * Gamma * Gamma - gammaAdi + 1.0) / Gamma;
    }

    inline double get_dGammaEffdGamma(const double Gamma, const double gammaAdi){
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
        if (Eint2 == 0)
            return 0.;
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
    inline long double get_Gamma43(const double Gamma,
                                   const double Gamma0,
                                   const long double beta,
                                   const long double beta0){
        /*
         * Idea:
         *
         * Gamma43 = Gamma * Gamma0 * (1 - betaSh * beta0)
         * However, this is numerically... difficult, as 1-betaSh*beta0 is subjected to truncation error a lot
         * Using expansions we write
         *
         * Gamma43 = Gamma*Gamma0 - np.sqrt(Gamma0**2 - 1) * np.sqrt(Gamma**2 - 1) # -- wolfram. Creats jump
         * if betaSh < 0.999:
         * Gamma43 = Gamma * Gamma0 * (1 - betaSh * beta0)
         *
         */
        const long double beta_switch = 0.9999;
        long double Gamma43 = Gamma * Gamma0 * (1.0 - beta * beta0); // Eq.(B8) in Nava+13
        long double Gamma43_ = Gamma * Gamma0 * (1.0 / (Gamma0 * Gamma0)
                                                + 1.0 / (Gamma * Gamma)
                                                - 1.0 / (Gamma * Gamma) / (Gamma0 * Gamma0)) / (1.0 + beta * beta0);
//        if (Gamma43_ < 1. || Gamma43 < 1.)
//            int x = 1;

//        double gamma43_minus_one;
//        if (Gamma > beta_switch * Gamma0){
//            gamma43_minus_one = 1.0 - 1.0;
//        } else {
//            gamma43_minus_one = Gamma * Gamma0 * \
//                            (1.0 / (Gamma0 * Gamma0) + 1.0 / (Gamma * Gamma) - 1.0 / (Gamma * Gamma) / (Gamma0 * Gamma0)) / \
//                            (1 + beta * beta0) - 1.0;
//        }
//        return gamma43_minus_one;
        return Gamma43_ >= 1. ? Gamma43_ : 1.0;
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
                                       const long double &dgamma43dGamma,
                                       const long double &gamma43){
        return gammaAdi3 * (1.0 + 1.0 / (Gamma * Gamma))
              - dgamma43dGamma / 3.0 / gamma43*gamma43 * (Gamma - 1.0 / Gamma) - 1.0/(Gamma * Gamma);

    }

    /*
     * See Nava+2013 for the details (appendix on the reverse shock implementation)
     */
    double get_dGammadR_fs_rs(const double &Gamma, const double &Gamma0, const double &gammaAdi,
                              const double &dlnrho1dR, const double &M2, const double &dM2dR,
                              const double &dlnrho4dR, const double &M3, const long double &dM3dR,
                              const double &Eint2, const double &Eint3, const double &gammaAdi3,
                              const long double &Gamma43){

        // Using asymptotic approximation of betaSh
        // 0.5 / Gamma0 + 0.5 / Gamma ** 2 * (0.5 / Gamma0 - Gamma0 - 3. * Gamma0 / 8. / Gamma ** 2)
        long double dgamma43dGamma = get_dgamma43dGamma(Gamma0, Gamma);

        // (gammaAdi3 * Gamma ** 2 - gammaAdi3 + 1) / Gamma # (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
        long double GammaEff3 = get_GammaEff(Gamma, gammaAdi3);

        // gammaAdi3 * (1. + Gamma ** -2) - dgamma43dGamma / 3. / (Gamma43 + 1.) ** 2 * (Gamma - 1. / Gamma) - Gamma ** -2
        long double dGammaEff3dGamma = get_dGammaEff3dGamma(Gamma, gammaAdi3, dgamma43dGamma, Gamma43);

        // 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.
        double dGammaEffdGamma = get_dGammaEffdGamma(Gamma, gammaAdi);
        double GammaEff = get_GammaEff(Gamma, gammaAdi);  // (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

        double f_2 = GammaEff * (gammaAdi - 1.0) * Eint2 / Gamma;
        double h_2 = GammaEff * (gammaAdi - 1.0) * Eint2 * (dM2dR / M2 - dlnrho1dR);

        long double fh_factor3 = GammaEff3 * (gammaAdi3 - 1.0) * Eint3;
        long double f_3 = fh_factor3 / (Gamma43) * dgamma43dGamma;
        long double h_3 = 0.0;
        if (Eint3 != 0)
            h_3 = fh_factor3 * (dM3dR / M3 - dlnrho4dR);

        long double dGammadR = -((Gamma - 1.0) * (GammaEff + 1.0) * dM2dR
                                + (Gamma - Gamma0 + GammaEff3 * (Gamma43 - 1.0)) * dM3dR - h_2 - h_3)
                               / ((M2 + M3) + Eint2 * dGammaEffdGamma + Eint3 * dGammaEff3dGamma + f_2 + f_3);
//        if (dGammadR > 0 and Gamma > .95 * Gamma0){
//            std::cerr << " error dGammadR > 0 and Gamma > .59 * Gamma0 \n";
//            exit(1);
//        }
        return (double)dGammadR;

    }
    /*
     * Density in the unshocked region of the ejecta (region 4 in Nava+2013)
     */
    double rho4( const double &R, const double &deltaR4, const double &tprompt, const double &beta0,
                 const double &M0, const double &theta_a, const double &theta_b0, bool use_exp){
        double alpha_of = tprompt * beta0 * CGS::c;
        double rho4_fac_1 = M0 / (2.0 * alpha_of * CGS::pi * (cos(theta_a) - cos(theta_b0)));
        double rho4_fac = rho4_fac_1 / (R * R);
        if (use_exp)
            return rho4_fac * exp(-deltaR4 / alpha_of); // assuming exponential cutoff
        else
            return rho4_fac;
    }

    double get_U_e_rs(const double &Gamma, const double &M3, const double &Eint3, const double &rho4){
        double rho3 = 4.0 * Gamma * rho4; // TODO is it correct iwth Gamma and not Gamma43?
        double V3 = M3 / rho3;            // comoving volume
        double U_e = Eint3 / V3;          // comoving energy density
        return U_e;
    }


}

namespace MEQS{
    /// Binding energy of the NS
    //  Prescription from Lattimer and Prakash (2001)
    double e_bind( const double ns_mass, const double ns_radius){
        double num = CGS::gravconst*ns_mass;
        double den = ns_radius * CGS::c*CGS::c - 0.5*CGS::gravconst * ns_mass;
        double out = 0.6 * ns_mass * CGS::c*CGS::c * num / den;
        return out;
    }
    double e_rot( const double ns_inertia, const double omega ){
        return 0.5 * ns_inertia * omega * omega;
    }
}

/// Blandford-McKee self-simialr solution
class BlandfordMcKee2{
    struct Pars{
        double Etot = -1;
        double k = 0;
        double Rscale = -1;
        double gamAdi = 4/3.;
        double eta     = 1.e-5;         //          eta = p/(rho*c^2)
        double n_ext = 1.e0;            // external medium number density

        double n0      = 1.;            // cm-3:    CBM number density
        double rho0    = n0*CGS::mp;    // g.cm-3:  comoving CBM mass density
        double rhoNorm = rho0;          // n0*mp_;       // g.cm-3:  comoving CBM mass density
        double lNorm = CGS::c;          // distance normalised to c
        double vNorm = CGS::c;          // velocity normalised to c
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
//        _p.n0 = n_ism ;
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
        Vector f{};//, 1e5);
        Vector eta{};
        Vector d{};
        Vector p{};
        Vector vv{};
    };
    Pars * p_pars{};
public:
    SedovTaylor() {
        p_pars = new Pars();
    }
    ~SedovTaylor() { delete p_pars; }
    void setPars(double gamma, size_t nu, double w, int size){
        std::cout << " allocating n="<<size<<" for Sedov-Taylor profile!\n";
        p_pars->f.resize(size);
        p_pars->eta.resize( p_pars->f.size() );

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

        p_pars->f = TOOLS::MakeLogspaceVec(log10(f_min), 0., (int)p_pars->f.size()); // log10(f_min)
        print_x_as_numpy(p_pars->f, 100, "f");


//                p_pars->c2 ? p_pars->w1 > w : p_pars->c6
        bool use_uniform_log_grid = false;
        if (use_uniform_log_grid){
            p_pars->f = TOOLS::MakeLogspaceVec(log10(f_min), 0., (int)p_pars->f.size());
        }
        else {
            double tmp_fmin = log10(f_min);
            double tmp_fmax = tmp_fmin;
            std::vector<size_t> segment_lengths = {10000, 1000, 100, 100};
            VecVector segments{};
            size_t tot_size = 0;
            for (size_t i_segment = 0; i_segment < segment_lengths.size() - 1; i_segment++) {
                double _p = (double) i_segment - (double) segment_lengths.size() - 1;
                double step = std::pow(10, (double) _p);
                if (tmp_fmax == 1) tmp_fmax = 0;
                //            std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<tmp_fmin + step<<"\n";
                segments.emplace_back(TOOLS::MakeLogspaceVec(tmp_fmin, tmp_fmin + step, (int) segment_lengths[i_segment]));
                //            std::cout << segments[i_segment] << "\n";
                tmp_fmin = tmp_fmin + step;
                tot_size = tot_size + segment_lengths[i_segment];
            }
            segments.emplace_back(TOOLS::MakeLogspaceVec(tmp_fmin, 0, (int) segment_lengths[segment_lengths.size() - 1]));
//            std::cout << segments[segments.size() - 1] << "\n";
            tot_size = tot_size + segment_lengths[segment_lengths.size() - 1];
            //        std::cout << " tmp_min="<<tmp_fmin<<" tmp_max="<<0<<"\n";
            Vector tmp(tot_size);
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
        std::vector<size_t> idx = sort_indexes(p_pars->eta);
        p_pars->f = sort_by_indexes(p_pars->f, idx);
        p_pars->eta = sort_by_indexes(p_pars->eta, idx);

        p_pars->d = parametrized_d(p_pars->f);
        p_pars->p = parametrized_p(p_pars->f);
        p_pars->vv = parametrized_v(p_pars->f);

        print_x_as_numpy(p_pars->eta,100,"eta");
        print_x_as_numpy(p_pars->d,100,"d");

        if (p_pars->eta[0] > 0.){
            std::cerr << " in Sedov Taylow, the extend to eta is too short. Addd more 0.00000..." << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        // Finally Calculate the normalization of R_s:
//        auto tmp = std::pow(p_pars->vv, 2);
        Vector integral_(p_pars->eta.size());// = std::pow(p_pars->eta, nu - 1) * (p_pars->d * std::pow(p_pars->vv, 2.) + p_pars->p);
        for (size_t i = 0; i < p_pars->eta.size(); i++)
            integral_[i] = std::pow(p_pars->eta[i], nu - 1) * (p_pars->d[i] * std::pow(p_pars->vv[i], 2.) + p_pars->p[i]);

        Vector integral( integral_.size() - 1 );
        Vector deta (integral_.size() - 1);
        double integ = 0;
        for (size_t i = 0; i < integral.size(); i++){
            integral[i] = 0.5 * (integral_[i+1] + integral_[i]);
            deta[i] = (p_pars->eta[i+1] - p_pars->eta[i]);
            integ += integral[i] * deta[i];
        }
        double alpha = integ * (8 * p_pars->c0) / ((std::pow(gamma, 2) - 1.) * std::pow(nu + 2. - w, 2));
//        p_pars->_c = pow(1. / alpha, 1. / (nu + 2 - w))


    }

    Vector parametrized_eta(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], -p_pars->b6)
                    * std::pow((p_pars->c1 * (var[i] - p_pars->c2)), p_pars->b2)
                     * std::pow( (p_pars->c3 * (p_pars->c4 - var[i])), -p_pars->b1);
        return res;
//        return std::pow(var, -p_pars->b6) * std::pow((p_pars->c1 * (var - p_pars->c2)), p_pars->b2)
//               * std::pow( (p_pars->c3 * (p_pars->c4 - var)), -p_pars->b1);
    }

    Vector parametrized_d(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], -p_pars->b7)
                     * ( std::pow(p_pars->c1 * (var[i] - p_pars->c2), p_pars->b3 - p_pars->w * p_pars->b2) )
                     * ( std::pow(p_pars->c3 * (p_pars->c4 - var[i]), p_pars->b4 + p_pars->w * p_pars->b1) )
                     * std::pow(p_pars->c5 * (p_pars->c6 - var[i]),  -p_pars->b5 );
        return res;
//        return std::pow(var, -p_pars->b7)
//               * ( std::pow(p_pars->c1 * (var - p_pars->c2), p_pars->b3 - p_pars->w * p_pars->b2) )
//               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + p_pars->w * p_pars->b1) )
//               * std::pow(p_pars->c5 * (p_pars->c6 - var),  -p_pars->b5 );
    }

    Vector parametrized_p(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < var.size(); i++)
            res[i] = std::pow(var[i], p_pars->b8)
                     * ( std::pow(p_pars->c3 * (p_pars->c4 - var[i]), p_pars->b4 + (p_pars->w - 2) * p_pars->b1))
                     * ( std::pow(p_pars->c5 * (p_pars->c6 - var[i]), 1 - p_pars->b5));
        return res;
//        return std::pow(var, p_pars->b8)
//               * ( std::pow(p_pars->c3 * (p_pars->c4 - var), p_pars->b4 + (p_pars->w - 2) * p_pars->b1))
//               * ( std::pow(p_pars->c5 * (p_pars->c6 - var), 1 - p_pars->b5));
    }

    Vector parametrized_v(Vector & var) {
        Vector res(var.size());
        for (size_t i = 0; i < p_pars->eta.size(); i++)
            res[i] = parametrized_eta(var)[i] * var[i];
        return res;
//        return parametrized_eta(var) * var;
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
        return EQS::GammaFromBeta( vint * vshock * (2. / (p_pars->gamma + 1.)) );
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
    std::unique_ptr<logger> p_log;
public:
    RhoISM(int loglevel){
        p_log=std::make_unique<logger>(std::cout, std::cerr, loglevel, "RhoISM");
    }
    void setPars(const double nn, const double A0, const double s,
                 const double R_EJ, const double R_ISM, bool set_static_ism){
        m_nn = nn;
        m_A0 = A0;// * CGS::mppme;
        m_s = (double)s;
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
        if (!is_overriden){

        }

        /// uniform ISM density
        if (m_nn > 0.0)  {
            m_rho_def = m_nn * CGS::mppme;
            m_drhodr_def = 0.0;
        }
            /// step function as a ISM density
        else {
            if ((m_s < 0)or(m_A0<0)){
                (*p_log)(LOG_ERR, AT) << "if n_ism < 0; positive A0 and s are expected (wind). Given A0="<<m_A0<<" s="<<m_s<<"\n";
                exit(1);
            }

            // Inner part of the CBM :: inside the merger ejecta (high density, constant)
            if (m_R < m_R_EJ) {
                m_rho_def = m_A0 * std::pow(m_R_EJ, -m_s) * CGS::mppme;
                m_drhodr_def = 0.0;
            }

            // middle part :: zone of decreasing density
            if ((m_R_EJ <= m_R) && (m_R < m_R_ISM)) {
                m_rho_def = m_A0 * std::pow(m_R, -m_s) * CGS::mppme,
                m_drhodr_def = m_rho_def * (-m_s / m_R);
            }

            // outer part, constant low density
            if (m_R_ISM <= m_R) {
                m_rho_def = m_A0 * std::pow(m_R_ISM, -m_s) * CGS::mppme,
                m_drhodr_def = 0.0;
            }
        }
        if(m_rho_def < 0){
            (*p_log)(LOG_ERR, AT)<<" evaluateRhoDrhoDrDefault failed m_rho_def < 0\n";
            exit(1);
        }
        if(m_rho_def > 1.){
            (*p_log)(LOG_ERR, AT) << "wrong walue of rho_ism="<<m_rho_def<<" for n_ism="<<m_nn<<" A0="<<m_A0<<" s="<<m_s<<"\n";
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
    double getGammaAdi(const double Gamma, const long double beta){
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
    static double AdiabaticIndexPeer(const double Gamma, const long double beta){
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
    enum METHODS { iNULL, iAdi, iAA, iAFGPY, iOUR };
    void setPars(const double aa, const double theta_max, const double thetaC, const double thetaW, METHODS method ){
        m_aa = aa; m_theta_max=theta_max; m_thetaC = thetaC; m_thetaW=thetaW; m_method=method;
        if ((thetaC < 0) || (thetaW < 0)){
            std::cerr << AT << " thetaC="<<thetaC<<" thetaW="<<thetaW<<"\n";
            exit(1);
        }
    }
    double getDthetaDr(double Gamma, double GammaSh, double R, double gammaAdi, double theta){
        double dthetadr;
        switch (m_method) {
            case iNULL:
                dthetadr = 0.;
                break;
            case iAdi:
                dthetadr = dthetadr_Adi(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iAA:
                dthetadr = dthetadr_AA(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iAFGPY:
                dthetadr = dthetadr_afterglopy(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iOUR:
                dthetadr = dthetadr_our(Gamma, GammaSh, R, gammaAdi, theta);
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
    [[nodiscard]] double dthetadr_Adi(double Gamma, double GammaSh, double R, double gammaAdi, double theta) const {
        //    vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
        double vperp = sqrt(gammaAdi * (gammaAdi - 1.0) * (Gamma - 1.0) / (1.0 + gammaAdi * (Gamma - 1.0))) * CGS::c;
        //        return vperp / R / Gamma / betaSh / cgs.c if (theta < thetamax) & useSpread else 0.
        double beta = EQS::BetaFromGamma(Gamma);
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
    [[nodiscard]] double dthetadr_AA(double Gamma, double GammaSh, double R, double gammaAdi, double theta) const{
        //1. / (R * Gamma ** (1. + aa) * theta ** (aa)) if (theta < thetamax) & useSpread else 0.
        double Q0 = 2.;
        double g = Gamma;
        double beta = EQS::BetaFromGamma(Gamma);
        double u = Gamma * beta;
        double thC = m_thetaC;
        double thW = m_thetaW;
        double dthetadr = 0;
//        if ((theta < 0.5 * M_PI) )//&& 2. * EQS::Beta(Gamma) * Gamma * thetaC < 1.)
//            dthetadr = 1 / ( R * pow( Gamma, 1. + m_aa) * pow(theta, m_aa));
        //* pow(std::tan(theta),-m_aa) * EQS::Beta(Gamma); // This did not seem to be working
        return 1. / ( R * std::pow( Gamma, 1. + m_aa) * std::pow(theta, m_aa));
    }

    /**
     * This is option 7 from 'afterglowpy' package
     * There the eq. is given as dtheta/dt
     * we thus return dtheta/dr = dtheta/dt * (dR/dt)^-1, both of
     * which are taken from the 'afterglopy' code.
     * Works only for 'adaptive" eats
     */
    [[nodiscard]] double dthetadr_afterglopy(double Gamma, double GammaSh, double R, double gammaAdi, double theta) const{
        if (m_thetaC < 0){
            std::cerr << " thetaC is not set!" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        double dThdt = 0.;
        double Q0 = 2.0;
        double Q = sqrt(2.0)*3.0;
        double th = theta;
        double u = EQS::MomFromGam(Gamma);
        double thC = m_thetaC;
        double th0 = m_theta_b0;
        double bes = 4*u*Gamma/(4*u*u+3);
        double val = sqrt((2*u*u+3)/(4*u*u+3));
        double bew = 0.5*val*bes/Gamma;

        if((th < 0.5 * CGS::pi) && (Q0 * u * thC < 1)){
            double fac = u*thC*Q < 1.0 ? 1.0 : Q*(1-Q0*u*thC) / (Q-Q0);
            if (th0 < thC)
                fac *= tan(0.5*th0)/tan(0.5*thC); //th0/thC;
            dThdt = fac * bew * CGS::c / R;
        }
        if (dThdt < 0.) {
            std::cerr << " dtheta/dr < 0\n ";
            std::cerr << "\n ";
            exit(1);
        }
        return dThdt * (1. / bes / CGS::c );
    }

    [[nodiscard]] double dthetadr_our(double Gamma, double GammaSh, double R, double gammaAdi, double theta) const{
        double Q0 = 2.0;
        double Q = sqrt(2.0)*3.0;
        double mom = EQS::MomFromGam(GammaSh);
        double vperp = sqrt(gammaAdi * (gammaAdi - 1.0) * (Gamma - 1.0)
                     / (1.0 + gammaAdi * (Gamma - 1.0))) * CGS::c; // Huang+2000
        if((theta < 0.5 * CGS::pi) && (Q0 * mom * m_thetaC < 1)){
            double fac = 1.;
            if (mom * m_thetaC * Q >= 1.0)
                fac = Q * (1. - Q0 * mom * m_thetaC) / (Q - Q0);  // Ryan+2019
            if (m_theta_b0 < m_thetaC)
                fac *= tan(0.5*m_theta_b0)/tan(0.5*m_thetaC); //th0/thC; // Ryan+2019
            double beta = EQS::BetaFromGamma(GammaSh);
            return vperp / R / GammaSh / beta / CGS::c * fac;
        }
        return 0.;
    }

    [[nodiscard]] double dthetadr_afterglopy_old(double Gamma, double GammaSh, double R, double gammaAdi, double theta) const{

        if (m_thetaC < 0){
            std::cerr << " thetaC is not set!" << "\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        double Q0 = 2.0;
        double Q = sqrt(2.0)*3.0;
        double beta = EQS::BetaFromGamma(Gamma);
        double g = Gamma;
        double u = Gamma * beta;
        double thC = m_thetaC;
        double thW = m_thetaW;
//        thC = thW;
        double th = theta;
        double th0 = m_theta_b0;
        double bes2 = 4*u*g/(4*u*u+3); // shock velocity
        double bes = std::sqrt(gammaAdi * (gammaAdi - 1.0) * (Gamma - 1.0) / (1.0 + gammaAdi * (Gamma - 1.0)));
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

#endif //SRC_BLASTWAVE_COMPONENTS_H
