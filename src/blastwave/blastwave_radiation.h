//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLASTWAVE_RADIATION_H
#define SRC_BLASTWAVE_RADIATION_H

//#include "blastwave_components.h"
//#include "blastwave_base.h"
#include "blastwave_base.h"


/// use precomputed emissivity and absorption for Piece Wise EATS and blastwave structure
static void fluxDensPieceWiseWithComov(
        double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
        double r, double & ctheta, double theta, double phi,
        size_t ia, size_t ib, double ta, double tb, double mu,
        double t_e, double nu_obs, void * params){

    auto * p_pars = (struct Pars *) params;
    auto & m_data = p_pars->m_data;

    /// if blastwave is not evolved
    if (p_pars->i_end_r==0)
        return;

    double Gamma = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iGamma]);
    double beta = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::ibeta]);
    double GammaShock = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iGammaFsh]);
    double dr = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::ithickness]);
    double nprime = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::irho2]) / CGS::mp;
    double nprotons = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iM2]) / CGS::mp;
    double rsh = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iRsh]);
    double acc_fac = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iacc_frac]);

    /// compute Doppler factor
    double a = 1.0 - beta * mu; // beaming factor
    double delta_D = Gamma * a; // doppler factor

    /// compute the comoving obs. frequency from given one in obs. frame
    double nuprime = (1.0 + p_pars->z) * nu_obs * delta_D;

    flux_dens = p_pars->p_mphys->fluxDens(
            EjectaID2::STUCT_TYPE::ipiecewise, theta, p_pars->ncells,
            ia, ib, nuprime, Gamma, GammaShock, acc_fac,
            mu, r, dr, rsh, nprime, nprotons, m_data[BW::Q::iR]
            );

    /// save the result in image
    ctheta = interpSegLin(ia, ib, ta, tb, t_e, m_data[BW::Q::ictheta]);

    if (p_pars->m_type == BW_TYPES::iFSRS && p_pars->do_rs_radiation) {
        if (m_data[BW::Q::ithickness_rs][ia] > 0 and m_data[BW::Q::ithickness_rs][ib] > 0) {

            double GammaShock_rs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iGammaRsh]);
            double dr_rs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::ithickness_rs]);
            double nprime_rs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::irho3]) / CGS::mp;
            double ne_rs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iM3]) / CGS::mp;
            double r_rsh = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iRrsh]);
            double acc_fac_rs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iacc_frac_rs]);

            double flux_dens_rs = p_pars->p_mphys->fluxDens(
                    EjectaID2::STUCT_TYPE::ipiecewise, theta, p_pars->ncells,
                    ia, ib, nuprime, Gamma, GammaShock_rs, acc_fac_rs,
                    mu, r, r_rsh, dr_rs, nprime_rs, ne_rs, m_data[BW::Q::iR]
                    );

            flux_dens += flux_dens_rs;
        }
    }
    flux_dens *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
}

/// use precomputed emissivity and absorption for adaptive EATS and blastwave structure
static void fluxDensAdaptiveWithComov(
        double & flux_dens, double & r, double & ctheta, double theta, double phi,
        size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs, void * params){

    auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
    auto & p_syna = p_pars->p_mphys;//->getAnSynch();
    auto & Dt = p_pars->m_data;
    auto & tburst = Dt[BW::Q::itburst];
    auto & r_arr = Dt[BW::Q::iR];

    /// interpolate the exact radial position of the blast that corresponds to the req. obs time
    r = interpSegLog(ia, ib, t_e, tburst, r_arr);

    if ((r <= 0.0) || !std::isfinite(r)) {
        (*p_pars->p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                     << " Current R grid us ["
                                     << r_arr[0] << ", "
                                     << r_arr[tburst.size() - 1] << "] "
                                     << "and tburst arr ["
                                     << tburst[0] << ", " << tburst[p_pars->nr - 1]
                                     << "] while the requried obs_time=" << t_obs
                                     << "\n";
        return;
    }
    double Gamma = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iGamma]);
    double beta = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::ibeta]);
    double GammaShock = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iGammaFsh]);
    double dr = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::ithickness]);
    double ne = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iM2]) / CGS::mp;
    double n_prime = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::irho2]) / CGS::mp;
    double rsh = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iRsh]);
    double acc_fac = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iacc_frac]);

    /// compute Doppler factor
    double a = 1.0 - beta * mu; // beaming factor
    double delta_D = Gamma * a; // doppler factor

    /// computeSynchrotronEmissivityAbsorptionAnalytic the comoving obs. frequency from given one in obs. frame
    double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;

    /// save the result in image
    ctheta = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::ictheta]);
    theta = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::itheta]);

    flux_dens = p_pars->p_mphys->fluxDens(
            EjectaID2::STUCT_TYPE::iadaptive, theta, p_pars->ncells,
            ia, ib, nuprime, Gamma, GammaShock, acc_fac,
            mu, r, rsh, dr, n_prime, ne, r_arr);

    if (p_pars->m_type == BW_TYPES::iFSRS && p_pars->do_rs_radiation) {
        if (Dt[BW::Q::ithickness_rs][ia] > 0 and Dt[BW::Q::ithickness_rs][ib] > 0) {
            double GammaShock_rs = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iGammaRsh]);
            double dr_rs = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::ithickness_rs]);
            double ne_rs = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iM3]) / CGS::mp;
            double n_prime_rs = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::irho3]) / CGS::mp;
            double r_rsh = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iRrsh]);
            double acc_fac_rs = interpSegLog(ia, ib, r, r_arr, Dt[BW::Q::iacc_frac_rs]);

            double flux_dens_rs = p_pars->p_mphys_rs->fluxDens(
                    EjectaID2::STUCT_TYPE::iadaptive, theta, p_pars->ncells,
                    ia, ib, nuprime, Gamma, GammaShock_rs, acc_fac_rs,
                    mu, r, r_rsh, dr_rs, n_prime_rs, ne_rs, r_arr);

            flux_dens += flux_dens_rs;
        }
    }
}


/// compute local emissivity and absoption from interpolated shock conditions
static double shock_synchrotron_flux_density(
        EjectaID2::STUCT_TYPE method_eats,
        double theta, double Gamma, double GammaShock, double acc_fac,
        double nprotons, double nprime, double acc_frac, double B,
        double gm, double gM, double gc, double Theta, double z_cool,
        double t_e, double mu, double R, double Rsh, double dr, double dr_tau,
        double nu_obs, void * params){

    auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
    auto & p_syna = p_pars->p_mphys;//->getAnSynch();

    /// relativistic effects on the emitting region
    double beta = EQS::BetaFromGamma(Gamma);
    double a = 1.0 - beta * mu; // beaming factor
    double delta_D = Gamma * a; // doppler factor

    /// properties of the emitting electrons
    double nuprime = (1.0 + p_pars->z) * nu_obs * delta_D;

    double em_prime,abs_prime;
    p_syna->setShockElectronParameters(nprime, nprotons, acc_frac,
                                       B, gm, gM, gc, Theta, z_cool);
    p_syna->computeSynchrotronEmissivityAbsorptionAnalytic(nuprime, em_prime, abs_prime);
    return p_syna->computeFluxDensity(
            method_eats,em_prime, abs_prime, Gamma, GammaShock, acc_fac, mu, R, Rsh, dr,
            nprime, nprotons, theta, p_pars->ncells);
}
#if 0
//    double em_prime, abs_prime;
    switch (p_syna->m_method_ne) {
        // default (presumably more correct version)
        case iusenprime:
            /// this is informed by van Earten et al. and afterglowpy
            p_syna->setShockElectronParameters(R,dr,nprime, Ne, acc_frac,
                                               B, gm, gM, gc, Theta, z_cool);
            p_syna->computeSynchrotronEmissivityAbsorptionAnalytic(nuprime, em_prime, abs_prime);
//            em_prime = p_syna->em;
            em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
//            abs_prime = p_syna->abs;
            abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
            dtau = optical_depth(abs_lab,dr_tau, mu, beta_shock);
            intensity = computeIntensity(em_lab, dtau, p_syna->method_tau);
//                flux_dens = (intensity * R * R * dr) * (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);

            flux_dens = (intensity * R * R * dr) ;
            break;
        case iuseNe: // TODO THIS SHOULD NOT BE USED (tau is incorrectly estimated)
            /// This is informed by the G. Lamb and Fernandez et al.
            p_syna->setShockElectronParameters(R,dr,nprime, Ne, acc_frac,
                                               B, gm, gM, gc, Theta, z_cool);
            p_syna->computeSynchrotronEmissivityAbsorptionAnalytic(nuprime, em_prime, abs_prime);
//            em_prime = p_syna->em;
            em_lab = em_prime / (delta_D * delta_D);
//            em_lab /= (delta_D); // TODO this should be from 'dr'...
//            em_lab /= (delta_D); // TODO this should be from 'dr'...
            abs_lab = abs_prime * delta_D; // TODO with Ne this might not work as we do not use 'dr' of the shock...
            dtau = optical_depth(abs_lab,dr_tau,mu,beta_shock);
            intensity = computeIntensity(em_lab, dtau, p_syna->method_tau);
//                flux_dens = intensity * (1.0 + p_eats->z) / (p_eats->d_l * p_eats->d_l) / 10;
            flux_dens = intensity * R * R * dr / Ne * nprime; // TODO why no '2'?? // why this need /10 to fit the upper result?
            break;
    }
    if ((!std::isfinite(flux_dens)) || (flux_dens < 0.)){
        (*p_pars->p_log)(LOG_ERR,AT) << " flux_dens="<<flux_dens
                        <<" em_lab="<<em_lab<<" abs_lab="<<abs_lab <<'\n';
        (*p_pars->p_log)(LOG_ERR,AT) << " wrong value in interpolation to EATS surface  \n"
                  << " nu_obs = " << nu_obs
                  << " nuprime = " << nuprime
                  << " R = " << R
                  << " mu = " << mu
                  << " dr = " << dr
                  << " m2 = " << m2
                  << " rno2 = " << rho2
                  << " B = " << B
                  << " Ne = " << Ne
                  << " nprime = " << nprime
                  << " gm = " << gm
                  << " gM = " << gM
                  << " gc = " << gc
                  << " Gamma = " << Gamma
                  << " theta = " << Theta
                  << " z_cool = " << z_cool
                  << " t_e = " << t_e
                  << " Exiting...\n";
        exit(1);
    }
#endif
//    return flux_dens;


/// use evaluated local emissivity and absoprtion to get observed flux density for Piece Wise EATS and blastwave structure
static void fluxDensPieceWiseWithObs(
        double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
        double r, double & ctheta, double theta, double phi,
        size_t ia, size_t ib, double ta, double tb, double mu,
        double t_e, double nu_obs, void * params){
    auto * p_pars = (struct Pars *) params;
    auto & m_data = p_pars->m_data;
    if (p_pars->i_end_r==0)
        return;

    double Gamma = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iGamma]);
    double GammaSh = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iGammaFsh]);
    double nprime = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::irho2]) / CGS::mp;
    double ne = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iM2]) / CGS::mp;
    double frac = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iacc_frac]);
    double thick = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::ithickness]);
    theta = interpSegLin(ia, ib, ta, tb, t_e, m_data[BW::Q::itheta]);
    ctheta = interpSegLin(ia, ib, ta, tb, t_e, m_data[BW::Q::ictheta]);
    double B = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iB]);
    double gm = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::igm]);
    double gM = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::igM]);
    double gc = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::igc]);
    double Theta = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iTheta]);
    double z_cool = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iz_cool]);
    double tburst = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::itburst]);
    double tt = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::itt]);
    double dr = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::ithickness]);
    double rsh = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iRsh]);
    double acc_fac = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iacc_frac]);
//            double cs = interpSegLog(ia, ib, ta, tb, t_e, m_data[BW::Q::iCSCBM]);

    if ((!std::isfinite(gm))||(!std::isfinite(B))||(!std::isfinite(ne))){
        (*p_pars->p_log)(LOG_ERR,AT)
                << "[ish=" << p_pars->ishell << ", " << "il=" << p_pars->ilayer << "] "
                << " nanss {"
                << " ia=" << ia << " ib=" << ib << " ta=" << ta << " tb=" << tb << " r=" << r << " mu=" << mu
                << " nu_obs=" << nu_obs << " t_e=" << t_e
                <<" phi="<<phi<<" theta="<<theta<<" ctheta="<<ctheta<<" flux_dens="<<flux_dens
                << " | " << " Gamma="<<Gamma<<" nprime="<<nprime<<" ne="<<ne<<" B="<<B
                <<"\n";
        ;
        exit(1);
    }

    double dFnu = 0.;
    if ((m_data[BW::Q::iB][ia] == 0.) || (m_data[BW::Q::iB][ib] == 0.))
        return;

#if 0
    if ((Gamma < 1. || !std::isfinite(Gamma))
        || (gm < 0.) || (!std::isfinite(gm))
        || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
        || (theta <= 0.)
        || (rho2 < 0.) || (!std::isfinite(rho2))
        || (thick <= 0.) || (GammaSh <= 1.)) {
        (*p_log)(LOG_ERR,AT)<< " Error in interpolation to EATS surface: \n"
                            << " R = " << r << "\n"
                            << " Gamma = " << Gamma << "\n"
                            << " GammaSh = " << GammaSh << "\n"
                            << " gm = " << gm << "\n"
                            << " gM = " << gM << "\n"
                            << " gc = " << gc << "\n"
                            << " B = " << B << "\n"
                            << " Theta = " << Theta << "\n"
                            << " z_cool = " << z_cool << "\n"
                            << " theta = " << theta << "\n"
                            << " rho2 = " << rho2 << "\n"
                            << " thick = " << thick << "\n"
                            << " t_e = " << t_e << "\n";
//                    exit(1);
    }
    if ((B != 0.) && (!std::isfinite(rho2))) {
        (*p_log)(LOG_ERR,AT)<< " B!=0 and rho2 is NAN \n"
                            << " Error in data \n"
                            << " R = " << r << "\n"
                            << " Gamma = " << Gamma << "\n"
                            << " gm = " << gm << "\n"
                            << " gM = " << gM << "\n"
                            << " gc = " << gc << "\n"
                            << " B = " << B << "\n"
                            << " Theta = " << Theta << "\n"
                            << " z_cool = " << z_cool << "\n"
                            << " theta = " << theta << "\n"
                            << " rho2 = " << rho2 << "\n"
                            << " thick = " << thick << "\n"
                            << " t_e = " << t_e << "\n"
                            << " Exiting...\n";
        exit(1);
    }
#endif

    double thick_tau = dr;
    flux_dens = shock_synchrotron_flux_density(
            EjectaID2::STUCT_TYPE::ipiecewise, theta,
            Gamma, GammaSh, acc_fac, ne, nprime, frac, B, gm, gM, gc,
            Theta, z_cool, t_e, mu, r, rsh, thick, thick_tau, nu_obs, p_pars);

    flux_dens *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);

    if (p_pars->m_type == BW_TYPES::iFSRS) {
        (*p_pars->p_log)(LOG_ERR,AT) << " not supported her...";
        exit(1);
    }
//    Vector & freq_arr = p_pars->p_mphys->m_freq_arr;
//    Vector & synch_em = p_pars->p_mphys->out_spectrum;
//    Vector & synch_abs = p_pars->p_mphys->out_spectrum;
//    size_t nfreqs = freq_arr.size();
}

#if 0
static void fluxDensPieceWiseWithObs(
        double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
        double r, double & ctheta, double theta, double phi,
        size_t ia, size_t ib, double ta, double tb, double mu,
        double t_obs, double nu_obs, void * params){

    auto * p_pars = (struct Pars *) params;
    auto & m_data = p_pars->m_data;
    if (p_pars->i_end_r==0)
        return;

    double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
    double GammaSh = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGammaFsh]);
    double rho2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::irho2]);
    double m2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iM2]);
    double frac = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iacc_frac]);
    double thick = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness]);
    theta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::itheta]);
    ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
    double B = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iB]);
    double gm = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igm]);
    double gM = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igM]);
    double gc = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igc]);
    double Theta = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iTheta]);
    double z_cool = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iz_cool]);
    double tburst = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itburst]);
    double tt = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itt]);
//            double cs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iCSCBM]);

    if ((!std::isfinite(gm))||(!std::isfinite(B))||(!std::isfinite(m2))){
        (*p_pars->p_log)(LOG_ERR,AT)
                <<"[ish="<<p_pars->ishell<<", "<<"il="<<p_pars->ilayer<<"] "
                <<" nanss {"
                <<" ia="<<ia<<" ib="<<ib<<" ta="<<ta<<" tb="<<tb<<" r="<<r<<" mu="<<mu
                <<" nu_obs="<<nu_obs<<" t_obs="<<t_obs
                <<" phi="<<phi<<" theta="<<theta<<" ctheta="<<ctheta<<" flux_dens="<<flux_dens
                << " | " << " Gamma="<<Gamma<<" rho2="<<rho2<<" m2="<<m2<<" B="<<B
                <<"\n";
        ;
        exit(1);
    }

    double dFnu = 0.;
    if ((m_data[BW::Q::iB][ia] == 0.) || (m_data[BW::Q::iB][ib] == 0.)){
        dFnu = 0.;
    }
    else{
#if 0
        if ((Gamma < 1. || !std::isfinite(Gamma))
                        || (gm < 0.) || (!std::isfinite(gm))
                        || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
                        || (theta <= 0.)
                        || (rho2 < 0.) || (!std::isfinite(rho2))
                        || (thick <= 0.) || (GammaSh <= 1.)) {
                        (*p_log)(LOG_ERR,AT)<< " Error in interpolation to EATS surface: \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " GammaSh = " << GammaSh << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n";
//                    exit(1);
                    }
                    if ((B != 0.) && (!std::isfinite(rho2))) {
                        (*p_log)(LOG_ERR,AT)<< " B!=0 and rho2 is NAN \n"
                                            << " Error in data \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n"
                                            << " Exiting...\n";
                        exit(1);
                    }
#endif
        double thick_tau = EQS::shock_delta(r,GammaSh); // TODO this is added becasue in Johanneson Eq. I use ncells
        flux_dens = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                   Theta, z_cool, t_obs, mu,
                                                   r, thick,  thick_tau, nu_obs, p_pars);
    }

    if (p_pars->m_type == BW_TYPES::iFSRS) {
        (*p_pars->p_log)(LOG_ERR,AT) << " not supported her...";
        exit(1);
    }

    flux_dens *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);

    Vector & freq_arr = p_pars->p_mphys->m_freq_arr;
    Vector & synch_em = p_pars->p_mphys->out_spectrum;
    Vector & synch_abs = p_pars->p_mphys->out_spectrum;
    size_t nfreqs = freq_arr.size();

    if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
        Interp2d int_em(freq_arr, m_data[BW::Q::iR], synch_em);
        Interp2d int_abs(freq_arr, m_data[BW::Q::iR], synch_abs);
        Interp1d::METHODS mth = Interp1d::iLagrangeLinear;

        double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
        double betaSh = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ibeta]);
        // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
        /// computeSynchrotronEmissivityAbsorptionAnalytic Doppler factor
        double a = 1.0 - betaSh * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        /// computeSynchrotronEmissivityAbsorptionAnalytic the comoving obs. frequency from given one in obs. frame
        double nuprime = (1.0 + p_pars->z) * nu_obs * delta_D;
        if (nuprime >= freq_arr[nfreqs-1]){
            (*p_pars->p_log)(LOG_ERR,AT)<<" freq_prime="<<nuprime
                                        <<" > grid freq[-1]="<<freq_arr[nfreqs-1]
                                        <<" increase 'freq2' parameter\n";
            exit(1);
        }
        if (nuprime <= freq_arr[0]){
            (*p_pars->p_log)(LOG_ERR,AT)<<" freq_prime="<<nuprime
                                        <<" < grid freq[0]="<<freq_arr[0]
                                        <<" decrease 'freq1' parameter\n";
            exit(1);
        }
        size_t ia_nu = findIndex(nuprime, freq_arr, freq_arr.size());
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
        double GammaShock = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGammaFsh]);
        double dr = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness]);
        double dr_tau = EQS::shock_delta(r,
                                         GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

        double beta_shock;
        switch (p_pars->method_shock_vel) {

            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> computeSynchrotronEmissivityAbsorptionAnalytic shock velocity
                beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                break;
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
        dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
        dr_tau /= ashock;
        double dtau = optical_depth(abs_lab, dr_tau, mu, beta_shock);
        double intensity = computeIntensity(em_lab, dtau, p_pars->p_mphys->method_tau);
        flux_dens = (intensity * r * r * dr) * (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
//        flux += flux_dens;
        /// save the result in image
        ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
        //  double ctheta = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::ictheta] ) ).Interpolate(r, mth );
//        image(Image::iintens, i) =
//                flux_dens / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
//        image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::imu, i) = mu_arr[i];

        if (p_pars->m_type == BW_TYPES::iFSRS) {
            if (p_pars->do_rs &&
                !(m_data[BW::Q::ithickness_rs][ia] == 0 || m_data[BW::Q::ithickness_rs][ib] == 0)) {

                Vector & freq_arr_rs = p_pars->p_mphys_rs->m_freq_arr;
                Vector & synch_em_rs = p_pars->p_mphys_rs->out_spectrum;
                Vector & synch_abs_rs = p_pars->p_mphys_rs->out_spectrum;

                Interp2d int_em_rs(freq_arr_rs, m_data[BW::Q::iR], synch_em_rs);
                Interp2d int_abs_rs(freq_arr_rs, m_data[BW::Q::iR], synch_abs_rs);
                double em_prime_rs = int_em_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                double abs_prime_rs = int_abs_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                double em_lab_rs = em_prime_rs / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                double abs_lab_rs = abs_prime_rs * delta_D; // conversion of absorption (see vanEerten+2010)

                double GammaShock_rs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma43]);
                double dr_rs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness_rs]);

                double dr_tau_rs = EQS::shock_delta(r, GammaShock_rs); // TODO this is added becasue in Johanneson Eq. I use ncells

                double beta_shock_rs;
                switch (p_pars->method_shock_vel) {
                    case isameAsBW:
                        beta_shock_rs = EQS::Beta(Gamma);
                        break;
                    case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> computeSynchrotronEmissivityAbsorptionAnalytic shock velocity
                        beta_shock_rs = EQS::Beta(GammaShock_rs);//us / sqrt(1. + us * us);
                        break;
                }
                double ashock_rs = (1.0 - mu * beta_shock_rs); // shock velocity beaming factor
                dr_rs /= ashock_rs; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
                dr_tau_rs /= ashock_rs;
                double dtau_rs = optical_depth(abs_lab_rs, dr_tau_rs, mu, beta_shock_rs);
                double intensity_rs = computeIntensity(em_lab_rs, dtau_rs,p_pars->p_mphys->method_tau);
                if (intensity_rs < 0 || !std::isfinite(intensity_rs)) {
                    (*p_pars->p_log)(LOG_ERR, AT) << "intensity_rs = " << intensity_rs << "\n";
                    exit(1);
                }
                double flux_dens_rs = intensity_rs * r * r * dr_rs;
                if (flux_dens_rs < 0 || !std::isfinite(flux_dens_rs)) {
                    (*p_pars->p_log)(LOG_ERR, AT) << "flux_dens_rs = " << flux_dens_rs << "\n";
                    exit(1);
                }
                flux_dens += flux_dens_rs * (1.0 + p_pars->z) / (2.0 * p_pars->d_l *
                                                                 p_pars->d_l);; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
            }
        }
    }
    else{
        double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
        double GammaSh = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGammaFsh]);
        double rho2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::irho2]);
        double m2 = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iM2]);
        double frac = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iacc_frac]);
        double thick = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithickness]);
        theta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::itheta]);
        ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
        double B = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iB]);
        double gm = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igm]);
        double gM = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igM]);
        double gc = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::igc]);
        double Theta = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iTheta]);
        double z_cool = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iz_cool]);
        double tburst = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itburst]);
        double tt = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::itt]);
//            double cs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iCSCBM]);

        if ((!std::isfinite(gm))||(!std::isfinite(B))||(!std::isfinite(m2))){

            (*p_pars->p_log)(LOG_ERR,AT)
                    <<"[ish="<<p_pars->ishell<<", "<<"il="<<p_pars->ilayer<<"] "
                    <<" nanss {"
                    <<" ia="<<ia<<" ib="<<ib<<" ta="<<ta<<" tb="<<tb<<" r="<<r<<" mu="<<mu
                    <<" nu_obs="<<nu_obs<<" t_obs="<<t_obs
                    <<" phi="<<phi<<" theta="<<theta<<" ctheta="<<ctheta<<" flux_dens="<<flux_dens
                    << " | " << " Gamma="<<Gamma<<" rho2="<<rho2<<" m2="<<m2<<" B="<<B
                    <<"\n";
            ;
            exit(1);
        }

        double dFnu = 0.;
        if ((m_data[BW::Q::iB][ia] == 0.) || (m_data[BW::Q::iB][ib] == 0.)){
            dFnu = 0.;
        }
        else{
#if 0
            if ((Gamma < 1. || !std::isfinite(Gamma))
                        || (gm < 0.) || (!std::isfinite(gm))
                        || (gM < 0.) || (gc < 0.) || (B < 0.) || (!std::isfinite(B))
                        || (theta <= 0.)
                        || (rho2 < 0.) || (!std::isfinite(rho2))
                        || (thick <= 0.) || (GammaSh <= 1.)) {
                        (*p_log)(LOG_ERR,AT)<< " Error in interpolation to EATS surface: \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " GammaSh = " << GammaSh << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n";
//                    exit(1);
                    }
                    if ((B != 0.) && (!std::isfinite(rho2))) {
                        (*p_log)(LOG_ERR,AT)<< " B!=0 and rho2 is NAN \n"
                                            << " Error in data \n"
                                            << " R = " << r << "\n"
                                            << " Gamma = " << Gamma << "\n"
                                            << " gm = " << gm << "\n"
                                            << " gM = " << gM << "\n"
                                            << " gc = " << gc << "\n"
                                            << " B = " << B << "\n"
                                            << " Theta = " << Theta << "\n"
                                            << " z_cool = " << z_cool << "\n"
                                            << " theta = " << theta << "\n"
                                            << " rho2 = " << rho2 << "\n"
                                            << " thick = " << thick << "\n"
                                            << " t_obs = " << t_obs << "\n"
                                            << " Exiting...\n";
                        exit(1);
                    }
#endif
            double thick_tau = EQS::shock_delta(r,GammaSh); // TODO this is added becasue in Johanneson Eq. I use ncells
            flux_dens = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2, frac, B, gm, gM, gc,
                                                       Theta, z_cool, t_obs, mu,
                                                       r, thick,  thick_tau, nu_obs, p_pars);
        }

        if (p_pars->m_type == BW_TYPES::iFSRS) {
            (*p_pars->p_log)(LOG_ERR,AT) << " not supported her...";
            exit(1);
        }

        flux_dens *= (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
    }
}
#endif

/// use evaluated local emissivity and absoprtion to get observed flux density for Adaptive EATS and blastwave structure
static void fluxDensAdaptiveWithObs(
        double & flux_dens, double & r, double & ctheta, double theta, double phi,
        size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs,
        void * params){

    auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
    auto & p_syna = p_pars->p_mphys;//->getAnSynch();
    auto & Dt = p_pars->m_data;
    auto & times = Dt[BW::Q::itburst];
    auto & r_arr = Dt[BW::Q::iR];

    r = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iR]);
    if (!std::isfinite(r)) {
        (*p_pars->p_log)(LOG_ERR,AT) << " R is NAN in integrand for microphysics"
                                     << " t_e=" << t_e
                                     << " m_data[BW::Q::iR][ia]=" << Dt[BW::Q::iR][ia]
                                     << " m_data[BW::Q::iR][ib]=" << Dt[BW::Q::iR][ib]
                                     << " m_data[BW::Q::itburst][ia]=" << Dt[BW::Q::itburst][ia]
                                     << " m_data[BW::Q::itburst][ib]=" << Dt[BW::Q::itburst][ib]
                                     << "\n";
        // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << m_data[BW::Q::iR] << "\n";
//            std::cout << " Gamma= " << m_data[BW::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
//                return 0.;
        exit(1);
    }

    double rho = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::irho]);
    double Gamma = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iGamma]);
    double GammaSh = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iGammaFsh]);
    double beta = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::ibeta]);
    double U_p = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM2));
    theta = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::itheta]);
    double nprime = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::irho2])/CGS::mp;
    double nprotons = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iM2])/CGS::mp;
    double frac = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iacc_frac]);
    double thick = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::ithickness]);
    double gm = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::igm]);
    double gM = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::igM]);
    double gc = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::igc]);
    double B = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iB]);
    double Theta = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iTheta]);
    double z_cool = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iz_cool]);
    double rsh = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iRsh]);
    double acc_fac = interpSegLog(ia, ib, t_e, times, Dt[BW::Q::iacc_frac]);

    if (rho < 0. || Gamma < 1. || !std::isfinite(Gamma)
        || U_p < 0. || theta <= 0. || nprime < 0. || thick <= 0.) {
        std::cerr << " wrong value in interpolation to EATS surface  \n"
                  << " r = " << r << "\n"
                  << " rho = " << rho << "\n"
                  << " Gamma = " << Gamma << "\n"
                  << " U_p = " << U_p << "\n"
                  << " theta = " << theta << "\n"
                  << " nprime = " << nprime << "\n"
                  << " thick = " << thick << "\n"
                  << " t_e = " << t_e << "\n"
                  << " Exiting...\n";
        exit(1);
    }

    flux_dens = shock_synchrotron_flux_density(
            EjectaID2::STUCT_TYPE::iadaptive, theta,
            Gamma, GammaSh, acc_fac, nprotons, nprime,
            frac, B, gm, gM, gc, Theta, z_cool,
            t_e, mu, r, rsh, thick, thick, nu_obs, params);
//            flux_dens*=(p_pars->d_l*p_pars->d_l*2.);
#if 0
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
                /// computeSynchrotronEmissivityAbsorptionAnalytic the 'thickness' of the shock (emitting region)
    //            double dr_rs = thick3;

    //          // TODO check if for the shock velocity gamma43 has to be used!!! I think it is gam43! See gammaAdi calc.
                dFnu_rs = Observables::shock_synchrotron_flux_density(
                        Gamma, gamma43, rho3, U_e3, t_e, mu, R, thick3, thick3, params );
            }
            dFnu += dFnu_rs;
#endif

    if (p_pars->m_type == BW_TYPES::iFSRS) {
        (*p_pars->p_log)(LOG_ERR,AT) << " BW type FSRS is not supported for mode of "
                                        "computing observer radiation without comoving radiation first\n";
        exit(1);
    }

//            p_pars->o_gam = Gamma;
//            p_pars->o_r = R;
//            p_pars->o_mu = mu;
//            p_pars->o_flux = dFnu;
//            p_pars->o_theta_j = theta;

}


class BlastWaveRadiation : public BlastWaveBase{

    std::unique_ptr<logger> p_log;
    std::unique_ptr<EATS> p_eats_fs = nullptr;
    std::unique_ptr<EATS> p_eats_opt_depth = nullptr;


    //    std::unique_ptr<Source> m_source = nullptr;
//    std::unique_ptr<ShockMicrophysicsNumeric> m_mfphys = nullptr;
public:
    BlastWaveRadiation(Vector & tb_arr, size_t ishell, size_t ilayer,size_t n_substeps,
                       BW_TYPES type, int loglevel)
            : BlastWaveBase(tb_arr,ishell,ilayer,n_substeps, type, loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW_radiation");

        /// init radiation
        p_pars->p_mphys = std::make_unique<ElectronAndRadiation>(loglevel, false);
        p_pars->p_mphys_rs = std::make_unique<ElectronAndRadiation>(loglevel, true);

        /// EATS integrator for forward shock
//        p_eats_pars = new EatsPars(mD); /// for static methods (data link)
        p_eats_fs = std::make_unique<EATS>(m_data[BW::Q::itburst],
                                           m_data[BW::Q::itt], m_data[BW::Q::iR], m_data[BW::Q::itheta],
                                           m_data[BW::Q::iGamma], m_data[BW::Q::ibeta],
//                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
//        p_eats_opt_depth = std::make_unique<EATS>(mD[BW::Q::itburst],
//                                           mD[BW::Q::itt], mD[BW::Q::iR], mD[BW::Q::itheta],
//                                           mD[BW::Q::iGamma],mD[BW::Q::ibeta],
//                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
//                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
//        if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
//            p_eats_fs->setFluxFunc(fluxDensPieceWiseWithComov);
//            p_eats_fs->setFluxFuncA(fluxDensAdaptiveWithComov);
//        }
//        else {
//            p_eats_fs->setFluxFunc(fluxDensPieceWiseWithObs);
//            p_eats_fs->setFluxFuncA(fluxDensAdaptiveWithObs);
//        }

        /// initialize the source
//        m_source = std::make_unique<Source>();
//        m_mfphys = std::make_unique<ShockMicrophysicsNumeric>(loglevel);

    }

    void setParams(std::unique_ptr<EjectaID2> & id,
                   StrDbMap & pars, StrStrMap & opts, size_t ilayer, size_t ii_eq){

        setBaseParams(id, pars, opts, ilayer, ii_eq);

        size_t & ish = p_pars->ishell;
        size_t & il = p_pars->ilayer;

        p_pars->theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        p_pars->d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        p_pars->z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
        p_pars->theta_c = id->get(ish,il,EjectaID2::Q::itheta_c);
        p_pars->theta_c_l = id->get(ish,il,EjectaID2::Q::itheta_c_l);
        p_pars->theta_c_h = id->get(ish,il,EjectaID2::Q::itheta_c_h);
        p_pars->theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);

        p_eats_fs->setEatsPars(
                pars,opts,id->nlayers, id->ncells,id->get(ish,il,EjectaID2::Q::ictheta),
                id->get(ish,il,EjectaID2::Q::itheta_c_l),
                id->get(ish,il,EjectaID2::Q::itheta_c_h),
                id->theta_wing,
                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));


        /// Set Electron And Radiation Model
        p_pars->p_mphys->setPars(pars, opts, p_pars->nr, p_pars->theta_c_h,
                                 EQS::initTComov(p_pars->R0, p_pars->beta0, p_pars->Gamma0));
        if (p_pars->do_rs_radiation)
            p_pars->p_mphys_rs->setPars(pars, opts, p_pars->nr, p_pars->theta_c_h,
                                        EQS::initTComov(p_pars->R0, p_pars->beta0, p_pars->Gamma0));

        /// Set EATS functions (interpolator functions)
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
            p_eats_fs->fluxFuncPW = fluxDensPieceWiseWithComov;
            p_eats_fs->fluxFuncA = fluxDensAdaptiveWithComov;
        }
        else {
            p_eats_fs->fluxFuncPW = fluxDensPieceWiseWithObs;
            p_eats_fs->fluxFuncA = fluxDensAdaptiveWithObs;
        }

    }

    std::unique_ptr<EATS> & getFsEATS(){ return p_eats_fs; }


    bool isBlastWaveValidForElectronCalc(size_t it){

        /// check if BW is evolved
        if (p_pars->R0 < 0)
            return false;

        /// check memory allocation
        if (it > m_data[BW::Q::itburst].size()){
            (*p_log)(LOG_ERR,AT) << " it=" << it << " out of range=m_data[BW::Q::itburst].size()="
                                 << m_data[BW::Q::itburst].size() << " mtburst.size()=" << p_pars->nr << "\n";
            exit(1);
        }

        auto & p_syna = p_pars->p_mphys;

        /// update the i_end_r with respect to minimum for radiation
        if ((m_data[BW::Q::ibeta][it] < p_syna->beta_min)){
            if (p_pars->i_end_r > it-1)
                p_pars->i_end_r = it-1;
            return false;
        }

        if (p_pars->i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "[ish="<<p_pars->ishell<<" il="<<p_pars->ilayer<<"] "
                                 <<"beta0="<<p_pars->beta0<< " i_end_r = 0"<<"\n";
            exit(1);
        }

        double Gamma_ = m_data[BW::Q::iGamma][it]; // TODO should be GammaSh

        /// check if the flow is subsonic # TODO this may not work anymore
        if (p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            if ((m_data[BW::Q::iCSCBM][it] >= m_data[BW::Q::ibeta][it])) {
                m_data[BW::Q::iB][it] = 0.;
                if (p_pars->i_end_r > it - 1)
                    p_pars->i_end_r = it - 1;
                return false;
            }
            if (m_data[BW::Q::iGammaREL][it] > 0.)
                Gamma_ = m_data[BW::Q::iGammaREL][it];
            else
                Gamma_ = m_data[BW::Q::iGamma][it];
        }

        // no ISM case -> no electron acceleration
        if ((m_data[BW::Q::iU_p][it] <= 0) || (m_data[BW::Q::irho2][it] <= 0)){
            if (p_pars->i0_failed_elecctrons == 0)
                p_pars->i0_failed_elecctrons = it;
            p_pars->n_fialed_electrons += 1;
//            std::cerr << AT << " at it="<<it<<" U_e="<<m_data[Q::iU_e][it]<<" and rho2="<<m_data[Q::irho2][it]<<" skipping electron calc.\n";
            return false;
        }

//        /// NOTE here we implicitely allow only comoving
//        if (p_pars->m_method_rad != METHODS_RAD::icomovspec || p_pars->R0 < 0)
//            return false;

        /// -- check if there are any data first
        double beta_ = beta_ = m_data[BW::Q::ibeta][it]; // TODO this should be beta_sh?

        if (p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            if (m_data[BW::Q::iGammaREL][it] > 0.)
                beta_ = EQS::BetaFromGamma(m_data[BW::Q::iGammaREL][it]);
            else
                beta_ = m_data[BW::Q::ibeta][it];
        }

        /// if BW is not evolved to this 'it' or velocity is smaller than minimum
        if ((m_data[BW::Q::iR][it] < 1) || (beta_ < p_syna->beta_min))
            return false;

        return true;
    }

    void storeShockPropertiesAndElectronDistributionLimits(size_t it, ElectronAndRadiaionBase * p_rad){

//        auto & p_rad = p_pars->p_rad;

        /// Save Results
        m_data[BW::Q::iB][it]      = p_rad->B;
        m_data[BW::Q::igm][it]     = p_rad->gamma_min;
        m_data[BW::Q::igM][it]     = p_rad->gamma_max;
        m_data[BW::Q::igc][it]     = p_rad->gamma_c;
        m_data[BW::Q::iTheta][it]  = p_rad->Theta;
        m_data[BW::Q::iz_cool][it] = p_rad->z_cool;
        m_data[BW::Q::inprime][it] = p_rad->n_prime;
        m_data[BW::Q::iacc_frac][it] = p_rad->accel_frac;

        if (m_data[BW::Q::igm][it] == 0){
            (*p_log)(LOG_ERR,AT)<< " in evolved blast wave, found gm = 0" << "\n";
            exit(1);
        }
        /// Check results
        if ((!std::isfinite(m_data[BW::Q::igm][it] )) || (m_data[BW::Q::iB][it] < 0.)) {
            (*p_log)(LOG_ERR,AT) << " Wrong value at i=" << it << " tb=" << m_data[BW::Q::itburst][it]
                                 << " iR=" << m_data[BW::Q::iR][it] << "\n"
                                 << " iGamma=" << m_data[BW::Q::iGamma][it] << "\n"
                                 << " ibeta=" << m_data[BW::Q::ibeta][it] << "\n"
                                 << " iM2=" << m_data[BW::Q::iM2][it] << "\n"
                                 << " iEint2=" << m_data[BW::Q::iEint2][it] << "\n"
                                 << " iU_p=" << m_data[BW::Q::iU_p][it] << "\n"
                                 << " irho=" << m_data[BW::Q::irho][it] << "\n"
                                 << " irho2=" << m_data[BW::Q::irho2][it] << "\n"
                                 << " iB=" << m_data[BW::Q::iB][it] << "\n"
                                 << " igm=" << m_data[BW::Q::igm][it] << "\n"
                                 << " igM=" << m_data[BW::Q::igM][it] << "\n"
                                 << " igc=" << m_data[BW::Q::igc][it] << "\n"
                                 << " iTheta=" << m_data[BW::Q::iTheta][it] << "\n"
                                 << " iz_cool=" << m_data[BW::Q::iz_cool][it] << "\n"
                                 << " inprime=" << m_data[BW::Q::inprime][it]
                                 << "\n";
            exit(1);
        }
    }

    void storeReverseShockPropertiesAndElectronDistributionLimits(size_t it, ElectronAndRadiaionBase * p_rad_rs){

//        auto & p_rad_rs = p_pars->p_rad_rs;

        // adding
        m_data[BW::Q::iB3][it] = p_rad_rs->B;
        m_data[BW::Q::igm_rs][it] = p_rad_rs->gamma_min;
        m_data[BW::Q::igM_rs][it] = p_rad_rs->gamma_max;
        m_data[BW::Q::igc_rs][it] = p_rad_rs->gamma_c;
        m_data[BW::Q::iTheta_rs][it] = p_rad_rs->Theta;
        m_data[BW::Q::iz_cool_rs][it] = p_rad_rs->z_cool;
        m_data[BW::Q::inprime_rs][it] = p_rad_rs->n_prime;
        m_data[BW::Q::iacc_frac_rs][it] = p_rad_rs->accel_frac;

        if (m_data[BW::Q::igm_rs][it] == 0){
            (*p_log)(LOG_ERR,AT)<< " in evolved blast wave, found gm = 0" << "\n";
            exit(1);
        }
        /// Check results
        if ((!std::isfinite(m_data[BW::Q::igm][it] )) || (m_data[BW::Q::iB][it] < 0.)) {
            (*p_log)(LOG_ERR,AT) << " Wrong value at i=" << it << " tb=" << m_data[BW::Q::itburst][it]
                                 << " iR=" << m_data[BW::Q::iR][it] << "\n"
                                 << " iGamma=" << m_data[BW::Q::iGamma][it] << "\n"
                                 << " iGamma43=" << m_data[BW::Q::iGamma43][it] << "\n"
                                 << " ibeta=" << m_data[BW::Q::ibeta][it] << "\n"
                                 << " iM2=" << m_data[BW::Q::iM2][it] << "\n"
                                 << " iM3=" << m_data[BW::Q::iM3][it] << "\n"
                                 << " iEint2=" << m_data[BW::Q::iEint2][it] << "\n"
                                 << " iEint3=" << m_data[BW::Q::iEint3][it] << "\n"
                                 << " iU_p=" << m_data[BW::Q::iU_p][it] << "\n"
                                 << " iU_p3=" << m_data[BW::Q::iU_p3][it] << "\n"
                                 << " irho=" << m_data[BW::Q::irho][it] << "\n"
                                 << " irho2=" << m_data[BW::Q::irho2][it] << "\n"
                                 << " irho3=" << m_data[BW::Q::irho3][it] << "\n"
                                 << " iB=" << m_data[BW::Q::iB][it] << "\n"
                                 << " iB3=" << m_data[BW::Q::iB3][it] << "\n"
                                 << " igm=" << m_data[BW::Q::igm][it] << "\n"
                                 << " igm_rs=" << m_data[BW::Q::igm_rs][it] << "\n"
                                 << " igM=" << m_data[BW::Q::igM][it] << "\n"
                                 << " igc=" << m_data[BW::Q::igc][it] << "\n"
                                 << " igc_rs=" << m_data[BW::Q::igc_rs][it] << "\n"
                                 << " iTheta=" << m_data[BW::Q::iTheta][it] << "\n"
                                 << " iz_cool=" << m_data[BW::Q::iz_cool][it] << "\n"
                                 << " inprime=" << m_data[BW::Q::inprime][it] << "\n"
                                 << " inprime3=" << m_data[BW::Q::inprime_rs][it]
                                 << "\n";
            exit(1);
        }
    }

    bool considerReverseShock(size_t it){
        /// compute electron distribution in reverse shock TODO remove "comov" req. add method to EATS integrator
        if ( ( p_pars->m_type == BW_TYPES::iFSRS
            && p_pars->m_method_rad == METHODS_RAD::icomovspec)
            && (p_pars->do_rs_radiation
                && (it > 0)
                && (m_data[BW::Q::ithickness_rs][it - 1] > 0) // for numerical electron evolution (otherwise it fails)
            && (m_data[BW::Q::iGammaRsh][it] > 0)
                && (m_data[BW::Q::iU_p3][it] > 0)) )
                return true;
        else
            return false;

    }

    void evolveElectronDistAndComputeRadiation(){

        (*p_log)(LOG_INFO,AT)
            << " computing comoving spectrum [ish="<<p_pars->ishell<<" il="<<p_pars->ilayer<<"] \n";

        auto & p_mphys = p_pars->p_mphys;
        auto & p_mphys_rs = p_pars->p_mphys_rs;

        /// evolve
        for (size_t it = 0; it < p_pars->nr; it++) {

            /// check if to consider this timestep
            if ( not isBlastWaveValidForElectronCalc( it ) )
                return;

            /// compute electron distribution in reverse shock
            if ( considerReverseShock(it) ){

                p_mphys_rs->updateSockProperties(//m_data[BW::Q::iR][it],
                        //m_data[BW::Q::ithickness_rs][it],
                        m_data[BW::Q::iU_p3][it],
                        m_data[BW::Q::iGamma][it],
//                               m_data[Q::iGamma][it],
                        m_data[BW::Q::iGammaRsh][it],
                        m_data[BW::Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
                        m_data[BW::Q::itcomov][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               m_data[Q::itburst][it], // emission time (TT)
                        m_data[BW::Q::irho3][it] / CGS::mp,
                        m_data[BW::Q::iM3][it] / CGS::mp
                );

                p_mphys_rs->evaluateElectronDistributionAnalytic();

                storeReverseShockPropertiesAndElectronDistributionLimits(
                        it,const_cast<ElectronAndRadiaionBase *>(p_mphys_rs->getThis()));

                /// compute comoving spectra
                if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
                    if (p_mphys_rs->m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt)
                        p_mphys_rs->computeSynchrotronSpectrumAnalytic(it
//                                                                       m_data[BW::Q::iR][it],
//                                                                       m_data[BW::Q::ithickness_rs][it]
                                                                       );
                    else{
                        if (not p_mphys_rs->is_distribution_initialized)
//                            p_mphys_rs->is_distribution_initialized = true;
                            p_mphys_rs->initializeElectronDistribution(
                                    m_data[BW::Q::itcomov][it],m_data[BW::Q::iM3][it]);
                        else
                            p_mphys_rs->evaluateElectronDistributionNumericMixed(
                                    m_data[BW::Q::itcomov][it - 1],
                                    m_data[BW::Q::itcomov][it],
                                    m_data[BW::Q::iM3][it - 1],
                                    m_data[BW::Q::iM3][it],
                                    m_data[BW::Q::irho3][it - 1],
                                    m_data[BW::Q::irho3][it],
                                    m_data[BW::Q::iR][it - 1],
                                    m_data[BW::Q::iR][it],
                                    m_data[BW::Q::ithickness_rs][it - 1] *
                                    m_data[BW::Q::iGammaRsh][it - 1],  // comoving shock thickness
                                    m_data[BW::Q::ithickness_rs][it] * m_data[BW::Q::iGammaRsh][it - 1]);
                        if (p_mphys_rs->is_distribution_initialized)
                            p_mphys_rs->storeSynchrotronSpectrumNumericMixed(it);
                    }
                }
            }


            /// compute electron distribution in forward shock
            p_mphys->updateSockProperties(//m_data[BW::Q::iR][it],
                                          //m_data[BW::Q::ithickness][it],
                                          m_data[BW::Q::iU_p][it],
                                          m_data[BW::Q::iGamma][it],
//                               m_data[Q::iGamma][it],
                                          m_data[BW::Q::iGammaFsh][it],
                                          m_data[BW::Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
                                          m_data[BW::Q::itcomov][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               m_data[Q::itburst][it], // emission time (TT)
                                          m_data[BW::Q::irho2][it] / CGS::mp,
//                                          (m_data[BW::Q::iM2][it+1]-m_data[BW::Q::iM2][it]) / CGS::mp
                                          m_data[BW::Q::iM2][it] / CGS::mp//(m_data[BW::Q::iM2][it+1]-m_data[BW::Q::iM2][it]) / CGS::mp
                                          );

            /// compute electron injection function / electron spectrum (analytically)
            p_mphys->evaluateElectronDistributionAnalytic();

            /// store the result in the main storage
            storeShockPropertiesAndElectronDistributionLimits(
                    it, const_cast<ElectronAndRadiaionBase *>(p_mphys->getThis()));

//            std::cout << m_data[BW::Q::ithickness][it-1]<< " " <<  m_data[BW::Q::ithickness][it] << "\n";
            /// compute comoving spectra
            if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
                if (p_mphys->m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt)
                    p_mphys->computeSynchrotronSpectrumAnalytic(it
//                                                                m_data[BW::Q::iR][it],
//                                                                m_data[BW::Q::ithickness][it]
                                                                );
                else {
                    if (not p_mphys->is_distribution_initialized) // Initialize electron distribution analytically
                        p_mphys->initializeElectronDistribution(
                                m_data[BW::Q::itcomov][it],m_data[BW::Q::iM2][it]);
                    else // Evolve electron distribution numerically (or analytically)
                        p_mphys->evaluateElectronDistributionNumericMixed(
                                m_data[BW::Q::itcomov][it - 1],
                                m_data[BW::Q::itcomov][it],
                                m_data[BW::Q::iM2][it - 1],
                                m_data[BW::Q::iM2][it],
                                m_data[BW::Q::irho2][it - 1],
                                m_data[BW::Q::irho2][it],
                                m_data[BW::Q::iR][it - 1],
                                m_data[BW::Q::iR][it],
                                m_data[BW::Q::ithickness][it - 1] * m_data[BW::Q::iGammaFsh][it -
                                                                                             1], // comoving shock thickness dR' = dR * Gamma... (Granot + 1999)
                                m_data[BW::Q::ithickness][it] * m_data[BW::Q::iGammaFsh][it]);
                    if (p_mphys->is_distribution_initialized)
                        p_mphys->storeSynchrotronSpectrumNumericMixed(it);
                }
            }
        }

        /// check if electron spectrum failed for any reason
        if ((m_data[BW::Q::iR][0] > 0) && (p_pars->n_fialed_electrons == p_pars->nr) && (!p_pars->end_evolution)){
            (*p_log)(LOG_ERR,AT)
                    << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                    << " Electron calculation failed for all iterations. Exiting...\n";
            exit(1);
        }
        else if ((m_data[BW::Q::iR][0] > 0) && (p_pars->n_fialed_electrons > 0) && (p_pars->n_fialed_electrons < p_pars->nr)){
            (*p_log)(LOG_ERR,AT)
                    << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                    <<" Electron calculation failed for n=" << p_pars->n_fialed_electrons
                    << " iterations starting from it=" << p_pars->i0_failed_elecctrons<<"\n";
        }

        if (p_pars->loglevel >= LOG_WARN) {
            double _gm_max = *std::max_element(m_data[BW::Q::igm].begin(), m_data[BW::Q::igm].end());
            double _gm_min = *std::min_element(m_data[BW::Q::igm].begin(), m_data[BW::Q::igm].end());
            if (_gm_max == _gm_min) {
                (*p_log)(LOG_WARN, AT) << " min(gm) == max(gm) = " << _gm_min << "\n";
                exit(1);
            }
        }
    }

#if 0
    static void optDepthPW(double & tau_Compton, double & tau_BH, double & tau_bf, double & r, double & ctheta,
                           size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                           Vector & ttobs, void * params){

        auto * p_pars = (struct PWNPars *) params;
        auto & m_data = p_pars->m_data;
        if (p_pars->i_end_r==0)
            return;

        /// interpolate ejecta properties for a given time
        double delta_ej = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iEJdelta]);
        double rho_ej = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iEJrho]);
        double Gamma = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::iGamma]);
        double betaSh = interpSegLog(ia, ib, t_obs, ttobs, m_data[BW::Q::ibeta]);
        ctheta = interpSegLin(ia, ib, t_obs, ttobs, m_data[BW::Q::ictheta]);

        /// account for relativistic motion of the shell
        double a = 1.0 - betaSh * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
        double nu_erg = nu_obs*4.1356655385381E-15*CGS::EV_TO_ERG;

        /// evaluate optical depths
        double e_gamma = nu_erg; //! TODO CHECK!;
        double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/p_pars->mu_e/CGS::M_PRO;
        tau_Compton = rho_ej*delta_ej*Kcomp;
        /// optical depth of BH pair production
//        double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/CGS::M_PRO/r_ej/r_ej;
        double KBH = (1.0+p_pars->Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/p_pars->mu_e/CGS::M_PRO;
        tau_BH = rho_ej*delta_ej*KBH;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
//        double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma, Z_eff, opacitymode)/r_ej/r_ej;
        double Kbf = (1.0-p_pars->albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, p_pars->Z_eff, p_pars->opacitymode);
        tau_bf = rho_ej*delta_ej*Kbf;

    }
#endif

    bool evalEATSindexes(size_t &ia, size_t &ib, double t_obs, double z, //size_t m_i_end_r,
                         double ctheta, double cphi, double theta_obs,
                         double (*obs_angle_func)( const double &, const double &, const double & )){

//        auto * p_pars = (struct EatsPars *) params; // removing EATS_pars for simplicity
        auto & m_mu =  p_pars->m_mu;
        auto & ttobs =  p_pars->ttobs;
        size_t m_i_end_r= p_pars->i_end_r;
        if (m_mu.size() < 1) {
            m_mu.resize(p_pars->i_end_r);
            ttobs.resize(p_pars->i_end_r);
        }
        /// evaluate mu array
        for (size_t i = 0; i < m_i_end_r; i++)
            m_mu[i] = (m_data[BW::Q::itburst][i] - t_obs / (1.0 + z) ) / m_data[BW::Q::iR][i] * CGS::c;
        double mu = obs_angle_func(ctheta, cphi, theta_obs);
        for (size_t i_ = 0; i_ < m_i_end_r; i_++) {
            ttobs[i_] = m_data[BW::Q::itt][i_] + m_data[BW::Q::iR][i_] / CGS::c * (1.0 - mu);
        }
        if (t_obs < ttobs[0]) {
            std::cerr << AT << "t_obs=" << t_obs << " < ttobs_arr[0]=" << ttobs[0] << ". Extend ttobs.\n";
            return false;
        }
        if ((t_obs > ttobs[m_i_end_r - 1])) {
            if (m_i_end_r == m_mu.size()-1)
                std::cerr <<AT << " t_obs="<<t_obs<<" > ttobs_arr[i_end="<<m_i_end_r - 1<<"]="
                          <<ttobs[m_i_end_r - 1]<<". Extend ttobs.\n";
//
//            << " time grid ends too early. "
//                                  << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                  << " while requested obs.time=" << t_obs
//                                  << " extend the grid to later time or request tobs at earlier times\n";
//                    std::cout << ttobs << std::endl;
//            exit(1);
            return false;
        }
        /// locate closest evolution points to the requested obs. time
        ia = findIndex(t_obs, ttobs, ttobs.size());
        if (ia >= m_i_end_r - 1)
            return false; // ??
        ib = ia + 1;
        return true;
    }

    void evalOpticalDepths(double & tau_comp, double & tau_BH, double & tau_bf,
                           size_t ia, size_t ib, double t_obs, double freqprime){

        double rej = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iR]);
        double rho_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iEJrho]);
        double delta_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, m_data[BW::Q::iEJdelta]);

        /// TODO this freq. should be doppler shifted
        double e_gamma = freqprime * 6.62606957030463E-27;; // Hz -> erg //* 6.62606957030463E-27;//* 4.1356655385381E-15 * CGS::EV_TO_ERG;
        double mu_e = p_pars->mu_e;
        double Z_eff = p_pars->Z_eff;
        int opacitymode = p_pars->opacitymode;
        double albd_fac = p_pars->albd_fac;

        /// --- optical depth due to compton scattering
        double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
        double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
        tau_comp=tau_comp_;
        /// optical depth of BH pair production
        double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
        double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
        tau_BH=tau_BH_;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
        double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
        double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
        tau_bf=tau_bf_;
    }

    /// Fraction of PWN emission that thermalises in ejecta
    double facPSRdep(const double rho_ej, const double delta_ej,
                     const double T_ej, const int opacitymode){
        if (p_pars->curr_b_pwn < 0. || !std::isfinite(p_pars->curr_b_pwn)){
            (*p_log)(LOG_ERR,AT)<<" b_pwn is not set\n";
            exit(0);
        }

        double e_gamma_min;
        double e_gamma_max_tmp;
        double e_gamma_gamma_ani_tmp;
        double e_gamma_syn_b_tmp;

        e_gamma_min = 1.0*EV_TO_ERG; // eV
        e_gamma_max_tmp = PWNradiationMurase::e_gamma_max(p_pars->curr_b_pwn);
        e_gamma_gamma_ani_tmp = PWNradiationMurase::e_gamma_gamma_ani(T_ej);
        e_gamma_syn_b_tmp = PWNradiationMurase::e_gamma_syn_b(
                p_pars->curr_b_pwn,p_pars->gamma_b_w);

        if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
            e_gamma_max_tmp = e_gamma_gamma_ani_tmp;

//        const int i_max = 1000;
//        double e_tmp = e_gamma_min;
//        double del_ln_e = 0.0;

        const double albd_fac = p_pars->albd_fac;
//        int opacitymode = 0;

//        tmp(rho_ej,delta_ej,T_ej,p_pars->albd_fac,opacitymode,e_gamma_max_tmp,e_gamma_syn_b_tmp, e_gamma_min);
//        std::cout << " 11 \n";
        double frac_psr_dep_tmp = 0.0; int i = 0;
//        std::cout << " frac_psr_dep_.size(0" << frac_psr_dep_.size()<<"\n";
//        for (size_t i = 0; i <= p_pars->iterations; i++)
//            frac_psr_dep_tmp += frac_psr_dep_[i];
        int nthreads = 6;
        ///
        if (e_gamma_max_tmp > e_gamma_syn_b_tmp){
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(iters+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=iters;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(
                        e_tmp,p_pars->curr_b_pwn,p_pars->gamma_b_w,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==iters)
                    frac_psr_dep_tmp_tmp = (1.0/2.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i] = frac_psr_dep_tmp_tmp;
            }
        }
        else {
            e_gamma_max_tmp = e_gamma_syn_b_tmp;
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(iters+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=iters;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(
                        e_tmp,p_pars->curr_b_pwn,p_pars->gamma_b_w,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==iters)
                    frac_psr_dep_tmp_tmp = (1.0/3.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i]= frac_psr_dep_tmp_tmp;
            }
        }

        for (size_t i = 0; i <= iters; ++i)
            frac_psr_dep_tmp += frac_psr_dep_[i];

        if (frac_psr_dep_tmp > 1.)
            frac_psr_dep_tmp = 1.;
        if (!std::isfinite(frac_psr_dep_tmp)||frac_psr_dep_tmp < 0){
            (*p_log)(LOG_ERR,AT) << "frac_psr_dep_tmp="<<frac_psr_dep_tmp<<"\n";
            exit(1);
        }
        return frac_psr_dep_tmp;
    }

    double getFacPWNdep( const double rho_ej, const double delta_ej, const double T_ej, const double Ye){
        if(!std::isfinite(rho_ej) || rho_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: rho_ej="<<rho_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(delta_ej) || delta_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: delta_ej="<<delta_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(T_ej) || T_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: T_ej="<<T_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(Ye) || Ye < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: Ye="<<Ye<<"\n"; exit(1);
        }

        int opacitymode=0; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
        if (Ye <= 0.2)
            opacitymode = 2; // r-process heavy
        else if ((Ye > 0.2) && (Ye <= 0.3))
            opacitymode = 1; // r-process light
        else if ((Ye > 0.3) && (Ye <= 0.5))
            opacitymode = 0;//0=iron-rich
        else if (Ye > 0.5)
            opacitymode = 3; // CO
        else{
            (*p_log)(LOG_ERR,AT) << " error \n";
            exit(1);
        }
        return facPSRdep(rho_ej, delta_ej, T_ej, opacitymode);
    }

};



#endif //SRC_BLASTWAVE_RADIATION_H
