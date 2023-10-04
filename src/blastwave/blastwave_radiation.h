//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLASTWAVE_RADIATION_H
#define SRC_BLASTWAVE_RADIATION_H

//#include "blastwave_components.h"
//#include "blastwave_base.h"
#include "blastwave_base.h"

class BlastWaveRadiation : public BlastWaveBase{
    std::unique_ptr<logger> p_log;
    std::unique_ptr<EATS> p_eats_fs = nullptr;
    std::unique_ptr<EATS> p_eats_opt_depth = nullptr;
public:
    BlastWaveRadiation(Vector & tb_arr, size_t ishell, size_t ilayer,size_t n_substeps, BW_TYPES type, int loglevel)
            : BlastWaveBase(tb_arr,ishell,ilayer,n_substeps, type, loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW_radiation");
        //        p_rad = std::make_unique<BlastWaveRadiation>(mD, ishell, ilayer, loglevel);
        p_pars->p_syna = std::make_unique<SynchrotronAnalytic>(loglevel, false);
        p_pars->p_syna_rs = std::make_unique<SynchrotronAnalytic>(loglevel, true);
        /// EATS integrator for forward shock
//        p_eats_pars = new EatsPars(mD); /// for static methods (data link)
        p_eats_fs = std::make_unique<EATS>(mD[BW::Q::itburst],
                                           mD[BW::Q::itt], mD[BW::Q::iR], mD[BW::Q::itheta],
                                           mD[BW::Q::iGamma], mD[BW::Q::ibeta],
                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
//        p_eats_opt_depth = std::make_unique<EATS>(mD[BW::Q::itburst],
//                                           mD[BW::Q::itt], mD[BW::Q::iR], mD[BW::Q::itheta],
//                                           mD[BW::Q::iGamma],mD[BW::Q::ibeta],
//                                           p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
//                                           p_pars->i_end_r, ishell, ilayer, loglevel, p_pars);
        p_eats_fs->setFluxFunc(fluxDensPW);
        p_eats_fs->setFluxFuncA(fluxDensA);
    }
    void setParams(std::unique_ptr<EjectaID2> & id, StrDbMap & pars, StrStrMap & opts, size_t ilayer, size_t ii_eq){

        setBaseParams(id, pars, opts, ilayer, ii_eq);

        size_t & ish = p_pars->ishell;
        size_t & il = p_pars->ilayer;

        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e14, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,100, true);//pars.at("nfreq");

        /// allocate space for the comoving spectrum
        p_pars->m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                                  << " freqs="<<p_pars->m_freq_arr.size() << " by radii=" << p_pars->nr << " Spec. grid="
                                  << p_pars->m_freq_arr.size() * p_pars->nr << "\n";
            p_pars->m_synch_em.resize( p_pars->m_freq_arr.size() * p_pars->nr );
            p_pars->m_synch_abs.resize( p_pars->m_freq_arr.size() *p_pars->nr );
            if (p_pars->do_rs){
                p_pars->m_synch_em_rs.resize( p_pars->m_freq_arr.size() * p_pars->nr );
                p_pars->m_synch_abs_rs.resize( p_pars->m_freq_arr.size() * p_pars->nr );
            }
        }

        p_pars->theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        p_pars->d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        p_pars->z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");
        p_pars->theta_c_l = id->get(ish,il,EjectaID2::Q::itheta_c_l);
        p_pars->theta_c_h = id->get(ish,il,EjectaID2::Q::itheta_c_h);
        p_pars->theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);

        p_eats_fs->setEatsPars(
                pars,opts,id->nlayers, id->ncells,id->get(ish,il,EjectaID2::Q::ictheta),
                id->get(ish,il,EjectaID2::Q::itheta_c_l),
                id->get(ish,il,EjectaID2::Q::itheta_c_h),
                id->theta_wing,
                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

        p_pars->p_syna->setPars(pars, opts);
        if (p_pars->do_rs)
            p_pars->p_syna_rs->setPars(pars, opts);
    }
    std::unique_ptr<EATS> & getFsEATS(){ return p_eats_fs; }
    /// ----------------- BLAST WAVE RADIATION ------------

    void addComputeForwardShockMicrophysics(size_t it){

        if (it > mD[BW::Q::itburst].size()){
            (*p_log)(LOG_ERR,AT) << " it=" << it << " out of range=mD[BW::Q::itburst].size()="
                                 << mD[BW::Q::itburst].size() << " mtburst.size()=" << p_pars->nr << "\n";
            exit(1);
        }

        auto & p_syna = p_pars->p_syna;
        /// update the i_end_r with respect to minimum for radiation
        if ((mD[BW::Q::ibeta][it] < p_syna->getPars()->beta_min)){
            if (p_pars->i_end_r > it-1)
                p_pars->i_end_r = it-1;
            return;
        }
        if (p_pars->i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "[ish="<<p_pars->ishell<<" il="<<p_pars->ilayer<<"] "
                                 <<"beta0="<<p_pars->beta0<< " i_end_r = 0"<<"\n";
            exit(1);
        }

        double Gamma_ = mD[BW::Q::iGamma][it]; // TODO should be GammaSh

        /// check if the flow is subsonic # TODO this may not work anymore
        if (p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            if ((mD[BW::Q::iCSCBM][it] >= mD[BW::Q::ibeta][it])) {
                mD[BW::Q::iB][it] = 0.;
                if (p_pars->i_end_r > it - 1)
                    p_pars->i_end_r = it - 1;
                return;
            }
            if (mD[BW::Q::iGammaREL][it] > 0.)
                Gamma_ = mD[BW::Q::iGammaREL][it];
            else
                Gamma_ = mD[BW::Q::iGamma][it];
        }

        // no ISM case -> no electron acceleration
        if ((mD[BW::Q::iU_p][it] <= 0) || (mD[BW::Q::irho2][it] <= 0)){
            if (p_pars->i0_failed_elecctrons == 0)
                p_pars->i0_failed_elecctrons = it;
            p_pars->n_fialed_electrons += 1;
//            std::cerr << AT << " at it="<<it<<" U_e="<<mD[Q::iU_e][it]<<" and rho2="<<mD[Q::irho2][it]<<" skipping electron calc.\n";
            return;
        }

        p_syna->evaluateElectronDistribution(mD[BW::Q::iU_p][it],
                                             mD[BW::Q::iGamma][it],
//                               mD[Q::iGamma][it],
                                             mD[BW::Q::iGammaFsh][it],
                                             mD[BW::Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               mD[Q::itburst][it], // emission time (TT)
                                             mD[BW::Q::irho2][it] / CGS::mp);
        /// Save Results
        mD[BW::Q::iB][it]      = p_syna->getPars()->B;
        mD[BW::Q::igm][it]     = p_syna->getPars()->gamma_min;
        mD[BW::Q::igM][it]     = p_syna->getPars()->gamma_max;
        mD[BW::Q::igc][it]     = p_syna->getPars()->gamma_c;
        mD[BW::Q::iTheta][it]  = p_syna->getPars()->Theta;
        mD[BW::Q::iz_cool][it] = p_syna->getPars()->z_cool;
        mD[BW::Q::inprime][it] = p_syna->getPars()->n_prime;
        mD[BW::Q::iacc_frac][it] = p_syna->getPars()->accel_frac;

        /// Check results
        if ((!std::isfinite(mD[BW::Q::igm][it] )) || (mD[BW::Q::iB][it] < 0.)) {
            (*p_log)(LOG_ERR,AT) << " Wrong value at i=" << it << " tb=" << mD[BW::Q::itburst][it]
                                 << " iR=" << mD[BW::Q::iR][it] << "\n"
                                 << " iGamma=" << mD[BW::Q::iGamma][it] << "\n"
                                 << " ibeta=" << mD[BW::Q::ibeta][it] << "\n"
                                 << " iM2=" << mD[BW::Q::iM2][it] << "\n"
                                 << " iEint2=" << mD[BW::Q::iEint2][it] << "\n"
                                 << " iU_p=" << mD[BW::Q::iU_p][it] << "\n"
                                 << " irho=" << mD[BW::Q::irho][it] << "\n"
                                 << " irho2=" << mD[BW::Q::irho2][it] << "\n"
                                 << " iB=" << mD[BW::Q::iB][it] << "\n"
                                 << " igm=" << mD[BW::Q::igm][it] << "\n"
                                 << " igM=" << mD[BW::Q::igM][it] << "\n"
                                 << " igc=" << mD[BW::Q::igc][it] << "\n"
                                 << " iTheta=" << mD[BW::Q::iTheta][it] << "\n"
                                 << " iz_cool=" << mD[BW::Q::iz_cool][it] << "\n"
                                 << " inprime=" << mD[BW::Q::inprime][it]
                                 << "\n";
            exit(1);
        }


        /// compute electron distribution in reverse shock
        if (p_pars->m_type == BW_TYPES::iFSRS) {
            if (p_pars->do_rs && (it > 0) && (mD[BW::Q::iGammaRsh][it] > 0)) {
                auto &p_syna_rs = p_pars->p_syna_rs;
                p_syna_rs->evaluateElectronDistribution(mD[BW::Q::iU_p3][it],
                                                        mD[BW::Q::iGamma][it],
//                               mD[Q::iGamma][it],
                                                        mD[BW::Q::iGammaRsh][it],
                                                        mD[BW::Q::itt][it], // TODO WHICH TIME IS HERE? tbirst? tcomov? observer time (TT)
//                               mD[Q::itburst][it], // emission time (TT)
                                                        mD[BW::Q::irho3][it] / CGS::mp);
                // adding
                mD[BW::Q::iB3][it] = p_syna_rs->getPars()->B;
                mD[BW::Q::igm_rs][it] = p_syna_rs->getPars()->gamma_min;
                mD[BW::Q::igM_rs][it] = p_syna_rs->getPars()->gamma_max;
                mD[BW::Q::igc_rs][it] = p_syna_rs->getPars()->gamma_c;
                mD[BW::Q::iTheta_rs][it] = p_syna_rs->getPars()->Theta;
//        mD[Q::ix][it]      = p_syn->getPars()->x; // IT is not evaluated/ It depends in freq
                mD[BW::Q::iz_cool_rs][it] = p_syna_rs->getPars()->z_cool;
                mD[BW::Q::inprime_rs][it] = p_syna_rs->getPars()->n_prime;
                mD[BW::Q::iacc_frac_rs][it] = p_syna_rs->getPars()->accel_frac;
            }
        }
    }

    void computeForwardShockElectronAnalyticVars(){

        if (p_pars->end_evolution && (p_pars->E0 < 0))
            return;

        for (size_t it = 0; it < p_pars->nr; it++)
            addComputeForwardShockMicrophysics(it);

        /// check if electron spectrum failed for any reason
        if ((mD[BW::Q::iR][0] > 0) && (p_pars->n_fialed_electrons == p_pars->nr) && (!p_pars->end_evolution)){
            (*p_log)(LOG_ERR,AT)
                    << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                    << " Electron calculation failed for all iterations. Exiting...\n";
            exit(1);
        }
        else if ((mD[BW::Q::iR][0] > 0) && (p_pars->n_fialed_electrons > 0) && (p_pars->n_fialed_electrons < p_pars->nr)){
            (*p_log)(LOG_ERR,AT)
                    << "[il="<<p_pars->ilayer<<", ish="<<p_pars->ishell<<"] "
                    <<" Electron calculation failed for n=" << p_pars->n_fialed_electrons
                    << " iterations starting from it=" << p_pars->i0_failed_elecctrons<<"\n";
        }
    }

    /// comoving spectrum
    void computeForwardShockSynchrotronAnalyticSpectrum(){
        if (p_pars->m_method_rad != METHODS_RAD::icomovspec)
            return;

        if (p_pars->R0 < 0)
            return;
        (*p_log)(LOG_INFO,AT) << " computing analytic comoving spectrum "
                                 " [ish="<<p_pars->ishell<<" il="<<p_pars->ilayer<<"] \n";

        if (p_pars->m_freq_arr.size() < 1){
            (*p_log)(LOG_ERR,AT) << " array for comoving spectrum is not initialized Exiting...\n";
            exit(1);
        }
        if (p_pars->m_synch_em.size() < 1){
            (*p_log)(LOG_ERR,AT)<< " array for comoving spectrum frequencies is not initialized Exiting...\n";
            exit(1);
        }
        if (p_pars->do_rs){
            if (p_pars->m_synch_em_rs.size() < 1){
                (*p_log)(LOG_ERR,AT)<< " array for comoving spectrum frequencies for RS is not initialized Exiting...\n";
                exit(1);
            }
        }

        for (size_t it = 0; it < p_pars->nr; ++it) {
            /// exit if the obs. radiation method of choice does not need comoving spectrum
            size_t nfreq = p_pars->m_freq_arr.size();
            /// -- check if there are any data first
            double beta_ = beta_ = mD[BW::Q::ibeta][it]; // TODO this should be beta_sh?
            if (p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
                if (mD[BW::Q::iGammaREL][it] > 0.)
                    beta_ = EQS::Beta(mD[BW::Q::iGammaREL][it]);
                else
                    beta_ = mD[BW::Q::ibeta][it];
            }

            auto & p_syna = p_pars->p_syna;
            /// if BW is not evolved to this 'it' or velocity is smaller than minimum
            if ((mD[BW::Q::iR][it] < 1) || (beta_ < p_syna->getPars()->beta_min))
                return;

            if (mD[BW::Q::igm][it] == 0){
                (*p_log)(LOG_ERR,AT)<< " in evolved blast wave, found gm = 0" << "\n";
                exit(1);
            }

            /// evaluateShycnhrotronSpectrum emissivity and absorption for each frequency
            for (size_t ifreq = 0; ifreq < p_pars->m_freq_arr.size(); ++ifreq){
                /// evaluateShycnhrotronSpectrum all types of emissivities and absoprtions
                p_syna->evaluateShycnhrotronSpectrum(
                        mD[BW::Q::irho2][it] / CGS::mp,//mD[Q::iM2][it] / CGS::mp,
                        mD[BW::Q::iM2][it] / CGS::mp,//mD[Q::iM2][it] / CGS::mp,
                        mD[BW::Q::iacc_frac][it],
                        mD[BW::Q::iB][it],
                        mD[BW::Q::igm][it],
                        mD[BW::Q::igM][it],
                        mD[BW::Q::igc][it],
                        mD[BW::Q::iTheta][it],
                        mD[BW::Q::iz_cool][it],
//                                       mD[Q::irho2][it] / CGS::mp,
                        p_pars->m_freq_arr[ifreq]
                );
                /// add evaluated data to the storage
//            double thick_tau = EQS::shock_delta(mD[Q::iRsh][it],mD[Q::iGammaFsh][it]);
//            p_syna->addIntensity(thick_tau, 0., 1.);
                p_pars->m_synch_em[ifreq + nfreq * it] = p_syna->getPars()->em;
                p_pars->m_synch_abs[ifreq + nfreq * it] = p_syna->getPars()->abs;
            }

            /// Reverse shock
            if (p_pars->m_type == BW_TYPES::iFSRS) {
                if (p_pars->do_rs && (it > 0) && (mD[BW::Q::iGammaRsh][it] > 0) && (mD[BW::Q::iU_p3][it] > 0)) {
                    auto &p_syna_rs = p_pars->p_syna_rs;
                    if (mD[BW::Q::igm_rs][it] == 0) {
                        (*p_log)(LOG_ERR, AT) << " gm_rs = 0" << "\n";
                        exit(1);
                    }
                    if (mD[BW::Q::iB3][it] <= 0) {
                        (*p_log)(LOG_ERR, AT) << " B3 <= 0" << "\n";
                        exit(1);
                    }

                    /// evaluateShycnhrotronSpectrum emissivity and absorption for each frequency
                    for (size_t ifreq = 0; ifreq < p_pars->m_freq_arr.size(); ++ifreq) {
                        /// evaluateShycnhrotronSpectrum all types of emissivities and absoprtions
                        p_syna_rs->evaluateShycnhrotronSpectrum(
                                mD[BW::Q::irho3][it] / CGS::mp,//mD[Q::iM2][it] / CGS::mp,
                                mD[BW::Q::iM3][it] / CGS::mp,//mD[Q::iM2][it] / CGS::mp,
                                mD[BW::Q::iacc_frac][it],
                                mD[BW::Q::iB3][it],
                                mD[BW::Q::igm_rs][it],
                                mD[BW::Q::igM_rs][it],
                                mD[BW::Q::igc_rs][it],
                                mD[BW::Q::iTheta_rs][it],
                                mD[BW::Q::iz_cool_rs][it],
//                                       mD[Q::irho2][it] / CGS::mp,
                                p_pars->m_freq_arr[ifreq]
                        );
                        /// add evaluated data to the storage
//            double thick_tau = EQS::shock_delta(mD[Q::iRsh][it],mD[Q::iGammaFsh][it]);
//            p_syna->addIntensity(thick_tau, 0., 1.);
                        p_pars->m_synch_em_rs[ifreq + nfreq * it] = p_syna_rs->getPars()->em;
                        p_pars->m_synch_abs_rs[ifreq + nfreq * it] = p_syna_rs->getPars()->abs;
                    }

                }
            }
        }
    }

#if 0
    static void optDepthPW(double & tau_Compton, double & tau_BH, double & tau_bf, double & r, double & ctheta,
                           size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                           Vector & ttobs, void * params){

        auto * p_pars = (struct PWNPars *) params;
        auto & mD = p_pars->mD;
        if (p_pars->i_end_r==0)
            return;

        /// interpolate ejecta properties for a given time
        double delta_ej = interpSegLog(ia, ib, t_obs, ttobs, mD[BW::Q::iEJdelta]);
        double rho_ej = interpSegLog(ia, ib, t_obs, ttobs, mD[BW::Q::iEJrho]);
        double Gamma = interpSegLog(ia, ib, t_obs, ttobs, mD[BW::Q::iGamma]);
        double beta = interpSegLog(ia, ib, t_obs, ttobs, mD[BW::Q::ibeta]);
        ctheta = interpSegLin(ia, ib, t_obs, ttobs, mD[BW::Q::ictheta]);

        /// account for relativistic motion of the shell
        double a = 1.0 - beta * mu; // beaming factor
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

    static void fluxDensPW(double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
                           double r, double & ctheta, double theta, double phi,
                           size_t ia, size_t ib, double ta, double tb, double mu,
                           double t_obs, double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params;
        auto & m_data = p_pars->m_data;
        if (p_pars->i_end_r==0)
            return;

        if (p_pars->m_method_rad == METHODS_RAD::icomovspec) {
            Interp2d int_em(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_em);
            Interp2d int_abs(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;

            double Gamma = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma]);
            double beta = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ibeta]);
            // double GammaSh = ( Interp1d(mD[BW::Q::iR], mD[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor
            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z) * nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_pars->m_freq_arr, p_pars->m_freq_arr.size());
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
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = RadiationBase::optical_depth(abs_lab, dr_tau, mu, beta_shock);
            double intensity = RadiationBase::computeIntensity(em_lab, dtau, p_pars->p_syna->getPars()->method_tau);
            flux_dens = (intensity * r * r * dr) * (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
//        flux += flux_dens;
            /// save the result in image
            ctheta = interpSegLin(ia, ib, ta, tb, t_obs, m_data[BW::Q::ictheta]);
            //  double ctheta = ( Interp1d(mD[BW::Q::iR], mD[BW::Q::ictheta] ) ).Interpolate(r, mth );
//        image(Image::iintens, i) =
//                flux_dens / (r * r * std::abs(mu)) * CGS::cgs2mJy; //(obs_flux / (delta_D * delta_D * delta_D)); //* tmp;
//        image(Image::ixr, i) = r * im_xxs(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::iyr, i) = r * im_yys(ctheta, phi_cell, p_pars->theta_obs);
//        image(Image::imu, i) = mu_arr[i];

            if (p_pars->m_type == BW_TYPES::iFSRS) {
                if (p_pars->do_rs &&
                    !(m_data[BW::Q::ithichness_rs][ia] == 0 || m_data[BW::Q::ithichness_rs][ib] == 0)) {
                    Interp2d int_em_rs(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_em_rs);
                    Interp2d int_abs_rs(p_pars->m_freq_arr, m_data[BW::Q::iR], p_pars->m_synch_abs_rs);
                    double em_prime_rs = int_em_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                    double abs_prime_rs = int_abs_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                    double em_lab_rs = em_prime_rs / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                    double abs_lab_rs = abs_prime_rs * delta_D; // conversion of absorption (see vanEerten+2010)

                    double GammaShock_rs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::iGamma43]);
                    double dr_rs = interpSegLog(ia, ib, ta, tb, t_obs, m_data[BW::Q::ithichness_rs]);

                    double dr_tau_rs = EQS::shock_delta(r, GammaShock_rs); // TODO this is added becasue in Johanneson Eq. I use ncells

                    double beta_shock_rs;
                    switch (p_pars->method_shock_vel) {
                        case isameAsBW:
                            beta_shock_rs = EQS::Beta(Gamma);
                            break;
                        case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                            beta_shock_rs = EQS::Beta(GammaShock_rs);//us / sqrt(1. + us * us);
                            break;
                    }
                    double ashock_rs = (1.0 - mu * beta_shock_rs); // shock velocity beaming factor
                    dr_rs /= ashock_rs; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
                    dr_tau_rs /= ashock_rs;
                    double dtau_rs = RadiationBase::optical_depth(abs_lab_rs, dr_tau_rs, mu, beta_shock_rs);
                    double intensity_rs = RadiationBase::computeIntensity(em_lab_rs, dtau_rs,
                                                                          p_pars->p_syna->getPars()->method_tau);
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

    static void fluxDensA(double & flux_dens, double & r, double & ctheta, double theta, double phi,
                          size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;//->getAnSynch();
        auto & Dt = p_pars->m_data;
        auto & tburst = Dt[BW::Q::itburst];
        auto & r_arr = Dt[BW::Q::iR];

        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            Interp2d int_em(p_pars->m_freq_arr, r_arr, p_pars->m_synch_em);
            Interp2d int_abs(p_pars->m_freq_arr, r_arr, p_pars->m_synch_abs);
            Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r = interpSegLog(ia, ib, t_e, tburst, r_arr);
            //  double r = ( Interp1d(ttobs, mD[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                             << " Current R grid us ["
                                             << r_arr[0] << ", "
                                             << r_arr[tburst.size() - 1] << "] "
                                             << "and tburst arr ["
                                             << tburst[0] << ", " << tburst[p_pars->nr - 1]
                                             << "] while the requried obs_time=" << t_obs
                                             << "\n";
//                std::cerr << AT << "\n";
                return;
            }
            double Gamma = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::iGamma]);
            double beta = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::ibeta]);
            // double GammaSh = ( Interp1d(mD[BW::Q::iR], mD[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - beta * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor
            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_pars->m_freq_arr, p_pars->m_freq_arr.size());
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
            double GammaShock = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::iGammaFsh]);
            double dr = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::ithickness]);
            double dr_tau = EQS::shock_delta(r, GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

            double beta_shock;
            switch (p_pars->method_shock_vel) {
                case isameAsBW:
                    beta_shock = EQS::Beta(Gamma);
                    break;
                case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
            double intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                               p_syna->getPars()->method_tau);
            flux_dens = (intensity * r * r * dr); //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
//            dFnu+=flux_dens;
            /// save the result in image
            ctheta = interpSegLin(ia, ib, t_e, tburst, Dt[BW::Q::ictheta]);
            double theta = interpSegLin(ia, ib, t_e, tburst, Dt[BW::Q::itheta]);

            if (p_pars->m_type == BW_TYPES::iFSRS) {
                if (p_pars->do_rs &&
                    !(Dt[BW::Q::ithichness_rs][ia] == 0 || Dt[BW::Q::ithichness_rs][ib] == 0)) {
                    Interp2d int_em_rs(p_pars->m_freq_arr, r_arr, p_pars->m_synch_em_rs);
                    Interp2d int_abs_rs(p_pars->m_freq_arr, r_arr, p_pars->m_synch_abs_rs);
                    double em_prime_rs = int_em_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                    double abs_prime_rs = int_abs_rs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
                    double em_lab_rs = em_prime_rs / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                    double abs_lab_rs = abs_prime_rs * delta_D; // conversion of absorption (see vanEerten+2010)

                    double GammaShock_rs = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::iGammaRsh]);
                    double dr_rs = interpSegLog(ia, ib, t_e, tburst, Dt[BW::Q::ithichness_rs]);

                    double dr_tau_rs = EQS::shock_delta(r,GammaShock_rs); // TODO this is added becasue in Johanneson Eq. I use ncells

                    double beta_shock_rs;
                    switch (p_pars->method_shock_vel_rs) {

                        case isameAsBW:
                            beta_shock_rs = EQS::Beta(Gamma);
                            break;
                        case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                            beta_shock_rs = EQS::Beta(GammaShock_rs);//us / sqrt(1. + us * us);
                            break;
                    }
                    double ashock_rs = (1.0 - mu * beta_shock_rs); // shock velocity beaming factor
                    dr_rs /= ashock_rs; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
                    dr_tau_rs /= ashock_rs;
                    double dtau_rs = RadiationBase::optical_depth(abs_lab_rs, dr_tau_rs, mu, beta_shock_rs);
                    double intensity_rs = RadiationBase::computeIntensity(em_lab_rs, dtau_rs,
                                                                          p_syna->getPars()->method_tau);
                    if (intensity_rs < 0 || !std::isfinite(intensity_rs)) {
                        (*p_pars->p_log)(LOG_ERR, AT) << "intensity_rs = " << intensity_rs << "\n";
                        exit(1);
                    }
                    double flux_dens_rs = intensity_rs * r * r * dr_rs;
                    if (flux_dens_rs < 0 || !std::isfinite(flux_dens_rs)) {
                        (*p_pars->p_log)(LOG_ERR, AT) << "flux_dens_rs = " << flux_dens_rs << "\n";
                        exit(1);
                    }
                    flux_dens += flux_dens_rs; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                }
            }

        }
        else{
            r = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iR]);
            if (!std::isfinite(r)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R is NAN in integrand for radiation"
                                             << " t_e=" << t_e
                                             << " mD[BW::Q::iR][ia]=" << Dt[BW::Q::iR][ia]
                                             << " mD[BW::Q::iR][ib]=" << Dt[BW::Q::iR][ib]
                                             << " mD[BW::Q::itburst][ia]=" << Dt[BW::Q::itburst][ia]
                                             << " mD[BW::Q::itburst][ib]=" << Dt[BW::Q::itburst][ib]
                    << "\n";
                // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << mD[BW::Q::iR] << "\n";
//            std::cout << " Gamma= " << mD[BW::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
//                return 0.;
                exit(1);
            }

            double rho = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::irho]);
            double Gamma = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iGamma]);
            double GammaSh = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iGammaFsh]);
            double beta = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::ibeta]);
            double U_p = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM2));
            double theta = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::itheta]);
            double rho2 = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::irho2]);
            double m2 = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iM2]);
            double frac = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iacc_frac]);
            double thick = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::ithickness]);
            double gm = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::igm]);
            double gM = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::igM]);
            double gc = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::igc]);
            double B = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iB]);
            double Theta = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iTheta]);
            double z_cool = interpSegLog(ia, ib, t_e, Dt[BW::Q::itburst], Dt[BW::Q::iz_cool]);

            if (rho < 0. || Gamma < 1. || !std::isfinite(Gamma)
                || U_p < 0. || theta <= 0. || rho2 < 0. || thick <= 0.) {
                std::cerr << " wrong value in interpolation to EATS surface  \n"
                          << " r = " << r << "\n"
                          << " rho = " << rho << "\n"
                          << " Gamma = " << Gamma << "\n"
                          << " U_p = " << U_p << "\n"
                          << " theta = " << theta << "\n"
                          << " rho2 = " << rho2 << "\n"
                          << " thick = " << thick << "\n"
                          << " t_e = " << t_e << "\n"
                          << " Exiting...\n";
                exit(1);
            }

            flux_dens = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2,
                                                       frac, B, gm, gM, gc, Theta, z_cool,
                                                       t_e, mu, r, thick, thick, nu_obs, params);
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
                /// evaluateShycnhrotronSpectrum the 'thickness' of the shock (emitting region)
    //            double dr_rs = thick3;

    //          // TODO check if for the shock velocity gamma43 has to be used!!! I think it is gam43! See gammaAdi calc.
                dFnu_rs = Observables::shock_synchrotron_flux_density(
                        Gamma, gamma43, rho3, U_e3, t_e, mu, R, thick3, thick3, params );
            }
            dFnu += dFnu_rs;
#endif
            if (flux_dens == 0 || !std::isfinite(flux_dens)) {
                // REMOVING LOGGER
                std::cerr << " flux density is zero ( dFnu = 0 )" << "\n";
            }

            if (p_pars->m_type == BW_TYPES::iFSRS) {
                (*p_pars->p_log)(LOG_ERR,AT) << " not supported her...";
                exit(1);
            }

//            p_pars->o_gam = Gamma;
//            p_pars->o_r = R;
//            p_pars->o_mu = mu;
//            p_pars->o_flux = dFnu;
//            p_pars->o_theta_j = theta;
        }
    }

    static double shock_synchrotron_flux_density(
            double Gamma, double GammaShock, double m2, double rho2, double acc_frac, double B,
            double gm, double gM, double gc, double Theta, double z_cool,
            double t_e, double mu, double R, double dr, double dr_tau,
            double nu_obs, void * params){

        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & p_syna = p_pars->p_syna;//->getAnSynch();

        /// relativistic effects on the emitting region
        double beta = EQS::Beta(Gamma);
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        double beta_shock;
        switch (p_pars->method_shock_vel) {

            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                break;
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
        dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
        dr_tau /= ashock;

        /// properties of the emitting electrons
        double nuprime = (1.0 + p_pars->z ) * nu_obs * delta_D;
        double Ne = m2 / CGS::mp; // numer of protons/electrons
        double nprime = rho2 / CGS::mp; // number density of protons/electrons
        double em_prime,em_lab,abs_prime,abs_lab,intensity,flux_dens,dtau;


#if 1
        switch (p_pars->m_method_ne) {
            // default (presumably more correct version)
            case iusenprime:
                /// this is informed by van Earten et al. and afterglowpy
                p_syna->evaluateShycnhrotronSpectrum(nprime, Ne, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime);
                em_prime = p_syna->getPars()->em;
                em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
                abs_prime = p_syna->getPars()->abs;
                abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_syna->getPars()->method_tau);
//                flux_dens = (intensity * R * R * dr) * (1.0 + p_eats->z) / (2.0 * p_eats->d_l * p_eats->d_l);

                flux_dens = (intensity * R * R * dr) ;
                break;
            case iuseNe: // TODO THIS SHOULD NOT BE USED (tau is incorrectly estimated)
                /// This is informed by the G. Lamb and Fernandez et al.
                p_syna->evaluateShycnhrotronSpectrum(Ne, Ne, acc_frac, B, gm, gM, gc, Theta, z_cool, nuprime);
                em_prime = p_syna->getPars()->em;
                em_lab = em_prime / (delta_D * delta_D);
                em_lab /= delta_D; // TODO this should be from 'dr'...
                abs_lab = abs_prime * delta_D; // TODO with Ne this might not work as we do not use 'dr' of the shock...
                dtau = RadiationBase::optical_depth(abs_lab,dr_tau,mu,beta_shock);
                intensity = RadiationBase::computeIntensity(em_lab, dtau,
                                                            p_syna->getPars()->method_tau);
//                flux_dens = intensity * (1.0 + p_eats->z) / (p_eats->d_l * p_eats->d_l) / 10;
                flux_dens = intensity / 5.; // TODO why no '2'?? // why this need /10 to fit the upper result?
                break;
        }
#endif
        return flux_dens;
    }

    bool evalEATSindexes(size_t &ia, size_t &ib, double t_obs, double z, //size_t m_i_end_r,
                         double ctheta, double cphi, double theta_obs,
                         double (*obs_angle_func)( const double &, const double &, const double & )){

//        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto & m_mu =  p_pars->m_mu;
        auto & ttobs =  p_pars->ttobs;
        size_t m_i_end_r= p_pars->i_end_r;
        if (m_mu.size() < 1) {
            m_mu.resize(p_pars->i_end_r);
            ttobs.resize(p_pars->i_end_r);
        }
        /// evaluate mu array
        for (size_t i = 0; i < m_i_end_r; i++)
            m_mu[i] = (mD[BW::Q::itburst][i] - t_obs / (1.0 + z) ) / mD[BW::Q::iR][i] * CGS::c;
        double mu = obs_angle_func(ctheta, cphi, theta_obs);
        for (size_t i_ = 0; i_ < m_i_end_r; i_++) {
            ttobs[i_] = mD[BW::Q::itt][i_] + mD[BW::Q::iR][i_] / CGS::c * (1.0 - mu);
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

        double rej = interpSegLog(ia, ib, t_obs, p_pars->ttobs, mD[BW::Q::iR]);
        double rho_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, mD[BW::Q::iEJrho]);
        double delta_ej_cell = interpSegLog(ia, ib, t_obs, p_pars->ttobs, mD[BW::Q::iEJdelta]);

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
