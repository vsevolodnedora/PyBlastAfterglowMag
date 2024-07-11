//
// Created by vsevolod on 05/04/23.
//

#ifndef SRC_BLASTWAVE_H
#define SRC_BLASTWAVE_H

//#include "../utilitites/pch.h"
//#include "../utilitites/utils.h"
//#include "../utilitites/interpolators.h"
//#include "../utilitites/ode_solvers.h"
//#include "../utilitites/quadratures.h"
//#include "../utilitites/rootfinders.h"
//#include "../image.h"
//#include "synchrotron_an.h"
//#include "../composition.h"
//#include "blastwave_components.h"
//#include "blastwave_pars.h"
//#include "blastwave_base.h"
//#include "blastwave_radiation.h"
#include "blastwave_radiation.h"
//#include "eats.h"

/// Main Blastwave class
class BlastWave : public BlastWaveRadiation{
    std::unique_ptr<logger> p_log;
public:
    enum CASES { i_INSIDE_BEHIND, i_OUTSIDE_BEHIND, i_INSIDE_ABOVE, i_OUTSIDE_ABOVE, i_AT_ZERO_INSIDE };
    BlastWave(Vector & tb_arr, size_t ishell, size_t ilayer, size_t n_substeps,
              BW_TYPES type, CommonTables & commonTables, int loglevel )
        : BlastWaveRadiation(tb_arr, ishell, ilayer, n_substeps, type, commonTables, loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW");
    }
    // ------------------------------------------------------
    ~BlastWave() = default;
    // --------------------------------------------------------
    void updateNucAtomic( const double * sol, const double t ){
        p_nuc->update(
                p_pars->Ye0,
                p_pars->M0 * (double)EjectaID2::CellsInLayer(p_pars->ilayer), // m_iso ( for eps_th )
//                EQS::BetFromMom( sol[ p_pars->ii_eq + SOL::QS::imom ] ),
                EQS::BetaFromGamma( sol[ p_pars->ii_eq + SOL::QS::iGamma ] ),
                t,
                sol[ p_pars->ii_eq + SOL::QS::iR ],
                p_pars->s0
        );
        p_pars->kappa = p_nuc->getPars()->kappa;
        p_pars->dEnuc = p_pars->M0 * p_nuc->getPars()->eps_nuc_thermalized;
    }
    void updateCurrentBpwn( const double * sol ){
        double r_w = sol[p_pars->ii_eq + SOL::iRw];
        if (r_w < 1 || !std::isfinite(r_w)){
            (*p_log)(LOG_ERR,AT) << " r_w is not set or nan \n";
            exit(1);
        }
        double u_b_pwn = 3.0*sol[p_pars->ii_eq + SOL::iWepwn]/4.0/M_PI/r_w/r_w/r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
        double b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
        p_pars->curr_b_pwn = b_pwn;
    }
    void updateEnergyInjection( double ldip, double lacc ){
        if (ldip < 0 or !std::isfinite(ldip)){
            (*p_log)(LOG_ERR,AT) << " ldip < 0 or nan; ldip="<<ldip<<"\n";
        }
        if (lacc < 0 or !std::isfinite(lacc)){
            (*p_log)(LOG_ERR,AT) << " lacc < 0 or nan; lacc="<<lacc<<"\n";
        }
        // ----------------------
        p_pars->curr_ldip = ldip;
        p_pars->curr_lacc = lacc;
    }

    /// set initial condition 'ic_arr' for this blast wave using data from PWNPars{} struct
    void setInitConditions( double * ic_arr, size_t i ) {
        if (p_pars->end_evolution)
            return;
        /// if layer does not have energy / mass -- do not evolve it
        if ((p_pars->M0 == 0.) && (p_pars->E0 == 0.)){
//            std::cerr << AT << "\n "
//                      << "[ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
//                      <<" M0=0 and E0=0 -> Ignoring this layer.\n";
            p_pars->end_evolution = true;
            for (size_t v = 0; v < SOL::neq; ++v){
                ic_arr[i+v] = 0.;
            }
            return;
        }
//        p_pars->Gamma0 = EQS::GamFromMom(p_pars->mom0);
        p_pars->beta0 = EQS::BetaFromGamma(p_pars->Gamma0);
        // ****************************************
//        p_pars->R0    = p_pars->tb0 * beta0 * CGS::c;

        p_dens->evaluateRhoDrhoDrDefault(p_pars->R0, p_pars->ctheta0);
        if (p_pars->is_within0){ //p_pars->j_i0!=123456789
            // behind a jet BW
            p_dens->m_rho_ = p_dens->m_rho_floor_val * p_dens->m_rho_def;
            p_dens->m_drhodr_ = 0.;
        }
        else{
            // outside jet BW
            p_dens->m_rho_ = p_dens->m_rho_def;
            p_dens->m_drhodr_=p_dens->m_drhodr_def;
        }

        double m_M20 = (2.0 / 3.0) * CGS::pi
                       * (std::cos(p_pars->theta_a) - std::cos(p_pars->theta_b0))
                       * p_dens->m_rho_
                       * std::pow(p_pars->R0, 3)
                       / p_pars->ncells; // mass accreted by the shock by the time tburst = tburst[0]
        double adi0 = p_eos->getGammaAdi(p_pars->Gamma0,p_pars->beta0);
        double GammaSh0 = EQS::GammaSh(p_pars->Gamma0,adi0);
        // ****************************************
        if ((p_pars->mom0 <= 0.) || (!std::isfinite(p_pars->mom0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " Mom0 < 0 (Mom0=" <<p_pars->Gamma0 << ") "
                                   << "Mom0="<<p_pars->mom0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")"
                                   << " \n";
            //std::cout << "[ Error ] " << "Gamma0 < 0 (Gamma0=" <<Gamma0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->M0 <= 0.) || (!std::isfinite(p_pars->M0))){
//            std::cout << "[ WARNING ] " << "M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            // REMOVING LOGGER
            (*p_log)(LOG_WARN, AT)  << " M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            p_pars->M0 = p_pars->E0 / (p_pars->Gamma0 * CGS::c * CGS::c);
        }
        if ((p_pars->R0 <= 1.) || (!std::isfinite(p_pars->R0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " R0 <= 0 (R0=" <<p_pars->R0 << ") " << "G0="<<p_pars->Gamma0
                                   << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<") "
                                   << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "R0 <= 0 (R0=" <<R0 << ")\n";
            exit(1);
        }
        if ((p_dens->m_rho_) <= 0.|| (!std::isfinite(p_dens->m_rho_))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " rho0 < 0 (rho0=" <<p_dens->m_rho_ << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                                   << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")"
                                   << " \n";
//            std::cerr << AT  << "\n";

            //std::cout << "[ Error ] " << "rho0 < 0 (rho0=" <<rho0 << ")\n";
            exit(1);
        }
        if ((p_pars->E0 <= 0.) || (!std::isfinite(p_pars->E0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << "E0 <= 0 (E0=" <<p_pars->E0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                                   << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")"
                                   << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "E0 < 0 (E0=" <<E0 << ")\n";
            exit(1);
        }
        if ((p_pars->Gamma0 < 1.) || (!std::isfinite(p_pars->Gamma0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " Gamma0 < 1 (Gamma0=" <<p_pars->Gamma0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")"
                                   << " \n";
            //std::cout << "[ Error ] " << "Gamma0 < 0 (Gamma0=" <<Gamma0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->theta_b0 < 0.) || (!std::isfinite(p_pars->theta_b0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " theta_b0 < 0 (theta_b0=" <<p_pars->theta_b0 << ") "
                                   << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")"
                                   << " \n";
            //std::cout << "[ Error ] " << "theta_b0 < 0 (theta_b0=" <<theta_b0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->ncells < 1 )|| (!std::isfinite(p_pars->ncells))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " ncells < 1 (ncells=" <<p_pars->ncells << ")"
                                   << " \n";
            //std::cout << "[ Error ] " << "ncells < 1 (ncells=" <<ncells << ")\n";
            exit(1);
        }
        if (cos(p_pars->theta_a) <= cos(p_pars->theta_b0)){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " cos(theta_a) < cos(theta_b0) (theta_a="<<p_pars->theta_a
                                   << ", theta_b0="<<p_pars->theta_b0<<")"
                                   << " \n";
            //std::cout << "[ Error ] " <<" cos(theta_a) < cos(theta_b0) (theta_a="<<theta_a<<", theta_b0="<<theta_b0<<")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        bool use_spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        p_pars->Rd = std::pow(3./(4.*CGS::pi)
                      * 1./(CGS::c*CGS::c*CGS::mp) *
                              p_pars->E0/( p_dens->m_rho_ / CGS::mp * p_pars->Gamma0*p_pars->Gamma0), 1./3.);

        double x = p_pars->Rd / p_pars->R0;
        // ***************************************
        // -------------- DYNAMICS ---------------
        ic_arr[i + SOL::QS::iR]      = p_pars->R0;//m_tb_arr[0] * beta0 * CGS::c; CO
//        ic_arr[i + SOL::QS::iRsh]    = p_pars->R0;//m_tb_arr[0] * EQS::Beta(GammaSh0) * CGS::c; # TODO change to Rsh
        ic_arr[i + SOL::QS::itt]     = EQS::init_elapsed_time(p_pars->R0, p_pars->mom0, use_spread);
        ic_arr[i + SOL::QS::itcomov] = EQS::initTComov(p_pars->R0, p_pars->beta0, p_pars->Gamma0);
        ic_arr[i + SOL::QS::iEint2]  = (p_pars->Gamma0 - 1. ) * m_M20 / p_pars->M0;  //TODO Isnt it just E0 / m_M0 ???? As M0 = E0 * cgs.c ** -2 / Gamma0
        ic_arr[i + SOL::QS::iEint2]  += p_pars->Eint0 / p_pars->M0 / CGS::c / CGS::c; // add initial internal energy
//        ic_arr[i + SOL::QS::imom]    = MomFromGam(p_pars->Gamma0);//std::log( p_pars->Gamma0 );
        ic_arr[i + SOL::QS::iGamma]  = p_pars->Gamma0;//std::log( p_pars->Gamma0 );
        ic_arr[i + SOL::QS::itheta]  = p_pars->theta_b0;
        ic_arr[i + SOL::QS::iErad2]  = 0.0;
        ic_arr[i + SOL::QS::iEsh2]   = 0.0;
        ic_arr[i + SOL::QS::iEad2]   = 0.0;
        ic_arr[i + SOL::QS::iM2]     = m_M20 / p_pars->M0;
        // ------------ PWN -------------------
        if ((p_pars->curr_ldip < 0 || p_pars->curr_lacc < 1) && p_pars->m_type==iFS_PWN_DENSE){
            (*p_log)(LOG_ERR, AT)  << "p_pars->curr_ldip = " << p_pars->curr_ldip
                                   << " p_pars->curr_lacc = "<< p_pars->curr_lacc
                                   << "\n";
            exit(1);
        }
        ic_arr[i + SOL::QS::iRw]   = .5 * ic_arr[i + SOL::QS::iR]; // TODO this is temporary!!!
        ic_arr[i + SOL::QS::iWenb] = p_pars->eps_e_w * p_pars->curr_ldip
                                   + p_pars->epsth_w * p_pars->curr_lacc;
        ic_arr[i + SOL::QS::iWenb] = p_pars->eps_mag_w * p_pars->curr_ldip;
        ic_arr[i + SOL::QS::iWtt]  = ic_arr[i + SOL::QS::itt]; // we assume that also time is the same
        /// reverse shock
        if (p_pars->m_type==BW_TYPES::iFSRS){
            ic_arr[i + SOL::QS::iEint3] = 0.; // assume shell is initially cold #1e-3*ic_arr[i + SOL::QS::iEint2];
            if (p_pars->init_deltaR4)
                ic_arr[i + SOL::QS::ideltaR4] = CGS::c * p_pars->tprompt * p_pars->beta0;
        }
        // ***************************************
        for (size_t v = 0; v < SOL::neq; ++v){
            if (!std::isfinite(ic_arr[i + v])){
                (*p_log)(LOG_ERR, AT)  << " NAN in initial data for shell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                                       << " v_n="<<SOL::vars[v]<<" val="<<ic_arr[i + v]<<" Exiting...\n";
//                std::cerr << AT  << "\n";
                exit(1);
            }
        }
        p_spread->m_theta_b0 = p_pars->theta_b0;
        // ***************************************
    }
    /// check the evolution result
    void checkEvolutionEnd(){
        /// limit the evaluation to the latest 'R' that is not 0 (before termination)
        if (p_pars->end_evolution)
            return;
        size_t nr = m_data[BW::Q::iR].size();
        size_t i_end_r = nr;
        for(size_t ir = 0; ir < nr; ++ir){
            if (m_data[BW::Q::iR][ir] == 0.) {
                i_end_r = ir;
                break;
            }
        }
//        if (i_end_r == 0){
//            (*p_log)(LOG_ERR,AT) << "[il="<< p_pars->ilayer << ", ish="<< p_pars->ishell
//                                            << "] Blastwave was not evolved: i_end_r = " << i_end_r << "\n";
//            exit(1);
//        }
        if (i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "[ish="<<p_pars->ishell<<" il="<<p_pars->ilayer<<"] "
                <<"beta0="<<p_pars->beta0<< " i_end_r = 0"<<"\n";
            exit(1);
        }
        p_pars->i_end_r = i_end_r;
    }
    /// add the current solution 'sol' to the 'm_data' which is Vector of Arrays (for all variables)
    void insertSolution( const double * sol, double t, size_t it, size_t i, VecVector & Dat ) {

        Dat[BW::Q::itburst][it]   = t;
        Dat[BW::Q::iR][it]        = sol[i + SOL::QS::iR]; // TODO you do not need 'i' -- this is p_pars->ii_eq
        Dat[BW::Q::iRsh][it]      = sol[i + SOL::QS::iRsh];
        Dat[BW::Q::itt][it]       = sol[i + SOL::QS::itt];
        Dat[BW::Q::imom][it]      = EQS::MomFromGamma(sol[i + SOL::QS::iGamma]);
        Dat[BW::Q::iGamma][it]    = sol[i + SOL::QS::iGamma];
        Dat[BW::Q::ibeta][it]     = EQS::BetaFromGamma(sol[i + SOL::QS::iGamma]);//sol[i+QS::iGamma];
        Dat[BW::Q::itheta][it]    = sol[i + SOL::QS::itheta];
        Dat[BW::Q::iM2][it]       = sol[i + SOL::QS::iM2];
        Dat[BW::Q::itcomov][it]   = sol[i + SOL::QS::itcomov];
        Dat[BW::Q::iEad2][it]     = sol[i + SOL::QS::iEad2];
        Dat[BW::Q::iEint2][it]    = sol[i + SOL::QS::iEint2];
        Dat[BW::Q::iEsh2][it]     = sol[i + SOL::QS::iEsh2];
        Dat[BW::Q::iErad2][it]    = sol[i + SOL::QS::iErad2];

        /// reverse shock
        if (p_pars->m_type==BW_TYPES::iFSRS){
//            double rrsh = sol[i + SOL::QS::iR] - sol[i + SOL::QS::iRrsh];
            Dat[BW::Q::iRrsh][it]     = sol[i + SOL::QS::iRrsh];
            Dat[BW::Q::iEint3][it]    = sol[i + SOL::QS::iEint3];
            Dat[BW::Q::iEad3][it]     = sol[i + SOL::QS::iEad3];
            Dat[BW::Q::iErad3][it]    = sol[i + SOL::QS::iErad3];
            Dat[BW::Q::iEsh3][it]     = sol[i + SOL::QS::iEsh3];
            Dat[BW::Q::iM3][it]       = sol[i + SOL::QS::iM3];
            Dat[BW::Q::ideltaR4][it]  = sol[i + SOL::QS::ideltaR4];
        }

        /// density upstream
        if (p_pars->m_type==BW_TYPES::iFS_DENSE){
            /// pass
        }

        /// PWN
        if (p_pars->m_type==BW_TYPES::iFS_PWN_DENSE){
            Dat[BW::Q::i_Wtt][it]     = sol[i + SOL::QS::iWtt];
            Dat[BW::Q::i_Wmom][it]    = sol[i + SOL::QS::iWmom];
            Dat[BW::Q::i_Wepwn][it]   = sol[i + SOL::QS::iWepwn];
            Dat[BW::Q::i_Wenb][it]    = sol[i + SOL::QS::iWenb];
            Dat[BW::Q::i_WGamma][it]  = EQS::GammaFromMom(sol[i + SOL::QS::iWmom]);
        }

        if (sol[i+SOL::QS::iR] < 1. || Dat[BW::Q::iGamma][it] < 1. || sol[i + SOL::QS::iGamma] < 1) {
            (*p_log)(LOG_ERR, AT)  << "Wrong value at i=" << it << " tb=" << sol[i + SOL::QS::iR]
                                   << " iR="      << sol[i + SOL::QS::iR]
//                                   << " iRsh="    << sol[i + SOL::QS::iRsh]
                                   << " iGamma="  << sol[i + SOL::QS::iGamma]
                                   << " itheta="  << sol[i + SOL::QS::itheta]
                                   << " iM2="     << sol[i + SOL::QS::iM2]
                                   << " itcomov=" << sol[i + SOL::QS::itcomov]
                                   << " iEad2="   << sol[i + SOL::QS::iEad2]
                                   << " iEint2="  << sol[i + SOL::QS::iEint2]
                                   << " iEsh2="   << sol[i + SOL::QS::iEsh2]
                                   << " iErad2="  << sol[i + SOL::QS::iErad2]
                                   << " "
                                   << " Exiting...\n";
            std::cerr << AT  << "\n";
            exit(1);
        }
    }
    /// add the current solution 'sol' to the 'm_data' which is Vector of Arrays (for all variables)
    void insertSolution( const double * sol, size_t it, size_t i ) {
        if (p_pars->end_evolution)
            return;
        double t = m_tb_arr[it];
        insertSolution(sol, t, it, i, m_data);
    }
    /// Mass and energy are evolved in units of M0 and M0c2 respectively
    void applyUnits( double * sol, size_t i ) {
        /// quantities for forward shock
        sol[i + SOL::QS::iM2]    *= p_pars->M0;
        sol[i + SOL::QS::iEint2] *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEad2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEsh2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iErad2] *= (p_pars->M0 * CGS::c * CGS::c);
        /// quantities for reverse shock
        sol[i + SOL::QS::iM3]    *= p_pars->M0;
        sol[i + SOL::QS::iEint3] *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEad3]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iEsh3]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + SOL::QS::iErad3] *= (p_pars->M0 * CGS::c * CGS::c);
    }
    /// insert the solution at every substep used in derivatives inside RHS
    void insertSolutionSustep(const double * sol, double t, size_t it, size_t i){
        double tmp[p_pars->n_substeps];
        /// Move previous solutions down in the tmp contrainer
        for (auto & j : BW::VARS.at(p_pars->m_type)) {
            for (size_t i = 0; i < p_pars->n_substeps-1; i++) {
                tmp[i] = mDtmp[j][i];
            }
            for (size_t i = 0; i < p_pars->n_substeps-1; i++) {
                mDtmp[j][i + 1] =tmp[i];
            }
            mDtmp[j][0] = 0.; // to be filled by the current solution
        }
        if (!p_pars->end_evolution)
        /// Store the current solution in the 'tmp' container
        insertSolution(sol, t, 0, i, mDtmp);
        int x = 1;

    }
    /// check if to terminate the evolution
    bool isToTerminate( double * sol, size_t i, size_t ix ) {
//        double mom = sol[i + SOL::QS::imom];
        double igamma = sol[i+SOL::QS::iGamma];//EQS::GamFromMom(mom);
//        double igamma = EQS::GamFromMom(mom);
        double ibeta = EQS::BetaFromGamma(igamma);//EQS::BetFromMom(mom);
        double iEint2 = sol[i + SOL::QS::iEint2];
        bool do_terminate = false;
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
//            double betaSh = EQS::Beta(igamma);
            if ((iEint2 <= 0.)||(!std::isfinite(iEint2))||(igamma < 1.)||(ibeta < 0)){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution at ix="<<ix
                        << " [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: Eint2 < 0 or NAN ("<<iEint2<<")\n";
            }
            if ((std::abs(ibeta) < p_pars->min_beta_terminate)||(!std::isfinite(ibeta))){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution at ix="<<ix
                        << " [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: betaSh < beta_min or NAN ("<<ibeta<<") with beta_min= "
                        << p_pars->min_beta_terminate<<"\n";
            }
            if (do_terminate) {
                (*p_log)(LOG_WARN, AT) << " TERMINATION Previous iterations: "
                                       << " E0=" << string_format("%.9e", p_pars->E0)
                                       << " Gamma0=" << string_format("%.9e", p_pars->Gamma0)
                                       << " M0=" << string_format("%.9e", p_pars->M0)
                                       << " \n";
                (*p_log)(LOG_WARN, AT) << " i=Current" << " Eint2=" << string_format("%.9e", sol[i + SOL::QS::iEint2])
                                       << " Gamma=" << string_format("%.9e", igamma)
                                       << " betaSh=" << string_format("%.9f", ibeta)
                                       << " M2=" << string_format("%.9e", sol[i + SOL::QS::iM2])
                                       << " rho=" << string_format("%.9e", p_dens->m_rho_)
                                       << " \n";
                for (size_t j = 0; j < p_pars->n_substeps; j++){
                    (*p_log)(LOG_WARN, AT) << " ix=" << ix - j << " Eint2="
                                           << string_format("%.9e", mDtmp[BW::Q::iEint2][j])
                                           << " Gamma=" << string_format("%.9e", mDtmp[BW::Q::iGamma][j])
                                           << " betaSh=" << string_format("%.9f", mDtmp[BW::Q::ibeta][j])
                                           << " M2=" << string_format("%.9e", mDtmp[BW::Q::iM2][j])
                                           << " rho=" << string_format("%.9e", mDtmp[BW::Q::irho][j])
                                           << " \n";
                }
            }
        }
        if (do_terminate && (!p_pars->allow_termination)){
            (*p_log)(LOG_ERR,AT)<<" termination is not allowed. Exiting.\n";
            exit(1);
        }
        return do_terminate;
    }
    /// check if to terminate lateral spreading
    bool isToStopLateralExpansion( double * sol, size_t i ) {
        double itheta = sol[i + SOL::QS::itheta];
        double iEint2 = sol[i + SOL::QS::iEint2];
        bool do_terminate = false;
        if (!p_pars->end_spreading) {
            if (itheta > 0.99 * p_pars->theta_max){
//                std::cerr << AT << " \n"
//                          << " Stopping Jet Lateral Expansion  ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer
//                          << " with Gamma0=" << p_pars->Gamma0 << " [iteration=" << p_pars->comp_ix << "] REASON :: theta > 0.99 theta_max ("<<string_format("%.3e", p_pars->theta_max)<<") \n";
//                std::cerr << " Last Gammas: " << "\n"
//                          << " Gamma[" << p_pars->comp_ix << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix)) << "\n"
//                          << " theta[" << p_pars->comp_ix << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix))
//                          << " Gamma[" << p_pars->comp_ix - 1 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 1)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 1 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 1))
//                          << " Gamma[" << p_pars->comp_ix - 2 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 2)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 2 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 2))
//                          << " Gamma[" << p_pars->comp_ix - 3 << "]=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 3)) << "\n"
//                          << " theta[" << p_pars->comp_ix - 3 << "]=" << string_format("%.4f", getVal(Q::itheta, p_pars->comp_ix - 3))
//                          << "\n";
                do_terminate = true;
            }
        }
        return do_terminate;
    }
    /// check if RS needs to be terminated
    bool isToStopReverseShock( double * sol, size_t i ) {
        if (p_pars->shutOff || (!p_pars->do_rs))
            return false;
        double iGamma = sol[i + SOL::QS::iGamma];
        double M3 = sol[i + SOL::QS::iM3];
        double rho4 = EQS::rho4(sol[i + SOL::QS::iR],
                                sol[i + SOL::QS::ideltaR4],
                                p_pars->tprompt,
                                p_pars->beta0,
                                p_pars->M0,
                                p_pars->theta_a,
                                p_pars->theta_b0,
                                p_pars->exponential_rho4);
        if ((iGamma < p_pars->Gamma0) &&
            (rho4 < p_pars->rs_shutOff_criterion_rho || M3 >= 1.)){
            p_pars->shutOff = true;
        }
        return p_pars->shutOff;
    }
    /// check if the solution 'makes sense'
    bool isSolutionOk( double * sol, double x, size_t ix, size_t i ) {

        bool no_issues = true;
//        double mom = sol[i + SOL::QS::imom];
//        double mom0 = p_pars->Gamma0 * EQS::Beta(p_pars->Gamma0);
        double igamma = sol[i+SOL::QS::iGamma];
//        double igamma = EQS::GamFromMom(mom);
        double ibeta = EQS::BetaFromGamma(igamma);//EQS::BetFromMom(mom);
//        double ibeta = EQS::BetFromMom(mom);
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
            double beta   = EQS::BetaFromGamma(igamma);
            if ((igamma < 1.) || (beta < 0)) {
                (*p_log)(LOG_ERR,AT)  << AT << " \n"
                                      << " Gamma < 1. ishell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                                      << " with Gamma0="<<p_pars->Gamma0 << " [iteration="<<p_pars->prev_idx_x<<"] \n";
                (*p_log)(LOG_ERR,AT) << " Last solutions: " << "\n";
                for (size_t j = 0; j < p_pars->n_substeps; j++){
                    (*p_log)(LOG_WARN, AT) << " ix=" << ix - j << " Eint2="
                                           << string_format("%.9e", mDtmp[BW::Q::iEint2][j])
                                           << " Gamma=" << string_format("%.9e", mDtmp[BW::Q::iGamma][j])
                                           << " betaSh=" << string_format("%.9f", mDtmp[BW::Q::ibeta][j])
                                           << " M2=" << string_format("%.9e", mDtmp[BW::Q::iM2][j])
                                           << " rho=" << string_format("%.9e", mDtmp[BW::Q::irho][j])
                                           << " \n";
                }

//                std::cerr << "Gamma:" << getData(Q::iGamma) << "\n";
//                std::cerr << "R:" << getData(Q::iR) << "\n";
                (*p_log)(LOG_ERR,AT)  << " Gamma cannot be less than 1. or too large. "
                                         "Found: Gamma(" << igamma << "). Gamma0="<<p_pars->Gamma0<<" \n";
                no_issues = false;
            }
            if (sol[i + SOL::QS::itheta] < p_pars->theta_b0 * 0.99) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical decrease in theta(" << sol[i + SOL::QS::itheta]
                                     << ") < theta_b0(" << p_pars->theta_b0 << ")\n";
                no_issues = false;
            }
            if (sol[i + SOL::QS::itt] <= 0.) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical value of observer time (on axis) itt("
                                     << sol[i + SOL::QS::itt] << ") < 0 \n";
                no_issues = false;
            }
        }
        p_pars->prev_x = x;
        p_pars->prev_idx_x = ix;
        return no_issues;
    }
    /// compute other quantities for BW and add them to the contaienr
    void addOtherVars(size_t it){ addOtherVars(it, m_data); }
    void addOtherVars(size_t it, VecVector & Dat){

        if (p_pars->end_evolution)
            return;

        Dat[BW::Q::ictheta][it] = ctheta(Dat[BW::Q::itheta][it]);//p_pars->ctheta0 + 0.5 * (2. * m_data[Q::itheta][it] - 2. * p_pars->theta_w);

        double rho_prev = Dat[BW::Q::irho][it - 1];
        double rho = Dat[BW::Q::irho][it];
        if ((rho < 0)||(!std::isfinite(rho))){
            (*p_log)(LOG_ERR,AT)<<" negative density!\n";
            exit(1);
        }

        /// related to the jet BW density profile
        Dat[BW::Q::irho][it] = p_dens->m_rho_;
        Dat[BW::Q::idrhodr][it] = p_dens->m_drhodr_;

        /// add density properties
        if (p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            Dat[BW::Q::iGammaCBM][it] = p_dens->m_GammaRho;
            Dat[BW::Q::iGammaREL][it] = p_dens->m_GammaRel;
            Dat[BW::Q::idGammaCBMdr][it] = p_dens->m_dGammaRhodR;
            Dat[BW::Q::idGammaRELdGamma][it] = p_dens->m_dGammaReldGamma;
            Dat[BW::Q::idPCBMdrho][it] = p_dens->m_dPCBMdrho;
            Dat[BW::Q::iPcbm][it] = p_dens->m_P_cbm;
            Dat[BW::Q::iMCBM][it] = p_dens->m_M_cbm;
            Dat[BW::Q::iCSCBM][it] = p_dens->m_CS_CBM;
            if (Dat[BW::Q::iGammaREL][it] < 1.) {
                Dat[BW::Q::iGammaREL][it] = Dat[BW::Q::iGamma][it];
                //            std::cerr << AT << "\n GammaRel="<<m_data[Q::iGammaREL][it]<<"; Exiting...\n";
                //            exit(1);
            }
        }
        // other parameters

//        if (p_pars->end_evolution)
//            return;

        if ((it>1)&&(rho_prev < 1e-10 * m_data[BW::Q::irho][it])){
            (*p_log)(LOG_ERR,AT) << " it="<<it<<" density gradient >10 orders of magnitude\n";
//            exit(1);
        }

        Dat[BW::Q::iadi][it] = p_eos->getGammaAdi(
                Dat[BW::Q::iGamma][it], // TODO ! is it adi or adi21 (using GammaRel)??
                Dat[BW::Q::ibeta][it]);
        Dat[BW::Q::irho2][it] = EQS::rho2t(
                Dat[BW::Q::iGamma][it], // TODO should there be a gammaRel?? with adi43??..
                Dat[BW::Q::iadi][it],
                Dat[BW::Q::irho][it]);

        /// FS: shock front velocity
        Dat[BW::Q::iGammaFsh][it] = EQS::GammaSh(Dat[BW::Q::iGamma][it], Dat[BW::Q::iadi][it]);
//        switch (p_pars->p_mphys->m_method_gamma_fsh) {
//
//            case iuseGammaShock:
//                Dat[BW::Q::iGammaFsh][it] = EQS::GammaSh(
//                        Dat[BW::Q::iGamma][it], Dat[BW::Q::iadi][it]);
//                break;
//            case iuseJustGamma:
//                Dat[BW::Q::iGammaFsh][it] = Dat[BW::Q::iGamma][it];
//                break;
//            case iuseJustGammaRel:
//                Dat[BW::Q::iGammaFsh][it] = Dat[BW::Q::iGammaREL][it];
//                break;
//            case iuseGammaRelShock:
//                Dat[BW::Q::iGammaFsh][it] = EQS::GammaSh(
//                        Dat[BW::Q::iGammaREL][it], Dat[BW::Q::iadi][it]);
//                break;
//        }

        /// FS: shock front radius
//        switch (p_pars->p_mphys->m_method_r_sh) {
//            case isameAsR:
//                Dat[BW::Q::iRsh][it] = Dat[BW::Q::iR][it];
//                break;
//            case iuseGammaSh:
//                break;
//        }

        /// FS: shock thickness
        Dat[BW::Q::ithickness][it] = p_pars->p_mphys->get_shock_thickness(
                Dat[BW::Q::iR][it],Dat[BW::Q::iM2][it],Dat[BW::Q::itheta][it],Dat[BW::Q::iGamma][it],
                Dat[BW::Q::irho2][it], p_pars->ncells);
//        if (Dat[BW::Q::ithickness][it] == 0)
//            std::cout << AT << " " << Dat[BW::Q::ithickness][it] << "\n";

        /// FS: shock downstream energy density
        Dat[BW::Q::iU_p][it] = p_pars->p_mphys->get_shock_Up(
                Dat[BW::Q::iGammaFsh][it],Dat[BW::Q::irho2][it],Dat[BW::Q::iM2][it],
                Dat[BW::Q::iEint2][it]);


        /// For reverse shock
        if ((it>0) && (p_pars->do_rs) && (!p_pars->shutOff) && (p_pars->m_type==BW_TYPES::iFSRS)){
            Dat[BW::Q::iGamma43][it] = (double)EQS::get_Gamma43(
                    Dat[BW::Q::iGamma][it],
                    p_pars->Gamma0,
                    EQS::BetaFromMom(Dat[BW::Q::imom][it]),
                    p_pars->beta0);

            Dat[BW::Q::irho4][it] = EQS::rho4(Dat[BW::Q::iR][it],
                                                 Dat[BW::Q::ideltaR4][it],
                                                 p_pars->tprompt,
                                                 p_pars->beta0,
                                                 p_pars->M0,
                                                 p_pars->theta_a,
                                                 p_pars->theta_b0,
                                                 p_pars->exponential_rho4);

            Dat[BW::Q::iadi3][it] = p_eos->getGammaAdi(
                    Dat[BW::Q::iGamma43][it],
                    EQS::BetaFromGamma(Dat[BW::Q::iGamma43][it]));

            Dat[BW::Q::irho3][it] = EQS::rho2t(
                    Dat[BW::Q::iGamma][it],
                    Dat[BW::Q::iadi3][it],
                    Dat[BW::Q::irho4][it] // /Dat[BW::Q::iGamma][0] // rho4 -> rho4prime
                    ); // TODO Check if here Gamma and adi3 are used
//            double rho3prim = 4 * Dat[BW::Q::irho4][it] * Dat[BW::Q::iGamma][it];
//            std::cout << " " << Dat[BW::Q::irho3][it] << " " << rho3prim << "\n";

            /// shock front velocity
            Dat[BW::Q::iGammaRsh][it] = EQS::GammaSh(Dat[BW::Q::iGamma43][it], Dat[BW::Q::iadi3][it]);
//            switch (p_pars->p_mphys_rs->m_method_gamma_rsh) {
//
//                case iuseGammaShock:
//                    Dat[BW::Q::iGammaRsh][it] = EQS::GammaSh(
//                            Dat[BW::Q::iGamma43][it], Dat[BW::Q::iadi3][it]
//                            );
//                    break;
//                case iuseJustGamma:
//                    Dat[BW::Q::iGammaRsh][it] = Dat[BW::Q::iGamma43][it];
//                    break;
//                case iuseJustGammaRel:
//                    (*p_log)(LOG_ERR,AT)<<" not implemented method_gamma_rsh = useJustGammaRel \n";
//                    exit(1);
//                    break;
//                case iuseGammaRelShock:
//                    (*p_log)(LOG_ERR,AT)<<" not implemented method_gamma_rsh = useGammaRelShock \n";
//                    exit(1);
//                    break;
//            }

            /// FS: shock front radius
//            switch (p_pars->p_mphys_rs->m_method_r_sh) {
//                case isameAsR:
//                    Dat[BW::Q::iRrsh][it] = Dat[BW::Q::iR][it];
//                    break;
//                case iuseGammaSh:
//                    break;
//            }

            Dat[BW::Q::ithickness_rs][it] = p_pars->p_mphys_rs->get_shock_thickness(
                    Dat[BW::Q::iR][it],
                    Dat[BW::Q::iM3][it],
                    Dat[BW::Q::itheta][it],
                    Dat[BW::Q::iGamma43][it],
                    Dat[BW::Q::irho3][it], p_pars->ncells);

            Dat[BW::Q::iU_p3][it] = p_pars->p_mphys_rs->get_shock_Up(
                    Dat[BW::Q::iGamma43][it],
                    Dat[BW::Q::irho3][it],
                    Dat[BW::Q::iM3][it],
                    Dat[BW::Q::iEint3][it]);
        }

        /// check values
        if ((Dat[BW::Q::iadi][it] < 1.) || (Dat[BW::Q::iR][it] < 1.) || (Dat[BW::Q::irho2][it] < 0.) ||
            (Dat[BW::Q::iU_p][it] < 0) || (Dat[BW::Q::iGammaFsh][it] < 1.) ) {
            std::cerr << AT << " \n Wrong value at i=" << it << " tb=" << Dat[BW::Q::itburst][it] << "\n"
                      << " ----------------------------------------------------------- \n"
                      << " end_evolution=" << p_pars->end_evolution << "\n"
                      << " which_jet_layer_to_use=" << p_pars->which_jet_layer_to_use << "\n"
                      << " first_entry_r=" << p_pars->first_entry_r << "\n"
                      << " use_flat_dens_floor=" << p_pars->use_flat_dens_floor << "\n"
                      << " use_exp_rho_decay_as_floor=" << p_pars->use_exp_rho_decay_as_floor << "\n"
                      << " use_st_dens_profile=" << p_pars->use_st_dens_profile << "\n"
                      << " use_bm_dens_profile=" << p_pars->use_bm_dens_profile << "\n"
                      << " is_within0=  " << p_pars->is_within0 << "\n"
                      << " is_within=   " << p_pars->is_within << "\n"
                      << " ----------------------------------------------------------- \n"
                      << " ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << " ii_eq=" << p_pars->ii_eq
                      << " Mom0=" << p_pars->mom0 << " E0=" << p_pars->E0 << " theta0=" << p_pars->theta_b0
                      << " theta_max=" << p_pars->theta_max << " ncells=" << p_pars->ncells << "\n"
                      << " ----------------------------------------------------------- \n"
                      << " iR=          " << Dat[BW::Q::iR][it] << "\n"
                      << " itt=         " << Dat[BW::Q::itt][it] << "\n"
                      << " imom=        " << Dat[BW::Q::imom][it] << "\n"
                      << " iGamma=      " << Dat[BW::Q::iGamma][it] << "\n"
                      << " iGammaFsh=   " << Dat[BW::Q::iGammaFsh][it] << "\n"
                      << " iEint2=      " << Dat[BW::Q::iEint2][it] << "\n"
                      << " iEad2=       " << Dat[BW::Q::iEad2][it] << "\n"
                      << " iEsh2=       " << Dat[BW::Q::iEsh2][it] << "\n"
                      << " ictheta=     " << Dat[BW::Q::ictheta][it] << "\n"
                      << " irho=        " << Dat[BW::Q::irho][it] << "\n"
                      << " iEinj=       " << Dat[BW::Q::iEint2][it] << "\n"
                      << " idrhodr=     " << Dat[BW::Q::idrhodr][it] << "\n"
                      << " ibeta=       " << Dat[BW::Q::ibeta][it] << "\n"
                      << " iadi=        " << Dat[BW::Q::iadi][it] << "\n"
                      << " irho2=       " << Dat[BW::Q::irho2][it] << "\n"
                      << " ithickness=  " << Dat[BW::Q::ithickness][it] << "\n"
                      << " iU_p=        " << Dat[BW::Q::iU_p][it] << "\n"
                      << " imom=        " << Dat[BW::Q::imom][it] << "\n";
            exit(1);
        }
        if ((p_pars->m_type == BW_TYPES::iFS_DENSE) || (p_pars->m_type == BW_TYPES::iFS_PWN_DENSE)) {
            if ((Dat[BW::Q::iGammaCBM][it] < 1.)) {
                std::cerr << AT << " \n Wrong value at i=" << it << " tb=" << Dat[BW::Q::itburst][it] << "\n"
                     << " iGammaCBM=   " << Dat[BW::Q::iGammaCBM][it] << "\n"
                     << " iGammaREL=   " << Dat[BW::Q::iGammaREL][it] << "\n"
                     << " idGammaCBMdr=" << Dat[BW::Q::idGammaCBMdr][it] << "\n"
                     << " idGammaRELdGamma=" << Dat[BW::Q::idGammaRELdGamma][it] << "\n";
                exit(1);
            }
        }


        /// -------- energy injections
        if (p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            Dat[BW::Q::ipsrFrac][it] = p_pars->facPSRdep;
            Dat[BW::Q::iLmag][it] = p_pars->dEinjdt;
            Dat[BW::Q::iLnuc][it] = p_pars->dEnuc;
            ///
            double r_w = Dat[BW::Q::i_Rw][it];
            double u_b_pwn = 3.0 * Dat[BW::Q::i_Wepwn][it] / 4.0 / M_PI / r_w / r_w / r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
            double b_pwn = pow(u_b_pwn * 8.0 * M_PI, 0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
            Dat[BW::Q::i_Wb][it] = b_pwn;
            Dat[BW::Q::i_Wdr][it] = it > 0 ? Dat[BW::Q::i_Rw][it] - Dat[BW::Q::i_Rw][it - 1] : 0.;
        }
    }

/// ---------------------------------------------------------------------------------------------
public:
    void rhs_dispatcher(double * out_Y, size_t i, double x, double const * Y );

    /// right hand side of the blast wave evolution equation (ODE)
    void rhs_fs(double * out_Y, size_t i, double x, double const * Y );

    /// right hand side of the blast wave evolution equation (ODE)
    void rhs_fsrs(double * out_Y, size_t i, double x, double const * Y );

    /// set standard ISM for a BW
    void setStandardISM(double * out_Y, size_t i, double x, double const * Y){
//        double ej_Gamma  = EQS::GamFromMom( Y[i + SOL::QS::imom] );
        double ej_Gamma  = Y[i + SOL::QS::iGamma];
        if (ej_Gamma < 1.) {
            (*p_log)(LOG_ERR,AT) << "Gamma < 1\n";
            ej_Gamma = 1. + 1e-5;
        }
        double ej_R      = Y[i + SOL::QS::iR];
        double theta_b0  = p_pars->theta_b0;
        double ej_theta  = Y[i + SOL::QS::itheta];
        double ej_ctheta = EjectaID2::ctheta(ej_theta,p_pars->ilayer,p_pars->nlayers);//ctheta(ej_theta);
        set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
    }
    /// RHS with ODEs for dGammadR modified for pre-accelerated ISM
    void rhs_fs_dense(double * out_Y, size_t i, double x, double const * Y );
    /// RHS with ODEs for dGammadR modified for pre-accelerated ISM with energy injection
    void rhs_fs_dense_pwn(double * out_Y, size_t i, double x, double const * Y );

    void evalDensProfileInsideBWset(double * out_Y, size_t i, double x, double const * Y,
                                    std::vector<std::unique_ptr<BlastWave>> & others){

        if (m_data_shells[BW::QSH::iRs].size() < others.size()){
            for (auto &arr: m_data_shells)
                arr.resize(others.size());
        }

//        auto * p_others = (std::vector<std::unique_ptr<BlastWave>> *) others;
//        auto & othes = * p_others;
        double r = Y[p_pars->ii_eq + SOL::QS::iR];
//        double gam  = EQS::GamFromMom( Y[p_pars->ii_eq + SOL::QS::imom] );
        double gam  = Y[p_pars->ii_eq + SOL::QS::iGamma];
        double theta  = Y[p_pars->ii_eq + SOL::QS::itheta];
        double ctheta = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(ej_theta);
        /// Set default values
        set_standard_ism(r, theta, gam);

        /// Check if BW is behind or ahead of the BW system: find min/max of the whole BW system
        double rmin=std::numeric_limits<double>::max();
        double rmax=std::numeric_limits<double>::min();
        size_t nactive = 0;
        double _r=0.,rtmp = 0.;
#if 1
        for (size_t j = 0; (j < others.size()) && (!others[j]->getPars()->end_evolution); j++) {
            _r = Y[others[j]->getPars()->ii_eq + SOL::QS::iR];
            m_data_shells[BW::QSH::iRs][j] = Y[others[j]->getPars()->ii_eq + SOL::QS::iR];
            m_data_shells[BW::QSH::iGammas][j] = Y[others[j]->getPars()->ii_eq + SOL::QS::iGamma];//EQS::GamFromMom(Y[others[j]->getPars()->ii_eq + SOL::QS::imom]);

            if (_r > rmax) { rmax = _r; } // get maximum
            if (_r < rmin) { rmin = _r; } // get minimum
            if (rtmp < _r) { rtmp = _r; } // check if radii are ordered
            else{
                if (r <= _r) {
                    (*p_log)(LOG_ERR, AT) << "Radii not ordered! (Background ejecta started to decelerate)" << "r=" << _r << "\n";
                    exit(1);
                }
            }
            nactive+=1;
        }
        if (nactive < 3){
            (*p_log)(LOG_ERR,AT)<<"nactive < 3"<<"\n"; exit(1); }
        if (r < rmin){
            // BW is behind the ejecta system: Use exponential decay?
            return;
        }
        else if (r > rmax){
            // BW is ahead the ejecta system: Use exponential decay?
            return;
        }
        else{


            /// Find between which two BWs the given BW is
            size_t i_m=0;//index of the one on the left
            size_t i_p=0;//index of the one on the right
            for (size_t j=0;j<nactive-1;j++) {
                m_data_shells[BW::QSH::iDeltas][j] = m_data_shells[BW::QSH::iRs][j + 1] - m_data_shells[BW::QSH::iRs][j];
                /// find the BW from the set near current BW
                double _r_m = Y[others[j]->getPars()->ii_eq + SOL::QS::iR];
                auto & _x = others[j+1]->getPars();
                double _r_p = Y[others[j+1]->getPars()->ii_eq + SOL::QS::iR];
                if ((r>=_r_m)&&(r<=_r_p)){
                    i_m = j; i_p=j+1;
//                    break;
                } // found the two BW in between which this one lies
            }
            /// fill the ast point in the grid TODO: do this with interpolation...
            m_data_shells[BW::QSH::iDeltas][nactive-1] = m_data_shells[BW::QSH::iDeltas][nactive-2]*2.;
            /// check if we are in the grid
            if ((i_m==0)&&(i_p==0)){
                (*p_log)(LOG_ERR,AT)<<"No i_m and i_p found for BW with r="<<r<<" while rmin="<<rmin<<" and rmax="<<rmax<<"\n";
                exit(1);
            }

            for (size_t j=0;j<nactive;j++){
                m_data_shells[BW::QSH::iVols][j] = (4./3.)*CGS::pi*m_data_shells[BW::QSH::iRs][j]*m_data_shells[BW::QSH::iRs][j]*m_data_shells[BW::QSH::iDeltas][j]/(others[j]->getPars()->ncells);
                m_data_shells[BW::QSH::irhos][j] = (others[j]->getPars()->M0)/m_data_shells[BW::QSH::iVols][j];
                m_data_shells[BW::QSH::iPs][j] = Y[others[j]->getPars()->ii_eq + SOL::QS::iEint2]/m_data_shells[BW::QSH::iVols][j];
                if (m_data_shells[BW::QSH::irhos][j] <=0){
                    std::cerr<<m_data_shells[BW::QSH::iRs]<<"\n";
                    (*p_log)(LOG_ERR,AT)<<"j="<<j<<" vol="<<m_data_shells[BW::QSH::iVols][j]
                        <<" delta="<<m_data_shells[BW::QSH::iDeltas][j]
                        <<" rho="<<m_data_shells[BW::QSH::irhos][j]
                        <<"\n";
                    exit(1);
                }
            }


            /// interpolate density, gamma, pressure at the location of the current BW
            double rho = interpSegLog(i_m,i_p,r,m_data_shells[BW::QSH::iRs],m_data_shells[BW::QSH::irhos]);
            double GammaCMB = interpSegLog(i_m,i_p,r,m_data_shells[BW::QSH::iRs],m_data_shells[BW::QSH::iGammas]);
            double P = interpSegLog(i_m,i_p,r,m_data_shells[BW::QSH::iRs],m_data_shells[BW::QSH::iPs]);

            if ((GammaCMB < 1.)||(rho<=0)){
                (*p_log)(LOG_ERR,AT)<<"No i_m and i_p found for BW with GammaCMB="<<GammaCMB<<" rho="<<rho<<"\n";
                exit(1);
            }

            /// update these values inside the storage of the current BW
            double rho_def = p_dens->m_rho_def;
            double drhodr_def = p_dens->m_drhodr_def;

            double rho_prev = mDtmp[BW::Q::irho][0];
            double r_prev = mDtmp[BW::Q::iR][0];
            if (rho_prev == 0){ (*p_log)(LOG_ERR,AT)  << " rho[i] = 0" << "Exiting...\n";
                exit(1); }
            if (r == r_prev){ (*p_log)(LOG_ERR,AT) << AT << " R = R_i-1 \n";
                exit(1);
            }
            double dr = r - r_prev;
            double drhodr = (rho - rho_prev) / (r - r_prev);
//            m_tmp_data[BW::Q::iGamma][prev_ix+1] = gam;
//            m_tmp_data[BW::Q::irho][prev_ix+1] = rho;
//            m_tmp_data[BW::Q::iR][prev_ix+1] = r;
//            m_tmp_data[BW::Q::iPcbm][prev_ix + 1] = P;
//            double drhodr = dydx(m_tmp_data[BW::Q::iR], m_tmp_data[BW::Q::irho], m_tmp_data[BW::Q::idrhodr], dr, prev_ix+1, r, true);
            if (drhodr < 0){
//                (*p_log)(LOG_ERR,AT)<<" drhodr="<<drhodr<<"\n";
//                exit(1);
            }
            // if density is too small, for numerical reasons, kill the profile
//            if (rho < p_dens->m_rho_floor_val*rho_def){ rho = p_dens->m_rho_floor_val*rho_def; drhodr = 0.; }

            /// Compute the CBM density derivative (with radius)
            double dGammaRhodR = p_dens->m_dGammaRhodR;
            double GammaCBM_prev = mDtmp[BW::Q::iGammaCBM][0];//;getVal(BW::Q::iGammaCBM, (int)prev_ix - 1);
            if (GammaCMB == GammaCBM_prev ){
                dGammaRhodR = 0;
            }
            else {
                dGammaRhodR = (GammaCMB - GammaCBM_prev) / (r - r_prev);
//                m_tmp_data[BW::Q::iGammaCBM][prev_ix + 1] = GammaCMB;
//                dGammaRhodR = dydx(m_tmp_data[BW::Q::iR], m_tmp_data[BW::Q::iGammaCBM], m_tmp_data[BW::Q::idGammaCBMdr],
//                                   dr, prev_ix + 1, r, false);
            }
            if (!std::isfinite(dGammaRhodR)){
                (*p_log)(LOG_ERR,AT)  << " Nan dGammaRhodR . setting to 0 ...";
                dGammaRhodR = 0.;
            }

            /// relative velocity of this BW and BW system
            double GammaREL = GammaCMB;
            GammaREL = EQS::GammaRel(gam,GammaCMB);
            if (GammaREL < 1.) GammaREL = 1.;
            if ((GammaREL < 1.) || (!std::isfinite(GammaREL)) ){
                (*p_log)(LOG_ERR,AT)  << " bad GammaREL=" << GammaREL << " (too big for ejecta or nan)"
                                      << " Gamma0=" << p_pars->Gamma0
                                      << " E0=" << p_pars->E0
                                      << " M0=" << p_pars->M0
                                      << " rho=" << rho
                                      << " rho_prev=" << rho_prev
                                      << " GammaCMB=" << GammaCMB
                                      << " GammaCBM_prev=" << GammaCBM_prev
                                      << " Exiting..." << "\n";
//            std::cerr << "| Gamma evol: \n";
//            std::cerr << m_data[Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
                exit(1);
            }

            double adi = p_eos->getGammaAdi(GammaCMB,EQS::BetaFromGamma(GammaCMB));//5/3.;
            double cscbm = sqrt(adi * P / rho) / CGS::c; // sound speed
            double mcbm = sqrt(rho * EQS::BetaFromGamma(GammaREL) * EQS::BetaFromGamma(GammaREL) * CGS::c * CGS::c / adi / P); // mach number

//            double ej_Gamma_prev = getVal(BW::Q::iGamma, (int)prev_ix - 1);
//            double GammaCMB_prev = getVal(BW::Q::iGammaCBM, (int)prev_ix - 1);
//            double GammaREL_prev = EQS::GammaRel(ej_Gamma_prev, GammaCMB_prev);
//            if ((gam == ej_Gamma_prev) || (GammaREL_prev < 1.) || (GammaCMB_prev == 1.)){
//                dGammaRELdGamma = 1.;
//            }
//            else{
//            m_tmp_data[BW::Q::iGammaREL][prev_ix+1] = GammaREL;
            double GammaRel_prev = mDtmp[BW::Q::iGammaREL][0];
            double Gamma_prev = mDtmp[BW::Q::iGamma][0];
            double dGammaRELdGamma = (GammaREL - GammaRel_prev) / (gam - Gamma_prev);
//            double dGammaRel = m_tmp_data[BW::Q::iGammaREL][prev_ix+1] - m_tmp_data[BW::Q::iGammaREL][prev_ix];
//            if (m_tmp_data[BW::Q::iGammaREL][prev_ix]<0)
//                dGammaRel = 0;
//            dGammaRELdGamma = dydx(m_tmp_data[BW::Q::iGammaREL], m_tmp_data[BW::Q::iGamma], m_tmp_data[BW::Q::idGammaRELdGamma],
//                                       dr, prev_ix+1, GammaREL, false);
//            dPdrho = dydx(m_tmp_data[BW::Q::iPcbm], m_tmp_data[BW::Q::irho], m_tmp_data[BW::Q::idPCBMdrho],
//                              dr, prev_ix+1, P, false);
            double P_prev = mDtmp[BW::Q::iPcbm][0];
            double dPdrho = (P - P_prev) / (rho - rho_prev);
            double tmp = sqrt(P / rho) / CGS::c;
//                int x = 1;
//            dGammaRELdGamma1 = (GammaREL - GammaREL_prev) / (ej_Gamma - ej_Gamma_prev);
//            }
            if ((!std::isfinite(dGammaRELdGamma)) || (dGammaRELdGamma > 1.e2)){
                (*p_log)(LOG_ERR,AT) << "HUGE dGammaRELdGamma="<<dGammaRELdGamma<<"\n";
                dGammaRELdGamma = 1.;
            }
            if (dGammaRELdGamma < 0){
//                (*p_log)(LOG_ERR,AT) << "NEGATIVE dGammaRELdGamma="<<dGammaRELdGamma<<"\n";
//                dGammaRELdGamma = 1.;
            }

            // Finally, apply the computed density profile to be used in RHS
//            if( (GammaCMB >= gam) ){
//            std::cerr << AT << " GammaCMB="<<GammaCMB<<" > ej_Gamma="<<ej_Gamma<<"\n";
//                GammaCMB = 1.; GammaREL = gam; dGammaRhodR = 0; dGammaRELdGamma = 1.; //rho = 1e-70; m_drhodr = 0.;
//            }

            p_dens->m_rho_ = rho;
            p_dens->m_drhodr_ = drhodr;
            p_dens->m_GammaRho = GammaCMB;
            p_dens->m_GammaRel = GammaREL;
            p_dens->m_dGammaRhodR = dGammaRhodR;
            p_dens->m_dGammaReldGamma = dGammaRELdGamma;
            p_dens->m_dPCBMdrho = dPdrho;
            p_dens->m_P_cbm = P;
            p_dens->m_M_cbm = mcbm; // 4.47 Zhang:2018
            p_dens->m_CS_CBM = cscbm; // 4.43 Zhang:2018
//#endif
        }
#endif
    }

    /// --- Evaluate the density (and its velocity) at point ej_R left by blast wave at a point j_R
    void evalDensAndItsVelocityBehindBlastWave(double & rho, double & GammaCMB, double & p_cbm,
                                               double rho_def, double rho_prev,
                                               double ej_R, double j_R, double j_rho, double j_rho2,
                                               double j_Gamma, double j_Gamma0, double j_P2){

        (*p_log)(LOG_ERR,AT)<<" IMPLEMENTATIONS IS NOT READY \n";
        exit(1);
#if 0
        GammaCMB = 1.;
        /// exponential decay from the point of entry
        double rho0 = getVal(BW::Q::irho, 0);
        /// if inside from the it=0 -- keep floor dens untill BM and ST profiles
        if ((rho0==(p_dens->m_rho_floor_val*rho_def))
            &&(rho_prev==(p_dens->m_rho_floor_val*rho_def))){
            rho = rho_prev;
        }
        else if (p_pars->use_exp_rho_decay_as_floor) {
            double coeff = p_pars->steepnes_of_exp_decay;
            double scale_rho = exp(1. - coeff * ej_R / getPars()->first_entry_r);
            double scale_dlnrhodr = -1. * coeff * scale_rho / getPars()->first_entry_r;
            if (scale_rho > 1.) {
                (*p_log)(LOG_ERR,AT)  << " scale_rho > 1 for ej_R=" << ej_R
                                      << "  entry_r=" << getPars()->first_entry_r << "\n"
                                      << " Exiting...";
//                std::cerr << AT  << "\n";
                exit(1);
            }
            rho = rho_def * scale_rho;
//            if (rho_prev < rho){
//                std::cerr <<AT<<"\n";
//                exit(1);
//            }
//            dlnrhodr = rho_def * scale_dlnrhodr;// analytc. derivative
        }
        else if ((p_pars->use_flat_dens_floor)){
            rho = rho_def * p_dens->m_rho_floor_val;
//            dlnrhodr = 0.;
        }

        /// check
        if (!std::isfinite(rho) || rho < 0. || rho > 1. || rho > rho_def ){
            (*p_log)(LOG_ERR,AT) << " wrong floor rho="<<rho<<"\n" << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        /// Blandford-Mckee tail from the jet blast wave (j_Gamma > GammaST and j_Gamma < GammaBM)
        double rho_bm = -1;
        double gam_cbm_bm = -1;
        if((ej_R > j_R * 0.01)&&(p_pars->use_bm_dens_profile)
           &&(j_Gamma > p_pars->Gamma_when_st_starts)
           &&(j_Gamma < (j_Gamma0 * p_pars->fraction_of_Gamma0_when_bm_for_bm))) {
            rho_bm = getBM()->rho_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2) * CGS::mp;
            gam_cbm_bm = getBM()->gamma_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2);
//            p_cbm_s = getSedov()->pressure_profile_int(ej_R, j_R, EQS::Beta(j_P2));
            (*p_log)(LOG_ERR,AT)  << " pressure is not implemented\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
            /// apply only of above the floor
            if (rho_bm > rho) { rho = rho_bm; GammaCMB = gam_cbm_bm; }
        }

        /// Sedov-Taylor tail from the jet blast wave (j_Gamma < GammaST)
        double rho_s = -1;
        double gam_cbm_s = -1;
        double p_cbm_s = -1.;
        if ((ej_R > j_R * 0.01) && // TODO this is needed to activate ST when kB BW is close to jet BW
            (p_pars->use_st_dens_profile)
            &&(j_Gamma <= p_pars->Gamma_when_st_starts)) {
            rho_s = getSedov()->rho_profile_int(ej_R, j_R, j_rho2);
            gam_cbm_s = getSedov()->Gamma_profile_int(ej_R, j_R, EQS::Beta(j_Gamma));
            p_cbm_s = getSedov()->pressure_profile_int(ej_R, j_R, j_P2);
            // apply only of above the floor
//            if ((rho_bm > rho_s)&&(rho_bm > rho)){
//                rho = rho_bm; GammaCMB = gam_cbm_bm;
//            }
//            else
            if (rho_s > rho) { // # THis is needed for some reason :)
                p_pars->is_using_st_prof = true;
                rho = rho_s; GammaCMB = gam_cbm_s; p_cbm = p_cbm_s;
            }


//            if (rho_prev < 1e-10 * rho){
//                std::cerr << AT << " rho_prev="<<rho_prev<<" rho="<<rho<<" rho0="<<rho0<<" density gradient >10 orders of magnitude\n";
////                exit(1);
//            }
        }



        if ((GammaCMB < 1.) || (!std::isfinite(GammaCMB))){
            (*p_log)(LOG_ERR,AT)  << " bad GammaCMB=" << GammaCMB << "\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        if ((((std::abs(rho_prev/rho) < 1e-4)||((std::abs(rho_prev/rho) > 1e4))))){
            //std::cerr << AT << " rho_prev="<<rho_prev<<" rho="<<rho<<" rho0="<<rho0<<" density gradient >4 orders of magnitude\n";
//            exit(1);
        }
#endif

    }

    // ---------------------------------------------------------

    /// Density profiles application based on the jet BW position and velocity
    void evalDensAndItsVelocityBehindBlastWave_Case1(
            double j_R, double j_Gamma, double j_ctheta, double j_rho, double j_rho2, double j_Gamma0, double j_P2,
            double ej_R, double ej_Gamma, double ej_theta ){
#if 0
        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);

        if (ej_Gamma <= 1.) { ej_Gamma = 1. + 5e-5; } // TODO REMOVE

//        size_t prev_ix = p_pars->comp_ix;

        // --- evaluate what density profile the ejecta blast wave would experience
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta); // set default values for density
        double rho_prev = getVal(BW::Q::irho, prev_ix-1);
        double rho_def = p_dens->m_rho_def;
        double drhodr_def = p_dens->m_drhodr_def;
        if (rho_def<(p_dens->m_rho_floor_val*rho_def)){
            (*p_log)(LOG_ERR,AT)  << " rho_def is not initialized\n Exiting...\n";
//            std::cerr<<AT<<"\n";
            exit(1);
        }
//        if ((rho_prev/p_dens->m_rho_) >= p_dens->m_rho_floor_val){
//            std::cerr << AT << " " << p_pars->is_within0 << " and "<< p_pars->is_within <<"\n";
//            std::cerr<<AT<<" rho_def is not initialized\n";
//            exit(1);
//        }
        double rho=rho_def, GammaCMB=1, P=0, GammaREL=ej_Gamma, drhodr=drhodr_def, dGammaRELdGamma=1, dGammaRhodR=0, dPdrho=0., cscbm=1., mcbm=0.;
        if (rho_prev<(p_dens->m_rho_floor_val*rho_def)){
            (*p_log)(LOG_ERR,AT) <<" rho_prev is not initialized\n Exiting...\n";
//            std::cerr<<AT<<"\n";
            exit(1);
        }
        if (rho_prev == 0){
            (*p_log)(LOG_ERR,AT)  << " rho[i] = 0" << "\n Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        /// eval dens seen by ejecta BW
        evalDensAndItsVelocityBehindBlastWave( rho, GammaCMB, P, rho_def, rho_prev, ej_R,  j_R,
                                               j_rho,  j_rho2,  j_Gamma, j_Gamma0, j_P2);

        /// Compute the density derivative (with radius)
        double r_prev = getVal(BW::Q::iR, prev_ix - 1);
//        double rho_prev = getVal(Q::irho, evaled_ix - 1);
        if (ej_R == r_prev){ (*p_log)(LOG_ERR,AT) << AT << " R = R_i-1 \n"; exit(1); }
        m_data[BW::Q::irho][prev_ix+1] = rho;
        m_data[BW::Q::iR][prev_ix+1] = ej_R;
        m_data[BW::Q::iPcbm][prev_ix + 1] = P;
        double dr = m_data[BW::Q::iR][prev_ix] - m_data[BW::Q::iR][prev_ix-1];
//            double dr1 = m_data[Q::iR][0] - m_data[Q::iR][1];
        drhodr = dydx(m_data[BW::Q::iR], m_data[BW::Q::irho], m_data[BW::Q::idrhodr],
                      dr, prev_ix+1, ej_R, true);
//        double drhodr2 = (rho - rho_prev) / (ej_R - r_prev);
        if (rho < p_dens->m_rho_floor_val*rho_def){ rho = p_dens->m_rho_floor_val*rho_def; drhodr = 0.; }
        if (!std::isfinite(rho) || !std::isfinite(drhodr)) {
            (*p_log)(LOG_ERR,AT)  <<" bad rho="<<rho<<" or m_drhodr="<<drhodr<<" \n Exiting...\n";
//            std::cerr << AT <<" \n";
            exit(1);
        }

        /// Compute the CBM density derivative (with radius)
        double GammaCBM_prev = getVal(BW::Q::iGammaCBM, prev_ix - 1);
        if (GammaCMB == GammaCBM_prev ){
            dGammaRhodR = 0;
        }
        else{
            m_data[BW::Q::iGammaCBM][prev_ix+1] = GammaCMB;
            dGammaRhodR = dydx(m_data[BW::Q::iR], m_data[BW::Q::iGammaCBM], m_data[BW::Q::idGammaCBMdr],
                               dr, prev_ix+1, ej_R, false);

//            double dGammaRhodR1 = (GammaCMB - GammaCBM_prev) / (ej_R - r_prev);
//            if (dGammaRhodR != dGammaRhodR1){
//                exit(1);
//            }
        }
        if (!std::isfinite(dGammaRhodR)){
            (*p_log)(LOG_ERR,AT)  << " Nan dGammaRhodR . setting to 0\n Exiting...";
//            std::cerr << AT << "\n";
            dGammaRhodR = 0.;
        }

        /// computeSynchrotronEmissivityAbsorptionAnalytic the relative velocity of the ISM (with respect to ejecta) and its derivative
//            double betaREL = EQS::BetaRel(EQS::Beta(ej_Gamma),EQS::Beta(GammaCMB));
//            GammaREL = EQS::Gamma(betaREL);
        GammaREL = EQS::GammaRel(ej_Gamma,GammaCMB);
        if (GammaREL < 1.) GammaREL = 1.;
        if ((GammaREL < 1.) || (!std::isfinite(GammaREL)) ){
            (*p_log)(LOG_ERR,AT)  << " bad GammaREL=" << GammaREL << " (too big for ejecta or nan)" << "\n"
                                  << " j_R=" << j_R << "\n"
                                  << " ej_R=" << ej_R << "\n"
                                  << " ej_theta=" << ej_theta << "\n"
                                  << " j_Gamma=" << j_Gamma << "\n"
                                  << " ej_Gamma=" << ej_Gamma << "\n"
                                  << " Gamma0=" << p_pars->Gamma0 << "\n"
                                  << " E0=" << p_pars->E0 << "\n"
                                  << " M0=" << p_pars->M0 << "\n"
                                  << " rho=" << rho << "\n"
                                  << " rho_prev=" << rho_prev << "\n"
                                  << " GammaCMB=" << GammaCMB << "\n"
                                  << " GammaCBM_prev=" << GammaCBM_prev << "\n"
                                  << " Exiting...";
//            std::cerr << "| Gamma evol: \n";
//            std::cerr << m_data[Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        double adi = p_eos->getGammaAdi(GammaCMB,EQS::Beta(GammaCMB));//5/3.;
        cscbm = sqrt(adi * P / rho) / CGS::c; // sound speed
        mcbm = sqrt(rho * EQS::Beta(GammaREL) * EQS::Beta(GammaREL) * CGS::c * CGS::c / adi / P); // mach number

        double ej_Gamma_prev = getVal(BW::Q::iGamma, prev_ix - 1);
        double GammaCMB_prev = getVal(BW::Q::iGammaCBM, prev_ix - 1);
        double GammaREL_prev = EQS::GammaRel(ej_Gamma_prev, GammaCMB_prev);
        if ((ej_Gamma == ej_Gamma_prev) || (GammaREL_prev < 1.) || (GammaCMB_prev == 1.)){
            dGammaRELdGamma = 1.;
        }
        else{
            m_data[BW::Q::iGammaREL][prev_ix+1] = GammaREL;
            double dGammaRel = m_data[BW::Q::iGammaREL][prev_ix] - m_data[BW::Q::iGammaREL][prev_ix-1];
            dGammaRELdGamma = dydx(m_data[BW::Q::iGammaREL], m_data[BW::Q::iGamma], m_data[BW::Q::idGammaRELdGamma],
                                   dr, prev_ix+1, GammaREL, false);
            dPdrho = dydx(m_data[BW::Q::iPcbm], m_data[BW::Q::irho], m_data[BW::Q::idPCBMdrho],
                          dr, prev_ix+1, P, false);
            double tmp = sqrt(P / rho) / CGS::c;
            int x = 1;
//            dGammaRELdGamma1 = (GammaREL - GammaREL_prev) / (ej_Gamma - ej_Gamma_prev);
        }
        if ((!std::isfinite(dGammaRELdGamma)) || (dGammaRELdGamma > 1.e2)){
//                std::cerr << AT << " Nan dGammaRELdGamma . setting to 1\n";
            dGammaRELdGamma = 1.;
        }

        // Finally, apply the computed density profile to be used in RHS
        if( (GammaCMB >= ej_Gamma) ){
//            std::cerr << AT << " GammaCMB="<<GammaCMB<<" > ej_Gamma="<<ej_Gamma<<"\n";
            GammaCMB = 1.; GammaREL = ej_Gamma; dGammaRhodR = 0; dGammaRELdGamma = 1.; //rho = 1e-70; m_drhodr = 0.;
        }

        p_dens->m_rho_ = rho;
        p_dens->m_drhodr_ = drhodr;
        p_dens->m_GammaRho = GammaCMB;
        p_dens->m_GammaRel = GammaREL;
        p_dens->m_dGammaRhodR = dGammaRhodR;
        p_dens->m_dGammaReldGamma = dGammaRELdGamma;
        p_dens->m_dPCBMdrho = dPdrho;
        p_dens->m_P_cbm = P;
        p_dens->m_M_cbm = mcbm; // 4.47 Zhang:2018
        p_dens->m_CS_CBM = cscbm; // 4.43 Zhang:2018
#endif

    }
    void set_no_ism(double ej_R, double ej_ctheta, double ej_Gamma_rel){
//        if (ej_Gamma_rel>10){ std::cerr << AT << "\n"; exit(1); }
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta);
        p_dens->m_rho_ = p_dens->m_rho_floor_val * p_dens->m_rho_def;
        p_dens->m_drhodr_ = 0.;

        p_dens->m_GammaRho = 1;
        p_dens->m_GammaRel = ej_Gamma_rel;
        p_dens->m_dGammaRhodR = 0;
        p_dens->m_dGammaReldGamma = 1;
        p_dens->m_dPCBMdrho = 0.;
        p_dens->m_P_cbm = 0.;
        p_dens->m_M_cbm = 0.;
        p_dens->m_CS_CBM = 0;
    }
    void set_standard_ism(double ej_R, double ej_ctheta, double ej_Gamma_rel){
//        if (ej_Gamma_rel>10){ std::cerr << AT << "\n"; exit(1); }
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta);
        p_dens->m_rho_ = p_dens->m_rho_def;
        p_dens->m_drhodr_ = p_dens->m_drhodr_def;
        p_dens->m_GammaRho = 1;
        p_dens->m_GammaRel = ej_Gamma_rel;
        p_dens->m_dGammaRhodR = 0;
        p_dens->m_dGammaReldGamma = 1;
        p_dens->m_dPCBMdrho = 0.;
        p_dens->m_P_cbm = 0.;
        p_dens->m_M_cbm = 0.;
        p_dens->m_CS_CBM = 0;
    }
    void prepareDensProfileFromJet(double * out_Y, size_t i, double x, double const * Y,
                                   std::vector<std::unique_ptr<BlastWave>> & others){

        (*p_log)(LOG_ERR,AT)<<" IMPLEMENTATIONS IS NOT READY \n";
        exit(1);
#if 0

//        auto * p_pars = (struct PWNPars *) params; // removing EATS_pars for simplicity
//        auto * p_others = (std::vector<std::unique_ptr<BlastWave>> *) _others;
//        auto & others = * p_others;
        // --- current ejecta bw
//        double ej_Gamma  = Y[p_pars->ii_eq + DynRadBlastWave::QS::iGamma];
        double ej_mom = Y[p_pars->ii_eq + SOL::QS::imom];
        if (ej_mom < 0.){
            (*p_log)(LOG_ERR,AT)<<" ej_mom < 0\n";
            exit(1);
        }
        double ej_Gamma  = EQS::GamFromMom(ej_mom);
//        double ej_Gamma = std::exp(ej_lnGamma);
        if (ej_Gamma < 1.) ej_Gamma = 1.0000000099999;
//        if (ej_Gamma > 10.) {
//            std::cerr << AT << " ej_Gamma="<<ej_Gamma<<"\n";
//            exit(1);
//        }
        double ej_R      = Y[p_pars->ii_eq + SOL::QS::iR];
        double theta_b0  = p_pars->theta_b0;
        double ej_theta  = Y[p_pars->ii_eq + SOL::QS::itheta];
        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);
        // -- loop over jets

        int i_ej_l = p_pars->which_jet_layer_to_use;//others.size()-1;
        bool is_within = false;
//        for (size_t ij = 0; ij < others.size(); ij++){
//            double j_R      = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            double j_Gamma  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//            double j_theta  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::QS::itheta];
//            double j_ctheta = others[ij]->ctheta(j_theta);
//            double j_theta0 = others[ij]->theta0(j_theta); //
//            double j_theta1 = others[ij]->theta1(j_theta);
//            if ((ej_ctheta > j_theta0) && (ej_ctheta <= j_theta1)){
////            if ((ej_ctheta < j_ctheta)){
//                i_ej_l = ij; is_within = true;
////                if(i_ej_l!=p_pars->j_i0){
////                    std::cerr << AT<< " \nj_theta0"<<p_pars->j_theta0<<" j_theta1="<<p_pars->j_theta1<<" prev_layer="<<p_pars->j_i0<<"\n";
////                    std::cerr << "ctheta="<<j_theta0<<" ctheta0="<< p_pars->ctheta0 <<"\n";
////                    exit(1);
////                }
//            }
//        }
//        if (!is_within){
//            double j_R      = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            double j_Gamma  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//            double j_theta  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::QS::itheta];
//            double j_ctheta = others[0]->ctheta(j_theta);
//            double j_theta0 = others[0]->theta0(j_theta); //
//            double j_theta1 = others[0]->theta1(j_theta);
//            if (ej_ctheta < j_theta0){
//                i_ej_l = 0.;
//                is_within = true;
//            }
//        }

//        i_ej_l = 0;//others.size()-1;
        double j_mom = Y[others[i_ej_l]->getPars()->ii_eq + SOL::QS::imom];
        double j_Gamma = EQS::GamFromMom(j_mom);
//        double j_Gamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::iGamma];
//        double ibeta = EQS::BetFromMom(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::imom]);
        double j_theta = Y[others[i_ej_l]->getPars()->ii_eq + SOL::QS::itheta];
//        double j_lnGamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::ilnGamma];
//        double j_Gamma = std::exp(j_lnGamma);
        double j_Gamma0 = others[i_ej_l]->getPars()->Gamma0;
//        double j_ctheta = others[i_ej_l]->ctheta(j_theta);
        double j_ctheta = EjectaID2::ctheta(j_theta,
                                            others[i_ej_l]->getPars()->ilayer,
                                            others[i_ej_l]->getPars()->nlayers);//others[i_ej_l]->ctheta(j_theta);
        if (ej_ctheta < j_ctheta){ is_within = true; }

        auto & other = others[i_ej_l];
        double j_R = Y[other->getPars()->ii_eq + SOL::QS::iR];

        if ((ej_R == p_pars->R0)||(evaled_ix == 0)){
//            if (evaled_ix > 1){
//                std::cerr<<AT<<" \n";exit(1);
//            }
            if (is_within){
                // inside -- NO ISM
                set_no_ism(ej_R, ej_ctheta, ej_Gamma);
//                p_pars->is_within0 = is_within;
            }
            else {
                // outside -- standard ISM
                set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
//                p_pars->is_within0 = is_within;
            }
//            if ((p_pars->ilayer==1)&&(p_dens->m_rho_!=(p_dens->m_rho_floor_val*p_dens->m_rho_def))){
//                exit(1);
//            }
        }
        else if ((ej_R < j_R)&&(ej_R > p_pars->R0)){
            if (is_within && (evaled_ix > 0)){
                // inside and behind -- complex ISM
                if ((p_pars->first_entry_r < 0) || (p_pars->first_entry_r > ej_R))
                    p_pars->first_entry_r = ej_R;
                // if entering a new jet layer
                if (i_ej_l != p_pars->ijl) {
                    p_pars->prev_ijl = p_pars->ijl; // -1 0 1 2 3 4
                    p_pars->ijl = i_ej_l; // 0 2 3 4 5
                }

                size_t other_i = others[i_ej_l]->getPars()->ii_eq;//other->getPars()->ii_eq;
                double other_Gamma = EQS::GamFromMom(Y[other_i + SOL::QS::imom]);//std::exp( Y[other_i + DynRadBlastWave::QS::ilnGamma] );
//                double other_Gamma = Y[other_i + DynRadBlastWave::QS::iGamma];//std::exp( Y[other_i + DynRadBlastWave::QS::ilnGamma] );
                // --- density experienced by the jet blast wave and the density in the downstream (rho2)
                other->getDensIsm()->evaluateRhoDrhoDrDefault(j_R, INFINITY);
                double j_rho = other->getDensIsm()->m_rho_def;
                double j_drhodr = other->getDensIsm()->m_drhodr_def;
                double j_adi = other->getEos()->getGammaAdi(other_Gamma, EQS::Beta(other_Gamma));
                double j_rho2 = EQS::rho2t(j_Gamma, j_adi, j_rho);
                double j_V2 = Y[other_i + SOL::QS::iM2] / j_rho2; // Units -> c^2 for energy
                double j_P2 = (j_adi - 1.) * Y[other_i + SOL::QS::iEint2] / j_V2 * CGS::c * CGS::c;//# * CGS::c * CGS::c; // Units -> c^2 for energy
                double cs = sqrt(j_adi * j_P2 / j_rho2) / CGS::c;
                double cs2 = sqrt(j_adi*(j_adi-1)*(j_Gamma-1)/(1+j_adi*(j_Gamma-1)));
                // ---
//                size_t i = p_pars->ii_eq;
//                double ej_Gamma = Y[i + DynRadBlastWave::QS::iGamma];
//                double ej_R = Y[i + DynRadBlastWave::QS::iR];
//                double theta_b0 = p_pars->theta_b0;
//                double ej_theta = Y[i + DynRadBlastWave::QS::itheta];
//                double ej_ctheta = ctheta(ej_theta);
                if (ej_Gamma <= 1.) { ej_Gamma = 1.0000000000999; } // --------------------------- [ kostil ]

                evalDensAndItsVelocityBehindBlastWave_Case1(j_R, j_Gamma, j_ctheta, j_rho, j_rho2, j_Gamma0, j_P2, ej_R, ej_Gamma, ej_theta);
//                evalDensAndItsVelocityBehindBlastWave_Case1( Y, other );
            }
            else{
                // outside -- standard ISM
                set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
            }
        }
        else if (ej_R >= j_R){
            // outside -- standard ISM
            set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
        }
        else{
            (*p_log)(LOG_ERR,AT)  << " should not be entered. Exiting..."<<"\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

//        if ((p_pars->is_within)&&(!is_within)){
//            exit(1);
//        }
//        double rho =  p_dens->m_rho_;
//        p_dens->evaluateRhoDrhoDr(ej_R, ej_ctheta); // set default values for density
//        double rho_def = p_dens->m_rho;
//        double drhodr_def = p_dens->m_drhodr;
        double rho_prev    = getVal(BW::Q::irho, evaled_ix-1);
        double drhodr_prev = getVal(BW::Q::idrhodr, evaled_ix-1);
        double prev_ej_R   = getVal(BW::Q::iR, evaled_ix-1);
//        if ((rho<1e-30)&&(evaled_ix>1)&&(rho_prev<rho)){
//            std::cerr<<AT<<" rho_prev is not initialized\n";
//            exit(1);
//        }
        p_pars->is_within = is_within;
        p_pars->r_dist = j_R - ej_R;

//        if (p_pars-> < pfsolvePars->dens_floor_frac*rho_def)
        /// evaluate the RHS

//        if ((p_pars->prev_ijl!=i_ej_l)&&(evaled_ix>20)){
//            std::cerr << " ishell="<<p_pars->ishell
//                      << " ilayer="<<p_pars->ilayer
//                      << " Gamma0="<<p_pars->Gamma0
//                      << " Gamma="<<ej_Gamma
//                      << " rho_prev"<<rho_prev
//                      << " drhodr_prev"<<drhodr_prev
//                      << " rho="<<p_dens->m_rho
//                      << " drhodr="<<p_dens->m_drhodr
//                      << " ej_ctheta="<<ej_ctheta
//                      << " ej_ctheta_prev="<<getVal(Q::ictheta, evaled_ix-1)
//                      << " j_theta0="<<others[ijl]->theta0(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::itheta])
//                      << " j_theta1="<<others[ijl]->theta1(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::QS::itheta])
//                      << "\n";
//        }

//        if (((rho_prev/p_dens->m_rho_) >= p_dens->m_rho_floor_val)){
//            exit(1);
//        }
//        if ((p_pars->ilayer==1)&&(p_dens->m_rho_!=(p_dens->m_rho_floor_val*p_dens->m_rho_def))){
//            exit(1);
//        }
//        if (p_pars->is_within0 != p_pars->is_within){
//            exit(1);
//        }
#endif
    }
};


/// select an appropriate RHS
void BlastWave::rhs_dispatcher(double * out_Y, size_t i, double x, double const * Y ){
    switch (p_pars->m_type) {
        case iFS:
            rhs_fs(out_Y, i, x, Y);
            break;
        case iFSRS:
            rhs_fsrs(out_Y, i, x, Y);
            break;
        case iFS_DENSE:
            rhs_fs_dense(out_Y, i, x, Y);
            break;
        case iFS_PWN_DENSE:
            rhs_fs_dense_pwn(out_Y, i, x, Y);
            break;
    }
}


/// right hand side of the blast wave evolution equation (ODE) with FS only
void BlastWave::rhs_fs(double * out_Y, size_t i, double x, double const * Y ) {
    // ****************************************
    double R      = Y[i+ SOL::QS::iR];
    double Rsh    = Y[i+ SOL::QS::iRsh];
    double tcomov = Y[i+ SOL::QS::itcomov];
//        double mom    = Y[i+ SOL::QS::imom];
    double Gamma  = Y[i+ SOL::QS::iGamma];
    double Eint2  = Y[i+ SOL::QS::iEint2];
    double theta  = Y[i+ SOL::QS::itheta];
    double M2     = Y[i+ SOL::QS::iM2];
    // ****************************************
    if (Gamma < 1 || M2 < 0 || theta < 0) {
        (*p_log)(LOG_ERR, AT) << "Wrong value in RHS: Gamma="<<Gamma
            << " Gamma0="<<p_pars->Gamma0
            << " E0=" <<p_pars->E0
            << " M0=" <<p_pars->M0
            << " R=" <<R
//            << " Rsh=" <<Rsh
            << " Eint2=" <<Eint2
            << " theta=" <<theta
            << " M2=" <<M2
            << " rho=" <<p_dens->m_rho_def / p_pars->M0
            << " drhodr=" <<p_dens->m_drhodr_def / p_pars->M0
            << " M2=" <<M2
            << " dthetadr_prev=" << p_pars->dthetadr_prev
            <<"\n";
        exit(1);
    }
//        double Gamma = EQS::GamFromMom(mom);
//        if (Gamma <= 1.)
//            Gamma = EQS::GamFromMom(mom);
//            Gamma = 1.0000000001;
//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }

    // ****************************************
//        auto *_pars = (struct RHS_pars *) rhs_pars;
    double ctheta_ = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);//p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
    p_dens->evaluateRhoDrhoDrDefault(R, ctheta_);
    double rho = p_dens->m_rho_def / p_pars->M0;
    double drhodr = p_dens->m_drhodr_def / p_pars->M0;

//        dlnrho1dR /= p_pars->M0;
//        double _rhoi = p_pars->rho/p_pars->M0;
    // ****************************************
    double beta   = EQS::BetaFromGamma(Gamma);
//    double dRdt = betaSh * CGS::c; // TODO this is beta_sh ! ...
    double gammaAdi  = p_eos->getGammaAdi(Gamma, beta);
    double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
//    long double betaSh = EQS::Beta2(GammaSh);
//    long double dRshdt = betaSh * CGS::c;
////    double dRdtcomov = 1.0 / betaSh / Gamma / CGS::c;
//    double dRshdt_ = betaSh / (Gamma * Gamma * (1. - betaSh * betaSh)) * CGS::c;// / Gamma;
//    std::cout << "dr= " << dRshdt << " !/ " << dRshdt_ << "\n";
    double dRshdt = EQS::dRsh_dt(Gamma, beta, GammaSh);
    double dRdt = beta * CGS::c;

//    std::cout << " drdt = " << dRdt << " drsh_dt = " << dRshdt << "\n";

    double dthetadr = 0.0;
    if ((theta < p_pars->theta_max) && (!p_pars->end_spreading)) {
        if (theta < 0){
            (*p_log)(LOG_ERR,AT) << " theta="<<theta<<"\n";
            exit(1);
        }
        switch (p_pars->method_limit_spread) {

            case iNone:
                dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iGamma0Frac:
                if (Gamma*beta < p_pars->fraction_of_mom0_when_spread_start*p_pars->mom0)
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iGammaVal:
                if (Gamma*beta < std::max(p_pars->value_of_mom_when_spread_start,
                                          p_pars->fraction_of_mom0_when_spread_start*p_pars->mom0))
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iRd:
                if (R > p_pars->Rd)
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
        }
    }
    p_pars->dthetadr_prev = dthetadr; // for debugging


    double dM2dR = 0;
    switch (p_pars->m_method_dmdr) {

        case iusingA:
            dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
            break;
        case iusingdthdR:
            dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
            break;
        case iNodmdr:
            break;
    }
    if (dM2dR < 0.){
        (*p_log)(LOG_ERR, AT) << " dM/dR < 0\n";
        exit(1);
    }

    // --- dGammadR ---
    double dGammadR = 0.,GammaEff=0.,dGammaEffdGamma=0., num=0., denum=0.;
    switch (p_pars->m_method_dgdr) {

        case iour:
            GammaEff = get_GammaEff(Gamma, gammaAdi); // gammaAdi) #(gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
            dGammaEffdGamma = get_dGammaEffdGamma(Gamma, gammaAdi); // 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.
            num = ((Gamma - 1.) * (GammaEff + 1.) * dM2dR - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR / M2 - drhodr / rho));
            denum = ((1. + M2) + Eint2 * dGammaEffdGamma + GammaEff * (gammaAdi - 1.) * Eint2 / Gamma);
            dGammadR = -num/denum;
            break;
        case ipeer:
            dGammadR = EQS::dgdr(1, Gamma, beta, M2, gammaAdi, dM2dR);
            break;
    }
//        double dmomdR = dGammadR * (Gamma / std::sqrt(Gamma * Gamma - 1.0));
//        double dmomdR = dGammadR / betaSh;


//        if (dGammadR > 0){
//            if ( Gamma > 0.95 * p_pars->Gamma0 ){
//                dGammadR = 0.;
//            }
//        }
    double dtcomov_dR = 1.0 / EQS::MomFromGamma(Gamma) / CGS::c;


    // -- Energies --
    double dEsh2dR  = (Gamma - 1.0) * dM2dR; // Shocked energy;
    double dlnV2dR  = dM2dR / M2 - drhodr / rho - dGammadR / Gamma;
    double dEad2dR  = 0.0;
    if ( p_pars->adiabLoss )
        dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;

    // -- Radiative losses
    double eps_rad = p_pars->eps_rad;
    if (p_pars->eps_rad < 0){
        if (Eint2 < 0 || M2 < 0 || Gamma < 1){
            (*p_log)(LOG_ERR,AT) << " Eint="<<Eint2<<" M2="<<M2<<" Gamma="<<Gamma<<"\n";
            exit(1);
        }
        eps_rad = p_pars->p_mphys->computeRadiationLoss(
                dM2dR*p_pars->M0, dEsh2dR*p_pars->M0*CGS::c*CGS::c,
                R, Gamma, GammaSh, gammaAdi, M2*p_pars->M0, Eint2*p_pars->M0*CGS::c*CGS::c,
                rho*p_pars->M0, theta, p_pars->ncells
                );
//        /// compute radiative losses by integrating synchrotron spectrum (analytic synchrotron spectrum)
//        double rho2 = EQS::rho2t(Gamma,gammaAdi,rho)*p_pars->M0;
//        double m = M2*p_pars->M0;
//        double U_p = p_pars->p_mphys->get_shock_Up(GammaSh,rho2,m,
//                                                   Eint2*CGS::c*CGS::c*p_pars->M0);
//        p_pars->p_mphys->updateSockProperties(U_p,Gamma,GammaSh,
//                                              x,tcomov,rho2/CGS::mp,
//                                              m/CGS::mp);
//        thickness = EQS::shock_thickness(m,rho2,
//                                         0.,theta,Gamma,R,p_pars->ncells);
//        double dm = dM2dR*p_pars->M0;
//        p_tot = p_pars->p_mphys->computePtot(dm);
//        p_tot = std::max(p_tot,0.);
//        /// compute total amount of energy radaited
//        double drcomov = thickness * GammaSh;
//        double dt_esc = drcomov / CGS::c;
//
////        dErad2dR = p_tot / beta / Gamma;
//        dErad2dR = p_tot * dt_esc;
//        dErad2dR = dErad2dR / (CGS::c*CGS::c*p_pars->M0);
//        eps_rad = dErad2dR / (dEsh2dR*p_pars->p_mphys->eps_e);
    }
    double dErad2dR = eps_rad * p_pars->p_mphys->eps_e * dEsh2dR;


    // -- Energy equation
    double dEint2dR = dEsh2dR + dEad2dR - dErad2dR; // / (m_pars.M0 * c ** 2)
#if 0
    double dttdr = 0.;
        if (p_spread->m_method != LatSpread::METHODS_SYNCH::iNULL)
            dttdr = 1. / (CGS::c * betaSh) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
        else
            dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
        dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
#endif
    bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
    double dttdr = EQS::evalElapsedTime(R, EQS::MomFromGamma(Gamma),dthetadr,spread);
    // ****************************************

    if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) ||
        !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr) ||
        !std::isfinite(dEint2dR) || !std::isfinite(dthetadr)) {
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives. Exiting..." << "\n";
        exit(1);
    }
    // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
    if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0. || dthetadr < 0
        || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
        (*p_log)(LOG_ERR,AT) << " nan in derivatives. "
                             << " dRdt="<<dRdt<<"\n"
                             << " dM2dR="<<dM2dR<<"\n"
                             << " dthetadr=" << dthetadr << "\n"
                             << " dGammadR=" << dGammadR << "\n"
                             << " dthetadr=" << dRdt << "\n"
                             << " dEsh2dR=" << dEsh2dR << "\n"
                             << " dlnV2dR=" << dlnV2dR << "\n"
                             << " dEad2dR=" << dEad2dR << "\n"
                             << " dErad2dR=" << dErad2dR << "\n"
                             << " dEint2dR=" << dEint2dR << "\n"
                             << " dttdr=" << dttdr << "\n"
                             << "  \n";
//            std::cerr << AT << "\n";
        exit(1);
    }
    if (std::abs(dRdt * dGammadR) > p_pars->Gamma0){
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives. "
                              << "i="<< i << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"] " << "\n"
                              << " dGamma/dr > Gamma0 "
                              <<  "dG/dr="<< dRdt * dGammadR
                              << " Gamma0=" <<p_pars->Gamma0 << "\n"
                              << " E0=" <<p_pars->E0 << "\n"
                              << " M0=" <<p_pars->M0 << "\n"
                              << " R=" <<R << " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
                              << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
        (*p_log)(LOG_ERR,AT)  << " maybe timestep is too large. Exiting...";
//            std::cerr << AT << "\n";
        exit(1);
    }

    out_Y[i + SOL::QS::iR] = dRdt;//1.0 / betaSh / CGS::c;
    out_Y[i + SOL::QS::iRsh] = dRshdt;//1.0 / betaSh / CGS::c;
    out_Y[i + SOL::QS::itt] = dRdt * dttdr;
    out_Y[i + SOL::QS::itcomov] = dRdt * dtcomov_dR;
    out_Y[i + SOL::QS::iGamma] = dRdt * dGammadR;//dRdt * dGammadR/Gamma;
//        out_Y[i + SOL::QS::imom] = dRdt * dmomdR;//dRdt * dGammadR/Gamma;
    out_Y[i + SOL::QS::iEint2] = dRdt * dEint2dR;
    out_Y[i + SOL::QS::itheta] = dRdt * dthetadr;
    out_Y[i + SOL::QS::iErad2] = dRdt * dErad2dR;
    out_Y[i + SOL::QS::iEsh2] = dRdt * dEsh2dR;
    out_Y[i + SOL::QS::iEad2] = dRdt * dEad2dR;
    out_Y[i + SOL::QS::iM2] = dRdt * dM2dR;
//        if (betaSh>1e-5)  p_pars->end_evolution = true;
    // ****************************************
//        if (mom+out_Y[i + QS::imom])
}

/// right hand side of the blast wave evolution equation (ODE) with FS and RS
void BlastWave::rhs_fsrs(double * out_Y, size_t i, double x, double const * Y ) {
    // ****************************************
    double R      = Y[i+ SOL::QS::iR];
    double Rsh    = Y[i+ SOL::QS::iRsh];
    double tcomov = Y[i+ SOL::QS::itcomov];
//        double mom    = Y[i+ SOL::QS::imom];
    double Gamma  = Y[i+ SOL::QS::iGamma];
    double Eint2  = Y[i+ SOL::QS::iEint2];
    double Eint3  = Y[i+ SOL::QS::iEint3];
    double Ead3   = Y[i+ SOL::QS::iEad3];
    double Erad3  = Y[i+ SOL::QS::iErad3];
    double theta  = Y[i+ SOL::QS::itheta];
    double M2     = Y[i+ SOL::QS::iM2];
    double M3     = Y[i+ SOL::QS::iM3];
    double deltaR4= Y[i+ SOL::QS::ideltaR4];
    // ****************************************
    if (Gamma < 1.) {
        (*p_log)(LOG_ERR, AT) << " Gamma < 1 Gamma="<<Gamma<<" -> Gamma0="<<p_pars->Gamma0<<"\n";
        if (p_pars->prev_idx_x == 0)
            Gamma = p_pars->Gamma0;
        else
            Gamma = 1.+1.e-5; //
//            exit(1);
    }
//    if (deltaR4 < 0){
//        (*p_log)(LOG_ERR, AT) << " deltaR4 < 1 deltaR4="<<deltaR4<<" -> 0\n";
//        if (p_pars->prev_idx_x == 0)
//            deltaR4 = 0;
//    }
//        if(Gamma>p_pars->Gamma0){
//            (*p_log)(LOG_ERR, AT) << " Gamma > Gamma0; Gamma="<<Gamma
//                <<" Gamma0="<<p_pars->Gamma0<<" dGammadt_last="<<p_pars->prev_dGammadR<<"\n";
//            exit(1);
//        }
    if(Eint3 < 0){
        std::cerr << " p_pars->prev_dM3dR = "<<p_pars->prev_dM3dR<<" Eint3="<<" "<<Eint3<<" -> 0\n";
        Eint3 = 0.0;
    }
    if(M3 < 0){
        std::cerr << " p_pars->prev_dM3dR = "<<p_pars->prev_dM3dR<<" M3="<<" "<<p_pars->prev_idx_x<<" -> 0\n";
        if (p_pars->prev_idx_x == 0) {
            M3 = 0;
            Eint3 = 0;
            deltaR4 = 0;
        }
//            M3 = 0.0;
    }
//        if(Eint2 < 0){
//            Eint2 = 0.0;
//        }
//        double Gamma = EQS::GamFromMom(mom);
//        if (Gamma <= 1.)
//            Gamma = EQS::GamFromMom(mom);
//            Gamma = 1.0000000001;
//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }
//        double betaSh   = BetFromMom(mom);
    double beta   = EQS::BetaFromGamma(Gamma);

    long double Gamma43 = EQS::get_Gamma43(Gamma, p_pars->Gamma0, beta, p_pars->beta0);
    if (Gamma43 < 1.l)
        Gamma43 = 1.l;
    long double beta43 = EQS::BetaFromGamma(Gamma43);
//    if (Gamma43*beta43 < 2e-3)
//        int z = 1;

    double gammaAdi3 = p_eos->getGammaAdi((double)Gamma43, beta43);

    double GammaRsh = EQS::GammaSh((double)Gamma43,gammaAdi3);
    double dRshRshdt = EQS::dRsh_dt(Gamma, beta, GammaRsh);

    // ****************************************
//        auto *_pars = (struct RHS_pars *) rhs_pars;
//    if (dRshRshdt > 0.){
//        int z = 1;
//    }
    double ctheta_ = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);//p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
    p_dens->evaluateRhoDrhoDrDefault(R, ctheta_);
    double rho = p_dens->m_rho_def / p_pars->M0;
    double drhodr = p_dens->m_drhodr_def / p_pars->M0;

//        dlnrho1dR /= p_pars->M0;
//        double _rhoi = p_pars->rho/p_pars->M0;
    // ****************************************
    double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...


    double gammaAdi  = p_eos->getGammaAdi(Gamma, beta);
    double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
//    double dRshdt = EQS::Beta(GammaSh) * CGS::c;
    double dRshdt = EQS::dRsh_dt(Gamma, beta, GammaSh);
//    std::cout << " drdt = " << dRdt << " drsh_dt = " << dRshdt << "\n";

    double dthetadr = 0.0;
    if ((theta < p_pars->theta_max) && (!p_pars->end_spreading)) {
        switch (p_pars->method_limit_spread) {

            case iNone:
                dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iGamma0Frac:
                if (Gamma*beta < p_pars->fraction_of_mom0_when_spread_start*p_pars->mom0)
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iGammaVal:
                if (Gamma*beta < std::max(p_pars->value_of_mom_when_spread_start,
                                          p_pars->fraction_of_mom0_when_spread_start*p_pars->mom0))
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
            case iRd:
                if (R > p_pars->Rd)
                    dthetadr = p_spread->getDthetaDr(Gamma, GammaSh, R, gammaAdi, theta);
                break;
        }
    }


    double dM2dR = 0;
    switch (p_pars->m_method_dmdr) {

        case iusingA:
            dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
            break;
        case iusingdthdR:
            dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
            break;
        case iNodmdr:
            break;
    }
    if (dM2dR < 0.){
        (*p_log)(LOG_ERR, AT) << " dM/dR < 0\n";
        exit(1);
    }

    // Reverse Shock [DENSITY] Compute the evolution of the Ejecta density profile
    double ddeltaR4dR = 0.;
    long double dM3dR = 0.;
    double rho4 = 0.;
    double dlnrho4dR = 0.;
    // (Gamma < 0.9999*p_pars->Gamma0) (x > p_pars->tb0*1.e-5)
    if ((!p_pars->shutOff) && (Gamma < 0.9999 * p_pars->Gamma0) ){ // and (not M3 > 1.):# and (deltaR4 > 0): # the last assures that jet is decelerating
        double alpha_of = p_pars->tprompt * p_pars->beta0 * CGS::c;// * EjectaID2::CellsInLayer(p_pars->ilayer);
        ddeltaR4dR = (1.0 / (beta*beta*beta*beta) - 1.0 / (p_pars->beta0*p_pars->beta0*p_pars->beta0*p_pars->beta0)) /
                     (1.0 / beta + 1.0 / p_pars->beta0)
                     / (1.0 / (beta*beta) + 1.0 / (p_pars->beta0*p_pars->beta0)) * p_pars->beta0;
        long double rho4_scale_factor = -deltaR4 / alpha_of;
        // rho4 = rho4_factor * np.exp(rho4_scale_factor)  # Here exp() does not really make a difference
        if (p_pars->exponential_rho4) {
            dM3dR = (1.0 / alpha_of) * ddeltaR4dR * std::exp(rho4_scale_factor);
            // rho3prim = 4 * Gamma * rho4
            dlnrho4dR = -2.0 / R - ddeltaR4dR / alpha_of;
        }
        else{
            dM3dR = (1.0 / alpha_of) * ddeltaR4dR;
            dlnrho4dR = -2.0 / R;
        }

//        double bb0m1 = 0.5 / (Gamma*Gamma)
//                + 0.5 / (p_pars->Gamma0*p_pars->Gamma0)
//                + (1./(Gamma*Gamma*Gamma*Gamma)) / 8
//                + (1/(p_pars->Gamma0*p_pars->Gamma0*p_pars->Gamma0*p_pars->Gamma0)) / 8
//                - 0.25 / (p_pars->Gamma0*p_pars->Gamma0) / (Gamma*Gamma);
//        dM3dR = 1 * (p_pars->beta0 *p_pars->beta0 - beta *beta)
//                 / (p_pars->beta0 + beta) / p_pars->beta0 / p_pars->tprompt / beta / c / bb0m1 ;
//        dlnrho4dR = -2 / R;
    }

    // --- dGammadR ---
    double dlnrhodr = drhodr / rho;
    double dGammadR = EQS::get_dGammadR_fs_rs(
            Gamma, p_pars->Gamma0, gammaAdi, dlnrhodr, M2, dM2dR,
            dlnrho4dR, M3, dM3dR, Eint2, Eint3, gammaAdi3, Gamma43);
    if(dGammadR > 0){
        if (Gamma > p_pars->rs_Gamma0_frac_no_exceed*p_pars->Gamma0) {
            dGammadR = 0.0;
        }
    }
//        if (dRdt * dGammadR > 0){
//            std::cerr << AT <<"\n [ Error ] dRdt * dGammadR > 0"
//                              " Eint2="<<Eint2<<" Gamma="<<Gamma<<" betaSh="<<EQS::Beta(Gamma)<<" Eint2="<<Eint2
//                      <<" Eint3"<<Eint3<<" M2"<<M2<<" M3"<<M3<<"Exiting...\n";
//            exit(1);
//        }
//        if (!std::isfinite(dGammadR)){
//            std::cerr << AT <<"\n [ Error ] dGammadR = nan"
//                              " Eint2="<<Eint2<<" Gamma="<<Gamma<<" betaSh="<<EQS::Beta(Gamma)<<" Eint2="<<Eint2
//                      <<" Eint3"<<Eint3<<" M2"<<M2<<" M3"<<M3<<" Exiting...\n";
//            exit(1);
//        }

////        double dmomdR = dGammadR * (Gamma / std::sqrt(Gamma * Gamma - 1.0));
//        double dmomdR = dGammadR / betaSh;


    // -- Energies --
    double dEsh2dR  = (Gamma - 1.0) * dM2dR; // Shocked energy;

    long double dEsh3dR = (Gamma43 - 1.0) * dM3dR; // Here it should be Gamma43 I presume...

    double dEad2dR  = 0.0;
    double dlnV2dR  = dM2dR / M2 - drhodr / rho - dGammadR / Gamma;
    if ( p_pars->adiabLoss )
        dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;

    double dEad3dR = 0.;
    if((!p_pars->shutOff) && (M3 > 0.0)){
        // Shocked energy
//        dEsh3dR = gamma43_minus_one * (double)dM3dR;
        // Expansion energy
        double dlnV3dR = (double)dM3dR / M3 - dlnrho4dR - dGammadR / Gamma;
        dEad3dR = 0.;
        if(p_pars->adiabLoss_rs){
            dEad3dR = -(gammaAdi3 - 1) * Eint3 * dlnV3dR;
        }
    }

    // -- Radiative losses
    double eps_rad = p_pars->eps_rad;
    if (p_pars->eps_rad < 0)
        eps_rad = p_pars->p_mphys->computeRadiationLoss(
                dM2dR*p_pars->M0, dEsh2dR*p_pars->M0*CGS::c*CGS::c,
                R, Gamma, GammaSh, gammaAdi, M2*p_pars->M0, Eint2*p_pars->M0*CGS::c*CGS::c,
                rho*p_pars->M0, theta, p_pars->ncells
        );
    double dErad2dR = eps_rad * p_pars->p_mphys->eps_e * dEsh2dR;

    double eps_rad_rs = p_pars->epsilon_rad_rs;
    if (p_pars->epsilon_rad_rs < 0)
        eps_rad_rs = p_pars->p_mphys_rs->computeRadiationLoss(
                (double)dM3dR*p_pars->M0, (double)dEsh3dR*p_pars->M0*CGS::c*CGS::c,
                R, Gamma, GammaRsh, gammaAdi3, M3*p_pars->M0, Eint3*p_pars->M0*CGS::c*CGS::c,
                rho4*p_pars->M0, theta, p_pars->ncells
        );
    double dErad3dR = eps_rad_rs * p_pars->p_mphys_rs->eps_e * (double)dEsh3dR;

    // -- Energy equation
    double dEint2dR = dEsh2dR + dEad2dR - dErad2dR; // / (m_pars.M0 * c ** 2)
    double dEint3dR = dEsh3dR + dEad3dR - dErad3dR;  // / (m_pars.M0 * c ** 2)
    double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
#if 0
    double dttdr = 0.;
        if (p_spread->m_method != LatSpread::METHODS_SYNCH::iNULL)
            dttdr = 1. / (CGS::c * betaSh) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
        else
            dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
        dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
#endif
    bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
    double dttdr = EQS::evalElapsedTime(R,Gamma*beta,dthetadr,spread);
    // ****************************************

    // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
    if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0.
        || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
        (*p_log)(LOG_ERR,AT) << " nan in derivatives. Try lowering 'rs_Gamma0_frac_no_exceed' or increasing 'ode_rtol' "
                             << " rs_Gamma0_frac_no_exceed="<<p_pars->rs_Gamma0_frac_no_exceed
                             << " dRdt="<<dRdt<<"\n"
                             << " dM2dR="<<dM2dR<<"\n"
                             << " dthetadr=" << dthetadr << "\n"
                             << " dGammadR=" << dGammadR << "\n"
                             << " dthetadr=" << dRdt << "\n"
                             << " dEsh2dR=" << dEsh2dR << "\n"
                             << " dlnV2dR=" << dlnV2dR << "\n"
                             << " dEad2dR=" << dEad2dR << "\n"
                             << " dErad2dR=" << dErad2dR << "\n"
                             << " dEint2dR=" << dEint2dR << "\n"
                             << " dttdr=" << dttdr << "\n"
                             << "  \n";
//            std::cerr << AT << "\n";
        exit(1);
    }
//        if (!std::isfinite(dM3dR) || !std::isfinite(dEad3dR) || dM3dR < 0.
//            || !std::isfinite(dlnrho4dR) || !std::isfinite(dEint3dR)) {
//            (*p_log)(LOG_ERR,AT) << " nan in derivatives. "
//                                 << " dRdt="<<dRdt<<"\n"
//                                 << " dM2dR="<<dM2dR<<"\n"
//                                 << " dM3dR="<<dM3dR<<"\n"
//                                 << " dthetadr=" << dthetadr << "\n"
//                                 << " dGammadR=" << dGammadR << "\n"
//                                 << " dthetadr=" << dRdt << "\n"
//                                 << " dEsh2dR=" << dEsh2dR << "\n"
//                                 << " dEsh3dR=" << dEsh3dR << "\n"
//                                 << " dlnV2dR=" << dlnV2dR << "\n"
//                                 << " dEad2dR=" << dEad2dR << "\n"
//                                 << " dEad3dR=" << dEad3dR << "\n"
//                                 << " dErad2dR=" << dErad2dR << "\n"
//                                 << " dErad3dR=" << dErad3dR << "\n"
//                                 << " dEint2dR=" << dEint2dR << "\n"
//                                 << " dEint3dR=" << dEint3dR << "\n"
//                                 << " ddeltaR4dR=" << ddeltaR4dR << "\n"
//                                 << " dttdr=" << dttdr << "\n"
//                                 << "  \n";
////            std::cerr << AT << "\n";
//            exit(1);
//        }
//        if (std::abs(dRdt * dGammadR) > p_pars->Gamma0){
//            (*p_log)(LOG_ERR,AT)  << " nan in derivatives. "
//                                  << "i="<< i << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"] "
//                                  << " dGamma/dr > Gamma0 "
//                                  <<  "dG/dr="<< dRdt * dGammadR
//                                  << " Gamma0=" <<p_pars->Gamma0
//                                  << " E0=" <<p_pars->E0
//                                  << " M0=" <<p_pars->M0
//                                  << " R=" <<R << " Rsh=" <<Rsh<< " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
//                                  << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
//            (*p_log)(LOG_ERR,AT)  << " maybe timestep is too large. Exiting...";
////            std::cerr << AT << "\n";
//            exit(1);
//        }

    out_Y[i + SOL::QS::iR] = dRdt;//1.0 / betaSh / CGS::c;
    out_Y[i + SOL::QS::iRsh] = dRshdt;//1.0 / betaSh / CGS::c;
    out_Y[i + SOL::QS::itt] = dRdt * dttdr;
    out_Y[i + SOL::QS::itcomov] = dRdt * dtcomov_dR;
    out_Y[i + SOL::QS::iGamma] = dRdt * dGammadR;//dRdt * dGammadR/Gamma;
//        out_Y[i + SOL::QS::imom] = dRdt * dmomdR;//dRdt * dGammadR/Gamma;
    out_Y[i + SOL::QS::iRrsh] = dRshRshdt;//1.0 / betaSh / CGS::c;
    out_Y[i + SOL::QS::iEint2] = dRdt * dEint2dR;
    out_Y[i + SOL::QS::iEint3] = dRdt * dEint3dR;
    out_Y[i + SOL::QS::itheta] = dRdt * dthetadr;
    out_Y[i + SOL::QS::iErad2] = dRdt * dErad2dR;
    out_Y[i + SOL::QS::iErad3] = dRdt * dErad3dR;
    out_Y[i + SOL::QS::iEsh2] = dRdt * dEsh2dR;
    out_Y[i + SOL::QS::iEsh3] = dRdt * dEsh3dR;
    out_Y[i + SOL::QS::iEad2] = dRdt * dEad2dR;
    out_Y[i + SOL::QS::iEad3] = dRdt * dEad3dR;
    out_Y[i + SOL::QS::iM2] = dRdt * dM2dR;
    out_Y[i + SOL::QS::iM3] = dRdt * dM3dR;
    out_Y[i + SOL::QS::ideltaR4] = dRdt * ddeltaR4dR;
//        if (betaSh>1e-5)  p_pars->end_evolution = true;
    p_pars->prev_dM3dR = dRdt * dM3dR;
    p_pars->prev_dGammadR = dRdt * dGammadR;
    p_pars->prev_dRdt = dRdt;


    // ****************************************
//        if (mom+out_Y[i + QS::imom])
}


/// RHS with ODEs for dGammadR modified for pre-accelerated ISM
void BlastWave::rhs_fs_dense(double * out_Y, size_t i, double x, double const * Y ) {
    double Gamma = Y[i+SOL::QS::iGamma];//EQS::GamFromMom(Y[i+QS::imom]);
//        double mom = Y[i+ SOL::QS::imom];
    // ****************************************
    double R      = Y[i+ SOL::QS::iR];
//    double Rsh    = Y[i+ SOL::QS::iRsh];
//        double tcomov = Y[i+Q::itcomov];
//        double Gamma  = std::exp(Y[i+QS::ilnGamma]);
    double Eint2  = Y[i+ SOL::QS::iEint2];
    double theta  = Y[i+ SOL::QS::itheta];
    double M2     = Y[i+ SOL::QS::iM2];
    // *****************************************
//        double Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
    if (Gamma < 1.){
        (*p_log)(LOG_ERR,AT) << "Error Gamma < 1.\n";
        Gamma = 1. + 1.e-5;
//            Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
//            betaSh = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
    }
    double beta = EQS::BetaFromGamma(Y[i+ SOL::QS::iGamma]);//EQS::BetFromMom(Y[i+QS::imom]);
//        double betaSh = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
    // ****************************************
    if (!std::isfinite(R) || !std::isfinite(Gamma) || M2 < 0.
        || !std::isfinite(M2) || !std::isfinite(Eint2) || Eint2<0) {
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " << "\n"
                              << " R="<<R<<"\n"
                              << " M2="<<M2<<"\n"
                              << " Gamma=" << Gamma << "\n"
                              //                                  << " Mom=" << mom << "\n"
                              << " Eint2=" << Eint2
                              << " \n";
        exit(1);
    }
//        if (Gamma <= 1.) { // TODO to be removed
//            Gamma = 1.0001;
//            (*p_log)(LOG_ERR, AT) << " Gamma < 1 in RHS for kN Ejecta\n";
//        }
//        if (Gamma > p_pars->Gamma0){ // TODO to be removed
//            Gamma = p_pars->Gamma0;
//            (*p_log)(LOG_ERR, AT) << " Gamma > Gamma0 in RHS for kN Ejecta\n";
//        }
//        if (mom < 0.){
//            (*p_log)(LOG_ERR, AT) << " mom < 0 = "<< mom << " in kN RHS dynamics\n";
//            exit(1);
//        }



//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }
//        double betaSh   = EQS::Beta(Gamma);


    // ****************************************
    // Get ISM density and its velocity
//        double rho, m_drhodr;
    double ctheta_ = EjectaID2::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);// = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
//        p_dens->getDrhoDr( rho, m_drhodr, R,ctheta );
//        rho /= p_pars->M0;
//        m_drhodr /= p_pars->M0;
//        p_dens->evaluate(R, ctheta);
    if (p_dens->m_GammaRel < 0){
        (*p_log)(LOG_ERR,AT) << AT << " GammaRel is not set or incorrect \n";
        exit(1);
    }
    double rho = p_dens->m_rho_ / p_pars->M0;
    double drhodr = p_dens->m_drhodr_ / p_pars->M0;
    double GammaRho = p_dens->m_GammaRho;
    double GammaRel = p_dens->m_GammaRel; // Gamma
    double dGammaRelDGamma = p_dens->m_dGammaReldGamma; // 1.
    double dGammaRhodR = p_dens->m_dGammaRhodR; // 0.
    double cs_cbm = p_dens->m_CS_CBM;


    if (GammaRho < 1.) {
        (*p_log)(LOG_ERR,AT)  << "GammaRho=" << GammaRho << "\n" << " Exiting...";
//            std::cerr << AT<< "\n";
        exit(1);
    }
    if (GammaRel < 1.) {
        (*p_log)(LOG_ERR,AT)  << "GammaRel=" << GammaRel << "\n"; GammaRel = Gamma;
//            std::cerr << AT<< "\n";
        exit(1);
    }
    if (!std::isfinite(dGammaRelDGamma)) { (*p_log)(LOG_ERR,AT) << "dGammaRelDGamma=" << dGammaRelDGamma << "\n"; exit(1); }
    if (!std::isfinite(dGammaRhodR)) { (*p_log)(LOG_ERR,AT)  << "dlnGammaCBMdR=" << dGammaRhodR << "\n"; exit(1); }
    // ****************************************
    double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...
    double Gamma_ = GammaRel;
    double beta_ = EQS::BetaFromGamma(Gamma_);
    double gammaAdi  = p_eos->getGammaAdi(Gamma_, beta_);//p_eos->getGammaAdi(Gamma, betaSh);

    double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
    double dRshdt = EQS::BetaFromGamma(GammaSh) * CGS::c;

    double dthetadr = 0.0;
//        if (theta < (p_pars->theta_max * 0.999999) ) {
//            dthetadr = p_spread->getDthetaDr(Gamma, R, gammaAdi, theta);
//        }
//        double dM2dR     = EQS::dmdr(Gamma, R, p_pars->theta_a,
//                                     theta, rho, p_spread->m_aa) / p_pars->ncells;

    double dM2dR = 0.;

    switch (p_pars->m_method_dmdr) {

        case iusingA:
            /// Old equation from Gavin P. Lamb work
            dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
            break;
        case iusingdthdR:
            /// Equation motivated by Gavin P. Lamb paper
            dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
            break;
        case iNodmdr:
            /// if no mass needs to be swept up
            break;
    }
//        dM2dR*=0;
    if (dM2dR < 0.){
        (*p_log)(LOG_ERR,AT) << " dMdR < 0 in RHS dyn for kN ejecta\n";
        exit(1);
    }

    // --- Energy injection --- ||
    double xi_inj = 1.; // TODO put in parfile
    double dEindt = p_pars->dEinjdt;
    if (!std::isfinite( dEindt) || dEindt < 0){
        (*p_log)(LOG_ERR,AT) << " dEindt = "<<dEindt << "\n";
        exit(1);
    }

    double dEinjdt = dEindt / (p_pars->M0 * CGS::c * CGS::c) / p_pars->ncells;
    double dEinjdR = dEinjdt / dRdt;
    double theta_ej = 0.; // assume that ejecta is alinged with magnetar emission?..
    double Doppler = Gamma / (1. - beta * std::cos(theta_ej));
    double dEinjdR_dop = dEinjdR * Doppler;
    double dEingdR_abs = dEinjdR;// * ( 1. - std::exp(-1.*p_pars->dtau) ) * std::exp(-1.*p_pars->tau_to0);
    double dEingdR_abs_dop = dEingdR_abs / Doppler / Doppler;

    double dEnuc = p_pars->dEnuc;//->getPars()->eps_nuc_thermalized;
    double dEnuc_norm = dEnuc / (p_pars->M0 * CGS::c * CGS::c);
    double dEnucdR = dEnuc_norm / dRdt;

    double dElum = p_pars->dElum;
    double dElum_norm = dElum / (p_pars->M0 * CGS::c * CGS::c);
    double dElumdR = dElum_norm / dRdt;
    if ((std::abs(dElumdR) > std::abs(dEnucdR*1.5)) and (dEnuc > 0.))
        dElumdR = dEnucdR;

    // --- dGammadR ---
    double dGammadR = 0., GammaEff=0.,dGammaEffdGamma=0.,num1=0.,num2=0.,num3=0.,denum1=0.,denum2=0.,denom3=0.;
//        double _tmp = (1. - GammaEff / Gamma * xi_inj) * dEingdR_abs;
    switch (p_pars->m_method_dgdr) {

        case iour:
            GammaEff = get_GammaEff(Gamma_, gammaAdi); // TODO check if GammaRel should be used for this!!!
            dGammaEffdGamma = get_dGammaEffdGamma(Gamma_, gammaAdi);
            num1 = (Gamma - GammaRho + GammaEff * (GammaRel - 1.)) * dM2dR;
            num2 = - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR/M2 - drhodr/rho - dGammaRhodR / GammaRho); // - 3.*Eint2/R
            num3 = - (1. + GammaEff / Gamma * xi_inj) * dEingdR_abs;
            denum1 = (1.+M2);
            denum2 = Eint2 * dGammaEffdGamma;
            denom3 = GammaEff * (gammaAdi - 1.) * Eint2 * dGammaRelDGamma / GammaRel;
            dGammadR = -1. * (num1 + num2 + num3) / (denum1+denum2+denom3);
            break;
        case ipeer:
            dGammadR = EQS::dgdr(1, Gamma, beta, M2, gammaAdi, dM2dR);
            break;
    }



//        if ((dGammadR > 0) && (Gamma > p_pars->Gamma0*0.99)){
//            std::cerr << AT << " \n"
//                      << "dGammadR>0 ("<<dGammadR<<") and Gamma="<<Gamma<<" > Gamma0="<<p_pars->Gamma0<<"\n";
//            std::cerr << " p_dens->m_GammaRho="<<p_dens->m_GammaRho<<"\n"
//                      << " p_dens->m_GammaRel="<<p_dens->m_GammaRel<<"\n"
//                      << " p_dens->m_dGammaReldGamma="<<p_dens->m_dGammaReldGamma<<"\n"
//                      << " p_dens->m_dGammaRhodR="<<p_dens->m_dGammaRhodR<<"\n";
//            dGammadR = 0.;
//        }
    //        double dmomdR = dGammadR * (Gamma / std::sqrt(Gamma * Gamma - 1.0));
//        double dmomdR = dGammadR / betaSh;


//        drhodr = 0.; // TODO why was this set?! This overrides the derivative
    // -- Energies --

//        double dlnV2dR  = dM2dR / M2 - m_drhodr - dGammadR / Gamma;
    double dlnV2dR  = dM2dR / M2 - drhodr / rho - (1./GammaRel)*dGammaRelDGamma*dGammadR + dGammaRhodR / GammaRho;// + 3./R;
    double dlnVdR = (3./R) - (1/Gamma * dGammadR);
    double dEad2dR  = 0.0;
    if ( p_pars->adiabLoss )
        dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;
//        double dx = dRdt * (x - p_pars->prev_x);

//        double mom_ = mom + dmomdR * dRdt * (x - p_pars->x);
//        if (mom_ < 0){
//            int x = 1;
//        }

    // --- shock energy
    double dEsh2dR = (GammaRel - 1.0) * dM2dR;
//        if (cs_cbm > 100){
//            int y = 1;
//        }
    if ((cs_cbm > EQS::BetaFromGamma(GammaRel))){//&&(p_pars->comp_ix > p_pars->nr*1e-2))){ // TODO Subsonic flow -- no shock
        dEsh2dR *= 1e-10;
    }
//        double dEsh2dR = 0.;

//        double dEsh2dR  = (GammaRel - 1.0) * dM2dR; // Shocked energy;
    // -- Radiative losses
    double dErad2dR = p_pars->eps_rad * dEsh2dR;
    // -- Energy equation
    double dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEnucdR - dElumdR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)
//        double _x = dEint2dR/Eint2;
//        double Eint_ = Eint2 + dEint2dR * dRdt * (x-p_pars->prev_x);
//        if (Eint_ < 0){
//            dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)
//            double Eint__ = Eint2 + dEint2dR * dRdt * (x-p_pars->prev_x);
//            if (Eint__ < 0){
//                int z = 1;
//            }
//        }

    double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
//        double dttdr;
//        if (p_spread->m_method != LatSpread::METHODS_SYNCH::iNULL)
//            dttdr = 1. / (CGS::c * betaSh) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
//        else
//            dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
//        dttdr = 1. / (CGS::c * Gamma * Gamma * betaSh * (1. + betaSh));
    bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
    double dttdr = EQS::evalElapsedTime(R,Gamma*beta,dthetadr,spread);
    // ****************************************
//        if (!std::isfinite(Gamma) ||
//        !std::isfinite(betaSh) //||
////        !std::isfinite(Eint2) ||
////        !std::isfinite(theta) ||
////        !std::isfinite(M2)
////        || Eint2 < 0.
////        || M2 < 0.
//          ) {
//            std::cerr << AT << " wrong data at the beginning of RHS t_b=" << x << "\n";
//            std::cerr << " Gamma="<<Gamma<<"\n"
//                      << " theta="<<theta<<"\n"
//                      << " M2="<<M2<<"\n"
//                      << " Eint2="<<Eint2<<"\n"
//                      << " p_dens->m_GammaRho="<<p_dens->m_GammaRho<<"\n"
//                      << " p_dens->Gamma_rel="<<p_dens->m_GammaRel<<"\n"
//                      << " p_dens->m_dGammaRhodR="<<p_dens->m_dGammaRhodR<<"\n"
//                      << " p_dens->m_dGammaReldGamma="<<p_dens->m_dGammaReldGamma<<"\n"
//                      << " p_dens->m_CS_CBM="<<p_dens->m_CS_CBM<<"\n"
//                      << " p_dens->m_P_cbm="<<p_dens->m_P_cbm<<"\n";
//            exit(1);
//        }
    if (Gamma < 1.){
        (*p_log)(LOG_ERR, AT) << " Gamma < 1 = "<< Gamma << " in kN RHS dynamics\n";
//            exit(1);
    }
    if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0.
        || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " //<< "\n"
                              << " dRdt="<<dRdt//<<"\n"
                              << " dM2dR="<<dM2dR//<<"\n"
                              << " dthetadr=" << dthetadr// << "\n"
                              << " dGammadR=" << dGammadR //<< "\n"
                              << " dthetadr=" << dRdt //<< "\n"
                              << " dEsh2dR=" << dEsh2dR //<< "\n"
                              << " dlnV2dR=" << dlnV2dR //<< "\n"
                              << " dEad2dR=" << dEad2dR //<< "\n"
                              << " dErad2dR=" << dErad2dR //<< "\n"
                              << " dEint2dR=" << dEint2dR //<< "\n"
                              << " dttdr=" << dttdr// << "\n"
                              << " dGammaRelDGamma=" << dGammaRelDGamma //<< "\n"
                              << " dGammaRhodR=" << dGammaRhodR //<< "\n"
                              << " drhodr=" << drhodr //<< "\n"
                              << " Gamma="<<Gamma//<< "\n"
                              << " betaSh="<<beta//<< "\n"
                              << " theta"<<theta//<< "\n"
                              << " M2="<<M2//<< "\n"
                              << " Eint2="<<Eint2//<< "\n"
                              << " GammaRel"<<GammaRel//<< "\n"
                              << " GammaRho"<<GammaRho//<< "\n"
                              << " \n";
        exit(1);
    }
    // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
//        if (dGammadR > 0){
//            std::cout << " Gamma="<<Gamma<< " Gamma0="<<p_pars->Gamma0<<" Gamma/Gamma0="<<Gamma/p_pars->Gamma0<<"\n";
//        }
//        p_pars->x = x; // update
//        if (Gamma > 2.*p_pars->Gamma0){
//            int x = 1;
//        }
//        if ( std::abs(dRdt * dGammadR ) > 2. * p_pars->Gamma0 )+{
//            dGammadR = 0.;
//        }
//        if ( std::abs(dRdt * dGammadR ) > p_pars->Gamma0 ){
//            std::cerr << AT
//                      << "i="<< i << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"] " << "\n"
//                      << " dGamma/dr > Gamma0 "
//                      <<  "dG/dr="<< dRdt * dGammadR
//                      << " Gamma0=" <<p_pars->Gamma0 << "\n"
//                      << " theta_b0=" <<p_pars->theta_b0 << "\n"
//                      << " GammaRho=" <<GammaRho << "\n"
//                      << " GammaRel=" <<GammaRel << "\n"
//                      << " dGammaRelDGamma=" <<dGammaRelDGamma << "\n"
//                      << " dGammaRhodR=" <<dGammaRhodR << "\n"
//                      << " R=" <<R << " Rsh=" <<Rsh<< " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
//                      << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
//            std::cerr << " maybe timestep is too large \n";
////            exit(1);
//            dGammadR = 0.;
//            out_Y[i+QS::iR]      = dRdt;
//            out_Y[i+QS::iRsh]    = dRshdt;
//            out_Y[i+QS::itt]     = dRdt * dttdr;
//            out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+QS::iGamma]  = 0 * dGammadR;
//            out_Y[i+QS::iEint2]  = 0 * dEint2dR;
//            out_Y[i+QS::itheta]  = 0 * dthetadr;
//            out_Y[i+QS::iErad2]  = 0 * dErad2dR;
//            out_Y[i+QS::iEsh2]   = 0 * dEsh2dR;
//            out_Y[i+QS::iEad2]   = 0 * dEad2dR;
//            out_Y[i+QS::iM2]     = 0 * dM2dR;
////            Gamma = p_pars->Gamma0;
//        }
//        else{
//            out_Y[i+QS::iR]      = dRdt;//1.0 / betaSh / CGS::c;
//            out_Y[i+QS::iRsh]    = dRshdt;//1.0 / betaSh / CGS::c;
//            out_Y[i+QS::itt]     = dRdt * dttdr;
//            out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+QS::iGamma]  = dRdt * dGammadR;
//            out_Y[i+QS::iEint2]  = dRdt * dEint2dR;
//            out_Y[i+QS::itheta]  = dRdt * dthetadr;
//            out_Y[i+QS::iErad2]  = dRdt * dErad2dR;
//            out_Y[i+QS::iEsh2]   = dRdt * dEsh2dR;
//            out_Y[i+QS::iEad2]   = dRdt * dEad2dR;
//            out_Y[i+QS::iM2]     = dRdt * dM2dR;
//        }
    out_Y[i+ SOL::QS::iR]      = dRdt;//1.0 / betaSh / CGS::c;
//    out_Y[i+ SOL::QS::iRsh]    = dRshdt;//1.0 / betaSh / CGS::c;
    out_Y[i+ SOL::QS::itt]     = dRdt * dttdr;
    out_Y[i+ SOL::QS::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i+ SOL::QS::imom]    = dRdt * dmomdR; //dRdt * dGammadR / Gamma;
    out_Y[i+ SOL::QS::iGamma]  = dRdt * dGammadR; //dRdt * dGammadR / Gamma;
    out_Y[i+ SOL::QS::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+QS::iEinj]   = dRdt * dEingdR_abs;
    out_Y[i+ SOL::QS::itheta]  = dRdt * dthetadr;
    out_Y[i+ SOL::QS::iErad2]  = dRdt * dErad2dR;
    out_Y[i+ SOL::QS::iEsh2]   = dRdt * dEsh2dR;
    out_Y[i+ SOL::QS::iEad2]   = dRdt * dEad2dR;
    out_Y[i+ SOL::QS::iM2]     = dRdt * dM2dR;
//        if (Gamma)
//        out_Y[i+QS::iR]      = dRdt;//1.0 / betaSh / CGS::c;
//        out_Y[i+QS::iRsh]    = dRshdt;//1.0 / betaSh / CGS::c;
//        out_Y[i+QS::itt]     = dRdt * dttdr;
//        out_Y[i+QS::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i+QS::iGamma]  = dRdt * dGammadR;
//        out_Y[i+QS::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+QS::itheta]  = dRdt * dthetadr;
//        out_Y[i+QS::iErad2]  = dRdt * dErad2dR;
//        out_Y[i+QS::iEsh2]   = dRdt * dEsh2dR;
//        out_Y[i+QS::iEad2]   = dRdt * dEad2dR;
//        out_Y[i+QS::iM2]     = dRdt * dM2dR;
    // ****************************************
//        if (betaSh<1e-5)  p_pars->end_evolution = true;
}

/// RHS with ODEs for dGammadR modified for pre-accelerated ISM with energy injection
void BlastWave::rhs_fs_dense_pwn(double * out_Y, size_t i, double x, double const * Y ) {
    // *************| Extract Values |***************************
//        double mom    = Y[i+ SOL::QS::imom];
    double Gamma  = Y[i+ SOL::QS::iGamma];
    double R      = Y[i+ SOL::QS::iR];
//    double Rsh    = Y[i+ SOL::QS::iRsh];
    double tcomov = Y[i+ SOL::QS::itcomov];
    double Eint2  = Y[i+ SOL::QS::iEint2];
    double theta  = Y[i+ SOL::QS::itheta];
    double M2     = Y[i+ SOL::QS::iM2];
//        double Gamma  = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
//        double betaSh   = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
    double beta   = BetaFromGamma(Gamma);
    /// ---- PWN ---
    double r_w   = Y[i+ SOL::QS::iRw];
    double e_nb  = Y[i+ SOL::QS::iWenb];
    double e_pwn = Y[i+ SOL::QS::iWepwn];


    // *************| Check Values |****************************
    if (Gamma < 1.){
        (*p_log)(LOG_ERR,AT) << "Error mom < 0! Resetting to default \n";
        Gamma = 1. + 1e-5;
//            Gamma = EQS::GamFromMom(Y[i+ SOL::QS::imom]);
//            betaSh = EQS::BetFromMom(Y[i+ SOL::QS::imom]);
        beta = BetaFromGamma(Gamma);
    }
    if ((!std::isfinite(R)) || (!std::isfinite(Gamma)) ||
        (M2 < 0.) || (!std::isfinite(M2)) || (!std::isfinite(Eint2)) || (Eint2<0)) {
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " << "\n"
                              << " R="<<R<<"\n"
                              << " M2="<<M2<<"\n"
                              << " Gamma=" << Gamma << "\n"
                              //                                  << " Mom=" << mom << "\n"
                              << " Eint2=" << Eint2
                              << " \n";
        exit(1);
    }

    // ****| Compute the Downstream properties |*********************
    double bw_thickness = -1; // TODO interpolate it from BM-ST self-similar solution
    double bw_rho = -1;
    double bw_temp = -1;
    double bw_tau = -1;
    double bw_ye = -1;
    if (p_pars->method_thick_for_rho == iFromRadius){
        /// simple estimation of the properties of the blastwave
        bw_thickness = R;
        double vol = (4./3.) * CGS::pi * R * R * R;
        bw_rho = p_pars->M0 * (1. + M2) / vol;
        double ene_th = Eint2 * p_pars->M0 * CGS::c * CGS::c;
        bw_temp = std::pow(ene_th / A_RAD / vol, 0.25); // Stephan-Boltzman law
        double bw_kappa = 5.;//p_pars->kappa;
        bw_tau = bw_kappa * bw_rho * bw_thickness;
        bw_ye = p_pars->Ye0;
    }

    // ****| Check the downstream conditions |******************
    if ((!std::isfinite(bw_thickness)) || (!std::isfinite(bw_rho)) || (bw_thickness < 0.)
        || (!std::isfinite(bw_temp)) || (bw_tau < 0) || ((bw_ye < 0))) {
        (*p_log)(LOG_ERR,AT) << " bad value downstream " << "\n"
                             << " bw_thickness="<<bw_thickness<<"\n"
                             << " bw_rho="<<bw_rho<<"\n"
                             << " bw_temp=" << bw_temp << "\n"
                             << " bw_tau=" << bw_tau << "\n"
                             << " \n";
        exit(1);
    }

    // ****| Update the PWN state |***************************
    double ldip = p_pars->curr_ldip;
    double lacc = p_pars->curr_lacc;
    if (r_w > R){ r_w = R; } // PWN should be bounded by the blastwave
    double v_w = std::sqrt( (7./6.) * e_nb / (4. * CGS::pi * r_w*r_w*r_w * bw_rho) );// (see Eq. 28 in Kashiyama+16)
    if (v_w > beta*CGS::c){
        v_w = beta*CGS::c; } // velocity should also be bounded by BW TODO see if not
    double mom_wind = EQS::MomFromBeta(v_w/CGS::c);
    double dEnbdt = 0; // computeSynchrotronEmissivityAbsorptionAnalytic nebula energy \int(Lem * min(1, tau_T^ej * V_ej / c))dt
    if (bw_tau * (r_w / R) > CGS::c / v_w){ // Eq.[28] in Eq. 28 in Kashiyama+16
        dEnbdt = (p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
    }
    else{
        dEnbdt = (bw_tau * (r_w / R) * (v_w/CGS::c) * p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
    }
    double tdyn = R / beta / CGS::c; // dyn. timescale
    double dEpwndt = p_pars->eps_mag_w * ldip - e_pwn / tdyn; // adiabatic loss inclided
    double dttdr_w = EQS::evalElapsedTime(r_w, mom_wind, 0., false);
    if (!std::isfinite(dttdr_w)||(dttdr_w<=0)){
        (*p_log)(LOG_ERR,AT)<<"dttdr_w="<<dttdr_w<<"\n";
        exit(1);
    }
    // using the downstream properties, compute the fraction of PWN energy thermlized in BW
    double fac_psr_dep_tmp = getFacPWNdep(bw_rho,bw_thickness,bw_temp,bw_ye); // double rho_ej, double delta_ej, double T_ej, double Ye
    double dEindt = fac_psr_dep_tmp * (p_pars->eps_e_w * ldip + p_pars->epsth_w * lacc);
    if (!std::isfinite(dEindt)||(dEindt<0)){
        (*p_log)(LOG_ERR,AT)<<"dEindt="<<dEindt<<"\n";
        exit(1);
    }

    // ****| Get the Upstream properties |******************

    double rho = p_dens->m_rho_ / p_pars->M0; // TODO interpolate from ejecta structure as it moves freely through ISM
    double drhodr = p_dens->m_drhodr_ / p_pars->M0;
    double GammaRho = p_dens->m_GammaRho;
    double GammaRel = p_dens->m_GammaRel; // Gamma
    double dGammaRelDGamma = p_dens->m_dGammaReldGamma; // 1.
    double dGammaRhodR = p_dens->m_dGammaRhodR; // 0.
    double cs_cbm = p_dens->m_CS_CBM;
//        drhodr = 0;
//        GammaRho = 1;
//        GammaRel = Gamma; // Gamma
//        dGammaRelDGamma = 1; // 1.
//        dGammaRhodR = 0.; // 0.
    if ((rho <= 0)||(!std::isfinite(drhodr))||(GammaRho<1)||(GammaRel<0)||(!std::isfinite(dGammaRhodR))){
        (*p_log)(LOG_ERR,AT) << "Wrong value upstream"
                             << " rho="<<rho<< " drhodr="<<drhodr<< " GammaRho="<<GammaRho<< " GammaRel="<<GammaRel
                             << " dGammaRelDGamma="<<dGammaRelDGamma<< " dGammaRhodR="<<dGammaRhodR<<"\n";
        exit(1);
    }

    // ***| Check the upstream conditions |******************
    if (GammaRho < 1.) {
        (*p_log)(LOG_ERR,AT)  << "GammaRho=" << GammaRho  << "\n";
//            std::cerr << AT<< "\n";
        exit(1);
    }
    if (GammaRel < 1.) {
        (*p_log)(LOG_ERR,AT)  << "GammaRel=" << GammaRel << "\n"; GammaRel = Gamma;
//            std::cerr << AT<< "\n";
        exit(1);
    }
    if (!std::isfinite(dGammaRelDGamma)) { (*p_log)(LOG_ERR,AT) << "dGammaRelDGamma=" << dGammaRelDGamma << "\n"; exit(1); }
    if (!std::isfinite(dGammaRhodR)) { (*p_log)(LOG_ERR,AT)  << "dlnGammaCBMdR=" << dGammaRhodR << "\n"; exit(1); }

    // ***| Update BW state |********************************
    double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...
    double Gamma_ = GammaRel;
    double beta_ = EQS::BetaFromGamma(Gamma_);
    double gammaAdi  = p_eos->getGammaAdi(Gamma_, beta_);//p_eos->getGammaAdi(Gamma, betaSh);
    double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
    double dRshdt = EQS::BetaFromGamma(GammaSh) * CGS::c;
    double dthetadr = 0.0; // NO lateral expansion for this BW

    double dM2dR = 0.;
    switch (p_pars->m_method_dmdr) {
        case iusingA:
            /// Old equation from Gavin P. Lamb work
            dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
            break;
        case iusingdthdR:
            /// Equation motivated by Gavin P. Lamb paper
            dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
            break;
        case iNodmdr:
            /// if no mass needs to be swept up
            break;
    }
    if (dM2dR < 0.){
        (*p_log)(LOG_ERR,AT) << " dMdR < 0 in RHS dyn for kN ejecta\n";
        exit(1);
    }

    // PWN energy injection
    double xi_inj = 1.; // Additional free parameters (used in dynamics eq.)
    double dEinjdt = dEindt / (p_pars->M0 * CGS::c * CGS::c) / p_pars->ncells;
    double dEinjdR = dEinjdt / dRdt;
    double theta_ej = 0.; // assume that ejecta is alinged with magnetar emission?..
    double Doppler = Gamma / (1. - beta * std::cos(theta_ej));
    double dEinjdR_dop = dEinjdR * Doppler;
    double dEingdR_abs = dEinjdR;// * ( 1. - std::exp(-1.*p_pars->dtau) ) * std::exp(-1.*p_pars->tau_to0);
    double dEingdR_abs_dop = dEingdR_abs / Doppler / Doppler; // in comoving frame

    // Nuclear heating by r-process
    double dEnuc = p_pars->dEnuc;//->getPars()->eps_nuc_thermalized;
    double dEnuc_norm = dEnuc / (p_pars->M0 * CGS::c * CGS::c);
    double dEnucdR = dEnuc_norm / dRdt;

    // Energy loss to thermal microphysics
    double dElum = p_pars->dElum;
    double dElum_norm = dElum / (p_pars->M0 * CGS::c * CGS::c);
    double dElumdR = dElum_norm / dRdt;
    if ((std::abs(dElumdR) > std::abs(dEnucdR*1.5)) and (dEnuc > 0.)) // TODO These are artificial limiters
        dElumdR = dEnucdR;

    // Dynamics equation
    double GammaEff = get_GammaEff(Gamma_, gammaAdi); // TODO check if GammaRel should be used for this!!!
    double dGammaEffdGamma = get_dGammaEffdGamma(Gamma_, gammaAdi);
    double num1 = (Gamma - GammaRho + GammaEff * (GammaRel - 1.)) * dM2dR;
    double num2 = - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR/M2 - drhodr/rho - dGammaRhodR / GammaRho); // - 3.*Eint2/R
    double num3 = - (1. + GammaEff / Gamma * xi_inj) * dEingdR_abs; // injection term
    double denum1 = (1. + M2);
    double denum2 = Eint2 * dGammaEffdGamma;
    double denom3 = GammaEff * (gammaAdi - 1.) * Eint2 * dGammaRelDGamma / GammaRel;
    double dGammadR = -1. * (num1 + num2 + num3) / (denum1+denum2+denom3);
//        double dmomdR = dGammadR / betaSh;

    // Adiabatic energy loss
    double dlnV2dR  = dM2dR / M2 - drhodr / rho - (1. / GammaRel) * dGammaRelDGamma * dGammadR + dGammaRhodR / GammaRho;// + 3./R;
    double dlnVdR = (3./R) - (1/Gamma * dGammadR);
    double dEad2dR  = 0.0;
    if ( p_pars->adiabLoss )
        dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;

    // Energy generation at the shock
    double dEsh2dR = (GammaRel - 1.0) * dM2dR;
    if ((cs_cbm > EQS::BetaFromGamma(GammaRel))){//&&(p_pars->comp_ix > p_pars->nr*1e-2))){ // TODO Subsonic flow -- no shock
        dEsh2dR *= 1e-10;
    }

    // -- Radiative losses
    double dErad2dR = p_pars->eps_rad * dEsh2dR;

    // -- Energy equation
    double dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEnucdR - dElumdR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)

    double dtcomov_dR = 1.0 / beta / Gamma / CGS::c; // comoving time

    // -- time in a observer frame
    bool spread = p_spread->m_method != LatSpread::METHODS::iNULL;
    double dttdr = EQS::evalElapsedTime(R,Gamma*beta,dthetadr,spread);

    // ***| final checks |*************************

    if (Gamma < 1.){
        (*p_log)(LOG_ERR, AT) << " Gamma < 1 = "<< Gamma << " in kN RHS dynamics\n";
//            exit(1);
    }
    if ((!std::isfinite(dRdt)) || (!std::isfinite(dGammadR)) || (dM2dR < 0.)
        || (!std::isfinite(dlnV2dR)) || (!std::isfinite(dthetadr))) {
        (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " //<< "\n"
                              << " dRdt="<<dRdt//<<"\n"
                              << " dM2dR="<<dM2dR//<<"\n"
                              << " dthetadr=" << dthetadr// << "\n"
                              << " dGammadR=" << dGammadR //<< "\n"
                              << " dthetadr=" << dRdt //<< "\n"
                              << " dEsh2dR=" << dEsh2dR //<< "\n"
                              << " dlnV2dR=" << dlnV2dR //<< "\n"
                              << " dEad2dR=" << dEad2dR //<< "\n"
                              << " dErad2dR=" << dErad2dR //<< "\n"
                              << " dEint2dR=" << dEint2dR //<< "\n"
                              << " dttdr=" << dttdr// << "\n"
                              << " dGammaRelDGamma=" << dGammaRelDGamma //<< "\n"
                              << " dGammaRhodR=" << dGammaRhodR //<< "\n"
                              << " drhodr=" << drhodr //<< "\n"
                              << " Gamma="<<Gamma//<< "\n"
                              << " betaSh="<<beta//<< "\n"
                              << " theta"<<theta//<< "\n"
                              << " M2="<<M2//<< "\n"
                              << " Eint2="<<Eint2//<< "\n"
                              << " GammaRel"<<GammaRel//<< "\n"
                              << " GammaRho"<<GammaRho//<< "\n"
                              << " \n";
        exit(1);
    }

    // ***| Save the state |***********************
    p_pars->prev_x = x;
    out_Y[i+ SOL::QS::iR]      = dRdt;//1.0 / betaSh / CGS::c;
//    out_Y[i+ SOL::QS::iRsh]    = dRshdt;//1.0 / betaSh / CGS::c;
    out_Y[i+ SOL::QS::itt]     = dRdt * dttdr;
    out_Y[i+ SOL::QS::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i+ SOL::QS::imom]    = dRdt * dmomdR; //?dRdt * dGammadR / Gamma;
    out_Y[i+ SOL::QS::iGamma]  = dRdt * dGammadR; //dRdt * dGammadR / Gamma;
    out_Y[i+ SOL::QS::iEint2]  = dRdt * dEint2dR;
    out_Y[i+ SOL::QS::itheta]  = dRdt * dthetadr;
    out_Y[i+ SOL::QS::iErad2]  = dRdt * dErad2dR;
    out_Y[i+ SOL::QS::iEsh2]   = dRdt * dEsh2dR;
    out_Y[i+ SOL::QS::iEad2]   = dRdt * dEad2dR;
    out_Y[i+ SOL::QS::iM2]     = dRdt * dM2dR;
    // --- PWN ---
    double drdtw = v_w;
    out_Y[i+ SOL::QS::iRw]     = drdtw;
    out_Y[i+ SOL::QS::iWtt]    = dttdr_w*drdtw;
    out_Y[i+ SOL::QS::iWenb]   = dEnbdt;
    out_Y[i+ SOL::QS::iWepwn]  = dEpwndt;
}

#endif //SRC_BLASTWAVE_H
