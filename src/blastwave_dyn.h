//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_DYN_H
#define SRC_BLASTWAVE_DYN_H

#include "pch.h"
#include "utils.h"
#include "base_equations.h"
#include "interpolators.h"
#include "ode_solvers.h"
#include "quadratures.h"
#include "rootfinders.h"
#include "observer.h"
#include "synchrotron_an.h"

#include "blastwave_base.h"
#include "blastwave_rad.h"



/// specific set of equations for a blast wave
class DynRadBlastWave : public RadBlastWave{
    std::unique_ptr<logger> p_log;
public:
    explicit DynRadBlastWave(Vector & tb_arr, size_t ishell, size_t ilayer, int loglevel )
            : RadBlastWave(tb_arr, ishell, ilayer, loglevel) {
//        p_log = std::make_unique<logger>(std::cout, loglevel, "RadBlastWave");
//        pfsolvePars = new FsolvePars;
//        pfsolvePars->p_eos = p_eos;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "DynRadBlastWave");
    }
    ~DynRadBlastWave(){ }
//    FsolvePars *& getCombPars(){ return pfsolvePars; }
    // RHS settings
    const static int neq = 11;
    std::vector<std::string> vars { "R", "Rsh", "tt", "tcomov", "mom", "Eint2", "theta", "Erad2", "Esh2", "Ead2", "M2" };
    enum Q_SOL { iR, iRsh, itt, itcomov, imom, iEint2, itheta, iErad2, iEsh2, iEad2, iM2 };
    size_t getNeq() override { return neq; }
    enum CASES { i_INSIDE_BEHIND, i_OUTSIDE_BEHIND, i_INSIDE_ABOVE, i_OUTSIDE_ABOVE, i_AT_ZERO_INSIDE };
    /// set initial condition 'ic_arr' for this blast wave using data from Pars{} struct
    void setInitConditions( double * ic_arr, size_t i ) override {


        /// if layer does not have energy / mass -- do not evolve it
        if ((p_pars->M0 == 0.) && (p_pars->E0 == 0.)){
//            std::cerr << AT << "\n "
//                      << "[ishell=" << p_pars->ishell << " ilayer="<<p_pars->ilayer << "] "
//                      <<" M0=0 and E0=0 -> Ignoring this layer.\n";
            p_pars->end_evolution = true;
            for (size_t v = 0; v < neq; ++v){
                ic_arr[i+v] = 0.;
            }
            return;
        }

        // ****************************************
        double beta0 = EQS::Beta(p_pars->Gamma0);
        p_pars->R0    = p_pars->tb0 * beta0 * CGS::c;

//        double _R0 = 1e8;
//        double _tb0 = _R0 / (beta0 * CGS::c);
//        m_mag_time = m_mag_time - m_mag_time[0] + _tb0;
//        p_pars->tb0 = m_mag_time[0];
//        p_pars->x = p_pars->tb0;
//        p_pars->R0    = p_pars->tb0 * beta0 * CGS::c;


//        double j_theta  = Y[other_i + DynRadBlastWave::Q_SOL::itheta];
//        double j_ctheta = other->ctheta(j_theta);

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

        double m_M20 = (2.0 / 3.0)
                     * CGS::pi
                     * (std::cos(p_pars->theta_a) - std::cos(p_pars->theta_b0))
                     * p_dens->m_rho_
                     * std::pow(p_pars->R0, 3)
                     / p_pars->ncells;
        double adi0 = p_eos->getGammaAdi(p_pars->Gamma0,EQS::Beta(p_pars->Gamma0));
        double GammaSh0 = EQS::GammaSh(p_pars->Gamma0,adi0);
        // ****************************************
        if ((p_pars->M0 <= 0.) || (!std::isfinite(p_pars->M0))){
//            std::cout << "[ WARNING ] " << "M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            // REMOVING LOGGER
            (*p_log)(LOG_WARN, AT)  << " M0 < 0 Setting M0=E0/(Gamma0 c^2)\n";
            p_pars->M0 = p_pars->E0 / (p_pars->Gamma0 * CGS::c * CGS::c);
        }
        if ((p_pars->R0 <= 1.) || (!std::isfinite(p_pars->R0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " R0 <= 0 (R0=" <<p_pars->R0 << ") " << "G0="<<p_pars->Gamma0
                       << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<") \n"
                       << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "R0 <= 0 (R0=" <<R0 << ")\n";
            exit(1);
        }
        if ((p_dens->m_rho_) <= 0.|| (!std::isfinite(p_dens->m_rho_))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " rho0 < 0 (rho0=" <<p_dens->m_rho_ << ") "
                       << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                       << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                       << " \n";
//            std::cerr << AT  << "\n";

            //std::cout << "[ Error ] " << "rho0 < 0 (rho0=" <<rho0 << ")\n";
            exit(1);
        }
        if ((p_pars->E0 <= 0.) || (!std::isfinite(p_pars->E0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << "E0 <= 0 (E0=" <<p_pars->E0 << ") "
                       << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0
                       << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                       << " \n";
//            std::cerr << AT  << "\n";
            //std::cout << "[ Error ] " << "E0 < 0 (E0=" <<E0 << ")\n";
            exit(1);
        }
        if ((p_pars->Gamma0 < 1.) || (!std::isfinite(p_pars->Gamma0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " Gamma0 < 1 (Gamma0=" <<p_pars->Gamma0 << ") "
                       << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                       << " \n";
            //std::cout << "[ Error ] " << "Gamma0 < 0 (Gamma0=" <<Gamma0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->theta_b0 < 0.) || (!std::isfinite(p_pars->theta_b0))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " theta_b0 < 0 (theta_b0=" <<p_pars->theta_b0 << ") "
                       << "G0="<<p_pars->Gamma0 << " E0="<<p_pars->E0 << " tb0="<<p_pars->tb0 << " (offset i="<<i<<")\n"
                       << " \n";
            //std::cout << "[ Error ] " << "theta_b0 < 0 (theta_b0=" <<theta_b0 << ")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        if ((p_pars->ncells < 1 )|| (!std::isfinite(p_pars->ncells))){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " ncells < 1 (ncells=" <<p_pars->ncells << ")\n"
                       << " \n";
            //std::cout << "[ Error ] " << "ncells < 1 (ncells=" <<ncells << ")\n";
            exit(1);
        }
        if (cos(p_pars->theta_a) <= cos(p_pars->theta_b0)){
            // REMOVING LOGGER
            (*p_log)(LOG_ERR, AT)  << " cos(theta_a) < cos(theta_b0) (theta_a="<<p_pars->theta_a
                       << ", theta_b0="<<p_pars->theta_b0<<")\n"
                       << " \n";
            //std::cout << "[ Error ] " <<" cos(theta_a) < cos(theta_b0) (theta_a="<<theta_a<<", theta_b0="<<theta_b0<<")\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
        bool use_spread = p_spread->m_method != LatSpread::METHODS::iNULL;
        p_pars->Rd = std::pow(3./(4.*CGS::pi) * 1./(CGS::c*CGS::c*CGS::mp) *
                              p_pars->E0/( p_dens->m_rho_ / CGS::mp * p_pars->Gamma0*p_pars->Gamma0), 1./3.);
//        if (p_pars->Rd < p_pars->R0){
//            std::cerr << AT << " error with Rd\n";
//            exit(1);
//        }
        double x = p_pars->Rd / p_pars->R0;
        // ***************************************
        ic_arr[i + Q_SOL::iR]      = m_tb_arr[0] * beta0 * CGS::c;
        ic_arr[i + Q_SOL::iRsh]    = m_tb_arr[0] * EQS::Beta(GammaSh0) * CGS::c;
        ic_arr[i + Q_SOL::itt]     = EQS::init_elapsed_time(p_pars->R0, p_pars->Gamma0, use_spread);
        ic_arr[i + Q_SOL::itcomov] = p_pars->R0 / (beta0 * p_pars->Gamma0 * CGS::c);
        ic_arr[i + Q_SOL::iEint2]  = (p_pars->Gamma0 - 1 ) * m_M20 / p_pars->M0;  //TODO Isnt it just E0 / m_M0 ???? As M0 = E0 * cgs.c ** -2 / Gamma0
        ic_arr[i + Q_SOL::imom]    = p_pars->Gamma0 * EQS::Beta(p_pars->Gamma0);//std::log( p_pars->Gamma0 );
//        ic_arr[i + Q_SOL::iGamma]  = p_pars->Gamma0;//std::log( p_pars->Gamma0 );
        ic_arr[i + Q_SOL::itheta]  = p_pars->theta_b0;
        ic_arr[i + Q_SOL::iErad2]  = 0.0;
        ic_arr[i + Q_SOL::iEsh2]   = 0.0;
        ic_arr[i + Q_SOL::iEad2]   = 0.0;
        ic_arr[i + Q_SOL::iM2]     = m_M20 / p_pars->M0;
        // ***************************************
        for (size_t v = 0; v < neq; ++v){
            if (!std::isfinite(ic_arr[i + v])){
                (*p_log)(LOG_ERR, AT)  << " NAN in initial data for shell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                           << " v_n="<<vars[v]<<" val="<<ic_arr[i + v]<<" Exiting...\n";
//                std::cerr << AT  << "\n";
                exit(1);
            }
        }
//        p_spread->m_theta_b0 = p_pars->theta_b0;
//        p_pars->x = p_pars->tb0;
        // ***************************************
    }
    /// add the current solution 'sol' to the 'm_data' which is Vector of Arrays (for all variables)
    void insertSolution( const double * sol, size_t it, size_t i ) override {
        if (p_pars->end_evolution)
            return;
//        double mom = sol[i+Q_SOL::imom];
//        double gam = sol[i+Q_SOL::iGamma];
        p_pars->comp_ix = it; // latest computed iteration
        m_data[Q::itburst][it]   = m_tb_arr[it];
        m_data[Q::iR][it]        = sol[i+Q_SOL::iR]; // TODO you do not need 'i' -- this is p_pars->ii_eq
        m_data[Q::iRsh][it]      = sol[i+Q_SOL::iRsh];
        m_data[Q::itt][it]       = sol[i+Q_SOL::itt];
        m_data[Q::imom][it]      = sol[i+Q_SOL::imom];
        m_data[Q::iGamma][it]    = EQS::GamFromMom(sol[i+Q_SOL::imom]);//sol[i+Q_SOL::iGamma];
        m_data[Q::itheta][it]    = sol[i+Q_SOL::itheta];
        m_data[Q::iM2][it]       = sol[i+Q_SOL::iM2];
        m_data[Q::itcomov][it]   = sol[i+Q_SOL::itcomov];
        m_data[Q::iEad2][it]     = sol[i+Q_SOL::iEad2];
        m_data[Q::iEint2][it]    = sol[i+Q_SOL::iEint2];
        m_data[Q::iEsh2][it]     = sol[i+Q_SOL::iEsh2];
        m_data[Q::iErad2][it]    = sol[i+Q_SOL::iErad2];
        if (sol[i+Q_SOL::iR] < 1. || m_data[Q::iGamma][it] < 1. || sol[i+Q_SOL::imom] < 0) {
            (*p_log)(LOG_ERR, AT)  << "Wrong value at i=" << it << " tb=" << sol[i + Q_SOL::iR]
                       << " iR="      << sol[i + Q_SOL::iR]
                       << " iRsh="    << sol[i + Q_SOL::iRsh]
//                       << " imom="  << m_data[Q::iGamma][it]
                       << " iMom="    << sol[i+Q_SOL::imom]
                       << " itheta="  << sol[i + Q_SOL::itheta]
                       << " iM2="     << sol[i + Q_SOL::iM2]
                       << " itcomov=" << sol[i + Q_SOL::itcomov]
                       << " iEad2="   << sol[i + Q_SOL::iEad2]
                       << " iEint2="  << sol[i + Q_SOL::iEint2]
                       << " iEsh2="   << sol[i + Q_SOL::iEsh2]
                       << " iErad2="  << sol[i + Q_SOL::iErad2]
                       << " "
                       << " Exiting...\n";
//            std::cerr << AT  << "\n";
            exit(1);
        }
    }
    /// Mass and energy are evolved in units of M0 and M0c2 respectively
    void applyUnits( double * sol, size_t i ) override {
        sol[i + Q_SOL::iM2]    *= p_pars->M0;
        sol[i + Q_SOL::iEint2] *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + Q_SOL::iEad2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + Q_SOL::iEsh2]  *= (p_pars->M0 * CGS::c * CGS::c);
        sol[i + Q_SOL::iEad2]  *= (p_pars->M0 * CGS::c * CGS::c);
    }
    /// check if to terminate the evolution
    bool isToTerminate( double * sol, size_t i ) override {
        double mom = sol[i+Q_SOL::imom];
//        double igamma = sol[i+Q_SOL::iGamma];//EQS::GamFromMom(mom);
        double igamma = EQS::GamFromMom(mom);
        double ibeta = EQS::Beta(igamma);//EQS::BetFromMom(mom);
        double iEint2 = sol[i + Q_SOL::iEint2];
        bool do_terminate = false;
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
//            double beta = EQS::Beta(igamma);
            if ((iEint2 <= 0.)||(!std::isfinite(iEint2))||(mom < 0)||(ibeta < 0)){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: Eint2 < 0 or NAN ("<<iEint2<<")\n";
            }
            if ((std::abs(ibeta) < p_pars->min_beta_terminate)||(!std::isfinite(ibeta))){
                do_terminate = true;
                (*p_log)(LOG_ERR,AT)
                        << " Terminating evolution [ishell=" << p_pars->ishell << " ilayer=" << p_pars->ilayer << "] "
                        << " REASON :: beta < beta_min or NAN ("<<ibeta<<") with beta_min= "
                        << p_pars->min_beta_terminate<<"\n";
            }
            if (do_terminate) {
                (*p_log)(LOG_WARN, AT) << " TERMINATION Previous iterations: \n"
                                      << " E0=" << string_format("%.9e", p_pars->E0)
                                      << " Gamma0=" << string_format("%.9e", p_pars->Gamma0)
                                      << " M0=" << string_format("%.9e", p_pars->M0)
                                      << " \n";
                (*p_log)(LOG_WARN, AT) << " i=Current" << " Eint2=" << string_format("%.9e", sol[i + Q_SOL::iEint2])
                                      << " Gamma=" << string_format("%.9e", igamma)
                                      << " beta=" << string_format("%.9f", ibeta)
                                      << " M2=" << string_format("%.9e", sol[i + Q_SOL::iM2])
                                      << " rho=" << string_format("%.9e", p_dens->m_rho_)
                                      << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix << " Eint2="
                                      << string_format("%.9e", getVal(Q::iEint2, p_pars->comp_ix))
                                      << " Gamma=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix))
                                      << " beta=" << string_format("%.9f", getVal(Q::ibeta, p_pars->comp_ix))
                                      << " M2=" << string_format("%.9e", getVal(Q::iM2, p_pars->comp_ix))
                                      << " rho=" << string_format("%.9e", getVal(Q::irho, p_pars->comp_ix))
                                      << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix - 1 << " Eint2="
                                      << string_format("%.9e", getVal(Q::iEint2, p_pars->comp_ix - 1))
                                      << " Gamma=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 1))
                                      << " beta=" << string_format("%.9f", getVal(Q::ibeta, p_pars->comp_ix - 1))
                                      << " M2=" << string_format("%.9e", getVal(Q::iM2, p_pars->comp_ix - 1))
                                      << " rho=" << string_format("%.9e", getVal(Q::irho, p_pars->comp_ix - 1))
                                      << " \n";
                (*p_log)(LOG_WARN, AT) << " i=" << p_pars->comp_ix - 2 << " Eint2="
                                      << string_format("%.9e", getVal(Q::iEint2, p_pars->comp_ix - 2))
                                      << " Gamma=" << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 2))
                                      << " beta=" << string_format("%.9f", getVal(Q::ibeta, p_pars->comp_ix - 2))
                                      << " M2=" << string_format("%.9e", getVal(Q::iM2, p_pars->comp_ix - 2))
                                      << " rho=" << string_format("%.9e", getVal(Q::irho, p_pars->comp_ix - 2))
                                      << " \n";
            }
        }
        return do_terminate;
    }
    /// check if to terminate the evolution
    bool isToStopLateralExpansion( double * sol, size_t i ) override {
        double itheta = sol[i + Q_SOL::itheta];
        double iEint2 = sol[i + Q_SOL::iEint2];
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
    /// check if the solution 'makes sense'
    bool isSolutionOk( double * sol, size_t i ) override {

        bool no_issues = true;
        double mom = sol[i+Q_SOL::imom];
        double mom0 = p_pars->Gamma0 * EQS::Beta(p_pars->Gamma0);
//        double igamma = sol[i+Q_SOL::iGamma];
        double igamma = EQS::GamFromMom(mom);
//        double ibeta = EQS::Beta(igamma);//EQS::BetFromMom(mom);
        double ibeta = EQS::BetFromMom(mom);
        if (!p_pars->end_evolution) {
            /// if BW is too slow numerical issues arise (Gamma goes below 1) # TODO rewrite ALL eqs in terms of GammaBeta
            double beta   = EQS::Beta(igamma);

            if ((igamma < 1.) || (mom < 0)) {
                // REMOVING LOGGER

                //            std::cout << "theta:" << getData(Q::iR)<<"\n";
                (*p_log)(LOG_ERR,AT)  << AT << " \n"
                           << " Gamma < 1. ishell="<<p_pars->ishell<<" ilayer="<<p_pars->ilayer
                           << " with Gamma0="<<p_pars->Gamma0 << " [iteration="<<p_pars->comp_ix<<"] \n";
                (*p_log)(LOG_ERR,AT) << " Last Gammas: " << "\n";
                (*p_log)(LOG_ERR,AT)<< " Gamma[" << p_pars->comp_ix << "]="
                           << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix << "]="
                           << string_format("%.9e", getVal(Q::ibeta, p_pars->comp_ix))
                           << " Gamma[" << p_pars->comp_ix - 1 << "]="
                           << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 1)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 1 << "]="
                           << string_format("%.9e", getVal(Q::ibeta, p_pars->comp_ix - 1))
                           << " Gamma[" << p_pars->comp_ix - 2 << "]="
                           << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 2)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 2 << "]="
                           << string_format("%.9e", getVal(Q::ibeta, p_pars->comp_ix - 2))
                           << " Gamma[" << p_pars->comp_ix - 3 << "]="
                           << string_format("%.9e", getVal(Q::iGamma, p_pars->comp_ix - 3)) << "\n";
                (*p_log)(LOG_ERR,AT)<< " beta[" << p_pars->comp_ix - 3 << "]="
                           << string_format("%.9e", getVal(Q::ibeta, p_pars->comp_ix - 3))
                           << "\n";
//                std::cerr << "Gamma:" << getData(Q::iGamma) << "\n";
//                std::cerr << "R:" << getData(Q::iR) << "\n";
                (*p_log)(LOG_ERR,AT)  << " Gamma cannot be less than 1. or too large. "
                                         "Found: Gamma(" << igamma << "). Gamma0="<<p_pars->Gamma0<<" \n";
                no_issues = false;
            }
//            if (igamma > p_pars->Gamma0 * 2.) {
//                // REMOVING LOGGER
//                (*p_log)(LOG_ERR,AT)  << " too large value of Gamma(" << igamma
//                           << ") > 2*Gamma0" << p_pars->Gamma0 * 2 << ")\n";
//                no_issues = false;
//            }
            if (sol[i + Q_SOL::itheta] < p_pars->theta_b0 * 0.99) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical decrease in theta(" << sol[i + Q_SOL::itheta]
                          << ") < theta_b0(" << p_pars->theta_b0 << ")\n";
                no_issues = false;
            }
            if (sol[i + Q_SOL::itt] <= 0.) {
                // REMOVING LOGGER
                (*p_log)(LOG_ERR,AT) << " unphysical value of observer time (on axis) itt("
                          << sol[i + Q_SOL::itt] << ") < 0 \n";
                no_issues = false;
            }
        }
        return no_issues;
    }
    /// right hand side of the blast wave evolution equation (ODE)
    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ) override {
        // ****************************************
        double R      = Y[i+Q_SOL::iR];
        double Rsh    = Y[i+Q_SOL::iRsh];
//        double tcomov = Y[i+Q::itcomov];
        double mom    = Y[i+Q_SOL::imom];
//        double Gamma  = Y[i+Q_SOL::iGamma];
        double Eint2  = Y[i+Q_SOL::iEint2];
        double theta  = Y[i+Q_SOL::itheta];
        double M2     = Y[i+Q_SOL::iM2];
        // ****************************************
        if (mom < 0) {
            (*p_log)(LOG_ERR, AT) << " negative momentum\n";
            exit(1);
        }
        double Gamma = EQS::GamFromMom(mom);
//        if (Gamma <= 1.)
//            Gamma = EQS::GamFromMom(mom);
//            Gamma = 1.0000000001;
//        if ((theta < p_pars->theta_b0) || (theta > p_pars->theta_max)){
////            std::cerr << "theta < theta_b0 || theta > theta_max\n";
////            exit(1);
//            theta = p_pars->theta_max;
//        }
        double beta   = EQS::Beta(Gamma);

        // ****************************************
//        auto *_pars = (struct RHS_pars *) rhs_pars;
        double ctheta_ = LatStruct::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);//p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
        p_dens->evaluateRhoDrhoDrDefault(R, ctheta_);
        double rho = p_dens->m_rho_def / p_pars->M0;
        double drhodr = p_dens->m_drhodr_def / p_pars->M0;

//        dlnrho1dR /= p_pars->M0;
//        double _rhoi = p_pars->rho/p_pars->M0;
        // ****************************************
        double dRdt = beta * CGS::c; // TODO this is beta_sh ! ...

        double gammaAdi  = p_eos->getGammaAdi(Gamma, beta);
        double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
        double dRshdt = EQS::Beta(GammaSh) * CGS::c;

        double dthetadr = 0.0;
        if ( (theta < 2*p_pars->theta_max)
             && (R > p_pars->Rd)
             && (Gamma < std::max(2., p_pars->Gamma0*p_pars->fraction_of_Gamma0_when_spread))
             && (!p_pars->end_spreading) ) { // &&(Gamma < p_pars->Gamma0*.95)&&
            dthetadr = p_spread->getDthetaDr(Gamma, R, gammaAdi, theta);
//            if (dthetadr > 0.){
//                int x = 1;
//            }
        }


        double dM2dR = 0;
        switch (p_pars->m_method_dmdr) {

            case iusingA:
                dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
                break;
            case iusingdthdR:
                dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
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
        double dmomdR = dGammadR / beta;


//        if (dGammadR > 0){
//            if ( Gamma > 0.95 * p_pars->Gamma0 ){
//                dGammadR = 0.;
//            }
//        }
        // -- Energies --
        double dEsh2dR  = (Gamma - 1.0) * dM2dR; // Shocked energy;
        double dlnV2dR  = dM2dR / M2 - drhodr / rho - dGammadR / Gamma;
        double dEad2dR  = 0.0;
        if ( p_pars->adiabLoss )
            dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;
        // -- Radiative losses
        double dErad2dR = p_pars->eps_rad * dEsh2dR;
        // -- Energy equation
        double dEint2dR = dEsh2dR + dEad2dR - dErad2dR; // / (m_pars.M0 * c ** 2)
        double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
        double dttdr = 0.;
        if (p_spread->m_method != LatSpread::METHODS::iNULL)
            dttdr = 1. / (CGS::c * beta) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
        else
            dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        // ****************************************

        if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || !std::isfinite(dlnV2dR) || !std::isfinite(dthetadr)) {
            (*p_log)(LOG_ERR,AT)  << " nan in derivatives. Exiting..." << "\n";
            exit(1);
        }
        // ****************************************
//        double theta_c_h = dthetadr * dRdt * (x - p_pars->x);
//        if (theta + theta_c_h >= p_pars->theta_max ){
////            exit(1);
//            dthetadr = 0.;
//        }
        if (!std::isfinite(dRdt) || !std::isfinite(dGammadR) || dM2dR < 0.
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
                       << " R=" <<R << " Rsh=" <<Rsh<< " Gamma=" <<Gamma<< " Eint2=" <<Eint2<< " theta=" <<theta
                       << " M2=" <<M2<< " rho=" <<rho<<" drhodr=" <<drhodr<<" dM2dR=" <<dM2dR<< "\n";
            (*p_log)(LOG_ERR,AT)  << " maybe timestep is too large \n Exiting...";
//            std::cerr << AT << "\n";
            exit(1);
        }
        p_pars->x = x; // update
        out_Y[i + Q_SOL::iR] = dRdt;//1.0 / beta / CGS::c;
        out_Y[i + Q_SOL::iRsh] = dRshdt;//1.0 / beta / CGS::c;
        out_Y[i + Q_SOL::itt] = dRdt * dttdr;
        out_Y[i + Q_SOL::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i + Q_SOL::iGamma] = 0.;//dRdt * dGammadR;//dRdt * dGammadR/Gamma;
        out_Y[i + Q_SOL::imom] = dRdt * dmomdR;//dRdt * dGammadR/Gamma;
        out_Y[i + Q_SOL::iEint2] = dRdt * dEint2dR;
        out_Y[i + Q_SOL::itheta] = dRdt * dthetadr;
        out_Y[i + Q_SOL::iErad2] = dRdt * dErad2dR;
        out_Y[i + Q_SOL::iEsh2] = dRdt * dEsh2dR;
        out_Y[i + Q_SOL::iEad2] = dRdt * dEad2dR;
        out_Y[i + Q_SOL::iM2] = dRdt * dM2dR;
//        if (beta>1e-5)  p_pars->end_evolution = true;
        // ****************************************
//        if (mom+out_Y[i + Q_SOL::imom])
    }
    /// RHS with ODEs for dGammadR modified for pre-accelerated ISM
    void evaluateRhsDens( double * out_Y, size_t i, double x, double const * Y ) override {
//        double Gamma = Y[i+Q_SOL::iGamma];//EQS::GamFromMom(Y[i+Q_SOL::imom]);
        double mom = Y[i+Q_SOL::imom];
        // ****************************************
        double R      = Y[i+Q_SOL::iR];
        double Rsh    = Y[i+Q_SOL::iRsh];
//        double tcomov = Y[i+Q::itcomov];
//        double Gamma  = std::exp(Y[i+Q_SOL::ilnGamma]);
        double Eint2  = Y[i+Q_SOL::iEint2];
        double theta  = Y[i+Q_SOL::itheta];
        double M2     = Y[i+Q_SOL::iM2];
        // *****************************************
        double Gamma = EQS::GamFromMom(Y[i+Q_SOL::imom]);
//        double beta = EQS::Beta(Y[i+Q_SOL::iGamma]);//EQS::BetFromMom(Y[i+Q_SOL::imom]);
        double beta = EQS::BetFromMom(Y[i+Q_SOL::imom]);
        if (mom < 0){
            (*p_log)(LOG_ERR,AT) << "Error\n";
            mom = 1e-5;
            Gamma = EQS::GamFromMom(Y[i+Q_SOL::imom]);
            beta = EQS::BetFromMom(Y[i+Q_SOL::imom]);
        }
        // ****************************************
//        if (!std::isfinite(R) || !std::isfinite(Gamma) || M2 < 0.
//            || !std::isfinite(M2) || !std::isfinite(Eint2) || Eint2<0) {
//            (*p_log)(LOG_ERR,AT)  << " nan in derivatives (may lead to crash) " << "\n"
//                                  << " R="<<R<<"\n"
//                                  << " M2="<<M2<<"\n"
//                                  << " Gamma=" << Gamma << "\n"
//                                  << " Mom=" << mom << "\n"
//                                  << " Eint2=" << Eint2
//                                  << " \n";
//            exit(1);
//        }
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
//        double beta   = EQS::Beta(Gamma);


        // ****************************************
        // Get ISM density and its velocity
//        double rho, m_drhodr;
        double ctheta_ = LatStruct::ctheta(theta,p_pars->ilayer,p_pars->nlayers);//ctheta(theta);// = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w);
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
        double beta_ = EQS::Beta(Gamma_);
        double gammaAdi  = p_eos->getGammaAdi(Gamma_, beta_);//p_eos->getGammaAdi(Gamma, beta);

        double GammaSh = EQS::GammaSh(Gamma,gammaAdi);
        double dRshdt = EQS::Beta(GammaSh) * CGS::c;

        double dthetadr = 0.0;
//        if (theta < (p_pars->theta_max * 0.999999) ) {
//            dthetadr = p_spread->getDthetaDr(Gamma, R, gammaAdi, theta);
//        }
//        double dM2dR     = EQS::dmdr(Gamma, R, p_pars->theta_a,
//                                     theta, rho, p_spread->m_aa) / p_pars->ncells;

        double dM2dR = 0;
        switch (p_pars->m_method_dmdr) {

            case iusingA:
                dM2dR = EQS::dmdr(Gamma, R, p_pars->theta_a, theta, rho, p_spread->m_aa) / p_pars->ncells;
                break;
            case iusingdthdR:
                dM2dR = EQS::dmdr(Gamma, R, dthetadr, theta, rho) / p_pars->ncells;
                break;
        }
//        dM2dR*=0;
        if (dM2dR < 0.){
            (*p_log)(LOG_ERR,AT) << " dMdR < 0 in RHS dyn for kN ejecta\n";
            exit(1);
        }

        // --- Energy injection --- ||
        double xi_inj = 1.;
        double dEindt = p_pars->dEinjdt;

        double dEinjdt = dEindt / (p_pars->M0 * CGS::c * CGS::c) / p_pars->ncells;
        double dEinjdR = dEinjdt / dRdt;
        double theta_ej = 0.; // assume that ejecta is alinged with magnetar emission?..
        double Doppler = Gamma / (1. - beta * std::cos(theta_ej));
        double dEinjdR_dop = dEinjdR * Doppler;
        double dEingdR_abs = dEinjdR;// * ( 1. - std::exp(-1.*p_pars->dtau) ) * std::exp(-1.*p_pars->tau_to0);
        double dEingdR_abs_dop = dEingdR_abs / Doppler / Doppler;


        // --- dGammadR ---
        double dGammadR = 0., GammaEff=0.,dGammaEffdGamma=0.,num1=0.,num2=0.,num3=0.,denum1=0.,denum2=0.,denom3=0.;
        double _tmp = (1. - GammaEff / Gamma * xi_inj) * dEingdR_abs;
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
        double dmomdR = dGammadR / beta;


//        drhodr = 0.; // TODO why was this set?! This overrides the derivative
        // -- Energies --

//        double dlnV2dR  = dM2dR / M2 - m_drhodr - dGammadR / Gamma;
        double dlnV2dR  = dM2dR / M2 - drhodr / rho - (1./GammaRel)*dGammaRelDGamma*dGammadR + dGammaRhodR / GammaRho;// + 3./R;
        double dlnVdR = (3./R) - (1/Gamma * dGammadR);
        double dEad2dR  = 0.0;
        if ( p_pars->adiabLoss )
            dEad2dR = -(gammaAdi - 1.0) * Eint2 * dlnV2dR;
        double dx = dRdt * (x - p_pars->x);
//        double Eint_ = Eint2 + dEad2dR * dRdt * (x - p_pars->x);
//        if (Eint_ < 0){
//            int x = 1;
//        }
//        double mom_ = mom + dmomdR * dRdt * (x - p_pars->x);
//        if (mom_ < 0){
//            int x = 1;
//        }

        // --- shock energy
        double dEsh2dR = (GammaRel - 1.0) * dM2dR;
//        if (cs_cbm > 100){
//            int y = 1;
//        }
        if ((cs_cbm > EQS::Beta(GammaRel)&&(p_pars->comp_ix > p_pars->nr*1e-2))){ // TODO Subsonic flow -- no shock
            dEsh2dR *= 1e-10;
        }
//        double dEsh2dR = 0.;

//        double dEsh2dR  = (GammaRel - 1.0) * dM2dR; // Shocked energy;
        // -- Radiative losses
        double dErad2dR = p_pars->eps_rad * dEsh2dR;
        // -- Energy equation
        double dEint2dR = dEsh2dR + dEad2dR - dErad2dR + dEingdR_abs_dop;// dEingdR_abs_dop; // / (m_pars.M0 * c ** 2)
        double _x = dEint2dR/Eint2;

        double dtcomov_dR = 1.0 / beta / Gamma / CGS::c;
        double dttdr;
//        if (p_spread->m_method != LatSpread::METHODS::iNULL)
//            dttdr = 1. / (CGS::c * beta) * sqrt(1. + R * R * dthetadr * dthetadr) - (1. / CGS::c);
//        else
//            dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        dttdr = 1. / (CGS::c * Gamma * Gamma * beta * (1. + beta));
        // ****************************************
//        if (!std::isfinite(Gamma) ||
//        !std::isfinite(beta) //||
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
        if (mom < 0.){
            (*p_log)(LOG_ERR, AT) << " mom < 0 = "<< mom << " in kN RHS dynamics\n";
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
                       << " beta="<<beta//<< "\n"
                       << " theta"<<theta//<< "\n"
                       << " M2="<<M2//<< "\n"
                       << " Eint2="<<Eint2//<< "\n"
                       << " GammaRel"<<GammaRel//<< "\n"
                       << " GammaRho"<<GammaRho//<< "\n"
                       << " \n";
//            exit(1);
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
//            out_Y[i+Q_SOL::iR]      = dRdt;
//            out_Y[i+Q_SOL::iRsh]    = dRshdt;
//            out_Y[i+Q_SOL::itt]     = dRdt * dttdr;
//            out_Y[i+Q_SOL::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+Q_SOL::iGamma]  = 0 * dGammadR;
//            out_Y[i+Q_SOL::iEint2]  = 0 * dEint2dR;
//            out_Y[i+Q_SOL::itheta]  = 0 * dthetadr;
//            out_Y[i+Q_SOL::iErad2]  = 0 * dErad2dR;
//            out_Y[i+Q_SOL::iEsh2]   = 0 * dEsh2dR;
//            out_Y[i+Q_SOL::iEad2]   = 0 * dEad2dR;
//            out_Y[i+Q_SOL::iM2]     = 0 * dM2dR;
////            Gamma = p_pars->Gamma0;
//        }
//        else{
//            out_Y[i+Q_SOL::iR]      = dRdt;//1.0 / beta / CGS::c;
//            out_Y[i+Q_SOL::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
//            out_Y[i+Q_SOL::itt]     = dRdt * dttdr;
//            out_Y[i+Q_SOL::itcomov] = dRdt * dtcomov_dR;
//            out_Y[i+Q_SOL::iGamma]  = dRdt * dGammadR;
//            out_Y[i+Q_SOL::iEint2]  = dRdt * dEint2dR;
//            out_Y[i+Q_SOL::itheta]  = dRdt * dthetadr;
//            out_Y[i+Q_SOL::iErad2]  = dRdt * dErad2dR;
//            out_Y[i+Q_SOL::iEsh2]   = dRdt * dEsh2dR;
//            out_Y[i+Q_SOL::iEad2]   = dRdt * dEad2dR;
//            out_Y[i+Q_SOL::iM2]     = dRdt * dM2dR;
//        }
        out_Y[i+Q_SOL::iR]      = dRdt;//1.0 / beta / CGS::c;
        out_Y[i+Q_SOL::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
        out_Y[i+Q_SOL::itt]     = dRdt * dttdr;
        out_Y[i+Q_SOL::itcomov] = dRdt * dtcomov_dR;
        out_Y[i+Q_SOL::imom]    = dRdt * dmomdR; //dRdt * dGammadR / Gamma;
//        out_Y[i+Q_SOL::iGamma]  = 0.;//dRdt * dGammadR; //dRdt * dGammadR / Gamma;
        out_Y[i+Q_SOL::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+Q_SOL::iEinj]   = dRdt * dEingdR_abs;
        out_Y[i+Q_SOL::itheta]  = dRdt * dthetadr;
        out_Y[i+Q_SOL::iErad2]  = dRdt * dErad2dR;
        out_Y[i+Q_SOL::iEsh2]   = dRdt * dEsh2dR;
        out_Y[i+Q_SOL::iEad2]   = dRdt * dEad2dR;
        out_Y[i+Q_SOL::iM2]     = dRdt * dM2dR;
//        if (Gamma)
//        out_Y[i+Q_SOL::iR]      = dRdt;//1.0 / beta / CGS::c;
//        out_Y[i+Q_SOL::iRsh]    = dRshdt;//1.0 / beta / CGS::c;
//        out_Y[i+Q_SOL::itt]     = dRdt * dttdr;
//        out_Y[i+Q_SOL::itcomov] = dRdt * dtcomov_dR;
//        out_Y[i+Q_SOL::iGamma]  = dRdt * dGammadR;
//        out_Y[i+Q_SOL::iEint2]  = dRdt * dEint2dR;
//        out_Y[i+Q_SOL::itheta]  = dRdt * dthetadr;
//        out_Y[i+Q_SOL::iErad2]  = dRdt * dErad2dR;
//        out_Y[i+Q_SOL::iEsh2]   = dRdt * dEsh2dR;
//        out_Y[i+Q_SOL::iEad2]   = dRdt * dEad2dR;
//        out_Y[i+Q_SOL::iM2]     = dRdt * dM2dR;
        // ****************************************
//        if (beta<1e-5)  p_pars->end_evolution = true;
    }
    /// --- Evaluate the density (and its velocity) at point ej_R left by blast wave at a point j_R
    void evalDensAndItsVelocityBehindBlastWave(double & rho, double & GammaCMB, double & p_cbm,
                                               double rho_def, double rho_prev,
                                               double ej_R, double j_R, double j_rho, double j_rho2, double j_Gamma, double j_Gamma0, double j_P2){

        GammaCMB = 1.;
        /// exponential decay from the point of entry
        double rho0 = getVal(Q::irho,0);
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
    }

    /// Density profiles application based on the jet BW position and velocity
    void evalDensAndItsVelocityBehindBlastWave_Case1(
            double j_R, double j_Gamma, double j_ctheta, double j_rho, double j_rho2, double j_Gamma0, double j_P2,
            double ej_R, double ej_Gamma, double ej_theta ){

        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);

        if (ej_Gamma <= 1.) { ej_Gamma = 1. + 5e-5; } // TODO REMOVE

        size_t prev_ix = p_pars->comp_ix;

        // --- evaluate what density profile the ejecta blast wave would experience
        p_dens->evaluateRhoDrhoDrDefault(ej_R, ej_ctheta); // set default values for density
        double rho_prev = getVal(Q::irho, prev_ix-1);
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
        double r_prev = getVal(Q::iR, prev_ix - 1);
//        double rho_prev = getVal(Q::irho, evaled_ix - 1);
        if (ej_R == r_prev){ (*p_log)(LOG_ERR,AT) << AT << " R = R_i-1 \n"; exit(1); }
        m_data[Q::irho][prev_ix+1] = rho;
        m_data[Q::iR][prev_ix+1] = ej_R;
        m_data[Q::iPcbm][prev_ix + 1] = P;
        double dr = m_data[Q::iR][prev_ix] - m_data[Q::iR][prev_ix-1];
//            double dr1 = m_data[Q::iR][0] - m_data[Q::iR][1];
        drhodr = dydx(m_data[Q::iR], m_data[Q::irho], m_data[Q::idrhodr],
                      dr, prev_ix+1, ej_R, true);
//        double drhodr2 = (rho - rho_prev) / (ej_R - r_prev);
        if (rho < p_dens->m_rho_floor_val*rho_def){ rho = p_dens->m_rho_floor_val*rho_def; drhodr = 0.; }
        if (!std::isfinite(rho) || !std::isfinite(drhodr)) {
            (*p_log)(LOG_ERR,AT)  <<" bad rho="<<rho<<" or m_drhodr="<<drhodr<<" \n Exiting...\n";
//            std::cerr << AT <<" \n";
            exit(1);
        }

        /// Compute the CBM density derivative (with radius)
        double GammaCBM_prev = getVal(Q::iGammaCBM, prev_ix - 1);
        if (GammaCMB == GammaCBM_prev ){
            dGammaRhodR = 0;
        }
        else{
            m_data[Q::iGammaCBM][prev_ix+1] = GammaCMB;
            dGammaRhodR = dydx(m_data[Q::iR], m_data[Q::iGammaCBM], m_data[Q::idGammaCBMdr],
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

        /// compute the relative velocity of the ISM (with respect to ejecta) and its derivative
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

        double ej_Gamma_prev = getVal(Q::iGamma, prev_ix - 1);
        double GammaCMB_prev = getVal(Q::iGammaCBM, prev_ix - 1);
        double GammaREL_prev = EQS::GammaRel(ej_Gamma_prev, GammaCMB_prev);
        if ((ej_Gamma == ej_Gamma_prev) || (GammaREL_prev < 1.) || (GammaCMB_prev == 1.)){
            dGammaRELdGamma = 1.;
        }
        else{
            m_data[Q::iGammaREL][prev_ix+1] = GammaREL;
            double dGammaRel = m_data[Q::iGammaREL][prev_ix] - m_data[Q::iGammaREL][prev_ix-1];
            dGammaRELdGamma = dydx(m_data[Q::iGammaREL], m_data[Q::iGamma], m_data[Q::idGammaRELdGamma],
                                   dr, prev_ix+1, GammaREL, false);
            dPdrho = dydx(m_data[Q::iPcbm], m_data[Q::irho], m_data[Q::idPCBMdrho],
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
                                   void * _others, size_t evaled_ix){
//        auto * p_pars = (struct Pars *) params; // removing EATS_pars for simplicity
        auto * p_others = (std::vector<std::unique_ptr<BlastWaveBase>> *) _others;
        auto & others = * p_others;
        // --- current ejecta bw
//        double ej_Gamma  = Y[p_pars->ii_eq + DynRadBlastWave::Q_SOL::iGamma];
        double ej_mom = Y[p_pars->ii_eq + DynRadBlastWave::Q_SOL::imom];
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
        double ej_R      = Y[p_pars->ii_eq + DynRadBlastWave::Q_SOL::iR];
        double theta_b0  = p_pars->theta_b0;
        double ej_theta  = Y[p_pars->ii_eq + DynRadBlastWave::Q_SOL::itheta];
        double ej_ctheta = p_pars->ctheta0;//ctheta(ej_theta);
        // -- loop over jets

        int i_ej_l = p_pars->which_jet_layer_to_use;//others.size()-1;
        bool is_within = false;
//        for (size_t ij = 0; ij < others.size(); ij++){
//            double j_R      = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
//            double j_Gamma  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iGamma];
//            double j_theta  = Y[others[ij]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::itheta];
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
//            double j_R      = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
//            double j_Gamma  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iGamma];
//            double j_theta  = Y[others[0]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::itheta];
//            double j_ctheta = others[0]->ctheta(j_theta);
//            double j_theta0 = others[0]->theta0(j_theta); //
//            double j_theta1 = others[0]->theta1(j_theta);
//            if (ej_ctheta < j_theta0){
//                i_ej_l = 0.;
//                is_within = true;
//            }
//        }

//        i_ej_l = 0;//others.size()-1;
        double j_mom = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom];
        double j_Gamma = EQS::GamFromMom(j_mom);
//        double j_Gamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iGamma];
//        double ibeta = EQS::BetFromMom(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom]);
        double j_theta = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::itheta];
//        double j_lnGamma = Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::ilnGamma];
//        double j_Gamma = std::exp(j_lnGamma);
        double j_Gamma0 = others[i_ej_l]->getPars()->Gamma0;
//        double j_ctheta = others[i_ej_l]->ctheta(j_theta);
        double j_ctheta = LatStruct::ctheta(j_theta,
                                            others[i_ej_l]->getPars()->ilayer,
                                            others[i_ej_l]->getPars()->nlayers);//others[i_ej_l]->ctheta(j_theta);
        if (ej_ctheta < j_ctheta){ is_within = true; }

        auto & other = others[i_ej_l];
        double j_R = Y[other->getPars()->ii_eq + Q_SOL::iR];

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
                double other_Gamma = EQS::GamFromMom(Y[other_i + DynRadBlastWave::Q_SOL::imom]);//std::exp( Y[other_i + DynRadBlastWave::Q_SOL::ilnGamma] );
//                double other_Gamma = Y[other_i + DynRadBlastWave::Q_SOL::iGamma];//std::exp( Y[other_i + DynRadBlastWave::Q_SOL::ilnGamma] );
                // --- density experienced by the jet blast wave and the density in the downstream (rho2)
                other->getDensIsm()->evaluateRhoDrhoDrDefault(j_R, INFINITY);
                double j_rho = other->getDensIsm()->m_rho_def;
                double j_drhodr = other->getDensIsm()->m_drhodr_def;
                double j_adi = other->getEos()->getGammaAdi(other_Gamma, EQS::Beta(other_Gamma));
                double j_rho2 = EQS::rho2t(j_Gamma, j_adi, j_rho);
                double j_V2 = Y[other_i + DynRadBlastWave::Q_SOL::iM2] / j_rho2; // Units -> c^2 for energy
                double j_P2 = (j_adi - 1.) * Y[other_i + DynRadBlastWave::Q_SOL::iEint2] / j_V2 * CGS::c * CGS::c;//# * CGS::c * CGS::c; // Units -> c^2 for energy
                double cs = sqrt(j_adi * j_P2 / j_rho2) / CGS::c;
                double cs2 = sqrt(j_adi*(j_adi-1)*(j_Gamma-1)/(1+j_adi*(j_Gamma-1)));
                // ---
//                size_t i = p_pars->ii_eq;
//                double ej_Gamma = Y[i + DynRadBlastWave::Q_SOL::iGamma];
//                double ej_R = Y[i + DynRadBlastWave::Q_SOL::iR];
//                double theta_b0 = p_pars->theta_b0;
//                double ej_theta = Y[i + DynRadBlastWave::Q_SOL::itheta];
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
        double rho_prev    = getVal(Q::irho, evaled_ix-1);
        double drhodr_prev = getVal(Q::idrhodr, evaled_ix-1);
        double prev_ej_R   = getVal(Q::iR, evaled_ix-1);
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
//                      << " j_theta0="<<others[ijl]->theta0(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::itheta])
//                      << " j_theta1="<<others[ijl]->theta1(Y[others[i_ej_l]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::itheta])
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
    }
    void evaluateRhsDensModel2(double * out_Y, size_t i, double x, double const * Y, void * others, size_t evaled_ix) override {

        /// do not evaluate RHS if the evolution was terminated
        if (p_pars->end_evolution) {
            return;
        }

        /// compute density profile in front of the kN BW
        if (p_pars->use_dens_prof_behind_jet_for_ejecta){
            prepareDensProfileFromJet(out_Y,i,x,Y,others,evaled_ix);
        }
        else{
            double ej_Gamma  = EQS::GamFromMom( Y[i + DynRadBlastWave::Q_SOL::imom] );
//            double ej_Gamma  = Y[i + DynRadBlastWave::Q_SOL::iGamma];
            if (ej_Gamma < 1.) {
                (*p_log)(LOG_ERR,AT) << "Gamma < 1\n";
                ej_Gamma = 1. + 1e-5;
            }
            double ej_R      = Y[i + DynRadBlastWave::Q_SOL::iR];
            double theta_b0  = p_pars->theta_b0;
            double ej_theta  = Y[i + DynRadBlastWave::Q_SOL::itheta];
            double ej_ctheta = LatStruct::ctheta(ej_theta,p_pars->ilayer,p_pars->nlayers);//ctheta(ej_theta);
//            if (ej_Gamma>10.){
//                (*p_log)(LOG_ERR,AT)
//                        << "["<<p_pars->ishell<<", "<<p_pars->ilayer<<"]"
//                        << " ii_eq = "<<p_pars->ii_eq
//                        << " i = " << i
//                        << " Gamma="<<ej_Gamma
//                        << " R = " << ej_R
//                        << " theta_b0 = " <<theta_b0
//                        << " ctheta =" <<ej_ctheta
//                        << "\n"
//                        << " Exiting...\n";
////                std::cerr << AT << "\n";
//                exit(1);
//            }
            set_standard_ism(ej_R, ej_ctheta, ej_Gamma);
        }

        /// evaluate actual RHS
        evaluateRhsDens(out_Y, i, x, Y);

    }

    /// ---


};


#endif //SRC_BLASTWAVE_DYN_H
