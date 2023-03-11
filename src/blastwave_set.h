//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_SET_H
#define SRC_BLASTWAVE_SET_H

#include "pch.h"
#include "utils.h"
#include "base_equations.h"
#include "interpolators.h"
#include "ode_solvers.h"
#include "quadratures.h"
#include "rootfinders.h"
#include "observer.h"
#include "synchrotron_an.h"

#include "magnetar.h"
#include "blastwave_base.h"
#include "blastwave_rad.h"
#include "blastwave_dyn.h"

class BlastWaveCollision{
    struct FsolvePars{
        EOSadi * p_eos = nullptr;
        FsolvePars(EOSadi * eos_pointer) : p_eos(eos_pointer){ }
        double m_gam1{}, m_beta1{}, m_mass1{}, m_eint1{}, m_adi1{};
        double m_gam2{}, m_beta2{}, m_mass2{}, m_eint2{}, m_adi2{};
        int getTheMostMassive() const{ return m_mass1 > m_mass2 ? 1 : 2; }
        int getTheFastest() const{ return m_gam1 > m_gam2 ? 1 : 2; }
        int getTheMostEnergetic() const{ return m_eint1 > m_eint2 ? 1 : 2; }
        /// chose which shell to delete: the slowest/less energetic/less massive
        int choseShell(){
            int idx_massive = getTheMostMassive();
            int idx_energetic = getTheMostEnergetic();
            int idx_fastest = getTheFastest();
            if (idx_massive == idx_energetic && idx_massive == idx_fastest){ return idx_massive;}
            if (idx_massive != idx_energetic && idx_massive == idx_fastest){ return idx_massive;}
            if (idx_massive == idx_energetic && idx_massive != idx_fastest){ return idx_fastest;}
            if (idx_massive != idx_energetic && idx_massive != idx_fastest){ return idx_fastest;}
            std::cerr<<AT<<" ERROR\n";
            exit(1);
        }
        void set(double gam1, double beta1, double mass1, double eint1, double adi1,
                 double gam2, double beta2, double mass2, double eint2, double adi2){
            m_gam1=gam1; m_beta1=beta1; m_mass1=mass1; m_eint1=eint1; m_adi1=adi1;
            m_gam2=gam2; m_beta2=beta2; m_mass2=mass2; m_eint2=eint2; m_adi2=adi2;
        }
    };
    FsolvePars * p_colsolve;
    std::unique_ptr<logger> p_log;
    EOSadi * p_eos = nullptr;
//    std::unique_ptr<RadBlastWave> & p_bw1;
//    std::unique_ptr<RadBlastWave> & p_bw2;
public:
//    BlastWaveCollision(std::unique_ptr<RadBlastWave> & p_bw1,
//                       std::unique_ptr<RadBlastWave> & p_bw2,
//                       std::unique_ptr<logger> & p_log)
//        : p_bw1(p_bw1), p_bw2(p_bw2), p_log(p_log){
//
//    }
    BlastWaveCollision(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BlastWaveCollision");
        p_eos = new EOSadi();
        p_colsolve = new FsolvePars(p_eos);
    }
    ~BlastWaveCollision(){ delete p_colsolve; delete p_eos; }

    void collideBlastWaves(std::unique_ptr<RadBlastWave> & bw1,
                           std::unique_ptr<RadBlastWave> & bw2,
                           double * Y, double rcoll, size_t il){
        if (p_colsolve->p_eos == nullptr){
            (*p_log)(LOG_ERR,AT) << " eos pointer is not set\n;";
            exit(1);
        }
        if ((bw1->getPars()->end_evolution) || (bw2->getPars()->end_evolution)){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active\n";
            exit(1);
        }
        /// get relevant parameters
        double gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double m0_1 = bw1->getPars()->M0;
        double m2_1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        double m2_1_ = m2_1 * m0_1;
        double eint2_1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        double eint2_1_ = eint2_1 * bw1->getPars()->M0 * CGS::c * CGS::c;
        double adi1 = bw1->getEos()->getGammaAdi(gam1, beta1);
        /// --------------------
        double gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double m0_2 = bw2->getPars()->M0;
        double m2_2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        double m2_2_ = m2_2 * m0_2;
        double eint2_2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        double eint2_2_ = eint2_1 * bw2->getPars()->M0 * CGS::c * CGS::c;
        double adi2 = bw2->getEos()->getGammaAdi(gam2, beta2);
        /// -------------------
        p_colsolve->set(gam1, beta1, m2_1_ + m0_1, eint2_1_, adi1,
                        gam2, beta2, m2_2_ + m0_2, eint2_2_, adi2);
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll=" << il << "] Colliding shells: "
                              << " Masses=(" << p_colsolve->m_mass1 << ", " << p_colsolve->m_mass2 << ")"
                              << " Gammas=(" << p_colsolve->m_gam1 << ", " << p_colsolve->m_gam2 << ")"
                              << " betas=(" << p_colsolve->m_beta1 << ", " << p_colsolve->m_beta2 << ")"
                              << " Eint2=(" << p_colsolve->m_eint1 << ", " << p_colsolve->m_eint2 << ")"
                              << "\n";
        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1 * p_colsolve->m_mass1 + p_colsolve->m_gam2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = p_colsolve->m_eint1 + p_colsolve->m_eint2;
                //(p_colsolve->m_eint1 * p_colsolve->m_mass1 + p_colsolve->m_eint2 * p_colsolve->m_mass2)
                      /// (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        /// solve the collision
        int result = solve(i_gM, i_eM, i_mM);
        double m0_c = m0_2 + m0_1;
        double eint2_c = i_eM / (m0_c * CGS::c * CGS::c);
        double m2_c = (1 * (m2_1 + 1) / (m0_2/m0_1 + 1)) + (1 * (m2_2 + 1) / (1 + m0_1/m0_2));
        if (m2_c==1)
            m2_c = m2_1 + m2_2;
        double mom_c = i_gM * EQS::Beta(i_gM);
        /// chose which shell to update/delete
        int ish = p_colsolve->choseShell();
        double eint_before,eint_after,m2_before,m2_after,m0_before,m0_after;
        if (ish == 1){
            bw2->getPars()->end_evolution = true;
            bw1->getPars()->M0 = m0_c;
            Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] = mom_c;
            Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2] = eint2_c;
            Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] = m2_c;
        }
        else {
            bw1->getPars()->end_evolution = true;
            bw2->getPars()->M0 = m0_c;
            Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] = mom_c;
            Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2] = eint2_c;
            Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] = m2_c;
        }
        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<il<<"] Outcome for"
                              << " Eint2/M0c^2: ["<<eint2_1<<", "<<eint2_2<<"] -> "<<eint2_c
                              << " M2/M0: ["<<m2_1<<", "<<m2_2<<"] -> "<<m2_c
                              << " M0: ["<<m0_1<<", "<<m0_2<<"] -> "<<m0_c
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }


#if 0
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double m2_1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] * bw1->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2]
                <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        double m0_1 = bw1->getPars()->M0;
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        p_colsolve->m_mass1 = m0_1 + m2_1;

        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        double m2_2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] * bw2->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2]
                  <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi2 = bw2->getEos()->getGammaAdi(p_colsolve->m_gam2, p_colsolve->m_beta2);
        /// apply units for the second shell (mass and energy)
        double m0_2 = bw2->getPars()->M0;
        p_colsolve->m_eint2 = p_colsolve->m_eint2 * (bw2->getPars()->M0 * CGS::c * CGS::c);
        p_colsolve->m_mass2 = m2_2 + m0_2;

        /// log the data
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll=" << il << "] Colliding shells: "
//                              << "idx1="<<idx1<<", idx2="<<idx2
                              << " Masses=(" << p_colsolve->m_mass1 << ", " << p_colsolve->m_mass2 << ")"
                              << " Gammas=(" << p_colsolve->m_gam1 << ", " << p_colsolve->m_gam2 << ")"
                              << " betas=(" << p_colsolve->m_beta1 << ", " << p_colsolve->m_beta2 << ")"
                              << " Eint2=(" << p_colsolve->m_eint1 << ", " << p_colsolve->m_eint2 << ")"
                              << "\n";

        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1 * p_colsolve->m_mass1 + p_colsolve->m_gam2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = (p_colsolve->m_eint1 * p_colsolve->m_mass1 + p_colsolve->m_eint2 * p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1 + p_colsolve->m_mass2);

        /// solve the collision
        int result = solve(i_gM, i_eM);
//        i_mM =

        /// chose which shell to update/delete
        int ish = p_colsolve->choseShell();
        size_t iieq; size_t iieq_other;
        double eint_before=0., eint_after=0.;
        double m2_before=0., m2_after=0.;
        double m0_before=0., m0_after=0.;
        if (ish == 1){
            iieq = bw1->getPars()->ii_eq;
            iieq_other = bw2->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            m0_before = bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            //
            bw1->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM / bw1->getPars()->M0;//Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;
            eint_after = i_eM / ( bw1->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::Q_SOL::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw1->getPars()->M0;
        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            m0_before = bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM / bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;
            eint_after = i_eM / ( bw2->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::Q_SOL::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw2->getPars()->M0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
//        double _mom = i_gM * EQS::Beta(i_gM);
//        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
//        Y[iieq + DynRadBlastWave::Q_SOL::imom] = _mom;
//        Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = _eint2; // use ODE units

        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<il<<"] Outcome for"
//                              << " idx1="<<idx1<<", idx2="<<idx2 << " collision:"
                              << " Eint2/M0c^2: "<<eint_before<<" -> "<<eint_after
                              << " M2/M0: "<<m2_before<<" -> "<<m2_after
                              << " M0: "<<m0_before<<" -> "<<m0_after
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }
#endif
    }
private:
    static void func( int n, double x[], double fx[], int &iflag, void * pars ){
        auto pp = (struct FsolvePars *) pars;
        double gM = x[0];
        double eM = x[1];
//        double em = x[1];
        double gAdiM = pp->p_eos->getGammaAdi(gM, EQS::Beta(gM));
        /// total mass conservation
        double mM = pp->m_mass1 + pp->m_mass2;
        /// total energy consercation (Ek + Eint)
        double fx1 = (pp->m_gam1 * pp->m_mass1 + EQS::get_GammaEff(pp->m_gam1, pp->m_adi1) * pp->m_eint1)
                     + (pp->m_gam2 * pp->m_mass2 + EQS::get_GammaEff(pp->m_gam2, pp->m_adi2) * pp->m_eint2)
                     - (gM * mM + EQS::get_GammaEff(gM,gAdiM) * eM );
        /// total momentum conservation
        double fx2 = sqrt(pp->m_gam1 * pp->m_gam1 - 1) * (pp->m_mass1 + pp->m_adi1 * pp->m_eint1 )
                     + sqrt(pp->m_gam2 * pp->m_gam2 - 1) * (pp->m_mass2 + pp->m_adi2 * pp->m_eint2 )
                     - sqrt(gM * gM - 1) * (mM + gAdiM * eM );
//        double fx3 = em
//                     - pp->m_mass1 + pp->m_mass2;
        if (!std::isfinite(fx1) || !std::isfinite(fx2))
            iflag = -1;

        if (!std::isfinite(sqrt(pp->m_gam1 * pp->m_gam1 - 1)))
            iflag = -1;

        fx[0] = fx1;//std::log10(fx1);// * fx1; TODO THis is wrong. It should not be sqred!
        fx[1] = fx2;//std::log10(fx2);// * fx2;
//        fx[2] = fx3;
    };
    int solve(double & iGamma, double & iEint, double &im){
        double *fx;
        int iflag;
        int info;
        int lwa;
        int n = 2;
        double tol = 1e-8;
        double *wa;
        double *x;

        lwa = ( n * ( 3 * n + 13 ) ) / 2;
        fx = new double[n];
        wa = new double[lwa];
        x = new double[n];

//            (*p_log)(LOG_INFO,AT) << "\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve_test2():\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve() solves a nonlinear system of 2 equations.\n";

        x[0] = iGamma;
        x[1] = iEint;
//        x[2] = im;
        iflag = 0;
        func ( n, x, fx, iflag, p_colsolve );
//            r8vec2_print ( n, x, fx, "  Initial X, F(X)" );
        info = fsolve ( func, n, x, p_colsolve, fx, tol, wa, lwa );

        if (info<0){
            (*p_log)(LOG_ERR,AT)<< "\tFsolve failed (try setting 'tol' lower). "
                                   "Using initial guess values. New shell has "
                                <<"Gamma="<<iGamma<<" beta="<<EQS::Beta(iGamma)<<" Eint="<<iEint<<" im="<<im<<"\n";
//                i_gM = x[0];
//                i_eM = x[1];
        }
        else{
            (*p_log)(LOG_INFO,AT)<< "Fsolve successful. New shell has "
                                 <<"Gamma="<<x[0]
                                 <<" beta="<<EQS::Beta(x[0])
                                 <<" Eint="<<x[1]<<"\n";
            iGamma = x[0];
            iEint = x[1];
//            im = x[2];
        }
//            std::cout << "\n";
//            std::cout << "  Returned value of INFO = " << info << "\n";
//            r8vec2_print ( n, x, fx, "  Final X, FX" );
        delete [] fx;
        delete [] wa;
        delete [] x;
        return info;
    };
};


/*
 * Class to hold blastwaves with the same angular coordinate (shells of the same layer)
 * Here the velocity distribution, relative position and photosphere can be found.
 */
class CumulativeShell{

//    enum Q { ikapp};
    std::unique_ptr<logger> p_log;
    std::unique_ptr<BlastWaveCollision> p_coll;
    std::vector<std::unique_ptr<RadBlastWave>> p_bws_ej;
    std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws_ej;
    std::vector<std::vector<size_t>> relative_position;
    std::vector<size_t> m_idxs;
    std::vector<size_t> m_idxs0;
    std::vector<size_t> idx_tau_eq1;
    Vector m_rho;
    Vector m_radii_init;
    Vector m_radii_sorted;
    Vector m_dtau;
    Vector m_dtau_cum;
    VecVector m_shell_data;
    size_t m_tarr_size;
    const Array m_t_grid;
    bool do_collision;
public:
    size_t m_nshells;
    size_t n_active_shells;
    size_t mylayer{};
    CumulativeShell( Array t_grid, size_t nshells, int ilayer, int loglevel )
        : m_t_grid(t_grid) {
        m_nshells = nshells;
        n_active_shells = nshells; // as shells collide this number will decrease

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CumulativeShell");
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws_ej.emplace_back( std::make_unique<DynRadBlastWave>(t_grid, ishell, ilayer, loglevel ) );

        idx_tau_eq1.resize(t_grid.size(), 0);
        relative_position.resize( nshells );
        for (auto & arr : relative_position)
            arr.resize( t_grid.size() );

        m_rho.resize(nshells, 0.0);
        m_radii_init.resize(nshells, std::numeric_limits<double>::max());
        m_radii_sorted.resize(nshells, std::numeric_limits<double>::max());
        m_dtau.resize(nshells, 0.0);
        m_dtau_cum.resize(nshells, 0.0);
//        m_shell_data.resize(m_nshells);
//        for (auto arr:)
        m_idxs.resize(nshells, nshells-1); // fill with last value (if not active)
        /// initial shell indexes
        m_idxs0.resize(nshells, nshells-1); // fill with last value (if not active)
        for (size_t i = 0; i<nshells; i++)
            m_idxs0[i] = i;
        mylayer = ilayer;
        m_tarr_size = t_grid.size();
//        p_colsolve = std::make_unique<FsolvePars>(p_bws_ej[0]->getEos());
//        p_colsolve = new FsolvePars(p_bws_ej[0]->getEos());
//        p_colsolve->p_eos;
        p_coll = std::make_unique<BlastWaveCollision>(loglevel);


    }
//    ~CumulativeShell(){ delete p_colsolve; }
    Vector & getSortedRadii(){return m_radii_init;}
//    std::vector<size_t> & getIdx(){return m_idxs;}
    inline double getR(size_t i){return m_radii_init[m_idxs[i]];}
    std::unique_ptr<RadBlastWave> & getBW(size_t ish){
        if (p_bws_ej.empty()){
            (*p_log)(LOG_ERR, AT) << " shell does not contain blast waves\n";
            exit(1);
        }
        if (ish > p_bws_ej.size()){
            (*p_log)(LOG_ERR, AT) << "invalid memory accessed\n"; exit(1);
        }
        return p_bws_ej[ish];
    }
    inline size_t nBWs() const { return p_bws_ej.size();}
    inline std::vector<std::unique_ptr<RadBlastWave>> & getBWs() { return p_bws_ej; }
    /// update the number of active shells
    void updateActiveShells(){
        std::fill(m_idxs.begin(), m_idxs.end(),std::numeric_limits<size_t>::max());
        n_active_shells = m_nshells;
        size_t idx = 0;
        for (size_t i=0; i<n_active_shells; i++) {
            /// loop over all blastwaves and collect the active ones
            if (p_bws_ej[i]->getPars()->end_evolution) {
                continue; }
            /// if shell is active record its current index (order)
            m_idxs[idx] = i; // add only active shells to the list
            idx++;
        }
        n_active_shells = idx;
    }
    /// check if active shells are ordered according to their radius
    bool checkIfActiveShellsOrdered(const double * Y, size_t & idx0, size_t & idx1 ){
        idx0 = std::numeric_limits<size_t>::max();
        idx1 = std::numeric_limits<size_t>::max();
        std::fill(m_radii_init.begin(), m_radii_init.end(),std::numeric_limits<double>::max());
        bool is_sorted = true;
        for (size_t i=0; i<n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            m_radii_init[i] = Y[ p_bws_ej[idx]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ];
        }
        for (size_t i=1; i<n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            double rim1 = m_radii_init[i-1];
            /// if the next shells in the list has a radius > than the previous; shells not ordered
            if (m_radii_init[i] > rim1){
                rim1 = m_radii_init[i];
            }
            else {
                is_sorted = false;
//                (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<mylayer<<"] active shells are not ordered at "
//                                      << " (r[i] - r[i-1]) = "<< m_radii_init[i] - m_radii_init[i-1]
//                                      << " is negative (overrun)"
//                                      << " \n";
                /// store where the overrun occured
                if (idx0 == std::numeric_limits<size_t>::max()){
                    idx0 = m_idxs[i-1];
                    idx1 = m_idxs[i];
                }
                break;
            }
        }
        return is_sorted;
    }
    ///
    void collide(size_t idx1, size_t idx2, double * Y, double rcoll){
        if (idx1 == idx2){
            (*p_log)(LOG_ERR,AT) << " shells staged for collision are the same ish=idx2=" << idx1 << "\n";
            exit(1);
        }
        auto & bw1 = p_bws_ej[idx1];
        auto & bw2 = p_bws_ej[idx2];
        p_coll->collideBlastWaves(bw1, bw2, Y, rcoll, mylayer);
    }

#if 0
    /// check if shells are ordered radially
    bool checkIfShellsOrderedSetIndexes( const double * Y, size_t & idx0, size_t & idx1 ){
        /// clear arrays to avoid keeping inacitve shells here
        std::fill(m_radii_init.begin(), m_radii_init.end(),std::numeric_limits<double>::max());
        std::fill(m_idxs.begin(), m_idxs.end(),std::numeric_limits<size_t>::max());
        /// collect radii of the blast wave
        double r = -1.;
        std::vector<size_t> idxs_to_collide;
        size_t idx = 0; double ri = 0; bool is_sorted = true;
        idx0 = std::numeric_limits<size_t>::max();
        idx1 = std::numeric_limits<size_t>::max();
        for (size_t i=0; i<n_active_shells; i++) {
            /// if shell shell is deactivated; skip
            if (p_bws_ej[i]->getPars()->end_evolution) {
                continue;
            }
            /// if shell is active record its current index (order)
            m_idxs[idx] = i; // add only active shells to the list
//            std::cerr << "i="<<i<<" R="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ]<<"\n";
            m_radii_init[idx] = Y[p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ];
            /// if the next shells in the list has a radius > than the previous; shells not ordered
            if (m_radii_init[idx] > ri){
                ri = m_radii_init[idx];
            }
            else {
                is_sorted = false;
                (*p_log)(LOG_INFO,AT) << "Active shells are not ordered at "
                                   << " ilayer=" << mylayer
                                   << " i[idx-1]=" << idx - 1 << " r[idx-1]=" << m_radii_init[idx-1]
                                   << " i[idx]=" << idx << " r[idx]="<<m_radii_init[idx]
                                   << " \n";
                /// store where the overrun occured
                if (idx0 == std::numeric_limits<size_t>::max()){
                    idx0 = idx-1;
                    idx1 = idx;
                }
            }
            idx++;
        }
        n_active_shells = idx;
//        bool is_sorted = std::is_sorted(m_radii_init.begin(),m_radii_init.end());
//        if (!is_sorted){
//            std::cout << m_radii_init << "\n";
//        }
        return is_sorted;
    }

    /// set indexes 'm_idxs' of shells based on their radii
//    void oderShells( const double * Y ){
//        areShellsOrderedRadially(Y);
//        sort_indexes(m_radii_init, m_idxs );
//    }
#endif


    /// get indexes of shells that will collide next
    void evalWhichShellsCollideNext(size_t & ish1, size_t & ish2, double & tcoll, double & rcoll,
                                    double tim1, double ti, const double * Ym1_, const double * Y_){
        Vector _tcoll{};
        Vector _rcoll{};
        std::vector<size_t> _idx1s{};
        std::vector<size_t> _idx2s{};
        tcoll = std::numeric_limits<double>::max();
        if (tim1 > ti){
            (*p_log)(LOG_ERR,AT)<<"tim1="<<tim1<<" is larger than ti="<<ti<<"\n";
            exit(1);
        }
        (*p_log)(LOG_INFO, AT)
                << "\tLayer [il="<<mylayer<<"] Checking which shells will collide in"
                <<" between tim1="<<tim1<<" and ti="<<ti<<"\n";
        size_t n_collisions = 0;
        for (size_t ii=0; ii<n_active_shells; ++ii) {
            n_collisions = 0;
            size_t i_idx = m_idxs[ii];
            double r_i = m_radii_init[ii];
            if (r_i == std::numeric_limits<double>::max()){
                continue;
            }
            std::vector<size_t> overrun_shells;
            Vector t_at_which_overrun;
            /// now check which blast-waves has this one overrun

            for (size_t jj = ii; jj < n_active_shells; ++jj) {
                size_t j_idx = m_idxs[jj];
                double r_j = m_radii_init[jj];
                if (r_i > r_j) {
                    n_collisions += 1;
                    overrun_shells.push_back(j_idx);
                    /// interpolate the time at which shells collided
//                    double tim1 = m_t_grid[ix - 1];
//                    double ti = m_t_grid[ix];
                    double dt = ti - tim1;
                    double r0im1 = Ym1_[p_bws_ej[i_idx]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];//p_bws_ej[i]->getVal(BlastWaveBase::iR, (int) (ix - 1));
                    double r0i   = Y_[p_bws_ej[i_idx]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
                    double r1im1 = Ym1_[p_bws_ej[j_idx]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];//p_bws_ej[j]->getVal(BlastWaveBase::iR, (int) (ix - 1));
                    double r1i   = Y_[p_bws_ej[j_idx]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
                    double slope0 = (r0i - r0im1) / dt;
                    double slope1 = (r1i - r1im1) / dt;
                    double r0int = r0im1 - slope0 * tim1;
                    double r1int = r1im1 - slope1 * tim1;
                    double tc = (r1int - r0int) / (slope0 - slope1);
                    double rc = slope1 * tc + r0int;
                    if (!std::isfinite(tc)){
                        (*p_log)(LOG_ERR, AT)<< " nan in tcoll evaluation tc=" << tc << "\n";
                        exit(1);
                    }
                    if ((tc < tim1) or (tc > ti)) {
                        (*p_log)(LOG_ERR, AT)<< " inferred tcoll=" << tc
                                              << " is not in tim1-ti range [" << tim1 << ", " << ti << "]"
                                              << " ish1="<<i_idx<<" ish2="<<j_idx
                                              <<" len(overruns)="<<overrun_shells.size()<<"\n";
                        exit(1);
                    }

                    /// log the result
//                    (*p_log)(LOG_WARN, AT)
//                            << " Interpolating shell collision time" << " [ilayer=" << mylayer
//                            << " ishell=" << i << "] "
//                            << " mom=" << Y[p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom]
//                            << " (mom0=" << p_bws_ej[i]->getPars()->mom0 << ") "
//                            << " has overrun shell j=" << overrun_shells[0]
//                            << " with momenta "
//                            << Y[p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom]
//                            << " (mom0=" << p_bws_ej[overrun_shells[0]]->getPars()->mom0 << ") "
//                            << "\n";
//                    if (x < tcoll)
                    /// record the collision shells
//                    ish1 = i_idx;
//                    ish2 = j_idx;
//                    tcoll = tc;
                    _idx1s.push_back(i_idx);
                    _idx2s.push_back(j_idx);
                    _tcoll.push_back(tc);
                    _rcoll.push_back(rc);
//                    dt_from_ix_to_coll = tc - m_t_grid[ix];
//                    return;
                }
            }
        }
        size_t idx_min = indexOfMinimumElement(_tcoll);
        tcoll = _tcoll[idx_min];
        rcoll = _rcoll[idx_min];
        ish1 = _idx1s[idx_min];
        ish2 = _idx2s[idx_min];
        (*p_log)(LOG_INFO, AT) << "\tLayer [il="<<mylayer<<"] interpolated tcoll="<<tcoll
                                         <<" (idx1="<<ish1<<" idx2="<<ish2<<")\n";
        if ((tcoll == std::numeric_limits<double>::max())||(!std::isfinite(tcoll))){
            (*p_log)(LOG_ERR,AT)<<" Failed to find tcoll in layer="<<mylayer<<"\n";
            exit(1);
        }
    }

#if 0
    /// solve energy, momentum and mass conservation to get properties of the shell after collision
    void collideBlastWaves_old(size_t idx1, size_t idx2, double * Y, double rcoll){
        if (idx1 == idx2){
            (*p_log)(LOG_ERR,AT) << " shells staged for collision are the same ish=idx2=" << idx1 << "\n";
            exit(1);
        }
        if (p_colsolve->p_eos == nullptr){
            (*p_log)(LOG_ERR,AT) << " eos pointer is not set\n;";
            exit(1);
        }
        auto & bw1 = p_bws_ej[idx1];
        auto & bw2 = p_bws_ej[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                    "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi2 = bw2->getEos()->getGammaAdi(p_colsolve->m_gam2, p_colsolve->m_beta2);
        /// apply units for the second shell (mass and energy)
        p_colsolve->m_mass2 = bw2->getPars()->M0 + (p_colsolve->m_mass2 * bw2->getPars()->M0);
        p_colsolve->m_eint2 = p_colsolve->m_eint2 * (bw2->getPars()->M0 * CGS::c * CGS::c);
        /// log the data
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll="<<mylayer<<"] Colliding shells: "
                              << "idx1="<<idx1<<", idx2="<<idx2
                              << " Masses=("<<p_colsolve->m_mass1<<", "<<p_colsolve->m_mass2<<")"
                              << " Gammas=("<<p_colsolve->m_gam1<<", "<<p_colsolve->m_gam2<<")"
                              << " betas=("<<p_colsolve->m_beta1<<", "<<p_colsolve->m_beta2<<")"
                              << " Eint2=("<<p_colsolve->m_eint1<<", "<<p_colsolve->m_eint2<<")"
                              << "\n";
        /// solve equation for the shell collision
        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1*p_colsolve->m_mass1 + p_colsolve->m_gam2*p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = (p_colsolve->m_eint1*p_colsolve->m_mass1+p_colsolve->m_eint2*p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        int x_ = 1;
        // define a system of non-linear equations to solve
        auto func = []( int n, double x[], double fx[],
                        int &iflag, void * pars ){
            auto pp = (struct FsolvePars *) pars;
            double gM = x[0];
            double eM = x[1];
            double gAdiM = pp->p_eos->getGammaAdi(gM, EQS::Beta(gM));
            /// total mass conservation
            double mM = pp->m_mass1 + pp->m_mass2;
            /// total energy consercation (Ek + Eint)
            double fx1 = (pp->m_gam1 * pp->m_mass1 + EQS::get_GammaEff(pp->m_gam1, pp->m_adi1) * pp->m_eint1)
                         + (pp->m_gam2 * pp->m_mass2 + EQS::get_GammaEff(pp->m_gam2, pp->m_adi2) * pp->m_eint2)
                         - (gM * mM + EQS::get_GammaEff(gM,gAdiM) * eM );
            /// total momentum conservation
            double fx2 = sqrt(pp->m_gam1 * pp->m_gam1 - 1) * (pp->m_mass1 + pp->m_adi1 * pp->m_eint1 )
                         + sqrt(pp->m_gam2 * pp->m_gam2 - 1) * (pp->m_mass2 + pp->m_adi2 * pp->m_eint2 )
                         - sqrt(gM * gM - 1) * (mM + gAdiM * eM );
            if (!std::isfinite(fx1) || !std::isfinite(fx2))
                iflag = -1;

            if (!std::isfinite(sqrt(pp->m_gam1 * pp->m_gam1 - 1)))
                iflag = -1;

            fx[0] = fx1;//std::log10(fx1);// * fx1; TODO THis is wrong. It should not be sqred!
            fx[1] = fx2;//std::log10(fx2);// * fx2;
        };
        using func_ptr = void(*)(int, double *, double *, int&, void*);
        func_ptr p_func = func;

        /// solve the system using a 'fsolve'
        auto solve = [&]( ){
            double *fx;
            int iflag;
            int info;
            int lwa;
            int n = 2;
            double tol = 1e-3;
            double *wa;
            double *x;

            lwa = ( n * ( 3 * n + 13 ) ) / 2;
            fx = new double[n];
            wa = new double[lwa];
            x = new double[n];

//            (*p_log)(LOG_INFO,AT) << "\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve_test2():\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve() solves a nonlinear system of 2 equations.\n";

            x[0] = i_gM;
            x[1] = i_eM;
            iflag = 0;
            func ( n, x, fx, iflag, p_colsolve );
//            r8vec2_print ( n, x, fx, "  Initial X, F(X)" );
            info = fsolve ( func, n, x, p_colsolve, fx, tol, wa, lwa );

            if (info<0){
                (*p_log)(LOG_ERR,AT)<< "Fsolve failed (try setting 'tol' lower). "
                                       "Using initial guess values. New shell has "
                                    <<"Gamma="<<i_gM<<" beta="<<EQS::Beta(i_gM)<<" Eint="<<i_eM<<"\n";
//                i_gM = x[0];
//                i_eM = x[1];
            }
            else{
                (*p_log)(LOG_INFO,AT)<< "Fsolve successful. New shell has "
                                     <<"Gamma="<<x[0]
                                     <<" beta="<<EQS::Beta(x[0])
                                     <<" Eint="<<x[1]<<"\n";
                i_gM = x[0];
                i_eM = x[1];
            }
//            std::cout << "\n";
//            std::cout << "  Returned value of INFO = " << info << "\n";
//            r8vec2_print ( n, x, fx, "  Final X, FX" );
            delete [] fx;
            delete [] wa;
            delete [] x;
        };
        solve();
        /// chose which shell to terminate (the slower one)
        int ish = p_colsolve->choseShell();
//        ish = 2;
        size_t iieq; size_t iieq_other;
        double eint_before, eint_after;
        double m2_before, m2_after;
        if (ish == 1){
            iieq = bw1->getPars()->ii_eq;
            iieq_other = bw2->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::Q_SOL::imom] = _mom;
        Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = _eint2; // use ODE units

        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<mylayer<<"] Outcome for"
                              << " idx1="<<idx1<<", idx2="<<idx2 << " collision:"
                              << " Eint2/M0c^2: "<<eint_before<<" -> "<<_eint2
                              << " M2/M0: "<<m2_before<<" -> "<<m2_after
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }
//


//
//        Y[iieq + DynRadBlastWave::Q_SOL::i] = i_gM * EQS::Beta(i_gM);


    }
    void collideBlastWaves(size_t idx1, size_t idx2, double * Y, double rcoll){
        if (idx1 == idx2){
            (*p_log)(LOG_ERR,AT) << " shells staged for collision are the same ish=idx2=" << idx1 << "\n";
            exit(1);
        }
        if (p_colsolve->p_eos == nullptr){
            (*p_log)(LOG_ERR,AT) << " eos pointer is not set\n;";
            exit(1);
        }
        auto & bw1 = p_bws_ej[idx1];
        auto & bw2 = p_bws_ej[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                    "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi2 = bw2->getEos()->getGammaAdi(p_colsolve->m_gam2, p_colsolve->m_beta2);
        /// apply units for the second shell (mass and energy)
        p_colsolve->m_mass2 = bw2->getPars()->M0 + (p_colsolve->m_mass2 * bw2->getPars()->M0);
        p_colsolve->m_eint2 = p_colsolve->m_eint2 * (bw2->getPars()->M0 * CGS::c * CGS::c);
        /// log the data
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll="<<mylayer<<"] Colliding shells: "
                              << "idx1="<<idx1<<", idx2="<<idx2
                              << " Masses=("<<p_colsolve->m_mass1<<", "<<p_colsolve->m_mass2<<")"
                              << " Gammas=("<<p_colsolve->m_gam1<<", "<<p_colsolve->m_gam2<<")"
                              << " betas=("<<p_colsolve->m_beta1<<", "<<p_colsolve->m_beta2<<")"
                              << " Eint2=("<<p_colsolve->m_eint1<<", "<<p_colsolve->m_eint2<<")"
                              << "\n";
        /// solve equation for the shell collision
        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1*p_colsolve->m_mass1 + p_colsolve->m_gam2*p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = (p_colsolve->m_eint1*p_colsolve->m_mass1+p_colsolve->m_eint2*p_colsolve->m_mass2)
                      / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        int x_ = 1;
        // define a system of non-linear equations to solve
        auto func = []( int n, double x[], double fx[],
                        int &iflag, void * pars ){
            auto pp = (struct FsolvePars *) pars;
            double gM = x[0];
            double eM = x[1];
            double gAdiM = pp->p_eos->getGammaAdi(gM, EQS::Beta(gM));
            /// total mass conservation
            double mM = pp->m_mass1 + pp->m_mass2;
            /// total energy consercation (Ek + Eint)
            double fx1 = (pp->m_gam1 * pp->m_mass1 + EQS::get_GammaEff(pp->m_gam1, pp->m_adi1) * pp->m_eint1)
                         + (pp->m_gam2 * pp->m_mass2 + EQS::get_GammaEff(pp->m_gam2, pp->m_adi2) * pp->m_eint2)
                         - (gM * mM + EQS::get_GammaEff(gM,gAdiM) * eM );
            /// total momentum conservation
            double fx2 = sqrt(pp->m_gam1 * pp->m_gam1 - 1) * (pp->m_mass1 + pp->m_adi1 * pp->m_eint1 )
                         + sqrt(pp->m_gam2 * pp->m_gam2 - 1) * (pp->m_mass2 + pp->m_adi2 * pp->m_eint2 )
                         - sqrt(gM * gM - 1) * (mM + gAdiM * eM );
            if (!std::isfinite(fx1) || !std::isfinite(fx2))
                iflag = -1;

            if (!std::isfinite(sqrt(pp->m_gam1 * pp->m_gam1 - 1)))
                iflag = -1;

            fx[0] = fx1;//std::log10(fx1);// * fx1; TODO THis is wrong. It should not be sqred!
            fx[1] = fx2;//std::log10(fx2);// * fx2;
        };
        using func_ptr = void(*)(int, double *, double *, int&, void*);
        func_ptr p_func = func;

        /// solve the system using a 'fsolve'
        auto solve = [&]( ){
            double *fx;
            int iflag;
            int info;
            int lwa;
            int n = 2;
            double tol = 1e-3;
            double *wa;
            double *x;

            lwa = ( n * ( 3 * n + 13 ) ) / 2;
            fx = new double[n];
            wa = new double[lwa];
            x = new double[n];

//            (*p_log)(LOG_INFO,AT) << "\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve_test2():\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve() solves a nonlinear system of 2 equations.\n";

            x[0] = i_gM;
            x[1] = i_eM;
            iflag = 0;
            func ( n, x, fx, iflag, p_colsolve );
//            r8vec2_print ( n, x, fx, "  Initial X, F(X)" );
            info = fsolve ( func, n, x, p_colsolve, fx, tol, wa, lwa );

            if (info<0){
                (*p_log)(LOG_ERR,AT)<< "Fsolve failed (try setting 'tol' lower). "
                                       "Using initial guess values. New shell has "
                                    <<"Gamma="<<i_gM<<" beta="<<EQS::Beta(i_gM)<<" Eint="<<i_eM<<"\n";
//                i_gM = x[0];
//                i_eM = x[1];
            }
            else{
                (*p_log)(LOG_INFO,AT)<< "Fsolve successful. New shell has "
                                     <<"Gamma="<<x[0]
                                     <<" beta="<<EQS::Beta(x[0])
                                     <<" Eint="<<x[1]<<"\n";
                i_gM = x[0];
                i_eM = x[1];
            }
//            std::cout << "\n";
//            std::cout << "  Returned value of INFO = " << info << "\n";
//            r8vec2_print ( n, x, fx, "  Final X, FX" );
            delete [] fx;
            delete [] wa;
            delete [] x;
        };
        solve();
        /// chose which shell to terminate (the slower one)
        int ish = p_colsolve->choseShell();
//        ish = 2;
        size_t iieq; size_t iieq_other;
        double eint_before, eint_after;
        double m2_before, m2_after;
        if (ish == 1){
            iieq = bw1->getPars()->ii_eq;
            iieq_other = bw2->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::Q_SOL::imom] = _mom;
        Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = _eint2; // use ODE units

        (*p_log)(LOG_INFO,AT) << "\tLayer [il="<<mylayer<<"] Outcome for"
                              << " idx1="<<idx1<<", idx2="<<idx2 << " collision:"
                              << " Eint2/M0c^2: "<<eint_before<<" -> "<<_eint2
                              << " M2/M0: "<<m2_before<<" -> "<<m2_after
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }
//


//
//        Y[iieq + DynRadBlastWave::Q_SOL::i] = i_gM * EQS::Beta(i_gM);


    }
#endif

#if 0
    void collideBlastWaves(size_t idx1, size_t idx2, double * Y, double rcoll){
        if (idx1 == idx2){
            (*p_log)(LOG_ERR,AT) << " shells staged for collision are the same ish=idx2=" << idx1 << "\n";
            exit(1);
        }
        if (p_colsolve->p_eos == nullptr){
            (*p_log)(LOG_ERR,AT) << " eos pointer is not set\n;";
            exit(1);
        }
        auto & bw1 = p_bws_ej[idx1];
        auto & bw2 = p_bws_ej[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                   "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
//        double m_gam1, m_beta1, m_mass1, m_eint1, m_adi1;
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
//        double m_gam2, m_beta2, m_mass2, m_eint2, m_adi2;
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
        p_colsolve->m_adi2 = bw2->getEos()->getGammaAdi(p_colsolve->m_gam2, p_colsolve->m_beta2);
        /// apply units for the second shell (mass and energy)
        p_colsolve->m_mass2 = bw2->getPars()->M0 + (p_colsolve->m_mass2 * bw2->getPars()->M0);
        p_colsolve->m_eint2 = p_colsolve->m_eint2 * (bw2->getPars()->M0 * CGS::c * CGS::c);
        /// log the data
        (*p_log)(LOG_INFO,AT) << "\tLayer [ll="<<mylayer<<"] Colliding shells: "
            << "idx1="<<idx1<<", idx2="<<idx2
            << " Masses=("<<p_colsolve->m_mass1<<", "<<p_colsolve->m_mass2<<")"
            << " Gammas=("<<p_colsolve->m_gam1<<", "<<p_colsolve->m_gam2<<")"
            << " betas=("<<p_colsolve->m_beta1<<", "<<p_colsolve->m_beta2<<")"
            << " Eint2=("<<p_colsolve->m_eint1<<", "<<p_colsolve->m_eint2<<")"
            << "\n";
        /// solve equation for the shell collision
        /// set initial guess for the solver
        double i_gM = (p_colsolve->m_gam1*p_colsolve->m_mass1 + p_colsolve->m_gam2*p_colsolve->m_mass2)
                    / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        double i_mM = p_colsolve->m_mass1 + p_colsolve->m_mass2;
        double i_eM = (p_colsolve->m_eint1*p_colsolve->m_mass1+p_colsolve->m_eint2*p_colsolve->m_mass2)
                    / (p_colsolve->m_mass1+p_colsolve->m_mass2);
        int x_ = 1;
        // define a system of non-linear equations to solve
        auto func = []( int n, double x[], double fx[],
                                                               int &iflag, void * pars ){
            auto pp = (struct FsolvePars *) pars;
            double gM = x[0];
            double eM = x[1];
            double gAdiM = pp->p_eos->getGammaAdi(gM, EQS::Beta(gM));
            /// total mass conservation
            double mM = pp->m_mass1 + pp->m_mass2;
            /// total energy consercation (Ek + Eint)
            double fx1 = (pp->m_gam1 * pp->m_mass1 + EQS::get_GammaEff(pp->m_gam1, pp->m_adi1) * pp->m_eint1)
                       + (pp->m_gam2 * pp->m_mass2 + EQS::get_GammaEff(pp->m_gam2, pp->m_adi2) * pp->m_eint2)
                       - (gM * mM + EQS::get_GammaEff(gM,gAdiM) * eM );
            /// total momentum conservation
            double fx2 = sqrt(pp->m_gam1 * pp->m_gam1 - 1) * (pp->m_mass1 + pp->m_adi1 * pp->m_eint1 )
                       + sqrt(pp->m_gam2 * pp->m_gam2 - 1) * (pp->m_mass2 + pp->m_adi2 * pp->m_eint2 )
                       - sqrt(gM * gM - 1) * (mM + gAdiM * eM );
            if (!std::isfinite(fx1) || !std::isfinite(fx2))
                iflag = -1;

            if (!std::isfinite(sqrt(pp->m_gam1 * pp->m_gam1 - 1)))
                iflag = -1;

            fx[0] = fx1;//std::log10(fx1);// * fx1; TODO THis is wrong. It should not be sqred!
            fx[1] = fx2;//std::log10(fx2);// * fx2;
        };
        using func_ptr = void(*)(int, double *, double *, int&, void*);
        func_ptr p_func = func;

        /// solve the system using a 'fsolve'
        auto solve = [&]( ){
            double *fx;
            int iflag;
            int info;
            int lwa;
            int n = 2;
            double tol = 1e-3;
            double *wa;
            double *x;

            lwa = ( n * ( 3 * n + 13 ) ) / 2;
            fx = new double[n];
            wa = new double[lwa];
            x = new double[n];

//            (*p_log)(LOG_INFO,AT) << "\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve_test2():\n";
//            (*p_log)(LOG_INFO,AT) << "fsolve() solves a nonlinear system of 2 equations.\n";

            x[0] = i_gM;
            x[1] = i_eM;
            iflag = 0;
            func ( n, x, fx, iflag, p_colsolve );
//            r8vec2_print ( n, x, fx, "  Initial X, F(X)" );
            info = fsolve ( func, n, x, p_colsolve, fx, tol, wa, lwa );

            if (info<0){
                (*p_log)(LOG_ERR,AT)<< "Fsolve failed (try setting 'tol' lower). "
                                                   "Using initial guess values. New shell has "
                          <<"Gamma="<<i_gM<<" beta="<<EQS::Beta(i_gM)<<" Eint="<<i_eM<<"\n";
//                i_gM = x[0];
//                i_eM = x[1];
            }
            else{
                (*p_log)(LOG_INFO,AT)<< "Fsolve successful. New shell has "
                                             <<"Gamma="<<x[0]
                                             <<" beta="<<EQS::Beta(x[0])
                                             <<" Eint="<<x[1]<<"\n";
                i_gM = x[0];
                i_eM = x[1];
            }
//            std::cout << "\n";
//            std::cout << "  Returned value of INFO = " << info << "\n";
//            r8vec2_print ( n, x, fx, "  Final X, FX" );
            delete [] fx;
            delete [] wa;
            delete [] x;
        };
        solve();
        /// chose which shell to terminate (the slower one)
        int ish = p_colsolve->choseShell();
//        ish = 2;
        size_t iieq; size_t iieq_other;
        double eint_before, eint_after;
        double m2_before, m2_after;
        if (ish == 1){
            iieq = bw1->getPars()->ii_eq;
            iieq_other = bw2->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iEint2];
            eint_before = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            Y[iieq + DynRadBlastWave::Q_SOL::iEad2] += Y[iieq + DynRadBlastWave::Q_SOL::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] += Y[iieq + DynRadBlastWave::Q_SOL::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::Q_SOL::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::Q_SOL::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::Q_SOL::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            m2_after = Y[iieq + DynRadBlastWave::Q_SOL::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::Q_SOL::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::Q_SOL::itt] = 0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::Q_SOL::imom] = _mom;
        Y[iieq + DynRadBlastWave::Q_SOL::iEint2] = _eint2; // use ODE units

        (*p_log)(LOG_INFO,AT) << "\tLayer [ll="<<mylayer<<"] Outcome for"
                              << " idx1="<<idx1<<", idx2="<<idx2 << " collision:"
                              << " Eint2/M0c^2: "<<eint_before<<" -> "<<_eint2
                              << " M2/M0: "<<m2_before<<" -> "<<m2_after
                              <<"\n";
        if (!std::isfinite(i_gM)||(!std::isfinite(i_eM)||(!std::isfinite(i_mM)))){
            (*p_log)(LOG_ERR,AT)<<"Nan in collision result\n";
            exit(1);
        }
//


//
//        Y[iieq + DynRadBlastWave::Q_SOL::i] = i_gM * EQS::Beta(i_gM);


    }

    /// get the time at which shells collided
    double evalShellCollisionTime( double x, size_t ix, const double * Y ){
        /// if not sorted however, we need to check what blast wave overrun what blastwave(s)...
        size_t n_overruns = 0;
        Vector min_t_at_which_overrun;

        for (size_t i=0; i<m_nshells; i++) {
            double r_i = m_radii_init[i];
            std::vector<size_t> overrun_shells;
            Vector t_at_which_overrun;
            /// now check which blast-waves has this one overrun
            for (size_t j=i; j<m_nshells; j++){
                double r_j = m_radii_init[j];
                if (r_i > r_j){
                    overrun_shells.push_back(j);
                    /// interpolate the time at which shells collided
                    double tim1 = m_t_grid[ix-1];
                    double ti = m_t_grid[ix];
                    double r0im1 = p_bws_ej[i]->getVal(BlastWaveBase::iR,(int)(ix-1));
                    double r0i = Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ];
                    double r1im1 = p_bws_ej[j]->getVal(BlastWaveBase::iR,(int)(ix-1));
                    double r1i = Y[ p_bws_ej[j]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ];
                    double slope0 = (r0i - r0im1) / (ti - tim1);
                    double slope1 = (r1i - r1im1) / (ti - tim1);
                    double r0int = r0im1 - slope0 * tim1;
                    double r1int = r1im1 - slope1 * tim1;
                    double tcoll = (r1int - r0int) / (slope0 - slope1);
                    double rcoll = slope1 * tcoll + r0int;
                    if ((tcoll < tim1)or(tcoll>ti)){
                        (*p_log)(LOG_ERR,AT) << " inferred tcoll="<<tcoll
                        <<" is not in tim1-ti range ["<<tim1<<", "<<ti<<"]\n";
                        exit(1);
                    }
                    t_at_which_overrun.push_back( tcoll );
                }
//                if (overrun_shells.empty()){
//                    std::cerr << AT << " error.\n";
//                    exit(1);
//                }
            }
            if (overrun_shells.size() > 1){
                (*p_log)(LOG_ERR,AT)
                        <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                        <<" mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                        <<" (mom0="<<p_bws_ej[i]->getPars()->mom0<<") "
                        <<" has overrun n_shells="<<overrun_shells.size()
                        <<" with momenta ["<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                        <<" (mom0="<<p_bws_ej[overrun_shells[0]]->getPars()->mom0<<") "
                        <<" ... " << Y[ p_bws_ej[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                        <<"] \n";

//                (*p_log)(LOG_ERR,AT) << " shell i="<<i<<" has overrun n="<<overrun_shells.size()<<" other shells.\n";
                return -1;
            }
            if (overrun_shells.size() == 1){
                (*p_log)(LOG_WARN,AT)
                        <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                        <<" mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                        <<" (mom0="<<p_bws_ej[i]->getPars()->mom0<<") "
                        <<" has overrun shell j="<<overrun_shells[0]
                        <<" with momenta "<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                        <<" (mom0="<<p_bws_ej[overrun_shells[0]]->getPars()->mom0<<") "
                        <<"\n";
                n_overruns += 1;
            }
            /// get time till the nearest collision
            min_t_at_which_overrun.push_back()
        }
        if (n_overruns > 1){
            (*p_log)(LOG_ERR,AT)

                    << " at t="<<x << " multimple overruns have occured n_overruns="<<n_overruns<<"\n";
//                    <<" with mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
//                    <<" has overrun n_shells="<<overrun_shells.size()
//                    <<" with momenta ["<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ]
//                    <<" ... " << Y[ p_bws_ej[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ]
//                    <<" \n";

//            (*p_log)(LOG_ERR,AT) << "n_overruns="<<n_overruns<<"\n";
            return -1;
        }
        return 1;
    }
    int areShellsOrderedRadially( double x, size_t ix, const double * Y ){
        /// collect radii of the blast wave
        double r = -1.;
        std::vector<size_t> idxs_to_collide;
        for (size_t i=0; i<m_nshells; i++) {
            m_idxs[i] = i;
            if (p_bws_ej[i]->getPars()->end_evolution)
                continue;
            m_radii_init[i] = Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ];
        }

        /// check if radii are sorted

        bool is_sorted = std::is_sorted(m_radii_init.begin(),m_radii_init.end());
        /// if sorted returm 0 as there is no collision and no problem)
        if (is_sorted)
            return 0;

        /// if not sorted however, we need to check what blast wave overrun what blastwave(s)...
        size_t n_overruns = 0;
        for (size_t i=0; i<m_nshells; i++) {
            double r_i = m_radii_init[i];
            std::vector<size_t> overrun_shells;
            /// now check which blast-waves has this one overrun
            for (size_t j=i; j<m_nshells; j++){
                double r_j = m_radii_init[j];
                if (r_i > r_j){
                    overrun_shells.push_back(j);
                }
//                if (overrun_shells.empty()){
//                    std::cerr << AT << " error.\n";
//                    exit(1);
//                }
            }
            if (overrun_shells.size() > 1){
                (*p_log)(LOG_ERR,AT)
                    <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                    <<" mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                    <<" (mom0="<<p_bws_ej[i]->getPars()->mom0<<") "
                    <<" has overrun n_shells="<<overrun_shells.size()
                    <<" with momenta ["<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                    <<" (mom0="<<p_bws_ej[overrun_shells[0]]->getPars()->mom0<<") "
                    <<" ... " << Y[ p_bws_ej[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                    <<"] \n";

//                (*p_log)(LOG_ERR,AT) << " shell i="<<i<<" has overrun n="<<overrun_shells.size()<<" other shells.\n";
                return -1;
            }
            if (overrun_shells.size() == 1){
                (*p_log)(LOG_WARN,AT)
                <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                <<" mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                <<" (mom0="<<p_bws_ej[i]->getPars()->mom0<<") "
                <<" has overrun shell j="<<overrun_shells[0]
                <<" with momenta "<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
                <<" (mom0="<<p_bws_ej[overrun_shells[0]]->getPars()->mom0<<") "
                <<"\n";
                n_overruns += 1;
            }
        }
        if (n_overruns > 1){
            (*p_log)(LOG_ERR,AT)

                    << " at t="<<x << " multimple overruns have occured n_overruns="<<n_overruns<<"\n";
//                    <<" with mom="<<Y[ p_bws_ej[i]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom ]
//                    <<" has overrun n_shells="<<overrun_shells.size()
//                    <<" with momenta ["<<Y[ p_bws_ej[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ]
//                    <<" ... " << Y[ p_bws_ej[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR ]
//                    <<" \n";

//            (*p_log)(LOG_ERR,AT) << "n_overruns="<<n_overruns<<"\n";
            return -1;
        }
        return 1;
    }

#endif

    void evalShellThicknessIsolated(size_t idx, const double * Y){

        double r_i=0., r_ip1=0., dr_i = 0., vol_i=0.;
        double frac = 0.1; // fraction of the shell volume to be used as its width. Inaccurate as we need adjacent shells to get the volume...
        auto & bw = p_bws_ej[idx];
        r_i =  Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
        dr_i = frac * r_i;
        vol_i = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
        bw->getPars()->thickness = dr_i;
        bw->getPars()->volume = vol_i;
    }
    /// Evaluate shell thickness using various methods. Assumes sorted shells (no overruns)
    void evalShellThickness( const double * Y ){
        double r_i=0., r_ip1=0., dr_i = 0., vol_i=0.;
        /// if there is one shell we cannot have the shell width that comes from shell separation.
        if (n_active_shells == 1){ evalShellThicknessIsolated(0, Y); }
        else{
            for (size_t ii=0; ii<n_active_shells-1; ii++) {
                r_i = 0., r_ip1 = 0., dr_i = 0., vol_i = 0.;
                ///
                size_t idx = m_idxs[ii];
                size_t nextidx = m_idxs[ii + 1];
                auto &bw = p_bws_ej[idx];
                auto &nextbw = p_bws_ej[nextidx];
                if ((bw->getPars()->end_evolution) || (nextbw->getPars()->end_evolution)){
                    evalShellThicknessIsolated(idx, Y);
                }
//            double r_cur = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
//            double r_cur = bw->getVal(DynRadBlastWave::Q::iR, 0);
                r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
                r_ip1 = Y[nextbw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
                if ((r_i == 0.) || (r_ip1 == 0.)) {
//                (*p_log)(LOG_WARN,AT)<<" shell="<<idx<<" has r_i="<<r_i<<" and r_ip1="<<r_ip1<<"\n";
                    continue;
                }
                dr_i = r_ip1 - r_i;
                /// evaluate the volume of the shell (fraction of the 4pi)
                vol_i = (4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i) / bw->getPars()->ncells;
                /// --------------------------- |
                bw->getPars()->thickness = dr_i;
                bw->getPars()->volume = vol_i;

            }
        }
    }

    /// Evaluate the radial extend of a velocity shell. Assume ordered shells. Assumes sorted shells
    void evalShellOptDepth(  const double * Y ){

        double kappa_i = 5.e50; // bw->getPars()->kappa0; // TODO!

        double r_i=0., r_ip1=0., dr_i = 0., m_i=0., m2_i=0., m_ip1=0.,
                m2_ip1=0., vol_i=0., rho_i=0., dtau_i=0., tau_tot=0.;

        /// if there is one shell we cannot have the shell width that comes from shell separation.
        if (n_active_shells == 1){
            double frac = 0.1; // fraction of the shell volume to be used as its width. Inaccurate as we need adjacent shells to get the volume...
            auto & bw = p_bws_ej[0];
//            r_i =  Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
            dr_i = bw->getPars()->thickness;//frac * r_i; // TODO add other methods 1/Gamma...
            vol_i = bw->getPars()->volume;// = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
            m_i = bw->getPars()->M0;
            m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] * m_i;
            m_i += m2_i;
            rho_i = m_i / vol_i;
            m_dtau[0] = kappa_i * rho_i * dr_i;
            m_dtau_cum[0] = 0.; // opt. depth to this shell
        }

        /// compute shell width from the shell separation
        for (size_t ii=0; ii<n_active_shells-1; ii++){
            r_i=0., r_ip1=0., dr_i = 0., m_i=0., m2_ip1=0., m_ip1=0., vol_i=0., rho_i=0., dtau_i=0.;
            ///
            size_t idx = m_idxs[ii];
            size_t nextidx = m_idxs[ii+1];
            auto & bw = p_bws_ej[idx];
            auto & nextbw = p_bws_ej[nextidx];

            dr_i = bw->getPars()->thickness;//r_ip1 - r_i;

            /// evaluate mass within the shell
            m_i = bw->getPars()->M0;
            m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] * m_i;
            m_i += m2_i;
            m_ip1 = nextbw->getPars()->M0;
            m2_ip1 = Y[nextbw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iM2] * m_ip1;
            m_ip1 += m2_ip1;
            /// evaluate the volume of the shell (fraction of the 4pi)
            vol_i = bw->getPars()->volume;//(4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i) / bw->getPars()->ncells;
            /// evaluate density within a shell (neglecting the accreted by the shock!)
            rho_i = m_i / vol_i;
            m_rho[ii] = rho_i;

            /// evaluate optical depth
            dtau_i = kappa_i * rho_i * dr_i;
            tau_tot += dtau_i;
            if ((!std::isfinite(dtau_i)) || (dtau_i < 0)){
                (*p_log)(LOG_ERR,AT) << "dtau is nan or < 0 and dtau="<<dtau_i<<"\n";
                exit(1);
            }
            m_dtau[ii] = dtau_i;
        }
//        (*p_log)(LOG_INFO,AT)<<" optical depth is evaluated. Total tau="<<tau_tot<<"\n";


        /// based on where tau = 1 get the photosphere radius, i.e., in which shell is the photosphere
        double tmp_tau = 0;
        size_t idx_photo = n_active_shells;
        for (size_t i = n_active_shells-1; i > 0; i--){
            size_t idx = m_idxs[i];
            auto & bw = p_bws_ej[idx];
//            if (bw->getPars()->M0 == 0)
//                continue;
            tmp_tau+=m_dtau[idx];
            if (tmp_tau > 1.){
                idx_photo = (int)idx+1;
                break;
            }
        }
//        (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<mylayer
//                                        <<"] Photosphere located at idx="<<idx_photo<<"\n";


        /// Compute the optical depth from 0 to a given shell
//        for (auto & cur_idx : m_idxs){
        for (size_t ii = 0 ; ii < n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            double tau_cum = 0.;
//            for (auto & other_idx : m_idxs){
            for (size_t jj = 0 ; jj < n_active_shells; jj++){
                size_t other_idx = m_idxs[jj];
                if (other_idx < cur_idx)
                    tau_cum += m_dtau[other_idx];
            }
            m_dtau_cum[cur_idx] = tau_cum;
        }
//        (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<mylayer
//                                        <<"] Cumulative optical depth to idx0"
//                                        <<" tau0="<<m_dtau_cum[0]
//                                        <<" tau[-1]="<<m_dtau_cum[n_active_shells-1]<<"\n";

        /// save results for the use in ODE when solving Energy equation
//        for (auto & cur_idx : m_idxs){
        for (size_t ii = 0 ; ii < n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            auto & bw = p_bws_ej[cur_idx];
            bw->getPars()->tau_to0 = m_dtau_cum[cur_idx];
            bw->getPars()->dtau = m_dtau[cur_idx];
            bw->getPars()->is_above_tau1 = idx_photo > cur_idx ? false : true;
        }
//        std::cout << m_dtau_cum << "\n";
//        std::cout << m_dtau << "\n";
//        std::cout << m_rho << "\n";
        int x = 1;

//        if ((idx_photo < 10) && (idx_photo > 0)){
//            int a = 1;
//            std::cout << " i_photo="<<idx_photo << " TauToIt="<<m_dtau_cum[idx_photo]<< " remaining tau="<<tau_tot-m_dtau_cum[idx_photo]<<'\n';
//        }


//        for (size_t i = 0; i<m_nshells; i++){
//            double tau_cum = 0.;
//            for (size_t j = 0; j<m_nshells; j++) {
//                if (i > j)
//                    continue;
//                else
//                    tau_cum += m_dtau[i];
//            }
//            m_dtau_cum[i] = tau_cum;
//        }
//        if (idx_photo>0)
//            std::cout << " i_photo="<<idx_photo << " TauToIt="<<m_dtau_cum[idx_photo]<< " remaining tau="<<tau_tot-m_dtau_cum[idx_photo]<<'\n';
        /// for the last shell i don't have difference, so... just the same as previous one.
//        size_t ish = m_idxs.back();
//        p_bws_ej[ish]->getLastVal(DynRadBlastWave::Q::idR) = dr_i;
//        p_bws_ej[ish]->getLastVal(DynRadBlastWave::Q::iRho) = rho_i;
//        p_bws_ej[ish]->getLastVal(DynRadBlastWave::Q::iTau) = tau_i; // update blastwave thickness

//        if (idx_photo == m_idxs.back()-1)
//            idx_photo = (int)m_idxs.back();
//        for (auto idx : m_idxs)
//            std::cout << "i_photo=" <<idx_photo<< "; "<< p_bws_ej[idx]->getLastVal(DynRadBlastWave::Q::iTau) << ", ";
//        std::cout <<"\n";

//        double tau_out = 0, tau_out_i = 0;
//        for (auto idx_curr : m_idxs){
//            tau_out = 0.;
//            for (auto idx = m_idxs.back(); idx > idx_curr; idx--){
//                tau_out_i = p_bws_ej[idx]->getLastVal(DynRadBlastWave::Q::iTau);
//                tau_out+=tau_out_i ;


    }

#if 0
    inline std::vector<size_t> & getCurrentIndexes( ) { return m_idxs; }
    void evaluateInitialShellLength(){

    }

    void evaluateOptDepthAtCurrR(size_t i_cur, double x, const double * Y){

        double kappa = 10.;

        /// evaluate relative positions of shells at this time
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws_ej[i]->getLastVal(std::static_pointer_cast<BlastWaveBase::Q>())
        }
        /// get indexes of sorted shells by radius
        sort_indexes(m_radii, m_idxs);
        /// evaluate opacity of each shell up to current
        double tau = 0.;
        for (auto & idx : m_idxs){
            auto & bw = p_bws_ej[idx];
            double dr = 0.;
            if (idx == 0)
                dr =0.5 * (m_radii[2]-m_radii[1]);
            else if (idx == m_idxs.back())
                dr =(m_radii[idx]-m_radii[idx-1]);
            else
                dr = 0.5 * (m_radii[idx+1]-m_radii[idx-1]);
//            if (idx == 0) {
//                double dr1 = 0.5 * (m_radii[0] - m_radii[1]);
//                double dr2 = 0.5 * (m_radii[1] - m_radii[2]);
//                dr = linearExtrapolate(m_radii[1], m_radii[2], dr1, dr2, m_radii[0]);
//            }
//            else if (idx==m_idxs.back()){
//                double drnm2 = 0.5 * (m_radii[m_radii.size()-2] - m_radii[m_radii.size()-3]);
//                double drnm1 = 0.5 * (m_radii[m_radii.size()-1] - m_radii[m_radii.size()-2]);
//                dr = linearExtrapolate(m_radii[1], m_radii[2], drnm1, drnm2, m_radii[0]);
            double dvolume = (4./3.) * CGS::pi * dr * dr * dr;
            double drho = bw->getPars()->M0 / dvolume;
            double dtau = kappa * drho * dvolume;
        }


        /// optical depth at a given radius
        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws_ej[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

    }





    inline void evalRelativePosition(double x, const double * Y){
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws_ej[i]->getLastVal(std::static_pointer_cast<BlastWaveBase::Q>())
        }
        sort_indexes(m_radii, m_idxs);

        iota(current_indexes.begin(), current_indexes.end(), 0);
        stable_sort(current_indexes.begin(), current_indexes.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
//        std::cout << radii << "\n";
        current_indexes = sort_indexes(m_radii);
//        std::cout << idxes << "\n";
//        if (idxes[0]!=0){
//            std::cout << x << " " << x / CGS::day << " " << x / CGS::year;
//            std::cout << radii << "\n" ;
//            std::cout << idxes << "\n" ;
//            exit(1);
//        }
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            bw->getPars()->i_position_in_layer = current_indexes[i];
        }

        std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws;
        for (size_t & idx : current_indexes){
            p_sorted_bws.emplace_back(p_bws_ej[idx]);
        }

    }
    ///
    void evaluateDensityProfileFromBWpositions(){
        for (size_t & idx : current_indexes) {
            auto & bw = p_bws_ej[idx];

            double rho_i = bw->getPars()->M0 / (4./3.*CGS::pi*())
        }

        Vector rho ( m_nshells, 0.0 );
        for (size_t i = 0; i < m_nshells; i++){
            size_t idx = current_indexes[i];

        }
    }

    void evaluateOpticalDepthToMagnetarEmission(size_t i_cur, double x, const double * Y){
        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws_ej[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            if (r_cur > bw->)
        }

        for (size_t i = 0; i < i_position; i++){
            tau +=
        }
    }

#endif
};

class GRB{
    std::unique_ptr<logger> p_log;
//    bool run_jet_bws = false;
    bool is_jBW_init = false;
    bool is_jet_obsrad_pars_set = false;
    bool is_jet_struct_set = false;
    int n_ode_eq;
    LatStruct::METHOD_eats jet_eats_method{};
    int m_loglevel{};
    Array & t_arr;
    double jet_layer_fnu_stop_frac;
protected:
    LatStruct jetStruct{};
    std::vector<std::unique_ptr<RadBlastWave>> p_bws_jet;
public:
    double jet_rtol = 1e-5;
    bool run_jet_bws = false;
    bool is_jet_obs_pars_set = false;
    bool is_jet_anal_synch_computed = false;
    GRB(Array & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "GRB");
    }
    size_t getNeq() const {
        if (!run_jet_bws)
            return 0;
        else
            return (p_bws_jet.size() * p_bws_jet[0]->getNeq());
    }
    void setJetStructAnalytic(StrDbMap pars, StrStrMap opts){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true, true);
        if (!run_jet_bws)
            return;
        std::string _jet_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        jet_eats_method = LatStruct::setEatsMethod(_jet_eats_method);
        jetStruct.initAnalytic( pars, opts, _jet_eats_method, m_loglevel );
        is_jet_struct_set = true;
    }
    std::vector<std::unique_ptr<RadBlastWave>> & getBWs(){return p_bws_jet;}
//    static auto listParsNumericJetStruct(){ return LatStruct::listParsCustomStruct(); }
    void setJetStructNumeric( Vector & dist_thetas, Vector & dist_EEs, Vector & dist_Gam0s, Vector & dist_MM0s,
                              bool force_grid, std::string eats_method ){
        (*p_log)(LOG_ERR,AT) << " not finished...\n"; exit(1);
        jetStruct.initCustom( dist_thetas, dist_EEs, dist_Gam0s, dist_MM0s, force_grid, eats_method,
                              m_loglevel);
        is_jet_struct_set = true;
    }
    void setJetBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true, true);
        if (!run_jet_bws)
            return;

        size_t n_layers = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
        for(size_t i = 0; i < n_layers; i++){
//            DynRadBlastWave x(t_grid, 0, i, p_pars->loglevel);
            p_bws_jet.emplace_back( std::make_unique<DynRadBlastWave>(t_arr, 0, i, m_loglevel) );
            p_bws_jet[i]->setAllParametersForOneLayer(jetStruct, pars, opts, i, ii_eq);
            p_bws_jet[i]->setEatsPars(pars, opts);
            p_bws_jet[i]->getSynchAnPtr()->setPars( pars, opts );
            ii_eq += p_bws_jet[i]->getNeq();
        }

        /// parameters for EATS methods
//        p_pars->jet_eats_method = LatStruct::setEatsMethod(
//                getStrOpt("method_eats",opts,AT,p_log,"",true));
        jet_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);
        jet_layer_fnu_stop_frac= getDoublePar("fnu_min_frac",pars,AT,p_log,-1, true);

        (*p_log)(LOG_INFO,AT) << "finished initializing jet...\n";
        is_jet_obs_pars_set = true;
        is_jBW_init = true;
    }
    void infoFastestLayer(size_t it, logger &sstream){
        size_t fastest = 0;
        double mom = 0;
        for (size_t il = 0; il < jetStruct.nlayers; il++){
            if (p_bws_jet[il]->getVal(RadBlastWave::Q::imom,(int)it) > mom) {
                mom = p_bws_jet[il]->getVal(RadBlastWave::Q::imom, (int) it);
                fastest = il;
            }
        }
        sstream << "[Jet: "<<"[l="<<fastest<<"]"
                << " Mom=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::imom, (int) it))
                << " R=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::iR, (int) it))
                << " Eint=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::iEint2, (int) it))
                << "] ";
    }
};

class RadGRB : protected GRB{
    std::unique_ptr<logger> p_log;
    Image tmp_im_pj;
    Image tmp_im_cj;
    Image tmp_im;
    RadGRB(Array & t_arr, int loglevel) : GRB(t_arr, loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "RadGRB");
        auto ncells = p_bws_jet[0]->getPars()->ncells;
        (*p_log)(LOG_INFO, AT) << " allocating memory for images with ncells="<<ncells<<"\n";
        tmp_im_pj.resize((int)ncells); tmp_im_cj.resize((int)ncells); tmp_im.resize((int)ncells);
    }


private:

public:
    void getJetLightCurve( VecVector & lc, Vector & obs_times, Vector & obs_freqs ){
        (*p_log)(LOG_INFO,AT)<<" starting grb light curve calculation\n";
        size_t nlayers = p_bws_jet.size();
        auto ncells = p_bws_jet[0]->getPars()->ncells;
        VecVector light_curves(nlayers); // [i_layer][i_time]
        light_curves.resize(nlayers);
        for (auto & arr : light_curves)
            arr.resize(obs_times.size(), 0.);
        double flux;
        double rtol = jet_rtol;
        for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
            auto & model = p_bws_jet[ilayer];
            (*p_log)(LOG_INFO,AT)<<" GRB LC ntimes="<<obs_times.size()<<" theta_layer="<<ilayer<<"/"<<nlayers<<
                                 " phi_cells="<<LatStruct::CellsInLayer(model->getPars()->ilayer)<<"\n";
            model->evalForwardShockLightCurve(tmp_im, tmp_im_pj, tmp_im_cj, light_curves[ilayer], obs_times, obs_freqs);
        }
        (*p_log)(LOG_INFO,AT)<<" grb light curve is computed\n";
    }

};

class Ejecta{
    VelocityAngularStruct ejectaStructs{};
    std::vector<std::unique_ptr<CumulativeShell>> p_ej;
    std::unique_ptr<logger> p_log;
    bool is_ejBW_init = false;
    bool is_ejecta_obsrad_pars_set = false;
    bool is_ejecta_struct_set = false;
    double jet_layer_fnu_stop_frac=1e-5;
    int n_ode_eq;
    int m_loglevel;
    LatStruct::METHOD_eats ejecta_eats_method{};
    Array & t_arr;
public:
    bool run_ej_bws = false;
    bool do_collision = false;
    bool is_ejecta_obs_pars_set = false;
    bool is_ejecta_anal_synch_computed = false;
    Ejecta(Array & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Ejecta");
    }
    size_t getNeq() const {
        if (!run_ej_bws)
            return 0;
        else
            return (p_ej.size() * p_ej[0]->nBWs() * p_ej[0]->getBW(0)->getNeq());
    }
    size_t nlayers() const { return p_ej.size(); }
    size_t nshells() const { return p_ej[0]->getBWs().size(); }
    int ncells() const { return (int)p_ej[0]->getBW(0)->getPars()->ncells; }
    std::vector<std::unique_ptr<CumulativeShell>> & getShells(){ return p_ej;}
    /// set ejecta lateral & velocity structure
//    static auto listParsAnalyticEjectaStruct(){ std::cerr << AT << " not implemented\n"; exit(1); }
    void setEjectaStructAnalytic(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_ERR,AT) << " not implimeneted\n Exiting...";
        exit(1);
    }
//    static std::vector<std::string> listParsNumericEjectaStruct(){ return VelocityAngularStruct::list_pars_v_ns(); }
    void setEjectaStructNumericUniformInTheta(Vector & dist_thetas0, Vector & dist_betas, Vector & dist_ek, size_t nlayers, double mfac, StrStrMap & opts){
        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log, false, true);
        if (!run_ej_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initUniform(dist_thetas0,dist_betas,dist_ek, nlayers,mfac,
                                  ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaStructNumeric(Vector dist_thetas, Vector dist_betas, VecVector dist_ek, bool force_grid, StrStrMap & opts){
        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log,false, true);
        if (!run_ej_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initCustom(dist_thetas, dist_betas, dist_ek, force_grid,
                                 ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq, size_t n_layers_jet){

        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log, false, true);
        if (!run_ej_bws)
            return;
        bool is_within = false;
        std::vector<size_t> which_within{};
        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t nshells = ejectaStructs.nshells;
        size_t n_layers_ej = ejectaStructs.structs[0].nlayers;
        std::vector<std::vector<size_t>> n_empty_images_layer_shell;
        for (auto & n_empty_images_layer : n_empty_images_layer_shell)
            n_empty_images_layer.resize(nshells);
        /// include blastwave collision between velocioty shells into the run
        for(size_t il = 0; il < n_layers_ej; il++){
            p_ej.push_back( std::make_unique<CumulativeShell>(t_arr, nshells, il,
                                                              p_log->getLogLevel()) );
            for (size_t ish = 0; ish < nshells; ish++){
                auto & bw = p_ej[il]->getBW(ish);
                auto & struc = ejectaStructs.structs[ish];
                bw->setAllParametersForOneLayer(struc, pars, opts, il, ii_eq);
/// Override the layer-to-use
                if (bw->getPars()->which_jet_layer_to_use == 0){
                    bw->getPars()->which_jet_layer_to_use = 0; // the fastest
                }
                else if(n_layers_ej==0){
                    n_ejecta_empty_images += 1;
                    n_empty_images_layer_shell[ish].emplace_back(il);
//                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
                }
                else if(n_layers_ej == 0){
                    // NO jet structure was set, so exiting I guess... :)
                    // TODO THIS MIGHT BE WRONG -- why 'n_layers_i'
                }
                else if(n_layers_jet == 0){
                    // NO jet structure was set, so exiting I guess... :)
                }
                else if ((bw->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
                    bw->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
                }
                else if ((bw->getPars()->which_jet_layer_to_use < n_layers_jet) &&
                         (bw->getPars()->which_jet_layer_to_use > -1)){
                    //
                }
                else{
                    (*p_log)(LOG_ERR,AT) << " which_jet_layer_to_use="<<bw->getPars()->which_jet_layer_to_use
                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
                              <<" for the slowest"
                              <<" or any N in between the two for a specific jet layer \n"
                              << "Exiting..."
                              << "\n";
                    exit(1);
                }
                ii_eq += bw->getNeq();

                bw->setEatsPars(pars, opts);
                bw->getSynchAnPtr()->setPars( pars, opts );

//                ii++;
            }
        }
# if 0
        size_t ii = 0;
        bool is_within = false;
        std::vector<size_t> which_within{};

        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t n_layers_i = ejectaStructs.structs[0].nlayers;

        for(size_t i = 0; i < ejectaStructs.nshells; i++){
            std::vector<size_t> n_empty_images_layer;
            auto & struc = ejectaStructs.structs[i];
            if (n_layers_i != struc.nlayers){
                (*p_log)(LOG_ERR,AT)<<" expected nlayers="<<n_layers_i<<" (from 0th shell) got="<<struc.nlayers<<"\n";
                exit(1);
            }
//            size_t n_layers_i = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for(size_t j = 0; j < n_layers_i; j++) { // ii = il + nlayers * ish
                p_bws_ej.emplace_back(std::make_unique<DynRadBlastWave>(t_grid, i, j, p_pars->loglevel));
                setAllParametersForOneLayer(struc, *(p_bws_ej[ii]), pars, opts, j, ii_eq);
                /// Override the layer-to-use
                if (p_bws_ej[ii]->getPars()->which_jet_layer_to_use == 0){
                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = 0; // the fastest
                }
                else if(n_layers_i==0){
                    n_ejecta_empty_images += 1;
                    n_empty_images_layer.emplace_back(j);
//                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
                }
                else if(n_layers_i == 0){
                    // NO jet structure was set, so exiting I guess... :)
                    // TODO THIS MIGHT BE WRONG -- why 'n_layers_i'
                }
                else if(n_layers_jet == 0){
                    // NO jet structure was set, so exiting I guess... :)
                }
                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
                }
                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use < n_layers_jet) &&
                         (p_bws_ej[ii]->getPars()->which_jet_layer_to_use > -1)){
                    //
                }
                else{
                    std::cerr << " which_jet_layer_to_use="<<p_bws_ej[ii]->getPars()->which_jet_layer_to_use
                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
                              <<" for the slowest"
                              <<" or any N in between the two for a specific jet layer \n"
                              << "Exiting..."
                              << "\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }
                ii_eq += p_bws_ej[ii]->getNeq();

                p_bws_ej[ii]->setEatsPars(pars, opts);
                p_bws_ej[ii]->getSynchAnPtr()->setPars( pars, opts );

                ii++;
            }

            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(i);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
        }
#endif
        is_ejBW_init = true;
        is_ejecta_obs_pars_set = true;

        if ((p_log->getLogLevel() > LOG_WARN)) {
            if (n_ejecta_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta blastwave is NOT initialized for total n="
                      << n_ejecta_empty_images << " layers. Specifically:\n";
                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
                    auto &ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
                    size_t n_layers_i = ejectaStruct.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [";
                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
                        ccerr << n_empty_images[ish][il] << " ";
                    }
                    ccerr << "] / (" << n_layers_i << " total layers) \n";
                }
            }
        }

        ej_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);

        std::string method_collision = getStrOpt("method_collision", opts, AT,p_log, "none", true);
        if (method_collision != "none") {
            do_collision = true;
            (*p_log)(LOG_INFO,AT)<<"Including shell collision into kN ejecta dynamics\n";
        }

        (*p_log)(LOG_INFO,AT) << "finished initializing ejecta...\n";
    }
    double ej_rtol = 1e-5;
    void infoFastestShell(size_t it, const double * Ym1, const double * Y, logger sstream){
        size_t n_active_min = std::numeric_limits<size_t>::max(); int il_with_min_nact = -1;
        size_t n_active_max = 0; int il_with_max_nact = -1;
        size_t n_accel_max = 0; int il_with_n_accel_max = -1;
        size_t n_decel_max = 0; int il_with_n_decel_max = -1;

        double Mom_max_over_Gamma0 = 0;
        int il_wich_fastest = -1; int ish_with_fastest = -1;
        double Eint2max = 0.;
        int il_wich_energetic = -1; int ish_with_energetic = -1;
        double Mom_min_over_Gamma0 = std::numeric_limits<double>::max();
        int il_with_slowest = -1; int ish_with_slowest = -1;
        /// collect info in active shells
        for (size_t il = 0; il < nlayers(); il++ ){
            if (p_ej[il]->n_active_shells > n_active_max){
                n_active_max = p_ej[il]->n_active_shells;
                il_with_min_nact = (int)il;
            }
            if (p_ej[il]->n_active_shells < n_active_min){
                n_active_min = p_ej[il]->n_active_shells;
                il_with_max_nact = (int)il;
            }
            /// find number of shells in each layer that (i) accelerating (ii) decelerating
            size_t n_accel = 0;
            size_t n_decel = 0;
            for (size_t ish = 0; ish < p_ej[il]->n_active_shells; ish++){
                auto & bws = p_ej[il]->getBW(ish);
                double MomIm1 = Ym1[bws->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom];
                double MomI = Y[bws->getPars()->ii_eq + DynRadBlastWave::Q_SOL::imom];
                double Eint2I = Y[bws->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2];
                double Mom0 = bws->getPars()->mom0;
                if (MomI > MomIm1){
                    /// acceleration
                    n_accel += 1;
                }
                else if (MomI < MomIm1){
                    /// deceleration
                    n_decel += 1;
                }
                /// find fastest
                if (MomI/Mom0 > Mom_max_over_Gamma0){
                    Mom_max_over_Gamma0 = MomI/Mom0;
                    il_wich_fastest = (int)il;
                    ish_with_fastest = (int)ish;
                }
                /// find slowest
                if (MomI/Mom0 < Mom_min_over_Gamma0){
                    Mom_min_over_Gamma0 = MomI/Mom0;
                    il_with_slowest = (int)il;
                    ish_with_slowest = (int)ish;
                }
                /// most energetic
                if (Eint2max < Eint2I){
                    Eint2max = Eint2I;
                    il_wich_energetic= (int)il;
                    ish_with_energetic = (int)ish;
                }
            }
            /// record layer with the maximum number of accelerating and decelerating shells
            if (n_accel > n_accel_max){
                n_accel_max = n_accel;
                il_with_n_accel_max = (int)il;
                int x = 1;
            }
            if (n_decel > n_decel_max){
                n_decel_max = n_decel;
                il_with_n_decel_max = (int)il;
            }
        }

        sstream << "Ej:"
            <<" [NActive max/min="<<n_active_max<<"/"<<n_active_min<<" (il="<<il_with_max_nact<<"/"<<il_with_min_nact<<")"
            <<" MAX_Nacc/dec="<<n_accel_max<<"/"<<n_decel_max<<" (il="<<il_with_n_accel_max<<"/"<<il_with_n_decel_max<<")"
            <<" Mmax/M0="<<Mom_max_over_Gamma0<<" (il="<<il_wich_fastest<<", ish="<<ish_with_fastest<<")"
            <<" Mmin/M0="<<Mom_min_over_Gamma0<<" (il="<<il_with_slowest<<", ish="<<ish_with_slowest<<")"
            <<" Eint2max="<<Eint2max<<" (il="<<il_wich_energetic<<", ish="<<ish_with_energetic<<")"
            <<"]";

//
//
//        size_t fastest_sh = 0;
//        size_t fastest_l = 0;
//        double mom = 0;
//        double n_active_min = 0;
//        double layer
//        size_t n_acc = 0;
//        for (size_t ish = 0; ish < ejectaStructs.nshells; ish++){
//            for (size_t il = 0; il < ejectaStructs.structs[0].nlayers; il++){
//                if (p_ej[il]->getBW(ish)->getVal(RadBlastWave::Q::imom,(int)it) > mom) {
//                    mom = p_ej[il]->getBW(ish)->getVal(RadBlastWave::Q::imom,(int)it) > mom;
//                    fastest_l = il;
//                    fastest_sh = ish;
//                }
//            }
//        }
//        sstream << "[Ej: "<<"[l="<<fastest_l<<", sh="<<fastest_sh<<"]"
//                << " Mom=" << string_format("%.2e",p_ej[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::imom,(int)it))
//                << " R=" << string_format("%.2e",p_ej[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::iR,(int)it))
//                << " Eint=" << string_format("%.2e",p_ej[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::iEint2,(int)it))
//                << "] ";
    }
};

class EvolveODEsystem{
    struct Pars{
        Pars(
                std::unique_ptr<Magnetar> & p_magnetar,
                std::unique_ptr<GRB> & p_grb,
                std::unique_ptr<Ejecta> & p_ej,
//                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
//                std::vector<std::unique_ptr<CumulativeShell>> & p_ej,
//                bool run_magnetar, bool run_jet_bws, bool run_ej_bws,
                Array & t_grid
        ) :
        p_magnetar(p_magnetar), p_grb(p_grb), p_ej(p_ej),
//            run_magnetar(run_magnetar), run_jet_bws(run_jet_bws), run_ej_bws(run_ej_bws),
            t_grid(t_grid){
//            if (!p_bws_jet.empty()) {
//                n_j_bws = p_bws_jet.size();
//                n_eq_j_bws = p_bws_jet[0]->getNeq();
//            }
//            else {
//                n_j_bws = 0;
//                n_eq_j_bws = 0;
//            }
//            if(!p_ej.empty()) {
//                n_ej_bws = p_ej.size()*p_ej[0]->nBWs();//p_bws_ej.size();
//                n_eq_ej_bws= p_ej[0]->getBW(0)->getNeq();//p_bws_ej[0]->getNeq();
//            }
//            else {
//                n_ej_bws = 0;
//                n_eq_ej_bws = 0;
//            }
//            if (run_magnetar){
//                n_magnetar_eqs = p_magnetar->getNeq();
//            }
//            else{
//                n_magnetar_eqs = 0;
//            }
////            n_eq_j_bws = p_bws_jet[0]->getNeq();
////            n_eq_ej_bws= p_bws_ej[0]->getNeq();
//            n_tot_eqs  = n_magnetar_eqs
//                       + n_j_bws // number of 'jet' blast waves
//                       * n_eq_j_bws // number of equations in each
//                       + n_ej_bws // number of 'ejecta' blast waves
//                       * n_eq_ej_bws; // number of equations in each;
//            if (n_tot_eqs < 1){
//                std::cerr << AT << "\n No equations to evolve; Both, jet and ejecta seems to be not given. Exiting...\n";
//                exit(1);
//            }
        }
        Array & t_grid;
        std::unique_ptr<Magnetar> & p_magnetar;
//        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet;
//        std::vector<std::unique_ptr<CumulativeShell>> & p_ej;
        std::unique_ptr<GRB> & p_grb;
        std::unique_ptr<Ejecta> & p_ej;
//        bool run_magnetar{};
//        bool run_jet_bws{};
//        bool run_ej_bws{};
//        size_t n_magnetar_eqs = 0;
//        size_t n_j_bws = 0;
//        size_t n_ej_bws = 0;
//        size_t n_eq_j_bws = 0;
//        size_t n_eq_ej_bws = 0;
//        size_t n_tot_eqs = 0;
//        // ---
//        size_t n_layers_j = 0;
//        size_t n_layers_ej = 0;
//        size_t n_shells_j = 0;
//        size_t n_shells_ej = 0;
        // ---
//        bool use_cthetajet_lim_for_rho = false;
//        bool use_rjet_lim_for_rho = false;
        // ---
//        bool is_entered = false;
//        double r_ej_ent = -1;
        // ----

//        inline size_t idx_j(const size_t i){
////            auto * p_pars = (struct Pars *) pars;
//            auto & bwPtr = p_bws_jet[i];
//            return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_j * bwPtr->getPars()->ishell);
//        }
//        inline size_t idx_ej(const size_t i){
////            auto * p_pars = (struct Pars *) pars;
//            auto & bwPtr = p_bws_ej[i];
//            size_t _idx_j = n_eq_j_bws * n_j_bws;
//            return _idx_j + bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_ej * bwPtr->getPars()->ishell);
//        }
        // ---
        size_t ix = 0;  // index latest solution
        double dx = 0;
        double x = 0;
//        void collision(double const * Y, size_t ij, size_t iej){
//            double gJ=Y[idx_j(ij)+DynRadBlastWave::Q_SOL::iGamma];
//            mJ=M2*M0 + M0; eJ=Eint2*M0; gAdiJ=gammaAdi
//            gEJ=Gamma_ej; mEJ=M2_ej*M0_ej+M0_ej; eEJ=Eint2_ej*M0_ej; gAdiEJ=gammaAdi_ej
//            i_gM=Gamma_ej; i_mM=M2_ej*M0_ej+M2*M0; i_eM=Eint2_ej*M0_ej
//        }
        size_t i_restarts = 0;
        int n_tot_eqs = 0;
//        int n_j_bws = 0;
//        int n_ej_bws = 0;
    };
    Pars * p_pars;
    std::unique_ptr<logger> p_log;
public:
    EvolveODEsystem(std::unique_ptr<Magnetar> & p_mag,
                    std::unique_ptr<GRB> & p_grb,//std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                    std::unique_ptr<Ejecta> & p_ej,//std::vector<std::unique_ptr<CumulativeShell>> & p_ej,
//                    bool run_magnetar,
                    Array & t_grid,
//                  size_t n_shells_j, size_t n_shells_ej, size_t n_layers_j, size_t n_layers_ej,
                    const Integrators::METHODS integrator,
                    int loglevel
    ){
        p_pars = new Pars(p_mag, p_grb, p_ej, t_grid);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EvolveODEsystem");
//        p_pars->n_layers_j  = run_jet_bws ? n_layers_j : 0;
//        p_pars->n_layers_ej = run_ej_bws ? n_layers_ej : 0;
//        p_pars->n_shells_j  = run_jet_bws ? n_shells_j : 0;
//        p_pars->n_shells_ej = run_ej_bws ? n_shells_ej : 0;
        // ----

        // allocate memory for the IC and solution
//        p_pars->n_j_bws =
        p_pars->n_tot_eqs = (int)p_mag->getNeq() + (int)p_grb->getNeq() + (int)p_ej->getNeq();
        (*p_log)(LOG_INFO,AT) << " ODE will solve"
                                        << " N_mag="<<p_mag->getNeq()
                                        << " N_grb="<<p_grb->getNeq()
                                        << " N_ej="<<p_ej->getNeq()
                                        << " (total " << p_pars->n_tot_eqs << ") equations. \n";
        m_InitData = new double [ p_pars->n_tot_eqs ];
        m_CurSol   = new double [ p_pars->n_tot_eqs ];
        m_TmpSol   = new double [ p_pars->n_tot_eqs ]; // for shell collision
        // checks for the settings of BWs interaction prescription

        // chose the integrator for the system
        m_loglevel = loglevel;
        m_Method = integrator;
        p_Integrator = nullptr;
        switch (m_Method) {
            case Integrators::RK4 :
                p_Integrator = new IntegratorStatic<Integrators::rk4_integrator>
                        (Integrators::rk4_integrator(), loglevel);
                break;
            case Integrators::EULER:
                p_Integrator = new IntegratorStatic<Integrators::euler_integrator>
                        (Integrators::euler_integrator(), loglevel);
                break;
            case Integrators::MIDPOINT:
                p_Integrator = new IntegratorStatic<Integrators::midpoint_integrator>
                        (Integrators::midpoint_integrator(), loglevel);
                break;
            case Integrators::HEUN:
                p_Integrator = new IntegratorStatic<Integrators::heun3_integrator>
                        (Integrators::heun3_integrator(), loglevel);
                break;
            case Integrators::RALSTON:
                p_Integrator = new IntegratorStatic<Integrators::ralston3_integrator>
                        (Integrators::ralston3_integrator(), loglevel);
                break;
            case Integrators::DOP853:
                p_Integrator = new IntegratorStatic<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), loglevel);
                break;
            case Integrators::DOP853E:
                p_Integrator = new IntegratorEmbedded<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), loglevel);
                break;
        }
    }
    Pars *& getPars(){ return p_pars; }
    ~EvolveODEsystem(){
        delete [] m_InitData;
        delete [] m_CurSol;
        delete [] m_TmpSol;
        delete p_pars;
        delete p_Integrator;
    }

    void setRelativePositions(){

        auto & jet_bws = p_pars->p_grb->getBWs();
        auto & ej_bws = p_pars->p_ej->getShells();

        size_t i_j_l;// = which_layer_to_use;//jet_bws.size() - 1;

        for (size_t il=0; il<ej_bws.size(); il++){
            for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++){
                auto & bw = ej_bws[il]->getBW(ish);
                i_j_l = bw->getPars()->which_jet_layer_to_use;
                double ej_ctheta = bw->getPars()->ctheta0;
                double j_ctheta = jet_bws[i_j_l]->getPars()->ctheta0;
                if (ej_ctheta < j_ctheta){
                    bw->getPars()->is_within0 = true;
                    bw->getPars()->ijl = i_j_l;
                }
            }
        }
#if 0
        for (size_t iej = 0; iej < ej_bws.size(); iej++){
            i_j_l = ej_bws[iej]->getPars()->which_jet_layer_to_use;
            double ej_ctheta = ej_bws[iej]->getPars()->ctheta0;
            double j_ctheta = jet_bws[i_j_l]->getPars()->ctheta0;
            if (ej_ctheta < j_ctheta){
                ej_bws[iej]->getPars()->is_within0 = true;
                ej_bws[iej]->getPars()->ijl = i_j_l;
            }
        }
#endif
    }

    /// fill the 'm_InitData' - vector of initial conditions for all blast waves
    void setInitialConditions( const double tb0 ){
        size_t ii = 0;
        // ***************************************
        if (p_pars->p_magnetar->run_magnetar){
            auto & magnetar = p_pars->p_magnetar;
            magnetar->setInitConditions(m_InitData, ii);
            ii += magnetar->getNeq();
        }
        // ***************************************
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t il = 0; il < jet_bws.size(); il++) {
//                auto &j_pars = jet_bws[il]->getPars();
//            double beta0 = EQS::Beta(j_pars->Gamma0);
//            double R0    = j_pars->tb0 * beta0 * CGS::c;
//            jet_bws[il]->getDensIsm()->evaluateRhoDrhoDrDefault(R0, j_pars->ctheta0);
//                std::cout << "Initializing blast wave (jet) layer=" << il << " (ii=" << ii << ")" << "\n";
                jet_bws[il]->setInitConditions(m_InitData, ii);
                ii += jet_bws[il]->getNeq();
            }
            double jet_extend = jet_bws[jet_bws.size()-1]->getPars()->ctheta0; // UNUSED
        }
        if ((p_pars->p_grb->run_jet_bws)&&(p_pars->p_ej->run_ej_bws))
            setRelativePositions();
        /// get indexes of sorted (with respect to the radius) shells
        // ***************************************
        if (p_pars->p_ej->run_ej_bws) {
            auto & ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++){
                    ej_bws[il]->getBW(ish)->setInitConditions(m_InitData, ii);
                    ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
//            for (size_t il = 0; il < p_pars->n_ej_bws; il++) {
//                ej_bws[il]->setInitConditions(m_InitData, ii);
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        // ***************************************

        // ***************************************
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            m_TmpSol[i] = m_InitData[i]; /// will be used for ix-1 solutions for restarts.
            m_CurSol[i] = m_InitData[i];
        }
        // **************************************
        if ( !isThereATermination() ){
            (*p_log)(LOG_ERR,AT)  <<" termination at initialization. Evolution canceled\n";
            exit(1);
        }
        // **************************************
        if ( !isSolutionOk() ) {
            (*p_log)(LOG_ERR,AT)   << " Unphysical value in the initial data for evolution \n Exiting...";
            exit(1);
        }
        // ***************************************
        if (p_pars->p_ej->run_ej_bws) //  sort; get thickness; optical depth
            for (auto & cumShell : p_pars->p_ej->getShells()) {
                size_t _i, _j;
                cumShell->updateActiveShells();
                if (!cumShell->checkIfActiveShellsOrdered(m_InitData, _i, _j)){
                    (*p_log)(LOG_ERR,AT) << " shells are initially not radially ordered: "
                                                    << " shell idx="<<_i<<"overruns shell idx="<<_j<<"\n";
                    exit(1);
//                    cumShell->updateSortActiveShells(m_InitData);
                }
                cumShell->evalShellThickness(m_InitData);
                cumShell->evalShellOptDepth(m_InitData);
            }
        // **************************************
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            if (!std::isfinite(m_InitData[i])){
                (*p_log)(LOG_ERR,AT) << AT << "\n Nan in initial data: i="<<i<<" val="<<m_InitData[i]<<"\n";
//                exit(1);
            }
        }
        // **************************************
        p_Integrator->SetRHS( RHS_comb, p_pars->n_tot_eqs, p_pars );
        p_Integrator->SetInitCond( tb0, m_InitData );
        applyUnits();
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(0);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(0);
        // add electron properties (needed for synchron calculation)

//        addComputeForwardShockMicrophysics(0);

        is_initialized = true;

    }
    // evolve all blast waves
    void evolve( const double dx, const size_t ix ){
        if (!is_initialized){
            (*p_log)(LOG_ERR,AT)<<" bw set is not initialized. Cannot evolve\n";
            exit(1);
        }
        Timer timer;
        Array & t_grid = p_pars->t_grid;
        // eval relative geometry (which BW is behind what BW)
//        check_layer_relative_position( ix );
/// also, update the shell thickness and optical depth after the substap
        for (size_t il = 0; il < p_pars->p_ej->nlayers(); il++) {
            auto &cumShell = p_pars->p_ej->getShells()[il];
            cumShell->evalShellThickness(m_TmpSol);
            cumShell->evalShellOptDepth(m_TmpSol);
        }
        // solve the ODE system for x_i = x_{i-1} + dx
        p_Integrator->Integrate( dx );
        // extract the solution vector from the ODE solver
        p_Integrator->GetY( m_CurSol );

        /// check if there one of the layers is terminated
        if ( !isThereATermination() ){
            if(p_pars->i_restarts > 1){
                (*p_log)(LOG_ERR,AT)  << " second restart. Should not occure. Exiting...";
                exit(1);
            }
            (*p_log)(LOG_ERR,AT)  << " restarting iteration as there was a termination\n";
            evolve( dx, ix );
            p_pars->i_restarts += 1;
        }
        /// check if any of the kilonova blastwaves have collided
//        std::cout << m_CurSol[p_pars->p_ej->getShells()[0]->getBW(0)->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2]<<"\n";
        /// TODO DEBUG: this helps when collisions are extremenly close to ofset the trunning back a bit
        double col_prec_fac = 1e-11;
        if ((p_pars->p_ej->run_ej_bws) && (p_pars->p_ej->do_collision)){
            size_t n_unordered_layers = 0;
            std::vector<size_t> il_unordered{};
            double tcoll_min, rcoll_min;// = std::numeric_limits<double>::max();
            /// loop over all layers(cumShells) and see if there is a collision in any, and when
            bool are_shells_sorted = true;
            size_t _i, _j;
            for (size_t il = 0; il < p_pars->p_ej->nlayers(); il++){
                auto & cumShell = p_pars->p_ej->getShells()[il];
                /// update the active shells (indexes array)
                cumShell->updateActiveShells();
                /// check if active shells are in order
                are_shells_sorted = cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j);
                if (!are_shells_sorted){
                    n_unordered_layers+=1;
                    il_unordered.push_back(il);
                }
            }
            /// if the active shells are not in order, it means there was a collision; step back; collide; evolve
            if (n_unordered_layers > 0){
                are_shells_sorted = false;
                (*p_log)(LOG_INFO,AT) << "============================================================== \n";
                (*p_log)(LOG_INFO,AT) << "Shell collision detected at it="<<ix
                                                << " between t[i-1]="<<t_grid[ix-1]<< " and t[i]="<<t_grid[ix]
                                                << " in Nlayer="<<n_unordered_layers
                                                << " e.g., i="<<_i<<" and j="<<_j
                                                << " \n";
//                double target_dx = 0;
                double trunning = t_grid[ix - 1];
                size_t i_substeps = 0;
                size_t ish1, ish2;
                size_t _ish1, _ish2;
                double tcoll, rcoll;
                /// do substepping until all shells, staged for collision in this step has collided


                while (!are_shells_sorted) {
//                    std::cout << m_CurSol[p_pars->p_ej->getShells()[0]->getBW(0)->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iEint2]<<"\n";
                    i_substeps += 1;
                    /// check if there is anything left to collide
                    if (i_substeps > p_pars->p_ej->nshells()-2){
                        (*p_log)(LOG_ERR,AT) << "N substeps exceeded n_shells.\n";
                        exit(1);
                    }
                    (*p_log)(LOG_INFO,AT) << "Starting i="<<i_substeps<<" substep.\n";
                    /// set lists of shells for a each layer to collide simultaneously
                    size_t il_in_which_collision = 0;
                    std::vector<size_t> multi_collision{};
                    std::vector<size_t> multi_collision_ish1_arr{};
                    std::vector<size_t> multi_collision_ish2_arr{};
                    ish1 = 0; ish2 = 0;
                    /// locate which shell will collide first (in different layers it may happen simultaneously)
                    tcoll_min = std::numeric_limits<double>::max();
                    Vector tcolls{};
                    bool _is_actually_sorted = true;
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il) {
                        auto &cumShell = p_pars->p_ej->getShells()[il];
                        if (!cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j)) {
                            (*p_log)(LOG_INFO,AT)<< "\tLayer [il="<<il<<"] with active shells"
                                                           <<" N="<<cumShell->n_active_shells<<"/"<<cumShell->m_nshells
                                                           <<" e.g. (i="<<_i<<" j="<<_j<<" collision)" <<" \n";
                            tcoll = 0.; _ish1 = 0; _ish2 = 0; rcoll = 0;
                            /// evaluate time between trunning and t_grid[ix] at which shells collide
                            cumShell->evalWhichShellsCollideNext(_ish1, _ish2, tcoll, rcoll,
                                                                 trunning, t_grid[ix],
                                                                 m_TmpSol, m_CurSol);
                            tcolls.push_back(tcoll);
                            if (tcoll < tcoll_min) {
                                if (std::abs(tcoll - tcoll_min) < 1e-5*trunning) { // TODO ADD epsilon window within which shell collision is simultanous!
                                    multi_collision.push_back(il); // this can only be entered second time
                                    multi_collision_ish1_arr.push_back(_ish1);
                                    multi_collision_ish2_arr.push_back(_ish2);
                                }
                                else {
//                                    (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<il<<"] tcoll="<<tcoll<<"\n";
                                    tcoll_min = tcoll;
                                    rcoll_min = rcoll;
                                    il_in_which_collision = il;
                                    ish1 = _ish1;
                                    ish2 = _ish2;
                                }
                            }
                            _is_actually_sorted = false;
                        }
                    }
//                    tcoll_min = tcolls[ indexOfMinimumElement(tcolls)  ];
                    (*p_log)(LOG_INFO,AT)<<"Nearest collision is in"
                                         <<" il="<<il_in_which_collision<<"/"<<n_unordered_layers
                                         <<" (idx1="<<ish1<<", idx2="<<ish2<<") "
                                         <<" tcoll="<<tcoll_min
                                         <<" while timestep is=["<<trunning<<", " << t_grid[ix]<<"]"
                                         <<"\n";


                    if (_is_actually_sorted){
                        (*p_log)(LOG_ERR,AT)<<" _is_actually_sorted\n";
                        exit(1);
                    }
                    if (tcoll_min == std::numeric_limits<double>::max()){
                        (*p_log)(LOG_ERR,AT)<<" failed to fined tcoll_min\n";
                        exit(1);
                    }
                    if ((tcoll_min <= t_grid[ix-1]) || (tcoll_min >= t_grid[ix])){
                        (*p_log)(LOG_ERR, AT) << " tcoll_min is outside t[i-1], t[i] range; "
                                              << " tcol_min="<<tcoll_min
                                              <<" t_grid[ix-1]="<<t_grid[ix-1]
                                              <<" t_grid[ix]="<<t_grid[ix]
                                              <<" \n";
                        exit(1);
                    }

                    /// insert the first layer in which there is a collision to a vector of layer with a collision
                    multi_collision.insert(multi_collision.begin(), il_in_which_collision);
                    multi_collision_ish1_arr.insert(multi_collision_ish1_arr.begin(), ish1);
                    multi_collision_ish2_arr.insert(multi_collision_ish2_arr.begin(), ish2);
                    /// log info
//                    (*p_log)(LOG_INFO, AT) << "Nearest collision is at"
//                                           << " tcoll_min=" << tcoll_min
//                                           << " while timestep is=["<<trunning<<", " << t_grid[ix]<<"]"
//                                           << " t_to_coll="<<trunning-tcoll_min<< "\n";
                    if (trunning - tcoll_min == 0){
                        (*p_log)(LOG_ERR,AT)<<" tcoll=trunning Multiple collisions?\n";
                        exit(1);
                    }
                    if (multi_collision.size() > 1) {
                        (*p_log)(LOG_INFO, AT) << "Multi-collision at tcoll_min=" << tcoll_min
                                               << " N=" << multi_collision.size()
                                               << " collisions (in N layers) occurs simultaneously at\n";
                        (*p_log)(LOG_INFO, AT) << "Performing i=" << i_substeps
                                               << " substep from trunning=" << trunning
                                               << " to tcoll_min=" << tcoll_min << "\n";
                    }
                    if (multi_collision.empty()){
                        (*p_log)(LOG_ERR,AT)<<" no collisions?\n";
                        exit(1);
                    }
                    (*p_log)(LOG_INFO,AT)
                            <<"Trying to integrate to collision from "
                            <<"trunning="<<trunning<<" to tcoll="<<tcoll_min<<" Precision is set="<<col_prec_fac<<"\n";
                    /// advance the ODE to time of the collision. Use previous solution [ix-1] as a starting point
                    p_Integrator->Integrate(trunning, tcoll_min*(1.-col_prec_fac), m_TmpSol);
                    /// extract the solution vector from the ODE solver into 'm_CurSol'
                    p_Integrator->GetY(m_CurSol);

                    /// in case of 'simultaneous' collision, collide shells in each layer
                    for (size_t j = 0; j < multi_collision.size(); ++j) {
                        auto &cumShell = p_pars->p_ej->getShells()[multi_collision[j]];
                        /// collide shells within each layer
                        cumShell->collide(multi_collision_ish1_arr[j],multi_collision_ish2_arr[j],
                                         m_CurSol, rcoll_min);
                        /// update the active shells (as one was deleted)
                        cumShell->updateActiveShells();
                    }
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il) {
                        auto &cumShell = p_pars->p_ej->getShells()[il];
                        if (!cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j)){
                            (*p_log)(LOG_ERR, AT) << "\tLayer [il="<< cumShell->mylayer << "] Not ordered after Nsteps="
                                                  <<i_substeps<<" Collision e.g. (i=" << _i << " j=" << _j << "), "
                                                  << " (R[i]=" <<cumShell->getR(_i)
                                                  << " R[j]=" <<cumShell->getR(_j)
                                                  << ") dR[j-i]="<<cumShell->getR(_j)-cumShell->getR(_i)
                                                  << " \n";
                            exit(1);
                        }
                    }
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il){
                        auto &cumShell = p_pars->p_ej->getShells()[il];
                        for (size_t ish = 0; ish < cumShell->n_active_shells; ish++){
                            size_t _ieq0 = cumShell->getBW(ish)->getPars()->ii_eq;
                            size_t _ieq1 = cumShell->getBW(ish)->getNeq();
                            for (size_t ii = _ieq0; ii < _ieq0 + _ieq1; ++ii){
                                if (!std::isfinite(m_CurSol[ii])||(m_CurSol[ii]>1e60)){

                                    std::cerr << ii << " " << m_CurSol[ii] << "\n";
                                    (*p_log)(LOG_ERR,AT)<<"nans...\n";
                                    exit(1);
                                }
                            }
                        }
                    }
                    /// after shells collided put the outcome of the shell collision into ODE solver internal storage
                    p_Integrator->SetY(m_CurSol);
                    /// copy this solution into temporary storage for restarts ' m_TmpSol '
                    for (size_t i = 0; i < p_pars->n_tot_eqs; ++i)
                        m_TmpSol[i] = m_CurSol[i];
                    /// also, update the shell thickness and optical depth after the substap
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il) {
                        auto & cumShell = p_pars->p_ej->getShells()[il];
                        /// check if after the collision the shells are now in order (they should be...)
                        are_shells_sorted = cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j);
                        if (!are_shells_sorted) {
                            (*p_log)(LOG_ERR, AT) << "\tLayer [il="<< il << "] Not ordered after Nsteps="
                                                  <<i_substeps<<" Collision e.g. (i=" << _i << " j=" << _j << "), "
                                                  << " (R[i]=" <<cumShell->getR(_i)
                                                  << " R[j]=" <<cumShell->getR(_j)
                                                  << ") dR[j-i]="<<cumShell->getR(_j)-cumShell->getR(_i)
                                                  << " \n";
                            exit(1);
                        }
                    }
                    /// advance the target timestep:
                    trunning = tcoll_min*(1.-col_prec_fac);
                    if (trunning == t_grid[ix]){
                        exit(1);
                    }
                    if (trunning > t_grid[ix]) {
                        (*p_log)(LOG_ERR, AT) << " trunning=" << trunning << " is larger than t_grid[ix]="
                                              << t_grid[ix] << "\n";
                        exit(1);
                    }
                    /// update the shell thickness and optical depth after the substap
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); il++) {
                        auto &cumShell = p_pars->p_ej->getShells()[il];
                        cumShell->evalShellThickness(m_CurSol);
                        cumShell->evalShellOptDepth(m_CurSol);
                    }
                    (*p_log)(LOG_INFO,AT)
                        <<"Trying to integrate after collision from trunning="<<trunning<<" to t[ix]="<<t_grid[ix]<<"\n";
                    /// try to complete the timestep, itegrating from trunning to t_grid[ix]
                    p_Integrator->Integrate(trunning, t_grid[ix], m_CurSol);
                    /// extract the solution vector from the ODE solver into 'm_CurSol'
                    p_Integrator->GetY(m_CurSol);
                    /// check if now shells are ordered:
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il){
                        auto & cumShell = p_pars->p_ej->getShells()[il];
                        are_shells_sorted = cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j);
                        if (!are_shells_sorted){
                            (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<il<<"]"
                                                            <<" Collisions in a substep N="
                                                            << i_substeps
                                                            <<" e.g., (i="<<_i<<" j="<<_j<<")"
                                                            <<" Continuing substep loop.\n";
                            (*p_log)(LOG_INFO,AT) << "--------------------------------------------------- \n";
                            break;
                        }
                    }
                    if (are_shells_sorted){
                        (*p_log)(LOG_INFO,AT)<<"Substepping is complete successfully after "
                                             <<"N="<<i_substeps<< " t[i]="<<t_grid[ix]<<" reached.\n";
                        (*p_log)(LOG_INFO,AT) << "============================================================== \n";
                    }
                    for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il){
                        auto &cumShell = p_pars->p_ej->getShells()[il];
                        for (size_t ish = 0; ish < cumShell->n_active_shells; ish++){
                            size_t _ieq0 = cumShell->getBW(ish)->getPars()->ii_eq;
                            size_t _ieq1 = cumShell->getBW(ish)->getNeq();
                            for (size_t ii = _ieq0; ii < _ieq0 + _ieq1; ++ii){
                                if (!std::isfinite(m_CurSol[ii])||(m_CurSol[ii]>1e60)){

                                    std::cerr << ii << " " << m_CurSol[ii] << "\n";
                                    (*p_log)(LOG_ERR,AT)<<"nans...\n";
                                    exit(1);
                                }
                            }
                        }
                    }
                }
            }
        }
        /// --- Log main results of the integration
        double time_took = timer.checkPoint();
        if ( ix % 10 == 0 ) {
            if ((*p_log).getLogLevel()>LOG_WARN){
                auto & sstream = (*p_log)(LOG_INFO,AT);
                sstream << "it=" << ix << "/" << t_grid.size()
                        << " t=" << string_format("%.2e",t_grid[ix]) << " s, "
                        << string_format("%.2e",t_grid[ix]/CGS::day) << " d "
                        << "[Solver="<< string_format("%.2e",time_took)<<" s] ";
//                if (p_pars->p_grb->run_jet_bws){
//                    p_pars->p_grb->infoFastestLayer(ix, sstream);
//                }
                if (p_pars->p_ej->run_ej_bws){
                    p_pars->p_ej->infoFastestShell(ix, m_TmpSol, m_CurSol, sstream);
                }
                sstream << "\n";
            }


//        for(size_t i_ej = 0; i_ej < p_pars->p_bws_ej.size(); i_ej++){
//            auto & ej = p_pars->p_bws_ej[i_ej];
//            std::cout << " ish=" << ej->getPars()->ishell
//                      << " il="  << ej->getPars()->ishell
//                      << " G=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iGamma, ix))
//                      << " M2=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iM2, ix))
//                      << " R=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iR, ix))
//                      << " rho=" << string_format("%.2e", ej->getDensIsm()->m_rho)
//                      << " drho=" << string_format("%.2e", ej->getDensIsm()->m_drhodr)
//                      << " Grho=" << string_format("%.2e", ej->getDensIsm()->m_drhodr)
//                      << " Grel=" << string_format("%.2e", ej->getDensIsm()->m_GammaRel)
//                      << " dGrho=" << string_format("%.2e", ej->getDensIsm()->m_dGammaRhodR)
//                      << " dGrhodG=" << string_format("%.2e", ej->getDensIsm()->m_dGammaReldGamma)
//                      << " Grel=" << string_format("%.2e", ej->getDensIsm()->m_GammaRel)
//        }
//             std::cout << p_pars->p_ej[0]->getCurrentIndexes() << "\n";
//        (*p_log)(LOG_INFO,AT) << "it=" << ix << "/"
//                                            << t_grid.size() << " t=" << t_grid[ix] << ", " << t_grid[ix]/CGS::day << "\n";
        }
        /// save previous step in a temporary solution for restarts
        for (size_t i = 0; i < p_pars->n_tot_eqs; ++i)
            m_TmpSol[i] = m_CurSol[i];
        // apply units, e.g., energy is usually evolved in E/E0
        applyUnits();
        // --- Log main results of the integration


//        isThereACollisionBetweenKnBlastWaves( dx, ix );
        /// check if blast wave has fully expanded
        isThereLateralExpansionTermiantion();
        // check if there are no nans/unphysical vals in solution
        if ( !isSolutionOk() ) {
            (*p_log)(LOG_ERR,AT)  << " Unphysical value in the solution \n";
            exit(1);
        }
        p_pars->ix = ix;// latest solution
        p_pars->dx = dx;// latest solution
        p_pars->x = t_grid[ix];
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(ix);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(ix);
        // add electron properties (needed for synchron calculation)
//        addComputeForwardShockMicrophysics(ix);
        ///


//        (*p_log)(LOG_INFO, AT) << "Initialization finished [" << timer.checkPoint() << " s]" << "\n";


    }
    inline auto * pIntegrator() { return p_Integrator; }

private:
    void collisions(){

    }
    void applyUnits(){
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->applyUnits(m_CurSol, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->applyUnits(m_CurSol, ii);
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->applyUnits(m_CurSol, ii);
                    ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }
//        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
//            ej_bws[i]->applyUnits( m_CurSol, ii);
//            ii += p_pars->n_eq_ej_bws;
//        }
    }
    bool isThereATermination(){
//        auto & magnetar = p_pars->p_magnetar;

        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToTerminate(m_CurSol, ii)) {
                    is_ok = false;
                    (*p_log)(LOG_ERR,AT) << " Terminating jet BW layer=" << i << " (of all) failed "
                              << " [ishell=" << jet_bws[i]->getPars()->ishell
                              << " ilayer=" << jet_bws[i]->getPars()->ilayer
                              << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                              << " ] \n";
                    jet_bws[i]->getPars()->end_evolution = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();
            }
        }

        if (p_pars->p_ej->run_ej_bws) {
            auto & ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if (bw->isToTerminate(m_CurSol, ii)) {
                        is_ok = false;
                        bw->getPars()->end_evolution = true; // SET TO END
                        (*p_log)(LOG_ERR,AT) << " Terminating ejecta BW [ish=" << ish << " il=" << il << " "
                                  << " ii_eq=" << bw->getPars()->ii_eq
                                  << " ] \n";
                    }
                    ii += bw->getNeq();//p_pars->n_eq_ej_bws;
                }
            }


//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                if (ej_bws[i]->isToTerminate(m_CurSol, ii)) {
//                    is_ok = false;
//                    std::cerr  << " Terminating ejecta BW layer=" << i << " (of all) failed "
//                               << " [ishell=" << ej_bws[i]->getPars()->ishell
//                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
//                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
//                               << " ] \n";
//                    ej_bws[i]->getPars()->end_evolution = true; // SET TO END
//                }
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        return is_ok;
    }
    bool isThereLateralExpansionTermiantion(){
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    jet_bws[i]->getPars()->end_spreading = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        return is_ok;
    }
#if 0
    int isThereACollisionBetweenKnBlastWaves(double x, size_t ix){
        if (p_pars->p_ej->run_ej_bws){
            /// for one shell there is no collision possible
            if (p_pars->p_ej->nshells() == 1)
                return true;
//            std::cout << "nshells="<<p_pars->p_ej->nshells()<<"\n"; exit(1);
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                int order_info = ej_bws[il]->areShellsOrderedRadially( m_CurSol );
                (*p_log)(LOG_INFO,AT)<< " Shell Order Infor " << order_info << "\n";
            }
        }
        return true;
        //        size_t ii = 0; bool is_ok = true;
//        /// move the ii beyong GRB
//        if (p_pars->p_grb->run_jet_bws) {
//            auto & jet_bws = p_pars->p_grb->getBWs();
//            for (size_t i = 0; i < jet_bws.size(); i++) {
//                ii += jet_bws[i]->getNeq();
//            }
//        }

//        if (p_pars->p_ej->run_ej_bws) {
//            auto & ej_bws = p_pars->p_ej->getShells();
//            for (size_t il=0; il<ej_bws.size(); il++){
//                ej_bws[il]->setShellOrder(ix, m_InitData);
//            }
//        }


//        if (p_pars->p_ej->run_ej_bws) {
//            auto & ej_bws = p_pars->p_ej->getShells();
//            for (size_t il=0; il<ej_bws.size(); il++){
//                Vector rs (ej_bws[il]->nBWs(), 0.0);
//                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
//                    auto & bw = ej_bws[il]->getBW(ish);
//
//
//                    if (bw->isToTerminate(m_CurSol, ii)) {
//                        is_ok = false;
//                        bw->getPars()->end_evolution = true; // SET TO END
//                        (*p_log)(LOG_ERR,AT) << " Terminating ejecta BW [ish=" << ish << " il=" << il << " "
//                                             << " ii_eq=" << bw->getPars()->ii_eq
//                                             << " ] \n";
//                    }
//                    ii += bw->getNeq();//p_pars->n_eq_ej_bws;
//                }
//            }
//        }
    }
#endif
    bool isSolutionOk(){
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            if (!magnetar->isSolutionOk(m_CurSol, ii)){
                is_ok = false;
                (*p_log)(LOG_ERR,AT)  << " Magnetar evolution failed "
//                           << " [ishell=" << jet_bws[i]->getPars()->ishell
//                           << " ilayer=" << jet_bws[i]->getPars()->ilayer
//                           << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                           << " \n";
            }
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (( !jet_bws[i]->getPars()->end_evolution ) && (!jet_bws[i]->isSolutionOk(m_CurSol, ii))) {
                    is_ok = false;
                    (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of jet BW layer=" << i << " (of all) failed "
                               << " [ishell=" << jet_bws[i]->getPars()->ishell
                               << " ilayer=" << jet_bws[i]->getPars()->ilayer
                               << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                               << " ] \n";
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if ((!bw->getPars()->end_evolution) && (!bw->isSolutionOk(m_CurSol, ii))) {
                        is_ok = false;
                        (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of ejecta BW failed "
                                   << " [ishell=" << ish
                                   << " ilayer=" << il
                                   << " ii_eq=" << bw->getPars()->ii_eq
                                   << " ] \n";
                    }
                    ii += bw->getNeq();
                }
            }
//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                if ((!ej_bws[i]->getPars()->end_evolution) && (!ej_bws[i]->isSolutionOk(m_CurSol, ii))) {
//                    is_ok = false;
//                    std::cerr  << " Dyn. evol. of ejecta BW layer=" << i << " (of all) failed "
//                               << " [ishell=" << ej_bws[i]->getPars()->ishell
//                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
//                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
//                               << " ] \n";
//                }
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        return is_ok;
    }
    void insertSolution(size_t it){
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->insertSolution(m_CurSol, it, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->insertSolution(m_CurSol, it, ii);
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolution(m_CurSol, it, ii);
                    ii += bw->getNeq();
                }
            }
//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                ej_bws[i]->insertSolution(m_CurSol, it, ii);
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
    }
    void addOtherVariables(size_t it){
        if (p_pars->p_magnetar->run_magnetar) {
            p_pars->p_magnetar->addOtherVariables(it);
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->addOtherVars( it );
            }
        }
        auto & ej_bws = p_pars->p_ej;
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->addOtherVars(it);
                }
            }
        }


//        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
//            p_pars->p_bws_ej[i]->addOtherVars( it );
//        }
    }
    /// get offset of the equations in the combined solution vector for a given blast wave from its ilayer and ishell
    static void RHS_comb(double * out_Y,
                         size_t n_eq,
                         double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;

        /// advance magnetar to the next timestep
        size_t ii = 0;
//        double dEsddt, dEpropdt;
        double dEinj = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->evaluateRhs(out_Y, ii, x, Y);
//            dEsddt = p_pars->p_magnetar->getLdip(x, Y, ii);
//            dEpropdt = p_pars->p_magnetar->getLprop(x, Y, ii);
            dEinj = p_pars->p_magnetar->getValInt(Magnetar::Q::iLtot, x);
            ii += magnetar->getNeq();
        }
        else if (p_pars->p_magnetar->load_magnetar){
            dEinj = p_pars->p_magnetar->getValInt(Magnetar::Q::iLtot, x);
        }

        /// evaluate RHS for the jet (advance it to the next sub-step)
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); ++i) {
                jet_bws[i]->evaluateRhs(out_Y, ii, x, Y);
//                ii += p_pars->n_eq_j_bws;
                ii += jet_bws[i]->getNeq();//n_eq_j_bws;
            }
        }
//        else
//            std::vector<std::unique_ptr<BlastWaveBase>> jet_bws ;



        /// evaluate RHS for the ejecta (advance it to the next sub-step)
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_layers = p_pars->p_ej->getShells();
            for (size_t il=0; il < ej_layers.size(); il++){
//                ej_layers[il]->setShellOrder(x, Y);
//                ej_layers[il]->evalShellOptDepth(x, Y); // dEinjdt
                for(size_t ish=0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
//                    if (ej_bw->getPars()->end_evolution)
                    ej_bw->getPars()->dEinjdt = dEinj;
                    if (ej_bw->getPars()->dEinjdt < 0){
                        ej_bw->getPars()->dEinjdt = 0.;
                        (*ej_bw->getPars()->p_log)(LOG_ERR,AT) << " wrong value of dEinjdt="<<ej_bw->getPars()->dEinjdt<<"\n";
                    }
                    if (p_pars->p_grb->run_jet_bws) {
                        auto & jet_bws = p_pars->p_grb->getBWs();
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     & reinterpret_cast<std::vector<std::unique_ptr<BlastWaveBase>> &>(jet_bws),
                                                     p_pars->ix);
                    }
                    else
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     NULL,
                                                     p_pars->ix);
                    ii += ej_bw->getNeq();
                }
            }


//            for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++) {
////            if (ii!=ej_layers[iej]->getPars()->ii_eq){ // -------------- TO BE REMOVED
////                std::cerr<<AT << ii << " != "<< ej_layers[iej]->getPars()->ii_eq << "\n";
////                exit(1);
////            }
////            if (ej_layers[iej]->getPars()->end_evolution)
////                continue;
//                ej_layers[iej]->evaluateRhsDensModel2(out_Y, ii, x, Y,
//                     reinterpret_cast<std::vector<std::unique_ptr<BlastWaveBase>> &>(jet_bws), p_pars->ix);
////                ii += p_pars->n_eq_ej_bws;
//                ii += ej_layers[iej]->getNeq();
//            }
        }

        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************

    }
    /// TODO add method for BW collision
//    static void collision(const double * Y, std::unique_ptr<RadBlastWave> & jetPtr,
//                          std::unique_ptr<RadBlastWave> & ejPtr, void *pars){
//
//        auto * p_pars = (struct Pars *) pars;
//        // ii = il + nlayers * ish
//        double gJ = Y[idx(jetPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ]; //Gamma;
////        double mJ=M2*M0 + M0;  double eJ=Eint2*M0;  double gAdiJ=gammaAdi;
//        double gEJ= Y[idx(ejPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ];//Gamma_ej;
////        double mEJ=M2_ej*M0_ej+M0_ej; double eEJ=Eint2_ej*M0_ej; double gAdiEJ=gammaAdi_ej
////        double i_gM=Gamma_ej; double i_mM=M2_ej*M0_ej+M2*M0; double i_eM=Eint2_ej*M0_ej;
//        exit(1);
//    }
private:
    double * m_InitData{};
    double * m_CurSol{};
    double * m_TmpSol{};
    Integrators::METHODS m_Method{};
    int m_loglevel{};
    IntegratorBase * p_Integrator;
    bool is_initialized = false;
//    Array m_t_grid;
};

#endif //SRC_BLASTWAVE_SET_H
