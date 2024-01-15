//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLAST_WAVE_COLLISION_H
#define SRC_BLAST_WAVE_COLLISION_H

//#include "../eats.h"

/// how to collide two blastwaves
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
public:
    BlastWaveCollision(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BlastWaveCollision");
        p_eos = new EOSadi();
        p_colsolve = new FsolvePars(p_eos);
    }

    ~BlastWaveCollision(){ delete p_colsolve; delete p_eos; }

    void collideBlastWaves(std::unique_ptr<BlastWave> & bw1,
                           std::unique_ptr<BlastWave> & bw2,
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
//        double gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + SOL::QS::imom] );
        double gam1 = Y[bw1->getPars()->ii_eq + SOL::QS::iGamma];
//        double beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + SOL::QS::imom] );
        double beta1 = Beta( Y[bw1->getPars()->ii_eq + SOL::QS::iGamma] );
        double m0_1 = bw1->getPars()->M0;
        double m2_1 = Y[bw1->getPars()->ii_eq + SOL::QS::iM2];
        double m2_1_ = m2_1 * m0_1;
        double eint2_1 = Y[bw1->getPars()->ii_eq + SOL::QS::iEint2];
        double eint2_1_ = eint2_1 * bw1->getPars()->M0 * CGS::c * CGS::c;
        double adi1 = bw1->getEos()->getGammaAdi(gam1, beta1);
        /// --------------------
//        double gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + SOL::QS::imom] );
        double gam2 = Y[bw2->getPars()->ii_eq + SOL::QS::iGamma];
//        double beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + SOL::QS::imom] );
        double beta2 = Beta( Y[bw2->getPars()->ii_eq + SOL::QS::iGamma] );
        double m0_2 = bw2->getPars()->M0;
        double m2_2 = Y[bw2->getPars()->ii_eq + SOL::QS::iM2];
        double m2_2_ = m2_2 * m0_2;
        double eint2_2 = Y[bw2->getPars()->ii_eq + SOL::QS::iEint2];
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
        double gam_c = i_gM;
        /// update the shell composition (mass averaged)
        double ye_c = (bw1->getPars()->Ye0 * bw1->getPars()->M0 + bw2->getPars()->Ye0 * bw2->getPars()->M0)
                      / (bw1->getPars()->M0 + bw2->getPars()->M0);
        if ((ye_c > bw1->getPars()->Ye0) && (ye_c > bw2->getPars()->Ye0)){
            (*p_log)(LOG_ERR,AT)<< "after collision ye="<<ye_c<<" while BWs had: 1 : "
                                << bw1->getPars()->Ye0<< " 2 : " << bw2->getPars()->Ye0 << "\n";
            exit(1);
        }
        /// chose which shell to update/delete
        int ish = p_colsolve->choseShell();
        double eint_before,eint_after,m2_before,m2_after,m0_before,m0_after;
        if (ish == 1){
            bw2->getPars()->end_evolution = true;
            bw1->getPars()->M0 = m0_c;
            bw1->getPars()->Ye0 = ye_c;
//            Y[bw1->getPars()->ii_eq + SOL::QS::imom] = mom_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::iGamma] = gam_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::iEint2] = eint2_c;
            Y[bw1->getPars()->ii_eq + SOL::QS::iM2] = m2_c;
        }
        else {
            bw1->getPars()->end_evolution = true;
            bw2->getPars()->M0 = m0_c;
            bw2->getPars()->Ye0 = ye_c;
//            Y[bw2->getPars()->ii_eq + SOL::QS::imom] = mom_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::iGamma] = gam_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::iEint2] = eint2_c;
            Y[bw2->getPars()->ii_eq + SOL::QS::iM2] = m2_c;
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
        p_colsolve->m_gam1 = EQS::GamFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom(Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        double m2_1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2]
                <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        double m0_1 = bw1->getPars()->M0;
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        p_colsolve->m_mass1 = m0_1 + m2_1;

        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom(Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        double m2_2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0;
        std::cout << " m2="<<Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2]
                  <<" m0="<<bw1->getPars()->M0<<" m2*m0="<<m2_1<<"\n";
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
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
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            m0_before = bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            //
            bw1->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw1->getPars()->M0;//Y[iieq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::QS::iM2];
            // terminated collided shell
            bw2->getPars()->end_evolution = true;
//            Y[iieq_other + DynRadBlastWave::QS::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itt] = 0;
            eint_after = i_eM / ( bw1->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::QS::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::QS::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw1->getPars()->M0;
        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            m0_before = bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0 + bw1->getPars()->M0;//i_mM; // update the total mass of the shell
//            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iM2] = (m0_1 + m2_1 + m2_2 + m0_2) / (m0_1 + m0_2);
            m2_after = Y[iieq + DynRadBlastWave::QS::iM2];
            // terminated collided shell
            bw1->getPars()->end_evolution = true;

//            Y[iieq_other + DynRadBlastWave::QS::iR] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iM2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEint2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEsh2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iErad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iEad2] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::iRsh] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itcomov] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itheta] = 0;
//            Y[iieq_other + DynRadBlastWave::QS::itt] = 0;
            eint_after = i_eM / ( bw2->getPars()->M0 * CGS::c * CGS::c );
            Y[iieq + DynRadBlastWave::QS::iEint2] = eint_after;
            Y[iieq + DynRadBlastWave::QS::imom] = i_gM * EQS::Beta(i_gM);
            m0_after = bw2->getPars()->M0;
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
//        double _mom = i_gM * EQS::Beta(i_gM);
//        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
//        Y[iieq + DynRadBlastWave::QS::imom] = _mom;
//        Y[iieq + DynRadBlastWave::QS::iEint2] = _eint2; // use ODE units

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

#endif //SRC_BLAST_WAVE_COLLISION_H
