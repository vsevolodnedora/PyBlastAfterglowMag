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

#include "blastwave_base.h"
#include "blastwave_rad.h"
#include "blastwave_dyn.h"

class SetDynRadBlastWaves{
    struct Pars{
        Pars(
                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej,
                bool run_jet_bws, bool run_ej_bws,
                Array & t_grid
        ) : p_bws_jet(p_bws_jet), p_bws_ej(p_bws_ej),
            run_jet_bws(run_jet_bws), run_ej_bws(run_ej_bws), t_grid(t_grid){
            if (!p_bws_jet.empty()) {
                n_j_bws = p_bws_jet.size();
                n_eq_j_bws = p_bws_jet[0]->getNeq();
            }
            else {
                n_j_bws = 0;
                n_eq_j_bws = 0;
            }
            if(!p_bws_ej.empty()) {
                n_ej_bws = p_bws_ej.size();
                n_eq_ej_bws= p_bws_ej[0]->getNeq();
            }
            else {
                n_ej_bws = 0;
                n_eq_ej_bws = 0;
            }
//            n_eq_j_bws = p_bws_jet[0]->getNeq();
//            n_eq_ej_bws= p_bws_ej[0]->getNeq();
            n_tot_eqs  = n_j_bws // number of 'jet' blast waves
                         * n_eq_j_bws // number of equations in each
                         + n_ej_bws // number of 'ejecta' blast waves
                           * n_eq_ej_bws; // number of equations in each;
            if (n_tot_eqs < 1){
                std::cerr << AT << "\n No equations to evolve; Both, jet and ejecta seems to be not given. Exiting...\n";
                exit(1);
            }
        }
        Array & t_grid;
        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet;
        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej;
        bool run_jet_bws{};
        bool run_ej_bws{};
        size_t n_j_bws = 0;
        size_t n_ej_bws = 0;
        size_t n_eq_j_bws = 0;
        size_t n_eq_ej_bws = 0;
        size_t n_tot_eqs = 0;
        // ---
        size_t n_layers_j = 0;
        size_t n_layers_ej = 0;
        size_t n_shells_j = 0;
        size_t n_shells_ej = 0;
        // ---
//        bool use_cthetajet_lim_for_rho = false;
//        bool use_rjet_lim_for_rho = false;
        // ---
        bool is_entered = false;
        double r_ej_ent = -1;
        // ----

        inline size_t idx_j(const size_t i){
//            auto * p_pars = (struct Pars *) pars;
            auto & bwPtr = p_bws_jet[i];
            return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_j * bwPtr->getPars()->ishell);
        }
        inline size_t idx_ej(const size_t i){
//            auto * p_pars = (struct Pars *) pars;
            auto & bwPtr = p_bws_ej[i];
            size_t _idx_j = n_eq_j_bws * n_j_bws;
            return _idx_j + bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_ej * bwPtr->getPars()->ishell);
        }
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

    };
    Pars * p_pars;
//    std::unique_ptr<logger> p_log;
public:
    SetDynRadBlastWaves(std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej,
                        bool run_jet_bws, bool run_ej_bws, Array & t_grid,
                        size_t n_shells_j, size_t n_shells_ej, size_t n_layers_j, size_t n_layers_ej,
                        const Integrators::METHODS integrator = Integrators::METHODS::RK4,
                        int loglevel = CurrLogLevel
    ){
        p_pars = new Pars(p_bws_jet, p_bws_ej,  run_jet_bws, run_ej_bws, t_grid);
//        p_log = std::make_unique<logger>(std::cout, loglevel, "SetDynRadBlastWaves");
        p_pars->n_layers_j  = run_jet_bws ? n_layers_j : 0;
        p_pars->n_layers_ej = run_ej_bws ? n_layers_ej : 0;
        p_pars->n_shells_j  = run_jet_bws ? n_shells_j : 0;
        p_pars->n_shells_ej = run_ej_bws ? n_shells_ej : 0;
        // ----

        // allocate memory for the IC and solution
        m_InitData = new double [ p_pars->n_tot_eqs ];
        m_CurSol   = new double [ p_pars->n_tot_eqs ];
        // checks for the settings of BWs interaction prescription

        // chose the integrator for the system
        m_loglevel = loglevel;
        m_Method = integrator;
        p_Integrator = nullptr;
        switch (m_Method) {
            case Integrators::RK4 :
                p_Integrator = new IntegratorStatic<Integrators::rk4_integrator>
                        (Integrators::rk4_integrator(), CurrLogLevel);
                break;
            case Integrators::EULER:
                p_Integrator = new IntegratorStatic<Integrators::euler_integrator>
                        (Integrators::euler_integrator(), CurrLogLevel);
                break;
            case Integrators::MIDPOINT:
                p_Integrator = new IntegratorStatic<Integrators::midpoint_integrator>
                        (Integrators::midpoint_integrator(), CurrLogLevel);
                break;
            case Integrators::HEUN:
                p_Integrator = new IntegratorStatic<Integrators::heun3_integrator>
                        (Integrators::heun3_integrator(), CurrLogLevel);
                break;
            case Integrators::RALSTON:
                p_Integrator = new IntegratorStatic<Integrators::ralston3_integrator>
                        (Integrators::ralston3_integrator(), CurrLogLevel);
                break;
            case Integrators::DOP853:
                p_Integrator = new IntegratorStatic<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), CurrLogLevel);
                break;
            case Integrators::DOP853E:
                p_Integrator = new IntegratorEmbedded<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), CurrLogLevel);
                break;
        }
    }
    Pars *& getPars(){ return p_pars; }
    ~SetDynRadBlastWaves(){
        delete [] m_InitData;
        delete [] m_CurSol;
        delete p_pars;
        delete p_Integrator;
    }

    void setRelativePositions(){

        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;

        size_t i_j_l;// = which_layer_to_use;//jet_bws.size() - 1;

        for (size_t iej = 0; iej < ej_bws.size(); iej++){
            i_j_l = ej_bws[iej]->getPars()->which_jet_layer_to_use;
            double ej_ctheta = ej_bws[iej]->getPars()->ctheta0;
            double j_ctheta = jet_bws[i_j_l]->getPars()->ctheta0;
            if (ej_ctheta < j_ctheta){
                ej_bws[iej]->getPars()->is_within0 = true;
                ej_bws[iej]->getPars()->ijl = i_j_l;
            }
        }
    }

    /// fill the 'm_InitData' - vector of initial conditions for all blast waves
    void setInitialConditions( const double tb0 ){
        size_t ii = 0;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        // ***************************************
        if (p_pars->n_layers_j > 0) {
            for (size_t il = 0; il < p_pars->n_j_bws; il++) {
                auto &j_pars = jet_bws[il]->getPars();
//            double beta0 = EQS::Beta(j_pars->Gamma0);
//            double R0    = j_pars->tb0 * beta0 * CGS::c;
//            jet_bws[il]->getDensIsm()->evaluateRhoDrhoDrDefault(R0, j_pars->ctheta0);
//                std::cout << "Initializing blast wave (jet) layer=" << il << " (ii=" << ii << ")" << "\n";
                jet_bws[il]->setInitConditions(m_InitData, ii);
                ii += p_pars->n_eq_j_bws;
            }
            double jet_extend = jet_bws[jet_bws.size()-1]->getPars()->ctheta0; // UNUSED
        }
        if ((p_pars->n_layers_j > 0)&&(p_pars->n_layers_ej > 0)) setRelativePositions();
        // ***************************************
        if (p_pars->n_layers_ej > 0) {
            for (size_t il = 0; il < p_pars->n_ej_bws; il++) {
//            auto & ej_pars = ej_bws[il]->getPars();
//            if (ej_pars->ctheta0 < jet_extend){
//                double beta0 = EQS::Beta(ej_pars->Gamma0);
//                double R0    = ej_pars->tb0 * beta0 * CGS::c;
//                ej_bws[il]->getDensIsm()->evaluateRhoDrhoDr(R0, ej_pars->ctheta0);
//                ej_bws[il]->getDensIsm()->m_rho *= 1e-20;
//                ej_bws[il]->getDensIsm()->m_drhodr *= 0.;
//            }
//            else {
//                double beta0 = EQS::Beta(ej_bws[il]->getPars()->Gamma0);
//                double R0 = ej_bws[il]->getPars()->tb0 * beta0 * CGS::c;
//                ej_bws[il]->getDensIsm()->evaluateRhoDrhoDr(R0, ej_pars->ctheta0);
//            }
//                std::cout << "Initializing blast wave (ejecta) layer=" << il << " (ii=" << ii << ")" << "\n";
                ej_bws[il]->setInitConditions(m_InitData, ii);
                ii += p_pars->n_eq_ej_bws;
            }
        }
        // ***************************************
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){ m_CurSol[i] = m_InitData[i]; }
        // **************************************
        if ( !isThereATermination() ){
            std::cerr  << " termination at initialization\n Exiting...";
            std::cerr << AT << "\n ";
            exit(1);
        }
        // **************************************
        if ( !isSolutionOk() ) {
            // REMOVING LOGGER
            std::cerr   << " Unphysical value in the initial data for evolution \n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
        // **************************************
//        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
//            if (!std::isfinite(m_InitData[i])){
//                std::cerr << AT << "\n Nan in initial data: i="<<i<<" val="<<m_InitData[i]<<"\n";
////                exit(1);
//            }
//        }
        // **************************************
        p_Integrator->SetRHS( RHS_comb, p_pars->n_tot_eqs, p_pars );
        p_Integrator->SetInitCond( tb0, m_InitData );
        applyUnits();
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(0);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(0);
        // add electron properties (needed for synchron calculation)
//        addElectronAnalyticVars(0);
        is_initialized = true;

    }
    // evolve all blast waves
    void evolve( const double dx, const size_t ix ){
        // eval relative geometry (which BW is behind what BW)
//        check_layer_relative_position( ix );
        // solve the ODE system for x_i = x_{i-1} + dx
        p_Integrator->Integrate( dx );
        // extract the solution vector from the ODE solver
        p_Integrator->GetY( m_CurSol );
        // apply units, e.g., energy is usually evolved in E/E0
        applyUnits();
        // ---
        Array & t_grid = p_pars->t_grid;
        if ( ix % 10 == 0 ) {
            std::cout << "it=" << ix << "/" << t_grid.size() << " t=" << t_grid[ix] << "\n";
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

        // check if there one of the layers is terminated
        if ( !isThereATermination() ){
            if(p_pars->i_restarts > 1){
                std::cerr  << " second restart. Should not occure. Exiting...";
                std::cerr << AT << "\n";
                exit(1);
            }
            std::cerr  << " restarting iteration as there was a termination\n";
            evolve( dx, ix );
            p_pars->i_restarts += 1;
        }
        isThereLateralExpansionTermiantion();
        // check if there are no nans/unphysical vals in solution
        if ( !isSolutionOk() ) {
            // REMOVING LOGGER
            std::cerr  << " Unphysical value in the solution \n";
            std::cerr << AT << "\n ";
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
//        addElectronAnalyticVars(ix);
    }
    inline auto * pIntegrator() { return p_Integrator; }

private:
    void applyUnits(){
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0;
        for(size_t i = 0; i < p_pars->n_j_bws; i++){
            jet_bws[i]->applyUnits( m_CurSol, ii );
            ii += p_pars->n_eq_j_bws;
        }
        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
            ej_bws[i]->applyUnits( m_CurSol, ii);
            ii += p_pars->n_eq_ej_bws;
        }
    }
    bool isThereATermination(){
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; i++) {
                if (jet_bws[i]->isToTerminate(m_CurSol, ii)) {
                    is_ok = false;
                    std::cerr << " Terminating jet BW layer=" << i << " (of all) failed "
                              << " [ishell=" << jet_bws[i]->getPars()->ishell
                              << " ilayer=" << jet_bws[i]->getPars()->ilayer
                              << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                              << " ] \n";
                    jet_bws[i]->getPars()->end_evolution = true; // SET TO END
                }
                ii += p_pars->n_eq_j_bws;

            }
        }
        if (p_pars->run_ej_bws) {
            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
                if (ej_bws[i]->isToTerminate(m_CurSol, ii)) {
                    is_ok = false;
                    std::cerr  << " Terminating ejecta BW layer=" << i << " (of all) failed "
                               << " [ishell=" << ej_bws[i]->getPars()->ishell
                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
                               << " ] \n";
                    ej_bws[i]->getPars()->end_evolution = true; // SET TO END
                }
                ii += p_pars->n_eq_ej_bws;
            }
        }
        return is_ok;
    }
    bool isThereLateralExpansionTermiantion(){
        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; i++) {
                if (jet_bws[i]->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    jet_bws[i]->getPars()->end_spreading = true; // SET TO END
                }
                ii += p_pars->n_eq_j_bws;
            }
        }
        return is_ok;
    }
    bool isSolutionOk(){
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if(p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; i++) {
                if ((!jet_bws[i]->getPars()->end_evolution) && (!jet_bws[i]->isSolutionOk(m_CurSol, ii))) {
                    is_ok = false;
                    std::cerr  << " Dyn. evol. of jet BW layer=" << i << " (of all) failed "
                               << " [ishell=" << jet_bws[i]->getPars()->ishell
                               << " ilayer=" << jet_bws[i]->getPars()->ilayer
                               << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                               << " ] \n";
                }
                ii += p_pars->n_eq_j_bws;
            }
        }
        if (p_pars->run_ej_bws) {
            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
                if ((!ej_bws[i]->getPars()->end_evolution) && (!ej_bws[i]->isSolutionOk(m_CurSol, ii))) {
                    is_ok = false;
                    std::cerr  << " Dyn. evol. of ejecta BW layer=" << i << " (of all) failed "
                               << " [ishell=" << ej_bws[i]->getPars()->ishell
                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
                               << " ] \n";
                }
                ii += p_pars->n_eq_ej_bws;
            }
        }
        return is_ok;
    }
    void insertSolution(size_t it){
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0;
        if (p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; i++) {
                jet_bws[i]->insertSolution(m_CurSol, it, ii);
                ii += p_pars->n_eq_j_bws;
            }
        }
        if (p_pars->run_ej_bws) {
            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
                ej_bws[i]->insertSolution(m_CurSol, it, ii);
                ii += p_pars->n_eq_ej_bws;
            }
        }
    }
    void addOtherVariables(size_t it){
        for(size_t i = 0; i < p_pars->n_j_bws; i++){
            p_pars->p_bws_jet[i]->addOtherVars( it );
        }
        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
            p_pars->p_bws_ej[i]->addOtherVars(it );
        }
    }
    /// get offset of the equations in the combined solution vector for a given blast wave from its ilayer and ishell
    static inline size_t idx(std::unique_ptr<RadBlastWave> & bwPtr, void *pars){
        auto * p_pars = (struct Pars *) pars;
        return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + p_pars->n_layers_j * bwPtr->getPars()->ishell);
    }
    static void collision(const double * Y, std::unique_ptr<RadBlastWave> & jetPtr,
                          std::unique_ptr<RadBlastWave> & ejPtr, void *pars){

        auto * p_pars = (struct Pars *) pars;
        // ii = il + nlayers * ish
        double gJ = Y[idx(jetPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ]; //Gamma;
//        double mJ=M2*M0 + M0;  double eJ=Eint2*M0;  double gAdiJ=gammaAdi;
        double gEJ= Y[idx(ejPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ];//Gamma_ej;
//        double mEJ=M2_ej*M0_ej+M0_ej; double eEJ=Eint2_ej*M0_ej; double gAdiEJ=gammaAdi_ej
//        double i_gM=Gamma_ej; double i_mM=M2_ej*M0_ej+M2*M0; double i_eM=Eint2_ej*M0_ej;
        exit(1);
    }

    static void RHS_comb(double * out_Y,
                         size_t n_eq,
                         double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;

        /// evaluate RHS for the jet (advance it to the next sub-step)
        double x_grb = x;
        size_t ii = 0;
        if (p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; ++i) {
                jet_bws[i]->evaluateRhs(out_Y, ii, x, Y);
                ii += p_pars->n_eq_j_bws;
            }
        }
        if (p_pars->run_ej_bws){
            for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++) {
//            if (ii!=ej_bws[iej]->getPars()->ii_eq){ // -------------- TO BE REMOVED
//                std::cerr<<AT << ii << " != "<< ej_bws[iej]->getPars()->ii_eq << "\n";
//                exit(1);
//            }
//            if (ej_bws[iej]->getPars()->end_evolution)
//                continue;
                ej_bws[iej]->evaluateRhsDensModel2(out_Y, ii, x, Y, jet_bws, p_pars->ix);
                ii += p_pars->n_eq_ej_bws;
            }
        }
#if 0
        for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++){


            auto & ej_bw = ej_bws[iej];
            double ej_ctheta = ej_bw->ctheta( Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::itheta] );
            double ej_Gamma = Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::iGamma];
            double ej_R = Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::iR];

            size_t ijl = jet_bws.size()-1;;

            for (size_t ij = 0; ij < p_pars->n_j_bws; ++ij) {
                auto & j_bw = jet_bws[ij];
                size_t ii_j = p_pars->idx_j(ij);
                double j_theta = Y[ii_j + DynRadBlastWave::Q_SOL::itheta];
                double j_theta0 = j_bw->theta0(j_theta); //
                double j_theta1 = j_bw->theta1(j_theta);
                double j_R = Y[ii_j + DynRadBlastWave::Q_SOL::iR];

                /// should be entered for a SINGLE layer of a jet
                if ((ej_ctheta > j_theta0) && (ej_ctheta <= j_theta1)){
                    ijl = ij;
                    break;
                }
            }
//            if (ijl == 123456789){
//                exit(1);
//            }

            size_t idx_closest_jet= ijl;//jet_bws.size()-1;
            ej_bws[iej]->evaluateRhsDensModel(out_Y, ii, x, Y, jet_bws[idx_closest_jet],
                                            p_pars->idx_j(idx_closest_jet), p_pars->ix);
            ii += p_pars->n_eq_ej_bws;

        }
#endif
#if 0
        /// evaluate RHS for the ejecta (advance it to the next sub-step)
        for(size_t iej = 0; iej < p_pars->n_ej_bws; ++iej) {
            auto & ej_bw = ej_bws[iej];
            double ej_ctheta = ej_bw->ctheta( Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::itheta] );
            double ej_Gamma = Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::iGamma];
            double ej_R = Y[p_pars->idx_ej(iej) + DynRadBlastWave::Q_SOL::iR];

            bool is_inside = false;
            bool is_dens_evald = false;
            for (size_t ij = 0; ij < p_pars->n_j_bws; ++ij) {
                auto & j_bw = jet_bws[ij];
                double j_theta = Y[p_pars->idx_j(ij) + DynRadBlastWave::Q_SOL::itheta];
                double j_theta0 = j_bw->theta0(j_theta); //
                double j_theta1 = j_bw->theta1(j_theta);
                double j_R = Y[p_pars->idx_j(ij) + DynRadBlastWave::Q_SOL::iR];

                /// should be entered for a SINGLE layer of a jet
                if ((ej_ctheta > j_theta0) && (ej_ctheta <= j_theta1)){
                    /// check if alrady computed
                    if (is_dens_evald){
                        std::cerr << AT << " second evaluation of dens for the same jet BW\n";
                        exit(1);
                    }
                    /// at R0 (IC)
                    if (ej_R == ej_bw->getPars()->R0){
                        // inside -- NO ISM
                        ej_bw->getDensIsm()->evaluateRhoDrhoDr(ej_R, ej_ctheta);
                        ej_bw->getDensIsm()->m_rho = 1e-20 * ej_bw->getDensIsm()->m_rho;
                        ej_bw->getDensIsm()->m_drhodr = 0.;

                        ej_bw->getDensIsm()->m_GammaRho = 1;
                        ej_bw->getDensIsm()->m_GammaRel = ej_Gamma;
                        ej_bw->getDensIsm()->m_dGammaRhodR = 0;
                        ej_bw->getDensIsm()->m_dGammaReldGamma = 1;
                    }
                    /// behind jet
                    else if (ej_R <= j_R){
                        // first entry (entry radius) (r_entry < 0) || (r_entry > ej_R)
                        if ((ej_bw->getPars()->first_entry_r < 0) || (ej_bw->getPars()->first_entry_r > ej_R))
                            ej_bw->getPars()->first_entry_r = ej_R;
                        // if entering a new jet layer
                        if (ij != ej_bw->getPars()->ijl) {
                            ej_bw->getPars()->prev_ijl = ej_bw->getPars()->ijl; // -1 0 1 2 3 4
                            ej_bw->getPars()->ijl = ij; // 0 2 3 4 5
                        }
                        // set density parameters (no RHS calculation)
                        ij = p_pars->n_j_bws-1;
                        ej_bw->evalDensAndItsVelocityBehindBlastWave_Case1(ii, x, Y, jet_bws[ij],
                                                                           p_pars->idx_j(ij), p_pars->ix);
                        is_inside = true;
                    }
                    /// outside the jet
                    else{
                        // outside -- standard ISM
                        ej_bw->getDensIsm()->evaluateRhoDrhoDr(ej_R, ej_ctheta);

                        ej_bw->getDensIsm()->m_GammaRho = 1;
                        ej_bw->getDensIsm()->m_GammaRel = ej_Gamma;
                        ej_bw->getDensIsm()->m_dGammaRhodR = 0;
                        ej_bw->getDensIsm()->m_dGammaReldGamma = 1;
                    }
                    is_dens_evald = true;
                    break;
                }
            }
            /// this ejecta layer is outside the jet BWs
            if (!is_dens_evald){
                // outside -- standard ISM
                ej_bw->getDensIsm()->evaluateRhoDrhoDr(ej_R, ej_ctheta);

                ej_bw->getDensIsm()->m_GammaRho = 1;
                ej_bw->getDensIsm()->m_GammaRel = ej_Gamma;
                ej_bw->getDensIsm()->m_dGammaRhodR = 0;
                ej_bw->getDensIsm()->m_dGammaReldGamma = 1;
//                std::cerr << AT << " density is not evaluated\n";
//                exit(1);
            }
            // evaluate RHS itself
//            if (ej_Gamma <= 1.+1.e-6)
//                continue;
//            if (x == (double)p_pars->ix){
//                std::cerr << " failed " << '\n';
//                exit(1);
//            }

            double rho_i = ej_bw->getDensIsm()->m_rho;
            double drhodr_i = ej_bw->getDensIsm()->m_drhodr;

            double Gammarho_i = ej_bw->getDensIsm()->m_GammaRho;
            double GammaRel_i = ej_bw->getDensIsm()->m_GammaRel;
            double dGammaRhodR_i = ej_bw->getDensIsm()->m_dGammaRhodR;
            double dGammaReldGamma_i = ej_bw->getDensIsm()->m_dGammaReldGamma;

//            std::cout << " rho="<<rho_i
//                      << " drhodr="<<drhodr_i
//                      << " Gammarho="<<Gammarho_i
//                      << " GammaRel="<<GammaRel_i
//                      << " dGammaReldGamma="<<dGammaReldGamma_i
//                      <<"\n";
//            if (is_inside){
//                std::cerr << '\n';
//            }
            ej_bw->evaluateRhsDens(out_Y, ii, x, Y);

        }
#endif
#if 0
        //        for (size_t ij = 0; ij < p_pars->n_j_bws; ij++){
//            auto & j_bw = jet_bws[ij];
//            double j_theta  = Y[p_pars->idx_j(ij) + DynRadBlastWave::Q_SOL::itheta];
//            double j_theta0 = j_bw->theta0(j_theta); //
//            double j_theta1 = j_bw->theta1(j_theta);
//            double j_R      = Y[p_pars->idx_j(ij) + DynRadBlastWave::Q_SOL::iR];
////            std::cout << "ij="<<ij<<" theta=["<<j_theta0<<", "<<j_theta1<<"]" << '\n';
//
//            ///
//            for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++){
//                double ej_ctheta = Y[p_pars->idx_ej(ij) + DynRadBlastWave::Q_SOL::itheta];
//                double ej_Gamma = Y[p_pars->idx_ej(ij) + DynRadBlastWave::Q_SOL::iGamma];
//                double ej_R = Y[p_pars->idx_ej(ij) + DynRadBlastWave::Q_SOL::iR];
//                auto & ej_bw = ej_bws[iej];
//
//                // ejecta behind jet & inside [ vv ]
//                if ((ej_R < j_R)&&(ej_ctheta>j_theta0)&&(ej_ctheta<=j_theta1)){
//                    // first entry (entry radius)
//                    if ((ej_bw->getPars()->first_entry_r < -1) && (ej_bw->getPars()->first_entry_r > ej_R))
//                        ej_bw->getPars()->first_entry_r = ej_R;
//                    // if entering a new jet layer
//                    if (ij != ej_bw->getPars()->ijl){
//                        ej_bw->getPars()->prev_ijl = ej_bw->getPars()->ijl; // -1 0 1 2 3 4
//                        ej_bw->getPars()->ijl = ij; // 0 2 3 4 5
//                    }
//                    // set density parameters (no RHS calculation)
//                    ej_bw->evalDensAndItsVelocityBehindBlastWave_Case1(ii, x, Y, jet_bws[ij],
//                                                                       p_pars->idx_j(ij), p_pars->ix);
//
//                }
//                else{
//                    // set defaults for the ISM density
//                    ej_bw->getDensIsm()->evaluateRhoDrhoDr(ej_R, ej_ctheta);
////                ej_bw->getDensIsm()->m_rho = pfsolvePars->dens_floor_frac * p_dens->m_rho;
////                ej_bw->getDensIsm()->m_drhodr = 0.;
//
//                    ej_bw->getDensIsm()->m_GammaRho = 1;
//                    ej_bw->getDensIsm()->m_GammaRel = ej_Gamma;
//                    ej_bw->getDensIsm()->m_dGammaRhodR = 0;
//                    ej_bw->getDensIsm()->m_dGammaReldGamma = 1;
//
//                }
//                // evaluate RHS itself
//                ej_bw->evaluateRhsDens(out_Y, ii, x, Y);
//            }
//        }


//
//        /// create current grid of cthetas of the jet bws
//        Array cthetas_bound (p_pars->n_j_bws+2);
//        cthetas_bound[0] = 0.;
//        for (size_t ij = 0; ij < p_pars->n_j_bws; ij++){
//            auto & j_bw = jet_bws[ij];
//            double theta = Y[p_pars->idx_j(ij) + DynRadBlastWave::Q_SOL::itheta];
//            cthetas_bound[ij+1] = jet_bws[ij]->ctheta(theta);
//        }
//


        /// original system
//        for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++){
//
//            size_t idx_closest_jet= jet_bws.size()-1;
//            ej_bws[iej]->evaluateRhsDensModel(out_Y, ii, x, Y, jet_bws[idx_closest_jet],
//                                            p_pars->idx_j(idx_closest_jet), p_pars->ix);
//            ii += p_pars->n_eq_ej_bws;
//
//        }
#endif
#if 0

        // Place to add interaction between blast waves
        double j_rho, j_dlnrhodr;
        size_t ijl = jet_bws.size();
        double j_R = Y[iijl(jet_bws,ijl) + DynRadBlastWave::Q_SOL::iR];
        double j_theta = Y[iijl(jet_bws,ijl) + DynRadBlastWave::Q_SOL::itheta];
        double j_ctheta = jet_bws[ijl-1]->ctheta(j_theta);
        jet_bws[ijl-1]->getDensIsm()->getRhoDlnrhdR(j_rho, j_dlnrhodr, j_R, INFINITY);// get default
        double j_adi = jet_bws[ijl-1]->getEos()->getGammaAdi(
                Y[iijl(jet_bws,ijl) + DynRadBlastWave::Q_SOL::iGamma],
                EQS::Beta(Y[iijl(jet_bws,ijl) + DynRadBlastWave::Q_SOL::iGamma]));
        double j_rho2 = EQS::rho2t(Y[iijl(jet_bws,jet_bws.size()) + DynRadBlastWave::Q_SOL::iGamma],
                                   j_adi,
                                   j_rho);


        // ****************************************
        for(size_t i = 0; i < p_pars->n_ej_bws; i++){

            double ej_R = Y[ii + DynRadBlastWave::Q_SOL::iR];
            double ej_theta = Y[ii + DynRadBlastWave::Q_SOL::itheta];
            double ej_ctheta = ej_bws[i]->ctheta(ej_theta);

            // 3 cases:
            double rho_def, dlnrhodr_def, rho, dlnrhodr; // set default values for density
            ej_bws[i]->getDensIsm()->overrideRhoDlnRho(-1., -1.); // reset previous override
            ej_bws[i]->getDensIsm()->getRhoDlnrhdR(rho_def, dlnrhodr_def, ej_R, ej_ctheta);// get default

            if ( (j_ctheta > ej_ctheta) && (j_R > ej_R) ){
                // inside jet && behind head -> NO ISM (or exponential decay)
                if (ej_bws[i]->getPars()->entry_r < 0) {
                    ej_bws[i]->getPars()->entry_r = ej_R;
                    ej_bws[i]->getPars()->entry_time = x;
                }
//                if(ej_R < ej_bws[i]->getPars()->entry_r){
//                    ej_bws[i]->getPars()->entry_r = -1;
//                    std::cerr << " impossible";
//                    exit(1);
//                }

                if (ej_R != ej_bws[i]->getPars()->R0) { // it is 0 if not initialized (at it = 0)
                    // set density exponentially decaying from the point of entry
                    double coeff = 1;
                    double scale_rho = exp( 1. - coeff * ej_R / ej_bws[i]->getPars()->entry_r );
                    double scale_dlnrhodr = -1. * coeff * scale_rho / ej_bws[i]->getPars()->entry_r;

                    if (scale_rho > 1.2){
                        std::cout << ej_R << "\n";
                        std::cout << ej_bws[i]->getPars()->entry_r << "\n";
                        exit(1);
                    }

                    rho = rho_def * scale_rho;
                    dlnrhodr = rho_def * scale_dlnrhodr;// analytc. derivative
                }
                else{
                    rho = 0.;
                    dlnrhodr = 0.;
                }
                if (!std::isfinite(rho) || !std::isfinite(dlnrhodr) || rho < 0. || rho > 1. || rho > rho_def ){
                    std::cerr << AT << " wrong rho="<<rho<<" or dlnrhodr="<<dlnrhodr<<"\n";
                    exit(1);
                }
                if (rho < 1e-70){ rho = 1e-70; dlnrhodr = 0.; }
                // from Sedov taylor
                double rho_s = ej_bws[i]->getSedov()->rho_profile_int(ej_R, j_R, j_rho2);
                double rho_s2 = ej_bws[i]->getSedov()->rho_profile_int(ej_R - ej_R * 1e-4, j_R, j_rho2);
                double dlnrhodr_s = (1. / ej_R) * (rho_s - rho_s2) / (ej_R - ej_R * 1e-4);

                if(rho_s > rho) { rho = rho_s; dlnrhodr = dlnrhodr_s; }
                if(rho > 1.){ std::cerr << AT << " rho > 1" << "\n"; exit(1); }
                if(rho > j_rho2){
                    std::cerr << AT<< " rho > j_rho2"; exit(1);
                }

                ej_bws[i]->getDensIsm()->overrideRhoDlnRho(rho, dlnrhodr );
//                ej_bws[i]->getDensIsm()->overrideRhoDlnRho(1e-70, 0.);
            }
            else if ((j_ctheta > ej_ctheta) && (j_R < ej_R)) {
//                jet_bws[jet_bws.size()-1]->evaluateCollision(out_Y,
//                                                             p_pars->idx_j(jet_bws.size()-1),
//                                                             x, Y,
//                                                             ej_bws[i],
//                                                             p_pars->idx_ej(i));
                // inside the jet & above -> use default density
//                ej_bws[i]->getDensIsm()->overrideRhoDlnRho(rho, dlnrhodr);
//                p_pars->collision( Y, jet_bws.size()-1, i );
//                collision(Y,jet_bws[jet_bws.size()-1],ej_bws[i],pars);

            }
            else if((j_ctheta < ej_ctheta)) {
                // outside the jet cone -> use default value
//                ej_bws[i]->getDensIsm()->overridePars(rho, dlnrhodr);
            }
            else {
                std::cerr << " should not be entered\n";
                exit(1);
            }

//            for (size_t ival = 0; ival < p_pars->n_tot_eqs; ival ++){
//                if (!std::isfinite(Y[ival])){
//                    std::cerr << " nan in solution for ejecta\n";
//                    exit(1);
//                }
//            }

            ej_bws[i]->evaluateRhs( out_Y, ii, x, Y );
            ii += p_pars->n_eq_ej_bws;

//            // check for nans
//            for (size_t ival = 0; ival < p_pars->n_tot_eqs; ival ++){
//                if (!std::isfinite(out_Y[ival])){
//                    std::cerr << " nan in solution for ejecta\n";
//                    exit(1);
//                }
//            }

        }
#endif
        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************
    }
private:
    double * m_InitData{};
    double * m_CurSol{};
    Integrators::METHODS m_Method{};
    int m_loglevel{};
    IntegratorBase * p_Integrator;
    bool is_initialized = false;
    Array m_t_grid;
};

#endif //SRC_BLASTWAVE_SET_H
