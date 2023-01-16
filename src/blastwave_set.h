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

class SetDynRadBlastWaves{
    struct Pars{
        Pars(
                std::unique_ptr<Magnetar> & p_magnetar,
                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej,
                bool run_magnetar, bool run_jet_bws, bool run_ej_bws,
                Array & t_grid
        ) : p_magnetar(p_magnetar), p_bws_jet(p_bws_jet), p_bws_ej(p_bws_ej),
            run_magnetar(run_magnetar), run_jet_bws(run_jet_bws), run_ej_bws(run_ej_bws), t_grid(t_grid){
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
            if (run_magnetar){
                n_magnetar_eqs = p_magnetar->getNeq();
            }
            else{
                n_magnetar_eqs = 0;
            }
//            n_eq_j_bws = p_bws_jet[0]->getNeq();
//            n_eq_ej_bws= p_bws_ej[0]->getNeq();
            n_tot_eqs  = n_magnetar_eqs
                       + n_j_bws // number of 'jet' blast waves
                       * n_eq_j_bws // number of equations in each
                       + n_ej_bws // number of 'ejecta' blast waves
                       * n_eq_ej_bws; // number of equations in each;
            if (n_tot_eqs < 1){
                std::cerr << AT << "\n No equations to evolve; Both, jet and ejecta seems to be not given. Exiting...\n";
                exit(1);
            }
        }
        Array & t_grid;
        std::unique_ptr<Magnetar> & p_magnetar;
        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet;
        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej;
        bool run_magnetar{};
        bool run_jet_bws{};
        bool run_ej_bws{};
        size_t n_magnetar_eqs = 0;
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
    std::unique_ptr<logger> p_log;
public:
    SetDynRadBlastWaves(std::unique_ptr<Magnetar> & p_magnetar,
                        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_ej,
                        bool run_magnetar, bool run_jet_bws, bool run_ej_bws, Array & t_grid,
                        size_t n_shells_j, size_t n_shells_ej, size_t n_layers_j, size_t n_layers_ej,
                        const Integrators::METHODS integrator = Integrators::METHODS::RK4,
                        int loglevel = CurrLogLevel
    ){
        p_pars = new Pars(p_magnetar, p_bws_jet, p_bws_ej, run_magnetar, run_jet_bws, run_ej_bws, t_grid);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SetDynRadBlastWaves");
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
        auto & magnetar = p_pars->p_magnetar;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        // ***************************************
        if (p_pars->run_magnetar){
            magnetar->setInitConditions(m_InitData, ii);
        }
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
        if ((p_pars->n_layers_j > 0)&&(p_pars->n_layers_ej > 0))
            setRelativePositions();
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
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            if (!std::isfinite(m_InitData[i])){
                std::cerr << AT << "\n Nan in initial data: i="<<i<<" val="<<m_InitData[i]<<"\n";
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
//        addElectronAnalyticVars(0);
        is_initialized = true;

    }
    // evolve all blast waves
    void evolve( const double dx, const size_t ix ){
        if (!is_initialized){
            (*p_log)(LOG_ERR,AT)<<" bw set is not initialized. Cannot evolve\n";
            exit(1);
        }
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
        auto & magnetar = p_pars->p_magnetar;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0;
        if (p_pars->run_magnetar) {
            magnetar->applyUnits(m_CurSol, ii);
            ii += magnetar->getNeq();
        }
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
        auto & magnetar = p_pars->p_magnetar;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if(p_pars->run_magnetar){
            if (!magnetar->isSolutionOk(m_CurSol, ii)){
                is_ok = false;
                std::cerr  << " Magnetar evolution failed "
//                           << " [ishell=" << jet_bws[i]->getPars()->ishell
//                           << " ilayer=" << jet_bws[i]->getPars()->ilayer
//                           << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                           << " \n";
            }
            ii += magnetar->getNeq();
        }
        if(p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; i++) {
                if (( !jet_bws[i]->getPars()->end_evolution ) && (!jet_bws[i]->isSolutionOk(m_CurSol, ii))) {
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
        auto & magnetar = p_pars->p_magnetar;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0;
        if (p_pars->run_magnetar){
            magnetar->insertSolution(m_CurSol, it, ii);
        }
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
        if(p_pars->run_magnetar)
            p_pars->p_magnetar->addOtherVariables( it );
        for(size_t i = 0; i < p_pars->n_j_bws; i++){
            p_pars->p_bws_jet[i]->addOtherVars( it );
        }
        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
            p_pars->p_bws_ej[i]->addOtherVars( it );
        }
    }
    /// get offset of the equations in the combined solution vector for a given blast wave from its ilayer and ishell
    static inline size_t idx(std::unique_ptr<RadBlastWave> & bwPtr, void *pars){
        auto * p_pars = (struct Pars *) pars;
        return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + p_pars->n_layers_j * bwPtr->getPars()->ishell);
    }
    static void RHS_comb(double * out_Y,
                         size_t n_eq,
                         double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
        auto & magnetar = p_pars->p_magnetar;
        auto & jet_bws = p_pars->p_bws_jet;
        auto & ej_bws = p_pars->p_bws_ej;

        /// advance magnetar to the next timestep
        size_t ii = 0;
        if (p_pars->run_magnetar){
            magnetar->evaluateRhs(out_Y, ii, x, Y);
            ii += magnetar->getNeq();
        }

        /// evaluate RHS for the jet (advance it to the next sub-step)
        if (p_pars->run_jet_bws) {
            for (size_t i = 0; i < p_pars->n_j_bws; ++i) {
                jet_bws[i]->evaluateRhs(out_Y, ii, x, Y);
//                ii += p_pars->n_eq_j_bws;
                ii += jet_bws[i]->getNeq();//n_eq_j_bws;
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
                ej_bws[iej]->evaluateRhsDensModel2(out_Y, ii, x, Y,
                     reinterpret_cast<std::vector<std::unique_ptr<BlastWaveBase>> &>(jet_bws), p_pars->ix);
//                ii += p_pars->n_eq_ej_bws;
                ii += ej_bws[iej]->getNeq();
            }
        }

        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************

    }
    /// TODO add method for BW collision
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
