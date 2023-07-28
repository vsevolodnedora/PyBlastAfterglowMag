//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_MODEL_EVOLVE_H
#define SRC_MODEL_EVOLVE_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"
#include "blastwave/blastwave_components.h"
#include "blastwave/blastwave_collision.h"
#include "model_magnetar.h"
#include "ejecta/ejecta.h"

class EvolveODEsystem{
    struct Pars{
        Pars(
                std::unique_ptr<Magnetar> & p_magnetar,
                std::unique_ptr<Ejecta> & p_grb,
                std::unique_ptr<Ejecta> & p_ej,
                std::unique_ptr<Ejecta> & p_ej_pwn2,
                std::unique_ptr<PWNset> & p_ej_pwn,
                Vector & t_grid,
                Vector & _t_grid
        ) :
        p_magnetar(p_magnetar), p_grb(p_grb), p_ej(p_ej), p_ej_pwn2(p_ej_pwn2), p_ej_pwn(p_ej_pwn),
            t_grid(t_grid), _t_grid(_t_grid){
        }
        Vector & t_grid;
        Vector & _t_grid;
        std::unique_ptr<Magnetar> & p_magnetar;
        std::unique_ptr<PWNset> & p_ej_pwn;
        std::unique_ptr<Ejecta> & p_grb;
        std::unique_ptr<Ejecta> & p_ej;
        std::unique_ptr<Ejecta> & p_ej_pwn2;
        // ---
//        size_t ix = 0;  // index latest solution
        double dx = 0;
        double x = 0;
        size_t prev_ix = 0;
        double prev_x = 0.;

        size_t i_restarts = 0;
        int n_tot_eqs = 0;
        size_t n_substeps = 0;
        Vector m_tmp_tsubstes{};
        bool do_average_solution = true;
    };
    Pars * p_pars;
    std::unique_ptr<logger> p_log;
public:
    EvolveODEsystem(std::unique_ptr<Magnetar> & p_mag,
                    std::unique_ptr<Ejecta> & p_grb,
                    std::unique_ptr<Ejecta> & p_ej,
                    std::unique_ptr<Ejecta> & p_ej_pwn2,
                    std::unique_ptr<PWNset> & p_ej_pwn,
                    Vector & t_grid,
                    Vector & _t_grid,
                    const Integrators::METHODS integrator,
                    bool do_average_solution,
                    int loglevel
    ){
        p_pars = new Pars(p_mag, p_grb, p_ej, p_ej_pwn2, p_ej_pwn, t_grid, _t_grid);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EvolveODEsystem");
        /// ---------------------------------------------------------------------------------------------
        p_pars->n_tot_eqs = (int)p_mag->getNeq() + (int)p_grb->getNeq() + (int)p_ej->getNeq() + (int)p_ej_pwn2->getNeq() + (int)p_ej_pwn->getNeq();
        (*p_log)(LOG_INFO,AT) << " ODE will solve"
                                        << " N_mag="<<p_mag->getNeq()
                                        << " N_grb="<<p_grb->getNeq()
                                        << " N_ej="<<p_ej->getNeq()
                                        << " N_ej_pwn2="<<p_ej_pwn2->getNeq()
                                        << " N_ej_pwn="<<p_ej_pwn->getNeq()
                                        << " (total " << p_pars->n_tot_eqs << ") equations. \n";
        p_pars->n_substeps = (size_t)_t_grid.size() / t_grid.size();
        p_pars->m_tmp_tsubstes.resize(p_pars->n_substeps, 0.);
        m_InitData = new double [ p_pars->n_tot_eqs ];
        m_CurSol   = new double [ p_pars->n_tot_eqs ];
        m_TmpSol   = new double [ p_pars->n_tot_eqs ]; // for shell collision
        /// init all arrays with zeroes
        for (size_t k = 0; k < p_pars->n_tot_eqs; k++){
            m_InitData[k] = 0.;m_CurSol[k]=0.;m_TmpSol[k]=0.;
        }

        m_CurSols  = new double * [ p_pars->n_substeps ];
        for (size_t j = 0; j < p_pars->n_substeps; ++j) // store result after each substep
            m_CurSols[j] = new double[p_pars->n_tot_eqs];
        // checks for the settings of BWs interaction prescription

        // chose the integrator for the system
        m_loglevel = loglevel;
        m_Method = integrator;
        p_Integrator = nullptr;
        p_pars->do_average_solution = do_average_solution;
        if (do_average_solution and p_pars->n_substeps == 1){
            (*p_log)(LOG_ERR,AT)<<" if do_average_solution = yes, iout must be > 1 \n";
            exit(1);
        }
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
        for (size_t i = 0; i < p_pars->n_substeps; ++i)
            delete[] m_CurSols[i];
        delete [] m_CurSols;
        delete p_pars;
        delete p_Integrator;
    }

    void setRelativePositions(){

        auto & ej_bws = p_pars->p_ej->getShells();

        if (p_pars->p_grb->nshells() > 1){
            (*p_log)(LOG_ERR,AT)<<"not implemented\n";
            exit(1);
        }

        auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();


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
        // ***********| M A G N E T A R |***********
        if (p_pars->p_magnetar->run_magnetar){
            auto & magnetar = p_pars->p_magnetar;
            magnetar->setInitConditions(m_InitData, ii);
            ii += magnetar->getNeq();
        }
#if 0
        // ***********| G R B |********
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
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
#endif
        // ***********| G R B |***********
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            (*p_log)(LOG_INFO,AT)<<"Initializing GRB blastwaves with "
                                <<" nlayers="<<ej_bws.size()<<" nshells="<<ej_bws[0]->nBWs()<<"\n";
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->setInitConditions(m_InitData, ii);
//                    ej_bws[il]->getBW(ish)->allocateSpaceForEachStepArrays(p_pars->_t_grid.size());
//                    ej_bws[il]->getBW(ish)->updateNucAtomic( m_InitData, tb0 );
                    ii += SOL::neq;//ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }

        // ***********| ????? |********
        /// get indexes of sorted (with respect to the radius) shells
        if ((p_pars->p_grb->run_bws)&&(p_pars->p_ej->run_bws))
            setRelativePositions();


        // ***********| E J E C T A |***********
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            (*p_log)(LOG_INFO,AT)<<"Initializing kn ejecta blastwaves with "
                                 <<" nlayers="<<ej_bws.size()<<" nshells="<<ej_bws[0]->nBWs()<<"\n";
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->setInitConditions(m_InitData, ii);
//                    ej_bws[il]->getBW(ish)->allocateSpaceForEachStepArrays(p_pars->_t_grid.size());
//                    ej_bws[il]->getBW(ish)->updateNucAtomic( m_InitData, tb0 );
                    ii += SOL::neq;//ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
            ejectaUpdate(tb0,0, m_InitData);
        }

        // ***********| E J E C T A    I N   E J E C T A |***********
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            (*p_log)(LOG_INFO,AT)<<"Initializing blastwaves inside kn ejecta with "
                                 <<" nlayers="<<ej_pwnbw.size()<<" nshells="<<ej_pwnbw[0]->nBWs()<<"\n";
            for (size_t il = 0; il < ej_pwnbw.size(); il++) {
                for (size_t ish = 0; ish < ej_pwnbw[il]->nBWs(); ish++) {
                    ej_pwnbw[il]->getBW(ish)->getPars()->curr_ldip = p_pars->p_magnetar->getMagValInt(MAG::ildip,tb0);
                    ej_pwnbw[il]->getBW(ish)->getPars()->curr_lacc = p_pars->p_magnetar->getMagValInt(MAG::ilacc,tb0);
                    ej_pwnbw[il]->getBW(ish)->setInitConditions(m_InitData, ii);
//                    ej_pwnbw[il]->getBW(ish)->allocateSpaceForEachStepArrays(p_pars->_t_grid.size());
//                    ej_bws[il]->getBW(ish)->updateNucAtomic( m_InitData, tb0 );
                    ii += SOL::neq;//ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }

        //************| S E T  I N I T I A L  V E C T O R |******
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            m_TmpSol[i] = m_InitData[i]; /// will be used for ix-1 solutions for restarts.
            m_CurSol[i] = m_InitData[i];
        }
        // **************************************
        if ( !isThereATermination(0) ){
            (*p_log)(LOG_ERR,AT)  <<" termination at initialization. Evolution canceled\n";
            exit(1);
        }
        // **************************************

        // **************************************
        if ( !isSolutionOk(tb0, 0) ) {
            (*p_log)(LOG_ERR,AT)   << " Unphysical value in the initial data for evolution \n Exiting...";
            exit(1);
        }
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
        // do the same for a smaller container of variables used in each stubstep
        insertSolutionSubstep(0, p_pars->_t_grid[0]);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(0);
        // add electron properties (needed for synchron calculation)
//        addComputeForwardShockMicrophysics(0);
        p_pars->prev_x = p_pars->_t_grid[0];
        p_pars->prev_ix = 0;
        is_initialized = true;

    }

    // evolve all blast waves
//    void evolve( const double dx, const size_t ix, const size_t i_x){ advanceTimeSubStep( dx, ix, i_x ); }

    void averageSolution(const size_t ix){
//        Vector m_tmp_time(p_pars->n_substeps);
        Vector m_tmp(p_pars->n_substeps);
        auto method = InterpBase::METHODS::iLagrangeUnivariate3m;
        /// compute the middle point between the start and end
        double t = p_pars->m_tmp_tsubstes[0] + (p_pars->m_tmp_tsubstes[p_pars->n_substeps-1] - p_pars->m_tmp_tsubstes[0]) / 2.;
        if (t < 0 || !std::isfinite(t)){
            (*p_log)(LOG_ERR,AT)<< "Interpolating between solutions: t="<<t<<"\n";
            exit(1);
        }
        auto intp = Interp1d(p_pars->m_tmp_tsubstes, m_tmp);
        if ((p_pars->do_average_solution) && (p_pars->n_substeps > 1)){
            for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
                for (size_t j = 0; j < p_pars->n_substeps; j++){
                    m_tmp[j] = m_CurSols[j][i];
                }
                double val = intp.Interpolate(t, method);
                if ((std::abs(val - m_CurSol[i]) > 1e-7) and
                    ((m_CurSol[i] < findMinimum(m_tmp)) || (m_CurSol[i] > findMaximum(m_tmp)))){
                    (*p_log)(LOG_ERR, AT) << " error in interpolating the solution: "
                        << " m_CurSol[i]="<<m_CurSol[i]<<" min="<<findMinimum(m_tmp)
                        << " max="<<findMaximum(m_tmp)<<"\n";
                    exit(1);
                }
                m_CurSol[i] = val;
            }
            /// overwrite the value of the time at which the solution is stored.
            p_pars->t_grid[ix] = t;
        }
    }

    void storeSolution(int ix){
        // average solution
        averageSolution(ix);
        // apply units, e.g., energy is usually evolved in E/E0
        applyUnits();
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(ix);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(ix);
        // add electron properties (needed for synchron calculation)
        //        addComputeForwardShockMicrophysics(ix);
        ///

    }

    inline auto * pIntegrator() { return p_Integrator; }
    void checkEvolution(){
        if (p_pars->p_magnetar->run_magnetar) {
//            p_pars->p_magnetar->addOtherVariables(it);
        }
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->checkEvolution();
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->checkEvolution();
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il = 0; il < ej_pwnbw.size(); il++) {
                for (size_t ish = 0; ish < ej_pwnbw[il]->nBWs(); ish++) {
                    ej_pwnbw[il]->getBW(ish)->checkEvolution();
                }
            }
        }
        if (p_pars->p_ej_pwn->run_pwn) {
            auto & pwns = p_pars->p_ej_pwn->getPWNs();
            for (auto & pwn : pwns){
//                pwn->addOtherVariables(it);
            }
        }
    }
    void insertSolutionSubstep(size_t it, double t){
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
//            magnetar->insertSolution(m_CurSol, it, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolutionSustep(m_CurSol, t, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolutionSustep(m_CurSol, t, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il=0; il<ej_pwnbw.size(); il++){
                for(size_t ish=0; ish<ej_pwnbw[il]->nBWs(); ish++) {
                    auto &bw = ej_pwnbw[il]->getBW(ish);
                    bw->insertSolutionSustep(m_CurSol, t, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn->run_pwn) {
            auto & pwns = p_pars->p_ej_pwn->getPWNs();
            for (auto & pwn : pwns){
//                pwn->insertSolutionS(m_CurSol, it, ii);
                ii += pwn->getNeq();
            }
        }
    }
    void advanceTimeSubStep(const double dx, const size_t ix, const size_t i_x){
        if (!is_initialized){
            (*p_log)(LOG_ERR,AT)<<" bw set is not initialized. Cannot evolve\n";
            exit(1);
        }
        Timer timer;
        Vector & t_grid = p_pars->_t_grid;

        /// update cumulative shell properties
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej->do_collision)) {
            for (auto &cumShell: p_pars->p_ej->getShells()) {
                size_t _i, _j;
                cumShell->updateActiveShells();
                if (!cumShell->checkIfActiveShellsOrdered(m_CurSol, _i, _j)) {
                    (*p_log)(LOG_ERR, AT) << " shells are not radially ordered: "
                                          << " shell idx=" << _i << "overruns shell idx=" << _j << "\n";
                    exit(1);
//                    cumShell->updateSortActiveShells(m_InitData);
                }
            }
        }

        /// solve the ODE system for x_i = x_{i-1} + dx
        p_Integrator->Integrate( dx );

        /// extract the solution vector from the ODE solver
        p_Integrator->GetY( m_CurSol );

        /// check if there one of the layers is terminated #TODO change this to not to use recursive call
        if ( !isThereATermination( ix ) ){
            if(p_pars->i_restarts > 1){
                (*p_log)(LOG_ERR,AT)  << " second restart. Should not occure. Exiting...";
                exit(1);
            }
            (*p_log)(LOG_ERR,AT)  << " restarting iteration as there was a termination\n";
            advanceTimeSubStep(dx, ix, i_x);
            p_pars->i_restarts += 1;
        }


        /// if blastwave collide, substuep through collision
        bool is_updated = false;
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej->do_collision) && (p_pars->p_ej->nMaxActiveShells() > 1))
            is_updated = doCollisionSubSteps( ix );

        if (p_pars->p_ej->run_bws && (!is_updated))
            ejectaUpdate(t_grid[ix], ix, m_CurSol);

//        bool is_laterly_sorted = true;
//        if ((p_pars->p_ej->run_bws)&&(true)&&(p_pars->p_ej->nMaxActiveShells())>1)
//            is_laterly_sorted = checkIfLateralShiftIsNeeded();

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
                if (p_pars->p_ej->run_bws){
                    p_pars->p_ej->infoFastestShell(ix, m_TmpSol, m_CurSol, sstream);
                }
                sstream << "\n";
            }


//        for(size_t i_ej = 0; i_ej < p_pars->p_bws.size(); i_ej++){
//            auto & ej = p_pars->p_bws[i_ej];
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
//             std::cout << p_pars->p_cumShells[0]->getCurrentIndexes() << "\n";
//        (*p_log)(LOG_INFO,AT) << "it=" << ix << "/"
//                                            << t_grid.size() << " t=" << t_grid[ix] << ", " << t_grid[ix]/CGS::day << "\n";
        }

        /// save previous step in a temporary solution for restarts
        for (size_t i = 0; i < p_pars->n_tot_eqs; ++i)
            m_TmpSol[i] = m_CurSol[i];

        /// store the solution for the timestep in a container for averaging
        p_pars->m_tmp_tsubstes[i_x] = t_grid[ix];
        for (size_t i = 0; i < p_pars->n_tot_eqs; ++i)
            m_CurSols[i_x][i] = m_CurSol[i];

//        // apply units, e.g., energy is usually evolved in E/E0
//        applyUnits();

        /// check if blast wave has fully expanded
        if (p_pars->p_grb->run_bws)
            isThereLateralExpansionTermiantion();
        /// check if blast wave has fully expanded
        if (p_pars->p_grb->run_bws)
            isThereReverseShockTermiantion();
        // check if there are no nans/unphysical vals in solution
        if ( !isSolutionOk( t_grid[ix], ix ) ) {
            (*p_log)(LOG_ERR,AT)  << " Unphysical value in the solution \n";
            exit(1);
        }
//        p_pars->ix = ix;// latest solution
//        p_pars->dx = dx;// latest solution
//        p_pars->x = t_grid[ix];

        p_pars->prev_x = t_grid[ix];
        p_pars->prev_ix = ix;

    }

private:
    void applyUnits(){
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->applyUnits(m_CurSol, ii);
            ii += magnetar->getNeq();
        }
#if 0
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->applyUnits(m_CurSol, ii);
                ii += jet_bws[i]->getNeq();
            }
        }
#endif
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->applyUnits(m_CurSol, ii);
                    ii += SOL::neq;;//ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->applyUnits(m_CurSol, ii);
                    ii += SOL::neq;//ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il = 0; il < ej_pwnbw.size(); il++) {
                for (size_t ish = 0; ish < ej_pwnbw[il]->nBWs(); ish++) {
                    ej_pwnbw[il]->getBW(ish)->applyUnits(m_CurSol, ii);
                    ii += SOL::neq;//ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }
    }
    bool isThereATermination(size_t ix){
//        auto & magnetar = p_pars->p_magnetar;
        size_t ii = 0; bool is_ok = true;
#if 0
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
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
#endif
        if (p_pars->p_grb->run_bws) {
            auto & ej_bws = p_pars->p_grb->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if (bw->isToTerminate(m_CurSol, ii, ix)) {
                        is_ok = false;
                        bw->getPars()->end_evolution = true; // SET TO END
                        (*p_log)(LOG_ERR,AT) << " Terminating grb BW [ish=" << ish << " il=" << il << " "
                                             << " ii_eq=" << bw->getPars()->ii_eq
                                             << " ] \n";
                    }
                    ii += SOL::neq;;//ii += bw->getNeq();//p_pars->n_eq_ej_bws;
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto & ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if (bw->isToTerminate(m_CurSol, ii, ix)) {
                        is_ok = false;
                        bw->getPars()->end_evolution = true; // SET TO END
                        (*p_log)(LOG_ERR,AT) << " Terminating ejecta BW [ish=" << ish << " il=" << il << " "
                                  << " ii_eq=" << bw->getPars()->ii_eq
                                  << " ] \n";
                    }
                    ii += SOL::neq;;//ii += bw->getNeq();//p_pars->n_eq_ej_bws;
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto & ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il=0; il<ej_pwnbw.size(); il++){
                for(size_t ish=0; ish<ej_pwnbw[il]->nBWs(); ish++) {
                    auto & bw = ej_pwnbw[il]->getBW(ish);
                    if (bw->isToTerminate(m_CurSol, ii, ix)) {
                        is_ok = false;
                        bw->getPars()->end_evolution = true; // SET TO END
                        (*p_log)(LOG_ERR,AT) << " Terminating ejecta BW [ish=" << ish << " il=" << il << " "
                                             << " ii_eq=" << bw->getPars()->ii_eq
                                             << " ] \n";
                    }
                    ii += SOL::neq;;//ii += bw->getNeq();//p_pars->n_eq_ej_bws;
                }
            }
        }
        return is_ok;
    }
    bool isThereLateralExpansionTermiantion(){
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_bws;
#if 0
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    jet_bws[i]->getPars()->end_spreading = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        return is_ok;
#endif
        size_t ii = 0; bool is_ok = true;
        auto & ej_bws = p_pars->p_grb->getShells();
        for (size_t il=0; il<ej_bws.size(); il++){
            for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                auto & bw = ej_bws[il]->getBW(ish);
                if (bw->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    bw->getPars()->end_spreading = true; // SET TO END
                }
                ii += SOL::neq;//ii += bw->getNeq();
            }
        }
        return is_ok;

    }
    bool isThereReverseShockTermiantion(){
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_bws;
#if 0
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    jet_bws[i]->getPars()->end_spreading = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        return is_ok;
#endif
        size_t ii = 0; bool is_ok = true;
        auto & ej_bws = p_pars->p_grb->getShells();
        for (size_t il=0; il<ej_bws.size(); il++){
            for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                auto & bw = ej_bws[il]->getBW(ish);
                if (bw->isToStopReverseShock(m_CurSol, ii)) {
                    is_ok = false;
                }
                ii += SOL::neq;//ii += bw->getNeq();
            }
        }
        return is_ok;
    }
    bool isSolutionOk(double x, size_t ix){
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_cumShells;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            if (!magnetar->isSolutionOk(m_CurSol, x, ix, ii)){
                is_ok = false;
                (*p_log)(LOG_ERR,AT)  << " Magnetar evolution failed "
//                           << " [ishell=" << jet_bws[i]->getPars()->ishell
//                           << " ilayer=" << jet_bws[i]->getPars()->ilayer
//                           << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                           << " \n";
            }
            ii += magnetar->getNeq();
        }
#if 0
        if (p_pars->p_grb->run_bws) {
            if (p_pars->p_grb->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            auto & jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//            auto & jet_bws = p_pars->p_grb->getBWs();
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
#endif
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    bw->getPars()->prev_x=x;
                    if ((!bw->getPars()->end_evolution) && (!bw->isSolutionOk(m_CurSol, x, ix, ii))) {
                        is_ok = false;
                        (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of GRB BW failed "
                                              << " [ishell=" << ish
                                              << " ilayer=" << il
                                              << " ii_eq=" << bw->getPars()->ii_eq
                                              << " ] \n";
                    }
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    bw->getPars()->prev_x=x;
                    if ((!bw->getPars()->end_evolution) && (!bw->isSolutionOk(m_CurSol, x, ix, ii))) {
                        is_ok = false;
                        (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of ejecta BW failed "
                                   << " [ishell=" << ish
                                   << " ilayer=" << il
                                   << " ii_eq=" << bw->getPars()->ii_eq
                                   << " ] \n";
                    }
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il=0; il<ej_pwnbw.size(); il++){
                for(size_t ish=0; ish<ej_pwnbw[il]->nBWs(); ish++) {
                    auto & bw = ej_pwnbw[il]->getBW(ish);
                    bw->getPars()->prev_x=x;
                    if ((!bw->getPars()->end_evolution) && (!bw->isSolutionOk(m_CurSol, x, ix, ii))) {
                        is_ok = false;
                        (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of ejecta PWN BW failed "
                                              << " [ishell=" << ish
                                              << " ilayer=" << il
                                              << " ii_eq=" << bw->getPars()->ii_eq
                                              << " ] \n";
                    }
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn->run_pwn) {
            auto & pwns = p_pars->p_ej_pwn->getPWNs();
            for (size_t il = 0; il < pwns.size(); ++il){
                if (!pwns[il]->isSolutionOk(m_CurSol, x, ix, ii)){
                    (*p_log)(LOG_ERR,AT)  << " PWN bound by ejecta (il="<<il<<") "
                                          //                           << " [ishell=" << jet_bws[i]->getPars()->ishell
                                          //                           << " ilayer=" << jet_bws[i]->getPars()->ilayer
                                          //                           << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                                          << " \n";
                    is_ok = false;
                }
                ii += pwns[il]->getNeq();
            }
        }

        return is_ok;
    }
    void insertSolution(size_t it){
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->insertSolution(m_CurSol, it, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolution(m_CurSol, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolution(m_CurSol, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_pwnbw = p_pars->p_ej_pwn2->getShells();
            for (size_t il=0; il<ej_pwnbw.size(); il++){
                for(size_t ish=0; ish<ej_pwnbw[il]->nBWs(); ish++) {
                    auto &bw = ej_pwnbw[il]->getBW(ish);
                    bw->insertSolution(m_CurSol, it, ii);
                    ii += SOL::neq;//ii += bw->getNeq();
                }
            }
        }
        if (p_pars->p_ej_pwn->run_pwn) {
            auto & pwns = p_pars->p_ej_pwn->getPWNs();
            for (auto & pwn : pwns){
                pwn->insertSolution(m_CurSol, it, ii);
                ii += pwn->getNeq();
            }
        }
    }
    void addOtherVariables(size_t it){
        if (p_pars->p_magnetar->run_magnetar) {
            p_pars->p_magnetar->addOtherVariables(it);
        }
        if (p_pars->p_grb->run_bws) {
            auto &ej_bws = p_pars->p_grb->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->addOtherVars(it);
                }
            }
        }
        if (p_pars->p_ej->run_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->addOtherVars(it);
                }
                /// save shell order; thickness; optical depth etc...
                if (p_pars->p_ej->do_collision)
                    ej_bws[il]->insertStatusInBWdata( it );
            }
        }
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_bws = p_pars->p_ej_pwn2->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->addOtherVars(it);
                }
            }
        }
        if (p_pars->p_ej_pwn->run_pwn) {
            auto & pwns = p_pars->p_ej_pwn->getPWNs();
            for (auto & pwn : pwns){
                pwn->addOtherVariables(it);
            }
        }
    }

// TODO this would be nice to have...
    bool checkIfLateralShiftIsNeeded(){
        for (size_t il = 0; il < p_pars->p_ej->nlayers()-1; il++){
//            auto & cumshell_jm1 = p_pars->p_ej->getShells()[il];
            auto & cumshell_j = p_pars->p_ej->getShells()[il];
            auto & cumshell_jp1 = p_pars->p_ej->getShells()[il+1];
            if (cumshell_j->getPars()->n_active_shells > 2 && cumshell_jp1->getPars()->n_active_shells > 2){
                if( (cumshell_j->getRvec()[0] > cumshell_jp1->getRvec()[0]) &&
                    (cumshell_j->getRvec()[0] > cumshell_jp1->getRvec()[1]) ){
                    double p_int = cumshell_j->getEthVec()[0] / cumshell_j->getVolVec()[0];
                    double lum = cumshell_j->getBWs()[cumshell_j->getIdx()[0]]->getPars()->dEinjdt;
                    double r = cumshell_j->getRvec()[0];
                    double p_rad = lum / (4. * CGS::pi * r * r * CGS::c);
                    std::cout << "p_rad/p_int = " << p_rad/p_int << "\n";
                    std::cout << "cumshell_j->getRvec()[0]="<<cumshell_j->getRvec()[0]
                        <<" cumshell_jp1->getRvec()[0]="<<cumshell_jp1->getRvec()[0]
                        <<" cumshell_jp1->getRvec()[1]="<<cumshell_jp1->getRvec()[1]
                        <<"\n";
                    int x = 1;
                    }
            }
        }
        return true;
    }

private:
    static void updateEnergyInjectionToEjectaBWs(double ldip, double lacc, double * Y, void * pars){
        auto * p_pars = (struct Pars *) pars;
        if (!p_pars->p_ej->run_bws) {
            return;
        }
        /// ********************| Ejecta Energy Injection |****************
        auto &ej_layers = p_pars->p_ej->getShells();
        for (size_t il = 0; il < ej_layers.size(); il++) {
            auto & pwn = p_pars->p_ej_pwn->getPWN(il);
            auto & cumShell = p_pars->p_ej->getShells()[il];
            double total_sd = 0.;
            if (p_pars->p_ej_pwn->run_pwn) {
                pwn->updateMagnetar(ldip, lacc);
                pwn->evalCurrBpwn(Y);
                if(!pwn->getPars()->is_init)
                    pwn->setInitConditions( Y, pwn->getPars()->iieq);
                /// apply eps_e and other coefficients from parfile
                total_sd = pwn->getAbsobedMagnetarLum(ldip, lacc, 1.);
            }
            /// update energy injection into blast waves
            for (size_t ish = 0; ish < cumShell->getPars()->n_active_shells; ish++) {
                auto & ej_bw = ej_layers[il]->getBW(ish);
                cumShell->getBW(ish)->getPars()->dEinjdt = 0.;
                if ((total_sd > 0.)) { // and(ish<10)
                    double fac_psr_dep_tmp = pwn->getFacPWNdep( // double rho_ej, double delta_ej, double T_ej, double Ye
                            cumShell->getRhoVec()[cumShell->getIdx()[ish]],
                            cumShell->getDeltaVec()[cumShell->getIdx()[ish]],
                            cumShell->getTempVec()[cumShell->getIdx()[ish]],
                            cumShell->getBW(ish)->getPars()->Ye0
                    );
                    ej_bw->getPars()->facPSRdep = fac_psr_dep_tmp;
                    ej_bw->getPars()->dEinjdt = fac_psr_dep_tmp * total_sd;
                    total_sd = total_sd - (fac_psr_dep_tmp * total_sd);
                    int x = 1;
                }
            }
        }
    }
    bool doCollisionSubSteps(size_t ix){

        double col_prec_fac = 1e-11; // TODO DEBUG: this helps when collisions are extremenly close to ofset the trunning back a bit
        Vector & t_grid = p_pars->_t_grid;
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
//                    std::cout << m_CurSol[p_pars->p_cumShells->getShells()[0]->getBW(0)->getPars()->ii_eq + DynRadBlastWave::QS::iEint2]<<"\n";
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
                                             <<" N="<<cumShell->getPars()->n_active_shells<<"/"
                                             <<cumShell->getPars()->nshells
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
#if 0
                /// update Magnetar energy injection
                if (!p_pars->p_cumShells->do_eninj_inside_rhs)
                    updateEnergyInjectionToEjectaBWs(
                            p_pars->p_magnetar->getMagValInt(Magnetar::Q::ildip, tcoll_min*(1.-col_prec_fac)),
                            p_pars->p_magnetar->getMagValInt(Magnetar::Q::ilacc, tcoll_min*(1.-col_prec_fac)),
                            m_CurSol,p_pars);

                /// update opacity and nuclear heating
                auto &ej_bws = p_pars->p_cumShells->getShells();
                for (size_t il = 0; il < ej_bws.size(); il++) {
                    for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                        ej_bws[il]->getBW(ish)->updateNucAtomic( m_CurSol, tcoll_min*(1.-col_prec_fac) );
                    }
                }
#endif
                /// advance the ODE to time of the collision. Use previous solution [ix-1] as a starting point
                p_Integrator->Integrate(trunning, tcoll_min*(1.-col_prec_fac), m_TmpSol);
                ejectaUpdate(tcoll_min*(1.-col_prec_fac), ix, m_TmpSol);
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
                        (*p_log)(LOG_ERR, AT) << "\tLayer [il="<< cumShell->getPars()->ilayer << "] Not ordered after Nsteps="
                                              <<i_substeps<<" Collision e.g. (i=" << _i << " j=" << _j << "), "
                                              << " (R[i]=" <<cumShell->getR(_i)
                                              << " R[j]=" <<cumShell->getR(_j)
                                              << ") dR[j-i]="<<cumShell->getR(_j)-cumShell->getR(_i)
                                              << " \n";
                        exit(1);
                    }
                }
#if 1
                /// check for nans anywhere
                for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il){
                    auto &cumShell = p_pars->p_ej->getShells()[il];
                    for (size_t ish = 0; ish < cumShell->getPars()->n_active_shells; ish++){
                        size_t _ieq0 = cumShell->getBW(ish)->getPars()->ii_eq;
                        size_t _ieq1 = SOL::neq;//cumShell->getBW(ish)->getNeq();
                        for (size_t ii = _ieq0; ii < _ieq0 + _ieq1; ++ii){
                            if (!std::isfinite(m_CurSol[ii])||(m_CurSol[ii]>1e60)){

                                std::cerr << ii << " " << m_CurSol[ii] << "\n";
                                (*p_log)(LOG_ERR,AT)<<"nans...\n";
                                exit(1);
                            }
                        }
                    }
                }
#endif
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
//                    updateEnergyInjection(m_CurSol, tcoll_min*(1.-col_prec_fac));

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
                    cumShell->updateSortedShellWidth(m_CurSol);
                    cumShell->updateSortedShellProperties(m_CurSol);
//                    cumShell->insertStatusInBWdata(ix);
                    /// update the shell properties for the PWN if needed
                    if (p_pars->p_ej_pwn->run_pwn){
                        auto & pwn = p_pars->p_ej_pwn->getPWN(il);
//                        pwn->updateOuterBoundary(cumShell->getRvec(),
//                                                 cumShell->getBetaVec(),
//                                                 cumShell->getRhoVec(),
//                                                 cumShell->getTauVec(),
//                                                 cumShell->getTempVec());
                        pwn->updateOuterBoundary(cumShell->getRvec()[0],
                                                 cumShell->getBetaVec()[0],
                                                 cumShell->getShellRho(m_CurSol),
                                                 cumShell->getShellOptDepth(),
                                                 cumShell->getTempVec()[0]);
                    }
                }
                (*p_log)(LOG_INFO,AT)
                        <<"Trying to integrate after collision from trunning="<<trunning<<" to t[ix]="<<t_grid[ix]<<"\n";
#if 0
                if (!p_pars->p_cumShells->do_eninj_inside_rhs)
                    updateEnergyInjectionToEjectaBWs(
                            p_pars->p_magnetar->getMagValInt(Magnetar::Q::ildip, t_grid[ix]),
                            p_pars->p_magnetar->getMagValInt(Magnetar::Q::ilacc, t_grid[ix]),
                            m_CurSol,p_pars);

                /// update opacity and nuclear heating
                for (size_t il = 0; il < p_pars->p_cumShells->m_nlayers(); il++) {
                    for (size_t ish = 0; ish < p_pars->p_cumShells->nshells(); ish++) {
                        auto & bw = p_pars->p_cumShells->getShells()[il]->getBW(ish);
                        bw->updateNucAtomic( m_CurSol, t_grid[ix] );
                    }
                }
#endif
                /// try to complete the timestep, itegrating from trunning to t_grid[ix]
                p_Integrator->Integrate(trunning, t_grid[ix], m_CurSol);
                ejectaUpdate(t_grid[ix], ix, m_CurSol);

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
#if 1
                /// check fo nans
                for (size_t il = 0; il < p_pars->p_ej->nlayers(); ++il){
                    auto &cumShell = p_pars->p_ej->getShells()[il];
                    for (size_t ish = 0; ish < cumShell->getPars()->n_active_shells; ish++){
                        size_t _ieq0 = cumShell->getBW(ish)->getPars()->ii_eq;
                        size_t _ieq1 = SOL::neq;//cumShell->getBW(ish)->getNeq();
                        for (size_t ii = _ieq0; ii < _ieq0 + _ieq1; ++ii){
                            if (!std::isfinite(m_CurSol[ii])||(m_CurSol[ii]>1e60)){

                                std::cerr << ii << " " << m_CurSol[ii] << "\n";
                                (*p_log)(LOG_ERR,AT)<<"nans...\n";
                                exit(1);
                            }
                        }
                    }
                }
#endif
            }
            return true;
        }
        else{
            return false;
        }
    }
    void ejectaUpdate(double time, size_t it, double * sol){

        /// update cumulative shell properties
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej->do_collision)) {
            for (auto &cumShell: p_pars->p_ej->getShells()) {
                size_t _i, _j;
                cumShell->updateActiveShells();
                if (!cumShell->checkIfActiveShellsOrdered(sol, _i, _j)) {
                    (*p_log)(LOG_ERR, AT) << " shells are not radially ordered: "
                                          << " shell idx=" << _i << "overruns shell idx=" << _j << "\n";
                    exit(1);
//                    cumShell->updateSortActiveShells(m_InitData);
                }
            }
        }

        /// update ejecta opacity and nuclear heating
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej->do_nuc)) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (auto &cumShell: p_pars->p_ej->getShells()) {
                for (auto & bw : cumShell->getBWs()) {
                    bw->updateNucAtomic( sol, time ); //t_grid[ix] );
                }
            }
        }

        /// update shell width, density, nuclear heating
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej->do_collision)) {
            for (auto &cumShell: p_pars->p_ej->getShells()) {
                cumShell->updateSortedShellWidth(sol);
                cumShell->updateSortedShellProperties(sol);
//                cumShell->insertStatusInBWdata(it);
            }
        }

        /// update outer boundary for PWN from ejecta (Shells SHOULD be updated)
        if ((p_pars->p_ej->run_bws) && (p_pars->p_ej_pwn->run_pwn)){
            if (!(p_pars->p_ej->do_collision)){
                (*p_log)(LOG_ERR,AT)<<" not supported\n";
                exit(1);
            }
            for (size_t il = 0; il < p_pars->p_ej->nlayers(); il++){
                auto & cumShell = p_pars->p_ej->getShells()[il];
                auto & ej_pwn = p_pars->p_ej_pwn->getPWNs()[il];
//                ej_pwn->updateOuterBoundary(
//                        cumShell->getRvec(),
//                        cumShell->getBetaVec(),
//                        cumShell->getRhoVec(),
//                        cumShell->getTauVec(),
//                        cumShell->getTempVec()
//                );
                ej_pwn->updateOuterBoundary(
                        cumShell->getRvec()[0],
                        cumShell->getBetaVec()[0],
                        cumShell->getShellRho(sol),
                        cumShell->getShellOptDepth(),
                        cumShell->getTempVec()[0]
                );
//                ej_pwn->evalCurrBpwn(sol);sss
            }
        }

        /// update outer boundary for PWN from ejecta (Shells SHOULD be updated)
        if (p_pars->p_ej_pwn2->run_bws){
            if (!(p_pars->p_ej_pwn2->do_collision) && (p_pars->p_ej_pwn2->nshells() > 1)){
                (*p_log)(LOG_ERR,AT)<<" not supported\n";
                exit(1);
            }
//            for (size_t il = 0; il < p_pars->p_ej->nlayers(); il++){
//                auto & cumShell = p_pars->p_ej->getShells()[il];
//                auto & ej_pwn = p_pars->p_ej_pwn->getPWNs()[il];
////                ej_pwn->updateOuterBoundary(
////                        cumShell->getRvec(),
////                        cumShell->getBetaVec(),
////                        cumShell->getRhoVec(),
////                        cumShell->getTauVec(),
////                        cumShell->getTempVec()
////                );
//                ej_pwn->updateOuterBoundary(
//                        cumShell->getRvec()[0],
//                        cumShell->getBetaVec()[0],
//                        cumShell->getShellRho(sol),
//                        cumShell->getShellOptDepth(),
//                        cumShell->getTempVec()[0]
//                );
////                ej_pwn->evalCurrBpwn(sol);sss
//            }
        }

        /// update Magnetar energy injecton into kN blast waves
        if (!p_pars->p_ej->do_eninj_inside_rhs) {
            double tinj = time;
            double ldip = p_pars->p_magnetar->getMagValInt(MAG::Q::ildip, tinj);
            double lacc = p_pars->p_magnetar->getMagValInt(MAG::Q::ilacc, tinj);
            updateEnergyInjectionToEjectaBWs(ldip, lacc, sol, p_pars);
        }

    }

    void prepareDensProfInFrontOfBW(double * out_Y, size_t ii, double x, double * Y){

    }
    /// evaluate RHS
#if 0
    static void RHS_comb_old(double * out_Y,
                         size_t n_eq,
                         const double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
        if (!std::isfinite(x)){
            std::cerr << AT <<" nan in tb!"<<"\n";
            exit(1);
        }
        /// *********************| M A G N E T A R |**********************
        size_t ii = 0;
        double ldip = 0;
        double lacc = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->evaluateRhsMag(out_Y, ii, x, Y);
//            ldip = p_pars->p_magnetar->(x, Y, ii);
//            lacc = p_pars->p_magnetar->getLprop(x, Y, ii);
//            ldip = p_pars->p_magnetar->getMagValInt(Magnetar::Q::ildip, x);
//            lacc = p_pars->p_magnetar->getMagValInt(Magnetar::Q::ilacc, x);
            ii += magnetar->getNeq();
        }
        else if (p_pars->p_magnetar->load_magnetar){
            ldip = p_pars->p_magnetar->getMagValInt(MAG::Q::ildip, x);
            lacc = p_pars->p_magnetar->getMagValInt(MAG::Q::ilacc, x);
        }

        /// ***************************| G R B |**************************
        if (p_pars->p_grb->run_bws) {
            auto &ej_layers = p_pars->p_grb->getShells();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                /// evaluate ejecta blast waves rhs
                for (size_t ish = 0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto &ej_bw = ej_layers[il]->getBW(ish);
                    ej_bw->evaluateRhs(out_Y, ii, x, Y);
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
        }

        /// ********************| Ejecta Energy Injection |****************
        if (p_pars->p_ej->do_eninj_inside_rhs)
            updateEnergyInjectionToEjectaBWs(ldip, lacc, const_cast<double *>(Y), pars);

        /// *************************| E J E C T A |***********************
        if (p_pars->p_ej->run_bws) {
            auto &ej_layers = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                /// evaluate ejecta blast waves rhs
                for(size_t ish=0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
                    if (p_pars->p_grb->run_bws) {
                        if (p_pars->p_grb->nshells() > 1){
                            std::cerr <<AT << " not implemented\n";
                            exit(1);
                        }
                        auto &jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
//                        auto & jet_bws = p_pars->p_grb->getBWs();
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
//                                                     & reinterpret_cast<std::vector<std::unique_ptr<BlastWave>> &>(jet_bws),
                                                     jet_bws,
                                                     p_pars->ix);
                    }
                    else {
//                        if (p_pars->use_dens_prof_behind_ejecta_for_ejecta){
//                            prepareDensProfInFrontOfBW(out_Y, ii, x, Y);
//                        }
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     std::vector<std::unique_ptr<BlastWave>>{},
                                                     p_pars->ix);
                    }
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
//            /// update density, velocity, int. energy of shells
//            if (p_pars->p_ej->do_shell_upd_in_rhs) {
//                for (size_t il = 0; il < ej_layers.size(); il++) {
//                    ej_layers[il]->updateActiveShells();
//                    ej_layers[il]->updateSortedShellWidth(Y);
//                    ej_layers[il]->updateSortedShellProperties(Y);
//                    for (size_t ish = 0; ish < ej_layers[il]->nBWs(); ish++) {
//                        auto &ej_bw = ej_layers[il]->getBW(ish);
//                        ej_bw->updateEnergyInjection(ldip, lacc);
//                    }
//                }
//            }
//            /// ****************************| P W N |***************************
//            if ((p_pars->p_ej_pwn->run_pwn)&&(p_pars->p_ej->do_eninj_inside_rhs)){
//                auto & pwn = p_pars->p_ej_pwn;
//                auto & ej_pwns = p_pars->p_ej_pwn->getPWNs();
//                auto & cumShells = p_pars->p_ej->getShells();
//                for (size_t il=0; il<pwn->nlayers(); il++){
//                    /// Set PWN ODE ICs
//                    ej_pwns[il]->updateMagnetar(ldip, lacc);
//                    ej_pwns[il]->evaluateRhs(out_Y, ii, x, Y);
//                    ii += ej_pwns[il]->getNeq();
//                }
//            }
        }

        /// *******************| E J E C T A  P W N |**********************
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_layers = p_pars->p_ej_pwn2->getShells();
            size_t nbws = ej_layers[0]->nBWs();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                /// evaluate ejecta blast waves rhs
                for(size_t ish=0; ish < nbws; ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
                    if (p_pars->p_ej->run_bws) {
                        if (p_pars->p_ej->nshells() < 2){
                            std::cerr <<AT << " not implemented\n";
                            exit(1);
                        }
                        /// for computing the upstram profile
                        auto &others = p_pars->p_ej->getShells()[il]->getBWs();
//                        auto & jet_bws = p_pars->p_grb->getBWs();
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     & static_cast<std::vector<std::unique_ptr<BlastWave>> &>(others),
                                                     p_pars->ix);
                    }
                    else {
//                        if (p_pars->use_dens_prof_behind_ejecta_for_ejecta){
//                            prepareDensProfInFrontOfBW(out_Y, ii, x, Y);
//                        }
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     NULL,
                                                     p_pars->ix);
                    }
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
        }

        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************

    }
#endif
    static void RHS_comb(double * out_Y,
                         size_t n_eq,
                         const double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
        if (!std::isfinite(x)){
            std::cerr << AT <<" nan in tb!"<<"\n";
            exit(1);
        }
        /// *********************| M A G N E T A R |**********************
        size_t ii = 0;
        double ldip = 0;
        double lacc = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->evaluateRhsMag(out_Y, ii, x, Y);
            ii += magnetar->getNeq();
        }
        else if (p_pars->p_magnetar->load_magnetar){
            ldip = p_pars->p_magnetar->getMagValInt(MAG::Q::ildip, x);
            lacc = p_pars->p_magnetar->getMagValInt(MAG::Q::ilacc, x);
        }

        /// ***************************| G R B |**************************
        if (p_pars->p_grb->run_bws) {
            auto &ej_layers = p_pars->p_grb->getShells();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                /// evaluate ejecta blast waves rhs
                for (size_t ish = 0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto &ej_bw = ej_layers[il]->getBW(ish);
                    if (ej_bw->getPars()->end_evolution) { ii += SOL::neq; continue; }
                    ej_bw->evaluateGRBRhs(out_Y, ii, x, Y);
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
        }

        /// ********************| Ejecta Energy Injection |****************
        if (p_pars->p_ej->do_eninj_inside_rhs)
            updateEnergyInjectionToEjectaBWs(ldip, lacc, const_cast<double *>(Y), pars);

        /// *************************| E J E C T A |***********************
        if (p_pars->p_ej->run_bws) {
            auto &ej_layers = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                for(size_t ish=0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
                    if (ej_bw->getPars()->end_evolution) { ii += SOL::neq; continue; }
                    if (p_pars->p_grb->run_bws) {
                        /// Run ejecta behind GRB accounting for GRB presence (pass GRB bws through)
                        if (p_pars->p_grb->nshells() > 1){ std::cerr <<AT << " not implemented\n"; exit(1); }
                        auto &jet_bws = p_pars->p_grb->getShells()[0]->getBWs();
                        ej_bw->prepareDensProfileFromJet(out_Y, ii, x, Y, jet_bws); // updates p_dens()
                        ej_bw->evaluateEjectaRhsDens(out_Y, ii, x, Y);
                    }
                    else {
                        /// Run ejecta without GRB presence
                        ej_bw->setStandardISM(out_Y, ii, x, Y); // updates p_dens()
                        ej_bw->evaluateEjectaRhsDens(out_Y, ii, x, Y);
                    }
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
        }

        /// *******************| E J E C T A  R A D  S H O C K |**********************
        if (p_pars->p_ej_pwn2->run_bws) {
            auto &ej_layers = p_pars->p_ej_pwn2->getShells();
            size_t nbws = ej_layers[0]->nBWs();
            for (size_t il = 0; il < ej_layers.size(); il++) {
                /// evaluate ejecta blast waves rhs
                for(size_t ish=0; ish < nbws; ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
                    if (ej_bw->getPars()->end_evolution) { ii += SOL::neq; continue; }
                    if (!p_pars->p_ej->run_bws) {
                        std::cerr <<AT << " for evolving BW inside ejecta, ejecta has to be evolved as well...\n";
                        exit(1);
                    }
                    if (p_pars->p_ej->nshells() < 2){
                        std::cerr <<AT << " for evolving BW inside ejecta, ejecta has to have > 2 shells \n";
                        exit(1);
                    }
                    auto &others = p_pars->p_ej->getShells()[il]->getBWs();
                    ej_bw->getPars()->curr_ldip = p_pars->p_magnetar->getMagValInt(MAG::ildip,x);
                    ej_bw->getPars()->curr_lacc = p_pars->p_magnetar->getMagValInt(MAG::ilacc,x);
                    ej_bw->evalDensProfileInsideBWset(out_Y, ii, x, Y, others); // update p_dens()
                    ej_bw->updateNucAtomic(Y,x);
                    ej_bw->updateCurrentBpwn(Y);
                    ej_bw->evaluateRhsDensPWN(out_Y, ii, x, Y);
                    ii += SOL::neq;//ii += ej_bw->getNeq();
                }
            }
        }



        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************

    }

private:
    double * m_InitData{};
    double * m_CurSol{};
    double * m_TmpSol{};
    double ** m_CurSols{};
    Integrators::METHODS m_Method{};
    int m_loglevel{};
    IntegratorBase * p_Integrator;
    bool is_initialized = false;
//    Array m_t_grid;
};

#endif //SRC_MODEL_EVOLVE_H
