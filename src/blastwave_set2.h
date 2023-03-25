//
// Created by vsevolod on 24/03/23.
//

#ifndef SRC_BLASTWAVE_SET2_H
#define SRC_BLASTWAVE_SET2_H



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
//#include "pulsar_wind_nebula.h"
#include "blastwave_base.h"
#include "blastwave_rad.h"
#include "blastwave_dyn.h"
#include "blastwave_set.h"


class EjectaLayer{
    struct Pars{
        Pars(Vector & _t_grid, std::unique_ptr<PWNmodel> & _p_pwn,
             std::unique_ptr<CumulativeShell> & _p_ej)
             : t_grid(_t_grid), p_pwn(_p_pwn), p_ej(_p_ej) { }
        Vector & t_grid;
        std::unique_ptr<PWNmodel> & p_pwn;
        std::unique_ptr<CumulativeShell> & p_ej;
        size_t ix = 0;  // index latest solution
        double x = 0;
        size_t i_restarts = 0;
        int n_tot_eqs = 0;
    };
    Pars * p_pars;
    std::unique_ptr<logger> p_log{};
    size_t m_ilayer{};
public:
    EjectaLayer( Vector & t_grid, std::unique_ptr<CumulativeShell> & p_ej,
                 std::unique_ptr<PWNmodel> & p_pwn,
                 const Integrators::METHODS integrator,
                 int ilayer,
                 int loglevel ){
        m_ilayer = ilayer;
        p_pars = new Pars(t_grid,p_pwn,p_ej);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EjectaLayer");
        p_pars->n_tot_eqs = (int)p_ej->getNeq() + (int)p_pwn->getNeq();
        (*p_log)(LOG_INFO,AT) << " ODE will solve"
                              << " N_mag="<<p_mag->getNeq()
                              << " N_grb="<<p_grb->getNeq()
                              << " N_ej="<<p_ej->getNeq()
                              << " N_ej_pwn="<<p_ej_pwn->getNeq()
                              << " (total " << p_pars->n_tot_eqs << ") equations. \n";
        m_InitData = new double [ p_pars->n_tot_eqs ];
        m_CurSol   = new double [ p_pars->n_tot_eqs ];
        m_TmpSol   = new double [ p_pars->n_tot_eqs ]; // for shell collision
    }



}

#endif //SRC_BLASTWAVE_SET2_H
