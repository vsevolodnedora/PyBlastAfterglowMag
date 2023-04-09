//
// Created by vsevolod on 02/04/23.
//

#ifndef SRC_MODEL_EJECTA_H
#define SRC_MODEL_EJECTA_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "blastwave_components.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"

#include "model_magnetar.h"
#include "blastwave.h"

/// Radially structured blastwave collection
class CumulativeShell{
    struct Pars{
        size_t ilayer = 0;
        size_t nshells = 0;
        size_t n_active_shells = 0;
        double ltot = 0.;
        double tautot = 0.;
        double rphot = 0.;
        size_t idxphoto = 0;
        double tphoto = 0.;
        bool do_thermrad_loss = false;
        bool thermradloss_at_photo = false;
    };
    std::vector<size_t> m_idxs{};
    std::unique_ptr<logger> p_log = nullptr;
    std::unique_ptr<Pars> p_pars = nullptr;
    std::unique_ptr<BlastWaveCollision> p_coll = nullptr;
    std::vector<std::unique_ptr<BlastWave>> p_bws{};
    VecVector m_data{};
public:
    enum Q { ir, irho, ibeta, idelta, ivol, idtau, itaucum, itaucum0, ieth, itemp, ilum, itdiff };
    std::vector<std::string> m_vnames{
            "r", "rho", "beta", "delta", "vol", "dtau", "taucum", "taucum0", "eth", "temp", "lum", "tdiff"
    };
    CumulativeShell(Vector t_grid, size_t nshells, int ilayer, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CumulativeShell");
        p_coll = std::make_unique<BlastWaveCollision>(loglevel);
        p_pars = std::make_unique<Pars>();
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws.emplace_back(std::make_unique<BlastWave>(t_grid, ishell, ilayer, loglevel ) );
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize(nshells);
        m_idxs.resize(nshells);
        std::fill(m_data[Q::ir].begin(), m_data[Q::ir].end(), std::numeric_limits<double>::max());
        std::fill(m_idxs.begin(), m_idxs.end(), std::numeric_limits<size_t>::max());
        p_pars->ilayer=ilayer;
        p_pars->nshells=nshells;
        p_pars->n_active_shells=nshells;
    }
    void setPars(StrDbMap & pars, StrStrMap & opts){
        p_pars->do_thermrad_loss= getBoolOpt("do_thermrad_loss",opts,AT,p_log,false,true);
        if (!p_pars->do_thermrad_loss)
            return;
        p_pars->thermradloss_at_photo= getBoolOpt("thermradloss_at_photo",opts,AT,p_log,false,true);
    }
    // ------------------------
    inline double getR(size_t i){return m_data[Q::ir][i];}
    inline Vector & getRvec(){return m_data[Q::ir]; }
    inline std::vector<size_t> & getIdx(){ return m_idxs; }
    inline Vector & getBetaVec(){return m_data[Q::ibeta]; }
    inline Vector & getRhoVec(){return m_data[Q::irho]; }
    inline Vector & getTauVec(){return m_data[Q::idtau]; }
    inline Vector & getTempVec(){return m_data[Q::itemp]; }
    inline Vector & getDeltaVec(){return m_data[Q::idelta]; }
    // -----------------------
    std::unique_ptr<BlastWave> & getBW(size_t ish){
        if (p_bws.empty()){
            (*p_log)(LOG_ERR, AT) << " shell does not contain blast waves\n";
            exit(1);
        }
        if (ish > p_bws.size()){
            (*p_log)(LOG_ERR, AT) << "invalid memory accessed\n"; exit(1);
        }
        return p_bws[ish];
    }
    std::unique_ptr<Pars> & getPars(){ return p_pars; }
    inline size_t nBWs() const { return p_bws.size();}
    inline std::vector<std::unique_ptr<BlastWave>> & getBWs() {
        return p_bws;
    }
    // -----------------------
    /// update the number of active shells
    void updateActiveShells(){
        std::fill(m_idxs.begin(), m_idxs.end(),std::numeric_limits<size_t>::max());
        p_pars->n_active_shells = p_pars->nshells;
        size_t idx = 0;
        for (size_t i=0; i<p_pars->n_active_shells; i++) {
            /// loop over all blastwaves and collect the active ones
            if (p_bws[i]->getPars()->end_evolution) {
                continue; }
            /// if shell is active record its current index (order)
            m_idxs[idx] = i; // add only active shells to the list
            idx++;
        }
        p_pars->n_active_shells = idx;
    }
    /// check if active shells are ordered according to their radius
    bool checkIfActiveShellsOrdered(const double * Y, size_t & idx0, size_t & idx1 ){
        idx0 = std::numeric_limits<size_t>::max();
        idx1 = std::numeric_limits<size_t>::max();
        std::fill(m_data[Q::ir].begin(), m_data[Q::ir].end(),std::numeric_limits<double>::max());
        bool is_sorted = true;
        for (size_t i=0; i<p_pars->n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            m_data[Q::ir][i] = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR ];
        }
        for (size_t i=1; i<p_pars->n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            double rim1 = m_data[Q::ir][i-1];
            /// if the next shells in the list has a radius > than the previous; shells not ordered
            if (m_data[Q::ir][i] > rim1){
                rim1 = m_data[Q::ir][i];
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
        auto & bw1 = p_bws[idx1];
        auto & bw2 = p_bws[idx2];
        p_coll->collideBlastWaves(bw1, bw2, Y, rcoll, p_pars->ilayer);
    }
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
                << "\tLayer [il="<<p_pars->ilayer<<"] Checking which shells will collide in"
                <<" between tim1="<<tim1<<" and ti="<<ti<<"\n";
        size_t n_collisions = 0;
        for (size_t ii=0; ii<p_pars->n_active_shells; ++ii) {
            n_collisions = 0;
            size_t i_idx = m_idxs[ii];
            double r_i = m_data[Q::ir][ii];
            if (r_i == std::numeric_limits<double>::max()){
                continue;
            }
            std::vector<size_t> overrun_shells;
            Vector t_at_which_overrun;
            /// now check which blast-waves has this one overrun

            for (size_t jj = ii; jj < p_pars->n_active_shells; ++jj) {
                size_t j_idx = m_idxs[jj];
                double r_j = m_data[Q::ir][jj];
                if (r_i > r_j) {
                    n_collisions += 1;
                    overrun_shells.push_back(j_idx);
                    /// interpolate the time at which shells collided
//                    double tim1 = m_t_grid[ix - 1];
//                    double ti = m_t_grid[ix];
                    double dt = ti - tim1;
                    double r0im1 = Ym1_[p_bws[i_idx]->getPars()->ii_eq + SOL::QS::iR];//p_bws[i]->getVal(BW::iR, (int) (ix - 1));
                    double r0i   = Y_[p_bws[i_idx]->getPars()->ii_eq + SOL::QS::iR];
                    double r1im1 = Ym1_[p_bws[j_idx]->getPars()->ii_eq + SOL::QS::iR];//p_bws[j]->getVal(BW::iR, (int) (ix - 1));
                    double r1i   = Y_[p_bws[j_idx]->getPars()->ii_eq + SOL::QS::iR];
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
//                            << " mom=" << Y[p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom]
//                            << " (mom0=" << p_bws[i]->getPars()->mom0 << ") "
//                            << " has overrun shell j=" << overrun_shells[0]
//                            << " with momenta "
//                            << Y[p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom]
//                            << " (mom0=" << p_bws[overrun_shells[0]]->getPars()->mom0 << ") "
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
        (*p_log)(LOG_INFO, AT) << "\tLayer [il="<<p_pars->ilayer<<"] interpolated tcoll="<<tcoll
                               <<" (idx1="<<ish1<<" idx2="<<ish2<<")\n";
        if ((tcoll == std::numeric_limits<double>::max())||(!std::isfinite(tcoll))){
            (*p_log)(LOG_ERR,AT)<<" Failed to find tcoll in layer="<<p_pars->ilayer<<"\n";
            exit(1);
        }
    }
    /// Evaluate shell thickness using various methods. Assumes sorted shells (no overruns)
    void updateSortedShellWidth( const double * Y ){
        /// if there is one shell we cannot have the shell width that comes from shell separation.
        if (p_pars->n_active_shells == 1){
//            double r_i=0., r_ip1=0., dr_i = 0., vol_i=0.;
            size_t idx = 0;
            double frac = 0.1; // TODO fraction of the shell volume to be used as its width. !!! Inaccurate as we need adjacent shells to get the volume...
            auto & bw = p_bws[idx];
            double r_i =  Y[bw->getPars()->ii_eq + SOL::QS::iR];
            double dr_i = frac * r_i;
            double vol_i = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
            if (dr_i < 0 or !std::isfinite(dr_i)){
                (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
                exit(1);
            }
            m_data[Q::idelta][idx] = dr_i;
            m_data[Q::ivol][idx] = vol_i;
        }
        else{
            for (size_t ii=0; ii<p_pars->n_active_shells-1; ii++) {
                ///
                size_t idx = m_idxs[ii];
                size_t nextidx = m_idxs[ii + 1];
                auto &bw = p_bws[idx];
                auto &nextbw = p_bws[nextidx];
                if ((bw->getPars()->end_evolution) || (nextbw->getPars()->end_evolution)){
//                    evalShellThicknessIsolated(idx, Y);
                    (*p_log)(LOG_ERR,AT) << "|error|\n";
                    exit(1);
                }
                double r_i = Y[bw->getPars()->ii_eq + SOL::QS::iR];
                double r_ip1 = Y[nextbw->getPars()->ii_eq + SOL::QS::iR];
                if ((r_i == 0.) || (r_ip1 == 0.)) {
//                (*p_log)(LOG_WARN,AT)<<" shell="<<idx<<" has r_i="<<r_i<<" and r_ip1="<<r_ip1<<"\n";
                    continue;
                }
                double dr_i = r_ip1 - r_i;
                /// evaluate the volume of the shell (fraction of the 4pi)
                double vol_i = (4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i) / bw->getPars()->ncells;
                if (dr_i < 0 or !std::isfinite(dr_i)){
                    (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
                    exit(1);
                }
                /// --------------------------- |
                m_data[Q::idelta][idx] = dr_i;
                m_data[Q::ivol][idx] = vol_i;
            }
            /// for the last shell we have to assume it width
            size_t idx = p_pars->n_active_shells-1;
            double frac = 0.1; // TODO fraction of the shell volume to be used as its width. !!! Inaccurate as we need adjacent shells to get the volume...
            auto & bw = p_bws[idx];
            double r_i =  Y[bw->getPars()->ii_eq + SOL::QS::iR];
            double dr_i = frac * r_i;
            double vol_i = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
            if (dr_i < 0 or !std::isfinite(dr_i)){
                (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
                exit(1);
            }
            m_data[Q::idelta][idx] = dr_i;
            m_data[Q::ivol][idx] = vol_i;
        }
    }
    /// Evaluate the radial extend of a velocity shell. Assume ordered shells. Assumes sorted shells. Assume update kappa
    void updateSortedShellProperties( const double * Y ){
        /// Get thermal energy and temperature of each of the active shells
        double ltot = 0.; double tau_tot=0.;
        for (size_t i = 0; i < p_pars->n_active_shells; i++) {

            size_t idx = m_idxs[i];
            auto & bw = p_bws[i];
            double dr_i = m_data[Q::idelta][i];// TODO add other methods 1/Gamma...
            double vol_i = m_data[Q::ivol][i];
            if (!std::isfinite(vol_i)||vol_i==0.){
                (*p_log)(LOG_ERR,AT)<<"Error.\n";
                exit(1);
            }
            ///store also velocity
            double m_beta = EQS::BetFromMom( Y[bw->getPars()->ii_eq + SOL::QS::imom] );
            double eint2_i = Y[bw->getPars()->ii_eq + SOL::QS::iEint2];
            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
            double ene_th = eint2_i; // TODO is this shock??
            double temp = std::pow(ene_th / A_RAD / vol_i, 0.25); // Stephan-Boltzman law
            m_data[Q::ieth][i] = ene_th;
            m_data[Q::itemp][i] = temp;
            m_data[Q::ibeta][i] = m_beta;

            double m_i = bw->getPars()->M0;
            double m2_i = Y[bw->getPars()->ii_eq + SOL::QS::iM2] * m_i;
//            m_i += m2_i; // TODO should I include this?/
            double rho_i = m_i / vol_i;
            double kappa_i = bw->getPars()->kappa;
            m_data[Q::irho][i] = rho_i;
            m_data[Q::idtau][i] = kappa_i * rho_i * dr_i;
            tau_tot += m_data[Q::idtau][i];

            double r_i = Y[bw->getPars()->ii_eq + SOL::QS::iR];
            double m_tlc = r_i / CGS::c;
            m_data[Q::itdiff][i] = std::max( kappa_i * m_i / m_data[Q::ibeta][i] / r_i / CGS::c, m_tlc); // TODO M0 or M0+M2 )

            m_data[Q::ilum][i] = ene_th / m_data[Q::itdiff][i];
            ltot += m_data[Q::ilum][i];
        }
        p_pars->ltot = ltot;
        p_pars->tautot = tau_tot;

        /// Compute the optical depth from 0 to a given shell
        int idx_photo=-1;
        for (size_t ii = 0 ; ii < p_pars->n_active_shells; ii++){
            /// Optcial depth from Outer Radius Boundary
            double taucum = 0.;
            double tdiff_cum = 0.;
            for (size_t jj = ii ; jj < p_pars->n_active_shells; jj++){
                taucum += m_data[Q::idtau][jj];
                tdiff_cum += m_data[Q::itdiff][jj];
            }
            m_data[Q::itaucum][ii] = taucum;
            m_data[Q::itdiff][ii] = tdiff_cum;
            if ((idx_photo==-1) and (m_data[Q::itaucum][ii] < 1.)){
                idx_photo = (int)ii+1;
            }
            /// optical depth from Inner Radius Boundary
            double tcaum0 = 0.;
            for (size_t jj = 0 ; jj < ii; jj++){
                tcaum0 += m_data[Q::idtau][jj];
            }
            m_data[itaucum0][ii] = tcaum0;

        }

        /// save results for the use in ODE when solving Energy equation
        for (size_t ii = 0 ; ii < p_pars->n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            auto & bw = p_bws[ii];
            if (p_pars->do_thermrad_loss){
                if((cur_idx >= idx_photo && p_pars->thermradloss_at_photo) || (!p_pars->thermradloss_at_photo))
                    bw->getPars()->dElum = m_data[Q::ilum][ii];
            }
//            bw->getPars()->tau_to0 = m_data[Q::itaucum][ii];
//            bw->getPars()->dtau = m_data[Q::idtau][ii];
//            bw->getPars()->is_above_tau1 = idx_photo > cur_idx ? false : true;
        }


#if 0
        /// if there is one shell we cannot have the shell width that comes from shell separation.
        double tau_tot=0.; double ltot = 0.;
        if (p_pars->n_active_shells == 1){
            size_t idx = 0;
            double frac = 0.1; // fraction of the shell volume to be used as its width. Inaccurate as we need adjacent shells to get the volume...
            auto & bw = p_bws[0];
            double kappa_i = bw->getPars()->kappa;
//            r_i =  Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            double dr_i = m_data[Q::idelta][idx];//bw->getPars()->thickness;//frac * r_i; // TODO add other methods 1/Gamma...
            double vol_i = m_data[Q::ivol][idx];//bw->getPars()->volume;// = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
            double m_i = bw->getPars()->M0;
            double m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_i;
//            m_i += m2_i; // TODO should I include this?/
            double rho_i = m_i / vol_i;
            //
            double r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
//            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
//            ene_th = eint2_i; // TODO is this shock??
            //
            m_data[Q::idtau][idx] = kappa_i * rho_i * dr_i;
            tau_tot += m_data[Q::idtau][idx];

            m_data[Q::itaucum][idx] = 0.; // opt. depth to this shell
            double m_tlc = r_i / CGS::c;
            m_data[Q::itdiff][idx] = std::max( kappa_i * m_i / m_data[Q::ibeta][idx] / r_i / CGS::c, m_tlc); // TODO M0 or M0+M2 )

        }
#endif
#if 0
        /// compute shell width from the shell separation
        double tau_tot=0.;
        for (size_t ii=0; ii<p_pars->n_active_shells-1; ii++){
            ///
            size_t idx = m_idxs[ii];
            size_t nextidx = m_idxs[ii+1];
            auto & bw = p_bws[idx];
            auto & nextbw = p_bws[nextidx];

            double dr_i = m_data[Q::idelta][ii];//bw->getPars()->thickness;//r_ip1 - r_i;
            double kappa_i = bw->getPars()->kappa;

            /// evaluate mass within the shell
            double m_i = bw->getPars()->M0;
            double m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_i;
            m_i += m2_i;
            double m_ip1 = nextbw->getPars()->M0;
            double m2_ip1 = Y[nextbw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_ip1;
            m_ip1 += m2_ip1;
            double r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            /// evaluate the volume of the shell (fraction of the 4pi)
            double vol_i = m_data[Q::ivol][ii];//bw->getPars()->volume;//(4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i) / bw->getPars()->ncells;
            /// evaluate density within a shell (neglecting the accreted by the shock!)
            double rho_i = m_i / vol_i;
            m_data[Q::irho][ii] = rho_i;//m_rho[ii] = rho_i;

//            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
//            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different

            /// evaluate optical depth
            double dtau_i = kappa_i * rho_i * dr_i;
            tau_tot += dtau_i;
            if ((!std::isfinite(dtau_i)) || (dtau_i < 0)){
                (*p_log)(LOG_ERR,AT) << "dtau is nan or < 0 and dtau="<<dtau_i<<"\n";
                exit(1);
            }
            m_data[Q::idtau][ii] = dtau_i;//m_dtau[ii] = dtau_i;

            double m_tlc = r_i / CGS::c;
            m_data[Q::itdiff][ii] = std::max(
                    kappa_i * bw->getPars()->M0 / m_data[Q::ibeta][ii] / r_i / CGS::c, m_tlc); // TODO M0 or M0+M2 )
//            tdiff = m_dtau[ii] * Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR] / CGS::c;

        }
#endif
#if 0
        /// Compute the optical depth from 0 to a given shell
        for (size_t ii = 0 ; ii < p_pars->n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            double tau_cum = 0.;
            for (size_t jj = 0 ; jj < p_pars->n_active_shells; jj++){
                size_t other_idx = m_idxs[jj];
                if (other_idx < cur_idx) {
                    tau_cum += m_data[Q::idtau][jj];
                }
            }
            m_data[Q::itaucum][cur_idx] = tau_cum;
        }

        /// Compute cumulative diffusive timescale for all shells
        for (size_t ii = 0 ; ii < p_pars->n_active_shells; ii++){
            double tdiff_cum = 0.;
            for (size_t jj = ii ; jj < p_pars->n_active_shells; jj++){
                tdiff_cum += m_data[Q::itdiff][jj];
            }
//            m_tdiff_out[ii] = std::max( t/diff_cum, m_tlc[ii] ); // Eq. 19 in Ren+19
        }
#endif
#if 0
        /// Compute luminocity from each shell
        double ltot = 0.;
        for (size_t ii=0; ii<p_pars->n_active_shells-1; ii++){
            size_t idx = m_idxs[ii];
            auto & bw = p_bws[idx];
            double eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
            double ene_th = eint2_i; // TODO is this shock??
            m_data[Q::ilum][ii] = ene_th / m_data[Q::itdiff][ii]; // m_dtau[ii] * ene_th * bw->getPars()->M0 * CGS::c * CGS::c / tdiff;
            bw->getPars()->dElum = m_data[Q::ilum][ii];
            ltot += m_data[Q::ilum][ii];
        }
        p_pars->ltot = ltot;
#endif
#if 0
        /// based on where tau = 1 get the photosphere radius, i.e., in which shell is the photosphere
        double tmp_tau = 0; double r_ph;
        size_t idx_photo = p_pars->n_active_shells;
        for (size_t i = p_pars->n_active_shells-1; i > 0; --i){
            size_t idx = m_idxs[i];
            auto & bw = p_bws[idx];
            tmp_tau+=m_data[Q::idtau][i];
            if (tmp_tau > 1.){
                idx_photo = (int)idx+1;
                break;
            }
        }
        if (idx_photo == 0) {
            /// Ejecta totally transparent
            int y = 1;
        }
        else if ((idx_photo == p_pars->n_active_shells-1) && (m_data[Q::idtau][idx_photo-1] > 1)) {
            /// Ejecta totally opaque
            int x = 1;
        }
        else{
            p_pars->rphot = Y[p_bws[idx_photo-1]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            p_pars->tphoto = 0.;
        }
#endif



//        std::cout << "m_rho      ="<< m_data[Q::irho] << "\n";
//        std::cout << "m_lum      ="<< m_data[Q::ilum] << "\n";
//        std::cout << "m_temp     ="<< m_data[Q::itemp] << "\n";
//        std::cout << "m_dtau     ="<< m_data[Q::idtau] << "\n";
//        std::cout << "m_tdiff    ="<< m_data[Q::itdiff] << "\n";
//        std::cout << "m_tau_cum  ="<< m_data[Q::itaucum] << "\n";
//        std::cout << "m_tau_cum0 ="<< m_data[Q::itaucum0] << "\n";



//        (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<mylayer
//                                        <<"] Cumulative optical depth to idx0"
//                                        <<" tau0="<<m_dtau_cum[0]
//                                        <<" tau[-1]="<<m_dtau_cum[n_active_shells-1]<<"\n";


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
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::idR) = dr_i;
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::iRho) = rho_i;
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::iTau) = tau_i; // update blastwave thickness

//        if (idx_photo == m_idxs.back()-1)
//            idx_photo = (int)m_idxs.back();
//        for (auto idx : m_idxs)
//            std::cout << "i_photo=" <<idx_photo<< "; "<< p_bws[idx]->getLastVal(DynRadBlastWave::Q::iTau) << ", ";
//        std::cout <<"\n";

//        double tau_out = 0, tau_out_i = 0;
//        for (auto idx_curr : m_idxs){
//            tau_out = 0.;
//            for (auto idx = m_idxs.back(); idx > idx_curr; idx--){
//                tau_out_i = p_bws[idx]->getLastVal(DynRadBlastWave::Q::iTau);
//                tau_out+=tau_out_i ;


    }
    /// Evaluate shell total mass, volume, density
    /// next three functions compute mass, volume and average density of the entire shell
    double getShellMass(const double * Y_){
        double mtot = 0.;
        for (size_t ii=0; ii<p_pars->n_active_shells; ++ii){
            size_t i_idx = m_idxs[ii];
            double m2 = Y_[p_bws[i_idx]->getPars()->ii_eq + SOL::QS::iM2];
            double m0 = p_bws[i_idx]->getPars()->M0;
            double m2plus0 = (1. + m2) * m0;
            mtot+=m2plus0;
        }
        if (!std::isfinite(mtot) || (mtot < 0)){
            (*p_log)(LOG_ERR,AT) << "mtot is nan or < 0; mtot="<<mtot<<"\n";
            exit(1);
        }
        return mtot;
    }
    double getShellVolume(const double * Y){
        double r0 = Y[p_bws[ m_idxs[0] ]->getPars()->ii_eq + SOL::QS::iR];
        double r1 = Y[p_bws[ m_idxs[p_pars->n_active_shells - 1] ]->getPars()->ii_eq + SOL::QS::iR];
        if ((r0 >= r1)||(r0==0)||(r1==0)){
            (*p_log)(LOG_ERR,AT)<<" r0 > r1. in the shell; r0="<<r0<<" r1="<<r1<<"\n";
            exit(1);
        }
        double delta = r1-r0;
        double volume = (4./3.) * CGS::pi * (r1*r1*r1 - r0*r0*r0) / p_bws[ m_idxs[0] ]->getPars()->ncells;
    }
    double getShellRho(const double * Y){
        double mtot = getShellMass(Y);
        double volume = getShellVolume(Y);
        return mtot / volume;
    }
};

/// Radially/Angular structured Blastwave collection
class Ejecta{
//    VelocityAngularStruct ejectaStructs{};
    std::unique_ptr<EjectaID2> id = nullptr;
    std::unique_ptr<Output> p_out = nullptr;
    std::vector<std::unique_ptr<CumulativeShell>> p_cumShells {};
    std::unique_ptr<logger> p_log = nullptr;
    bool is_ejBW_init = false;
    bool is_ejecta_obsrad_pars_set = false;
    bool is_ejecta_struct_set = false;
    double jet_layer_fnu_stop_frac=1e-5;
    int n_ode_eq{};
    int m_loglevel{};
//    LatStruct::METHOD_eats ejecta_eats_method{};
    Vector & t_arr;
//    size_t m_nshells = 0;
//    size_t m_nlayers = 0;
public:
    bool do_eninj_inside_rhs = false;
    bool run_bws=false, save_dyn=false, do_ele=false, do_spec=false, do_lc=false, do_skymap=false;
    bool do_collision = false;
    bool is_ejecta_obs_pars_set = false;
    bool is_ejecta_anal_synch_computed = false;
    StrDbMap grb_pars; StrStrMap grb_opts;
    std::string working_dir{}; std::string parfilename{};
    Ejecta(Vector & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Ejecta");
        p_out = std::make_unique<Output>(loglevel);
    }
    size_t getNeq() const {
        if (!run_bws)
            return 0;
        else {
            if (p_cumShells.empty()){
                (*p_log)(LOG_ERR,AT)<<" error\n";
                exit(1);
            }
            return (p_cumShells.size() * p_cumShells[0]->nBWs() * SOL::neq);//p_cumShells[0]->getBW(0)->getNeq());
        }
    }
    VecArray & getData(size_t il, size_t ish){ return getShells()[il]->getBW(ish)->getData(); }
    Vector & getTbGrid(){ return t_arr; }
    Vector getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0))
            return t_arr;
        Vector tmp{};
        for (size_t it = 0; it < t_arr.size(); it = it + every_it){
            tmp.push_back(t_arr[it]);
        }
//        Vector tmp2 (tmp.data(), tmp.size());
        return std::move(tmp);
    }
    size_t nlayers() const {
//        return run_bws ? ejectaStructs.structs[0].nlayers : 0;
        return run_bws ? id->nlayers : 0;
    }
    size_t nshells() const {
//        return run_bws ? ejectaStructs.nshells : 0;
        return run_bws ? id->nshells : 0;
    }
    size_t nMaxActiveShells() {
        size_t nsh = 0;
        for (auto & cumShell : getShells() )
            if (nsh < cumShell->getPars()->n_active_shells)
                nsh = cumShell->getPars()->n_active_shells;
        return nsh;
    }
    int ncells() const {
//        return run_bws ? (int)ejectaStructs.structs[0].ncells : 0;
        return run_bws ? (int)id->ncells : 0;
    }
    std::vector<std::unique_ptr<CumulativeShell>> & getShells(){
        if (p_cumShells.empty()){
            (*p_log)(LOG_ERR,AT)<<" ejecta not initialized\n";
            exit(1);
        }
        return p_cumShells;
    }

    void setPars(StrDbMap & pars, StrStrMap & opts,
                 std::string working_dir_, std::string parfilename_,
                 StrDbMap main_pars, size_t ii_eq, size_t iljet){
        working_dir = working_dir_;
        parfilename = parfilename_;
        grb_pars = pars;
        grb_opts = opts;
        /// read GRB afterglow parameters
//        StrDbMap grb_pars; StrStrMap grb_opts;
//        run_bws=false; bool save_dyn=false, do_ele=false, do_spec=false, do_lc=false, do_skymap=false;
        if ((!grb_pars.empty()) || (!grb_opts.empty())) {
            run_bws = getBoolOpt("run_bws", grb_opts, AT, p_log, false, true);
            save_dyn = getBoolOpt("save_dynamics", grb_opts, AT, p_log, false, true);
            do_ele = getBoolOpt("do_ele", grb_opts, AT, p_log, false, true);
            do_spec = getBoolOpt("do_spec", grb_opts, AT, p_log, false, true);
            do_lc = getBoolOpt("do_lc", grb_opts, AT, p_log, false, true);
            do_skymap = getBoolOpt("do_skymap", grb_opts, AT, p_log, false, true);
            for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
                if (main_pars.find(key) == main_pars.end()) {
                    (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                    exit(1);
                }
                grb_pars[key] = main_pars.at(key);
            }
            grb_opts["workingdir"] = working_dir; // For loading Nuclear Heating table
            if (run_bws) {
                std::string fname_ejecta_id = getStrOpt("fname_ejecta_id", grb_opts, AT, p_log, "", true);
                bool use_1d_id = getBoolOpt("use_1d_id", grb_opts, AT, p_log, false, true);
                if (!std::experimental::filesystem::exists(working_dir + fname_ejecta_id)) {
                    (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_ejecta_id << "\n";
                    exit(1);
                }
#if 0
                EjectaID2 ID(working_dir + fname_ejecta_id, m_loglevel);
                if (ID.idtype == EjectaID::IDTYPE::i_id_corr && (!use_1d_id))
                    setEjectaStructNumeric(
                            ID,
                            getBoolOpt("enforce_angular_grid", grb_opts, AT, p_log, false, true),
                            grb_opts);
                else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (use_1d_id))
                    setEjectaStructNumeric(ID, grb_opts);
                else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (!use_1d_id))
                    setEjectaStructNumericUniformInTheta(
                            ID,
                            (size_t) getDoublePar("nlayers", grb_pars, AT, p_log, 30, true),
                            getDoublePar("mfac", grb_pars, AT, p_log, 1.0, true),
                            grb_opts);
                else {
                    (*p_log)(LOG_ERR, AT) << " no valid grb ejecta structure given\n";
                    exit(1);
                }
#endif
                id = std::make_unique<EjectaID2>(
                        working_dir + fname_ejecta_id,
                        getStrOpt("method_eats",grb_opts,AT,p_log,"", true),
                        getBoolOpt("use_1d_id", grb_opts, AT, p_log, false, true),
                        m_loglevel
                        );


//                size_t ii_eq = pba.getMag()->getNeq();
//            size_t nlayers_jet = pba.getGRB()->getBWs().size();
                setEjectaBwPars(grb_pars, grb_opts, ii_eq, iljet);

                /// initialize Ejecta Bound PWN
//                if (run_pwn) {
//                    (*p_log)(LOG_ERR,AT) << " not implemented\n";
//                    exit(1);
////                size_t ii_eq_ = pba.getMag()->getNeq() + pba.getGRB()->getNeq() + pba.getEj()->getNeq();
//                    /// init ejecta-bound PWN
////                pba.getEjPWN()->setPWNpars(pba.getTburst(), pwn_pars, pwn_opts, ii_eq_,
////                                           pba.getEj()->getShells().size());
//                }
//                std::cout<<"nshells="<<nshells()<<" nlayers="<<nlayers()<<"\n";

            }
        }
        else{
            (*p_log)(LOG_INFO, AT) << "ejecta is not initialized and will not be considered.\n";
        }
    }

    void computeAndOutputObservables(StrDbMap main_pars, StrStrMap main_opts){
        /// work on GRB afterglow
        if (run_bws){

            bool lc_freq_to_time = getBoolOpt("lc_use_freq_to_time",main_opts,AT,p_log,false,true);
            Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
            Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

            if (save_dyn)
                saveEjectaBWsDynamics(
                        working_dir,
                        getStrOpt("fname_dyn", grb_opts, AT, p_log, "", true),
                        (int)getDoublePar("save_dyn_every_it", grb_pars, AT, p_log, 1, true),
                        main_pars, grb_pars);
            if (do_ele) {
                setPreComputeEjectaAnalyticElectronsPars();
//                (*p_log)(LOG_INFO, AT) << "jet analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_lc or do_skymap) {
                if (do_lc) {
                    computeSaveEjectaLightCurveAnalytic(
                            working_dir,
                            getStrOpt("fname_light_curve", grb_opts, AT, p_log, "", true),
                            getStrOpt("fname_light_curve_layers", grb_opts, AT, p_log, "", true),
                            lc_times, lc_freqs, main_pars, grb_pars, lc_freq_to_time);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
                }
                if (do_skymap) {
                    computeSaveEjectaSkyImagesAnalytic(
                            working_dir,
                            getStrOpt("fname_sky_map", grb_opts, AT, p_log, "", true),
                            skymap_times, skymap_freqs, main_pars, grb_pars);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";
                }
            }
        }
    }

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
            if (p_cumShells[il]->getPars()->n_active_shells > n_active_max){
                n_active_max = p_cumShells[il]->getPars()->n_active_shells;
                il_with_min_nact = (int)il;
            }
            if (p_cumShells[il]->getPars()->n_active_shells < n_active_min){
                n_active_min = p_cumShells[il]->getPars()->n_active_shells;
                il_with_max_nact = (int)il;
            }
            /// find number of shells in each layer that (i) accelerating (ii) decelerating
            size_t n_accel = 0;
            size_t n_decel = 0;
            for (size_t ish = 0; ish < p_cumShells[il]->getPars()->n_active_shells; ish++){
                auto & bws = p_cumShells[il]->getBW(ish);
                double MomIm1 = Ym1[bws->getPars()->ii_eq + SOL::QS::imom];
                double MomI = Y[bws->getPars()->ii_eq + SOL::QS::imom];
                double Eint2I = Y[bws->getPars()->ii_eq + SOL::QS::iEint2];
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
//                if (p_cumShells[il]->getBW(ish)->getVal(RadBlastWave::Q::imom,(int)it) > mom) {
//                    mom = p_cumShells[il]->getBW(ish)->getVal(RadBlastWave::Q::imom,(int)it) > mom;
//                    fastest_l = il;
//                    fastest_sh = ish;
//                }
//            }
//        }
//        sstream << "[Ej: "<<"[l="<<fastest_l<<", sh="<<fastest_sh<<"]"
//                << " Mom=" << string_format("%.2e",p_cumShells[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::imom,(int)it))
//                << " R=" << string_format("%.2e",p_cumShells[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::iR,(int)it))
//                << " Eint=" << string_format("%.2e",p_cumShells[fastest_l]->getBW(fastest_sh)->getVal(RadBlastWave::Q::iEint2,(int)it))
//                << "] ";
    }

private:
#if 0
    void setEjectaStructNumericUniformInTheta(EjectaID & id, size_t nlayers, double mfac, StrStrMap & opts){
        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        if (!run_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initUniform(id, nlayers, mfac, ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaStructNumeric(EjectaID & id, StrStrMap & opts){
        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        if (!run_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initCustom(id, ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaStructNumeric(EjectaID & id, bool force_grid, StrStrMap & opts){
        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        if (!run_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initCustom(id, force_grid, ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
#endif
    void setEjectaBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq, size_t n_layers_jet){

        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        if (!run_bws)
            return;
        bool is_within = false;
        std::vector<size_t> which_within{};
        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t nshells_ = nshells();//ejectaStructs.nshells;
        size_t n_layers_ej_ = nlayers();//ejectaStructs.structs[0].nlayers;
        if (n_layers_ej_ == 0){
            (*p_log)(LOG_ERR,AT)<<" no layers found to evolve!\n";
            exit(1);
        }
        std::vector<std::vector<size_t>> n_empty_images_layer_shell;
        for (auto & n_empty_images_layer : n_empty_images_layer_shell)
            n_empty_images_layer.resize(nshells_);
        /// include blastwave collision between velocioty shells into the run
        for(size_t il = 0; il < n_layers_ej_; il++){
            p_cumShells.push_back(
                    std::make_unique<CumulativeShell>(t_arr, nshells_, il,
                                                      p_log->getLogLevel()) );
            p_cumShells[il]->setPars(pars, opts);
            for (size_t ish = 0; ish < nshells_; ish++){
                auto & bw = p_cumShells[il]->getBW(ish);
//                auto & struc = ejectaStructs.structs[ish];
                bw->setAllParametersForOneLayer(id, pars, opts, il, ii_eq);
/// Override the layer-to-use
                if (bw->getPars()->which_jet_layer_to_use == 0){
                    bw->getPars()->which_jet_layer_to_use = 0; // the fastest
                }
                else if(n_layers_ej_ == 0){
                    n_ejecta_empty_images += 1;
                    n_empty_images_layer_shell[ish].emplace_back(il);
//                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
                }
                else if(n_layers_ej_ == 0){
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
                ii_eq += SOL::neq;//bw->getNeq();

//                bw->getRad()->setEatsPars(pars, opts);
//                bw->getSynchAnPtr()->setPars( pars, opts );

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

        for(size_t i = 0; i < ejectaStructs.nshells_; i++){
            std::vector<size_t> n_empty_images_layer;
            auto & struc = ejectaStructs.structs[i];
            if (n_layers_i != struc.nlayers){
                (*p_log)(LOG_ERR,AT)<<" expected nlayers="<<n_layers_i<<" (from 0th shell) got="<<struc.nlayers<<"\n";
                exit(1);
            }
//            size_t n_layers_i = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for(size_t j = 0; j < n_layers_i; j++) { // ii = il + nlayers * ish
                p_bws.emplace_back(std::make_unique<DynRadBlastWave>(t_grid, i, j, p_pars->loglevel));
                setAllParametersForOneLayer(struc, *(p_bws[ii]), pars, opts, j, ii_eq);
                /// Override the layer-to-use
                if (p_bws[ii]->getPars()->which_jet_layer_to_use == 0){
                    p_bws[ii]->getPars()->which_jet_layer_to_use = 0; // the fastest
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
                else if ((p_bws[ii]->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
                    p_bws[ii]->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
                }
                else if ((p_bws[ii]->getPars()->which_jet_layer_to_use < n_layers_jet) &&
                         (p_bws[ii]->getPars()->which_jet_layer_to_use > -1)){
                    //
                }
                else{
                    std::cerr << " which_jet_layer_to_use="<<p_bws[ii]->getPars()->which_jet_layer_to_use
                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
                              <<" for the slowest"
                              <<" or any N in between the two for a specific jet layer \n"
                              << "Exiting..."
                              << "\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }
                ii_eq += p_bws[ii]->getNeq();

                p_bws[ii]->setEatsPars(pars, opts);
                p_bws[ii]->getSynchAnPtr()->setPars( pars, opts );

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
//                    auto &ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
                    size_t n_layers_i = nlayers();//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
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
        do_eninj_inside_rhs = getBoolOpt("do_eninj_inside_rhs", opts, AT, p_log, "no", false);
//        std::cout<<"nshells_="<<nshells()<<" nlayers="<<nlayers()<<"\n";
//        std::cout<<"nshells_="<<nshells_<<" nlayers="<<n_layers_ej_<<"\n";

        (*p_log)(LOG_INFO,AT) << "finished initializing ejecta. "
                                 "nshells="<<nshells()<<" nlayers="<<nlayers()<<"\n";
//        for(size_t il = 0; il < n_layers_ej_; il++){
//            for (size_t ish = 0; ish < nshells_; ish++) {
//                p_cumShells[il]->getBW(ish)->getRad();
//            }
//        }
    }

    double ej_rtol = 1e-5;

    /// OUTPUT
    void saveEjectaBWsDynamics(std::string workingdir, std::string fname, size_t every_it,
                               StrDbMap & main_pars, StrDbMap & ej_pars){
        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

//        size_t nshells = p_cumShells->nshells();
//        size_t nlayers = p_cumShells->nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        auto & models = getShells();

        std::vector<std::string> table_names;
        std::vector<std::vector<std::vector<double>>> tot_dyn_out ( nshells() * nlayers() );
        size_t i = 0;
        VecVector other_data;
        std::vector<std::string> other_names;
        auto & arr_names = BW::m_vnames;//models[0]->getBW(0)->m_vnames;//models[0][0].getBWdynPtr()->m_vnames;
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            for(size_t ilayer = 0; ilayer < nlayers(); ilayer++){
                table_names.push_back("shell="+std::to_string(ishell)+" layer="+std::to_string(ilayer));
                tot_dyn_out[i].resize( arr_names.size() );
                auto & bw = models[ilayer]->getBW(ishell);

                std::unordered_map<std::string,double> group_attr{
                        {"Gamma0",bw->getPars()->Gamma0},
                        {"M0",bw->getPars()->M0},
                        {"R0",bw->getPars()->R0},
                        {"theta0",bw->getPars()->theta_b0},
                        {"theta_max",bw->getPars()->theta_max},
                        {"tb0",bw->getPars()->tb0},
                        {"ijl",bw->getPars()->ijl},
                        {"ncells",bw->getPars()->ncells},
                        {"ilayer",bw->getPars()->ilayer},
                        {"ishell",bw->getPars()->ishell},
                        {"ctheta0",bw->getPars()->ctheta0},
                        {"E0",bw->getPars()->E0},
                        {"theta_c_l",bw->getPars()->theta_c_l},
                        {"theta_c_h",bw->getPars()->theta_c_h},
                        {"eps_rad",bw->getPars()->eps_rad},
                        {"entry_time",bw->getPars()->entry_time},
                        {"entry_r",bw->getPars()->entry_r},
                        {"first_entry_r",bw->getPars()->first_entry_r},
                        {"min_beta_terminate",bw->getPars()->min_beta_terminate}
                };

                group_attrs.emplace_back( group_attr );

                for (size_t ivar = 0; ivar < arr_names.size(); ivar++) {
                    for (size_t it = 0; it < bw->getTbGrid().size(); it = it + every_it)
                        tot_dyn_out[i][ivar].emplace_back( (*bw)[ static_cast<BW::Q>(ivar) ][it] );
                }
                i++;
            }
        }
        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells() },
                {"nlayers", nlayers() }
        };
//        attrs.insert(ej_pars.begin(),ej_pars.end());
//        attrs.insert(main_pars.begin(),main_pars.end());

        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VecVectorOfVectorsAsGroupsH5(tot_dyn_out, table_names, arr_names,
                                            workingdir+fname, attrs, group_attrs);
    }
private:
    void computeEjectaSkyMapPieceWise( std::vector<Image> & images, double obs_time, double obs_freq ){

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();

        for(size_t il = 0; il < nlayers_; il++){
            for (size_t ish = 0; ish < nshells_; ish++) {
                p_cumShells[il]->getBW(ish)->getRad();
            }
        }

//        std::vector<Image> images;
//        images.resize(ejectaStructs.nshells);
        if (images.empty()){
            (*p_log)(LOG_ERR,AT) << " empty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nshells_){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }
        // ejectaStructs.nshells);
//        size_t ii = 0;
        size_t n_jet_empty_images = 0;


//        std::vector<size_t> n_empty_images_layer;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        std::vector<Image> tmp (nlayers_);
        for (auto & _tmp : tmp)
            _tmp.resize( ncells_ );
        Image tmp_pj( ncells_); // std::vector<Image> tmp_pj (ejectaStructs.structs[0].nlayers);
        Image tmp_cj( ncells_); // std::vector<Image> tmp_cj (ejectaStructs.structs[0].nlayers);
        for (size_t ishell = 0; ishell < nshells_; ishell++){
            for (auto & _tmp : tmp)
                _tmp.clearData();
            tmp_pj.clearData(); tmp_cj.clearData();
//            auto & shells = getShells();
//            if (shells.size() < nlayers_-1){
//                std::cerr << AT << "\n"; exit(1);
//            }
            std::vector<size_t> n_empty_images_layer;
            for (size_t ilayer = 0; ilayer < nlayers_; ilayer++){
//                auto & model = shells[ilayer];
//                tmp.emplace_back( model->evalImagePW(obs_time, obs_freq) );
//                auto & x1 = model->getBW(ishell);
//                auto & x2 = x1->getRad();
                auto & bw_rad = p_cumShells[ilayer]->getBW(ishell)->getRad();
                bw_rad->evalImagePW(tmp[ilayer], tmp_pj, tmp_cj, obs_time, obs_freq);

//                std::cout<<"ish="<<ishell<<" il="<<ilayer<<" pj="<<tmp_pj.m_f_tot<<" cj="<<tmp_cj.m_f_tot<<"\n";
                if (tmp[ilayer].m_f_tot == 0){
                    n_jet_empty_images += 1;
                    n_empty_images_layer.emplace_back(ilayer);
                }
//                ii ++ ;
            }
            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
            int x = 1;

//            images.emplace_back( combineImages(ejectaStruct,tmp) );
            combineImages(images[ishell], ncells_, nlayers_, tmp) ;
        }

        /// print which layers/shells gave empty image
        if (p_log->getLogLevel() == LOG_INFO) {
            if (n_jet_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an empty images for total n="
                      << n_jet_empty_images << " layers. Specifically:\n";
                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
//                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                size_t n_layers_ej = nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [ ";
                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
                        ccerr << n_empty_images[ish][il] << " ";
                    }
                    ccerr << "] / (" << nlayers_ << " total layers) \n";
                }
            }
        }

//        return std::move( images );
    }

    std::vector<VecVector> evalEjectaLightCurves( Vector & obs_times, Vector & obs_freqs){
        (*p_log)(LOG_INFO,AT)<<" starting ejecta light curve calculation\n";
//        size_t nshells = p_cumShells->nshells();
//        size_t nlayers = p_cumShells->nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();
        std::vector<VecVector> light_curves(nshells()); // [ishell][i_layer][i_time]
        for (auto & arr : light_curves){
            size_t n_layers_ej = nlayers();//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
            arr.resize(n_layers_ej);
            for (auto & arrr : arr){
                arrr.resize( obs_times.size(), 0. );
            }
        }
        double flux_pj, flux_cj; size_t ii = 0;
//        Image image;
        double rtol = ej_rtol;
        Image image_i ( ncells() );
        Image im_pj ( ncells() );
        Image im_cj ( ncells() );
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            image_i.clearData();
            im_pj.clearData();
            im_cj.clearData();
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStructs.structs[ishell].nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer];//ejectaModels[ishell][ilayer];
//                model->setEatsPars( pars, opts );
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC ntimes="<<obs_times.size()
                        << " vel_shell="<<ishell<<"/"<<nshells()-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers()
                        << " phi_cells="<<EjectaID2::CellsInLayer(ilayer)<<"\n";
                model->getBW(ishell)->getRad()->evalForwardShockLightCurve(
                        id->method_eats,
                        image_i, im_pj, im_cj, light_curves[ishell][ilayer], obs_times, obs_freqs);
                ii ++;
            }
        }
        return std::move( light_curves );
    }

public:
    /// electrons
    void setPreComputeEjectaAnalyticElectronsPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Ejecta analytic electron pars...\n";

        if (!run_bws){
            std::cerr << " ejecta BWs were not evolved. Cannot compute electrons (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = getShells();
        for (auto & model : models) {
            for (auto & bw : model->getBWs()) {
                bw->getRad()->computeForwardShockElectronAnalyticVars();
                bw->getRad()->computeForwardShockSynchrotronAnalyticSpectrum();
            }
        }
        is_ejecta_anal_synch_computed = true;
    }

    void updateEjectaObsPars(StrDbMap pars) {

        auto & models = getShells();

//        size_t nshells = nshells();
//        size_t nlayers = p_cumShells->nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        (*p_log)(LOG_ERR,AT) << "Updating Ejecta observer pars...\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++) {
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer]->getBW(ishell);//ejectaModels[ishell][ilayer];
                model->getRad()->updateObsPars(pars);
                ii++;
            }
        }
//        p_pars->is_ejecta_obs_pars_set = true;
    }

    void computeSaveEjectaSkyImagesAnalytic(std::string workingdir, std::string fname, Vector times, Vector freqs,
                                            StrDbMap & main_pars, StrDbMap & ej_pars){
        if (!run_bws)
            return;

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta sky image with analytic synchrotron...\n";

        if (!is_ejecta_anal_synch_computed){
            std::cerr  << "ejecta analytic electrons were not evolved. Cannot compute images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            std::cerr<< "ejecta observer parameters are not set. Cannot compute image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();

//        auto & ejectaStruct = ejectaStructs;
//        size_t nshells = ejectaStruct.nshells;
        Image dummy(1,0);

        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};

        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(dummy.m_names.size());
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); ++i_vn) {
                    out[ii][i_vn].resize(nshells_);
                }
                ii++;
            }
        }

        VecVector other_data{
                times,
                freqs
        };
        ii = 0;
        std::vector<std::string> other_names { "times", "freqs" };
        std::vector<Image> images(nshells_);

        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            VecVector total_flux_shell( nshells_ );
            for (auto & total_flux_shel : total_flux_shell)
                total_flux_shel.resize( times.size(), 0.0 );
            for (size_t it = 0; it < times.size(); ++it){
//                auto images = computeEjectaSkyMapPieceWise( times[it],freqs[ifreq]);
//                std::vector<Image> images(nshells_);
                computeEjectaSkyMapPieceWise( images, times[it],freqs[ifreq]);
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells_; ish++) {
                        out[ii][i_vn][ish] = images[ish].m_data[i_vn];//arrToVec(images[ish].m_data[i_vn]);
                    }
                }
                for (size_t ish = 0; ish < nshells_; ish++) {
                    tota_flux[it] += images[ish].m_f_tot;
                    total_flux_shell[ish][it] = images[ish].m_f_tot;
                }
                ii++;
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );

            for (size_t ish = 0; ish < nshells_; ish++){
                other_data.emplace_back( total_flux_shell[ish] );
                other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) +
                                          " shell=" + string_format("%d", ish));
            }

        }

        std::vector<std::string> group_names{};
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            for (size_t it = 0; it < times.size(); it++){
                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
            }
        }


//        for (size_t it = 0; it < times.size(); it++){
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
//                                       + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
//            }
//        }

        auto in_group_names = dummy.m_names;

        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells_},
                {"nshells", nlayers_}
        };
//                {"thetaObs", p_cumShells->getShells()[0]->getBW(0)->getEatsPars()->theta_obs },
//                {"d_L", p_cumShells->getShells()[0]->getBW(0)->getEatsPars()->d_l },
//                {"z",  p_cumShells->getShells()[0]->getBW(0)->getEatsPars()->z },
//                {"eps_e",  p_cumShells->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
//                {"eps_b",  p_cumShells->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
//                {"eps_t",  p_cumShells->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
//                {"p",  p_cumShells->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
//                {"ksi_n",  p_cumShells->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
//        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(workingdir+fname,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    void computeSaveEjectaLightCurveAnalytic(std::string workingdir,std::string fname, std::string fname_shells_layers,
                                             Vector lc_times, Vector lc_freqs, StrDbMap & main_pars, StrDbMap & ej_pars,
                                             bool lc_freq_to_time){

        Vector _times, _freqs;
        cast_times_freqs(lc_times,lc_freqs,_times,_freqs,lc_freq_to_time,p_log);

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta light curve with analytic synchrotron...\n";

//        size_t nshells = p_cumShells->nshells();
//        size_t nlayers = p_cumShells->nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        if (!is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            std::cerr << " ejecta observer parameters are not set. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

//        auto & tmp = getShells()[0]->getBW(0)->getSynchAnPtr();

        std::vector< // layers / shells
                std::vector< // options
                        std::vector<double>>> // freqs*times
        out {};

        /// evaluate light curve
        auto light_curve = evalEjectaLightCurves( _times, _freqs);

        /// save total lightcurve
        size_t n = _times.size();
        Vector total_fluxes (n, 0.0);
        for (size_t itnu = 0; itnu < n; ++itnu) {
            size_t ishil = 0;
            for (size_t ishell = 0; ishell < nshells(); ++ishell) {
                for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                    total_fluxes[itnu] += light_curve[ishell][ilayer][itnu];
                    ishil++;
                }
            }
        }
        std::vector<std::string> other_names { "times", "freqs", "total_fluxes" };
        VecVector out_data {_times, _freqs, total_fluxes};

        std::unordered_map<std::string,double> attrs{ {"nshells", nshells()}, {"nlayers", nlayers()} };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);


        /// save light curve for each shell and layer
        if (fname_shells_layers == "none")
            return;
        std::vector<std::string> group_names;
        VecVector total_fluxes_shell_layer(nshells()*nlayers());
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ++ishell) {
            for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                total_fluxes_shell_layer[ii].resize(n,0.);
                for (size_t ifnu = 0; ifnu < n; ifnu++){
                    total_fluxes_shell_layer[ii][ifnu] = light_curve[ishell][ilayer][ifnu];
                }
                ii++;
            }
        }
        total_fluxes_shell_layer.emplace_back(_times);
        total_fluxes_shell_layer.emplace_back(_freqs);
        total_fluxes_shell_layer.emplace_back(total_fluxes);

        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_fluxes");
        p_out->VectorOfVectorsH5(total_fluxes_shell_layer, group_names, workingdir+fname,  attrs);
    }

};




#if 0
/*
 * Class to hold blastwaves with the same angular coordinate (shells of the same layer)
 * Here the velocity distribution, relative position and photosphere can be found.
 */
class CumulativeShell_{

//    enum Q { ikapp};
    std::unique_ptr<logger> p_log = nullptr;
    std::unique_ptr<BlastWaveCollision> p_coll = nullptr;
    std::vector<std::unique_ptr<RadBlastWave>> p_bws{};
    std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws_ej{};
    std::vector<std::vector<size_t>> relative_position{};
    std::vector<size_t> m_idxs{};
    std::vector<size_t> m_idxs0{};
    std::vector<size_t> idx_tau_eq1{};
    Vector m_rho{};
    Vector m_beta{};
    Vector m_temp{};
    Vector m_radii_init{};
    Vector m_delta{};
    Vector m_lum{};
    Vector m_tdiff{};
    Vector m_tdiff_out{};
    Vector m_tlc{};
//    Vector m_radii_sorted;
    Vector m_dtau{};
    Vector m_frac{};
    Vector m_dtau_cum{};
    VecVector m_shell_data{};
    size_t m_tarr_size{};
    double m_tot_lum = 0.;
    double m_photo_r = 0.;
    double m_photo_teff = 0.;
    const Vector m_t_grid{};
    bool do_collision{};
    double total_mass{};
    double mass_averaged_beta{};
public:
    size_t m_nshells;
    size_t n_active_shells;
    size_t mylayer{};
    CumulativeShell_( Vector t_grid, size_t nshells, int ilayer, int loglevel ) : m_t_grid(t_grid) {
        m_nshells = nshells;
        n_active_shells = nshells; // as shells collide this number will decrease

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CumulativeShell");
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws.emplace_back( std::make_unique<DynRadBlastWave>(t_grid, ishell, ilayer, loglevel ) );

        idx_tau_eq1.resize(t_grid.size(), 0);
        relative_position.resize( nshells );
        for (auto & arr : relative_position)
            arr.resize( t_grid.size() );

        m_rho.resize(nshells, 0.0);
        m_frac.resize(nshells, 0.0);
        m_temp.resize(nshells, 0.0);
        m_delta.resize(nshells, 0.0);
        m_radii_init.resize(nshells, std::numeric_limits<double>::max());
//        m_radii_sorted.resize(nshells, std::numeric_limits<double>::max());
        m_dtau.resize(nshells, 0.0);
        m_dtau_cum.resize(nshells, 0.0);
        m_beta.resize(nshells, 0.0);
        m_lum.resize(nshells, 0.0);
        m_tdiff.resize(nshells, 0.0);
        m_tdiff_out.resize(nshells, 0.0);
        m_tlc.resize(nshells, 0.0);
//        m_shell_data.resize(m_nshells);
//        for (auto arr:)
        m_idxs.resize(nshells, nshells-1); // fill with last value (if not active)
        /// initial shell indexes
        m_idxs0.resize(nshells, nshells-1); // fill with last value (if not active)
        for (size_t i = 0; i<nshells; i++)
            m_idxs0[i] = i;
        mylayer = ilayer;
        m_tarr_size = t_grid.size();
//        p_colsolve = std::make_unique<FsolvePars>(p_bws[0]->getEos());
//        p_colsolve = new FsolvePars(p_bws[0]->getEos());
//        p_colsolve->p_eos;
        p_coll = std::make_unique<BlastWaveCollision>(loglevel);


    }

    Vector & getSortedRadii(){return m_radii_init;}

    inline double getR(size_t i){return m_radii_init[m_idxs[i]];}
    inline Vector & getRvec(){return m_radii_init; }
    inline std::vector<size_t> & getIdx(){return m_idxs; }
    inline Vector & getBetaVec(){return m_beta; }
    inline Vector & getRhoVec(){return m_rho; }
    inline Vector & getTauVec(){return m_dtau; }
    inline Vector & getTempVec(){return m_temp; }
    inline Vector & getDeltaVec(){return m_delta; }
    std::unique_ptr<RadBlastWave> & getBW(size_t ish){
        if (p_bws.empty()){
            (*p_log)(LOG_ERR, AT) << " shell does not contain blast waves\n";
            exit(1);
        }
        if (ish > p_bws.size()){
            (*p_log)(LOG_ERR, AT) << "invalid memory accessed\n"; exit(1);
        }
        return p_bws[ish];
    }
    inline size_t nBWs() const { return p_bws.size();}
    inline std::vector<std::unique_ptr<RadBlastWave>> & getBWs() { return p_bws; }
    /// update the number of active shells
    void updateActiveShells(){
        std::fill(m_idxs.begin(), m_idxs.end(),std::numeric_limits<size_t>::max());
        n_active_shells = m_nshells;
        size_t idx = 0;
        for (size_t i=0; i<n_active_shells; i++) {
            /// loop over all blastwaves and collect the active ones
            if (p_bws[i]->getPars()->end_evolution) {
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
            m_radii_init[i] = Y[ p_bws[idx]->getPars()->ii_eq + DynRadBlastWave::QS::iR ];
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
        auto & bw1 = p_bws[idx1];
        auto & bw2 = p_bws[idx2];
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
            if (p_bws[i]->getPars()->end_evolution) {
                continue;
            }
            /// if shell is active record its current index (order)
            m_idxs[idx] = i; // add only active shells to the list
//            std::cerr << "i="<<i<<" R="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::iR ]<<"\n";
            m_radii_init[idx] = Y[p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::iR ];
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
                    double r0im1 = Ym1_[p_bws[i_idx]->getPars()->ii_eq + DynRadBlastWave::QS::iR];//p_bws[i]->getVal(BW::iR, (int) (ix - 1));
                    double r0i   = Y_[p_bws[i_idx]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
                    double r1im1 = Ym1_[p_bws[j_idx]->getPars()->ii_eq + DynRadBlastWave::QS::iR];//p_bws[j]->getVal(BW::iR, (int) (ix - 1));
                    double r1i   = Y_[p_bws[j_idx]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
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
//                            << " mom=" << Y[p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom]
//                            << " (mom0=" << p_bws[i]->getPars()->mom0 << ") "
//                            << " has overrun shell j=" << overrun_shells[0]
//                            << " with momenta "
//                            << Y[p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom]
//                            << " (mom0=" << p_bws[overrun_shells[0]]->getPars()->mom0 << ") "
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
        auto & bw1 = p_bws[idx1];
        auto & bw2 = p_bws[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                    "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
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
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
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

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
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
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::QS::imom] = _mom;
        Y[iieq + DynRadBlastWave::QS::iEint2] = _eint2; // use ODE units

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
//        Y[iieq + DynRadBlastWave::QS::i] = i_gM * EQS::Beta(i_gM);


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
        auto & bw1 = p_bws[idx1];
        auto & bw2 = p_bws[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                    "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
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
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
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

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            m2_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
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
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::QS::imom] = _mom;
        Y[iieq + DynRadBlastWave::QS::iEint2] = _eint2; // use ODE units

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
//        Y[iieq + DynRadBlastWave::QS::i] = i_gM * EQS::Beta(i_gM);


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
        auto & bw1 = p_bws[idx1];
        auto & bw2 = p_bws[idx2];
        if (bw1->getPars()->end_evolution || bw2->getPars()->end_evolution){
            (*p_log)(LOG_ERR,AT) << " of the shells staged for collision is not active: "
                                   "idx1=" << idx1 << " idx2=" << idx2 << "\n";
            exit(1);
        }
        /// Extract data for first shell
//        double m_gam1, m_beta1, m_mass1, m_eint1, m_adi1;
        p_colsolve->m_gam1 = EQS::GamFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta1 = EQS::BetFromMom( Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint1 = Y[bw1->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
        p_colsolve->m_adi1 = bw1->getEos()->getGammaAdi(p_colsolve->m_gam1, p_colsolve->m_beta1);
        /// apply units for the first shell (mass and energy)
        p_colsolve->m_mass1 = bw1->getPars()->M0 + (p_colsolve->m_mass1 * bw1->getPars()->M0);
        p_colsolve->m_eint1 = p_colsolve->m_eint1 * (bw1->getPars()->M0 * CGS::c * CGS::c);
        /// extract data for the second shell
//        double m_gam2, m_beta2, m_mass2, m_eint2, m_adi2;
        p_colsolve->m_gam2 = EQS::GamFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_beta2 = EQS::BetFromMom( Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
        p_colsolve->m_mass2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
        p_colsolve->m_eint2 = Y[bw2->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
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
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            eint_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw1->getPars()->M0 / bw2->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw1->getPars()->M0 / bw2->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            //
            //
            bw1->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM / bw2->getPars()->M0;//Y[iieq + DynRadBlastWave::QS::iM2] * bw2->getPars()->M0 / bw1->getPars()->M0;
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

        }
        else{
            iieq = bw2->getPars()->ii_eq;
            iieq_other = bw1->getPars()->ii_eq;
            eint_before = Y[iieq + DynRadBlastWave::QS::iEint2];
            eint_before = Y[iieq + DynRadBlastWave::QS::iM2];
            Y[iieq + DynRadBlastWave::QS::iEad2] += Y[iieq + DynRadBlastWave::QS::iEint2] * bw2->getPars()->M0 / bw1->getPars()->M0;
            Y[iieq + DynRadBlastWave::QS::iEsh2] += Y[iieq + DynRadBlastWave::QS::iEsh2] * bw2->getPars()->M0 / bw1->getPars()->M0;
//            Y[iieq + DynRadBlastWave::QS::iR] = rcoll;
            ///
            ///
            bw2->getPars()->M0 = bw2->getPars()->M0+bw1->getPars()->M0;//i_mM; // update the total mass of the shell
            Y[iieq + DynRadBlastWave::QS::iM2] = i_mM/ bw2->getPars()->M0;// Y[iieq + DynRadBlastWave::QS::iM2] * bw1->getPars()->M0 / bw2->getPars()->M0;
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
        }
        /// using the solution (mass, lorentz factor, energy) update the state vector
        double _mom = i_gM * EQS::Beta(i_gM);
        double _eint2 = i_eM / ( i_mM * CGS::c * CGS::c );
        Y[iieq + DynRadBlastWave::QS::imom] = _mom;
        Y[iieq + DynRadBlastWave::QS::iEint2] = _eint2; // use ODE units

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
//        Y[iieq + DynRadBlastWave::QS::i] = i_gM * EQS::Beta(i_gM);


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
                    double r0im1 = p_bws[i]->getVal(BW::iR,(int)(ix-1));
                    double r0i = Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::iR ];
                    double r1im1 = p_bws[j]->getVal(BW::iR,(int)(ix-1));
                    double r1i = Y[ p_bws[j]->getPars()->ii_eq + DynRadBlastWave::QS::iR ];
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
                        <<" mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                        <<" (mom0="<<p_bws[i]->getPars()->mom0<<") "
                        <<" has overrun n_shells="<<overrun_shells.size()
                        <<" with momenta ["<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                        <<" (mom0="<<p_bws[overrun_shells[0]]->getPars()->mom0<<") "
                        <<" ... " << Y[ p_bws[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                        <<"] \n";

//                (*p_log)(LOG_ERR,AT) << " shell i="<<i<<" has overrun n="<<overrun_shells.size()<<" other shells.\n";
                return -1;
            }
            if (overrun_shells.size() == 1){
                (*p_log)(LOG_WARN,AT)
                        <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                        <<" mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                        <<" (mom0="<<p_bws[i]->getPars()->mom0<<") "
                        <<" has overrun shell j="<<overrun_shells[0]
                        <<" with momenta "<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                        <<" (mom0="<<p_bws[overrun_shells[0]]->getPars()->mom0<<") "
                        <<"\n";
                n_overruns += 1;
            }
            /// get time till the nearest collision
            min_t_at_which_overrun.push_back()
        }
        if (n_overruns > 1){
            (*p_log)(LOG_ERR,AT)

                    << " at t="<<x << " multimple overruns have occured n_overruns="<<n_overruns<<"\n";
//                    <<" with mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
//                    <<" has overrun n_shells="<<overrun_shells.size()
//                    <<" with momenta ["<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::iR ]
//                    <<" ... " << Y[ p_bws[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::QS::iR ]
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
            if (p_bws[i]->getPars()->end_evolution)
                continue;
            m_radii_init[i] = Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::iR ];
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
                    <<" mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                    <<" (mom0="<<p_bws[i]->getPars()->mom0<<") "
                    <<" has overrun n_shells="<<overrun_shells.size()
                    <<" with momenta ["<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                    <<" (mom0="<<p_bws[overrun_shells[0]]->getPars()->mom0<<") "
                    <<" ... " << Y[ p_bws[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                    <<"] \n";

//                (*p_log)(LOG_ERR,AT) << " shell i="<<i<<" has overrun n="<<overrun_shells.size()<<" other shells.\n";
                return -1;
            }
            if (overrun_shells.size() == 1){
                (*p_log)(LOG_WARN,AT)
                <<" at t="<<x << " ("<<ix<<"/"<<m_tarr_size << ")" <<" [ilayer="<<mylayer<<" ishell="<<i<<"] "
                <<" mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                <<" (mom0="<<p_bws[i]->getPars()->mom0<<") "
                <<" has overrun shell j="<<overrun_shells[0]
                <<" with momenta "<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
                <<" (mom0="<<p_bws[overrun_shells[0]]->getPars()->mom0<<") "
                <<"\n";
                n_overruns += 1;
            }
        }
        if (n_overruns > 1){
            (*p_log)(LOG_ERR,AT)

                    << " at t="<<x << " multimple overruns have occured n_overruns="<<n_overruns<<"\n";
//                    <<" with mom="<<Y[ p_bws[i]->getPars()->ii_eq + DynRadBlastWave::QS::imom ]
//                    <<" has overrun n_shells="<<overrun_shells.size()
//                    <<" with momenta ["<<Y[ p_bws[overrun_shells[0]]->getPars()->ii_eq + DynRadBlastWave::QS::iR ]
//                    <<" ... " << Y[ p_bws[overrun_shells.back()]->getPars()->ii_eq + DynRadBlastWave::QS::iR ]
//                    <<" \n";

//            (*p_log)(LOG_ERR,AT) << "n_overruns="<<n_overruns<<"\n";
            return -1;
        }
        return 1;
    }

#endif

    void evalShellThicknessIsolated(size_t idx, const double * Y){

        double r_i=0., r_ip1=0., dr_i = 0., vol_i=0.;
        double frac = 0.1; // TODO fraction of the shell volume to be used as its width. !!! Inaccurate as we need adjacent shells to get the volume...
        auto & bw = p_bws[idx];
        r_i =  Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
        dr_i = frac * r_i;
        vol_i = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
        bw->getPars()->thickness = dr_i;
        bw->getPars()->volume = vol_i;
    }
    /// Evaluate shell thickness using various methods. Assumes sorted shells (no overruns)
    void evalShellThickness( const double * Y ){
        double r_i=0., r_ip1=0., dr_i = 0., vol_i=0.;
        /// if there is one shell we cannot have the shell width that comes from shell separation.
        if (n_active_shells == 1){
            evalShellThicknessIsolated(0, Y);
        }
        else{
            for (size_t ii=0; ii<n_active_shells-1; ii++) {
                r_i = 0., r_ip1 = 0., dr_i = 0., vol_i = 0.;
                ///
                size_t idx = m_idxs[ii];
                size_t nextidx = m_idxs[ii + 1];
                auto &bw = p_bws[idx];
                auto &nextbw = p_bws[nextidx];
                if ((bw->getPars()->end_evolution) || (nextbw->getPars()->end_evolution)){
                    evalShellThicknessIsolated(idx, Y);
                }
//            double r_cur = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            double r_cur = bw->getVal(DynRadBlastWave::Q::iR, 0);
                r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
                r_ip1 = Y[nextbw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
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
                m_delta[ii] = dr_i;
            }
        }
    }

    /// Evaluate the radial extend of a velocity shell. Assume ordered shells. Assumes sorted shells
    void evalShellOptDepth(  const double * Y ){

        double r_i=0., r_ip1=0., dr_i = 0., m_i=0., m2_i=0., m_ip1=0.,
                m2_ip1=0., vol_i=0., rho_i=0., dtau_i=0., tau_tot=0., eint2_i=0., ene_th=0.;
        double tdiff=0., lum=0.;
        double kappa_i = 0.;/// = bw->; // bw->getPars()->kappa0; // TODO!

        /// Compute temperature of each of the active shells
        if (n_active_shells == 1){
            auto & bw = p_bws[0];
            dr_i = bw->getPars()->thickness;//frac * r_i; // TODO add other methods 1/Gamma...
//            m_i = bw->getPars()->M0;//frac * r_i; // TODO add other methods 1/Gamma...
            vol_i = bw->getPars()->volume;// = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
//            kappa_i = bw->getPars()->kappa;
            ///store also velocity
            m_beta[0] = EQS::BetFromMom( Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::imom] );
//            r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
            ene_th = eint2_i; // TODO is this shock??
            m_temp[0] = std::pow(ene_th/A_RAD/vol_i,0.25); // Stephan-Boltzman law
            int x = 1;
        }
        else {
            for (size_t i = 0; i < n_active_shells; i++) {
                size_t idx = m_idxs[i];
                auto &bw = p_bws[idx];
                dr_i = bw->getPars()->thickness;//frac * r_i; // TODO add other methods 1/Gamma...
                vol_i = bw->getPars()->volume;// = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
                ///store also velocity
                m_beta[i] = EQS::BetFromMom(Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::imom]);
                eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
                eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
                ene_th = eint2_i; // TODO is this shock??
                m_temp[i] = std::pow(ene_th / A_RAD / vol_i, 0.25); // Stephan-Boltzman law
                int x = 1;
            }
        }



        /// if there is one shell we cannot have the shell width that comes from shell separation.
        if (n_active_shells == 1){
            double frac = 0.1; // fraction of the shell volume to be used as its width. Inaccurate as we need adjacent shells to get the volume...
            auto & bw = p_bws[0];
            kappa_i = bw->getPars()->kappa;
//            r_i =  Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            dr_i = bw->getPars()->thickness;//frac * r_i; // TODO add other methods 1/Gamma...
            vol_i = bw->getPars()->volume;// = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) / bw->getPars()->ncells;
            m_i = bw->getPars()->M0;
            m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_i;
            m_i += m2_i;
            rho_i = m_i / vol_i;
            //
            r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
//            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
//            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
//            ene_th = eint2_i; // TODO is this shock??
            //
            m_dtau[0] = kappa_i * rho_i * dr_i;
            m_dtau_cum[0] = 0.; // opt. depth to this shell
            m_tdiff[0] = kappa_i * m_i / m_beta[0] / r_i / CGS::c; // TODO M0 or M0+M2 )
            m_tlc[0] = r_i / CGS::c;

        }

        /// compute shell width from the shell separation
        for (size_t ii=0; ii<n_active_shells-1; ii++){
            r_i=0., r_ip1=0., dr_i = 0., m_i=0., m2_ip1=0., m_ip1=0., vol_i=0., rho_i=0., dtau_i=0., eint2_i=0;
            ///
            size_t idx = m_idxs[ii];
            size_t nextidx = m_idxs[ii+1];
            auto & bw = p_bws[idx];
            auto & nextbw = p_bws[nextidx];

            dr_i = bw->getPars()->thickness;//r_ip1 - r_i;
            kappa_i = bw->getPars()->kappa;

            /// evaluate mass within the shell
            m_i = bw->getPars()->M0;
            m2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_i;
            m_i += m2_i;
            m_ip1 = nextbw->getPars()->M0;
            m2_ip1 = Y[nextbw->getPars()->ii_eq + DynRadBlastWave::QS::iM2] * m_ip1;
            m_ip1 += m2_ip1;
            r_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            /// evaluate the volume of the shell (fraction of the 4pi)
            vol_i = bw->getPars()->volume;//(4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i) / bw->getPars()->ncells;
            /// evaluate density within a shell (neglecting the accreted by the shock!)
            rho_i = m_i / vol_i;
            m_rho[ii] = rho_i;

//            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
//            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different

            /// evaluate optical depth
            dtau_i = kappa_i * rho_i * dr_i;
            tau_tot += dtau_i;
            if ((!std::isfinite(dtau_i)) || (dtau_i < 0)){
                (*p_log)(LOG_ERR,AT) << "dtau is nan or < 0 and dtau="<<dtau_i<<"\n";
                exit(1);
            }
            m_dtau[ii] = dtau_i;

            m_tlc[ii] = r_i / CGS::c;
            m_tdiff[ii] = kappa_i * bw->getPars()->M0 / m_beta[ii] / r_i / CGS::c; // TODO M0 or M0+M2 )
//            tdiff = m_dtau[ii] * Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR] / CGS::c;

        }

        /// compute diffusion out
//        for (size_t ii = 0 ; ii < n_active_shells; ii++){
//            m_tdiff_out[i]
//        }


//        (*p_log)(LOG_INFO,AT)<<" optical depth is evaluated. Total tau="<<tau_tot<<"\n";
//        print_xy_as_numpy<Vector>(m_radii_init,m_rho,m_rho.size(),1);
//        std::cout << " rho_ave="<<getShellRho(Y)<<"\n";
//        std::cout << " m_tot="<<getShellMass(Y)<<"\n";


//        (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<mylayer
//                                        <<"] Photosphere located at idx="<<idx_photo<<"\n";


        /// Compute the optical depth from 0 to a given shell
        for (size_t ii = 0 ; ii < n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            double tau_cum = 0.;
            for (size_t jj = 0 ; jj < n_active_shells; jj++){
                size_t other_idx = m_idxs[jj];
                if (other_idx < cur_idx) {
                    tau_cum += m_dtau[other_idx];
                }
            }
            m_dtau_cum[cur_idx] = tau_cum;
        }

        /// Compute cumulative diffusive timescale for all shells
        for (size_t ii = 0 ; ii < n_active_shells; ii++){
            double tdiff_cum = 0.;
            for (size_t jj = ii ; jj < n_active_shells; jj++){
                tdiff_cum += m_tdiff[jj];
            }
            m_tdiff_out[ii] = std::max( tdiff_cum, m_tlc[ii] ); // Eq. 19 in Ren+19
        }

        /// Compute luminocity from each shell
        double ltot = 0.;
        for (size_t ii=0; ii<n_active_shells-1; ii++){
            size_t idx = m_idxs[ii];
            auto & bw = p_bws[idx];
            eint2_i = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iEint2];
            eint2_i *= bw->getPars()->M0 * CGS::c * CGS::c; /// remember the units in ODE solver are different
            ene_th = eint2_i; // TODO is this shock??
            m_lum[ii] = ene_th / m_tdiff_out[ii]; // m_dtau[ii] * ene_th * bw->getPars()->M0 * CGS::c * CGS::c / tdiff;
            ltot += m_lum[ii];
        }
        m_tot_lum = ltot;


        /// based on where tau = 1 get the photosphere radius, i.e., in which shell is the photosphere
        double tmp_tau = 0; double r_ph;
        size_t idx_photo = n_active_shells;
        for (size_t i = n_active_shells-1; i >= 0; --i){
            size_t idx = m_idxs[i];
            auto & bw = p_bws[idx];
            tmp_tau+=m_dtau[idx];
            if (tmp_tau > 1.){
                idx_photo = (int)idx+1;
                break;
            }
        }
        if (idx_photo == 0) {
            /// Ejecta totally transparent
            int y = 1;
        }
        else if ((idx_photo == n_active_shells-1) && (m_dtau[idx_photo-1] > 1)) {
            /// Ejecta totally opaque
            int x = 1;
        }
        else{
            m_photo_r = Y[p_bws[idx_photo-1]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            m_photo_teff = 0.;
        }




        std::cout << "m_rho      ="<< m_rho << "\n";
        std::cout << "m_lum      ="<< m_lum << "\n";
        std::cout << "m_temp     ="<< m_temp << "\n";
        std::cout << "m_dtau     ="<< m_dtau << "\n";
        std::cout << "m_tdiff    ="<< m_tdiff << "\n";
        std::cout << "m_tdiff_out="<< m_tdiff_out << "\n";
        std::cout << "m_dtau_cum ="<< m_dtau_cum << "\n";


//        (*p_log)(LOG_INFO,AT)<<"\tLayer [il="<<mylayer
//                                        <<"] Cumulative optical depth to idx0"
//                                        <<" tau0="<<m_dtau_cum[0]
//                                        <<" tau[-1]="<<m_dtau_cum[n_active_shells-1]<<"\n";

        /// save results for the use in ODE when solving Energy equation
        for (size_t ii = 0 ; ii < n_active_shells; ii++){
            size_t cur_idx = m_idxs[ii];
            auto & bw = p_bws[cur_idx];
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
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::idR) = dr_i;
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::iRho) = rho_i;
//        p_bws[ish]->getLastVal(DynRadBlastWave::Q::iTau) = tau_i; // update blastwave thickness

//        if (idx_photo == m_idxs.back()-1)
//            idx_photo = (int)m_idxs.back();
//        for (auto idx : m_idxs)
//            std::cout << "i_photo=" <<idx_photo<< "; "<< p_bws[idx]->getLastVal(DynRadBlastWave::Q::iTau) << ", ";
//        std::cout <<"\n";

//        double tau_out = 0, tau_out_i = 0;
//        for (auto idx_curr : m_idxs){
//            tau_out = 0.;
//            for (auto idx = m_idxs.back(); idx > idx_curr; idx--){
//                tau_out_i = p_bws[idx]->getLastVal(DynRadBlastWave::Q::iTau);
//                tau_out+=tau_out_i ;


    }
    /// Evaluate shell total mass, volume, density
    /// next three functions compute mass, volume and average density of the entire shell
    double getShellMass(const double * Y_){
        double mtot = 0.;
        for (size_t ii=0; ii<n_active_shells; ++ii){
            size_t i_idx = m_idxs[ii];
            double m2 = Y_[p_bws[i_idx]->getPars()->ii_eq + DynRadBlastWave::QS::iM2];
            double m0 = p_bws[i_idx]->getPars()->M0;
            double m2plus0 = (1. + m2) * m0;
            mtot+=m2plus0;
        }
        if (!std::isfinite(mtot) || (mtot < 0)){
            (*p_log)(LOG_ERR,AT) << "mtot is nan or < 0; mtot="<<mtot<<"\n";
            exit(1);
        }
        return mtot;
    }
    double getShellVolume(const double * Y){
        double r0 = Y[p_bws[ m_idxs[0] ]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
        double r1 = Y[p_bws[ m_idxs[n_active_shells-1] ]->getPars()->ii_eq + DynRadBlastWave::QS::iR];
        if ((r0 >= r1)||(r0==0)||(r1==0)){
            (*p_log)(LOG_ERR,AT)<<" r0 > r1. in the shell; r0="<<r0<<" r1="<<r1<<"\n";
            exit(1);
        }
        double delta = r1-r0;
        double volume = (4./3.) * CGS::pi * (r1*r1*r1 - r0*r0*r0) / p_bws[ m_idxs[0] ]->getPars()->ncells;
    }
    double getShellRho(const double * Y){
        double mtot = getShellMass(Y);
        double volume = getShellVolume(Y);
        return mtot / volume;
    }
    /// evaluate how much each shell in the cumShell gets based on cumulative optical depth

#if 0
    inline std::vector<size_t> & getCurrentIndexes( ) { return m_idxs; }
    void evaluateInitialShellLength(){

    }

    void evaluateOptDepthAtCurrR(size_t i_cur, double x, const double * Y){

        double kappa = 10.;

        /// evaluate relative positions of shells at this time
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws[i]->getLastVal(std::static_pointer_cast<BW::Q>())
        }
        /// get indexes of sorted shells by radius
        sort_indexes(m_radii, m_idxs);
        /// evaluate opacity of each shell up to current
        double tau = 0.;
        for (auto & idx : m_idxs){
            auto & bw = p_bws[idx];
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
        auto & bw_cur = p_bws[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::QS::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

    }





    inline void evalRelativePosition(double x, const double * Y){
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::QS::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws[i]->getLastVal(std::static_pointer_cast<BW::Q>())
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
            auto & bw = p_bws[i];
            bw->getPars()->i_position_in_layer = current_indexes[i];
        }

        std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws;
        for (size_t & idx : current_indexes){
            p_sorted_bws.emplace_back(p_bws[idx]);
        }

    }
    ///
    void evaluateDensityProfileFromBWpositions(){
        for (size_t & idx : current_indexes) {
            auto & bw = p_bws[idx];

            double rho_i = bw->getPars()->M0 / (4./3.*CGS::pi*())
        }

        Vector rho ( m_nshells, 0.0 );
        for (size_t i = 0; i < m_nshells; i++){
            size_t idx = current_indexes[i];

        }
    }

    void evaluateOpticalDepthToMagnetarEmission(size_t i_cur, double x, const double * Y){
        double tau = 0;
        auto & bw_cur = p_bws[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::QS::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

        double tau = 0;
        auto & bw_cur = p_bws[i_cur];
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws[i];
            if (r_cur > bw->)
        }

        for (size_t i = 0; i < i_position; i++){
            tau +=
        }
    }

#endif
};
#endif

#endif //SRC_MODEL_EJECTA_H
