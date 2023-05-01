//
// Created by vsevolod on 02/04/23.
//

#ifndef SRC_MODEL_EJECTA_H
#define SRC_MODEL_EJECTA_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"
#include "synchrotron_an.h"

#include "blastwave_components.h"
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
        unsigned m_loglevel = 0;
    };
    std::vector<size_t> m_idxs{};
    std::unique_ptr<logger> p_log = nullptr;
    std::unique_ptr<Pars> p_pars = nullptr;
    std::unique_ptr<BlastWaveCollision> p_coll = nullptr;
    std::vector<std::unique_ptr<BlastWave>> p_bws{};
    VecVector m_data{};
//    VecVector m_data_prev{};
public:
    enum Q { irEj, irhoEj, ibetaEj, idelta, ivol, idtau, itaucum, itaucum0, ieth, itemp, ilum, itdiff };
    std::vector<std::string> m_vnames{
            "r", "rho", "beta", "delta", "vol", "dtau", "taucum", "taucum0", "eth", "temp", "lum", "tdiff"
    };
    CumulativeShell(Vector t_grid, size_t nshells, int ilayer, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CumulativeShell");
        p_coll = std::make_unique<BlastWaveCollision>(loglevel);
        p_pars = std::make_unique<Pars>();
        p_pars->ilayer=ilayer;
        p_pars->nshells=nshells;
        p_pars->n_active_shells=nshells;
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws.emplace_back(std::make_unique<BlastWave>(t_grid, ishell, ilayer, loglevel ) );
        if (t_grid.empty())
            return;
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize(nshells);
//        for (auto & arr : m_data_prev)
//            arr.resize(nshells);
        m_idxs.resize(nshells);
        std::fill(m_data[Q::irEj].begin(), m_data[Q::irEj].end(), std::numeric_limits<double>::max());
        std::fill(m_idxs.begin(), m_idxs.end(), std::numeric_limits<size_t>::max());
        p_pars->m_loglevel=loglevel;
    }
    void setPars(StrDbMap & pars, StrStrMap & opts){
        p_pars->do_thermrad_loss= getBoolOpt("do_thermrad_loss",opts,AT,p_log,false,true);
        if (!p_pars->do_thermrad_loss)
            return;
        p_pars->thermradloss_at_photo= getBoolOpt("thermradloss_at_photo",opts,AT,p_log,false,true);
    }
    /// copy current shell properties to each blastwave "data container"
    void insertStatusInBWdata(size_t it){
//        if (p_pars->)

        for (size_t i=0; i<p_pars->n_active_shells; i++){
            size_t idx = m_idxs[i];
            if (p_bws[idx]->getPars()->end_evolution) {
                continue;
            }
            auto & bw1 = p_bws[idx];
//            for (size_t key = 0; key < m_vnames.size(); key++)
//                bw1->getData(static_cast<BW::Q>(key))[it] = m_data[key][i];

            bw1->getData()[BW::Q::iEJr][it]       = m_data[Q::irEj][i];
            bw1->getData()[BW::Q::iEJbeta][it]    = m_data[Q::ibetaEj][i];
            bw1->getData()[BW::Q::iEJdelta][it]   = m_data[Q::idelta][i];
            bw1->getData()[BW::Q::iEJdtau][it]    = m_data[Q::idtau][i];
            bw1->getData()[BW::Q::iEJeth][it]     = m_data[Q::ieth][i];
            bw1->getData()[BW::Q::iEJlum][it]     = m_data[Q::ilum][i];
            bw1->getData()[BW::Q::iEJrho][it]     = m_data[Q::irhoEj][i];
            bw1->getData()[BW::Q::iEJtaucum][it]  = m_data[Q::itaucum][i];
            bw1->getData()[BW::Q::iEJtdiff][it]   = m_data[Q::itdiff][i];
            bw1->getData()[BW::Q::iEJtemp][it]    = m_data[Q::itemp][i];
            bw1->getData()[BW::Q::iEJtaucum0][it] = m_data[Q::itaucum0][i];
            bw1->getData()[BW::Q::iEJvol][it]     = m_data[Q::ivol][i];
        }
    }
    // -------------------------
    inline double getR(size_t i){return m_data[Q::irEj][i];}
    inline Vector & getRvec(){return m_data[Q::irEj]; }
    inline std::vector<size_t> & getIdx(){ return m_idxs; }
    inline Vector & getBetaVec(){return m_data[Q::ibetaEj]; }
    inline Vector & getRhoVec(){return m_data[Q::irhoEj]; }
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
        std::fill(m_data[Q::irEj].begin(), m_data[Q::irEj].end(),std::numeric_limits<double>::max());
        bool is_sorted = true;
        for (size_t i=0; i<p_pars->n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            m_data[Q::irEj][i] = Y[ p_bws[idx]->getPars()->ii_eq + SOL::QS::iR ];
        }
        for (size_t i=1; i<p_pars->n_active_shells; ++i) {
            size_t idx = m_idxs[i];
            double rim1 = m_data[Q::irEj][i-1];
            /// if the next shells in the list has a radius > than the previous; shells not ordered
            if (m_data[Q::irEj][i] > rim1){
                rim1 = m_data[Q::irEj][i];
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
            double r_i = m_data[Q::irEj][ii];
            if (r_i == std::numeric_limits<double>::max()){
                continue;
            }
            std::vector<size_t> overrun_shells;
            Vector t_at_which_overrun;
            /// now check which blast-waves has this one overrun

            for (size_t jj = ii; jj < p_pars->n_active_shells; ++jj) {
                size_t j_idx = m_idxs[jj];
                double r_j = m_data[Q::irEj][jj];
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
#if 0
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
#endif
    void updateSortedShellWidth_( const double * Y ){

        /// for the first shell we have to assume it width
#if 0
        size_t idx = m_idxs[0];//p_pars->n_active_shells-1;
        double frac = 0.5; // TODO fraction of the shell volume to be used as its width. !!! Inaccurate as we need adjacent shells to get the volume...
        auto & bw = p_bws[idx];
        double r_i =  Y[bw->getPars()->ii_eq + SOL::QS::iR];
        double dr_i = frac * r_i;
        double vol_i = frac * (4./3.) * CGS::pi * (r_i*r_i*r_i) ;/// bw->getPars()->ncells;
        if (dr_i <= 0 or !std::isfinite(dr_i)){
            (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
            exit(1);
        }
        m_data[Q::idelta][0] = dr_i;
        m_data[Q::ivol][0] = vol_i;
#endif
        if (p_pars->n_active_shells > 1) {
            size_t idx = m_idxs[0];//p_pars->n_active_shells-1;
            size_t nextidx = m_idxs[1];//p_pars->n_active_shells-1;
            auto & bw = p_bws[idx];
            auto & nextbw = p_bws[nextidx];

            double r_i = Y[bw->getPars()->ii_eq + SOL::QS::iR];
            double r_ip1 = Y[nextbw->getPars()->ii_eq + SOL::QS::iR];
            double dr_i = r_ip1 - r_i;
            double vol_i =
                    (4. / 3.) * CGS::pi * (r_ip1 * r_ip1 * r_ip1 - r_i * r_i * r_i);// / currbw->getPars()->ncells;
            if (dr_i <= 0 or !std::isfinite(dr_i)) {
                (*p_log)(LOG_ERR, AT) << " dr_i = " << dr_i << "\n";
                exit(1);
            }
            m_data[Q::idelta][0] = dr_i;
            m_data[Q::ivol][0] = vol_i;
            bw->getPars()->delta = dr_i;
            bw->getPars()->vol = vol_i;
        }
        else{
            double koeff = 0.1;
            size_t idx = m_idxs[0];//p_pars->n_active_shells-1;
            auto &bw = p_bws[idx];
            double r_i = Y[bw->getPars()->ii_eq + SOL::QS::iR];

            double _r_i = r_i * koeff;
            double _vol_i = koeff * (4./3.) * CGS::pi * (r_i*r_i*r_i);
            double _r_i_pref = bw->getPars()->delta;
            double _vol_i_prev = bw->getPars()->vol;
//            double k_r =
//            size_t idx = m_idxs[0];
            double _r_ip1 = 1.01 * r_i;
            m_data[Q::idelta][0] = _r_ip1 - r_i;// p_bws[idx]->getPars()->delta;
            m_data[Q::ivol][0] = (4./3.) * CGS::pi * (_r_ip1*_r_ip1*_r_ip1 - r_i*r_i*r_i); //p_bws[idx]->getPars()->vol;
        }
        for (size_t ii=1; ii<p_pars->n_active_shells; ii++) {
            ///
            size_t previdx = m_idxs[ii-1];
            size_t curidx = m_idxs[ii];
            auto & prevbw = p_bws[previdx];
            auto & currbw = p_bws[curidx];
            if ((prevbw->getPars()->end_evolution) || (currbw->getPars()->end_evolution)){
//                    evalShellThicknessIsolated(idx, Y);
                (*p_log)(LOG_ERR,AT) << "|error|\n";
                exit(1);
            }
            double rm1_i = Y[prevbw->getPars()->ii_eq + SOL::QS::iR];
            double r_i = Y[currbw->getPars()->ii_eq + SOL::QS::iR];
            if ((rm1_i == 0.) || (r_i == 0.)) {
//                (*p_log)(LOG_WARN,AT)<<" shell="<<idx<<" has r_i="<<r_i<<" and r_ip1="<<r_ip1<<"\n";
                continue;
            }
            double dr_i = r_i - rm1_i;
            /// evaluate the volume of the shell (fraction of the 4pi)
            double vol_i = (4./3.) * CGS::pi * (r_i*r_i*r_i - rm1_i*rm1_i*rm1_i);// / currbw->getPars()->ncells;
            if ((dr_i <= 0) or (!std::isfinite(dr_i))){
                (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
                exit(1);
            }
            /// --------------------------- |
            m_data[Q::idelta][ii] = dr_i;
            m_data[Q::ivol][ii] = vol_i;
            prevbw->getPars()->delta = dr_i;
            prevbw->getPars()->vol = vol_i;
        }
    }
    void updateSortedShellWidth( const double * Y ){
        /// -------------------------------------------------------
        for (size_t ii=0; ii<p_pars->n_active_shells-1; ii++) {
            ///
            size_t nextidx = m_idxs[ii+1];
            size_t curidx = m_idxs[ii];
            auto & nextbw = p_bws[nextidx];
            auto & currbw = p_bws[curidx];
            if ((nextbw->getPars()->end_evolution) || (currbw->getPars()->end_evolution)){
//                    evalShellThicknessIsolated(idx, Y);
                (*p_log)(LOG_ERR,AT) << "|error|\n";
                exit(1);
            }
            double r_ip1 = Y[nextbw->getPars()->ii_eq + SOL::QS::iR];
            double r_i = Y[currbw->getPars()->ii_eq + SOL::QS::iR];
            if ((r_ip1 == 0.) || (r_i == 0.)) {
//                (*p_log)(LOG_WARN,AT)<<" shell="<<idx<<" has r_i="<<r_i<<" and r_ip1="<<r_ip1<<"\n";
                continue;
            }
            double dr_i = r_ip1 - r_i;
            /// evaluate the volume of the shell (fraction of the 4pi)
            double vol_i = (4./3.) * CGS::pi * (r_ip1*r_ip1*r_ip1 - r_i*r_i*r_i);// / currbw->getPars()->ncells;
            if ((dr_i <= 0) or (!std::isfinite(dr_i))){
                (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
                exit(1);
            }
            /// --------------------------- |
            m_data[Q::idelta][ii] = dr_i;
            m_data[Q::ivol][ii] = vol_i;
            currbw->getPars()->delta = dr_i;
            currbw->getPars()->vol = vol_i;
            /// ---------------------------- | for when this shell is the only one left |------
            double _r_i = Y[currbw->getPars()->ii_eq + SOL::QS::iR];
            double _vol_i = (4./3.) * CGS::pi * (_r_i*_r_i*_r_i);
            currbw->getPars()->_last_frac = _r_i / m_data[Q::idelta][ii];
            currbw->getPars()->_last_frac_vol = _vol_i / m_data[Q::ivol][ii];
        }

        /// ----------------------------
        if (p_pars->n_active_shells>1){
            /// copy from previos shell
            size_t ii = p_pars->n_active_shells-1;
            size_t idx = m_idxs[ii];
            p_bws[idx]->getPars()->delta = 1.01 * p_bws[idx-1]->getPars()->delta;
            p_bws[idx]->getPars()->vol = 1.01 * p_bws[idx-1]->getPars()->vol;
            m_data[Q::idelta][ii] = p_bws[idx]->getPars()->delta;//dr_i;
            m_data[Q::ivol][ii] = p_bws[idx]->getPars()->vol;
            double _r_i = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR];
            double _vol_i = (4./3.) * CGS::pi * (_r_i*_r_i*_r_i);
            p_bws[idx]->getPars()->_last_frac = _r_i / m_data[Q::idelta][ii];
            p_bws[idx]->getPars()->_last_frac_vol = _vol_i / m_data[Q::ivol][ii];
        }
        else {
            // use same fraction as at the last time there were many shells
            size_t ii = p_pars->n_active_shells-1;
            size_t idx = m_idxs[ii];
//            auto & curbw = p_bws[p_pars->n_active_shells-1];
//            if (curbw->getPars()->delta<=0){
//                (*p_log)(LOG_ERR,AT)<<" previous dr <= 0 dr="<<curbw->getPars()->delta<<"\n";
//                exit(1);
//            }
//            if (p_bws[idx]->getPars()->_last_frac <= 0. || p_bws[idx]->getPars()->_last_frac_vol <= 0.){
//                (*p_log)(LOG_ERR,AT)<<" frac <= 0 frac="<<curbw->getPars()->_last_frac
//                                               <<" frac_vol <= 0 frac="<<curbw->getPars()->_last_frac_vol<<"\n";
//                exit(1);
//            }
//            double _r_i = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR];
//            double _vol_i = (4./3.) * CGS::pi * (_r_i*_r_i*_r_i);
//            p_bws[idx]->getPars()->delta = _r_i * p_bws[idx]->getPars()->_last_frac;
//            p_bws[idx]->getPars()->vol = _vol_i * p_bws[idx]->getPars()->_last_frac_vol;
            m_data[Q::idelta][ii] = p_bws[idx]->getPars()->delta;//dr_i;
            m_data[Q::ivol][ii] = p_bws[idx]->getPars()->vol;
        }

//
//        if (p_pars->n_active_shells>1){
//            double frac = 0.5;
//            size_t ii = p_pars->n_active_shells-1;
//            size_t idx = m_idxs[ii];
//            size_t previdx = m_idxs[ii-1];
//            auto & prevbw = p_bws[previdx];
//            auto & curbw = p_bws[idx];
//            if ((prevbw->getPars()->end_evolution) || (curbw->getPars()->end_evolution)){
////                    evalShellThicknessIsolated(idx, Y);
//                (*p_log)(LOG_ERR,AT) << "|error|\n";
//                exit(1);
//            }
//            double prev_r = Y[prevbw->getPars()->ii_eq + SOL::QS::iR];
//            double curbw_r = Y[curbw->getPars()->ii_eq + SOL::QS::iR];
//            double dr_i = curbw_r - prev_r;
//            if ((dr_i <= 0) or (!std::isfinite(dr_i))){
//                (*p_log)(LOG_ERR,AT)<<" dr_i = "<<dr_i<<"\n";
//                exit(1);
//            }
//            double vol_i = (4./3.) * CGS::pi * (curbw_r*curbw_r*curbw_r - prev_r*prev_r*prev_r);
//            m_data[Q::idelta][ii] = dr_i;
//            m_data[Q::ivol][ii] = vol_i;
//            curbw->getPars()->delta = dr_i;
//            curbw->getPars()->vol = vol_i;
//        }
//        else{
//            size_t ii = p_pars->n_active_shells-1;
//            auto & curbw = p_bws[p_pars->n_active_shells-1];
//            if (curbw->getPars()->delta<=0){
//                (*p_log)(LOG_ERR,AT)<<" previous dr <= 0 dr="<<curbw->getPars()->delta<<"\n";
//                exit(1);
//            }
////            double r_i = Y[curbw->getPars()->ii_eq + SOL::QS::iR];
////            double vol_i = (4./3.) * CGS::pi * (r_i*r_i*r_i);
////            double koeff = vol_i / curbw->getPars()->vol;
////            double koeff_dr = std::pow(koeff, 1./3.);
////            if (koeff < 1.){
////                (*p_log)(LOG_ERR,AT)<<" koeff < 1; koeff="<<koeff<<"\n";
////                exit(1);
////            }
//
//
//        }
        int x = 1;
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
            m_data[Q::ibetaEj][i] = m_beta;

            double m_i = bw->getPars()->M0;
            double m2_i = Y[bw->getPars()->ii_eq + SOL::QS::iM2] * m_i;
//            m_i += m2_i; // TODO should I include this?/
            double rho_i = m_i / vol_i;
            double kappa_i = bw->getPars()->kappa;
            m_data[Q::irhoEj][i] = rho_i;
            m_data[Q::idtau][i] = kappa_i * rho_i * dr_i;
            tau_tot += m_data[Q::idtau][i];

            double r_i = Y[bw->getPars()->ii_eq + SOL::QS::iR];
            double m_tlc = r_i / CGS::c;
            m_data[Q::itdiff][i] = std::max( kappa_i * m_i / m_data[Q::ibetaEj][i] / r_i / CGS::c, m_tlc); // TODO M0 or M0+M2 )

            m_data[Q::ilum][i] = ene_th / m_data[Q::itdiff][i];
            if (m_data[Q::itdiff][i] < 1){
                (*p_log)(LOG_WARN,AT) << " small tdiff\n";
            }
            ltot += m_data[Q::ilum][i];
        }
        p_pars->ltot = ltot;
        p_pars->tautot = tau_tot;
        if (tau_tot < 1.){
            int x = 1;
        }

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

        int x = 1;


    }
    /// Evaluate shell total mass, volume, density
    /// next three functions evaluateShycnhrotronSpectrum mass, volume and average density of the entire shell
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
public:
    bool do_eninj_inside_rhs = false;
    bool run_bws=false, save_dyn=false, load_dyn=false, do_ele=false, do_spec=false, do_lc=false, do_skymap=false;
    bool do_collision = false;
    bool do_nuc = false;
    bool is_ejecta_obs_pars_set = false;
    bool is_ejecta_anal_synch_computed = false;
    StrDbMap m_pars; StrStrMap m_opts;
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
    VecVector & getData(size_t il, size_t ish){ return getShells()[il]->getBW(ish)->getData(); }
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
//        return run_bws ? ejectaStructs.structs[0].m_nlayers : 0;
        return (run_bws||load_dyn) ? id->nlayers : 0;
    }
    size_t nshells() const {
//        return run_bws ? ejectaStructs.nshells : 0;
        return (run_bws||load_dyn) ? id->nshells : 0;
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
        return (run_bws || load_dyn) ? (int)id->ncells : 0;
        int x = 1;
    }
    std::vector<std::unique_ptr<CumulativeShell>> & getShells(){
        if (p_cumShells.empty()){
            (*p_log)(LOG_ERR,AT)<<" ejecta not initialized\n";
            exit(1);
        }
        return p_cumShells;
    }
    std::unique_ptr<EjectaID2> & getId(){ return id; }

    void setPars(StrDbMap & pars, StrStrMap & opts,
                 std::string working_dir_, std::string parfilename_,
                 StrDbMap & main_pars, size_t ii_eq, size_t iljet){
        working_dir = working_dir_;
        parfilename = parfilename_;
        m_pars = pars;
        m_opts = opts;
        /// read GRB afterglow parameters
//        StrDbMap m_pars; StrStrMap m_opts;
//        run_bws=false; bool save_dyn=false, do_ele=false, do_spec=false, do_lc=false, do_skymap=false;
        if ((!m_pars.empty()) || (!m_opts.empty())) {
            run_bws = getBoolOpt("run_bws", m_opts, AT, p_log, false, true);
            save_dyn = getBoolOpt("save_dynamics", m_opts, AT, p_log, false, true);
            load_dyn = getBoolOpt("load_dynamics", m_opts, AT, p_log, false, true);
            do_ele = getBoolOpt("do_ele", m_opts, AT, p_log, false, true);
            do_spec = getBoolOpt("do_spec", m_opts, AT, p_log, false, true);
            do_lc = getBoolOpt("do_lc", m_opts, AT, p_log, false, true);
            do_skymap = getBoolOpt("do_skymap", m_opts, AT, p_log, false, true);
            for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
                if (main_pars.find(key) == main_pars.end()) {
                    (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                    exit(1);
                }
                m_pars[key] = main_pars.at(key);
            }
            m_opts["workingdir"] = working_dir; // For loading Nuclear Heating table
            if (run_bws || load_dyn) {
                std::string fname_ejecta_id = getStrOpt("fname_ejecta_id", m_opts, AT, p_log, "", true);
                bool use_1d_id = getBoolOpt("use_1d_id", m_opts, AT, p_log, false, true);
                if (!std::experimental::filesystem::exists(working_dir + fname_ejecta_id)) {
                    (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_ejecta_id << "\n";
                    exit(1);
                }
                id = std::make_unique<EjectaID2>(
                        working_dir + fname_ejecta_id,
                        getStrOpt("method_eats", m_opts, AT, p_log, "", true),
                        getBoolOpt("use_1d_id", m_opts, AT, p_log, false, true),
                        getBoolOpt("load_r0", m_opts, AT, p_log, false, true),
                        t_arr[0],
                        m_loglevel );
                setEjectaBwPars(m_pars, m_opts, ii_eq, iljet);
            }
        }
        else{
            (*p_log)(LOG_INFO, AT) << "ejecta is not initialized and will not be considered.\n";
        }
    }

    void computeAndOutputObservables(StrDbMap & main_pars, StrStrMap & main_opts){
        /// work on GRB afterglow
        if (run_bws || load_dyn){
            bool lc_freq_to_time = getBoolOpt("lc_use_freq_to_time",main_opts,AT,p_log,false,true);
            Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
            Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

            if (save_dyn)
                saveEjectaBWsDynamics(
                        working_dir,
                        getStrOpt("fname_dyn", m_opts, AT, p_log, "", true),
                        (int)getDoublePar("save_dyn_every_it", m_pars, AT, p_log, 1, true),
                        main_pars, m_pars);
            if (load_dyn)
                loadEjectaBWDynamics(working_dir,
                                     getStrOpt("fname_dyn", m_opts, AT, p_log, "", true));

            if (do_ele)
                setPreComputeEjectaAnalyticElectronsPars();
//                (*p_log)(LOG_INFO, AT) << "jet analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";

            if (do_lc) {
                computeSaveEjectaLightCurveAnalytic(
                        working_dir,
                        getStrOpt("fname_light_curve", m_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", m_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, m_pars, lc_freq_to_time);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_skymap)
                computeSaveEjectaSkyImagesAnalytic(
                        working_dir,
                        getStrOpt("fname_sky_map", m_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, m_pars);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";

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
//            for (size_t il = 0; il < ejectaStructs.structs[0].m_nlayers; il++){
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
    void setEjectaBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq, size_t n_layers_jet){

        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        load_dyn = getBoolOpt("load_dynamics", opts, AT, p_log, false, true);
        if ((!run_bws) && (!load_dyn))
            return;

        if ((!run_bws) && (load_dyn)){
            is_ejecta_obs_pars_set = true;
            for(size_t il = 0; il < nlayers(); il++) {
                p_cumShells.push_back(
                        std::make_unique<CumulativeShell>(Vector {}, nshells(), il,
                                                          p_log->getLogLevel()));
                p_cumShells[il]->setPars(pars, opts);
                for (size_t ish = 0; ish < nshells(); ish++){
                    auto & bw = p_cumShells[il]->getBW(ish);
                    bw->setAllParametersForOneLayer(id, pars, opts, il, ii_eq);
#if 0
                    switch (id->method_eats) {
                        case EjectaID2::iadaptive:
                            bw->getFsEATS()->setEatsPars(
                                    pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                                    0.,0.,id->theta_wing,
                                    getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                            break;
                        case EjectaID2::ipiecewise:
                            bw->getFsEATS()->setEatsPars(
                                    pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                                    id->get(ish,il,EjectaID2::Q::itheta_c_l),
                                    id->get(ish,il,EjectaID2::Q::itheta_c_h),0.,
                                    getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                            break;
                    }
#endif
                    ii_eq += SOL::neq;//bw->getNeq();
                }
            }
            return;
        }


        bool is_within = false;
        std::vector<size_t> which_within{};
        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t nshells_ = nshells();//ejectaStructs.nshells;
        size_t n_layers_ej_ = nlayers();//ejectaStructs.structs[0].m_nlayers;
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
#if 0
                switch (id->method_eats) {
                    case EjectaID2::iadaptive:
                        bw->getRad()->setEatsPars(
                                pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                                0.,0.,id->theta_wing,
                                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                        break;
                    case EjectaID2::ipiecewise:
                        bw->getRad()->setEatsPars(
                                pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                                id->get(ish,il,EjectaID2::Q::itheta_c_l),
                                id->get(ish,il,EjectaID2::Q::itheta_c_h),0.,
                                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                        break;
                }
#endif

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
        do_nuc = getBoolOpt("do_nucinj", opts, AT,p_log, false, true);

        do_eninj_inside_rhs = getBoolOpt("do_eninj_inside_rhs", opts, AT, p_log, "no", false);

        (*p_log)(LOG_INFO,AT) << "finished initializing ejecta. "
                                 "nshells="<<nshells()<<" nlayers="<<nlayers()<<"\n";
    }

    double ej_rtol = 1e-5;

    /// OUTPUT
    void saveEjectaBWsDynamics_old(std::string workingdir, std::string fname, size_t every_it,
                               StrDbMap & main_pars, StrDbMap & ej_pars){
        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

//        size_t nshells = p_cumShells->nshells();
//        size_t m_nlayers = p_cumShells->m_nlayers();
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
//                        {"theta_c_l",bw->getPars()->theta_c_l},
//                        {"theta_c_h",bw->getPars()->theta_c_h},
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
    void saveEjectaBWsDynamics(std::string workingdir, std::string fname, size_t every_it,
                               StrDbMap & main_pars, StrDbMap & ej_pars){
        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = getShells();
        std::vector<std::vector<double>> tot_dyn_out ( nshells() * nlayers() * BW::m_vnames.size() );
        for (auto & arr : tot_dyn_out)
            arr.resize(t_arr.size());
        std::vector<std::string> arr_names{};
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            for(size_t ilayer = 0; ilayer < nlayers(); ilayer++){
                for (size_t ivar = 0; ivar < BW::m_vnames.size(); ivar++) {
                    arr_names.push_back("shell="+std::to_string(ishell)
                                       +" layer="+std::to_string(ilayer)
                                       +" key="+BW::m_vnames[ivar]);
                    auto & bw = models[ilayer]->getBW(ishell);
                    for (size_t it = 0; it < t_arr.size(); it++)
                        tot_dyn_out[ii][it] = bw->getData(static_cast<BW::Q>(ivar))[it];
                    size_t size = tot_dyn_out[ii].size();
                    auto & x = tot_dyn_out[ii];
                    ii++;
                }
            }
        }

        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells() },
                {"nlayers", nlayers() },
                {"ntimes", t_arr.size() },
                {"ncells", ncells() }
        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(tot_dyn_out,arr_names,workingdir+fname,attrs);
    }

    /// INPUT
    void loadEjectaBWDynamics(std::string workingdir, std::string fname){
        if (!std::experimental::filesystem::exists(working_dir+fname))
            throw std::runtime_error("File not found. " + workingdir+fname);

        Exception::dontPrint();
        H5std_string FILE_NAME(workingdir+fname);
        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        size_t nshells_ = (size_t)getDoubleAttr(file,"nshells");
        size_t nlayers_ = (size_t)getDoubleAttr(file, "nlayers");
        size_t ntimes_ = (size_t)getDoubleAttr(file, "ntimes");
        if (nshells_ != nshells()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nshells_="<<nshells_<<" expected nshells="<<nshells()<<"\n";
//            exit(1);
        }
        if (nlayers_ != nlayers()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nlayers_="<<nlayers_<<" expected nlayers_="<<nlayers()<<"\n";
//            exit(1);
        }
//        if (ntimes_ != ()){
//            (*p_log)(LOG_ERR,AT) << "Wring attribute: nlayers_="<<nlayers_<<" expected nlayers_="<<m_nlayers()<<"\n";
//            exit(1);
//        }
        //        double ntimes = getDoubleAttr(file, "ntimes");

        auto & models = getShells();
        for (size_t il = 0; il < nlayers(); il++) {
            size_t n_empty_shells = 0;
            for (size_t ish = 0; ish < nshells(); ish++) {
                auto &bw = models[il]->getBW(ish);
                for (size_t ivar = 0; ivar < BW::m_vnames.size(); ivar++) {
                    std::string key = "shell=" + std::to_string(ish)
                                      + " layer=" + std::to_string(il)
                                      + " key=" + BW::m_vnames[ivar];
                    auto & vec = bw->getData()[static_cast<BW::Q>(ivar)];
                    if (!vec.empty()){
                        (*p_log)(LOG_ERR,AT) << " container is not empty\n";
                    }

                    DataSet dataset = file.openDataSet(key);
                    DataType datatype = dataset.getDataType();
                    DataSpace dataspace = dataset.getSpace();
                    const int npts = dataspace.getSimpleExtentNpoints();

                    H5T_class_t classt = datatype.getClass();
                    if ( classt != 1 )
                    {
                        std::cout << key << " is not a float... you can't save this as a float." << std::endl;
                        exit(1);
                    }
                    FloatType ftype = dataset.getFloatType();
                    H5std_string order_string;
                    H5T_order_t order = ftype.getOrder( order_string);
                    size_t size = ftype.getSize();
//                    vec.resize(1);
                    double * data = new double[npts];
                    if ( order==0 && size == 4 )
                    {
                        std::cout << "NOTE: This is actually float data. We are casting to double" << std:: endl;
                        dataset.read((double*)data, PredType::IEEE_F32LE); // Our standard integer
                    }
                    else if ( order == 0 && size == 8 )
                        dataset.read(data, PredType::IEEE_F64LE);
                    else if ( order == 1 && size == 4 )
                    {
                        std::cout << "NOTE: This is actually float data. We are casting to double" << std:: endl;
                        dataset.read((double*)data, PredType::IEEE_F32BE);
                    }
                    else if ( order ==1 && size == 8 )
                        dataset.read((double*)data, PredType::IEEE_F64BE);
                    else
                        std::cout << "Did not find data type" << std::endl;
                    std::vector<double> v(data, data + npts);
                    vec = std::move( v );
//                    delete[] data;
                    dataspace.close();
                    datatype.close();
                    dataset.close();

                    if ( bw->getData()[static_cast<BW::Q>(ivar)].empty() ){
                        std::cout << key << " faild" << std::endl;
                        exit(1);
                    }
                }
                if (bw->getData()[BW::iR][0] == 0){
//                    (*p_log)(LOG_WARN,AT) << "Loaded not evolved shell [il="<<il<<", "<<"ish="<<ish<<"] \n";
                    n_empty_shells+=1;
                }
                bw->checkEvolution();
                bw->getPars()->nr = bw->getData()[BW::iR].size();
            }
            (*p_log)(LOG_INFO,AT) << "Loaded [il="<<il<<"] N empty shells ="<<n_empty_shells<<"\n";
        }
        file.close();
//        if ( p_cumShells[0]->getBW(0)->getData()[BW::iR][0] == 0 ){
//            std::cout << p_cumShells[0]->getBW(0)->getData()[BW::iR] << "\n";
//            std::cout << " faild" << std::endl;
//            exit(1);
//        }

#if 0
        LoadH5 ldata;
        ldata.setFileName(workingdir+fname);
        ldata.setVarName("nshells");
        double nshells = ldata.getDoubleAttr("nshells");
        double m_nlayers = ldata.getDoubleAttr("m_nlayers");
        auto & models = getShells();
        for (size_t ish = 0; ish < nshells-1; ish++){
            for (size_t il = 0; il < m_nlayers-1; il++){
                auto & bw = models[il]->getBW(ish);
                for (size_t ivar = 0; ivar < BW::m_vnames.size(); ivar++) {
                    std::string key = "shell=" + std::to_string(ish)
                                    + " layer=" + std::to_string(il)
                                    + " key=" + BW::m_vnames[ivar];
                    ldata.setVarName(BW::m_vnames[ivar]);
                    bw->getData().emplace_back( std::move( ldata.getDataVDouble() ) );
                    (*p_log)(LOG_INFO,AT) << "Reading: "<<key<<"\n";
//                    bw->getData(static_cast<BW::Q>(ivar))
//                        = std::move( ldata.getDataVDouble() );
                }

//                auto & bw = models[il]->getBW(ish);
//                std::string
//                bw->getData()[]
            }
        }
#endif
        (*p_log)(LOG_INFO,AT)<<" dynamics loaded successfully\n";

    }

private:

    void computeEjectaSkyMapPieceWise( Images & images, double obs_time, double obs_freq ){

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();

        if (images.empty()){
            (*p_log)(LOG_ERR,AT) << " empty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nshells_){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }

        size_t n_jet_empty_images = 0;

        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        const std::vector<std::string> x {};
        Images tmp (nlayers_, IMG::m_names.size());
        tmp.resize(ncells_);
//        for (auto & _tmp : tmp)
//            _tmp.resize( ncells_ );
        Image tmp_pj( ncells_, IMG::m_names.size(), 0, m_loglevel);
        Image tmp_cj( ncells_, IMG::m_names.size(), 0, m_loglevel);
        for (size_t ishell = 0; ishell < nshells_; ishell++){
//            for (auto & _tmp : tmp)
//                _tmp.clearData();
            tmp.clear();
            tmp_pj.clearData(); tmp_cj.clearData();
            std::vector<size_t> n_empty_images_layer;
            for (size_t ilayer = 0; ilayer < nlayers_; ilayer++){
                /// Evaluate a given image --------------------------------------
                auto & bw_rad = p_cumShells[ilayer]->getBW(ishell)->getFsEATS();
                bw_rad->evalImagePW(tmp.getImgRef(ilayer), tmp_pj, tmp_cj, obs_time, obs_freq);
                /// -------------------------------------------------------------
                if (tmp.getImgRef(ilayer).m_f_tot == 0){
                    n_jet_empty_images += 1;
                    n_empty_images_layer.emplace_back(ilayer);
                }
            }
            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
            combineImages(images.getImgRef(ishell), ncells_, nlayers_, tmp) ;
        }

        /// print which layers/shells gave empty image
        if (p_log->getLogLevel() == LOG_INFO) {
            if (n_jet_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an empty images for total n="
                      << n_jet_empty_images << " layers. Specifically:\n";
                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
//                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                size_t n_layers_ej = m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
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
//        size_t m_nlayers = p_cumShells->m_nlayers();
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
        Image image_i ( ncells(), IMG::m_names.size(), 0, m_loglevel );
        Image im_pj ( ncells(),IMG::m_names.size(), 0, m_loglevel  );
        Image im_cj ( ncells(),IMG::m_names.size(), 0, m_loglevel  );
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            image_i.clearData();
            im_pj.clearData();
            im_cj.clearData();
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStructs.structs[ishell].m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer];//ejectaModels[ishell][ilayer];
//                model->setEatsPars( pars, opts );
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC ntimes="<<obs_times.size()
                        << " vel_shell="<<ishell<<"/"<<nshells()-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers()
                        << " phi_cells="<<EjectaID2::CellsInLayer(ilayer)<<"\n";
                model->getBW(ishell)->getFsEATS()->evalLC(
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

        if ((!run_bws)&&(!load_dyn)){
            std::cerr << " ejecta BWs were not evolved. Cannot evaluateShycnhrotronSpectrum electrons (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = getShells();
        for (auto & model : models) {
            for (auto & bw : model->getBWs()) {
                bw->computeForwardShockElectronAnalyticVars();
                bw->computeForwardShockSynchrotronAnalyticSpectrum();
            }
        }
        is_ejecta_anal_synch_computed = true;
    }

    void updateEjectaObsPars(StrDbMap pars) {

        auto & models = getShells();

//        size_t nshells = nshells();
//        size_t m_nlayers = p_cumShells->m_nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        (*p_log)(LOG_ERR,AT) << "Updating Ejecta observer pars...\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++) {
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = struc.m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer]->getBW(ishell);//ejectaModels[ishell][ilayer];
                model->getFsEATS()->updateObsPars(pars);
                ii++;
            }
        }
//        p_pars->is_ejecta_obs_pars_set = true;
    }

    void evalOptDepthsAlongLineOfSight(double & frac, double mu, double time, double freq){
        size_t nshells_ = nshells();
        Images images(nshells_, IMG_TAU::m_names.size());
        double tau_comp=0.;
        double tau_bh=0.;
        double tau_pf=0.;
        for (size_t il = 0; il < nlayers(); il++){
            auto & cumshell = p_cumShells[il];
            for (size_t ish = 0; ish < cumshell->getPars()->n_active_shells; ish++){

            }
        }

        computeEjectaSkyMapPieceWise( images, time, freq );
        for (size_t i = 0; i < images.size(); i++){
            /// find the layer within a shell that has the smallest 'distance' to the line of sight mu = 0

        }


//        std::vector<double> tau_comp(images.size());
//        std::vector<double> tau_bh(images.size());
//        std::vector<double> tau_comp(images.size());


    }

    void computeSaveEjectaSkyImagesAnalytic(std::string workingdir, std::string fname, Vector times, Vector freqs,
                                            StrDbMap & main_pars, StrDbMap & ej_pars){
        if ((!run_bws)&&(!load_dyn))
            return;

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta sky image with analytic synchrotron...\n";

        if (!is_ejecta_anal_synch_computed){
            std::cerr  << "ejecta analytic electrons were not evolved. Cannot evaluateShycnhrotronSpectrum images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            std::cerr<< "ejecta observer parameters are not set. Cannot evaluateShycnhrotronSpectrum image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
//        size_t ncells_ =  (int)ncells();

//        Image dummy(1,0,);

        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};

        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(IMG::m_names.size());
                for (size_t i_vn = 0; i_vn < IMG::m_names.size(); ++i_vn) {
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
        Images images(nshells_,IMG::m_names.size());

        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            Vector tota_flux(times.size(), 0.0);
            VecVector total_flux_shell( nshells_ );
            for (auto & total_flux_shel : total_flux_shell)
                total_flux_shel.resize( times.size(), 0.0 );
            for (size_t it = 0; it < times.size(); ++it){
//                auto images = computeEjectaSkyMapPieceWise( times[it],freqs[ifreq]);
//                std::vector<Image> images(nshells_);
                computeEjectaSkyMapPieceWise( images, times[it],freqs[ifreq]);
                for (size_t i_vn = 0; i_vn < IMG::m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells_; ish++) {
                        out[ii][i_vn][ish] = images.getImgRef(ish).m_data[i_vn];//arrToVec(images[ish].m_data[i_vn]);
                    }
                }
                for (size_t ish = 0; ish < nshells_; ish++) {
                    tota_flux[it] += images.getImgRef(ish).m_f_tot;
                    total_flux_shell[ish][it] = images.getImgRef(ish).m_f_tot;
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

        auto in_group_names = IMG::m_names;//dummy.m_names;

        /// add attributes from model parameters
        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells_},
                {"nshells", nlayers_}
        };

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
//        size_t m_nlayers = p_cumShells->m_nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        if (!is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot evaluateShycnhrotronSpectrum light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            std::cerr << " ejecta observer parameters are not set. Cannot evaluateShycnhrotronSpectrum light curve (analytic) exiting...\n";
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

#endif //SRC_MODEL_EJECTA_H
