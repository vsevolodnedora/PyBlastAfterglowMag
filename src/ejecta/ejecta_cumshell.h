//
// Created by vsevolod on 27/07/23.
//

#ifndef SRC_EJECTA_CUMSHELL_H
#define SRC_EJECTA_CUMSHELL_H

//#include "../utilitites/pch.h"
//#include "../utilitites/utils.h"
//#include "../utilitites/interpolators.h"
//#include "../utilitites/ode_solvers.h"
//#include "../utilitites/quadratures.h"
//#include "../utilitites/rootfinders.h"
//#include "../image.h"
//#include "../synchrotron_an.h"

//#include "model_magnetar.h"
//#include "blastwave/blastwave_components.h"
//#include "blastwave/blastwave.h"
//#include "blastwave/blastwave_collision.h"

/// Radially structured blastwave collection
class CumulativeShell{
    struct Pars{
        size_t ilayer = 0;
        size_t nshells = 0;
        size_t n_active_shells = 0;
//        size_t n_active_shells_im1 = 0;
        bool first_entry_as_single_shell = true;
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
    CumulativeShell(Vector t_grid, size_t nshells, int ilayer, size_t n_substeps, BW_TYPES type, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CumulativeShell");
        p_coll = std::make_unique<BlastWaveCollision>(loglevel);
        p_pars = std::make_unique<Pars>();
        p_pars->ilayer=ilayer;
        p_pars->nshells=nshells;
        p_pars->n_active_shells=nshells;
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws.emplace_back(std::make_unique<BlastWave>(t_grid, ishell, ilayer, n_substeps, type, loglevel ) );
        if (t_grid.empty())
            return;
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize(nshells);
//        for (auto & arr : m_data_prev)
//            arr.resizeEachImage(nshells);
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
//                bw1->getData(static_cast<BW::Q>(key))[it] = mD[key][i];

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
    inline Vector & getEthVec(){return m_data[Q::ieth]; }
    inline Vector & getVolVec(){return m_data[Q::ivol]; }
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
    inline size_t nBWs() const {
        if (p_bws.empty() || p_bws.size() < 1){
            (*p_log)(LOG_ERR,AT) << " bws are not initizlized\n";
            exit(1);
        }
        return p_bws.size();
    }
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
//        p_pars->n_active_shells_im1 = p_pars->n_active_shells;
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
            mD[Q::idelta][idx] = dr_i;
            mD[Q::ivol][idx] = vol_i;
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
                mD[Q::idelta][idx] = dr_i;
                mD[Q::ivol][idx] = vol_i;
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
            mD[Q::idelta][idx] = dr_i;
            mD[Q::ivol][idx] = vol_i;
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
        mD[Q::idelta][0] = dr_i;
        mD[Q::ivol][0] = vol_i;
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

        /// --------------------------------------------
        double _r_i,_vol_i;
        if (p_pars->n_active_shells>1){
            /// copy from previos shell
            size_t ii = p_pars->n_active_shells-1;
            size_t idx = m_idxs[ii];
            p_bws[idx]->getPars()->delta = 1.01 * p_bws[idx-1]->getPars()->delta;
            p_bws[idx]->getPars()->vol = 1.01 * p_bws[idx-1]->getPars()->vol;
            m_data[Q::idelta][ii] = p_bws[idx]->getPars()->delta;//dr_i;
            m_data[Q::ivol][ii] = p_bws[idx]->getPars()->vol;
            _r_i = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR];
            _vol_i = (4./3.) * CGS::pi * (_r_i*_r_i*_r_i);
            p_bws[idx]->getPars()->_last_frac = _r_i / m_data[Q::idelta][ii];//maxValue(mD[Q::idelta]);//mD[Q::idelta][ii];
            p_bws[idx]->getPars()->_last_frac_vol = _vol_i / m_data[Q::ivol][ii];//maxValue(mD[Q::ivol]); // mD[Q::ivol][ii];
        }
            /// if only one shell exists
        else {
            size_t idx = m_idxs[0];
            auto & p_bwpars = p_bws[idx]->getPars();
            auto & p_bw = p_bws[idx];
            auto & bw_data = p_bw->getData();

            // use same fraction as at the last time there were many shells
            size_t ii = p_pars->n_active_shells-1;
//            size_t idx = m_idxs[ii];
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
//            mD[Q::idelta][ii] = p_bws[idx]->getPars()->delta;//dr_i;
//            mD[Q::ivol][ii] = p_bws[idx]->getPars()->vol;

            switch (p_bwpars->method_single_bw_delta) {
                /// keep thickness constant
                case iconst:
                    break;
                    /// use the radius and apply a fraction to it from last time there were > 1 BW
                case ifrac_last:
                    if ((p_pars->n_active_shells == 1) and (p_pars->first_entry_as_single_shell)) {
                        (*p_log)(LOG_INFO, AT) << "[il=" << p_pars->ilayer << ", ish=" << idx << "] "
                                                                                                 "Only one shell left. Using MAX delta & volume...\n";

                        auto iter = std::max_element(bw_data[BW::Q::iEJdelta].rbegin(), bw_data[BW::Q::iEJdelta].rend()).base();
                        size_t idx_max = std::distance(bw_data[BW::Q::iEJdelta].begin(), std::prev(iter));
                        double max_delta = bw_data[BW::Q::iEJdelta][idx_max];
                        double r_at_max = bw_data[BW::Q::iR][idx_max];

                        p_bws[idx]->getPars()->_last_frac = r_at_max / max_delta;//maxValue(mD[Q::idelta]);//mD[Q::idelta][ii];

                        auto iter_ = std::max_element(bw_data[BW::Q::iEJvol].rbegin(), bw_data[BW::Q::iEJvol].rend()).base();
                        size_t idx_max_ = std::distance(bw_data[BW::Q::iEJvol].begin(), std::prev(iter));
                        double max_vol_ = bw_data[BW::Q::iEJvol][idx_max];
                        double r_at_max_ = bw_data[BW::Q::iR][idx_max];
                        _vol_i = (4./3.) * CGS::pi * (r_at_max*r_at_max*r_at_max);
                        p_bws[idx]->getPars()->_last_frac_vol = _vol_i / max_vol_;//maxValue(mD[Q::ivol]); // mD[Q::ivol][ii];

                        p_pars->first_entry_as_single_shell = false;
                    }

                    _r_i = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR];
                    _vol_i = (4./3.) * CGS::pi * (_r_i*_r_i*_r_i);
                    p_bws[idx]->getPars()->delta = _r_i / p_bws[idx]->getPars()->_last_frac;
                    p_bws[idx]->getPars()->vol = _vol_i / p_bws[idx]->getPars()->_last_frac_vol;
                    break;
                case ibm:
                    (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                    exit(1);
//                    rho_bm = getBM()->rho_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2) * CGS::mp;
//                    gam_cbm_bm = getBM()->gamma_downstream(ej_R, j_R, j_Gamma, j_rho / CGS::mp, j_rho2);
                    break;
                case ist:
                    (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                    exit(1);
                    break;
                case ibm_st:
                    (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                    exit(1);
                    break;
                case ilr:
                    if ((p_pars->n_active_shells == 1) and (p_pars->first_entry_as_single_shell)){
                        (*p_log)(LOG_INFO,AT) << "[il="<<p_pars->ilayer<<", ish="<<idx <<"] "
                                                                                         "Only one shell left. Interpolating delta & volume histories...\n";
                        if (p_bw->getLRforDelta()->isTrained()){
                            (*p_log)(LOG_ERR,AT) << " Parameters are already set for regression... (they should be set now, npt before...)!\n";
                            exit(1);
                        }

                        p_bw->getLRforDelta()->commputeInput();
                        p_bw->getLRforDelta()->calculateCoefficient();
                        p_bw->getLRforDelta()->calculateConstantTerm();

                        p_bw->getLRforVol()->commputeInput();
                        p_bw->getLRforVol()->calculateCoefficient();
                        p_bw->getLRforVol()->calculateConstantTerm();

                        std::cout << " ----------------------------------------------------------------------------- \n";
                        (*p_log)(LOG_INFO,AT) << " Result: Delta coeffs=(w="<<p_bw->getLRforDelta()->getCoeff()
                                              <<", b="<<p_bw->getLRforDelta()->getConst()<<"\n";
                        print_xy_as_numpy(p_bw->getData()[BW::Q::iR],p_bw->getData()[BW::Q::iEJdelta],
                                          p_bw->getData()[BW::Q::iR].size(),100);
                        std::cout << " ----------------------------------------------------------------------------- \n";
                        std::cout << p_bw->getData()[BW::Q::iR]<<"\n";
                        (*p_log)(LOG_INFO,AT) << " Result: Volume coeffs=(w="<<p_bw->getLRforDelta()->getCoeff()
                                              <<", b="<<p_bw->getLRforDelta()->getConst()<<"\n";
                        print_xy_as_numpy(p_bw->getData()[BW::Q::iR],p_bw->getData()[BW::Q::iEJvol],p_bw->getData()[BW::Q::iR].size(),100);
                        std::cout << " ----------------------------------------------------------------------------- \n";
                        p_pars->first_entry_as_single_shell = false;
                    }
                    _r_i = Y[p_bws[idx]->getPars()->ii_eq + SOL::QS::iR];
                    p_bws[idx]->getPars()->delta = p_bw->getLRforDelta()->predict(_r_i);
                    p_bws[idx]->getPars()->vol = p_bw->getLRforVol()->predict(_r_i);
                    break;
            }
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
//            mD[Q::idelta][ii] = dr_i;
//            mD[Q::ivol][ii] = vol_i;
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
//            double m_beta = EQS::BetFromMom( Y[bw->getPars()->ii_eq + SOL::QS::imom] );
            double m_beta = Beta( Y[bw->getPars()->ii_eq + SOL::QS::iGamma] );
            if (m_beta > 1){
                (*p_log)(LOG_ERR,AT) << " m_beta="<<m_beta<<"\n";
                exit(1);
            }
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
//            bw->getPars()->tau_to0 = mD[Q::itaucum][ii];
//            bw->getPars()->dtau = mD[Q::idtau][ii];
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
        if (!std::isfinite(mtot) || (mtot <= 0)){
            (*p_log)(LOG_ERR,AT) << "mtot is nan or < 0; mtot="<<mtot<<"\n";
            exit(1);
        }
        return mtot;
    }
    double getShellVolume(const double * Y){
        double total_volume = 0.;
        double mtot = getShellMass(Y);
        double cummass = 0.;
        for (size_t i = 0; i < p_pars->n_active_shells; i++) {

            /// try to limit the volume to the dense part of the ejecta...
            size_t i_idx = m_idxs[i];
            double m2 = Y[p_bws[i_idx]->getPars()->ii_eq + SOL::QS::iM2];
            double m0 = p_bws[i_idx]->getPars()->M0;
            double m2plus0 = (1. + m2) * m0;

            cummass += m2plus0;

            if ((cummass >= .9 * mtot)&&(p_pars->n_active_shells > 1)&&(total_volume>0))
                break;

            total_volume += m_data[Q::ivol][i];
        }
//        double r0 = Y[p_bws[ m_idxs[0] ]->getPars()->ii_eq + SOL::QS::iR];
//        double r1 = Y[p_bws[ m_idxs[p_pars->n_active_shells - 1] ]->getPars()->ii_eq + SOL::QS::iR];
//        if ((r0 >= r1)||(r0==0)||(r1==0)){
//            (*p_log)(LOG_ERR,AT)<<" r0 > r1. in the shell; r0="<<r0<<" r1="<<r1<<"\n";
//            exit(1);
//        }
//        double delta = r1-r0;
//        double total_volume = (4./3.) * CGS::pi * (r1*r1*r1 - r0*r0*r0) / p_bws[ m_idxs[0] ]->getPars()->ncells;
        if (!std::isfinite(total_volume) || (total_volume <= 0)){
            (*p_log)(LOG_ERR,AT) << "volume is nan or < 0; volume="<<total_volume<<"\n";
            exit(1);
        }
        return total_volume;
    }
    double getShellRho(const double * Y){
        double mtot = getShellMass(Y);
        double volume = getShellVolume(Y);
        return (mtot / volume);
    }
    double getShellOptDepth(){
        if (!std::isfinite(p_pars->tautot) || (p_pars->tautot <= 0)){
            (*p_log)(LOG_ERR,AT) << "p_pars->tautot is nan or < 0; p_pars->tautot="<<p_pars->tautot<<"\n";
            exit(1);
        }
        return p_pars->tautot;
    }
};

#endif //SRC_EJECTA_CUMSHELL_H
