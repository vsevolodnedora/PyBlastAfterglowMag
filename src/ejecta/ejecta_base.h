//
// Created by vsevolod on 28/07/23.
//

#ifndef SRC_EJECTA_BASE_H
#define SRC_EJECTA_BASE_H

//#include "../utilitites/pch.h"
//#include "../utilitites/utils.h"
//#include "../utilitites/interpolators.h"
//#include "../utilitites/ode_solvers.h"
//#include "../utilitites/quadratures.h"
//#include "../utilitites/rootfinders.h"
//#include "../image.h"
//#include "../synchrotron_an.h"


//#include "../blastwave/blastwave_components.h"
//#include "../blastwave/blastwave.h"
//#include "../blastwave/blastwave_collision.h"
#include "ejecta_cumshell.h"

struct ImageExtend{
    double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min(), ymax = std::numeric_limits<double>::min();
};

/// Radially/Angular structured Blastwave collection
class EjectaBase{
    std::unique_ptr<logger> p_log = nullptr;
public:
    std::unique_ptr<EjectaID2> id = nullptr;
    std::vector<std::unique_ptr<CumulativeShell>> p_cumShells {};
    int m_loglevel{};
//    LatStruct::METHOD_eats ejecta_eats_method{};
    Vector & t_arr;
    /// -------------------------------------
    std::vector<std::vector<char>> status{};
    Vector tobs_ej_max{};
    /// -------------------------------------

    bool do_eninj_inside_rhs = false;
    bool run_bws=false, save_dyn=false, load_dyn=false, do_ele=false, do_spec=false, save_spec=false, do_lc=false, do_skymap=false;
    bool do_collision = false;
    bool do_nuc = false;
    bool is_ejecta_obs_pars_set = false;
    bool is_ejecta_anal_ele_computed = false;
    bool is_ejecta_anal_synch_computed = false;
    double im_max_theta;
    StrDbMap m_pars; StrStrMap m_opts;
    std::string working_dir{}; std::string parfilename{};
    EjectaBase(Vector & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Ejecta");
//        p_out = std::make_unique<Output>(loglevel);
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

    /// Parameters
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
            run_bws   = getBoolOpt("run_bws", m_opts, AT, p_log, false, true);
            save_dyn  = getBoolOpt("save_dynamics", m_opts, AT, p_log, false, true);
            load_dyn  = getBoolOpt("load_dynamics", m_opts, AT, p_log, false, true);
            do_ele    = getBoolOpt("do_ele", m_opts, AT, p_log, false, true);
            do_spec   = getBoolOpt("do_spec", m_opts, AT, p_log, false, true);
            save_spec = getBoolOpt("save_spec", m_opts, AT, p_log, false, true);
            do_lc     = getBoolOpt("do_lc", m_opts, AT, p_log, false, true);
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
//                bool use_1d_id = getBoolOpt("use_1d_id", m_opts, AT, p_log, false, true);
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
//                std::cerr << "R0=" << id->get(0,0,EjectaID2::Q::ir) << "\n";
                setEjectaBwPars(m_pars, m_opts, ii_eq, iljet);
            }
            /// --- For optical depth verbosity ....
            status.resize(nlayers());
            for (auto & arr : status)
                arr.resize(nshells(), '-');
            tobs_ej_max.resize(nlayers(),0.);
            im_max_theta = getDoublePar("im_max_theta", m_pars, AT, p_log, CGS::pi/2., true);
            if (im_max_theta <= 0.)
                im_max_theta = id->theta_wing;
        }
        else{
            (*p_log)(LOG_INFO, AT) << "ejecta is not initialized and will not be considered.\n";
        }
    }

    void setEjectaBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq, size_t n_layers_jet){

        run_bws = getBoolOpt("run_bws", opts, AT, p_log, false, true);
        load_dyn = getBoolOpt("load_dynamics", opts, AT, p_log, false, true);
        size_t n_substeps = (size_t)getDoublePar("n_store_substeps",pars,AT,p_log,10, true);

        if ((!run_bws) && (!load_dyn))
            return;

        if ((!run_bws) && (load_dyn)){
            is_ejecta_obs_pars_set = true;
            for(size_t il = 0; il < nlayers(); il++) {
                p_cumShells.push_back(
                        std::make_unique<CumulativeShell>(Vector {}, nshells(), il, n_substeps,
                                                          p_log->getLogLevel()));
                p_cumShells[il]->setPars(pars, opts);
                for (size_t ish = 0; ish < nshells(); ish++){
                    auto & bw = p_cumShells[il]->getBW(ish);
                    bw->setParams(id, pars, opts, il, ii_eq);
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
//        size_t n_ejecta_empty_images = 0;
//        std::vector<std::vector<size_t>> n_empty_images;
//        std::vector<size_t> n_empty_images_shells;
        size_t nshells_ = nshells();//ejectaStructs.nshells;
        size_t n_layers_ej_ = nlayers();//ejectaStructs.structs[0].m_nlayers;
        if (n_layers_ej_ == 0){
            (*p_log)(LOG_ERR,AT)<<" no layers found to evolve!\n";
            exit(1);
        }
//        std::vector<std::vector<size_t>> n_empty_images_layer_shell;
//        for (auto & n_empty_images_layer : n_empty_images_layer_shell)
//            n_empty_images_layer.resizeEachImage(nshells_);
        /// include blastwave collision between velocioty shells into the run
        std::vector<std::string> empty_bws{};
        size_t n_unitinitilized_shells=0;
        for(size_t il = 0; il < n_layers_ej_; il++){
            p_cumShells.push_back(
                    std::make_unique<CumulativeShell>(t_arr, nshells_, il, n_substeps,
                                                      p_log->getLogLevel()) );
            p_cumShells[il]->setPars(pars, opts);
            empty_bws.emplace_back("il="+std::to_string(il)+" | shells ");
            for (size_t ish = 0; ish < nshells_; ish++){
                auto & bw = p_cumShells[il]->getBW(ish);

                if (id->get(ish,il,EjectaID2::Q::ir)<=0||
                    id->get(ish,il,EjectaID2::Q::iek)<=0||
                    id->get(ish,il,EjectaID2::Q::imass)<=0){
                    empty_bws[il]+=" "+std::to_string(ish) + " ENDING evolution";
                    n_unitinitilized_shells++;
                    bw->getPars()->end_evolution = true;
                    ii_eq += SOL::neq;//bw->getNeq();
                    continue;
                }
                else{
                    bw->setParams(id, pars, opts, il, ii_eq);
                    ii_eq += SOL::neq;//bw->getNeq();
                }

//                auto & struc = ejectaStructs.structs[ish];

#if 0
                switch (id->method_eats) {
                    case EjectaID2::iadaptive:
                        bw->getRad()->setEatsPars(
                                pars,opts,id->nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
                                0.,0.,id->theta_wing,
                                getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                        break;
                    case EjectaID2::ipiecewise:
                        bw->getEATS()->setEatsPars(
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
//                else if(n_layers_ej_ == 0){
////                    n_ejecta_empty_images += 1;
////                    n_empty_images_layer_shell[ish].emplace_back(il);
////                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
//                }
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


//                bw->getEATS()->setEatsPars(pars, opts);
//                bw->getSynchAnPtr()->setPars( pars, opts );

//                ii++;
            }
        }

//        is_ejBW_init = true;
        is_ejecta_obs_pars_set = true;

        if ((n_unitinitilized_shells > 0)){
            (*p_log)(LOG_WARN,AT)<<"-------------- NO ID FOR ------------"<<"\n";
            for (const auto & empty_bw : empty_bws){
                (*p_log)(LOG_WARN,AT)<<empty_bw<<"\n";
            }
            (*p_log)(LOG_WARN,AT)<<"-------------------------------------"<<"\n";
        }
        if (n_unitinitilized_shells == nshells_*n_layers_ej_){
            (*p_log)(LOG_ERR,AT)<<" No BW from ID file were intilazied."<<"\n";
            exit(1);
        }

//        if ((p_log->getLogLevel() > LOG_WARN)) {
//            if (n_ejecta_empty_images > 0) {
//                auto &ccerr = std::cout;
//                ccerr << "Ejecta blastwave is NOT initialized for total n="
//                      << n_ejecta_empty_images << " layers. Specifically:\n";
//                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
////                    auto &ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                    size_t n_layers_i = nlayers();//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
//                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [";
//                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
//                        ccerr << n_empty_images[ish][il] << " ";
//                    }
//                    ccerr << "] / (" << n_layers_i << " total layers) \n";
//                }
//            }
//        }

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
//                double MomIm1 = Ym1[bws->getPars()->ii_eq + SOL::QS::imom];
                double GamIm1 = Ym1[bws->getPars()->ii_eq + SOL::QS::iGamma];
//                double MomI = Y[bws->getPars()->ii_eq + SOL::QS::imom];
                double GamI = Y[bws->getPars()->ii_eq + SOL::QS::iGamma];
                double Eint2I = Y[bws->getPars()->ii_eq + SOL::QS::iEint2];
                double Mom0 = bws->getPars()->mom0;
                double Gam0 = bws->getPars()->Gamma0;
//                if (MomI > MomIm1){
//                    /// acceleration
//                    n_accel += 1;
//                }
//                else if (MomI < MomIm1){
//                    /// deceleration
//                    n_decel += 1;
//                }
//                /// find fastest
//                if (MomI/Mom0 > Mom_max_over_Gamma0){
//                    Mom_max_over_Gamma0 = MomI/Mom0;
//                    il_wich_fastest = (int)il;
//                    ish_with_fastest = (int)ish;
//                }
//                /// find slowest
//                if (MomI/Mom0 < Mom_min_over_Gamma0){
//                    Mom_min_over_Gamma0 = MomI/Mom0;
//                    il_with_slowest = (int)il;
//                    ish_with_slowest = (int)ish;
//                }
                if (GamI > GamIm1){
                    /// acceleration
                    n_accel += 1;
                }
                else if (GamI < GamIm1){
                    /// deceleration
                    n_decel += 1;
                }
                /// find fastest
                if (GamI/Gam0 > Mom_max_over_Gamma0){
                    Mom_max_over_Gamma0 = GamI/Gam0;
                    il_wich_fastest = (int)il;
                    ish_with_fastest = (int)ish;
                }
                /// find slowest
                if (GamI/Gam0 < Mom_min_over_Gamma0){
                    Mom_min_over_Gamma0 = GamI/Gam0;
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

    /// Compute Skymap for piece-wise EATS
    ImageExtend computeEjectaSkyMapPW_new(std::vector<VecVector> & out, Vector & fluxes, double obs_time, double obs_freq ){
        /// out is [i_vn][ish][itheta_iphi]
        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();
        if (out.empty()){
            (*p_log)(LOG_ERR,AT) << " isEmpty image passed. Exiting...\n";
            exit(1);
        }
        if (out[0].size() != nshells_){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }
        /// Allocate memory for ctheta and cphi grids ('2' is for principle and counter jet)
        for (size_t ivn = 0; ivn < IMG::m_names.size(); ivn++)
            for (size_t ish = 0; ish < nshells_; ish++)
                out[ivn][ish].resize(2 * ncells_, 0. );


        ImageExtend im{};
        size_t n_jet_empty_images = 0;
        for (size_t ishell = 0; ishell < nshells_; ishell++){
            std::vector<size_t> n_empty_images_layer;
            size_t offset = 0; // offset the data from this layer in the container
            for (size_t ilayer = 0; ilayer < nlayers_; ilayer++){
                /// Evaluate a given image
                size_t cil = EjectaID2::CellsInLayer(ilayer);
                auto & bw_rad = p_cumShells[ilayer]->getBW(ishell)->getFsEATS();
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA SkyMap obs_time="<<obs_time<<" obs_freq="<<obs_freq
                        << " vel_shell="<<ishell<<"/"<<nshells()-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers()
                        << " phi_cells="<<cil<<"\n";
                double flux = bw_rad->evalSkyMapPW(out, obs_time, obs_freq, offset);
                fluxes[ishell] += flux;
                offset += cil;
            }

            /// find the image extend for this shell
            double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
            double imin = std::numeric_limits<double>::max(), imax = std::numeric_limits<double>::min();
            for (size_t i = 0; i < 2 * ncells_; i++){
                if (out[IMG::Q::ixr][ishell][i] < xmin) xmin = out[IMG::Q::ixr][ishell][i];
                if (out[IMG::Q::ixr][ishell][i] > xmax) xmax = out[IMG::Q::ixr][ishell][i];
                if (out[IMG::Q::iyr][ishell][i] < ymin) ymin = out[IMG::Q::iyr][ishell][i];
                if (out[IMG::Q::iyr][ishell][i] > ymax) ymax = out[IMG::Q::iyr][ishell][i];
                if (out[IMG::Q::iintens][ishell][i] < imin) imin = out[IMG::Q::iintens][ishell][i];
                if (out[IMG::Q::iintens][ishell][i] > imax) imax = out[IMG::Q::iintens][ishell][i];
            }

            /// normalize image to be mJy/mas^2
            double delta_x = xmax - xmin;
            double delta_y = ymax - ymin;
            double dfnu = fluxes[ishell] / (delta_x * delta_y);
            for (size_t i = 0; i < 2 * ncells_; i++){
                out[IMG::Q::iintens][ishell][i] *= (dfnu / imax);
            }

            /// update image boundaries for overall image
            if (xmin < im.xmin) im.xmin = xmin;
            if (ymin < im.ymin) im.ymin = ymin;
            if (xmax > im.xmax) im.xmax = xmax;
            if (ymax > im.ymax) im.ymax = ymax;
        }
        return im;
    }

    /// Compute Skymap for adaptive EATS
    ImageExtend computeEjectaSkyMapA_new(std::vector<VecVector> & out, Vector & fluxes, double obs_time, double obs_freq, size_t nsublayers, size_t nphi ){

        /// out is [i_vn][ish][itheta_iphi]

        size_t nlayers_ = nlayers();
        // ------------------------
        if (out.empty()){
            (*p_log)(LOG_ERR,AT) << " isEmpty image passed. Exiting...\n";
            exit(1);
        }

        double rtol = 1.e-6;

        /// Allocate memory for ctheta and cphi grids ('2' is for principle and counter jet)
        size_t ncells = 0;
        for (size_t ill = 0; ill < nsublayers; ++ill){
            ncells += EjectaID2::CellsInLayer_(ill);
        }
        (*p_log)(LOG_INFO,AT)
                << " With nlayers=" << nlayers_ << " nsublayers=" << nsublayers << " ncells=" << ncells << "\n";

        /// resize or empty the contaner for temporary solution
        for (size_t ivn = 0; ivn < IMG::m_names.size(); ivn++) {
            for (size_t ish = 0; ish < nlayers_; ish++) {
//                if (out[ivn][ish].size() != 2 * ncells)
                out[ivn][ish].resize(2 * ncells, 0.);
//                else{
//                    std::fill(out[ivn][ish].begin(), out[ivn][ish].end(), 0.0);
//                }
            }
        }


        /// locate the extend of the jet at a given time (for further grid creation)
        double th_l_prev = 0.;
        double th_jet_min = 0., th_jet_max = 0.;
        Vector th_b_layers (nlayers_);
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer){
            double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;
            if (ilayer > 0)
                th_l_prev = p_cumShells[ilayer-1]->getBW(0)->getPars()->theta_c_l;
            auto & bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
            double th_b_min, th_b_max;
            bw_rad->findJetEdgeA(th_b_min, th_b_max, obs_time, obs_freq, th_l_prev);
            th_b_layers[ilayer] = th_b_max;
            if (th_b_max > th_jet_max)
                th_jet_max = th_b_max;
        }
        (*p_log)(LOG_INFO,AT)
                << " Jet Edge Found theta="<<th_jet_max<<" at obs_time="<<obs_time<<" obs_freq="<<obs_freq << "\n";


        /// compute flux density for each layer ( for normalization )
        double tot_flux = 0.;
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer) {
            auto &bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
//            double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;
            double atol = tot_flux * rtol / (double) nlayers_;
            double layer_flux = bw_rad->evalFluxDensA(obs_time, obs_freq, atol);
            tot_flux += layer_flux;
            fluxes[ilayer] = layer_flux;
        }
        /// make an interpolator for sublayers (for a given ctheta -> get flux)
//        auto interp = Interp1d(id->getVec(0,EjectaID2::Q::itheta_c_l), fluxes);

        Vector theta_pw (nsublayers + 1 );
        for (size_t i = 0; i < nsublayers + 1; i++){
            double fac = (double) i / (double)nsublayers;
            theta_pw[i] = 2.0 * std::asin( fac * std::sin(th_jet_max / 2.0 ) );
        }
        Vector cthetas (nsublayers );
        for (size_t i = 0; i < nsublayers; i++)
            cthetas[i] = theta_pw[i] + (theta_pw[i+1] - theta_pw[i]) * 0.5;

        ImageExtend im{};

        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer){

            auto & bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();

            size_t iii = 0;
            double layer_summed_intensity = 0;
            /// loop for all thetas from 0 to overall theta_edge # TODO see if you can just do it until current theta_edge



            for (size_t ith = 0; ith < nsublayers; ith++){
                size_t cil = EjectaID2::CellsInLayer_(ith);

                /// fill the angular grid for intensity evaluations
                for (size_t iphi = 0;  iphi < cil; iphi++){
                    double cphi = (double)iphi * 2.0 * CGS::pi / (double)cil;
                    /// principle jet
                    out[IMG::Q::icphi][ilayer][iii+iphi] = cphi;
                    out[IMG::Q::ictheta][ilayer][iii+iphi] = cthetas[ith];
                    /// counter jet
                    out[IMG::Q::icphi][ilayer][ncells+iii+iphi] = cphi;
                    out[IMG::Q::ictheta][ilayer][ncells+iii+iphi] = cthetas[ith];
                }

                /// evaluate intensity
                layer_summed_intensity += bw_rad->evalSkyMapA(out, obs_time, obs_freq, ilayer, iii, cil, ncells);
                iii = iii + cil;
            }

            /// normalize the skymap
            double total_int = 0.;
            double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
            double max_int = std::numeric_limits<double>::min();
            for (size_t i = 0; i < 2 * ncells; i++){
                if (out[IMG::Q::ixr][ilayer][i] < xmin) xmin = out[IMG::Q::ixr][ilayer][i];
                if (out[IMG::Q::ixr][ilayer][i] > xmax) xmax = out[IMG::Q::ixr][ilayer][i];
                if (out[IMG::Q::iyr][ilayer][i] < ymin) ymin = out[IMG::Q::iyr][ilayer][i];
                if (out[IMG::Q::iyr][ilayer][i] > ymax) ymax = out[IMG::Q::iyr][ilayer][i];

                total_int += out[IMG::Q::iintens][ilayer][i];
                if (max_int < out[IMG::Q::iintens][ilayer][i])
                    max_int = out[IMG::Q::iintens][ilayer][i];
            }

            if (xmax == xmin){
                (*p_log)(LOG_ERR,AT)<<" xmin=xmax="<<xmin<<"\n";
                exit(1);
            }
            if (ymax == ymin){
                (*p_log)(LOG_ERR,AT)<<" ymin=ymax="<<ymin<<"\n";
                exit(1);
            }
            double delta_x = xmax - xmin;
            double delta_y = ymax - ymin;

            double d_l = p_cumShells[ilayer]->getBW(0)->getPars()->d_l;
            for (size_t i = 0; i < 2 * ncells; i++) {
                out[IMG::Q::iintens][ilayer][i] = out[IMG::Q::iintens][ilayer][i] / max_int * fluxes[ilayer] / delta_x / delta_y;// / (delta_x * CGS::rad2mas * delta_y * CGS::rad2mas) * d_l * d_l;
            }

            /// update image boundaries for overall image
//            if (xmin < im.xmin) im.xmin = xmin;
//            if (ymin < im.ymin) im.ymin = ymin;
//            if (xmax > im.xmax) im.xmax = xmax;
//            if (ymax > im.ymax) im.ymax = ymax;
        }
        return im;
    }

    /// Compute lightcurve for piece-wise EATS
    void evalEjectaLightCurves_new(VecVector & out, Vector & obs_times, Vector & obs_freqs) {
        /// out // [ish*il][it*inu]
        (*p_log)(LOG_INFO, AT) << " starting ejecta light curve calculation\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++) {
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                (*p_log)(LOG_INFO, AT)
                        << " EJECTA LC ntimes=" << obs_times.size()
                        << " vel_shell=" << ishell << "/" << nshells() - 1
                        << " theta_layer=" << ilayer << "/" << nlayers()
                        << "\n";
                auto &model = getShells()[ilayer];//ejectaModels[ishell][ilayer];
                model->getBW(ishell)->getFsEATS()->evalLightCurve(out[ii], id->method_eats, obs_times, obs_freqs);
                ii++;
            }
        }
    }
};


#endif //SRC_EJECTA_BASE_H
