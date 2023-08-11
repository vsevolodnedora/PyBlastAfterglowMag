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

/// Radially/Angular structured Blastwave collection
class EjectaBase{
//    VelocityAngularStruct ejectaStructs{};
//    std::unique_ptr<Output> p_out = nullptr;
    std::unique_ptr<logger> p_log = nullptr;
public:
    std::unique_ptr<EjectaID2> id = nullptr;
    std::vector<std::unique_ptr<CumulativeShell>> p_cumShells {};
//    bool is_ejBW_init = false;
//    bool is_ejecta_obsrad_pars_set = false;
//    bool is_ejecta_struct_set = false;
//    double jet_layer_fnu_stop_frac=1e-5;
//    int n_ode_eq{};
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
    void computeEjectaSkyMapPW(Images & images, double obs_time, double obs_freq ){

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();

        if (images.isEmpty()){
            (*p_log)(LOG_ERR,AT) << " isEmpty image passed. Exiting...\n";
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
        Images tmpImages (nlayers_, IMG::m_names.size());
        tmpImages.resizeEachImage(ncells_);
//        for (auto & _tmp : tmpImages)
//            _tmp.resizeEachImage( ncells_ );
        Image tmp_pj( ncells_, IMG::m_names.size(), 0, m_loglevel);
        Image tmp_cj( ncells_, IMG::m_names.size(), 0, m_loglevel);
        for (size_t ishell = 0; ishell < nshells_; ishell++){
//            for (auto & _tmp : tmpImages)
//                _tmp.clearData();
            tmpImages.clearEachImage();
            std::vector<size_t> n_empty_images_layer;
            double atol=0; // TODO make it depend on the layer flux density
            Image & layerImage = images.getReferenceToTheImage(ishell);
            for (size_t ilayer = 0; ilayer < nlayers_; ilayer++){
                tmp_pj.clearData(); tmp_cj.clearData();
                /// Evaluate a given image --------------------------------------
                auto & bw_rad = p_cumShells[ilayer]->getBW(ishell)->getFsEATS();
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC obs_time="<<obs_time<<" obs_freq="<<obs_freq
                        << " vel_shell="<<ishell<<"/"<<nshells()-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers()
                        << " phi_cells="<<EjectaID2::CellsInLayer(ilayer)<<"\n";
                bw_rad->evalImagePW(tmpImages.getReferenceToTheImage(ilayer), tmp_pj, tmp_cj, obs_time, obs_freq);

//                bw_rad->evalImageA(tmpImages.getReferenceToTheImage(ilayer), tmp_pj, tmp_cj, obs_time, obs_freq);
                /// -------------------------------------------------------------
                if (tmpImages.getReferenceToTheImage(ilayer).m_f_tot == 0){
                    n_jet_empty_images += 1;
                    n_empty_images_layer.emplace_back(ilayer);
                }
            }
            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
            combineImages(layerImage, ncells_, nlayers_, tmpImages) ;

//            Vector times {obs_time}; Vector freqs {obs_freq};
//            auto lc = evalEjectaLightCurves(times, freqs);
//            double totflux = 0;
//            for (size_t il = 0; il < nlayers_; il++)
//                totflux += lc[0][il][0];
//            std::cout <<" Lc: "<< totflux << " Im: "<<images.getReferenceToTheImage(0).m_f_tot<<" t="<<obs_time<<" nu="<<obs_freq<<"\n";
        }

        /// print which layers/shells gave isEmpty image
        if (p_log->getLogLevel() == LOG_INFO) {
            if (n_jet_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an isEmpty images for total n="
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

    /// Compute Skymap for adaptive EATS
    void computeEjectaSkyMapA_old(Images & images, double obs_time, double obs_freq, size_t nsublayers ){

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
//        size_t ncells_ =  (int)ncells();

//        size_t nsublayers = ;

        /// for skymap we need an extended grid of layers so that there is enough resolution for the image
        Vector theta_c{};
        Vector theta_c_l{};
        Vector theta_c_h{};
        std::vector<size_t> cils{};
//        EjectaID2::_init_a_grid(theta_c_l, theta_c_h, theta_c, nlayers_ * nsublayers, CGS::pi/2.);
        EjectaID2::_init_a_grid(theta_c_l, theta_c_h, theta_c, nlayers_ * nsublayers, im_max_theta);//id->theta_wing);
        EjectaID2::_evalCellsInLayer(nlayers_ * nsublayers, cils);
        size_t ncells = EjectaID2::_evalTotalNcells(nlayers_ * nsublayers);
        if (ncells > 1e6){
            (*p_log)(LOG_WARN,AT) << " for adaptive EATS image calc. large ncells="<<ncells<<"\n";
        }




//        size_t nphi = 100;
//        ncells = nlayers_ * nsublayers * nphi;


//        size_t ntheta = nlayers_*nsublayers;
//        Vector cthetas0;
//        Vector thetas ( ntheta + 1 );
//        cthetas0.resize( ntheta );
//        for (size_t i = 0; i < ntheta + 1; i++){
//            double fac = (double)i / (double)ntheta;
//            thetas[i] = 2.0 * asin( fac * sin(M_PI / 2.0 ) );
//        }
//        for (size_t i = 0; i < ntheta; ++i){
//            cthetas0[i] = 0.5 * ( thetas[i+1] + thetas[i] );
//        }
//        std::vector<size_t> cil;
//        cil.resize(ntheta );
//        for (size_t i = 0; i < ntheta; i++)
//            cil[i] = EjectaID2::CellsInLayer(i);
//        size_t ncells = (size_t)std::accumulate(cil.begin(), cil.end(),0); /// total number of cells
//
//        Vector all_cthetas(ncells );
//        Vector all_cphis(ncells );
//        int k = 0;
//        for (size_t i = 0; i < ntheta; i++){
//            for (size_t j = 0; j < cil[i]; j++){
//                all_cthetas[k] = cthetas0[i];
//                all_cphis[k] = (double)j * 2.0 * CGS::pi / (double)cil[i];
//                k++;
//            }
////                std::cout << all_cphis << "\ncells";
//        }


        if (images.isEmpty()){
            (*p_log)(LOG_ERR,AT) << " isEmpty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nshells_){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }

        size_t n_jet_empty_images = 0;

//        size_t ii = 0;
//        size_t tot_ncells = 0;
//        std::vector<size_t> cells_in_layers{};
//        Vector tot_theta_l{};
//        Vector tot_theta_h{};
//        for (size_t il = 0; il < nlayers_; il++){
//            double theta_l = id->get(0,il,EjectaID2::itheta_c_l);
//            double theta_h = id->get(0,il,EjectaID2::itheta_c_h);
//            tot_theta_l.push_back(theta_l);
//            tot_theta_h.push_back(theta_h);
//            double dtheta = (theta_h - theta_l) / (double)nsublayers;
//            for (size_t isubl = 0; isubl < nsublayers; isubl++){
//                tot_theta_l.push_back(theta_l + dtheta*(double)isubl);
//                tot_theta_l.push_back(theta_h + dtheta*(double)isubl);
//            }
//        }

        double rtol = 1e-6;

        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        const std::vector<std::string> x {};
        Images tmpImagesSet (nlayers_*nsublayers, IMG::m_names.size());
        tmpImagesSet.resizeEachImage(ncells*2);
//        for (auto & _tmp : tmp)
//            _tmp.resizeEachImage( ncells_ );
        Image tmp_pj( ncells, IMG::m_names.size(), 0, m_loglevel);
        Image tmp_cj( ncells, IMG::m_names.size(), 0, m_loglevel);
        for (size_t ishell = 0; ishell < nshells_; ishell++){
            tmpImagesSet.clearEachImage();
            std::vector<size_t> n_empty_images_layer ;

            size_t ii = 0;
            double tot_flux = 0.;
            for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer){
                auto & bw_rad = p_cumShells[ilayer]->getBW(ishell)->getFsEATS();

                /// eval. total flux.dens from the layer (one BW)
                double atol = tot_flux * rtol / (double)nlayers();
                double layer_flux = bw_rad->evalFluxDensA(obs_time,obs_freq, atol);
//                layer_flux /= (double)nsublayers;
                tot_flux += layer_flux;

                /// clear emages
//                tmp_pj.clearData(); tmp_cj.clearData();
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC obs_time="<<obs_time<<" obs_freq="<<obs_freq
                        << " vel_shell="<<ishell<<"/"<<nshells()-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers()
                        <<" Fnu="<<layer_flux<<" mJy \n";

                /// loop over sublayer and evaluate the intensity distribution
                for(size_t iilayer = 0; iilayer < nsublayers; iilayer++){
//                    size_t cil = EjectaID2::CellsInLayer(ii);
                    double ctheta = theta_c_l[ii]+0.5*(theta_c_h[ii]-theta_c_l[ii]);
                    size_t nphi = cils[ii];
                    for (size_t j = 0; j < nphi; j++){
                        double cphis = (double)j * 2.0 * CGS::pi / (double)nphi;
                        tmp_pj.gerArr(IMG::Q::icphi)[j] = cphis;
                        tmp_cj.gerArr(IMG::Q::icphi)[j] = cphis;
                        tmp_pj.gerArr(IMG::Q::ictheta)[j] = ctheta;
                        tmp_cj.gerArr(IMG::Q::ictheta)[j] = ctheta;
                    }
                    if (tmp_pj.m_size==0||tmp_cj.m_size==0){
                        (*p_log)(LOG_ERR,AT)<<"error\n";
                        exit(1);
                    }
                    tmpImagesSet.getReferenceToTheImage(ii).m_size_active = nphi;
                    bw_rad->evalImageA(tmpImagesSet.getReferenceToTheImage(ii), tmp_pj, tmp_cj,
                                       theta_c_l[ii], theta_c_h[ii], nphi,
                                       obs_time, obs_freq, atol);
                    /// [debug] trying to account for if thetamax > theta_w hist image does not work...
//                    auto & ints = tmpImagesSet.getReferenceToTheImage(ii).gerArr(IMG::iintens);
//                    for (auto & int_ : ints) int_ /=(double)cils[ii];
//                    for (auto & int_ : ints) int_ /=(im_max_theta/(theta_c_h[ii]-theta_c_l[ii]));


//                    if (tmpImagesSet.getReferenceToTheImage(ilayer).m_f_tot == 0){
//                        n_jet_empty_images += 1; n_empty_images_layer.emplace_back(ilayer);
//                    }
                    /// override the image total flux density (normalize to the number of sublayers)
//                    if (tmpImagesSet.getReferenceToTheImage(ii).m_f_tot==0){
//                        (*p_log)(LOG_WARN,AT) << "image summed intensity = 0 il="<<ilayer<<" sub_il="<<iilayer<<" nphis="<<nphi<<"\n";
//                    }
                    tmpImagesSet.getReferenceToTheImage(ii).m_f_tot = layer_flux / (double)(nsublayers);
                    ii++;
                }
            }

            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
            images.getReferenceToTheImage(ishell).m_f_tot = 0.;
            auto & imageForEntireShell = images.getReferenceToTheImage(ishell);
//            if (std::accumulate(imageForEntireShell.gerArr(IMG::Q::iintens).begin(),
//                                imageForEntireShell.gerArr(IMG::Q::iintens).end(),0.)==0.){
//                (*p_log)(LOG_WARN,AT) << "image summed intensity = 0\n";
//                exit(1);
//            }
            imageForEntireShell.resize(2 * ncells, 0. );
            combineImagesA(imageForEntireShell, ncells, nlayers_ * nsublayers, tmpImagesSet) ;

            if (std::accumulate(imageForEntireShell.gerArr(IMG::Q::iintens).begin(),
                                imageForEntireShell.gerArr(IMG::Q::iintens).end(),0.)==0.){
                (*p_log)(LOG_WARN,AT) << "image summed intensity = 0\n";
                exit(1);
            }
            if (imageForEntireShell.m_f_tot == 0){
                (*p_log)(LOG_WARN,AT) << "image summed flux = 0\n";
                exit(1);
            }
            int i = 0;

//            if (imageForEntireShell.m_f_tot != tot_flux){
//                (*p_log)(LOG_ERR,AT)<<" imageForEntireShell.m_f_tot="<<imageForEntireShell.m_f_tot
//                    <<" tot_flux="<<tot_flux<<"\n";
//                exit(1);
//            }
//            double atol = tot_flux * rtol / (double)nlayers();
//            double layer_flux = bw_rad->evalFluxDensA(obs_time,obs_freq, atol);
//            imageForEntireShell.m_f_tot = layer_flux;
        }

        /// print which layers/shells gave isEmpty image
//        if (p_log->getLogLevel() == LOG_INFO) {
//            if (n_jet_empty_images > 0) {
//                auto &ccerr = std::cout;
//                ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an isEmpty images for total n="
//                      << n_jet_empty_images << " layers. Specifically:\n";
//                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
////                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
////                size_t n_layers_ej = m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
//                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [ ";
//                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
//                        ccerr << n_empty_images[ish][il] << " ";
//                    }
//                    ccerr << "] / (" << nlayers_ << " total layers) \n";
//                }
//            }
//        }

//        return std::move( images );
    }

    /// Compute Skymap for adaptive EATS
    void computeEjectaSkyMapA(Images & images, double obs_time, double obs_freq, size_t nsublayers ){
        size_t nlayers_ = nlayers();
        // ------------------------
        if (images.isEmpty()){
            (*p_log)(LOG_ERR,AT) << " isEmpty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nlayers_){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }
        // ----------------------
        double rtol = 1e-6;
        Images tmpImagesSet (nsublayers, IMG::m_names.size());
        size_t ii = 0;
        double tot_flux = 0.;
        std::vector<size_t> cils(nlayers_ * nsublayers, 100);
        EjectaID2::_evalCellsInLayer(nlayers_ * nsublayers, cils);

//        size_t ncells = EjectaID2::_evalTotalNcells(nlayers_ * nsublayers);
        size_t ncells = std::accumulate(cils.begin(), cils.end(), decltype(cils)::value_type(0));
        tmpImagesSet.resizeEachImage(ncells*2);
        Image tmp_pj( ncells, IMG::m_names.size(), 0, m_loglevel);
        Image tmp_cj( ncells, IMG::m_names.size(), 0, m_loglevel);

        Vector theta_c{};
        Vector theta_c_l{};
        Vector theta_c_h{};

        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer){

            auto & bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
            double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;

//        EjectaID2::_init_a_grid(theta_c_l, theta_c_h, theta_c, nlayers_ * nsublayers, CGS::pi/2.);
            EjectaID2::_init_pw_grid_(theta_c_l, theta_c_h, theta_c, nsublayers, theta_l, im_max_theta);//id->theta_wing);

            double atol = tot_flux * rtol / (double)nlayers();
            double layer_flux = bw_rad->evalFluxDensA(obs_time,obs_freq, atol);
            tot_flux += layer_flux;
            (*p_log)(LOG_INFO,AT)
                    << " EJECTA LC obs_time="<<obs_time<<" obs_freq="<<obs_freq
                    << " theta_layer="<<ilayer<<"/"<<nlayers()
                    << " theta_l="<<theta_l
                    << " Fnu="<<layer_flux<<" mJy \n";
            for(size_t iilayer = 0; iilayer < nsublayers; iilayer++){
                double ctheta = theta_c_l[iilayer];// + 0.5 * (theta_c_h[iilayer] - theta_c_l[iilayer]);
                size_t nphi = cils[ii];
                for (size_t j = 0; j < nphi; j++){
                    double cphis = (double)j * 2.0 * CGS::pi / (double)nphi;
                    tmp_pj.gerArr(IMG::Q::icphi)[j] = cphis;
                    tmp_cj.gerArr(IMG::Q::icphi)[j] = cphis;
                    tmp_pj.gerArr(IMG::Q::ictheta)[j] = ctheta;
                    tmp_cj.gerArr(IMG::Q::ictheta)[j] = ctheta;
                }
                if (tmp_pj.m_size==0||tmp_cj.m_size==0){
                    (*p_log)(LOG_ERR,AT)<<"error\n";
                    exit(1);
                }
                tmpImagesSet.getReferenceToTheImage(iilayer).m_size_active = nphi;
                bw_rad->evalImageA(tmpImagesSet.getReferenceToTheImage(iilayer), tmp_pj, tmp_cj,
                                   theta_c_l[iilayer], theta_c_h[iilayer], nphi,
                                   obs_time, obs_freq, atol);
                tmpImagesSet.getReferenceToTheImage(iilayer).m_f_tot = layer_flux / (double)(nsublayers);
                ii++;
            }
            images.getReferenceToTheImage(ilayer).m_f_tot = 0.;
            auto & imageForEntireShell = images.getReferenceToTheImage(ilayer);

//            size_t ncells_max = 0;
//            for (size_t iim=0; iim<tmpImagesSet.size(); iim++){
//                size_t ncells_ = 0;
//                auto & im = tmpImagesSet.getReferenceToTheImage(iim);
//                for (size_t ii_ = 0; ii_ < im.m_size; ii_++){
//                    if (im.gerArr(IMG::Q::iintens)[ii_]>0){
//                        ncells_++;
//                    }
//                }
//                if (ncells_max<ncells_)
//                    ncells_max=ncells_;
//            }
//            if (ncells_max < ncells){
//                (*p_log)(LOG_INFO,AT)<<" ncells_max="<<ncells_max<<" ncells="<<ncells<<"\n";
//            }

            imageForEntireShell.resize(2 * ncells, 0. );
            combineImagesA(imageForEntireShell, ncells, nsublayers, tmpImagesSet) ;

//            if (std::accumulate(imageForEntireShell.gerArr(IMG::Q::iintens).begin(),
//                                imageForEntireShell.gerArr(IMG::Q::iintens).end(),0.)==0.){
//                (*p_log)(LOG_WARN,AT) << "image summed intensity = 0\n";
//                exit(1);
//            }
            if (imageForEntireShell.m_f_tot == 0){
                (*p_log)(LOG_WARN,AT) << "image summed flux = 0\n";
                exit(1);
            }
        }
    }

    /// Compute lightcurve for piece-wise EATS
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
//        double flux_pj, flux_cj; size_t ii = 0;
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
//                ii ++;
            }
        }
        return std::move( light_curves );
    }


};


#endif //SRC_EJECTA_BASE_H
