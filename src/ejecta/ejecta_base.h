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

/// select type for the BW (controls which RHS to use)
BW_TYPES select_bw_type(StrStrMap opts, std::unique_ptr<logger> & p_log){

    std::string opt = "bw_type";
    BW_TYPES rhs_type;
    if ( opts.find(opt) == opts.end() ) {
        (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. No default option set\n";
        exit(1);
//            rhs_type = RHS_TYPES::iGRG_FS;
    }
    else{
        if(opts.at(opt) == "fs")
            rhs_type = BW_TYPES::iFS;
        else if(opts.at(opt) == "fsrs")
            rhs_type = BW_TYPES::iFSRS;
        else if(opts.at(opt) == "fs_dense")
            rhs_type = BW_TYPES::iFS_DENSE;
        else if(opts.at(opt) == "fs_dense_pwn")
            rhs_type = BW_TYPES::iFS_PWN_DENSE;
        else{
            (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                 <<" given: " << opts.at(opt)
                                 << " is not recognized. "
                                 << "Possible options: "
                                 << " grb_fs "<< " grb_fsrs " << " ej " << " ej_pwn " << "\n";
//                std::cerr << AT << "\n";
            exit(1);
        }
    }
    return rhs_type;
}


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
        if ((!m_pars.empty()) || (!m_opts.empty())) {
            run_bws   = getBoolOpt("run_bws", m_opts, AT, p_log, false, true);
            save_dyn  = getBoolOpt("save_dynamics", m_opts, AT, p_log, false, true);
            load_dyn  = getBoolOpt("load_dynamics", m_opts, AT, p_log, false, true);
            do_ele    = getBoolOpt("do_ele", m_opts, AT, p_log, false, true);
            do_spec   = getBoolOpt("do_spec", m_opts, AT, p_log, false, true);
            save_spec = getBoolOpt("save_spec", m_opts, AT, p_log, false, true);
            do_lc     = getBoolOpt("do_lc", m_opts, AT, p_log, false, true);
            do_skymap = getBoolOpt("do_skymap", m_opts, AT, p_log, false, true);
            /// copy parameters from the main to ejecta (same for all ejecta types)
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
        else
            (*p_log)(LOG_INFO, AT) << "ejecta is not initialized and will not be considered.\n";

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
                BW_TYPES m_type = select_bw_type(opts,p_log);
                p_cumShells.push_back(
                        std::make_unique<CumulativeShell>(Vector {}, nshells(),
                                                          il, n_substeps,m_type,
                                                          p_log->getLogLevel())
                                                          );
                /// set parameters for the shell (concerns BW interaction, structure)
                p_cumShells[il]->setPars(pars, opts);
                for (size_t ish = 0; ish < nshells(); ish++){
                    auto & bw = p_cumShells[il]->getBW(ish);
                    /// set parameters for each blast wave within the shell
                    bw->setParams(id, pars, opts, il, ii_eq);
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
            BW_TYPES type = select_bw_type(opts, p_log);
            p_cumShells.push_back(
                    std::make_unique<CumulativeShell>(t_arr, nshells_, il, n_substeps,
                                                      type, p_log->getLogLevel()) );
            p_cumShells[il]->setPars(pars, opts);
            empty_bws.emplace_back("il="+std::to_string(il)+" | shells:");
            for (size_t ish = 0; ish < nshells_; ish++){
                auto & bw = p_cumShells[il]->getBW(ish);

                /// If ID is not correct, skip this blastwave initiialization (will not be evolved)
                if (id->get(ish,il,EjectaID2::Q::ir) <= 0 ||
                    id->get(ish,il,EjectaID2::Q::iek) <= 0 ||
                    EQS::Beta(EQS::GamFromMom(id->get(ish,il,EjectaID2::Q::imom))) <= 1.e-6 ||
                    id->get(ish,il,EjectaID2::Q::imass) <= 0){
                    empty_bws[il]+=" "+std::to_string(ish);
                    n_unitinitilized_shells++;
                    bw->getPars()->end_evolution = true;
                    ii_eq += SOL::neq;//bw->getNeq();
                    continue;
                }
                /// init parameters for blastwave (incl. initial conditions)
                else{
                    bw->setParams(id, pars, opts, il, ii_eq);
                    ii_eq += SOL::neq;//bw->getNeq();
                }

                if (bw->getPars()->which_jet_layer_to_use == 0){
                    bw->getPars()->which_jet_layer_to_use = 0; // the fastest
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
            }
        }

//        is_ejBW_init = true;
        is_ejecta_obs_pars_set = true;

        if ((n_unitinitilized_shells > 0)){
            (*p_log)(LOG_WARN,AT)<<"-------------- NO ID FOR ------------"<<"\n";
            for (const auto & empty_bw : empty_bws){
                (*p_log)(LOG_WARN,AT)<<"\t NOT EVOLVING: "<<empty_bw<<"\n";
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

//        ej_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);

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
    void computeEjectaSkyMapPW(ImageExtend & im, double obs_time, double obs_freq){
        /// out is [i_vn][ish][itheta_iphi]
        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();
        size_t ncells_ =  (int)ncells();


        size_t n_jet_empty_images = 0;
        for (size_t ishell = 0; ishell < nshells_; ishell++){
            auto & out = im.raw_data[ishell];
            /// Allocate memory for ctheta and cphi grids ('2' is for principle and counter jet)
            for (size_t ivn = 0; ivn < IMG::m_names.size(); ivn++)
                out[ivn].resize(2 * ncells_, 0.);
            /// compute the skymap for each layer (stacking cells batch after batch along one axis)
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
                im.fluxes_shells[ishell] += flux;
                offset += cil;
            }
            if ((im.fluxes_shells[ishell] <=0) || (!std::isfinite(im.fluxes_shells[ishell]))){

                for (size_t il = 0; il < nlayers_; il++){
                    auto & pars = p_cumShells[il]->getBW(ishell)->getPars();
                    auto & bw = p_cumShells[il]->getBW(ishell);
                    (*p_log)(LOG_WARN,AT)<<"\t il="<<il<<" G0="<<pars->Gamma0<<" E0="<<pars->E0<<" i_end="
                            <<pars->i_end_r<<" tb[iend]="<<bw->getData(BW::itburst, pars->i_end_r-1)<<"\n";
                }

                (*p_log)(LOG_WARN,AT)<<" im.fluxes_shells[ishell]="<<im.fluxes_shells[ishell]
                    <<" ish="<<ishell<< " obs_time="<<obs_time<<"  obs_freq="<<obs_freq<<""<<"\n";
                im.total_flux += im.fluxes_shells[ishell];
                return;
            }

            /// find the image extend for this shell
            double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
            double imin = std::numeric_limits<double>::max(), imax = std::numeric_limits<double>::min();
            /// TODO replace with algorithmic search
            for (size_t i = 0; i < 2 * ncells_; i++){
                if (out[IMG::Q::ixr][i] < xmin) xmin = out[IMG::Q::ixr][i];
                if (out[IMG::Q::ixr][i] > xmax) xmax = out[IMG::Q::ixr][i];
                if (out[IMG::Q::iyr][i] < ymin) ymin = out[IMG::Q::iyr][i];
                if (out[IMG::Q::iyr][i] > ymax) ymax = out[IMG::Q::iyr][i];
                if (out[IMG::Q::iintens][i] < imin) imin = out[IMG::Q::iintens][i];
                if (out[IMG::Q::iintens][i] > imax) imax = out[IMG::Q::iintens][i];
            }
            ///
            if ((imax <= 0) || (!std::isfinite(imax))){
                (*p_log)(LOG_ERR,AT) << "  imax="<<imax<<"\n";
                exit(1);
            }

            /// normalize image to be mJy/mas^2
            double delta_x = xmax - xmin;
            double delta_y = ymax - ymin;
            double dfnu = im.fluxes_shells[ishell] / (delta_x * delta_y);
            for (size_t i = 0; i < 2 * ncells_; i++){
                out[IMG::Q::iintens][i] *= (dfnu / imax);
            }

            /// update image boundaries for overall image
            if (xmin < im.xmin) im.xmin = xmin;
            if (ymin < im.ymin) im.ymin = ymin;
            if (xmax > im.xmax) im.xmax = xmax;
            if (ymax > im.ymax) im.ymax = ymax;
            im.total_flux += im.fluxes_shells[ishell];
        }
    }

    bool _computeEjectaSkyMapA(std::vector<VecVector> & row_data,
                               double obs_time, double obs_freq,
                               size_t & nsublayers, size_t & nrestarts,
                               Vector & th_b_layers, double th_jet_max){

        size_t min_sublayers = 3;
        double frac_to_increase = 1.5;
        size_t max_restarts = 10;
        size_t min_non_zero_cells = 3;

        /// allocate memory for angular grid for this shell
        Vector theta_pw(nsublayers + 1);
        for (size_t i = 0; i < nsublayers + 1; i++) {
            double fac = (double) i / (double) nsublayers;
            theta_pw[i] = 2.0 * std::asin(fac * std::sin(th_jet_max / 2.0));
        }
        Vector cthetas(nsublayers, 0.);

        size_t nlayers_ = nlayers();
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer) {
            /// get data conter for this time/freq/shell
            auto & out = row_data[ilayer];
            /// rotate index to the next time/freq/shell
            size_t ncells = 0;

            /// find how many sublayers lie withing jet openning angle and how many phi-cells are there in total
            std::vector<size_t> sublayers_l{};
            auto &bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
            double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;
            size_t ncells_il = 0;
            for (size_t ith = 0; ith < nsublayers; ith++) {
                cthetas[ith] = theta_pw[ith] + (theta_pw[ith + 1] - theta_pw[ith]) * 0.5;
                size_t cil = EjectaID2::CellsInLayer_(ith);
                if ((cthetas[ith] >= theta_l) // check that BW lower angle is above the ctheta[il]
                    &&
                    (cthetas[ith] < th_b_layers[ilayer]) // check if BW maximum openning angle is below ctheta[il]
                        ) {
                    ncells += cil;
                    sublayers_l.push_back(ith);
                }
            }
            if (sublayers_l.size() < min_sublayers) {
                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"] "
                                       <<" Reason: nsublayers="<<sublayers_l.size()<<" < "<<"min_sublayers="<<min_sublayers<<" [restarts="<<nrestarts<<"]"
                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                if (nrestarts>max_restarts){
                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                    exit(1);
                }
                nsublayers = nsublayers_new;
                nrestarts+=1;
                return false;
            }

            /// resize the output according to the number of cell to compute
            for (size_t ivn = 0; ivn < IMG::m_names.size(); ivn++) {
//                if (out[ivn].size() < 2 * ncells)
//                    out[ivn].clear();
//                std::cout << out[ivn].size() << "\t" << ncells<< "\n";
                out[ivn].resize(2 * ncells);
            }

            /// compute the intensty for these cells
            size_t iii = 0;
            double layer_summed_intensity = 0;
            for (auto &ith: sublayers_l) {
                size_t cil = EjectaID2::CellsInLayer_(ith);
                /// fill the angular grid for intensity evaluations
                for (size_t iphi = 0; iphi < cil; iphi++) {
                    double cphi = (double) iphi * 2.0 * CGS::pi / (double) cil;
                    /// principle jet
                    out[IMG::Q::icphi][iii + iphi] = cphi;
                    out[IMG::Q::ictheta][iii + iphi] = cthetas[ith];
                    /// counter jet
                    out[IMG::Q::icphi][ncells + iii + iphi] = cphi;
                    out[IMG::Q::ictheta][ncells + iii + iphi] = cthetas[ith];
                }
                /// evaluate intensity
                layer_summed_intensity += bw_rad->evalSkyMapA(out, obs_time, obs_freq, ilayer, iii, cil, ncells);
                iii = iii + cil;
            }
            /// check if the principle jet an counter jets are resolved (intensity found)
            double summed_int_pj = 0.;
            double summed_int_cj = 0.;
            int nx_pj=0, nx_cj=0;
            for (size_t i = 0; i < ncells; i++) {
                summed_int_pj += out[IMG::Q::iintens][i];
                summed_int_cj += out[IMG::Q::iintens][ncells + i];
                if (out[IMG::Q::iintens][i] > 0) nx_pj++;
                if (out[IMG::Q::iintens][ncells + i] > 0) nx_cj++;
            }

            if ((summed_int_pj == 0) || (summed_int_cj == 0)) {
//                    (*p_log)(LOG_ERR, AT)
//                            << " sum(summed_int_pj)=" << summed_int_pj << " sum(summed_int_cj)=" << summed_int_cj
//                            << "  for layer=" << ilayer << " / " << nlayers_ << " with nsublayers=" << nsublayers
//                            << "  layer_flux=" << fluxes[ilayer]
//                            << "  ID layer theta_c_l=" << id->get(0, ilayer, EjectaID2::Q::itheta_c_l)
//                            << "  sublayer grid ctheta=[" << cthetas[0] << " " << cthetas[nsublayers - 1]
//                            << "  jet edge th_b=" << th_b_layers[ilayer]
//                            << "\n";
                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                                       <<" Reason: summed_int_pj=0 or summed_int_cj=0 [restarts="<<nrestarts<<"]"
                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                if (nrestarts>max_restarts){
                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                    exit(1);
                }
                nsublayers = nsublayers_new;
                nrestarts+=1;
                return false;
//                    exit(1);
            }
            if (nx_pj < min_non_zero_cells){
                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                                       <<" Reason: nx_pj="<<nx_pj<<" < min_non_zero_cells="<<min_non_zero_cells<<" [restarts="<<nrestarts<<"]"
                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                if (nrestarts>max_restarts){
                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                    exit(1);
                }
                nsublayers = nsublayers_new;
                nrestarts+=1;
                return false;
            }
            if (nx_cj < min_non_zero_cells){
                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                                       <<" Reason: nx_cj="<<nx_cj<<" < min_non_zero_cells="<<min_non_zero_cells<<" [restarts="<<nrestarts<<"]"
                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                if (nrestarts>max_restarts){
                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                    exit(1);
                }
                nsublayers = nsublayers_new;
                nrestarts+=1;
                return false;
            }

        }
        return true;
    }

    /// Compute Skymap for adaptive EATS
    void computeEjectaSkyMapA(ImageExtend & im, double obs_time, double obs_freq, size_t nsublayers_, size_t nphi ){

        /// out is [i_vn][ish][itheta_iphi]
        size_t nlayers_ = nlayers();
        // ------------------------

        double rtol = 1.e-6;

        (*p_log)(LOG_INFO,AT)  << " Computing skymap with N layers=" << nlayers_ << " and N sublayers=" << nsublayers_ << "\n";

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
        (*p_log)(LOG_INFO,AT)  << " Found outermost theta_edge="<<th_jet_max<<" at obs_time="<<obs_time<<" obs_freq="<<obs_freq << "\n";

        /// compute flux density for each layer ( for normalization )
        double tot_flux = 0.;
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer) {
            auto &bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
//            double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;
            double atol = tot_flux * rtol / (double) nlayers_;
            double layer_flux = bw_rad->evalFluxDensA(obs_time, obs_freq, atol);
            tot_flux += layer_flux;
            im.fluxes_shells[ilayer] = layer_flux;
        }
        /// make an interpolator for sublayers (for a given ctheta -> get flux)
        auto interp = Interp1d(id->getVec(0,EjectaID2::Q::itheta_c_l), im.fluxes_shells);

        /// find non-zero sublayers in each layer and corresponding cells
        size_t nrestarts = 0;

        /// Depending on the jet layer opening angle (which depends on spreading) it is possible that
        /// no jet sublayer will fall within the grid resolution for computing intensity.
        /// So the entire layer will not be resolved. In order to avoid this
        /// we iteratively increase the resolution until the layer is resolved.
        /// this value can be different for different shells
        /// TODO optimize memory allocation

        bool is_ok = _computeEjectaSkyMapA(im.raw_data, obs_time, obs_freq,
                                           nsublayers_, nrestarts, th_b_layers, th_jet_max);
        while (!is_ok){
            is_ok = _computeEjectaSkyMapA(im.raw_data, obs_time, obs_freq,
                                          nsublayers_, nrestarts, th_b_layers, th_jet_max);
        }

        /// normalize the image
        for (size_t il = 0; il < nlayers_; il++){
            auto & out = im.raw_data[il];

            /// normalize the intensity (divide by the size of the cell)
            double max_int = std::numeric_limits<double>::min();
            double total_int = 0.;
            /// find overall maximum (in principle and counter jet)
            for (size_t i = 0; i < out[IMG::Q::iintens].size(); i++) {
                total_int += out[IMG::Q::iintens][i];
                if (max_int < out[IMG::Q::iintens][i])
                    max_int = out[IMG::Q::iintens][i];

            }
            /// size of the principle jet
            double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
            for (size_t i = 0; i < out[IMG::Q::ixr].size(); i++) {
                if (out[IMG::Q::ixr][i] < xmin) xmin = out[IMG::Q::ixr][i];
                if (out[IMG::Q::ixr][i] > xmax) xmax = out[IMG::Q::ixr][i];
                if (out[IMG::Q::iyr][i] < ymin) ymin = out[IMG::Q::iyr][i];
                if (out[IMG::Q::iyr][i] > ymax) ymax = out[IMG::Q::iyr][i];
            }
            if (xmax == xmin) {
                (*p_log)(LOG_ERR, AT) << " xmin=xmax="<<xmin<<"\n";
                exit(1);
//                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
//                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
//                                       <<" Reason: xmax=xmin [restarts="<<nrestarts<<"]"
//                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
//                if (nrestarts>max_restarts){
//                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
//                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
//                    exit(1);
//                }
//                nsublayers = nsublayers_new;
//                nrestarts+=1;
//                return false;
            }
            if (ymax == ymin) {
                (*p_log)(LOG_ERR, AT) << " ymin=ymax="<<ymin<<"\n";
                exit(1);
//                auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
//                (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
//                                       <<" Reason: ymax=ymin [restarts="<<nrestarts<<"]"
//                                       <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
//                if (nrestarts>max_restarts){
//                    (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
//                                         << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
//                    exit(1);
//                }
//                nsublayers = nsublayers_new;
//                nrestarts+=1;
//                return false;
            }
            double delta_x = xmax - xmin;
            double delta_y = ymax - ymin;

            /// normalize to mJy/mas^2
            for (size_t i = 0; i < out[IMG::Q::iintens].size(); i++) {
                out[IMG::Q::iintens][i] = out[IMG::Q::iintens][i] / max_int * im.fluxes_shells[il] / delta_x / delta_y;
            }

            /// save the maximum extend
            if (im.xmin > xmin) im.xmin = xmin;
            if (im.ymin > ymin) im.ymin = ymin;
            if (im.xmax < xmax) im.xmax = xmax;
            if (im.ymax < xmax) im.ymax = ymax;
            im.total_flux += im.fluxes_shells[il];
        }
#if 0
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer) {

            /// get data conter for this time/freq/shell
            auto &out = out_[itnush];
            /// rotate index to the next time/freq/shell
            itnush += 1;
            /// Depending on the jet layer opening angle (which depends on spreading) it is possible that
            /// no jet layer will fall within the grid resolution for computing intensity. In order to avoid this
            /// we iteratively increase the resolution until the layer is resolved.
            size_t nsublayers = nsublayers_; // initial guess for the sufficient resolution
            size_t nrestarts = 0;
            bool found = false;
            while (!found) {

                /// this value can be different for different shells
                size_t ncells = 0;

                /// allocate memory for angular grid for this shell
                Vector theta_pw(nsublayers + 1);
                for (size_t i = 0; i < nsublayers + 1; i++) {
                    double fac = (double) i / (double) nsublayers;
                    theta_pw[i] = 2.0 * std::asin(fac * std::sin(th_jet_max / 2.0));
                }
                Vector cthetas(nsublayers, 0.);

                /// find how many sublayers lie withing jet openning angle and how many phi-cells are there in total
                std::vector<size_t> sublayers_l{};
                auto &bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();
                double theta_l = p_cumShells[ilayer]->getBW(0)->getPars()->theta_c_l;
                size_t ncells_il = 0;
                for (size_t ith = 0; ith < nsublayers; ith++) {
                    cthetas[ith] = theta_pw[ith] + (theta_pw[ith + 1] - theta_pw[ith]) * 0.5;
                    size_t cil = EjectaID2::CellsInLayer_(ith);
                    if ((cthetas[ith] >= theta_l) // check that BW lower angle is above the ctheta[il]
                        &&
                        (cthetas[ith] < th_b_layers[ilayer]) // check if BW maximum openning angle is below ctheta[il]
                            ) {
                        ncells += cil;
                        sublayers_l.push_back(ith);
                    }
                }
                if (sublayers_l.size() < min_sublayers) {
                    auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                    (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"] "
                        <<" Reason: nsublayers="<<sublayers_l.size()<<" < "<<"min_sublayers="<<min_sublayers<<" [restarts="<<nrestarts<<"]"
                        <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                    if (nrestarts>max_restarts){
                        (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                            << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                        exit(1);
                    }
                    nsublayers = nsublayers_new;
                    nrestarts+=1;
                    continue;
                }

                /// resize the output according to the number of cell to compute
                for (size_t ivn = 0; ivn < IMG::m_names.size(); ivn++) {
                    out[ivn].resize(2 * ncells, 0.);
                }

                /// compute the intensty for these cells
                size_t iii = 0;
                double layer_summed_intensity = 0;
                for (auto &ith: sublayers_l) {
                    size_t cil = EjectaID2::CellsInLayer_(ith);
                    /// fill the angular grid for intensity evaluations
                    for (size_t iphi = 0; iphi < cil; iphi++) {
                        double cphi = (double) iphi * 2.0 * CGS::pi / (double) cil;
                        /// principle jet
                        out[IMG::Q::icphi][iii + iphi] = cphi;
                        out[IMG::Q::ictheta][iii + iphi] = cthetas[ith];
                        /// counter jet
                        out[IMG::Q::icphi][ncells + iii + iphi] = cphi;
                        out[IMG::Q::ictheta][ncells + iii + iphi] = cthetas[ith];
                    }
                    /// evaluate intensity
                    layer_summed_intensity += bw_rad->evalSkyMapA(out, obs_time, obs_freq, ilayer, iii, cil, ncells);
                    iii = iii + cil;

                }
                /// check if the principle jet an counter jets are resolved (intensity found)
                double summed_int_pj = 0.;
                double summed_int_cj = 0.;
                for (size_t i = 0; i < ncells; i++) {
                    summed_int_pj += out[IMG::Q::iintens][i];
                    summed_int_cj += out[IMG::Q::iintens][ncells + i];
                }
                if ((summed_int_pj == 0) || (summed_int_cj == 0)) {
//                    (*p_log)(LOG_ERR, AT)
//                            << " sum(summed_int_pj)=" << summed_int_pj << " sum(summed_int_cj)=" << summed_int_cj
//                            << "  for layer=" << ilayer << " / " << nlayers_ << " with nsublayers=" << nsublayers
//                            << "  layer_flux=" << fluxes[ilayer]
//                            << "  ID layer theta_c_l=" << id->get(0, ilayer, EjectaID2::Q::itheta_c_l)
//                            << "  sublayer grid ctheta=[" << cthetas[0] << " " << cthetas[nsublayers - 1]
//                            << "  jet edge th_b=" << th_b_layers[ilayer]
//                            << "\n";
                    auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                    (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                        <<" Reason: summed_int_pj=0 or summed_int_cj=0 [restarts="<<nrestarts<<"]"
                        <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                    if (nrestarts>max_restarts){
                        (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                             << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                        exit(1);
                    }
                    nsublayers = nsublayers_new;
                    nrestarts+=1;
                    continue;
//                    exit(1);
                }

                /// normalize the intensity (divide by the size of the cell)
                double max_int = std::numeric_limits<double>::min();
                double total_int = 0.;
                /// find overall maximum (in principle and counter jet)
                for (size_t i = 0; i < 2 * ncells; i++) {
                    total_int += out[IMG::Q::iintens][i];
                    if (max_int < out[IMG::Q::iintens][i])
                        max_int = out[IMG::Q::iintens][i];

                }
                /// size of the principle jet
                double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
                double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
                for (size_t i = 0; i < 2 * ncells; i++) {
                    if (out[IMG::Q::ixr][i] < xmin) xmin = out[IMG::Q::ixr][i];
                    if (out[IMG::Q::ixr][i] > xmax) xmax = out[IMG::Q::ixr][i];
                    if (out[IMG::Q::iyr][i] < ymin) ymin = out[IMG::Q::iyr][i];
                    if (out[IMG::Q::iyr][i] > ymax) ymax = out[IMG::Q::iyr][i];
                }
                if (xmax == xmin) {
                    auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                    (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                                           <<" Reason: xmax=xmin [restarts="<<nrestarts<<"]"
                                           <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                    if (nrestarts>max_restarts){
                        (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                             << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                        exit(1);
                    }
                    nsublayers = nsublayers_new;
                    nrestarts+=1;
                    continue;
                }
                if (ymax == ymin) {
                    auto nsublayers_new = size_t((double)nsublayers * frac_to_increase);
                    (*p_log)(LOG_WARN, AT) << "\t restarting skymap [il="<<ilayer<<"]"
                                           <<" Reason: ymax=ymin [restarts="<<nrestarts<<"]"
                                           <<" Increasing nsublayers="<<nsublayers<<"->"<<nsublayers_new<<"\n";
                    if (nrestarts>max_restarts){
                        (*p_log)(LOG_ERR, AT)<<" maximum number of restarts exceeded for [il="<<ilayer<<"] "
                                             << " nrestarts=" << nrestarts << " sublayers=" << nsublayers <<"\n";
                        exit(1);
                    }
                    nsublayers = nsublayers_new;
                    nrestarts+=1;
                    continue;
                }
                double delta_x = xmax - xmin;
                double delta_y = ymax - ymin;

                /// normalize jet to its size
                double d_l = p_cumShells[ilayer]->getBW(0)->getPars()->d_l;
                for (size_t i = 0; i < 2 * ncells; i++) {
                    out[IMG::Q::iintens][i] = out[IMG::Q::iintens][i] / max_int * fluxes[ilayer] / delta_x / delta_y;
                }
                found = true;
            }
        }

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

        /// evaluate intensity distributions
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer) {

            auto &bw_rad = p_cumShells[ilayer]->getBW(0)->getFsEATS();

            size_t iii = 0;
            double layer_summed_intensity = 0;
            /// loop for all thetas from 0 to overall theta_edge # TODO see if you can just do it until current theta_edge

            for (size_t ith = 0; ith < nsublayers; ith++) {
                size_t cil = EjectaID2::CellsInLayer_(ith);

                /// fill the angular grid for intensity evaluations
                for (size_t iphi = 0; iphi < cil; iphi++) {
                    double cphi = (double) iphi * 2.0 * CGS::pi / (double) cil;
                    /// principle jet
                    out[IMG::Q::icphi][ilayer][iii + iphi] = cphi;
                    out[IMG::Q::ictheta][ilayer][iii + iphi] = cthetas[ith];
                    /// counter jet
                    out[IMG::Q::icphi][ilayer][ncells + iii + iphi] = cphi;
                    out[IMG::Q::ictheta][ilayer][ncells + iii + iphi] = cthetas[ith];
                }

                /// evaluate intensity
                layer_summed_intensity += bw_rad->evalSkyMapA(out, obs_time, obs_freq, ilayer, iii, cil, ncells);
                iii = iii + cil;
            }

            double summed_int_pj = 0.;
            double summed_int_cj = 0.;
            for (size_t i = 0; i < ncells; i++) {
                summed_int_pj+=out[IMG::Q::iintens][ilayer][i];
                summed_int_cj+=out[IMG::Q::iintens][ilayer][ncells + i];
            }
            if ((summed_int_pj == 0) || (summed_int_cj == 0)){
                (*p_log)(LOG_ERR,AT)
                    <<" sum(summed_int_pj)="<<summed_int_pj<<" sum(summed_int_cj)="<<summed_int_cj
                    <<"  for layer="<<ilayer<<" / "<<nlayers_ << " with nsublayers="<<nsublayers
                    <<"  layer_flux="<<fluxes[ilayer]
                    <<"  ID layer theta_c_l="<<id->get(0,ilayer,EjectaID2::Q::itheta_c_l)
                    <<"  sublayer grid ctheta=["<<cthetas[0]<<" "<<cthetas[nsublayers-1]
                    <<"  jet edge th_b="<<th_b_layers[ilayer]
                <<"\n";
                exit(1);
            }

        }


        ImageExtend im{};
        /// normalize the skymap
        for (size_t ilayer = 0; ilayer < nlayers_; ++ilayer){

            /// find overall maximum
            double max_int = std::numeric_limits<double>::min();
            double total_int = 0.;
            for (size_t i = 0; i < 2 * ncells; i++) {
                total_int += out[IMG::Q::iintens][ilayer][i];
                if (max_int < out[IMG::Q::iintens][ilayer][i])
                    max_int = out[IMG::Q::iintens][ilayer][i];

            }

            /// size of the principle jet
            double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::min();
            double ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::min();
            for (size_t i = 0; i < 2 * ncells; i++){
                if (out[IMG::Q::ixr][ilayer][i] < xmin) xmin = out[IMG::Q::ixr][ilayer][i];
                if (out[IMG::Q::ixr][ilayer][i] > xmax) xmax = out[IMG::Q::ixr][ilayer][i];
                if (out[IMG::Q::iyr][ilayer][i] < ymin) ymin = out[IMG::Q::iyr][ilayer][i];
                if (out[IMG::Q::iyr][ilayer][i] > ymax) ymax = out[IMG::Q::iyr][ilayer][i];
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

            /// normalize principle jet to its size
            double d_l = p_cumShells[ilayer]->getBW(0)->getPars()->d_l;
            for (size_t i = 0; i < 2 * ncells; i++) {
                out[IMG::Q::iintens][ilayer][i] = out[IMG::Q::iintens][ilayer][i] / max_int * fluxes[ilayer] / delta_x / delta_y;
            }



//            size_t iii = 0;
//            for (size_t ith = 0; ith < nsublayers; ith++){
//                size_t cil = EjectaID2::CellsInLayer_(ith);
//                double flux_ctheta = 0;
//                if (cthetas[ith] < id->getVec(0,EjectaID2::Q::itheta_c_l)[0])
//                    flux_ctheta = fluxes[0];
//                else if (cthetas[ith] > id->getVec(0,EjectaID2::Q::itheta_c_l)[id->nlayers-1])
//                    flux_ctheta = fluxes[nlayers_-1];
//                else
//                    flux_ctheta = interp.Interpolate(cthetas[ith], Interp1d::METHODS_SYNCH::iLinear);
//
//                if (flux_ctheta <= 0){
//                    std::cerr << AT << "error \n";
//                    exit(1);
//                }
//                for (size_t iphi = 0;  iphi < cil; iphi++){
//                    out[IMG::Q::iintens][ilayer][iii+iphi] = out[IMG::Q::iintens][ilayer][iii+iphi] / max_int * flux_ctheta / delta_x / delta_y;
//                    out[IMG::Q::iintens][ilayer][ncells+iii+iphi] = out[IMG::Q::iintens][ilayer][ncells+iii+iphi] / max_int * flux_ctheta / delta_x / delta_y;
//                }
//                iii = iii + cil;
//            }

            /// update image boundaries for overall image
//            if (xmin < im.xmin) im.xmin = xmin;
//            if (ymin < im.ymin) im.ymin = ymin;
//            if (xmax > im.xmax) im.xmax = xmax;
//            if (ymax > im.ymax) im.ymax = ymax;
        }
#endif
    }

    /// Compute lightcurve for piece-wise EATS
    void evalEjectaLightCurves(VecVector & out, Vector & obs_times, Vector & obs_freqs) {
        /// out // [ish*il][it*inu]
        (*p_log)(LOG_INFO, AT) << " starting ejecta light curve calculation\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++) {
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer];//ejectaModels[ishell][ilayer];
                auto & bw = model->getBW(ishell);
                if (bw->getPars()->i_end_r == 0) {
                    (*p_log)(LOG_WARN, AT)
                            << " NOT EVOLVED EJECTA BW"
                            << " vel_shell=" << ishell << "/" << nshells() - 1
                            << " theta_layer=" << ilayer << "/" << nlayers()
                            << " beta0=" << bw->getPars()->beta0
                            << "\n";
                }
                else {
                    (*p_log)(LOG_INFO, AT)
                            << " EJECTA LC ntimes=" << obs_times.size()
                            << " vel_shell=" << ishell << "/" << nshells() - 1
                            << " theta_layer=" << ilayer << "/" << nlayers()
                            << "\n";
                    bw->getFsEATS()->evalLightCurve(out[ii], id->method_eats, obs_times, obs_freqs);
                }
                ii++;
            }
        }
    }
};


#endif //SRC_EJECTA_BASE_H
