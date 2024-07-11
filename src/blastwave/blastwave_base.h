//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLASTWAVE_BASE_H
#define SRC_BLASTWAVE_BASE_H

#include "../composition.h"
//#include "blastwave_components.h"
#include "blastwave_pars.h"
#include "../microphysics/common.h"

class BlastWaveBase{
    std::unique_ptr<logger> p_log;
protected:
    Vector m_tb_arr;
    VecVector m_data{}; // container for the solution of the evolution
//    VecVector m_tmp_data{}; // container for the solution for each evolved step (for derivatives)
    VecVector mDtmp{};
    VecVector m_data_shells{};
    /// -----------------------------------------------------
    Pars * p_pars = nullptr;
    std::unique_ptr<LatSpread> p_spread = nullptr;
    std::unique_ptr<EOSadi> p_eos = nullptr;
    std::unique_ptr<RhoISM> p_dens = nullptr;
    std::unique_ptr<SedovTaylor> p_sedov = nullptr;
    std::unique_ptr<BlandfordMcKee2> p_bm = nullptr;
    std::unique_ptr<NuclearAtomic> p_nuc = nullptr;
    std::unique_ptr<LinearRegression> p_lr_delta = nullptr;
    std::unique_ptr<LinearRegression> p_lr_vol = nullptr;
//    std::unique_ptr<ShockMicrophysics> p_fs;
//    std::unique_ptr<ShockMicrophysics> p_rs;
    static constexpr int iters=1000; // for PWN, resolution of frac_psr_dep_
    Vector frac_psr_dep_{}; // for PWN, fraction of rad. absorbed by BW f(opacity)
    CommonTables & commonTables;

public:
    bool is_initialized = false;

    BlastWaveBase(Vector & tb_arr, size_t ishell, size_t ilayer, size_t n_substeps,
                  BW_TYPES type, CommonTables & commonTables, int loglevel)
        : m_tb_arr(tb_arr), commonTables(commonTables) {

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW_Base");

        /// parameters (pass p_fs, p_rs for later use in EATS()
        p_pars = new Pars(m_data, mDtmp, loglevel); //
        p_pars->m_type = type;

        /// the container for the final solution
        if (m_data.empty())
            m_data.resize(BW::TOTAL_VARS );

        /// the container for the last N substeps
        if (mDtmp.empty())
            mDtmp.resize( BW::TOTAL_VARS );

        /// Check if contener will be filled by evolving or loading
        if (m_tb_arr.empty())
            (*p_log)(LOG_WARN,AT) << " Time grid was not initialized\n";

        /// if no evolution required; do not allocate memory for each variable (only for evolved)
        if (m_data[BW::Q::itburst].size() < 1)
            for (auto & ivn : BW::VARS.at(p_pars->m_type))
                m_data[ivn].resize(tb_arr.size(), 0.0);

        // ---------------------- Methods ---------------
        p_lr_delta = std::make_unique<LinearRegression>(m_data[BW::Q::iR], m_data[BW::Q::iEJdelta]);
        p_lr_vol = std::make_unique<LinearRegression>(m_data[BW::Q::iR], m_data[BW::Q::iEJvol]);
        p_spread = std::make_unique<LatSpread>();
        p_eos = std::make_unique<EOSadi>();
        p_dens = std::make_unique<RhoISM>(loglevel);
        p_sedov = std::make_unique<SedovTaylor>();
        p_bm = std::make_unique<BlandfordMcKee2>();
        p_nuc = std::make_unique<NuclearAtomic>(loglevel);
//        p_fs = std::make_unique<ShockMicrophysics>(m_data, loglevel);
//        p_rs = std::make_unique<ShockMicrophysics>(m_data, loglevel);

        /// if no evolution required; do not allocate memory for each variable
        p_pars->n_substeps = n_substeps;
        if (mDtmp[BW::Q::itburst].size() < 1)
            for (auto & ivn : BW::VARS.at(p_pars->m_type))
                mDtmp[ivn].resize(p_pars->n_substeps, 0.0);

        /// ----------------------
        p_pars->loglevel = loglevel;
        p_pars->nr = m_tb_arr.size();
        p_pars->ilayer = ilayer;
        p_pars->ishell = ishell;
        p_pars->n_substeps = n_substeps;
        is_initialized = true;
    }

    void setBaseParams(std::unique_ptr<EjectaID2> & id, StrDbMap & pars, StrStrMap & opts, size_t ilayer, size_t ii_eq){

        // set parameters for ISM density
        double nism = getDoublePar("n_ism", pars, AT, p_log, -1, true);//pars.at("nism");
        double A0 = getDoublePar("A0", pars, AT,p_log,-1,false);//pars.at("A0");
        double s = getDoublePar("s", pars, AT,p_log,-1,false);//pars.at("s");
        double r_ej = getDoublePar("r_ej", pars, AT,p_log,-1,false);//pars.at("r_ej");
        double r_ism = getDoublePar("r_ism", pars, AT,p_log,-1,false);//pars.at("r_ism");
        p_dens->setPars(nism, A0, s, r_ej, r_ism, true);

        // spreading
        double a = getDoublePar("a", pars, AT,p_log,-1,false);//pars.at("a");
        double theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);//pars.at("theta_max");

        // radiative losses
        p_pars->eps_rad   = getDoublePar("epsilon_e_rad", pars, AT,p_log,0.,false);// pars.at("epsilon_e_rad");

        // Mass accretion from ISM; set options
        std::string opt;
        opt = "method_dmdr";
        METHOD_dmdr method_dmdr;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_dmdr = METHOD_dmdr::iusingdthdR;
        }
        else{
            if(opts.at(opt) == "usingA")
                method_dmdr = METHOD_dmdr::iusingA;
            else if(opts.at(opt) == "usingdthdr")
                method_dmdr = METHOD_dmdr::iusingdthdR;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                      << " given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " usingA " <<" usingdthdr "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_dmdr = method_dmdr;

        /// evolution eq
        opt = "method_dgdr";
        METHOD_dgdr method_dgdr;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_dgdr = METHOD_dgdr::iour;
        }
        else{
            if(opts.at(opt) == "our")
                method_dgdr = METHOD_dgdr::iour;
            else if(opts.at(opt) == "peer")
                method_dgdr = METHOD_dgdr::ipeer;
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                      << " given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " our " <<" peer "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_dgdr = method_dgdr;

        /// set parameters for lateral expanding
        opt = "method_spread";
        LatSpread::METHODS method_spread;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_spread = LatSpread::iNULL;
        }
        else{
            if(opts.at(opt) == "None")
                method_spread = LatSpread::iNULL;
            else if(opts.at(opt) == "AFGPY")
                method_spread = LatSpread::iAFGPY;
            else if(opts.at(opt) == "Adi")
                method_spread = LatSpread::iAdi;
            else if(opts.at(opt) == "AA")
                method_spread = LatSpread::iAA;
            else if(opts.at(opt) == "our")
                method_spread = LatSpread::iOUR;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " None " <<" AFGPY " << " Adi " << " AA " << " our "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_spread->setPars(a, theta_max,
                          id->theta_core,
                          id->theta_wing, method_spread);

        // set parameters for EOS
        opt = "method_eos";
        EOSadi::METHODS method_eos;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_eos = EOSadi::iNava13;
        }
        else{
            if(opts.at(opt) == "Nava13")
                method_eos = EOSadi::iNava13;
            else if(opts.at(opt) == "Peer12")
                method_eos = EOSadi::iPeer12;
            else{
                (*p_log)(LOG_ERR,AT)<< " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " Nava13 " << " Peer12 "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_eos->setPars(method_eos);


        /// set Nuclear Atomic pars
        p_nuc->setPars(pars, opts);

        /// limit the spreading
        opt = "method_limit_spread";
        METHOD_LIMIT_SPREAD method_limit_spread;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_limit_spread = METHOD_LIMIT_SPREAD::iNone;
        }
        else{
            if(opts.at(opt) == "None") {
                method_limit_spread = METHOD_LIMIT_SPREAD::iNone;
            }
            else if(opts.at(opt) == "Mom0Frac") {
                p_pars->fraction_of_mom0_when_spread_start =
                        (double)getDoublePar("mom0_frac_when_start_spread",pars, AT,p_log,.75, true);
                if (p_pars->fraction_of_mom0_when_spread_start <=0. or p_pars->fraction_of_mom0_when_spread_start > 1.){
                    (*p_log)(LOG_ERR,AT) << " Bad Value. mom0_frac_when_start_spread = "
                        << p_pars->fraction_of_mom0_when_spread_start
                        << "\n";
                    exit(1);
                }
                method_limit_spread = METHOD_LIMIT_SPREAD::iGamma0Frac;
            }
            else if(opts.at(opt) == "MomValAndFrac") {
                p_pars->value_of_mom_when_spread_start =
                        (double)getDoublePar("mom_when_start_spread",pars, AT,p_log,2., true);
                method_limit_spread = METHOD_LIMIT_SPREAD::iGammaVal;
            }
            else if(opts.at(opt) == "Rd") {
                method_limit_spread = METHOD_LIMIT_SPREAD::iRd;
            }
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " None "
                                      << " Rd "
                                      << " MomValAndFrac "
                                      << " Mom0Frac "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->method_limit_spread = method_limit_spread;

        /// EATS method guverns the BW discretization (piece-wise versus adaptive)
        p_pars->m_method_eats = id->method_eats;
        p_pars->nlayers = id->nlayers;
        p_pars->nshells = id->nshells;

        /// set boolean pars
        p_pars->allow_termination =
                getBoolOpt("allow_termination", opts, AT,p_log,false, true);

        p_pars->adiabLoss =
                getBoolOpt("use_adiabLoss", opts, AT,p_log,true, false);

        p_pars->do_mphys_in_ppr =
                getBoolOpt("do_mphys_in_ppr", opts, AT,p_log,true, false);
        p_pars->do_mphys_in_situ =
                getBoolOpt("do_mphys_in_situ", opts, AT,p_log,true, false);
        if (p_pars->do_mphys_in_ppr && p_pars->do_mphys_in_situ){
            (*p_log)(LOG_ERR,AT) << " only 'do_mphys_in_ppr' or 'do_mphys_in_situ' should be 'yes' "<<'\n';
            exit(1);
        }

        /// -----------  set initials and constants for the blast wave ----------------------
        size_t & ish = p_pars->ishell;
        size_t & il = p_pars->ilayer;

        p_pars->E0        = (double)id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_pw[ilayer];
        p_pars->Ye0       = (double)id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_pw[ilayer];
        p_pars->M0        = (double)id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_pw[ilayer];
        p_pars->R0        = (double)id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_pw[ilayer];
        p_pars->mom0      = (double)id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_pw[ilayer];
        p_pars->Gamma0    = EQS::GammaFromMom(id->get(ish,il,EjectaID2::Q::imom));//latStruct.dist_Mom0_pw[ilayer];
        p_pars->beta0     = EQS::BetaFromMom(id->get(ish,il,EjectaID2::Q::imom));
        p_pars->s0        = (double)id->get(ish,il,EjectaID2::Q::ientr);//latStruct.dist_s_pw[ilayer];
        p_pars->tb0       = m_tb_arr.empty() ? 0 : m_tb_arr[0];
        p_pars->theta_a   = 0.;
        p_pars->theta_b0  = ((id->method_eats) == EjectaID2::ipiecewise)
                            ? id->theta_wing : id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.m_theta_w; // theta_b0
        p_pars->ctheta0   = id->get(ish,il,EjectaID2::Q::ictheta); // TODO 0.5 * (latStruct.thetas_c_l[ilayer] + latStruct.thetas_c_h[ilayer]);
        p_pars->theta_w   = id->theta_wing;//latStruct.m_theta_w; //
        p_pars->theta_max = theta_max;
        p_pars->ncells    = ((id->method_eats) == EjectaID2::ipiecewise)
                            ? (double)id->ncells : 1.;//(double) latStruct.ncells;

        /// quick check
        if (p_pars->beta0 < p_pars->min_beta_0){
            (*p_log)(LOG_ERR,AT) << "ish="<<ish<<" il="<<il
                <<" Gamma0="<<p_pars->Gamma0<<" beta0="<<p_pars->beta0<<" < min="<<p_pars->min_beta_0<<"\n";
            exit(1);
        }
        if (p_pars->E0 < 1e10){
            (*p_log)(LOG_ERR,AT) << "ish="<<ish<<" il="<<il<<" E0="<<p_pars->E0<<" < min="<<1e10<<"\n";
            exit(1);
        }

        p_spread->m_theta_b0 = p_pars->theta_b0;
        p_pars->prev_x = p_pars->tb0;
        p_pars->prev_idx_x = 0;
        p_pars->ii_eq  = ii_eq;

        /// how to compute obs.flux via interpolation of comoving intensity, or on the fly
        opt = "method_comp_mode";
        METHODS_RAD methodCompMode;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodCompMode = METHODS_RAD::iobservflux;
        }
        else{
            if(opts.at(opt) == "observFlux")
                methodCompMode = METHODS_RAD::iobservflux;
            else if(opts.at(opt) == "comovSpec")
                methodCompMode = METHODS_RAD::icomovspec;
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " observFlux " << " comovSpec " << "\n";
                exit(1);
            }
        }
        p_pars->m_method_rad = methodCompMode;

        if (p_pars->m_type == BW_TYPES::iFSRS)
            setParamsRS(pars,opts);

        if (p_pars->m_type == BW_TYPES::iFS_DENSE)
            setParamsDense(pars,opts);

        if (p_pars->m_type == BW_TYPES::iFS_PWN_DENSE) {
            setParamsDense(pars, opts);
            setParamsPWN(pars, opts);
        }
    }

    void setParamsRS(StrDbMap & pars, StrStrMap & opts){

        p_pars->rs_Gamma0_frac_no_exceed =
                (double)getDoublePar("rs_Gamma0_frac_no_exceed",pars, AT,p_log,.98,false);

        /// rs parameters
        p_pars->tprompt = getDoublePar("tprompt",pars,AT,p_log,1000.,false);
        p_pars->epsilon_rad_rs = getDoublePar("epsilon_e_rad_rs",pars,AT,p_log,.0,false);
        p_pars->rs_shutOff_criterion_rho = getDoublePar("rs_shutOff_criterion_rho",pars,AT,p_log,1.e-50,false);
        p_pars->init_deltaR4 = getBoolOpt("init_deltaR4",opts,AT,p_log,false, true);
        p_pars->exponential_rho4 = getBoolOpt("exponential_rho4",opts,AT,p_log, true, true);
//        /// RS: method for shock velocity
//        std::string opt = "method_GammaRsh";
//        METHOD_Gamma_sh m_method_gamma_rsh;
//        if ( opts.find(opt) == opts.end() ) {
//            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
//            m_method_gamma_rsh = METHOD_Gamma_sh::iuseGammaShock;
//        }
//        else{
//            if(opts.at(opt) == "useJustGamma")
//                m_method_gamma_rsh = METHOD_Gamma_sh::iuseJustGamma;
//            else if(opts.at(opt) == "useGammaShock")
//                m_method_gamma_rsh = METHOD_Gamma_sh::iuseGammaShock;
//            else{
//                (*p_log)(LOG_ERR,AT) << " option for: " << opt
//                                     <<" given: " << opts.at(opt)
//                                     << " is not recognized "
//                                     << " Possible options: "
//                                     << " useGammaShock " << " useJustGamma "
//                                     << " Exiting...\n";
////                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        p_pars->m_method_gamma_rsh = m_method_gamma_rsh;

//        opt = "method_shock_vel_rs";
//        METHODS_SHOCK_VEL methodsShockVel_rs;
//        if ( opts.find(opt) == opts.end() ) {
//            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
//            methodsShockVel_rs = METHODS_SHOCK_VEL::isameAsBW;
//        }
//        else{
//            if(opts.at(opt) == "sameAsBW")
//                methodsShockVel_rs = METHODS_SHOCK_VEL::isameAsBW;
//            else if(opts.at(opt) == "shockVel")
//                methodsShockVel_rs = METHODS_SHOCK_VEL::ishockVel;
//            else{
//                (*p_log)(LOG_ERR,AT) << " option for: " << opt
//                                     <<" given: " << opts.at(opt)
//                                     << " is not recognized. "
//                                     << " Possible options: "
//                                     << " sameAsBW " << " shockVel " << "\n";
////                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        p_pars->method_shock_vel_rs = methodsShockVel_rs;


        p_pars->do_rs = getBoolOpt("do_rs", opts, AT,p_log, false, true);
        if (p_pars->do_rs && p_pars->m_type!=BW_TYPES::iFSRS){
            (*p_log)(LOG_ERR,AT)<<" if do_rs = yes, the rhs_type should be 'fsrs' \n";
            exit(1);
        }
        p_pars->do_rs_radiation = getBoolOpt("do_rs_radiation", opts, AT,p_log, false, true);
        if (p_pars->do_rs_radiation && p_pars->m_type!=BW_TYPES::iFSRS){
            (*p_log)(LOG_ERR,AT)<<" if do_rs_radiation = yes, the rhs_type should be 'fsrs' \n";
            exit(1);
        }
        if (p_pars->do_rs_radiation and not p_pars->do_rs){
            (*p_log)(LOG_ERR,AT)<<" if do_rs_radiation = yes, the do_rs should be 'yes' \n";
            exit(1);
        }

        p_pars->adiabLoss_rs =
                getBoolOpt("use_adiabLoss_rs", opts, AT,p_log,true, false);

        if (p_pars->do_rs) {
            p_pars->min_Gamma0_for_rs =
                    getDoublePar("min_Gamma0_for_rs", pars, AT, p_log, 0., true);
        }
        else
            p_pars->min_Gamma0_for_rs = 0;
        if (p_pars->Gamma0 < p_pars->min_Gamma0_for_rs) {
            p_pars->do_rs = false;
            p_pars->do_rs_radiation = false;
            p_pars->shutOff = true;
            p_pars->m_type = BW_TYPES::iFS;
        }

    }

    void setParamsDense(StrDbMap & pars, StrStrMap & opts){

        // interaction parameters
        p_pars->which_jet_layer_to_use =
                (int)getDoublePar("which_jet_layer_to_use",pars,AT,p_log,1.e5,false);
        p_pars->steepnes_of_exp_decay =
                (double)getDoublePar("steepnes_of_exp_decay",pars,AT,p_log,1.,false);
        p_pars->Gamma_when_st_starts =
                (double)getDoublePar("Gamma_when_st_starts",pars,AT,p_log,2.,false);
        p_pars->fraction_of_Gamma0_when_bm_for_bm =
                (double)getDoublePar("fraction_of_Gamma0_when_bm_for_bm",pars, AT,p_log,1.98,false);


        p_pars->use_dens_prof_behind_jet_for_ejecta =
                getBoolOpt("use_dens_prof_behind_jet_for_ejecta", opts, AT,p_log,false, false);

        p_pars->use_dens_prof_inside_ejecta =
                getBoolOpt("use_dens_prof_inside_ejecta", opts, AT,p_log,false, false);

        p_pars->use_dens_prof_inside_ejecta =
                getBoolOpt("use_dens_prof_inside_ejecta", opts, AT,p_log,false, false);

        p_pars->use_exp_rho_decay_as_floor =
                getBoolOpt("use_exp_rho_decay_as_floor", opts, AT,p_log,false, false);

        p_pars->use_flat_dens_floor =
                getBoolOpt("use_flat_dens_floor", opts, AT,p_log,false, false);

        p_pars->use_st_dens_profile =
                getBoolOpt("use_st_dens_profile", opts, AT,p_log,false, false);

        p_pars->use_bm_dens_profile =
                getBoolOpt("use_bm_dens_profile", opts, AT,p_log,false, false);

        /// set sedov-taylor profile (for jet to be seen by ejecta as it moves behind)
        if (p_pars->use_st_dens_profile) {
            p_sedov->setPars(1.5, 3, 0., 1e5); // TODO this should not be here and evaluated for EVERY bw...
            p_sedov->evaluate();
        }
    }

    void setParamsPWN(StrDbMap & pars, StrStrMap & opts){
        /// if evolve BW with energy injection; method to eval. its vol;rho;thick;tau
        std::string opt = "method_thick_for_rho";
        METHOD_THICK_FOR_RHO m_method_thick_for_rho;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_thick_for_rho = METHOD_THICK_FOR_RHO::iFromRadius;
        }
        else{
            if(opts.at(opt) == "bw_radius")
                m_method_thick_for_rho = METHOD_THICK_FOR_RHO::iFromRadius;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " bw_radius "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->method_thick_for_rho = m_method_thick_for_rho;

        //        int p_pars->opacitymode=0; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
        if (p_pars->Ye0 <= 0.2)
            p_pars->opacitymode = 2;
        else if ((p_pars->Ye0 > 0.2) or (p_pars->Ye0 < 0.3))
            p_pars->opacitymode = 1;
        else if (p_pars->Ye0 >= 0.3)
            p_pars->opacitymode = 0;
        else if (p_pars->Ye0 > 0.5)
            p_pars->opacitymode = 3;
        else{
            (*p_log)(LOG_ERR,AT) << " error \n";
            exit(1);
        }

        if(p_pars->opacitymode==0){
            p_pars->Z_eff = 24.21; /* effective nuclear weight */
            p_pars->mu_e = 2.148;
        }
        else if(p_pars->opacitymode==1){
            p_pars->Z_eff = 26.74;
            p_pars->mu_e = 2.2353;
        }
        else if(p_pars->opacitymode==2){
            p_pars->Z_eff = 53.90;
            p_pars->mu_e= 2.4622;
        }
        else if(p_pars->opacitymode==3){
            p_pars->Z_eff = 7.0;
            p_pars->mu_e = 2.0;
        }

        /// ionization fraction in ejecta
        p_pars->albd_fac = getDoublePar("albd_fac", pars, AT, p_log,0.5, false);//pars.at("freq1");
        p_pars->eps_e_w = getDoublePar("eps_e_w",pars,AT,p_log,-1,false); // electron acceleration efficiency
        p_pars->eps_mag_w = getDoublePar("eps_mag_w",pars,AT,p_log,-1,false); // magnetic field efficiency
        p_pars->epsth_w = getDoublePar("eps_th_w",pars,AT,p_log,-1,false); // initial absorption fraction
        p_pars->gamma_b_w = getDoublePar("gamma_b_w",pars,AT,p_log,-1,false); //break Lorentz factor of electron injection spectrum

        /// for shell structured, if shells collide, how to treat last shell properties
        opt = "method_single_bw_delta";
        METHOD_SINGLE_BW_DELTA method_single_bw_delta;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_single_bw_delta = METHOD_SINGLE_BW_DELTA::iconst;
        }
        else{
            if(opts.at(opt) == "frac_last") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ifrac_last;
//                p_pars->fraction_last_delta =
//                        (double)getDoublePar("fraction_last_delta",pars, AT,p_log,.75,true);
            }
            else if(opts.at(opt) == "const")
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::iconst;
            else if(opts.at(opt) == "lr")
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ilr;
            else if(opts.at(opt) == "bm") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ibm;
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
            }
            else if(opts.at(opt) == "st") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ist;
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
            }
            else if(opts.at(opt) == "bm_st") {
                method_single_bw_delta = METHOD_SINGLE_BW_DELTA::ibm_st; // TODO add an option that smoothly connects BM and ST profiles
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
                p_pars->mom_when_bm_to_st =
                        (double)getDoublePar("mom_when_bm_to_st",pars, AT,p_log,.75,true);
            }
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " frac_last " << " bm " << " st "<< " bm "<< " bm_st "<< " bm "<< "\n";
                exit(1);
            }
        }
        p_pars->method_single_bw_delta = method_single_bw_delta;

        if (p_pars->use_dens_prof_inside_ejecta) {
            m_data_shells.resize(BW::NVALS_SH);
            for (auto &arr: m_data_shells)
                arr.resize(p_pars->nshells);
        }

        frac_psr_dep_.resize(iters);


    }

    /// --------------------------------------------------------
    Pars *& getPars(){ return p_pars; }
    std::unique_ptr<EOSadi> & getEos(){ return p_eos; }
    std::unique_ptr<LatSpread> & getSpread(){ return p_spread; }
    std::unique_ptr<RhoISM> & getDensIsm(){ return p_dens; }
    std::unique_ptr<SedovTaylor> & getSedov(){ return p_sedov; }
    std::unique_ptr<BlandfordMcKee2> & getBM(){ return p_bm; }
    std::unique_ptr<LinearRegression> & getLRforDelta(){ return p_lr_delta; }
    std::unique_ptr<LinearRegression> & getLRforVol(){ return p_lr_vol; }
    /// --------------------------------------------------------
    inline Vector & operator[](unsigned ll){ return this->m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return this->m_data[ivn][ir]; }
    inline double ctheta(double theta){
        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        if (theta > p_pars->theta_max ){
//            std::cerr << AT << " theta="<<theta<<" > theta_max=" << p_pars->theta_max << "\n";
//        }
//        if (std::fabs( theta - p_pars->theta_b0) > 1e-2){
//            std::cerr << AT << " theta="<<theta<<" < theta_b0=" << p_pars->theta_b0 <<"\n";
//            exit(1);
//        }
//        double ctheta = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w); // TODO WROOOONG

        double ctheta = 0.;
        if (p_pars->ilayer > 0) {
            //
            double fac0 = (double)p_pars->ilayer/(double)p_pars->nlayers;
            double fac1 = (double)(p_pars->ilayer+1)/(double)p_pars->nlayers;
//            std::cout << std::asin(CGS::pi*3/4.) << "\n";
            if (!std::isfinite(std::sin(theta))){
                (*p_log)(LOG_ERR,AT) << " sin(theta= "<<theta<<") is not finite... Exiting..." << "\n";
                exit(1);
            }

            double x2 = fac1*std::sin(theta / 2.);
            double xx2 = 2.*std::asin(x2);
            double x1 = fac0*std::sin(theta / 2.);
            double xx1 = 2.*std::asin(x1);

            ctheta = 0.5 * (xx1 + xx2);
            if (!std::isfinite(ctheta)){
                (*p_log)(LOG_ERR,AT) << "ctheta is not finite. ctheta="<<ctheta<<" Exiting..." << "\n";
                exit(1);
            }
        }
        return ctheta;
    }
    inline VecVector & getData(){ return m_data; }
    inline Vector & getData(BW::Q var){ return m_data[ var ]; }
    inline double & getData(BW::Q var, size_t i){ return m_data[ var ][i]; }
    inline Vector & get_tburst(){return m_tb_arr;}
    ~BlastWaveBase(){ delete p_pars; }
};


#endif //SRC_BLASTWAVE_BASE_H
