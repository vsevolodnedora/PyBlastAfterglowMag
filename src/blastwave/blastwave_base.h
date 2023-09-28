//
// Created by vsevolod on 14/07/23.
//

#ifndef SRC_BLASTWAVE_BASE_H
#define SRC_BLASTWAVE_BASE_H

#include "../composition.h"
//#include "blastwave_components.h"
#include "blastwave_pars.h"

class BlastWaveBase{
    std::unique_ptr<logger> p_log;
protected:
    Vector m_tb_arr;
    VecVector mD{}; // container for the solution of the evolution
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
    static constexpr int iters=1000;
    Vector frac_psr_dep_{};
public:
    bool is_initialized = false;
    BlastWaveBase(Vector & tb_arr, size_t ishell, size_t ilayer, size_t n_substeps, BW_TYPES type, int loglevel)
        : m_tb_arr(tb_arr){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "BW_Base");

        /// parameters
        p_pars = new Pars(mD, mDtmp, loglevel); //
        p_pars->m_type = type;

        /// the container for the final solution
        if (mD.empty()){
            mD.resize( BW::TOTAL_VARS );
        }
        /// the container for the last N substeps
        if (mDtmp.empty()){
            mDtmp.resize( BW::TOTAL_VARS );
        }
        /// Check if contener will be filled by evolving or loading
        if (m_tb_arr.empty()){
            (*p_log)(LOG_WARN,AT) << " Time grid was not initialized\n";
        }
        /// if no evolution required; do not allocate memory for each variable (only for evolved)
        if (mD[BW::Q::itburst].size() < 1)
            for (auto & ivn : BW::VARS.at(p_pars->m_type))
                mD[ivn].resize(tb_arr.size(), 0.0);

        // ---------------------- Methods
        p_lr_delta = std::make_unique<LinearRegression>(mD[BW::Q::iR], mD[BW::Q::iEJdelta]);
        p_lr_vol = std::make_unique<LinearRegression>(mD[BW::Q::iR], mD[BW::Q::iEJvol]);
        p_spread = std::make_unique<LatSpread>();
        p_eos = std::make_unique<EOSadi>();
        p_dens = std::make_unique<RhoISM>(loglevel);
        p_sedov = std::make_unique<SedovTaylor>();
        p_bm = std::make_unique< BlandfordMcKee2>();
        p_nuc = std::make_unique<NuclearAtomic>(loglevel);
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

        double nism, A0, s, r_ej, r_ism,  a, theta_max, epsilon_e_rad;

        // set parameters for ISM density
        nism = getDoublePar("n_ism", pars, AT, p_log, -1, true);//pars.at("nism");
        A0 = getDoublePar("A0", pars, AT,p_log,-1,false);//pars.at("A0");
        s = getDoublePar("s", pars, AT,p_log,-1,false);//pars.at("s");
        r_ej = getDoublePar("r_ej", pars, AT,p_log,-1,false);//pars.at("r_ej");
        r_ism = getDoublePar("r_ism", pars, AT,p_log,-1,false);//pars.at("r_ism");
        p_dens->setPars(nism, A0, s, r_ej, r_ism, true);

        // spreading
        a = getDoublePar("a", pars, AT,p_log,-1,false);//pars.at("a");
        theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);//pars.at("theta_max");

        // radiative losses
        epsilon_e_rad = getDoublePar("epsilon_e_rad", pars, AT,p_log,0.,false);// pars.at("epsilon_e_rad");

        // interaction parameters
        p_pars->which_jet_layer_to_use =
                (int)getDoublePar("which_jet_layer_to_use",pars,AT,p_log,1.e5,false);
        p_pars->steepnes_of_exp_decay =
                (double)getDoublePar("steepnes_of_exp_decay",pars,AT,p_log,1.,false);
        p_pars->Gamma_when_st_starts =
                (double)getDoublePar("Gamma_when_st_starts",pars,AT,p_log,2.,false);
        p_pars->fraction_of_Gamma0_when_bm_for_bm =
                (double)getDoublePar("fraction_of_Gamma0_when_bm_for_bm",pars, AT,p_log,1.98,false);
        p_pars->rs_Gamma0_frac_no_exceed =
                (double)getDoublePar("rs_Gamma0_frac_no_exceed",pars, AT,p_log,.98,false);

        /// rs parameters
        p_pars->tprompt = getDoublePar("tprompt",pars,AT,p_log,1000.,false);
        p_pars->epsilon_rad_rs = getDoublePar("epsilon_rad_rs",pars,AT,p_log,.0,false);
        p_pars->rs_shutOff_criterion_rho = getDoublePar("rs_shutOff_criterion_rho",pars,AT,p_log,1.e-50,false);

        // set options
        std::string opt;
        // mass accretion from ISM
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
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
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
//        p_pars->only_last_shell_dmdr
//            = getBoolOpt("only_last_shell_dmdr", opts, AT,p_log,false, false);


        // evolution eq
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
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
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

        // set parameters for lateral expanding
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
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " None " <<" AFGPY " << " Adi " << " AA "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_spread->setPars(a,theta_max,id->theta_core,
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
                (*p_log)(LOG_WARN,AT)<< " option for: " << opt
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

        // set Nuclear Atomic pars
        p_nuc->setPars(pars, opts);

        /// method for shock radius
        opt = "method_Rsh";
        METHOD_RSh m_method_r_sh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_r_sh = METHOD_RSh::isameAsR;
        }
        else{
            if(opts.at(opt) == "sameAsR")
                m_method_r_sh = METHOD_RSh::isameAsR;
            else if(opts.at(opt) == "useGammaSh")
                m_method_r_sh = METHOD_RSh::iuseGammaSh;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " sameAsR " << " useGammaSh "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_r_sh = m_method_r_sh;

        /// method for shock velocity
        opt = "method_GammaSh";
        METHOD_GammaSh m_method_gamma_sh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_gamma_sh = METHOD_GammaSh::iuseGammaShock;
        }
        else{
            if(opts.at(opt) == "useJustGamma")
                m_method_gamma_sh = METHOD_GammaSh::iuseJustGamma;
            else if(opts.at(opt) == "useGammaShock")
                m_method_gamma_sh = METHOD_GammaSh::iuseGammaShock;
            else if(opts.at(opt) == "useJustGammaRel") {
                if (!(p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE)){
                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useJustGammaRel "<< " if bw_type ="<<p_pars->m_type<<"\n";
                    exit(1);
                }
                m_method_gamma_sh = METHOD_GammaSh::iuseJustGammaRel;
            }
            else if(opts.at(opt) == "useGammaRelShock") {
                if (!(p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE)){
                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useGammaRelShock "<< " if bw_type ="<<p_pars->m_type<<"\n";
                    exit(1);
                }
                m_method_gamma_sh = METHOD_GammaSh::iuseGammaRelShock;
            }
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " useGammaShock " << " useJustGammaRel " << " useJustGammaRel "<< " useGammaRelShock "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_gamma_sh = m_method_gamma_sh;

        /// method for energy density behind shock
        opt = "method_Up";
        METHODS_Up m_method_up;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_up = METHODS_Up::iuseEint2;
        }
        else{
            if(opts.at(opt) == "useEint2")
                m_method_up = METHODS_Up::iuseEint2;
            else if(opts.at(opt) == "useGamma")
                m_method_up = METHODS_Up::iuseGamma;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " useEint2 " << " useGamma "
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_up = m_method_up;

        /// method for shock thickness
        opt = "method_Delta";
        METHOD_Delta m_method_delta;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_delta = METHOD_Delta::iuseJoh06;
        }
        else{
            if(opts.at(opt) == "useJoh06")
                m_method_delta = METHOD_Delta::iuseJoh06;
            else if(opts.at(opt) == "useVE12")
                m_method_delta = METHOD_Delta::iuseVE12;
            else if(opts.at(opt) == "None")
                m_method_delta = METHOD_Delta::iNoDelta;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " useJoh06 " << " useVE12 " << "None"
                                      << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_Delta = m_method_delta;

        /// if evolve BW with energy injection; method to eval. its vol;rho;thick;tau
        opt = "method_thick_for_rho";
        METHOD_THICK_FOR_RHO m_method_thick_for_rho;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            m_method_thick_for_rho = METHOD_THICK_FOR_RHO::iFromRadius;
        }
        else{
            if(opts.at(opt) == "bw_radius")
                m_method_thick_for_rho = METHOD_THICK_FOR_RHO::iFromRadius;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
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
            else if(opts.at(opt) == "Gamma0Frac") {
                p_pars->fraction_of_mom0_when_spread_start =
                        (double)getDoublePar("mom0_frac_when_start_spread",pars, AT,p_log,.75, true);
                method_limit_spread = METHOD_LIMIT_SPREAD::iGamma0Frac;
            }
            else if(opts.at(opt) == "GammaVal") {
                p_pars->value_of_mom_when_spread_start =
                        (double)getDoublePar("mom_when_start_spread",pars, AT,p_log,.75, true);
                method_limit_spread = METHOD_LIMIT_SPREAD::iGammaVal;
            }
            else if(opts.at(opt) == "Rd") {
                method_limit_spread = METHOD_LIMIT_SPREAD::iRd;
            }
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized "
                                      << " Possible options: "
                                      << " None "
                                      << " Rd "
                                      << " GammaVal "
                                      << " Gamma0Frac "
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

        p_pars->adiabLoss =
                getBoolOpt("use_adiabLoss", opts, AT,p_log,true, false);



        /// set sedov-taylor profile (for jet to be seen by ejecta as it moves behind)
        if (p_pars->use_st_dens_profile) {
            p_sedov->setPars(1.5, 3, 0., 1e5); // TODO this should not be here and evaluated for EVERY bw...
            p_sedov->evaluate();
        }

        /// set parameters for computing observed emission

        /// -----------  set initials and constants for the blast wave ----------------------
        size_t & ish = p_pars->ishell;
        size_t & il = p_pars->ilayer;

        p_pars->E0        = (double)id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_pw[ilayer];
        p_pars->Ye0       = (double)id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_pw[ilayer];
        p_pars->M0        = (double)id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_pw[ilayer];
        p_pars->R0        = (double)id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_pw[ilayer];
        p_pars->mom0      = (double)id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_pw[ilayer];
        p_pars->Gamma0    = GamFromMom(id->get(ish,il,EjectaID2::Q::imom));//latStruct.dist_Mom0_pw[ilayer];
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
        p_pars->eps_rad   = epsilon_e_rad;

        p_spread->m_theta_b0 = p_pars->theta_b0;
        p_pars->prev_x = p_pars->tb0;
        p_pars->prev_idx_x = 0;
        p_pars->ii_eq  = ii_eq;
#if 0
        switch (id->method_eats) {

            case EjectaID2::ipiecewise:
#if 0
                if ( latStruct.dist_E0_pw.empty() || latStruct.dist_M0_pw.empty() || latStruct.dist_Mom0_pw.empty()
                     || latStruct.dist_Ye_pw.empty() || latStruct.theta_pw.empty()){
                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
                    exit(1);
                }
                (*p_log)(LOG_INFO,AT) << " Init. [pw] "
                                      << " E0="<<latStruct.dist_E0_pw[ilayer]
                                      << " Ye="<<latStruct.dist_Ye_pw[ilayer]
                                      << " s="<<latStruct.dist_s_pw[ilayer]
                                      << " M0="<<latStruct.dist_M0_pw[ilayer]
                                      << " M0m0="<<latStruct.dist_Mom0_pw[ilayer]
                                      << " beta0="<<EQS::BetFromMom(latStruct.dist_Mom0_pw[ilayer])
                                      << " tb0="<<m_tb_arr[0]
                                      << " thetab0="<<latStruct.m_theta_w
                                      << " theta0="<<latStruct.theta_pw[ilayer]
                                      << " theta1="<<latStruct.theta_pw[ilayer]
                                      << " theta_w="<<latStruct.m_theta_w
                                      << " ii_eq="<<ii_eq
                                      << " ncells="<<latStruct.ncells
                                      << "\n";

                if (latStruct.m_theta_w > theta_max){
                    (*p_log)(LOG_ERR,AT) << " theta_b0="<<latStruct.m_theta_w<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
#endif

                p_pars->E0        = id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_pw[ilayer];
                p_pars->Ye0       = id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_pw[ilayer];
                p_pars->s0        = id->get(ish,il,EjectaID2::Q::is);//latStruct.dist_s_pw[ilayer];
                p_pars->M0        = id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_pw[ilayer];
                p_pars->R0        = id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_pw[ilayer];
                p_pars->mom0      = id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_pw[ilayer];
//                p_pars->Eint0     = id->get(ish,il,EjectaID2::Q::ieint);
                p_pars->tb0       = m_tb_arr.empty() ? 0 : m_tb_arr[0];
                p_pars->theta_a   = 0.; // theta_a
                p_pars->theta_b0  = id->theta_wing;//latStruct.m_theta_w; // theta_b0
                p_pars->ctheta0   = id->get(ish,il,EjectaID2::Q::ictheta); //TODO !! 0.5 * (latStruct.theta_pw[ilayer] + latStruct.theta_pw[ilayer]);
//        p_pars->theta_h0= theta_c_h;
//                p_pars->theta_c_l = 0.;//id->get(ish,il,EjectaID2::Q::itheta_c_l);//latStruct.theta_pw[ilayer];//theta_c_l;
//                p_pars->theta_c_h = 0.,//id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.theta_pw[ilayer];
                p_pars->theta_w   = id->theta_wing;//latStruct.m_theta_w; //
                p_pars->theta_max = theta_max;
                p_pars->ncells    = (double)id->ncells;//(double) latStruct.ncells;
                p_pars->eps_rad   = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->prev_x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;

                /// initialize non-thermal radiation module
//                p_rad->setEatsPars(pars,opts,id->m_nlayers,
//                                   id->get(ish,il,EjectaID2::Q::ictheta),
//                                   0.,0.,id->theta_wing,
//                                   getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                break;

            case EjectaID2::iadaptive:
#if 0
                if ( latStruct.dist_E0_a.empty() || latStruct.dist_M0_a.empty() || latStruct.dist_Mom0_a.empty()
                     || latStruct.dist_Ye_a.empty() || latStruct.thetas_c_h.empty()){
                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
                    exit(1);
                }
                (*p_log)(LOG_INFO,AT)<<"Init. [a] "
                                     << " E0="<<latStruct.dist_E0_a[ilayer]
                                     << " Ye="<<latStruct.dist_Ye_a[ilayer]
                                     << " s="<<latStruct.dist_s_a[ilayer]
                                     << " M0="<<latStruct.dist_M0_a[ilayer]
                                     << " G0="<<latStruct.dist_Mom0_a[ilayer]
                                     << " beta0="<<EQS::BetFromMom(latStruct.dist_Mom0_a[ilayer])
                                     << " tb0="<<m_tb_arr[0]
                                     << " thetab0="<<latStruct.thetas_c_h[ilayer]
                                     << " theta0="<<latStruct.thetas_c_l[ilayer]
                                     << " theta1="<<latStruct.thetas_c_h[ilayer]
                                     << " theta_w="<<latStruct.m_theta_w
                                     << " ii_eq="<<ii_eq
                                     << " ncells="<<1.
                                     << "\n";
                double fac = 2 * std::sin(0.5 * latStruct.thetas_c_h[ilayer]) * std::sin(0.5 * latStruct.thetas_c_h[ilayer]);
                if (latStruct.thetas_c_h[ilayer] > theta_max){
                    (*p_log)(LOG_ERR,AT) << " theta_b0="<<latStruct.thetas_c_h[ilayer]<<" exceeds theta_max="<<theta_max<<" \n Exiting... \n";
//                    std::cerr << AT << "\n";
                    exit(1);
                }
#endif
                p_pars->E0      = id->get(ish,il,EjectaID2::Q::iek);//latStruct.dist_E0_a[ilayer];
                p_pars->Ye0     = id->get(ish,il,EjectaID2::Q::iye);//latStruct.dist_Ye_a[ilayer];
                p_pars->s0      = id->get(ish,il,EjectaID2::Q::is);//latStruct.dist_s_a[ilayer];
                p_pars->M0      = id->get(ish,il,EjectaID2::Q::imass);//latStruct.dist_M0_a[ilayer];
                p_pars->R0      = id->get(ish,il,EjectaID2::Q::ir);//latStruct.dist_M0_a[ilayer];
                p_pars->mom0    = id->get(ish,il,EjectaID2::Q::imom);//latStruct.dist_Mom0_a[ilayer];
                p_pars->tb0     = m_tb_arr.empty() ? 0 : m_tb_arr[0];
                p_pars->theta_a = 0.;
                p_pars->theta_b0= id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.thetas_c_h[ilayer];
                p_pars->ctheta0 = id->get(ish,il,EjectaID2::Q::ictheta); // TODO 0.5 * (latStruct.thetas_c_l[ilayer] + latStruct.thetas_c_h[ilayer]);
//        p_pars->theta_h0= theta_c_h;
//                p_pars->theta_c_l = id->get(ish,il,EjectaID2::Q::itheta_c_l);//latStruct.thetas_c_l[ilayer];
//                p_pars->theta_c_h = id->get(ish,il,EjectaID2::Q::itheta_c_h);//latStruct.thetas_c_h[ilayer];
                p_pars->theta_w = 0.; //
                p_pars->theta_max = theta_max;
                p_pars->ncells  = 1.;
                p_pars->eps_rad = epsilon_e_rad;

                p_spread->m_theta_b0 = p_pars->theta_b0;
                p_pars->prev_x = p_pars->tb0;

                p_pars->ii_eq  = ii_eq;

                /// initialize non-thermal radiation module
//                p_rad->setEatsPars(pars,opts,id->m_nlayers,id->get(ish,il,EjectaID2::Q::ictheta),
//                                   id->get(ish,il,EjectaID2::Q::itheta_c_l),
//                                   id->get(ish,il,EjectaID2::Q::itheta_c_h),0.,
//                                   getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false));

                // double E0, double M0, double Gamma0, double tb0, double theta_a, double theta_b0,
                // double theta_c_l, double theta_c_h, double theta_w, double theta_max, double epsilon_e_rad,
                // size_t ii_eq,
                // double ncells
//                bw_obj.setMagPars(latStruct.dist_E0_a[ilayer], = double E0,
//                               latStruct.dist_M0_a[ilayer], = double M0,
//                               latStruct.dist_G0_a[ilayer], = double Gamma0,
//                               t_grid[0],                   = double tb0,
//                               0.,                          = double theta_a,
//                               latStruct.thetas_c_h[ilayer],= double theta_b0,
//                               latStruct.thetas_c_l[ilayer],= double theta_c_l,
//                               latStruct.thetas_c_h[ilayer],= double theta_c_h,
//                               0.,                          = double theta_w,
//                               theta_max,                   = double theta_max,
////                             latStruct.cthetas0[ilayer],  = double epsilon_e_rad,
//                               epsilon_e_rad,               = size_t ii_eq,
//                               ii_eq,
//                               1.);
                break;
        }
#endif

        /// ----------------------- set options ------------------------------


        opt = "method_ne";
        METHOD_NE methodNe;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodNe = METHOD_NE::iuseNe;
        }
        else{
            if(opts.at(opt) == "useNe")
                methodNe = METHOD_NE::iuseNe;
            else if(opts.at(opt) == "usenprime")
                methodNe = METHOD_NE::iusenprime;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useNe " << " usenprime " << "\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->m_method_ne = methodNe;

        opt = "method_shock_vel";
        METHODS_SHOCK_VEL methodsShockVel;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
        }
        else{
            if(opts.at(opt) == "sameAsBW")
                methodsShockVel = METHODS_SHOCK_VEL::isameAsBW;
            else if(opts.at(opt) == "shockVel")
                methodsShockVel = METHODS_SHOCK_VEL::ishockVel;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << " Possible options: "
                                     << " sameAsBW " << " shockVel " << "\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->method_shock_vel = methodsShockVel;

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

        /// -------------------------------------

        p_pars->do_rs = getBoolOpt("do_rs", opts, AT,p_log, false, true);
        if (p_pars->do_rs && p_pars->m_type!=BW_TYPES::iFSRS){
            (*p_log)(LOG_ERR,AT)<<" if do_rs = yes, the rhs_type should be 'fsrs' \n";
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
            p_pars->shutOff = true;
            p_pars->m_type = BW_TYPES::iFS;
        }

        /// ---------------------------------------
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

        /// PWN wind initial radius
//        p_pars->pwn_Rw0 = getDoublePar("Rw0", pars, AT, p_log, -1, false); // PWN radius at t=0; [km]



//        std::cout << " ["<<bw_obj.getPars()->ishell<<", "<<bw_obj.getPars()->ilayer<<"] "
//                  <<" G0="<<bw_obj.getPars()->Gamma0
//                  <<" E0="<<bw_obj.getPars()->E0
//                  <<" M0="<<bw_obj.getPars()->M0
//                  <<" ctheta0="<<bw_obj.getPars()->ctheta0
//                  <<"\n";
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
//    size_t ntb() const { return m_tb_arr.size(); }
//    Vector & getTbGrid() {return m_tb_arr;}
//    Vector getTbGrid(size_t every_it) {
//        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
//        Vector tmp{};
//        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
//            tmp.push_back(m_tb_arr[it]);
//        }
////        Vector tmp2 (tmp.data(), tmp.size());
//        return std::move(tmp);
//    }
    inline Vector & operator[](unsigned ll){ return this->mD[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return this->mD[ivn][ir]; }
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
    inline VecVector & getData(){ return mD; }
    inline Vector & getData(BW::Q var){ return mD[ var ]; }
    inline Vector & get_tburst(){return m_tb_arr;}
    ~BlastWaveBase(){ delete p_pars; }
};


#endif //SRC_BLASTWAVE_BASE_H
