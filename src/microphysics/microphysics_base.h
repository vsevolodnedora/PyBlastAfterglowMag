//
// Created by vsevolod on 4/29/24.
//

#ifndef SRC_MICROPHYSICS_BASE_H
#define SRC_MICROPHYSICS_BASE_H

#include "../utilitites/pch.h"
#include "../utilitites/logger.h"
#include "../utilitites/utils.h"
#include "../utilitites/interpolators.h"
#include "../utilitites/quadratures.h"
#include "../utilitites/rootfinders.h"

#include "../blastwave/blastwave_components.h"

#include "analytic_models.h"
#include "numeric_model.h"

enum METHODS_Up_sh { iuseEint2, iuseGamma }; // energy density behind the shock

enum METHOD_thickness_sh { iuseJoh06, iuseVE12, iNoDelta }; // thickness of the shock

enum METHOD_Gamma_sh { iuseGammaShock, iuseJustGamma, iuseJustGammaRel, iuseGammaRelShock };

enum METHOD_R_sh { isameAsR, iuseGammaSh };

enum METHODS_SHOCK_VEL { isameAsBW, ishockVel };

enum METHODS_SHOCK_ELE { iShockEleAnalyt, iShockEleNum, iShockEleMix };

enum METHOD_TAU { iAPPROX, iTHICK, iSMOOTH, iSHARP };

enum METHODS_SYNCH { iWSPN99, iJOH06, iDER06, iMARG21, iNumeric, iGSL, iCSYN, iBessel };

enum METHOD_SSC { inoSSC, iNumSSC };

enum METHOD_ELE { iEleAnalytic, iEleNumeric };

enum METHODS_LFMIN { igmNum, igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };

enum METHODS_B { iBuseUb, iBasMAG21, iBuseGammaSh };

enum METHOD_LFCOOL{ iuseConst, iuseTe, iuseTcomov };

enum METHOD_LFMAX { iConst, iuseB };

enum METHODS_SSA { iSSAoff, iSSAon };

enum METHOD_PP { iPPoff, iPPnum };

enum METHOD_NONRELDIST{ inone, iuseGm, iuseAyache, iuseSironi };

enum METHOD_NE{ iusenprime, iuseNe };




class ElectronAndRadiaionBase{
public:
    Vector & freqs_grid;
    Vector & gams_grid;
    Electrons ele = Electrons(gams_grid); // std::unique_ptr<State> ele = nullptr;
    Photons syn = Photons(freqs_grid); //std::unique_ptr<State> syn = nullptr;
    Photons ssc = Photons(freqs_grid); //std::unique_ptr<State> ssc = nullptr;
    Photons total_rad = Photons(freqs_grid);
    /// ------------------------------------
    double B=-1, gamma_min=-1, gamma_max=-1, gamma_c=-1;


    double n_prime=-1; // used in B and gamma_min
    double n_protons=-1; // conditionally used on synchrotron emissivity
    double nn=-1; // Ne or nprime depending on the setting 'm_method_ne'
    double eprime=-1.,Gamma=-1,GammaSh=-1, betaSh=-1.,t_e=-1.,t_b=-1;
    double accel_frac=1.;
    double Theta=-1, z_cool=-1, x=-1;
    size_t max_substeps=0; // limit to the number of substeps to evolve electrons
    METHOD_TAU method_tau{}; double beta_min = -1;
    METHODS_SHOCK_ELE m_eleMethod{}; bool num_ele_use_adi_loss = false;
    METHOD_NE m_method_ne{};
//    bool do_ssa{};
    METHODS_SHOCK_VEL method_shock_vel{};
    METHODS_Up_sh m_method_up{};
    METHOD_thickness_sh m_method_Delta{};
//    METHOD_Gamma_sh m_method_gamma_fsh{};
//    METHOD_Gamma_sh m_method_gamma_rsh{};
//    METHOD_R_sh m_method_r_sh{};
    /// --------------------------------------
    void setBasePars( StrDbMap & pars, StrStrMap & opts ){
//        n_layers = nlayers;
        // set parameters
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";
        else
            fs_or_rs += "_fs";

//        ksi_n = getDoublePar("ksi_n" + fs_or_rs, pars, AT, p_log, 1., false);//pars.at("ksi_n");
        eps_e = getDoublePar("eps_e" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_e");
        eps_b = getDoublePar("eps_b" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_b");
        eps_t = getDoublePar("eps_t" + fs_or_rs, pars, AT, p_log, 0., true);//pars.at("eps_t");
        p = getDoublePar("p" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("p");
        mu = 0.62;//getDoublePar("mu" + fs_or_rs, pars, AT, p_log, 0.62, false);//pars.at("mu");
        mu_e = 1.18;//getDoublePar("mu_e" + fs_or_rs, pars, AT, p_log, 1.18, false);//pars.at("mu_e");
        beta_min = getDoublePar("beta_min" + fs_or_rs, pars, AT, p_log, 1.e-5, false);//pars.at("beta_min");
        gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars, AT, p_log, 1.e7, false);//pars.at("beta_min");
        max_substeps = (size_t)getDoublePar("max_substeps" + fs_or_rs, pars, AT, p_log, 1000, false);//pars.at("beta_min");

//        lim_gm_to_1 = getBoolOpt("limit_lf_min_to1" + fs_or_rs, opts, AT, p_log, false, false);//pars.at("beta_min");

        // set options
        std::string opt;

        /// method for shock radius
//        opt = "method_radius" + fs_or_rs;
//        METHOD_R_sh method_r_sh;
//        if ( opts.find(opt) == opts.end() ) {
//            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_r_sh = METHOD_R_sh::isameAsR;
//        }
//        else{
//            if(opts.at(opt) == "sameAsR")
//                method_r_sh = METHOD_R_sh::isameAsR;
////            else if(opts.at(opt) == "useGammaSh")
////                method_r_sh = METHOD_R_sh::iuseGammaSh;
//            else{
//                (*p_log)(LOG_ERR,AT) << " option for: " << opt
//                                     <<" given: " << opts.at(opt)
//                                     << " is not recognized "
//                                     << " Possible options: "
//                                     << " sameAsR " //<< " useGammaSh "
//                                     << " Exiting...\n";
////                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        m_method_r_sh = method_r_sh;

        /// method for shock velocity
//        opt = "method_Gamma" + fs_or_rs;
//        METHOD_Gamma_sh method_gamma_fsh;
//        if ( opts.find(opt) == opts.end() ) {
//            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_gamma_fsh = METHOD_Gamma_sh::iuseGammaShock;
//        }
//        else{
//            if(opts.at(opt) == "useJustGamma")
//                method_gamma_fsh = METHOD_Gamma_sh::iuseJustGamma;
//            else if(opts.at(opt) == "useGammaShock")
//                method_gamma_fsh = METHOD_Gamma_sh::iuseGammaShock;
//            else if(opts.at(opt) == "useJustGammaRel") {
////                if (!(m_type == BW_TYPES::iFS_DENSE || m_type == BW_TYPES::iFS_PWN_DENSE)){
////                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useJustGammaRel "<< " if bw_type ="<<p_pars->m_type<<"\n";
////                    exit(1);
////                }
//                method_gamma_fsh = METHOD_Gamma_sh::iuseJustGammaRel;
//            }
//            else if(opts.at(opt) == "useGammaRelShock") {
////                if (!(p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE)){
////                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useGammaRelShock "<< " if bw_type ="<<p_pars->m_type<<"\n";
////                    exit(1);
////                }
//                method_gamma_fsh = METHOD_Gamma_sh::iuseGammaRelShock;
//            }
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
//        m_method_gamma_fsh = method_gamma_fsh;

        /// method for energy density behind shock
        opt = "method_Up" + fs_or_rs;
        METHODS_Up_sh method_up;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_up = METHODS_Up_sh::iuseEint2;
        }
        else{
            if(opts.at(opt) == "useEint2")
                method_up = METHODS_Up_sh::iuseEint2;
            else if(opts.at(opt) == "useGamma")
                method_up = METHODS_Up_sh::iuseGamma;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " useEint2 " << " useGamma "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        m_method_up = method_up;

        /// method for shock thickness
        opt = "method_thickness" + fs_or_rs;
        METHOD_thickness_sh method_delta;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_delta = METHOD_thickness_sh::iuseJoh06;
        }
        else{
            if(opts.at(opt) == "useJoh06")
                method_delta = METHOD_thickness_sh::iuseJoh06;
            else if(opts.at(opt) == "useVE12")
                method_delta = METHOD_thickness_sh::iuseVE12;
            else if(opts.at(opt) == "None")
                method_delta = METHOD_thickness_sh::iNoDelta;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " useJoh06 " << " useVE12 " << " None "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        m_method_Delta = method_delta;

        /// should the shock gamma be the bw gamma, or from RH jump conditions
        opt = "method_vel" + fs_or_rs;
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
        method_shock_vel = methodsShockVel;

        opt = "method_ele" + fs_or_rs;
        METHODS_SHOCK_ELE val_ele;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_ele = METHODS_SHOCK_ELE::iShockEleAnalyt;
        }
        else{
            if(opts.at(opt) == "analytic")
                val_ele = METHODS_SHOCK_ELE::iShockEleAnalyt;
            else if(opts.at(opt) == "numeric")
                val_ele = METHODS_SHOCK_ELE::iShockEleNum;
            else if(opts.at(opt) == "mix")
                val_ele = METHODS_SHOCK_ELE::iShockEleMix;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " analytic " << " numeric " << " mix "<< "\n";
                exit(1);
            }
        }
        m_eleMethod = val_ele;

        num_ele_use_adi_loss = getBoolOpt("num_ele_use_adi_loss"+fs_or_rs,opts,AT,p_log,true,false);

        /// how to compute number of electrons that are emitting
        opt = "method_ne" + fs_or_rs;
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
        m_method_ne = methodNe;

        opt = "method_synchrotron" + fs_or_rs;
        METHODS_SYNCH val_synch;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_synch = METHODS_SYNCH::iJOH06;
        }
        else{
            if(opts.at(opt) == "Joh06")
                val_synch = METHODS_SYNCH::iJOH06;
            else if(opts.at(opt) == "WSPN99")
                val_synch = METHODS_SYNCH::iWSPN99;
            else if(opts.at(opt) == "Marg21")
                val_synch = METHODS_SYNCH::iMARG21;
            else if(opts.at(opt) == "Dermer09")
                val_synch = METHODS_SYNCH::iDER06;
            else if(opts.at(opt) == "GSL") {
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt){
                    (*p_log)(LOG_ERR,AT) << " options synchrotron GSL and analytic electron evol. are incompatible\n";
                    exit(1);
                }
                val_synch = METHODS_SYNCH::iGSL;
            }
            else if(opts.at(opt) == "Bessel") {
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt){
                    (*p_log)(LOG_ERR,AT) << " options synchrotron Bessel and analytic electron evol. are incompatible\n";
                    exit(1);
                }
                val_synch = METHODS_SYNCH::iBessel;
            }
            else if(opts.at(opt) == "CSYN") {
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt){
                    (*p_log)(LOG_ERR,AT) << " options synchrotron CSYN and analytic electron evol. are incompatible\n";
                    exit(1);
                }
                val_synch = METHODS_SYNCH::iCSYN;
            }
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " Joh06 " << " WSPN99 " << " Marg21 " << " Dermer09 " << " GSL " << " CSYN " << "\n";
                exit(1);
            }
        }
        m_sychMethod = val_synch;
        do_th_marg21 = (bool) getBoolOpt("incl_th_in_marg21"+fs_or_rs,opts,AT,p_log, true, false);

        opt = "method_ssc" + fs_or_rs;
        METHOD_SSC val_ssc;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_ssc = METHOD_SSC::inoSSC;
//            do_ssa = false;
        }
        else{
            if(opts.at(opt) == "none") {
                val_ssc = METHOD_SSC::inoSSC;
//                do_ssa = false;
            }
            else if(opts.at(opt) == "numeric") {
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleAnalyt){
                    (*p_log)(LOG_ERR,AT) << " SSC is not supported for analytic electron dsitrib. model\n";
                    exit(1);
                }
                val_ssc = METHOD_SSC::iNumSSC;
//                do_ssa = true;
            }
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " none " << " numeric " << "\n";
                exit(1);
            }
        }
        m_methods_ssc = val_ssc;

        opt = "method_pp" + fs_or_rs;
        METHOD_PP val_pp;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_pp = METHOD_PP::iPPoff;
        }
        else{
            if(opts.at(opt) == "none") {
                val_pp = METHOD_PP::iPPoff;
            }
            else if(opts.at(opt) == "numeric") {
                if (m_methods_ssc == METHOD_SSC::inoSSC){
                    (*p_log)(LOG_ERR,AT) << " PP requires SSC to be 'numeric'\n";
                    exit(1);
                }
                val_pp = METHOD_PP::iPPnum;
//                do_ssa = true;
            }
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " none " << " numeric " << "\n";
                exit(1);
            }
        }
        m_methods_pp = val_pp;

        opt = "method_nonrel_dist" + fs_or_rs;
        METHOD_NONRELDIST val_monreldist;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for " << opt << " is not set. Using default value.\n";
            val_monreldist = METHOD_NONRELDIST::inone;
        }
        else{
            if(opts.at(opt) == "none")
                val_monreldist = METHOD_NONRELDIST::inone;
            else if(opts.at(opt) == "use_gamma_min")
                val_monreldist = METHOD_NONRELDIST::iuseGm;
            else if(opts.at(opt) == "use_Ayache")
                val_monreldist = METHOD_NONRELDIST::iuseAyache;
            else if(opts.at(opt) == "use_Sironi")
                val_monreldist = METHOD_NONRELDIST::iuseSironi;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " none " << " use_Ayache " << " use_Sironi " << " use_gamma_min " << "\n";
                exit(1);
            }
        }
        m_method_nonreldist = val_monreldist;

        opt = "method_gamma_min" + fs_or_rs;
        METHODS_LFMIN val_lfmin;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_lfmin = METHODS_LFMIN::igmUprime;
        }
        else{
            if(opts.at(opt) == "useU_e")
                val_lfmin = METHODS_LFMIN::igmUprime;
            else if(opts.at(opt) == "useNumeric")
                val_lfmin = METHODS_LFMIN::igmNum;
            else if(opts.at(opt) == "useBeta")
                val_lfmin = METHODS_LFMIN::igmNakarPiran;
            else if(opts.at(opt) == "useGamma")
                val_lfmin = METHODS_LFMIN::igmJoh06;
            else if(opts.at(opt) == "useTheta")
                val_lfmin = METHODS_LFMIN::igmMAG21;
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useU_e " << " useTheta "<< " useBeta " << " useGamma " << "useNumericGamma" << "\n";
                exit(1);
            }
        }
        m_methodsLfmin = val_lfmin;

        opt = "method_B" + fs_or_rs;
        METHODS_B methodsB;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsB = METHODS_B::iBuseUb;
        }
        else{
            if(opts.at(opt) == "useU_b")
                methodsB = METHODS_B::iBuseUb;
            else if(opts.at(opt) == "useMAG21")
                methodsB = METHODS_B::iBasMAG21;
            else if(opts.at(opt) == "useGammaSh")
                methodsB = METHODS_B::iBuseGammaSh;
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useUb " << " useMAG21 " << " useGammaSh " << "\n";
                exit(1);
            }
        }
        m_methodsB = methodsB;

        opt = "method_gamma_max" + fs_or_rs;
        METHOD_LFMAX methodLfmax;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodLfmax = METHOD_LFMAX::iuseB;
        }
        else{
            if(opts.at(opt) == "useB")
                methodLfmax = METHOD_LFMAX::iuseB;
            else if(opts.at(opt) == "useConst") {
                gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars,
                                         AT, p_log, 1.e7, true);//pars.at("beta_min");
                methodLfmax = METHOD_LFMAX::iConst;
            }
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useB " << " useConst " << "\n";
                exit(1);
            }
        }
        m_methodsLfmax = methodLfmax;

        opt = "method_gamma_c" + fs_or_rs;
        METHOD_LFCOOL methodLfcool;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodLfcool = METHOD_LFCOOL::iuseTe;
        }
        else{
            if(opts.at(opt) == "useTe")
                methodLfcool = METHOD_LFCOOL::iuseTe;
            else if(opts.at(opt) == "useTcomov")
                methodLfcool = METHOD_LFCOOL::iuseTcomov;
            else if(opts.at(opt) == "useConst") {
                gamma_c = getDoublePar("gamma_c" + fs_or_rs, pars,
                                       AT, p_log, 1.e5, true);//pars.at("beta_min");
                methodLfcool = METHOD_LFCOOL::iuseConst;
            }
            else{
                (*p_log)(LOG_ERR,AT) << AT << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " useTe " << " useTcomov " << " useConst "<< "\n";
                exit(1);
            }
        }
        m_methodsLfcool = methodLfcool;

        bool tmp = getBoolOpt("use_ssa" + fs_or_rs, opts, AT, p_log, false, true);
        if (tmp) m_methods_ssa = METHODS_SSA::iSSAon;
        else m_methods_ssa = METHODS_SSA::iSSAoff;

        opt = "method_tau" + fs_or_rs;
        METHOD_TAU methodTau;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodTau = METHOD_TAU::iSMOOTH;
        }
        else{
            if(opts.at(opt) == "thick")
                methodTau = METHOD_TAU::iTHICK;
            else if(opts.at(opt) == "sharp")
                methodTau = METHOD_TAU::iSHARP;
            else if(opts.at(opt) == "smooth")
                methodTau = METHOD_TAU::iSMOOTH;
            else if(opts.at(opt) == "approx")
                methodTau = METHOD_TAU::iAPPROX;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << " Possible options: "
                                     << " approx " << " thick " << " sharp " << " smooth " << "\n";
                exit(1);
            }
        }
        method_tau = methodTau;

    }

    /// store current shock properties
    void updateSockProperties(double e_prime, double Gamma_, double Gamma_shock, double t_e_, double t_b_,
                              double tcomov_, double n_prime_, double n_protons_){
        /// Store current parameters

        GammaSh = Gamma_shock;
        eprime = e_prime;
        Gamma = Gamma_;
        betaSh = EQS::BetaFromGamma(Gamma_shock);
        t_e = t_e_;
        t_b = t_b_;
        tcomov = tcomov_;
        n_prime = n_prime_;
        n_protons = n_protons_;
//        r = r_;
//        dr = dr_;

        switch (m_method_ne) {
            case iusenprime:
                nn = n_prime;
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT)<<" cannot use shock electron density for "
                                          "numeric electron evolution. Must use actaul electron number, Ne.\n";
                    exit(1);
                }
                break;
            case iuseNe:
                nn = n_protons;
                break;
        }

    }
    /// In case the electrons were computed elsewhere E.g., if interpolated for EATS plane
    void setShockElectronParameters(double n_prime_, double n_protons_, double acc_frac,
                                    double B_, double gm, double gM, double gc,
                                    double Theta_, double z_cool_){

        if (!std::isfinite(gm) || !std::isfinite(gc) || !std::isfinite(n_prime_)) {
            (*p_log)(LOG_ERR, AT) << " nans is synchrotron spectum\n";
            exit(1);
        }

        gamma_min = gm;
        gamma_max = gM;
        n_prime = n_prime_;
        n_protons = n_protons_;
        B = B_;
        accel_frac = acc_frac;
        gamma_c = gc;
        z_cool = z_cool_;
        Theta = Theta_;
//        r = r_;
//        dr = dr_;
        switch (m_method_ne) {
            case iusenprime:
                nn = n_prime_;
                if (m_eleMethod == METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT)<<" cannot use shock electron density for "
                                          "numeric electron evolution. Must use actaul electron number, Ne.\n";
                    exit(1);
                }
                break;
            case iuseNe:
                nn = n_protons_;
                break;
        }
    }
    /// evaluate frequency independent quantities (critical LFs, Bfield, etc)
    void computeEleSpectrumLimits() {

        if (not checkParams())
            return;
        /// compute comoving magnetic field
        computeMagneticField(m_methodsB);
        /// compute limits of electron spectrum
        computeGammaMax(m_methodsLfmax);
        computeGammaMin(m_methodsLfmin);
        computeGammaCool(m_methodsLfcool);
    }

    double inline get_shock_thickness(double R, double M2, double theta, double Gamma, double rho2, double ncells) const {
        double thickness;
        double delta_joh06 = shock_delta_joh06(R,M2,theta,Gamma,rho2,ncells);
        double delta_ve12 = shock_delta(R, Gamma);
        switch (m_method_Delta) {
            case iuseJoh06:
                thickness = delta_joh06;
                break;
            case iuseVE12:
                thickness = delta_ve12;
                break;
            case iNoDelta:
                thickness = 1.;
                break;
        }
        if (!std::isfinite(thickness)){
            (*p_log)(LOG_ERR,AT) << " shock thickness = " << thickness << "\n";
            exit(1);
        }
        return thickness;
    }

    double inline get_shock_Up(double GammaShock, double rho2, double m2, double Eint2){
        if (Eint2 == 0 || m2 == 0)
            return 0.;
        double up; //comoving energy density (MF)%
        //        double rhoprim = 4. * rho * Gamma ;     // comoving density
        double V2 = m2 / rho2 ; // comoving volume
        double U_p = Eint2 / V2 ; // comoving energy density (electrons)
        double up_from_internal_energy = U_p;
        double up_from_lorentz_factor = (GammaShock - 1.) * rho2 * CGS::c * CGS::c;;
        switch(m_method_up){
            case iuseEint2:
                up = up_from_internal_energy;
                break;
            case iuseGamma:
                up = up_from_lorentz_factor;
                break;
        }
        if (up < 0 || !std::isfinite(up)){
            (*p_log)(LOG_ERR,AT) << " wrong value Up = "<<up<<"\n";
            exit(1);
        }
        return up;
    }

    /// ---
    [[nodiscard]] ElectronAndRadiaionBase const * getThis() const { return this; }

public:
    /// --------------------------------------
    METHODS_SYNCH m_sychMethod{};
    METHODS_LFMIN m_methodsLfmin{};
    METHOD_LFCOOL m_methodsLfcool{};
    METHOD_LFMAX m_methodsLfmax{};
    METHODS_B m_methodsB{};
    METHOD_NONRELDIST m_method_nonreldist{};
    METHODS_SSA m_methods_ssa{};
    METHOD_SSC m_methods_ssc{};
    METHOD_PP m_methods_pp{};

    double eps_e=-1;
protected:
    int m_loglevel = -1;
    std::unique_ptr<logger> p_log = nullptr;
    bool is_rs = false;
    bool do_th_marg21 = true;
    double tcomov0=-1.;
    double tcomov=-1;
    /// --------------------------------------
    double eps_b=-1, eps_t=-1, p=-1;// ksi_n=-1;
    double mu=-1, mu_e=-1;
//    bool lim_gm_to_1= true;
    double gam1=-1,gam2=-1.,freq1=-1.,freq2=-1.;
    size_t ngam=0,nfreq=0;
    /// --------------------------------------
    ElectronAndRadiaionBase(Vector & gams_grid, Vector & freqs_grid, int loglevel, bool _is_rs)
            : gams_grid(gams_grid), freqs_grid(freqs_grid) {
        m_loglevel = loglevel; is_rs = _is_rs;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronAnalytic");
    }
    /// --------------------------------------
    bool checkParams(){

        /// if velocity is too small shock may not be possible TODO update it with proper sound speed check
        if (betaSh < beta_min)
            return false;


        /// in Margalit+21 model thermal electrons are present. Different treatment :: Assert:
        if (m_sychMethod == METHODS_SYNCH::iMARG21) {
            /* Margalit+21 arXiv:2111.00012 */
            if (m_methodsLfmin != igmMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsLfmin=gmMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
            if (m_methodsB != iBasMAG21) {
                (*p_log)(LOG_ERR, AT) << "use m_methodsB=asMAG21 when using m_sychMethod=MARG21\n";
                exit(1);
            }
        }

        // check
        if (eps_e <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<eps_e<<")\n";
            exit(1);
        }
        if (eps_b <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_b is not set (eps_e="<<eps_b<<")\n";
            exit(1);
        }
        if (eps_t < 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<eps_t<<") (even through it is only needed in Marg21 model)\n";
            exit(1);
        }
        if (p <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<p<<")\n";
            exit(1);
        }
//        if (ksi_n <= 0){
//            (*p_log)(LOG_ERR,AT)<< " ksi_n is not set (ksi_n="<<ksi_n<<")\n";
//            exit(1);
//        }
        if (n_prime <= 0){
            (*p_log)(LOG_ERR,AT) << " n_prime is not set (n_prime="<<n_prime<<")\n";
            exit(1);
        }
        if ((eprime < 0) || !std::isfinite(eprime)){
            (*p_log)(LOG_ERR,AT) << " eprime is not set (eprime=" << eprime << ")\n";
            exit(1);
        }

        return true;
    }
    /// --------------------------------------
    /*
     * Minimum lorentz factor of the injection electron distribution (See Eq.8.76 - 8.85 in Zhang+2018 book)
     */
    static double gammaMinFunc(const double &x, void * pars){
        auto * pp = (struct ElectronAndRadiaionBase *) pars;
        double gamma_ave = pp->eps_e * pp->eprime / (pp->n_prime * CGS::me * CGS::c * CGS::c);//* mp / me * (pp->GammaSh - 1)
        if (pp->p == 1)
            return (pp->gamma_max - pp->gamma_min) / (std::log(pp->gamma_max) - std::log(pp->gamma_min))
                   - gamma_ave;
        else if (pp->p == 2)
            return (std::log(pp->gamma_max) - std::log(pp->gamma_min)) / (-1./pp->gamma_max + 1./pp->gamma_min)
                   - gamma_ave;
        else
            return (pp->p - 1) / (pp->p - 2)
                   * (std::pow(pp->gamma_max, -pp->p + 2) - std::pow(x, -pp->p + 2))
                   / (std::pow(pp->gamma_max, -pp->p + 1) - std::pow(x, -pp->p + 1))
                   - gamma_ave;
    }
    /// --------------------------------------
    void computeMagneticField(METHODS_B methodsB){
        /// fraction of the shock energy density in electrons
        double U_b_prime = eprime * eps_b;
        double B1 = sqrt(8. * CGS::pi * U_b_prime);
        /// See Nava 2013 paper or others; classical equation
        double B2 = std::sqrt( 32. * M_PI * eps_b * (GammaSh - 1.)
                               * (GammaSh + 3. / 4) * n_prime * CGS::mp * CGS::c2);
        double B3 = sqrt(9.0 * M_PI * eps_b * n_prime * mu * CGS::mp)
                    * (betaSh * CGS::c);

        switch (methodsB) {
            case iBuseUb:       B=B1; break;
            case iBuseGammaSh:  B=B2; break;
            case iBasMAG21:     B=B3; break;
        }
    }
    void computeGammaMax(METHOD_LFMAX methodLfmax){
        /// get electron distribution boundaries (has ot be set BEFORE gamma_min)
        double gamma_max_1 = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / B); // Kumar+14
        switch (methodLfmax) {
            case iConst:
                break;
            case iuseB:
                gamma_max = gamma_max_1;
                break;
        }
        if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleAnalyt)
            if (gamma_max > 0.99 * ele.e[ele.numbins - 1])
                gamma_max = 0.99 * ele.e[ele.numbins - 1];
    }
    void computeGammaMin(METHODS_LFMIN methodsLfmin){
        /// compute injection LF (lower bound of source injection gammaMinFunc)
        int status = 0; double U_e_prime = eprime * eps_e;
        /// Eq. A18 vanEarten+10 (and Sironi+13)
        double gamma_min_1 = (p - 2.) / (p - 1.) * U_e_prime / (n_prime * CGS::me * CGS::c * CGS::c);
        /// Eq.(before 4.) in Nakar&Piran 1102.1020
        double gamma_min_2 = (p - 2.) / (p - 1.) * (CGS::mp / CGS::me) * eps_e * EQS::BetaFromGamma(GammaSh);
        /// Eq. A3 in J+06
        double gamma_min_3 = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (GammaSh - 1.) + 1.);
        /// downstream electron temperature:
        Theta = Margalit21::Theta_fun(betaSh, mu, mu_e, eps_t);
        /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21): NOTE we use mean dynamical time:
        z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * B * B * t_e)) / Theta;
        /// Margalit Piran 2022
        double gamma_min_4 = Margalit21::gamma_m_fun(Theta);
        /// switch
        switch (methodsLfmin) {
            case igmUprime:
                gamma_min = gamma_min_1;
                break;
            case igmNakarPiran:
                gamma_min = gamma_min_2;
                break;
            case igmJoh06:
                gamma_min = gamma_min_3;
                break;
            case igmNum:
                /// solve gamma_min fun numerically; use fixed limits and number of iterations
                gamma_min = Bisect(ElectronAndRadiaionBase::gammaMinFunc,
                                   1., 1.e8, 0., .001, 100, this, status);
                /// If numerical solution failed for expected large gamma_min, use simple analytical solution
                if ((status < 0) && (gamma_min_1 > 1.5) ) {
                    (*p_log)(LOG_ERR,AT) << "Bisect for gamma_min failed. Obtained = " << gamma_min
                                         << " while gamma_min_1 = "<< gamma_min_1
                                         << " \n";
                    gamma_min = gamma_min_1;
                }
                break;
            case igmMAG21:
                gamma_min = gamma_min_4;
                break;
        }



        // Deep Newtonian Phase
        double accel_frac_ayache = (std::pow(gamma_max, 2 - p) - std::pow(gamma_min - 1., 2. - p))
                            / (std::pow(gamma_max, 2. - p) - 1.)
                            * (std::pow(gamma_max, 1. - p) - 1.) / (std::pow(gamma_max, 1. - p)
                                                                    - std::pow(gamma_min - 1., 1. - p));
        double accel_frac_sironi = (p - 2.0) / (p - 1.0) * eps_e * eprime / (n_prime * CGS::me * CGS::c * CGS::c);


        if ((gamma_min < 2.) || (status < 0)) {
            switch (m_method_nonreldist) {
                case inone:
                    accel_frac = 1;
                    break;
                case iuseGm:
                    accel_frac = gamma_min;
                    break;
                case iuseAyache:
                    if (p < 2.){
                        (*p_log)(LOG_ERR,AT) << "When setting p < 2, Newtonain regime is not available." <<"\n";
                        exit(1);
                    }
                    accel_frac = accel_frac_ayache;
                    break;
                case iuseSironi: // Sironi et al 2013
                    if (p < 2.){
                        (*p_log)(LOG_ERR,AT) << "When setting p < 2, Newtonain regime is not available." <<"\n";
                        exit(1);
                    }
                    accel_frac = accel_frac_sironi;
                    break;
            }
            accel_frac = accel_frac > 1. ? 1 : accel_frac;
            gamma_min = gamma_min < 1. ? 1 : gamma_min;
        }

        /// check
        if (!std::isfinite(gamma_min) || (!std::isfinite(accel_frac)  || accel_frac < 0)){
            (*p_log)(LOG_ERR,AT) << " Wrong value gamma_min="<<gamma_min<<" accel_frac="<<accel_frac<<"\n";
            exit(1);
        }

//        printf("gamma_min to %.2e GammaSh=%.2e eprime=%.2e p=%.2f gamma_max=%.2e ", gamma_min, GammaSh, eprime, p, gamma_max);

        /// prevent gamma_max to go beyond the electron grid
        if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleAnalyt)
            if (gamma_min < ele.e[0]) {
                gamma_min = ele.e[0];
                if (gamma_min > 1e8){
                    (*p_log)(LOG_WARN,AT)<< " uphysical value in gamma_min="<<gamma_min<<"\n";
                    exit(1);
                }
            }
//        printf("gamma_min to %.2e GammaSh=%.2e eprime=%.2e p=%.2f gamma_max=%.2e ", gamma_min, GammaSh, eprime, p, gamma_max);

    }
    void computeGammaCool(METHOD_LFCOOL methodLfcool){
        // Eq. A19 in vanEarten+10
        double gamma_c_1 = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * B * B) * t_e / Gamma;

        double gamma_c_2 = ((6. * CGS::pi * me * c) / (CGS::sigmaT * B * B * (tcomov))); // (tcomov0-1.e-6)

        switch (methodLfcool) {
            case iuseConst:
                break;
            case iuseTe:
                gamma_c = gamma_c_1;
                break;
            case iuseTcomov:
                gamma_c = gamma_c_2;
                break;
        }
//        if (gamma_c < 1. || !std::isfinite(gamma_c)){
//            (*p_log)(LOG_ERR,AT)<<" gamma_c="<<gamma_c<<"\n";
//            exit(1);
//        }
//        if (gamma_c >= gamma_max)
//            gamma_c = 0.999 * gamma_max;
        if (gamma_c <= 1.)
            gamma_c = 1.001;
        if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleAnalyt)
            if (gamma_c > gamma_max)
                gamma_c = gamma_max;

    }

};

#endif //SRC_MICROPHYSICS_BASE_H
