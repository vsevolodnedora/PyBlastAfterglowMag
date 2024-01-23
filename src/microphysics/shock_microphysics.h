//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_SHOCK_MICROPHYSICS_H
#define SRC_SHOCK_MICROPHYSICS_H

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

enum METHODS_SHOCK_ELE { iShockEleAnalyt, iShockEleNum };

enum METHOD_TAU { iAPPROX, iTHICK, iSMOOTH, iSHARP };

enum METHODS_SYNCH { iWSPN99, iJOH06, iDER06, iMARG21, iNumeric, iGSL };

enum METHOD_SSC { inoSSC, iNumSSC };

enum METHOD_ELE { iEleAnalytic, iEleNumeric };

enum METHODS_LFMIN { igmNumGamma, igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };

enum METHODS_B { iBuseUb, iBasMAG21, iBuseGammaSh };

enum METHOD_LFMAX { iConst, iuseB };

enum METHODS_SSA { iSSAoff, iSSAon };

enum METHOD_NONRELDIST{ inone, iuseGm };

enum METHOD_NE{ iusenprime, iuseNe };

/// Evaluate optical depth
inline double optical_depth(const double abs_lab, const double dr,
                            const double mu, const double beta_shock){
    if (abs_lab <= 0.)
        return 0.;
    double dtau;
    if(mu == beta_shock)
        dtau = 1.0e100; // HUGE VAL, just in case
    else
        dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // basic dtau = abs_lab * dr;
    return dtau;
}

/// Compute the I = emission/absorption * ( 1 - tau^{-absorption * thickness})
double computeIntensity(const double em_lab, double dtau, METHOD_TAU m_tau){

    if ( dtau == 0 )
        return em_lab;

    double intensity = em_lab;

    // (Signed) Optical depth through this shell.
    // if negative, face is oriented away from observer.
//        double dtau;
//        if(mu == beta_shock)
//            dtau = 1.0e100; // HUGE VAL, just in case
//        else
//            dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // basic dtau = abs_lab * dr;

    // Now that we know the optical depth, we apply it in a way
    // according to the given specType
    if(m_tau == iTHICK) {
        // Special case: use the optically thick limit *everywhere*
        if(dtau <= 0.0)
            intensity = 0.0;
        else
            intensity /= dtau;
    }
    else if (m_tau == iAPPROX){
        double u = 1. / 2. + exp(-dtau) / dtau - (1 - exp(-dtau)) / (dtau * dtau);
        if (dtau < 1e-3) {
            intensity *= 1.;
        }
        else {
            intensity *= (3. * u / dtau);
        }
//            return np.where(tau < 1e-3, 1., 3 * u / tau)
    }
    else if(m_tau == iSMOOTH) {
        // Apply self-absorption "properly"
        //
        // correction factor to emissivity from absorption
        // ( 1 - e^(-tau) ) / tau  (on front face)
        //
        // back face has extra factor ~ e^-beta_shock/(mu-beta_shock)
        //
        // for now ignoring shadowing by the front face.
        double abs_fac;
        if(dtau == 0.0)
            abs_fac = 1.0;
        else if(dtau > 0.0)
            abs_fac = -expm1(-dtau) / dtau;
        else {
            abs_fac = expm1(dtau) / dtau; //* exp(
            //abs * DR * beta_shock*mu / (mu - beta_shock));
        }

        intensity *= abs_fac;
    }
    else if(m_tau == iSHARP){
        // Apply self-absorption "simply".
        //
        // Compute flux in optically thick limit,
        // use in final result if less than optically thin calculation.
        //
        // e.g. use tau->infty limit if tau > 1.0

        // "Forward" face
        if(dtau > 1.0)
            intensity /= dtau;
            // "Back" face --> assume shadowed by front
        else if(dtau < -1.0)
            intensity = 0.0;
    }

    return intensity;
}


class ElectronAndRadiaionBase{
public:
    Vector out_spectrum{}; Vector out_specturm_ssa{};
    Vector m_freq_arr{};
    /// ------------------------------------
    double B=-1, gamma_min=-1, gamma_max=-1, gamma_c=-1;
    double n_prime=-1; // used in B and gamma_min
    double n_protons=-1; // conditionally used on synchrotron emissivity
    double nn=-1; // Ne or nprime depending on the setting 'm_method_ne'
    double eprime=-1.,Gamma=-1,GammaSh=-1, beta=-1.,t_e=-1.,r=-1,dr=-1;
    double accel_frac=-1.;
    double Theta=-1, z_cool=-1, x=-1;
    size_t max_substeps=0; // limit to the number of substeps to evolve electrons
    METHOD_TAU method_tau{}; double beta_min = -1;
    METHODS_SHOCK_ELE m_eleMethod{};
    METHOD_NE m_method_ne{};
    bool do_ssa{};
    METHODS_SHOCK_VEL method_shock_vel{};
    METHODS_Up_sh m_method_up{};
    METHOD_thickness_sh m_method_Delta{};
    METHOD_Gamma_sh m_method_gamma_fsh{};
    METHOD_Gamma_sh m_method_gamma_rsh{};
    METHOD_R_sh m_method_r_sh{};
    /// --------------------------------------
    void setBasePars( StrDbMap & pars, StrStrMap & opts ){
//        n_layers = nlayers;
        // set parameters
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";
        else
            fs_or_rs += "_fs";

        ksi_n = getDoublePar("ksi_n" + fs_or_rs, pars, AT, p_log, 1., false);//pars.at("ksi_n");
        eps_e = getDoublePar("eps_e" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_e");
        eps_b = getDoublePar("eps_b" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_b");
        eps_t = getDoublePar("eps_t" + fs_or_rs, pars, AT, p_log, 0., true);//pars.at("eps_t");
        p = getDoublePar("p" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("p");
        mu = getDoublePar("mu" + fs_or_rs, pars, AT, p_log, 0.62, false);//pars.at("mu");
        mu_e = getDoublePar("mu_e" + fs_or_rs, pars, AT, p_log, 1.18, false);//pars.at("mu_e");
        beta_min = getDoublePar("beta_min" + fs_or_rs, pars, AT, p_log, 1.e-5, false);//pars.at("beta_min");
        gamma_max = getDoublePar("gamma_max" + fs_or_rs, pars, AT, p_log, 1.e7, false);//pars.at("beta_min");
        max_substeps = (size_t)getDoublePar("max_substeps" + fs_or_rs, pars, AT, p_log, 1000, false);//pars.at("beta_min");

        lim_gm_to_1 = getBoolOpt("limit_lf_min_to1" + fs_or_rs, opts, AT, p_log, false, false);//pars.at("beta_min");

        // set options
        std::string opt;

        /// method for shock radius
        opt = "method_radius" + fs_or_rs;
        METHOD_R_sh method_r_sh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_r_sh = METHOD_R_sh::isameAsR;
        }
        else{
            if(opts.at(opt) == "sameAsR")
                method_r_sh = METHOD_R_sh::isameAsR;
            else if(opts.at(opt) == "useGammaSh")
                method_r_sh = METHOD_R_sh::iuseGammaSh;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " sameAsR " << " useGammaSh "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        m_method_r_sh = method_r_sh;

        /// method for shock velocity
        opt = "method_Gamma" + fs_or_rs;
        METHOD_Gamma_sh method_gamma_fsh;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            method_gamma_fsh = METHOD_Gamma_sh::iuseGammaShock;
        }
        else{
            if(opts.at(opt) == "useJustGamma")
                method_gamma_fsh = METHOD_Gamma_sh::iuseJustGamma;
            else if(opts.at(opt) == "useGammaShock")
                method_gamma_fsh = METHOD_Gamma_sh::iuseGammaShock;
            else if(opts.at(opt) == "useJustGammaRel") {
//                if (!(m_type == BW_TYPES::iFS_DENSE || m_type == BW_TYPES::iFS_PWN_DENSE)){
//                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useJustGammaRel "<< " if bw_type ="<<p_pars->m_type<<"\n";
//                    exit(1);
//                }
                method_gamma_fsh = METHOD_Gamma_sh::iuseJustGammaRel;
            }
            else if(opts.at(opt) == "useGammaRelShock") {
//                if (!(p_pars->m_type == BW_TYPES::iFS_DENSE || p_pars->m_type == BW_TYPES::iFS_PWN_DENSE)){
//                    (*p_log)(LOG_ERR,AT)<<" Cannot use "<<opt<<" = useGammaRelShock "<< " if bw_type ="<<p_pars->m_type<<"\n";
//                    exit(1);
//                }
                method_gamma_fsh = METHOD_Gamma_sh::iuseGammaRelShock;
            }
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized "
                                     << " Possible options: "
                                     << " useGammaShock " << " useJustGamma "
                                     << " Exiting...\n";
//                std::cerr << AT << "\n";
                exit(1);
            }
        }
        m_method_gamma_fsh = method_gamma_fsh;

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
                                     << " useJoh06 " << " useVE12 " << "None"
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
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " analytic " << " numeric " << "\n";
                exit(1);
            }
        }
        m_eleMethod = val_ele;

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
                if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT) << " options synchrotron GSL and analytic electron evol. are incompatible\n";
                    exit(1);
                }
                val_synch = METHODS_SYNCH::iGSL;

            }
//            else if(opts.at(opt) == "Bretta")
//                val_synch = SynchrotronAnalytic::METHODS_SYNCH::iBerrettaSynch;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " Joh06 " << " WSPN99 " << " Marg21 " << " Dermer09 " << " GSL " << "\n";
                exit(1);
            }
        }
        m_sychMethod = val_synch;

        opt = "method_ssc" + fs_or_rs;
        METHOD_SSC val_ssc;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_ssc = METHOD_SSC::inoSSC;
            do_ssa = false;
        }
        else{
            if(opts.at(opt) == "none") {
                val_ssc = METHOD_SSC::inoSSC;
                do_ssa = false;
            }
            else if(opts.at(opt) == "numeric") {
                if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleNum){
                    (*p_log)(LOG_ERR,AT) << " SSC is not supported for analytic electron dsitrib. model\n";
                    exit(1);
                }
                val_ssc = METHOD_SSC::iNumSSC;
                do_ssa = true;
            }
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " none " << " numeric " << "\n";
                exit(1);
            }
        }
        m_methods_ssc = val_ssc;

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
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " none " << " use_gamma_min " << "\n";
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
            else if(opts.at(opt) == "useNumericGamma")
                val_lfmin = METHODS_LFMIN::igmNumGamma;
            else if(opts.at(opt) == "useBeta")
                val_lfmin = METHODS_LFMIN::igmNakarPiran;
            else if(opts.at(opt) == "useGamma")
                val_lfmin = METHODS_LFMIN::igmJoh06;
            else if(opts.at(opt) == "useTheta")
                val_lfmin = METHODS_LFMIN::igmMAG21;
            else{
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useU_e " << " useTheta "<< " useBeta " << " useGamma " << "\n";
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
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
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
                (*p_log)(LOG_WARN,AT) << AT << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << "Possible options: "
                                      << " useB " << " useConst " << "\n";
                exit(1);
            }
        }
        m_methodsLfmax = methodLfmax;

        bool tmp = getBoolOpt("use_ssa" + fs_or_rs, opts, AT, p_log, false, true);
        if (tmp) m_methods_ssa = METHODS_SSA::iSSAon;
        else m_methods_ssa = METHODS_SSA::iSSAoff;

#if 0
        opt = "emissivity";
        SynchrotronAnalytic::QQ val_em;
        if ( opts.find(opt) == opts.end() ) {
            val_em = SynchrotronAnalytic::QQ::i_em;
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
        }
        else{
            if(opts.at(opt) == "em_pl")
                val_em = SynchrotronAnalytic::QQ::i_em_pl;
            else if(opts.at(opt) == "em_th")
                val_em = SynchrotronAnalytic::QQ::i_em_th;
            else if(opts.at(opt) == "em")
                val_em = SynchrotronAnalytic::QQ::i_em;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " em_pl " << " em_th " << " em " << "\n";
                exit(1);
            }
        }
        p_pars->m_marg21opt_em = val_em;

        opt = "absorption";
        SynchrotronAnalytic::QQ val_abs;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_WARN,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            val_abs = SynchrotronAnalytic::QQ::i_abs;
        }
        else{
            if(opts.at(opt) == "abs_pl")
                val_abs = SynchrotronAnalytic::QQ::i_abs_pl;
            else if(opts.at(opt) == "abs_th")
                val_abs = SynchrotronAnalytic::QQ::i_abs_th;
            else if(opts.at(opt) == "abs")
                val_abs = SynchrotronAnalytic::QQ::i_abs;
            else{
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                          <<" given: " << opts.at(opt)
                          << " is not recognized. "
                          << "Possible options: "
                          << " abs_pl " << " abs_th " << " abs " << "\n";
                exit(1);
            }
        }
        p_pars->m_marg21opt_abs = val_abs;
        // ---
#endif
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
                (*p_log)(LOG_WARN,AT) << " option for: " << opt
                                      <<" given: " << opts.at(opt)
                                      << " is not recognized. "
                                      << " Possible options: "
                                      << " approx " << " thick " << " sharp " << " smooth " << "\n";
                exit(1);
            }
        }
        method_tau = methodTau;
    }
    void allocateStorageForOutputSpectra(size_t nr, double freq1, double freq2,  size_t nfreq){
        m_freq_arr = TOOLS::MakeLogspaceVec(log10(freq1), log10(freq2),(int)nfreq);
        /// allocate space for the comoving spectrum
        out_spectrum.resize(nfreq * nr );
        out_specturm_ssa.resize(nfreq * nr );

    }
    /// store current shock properties
    void updateSockProperties(double e_prime, double Gamma_, double Gamma_shock, double t_e_,
                              double n_prime_, double n_protons_){
        /// Store current parameters

        GammaSh = Gamma_shock; // for bisect solver
        eprime = e_prime; // for bisect solver
        Gamma = Gamma_; // for bisect solver
        beta = EQS::Beta(Gamma_shock); // for bisect solver
        t_e = t_e_; // for bisect solver
        n_prime = ksi_n * n_prime_; //  for bisect solver
        n_protons = ksi_n * n_protons_;
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
    void evaluateElectronDistributionAnalytic() {

        if (not checkParams())
            return;

        computeMagneticField();

        computeGammaMax();

        computeGammaMin();

        computeNonRelativisticFlattening();

        /// compute cooling lorentz factor
        gamma_c = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * t_e * B * B) / Gamma; // Eq. A19 in vanEarten+10
    }
    /// ---
    [[nodiscard]] ElectronAndRadiaionBase const * getThis() const { return this; }
protected:
    int m_loglevel = -1;
    std::unique_ptr<logger> p_log = nullptr;
    bool is_rs = false;
    /// --------------------------------------
    double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
    double mu=-1, mu_e=-1;
    bool lim_gm_to_1= true;
    /// --------------------------------------
    METHODS_SYNCH m_sychMethod{};
    METHODS_LFMIN m_methodsLfmin{};
    METHODS_B m_methodsB{};
    METHOD_LFMAX m_methodsLfmax{};
    METHOD_NONRELDIST m_method_nonreldist{};
    METHODS_SSA m_methods_ssa{};
    METHOD_SSC m_methods_ssc{};
    /// --------------------------------------
    ElectronAndRadiaionBase(int loglevel, bool _is_rs){
        m_loglevel = loglevel; is_rs = _is_rs;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronAnalytic");
    }
    /// --------------------------------------
    bool checkParams(){

        /// if velocity is too small shock may not be possible TODO update it with proper sound speed check
        if (beta < beta_min)
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
        if (ksi_n <= 0){
            (*p_log)(LOG_ERR,AT)<< " ksi_n is not set (ksi_n="<<ksi_n<<")\n";
            exit(1);
        }
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
    static double gammaMinFunc(const double &x, void * pars){
        auto * pp = (struct ElectronAndRadiaionBase *) pars;
        return (pp->p - 1) / (pp->p - 2)
            * (std::pow(pp->gamma_max, -pp->p + 2) - std::pow(x, -pp->p + 2))
            / (std::pow(pp->gamma_max, -pp->p + 1) - std::pow(x, -pp->p + 1))
            - pp->eps_e * mp / me * (pp->GammaSh - 1);
    }
    /// --------------------------------------
    void computeMagneticField(){
        /// fraction of the shock energy density in electrons
        double U_b_prime;
        switch (m_methodsB) {
            case iBuseUb:
                U_b_prime = eprime * eps_b;
                B = sqrt(8. * CGS::pi * U_b_prime);
                break;
            case iBuseGammaSh:
                /// See Nava 2013 paper or others; classical equation
                B = std::sqrt( 32. * M_PI * eps_b * (GammaSh - 1.)
                  * (GammaSh + 3. / 4) * n_prime * CGS::mp * CGS::c2);
                break;
            case iBasMAG21:
                B = sqrt(9.0 * M_PI * eps_b * n_prime * mu * CGS::mp)
                  * (beta * CGS::c);
                break;
        }
    }
    void computeGammaMax(){
        /// get electron distribution boundaries (has ot be set BEFORE gamma_min)
        switch (m_methodsLfmax) {
            case iConst:
                break;
            case iuseB:
                gamma_max = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / B); // Kumar+14
                break;
        }
    }
    void computeGammaMin(){
        /// compute injection LF (lower bound of source injection gammaMinFunc)
        int status = 0; double U_e_prime = -1;
//        void * pars = (struct ElectronAndRadiaionBase *) this; // this "this" i do not understand myself...
        switch (m_methodsLfmin) {
            case igmUprime:
                U_e_prime = eprime * eps_e;
                gamma_min = (p - 2.) / (p - 1.) * U_e_prime / (n_prime * CGS::me * CGS::c * CGS::c); // Eq. A18 vanEarten+10 (and Sironi+13)
                break;
            case igmNakarPiran:
                gamma_min = (p - 2.) / (p - 1.) * (CGS::mp / CGS::me) * eps_e * EQS::Beta2(GammaSh); // Eq.(before 4.) in Nakar&Piran 1102.1020
                break;
            case igmJoh06:
                gamma_min = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (GammaSh - 1.) + 1.); // Eq. A3 in J+06
                break;
            case igmNumGamma:
                /// solve gamma_min fun numerically; use fixed limits and number of iterations
                gamma_min = Bisect(ElectronAndRadiaionBase::gammaMinFunc,
                                   1, 1e8, 0, .001, 100, this, status);
                /// If numerical solution failed, use simple analytical solution
                if (status < 0)
                    gamma_min = (p - 2.) / (p - 1.) * (eps_e * CGS::mp / CGS::me * (GammaSh - 1.) + 1.);
                break;
            case igmMAG21:
                /// downstream electron temperature:
                Theta = Margalit21::Theta_fun(beta, mu, mu_e, eps_t);
                gamma_min = Margalit21::gamma_m_fun(Theta);
                break;
        }

        /// check
        if (!std::isfinite(gamma_min)){
            (*p_log)(LOG_ERR,AT) << " error gm nan \n";
            exit(1);
        }

        /// limit the min lf to 1
        if ((lim_gm_to_1) && (gamma_min < 1.))
            gamma_min = 1.; // Sironi et al 2013 suggestion to limi gm=1. for Deep Newtinoan regime # TODO to be removed. I did no understand it

        /// calculate the (normalized) cooling Lorentz factor (eq. 18, MQ21): NOTE we use mean dynamical time:
        z_cool = (6.0 * M_PI * CGS::me * CGS::c / (CGS::sigmaT * B * B * t_e)) / Theta;

        if (!std::isfinite(z_cool)){
            (*p_log)(LOG_ERR,AT) << AT << " Theta = " << Theta << " z_cool="<<z_cool<<"\n";
            exit(1);
        }
    }
    /// for semi-neutonian regime, where gm ->
    void computeNonRelativisticFlattening(){
        switch (m_method_nonreldist) {
            case inone:
                break;
            case iuseGm:
                if (lim_gm_to_1) {
                    (*p_log)(LOG_ERR, AT) << " 'method_nonreldist' is incopatible with 'lim_gm_to_1' \n";
                    exit(1);
                }
                accel_frac = (std::pow(gamma_max, 2-p) - std::pow(gamma_min-1, 2-p))
                           / (std::pow(gamma_max,2-p) - 1.)
                           * (std::pow(gamma_max,1-p)-1.) / (std::pow(gamma_max,1-p)
                           - std::pow(gamma_min-1, 1-p));
                if (accel_frac > 1.)
                    accel_frac = 1.;
                break;
        }
    }

};


class ElectronAndRadiation : public ElectronAndRadiaionBase{

    /// --------------------------------------
    Source source{};
    State ele{}; // std::unique_ptr<State> ele = nullptr;
    State syn{}; //std::unique_ptr<State> syn = nullptr;
    State ssc{}; //std::unique_ptr<State> ssc = nullptr;
    SynKernel syn_kernel{};//std::unique_ptr<SSCKernel> ssc_kernel = nullptr;
    SSCKernel ssc_kernel{};//std::unique_ptr<SynKernel> syn_kernel = nullptr;
    double vol=-1.;
    double vol_p1=-1.;
    double dm=-1.;

    /// Analytical Synchrotron Sectrum; BPL;
    void checkEmssivityAbsorption(double em, double abs, double nuprime){
        if (( em < 0.) || (!std::isfinite( em )) ){

            (*p_log)(LOG_ERR,AT) << " em_pl_prime < 0 or nan ("<< em<<") or \n";
            (*p_log)(LOG_ERR,AT) << " abs_pl_prime < 0 or nan ("<< abs<<")\n";
            (*p_log)(LOG_ERR,AT) << " Error in data "
                                 << " eps_e = " << eps_e //<< "\n"
                                 << " eps_t = " << eps_t //<< "\n"
                                 << " n_prime = " << n_prime //<< "\n"
                                 << " gamma_min = " << gamma_min //<< "\n"
                                 << " gamma_max = " << gamma_max //<< "\n"
                                 << " gamma_c = " << gamma_c //<< "\n"
                                 << " B = " << B //<< "\n"
                                 << " Theta = " << Theta //<< "\n"
                                 << " z_cool = " << z_cool //<< "\n"
                                 << " nuprime = " << nuprime << "\n";
            exit(1);
        }
    }

public: // ---------------- ANALYTIC -------------------------- //

    ElectronAndRadiation(int loglevel, bool _is_rs) : ElectronAndRadiaionBase(loglevel, _is_rs){}

    void setPars( StrDbMap & pars, StrStrMap & opts, size_t nr ){
        setBasePars(pars, opts);
        // ---
        switch (m_eleMethod) {
            case iShockEleAnalyt:
                allocateStorageForAnalyticSpectra(pars, nr);
                break;
            case iShockEleNum:
                allocateStorageForNumericSpectra(pars, nr);
                break;
        }
    }

    void allocateStorageForAnalyticSpectra(StrDbMap & pars, size_t nr){
        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e28, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,200, true);//pars.at("nfreq");

        (*p_log)(LOG_INFO,AT) << " allocating comoving spectrum array (fs) "
                              << " freqs="<<nfreq << " by radii=" << nr << " Spec. grid="
                              << nfreq * nr << "\n";

        /// allocate space for the comoving spectrum
        out_spectrum.resize(m_freq_arr.size() * nr );
        out_specturm_ssa.resize(m_freq_arr.size() * nr );

        allocateStorageForOutputSpectra(nr, freq1, freq2, nfreq);
    }

    /// evaluate the comoving emissivity and absorption (frequency dependent)
    void computeSynchrotronEmissivityAbsorptionAnalytic( double nuprime, double & em, double & abs ) {

        // TODO WARNING I did replace n_prime with ne is absorption, but this might not be correct!!!

        em=0., abs=0.;

        nn *= ksi_n;// TODO this is not included in normalization obs.flux

        // defined in Margalit+21 arXiv:2111.00012
        double delta = eps_e / eps_t;
        double ne_ = mu_e * nn;

        SynchrotronAnalytic syn_an =
                SynchrotronAnalytic(gamma_min,gamma_c,gamma_max,B,p,do_ssa,p_log);

        if (m_sychMethod == METHODS_SYNCH::iJOH06)
            syn_an.computeAnalyticSynchJOH06(em, abs, nuprime, nn);
        else if (m_sychMethod == METHODS_SYNCH::iWSPN99)
            syn_an.computeAnalyticSynchWSPN99(em, abs, nuprime, nn);
        else if (m_sychMethod == METHODS_SYNCH::iDER06)
            syn_an.computeAnalyticSynchDER06(em, abs, nuprime, nn);
        else if (m_sychMethod == METHODS_SYNCH::iMARG21)
            syn_an.computeAnalyticSynchMARG21(em, abs, nuprime, ne_,
                                              delta,Theta,z_cool,accel_frac);
        else{
            (*p_log)(LOG_ERR,AT)<<" analytic synchrotron method is not supported \n";
            exit(1);
        }

        checkEmssivityAbsorption(em, abs, nuprime);


//        em /= (r * r * dr);
//        abs /= (r * r * dr);



        /// in piece-wise EATS you also need to normalize by n_layers
//        if ((m_method_ne == METHOD_NE::iuseNe) && (n_layers > 0)){
//            em /= (double)n_layers;
//            abs /= (double)n_layers;
//        }

    }

    /// compute spectrum for all freqs and add it to the container
    void computeSynchrotronSpectrumAnalytic(size_t it){
        size_t nfreq = m_freq_arr.size();
        double em, abs;
        /// computeSynchrotronEmissivityAbsorptionAnalytic emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
            computeSynchrotronEmissivityAbsorptionAnalytic(m_freq_arr[ifreq], em, abs );
            out_spectrum[ifreq + nfreq * it] = em;
            out_specturm_ssa[ifreq + nfreq * it] = abs;
            /// compute emissivity density
            if (m_method_ne == METHOD_NE::iuseNe){
                out_spectrum[ifreq + nfreq * it] /= n_protons;// /= (r * r * dr);///= (r * r * dr);
                out_specturm_ssa[ifreq + nfreq * it] /= n_protons;// /= (r * r * dr);//  /= n_protons;///= (r * r * dr);
            }
        }
    }

public: // -------------------- NUMERIC -------------------------------- //

    void allocateStorageForNumericSpectra(StrDbMap & pars, size_t nr){

        /// get freq. boundaries for calculation of the comoving spectrum
        double freq1 = getDoublePar("freq1", pars, AT, p_log,1.e7, true);//pars.at("freq1");
        double freq2 = getDoublePar("freq2", pars, AT, p_log,1.e28, true);//pars.at("freq2");
        size_t nfreq = (size_t)getDoublePar("nfreq", pars, AT, p_log,200, true);//pars.at("nfreq");

        double gam1 = getDoublePar("gam1", pars, AT, p_log,1, true);//pars.at("freq1");
        double gam2 = getDoublePar("gam2", pars, AT, p_log,1.e8, true);//pars.at("freq2");
        size_t ngams = (size_t)getDoublePar("ngam", pars, AT, p_log,250, true);//pars.at("nfreq");

        ele.allocate(gam1, gam2, ngams);
        ele.build_grid_chang_cooper(); // build the grid for the implicit solver

        syn.allocate(freq1*CGS::freqToErg, freq2*CGS::freqToErg, nfreq);

        if ((m_sychMethod != METHODS_SYNCH::iGSL) and (m_sychMethod != METHODS_SYNCH::iDER06)){
            (*p_log)(LOG_ERR,AT) << " only GSL or DER06 synchrotron options are "
                                    "avaialble when evolving electrons numerically. \n";
            exit(1);
        }

        /// allocate memory for synchrotron kernel
        if (m_sychMethod == METHODS_SYNCH::iGSL)
            syn_kernel.allocate(ele, syn, SynKernel::synchGSL);
        else
            syn_kernel.allocate(ele, syn, SynKernel::synchDermer);

        /// allocate memory for SSC structure and kernel
        if (m_methods_ssc == METHOD_SSC::iNumSSC) {
            ssc.allocate(freq1 * CGS::freqToErg, freq2 * CGS::freqToErg, nfreq);
            ssc_kernel.allocate(ele, syn, ssc, SSCKernel::sscNava);
            ssc_kernel.evalSSCkernel(ele, syn, ssc);
            ssc_kernel.evalSSCgridIntergral(ele, syn, ssc);
        }

        allocateStorageForOutputSpectra(nr, freq1, freq2, nfreq);

    }

    void evaluateElectronDistributionNumeric(double dt, double d_m, double r, double dr, double rp1, double drp1){

        /// check parameters
        if (ele.e.size() < 1){
            (*p_log)(LOG_ERR,AT) << " electon grid is not initialized\n";
            exit(1);
        }
        if (B <= 0 || gamma_min <= 0 || gamma_max <= 0){
            (*p_log)(LOG_ERR,AT) << "wrong value, cannot evolve electrons: "
                                    "B="<<B<<" gm="<<gamma_min<<" gM="<<gamma_max<<"\n";
            exit(1);
        }

        /// for adiabatic cooling of electron distribution
        vol = r * r * dr;
        vol_p1 = rp1 * rp1 * drp1;
        dm = d_m;
        double dlnVdt = 1. / dt * (1. - vol / vol_p1);

        /// number of injected electrons
        double N = dm / CGS::mp / dt;

        /// compute substepping to account for electron cooling timescales
        double dlog_gam = ((ele.e[1]-ele.e[0])/ele.e[0]);
        double delta_t_syn = CGS::sigmaT * gamma_max * B * B / (6. * M_PI * CGS::me * CGS::c);
        double delta_t_adi = (1.-gamma_max*gamma_max)/(3.*gamma_max*gamma_max)*dlnVdt;
        double delta_t = dlog_gam / (delta_t_syn + delta_t_adi);

        /// if cooling is too slow, we still need to evolve distribution
        if (delta_t <= dt)
            delta_t = dt;

        /// limit the max number of substeps
        else if (delta_t > dt / (double)max_substeps)
            delta_t = dt / (double) max_substeps;

        auto n_substeps = (size_t) (dt / delta_t); // number of substeps for electron evolution

        delta_t = std::max(delta_t, dt / (double)max_substeps);

        /// update source properties
        source.B = B;
        source.dlnVdt = dlnVdt;
        source.r = r;
        source.dr = dr;
        source.N = N;

        ChangCooper model = ChangCooper(source, ele, syn, ssc, syn_kernel, ssc_kernel);

        /// Assume source/escape do not change during substepping
        model.setSourceFunction(gamma_min, gamma_max, -p, N);
        model.setEscapeFunction(gamma_min, gamma_max);

        /// Init photon field
        Vector photon_field (syn.numbins, 0);

        for (size_t i = 0; i < (size_t)n_substeps; i++){

            /// Set gamma_dot terms (uses syn.n[] for photon density for SSC from previous solution)
            model.setHeatCoolTerms(2);

            /// Assuming sources change during substepping
//            model.setSourceFunction(m_data[Q::igm][it], m_data[Q::igM][it], index, N);
//            model.setEscapeFunction(m_data[Q::igm][it], m_data[Q::igM][it]);

            /// setup grid of differences for escape term
            model.computeDeltaJ();

            /// setup tridiagonal solver arrays
            model.setupVectors(delta_t);

            /// solve system
            model.solve(delta_t);

            /// Update photon field during electron evolution
//            model.update_radiation(); // Too long, if inside

        }

        /// Update photon field during electron evolution
        model.update_radiation(); // saves photon density in syn.n[] -> used in ssc cooling next iteration

    }

    void computeSynchrotronSpectrumNumeric(size_t it){

//        m_sychMethod = METHODS_SYNCH::iDER06;

        size_t nfreq = m_freq_arr.size();
        /// store emissivity and absorption for each frequency
        for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {

//            computeSynchrotronEmissivityAbsorptionAnalytic(m_freq_arr[ifreq]);
//
//            out_spectrum[ifreq + nfreq * it] = syn.j[ifreq] / n_or_nprime * n_prime; // M2 / mp * (4. * Gamma * rho_ism)
            out_spectrum[ifreq + nfreq * it] = syn.j[ifreq]; // M2 / mp * (4. * Gamma * rho_ism)
//            out_specturm_ssa[ifreq + nfreq * it] = syn.a[ifreq] / n_or_nprime * n_prime;
            out_specturm_ssa[ifreq + nfreq * it] = syn.a[ifreq];

            if (m_method_ne == METHOD_NE::iuseNe){
                out_spectrum[ifreq + nfreq * it] /= n_protons; //;/= (r * r * dr);
                out_specturm_ssa[ifreq + nfreq * it] /= n_protons; // /= (r * r * dr);
            }
//            std::cout<<ifreq<<" em="<<em / (r * r * dr)<<" j="<<out_spectrum[ifreq + nfreq * it]<<"\n";
//            std::cout<<ifreq<<" em="<<abs / (r * r * dr)<<" a="<<out_specturm_ssa[ifreq + nfreq * it]<<"\n";
//            int x = 1;

        }

        /// implicitely assume that SSC and Syn grids are the same. TODO generalize
        if (m_methods_ssc != METHOD_SSC::inoSSC)
            for (size_t ifreq = 0; ifreq < m_freq_arr.size(); ++ifreq) {
                out_spectrum[ifreq + nfreq * it] += ssc.j[ifreq];
            }
    }

public: // ---------------------- EATS -------------------------------- //

    double fixMe(double & em_prime, double & abs_prime, double Gamma, double GammaSh,
                 double mu, double r, double dr, double n_prime, double ne){
//    double flux_dens = 0.;
//        auto * p_pars = (struct Pars *) params;

        double beta = EQS::Beta(Gamma);
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor

        /// convert to the laboratory frame
        double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
        double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

        double beta_shock;
        switch (method_shock_vel) {

            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> computeSynchrotronEmissivityAbsorptionAnalytic shock velocity
                beta_shock = EQS::Beta(GammaSh);//us / sqrt(1. + us * us);
                break;
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
        dr /= ashock;

        double dr_tau = EQS::shock_delta(r,GammaSh);
        dr_tau /= ashock;
        double dtau = optical_depth(abs_lab, dr_tau, mu, beta_shock);
        double intensity = computeIntensity(em_lab, dtau, method_tau);

        if ((intensity < 0) || (!std::isfinite(intensity))) {
            (*p_log)(LOG_ERR, AT) << "intensity = " << intensity << "\n";
            exit(1);
        }

        double flux_dens=0.;
        switch (m_method_ne) {
            case iusenprime:
                flux_dens = intensity * r * r * dr;
                break;
            case iuseNe:
                flux_dens = intensity * r * r * dr * n_prime;
                break;
        }

        if (flux_dens < 0 || !std::isfinite(flux_dens)) {
            (*p_log)(LOG_ERR, AT) << "flux_dens_rs = " << flux_dens << "\n";
            exit(1);
        }

        return flux_dens;
#if 0
        switch (p_pars->m_method_eats) {
        case EjectaID2::iadaptive:
            switch (p_pars->m_method_rad) {
                // fluxDensAdaptiveWithComov()
                case icomovspec:
                    switch (p_pars->p_syn_a->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:
                    switch (p_pars->p_syn_a->m_method_ne){
                        case iusenprime:
                            flux_dens = (intensity * r * r * dr);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr / ne * n_prime;
                            break;
                    }
                    break;
            }
            break;
        case EjectaID2::ipiecewise:
            switch (p_pars->m_method_rad) {
                case icomovspec:
                    switch (p_pars->p_syn_a->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:

                    switch (p_pars->p_syn_a->m_method_ne){
                        case iusenprime:
                            flux_dens = (intensity * r * r * dr);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr / ne * n_prime;
                            break;
                    }
                    break;
            }
            break;
    }
#endif
//        return flux_dens;
    }

    double fluxDens(size_t ia, size_t ib, double freq, double Gamma, double GammaSh,
                    double mu, double r, double dr, double n_prime, double ne, Vector & r_arr){


        if (freq >= m_freq_arr[m_freq_arr.size()-1]){
            (*p_log)(LOG_ERR,AT)<<" freq_prime="<<freq
                                <<" > grid freq[-1]="<<m_freq_arr[m_freq_arr.size()-1]
                                <<" increase 'freq2' parameter\n";
            exit(1);
        }
        if (freq <= m_freq_arr[0]){
            (*p_log)(LOG_ERR,AT)<<" freq_prime="<<freq
                                <<" < grid freq[0]="<<m_freq_arr[0]
                                <<" decrease 'freq1' parameter\n";
            exit(1);
        }

        size_t ia_nu = findIndex(freq, m_freq_arr, m_freq_arr.size());
        size_t ib_nu = ia_nu + 1;

        /// interpolate the emissivity and absorption coefficines
        Interp2d int_em(m_freq_arr, r_arr, out_spectrum);
        double em_prime = int_em.InterpolateBilinear(freq, r, ia_nu, ib_nu, ia, ib);

        Interp2d int_abs(m_freq_arr, r_arr, out_specturm_ssa);
        double abs_prime = int_abs.InterpolateBilinear(freq, r, ia_nu, ib_nu, ia, ib);

        /// compute flux density
        return fixMe(em_prime,abs_prime, Gamma,GammaSh,mu,r,dr,n_prime,ne);
    }

};


// ---------------
enum METHOD_PWN_SPEC { inumBPL, ianaBPL };
static void initialize_e_dis(double *gam, double *dgam, double *dN_dgam, double *dN_dgam_dt,
                             double *dgam_dt, double *tad, double *tsyn, double gam_max, int Nbin_e){
    double dln_gam = log(gam_max)/(double)(Nbin_e-1);
    int i;
    for (i=0;i<Nbin_e;i++){
        gam[i] = exp(dln_gam*i);
        dgam[i] = gam[i]*(exp(dln_gam)-1.);
        dN_dgam[i] = 0.;
        dN_dgam_dt[i] = 0.;
        dgam_dt[i] = 0.;
        tad[i] = 0.;
        tsyn[i] = 0.;
    }
}
static void initialize_ph_dis(Vector & freqs, double *gam_ph, double *P_nu_syn, double *alpha_nu_syn){
    int i;
//    double del_ln_gam_ph = (log(gam_ph_max)-log(gam_ph_min))/(double)(Nbin_ph-1);
    for (i=0;i<freqs.size();i++){
        gam_ph[i] = freqs[i] / (CGS::MeC2 / CGS::H); // Hz->erg //gam_ph_min*exp(del_ln_gam_ph*(double)i);
        P_nu_syn[i] = 0.;
        alpha_nu_syn[i] = 0.;
    }
}
static double syn_func_fit(double x)
{
    /* analytical fitting of synchrotron function F(x) */
    /* see http://arxiv.org/pdf/1301.6908.pdf */

    double F1 = M_PI*pow(2.0,5.0/3.0)/sqrt(3.0)/GAMMA13*pow(x,1.0/3.0);
    double F2 = sqrt(M_PI/2.0)*exp(-x)*pow(x,1.0/2.0);

    double a1_1 = -0.97947838884478688;
    double a1_2 = -0.83333239129525072;
    double a1_3 = 0.1554179602681624;
    double H_1 = a1_1*pow(x,1.0)+a1_2*pow(x,1.0/2.0)+a1_3*pow(x,1.0/3.0);
    double delta_1 = exp(H_1);

    double a2_1 = -0.0469247165562628882;
    double a2_2 = -0.70055018056462881;
    double a2_3 = 0.0103876297841949544;

    double H_2 = a2_1*pow(x,1.0)+a2_2*pow(x,1.0/2.0)+a2_3*pow(x,1.0/3.0);
    double delta_2 = 1.0-exp(H_2);

    return F1*delta_1+F2*delta_2;
}
static void calc_syn_spec(const double B, const double dr, const double vol,
                          const double *gam, const double *dgam, const double *dN_dgam,
                          const double gam_max, const int Nbin_e, const double *gam_ph,
                          double *P_nu_syn, double *alpha_nu_syn, const int Nbin_ph){

//    double nu=0.,x=0.,sin_alpha=2./3.,tau_sa=0.;
//    double integ=0.,integ_alpha=0.;
//#pragma omp parallel for default(shared) reduction (P_nu_syn,alpha_nu_syn) firstprivate(nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol)// private(nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol) shared(std::cout,std::cerr,B,r,dr,gam,dgam,dN_dgam,gam_max,Nbin_e,Nbin_ph,gam_ph,P_nu_syn,alpha_nu_syn,CGS::c,CGS::MeC2,CGS::H,CGS::ELEC,CGS::me) num_threads( 6 ) // private(i,nu,x,sin_alpha,tau_sa,integ,integ_alpha,vol) shared(B,r,dr,gam,dgam,dN_dgam,gam_max,Nbin_e,gam_ph,P_nu_syn,alpha_nu_syn,Nbin_ph) num_threads( 6 )
#pragma omp parallel for num_threads( 6 )
    for (size_t k=0;k<Nbin_ph;k++) {
        double integ = 0.0;
        double sin_alpha=2./3.;
        double integ_alpha = 0.0;
        double nu = gam_ph[k]*CGS::MeC2/CGS::H;
//#pragma omp parallel for firstprivate(integ,integ_alpha,nu,x) num_threads( 10 )
        for (size_t i=0;i<Nbin_e;i++) {
            double x = (2.0*M_PI*nu)/(3.0*CGS::ELEC*gam[i]*gam[i]*B/2.0/CGS::me/CGS::c*sin_alpha); /* Eq. (6.17c) of Rybicki & Lightman */

            //            std::cout << x << ", ";
            if ((i==0) || (i==Nbin_e-1)) {
                integ += 0.5*dN_dgam[i]*dgam[i]*syn_func_fit(x);
                integ_alpha += -0.5*sin_alpha*std::pow(gam[i],2.0)
                               * (-dN_dgam[i]/std::pow(gam[i],2.0))
                               / dgam[i]*syn_func_fit(x)
                               * dgam[i] / CGS::MeC2;
            }
            else {
                integ += dN_dgam[i]*dgam[i]*syn_func_fit(x);
                integ_alpha += -sin_alpha*std::pow(gam[i],2.0)*
                               (dN_dgam[i+1]/std::pow(gam[i+1],2.0) - dN_dgam[i]/std::pow(gam[i],2.0))
                               / dgam[i] * syn_func_fit(x) * dgam[i] / CGS::MeC2;
            }
//                if ((!std::isfinite(P_nu_syn[k]))||(P_nu_syn[k]<0)) {
//                    std::cerr << AT << " nan in pwn spec P_nu_syn[i]=" << P_nu_syn[k] << "\n";
//                    std::cout << " x=" << x << " integ=" << integ << " nu=" << nu << " dN_dgam[0]=" << dN_dgam[0]
//                              << " dgam[0]=" << dgam[0] << " gam[0]=" << gam[0] << " B=" << B
//                              << " vol=" << vol << " tau_sa=" << tau_sa << " alpha_nu_sym[0]=" << alpha_nu_syn[0] << "\n";
//                    exit(1);
//                }
        }

        P_nu_syn[k] = sqrt(3.0)*pow(CGS::ELEC,3.0)*B*sin_alpha/CGS::MeC2*integ; /* Eq. (6.33) x (2 pi) of Rybicki & Lightman */
        alpha_nu_syn[k] = CGS::c*CGS::c/8.0/M_PI/nu/nu*sqrt(3.0)*pow(CGS::ELEC,3.0)*B*integ_alpha/vol; /* Eq. (6.52) of Rybicki & Lightman */

        double tau_sa = alpha_nu_syn[k] * dr;
        if (tau_sa > 1.0e-6){
            P_nu_syn[k] = (1.0-exp(-tau_sa)) * P_nu_syn[k] / tau_sa;
        }

        integ = 0.0;
        integ_alpha = 0.0;


    }
}
static double Ntot(double *dgam, double *dN_dgam, int Nbin_e) {
    int i;
    double tmp=0.;
    for (i=0;i<Nbin_e-1;i++){
        tmp += dN_dgam[i]*dgam[i];
    }
    tmp += 0.5*dN_dgam[Nbin_e-1]*dgam[Nbin_e-1];
    return tmp;
}
static double dgam_dt_ad(double gam, double t) {
    return gam/t;
}
static double dgam_dt_syn(double gam, double B) {
    // electron synchrotron energy loss rate (see e.g., Eq. 7.13 of Dermer & Menon)
    // double sin2phi = 2.0/3.0; /* averaging pitch angle */
    // double beta_par = 1.0; /* assuming that particles are relativistic */

    return 4.0/3.0*CGS::c*CGS::SIGMA_T*(B*B/8.0/M_PI)*gam*gam/CGS::MeC2;
}
static void cooling(double t, double B, double *dgam_dt, double *gam, double *tad, double *tsyn, int Nbin_e){
    int i;
    for (i=0;i<Nbin_e;i++) {
        dgam_dt[i] = dgam_dt_ad(gam[i],t)+dgam_dt_syn(gam[i],B);
        tad[i] = gam[i]/dgam_dt_ad(gam[i],t);
        tsyn[i]= gam[i]/dgam_dt_syn(gam[i],B);
    }
}
static double dN_dgam_dt_inj(double gam, double Lpsr, double epse, double gam_b, double gam_max, double p1, double p2){
    double fac_bol = 1./(2. - p1) * (1. - pow(gam_b, -2. + p1))
                   + 1./(2. - p2) * (pow(gam_max/gam_b, 2. - p2) - 1.);
    if (gam < gam_b){
        return epse * Lpsr/(MeC2*gam_b*gam_b) / fac_bol * pow(gam/gam_b,-p1);
    }
    else {
        return epse * Lpsr/(MeC2*gam_b*gam_b) / fac_bol * pow(gam/gam_b,-p2);
    }
}
static void time_evolution_e(double dt, double *gam, double *dgam, double *dN_dgam,
                             double *dN_dgam_dt, double *dgam_dt, int Nbin_e)
{
    int i;
    double dN_dgam_old[Nbin_e];
    double N_cool = 0.0;
    double dN_dgam_cool = 0.0;

    for (i=0; i<Nbin_e; i++){
        dN_dgam_old[i] = dN_dgam[i];
    }

    dN_dgam[Nbin_e-1] = (dN_dgam_old[Nbin_e-1])
                      / (1.0+dt/dgam[Nbin_e-1]*dgam_dt[Nbin_e-1])
                      + dN_dgam_dt[Nbin_e-1]*dt;
    for(i = Nbin_e-2; i>0; i--){
        dN_dgam[i] = (dN_dgam_old[i]+dN_dgam[i+1]*dt/dgam[i]*dgam_dt[i+1]) / (1.0+dt/dgam[i]*dgam_dt[i])
                   + dN_dgam_dt[i] * dt;
    }
    dN_dgam[0] = dN_dgam_old[0]
               + dN_dgam_dt[1] * dt / dgam[0] * dgam_dt[1]
               + dN_dgam_dt[0] * dt;
}
static void injection(double *gam, double *dgam, double *dN_dgam_dt, double *dgam_dt, double Lpsr, double dt,
                      double *N_inj_tot, double epse, double gam_b, double gam_max,
                      double p1, double p2, int Nbin_e){
    int i,j;
    double tmp=0.0,total_frac_depo=0.0,frac_resi=0.0,frac_depo=0.0;
    double dN_dgam_dt_eff[Nbin_e];

    /// standard integration
    for (i=0;i<Nbin_e;i++){
        dN_dgam_dt[i] = dN_dgam_dt_inj(gam[i],Lpsr,epse,gam_b,gam_max,p1,p2);
        dN_dgam_dt_eff[i] = 0.0;
        tmp += dN_dgam_dt[i] * dt * dgam[i];
    }
    tmp -= 0.5*(dN_dgam_dt[0] * dt * dgam[0] + dN_dgam_dt[Nbin_e-1 ]* dt * dgam[Nbin_e-1]);
    *N_inj_tot += tmp;

    /// ???
    for (i=Nbin_e-1; i>0; i--){
        total_frac_depo = 0.0;
        frac_resi = 1.0;
        frac_depo = 0.0;
        for (j=i; j>0; j--){
            frac_depo = frac_resi/(1.0 + dt/dgam[i] * dgam_dt[i]);
            dN_dgam_dt_eff[j] += frac_depo * dN_dgam_dt[i] * (dgam[i]/dgam[j]); // eff?

            total_frac_depo += frac_depo;
            if(total_frac_depo > 1.0)
                break;

            frac_resi = (1.-total_frac_depo);
        }
        dN_dgam_dt_eff[0] += frac_resi*dN_dgam_dt[i]*(dgam[i]/dgam[0]);
    }
    dN_dgam_dt_eff[0] += dN_dgam_dt[0];

    for (i=0;i<Nbin_e;i++){
        dN_dgam_dt[i] = dN_dgam_dt_eff[i];
    }

}
// ------------------
class MagnetarSynchrotron{
    struct Pars{
        double gam_b{};
        double p1{};
        double p2{};
        double eps_e{};
        double gam_max{};
//        int Nbin_e{};
//        int Nbin_ph{};
    };
    std::unique_ptr<Pars> p_pars = nullptr;
    std::unique_ptr<logger> p_log = nullptr;
//    VecVector & m_data;
//    Vector spectrum{}; Vector spec_gams{}; Vector spec_freqs{};
//    Vector emissivity{}; Vector absorption{};
public:
    MagnetarSynchrotron(unsigned loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "MagnetarSynchrotron");
        p_pars = std::make_unique<Pars>();
    }
    std::unique_ptr<Pars> & getPars(){ return p_pars; }
//    Vector & getSpectrum(){return spectrum;}
//    Vector & getEmissivity(){return emissivity;}
//    Vector & getAbsorption(){return absorption;}

    void setPars(StrDbMap & pars, StrStrMap & opts){
        p_pars->gam_b = getDoublePar("gam_b",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->p1 = getDoublePar("p1",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->p2 = getDoublePar("p2",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->eps_e = getDoublePar("eps_e_spec",pars,AT,p_log,1., true);//pars.at("ksi_n");
        p_pars->gam_max = getDoublePar("gam_max",pars,AT,p_log,1.e6, false);//pars.at("ksi_n");
        /// ---
//        spec_gams = makeVecFromString(getStrOpt("spec_gams",opts,AT,p_log,"",true),p_log);
//        spec_freqs = makeVecFromString(getStrOpt("spec_freqs",opts,AT,p_log,"",true),p_log);
//        p_pars->Nbin_e = (int)spec_gams.size();
//        p_pars->Nbin_ph = (int)spec_freqs.size();

    }

    void computeComovingSpectrum(Vector & spectrum, Vector & emissivity, Vector & absorption,
                                 Vector & times, Vector & freqs, Vector & gams, Vector & Bnb,
                                 Vector & vol, Vector & drnb, Vector & Lpsr){

        if ((spectrum.size() < 1)||(emissivity.empty() < 1)||(absorption.empty() < 1)){
            spectrum.resize(times.size()*freqs.size(), 0.);
            emissivity.resize(times.size()*freqs.size(), 0.);
            absorption.resize(times.size()*freqs.size(), 0.);
        }
        /// initialize storage
        int Nbin_e = (int)gams.size();
        double gam[Nbin_e],dN_dgam[Nbin_e],dgam[Nbin_e],dN_dgam_dt[Nbin_e],dgam_dt[Nbin_e],tad[Nbin_e],tsyn[Nbin_e];
        double N_inj_tot=0.;
        initialize_e_dis(gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,tad,tsyn,p_pars->gam_max,Nbin_e);
        int Nbin_ph = (int)freqs.size();
        double gam_ph[Nbin_ph],P_nu_syn[Nbin_ph],alpha_nu_syn[Nbin_ph];
        initialize_ph_dis(freqs,gam_ph,P_nu_syn,alpha_nu_syn);
        /// ----------------
        size_t ii = 0;
        for (size_t it = 0; it < times.size(); it++) {
//            double Bnb = 0., rnb = 0., dr = 0., Lpsr = 0., dt = 0.;
            double t = times[it];
            double dt = (times[it] - times[it - 1]);
            if (it > 1)
                calc_syn_spec(Bnb[it],drnb[it],vol[it],gam,dgam,dN_dgam,p_pars->gam_max,
                              Nbin_e,gam_ph,P_nu_syn,alpha_nu_syn,Nbin_ph);
            double number_conservation = Ntot(dgam,dN_dgam,Nbin_e)/N_inj_tot;
            if (number_conservation > 1.2 or number_conservation < 0.9){
                std::cerr << AT << " conserv="<<number_conservation<<"\n";
            }
            cooling(t,Bnb[it],dgam_dt,gam,tad,tsyn,Nbin_e);
            injection(gam, dgam, dN_dgam_dt, dgam_dt, Lpsr[it], dt, &N_inj_tot,
                      p_pars->eps_e, p_pars->gam_b, p_pars->gam_max, p_pars->p1, p_pars->p2, Nbin_e);
            time_evolution_e(dt,gam,dgam,dN_dgam,dN_dgam_dt,dgam_dt,Nbin_e);
            /// --------------------------
            if ((!std::isfinite(P_nu_syn[0]))||(!std::isfinite(P_nu_syn[Nbin_ph-1])) || (P_nu_syn[0] < 0)||(P_nu_syn[Nbin_ph-1] < 0)){
                (*p_log)(LOG_ERR,AT) << " nan in pwn spec P_nu_syn[0]=" << P_nu_syn[0] << "\n";
                std::cout << "Bnb="<<Bnb<<" gam[0]="<<gam[0]<<" dgam[0]="<<dgam[0]<<" dN_dgam[0]="<<dN_dgam[0]<<" gam_max="<<p_pars->gam_max
                          <<" Nbin_e="<<Nbin_e<<" Lpsr="<<Lpsr<<" n_ing_tot="<<N_inj_tot<<" number_conserv="<<number_conservation<<"\n";
                exit(1);
            }
//                spec[ii] = gam_ph[inu] * CGS::MeC2 / CGS::H * P_nu_syn[inu]; // [erg/s/Hz]
            for (size_t inu = 0; inu < freqs.size();  inu++) {
                emissivity[inu + freqs.size() * it] = P_nu_syn[inu]; // TODO ? IS IT TRUE?
                absorption[inu + freqs.size() * it] = alpha_nu_syn[inu];
                spectrum[inu + freqs.size() * it] = gam_ph[inu] * CGS::MeC2 / CGS::H * P_nu_syn[inu]; // [erg/s/Hz]

                double abs_fac=0.;
                double tau_sa = absorption[inu + freqs.size() * it] * drnb[it];
                if(tau_sa == 0.0)
                    abs_fac = 1.0;
                else if(tau_sa > 0.0)
                    abs_fac = -expm1(-tau_sa) / tau_sa;
                else {
                    abs_fac = expm1(tau_sa) / tau_sa; //* exp(
                    //abs * DR * beta_shock*mu / (mu - beta_shock));
                }
                spectrum[inu + freqs.size() * it] = tau_sa > 1e-6 ? emissivity[inu + freqs.size() * it] * abs_fac : emissivity[inu + freqs.size() * it];

                ii++;
            }
        }
    }

};


//    void allocateSpaceForComovingSpectrum(size_t nr){
//        /// allocate memory for comoving spectra to be evaluated
//        m_freq_grid = TOOLS::MakeLogspace(log10(p_pars->freq1),
//                                          log10(p_pars->freq2),(int)p_pars->nfreq);
//        if (p_pars->method_comp_mode == SynchrotronAnalytic::METHODS_RAD::icomovspec){
//            std::cerr << " allocating comoving spectrum array (fs) "
//                               << " freqs="<<m_freq_grid.size() << " by radii=" << nr << " Spec. grid="
//                               << m_freq_grid.size() * nr << "\n";
//            m_spectrum.resizeEachImage( m_freq_grid.size() * nr );
//        }
//    }

//class SynchrotronNumeric{
//    std::unique_ptr<logger> p_log;
//    Vector farr_sync;
//    Vector freq_ssc;
//    struct PWNPars{
//        double f1_sync, f2_sync; size_t nf_sync;
//        double f1_ssc, f2_ssc; size_t nf_ssc;
//
//    };
//    std::unique_ptr<PWNPars> p_pars;
//public:
//    SynchrotronNumeric(int loglevel){
//        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "SynchrotronNumeric");
//        p_pars = std::make_unique<PWNPars>();
//    }
//    void setMagPars(StrDbMap pars, StrStrMap sets){
////        p_pars->f1_ssc =
//    }
//};


//class SynchrotronAnalyticComoving{
//    enum METHODS_SYNCH { iWSPN99, iJOH06, iDER06, iMARG21 };
//    enum METHODS_LFMIN { igmUprime, igmNakarPiran, igmJoh06, igmMAG21 };
//    enum METHODS_SSA { iSSAoff, iSSAon };
//    enum METHOD_TAU { iTHICK, iSMOOTH, iSHARP };
//    enum QQ { i_em_pl, i_em_th, i_abs_pl, i_abs_th, i_em, i_abs };
//    struct PWNPars{
//        // --- in
//        double eps_e=-1, eps_b=-1, eps_t=-1, p=-1, ksi_n=-1;
//        double mu=-1, mu_e=-1;
//        bool lim_gm_to_1= true;
//        double beta_min = -1;
//        // --- methods
//        METHODS_SYNCH m_sychMethod{};
//        METHODS_LFMIN m_methodsLfmin{};
//
//        METHODS_SSA m_methods_ssa{};
//        QQ m_marg21opt_em = i_em_pl;
//        QQ m_marg21opt_abs = i_abs_pl;
//    };
//    PWNPars * p_pars = nullptr;
//
//    static constexpr size_t n_vars = 5;
//    enum QS { igm, igc };
//    std::vector<std::string> m_vars{ "gm", "gc", "num", "nuc", "pmax" };
//    VecArray m_data{};// ( n_vars );
//public:
//    SynchrotronAnalyticComoving(size_t nt){
//        p_pars = new PWNPars();
//        allocateSpace(nt);
//    }
//    ~SynchrotronAnalyticComoving(){ delete p_pars; }
//    void allocateSpace(size_t nt){
//        m_data.resizeEachImage(n_vars);
//        for (auto & arr : m_data){
//            arr.resizeEachImage( nt, 0. );
//        }
//    }
//
//    void setMagPars(std::unordered_map<std::string,double> & pars, std::unordered_map<std::string,std::string> & opts){
//        // set parameters
//
//    }
//};

#endif //SRC_SHOCK_MICROPHYSICS_H
