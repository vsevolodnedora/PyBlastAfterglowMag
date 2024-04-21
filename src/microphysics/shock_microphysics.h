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

enum METHOD_NONRELDIST{ inone, iuseGm, iuseAyache, iuseSironi };

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
double computeIntensity(const double em_lab, const double dtau, METHOD_TAU m_tau) {

    if (dtau == 0)
        return em_lab;

    /// thick: Special case: use the optically thick limit *everywhere*
    double intensity_thick = dtau > 0 ? em_lab / dtau : 0.;

    /// from Dermer+09 book
    double u = 1. / 2. + exp(-dtau) / dtau - (1 - exp(-dtau)) / (dtau * dtau);
    double intensity_approx = dtau > 1.e-3 ? em_lab * (3. * u / dtau) : em_lab;
//    if (dtau > 1)
//        int z = 1;
    // correction factor to emissivity from absorption
    // ( 1 - e^(-tau) ) / tau  (on front face)
    // back face has extra factor ~ e^-beta_shock/(mu-beta_shock)
    // for now ignoring shadowing by the front face.
    double abs_fac;
    if (dtau > 0.0) abs_fac = -std::expm1(-dtau) / dtau; // exp(x) - 1.
    else abs_fac = std::expm1(dtau) / dtau; //* exp(//abs * DR * beta_shock*mu / (mu - beta_shock));
    double intensity_smooth = em_lab * abs_fac;

    // Apply self-absorption "simply".
    // Compute flux in optically thick limit,
    // use in final result if less than optically thin calculation.
    // e.g. use tau->infty limit if tau > 1.0
    double intensity_sharp = em_lab;
    if (dtau > 1.0) intensity_sharp = em_lab / dtau;   // "Forward" face
    else if (dtau < -1.0) intensity_sharp = 0.0;  // "Back" face --> assume shadowed by front

    double intensity = em_lab;
    switch (m_tau) {

        case iAPPROX:
            intensity = intensity_approx;
            break;
        case iTHICK:
            intensity = intensity_thick;
            break;
        case iSMOOTH:
            intensity = intensity_smooth;
            break;
        case iSHARP:
            intensity = intensity_sharp;
            break;
    }
    return intensity;
}
double computeIntensity_OLD(const double em_lab, double dtau, METHOD_TAU m_tau){

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
            abs_fac = -std::expm1(-dtau) / dtau; // exp(x) - 1.
        else {
            abs_fac = std::expm1(dtau) / dtau; //* exp(//abs * DR * beta_shock*mu / (mu - beta_shock));
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
    double eprime=-1.,Gamma=-1,GammaSh=-1, betaSh=-1.,t_e=-1.;
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
    void updateSockProperties(double e_prime, double Gamma_, double Gamma_shock, double t_e_, double tcomov_,
                              double n_prime_, double n_protons_){
        /// Store current parameters

        GammaSh = Gamma_shock; // for bisect solver
        eprime = e_prime; // for bisect solver
        Gamma = Gamma_; // for bisect solver
        betaSh = EQS::BetaFromGamma(Gamma_shock); // for bisect solver
        t_e = t_e_; // for bisect solver
        tcomov = tcomov_; // for bisect solver
        n_prime = n_prime_; //  for bisect solver
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
    void evaluateElectronDistributionAnalytic() {

        if (not checkParams())
            return;
        /// compute comoving magnetic field
        computeMagneticField();
        /// compute limits of electron spectrum
        computeGammaMax();
        computeGammaMin();
        computeGammaCool();
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
protected:
    int m_loglevel = -1;
    std::unique_ptr<logger> p_log = nullptr;
    bool is_rs = false;
    double tcomov0=-1.;
    double tcomov=-1;
    /// --------------------------------------
    double eps_e=-1, eps_b=-1, eps_t=-1, p=-1;// ksi_n=-1;
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
    void computeMagneticField(){
        /// fraction of the shock energy density in electrons
        double U_b_prime = eprime * eps_b;
        double B1 = sqrt(8. * CGS::pi * U_b_prime);
        /// See Nava 2013 paper or others; classical equation
        double B2 = std::sqrt( 32. * M_PI * eps_b * (GammaSh - 1.)
                       * (GammaSh + 3. / 4) * n_prime * CGS::mp * CGS::c2);
        double B3 = sqrt(9.0 * M_PI * eps_b * n_prime * mu * CGS::mp)
            * (betaSh * CGS::c);

        switch (m_methodsB) {
            case iBuseUb:       B=B1; break;
            case iBuseGammaSh:  B=B2; break;
            case iBasMAG21:     B=B3; break;
        }
    }
    void computeGammaMax(){
        /// get electron distribution boundaries (has ot be set BEFORE gamma_min)
        double gamma_max_1 = sqrt(6.0 * CGS::pi * CGS::qe / CGS::sigmaT / B); // Kumar+14
        switch (m_methodsLfmax) {
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
    void computeGammaMin(){
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
        switch (m_methodsLfmin) {
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
        if ((gamma_min < 1.) || (status < 0)) {
            switch (m_method_nonreldist) {
                case inone:
                    accel_frac = 1;
                    break;
                case iuseGm:
                    accel_frac = gamma_min;
                    break;
                case iuseAyache:
                    accel_frac = (std::pow(gamma_max, 2 - p) - std::pow(gamma_min - 1., 2. - p))
                                 / (std::pow(gamma_max, 2. - p) - 1.)
                                 * (std::pow(gamma_max, 1. - p) - 1.) / (std::pow(gamma_max, 1. - p)
                                                                         - std::pow(gamma_min - 1., 1. - p));
                    break;
                case iuseSironi: // Sironi et al 2013
//                    accel_frac = (p - 2.0) / (p - 1.0) * eps_e * CGS::mp / CGS::me * (GammaSh - 1.0) / 1.;
                    accel_frac = (p - 2.0) / (p - 1.0) * eps_e * eprime / (n_prime * CGS::me * CGS::c * CGS::c);
                    break;
            }
            gamma_min = 1.;
            accel_frac = accel_frac > 1. ? 1 : accel_frac;
        }

        /// check
        if (!std::isfinite(gamma_min) || (!std::isfinite(accel_frac))){
            (*p_log)(LOG_ERR,AT) << " Wrong value gamma_min="<<gamma_min<<" accel_frac="<<accel_frac<<"\n";
            exit(1);
        }

//        printf("gamma_min to %.2e GammaSh=%.2e eprime=%.2e p=%.2f gamma_max=%.2e ", gamma_min, GammaSh, eprime, p, gamma_max);

        /// prevent gamma_max to go beyond the electron grid
        if (m_eleMethod != METHODS_SHOCK_ELE::iShockEleAnalyt)
            if (gamma_min < ele.e[0]) {
                gamma_min = ele.e[0];

        }
//        printf("gamma_min to %.2e GammaSh=%.2e eprime=%.2e p=%.2f gamma_max=%.2e ", gamma_min, GammaSh, eprime, p, gamma_max);

    }
    void computeGammaCool(){
        // Eq. A19 in vanEarten+10
        double gamma_c_1 = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * B * B) * t_e / Gamma;

        double gamma_c_2 = ((6. * CGS::pi * me * c) / (CGS::sigmaT * B * B * (tcomov))); // (tcomov0-1.e-6)

        switch (m_methodsLfcool) {
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


class ElectronAndRadiation : public ElectronAndRadiaionBase{

    Source source{};
    SynKernel & syn_kernel;//std::unique_ptr<SSCKernel> ssc_kernel = nullptr;
    SSCKernel & ssc_kernel;//std::unique_ptr<SynKernel> syn_kernel = nullptr;
    ChangCooper model = ChangCooper(source, ele, syn, ssc, syn_kernel, ssc_kernel);
    double theta_h = 0;
    double dr_comov;
    size_t n_substeps=0;

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

    ElectronAndRadiation(int loglevel, Vector & gam_grid, Vector & freq_grid,
                         SynKernel & syn_kernel, SSCKernel & ssc_kernel, bool _is_rs) :
        syn_kernel(syn_kernel), ssc_kernel(ssc_kernel),
        ElectronAndRadiaionBase(gam_grid, freq_grid, loglevel, _is_rs){
    }

    void setPars( StrDbMap & pars, StrStrMap & opts, size_t nr, size_t ish, size_t il,
                  double theta_h_, double tcomov0_ ){
        setBasePars(pars, opts);
        if (tcomov0_<=0.){
            (*p_log)(LOG_ERR,AT) << "Bad value: tcomov="<<tcomov0_<<"\n";
            exit(1);
        }
        tcomov0 = tcomov0_;
        theta_h = theta_h_;

        /// allocate space for spectra
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";
        else
            fs_or_rs += "_fs";
        freq1 = getDoublePar("freq1"+fs_or_rs, pars, AT, p_log, 1.e7, true);//pars.at("freq1");
        freq2 = getDoublePar("freq2"+fs_or_rs, pars, AT, p_log, 1.e28, true);//pars.at("freq2");
        nfreq = (size_t) getDoublePar("nfreq"+fs_or_rs, pars, AT, p_log, 200, true);//pars.at("nfreq");

        is_dense_output = getBoolOpt("save_spec", opts, AT, p_log, false, true);

        /// allocate memory for spectra to be used in EATS interpolation
        if (m_eleMethod==METHODS_SHOCK_ELE::iShockEleAnalyt)
            allocateForAnalyticSpectra(nr, ish, il);
        else
            allocateForNumericSpectra(pars, nr, ish, il);

    }
    void allocateForAnalyticSpectra(size_t nr, size_t ish, size_t il){
        /// allocate space for spectra
        (*p_log)(LOG_INFO, AT)
                << "[ish="<<ish<<" il="<<il<<"] AN "
                << " nfreq="<<nfreq
                << (m_methods_ssa==METHODS_SSA::iSSAon ? " SSA " : " ")
                << " ntimes="<<nr
                << " ALLOCATING\n";

        total_rad.allocate(freq1, freq2, nfreq, nr, true);// todo always dense
        syn.allocate(freq1, freq2, nfreq, nr, true);// todo add if statment if to save dense output
        if (m_methods_ssc!=METHOD_SSC::inoSSC)
            ssc.allocate(freq1, freq2, nfreq, nr, true);
    }


    /// evaluate the comoving emissivity and absorption (frequency dependent)
    void computeSynchrotronEmissivityAbsorptionAnalytic( double nuprime, double & em, double & abs ) {

        // TODO WARNING I did replace n_prime with ne is absorption, but this might not be correct!!!

        SynchrotronAnalytic syn_an =
                SynchrotronAnalytic(gamma_min, gamma_c, gamma_max, B, p,
                                    m_methods_ssa == METHODS_SSA::iSSAon,
                                    p_log);

        em=0., abs=0.;
        if (m_sychMethod == METHODS_SYNCH::iJOH06)
            syn_an.computeAnalyticSynchJOH06(em, abs, nuprime, nn*accel_frac);
        else if (m_sychMethod == METHODS_SYNCH::iWSPN99)
            syn_an.computeAnalyticSynchWSPN99(em, abs, nuprime, nn*accel_frac);
        else if (m_sychMethod == METHODS_SYNCH::iDER06)
            syn_an.computeAnalyticSynchDER06(em, abs, nuprime, nn*accel_frac);
        else if (m_sychMethod == METHODS_SYNCH::iMARG21) {
            // defined in Margalit+21 arXiv:2111.00012
            syn_an.computeAnalyticSynchMARG21(em, abs, nuprime, mu_e * nn,
                                              eps_e / eps_t, Theta, z_cool, accel_frac);
        }
        else{
            (*p_log)(LOG_ERR,AT)<<" analytic synchrotron method is not supported \n";
            exit(1);
        }

        checkEmssivityAbsorption(em, abs, nuprime);

        // get emissivity and absorption per particle (for EATS)
        if (m_method_ne == METHOD_NE::iuseNe){
            em /= n_protons; // /= (source.r * source.r * source.dr); // /= n_protons;
            abs = abs / n_protons * n_prime; // /= (source.r * source.r * source.dr);// /= n_protons;
        }
    }

    /// compute spectrum for all freqs and add it to the container
    void computeSynchrotronSpectrumAnalytic(size_t it) {
        for (size_t i = 0; i < syn.numbins; ++i) {
            computeSynchrotronEmissivityAbsorptionAnalytic(
                    syn.e[i], syn.j[i], syn.a[i]
            );
            syn.intensity[i] = computeIntensity(syn.j[i], syn.a[i] * dr_comov, METHOD_TAU::iAPPROX);
        }

        /// store the result
        syn.save_to_all(it);
        total_rad.add_to_all(syn, it);
        if (m_methods_ssc!=METHOD_SSC::inoSSC) {
            ssc.save_to_all(it);
            total_rad.add_to_all(ssc, it);
        }

    }

public: // -------------------- NUMERIC -------------------------------- //

    bool is_distribution_initialized = false; // if the analytic profile for the first iteration is set
    bool is_dense_output = false;

    /// Analytic electron spectrum, power-law
    void powerLawElectronDistributionAnalytic(double tcomov0_, double n_ele_inj, Vector & N_ele){
//        /// prevent gamma_max to go beyond the electron grid
//        if (gamma_max > 0.99 * ele.e[ele.numbins-1])
//            gamma_max = 0.99 * ele.e[ele.numbins-1];

//        double gamma_c_ = 6. * CGS::pi * CGS::me * CGS::c / (CGS::sigmaT * dt * B * B) / Gamma; // Eq. A19 in vanEarten+10
//        gamma_c = ((6. * CGS::pi * me * c) / (CGS::sigmaT * B * B * (tcomov0_ - (tcomov0-1.e-6))));
//        gamma_c = std::min(gamma_c, gamma_max);
        double K_s=0.,K_f=0.;
        if (gamma_min > gamma_c) { // fast cooling regime
            K_f = n_ele_inj / (
                    pow(gamma_min, 1. - p) * (1. / gamma_c - 1. / gamma_min)
                    + 1. / p * (pow(gamma_min, -p) - pow(gamma_max, -p))
                    );
            for (size_t i = 0; i < ele.numbins; i++) {
                double gamma = ele.e[i];
                double N = 0.0;
                if ((gamma_c <= gamma) && (gamma < gamma_min))
                    N_ele[i] = K_f * pow(gamma_min, 1. - p) * pow(gamma, -2.);

                else if ((gamma_min < gamma) && (gamma <= gamma_max))
                    N_ele[i] = K_f * pow(gamma, -p - 1.);

                if (!std::isfinite(N_ele[i]) || N_ele[i] < 0.){
                    (*p_log)(LOG_ERR,AT)<< " N_ele() = "<< N_ele[i]<< "\n";
                    exit(1);
                }
            }
        }
        else{//else if (gamma_c > gamma_min) { // slow cooling regime
            K_s = n_ele_inj /
                    (
                    (1. / (1. - p)) * (pow(gamma_c, 1. - p)
                    - pow(gamma_min, 1. - p))
                    - gamma_c / p * (pow(gamma_max, -p) - pow(gamma_c, -p))
                    );
            for (size_t i = 0; i < ele.numbins; ++i) {
                double gamma = ele.e[i];
                double N = 0.0;
                if ((gamma_min <= gamma) && (gamma < gamma_c))
                    N_ele[i] = K_s * pow(gamma, -p);

                else if ((gamma_c <= gamma) && (gamma <= gamma_max))
                    N_ele[i] = K_s * gamma_c * pow(gamma, -p - 1.);

                if (!std::isfinite(N_ele[i]) || N_ele[i] < 0.){
                    (*p_log)(LOG_ERR,AT)<< " N_ele() = "<< N_ele[i]<< "\n";
                    exit(1);
                }
            }
        }
//        int x = 0;
    }

    void allocateForNumericSpectra(StrDbMap & pars, size_t nr, size_t ish, size_t il){
        /// get freq. boundaries for calculation of the comoving spectrum
        /// allocate space for spectra
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";
        else
            fs_or_rs += "_fs";

        gam1 = getDoublePar("gam1"+fs_or_rs, pars, AT, p_log,1, true);//pars.at("freq1");
        gam2 = getDoublePar("gam2"+fs_or_rs, pars, AT, p_log,1.e8, true);//pars.at("freq2");
        ngam = (size_t)getDoublePar("ngam"+fs_or_rs, pars, AT, p_log,250, true);//pars.at("nfreq");

        (*p_log)(LOG_INFO, AT)
                << "[ish="<<ish<<" il="<<il<<"] NUM "
                << "ngam="<<ngam<<" nfreq="<<nfreq
                << (m_methods_ssa==METHODS_SSA::iSSAon ? " SSA " : " ")
                << (m_methods_ssc==METHOD_SSC::iNumSSC ? " SSC " : " ")
                << " ntimes="<<nr
                << " ALLOCATING\n";

        ele.allocate(gam1, gam2, ngam, nr, is_dense_output); // todo add if statement if dense
        ele.build_grid_chang_cooper(); // build the grid for the implicit solver
        syn.allocate(freq1, freq2, nfreq, nr, is_dense_output); // todo add if statement if dense

        if ((m_sychMethod != METHODS_SYNCH::iGSL)
            and (m_sychMethod != METHODS_SYNCH::iDER06)
            and (m_sychMethod != METHODS_SYNCH::iBessel)
            and (m_sychMethod != METHODS_SYNCH::iCSYN)){
            (*p_log)(LOG_ERR,AT) << " only GSL DER06 Bessel or CSYN synchrotron options are "
                                    "avaialble when evolving electrons numerically / or mixed. \n";
            exit(1);
        }

        total_rad.allocate(freq1, freq2, nfreq, nr, true); // always true

        /// allocate memory for synchrotron kernel
        syn_kernel.allocate(ele.numbins, syn.numbins);

        if (m_sychMethod == METHODS_SYNCH::iGSL)
            syn_kernel.setFunc(SynKernel::synchGSL);
        else if (m_sychMethod == METHODS_SYNCH::iDER06)
            syn_kernel.setFunc(SynKernel::synchDermer);
        else if (m_sychMethod == METHODS_SYNCH::iBessel)
            syn_kernel.setFunc(SynKernel::synchBessel);
        else if (m_sychMethod == METHODS_SYNCH::iCSYN)
            syn_kernel.setFunc(SynKernel::cycSynch);
        else {
            (*p_log)(LOG_ERR, AT) << " method for synchrotron is not recognized \n";
            exit(1);
        }


        /// allocate memory for SSC structure and kernel
        if (m_methods_ssc == METHOD_SSC::iNumSSC) {
            ssc.allocate(freq1, freq2, nfreq, nr, is_dense_output); // todo add if statement if dense
            ssc_kernel.allocate(ele.numbins, syn.numbins, ssc.numbins,
                                SSCKernel::sscNava);
            ssc_kernel.evalSSCkernel(ele.e, syn.e, ssc.e);
        }

        /// ---
        model.setSolver();
    }

    void initializeElectronDistribution(double tcomov, double m2, double rho2, double u_p){
        source.N_ele_tot = m2 / CGS::mp;
        source.vol = m2 / rho2;
        source.u_e = u_p * eps_e;
        powerLawElectronDistributionAnalytic(tcomov, source.N_ele_tot, ele.f);
        is_distribution_initialized = true;
        /// apply Deep Newtonian limit
        for (size_t i = 0; i < ele.numbins-1; i++)
            ele.f[i] *= accel_frac;
    }

    void evaluateElectronDistributionNumericMixed(
            double tc, double tc_p1,
            double m, double m_p1,
            double rho_prime, double rho_prime_p1,
            double r, double rp1,
            double dr_comov_, double drp1_comov,
            double u_p, double u_p_p1){

        /// check parameters
        if (ele.e.size() < 1){
            (*p_log)(LOG_ERR,AT) << " electon grid is not initialized\n";
            exit(1);
        }
        if ((B <= 0) || (gamma_min < 1.) || (gamma_max <= 1.) || (gamma_min >= gamma_max) || (accel_frac > 1.)){
            (*p_log)(LOG_ERR,AT)
                << "wrong value, cannot evolve electrons: "
                << "B="<<B<<" gm="<<gamma_min<<" gM="<<gamma_max<<" accel_frac="<<accel_frac<<"\n";
            exit(1);
        }
        dr_comov = dr_comov_;
        source.N_ele_tot = m_p1 / CGS::mp;
        source.u_e = u_p_p1 * eps_e;

//        accel_frac = std::max(accel_frac,1.);
        /// compute analytical electron spectrum
//        double n_inj = (m_p1 - m) / CGS::mp;
//        n_inj *= accel_frac;
//        if (accel_frac < 1.)
//            int z = 1;
        double n_ele_an=0.;
        Vector tmp (ele.numbins,0.);
        powerLawElectronDistributionAnalytic(tc, source.N_ele_tot, tmp);

        /// check continouity
        for (size_t i = 0; i < ele.numbins-1; i++)
            n_ele_an += tmp[i] * (ele.e[i+1]-ele.e[i]);
//        n_ele_out = n_ele_an;
//        ratio_an = source.N_ele_tot / n_ele_an;

        double dt = tc_p1 - tc;
        source.dt = dt;
        /// for adiabatic cooling of electron distribution (Note it may be turned off)
        double volume = m / rho_prime;//r * r * dr_comov;
//        double volume = r * r * dr_comov;
        double volume_p1 = m_p1 / rho_prime_p1;//rp1 * rp1 * drp1_comov;
//        double volume_p1 = rp1 * rp1 * drp1_comov;
        double dlnVdt = num_ele_use_adi_loss ? 1. / dt * (1. - volume / volume_p1) : 0.;
        if ((!std::isfinite(dlnVdt)) || (dlnVdt < 0)){
            (*p_log)(LOG_ERR,AT) << " dlnVdt= "<<dlnVdt<<"\n";
            exit(1);
        }


        /// compute substepping to account for electron cooling timescales
        double dlog_gam = (ele.e[1]-ele.e[0])/ele.e[0];
        double delta_t_syn = CGS::sigmaT * gamma_max * B * B / (6. * M_PI * CGS::me * CGS::c);
        double delta_t_adi = (gamma_max*gamma_max-1.)/(3.*gamma_max*gamma_max)*dlnVdt;
        double delta_t = dlog_gam / (delta_t_syn + delta_t_adi);
//        size_t n_substeps = 0;

        /// if cooling is too slow, we still need to evolve distribution
        size_t min_substeps = 1;
//        max_substeps = 1;
        if (delta_t >= dt/(double)min_substeps) {
            delta_t = dt/(double)min_substeps;
            n_substeps = min_substeps;
        }
        else if (delta_t < dt/(double)max_substeps){
            delta_t = dt/(double)max_substeps;
            n_substeps = max_substeps;
        }
        else
            n_substeps = (size_t)(dt/delta_t);
        if (n_substeps < min_substeps){
            (*p_log)(LOG_ERR,AT) << " n_substeps="<<n_substeps<<"\n";
            exit(1);
        }

        /// update source properties
        source.B = B;
        source.dlnVdt = dlnVdt;
        source.dr = dr_comov_;
        source.vol = volume;
//        double n_ele_num = 0.;
//        ChangCooper model = ChangCooper(source, ele, syn, ssc, syn_kernel, ssc_kernel);

        /// compute radiation numerically from ANALYTIC electron distribution or evolve electron distribution
        if (m_eleMethod==METHODS_SHOCK_ELE::iShockEleMix ) // || gamma_min == 1.
            for (size_t i = 0; i < ele.numbins; i++)
                ele.f[i] = tmp[i];
        /// evolve electron distribution via Chang Cooper Scheme
        else{
            source.N = (m_p1 - m) / CGS::mp / dt; //n_inj / dt; // number of injected electrons per timestep

            /// Assume source/escape do not change during substepping
            model.setSourceFunction(gamma_min, gamma_max, -p);
            model.setEscapeFunction(gamma_min, gamma_max);

            /// add Pair production as a source term
            if (m_methods_ssc == METHOD_SSC::iNumSSC)
                model.addPPSource(total_rad.f);

            for (size_t i = 0; i < (size_t) n_substeps; i++) {

                /// Set gamma_dot terms (uses syn.n[] for photon density for SSC from previous solution)
                model.setHeatCoolTerms(2, total_rad.f);

                /// Assuming sources change during substepping
//            model.setSourceFunction(m_data[Q::igm][it], m_data[Q::igM][it], index, N);
//            model.setEscapeFunction(m_data[Q::igm][it], m_data[Q::igM][it]);

                /// setup grid of differences for escape term
                model.computeDeltaJ();

                /// setup tridiagonal solver arrays
                model.setupVectors(delta_t);

                /// solve system
                model.solve(delta_t);

//                model.resetSolver();
            }
            /// check continouity
//            for (size_t i = 0; i < ele.numbins-1; i++)
//                n_ele_num += ele.f[i] * ele.de[i];
//            n_ele_out = n_ele_num;
        }

        /// clear array (resets everything!)
        model.resetSolver();
    }

    void computeRadiationNumericMixed(size_t it){
        /// apply Deep Newtonian limit TODO !!!
//        for (size_t i = 0; i < ele.numbins-1; i++)
//            ele.f[i] *= accel_frac;
//        n_ele_num*=accel_frac;

        /// 1. Update synchrotron kernel for current B
        syn_kernel.evalSynKernel(ele.e, syn.e, B);

        /// 2. Compute synchrotron emission spectrum
        model.computeSynSpectrum(); // compute syn.j[] syn.f[]

        /// 3. compute SSA
        if (m_methods_ssa == METHODS_SSA::iSSAon)
            model.computeSSA(); // compute syn.a[]

        /// 4. Compute SSC spectrum (using PREVIOUS step photon density 'total_rad.f')
        if (m_methods_ssc == METHOD_SSC::iNumSSC)
            model.computeSSCSpectrum(total_rad.f); // using syn.f[]+ssc.f[] -> compute -> ssc.j[], ssc.f[]

        /// 5. compute new photon density
        double T = dr_comov / CGS::c; // escape time (See Huang+2022)
        double fac = T / source.vol; // volume = 4. * M_PI * r * r * dr_comov;
        for (size_t i = 0; i < syn.numbins; i++)
            total_rad.f[i] = syn.j[i] * fac / (CGS::h * syn.e[i]);
        if ((m_eleMethod==METHODS_SHOCK_ELE::iShockEleNum) && (m_methods_ssc==METHOD_SSC::iNumSSC))
            for (size_t i = 0; i < ssc.numbins; i++)
                total_rad.f[i] += ssc.j[i] * fac / (CGS::h * ssc.e[i]);

        /// 5. Substract photons due to pair-production from new photon density
        for (size_t i = 0; i < ssc.numbins; i++) {
            double pp_loss = model.computePP(ssc.e[i], total_rad.f);
            total_rad.f[i] = std::max(0., total_rad.f[i] - pp_loss);
//            double val1 = ssc.f[i];
//            double val2 = model.computePP(syn.e[i]);
//            std::cout << val1 << " | " << val2 << " -> " << std::max(0., ssc.f[i] - model.computePP(ssc.e[i])) << "\n";
        }

        /// 6. convert photon field back to emissivity (to account for photons lost to pp-production)
        for (size_t i=0; i < syn.numbins; i++) {
            total_rad.j[i] = total_rad.f[i] / fac * (CGS::h * ssc.e[i]);
            total_rad.a[i] = syn.a[i];
        }
        if ((m_eleMethod==METHODS_SHOCK_ELE::iShockEleNum) && (m_methods_ssc==METHOD_SSC::iNumSSC))
            for (size_t i = 0; i < syn.numbins; i++)
                total_rad.a[i] += ssc.a[i];

        /// 7. compute total numer of electrons (solution)
        double n_ele=0.;
        for (size_t i = 0; i < ele.numbins-1; i++)
            n_ele+=ele.f[i]*(ele.e[i+1]-ele.e[i]); //total number of electrons
//        n_ele *= accel_frac;// Deep Newtonian corretion
        if (n_ele <= 0){
            (*p_log)(LOG_ERR,AT) << " n_ele = "<<n_ele<<"\n";
            exit(1);
        }

        /// 8. compute final required emissivity and absorption (normalized)
        for (size_t i = 0; i < syn.numbins; i++) {
            total_rad.j[i] = total_rad.j[i] / n_ele * accel_frac;
            total_rad.a[i] = total_rad.a[i] / n_ele * accel_frac * n_prime; // absorption (depth) requires the comoving particle density
            total_rad.intensity[i] = computeIntensity(
                    total_rad.j[i], total_rad.a[i]*dr_comov, METHOD_TAU::iAPPROX);
        }

        /// store the result in long arrays for entire evolution
        total_rad.save_to_all(it); // store j and a

        /// process synchrotron, ssa, ssc separately for comoving spectra analysis
        if (is_dense_output)
            saveDenseOutput(it, n_ele);

        /// log the result
        if (m_loglevel>=LOG_INFO)
            logResultToConcole(it, n_ele);

    }

    void saveDenseOutput(size_t it, double n_ele){
        // Normalize to get emissivity per particle (for EATS)
        for (size_t i = 0; i < syn.numbins; i++) {
            syn.j[i] = syn.j[i] / n_ele * accel_frac;
            syn.a[i] = syn.a[i] / n_ele * accel_frac * n_prime; // absorption (depth) requires the comoving particle density
            double dtau = syn.a[i]*dr_comov;
            syn.intensity[i] = computeIntensity(syn.j[i], dtau, METHOD_TAU::iAPPROX);
        }

        if (m_methods_ssc != METHOD_SSC::inoSSC) {
            for (size_t i = 0; i < ssc.numbins; i++) {
                ssc.j[i] = ssc.j[i] / n_ele * accel_frac;
                ssc.a[i] = ssc.a[i] / n_ele * accel_frac * n_prime; // absorption (depth) requires the comoving particle density
                double dtau = ssc.a[i]*dr_comov;
                ssc.intensity[i] = computeIntensity(ssc.j[i], dtau, METHOD_TAU::iAPPROX);
            }
        }
        /// store comoving spectra in long arrays for the entire evolution
        ele.save_to_all(it);
        syn.save_to_all(it);
//        total_rad.add_to_all(syn, it);
        if (m_methods_ssc != METHOD_SSC::inoSSC) {
            ssc.save_to_all(it);
//            total_rad.add_to_all(ssc, it);
        }
        double val = std::accumulate(syn.j_all.begin(),syn.j_all.end(),0.);
        if (val == 0. and it > 0){
            (*p_log)(LOG_ERR,AT)<<" sum(syn.j)=0\n";
            exit(1);
        }

    }

    void logResultToConcole(size_t it, double n_ele){

        double nu_tau = -1;
        double _min = std::numeric_limits<double>::max();
        for (size_t i = 0; i < syn.numbins; i++) {
            double dtau = syn.a[i]*dr_comov;
            if (dtau-1. < _min and dtau-1.>0) {
                _min = dtau - 1.;
                nu_tau = syn.e[i];
            }
        }

        /// compute energy conservation (https://www.mpi-hd.mpg.de/personalhomes/frieger/HEA6.pdf)
        double u_ele=0.;
        for (size_t i = 0; i < ele.numbins-1; i++)
            u_ele+=ele.f[i]*ele.e[i]*(ele.e[i+1]-ele.e[i])*CGS::me*CGS::c*CGS::c/source.vol; //total energy density in electrons
        u_ele*=accel_frac;


        /// find numerical gamma_c and compare with analytic one
        if(m_eleMethod==METHODS_SHOCK_ELE::iShockEleNum) {
            double gamma_c_num = 0.;
            for (size_t i = 1; i < ele.numbins; i++) {
                if ((model.gamma_dot_syn(ele.e[i]) > model.gamma_dot_adi(ele.e[i])) && \
                (model.gamma_dot_syn(ele.e[i - 1]) < model.gamma_dot_adi(ele.e[i - 1]))) {
                    gamma_c_num = ele.e[i];
                    break;
                }
            }
            gamma_c = gamma_c_num;
        }

        /// locate spectral limits (synchrotron)
        size_t max_idx;
        auto [min_iter, max_iter] = std::minmax_element(std::begin(syn.j), std::end(syn.j));
        max_idx = std::distance(std::begin(syn.j), max_iter);
        double synch_peak = syn.e[max_idx];
        double synch_peak_flux = syn.e[max_idx] * syn.j[max_idx];

        /// locate spectral limits (ssc)
        size_t max_idx_c = 0;
        double ssc_peak = 0.;
        double ssc_peak_flux = 0.;
        if ((m_eleMethod == METHODS_SHOCK_ELE::iShockEleNum) && (m_methods_ssc == METHOD_SSC::iNumSSC)) {
            auto [min_iter_c, max_iter_c] =
                    std::minmax_element(std::begin(ssc.j), std::end(ssc.j));
            max_idx_c = std::distance(std::begin(ssc.j), max_iter_c);
            ssc_peak = ssc.e[max_idx_c];
            ssc_peak_flux = ssc.e[max_idx_c] * ssc.j[max_idx_c];
        }

        /// energy density of the synchrotron radiation
        double u_syn = 0;
        for (size_t i = 0; i < syn.numbins-1; i++)
//                u_syn += syn.f[i] * (CGS::h * syn.e[i]);
            u_syn += syn.j[i] * (syn.e[i+1]-syn.e[i]);
        u_syn *= (4. * CGS::pi / CGS::c);
//            double U_syn = u_syn * source.vol;

        double u_ssc = 0;
        if (m_methods_ssc == METHOD_SSC::iNumSSC)
            for (size_t i = 0; i < syn.numbins-1; i++)
//                    u_ssc += ssc.f[i] * (CGS::h * syn.e[i]);
                u_ssc += ssc.j[i] * (ssc.e[i+1]-ssc.e[i]);
        u_ssc *= (4. * CGS::pi / CGS::c);
//            double U_ssc = u_ssc * source.vol;

//            std::cout << (u_syn + u_ssc) * source.dt / source.u_e << "\n";
//
        double y = u_ssc/u_syn;

        /// log result
        (*p_log)(LOG_INFO, AT)
                << (is_rs ? "RS " : "FS " )
                << (n_substeps == 0 ? "I " : "E ")
                << "it=" << it
                << " n=" << n_substeps
                << " Ne=" << n_ele / (source.N_ele_tot) // should be 1
                << " Ue=" << u_ele / (source.u_e) // should be 1
                << " | "
                << " gm=" << gamma_min
                << " gc=" << gamma_c
                << " gM=" << gamma_max
                << " | "
                << " y="<<y
                << " syn[" << synch_peak << "]="
                << synch_peak_flux//<< std::accumulate(syn.j.begin(), syn.j.end(),0.)
                << " nu_tau="<<nu_tau
                << " ssc[" << ssc_peak << "]="
                << ssc_peak_flux//<<" ssc="<< std::accumulate(ssc.j.begin(), ssc.j.end(),0.)
                << "\n";

//                << "Ne="<<n_ele/source.u_e
//
//                                   << " N/Nan=" << ratio_an
//                                   << " N/Num=" << ratio_num
//                                   << " Nan/Num=" << ratio_an_num
//                                   << " gm=" << gamma_min
//                                   << " gc=" << gamma_c
//                                   << " gM=" << gamma_max
//                                   << " u_syn=" << u_syn
//                                   << " u_ssc=" << u_ssc
////                                   << " y="<<y
//                                   << " syn[" << synch_peak << "]="
//                                   << synch_peak_flux//<< std::accumulate(syn.j.begin(), syn.j.end(),0.)
//                                   << " nu_tau="<<nu_tau
//                                   << " ssc[" << ssc_peak << "]="
//                                   << ssc_peak_flux//<<" ssc="<< std::accumulate(ssc.j.begin(), ssc.j.end(),0.)
//                                   << "\n";
    }

public: // ---------------------- FOR EATS -------------------------------- //

    double computeFluxDensity(
            EjectaID2::STUCT_TYPE eats_method,
            double & em_prime, double & abs_prime, double Gamma, double GammaSh, double acc_fac,
            double mu, double r, double rsh, double dr, double n_prime, double ne, double theta, double ncells) {

        /// compute beaming and doppler factor
        double beta = EQS::BetaFromGamma(Gamma);
        double a = 1.0 - beta * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        /// convert to the laboratory frame
        double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
        double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

        /// compute shock velocity
        double beta_shock;
        switch (method_shock_vel) {
            case isameAsBW:
                beta_shock = EQS::BetaFromGamma(Gamma);
                break;
            case ishockVel:
                beta_shock = EQS::BetaFromGamma(GammaSh);
                break;
        }
//        double dr_tau = dr; //EQS::shock_delta(r,GammaSh);
        /// mass in the shock downstream
        double m2 = ne * CGS::mp;
        /// Using Johannesson paper 2006 Eq.33 (shock thickness in a burster frame)
        dr = m2 / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
//        double dr = r / (12. * Gamma * Gamma); // for optical depth; (see vanEerten+2010)
//        dr_/=(1. - cos(theta) / ncells);
//        std::cout << "dr="<<dr<<" dr_="<<dr_<<" dr/dr_="<<dr/dr_<<"\n";
        /// compute shock optical depth along the line of sight
        // (Signed) Optical depth through the shell.
        // if negative, face is oriented away from observer.
        double dtau;
        if (mu == beta_shock) dtau = 1.0e100; // HUGE VAL, just in case
        else dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // correct for abberation
        /// compute intensity
        double intensity = computeIntensity(em_lab, dtau, method_tau);
        if ((intensity < 0) || (!std::isfinite(intensity))) {
            (*p_log)(LOG_ERR, AT) << "intensity = " << intensity << "\n";
            exit(1);
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor (See vanEerten+2010 paper, Eq.A9)
        double luminocity = 0.;
        switch (m_method_ne) {
            case iusenprime:
                luminocity = intensity * r * r * (dr / ashock);
                break;
            case iuseNe:
                luminocity = intensity
                             * ne
                             / (2. * CGS::pi * (1. - std::cos(theta) / ncells))
                             / Gamma
                             / (1.0 - mu * beta_shock);
//                if (eats_method==EjectaID2::STUCT_TYPE::iadaptive)
//                    luminocity /= (2. * CGS::pi * (1. - std::cos(theta_h)/ncells));
//                else
//                    luminocity /= (2. * CGS::pi * (1. - std::cos(theta)/ncells));
                break;
        }

        if (luminocity < 0 || !std::isfinite(luminocity)) {
            (*p_log)(LOG_ERR, AT) << "flux_dens_rs = " << luminocity << "\n";
            exit(1);
        }

        return luminocity;
    }
#if 0
    double computeFluxDensity_new(
            EjectaID2::STUCT_TYPE eats_method,
            double & em_prime, double & abs_prime, double Gamma, double GammaSh,
            double mu, double r, double rsh, double dr, double n_prime, double ne, double theta, double ncells) {
        /// compute beaming and doppler factor
        double betaSh = (double)EQS::Beta2(Gamma);
        double a = 1.0 - betaSh * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        /// convert to the laboratory frame
        double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
        double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
        /// compute shock velocity
        long double beta_shock = EQS::Beta2(GammaSh);
//        double beta_shock = EQS::Beta2(GammaSh);

        /// Shock thicnkness using Johannesson paper 2006 Eq.33
        dr = (ne * CGS::mp) / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
//        dr = r / (12. * Gamma * Gamma); // Van Earten
        /// compute optical depth along line of sight
        double dtau;
        if(mu == beta_shock) dtau = 1.0e100; // HUGE VAL, just in case
        else dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock);
        /// compute observer intensity
        double intensity = computeIntensity(em_lab, dtau, method_tau);
        /// compute flux density
        double flux_dens=0., ashock=0.;
        switch (m_method_ne) {
            case iusenprime:
//                m2 = ne * CGS::mp;
//                dr = m2 / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
                ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor (See vanEerten+2010 paper, Eq.A9)
//                dr_scaled = dr / ashock; // / Gamma; // moved from thicknesss
                flux_dens = intensity * r * r * dr / ashock;
                break;
            case iuseNe:
                /// combine dr equation with F = I * r * r * dr / ashock
                flux_dens = intensity
                            * ne
                            / Gamma
                            / (1.0 - mu * beta_shock);
                if (eats_method == EjectaID2::STUCT_TYPE::ipiecewise)
                    /// from Johanesson dr equation
                    flux_dens /= (2. * CGS::pi * (1. - std::cos(theta)/ncells));
                break;
        }
        return flux_dens;
    }


    double computeFluxDensity_(
            EjectaID2::STUCT_TYPE eats_method,
            double & em_prime, double & abs_prime, double Gamma, double GammaSh,
            double mu, double r, double rsh, double dr, double n_prime, double ne, double theta, double ncells){
        /// compute beaming and doppler factor
        double betaSh = (double)EQS::Beta2(Gamma);
        double a = 1.0 - betaSh * mu; // beaming factor
        double delta_D = Gamma * a; // doppler factor
        /// convert to the laboratory frame
        double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
        double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)
        /// compute shock velocity
        double beta_shock;
        switch (method_shock_vel) {
            case isameAsBW:
                beta_shock = EQS::Beta(Gamma);
                break;
            case ishockVel:
                beta_shock = EQS::Beta(GammaSh);
                break;
        }
        double dr_tau = dr; //EQS::shock_delta(r,GammaSh);
        /// mass in the shock downstream
        double m2 = ne * CGS::mp;
        /// Using Johannesson paper 2006 Eq.33
        dr = m2 / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
//        double dr_ = r / (12. * Gamma * Gamma);
//        dr_/=(1. - cos(theta) / ncells);
//        std::cout << "dr="<<dr<<" dr_="<<dr_<<" dr/dr_="<<dr/dr_<<"\n";
        /// compute shock optical depth along the line of sight
        // (Signed) Optical depth through the shell.
        // if negative, face is oriented away from observer.
        double dtau;
        if(mu == beta_shock) dtau = 1.0e100; // HUGE VAL, just in case
        else dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock); // correct for abberation
        /// compute intensity
        double intensity = computeIntensity(em_lab, dtau, method_tau);
        if ((intensity < 0) || (!std::isfinite(intensity))) {
            (*p_log)(LOG_ERR, AT) << "intensity = " << intensity << "\n";
            exit(1);
        }
        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor (See vanEerten+2010 paper, Eq.A9)
        double flux_dens=0.;
        switch (m_method_ne) {
            case iusenprime:
                flux_dens = intensity * r * r * dr / ashock;
                break;
            case iuseNe:
                flux_dens = intensity
                          * ne
                          / (2. * CGS::pi * (1. - std::cos(theta)/ncells))
                          / Gamma
                          / (1.0 - mu * beta_shock);
//                if (eats_method==EjectaID2::STUCT_TYPE::iadaptive)
//                    flux_dens /= (2. * CGS::pi * (1. - std::cos(theta_h)/ncells));
//                else
//                    flux_dens /= (2. * CGS::pi * (1. - std::cos(theta)/ncells));
                break;
        }

        if (flux_dens < 0 || !std::isfinite(flux_dens)) {
            (*p_log)(LOG_ERR, AT) << "flux_dens_rs = " << flux_dens << "\n";
            exit(1);
        }

        return flux_dens;

        switch (p_pars->m_method_eats) {
        case EjectaID2::iadaptive:
            switch (p_pars->m_method_rad) {
                // fluxDensAdaptiveWithComov()
                case icomovspec:
                    switch (p_pars->p_mphys->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:
                    switch (p_pars->p_mphys->m_method_ne){
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
                    switch (p_pars->p_mphys->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:

                    switch (p_pars->p_mphys->m_method_ne){
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

//        return flux_dens;
    }


    double computeFluxDensity_old(
            EjectaID2::STUCT_TYPE eats_method,
            double & em_prime, double & abs_prime, double Gamma, double GammaSh,
            double mu, double r, double rsh, double dr, double n_prime, double ne, double theta, double ncells){


        /// Using Johannesson paper 2006 Eq.33
//        double m2 = ne * CGS::mp;
//        std::cout << " dr = "<< dr << " r-rsh =" << rsh - r << "\n";
//        theta = 0.1;
//        dr = std::abs(rsh - r);


        /// compute beaming and doppler factor
        double betaSh = (double)EQS::Beta2(Gamma);
        double a = 1.0 - betaSh * mu; // beaming factor
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
//        beta_shock = EQS::Beta(Gamma);
        double dr_tau = dr; //EQS::shock_delta(r,GammaSh);
//        double tmp_ = dr / dr_tau;
//        dr_tau /= ashock;
        double m2 = ne * CGS::mp;
        dr = m2 / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
        double dtau;
        if(mu == beta_shock)
            dtau = 1.0e100; // HUGE VAL, just in case
        else
            dtau = abs_lab * dr * (1. - mu * beta_shock) / (mu - beta_shock);
        double intensity = computeIntensity(em_lab, dtau, method_tau);

        if ((intensity < 0) || (!std::isfinite(intensity))) {
            (*p_log)(LOG_ERR, AT) << "intensity = " << intensity << "\n";
            exit(1);
        }

//        double flux_dens = intensity * r * r * dr;

//        if (eats_method == EjectaID2::STUCT_TYPE::iadaptive)
//            intensity = std::abs(mu * beta_shock) * intensity / (1.0 - mu * beta_shock);

        double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor (See vanEerten+2010 paper, Eq.A9)
        double flux_dens=0.;
        switch (m_method_ne) {
            case iusenprime:
                flux_dens = intensity * r * r * dr / ashock;
                break;
            case iuseNe:
                // TODO this nprime / ne is so that results afree with when we use nprime not Ne option
//                double m2 = ne * CGS::mp;
//                dr = m2 / (2. * CGS::pi * CGS::mp * r * r * (1. - cos(theta) / ncells) * Gamma * n_prime);
//                flux_dens = intensity
//                          * r * r * dr / ashock
//                          * n_prime;// / ne;
                flux_dens = intensity
                          * ne
                          / (2. * CGS::pi * (1. - std::cos(theta)/ncells))
                          / Gamma
                          / (1.0 - mu * beta_shock);
                break;
        }
//        flux_dens = intensity / ashock;

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
                    switch (p_pars->p_mphys->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:
                    switch (p_pars->p_mphys->m_method_ne){
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
                    switch (p_pars->p_mphys->m_method_ne) {
                        case iusenprime:
                            flux_dens = intensity * r * r * dr; //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
                            break;
                        case iuseNe:
                            flux_dens = intensity * r * r * dr * n_prime;
                            break;
                    }
                    break;
                case iobservflux:

                    switch (p_pars->p_mphys->m_method_ne){
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
#endif


    double fluxDens(
            EjectaID2::STUCT_TYPE eats_method, double theta, double ncells,
            size_t ia, size_t ib, double freq, double Gamma, double GammaSh, double acc_fac,
            double mu, double r, double rsh, double dr, double n_prime, double ne, Vector & r_arr){

        Vector & freq_arr = total_rad.e;
        if (freq >= freq_arr[total_rad.numbins-1]){
            (*p_log)(LOG_ERR,AT)<<" freq_prime="<<freq
                                <<" > grid freq[-1]="<<freq_arr[total_rad.numbins-1]
                                <<" increase 'freq2' parameter\n";
            exit(1);
        }
        if (freq <= freq_arr[0]){
            (*p_log)(LOG_ERR,AT)<<" freq_prime="<<freq
                                <<" < grid freq[0]="<<freq_arr[0]
                                <<" decrease 'freq1' parameter\n";
            exit(1);
        }

        size_t ia_nu = findIndex(freq, freq_arr, freq_arr.size());
        size_t ib_nu = ia_nu + 1;

        /// interpolate the emissivity and absorption coefficines
        Interp2d int_em(freq_arr, r_arr, total_rad.j_all);
        double em_prime = int_em.InterpolateBilinear(freq, r, ia_nu, ib_nu, ia, ib);

        Interp2d int_abs(freq_arr, r_arr, total_rad.a_all);
        double abs_prime = int_abs.InterpolateBilinear(freq, r, ia_nu, ib_nu, ia, ib);

        /// compute flux density
        return computeFluxDensity(eats_method,em_prime, abs_prime, Gamma, GammaSh, acc_fac,
                                  mu, r, rsh, dr, n_prime, ne, theta, ncells);
    }
};


// --------------- OLD code to compute electron dist. at shock when Magnetar wind injection present
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
// #pragma omp parallel for num_threads( 6 )
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
