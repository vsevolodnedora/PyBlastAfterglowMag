//
// Created by vsevolod on 02/04/23.
//

#ifndef SRC_EJECTA_GRB_H
#define SRC_EJECTA_GRB_H

#include "pch.h"
#include "utils.h"
#include "base_equations.h"
#include "interpolators.h"
#include "ode_solvers.h"
#include "quadratures.h"
#include "rootfinders.h"
#include "observer.h"
#include "synchrotron_an.h"

#include "magnetar.h"
//#include "pulsar_wind_nebula.h"
#include "blastwave_base.h"
#include "blastwave_rad.h"
#include "blastwave_dyn.h"

class GRB{
    std::unique_ptr<logger> p_log;
//    bool run_jet_bws = false;
    bool is_jBW_init = false;
    bool is_jet_obsrad_pars_set = false;
    bool is_jet_struct_set = false;
    int n_ode_eq;
    LatStruct::METHOD_eats jet_eats_method{};
    int m_loglevel{};
    Vector & t_arr;
    double jet_layer_fnu_stop_frac;
protected:
    LatStruct jetStruct{};
    std::vector<std::unique_ptr<RadBlastWave>> p_bws_jet;
public:
    double jet_rtol = 1e-5;
    bool run_jet_bws = false;
    bool is_jet_obs_pars_set = false;
    bool is_jet_anal_synch_computed = false;
    GRB(Vector & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "GRB");
    }
    size_t getNeq() const {
        if (!run_jet_bws)
            return 0;
        else
            return (p_bws_jet.size() * p_bws_jet[0]->getNeq());
    }
    void setJetStructAnalytic(StrDbMap pars, StrStrMap opts){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true, true);
        if (!run_jet_bws)
            return;
        std::string _jet_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        jet_eats_method = LatStruct::setEatsMethod(_jet_eats_method);
        jetStruct.initAnalytic( pars, opts, _jet_eats_method, m_loglevel );
        is_jet_struct_set = true;
    }
    std::vector<std::unique_ptr<RadBlastWave>> & getBWs(){return p_bws_jet;}
//    static auto listParsNumericJetStruct(){ return LatStruct::listParsCustomStruct(); }
    void setJetStructNumeric( Vector & dist_thetas, Vector & dist_EEs, Vector & dist_Yes, Vector & dist_s, Vector & dist_Gam0s, Vector & dist_MM0s,
                              bool force_grid, std::string eats_method ){
        (*p_log)(LOG_ERR,AT) << " not finished...\n"; exit(1);
        jetStruct.initCustom( dist_thetas, dist_EEs, dist_Yes, dist_s, dist_Gam0s, dist_MM0s,
                              force_grid, eats_method, m_loglevel);
        is_jet_struct_set = true;
    }
    void setJetBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true, true);
        if (!run_jet_bws)
            return;

        size_t n_layers = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
        for(size_t i = 0; i < n_layers; i++){
//            DynRadBlastWave x(t_grid, 0, i, p_pars->loglevel);
            p_bws_jet.emplace_back( std::make_unique<DynRadBlastWave>(t_arr, 0, i, m_loglevel) );
            p_bws_jet[i]->setAllParametersForOneLayer(jetStruct, pars, opts, i, ii_eq);
            p_bws_jet[i]->setEatsPars(pars, opts);
            p_bws_jet[i]->getSynchAnPtr()->setPars( pars, opts );
            ii_eq += p_bws_jet[i]->getNeq();
        }

        /// parameters for EATS methods
//        p_pars->jet_eats_method = LatStruct::setEatsMethod(
//                getStrOpt("method_eats",opts,AT,p_log,"",true));
        jet_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);
        jet_layer_fnu_stop_frac= getDoublePar("fnu_min_frac",pars,AT,p_log,-1, true);

        (*p_log)(LOG_INFO,AT) << "finished initializing jet...\n";
        is_jet_obs_pars_set = true;
        is_jBW_init = true;
    }
    void infoFastestLayer(size_t it, logger &sstream){
        size_t fastest = 0;
        double mom = 0;
        for (size_t il = 0; il < jetStruct.nlayers; il++){
            if (p_bws_jet[il]->getVal(RadBlastWave::Q::imom,(int)it) > mom) {
                mom = p_bws_jet[il]->getVal(RadBlastWave::Q::imom, (int) it);
                fastest = il;
            }
        }
        sstream << "[Jet: "<<"[l="<<fastest<<"]"
                << " Mom=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::imom, (int) it))
                << " R=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::iR, (int) it))
                << " Eint=" << string_format("%.2e",p_bws_jet[fastest]->getVal(RadBlastWave::Q::iEint2, (int) it))
                << "] ";
    }
};

#endif //SRC_EJECTA_GRB_H
