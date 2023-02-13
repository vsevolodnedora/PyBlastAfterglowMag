//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_BLASTWAVE_SET_H
#define SRC_BLASTWAVE_SET_H

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
#include "blastwave_base.h"
#include "blastwave_rad.h"
#include "blastwave_dyn.h"

/*
 * Class to hold blastwaves with the same angular coordinate (shells of the same layer)
 * Here the velocity distribution, relative position and photosphere can be found.
 */
class CumulativeShell{
    size_t m_nshells;
    std::unique_ptr<logger> p_log;
    std::vector<std::unique_ptr<RadBlastWave>> p_bws_ej;
    std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws_ej;
    std::vector<std::vector<size_t>> relative_position;
    std::vector<size_t> m_idxs;
    Vector m_rho;
    Vector m_radii;
    Vector m_kappas;
public:

    CumulativeShell( Array t_grid, size_t nshells, int ilayer, int loglevel ){
        m_nshells = nshells;

        p_log = std::make_unique<logger>(std::cout, std::cerr, CurrLogLevel, "LatStruct");
        for (size_t ishell = 0; ishell < nshells; ishell++)
            p_bws_ej.emplace_back( std::make_unique<DynRadBlastWave>(t_grid, ishell, ilayer, loglevel ) );

        relative_position.resize( nshells );
        for (auto & arr : relative_position)
            arr.resize( t_grid.size() );

        m_rho.resize(nshells, 0.0);
        m_radii.resize(nshells, 0.0);
    }
    std::unique_ptr<RadBlastWave> & getBW(size_t ish){
        if (p_bws_ej.empty()){
            (*p_log)(LOG_ERR, AT) << " shell does not contain blast waves\n"; exit(1);
        }
        if (ish > p_bws_ej.size()){
            (*p_log)(LOG_ERR, AT) << "invalid memory accessed\n"; exit(1);
        }
        return p_bws_ej[ish];
    }
    inline size_t nBWs() const { return p_bws_ej.size();}
    inline std::vector<std::unique_ptr<RadBlastWave>> & getBWs() {return p_bws_ej; }
    inline std::vector<size_t> & getCurrentIndexes( ) { return m_idxs; }

#if 0
    void evaluateOptDepthAtCurrR(size_t i_cur, double x, const double * Y){
        /// evaluate relative positions of shells at this time
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws_ej[i]->getLastVal(std::static_pointer_cast<BlastWaveBase::Q>())
        }
        /// get indexes of sorted shells by radius
        sort_indexes(m_radii, m_idxs);
        /// evaluate opacity of each shell up to current
        double tau = 0.;
        for (auto & i : m_idxs){
            double kappa = 10.;
            double dr = 0.;
            if (i == 0) {
                double dr1 = 0.5 * (m_radii[0] - m_radii[1]);
                double dr2 = 0.5 * (m_radii[1] - m_radii[2]);
                dr = linearExtrapolate(m_radii[1], m_radii[2], dr1, dr2, m_radii[0]);
            }
            else if (i==m_idxs.back()){
                size_t
                double drnm2 = 0.5 * (m_radii[0] - m_radii[1]);
                double drnm1 = 0.5 * (m_radii[1] - m_radii[2]);
                dr = linearExtrapolate(m_radii[1], m_radii[2], dr1, dr2, m_radii[0]);
            }
            else
                dr = 0.5 * (m_radii[i+1]-m_radii[i-1])
            tau +=
        }


        /// optical depth at a given radius
        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws_ej[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

    }





    inline void evalRelativePosition(double x, const double * Y){
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            double r = Y[bw->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
            m_radii[i] = r == 0. ? 1e90 : r;
//            radii[i] = p_bws_ej[i]->getLastVal(std::static_pointer_cast<BlastWaveBase::Q>())
        }
        sort_indexes(m_radii, m_idxs);

        iota(current_indexes.begin(), current_indexes.end(), 0);
        stable_sort(current_indexes.begin(), current_indexes.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
//        std::cout << radii << "\n";
        current_indexes = sort_indexes(m_radii);
//        std::cout << idxes << "\n";
//        if (idxes[0]!=0){
//            std::cout << x << " " << x / CGS::day << " " << x / CGS::year;
//            std::cout << radii << "\n" ;
//            std::cout << idxes << "\n" ;
//            exit(1);
//        }
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            bw->getPars()->i_position_in_layer = current_indexes[i];
        }

        std::vector<std::unique_ptr<RadBlastWave>> p_sorted_bws;
        for (size_t & idx : current_indexes){
            p_sorted_bws.emplace_back(p_bws_ej[idx]);
        }

    }
    ///
    void evaluateDensityProfileFromBWpositions(){
        for (size_t & idx : current_indexes) {
            auto & bw = p_bws_ej[idx];

            double rho_i = bw->getPars()->M0 / (4./3.*CGS::pi*())
        }

        Vector rho ( m_nshells, 0.0 );
        for (size_t i = 0; i < m_nshells; i++){
            size_t idx = current_indexes[i];

        }
    }

    void evaluateOpticalDepthToMagnetarEmission(size_t i_cur, double x, const double * Y){
        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        double r_cur = Y[bw_cur->getPars()->ii_eq + DynRadBlastWave::Q_SOL::iR];
        for (size_t i=0; i<bw_cur->getPars()->i_position_in_layer; i++) {
            auto &bw = p_bws_ej[i];
            tau += bw->getVal(RadBlastWave::Q::ikappa, bw->getPars()->)
        }

        double tau = 0;
        auto & bw_cur = p_bws_ej[i_cur];
        for (size_t i=0; i<m_nshells; i++){
            auto & bw = p_bws_ej[i];
            if (r_cur > bw->)
        }

        for (size_t i = 0; i < i_position; i++){
            tau +=
        }
    }

#endif
};

class GRB{
    std::unique_ptr<logger> p_log;
//    bool run_jet_bws = false;
    bool is_jBW_init = false;
    bool is_jet_obsrad_pars_set = false;
    bool is_jet_struct_set = false;
    int n_ode_eq;
    LatStruct::METHOD_eats jet_eats_method{};
    int m_loglevel;
    Array & t_arr;
    double jet_layer_fnu_stop_frac;
protected:
    LatStruct jetStruct{};
    std::vector<std::unique_ptr<RadBlastWave>> p_bws_jet;
public:
    double jet_rtol = 1e-5;
    bool run_jet_bws = false;
    bool is_jet_obs_pars_set = false;
    bool is_jet_anal_synch_computed = false;
    GRB(Array & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "GRB");
    }
    size_t getNeq() const {
        if (!run_jet_bws)
            return 0;
        else
            return (p_bws_jet.size() * p_bws_jet[0]->getNeq());
    }
    void setJetStructAnalytic(StrDbMap pars, StrStrMap opts){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true);
        if (!run_jet_bws)
            return;
        std::string _jet_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        jet_eats_method = LatStruct::setEatsMethod(_jet_eats_method);
        jetStruct.initAnalytic( pars, opts, _jet_eats_method, m_loglevel );
        is_jet_struct_set = true;
    }
    std::vector<std::unique_ptr<RadBlastWave>> & getBWs(){return p_bws_jet;}
//    static auto listParsNumericJetStruct(){ return LatStruct::listParsCustomStruct(); }
    void setJetStructNumeric( Vector & dist_thetas, Vector & dist_EEs, Vector & dist_Gam0s, Vector & dist_MM0s,
                              bool force_grid, std::string eats_method ){
        (*p_log)(LOG_ERR,AT) << " not finished...\n"; exit(1);
        jetStruct.initCustom( dist_thetas, dist_EEs, dist_Gam0s, dist_MM0s, force_grid, eats_method,
                              m_loglevel);
        is_jet_struct_set = true;
    }
    void setJetBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq){
        run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true);
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
};

class RadGRB : protected GRB{
    std::unique_ptr<logger> p_log;
    Image tmp_im_pj;
    Image tmp_im_cj;
    Image tmp_im;
    RadGRB(Array & t_arr, int loglevel) : GRB(t_arr, loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "RadGRB");
        auto ncells = p_bws_jet[0]->getPars()->ncells;
        (*p_log)(LOG_INFO, AT) << " allocating memory for images with ncells="<<ncells<<"\n";
        tmp_im_pj.resize((int)ncells); tmp_im_cj.resize((int)ncells); tmp_im.resize((int)ncells);
    }


private:

public:
    void getJetLightCurve( VecVector & lc, Vector & obs_times, Vector & obs_freqs ){
        (*p_log)(LOG_INFO,AT)<<" starting grb light curve calculation\n";
        size_t nlayers = p_bws_jet.size();
        auto ncells = p_bws_jet[0]->getPars()->ncells;
        VecVector light_curves(nlayers); // [i_layer][i_time]
        light_curves.resize(nlayers);
        for (auto & arr : light_curves)
            arr.resize(obs_times.size(), 0.);
        double flux;
        double rtol = jet_rtol;
        for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
            auto & model = p_bws_jet[ilayer];
            (*p_log)(LOG_INFO,AT)<<" GRB LC ntimes="<<obs_times.size()<<" theta_layer="<<ilayer<<"/"<<nlayers<<
                                 " phi_cells="<<LatStruct::CellsInLayer(model->getPars()->ilayer)<<"\n";
            model->evalForwardShockLightCurve(tmp_im, tmp_im_pj, tmp_im_cj, light_curves[ilayer], obs_times, obs_freqs);
        }
        (*p_log)(LOG_INFO,AT)<<" grb light curve is computed\n";
    }

};

class Ejecta{
    VelocityAngularStruct ejectaStructs{};
    std::vector<std::unique_ptr<CumulativeShell>> p_ej;
    std::unique_ptr<logger> p_log;
    bool is_ejBW_init = false;
    bool is_ejecta_obsrad_pars_set = false;
    bool is_ejecta_struct_set = false;
    double jet_layer_fnu_stop_frac=1e-5;
    int n_ode_eq;
    int m_loglevel;
    LatStruct::METHOD_eats ejecta_eats_method{};
    Array & t_arr;
public:
    bool run_ej_bws = false;
    bool is_ejecta_obs_pars_set = false;
    bool is_ejecta_anal_synch_computed = false;
    Ejecta(Array & t_arr, int loglevel) : t_arr(t_arr), m_loglevel(loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Ejecta");
    }
    size_t getNeq() const {
        if (!run_ej_bws)
            return 0;
        else
            return (p_ej.size() * p_ej[0]->nBWs() * p_ej[0]->getBW(0)->getNeq());
    }
    size_t nlayers() const { return p_ej.size(); }
    size_t nshells() const { return p_ej[0]->getBWs().size(); }
    int ncells() const { return (int)p_ej[0]->getBW(0)->getPars()->ncells; }
    std::vector<std::unique_ptr<CumulativeShell>> & getShells(){ return p_ej;}
    /// set ejecta lateral & velocity structure
//    static auto listParsAnalyticEjectaStruct(){ std::cerr << AT << " not implemented\n"; exit(1); }
    void setEjectaStructAnalytic(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_ERR,AT) << " not implimeneted\n Exiting...";
        exit(1);
    }
//    static std::vector<std::string> listParsNumericEjectaStruct(){ return VelocityAngularStruct::list_pars_v_ns(); }
    void setEjectaStructNumericUniformInTheta(Vector & dist_thetas0, Vector & dist_betas, Vector & dist_ek, size_t nlayers, double mfac, StrStrMap & opts){
        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log,true);
        if (!run_ej_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initUniform(dist_thetas0,dist_betas,dist_ek, nlayers,mfac,
                                  ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaStructNumeric(Vector dist_thetas, Vector dist_betas, VecVector dist_ek, bool force_grid, StrStrMap & opts){
        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log,true);
        if (!run_ej_bws)
            return;
        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
        ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
        ejectaStructs.initCustom(dist_thetas, dist_betas, dist_ek, force_grid,
                                 ej_eats_method, m_loglevel);
        is_ejecta_struct_set = true;
    }
    void setEjectaBwPars(StrDbMap pars, StrStrMap opts, size_t ii_eq, size_t n_layers_jet){

        run_ej_bws = getBoolOpt("run_ej_bws", opts, AT,p_log,true);
        if (!run_ej_bws)
            return;
        bool is_within = false;
        std::vector<size_t> which_within{};
        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t nshells = ejectaStructs.nshells;
        size_t n_layers_ej = ejectaStructs.structs[0].nlayers;
        std::vector<std::vector<size_t>> n_empty_images_layer_shell;
        for (auto & n_empty_images_layer : n_empty_images_layer_shell)
            n_empty_images_layer.resize(nshells);
        for(size_t il = 0; il < n_layers_ej; il++){
            p_ej.push_back( std::make_unique<CumulativeShell>(t_arr, nshells, il, p_log->getLogLevel()) );
            for (size_t ish = 0; ish < nshells; ish++){
                auto & bw = p_ej[il]->getBW(ish);
                auto & struc = ejectaStructs.structs[ish];
                bw->setAllParametersForOneLayer(struc, pars, opts, il, ii_eq);
/// Override the layer-to-use
                if (bw->getPars()->which_jet_layer_to_use == 0){
                    bw->getPars()->which_jet_layer_to_use = 0; // the fastest
                }
                else if(n_layers_ej==0){
                    n_ejecta_empty_images += 1;
                    n_empty_images_layer_shell[ish].emplace_back(il);
//                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
                }
                else if(n_layers_ej == 0){
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
                ii_eq += bw->getNeq();

                bw->setEatsPars(pars, opts);
                bw->getSynchAnPtr()->setPars( pars, opts );

//                ii++;
            }
        }
# if 0
        size_t ii = 0;
        bool is_within = false;
        std::vector<size_t> which_within{};

        size_t n_ejecta_empty_images = 0;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        size_t n_layers_i = ejectaStructs.structs[0].nlayers;

        for(size_t i = 0; i < ejectaStructs.nshells; i++){
            std::vector<size_t> n_empty_images_layer;
            auto & struc = ejectaStructs.structs[i];
            if (n_layers_i != struc.nlayers){
                (*p_log)(LOG_ERR,AT)<<" expected nlayers="<<n_layers_i<<" (from 0th shell) got="<<struc.nlayers<<"\n";
                exit(1);
            }
//            size_t n_layers_i = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for(size_t j = 0; j < n_layers_i; j++) { // ii = il + nlayers * ish
                p_bws_ej.emplace_back(std::make_unique<DynRadBlastWave>(t_grid, i, j, p_pars->loglevel));
                setAllParametersForOneLayer(struc, *(p_bws_ej[ii]), pars, opts, j, ii_eq);
                /// Override the layer-to-use
                if (p_bws_ej[ii]->getPars()->which_jet_layer_to_use == 0){
                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = 0; // the fastest
                }
                else if(n_layers_i==0){
                    n_ejecta_empty_images += 1;
                    n_empty_images_layer.emplace_back(j);
//                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
                }
                else if(n_layers_i == 0){
                    // NO jet structure was set, so exiting I guess... :)
                    // TODO THIS MIGHT BE WRONG -- why 'n_layers_i'
                }
                else if(n_layers_jet == 0){
                    // NO jet structure was set, so exiting I guess... :)
                }
                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
                }
                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use < n_layers_jet) &&
                         (p_bws_ej[ii]->getPars()->which_jet_layer_to_use > -1)){
                    //
                }
                else{
                    std::cerr << " which_jet_layer_to_use="<<p_bws_ej[ii]->getPars()->which_jet_layer_to_use
                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
                              <<" for the slowest"
                              <<" or any N in between the two for a specific jet layer \n"
                              << "Exiting..."
                              << "\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }
                ii_eq += p_bws_ej[ii]->getNeq();

                p_bws_ej[ii]->setEatsPars(pars, opts);
                p_bws_ej[ii]->getSynchAnPtr()->setPars( pars, opts );

                ii++;
            }

            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(i);
                n_empty_images.emplace_back(n_empty_images_layer);
            }
        }
#endif
        is_ejBW_init = true;
        is_ejecta_obs_pars_set = true;

        if ((p_log->getLogLevel() > LOG_WARN)) {
            if (n_ejecta_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta blastwave is NOT initialized for total n="
                      << n_ejecta_empty_images << " layers. Specifically:\n";
                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
                    auto &ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
                    size_t n_layers_i = ejectaStruct.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [";
                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
                        ccerr << n_empty_images[ish][il] << " ";
                    }
                    ccerr << "] / (" << n_layers_i << " total layers) \n";
                }
            }
        }

        ej_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);

        (*p_log)(LOG_INFO,AT) << "finished initializing ejecta...\n";
    }
    double ej_rtol = 1e-5;
};

class EvolveODEsystem{
    struct Pars{
        Pars(
                std::unique_ptr<Magnetar> & p_magnetar,
                std::unique_ptr<GRB> & p_grb,
                std::unique_ptr<Ejecta> & p_ej,
//                std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
//                std::vector<std::unique_ptr<CumulativeShell>> & p_ej,
//                bool run_magnetar, bool run_jet_bws, bool run_ej_bws,
                Array & t_grid
        ) :
        p_magnetar(p_magnetar), p_grb(p_grb), p_ej(p_ej),
//            run_magnetar(run_magnetar), run_jet_bws(run_jet_bws), run_ej_bws(run_ej_bws),
            t_grid(t_grid){
//            if (!p_bws_jet.empty()) {
//                n_j_bws = p_bws_jet.size();
//                n_eq_j_bws = p_bws_jet[0]->getNeq();
//            }
//            else {
//                n_j_bws = 0;
//                n_eq_j_bws = 0;
//            }
//            if(!p_ej.empty()) {
//                n_ej_bws = p_ej.size()*p_ej[0]->nBWs();//p_bws_ej.size();
//                n_eq_ej_bws= p_ej[0]->getBW(0)->getNeq();//p_bws_ej[0]->getNeq();
//            }
//            else {
//                n_ej_bws = 0;
//                n_eq_ej_bws = 0;
//            }
//            if (run_magnetar){
//                n_magnetar_eqs = p_magnetar->getNeq();
//            }
//            else{
//                n_magnetar_eqs = 0;
//            }
////            n_eq_j_bws = p_bws_jet[0]->getNeq();
////            n_eq_ej_bws= p_bws_ej[0]->getNeq();
//            n_tot_eqs  = n_magnetar_eqs
//                       + n_j_bws // number of 'jet' blast waves
//                       * n_eq_j_bws // number of equations in each
//                       + n_ej_bws // number of 'ejecta' blast waves
//                       * n_eq_ej_bws; // number of equations in each;
//            if (n_tot_eqs < 1){
//                std::cerr << AT << "\n No equations to evolve; Both, jet and ejecta seems to be not given. Exiting...\n";
//                exit(1);
//            }
        }
        Array & t_grid;
        std::unique_ptr<Magnetar> & p_magnetar;
//        std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet;
//        std::vector<std::unique_ptr<CumulativeShell>> & p_ej;
        std::unique_ptr<GRB> & p_grb;
        std::unique_ptr<Ejecta> & p_ej;
//        bool run_magnetar{};
//        bool run_jet_bws{};
//        bool run_ej_bws{};
//        size_t n_magnetar_eqs = 0;
//        size_t n_j_bws = 0;
//        size_t n_ej_bws = 0;
//        size_t n_eq_j_bws = 0;
//        size_t n_eq_ej_bws = 0;
//        size_t n_tot_eqs = 0;
//        // ---
//        size_t n_layers_j = 0;
//        size_t n_layers_ej = 0;
//        size_t n_shells_j = 0;
//        size_t n_shells_ej = 0;
        // ---
//        bool use_cthetajet_lim_for_rho = false;
//        bool use_rjet_lim_for_rho = false;
        // ---
//        bool is_entered = false;
//        double r_ej_ent = -1;
        // ----

//        inline size_t idx_j(const size_t i){
////            auto * p_pars = (struct Pars *) pars;
//            auto & bwPtr = p_bws_jet[i];
//            return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_j * bwPtr->getPars()->ishell);
//        }
//        inline size_t idx_ej(const size_t i){
////            auto * p_pars = (struct Pars *) pars;
//            auto & bwPtr = p_bws_ej[i];
//            size_t _idx_j = n_eq_j_bws * n_j_bws;
//            return _idx_j + bwPtr->getNeq() * (bwPtr->getPars()->ilayer + n_layers_ej * bwPtr->getPars()->ishell);
//        }
        // ---
        size_t ix = 0;  // index latest solution
        double dx = 0;
        double x = 0;
//        void collision(double const * Y, size_t ij, size_t iej){
//            double gJ=Y[idx_j(ij)+DynRadBlastWave::Q_SOL::iGamma];
//            mJ=M2*M0 + M0; eJ=Eint2*M0; gAdiJ=gammaAdi
//            gEJ=Gamma_ej; mEJ=M2_ej*M0_ej+M0_ej; eEJ=Eint2_ej*M0_ej; gAdiEJ=gammaAdi_ej
//            i_gM=Gamma_ej; i_mM=M2_ej*M0_ej+M2*M0; i_eM=Eint2_ej*M0_ej
//        }
        size_t i_restarts = 0;
        int n_tot_eqs = 0;
//        int n_j_bws = 0;
//        int n_ej_bws = 0;
    };
    Pars * p_pars;
    std::unique_ptr<logger> p_log;
public:
    EvolveODEsystem(std::unique_ptr<Magnetar> & p_mag,
                    std::unique_ptr<GRB> & p_grb,//std::vector<std::unique_ptr<RadBlastWave>> & p_bws_jet,
                    std::unique_ptr<Ejecta> & p_ej,//std::vector<std::unique_ptr<CumulativeShell>> & p_ej,
//                    bool run_magnetar,
                    Array & t_grid,
//                  size_t n_shells_j, size_t n_shells_ej, size_t n_layers_j, size_t n_layers_ej,
                    const Integrators::METHODS integrator,
                    int loglevel
    ){
        p_pars = new Pars(p_mag, p_grb, p_ej, t_grid);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EvolveODEsystem");
//        p_pars->n_layers_j  = run_jet_bws ? n_layers_j : 0;
//        p_pars->n_layers_ej = run_ej_bws ? n_layers_ej : 0;
//        p_pars->n_shells_j  = run_jet_bws ? n_shells_j : 0;
//        p_pars->n_shells_ej = run_ej_bws ? n_shells_ej : 0;
        // ----

        // allocate memory for the IC and solution
//        p_pars->n_j_bws =
        p_pars->n_tot_eqs = (int)p_mag->getNeq() + (int)p_grb->getNeq() + (int)p_ej->getNeq();
        (*p_log)(LOG_INFO,AT) << " ODE will solve"
                                        << " N_mag="<<p_mag->getNeq()
                                        << " N_grb="<<p_grb->getNeq()
                                        << " N_ej="<<p_ej->getNeq()
                                        << " (total " << p_pars->n_tot_eqs << ") equations. \n";
        m_InitData = new double [ p_pars->n_tot_eqs ];
        m_CurSol   = new double [ p_pars->n_tot_eqs ];
        // checks for the settings of BWs interaction prescription

        // chose the integrator for the system
        m_loglevel = loglevel;
        m_Method = integrator;
        p_Integrator = nullptr;
        switch (m_Method) {
            case Integrators::RK4 :
                p_Integrator = new IntegratorStatic<Integrators::rk4_integrator>
                        (Integrators::rk4_integrator(), loglevel);
                break;
            case Integrators::EULER:
                p_Integrator = new IntegratorStatic<Integrators::euler_integrator>
                        (Integrators::euler_integrator(), loglevel);
                break;
            case Integrators::MIDPOINT:
                p_Integrator = new IntegratorStatic<Integrators::midpoint_integrator>
                        (Integrators::midpoint_integrator(), loglevel);
                break;
            case Integrators::HEUN:
                p_Integrator = new IntegratorStatic<Integrators::heun3_integrator>
                        (Integrators::heun3_integrator(), loglevel);
                break;
            case Integrators::RALSTON:
                p_Integrator = new IntegratorStatic<Integrators::ralston3_integrator>
                        (Integrators::ralston3_integrator(), loglevel);
                break;
            case Integrators::DOP853:
                p_Integrator = new IntegratorStatic<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), loglevel);
                break;
            case Integrators::DOP853E:
                p_Integrator = new IntegratorEmbedded<Integrators::dop853_integrator>
                        (Integrators::dop853_integrator(), loglevel);
                break;
        }
    }
    Pars *& getPars(){ return p_pars; }
    ~EvolveODEsystem(){
        delete [] m_InitData;
        delete [] m_CurSol;
        delete p_pars;
        delete p_Integrator;
    }

    void setRelativePositions(){

        auto & jet_bws = p_pars->p_grb->getBWs();
        auto & ej_bws = p_pars->p_ej->getShells();

        size_t i_j_l;// = which_layer_to_use;//jet_bws.size() - 1;

        for (size_t il=0; il<ej_bws.size(); il++){
            for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++){
                auto & bw = ej_bws[il]->getBW(ish);
                i_j_l = bw->getPars()->which_jet_layer_to_use;
                double ej_ctheta = bw->getPars()->ctheta0;
                double j_ctheta = jet_bws[i_j_l]->getPars()->ctheta0;
                if (ej_ctheta < j_ctheta){
                    bw->getPars()->is_within0 = true;
                    bw->getPars()->ijl = i_j_l;
                }
            }
        }
#if 0
        for (size_t iej = 0; iej < ej_bws.size(); iej++){
            i_j_l = ej_bws[iej]->getPars()->which_jet_layer_to_use;
            double ej_ctheta = ej_bws[iej]->getPars()->ctheta0;
            double j_ctheta = jet_bws[i_j_l]->getPars()->ctheta0;
            if (ej_ctheta < j_ctheta){
                ej_bws[iej]->getPars()->is_within0 = true;
                ej_bws[iej]->getPars()->ijl = i_j_l;
            }
        }
#endif
    }

    /// fill the 'm_InitData' - vector of initial conditions for all blast waves
    void setInitialConditions( const double tb0 ){
        size_t ii = 0;
        // ***************************************
        if (p_pars->p_magnetar->run_magnetar){
            auto & magnetar = p_pars->p_magnetar;
            magnetar->setInitConditions(m_InitData, ii);
            ii += magnetar->getNeq();
        }
        // ***************************************
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t il = 0; il < jet_bws.size(); il++) {
//                auto &j_pars = jet_bws[il]->getPars();
//            double beta0 = EQS::Beta(j_pars->Gamma0);
//            double R0    = j_pars->tb0 * beta0 * CGS::c;
//            jet_bws[il]->getDensIsm()->evaluateRhoDrhoDrDefault(R0, j_pars->ctheta0);
//                std::cout << "Initializing blast wave (jet) layer=" << il << " (ii=" << ii << ")" << "\n";
                jet_bws[il]->setInitConditions(m_InitData, ii);
                ii += jet_bws[il]->getNeq();
            }
            double jet_extend = jet_bws[jet_bws.size()-1]->getPars()->ctheta0; // UNUSED
        }
        if ((p_pars->p_grb->run_jet_bws)&&(p_pars->p_ej->run_ej_bws))
            setRelativePositions();
        // ***************************************
        if (p_pars->p_ej->run_ej_bws) {
            auto & ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++){
                    ej_bws[il]->getBW(ish)->setInitConditions(m_InitData, ii);
                    ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
//            for (size_t il = 0; il < p_pars->n_ej_bws; il++) {
//                ej_bws[il]->setInitConditions(m_InitData, ii);
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        // ***************************************

        // ***************************************
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            m_CurSol[i] = m_InitData[i];
        }
        // **************************************
        if ( !isThereATermination() ){
            (*p_log)(LOG_ERR,AT)  <<" termination at initialization. Evolution canceled\n";
//            exit(1);
        }
        // **************************************
        if ( !isSolutionOk() ) {
            (*p_log)(LOG_ERR,AT)   << " Unphysical value in the initial data for evolution \n Exiting...";
            exit(1);
        }
        // **************************************
        for (size_t i = 0; i < p_pars->n_tot_eqs; i++){
            if (!std::isfinite(m_InitData[i])){
                (*p_log)(LOG_ERR,AT) << AT << "\n Nan in initial data: i="<<i<<" val="<<m_InitData[i]<<"\n";
//                exit(1);
            }
        }
        // **************************************
        p_Integrator->SetRHS( RHS_comb, p_pars->n_tot_eqs, p_pars );
        p_Integrator->SetInitCond( tb0, m_InitData );
        applyUnits();
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(0);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(0);
        // add electron properties (needed for synchron calculation)
//        addComputeForwardShockMicrophysics(0);
        is_initialized = true;

    }
    // evolve all blast waves
    void evolve( const double dx, const size_t ix ){
        if (!is_initialized){
            (*p_log)(LOG_ERR,AT)<<" bw set is not initialized. Cannot evolve\n";
            exit(1);
        }
        // eval relative geometry (which BW is behind what BW)
//        check_layer_relative_position( ix );
        // solve the ODE system for x_i = x_{i-1} + dx
        p_Integrator->Integrate( dx );
        // extract the solution vector from the ODE solver
        p_Integrator->GetY( m_CurSol );
        // apply units, e.g., energy is usually evolved in E/E0
        applyUnits();
        // ---
        Array & t_grid = p_pars->t_grid;
        if ( ix % 10 == 0 ) {
            // std::cout << p_pars->p_ej[0]->getCurrentIndexes() << "\n";
            (*p_log)(LOG_INFO,AT) << "it=" << ix << "/" << t_grid.size() << " t=" << t_grid[ix] << "\n";
        }
//        for(size_t i_ej = 0; i_ej < p_pars->p_bws_ej.size(); i_ej++){
//            auto & ej = p_pars->p_bws_ej[i_ej];
//            std::cout << " ish=" << ej->getPars()->ishell
//                      << " il="  << ej->getPars()->ishell
//                      << " G=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iGamma, ix))
//                      << " M2=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iM2, ix))
//                      << " R=" << string_format("%.2e", ej->getVal(DynRadBlastWave::Q::iR, ix))
//                      << " rho=" << string_format("%.2e", ej->getDensIsm()->m_rho)
//                      << " drho=" << string_format("%.2e", ej->getDensIsm()->m_drhodr)
//                      << " Grho=" << string_format("%.2e", ej->getDensIsm()->m_drhodr)
//                      << " Grel=" << string_format("%.2e", ej->getDensIsm()->m_GammaRel)
//                      << " dGrho=" << string_format("%.2e", ej->getDensIsm()->m_dGammaRhodR)
//                      << " dGrhodG=" << string_format("%.2e", ej->getDensIsm()->m_dGammaReldGamma)
//                      << " Grel=" << string_format("%.2e", ej->getDensIsm()->m_GammaRel)
//        }

        // check if there one of the layers is terminated
        if ( !isThereATermination() ){
            if(p_pars->i_restarts > 1){
                (*p_log)(LOG_ERR,AT)  << " second restart. Should not occure. Exiting...";
                exit(1);
            }
            (*p_log)(LOG_ERR,AT)  << " restarting iteration as there was a termination\n";
            evolve( dx, ix );
            p_pars->i_restarts += 1;
        }
        isThereLateralExpansionTermiantion();
        // check if there are no nans/unphysical vals in solution
        if ( !isSolutionOk() ) {
            (*p_log)(LOG_ERR,AT)  << " Unphysical value in the solution \n";
            exit(1);
        }
        p_pars->ix = ix;// latest solution
        p_pars->dx = dx;// latest solution
        p_pars->x = t_grid[ix];
        // plug the solution vector to the solution container of a given blast wave
        insertSolution(ix);
        // add other variables (that are not part of ODE but still needed)
        addOtherVariables(ix);
        // add electron properties (needed for synchron calculation)
//        addComputeForwardShockMicrophysics(ix);
    }
    inline auto * pIntegrator() { return p_Integrator; }

private:
    void applyUnits(){
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->applyUnits(m_CurSol, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->applyUnits(m_CurSol, ii);
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->applyUnits(m_CurSol, ii);
                    ii += ej_bws[il]->getBW(ish)->getNeq();
                }
            }
        }
//        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
//            ej_bws[i]->applyUnits( m_CurSol, ii);
//            ii += p_pars->n_eq_ej_bws;
//        }
    }
    bool isThereATermination(){
//        auto & magnetar = p_pars->p_magnetar;

        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToTerminate(m_CurSol, ii)) {
                    is_ok = false;
                    (*p_log)(LOG_ERR,AT) << " Terminating jet BW layer=" << i << " (of all) failed "
                              << " [ishell=" << jet_bws[i]->getPars()->ishell
                              << " ilayer=" << jet_bws[i]->getPars()->ilayer
                              << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                              << " ] \n";
                    jet_bws[i]->getPars()->end_evolution = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();

            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto & ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if (bw->isToTerminate(m_CurSol, ii)) {
                        is_ok = false;
                        bw->getPars()->end_evolution = true; // SET TO END
                        (*p_log)(LOG_ERR,AT) << " Terminating ejecta BW [ish=" << ish << " il=" << il << " "
                                  << " ii_eq=" << bw->getPars()->ii_eq
                                  << " ] \n";
                    }
                    ii += bw->getNeq();//p_pars->n_eq_ej_bws;
                }
            }


//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                if (ej_bws[i]->isToTerminate(m_CurSol, ii)) {
//                    is_ok = false;
//                    std::cerr  << " Terminating ejecta BW layer=" << i << " (of all) failed "
//                               << " [ishell=" << ej_bws[i]->getPars()->ishell
//                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
//                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
//                               << " ] \n";
//                    ej_bws[i]->getPars()->end_evolution = true; // SET TO END
//                }
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        return is_ok;
    }
    bool isThereLateralExpansionTermiantion(){
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_bws_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (jet_bws[i]->isToStopLateralExpansion(m_CurSol, ii)) {
                    is_ok = false;
                    jet_bws[i]->getPars()->end_spreading = true; // SET TO END
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        return is_ok;
    }
    bool isSolutionOk(){
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;
        size_t ii = 0; bool is_ok = true;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            if (!magnetar->isSolutionOk(m_CurSol, ii)){
                is_ok = false;
                (*p_log)(LOG_ERR,AT)  << " Magnetar evolution failed "
//                           << " [ishell=" << jet_bws[i]->getPars()->ishell
//                           << " ilayer=" << jet_bws[i]->getPars()->ilayer
//                           << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                           << " \n";
            }
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                if (( !jet_bws[i]->getPars()->end_evolution ) && (!jet_bws[i]->isSolutionOk(m_CurSol, ii))) {
                    is_ok = false;
                    (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of jet BW layer=" << i << " (of all) failed "
                               << " [ishell=" << jet_bws[i]->getPars()->ishell
                               << " ilayer=" << jet_bws[i]->getPars()->ilayer
                               << " ii_eq=" << jet_bws[i]->getPars()->ii_eq
                               << " ] \n";
                }
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto & bw = ej_bws[il]->getBW(ish);
                    if ((!bw->getPars()->end_evolution) && (!bw->isSolutionOk(m_CurSol, ii))) {
                        is_ok = false;
                        (*p_log)(LOG_ERR,AT)  << " Dyn. evol. of ejecta BW failed "
                                   << " [ishell=" << ish
                                   << " ilayer=" << il
                                   << " ii_eq=" << bw->getPars()->ii_eq
                                   << " ] \n";
                    }
                    ii += bw->getNeq();
                }
            }
//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                if ((!ej_bws[i]->getPars()->end_evolution) && (!ej_bws[i]->isSolutionOk(m_CurSol, ii))) {
//                    is_ok = false;
//                    std::cerr  << " Dyn. evol. of ejecta BW layer=" << i << " (of all) failed "
//                               << " [ishell=" << ej_bws[i]->getPars()->ishell
//                               << " ilayer=" << ej_bws[i]->getPars()->ilayer
//                               << " ii_eq=" << ej_bws[i]->getPars()->ii_eq
//                               << " ] \n";
//                }
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
        return is_ok;
    }
    void insertSolution(size_t it){
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->insertSolution(m_CurSol, it, ii);
            ii += magnetar->getNeq();
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->insertSolution(m_CurSol, it, ii);
                ii += jet_bws[i]->getNeq();
            }
        }
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il=0; il<ej_bws.size(); il++){
                for(size_t ish=0; ish<ej_bws[il]->nBWs(); ish++) {
                    auto &bw = ej_bws[il]->getBW(ish);
                    bw->insertSolution(m_CurSol, it, ii);
                    ii += bw->getNeq();
                }
            }
//            for (size_t i = 0; i < p_pars->n_ej_bws; i++) {
//                ej_bws[i]->insertSolution(m_CurSol, it, ii);
//                ii += p_pars->n_eq_ej_bws;
//            }
        }
    }
    void addOtherVariables(size_t it){
        if (p_pars->p_magnetar->run_magnetar) {
            p_pars->p_magnetar->addOtherVariables(it);
        }
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); i++) {
                jet_bws[i]->addOtherVars( it );
            }
        }
        auto & ej_bws = p_pars->p_ej;
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_bws = p_pars->p_ej->getShells();
            for (size_t il = 0; il < ej_bws.size(); il++) {
                for (size_t ish = 0; ish < ej_bws[il]->nBWs(); ish++) {
                    ej_bws[il]->getBW(ish)->addOtherVars(it);
                }
            }
        }


//        for(size_t i = 0; i < p_pars->n_ej_bws; i++){
//            p_pars->p_bws_ej[i]->addOtherVars( it );
//        }
    }
    /// get offset of the equations in the combined solution vector for a given blast wave from its ilayer and ishell
//    static inline size_t idx(std::unique_ptr<RadBlastWave> & bwPtr, void *pars){
//        auto * p_pars = (struct Pars *) pars;
//        return bwPtr->getNeq() * (bwPtr->getPars()->ilayer + p_pars->n_layers_j * bwPtr->getPars()->ishell);
//    }
    static void RHS_comb(double * out_Y,
                         size_t n_eq,
                         double x,
                         double const * Y,
                         void *pars){

        auto * p_pars = (struct Pars *) pars;
//        auto & magnetar = p_pars->p_magnetar;
//        auto & jet_bws = p_pars->p_bws_jet;
//        auto & ej_bws = p_pars->p_ej;

        /// advance magnetar to the next timestep
        size_t ii = 0;
        if (p_pars->p_magnetar->run_magnetar) {
            auto & magnetar = p_pars->p_magnetar;
            magnetar->evaluateRhs(out_Y, ii, x, Y);
            ii += magnetar->getNeq();
        }

        /// evaluate RHS for the jet (advance it to the next sub-step)
        if (p_pars->p_grb->run_jet_bws) {
            auto & jet_bws = p_pars->p_grb->getBWs();
            for (size_t i = 0; i < jet_bws.size(); ++i) {
                jet_bws[i]->evaluateRhs(out_Y, ii, x, Y);
//                ii += p_pars->n_eq_j_bws;
                ii += jet_bws[i]->getNeq();//n_eq_j_bws;
            }
        }
//        else
//            std::vector<std::unique_ptr<BlastWaveBase>> jet_bws ;


        /// evaluate RHS for the ejecta (advance it to the next sub-step)
        if (p_pars->p_ej->run_ej_bws) {
            auto &ej_layers = p_pars->p_ej->getShells();
            for (size_t il=0; il < ej_layers.size(); il++){
//                ej_layers[il]->evalRelativePosition(x, Y); // dEinjdt
                for(size_t ish=0; ish < ej_layers[il]->nBWs(); ish++) {
                    auto & ej_bw = ej_layers[il]->getBW(ish);
                    if (p_pars->p_grb->run_jet_bws) {
                        auto & jet_bws = p_pars->p_grb->getBWs();
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     & reinterpret_cast<std::vector<std::unique_ptr<BlastWaveBase>> &>(jet_bws),
                                                     p_pars->ix);
                    }
                    else
                        ej_bw->evaluateRhsDensModel2(out_Y, ii, x, Y,
                                                     NULL,
                                                     p_pars->ix);
                    ii += ej_bw->getNeq();
                }
            }


//            for(size_t iej = 0; iej < p_pars->n_ej_bws; iej++) {
////            if (ii!=ej_layers[iej]->getPars()->ii_eq){ // -------------- TO BE REMOVED
////                std::cerr<<AT << ii << " != "<< ej_layers[iej]->getPars()->ii_eq << "\n";
////                exit(1);
////            }
////            if (ej_layers[iej]->getPars()->end_evolution)
////                continue;
//                ej_layers[iej]->evaluateRhsDensModel2(out_Y, ii, x, Y,
//                     reinterpret_cast<std::vector<std::unique_ptr<BlastWaveBase>> &>(jet_bws), p_pars->ix);
////                ii += p_pars->n_eq_ej_bws;
//                ii += ej_layers[iej]->getNeq();
//            }
        }

        // ****************************************
        // Place to add interaction between blast waves
        // ****************************************

    }
    /// TODO add method for BW collision
//    static void collision(const double * Y, std::unique_ptr<RadBlastWave> & jetPtr,
//                          std::unique_ptr<RadBlastWave> & ejPtr, void *pars){
//
//        auto * p_pars = (struct Pars *) pars;
//        // ii = il + nlayers * ish
//        double gJ = Y[idx(jetPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ]; //Gamma;
////        double mJ=M2*M0 + M0;  double eJ=Eint2*M0;  double gAdiJ=gammaAdi;
//        double gEJ= Y[idx(ejPtr, pars) + DynRadBlastWave::Q_SOL::iGamma ];//Gamma_ej;
////        double mEJ=M2_ej*M0_ej+M0_ej; double eEJ=Eint2_ej*M0_ej; double gAdiEJ=gammaAdi_ej
////        double i_gM=Gamma_ej; double i_mM=M2_ej*M0_ej+M2*M0; double i_eM=Eint2_ej*M0_ej;
//        exit(1);
//    }
private:
    double * m_InitData{};
    double * m_CurSol{};
    Integrators::METHODS m_Method{};
    int m_loglevel{};
    IntegratorBase * p_Integrator;
    bool is_initialized = false;
//    Array m_t_grid;
};

#endif //SRC_BLASTWAVE_SET_H
