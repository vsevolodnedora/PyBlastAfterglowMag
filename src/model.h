//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_MODEL_H
#define SRC_MODEL_H

#include "pch.h"
#include "utils.h"
#include "ejecta_combined.h"
#include "observer.h"
#include "synchrotron_an.h"


void cast_times_freqs(Vector& lc_times, Vector& lc_freqs,
                      Vector& _times, Vector& _freqs,
                      bool is_one_to_one_already, std::unique_ptr<logger> & p_log){
    if (lc_times.empty() || lc_freqs.empty()){
        (*p_log)(LOG_ERR,AT)<<" empty time or freq arr.\n";
        exit(1);
    }
    if (is_one_to_one_already){
        if (lc_times.size()!=lc_freqs.size()){
            (*p_log)(LOG_ERR,AT)<<" size mismatch between arrays time and freq (for one-to-one freq-to-time)\n";
            exit(1);
        }

        _times = lc_times;
        _freqs = lc_freqs;

    }
    else {
        _times.resize(lc_freqs.size() * lc_times.size(), 0.0);
        _freqs.resize(lc_freqs.size() * lc_times.size(), 0.0);
        size_t ii = 0;
        for (double freq: lc_freqs) {
            for (double time: lc_times) {
                _times[ii] = time;
                _freqs[ii] = freq;
                ii++;
            }
        }
    }
}

class PyBlastAfterglow{
    struct Pars{
        double tb0{}; double tb1{}; int ntb{}; int iout{};
        Integrators::METHODS integrator = Integrators::METHODS::RK4;
        double rtol = 1e-5;
//        double jet_rtol = 1e-5;
//        double ej_rtol = 1e-5;
//        double jet_layer_fnu_stop_frac=1e-5;
        int nmax = 100000;
        int loglevel{};
        // ---
//        std::string method_eats = "";
//        LatStruct::METHOD_eats method_eats{};
        // ---
//        bool run_magnetar = false;
//        bool run_jet_bws = false;
//        bool run_ejecta_bws = false;
        // ---
//        bool is_jBW_init = false;
//        bool is_ejBW_init = false;
//        LatStruct::METHOD_eats jet_eats_method{};
//        LatStruct::METHOD_eats ejecta_eats_method{};
        // ---
//        bool is_jet_anal_synch_computed = false;
//        bool is_ejecta_anal_synch_computed = false;
//        bool is_jet_obs_pars_set = false;
//        bool is_ejecta_obs_pars_set = false;
//        bool is_jet_obsrad_pars_set = false;
//        bool is_ejecta_obsrad_pars_set = false;
//        bool is_jet_struct_set = false;
//        bool is_ejecta_struct_set = false;
        bool is_main_pars_set = false;
//        bool is_mag_pars_set = false;
        // ---

    };
    std::unique_ptr<logger> p_log;
    Pars * p_pars;
    std::unique_ptr<EvolveODEsystem> p_model;
//    LatStruct jetStruct{};
//    VelocityAngularStruct ejectaStructs{};
    std::unique_ptr<Output> p_out;
    std::unique_ptr<Magnetar> p_mag;
    std::unique_ptr<GRB> p_grb;
    std::unique_ptr<Ejecta> p_ej;
    std::unique_ptr<PWNset> p_ej_pwn;
    Vector _t_grid;
    Vector t_grid;
    int m_loglevel;
public:
    std::unique_ptr<Magnetar> & getMag(){return p_mag;}
    std::unique_ptr<PWNset> & getEjPWN(){return p_ej_pwn;}
    std::unique_ptr<GRB> & getGRB(){return p_grb;}
    std::unique_ptr<Ejecta> & getEj(){return p_ej;}
    Vector & getTburst(){return t_grid;}
    PyBlastAfterglow(int loglevel){
        m_loglevel = loglevel;
        p_pars = new Pars;
//        p_pars->loglevel = loglevel;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PyBlastAfterglow");
//        p_out = std::make_unique<Output>(loglevel);
        p_mag = std::make_unique<Magnetar>(loglevel);
        p_ej_pwn = std::make_unique<PWNset>(loglevel);
        p_grb = std::make_unique<GRB>(t_grid, loglevel);
        p_ej  = std::make_unique<Ejecta>(t_grid, loglevel);
        p_out = std::make_unique<Output>(loglevel);

    }
    ~PyBlastAfterglow() {
        (*p_log)(LOG_INFO,AT) << "Deleting PyBlastAfterglow instance...\n";
//        delete p_log;
        delete p_pars;
    }

    void setModelPars(StrDbMap pars, StrStrMap opts) {

        /// check if parameters present
        p_pars->tb0 = getDoublePar("tb0",pars,AT,p_log,-1,true);//pars.at("tb0");
        p_pars->tb1 = getDoublePar("tb1",pars,AT,p_log,-1,true);//(double) pars.at("tb1");
        p_pars->ntb = (int) getDoublePar("ntb",pars,AT,p_log,-1,true);//pars.at("ntb");
        p_pars->iout = (int) getDoublePar("iout",pars,AT,p_log,-1,true);//pars.at("ntb");
//        p_pars->ntb -= 1;

//        p_pars->loglevel = (int) getDoublePar("loglevel",pars,0,false);//pars.at("loglevel");
        p_pars->rtol = getDoublePar("rtol",pars,AT,p_log,1e-13,true);//(double) pars.at("rtol");
        p_pars->nmax = (int)getDoublePar("nmax",pars,AT,p_log,100000,false);//(double) pars.at("rtol");
//        p_pars->method_eats = getStrOpt("method_eats", opts, AT, p_log, "", true);
//        p_pars->jet_eats_method = LatStruct::setEatsMethod(p_pars->method_eats);

        // set options
        std::string opt = "integrator";
        Integrators::METHODS val;
        if (opts.find(opt) == opts.end()) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            val = Integrators::METHODS::RK4;
        }
        else {
            if (opts.at(opt) == "RK4")
                val = Integrators::METHODS::RK4;
            else if (opts.at(opt) == "DOP853")
                val = Integrators::METHODS::DOP853;
            else if (opts.at(opt) == "DOP853E")
                val = Integrators::METHODS::DOP853E;
            else {
                std::cerr << " option for: " << opt
                          << " given: " << opts.at(opt)
                          << " is not recognized \n" << " Possible options: "
                          << " RK4 " << " DOP853 " << " DOP853E " << "\n Exititng...";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->integrator = val;

//        p_pars->run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log,true);
//        p_pars->run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true);
//        p_pars->run_ejecta_bws = getBoolOpt("run_ejecta_bws", opts, AT,p_log,true);

        // ---------------------------------------------------------------
        /// Make a grid for ODE solver to follow
        _t_grid = TOOLS::MakeLogspaceVec(std::log10(p_pars->tb0),
                                         std::log10(p_pars->tb1),
                                         p_pars->ntb);
        /// Make a grid for solutions to be sotred
        if (p_pars->iout == 1)
            t_grid = _t_grid;
        else {
            t_grid.resize((int) p_pars->ntb / p_pars->iout);
            int k = 0;
            for (size_t i = 0; i < p_pars->ntb; i++)
                if (i % p_pars->iout == 0) {
                    t_grid[k] = _t_grid[i];
                    k++;
                }
        }
        p_pars->ntb = (int)t_grid.size();
//        std::cout << "_t_grd="<<_t_grid<<"\n";
//        std::cout << "t_grid="<<t_grid<<"\n";

        (*p_log)(LOG_INFO,AT) << "Computation tgrid = ["<<_t_grid[0]<<", "<<_t_grid[_t_grid.size()-1]<<"] n="<<_t_grid.size()<<"\n";
        (*p_log)(LOG_INFO,AT) << "Output      tgrid = ["<<t_grid[0]<<", "<<t_grid[t_grid.size()-1]<<"] n="<<t_grid.size()<<"\n";

        // -------------------------------------------------------------
        p_pars->is_main_pars_set = true;
//        std::cout << "setting model pars...\n";
        //std::cout << pars << "\n";
    }

//    void load

    /// run the time-evolution
    void run(){
        p_model = std::make_unique<EvolveODEsystem>( p_mag, p_grb, p_ej, p_ej_pwn,
                                                     t_grid, _t_grid, p_pars->integrator, m_loglevel );
        p_model->pIntegrator()->pPars()->rtol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->atol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->nmax = p_pars->nmax;
        p_model->setInitialConditions(p_pars->tb0);
        (*p_log)(LOG_INFO, AT) << "all initial conditions set. Starting evolution\n";
        // evolve
        double dt;
        int ixx = 1;
        for (size_t it = 1; it < _t_grid.size(); it++){
            p_model->getPars()->i_restarts = 0;
            dt = _t_grid[it] - _t_grid[it - 1];
            p_model->evolve(dt, it);
            if ((it % p_pars->iout == 0)) {
                p_model->storeSolution(ixx);
//                (*p_log)(LOG_INFO,AT)<<"Storing solution at i="<<ixx<<" t="<<t_grid[ixx]<<"\n";
                ixx++;
            }
        }
        (*p_log)(LOG_INFO, AT) << "evolution is completed\n";
    }

    /// save magnetar evolution
    void saveMagnetarEvolution(std::string workingdir, std::string fname, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Saving magnetar evolution...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"\n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & magnetar = p_mag;

        std::vector<std::vector<double>>  tot_mag_out{};
        std::vector<std::string> tot_names {};
        std::unordered_map<std::string,double> group_attrs{};

        auto & mag_v_ns = magnetar->m_vnames;
        auto t_arr = magnetar->getTbGrid(every_it);

        tot_mag_out.resize( mag_v_ns.size() );
        for (size_t ivar = 0; ivar < mag_v_ns.size(); ivar++){
            for (size_t it = 0; it < magnetar->getTbGrid().size(); it = it + every_it)
                tot_mag_out[ivar].emplace_back( (*magnetar)[ static_cast<Magnetar::Q>(ivar) ][it] );
        }

        std::unordered_map<std::string,double> attrs {};

        p_out->VectorOfVectorsH5(tot_mag_out, mag_v_ns, workingdir+fname, attrs);
    }

    /// save magnetar evolution
    void savePWNEvolution(std::string workingdir, std::string fname, size_t every_it, StrDbMap & main_pars, StrDbMap & pwn_pars){

        (*p_log)(LOG_INFO,AT) << "Saving PWN ...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & pwn = p_ej_pwn->getPWNs();
        std::vector<
                std::vector<
                        std::vector<double>>> tot_mag_out (pwn.size() );
        std::vector<std::string> tot_names {};
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        auto & dyn_v_ns = pwn[0]->m_vnames;
        auto t_arr = pwn[0]->getTbGrid(every_it);
        for (size_t i = 0; i < pwn.size(); i++) {
            tot_names.push_back("layer="+std::to_string(i));
            tot_mag_out[i].resize(dyn_v_ns.size() );
            for (size_t ivar = 0; ivar < dyn_v_ns.size(); ivar++){
                for (size_t it = 0; it < pwn[i]->getTbGrid().size(); it = it + every_it)
                    tot_mag_out[i][ivar].emplace_back((*pwn[i])[ static_cast<PWNmodel::Q>(ivar) ][it] );
            }
            ///write attributes
            auto & model = pwn[i];
            std::unordered_map<std::string,double> group_attr{
                    {"radius_w0",model->getPars()->radius_w0},
                    {"every_it",every_it}
            };
            group_attrs.emplace_back( group_attr );
        }
//    VecVector other_data { latStruct.cthetas0 };
//    std::vector<std::string> other_names { "cthetas0" };

        std::unordered_map<std::string, double> attrs{
                {"nlayers", pwn.size()}
        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: pwn_pars) { attrs[key] = value; }


        p_out->VecVectorOfVectorsAsGroupsH5(tot_mag_out, tot_names, dyn_v_ns,
                                            workingdir+fname, attrs, group_attrs);

    }

    /// save dynamical evolution of the jet blast-waves
    void saveJetBWsDynamics(std::string workingdir, std::string fname, size_t every_it, StrDbMap & main_pars, StrDbMap & grb_pars){
        (*p_log)(LOG_INFO,AT) << "Saving jet BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"\n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & models = p_grb->getBWs();
        std::vector<
                std::vector<
                        std::vector<double>>> tot_dyn_out (models.size() );
        std::vector<std::string> tot_names {};
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        auto & dyn_v_ns = models[0]->m_vnames;
        auto t_arr = models[0]->getTbGrid(every_it);
        for (size_t i = 0; i < models.size(); i++) {
            tot_names.push_back("layer="+std::to_string(i));
            tot_dyn_out[i].resize( dyn_v_ns.size() );
            for (size_t ivar = 0; ivar < dyn_v_ns.size(); ivar++){
                for (size_t it = 0; it < models[i]->getTbGrid().size(); it = it + every_it)
                    tot_dyn_out[i][ivar].emplace_back( (*models[i])[ static_cast<BlastWaveBase::Q>(ivar) ][it] );
            }
            ///write attributes
            auto & model = models[i];
            std::unordered_map<std::string,double> group_attr{
                    {"Gamma0",model->getPars()->Gamma0},
                    {"M0",model->getPars()->M0},
                    {"R0",model->getPars()->R0},
                    {"theta0",model->getPars()->theta_b0},
                    {"theta_max",model->getPars()->theta_max},
                    {"tb0",model->getPars()->tb0},
                    {"ijl",model->getPars()->ijl},
                    {"ncells",model->getPars()->ncells},
                    {"ilayer",model->getPars()->ilayer},
                    {"ishell",model->getPars()->ishell},
                    {"ctheta0",model->getPars()->ctheta0},
                    {"E0",model->getPars()->E0},
                    {"theta_c_l",model->getPars()->theta_c_l},
                    {"theta_c_h",model->getPars()->theta_c_h},
                    {"eps_rad",model->getPars()->eps_rad},
                    {"entry_time",model->getPars()->entry_time},
                    {"entry_r",model->getPars()->entry_r},
                    {"first_entry_r",model->getPars()->first_entry_r},
                    {"min_beta_terminate",model->getPars()->min_beta_terminate},
                    {"every_it",every_it}
            };
            group_attrs.emplace_back( group_attr );
        }
//    VecVector other_data { latStruct.cthetas0 };
//    std::vector<std::string> other_names { "cthetas0" };

        std::unordered_map<std::string, double> attrs{
                {"nlayers", models.size()}
        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: grb_pars) { attrs[key] = value; }


        p_out->VecVectorOfVectorsAsGroupsH5(tot_dyn_out, tot_names, dyn_v_ns,
                                            workingdir+fname, attrs, group_attrs);
    }

    /// save dynamical evolution of the ejecta blast-waves
    void saveEjectaBWsDynamics(std::string workingdir, std::string fname, size_t every_it, StrDbMap & main_pars, StrDbMap & ej_pars){
        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        auto & models = p_ej->getShells();

        std::vector<std::string> table_names;
        std::vector<std::vector<std::vector<double>>> tot_dyn_out ( nshells * nlayers );
        size_t i = 0;
        VecVector other_data;
        std::vector<std::string> other_names;
        auto & arr_names = models[0]->getBW(0)->m_vnames;//models[0][0].getBWdynPtr()->m_vnames;
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        for (size_t ishell = 0; ishell < nshells; ishell++){
            for(size_t ilayer = 0; ilayer < nlayers; ilayer++){
                table_names.push_back("shell="+std::to_string(ishell)+" layer="+std::to_string(ilayer));
                tot_dyn_out[i].resize( arr_names.size() );
                auto & bw = models[ilayer]->getBW(ishell);

                std::unordered_map<std::string,double> group_attr{
                        {"Gamma0",bw->getPars()->Gamma0},
                        {"M0",bw->getPars()->M0},
                        {"R0",bw->getPars()->R0},
                        {"theta0",bw->getPars()->theta_b0},
                        {"theta_max",bw->getPars()->theta_max},
                        {"tb0",bw->getPars()->tb0},
                        {"ijl",bw->getPars()->ijl},
                        {"ncells",bw->getPars()->ncells},
                        {"ilayer",bw->getPars()->ilayer},
                        {"ishell",bw->getPars()->ishell},
                        {"ctheta0",bw->getPars()->ctheta0},
                        {"E0",bw->getPars()->E0},
                        {"theta_c_l",bw->getPars()->theta_c_l},
                        {"theta_c_h",bw->getPars()->theta_c_h},
                        {"eps_rad",bw->getPars()->eps_rad},
                        {"entry_time",bw->getPars()->entry_time},
                        {"entry_r",bw->getPars()->entry_r},
                        {"first_entry_r",bw->getPars()->first_entry_r},
                        {"min_beta_terminate",bw->getPars()->min_beta_terminate}
                };
                group_attrs.emplace_back( group_attr );

                for (size_t ivar = 0; ivar < arr_names.size(); ivar++) {
                    for (size_t it = 0; it < bw->getTbGrid().size(); it = it + every_it)
                        tot_dyn_out[i][ivar].emplace_back( (*bw)[ static_cast<BlastWaveBase::Q>(ivar) ][it] );
                }
                i++;
            }
        }
        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells },
                {"nlayers", nlayers }
        };

        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VecVectorOfVectorsAsGroupsH5(tot_dyn_out, table_names, arr_names,
                                            workingdir+fname, attrs, group_attrs);
//    OutputVecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(fpath,tot_dyn_out, table_names,
//                                                         arr_names,other_data,other_names,attrs);
    }

    /// set parameters for analytical electron/synchrotron calculations
    void setPreComputeJetAnalyticElectronsPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Jet analytic electron pars...\n";

        if (!p_grb->run_jet_bws){
            std::cerr << " jet BWs were not evolved. Cannot compute electrons (analytic) \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_grb->getBWs();
        for (auto & model : models) {
//            model->setEatsPars(pars, opts);
//            model->getSynchAnPtr()->setMagPars( pars, opts );
            model->computeForwardShockElectronAnalyticVars();
            model->computeForwardShockSynchrotronAnalyticSpectrum();
        }
        p_grb->is_jet_anal_synch_computed = true;
    }
    void setPreComputeEjectaAnalyticElectronsPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Ejecta analytic electron pars...\n";

        if (!p_ej->run_ej_bws){
            std::cerr << " ejecta BWs were not evolved. Cannot compute electrons (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_ej->getShells();
        for (auto & model : models) {
            for (auto & bw : model->getBWs()) {
//            model->setEatsPars( pars, opts );
//            model->getSynchAnPtr()->setMagPars( pars, opts );
                bw->computeForwardShockElectronAnalyticVars();
                bw->computeForwardShockSynchrotronAnalyticSpectrum();
            }
        }
        p_ej->is_ejecta_anal_synch_computed = true;
    }

    /// compute and save comoving spectrum for all layers of the jet blast wave
    void computeSaveJetAnalyticSynchrotronSpectrum(std::string fpath, Vector & freq_array, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet analytic synchrotron spectrum...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr<< " jet analytic electrons were not evolved. Cannot compute spectrum (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_grb->getBWs();

        std::vector< // layers
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out{};
        auto t_arr = models[0]->getTbGrid(every_it);
        auto in_group_names = models[0]->getSynchAnPtr()->m_names_;
        std::vector<std::string> group_names{};
        for (size_t il = 0; il < models.size(); ++il ){
            group_names.emplace_back( "layer=" + std::to_string(il) );
            models[il]->computeForwardShockElectronAnalyticVars();
            out.emplace_back(models[il]->evalForwardShockComovingSynchrotron(freq_array, every_it) );
        }
        VecVector other_data{
                t_arr,
                freq_array
        };
        std::vector<std::string> other_names { "times", "freqs" };
        std::unordered_map<std::string,double> attrs{
                {"eps_e", models[0]->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", models[0]->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", models[0]->getSynchAnPtr()->getPars()->eps_t },
                {"p", models[0]->getSynchAnPtr()->getPars()->p },
                {"ksi_n", models[0]->getSynchAnPtr()->getPars()->ksi_n }
        };
//    auto out = { emissivity, absorption, em_th, em_pl, abs_th, abs_pl };
        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);
    }
    void computeSaveEjectaAnalyticSynchrotronSpectrum(std::string fpath, Vector & freq_array, size_t every_it){

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta analytic synchrotron spectrum...\n";

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot compute spectrum (analytic) \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_ej->getShells();

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        std::vector< // layers
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out {};
        auto t_arr = models[0]->getBW(0)->getTbGrid(every_it);
        auto in_group_names = models[0]->getBW(0)->getSynchAnPtr()->m_names_;
        std::vector<std::string> group_names{};
        size_t ish_il = 0;
        for (size_t ishell = 0; ishell < nshells; ishell++) {
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                out.emplace_back(
                        models[ilayer]->getBW(ishell)->evalForwardShockComovingSynchrotron(freq_array, every_it) );
                ish_il++;
            }
        }
        VecVector other_data{
                t_arr,
                freq_array
        };
        std::vector<std::string> other_names { "times", "freqs" };
        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells },
                {"nlayers", nlayers },
                {"eps_e", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
                {"p", models[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
                {"ksi_n", models[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
        };
//    auto out = { emissivity, absorption, em_th, em_pl, abs_th, abs_pl };
        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);
    }

    void updateJetObsPars(StrDbMap pars) {
        (*p_log)(LOG_INFO,AT) << " updating Jet observer pars...\n";
        size_t n_layers_j = p_grb->getBWs().size();//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        size_t n_layers_ej = (p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
        for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++) {
            auto &model = p_grb->getBWs()[ilayer];
            model->updateObsPars(pars);
        }
//        p_pars->is_jet_obs_pars_set = true;
    }
    void updateEjectaObsPars(StrDbMap pars) {

        auto & models = p_ej->getShells();

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        (*p_log)(LOG_ERR,AT) << "Updating Ejecta observer pars...\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells; ishell++) {
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                auto & model = p_ej->getShells()[ilayer]->getBW(ishell);//ejectaModels[ishell][ilayer];
                model->updateObsPars(pars);
                ii++;
            }
        }
//        p_pars->is_ejecta_obs_pars_set = true;
    }

    void computeSaveJetSkyImagesAnalytic(std::string workingdir,std::string fname, Vector times, Vector freqs,
                                         StrDbMap & main_pars, StrDbMap & grb_pars){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet sky image with analytic synchrotron...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr << "jet analytic electrons were not evolved. Cannot compute images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_grb->is_jet_obs_pars_set){
            std::cerr << "jet observer parameters are not set. Cannot compute image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells = 1; // jet has one shell only
        Image dummy(1);

        /// allocate memory for output data structure
        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};
        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(dummy.m_names.size());
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); ++i_vn) {
                    out[ii][i_vn].resize(nshells);
                }
                ii++;
            }
        }
        VecVector other_data{ times, freqs };
        std::vector<std::string> other_names { "times", "freqs" };

        /// allocate memory for image calculation
        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;
        Image image ( 2 * ncells );
        std::vector<Image> images ( nlayers );
        for (auto & im : images)
            im.resize( ncells );
        Image im_pj ( ncells );
        Image im_cj ( ncells );

        /// compute image for each time/freq
        ii = 0;
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            for (size_t it = 0; it < times.size(); ++it){
                im_pj.clearData();
                im_cj.clearData();
                for (auto & _image : images)
                    _image.clearData();

                (*p_log)(LOG_INFO,AT)<<" GRB SKYMAP "
                                     << " freq="<<freqs[ifreq] << " (" << ifreq << "/" << freqs.size()-1 << ") "
                                     << " time="<<times[it] << " (" << it << "/" << times.size()-1 << ") "
                                     << " theta_layers="<<nlayers
                                     << " ncells="<<ncells
                                     << "\n";


                size_t n_jet_empty_images = 0;
                std::vector<size_t> n_empty_images;
                size_t n_layers_j = nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
                for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++){
                    auto & model = p_grb->getBWs()[ilayer];
                    model->evalImagePW(images[ilayer], im_pj, im_cj, times[it], freqs[ifreq]);
                    if (images[ilayer].m_f_tot == 0){
                        n_jet_empty_images += 1;
                        n_empty_images.emplace_back(ilayer);
                    }
                }
                if (n_jet_empty_images > 0){
                    (*p_log)(LOG_INFO,AT)
                            << " grb empty images n="<<n_jet_empty_images<<"/"<<n_layers_j<<" layers.\t"
                            << " Layers with empty image: " << n_jet_empty_images << "\n";
                }

                /// combine all images from all layers into one image
                combineImages(image, ncells, nlayers, images);

                /// copy data from images into output structure
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells; ish++) {
                        out[ii][i_vn][ish] = image.m_data[i_vn];//arrToVec(image.m_data[i_vn]);
                    }
                }

                /// separately output total fluxes
                for (size_t ish = 0; ish < nshells; ish++) {
                    tota_flux[it] += image.m_f_tot;
                }
                ii++;
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );
        }

        std::vector<std::string> group_names{};
//        for (size_t it = 0; it < times.size(); it++){
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
//                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
//            }
//        }
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            for (size_t it = 0; it < times.size(); it++){
                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
            }
        }

        auto in_group_names = dummy.m_names;

        std::unordered_map<std::string,double> attrs{
                {"nlayers", nlayers}
        };
//                {"E0", p_grb->getBWs()[0]->getPars()->E0 },
//                {"Gamma0", p_grb->getBWs()[0]->getPars()->Gamma0 },
//                {"M0", p_grb->getBWs()[0]->getPars()->M0 },
//                {"R0", p_grb->getBWs()[0]->getPars()->R0 },
//                {"ctheta0", p_grb->getBWs()[0]->getPars()->ctheta0 },
//                {"theta_b0", p_grb->getBWs()[0]->getPars()->theta_b0 },
//                {"theta_w", p_grb->getBWs()[0]->getPars()->theta_w },
//                {"thetaObs", p_grb->getBWs()[0]->getEatsPars()->theta_obs },
//                {"d_L", p_grb->getBWs()[0]->getEatsPars()->d_l },
//                {"z", p_grb->getBWs()[0]->getEatsPars()->z },
//                {"eps_e", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_e },
//                {"eps_b", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_b },
//                {"eps_t", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_t },
//                {"p", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->p },
//                {"ksi_n", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->ksi_n }
//        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: grb_pars) { attrs[key] = value; }

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(workingdir+fname,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    void computeSaveEjectaSkyImagesAnalytic(std::string workingdir, std::string fname, Vector times, Vector freqs,
                                            StrDbMap & main_pars, StrDbMap & ej_pars){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta sky image with analytic synchrotron...\n";

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr  << "ejecta analytic electrons were not evolved. Cannot compute images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_ej->is_ejecta_obs_pars_set){
            std::cerr<< "ejecta observer parameters are not set. Cannot compute image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

//        auto & ejectaStruct = ejectaStructs;
//        size_t nshells = ejectaStruct.nshells;
        Image dummy(1,0);

        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};

        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(dummy.m_names.size());
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); ++i_vn) {
                    out[ii][i_vn].resize(nshells);
                }
                ii++;
            }
        }

        VecVector other_data{
                times,
                freqs
        };
        ii = 0;
        std::vector<std::string> other_names { "times", "freqs" };
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            VecVector total_flux_shell( nshells );
            for (auto & total_flux_shel : total_flux_shell)
                total_flux_shel.resize( times.size(), 0.0 );
            for (size_t it = 0; it < times.size(); ++it){
//                auto images = computeEjectaSkyMapPieceWise( times[it],freqs[ifreq]);
                std::vector<Image> images(nshells);
                computeEjectaSkyMapPieceWise( images, times[it],freqs[ifreq]);
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells; ish++) {
                        out[ii][i_vn][ish] = images[ish].m_data[i_vn];//arrToVec(images[ish].m_data[i_vn]);
                    }
                }
                for (size_t ish = 0; ish < nshells; ish++) {
                    tota_flux[it] += images[ish].m_f_tot;
                    total_flux_shell[ish][it] = images[ish].m_f_tot;
                }
                ii++;
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );

            for (size_t ish = 0; ish < nshells; ish++){
                other_data.emplace_back( total_flux_shell[ish] );
                other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) +
                                          " shell=" + string_format("%d", ish));
            }

        }

        std::vector<std::string> group_names{};
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            for (size_t it = 0; it < times.size(); it++){
                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
            }
        }


//        for (size_t it = 0; it < times.size(); it++){
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
//                                       + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
//            }
//        }

        auto in_group_names = dummy.m_names;

        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells},
                {"nshells", nlayers}
        };
//                {"thetaObs", p_ej->getShells()[0]->getBW(0)->getEatsPars()->theta_obs },
//                {"d_L", p_ej->getShells()[0]->getBW(0)->getEatsPars()->d_l },
//                {"z",  p_ej->getShells()[0]->getBW(0)->getEatsPars()->z },
//                {"eps_e",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
//                {"eps_b",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
//                {"eps_t",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
//                {"p",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
//                {"ksi_n",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
//        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(workingdir+fname,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    /*   */
    void computeSaveJetLightCurveAnalytic(std::string workingdir, std::string fname, std::string fname_layers,
                                          Vector lc_times, Vector lc_freqs, StrDbMap &main_pars, StrDbMap &grb_pars,
                                          bool lc_freq_to_time){

        Vector _times, _freqs;
        cast_times_freqs(lc_times,lc_freqs,_times,_freqs,lc_freq_to_time,p_log);

        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet light curve with analytic synchrotron...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr << "jet analytic electrons were not evolved. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_grb->is_jet_obs_pars_set){
            std::cerr << " jet observer parameters are not set. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;

//        auto & structure_hist = jetStruct;
        auto & tmp = p_grb->getBWs()[0]->getSynchAnPtr();
//        size_t nshells =structure_hist.nshells;
//        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;

        std::vector< // layers
                std::vector< // options
                        std::vector<double>>> // freqs * times
        out {};

        /// compute actual light curve for each layer
        auto light_curve = evalJetLightCurves( _times,_freqs);//[ilayer][it+ifreq]

        /// save total light curve for freqs
        int ii = 0;
        size_t n = _times.size();
        Vector total_fluxes (n,0.0);
        for(size_t ifnu = 0; ifnu < n; ifnu++){
            for (size_t ishil = 0; ishil < nlayers; ++ishil) {
                total_fluxes[ifnu] += light_curve[ishil][ifnu];
            }
        }
        std::vector<std::string> other_names { "times", "freqs", "total_fluxes" };
        VecVector out_data {_times, _freqs, total_fluxes};

        StrDbMap attrs {{"nlayers", nlayers}};
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: grb_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);

        /// save layer-by-layer data
        if (fname_layers == "none")
            return;

        light_curve.emplace_back(_times);
        light_curve.emplace_back(_freqs);
        light_curve.emplace_back(total_fluxes);

        std::vector<std::string> group_names;
        for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
            group_names.emplace_back("layer=" + std::to_string(ilayer));
        }
        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_fluxes");
        p_out->VectorOfVectorsH5(light_curve, group_names, workingdir+fname,  attrs);
    }

    void computeSaveEjectaLightCurveAnalytic(std::string workingdir,std::string fname, std::string fname_shells_layers,
                                             Vector lc_times, Vector lc_freqs, StrDbMap & main_pars, StrDbMap & ej_pars,
                                             bool lc_freq_to_time){

        Vector _times, _freqs;
        cast_times_freqs(lc_times,lc_freqs,_times,_freqs,lc_freq_to_time,p_log);

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta light curve with analytic synchrotron...\n";

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_ej->is_ejecta_obs_pars_set){
            std::cerr << " ejecta observer parameters are not set. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        auto & tmp = p_ej->getShells()[0]->getBW(0)->getSynchAnPtr();

        std::vector< // layers / shells
                std::vector< // options
                        std::vector<double>>> // freqs*times
        out {};

        /// evaluate light curve
        auto light_curve = evalEjectaLightCurves( _times, _freqs);

        /// save total lightcurve
        size_t n = _times.size();
        Vector total_fluxes (n, 0.0);
        for (size_t itnu = 0; itnu < n; ++itnu) {
            size_t ishil = 0;
            for (size_t ishell = 0; ishell < nshells; ++ishell) {
                for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
                    total_fluxes[itnu] += light_curve[ishell][ilayer][itnu];
                    ishil++;
                }
            }
        }
        std::vector<std::string> other_names { "times", "freqs", "total_fluxes" };
        VecVector out_data {_times, _freqs, total_fluxes};

        std::unordered_map<std::string,double> attrs{ {"nshells", nshells}, {"nlayers", nlayers} };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);


        /// save light curve for each shell and layer
        if (fname_shells_layers == "none")
            return;
        std::vector<std::string> group_names;
        VecVector total_fluxes_shell_layer(nshells*nlayers);
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells; ++ishell) {
            for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                for (size_t ifnu = 0; ifnu < n; ifnu++){
                    total_fluxes_shell_layer[ii][ifnu] = light_curve[ishell][ilayer][ifnu];
                }
                ii++;
            }
        }
        total_fluxes_shell_layer.emplace_back(_times);
        total_fluxes_shell_layer.emplace_back(_freqs);
        total_fluxes_shell_layer.emplace_back(total_fluxes);

        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_fluxes");
        p_out->VectorOfVectorsH5(total_fluxes_shell_layer, group_names, workingdir+fname,  attrs);
    }

private:

    void computeEjectaSkyMapPieceWise( std::vector<Image> & images, double obs_time, double obs_freq ){

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

//        std::vector<Image> images;
//        images.resize(ejectaStructs.nshells);
        if (images.empty()){
            (*p_log)(LOG_ERR,AT) << " empty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nshells){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }
        // ejectaStructs.nshells);
        size_t ii = 0;
        size_t n_jet_empty_images = 0;
//        std::vector<size_t> n_empty_images_layer;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        std::vector<Image> tmp (nlayers);
        for (auto & _tmp : tmp)
            _tmp.resize( ncells );
        Image tmp_pj( ncells); // std::vector<Image> tmp_pj (ejectaStructs.structs[0].nlayers);
        Image tmp_cj( ncells); // std::vector<Image> tmp_cj (ejectaStructs.structs[0].nlayers);
        for (size_t ishell = 0; ishell <  nshells; ishell++){
            for (auto & _tmp : tmp)
                _tmp.clearData();
//            tmp[].clearData();
            tmp_pj.clearData(); tmp_cj.clearData();
//            image_i.clearData();
//            im_pj.clearData();
//            im_cj.clearData();
//            auto & ejectaStruct = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStruct.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
//            std::vector<Image> tmp;
//            std::vector<Image> tmp (n_layers_ej);
            std::vector<size_t> n_empty_images_layer;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
                auto & model = p_ej->getShells()[ ilayer ];
//                tmp.emplace_back( model->evalImagePW(obs_time, obs_freq) );
                model->getBW(ishell)->evalImagePW(tmp[ilayer], tmp_pj, tmp_cj, obs_time, obs_freq);
                if (tmp[ilayer].m_f_tot == 0){
                    n_jet_empty_images += 1;
                    n_empty_images_layer.emplace_back(ilayer);
                }
                ii ++ ;
            }
            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }

//            images.emplace_back( combineImages(ejectaStruct,tmp) );
            combineImages(images[ishell], ncells, nlayers, tmp) ;
        }

        /// print which layers/shells gave empty image
        if (p_log->getLogLevel() == LOG_INFO) {
            if (n_jet_empty_images > 0) {
                auto &ccerr = std::cout;
                ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an empty images for total n="
                      << n_jet_empty_images << " layers. Specifically:\n";
                for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++) {
//                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                size_t n_layers_ej = nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
                    ccerr << "\t [ishell=" << n_empty_images_shells[ish] << " ilayer] = [ ";
                    for (size_t il = 0; il < n_empty_images[ish].size(); il++) {
                        ccerr << n_empty_images[ish][il] << " ";
                    }
                    ccerr << "] / (" << nlayers << " total layers) \n";
                }
            }
        }

//        return std::move( images );
    }
    VecVector evalJetLightCurves( Vector & obs_times, Vector & obs_freqs ){
        (*p_log)(LOG_INFO,AT)<<" starting grb light curve calculation\n";
        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;
//        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
        VecVector light_curves(nlayers); // [i_layer][i_time]
        for (auto & arr : light_curves)
            arr.resize(obs_times.size(), 0.);
        double flux;
        double rtol = p_grb->jet_rtol;
//        Image image;
        Image image_i ( ncells ); Image im_pj ( ncells ); Image im_cj ( ncells );
        for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
            auto & model = p_grb->getBWs()[ilayer];
            (*p_log)(LOG_INFO,AT)<<" GRB LC ntimes="<<obs_times.size()<<" theta_layer="<<ilayer<<"/"<<nlayers<<
                                 " phi_cells="<<LatStruct::CellsInLayer(model->getPars()->ilayer)<<"\n";
            model->evalForwardShockLightCurve(image_i, im_pj, im_cj, light_curves[ilayer], obs_times, obs_freqs);
//            for (size_t it = 0; it < obs_times.size(); it++) {
//                if (model->getPars()->m_method_eats == LatStruct::i_pw) {
////                    auto image = model->evalImagePW(obs_times[it], obs_freqs[it]);
//                    model->evalImagePW(image, obs_times[it], obs_freqs[it]);
//                    light_curves[ilayer][it] += image.m_f_tot;
//                }
//                else{
//                    double atol = light_curves[ilayer][it] * rtol / (double)n_layers_j;
//                    light_curves[ilayer][it] += model->evalFluxDensA(obs_times[it], obs_freqs[it], atol);
//                }
//            }
        }
        return std::move( light_curves );
    }
    std::vector<VecVector> evalEjectaLightCurves( Vector & obs_times, Vector & obs_freqs){
        (*p_log)(LOG_INFO,AT)<<" starting ejecta light curve calculation\n";
        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();
        std::vector<VecVector> light_curves(nshells); // [ishell][i_layer][i_time]
        for (auto & arr : light_curves){
            size_t n_layers_ej = nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
            arr.resize(n_layers_ej);
            for (auto & arrr : arr){
                arrr.resize( obs_times.size(), 0. );
            }
        }
        double flux_pj, flux_cj; size_t ii = 0;
//        Image image;
        double rtol = p_ej->ej_rtol;
        Image image_i ( ncells );
        Image im_pj ( ncells );
        Image im_cj ( ncells );
        for (size_t ishell = 0; ishell < nshells; ishell++){
            image_i.clearData();
            im_pj.clearData();
            im_cj.clearData();
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStructs.structs[ishell].nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                auto & model = p_ej->getShells()[ilayer];//ejectaModels[ishell][ilayer];
//                model->setEatsPars( pars, opts );
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC ntimes="<<obs_times.size()
                        << " vel_shell="<<ishell<<"/"<<nshells-1
                        << " theta_layer="<<ilayer<<"/"<<nlayers
                        << " phi_cells="<<LatStruct::CellsInLayer(ilayer)<<"\n";
                model->getBW(ishell)->evalForwardShockLightCurve(image_i, im_pj, im_cj,
                                                                     light_curves[ishell][ilayer], obs_times, obs_freqs);

//                for (size_t it = 0; it < obs_times.size(); it++) {
//                    if (model->getPars()->m_method_eats == LatStruct::i_pw) {
//                        //                    auto tmp = model->evalImagePW(obs_times[it], obs_freqs[it] );
//                        model->evalImagePW(image, obs_times[it], obs_freqs[it] );
//                        light_curves[ishell][ilayer][it] += image.m_f_tot;
//                    }
//                    else{
//                        double atol = light_curves[ishell][ilayer][it] * rtol / (double)n_layers_ej;
//                        light_curves[ishell][ilayer][it] += model->evalFluxDensA(obs_times[it], obs_freqs[it], atol);
//                    }
//                }
                ii ++;
            }
        }
        return std::move( light_curves );
    }
private:
};




















#if 0
class PyBlastAfterglow{
    struct Pars{
        double tb0{}; double tb1{}; int ntb{};
        Integrators::METHODS integrator = Integrators::METHODS::RK4;
        double rtol = 1e-5;
//        double jet_rtol = 1e-5;
//        double ej_rtol = 1e-5;
//        double jet_layer_fnu_stop_frac=1e-5;
        int nmax = 100000;
        int loglevel = CurrLogLevel;
        // ---
//        std::string method_eats = "";
//        LatStruct::METHOD_eats method_eats{};
        // ---
//        bool run_magnetar = false;
//        bool run_jet_bws = false;
//        bool run_ejecta_bws = false;
        // ---
//        bool is_jBW_init = false;
//        bool is_ejBW_init = false;
//        LatStruct::METHOD_eats jet_eats_method{};
//        LatStruct::METHOD_eats ejecta_eats_method{};
        // ---
//        bool is_jet_anal_synch_computed = false;
//        bool is_ejecta_anal_synch_computed = false;
//        bool is_jet_obs_pars_set = false;
//        bool is_ejecta_obs_pars_set = false;
//        bool is_jet_obsrad_pars_set = false;
//        bool is_ejecta_obsrad_pars_set = false;
//        bool is_jet_struct_set = false;
//        bool is_ejecta_struct_set = false;
        bool is_main_pars_set = false;
//        bool is_mag_pars_set = false;
        // ---

    };
    std::unique_ptr<logger> p_log;
    Pars * p_pars;
    std::unique_ptr<EvolveODEsystem> p_model;
//    LatStruct jetStruct{};
//    VelocityAngularStruct ejectaStructs{};
    std::unique_ptr<Output> p_out;
    std::unique_ptr<Magnetar> p_mag;
    std::unique_ptr<GRB> p_grb;
    std::unique_ptr<Ejecta> p_ej;
    Array t_grid;
public:
    std::unique_ptr<Magnetar> & getMag(){return p_mag;}
    std::unique_ptr<GRB> & getGRB(){return p_grb;}
    std::unique_ptr<Ejecta> & getEj(){return p_ej;}

    PyBlastAfterglow(int loglevel){
        p_pars = new Pars;
//        p_pars->loglevel = loglevel;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PyBlastAfterglow");
//        p_out = std::make_unique<Output>(loglevel);
        p_mag = std::make_unique<Magnetar>(t_grid, loglevel);
        p_grb = std::make_unique<GRB>(t_grid, loglevel);
        p_ej = std::make_unique<Ejecta>(t_grid, loglevel);

    }
    ~PyBlastAfterglow() {
        std::cout << "Deleting PyBlastAfterglow instance...\n";
//        delete p_log;
        delete p_pars;
    }

    /// set jet lateral structure
//    static auto listParsUniformJetStruct(){ return LatStruct::listParametersAnalyticBlastWave(); }
//    static auto listParametersAnalyticBlastWave(){ return LatStruct::listParametersAnalyticBlastWave(); }
//    static auto listJetStructOpts(){ return LatStruct::listStructOpts(); }
//    void setJetStructAnalytic(StrDbMap pars, StrStrMap opts){
//        std::string jet_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
//        p_pars->jet_eats_method = LatStruct::setEatsMethod(jet_eats_method);
//        jetStruct.initAnalytic( pars, opts, jet_eats_method, p_pars->loglevel );
//        p_pars->is_jet_struct_set = true;
//    }
//
////    static auto listParsNumericJetStruct(){ return LatStruct::listParsCustomStruct(); }
//    void setJetStructNumeric( Vector & dist_thetas, Vector & dist_EEs, Vector & dist_Gam0s, Vector & dist_MM0s,
//                              bool force_grid, std::string eats_method ){
//        jetStruct.initCustom( dist_thetas, dist_EEs, dist_Gam0s, dist_MM0s, force_grid, eats_method,
//                              p_pars->loglevel);
//        p_pars->is_jet_struct_set = true;
//    }


//    /// set ejecta lateral & velocity structure
////    static auto listParsAnalyticEjectaStruct(){ std::cerr << AT << " not implemented\n"; exit(1); }
//    void setEjectaStructAnalytic(StrDbMap pars, StrStrMap opts){
//        std::cerr << " not implimeneted\n Exiting...";
//        std::cerr << AT << "\n";
//        exit(1);
//    }
////    static std::vector<std::string> listParsNumericEjectaStruct(){ return VelocityAngularStruct::list_pars_v_ns(); }
//    void setEjectaStructNumericUniformInTheta(Vector & dist_thetas0, Vector & dist_betas, Vector & dist_ek,
//                                              size_t nlayers, double mfac, StrStrMap & opts){
//        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
//        p_pars->ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
//        ejectaStructs.initUniform(dist_thetas0,dist_betas,dist_ek, nlayers,mfac,
//                                  ej_eats_method, p_pars->loglevel);
//        p_pars->is_ejecta_struct_set = true;
//    }
//    void setEjectaStructNumeric(Vector dist_thetas, Vector dist_betas, VecVector dist_ek, bool force_grid,
//                                StrStrMap & opts){
//        std::string ej_eats_method = getStrOpt("method_eats",opts,AT,p_log,"", true);
//        p_pars->ejecta_eats_method = LatStruct::setEatsMethod(ej_eats_method);
//        ejectaStructs.initCustom(dist_thetas, dist_betas, dist_ek, force_grid,
//                                 ej_eats_method, p_pars->loglevel);
//        p_pars->is_ejecta_struct_set = true;
//    }


    /// list and set parameters that the whole model has (time grid, solver, general settings)
//    static std::vector<std::string> list_pars(){ return { "tb0", "tb1", "ntb", "rtol"}; }
//    static std::vector<std::string> list_opts(){ return { "integrator", "run_jet_bws", "run_ejecta_bws" }; }; //, "jet_eats", "ejecta_eats"
    void setModelPars(StrDbMap pars, StrStrMap opts) {

        /// check if parameters present
        p_pars->tb0 = getDoublePar("tb0",pars,AT,p_log,-1,true);//pars.at("tb0");
        p_pars->tb1 = getDoublePar("tb1",pars,AT,p_log,-1,true);//(double) pars.at("tb1");
        p_pars->ntb = (int) getDoublePar("ntb",pars,AT,p_log,-1,true);//pars.at("ntb");
//        p_pars->loglevel = (int) getDoublePar("loglevel",pars,0,false);//pars.at("loglevel");
        p_pars->rtol = getDoublePar("rtol",pars,AT,p_log,1e-13,true);//(double) pars.at("rtol");
        p_pars->nmax = (int)getDoublePar("nmax",pars,AT,p_log,100000,false);//(double) pars.at("rtol");
//        p_pars->method_eats = getStrOpt("method_eats", opts, AT, p_log, "", true);
//        p_pars->jet_eats_method = LatStruct::setEatsMethod(p_pars->method_eats);

        // set options
        std::string opt = "integrator";
        Integrators::METHODS val;
        if (opts.find(opt) == opts.end()) {
            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
            val = Integrators::METHODS::RK4;
        }
        else {
            if (opts.at(opt) == "RK4")
                val = Integrators::METHODS::RK4;
            else if (opts.at(opt) == "DOP853")
                val = Integrators::METHODS::DOP853;
            else if (opts.at(opt) == "DOP853E")
                val = Integrators::METHODS::DOP853E;
            else {
                std::cerr << " option for: " << opt
                          << " given: " << opts.at(opt)
                          << " is not recognized \n" << " Possible options: "
                          << " RK4 " << " DOP853 " << " DOP853E " << "\n Exititng...";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        p_pars->integrator = val;

//        p_pars->run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log,true);
//        p_pars->run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true);
//        p_pars->run_ejecta_bws = getBoolOpt("run_ejecta_bws", opts, AT,p_log,true);

        // ---------------------------------------------------------------
        // place all blast wave into the "evolver"
        t_grid = TOOLS::MakeLogspace(log10(p_pars->tb0), log10(p_pars->tb1),p_pars->ntb);

        // -------------------------------------------------------------
        p_pars->is_main_pars_set = true;
//        std::cout << "setting model pars...\n";
        //std::cout << pars << "\n";
    }

    ///
//    void setMagnetarPars(StrDbMap pars, StrStrMap opts){
////        bool run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log,true);
//        if (!p_pars->is_main_pars_set){
//            std::cerr << "\n model parameters were not set. Cannot run model. \n Exiting...\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
//        p_magnetar = std::make_unique<Magnetar>(t_grid, p_pars->loglevel);
////        p_model->getPars()->p_magnetar->setPars(pars, opts);
//        p_magnetar->setPars(pars, opts);
////        p_magnetar->is_mag_pars_set = true;
//    }

    /// list and set parameters that each blast wave needs (dynamics); interaction with others
//    static auto listBwPars(){ return listBwPars_(); }
//    static auto listBwOpts(){ return listBwOpts_(); }
//    void setJetBwPars(StrDbMap pars, StrStrMap opts){
//        bool run_jet_bws = getBoolOpt("run_jet_bws", opts, AT,p_log,true);
//        if (!p_pars->is_main_pars_set){
//            std::cerr << "\n model parameters were not set. Cannot run model. \n Exiting...\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
//
//        size_t ii_eq = 0;
//
//        if (p_pars->run_magnetar)
//            ii_eq += p_magnetar->getNeq();
//
//        size_t n_layers = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        for(size_t i = 0; i < n_layers; i++){
////            DynRadBlastWave x(t_grid, 0, i, p_pars->loglevel);
//            p_bws_jet.emplace_back( std::make_unique<DynRadBlastWave>(t_grid, 0, i, p_pars->loglevel) );
//            setAllParametersForOneLayer(jetStruct, *(p_bws_jet[i]), pars, opts, i, ii_eq);
//            p_bws_jet[i]->setEatsPars(pars, opts);
//            p_bws_jet[i]->getSynchAnPtr()->setPars( pars, opts );
//            ii_eq += p_bws_jet[i]->getNeq();
//        }
//
//        /// parameters for EATS methods
////        p_pars->jet_eats_method = LatStruct::setEatsMethod(
////                getStrOpt("method_eats",opts,AT,p_log,"",true));
//        p_pars->jet_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);
//        p_pars->jet_layer_fnu_stop_frac= getDoublePar("fnu_min_frac",pars,AT,p_log,-1, true);
//
//        (*p_log)(LOG_INFO,AT) << "finished initializing jet...\n";
//        p_pars->is_jet_obs_pars_set = true;
//        p_pars->is_jBW_init = true;
//    }

//    void setEjectaBwPars(StrDbMap pars, StrStrMap opts){
//        bool run_ejecta_bws = getBoolOpt("run_ejecta_bws", opts, AT,p_log,true);
//        if (!p_pars->is_main_pars_set){
//            std::cerr<< "\n model parameters were not set. Cannot run model. \n Exiting...\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
//
//        /// counter for total equations to be solved with ODE solver
//        size_t ii_eq = 0;
//
//        /// offset counter for magnetar
//        if (p_pars->run_magnetar)
//            ii_eq += p_magnetar->getNeq();
//
//        /// offset counter for grb
//        size_t n_layers_jet = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        if (p_pars->run_jet_bws)
//            for(size_t i = 0; i < n_layers_jet; i++){
//                ii_eq += p_bws_jet[i]->getNeq();
//            }
//
////        size_t n_layers_i = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
////        size_t ii = 0;
//        bool is_within = false;
//        std::vector<size_t> which_within{};
//        size_t n_ejecta_empty_images = 0;
//        std::vector<std::vector<size_t>> n_empty_images;
//        std::vector<size_t> n_empty_images_shells;
//        size_t nshells = ejectaStructs.nshells;
//        size_t n_layers_ej = ejectaStructs.structs[0].nlayers;
//        std::vector<std::vector<size_t>> n_empty_images_layer_shell;
//        for (auto & n_empty_images_layer : n_empty_images_layer_shell)
//            n_empty_images_layer.resize(nshells);
//        for(size_t il = 0; il < n_layers_ej; il++){
//            p_ej.push_back( std::make_unique<CumulativeShell>(t_grid, nshells, il, p_log->getLogLevel()) );
//            for (size_t ish = 0; ish < nshells; ish++){
//                auto & bw = p_ej[il]->getBW(ish);
//                auto & struc = ejectaStructs.structs[ish];
//                setAllParametersForOneLayer(struc, * bw, pars, opts, il, ii_eq);
///// Override the layer-to-use
//                if (bw->getPars()->which_jet_layer_to_use == 0){
//                    bw->getPars()->which_jet_layer_to_use = 0; // the fastest
//                }
//                else if(n_layers_ej==0){
//                    n_ejecta_empty_images += 1;
//                    n_empty_images_layer_shell[ish].emplace_back(il);
////                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
//                }
//                else if(n_layers_ej == 0){
//                    // NO jet structure was set, so exiting I guess... :)
//                    // TODO THIS MIGHT BE WRONG -- why 'n_layers_i'
//                }
//                else if(n_layers_jet == 0){
//                    // NO jet structure was set, so exiting I guess... :)
//                }
//                else if ((bw->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
//                    bw->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
//                }
//                else if ((bw->getPars()->which_jet_layer_to_use < n_layers_jet) &&
//                         (bw->getPars()->which_jet_layer_to_use > -1)){
//                    //
//                }
//                else{
//                    std::cerr << " which_jet_layer_to_use="<<bw->getPars()->which_jet_layer_to_use
//                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
//                              <<" for the slowest"
//                              <<" or any N in between the two for a specific jet layer \n"
//                              << "Exiting..."
//                              << "\n";
//                    std::cerr << AT << "\n";
//                    exit(1);
//                }
//                ii_eq += bw->getNeq();
//
//                bw->setEatsPars(pars, opts);
//                bw->getSynchAnPtr()->setPars( pars, opts );
//
////                ii++;
//            }
//        }
//
//
//
//# if 0
//        size_t ii = 0;
//        bool is_within = false;
//        std::vector<size_t> which_within{};
//
//        size_t n_ejecta_empty_images = 0;
//        std::vector<std::vector<size_t>> n_empty_images;
//        std::vector<size_t> n_empty_images_shells;
//        size_t n_layers_i = ejectaStructs.structs[0].nlayers;
//
//        for(size_t i = 0; i < ejectaStructs.nshells; i++){
//            std::vector<size_t> n_empty_images_layer;
//            auto & struc = ejectaStructs.structs[i];
//            if (n_layers_i != struc.nlayers){
//                (*p_log)(LOG_ERR,AT)<<" expected nlayers="<<n_layers_i<<" (from 0th shell) got="<<struc.nlayers<<"\n";
//                exit(1);
//            }
////            size_t n_layers_i = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
//            for(size_t j = 0; j < n_layers_i; j++) { // ii = il + nlayers * ish
//                p_bws_ej.emplace_back(std::make_unique<DynRadBlastWave>(t_grid, i, j, p_pars->loglevel));
//                setAllParametersForOneLayer(struc, *(p_bws_ej[ii]), pars, opts, j, ii_eq);
//                /// Override the layer-to-use
//                if (p_bws_ej[ii]->getPars()->which_jet_layer_to_use == 0){
//                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = 0; // the fastest
//                }
//                else if(n_layers_i==0){
//                    n_ejecta_empty_images += 1;
//                    n_empty_images_layer.emplace_back(j);
////                    std::cerr << AT << "\n jet structure was NOT initialized. No layer selected for ejecta to propagate through.\n";
//                }
//                else if(n_layers_i == 0){
//                    // NO jet structure was set, so exiting I guess... :)
//                    // TODO THIS MIGHT BE WRONG -- why 'n_layers_i'
//                }
//                else if(n_layers_jet == 0){
//                    // NO jet structure was set, so exiting I guess... :)
//                }
//                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use > n_layers_jet - 1)){
//                    p_bws_ej[ii]->getPars()->which_jet_layer_to_use = (int)n_layers_jet - 1;
//                }
//                else if ((p_bws_ej[ii]->getPars()->which_jet_layer_to_use < n_layers_jet) &&
//                         (p_bws_ej[ii]->getPars()->which_jet_layer_to_use > -1)){
//                    //
//                }
//                else{
//                    std::cerr << " which_jet_layer_to_use="<<p_bws_ej[ii]->getPars()->which_jet_layer_to_use
//                              << "\n" << " expected 0 (for fasterst) or any N larger than n_layers_jet=" << (int)n_layers_jet-1
//                              <<" for the slowest"
//                              <<" or any N in between the two for a specific jet layer \n"
//                              << "Exiting..."
//                              << "\n";
//                    std::cerr << AT << "\n";
//                    exit(1);
//                }
//                ii_eq += p_bws_ej[ii]->getNeq();
//
//                p_bws_ej[ii]->setEatsPars(pars, opts);
//                p_bws_ej[ii]->getSynchAnPtr()->setPars( pars, opts );
//
//                ii++;
//            }
//
//            if(!n_empty_images_layer.empty()){
//                n_empty_images_shells.emplace_back(i);
//                n_empty_images.emplace_back(n_empty_images_layer);
//            }
//        }
//#endif
//        p_pars->is_ejBW_init = true;
//        p_pars->is_ejecta_obs_pars_set = true;
//
//        if (n_ejecta_empty_images > 0){
//            auto & ccerr = std::cout;
//            ccerr << "Ejecta blastwave is NOT initialized for total n="
//                  << n_ejecta_empty_images << " layers. Specifically:\n";
//            for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++){
//                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                size_t n_layers_i = ejectaStruct.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
//                ccerr << "\t [ishell="<<n_empty_images_shells[ish] << " ilayer] = [";
//                for (size_t il = 0; il < n_empty_images[ish].size(); il++){
//                    ccerr << n_empty_images[ish][il] << " ";
//                }
//                ccerr << "] / (" << n_layers_i << " total layers) \n";
//            }
//        }
//
//        p_pars->jet_rtol = getDoublePar("rtol_adapt",pars,AT,p_log,-1, true);
//
//        (*p_log)(LOG_INFO,AT) << "finished initializing ejecta...\n";
//    }

    /// run the time-evolution
    void run(){
//        size_t n_layers_j = 0;
//        if ((!p_pars->run_jet_bws)&&(p_pars->is_jBW_init)){
//            p_pars->is_jBW_init = false;
//            p_bws_jet.clear();
//        }
//        if(p_pars->is_jBW_init) {
//            n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        }
//
//        size_t n_layers_ej = 0;
//        if ((!p_pars->run_ejecta_bws)&&(p_pars->is_ejBW_init)){
//            p_pars->is_ejBW_init = false;
////            p_bws_ej.clear();
//            p_ej.clear();
//        }
//        if(p_pars->is_ejBW_init){
//            n_layers_ej = ejectaStructs.structs[0].nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
//        }
//
//        if (!p_pars->is_main_pars_set){
//            (*p_log)(LOG_ERR,AT) << " model parameters were not set. Cannot run model. \n";
////            std::cerr << AT << "\n";
//            exit(1);
//        }
//        if ((!p_pars->is_mag_pars_set)&&(p_pars->run_magnetar)){
//            (*p_log)(LOG_ERR,AT)  << "magnetar pars are not set\n";
////            std::cerr << AT << "\n";
//            exit(1);
//        }
//        if ((!p_pars->is_jBW_init)&&(p_pars->run_jet_bws)){
//            (*p_log)(LOG_ERR,AT)  << "jet BWs are not set \n";
////            std::cerr << AT << "\n";
//            exit(1);
//        }
//        if ((!p_pars->is_ejBW_init)&&(p_pars->run_ejecta_bws)){
//            std::cerr << " ejecta BWs are not set\n Exiting...\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
//        if ((!p_pars->is_mag_pars_set)&&(!p_pars->is_jBW_init)&&(!p_pars->is_ejBW_init)){
//            (*p_log)(LOG_ERR,AT) << " magnetar AND jet AND ejecta BWs are not set\n";
////            std::cerr << AT << "\n";
//            exit(1);
//        }

        p_model = std::make_unique<EvolveODEsystem>( p_mag, p_grb, p_ej, t_grid );
//                p_magnetar, p_bws_jet, p_ej,
//                p_pars->run_magnetar, p_pars->run_jet_bws, p_pars->run_ejecta_bws,
//                t_grid, 0, ejectaStructs.nshells,
//                n_layers_j, n_layers_ej, p_pars->integrator, p_pars->loglevel);
        p_model->pIntegrator()->pPars()->rtol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->atol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->nmax = p_pars->nmax;
        p_model->setInitialConditions(p_pars->tb0);
        (*p_log)(LOG_INFO, AT) << "all initial conditions set. Starting evolution\n";
        // evolve
        double dt;
        for (size_t it = 1; it < t_grid.size(); it++){
            p_model->getPars()->i_restarts = 0;
            dt = t_grid[it] - t_grid[it - 1];
            p_model->evolve(dt, it);
        }
        (*p_log)(LOG_INFO, AT) << "evolution is completed\n";
    }

    /// save magnetar evolution
    void saveMagnetarEvolution(std::string fpath, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Saving magnetar evolution...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"\n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & magnetar = p_mag;

        std::vector<std::vector<double>>  tot_mag_out{};
        std::vector<std::string> tot_names {};
        std::unordered_map<std::string,double> group_attrs{};

        auto & mag_v_ns = magnetar->m_vnames;
        auto t_arr = magnetar->getTbGrid(every_it);

        tot_mag_out.resize( mag_v_ns.size() );
        for (size_t ivar = 0; ivar < mag_v_ns.size(); ivar++){
            for (size_t it = 0; it < magnetar->getTbGrid().size(); it = it + every_it)
                tot_mag_out[ivar].emplace_back( (*magnetar)[ static_cast<Magnetar::Q>(ivar) ][it] );
        }

        std::unordered_map<std::string,double> attrs {};

        p_out->VectorOfVectorsH5(tot_mag_out, mag_v_ns, fpath, attrs);
    }

    /// save dynamical evolution of the jet blast-waves
    void saveJetBWsDynamics(std::string fpath, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Saving jet BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"\n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & models = p_grb->getBWs();
        std::vector<
                std::vector<
                        std::vector<double>>> tot_dyn_out (models.size() );
        std::vector<std::string> tot_names {};
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        auto & dyn_v_ns = models[0]->m_vnames;
        auto t_arr = models[0]->getTbGrid(every_it);
        for (size_t i = 0; i < models.size(); i++) {
            tot_names.push_back("layer="+std::to_string(i));
            tot_dyn_out[i].resize( dyn_v_ns.size() );
            for (size_t ivar = 0; ivar < dyn_v_ns.size(); ivar++){
                for (size_t it = 0; it < models[i]->getTbGrid().size(); it = it + every_it)
                    tot_dyn_out[i][ivar].emplace_back( (*models[i])[ static_cast<BlastWaveBase::Q>(ivar) ][it] );
            }
            ///write attributes
            auto & model = models[i];
            std::unordered_map<std::string,double> group_attr{
                    {"Gamma0",model->getPars()->Gamma0},
                    {"M0",model->getPars()->M0},
                    {"R0",model->getPars()->R0},
                    {"theta0",model->getPars()->theta_b0},
                    {"theta_max",model->getPars()->theta_max},
                    {"tb0",model->getPars()->tb0},
                    {"ijl",model->getPars()->ijl},
                    {"ncells",model->getPars()->ncells},
                    {"ilayer",model->getPars()->ilayer},
                    {"ishell",model->getPars()->ishell},
                    {"ctheta0",model->getPars()->ctheta0},
                    {"E0",model->getPars()->E0},
                    {"theta_c_l",model->getPars()->theta_c_l},
                    {"theta_c_h",model->getPars()->theta_c_h},
                    {"eps_rad",model->getPars()->eps_rad},
                    {"entry_time",model->getPars()->entry_time},
                    {"entry_r",model->getPars()->entry_r},
                    {"first_entry_r",model->getPars()->first_entry_r},
                    {"min_beta_terminate",model->getPars()->min_beta_terminate},
                    {"every_it",every_it}
            };
            group_attrs.emplace_back( group_attr );
        }
//    VecVector other_data { latStruct.cthetas0 };
//    std::vector<std::string> other_names { "cthetas0" };

        std::unordered_map<std::string, double> attrs{
                {"nlayers", models.size()}
        };

        p_out->VecVectorOfVectorsAsGroupsH5(tot_dyn_out, tot_names, dyn_v_ns,
                                            fpath, attrs, group_attrs);
    }

    /// save dynamical evolution of the ejecta blast-waves
    void saveEjectaBWsDynamics(std::string fpath, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        auto & models = p_ej->getShells();

        std::vector<std::string> table_names;
        std::vector<std::vector<std::vector<double>>> tot_dyn_out ( nshells * nlayers );
        size_t i = 0;
        VecVector other_data;
        std::vector<std::string> other_names;
        auto & arr_names = models[0]->getBW(0)->m_vnames;//models[0][0].getBWdynPtr()->m_vnames;
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        for (size_t ishell = 0; ishell < nshells; ishell++){
            for(size_t ilayer = 0; ilayer < nlayers; ilayer++){
                table_names.push_back("shell="+std::to_string(ishell)+" layer="+std::to_string(ilayer));
                tot_dyn_out[i].resize( arr_names.size() );
                auto & bw = models[i]->getBW(ishell);

                std::unordered_map<std::string,double> group_attr{
                        {"Gamma0",bw->getPars()->Gamma0},
                        {"M0",bw->getPars()->M0},
                        {"R0",bw->getPars()->R0},
                        {"theta0",bw->getPars()->theta_b0},
                        {"theta_max",bw->getPars()->theta_max},
                        {"tb0",bw->getPars()->tb0},
                        {"ijl",bw->getPars()->ijl},
                        {"ncells",bw->getPars()->ncells},
                        {"ilayer",bw->getPars()->ilayer},
                        {"ishell",bw->getPars()->ishell},
                        {"ctheta0",bw->getPars()->ctheta0},
                        {"E0",bw->getPars()->E0},
                        {"theta_c_l",bw->getPars()->theta_c_l},
                        {"theta_c_h",bw->getPars()->theta_c_h},
                        {"eps_rad",bw->getPars()->eps_rad},
                        {"entry_time",bw->getPars()->entry_time},
                        {"entry_r",bw->getPars()->entry_r},
                        {"first_entry_r",bw->getPars()->first_entry_r},
                        {"min_beta_terminate",bw->getPars()->min_beta_terminate}
                };
                group_attrs.emplace_back( group_attr );

                for (size_t ivar = 0; ivar < arr_names.size(); ivar++) {
                    for (size_t it = 0; it < bw->getTbGrid().size(); it = it + every_it)
                        tot_dyn_out[i][ivar].emplace_back( (*bw)[ static_cast<BlastWaveBase::Q>(ivar) ][it] );
                }
                i++;
            }
        }
        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells },
                {"nlayers", nlayers }
        };
        p_out->VecVectorOfVectorsAsGroupsH5(tot_dyn_out, table_names, arr_names,
                                            fpath, attrs, group_attrs);
//    OutputVecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(fpath,tot_dyn_out, table_names,
//                                                         arr_names,other_data,other_names,attrs);
    }

    /// list parameters for analytic synchrotron & observer
//    static auto listAnalyticSynchrotronPars(){ return SynchrotronAnalytic::listPars(); }
//    static auto listAnalyticSynchrotronOpts(){ return SynchrotronAnalytic::listOpts(); }
//    static auto listObsPars(){ return RadBlastWave::listPars(); }
//    static auto listObsOpts(){ return RadBlastWave::listOpts(); }
    /// set parameters for analytical electron/synchrotron calculations
    void setPreComputeJetAnalyticElectronsPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Jet analytic electron pars...\n";

        if (!p_grb->run_jet_bws){
            std::cerr << " jet BWs were not evolved. Cannot compute electrons (analytic) \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_grb->getBWs();
        for (auto & model : models) {
//            model->setEatsPars(pars, opts);
//            model->getSynchAnPtr()->setPars( pars, opts );
            model->computeForwardShockElectronAnalyticVars();
            model->computeForwardShockSynchrotronAnalyticSpectrum();
        }
        p_grb->is_jet_anal_synch_computed = true;
    }
    void setPreComputeEjectaAnalyticElectronsPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Ejecta analytic electron pars...\n";

        if (!p_ej->run_ej_bws){
            std::cerr << " ejecta BWs were not evolved. Cannot compute electrons (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_ej->getShells();
        for (auto & model : models) {
            for (auto & bw : model->getBWs()) {
//            model->setEatsPars( pars, opts );
//            model->getSynchAnPtr()->setPars( pars, opts );
                bw->computeForwardShockElectronAnalyticVars();
                bw->computeForwardShockSynchrotronAnalyticSpectrum();
            }
        }
        p_ej->is_ejecta_anal_synch_computed = true;
    }

    /// compute and save comoving spectrum for all layers of the jet blast wave
    void computeSaveJetAnalyticSynchrotronSpectrum(std::string fpath, Vector & freq_array, size_t every_it){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet analytic synchrotron spectrum...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr<< " jet analytic electrons were not evolved. Cannot compute spectrum (analytic) exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_grb->getBWs();

        std::vector< // layers
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out{};
        auto t_arr = models[0]->getTbGrid(every_it);
        auto in_group_names = models[0]->getSynchAnPtr()->m_names_;
        std::vector<std::string> group_names{};
        for (size_t il = 0; il < models.size(); ++il ){
            group_names.emplace_back( "layer=" + std::to_string(il) );
            models[il]->computeForwardShockElectronAnalyticVars();
            out.emplace_back(models[il]->evalForwardShockComovingSynchrotron(freq_array, every_it) );
        }
        VecVector other_data{
                arrToVec( t_arr ),
                freq_array
        };
        std::vector<std::string> other_names { "times", "freqs" };
        std::unordered_map<std::string,double> attrs{
                {"eps_e", models[0]->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", models[0]->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", models[0]->getSynchAnPtr()->getPars()->eps_t },
                {"p", models[0]->getSynchAnPtr()->getPars()->p },
                {"ksi_n", models[0]->getSynchAnPtr()->getPars()->ksi_n }
        };
//    auto out = { emissivity, absorption, em_th, em_pl, abs_th, abs_pl };
        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);
    }
    void computeSaveEjectaAnalyticSynchrotronSpectrum(std::string fpath, Vector & freq_array, size_t every_it){

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta analytic synchrotron spectrum...\n";

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot compute spectrum (analytic) \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = p_ej->getShells();

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        std::vector< // layers
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out {};
        auto t_arr = models[0]->getBW(0)->getTbGrid(every_it);
        auto in_group_names = models[0]->getBW(0)->getSynchAnPtr()->m_names_;
        std::vector<std::string> group_names{};
        size_t ish_il = 0;
        for (size_t ishell = 0; ishell < nshells; ishell++) {
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                out.emplace_back(
                        models[ilayer]->getBW(ishell)->evalForwardShockComovingSynchrotron(freq_array, every_it) );
                ish_il++;
            }
        }
        VecVector other_data{
                arrToVec( t_arr ),
                freq_array
        };
        std::vector<std::string> other_names { "times", "freqs" };
        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells },
                {"nlayers", nlayers },
                {"eps_e", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", models[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
                {"p", models[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
                {"ksi_n", models[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
        };
//    auto out = { emissivity, absorption, em_th, em_pl, abs_th, abs_pl };
        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);
    }

    void updateJetObsPars(StrDbMap pars) {
        std::cout << "Updating Jet observer pars...\n";
        size_t n_layers_j = p_grb->getBWs().size();//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        size_t n_layers_ej = (p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
        for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++) {
            auto &model = p_grb->getBWs()[ilayer];
            model->updateObsPars(pars);
        }
//        p_pars->is_jet_obs_pars_set = true;
    }
    void updateEjectaObsPars(StrDbMap pars) {

        auto & models = p_ej->getShells();

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        std::cout << "Updating Ejecta observer pars...\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells; ishell++) {
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = struc.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                auto & model = p_ej->getShells()[ilayer]->getBW(ishell);//ejectaModels[ishell][ilayer];
                model->updateObsPars(pars);
                ii++;
            }
        }
//        p_pars->is_ejecta_obs_pars_set = true;
    }

    void computeSaveJetSkyImagesAnalytic(std::string fpath, Vector times, Vector freqs){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet sky image with analytic synchrotron...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr << "jet analytic electrons were not evolved. Cannot compute images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_grb->is_jet_obs_pars_set){
            std::cerr << "jet observer parameters are not set. Cannot compute image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells = 1; // jet has one shell only
        Image dummy(1);

        /// allocate memory for output data structure
        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};
        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(dummy.m_names.size());
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); ++i_vn) {
                    out[ii][i_vn].resize(nshells);
                }
                ii++;
            }
        }
        VecVector other_data{ times, freqs };
        std::vector<std::string> other_names { "times", "freqs" };

        /// allocate memory for image calculation
        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;
        Image image ( 2 * ncells );
        std::vector<Image> images ( nlayers );
        for (auto & im : images)
            im.resize( ncells );
        Image im_pj ( ncells );
        Image im_cj ( ncells );

        /// compute image for each time/freq
        ii = 0;
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            for (size_t it = 0; it < times.size(); ++it){
                im_pj.clearData();
                im_cj.clearData();
                for (auto & _image : images)
                    _image.clearData();

                (*p_log)(LOG_INFO,AT)<<" GRB SKYMAP "
                    << " freq="<<freqs[ifreq] << " (" << ifreq << "/" << freqs.size()-1 << ") "
                    << " time="<<times[it] << " (" << it << "/" << times.size()-1 << ") "
                    << " theta_layers="<<nlayers
                    << " ncells="<<ncells
                    << "\n";


                size_t n_jet_empty_images = 0;
                std::vector<size_t> n_empty_images;
                size_t n_layers_j = nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
                for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++){
                    auto & model = p_grb->getBWs()[ilayer];
                    model->evalImagePW(images[ilayer], im_pj, im_cj, times[it], freqs[ifreq]);
                    if (images[ilayer].m_f_tot == 0){
                        n_jet_empty_images += 1;
                        n_empty_images.emplace_back(ilayer);
                    }
                }
                if (n_jet_empty_images > 0){
                    (*p_log)(LOG_INFO,AT)
                        << " grb empty images n="<<n_jet_empty_images<<"/"<<n_layers_j<<" layers.\t"
                              << " Layers with empty image: " << n_jet_empty_images << "\n";
                }

                /// combine all images from all layers into one image
                combineImages(image, ncells, nlayers, images);

                /// copy data from images into output structure
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells; ish++) {
                        out[ii][i_vn][ish] = image.m_data[i_vn];//arrToVec(image.m_data[i_vn]);
                    }
                }

                /// separately output total fluxes
                for (size_t ish = 0; ish < nshells; ish++) {
                    tota_flux[it] += image.m_f_tot;
                }
                ii++;
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );
        }

        std::vector<std::string> group_names{};
//        for (size_t it = 0; it < times.size(); it++){
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
//                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
//            }
//        }
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            for (size_t it = 0; it < times.size(); it++){
                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
            }
        }

        auto in_group_names = dummy.m_names;

        std::unordered_map<std::string,double> attrs{
//                {"nlayers", lay },
                {"E0", p_grb->getBWs()[0]->getPars()->E0 },
                {"Gamma0", p_grb->getBWs()[0]->getPars()->Gamma0 },
                {"M0", p_grb->getBWs()[0]->getPars()->M0 },
                {"R0", p_grb->getBWs()[0]->getPars()->R0 },
                {"ctheta0", p_grb->getBWs()[0]->getPars()->ctheta0 },
                {"theta_b0", p_grb->getBWs()[0]->getPars()->theta_b0 },
                {"theta_w", p_grb->getBWs()[0]->getPars()->theta_w },
                {"thetaObs", p_grb->getBWs()[0]->getEatsPars()->theta_obs },
                {"d_L", p_grb->getBWs()[0]->getEatsPars()->d_l },
                {"z", p_grb->getBWs()[0]->getEatsPars()->z },
                {"eps_e", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_t },
                {"p", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->p },
                {"ksi_n", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->ksi_n }
        };

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    void computeSaveEjectaSkyImagesAnalytic(std::string fpath, Vector times, Vector freqs){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta sky image with analytic synchrotron...\n";

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr  << "ejecta analytic electrons were not evolved. Cannot compute images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_ej->is_ejecta_obs_pars_set){
            std::cerr<< "ejecta observer parameters are not set. Cannot compute image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

//        auto & ejectaStruct = ejectaStructs;
//        size_t nshells = ejectaStruct.nshells;
        Image dummy(1,0);

        std::vector< // times & freqs
                std::vector< // v_ns
                        std::vector< // shells
                                std::vector<double>>>> // data
        out {};

        size_t ii = 0;
        out.resize(times.size() * freqs.size());
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            for (size_t it = 0; it < times.size(); ++it){
                out[ii].resize(dummy.m_names.size());
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); ++i_vn) {
                    out[ii][i_vn].resize(nshells);
                }
                ii++;
            }
        }

        VecVector other_data{
                times,
                freqs
        };
        ii = 0;
        std::vector<std::string> other_names { "times", "freqs" };
        for (size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            VecVector total_flux_shell( nshells );
            for (auto & total_flux_shel : total_flux_shell)
                total_flux_shel.resize( times.size(), 0.0 );
            for (size_t it = 0; it < times.size(); ++it){
//                auto images = computeEjectaSkyMapPieceWise( times[it],freqs[ifreq]);
                std::vector<Image> images(nshells);
                computeEjectaSkyMapPieceWise( images, times[it],freqs[ifreq]);
                for (size_t i_vn = 0; i_vn < dummy.m_names.size(); i_vn++) {
                    for (size_t ish = 0; ish < nshells; ish++) {
                        out[ii][i_vn][ish] = images[ish].m_data[i_vn];//arrToVec(images[ish].m_data[i_vn]);
                    }
                }
                for (size_t ish = 0; ish < nshells; ish++) {
                    tota_flux[it] += images[ish].m_f_tot;
                    total_flux_shell[ish][it] = images[ish].m_f_tot;
                }
                ii++;
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );

            for (size_t ish = 0; ish < nshells; ish++){
                other_data.emplace_back( total_flux_shell[ish] );
                other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) +
                                          " shell=" + string_format("%d", ish));
            }

        }

        std::vector<std::string> group_names{};
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
            for (size_t it = 0; it < times.size(); it++){
                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
                                         + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
            }
        }


//        for (size_t it = 0; it < times.size(); it++){
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                group_names.emplace_back("time=" +  string_format("%.4e",times[it])  //std::to_string(times[it])
//                                       + " freq=" + string_format("%.4e",freqs[ifreq])); //std::to_string(freqs[ifreq]));
//            }
//        }

        auto in_group_names = dummy.m_names;

        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells },
                {"thetaObs", p_ej->getShells()[0]->getBW(0)->getEatsPars()->theta_obs },
                {"d_L", p_ej->getShells()[0]->getBW(0)->getEatsPars()->d_l },
                {"z",  p_ej->getShells()[0]->getBW(0)->getEatsPars()->z },
                {"eps_e",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
                {"p",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
                {"ksi_n",  p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
        };

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    /*   */
    void computeSaveJetLightCurveAnalytic(std::string fpath, Vector times, Vector freqs){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Jet light curve with analytic synchrotron...\n";

        if (!p_grb->is_jet_anal_synch_computed){
            std::cerr << "jet analytic electrons were not evolved. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_grb->is_jet_obs_pars_set){
            std::cerr << " jet observer parameters are not set. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;

//        auto & structure_hist = jetStruct;
        auto & tmp = p_grb->getBWs()[0]->getSynchAnPtr();
//        size_t nshells =structure_hist.nshells;
//        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;

        std::vector< // layers
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out {};

        Vector _times ( freqs.size() * times.size(), 0.0 );
        Vector _freqs ( freqs.size() * times.size(), 0.0 );
        size_t ii = 0;
        for (double freq : freqs){
            for (double time : times){
                _times[ii] = time;
                _freqs[ii] = freq;
                ii ++ ;
            }
        }

        auto light_curve = evalJetLightCurves( _times,_freqs);//[ilayer][it+ifreq]

        std::vector<std::string> group_names;
        std::vector<std::string> in_group_names{  // single entry
                tmp->m_names_[tmp->getPars()->m_marg21opt_em] + " " +
                tmp->m_names_[tmp->getPars()->m_marg21opt_abs]
        };
        size_t i_v_n = 0;
        size_t ishil = 0;
        out.resize(nlayers);
        for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {

            out[ishil].resize(1);
            ii = 0;
            out[ishil][i_v_n].resize(freqs.size());
            for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
                out[ishil][i_v_n][ifreq].resize(times.size());
                for (size_t it = 0; it < times.size(); ++it){
                    out[ishil][i_v_n][ifreq][it] = light_curve[ilayer][ii];
                    ii++;
                }
            }
            ishil++;
            group_names.emplace_back("layer=" + std::to_string(ilayer));
        }


        VecVector other_data{
                times,
                freqs
        };
        std::vector<std::string> other_names { "times", "freqs" };

        for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            for (size_t ishil = 0; ishil < nlayers; ++ishil){
                for (size_t it = 0; it < times.size(); ++it){
                    tota_flux[it] += out[ishil][i_v_n][ifreq][it];
                }
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e",freqs[ifreq]) );
        }

        StrDbMap attrs{
                {"E0", p_grb->getBWs()[0]->getPars()->E0 },
                {"Gamma0", p_grb->getBWs()[0]->getPars()->Gamma0 },
                {"M0", p_grb->getBWs()[0]->getPars()->M0 },
                {"R0", p_grb->getBWs()[0]->getPars()->R0 },
                {"ctheta0", p_grb->getBWs()[0]->getPars()->ctheta0 },
                {"theta_b0", p_grb->getBWs()[0]->getPars()->theta_b0 },
                {"theta_w", p_grb->getBWs()[0]->getPars()->theta_w },
                {"thetaObs", p_grb->getBWs()[0]->getEatsPars()->theta_obs },
                {"nlayers", nlayers },
                {"thetaObs", p_grb->getBWs()[0]->getEatsPars()->theta_obs },
                {"d_L", p_grb->getBWs()[0]->getEatsPars()->d_l },
                {"z", p_grb->getBWs()[0]->getEatsPars()->z },
                {"eps_e", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->eps_t },
                {"p", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->p },
                {"ksi_n", p_grb->getBWs()[0]->getSynchAnPtr()->getPars()->ksi_n }
        };

        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);

    }

    void computeSaveEjectaLightCurveAnalytic(std::string fpath, Vector times, Vector freqs){
        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta light curve with analytic synchrotron...\n";

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

        if (!p_ej->is_ejecta_anal_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!p_ej->is_ejecta_obs_pars_set){
            std::cerr << " ejecta observer parameters are not set. Cannot compute light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

//        auto & structure_hist = ejectaStructs;
        auto & tmp = p_ej->getShells()[0]->getBW(0)->getSynchAnPtr();
//        size_t nshells =structure_hist.nshells;
//        size_t n_layers_j = (p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        size_t n_layers_ej = structure_hist.structs[0].nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? structure_hist.structs[0].nlayers_pw : structure_hist.structs[0].nlayers_a ;
//        size_t nlayers = structure_hist.structs[0].nlayers_pw;

//        std::vector<std::string> group_names{};
//        std::vector< // freq
//                std::vector< // layers+shells
//                        std::vector<double>>> // times
//        out{};

        std::vector< // layers / shells
                std::vector< // options
                        std::vector< // freqs
                                std::vector<double>>>> // times
        out {};

        Vector _times ( freqs.size() * times.size(), 0.0 );
        Vector _freqs ( freqs.size() * times.size(), 0.0 );
        size_t ii = 0;
        for (double freq : freqs){
            for (double time : times){
                _times[ii] = time;
                _freqs[ii] = freq;
                ii ++ ;
            }
        }



        auto light_curve = evalEjectaLightCurves( _times,_freqs);//[ishll][ilayer][it+ifreq]

        std::vector<std::string> group_names;
        std::vector<std::string> in_group_names{  // single entry
                tmp->m_names_[tmp->getPars()->m_marg21opt_em] + " " +
                tmp->m_names_[tmp->getPars()->m_marg21opt_abs]
        };
        size_t i_v_n = 0;
        size_t ishil = 0;
        out.resize(nshells*nlayers);
        for (size_t ishell = 0; ishell < nshells; ++ishell) {
            for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
//                (*p_log)(LOG_INFO,AT)
//                         <<" EJECTA LC shell="<<ishell<<"/"<<nshells<<" layer="<<ilayer<<"/"<<n_layers_ej<<"\n";
                out[ishil].resize(1);
                ii = 0;
                out[ishil][i_v_n].resize(freqs.size());
                for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
                    out[ishil][i_v_n][ifreq].resize(times.size());
                    for (size_t it = 0; it < times.size(); ++it){
                        out[ishil][i_v_n][ifreq][it] = light_curve[ishell][ilayer][ii];
                        ii++;
                    }
                }
                ishil++;
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
            }
        }

        VecVector other_data{
                times,
                freqs
        };
        std::vector<std::string> other_names { "times", "freqs" };

        for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
            Vector tota_flux(times.size(), 0.0);
            for (size_t ishil = 0; ishil < nshells * nlayers; ++ishil){
                for (size_t it = 0; it < times.size(); ++it){
                    tota_flux[it] += out[ishil][i_v_n][ifreq][it];
                }
            }
            other_data.emplace_back( tota_flux );
            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e",freqs[ifreq]) );
        }

        std::unordered_map<std::string,double> attrs{
                {"nshells", nshells },
                {"nlayers", nlayers },
                {"thetaObs", p_ej->getShells()[0]->getBW(0)->getEatsPars()->theta_obs },
                {"d_L", p_ej->getShells()[0]->getBW(0)->getEatsPars()->d_l },
                {"z", p_ej->getShells()[0]->getBW(0)->getEatsPars()->z },
                {"eps_e", p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_e },
                {"eps_b", p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_b },
                {"eps_t", p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->eps_t },
                {"p", p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->p },
                {"ksi_n", p_ej->getShells()[0]->getBW(0)->getSynchAnPtr()->getPars()->ksi_n }
        };


        p_out->VectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath,out,group_names,
                                                          in_group_names,
                                                          other_data,other_names,attrs);


//        out.resize(freqs.size());
//        size_t ishil = 0;
//        for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
//            out[ifreq].resize(nshells*nlayers);
//            ishil = 0;
//            for (size_t ishell = 0; ishell < nshells; ++ishell) {
//                for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
//                    out[ifreq][ishil].resize( times.size() );
//                }
//                ishil ++;
//            }
//            ii++;
//        }
//
//        ishil = 0;
//        for (size_t ishell = 0; ishell < nshells; ++ishell) {
//            for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
//                ii = 0;
//                for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq){
//                    for (size_t it = 0; it < times.size(); ++it){
//                        out[ifreq][ishil][it] = light_curve[ishell][ilayer][ii];
//                        ii ++ ;
//                    }
//                }
//                ishil++;
//            }
//        }




//
////        std::cout << _times << "\n";
////        std::cout << _freqs << "\n";
//
//
//        {
////            out.resize(1);
////            group_names.emplace_back(
////                    tmp->m_names_[tmp->getPars()->m_marg21opt_em] + " " +
////                    tmp->m_names_[tmp->getPars()->m_marg21opt_abs]
////            );
//            /// compute light curves
//            auto light_curve = evalEjectaLightCurves( _times,_freqs);
//            size_t ish_il = 0;
//            ii = 0;
//            /// place light curves into storage
//            out.resize( nshells * nlayers );
//            for (size_t ishell = 0; ishell < nshells; ++ishell) {
//                for (size_t ilayer = 0; ilayer < nlayers; ++ilayer) {
//
//                    group_names.emplace_back("shell="+std::to_string(ishell)+" layer="+std::to_string(ilayer));
//
//                    ii = 0;
//                    out[ish_il].resize(freqs.size() );
//                    for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq) {
//                        out[ish_il][ifreq].resize(times.size());
//                        for (size_t it = 0; it < times.size(); ++it){
//                            out[ish_il][ifreq][it] = light_curve[ishell][ilayer][ii];//.resize( times.size() );
//                            ii++;
//                        }
//                        //out[o].emplace_back(light_curve[ishell][ilayer][ii]);
//                    }
//                    ish_il++;
//                }
//            }
//            VecVector tot_fluxes(freqs.size());
//            for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq) {
//                tot_fluxes[ifreq].resize(times.size(), 0.);
//                for (size_t it = 0; it < times.size(); ++it){
//                    for (size_t ish_il = 0; ish_il < nshells * nlayers; ++ish_il) {
//                        tot_fluxes[ifreq][it] += out[ish_il][ifreq][it];
//                    }
//                }
//            }
//
////            Vector tot_flux_ejetca (times.size(), 0);
////            for (size_t k = 0; k < light_curve.size(); k++) {
////                for (size_t l = 0; l < light_curve[0].size(); l++) { // layer
////                    for (size_t m = 0; m < light_curve[0][0].size(); m++) { // time
////                        tot_flux_ejetca[m] += light_curve[k][l][m];
////                    }
////                }
////            }
//            out.emplace_back(tot_fluxes);
//        }
//
//
//        /// names for arrays within a group
////        std::vector<std::string> out_names_ej{};
////        for (size_t ishell = 0; ishell < nshells; ishell++) {
////            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
////                out_names_ej.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
////            }
////        }
////        out_names_ej.emplace_back("totalflux");
//        std::vector<std::string> out_names_ej{};
//        for(size_t ifreq = 0; ifreq < freqs.size(); ++ifreq) {
//            out_names_ej.emplace_back("freq="+std::to_string(freqs[ifreq]));
//        }
//
//        std::unordered_map<std::string,double> attrs{
//                {"nshells", nshells },
//                {"nlayers", nlayers },
//                {"thetaObs", p_bws_ej[0]->getEatsPars()->theta_obs },
//                {"d_L", p_bws_ej[0]->getEatsPars()->d_l },
//                {"z", p_bws_ej[0]->getEatsPars()->z },
//                {"eps_e", p_bws_ej[0]->getSynchAnPtr()->getPars()->eps_e },
//                {"eps_b", p_bws_ej[0]->getSynchAnPtr()->getPars()->eps_b },
//                {"eps_t", p_bws_ej[0]->getSynchAnPtr()->getPars()->eps_t },
//                {"p", p_bws_ej[0]->getSynchAnPtr()->getPars()->p },
//                {"ksi_n", p_bws_ej[0]->getSynchAnPtr()->getPars()->ksi_n }
//        };
//
//        VecVector common_data{
//                times,
//                freqs
//        };
//        std::vector<std::string> common_names{
//                "times",
//                "freqs"
//        };
//
//        OutputVecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(fpath ,
//                                                             out,
//                                                             group_names,
//                                                             out_names_ej,
//                                                             common_data,
//                                                             common_names,
//                                                             attrs);
//
////        OutputVectorOfTablesAsGroupsAndVectorOfVectorsH5(fpath, out, group_names,
////                                           fpath, attrs, group_attrs);
////
////        OutputVecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(
////                fpath ,
////                out,
////                group_names,
////                out_names_ej,
////                common_data,
////                common_names,
////                attrs
////        );
    }

private:
//    void computeJetSkyMapPieceWise(Image & image, std::vector<Image> & images, Image & im_pj, Image & im_cj,
//                                   double obs_time, double obs_freq){
////        std::vector<Image> images (jetStruct.nlayers);
////        for (auto & im : images)
////            im.resize( jetStruct.ncells );
//////        std::vector<Image> im_pjs (jetStruct.nlayers);
//////        std::vector<Image> im_cjs (jetStruct.nlayers);
////        Image im_pj ( jetStruct.ncells );
////        Image im_cj ( jetStruct.ncells );
//        size_t n_jet_empty_images = 0;
//        std::vector<size_t> n_empty_images;
//        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
//        for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++){
//
//            auto & model = p_bws_jet[ilayer];
////            model->setEatsPars(pars, opts);
////            images.emplace_back( model->evalImagePW(obs_time, obs_freq) );
//            model->evalImagePW(images[ilayer], im_pj, im_cj, obs_time, obs_freq);
//            if (images[ilayer].m_f_tot == 0){
//                n_jet_empty_images += 1;
//                n_empty_images.emplace_back(ilayer);
//            }
//        }
//        if (n_jet_empty_images > 0){
//            std::cout << "Jet Empty images n="<<n_jet_empty_images<<"/"<<n_layers_j<<" layers.\t"
//                      << " Layers with empty image: " << n_jet_empty_images << "\n";
//        }
//        combineImages(image, jetStruct, images);
////        return std::move( combineImages(jetStruct, images) );
//    }


    void computeEjectaSkyMapPieceWise( std::vector<Image> & images, double obs_time, double obs_freq ){

        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();

//        std::vector<Image> images;
//        images.resize(ejectaStructs.nshells);
        if (images.empty()){
            (*p_log)(LOG_ERR,AT) << " empty image passed. Exiting...\n";
            exit(1);
        }
        if (images.size() != nshells){
            (*p_log)(LOG_ERR,AT) << " number of images does not equal to the number of shells. Exiting...\n";
            exit(1);
        }
        // ejectaStructs.nshells);
        size_t ii = 0;
        size_t n_jet_empty_images = 0;
//        std::vector<size_t> n_empty_images_layer;
        std::vector<std::vector<size_t>> n_empty_images;
        std::vector<size_t> n_empty_images_shells;
        std::vector<Image> tmp (nlayers);
        for (auto & _tmp : tmp)
            _tmp.resize( ncells );
        Image tmp_pj( ncells); // std::vector<Image> tmp_pj (ejectaStructs.structs[0].nlayers);
        Image tmp_cj( ncells); // std::vector<Image> tmp_cj (ejectaStructs.structs[0].nlayers);
        for (size_t ishell = 0; ishell <  nshells; ishell++){
            for (auto & _tmp : tmp)
                _tmp.clearData();
//            tmp[].clearData();
            tmp_pj.clearData(); tmp_cj.clearData();
//            image_i.clearData();
//            im_pj.clearData();
//            im_cj.clearData();
//            auto & ejectaStruct = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStruct.nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
//            std::vector<Image> tmp;
//            std::vector<Image> tmp (n_layers_ej);
            std::vector<size_t> n_empty_images_layer;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
                auto & model = p_ej->getShells()[ ilayer ];
//                tmp.emplace_back( model->evalImagePW(obs_time, obs_freq) );
                model->getBW(ishell)->evalImagePW(tmp[ilayer], tmp_pj, tmp_cj, obs_time, obs_freq);
                if (tmp[ilayer].m_f_tot == 0){
                    n_jet_empty_images += 1;
                    n_empty_images_layer.emplace_back(ilayer);
                }
                ii ++ ;
            }
            if(!n_empty_images_layer.empty()){
                n_empty_images_shells.emplace_back(ishell);
                n_empty_images.emplace_back(n_empty_images_layer);
            }

//            images.emplace_back( combineImages(ejectaStruct,tmp) );
            combineImages(images[ishell], ncells, nlayers, tmp) ;
        }
        if (n_jet_empty_images > 0){
            auto & ccerr =std::cout;
            ccerr << "Ejecta at tobs=" << obs_time << " freq=" << obs_freq << " gave an empty images for total n="
                  <<n_jet_empty_images<<" layers. Specifically:\n";
            for (size_t ish = 0; ish < n_empty_images_shells.size(); ish++){
//                auto & ejectaStruct = ejectaStructs.structs[n_empty_images_shells[ish]];
//                size_t n_layers_ej = nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStruct.nlayers_pw : ejectaStruct.nlayers_a ;
                ccerr << "\t [ishell="<<n_empty_images_shells[ish] << " ilayer] = [ ";
                for (size_t il = 0; il < n_empty_images[ish].size(); il++){
                    ccerr << n_empty_images[ish][il] << " ";
                }
                ccerr << "] / (" << nlayers << " total layers) \n";
            }
        }

//        return std::move( images );
    }
#if 0
    VecVector evalJetLightCurves_NEW( Vector & obs_times, Vector & obs_freqs ){
        double stop_frac = p_pars->jet_layer_fnu_stop_frac; /// WARNING! Requires monotonic decrease in Flux with layer angle!
        double rtol = p_pars->jet_rtol;
//        std::vector<size_t> sorted_layers ( jetStruct.nlayers, 0 );
//        Vector theta_c_minus_theta_obs ( jetStruct.nlayers, 0. );
//        for (size_t il = 0; il < jetStruct.nlayers; il++){
//            auto & layer = p_bws_jet[il];
//            theta_c_minus_theta_obs[il] = mu_max()
//        }

        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
        VecVector light_curves(n_layers_j); // [i_layer][i_time]
        for (auto & arr : light_curves)
            arr.resize(obs_times.size(), 0.);
        /// evaluate flux density for each time/frequency
        Image image; Image im_pj; Image im_cj;
        for (size_t it = 0; it < obs_times.size(); it++) {
            Vector layer_fluxes ( n_layers_j, 0.0 );
            double peak_flux_for_layers = 0.0;
            size_t last_layer = 0;
//            std::fill(layer_fluxes.begin(), layer_fluxes.end(), 0.0);
            /// evaluate flux density for each lateral layer
            for (size_t ilayer = 0; ilayer < n_layers_j; ilayer++){
                auto & model = p_bws_jet[ilayer];
                /// check if blast wave for this layer was evolved
                if ((model->getVal(BlastWaveBase::Q::iR,0) == 0.0)
                    && (model->getVal(BlastWaveBase::Q::iR,-1) == 0.0)){
                    (*p_log)(LOG_WARN, AT)
                            << " blast wave not evolved, flux=0 [ishell="
                            << model->getPars()->ishell<<", ilayer="<<model->getPars()->ilayer<<"]\n";
                }
                /// piece-wise
                if (p_pars->jet_eats_method == LatStruct::METHOD_eats::i_pw) {
                    model->evalImagePW(image, im_pj, im_cj, obs_times[it], obs_freqs[it]);
                    layer_fluxes[ilayer] = image.m_f_tot;
                    light_curves[ilayer][it] += image.m_f_tot;

                    /// try to spip more equatorial layers where there is not much flux for a collimated grb...
                    if ((ilayer > 0) && (peak_flux_for_layers == 0.0) &&
                        (layer_fluxes[ilayer] < layer_fluxes[ilayer - 1])) {
                        peak_flux_for_layers = layer_fluxes[ilayer];
                    }
                    if (layer_fluxes[ilayer] < (stop_frac * peak_flux_for_layers)) {
                        (*p_log)(LOG_INFO, AT)
//                              << model->getPars()->ishell<<", ilayer="<<model->getPars()->ilayer
                                << " t=" << obs_times[it]
                                << " freq="<< obs_freqs[it]
                                << " flux decrease after peak="<<peak_flux_for_layers
                                << " stopping at layer="<<model->getPars()->ilayer<<"\n";
//                                << " skipping other layers as Fnu < frac*Fnu_max; "
//                                << " Fnu=" << layer_fluxes[ilayer]
//                                << " Fnu_max=" << peak_flux_for_layers
//                                << " fraction set="<<stop_frac
//                                << " [ishell="<< model->getPars()->ishell<<", ilayer="<<model->getPars()->ilayer<<"]\n";
                        break;
                    }
                }
                /// adaptive
                if (p_pars->jet_eats_method == LatStruct::METHOD_eats::i_adap) {
                    double atol = light_curves[ilayer][it] * rtol / (double) model->getPars()->nlayers;
                    light_curves[ilayer][it] += model->evalFluxDensA(obs_times[it], obs_freqs[it], atol);
                }
                last_layer = ilayer;
            }
            (*p_log)(LOG_INFO, AT)
                << " t=" << obs_times[it]
                << " freq="<< obs_freqs[it]
                << " done"
                << " (it="<<it<<'/'<<obs_times.size()<<")"
                << " Last layer="<<last_layer<<'/'<<n_layers_j-1
                << " (peak flux="<<peak_flux_for_layers
                << " last layer flux ="<<layer_fluxes[last_layer]
                << ") "
                <<'\n';
        }
        return std::move( light_curves );
    }
#endif
    VecVector evalJetLightCurves( Vector & obs_times, Vector & obs_freqs ){
        (*p_log)(LOG_INFO,AT)<<" starting grb light curve calculation\n";
        size_t nlayers = p_grb->getBWs().size();
        auto ncells = (size_t)p_grb->getBWs()[0]->getPars()->ncells;
//        size_t n_layers_j = jetStruct.nlayers;//(p_pars->jet_method_eats == LatStruct::i_pw) ? jetStruct.nlayers_pw : jetStruct.nlayers_a ;
        VecVector light_curves(nlayers); // [i_layer][i_time]
        for (auto & arr : light_curves)
            arr.resize(obs_times.size(), 0.);
        double flux;
        double rtol = p_grb->jet_rtol;
//        Image image;
        Image image_i ( ncells ); Image im_pj ( ncells ); Image im_cj ( ncells );
        for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
            auto & model = p_grb->getBWs()[ilayer];
            (*p_log)(LOG_INFO,AT)<<" GRB LC ntimes="<<obs_times.size()<<" theta_layer="<<ilayer<<"/"<<nlayers<<
            " phi_cells="<<LatStruct::CellsInLayer(model->getPars()->ilayer)<<"\n";
            model->evalForwardShockLightCurve(image_i, im_pj, im_cj, light_curves[ilayer], obs_times, obs_freqs);
//            for (size_t it = 0; it < obs_times.size(); it++) {
//                if (model->getPars()->m_method_eats == LatStruct::i_pw) {
////                    auto image = model->evalImagePW(obs_times[it], obs_freqs[it]);
//                    model->evalImagePW(image, obs_times[it], obs_freqs[it]);
//                    light_curves[ilayer][it] += image.m_f_tot;
//                }
//                else{
//                    double atol = light_curves[ilayer][it] * rtol / (double)n_layers_j;
//                    light_curves[ilayer][it] += model->evalFluxDensA(obs_times[it], obs_freqs[it], atol);
//                }
//            }
        }
        return std::move( light_curves );
    }
    std::vector<VecVector> evalEjectaLightCurves( Vector & obs_times, Vector & obs_freqs){
        (*p_log)(LOG_INFO,AT)<<" starting ejecta light curve calculation\n";
        size_t nshells = p_ej->nshells();
        size_t nlayers = p_ej->nlayers();
        size_t ncells =  (int)p_ej->ncells();
        std::vector<VecVector> light_curves(nshells); // [ishell][i_layer][i_time]
        for (auto & arr : light_curves){
            size_t n_layers_ej = nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
            arr.resize(n_layers_ej);
            for (auto & arrr : arr){
                arrr.resize( obs_times.size(), 0. );
            }
        }
        double flux_pj, flux_cj; size_t ii = 0;
        Image image;
        double rtol = p_ej->ej_rtol;
        Image image_i ( ncells );
        Image im_pj ( ncells );
        Image im_cj ( ncells );
        for (size_t ishell = 0; ishell < nshells; ishell++){
            image_i.clearData();
            im_pj.clearData();
            im_cj.clearData();
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStructs.structs[ishell].nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers; ilayer++) {
                auto & model = p_ej->getShells()[ilayer];//ejectaModels[ishell][ilayer];
//                model->setEatsPars( pars, opts );
                (*p_log)(LOG_INFO,AT)
                    << " EJECTA LC ntimes="<<obs_times.size()
                    << " vel_shell="<<ishell<<"/"<<nshells-1
                    << " theta_layer="<<ilayer<<"/"<<nlayers
                    << " phi_cells="<<LatStruct::CellsInLayer(ilayer)<<"\n";
                model->getBW(ishell)->evalForwardShockLightCurve(image_i, im_pj, im_cj,
                                                                 light_curves[ishell][ilayer], obs_times, obs_freqs);

//                for (size_t it = 0; it < obs_times.size(); it++) {
//                    if (model->getPars()->m_method_eats == LatStruct::i_pw) {
//                        //                    auto tmp = model->evalImagePW(obs_times[it], obs_freqs[it] );
//                        model->evalImagePW(image, obs_times[it], obs_freqs[it] );
//                        light_curves[ishell][ilayer][it] += image.m_f_tot;
//                    }
//                    else{
//                        double atol = light_curves[ishell][ilayer][it] * rtol / (double)n_layers_ej;
//                        light_curves[ishell][ilayer][it] += model->evalFluxDensA(obs_times[it], obs_freqs[it], atol);
//                    }
//                }
                ii ++;
            }
        }
        return std::move( light_curves );
    }

    /// getters
//    auto & getBWjetPtrs(){ return p_bws_jet; }
//    auto & getBWejectaPtrs(){ return p_bws_ej; }
//    auto & getStructJet(){ return jetStruct; }
//    auto & getStructEjecta(){ return ejectaStructs; }
//    auto & getBWjetUptr( size_t ilayer ){ return p_bws_jet[ilayer]; }
//    auto & getBWejUptr( size_t ishell, size_t ilayer ){ //x + m_nX * y
//        return p_bws_ej[ishell + ejectaStructs.nshells * ilayer];
//    }
    /// main methods

private:
    /// List of the parameters for each blast wave
//    static std::vector<std::string> listBwPars_(){ return { "nn","A0","s","r_ej","r_ism",
//                                                            "a", "theta_max", "epsilon_e_rad",
//                                                            "which_jet_layer_to_use",
//                                                            "steepnes_of_exp_decay",
//                                                            "Gamma_when_st_starts",
//                                                            "fraction_of_Gamma0_when_bm_for_bm"
//        };
//    }
//    static std::vector<std::string> listBwOpts_(){ return { "method_spread","method_eos",
//                                                            "use_dens_prof_behind_jet_for_ejecta",
//                                                            "use_exp_rho_decay_as_floor",
//                                                            "use_flat_dens_floor",
//                                                            "use_st_dens_profile",
//                                                            "use_bm_dens_profile",
//                                                            "method_GammaSh",
//                                                            "method_Up",
//                                                            "method_Delta",
//                                                            "use_adiabLoss",
//                                                            "method_Rsh",
//                                                            "method_dmdr",
//                                                            "method_dgdr"
//        };
//    }
//    void setAllParametersForOneLayer(LatStruct & latStruct, RadBlastWave & bw_obj,
//                                     StrDbMap & pars, StrStrMap & opts,
//                                     size_t ilayer, size_t ii_eq){
//
//        if (!std::isfinite(latStruct.dist_E0_pw[ilayer])||
//            !std::isfinite(latStruct.dist_M0_pw[ilayer])||
//            !std::isfinite(latStruct.dist_G0_pw[ilayer])){
//            std::cerr << AT << "\n E0="<<latStruct.dist_E0_pw[ilayer]
//                      << " M0=" << latStruct.dist_M0_pw[ilayer]
//                      << " G0=" << latStruct.dist_G0_pw[ilayer]
//                      << " wrong value. exiting...\n";
////            std::cerr << AT << " error." << "\n";
//            exit(1);
//        }
//
//        double nism, A0, s, r_ej, r_ism,  a, theta_max, epsilon_e_rad;
//
//        // set parameters for ISM density
//        nism = getDoublePar("n_ism", pars, AT, p_log, -1, true);//pars.at("nism");
//        A0 = getDoublePar("A0", pars, AT,p_log,-1,false);//pars.at("A0");
//        s = getDoublePar("s", pars, AT,p_log,-1,false);//pars.at("s");
//        r_ej = getDoublePar("r_ej", pars, AT,p_log,-1,false);//pars.at("r_ej");
//        r_ism = getDoublePar("r_ism", pars, AT,p_log,-1,false);//pars.at("r_ism");
//        bw_obj.getDensIsm()->setPars(nism, A0, s, r_ej, r_ism, true);
//
//        // spreading
//        a = getDoublePar("a", pars, AT,p_log,-1,false);//pars.at("a");
//        theta_max = getDoublePar("theta_max", pars, AT,p_log,CGS::pi/2.,false);//pars.at("theta_max");
//
//        // radiative losses
//        epsilon_e_rad = getDoublePar("epsilon_e_rad", pars, AT,p_log,0.,false);// pars.at("epsilon_e_rad");
//
//        // interaction parameters
//        bw_obj.getPars()->which_jet_layer_to_use =
//                (int)getDoublePar("which_jet_layer_to_use",pars,AT,p_log,1.e5,false);
//        bw_obj.getPars()->steepnes_of_exp_decay =
//                (double)getDoublePar("steepnes_of_exp_decay",pars,AT,p_log,1.,false);
//        bw_obj.getPars()->Gamma_when_st_starts =
//                (double)getDoublePar("Gamma_when_st_starts",pars,AT,p_log,2.,false);
//        bw_obj.getPars()->fraction_of_Gamma0_when_bm_for_bm =
//                (double)getDoublePar("fraction_of_Gamma0_when_bm_for_bm",pars, AT,p_log,1.98,false);
//        bw_obj.getPars()->fraction_of_Gamma0_when_spread =
//                (double)getDoublePar("fraction_of_Gamma0_when_spread",pars, AT,p_log,.75,false);
//
//
//        // check if unknown option is given
////        check_for_unexpected_opt(opts, listBwOpts(), "listBwOpts()");
//
//        // set options
//        std::string opt;
//
//        // mass accretion from ISM
//        opt = "method_dmdr";
//        METHOD_dmdr method_dmdr;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_dmdr = METHOD_dmdr::iusingdthdR;
//        }
//        else{
//            if(opts.at(opt) == "usingA")
//                method_dmdr = METHOD_dmdr::iusingA;
//            else if(opts.at(opt) == "usingdthdr")
//                method_dmdr = METHOD_dmdr::iusingdthdR;
//            else{
//                std::cerr << " option for: " << opt
//                          << " given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " usingA " <<" usingdthdr " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_dmdr = method_dmdr;
//
//        // evolution eq
//        opt = "method_dgdr";
//        METHOD_dgdr method_dgdr;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_dgdr = METHOD_dgdr::iour;
//        }
//        else{
//            if(opts.at(opt) == "our")
//                method_dgdr = METHOD_dgdr::iour;
//            else if(opts.at(opt) == "peer")
//                method_dgdr = METHOD_dgdr::ipeer;
//            else{
//                std::cerr << AT << " option for: " << opt
//                          << " given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " our " <<" peer " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_dgdr = method_dgdr;
//
//        // set parameters for lateral expanding
//        opt = "method_spread";
//        LatSpread::METHODS method_spread;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_spread = LatSpread::iNULL;
//        }
//        else{
//            if(opts.at(opt) == "None")
//                method_spread = LatSpread::iNULL;
//            else if(opts.at(opt) == "AFGPY")
//                method_spread = LatSpread::iAFGPY;
//            else if(opts.at(opt) == "Adi")
//                method_spread = LatSpread::iAdi;
//            else if(opts.at(opt) == "AA")
//                method_spread = LatSpread::iAA;
//            else{
//                std::cerr << " option for: " << opt
//                          <<" given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " None " <<" AFGPY " << " Adi " << " AA " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getSpread()->setPars(a,theta_max,latStruct.m_theta_c,
//                                    latStruct.m_theta_w, method_spread);
//
//
//        // set parameters for EOS
//        opt = "method_eos";
//        EOSadi::METHODS method_eos;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            method_eos = EOSadi::iNava13;
//        }
//        else{
//            if(opts.at(opt) == "Nava13")
//                method_eos = EOSadi::iNava13;
//            else if(opts.at(opt) == "Peer12")
//                method_eos = EOSadi::iPeer12;
//            else{
//                std::cerr<< " option for: " << opt
//                         <<" given: " << opts.at(opt)
//                         << " is not recognized \n"
//                         << " Possible options: "
//                         << " Nava13 " << " Peer12 " << "\n"
//                         << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getEos()->setMagPars(method_eos);
//
//        /// method for shock radius
//        opt = "method_Rsh";
//        METHOD_RSh m_method_r_sh;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            m_method_r_sh = METHOD_RSh::isameAsR;
//        }
//        else{
//            if(opts.at(opt) == "sameAsR")
//                m_method_r_sh = METHOD_RSh::isameAsR;
//            else if(opts.at(opt) == "useGammaSh")
//                m_method_r_sh = METHOD_RSh::iuseGammaSh;
//            else{
//                std::cerr << " option for: " << opt
//                          <<" given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " sameAsR " << " useGammaSh " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_r_sh = m_method_r_sh;
//
//        /// method for shock velocity
//        opt = "method_GammaSh";
//        METHOD_GammaSh m_method_gamma_sh;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            m_method_gamma_sh = METHOD_GammaSh::isameAsGamma;
//        }
//        else{
//            if(opts.at(opt) == "useGamma")
//                m_method_gamma_sh = METHOD_GammaSh::isameAsGamma;
//            else if(opts.at(opt) == "useGammaRel")
//                m_method_gamma_sh = METHOD_GammaSh::iuseGammaRel;
//            else if(opts.at(opt) == "useJK")
//                m_method_gamma_sh = METHOD_GammaSh::iuseJK;
//            else if(opts.at(opt) == "useJKwithGammaRel")
//                m_method_gamma_sh = METHOD_GammaSh::iuseJKwithGammaRel;
//            else{
//                std::cerr << " option for: " << opt
//                          <<" given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " useGamma " << " useGammaRel " << " useJK "<< " useJKwithGammaRel " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_gamma_sh = m_method_gamma_sh;
//
//        /// method for energy density behind shock
//        opt = "method_Up";
//        METHODS_Up m_method_up;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            m_method_up = METHODS_Up::iuseEint2;
//        }
//        else{
//            if(opts.at(opt) == "useEint2")
//                m_method_up = METHODS_Up::iuseEint2;
//            else if(opts.at(opt) == "useGamma")
//                m_method_up = METHODS_Up::iuseGamma;
//            else{
//                std::cerr << " option for: " << opt
//                          <<" given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " useEint2 " << " useGamma " << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_up = m_method_up;
//
//        /// method for shock thickness
//        opt = "method_Delta";
//        METHOD_Delta m_method_delta;
//        if ( opts.find(opt) == opts.end() ) {
//            std::cerr << " Option for '" << opt << "' is not set. Using default value.\n";
//            m_method_delta = METHOD_Delta::iuseJoh06;
//        }
//        else{
//            if(opts.at(opt) == "useJoh06")
//                m_method_delta = METHOD_Delta::iuseJoh06;
//            else if(opts.at(opt) == "useVE12")
//                m_method_delta = METHOD_Delta::iuseVE12;
//            else if(opts.at(opt) == "None")
//                m_method_delta = METHOD_Delta::iNoDelta;
//            else{
//                std::cerr << " option for: " << opt
//                          <<" given: " << opts.at(opt)
//                          << " is not recognized \n"
//                          << " Possible options: "
//                          << " useJoh06 " << " useVE12 " << "None" << "\n"
//                          << " Exiting...\n";
//                std::cerr << AT << "\n";
//                exit(1);
//            }
//        }
//        bw_obj.getPars()->m_method_Delta = m_method_delta;
//
//        /// EATS method guverns the BW discretization (piece-wise versus adaptive)
//        bw_obj.getPars()->m_method_eats = latStruct.m_method_eats;
//        bw_obj.getPars()->nlayers = latStruct.nlayers;
//
//        /// set boolean pars
//        bw_obj.getPars()->use_dens_prof_behind_jet_for_ejecta =
//                getBoolOpt("use_dens_prof_behind_jet_for_ejecta", opts, AT,p_log,false);
//
//        bw_obj.getPars()->use_exp_rho_decay_as_floor =
//                getBoolOpt("use_exp_rho_decay_as_floor", opts, AT,p_log,false);
//
//        bw_obj.getPars()->use_flat_dens_floor =
//                getBoolOpt("use_flat_dens_floor", opts, AT,p_log,false);
//
//        bw_obj.getPars()->use_st_dens_profile =
//                getBoolOpt("use_st_dens_profile", opts, AT,p_log,false);
//
//        bw_obj.getPars()->use_bm_dens_profile =
//                getBoolOpt("use_bm_dens_profile", opts, AT,p_log,false);
//
//        bw_obj.getPars()->adiabLoss =
//                getBoolOpt("use_adiabLoss", opts, AT,p_log,true);
//
//        /// set sedov-taylor profile (for jet to be seen by ejecta as it moves behind)
//        if (bw_obj.getPars()->use_st_dens_profile) {
//            bw_obj.getSedov()->setPars(1.5, 3, 0.); // TODO this should not be here and evaluated for EVERY bw...
//            bw_obj.getSedov()->evaluate();
//        }
//
//        /// set parameters for computing observed emission
//
//        // set initials and costants for the blast wave
//        switch (latStruct.m_method_eats) {
//
//            case LatStruct::i_pw:
//                if ( latStruct.dist_E0_pw.empty() || latStruct.dist_M0_pw.empty() || latStruct.dist_G0_pw.empty() ||
//                     latStruct.theta_pw.empty()){
//                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
//                    exit(1);
//                }
//                (*p_log)(LOG_INFO,AT)<<"Init. [pw] "
//                    << " E0="<<latStruct.dist_E0_pw[ilayer]
//                    << " M0="<<latStruct.dist_M0_pw[ilayer]
//                    << " G0="<<latStruct.dist_G0_pw[ilayer]
//                    << " beta0="<<EQS::Beta(latStruct.dist_G0_pw[ilayer])
//                    << " tb0="<<t_grid[0]
//                    << " thetab0="<<latStruct.m_theta_w
//                    << " theta0="<<latStruct.theta_pw[ilayer]
//                    << " theta1="<<latStruct.theta_pw[ilayer]
//                    << " theta_w="<<latStruct.m_theta_w
//                    << " ii_eq="<<ii_eq
//                    << " ncells="<<latStruct.ncells
//                    << "\n";
//                bw_obj.setPars(latStruct.dist_E0_pw[ilayer], latStruct.dist_M0_pw[ilayer],
//                               latStruct.dist_G0_pw[ilayer], t_grid[0], 0.,
//                               latStruct.m_theta_w, latStruct.theta_pw[ilayer],
//                               latStruct.theta_pw[ilayer],//latStruct.cthetas0[ilayer],
//                               latStruct.m_theta_w, theta_max, //latStruct.cthetas0[ilayer],
//                               epsilon_e_rad, ii_eq, (double) latStruct.ncells);
//                break;
//            case LatStruct::i_adap:
//                if ( latStruct.dist_E0_a.empty() || latStruct.dist_M0_a.empty() || latStruct.dist_G0_a.empty() ||
//                    latStruct.thetas_c_h.empty()){
//                    (*p_log)(LOG_ERR, AT) << "one of the blast-wave initial data arrays is empty. \n";
//                    exit(1);
//                }
//                (*p_log)(LOG_INFO,AT)<<"Init. [a] "
//                                     << " E0="<<latStruct.dist_E0_a[ilayer]
//                                     << " M0="<<latStruct.dist_M0_a[ilayer]
//                                     << " G0="<<latStruct.dist_G0_a[ilayer]
//                                     << " beta0="<<EQS::Beta(latStruct.dist_G0_a[ilayer])
//                                     << " tb0="<<t_grid[0]
//                                     << " thetab0="<<latStruct.thetas_c_h[ilayer]
//                                     << " theta0="<<latStruct.thetas_c_l[ilayer]
//                                     << " theta1="<<latStruct.thetas_c_h[ilayer]
//                                     << " theta_w="<<latStruct.m_theta_w
//                                     << " ii_eq="<<ii_eq
//                                     << " ncells="<<1.
//                                     << "\n";
//                double fac = 2 * std::sin(0.5 * latStruct.thetas_c_h[ilayer]) * std::sin(0.5 * latStruct.thetas_c_h[ilayer]);
//                bw_obj.setPars(latStruct.dist_E0_a[ilayer],
//                               latStruct.dist_M0_a[ilayer],
//                               latStruct.dist_G0_a[ilayer], t_grid[0],
//                               0.,
//                               latStruct.thetas_c_h[ilayer],
//                               latStruct.thetas_c_l[ilayer],
//                               latStruct.thetas_c_h[ilayer],//latStruct.cthetas0[ilayer],
//                               0., theta_max,
////                               latStruct.cthetas0[ilayer],
//                               epsilon_e_rad, ii_eq, 1.);
//                break;
//        }
//
////        std::cout << " ["<<bw_obj.getPars()->ishell<<", "<<bw_obj.getPars()->ilayer<<"] "
////                  <<" G0="<<bw_obj.getPars()->Gamma0
////                  <<" E0="<<bw_obj.getPars()->E0
////                  <<" M0="<<bw_obj.getPars()->M0
////                  <<" ctheta0="<<bw_obj.getPars()->ctheta0
////                  <<"\n";
//    }

//private:
//    Array t_grid;
//    std::unique_ptr<Magnetar> p_magnetar;
//    std::vector<std::unique_ptr<RadBlastWave>> p_bws_jet;
//    std::vector<std::unique_ptr<RadBlastWave>> p_bws_ej;
//    std::vector<std::unique_ptr<CumulativeShell>> p_ej;
};
#endif
#endif //SRC_MODEL_H
