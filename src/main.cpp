///////////////////////////////////////////////////////////////////////////////
// PyBlastAfterglowMag
// Copyright (C) 2019-2028, Vsevolod Nedora <vsevolod.nedora@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////
// This package contains utilities for
// TODO describe
//
// Compiling
// . $EDITOR main.cfg
// . g++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5
// . h5c++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5
// . gcc main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5 -lstdc++
// | g++ src/main.cpp -o src/pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp -vv
// | g++ `pkg-config --cflags hdf5-serial` src/main.cpp `pkg-config --libs hdf5-serial` -lhdf5_cpp -o src/pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp
//
// Usage
// . ./pba.out /path/to/working/dir
//
// Requirements :
// Following files need to exists in working dir:
// - parfile.par
// - ejecta_id.h5 (if run_ejecta is set 'yes' in parfile
//
// Assumptions
// . Many...
///////////////////////////////////////////////////////////////////////////////

// 12/2022
// Version 1 (base) Adding the content of PyBlastAfterglow lib here
// Added parfile reader and setupe the project structure
// Issues:


// include necessary files
#include "utilitites/utils.h"
#include "utilitites/H5Easy.h"
#include "blastwave.h"
#include "model_ejecta.h"
#include "model_evolve.h"
#include "model_magnetar.h"

class PyBlastAfterglow{
    struct Pars{
        double tb0{}; double tb1{}; int ntb{}; int iout{};
        Integrators::METHODS integrator = Integrators::METHODS::RK4;
        double rtol = 1e-5;
        int nmax = 100000;
        int loglevel{};
    };
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    std::unique_ptr<EvolveODEsystem> p_model;
    std::unique_ptr<Output> p_out;
    std::unique_ptr<Magnetar> p_mag;
    std::unique_ptr<Ejecta> p_grb;
    std::unique_ptr<Ejecta> p_ej;
    std::unique_ptr<PWNset> p_ej_pwn;
    Vector _t_grid;
    Vector t_grid;
    int m_loglevel;
public:
    std::unique_ptr<Magnetar> & getMag(){return p_mag;}
    std::unique_ptr<PWNset> & getEjPWN(){return p_ej_pwn;}
    std::unique_ptr<Ejecta> & getGRB(){return p_grb;}
    std::unique_ptr<Ejecta> & getEj(){return p_ej;}
    Vector & getTburst(){return t_grid;}
    PyBlastAfterglow(int loglevel){
        m_loglevel = loglevel;
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PyBlastAfterglow");
        p_mag = std::make_unique<Magnetar>(loglevel);
        p_ej_pwn = std::make_unique<PWNset>(loglevel);
        p_grb = std::make_unique<Ejecta>(t_grid, loglevel);
        p_ej  = std::make_unique<Ejecta>(t_grid, loglevel);
    }

    void setModelPars(StrDbMap pars, StrStrMap opts) {

        /// check if parameters present
        p_pars->tb0 = getDoublePar("tb0",pars,AT,p_log,-1,true);//pars.at("tb0");
        p_pars->tb1 = getDoublePar("tb1",pars,AT,p_log,-1,true);//(double) pars.at("tb1");
        p_pars->ntb = (int) getDoublePar("ntb",pars,AT,p_log,-1,true);//pars.at("ntb");
        p_pars->iout = (int) getDoublePar("iout",pars,AT,p_log,-1,true);//pars.at("ntb");

        p_pars->rtol = getDoublePar("rtol",pars,AT,p_log,1e-13,true);//(double) pars.at("rtol");
        p_pars->nmax = (int)getDoublePar("nmax",pars,AT,p_log,100000,false);//(double) pars.at("rtol");

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

        (*p_log)(LOG_INFO,AT) << "Computation tgrid = ["<<_t_grid[0]<<", "<<_t_grid[_t_grid.size()-1]<<"] n="<<_t_grid.size()<<"\n";
        (*p_log)(LOG_INFO,AT) << "Output      tgrid = ["<<t_grid[0]<<", "<<t_grid[t_grid.size()-1]<<"] n="<<t_grid.size()<<"\n";

    }

    /// run the time-evolution
    void run(){

        bool dorun = p_mag->run_magnetar || p_grb->run_bws || p_ej->run_bws || p_ej_pwn->run_pwn;
        if (!dorun){
            (*p_log)(LOG_INFO, AT) << "Nothing to evolve. Skipping ODE integration\n";
            return;
        }

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

};


/// -------------- Code configuration ------------------
struct {
    struct {
        double n_vals_ode;
    } dynamics;
} const config =
#include "main.cfg"

/// -------------- Read H5 table with ID ------------------

int main(int argc, char** argv) {
    int loglevel;
    std::string working_dir; std::string parfilename;
    /// ------------------------------------------------------
    if (argc<4){
//        working_dir = "../tst/grbafg_gauss_offaxis/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_gauss_offaxis_io/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_skymap_io/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_tophat_afgpy/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
        working_dir = "../tst/knafg_nrinformed/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/knafg_nrinformed_eats/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/magnetar/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_tophat_wind/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/knafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbtophat_parallel/"; parfilename="tophat_EisoC500_Gamma0c1000_thetaC50_thetaW50_theta00_nism10_p22_epse05_epsb05_parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbgauss_mcmc/working/"; parfilename="tophat_7549a8d74ce86fc502b087d8eb0e341656ee536a.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="tst.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="failed_tophat_4d4f9670289d90ea734b93aeb3ba05795defeca9.par"; loglevel=LOG_INFO;
        std::cerr << " working directory and parfile are not given. Using default: "
                  << " workdir=" << working_dir << " parfile="<<parfilename << " loglevel="<< loglevel<<"\n";
    }
    else if (argc>4){
        std::cerr << "args="<<argc<<"\n";
        std::cerr << argv[0] << " " << argv[1] << " "<< argv[2] << " " << argv[3] <<"\n";
        std::cerr << "Code requires 3 argument (path to working dir, name of the parfile, and loglevel). Given: " << argc << " \n";
        exit(1);
//        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
    }
    else{
        working_dir = argv[1];
        parfilename = argv[2]; // "../../tst/dynamics/parfile.par"
        std::string _loglevel = argv[3]; // "../../tst/dynamics/parfile.par"
        if (_loglevel=="debug"){
            loglevel = LOG_DEBUG;
        }
        else if (_loglevel == "info"){
            loglevel = LOG_INFO;
        }
        else if (_loglevel == "warn"){
            loglevel = LOG_WARN;
        }
        else if (_loglevel == "err"){
            loglevel = LOG_ERR;
        }
        else if (_loglevel == "silent"){
            loglevel = LOG_SILENT;
        }
        else {
            std::cerr << " loglevel is not recognzed. Use one of: 'debug' 'info' 'warn' 'err' 'silent' \n";
            exit(1);
        }
//        parfile_arrs_path = argv[2]; // "../../tst/dynamics/parfile_arrs.par"
    }
    /// ------------------------------------------------------
    std::unique_ptr<logger>(p_log);
    p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "main");
    Timer timer;
    /// ------------------------------------------------------
    PyBlastAfterglow pba(loglevel);
    /// read main parameters of the model
    StrDbMap main_pars; StrStrMap main_opts;
    readParFile2(main_pars, main_opts, p_log, working_dir + parfilename,
                 "# -------------------------- main ---------------------------",
                 "# --------------------------- END ---------------------------");
    pba.setModelPars(main_pars, main_opts);


    StrDbMap mag_pars; StrStrMap mag_opts;
    readParFile2(mag_pars, mag_opts, p_log, working_dir + parfilename,
                 "# ------------------------ Magnetar -------------------------",
                 "# --------------------------- END ---------------------------");
    pba.getMag()->setPars(mag_pars,mag_opts,working_dir,parfilename,pba.getTburst());


    StrDbMap grb_pars; StrStrMap grb_opts;
    readParFile2(grb_pars, grb_opts, p_log, working_dir + parfilename,
                 "# ---------------------- GRB afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    pba.getGRB()->setPars(grb_pars,grb_opts,working_dir,
                          parfilename,main_pars,pba.getMag()->getNeq(),0);


    StrDbMap kn_pars; StrStrMap kn_opts;
    readParFile2(kn_pars, kn_opts, p_log, working_dir + parfilename,
                 "# ----------------------- kN afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    size_t ii_eq = pba.getMag()->getNeq() + pba.getGRB()->getNeq();
    size_t nlayers_jet = pba.getGRB()->nlayers();
    pba.getEj()->setPars(kn_pars,kn_opts,working_dir,parfilename,main_pars,ii_eq,nlayers_jet);


    StrDbMap pwn_pars; StrStrMap pwn_opts;
    readParFile2(pwn_pars, pwn_opts, p_log, working_dir + parfilename,
                 "# --------------------------- PWN ---------------------------",
                 "# --------------------------- END ---------------------------");
    size_t ii_eq_ = pba.getMag()->getNeq() + pba.getGRB()->getNeq() + pba.getEj()->getNeq();
    pba.getEjPWN()->setPars(pwn_pars,pwn_opts,working_dir,parfilename,
                            pba.getTburst(),ii_eq_,pba.getEj()->nlayers());

    (*p_log)(LOG_INFO, AT) << "Initialization finished [" << timer.checkPoint() << " s]" << "\n";

    pba.run();

    (*p_log)(LOG_INFO, AT) << "evolution finished [" << timer.checkPoint() << " s]" << "\n";

    pba.getMag()->saveMagnetarEvolution();

    pba.getEjPWN()->savePWNEvolution(main_pars);

    pba.getGRB()->computeAndOutputObservables(main_pars, main_opts);

    pba.getEj()->computeAndOutputObservables(main_pars, main_opts);

}

#if 0

// driver function
int main(int argc, char** argv) {


    std::string working_dir; std::string parfilename;
//    std::string parfile_arrs_path;
//    if (argc==2){
//        std::cerr << "args="<<argc<<"\n";
//        std::cerr << argv[0] << " " << argv[1] << "\n";
//        (*p_log)(LOG_WARN,AT) << "Code requires 2 argument (paths to parfile, arrs). Given: " << argc << " pr\n";
////        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
//    }
    int loglevel;
    if (argc<4){
//        working_dir = "../tst/grbafg_gauss_offaxis/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
        working_dir = "../tst/grbafg_gauss_offaxis_io/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_tophat_afgpy/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/knafg_nrinformed/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/magnetar/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_tophat_wind/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/knafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbtophat_parallel/"; parfilename="tophat_EisoC500_Gamma0c1000_thetaC50_thetaW50_theta00_nism10_p22_epse05_epsb05_parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbgauss_mcmc/working/"; parfilename="tophat_7549a8d74ce86fc502b087d8eb0e341656ee536a.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="tst.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="failed_tophat_4d4f9670289d90ea734b93aeb3ba05795defeca9.par"; loglevel=LOG_INFO;
        std::cerr << " working directory and parfile are not given. Using default: "
                  << " workdir=" << working_dir << " parfile="<<parfilename << " loglevel="<< loglevel<<"\n";
    }
    else if (argc>4){
        std::cerr << "args="<<argc<<"\n";
        std::cerr << argv[0] << " " << argv[1] << " "<< argv[2] << " " << argv[3] <<"\n";
        std::cerr << "Code requires 3 argument (path to working dir, name of the parfile, and loglevel). Given: " << argc << " \n";
        exit(1);
//        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
    }
    else{
        working_dir = argv[1];
        parfilename = argv[2]; // "../../tst/dynamics/parfile.par"
        std::string _loglevel = argv[3]; // "../../tst/dynamics/parfile.par"
        if (_loglevel=="debug"){
            loglevel = LOG_DEBUG;
        }
        else if (_loglevel == "info"){
            loglevel = LOG_INFO;
        }
        else if (_loglevel == "warn"){
            loglevel = LOG_WARN;
        }
        else if (_loglevel == "err"){
            loglevel = LOG_ERR;
        }
        else if (_loglevel == "silent"){
            loglevel = LOG_SILENT;
        }
        else {
            std::cerr << " loglevel is not recognzed. Use one of: 'debug' 'info' 'warn' 'err' 'silent' \n";
            exit(1);
        }
//        parfile_arrs_path = argv[2]; // "../../tst/dynamics/parfile_arrs.par"
    }
//    std::cout << " Starting...\n";
    std::unique_ptr<logger>(p_log);
    p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "main");
    Timer timer;
//    std::string parfilename = "parfile.par";

    /// initialize the model
    PyBlastAfterglow pba(loglevel);

    /// read main parameters of the model
    StrDbMap main_pars; StrStrMap main_opts;
    readParFile2(main_pars, main_opts, p_log, working_dir + parfilename,
                 "# -------------------------- main ---------------------------",
                 "# --------------------------- END ---------------------------");
    pba.setModelPars(main_pars, main_opts);

    /// observer times and frequencies
    bool lc_freq_to_time = getBoolOpt("lc_use_freq_to_time",main_opts,AT,p_log,false,true);
    Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
    Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
    Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
    Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

    /// read magnetar parameters of the magnetar
    StrDbMap mag_pars; StrStrMap mag_opts;
    bool run_magnetar = false, load_magnetar = false, save_magnetar = false;
    readParFile2(mag_pars, mag_opts, p_log, working_dir + parfilename,
                 "# ------------------------ Magnetar -------------------------",
                 "# --------------------------- END ---------------------------");
    if ((!mag_pars.empty()) || (!mag_opts.empty())) {
        run_magnetar = getBoolOpt("run_magnetar", mag_opts, AT, p_log, false, true);
        load_magnetar = getBoolOpt("load_magnetar", mag_opts, AT, p_log, false, true);
        save_magnetar = getBoolOpt("save_magnetar", mag_opts, AT, p_log, false, true);
        if (load_magnetar && run_magnetar) {
            (*p_log)(LOG_ERR, AT) << "Cannot run and load magnetar evolution at the same time. Chose one.\n";
            exit(1);
        }
        if (load_magnetar) {
            std::string fname_magnetar_ev = getStrOpt("fname_magnetar_evol", mag_opts, AT, p_log, "", true);
            if (!std::experimental::filesystem::exists(working_dir + fname_magnetar_ev)) {
                (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_magnetar_ev << "\n";
                exit(1);
            }
            /// load the magnetar evolution h5 file and parse the key quantity to magnetar class
            ReadMagnetarEvolutionFile dfile(working_dir + fname_magnetar_ev, loglevel);
            pba.getMag()->loadMagnetarEvolution(Magnetar::Q::itb, dfile.get("time"));
            pba.getMag()->loadMagnetarEvolution(Magnetar::Q::ildip, dfile.get("ldip"));
            pba.getMag()->loadMagnetarEvolution(Magnetar::Q::ilacc, dfile.get("lacc"));
            /// check if time grids mismatch
            if (pba.getMag()->getTbGrid()[0] > pba.getTburst()[0]) {
                (*p_log)(LOG_ERR, AT)
                        << "Magnetar timegrid is starts too late. Loaded magnetar times[0] > the one set in parfile "
                        << "(" << pba.getMag()->getTbGrid()[0] << " > " << pba.getTburst()[0] << ")\n";
                exit(1);
            }
        }
        if (run_magnetar) {
            /// pass the t_array manually as when loaded, the time grid may be different
            pba.getMag()->setPars(mag_pars, mag_opts);

        }
    }
    else {
        (*p_log)(LOG_INFO, AT) << "Magnetar is not initialized and will not be considered.\n";
    }

    /// read pwn parameters of the pwn
    StrDbMap pwn_pars; StrStrMap pwn_opts;
    bool run_pwn = false, save_pwn = false;
    readParFile2(pwn_pars, pwn_opts, p_log, working_dir + parfilename,
                 "# --------------------------- PWN ---------------------------",
                 "# --------------------------- END ---------------------------");
    if ((!pwn_pars.empty()) || (!pwn_pars.empty())) {
        run_pwn = getBoolOpt("run_pwn", pwn_opts, AT, p_log, false, true);
        save_pwn = getBoolOpt("save_pwn", pwn_opts, AT, p_log, false, true);
    }
    else {
        (*p_log)(LOG_INFO, AT) << "PWN is not initialized and will not be considered.\n";
    }
#if 0
    /// read grb afterglow parameters
    StrDbMap grb_pars; StrStrMap grb_opts;
    bool run_jet_bws=false, save_j_dynamics=false, do_j_ele=false, do_j_spec=false, do_j_lc=false, do_j_skymap=false;
    readParFile2(grb_pars, grb_opts, p_log, working_dir+parfilename,
                 "# ---------------------- GRB afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    if ((!grb_pars.empty()) || (!grb_opts.empty())) {
        run_jet_bws = getBoolOpt("run_jet_bws", grb_opts, AT, p_log, false, true);
        save_j_dynamics = getBoolOpt("save_j_dynamics", grb_opts, AT, p_log, false, true);
        do_j_ele = getBoolOpt("do_j_ele", grb_opts, AT, p_log, false, true);
        do_j_spec = getBoolOpt("do_j_spec", grb_opts, AT, p_log, false, true);
        do_j_lc = getBoolOpt("do_j_lc", grb_opts, AT, p_log, false, true);
        do_j_skymap = getBoolOpt("do_j_skymap", grb_opts, AT, p_log, false, true);
        for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
            if (main_pars.find(key) == main_pars.end()) {
                (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                exit(1);
            }
            grb_pars[key] = main_pars.at(key);
        }
        if (run_jet_bws) {
            pba.getGRB()->setJetStructAnalytic(grb_pars, grb_opts);
//        pba.setJetStructAnalytic(grb_pars, grb_opts);
            pba.getGRB()->setJetBwPars(grb_pars, grb_opts, pba.getMag()->getNeq());
        }
    }
    else {
        (*p_log)(LOG_INFO, AT) << "GRB is not initialized and will not be considered.\n";
    }
#endif

    /// read GRB afterglow parameters
    StrDbMap grb_pars; StrStrMap grb_opts;
    bool run_bws=false,save_dyn=false, do_ele=false, do_spec=false, do_lc=false, do_skymap=false;
    readParFile2(grb_pars, grb_opts, p_log, working_dir + parfilename,
                 "# ---------------------- GRB afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    if ((!grb_pars.empty()) || (!grb_opts.empty())) {
        run_bws = getBoolOpt("run_bws", grb_opts, AT, p_log, false, true);
        save_dyn = getBoolOpt("save_dynamics", grb_opts, AT, p_log, false, true);
        do_ele = getBoolOpt("do_ele", grb_opts, AT, p_log, false, true);
        do_spec = getBoolOpt("do_spec", grb_opts, AT, p_log, false, true);
        do_lc = getBoolOpt("do_lc", grb_opts, AT, p_log, false, true);
        do_skymap = getBoolOpt("do_skymap", grb_opts, AT, p_log, false, true);
        for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
            if (main_pars.find(key) == main_pars.end()) {
                (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                exit(1);
            }
            grb_pars[key] = main_pars.at(key);
        }
        grb_opts["workingdir"] = working_dir; // For loading Nuclear Heating table
        if (run_bws) {
            std::string fname_ejecta_id = getStrOpt("fname_ejecta_id", grb_opts, AT, p_log, "", true);
            bool use_1d_id = getBoolOpt("use_1d_id", grb_opts, AT, p_log, false, true);
            if (!std::experimental::filesystem::exists(working_dir + fname_ejecta_id)) {
                (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_ejecta_id << "\n";
                exit(1);
            }
            EjectaID ID(working_dir + fname_ejecta_id, loglevel);
            if (ID.idtype == EjectaID::IDTYPE::i_id_corr && (!use_1d_id))
                pba.getEj()->setEjectaStructNumeric(
                        ID,
                        getBoolOpt("enforce_angular_grid", grb_opts, AT, p_log, false, true),
                        grb_opts);
            else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (use_1d_id))
                pba.getEj()->setEjectaStructNumeric(ID, grb_opts);
            else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (!use_1d_id))
                pba.getEj()->setEjectaStructNumericUniformInTheta(
                        ID,
                        (size_t) getDoublePar("nlayers", grb_pars, AT, p_log, 30, true),
                        getDoublePar("mfac", grb_pars, AT, p_log, 1.0, true),
                        grb_opts);
            else {
                (*p_log)(LOG_ERR, AT) << " no valid grb ejecta structure given\n";
                exit(1);
            }
            size_t ii_eq = pba.getMag()->getNeq();
//            size_t nlayers_jet = pba.getGRB()->getBWs().size();
            pba.getEj()->setEjectaBwPars(grb_pars, grb_opts, ii_eq, 0);

            /// initialize Ejecta Bound PWN
            if (run_pwn) {
                (*p_log)(LOG_ERR,AT) << " not implemented\n";
                exit(1);
//                size_t ii_eq_ = pba.getMag()->getNeq() + pba.getGRB()->getNeq() + pba.getEj()->getNeq();
                /// init ejecta-bound PWN
//                pba.getEjPWN()->setPWNpars(pba.getTburst(), pwn_pars, pwn_opts, ii_eq_,
//                                           pba.getEj()->getShells().size());
            }
        }
    }
    else{
        (*p_log)(LOG_INFO, AT) << "kN is not initialized and will not be considered.\n";
    }

    /// read kn afterglow parameters
    StrDbMap kn_pars; StrStrMap kn_opts;
    bool run_ejecta_bws=false,save_ej_dynamics=false, do_ej_ele=false, do_ej_spec=false, do_ej_lc=false, do_ej_skymap=false;
    readParFile2(kn_pars, kn_opts, p_log, working_dir + parfilename,
                 "# ----------------------- kN afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    if ((!kn_pars.empty()) || (!kn_opts.empty())) {
        run_bws = getBoolOpt("run_bws", kn_opts, AT, p_log, false, true);
        save_dyn = getBoolOpt("save_dyn", kn_opts, AT, p_log, false, true);
        do_ele = getBoolOpt("do_ele", kn_opts, AT, p_log, false, true);
        do_spec = getBoolOpt("do_spec", kn_opts, AT, p_log, false, true);
        do_lc = getBoolOpt("do_lc", kn_opts, AT, p_log, false, true);
        do_skymap = getBoolOpt("do_skymap", kn_opts, AT, p_log, false, true);
        for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
            if (main_pars.find(key) == main_pars.end()) {
                (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                exit(1);
            }
            kn_pars[key] = main_pars.at(key);
        }
        kn_opts["workingdir"] = working_dir; // For loading Nuclear Heating table
        if (run_bws) {
            std::string fname_ejecta_id = getStrOpt("fname_ejecta_id", kn_opts, AT, p_log, "", true);
            bool use_1d_id = getBoolOpt("use_1d_id", kn_opts, AT, p_log, false, true);
            if (!std::experimental::filesystem::exists(working_dir + fname_ejecta_id)) {
                (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_ejecta_id << "\n";
                exit(1);
            }
            EjectaID ID(working_dir + fname_ejecta_id, loglevel);
            if (ID.idtype == EjectaID::IDTYPE::i_id_corr && (!use_1d_id))
                pba.getEj()->setEjectaStructNumeric(
                        ID,
                        getBoolOpt("enforce_angular_grid", kn_opts, AT, p_log, false, true),
                        kn_opts);
            else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (use_1d_id))
                pba.getEj()->setEjectaStructNumeric(ID, kn_opts);
            else if (ID.idtype == EjectaID::IDTYPE::i_id_hist && (!use_1d_id))
                pba.getEj()->setEjectaStructNumericUniformInTheta(
                        ID,
                        (size_t) getDoublePar("nlayers", kn_pars, AT, p_log, 30, true),
                        getDoublePar("mfac", kn_pars, AT, p_log, 1.0, true),
                        kn_opts);
            else {
                (*p_log)(LOG_ERR, AT) << " no valid ejecta structure given\n";
                exit(1);
            }
            size_t ii_eq = pba.getMag()->getNeq() + pba.getGRB()->getNeq();
            if (pba.getGRB()->nshells() > 1){
                (*p_log)(LOG_ERR,AT)<<"not implemented\n";
                exit(1);
            }
            size_t nlayers_jet =pba.getGRB()->nlayers();
            pba.getEj()->setEjectaBwPars(kn_pars, kn_opts, ii_eq, nlayers_jet);

            /// initialize Ejecta Bound PWN
            if (run_pwn) {
                size_t ii_eq_ = pba.getMag()->getNeq() + pba.getGRB()->getNeq() + pba.getEj()->getNeq();
                /// init ejecta-bound PWN
                pba.getEjPWN()->setPWNpars(pba.getTburst(), pwn_pars, pwn_opts, ii_eq_,
                                           pba.getEj()->getShells().size());
            }
        }
    }
    else{
        (*p_log)(LOG_INFO, AT) << "kN is not initialized and will not be considered.\n";
    }


    (*p_log)(LOG_INFO, AT) << "Initialization finished [" << timer.checkPoint() << " s]" << "\n";

    /// evolve ODEs
    pba.run();

    (*p_log)(LOG_INFO, AT) << "evolution finished [" << timer.checkPoint() << " s]" << "\n";

    /// save Magnetar data
    if (run_magnetar and save_magnetar){
        pba.getMag()->saveMagnetarEvolution(
                working_dir,
                getStrOpt("fname_mag", mag_opts, AT, p_log, "", true),
                (int)getDoublePar("save_mag_every_it", mag_pars, AT, p_log, 1, true));
    }

    /// save PWN
    if (run_pwn and save_pwn){
        pba.getEjPWN()->savePWNEvolution(
                working_dir,
                getStrOpt("fname_pwn", pwn_opts, AT, p_log, "", true),
                (int)getDoublePar("save_pwn_every_it", pwn_pars, AT, p_log, 1, true),
                main_pars, pwn_pars );
    }

    /// work on GRB afterglow
    if (run_bws){
        if (save_dyn)
            pba.getGRB()->saveEjectaBWsDynamics(
                    working_dir,
                    getStrOpt("fname_dyn", grb_opts, AT, p_log, "", true),
                    (int)getDoublePar("save_dyn_every_it", grb_pars, AT, p_log, 1, true),
                    main_pars, grb_pars);
        if (do_ele) {
            pba.getGRB()->setPreComputeEjectaAnalyticElectronsPars();
            (*p_log)(LOG_INFO, AT) << "jet analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";
        }
        if (do_lc or do_skymap) {
            if (do_lc) {
                pba.getGRB()->computeSaveEjectaLightCurveAnalytic(
                        working_dir,
                        getStrOpt("fname_light_curve", grb_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", grb_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, grb_pars, lc_freq_to_time);
                (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_skymap) {
                pba.getGRB()->computeSaveEjectaSkyImagesAnalytic(
                        working_dir,
                        getStrOpt("fname_sky_map", grb_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, grb_pars);
                (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";
            }
        }
    }

    /// work on kN afterglow
    if (run_bws){
        if (save_dyn)
            pba.getEj()->saveEjectaBWsDynamics(
                    working_dir,
                    getStrOpt("fname_dyn", kn_opts, AT, p_log, "", true),
                    (int) getDoublePar("save_dyn_every_it", kn_pars, AT, p_log, 1, true),
                    main_pars, grb_pars);
        if (do_ele) {
            pba.getEj()->setPreComputeEjectaAnalyticElectronsPars();
            (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";

        }
        if (do_lc or do_skymap) {
            if (do_lc) {
                pba.getEj()->computeSaveEjectaLightCurveAnalytic(
                        working_dir,
                        getStrOpt("fname_light_curve", kn_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", kn_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, grb_pars, lc_freq_to_time);
                (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_skymap) {
                pba.getEj()->computeSaveEjectaSkyImagesAnalytic(
                        working_dir,
                        getStrOpt("fname_sky_map", kn_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, grb_pars);
                (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";

            }
        }
    }
}
#endif
