///////////////////////////////////////////////////////////////////////////////
// PyBlastAfterglowMag
// 2019-2024, Vsevolod Nedora <vsevolod.nedora@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// MIT License for more details.
//
// You should have received a copy of the MIT License
// along with this program.
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


#define USE_OPENMP true // Switch to use openmp or not

// include necessary files
#include "utilitites/utils.h"
#include "utilitites/H5Easy.h"
#include "blastwave/blastwave.h"
#include "ejecta/ejecta.h"
#include "model_evolve.h"
#include "model_magnetar.h"
//#include "model_h"

class PyBlastAfterglow{
    struct Pars{
        double tb0{}; double tb1{}; int ntb{}; int iout{};
        Integrators::METHODS integrator = Integrators::METHODS::RK4;
        double rtol = 1e-5;
        int nmax = 100000;
        int loglevel{};
        bool do_average_solution = false;
    };
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    std::unique_ptr<EvolveODEsystem> p_model;
    std::unique_ptr<Output> p_out;
    std::unique_ptr<Magnetar> p_mag;
    std::unique_ptr<Ejecta> p_grb;
    std::unique_ptr<Ejecta> p_ej;
    std::unique_ptr<Ejecta> p_ej_pwn2;
    std::unique_ptr<PWNset> p_ej_pwn;
    Vector _t_grid;
    Vector t_grid;
    int m_loglevel;
    CommonTables commonTables;
public:
    std::unique_ptr<Magnetar> & getMag(){return p_mag;}
    std::unique_ptr<PWNset> & getEjPWN(){return p_ej_pwn;}
    std::unique_ptr<Ejecta> & getGRB(){return p_grb;}
    std::unique_ptr<Ejecta> & getEj(){return p_ej;}
    std::unique_ptr<Ejecta> & getEjPWN2(){return p_ej_pwn2;}
    Vector & getTburst(){return t_grid;}
    PyBlastAfterglow(int loglevel){
        m_loglevel = loglevel;
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PyBlastAfterglow");
        commonTables = CommonTables();
        p_mag = std::make_unique<Magnetar>(loglevel);
        p_grb = std::make_unique<Ejecta>(t_grid, commonTables, loglevel);
        p_ej  = std::make_unique<Ejecta>(t_grid, commonTables, loglevel);
        p_ej_pwn2 = std::make_unique<Ejecta>(t_grid, commonTables, loglevel);
        p_ej_pwn = std::make_unique<PWNset>(p_mag, p_ej, commonTables, loglevel); // depends on ang.struct. of the ejecta
    }

    void setModelPars(StrDbMap pars, StrStrMap opts) {

        /// check if parameters present
        p_pars->tb0 = getDoublePar("tb0",pars,AT,p_log,-1,true);//pars.at("tb0");
        p_pars->tb1 = getDoublePar("tb1",pars,AT,p_log,-1,true);//(double) pars.at("tb1");
        p_pars->ntb = (int) getDoublePar("ntb",pars,AT,p_log,-1,true);//pars.at("ntb");
        p_pars->iout = (int) getDoublePar("iout",pars,AT,p_log,-1,true);//pars.at("ntb");

        p_pars->rtol = getDoublePar("rtol",pars,AT,p_log,1e-13,true);//(double) pars.at("rtol");
        p_pars->nmax = (int)getDoublePar("nmax",pars,AT,p_log,100000,false);//(double) pars.at("rtol");

        p_pars->do_average_solution = getBoolOpt("do_average_solution",opts,AT,p_log, false, true);

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

        (*p_log)(LOG_INFO,AT) << "Computation tgrid = "
            <<"["<<_t_grid[0]<<", "<<_t_grid[_t_grid.size()-1]<<"] n="<<_t_grid.size()<<"\n";
        (*p_log)(LOG_INFO,AT) << "Output      tgrid = "
            <<"["<<t_grid[0]<<", "<<t_grid[t_grid.size()-1]<<"] n="<<t_grid.size()<<"\n";
    }

    /// run the time-evolution
    void run(){

        bool dorun = p_mag->run_magnetar || p_grb->run_bws || p_ej->run_bws || p_ej_pwn2->run_bws || p_ej_pwn->run_pwn;
        if (!dorun){
            (*p_log)(LOG_INFO, AT) << "Nothing to evolve. Skipping ODE integration\n";
            return;
        }

        p_model = std::make_unique<EvolveODEsystem>( p_mag, p_grb, p_ej, p_ej_pwn2,p_ej_pwn,
                                                     t_grid, _t_grid, p_pars->integrator,
                                                     p_pars->do_average_solution, m_loglevel );
        p_model->pIntegrator()->pPars()->rtol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->atol = p_pars->rtol;
        p_model->pIntegrator()->pPars()->nmax = p_pars->nmax;
        p_model->setInitialConditions(p_pars->tb0);
        (*p_log)(LOG_INFO, AT) << "all initial conditions set. Starting evolution\n";
        // evolve
        double dt;
        int ixx = 1; // [1 ... ... t_grid[-1]; where len(t_grid) <= len(_t_grid)
        int i_x = -1; // [-1, ... n_substeps]
        for (size_t it = 1; it < _t_grid.size(); it++){
            /// compute blastwave evolution using possible larger _t_grid[]
            i_x ++ ;
            p_model->getPars()->i_restarts = 0;
            dt = _t_grid[it] - _t_grid[it - 1];
            p_model->advanceTimeSubStep(dt, it, i_x);
            p_model->insertSolutionSubstep(it, _t_grid[it]);

            /// compute additional props. (e.g., microphsyics) for solutions of possible coarser grid t_grid[]
            if (it % p_pars->iout == 0) {
                p_model->processSolutionComputeAdditionalQuantities(ixx);
                ixx++;
                i_x = -1;
            }
        }
        p_model->checkEvolution();

        (*p_log)(LOG_INFO, AT) << "evolution is completed\n";
    }
};

/// -------------- Code configuration ------------------
struct { // TODO move here code options not to be changed by user
    struct {
        double n_vals_ode;
    } dynamics;
} const config =
#include "main.cfg"

/// -------------- Read H5 table with ID ------------------

int main(int argc, char** argv) {
    int loglevel = LOG_INFO;
    std::string working_dir = "../"; std::string parfilename = "parfile.par";
    /// ------------------------------------------------------
    if (argc<4){
        // test cases
//        working_dir = "../tst/grb/compare_jetsim_afgpy/tmp1/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/kn/knafg_nrinformed/workingdir_DD2_M135135_M0/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/kn/knafg_nrinformed/tmp1/"; parfilename = "parfile.par"; loglevel=LOG_INFO;

        // analysis
//        working_dir = "../analysis/grb/fsrs/tmp1/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
        working_dir = "../analysis/grb/fsrs/working_dirs/dyn_fs__rad_fs__num__ssa__ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../analysis/grb/fsrs/working_dirs/dyn_fs__rad_fs__num_ssa/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../analysis/grb/fsrs/working_dirs/dyn_fs__rad_fs__num__ssa__ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../analysis/grb/fsrs/working_dirs/dyn_fs__rad_rs__num__ssa__ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../analysis/grb/fsrs/working_dirs/dyn_spread_gp12/"; parfilename = "parfile.par"; loglevel=LOG_INFO;

        // test cases || GRG
//        working_dir = "../tst/grb/grbafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // gauss, skymap
//        working_dir = "../tst/grb/tophat_afgpy/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs
//        working_dir = "../tst/grb/tophat_spec/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs
//        working_dir = "../tst/grb/tophat_spec/working_dir/"; parfilename = "parfile_logn-30_th00_rsdyn-yes_rsrad-yes_ssa-yes_adi-yes_num.par"; loglevel=LOG_INFO; // tophat, lcs
//        working_dir = "../tst/grb/tophat_fsrs/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs

        // project test cases || GRB
//        working_dir = "../../PBA_projects/grbafg/tophat_ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs
//        working_dir = "../../PBA_projects/grbafg/tophat/output/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs

        // test cases || kn
//        working_dir = "../tst/kn/knafg_nrinformed/"; parfilename = "parfile.par"; loglevel=LOG_INFO; // tophat, lcs

        // From projects
//        working_dir = "../../PBA_projects/grbafg/tophat_ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../PBA_projects/grbafg/tophat/output/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../PBA_projects/grbafg/structured/output/"; parfilename = "parfile.par"; loglevel=LOG_INFO;



//        working_dir = "../projects/grb_prj/structured/output/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../projects/grb_prj/structured/output_rs/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../PBA_projects/grbafg/tophat_ssc/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_gauss_dyn_rs/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_gauss_dyn_spread/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_gauss_spec_rs/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grb/tophat/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/gauss_lcs_rs/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grb/grbafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_skymap_tophat_eats/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/grbafg_skymap_gauss_eats/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/tophat_afgpy/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/knafg_nrinformed/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/knafg_nrinformed_eats/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/magnetar/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/magnetar_bw/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_tophat_wind/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/kn/knafg_skymap/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbtophat_parallel/"; parfilename="tophat_EisoC500_Gamma0c1000_thetaC50_thetaW50_theta00_nism10_p22_epse05_epsb05_parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../projects/grbgauss_mcmc/working/"; parfilename="tophat_7549a8d74ce86fc502b087d8eb0e341656ee536a.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="tst.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/problems/"; parfilename="failed_tophat_4d4f9670289d90ea734b93aeb3ba05795defeca9.par"; loglevel=LOG_INFO;
//        working_dir="/media/vsevolod/T7/work/afterglowKentaProject/run1/afg/";parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../tst/crashes/"; parfilename="parfile.par"; loglevel=LOG_INFO;
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
    Timer timer{};
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
    pba.getEj()->setPars(kn_pars,kn_opts,working_dir,parfilename,
                         main_pars,ii_eq,nlayers_jet);


    StrDbMap pwn2_pars; StrStrMap pwn2_opts;
    readParFile2(pwn2_pars, pwn2_opts, p_log, working_dir + parfilename,
                 "# -------------------------- PWNBW --------------------------",
                 "# --------------------------- END ---------------------------");
    size_t ii_eq_ = pba.getMag()->getNeq() + pba.getGRB()->getNeq() + pba.getEj()->getNeq();
    pba.getEjPWN2()->setPars(pwn2_pars,pwn2_opts,working_dir,parfilename,
                             main_pars,ii_eq_,nlayers_jet);


    StrDbMap pwn_pars; StrStrMap pwn_opts;
    readParFile2(pwn_pars, pwn_opts, p_log, working_dir + parfilename,
                 "# --------------------------- PWN ---------------------------",
                 "# --------------------------- END ---------------------------");
    size_t ii_eq_t = pba.getMag()->getNeq() + pba.getGRB()->getNeq()
                   + pba.getEj()->getNeq() + pba.getEjPWN2()->getNeq();
    pba.getEjPWN()->setPars(pwn_pars,pwn_opts,working_dir,parfilename,
                            pba.getTburst(),main_pars,ii_eq_t);

    (*p_log)(LOG_INFO, AT) << "Finished reading parfile [" << timer.checkPoint() << " s]" << "\n";

    pba.run();

    (*p_log)(LOG_INFO, AT) << "evolution finished [" << timer.checkPoint() << " s]" << "\n";

//    pba.getMag()->saveMagnetarEvolution(); // TODO add magnetar here

    pba.getGRB()->processEvolved(main_pars, main_opts);

    pba.getEj()->processEvolved(main_pars, main_opts);

    pba.getEjPWN2()->processEvolved(main_pars, main_opts);

    pba.getEjPWN()->processEvolved(main_pars, main_opts);

}

