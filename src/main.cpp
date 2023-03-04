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
#include "utils.h"
#include "model.h"
#include "H5Easy.h"

/// -------------- Code configuration ------------------
struct {
    struct {
        double n_vals_ode;
    } dynamics;
} const config =
#include "main.cfg"

/// -------------- Read H5 table with ID ------------------
class ReadH5ThetaVinfCorrelationFile{
    Vector m_vel_inf;
    Vector m_theta;
    VecVector m_ek_corr;
    Vector m_ek_hist;
    std::unique_ptr<logger> p_log;
public:
    enum IDTYPE { i_id_hist, i_id_corr };
    IDTYPE idtype;
    ReadH5ThetaVinfCorrelationFile(std::string path_to_table, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "ReadH5ThetaVinfCorrelationFile");
//        auto path_to_table = pars.m_path_to_ejecta_id;
//        path_to_table = "../../tst/dynamics/corr_vel_inf_theta.h5";
        if (!std::experimental::filesystem::exists(path_to_table))
            throw std::runtime_error("File not found. " + path_to_table);

        LoadH5 ldata;
        ldata.setFileName(path_to_table);
        ldata.setVarName("vel_inf");
        m_vel_inf = ldata.getData();
        int beta_size = ldata.getSize();

        ldata.setVarName("theta");
        m_theta = ldata.getData();

        ldata.setVarName("ek");
        int ek_size = ldata.getSize();


        if (ek_size!=beta_size) {
            m_ek_corr = ldata.getData2Ddouble();
            m_ek_hist = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 2D loaded [vinf="<<m_vel_inf.size()<<", theta="
                    << m_theta.size()<<", ek="<<m_ek_corr.size()<< "x" << m_ek_corr[0].size()<<"] \n";
            idtype = i_id_corr;
        }
        else{
            m_ek_hist = ldata.getDataVDouble();
            m_ek_corr = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 1D loaded [vinf="<<m_vel_inf.size()<<", theta="
                    << m_theta.size()<<", ek="<<m_ek_hist.size()<<"] \n";
            idtype = i_id_hist;
        }

        (*p_log)(LOG_INFO, AT) << "Ejecta ID loaded\n";

    }
    Vector & getVelInf(){ return m_vel_inf; }
    Vector & getTheta(){ return m_theta; }
    VecVector & getEkCorr(){ return m_ek_corr; }
    Vector & getEkHist(){ return m_ek_hist; }
};

/// ----------- Read H5 file with magnetar table -------
class ReadMagnetarEvolutionFile{
    std::unique_ptr<logger> p_log;
    LoadH5 m_ldata;
    size_t data_size = 0; // track the size of arrays to avoid mismatching
public:
    ReadMagnetarEvolutionFile(std::string fapth, int loglevel) {
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "ReadMagnetarEvolutionFile");
//        auto path_to_table = pars.m_path_to_ejecta_id;
//        path_to_table = "../../tst/dynamics/corr_vel_inf_theta.h5";
        if (!std::experimental::filesystem::exists(fapth))
            throw std::runtime_error("File not found. " + fapth);
        /// loading the datafile
        m_ldata.setFileName(fapth);
        (*p_log)(LOG_INFO, AT) << "Ejecta ID loaded\n";
    }

    Vector get(std::string varname) {
        m_ldata.setVarName(varname);
        Vector data = m_ldata.getData();
        if (data_size == 0)
            data_size =  m_ldata.getSize();
        if (m_ldata.getSize() != data_size){
            (*p_log)(LOG_ERR,AT)<<"Input data size mismatch. All arrays should be the same size."
                                <<" Given array v_n="<<varname<<" has size="<<m_ldata.getSize()
                                <<" Expected size="<<data_size<<"\n";
            exit(1);
        }
        else
            data_size = m_ldata.getSize();
        return std::move(data);
    }
};


void readParFile2(std::unordered_map<std::string, double> & pars,
                 std::unordered_map<std::string, std::string> & opts,
                 std::unique_ptr<logger> & p_log,
                 std::string parfile_path,
                 std::string from_line, std::string until_line
                 ){

    /// settings for reading the parfile
    std::string key_after_which_to_look_for_parameters = "* Parameters";
    std::string key_after_which_to_look_for_settings = "* Settings";
    char char_that_separaters_name_and_value = '=';
    char char_that_separaters_value_and_comment = '#';
    std::vector<std::string> leave_spaces_for = {
            "lc_freqs", "lc_times", "skymap_freqs", "skymap_times"
    };

//    std::unique_ptr<logger> p_log;
//    p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "readParFile2");
    if (!std::experimental::filesystem::exists(parfile_path)) {
        (*p_log)(LOG_ERR, AT) << " Parfile not found. " + parfile_path << "\n";
        exit(1);

    }
    std::ifstream fin(parfile_path);
    std::string line;

    bool is_in_the_reqired_block = false;
    bool reading_pars = false;
    bool reading_opts = false;
    while (std::getline(fin, line)) {
        /// check if reading the required block of parfile
        if (line == from_line)
            is_in_the_reqired_block = true;
        if (line == until_line)
            is_in_the_reqired_block = false;
        if (!is_in_the_reqired_block)
            continue;

        /// read parameters (double) and settings (str) separately
        if (line == key_after_which_to_look_for_parameters) {
            reading_pars = true; reading_opts = false;
        }
        if (line == key_after_which_to_look_for_settings) {
            reading_pars = false; reading_opts = true;
        }
        if (line[0] == char_that_separaters_value_and_comment)
            continue;

        /// read pars (str, double)
        if (reading_pars and (line.length() > 1) and (line != key_after_which_to_look_for_parameters)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1), par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
//            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            if (val.find(char_that_separaters_value_and_comment) != std::string::npos){
                unsigned long _pos = val.find_first_of(char_that_separaters_value_and_comment);
                std::string _comment = val.substr(_pos + 1), _val = val.substr(0, _pos);
                val = _val;
            }
            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            double value = std::stod(val);
            pars.insert(std::pair<std::string, double>(par, value));
        }
        /// read opts (str, str)
        if (reading_opts and (line.length() > 1) and (line != key_after_which_to_look_for_settings)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1),
                    par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
//            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            if (val.find(char_that_separaters_value_and_comment) != std::string::npos){
                unsigned long _pos = val.find_first_of(char_that_separaters_value_and_comment);
                std::string _comment = val.substr(_pos + 1), _val = val.substr(0, _pos);
                val = _val;
            }
            if (std::find(leave_spaces_for.begin(), leave_spaces_for.end(), par) == leave_spaces_for.end())
                val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
//            double value = std::stod(val);
            opts.insert(std::pair<std::string, std::string>(par, val));
        }
    }

}

Vector makeVecFromString(std::string line, std::unique_ptr<logger> & p_log){
    char space_char = ' ';
//    char end_char = '#';
    std::vector<std::string> words{};
    std::stringstream sstream(line);
    std::string word;
    while (std::getline(sstream, word, space_char)){
        word.erase(std::remove_if(word.begin(), word.end(), ispunct), word.end());
        words.push_back(word);
    }
    /// remove first emtpy element left after '= ' this
    if ((words[0].empty())or(words[0] == " ")){
        words.erase(words.begin());
    }

    /// construct the vector
    Vector res;
    if (words[0] == "array"){
        if (words[1] == "logspace"){
            const double val1 = std::log10( std::stod(words[2]) );
            const double val2 = std::log10( std::stod(words[3]) );
            const int nvals = (int)std::stod(words[4]);
            res = TOOLS::MakeLogspaceVec(val1,val2,nvals,10);
        }
        else{
            words.erase(words.begin());
            for (const auto &str : words){
                res.push_back(std::stod(str));
            }
        }
    }
    else{
        (*p_log)(LOG_ERR, AT) << "incorrect first word in array line. Expected 'array' found " << words[0] << '\n';
        exit(1);
    }
//    std::cout << res << "\n";
    return std::move( res );
}



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
//        working_dir = "../../tst/grbafg_gauss_offaxis/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/grbafg_tophat_afgpy/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
//        working_dir = "../../tst/knafg_nrinformed/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
        working_dir = "../../tst/magnetar/"; parfilename = "parfile.par"; loglevel=LOG_INFO;
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

    /// read magnetar parameters of the magnetar # TODO
    StrDbMap mag_pars; StrStrMap mag_opts;
    readParFile2(mag_pars, mag_opts, p_log, working_dir + parfilename,
                 "# ------------------------ Magnetar -------------------------",
                 "# --------------------------- END ---------------------------");
    bool run_magnetar = getBoolOpt("run_magnetar", mag_opts, AT, p_log, false, true);
    bool load_magnetar = getBoolOpt("load_magnetar", mag_opts, AT, p_log, false, true);
    bool save_magnetar = getBoolOpt("save_magnetar", mag_opts, AT, p_log, false, true);
    if (load_magnetar && run_magnetar){
        (*p_log)(LOG_ERR,AT)<<"Cannot run and load magnetar evolution at the same time. Chose one.\n";
        exit(1);
    }
    if (load_magnetar){
        std::string fname_magnetar_ev = getStrOpt("fname_magnetar_evol", mag_opts, AT, p_log, "", true);
        if (!std::experimental::filesystem::exists(working_dir+fname_magnetar_ev)){
            (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir+fname_magnetar_ev << "\n";
            exit(1);
        }
        /// load the magnetar evolution h5 file and parse the key quantity to magnetar class
        ReadMagnetarEvolutionFile dfile(working_dir+fname_magnetar_ev, loglevel);
        pba.getMag()->loadEvolution(Magnetar::Q::itb, dfile.get("tburst"));
        pba.getMag()->loadEvolution(Magnetar::Q::iLtot, dfile.get("ltot"));
    }
    if (run_magnetar) {
        /// pass the t_array manually as when loaded, the time grid may be different
        pba.getMag()->setPars(pba.getTburst(),mag_pars, mag_opts);
    }


    /// read grb afterglow parameters
    StrDbMap grb_pars; StrStrMap grb_opts;
    readParFile2(grb_pars, grb_opts, p_log, working_dir+parfilename,
                 "# ---------------------- GRB afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    const bool run_jet_bws = getBoolOpt("run_jet_bws", grb_opts, AT, p_log, false, true);
    const bool save_j_dynamics = getBoolOpt("save_j_dynamics", grb_opts, AT, p_log, false, true);
    const bool do_j_ele = getBoolOpt("do_j_ele", grb_opts, AT, p_log, false, true);
    const bool do_j_spec = getBoolOpt("do_j_spec", grb_opts, AT, p_log, false, true);
    const bool do_j_lc = getBoolOpt("do_j_lc", grb_opts, AT, p_log, false, true);
    const bool do_j_skymap = getBoolOpt("do_j_skymap", grb_opts, AT, p_log, false, true);
    for (auto & key : {"n_ism","d_l","z","theta_obs","A0","s","r_ej","r_ism"}){
        if (main_pars.find(key) == main_pars.end()){
            (*p_log)(LOG_ERR,AT) << " keyword '"<<key<<"' is not found in main parameters. \n";
            exit(1);
        }
        grb_pars[key] = main_pars.at(key);
    }
    if (run_jet_bws) {
        pba.getGRB()->setJetStructAnalytic(grb_pars, grb_opts);
//        pba.setJetStructAnalytic(grb_pars, grb_opts);
        pba.getGRB()->setJetBwPars(grb_pars, grb_opts, pba.getMag()->getNeq());
    }


    /// read kn afterglow parameters
    StrDbMap kn_pars; StrStrMap kn_opts;
    readParFile2(kn_pars, kn_opts, p_log, working_dir + parfilename,
                 "# ----------------------- kN afterglow ----------------------",
                 "# --------------------------- END ---------------------------");
    const bool run_ejecta_bws = getBoolOpt("run_ej_bws", kn_opts, AT, p_log, false, true);
//    if (run_ejecta_bws) {
//        for (const auto& i : kn_opts)
//            std::cout << i.first << "=" << i.second << "|" << std::endl;
//        std::cerr << AT << "\n"; exit(1);
//    }
    const bool save_ej_dynamics = getBoolOpt("save_ej_dynamics", kn_opts, AT, p_log, false, true);
    const bool do_ej_ele = getBoolOpt("do_ej_ele", kn_opts, AT, p_log, false, true);
    const bool do_ej_spec = getBoolOpt("do_ej_spec", kn_opts, AT, p_log, false, true);
    const bool do_ej_lc = getBoolOpt("do_ej_lc", kn_opts, AT, p_log, false, true);
    const bool do_ej_skymap = getBoolOpt("do_ej_skymap", kn_opts, AT, p_log, false, true);
    for (auto & key : {"n_ism","d_l","z","theta_obs","A0","s","r_ej","r_ism"}){
        if (main_pars.find(key) == main_pars.end()){
            (*p_log)(LOG_ERR,AT) << " keyword '"<<key<<"' is not found in main parameters. \n";
            exit(1);
        }
        kn_pars[key] = main_pars.at(key);
    }
    if (run_ejecta_bws) {

        std::string fname_ejecta_id = getStrOpt("fname_ejecta_id", kn_opts, AT, p_log, "", true);
        if (!std::experimental::filesystem::exists(working_dir+fname_ejecta_id)){
            (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir+fname_ejecta_id << "\n";
            exit(1);
        }
        ReadH5ThetaVinfCorrelationFile dfile(working_dir+fname_ejecta_id, loglevel);
//        dfile.loadTable(working_dir+fname_ejecta_id, loglevel);
        if (dfile.idtype==ReadH5ThetaVinfCorrelationFile::IDTYPE::i_id_corr)
            pba.getEj()->setEjectaStructNumeric(dfile.getTheta(), dfile.getVelInf(), dfile.getEkCorr(),
                                                getBoolOpt("enforce_angular_grid", kn_opts, AT, p_log, false, true),
                                                kn_opts);
        else if (dfile.idtype==ReadH5ThetaVinfCorrelationFile::IDTYPE::i_id_hist)
            pba.getEj()->setEjectaStructNumericUniformInTheta(dfile.getTheta(),dfile.getVelInf(),dfile.getEkHist(),
                                                              (size_t)getDoublePar("nlayers", kn_pars, AT, p_log, 30, true),
                                                              getDoublePar("mfac", kn_pars, AT, p_log, 1.0, true),
                                                              kn_opts);
        else {
            (*p_log)(LOG_ERR, AT) << " no valid ejecta structure given\n";
            exit(1);
        }
        size_t ii_eq = pba.getMag()->getNeq() + pba.getGRB()->getNeq();
        size_t nlayers_jet = pba.getGRB()->getBWs().size();
        pba.getEj()->setEjectaBwPars(kn_pars, kn_opts, ii_eq, nlayers_jet);
    }

    (*p_log)(LOG_INFO, AT) << "Initialization finished [" << timer.checkPoint() << " s]" << "\n";

    /// evolve ODEs
    pba.run();

    (*p_log)(LOG_INFO, AT) << "evolution finished [" << timer.checkPoint() << " s]" << "\n";

    /// save Magnetar data
    if (run_magnetar and save_magnetar){
        pba.saveMagnetarEvolution(
                working_dir,
                getStrOpt("fname_mag", mag_opts, AT, p_log, "", true),
                (int)getDoublePar("save_mag_every_it", mag_pars, AT, p_log, 1, true) );
    }

    /// work on GRB afterglow
    if (run_jet_bws){
        if (save_j_dynamics)
            pba.saveJetBWsDynamics(
                    working_dir,
                    getStrOpt("fname_dyn", grb_opts, AT, p_log, "", true),
                    (int)getDoublePar("save_dyn_every_it", grb_pars, AT, p_log, 1, true),
                    main_pars, grb_pars);
        if (do_j_ele) {
            pba.setPreComputeJetAnalyticElectronsPars();
            (*p_log)(LOG_INFO, AT) << "jet analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";
        }
        if (do_j_lc or do_j_skymap) {
            if (do_j_lc) {
                pba.computeSaveJetLightCurveAnalytic(
                        working_dir,
                        getStrOpt("fname_light_curve", grb_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", grb_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, grb_pars, lc_freq_to_time);
                (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_j_skymap) {
                pba.computeSaveJetSkyImagesAnalytic(
                        working_dir,
                        getStrOpt("fname_sky_map", grb_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, grb_pars);
                (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";
            }
        }
    }

    /// work on kN afterglow
    if (run_ejecta_bws){
        if (save_ej_dynamics)
            pba.saveEjectaBWsDynamics(
                    working_dir,
                    getStrOpt("fname_dyn", kn_opts, AT, p_log, "", true),
                    (int) getDoublePar("save_dyn_every_it", kn_pars, AT, p_log, 1, true),
                    main_pars, grb_pars);
        if (do_ej_ele) {
            pba.setPreComputeEjectaAnalyticElectronsPars();
            (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";

        }
        if (do_ej_lc or do_ej_skymap) {
            if (do_ej_lc) {
                pba.computeSaveEjectaLightCurveAnalytic(
                        working_dir,
                        getStrOpt("fname_light_curve", kn_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", kn_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, grb_pars, lc_freq_to_time);
                (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_ej_skymap) {
                pba.computeSaveEjectaSkyImagesAnalytic(
                        working_dir,
                        getStrOpt("fname_sky_map", kn_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, grb_pars);
                (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";

            }
        }
    }

}


//
//    std::unordered_map<std::string, double> pars;
//    std::unordered_map<std::string, std::string> opts;
//    loadParFile(pars, opts, parfile_path);
//
//
//
//
//    std::unordered_map<std::string, double> pars;
//    std::unordered_map<std::string, std::string> opts;
//    loadParFile(pars, opts, parfile_path);
//
//    /// initialize the model
//    PyBlastAfterglow pba(loglevel);
//
//    bool do_grb_dyn = getBoolOpt("do_grb_dyn", opts, AT, p_log, false, true);
//    if (do_grb_dyn){
//        pba.setJetStructAnalytic(pars, opts);
//    }
//
//    pba.setModelPars(pars, opts);
//
//    if (do_grb_dyn){
//        pba.setJetBwPars(pars, opts);
//    }
//
//    bool do_kn_syn_an_lc = getBoolOpt("do_kn_syn_an_lc", opts, AT, p_log, false, true);
//    bool do_kn_syn_an_skymap = getBoolOpt("do_kn_syn_an_skymap", opts, AT, p_log, false, true);
//    if (do_kn_syn_an_lc or do_kn_syn_an_skymap){
//        pba.setPreComputeEjectaAnalyticElectronsPars();
//    }
//
//    LoadH5 pars_arrs;
//    pars_arrs.setFileName(parfile_arrs_path);
//    pars_arrs.setVarName("light_curve_times");
//    Vector times = pars_arrs.getDataVDouble();
//    pars_arrs.setVarName("light_curve_freqs");
//    Vector freqs = pars_arrs.getDataVDouble();
//    if (do_kn_syn_an_lc){
//        pba.computeSaveJetLightCurveAnalytic(
//                getStrOpt("fpath_grb_light_curve", opts, AT, p_log, "", true),
//                times, freqs );
//    }

//
//    pr.m_path_to_ejecta_id = getStrOpt("path_to_ejecta_id", opts, AT, p_log, "", true);
//    ReadH5ThetaVinfCorrelationFile ejecta_id;
//    ejecta_id.loadTable(pr, CurrLogLevel);
//
//
//    ///
//    pr.m_task = getStrOpt("task", opts, AT, p_log, "", true);
//    if (pr.m_task == "evolve_gaussian_bws"){
//        (*p_log)(LOG_INFO, AT) << "Starting task: '" << pr.m_task << "'\n";
//        pr.initGaussianStructure(pars, opts);
//    }





//    std::cout << config.dynamics.n_vals_ode << "\n";

//}
