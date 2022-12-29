///////////////////////////////////////////////////////////////////////////////
// PyBlastAfterglowMag
// Copyright (C) 202-2028, Vsevolod Nedora <vsevolod.nedora@gmail.com>
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
// | g++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp -vv
// | g++ `pkg-config --cflags hdf5-serial` main.cpp `pkg-config --libs hdf5-serial` -lhdf5_cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp
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

// -------------- Code configuration ------------------
struct {
    struct {
        double n_vals_ode;
    } dynamics;
} const config =
#include "main.cfg"

// -------------- Read H5 table with ID ------------------
class ReadH5ThetaVinfCorrelationFile{
    Vector m_vel_inf;
    Vector m_theta;
    VecVector m_ek;
    std::unique_ptr<logger> p_log;
public:
    void loadTable(std::string path_to_table, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "ReadH5ThetaVinfCorrelationFile");
//        auto path_to_table = pars.m_path_to_ejecta_id;
//        path_to_table = "../../tst/dynamics/corr_vel_inf_theta.h5";
        if (!std::experimental::filesystem::exists(path_to_table))
            throw std::runtime_error("File not found. " + path_to_table);

        LoadH5 ldata;
        ldata.setFileName(path_to_table);
        ldata.setVarName("vel_inf");
        m_vel_inf = ldata.getData();

        ldata.setVarName("theta");
        m_theta = ldata.getData();

        ldata.setVarName("ek");
        m_ek = ldata.getData2Ddouble();

        (*p_log)(LOG_INFO, AT)
            << "h5table loaded [vinf="<<m_vel_inf.size()<<", theta="
            << m_theta.size()<<", ek="<<m_ek.size()<< "x" << m_ek[0].size()<<"] \n";

        (*p_log)(LOG_INFO, AT) << "Ejecta ID loaded.\n";

    }
    Vector & getVelInf(){ return m_vel_inf; }
    Vector & getTheta(){ return m_theta; }
    VecVector & getEk(){ return m_ek; }
};

void readParFile2(std::unordered_map<std::string, double> & pars,
                 std::unordered_map<std::string, std::string> & opts,
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

    std::unique_ptr<logger> p_log;
    p_log = std::make_unique<logger>(std::cout, std::cerr, CurrLogLevel, "readParFile2");
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

struct Timer {

    //std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::system_clock::time_point start, end;
    std::chrono::duration<float> duration;

    Timer() {
        start = std::chrono::high_resolution_clock::now();
    }

    auto checkPoint() const{
        using namespace std::literals::chrono_literals;
        auto m_end = std::chrono::high_resolution_clock::now();
        auto m_duration = m_end-start; // duration in seconds
        auto duration_ms = m_duration.count()/1.e9;// * 1000.0f;
        return duration_ms; // s ???
//        std::cout << "Time took " << duration_ms << " ms" << std::endl;
    }

//    ~Timer()
//    {
//        end = std::chrono::high_resolution_clock::now();
//        duration = end-start; // duration in seconds
//        float duration_ms = duration.count() * 1000.0f;
//        std::cout << "Time took " << duration_ms << " ms" << std::endl;
//    }

};

// driver function
int main(int argc, char** argv) {



    int loglevel = CurrLogLevel;

    /// get path to parfile
    std::unique_ptr<logger>(p_log);
    p_log = std::make_unique<logger>(std::cout, std::cerr, CurrLogLevel, "main");
    Timer timer;

    std::string working_dir;
//    std::string parfile_arrs_path;
//    if (argc==2){
//        std::cerr << "args="<<argc<<"\n";
//        std::cerr << argv[0] << " " << argv[1] << "\n";
//        (*p_log)(LOG_WARN,AT) << "Code requires 2 argument (paths to parfile, arrs). Given: " << argc << " pr\n";
////        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
//    }
    if (argc==1){
        working_dir = "../../tst/grbafg_gauss_offaxis/";
//        parfile_arrs_path = "../../tst/dynamics/parfile_arrs.h5";
        (*p_log)(LOG_WARN,AT) << "Paths to working dir is not given. Using "
                                         << working_dir <<"\n";
    }
    else if (argc>2){
        std::cerr << "args="<<argc<<"\n";
        std::cerr << argv[0] << " " << argv[1] << "\n";
        (*p_log)(LOG_WARN,AT) << "Code requires 1 argument (path to working dir). Given: " << argc << " pr\n";
        exit(1);
//        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
    }
    else{
        working_dir = argv[1]; // "../../tst/dynamics/parfile.par"
//        parfile_arrs_path = argv[2]; // "../../tst/dynamics/parfile_arrs.par"
    }

    std::string parfilename = "parfile.par";

    /// initialize the model
    PyBlastAfterglow pba(loglevel);

    /// read main parameters of the model
    StrDbMap main_pars;
    StrStrMap main_opts;
    readParFile2(main_pars, main_opts, working_dir + parfilename,
                 "# -------------------------- main ---------------------------",
                 "# --------------------------- END ---------------------------");
    pba.setModelPars(main_pars, main_opts);

    /// observer times and frequencies
    Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
    Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
    Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
    Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

    /// read main parameters of the magnetar # TODO
//    StrDbMap mag_pars;
//    StrStrMap mag_opts;
//    readParFile2(mag_pars, mag_opts, working_dir + parfilename,
//                 "# ------------------------ Magnetar -------------------------",
//                 "# --------------------------- END ---------------------------");
//    pba.setMagnetarPars(mag_pars, mag_opts);

    /// read grb afterglow parameters
    bool run_jet_bws = getBoolOpt("run_jet_bws", main_opts, AT, p_log, false, true);
    bool save_j_dynamics = getBoolOpt("save_j_dynamics", main_opts, AT, p_log, false, true);
    bool do_j_ele = getBoolOpt("do_j_ele", main_opts, AT, p_log, false, true);
    bool do_j_spec = getBoolOpt("do_j_spec", main_opts, AT, p_log, false, true);
    bool do_j_lc = getBoolOpt("do_j_lc", main_opts, AT, p_log, false, true);
    bool do_j_skymap = getBoolOpt("do_j_skymap", main_opts, AT, p_log, false, true);
    StrDbMap grb_pars; StrStrMap grb_opts;
    for (auto & key : {"n_ism","d_l","z","theta_obs"})
        grb_pars[key] = main_pars.at(key);
    if (run_jet_bws) {
        readParFile2(grb_pars, grb_opts, working_dir+parfilename,
                     "# ---------------------- GRB afterglow ----------------------",
                     "# --------------------------- END ---------------------------");
        pba.setJetStructAnalytic(grb_pars, grb_opts);
        pba.setJetBwPars(grb_pars, grb_opts);
    }


    /// read kn afterglow parameters
    bool run_ejecta_bws = getBoolOpt("run_ejecta_bws", main_opts, AT, p_log, false, true);
    bool save_ej_dynamics = getBoolOpt("save_ej_dynamics", main_opts, AT, p_log, false, true);
    bool do_ej_ele = getBoolOpt("do_ej_ele", main_opts, AT, p_log, false, true);
    bool do_ej_spec = getBoolOpt("do_ej_spec", main_opts, AT, p_log, false, true);
    bool do_ej_lc = getBoolOpt("do_ej_lc", main_opts, AT, p_log, false, true);
    bool do_ej_skymap = getBoolOpt("do_ej_skymap", main_opts, AT, p_log, false, true);
    StrDbMap kn_pars; StrStrMap kn_opts;
    for (auto & key : {"n_ism","d_l","z","theta_obs"})
        grb_pars[key] = main_pars.at(key);
    if (run_ejecta_bws) {
        readParFile2(kn_pars, kn_opts, working_dir + parfilename,
                     "# ----------------------- kN afterglow ----------------------",
                     "# --------------------------- END ---------------------------");

        std::string path_to_arrs = getStrOpt("path_to_ejecta_id", kn_opts, AT, p_log, "", true);
        if (!std::experimental::filesystem::exists(path_to_arrs)){
            (*p_log)(LOG_ERR, AT) << " File not found. " + path_to_arrs << "\n";
            exit(1);
        }
        ReadH5ThetaVinfCorrelationFile dfile;
        dfile.loadTable(path_to_arrs, loglevel);
        pba.setEjectaStructNumeric(dfile.getTheta(), dfile.getVelInf(),
                                   dfile.getEk(), 2., false);
        pba.setEjectaBwPars(kn_pars, kn_opts);
    }

    (*p_log)(LOG_INFO, AT) << "Initialization finished [" << timer.checkPoint() << " s]" << "\n";

    /// evolve ODEs
    pba.run();

    (*p_log)(LOG_INFO, AT) << "evolution finished [" << timer.checkPoint() << " s]" << "\n";

    /// work on GRB afterglow
    if (save_j_dynamics)
        pba.saveJetBWsDynamics(
                working_dir+getStrOpt("fname_dyn", grb_opts, AT, p_log, "", true),
                (int)getDoublePar("save_dyn_every_it", grb_pars, AT, p_log, 1, true) );
    if (run_jet_bws and do_j_ele) {
        pba.setPreComputeJetAnalyticElectronsPars();
        (*p_log)(LOG_INFO, AT) << "jet analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";
    }
    if (run_jet_bws and (do_j_lc or do_j_skymap)) {

        if (do_j_lc) {
            pba.computeSaveJetLightCurveAnalytic(
                    working_dir + getStrOpt("fname_light_curve", grb_opts, AT, p_log, "", true),
                    lc_times, lc_freqs);
            (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
        }
        if (do_j_skymap) {
            pba.computeSaveJetSkyImagesAnalytic(
                    working_dir + getStrOpt("fname_sky_map", grb_opts, AT, p_log, "", true),
                    skymap_times, skymap_freqs);
            (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";
        }
    }


    /// work on kN afterglow
    if (save_ej_dynamics) {
        pba.saveEjectaBWsDynamics(
                working_dir + getStrOpt("fname_dyn", kn_opts, AT, p_log, "", true),
                (int) getDoublePar("save_dyn_every_it", kn_pars, AT, p_log, 1, true));
    }
    if (run_ejecta_bws and do_ej_ele) {
        pba.setPreComputeEjectaAnalyticElectronsPars();
        (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. electrons finished [" << timer.checkPoint() << " s]" << "\n";

    }
    if (run_ejecta_bws and (do_ej_lc or do_ej_skymap)) {
//        LoadH5 pars_arrs;
//        pars_arrs.setFileName(parfile_arrs_path);
//        pars_arrs.setVarName("light_curve_times");
//        Vector times = pars_arrs.getDataVDouble();
//        pars_arrs.setVarName("light_curve_freqs");
//        Vector freqs = pars_arrs.getDataVDouble();
        if (do_ej_lc) {
            pba.computeSaveJetLightCurveAnalytic(
                    working_dir + getStrOpt("fname_light_curve", kn_opts, AT, p_log, "", true),
                    lc_times, lc_freqs);
            (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";

        }
        if (do_ej_skymap) {
            pba.computeSaveEjectaSkyImagesAnalytic(
                    working_dir + getStrOpt("fname_sky_map", grb_opts, AT, p_log, "", true),
                    skymap_times, skymap_freqs);
            (*p_log)(LOG_INFO, AT) << "ejecta analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";

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
