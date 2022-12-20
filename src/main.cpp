///////////////////////////////////////////////////////////////////////////////
// PyBlastAfterglowMag
// Copyright (C) 2022, Vsevolod Nedora <vsevolod.nedora@gmail.com>
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
//
//
// Compiling
// . $EDITOR main.cfg
// . g++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5
// . h5c++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5
// . gcc main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_cpp -lhdf5 -lstdc++
// | g++ main.cpp -o pba.out -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp -v
//
// Usage
// . ./pba.out ../tst/dynamics/parfile.par
// All of the analysis output will be done in the current directory
//
// Assumptions
// . Many...
///////////////////////////////////////////////////////////////////////////////

// 12/2022
// Version 1 (base) Adding the content of PyBlastAfterglow lib here
// Issues:

// include necessary files
#include <algorithm> // for std::copy
#include <unordered_map>
#include <fstream>   // for std::ifstream
#include <iostream>  // for std::cout
#include <iterator>  // for std::ostream_iterator
#include <map>       // for std::map
#include <string>    // for std::string
#include <vector>    // for std::vector
#include <cctype>
#include "omp.h"
#include <stdexcept>
#include <experimental/filesystem>// for std::filesystem::exists

#include "H5Easy.h" // "external" lib with simple read/write funcs for hdf files
#include "Logger.h"

//#include "hdf5.h"
// #include "H5Cpp.h"
// #include "H5LTpublic.h"

// #include <hdf5.h>
// #include <hdf5_hl.h>


// -------------- Code configuration ------------------
struct {
    struct {
        double n_vals_ode;
    } dynamics;
} const config =
#include "main.cfg"

// ------------- Main Model Parameters -----------------
struct Parameters{
    std::string m_path_to_ejecta_id;
};

// -------------- Read H5 table with ID ------------------
#define MYH5CHECK(ierr) if (ierr < 0) { throw std::runtime_error("Error in the HDF5 library"); }
class ReadH5ThetaVinfCorrelationFile{
    std::vector<double> m_vel_inf;
    std::vector<double> m_theta;
    std::vector<double> m_mass;
    std::unique_ptr<Logger> p_log;
public:
    void loadTable(Parameters & pars, int loglevel){
        p_log = std::make_unique<Logger>(std::cout, std::cerr, loglevel, "ReadH5ThetaVinfCorrelationFile");
        auto path_to_table = pars.m_path_to_ejecta_id;
//        path_to_table = "../../tst/dynamics/corr_vel_inf_theta.h5";
        if (!std::experimental::filesystem::exists(path_to_table))
            throw std::runtime_error("File not found. " + path_to_table);

        LoadH5 ldata;
        ldata.setFileName(path_to_table);
        ldata.setVarName("vel_inf");
        std::vector<double> vinf = ldata.getData();

        ldata.setVarName("theta");
        std::vector<double> theta = ldata.getData();

        ldata.setVarName("mass");
        std::vector<std::vector<double>> mass = ldata.getData2Ddouble();

        (*p_log)(LOG_INFO, AT)
            << "h5table loaded [vinf="<<vinf.size()<<", theta="
            << theta.size()<<", mass="<<mass.size()<< "x" << mass[0].size()<<"] \n";

        std::cout << "[Info] Ejecta ID loaded.\n";
    }

};

// ---------- Read text file with model parameters -----
void readParFile(std::unordered_map<std::string, double> pars,
                 std::unordered_map<std::string, std::string> opts,
                 std::string parfile_path, Parameters & parameters){

    if (!std::experimental::filesystem::exists(parfile_path))
        throw std::runtime_error("File not found. " + parfile_path);
    else
        std::cout << "[Info] Parfile loaded.\n";

    std::string key_after_which_to_look_for_parameters = "* Model Parameters";
    std::string key_after_which_to_look_for_settings = "* Model Settings";
    char char_that_separaters_name_and_value = '=';

    std::ifstream fin(parfile_path);
    std::string line;

    bool reading_pars = false;
    bool reading_opts = false;
    while (std::getline(fin, line)) {
        if (line == key_after_which_to_look_for_parameters) {
            reading_pars = true; reading_opts = false;
        }
        if (line == key_after_which_to_look_for_settings) {
            reading_pars = false; reading_opts = true;
        }

        if (reading_pars and (line.length() > 1) and (line != key_after_which_to_look_for_parameters)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1),
                    par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
            double value = std::stod(val);
            pars.insert(std::pair<std::string, double>(par, value));
        }

        if (reading_opts and (line.length() > 1) and (line != key_after_which_to_look_for_settings)) {
            unsigned long pos = line.find_first_of(char_that_separaters_name_and_value);
            std::string val = line.substr(pos + 1),
                    par = line.substr(0, pos);
            par.erase(std::remove_if(par.begin(), par.end(), ::isspace), par.end());
            val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
//            double value = std::stod(val);
            opts.insert(std::pair<std::string, std::string>(par, val));
        }
    }
    std::cout << "[Info] --------------------------------------- \n";
    std::cout << "[Info] Model parameters:\n";
    for (auto &&[ch, v]: pars) {
        std::cout << '\t' << ch << char_that_separaters_name_and_value << v << "\n";
    }
    std::cout << "[Info] --------------------------------------- \n";
    std::cout << "[Info] Model options:\n";
    for (auto &&[ch, v]: opts) {
        std::cout << '\t' << ch << char_that_separaters_name_and_value << v << "\n";
    }
    std::cout << "[Info] --------------------------------------- \n";

    // -------------------------------------------------------------
    parameters.m_path_to_ejecta_id = opts.at("path_to_ejecta_id");
}



// driver function
int main(int argc, char** argv) {

    /// get path to parfile
    std::unique_ptr<Logger>(p_log);
    p_log = std::make_unique<Logger>(std::cout, std::cerr, CurrLogLevel, "main");

    std::string parfile_path;
    if (argc==1){
        parfile_path = "../../tst/dynamics/parfile.par";
        (*p_log)(LOG_WARN,AT) << "Path to parfile is not given. Using debug file: " <<parfile_path<<"\n";
    }
    else if (argc>2){
        std::cerr << "args="<<argc<<"\n";
        std::cerr << argv[0] << " " << argv[1] << "\n";
        (*p_log)(LOG_WARN,AT) << "Code requires 1 argument (path to parfile). Given: " << argc << " parameters\n";
//        throw std::invalid_argument("Code requires 1 argument (path to parfile)");
    }
    else{
        parfile_path = argv[1]; // "../../tst/dynamics/parfile.par"
    }


    /// read parfile
    std::unordered_map<std::string, double> pars;
    std::unordered_map<std::string, std::string> opts;
    Parameters parameters;
    readParFile(pars, opts, parfile_path, parameters);

    ReadH5ThetaVinfCorrelationFile ejecta_id;
    ejecta_id.loadTable(parameters, CurrLogLevel);

//    std::cout << config.dynamics.n_vals_ode << "\n";

}
