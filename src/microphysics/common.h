//
// Created by vsevolod on 4/22/24.
//

#ifndef SRC_COMMON_H
#define SRC_COMMON_H

#include "../utilitites/utils.h"
#include "../utilitites/interpolators.h"
#include "../utilitites/H5Easy.h"
#include "kernels.h"

struct EBL{
    Vector m_freq,m_z,m_tau;
    EBL() = default;
//    std::string table_fpath = "../data/EBL/Franceschini18/table.h5";
    std::string table_fpath = "../../../data/EBL/Franceschini18/table.h5";
    void load_h5_table(){
        if (!std::experimental::filesystem::exists(table_fpath)) {
            std::cerr << AT << " EBL data file not found: " + table_fpath << "\n";
            exit(1);
        }
        LoadH5 ldata;
        ldata.setFileName(table_fpath);

        ldata.setVarName("freq");
        m_freq = std::move ( ldata.getDataVDouble() ); // load 1D vec

        ldata.setVarName("z");
        m_z = std::move ( ldata.getDataVDouble() ); // load 1D vec

        ldata.setVarName("tau");
        VecVector tau_ = std::move( ldata.getData2Ddouble() ); // [freq][z]

//        m_freq.resize(freq_.size()*z_.size());
//        m_z.resize(freq_.size()*z_.size());
        m_tau.resize(m_freq.size()*m_z.size());

        if (tau_.size() != m_freq.size()){
            std::cerr << AT<< " size mismatch tau.size()="<<tau_.size()<<" != freq.size()="<<m_freq.size()<<"\n";
            exit(1);
        }
        if (tau_[0].size() != m_z.size()){
            std::cerr << AT<< " size mismatch tau[0].size()="<<tau_[0].size()<<" != z.size()="<<m_z.size()<<"\n";
            exit(1);
        }

        size_t ii = 0;
        for (size_t i = 0; i < m_freq.size(); i++){
            for (size_t j = 0; j < m_z.size(); j++){
//                m_freq[ii] = freq_[i];
//                m_z[ii] = z_[j];
                m_tau[ii] = tau_[i][j];
                ii++;
            }
        }
    }
    double interpolate(double i_freq, double i_z){
        /// check if loaded
        if ((m_freq.size() < 1) || (m_z.size() < 1) || (m_tau.size() < 1))
            load_h5_table();

        /// check if in bounds
        if (i_freq < m_freq[0]) return 0;
        if (i_freq > m_freq[m_freq.size()-1]) i_freq = m_freq[m_freq.size()-1];
        if (i_z < m_z[0]) return 0;
        if (i_z > m_z[m_z.size()-1]) i_z = m_z[m_z.size()-1];
        /// interpolate EBL table
        Interp2d interp2D = Interp2d( m_z,m_freq, m_tau);
        double tau = interp2D.InterpolateBilinear(i_z,i_freq);
        return tau;
    }

    void test_cases(){
        std::cout << interpolate(4.17852721e+25,0.5)<< " Expected 0.62492 "<< "\n";
        std::cout << interpolate(4.46885518e+27,1.) << " Expected 414.81 "<< "\n";
        std::cout << interpolate(4.46885518e+27,2.) << " Expected 1776.4800000000002 "<< "\n";
        exit(0);
    }
};

/// look-up tables that are general for all blast waves (pass by a reference)
struct CommonTables{
    std::unique_ptr<logger> p_log;
    CommonTables(){
//        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "CommonTables");
    }

    /// arrays
    Vector comov_gam_grid_fs{};
    Vector comov_gam_grid_rs{};
    Vector comov_freq_grid_fs{};
    Vector comov_freq_grid_rs{};

    /// kernels
    SynKernel synKernel_fs{};
    SSCKernel sscKernel_fs{};

    SynKernel synKernel_rs{};
    SSCKernel sscKernel_rs{};

    EBL ebl{};

    // todo nuclear heating table ...

    // magnetar table ...
};


#endif //SRC_COMMON_H
