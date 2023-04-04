//
// Created by vsevolod on 04/04/23.
//

#include "utils.h"
#include "logger.h"
#include "pch.h"
#include "H5Easy.h"

#ifndef SRC_INITIAL_DATA_H
#define SRC_INITIAL_DATA_H



class EjectaID{
    Vector m_mom{};
    Vector m_theta0{};
    Vector m_ctheta0{};
    VecVector m_ek_corr{};
    VecVector m_mass_corr{};
    VecVector m_ye_corr{};
    VecVector m_s_corr{};
    Vector m_ek_hist{};
    Vector m_mass_hist{};
    Vector m_s_hist{};
    Vector m_ye_hist{};
    std::unique_ptr<logger> p_log = nullptr;
public:
    enum IDTYPE { i_id_hist, i_id_corr };
    enum STUCT_TYPE { iadaptive, ipiecewise };
//    enum Q {iE, iMom, iMass, iYe, iS, };
private:
    struct Pars{

    };
    std::unique_ptr<Pars> p_pars = nullptr;
public:
    IDTYPE idtype{};
    STUCT_TYPE stuctType{};
    EjectaID(std::string path_to_table, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EjectaID");
        p_pars = std::make_unique<Pars>();
//        auto path_to_table = pars.m_path_to_ejecta_id;
//        path_to_table = "../../tst/dynamics/corr_vel_inf_theta.h5";
        if (!std::experimental::filesystem::exists(path_to_table))
            throw std::runtime_error("File not found. " + path_to_table);

        LoadH5 ldata;
        ldata.setFileName(path_to_table);
        ldata.setVarName("mom");
        m_mom = ldata.getData();
        int beta_size = ldata.getSize();

        ldata.setVarName("ctheta");
        m_ctheta0 = ldata.getData();

        ldata.setVarName("theta");
        m_theta0 = ldata.getData();

        if (m_theta0.size()!=m_ctheta0.size()){
            (*p_log)(LOG_ERR,AT) << " size mismatch. ctheta+1="<<m_ctheta0.size()+1<<" theta="<<m_theta0.size()<<"\n";
            exit(1);
        }

        ldata.setVarName("ek");
        int ek_size = ldata.getSize();
        if (ek_size!=beta_size) {
            m_ek_corr = ldata.getData2Ddouble();
            m_ek_hist = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 2D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", ek="<<m_ek_corr.size()<< "x" << m_ek_corr[0].size()<<"] \n";
            idtype = i_id_corr;
        }
        else{
            m_ek_hist = ldata.getDataVDouble();
            m_ek_corr = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 1D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", ek="<<m_ek_hist.size()<<"] \n";
            idtype = i_id_hist;
        }

        ldata.setVarName("mass");
        int mass_size = ldata.getSize();
        if (ek_size!=beta_size) {
            m_mass_corr = ldata.getData2Ddouble();
            m_mass_hist = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 2D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", mass="<<m_mass_corr.size()<< "x" << m_mass_corr[0].size()<<"] \n";
            idtype = i_id_corr;
        }
        else{
            m_mass_hist = ldata.getDataVDouble();
            m_mass_corr = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 1D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", mass="<<m_mass_hist.size()<<"] \n";
            idtype = i_id_hist;
        }

        ldata.setVarName("ye");
        int ye_size = ldata.getSize();
        if (ye_size!=beta_size) {
            m_ye_corr = ldata.getData2Ddouble();
            m_ye_hist = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 2D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", ye="<<m_ye_corr.size()<< "x" << m_ye_corr[0].size()<<"] \n";
            idtype = i_id_corr;
        }
        else{
            m_ye_hist = ldata.getDataVDouble();
            m_ye_corr = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 1D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", ye="<<m_ye_hist.size()<<"] \n";
            idtype = i_id_hist;
        }

        ldata.setVarName("s");
        int s_size = ldata.getSize();
        if (s_size!=beta_size) {
            m_s_corr = ldata.getData2Ddouble();
            m_s_hist = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 2D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", s="<<m_s_corr.size()<< "x" << m_s_corr[0].size()<<"] \n";
            idtype = i_id_corr;
        }
        else{
            m_s_hist = ldata.getDataVDouble();
            m_s_corr = {};
            (*p_log)(LOG_INFO, AT)
                    << "h5table 1D loaded [mom="<<m_mom.size()<<", m_ctheta0="
                    << m_ctheta0.size()<<", s="<<m_s_hist.size()<<"] \n";
            idtype = i_id_hist;
        }

        (*p_log)(LOG_INFO, AT) << "Ejecta ID loaded\n";

    }
    void setPars(StrDbMap & pars, StrStrMap & opts){
        /// method for EATS (piece-wise (AKA G. Lamb et al) or Adaptive (AKA vanEarten et al)
        std::string method = getStrOpt("eats_method",opts,AT,p_log,"",true);
        STUCT_TYPE method_eats;
        if(method == "piece-wise")
            method_eats = STUCT_TYPE::ipiecewise;
        else if(method == "adaptive")
            method_eats = STUCT_TYPE::iadaptive;
        else{
            std::cerr << " option for: " << "eats_method"
                      <<" given: " << method
                      << " is not recognized \n"
                      << "Possible options: "
                      << " piece-wise " << " adaptive " << "\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        stuctType = method_eats;
    }
    Vector & getMome(){ return m_mom; }
    Vector & getTheta0(){ return m_theta0; }
    Vector & getCtheta0(){ return m_ctheta0; }
    VecVector & getEkCorr(){ return m_ek_corr; }
    VecVector & getMassCorr(){ return m_mass_corr; }
    VecVector & getYeCorr(){ return m_ye_corr; }
    VecVector & getEntropyCorr(){ return m_s_corr; }
    Vector & getEkHist(){ return m_ek_hist; }
    Vector & getMassHist(){ return m_mass_hist; }
    Vector & getEntropyHist(){ return m_s_hist; }
    Vector & getYeHist(){ return m_ye_hist; }
    /// ---------------------------------------------
};


#endif //SRC_INITIAL_DATA_H
