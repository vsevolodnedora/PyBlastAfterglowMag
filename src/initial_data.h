//
// Created by vsevolod on 04/04/23.
//

#include "utilitites/utils.h"
#include "utilitites/logger.h"
#include "utilitites/pch.h"
#include "utilitites/H5Easy.h"

#ifndef SRC_INITIAL_DATA_H
#define SRC_INITIAL_DATA_H

class EjectaID2{
    enum IDTYPE { i_id_hist, i_id_corr };
    struct Pars{};
    std::vector< // v_ns
            std::vector< // shells
                    std::vector<double>>> // layers
    m_data{};
    std::unique_ptr<logger> p_log = nullptr;
    std::unique_ptr<Pars> p_pars = nullptr;
    unsigned m_loglevel{};
public:
    enum Q{ ir,imom,itheta,ictheta,iek,imass,iye,is,
            itheta_c_l, itheta_c_h, itheta_c, ieint    };
    std::vector<std::string> m_names{
            "r","mom","theta","ctheta","ek","mass","ye","s",
            "theta_c_l","theta_c_h","theta_c", "eint"
    };
    std::vector<std::string> m_v_ns {
            "r","mom", "theta", "ctheta", "ek", "mass", "ye", "s"
    };
    enum STUCT_TYPE { iadaptive, ipiecewise };
    IDTYPE idtype{};  STUCT_TYPE method_eats{};
    size_t nshells=0, nlayers=0, ncells=0;
    double theta_core{}; double theta_wing{};
    double theta_max{};
    EjectaID2(std::string path_to_table,
              std::string eats_method,
              bool use_1d_id,
              bool load_r0, double t0,
              unsigned loglevel){
        m_loglevel=loglevel;
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EjectaID");
        p_pars = std::make_unique<Pars>();
        m_data.resize(m_names.size());
        /// --------------------------------
        if(eats_method == "piece-wise")
            method_eats = STUCT_TYPE::ipiecewise;
        else if(eats_method == "adaptive")
            method_eats = STUCT_TYPE::iadaptive;
        else{
            std::cerr << " option for: " << "eats_method"
                      <<" given: " << eats_method
                      << " is not recognized \n"
                      << "Possible options: "
                      << " piece-wise " << " adaptive " << "\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        _load_id_file(path_to_table,use_1d_id, load_r0, t0);
        theta_max = CGS::pi/2.; // for dynamics
    }
    double get(size_t ish, size_t il, Q iv){
        return m_data[iv][ish][il];
    }
    inline static double ctheta(double theta, size_t ilayer, size_t nlayers_pw){
        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        if (theta > p_pars->theta_max ){
//            std::cerr << AT << " theta="<<theta<<" > theta_max=" << p_pars->theta_max << "\n";
//        }
//        if (std::fabs( theta - p_pars->theta_b0) > 1e-2){
//            std::cerr << AT << " theta="<<theta<<" < theta_b0=" << p_pars->theta_b0 <<"\n";
//            exit(1);
//        }
//        double ctheta = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w); // TODO WROOOONG

        double ctheta = 0.;
        if (ilayer > 0) {
            //
            double fac0 = (double)ilayer/(double)nlayers_pw;
            double fac1 = (double)(ilayer+1)/(double)nlayers_pw;
//            std::cout << std::asin(CGS::pi*3/4.) << "\n";
            if (!std::isfinite(std::sin(theta))){
                std::cerr << " sin(theta= "<<theta<<") is not finite... Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }

            double x2 = fac1*std::sin(theta / 2.);
            double xx2 = 2.*std::asin(x2);
            double x1 = fac0*std::sin(theta / 2.);
            double xx1 = 2.*std::asin(x1);

            ctheta = 0.5 * (xx1 + xx2);
            if (!std::isfinite(ctheta)){
                std::cerr << "ctheta is not finite. ctheta="<<ctheta<<" Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        return ctheta;
    }
    /// phi grid for a given layer (given number of phi cells)
    static Vector getCphiGridPW( size_t ilayer ){
        size_t cil = CellsInLayer(ilayer);
        Vector cphis ( cil );
        for (size_t j = 0; j < cil; j++){
            cphis[j] = (double)j * 2.0 * CGS::pi / (double)cil;
        }
        return std::move ( cphis );
    }
//    static Vector getCphiGridA( size_t ilayer){
//        getCphiGridA
//    }
    static void getCphiGridPW( Vector & cphis, size_t ilayer ){
        size_t cil = CellsInLayer(ilayer);
//        Vector cphis ( cil );
        for (size_t j = 0; j < cil; j++){
            cphis[j] = (double)j * 2.0 * CGS::pi / (double)cil;
        }
    }
    static size_t CellsInLayer(const size_t &i_layer){
        return 2 * i_layer + 1;
    }


    /// initial grid for [a] EATS method
    static void _init_a_grid(Vector & theta_c_l, Vector & theta_c_h, Vector & theta_c, size_t nlayers, double theta_w){
        theta_c_l.resize(nlayers,0.);
        theta_c_h.resize(nlayers,0.);
        theta_c.resize(nlayers,0.);
        double dtheta = theta_w / (double) nlayers;
        for (size_t i = 0; i < nlayers; i++) {
            /// account for geometry
            double theta_c_i   = (double) i * dtheta + dtheta / 2.;
            double i_theta_c_l = (double) i * dtheta;
            double i_theta_c_h = (double) (i + 1) * dtheta;
            theta_c[i]   = theta_c_i;//thetas_c[i] = theta_c_i ;
            theta_c_l[i] = i_theta_c_l;//thetas_c_l[i] = i_theta_c_l ;
            theta_c_h[i] = i_theta_c_h;//thetas_c_h[i] = i_theta_c_h ;
        }
    }

    static void _evalCellsInLayer(size_t nlayers, std::vector<size_t> & cil){
        cil.resize( nlayers );
        for (size_t i = 0; i < nlayers; i++)
            cil[i] = CellsInLayer(i);
    }

    static size_t _evalTotalNcells(size_t nlayers){
        /// evaluateShycnhrotronSpectrum the number of phi cells in each 'theta' layer
        std::vector<size_t> cil;
        _evalCellsInLayer(nlayers, cil);
        size_t ncells = std::accumulate(cil.begin(), cil.end(), decltype(cil)::value_type(0));
        return ncells;
    }
private:

    //    void _init_a_grid(size_t ish, size_t nlayers, double theta_w){
//        m_data[Q::itheta_c][ish].resize(nlayers,0.);
//        m_data[Q::itheta_c_h][ish].resize(nlayers,0.);
//        m_data[Q::itheta_c_l][ish].resize(nlayers,0.);
//        double dtheta = theta_w / (double) nlayers;
//        for (size_t i = 0; i < nlayers; i++) {
//            /// account for geometry
//            double theta_c_i   = (double) i * dtheta + dtheta / 2.;
//            double i_theta_c_l = (double) i * dtheta;
//            double i_theta_c_h = (double) (i + 1) * dtheta;
//            m_data[Q::itheta_c][ish][i]   = theta_c_i;//thetas_c[i] = theta_c_i ;
//            m_data[Q::itheta_c_l][ish][i] = i_theta_c_l;//thetas_c_l[i] = i_theta_c_l ;
//            m_data[Q::itheta_c_h][ish][i] = i_theta_c_h;//thetas_c_h[i] = i_theta_c_h ;
//        }
//    }

    /// initial grid for [a] EATS method
    void _init_a_grid(size_t ish, size_t nlayers, double theta_w) {
        _init_a_grid(m_data[Q::itheta_c_l][ish], m_data[Q::itheta_c_h][ish], m_data[Q::itheta_c][ish],
                     nlayers, theta_w);
    }

    void _init_pw_grid(size_t ish, size_t nlayers, double theta_w){
        m_data[ictheta][ish].resize(nlayers,0.);
        Vector theta_pw ( nlayers + 1 );
//        cthetas0.resizeEachImage( nlayers_pw );
        for (size_t i = 0; i < nlayers + 1; i++){
            double fac = (double)i / (double)nlayers;
            theta_pw[i] = 2.0 * std::asin( fac * std::sin(theta_w / 2.0 ) );

        }

        Vector thetas_h0_pw (nlayers );
        for (size_t i = 0; i < nlayers; ++i){
            m_data[ictheta][ish][i] = 0.5 * ( theta_pw[i+1] + theta_pw[i] );
            thetas_h0_pw[i] = theta_pw[i + 1];
            /// for tau_along_los # TODO replace the above 'thetas_h0_pw' with these two and  make sure they are correct
            m_data[Q::itheta_c_l][ish][i] = theta_pw[i];//thetas_c_l[i] = i_theta_c_l ;
            m_data[Q::itheta_c_h][ish][i] = theta_pw[i+1];//thetas_c_h[i] = i_theta_c_h ;
        }

        /// eval the number of phi cells in each 'theta' layer
        ncells = _evalTotalNcells(nlayers);
    }
    void _load_id_file(std::string & path_to_table, bool & use_1d_id, bool loadr0, double t0){
        if (!std::experimental::filesystem::exists(path_to_table))
            throw std::runtime_error("File not found. " + path_to_table);
        /// ---------------------------------------------
        LoadH5 ldata;
        ldata.setFileName(path_to_table);

//        size_t nshells=0, m_nlayers=0;
        if (use_1d_id){
            /// expect 1 shell with varying properties
            nshells = 1;
            for (auto & arr : m_data) arr.resize(m_names.size());
            for (size_t ish = 0; ish < nshells; ish++) {
                for (size_t i_v_n = 0; i_v_n < m_v_ns.size(); i_v_n++) {
                    ldata.setVarName(m_v_ns[i_v_n]);
                    m_data[i_v_n][ish] = std::move ( ldata.getDataVDouble() ); // Mode data
                }
            }
            nlayers = m_data[0][0].size();
        }
        else{
            /// expect N shells with varying properties [nshells; m_nlayers]
            ldata.setVarName("mass");
            VecVector tmp = ldata.getData2Ddouble();
            nshells = tmp.size();
            nlayers = tmp[0].size();
            for (auto & arr : m_data) arr.resize(nshells);
            for (size_t ish = 0; ish < nshells; ish++) {
                for (size_t i_v_n = 0; i_v_n < m_v_ns.size(); i_v_n++) {
                    ldata.setVarName(m_v_ns[i_v_n]);
                    m_data[i_v_n] = std::move( ldata.getData2Ddouble() ); //
                }
            }
        }
        /// ---------------------------
        (*p_log)(LOG_INFO,AT) << "Initial data loaded with nshells="<<nshells<<" m_nlayers="<<nlayers<<"\n";
        /// ---------------------------
        for (size_t ish = 0; ish < nshells; ish++){
            (*p_log)(LOG_ERR,AT)<<" theta_core is NOT given in new ID. Fix it by evaluating it FROM profile!\n";
            theta_wing = m_data[Q::itheta][ish][nlayers-1];
            theta_core = m_data[Q::itheta][ish][nlayers-1];
            _init_a_grid(ish, nlayers, theta_wing);
            _init_pw_grid(ish, nlayers, theta_wing);
        }
        /// ---------------------------
        (*p_log)(LOG_INFO,AT) << "Angular grids are initialized. nshells="<<nshells<<" m_nlayers="<<nlayers<<"\n";
        /// ---------------------------
        switch (method_eats) {
            case iadaptive:
                for (size_t ish = 0; ish < nshells; ish++) {
                    for (size_t i = 0; i < nlayers; i++) {
                        double frac_of_solid_ang = 2
                                                   * std::sin(0.5 * m_data[Q::itheta_c_h][ish][i])
                                                   * std::sin(0.5 * m_data[Q::itheta_c_h][ish][i]);
                        m_data[Q::iek][ish][i] *= (frac_of_solid_ang / 2.);// TODO remove it and include into ID maker
                        m_data[Q::imass][ish][i] *= (frac_of_solid_ang / 2.);
                    }
                }
                break ;
            case ipiecewise:
                for (size_t ish = 0; ish < nshells; ish++) {
                    for (size_t i = 0; i < nlayers; i++) {
                        m_data[Q::iek][ish][i] /= (double) CellsInLayer(i); // TODO remove it and make ID that includes it
                        m_data[Q::imass][ish][i] /= (double) CellsInLayer(i);
                    }
                }
                break ;
        }
        (*p_log)(LOG_INFO,AT) << "Energy and mass are rescaled."<<"\n";
        /// ---------------------------
        if (loadr0){
            /// pass

        }
        else{
            /// pass
            for (size_t ish = 0; ish < nshells; ish++) {
                for (size_t i = 0; i < nlayers; i++) {
                    m_data[Q::ir][ish][i] = BetFromMom(m_data[Q::imom][ish][i]) * CGS::c * t0;
                }
            }
        }
        /// ---------------------------
//        ldata.setVarName(m_v_ns[i_v_n]);
//        m_data[i_v_n][ish] = ldata.getDataVDouble();

    }

};

#if 0
class EjectaID {
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
    struct PWNPars{

    };
    std::unique_ptr<PWNPars> p_pars = nullptr;
public:
    IDTYPE idtype{};
    STUCT_TYPE stuctType{};
    EjectaID(std::string path_to_table, int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EjectaID");
        p_pars = std::make_unique<PWNPars>();
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
            (*p_log)(LOG_ERR,AT) << " size mismatch. ctheta="<<m_ctheta0.size()<<" theta="<<m_theta0.size()<<"\n";
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


struct LatStruct{

private:
//    logger * p_log;
    std::unique_ptr<logger> p_log;

private:

    /**
     * Computes the number of 'phi' cell in each 'theta' layer
     *
     * @param i_layer index of the 'theta' layer
     * @return
     */

    void setThetaGridPW(){
        theta_pw.resize( nlayers_pw + 1 );
        cthetas0.resize( nlayers_pw );
        for (size_t i = 0; i < nlayers_pw + 1; i++){
            double fac = (double)i / (double)nlayers_pw;
            theta_pw[i] = 2.0 * std::asin( fac * std::sin(m_theta_w / 2.0 ) );
        }

        thetas_h0_pw.resize(nlayers_pw );
        for (size_t i = 0; i < nlayers_pw; ++i){
            cthetas0[i] = 0.5 * ( theta_pw[i+1] + theta_pw[i] );
            thetas_h0_pw[i] = theta_pw[i + 1];
        }
//        ncells_pw = 0;

        /// compute the number of phi cells in each 'theta' layer
        cil.resize( nlayers_pw );
        for (size_t i = 0; i < nlayers_pw; i++)
            cil[i] = CellsInLayer(i);
        ncells = cil.sum(); /// total number of cells

//        std::cout << thetas << "\n";
//        std::cout << cthetas0 << "\n";
//        exit(1);


//        cil.resize( nlayers_pw );
//        for (size_t i = 0; i < nlayers_pw; i++){ cil[i] = CellsInLayer(i); }
//        for (size_t i = 1; i < nlayers_pw; i++){ cil[i] += cil[i-1]; }
//        ncells = cil.sum();

    }

    void setThetaGridA() {

        thetas_c_l.resize( nlayers_a );
        thetas_c_h.resize( nlayers_a );
        thetas_c.resize( nlayers_a );

        double dtheta = m_theta_w / (double) nlayers_a;
//        double _tmp=0;
        for (size_t i = 0; i < nlayers_a; i++) {

            /// account for geometry
            double theta_c_i = (double)i * dtheta + dtheta / 2.;
            double i_theta_c_l = (double) i * dtheta;
            double i_theta_c_h = (double) (i + 1) * dtheta;
//            double i_theta_h = i_theta_c_h;

            thetas_c[i] = theta_c_i ;
            thetas_c_l[i] = i_theta_c_l ;
            thetas_c_h[i] = i_theta_c_h ;
//            std::cout << "ilayer=" << i << "theta"
//            _tmp+=i_theta_c_h;
        }
//        std::cout << "theta_w="<<m_theta_w<<" thetas_c_h[-1]=" << thetas_c_h[thetas_c_h.m_size()-1] << "\n";
//        exit(0);
    }

public:

    enum METHOD_eats{i_pw, i_adap};
    static METHOD_eats setEatsMethod(std::string method){
        /// method for EATS (piece-wise (AKA G. Lamb et al) or Adaptive (AKA vanEarten et al)
        LatStruct::METHOD_eats method_eats;
        if(method == "piece-wise")
            method_eats = LatStruct::METHOD_eats::i_pw;
        else if(method == "adaptive")
            method_eats = LatStruct::METHOD_eats::i_adap;
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
        return method_eats;
    }

    enum METHODS { iUniform, iGaussian, iCustom };

    void initUniform(double E_iso, double Mom0, double theta_h, double M0, double Ye, double s,
                     size_t n_layers, std::string eats_method){
        method = iUniform;
        m_theta_w = theta_h;
        m_theta_c = theta_h;

        m_method_eats = setEatsMethod(eats_method);
        if ((m_method_eats==LatStruct::i_adap) && (n_layers > 1)){
            std::cerr<< " When using 'adaptive' eats and 'uniform' structure, "
                        "the n_layers should 1, Given={}"<<n_layers<<"\n";
            std::cerr << AT << "\n";
            exit(1);
        }

        // set piece-wise
        nlayers_pw = n_layers;
        dist_E0_pw.resize( nlayers_pw );
        dist_Mom0_pw.resizeEachImage( nlayers_pw );
        dist_M0_pw.resize( nlayers_pw );
        dist_Ye_pw.resize( nlayers_pw );
        dist_s_pw.resize( nlayers_pw );
        setThetaGridPW();
        if (nlayers_pw  < 1){
            std::cerr << " grid not initialized.\n Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
        double one_min_cos = 2 * std::sin(0.5 * theta_h) * std::sin(0.5 * theta_h); // for pi/2 -> 1.
        double ang_size_layer = 2.0 * CGS::pi * one_min_cos;
        for (size_t i = 0; i < cthetas0.size(); i++){
            dist_E0_pw[i] = E_iso * ang_size_layer / (4.0 * CGS::pi);
            dist_Mom0_pw[i] = Mom0;
            double Gamma = EQS::GamFromMom(dist_Mom0_pw[i]);
            dist_M0_pw[i] = M0 < 0 ? dist_E0_pw[i] / (Gamma * CGS::c * CGS::c) : M0 * ang_size_layer / (4 * CGS::pi);

            dist_E0_pw[i] /= (double)ncells;
            dist_M0_pw[i] /= (double)ncells;
            dist_Ye_pw[i] = Ye;
            dist_s_pw[i] = s;
        }

        // set adaptive
        nlayers_a = n_layers;
        setThetaGridA();
        dist_E0_a.resize( nlayers_a );
        dist_Mom0_a.resize( nlayers_a );
        dist_M0_a.resize( nlayers_a );
        dist_Ye_a.resize( nlayers_a );
        dist_s_a.resize( nlayers_a );
        double frac_of_solid_ang = 2 * sin(0.5 * theta_h) * sin(0.5 * theta_h); // for pi/2 -> 1.
        dist_E0_a[0] = E_iso * frac_of_solid_ang / 2.;
        dist_Mom0_a[0] = Mom0;
        double Gamma = EQS::GamFromMom(dist_Mom0_a[0]);
        dist_M0_a[0] = M0 < 0 ? E_iso / ((Gamma - 1.0) * CGS::c*CGS::c) * frac_of_solid_ang / 2. : M0 * frac_of_solid_ang / 2.;
        dist_Ye_a[0] = 0.;
        dist_s_a[0] = 0.;
//        std::cout<<m_theta_w<<"\n";
//        exit(1);

        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

    }

    void initGaussian( double E_iso_c, double Gamma0c, double theta_c, double theta_w, double M0c,
                       size_t n_layers_pw, size_t n_layers_a, bool gflat=false,
                       std::string eats_method="piece-wise"){

        method = iGaussian;
        m_method_eats = setEatsMethod(eats_method);
        m_theta_w = theta_w;
        m_theta_c = theta_c;

        // set piece-wise
        nlayers_pw = n_layers_pw;
        dist_E0_pw.resize( nlayers_pw );
        dist_Mom0_pw.resize( nlayers_pw );
        dist_M0_pw.resize( nlayers_pw );
        dist_Ye_pw.resize( nlayers_pw, 0. );
        dist_s_pw.resize( nlayers_pw, 0. );
        setThetaGridPW();
        double ang_size_layer = 2.0 * CGS::pi * ( 2.0 * std::sin(0.5 * theta_w) * std::sin(0.5 * theta_w) );
        for (size_t i = 0; i < cthetas0.size(); i++){
            dist_E0_pw[i] = E_iso_c * ang_size_layer / (4.0 * CGS::pi) * std::exp( -1. * cthetas0[i] * cthetas0[i] / (theta_c * theta_c) );
            double Gamma;
            if (gflat)
                Gamma = EQS::MomFromGam( Gamma0c );
            else {
                Gamma = 1. + (Gamma0c - 1.) * std::exp(-1. * cthetas0[i] * cthetas0[i] / (2. * theta_c * theta_c));
            }
            dist_Mom0_pw[i] = EQS::MomFromGam(Gamma);
            dist_M0_pw[i] = dist_E0_pw[i] / (Gamma * CGS::c * CGS::c);
            dist_E0_pw[i] /= (double)ncells;
            dist_M0_pw[i] /= (double)ncells;
        }

        // set adaptive
        nlayers_a = n_layers_a;
        setThetaGridA();
        dist_E0_a.resize( nlayers_a );
        dist_Mom0_a.resize( nlayers_a );
        dist_M0_a.resize( nlayers_a );
        dist_Ye_a.resize( nlayers_a, 0. );
        dist_s_a.resize( nlayers_a, 0. );
        for (size_t i = 0; i < nlayers_a; i++) {
            double frac_of_solid_ang = 2 * std::sin(0.5 * thetas_c_h[i]) * std::sin(0.5 * thetas_c_h[i]);
            dist_E0_a[i] = E_iso_c * std::exp(-0.5 * ( thetas_c[i] * thetas_c[i] / theta_c / theta_c ) );
            double Gamma;
            if (gflat)
                Gamma = Gamma0c;
            else
                Gamma = 1.0 + (Gamma0c - 1) * dist_E0_a[i] / E_iso_c;
            dist_Mom0_a[i] = EQS::MomFromGam( Gamma );
            dist_M0_a[i] = dist_E0_a[i] / (( Gamma - 1.0) * CGS::c*CGS::c );

            dist_E0_a[i] *= ( frac_of_solid_ang / 2. );
            dist_M0_a[i] *= ( frac_of_solid_ang / 2. );
        }

        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
//        int a = 1;

        std::cout << dist_E0_pw << "\n";
        std::cout << dist_Mom0_pw << "\n";
        std::cout << dist_M0_pw << "\n";
        std::cout << cthetas0 << "\n";
        exit(1);
    }
//    ~LatStruct(){ delete p_log; }

//    static std::vector<std::string> listParsUniformStruct() {
//        return {"Eiso", "Gamma0", "theta_h", "M0", "nlayers_pw"};}
    static std::vector<std::string> listParametersAnalyticBlastWave() {
        return {"Eiso_c", "Gamma0c", "theta_c", "theta_w", "M0c", "nlayers_pw", "nlayers_a"};
    }
//    static std::vector<std::string> listParsCustomStruct() {
//        return {"dist_thetas", "dist_EEs", "dist_Gam0s", "dist_MM0s"}; }
//    static std::vector<std::string> listStructOpts() {
//        return {"type", "gflat"}; }

    void initAnalytic(StrDbMap & pars, StrStrMap & opts, std::string method_eats,  unsigned loglevel){
//        p_log = new logger(std::cout, CurrLogLevel, "LatStruct");
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");
//        p_log->set_log_level(loglevel);
//        size_t n_layers = 0;
        nlayers_pw = (size_t)getDoublePar("nlayers_pw", pars, AT, p_log, 0, true);
        nlayers_a = (size_t)getDoublePar("nlayers_a", pars, AT, p_log, 0, true);
        (*p_log)(LOG_INFO,AT) << "Initializing analytic lateral structure ["
                              << method_eats<<", nlpw="<<nlayers_pw<<", nla="<<nlayers_a<<"] \n";
//        std::string opt = "eats_method";
//        if ( opts.find(opt) == opts.end() ){
//            std::cerr << " option for: " << opt
//                      << " is not given \n"
//                      << "Possible options: "
//                      << " piece-wise " << " adaptive " << "\n"
//                      << " Exiting...";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
        m_method_eats = setEatsMethod(method_eats);
        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

//        n_layers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;


        if (opts.find("type")==opts.end()){
            (*p_log)(LOG_ERR,AT) << " type of the lateral structure is not set. " << "Exiting...";
            exit(1);
        }
        if (opts.at("type") == "uniform"){
            // uniform structure
            for (auto & v_n : listParametersAnalyticBlastWave()){
                if (pars.find(v_n) == pars.end()){
                    (*p_log)(LOG_ERR,AT) << " not given parameter:" << v_n << " that is required for "
                                         << " uniform " << " structures. " << " Exiting...\n";
                    exit(1);
                }
            }
            double _theta_c = getDoublePar("theta_c", pars, AT, p_log, 0, true);
            double _theta_w = getDoublePar("theta_w", pars, AT, p_log, 0, true);
            if (nlayers_a!=1){
                (*p_log)(LOG_ERR, AT) << " nlayers_a must be 1 for uniform blastwave Given = "
                                      << nlayers_a <<"\n";
                exit(1);
            }
            if (_theta_c!=_theta_w){
                (*p_log)(LOG_ERR, AT) << " theta_c must be equal theta_w for uniform blast wave. Given: "
                                      << "theta_c = " << _theta_c << " theta_w = "<<_theta_w<<"\n";
                exit(1);
            }
            initUniform( getDoublePar("Eiso_c", pars, AT, p_log, 0, true),//(double)pars.at("Eiso"),
                         getDoublePar("Gamma0c", pars, AT, p_log, 0, true),//(double)pars.at("Gamma0"),
                         getDoublePar("theta_w", pars, AT, p_log, 0, true),//(double)pars.at("theta_h"),
                         getDoublePar("M0c", pars, AT, p_log, -1, false),//(double)pars.at("M0"),
                         0.,
                         0.,
                         m_nlayers,
                         method_eats//pars.at("nlayers_pw")
            );
        }
        else if (opts.at("type") == "gaussian"){
            // gaussian structure
            for (auto & v_n : listParametersAnalyticBlastWave()){
                if (pars.find(v_n) == pars.end()){
                    std::cerr << " not given parameter:" << v_n << " that is required for "
                              << " gaussian " << " structures " << "\n"
                              << " Exiting...\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }
            }
            std::string par = "gflat"; auto val = false;
            if (opts.find(par)!=opts.end()){
                if (opts.at(par) == "yes"){
                    val = true;
                }
                else if (opts.at(par) == "no"){
                    val = false;
                }
                else {
                    std::cerr << " option for: " << par
                              <<" given: " << opts.at(par)
                              << " is not recognized \n"
                              << "Possible options: "
                              << " yes " << " no " << "\n";
                    std::cerr << AT << "\n";
                    exit(1);
                }

            }
            initGaussian( getDoublePar("Eiso_c", pars, AT, p_log, -1., true),//pars.at("Eiso_c"),
                          getDoublePar("Gamma0c", pars, AT, p_log, -1., true),//pars.at("Gamma0c"),
                          getDoublePar("theta_c", pars, AT, p_log, -1., true),//pars.at("theta_c"),
                          getDoublePar("theta_w", pars, AT, p_log, -1., true),//pars.at("theta_w"),
                          getDoublePar("M0c", pars, AT, p_log, -1., false),//(double)pars.at("M0c"),
                          nlayers_pw,//pars.at("nlayers_pw"),
                          nlayers_a,//pars.at("nlayers_a"),
                          val,
                          method_eats
            );
        }
        else{
            std::cerr << " lateral structure type: " << opts.at("type")
                      << " is not recognized. \n"
                      << " Available options: "
                      << " uniform " << " gaussian " << "\n"
                      << " Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }

//        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
        (*p_log)(LOG_INFO,AT) << " setting " << opts.at("type") << " lateral structure\n";
    }

    void initCustom(Vector & dist_thetas, Vector & dist_cthetas, Vector & dist_EEs,
                    Vector & dist_Ye, Vector & dist_s, Vector & dist_Mom0s, Vector & dist_MM0s,
                    bool force_grid, std::string eats_method, unsigned loglevel){

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");

        m_method_eats = setEatsMethod(eats_method);
//        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

        if (dist_cthetas.empty() || dist_EEs.empty() || dist_Mom0s.empty()
            || dist_MM0s.empty() || dist_Ye.empty() || dist_s.empty()){
            std::cerr << " One of the input arrays is empty: "
                      << "dist_cthetas(" << dist_cthetas.size()
                      << ") dist_EEs(" << dist_EEs.size()
                      << ") dist_Yes(" << dist_Ye.size()
                      << ") dist_Ss(" << dist_s.size()
                      << ") dist_Mom0s(" << dist_Mom0s.size()
                      //                                  << ") dist_Beta0s(" << dist_Beta0s.size()
                      << ") dist_MM0s(" << dist_MM0s.size() << ")\n"
                      << " Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        if ( (dist_cthetas.size() != dist_EEs.size())
             || (dist_Mom0s.size() != dist_EEs.size())
             || (dist_Ye.size() != dist_EEs.size())
             || (dist_s.size() != dist_EEs.size())
             || (dist_MM0s.size() != dist_EEs.size()) ){
            std::cerr << " Mismatch in input array size "
                      << "dist_cthetas(" << dist_cthetas.size()
                      << ") dist_EEs(" << dist_EEs.size()
                      << ") dist_Yes(" << dist_Ye.size()
                      << ") dist_Ss(" << dist_s.size()
                      << ") dist_Mom0s(" << dist_Mom0s.size()
                      << ") dist_MM0s(" << dist_MM0s.size() << ")\n"
                      << " Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        if (dist_cthetas[dist_cthetas.size()-1] > CGS::pi/2.){
            (*p_log)(LOG_ERR,AT) << "wrong initial data dist_thetas[-1]="
                                 << dist_cthetas[dist_cthetas.size()-1]<<" > pi/2.="<<CGS::pi/2.<<"\n";
            exit(1);
        }
        if (dist_thetas[dist_thetas.size()-1] > CGS::pi/2.){
            (*p_log)(LOG_ERR,AT) << "wrong initial data dist_thetas[-1]="
                                 << dist_thetas[dist_thetas.size()-1]<<" > pi/2.="<<CGS::pi/2.<<"\n";
            exit(1);
        }
        if (dist_cthetas[0] < 0.){
            (*p_log)(LOG_ERR,AT) << "wrong initial data dist_thetas[]="
                                 << dist_cthetas[dist_cthetas.size()-1]<<" < 0.\n";
            exit(1);
        }
        if (dist_thetas[0] < 0.){
            (*p_log)(LOG_ERR,AT) << "wrong initial data dist_thetas[]="
                                 << dist_thetas[dist_thetas.size()-1]<<" < 0.\n";
            exit(1);
        }

        if (m_method_eats == METHOD_eats::i_pw) {
            method = iCustom;
            m_theta_w = dist_thetas[dist_thetas.size() - 1];
            m_theta_c = dist_thetas[dist_thetas.size() - 1];
            (*p_log)(LOG_ERR, AT) << " m_theta_c is inaccuate\n";

            // set piece-wise
            nlayers_pw = dist_cthetas.size();
            m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
            dist_E0_pw.resize(m_nlayers);
            dist_Ye_pw.resize(m_nlayers);
            dist_s_pw.resize(m_nlayers);
            dist_Mom0_pw.resize(m_nlayers);
            dist_M0_pw.resize(m_nlayers);
            setThetaGridPW();
            if ((!force_grid) &&
                (cthetas0[0] != dist_cthetas[0] || cthetas0[m_nlayers - 1] != dist_cthetas[m_nlayers - 1])) {
                dist_E0_pw.resize(m_nlayers);
                dist_Ye_pw.resize(m_nlayers);
                dist_s_pw.resize(m_nlayers);
                dist_Mom0_pw.resize(m_nlayers);
                dist_M0_pw.resize(m_nlayers);

                std::cerr << "force_grid=" << force_grid << "\n";
                std::cerr << "dist_thetas.size()=" << dist_cthetas.size() << " cthetas0.size()=" << cthetas0.size()
                          << "\n";
                std::cerr << "dist_thetas = " << dist_cthetas << "\n";
                std::cerr << "cthetas0    = " << cthetas0 << "\n";
                (*p_log)(LOG_ERR, AT) << "angular grid mismatch: "
                                      << "cthetas0[0]=" << cthetas0[0] << " != dist_thetas[0]=" << dist_cthetas[0]
                                      << " OR cthetas[m_nlayers-1]=" << cthetas0[m_nlayers - 1]
                                      << " != dist_thetas[nlayers_pw-1]=" << dist_cthetas[m_nlayers - 1]
                                      << " [m_nlayers=" << m_nlayers << ", " << eats_method << "]\n";
//            std::cerr << "Angular grid mismatch. Interpolating. No. Exititng...";
//            std::cerr << AT << "\n";
                exit(1);

                Array arr_cthetas(dist_cthetas.data(), dist_cthetas.size());
                Array arr_EEs(dist_EEs.data(), dist_EEs.size());
                Array arr_Yes(dist_EEs.data(), dist_EEs.size());
                Array arr_Ss(dist_EEs.data(), dist_EEs.size());
                Array arr_Mom0s(dist_Mom0s.data(), dist_Mom0s.size());
                Array arr_MM0s(dist_MM0s.data(), dist_MM0s.size());
                Array arr_cthetas0(cthetas0.data(), cthetas0.size());

                Interp1d::METHODS mth = Interp1d::iLagrangeLinear;
                Array int_E = Interp1d(arr_cthetas, arr_EEs).Interpolate(arr_cthetas0, mth);
                Array int_Ye = Interp1d(arr_cthetas, arr_Yes).Interpolate(arr_cthetas0, mth);
                Array int_s = Interp1d(arr_cthetas, arr_Ss).Interpolate(arr_cthetas0, mth);
                Array int_Mom = Interp1d(arr_cthetas, arr_Mom0s).Interpolate(arr_cthetas0, mth);
                Array int_M = Interp1d(arr_cthetas, arr_MM0s).Interpolate(arr_cthetas0, mth);

                dist_E0_pw = arrToVec(int_E);
                dist_Ye_pw = arrToVec(int_Ye);
                dist_s_pw = arrToVec(int_s);
                dist_Mom0_pw = arrToVec(int_Mom);
                dist_M0_pw = arrToVec(int_M);
            } else {
                dist_E0_pw = dist_EEs;
                dist_Mom0_pw = dist_Mom0s;
                dist_M0_pw = dist_MM0s;
                dist_Ye_pw = dist_Ye;
                dist_s_pw = dist_s;
            }

            for (size_t i = 0; i < m_nlayers; ++i) {
                dist_E0_pw[i] /= (double) CellsInLayer(i); // TODO remove it and make ID that includes it
                dist_M0_pw[i] /= (double) CellsInLayer(i);
            }

            /*
( 6.80844e+46, 6.75726e+46, 6.65604e+46, 6.50705e+46, 6.31357e+46, 6.07977e+46, 5.8106e+46, 5.51158e+46, 5.18861e+46, 4.84782e+46, 4.49531e+46, 4.13706e+46, 3.77868e+46, 3.42535e+46, 3.08165e+46, 2.75154e+46, 2.43827e+46, 2.14435e+46, 1.87163e+46, 1.62126e+46, 1.39377e+46, 1.18914e+46, 1.00688e+46, 8.46098e+45, 7.05607e+45, 5.83984e+45, 4.79659e+45, 3.90981e+45, 3.16277e+45, 2.53902e+45, 2.02277e+45, 1.59923e+45, 1.25473e+45, 9.76938e+44, 7.54841e+44, 5.78779e+44, 4.40389e+44, 3.32523e+44, 2.49154e+44, 1.85254e+44, 1.36685e+44, 1.00074e+44, 7.27058e+43, 5.24151e+43, 3.74956e+43, 2.66156e+43, 1.87466e+43, 1.31018e+43, 9.08577e+42, 6.25183e+42, )
( 299.857, 298.732, 296.494, 293.168, 288.791, 283.412, 277.09, 269.892, 261.895, 253.181, 243.839, 233.962, 223.643, 212.978, 202.062, 190.987, 179.845, 168.719, 157.691, 146.834, 136.215, 125.895, 115.925, 106.35, 97.2057, 88.5215, 80.3185, 72.6108, 65.4057, 58.7047, 52.5034, 46.7927, 41.5592, 36.7858, 32.4527, 28.5374, 25.0159, 21.8631, 19.0528, 16.5591, 14.3559, 12.4177, 10.7197, 9.23825, 7.95083, 6.83621, 5.87461, 5.0477, 4.33867, 3.73217, )
( 2.52632e+23, 2.51678e+23, 2.4978e+23, 2.46959e+23, 2.43247e+23, 2.38685e+23, 2.33322e+23, 2.27218e+23, 2.20435e+23, 2.13045e+23, 2.05122e+23, 1.96744e+23, 1.87992e+23, 1.78947e+23, 1.69689e+23, 1.60297e+23, 1.50846e+23, 1.41411e+23, 1.32058e+23, 1.2285e+23, 1.13844e+23, 1.05092e+23, 9.66366e+22, 8.85164e+22, 8.07619e+22, 7.33978e+22, 6.64418e+22, 5.99062e+22, 5.37971e+22, 4.81159e+22, 4.28588e+22, 3.80182e+22, 3.35828e+22, 2.95382e+22, 2.58677e+22, 2.25523e+22, 1.95718e+22, 1.6905e+22, 1.45301e+22, 1.24251e+22, 1.05681e+22, 8.93794e+21, 7.51387e+21, 6.27618e+21, 5.20617e+21, 4.2863e+21, 3.50025e+21, 2.83294e+21, 2.27052e+21, 1.80031e+21, )
( 0.00261053, 0.00783162, 0.0130528, 0.018274, 0.0234953, 0.0287168, 0.0339385, 0.0391605, 0.0443827, 0.0496052, 0.054828, 0.0600513, 0.0652749, 0.070499, 0.0757235, 0.0809486, 0.0861742, 0.0914004, 0.0966273, 0.101855, 0.107083, 0.112312, 0.117542, 0.122772, 0.128003, 0.133236, 0.138469, 0.143703, 0.148938, 0.154174, 0.159411, 0.164649, 0.169889, 0.175129, 0.180371, 0.185614, 0.190858, 0.196104, 0.201351, 0.206599, 0.211849, 0.2171, 0.222353, 0.227607, 0.232863, 0.238121, 0.24338, 0.248641, 0.253903, 0.259167, )

             */

//            std::cout << dist_E0_pw << "\n";
//            std::cout << dist_Mom0_pw << "\n";
//            std::cout << dist_M0_pw << "\n";
//            std::cout << cthetas0 << "\n";
//            exit(1);

        }
        else {
            // set adaptive
            // TODO IMPLEMENT ! But it might require re-inteprolation of angular dostributions as [a] and [pw] thetagrids differ
//        nlayers_a = n_layers_a;
            m_theta_w = dist_thetas[dist_thetas.size() - 1];
            m_theta_c = dist_thetas[dist_thetas.size() - 1];
            (*p_log)(LOG_ERR, AT) << " m_theta_c is inaccuate\n";
            nlayers_a = dist_cthetas.size();

            setThetaGridA();
            if ((!force_grid) &&
                (thetas_c[0] != dist_cthetas[0] || thetas_c[m_nlayers - 1] != dist_cthetas[m_nlayers - 1])) {
                std::cerr << "force_grid=" << force_grid << "\n";
                std::cerr << "thetas_c.size()=" << dist_cthetas.size() << " thetas_c.size()=" << cthetas0.size()
                          << "\n";
                std::cerr << "dist_thetas = " << dist_cthetas << "\n";
                std::cerr << "thetas_c    = " << thetas_c << "\n";
                (*p_log)(LOG_ERR, AT) << "angular grid mismatch: "
                                      << "thetas_c[0]=" << thetas_c[0] << " != dist_thetas[0]=" << dist_cthetas[0]
                                      << " OR thetas_c[m_nlayers-1]=" << thetas_c[m_nlayers - 1]
                                      << " != dist_thetas[nlayers_a-1]=" << dist_cthetas[m_nlayers - 1]
                                      << " [m_nlayers=" << m_nlayers << ", " << eats_method << "]\n";
//            std::cerr << "Angular grid mismatch. Interpolating. No. Exititng...";
//            std::cerr << AT << "\n";
                exit(1);
            }
            dist_E0_a = dist_EEs;
            dist_Ye_a = dist_Ye;
            dist_s_a = dist_s;
            dist_Mom0_a = dist_Mom0s;
            dist_M0_a = dist_MM0s;
            for (size_t i = 0; i < nlayers_a; i++) {
                double frac_of_solid_ang = 2 * std::sin(0.5 * thetas_c_h[i]) * std::sin(0.5 * thetas_c_h[i]);
//            dist_M0_a[i] = dist_E0_a[i] / (( dist_G0_a[i] - 1.0) * CGS::c*CGS::c );
                dist_E0_a[i] *= (frac_of_solid_ang / 2.); // TODO remove it and include into ID maker
                dist_M0_a[i] *= (frac_of_solid_ang / 2.);
//            ncells = (double)frac_of_solid_ang;
            }
        }

        m_nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

    }

    /// phi grid for a given layer (given number of phi cells)
    static Array getCphiGridPW( size_t ilayer ){
        size_t cil = CellsInLayer(ilayer);
        Array cphis ( cil );
        for (size_t j = 0; j < cil; j++){
            cphis[j] = (double)j * 2.0 * CGS::pi / (double)cil;
        }
        return std::move ( cphis );
    }

    /// center of theta (ctheta0) grid for a given layer and theta from dynamical evolution
    Array getCthetasGrid(size_t ilayer, Array & theta){
        return cthetas0[ilayer] + 0.5 * (2.0 * theta - 2.0 * m_theta_w);
    }

    static size_t CellsInLayer(const size_t &i_layer){
        return 2 * i_layer + 1;
    }

    static Vector getCthetaArr(size_t m_nlayers, double theta_max){
        Array thetas ( m_nlayers + 1 );
        Vector cthetas0( m_nlayers );
        for (size_t i = 0; i < m_nlayers + 1; i++){
            double fac = (double)i / (double)m_nlayers;
            thetas[i] = 2.0 * asin( fac * sin(theta_max / 2.0 ) );
        }
        for (size_t i = 0; i < m_nlayers; ++i){
            cthetas0[i] = 0.5 * ( thetas[i+1] + thetas[i] );
        }
        return std::move(cthetas0);
    }
    static Vector getThetaArr(size_t m_nlayers, double theta_max){
        Vector thetas ( m_nlayers + 1 );
        Vector cthetas0( m_nlayers );
        for (size_t i = 0; i < m_nlayers + 1; i++){
            double fac = (double)i / (double)m_nlayers;
            thetas[i] = 2.0 * asin( fac * sin(theta_max / 2.0 ) );
        }
        return std::move(thetas);
    }

    /**
     * Offset the center of the 'theta' cell accounting for the jet spreading
     * @param ilayer
     * @param theta
     * @return
     */
//    inline double getCurrCtheta(size_t ilayer, double theta){
//        return cthetas0[ilayer] + 0.5 * (2.0 * theta - 2.0 * m_theta_w);
//    }

    inline static double ctheta(double theta, size_t ilayer, size_t nlayers_pw){
        // cthetas = 0.5*(2.*arcsin(facs[0]*sin(self.joAngles[:,layer-1]/2.)) + 2.*arcsin(facs[1]*sin(self.joAngles[:,layer-1]/2.)))
//        if (theta > p_pars->theta_max ){
//            std::cerr << AT << " theta="<<theta<<" > theta_max=" << p_pars->theta_max << "\n";
//        }
//        if (std::fabs( theta - p_pars->theta_b0) > 1e-2){
//            std::cerr << AT << " theta="<<theta<<" < theta_b0=" << p_pars->theta_b0 <<"\n";
//            exit(1);
//        }
//        double ctheta = p_pars->ctheta0 + 0.5 * (2. * theta - 2. * p_pars->theta_w); // TODO WROOOONG

        double ctheta = 0.;
        if (ilayer > 0) {
            //
            double fac0 = (double)ilayer/(double)nlayers_pw;
            double fac1 = (double)(ilayer+1)/(double)nlayers_pw;
//            std::cout << std::asin(CGS::pi*3/4.) << "\n";
            if (!std::isfinite(std::sin(theta))){
                std::cerr << " sin(theta= "<<theta<<") is not finite... Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }

            double x2 = fac1*std::sin(theta / 2.);
            double xx2 = 2.*std::asin(x2);
            double x1 = fac0*std::sin(theta / 2.);
            double xx1 = 2.*std::asin(x1);

            ctheta = 0.5 * (xx1 + xx2);
            if (!std::isfinite(ctheta)){
                std::cerr << "ctheta is not finite. ctheta="<<ctheta<<" Exiting..." << "\n";
                std::cerr << AT << "\n";
                exit(1);
            }
        }
        return ctheta;
    }

public:

    METHODS method{};
    METHOD_eats m_method_eats{};
    size_t m_nlayers{};

    size_t nlayers_a{};
    Vector dist_E0_a{};
    Vector dist_Mom0_a{};
    Vector dist_M0_a{};
    Vector dist_Ye_a{};
    Vector dist_s_a{};
    Vector thetas_c_l{};
    Vector thetas_c_h{};
    Vector thetas_c{};

    size_t ncells = 0;
    size_t nlayers_pw = 0;
    Vector dist_E0_pw{};
    Vector dist_Mom0_pw{};
    Vector dist_M0_pw{};
    Vector dist_Ye_pw{};
    Vector dist_s_pw{};
    Vector cthetas0{};
    Vector theta_pw{};
    Vector thetas_h0_pw{};
    std::valarray<size_t> cil{};
//    std::valarray<size_t> cil_pw{};

    double m_theta_w = 0.;
    double m_theta_c = CGS::pi/2.;
};


struct VelocityAngularStruct{
private:
//    logger * p_log;
    std::unique_ptr<logger> p_log;
public:
    enum METHODS { iPoly22Ej, iUniform, iCustom };

    void initUniform(EjectaID & id, size_t m_nlayers, double mfac, std::string eats_method, unsigned loglevel){

        Vector & dist_thetas0 = id.getTheta0();
        Vector & dist_cthetas0 = id.getCtheta0();
        Vector & dist_moms = id.getMome();
        Vector & dist_ek = id.getEkHist();
        Vector & dist_mass = id.getMassHist();
        Vector & dist_ye = id.getYeHist();
        Vector & dist_s = id.getEntropyHist();

//        p_log = new logger(std::cout, CurrLogLevel, "LatStruct");
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");
//        p_log->set_log_level(loglevel);

        method = iUniform;
        if (dist_ek.empty() || dist_ye.empty() || dist_s.empty() || dist_moms.empty()){
            (*p_log)(LOG_ERR,AT) << " One of the input arrays is empty: "
                                 << "dist_thetas(" << dist_thetas0.size()
                                 << ") dist_moms(" << dist_moms.size()
                                 << ") dist_ye(" << dist_ye.size()
                                 << ") dist_s(" << dist_s.size()
                                 << ") dist_ek(" << dist_ek.size()
                                 << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }

        if ( (dist_moms.size() != dist_ek.size())
             || (dist_moms.size() != dist_ye.size())
             || (dist_moms.size() != dist_s.size())
             || (dist_moms.size() != dist_thetas0.size()) ){
            (*p_log)(LOG_ERR,AT) << " Mismatch in input array size "
                                 << "dist_thetas0(" << dist_thetas0.size()
                                 << ") dist_ek(" << dist_ek.size()
                                 << ") dist_ye(" << dist_ye.size()
                                 << ") dist_s(" << dist_s.size()
                                 << ") dist_moms(" << dist_moms.size()
                                 << ") Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        for (size_t i = 0; i < dist_moms.size(); i++){
            if (dist_moms[i] <= 0.0){
                (*p_log)(LOG_ERR,AT) << "bad value in initial data: beta[" << i << "] = " << dist_moms[i] << "\n";
                exit(1);
            }
        }

//        Vector dist_gams (dist_moms.size(), 0.0 );
//        for (size_t i = 0; i < dist_moms.size(); i++)
//            dist_gams[i] = EQS::Gamma(dist_moms[i]);

//        Vector dist_masses (dist_moms.size(), 0.0 );
//        for (size_t i = 0; i < dist_moms.size(); i++) {
//            double beta = EQS::BetFromMom(dist_moms[i]);
//            dist_masses[i] = dist_ek[i] / (beta * beta * CGS::c * CGS::c);
//
//        }
        nshells = dist_ek.size();
        for (size_t imom = 0; imom < nshells; imom++){
            structs.emplace_back( LatStruct() );
            structs[imom].initUniform(dist_ek[imom] * mfac,
                                      dist_moms[imom],
                                      dist_thetas0[imom],
                                      dist_mass[imom] * mfac,
                                      dist_ye[imom],
                                      dist_s[imom],
                                      m_nlayers,
                                      eats_method);
        }
    }

    /// Assuming that each shells has one layer
    void initCustom(EjectaID & id, std::string method_eats, unsigned loglevel){

        Vector & dist_thetas = id.getTheta0();
        Vector & dist_cthetas = id.getCtheta0();
        Vector & dist_moms = id.getMome();
        Vector & dist_ek = id.getEkHist();
        Vector & dist_mass = id.getMassHist();
        Vector & dist_ye = id.getYeHist();
        Vector & dist_s = id.getEntropyHist();

//        p_log = new logger(std::cout, CurrLogLevel, "LatStruct");
//        p_log = std::make_unique<logger>(std::cout, CurrLogLevel, "LatStruct");
//        p_log->set_log_level(loglevel);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");
        if (dist_moms.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs betas="
                                  << dist_moms.size() << " dist_ek=" << dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_cthetas.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs thetas="
                                  <<dist_cthetas.size() << " dist_ek="<<dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_ye.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_ye="
                                  <<dist_ye.size() << " dist_ek="<<dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_s.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_s="
                                  <<dist_s.size() << " dist_ek="<<dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }

//        Vector dist_mass(dist_cthetas.size());
//        for (size_t i = 0; i < dist_cthetas.size(); i++){
//            double beta = BetFromMom(dist_moms[i]);
//            dist_mass[i] = dist_ek[i] / (beta * beta * CGS::c * CGS::c);
//        }
        structs.emplace_back( LatStruct() );// TODO this is necessary to avoid segfault but should be replaced.
        structs[0].initCustom(dist_thetas, dist_cthetas, dist_ek, dist_ye, dist_s, dist_moms, dist_mass,
                              true, method_eats, loglevel);
        method = iCustom;
        nshells = 1.;

        (*p_log)(LOG_INFO,AT) << " setting " << "custom" << " lateral structure"
                              << " with ONE shell "<< " with beta[" << dist_moms[0] << ", " << dist_moms[nshells - 1]
                              << "] ncthetas=" << dist_cthetas.size() << " with ctheta[" << dist_cthetas[0] << ", "
                              << dist_cthetas[-1] << "]" << "\n";
    }


    /**
     *
     * @param dist_thetas  Vec[n_thetas]
     * @param dist_moms   Vec[n_betas]
     * @param dist_ek      VecVec[n_betas][n_thetas]
     * @param mfac
     */
    void initCustom(EjectaID & id, bool force_grid, std::string method_eats, unsigned loglevel){

//        ID.getTheta0(), ID.getCtheta0(), ID.getMome(),
//                ID.getEkHist(), ID.getMassHist(),  ID.getYeHist(), ID.getEntropyHist(),
//        id.getTheta0(), id.getCtheta0(), id.getMome(),
//                id.getEkCorr(), id.getMassCorr(), id.getYeCorr(), id.getEntropyCorr(),

        Vector & dist_thetas = id.getTheta0();
        Vector & dist_cthetas = id.getCtheta0();
        Vector & dist_moms = id.getMome();
        VecVector & dist_ek = id.getEkCorr();
        VecVector & dist_mass = id.getMassCorr();
        VecVector & dist_ye = id.getYeCorr();
        VecVector & dist_s = id.getEntropyCorr();

//        p_log = new logger(std::cout, CurrLogLevel, "LatStruct");
//        p_log = std::make_unique<logger>(std::cout, CurrLogLevel, "LatStruct");
//        p_log->set_log_level(loglevel);
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");
        if (dist_moms.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs betas="
                                  << dist_moms.size() << " dist_ek=" << dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_thetas.size() != dist_ek[0].size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs thetas="
                                  <<dist_thetas.size() << " dist_ek[0]="<<dist_ek[0].size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_ye[0].size() != dist_ek[0].size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_ye[0]="
                                  <<dist_ye[0].size() << " dist_ek[0]="<<dist_ek[0].size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_ye.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_ye="
                                  <<dist_ye.size() << " dist_ek="<<dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_s[0].size() != dist_ek[0].size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_s[0]="
                                  <<dist_s[0].size() << " dist_ek[0]="<<dist_ek[0].size() << "\n";
            exit(1);//throw std::runtime_error("");
        }
        if (dist_s.size() != dist_ek.size()){
            (*p_log)(LOG_ERR, AT) << "Size mismatch in ejecta distrib. arrs dist_s="
                                  <<dist_s.size() << " dist_ek="<<dist_ek.size() << "\n";
            exit(1);//throw std::runtime_error("");
        }

#if 0
        for (size_t ish = 0; ish < dist_moms.size(); ish++){
            auto nth = dist_cthetas.size();
            (*p_log)(LOG_INFO,AT)
                    << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\t"
                    << "\tThetas:  "
                    << dist_cthetas[0] << ", " << dist_cthetas[1] << ", " << dist_cthetas[2] << ", "
                    << dist_cthetas[3] << ", " << dist_cthetas[4] << ", " << dist_cthetas[5]
                    << "\t ...\t"
                    << dist_cthetas[nth-6] << ", " << dist_cthetas[nth-5] << ", " << dist_cthetas[nth-4] << ", "
                    << dist_cthetas[nth-3] << ", " << dist_cthetas[nth-2] << ", " << dist_cthetas[nth-1]
                    << "\n";
//                    << " Exiting...\n";
//            std::cout << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\n";
//            std::cout << "\tThetas" << dist_thetas << "\n";
//            std::cout << "\tEks   " << dist_ek[ish] << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
        }
        for (size_t ish = 0; ish < dist_moms.size(); ish++){
            auto nek = dist_ek[ish].size();
            (*p_log)(LOG_INFO,AT)
                    << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\t"
                    << "\tEks:  "
                    << dist_ek[ish][0] << "," << dist_ek[ish][1] << "," << dist_ek[ish][2] << ","
                    << dist_ek[ish][3] << "," << dist_ek[ish][4] << "," << dist_ek[ish][5]
                    << "\t ...\t"
                    << dist_ek[ish][nek-6] << "," << dist_ek[ish][nek-5] << "," << dist_ek[ish][nek-4] << ","
                    << dist_ek[ish][nek-3] << "," << dist_ek[ish][nek-2] << "," << dist_ek[ish][nek-1]
                    << "\n";
//                    << " Exiting...\n";
//            std::cout << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\n";
//            std::cout << "\tThetas" << dist_thetas << "\n";
//            std::cout << "\tEks   " << dist_ek[ish] << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
        }
        for (size_t ish = 0; ish < dist_moms.size(); ish++){
            auto nek = dist_ek[ish].size();
            (*p_log)(LOG_INFO,AT)
                    << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\t"
                    << "\tYe:  "
                    << dist_ye[ish][0] << "," << dist_ye[ish][1] << "," << dist_ye[ish][2] << ","
                    << dist_ye[ish][3] << "," << dist_ye[ish][4] << "," << dist_ye[ish][5]
                    << "\t ...\t"
                    << dist_ye[ish][nek-6] << "," << dist_ye[ish][nek-5] << "," << dist_ye[ish][nek-4] << ","
                    << dist_ye[ish][nek-3] << "," << dist_ye[ish][nek-2] << "," << dist_ye[ish][nek-1]
                    << "\n";
//                    << " Exiting...\n";
//            std::cout << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\n";
//            std::cout << "\tThetas" << dist_thetas << "\n";
//            std::cout << "\tEks   " << dist_ek[ish] << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
        }
        for (size_t ish = 0; ish < dist_moms.size(); ish++){
            auto nek = dist_ek[ish].size();
            (*p_log)(LOG_INFO,AT)
                    << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\t"
                    << "\tEntropy:  "
                    << dist_s[ish][0] << "," << dist_s[ish][1] << "," << dist_s[ish][2] << ","
                    << dist_s[ish][3] << "," << dist_s[ish][4] << "," << dist_s[ish][5]
                    << "\t ...\t"
                    << dist_s[ish][nek-6] << "," << dist_s[ish][nek-5] << "," << dist_s[ish][nek-4] << ","
                    << dist_s[ish][nek-3] << "," << dist_s[ish][nek-2] << "," << dist_s[ish][nek-1]
                    << "\n";
//                    << " Exiting...\n";
//            std::cout << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\n";
//            std::cout << "\tThetas" << dist_thetas << "\n";
//            std::cout << "\tEks   " << dist_ek[ish] << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
        }
#endif
        method = iCustom;
        if (dist_cthetas.empty() || dist_moms.empty() || dist_ek.empty() || dist_ye.empty()){
            std::cerr << "One of the input arrays is empty: "
                      << "dist_thetas(" << dist_cthetas.size()
                      << ") dist_ek(" << dist_ek.size()
                      << ") dist_ye(" << dist_ye.size()
                      << ") dist_s(" << dist_s.size()
                      << ") dist_moms(" << dist_moms.size() << ")\n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        size_t n_thetas = dist_cthetas.size();
        size_t n_betas = dist_moms.size();
        if ((dist_ek.size()!=n_betas) || (dist_ek[0].size()!=n_thetas)){
            std::cerr << "Input data mismatch. Expcted dist_ek[n_betas][n_thetas]" << " while got"
                      <<" dist_ek["<<dist_ek.size()<<"]["<<dist_ek[0].size()<<"] and n_betas["<<n_betas
                      <<"] n_thetas["<<dist_cthetas.size()<<"] exiting..." <<"\n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        size_t n_shells = 0;
        for (size_t imom = 0; imom < dist_moms.size(); ++imom){
            /// check if all data in the velocity shell are 0.
            double shell_sum = 0;
            for (size_t ith = 0; ith < dist_cthetas.size(); ++ith)
                shell_sum += dist_ek[imom][ith];
            if(shell_sum==0.){
                std::cerr << "shell beta=" << dist_moms[imom] << " [" << imom << "] is ignored as sum(Ek)=0.\n";
                continue;
            }
            if(dist_moms[imom] < 0.){ //  || dist_moms[imom] > 1.
                std::cerr << "shell beta=" << dist_moms[imom] << " [" << imom << "] is ignored (wrong value)\n";
                continue;
            }
//            // evaluate other needed vals for a layer
//            Vector mass_layer(dist_cthetas.size() );
//            Vector ek_layer( dist_cthetas.size() );
//            Vector ye_layer( dist_cthetas.size() );
//            Vector s_layer( dist_cthetas.size() );
            Vector mom_layer ( dist_cthetas.size(), dist_moms[imom] );//(dist_thetas.size(), EQS::Gamma(dist_moms[imom]) );
//            for (size_t ith = 0; ith < dist_cthetas.size(); ++ith){
//                ek_layer[ith] = dist_ek[imom][ith];
//                ye_layer[ith] = dist_ye[imom][ith];
//                s_layer[ith] = dist_s[imom][ith];
//                double beta = EQS::BetFromMom(dist_moms[imom]);
////                mass_layer[ith] = ek_layer[ith] / (beta * beta * CGS::c * CGS::c);
//                mass_layer[ith] = dist_mass[imom][ith];//= ek_layer[ith] / (beta * beta * CGS::c * CGS::c);
//            }
//            std::cout << dist_cthetas << "\n";
            /// emplace the layer grid
            structs.emplace_back( LatStruct() );// TODO this is necessary to avoid segfault but should be replaced.
            structs[structs.size()-1].initCustom(dist_thetas, dist_cthetas,
                                                 dist_ek[imom], dist_ye[imom],
                                                 dist_s[imom],mom_layer, dist_mass[imom],
                                                 force_grid, method_eats, loglevel);
            n_shells += 1;

#if 0
            //            Vector mom_layer ( dist_thetas.size(), EQS::Gamma(dist_moms[imom]) );
            Vector mass{};// ( dist_thetas.size() );
            Vector ek_layer{};// ( dist_thetas.size() );
            Vector _thetas{};
            for (size_t itheta = 0; itheta < n_thetas; itheta++) {
                if ( dist_ek[imom][itheta] > 0. ) {
//                    ek_layer[itheta] = dist_ek[imom][itheta];
//                    mass_layer[itheta] = ek_layer[itheta] / (dist_moms[imom] * dist_moms[imom] * CGS::c * CGS::c);
//                    ek_layer[itheta] *= mfac;
//                    mass_layer[itheta] *= mfac;
                    ek_layer.emplace_back( dist_ek[imom][itheta] );
                    mass_layer.emplace_back(ek_layer[itheta] / (dist_moms[imom] * dist_moms[imom] * CGS::c * CGS::c) );
                    ek_layer[ek_layer.size()-1] *= mfac;
                    mass_layer[ek_layer.size() - 1] *= mfac;
                    _thetas.emplace_back( dist_thetas[itheta] );
                }
            }
            if (mass_layer.size() == dist_thetas.size()) {
                n_shells += 1;
                Vector mom_layer ( _thetas.size(), EQS::Gamma(dist_moms[imom]) );
                structs.emplace_back( LatStruct() );
                structs[structs.size()-1].initCustom(_thetas, ek_layer, mom_layer, mass_layer, force_grid);
            }
            else{
                std::cout << dist_thetas << "\n";
                std::cout << ek_layer << "\n";
                std::cerr << AT << "\n shell beta=" << dist_moms[imom] << " [" << imom << "] is ignored as expected length["
                          << dist_thetas.size() << "]" << " != non-zero Ek slice[" << mass_layer.size() << "] \n";
                //std::cerr << AT << " shell="<<imom<<" with beta="<<dist_moms[imom]<<" ignored. Incomplete data.\n";
            }
#endif
        }
        if (n_shells == 0){
            std::cerr << "Zero velocity shells given. Nothing to evaluateShycnhrotronSpectrum. Exiting...\n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        nshells = n_shells;
        (*p_log)(LOG_INFO,AT) << " setting " << "custom" << " lateral structure"
                              << " nshells= " << n_shells << " with beta[" << dist_moms[0] << ", " << dist_moms[nshells - 1]
                              << "] nthetas=" << dist_thetas.size() << " with theta[" << dist_thetas[0] << ", "
                              << dist_thetas[-1] << "]" << "\n";


//        exit(1);
    }

    std::vector<LatStruct> structs{};
    size_t nshells=0;
    METHODS method{};
};
#endif





#endif //SRC_INITIAL_DATA_H
