//
// Created by vsevolod on 04/04/23.
//

//#include "utilitites/utils.h"
//#include "utilitites/logger.h"
//#include "utilitites/pch.h"
//#include "utilitites/H5Easy.h"

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
    enum Q{ ir,imom,itheta,ictheta,iek,imass,iye,ientr,ipress,ieps,itemp,
            itheta_c_l, itheta_c_h, itheta_c };/// additioanl quntities (will not be loaded)
    std::vector<std::string> m_names{
            "r","mom","theta","ctheta","ek","mass","ye","entr","press","eps","temp",
            "theta_c_l","theta_c_h","theta_c", /// additioanl quntities (will not be loaded)
    };
    std::vector<std::string> m_v_ns {
            "r","mom", "theta", "ctheta", "ek", "mass","ye","entr","press","eps","temp",
    };
    enum STUCT_TYPE { iadaptive, ipiecewise };
    IDTYPE idtype{};  STUCT_TYPE method_eats{};
    size_t nshells=0, nlayers=0, ncells=0;
    double theta_core=-1; double theta_wing=-1;
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
        if (m_data[iv].empty()){
            (*p_log)(LOG_ERR,AT) << " Out of bound: iv="<<iv<< " > m_data.size()="<<m_data.size()<<"\n";
            exit(1);
        }
        if (ish > nshells-1){
            (*p_log)(LOG_ERR,AT) << " Out of bound: ish="<<ish<< " > nshells-1"<<nshells-1<<"\n";
            exit(1);
        }
        if (il > nlayers-1){
            (*p_log)(LOG_ERR,AT) << " Out of bound: il="<<il<< " > nlayers-1"<<nlayers-1<<"\n";
            exit(1);
        }
        return m_data[iv][ish][il];
    }
    Vector & getVec(size_t ish, Q iv){ return m_data[iv][ish]; }
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

    static size_t CellsInLayer(const size_t &i_layer){
        return 2 * i_layer + 1;
    }
    static size_t CellsInLayer_(const size_t &i_layer){
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
        /// computeSynchrotronEmissivityAbsorptionAnalytic the number of phi cells in each 'theta' layer
        std::vector<size_t> cil;
        _evalCellsInLayer(nlayers, cil);
        size_t ncells = std::accumulate(cil.begin(), cil.end(), decltype(cil)::value_type(0));
        return ncells;
    }
private:
    /// initial grid for [a] EATS method
    void _init_a_grid(size_t ish, size_t nlayers, double theta_w) {
        _init_a_grid(m_data[Q::itheta_c_l][ish], m_data[Q::itheta_c_h][ish], m_data[Q::itheta_c][ish], nlayers, theta_w);
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
private:
    void _load_id_file(std::string & path_to_table, bool & use_1d_id, bool loadr0, double t0){
//        if (!std::experimental::filesystem::exists(path_to_table))
//            throw std::runtime_error("File not found. " + path_to_table);
        if (!std::experimental::filesystem::exists(path_to_table)) {
            (*p_log)(LOG_ERR, AT) << " Initial Data file not found: " + path_to_table << "\n";
            exit(1);
        }
        /// ---------------------------------------------
        LoadH5 ldata;
        ldata.setFileName(path_to_table);

//        size_t nshells=0, m_nlayers=0;
        if (use_1d_id){
            /// expect 1 shell with varying properties
            nshells = 1;
            for (auto & arr : m_data)
                arr.resize(m_names.size());
            for (size_t ish = 0; ish < nshells; ish++) {
                for (size_t i_v_n = 0; i_v_n < m_v_ns.size(); i_v_n++) {
                    ldata.setVarName(m_v_ns[i_v_n]);
                    m_data[i_v_n][ish] = std::move ( ldata.getDataVDouble() ); // load 1D vec
                }
            }
            nlayers = m_data[0][0].size();
            theta_wing = ldata.getAttr("theta_wing");
            theta_core = ldata.getAttr("theta_core");
            (*p_log)(LOG_INFO,AT) << " 1D ID has theta_wing="<<theta_wing<<" theta_core="<<theta_core<<"\n";
        }
        else{
            /// expect N shells with varying properties [nshells; m_nlayers]
            ldata.setVarName("mass");
            VecVector tmp = ldata.getData2Ddouble();
            nshells = tmp.size();
            nlayers = tmp[0].size();
            for (auto & arr : m_data)
                arr.resize(nshells);

            /// m_data [iv][ish][il]
            for (size_t i_v_n = 0; i_v_n < m_v_ns.size(); i_v_n++) {
                ldata.setVarName(m_v_ns[i_v_n]);
                VecVector vec = ldata.getData2Ddouble();
                if (vec.size() == nshells and vec[0].size() == nlayers)
                    m_data[i_v_n] = std::move(vec);
                else
                    m_data[i_v_n] = std::move( VecVector(nshells,Vector(nlayers,0.)) );
            }
            /// check main parameters (in case of mistakes in ID)
            for (size_t ish = 0; ish < nshells; ish++){
                for (size_t il = 0; il < nlayers; il++){
                    double mom = m_data[EjectaID2::Q::imom][ish][il];
                    double mass = m_data[EjectaID2::Q::imass][ish][il];
                    double ek = m_data[EjectaID2::Q::iek][ish][il];
                    double ctheta = m_data[EjectaID2::Q::ictheta][ish][il];
                    ///
                    if (mom < 0 || mom > 1e5 || !std::isfinite(mom)){
                        (*p_log)(LOG_ERR,AT)<<"Bad Value ID: ish="<<ish<<" il="<<il<<" mom="<<mom<<"\n";
                        exit(1);
                    }
                    if ((mass < 1 && mass != 0) || mass > 1e50 || !std::isfinite(mass)){
                        (*p_log)(LOG_ERR,AT)<<"Bad Value ID: ish="<<ish<<" il="<<il<<" mass="<<mass<<"\n";
                        exit(1);
                    }
                    if ((ek < 1e10 && ek != 0) || ek > 1e100 || !std::isfinite(ek)){
                        (*p_log)(LOG_ERR,AT)<<"Bad Value ID: ish="<<ish<<" il="<<il<<" ek="<<ek<<"\n";
                        exit(1);
                    }
                    if (ctheta < 0 || ctheta > 1.57 || !std::isfinite(ctheta)){
                        (*p_log)(LOG_ERR,AT)<<"Bad Value ID: ish="<<ish<<" il="<<il<<" ctheta="<<ctheta<<"\n";
                        exit(1);
                    }

                }
                /// we do not expect other values for structured ejecta
                theta_wing = m_data[Q::itheta][ish][nlayers-1];
                theta_core = theta_wing;
            }
        }
        /// ---------------------------
        (*p_log)(LOG_INFO,AT) << "Initial data loaded with nshells="<<nshells<<" m_nlayers="<<nlayers<<"\n";
        /// ---------------------------
        for (size_t ish = 0; ish < nshells; ish++){
//            (*p_log)(LOG_WARN,AT)<<" theta_core is NOT given in new ID. Fix it by evaluating it FROM profile!\n";
//            theta_wing = m_data[Q::itheta][ish][nlayers-1];
//            theta_core = m_data[Q::itheta][ish][nlayers-1]/4.;
//            double mom_max = std::numeric_limits<double>::max();
//            double mom_min = 0.;
//            for (size_t il = 0; il < nlayers; il++){
//                double mom = m_data[Q::imom][ish][il];
//                /// find
//                if (mom_max < mom)
//                    mom_max = mom;
//                if (mom_min > mom)
//                    mom_min = mom;
//            }


            _init_a_grid(ish, nlayers, theta_wing);
            _init_pw_grid(ish, nlayers, theta_wing);
        }
        /// ---------------------------
        (*p_log)(LOG_INFO,AT) << "Angular grids are initialized. nshells="<<nshells<<" m_nlayers="<<nlayers<<"\n";
        /// ---------------------------
        switch (method_eats) {
            case iadaptive:
//                for (size_t ish = 0; ish < nshells; ish++) {
//                    for (size_t i = 0; i < nlayers; i++) {
//                        double frac_of_solid_ang = 2
//                                                   * std::sin(0.5 * m_data[Q::itheta_c_h][ish][i])
//                                                   * std::sin(0.5 * m_data[Q::itheta_c_h][ish][i]);
//                        m_data[Q::iek][ish][i] *= (frac_of_solid_ang / 2.);// TODO remove it and include into ID maker
//                        m_data[Q::imass][ish][i] *= (frac_of_solid_ang / 2.);
//                    }
//                }
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
            (*p_log)(LOG_INFO,AT) << "Loading R0 for BWs from ID file"<<"\n";
        }
        else{
            /// pass
            (*p_log)(LOG_INFO,AT) << "Computing R0 for BWs from momenta and Global t0"<<"\n";
            for (size_t ish = 0; ish < nshells; ish++) {
                for (size_t i = 0; i < nlayers; i++) {
                    m_data[Q::ir][ish][i] = BetFromMom(m_data[Q::imom][ish][i]) * CGS::c * t0;
                }
            }
        }
        /// ---------------------------
    }
};

#endif //SRC_INITIAL_DATA_H
