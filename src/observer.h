//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_OBSERVER_H
#define SRC_OBSERVER_H

#include "pch.h"
#include "utils.h"
#include "interpolators.h"
#include "base_equations.h"
//#include "blastwave.h"


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
        dist_Mom0_pw.resize( nlayers_pw );
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

        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

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

        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
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
        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

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
                         nlayers,
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

//        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
        (*p_log)(LOG_INFO,AT) << " setting " << opts.at("type") << " lateral structure\n";
    }

    void initCustom(Vector & dist_thetas, Vector & dist_cthetas, Vector & dist_EEs,
                    Vector & dist_Ye, Vector & dist_s, Vector & dist_Mom0s, Vector & dist_MM0s,
                    bool force_grid, std::string eats_method, unsigned loglevel){

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "LatStruct");

        m_method_eats = setEatsMethod(eats_method);
//        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

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
                                 << dist_thetas[dist_thetas.size()-1]<<" > pi/2.="<<CGS::pi/2.<<"\n";
            exit(1);
        }
        if (dist_cthetas[0] < 0.){
            (*p_log)(LOG_ERR,AT) << "wrong initial data dist_thetas[]="
                                 << dist_cthetas[dist_cthetas.size()-1]<<" < 0.\n";
            exit(1);
        }

        if (m_method_eats == METHOD_eats::i_pw) {
            method = iCustom;
            m_theta_w = dist_thetas[dist_thetas.size() - 1];
            m_theta_c = dist_thetas[dist_thetas.size() - 1];
            (*p_log)(LOG_ERR, AT) << " m_theta_c is inaccuate\n";

            // set piece-wise
            nlayers_pw = dist_cthetas.size();
            nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;
            dist_E0_pw.resize(nlayers);
            dist_Ye_pw.resize(nlayers);
            dist_s_pw.resize(nlayers);
            dist_Mom0_pw.resize(nlayers);
            dist_M0_pw.resize(nlayers);
            setThetaGridPW();
            if ((!force_grid) &&
                (cthetas0[0] != dist_cthetas[0] || cthetas0[nlayers - 1] != dist_cthetas[nlayers - 1])) {
                std::cerr << "force_grid=" << force_grid << "\n";
                std::cerr << "dist_thetas.size()=" << dist_cthetas.size() << " cthetas0.size()=" << cthetas0.size()
                          << "\n";
                std::cerr << "dist_thetas = " << dist_cthetas << "\n";
                std::cerr << "cthetas0    = " << cthetas0 << "\n";
                (*p_log)(LOG_ERR, AT) << "angular grid mismatch: "
                                      << "cthetas0[0]=" << cthetas0[0] << " != dist_thetas[0]=" << dist_cthetas[0]
                                      << " OR cthetas[nlayers-1]=" << cthetas0[nlayers - 1]
                                      << " != dist_thetas[nlayers_pw-1]=" << dist_cthetas[nlayers - 1]
                                      << " [nlayers=" << nlayers << ", " << eats_method << "]\n";
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

            for (size_t i = 0; i < nlayers; ++i) {
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
                (thetas_c[0] != dist_cthetas[0] || thetas_c[nlayers - 1] != dist_cthetas[nlayers - 1])) {
                std::cerr << "force_grid=" << force_grid << "\n";
                std::cerr << "thetas_c.size()=" << dist_cthetas.size() << " thetas_c.size()=" << cthetas0.size()
                          << "\n";
                std::cerr << "dist_thetas = " << dist_cthetas << "\n";
                std::cerr << "thetas_c    = " << thetas_c << "\n";
                (*p_log)(LOG_ERR, AT) << "angular grid mismatch: "
                                      << "thetas_c[0]=" << thetas_c[0] << " != dist_thetas[0]=" << dist_cthetas[0]
                                      << " OR thetas_c[nlayers-1]=" << thetas_c[nlayers - 1]
                                      << " != dist_thetas[nlayers_a-1]=" << dist_cthetas[nlayers - 1]
                                      << " [nlayers=" << nlayers << ", " << eats_method << "]\n";
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

        nlayers = (m_method_eats == i_pw) ? nlayers_pw : nlayers_a;

    }

    /// phi grid for a given layer (given number of phi cells)
    static Array getCphiGridPW( size_t ilayer ){
        size_t cil = CellsInLayer(ilayer);
        Array cphis ( cil );
        for (size_t j = 0; j < cil; j++){
            cphis[j] = (double)j * 2.0 * CGS::pi / (double)cil;
        }
        return cphis;
    }

    /// center of theta (ctheta0) grid for a given layer and theta from dynamical evolution
    Array getCthetasGrid(size_t ilayer, Array & theta){
        return cthetas0[ilayer] + 0.5 * (2.0 * theta - 2.0 * m_theta_w);
    }

    static size_t CellsInLayer(const size_t &i_layer){
        return 2 * i_layer + 1;
    }

    static Vector getCthetaArr(size_t nlayers, double theta_max){
        Array thetas ( nlayers + 1 );
        Vector cthetas0( nlayers );
        for (size_t i = 0; i < nlayers + 1; i++){
            double fac = (double)i / (double)nlayers;
            thetas[i] = 2.0 * asin( fac * sin(theta_max / 2.0 ) );
        }
        for (size_t i = 0; i < nlayers; ++i){
            cthetas0[i] = 0.5 * ( thetas[i+1] + thetas[i] );
        }
        return std::move(cthetas0);
    }
    static Vector getThetaArr(size_t nlayers, double theta_max){
        Vector thetas ( nlayers + 1 );
        Vector cthetas0( nlayers );
        for (size_t i = 0; i < nlayers + 1; i++){
            double fac = (double)i / (double)nlayers;
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
    size_t nlayers{};

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

//    ~VelocityAngularStruct(){ delete p_log; }

    void initPoly( Vector coeffs, Vector coeffs2, size_t n_shells, size_t n_layers, std::string eats_method, unsigned loglevel){

//        p_log = new logger(std::cout, CurrLogLevel, "LatStruct");
//        p_log = std::make_unique<logger>(std::cout, CurrLogLevel, "LatStruct");
//        p_log->set_log_level(loglevel);

        method = iPoly22Ej;
        auto ekFromPoly22 = [](Vector& coeffs, double beta, double theta){
            double b0 = coeffs[0];
            double b1 = coeffs[1];
            double b2 = coeffs[2];
            double b3 = coeffs[3];
            double b4 = coeffs[4];
            double b5 = coeffs[5];
            double ek = b0 + (b1 * theta) + (b2 * beta) + (b3 * theta * theta) + (b4 * theta * beta) + (b5 * std::pow(beta, 1./4.));
            return std::pow(10, ek);
        };
        nshells = n_shells;
        Vector betas = TOOLS::linspace(0.001, 0.9, (int)nshells);
        for (size_t ibeta = 0; ibeta < betas.size(); ibeta++){
            structs.emplace_back( LatStruct() );
            Vector thetas = structs[ibeta].getThetaArr(n_layers, CGS::pi/2.);
            Vector cthetas = structs[ibeta].getCthetaArr(n_layers, CGS::pi/2.);
            Vector gammas (cthetas.size());
            Vector ek (cthetas.size());
            Vector ye (cthetas.size());
            Vector s (cthetas.size());
            Vector mass (cthetas.size());
            for (size_t itheta = 0; itheta < cthetas.size(); itheta++) {
                gammas[itheta] = EQS::Gamma(betas[ibeta]); // same for all layers (only mass changes, i.e., Ek)
                ek[itheta] = ekFromPoly22(coeffs, betas[ibeta], cthetas[itheta]);
                ye[itheta] = ekFromPoly22(coeffs2, betas[ibeta], cthetas[itheta]);
                s[itheta] = ekFromPoly22(coeffs2, betas[ibeta], cthetas[itheta]);
                mass[itheta] = ek[itheta] / ( betas[ibeta] * betas[ibeta] * CGS::c * CGS::c ); // in grams
            }
            structs[ibeta].initCustom( thetas, cthetas, ek, gammas, mass, ye, s, true, eats_method, loglevel);
        }
    }

//    static std::vector<std::string> list_arr_v_ns() { return {"dist_thetas", "dist_EEs", "dist_Gam0s", "dist_MM0s"}; }
//    static std::vector<std::string> list_pars_v_ns() { return {"nlayers", "mfac"}; };
//    static std::vector<std::string> list_opts_v_ns() { return {"force_grid"}; }
    void initUniform(Vector & dist_thetas0, Vector & dist_cthetas0, Vector & dist_moms,
                     Vector & dist_ek, Vector & dist_mass, Vector & dist_ye, Vector & dist_s,
                     size_t nlayers, double mfac, std::string eats_method, unsigned loglevel){

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
                                      nlayers,
                                      eats_method);
        }
    }

    /// Assuming that each shells has one layer
    void initCustom(Vector & dist_thetas, Vector & dist_cthetas, Vector & dist_moms,
                    Vector & dist_ek, Vector & dist_mass, Vector & dist_ye, Vector & dist_s,
                    std::string method_eats, unsigned loglevel){

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
    void initCustom(Vector & dist_thetas, Vector & dist_cthetas, Vector & dist_moms, VecVector & dist_ek, VecVector & dist_mass,
                    VecVector & dist_ye, VecVector & dist_s,
                    bool force_grid, std::string method_eats, unsigned loglevel){

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


        for (size_t ish = 0; ish < dist_moms.size(); ish++){
            auto nth = dist_thetas.size();
            (*p_log)(LOG_INFO,AT)
                    << "Shell[" << ish << "] beta=" << dist_moms[ish] << "\t"
                    << "\tThetas:  "
                    << dist_thetas[0] << ", " << dist_thetas[1] << ", " << dist_thetas[2] << ", "
                    << dist_thetas[3] << ", " << dist_thetas[4] << ", " << dist_thetas[5]
                    << "\t ...\t"
                    << dist_thetas[nth-6] << ", " << dist_thetas[nth-5] << ", " << dist_thetas[nth-4] << ", "
                    << dist_thetas[nth-3] << ", " << dist_thetas[nth-2] << ", " << dist_thetas[nth-1]
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
        method = iCustom;
        if (dist_thetas.empty() || dist_moms.empty() || dist_ek.empty() || dist_ye.empty()){
            std::cerr << "One of the input arrays is empty: "
                      << "dist_thetas(" << dist_thetas.size()
                      << ") dist_ek(" << dist_ek.size()
                      << ") dist_ye(" << dist_ye.size()
                      << ") dist_s(" << dist_s.size()
                      << ") dist_moms(" << dist_moms.size() << ")\n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        size_t n_thetas = dist_thetas.size();
        size_t n_betas = dist_moms.size();
        if ((dist_ek.size()!=n_betas) || (dist_ek[0].size()!=n_thetas)){
            std::cerr << "Input data mismatch. Expcted dist_ek[n_betas][n_thetas]" << " while got"
                      <<" dist_ek["<<dist_ek.size()<<"]["<<dist_ek[0].size()<<"] and n_betas["<<n_betas
                      <<"] n_thetas["<<dist_thetas.size()<<"] exiting..." <<"\n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        size_t n_shells = 0;
        for (size_t imom = 0; imom < dist_moms.size(); ++imom){
            /// check if all data in the velocity shell are 0.
            double shell_sum = 0;
            for (size_t ith = 0; ith < dist_thetas.size(); ++ith)
                shell_sum += dist_ek[imom][ith];
            if(shell_sum==0.){
                std::cerr << "shell beta=" << dist_moms[imom] << " [" << imom << "] is ignored as sum(Ek)=0.\n";
                continue;
            }
            if((dist_moms[imom] < 0.) || (dist_moms[imom] > 1.)){
                std::cerr << "shell beta=" << dist_moms[imom] << " [" << imom << "] is ignored (wrong value)\n";
                continue;
            }
            // evaluate other needed vals for a layer
            Vector mass_layer(dist_thetas.size() );
            Vector ek_layer( dist_thetas.size() );
            Vector ye_layer( dist_thetas.size() );
            Vector s_layer( dist_thetas.size() );
            Vector mom_layer ( dist_thetas.size(), dist_moms[imom] );//(dist_thetas.size(), EQS::Gamma(dist_moms[imom]) );
            for (size_t ith = 0; ith < dist_thetas.size(); ++ith){
                ek_layer[ith] = dist_ek[imom][ith];
                ye_layer[ith] = dist_ye[imom][ith];
                s_layer[ith] = dist_s[imom][ith];
//                double beta = EQS::BetFromMom(dist_moms[imom]);
                mass_layer[ith] = dist_mass[imom][ith];//= ek_layer[ith] / (beta * beta * CGS::c * CGS::c);
            }
            /// emplace the layer grid
            structs.emplace_back( LatStruct() );// TODO this is necessary to avoid segfault but should be replaced.
            structs[structs.size()-1].initCustom(dist_thetas, dist_cthetas, ek_layer, ye_layer, s_layer,
                                                 mom_layer, mass_layer,
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
            std::cerr << "Zero velocity shells given. Nothing to compute. Exiting...\n Exiting...\n";
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


struct Image {

//    std::vector<std::string> m_names{"theta", "phi", "r", "theta_j", "theta0", "mu", "xrs", "yrs", "gamma", "fluxes", "intensity", "gm", "gc", "B", "tburst", "tt"};
    std::vector<std::string> m_names{"mu", "xrs", "yrs", "intensity"};
//    enum Q { itheta, iphi, ir, itheta_j, itheta0, imu, ixr, iyr, igam, iflux, iintens, igm, igc, iB, itb, itt };
    enum Q {imu, ixr, iyr, iintens};


    VecVector m_data {};
    std::unique_ptr<logger> p_log;
    explicit Image( size_t size=1, double fill_value=0., unsigned loglevel=LOG_DEBUG) {
//        p_log = new logger(std::cout, loglevel, "Image");
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Image");
//        std::cerr << " creating image...\n";
        if (size < 1){
            (*p_log)(LOG_ERR,AT) << " required image size < 1 = "<< size << "\n";
//            std::throw_with_nested("error");
            throw std::runtime_error("error");
            exit(1);
        }
        m_size = size;
//        m_data.clear();
        m_data.resize(m_names.size());
        for (auto & arr : m_data){
            arr.resize(size, fill_value);
        }

//            m_size = size;
//            m_intens.resize(m_size, fill_value );
//            m_xrs.resize( m_size, fill_value );
//            m_yrs.resize( m_size, fill_value );
//            m_rrs.resize( m_size, fill_value );
//            m_grs.resize( m_size, fill_value );
//            m_mu.resize( m_size, fill_value );
//            m_theta_j.resize( m_size, fill_value );
//            m_thetas.resize( m_size, fill_value );
//            m_phis.resize( m_size, fill_value );
        m_f_tot = 0.0;
    }
    ~Image(){
//        std::cerr << AT << " deleting image...\n";
//        delete p_log;
    }
//    Image& operator=(const Image& other)
//    {
//        Image tmp(other);
////        swap(tmp);
//        std::cerr << AT << " copying imgae\n";
//        return *this;
//    }
//    void copy(Image & another, bool add_intensity=true){
//
//        for(size_t i = 0; i < m_names.size(); i++ )
//            m_data
//    }



    void resize(size_t size, double fill_value=0.){
        m_size = size;
        for (auto & arr : m_data){
            std::destroy(arr.begin(), arr.end());
        }
        m_data.resize(m_names.size());
        for (auto & arr : m_data){
            arr.resize(size, fill_value);
        }
    }

    void clearData(){
        if ((m_size == 0)||(m_data.empty())||(m_data[0].empty())){
            (*p_log)(LOG_ERR,AT) << "cannot clean empty image\n";
            exit(1);
        }
        for (auto & arr : m_data){
            std::fill(arr.begin(), arr.end(), 0.0);
        }
    }

    Vector & gerArr(size_t ivn){ return m_data[ivn]; }
    VecVector & getAllArrs(){ return m_data; }
    inline Vector & operator[](size_t iv_n){ return this->m_data[iv_n]; }
    inline double & operator()(size_t iv_n, size_t ii){
//        if(iv_n > m_names.size()-1){
//            if (USELOGGER){ (*p_log)(LOG_ERR) << " Access beyong memory index="<<ii<<" is above max="<<m_size-1<<"\n"; }
//            else{
//                std::cout << AT << " Access beyong memory index="<<iv_n<<" is above name_max="<<m_names.size()-1<<"\n";
//            }
//            exit(1);
//        }
//        if (ii > m_size-1){
//            if (USELOGGER){ (*p_log)(LOG_ERR) << " Access beyong memory index="<<ii<<" is above max="<<m_size-1<<"\n"; }
//            else{ std::cout << AT << " Access beyong memory index="<<ii<<" is above max="<<m_size-1<<"\n"; }
//            exit(1);
//        }
        return this->m_data[iv_n][ii];
    }
    VecVector getData(){
        if ((m_data.empty()) || (m_size == 0)){
            (*p_log)(LOG_ERR, AT) << " no data in the image. Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
        if(m_data.size() != m_names.size()){
            (*p_log)(LOG_ERR, AT) << " something is wrong with the image. Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
//        std::cout << " Reallocating [" << m_data.size() << ", " << m_data[0].size() << "]" << "\n";
        VecVector tmp;//(m_names.size(), Vector(m_size));
        tmp.resize(m_names.size());
        for(size_t i = 0; i < m_names.size(); i++){
            tmp[i].resize(m_size);
            for (size_t j = 0; j < m_size; j++){
                tmp[i][j] = m_data[i][j];
            }
        }
        return std::move( tmp );
    }
//    inline Image operator=(){
//        Image image(m_size);
//        for (size_t ivn = 0; ivn < image.m_names.size(); ivn++)
//            image(ivn) = this->gerArr(ivn);
//        return image;
//    }


    double m_f_tot{}; // total flux density in the image (in mJy) aka
    // J * (1.0 + z) / (2.0 * d_l * d_l) * CGS::cgs2mJy
    size_t m_size = 0;
//        Array m_thetas;
//        Array m_phis;
//        Array m_theta_j;
//        Array m_intens;
//        Array m_xrs;
//        Array m_yrs;
//        Array m_rrs;
//        Array m_grs;
//        Array m_mu;

};

void combineImages(Image & image, size_t ncells, size_t nlayers, std::vector<Image> & images){
    if (images.size() != nlayers){
        std::cerr << " nlayeyers="<<nlayers<<" != n_images="<<images.size()<<"\n";
        std::cerr << AT << "\n";
        exit(1);
    }
    size_t ii = 0;
    size_t icell = 0;
//    Image image(2 * struc.ncells, 0. );
    image.resize(2 * ncells, 0. );
    for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
        size_t ncells_in_layer = LatStruct::CellsInLayer(ilayer);//struc.cil[ilayer];
        auto & tmp = images[ilayer];
//        if ( tmp.m_size != 2 * ncells_in_layer ){
//            std::cerr <<  " Error !" << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
        for( size_t ipj = 0; ipj < ncells_in_layer; ipj++ ){
            for (size_t ivn = 0; ivn < image.m_names.size(); ivn++)
                (image)(ivn, ii + ipj) = tmp(ivn, ipj);
        }
        for( size_t icj = 0; icj < ncells_in_layer; icj++ ){
            for (size_t ivn = 0; ivn < image.m_names.size(); ivn++)
                (image)(ivn, ncells + ii + icj) = tmp(ivn, ncells_in_layer + icj);
        }
        ii += ncells_in_layer;
//
//        if (ilayer == 0) ii = 0;
//        else ii = struc.cil[ilayer-1];
//        for (size_t icell = 0; icell < struc.cil[ilayer]; icell++) {
////            if (ii+ncells > images[it].fluxes_layer.size()-1){ exit(1); }
//            for (size_t ivn = 0; ivn < image.m_names.size(); ivn++)
//                image(ivn, ii) = tmp(ivn, icell);
//            ii++;
//        }
        image.m_f_tot += tmp.m_f_tot;
    }
//    std::cout << image[Image::iintens].min() << ", " << image[Image::iintens].max() << "\n";
//    return image;//#std::move(image);
}



static inline double cosToSin(const double &cos_theta){
    return sqrt((1.0 - cos_theta) * (1.0 + cos_theta) );
}
static inline double arccos(const double &cos_theta){
    return 2.0 * asin( sqrt(0.5 * (1.0 - cos_theta)) );
}
static inline double obsAngle(const double &theta, const double &phi, const double &alpha_obs){
//    sin(alpha_obs) * sin(thetas) * sin(phis) + cos(alpha_obs) * cos(thetas)
    double val = sin(alpha_obs) * ( sin(theta) * sin(phi) ) + cos(alpha_obs) * cos(theta);
//    if (val < 1e-2){
//        std::cout<<AT<<"  sin(alpha_obs)="<< sin(alpha_obs)<<"\n"
//                 << " sin(theta)="<<sin(theta)<<"\n"
//                 << " sin(phi)="<<sin(phi)<<"\n"
//                 << " sin(alpha_obs) * ( sin(theta) * sin(phi) )="<<sin(alpha_obs) * ( sin(theta) * sin(phi) )<<"\n"
//                 << " cos(alpha_obs)="<<cos(alpha_obs)<<"\n"
//                 << " cos(theta)="<<cos(theta)<<"\n"
//                 << " cos(alpha_obs) * cos(theta)="<<cos(alpha_obs) * cos(theta)<<"\n";
//        std::cerr << AT << "\n";
//    }
    return val;
}
static inline double obsAngleCJ(const double &theta, const double &phi, const double &alpha_obs){
    return sin(alpha_obs) * (sin(CGS::pi - theta) * sin(phi)) + cos(alpha_obs) * cos(CGS::pi - theta);
}
static inline double imageXXs(const double &theta, const double &phi, const double &alpha_obs){
//    return -1.0 * cos(alpha_obs) * sin(theta) * sin(phi) + sin(alpha_obs) * cos(theta); // orinal
    return -1.0 * cos(alpha_obs) * ( sin(theta) * sin(phi) ) + sin(alpha_obs) * cos(theta);
}
static inline double imageYYs(const double &theta, const double &phi, const double &alpha_obs){
//    return sin(theta) * cos(phi); // original
    return +1. * sin(theta) * cos(phi);// * sin(alpha_obs);
}
static inline double imageXXsCJ(const double &theta, const double &phi, const double &alpha_obs){
    return -1.0*cos(alpha_obs)*sin(CGS::pi-theta)*sin(phi) + sin(alpha_obs)*cos(CGS::pi-theta);
}
static inline double imageYYsCJ(const double &theta, const double &phi, const double &alpha_obs){
    return +1. * sin(CGS::pi-theta) * cos(phi);
}


#endif //SRC_OBSERVER_H
