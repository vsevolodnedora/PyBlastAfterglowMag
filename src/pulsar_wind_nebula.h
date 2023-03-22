//
// Created by vsevolod on 15/03/23.
//

#ifndef SRC_PULSAR_WIND_NEBULA_H
#define SRC_PULSAR_WIND_NEBULA_H

#include "pch.h"
#include "logger.h"

class PulsarWindNebula{
    struct Pars{
        double radius_w0{};
        double vel_w0{};
        // ---------------
        double eps_e{};
        double eps_mag{};
        double eps_th{};
        // --------------
        double rho_ej_curr{};
        double tau_ej_curr{};
        double r_ej_curr{};

    };
    Array m_tb_arr;
    VecArray m_data{}; // container for the solution of the evolution
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    bool is_mag_pars_set = false;
    Interp1d::METHODS mth;
public:
    bool run_pwn = false;
    /// RHS pars
    const static int neq = 1;
    std::vector<std::string> vars {  };
    size_t getNeq() const {
        if (!run_pwn)
            return 0;
        else
            return neq;
    }
    enum Q_SOL {
        iRw, // Wind radius (PWN radius)
        iEnb // PWN total energy
    };

    /// All variables
    enum Q {
        // -- dynamics ---
        itb
    };
    std::vector<std::string> m_vnames{
            "tburst"
    };
//    static constexpr size_t NVALS = 1; // number of variables to save
    inline Array & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    PulsarWindNebula( int loglevel ){// : m_mag_time(t_grid) {
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PulsarWindNebula");
    }

    Array & getTbGrid() { return m_tb_arr; }
    Array getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
        Vector tmp{};
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            tmp.push_back(m_tb_arr[it]);
        }
        Array tmp2 (tmp.data(), tmp.size());
        return std::move(tmp2);
    }
    double getValInt(Q vname, double tb){
        if (tb < m_tb_arr[0] or tb > m_tb_arr[m_tb_arr.size()-1])
            return 0.;
        if (tb == m_tb_arr[0])
            return m_data[vname][0];
        if (tb == m_tb_arr[m_tb_arr.size()-1])
            return m_data[vname][m_tb_arr.size()-1];
        /// interpolate value
        size_t ia = findIndex(tb, m_tb_arr, m_tb_arr.size());
        size_t ib = ia + 1;
        /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
//        double val = interpSegLin(ia, ib, tb, m_data[vname], m_mag_time);
        double val = interpSegLog(ia, ib, tb, m_tb_arr, m_data[vname]);
        return val;
    }

    /// ------- EVOLVE PWN -------------

    void setPars(Array & t_arr, StrDbMap & pars, StrStrMap & opts){
        run_pwn = getBoolOpt("run_pwn", opts, AT,p_log, false, true);
        if (!run_pwn)
            return;
        m_tb_arr = t_arr; // TODO May Not Work. But mangetar needs its own time array...
        // *************************************
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize( m_tb_arr.size() );
        // **************************************
        double radius_wind_0 = getDoublePar("Rw0",pars,AT,p_log,-1,true); // PWN radius at t=0; [km]
        radius_wind_0 *= (1000. * 100); // km -> cm
        double vel_wind0 = radius_wind_0 / t_arr[0]; // initial wind velocity
        double eps_e = getDoublePar("eps_e",pars,AT,p_log,-1,true); // electron acceleration efficiency
        double eps_mag = getDoublePar("eps_mag",pars,AT,p_log,-1,true); // magnetic field efficiency
        double epsth0 = getDoublePar("eps_th",pars,AT,p_log,-1,true); // initial absorption fraction
        // **************************************
        p_pars->radius_w0 = radius_wind_0;
        p_pars->vel_w0 = vel_wind0;
        p_pars->eps_e = eps_e;
        p_pars->eps_mag = eps_mag;
        p_pars->eps_th = epsth0; //0=all thermal escape, 1.0=all radiation absorbed



        double ns_mass = getDoublePar("NS_mass0",pars,AT,p_log,-1,true);
        ns_mass *= CGS::solar_m; // Msun -> g
        double ns_radius = getDoublePar("NS_radius0",pars,AT,p_log,-1,true);
        ns_radius *= (1000. * 100); // km -> cm
        double ns_inertia = getDoublePar("NS_I0",pars,AT,p_log,
                                         (2./5.*ns_mass*std::pow(ns_radius,3.0)),false);
        double ns_eps = getDoublePar("NS_eps0",pars,AT,p_log, -1,true);
        double ns_crit_beta = getDoublePar("NS_crit_beta",pars,AT,p_log, -1,true);
        double ns_B0 = getDoublePar("NS_B0",pars,AT,p_log, -1,true);
        double ns_alpha0 = getDoublePar("NS_alpha0",pars,AT,p_log, -1,true);
        double ns_sigma = getDoublePar("NS_sigma_cond",pars,AT,p_log, -1,true);
        double ns_char_length = getDoublePar("NS_char_length",pars,AT,p_log, -1,true);
        double ns_char_edens = getDoublePar("NS_char_edens",pars,AT,p_log, -1,true);
        /// Auxiliary quantity beta as defined in eq. (72) of Pons & Vigano (2019).
        double beta_ax = 1.0 / 4.0 * std::pow(ns_radius, 6.) / (ns_inertia * CGS::c*CGS::c*CGS::c);
        /// Dimensionless coefficients k_0, k_1, k_2 for a force-free magnetosphere
        /// taken from Spitkovsky (2006) and Philippov et al. (2014).
        /// For comparison, in vacuum k_0 = 0 and k_1 = k_2 = 2/3.
        Vector k_coefficients{1.0, 1.0, 1.0};
        /// Eddington luminosity assuming Thompson scattering in [erg / s].
        double eddigton_lum = 4.0 * CGS::pi * CGS::gravconst * ns_mass * CGS::mp * CGS::c / CGS::sigmaT;
        // ----
        double disk_circ_radius = getDoublePar("Disk_circ_R0",pars,AT,p_log,-1,true);
        double disk_cent_temp = getDoublePar("Disk_cent_T0",pars,AT,p_log,-1,true);
        double disk_radius0 = getDoublePar("Disk_R0",pars,AT,p_log,-1,true);
        double disk_alpha = getDoublePar("Disk_alpha",pars,AT,p_log,-1,true);
        double disk_cs = getDoublePar("Disk_cs",pars,AT,p_log,-1,true);
        double disk_aspect0 = getDoublePar("Disk_H/R0",pars,AT,p_log,-1,true);
        disk_circ_radius *= (1000. * 100); // rm -> cm
        double disk_mdot0 = getDoublePar("Disk_mdot0",pars,AT,p_log,-1,true);
        double disk_mass0 = getDoublePar("Disk_mass0",pars,AT,p_log,-1,true);
        /// --- For Collapse Criterion
//        double eos_tov_mass = getDoublePar("eos_tov_mass",pars,AT,p_log,-1,true);
//        double eos_alpha = getDoublePar("eos_alpha",pars,AT,p_log,-1,true);
//        double eos_beta = getDoublePar("eos_beta",pars,AT,p_log,-1,true);


        std::string opt;
        /// ----
        double t_v;
        opt = "method_tvis";
        std::string method_tvis = getStrOpt(opt,opts,AT,p_log,"none",true);
        if (method_tvis=="menou"){
            /// typical initial viscousity timescale of the disk; see Menou et al. 2001.
            t_v = 2080.0
                  * std::pow(disk_cent_temp / 1.0e6,-1.0)
                  * std::pow(disk_circ_radius / std::pow(10, 8.),1.0 / 2.0); // s
        }
        else if (method_tvis == "gompertz"){
            /// prescription used in Gompertz 2014
            t_v = disk_radius0 / (3. * disk_alpha * disk_cs * disk_radius0); // s
        }
        else if (method_tvis == "rrayand"){
            double disk_H = disk_radius0 * disk_aspect0;
            double disk_omega_kep = std::sqrt(CGS::gravconst * ns_mass / std::pow(ns_radius,3.));
            if (disk_cs < 0)
                disk_cs = disk_H * disk_omega_kep * std::pow(ns_radius/disk_radius0,1.5);
            t_v = disk_radius0*disk_radius0 / (3. * disk_alpha * disk_cs * disk_H); // s
        }
        else{
            (*p_log)(LOG_ERR,AT)<<" option="<<opt<<" is not recognized."
                                << "Use: "<< "menou" << " , " << "gompertz" << " , "<< "rrayand" << "\n";
            exit(1);
        }
        /// ----
        double mdot0;
        opt = "method_mdot0";
        std::string method_mdot = getStrOpt(opt,opts,AT,p_log,"none",true);
        if (method_mdot=="given"){
            /// typical initial viscousity timescale of the disk; see Menou et al. 2001.
            mdot0 = getDoublePar("Disk_mdot0",pars,AT,p_log,-1, true);
        }
        else if (method_mdot == "disk0"){
            /// prescription used in Gompertz 2014
            mdot0 = disk_mass0 / t_v;
        }
        else{
            (*p_log)(LOG_ERR,AT)<<" option="<<opt<<" is not recognized."
                                << "Use: "<< "given" << " , " << "disk0" << "\n";
            exit(1);
        }
        /// ----
        opt = "method_disk_acc";
        double disk_acc_pl = -1;
        double disk_acc_exp = -1;
        std::string method_disk_acc = getStrOpt(opt,opts,AT,p_log,"none",true);
        if (method_mdot=="pl"){
            /// Power-law index for the disk accretion rate decay.
            disk_acc_pl = getDoublePar("Disk_acc_pl",pars,AT,p_log,-1,true);
        }
        else if (method_mdot == "exp"){
            /// Power-law index for the disk outer radius evolution.
            disk_acc_exp = getDoublePar("Disk_acc_exp",pars,AT,p_log,-1,true);
        }
        else{
            (*p_log)(LOG_ERR,AT)<<" option="<<opt<<" is not recognized."
                                << "Use: "<< "pl" << " , " << "exp" << "\n";
            exit(1);
        }
        /// ----
        double ns_P0;
        opt = "method_P0";
        if (method_mdot=="given"){
            /// Power-law index for the disk accretion rate decay.
            ns_P0 = getDoublePar("NS_P0",pars,AT,p_log,-1,true);
        }
        else if (method_mdot == "crit"){
            /// NS_critical_beta = 0.27 # bar-mode instability criterion
            double OmegaCrit = std::sqrt(2.*ns_crit_beta
                                         * MEQS::e_bind(ns_mass,ns_radius)/ns_inertia);
            ns_P0 = 2. * CGS::pi * OmegaCrit;
        }
        else{
            (*p_log)(LOG_ERR,AT)<<" option="<<opt<<" is not recognized."
                                << "Use: "<< "pl" << " , " << "exp" << "\n";
            exit(1);
        }
        // *************************************

        // *************************************
        is_mag_pars_set = true;
    }



private:

    void setInitConditions( double * arr, size_t i ) {
        arr[Q_SOL::iRw] = p_pars->radius_w0;
        arr[Q_SOL::iEnb] = 0.;
    }

    void insertSolution( const double * sol, size_t it, size_t i ){

    }

    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ){
        double r_w = Y[Q_SOL::iRw];
        double e_nb = Y[Q_SOL::iEnb];
        // ******************************************
        double rho_ej = p_pars->rho_ej_curr; // current ejecta density (shell?)
        double tau_ej = p_pars->tau_ej_curr; // current ejecta density (shell?)
        double r_ej = p_pars->r_ej_curr; // current ejecta density (shell?)
        // ******************************************
        double v_w = std::sqrt( (7./6.) * e_nb / (4. * CGS::pi * r_w*r_w*r_w * rho_ej) );

        double dEnbdt = 0;
        if (tau_ej * (r_w / r_ej) > CGS::c / v_w){
            ene_nb += (p_pars->eps_e * l_d(r_ns, b_p, ome,t) + eps_th * l_disk(t)) * t * del_ln_t;
        }else{
            ene_nb += (tau_ej_tmp * (r_w / r_ej) * (v_w/C)
                       * eps_e * l_d(r_ns,b_p,ome,t) + eps_th * l_disk(t)) * t * del_ln_t;
        }

        /// Using pressure equilibrium, Pmag = Prad; Following approach (see Eq. 28 in Kashiyama+16)
        out_Y[Q_SOL::iRw] = v_w + r_w / t;
        out_Y[Q_SOL::iEnb] = dEnbdt;
    }

    void addOtherVariables( size_t it ){

    }

    bool isSolutionOk( double * sol, size_t i ){
        return true;
    }

    void applyUnits( double * sol, size_t i ){

    }

};

#endif //SRC_PULSAR_WIND_NEBULA_H
