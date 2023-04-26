//
// Created by vsevolod on 08/01/23.
//

#ifndef SRC_MODEL_MAGNETAR_H
#define SRC_MODEL_MAGNETAR_H

#include "utilitites/pch.h"
#include "utilitites/logger.h"
#include "utilitites/utils.h"
#include "utilitites/interpolators.h"

#include "model_ejecta.h"

class Magnetar_OLD{
    struct Pars{
        size_t iieq = 0;
        bool useGompertz = false;

        double ns_mass = -1;
        double ns_radius = -1;
        double ns_ellipticity = -1;
        double e_bind = -1; // E_bind(self): Binding energy of the NS
        double eos_i = -1;
        double mdot0 = -1;
        double viscous_time = -1;
        double mu = -1; // G cm^3
        double omega_c = -1; // critical freq. (collapse)
        double ns_crit_beta = -1; // initial condition for EOS; beta = T/|W| parameter (Gompertz 2014)
        double omega_kep = -1;

        double omega0=-1;
    };
    Vector & m_tb_arr;
    VecVector m_data{}; // container for the solution of the evolution
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    bool is_mag_pars_set = false;
public:
    bool run_magnetar = false;
    /// RHS pars
    const static int neq = 1;
    std::vector<std::string> vars {  };
    size_t getNeq() const {
        if (!run_magnetar)
            return 0;
        else
            return neq;
    }
    enum Q_SOL { iomega };

    /// All variables
    enum Q {
        // -- dynamics ---
        itb, iOmega, iMdot, iRlc, iRmag, iRcorot, iFastness, iNdip, iNacc, iNgrav, iLdip, iLprop
    };
    std::vector<std::string> m_vnames{
            "tburst", "omega", "mdot", "r_lc", "r_mag", "r_corot", "fastness", "n_dip", "n_acc", "n_grav", "ldip", "lprop"
    };
//    static constexpr size_t NVALS = 1; // number of variables to save
    inline Vector & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    Magnetar_OLD( Vector & t_grid, int loglevel ) : m_tb_arr(t_grid) {
        ///
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Magnetar");
        /// allocate storage for all data

    }

    Vector & getTbGrid() { return m_tb_arr; }
    Vector getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
        Vector tmp{};
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            tmp.push_back(m_tb_arr[it]);
        }
//        Vector tmp2 (tmp.data(), tmp.size());
        return std::move(tmp);
    }

    void setPars(StrDbMap & pars, StrStrMap & opts){
        run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log, false, true);
        if (!run_magnetar)
            return;
        // *************************************
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize( m_tb_arr.size() );
        // **************************************
        double eos_tov_mass = getDoublePar("eos_tov_mass",pars,AT,p_log,-1,true);
        double eos_alpha = getDoublePar("eos_alpha",pars,AT,p_log,-1,true);
        double eos_beta = getDoublePar("eos_beta",pars,AT,p_log,-1,true);
        double eos_i = getDoublePar("eos_i",pars,AT,p_log,-1,true);
        double ns_b = getDoublePar("ns_b",pars,AT,p_log,-1,true);
        double ns_period = getDoublePar("ns_period",pars,AT,p_log,-1,true);
        double ns_crit_beta = getDoublePar("ns_crit_beta",pars,AT,p_log,-1,true);
        double ns_mass = getDoublePar("ns_mass",pars,AT,p_log,-1,true);
        double ns_radius = getDoublePar("ns_radius",pars,AT,p_log,-1,true);
        double ns_ellipticity = getDoublePar("ns_ellipticity",pars,AT,p_log,-1,true);
        double disk_radius = getDoublePar("disk_radius",pars,AT,p_log,-1,true);
        double disk_alpha = getDoublePar("disk_alpha",pars,AT,p_log,-1,true);
        double disk_mass0 = getDoublePar("disk_mass0",pars,AT,p_log,-1,true);
        double disk_aspect_ratio = getDoublePar("disk_aspect_ratio",pars,AT,p_log,-1,true);
        bool useGompertz = getBoolOpt("useGompertz",opts,AT,p_log, false,true);
        // **************************************

        ns_mass *= CGS::solar_m;
        disk_mass0 *= CGS::solar_m;
        eos_tov_mass *= CGS::solar_m;

        double e_bind = 0.6*ns_mass*CGS::c*CGS::c*(CGS::gravconst * ns_mass) / (ns_radius * CGS::c * CGS::c - 0.5*CGS::gravconst*ns_mass); // Lattimer and Prakash (2001)
        if (eos_i < 0){
            (*p_log)(LOG_WARN, AT) << "eos_i is not correct = "<<eos_i<<", computing from mass, preriod, radius\n";
            eos_i = 0.35 * ns_mass*ns_period*ns_radius; // SHOULD BE GIVEN! IF not: magnetar moment of inertia from Gompertz (2014) norm = 2./5
        }

        /// evaluateShycnhrotronSpectrum omega0
        double Omega0; // Gompertz et al. 2014, MNRAS 438, 240-250 ; eq. (10)
        if (!std::isfinite(ns_period) or ns_period < 0.)
            Omega0 = std::sqrt( 2. * ns_crit_beta * e_bind / eos_i );
        else
            Omega0 = 2.*CGS::pi/ns_period;
        double P0 = 2*CGS::pi/Omega0;
        /// evaluateShycnhrotronSpectrum Tem (dipole spin-down time) eq. (6) of Zhang & Meszaros (2001) [s]
        double time_spindown = (3.*CGS::c*CGS::c*CGS::c*eos_i) / (ns_b*ns_b*std::pow(ns_radius,6.)*Omega0*Omega0 );
        /// spin-down (plateu) luminocity; eq. (8) of Zhang & Meszaros (2001); [ergs/s]
        double L_em0 = (eos_i*Omega0*Omega0) / (2.*time_spindown);
        /// magnetar magnetic moment [G cm^3]
        double mu = ns_b * ns_radius*ns_radius*ns_radius;
        /// keplerian angular freq. [s^-1]
        double OmegaKep = std::sqrt(CGS::gravconst * ns_mass / (ns_radius*ns_radius*ns_radius));//**0.5
        /// evaluateShycnhrotronSpectrum viscous timescale (two versions: Gompertz 2014 and Rrayand) [s]
        double viscous_time = -1;
        if (useGompertz){
            viscous_time = disk_radius*disk_radius;
            double disk_cs = 1e7;
            viscous_time /= (3. * disk_alpha * disk_cs * disk_radius);
        } else {
            double h_disk = disk_radius * disk_aspect_ratio;
            double disk_cs = h_disk*OmegaKep*std::pow(ns_radius/disk_radius, 1.5);
            viscous_time = disk_radius*disk_radius / (3. * disk_alpha * disk_cs * h_disk);
        }
        /// evaluate initial mass accretion rate; eq (3) of King and Ritter 1998 [g/s]
        double Mdot0 = disk_mass0/viscous_time;
        /// critical_angular_velocity; Sun, Zhang & Gao (2017); eq. (25); (NS collapse for Omega < Omega_c (P>Pc); assuming ns_mass=const) [s]
        double crit = ns_mass - eos_tov_mass;
        double critical_period = -1;
        double omega_c = -1;
        if (crit > 0){
            critical_period = std::pow((ns_mass - eos_tov_mass)/(eos_alpha * eos_tov_mass),1./eos_beta);
            omega_c = 2.*CGS::pi/critical_period;
        }
        // **************************************
        p_pars->ns_mass = ns_mass;
        p_pars->ns_radius = ns_radius;
        p_pars->ns_ellipticity = ns_ellipticity;
        p_pars->eos_i = eos_i;
        p_pars->mu = mu;
        p_pars->omega_c = omega_c;
        p_pars->mdot0=Mdot0;
        p_pars->omega0=Omega0;
        p_pars->e_bind=e_bind;
        p_pars->ns_crit_beta=ns_crit_beta;
        p_pars->omega_kep=OmegaKep;
        p_pars->viscous_time=viscous_time;

        p_pars->useGompertz=useGompertz;
        // *************************************
        is_mag_pars_set = true;
    }

    void setInitConditions( double * arr, size_t i ) {
        if (p_pars->omega0 < 0){
            (*p_log)(LOG_ERR,AT)<<" omega0 is not set. Cannot set init. conditions\n";
            exit(1);
        }
        arr[i+Q_SOL::iomega] = p_pars->omega0;
    }

    void insertSolution( const double * sol, size_t it, size_t i ){
        m_data[Q::iOmega][it] = sol[i+Q_SOL::iomega];
    }

    inline double dmacc_dt(double t) {
        double mdot = p_pars->mdot0*std::exp(-t/p_pars->viscous_time);
        mdot = mdot < 1.e-10 ? 1.e-10 : mdot; // setting floor value (numerical reason)
        return mdot;
    }
    inline double radius_magnetospheric(double mdot, double r_lc){
        /// evaluateShycnhrotronSpectrum magnetospheric radius
        double r_mag = std::pow(p_pars->mu, 4./7.)
                       * std::pow(CGS::gravconst*p_pars->ns_mass, -1./7.)
                       * std::pow(mdot, -2./7.);
        r_mag = r_mag > 0.999*r_lc ? 0.999*r_lc : r_mag; // cannot exceed lc radius
        return r_mag;
    }
    inline double torque_dipol(double omega, double r_lc, double r_mag){
        /// Compute Dipole spindown torque. Eq (8) of Zhang and Meszaros 2001
        double n_dip = 0.;
        if (p_pars->useGompertz){
            // Gompertz uses the disk's alfven radius in the Bucciantini prescription, but it should actually be the alfven radius of the NS wind...
            n_dip = - 2./3. * p_pars->mu*p_pars->mu * std::pow(omega,3.)
                    / (CGS::c*CGS::c*CGS::c) * std::pow(r_lc/r_mag, 3.);
        }
        else{
            // Eq (2) of Bucciantini et al 2006 (not implemented)
            // Standard dipole spindown, no wind or disk
            n_dip = - 1./6. * p_pars->mu*p_pars->mu * std::pow(omega,3.) / (CGS::c*CGS::c*CGS::c);
        }
        n_dip = omega < p_pars->omega_c ? 0. : n_dip; // Check if NS has collapse
        return n_dip ;
    }
    inline double torque_propeller(double omega, double fastness, double r_mag, double mdot){
        double e_rot = 0.5*p_pars->eos_i*omega*omega; // Rotational energy of the NS
        double beta = e_rot / std::abs(p_pars->e_bind); // beta = T/|W| parameter (Gompertz 2014)
        /// evaluateShycnhrotronSpectrum accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = 0.;
        if ((omega > p_pars->omega_c) and (beta < p_pars->ns_crit_beta)){
            // if NS hasn't collapsed and bar-mode instability is still present
            if (r_mag > p_pars->ns_radius)
                n_acc = (1. - fastness) * std::sqrt(CGS::gravconst * p_pars->ns_mass * r_mag) * mdot; // Eq. (6)
            else
                n_acc = ((1. - omega/p_pars->omega_kep) * std::sqrt(CGS::gravconst*p_pars->ns_mass*r_mag) * mdot); // Eq. (7)
        }
        return n_acc ;
    }
    inline double torque_gws(double omega){
        /// evaluateShycnhrotronSpectrum Gravitational wave spindown torque (Zhang and Meszaros 2001)
        double n_grav = -32./5. * CGS::gravconst
                        * std::pow(CGS::pi,6.)
                        * p_pars->eos_i*p_pars->eos_i
                        * p_pars->ns_ellipticity*p_pars->ns_ellipticity
                        * std::pow(omega,5.) / std::pow(CGS::c,5.);
        n_grav = omega < 1.e4 ? 0. : n_grav; // TODO ? Why 1e4 is a limiting factor?!
        return n_grav ;
    }

    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ){
        // *********************************
        double omega = Y[i+Q_SOL::iomega]; /// dOmega/dt
        // *********************************
        double t = x;

        /// Compute accretion rate on a NS (not evolved, analytic)
        double mdot = dmacc_dt(t); // setting floor value (numerical reason)

        /// Compute light cylinder radius (for a given NS rotation)
        double r_lc = CGS::c/omega;

        /// evaluateShycnhrotronSpectrum magnetospheric radius
        double r_mag = radius_magnetospheric(mdot, r_lc);

        /// evaluateShycnhrotronSpectrum corotation radius (for a given NS mass and spin)
        double r_corot =  std::pow(CGS::gravconst * p_pars->ns_mass / (omega*omega), 1./3.);

        double fastness = std::pow(r_mag / r_corot, 1.5);

//        double e_rot = 0.5*p_pars->eos_i*omega*omega; // Rotational energy of the NS

//        double beta = e_rot / std::abs(p_pars->e_bind); // beta = T/|W| parameter (Gompertz 2014)

        /// Compute Dipole spindown torque. Eq (8) of Zhang and Meszaros 2001
        double n_dip = torque_dipol(omega, r_lc, r_mag);

        /// evaluateShycnhrotronSpectrum accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = torque_propeller(omega, fastness, r_mag, mdot);

        /// evaluateShycnhrotronSpectrum Gravitational wave spindown torque (Zhang and Meszaros 2001)
        double n_grav = torque_gws(omega);

        /// domega/dt
        double domegadt = (n_dip + n_acc + n_grav)/p_pars->eos_i;

        double dEsddt = n_dip*omega;
        double dEpropdt = n_acc*omega - CGS::gravconst*p_pars->ns_mass*mdot/r_mag;
        // **************************
        out_Y[iOmega] = domegadt;
//        out_Y[idEsddt] = -n_dip*omega;
//        out_Y[idEpropdt] = -n_acc*omega - CGS::gravconst*p_pars->ns_mass*mdot/r_mag;
        // **************************
    }

    /// TODO use one get func for all
    inline double getLdip(double tburst, const double * Y, size_t i){
        double omega = Y[i+Q_SOL::iomega]; /// dOmega/dt
        double mdot = dmacc_dt(tburst);
        double r_lc = CGS::c/omega;
        double r_mag = radius_magnetospheric(mdot, r_lc);
        /// Compute Dipole spindown torque. Eq (8) of Zhang and Meszaros 2001
        double n_dip = torque_dipol(omega, r_lc, r_mag);
        /// Dipole spindown luminosity
        double ldip = -n_dip*omega;
        return ldip;
    }
    inline double getLprop(double tburst, const double * Y, size_t i){
        double omega = Y[i+Q_SOL::iomega]; /// dOmega/dt
        double mdot = dmacc_dt(tburst);
        double r_lc = CGS::c/omega;
        double r_mag = radius_magnetospheric(mdot, r_lc);
        double r_corot =  std::pow(CGS::gravconst * p_pars->ns_mass / (omega*omega), 1./3.);
        double fastness = std::pow(r_mag / r_corot, 1.5);
        /// evaluateShycnhrotronSpectrum accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = torque_propeller(omega, fastness, r_mag, mdot);
        /// propeller luminocity (Gompertz et al. (2014))
        double lprop = - n_acc*omega - CGS::gravconst*p_pars->ns_mass*mdot/r_mag;
        return lprop;
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
        double val = interpSegLog(ia, ib, tb, m_data[vname], m_tb_arr);
        return val;
    }

    void addOtherVariables( size_t it ){
        double mdot = dmacc_dt(m_tb_arr[it]);
        double omega = m_data[Q::iOmega][it];
        double r_lc = CGS::c/omega;
        double r_mag = radius_magnetospheric(mdot, r_lc);
        double r_corot =  std::pow(CGS::gravconst * p_pars->ns_mass / (omega*omega), 1./3.);
        double fastness = std::pow(r_mag / r_corot, 1.5);
        /// Compute Dipole spindown torque. Eq (8) of Zhang and Meszaros 2001
        double n_dip = torque_dipol(omega, r_lc, r_mag);
        /// evaluateShycnhrotronSpectrum accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = torque_propeller(omega, fastness, r_mag, mdot);
        /// evaluateShycnhrotronSpectrum Gravitational wave spindown torque (Zhang and Meszaros 2001)
        double n_grav = torque_gws(omega);
        /// Dipole spindown luminosity
        double ldip = -n_dip*omega;
        /// propeller luminocity (Gompertz et al. (2014))
        double lprop = - n_acc*omega - CGS::gravconst*p_pars->ns_mass*mdot/r_mag;
        // ************************
        m_data[Q::itb][it] = m_tb_arr[it];
        m_data[Q::iMdot][it] = mdot;
        m_data[Q::iRlc][it] = r_lc;
        m_data[Q::iRmag][it] = r_mag;
        m_data[Q::iRcorot][it] = r_corot;
        m_data[Q::iFastness][iFastness] = fastness;
        m_data[Q::iNdip][it] = n_dip;
        m_data[Q::iNacc][it] = n_acc;
        m_data[Q::iNgrav][it] = n_grav;
        m_data[Q::iLdip][it] = ldip;
        m_data[Q::iLprop][it] = lprop;
        // ************************
    }

    bool isSolutionOk( double * sol, size_t i ){
        bool is_ok = true;
        if (sol[i+Q_SOL::iomega] < 0){
            (*p_log)(LOG_ERR,AT)<<" fatal. Wrong omega="<<sol[i+Q_SOL::iomega]<<"\n";
            is_ok = false;
        }
        return is_ok;
    }

    void applyUnits( double * sol, size_t i ){

    }
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

class Magnetar{
    struct Pars{

    };
    Vector m_mag_time;
    VecVector m_data{}; // container for the solution of the evolution
    std::unique_ptr<logger> p_log;
    std::unique_ptr<Pars> p_pars;
    bool is_mag_pars_set = false;
    Interp1d::METHODS mth;
    int m_loglevel{};
public:
    StrDbMap mag_pars; StrStrMap mag_opts;
    bool run_magnetar = false, load_magnetar = false, save_magnetar = false;
    /// RHS pars
    const static int neq = 1;
    std::vector<std::string> vars {  };
    size_t getNeq() const {
        if (!run_magnetar)
            return 0;
        else
            return neq;
    }
    enum Q_SOL {
        iomega,
    };

    /// All variables
    enum Q {
        // -- dynamics ---
        itb, iOmega, ildip, ilacc

    };
    std::vector<std::string> m_vnames{
            "tburst", "omega", "ildip", "ilacc"
    };
//    static constexpr size_t NVALS = 1; // number of variables to save
    inline Vector & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    Magnetar( int loglevel ){// : m_mag_time(t_grid) {
        m_loglevel = loglevel;
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Magnetar");
    }

    Vector & getTbGrid() { return m_mag_time; }
    Vector getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_mag_time;
        Vector tmp{};
        for (size_t it = 0; it < m_mag_time.size(); it = it + every_it){
            tmp.push_back(m_mag_time[it]);
        }
//        Vector tmp2 (tmp.data(), tmp.size());
        return std::move(tmp);
    }
    double getMagValInt(Q vname, double tb){
        if (!std::isfinite(tb)){
            (*p_log)(LOG_ERR,AT)<<" nan in tb!"<<"\n";
            exit(1);
        }
//        std::cerr << m_data[0][0] << " " << m_data[1][0] << " " << m_data[2][0] << " " << m_data[3][0] << "\n";
        if (!is_mag_pars_set || !load_magnetar){
            (*p_log)(LOG_ERR,AT) << " magnetar is not set/loaded\n";
            exit(1);
        }
        if (tb < m_mag_time[0] or tb > m_mag_time[m_mag_time.size() - 1])
            return 0.;
        if (tb == m_mag_time[0])
            return m_data[vname][0];
        if (tb == m_mag_time[m_mag_time.size() - 1])
            return m_data[vname][m_mag_time.size() - 1];
        /// interpolate value
        size_t ia = findIndex(tb, m_mag_time, m_mag_time.size());
        size_t ib = ia + 1;
        /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
//        double val = interpSegLin(ia, ib, tb, m_data[vname], m_mag_time);
        if ((m_data[vname][ia] == 0.) || (m_data[vname][ib]==0.))
            return 0.;
        double val = interpSegLog(ia, ib, tb, m_mag_time, m_data[vname]);
        if (!std::isfinite(val)){
            (*p_log)(LOG_ERR,AT)
                << " nan in interpolated value; vname="<<vname<<" tb="<<tb
                << " ia="<<ia<<" ib="<<ib<<" m_data["<<m_vnames[vname]<<"]["<<ia<<"]="<<m_data[vname][ia]
                << " m_data["<<m_vnames[vname]<<"]["<<ib<<"]="<<m_data[vname][ib]
                << "\n";
            std::cout << m_mag_time << "\n";
            std::cout << m_data[vname] << "\n";
            exit(1);
        }
        return val;
    }

    /// --------- LOAD MAGNETAR -------------

    void loadMagnetarEvolution(Q vname, Vector values){
        load_magnetar = true;
        // *************************************
        if ((m_mag_time.size() == 0) || (m_mag_time.size() != values.size()))
            m_mag_time.resize(values.size(), 0.0 );
        // *************************************
        if (m_data.empty()) {
            m_data.resize(m_vnames.size());
            for (auto &arr: m_data)
                arr.resize(m_mag_time.size());
        }
        // **************************************
        mth = InterpBase::select("linear",p_log);
        // **************************************
        if (vname == Q::itb)
            m_mag_time = values;
//            vecToArr(values, m_mag_time);
        else
            m_data[vname] = values;
//            vecToArr(values,m_data[vname]);
        is_mag_pars_set = true;
//        std::cerr << m_data[0][0] << " " << m_data[1][0] << " " << m_data[2][0] << " " << m_data[3][0] << "\n";
    }

    /// ------- EVOLVE MAGNETAR -------------
    void setPars(StrDbMap & pars, StrStrMap & opts, std::string working_dir, std::string parfilename, Vector ej_tarr){
        mag_pars = pars;
        mag_opts = opts;

        if ((!mag_pars.empty()) || (!mag_opts.empty())) {
            run_magnetar = getBoolOpt("run_magnetar", mag_opts, AT, p_log, false, true);
            load_magnetar = getBoolOpt("load_magnetar", mag_opts, AT, p_log, false, true);
            save_magnetar = getBoolOpt("save_magnetar", mag_opts, AT, p_log, false, true);
            if (load_magnetar && run_magnetar) {
                (*p_log)(LOG_ERR, AT) << "Cannot run and load magnetar evolution at the same time. Chose one.\n";
                exit(1);
            }
            if (load_magnetar) {
                std::string fname_magnetar_ev = getStrOpt("fname_magnetar_evol", mag_opts, AT, p_log, "", true);
                if (!std::experimental::filesystem::exists(working_dir + fname_magnetar_ev)) {
                    (*p_log)(LOG_ERR, AT) << " File not found. " + working_dir + fname_magnetar_ev << "\n";
                    exit(1);
                }
                /// load the magnetar evolution h5 file and parse the key quantity to magnetar class
                ReadMagnetarEvolutionFile dfile(working_dir + fname_magnetar_ev, m_loglevel);
                loadMagnetarEvolution(Magnetar::Q::itb, dfile.get("time"));
                loadMagnetarEvolution(Magnetar::Q::ildip, dfile.get("ldip"));
                loadMagnetarEvolution(Magnetar::Q::ilacc, dfile.get("lacc"));
                /// check if time grids mismatch
                if (getTbGrid()[0] > ej_tarr[0]) {
                    (*p_log)(LOG_ERR, AT)
                            << "Magnetar timegrid is starts too late. Loaded magnetar times[0] > the one set in parfile "
                            << "(" << getTbGrid()[0] << " > " << ej_tarr[0] << ")\n";
                    exit(1);
                }
            }
            if (run_magnetar) {
                /// pass the t_array manually as when loaded, the time grid may be different
//                setPars(mag_pars, mag_opts, ej_tarr);

            }
        }
        else {
            (*p_log)(LOG_INFO, AT) << "Magnetar is not initialized and will not be considered.\n";
        }
    }
#if 0
    void setPars(StrDbMap & pars, StrStrMap & opts){

        run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log, false, true);
        load_magnetar = getBoolOpt("load_magnetar", opts, AT,p_log, false, true);

        if (!run_magnetar)
            return;
        m_mag_time = t_arr; // TODO May Not Work. But mangetar needs its own time array...
        // *************************************
        m_data.resize(m_vnames.size());
        for (auto & arr : m_data)
            arr.resize( m_mag_time.size() );
        // **************************************
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
#endif
    void setInitConditions( double * arr, size_t i ) {
    }

    void insertSolution( const double * sol, size_t it, size_t i ){

    }
    void evaluateRhsMag( double * out_Y, size_t i, double x, double const * Y ){

    }
    void evaluateRhsPWN( double * out_Y, size_t i, double x, double const * Y ){

    }

    void addOtherVariables( size_t it ){

    }

    bool isSolutionOk( double * sol, size_t i ){
        return true;
    }

    void applyUnits( double * sol, size_t i ){

    }

    void saveMagnetarEvolution(){
        if (!run_magnetar)
            return;

        (*p_log)(LOG_ERR,AT) << "Not implemented\n";
        exit(1);

//        if (every_it < 1){
//            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"\n Exiting...";
//            std::cerr << AT << "\n";
//            exit(1);
//        }

//        auto & magnetar = p_mag;
//
//        std::vector<std::vector<double>>  tot_mag_out{};
//        std::vector<std::string> tot_names {};
//        std::unordered_map<std::string,double> group_attrs{};
//
//        auto & mag_v_ns = magnetar->m_vnames;
//        auto t_arr = magnetar->getTbGrid(every_it);
//
//        tot_mag_out.resize( mag_v_ns.size() );
//        for (size_t ivar = 0; ivar < mag_v_ns.size(); ivar++){
//            for (size_t it = 0; it < magnetar->getTbGrid().size(); it = it + every_it)
//                tot_mag_out[ivar].emplace_back( (*magnetar)[ static_cast<Magnetar::Q>(ivar) ][it] );
//        }
//
//        std::unordered_map<std::string,double> attrs {};
//
//        p_out->VectorOfVectorsH5(tot_mag_out, mag_v_ns, workingdir+fname, attrs);
    }

};

namespace PWNradiationMurase{
    /// maximum energy of electrons; Eq. (21) of Murase+15, but neglecting Y
    /// THe pair-limited lorentz -factor
    static double gamma_e_max(const double b_pwn) {
        double eta = 1;
        return sqrt(6.0*M_PI*CGS::ELEC/eta/CGS::SIGMA_T/b_pwn);
    }
    /// possible maxium energy of photons; Eq. (22) of Murase+15 in unit of [erg]
    static double e_gamma_max(const double b_pwn) {
        return gamma_e_max(b_pwn)*CGS::M_ELEC*CGS::c*CGS::c;
    }
    /// maximum energy of photons limited by gamma-gamma. Eq.42
    static double e_gamma_gamma_ani(const double T_ej) {
        return pow(CGS::M_ELEC*CGS::c*CGS::c,2.0)/2.0/CGS::K_B/T_ej;
    }
    /// The characterisitc synchrotron energy; Eq. (24) of Murase+15 in unit of [erg]
    static double e_gamma_syn_b(const double b_pwn, const double gamma_b) {
        return 3.0/2.0*CGS::H/(2.0*M_PI)*gamma_b*gamma_b*CGS::ELEC*b_pwn/CGS::M_ELEC/CGS::c;
    }

    /// The total KN cross section in unit of cm^2 (Eq.46 in Murase+15)
    static double sigma_kn(const double e_gamma) {
        double x = e_gamma/CGS::M_ELEC/CGS::c/CGS::c;
        if (x > 1.0e-3)
            return 3.0/4.0*CGS::SIGMA_T*((1.0+x)/x/x/x*(2.0*x*(1.0+x)/(1.0+2.0*x)
                                                        -log(1.0+2.0*x))+1.0/2.0/x*log(1.0+2.0*x)-(1.0+3.0*x)/pow(1.0+2.0*x,2.0));
        else
            return SIGMA_T;
    }

    /// The total BH cross section in unit of cm^2 Eq.49 in Murase+15
    static double sigma_BH_p(const double e_gamma) {
        /// Eq. (A1) and (A2) of Chodorowski+92
        double x = e_gamma/CGS::M_ELEC/CGS::c/CGS::c;
        double log2x = log(2.0*x);
        double eta = (x-2.0)/(x+2.0);
        double alpha = 1.0/137.0;
        double zeta3 = 1.2020569;

        if (x > 4.0)
            return 3.0*alpha/8.0/M_PI*SIGMA_T
                   * (28.0/9.0*log2x-218.0/27.0+pow(2.0/x,2.0)*(6.0*log2x-7.0/2.0+2.0/3.0*pow(log2x,3.0)
                                                                -pow(log2x,2.0)-1.0/3.0*M_PI*M_PI*log2x+2.0*zeta3+M_PI*M_PI/6.0)
                      -pow(2.0/x,4.0)*(3.0/16.0*log2x+1.0/8.0)
                      -pow(2.0/x,6.0)*(29.0/9.0/256.0*log2x-77.0/27.0/512.0));
        else if (x > 2.0)
            return 1.0/4.0*alpha*SIGMA_T*pow((x-2.0)/x,3.0)
                   *(1.0+eta/2.0+23.0/40.0*eta*eta+37.0/120.0*pow(eta,3.0)+61.0/192*pow(eta,4.0));
        else
            return 0.0;
    }

    /// opacity of boud-free emission
    static double kappa_bf(const double e_gamma, const double Z_eff, const int opacitymode) {
        double zeta = 1.0;//0.5 for OMK18, 1.0 for Murase+18; /* neutral fraction */

        /* See http://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html */
        //return 5.0*zeta*pow(e_gamma/EV_TO_ERG/1.0e4,-3.0)*pow(Z_eff/7.0,3.0);

        double ironopacity,seleopacity,xenonopacity,goldopacity;

        double A0,A1,A2,A3,A4;

        //iron
        if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.14801){
            A0 = -3.95919414261072;
            A1 = -2.64892215754265;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.695){
            A0 = -3.37291030805215;
            A1 = -2.75161208271434;
        }
        else{
            A0 = -1.59069232728045;
            A1 = -0.206813265289848;
        }
        ironopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));

        //selenium
        if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.84291){
            A0 = -3.78835348616654;
            A1 = -2.38423803432305;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.89764){
            A0 = -3.68604734441612;
            A1 = -2.66063041649055;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.55){
            A0 = -2.93083141712927;
            A1 = -2.60630263958148;
        }
        else{
            A0 = -1.58243193342158;
            A1 = -0.102384218895718;
        }
        seleopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));

        //xenon
        if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.32037){
            A0 = -3.07458863553159;
            A1 = -2.35739410975398;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.46141){
            A0 = -3.17731357386225;
            A1 = -2.68342346938979;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-0.282){
            A0 = -2.17345283895274;
            A1 = -2.26742402391864;
        }
        else{
            A0 = -1.55866825608716;
            A1 = -0.0467127630289143;
        }
        xenonopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));

        //gold
        if(log10(e_gamma/EV_TO_ERG/1.0e6)<-2.65645){
            A0 = -2.5113444206149;
            A1 = -2.06076550942316;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.92377){
            A0 = -2.45512766933321;
            A1 = -2.27191147638504;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<-1.09299){
            A0 = -2.32582461907071;
            A1 = -2.39063081204358;
        }
        else if(log10(e_gamma/EV_TO_ERG/1.0e6)<0.02){
            A0 = -1.55199838865465;
            A1 = -1.82076527957878;
        }
        else{
            A0 = -1.58507209319691;
            A1 = 0.0628004018301846;
        }
        goldopacity=zeta*pow(10.,A0+A1*log10(e_gamma/EV_TO_ERG/1.0e6));

        if(opacitymode==0){
            return ironopacity;
        }
        else if(opacitymode==1){
            return 0.5*ironopacity+0.5*seleopacity;
        }
        else if(opacitymode==2){
            return 0.5*xenonopacity+0.5*goldopacity;
        }
        else if(opacitymode==3){
            //approximate formula used in Murase et al. 2015
            return 5.0*zeta*pow(e_gamma/EV_TO_ERG/1.0e4,-3.0)*pow(Z_eff/7.0,3.0);
        }
        else{
            std::cerr << AT << " should not be entered: opacitymode="<<opacitymode<<"\n";
            exit(1);
        }

    }


    /// kappa_comp * sigma_comp; energy transfer coefficient from gamma rays to the thermal bath by Compton (Eq. 46 of Murase+15)
    static double gamma_ene_depo_frac_Compton(const double e_gamma) {
        double x = e_gamma/CGS::M_ELEC/CGS::c/CGS::c;

        if (x > 1.0e-3)
            return 3.0/4.0*CGS::SIGMA_T*(2.0*pow(1.0+x,2.0)/x/x/(1.0+2.0*x)
                                         -(1.0+3.0*x)/pow(1.0+2.0*x,2.0)
                                         -(1.0+x)*(2.0*x*x-2.0*x-1.0)/x/x/pow(1.0+2.0*x,2.0)
                                         -4*x*x/3.0/pow(1.0+2.0*x,3.0)
                                         -((1.0+x)/x/x/x-1.0/2.0/x+1.0/2.0/x/x/x)*log(1.0+2.0*x));
        else
            return CGS::SIGMA_T * x;

    }
    /// inelastisity of gamma rays in Compton; kappa_comp * sigma_comp / sigma_klein_nishina
    static double gamma_inelas_Compton(const double e_gamma) {
        double x = e_gamma/CGS::M_ELEC/CGS::c/CGS::c;
        return gamma_ene_depo_frac_Compton(e_gamma) / sigma_kn(e_gamma);
    }

    /// escape fraction of gamma rays interms of energy
    static double f_gamma_dep(const double e_gamma, const double rho_ej,
                              const double delta_ej, const double albd_fac, const int opacitymode) {
        double mu_e; /* electron mean molecular weight */
        //double Z_eff = 7.0; /* effective nuclear weight */
        /* this corresponds to C:O = 1:1 */
        double Z_eff;
        if(opacitymode==0){
            Z_eff = 24.21; /* effective nuclear weight */
            mu_e = 2.148;
        }
        else if(opacitymode==1){
            Z_eff = 26.74;
            mu_e = 2.2353;
        }
        else if(opacitymode==2){
            Z_eff = 53.90;
            mu_e= 2.4622;
        }
        else if(opacitymode==3){
            Z_eff = 2.0;
            mu_e= 7.0;
        }

        /// The Compton optical depth Kcomp * \rho_ej * Rej
//        double tau_Compton = (3.0-delta)/4.0/M_PI*m_ej*sigma_kn(e_gamma)/mu_e/CGS::M_PRO/r_ej/r_ej;
        double Kcomp = sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
        double tau_Compton = rho_ej*delta_ej*Kcomp;
        /// optical depth of BH pair production
//        double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/CGS::M_PRO/r_ej/r_ej;
        double KBH = (1.0+Z_eff)*sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
        double tau_BH = rho_ej*delta_ej*KBH;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
//        double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma, Z_eff, opacitymode)/r_ej/r_ej;
        double Kbf = (1.0-albd_fac)*kappa_bf(e_gamma, Z_eff, opacitymode);
        double tau_bf = rho_ej*delta_ej*Kbf;

        /// In the small inelasticity limit, a particle loses kg per
        /// interaction, so the survival fraction is (1 - kappa_gamma)^max[tau,tau^2] where
        /// max[tau,tau^2] represents the number of scatterings. In the large
        /// inelasticity limit, as in the attenuation case, the survival fraction
        /// is given by -t e .
        double power_Compton=0.0;
        if (tau_Compton > 1.0)
            power_Compton = tau_Compton*tau_Compton;
        else
            power_Compton = tau_Compton;

        /// inelasticity of Compton scattering K_comp = sigma_comp/(mu_e * mu); sigma_comp is the gamma-ray inelasticity.
        double k_comp = gamma_inelas_Compton(e_gamma);

        /// The contribution from the Compton scattering; Eq.38 in Kashiyama+2016
        double fdep_sc = (1.0 - pow(1.0 - k_comp,power_Compton));
        /// Eq.39 & 40 in Kashiyama+2016
        double fdep_ab = (1.0 - exp(-(tau_BH + tau_bf)));

        /// fdep = fdep,sc + fdep,ab
        double f_gamma_dep_tmp = fdep_sc + fdep_ab; // Eq. 39 & 40 in Kashiyama+2016

        if (f_gamma_dep_tmp > 1.0)
            return 1.0;
        else
            return f_gamma_dep_tmp;
    }

    static double spec_non_thermal(const double e_gamma, const double b_pwn,
                                   const double gamma_b, const double T_ej) {
        /* psr non-thermal emission injection spectrum "E*dF/dE/(eps_e*L_d) [erg^-1]" */
        /* We assume a broken power law with the low and high energy spectral
         * indices are -p_1 and -2 */
        /* This is motivated by more detailed calculation by Murase+15 */
        /// The source of this equations is unclear.
        /// seems to be based on https://iopscience.iop.org/article/10.1088/0004-637X/805/1/82/pdf
        /// but with no clear corresponding equation. A combination of Eq.23 & 26

        double p_1 = 1.5; //photon index -- modified
        double e_gamma_min;
//        if(diskwindmode==0){
//            e_gamma_min= 1.0*EV_TO_ERG; //0=no disk for pulsar
//        }else if(diskwindmode==2){
//            e_gamma_min = 1.0*1e-6*EV_TO_ERG; //2=disk wind
//        }

        e_gamma_min = 1.0*1e-6*EV_TO_ERG; //2=disk wind

        double e_gamma_max_tmp = e_gamma_max(b_pwn);
        double e_gamma_gamma_ani_tmp = e_gamma_gamma_ani(T_ej);
        double e_gamma_syn_b_tmp = e_gamma_syn_b(b_pwn,gamma_b);
        double norm_fac = 0.0;

        if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
            e_gamma_max_tmp = e_gamma_gamma_ani_tmp;

        if((e_gamma_max_tmp > e_gamma_syn_b_tmp)
           && (e_gamma_min < e_gamma_syn_b_tmp)){
            norm_fac = 1.0 / ( (1.0 / (2.0 - p_1) ) * (1.0 - pow(e_gamma_min/e_gamma_syn_b_tmp,2.0 - p_1)) + log(e_gamma_max_tmp/e_gamma_syn_b_tmp)) / e_gamma_syn_b_tmp;
            if((e_gamma < e_gamma_syn_b_tmp) && (e_gamma >= e_gamma_min))
                return norm_fac * pow(e_gamma/e_gamma_syn_b_tmp,-p_1+1.0);
            else if((e_gamma > e_gamma_syn_b_tmp) && (e_gamma <= e_gamma_max_tmp))
                return norm_fac * pow(e_gamma/e_gamma_syn_b_tmp,-1);
            else
                return 0.0;
        }
        else if(e_gamma_min > e_gamma_syn_b_tmp){
            norm_fac = 1.0 / log(e_gamma_max_tmp/e_gamma_min) / e_gamma_min;
            if ((e_gamma < e_gamma_max_tmp) && (e_gamma > e_gamma_min))
                return norm_fac * pow(e_gamma/e_gamma_min,-1);
            else
                return 0.0;
        }
        else{
            norm_fac = 1.0 / ((1.0/(2.0 - p_1)) * (1.0 - pow(e_gamma_min/e_gamma_max_tmp,2.0 - p_1))) / e_gamma_max_tmp;
            if ((e_gamma < e_gamma_max_tmp) && (e_gamma >= e_gamma_min))
                return norm_fac * pow(e_gamma/e_gamma_max_tmp,-p_1+1.0);
            else
                return 0.0;
        }
    }

    double f_gamma_esc(double e_gamma, double rho_ej, double delta_ej,
                       double r_ej, double albd_fac, const int opacitymode) {
        double mu_e; /* electron mean molecular weight */
        //double Z_eff = 7.0; /* effective nuclear weight */
        /* this corresponds to C:O = 1:1 */
        double Z_eff;
        if(opacitymode==0){
            Z_eff = 24.21; /* effective nuclear weight */
            mu_e = 2.148;
        }
        else if(opacitymode==1){
            Z_eff = 26.74;
            mu_e = 2.2353;
        }
        else if(opacitymode==2){
            Z_eff = 53.90;
            mu_e= 2.4622;
        }
        else if(opacitymode==3){
            Z_eff = 7.0;
            mu_e = 2.0;
        }

//        double tau_Compton = (3.0-delta)/4.0/M_PI*m_ej*sigma_kn(e_gamma)/mu_e/M_PRO/r_ej/r_ej;
//        double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/M_PRO/r_ej/r_ej;
        //double tau_bf = 2.0*(1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma,Z_eff)/r_ej/r_ej;
//        double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma,Z_eff,opacitymode)/r_ej/r_ej;

        /// The Compton optical depth Kcomp * \rho_ej * Rej
//        double tau_Compton = (3.0-delta)/4.0/M_PI*m_ej*sigma_kn(e_gamma)/mu_e/CGS::M_PRO/r_ej/r_ej;
        double Kcomp = sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
        double tau_Compton = rho_ej*delta_ej*Kcomp;
        /// optical depth of BH pair production
//        double tau_BH = (3.0-delta)/4.0/mu_e/M_PI*m_ej*(1.0+Z_eff)*sigma_BH_p(e_gamma)/CGS::M_PRO/r_ej/r_ej;
        double KBH = (1.0+Z_eff)*sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
        double tau_BH = rho_ej*delta_ej*KBH;
        /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
//        double tau_bf = (1.0-albd_fac)*(3.0-delta)/4.0/M_PI*m_ej*kappa_bf(e_gamma, Z_eff, opacitymode)/r_ej/r_ej;
        double Kbf = (1.0-albd_fac)*kappa_bf(e_gamma, Z_eff, opacitymode);
        double tau_bf = rho_ej*delta_ej*Kbf;

        double tau_abs = (1.0+gamma_inelas_Compton(e_gamma))*(tau_BH+tau_bf);
        double tau_eff = sqrt((tau_abs+tau_Compton)*tau_abs);

        double power_Compton=0.0;
        if (tau_Compton > 1.0)
            power_Compton = tau_Compton*tau_Compton;
        else
            power_Compton = tau_Compton;




        //return exp(-tau_eff);
        return exp(-(tau_BH+tau_bf))
               * (exp(-(tau_Compton)) + (1.0-exp(-(tau_Compton))) * pow(1.0-gamma_inelas_Compton(e_gamma),power_Compton));
        //return (exp(-(tau_Compton))+(1.0-exp(-(tau_Compton)))*pow(1.0-gamma_inelas_Compton(e_gamma),power_Compton));

    }


};

namespace PWN{
    /// All variables
    enum Q {
        // -- dynamics ---
        itburst, itt, iR, imom, itheta, iGamma, ibeta, iEnb, iEpwn,

    };
    std::vector<std::string> m_vnames{
            "tburst", "tt", "Rw", "mom", "theta", "Gamma", "beta", "Enebula", "Epwn"
    };
}



/// PWN model from Murase+2015
class PWNmodel{
    struct Pars{
        Pars(VecVector & data, unsigned loglevel) : m_data(data){
            p_log = std::make_unique<logger>(std::cout,std::cerr,loglevel,"PWN Pars");
        }
//        std::unique_ptr<SynchrotronAnalytic> p_syna = nullptr;
        std::unique_ptr<logger> p_log;
        Vector m_freq_arr{}; Vector m_synch_em{}; Vector m_synch_abs{};
        VecVector & m_data;
        double ctheta0=-1.;
        // --------------
        size_t iieq{};
        double Rw0{};
        double vel_w0{};
        // ---------------
        double eps_e{};
        double eps_mag{};
        double eps_th{};
        // --------------
        double rho_ej_curr{};
        double tau_ej_curr{};
        double r_ej_curr{};
        double v_ej_curr{};
        double temp_ej_curr{};
        // -------------
        double curr_ldip{};
        double curr_lacc{};
        // -------------
        double albd_fac{};
        double gamma_b{};
        int iterations{};
        int nthreads_for_frac{};
        // --------------
        double b_pwn=-1;
        bool is_init = false;
    };
//    Vector m_tb_arr;
    VecVector m_data{}; // container for the solution of the evolution
    std::unique_ptr<logger> p_log;
//    std::unique_ptr<Pars> p_pars;
    Pars * p_pars = nullptr;
    std::unique_ptr<EATS> p_eats;
    Vector frac_psr_dep_{};
    size_t m_ilayer=0;size_t m_ishell=0;size_t m_ncells=0;
    std::unique_ptr<BlastWave> & p_bw;
    unsigned m_loglevel{};
public:
    bool run_pwn = false;
    /// RHS pars
    const static int neq = 5;
    std::vector<std::string> vars {  };
    size_t getNeq() const {
        if (!run_pwn)
            return 0;
        else
            return neq;
    }
    enum Q_SOL {
        i_Rw, // Wind radius (PWN radius)
        i_mom,
        i_tt,
        i_Enb, // PWN total energy
        i_Epwn
    };

//    static constexpr size_t NVALS = 1; // number of variables to save
    inline Vector & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    PWNmodel( //Vector & tarr, double ctheta0, size_t ilayer, size_t ishell, size_t ncells, int loglevel,
              std::unique_ptr<BlastWave> & pp_bw)
            : p_bw(pp_bw) {// : m_mag_time(t_grid) {
        run_pwn = true;
        m_loglevel = p_bw->getPars()->loglevel;
        p_log = std::make_unique<logger>(std::cout, std::cerr, m_loglevel, "PWN");
        // ------------
        m_data.resize(PWN::m_vnames.size());
        for (auto & arr : m_data)
            arr.resize(p_bw->ntb(), 0.0);
        p_pars = new Pars(m_data, m_loglevel);//std::make_unique<Pars>();
        p_pars->ctheta0 = p_bw->getPars()->ctheta0;
        m_ilayer = p_bw->getPars()->ilayer;
        m_ishell = p_bw->getPars()->ishell;
        m_ncells = (size_t) p_bw->getPars()->ncells;
        /// ----------
        // Vector & tburst, Vector & tt, Vector & r, Vector & theta,Vector & m_gam, Vector & m_bet,
        // Vector & freq_arr, Vector & synch_em, Vector & synch_abs,
        // size_t & i_end_r, size_t ish, size_t il, int loglevel, void * params
        size_t i_end_r = m_data[PWN::Q::itburst].size();
        p_eats = std::make_unique<EATS>(m_data[PWN::Q::itburst],
                                        m_data[PWN::Q::itt], m_data[PWN::Q::iR], m_data[PWN::Q::itheta],
                                        m_data[PWN::Q::iGamma],m_data[PWN::Q::ibeta],
                                        p_pars->m_freq_arr, p_pars->m_synch_em, p_pars->m_synch_abs,
                                        i_end_r, 0, layer(), m_loglevel, p_pars);
        p_eats->setFluxFunc(fluxDensPW);
    }

    std::unique_ptr<EATS> & getRad(){ return p_eats; }

    Vector & getTbGrid() { return p_bw->getTbGrid(); }
    VecVector & getData(){ return m_data; }
    Vector & getData(PWN::Q q){ return m_data[q]; }

    size_t layer() const { return m_ilayer; }
    size_t shell() const { return m_ishell; }
    size_t ncells() const { return m_ncells; }

    void updateOuterBoundary(Vector & r, Vector & beta, Vector & rho, Vector & tau, Vector & temp){
        if ((r[0] < 0) || (beta[0] < 0) || (rho[0] < 0) || (tau[0] < 0) || (temp[0] < 0)){
            (*p_log)(LOG_ERR,AT) << " wrong value\n";
            exit(1);
        }
        p_pars->rho_ej_curr = rho[0];
        p_pars->tau_ej_curr = tau[0];
        p_pars->r_ej_curr = r[0];
        p_pars->v_ej_curr = beta[0] * CGS::c;
        p_pars->temp_ej_curr = temp[0];
//        std::cout << rho[0] << " "
//                  << tau[0] << " "
//                  << r[0] << " "
//                  << beta[0] * CGS::c
//                  << "\n";
//        int x = 1;
    }

    void updateMagnetar(double ldip, double lacc){
        if (ldip < 0 or !std::isfinite(ldip)){
            (*p_log)(LOG_ERR,AT) << " ldip < 0 or nan; ldip="<<ldip<<"\n";
        }
        if (lacc < 0 or !std::isfinite(lacc)){
            (*p_log)(LOG_ERR,AT) << " lacc < 0 or nan; lacc="<<lacc<<"\n";
        }
        // ----------------------
        p_pars->curr_ldip = ldip;
        p_pars->curr_lacc = lacc;
    }

    /// ------- PWN -------------
    void setPars(StrDbMap & pars, StrStrMap & opts, size_t iieq){
        double radius_wind_0 = getDoublePar("Rw0",pars,AT,p_log,-1,true); // PWN radius at t=0; [km]
//        radius_wind_0 *= (1000. * 100); // km -> cm
        double vel_wind0 = radius_wind_0 / p_bw->getTbGrid()[0]; // initial wind velocity
        double eps_e = getDoublePar("eps_e",pars,AT,p_log,-1,true); // electron acceleration efficiency
        double eps_mag = getDoublePar("eps_mag",pars,AT,p_log,-1,true); // magnetic field efficiency
        double epsth0 = getDoublePar("eps_th",pars,AT,p_log,-1,true); // initial absorption fraction
        double albd_fac = getDoublePar("albd_fac",pars,AT,p_log,-1,true); // initial albedo fraction
        double gamma_b = getDoublePar("gamma_b",pars,AT,p_log,-1,true); //break Lorentz factor of electron injection spectrum
        int iterations = (int)getDoublePar("iters",pars,AT,p_log,-1,true); //iterations for absorption spectrum calcs
        // **************************************
        p_pars->Rw0 = radius_wind_0;
        p_pars->vel_w0 = vel_wind0;
        p_pars->eps_e = eps_e;
        p_pars->eps_mag = eps_mag;
        p_pars->eps_th = epsth0;
        p_pars->albd_fac = albd_fac;
        p_pars->gamma_b = gamma_b;
        p_pars->iieq = iieq;
        p_pars->iterations = iterations;
//        p_pars->nthreads_for_frac = nthreads_frac;
        // *************************************
        frac_psr_dep_.resize(p_pars->iterations + 1, 0.0);

    }

//    std::unique_ptr<Pars> & getPars(){ return p_pars; }
    Pars *& getPars(){ return p_pars; }

    void setInitConditions( double * arr, size_t i ) {
        double beta0 = p_pars->Rw0 / p_bw->getTbGrid()[0] / CGS::c; // m_tb_arr[0] * beta0 * CGS::c;
        double mom0 = MomFromBeta(beta0);
        arr[i + Q_SOL::i_Rw] = p_pars->Rw0;
        arr[i + Q_SOL::i_mom] = mom0;
        arr[i + Q_SOL::i_Enb] = p_pars->eps_e * p_pars->curr_ldip
                              + p_pars->eps_th * p_pars->curr_lacc;
        arr[i + Q_SOL::i_Epwn] = p_pars->eps_mag * p_pars->curr_ldip;
        arr[i + Q_SOL::i_tt] = EQS::init_elapsed_time(p_pars->Rw0,
                               mom0, false);
        p_pars->b_pwn = evalCurrBpwn(arr);
        p_pars->is_init = true;
    }

    void insertSolution( const double * sol, size_t it, size_t i ){
        m_data[PWN::Q::itburst][it] = p_bw->getTbGrid()[it];
        m_data[PWN::Q::iR][it] = sol[i + Q_SOL::i_Rw];
        m_data[PWN::Q::iEnb][it] = sol[i+Q_SOL::i_Enb];
        m_data[PWN::Q::iEpwn][it] = sol[i+Q_SOL::i_Epwn];
    }

    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ){
        double r_w = Y[i + Q_SOL::i_Rw];
        double e_nb = Y[i + Q_SOL::i_Enb];
        double e_pwn = Y[i + Q_SOL::i_Epwn];
        double mom = Y[i + Q_SOL::i_mom];
        // ******************************************

        // ******************************************
        double rho_ej = p_pars->rho_ej_curr; // current ejecta density (shell?)
        double tau_ej = p_pars->tau_ej_curr; // current ejecta density (shell?)
        double r_ej = p_pars->r_ej_curr; // current ejecta density (shell?)
        double v_ej = p_pars->v_ej_curr; // current ejecta density (shell?)
        double ldip = p_pars->curr_ldip;
        double lacc = p_pars->curr_lacc;
        // ******************************************
        if (r_w > r_ej){
            r_w = r_ej;
        }
        // (see Eq. 28 in Kashiyama+16)
        double v_w = std::sqrt( (7./6.) * e_nb / (4. * CGS::pi * r_w*r_w*r_w * rho_ej) );
        // if (v_w > pow((2 * l_disk(t) * r_w / (3 - delta) / m_ej),1./3.)) # TODO finish this
        if (v_w > v_ej){
//            std::cerr << AT << "\t" << "v_w > v_ej\n";
            v_w = v_ej;
        }

        // evaluateShycnhrotronSpectrum nebula energy \int(Lem * min(1, tau_T^ej * V_ej / c))dt Eq.[28] in Eq. 28 in Kashiyama+16
        double dEnbdt = 0;
        if (tau_ej * (r_w / r_ej) > CGS::c / v_w){
            dEnbdt = (p_pars->eps_e * p_pars->curr_ldip
                   + p_pars->eps_th * p_pars->curr_lacc);
        }
        else{
            dEnbdt = (tau_ej * (r_w / r_ej) * (v_w/CGS::c)
                      * p_pars->eps_e * p_pars->curr_ldip
                      + p_pars->eps_th * p_pars->curr_lacc);
        }

        double tdyn = r_ej / v_ej;
        double dEpwndt = p_pars->eps_mag * p_pars->curr_ldip - e_pwn / tdyn;

        double dttdr = EQS::evalElapsedTime(r_w,mom,0., false);
        double drdt = v_w + r_w / x;
        /// Using pressure equilibrium, Pmag = Prad; Following approach (see Eq. 28 in Kashiyama+16)
        out_Y[i + Q_SOL::i_Rw] = drdt;
        out_Y[i + Q_SOL::i_mom] = EQS::MomFromBeta( drdt / CGS::c);
        out_Y[i + Q_SOL::i_tt] = dttdr*drdt;
        out_Y[i + Q_SOL::i_Enb] = dEnbdt;
        out_Y[i + Q_SOL::i_Epwn] = dEpwndt;
    }

    void addOtherVariables( size_t it ){
        double r_w = m_data[PWN::Q::iR][it];
        double u_b_pwn = 3.0*m_data[PWN::Q::iEpwn][it]/4.0/M_PI/r_w/r_w/r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
        double b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
        /// --------------------------------------------
        m_data[PWN::Q::iGamma][it] = EQS::GamFromMom(m_data[PWN::Q::imom][it]);
        m_data[PWN::Q::ibeta][it] = EQS::GamFromMom(m_data[PWN::Q::ibeta][it]);
        m_data[PWN::Q::ibeta][it] = p_pars->ctheta0;
    }

    /// Get current PWN magnetic field
    double evalCurrBpwn(const double * Y){
        if (!run_pwn)
            return 0.;
        double r_w = Y[p_pars->iieq + Q_SOL::i_Rw];
        double u_b_pwn = 3.0*Y[p_pars->iieq + Q_SOL::i_Epwn]/4.0/M_PI/r_w/r_w/r_w; // Eq.17 in Murase+15; Eq.34 in Kashiyama+16
        double b_pwn = pow(u_b_pwn*8.0*M_PI,0.5); //be careful: epsilon_B=epsilon_B^pre/8/PI for epsilon_B^pre used in Murase et al. 2018
        p_pars->b_pwn = b_pwn;
        return b_pwn;
    }

    double getFacPWNdep(double rho_ej, double delta_ej, double T_ej, double Ye){
        if(!std::isfinite(rho_ej) || rho_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: rho_ej="<<rho_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(delta_ej) || delta_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: delta_ej="<<delta_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(T_ej) || T_ej < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: T_ej="<<T_ej<<"\n"; exit(1);
        }
        if(!std::isfinite(Ye) || Ye < 0){
            (*p_log)(LOG_ERR,AT)<<"bad value in PWN frac calc: Ye="<<Ye<<"\n"; exit(1);
        }

        int opacitymode=0; //0=iron, 1=Y_e~0.3-0.5, 2=Y_e~0.1-0.2, 3=CO
        if (Ye <= 0.2)
            opacitymode = 2;
        else if ((Ye > 0.2) or (Ye < 0.3))
            opacitymode = 1;
        else if (Ye >= 0.3)
            opacitymode = 0;
        else if (Ye > 0.5)
            opacitymode = 3;
        else{
            (*p_log)(LOG_ERR,AT) << " error \n";
            exit(1);
        }
        return facPSRdep(rho_ej, delta_ej, T_ej, opacitymode);
    }

    double getAbsobedMagnetarLum(double lsd, double lacc, double fac_pwn_dep){
        return fac_pwn_dep * (p_pars->eps_e * lsd + p_pars->eps_e * lacc);
    }

    bool isSolutionOk( double * sol, size_t i ){
        bool is_ok = true;
        if (sol[i+Q_SOL::i_Enb] < 0){
            (*p_log)(LOG_ERR,AT)<<"Wrong Enb="<<sol[i+Q_SOL::i_Enb]<<"\n";
            is_ok = false;
        }
        return is_ok;
    }

    void applyUnits( double * sol, size_t i ){

    }

public:
    void evalElectronDistribution(){
        (*p_log)(LOG_ERR,AT)<<" not implemented\n";
        exit(1);
    }
    void computeSynchrotronSpectrum(){
        (*p_log)(LOG_ERR,AT)<<" not implemented\n";
        exit(1);
    }

private: // -------- RADIATION ----------

    double facPSRdep(const double rho_ej, const double delta_ej,
                     const double T_ej, const int opacitymode){
        if (p_pars->b_pwn < 0 || !std::isfinite(p_pars->b_pwn)){
            (*p_log)(LOG_ERR,AT)<<" b_pwn is not set\n";
            exit(1);
        }

        double e_gamma_min;
        double e_gamma_max_tmp;
        double e_gamma_gamma_ani_tmp;
        double e_gamma_syn_b_tmp;

        e_gamma_min = 1.0*EV_TO_ERG; // eV
        e_gamma_max_tmp = PWNradiationMurase::e_gamma_max(p_pars->b_pwn);
        e_gamma_gamma_ani_tmp = PWNradiationMurase::e_gamma_gamma_ani(T_ej);
        e_gamma_syn_b_tmp = PWNradiationMurase::e_gamma_syn_b(p_pars->b_pwn,p_pars->gamma_b);

        if (e_gamma_gamma_ani_tmp < e_gamma_max_tmp)
            e_gamma_max_tmp = e_gamma_gamma_ani_tmp;

        int i_max = p_pars->iterations;
//        double e_tmp = e_gamma_min;
//        double del_ln_e = 0.0;

        double albd_fac = p_pars->albd_fac;
//        int opacitymode = 0;

//        tmp(rho_ej,delta_ej,T_ej,p_pars->albd_fac,opacitymode,e_gamma_max_tmp,e_gamma_syn_b_tmp, e_gamma_min);
//        std::cout << " 11 \n";
        double frac_psr_dep_tmp = 0.0; int i = 0;
//        std::cout << " frac_psr_dep_.size(0" << frac_psr_dep_.size()<<"\n";
//        for (size_t i = 0; i <= p_pars->iterations; i++)
//            frac_psr_dep_tmp += frac_psr_dep_[i];
        int nthreads = 6;
        ///
        if (e_gamma_max_tmp > e_gamma_syn_b_tmp){
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(i_max+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=i_max;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(e_tmp,p_pars->b_pwn,p_pars->gamma_b,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==i_max)
                    frac_psr_dep_tmp_tmp = (1.0/2.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i] = frac_psr_dep_tmp_tmp;
            }
        }
        else {
            e_gamma_max_tmp = e_gamma_syn_b_tmp;
            double del_ln_e = (log(e_gamma_max_tmp)-log(e_gamma_min))/(double)(i_max+1);
#pragma omp parallel for num_threads( nthreads )
            for (i=0;i<=i_max;i++) {
                double e_tmp = e_gamma_min * exp(del_ln_e*(double)i);
                double f_gamma_dep_ = PWNradiationMurase::f_gamma_dep(e_tmp,rho_ej,delta_ej,albd_fac,opacitymode);
                double spec_non_thermal_ = PWNradiationMurase::spec_non_thermal(e_tmp,p_pars->b_pwn,p_pars->gamma_b,T_ej);
                double frac_psr_dep_tmp_ = f_gamma_dep_ * spec_non_thermal_ * e_tmp * del_ln_e;
                double frac_psr_dep_tmp_tmp = 0;
                if (i == 0 || i==i_max)
                    frac_psr_dep_tmp_tmp = (1.0/3.0) * frac_psr_dep_tmp_;
                else if (i % 2 == 0)
                    frac_psr_dep_tmp_tmp = (2.0/3.0) * frac_psr_dep_tmp_;
                else
                    frac_psr_dep_tmp_tmp = (4.0/3.0) * frac_psr_dep_tmp_;
                frac_psr_dep_[i]= frac_psr_dep_tmp_tmp;
            }
        }

        for (size_t i = 0; i <= i_max; ++i)
            frac_psr_dep_tmp += frac_psr_dep_[i];

        if (frac_psr_dep_tmp > 1.)
            frac_psr_dep_tmp = 1.;

        return frac_psr_dep_tmp;
    }

public:
    static void fluxDensPW(double & flux_dens, double & r, double & ctheta,
                           size_t ia, size_t ib, double mu, double t_obs, double nu_obs,
                           Vector ttobs, void * params){
        double nu_erg = nu_obs*4.1356655385381E-15*CGS::EV_TO_ERG;
        double b_
        fnu = PWNradiationMurase::spec_non_thermal(nu_erg, b_pwn, gamma_b, T_ej);
    }

};

/// Container for independent layers of PWN model
class PWNset{
    std::unique_ptr<Ejecta> & p_ej;
    std::unique_ptr<Output> p_out = nullptr;
    std::vector<std::unique_ptr<PWNmodel>> p_pwns{};
    std::unique_ptr<logger> p_log;
    int loglevel;
    size_t m_nlayers = 0;
    size_t m_nshells = 1;
    size_t m_ncells = 0;
    bool is_obs_pars_set = false;
    bool is_synch_computed = false;
public:
    PWNset(int loglevel, std::unique_ptr<Ejecta> & p_ej) : loglevel(loglevel), p_ej(p_ej){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "PWNset");
        p_out = std::make_unique<Output>(loglevel);
        m_nlayers = p_ej->nlayers();
        m_nshells = p_ej->nshells();
        m_ncells = p_ej->ncells();
    }
    StrDbMap pwn_pars{}; StrStrMap pwn_opts{};
    std::string workingdir{};
    bool run_pwn = false, save_pwn = false, load_pwn = false, do_ele=false, do_lc=false, do_skymap=false, do_spec=false;
    std::vector<std::unique_ptr<PWNmodel>> & getPWNs(){return p_pwns;}
    std::unique_ptr<PWNmodel> & getPWN(size_t i){return p_pwns[i];}
    size_t getNeq(){
        if (!run_pwn){
            return 0;
        }
        if (p_pwns.empty()){
            (*p_log)(LOG_ERR,AT) << "PWN is not initialized\n";
            exit(1);
        }
        size_t neq = m_nlayers * p_pwns[0]->getNeq();
        return neq;
    };
    size_t nlayers() const {return run_pwn ? m_nlayers : 0;}
    size_t nshells() const {return run_pwn ? m_nshells : 0;}
    int ncells() const {
//        return run_bws ? (int)ejectaStructs.structs[0].ncells : 0;
        return (run_pwn || load_pwn) ? (int)m_ncells : 0;
    }
    void setPars(StrDbMap & pars, StrStrMap & opts, std::string & working_dir, std::string parfilename,
                 Vector & tarr, StrDbMap & main_pars, size_t ii_eq){
        /// read pwn parameters of the pwn
        pwn_pars = pars;
        pwn_opts = opts;
        workingdir = working_dir;
        m_nlayers = p_ej->nlayers(); // as ejecta is  set later
        m_nshells = p_ej->nshells();
        if ((!pwn_pars.empty()) || (!pwn_pars.empty())) {
            run_pwn = getBoolOpt("run_pwn", pwn_opts, AT, p_log, false, true);
            save_pwn = getBoolOpt("save_pwn", pwn_opts, AT, p_log, false, true);
            load_pwn = getBoolOpt("load_pwn", pwn_opts, AT, p_log, false, true);
//            do_ele = getBoolOpt("do_ele", grb_opts, AT, p_log, false, true);
//            do_spec = getBoolOpt("do_spec", grb_opts, AT, p_log, false, true);
            do_lc = getBoolOpt("do_lc", pwn_opts, AT, p_log, false, true);
            do_skymap = getBoolOpt("do_skymap", pwn_opts, AT, p_log, false, true);
        }
        else {
            (*p_log)(LOG_INFO, AT) << "PWN is not initialized and will not be considered.\n";
        }
        if (!run_pwn)
            return;
        for(size_t il = 0; il < nlayers(); il++) {
//            auto & _x = p_ej->getShells()[il]->getBottomBW();
            p_pwns.push_back( std::make_unique<PWNmodel>( p_ej->getShells()[il]->getBW(0)) );
            p_pwns[il]->setPars(pars, opts, ii_eq);
            ii_eq += p_pwns[il]->getNeq();
        }
        // --------------------------------------------
        for (auto &key: {"n_ism", "d_l", "z", "theta_obs", "A0", "s", "r_ej", "r_ism"}) {
            if (main_pars.find(key) == main_pars.end()) {
                (*p_log)(LOG_ERR, AT) << " keyword '" << key << "' is not found in main parameters. \n";
                exit(1);
            }
            pwn_pars[key] = main_pars.at(key);
        }

    }
    void setInitConditions(double * m_InitData){
//        if (!run_pwn)
//            return;
        for (size_t il=0; il < m_nlayers; il++) {
            auto &ej_pwn = getPWNs()[il];
            ej_pwn->setInitConditions(m_InitData, ej_pwn->getPars()->iieq);
        }
    }
#if 0
    void savePWNEvolution_old(StrDbMap & main_pars){
        if (!run_pwn)
            return;

//        auto workingdir = ;
        auto fname = getStrOpt("fname_pwn", pwn_opts, AT, p_log, "", true);
        size_t every_it = (int)getDoublePar("save_pwn_every_it", pwn_pars, AT, p_log, 1, true);

        (*p_log)(LOG_INFO,AT) << "Saving PWN ...\n";

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<"Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }

        auto & pwn = getPWNs();
        std::vector<
                std::vector<
                        std::vector<double>>> tot_mag_out (pwn.size() );
        std::vector<std::string> tot_names {};
        std::vector<std::unordered_map<std::string,double>> group_attrs{};
        auto & dyn_v_ns = PWN::m_vnames;
        auto t_arr = pwn[0]->getTbGrid(every_it);
        for (size_t i = 0; i < pwn.size(); i++) {
            tot_names.push_back("layer="+std::to_string(i));
            tot_mag_out[i].resize(dyn_v_ns.size() );
            for (size_t ivar = 0; ivar < dyn_v_ns.size(); ivar++){
                for (size_t it = 0; it < pwn[i]->getTbGrid().size(); it = it + every_it)
                    tot_mag_out[i][ivar].emplace_back((*pwn[i])[ static_cast<PWN::Q>(ivar) ][it] );
            }
            ///write attributes
            auto & model = pwn[i];
            std::unordered_map<std::string,double> group_attr{
                    {"Rw0",model->getPars()->Rw0},
                    {"every_it",every_it}
            };
            group_attrs.emplace_back( group_attr );
        }
//    VecVector other_data { latStruct.cthetas0 };
//    std::vector<std::string> other_names { "cthetas0" };

        std::unordered_map<std::string, double> attrs{
                {"m_nlayers", pwn.size()}
        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: pwn_pars) { attrs[key] = value; }


        p_out->VecVectorOfVectorsAsGroupsH5(tot_mag_out, tot_names, dyn_v_ns,
                                            workingdir+fname, attrs, group_attrs);

    }
#endif
    /// output
    void savePWNdyn(StrDbMap & main_pars){
        (*p_log)(LOG_INFO,AT) << "Saving PWN BW dynamics...\n";
        auto fname = getStrOpt("fname_pwn", pwn_opts, AT, p_log, "", true);
        size_t every_it = (int)getDoublePar("save_pwn_every_it", pwn_pars, AT, p_log, 1, true);

        if (every_it < 1){
            std::cerr << " every_it must be > 1; Given every_it="<<every_it<<" \n Exiting...\n";
            std::cerr << AT << "\n";
            exit(1);
        }
        auto & models = getPWNs();
        auto & t_arr = models[0]->getTbGrid();
        std::vector<std::vector<double>> tot_dyn_out ( nshells() * nlayers() * PWN::m_vnames.size() );
        for (auto & arr : tot_dyn_out)
            arr.resize(t_arr.size());
        std::vector<std::string> arr_names{};
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            for(size_t ilayer = 0; ilayer < nlayers(); ilayer++){
                for (size_t ivar = 0; ivar < PWN::m_vnames.size(); ivar++) {
                    arr_names.push_back("shell="+std::to_string(ishell)
                                        +" layer="+std::to_string(ilayer)
                                        +" key="+BW::m_vnames[ivar]);
                    auto & pwn = models[ilayer];//->getBW(ishell);
                    for (size_t it = 0; it < t_arr.size(); it++)
                        tot_dyn_out[ii][it] = pwn->getData(static_cast<PWN::Q>(ivar))[it];
                    size_t size = tot_dyn_out[ii].size();
                    auto & x = tot_dyn_out[ii];
                    ii++;
                }
            }
        }

        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells() },
                {"m_nlayers", nlayers() },
                {"ntimes", t_arr.size() },
//                {"ncells", ncells() }
        };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: pwn_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(tot_dyn_out,arr_names,workingdir+fname,attrs);
    }

    /// INPUT
    void loadPWNDynamics(){
        auto fname = getStrOpt("fname_pwn", pwn_opts, AT, p_log, "", true);
        if (!std::experimental::filesystem::exists(workingdir+fname))
            throw std::runtime_error("File not found. " + workingdir+fname);
        Exception::dontPrint();
        H5std_string FILE_NAME(workingdir+fname);
        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        size_t nshells_ = (size_t)getDoubleAttr(file,"nshells");
        size_t nlayers_ = (size_t)getDoubleAttr(file, "nlayers");
        size_t ntimes_ = (size_t)getDoubleAttr(file, "ntimes");
        if (nshells_ != nshells()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nshells_="<<nshells_<<" expected nshells="<<nshells()<<"\n";
//            exit(1);
        }
        if (nlayers_ != nlayers()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nlayers_="<<nlayers_<<" expected nlayers_="<<nlayers()<<"\n";
//            exit(1);
        }
//        if (ntimes_ != ()){
//            (*p_log)(LOG_ERR,AT) << "Wring attribute: nlayers_="<<nlayers_<<" expected nlayers_="<<m_nlayers()<<"\n";
//            exit(1);
//        }
        //        double ntimes = getDoubleAttr(file, "ntimes");
        auto & models = getPWNs();
        for (size_t ish = 0; ish < nshells(); ish++) {
            for (size_t il = 0; il < nlayers(); il++) {
                auto & bw = models[il];//->getBW(ish);
                for (size_t ivar = 0; ivar < PWN::m_vnames.size(); ivar++) {
                    std::string key = "shell=" + std::to_string(ish)
                                      + " layer=" + std::to_string(il)
                                      + " key=" + PWN::m_vnames[ivar];
                    auto & vec = bw->getData()[static_cast<PWN::Q>(ivar)];
                    if (!vec.empty()){
                        (*p_log)(LOG_ERR,AT) << " container is not empty\n";
                    }

                    DataSet dataset = file.openDataSet(key);
                    DataType datatype = dataset.getDataType();
                    DataSpace dataspace = dataset.getSpace();
                    const int npts = dataspace.getSimpleExtentNpoints();

                    H5T_class_t classt = datatype.getClass();
                    if ( classt != 1 )
                    {
                        std::cout << key << " is not a float... you can't save this as a float." << std::endl;
                        exit(1);
                    }
                    FloatType ftype = dataset.getFloatType();
                    H5std_string order_string;
                    H5T_order_t order = ftype.getOrder( order_string);
                    size_t size = ftype.getSize();
//                    vec.resize(1);
                    double * data = new double[npts];
                    if ( order==0 && size == 4 )
                    {
                        std::cout << "NOTE: This is actually float data. We are casting to double" << std:: endl;
                        dataset.read((double*)data, PredType::IEEE_F32LE); // Our standard integer
                    }
                    else if ( order == 0 && size == 8 )
                        dataset.read(data, PredType::IEEE_F64LE);
                    else if ( order == 1 && size == 4 )
                    {
                        std::cout << "NOTE: This is actually float data. We are casting to double" << std:: endl;
                        dataset.read((double*)data, PredType::IEEE_F32BE);
                    }
                    else if ( order ==1 && size == 8 )
                        dataset.read((double*)data, PredType::IEEE_F64BE);
                    else
                        std::cout << "Did not find data type" << std::endl;
                    std::vector<double> v(data, data + npts);
                    vec = std::move( v );
//                    delete[] data;
                    dataspace.close();
                    datatype.close();
                    dataset.close();

                    if ( bw->getData()[static_cast<PWN::Q>(ivar)].empty() ){
                        std::cout << key << " faild" << std::endl;
                        exit(1);
                    }
                }
                if (bw->getData()[PWN::iR][0] == 0){
                    (*p_log)(LOG_WARN,AT) << "Loaded not evolved shell [il="<<il<<", "<<"ish="<<ish<<"] \n";
                }
            }
        }
        file.close();
//        if ( p_cumShells[0]->getBW(0)->getData()[BW::iR][0] == 0 ){
//            std::cout << p_cumShells[0]->getBW(0)->getData()[BW::iR] << "\n";
//            std::cout << " faild" << std::endl;
//            exit(1);
//        }

#if 0
        LoadH5 ldata;
        ldata.setFileName(workingdir+fname);
        ldata.setVarName("nshells");
        double nshells = ldata.getDoubleAttr("nshells");
        double m_nlayers = ldata.getDoubleAttr("m_nlayers");
        auto & models = getShells();
        for (size_t ish = 0; ish < nshells-1; ish++){
            for (size_t il = 0; il < m_nlayers-1; il++){
                auto & bw = models[il]->getBW(ish);
                for (size_t ivar = 0; ivar < BW::m_vnames.size(); ivar++) {
                    std::string key = "shell=" + std::to_string(ish)
                                    + " layer=" + std::to_string(il)
                                    + " key=" + BW::m_vnames[ivar];
                    ldata.setVarName(BW::m_vnames[ivar]);
                    bw->getData().emplace_back( std::move( ldata.getDataVDouble() ) );
                    (*p_log)(LOG_INFO,AT) << "Reading: "<<key<<"\n";
//                    bw->getData(static_cast<BW::Q>(ivar))
//                        = std::move( ldata.getDataVDouble() );
                }

//                auto & bw = models[il]->getBW(ish);
//                std::string
//                bw->getData()[]
            }
        }
#endif
        (*p_log)(LOG_INFO,AT)<<" dynamics loaded successfully\n";

    }

    void computeAndOutputObservables(StrDbMap & main_pars, StrStrMap & main_opts){
        if (run_pwn || load_pwn) {
            bool lc_freq_to_time = getBoolOpt("lc_use_freq_to_time",main_opts,AT,p_log,false,true);
            Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
            Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

            if (save_pwn)
                savePWNdyn(main_pars);

            if (load_pwn)
                loadPWNDynamics();

            if (do_ele)
                setPreComputePWNAnalyticElectronsPars(workingdir,
                    getStrOpt("fname_e_spec", pwn_opts, AT, p_log, "", true));

            if (do_spec)
                setPreComputePWNAnalyticSynchrotronPars(workingdir,
                    getStrOpt("fname_sync_spec", pwn_opts, AT, p_log, "", true));

            if (do_lc)
                computePWNlightcurve(
                        workingdir,
                        getStrOpt("fname_light_curve", pwn_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", pwn_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, pwn_pars, lc_freq_to_time);

            is_obs_pars_set = true;
        }
    }

private:
    void setPreComputePWNAnalyticElectronsPars(std::string workingdir,std::string fname){
        (*p_log)(LOG_INFO,AT) << "Computing PWN electron dists...\n";
        auto & models = getPWNs();
        for (auto & model : models) {
            model->evalElectronDistribution();
        }
    };
    void setPreComputePWNAnalyticSynchrotronPars(std::string workingdir,std::string fname){
        (*p_log)(LOG_INFO,AT) << "Computing PWN electron dists...\n";
        auto & models = getPWNs();
        for (auto & model : models) {
            model->computeSynchrotronSpectrum();
        }
        is_synch_computed = true;
    };
    void computePWNlightcurve(std::string workingdir,std::string fname, std::string fname_shells_layers,
                              Vector lc_times, Vector lc_freqs, StrDbMap & main_pars, StrDbMap & ej_pars,
                              bool lc_freq_to_time){

        Vector _times, _freqs;
        cast_times_freqs(lc_times,lc_freqs,_times,_freqs,lc_freq_to_time,p_log);

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta light curve with analytic synchrotron...\n";

        if (!is_synch_computed){
            std::cerr << " ejecta analytic electrons were not evolved. Cannot evaluateShycnhrotronSpectrum light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_obs_pars_set){
            std::cerr << " ejecta observer parameters are not set. Cannot evaluateShycnhrotronSpectrum light curve (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

//        auto & tmp = getShells()[0]->getBW(0)->getSynchAnPtr();

        std::vector< // layers / shells
                std::vector< // options
                        std::vector<double>>> // freqs*times
        out {};

        /// evaluate light curve
        auto light_curve = evalPWNLightCurves( _times, _freqs);

        /// save total lightcurve
        size_t n = _times.size();
        Vector total_fluxes (n, 0.0);
        for (size_t itnu = 0; itnu < n; ++itnu) {
            size_t ishil = 0;
            for (size_t ishell = 0; ishell < nshells(); ++ishell) {
                for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                    total_fluxes[itnu] += light_curve[ishell][ilayer][itnu];
                    ishil++;
                }
            }
        }
        std::vector<std::string> other_names { "times", "freqs", "total_fluxes" };
        VecVector out_data {_times, _freqs, total_fluxes};

        std::unordered_map<std::string,double> attrs{ {"nshells", nshells()}, {"m_nlayers", nlayers()} };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);


        /// save light curve for each shell and layer
        if (fname_shells_layers == "none")
            return;
        std::vector<std::string> group_names;
        VecVector total_fluxes_shell_layer(nshells() * nlayers());
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ++ishell) {
            for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                total_fluxes_shell_layer[ii].resize(n,0.);
                for (size_t ifnu = 0; ifnu < n; ifnu++){
                    total_fluxes_shell_layer[ii][ifnu] = light_curve[ishell][ilayer][ifnu];
                }
                ii++;
            }
        }
        total_fluxes_shell_layer.emplace_back(_times);
        total_fluxes_shell_layer.emplace_back(_freqs);
        total_fluxes_shell_layer.emplace_back(total_fluxes);

        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_fluxes");
        p_out->VectorOfVectorsH5(total_fluxes_shell_layer, group_names, workingdir+fname,  attrs);
    }

private:
    std::vector<VecVector> evalPWNLightCurves( Vector & obs_times, Vector & obs_freqs ){


        (*p_log)(LOG_INFO,AT)<<" starting ejecta light curve calculation\n";
        (*p_log)(LOG_ERR,AT)<<" not implemented\n";
//        std::vector<VecVector> light_curves{};
        exit(1);
#if 1

//        size_t nshells = p_cumShells->nshells();
//        size_t m_nlayers = p_cumShells->m_nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();
        std::vector<VecVector> light_curves(nshells()); // [ishell][i_layer][i_time]
        for (auto & arr : light_curves){
            size_t n_layers_ej = nlayers();//(p_pars->ej_method_eats == LatStruct::i_pw) ? ejectaStructs.structs[0].nlayers_pw : ejectaStructs.structs[0].nlayers_a ;
            arr.resize(n_layers_ej);
            for (auto & arrr : arr)
                arrr.resize( obs_times.size(), 0. );

        }
        double flux_pj, flux_cj; size_t ii = 0;
//        Image image;
        Image image_i ( ncells() );
        Image im_pj ( ncells() );
        Image im_cj ( ncells() );
        for (size_t ishell = 0; ishell < nshells(); ishell++){
            image_i.clearData();
            im_pj.clearData();
            im_cj.clearData();
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = ejectaStructs.structs[ishell].m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getPWNs()[ilayer];//ejectaModels[ishell][ilayer];
//                model->setEatsPars( pars, opts );
                (*p_log)(LOG_INFO,AT)
                        << " EJECTA LC ntimes=" << obs_times.size()
                        << " vel_shell=" << ishell << "/" <<nshells()-1
                        << " theta_layer=" << ilayer << "/" << nlayers()
                        << " phi_cells=" << EjectaID2::CellsInLayer(ilayer) << "\n";
                model->getRad()->evalForwardShockLightCurve(
                        p_ej->getId()->method_eats,
                        image_i, im_pj, im_cj, light_curves[ishell][ilayer], obs_times, obs_freqs);
                ii ++;
            }
        }
#endif
        return std::move( light_curves );
    }

};

#endif //SRC_MODEL_MAGNETAR_H
