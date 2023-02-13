//
// Created by vsevolod on 08/01/23.
//

#ifndef SRC_MAGNETAR_H
#define SRC_MAGNETAR_H

#include "pch.h"
#include "logger.h"

class Magnetar{
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
    Array & m_tb_arr;
    VecArray m_data{}; // container for the solution of the evolution
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
    inline Array & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    Magnetar( Array & t_grid, int loglevel ) : m_tb_arr(t_grid) {
        ///
        p_pars = std::make_unique<Pars>();
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Magnetar");
        /// allocate storage for all data

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

    void setPars(StrDbMap & pars, StrStrMap & opts){
        run_magnetar = getBoolOpt("run_magnetar", opts, AT,p_log,true);
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

        /// compute omega0
        double Omega0; // Gompertz et al. 2014, MNRAS 438, 240-250 ; eq. (10)
        if (!std::isfinite(ns_period) or ns_period < 0.)
            Omega0 = std::sqrt( 2. * ns_crit_beta * e_bind / eos_i );
        else
            Omega0 = 2.*CGS::pi/ns_period;
        double P0 = 2*CGS::pi/Omega0;
        /// compute Tem (dipole spin-down time) eq. (6) of Zhang & Meszaros (2001) [s]
        double time_spindown = (3.*CGS::c*CGS::c*CGS::c*eos_i) / (ns_b*ns_b*std::pow(ns_radius,6.)*Omega0*Omega0 );
        /// spin-down (plateu) luminocity; eq. (8) of Zhang & Meszaros (2001); [ergs/s]
        double L_em0 = (eos_i*Omega0*Omega0) / (2.*time_spindown);
        /// magnetar magnetic moment [G cm^3]
        double mu = ns_b * ns_radius*ns_radius*ns_radius;
        /// keplerian angular freq. [s^-1]
        double OmegaKep = std::sqrt(CGS::gravconst * ns_mass / (ns_radius*ns_radius*ns_radius));//**0.5
        /// compute viscous timescale (two versions: Gompertz 2014 and Rrayand) [s]
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
        /// compute magnetospheric radius
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
        /// compute accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
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
        /// compute Gravitational wave spindown torque (Zhang and Meszaros 2001)
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

        /// compute magnetospheric radius
        double r_mag = radius_magnetospheric(mdot, r_lc);

        /// compute corotation radius (for a given NS mass and spin)
        double r_corot =  std::pow(CGS::gravconst * p_pars->ns_mass / (omega*omega), 1./3.);

        double fastness = std::pow(r_mag / r_corot, 1.5);

//        double e_rot = 0.5*p_pars->eos_i*omega*omega; // Rotational energy of the NS

//        double beta = e_rot / std::abs(p_pars->e_bind); // beta = T/|W| parameter (Gompertz 2014)

        /// Compute Dipole spindown torque. Eq (8) of Zhang and Meszaros 2001
        double n_dip = torque_dipol(omega, r_lc, r_mag);

        /// compute accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = torque_propeller(omega, fastness, r_mag, mdot);

        /// compute Gravitational wave spindown torque (Zhang and Meszaros 2001)
        double n_grav = torque_gws(omega);

        /// domega/dt
        double domegadt = (n_dip + n_acc + n_grav)/p_pars->eos_i;

        // **************************
        out_Y[iOmega] = domegadt;
        // **************************
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
        /// compute accretion torque ( Accretion torque, taking into account the propeller model ) Eq (6-7) of Gompertz et al. 2014
        double n_acc = torque_propeller(omega, fastness, r_mag, mdot);
        /// compute Gravitational wave spindown torque (Zhang and Meszaros 2001)
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

#endif //SRC_MAGNETAR_H
