//
// Created by vsevolod on 21/04/23.
//

#ifndef SRC_EATS_H
#define SRC_EATS_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "blastwave/blastwave.h"
#include "blastwave/blastwave_components.h"
#include "utilitites/interpolators.h"
#include "utilitites/ode_solvers.h"
#include "utilitites/quadratures.h"
#include "utilitites/rootfinders.h"
#include "image.h"

//#include "blastwave.h"


static double check_emission_time( double t_e, double mu, double t_obs, Vector & mu_arr, int N ) {
    if(mu > mu_arr[N - 1]) {

        printf("mu >> 1? this should not have happened\n");
        printf("   t_obs=%.6lg t_e=%.6lg mu=%.6lg mu_arr[-1]=%.6lg\n", t_obs, t_e, mu, mu_arr[N - 1]);
        return -1;
    }
    if(mu_arr[0] >= mu) // happens only if t_e very small
    {
        printf("very small mu: mu=%.3lg, mu[0]=%.3lg\n", mu, mu_arr[0]);
        //return t_obs / (1.0 - mu); // so return small t_e limit
        return -1;
    }
    return t_e;
}

struct ImageExtend{
    ImageExtend(double time_, double freq_, size_t nshells_,size_t nlayers_, int loglevel){
        nshells = nshells_;
        nlayers = nlayers_;
        time = time_;
        freq = freq_;

        fluxes_shells.resize(nshells, 0.);
        raw_data.resize(nshells, VecVector(IMG::m_names.size()));

        p_log = std::make_unique<logger> (std::cout, std::cerr, loglevel, "image");

    }
    size_t nshells = 0;
    size_t nlayers = 0;
    double freq = 0.;
    double time = 0.;
    std::unique_ptr<logger> p_log;
    /// ---
    double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min(), ymax = std::numeric_limits<double>::min();
    Vector fluxes_shells{};
    Vector xmins_shells{};
    Vector xmaxs_shells{};
    Vector ymins_shells{};
    Vector ymaxs_shells{};
    double total_flux;
    double xc = 0.;
    double yc = 0.;
    std::vector<VecVector> raw_data{}; // [ish][Q][i]
    VecVector binned_data{};
    VecVector interp_data{};

    void binImage(size_t nx, size_t ny ) {

        // Initialize histogram with zeros
        binned_data.resize(nx, std::vector<double>(ny, 0.0));
        double binWidthX = (xmax - xmin) / (double)nx;
        double binWidthY = (ymax - ymin) / (double)ny;

        for (size_t ish = 0; ish < nshells; ish++) {

            auto & x = raw_data[ish][IMG::Q::ixr];
            auto & z = raw_data[ish][IMG::Q::iyr];
            auto & f = raw_data[ish][IMG::Q::iintens];

            for (size_t i = 0; i < raw_data[ish][0].size(); i++) {
                if ((x[i] >= xmin) && (x[i] <= xmax) &&
                    (z[i] >= ymin) && (z[i] <= ymax)) {

                    auto binIndexX = static_cast<size_t>((x[i] - xmin) / binWidthX);
                    auto binIndexY = static_cast<size_t>((z[i] - ymin) / binWidthY);

                    // Clamp to the last bin if value is exactly at the upper boundary
                    if (binIndexX == nx) {
                        binIndexX = nx - 1;
                    }
                    if (binIndexY == ny) {
                        binIndexY = ny - 1;
                    }

                    binned_data[binIndexX][binIndexY] += f[i];
                }
            }
        }
    }

    void interpImage(size_t nx, size_t ny){
        double dx = (xmax - xmin) / (double)(nx - 1);
        double dy = (ymax - ymin) / (double)(ny - 1);

        interp_data.resize(nx, std::vector<double>(ny));
        for (size_t ish = 0; ish < nshells; ish++) {
            auto & X = raw_data[ish][IMG::Q::ixr];
            auto & Y = raw_data[ish][IMG::Q::iyr];
            auto & F = raw_data[ish][IMG::Q::iintens];
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    double x = xmin + i * dx;
                    double y = ymin + j * dy;
                    interp_data[i][j] += interp2d(X, Y, F, x, y);
                }
            }
        }
    }

//    void interpImage(size_t nx, size_t ny){
//        interp_data.resize(nx, std::vector<double>(ny, 0.0));
//        for (size_t ish = 0; ish < nshells; ish++) {
//            auto & x = raw_data[ish][IMG::Q::ixr];
//            auto & y = raw_data[ish][IMG::Q::iyr];
//            auto & f = raw_data[ish][IMG::Q::iintens];
//            Interp2d interp2D(x, y, f);
//            for (size_t ix = 0; ix < nx; ix++){
//                double xp = xmin + (double)ix * (xmax - xmin) / (double)nx;
//                for (size_t iy = 0; iy < ny; iy++){
//                    double yp = ymin + (double)iy * (ymax - ymin) / (double)ny;
//                    double val = interp2D.InterpolateBilinear(xp, yp);
//                    interp_data[ix][iy] += val;
//                }
//            }
//        }
//    }

    void evalXcYc(){
        /// compute weighted average
        double _num_x = 0.;
        double _num_y = 0.;
        double _denum = 0.;
        for (size_t i = 0; i < nshells; i++){
            if (raw_data[i][IMG::Q::iintens].empty()){
                std::cerr << AT << " intensity array is 0\n";
                exit(1);
            }
            size_t nn = raw_data[i][IMG::ixr].size();
            for (size_t j = 0; j < nn; j++){
                _num_x += (raw_data[i][IMG::ixr][j] * raw_data[i][IMG::iintens][j]);
                _num_y += (raw_data[i][IMG::iyr][j] * raw_data[i][IMG::iintens][j]);
                _denum += raw_data[i][IMG::iintens][j];
            }
        }
        xc = _num_x / _denum;
        yc = _num_y / _denum;
    }
};

void saveRawImage(ImageExtend & im, size_t itinu, std::string workdir,
                  StrDbMap & main_pars, StrDbMap & ej_pars,
                  std::unique_ptr<logger> & p_log){

    if (im.raw_data.size() < 1){
        (*p_log)(LOG_ERR,AT) << " passed empty data. Cannot save\n";
        exit(1);
    }

    std::string fpath = workdir + "raw_skymap_"+std::to_string(itinu)+".h5";

    Output::remove_file_if_existis(fpath,p_log);

    H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

    /// add data for each shell
    size_t nshells_ = 0;
    for (size_t ish = 0; ish < im.nshells; ish++){
        std::string group_name = "shell=" + std::to_string(ish);
        if ((!im.raw_data[ish].empty()) && (!im.raw_data[ish][0].empty())) {
            Output::addGroupWith1Ddata(im.raw_data[ish], group_name, IMG::m_names, file);
            nshells_++;
        }
    }
    /// add fluxes for from shell
    std::string fname = "fluxs";
    Output::addVector(im.fluxes_shells, fname, file);

    /// add attributes from model parameters
    std::unordered_map<std::string,double> attrs{
            {"flux", im.total_flux},
            {"nshells", nshells_},
            {"nlayers", im.nlayers},
            {"time", im.time},
            {"freq", im.freq},
            {"xmin", im.xmin},
            {"ymin", im.ymin},
            {"xmax", im.xmax},
            {"ymax", im.ymax},
    };
    for (auto& [key, value]: main_pars) { attrs[key] = value; }
    for (auto& [key, value]: ej_pars) { attrs[key] = value; }

    Output::addStrDbMap(attrs, file);

    file.close();


}

/// To be used once I figure out inteproaltion of the image
void saveImages(std::vector<ImageExtend> & ims, Vector & times, Vector & freqs,
                StrDbMap & main_pars, StrDbMap & ej_pars,
                std::string fpath, std::unique_ptr<logger> & p_log){

//    std::cout << "S"

    // remove the old file (so to avoid data piling up)
    Output::remove_file_if_existis(fpath, p_log);

    H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

    size_t n_entries = ims.size();

    std::unordered_map<std::string,double> attrs{ {"nshells", ims[0].nshells},
                                                  {"nlayers", ims[0].nlayers} };
    for (auto& [key, value]: main_pars) { attrs[key] = value; }
    for (auto& [key, value]: ej_pars) { attrs[key] = value; }

    /// store binned images in h5
    std::vector<std::string> im_names {};
    for (auto & im : ims) {

        std::string name = " time=" + string_format("%.4e", im.time) +
                              " freq=" + string_format("%.4e", im.freq);
        Output::add2Dtable(im.binned_data, "image" + name, file);
//        Output::add2Dtable(im.interp_data, "interp_im" + name, file);
        Output::addVector(times, "times", file);
        Output::addVector(freqs, "freqs", file);

        attrs["flux" + name] = im.total_flux;
        attrs["xmin" + name] = im.xmin;
        attrs["xmax" + name] = im.xmax;
        attrs["ymin" + name] = im.ymin;
        attrs["ymax" + name] = im.ymax;
        attrs["xc" + name] = im.xc;
        attrs["yc" + name] = im.yc;
    }

    /// add atributes
    Output::addStrDbMap(attrs, file);

    file.close();
}

/// methods to compute radiation from a Blastwave


class EATS {

    void * m_params; /// parameters of the blast wave from which EATS is computed
    Vector & m_tburst; Vector & m_tt;
    Vector & m_r; Vector & m_theta; Vector & m_gam; Vector & m_bet;

//    Vector & m_freq_arr; Vector & out_spectrum; Vector & out_specturm_ssa;

//    VecVector & m_data;
    Vector m_mu{};
    double ctheta0 = -1.;
    std::unique_ptr<logger> p_log = nullptr;
    /// -------------------------------------------------------------------------------
    size_t ishell = 0, ilayer = 0;  size_t nlayers = 0; size_t & m_i_end_r; size_t ncells = 0;
    double theta_c_l = -1.;
    double theta_c_h = -1.;
    double theta_w = -1.;
    double theta_max=-1.;
    bool use_t_e = false;
    bool skymap_remove_mu = false;
    double z{}, d_l{}, nu_obs{}, t_obs{}, theta_obs{};
    double (* obsangle)(const double &, const double &, const double &) = nullptr;
    double (* im_xxs)( const double &, const double &, const double & ) = nullptr;
    double (* im_yys)( const double &, const double &, const double & ) = nullptr;
    Vector ttobs{}; // for PW eats only

    // --- for adaptive quad. method
    int nmax_phi=-1, nmax_theta=-1;
    double atol_theta=-1., atol_phi=-1.;
    int error = 0;
    double current_phi_hi=-1.;
    double current_phi_low=-1.;
    double current_theta_cone_hi=-1;
    double current_theta_cone_low=-1.;
    double cos_theta_obs=-1.;
    double sin_theta_obs=-1.;
    double phi = -1.;
//        double cos_phi = -1.;
    double theta=-1.;
//        double o_phi=-1, o_theta=-1., o_gam=-1, o_mu=-1, o_r=-1, o_flux=-1, o_theta_j=-1; // for map
//        double freq1=-1.0, freq2=-1.0;
//        size_t nfreq=0;

    METHODS_QUADRATURES method_quad{};

    //        METHOD_TAU method_tau{};
    bool counter_jet = true;
    int spread_method = 1; // TODO this is the only option!
    long nevals = 0; // counter for how many times the integrand was evaluated
    size_t n_tburst = 0;
    size_t i0_failed_elecctrons = 0;
    long n_fialed_electrons = 0;

    /// check if during the quadrature integration there was an error
    static int check_error(void *params) {
        auto *fp = (struct EATS *) params; // removing EATS_pars for simplicity
        return fp->error;
//        return 0;
    }
    /// find angle at which currently jet ends (after lateral expansion)
    static double find_jet_edge(double phi, double theta_obs, //double cos_th_obs, double sin_th_obs,
                                double th_con_hi, Vector & a_mu, Vector & a_thj, int N,
                                double (*obs_angle)( const double &, const double &, const double & )) {
        /*
         *
         * Return the (outer) edge of the jet section for a particular obs.
         *
         * phi: double
         *      phi-coordinate of jet along which to search
         * cos_theta_obs: double
         *      cos(theta_obs) cosine of observer angle
         * sin_theta_obs: double
         *      sin(theta_obs) sine of observer angle
         * th_con_hi: double
         *
         * a_mu: double array
         *      time-ordered array of mu values for this observation.
         *      mu = c * (t_em - t_obs) / R(t_em)
         *         = cos(th_obs)*cos(theta) + sin(theta_obs)*sin(theta)*cos(phi)
         * a_thj: double array
         *      time ordered array of jet-edge values. [from dynamic evolution]
         */

//        double cos_phi = cos(phi);
//        double mu = cos(th_con_hi) * cos_th_obs + sin(th_con_hi) * sin_th_obs * cos(phi); // cos(th_obs)*cos(theta) + sin(theta_obs)*sin(theta)*cos(phi)
        double mu = obs_angle( th_con_hi, phi, theta_obs );

        // in the mu_array -- find the index of the current 'mu'
        size_t ia = findIndex(mu, a_mu, N);
        // [ia,ib] contains the boundary of the shell, based on mu
//        std::cout <<AT << ' '<< a_mu[0]<< ", "<<a_mu[1]<<" ... "<<a_mu[N-2]<<", "<<a_mu[N-1]<<"\n";

        // *if the jet has not yet spread* = if theta[ia] from dynamic evolution < th_con_hi
        if(a_thj[ia] <= th_con_hi && th_con_hi <= a_thj[ia + 1])
            return th_con_hi;

        // *if the jet started to spread* start with the assumptions

        double theta_a, theta_b;
        // if current angle is larger then th_con_hi
        if(th_con_hi < a_thj[ia]) {
            //The jet is spreading
            theta_a = th_con_hi;   // set boundary
            theta_b = 0.5 * M_PI; // set boundary
        }
        else // if th_con_hi >= a_thj[ia]
        {
            //Guessed too far out!
            theta_a = 0.0;    // set boundary
            theta_b = th_con_hi; // set boundary
        }

        int i = 0;
        // similar to search sorted alborithim :: iteratively approach the 'theta' (theta_a) that correspods to mu
        while(theta_b - theta_a > 1.0e-5 && i < 100) {
            double theta_c_i = 0.5 * (theta_a + theta_b);          // middle of the two boundaries of the edge [theta_a ... edge ... ... theta_b]
//            mu = cos(theta_c_i) * cos_th_obs + sin(theta_c_i) * sin_th_obs * cos(phi);  // evaluate 'mu' at the current iteration of the search
            mu = obs_angle(theta_c_i, phi, theta_obs);

            ia = findIndex(mu, a_mu, N);     //
            if(theta_c_i < a_thj[ia])
                theta_a = theta_c_i;
            else
                theta_b = theta_c_i;
            i++;
        }

        //printf("iter: %d, th0=%.6lf, th=%.6lf\n", i, th_con_hi, theta_a);

        return theta_a;
    }
    /// function to be integrated for theta and phi
    static double integrand( double i_cos_theta, double i_phi, double & r, double & mu, double & gam, double & ctheta,
                             void* params ){

        auto * p_pars = (struct EATS *) params; // removing EATS_pars for simplicity
        if (p_pars->m_r[0] == 0.0 && p_pars->m_r[p_pars->m_i_end_r - 1] == 0.0){
            (*p_pars->p_log)(LOG_WARN, AT)
                    << " blast wave not evolved, flux=0 [ishell=" << p_pars->ishell << ", ilayer=" << p_pars->ilayer << "]\n";
            return 0.0;
        }

        /// current theta
        double a_theta = arccos(i_cos_theta);
        mu = p_pars->obsangle(a_theta, i_phi, p_pars->theta_obs);
        p_pars->theta = a_theta;

        /// find indeces of the BW dynamics that correspond to this EATS time
        size_t ia = findIndex(mu, p_pars->m_mu, p_pars->m_i_end_r);
        size_t ib = ia + 1;

        double t_e = interpSegLin(ia, ib, mu, p_pars->m_mu, p_pars->m_tburst);


        /// interpolate the time in burster frame that corresponds to the t_obs in observer frame
//        double t_e = interpSegLin(ia, ib, mu, p_pars->m_mu, p_pars->m_tburst);
        t_e = check_emission_time(t_e, mu, p_pars->t_obs, p_pars->m_mu, (int) p_pars->m_i_end_r);
        if ((t_e < 0.0)||(!std::isfinite(t_e))) {
            // REMOVING LOGGER
            (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
            return 0.;
        }
        /// -----------------------------------------
        double flux_dens = 0;
        p_pars->fluxFuncA(flux_dens, r, ctheta, a_theta, i_phi,
                          ia, ib, mu, t_e, p_pars->t_obs, p_pars->nu_obs, p_pars->m_params);

        /// ----------------------------------------
        gam = interpSegLog(ia, ib, t_e, p_pars->m_tburst, p_pars->m_gam);

#if 0
        /// Observed flux density evaluation (interpolate comoving spectrum)
        if (p_pars->m_method_rad == METHODS_RAD::icomovspec){
            Interp2d int_em(p_pars->m_freq_arr, r, p_pars->out_spectrum);
            Interp2d int_abs(p_pars->m_freq_arr, r, p_pars->out_specturm_ssa);
            Interp1d::METHODS_SYNCH mth = Interp1d::iLagrangeLinear;
            ///----
//            auto & mu_arr = m_data[BW::Q::imu];
            size_t ia = findIndex(mu, mu_arr, p_pars->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, mu_arr, tburst);
            t_e = check_emission_time(t_e, mu, p_pars->t_obs, mu_arr, (int) p_pars->nr);
            if (t_e < 0.0) {
                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }


            p_pars->fluxFuncA();


            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            double r = interpSegLog(ia, ib, t_e, tburst, r_arr);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r <= 0.0) || !std::isfinite(r)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R <= 0. Extend R grid (increasing R0, R1). "
                                             << " Current R grid us ["
                                             << r_arr[0] << ", "
                                             << r_arr[tburst.size() - 1] << "] "
                                             << "and tburst arr ["
                                             << tburst[0] << ", " << tburst[p_pars->nr - 1]
                                             << "] while the requried obs_time=" << p_pars->t_obs
                                             << "\n";
//                std::cerr << AT << "\n";
                return 0.;
            }
            double Gamma = interpSegLog(ia, ib, t_e, tburst, p_pars->m_gam);
            double betaSh = interpSegLog(ia, ib, t_e, tburst, p_pars->m_bet);
            // double GammaSh = ( Interp1d(m_data[BW::Q::iR], m_data[BW::Q::iGammaFsh] ) ).Interpolate(r, mth );
            /// evaluateShycnhrotronSpectrum Doppler factor
            double a = 1.0 - betaSh * mu; // beaming factor
            double delta_D = Gamma * a; // doppler factor

            /// evaluateShycnhrotronSpectrum the comoving obs. frequency from given one in obs. frame
            double nuprime = (1.0 + p_pars->z ) * p_pars->nu_obs * delta_D;
            size_t ia_nu = findIndex(nuprime, p_pars->m_freq_arr, p_pars->m_freq_arr.size());
            size_t ib_nu = ia_nu + 1;
            /// interpolate the emissivity and absorption coefficines
//                double em_prime = int_em.Interpolate(nuprime, r, mth);
            double em_prime = int_em.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
//                double abs_prime = int_abs.Interpolate(nuprime, r, mth);
            double abs_prime = int_abs.InterpolateBilinear(nuprime, r, ia_nu, ib_nu, ia, ib);
            /// convert to the laboratory frame
            double em_lab = em_prime / (delta_D * delta_D); // conversion of emissivity (see vanEerten+2010)
            double abs_lab = abs_prime * delta_D; // conversion of absorption (see vanEerten+2010)

            /// evaluateShycnhrotronSpectrum optical depth (for this shock radius and thickness are needed)
            double GammaShock = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::iGammaFsh]);
            double dr = interpSegLog(ia, ib, t_e, tburst, m_data[BW::Q::ithickness]);
            double dr_tau = EQS::shock_delta(r, GammaShock); // TODO this is added becasue in Johanneson Eq. I use ncells

            double beta_shock;
            switch (p_pars->method_shock_vel) {

                case isameAsBW:
                    beta_shock = EQS::Beta(Gamma);
                    break;
                case ishockVel:
//                double u = sqrt(GammaShock * GammaShock - 1.0);
//                double us = 4.0 * u * sqrt((u * u + 1) / (8. * u * u + 9.)); // from the blast wave velocity -> evaluateShycnhrotronSpectrum shock velocity
                    beta_shock = EQS::Beta(GammaShock);//us / sqrt(1. + us * us);
                    break;
            }
            double ashock = (1.0 - mu * beta_shock); // shock velocity beaming factor
            dr /= ashock; // TODO why is this here? What it means? Well.. Without it GRB LCs do not work!
            dr_tau /= ashock;
            double dtau = ElectronAndRadiaionBase::optical_depth(abs_lab,dr_tau, mu, beta_shock);
            double intensity = ElectronAndRadiaionBase::computeIntensity(em_lab, dtau,
                                                               p_mphys->getPars()->method_tau);
            double flux_dens = (intensity * r * r * dr); //* (1.0 + p_pars->z) / (2.0 * p_pars->d_l * p_pars->d_l);
            dFnu+=flux_dens;
            /// save the result in image
            double ctheta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::ictheta]);
            double theta = interpSegLin(ia, ib, t_e, tburst, m_data[BW::Q::itheta]);



        }
            /// Observed flux density evaluation (evaluateShycnhrotronSpectrum directly)
        else {

            size_t ia = findIndex(mu, p_pars->m_mu, p_pars->nr);
            size_t ib = ia + 1;
            /// interpolate the time in comobing frame that corresponds to the t_obs in observer frame
            double t_e = interpSegLin(ia, ib, mu, p_pars->m_mu, p_pars->m_tburst);
            t_e = check_emission_time(t_e, mu, p_pars->t_obs, p_pars->m_mu, (int) p_pars->nr);
            if (t_e < 0.0) {
//                auto & _r = m_data[BW::Q::iR];
//                auto & _mu = m_data[BW::Q::imu];
//                auto & _tt = m_data[BW::Q::itt];
//                auto & _tburst = m_data[BW::Q::itburst];
//                auto & _beta = m_data[BW::Q::ibeta];
//                (*p_pars->p_log)(LOG_ERR,AT) << "R      " << _r << "\n";
//                (*p_log)(LOG_ERR,AT) << "Mu     " << _mu << "\n";
//                (*p_log)(LOG_ERR,AT) << "tt     " << _tt << "\n";
//                (*p_log)(LOG_ERR,AT) << "tburst " << _tburst << "\n";
//                (*p_log)(LOG_ERR,AT) << "betaSh   " << _beta << "\n";


                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
//            std::cerr << AT  << "Error t_e < 0 = " << t_e << " Change R0/R1 parameters " << "\n";
                return 0.;
            }

            double flux_dens = 0;
            p_pars->fluxFuncA(flux_dens);


//        double R = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->r_arr, p_pars->nr);
            double R = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iR]);
            if (!std::isfinite(R)) {
                (*p_pars->p_log)(LOG_ERR,AT) << " R is NAN in integrand for microphysics" << "\n";
                // REMOVING LOGGER
//            std::cerr  << "R = " << R << "\n";
//            std::cout << " R = " << m_data[BW::Q::iR] << "\n";
//            std::cout << " Gamma= " << m_data[BW::Q::iGamma] << "\n";
//            std::cerr << AT << "\n";
                return 0.;
            }

            double rho = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::irho]);
            double Gamma = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::iGamma]);
            double GammaSh = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                          m_data[BW::Q::iGammaFsh]);
            double betaSh = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::ibeta]);
            double U_p = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iU_p]);
//        double M2    = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM2));
            double theta = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::itheta]);
            double rho2 = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::irho2]);
            double m2 = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iM2]);
            double frac = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                       m_data[BW::Q::iacc_frac]);
            double thick = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::ithickness]);
            double gm = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igm]);
            double gM = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igM]);
            double gc = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::igc]);
            double B = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst], m_data[BW::Q::iB]);
            double Theta = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                        m_data[BW::Q::iTheta]);
            double z_cool = interpSegLog(ia, ib, t_e, m_data[BW::Q::itburst],
                                         m_data[BW::Q::iz_cool]);

            if (rho < 0. || Gamma < 1. || U_p < 0. || theta <= 0. || rho2 < 0. || thick <= 0.) {
                (*p_pars->p_log)(LOG_ERR,AT) << " wrong value in interpolation to EATS surface  \n"
                                             << " R = " << R << "\n"
                                             << " rho = " << rho << "\n"
                                             << " Gamma = " << Gamma << "\n"
                                             << " U_p = " << U_p << "\n"
                                             << " theta = " << theta << "\n"
                                             << " rho2 = " << rho2 << "\n"
                                             << " thick = " << thick << "\n"
                                             << " t_e = " << t_e << "\n"
                                             << " Exiting...\n";
                exit(1);
            }

            dFnu = shock_synchrotron_flux_density(Gamma, GammaSh, m2, rho2,
                                                  frac, B, gm, gM, gc, Theta, z_cool,
                                                  t_e, mu, R, thick, thick, params);
//            dFnu*=(p_pars->d_l*p_pars->d_l*2.);
#if 1
            /* -- Reverse shock --- */
            double dFnu_rs = 0.0;
            if (p_pars->synch_rs){
                std::cout << AT << " Warning! Reverse shock is not finished!" << "\n";
                double rho4   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::irho4));
                double U_e3   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iU_e3));
                double rho3   = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::irho2));
                double thick3 = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::ithickness));
    //            double M3     = interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iM3));
                double gamma43= interpSegLog(ia, ib, t_e, p_pars->t_arr_burst, p_pars->dyn(BWDyn::iGamma43));
    //            double gammaAdi_rs = p_pars->p_eos->getGammaAdi(gamma43, EQS::Beta(gamma43));
    //            double nprime3 = 4.0 * Gamma * (rho4 / CGS::mppme);
    //            double nprime3 = p_pars->eq_rho2(Gamma, rho4 / CGS::mp, gammaAdi_rs); // TODO check if for RS here is Gamma!
                double nprime3 = rho3 / CGS::mp;
                /// computeSynchrotronEmissivityAbsorptionAnalytic the 'thickness' of the shock (emitting region)
    //            double dr_rs = thick3;

    //          // TODO check if for the shock velocity gamma43 has to be used!!! I think it is gam43! See gammaAdi calc.
                dFnu_rs = Observables::shock_synchrotron_flux_density(
                        Gamma, gamma43, rho3, U_e3, t_e, mu, R, thick3, thick3, params );
            }
            dFnu += dFnu_rs;
#endif
            if (dFnu == 0 || !std::isfinite(dFnu)) {
                // REMOVING LOGGER
                (*p_pars->p_log)(LOG_ERR,AT) << " flux density is zero ( dFnu = 0 )" << "\n";
            }



        }
#endif
//        p_pars->o_gam = gam;
//        p_pars->o_r = r;
//        p_pars->o_flux = flux_dens;
//        p_pars->o_theta_j = ctheta;
        return flux_dens;
    }
    static double costheta_integrand( double aomct, void* params ){
        auto * p_eats = (struct EATS *) params; // removing EATS_pars for simplicity
        p_eats->nevals+=1;
        if (p_eats->nevals > (long)1e7){
            (*p_eats->p_log)(LOG_ERR,AT) << " EATS adaptive integration frozen, "
                                                       " reached "<<p_eats->nevals<<" iterations\n";
            exit(1);
        }

        double act = 1 - aomct; // one minus cos theta 0 -> 'doppler_d' cos theta
        double r = 0., mu=0., gam=0., ctheta=0.;
        double integ = integrand(act, p_eats->phi, r, mu, gam, ctheta, params );
        return integ;
    }
    /// integral of the (costheta_integrand)dtheta
    static double phi_integrand( double a_phi, void* params ){
        double result;
        auto * p_eats = (struct EATS *) params; // removing EATS_pars for simplicity
//        auto & m_data = p_eats->m_data;
//        std::cout <<  p_eats->mu_arr[0]  << ' ' << p_eats->mu_arr[100] << "\n";
        p_eats->phi = a_phi;
//        p_eats->cos_phi = cos(a_phi);
        // INITIAL GUESSES FOR BOUNDARIES :: TO BE REFINED (based on spreading model)
        double theta_1 = p_eats->current_theta_cone_hi;  // shel upper theta boudnary
        double theta_0 = p_eats->current_theta_cone_low; // shel lower theta boundary
        double Dtheta = theta_1 - theta_0;             // m_size of the theta shell
        /// get the boundaries of the jet
        double th_0, th_1;
        if (p_eats->spread_method == 1){
            th_1 = find_jet_edge(a_phi, p_eats->theta_obs,// p_eats->cos_theta_obs, p_eats->sin_theta_obs,
                                 theta_1, p_eats->m_mu, p_eats->m_theta,//m_data[BW::Q::itheta],
                                 p_eats->m_i_end_r,
                                 p_eats->obsangle);
            double frac = theta_0 / theta_1;
            th_0 = frac * th_1;
            theta_0 = th_0;
            theta_1 = th_1;
            // assure that the inner and outer boundaries are within then 0.5Pi steradian
            if(theta_0 > 0.5*M_PI) theta_0 = 0.5*M_PI;
            if(theta_1 > 0.5*M_PI) theta_1 = 0.5*M_PI;
        } else {
            (*p_eats->p_log)(LOG_ERR,AT)  << " method not implemented. Exiting..." << "\n";
//            std :: cout << AT << "\n";
            exit(1);
        }
        /// For a given phi, integrate over 1 - cos(theta) (to avoid sin(theta->0) errors)
        double sin_half_theta_0 = std::sin(0.5 * theta_0);
        double sin_half_theta_1 = std::sin(0.5 * theta_1);
        double one_minus_cos_theta_0 = 2.0 * sin_half_theta_0 * sin_half_theta_0;
        double one_minus_cos_theta_1 = 2.0 * sin_half_theta_1 * sin_half_theta_1;

        int Neval;
        double err;
        switch (p_eats->method_quad) {

            case INT_TRAP_FIXED:
                result = trap(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              params, check_error);
                break;
            case INT_TRAP_ADAPT:
                result = trap_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                    check_error, nullptr, nullptr);
                break;
            case INT_SIMP_FIXED:
                result = simp(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              params, check_error);
                break;
            case INT_SIMP_ADAPT:
                result = simp_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                    check_error, nullptr, nullptr);
                break;
            case INT_ROMB_ADAPT:
                Neval = 0;
                err = 0;
                result = romb(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1, p_eats->nmax_theta,
                              p_eats->atol_theta, p_eats->rtol_theta, params,
                              &Neval, &err, 0, check_error, nullptr, nullptr);
                break;
            case INT_TRAP_NL:
                result = trapNL_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, nullptr, 0,
                                      check_error, nullptr, nullptr);
                break;
            case INT_HYBRID:
                result = hybrid_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                      check_error, nullptr, nullptr);
                break;
            case INT_CADRE:
                result = cadre_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                     p_eats->nmax_theta, p_eats->atol_theta,
                                     p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                     check_error, nullptr, nullptr);
                break;
            case INT_GK49_ADAPT:
                result = gk49_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                    p_eats->nmax_theta, p_eats->atol_theta,
                                    p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                    check_error);
                break;
            case INT_GK715_ADAPT:
                result = gk715_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                     p_eats->nmax_theta, p_eats->atol_theta,
                                     p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                     check_error);
                break;
            case INT_GK1021_ADAPT:
                result = gk1021_adapt(&costheta_integrand, one_minus_cos_theta_0, one_minus_cos_theta_1,
                                      p_eats->nmax_theta, p_eats->atol_theta,
                                      p_eats->rtol_theta, params, nullptr, nullptr, 0,
                                      check_error);
                break;
        }

        if(result != result || result < 0.0){
            // REMOVING LOGGER
            (*p_eats->p_log)(LOG_ERR,AT) << " phi_integrand failed. Nan or negative. Result=" << result
                                         << " for t_obs=" << p_eats->t_obs
                                         << " theta_lo=" << theta_0
                                         << " theta_hi=" << theta_1
                                         << " phi=" << p_eats->phi
                                         << "\n"
                                         << " Exiting...";
//            std::cerr << AT << "\n";
            return 0.;
        }
//        printf("   a_phi: %.6lf (%.6le)\n", a_phi, result);
        if (p_eats->nevals > (long)1e7){
            (*p_eats->p_log)(LOG_ERR,AT) << " EATS adaptive integration frozen, "
                                            " reached "<<p_eats->nevals<<" iterations\n";
            exit(1);
        }
        return result;

    }
    /// integral of the (phi_integrand)dphi
    double integrate_theta_phi(){

        auto * p_eats = (struct EATS *) this; // this "this" i do not understand myself...
        p_eats->nevals = 0;
        double atol = p_eats->atol_theta;
        // check if the parameters are set
        if(p_eats->nmax_theta < 0 || p_eats->nmax_phi < 0 || p_eats->rtol_phi < 0
           || p_eats->rtol_theta < 0 || p_eats->atol_theta < 0 || p_eats->atol_phi < 0 ){
            (*p_eats->p_log)(LOG_ERR,AT) << " one of the tolerance parameters for adaptive quad. is not set (< 0). Given:\n"
                                         << " nmax_theta=" << p_eats->nmax_theta
                                         << " nmax_phi=" << p_eats->nmax_phi
                                         << " rtol_phi=" << p_eats->rtol_phi
                                         << " rtol_theta=" << p_eats->rtol_theta
                                         << " atol_theta=" << p_eats->atol_theta
                                         << " atol_phi=" << p_eats->atol_phi << "\n"
                                         << " Exiting...\n";
//            std::cerr << AT << "\n";
            exit(1);
        }
        double phi_0 = p_eats->current_phi_low;//0.0;
        double phi_1 = p_eats->current_phi_hi;//M_PI;
        double Fcoeff = CGS::cgs2mJy / (4*M_PI * p_eats->d_l* p_eats->d_l);

        double result;
        double phi_a;
        switch (p_eats->method_quad) {

            case INT_TRAP_FIXED:
                result = 2 * Fcoeff * trap(&phi_integrand, phi_0, phi_1,
                                           p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_TRAP_ADAPT:
                result = 2 * Fcoeff * trap_adapt(&phi_integrand, phi_0, phi_1,
                                                 p_eats->nmax_phi, atol,
                                                 p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                 nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_SIMP_FIXED:
                result = 2 * Fcoeff * simp(&phi_integrand, phi_0, phi_1,
                                           p_eats->nmax_phi, p_eats, check_error);
                break;
            case INT_SIMP_ADAPT:
                result = 2 * Fcoeff * simp_adapt(&phi_integrand, phi_0, phi_1,
                                                 p_eats->nmax_phi, atol,
                                                 p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                 nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_ROMB_ADAPT:
                phi_a = phi_0 + 0.5*(phi_1-phi_0);
                result = 2 * Fcoeff * romb(&phi_integrand, phi_0, phi_a,
                                           p_eats->nmax_phi, atol,
                                           p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                                           check_error, nullptr, nullptr);
                result += 2 * Fcoeff * romb(&phi_integrand, phi_a, phi_1,
                                            p_eats->nmax_phi,
                                            (atol + p_eats->rtol_phi * result),
                                            p_eats->rtol_phi, p_eats, nullptr, nullptr, 0,
                                            check_error, nullptr, nullptr);
                break;
            case INT_TRAP_NL:
                result = 2 * Fcoeff * trapNL_adapt(&phi_integrand, phi_0, phi_1,
                                                   p_eats->nmax_phi, atol,
                                                   p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                   nullptr, 0, check_error, nullptr, nullptr);
                break;
            case INT_HYBRID:
                result = 2 * Fcoeff * hybrid_adapt(&phi_integrand, phi_0, phi_1,
                                                   p_eats->nmax_phi, atol,
                                                   p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                   0, check_error, nullptr, nullptr);
                break;
            case INT_CADRE:
                result = 2 * Fcoeff * cadre_adapt(&phi_integrand, phi_0, phi_1,
                                                  p_eats->nmax_phi, atol,
                                                  p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                  0, check_error, nullptr, nullptr);
                break;
            case INT_GK49_ADAPT:
                result = 2 * Fcoeff * gk49_adapt(&phi_integrand, phi_0, phi_1,
                                                 p_eats->nmax_phi, atol,
                                                 p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                 0, check_error);
                break;
            case INT_GK715_ADAPT:
                result = 2 * Fcoeff * gk715_adapt(&phi_integrand, phi_0, phi_1,
                                                  p_eats->nmax_phi, atol,
                                                  p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                  0, check_error);
                break;
            case INT_GK1021_ADAPT:
                result = 2 * Fcoeff * gk1021_adapt(&phi_integrand, phi_0, phi_1,
                                                   p_eats->nmax_phi, atol,
                                                   p_eats->rtol_phi, p_eats, nullptr, nullptr,
                                                   0, check_error);
                break;
        }
        return result;
    }

public:

    double rtol_theta=-1., rtol_phi=-1.;

    /// ---------------------------------------------------------------------------

    EATS(Vector & tburst, Vector & tt, Vector & r, Vector & theta, Vector & m_gam, Vector & m_bet,
//             Vector & freq_arr, Vector & synch_em, Vector & synch_abs,
             size_t & i_end_r, size_t ish, size_t il, int loglevel, void * params)
            : m_tburst(tburst), m_tt(tt),  m_r(r), m_theta(theta), m_gam(m_gam), m_bet(m_bet),
//              m_freq_arr(freq_arr), out_spectrum(synch_em), out_specturm_ssa(synch_abs),
              m_i_end_r(i_end_r),
              m_params(params) {

        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "EATS");

        ishell = ish;
        ilayer= il;
        m_i_end_r = i_end_r;
    }

    void (* fluxFuncPW)(
            double & flux_dens, double & tau_comp, double & tau_BH, double & tau_bf,
            double r, double & ctheta, double theta, double phi,
            size_t ia, size_t ib, double ta, double tb, double mu, double t_obs, double nu_obs, void * params
    ) = nullptr;
//        void (* funcOptDepth)(
//                double & frac, double ctheta, double r,
//                double phi, double theta,
//                double phi_obs, double theta_obs, double r_obs,
//                double mu, double time, double freq, void * params
//                ) = nullptr;
    void (* fluxFuncA)(
            double & flux_dens, double & r, double & ctheta, double theta, double phi,
            size_t ia, size_t ib, double mu, double t_e, double t_obs, double nu_obs, void * params
    ) = nullptr;

    /// -----------------------------------------------------------------------------

    void setEatsPars(StrDbMap & pars, StrStrMap & opts, size_t n_layers, size_t ncells_, double ctheta_0,
                     double theta_c_l_, double theta_c_h_, double theta_w_, double theta_max_){
        nlayers=n_layers;
        ctheta0=ctheta_0;
        theta_c_l=theta_c_l_;
        theta_c_h=theta_c_h_;
        theta_w=theta_w_;
        theta_max=theta_max_;
        ncells=ncells_;

        // set parameters
        nmax_phi = (int)getDoublePar("nmax_phi", pars, AT, p_log,1000, false);//pars.at("nmax_phi");
        nmax_theta = (int)getDoublePar("nmax_theta", pars, AT, p_log,1000, false);//pars.at("nmax_theta");
//        rtol_theta = getDoublePar("rtol_theta", pars, AT, p_log,1e-2, false);//pars.at("rtol_theta");
//        rtol_phi = getDoublePar("rtol_phi", pars, AT, p_log,1e-2, false);//pars.at("rtol_phi");
//        atol_theta = getDoublePar("atol_theta", pars, AT, p_log,1e-2, false);//pars.at("atol_theta");
//        atol_phi = getDoublePar("atol_phi", pars, AT, p_log,1e-2, false);//pars.at("atol_phi");
        rtol_phi = getDoublePar("rtol_phi", pars, AT, p_log,1e-15, true);;
        rtol_theta = getDoublePar("rtol_theta", pars, AT, p_log,1e-15, true);;
        theta_obs = getDoublePar("theta_obs", pars, AT, p_log,-1, true);//pars.at("theta_obs");
        d_l = getDoublePar("d_l", pars, AT, p_log,-1, true);//pars.at("d_l");
        z = getDoublePar("z", pars, AT, p_log,-1, true);//pars.at("z");

        counter_jet = getBoolOpt("counter_jet", opts, AT, p_log,true, false);

        // set options
        std::string opt;

        opt = "method_quad";
        METHODS_QUADRATURES methodsQuadratures;
        if ( opts.find(opt) == opts.end() ) {
            (*p_log)(LOG_ERR,AT) << " Option for '" << opt << "' is not set. Using default value.\n";
            methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
        }
        else{
            if(opts.at(opt) == "CADRE")
                methodsQuadratures = METHODS_QUADRATURES::INT_CADRE;
            else if(opts.at(opt) == "TRAP_FIXED")
                methodsQuadratures = METHODS_QUADRATURES::INT_TRAP_FIXED;
            else{
                (*p_log)(LOG_ERR,AT) << " option for: " << opt
                                     <<" given: " << opts.at(opt)
                                     << " is not recognized. "
                                     << "Possible options: "
                                     << " CADRE " << " TRAP_FIXED " << "\n";
                exit(1);
            }
        }
        method_quad = methodsQuadratures;

        /// set synchrotron parameters
//            p_mphys->setPars(pars, opts);
        skymap_remove_mu = getBoolOpt("skymap_remove_mu", opts, AT, p_log,true, true);

    }


    /// eval light curve using Adapitve or Piece-Wise EATS method
    void evalLightCurve(Vector & out, EjectaID2::STUCT_TYPE m_method_eats, Vector & times, Vector & freqs ){
        VecVector empty{};
        for (size_t it = 0; it < times.size(); it++) {
            if (m_method_eats == EjectaID2::STUCT_TYPE::ipiecewise)
                out[it] = evalSkyMapPW(empty, times[it], freqs[it], 0);
            else{
                double atol = out[it] * rtol_theta / (double)nlayers;
                out[it] += evalFluxDensA(times[it], freqs[it], atol);
            }
        }
    }

    /// ----------------------------------------------------------------------------------------------
    /// evaluate flux density using adaptive integrator
    double evalFluxDensA(double t_obs, double nu_obs, double atol) {
        /// check if given observer time is at least within burster time
        if (m_i_end_r > 0 && (t_obs < m_tburst[0] || t_obs> m_tburst[m_i_end_r-1])){
            (*p_log)(LOG_ERR,AT) << " t_obs="<<t_obs<<" is not in tburst["<<m_tburst[0]
                                 <<", "<<m_tburst[m_i_end_r-1]<<" Extend tburst grid or shorten tobs grid. \n";
            exit(1);
        }
        double fluxdens = 0.;
        parsPars(t_obs, nu_obs, theta_c_l, theta_c_h,
                 0., M_PI, obsAngle);
        check_pars();
        double Fcoeff = CGS::cgs2mJy / (4.0 * M_PI * d_l * d_l); // result will be in mJy
//        rtol_theta = rtol;
//        rtol_phi = rtol;
        atol_theta = atol/(2.*Fcoeff*M_PI);// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
        atol_phi = atol/(2.*Fcoeff*M_PI);//  / (2.0 * Fcoeff);
        fluxdens += integrate_theta_phi(); // 2. because Integ_0^pi (not 2pi)
        if (counter_jet){
            parsPars(t_obs, nu_obs, theta_c_l, theta_c_h,
                     0., M_PI, obsAngleCJ);
            check_pars();
//            atol_theta = atol/(2.*Fcoeff*M_PI);// / M_PI / (2.0 * Fcoeff * M_PI);  // correct the atol to the scale
//            atol_phi = atol/(2.*Fcoeff*M_PI);//  / (2.0 * Fcoeff);
            fluxdens += integrate_theta_phi();
        }
        return fluxdens;
    }

    /// evaluate intensity/flux density distribution using piece-wise summation
    double evalSkyMapPW(VecVector & out, double t_obs, double freq_obs, size_t offset){
        /// check if given observer time is at least within burster time
        if (m_i_end_r > 0 && (t_obs < m_tburst[0] || t_obs> m_tburst[m_i_end_r-1])){
            (*p_log)(LOG_ERR,AT) << " t_obs="<<t_obs<<" is not in tburst["<<m_tburst[0]
                                 <<", "<<m_tburst[m_i_end_r-1]<<" Extend tburst grid or shorten tobs grid. \n";
            exit(1);
        }
        double flux_pj=0., flux_cj=0.;
        flux_pj = evalSkyMapPW(out, t_obs, freq_obs, offset, obsAngle, imageXXs, imageYYs);
        if (counter_jet) // p_eats->counter_jet
            flux_cj = evalSkyMapPW(out, t_obs, freq_obs, ncells + offset, obsAngleCJ, imageXXsCJ, imageYYsCJ);
        return flux_pj + flux_cj;
    }

    /// evaluate intensity/flux density distribution using adaptive summation
    double evalSkyMapA(VecVector & out, double t_obs, double freq_obs, size_t il, size_t offset, size_t cil, size_t ncells){
        /// check if given observer time is at least within burster time
        if (m_i_end_r > 0 && (t_obs < m_tburst[0] || t_obs> m_tburst[m_i_end_r-1])){
            (*p_log)(LOG_ERR,AT) << " t_obs="<<t_obs<<" is not in tburst["<<m_tburst[0]
                                 <<", "<<m_tburst[m_i_end_r-1]<<" Extend tburst grid or shorten tobs grid. \n";
            exit(1);
        }
        /// computeSynchrotronEmissivityAbsorptionAnalytic image for primary jet and counter jet
        double int_pj=0., int_cj=0.;
        int_pj = evalSkyMapA(out, t_obs, freq_obs, il, offset, cil, obsAngle, imageXXs, imageYYs);
        if (counter_jet) // p_eats->counter_jet
            int_cj = evalSkyMapA(out, t_obs, freq_obs, il, offset + ncells, cil, obsAngleCJ, imageXXsCJ, imageYYsCJ);
        return (int_pj+int_cj);
    }

    /// evaluate jet extend for further grid allocation
    void findJetEdgeA(double & th_b_min, double & th_b_max, double t_obs, double freq_obs, double th_l_prev){
        /// check if given observer time is at least within burster time
        if (m_i_end_r > 0 && (t_obs < m_tburst[0] || t_obs> m_tburst[m_i_end_r-1])){
            (*p_log)(LOG_ERR,AT) << " t_obs="<<t_obs<<" is not in tburst["<<m_tburst[0]
                                 <<", "<<m_tburst[m_i_end_r-1]<<" Extend tburst grid or shorten tobs grid. \n";
            exit(1);
        }
        double th_b_min_pj=0., th_b_max_pj=0.;
        double th_b_min_cj=0., th_b_max_cj=0.;
        findJetEdge(th_b_min_pj, th_b_max_pj, t_obs, freq_obs, th_l_prev, obsAngle);
        if (counter_jet) // p_eats->counter_jet
            findJetEdge(th_b_min_cj, th_b_max_cj, t_obs, freq_obs, th_l_prev, obsAngleCJ);
        th_b_min = std::min(th_b_min_pj,th_b_min_cj);
        th_b_max = std::max(th_b_max_pj,th_b_max_cj);
        if (th_b_min < th_l_prev)
            th_b_min = th_l_prev;
        if (th_b_max > CGS::pi/2.){
            (*p_log)(LOG_ERR,AT) << " th_b_max > CGS::pi/2.\n";
            exit(1);
        }
    }

private:

    void parsPars(double t_obs_, double nu_obs_,
                  double theta_cone_low, double theta_cone_hi, double phi_low, double phi_hi,
                  double (*obs_angle)( const double &, const double &, const double & )){
//        auto & m_data = p_eats->m_data;
        nu_obs = nu_obs_;
        t_obs = t_obs_;
        // -- settings
        error = 0;
        current_phi_hi = phi_hi;
        current_phi_low = phi_low;
        current_theta_cone_hi = theta_cone_hi;
        current_theta_cone_low = theta_cone_low;
        cos_theta_obs = std::cos(theta_obs);
        sin_theta_obs = std::sin(theta_obs);
        obsangle = obs_angle;
        // ---
        if (m_tburst.empty()){
            (*p_log)(LOG_ERR,AT) << " m_tburst is empty\n";
            exit(1);
        }
//            nr = m_tburst.size();

        m_mu.resize(m_tburst.size());
//        p_eats->nr = p_pars->nr; // removing EATS_pars for simplicity
        for (size_t i = 0; i < m_i_end_r; i++) {
            m_mu[i] = (m_tburst[i] - t_obs / (1.0 + z)) / m_r[i] * CGS::c;
//                if (m_mu[i] == 0.){
//                    printf("error");
//                    exit(1);
//                }
        }
//            std::cout <<AT << ' '<< m_mu[0]<< ", "<<m_mu[1]<<" ... "<<m_mu[N-2]<<", "<<m_mu[N-1]<<"\n";

        if(m_mu[m_i_end_r - 1] < 1. ){
//                std::cout << " tburst = " << m_tburst << "\n";
//                std::cout << " r      = "<< m_r << "\n";
//                std::cout << " mu     = "<<m_mu<<"\n";
            if (m_i_end_r == m_mu.size()-1)
                (*p_log)(LOG_WARN,AT) << " mu[-1]=" <<m_mu[m_i_end_r - 1] << " < 1 (expected >1) "
                                      << " tobs=" << t_obs
                                      << " tbutst[-1]=" << m_tburst[m_i_end_r - 1]
                                      << " R[nr-1]=" << m_r[m_i_end_r - 1]
                                      << "\n";
            int x = 1;
        }
        if(m_mu[0] > 1. ){
            (*p_log)(LOG_WARN,AT) << " mu[0]=" << m_mu[0] << " > 1 (expected <1) "
                                  << " tobs=" << t_obs
                                  << " tbutst[0]=" << m_tburst[0]
                                  << " R[0]=" << m_r[0]
                                  << "\n";
        }
        if(m_mu[0]==m_mu[m_mu.size()-1] and m_mu[0]==0){
            (*p_log)(LOG_ERR,AT) << " m_mu[0]=m_mu[-1]=0 for [il="
                                 <<ilayer<<"] m_i_end_r="<<m_i_end_r<<"\n";
        }
        if (ttobs.empty())
            ttobs.resize( m_r.size(), std::numeric_limits<double>::max() );

        //        std::cerr << AT
//                  << "Observables for Gamma0=[" << (p_pars->p_dyn)->getData()[(*p_pars->p_dyn).Q::iGamma][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iGamma][p_pars->nr - 1] << "] "
//                  << " R0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::iR][p_pars->nr - 1] << "] "
//                  << " time0=[" << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][0] << ", " << (p_pars->p_dyn)->getData()[(p_pars->p_dyn)->Q::itburst][p_pars->nr - 1] << "] "
//                  << " mu0=[" << p_pars->mu_arr[0] << ", " << p_pars->mu_arr[p_pars->nr-1]<<"] "
//                  <<"\n";

    }

    void check_pars() const {
        if ((nu_obs <= 0.) || (z < 0) || (z > 10) || (d_l < 1) || (t_obs < 1) || (theta_c_h < 0) || (theta_c_l < 0)){
            (*p_log)(LOG_ERR,AT) << " error in input parameters"
                                 << " nu_obs="<<nu_obs
                                 << " z="<<z
                                 << " d_l="<<d_l
                                 << " t_obs="<<t_obs
                                 << " theta_obs="<<theta_obs
                                 << " theta_cone_low="<<current_theta_cone_low
                                 << " theta_cone_hi="<<current_theta_cone_hi
                                 << " phi_low="<<current_phi_low
                                 << " phi_hi="<<current_phi_hi
                                 << " theta_c_h="<<theta_c_h
                                 << " theta_c_l="<<theta_c_l
                                 << "\n Exiting... \n";
            exit(1);
        }
    }

    bool evalEATSindexes(size_t & ia, size_t & ib,
                         double t_obs, double obs_angle, double ctheta_cell, double phi_cell,
                         double (*obs_angle_func)( const double &, const double &, const double & )){
        check_pars();

//        if (p_pars->theta_obs < 0.){
//            (*p_log)(LOG_ERR,AT)<<" theta_obs="<<p_pars->theta_obs<<"\n";
//            exit(1);
//        }
        if (m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT)<<" p_pars->m_i_end_r="<<m_i_end_r<<"\n";
            return false;
        }

//        double ctheta_cell = p_pars->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
        double mu = obs_angle_func(ctheta_cell, phi_cell, obs_angle);
        for (size_t i_ = 0; i_ < m_i_end_r; i_++) {
            ttobs[i_] = m_tt[i_] + m_r[i_] / CGS::c * (1.0 - mu);
        }
        /// check if req. obs time is outside of the evolved times (throw error)
        if (t_obs < ttobs[0]) {
            (*p_log)(LOG_WARN, AT) << "t_obs="<<t_obs<<" < ttobs_arr[0]="<<ttobs[0]<<". Extend ttobs.\n";

//            " time grid starts too late "
//                                  << " t_grid[0]=" << p_pars->ttobs[0] << " while requested obs.time=" << t_obs << "\n"
//                                  << " extend the grid to earlier time or request tobs at later times\n"
//                                  << " Exiting...\n";
            exit(1);
        }
        if ((t_obs > ttobs[m_i_end_r - 1])) {
            if (m_i_end_r == m_mu.size()-1)
                (*p_log)(LOG_WARN, AT) << "t_obs="<<t_obs<<" > ttobs_arr[i_end="<<m_i_end_r - 1<<"]="
                                       <<ttobs[m_i_end_r - 1]<<". Extend ttobs.\n";
//
//            << " time grid ends too early. "
//                                  << " t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                  << " while requested obs.time=" << t_obs
//                                  << " extend the grid to later time or request tobs at earlier times\n";
//                    std::cout << ttobs << std::endl;
//            exit(1);
            return false;
        }
//        else if ((t_obs > p_pars->ttobs[p_pars->m_i_end_r - 1])) {
//            (*p_log)(LOG_WARN, AT) << " time grid was shorten to i=" << p_pars->m_i_end_r
//                                   << " from nr=" << p_pars->m_i_end_r
//                                   << " and now ends at t_grid[i_end_r-1]=" << p_pars->ttobs[p_pars->m_i_end_r - 1]
//                                   << " while t_obs=" << t_obs << "\n";
//            return;
//        }
        /// locate closest evolution points to the requested obs. time
        ia = findIndex(t_obs, ttobs, ttobs.size());
        if (ia >= m_i_end_r - 1)
            return false; // ??
        ib = ia + 1;

        return true;
        ///
//        return r;
    }

    double evalSkyMapPW(VecVector & out, double obs_time, double obs_freq, size_t offset,
                        double (*obs_angle)( const double &, const double &, const double & ),
                        double (*im_xxs)( const double &, const double &, const double & ),
                        double (*im_yys)( const double &, const double &, const double & )){

        /// out is [i_vn][ish][itheta_iphi]
        bool save_im = true;
        if ( out.empty() )
            save_im = false;

        if ((m_r[0] == 0.) && (m_gam[0] == 0.)){
            (*p_log)(LOG_WARN,AT) << " [ishell=" << ishell << " ilayer="<<ilayer << "] "
                                  << " R[0]=0. Seems not evolved -> returning isEmpty image." << "\n";
            return 0;
        }
        if (m_i_end_r == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_i_end_r = 0\n";
            exit(1);
        }
        if (m_tt[0] == 0 and m_tt[m_i_end_r] == 0){
            (*p_log)(LOG_ERR,AT) << "p_pars->m_tt = 0\n";
            exit(1);
        }

        parsPars(obs_time, obs_freq, 0., 0., 0., 0., obs_angle);
        check_pars();
        size_t cil = EjectaID2::CellsInLayer(ilayer);

        /// allocate memory for quantities that help
        Vector mu (m_i_end_r, std::numeric_limits<double>::max());
        std::vector<size_t> ia(m_i_end_r, 0);
        std::vector<size_t> ib(m_i_end_r, 0);
        Vector ta(m_i_end_r, 0);
        Vector tb(m_i_end_r, 0);
        Vector r(m_i_end_r, 0);
        Vector fluxes(m_i_end_r, 0);

        /// Find the region for EATS interpoaltion
        size_t nskipped_h = 0;
        size_t nskipped_r = 0;
        size_t nskipped_p = 0;
        for (size_t i = 0; i < cil; i++) {
            double cphi = (double)i * 2.0 * M_PI / (double)cil;
            double ctheta_cell = ctheta0;
            mu[i] = obs_angle(ctheta_cell, cphi, theta_obs);
            /// fill-in the time on EATS plane
            for (size_t i_ = 0; i_ < m_i_end_r; i_++) {
                ttobs[i_] = m_tt[i_] + m_r[i_] / CGS::c * (1.0 - mu[i]);
            }
            /// check if req. obs time is outside of the evolved times (throw error)
            if (obs_time < ttobs[0]) {
                (*p_log)(LOG_ERR, AT) << " time grid starts too late "
                                      << " t_grid[0]=" << ttobs[0] << " while requested obs.time=" << obs_time << "\n"
                                      << " extend the grid to earlier time or request tobs at later times\n"
                                      << " Exiting...\n";
                exit(1);
            }
            else if ((obs_time > ttobs[m_i_end_r - 1])) {
                nskipped_h++;
                mu[i] = std::numeric_limits<double>::max();
                continue;
            }
            ///
            int guess = i > 0 ? (int)ia[i-1] : (int)(ttobs.size()/2);
            ia[i] = findClosestIndex(obs_time, ttobs, guess);
            if (ia[i] >= m_i_end_r - 1) {
                nskipped_p++;
                mu[i] = std::numeric_limits<double>::max();
                continue; // ??
            }
            ib[i] = ia[i] + 1;
            ta[i] = ttobs[ia[i]];
            tb[i] = ttobs[ib[i]];
            /// interpolate the exact radial position of the blast that corresponds to the req. obs time
            r[i] = interpSegLog(ia[i], ib[i], obs_time, ttobs, m_r);
            //  double r = ( Interp1d(ttobs, m_data[BW::Q::iR] ) ).Interpolate(t_obs, mth );
            if ((r[i] <= 0.0) || (!std::isfinite(r[i]))) {
                nskipped_r++;
                mu[i] = std::numeric_limits<double>::max();
                continue;
            }
        }

        auto * _params = m_params;
        double phi_cell=0., ctheta_cell=0., flux_dens=0., ctheta=0.;
        double tau_comp=-1., tau_bf=-1., tau_BH=-1.;
        double tot_flux = 0.;
        /// Perform EATS integration
        for (size_t i = 0; i < cil; i++) {
            if (mu[i] == std::numeric_limits<double>::max())
                continue;

            phi_cell = (double) i * 2.0 * M_PI / (double) cil;
            ctheta_cell = ctheta0;

            fluxFuncPW(flux_dens, tau_comp, tau_BH, tau_bf, r[i], ctheta, ctheta_cell, phi_cell,
                       ia[i], ib[i], ta[i], tb[i], mu[i], obs_time, obs_freq, _params);

            if (!std::isfinite(flux_dens) || flux_dens < 0){
                (*p_log)(LOG_ERR,AT) << " flux_dens="<<flux_dens<<"\n";
                exit(1);
            }
            /// save data in container
            if (save_im) {
                out[IMG::Q::iintens][offset + i] = flux_dens / (r[i] * r[i]) * CGS::cgs2mJy;
                out[IMG::Q::ixr][offset + i] = r[i] * im_xxs(ctheta, phi_cell, theta_obs);
                out[IMG::Q::iyr][offset + i] = r[i] * im_yys(ctheta, phi_cell, theta_obs);
                out[IMG::Q::ir][offset + i] = r[i];
                out[IMG::Q::ictheta][offset + i] = ctheta_cell;
                out[IMG::Q::icphi][offset + i] = phi_cell;
                out[IMG::Q::imu][offset + i] = mu[i];

                if (!skymap_remove_mu)
                    out[IMG::Q::iintens][offset + i] /= std::abs(mu[i]);

                /// convert to mas
                out[IMG::Q::ixr][offset + i] *= (CGS::rad2mas / d_l);
                out[IMG::Q::iyr][offset + i] *= (CGS::rad2mas / d_l);
                out[IMG::Q::iintens][offset + i] *= (d_l * d_l / CGS::rad2mas / CGS::rad2mas);
            }
            tot_flux += flux_dens;
        }

        if (nskipped_h == cil || nskipped_p == cil || nskipped_r == cil || (nskipped_p+nskipped_r+nskipped_h) == cil){
            (*p_log)(LOG_WARN,AT)
                    << " N(obs_time > p_pars->ttobs[p_pars->m_i_end_r - 1]) = "<<nskipped_h<<"/"<<cil
                    << " N(ia[i] >= p_pars->m_i_end_r - 1) = "<<nskipped_p<<"/"<<cil
                    << " N((r[i] <= 0.0)||(isnan(r[i]))) = "<<nskipped_r<<"/"<<cil
                    << " \n";
//            (*p_log)(LOG_WARN,AT)<<" try extending tburst() grid\n";
        }

        return (tot_flux * CGS::cgs2mJy);
    }

    double evalSkyMapA(VecVector & out, double obs_time, double obs_freq,
                       size_t il, size_t offset, size_t cil,
                       double (*obs_angle)( const double &, const double &, const double & ),
                       double (*im_xxs)( const double &, const double &, const double & ),
                       double (*im_yys)( const double &, const double &, const double & )) {

        if (out.size() == 0) {
            (*p_log)(LOG_ERR, AT) << " image is empty\n";
            exit(1);
        }
        if (il < 0) {
            (*p_log)(LOG_ERR, AT) << " il < 0\n";
            exit(1);
        }
        parsPars(obs_time, obs_freq, theta_c_l, theta_c_h,
                 0., 2. * M_PI, obs_angle);
        check_pars();

        double summed_intensity = 0;
        double Fcoeff = cgs2mJy / (4. * M_PI * d_l * d_l);

//        size_t cill = out[IMG::Q::ictheta][il].size() / 2;
        size_t _i = 0;
        double th_b = 0; double th_a = 0;
        for (size_t iphi = 0; iphi < cil; iphi++) {
            double cth = out[IMG::Q::ictheta][offset + iphi];
            double cphi = out[IMG::Q::icphi][offset + iphi];

            // if jet is spreading, compute the upper boundary of the jet
            th_b = find_jet_edge(cphi, theta_obs, //p_pars->cos_theta_obs, p_pars->sin_theta_obs,
                                         theta_c_h, m_mu, m_theta,//m_data[BW::Q::itheta],
                                         (int) m_i_end_r, obsangle);
            th_a = (theta_c_l / theta_c_h) * th_b; // ???
            if (cth < th_a || cth > th_b)
                continue;
            _i += 1;

            // compute intensity
            double r = 0., mu = 0., gam = 0., ctheta_bw = 0.;
            double intensity = integrand(std::cos(cth), cphi, r, mu, gam, ctheta_bw, this);

            double x = r * im_xxs(cth, cphi, theta_obs);
            double y = r * im_yys(cth, cphi, theta_obs);

            out[IMG::Q::iintens][offset + iphi] = intensity;// / (r * r * std::abs(mu));//* CGS::cgs2mJy;
            out[IMG::Q::ir][offset + iphi] = r;
            out[IMG::Q::ixr][offset + iphi] = x;
            out[IMG::Q::iyr][offset + iphi] = y;
            out[IMG::Q::imu][offset + iphi] = mu;
            summed_intensity += out[IMG::Q::iintens][offset + iphi];

//            if (!p_pars->skymap_remove_mu)
//                out[IMG::Q::iintens][p_pars->ishell][offset + iphi] /= std::abs(mu);

            /// convert to mas
            out[IMG::Q::ixr][offset + iphi] *= (CGS::rad2mas / d_l);
            out[IMG::Q::iyr][offset + iphi] *= (CGS::rad2mas / d_l);
            out[IMG::Q::iintens][offset + iphi] *= (d_l * d_l / CGS::rad2mas / CGS::rad2mas);
        }
//        if (_i == 0){
//            (*p_log)(LOG_ERR,AT) << " error, no intensity collected from il="<<il
//            <<" ctheta1="<<out[IMG::Q::ictheta][il][offset + cil - 1]
//            <<" th_b_last="<<th_b
//            <<"\n";
//        }

        return summed_intensity;
    }

    void findJetEdge(double & th_b_min, double & th_b_max, double obs_time, double obs_freq, double theta_l_prev,
                     double (*obs_angle)( const double &, const double &, const double & )){
        parsPars(obs_time, obs_freq, theta_c_l, theta_c_h, 0., 2. * M_PI, obs_angle);
        check_pars();

        size_t nphi_ = 100;

        double phi0 = 0.;
        double phi1 = 2. * M_PI;

        for (size_t iphi = 0; iphi < nphi_; iphi++) {
            double dphi = (phi1 - phi0) / (double) nphi_;
            double cphi = phi0 + (double) iphi * dphi;
            /// if jet is spreading, compute the upper boundary of the jet
            double th_b = find_jet_edge(cphi, theta_obs, //p_pars->cos_theta_obs, p_pars->sin_theta_obs,
                                                theta_c_h, m_mu, m_theta,//m_data[BW::Q::itheta],
                                                (int) m_i_end_r, obsangle);
            if (th_b > th_b_max) th_b_max = th_b;
            double th_a = ( theta_c_l / theta_c_h ) * th_b;
            if (th_b < th_b_min) th_b_min = th_b;
        }
    }

};


#endif //SRC_EATS_H
