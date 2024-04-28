//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_IMAGE_H
#define SRC_IMAGE_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "utilitites/interpolators.h"
#include "blastwave/blastwave_components.h"
#include "ejecta/ejecta_id.h"

namespace IMG{
    //    std::vector<std::string> m_names{"theta", "phi", "r", "theta_j", "theta0", "mu", "xrs", "yrs", "gamma", "fluxes", "intensity", "gm", "gc", "B", "tburst", "tt"};
    std::vector<std::string> m_names{
            "r", "ctheta", "cphi",
            "mu", "xrs", "yrs", "intensity",
//      "tau_compton", "tau_bh", "tau_bf"
    };
//    enum Q { itheta, iphi, ir, itheta_j, itheta0, imu, ixr, iyr, igam, iflux, iintens, igm, igc, iB, itburst, itt };
    enum Q {
        ir, ictheta, icphi,
        imu, ixr, iyr, iintens,
//      itau_comp, itau_bh, itau_bf
    };
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
    std::unique_ptr<logger> p_log = nullptr;
    /// ---
//    Vector cthetas{};
//    Vector cphis{};
    /// ---
    double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min(), ymax = std::numeric_limits<double>::min();
    Vector fluxes_shells{};
    Vector xmins_shells{};
    Vector xmaxs_shells{};
    Vector ymins_shells{};
    Vector ymaxs_shells{};
    double total_flux=0.;
    double xc = 0.;
    double yc = 0.;
    std::vector<VecVector> raw_data{}; // [ish][Q][i]
    std::vector<VecVector> cleared_raw_data{};// [ish][Q][i]
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

    /// memory optimization; remove all non-zero data from 'raw_data'
    void removeZeroEntries(){
        double frac_max_intensity = 1.e-10; // or 0
        for (size_t ish = 0; ish < raw_data.size(); ish++) {
            if (raw_data[ish][IMG::Q::iintens].size() % 2 > 0){
                (*p_log)(LOG_ERR,AT) << " expected ncells to be even (prime & counter jets). Got: ncells="
                    <<raw_data[ish][IMG::Q::iintens].size()
                    <<" for ish="<<ish
                    <<"\n";
            }
            size_t ncells = raw_data[ish][IMG::Q::iintens].size() / 2;
            auto & data = raw_data[ish];

            /// find max value in principle and counter jets
            double max_pj = std::numeric_limits<double>::min();
            double max_cj = std::numeric_limits<double>::min();
            for (size_t i = 0; i < ncells; i++)
                if (data[IMG::Q::iintens][i] > max_pj) max_pj = data[IMG::Q::iintens][i];
            for (size_t i = ncells; i < 2*ncells; i++)
                if (data[IMG::Q::iintens][i] > max_cj) max_cj = data[IMG::Q::iintens][i];
            double max_int = std::max(max_pj,max_cj);

            /// collect number of non-zero cells in both jets
            size_t n_nonzero_pj = 0;
            size_t n_nonzero_cj = 0;
            for (size_t i = 0; i < ncells; i++)
                if (data[IMG::Q::iintens][i] > max_int * frac_max_intensity)
                    n_nonzero_pj+=1;
            for (size_t i = ncells; i < 2*ncells; i++)
                if (data[IMG::Q::iintens][i] >  max_int * frac_max_intensity)
                    n_nonzero_cj+=1;

            size_t ncells_new = 0.;
            bool copy = false;

            /// if all elements are to be used; just move the VecVector()
            if ((n_nonzero_pj == ncells) || (n_nonzero_cj == ncells)) {
                ncells_new = ncells;
                cleared_raw_data.emplace_back(std::move(raw_data[ish]));
            }

            /// if only some cells are non zero; copy them into a new dynamically allocated array
            if ((n_nonzero_pj < ncells) && (n_nonzero_cj < ncells)){
                std::vector<size_t> idx_non_zeo{};
                for (size_t i = 0; i < ncells; i++){
                    if ((data[IMG::Q::iintens][i] > max_int * frac_max_intensity) ||
                        (data[IMG::Q::iintens][ncells+i] > max_int * frac_max_intensity))
                        idx_non_zeo.push_back(i);
                }
                ncells_new = idx_non_zeo.size();
                copy = true;

                for (size_t i = 0; i < ncells_new; i++)
                    idx_non_zeo.push_back(ncells+idx_non_zeo[i]);

                cleared_raw_data.emplace_back(IMG::m_names.size(), Vector(2*ncells_new, 0));

                auto & tmp = cleared_raw_data[cleared_raw_data.size()-1];

                /// copy non-zero data to new storage
                for (size_t iv = 0; iv < IMG::m_names.size(); iv++)
                    for (size_t ic = 0; ic < ncells_new; ic ++) {
                        tmp[iv][ic] = raw_data[ish][iv][idx_non_zeo[ic]];
                        tmp[iv][ic + ncells_new] = raw_data[ish][iv][idx_non_zeo[ic + ncells_new]];
                    }
            }

            (*p_log)(LOG_INFO,AT)
                << "[ish="<<ish<<"] "
                <<" Removing zeroes for Skymap at"<<" t="<<time<<" freq="<<freq<<" Fnu="<<total_flux<<" mJy"
                <<" ncells = "<<ncells<<" -> " <<ncells_new
                << (copy ? " [copied]" : " [moved]") <<"\n";

        }
        /// de-allocate memory used for original array to save space
        raw_data.clear();
        if (cleared_raw_data.size() != fluxes_shells.size()){
            (*p_log)(LOG_ERR,AT)
                << " cleared_Data.size()= " <<cleared_raw_data.size()
                << " while fluxes_shells.size()="<<fluxes_shells.size()
                <<"\n";
            exit(1);
        }
    }

    void clear(){
        raw_data.clear();
        cleared_raw_data.clear();
        binned_data.clear();
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

    if (im.cleared_raw_data.size() < 1){
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
        if ((!im.cleared_raw_data[ish].empty()) && (!im.cleared_raw_data[ish][0].empty())) {
            Output::addGroupWith1Ddata(im.cleared_raw_data[ish], group_name, IMG::m_names, file);
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

namespace IMG_TAU{
    //    std::vector<std::string> m_names{"theta", "phi", "r", "theta_j", "theta0", "mu", "xrs", "yrs", "gamma", "fluxes", "intensity", "gm", "gc", "B", "tburst", "tt"};
    std::vector<std::string> m_names{"mu", "xrs", "yrs", "tau_compton", "tau_bh", "tau_bf"};
//    enum Q { itheta, iphi, ir, itheta_j, itheta0, imu, ixr, iyr, igam, iflux, iintens, igm, igc, iB, itburst, itt };
    enum Q {imu, ixr, iyr, itau_comp, itau_bh, itau_bf};
}

static inline double cosToSin(const double &cos_theta){
    return sqrt((1.0 - cos_theta) * (1.0 + cos_theta) );
}
static inline double arccos(const double &cos_theta){
//    return 2.0 * asin( sqrt(0.5 * (1.0 - cos_theta)) );
    return 2.0 * std::asin( std::sqrt(0.5 * (1.0 - cos_theta)) );
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



#endif //SRC_IMAGE_H
