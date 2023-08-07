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
    std::vector<std::string> m_names{"mu", "xrs", "yrs", "intensity", "r", "ctheta", "cphi",
                                     "tau_compton", "tau_bh", "tau_bf"};
//    enum Q { itheta, iphi, ir, itheta_j, itheta0, imu, ixr, iyr, igam, iflux, iintens, igm, igc, iB, itburst, itt };
    enum Q {imu, ixr, iyr, iintens, ir, ictheta, icphi,
            itau_comp, itau_bh, itau_bf};
}

namespace IMG_TAU{
    //    std::vector<std::string> m_names{"theta", "phi", "r", "theta_j", "theta0", "mu", "xrs", "yrs", "gamma", "fluxes", "intensity", "gm", "gc", "B", "tburst", "tt"};
    std::vector<std::string> m_names{"mu", "xrs", "yrs", "tau_compton", "tau_bh", "tau_bf"};
//    enum Q { itheta, iphi, ir, itheta_j, itheta0, imu, ixr, iyr, igam, iflux, iintens, igm, igc, iB, itburst, itt };
    enum Q {imu, ixr, iyr, itau_comp, itau_bh, itau_bf};
}

struct Image {

    size_t m_n_vn = 0;
    VecVector m_data {};
    std::unique_ptr<logger> p_log;
//    size_t ia, ib, nu_ia, nu_ib;
//    double theta=0.,phi=0.,r=0.,ctheta=0.,mu=0.;
//    double flux_dens=0.;
//    double t_obs, nu_obs;
    Image(size_t size, size_t n_vn, double fill_value, unsigned loglevel=LOG_DEBUG){
        m_n_vn=n_vn;
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
//        m_data.clearEachImage();
        m_data.resize(m_n_vn);
        for (auto & arr : m_data){
            arr.resize(size, fill_value);
        }

//            m_size = size;
//            m_intens.resizeEachImage(m_size, fill_value );
//            m_xrs.resizeEachImage( m_size, fill_value );
//            m_yrs.resizeEachImage( m_size, fill_value );
//            m_rrs.resizeEachImage( m_size, fill_value );
//            m_grs.resizeEachImage( m_size, fill_value );
//            m_mu.resizeEachImage( m_size, fill_value );
//            m_theta_j.resizeEachImage( m_size, fill_value );
//            m_thetas.resizeEachImage( m_size, fill_value );
//            m_phis.resizeEachImage( m_size, fill_value );
        m_f_tot = 0.0;
    }
    ~Image(){
//        std::cerr << AT << " deleting image... size ="<<m_size<<"\n";
    }

    void resize(size_t size, double fill_value=0.){
        m_size = size;
        for (auto & arr : m_data){
            std::destroy(arr.begin(), arr.end());
        }
        m_data.resize(m_n_vn);
        for (auto & arr : m_data){
            if (arr.size() != size) {
                arr.resize(size);
                std::fill(arr.begin(), arr.end(), fill_value);
            }
            else{
                std::fill(arr.begin(), arr.end(), fill_value);
            }
        }
    }

    void clearData(){
        if ((m_size == 0)||(m_data.empty())||(m_data[0].empty())){
            (*p_log)(LOG_ERR,AT) << "cannot clean isEmpty image\n";
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
//        if(iv_n > m_n_vn-1){
//            std::cerr << AT << " Access beyong memory index="<<iv_n<<" is above name_max="<<m_n_vn-1<<"\n";
//            exit(1);
//        }
//        if (ii > m_size-1){
//            std::cerr << AT << " Access beyong memory index="<<ii<<" is above max="<<m_size-1<<"\n";
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
        if(m_data.size() != m_n_vn){
            (*p_log)(LOG_ERR, AT) << " something is wrong with the image. Exiting...";
            std::cerr << AT << "\n";
            exit(1);
        }
//        std::cout << " Reallocating [" << m_data.size() << ", " << m_data[0].size() << "]" << "\n";
        VecVector tmp;//(m_names.size(), Vector(m_size));
        tmp.resize(m_n_vn);
        for(size_t i = 0; i < m_n_vn; i++){
            tmp[i].resize(m_size);
            for (size_t j = 0; j < m_size; j++){
                tmp[i][j] = m_data[i][j];
            }
        }
        return std::move( tmp );
    }

    double m_f_tot{}; // total flux density in the image (in mJy) aka
    // J * (1.0 + z) / (2.0 * d_l * d_l) * CGS::cgs2mJy
    size_t m_size = 0;
    size_t m_size_active = 0;
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

struct Images{
    size_t m_n=0;
    std::vector<std::unique_ptr<Image>> m_images;
    Images(size_t n, size_t n_vn, size_t size=1, double fill_value=0., unsigned loglevel=LOG_DEBUG){
        if (n_vn==0 || n_vn > 100){
            std::cerr << AT << " error in images initialization\n";
            exit(1);
        }
        m_n = n;
        for (size_t i = 0; i < n; i++)
            m_images.emplace_back(std::make_unique<Image>(size,n_vn,fill_value,loglevel));
    }
    void resizeEachImage(size_t new_size){
        for (size_t i = 0; i < m_n; i++)
            m_images[i]->resize(new_size);
    }
    void clearEachImage(){
        for (size_t i = 0; i < m_n; i++)
            m_images[i]->clearData();
    }
    bool isEmpty() const {
        if (m_n == 0 || m_images.empty()){
            return true;
        }
        return false;
    }
    size_t size() const {return m_images.size();}
    std::vector<std::unique_ptr<Image>> & getImgs(){return m_images;}
    std::unique_ptr<Image> & getImg(size_t i){return m_images[i];}
    Image & getReferenceToTheImage(size_t i){
        if (i>m_images.size()-1){
            std::cerr << AT << " index is out of boundary\n";
            exit(1);
        }
        return * m_images[i];
    }
};

void combineImages(Image & image, size_t ncells, size_t nlayers, Images & images){
    if (images.size() != nlayers){
        std::cerr << " nlayeyers="<<nlayers<<" != n_images="<<images.size()<<"\n";
        std::cerr << AT << "\n";
        exit(1);
    }
    size_t size = images.getReferenceToTheImage(0).m_size;
    for (auto & im : images.getImgs()){
        if (im->m_size != size){
            std::cerr << AT << " error...\n";
        }
    }
    image.m_f_tot = 0.;
    size_t ii = 0;
    if (image.m_size!=2*ncells)
        image.resize(2 * ncells, 0. );
    for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
        size_t ncells_in_layer = EjectaID2::CellsInLayer(ilayer);//struc.cil[ilayer];
        auto & tmp = images.getReferenceToTheImage(ilayer);
//        if ( tmp.m_size != 2 * ncells_in_layer ){
//            std::cerr <<  " Error !" << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
        for( size_t ipj = 0; ipj < ncells_in_layer; ipj++ ){
            for (size_t ivn = 0; ivn < image.m_n_vn; ivn++)
                (image)(ivn, ii + ipj) = tmp(ivn, ipj);
        }
        for( size_t icj = 0; icj < ncells_in_layer; icj++ ){
            for (size_t ivn = 0; ivn < image.m_n_vn; ivn++)
                (image)(ivn, ncells + ii + icj) = tmp(ivn, ncells_in_layer + icj);
        }
        ii += ncells_in_layer;

        image.m_f_tot += tmp.m_f_tot;
    }
}

void combineImagesA(Image & image, size_t ncells, size_t nlayers, Images & images){
    if (images.size() != nlayers){
        std::cerr << " nlayeyers="<<nlayers<<" != n_images="<<images.size()<<"\n";
        std::cerr << AT << "\n";
        exit(1);
    }
    size_t size = images.getReferenceToTheImage(0).m_size;
    for (auto & im : images.getImgs()){
        if (im->m_size != size){
            std::cerr << AT << " error...\n";
        }
    }
    image.m_f_tot = 0.;
    size_t ii = 0;
    if (image.m_size!=2*ncells)
        image.resize(2 * ncells, 0. );
    for (size_t ilayer = 0; ilayer < nlayers; ilayer++){
//        size_t ncells_in_layer = EjectaID2::CellsInLayer(ilayer);//struc.cil[ilayer];
        auto & tmp = images.getReferenceToTheImage(ilayer);
        size_t ncells_in_layer = tmp.m_size_active;
        if (ncells_in_layer == 0){
            std::cerr<<AT<<" image does not have active cells\n";
            exit(1);
        }
//        if ( tmp.m_size != 2 * ncells_in_layer ){
//            std::cerr <<  " Error !" << "\n";
//            std::cerr << AT << "\n";
//            exit(1);
//        }
        for( size_t ipj = 0; ipj < ncells_in_layer; ipj++ ){
            for (size_t ivn = 0; ivn < image.m_n_vn; ivn++)
                (image)(ivn, ii + ipj) = tmp(ivn, ipj);
        }
        for( size_t icj = 0; icj < ncells_in_layer; icj++ ){
            for (size_t ivn = 0; ivn < image.m_n_vn; ivn++)
                (image)(ivn, ncells + ii + icj) = tmp(ivn, ncells_in_layer + icj);
        }
        ii += ncells_in_layer;

        image.m_f_tot += tmp.m_f_tot;
    }
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


#endif //SRC_IMAGE_H
