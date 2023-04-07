//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_IMAGE_H
#define SRC_IMAGE_H

#include "utilitites/pch.h"
#include "utilitites/utils.h"
#include "utilitites/interpolators.h"
#include "blastwave_components.h"
#include "initial_data.h"

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
    size_t size = images[0].m_size;
    for (auto & im : images){
        if (im.m_size != size){
            std::cerr << AT << " error...\n";
        }
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


#endif //SRC_IMAGE_H
