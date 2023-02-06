//
// Created by vsevolod on 05/02/23.
//

#ifndef SRC_SSC_H
#define SRC_SSC_H

#include "utils.h"

/// Initilizing source specific variable dictionary
struct Source{
    double h0 = 0.0;
    double Omega_M = 0.0;
    double Omega_Lambda = 0.0;
    double distance = 0.0;
    double redshift = 0.0;
    double radius = 0.0;
    double volume = 0.0;
    double B = 0.0;
    double Gamma = 0.0;
    double theta = 0.0;
    double beta = 0.0;
    double delta = 0.0;
    double lum = 0.0;

    double ECflag = 0.0;
    double Mbh = 0.0;
    double Mdot = 0.0;
    double blobDist = 0.0;

    double E_min = 0.0;
    double E_max = 0.0;
    double E_break = 0.0;
    double pow1 = 0.0;
    double pow2 = 0.0;
    double gamma_max = 0.0;

    double Nu_max_sy = 0.0;
    double Nu_max_ic = 0.0;
    double Po_max_sy = 0.0;
    double Po_max_ic = 0.0;
    double od1 = 0.0;
    double pp1 = 0.0;

    double Nu_KN_lim1 = 0.0;
    double Nu_KN_lim2 = 0.0;
    double chi2 = 0.0;
    double dof = 0;

    std::string EBLModel = "";
};

/// Initialises beginning spectrum
struct InitSt{
    double x1;//=p1,
    double x2;//' : p2,
    int numbins;//' : numbins,
    double delta;//' : (p2 - p1)/numbins,
    Vector f;//' : [0]*numbins,
    Vector e;//' : [0]*numbins,
    Vector de;//' : [0]*numbins,
    Vector nu;//' : [0]*numbins,
    Vector j;//' : [0]*numbins,
    Vector n;//' : [0]*numbins,
    Vector a;//' : [0]*numbins,
    Vector tau_i;//' : [0]*numbins,
    VecVector tau_e;//' : np.zeros(shape = (NMODEL,MAX_BINS))
    InitSt(double p1, double p2, int _numbins, int NMODEL=5, int MAX_BINS=1001){
        x1=p1;
        x2=p2;
        numbins = _numbins;
        delta = (p2-p1)/numbins;
        f.resize(numbins,0.);
        e.resize(numbins,0.);
        de.resize(numbins,0.);
        nu.resize(numbins,0.);
        j.resize(numbins,0.);
        n.resize(numbins,0.);
        a.resize(numbins,0.);
//        tau_e.resize(numbins,0.);
        tau_e.resize(NMODEL); //
        for (auto & arr : tau_e)
            arr.resize(MAX_BINS);
        tau_i.resize(numbins,0.);
        ///
        for (int i=0; i<numbins; i++){
            e[i] = xen(x1, delta, i);
            de[i] = fxen(x1, delta, i+0.5) - fxen(x1,delta,i-0.5);
            nu[i] = xen(x1, delta, i);
        }
    }
    // Returns the center value of a given float bin index
    static inline double xen(double x1, double delta, int i){
        double val = x1 + (double(i) + 0.5) * delta;
        return std::pow(10.,val);
    }
    static inline double fxen(double x1, double delta, double i){
        double val = double(x1 + (i + 0.5) * delta);
        return std::pow(10.,val);
    }
};

double NedWright(double z, double h0, double Omega_M, double Omega_Lambda, double age_Gyr, double zage_Gyr,
               double DTT_Gyr, double DA_Mpc, double kpc_DA, double DL_Mpc, double DL_Gyr) {

    double n = 1000;    //number of integral points
    double c = 299792.458; //c in km / sec
    double Tyr = 977.8; // 1 / Hto Gyr

    double H0 = h0 * 100;    //Hubble constant

    /// 'Densities'
    double WM = Omega_M;  // Omega(matter)
    double WV = Omega_Lambda; // Omega(lambda)
    double h = H0 / 100;
    double WR = 4.165E-5 / (h * h);   //  Omega(radiation), includes 3 massless neutrino species, T0 = 2.72528
    double WK = 1 - WM - WR - WV;  //  Omega curvaturve = 1 - Omega(total)

    double a = 1.0; //    scale factor
    double az = 1.0 / (1 + 1.0 * z);

    double  age = 0; //
    for (size_t i=0; i<n; i++){
        a = az * (i + 0.5) / n;
        double adot = std::sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a));
        age = age + 1 / adot;
    }
    double zage = az * age / n;

    /// Correction for annihilations of particles not present not like e +/e -

    double lpz = std::log((1 + 1.0 * z)) / std::log(10.0);
    double dzage = 0;
    if (lpz > 7.500)
        dzage = 0.002 * (lpz - 7.500);
    if (lpz > 8.000)
        dzage = 0.014 * (lpz - 8.000) + 0.001;
    if (lpz > 8.500)
        dzage = 0.040 * (lpz - 8.500) + 0.008;
    if (lpz > 9.000)
        dzage = 0.020 * (lpz - 9.000) + 0.028;
    if (lpz > 9.500)
        dzage = 0.019 * (lpz - 9.500) + 0.039;
    if (lpz > 10.000)
        dzage = 0.048;
    if (lpz > 10.775)
        dzage = 0.035 * (lpz - 10.775) + 0.048;
    if (lpz > 11.851)
        dzage = 0.069 * (lpz - 11.851) + 0.086;
    if (lpz > 12.258)
        dzage = 0.461 * (lpz - 12.258) + 0.114;
    if (lpz > 12.382)
        dzage = 0.024 * (lpz - 12.382) + 0.171;
    if (lpz > 13.055)
        dzage = 0.013 * (lpz - 13.055) + 0.188;
    if (lpz > 14.081)
        dzage = 0.013 * (lpz - 14.081) + 0.201;
    if (lpz > 15.107)
        dzage = 0.214;
    zage = zage * std::pow(10.0, dzage);
    zage_Gyr = (Tyr / H0) * zage;

    double DTT = 0.0;
    double DCMR = 0.0;

    /// integral over a=1/(1+z) from az to 1 in n steps using midpoint rule
    for (size_t i = 0; i < n; i++) {
        a = az + (1. - az) * (i + 0.5) / n;
        double adot = std::sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a));
        DTT = DTT + 1. / adot;
        DCMR = DCMR + 1. / (a * adot);
    }

    DTT = (1 - az) * DTT / n;
    DCMR = (1 - az) * DCMR / n;

    age = DTT + zage;
    age_Gyr = age * (Tyr / H0);


    DTT_Gyr = (Tyr / H0) * DTT;


    // 'tangential comoving distance'
    double ratio = 1.00;

    double x = std::sqrt(std::abs(WK)) * DCMR;
    double y;
    if (x > 0.1) {
        if (WK > 0)
            ratio = 0.5 * (std::exp(x) - std::exp(-x)) / x;
        else
            ratio = std::sin(x) / x;
        y = ratio * DCMR;
    }
    else {
        y = x * x;
        if (WK < 0)
            y = -y;

        ratio = 1 + y / 6 + y * y / 120;
        y = ratio * DCMR;
    }
    double DCMT = y;

    double DA = az * DCMT;
    DA_Mpc = (c / H0) * DA;
    kpc_DA = DA_Mpc / 206.264806;
    double DL = DA / (az * az);
    DL_Mpc = (c / H0) * DL;

    printf("Distance parameters: age_Gyr %e zage_Gyr %e DTT_Gyr %e DA_Mpc %e kpc_DA %e DL_Mpc %e \n",
          age_Gyr, zage_Gyr, DTT_Gyr, DA_Mpc, kpc_DA, DL_Mpc);

    return DL_Mpc;
}

class SSC{
    std::unique_ptr<Source> p_source;
    int NUMEL =400;
    int NUMSYNC = 400;
    int NUMIC = 73;
    int NUMSUM = 100;
public:

    SSC(){
        p_source = std::make_unique<Source>();
    }
    void run(double z, double gamma, double theta, double B, double r, double w_p_soll,
             double Emin, double Emax, double Ebreak, double p1, double p2,
             double ECFlag, double Mbh, double Mdot, double blobDist, double EBLModel){
        p_source->h0 = 0.7;
        p_source->Omega_M = 0.3;
        p_source->Omega_Lambda = 0.7;
        double age_Gyr=0, zage_Gyr=0, DTT_Gyr=0, DA_Mpc=0, kpc_DA=0, DL_Mpc=0, DL_Gyr=0;
        /// Calling NedWright to calculate DL_Mpc (luminosity distance) for further calculations
        DL_Mpc = NedWright(z, p_source->h0, p_source->Omega_M, p_source->Omega_Lambda,
                           age_Gyr, zage_Gyr, DTT_Gyr, DA_Mpc, kpc_DA, DL_Mpc, DL_Gyr);
        /// Mpc to m
        double d_l = DL_Mpc * 1e6*CGS::SI_pc;
        printf("Luminosity distance %e Mpc \n", DL_Mpc);
        p_source->redshift = float(z);
        p_source->distance = d_l;
        //    Initilising variables for electrons, synchrotron, inverse Compton, disk, external Compton and sum spectra.
        //    arg 1 = Dictionary
        //    arg 2 = log10 of lower energy limit
        //    arg 3 = log10 of upper energy limit
        //    arg 4 = Number of bins. Greater value = more accuracy
        //-----
        //    Electron Spectrum from 1e4 to 1e15.
        //    Energy in eV = ple['e']
        //    Electron Density in m^-3 eV^-1 = ple['f']
        // ----
        InitSt ple(4, 15, NUMEL);
        //    Synchrotron Spectrum from 1e4 to 1e25
        //    Frequency in Hz = psync['e']
        //    Flux in  erg cm^-2 s^-2 = psync['f']
        InitSt psync (4.0, 25.0, NUMSYNC);
        //    Inverse Compton Spectrum from 1e10 to 1e27.5
        //    Frequency in Hz = pic['e']
        //    Flux in  erg cm^-2 s^-2 = pic['f']
        InitSt pic (10.0, 27.5, NUMIC);
        //    Disk Spectrum from 1e4 to 1e25
        //    Frequency in Hz = pbb['e']
        //    Flux in  erg cm^-2 s^-2 = pbb['f']
        InitSt pbb (4.0, 25.0, NUMSYNC);
        //    External Compton Spectrum from 1e10 to 1e28.5
        //    Frequency in Hz = pec['e']
        //    Flux in  erg cm^-2 s^-2 = pec['f']
        InitSt pec (10.0, 28.5, NUMIC);
        //    Sum Spectrum from 1e4 to 1e30
        //    Frequency in Hz = psum['e']
        //    Flux in  erg cm^-2 s^-2 = psum['f']
        InitSt psum(4.0, 30.0, NUMSUM);
        ///
        int idn = 1; // ???????
        p_source->lum = 1; // ?????
        p_source->Gamma = gamma;
        p_source->theta = theta/CGS::RtD;
        p_source->B = B;
        p_source->radius = r;
        w_p_soll = w_p_soll;
        p_source->beta = std::sqrt(1-1/(p_source->Gamma*p_source->Gamma));
        p_source->delta = 1./(p_source->Gamma*(1.-p_source->beta*std::cos(p_source->theta)));

        p_source->volume = 4./3.*CGS::pi*(p_source->radius*p_source->radius*p_source->radius);

        p_source->E_min =  Emin;
        p_source->E_max =  Emax;
        p_source->E_break =  Ebreak;

        if(p_source->E_min < ple.x1)
            p_source->E_min'] =  ple['x1'] + ple['delta']

        if(p_source->E_max > ple['x2']):
            p_source->E_max = ple['x2'] - ple['delta']

        if(p_source->E_break < source['E_min']):
            p_source->E_break = source['E_min']

        if(p_source->E_break'] > source['E_max'] ):
            p_sourcevE_break']  = source['E_max']

        p_source->chi2'] = 0.0

        p_source->pow1'] = -p1
        p_source->pow2'] = -p2

        p_source->ECflag'] = ECFlag
        p_source->Mbh'] = Mbh
        p_source->Mdot'] = Mdot
        p_source->blobDist'] = blobDist
        p_source->EBLModel'] = EBLModel

    }
};

#endif //SRC_SSC_H
