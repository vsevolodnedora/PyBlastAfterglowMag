//
// Created by vsevolod on 1/16/24.
//

#ifndef SRC_NUMERIC_MODEL_H
#define SRC_NUMERIC_MODEL_H

#include "base_state.h"
#include "kernels.h"

/**
 * Simple tridiagonal solver. O(n) time complexity O(n) space complexity on stack
 * @param d_j
 * @param d_j_plus_one
 * @param a
 * @param b
 * @param c
 */
static void tridiagonal_solver(
        Vector & d_j, Vector & d_j_plus_one, Vector & a, Vector & b, Vector & c_){
    size_t n_grid_points = a.size();
    double cprime[n_grid_points];
    double dprime[n_grid_points];
    /// This is the forward sweep of the tridiagonal solver
    for (size_t i = 0; i < n_grid_points; i++){
        cprime[i] = c_[i] / b[i];
        dprime[i] = d_j[i] / b[i];
    }
    /// forward_sweep if Non-zero
    double sum = 0;
    sum = std::accumulate(a.begin(), a.end(), sum);
    if (sum > 0)
        for (size_t i = 1; i < n_grid_points; i++){
            double b_minus_ac = b[i] - a[i] * cprime[i - 1];
            cprime[i] = c_[i] / b_minus_ac;
            dprime[i] = (d_j[i] - a[i] * dprime[i - 1]) / b_minus_ac;
        }

    /// This is the backwards substitution step of the tridiagonal solver.
    d_j_plus_one[n_grid_points-1] = dprime[n_grid_points-1];
    for (int j = (int)n_grid_points - 2; j > -1; j--)
        d_j_plus_one[j] = dprime[j] - cprime[j] * d_j_plus_one[j + 1];
}


/**
 * equation for the CC n_j-1 term
 * @param one_over_delta_grid  the total_rad change in energy
 * @param one_over_delta_grid_bar_backward  the forward change in energy for the second derivative
 * @param C_backward  the backward dispersion term
 * @param B_backward  the backward heating term
 * @param delta_j_minus_one  1 - delta_j
 * @return
 */
static inline double
compute_n_j_minus_one_term(
        double one_over_delta_grid, double one_over_delta_grid_bar_backward,
        double C_backward, double B_backward, double delta_j_minus_one){
    return one_over_delta_grid * (one_over_delta_grid_bar_backward * C_backward - delta_j_minus_one * B_backward);
}
/**
 * equation for the CC n_j term
 * @param one_over_delta_grid  the total_rad change in energy
 * @param one_over_delta_grid_bar_backward  the backward change in energy for the second derivative
 * @param one_over_delta_grid_bar_forward  the forward change in energy for the second derivative
 * @param C_backward  the forward dispersion term
 * @param C_forward  the backward dispersion term
 * @param B_backward  the forward heating term
 * @param B_forward  the backward heating term
 * @param one_minus_delta_j_minus_one  1 - delta_j-1
 * @param delta_j
 * @return
 */
static inline double
compute_n_j(double one_over_delta_grid, double one_over_delta_grid_bar_backward,
            double one_over_delta_grid_bar_forward,
            double C_backward, double C_forward, double B_backward, double B_forward,
            double one_minus_delta_j_minus_one, double delta_j){
    return -one_over_delta_grid * (
            (
                    one_over_delta_grid_bar_forward * C_forward
                    + one_over_delta_grid_bar_backward * C_backward
            )
            + one_minus_delta_j_minus_one * B_backward
            - delta_j * B_forward
    );
}
/**
 * equation for the CC n_j +1 term
 * @param one_over_delta_grid  the total_rad change in energy
 * @param one_over_delta_grid_bar_forward  the backward change in energy for the second derivative
 * @param C_forward  the forward dispersion term
 * @param B_forward  the forward heating term
 * @param one_minus_delta_j  1 - delta_j
 * @return
 */
static inline double
compute_n_j_plus_one(double one_over_delta_grid, double one_over_delta_grid_bar_forward,
                     double C_forward, double B_forward, double one_minus_delta_j){
    return one_over_delta_grid * (one_minus_delta_j * B_forward + one_over_delta_grid_bar_forward * C_forward);
}



class ElectronDistEvolutionBase {
protected:
    Source & source;
    Electrons & ele;
    Photons & syn;
    Photons & ssc;
    // =====================
    SynKernel & synKernel;
    SSCKernel & sscKernel;
public:
    ElectronDistEvolutionBase(Source & source, Electrons & ele, Photons & syn, Photons & ssc,
                              SynKernel & synKernel, SSCKernel & SSCKernel)
            : source(source), ele(ele), syn(syn), ssc(ssc), synKernel(synKernel), sscKernel(SSCKernel){
    }


    /**
 * Update radiation fields for current electron distribution
 */
    void update_radiation(Vector & photon_density, bool do_ssa, bool do_ssc){
        /// 1. Update synchrotron kernel for current B
        synKernel.evalSynKernel(ele.e, syn.e, source.B);

        /// 2. Compute synchrotron emission spectrum
        computeSynSpectrum();

        /// 3. compute SSA
        if (do_ssa)
            computeSSA();

        /// 4. Compute SSC spectrum (using PREVIOUS step syn & ssc spectra)
        if (do_ssc)
            computeSSCSpectrum(photon_density);

        /// 5. compute photon density
//        computePhotonDensity(photon_filed);
    }


    /**
     * Compute synchrotron emissivity assuming that electron distribution has been evolved
     * Complexity O(n*m)
     */
    void computeSynSpectrum(){
        auto & kernel = synKernel.getKernel();

        /// clear the array
        std::fill(syn.j.begin(), syn.j.end(), 0.);
//        double tmp[ele.numbins];
#ifdef USE_OPENMP
#pragma omp parallel for num_threads( 4 ) if (USE_OPENMP)
#endif
        for (size_t i = 0; i < syn.numbins - 1; i++) {
            // --- old integral
            double res = 0.;
            for (size_t j = 0; j < ele.numbins - 1; j++)
                res += kernel[i][j] * ele.f[j] * (ele.e[j+1]-ele.e[j]);
            // --- new integral
//            for (size_t j = 0; j < ele.numbins; j++)
//                tmp[j] = kernel[i][j] * ele.f[j];
//            double res_ = IntegrateSimpson38(ele.e,tmp);

//            std::cout << res << " " << res_ << " " << res / res_ << "\n";

            syn.j[i] = res;
//            syn.j[i] *= 2.3443791412546505e-22 * source.B; // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)

            if (!std::isfinite(syn.j[i]) || syn.j[i] < 0.){
                std::cerr << AT<< " nan in computeSynSpectrum()\n";
                exit(1);
            }
        }

//        /// store the number of photons per freq. bin
//        double T = source.dr / CGS::c; // escape time (See Huang+2022)
//        double fac = T / source.vol;
//        for (size_t i = 0; i < syn.numbins - 1; i++)
//            syn.f[i] = syn.j[i] * fac / (CGS::h * syn.e[i]); // (using T instead of dt) add synchrotron

    }

    /**
     * Compute synchrotron self-absorption by convolving
     * the derivative of the electron distribution
     * with emissivity kernel
     */
    void computeSSA(){
        /// compute derivative of electron distribution d (N / gam^2)
        Vector tmp (ele.numbins, 0.);
        for (size_t i = 0; i < ele.numbins; i++)
            tmp[i] = ele.f[i] / ele.e[i] / ele.e[i];
        dfdx(ele.dfde, ele.e, tmp, 1);

        /// clear the array
        std::fill(syn.a.begin(), syn.a.end(),0.);

        /// compute absorption
        auto & kernel = synKernel.getKernel();
#ifdef USE_OPENMP
#pragma omp parallel for num_threads( 4 ) if (USE_OPENMP)
#endif
        for (size_t i = 0; i < syn.numbins; i++){
            for (size_t j = 0; j < ele.numbins-1; j++)
                syn.a[i] += ele.e[j] * ele.e[j] * (ele.e[j+1]-ele.e[j]) * ele.dfde[j] * kernel[i][j];
//            syn.a[i] *= 2.3443791412546505e-22 * source.B; // np.sqrt(3) * np.power(e, 3) / h * (h/mec2)
//            syn.a[i] *= -1. / (8. * M_PI * CGS::me * std::pow(syn.e[i]/8.093440820813486e-21, 2));
            syn.a[i] *= -1. / (8. * M_PI * CGS::me * syn.e[i] * syn.e[i]); //std::pow(syn.e[i], 2));
//            if (syn.a[i] < 0){
//                std::cerr << AT<< " nan in computeSSA()\n";
//            }
            syn.a[i] = std::abs(syn.a[i]);
            if (!std::isfinite(syn.a[i])){
                std::cerr << AT << " nan in computeSSA()\n";
                syn.a[i] = 1e99;
//                exit(1);
            }
        }
//        double nu_tau = -1;
//        double tmp_[syn.numbins];
//        for (size_t i = 1; i < syn.numbins; i++) {
//            d
//            double u = 1. / 2. + exp(-dtau) / dtau - (1 - exp(-dtau)) / (dtau * dtau);
//            double intensity_approx = dtau > 1.e-3 ? em_lab * (3. * u / dtau) : em_lab;
//            syn.intensity[i] =
//            tmp_[i] = syn.intensity[i] / syn.j[i];
//        }
//        for (size_t i = syn.numbins-2; i > -1; i--){
//            if (tmp_[i] != 1){
//                nu_tau = syn.e[i];
//                break;
//            }
//        }
//        std::cout << nu_tau << ' ' <<std::accumulate(syn.a.begin(),syn.a.end(),0.) << "\n";
    }

    /**
     * Compute SSC emissivity by convolving SSC scattering
     * kernel with seed photon distribution and electron
     * distribution
     *
     * O(n*m*q) algorithm
     */
    void computeSSCSpectrum(Vector & photon_density){
        /// clear if old is present
        std::fill(ssc.j.begin(), ssc.j.end(),0.);

        /// save number of photons per freq. bin
//        double T = source.dr / CGS::c; // escape time (See Huang+2022)
//        double fac = T / source.vol;

        /// Compute SSC spectrum for each photon energy
#ifdef USE_OPENMP
#pragma omp parallel for num_threads( 4 ) if (USE_OPENMP)
#endif
        for (size_t i = 0; i < ssc.numbins; i++) {
            double scat_ph_spec_per_ele[ele.numbins - 1];

            // clean the buffer
            for (size_t i_ = 0; i_ < ele.numbins-1; i_++)
                scat_ph_spec_per_ele[i_] = 0.;

            /// Integrate seed photon spectrum [Inner integral]
            auto & kernel = sscKernel.getKernel();
            for (size_t j = 0; j < ele.numbins - 1; j++)
                for (size_t k = 0; k < syn.numbins; k++)
                    /// compute scattered photon spectrum per electron (Blumenthal 1970)
                    scat_ph_spec_per_ele[j] += 3. / 4. * CGS::sigmaT * CGS::c / (ele.e[j] * ele.e[j])
                                               * photon_density[k] / syn.e[k] * (syn.e[k+1]-syn.e[k]) * kernel[i][j][k];
//                    scat_ph_spec_per_ele[j] += syn.n[k] / syn.e[k] * syn.de[k] * kernel[i][j][k];

            /// Integrate the electron distribution [Outer integral]
            for (size_t j = 0; j < ele.numbins - 1; j++)
                ssc.j[i] += ele.f[j] * CGS::h * ssc.e[i] * scat_ph_spec_per_ele[j] * (ele.e[j+1]-ele.e[j]);

            /// add final terms
//            double ssc_freq = ssc.e[i] * CGS::ergToFreq;// / 8.093440820813486e-21;
//            ssc.j[i] *= constant * (ssc.e[i] * CGS::ergToFreq);
//            ssc.j[i] *= ssc.e[i] * CGS::h;
//            ssc.f[i] = ssc.j[i] * fac / (CGS::h * ssc.e[i]); // (using T instead of dt) add SSC
        }
    }

    /**
     * Annihilation rate
     * https://iopscience.iop.org/article/10.1088/0004-637X/732/2/77/pdf
     * https://arxiv.org/abs/2205.12146
     * @param nu_t
     * @param nu
     * @return
     */
    inline static double R_kernel(double nu, double nu_t){
        double x = nu * nu_t;
        double varbl = x * (CGS::h/CGS::me/CGS::c/CGS::c) * (CGS::h/CGS::me/CGS::c/CGS::c);
        double res = 0.;
        if (varbl - 1 > 0)
            res = 0.652 * CGS::c * CGS::sigmaT * (varbl * varbl - 1) * std::log(varbl) / (varbl*varbl*varbl);
        return res;
    }

    double computePPinjection(double nu, Vector & photon_density){
        /// find neares index to nu in syn.e grid
        size_t idx = syn.numbins-1;
        for (size_t i = 0; i < syn.numbins-1; i++)
            if ((syn.e[i] < nu) && (syn.e[i + 1] >= nu)){
                idx=i;
                break;
            }
        double integ = 0;
        for (size_t i=0; i<syn.numbins-1; i++){
            if (syn.e[i] * nu * CGS::h * CGS::h / CGS::me / CGS::me / CGS::c4 < 1.)
                continue;
            integ += photon_density[i] * R_kernel(nu, syn.e[i]) * (syn.e[i + 1] - syn.e[i]); //
        }
        double fac = 4. * CGS::me * CGS::c2 / CGS::h;
        double res = fac * photon_density[idx] * integ; // Eq.51 in Micelli&Nava
//
//        double integ = 0;
//        int k1 = -1;
//        if (nu >= syn.e[syn.numbins-1])
//            k1 = syn.numbins-1;
//        for (size_t i=0; i<syn.numbins-1; i++){
//            if (syn.e[i]*nu*CGS::h*CGS::h/CGS::me/CGS::me/CGS::c4 < 1.)
//                continue;
//            if (nu > ele.e[i])
//                k1 = i;
//            integ += photon_density[i] * R_kernel(nu * syn.e[i]) * (syn.e[i+1]-syn.e[i]);
//        }
//        double res = CGS::c * CGS::sigmaT * photon_density[idx] * integ;
//        if (!std::isfinite(res) || res < 0){
//            std::cerr << AT << " res = "<<res<< " in pp-roduction\n";
//            exit(1);
//        }
////        if (res!=0){
////            int z = 0;
////        }
        return res;
    }

    void computePPabsorption(Vector & photon_density){
        /// clear the array
        std::fill(ssc.a.begin(), ssc.a.end(),0.);
#ifdef USE_OPENMP
#pragma omp parallel for num_threads( 4 ) if (USE_OPENMP)
#endif
        for (size_t i = 0; i < ssc.numbins; i++) {
            double integ = 0;
            for (size_t j = 0; j < syn.numbins-1; j++) {
                /// Eq. 49 in Micelli & Nava
                integ += R_kernel(ssc.e[i], syn.e[j]) * photon_density[j] * (syn.e[j+1]-syn.e[j]);
            }
            ssc.a[i] = integ / CGS::c;
        }
//        std::cout << ssc.e << "\n";
//        std::cout << ssc.a << "\n";
//        print_xy_as_numpy(ssc.e,ssc.a,ssc.numbins,1);
//        int z =1;
    }

    /**
     * Define a factor that is dependent on the magnetic field B.
     * DT = 6 * pi * me*c / (sigma_T B^2 gamma)
     * DT = (1.29234E-9 B^2 gamma)^-1 seconds
     * DT = (cool * gamma)^-1 seconds
     * where B is in Gauss.
     * @param B
     * @return
     */
    double inline gamma_dot_syn(double gam) const{
        double const_factor = CGS::sigmaT / (6.*M_PI*CGS::me*CGS::c);
        return const_factor * source.B * source.B * gam * gam;
    }

    /**
     * Adiabatic cooling for electron
     * @param gam
     * @return
     */
    double inline gamma_dot_adi(double gam) const{
        return source.dlnVdt * (gam * gam - 1.) / (3. * gam);
    }

    /**
     * Integrate SSC kernel for electrons at index 'idx' O(n)
     * @param i
     * @param gam
     * @return
     */
    double inline gamma_dot_ssc(size_t idx, double gam, Vector & photon_density) const{
        auto & kernel = sscKernel.getKernelInteg(); // [i_gamma, i_energy_syn]
        double res = 0;
        /// integrate the seed photon spectrum (using pre-computed inner integral)
        for (size_t j = 0; j < syn.numbins-1; j++)
            res += photon_density[j] / syn.e[j] * (syn.e[j+1]-syn.e[j]) * kernel[idx][j];
        /// add constant from the paper
        double constant = 3./4. * CGS::h * CGS::sigmaT / CGS::me / CGS::c;
        return constant * res / gam / gam;
    }


    /**
     * Integrate SSC kernel for electrons at index 'idx' O(n)
     * @param idx
     * @return
     */
//    double sscIntegralTerm(size_t idx){
//
//        /// integrate photon density distribution over electron distribution
//        double res = 0;
//        auto & kernel = sscKernel.getKernelInteg(); // [i_gamma, i_energy_syn]
//
//        /// integrate the seed photon spectrum (using pre-computed inner integral)
//        for (size_t j = 0; j < syn.numbins-1; j++)
//            res += (syn.f[j]+ssc.f[j]) / syn.e[j] * syn.de[j] * kernel[idx][j];
//
//        /// full integral (very slow...)
////        auto & kernel_ = sscKernel.getKernel(); // [i_nu_ssc, i_gam, i_nu_syn]
////        double integ = 0;
////        for (size_t k = 0; k < ssc.numbins-1; k++) {
////            double inner = 0;
////            for (size_t j = 0; j < syn.numbins-1; j++){
////                inner+=syn.de[j]/syn.e[j]*photon_filed[j]*kernel_[k][idx][j];
////            }
////            integ+=inner*ssc.e[k]*ssc.de[k];
////        }
////        if (std::abs(res/integ) > 10 || std::abs(res/integ) < 0.1){
////            std::cerr << " res="<<res<<" integ="<<integ<<" ratio="<<res/integ<<"\n";
////        }
////        res= integ;
////        std::cout << std::accumulate
////        double tot = std::accumulate(syn.f.begin(),syn.f.end(),0.);
////        if (tot>0) {
////            std::cout << std::accumulate(syn.f.begin(), syn.f.end(), 0.) << "\n";
////            int z = 1;
////        }
//
//
//        /// add constant from the paper
//        double constant = 3./4. * CGS::h * CGS::sigmaT / CGS::me / CGS::c;
//        res *= constant;
//
//        return res;
//    }
};

class ChangCooper : public ElectronDistEvolutionBase{
    size_t n_grid_points = 0;
    Vector dispersion_term {};
    Vector heating_term {};
    Vector source_grid {};
    Vector escape_grid {};
    size_t iterations = 0;
    double current_time = 0;
    Vector delta_j {};
    Vector a{},b{},c{},d{}; // for solver (d is the solution)
    // =====================
    VecVector ssc_inner_integral {};
    std::vector<VecVector> ssc_kernel {};
    bool is_ssc = false;
    // -------------------
public:

    ChangCooper(Source &source, Electrons &ele, Photons &syn, Photons &ssc,
                SynKernel & synKernel, SSCKernel & SSCKernel)
            : ElectronDistEvolutionBase(source, ele, syn, ssc, synKernel, SSCKernel) {
    }

    void setSolver(){
        /// allocate memory for additional arrays
        n_grid_points = ele.numbins;
        dispersion_term.resize(n_grid_points, 0);
        heating_term.resize(n_grid_points, 0);
        source_grid.resize(n_grid_points, 0);
        escape_grid.resize(n_grid_points, 0);
        dispersion_term.resize(n_grid_points, 0);
        delta_j.resize(n_grid_points, 0);
        a.resize(n_grid_points, 0);
        b.resize(n_grid_points, 0);
        c.resize(n_grid_points, 0);
        d.resize(n_grid_points, 0);
//        gam_dot_syn.resize(n_grid_points, 0);
//        gam_dot_adi.resize(n_grid_points, 0);
//        gam_dot_ssc.resize(n_grid_points, 0);
        /// reset grids/terms for new calculation
        resetSolver();
        if (not sscKernel.getKernel().empty())
            is_ssc = true;
    }

    void resetSolver(){
        n_grid_points = ele.numbins;
        std::fill(dispersion_term.begin(), dispersion_term.end(),0.);
        std::fill(heating_term.begin(), heating_term.end(),0.);
        std::fill(source_grid.begin(), source_grid.end(),0.);
        std::fill(escape_grid.begin(), escape_grid.end(),0.);
        std::fill(dispersion_term.begin(), dispersion_term.end(),0.);
        std::fill(delta_j.begin(), delta_j.end(),0.);
        std::fill(a.begin(), a.end(),0.);
        std::fill(b.begin(), b.end(),0.);
        std::fill(c.begin(), c.end(),0.);
        std::fill(d.begin(), d.end(),0.);
//        std::fill(gam_dot_syn.begin(), gam_dot_syn.end(),0.);
//        std::fill(gam_dot_adi.begin(), gam_dot_adi.end(),0.);
//        std::fill(gam_dot_ssc.begin(), gam_dot_ssc.end(),0.);
        iterations = 0;
        current_time = 0.;
    }

    /**
     * delta_j controls where the differences are computed. If there are no dispersion
     * terms, then delta_j is zero
     */
    void computeDeltaJ(){
        for (size_t j = 0; j < n_grid_points; j++){
            delta_j[j] = 0.;
            // if the dispersion term is 0 => delta_j = 0
            if (dispersion_term[j] != 0){
                // w = ( self.ele._delta_grid[1:-1][j] * self._heating_term[j] ) / self._dispersion_term[j] # REPLACED
                double w = ( ele.delta_grid[j+1] * heating_term[j] ) / dispersion_term[j];
                // w asymptotically approaches 1/2, but we need to set it manually
                if (w == 0)
                    delta_j[j] = 0.5;
                    // otherwise, we use appropriate bounds
                else
                    delta_j[j] = (1.0 / w) - 1.0 / (std::exp(w) - 1.0);
            }
            // one_minus_delta_j[j] = 1 - delta_j[j]
        }
//        int x = 1;
    }

    /**
     * from the specified terms in the subclasses, setup the tridiagonal terms
     * @param delta_t
     */
    void setupVectors(double delta_t){

        double one_over_delta_grid_forward,one_over_delta_grid_backward,
                one_over_delta_grid_bar,B_forward,B_backward,C_forward,C_backward;

        /// clean to prevent leackage
        for (size_t k = 0; k < n_grid_points; k++){
            a[k] = 0.; b[k] = 0.; c[k] = 0.;
        }

        /// walk backwards in j starting from the second to last index; then set the end points
        for (size_t k = n_grid_points - 2; k > 0; k--){
            /// pre compute one over the delta of the grid
            /// this is the 1/delta_grid in front of the F_j +/- 1/2.
            one_over_delta_grid_forward = 1.0 / ele.delta_grid[k + 1];
            one_over_delta_grid_backward = 1.0 / ele.delta_grid[k];
            /// this is the delta grid in front of the full equation
            one_over_delta_grid_bar = 1.0 / ele.delta_grid_bar[k];
            /// The B_j +/- 1/2 from CC
            B_forward = heating_term[k];
            B_backward = heating_term[k - 1];
            /// The C_j +/- 1/2 from CC
            C_forward = dispersion_term[k];
            C_backward = dispersion_term[k - 1];
            /// in order to keep math errors at a minimum, the tridiagonal terms
            /// are computed in separate functions so that boundary conditions are set consistently.
            /// First we solve (N - N) = F
            /// then we will move the terms to form a tridiagonal equation
            /// n_j-1 term
            a[k] = compute_n_j_minus_one_term(
                    one_over_delta_grid_bar,
                    one_over_delta_grid_backward,
                    C_backward,
                    B_backward,
                    delta_j[k - 1]
            );
            /// n_j term
            b[k] = compute_n_j(
                    one_over_delta_grid_bar,
                    one_over_delta_grid_backward,
                    one_over_delta_grid_forward,
                    C_backward,
                    C_forward,
                    B_backward,
                    B_forward,
                    1. - delta_j[k - 1],
                    delta_j[k]
            );
            /// n_j+1 term
            c[k] = compute_n_j_plus_one(
                    one_over_delta_grid_bar,
                    one_over_delta_grid_forward,
                    C_forward,
                    B_forward,
                    1. - delta_j[k]
            );
        }

        /// now set the end points
        /// ----------- right boundary
        /// j+1/2 = 0
        one_over_delta_grid_forward = 0.0;
        one_over_delta_grid_backward = 1.0 / ele.delta_grid[ele.delta_grid.size()-1];

        one_over_delta_grid_bar = 1.0 / ele.delta_grid_bar[ele.delta_grid_bar.size()-1];

        /// n_j-1 term
        a[a.size()-1] = compute_n_j_minus_one_term(
                one_over_delta_grid_bar,
                one_over_delta_grid_backward,
                dispersion_term[dispersion_term.size()-1],
                heating_term[heating_term.size()-1],
                delta_j[delta_j.size()-1]
        );

        /// n_j term
        b[b.size()-1] = compute_n_j(
                one_over_delta_grid_bar,
                one_over_delta_grid_backward,
                one_over_delta_grid_forward,
                dispersion_term[dispersion_term.size()-1],
                0.,
                heating_term[heating_term.size()-1],
                0.,
                1.-delta_j[delta_j.size()-1],
                0.
        );

        /// n_j+1 term
        c[c.size()-1] = 0.;

        /// ----------- left boundary
        /// j-1/2 = 0
//        double one_over_delta_grid = 1.0 / (ele.half_e[0] - ele.e[0] / std::sqrt(ele.step));
        double one_over_delta_grid_bar_forward = 1.0 / ele.delta_grid_bar[0];
        double one_over_delta_grid_bar_backward = 0.0;

        one_over_delta_grid_forward = 1.0 / ele.delta_grid[0];
        one_over_delta_grid_backward = 0;

        one_over_delta_grid_bar = 1.0 / ele.delta_grid_bar[0];

        /// n_j-1 term
        a[0] = 0.0;

        /// n_j term
        b[0] = compute_n_j(
                one_over_delta_grid_bar,
                one_over_delta_grid_backward,
                one_over_delta_grid_forward,
                0,
                dispersion_term[0],
                0,
                heating_term[0],
                0,
                delta_j[0]
        );

        /// n_j+1 term
        c[0] = compute_n_j_plus_one(
                one_over_delta_grid_bar,
                one_over_delta_grid_forward,
                dispersion_term[0],
                heating_term[0],
                1.-delta_j[0]
        );

        /// carry terms to the other side to form a tridiagonal equation
        /// the escape term is added on but is zero unless created in other func.
        for (size_t i = 0; i < n_grid_points; i++){
            a[i] *= -delta_t;
            b[i] = (1. - b[i] * delta_t) + escape_grid[i] * delta_t;
            c[i] *= -delta_t;
        }


    }

    /**
     * Evolve over substep
     * @param delta_t
     * @return
     */
    void solve(double delta_t){
        /// set up the right side of the tridiagonal equation.
        /// This is the current distribution plus the source
        /// unless it is zero
        for (size_t i = 0; i < n_grid_points; i++) {
            d[i] = ele.f[i] + source_grid[i] * delta_t;
            if (!std::isfinite(ele.f[i]) || ele.f[i] < 0.){
                std::cerr << AT<< " nan in computeSynSpectrum()\n";
                exit(1);
            }
        }

        /// now make a tridiagonal_solver for these terms
//        TridiagonalSolver tridiagonalSolver = TridiagonalSolver(a,b,c);

        /// insert new solution to electron distribution object
//        tridiagonalSolver.solve(d, ele.f);

        tridiagonal_solver(d, ele.f, a, b, c);



        /// increase the run iterator and the current time
        iterations += 1;
        /// increase the time
        current_time += delta_t;
    }

    /**
     * Power law injection of electrons
     * @param x_min
     * @param x_max
     * @param index
     * @param N
     */
    void setSourceFunction(double x_min, double x_max, double index){
        double norm = (1 + index) / (std::pow(x_max, index + 1) - std::pow(x_min, index + 1));
        for (size_t i = 0; i < n_grid_points; i++)
            if ((ele.e[i] > x_min) && (ele.e[i] <= x_max))
                source_grid[i] = source.N * norm * std::pow(ele.e[i], index);
            else
                source_grid[i] = 0.;
    }


    void addPPSource(Vector & photon_density){
        for (size_t i = 0; i < n_grid_points; i++){
            double nu1 = 2. * ele.e[i] * CGS::me * CGS::c * CGS::c / CGS::h;
            double n_ele_from_pp = computePPinjection(nu1, photon_density) * source.vol;
            source_grid[i] += n_ele_from_pp;
        }
    }

    /**
     * Escape function of elections
     * @param x_min
     * @param x_max
     */
    void setEscapeFunction(double x_min, double x_max){
        for (size_t i = 0; i < n_grid_points; i++) {
            escape_grid[i] = 0.;
            if (ele.e[i] > x_min && ele.e[i] <= x_max)
                escape_grid[i] = 0;
        }

    }

    /**
     * Compute gamma_dot term in kinteic equation
     * @param cool_index
     */
    void setHeatCoolTerms(double cool_index, Vector & photon_density){
        /// Assume no dispersion so
        std::fill(dispersion_term.begin(), dispersion_term.end(),0.);

        /// set heating / cooling
#ifdef USE_OPENMP
#pragma omp parallel for num_threads( 4 ) if (USE_OPENMP)
#endif
        for (size_t i = 0; i < n_grid_points-1; i++){
            /// compute synchrotron cooling
            ele.gam_dot_syn[i] = gamma_dot_syn(ele.half_e[i]);
//                    synchCoolComponent(source.B) * ele.half_e[i] * ele.half_e[i];
//                            * std::pow(ele.half_e[i], cool_index);

            /// compute adiabatic cooling
            ele.gam_dot_adi[i] = gamma_dot_adi(ele.half_e[i]);
//                    source.dlnVdt \
//                            * (ele.half_e[i] * ele.half_e[i] - 1.) / (3. * ele.half_e[i]);

            /// compute SSC term
            if (is_ssc)
                ele.gam_dot_ssc[i] = gamma_dot_ssc(i, ele.half_e[i], photon_density);
//                        sscIntegralTerm(i) / ele.half_e[i] / ele.half_e[i];

            /// overall heating/cooling term
            heating_term[i] = ele.gam_dot_syn[i] + ele.gam_dot_adi[i] + ele.gam_dot_ssc[i];

//            if (ssc_term != 0 and ssc_term > syn_term){
//                int z = 1;
//            }
//            if (adi_term < syn_term){
//                int z = 1;
//            }
//            if (ele.gam_dot_ssc[i]>0) {
//                std::cout << ele.gam_dot_ssc[i] << "\n";
//                int z = 0;
//            }

            if (!std::isfinite(heating_term[i])){
                std::cerr << AT<< " nan in setHeatCoolTerms()\n";
                exit(1);
            }
        }
//        if (std::accumulate(ele.gam_dot_ssc.begin(),ele.gam_dot_ssc.end(),0.)>0.){
//            std::cout<<std::accumulate(syn.f.begin(),syn.f.end(),0.)<<"\n";
//            std::cout<<ele.gam_dot_ssc<<"\n";
//            int z = 1;
//        }
    }


};

#endif //SRC_NUMERIC_MODEL_H
