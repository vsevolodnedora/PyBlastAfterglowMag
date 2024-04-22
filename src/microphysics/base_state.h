//
// Created by vsevolod on 1/16/24.
//

#ifndef SRC_BASE_STATE_H
#define SRC_BASE_STATE_H

#include "../utilitites/utils.h"
#include "kernels.h"

struct Source{
    size_t it=0;
    double dt=-1, dlnVdt=-1, B=-1, N=-1;// r=-1, dr=-1;
    double dr=-1, vol=-1; //  comoving thickness and volume
    double N_ele_tot=-1;
    double u_e=-1;
};

#if 0
struct State{
    /// limits
    double x1=-1;
    double x2=-1;
    size_t numbins=0;
    /// arrays
    Vector e{}, de{}, half_e{}, j{}, a{}, yg{}, f{}, intensity{};
    Vector dfde{}; /// For SSA calcualtions

    Vector gam_dot_syn {}, gam_dot_syn_all {};
    Vector gam_dot_adi {}, gam_dot_adi_all {};
    Vector gam_dot_ssc {}, gam_dot_ssc_all {};

    /// for chang cooper scheme
    Vector delta_grid{}, delta_grid_bar{};

    /// output vectors
    Vector j_all{}, a_all{}, f_all{}, i_all{}, yg_all{};



    State() = default;

    void allocate(double x1_, double x2_, size_t numbins_){

        x1=x1_;
        x2=x2_;
        numbins=numbins_;

        e = TOOLS::MakeLogspaceVec(std::log10(x1),std::log10(x2), (int)numbins);

        de.resize(numbins-1);
        for (size_t i = 0; i < numbins-1; i++) // de = e[1:] - e[:-1]
            de[i] = e[i+1] - e[i];

        half_e.resize(numbins-1);
        j.resize(numbins);
        a.resize(numbins);
        yg.resize(numbins);
        f.resize(numbins);
        dfde.resize(numbins);
        intensity.resize(numbins);
    }
    /// allocate output vectors of the size numbins * num_timesteps
    void allocate_output(size_t nx, size_t nt, bool gam_dot){
        size_t nn = nx * nt; // space * time (total array size)
        f_all.resize(nn, 0.);
        j_all.resize(nn, 0.);
        a_all.resize(nn, 0.);
        i_all.resize(nn,0);
        yg_all.resize(nn,0);
        if (gam_dot){
            gam_dot_syn_all.resize(nn,0.);
            gam_dot_adi_all.resize(nn,0.);
            gam_dot_ssc_all.resize(nn,0.);
        }
    }
    /// store the current state for this timestep in total_rad vector
    void save_to_all(size_t it){
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] = f[i_];
            j_all[i_ + numbins * it] = j[i_];
            a_all[i_ + numbins * it] = a[i_];
            i_all[i_ + numbins * it] = intensity[i_];
            yg_all[i_ + numbins * it] = yg[i_];
        }
        /// store cooling terms (for electrons)
        if (gam_dot_syn.size() > 1)
            for (size_t i_ = 0; i_ < numbins; i_++) {
                gam_dot_syn_all[i_ + numbins * it] = gam_dot_syn[i_];
                gam_dot_adi_all[i_ + numbins * it] = gam_dot_adi[i_];
                gam_dot_ssc_all[i_ + numbins * it] = gam_dot_ssc[i_];
            }
    }

    void add_to_all(State & other, size_t it){
        if (numbins != other.numbins){
            std::cerr << " cannot sum states. Current nbins="<<numbins<<" other nbins="<<other.numbins<<"\n";
            exit(1);
        }
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] += other.f[i_];
            j_all[i_ + numbins * it] += other.j[i_];
            a_all[i_ + numbins * it] += other.a[i_];
            i_all[i_ + numbins * it] += other.intensity[i_];
            yg_all[i_ + numbins * it] += other.yg[i_];
        }
    }

    /**
     * Construct extra grids that are used by Chang Cooper scheme
     */
    void build_grid_chang_cooper(){

        /// cooling terms in Chang-Cooper solver (for saving later)
        gam_dot_syn.resize(numbins, 0);
        gam_dot_adi.resize(numbins, 0);
        gam_dot_ssc.resize(numbins, 0);

        if (numbins < 1){
            std::cout << " state is not initialized, empty \n";
            exit(1);
        }

        double step = std::exp((std::log(x2) / (double)numbins));
        double step_plus_one = step + 1.0;

        /// now build the grid and the half grid points
        /// we also make squared terms just incase
        for (size_t i = 0; i < numbins; i++){
            e[i] = std::pow(step, (double)i);
            if (i < numbins - 1)
                half_e[i] = .5 * e[i] * step_plus_one;
        }

        delta_grid.resize(numbins+1);
        for (size_t i = 0; i < numbins-1; i++)
            delta_grid[i+1] = e[i+1]-e[i];
        ///  we need to add extra end points to the grid
        ///  so that we can compute the delta at the boundaries
        delta_grid[0] = e[0] * (1 - 1.0 / step);
        delta_grid[numbins] = e[numbins-1] * (step - 1);

        delta_grid_bar.resize(numbins);
        for (size_t i = 0; i < numbins; i++)
            delta_grid_bar[i] = (delta_grid[i] + delta_grid[i+1])/2.;

//        std::cout<<delta_grid_bar<<"\n";
    }

};
#endif

struct BaseState{
    /// limits
    double x1=-1;
    double x2=-1;
    size_t numbins=0;
    Vector & e;
    bool is_dense= false;
    BaseState(Vector & e) : e(e) {}
    void allocateBase(double x_1, double x_2, size_t n_x, bool dense){
        x1=x_1; x2=x_2; numbins=n_x; is_dense=dense;
        if (e.size() != n_x)
            e = TOOLS::MakeLogspaceVec(std::log10(x1),std::log10(x2), (int)numbins);
    }
};

struct Electrons : BaseState{

    /// basic quantities
    Vector half_e{};
    Vector f{};
    Vector dfde{}; /// For SSA calcualtions
    /// for dense
    Vector f_all{};
    Vector gam_dot_syn {}, gam_dot_syn_all {};
    Vector gam_dot_adi {}, gam_dot_adi_all {};
    Vector gam_dot_ssc {}, gam_dot_ssc_all {};
    /// for chang cooper scheme
    Vector delta_grid{}, delta_grid_bar{};

    Electrons(Vector & e) : BaseState(e) {};

    void allocate(double x_1, double x_2, size_t n_x, size_t n_t, bool dense){

        allocateBase(x_1,x_2,n_x,dense);

        half_e.resize(numbins-1,0.);
        f.resize(numbins,0.);
        dfde.resize(numbins,0.);
        if (is_dense){
            f_all.resize(n_x*n_t, 0.);
            gam_dot_syn_all.resize(n_x*n_t,0.);
            gam_dot_adi_all.resize(n_x*n_t,0.);
            gam_dot_ssc_all.resize(n_x*n_t,0.);
        }
    }

    /// store the current state for this timestep in total_rad vector
    void save_to_all(size_t it){
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] = f[i_];
        }
        /// store cooling terms (for electrons)
        if (gam_dot_syn.size() > 1)
            for (size_t i_ = 0; i_ < numbins; i_++) {
                gam_dot_syn_all[i_ + numbins * it] = gam_dot_syn[i_];
                gam_dot_adi_all[i_ + numbins * it] = gam_dot_adi[i_];
                gam_dot_ssc_all[i_ + numbins * it] = gam_dot_ssc[i_];
            }
    }


    /**
     * Construct extra grids that are used by Chang Cooper scheme
     */
    void build_grid_chang_cooper(){

        /// cooling terms in Chang-Cooper solver (for saving later)
        gam_dot_syn.resize(numbins, 0);
        gam_dot_adi.resize(numbins, 0);
        gam_dot_ssc.resize(numbins, 0);

        if (numbins < 1){
            std::cout << " state is not initialized, empty \n";
            exit(1);
        }

        double step = std::exp((std::log(x2) / (double)numbins));
        double step_plus_one = step + 1.0;

        /// now build the grid and the half grid points
        /// we also make squared terms just incase
        for (size_t i = 0; i < numbins; i++){
            e[i] = std::pow(step, (double)i);
            if (i < numbins - 1)
                half_e[i] = .5 * e[i] * step_plus_one;
        }

        delta_grid.resize(numbins+1);
        for (size_t i = 0; i < numbins-1; i++)
            delta_grid[i+1] = e[i+1]-e[i];
        ///  we need to add extra end points to the grid
        ///  so that we can compute the delta at the boundaries
        delta_grid[0] = e[0] * (1 - 1.0 / step);
        delta_grid[numbins] = e[numbins-1] * (step - 1);

        delta_grid_bar.resize(numbins);
        for (size_t i = 0; i < numbins; i++)
            delta_grid_bar[i] = (delta_grid[i] + delta_grid[i+1])/2.;

//        std::cout<<delta_grid_bar<<"\n";
    }

};

struct Photons : BaseState{

    Vector f{};
    Vector f_all{};
    Vector j{};
    Vector j_all{};
    Vector a{};
    Vector a_all{};
    Vector intensity{};
    Vector i_all{};

    Photons(Vector & e) : BaseState(e){}

    void allocate(double x_1, double x_2, size_t n_x, size_t n_t, bool dense){

        allocateBase(x_1,x_2,n_x,dense);

        f.resize(numbins);
        j.resize(numbins);
        a.resize(numbins);
        intensity.resize(numbins);

        if (is_dense){
            f_all.resize(n_x*n_t, 0.);
            j_all.resize(n_x*n_t, 0.);
            a_all.resize(n_x*n_t, 0.);
            i_all.resize(n_x*n_t,0.);
        }
    }

    /// store the current state for this timestep in total_rad vector
    void save_to_all(size_t it){
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] = f[i_];
            j_all[i_ + numbins * it] = j[i_];
            a_all[i_ + numbins * it] = a[i_];
            i_all[i_ + numbins * it] = intensity[i_];
        }
    }

    void add_to_all(Photons & other, size_t it){
        if (numbins != other.numbins){
            std::cerr << " cannot sum states. Current nbins="<<numbins<<" other nbins="<<other.numbins<<"\n";
            exit(1);
        }
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] += other.f[i_];
            j_all[i_ + numbins * it] += other.j[i_];
            a_all[i_ + numbins * it] += other.a[i_];
            i_all[i_ + numbins * it] += other.intensity[i_];
        }
    }
};


#if 0 // todo remove if successfull
struct Electrons{
    size_t it=0;
    /// limits
    double x1=-1;
    double x2=-1;
    size_t numbins=0, nt=0;
    Vector e{}, de{}, half_e{};
    Vector dfde{}; /// For SSA calcualtions
    /// ---
    Vector o_f{}, o_gam_dot_syn {}, o_gam_dot_adi {}, o_gam_dot_ssc {};
    Vector f{}, gam_dot_syn {}, gam_dot_adi {}, gam_dot_ssc {};
    /// for chang cooper scheme
    Vector delta_grid{}, delta_grid_bar{};

    void allocate(double x1_, double x2_, size_t numbins_, size_t nt_){
        x1=x1_;
        x2=x2_;
        nt=nt_;
        numbins=numbins_;
        e = TOOLS::MakeLogspaceVec(std::log10(x1),std::log10(x2), (int)numbins);
        /// grid for derivatives
        de.resize(numbins-1);
        for (size_t i = 0; i < numbins-1; i++) // de = e[1:] - e[:-1]
            de[i] = e[i+1] - e[i];
        half_e.resize(numbins-1);
        dfde.resize(numbins,0);

        f.resize(numbins,0.);

        gam_dot_syn.resize(numbins * nt, 0);
        gam_dot_adi.resize(numbins * nt, 0);
        gam_dot_ssc.resize(numbins * nt, 0);

    }
    void allocate_output(size_t nt){
        /// resize output arrays
        o_f.resize(nt * numbins, 0.);
        o_gam_dot_syn.resize(numbins, 0);
        o_gam_dot_adi.resize(numbins, 0);
        o_gam_dot_ssc.resize(numbins, 0);

    }
    /**
     * Construct extra grids that are used by Chang Cooper scheme
     */
    void build_grid_chang_cooper(){

        /// cooling terms in Chang-Cooper solver (for saving later)
        if (numbins < 1){
            std::cout << " state is not initialized, empty \n";
            exit(1);
        }

        double step = std::exp((std::log(x2) / (double)numbins));
        double step_plus_one = step + 1.0;

        /// now build the grid and the half grid points
        /// we also make squared terms just incase
        for (size_t i = 0; i < numbins; i++){
            e[i] = std::pow(step, (double)i);
            if (i < numbins - 1)
                half_e[i] = .5 * e[i] * step_plus_one;
        }

        delta_grid.resize(numbins+1);
        for (size_t i = 0; i < numbins-1; i++)
            delta_grid[i+1] = e[i+1]-e[i];
        ///  we need to add extra end points to the grid
        ///  so that we can compute the delta at the boundaries
        delta_grid[0] = e[0] * (1 - 1.0 / step);
        delta_grid[numbins] = e[numbins-1] * (step - 1);

        delta_grid_bar.resize(numbins);
        for (size_t i = 0; i < numbins; i++)
            delta_grid_bar[i] = (delta_grid[i] + delta_grid[i+1])/2.;

//        std::cout<<delta_grid_bar<<"\n";
    }
};
struct Photons{
    size_t it=0;
    /// limits
    double x1=-1;
    double x2=-1;
    size_t numbins=0,nt=0;
    /// arrays
    Vector e{}, de{};
    /// arrays
    Vector j_syn{}, a_syn{}, f_syn{}, i_syn{};
    Vector j_ssc{}, f_ssc{}, i_ssc{};
    /// output arrays
    Vector o_j_syn{}, o_a_syn{}, o_f_syn{}, o_i_syn{};
    Vector o_j_ssc{}, o_f_ssc{}, o_i_ssc{};
    /// output arrays total
    Vector o_j{}, o_a{}, o_f{};
    void allocate(double x1_, double x2_, size_t numbins_, size_t nt_, bool is_ssc) {
        x1 = x1_;
        x2 = x2_;
        nt = nt_;
        numbins = numbins_;
        e = TOOLS::MakeLogspaceVec(std::log10(x1), std::log10(x2), (int) numbins);
        /// grid for derivatives
        de.resize(numbins - 1);
        for (size_t i = 0; i < numbins - 1; i++) // de = e[1:] - e[:-1]
            de[i] = e[i + 1] - e[i];

        /// resize arrays
        f_syn.resize(numbins, 0.);
        j_syn.resize(numbins, 0.);
        a_syn.resize(numbins, 0.);
        i_syn.resize(numbins, 0.);

        /// resize arrays
        if (is_ssc) {
            f_ssc.resize(numbins, 0.);
            j_ssc.resize(numbins, 0.);
            i_ssc.resize(numbins, 0.);
        }

        o_j.resize(nt * numbins, 0.);
        o_a.resize(nt * numbins, 0.);
    }
    void allocate_output(size_t numbins_, size_t nt_, bool is_ssc){
        nt = nt_;
        numbins = numbins_;
        /// resize output arrays
        o_f_syn.resize(nt * numbins, 0.);
        o_j_syn.resize(nt * numbins, 0.);
        o_a_syn.resize(nt * numbins, 0.);
        o_i_syn.resize(nt * numbins, 0.);

        if (is_ssc) {
            o_f_ssc.resize(nt * numbins, 0.);
            o_j_ssc.resize(nt * numbins, 0.);
            o_i_ssc.resize(nt * numbins, 0.);
        }
        /// total radiation field
        o_j.resize(nt * numbins, 0.);
        o_a.resize(nt * numbins, 0.);
        o_f.resize(nt * numbins, 0.);
    }
    /// store the current state for this timestep in total_rad vector
    void saveSynForIt(size_t it_){
        for (size_t i_ = 0; i_ < numbins; i_++){
            o_f_syn[i_ + numbins * it_] = f_syn[i_];
            o_j_syn[i_ + numbins * it_] = j_syn[i_];
            o_a_syn[i_ + numbins * it_] = a_syn[i_];
            o_i_syn[i_ + numbins * it_] = i_syn[i_];
        }
    }
    void addSynToTotal(size_t it_){
        for (size_t i_ = 0; i_ < numbins; i_++){
            o_j[i_ + numbins * it_] += j_syn[i_];
            o_a[i_ + numbins * it_] += a_syn[i_];
            o_f[i_ + numbins * it_] += f_syn[i_];
        }
    }
};
#endif

/**
 * Finite differencing stencil with boundaries using one-sided stencils
 * https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
 * @param df
 * @param x
 * @param f
 * @param acc
 */
void dfdx(Vector & df, Vector & x, Vector & f, int acc){
    double dx = 0; size_t i = 0;
//    for (i = 1; i < x.size(); i++){
//        dx = x[i] - x[i - 1]; // dx for backward or forward difference
//        df[i] = (f[i] - f[i - 1]) / dx; // backward difference
//    }
//    df[0] = df[1];

    for (i = 1; i < x.size()-1; i++){
        if (f[i] > 0) {
            dx = x[i] - x[i - 1];
            df[i] = 0.5 * (f[i + 1] - f[i - 1]) / dx;
        }
    }
    i = 0;
    if (f[i] > 0)
        df[i] = -1.5 * f[i] + 2. * f[i+1] - 0.5 * f[i+2];
    i = x.size()-1;
    if (f[i] > 0)
        df[i] = 1.5 * f[i] - 2. * f[i-1] + 0.5 * f[i-2];

//    i = 0;
//    double dx = 0;
//    size_t i = 0;
//    if (acc == 1){
//        // Assuming second-order central differencing in the comment was a mistake,
//        // and aiming for first-order accuracy in the interior points,
//        // the loop should correctly calculate differences for a first-order scheme.
//        for (i = 1; i < x.size(); i++){
//            dx = x[i] - x[i - 1]; // dx for backward or forward difference
//            df[i] = (f[i] - f[i - 1]) / dx; // backward difference
//        }
//        df[0] = df[1];
////        // one-sided stencils for boundaries
////        i = 0;
////        dx = x[i + 1] - x[i]; // Correct dx calculation for the first boundary
////        df[i] = (-3. * f[i] + 4. * f[i + 1] - f[i + 2]) / (2. * dx);
////        i = x.size() - 1;
////        dx = x[i] - x[i - 1]; // Correct dx calculation for the last boundary
////        df[i] = (3. * f[i] - 4. * f[i - 1] + f[i - 2]) / (2. * dx);
//    }
////    else if (acc == 2){
////        /// 2nd order accuracy
////        dx = x[i + 1] - x[i];
////        df[i] = (-f[i + 2] + 8. * f[i + 1] - 8. * f[i - 1] + f[i - 2]) / (12. * dx);
////        /// one-sided stencils for boundaries
////        for (i = 0; i < 2; i++)
////            df[i] = (-3. * f[i] + 4. * f[i + 1] - f[i + 2]) / (2. * (x[i + 1] - x[i]));
////        for (i = x.size()-3; i < x.size()-1; i++)
////            df[i] = (3. * f[i] - 4. * f[i - 1] + f[i - 2]) / (2. * (x[i] - x[i - 1]));
////    }
////    else
////        throw std::runtime_error("Not implemented");
////    int _x = 1;
}

#endif //SRC_BASE_STATE_H
