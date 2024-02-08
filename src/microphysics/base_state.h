//
// Created by vsevolod on 1/16/24.
//

#ifndef SRC_BASE_STATE_H
#define SRC_BASE_STATE_H

#include "../utilitites/utils.h"

struct Source{
    double dt=-1, dlnVdt=-1, B=-1, r=-1, dr=-1, N=-1;
};

struct State{
    /// limits
    double x1=-1;
    double x2=-1;
    size_t numbins=0;
    /// arrays
    Vector e{}, de{}, half_e{}, j{}, a{}, n{}, f{};
    Vector dfde{}; /// For SSA calcualtions

    /// for chang cooper scheme
    Vector delta_grid{}, delta_grid_bar{};

    /// output vectors
    Vector j_all{}, a_all{}, f_all{};

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
        n.resize(numbins);
        f.resize(numbins);
        dfde.resize(numbins);
    }
    /// allocate output vectors of the size numbins * num_timesteps
    void allocate_output(size_t nx, size_t nt){
        size_t nn = nx * nt;
//        Vector tmp = TOOLS::MakeLogspaceVec(std::log10(x1_),std::log10(x2_), (int)nx);
//        e_all.resize(nn, 0.);
//        size_t ii = 0;
//        for (size_t i_ = 0; i_ < nt; i_++) // time steps or radial steps
//            for (size_t j_ = 0; j_ < nx; j_++) {// freq or gam steps
//                e_all[ii] = tmp[j_];
//                ii++;
//            }
        f_all.resize(nn, 0.);
        j_all.resize(nn, 0.);
        a_all.resize(nn, 0.);
    }
    /// store the current state for this timestep in total_rad vector
    void save_to_all(size_t it){
        for (size_t i_ = 0; i_ < numbins; i_++){
            f_all[i_ + numbins * it] = f[i_];
            j_all[i_ + numbins * it] = j[i_];
            a_all[i_ + numbins * it] = a[i_];
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
        }
    }

    /**
     * Construct extra grids that are used by Chang Cooper scheme
     */
    void build_grid_chang_cooper(){

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

        std::cout<<delta_grid_bar<<"\n";
    }

};

/**
 * Finite differencing stencil with boundaries using one-sided stencils
 * https://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf
 * @param df
 * @param x
 * @param f
 * @param acc
 */
void dfdx(Vector & df, Vector & x, Vector & f, int acc){
    double dx=0;
    size_t i = 0;
    if (acc == 1){
        /// 1st order accuracy
        for (i = 1; i < x.size()-1; i++){
            dx = x[i + 1] - x[i];
            df[i] = (f[i + 1] - f[i - 1]) / (2. * dx);
        }
        /// one-sided stencils for boundaries
        i = 0;
        df[i] = (-3. * f[i] + 4. * f[i + 1] - f[i + 2]) / (2. * (x[i + 1] - x[i]));
        i = x.size() - 1;
        df[i] = (3. * f[i] - 4. * f[i - 1] + f[i - 2]) / (2. * (x[i] - x[i - 1]));
    }
    else if (acc == 2){
        /// 2nd order accuracy
        dx = x[i + 1] - x[i];
        df[i] = (-f[i + 2] + 8. * f[i + 1] - 8. * f[i - 1] + f[i - 2]) / (12. * dx);
        /// one-sided stencils for boundaries
        for (i = 0; i < 2; i++)
            df[i] = (-3. * f[i] + 4. * f[i + 1] - f[i + 2]) / (2. * (x[i + 1] - x[i]));
        for (i = x.size()-3; i < x.size()-1; i++)
            df[i] = (3. * f[i] - 4. * f[i - 1] + f[i - 2]) / (2. * (x[i] - x[i - 1]));
    }
    else
        throw std::runtime_error("Not implemented");
}

#endif //SRC_BASE_STATE_H
