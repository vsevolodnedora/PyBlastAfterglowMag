//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_ODE_SOLVERS_H
#define SRC_ODE_SOLVERS_H

#include "pch.h"
#include "utils.h"
#include "logger.h"

#ifndef _restrict
# if defined(_MSC_VER) && _MSC_VER >= 1400
#  define _restrict __restrict
# elif defined(__GNUC__)
# define _restrict __restrict__
# else
#  define _restrict
# endif
#endif

#ifndef unused
# define unused(var) ((void) var)
#endif

namespace Detail {
    double constexpr sqrtNewtonRaphson(double x, double curr, double prev) {
        return curr == prev
               ? curr
               : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
}
/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
*/
double constexpr sqrt_c(double x)
{
    return x >= 0 && x < std::numeric_limits<double>::infinity()
           ? Detail::sqrtNewtonRaphson(x, x, 0)
           : std::numeric_limits<double>::quiet_NaN();
}

// internal storage for the solver
struct working_buffer {
    std::vector<char> data;

public:
    template<typename T = char>
    void ensure(std::size_t minimalSize) {
        if(data.size() < sizeof(T) * minimalSize)
            data.resize(sizeof(T) * minimalSize, 0.0);
    }

    template<typename T = char> T      * ptr()       { return reinterpret_cast<T*>( &data[0] ); }
    template<typename T = char> T const* ptr() const { return reinterpret_cast<T const*>( &data[0] );  }

    std::size_t size() const { return this->data.size(); }

//    ~working_buffer(){
//        std::cout << data << "\n";
//        std::destroy(data.begin(), data.end());
//    }
};

namespace Integrators {

    // Euler (1768) Institutiones Calculi Integralis
    struct euler_integrator {
        static const int stage = 1;
        static const int order = 1;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(size);
            double* _restrict knode = buffer.ptr<double>();
            f(knode, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * knode[i];

            time += h;
        }
    };

    // Midpoint method, improved Euler method, Runge secondary formula (1895)
    // - Neither modified / improved is written on wikipedia.
    // - http://www.mymathlib.com/diffeq/runge-kutta/improved_euler_method.html
    //   According to the midpoint method, it is the improved Euler's method.
    // - According to the Spectral method book, this is an improved Euler method.
    struct midpoint_integrator {
        static const int stage = 2;
        static const int order = 2;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict knode = buffer.ptr<double>();
            double* _restrict xnode = buffer.ptr<double>() + size;

            f(knode, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + 0.5 * h * knode[i];

            f(knode, size, time + 0.5 * h, xnode, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * knode[i];

            time += h;
        }
    };

    //
    //   Heun's method, modified Euler's method
    // - en.wikipedia has both improved / modified names_lc.
    // - According to the Spectral method book, this is called the modified Euler method.
    // - http://pc-physics.com/syuseieuler1.html This is also the modified Euler method.
    // - http://detail.chiebukuro.yahoo.co.jp/qa/question_detail/q1091336470 This page too
    struct heun_integrator {
        static const int stage = 2;
        static const int order = 2;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + h * k[i];
                value[i] += (1.0 / 2.0) * h * k[i];
            }

            f(k, size, time + h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (1.0 / 2.0) * h * k[i];

            time += h;
        }
    };

    struct ralston_integrator {
        static const int stage = 2;
        static const int order = 2;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + (2.0 / 3.0) * h * k[i];
                value[i] += (1.0 / 4.0) * h * k[i];
            }

            f(k, size, time + (2.0 / 3.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (3.0 / 4.0) * h * k[i];

            time += h;
        }
    };

    //---------------------------------------------------------------------------
    // Third formula

    /* Runge's 3rd order formula (inefficient because it is a 4th step 3rd order formula)
     * - Hiler's book
     */
    struct runge3_integrator {
        static const int stage = 4;
        static const int order = 3;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(3 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;
            double* _restrict y = buffer.ptr<double>() + size * 2;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + (1.0 / 2.0) * h * k[i];
                y[i] = value[i] + (1.0 / 6.0) * h * k[i];
            }

            f(k, size, time + (1.0 / 2.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + h * k[i];
                y[i] += (2.0 / 3.0) * h * k[i];
            }

            f(k, size, time + h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * k[i];

            f(k, size, time + h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] = y[i] + (1.0 / 6.0) * h * k[i];

            time += h;
        }
    };

    /* Heun cubic method
     * - Hiler's book
     * - http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
     * - http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v2_3.html
     *   Is introduced as "Runge-Kutta cubic method v2".
     */
    struct heun3_integrator {
        static const int stage = 3;
        static const int order = 3;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + (1.0 / 3.0) * h * k[i];
                value[i] += (1.0 / 4.0) * h * k[i];
            }

            f(k, size, time + (1.0 / 3.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                x[i] = 4.0 * value[i] - 3.0 * x[i] + (2.0 / 3.0) * h * k[i];

            f(k, size, time + (2.0 / 3.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (3.0 / 4.0) * h * k[i];

            time += h;
        }
    };

    // Ralston 3rd method
    // * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods according to,
    //   "Kutta's third-order method" is
    struct ralston3_integrator {
        static const int stage = 3;
        static const int order = 3;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + (1.0 / 2.0) * h * k[i];
                value[i] += (2.0 / 9.0) * h * k[i];
            }

            f(k, size, time + (1.0 / 2.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = (-4.0 / 5.0) * x[i] + (9.0 / 5.0) * value[i] + (3.0 / 4.0) * h * k[i];
                value[i] += (1.0 / 3.0) * h * k[i];
            }

            f(k, size, time + (3.0 / 4.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (4.0 / 9.0) * h * k[i];

            time += h;
        }
    };


    // Kutta third method, Classical RK3
    // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v1_3.html
    //   According to it, it is simply "Third-order method v1".
    // * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods according to,
    //   "Kutta's third-order method" である。
    // * http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau according to,
    // It is written as Kutta cubic method or Classical RK3.
    struct kutta3_integrator {
        static const int stage = 3;
        static const int order = 3;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k = buffer.ptr<double>();
            double* _restrict x = buffer.ptr<double>() + size;

            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + (1.0 / 2.0) * h * k[i];
                value[i] += (1.0 / 6.0) * h * k[i];
            }

            f(k, size, time + (1.0 / 2.0) * h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x[i] = (-7.0 / 2.0) * x[i] + (9.0 / 2.0) * value[i] + 2.0 * h * k[i];
                value[i] += (2.0 / 3.0) * h * k[i];
            }

            f(k, size, time + h, x, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (1.0 / 6.0) * h * k[i];

            time += h;
        }
    };

    //---------------------------------------------------------------------------
    // 4th dan 4th formula

    // RK4 (classical Runge-Kutta method)
    struct rk4_integrator {
        static const int stage = 4;
        static const int order = 4;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(3 * size);
            double* _restrict knode = buffer.ptr<double>();
            double* _restrict xnode = buffer.ptr<double>() + size;
            double* _restrict delta = buffer.ptr<double>() + size * 2;

//            double knode[m_size];
//            double xnode[m_size];
//            double delta[m_size];

            // k1
            f(knode, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                delta[i] = (1.0 / 6.0) * h * knode[i];
                xnode[i] = value[i] + 0.5 * h * knode[i];
            }

            // k2
            f(knode, size, time + 0.5 * h, xnode, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                delta[i] += (2.0 / 6.0) * h * knode[i];
                xnode[i] = value[i] + 0.5 * h * knode[i];
            }

            // k3
            f(knode, size, time + 0.5 * h, xnode, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                delta[i] += (2.0 / 6.0) * h * knode[i];
                xnode[i] = value[i] + h * knode[i];
            }

            // k4
            f(knode, size, time + h, xnode, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                double const a = delta[i] + (1.0 / 6.0) * h * knode[i];
                value[i] += a;
            }

            time += h;
        }
//        ~rk4_integrator(){std::cout << "deleting rk4..."<<"\n"; }
    };

    // Kutta 3/8-rule
    struct kutta_3_8_integrator {
        static const int stage = 4;
        static const int order = 4;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(3 * size);
            double* _restrict k  = buffer.ptr<double>();
            double* _restrict xi = buffer.ptr<double>() + size;
            double* _restrict x4 = buffer.ptr<double>() + size * 2;

            // k1
            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                xi[i] = value[i] + (1.0 / 3.0) * h * k[i];
                x4[i] = value[i] + h * k[i];
                value[i] += (1.0 / 8.0) * h * k[i];
            }

            // k2
            f(k, size, time + (1.0 / 3.0) * h, xi, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                xi[i] = 2.0 * xi[i] - x4[i] + h * k[i];
                x4[i] -= h * k[i];
                value[i] += (3.0 / 8.0) * h * k[i];
            }

            // k3
            f(k, size, time + (2.0 / 3.0) * h, xi, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                x4[i] += h * k[i];
                value[i] += (3.0 / 8.0) * h * k[i];
            }

            // k4
            f(k, size, time + h, x4, rhs_pars);
            for (std::size_t i = 0; i < size; i++)
                value[i] += (1.0 / 8.0) * h * k[i];

            time += h;
        }
    };

    // Runge-Kutta Gill method
    //   Ref. gill.1 According to it, there is a calculation procedure to properly correct the rounding error.
    //   [gill.1] [[Runge-Kutta-Gill About the law-Keisuke Araki's miscellaneous notes> http://d.hatena.ne.jp/arakik10/20091004/p1]]
    struct gill_integrator {
        static const int stage = 4;
        static const int order = 4;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(2 * size);
            double* _restrict k  = buffer.ptr<double>();
            double* _restrict q  = buffer.ptr<double>() + size;

            static constexpr double alpha1 = 1.0 / 2.0;
            static constexpr double alpha2 = 1.0 - sqrt_c(1.0 / 2.0);
            static constexpr double alpha3 = 1.0 + sqrt_c(1.0 / 2.0);
            static constexpr double alpha4 = 1.0 / 2.0;

            // k1
            f(k, size, time, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                double const y = value[i] + alpha1 * (h * k[i] - 2.0 * q[i]);

                // ※丸め誤差をキャンセルする為に r = ynew - yold とする必要がある。
                //   先に r を計算してから y に足すのでは駄目らしい。
                //   http://d.hatena.ne.jp/arakik10/20091004/p1
                //   http://ci.nii.ac.jp/naid/110002718589/
                double const r = y - value[i];

                q[i] = q[i] + 3.0 * r - alpha1 * h * k[i];
                value[i] = y;
            }

            // k2
            f(k, size, time + 0.5 * h, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                double const y = value[i] + alpha2 * (h * k[i] - q[i]);
                double const r = y - value[i];
                q[i] = q[i] + 3.0 * r - alpha2 * h * k[i];
                value[i] = y;
            }

            // k3
            f(k, size, time + 0.5 * h, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                double const y = value[i] + alpha3 * (h * k[i] - q[i]);
                double const r = y - value[i];
                q[i] = q[i] + 3.0 * r - alpha3 * h * k[i];
                value[i] = y;
            }

            // k4
            f(k, size, time + h, value, rhs_pars);
            for (std::size_t i = 0; i < size; i++) {
                double const y = value[i] + (alpha4 / 3.0) * (h * k[i] - 2.0 * q[i]);
                double const r = y - value[i];
                q[i] = q[i] + 3.0 * r - alpha4 * h * k[i]; // ※次のステップで使う
                value[i] = y;
            }

            time += h;
        }
    };


    //---------------------------------------------------------------------------
    // 6th dan 5th formula

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct butcher5v1_integrator {
        static const int stage = 6;
        static const int order = 5;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(6 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict  k3 = buffer.ptr<double>() + size * 3;
            double* _restrict  k4 = buffer.ptr<double>() + size * 4;
            double* _restrict  k5 = buffer.ptr<double>() + size * 5;
            double* _restrict& k6 = k2;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = 1.0 / 8.0;
            static constexpr double c2  = 1.0 / 8.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a21 * h * k1[i];
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            static constexpr double a32 = 1.0 / 4.0;
            static constexpr double c3  = 1.0 / 4.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a32 * h * k2[i];
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            static constexpr double a41 = -1.0 / 2.0;
            static constexpr double a42 =  1.0;
            static constexpr double c4  =  1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i]);
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            static constexpr double a51 = 15.0 / 16.0;
            static constexpr double a52 = -3.0 /  2.0;
            static constexpr double a53 =  3.0 /  4.0;
            static constexpr double a54 =  9.0 / 16.0;
            static constexpr double c5  =  3.0 /  4.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            static constexpr double a61 = -17.0 / 7.0;
            static constexpr double a62 =   4.0;
            static constexpr double a64 = -12.0 / 7.0;
            static constexpr double a65 =   8.0 / 7.0;
            static constexpr double c6  =   1.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a64 * k4[i] + a65 * k5[i]);
            f(k6, size, time + c6 * h, x, rhs_pars);

            // increment
            static constexpr double b1  =  7.0 / 90.0;
            static constexpr double b3  = 16.0 / 45.0;
            static constexpr double b4  =  2.0 / 15.0;
            static constexpr double b5  = 16.0 / 45.0;
            static constexpr double b6  =  7.0 / 90.0;
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);

            time += h;
        }
    };

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct butcher5v2_integrator {
        static const int stage = 6;
        static const int order = 5;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(5 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  y  = buffer.ptr<double>() + size * 1;
            double* _restrict  k1 = buffer.ptr<double>() + size * 2;
            double* _restrict  k2 = buffer.ptr<double>() + size * 3;
            double* _restrict& k3 = k2;
            double* _restrict& k4 = k2;
            double* _restrict& k5 = k2;
            double* _restrict  k6 = buffer.ptr<double>() + size * 4;

            static constexpr double a21 = 1.0 / 4.0;
            static constexpr double c2  = 1.0 / 4.0;

            static constexpr double a31 = 1.0 / 8.0;
            static constexpr double a32 = 1.0 / 8.0;
            static constexpr double c3  = 1.0 / 4.0;

            static constexpr double a42 = -1.0 / 2.0;
            static constexpr double a43 =  1.0;
            static constexpr double c4  =  1.0 / 2.0;

            static constexpr double a51 = 3.0 / 16.0;
            static constexpr double a54 = 9.0 / 16.0;
            static constexpr double c5  = 3.0 /  4.0;

            static constexpr double a61 =  -3.0 / 7.0;
            static constexpr double a62 =   2.0 / 7.0;
            static constexpr double a63 =  12.0 / 7.0;
            static constexpr double a64 = -12.0 / 7.0;
            static constexpr double a65 =   8.0 / 7.0;
            static constexpr double c6  =   1.0;

            static constexpr double b1  =  7.0 / 90.0;
            static constexpr double b3  = 16.0 / 45.0;
            static constexpr double b4  =  2.0 / 15.0;
            static constexpr double b5  = 16.0 / 45.0;
            static constexpr double b6  =  7.0 / 90.0;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            for (std::size_t i = 0; i < size; i++) {
                x[i]  = value[i] + a21 * h * k1[i];
                k6[i] = value[i] + a61 * h * k1[i];
                y[i]  = value[i] + b1 * h * k1[i];
            }
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
                k6[i] += a62 * h * k2[i];
            }
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            for (std::size_t i = 0; i < size; i++) {
                x[i] = a42 / a32 * x[i] + (1.0 - a42 / a32) * value[i] - ((a42 / a32) * a31) * h * k1[i] + a43 * h * k3[i];
                k6[i] += a63 * h * k3[i];
                y[i] += b3 * h * k3[i];
            }
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            for (std::size_t i = 0; i < size; i++) {
                x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
                k6[i] += a64 * h * k4[i];
                y[i] += b4 * h * k4[i];
            }
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            for (std::size_t i = 0; i < size; i++) {
                x[i] = k6[i] + a65 * h * k5[i];
                y[i] += b5 * h * k5[i];
            }
            f(k6, size, time + c6 * h, x, rhs_pars);

            // increment
            for (std::size_t i = 0; i < size; i++)
                value[i] = y[i] + b6 * h * k6[i];

            time += h;

        }
    };

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct butcher5v3_integrator {
        static const int stage = 6;
        static const int order = 5;
        mutable working_buffer buffer;

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(6 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict  k3 = buffer.ptr<double>() + size * 3;
            double* _restrict  k4 = buffer.ptr<double>() + size * 4;
            double* _restrict  k5 = buffer.ptr<double>() + size * 5;
            double* _restrict& k6 = k2;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = -1.0 / 2.0;
            static constexpr double c2  = -1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a21 * h * k1[i];
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            static constexpr double a31 =  5.0 / 16.0;
            static constexpr double a32 = -1.0 / 16.0;
            static constexpr double c3  =  1.0 /  4.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            static constexpr double a41 = -3.0 / 4.0;
            static constexpr double a42 =  1.0 / 4.0;
            static constexpr double a43 =  1.0;
            static constexpr double c4  =  1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            static constexpr double a51 = 3.0 / 16.0;
            static constexpr double a54 = 9.0 / 16.0;
            static constexpr double c5  = 3.0 /  4.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            static constexpr double a62 =  -1.0 / 7.0;
            static constexpr double a63 =  12.0 / 7.0;
            static constexpr double a64 = -12.0 / 7.0;
            static constexpr double a65 =   8.0 / 7.0;
            static constexpr double c6  = 1.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            f(k6, size, time + c6 * h, x, rhs_pars);

            // increment
            static constexpr double b1  =  7.0 / 90.0;
            static constexpr double b3  = 16.0 / 45.0;
            static constexpr double b4  =  2.0 / 15.0;
            static constexpr double b5  = 16.0 / 45.0;
            static constexpr double b6  =  7.0 / 90.0;
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);

            time += h;
        }
    };

    //---------------------------------------------------------------------------
    // Higher official

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct hammud6_integrator {
        static const int stage = 7;
        static const int order = 6;
        mutable working_buffer buffer;

        static constexpr double sqrt21 = sqrt_c(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(7 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict  k3 = buffer.ptr<double>() + size * 3;
            double* _restrict  k4 = buffer.ptr<double>() + size * 4;
            double* _restrict  k5 = buffer.ptr<double>() + size * 5;
            double* _restrict  k6 = buffer.ptr<double>() + size * 6;
            double* _restrict& k7 = k2;

            static constexpr double sqrt5 = sqrt_c(5);
            static constexpr double b1 = 1.0 / 12.0;
            static constexpr double b5 = 5.0 / 12.0;
            static constexpr double b6 = 5.0 / 12.0;
            static constexpr double b7 = 1.0 / 12.0;


            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = 4.0 / 7.0;
            static constexpr double c2  = 4.0 / 7.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a21 * h * k1[i];
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            static constexpr double a31 = 115.0 / 112.0;
            static constexpr double a32 =  -5.0 /  16.0;
            static constexpr double c3  =   5.0 /   7.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            static constexpr double a41 = 589.0 / 630.0;
            static constexpr double a42 =   5.0 /  18.0;
            static constexpr double a43 = -16.0 /  45.0;
            static constexpr double c4  =   6.0 /   7.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            static constexpr double a51 = 229.0 / 1200.0 -  29.0 * sqrt5 / 6000.0;
            static constexpr double a52 = 119.0 /  240.0 - 187.0 * sqrt5 / 1200.0;
            static constexpr double a53 = -14.0 /   75.0 +  34.0 * sqrt5 /  375.0;
            static constexpr double a54 =                   -3.0 * sqrt5 /  100.0;
            static constexpr double c5  = (5.0 - sqrt5) / 10.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            static constexpr double a61 =  71.0 / 2400.0 - 587.0 * sqrt5 / 12000.0;
            static constexpr double a62 = 187.0 / 480.0  - 391.0 * sqrt5 /  2400.0;
            static constexpr double a63 = -38.0 / 75.0   +  26.0 * sqrt5 /   375.0;
            static constexpr double a64 =  27.0 / 80.0   -   3.0 * sqrt5 /   400.0;
            static constexpr double a65 =   1.0 / 4.0    +         sqrt5 /     4.0;
            static constexpr double c6  = (5.0 + sqrt5) / 10.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            f(k6, size, time + c6 * h, x, rhs_pars);

            // k7 <= k2
            static constexpr double a71 = -49.0 / 480.0 + 43.0 * sqrt5 / 160.0;
            static constexpr double a72 = -425.0 / 96.0 + 51.0 * sqrt5 /  32.0;
            static constexpr double a73 =   52.0 / 15.0 -  4.0 * sqrt5 /   5.0;
            static constexpr double a74 =  -27.0 / 16.0 +  3.0 * sqrt5 /  16.0;
            static constexpr double a75 =    5.0 / 4.0  -  3.0 * sqrt5 /   4.0;
            static constexpr double a76 =    5.0 / 2.0  -        sqrt5 /   2.0;
            static constexpr double c7 = 1.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            f(k7, size, time + c7 * h, x, rhs_pars);

            // increment
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * (b1 * k1[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]);

            time += h;
        }
    };

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct shanks7_integrator {
        static const int stage = 9;
        static const int order = 7;
        mutable working_buffer buffer;

        static constexpr double sqrt21 = sqrt_c(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(9 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict& k3 = k2;
            double* _restrict  k4 = buffer.ptr<double>() + size * 4;
            double* _restrict  k5 = buffer.ptr<double>() + size * 5;
            double* _restrict  k6 = buffer.ptr<double>() + size * 6;
            double* _restrict  k7 = buffer.ptr<double>() + size * 7;
            double* _restrict  k8 = buffer.ptr<double>() + size * 8;
            double* _restrict& k9 = k2;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = 2.0 / 9.0;
            static constexpr double c2  = 2.0 / 9.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a21 * h * k1[i];
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            static constexpr double a31 = 1.0 / 12.0;
            static constexpr double a32 = 1.0 / 4.0;
            static constexpr double c3  = 1.0 / 3.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            static constexpr double a41 = 1.0 / 8.0;
            static constexpr double a43 = 3.0 / 8.0;
            static constexpr double c4  = 1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            static constexpr double a51 = 23.0 / 216.0;
            static constexpr double a53 =  7.0 / 72.0;
            static constexpr double a54 = -1.0 / 27.0;
            static constexpr double c5  =  1.0 / 6.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            static constexpr double a61 = -4136.0 / 729.0;
            static constexpr double a63 = -4528.0 / 243.0;
            static constexpr double a64 =  5264.0 / 729.0;
            static constexpr double a65 =  1456.0 / 81.0;
            static constexpr double c6  =     8.0 / 9.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            f(k6, size, time + c6 * h, x, rhs_pars);

            // k7
            static constexpr double a71 = 8087.0 / 11664.0;
            static constexpr double a73 =  484.0 / 243.0;
            static constexpr double a74 = -518.0 / 729.0;
            static constexpr double a75 = -658.0 / 351.0;
            static constexpr double a76 =    7.0 / 624.0;
            static constexpr double c7  =    1.0 / 9.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            f(k7, size, time + c7 * h, x, rhs_pars);

            // k8
            static constexpr double a81 = -1217.0 / 2160.0;
            static constexpr double a83 =  -145.0 / 72.0;
            static constexpr double a84 =  8342.0 / 6615.0;
            static constexpr double a85 =   361.0 / 195.0;
            static constexpr double a86 =  3033.0 / 50960.0;
            static constexpr double a87 =   117.0 / 490.0;
            static constexpr double c8  =     5.0 / 6.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a81 * k1[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
            f(k8, size, time + c8 * h, x, rhs_pars);

            // k9
            static constexpr double a91 =    259.0 / 2768.0;
            static constexpr double a93 =    -84.0 / 173.0;
            static constexpr double a94 =    -14.0 / 173.0;
            static constexpr double a95 =   6210.0 / 2249.0;
            static constexpr double a96 = -99873.0 / 251888.0;
            static constexpr double a97 = -29160.0 / 15743.0;
            static constexpr double a98 =   2160.0 / 2249.0;
            static constexpr double c9  = 1.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
            f(k9, size, time + c9 * h, x, rhs_pars);

            // increment
            static constexpr double b1 =    173.0 / 3360.0;
            static constexpr double b4 =   1846.0 / 5145.0;
            static constexpr double b5 =     27.0 / 91.0;
            static constexpr double b6 = -19683.0 / 713440.0;
            static constexpr double b7 = -19683.0 / 713440.0;
            static constexpr double b8 =     27.0 / 91.0;
            static constexpr double b9 =    173.0 / 3360.0;
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * (b1 * k1[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i]);

            time += h;
        }
    };

    // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    struct cooper_verner7_integrator {
        static const int stage = 9;
        static const int order = 7;
        mutable working_buffer buffer;

        static constexpr double sqrt21 = sqrt_c(21.0);

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(9 * size);
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict& k3 = k2;
            double* _restrict  k4 = buffer.ptr<double>() + size * 4;
            double* _restrict  k5 = buffer.ptr<double>() + size * 5;
            double* _restrict  k6 = buffer.ptr<double>() + size * 6;
            double* _restrict  k7 = buffer.ptr<double>() + size * 7;
            double* _restrict  k8 = buffer.ptr<double>() + size * 8;
            double* _restrict& k9 = k2;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = (7.0 + 1.0 * sqrt21) / 42.0;
            static constexpr double c2  = (7.0 +       sqrt21) / 42.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a21 * h * k1[i];
            f(k2, size, time + c2 * h, x, rhs_pars);

            // k3
            static constexpr double a32 = (7.0 + 1.0 * sqrt21) / 21.0;
            static constexpr double c3  = (7.0 +       sqrt21) / 21.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + a32 * h * k2[i];
            f(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            static constexpr double a41 = ( 7.0 + 1.0 * sqrt21) / 56.0;
            static constexpr double a43 = (21.0 + 3.0 * sqrt21) / 56.0;
            static constexpr double c4  = ( 7.0 +       sqrt21) / 14.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
            f(k4, size, time + c4 * h, x, rhs_pars);

            // k5
            static constexpr double a51 = (  8.0 - 1.0 * sqrt21) / 16.0;
            static constexpr double a53 = (-21.0 + 6.0 * sqrt21) / 16.0;
            static constexpr double a54 = ( 21.0 - 5.0 * sqrt21) / 16.0;
            static constexpr double c5 = 1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
            f(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            static constexpr double a61 = (-1687.0 + 374.0 * sqrt21) / 196.0;
            static constexpr double a63 = (  969.0 - 210.0 * sqrt21) /  28.0;
            static constexpr double a64 = ( -381.0 +  83.0 * sqrt21) /  14.0;
            static constexpr double a65 = (   84.0 -  20.0 * sqrt21) /  49.0;
            static constexpr double c6  = (    7.0 -         sqrt21) /  14.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            f(k6, size, time + c6 * h, x, rhs_pars);

            // k7
            static constexpr double a71 = (  583.0 - 131.0 * sqrt21) / 128.0;
            static constexpr double a73 = (-2373.0 + 501.0 * sqrt21) / 128.0;
            static constexpr double a74 = ( 4221.0 - 914.0 * sqrt21) / 288.0;
            static constexpr double a75 = (   -9.0 +   4.0 * sqrt21) /  18.0;
            static constexpr double a76 = (  189.0 +  35.0 * sqrt21) / 576.0;
            static constexpr double c7 = 1.0 / 2.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            f(k7, size, time + c7 * h, x, rhs_pars);

            // k8
            static constexpr double a81 = ( -623.0 +  169.0 * sqrt21) /  392.0;
            static constexpr double a83 = (  435.0 -   81.0 * sqrt21) /   56.0;
            static constexpr double a84 = (-1437.0 +  307.0 * sqrt21) /  252.0;
            static constexpr double a85 = (-2028.0 - 1468.0 * sqrt21) / 7497.0;
            static constexpr double a86 = (  -21.0 -    4.0 * sqrt21) /  126.0;
            static constexpr double a87 = (  384.0 +   80.0 * sqrt21) /  833.0;
            static constexpr double c8  = (    7.0 +          sqrt21) /   14.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a81 * k1[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
            f(k8, size, time + c8 * h, x, rhs_pars);

            // k9
            static constexpr double a91 = (  579.0 -  131.0 * sqrt21) /  24.0;
            static constexpr double a93 = ( -791.0 +  167.0 * sqrt21) /   8.0;
            static constexpr double a94 = ( 8099.0 - 1765.0 * sqrt21) / 108.0;
            static constexpr double a95 = (-1976.0 +  784.0 * sqrt21) / 459.0;
            static constexpr double a96 = (   70.0 +    7.0 * sqrt21) /  54.0;
            static constexpr double a97 = (  160.0 -   80.0 * sqrt21) / 153.0;
            static constexpr double a98 = (   49.0 -    7.0 * sqrt21) /  18.0;
            static constexpr double c9 = 1.0;
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
            f(k9, size, time + c9 * h, x, rhs_pars);

            // increment
            static constexpr double b1 =  1.0 /  20.0;
            static constexpr double b6 = 49.0 / 180.0;
            static constexpr double b7 = 16.0 /  45.0;
            static constexpr double b8 = 49.0 / 180.0;
            static constexpr double b9 =  1.0 /  20.0;
            for (std::size_t i = 0; i < size; i++)
                value[i] += h * (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i]);

            time += h;
        }
    };

    // Cooper-Verner 8th order
    //   According to Ref. Cv8.3, it is named "Verner's 8th order method".
    //   In Ref. cv8.2, the sign before sqrt (21) is inverted. Either one seems to be fine.
    //   [cv8.1] E. Heiler, Taketomo Mitsui, "Numerical Solution of Common Differential Equations I", Maruzen Publishing (2012/7/17).
    //   [cv8.2] http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
    //   [cv8.3] http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_verner.html
    struct cooper_verner8_integrator {
        static const int stage = 11;
        static const int order = 8;
        mutable working_buffer buffer;

        static constexpr double sqrt21 = sqrt_c(21.0); // Ref [cv8.2] Then -sqrt(21.0). Either OK.

        template<typename F>
        void operator()(double& time, double* _restrict value, std::size_t size, F const& f, void * rhs_pars, double h) const {
            buffer.ensure<double>(8 * size);
            double* _restrict  delta = buffer.ptr<double>();
            double* _restrict  xnode = buffer.ptr<double>() + size;
            double* _restrict  k1    = buffer.ptr<double>() + size * 2;
            double* _restrict  k2    = buffer.ptr<double>() + size * 3;
            double* _restrict  k3    = buffer.ptr<double>() + size * 4;
            double* _restrict& k4    = k2;
            double* _restrict  k5    = buffer.ptr<double>() + size * 5;
            double* _restrict  k6    = buffer.ptr<double>() + size * 6;
            double* _restrict& k7    = k3;
            double* _restrict& k8    = k4;
            double* _restrict  k9    = buffer.ptr<double>() + size * 7;
            double* _restrict& kA    = k1;
            double* _restrict& kB    = k5;

            static constexpr double b1 =  1.0 /  20.0;
            static constexpr double b8 = 49.0 /  180.0;
            static constexpr double b9 = 16.0 /  45.0;
            static constexpr double bA = 49.0 / 180.0;
            static constexpr double bB =  1.0 /  20.0;

            // k1
            f(k1, size, time, value, rhs_pars);

            // k2
            static constexpr double a21 = 0.5;
            static constexpr double c20 = a21;
            for (std::size_t i = 0; i < size; i++) {
                delta[i] = b1 * h * k1[i];
                xnode[i] = value[i] + a21 * h * k1[i];
            }
            f(k2, size, time + c20 * h, xnode, rhs_pars);

            // k3
            static constexpr double a31 = 0.25;
            static constexpr double a32 = 0.25;
            static constexpr double c30 = a31 + a32;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a31 * h * k1[i] + a32 * h * k2[i];
            f(k3, size, time + c30 * h, xnode, rhs_pars);

            // k4 <= k2
            static constexpr double a41 = (1.0 /  7.0);
            static constexpr double a42 = (1.0 / 98.0) * (-7 - 3 * sqrt21);
            static constexpr double a43 = (1.0 / 49.0) * (21 + 5 * sqrt21);
            static constexpr double c40 = a41 + a42 + a43;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a41 * h * k1[i] + a42 * h * k2[i] + a43 * h * k3[i];
            f(k4, size, time + c40 * h, xnode, rhs_pars);

            // k5
            static constexpr double a51 = (1.0 /  84.0) * (11 + 1 * sqrt21);
            static constexpr double a53 = (1.0 /  63.0) * (18 + 4 * sqrt21);
            static constexpr double a54 = (1.0 / 252.0) * (21 - 1 * sqrt21);
            static constexpr double c50 = a51 + a53 + a54;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a51 * h * k1[i] + a53 * h * k3[i] + a54 * h * k4[i];
            f(k5, size, time + c50 * h, xnode, rhs_pars);

            // k6
            static constexpr double a61 = (1.0 /  48.0) * (   5 +  1 * sqrt21);
            static constexpr double a63 = (1.0 /  36.0) * (   9 +  1 * sqrt21);
            static constexpr double a64 = (1.0 / 360.0) * (-231 + 14 * sqrt21);
            static constexpr double a65 = (1.0 /  80.0) * (  63 -  7 * sqrt21);
            static constexpr double c60 = a61 + a63 + a64 + a65;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a61 * h * k1[i] + a63 * h * k3[i] + a64 * h * k4[i] + a65 * h * k5[i];
            f(k6, size, time + c60 * h, xnode, rhs_pars);

            // k7 <= k3
            static constexpr double a71 = (1.0 /  42.0) * (  10 -   1 * sqrt21);
            static constexpr double a73 = (1.0 / 315.0) * (-432 +  92 * sqrt21);
            static constexpr double a74 = (1.0 /  90.0) * ( 633 - 145 * sqrt21);
            static constexpr double a75 = (1.0 /  70.0) * (-504 + 115 * sqrt21);
            static constexpr double a76 = (1.0 /  35.0) * (  63 -  13 * sqrt21);
            static constexpr double c70 = a71 + a73 + a74 + a75 + a76;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a71 * h * k1[i] + a73 * h * k3[i] + a74 * h * k4[i] + a75 * h * k5[i] + a76 * h * k6[i];
            f(k7, size, time + c70 * h, xnode, rhs_pars);

            // k8 <= k4
            static constexpr double a81 =  1.0 /  14.0;
            static constexpr double a85 = (1.0 / 126.0) * (14 - 3 * sqrt21);
            static constexpr double a86 = (1.0 /  63.0) * (13 - 3 * sqrt21);
            static constexpr double a87 =  1.0 /   9.0;
            static constexpr double c80 = a81 + a85 + a86 + a87;
            for (std::size_t i = 0; i < size; i++)
                xnode[i] = value[i] + a81 * h * k1[i] + a85 * h * k5[i] + a86 * h * k6[i] + a87 * h * k7[i];
            f(k8, size, time + c80 * h, xnode, rhs_pars);

            // k9
            static constexpr double a91 =   1.0 /   32.0;
            static constexpr double a95 = ( 1.0 /  576.0) * (  91.0 - 21.0 * sqrt21);
            static constexpr double a96 =  11.0 /   72.0;
            static constexpr double a97 = ( 1.0 / 1152.0) * (-385.0 - 75.0 * sqrt21);
            static constexpr double a98 = ( 1.0 /  128.0) * (  63.0 + 13.0 * sqrt21);
            static constexpr double c90 = a91 + a95 + a96 + a97 + a98;
            for (std::size_t i = 0; i < size; i++) {
                delta[i] += b8 * h * k8[i];
                xnode[i] = value[i] + a91 * h * k1[i] + a95 * h * k5[i] + a96 * h * k6[i] + a97 * h * k7[i] + a98 * h * k8[i];
            }
            f(k9, size, time + c90 * h, xnode, rhs_pars);

            // kA <= k1
            static constexpr double aA1 =  1.0 /   14.0;
            static constexpr double aA5 =  1.0 /    9.0;
            static constexpr double aA6 = (1.0 / 2205.0) * (-733.0 - 147.0 * sqrt21);
            static constexpr double aA7 = (1.0 /  504.0) * ( 515.0 + 111.0 * sqrt21);
            static constexpr double aA8 = (1.0 /   56.0) * (- 51.0 -  11.0 * sqrt21);
            static constexpr double aA9 = (1.0 /  245.0) * ( 132.0 +  28.0 * sqrt21);
            static constexpr double cA0 = aA1 + aA5 + aA6 + aA7 + aA8 + aA9;
            for (std::size_t i = 0; i < size; i++) {
                delta[i] += b9 * h * k9[i];
                xnode[i] = value[i] + aA1 * h * k1[i] + aA5 * h * k5[i] + aA6 * h * k6[i] + aA7 * h * k7[i] + aA8 * h * k8[i] + aA9 * h * k9[i];
            }
            f(kA, size, time + cA0 * h, xnode, rhs_pars);

            // kB <= k5
            static constexpr double aB5 = (1.0 / 18.0) * (- 42.0 +  7.0 * sqrt21);
            static constexpr double aB6 = (1.0 / 45.0) * (- 18.0 + 28.0 * sqrt21);
            static constexpr double aB7 = (1.0 / 72.0) * (-273.0 - 53.0 * sqrt21);
            static constexpr double aB8 = (1.0 / 72.0) * ( 301.0 + 53.0 * sqrt21);
            static constexpr double aB9 = (1.0 / 45.0) * (  28.0 - 28.0 * sqrt21);
            static constexpr double aBA = (1.0 / 18.0) * (  49.0 -  7.0 * sqrt21);
            static constexpr double cB0 = aB5 + aB6 + aB7 + aB8 + aB9 + aBA;
            for (std::size_t i = 0; i < size; i++) {
                delta[i] += bA * h * kA[i];
                xnode[i] = value[i] + aB5 * h * k5[i] + aB6 * h * k6[i] + aB7 * h * k7[i] + aB8 * h * k8[i] + aB9 * h * k9[i] + aBA * h * kA[i];
            }
            f(kB, size, time + cB0 * h, xnode, rhs_pars);

            // increment
            for (std::size_t i = 0; i < size; i++) {
                double const a = delta[i] + bB * h * kB[i];
                value[i] += a;
            }

            time += h;
        }
    };

    // TODO
    // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_4.html Ralston's 4th
    // * Is it the same as Hyler's Ralson (1962), Hull (1967) (u, vv) = (0.4, 0.45)?
    // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_nystrom.html Nystrom's 5th
    // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_butcher.html Butcher's 6th


    /* ==================================================================================== */

    // one of the utils
    template<typename T>
    T const& clamp(T const& value, T const& lowerBound, T const& upperBound) {
        return value < lowerBound ? lowerBound : value > upperBound ? upperBound : value;
    }


    struct iequation_for_erk;

    struct stat_t {
        int nfcn   {0};
        int nstep  {0};
        int naccpt {0};
        int nrejct {0};

        iequation_for_erk* eq            {nullptr};
        double const*      previousValue {nullptr};
        std::size_t        previousSize  {0};
        double             previousTime  {0.0};
        double             previousStep  {0.0};
//            void *             rhs_pars      {nullptr};


        void clear() {
            nfcn   = 0;
            nstep  = 0;
            naccpt = 0;
            nrejct = 0;
        }
    };

    struct idense_output {

        virtual void get_values_at(
                double* _restrict interpolated, double time,
                std::size_t const* icomp, std::size_t ncomp
        ) = 0;

        void get_values_at(double* _restrict interpolated, double time) {
            return get_values_at(interpolated, time, nullptr, 0);
        }
    };

    /* used as a type for RHS function */
    struct iequation_for_erk {
        virtual ~iequation_for_erk() { }
        virtual void eval_derivative(double* _restrict slope,
                                     std::size_t size,
                                     double t,
                                     double const* _restrict value,
                                     void * rhs_pars) = 0;
        virtual void onstep() {}
        virtual void ondense( stat_t const&, idense_output& ) { }

    public:
        virtual bool is_stopping() { return false; }
    };


    // DOP853: Dormand-Prince 8(5, 3)
    //
    //   The coefficients are taken from the supplementary material of the Heirer's book:
    //   from http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
    //
    struct dop853_integrator {
        static const int stage = 12;
        static const int order = 8;
        mutable working_buffer buffer;

        //   buffer State
        //   before [  ? | k1 |  ? |  ? |  ? |  ? |  ? |  ? |  ? |  ? ]
        //   after  [  x | k1 | k5 | kC | x6 | x7 | x8 | k9 | kA | kB ]
        //   However, x is the value of the next step.
        //   xE, xF, xG are the variables used internally to generate the honey output data.
        //   This function can only be called once, as calling this function destroys internal information.
        void _integrate8(
                double& time, double* _restrict value, std::size_t size,
                iequation_for_erk& eq, void * rhs_pars, double h,
                double atol, double rtol, double& _err, double& _stf
        ) const {
            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k2 = buffer.ptr<double>() + size * 2;
            double* _restrict& k3 = k2;
            double* _restrict  k4 = buffer.ptr<double>() + size * 3;
            double* _restrict& k5 = k3;
            double* _restrict  k6 = buffer.ptr<double>() + size * 4;
            double* _restrict  k7 = buffer.ptr<double>() + size * 5;
            double* _restrict  k8 = buffer.ptr<double>() + size * 6;
            double* _restrict  k9 = buffer.ptr<double>() + size * 7;
            double* _restrict  kA = buffer.ptr<double>() + size * 8;
            double* _restrict  kB = buffer.ptr<double>() + size * 9;
            double* _restrict& kC = k4;

            static constexpr double c2 = 0.526001519587677318785587544488E-01;
            static constexpr double c3 = 0.789002279381515978178381316732E-01;
            static constexpr double c4 = 0.118350341907227396726757197510E+00;
            static constexpr double c5 = 0.281649658092772603273242802490E+00;
            static constexpr double c6 = 0.333333333333333333333333333333E+00;
            static constexpr double c7 = 0.25E+00;
            static constexpr double c8 = 0.307692307692307692307692307692E+00;
            static constexpr double c9 = 0.651282051282051282051282051282E+00;
            static constexpr double cA = 0.6E+00;
            static constexpr double cB = 0.857142857142857142857142857142E+00;
            static constexpr double cC = 1.0;

            static constexpr double a21 =  5.26001519587677318785587544488E-2;
            static constexpr double a31 =  1.97250569845378994544595329183E-2;
            static constexpr double a32 =  5.91751709536136983633785987549E-2;
            static constexpr double a41 =  2.95875854768068491816892993775E-2;
            static constexpr double a43 =  8.87627564304205475450678981324E-2;
            static constexpr double a51 =  2.41365134159266685502369798665E-1;
            static constexpr double a53 = -8.84549479328286085344864962717E-1;
            static constexpr double a54 =  9.24834003261792003115737966543E-1;
            static constexpr double a61 =  3.7037037037037037037037037037E-2;
            static constexpr double a64 =  1.70828608729473871279604482173E-1;
            static constexpr double a65 =  1.25467687566822425016691814123E-1;
            static constexpr double a71 =  3.7109375E-2;
            static constexpr double a74 =  1.70252211019544039314978060272E-1;
            static constexpr double a75 =  6.02165389804559606850219397283E-2;
            static constexpr double a76 = -1.7578125E-2;
            static constexpr double a81 =  3.70920001185047927108779319836E-2;
            static constexpr double a84 =  1.70383925712239993810214054705E-1;
            static constexpr double a85 =  1.07262030446373284651809199168E-1;
            static constexpr double a86 = -1.53194377486244017527936158236E-2;
            static constexpr double a87 =  8.27378916381402288758473766002E-3;
            static constexpr double a91 =  6.24110958716075717114429577812E-1;
            static constexpr double a94 = -3.36089262944694129406857109825;
            static constexpr double a95 = -8.68219346841726006818189891453E-1;
            static constexpr double a96 =  2.75920996994467083049415600797E+1;
            static constexpr double a97 =  2.01540675504778934086186788979E+1;
            static constexpr double a98 = -4.34898841810699588477366255144E+1;
            static constexpr double aA1 =  4.77662536438264365890433908527E-1;
            static constexpr double aA4 = -2.48811461997166764192642586468;
            static constexpr double aA5 = -5.90290826836842996371446475743E-1;
            static constexpr double aA6 =  2.12300514481811942347288949897E+1;
            static constexpr double aA7 =  1.52792336328824235832596922938E+1;
            static constexpr double aA8 = -3.32882109689848629194453265587E+1;
            static constexpr double aA9 = -2.03312017085086261358222928593E-2;
            static constexpr double aB1 = -9.3714243008598732571704021658E-1;
            static constexpr double aB4 =  5.18637242884406370830023853209;
            static constexpr double aB5 =  1.09143734899672957818500254654;
            static constexpr double aB6 = -8.14978701074692612513997267357;
            static constexpr double aB7 = -1.85200656599969598641566180701E+1;
            static constexpr double aB8 =  2.27394870993505042818970056734E+1;
            static constexpr double aB9 =  2.49360555267965238987089396762;
            static constexpr double aBA = -3.0467644718982195003823669022;
            static constexpr double aC1 =  2.27331014751653820792359768449;
            static constexpr double aC4 = -1.05344954667372501984066689879E+1;
            static constexpr double aC5 = -2.00087205822486249909675718444;
            static constexpr double aC6 = -1.79589318631187989172765950534E+1;
            static constexpr double aC7 =  2.79488845294199600508499808837E+1;
            static constexpr double aC8 = -2.85899827713502369474065508674;
            static constexpr double aC9 = -8.87285693353062954433549289258;
            static constexpr double aCA =  1.23605671757943030647266201528E+1;
            static constexpr double aCB =  6.43392746015763530355970484046E-1;

            static constexpr double b1 =  5.42937341165687622380535766363E-2;
            static constexpr double b6 =  4.45031289275240888144113950566;
            static constexpr double b7 =  1.89151789931450038304281599044;
            static constexpr double b8 = -5.8012039600105847814672114227;
            static constexpr double b9 =  3.1116436695781989440891606237E-1;
            static constexpr double bA = -1.52160949662516078556178806805E-1;
            static constexpr double bB =  2.01365400804030348374776537501E-1;
            static constexpr double bC =  4.47106157277725905176885569043E-2;

            static constexpr double er1 =  0.1312004499419488073250102996E-01;
            static constexpr double er6 = -0.1225156446376204440720569753E+01;
            static constexpr double er7 = -0.4957589496572501915214079952E+00;
            static constexpr double er8 =  0.1664377182454986536961530415E+01;
            static constexpr double er9 = -0.3503288487499736816886487290E+00;
            static constexpr double erA =  0.3341791187130174790297318841E+00;
            static constexpr double erB =  0.8192320648511571246570742613E-01;
            static constexpr double erC = -0.2235530786388629525884427845E-01;

            static constexpr double bhh1 = 0.244094488188976377952755905512E+00;
            static constexpr double bhh9 = 0.733846688281611857341361741547E+00;
            static constexpr double bhhC = 0.220588235294117647058823529412E-01;

            // k1: Set before calling this function (FSAL)

            // k2
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a21 * k1[i]);
            eq.eval_derivative(k2, size, time + c2 * h, x, rhs_pars);

            // k3 := k2
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
            eq.eval_derivative(k3, size, time + c3 * h, x, rhs_pars);

            // k4
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
            eq.eval_derivative(k4, size, time + c4 * h, x, rhs_pars);

            // k5 := k3
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
            eq.eval_derivative(k5, size, time + c5 * h, x, rhs_pars);

            // k6
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a61 * k1[i] + a64 * k4[i] + a65 * k5[i]);
            eq.eval_derivative(k6, size, time + c6 * h, x, rhs_pars);

            // k7
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            eq.eval_derivative(k7, size, time + c7 * h, x, rhs_pars);

            // k8
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
            eq.eval_derivative(k8, size, time + c8 * h, x, rhs_pars);

            // k9
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
            eq.eval_derivative(k9, size, time + c9 * h, x, rhs_pars);

            // kA
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (aA1 * k1[i] + aA4 * k4[i] + aA5 * k5[i] + aA6 * k6[i] + aA7 * k7[i] + aA8 * k8[i] + aA9 * k9[i]);
            eq.eval_derivative(kA, size, time + cA * h, x, rhs_pars);

            // kB
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (aB1 * k1[i] + aB4 * k4[i] + aB5 * k5[i] + aB6 * k6[i] + aB7 * k7[i] + aB8 * k8[i] + aB9 * k9[i] + aBA * kA[i]);
            eq.eval_derivative(kB, size, time + cB * h, x, rhs_pars);

            // kC
            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + h * (aC1 * k1[i] + aC4 * k4[i] + aC5 * k5[i] + aC6 * k6[i] + aC7 * k7[i] + aC8 * k8[i] + aC9 * k9[i] + aCA * kA[i] + aCB * kB[i]);
            eq.eval_derivative(kC, size, time + cC * h, x, rhs_pars);

            // increment
            double err1 = 0.0, err2 = 0.0, err3 = 0.0;
            for (std::size_t i = 0; i < size; i++) {
                double const slope = b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + bA * kA[i] + bB * kB[i] + bC * kC[i];
                double const _y = value[i] + h * slope;
                double const sk = atol + rtol * std::max(std::abs(value[i]), std::abs(_y));
                double const e1 = er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] + erA * kA[i] + erB * kB[i] + erC * kC[i];
                double const e2 = slope - bhh1 * k1[i] - bhh9 * k9[i] - bhhC * kC[i];

                err1 += (e1 / sk) * (e1 / sk);
                err2 += (e2 / sk) * (e2 / sk);
                err3 += slope * slope;
                x[i] = _y;
            }

            double deno = err1 + 0.01 * err2;
            if (deno <= 0.0) deno = 1.0;
            _err = std::abs(h) * err1 / std::sqrt(size * deno);
            _stf = h * h * err3;

            // std::printf("h = %g, tol = (%g %g), (|slope|, |err1|, |err2|) = (%g, %g, %g), _err=%g\n",
            //   h, atol, rtol, std::sqrt(err3 / m_size), std::sqrt(err1 / m_size), std::sqrt(err2 / m_size),_err);
        };
//        ~dop853_integrator(){ delete[] buffer; }

        struct param_t {
            double atol {1e-13};
            double rtol {1e-13};
            double beta {0.0};
            double fac1 {0.0};
            double fac2 {0.0};
            double safe {0.0};
            double hmax {0.0};
            double step {0.0};
            int    nstif {0};
            std::ptrdiff_t nmax {100000};
            void * rhs_pars {nullptr};
            int verbocity = -1;
        };


        void _dense_output_initialize(
                working_buffer& interpBuffer, stat_t& stat, std::size_t const* icomp, std::size_t nrd,
                void * rhs_pars) const {
            double* _restrict             cont  = interpBuffer.ptr<double>();
            iequation_for_erk&               eq    = *stat.eq;
            double const                     time  = stat.previousTime;
            double const* _restrict const value = stat.previousValue;
            std::size_t const                size  = stat.previousSize;
            double const                     h     = stat.previousStep;
//                void *                           rhs_ptr = stat.rhs_pars;

            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            double* _restrict  k6 = buffer.ptr<double>() + size * 4;
            double* _restrict  k7 = buffer.ptr<double>() + size * 5;
            double* _restrict  k8 = buffer.ptr<double>() + size * 6;
            double* _restrict  k9 = buffer.ptr<double>() + size * 7;
            double* _restrict  kA = buffer.ptr<double>() + size * 8;
            double* _restrict  kB = buffer.ptr<double>() + size * 9;
            double* _restrict  kC = buffer.ptr<double>() + size * 3;
            double* _restrict  kD = buffer.ptr<double>() + size * 2;

            double* _restrict& xE = k6;
            double* _restrict& xF = k7;
            double* _restrict& xG = k8;
            double* _restrict& kE = k9;
            double* _restrict& kF = kA;
            double* _restrict& kG = kB;

            static constexpr double cE  = 0.1E+00;
            static constexpr double cF  = 0.2E+00;
            static constexpr double cG  = 0.777777777777777777777777777778E+00;

            static constexpr double aE1 =  5.61675022830479523392909219681E-2;
            static constexpr double aE7 =  2.53500210216624811088794765333E-1;
            static constexpr double aE8 = -2.46239037470802489917441475441E-1;
            static constexpr double aE9 = -1.24191423263816360469010140626E-1;
            static constexpr double aEA =  1.5329179827876569731206322685E-1;
            static constexpr double aEB =  8.20105229563468988491666602057E-3;
            static constexpr double aEC =  7.56789766054569976138603589584E-3;
            static constexpr double aED = -8.298E-3;
            static constexpr double aF1 =  3.18346481635021405060768473261E-2;
            static constexpr double aF6 =  2.83009096723667755288322961402E-2;
            static constexpr double aF7 =  5.35419883074385676223797384372E-2;
            static constexpr double aF8 = -5.49237485713909884646569340306E-2;
            static constexpr double aFB = -1.08347328697249322858509316994E-4;
            static constexpr double aFC =  3.82571090835658412954920192323E-4;
            static constexpr double aFD = -3.40465008687404560802977114492E-4;
            static constexpr double aFE =  1.41312443674632500278074618366E-1;
            static constexpr double aG1 = -4.28896301583791923408573538692E-1;
            static constexpr double aG6 = -4.69762141536116384314449447206;
            static constexpr double aG7 =  7.68342119606259904184240953878;
            static constexpr double aG8 =  4.06898981839711007970213554331;
            static constexpr double aG9 =  3.56727187455281109270669543021E-1;
            static constexpr double aGD = -1.39902416515901462129418009734E-3;
            static constexpr double aGE =  2.9475147891527723389556272149;
            static constexpr double aGF = -9.15095847217987001081870187138;

            static constexpr double d41 = -0.84289382761090128651353491142E+01;
            static constexpr double d46 =  0.56671495351937776962531783590E+00;
            static constexpr double d47 = -0.30689499459498916912797304727E+01;
            static constexpr double d48 =  0.23846676565120698287728149680E+01;
            static constexpr double d49 =  0.21170345824450282767155149946E+01;
            static constexpr double d4A = -0.87139158377797299206789907490E+00;
            static constexpr double d4B =  0.22404374302607882758541771650E+01;
            static constexpr double d4C =  0.63157877876946881815570249290E+00;
            static constexpr double d4D = -0.88990336451333310820698117400E-01;
            static constexpr double d4E =  0.18148505520854727256656404962E+02;
            static constexpr double d4F = -0.91946323924783554000451984436E+01;
            static constexpr double d4G = -0.44360363875948939664310572000E+01;

            static constexpr double d51 =  0.10427508642579134603413151009E+02;
            static constexpr double d56 =  0.24228349177525818288430175319E+03;
            static constexpr double d57 =  0.16520045171727028198505394887E+03;
            static constexpr double d58 = -0.37454675472269020279518312152E+03;
            static constexpr double d59 = -0.22113666853125306036270938578E+02;
            static constexpr double d5A =  0.77334326684722638389603898808E+01;
            static constexpr double d5B = -0.30674084731089398182061213626E+02;
            static constexpr double d5C = -0.93321305264302278729567221706E+01;
            static constexpr double d5D =  0.15697238121770843886131091075E+02;
            static constexpr double d5E = -0.31139403219565177677282850411E+02;
            static constexpr double d5F = -0.93529243588444783865713862664E+01;
            static constexpr double d5G =  0.35816841486394083752465898540E+02;

            static constexpr double d61 =  0.19985053242002433820987653617E+02;
            static constexpr double d66 = -0.38703730874935176555105901742E+03;
            static constexpr double d67 = -0.18917813819516756882830838328E+03;
            static constexpr double d68 =  0.52780815920542364900561016686E+03;
            static constexpr double d69 = -0.11573902539959630126141871134E+02;
            static constexpr double d6A =  0.68812326946963000169666922661E+01;
            static constexpr double d6B = -0.10006050966910838403183860980E+01;
            static constexpr double d6C =  0.77771377980534432092869265740E+00;
            static constexpr double d6D = -0.27782057523535084065932004339E+01;
            static constexpr double d6E = -0.60196695231264120758267380846E+02;
            static constexpr double d6F =  0.84320405506677161018159903784E+02;
            static constexpr double d6G =  0.11992291136182789328035130030E+02;

            static constexpr double d71 = -0.25693933462703749003312586129E+02;
            static constexpr double d76 = -0.15418974869023643374053993627E+03;
            static constexpr double d77 = -0.23152937917604549567536039109E+03;
            static constexpr double d78 =  0.35763911791061412378285349910E+03;
            static constexpr double d79 =  0.93405324183624310003907691704E+02;
            static constexpr double d7A = -0.37458323136451633156875139351E+02;
            static constexpr double d7B =  0.10409964950896230045147246184E+03;
            static constexpr double d7C =  0.29840293426660503123344363579E+02;
            static constexpr double d7D = -0.43533456590011143754432175058E+02;
            static constexpr double d7E =  0.96324553959188282948394950600E+02;
            static constexpr double d7F = -0.39177261675615439165231486172E+02;
            static constexpr double d7G = -0.14972683625798562581422125276E+03;

            if (!icomp) nrd = size;

            for (std::size_t j = 0; j < nrd; j++) {
                std::size_t const i = icomp? icomp[j]: j; // index Of the variable to output
                cont[j] = value[i];
                double const ydiff = x[i] - value[i];
                cont[j + nrd] = ydiff;
                double const bspl = h * k1[i] - ydiff;
                cont[j + nrd * 2] = bspl;
                cont[j + nrd * 3] = ydiff - h * kD[i] - bspl;
                cont[j + nrd * 4] = d41 * k1[i] + d46 * k6[i] + d47 * k7[i] + d48 * k8[i] + d49 * k9[i] + d4A * kA[i] + d4B * kB[i] + d4C * kC[i];
                cont[j + nrd * 5] = d51 * k1[i] + d56 * k6[i] + d57 * k7[i] + d58 * k8[i] + d59 * k9[i] + d5A * kA[i] + d5B * kB[i] + d5C * kC[i];
                cont[j + nrd * 6] = d61 * k1[i] + d66 * k6[i] + d67 * k7[i] + d68 * k8[i] + d69 * k9[i] + d6A * kA[i] + d6B * kB[i] + d6C * kC[i];
                cont[j + nrd * 7] = d71 * k1[i] + d76 * k6[i] + d77 * k7[i] + d78 * k8[i] + d79 * k9[i] + d7A * kA[i] + d7B * kB[i] + d7C * kC[i];
            }

            for (std::size_t i = 0; i < size; i++) {
                double const _xE = value[i] + h * (aE1 * k1[i] + aE7 * k7[i] + aE8 * k8[i] + aE9 * k9[i] + aEA * kA[i] + aEB * kB[i] + aEC * kC[i] + aED * kD[i]);
                double const _xF = value[i] + h * (aF1 * k1[i] + aF6 * k6[i] + aF7 * k7[i] + aF8 * k8[i] + aFB * kB[i] + aFC * kC[i] + aFD * kD[i]);
                double const _xG = value[i] + h * (aG1 * k1[i] + aG6 * k6[i] + aG7 * k7[i] + aG8 * k8[i] + aG9 * k9[i] + aGD * kD[i]);
                /* xE == k6 */ k6[i] = _xE;
                /* xF == k7 */ k7[i] = _xF;
                /* xG == k8 */ k8[i] = _xG;
            }

            eq.eval_derivative(kE, size, time + cE * h, xE, rhs_pars);

            for (std::size_t i = 0; i < size; i++) xF[i] += h * aFE * kE[i];
            eq.eval_derivative(kF, size, time + cF * h, xF, rhs_pars);

            for (std::size_t i = 0; i < size; i++) xG[i] += h * (aGE * kE[i] + aGF * kF[i]);
            eq.eval_derivative(kG, size, time + cG * h, xG, rhs_pars);

            stat.nfcn += 3;
            for (std::size_t j = 0; j < nrd; j++) {
                std::size_t const i = icomp ? icomp[j]: j;
                cont[j + nrd * 4] = h * (cont[j + nrd * 4] + d4D * kD[i] + d4E * kE[i] + d4F * kF[i] + d4G * kG[i]);
                cont[j + nrd * 5] = h * (cont[j + nrd * 5] + d5D * kD[i] + d5E * kE[i] + d5F * kF[i] + d5G * kG[i]);
                cont[j + nrd * 6] = h * (cont[j + nrd * 6] + d6D * kD[i] + d6E * kE[i] + d6F * kF[i] + d6G * kG[i]);
                cont[j + nrd * 7] = h * (cont[j + nrd * 7] + d7D * kD[i] + d7E * kE[i] + d7F * kF[i] + d7G * kG[i]);
            }
        };


        double _determine_initial_step(
                double time, double* _restrict value, std::size_t size,
                iequation_for_erk& eq, void * rhs_pars,
                int bwd, double atol, double rtol, double hmax
        ) const {
            double* _restrict const x  = buffer.ptr<double>();
            double* _restrict const k1 = buffer.ptr<double>() + size * 1;
            double* _restrict const k2 = buffer.ptr<double>() + size * 2;

            // evaluateShycnhrotronSpectrum a first guess for explicit euler as:
            //   h1 = 0.01 * NORM(value) / NORM(k1).
            double dnf = 0.0, dny = 0.0;
            for (std::size_t i = 0; i < size; i++) {
                double const sk = atol + rtol * std::abs(value[i]);
                dnf += (k1[i] / sk) * (k1[i] / sk);
                dny += (value[i] / sk) * (value[i] / sk);
            }

            double h1;
            if (dnf <= 1e-10 || dny <= 1e-10)
                h1 = 1.0e-6;
            else
                h1 = std::sqrt(dny / dnf) * 0.01;
            h1 = std::min(h1, hmax);

            for (std::size_t i = 0; i < size; i++)
                x[i] = value[i] + bwd * h1 * k1[i];
            eq.eval_derivative(k2, size, time + bwd * h1, x, rhs_pars);

            // Second derivative norm | f'|
            double norm2 = 0.0;
            for (std::size_t i = 0; i < size; i++) {
                double const sk = atol + rtol * abs(value[i]);
                double const z =(k2[i] - k1[i]) / sk;
                norm2 += z * z;
            }
            norm2 = std::sqrt(norm2) / h1;

            // step m_size `h2' is computed such that:
            //   h2^order * max {norm(f), norm(f')} = 0.01
            double const der12 = std::max(norm2, std::sqrt(dnf));
            double h2;
            if (der12 <= 1e-15)
                h2 = std::max(1e-6, std::abs(h1) * 1e-3);
            else
                h2 = std::pow(0.01 / der12, 1.0 / order);

            return bwd * std::min(std::min(100 * h1, h2), hmax);
        };


    private:
        template<typename F>
        struct _equation_from_f : iequation_for_erk {
            F const& f;

            _equation_from_f(F const& f): f(f) { }

            virtual void eval_derivative(double* _restrict slope,
                                         std::size_t size,
                                         double t,
                                         double const* _restrict value,
                                         void * rhs_pars) override {
                f(slope, size, t, value, rhs_pars);
            }
        };

        template<typename F, typename CB>
        struct _equation_from_f_and_callback : iequation_for_erk {
            F const& f;
            CB const& stepCallback;

            _equation_from_f_and_callback(F const& f, CB const& stepCallback): f(f), stepCallback(stepCallback) {}
            virtual void eval_derivative(double* _restrict slope,
                                         std::size_t size,
                                         double t,
                                         double const* _restrict value,
                                         void * rhs_pars) override {
                f(slope, size, t, value, rhs_pars);
            }
            virtual void onstep() {
                stepCallback();
            }
        };


    public:
        /* --- main way to call the integrator --- */
        template<typename F>
        void operator()(double& time,
                        double* _restrict value,
                        std::size_t size,
                        F const& f,
                        void * rhs_pars,
                        double h) const {
            _equation_from_f<F> eq(f);

            buffer.ensure<double>(10 * size);

            double const atol = 1e-12;
            double const rtol = 1e-12;
            double err, stf;

            double* _restrict  x  = buffer.ptr<double>();
            double* _restrict  k1 = buffer.ptr<double>() + size * 1;
            eq.eval_derivative(k1, size, time, value, rhs_pars);
            this->_integrate8(
                    time, value, size, eq, rhs_pars, h,
                    atol, rtol, err, stf
            );

            for (std::size_t i = 0; i < size; i++)
                value[i] = x[i];
            time += h;
        }

        template<typename F>
        typename std::enable_if<!std::is_base_of<iequation_for_erk, F>::value, void>::type
        integrate(
                double& time, double* _restrict value, std::size_t size, F const& f,
                double timeN, stat_t& stat, param_t const& params
        ) const {
            _equation_from_f<F> eq(f);
            this->integrate(time, value, size, eq, timeN, stat, params);
        }

        template<typename F, typename CB>
        void integrate(
                double& time, double* _restrict value, std::size_t size, F const& f,
                double timeN, stat_t& stat, param_t const& params, CB const& stepCallback
        ) const {
            _equation_from_f_and_callback<F, CB> eq(f, stepCallback);
            this->integrate(time, value, size, eq, timeN, stat, params);
        }

        /* ----------------------------------- */
    public:
        /* --- interface method --- */
        void integrate(
                double& time,
                double* _restrict value,
                std::size_t size,
                iequation_for_erk& eq,
                double timeN,
                stat_t& stat,
                param_t const& params
        ) const {

            buffer.ensure<double>(10 * size);
            double const beta  = clamp(params.beta, 0.0, 0.2); // e.g. 0.04
            double const safe  = params.safe == 0.0? 0.9: clamp(params.safe, 1e-4, 1.0);
            double const facc1 = 1.0 / (params.fac1 == 0.0 ? 0.333 : params.fac1);
            double const facc2 = 1.0 / (params.fac2 == 0.0 ? 6.000 : params.fac2);
            double const expo1 = 1.0 / 8.0 - beta * 0.2;
            double const hmax  = std::abs(params.hmax == 0.0 ? timeN - time : params.hmax);
            int    const bwd   = time < timeN ? 1 : -1;
            double const rtol  = std::abs(params.rtol);
            double const atol  = std::abs(params.atol);
            int    const nstif = std::abs(params.nstif == 0 ? 10 : params.nstif);
            std::ptrdiff_t const nmax = params.nmax;
            void *    rhs_pars = params.rhs_pars;

            double* _restrict const x  = buffer.ptr<double>();
            double* _restrict const k1 = buffer.ptr<double>() + size * 1; // Where to put the first derivative (c = 0.0)
            double* _restrict const kD = buffer.ptr<double>() + size * 2; // Where to put the first derivative (FSAL) of the next step
            double* _restrict const kC = buffer.ptr<double>() + size * 3; // Where the final differential evaluation (c = 1.0) enters

            // first RHS evaluation
            eq.eval_derivative(k1, size, time, value, rhs_pars);
            stat.nfcn++;

            // evaluate what can be the initial step
            double h = bwd * std::abs(params.step);
            if (h == 0.0) {
                h = this->_determine_initial_step(time, value, size, eq, rhs_pars, bwd, atol, rtol, hmax);
                stat.nfcn++;
            }

            // initialize the struct for the dest data
            dense_output denseData(*this, stat, rhs_pars);

            double facold = 1e-4;
            bool reject = false, last = false;
            double hlamb = 0.0;
            int iasti = 0, nonsti = 0;
            for (/* none */ ; /* none */ ; stat.nstep++) {
//                mwg_check(nmax < 0 || stat.nstep < nmax, "Does not converge time = %g, h = %g at step#%d", time, h, stat.nstep);
                if ((nmax < 0 || stat.nstep < nmax) && params.verbocity > 0)
                    printf("\tDoes not converge time = %g, h = %g at step #%d \n", time, h, stat.nstep);
//                mwg_check(0.1 * std::abs(h) > std::abs(time) * std::numeric_limits<double>::epsilon(), "Time digit drop: time = %g, h = %g", time, h);
//                if (nmax < 0 || stat.nstep < nmax)
//                    std::cout << "\tDoes not converge x = "<<x<<" h = " << h << " d = " << stat.nstep << "\n";
                if ((0.1 * std::abs(h) > std::abs(time) * std::numeric_limits<double>::epsilon()) && params.verbocity > 0)
                    printf("\t\tTime digit drop: time = %g, h = %g \n", time, h);

                if ((time + 1.01 * h - timeN) * bwd > 0.0) {
                    h = timeN - time;
                    last = true;
                }

                //mwg_printd("time = %g, h = %g", time, h);
                double err, stf;
                this->_integrate8(
                        time, value, size, eq, rhs_pars, h, atol, rtol, err, stf
                );
                stat.nfcn += 11;

                double const fac11 = std::pow(err, expo1);
                double const fac = clamp(fac11 / (std::pow(facold, beta) * safe), facc2, facc1);
                double hnew = h / fac;
                if (err > 1.0) {
                    h /= std::min(facc1, fac11 / safe);
                    reject = true;
                    last = false;
                    if (stat.naccpt >= 1) stat.nrejct++;
                } else {
                    stat.naccpt++;
                    facold = std::max(err, 1e-4);
                    eq.eval_derivative(kD, size, time + h, x, rhs_pars);
                    stat.nfcn++;

                    // stiffness detection
                    if (stat.naccpt % nstif == 0 || iasti > 0) {
                        if (stf > 0.0) {
                            double stnum = 0.0;
                            for (std::size_t i = 0; i < size; i++)
                                stnum += (kD[i] - kC[i]) * (kD[i] - kC[i]);
                            hlamb = std::abs(h) * std::sqrt(stnum / stf);
                        }
                        if (hlamb > 6.1) {
                            nonsti = 0;
                            iasti++;
                            if (iasti != 15)
                                std::fprintf(stderr,"the problem seems to become stiff at time = %g\n", time);
                        } else {
                            nonsti++;
                            if (nonsti == 6) iasti = 0;
                        }
                    }

                    // Dense output

                    stat.eq = &eq;
                    stat.previousValue = value;
                    stat.previousSize  = size;
                    stat.previousTime  = time;
                    stat.previousStep  = h;
                    eq.ondense(stat, denseData);

                    for (std::size_t i = 0; i < size; i++) {
                        k1[i] = kD[i];
                        value[i] = x[i]; // overriding the solution array :: initial data -> solution
                    }

                    time += h;

                    eq.onstep();

                    if (eq.is_stopping() || last) return;

                    if (std::abs(hnew) > hmax) hnew = bwd * hmax;
                    if (reject) hnew = bwd * std::min(std::abs(hnew), std::abs(h));
                    reject = false;

                    h = hnew;
                }
            }
        };

        struct dense_output : idense_output {
            dop853_integrator const& integ;
            stat_t& stat;
            int bufferId {-1};
            working_buffer buffer;
            void * m_rhs_pars;

        public:
            /* --- constructor ---- */
            dense_output(dop853_integrator const& integ, stat_t& stat, void * rhs_pars)
                    : integ(integ), stat(stat), m_rhs_pars(rhs_pars) { }

        private:
            void initialize_data(std::size_t const* icomp, std::size_t ncomp) {
                if (this->bufferId == stat.nstep) return;

                this->bufferId = stat.nstep;

                std::size_t const size = stat.previousSize;
                buffer.ensure<double>(size * 8);
                integ._dense_output_initialize(buffer, stat, icomp, ncomp, m_rhs_pars);
            }

        public:
            virtual void get_values_at(
                    double* _restrict interpolated, double time,
                    std::size_t const* icomp, std::size_t _ncomp
            ) override {
//                mwg_check(stat.previousTime - 1e-14 <= time && time <= stat.previousTime + stat.previousStep + 1e-14,
//                          "time out of range: time=%g in [%g, %g]?",
//                          time, stat.previousTime, stat.previousTime + stat.previousStep);
                if (stat.previousTime - 1e-14 <= time && time <= stat.previousTime + stat.previousStep + 1e-14){
                    printf("time out of range: time=%g in [%g, %g]?", time, stat.previousTime, stat.previousTime + stat.previousStep);
                }

                this->initialize_data(icomp, _ncomp);

                double const s = (time - stat.previousTime) / stat.previousStep;
                double const t = 1.0 - s;

                double const* data = buffer.ptr<double>();
                std::size_t const ncomp = icomp ? _ncomp : stat.previousSize;
                for (std::size_t i = 0; i < ncomp; i++) {
                    int const j = icomp ? icomp[i]: i;

                    double a = data[j + 7 * ncomp];
                    a = a * s + data[j + 6 * ncomp];
                    a = a * t + data[j + 5 * ncomp];
                    a = a * s + data[j + 4 * ncomp];
                    a = a * t + data[j + 3 * ncomp];
                    a = a * s + data[j + 2 * ncomp];
                    a = a * t + data[j + 1 * ncomp];
                    a = a * s + data[j + 0 * ncomp];
                    interpolated[i] = a;
                }
            }
        };
    };

    /* ============================================================== */
    enum METHODS { EULER, MIDPOINT, HEUN, RALSTON, RK4, DOP853, DOP853E };
    static Integrators::METHODS selectOpt (const std::string & method, int loglevel ){
        // REMOVING LOGGER
//        logger p_log ( std::cerr, loglevel, "Integrators" );

        std::vector<std::string> str_METHODS { "Euler", "Midpoint", "Heun", "Ralston", "RK4", "DOP853", "DOP853E"};

        if (method == "Euler")
            return Integrators::METHODS::EULER;
        if (method == "Midpoint")
            return Integrators::METHODS::MIDPOINT;
        if (method == "Heun")
            return Integrators::METHODS::HEUN;
        if (method == "Ralston")
            return Integrators::METHODS::RALSTON;
        if (method == "RK4")
            return Integrators::METHODS::RK4;
        if (method == "DOP853")
            return Integrators::METHODS::DOP853;
        if (method == "DOP853E")
            return Integrators::METHODS::DOP853E;

        // REMOVING LOGGER
        std::cerr << AT  << "ODE integrator " << method << " is not recognized. Possible options: \n[ ";
        for (const auto& opt : str_METHODS){
            std::cout << opt << " ";
        }
        std::cout << "]\n";
        exit(1);
    }
}

class IntegratorRk4 {
    Integrators::rk4_integrator m_method;
protected:
    typedef void (*rhs)( double * ,        // solution array Y_1
                         std::size_t ,              // N_EQ of the system
                         double ,                   // X
                         double const * ,  // initial data Y_0
                         void *            // parameters for RHS
    );
    rhs m_Func      = nullptr; // pointer to the func to be integrated
    double m_X      = 0.0;
    double*  m_Y_i   = nullptr;
    double*  m_Y_o   = nullptr;
    void*  m_rhs_pars = nullptr;
    size_t m_n_eq   = 0;
    int m_loglevel{};
public:
    IntegratorRk4(int loglevel ) {
        m_method = Integrators::rk4_integrator();

        // REMOVING LOGGER
//        p_log = std::make_unique<logger>(std::cerr, m_loglevel, "ODE IntegratorStatic");

    }
    ~IntegratorRk4(){ std::cout << "deleing IntegratorRk4" << "\n"; }
//    void SetIntegrator( Method method ) { m_method = method; }
    void SetRHS( rhs Func, size_t neq, void *  pars ){
        m_Func = Func; m_n_eq = neq; m_rhs_pars = pars;
    }
    void SetInitCond( double x_0, double *  Y_0 ){ m_X = x_0; m_Y_o = Y_0; }
    void GetY( double * out_Y ) { for (size_t i = 0; i < m_n_eq; i++) out_Y[i] = m_Y_o[i]; }
    void Integrate( const double &dx ) {
        m_method( m_X, m_Y_o, m_n_eq, m_Func, m_rhs_pars, dx );
    }
//    void Integrate( const double &x0, const double &x1, double * Y0 ){
//        m_method( x0, Y0, m_n_eq, m_Func, m_rhs_pars, x1-x0 );
//    }
    void SetIntegratorPars( const std::unordered_map<std::string, double> &pars ) { unused(pars); }
    void SetIntegratorPars( Integrators::dop853_integrator::param_t & pars ) { unused(pars); }

};

class IntegratorBase{

protected:
    typedef void (*rhs)( double * ,        // solution array Y_1
                         std::size_t ,     // N_EQ of the system
                         double ,          // X
                         double const * ,  // initial data Y_0
                         void *            // parameters for RHS
    );
    rhs m_Func      = nullptr; // pointer to the func to be integrated
    double m_X      = 0.0;
    double*  m_Y_i   = nullptr;
    double*  m_Y_o   = nullptr;
    void*  m_rhs_pars = nullptr;
    size_t m_n_eq   = 0;
    int m_loglevel{};
    // REMOVING LOGGER
//    std::unique_ptr<logger> p_log;
    Integrators::dop853_integrator::param_t m_params{};
    Integrators::dop853_integrator::param_t * p_params{};
public:
    IntegratorBase() { m_params = Integrators::dop853_integrator::param_t(); p_params = &m_params; };
    virtual ~IntegratorBase() = default;
    void SetRHS( rhs Func, size_t neq, void * pars ){
        m_Func = Func;
        m_n_eq = neq;
        m_rhs_pars = pars;
    }
    void SetInitCond( double x_0, double * Y_0 ){ m_X = x_0; m_Y_o = Y_0; }
    void GetY( double * out_Y ) { for (size_t i = 0; i < m_n_eq; i++) out_Y[i] = m_Y_o[i]; }
    void SetY( const double * out_Y ) {
        for (size_t i = 0; i < m_n_eq; ++i)
            m_Y_o[i] = out_Y[i];
    }
    virtual void Integrate( const double x0, const double x1, double * Y0 ) = 0;
    virtual void Integrate( const double dx ) = 0;
    inline Integrators::dop853_integrator::param_t * pPars() { return p_params; }
};

template<typename Method>
class IntegratorStatic : public IntegratorBase {
    Method m_method;
public:
    IntegratorStatic( const Method method, int loglevel ) {
        m_method = method,
                m_loglevel = loglevel;
        // REMOVING LOGGER
//        p_log = std::make_unique<logger>(std::cerr, m_loglevel, "ODE IntegratorStatic");
    }
    void Integrate( const double dx ) override {
        m_method( m_X, m_Y_o, m_n_eq, m_Func, m_rhs_pars, dx );
    }
    /// updates the m_x with x0 and m_Y_o with Y0 and integrates till x1
    void Integrate( const double x0, const double x1, double * Y0 ) override {
        /// udate x and the state vector with new values:
        m_X = x0;
        for (size_t i = 0; i < m_n_eq; i++) { m_Y_o[i] = Y0[i]; }
        m_method( m_X, m_Y_o, m_n_eq, m_Func, m_rhs_pars, x1-x0 );
    }
};

template<typename Method>
class IntegratorEmbedded : public IntegratorBase {
    Method m_method;
public:

    IntegratorEmbedded( const Method method, int loglevel ) {
        m_method = method;
        m_loglevel = loglevel;
    }
    void SetIntegratorPars(const double atol=1e-11, const double rtol=1e-11, const double beta=0.0,
                           const double fac1=0.0, const double fac2=0.0, const double safe=0.0,
                           const double hmax=0.0, const double step=0.0, const int nstif=0,
                           const long nmax=100000, const int verpocity = -1) {
        m_params.atol = atol;
        m_params.rtol = rtol;
        m_params.beta = beta;
        m_params.fac1 = fac1;
        m_params.fac2 = fac2;
        m_params.safe = safe;
        m_params.hmax = hmax;
        m_params.step = step;
        m_params.nstif = nstif;
        m_params.nmax = nmax;
        m_params.verbocity = verpocity;
    }
    /// updates the m_x with x0 and m_Y_o with Y0 and integrates till x1
    void Integrate( const double x0, const double x1, double * Y0 ) override {
        Integrators::stat_t stat;
        m_params.rhs_pars = m_rhs_pars;
        /// udate x and the state vector with new values:
        m_X = x0;
        for (size_t i = 0; i < m_n_eq; i++) { m_Y_o[i] = Y0[i]; }
        m_method.integrate(m_X,        // beginning of the step 'x'
                           m_Y_o,      // pointer to double[] with initial conidtions to be overridden?
                           m_n_eq,     // number of equations in the system
                           m_Func,     // function
                           x1,//m_X + dx,   // the end of the step 'x + dx'
                           stat,       // statistics of the solver and internal vars
                           m_params    // options for the solver i.e tolerence, max iterations, struct rhs_pars
        );
    }
    void Integrate( const double dx  )  override  {
        Integrators::stat_t stat;
        m_params.rhs_pars = m_rhs_pars;
        m_method.integrate(m_X,        // beginning of the step 'x'
                           m_Y_o,      // pointer to double[] with initial conidtions to be overridden?
                           m_n_eq,     // number of equations in the system
                           m_Func,     // function
                           m_X + dx,   // the end of the step 'x + dx'
                           stat,       // statistics of the solver and internal vars
                           m_params    // options for the solver i.e tolerence, max iterations, struct rhs_pars
        );
    }
};


#endif //SRC_ODE_SOLVERS_H
