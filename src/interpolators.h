//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_INTERPOLATORS_H
#define SRC_INTERPOLATORS_H

#include "logger.h"
#include "pch.h"
#include "utils.h"

static size_t findIndex( const double & x, const Array & arr, size_t N ) {
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N-2;

    unsigned int i = ((unsigned int) N) >> 1;
    unsigned int a = 0;
    unsigned int b = N-1;

    // https://stackoverflow.com/questions/4192440/is-there-any-difference-between-1u-and-1-in-c/4192469
    // untill the m_size of b-a > 1 continue shrinking the array, approaching the 'x'
    while (b-a > 1u) // 1U is an unsigned value with the single bit 0 set
    {
        i = (b+a) >> 1; // ???
        if (arr[i] > x)
            b = i;
        else
            a = i;
    }

    return (int)a;
}
static size_t findIndex( const double & x, const Vector & arr, size_t N ) {
    if(x <= arr[0])
        return 0;
    else if(x >= arr[N-1])
        return N-2;

    unsigned int i = ((unsigned int) N) >> 1;
    unsigned int a = 0;
    unsigned int b = N-1;

    // https://stackoverflow.com/questions/4192440/is-there-any-difference-between-1u-and-1-in-c/4192469
    // untill the m_size of b-a > 1 continue shrinking the array, approaching the 'x'
    while (b-a > 1u) // 1U is an unsigned value with the single bit 0 set
    {
        i = (b+a) >> 1; // ???
        if (arr[i] > x)
            b = i;
        else
            a = i;
    }

    return (int)a;
}
static inline double interpSegLin( size_t & a, size_t & b, const double & x, Array & X, Array & Y) {
    // take two indexis, 'a' 'b' and two arrays 'X' and 'Y' and interpolate
    // between them 'Y' for at 'x" of "X'
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];
    return ya + (yb-ya) * (x-xa)/(xb-xa);
}
static inline double interpSegLog( size_t & a, size_t & b, double x, Array & X, Array & Y) {
//        std::cout << a << ' ' << b << ' '<< x << ' '<< X << ' '<< Y << ' '<< N << "\n";
    double xa = X[a];
    double xb = X[b];
    double ya = Y[a];
    double yb = Y[b];

    return ya * std::pow(yb/ya, log(x/xa)/log(xb/xa));
}


class InterpBase{
public:
    enum METHODS {
        iLagrangeLinear,  // Lagrange linear interpolation
        iLagrangeUnivariate3m, // Lagrange univariate three - point linear interpolation
        iSmoothCubicPolynomial,  // Smooth interpolation Cubic polynomial interpolation
    };
    static InterpBase::METHODS select( const std::string & method, std::unique_ptr<logger> & p_log ){

        // REMOVING LOGGER
//        logger log(std::cerr, loglevel, "InterpBase" );

        METHODS _opt;

        if (method == "linear")
            _opt = InterpBase::METHODS::iLagrangeLinear;
        else if (method == "univariate")
            _opt = InterpBase::METHODS::iLagrangeUnivariate3m;
        else if (method == "cubic")
            _opt = InterpBase::METHODS::iSmoothCubicPolynomial;
        else{
            // REMOVING LOGGER
            (*p_log)(LOG_ERR,AT)
                << " Interpolation option=" << method << " is not recognized. Possible options: "
                << "linear" << " , " << "univariate" << " , " << "cubic" << "\n";
            exit(1);
        }
        return _opt;
    }
};
/**
 *
 * 1D and 2D interpolation method_shock_vel
 * ispired by https://github.com/lanqianlong/Interp2CPP
 * Methods similar to 'matlab interp2'
 *
 */
class Interp1d : public InterpBase{
public:

    Interp1d(Array &x, Array &y) : m_X(x), m_Y(y)
    {
        m_nX = x.size();
        m_nY = y.size();

        if (m_nX != m_nY){
            std::cerr << "Error: arrays have different m_size. X[" << m_nX << "] Y" << m_nY << "]\n";
            std::cerr << AT << "\n";
            exit(1);
        }
    }

    double Interpolate(const double &ix, METHODS method) {

        switch (method) {

            case iLagrangeLinear:
                return lgr(ix);
                break;

            case iLagrangeUnivariate3m:
                return lg3(ix);
                break;

            case iSmoothCubicPolynomial:
                static double s[5];
                spl(-1, ix, s);
                return s[4];
                break;
        }
        std::cerr << "Unexpected error at " << AT << "\n";
        std::cerr << AT << "\n";
        exit(1);
    }

    //return multipoint value, point number is m
    Array Interpolate(Array &t, METHODS method) {
        Array result ( t.size() );
        size_t m = t.size();
//        if (t == NULL)
//            return;
        double tempVal = 0.0;;
        for (size_t k = 0; k < m; k++) {
            tempVal = Interpolate(t[k], method);
            result[k] = tempVal;
        }
        return std::move(result);
    }

    static int FindIdx(Array &x, const double &t, int n) {
        int kk, m, lc;
        /*if(t<x[0])
        return( 0);
        if(t>=x[n-1])
        return(n-1);*/

        if (t <= x[1])
            kk = 0;
        else if (t >= x[n - 1])
            kk = n - 2;
        else {
            kk = 1;
            m = n;
            while (((kk - m) != 1) && ((kk - m) != -1)) {
                lc = (kk + m) / 2;
                if (t < x[lc - 1])
                    m = lc;
                else
                    kk = lc;
            }
            kk = kk - 1;
        }
        return(kk);
    }

private:
    Array m_X, m_Y;
    size_t m_nX, m_nY;

protected:
    /// Lagrange linear interpolation
    double lgr(const double &t)
    {
        int i, j, k, m;
        double z, s;
        z = 0.0f;

        if (m_nX < 1)
            return(z);
        if (m_nX == 1) {
            z = m_Y[0];
            return(z);
        }

        if (m_nX == 2) {
            z = (m_Y[0] * (t - m_X[1]) - m_Y[1] * (t - m_X[0])) / (m_X[0] - m_X[1]);
            return(z);
        }
        if (m_nX > 1 && t > m_X[m_nX - 1]) {
            z = m_Y[m_nX - 1];
            return(z);
        }
        if (m_nX > 1 && t < m_X[0]) {
            z = m_Y[0];
            return(z);
        }


        i = 0;
        while ((m_X[i] < t) && (i < m_nX)) {
            i = i + 1;
        }


        k = i - 4;
        if (k < 0)
            k = 0;
        m = i + 3;
        if (m > m_nX - 1)
            m = (int)m_nX - 1;
        for (i = k; i <= m; i++) {
            s = 1.0;
            for (j = k; j <= m; j++) {
                if (j != i)
                    s = s * (t - m_X[j]) / (m_X[i] - m_X[j]);
            }
            z = z + s * m_Y[i];
        }
        return(z);
    }

    //Lagrange univariate three - point linear interpolation
    double lg3(const double &t)
    {
        int i, j, k, m;
        double z, s;
        z = 0.0;
        if (m_nX < 1) return (z);
        if (m_nX == 1) { z = m_Y[0]; return(z); }
        if (m_nX == 2)
        {
            z = (m_Y[0] * (t - m_X[1]) - m_Y[1] * (t - m_X[0])) / (m_X[0] - m_X[1]);
            return(z);
        }

        if (m_nX > 0 && t < m_X[0]) {
            z = m_Y[0];
            return(z);
        }

        if (m_nX > 0 && t > m_X[m_nX - 1]) {
            z = m_Y[m_nX - 1];
            return(z);
        }


        if (t <= m_X[1]) {
            k = 0;
            m = 2;
        }
        else if (t >= m_X[m_nX - 2]) {
            k = (int)m_nX - 3;
            m = (int)m_nX - 1;
        }
        else {
            k = 1;
            m = (int)m_nX;
            while ((m - k) != 1) {
                i = (k + m) / 2;
                if (t < m_X[i - 1])
                    m = i;
                else
                    k = i;
            }
            k = k - 1;
            m = m - 1;
            if (fabs(t - m_X[k]) < fabs(t - m_X[m]))
                k = k - 1;
            else
                m = m + 1;
        }

        z = 0.0;
        for (i = k; i <= m; i++) {
            s = 1.0;
            for (j = k; j <= m; j++) {
                if (j != i)
                    s = s * (t - m_X[j]) / (m_X[i] - m_X[j]);
            }
            z = z + s * m_Y[i];

        }
        return (z);
    }

    //Smooth interpolation Cubic polynomial interpolation
    void spl(int k, double t, double* s)
    {
        int kk, m, lc;
        double u[5], p, q;

        s[4] = 0.0;
        s[0] = 0.0;
        s[1] = 0.0;
        s[2] = 0.0;
        s[3] = 0.0;

        if (m_nX < 1)
            return;
        if (m_nX == 1) {
            s[0] = m_Y[0];
            s[4] = m_Y[0];
            return;
        }
        if (m_nX == 2) {
            s[0] = m_Y[0];
            s[1] = (m_Y[1] - m_Y[0]) / (m_X[1] - m_X[0]);
            if (k < 0)
                s[4] = (m_Y[0] * (t - m_X[1]) - m_Y[1] * (t - m_X[0])) / (m_X[0] - m_X[1]);
            return;
        }


        if (k < 0 && m_nX > 0 && t < m_X[0]) {
            s[4] = m_Y[0];
            return;
        }

        if (k < 0 && m_nX > 0 && t > m_X[m_nX - 1]) {
            s[4] = m_Y[m_nX - 1];
            return;
        }

        if (k < 0) {
            if (t <= m_X[1])
                kk = 0;
            else if (t >= m_X[m_nX - 1])
                kk = (int)m_nX - 2;
            else {
                kk = 1;
                m = (int)m_nX;
                while (((kk - m) != 1) && ((kk - m) != -1)) {
                    lc = (kk + m) / 2;
                    if (t < m_X[lc - 1])
                        m = lc;
                    else
                        kk = lc;
                }
                kk = kk - 1;
            }
        }
        else
            kk = k;

        if (kk > m_nX - 1)
            kk = (int)m_nX - 2;


        u[2] = (m_Y[kk + 1] - m_Y[kk]) / (m_X[kk + 1] - m_X[kk]);
        if (m_nX == 3) {
            if (kk == 0) {
                u[3] = (m_Y[2] - m_Y[1]) / (m_X[2] - m_X[1]);
                u[4] = 2.0 * u[3] - u[2];
                u[1] = 2.0 * u[2] - u[3];
                u[0] = 2.0 * u[1] - u[2];
            } else {
                u[1] = (m_Y[1] - m_Y[0]) / (m_X[1] - m_X[0]);
                u[0] = 2.0 * u[1] - u[2];
                u[3] = 2.0 * u[2] - u[1];
                u[4] = 2.0 * u[3] - u[2];
            }


        } else {
            if (kk <= 1) {
                u[3] = (m_Y[kk + 2] - m_Y[kk + 1]) / (m_X[kk + 2] - m_X[kk + 1]);
                if (kk == 1) {
                    u[1] = (m_Y[1] - m_Y[0]) / (m_X[1] - m_X[0]);
                    u[0] = 2.0 * u[1] - u[2];
                    if (m_nX == 4)
                        u[4] = 2.0 * u[3] - u[2];
                    else
                        u[4] = (m_Y[4] - m_Y[3]) / (m_X[4] - m_X[3]);
                } else {
                    u[1] = 2.0 * u[2] - u[3];
                    u[0] = 2.0 * u[1] - u[2];
                    u[4] = (m_Y[3] - m_Y[2]) / (m_X[3] - m_X[2]);
                }


            } else if (kk >= (m_nX - 3)) {
                u[1] = (m_Y[kk] - m_Y[kk - 1]) / (m_X[kk] - m_X[kk - 1]);
                if (kk == (m_nX - 3)) {
                    u[3] = (m_Y[m_nX - 1] - m_Y[m_nX - 2]) / (m_X[m_nX - 1] - m_X[m_nX - 2]);
                    u[4] = 2.0 * u[3] - u[2];
                    if (m_nX == 4)
                        u[0] = 2.0 * u[1] - u[2];
                    else
                        u[0] = (m_Y[kk - 1] - m_Y[kk - 2]) / (m_X[kk - 1] - m_X[kk - 2]);
                } else {
                    u[3] = 2.0 * u[2] - u[1];
                    u[4] = 2.0 * u[3] - u[2];
                    u[0] = (m_Y[kk - 1] - m_Y[kk - 2]) / (m_X[kk - 1] - m_X[kk - 2]);
                }


            } else {
                u[1] = (m_Y[kk] - m_Y[kk - 1]) / (m_X[kk] - m_X[kk - 1]);
                u[0] = (m_Y[kk - 1] - m_Y[kk - 2]) / (m_X[kk - 1] - m_X[kk - 2]);
                u[3] = (m_Y[kk + 2] - m_Y[kk + 1]) / (m_X[kk + 2] - m_X[kk + 1]);
                u[4] = (m_Y[kk + 3] - m_Y[kk + 2]) / (m_X[kk + 3] - m_X[kk + 2]);
            }
        }


        s[0] = fabs(u[3] - u[2]);
        s[1] = fabs(u[0] - u[1]);
        if ((s[0] + 1.0 == 1.0) && (s[1] + 1.0 == 1.0))
            p = (u[1] + u[2]) / 2.0;
        else
            p = (s[0] * u[1] + s[1] * u[2]) / (s[0] + s[1]);
        s[0] = fabs(u[3] - u[4]);
        s[1] = fabs(u[2] - u[1]);
        if ((s[0] + 1.0 == 1.0) && (s[1] + 1.0 == 1.0))
            q = (u[2] + u[3]) / 2.0;
        else
            q = (s[0] * u[2] + s[1] * u[3]) / (s[0] + s[1]);


        s[0] = m_Y[kk];
        s[1] = p;
        s[3] = m_X[kk + 1] - m_X[kk];
        s[2] = (3.0 * u[2] - 2.0 * p - q) / s[3];
        s[3] = (p + q - 2.0 * u[2]) / (s[3] * s[3]);
        if (k < 0) {
            p = t - m_X[kk];
            s[4] = s[0] + s[1] * p + s[2] * p * p + s[3] * p * p * p;
        }
//        return;
    }




};
/**
 * Simple 2 dimansional interpolation.
 * Inspired by https://github.com/lanqianlong/Interp2CPP
 */
class Interp2d : public InterpBase{

public:

    /**
     *
     * Note: 'z' is a 1-D array filled via
     * k = 0
     * for (size_t i = 0; i < x.m_size(); i++){
     *      for (size_t j = 0; j < y.m_size(); j++){
     *          z[k] = some_value;
     *          k++;
     *      }
     * }
     *
     * Function used to access the element is
     *      idx = ix + x.m_size() * iy
     *
     *
     * @param x array of the m_size 'm'
     * @param y array of rhe m_size 'n'
     * @param z array of the m_size 'm*n'
     *
     */
    Interp2d(Array &x, Array &y, Array &z)
            : m_X(x), m_Y(y), m_Z(z)
    {
        m_nX = x.size();
        m_nY = y.size();
        m_nZ = z.size();

//        if (m_nX != m_nY){
//            std::cout << "Error: arrays have different m_size. X[" << m_nX << "] Y" << m_nY << "]\n";
//            exit(1);
//        }

        if (m_nZ != (m_nX * m_nY)){
            std::cerr << "Error: Z.m_size() != X.m_size()*Y.m_size() "
                         "Z[" << m_nZ << "] X[" << m_nX << "] Y[" << m_nY << "]";
            std::cerr << AT << "\n";
            exit(1);
        }

    }

    double Interpolate(double a, double b, Interp1d::METHODS method){
        //find a，b
        int tempi, tempj;
        double w1, w2, w;
        //double* tempx(2), tempz(2);
//        double tempx[2] = { 0,0 };
//        double tempz[2] = { 0,0 };

        Array tempx { 0,0 };
        Array tempz { 0,0 };

        tempi = Interp1d::FindIdx(m_X, a, (int)m_nX);
        tempj = Interp1d::FindIdx(m_Y, b, (int)m_nY);
        if (tempi < 0)
            tempi = 0;
        if (tempi > m_nX - 2)
            tempi = (int)m_nX - 2;
        if (tempj < 0)
            tempj = 0;
        if (tempj > m_nY - 2)
            tempj = (int)m_nY - 2;


        tempx[0] = m_X[tempi];
        tempx[1] = m_X[tempi + 1];
        /*tempz[0] = m_Z[tempj * m_nX + tempi];
        tempz[1] = m_Z[tempj * m_nX + tempi + 1];*/
//        tempz[0] = m_Z[tempj][tempi];
//        tempz[1] = m_Z[tempj][tempi + 1];
//        tempz[0] = m_Z[Idx(tempj, tempi)];
//        tempz[1] = m_Z[Idx(tempj,tempi + 1)];
        tempz[0] = m_Z[Idx(tempi,tempj)];
        tempz[1] = m_Z[Idx(tempi + 1, tempj)];

        //tempz[0] = m_Z[tempi][tempj];
        //tempz[1] = m_Z[tempi + 1][tempj];

        Interp1d o1(tempx, tempz);
        w1 = o1.Interpolate(a, method);
        //update tempz
        /*tempz[0] = m_Z[(tempj + 1) * m_nX + tempi];
        tempz[1] = m_Z[(tempj + 1) * m_nX + tempi + 1];*/
//        tempz[0] = m_Z[(tempj + 1)][tempi];
//        tempz[1] = m_Z[(tempj + 1)][tempi + 1];
//        tempz[0] = m_Z[Idx(tempj + 1, tempi)];
//        tempz[1] = m_Z[Idx(tempj + 1, tempi + 1)];
        tempz[0] = m_Z[Idx(tempi, tempj + 1)];
        tempz[1] = m_Z[Idx(tempi + 1, tempj + 1)];

        //tempz[0] = m_Z[tempi][tempj+1];
        //tempz[1] = m_Z[tempi+1][tempj + 1];
        Interp1d o2(tempx, tempz);
        w2 = o2.Interpolate(a, method);


        //update tempx and tempz
        tempx[0] = m_Y[tempj];
        tempx[1] = m_Y[tempj + 1];


        tempz[0] = w1;
        tempz[1] = w2;

        Interp1d o3(tempx, tempz);
        w = o3.Interpolate(b, method);
        //}
        return(w);
    }

    Array Interpolate(Array &x, Array &y, Interp1d::METHODS method){
        double tempVal = 0.0;
        Array result ( x.size() * y.size() );
        size_t idx = 0;
        for (size_t i = 0; i < x.size(); i++){
            for(size_t j = 0; j < y.size(); j++){
                result[idx] = Interpolate(x[i], y[j], method);
                idx++;
            }
        }
        return result;
    }


    double InterpolateBilinear(const double &x, const double &y,
                               const int & ix1, const int &ix2, const int & iy1, const int & iy2){
        double z = BilinearInterpolation(
                m_Z[Idx(ix1, iy1)],
                m_Z[Idx(ix1, iy2)],
                m_Z[Idx(ix2, iy1)],
                m_Z[Idx(ix2, iy2)],
                m_X[ix1], m_X[ix2], m_Y[iy1], m_Y[iy2], x, y
        );
        return z;
    }

    double InterpolateBilinear(const double &x, const double &y){
        double z = 0.0;

        // check if in-bounds
        if (x > m_X[m_nX-1]){
            std::cerr << "X=" << x << " is above max(X)=" << m_X[m_nX-1] << std::endl;
            std::cerr << AT << "\n";
            exit(1);
        }
        if (x < m_X[0]){
            std::cerr << "X=" << x << " is below min(X)=" << m_X[0] << std::endl;
            std::cerr << AT << "\n";
            exit(1);
        }
        if (y > m_Y[m_nY-1]){
            std::cerr << "Y=" << y << " is above max(Y)=" << m_Y[m_nY-1] << std::endl;
            std::cerr << AT << "\n";
            exit(1);
        }
        if (y < m_Y[0]){
            std::cerr << "Y=" << y << " is below min(Y)=" << m_Y[0] << std::endl;
            std::cerr << AT << "\n";
            exit(1);
        }

        // find the i_X before the 'x' point
        int ix1 = 0;
        for (int i = 0; i < m_nX; i++){
            if (x > m_X[i])
                ix1 = i;
            else
                break;
        }
        int ix2 = ix1+1;

        int iy1 = 0;
        for (int i = 0; i < m_nY; i++){
            if (y > m_Y[i])
                iy1 = i;
            else
                break;
        }
        int iy2 = iy1+1;

        z = BilinearInterpolation(
                m_Z[Idx(ix1, iy1)],
                m_Z[Idx(ix1, iy2)],
                m_Z[Idx(ix2, iy1)],
                m_Z[Idx(ix2, iy2)],
                m_X[ix1], m_X[ix2], m_Y[iy1], m_Y[iy2], x, y
        );
        return z;
    }

private:
    Array m_X, m_Y, m_Z;
    size_t m_nX, m_nY, m_nZ;

    inline size_t Idx( int x, int y ) const { return x + m_nX * y; } // m_width

    static inline double BilinearInterpolation(double &q11, double &q12, double &q21, double &q22, double &x1,
                                               double &x2, double &y1, double &y2, const double &x, const double &y) {
        double x2x1, y2y1, x2x, y2y, yy1, xx1;
        x2x1 = x2 - x1;
        y2y1 = y2 - y1;
        x2x = x2 - x;
        y2y = y2 - y;
        yy1 = y - y1;
        xx1 = x - x1;
        return 1.0 / (x2x1 * y2y1) * (
                q11 * x2x * y2y +
                q21 * xx1 * y2y +
                q12 * x2x * yy1 +
                q22 * xx1 * yy1
        );
    }

};

class TriliniarInterpolation{
    Vector & x_arr;
    Vector & y_arr;
    Vector & z_arr;
    Vector & f_arr;
public:
    TriliniarInterpolation(Vector & x_arr, Vector & y_arr, Vector & z_arr, Vector & f_arr)
        : x_arr(x_arr), y_arr(y_arr), z_arr(z_arr), f_arr(f_arr){
        ///
    }
    double interp(double x_val, double y_val, double z_val){
        if (x_val < x_arr[0])
            x_val = x_arr[0];
        if (x_val > x_arr[x_arr.size()-1])
            x_val = x_arr[x_arr.size()-1];

        if (y_val < y_arr[0])
            y_val = y_arr[0];
        if (y_val > y_arr[y_arr.size()-1])
            y_val = y_arr[y_arr.size()-1];

        if (z_val < z_arr[0])
            z_val = z_arr[0];
        if (z_val > z_arr[z_arr.size()-1])
            z_val = z_arr[z_arr.size()-1];

        size_t _i_ye = findIndex(x_val, x_arr, x_arr.size());
        size_t _i_s = findIndex(y_val, y_arr, y_arr.size());
        size_t _i_tau = findIndex(z_val, z_arr, z_arr.size());

        size_t _ip1_ye = _i_ye + 1;
        if (_i_ye == x_arr.size()-1)
            _ip1_ye = _i_ye;

        size_t _ip1_s = _i_s + 1;
        if (_i_s == y_arr.size()-1)
            _ip1_s = _i_s;

        size_t _ip1_tau = _i_tau + 1;
        if (_i_tau == z_arr.size()-1)
            _ip1_tau = _i_tau;

        // _idx = k + len(new_s) * (j + len(new_ye) * i)
        double _i_val = f_arr[ _i_tau + y_arr.size() * (_i_s + x_arr.size() * _i_ye) ];
        double _ip1_val = f_arr[ _ip1_tau + y_arr.size() * (_ip1_s + z_arr.size() * _ip1_ye) ];

        double val = trilinear_interpolation( x_arr[_i_ye], y_arr[_i_s], z_arr[_i_tau], _i_val,
                                              x_arr[_ip1_ye], y_arr[_ip1_s], z_arr[_ip1_tau], _ip1_val,
                                              x_val, y_val, z_val);
        // print("my={:.4e} scipy={:.4e}".format(val,val1))
        // if (val < min(_i_val,_ip1_val) or val > max(_i_val, _ip1_val)):
        //     raise ValueError(f"val={val} < min(_i_val,_ip1_val)={min(_i_val,_ip1_val)} "
        //                 f"or val > max(_i_val, _ip1_val)={max(_i_val, _ip1_val)}")
        // return _ip1_val
        return val;

    }
private:
    static double trilinear_interpolation(double x1, double y1, double z1, double f1,
                                          double x2, double y2, double z2, double f2,
                                          double x3, double y3, double z3) {
        // Calculate interpolation factors
        double d = (x3 - x1) / (x2 - x1);
        double e = (y3 - y1) / (y2 - y1);
        double f = (z3 - z1) / (z2 - z1);

        // Get the values of f at the eight corners of the cube
        double f3 = f1 + d * (f2 - f1);
        double f4 = f1 + d * (f2 - f1) + e * (f2 - f1 - (f3 - f1));
        double f5 = f1 + d * (f2 - f1) + (1 - e) * (f3 - f1);
        double f6 = f1 + d * (f2 - f1) + e * (f2 - f1 - (f3 - f1)) + (1 - f) * (f3 - f1);
        double f7 = f1 + (1 - d) * (f2 - f1) + e * (f2 - f1 - (f3 - f1));
        double f8 = f1 + (1 - d) * (f2 - f1) + e * (f2 - f1 - (f3 - f1)) + (1 - f) * (f3 - f1);

        // Interpolate along the fourth dimension
        double res = (1 - d) * ((1 - e) * ((1 - f) * f1 + f * f2) + e * ((1 - f) * f3 + f * f4))
                     + d * ((1 - e) * ((1 - f) * f5 + f * f6) + e * ((1 - f) * f7 + f * f8));

        return res;
    }
};

inline double lagrangeInterpolation(Array & x, Array & y, size_t n, double xp){
    double intp = 0;
    for (size_t i = 0; i < n; i++){
        double m = 1;
        for (size_t j = 0; j < n; j++){
            if (i != j)
                m = m * (xp - x[j]) / (x[i] - x[j]);
        }
        m = m * y[i];
        intp = intp + m;
    }
    return intp;
}

inline double dydx(Array & x, Array & y, Array & dy, double dx, size_t idx, double curr_x, bool lagrange_interp=true){
    double dy_idx;
    if (idx == 2){
        // first order derivative
        dy_idx = (y[idx] - y[idx-1]) / dx;
    }
    if (idx == 3){
        // second order derivative
        dy_idx = (3.*y[idx] - 4.*y[idx-1] + y[idx-2]) / (2. * dx); //
    }
    if (idx == 4){
        // third order derivative
        dy_idx = (10. * y[idx] - 15.*y[idx-1] + 6.*y[idx-2] - y[idx-3]) / (6. * dx); // error O(dx^2)
    }
    if (idx >= 5){
        dy_idx = (35. * y[idx] - 56. * y[idx-1] + 28. * y[idx-2] - 8. * y[idx-3] + y[idx-4]) / (20. * dx);
    }

    /// interpolate solution if not at the grid point
    if (curr_x == x[idx]) return dy_idx;
    if ((curr_x < x[idx])&&(!lagrange_interp)){
        // two point lagrange interpolation
        return dy[idx-1] + ( dy[idx] - dy[idx-1] ) / ( x[idx] - x[idx-1] ) * ( curr_x - x[idx-1] );
    }
    if ((curr_x < x[idx])&&(lagrange_interp)){
        // full lagrange interpolation
        return lagrangeInterpolation(x, dy, idx+1, curr_x);
    }
    std::cerr << AT << " unexpected error\n";
    exit(1);
}


#endif //SRC_INTERPOLATORS_H
