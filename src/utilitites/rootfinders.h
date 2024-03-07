//
// Created by vsevolod on 21/12/22.
//

#ifndef SRC_ROOTFINDERS_H
#define SRC_ROOTFINDERS_H

#include "pch.h"

/*
 * Basic bisection routine to find a zero of the function `f` between the
 * arguments `a` and `b`. `f(a)` and `f(b)` cannot have the same signs.
 * Slow but sure.
 *
 * Source: https://github.com/scipy/scipy/blob/master/scipy/optimize/Zeros/bisect.c
 */
double Bisect(double (*f)(const double&, void*), double xa, double xb,
              double xtol, double rtol, int iter, void *func_data, int & status) {
    int i;
    double dm,xm,fm,fa,fb;
    status = 0;

    fa = f(xa, func_data); //fa = (*f)(xa, func_data);
    fb = f(xb, func_data); //(*f)(xb, func_data);

    if (fa*fb > 0) {
        status = -1;
//        std::cerr<< "[ Error ] rootfinder :: fa*fb=0\n";
        return 0.;
    }
    if (fa == 0) {
        status = -1;
//        std::cerr<< "[ Error ] rootfinder :: fa=0\n";
        return xa;
    }
    if (fb == 0) {
        status = -1;
//        std::cerr<< "[ Error ] rootfinder :: fb=0\n";
        return xb;
    }
    dm = xb - xa;
    for (i=0; i<iter; i++) {

        dm *= 0.5;
        xm = xa + dm;
        fm = f(xm, func_data); //(*f)(xm, func_data);

        if (fm*fa >= 0) {
            xa = xm;
        }
        if (fm == 0 || fabs(dm) < xtol + rtol*fabs(xm)) {
            return xm;
        }
    }
    status = -1;
    std::cerr << "[ Error ] rootfinder :: Root not found\n";
    return xa;
}


//****************************************************************************80
//  The following code is COPIED from open source:
//
//  https://people.sc.fsu.edu/~jburkardt/cpp_src/fsolve/fsolve.cpp
//****************************************************************************80

void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
              double delta, double x[], double wa1[], double wa2[] );
double enorm ( int n, double x[] );
void fdjac1 ( void fcn ( int n, double x[], double f[], int &iflag, void * pars ),
              int n, double x[], double fvec[], double fjac[], int ldfjac, int &iflag,
              int ml, int mu, double epsfcn, double wa1[], double wa2[] );
int fsolve ( void fcn ( int n, double x[], double fvec[], int &iflag, void * pars ), int n,
             double x[], void * pars, double fvec[], double tol, double wa[], int lwa );
int hybrd ( void fcn ( int n, double x[], double fvec[], int &iflag, void * pars ),
            int n, double x[], void * pars, double fvec[], double xtol, int maxfev, int ml,
            int mu, double epsfcn, double diag[], int mode, double factor, int nprint,
            int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[],
            double wa1[], double wa2[], double wa3[], double wa4[] );
void qform ( int m, int n, double q[], int ldq );
void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
             int lipvt, double rdiag[], double acnorm[] );
void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] );
bool r1updt ( int m, int n, double s[], int ls, double u[], double v[], double w[] );



//****************************************************************************80

void dogleg ( int n, double r[], int lr, double diag[], double qtb[],
              double delta, double x[], double wa1[], double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    dogleg() combines Gauss-Newton and gradient for a minimizing step.
//
//  Discussion:
//
//    Given an M by N matrix A, an n by n nonsingular diagonal
//    matrix d, an m-vector b, and a positive number delta, the
//    problem is to determine the convex combination x of the
//    gauss-newton and scaled gradient directions that minimizes
//    (a*x - b) in the least squares sense, subject to the
//    restriction that the euclidean norm of d*x be at most delta.
//
//    This function completes the solution of the problem
//    if it is provided with the necessary information from the
//    qr factorization of a.
//
//    That is, if a = q*r, where q has orthogonal columns and r is an upper
//    triangular matrix, then dogleg expects the full upper triangle of r and
//    the first n components of Q'*b.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int N, the order of R.
//
//    Input, double R[LR], the upper triangular matrix R stored by rows.
//
//    Input, int LR, the size of the storage for R, which should be at
//    least (n*(n+1))/2.
//
//    Input, double DIAG[N], the diagonal elements of the matrix D.
//
//    Input, double QTB[N], the first n elements of the vector
//    (q transpose)*b.
//
//    Input, double DELTA, an upper bound on the euclidean norm of d*x.
//
//    Output, double X[N], contains the desired convex combination of the
//    gauss-newton direction and the scaled gradient direction.
//
//    Workspace, WA1[N].
//
//    Workspace, WA2[N].
//
{
    double alpha;
    double bnorm;
    double epsmch;
    double gnorm;
    int i;
    int j;
    int jj;
    int jp1;
    int k;
    int l;
    double qnorm;
    double sgnorm;
    double sum;
    double temp;
//
//  EPSMCH is the machine precision.
//
    epsmch = std::numeric_limits<double>::epsilon();
//
//  Calculate the Gauss-Newton direction.
//
    jj = ( n * ( n + 1 ) ) / 2 + 1;

    for ( k = 1; k <= n; k++ )
    {
        j = n - k + 1;
        jp1 = j + 1;
        jj = jj - k;
        l = jj + 1;
        sum = 0.0;
        for ( i = jp1; i <= n; i++ )
        {
            sum = sum + r[l-1] * x[i-1];
            l = l + 1;
        }
        temp = r[jj-1];
        if ( temp == 0.0 )
        {
            l = j;
            for ( i = 1; i <= j; i++ )
            {
                temp = fmax ( temp, fabs ( r[l-1] ) );
                l = l + n - i;
            }
            temp = epsmch * temp;
            if ( temp == 0.0 )
            {
                temp = epsmch;
            }
        }
        x[j-1] = ( qtb[j-1] - sum ) / temp;
    }
//
//  Test whether the Gauss-Newton direction is acceptable.
//
    for ( j = 0; j < n; j++ )
    {
        wa1[j] = 0.0;
        wa2[j] = diag[j] * x[j];
    }
    qnorm = enorm ( n, wa2 );

    if ( qnorm <= delta )
    {
        return;
    }
//
//  The Gauss-Newton direction is not acceptable.
//  Calculate the scaled gradient direction.
//
    l = 0;
    for ( j = 0; j < n; j++ )
    {
        temp = qtb[j];
        for ( i = j; i < n; i++ )
        {
            wa1[i-1] = wa1[i-1] + r[l-1] * temp;
            l = l + 1;
        }
        wa1[j] = wa1[j] / diag[j];
    }
//
//  Calculate the norm of the scaled gradient and test for
//  the special case in which the scaled gradient is zero.
//
    gnorm = enorm ( n, wa1 );
    sgnorm = 0.0;
    alpha = delta / qnorm;
//
//  Calculate the point along the scaled gradient
//  at which the quadratic is minimized.
//
    if ( gnorm != 0.0 )
    {
        for ( j = 0; j < n; j++ )
        {
            wa1[j] = ( wa1[j] / gnorm ) / diag[j];
        }
        l = 0;
        for ( j = 0; j < n; j++ )
        {
            sum = 0.0;
            for ( i = j; i < n; i++ )
            {
                sum = sum + r[l] * wa1[i];
                l = l + 1;
            }
            wa2[j] = sum;
        }
        temp = enorm ( n, wa2 );
        sgnorm = ( gnorm / temp ) / temp;
        alpha = 0.0;
//
//  If the scaled gradient direction is not acceptable,
//  calculate the point along the dogleg at which the quadratic is minimized.
//
        if ( sgnorm < delta)
        {
            bnorm = enorm ( n, qtb );
            temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta );
            temp = temp - ( delta / qnorm ) * ( sgnorm / delta ) * ( sgnorm / delta )
                   + sqrt ( std::pow ( temp - ( delta / qnorm ), 2 )
                            + ( 1.0 - ( delta / qnorm ) * ( delta / qnorm ) )
                              * ( 1.0 - ( sgnorm / delta ) * ( sgnorm / delta ) ) );
            alpha = ( ( delta / qnorm )
                      * ( 1.0 - ( sgnorm / delta ) * ( sgnorm / delta ) ) ) / temp;
        }
    }
//
//  Form appropriate convex combination of the Gauss-Newton
//  direction and the scaled gradient direction.
//
    temp = ( 1.0 - alpha ) * fmin ( sgnorm, delta );
    for ( j = 0; j < n; j++ )
    {
        x[j] = temp * wa1[j] + alpha * x[j];
    }
    return;
}
//****************************************************************************80

double enorm ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    enorm() returns the Euclidean norm of a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    int N, the number of entries in A.
//
//    double X[N], the vector whose norm is desired.
//
//  Output:
//
//    double ENORM, the norm of X.
//
{
    int i;
    double value;

    value = 0.0;
    for ( i = 0; i < n; i++ )
    {
        value = value + x[i] * x[i];
    }
    value = sqrt ( value );
    return value;
}
//****************************************************************************80

void fdjac1 ( void fcn ( int n, double x[], double f[], int &iflag, void * pars ),
              int n, double x[], void * pars, double fvec[], double fjac[], int ldfjac, int &iflag,
              int ml, int mu, double epsfcn, double wa1[], double wa2[] )

//****************************************************************************80
//
//  Purpose:
//
//    fdjac1() estimates an N by N Jacobian matrix using forward differences.
//
//  Discussion:
//
//    This function computes a forward-difference approximation
//    to the N by N jacobian matrix associated with a specified
//    problem of N functions in N variables.
//
//    If the jacobian has a banded form, then function evaluations are saved
//    by only approximating the nonzero terms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//       fcn is the name of the user-supplied subroutine which
//         calculates the functions. fcn must be declared
//         in an external statement in the user calling
//         program, and should be written as follows.
//
//         subroutine fcn(n,x,fvec,iflag)
//         integer n,iflag
//         double precision x(n),fvec(n)
//         ----------
//         calculate the functions at x and
//         return this vector in fvec.
//         ----------
//         return
//         end
//
//         the value of iflag should not be changed by fcn unless
//         the user wants to terminate execution of fdjac1.
//         in this case set iflag to a negative integer.
//
//    Input, int N, the number of functions and variables.
//
//    Input, double X[N], the evaluation point.
//
//    Input, double FVEC[N], the functions evaluated at X.
//
//    Output, double FJAC[N*N], the approximate jacobian matrix at X.
//
//       ldfjac is a positive integer input variable not less than n
//         which specifies the leading dimension of the array fjac.
//
//       iflag is an integer variable which can be used to terminate
//         the execution of fdjac1. see description of fcn.
//
//       ml is a nonnegative integer input variable which specifies
//         the number of subdiagonals within the band of the
//         jacobian matrix. if the jacobian is not banded, set
//         ml to at least n - 1.
//
//       epsfcn is an input variable used in determining a suitable
//         step length for the forward-difference approximation. this
//         approximation assumes that the relative errors in the
//         functions are of the order of epsfcn. if epsfcn is less
//         than the machine precision, it is assumed that the relative
//         errors in the functions are of the order of the machine
//         precision.
//
//       mu is a nonnegative integer input variable which specifies
//         the number of superdiagonals within the band of the
//         jacobian matrix. if the jacobian is not banded, set
//         mu to at least n - 1.
//
//       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
//         least n, then the jacobian is considered dense, and wa2 is
//         not referenced.
{
    double eps;
    double epsmch;
    double h;
    int i;
    int j;
    int k;
    int msum;
    double temp;
//
//  EPSMCH is the machine precision.
//
    epsmch = std::numeric_limits<double>::epsilon();

    eps = sqrt ( fmax ( epsfcn, epsmch ) );
    msum = ml + mu + 1;
//
//  Computation of dense approximate jacobian.
//
    if ( n <= msum )
    {
        for ( j = 0; j < n; j++ )
        {
            temp = x[j];
            h = eps * fabs ( temp );
            if ( h == 0.0 )
            {
                h = eps;
            }
            x[j] = temp + h;
            fcn ( n, x, wa1, iflag, pars );
            if ( iflag < 0 )
            {
                break;
            }
            x[j] = temp;
            for ( i = 0; i < n; i++ )
            {
                fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
            }
        }
    }
//
//  Computation of a banded approximate jacobian.
//
    else
    {
        for ( k = 0; k < msum; k++ )
        {
            for ( j = k; j < n; j = j + msum )
            {
                wa2[j] = x[j];
                h = eps * fabs ( wa2[j] );
                if ( h == 0.0 )
                {
                    h = eps;
                }
                x[j] = wa2[j] + h;
            }
            fcn ( n, x, wa1, iflag, pars );
            if ( iflag < 0 )
            {
                break;
            }
            for ( j = k; j < n; j = j + msum )
            {
                x[j] = wa2[j];
                h = eps * fabs ( wa2[j] );
                if ( h == 0.0 )
                {
                    h = eps;
                }
                for ( i = 0; i < n; i++ )
                {
                    if ( j - mu <= i && i <= j + ml )
                    {
                        fjac[i+j*ldfjac] = ( wa1[i] - fvec[i] ) / h;
                    }
                    else
                    {
                        fjac[i+j*ldfjac] = 0.0;
                    }
                }
            }
        }
    }
    return;
}
//****************************************************************************80

int fsolve ( void fcn ( int n, double x[], double fvec[], int &iflag, void * pars ), int n,
             double x[], void * pars, double fvec[], double tol, double wa[], int lwa )

//****************************************************************************80
//
//  Purpose:
//
//    fsolve() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The purpose is to find a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.
//
//    This is done by using the more general nonlinear equation solver HYBRD.
//
//    The user must provide FCN, which calculates the functions.
//
//    The jacobian is calculated by a forward-difference approximation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//       fcn is the name of the user-supplied subroutine which
//         calculates the functions. fcn must be declared
//         in an external statement in the user calling
//         program, and should be written as follows.
//
//         subroutine fcn(n,x,fvec,iflag)
//         integer n,iflag
//         double precision x(n),fvec(n)
//         ----------
//         calculate the functions at x and
//         return this vector in fvec.
//         ---------
//         return
//         end
//
//         the value of iflag should not be changed by fcn unless
//         the user wants to terminate execution of hybrd1.
//         in this case set iflag to a negative integer.
//
//    Input, int N, the number of functions and variables.
//
//    Input/output, double X[N].  On input, an initial estimate of the solution.
//    On output, the final estimate of the solution.
//
//    Output, double FVEC[N], the functions evaluated at the output X.
//
//    Input, double TOL, a nonnegative variable. tTermination occurs when the
//    algorithm estimates that the relative error between X and the solution
//    is at most TOL.
//
//       info is an integer output variable. if the user has
//         terminated execution, info is set to the (negative)
//         value of iflag. see description of fcn. otherwise,
//         info is set as follows.
//         info = 0   improper input parameters.
//         info = 1   algorithm estimates that the relative error
//                    between x and the solution is at most tol.
//         info = 2   number of calls to fcn has reached or exceeded
//                    200*(n+1).
//         info = 3   tol is too small. no further improvement in
//                    the approximate solution x is possible.
//         info = 4   iteration is not making good progress.
//
//    Workspace, double WA[LWA].
//
//       lwa is a positive integer input variable not less than
//         (n*(3*n+13))/2.
//
{
    double epsfcn;
    double factor;
    int index;
    int info;
    int j;
    int lr;
    int maxfev;
    int ml;
    int mode;
    int mu;
    int nfev;
    int nprint;
    double xtol;

    info = 0;
//
//  Check the input.
//
    if ( n <= 0 )
    {
        return info;
    }
    if ( tol <= 0.0 )
    {
        return info;
    }
    if ( lwa < ( n * ( 3 * n + 13 ) ) / 2 )
    {
        return info;
    }
//
//  Call HYBRD.
//
    xtol = tol;
    maxfev = 200 * ( n + 1 );
    ml = n - 1;
    mu = n - 1;
    epsfcn = 0.0;
    for ( j = 0; j < n; j++ )
    {
        wa[j] = 1.0;
    }
    mode = 2;
    factor = 100.0;
    nprint = 0;
    nfev = 0;
    lr = ( n * ( n + 1 ) ) / 2;
    index = 6 * n + lr;

    info = hybrd ( fcn, n, x, pars, fvec, xtol, maxfev, ml, mu, epsfcn, wa, mode,
                   factor, nprint, nfev, wa+index, n, wa+6*n, lr,
                   wa+n, wa+2*n, wa+3*n, wa+4*n, wa+5*n );

    if ( info == 5 )
    {
        info = 4;
    }
    return info;
}
//****************************************************************************80

int hybrd ( void fcn ( int n, double x[], double fvec[], int &iflag, void * pars ),
            int n, double x[], void * pars,
            double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn,
            double diag[], int mode, double factor, int nprint, int nfev,
            double fjac[], int ldfjac, double r[], int lr, double qtf[], double wa1[],
            double wa2[], double wa3[], double wa4[] )

//****************************************************************************80
//
//  Purpose:
//
//    hybrd() finds a zero of a system of N nonlinear equations.
//
//  Discussion:
//
//    The purpose of HYBRD is to find a zero of a system of
//    N nonlinear functions in N variables by a modification
//    of the Powell hybrid method.
//
//    The user must provide FCN, which calculates the functions.
//
//    The jacobian is calculated by a forward-difference approximation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//       fcn is the name of the user-supplied subroutine which
//         calculates the functions. fcn must be declared
//         in an external statement in the user calling
//         program, and should be written as follows.
//
//         subroutine fcn(n,x,fvec,iflag)
//         integer n,iflag
//         double precision x(n),fvec(n)
//         ----------
//         calculate the functions at x and
//         return this vector in fvec.
//         ---------
//         return
//         end
//
//         the value of iflag should not be changed by fcn unless
//         the user wants to terminate execution of hybrd.
//         in this case set iflag to a negative integer.
//
//    Input, int N, the number of functions and variables.
//
//    Input/output, double X[N].  On input an initial estimate of the solution.
//    On output, the final estimate of the solution.
//
//    Output, double FVEC[N], the functions evaluated at the output value of X.
//
//    Input, double XTOL, a nonnegative value.  Termination occurs when the
//    relative error between two consecutive iterates is at most XTOL.
//
//    Input, int MAXFEV.  Termination occurs when the number of calls to FCN
//    is at least MAXFEV by the end of an iteration.
//
//    Input, int ML, specifies the number of subdiagonals within the band of
//    the jacobian matrix.  If the jacobian is not banded, set
//    ml to at least n - 1.
//
//       mu is a nonnegative integer input variable which specifies
//         the number of superdiagonals within the band of the
//         jacobian matrix. if the jacobian is not banded, set
//         mu to at least n - 1.
//
//       epsfcn is an input variable used in determining a suitable
//         step length for the forward-difference approximation. this
//         approximation assumes that the relative errors in the
//         functions are of the order of epsfcn. if epsfcn is less
//         than the machine precision, it is assumed that the relative
//         errors in the functions are of the order of the machine
//         precision.
//
//       diag is an array of length n. if mode = 1 (see
//         below), diag is internally set. if mode = 2, diag
//         must contain positive entries that serve as
//         multiplicative scale factors for the variables.
//
//       mode is an integer input variable. if mode = 1, the
//         variables will be scaled internally. if mode = 2,
//         the scaling is specified by the input diag. other
//         values of mode are equivalent to mode = 1.
//
//       factor is a positive input variable used in determining the
//         initial step bound. this bound is set to the product of
//         factor and the euclidean norm of diag*x if nonzero, or else
//         to factor itself. in most cases factor should lie in the
//         interval (.1,100.). 100. is a generally recommended value.
//
//       nprint is an integer input variable that enables controlled
//         printing of iterates if it is positive. in this case,
//         fcn is called with iflag = 0 at the beginning of the first
//         iteration and every nprint iterations thereafter and
//         immediately prior to return, with x and fvec available
//         for printing. if nprint is not positive, no special calls
//         of fcn with iflag = 0 are made.
//
//       info is an integer output variable. if the user has
//         terminated execution, info is set to the (negative)
//         value of iflag. see description of fcn. otherwise,
//         info is set as follows.
//
//         info = 0   improper input parameters.
//
//         info = 1   relative error between two consecutive iterates
//                    is at most xtol.
//
//         info = 2   number of calls to fcn has reached or exceeded
//                    maxfev.
//
//         info = 3   xtol is too small. no further improvement in
//                    the approximate solution x is possible.
//
//         info = 4   iteration is not making good progress, as
//                    measured by the improvement from the last
//                    five jacobian evaluations.
//
//         info = 5   iteration is not making good progress, as
//                    measured by the improvement from the last
//                    ten iterations.
//
//       nfev is an integer output variable set to the number of
//         calls to fcn.
//
//       fjac is an output n by n array which contains the
//         orthogonal matrix q produced by the qr factorization
//         of the final approximate jacobian.
//
//       ldfjac is a positive integer input variable not less than n
//         which specifies the leading dimension of the array fjac.
//
//       r is an output array of length lr which contains the
//         upper triangular matrix produced by the qr factorization
//         of the final approximate jacobian, stored rowwise.
//
//       lr is a positive integer input variable not less than
//         (n*(n+1))/2.
//
//       qtf is an output array of length n which contains
//         the vector (q transpose)*fvec.
//
//       wa1, wa2, wa3, and wa4 are work arrays of length n.
//
{
    double actred;
    double delta;
    double epsmch;
    double fnorm;
    double fnorm1;
    int i;
    int iflag;
    int info;
    int iter;
    int iwa[1];
    int j;
    bool jeval;
    int l;
    int msum;
    int ncfail;
    int ncsuc;
    int nslow1;
    int nslow2;
    const double p001 = 0.001;
    const double p0001 = 0.0001;
    const double p1 = 0.1;
    const double p5 = 0.5;
    double pnorm;
    double prered;
    double ratio;
    double sum;
    double temp;
    double xnorm;
//
//  Certain loops in this function were kept closer to their original FORTRAN77
//  format, to avoid confusing issues with the array index L.  These loops are
//  marked "DO NOT ADJUST", although they certainly could be adjusted (carefully)
//  once the initial translated code is tested.
//

//
//  EPSMCH is the machine precision.
//
    epsmch = std::numeric_limits<double>::epsilon();

    info = 0;
    iflag = 0;
    nfev = 0;
//
//  Check the input parameters.
//
    if ( n <= 0 )
    {
        info = 0;
        return info;
    }
    if ( xtol < 0.0 )
    {
        info = 0;
        return info;
    }
    if ( maxfev <= 0 )
    {
        info = 0;
        return info;
    }
    if ( ml < 0 )
    {
        info = 0;
        return info;
    }
    if ( mu < 0 )
    {
        info = 0;
        return info;
    }
    if ( factor <= 0.0 )
    {
        info = 0;
        return info;
    }
    if ( ldfjac < n )
    {
        info = 0;
        return info;
    }
    if ( lr < ( n * ( n + 1 ) ) / 2 )
    {
        info = 0;
        return info;
    }
    if ( mode == 2 )
    {
        for ( j = 0; j < n; j++ )
        {
            if ( diag[j] <= 0.0 )
            {
                info = 0;
                return info;
            }
        }
    }
//
//  Evaluate the function at the starting point and calculate its norm.
//
    iflag = 1;
    fcn ( n, x, fvec, iflag, pars );
    nfev = 1;
    if ( iflag < 0 )
    {
        info = iflag;
        return info;
    }

    fnorm = enorm ( n, fvec );
//
//  Determine the number of calls to FCN needed to computeSynchrotronEmissivityAbsorptionAnalytic the jacobian matrix.
//
    if ( ml + mu + 1 < n )
    {
        msum = ml + mu + 1;
    }
    else
    {
        msum = n;
    }
//
//  Initialize iteration counter and monitors.
//
    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;
//
//  Beginning of the outer loop.
//
    for ( ; ; )
    {
        jeval = true;
//
//  Calculate the jacobian matrix.
//
        iflag = 2;
        fdjac1 ( fcn, n, x, pars, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn, wa1, wa2 );

        nfev = nfev + msum;
        if ( iflag < 0 )
        {
            info = iflag;
            return info;
        }
//
//  Compute the QR factorization of the jacobian.
//
        qrfac ( n, n, fjac, ldfjac, false, iwa, 1, wa1, wa2 );
//
//  On the first iteration and if MODE is 1, scale according
//  to the norms of the columns of the initial jacobian.
//
        if ( iter == 1 )
        {
            if ( mode == 1 )
            {
                for ( j = 0; j < n; j++ )
                {
                    if ( wa2[j] != 0.0 )
                    {
                        diag[j] = wa2[j];
                    }
                    else
                    {
                        diag[j] = 1.0;
                    }
                }
            }
//
//  On the first iteration, calculate the norm of the scaled X
//  and initialize the step bound DELTA.
//
            for ( j = 0; j < n; j++ )
            {
                wa3[j] = diag[j] * x[j];
            }
            xnorm = enorm ( n, wa3 );

            if ( xnorm == 0.0 )
            {
                delta = factor;
            }
            else
            {
                delta = factor * xnorm;
            }
        }
//
//  Form Q' * FVEC and store in QTF.
//
        for ( i = 0; i < n; i++ )
        {
            qtf[i] = fvec[i];
        }
        for ( j = 0; j < n; j++ )
        {
            if ( fjac[j+j*ldfjac] != 0.0 )
            {
                sum = 0.0;
                for ( i = j; i < n; i++ )
                {
                    sum = sum + fjac[i+j*ldfjac] * qtf[i];
                }
                temp = - sum / fjac[j+j*ldfjac];
                for ( i = j; i < n; i++ )
                {
                    qtf[i] = qtf[i] + fjac[i+j*ldfjac] * temp;
                }
            }
        }
//
//  Copy the triangular factor of the QR factorization into R.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
        for ( j = 1; j <= n; j++ )
        {
            l = j;
            for ( i = 1; i <= j - 1; i++ )
            {
                r[l-1] = fjac[(i-1)+(j-1)*ldfjac];
                l = l + n - i;
            }
            r[l-1] = wa1[j-1];
            if ( wa1[j-1] == 0.0 )
            {
                std::cout << "  Matrix is singular.\n";
            }
        }
//
//  Accumulate the orthogonal factor in FJAC.
//
        qform ( n, n, fjac, ldfjac );
//
//  Rescale if necessary.
//
        if ( mode == 1 )
        {
            for ( j = 0; j < n; j++ )
            {
                diag[j] = fmax ( diag[j], wa2[j] );
            }
        }
//
//  Beginning of the inner loop.
//
        for ( ; ; )
        {
//
//  If requested, call FCN to enable printing of iterates.
//
            if ( 0 < nprint )
            {
                if ( ( iter - 1 ) % nprint == 0 )
                {
                    iflag = 0;
                    fcn ( n, x, fvec, iflag, pars );
                    if ( iflag < 0 )
                    {
                        info = iflag;
                        return info;
                    }
                }
            }
//
//  Determine the direction P.
//
            dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 );
//
//  Store the direction P and X + P.  Calculate the norm of P.
//
            for ( j = 0; j < n; j++ )
            {
                wa1[j] = - wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            pnorm = enorm ( n, wa3 );
//
//  On the first iteration, adjust the initial step bound.
//
            if ( iter == 1 )
            {
                delta = fmin ( delta, pnorm );
            }
//
//  Evaluate the function at X + P and calculate its norm.
//
            iflag = 1;
            fcn ( n, wa2, wa4, iflag, pars );
            nfev = nfev + 1;
            if ( iflag < 0 )
            {
                info = iflag;
                return info;
            }
            fnorm1 = enorm ( n, wa4 );
//
//  Compute the scaled actual reduction.
//
            if ( fnorm1 < fnorm )
            {
                actred = 1.0 - ( fnorm1 / fnorm ) * ( fnorm1 / fnorm );
            }
            else
            {
                actred = - 1.0;
            }
//
//  Compute the scaled predicted reduction.
//
//  DO NOT ADJUST THIS LOOP, BECAUSE OF L.
//
            l = 1;
            for ( i = 1; i <= n; i++ )
            {
                sum = 0.0;
                for ( j = i; j <= n; j++ )
                {
                    sum = sum + r[l-1] * wa1[j-1];
                    l = l + 1;
                }
                wa3[i-1] = qtf[i-1] + sum;
            }
            temp = enorm ( n, wa3 );

            if ( temp < fnorm )
            {
                prered = 1.0 - ( temp / fnorm ) * ( temp / fnorm );
            }
            else
            {
                prered = 0.0;
            }
//
//  Compute the ratio of the actual to the predicted reduction.
//
            if ( 0.0 < prered )
            {
                ratio = actred / prered;
            }
            else
            {
                ratio = 0.0;
            }
//
//  Update the step bound.
//
            if ( ratio < p1 )
            {
                ncsuc = 0;
                ncfail = ncfail + 1;
                delta = p5 * delta;
            }
            else
            {
                ncfail = 0;
                ncsuc = ncsuc + 1;
                if ( p5 <= ratio || 1 < ncsuc )
                {
                    delta = fmax ( delta, pnorm / p5 );
                }
                if ( fabs ( ratio - 1.0 ) <= p1 )
                {
                    delta = pnorm / p5;
                }
            }
//
//  On successful iteration, update X, FVEC, and their norms.
//
            if ( p0001 <= ratio )
            {
                for ( j = 0; j < n; j++ )
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                xnorm = enorm ( n, wa2 );
                fnorm = fnorm1;
                iter = iter + 1;
            }
//
//  Determine the progress of the iteration.
//
            nslow1 = nslow1 + 1;
            if ( p001 <= actred )
            {
                nslow1 = 0;
            }
            if ( jeval )
            {
                nslow2 = nslow2 + 1;
            }
            if ( p1 <= actred )
            {
                nslow2 = 0;
            }
//
//  Test for convergence.
//
            if ( delta <= xtol * xnorm || fnorm == 0.0 )
            {
                info = 1;
                return info;
            }
//
//  Tests for termination and stringent tolerances.
//
            if ( maxfev <= nfev )
            {
                info = 2;
                return info;
            }
            if ( p1 * fmax ( p1 * delta, pnorm ) <= epsmch * xnorm )
            {
                info = 3;
                return info;
            }
            if ( nslow2 == 5 )
            {
                info = 4;
                return info;
            }
            if ( nslow1 == 10 )
            {
                info = 5;
                return info;
            }
//
//  Criterion for recalculating jacobian approximation by forward differences.
//
            if ( ncfail == 2 )
            {
                break;
            }
//
//  Calculate the rank one modification to the jacobian
//  and update QTF if necessary.
//
            for ( j = 0; j < n; j++ )
            {
                sum = 0.0;
                for ( i = 0; i < n; i++ )
                {
                    sum = sum + fjac[i+j*ldfjac] * wa4[i];
                }
                wa2[j] = ( sum - wa3[j] ) / pnorm;
                wa1[j] = diag[j] * ( ( diag[j] * wa1[j] ) / pnorm );
                if ( p0001 <= ratio )
                {
                    qtf[j] = sum;
                }
            }
//
//  Compute the QR factorization of the updated jacobian.
//
            r1updt ( n, n, r, lr, wa1, wa2, wa3 );
            r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 );
            r1mpyq ( 1, n, qtf, 1, wa2, wa3 );

            jeval = false;
        }
//
//  End of the inner loop.
//
    }
//
//  End of the outer loop.
//
}
//****************************************************************************80

void qform ( int m, int n, double q[], int ldq )

//****************************************************************************80
//
//  Purpose:
//
//    qform() constructs the standard form of Q from its factored form.
//
//  Discussion:
//
//    This function proceeds from the computed QR factorization of
//    an M by N matrix A to accumulate the M by M orthogonal matrix
//    Q from its factored form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2018
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of A, and the order of Q.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, double Q[LDQ*N].  On input, the full lower trapezoid in
//    the first min(M,N) columns of Q contains the factored form.
//    On output Q has been accumulated into a square matrix.
//
//    Input, int LDQ, the leading dimension of the array Q.
//
{
    int i;
    int j;
    int k;
    int minmn;
    double sum;
    double temp;
    double *wa;
//
//  Zero out the upper triangle of Q in the first min(M,N) columns.
//
    if ( m < n )
    {
        minmn = m;
    }
    else
    {
        minmn = n;
    }

    for ( j = 1; j < minmn; j++ )
    {
        for ( i = 0; i <= j - 1; i++ )
        {
            q[i+j*ldq] = 0.0;
        }
    }
//
//  Initialize remaining columns to those of the identity matrix.
//
    for ( j = n; j < m; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            q[i+j*ldq] = 0.0;
        }
        q[j+j*ldq] = 1.0;
    }
//
//  Accumulate Q from its factored form.
//
    wa = new double[m];

    for ( k = minmn - 1; 0 <= k; k-- )
    {
        for ( i = k; i < m; i++ )
        {
            wa[i] = q[i+k*ldq];
            q[i+k*ldq] = 0.0;
        }
        q[k+k*ldq] = 1.0;

        if ( wa[k] != 0.0 )
        {
            for ( j = k; j < m; j++ )
            {
                sum = 0.0;
                for ( i = k; i < m; i++ )
                {
                    sum = sum + q[i+j*ldq] * wa[i];
                }
                temp = sum / wa[k];
                for ( i = k; i < m; i++ )
                {
                    q[i+j*ldq] = q[i+j*ldq] - temp * wa[i];
                }
            }
        }
    }
//
//  Free memory.
//
    delete [] wa;

    return;
}
//****************************************************************************80

void qrfac ( int m, int n, double a[], int lda, bool pivot, int ipvt[],
             int lipvt, double rdiag[], double acnorm[] )

//****************************************************************************80
//
//  Purpose:
//
//    qrfac() computes the QR factorization of an M by N matrix.
//
//  Discussion:
//
//    This function uses Householder transformations with optional column
//    pivoting to computeSynchrotronEmissivityAbsorptionAnalytic a QR factorization of the M by N matrix A.
//
//    That is, QRFAC determines an orthogonal
//    matrix Q, a permutation matrix P, and an upper trapezoidal
//    matrix R with diagonal elements of nonincreasing magnitude,
//    such that A*P = Q*R.
//
//    The Householder transformation for
//    column k, k = 1,2,...,min(m,n), is of the form
//
//      i - (1/u(k))*u*u'
//
//    where U has zeros in the first K-1 positions.
//
//    The form of this transformation and the method of pivoting first
//    appeared in the corresponding LINPACK function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2017
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, double A[M*N].  On input, the matrix for which the QR
//    factorization is to be computed.  On output, the strict upper trapezoidal
//    part contains the strict upper trapezoidal part of the R factor, and
//    the lower trapezoidal part contains a factored form of Q, the non-trivial
//    elements of the U vectors described above.
//
//    Input, int LDA, a positive value not less than M which specifies the
//    leading dimension of the array A.
//
//    Input, bool PIVOT.  If true, then column pivoting is enforced.
//
//    Output, integer IPVT[LIPVT].  If PIVOT is true, then on output IPVT
//    defines the permutation matrix P such that A*P = Q*R.  Column J of P
//    is column IPVT[J] of the identity matrix.
//
//       lipvt is a positive integer input variable. if pivot is false,
//         then lipvt may be as small as 1. if pivot is true, then
//         lipvt must be at least n.
//
//       rdiag is an output array of length n which contains the
//         diagonal elements of r.
//
//       acnorm is an output array of length n which contains the
//         norms of the corresponding columns of the input matrix a.
//         if this information is not needed, then acnorm can coincide
//         with rdiag.
//
{
    double ajnorm;
    double epsmch;
    int i;
    int j;
    int k;
    int kmax;
    int minmn;
    const double p05 = 0.05;
    double sum;
    double temp;
    double *wa;
//
//  EPSMCH is the machine precision.
//
    epsmch = std::numeric_limits<double>::epsilon();
//
//  Compute the initial column norms and initialize several arrays.
//
    wa = new double[n];

    for ( j = 0; j < n; j++ )
    {
        acnorm[j] = enorm ( m, a+j*lda );
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        if ( pivot )
        {
            ipvt[j] = j;
        }
    }
//
//  Reduce A to R with Householder transformations.
//
    if ( m < n )
    {
        minmn = m;
    }
    else
    {
        minmn = n;
    }

    for ( j = 0; j < minmn; j++ )
    {
        if ( pivot )
        {
//
//  Bring the column of largest norm into the pivot position.
//
            kmax = j;
            for ( k = j; k < n; k++ )
            {
                if ( rdiag[kmax] < rdiag[k] )
                {
                    kmax = k;
                }
            }
            if ( kmax != j )
            {
                for ( i = 0; i < m; i++ )
                {
                    temp          = a[i+j*lda];
                    a[i+j*lda]    = a[i+kmax*lda];
                    a[i+kmax*lda] = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax]    = wa[j];
                k          = ipvt[j];
                ipvt[j]    = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }
//
//  Compute the Householder transformation to reduce the
//  J-th column of A to a multiple of the J-th unit vector.
//
        ajnorm = enorm ( m - j, a+j+j*lda );

        if ( ajnorm != 0.0 )
        {
            if ( a[j+j*lda] < 0.0 )
            {
                ajnorm = - ajnorm;
            }
            for ( i = j; i < m; i++ )
            {
                a[i+j*lda] = a[i+j*lda] / ajnorm;
            }
            a[j+j*lda] = a[j+j*lda] + 1.0;
//
//  Apply the transformation to the remaining columns and update the norms.
//
            for ( k = j + 1; k < n; k++ )
            {
                sum = 0.0;
                for ( i = j; i < m; i++ )
                {
                    sum = sum + a[i+j*lda] * a[i+k*lda];
                }
                temp = sum / a[j+j*lda];
                for ( i = j; i < m; i++ )
                {
                    a[i+k*lda] = a[i+k*lda] - temp * a[i+j*lda];
                }
                if ( pivot && rdiag[k] != 0.0 )
                {
                    temp = a[j+k*lda] / rdiag[k];
                    rdiag[k] = rdiag[k] * sqrt ( fmax ( 0.0, 1.0 - temp * temp ) );
                    if ( p05 * ( rdiag[k] / wa[k] ) * ( rdiag[k] / wa[k] ) <= epsmch )
                    {
                        rdiag[k] = enorm ( m - 1 - j, a+(j+1)+k*lda );
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = - ajnorm;
    }
//
//  Free memory.
//
    delete [] wa;

    return;
}
//****************************************************************************80

void r1mpyq ( int m, int n, double a[], int lda, double v[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    r1mpyq() multiplies an M by N matrix A by the Q factor.
//
//  Discussion:
//
//    Given an M by N matrix A, this function computes a*q where
//    q is the product of 2*(n - 1) transformations
//
//      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
//    and gv(i), gw(i) are givens rotations in the (i,n) plane which
//    eliminate elements in the i-th and n-th planes, respectively.
//
//    Q itself is not given, rather the information to recover the
//    GV and GW rotations is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, double A[M*N].  On input, the matrix to be postmultiplied
//    by the orthogonal matrix Q described above.  On output, the value of A*Q.
//
//    Input, int LDA, a positive value not less than M
//    which specifies the leading dimension of the array A.
//
//    Input, double V[N].  V(I) must contain the information necessary to
//    recover the givens rotation GV(I) described above.
//
//    Input, double W[N], contains the information necessary to recover the
//    Givens rotation gw(i) described above.
//
{
    double c;
    int i;
    int j;
    double s;
    double temp;
//
//  Apply the first set of Givens rotations to A.
//
    for ( j = n - 2; 0 <= j; j-- )
    {
        if ( 1.0 < fabs ( v[j] ) )
        {
            c = 1.0 / v[j];
            s = sqrt ( 1.0 - c * c );
        }
        else
        {
            s = v[j];
            c = sqrt ( 1.0 - s * s );
        }
        for ( i = 0; i < m; i++ )
        {
            temp           = c * a[i+j*lda] - s * a[i+(n-1)*lda];
            a[i+(n-1)*lda] = s * a[i+j*lda] + c * a[i+(n-1)*lda];
            a[i+j*lda]     = temp;
        }
    }
//
//  Apply the second set of Givens rotations to A.
//
    for ( j = 0; j < n - 1; j++ )
    {
        if ( 1.0 < fabs ( w[j] ) )
        {
            c = 1.0 / w[j];
            s = sqrt ( 1.0 - c * c );
        }
        else
        {
            s = w[j];
            c = sqrt ( 1.0 - s * s );
        }
        for ( i = 0; i < m; i++ )
        {
            temp           =   c * a[i+j*lda] + s * a[i+(n-1)*lda];
            a[i+(n-1)*lda] = - s * a[i+j*lda] + c * a[i+(n-1)*lda];
            a[i+j*lda]     = temp;
        }
    }

    return;
}
//****************************************************************************80

bool r1updt ( int m, int n, double s[], int ls, double u[], double v[],
              double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    r1updt() updates the Q factor after a rank one update of the matrix.
//
//  Discussion:
//
//    Given an M by N lower trapezoidal matrix S, an M-vector U,
//    and an N-vector V, the problem is to determine an
//    orthogonal matrix Q such that
//
//      (S + U*V') * Q
//
//    is again lower trapezoidal.
//
//    This function determines q as the product of 2*(n - 1) transformations
//
//      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
//    where gv(i), gw(i) are givens rotations in the (i,n) plane
//    which eliminate elements in the i-th and n-th planes,
//    respectively.
//
//    Q itself is not accumulated, rather the
//    information to recover the gv, gw rotations is returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2010
//
//  Author:
//
//    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    User Guide for MINPACK-1,
//    Technical Report ANL-80-74,
//    Argonne National Laboratory, 1980.
//
//  Parameters:
//
//    Input, int M, the number of rows of S.
//
//    Input, int N, the number of columns of S.  N must not exceed M.
//
//    Input/output, double S[LS].  On input, the lower trapezoidal matrix S
//    stored by columns.  On output, the lower trapezoidal matrix produced
//    as described above.
//
//    Input, int LS, a positive value not less than (N*(2*M-N+1))/2.
//
//    Input, double U[M], contains the vector U.
//
//    Input/output, double V[N].  On input, the vector vv.  On output,
//    information necessary to recover the givens rotation gv(i) described above.
//
//       w is an output array of length m. w(i) contains information
//         necessary to recover the givens rotation gw(i) described
//         above.
//
//       sing is a logical output variable. sing is set true if any
//         of the diagonal elements of the output s are zero. otherwise
//         sing is set false.
//
{
    double cotan;
    double cs;
    double giant;
    int i;
    int j;
    int jj;
    int l;
    int nm1;
    const double p25 = 0.25;
    const double p5 = 0.5;
    double sn;
    bool sing;
    double tan;
    double tau;
    double temp;
//
//  Because of the computation of the pointer JJ, this function was
//  converted from FORTRAN77 to C++ in a conservative way.  All computations
//  are the same, and only array indexing is adjusted.
//
//  GIANT is the largest magnitude.
//
    giant = std::numeric_limits<double>::max();
//
//  Initialize the diagonal element pointer.
//
    jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n );
//
//  Move the nontrivial part of the last column of S into W.
//
    l = jj;
    for ( i = n; i <= m; i++ )
    {
        w[i-1] = s[l-1];
        l = l + 1;
    }
//
//  Rotate the vector V into a multiple of the N-th unit vector
//  in such a way that a spike is introduced into W.
//
    nm1 = n - 1;

    for ( j = n - 1; 1 <= j; j-- )
    {
        jj = jj - ( m - j + 1 );
        w[j-1] = 0.0;

        if ( v[j-1] != 0.0 )
        {
//
//  Determine a Givens rotation which eliminates the J-th element of V.
//
            if ( fabs ( v[n-1] ) < fabs ( v[j-1] ) )
            {
                cotan = v[n-1] / v[j-1];
                sn = p5 / sqrt ( p25 + p25 * cotan * cotan );
                cs = sn * cotan;
                tau = 1.0;
                if ( 1.0 < fabs ( cs ) * giant )
                {
                    tau = 1.0 / cs;
                }
            }
            else
            {
                tan = v[j-1] / v[n-1];
                cs = p5 / sqrt ( p25 + p25 * tan * tan );
                sn = cs * tan;
                tau = sn;
            }
//
//  Apply the transformation to V and store the information
//  necessary to recover the Givens rotation.
//
            v[n-1] = sn * v[j-1] + cs * v[n-1];
            v[j-1] = tau;
//
//  Apply the transformation to S and extend the spike in W.
//
            l = jj;
            for ( i = j; i <= m; i++ )
            {
                temp   = cs * s[l-1] - sn * w[i-1];
                w[i-1] = sn * s[l-1] + cs * w[i-1];
                s[l-1] = temp;
                l = l + 1;
            }
        }
    }
//
//  Add the spike from the rank 1 update to W.
//
    for ( i = 1; i <= m; i++ )
    {
        w[i-1] = w[i-1] + v[n-1] * u[i-1];
    }
//
//  Eliminate the spike.
//
    sing = false;

    for ( j = 1; j <= nm1; j++ )
    {
//
//  Determine a Givens rotation which eliminates the
//  J-th element of the spike.
//
        if ( w[j-1] != 0.0 )
        {

            if ( fabs ( s[jj-1] ) < fabs ( w[j-1] ) )
            {
                cotan = s[jj-1] / w[j-1];
                sn = p5 / sqrt ( p25 + p25 * cotan * cotan );
                cs = sn * cotan;
                tau = 1.0;
                if ( 1.0 < fabs ( cs ) * giant )
                {
                    tau = 1.0 / cs;
                }
            }
            else
            {
                tan = w[j-1] / s[jj-1];
                cs = p5 / sqrt ( p25 + p25 * tan * tan );
                sn = cs * tan;
                tau = sn;
            }
//
//  Apply the transformation to s and reduce the spike in w.
//
            l = jj;

            for ( i = j; i <= m; i++ )
            {
                temp   =   cs * s[l-1] + sn * w[i-1];
                w[i-1] = - sn * s[l-1] + cs * w[i-1];
                s[l-1] = temp;
                l = l + 1;
            }
//
//  Store the information necessary to recover the givens rotation.
//
            w[j-1] = tau;
        }
//
//  Test for zero diagonal elements in the output s.
//
        if ( s[jj-1] == 0.0 )
        {
            sing = true;
        }
        jj = jj + ( m - j + 1 );
    }
//
//  Move W back into the last column of the output S.
//
    l = jj;
    for ( i = n; i <= m; i++ )
    {
        s[l-1] = w[i-1];
        l = l + 1;
    }
    if ( s[jj-1] == 0.0 )
    {
        sing = true;
    }
    return sing;
}

//****************************************************************************80

void r8vec2_print ( int n, double a1[], double a2[], std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec2_print() prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    int N, the number of components of the vector.
//
//    double A1[N], double A2[N], the vectors to be printed.
//
//    string TITLE, a title.
//
{
    int i;

    std::cout << "\n";
    std::cout << title << "\n";
    std::cout << "\n";
    for ( i = 0; i <= n - 1; i++ )
    {
        std::cout << std::setw(6)  << i
                  << ": " << std::setw(14) << a1[i]
                  << "  " << std::setw(14) << a2[i] << "\n";
    }

    return;
}

#endif //SRC_ROOTFINDERS_H
