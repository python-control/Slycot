/* MD03BB.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int md03bb_(char *cond, integer *n, integer *ipar, integer *
	lipar, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag,
	 doublereal *qtb, doublereal *delta, doublereal *par, integer *ranks, 
	doublereal *x, doublereal *rx, doublereal *tol, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset;

    /* Local variables */
    extern /* Subroutine */ int md03by_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen);


/*     SLICOT RELEASE 5.0. */

/*     Copyright (c) 2002-2009 NICONET e.V. */

/*     This program is free software: you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation, either version 2 of */
/*     the License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see */
/*     <http://www.gnu.org/licenses/>. */

/*     PURPOSE */

/*     To determine a value for the parameter PAR such that if x solves */
/*     the system */

/*           A*x = b ,     sqrt(PAR)*D*x = 0 , */

/*     in the least squares sense, where A is an m-by-n matrix, D is an */
/*     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if */
/*     DELTA is a positive number, DXNORM is the Euclidean norm of D*x, */
/*     then either PAR is zero and */

/*           ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*     It is assumed that a QR factorization, with column pivoting, of A */
/*     is available, that is, A*P = Q*R, where P is a permutation matrix, */
/*     Q has orthogonal columns, and R is an upper triangular matrix */
/*     with diagonal elements of nonincreasing magnitude. */
/*     The routine needs the full upper triangle of R, the permutation */
/*     matrix P, and the first n components of Q'*b (' denotes the */
/*     transpose). On output, MD03BB also provides an upper triangular */
/*     matrix S such that */

/*           P'*(A'*A + PAR*D*D)*P = S'*S . */

/*     Matrix S is used in the solution process. */

/*     This routine is an interface to SLICOT Library routine MD03BY, */
/*     for solving standard nonlinear least squares problems using SLICOT */
/*     routine MD03BD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the matrices R and S */
/*             should be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation for R and S; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R and S for zero values; */
/*             = 'U' :  use the rank already stored in RANKS (for R). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix R. IPAR and LIPAR are not used by this routine, */
/*             but are provided for compatibility with SLICOT Library */
/*             routine MD03BD. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR, N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the full upper triangle is unaltered, and the */
/*             strict lower triangle contains the strict upper triangle */
/*             (transposed) of the upper triangular matrix S. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             A*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             An upper bound on the Euclidean norm of D*x.  DELTA > 0. */

/*     PAR     (input/output) DOUBLE PRECISION */
/*             On entry, PAR must contain an initial estimate of the */
/*             Levenberg-Marquardt parameter.  PAR >= 0. */
/*             On exit, it contains the final estimate of this parameter. */

/*     RANKS   (input or output) INTEGER array, dimension (1) */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical rank of the matrix R. */
/*             On exit, this array contains the numerical rank of the */
/*             matrix S. */
/*             RANKS is defined as an array for compatibility with SLICOT */
/*             Library routine MD03BD. */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system A*x = b, sqrt(PAR)*D*x = 0. */

/*     RX      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the matrix-vector product -R*P'*x. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             rank of the matrices R and S. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 4*N, if COND =  'E'; */
/*             LDWORK >= 2*N, if COND <> 'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine calls SLICOT Library routine MD03BY to perform the */
/*     calculations. */

/*     FURTHER COMMENTS */

/*     For efficiency, the arguments are not checked. This is done in */
/*     the routine MD03BY (except for LIPAR). */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipar;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --ipvt;
    --diag;
    --qtb;
    --ranks;
    --x;
    --rx;
    --dwork;

    /* Function Body */
    md03by_(cond, n, &r__[r_offset], ldr, &ipvt[1], &diag[1], &qtb[1], delta, 
	    par, &ranks[1], &x[1], &rx[1], tol, &dwork[1], ldwork, info, (
	    ftnlen)1);
    return 0;

/* *** Last line of MD03BB *** */
} /* md03bb_ */

