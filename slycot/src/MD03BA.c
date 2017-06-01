/* MD03BA.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int md03ba_(integer *n, integer *ipar, integer *lipar, 
	doublereal *fnorm, doublereal *j, integer *ldj, doublereal *e, 
	doublereal *jnorms, doublereal *gnorm, integer *ipvt, doublereal *
	dwork, integer *ldwork, integer *info)
{
    extern /* Subroutine */ int md03bx_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);


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

/*     To compute the QR factorization with column pivoting of an */
/*     m-by-n Jacobian matrix J (m >= n), that is, J*P = Q*R, where Q is */
/*     a matrix with orthogonal columns, P a permutation matrix, and */
/*     R an upper trapezoidal matrix with diagonal elements of */
/*     nonincreasing magnitude, and to apply the transformation Q' on */
/*     the error vector e (in-situ). The 1-norm of the scaled gradient */
/*     is also returned. */

/*     This routine is an interface to SLICOT Library routine MD03BX, */
/*     for solving standard nonlinear least squares problems using SLICOT */
/*     routine MD03BD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J.  N >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain the number of rows M of the Jacobian */
/*                     matrix J.  M >= N. */
/*             IPAR is provided for compatibility with SLICOT Library */
/*             routine MD03BD. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 1. */

/*     FNORM   (input) DOUBLE PRECISION */
/*             The Euclidean norm of the vector e.  FNORM >= 0. */

/*     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the Jacobian matrix J. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*             array contains the upper triangular factor R of the */
/*             Jacobian matrix. Note that for efficiency of the later */
/*             calculations, the matrix R is delivered with the leading */
/*             dimension MAX(1,N), possibly much smaller than the value */
/*             of LDJ on entry. */

/*     LDJ     (input/output) INTEGER */
/*             The leading dimension of array J. */
/*             On entry, LDJ >= MAX(1,M). */
/*             On exit,  LDJ >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (M) */
/*             On entry, this array must contain the error vector e. */
/*             On exit, this array contains the updated vector Q'*e. */

/*     JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the Euclidean norms of the columns */
/*             of the Jacobian matrix, considered in the initial order. */

/*     GNORM   (output) DOUBLE PRECISION */
/*             If FNORM > 0, the 1-norm of the scaled vector */
/*             J'*Q'*e/FNORM, with each element i further divided */
/*             by JNORMS(i) (if JNORMS(i) is nonzero). */
/*             If FNORM = 0, the returned value of GNORM is 0. */

/*     IPVT    (output) INTEGER array, dimension (N) */
/*             This array defines the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1,      if N = 0 or  M = 1; */
/*             LDWORK >= 4*N+1,  if N > 1. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine calls SLICOT Library routine MD03BX to perform the */
/*     calculations. */

/*     FURTHER COMMENTS */

/*     For efficiency, the arguments are not checked. This is done in */
/*     the routine MD03BX (except for LIPAR). */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, Jacobian matrix, matrix algebra, */
/*     matrix operations. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --dwork;
    --ipvt;
    --jnorms;
    --e;
    --j;
    --ipar;

    /* Function Body */
    md03bx_(&ipar[1], n, fnorm, &j[1], ldj, &e[1], &jnorms[1], gnorm, &ipvt[1]
	    , &dwork[1], ldwork, info);
    return 0;

/* *** Last line of MD03BA *** */
} /* md03ba_ */

