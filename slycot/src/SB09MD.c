/* SB09MD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sb09md_(integer *n, integer *nc, integer *nb, doublereal 
	*h1, integer *ldh1, doublereal *h2, integer *ldh2, doublereal *ss, 
	integer *ldss, doublereal *se, integer *ldse, doublereal *pre, 
	integer *ldpre, doublereal *tol, integer *info)
{
    /* System generated locals */
    integer h1_dim1, h1_offset, h2_dim1, h2_offset, pre_dim1, pre_offset, 
	    se_dim1, se_offset, ss_dim1, ss_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal var, sse, sss, vare, epso, toler;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noflow;


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

/*     To compare two multivariable sequences M1(k) and M2(k) for */
/*     k = 1,2,...,N, and evaluate their closeness. Each of the */
/*     parameters M1(k) and M2(k) is an NC by NB matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of parameters.  N >= 0. */

/*     NC      (input) INTEGER */
/*             The number of rows in M1(k) and M2(k).  NC >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns in M1(k) and M2(k).  NB >= 0. */

/*     H1      (input) DOUBLE PRECISION array, dimension (LDH1,N*NB) */
/*             The leading NC-by-N*NB part of this array must contain */
/*             the multivariable sequence M1(k), where k = 1,2,...,N. */
/*             Each parameter M1(k) is an NC-by-NB matrix, whose */
/*             (i,j)-th element must be stored in H1(i,(k-1)*NB+j) for */
/*             i = 1,2,...,NC and j = 1,2,...,NB. */

/*     LDH1    INTEGER */
/*             The leading dimension of array H1.  LDH1 >= MAX(1,NC). */

/*     H2      (input) DOUBLE PRECISION array, dimension (LDH2,N*NB) */
/*             The leading NC-by-N*NB part of this array must contain */
/*             the multivariable sequence M2(k), where k = 1,2,...,N. */
/*             Each parameter M2(k) is an NC-by-NB matrix, whose */
/*             (i,j)-th element must be stored in H2(i,(k-1)*NB+j) for */
/*             i = 1,2,...,NC and j = 1,2,...,NB. */

/*     LDH2    INTEGER */
/*             The leading dimension of array H2.  LDH2 >= MAX(1,NC). */

/*     SS      (output) DOUBLE PRECISION array, dimension (LDSS,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             matrix SS. */

/*     LDSS    INTEGER */
/*             The leading dimension of array SS.  LDSS >= MAX(1,NC). */

/*     SE      (output) DOUBLE PRECISION array, dimension (LDSE,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             quadratic error matrix SE. */

/*     LDSE    INTEGER */
/*             The leading dimension of array SE.  LDSE >= MAX(1,NC). */

/*     PRE     (output) DOUBLE PRECISION array, dimension (LDPRE,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             percentage relative error matrix PRE. */

/*     LDPRE   INTEGER */
/*             The leading dimension of array PRE.  LDPRE >= MAX(1,NC). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in the computation of the error */
/*             matrices SE and PRE. If the user sets TOL to be less than */
/*             EPS then the tolerance is taken as EPS, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The (i,j)-th element of the matrix SS is defined by: */
/*                        N          2 */
/*               SS    = SUM  M1  (k) .                            (1) */
/*                 ij    k=1    ij */

/*     The (i,j)-th element of the quadratic error matrix SE is defined */
/*     by: */
/*                        N                      2 */
/*               SE    = SUM  (M1  (k) - M2  (k)) .                (2) */
/*                 ij    k=1     ij        ij */

/*     The (i,j)-th element of the percentage relative error matrix PRE */
/*     is defined by: */

/*               PRE   = 100 x SQRT( SE  / SS  ).                  (3) */
/*                  ij                 ij    ij */

/*     The following precautions are taken by the routine to guard */
/*     against underflow and overflow: */

/*     (i) if ABS( M1  (k) ) > 1/TOL or ABS( M1  (k) - M2  (k) ) > 1/TOL, */
/*                   ij                        ij        ij */

/*         then SE   and SS   are set to 1/TOL and PRE   is set to 1; and */
/*                ij       ij                         ij */

/*     (ii) if ABS( SS  ) <= TOL, then PRE   is set to 100. */
/*                    ij                  ij */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */
/*        2xNBxNCx(N+1) multiplications/divisions, */
/*        4xNBxNCxN     additions/subtractions and */
/*          NBxNC       square roots. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB09AD by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Closeness multivariable sequences, elementary matrix operations, */
/*     real signals, system response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    h1_dim1 = *ldh1;
    h1_offset = 1 + h1_dim1;
    h1 -= h1_offset;
    h2_dim1 = *ldh2;
    h2_offset = 1 + h2_dim1;
    h2 -= h2_offset;
    ss_dim1 = *ldss;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    se_dim1 = *ldse;
    se_offset = 1 + se_dim1;
    se -= se_offset;
    pre_dim1 = *ldpre;
    pre_offset = 1 + pre_dim1;
    pre -= pre_offset;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*nc < 0) {
	*info = -2;
    } else if (*nb < 0) {
	*info = -3;
    } else if (*ldh1 < max(1,*nc)) {
	*info = -5;
    } else if (*ldh2 < max(1,*nc)) {
	*info = -7;
    } else if (*ldss < max(1,*nc)) {
	*info = -9;
    } else if (*ldse < max(1,*nc)) {
	*info = -11;
    } else if (*ldpre < max(1,*nc)) {
	*info = -13;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB09MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *nc == 0 || *nb == 0) {
	return 0;
    }

/* Computing MAX */
    d__1 = *tol, d__2 = dlamch_("Epsilon", (ftnlen)7);
    toler = max(d__1,d__2);
    epso = 1. / toler;

    i__1 = *nb;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *nc;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sse = 0.;
	    sss = 0.;
	    noflow = TRUE_;
	    k = 0;

/*           WHILE ( ( NOFLOW .AND. ( K .LT. N*NB ) ) DO */
L20:
	    if (noflow && k < *n * *nb) {
		var = h1[i__ + (k + j) * h1_dim1];
		vare = h2[i__ + (k + j) * h2_dim1] - var;
		if (abs(var) > epso || abs(vare) > epso) {
		    se[i__ + j * se_dim1] = epso;
		    ss[i__ + j * ss_dim1] = epso;
		    pre[i__ + j * pre_dim1] = 1.;
		    noflow = FALSE_;
		} else {
		    if (abs(vare) > toler) {
			sse += vare * vare;
		    }
		    if (abs(var) > toler) {
			sss += var * var;
		    }
		    k += *nb;
		}
		goto L20;
	    }
/*           END WHILE 20 */

	    if (noflow) {
		se[i__ + j * se_dim1] = sse;
		ss[i__ + j * ss_dim1] = sss;
		pre[i__ + j * pre_dim1] = 100.;
		if (sss > toler) {
		    pre[i__ + j * pre_dim1] = sqrt(sse / sss) * 100.;
		}
	    }
/* L40: */
	}

/* L60: */
    }

    return 0;
/* *** Last line of SB09MD *** */
} /* sb09md_ */

