/* MB04YD.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static doublereal c_b13 = -.125;
static integer c__1 = 1;
static doublereal c_b18 = 0.;
static doublereal c_b19 = 1.;

/* Subroutine */ int mb04yd_(char *jobu, char *jobv, integer *m, integer *n, 
	integer *rank, doublereal *theta, doublereal *q, doublereal *e, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, logical *
	inul, doublereal *tol, doublereal *reltol, doublereal *dwork, integer 
	*ldwork, integer *iwarn, integer *info, ftnlen jobu_len, ftnlen 
	jobv_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p, r__;
    static doublereal x;
    static integer i1;
    static doublereal eps;
    static logical noc12;
    static integer oldi, oldk;
    static doublereal cosl;
    static integer iter;
    static doublereal rmin, cosr, rmax, sinl, smax, sinr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical qrit;
    static integer info1;
    extern /* Subroutine */ int mb03md_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern integer mb03nd_(integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iascl;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb02ny_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *);
    static doublereal shift, sigmn;
    static integer maxit;
    extern /* Subroutine */ int mb04yw_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *);
    static doublereal sigmx;
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal thetac, safemn;
    static logical ljobua, ljobva;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobui, ljobvi;
    static integer numeig;
    static doublereal tolabs, thresh, pivmin, tolrel, smlnum;


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

/*     To partially diagonalize the bidiagonal matrix */

/*               |q(1) e(1)  0    ...       0      | */
/*               | 0   q(2) e(2)            .      | */
/*           J = | .                        .      |                  (1) */
/*               | .                  e(MIN(M,N)-1)| */
/*               | 0   ...        ...  q(MIN(M,N)) | */

/*     using QR or QL iterations in such a way that J is split into */
/*     unreduced bidiagonal submatrices whose singular values are either */
/*     all larger than a given bound or are all smaller than (or equal */
/*     to) this bound. The left- and right-hand Givens rotations */
/*     performed on J (corresponding to each QR or QL iteration step) may */
/*     be optionally accumulated in the arrays U and V. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the left-hand Givens rotations, as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is initialized to the M-by-MIN(M,N) submatrix of */
/*                     the unit matrix and the left-hand Givens rotations */
/*                     are accumulated in U; */
/*             = 'U':  The given matrix U is updated by the left-hand */
/*                     Givens rotations used in the calculation. */

/*     JOBV    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the right-hand Givens rotations, as follows: */
/*             = 'N':  Do not form V; */
/*             = 'I':  V is initialized to the N-by-MIN(M,N) submatrix of */
/*                     the unit matrix and the right-hand Givens */
/*                     rotations are accumulated in V; */
/*             = 'U':  The given matrix V is updated by the right-hand */
/*                     Givens rotations used in the calculation. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in matrix U.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows in matrix V.  N >= 0. */

/*     RANK    (input/output) INTEGER */
/*             On entry, if RANK < 0, then the rank of matrix J is */
/*             computed by the routine as the number of singular values */
/*             larger than THETA. */
/*             Otherwise, RANK must specify the rank of matrix J. */
/*             RANK <= MIN(M,N). */
/*             On exit, if RANK < 0 on entry, then RANK contains the */
/*             computed rank of J. That is, the number of singular */
/*             values of J larger than THETA. */
/*             Otherwise, the user-supplied value of RANK may be */
/*             changed by the routine on exit if the RANK-th and the */
/*             (RANK+1)-th singular values of J are considered to be */
/*             equal. See also the parameter TOL. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, if RANK < 0, then THETA must specify an upper */
/*             bound on the smallest singular values of J. THETA >= 0.0. */
/*             Otherwise, THETA must specify an initial estimate (t say) */
/*             for computing an upper bound such that precisely RANK */
/*             singular values are greater than this bound. */
/*             If THETA < 0.0, then t is computed by the routine. */
/*             On exit, if RANK >= 0 on entry, then THETA contains the */
/*             computed upper bound such that precisely RANK singular */
/*             values of J are greater than THETA + TOL. */
/*             Otherwise, THETA is unchanged. */

/*     Q       (input/output) DOUBLE PRECISION array, dimension */
/*             (MIN(M,N)) */
/*             On entry, this array must contain the diagonal elements */
/*             q(1),q(2),...,q(MIN(M,N)) of the bidiagonal matrix J. That */
/*             is, Q(i) = J(i,i) for i = 1,2,...,MIN(M,N). */
/*             On exit, this array contains the leading diagonal of the */
/*             transformed bidiagonal matrix J. */

/*     E       (input/output) DOUBLE PRECISION array, dimension */
/*             (MIN(M,N)-1) */
/*             On entry, this array must contain the superdiagonal */
/*             elements e(1),e(2),...,e(MIN(M,N)-1) of the bidiagonal */
/*             matrix J. That is, E(k) = J(k,k+1) for k = 1,2,..., */
/*             MIN(M,N)-1. */
/*             On exit, this array contains the superdiagonal of the */
/*             transformed bidiagonal matrix J. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, if JOBU = 'U', the leading M-by-MIN(M,N) part */
/*             of this array must contain a left transformation matrix */
/*             applied to the original matrix of the problem, and */
/*             on exit, the leading M-by-MIN(M,N) part of this array */
/*             contains the product of the input matrix U and the */
/*             left-hand Givens rotations. */
/*             On exit, if JOBU = 'I', then the leading M-by-MIN(M,N) */
/*             part of this array contains the matrix of accumulated */
/*             left-hand Givens rotations used. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. If JOBU = 'U' or */
/*             JOBU = 'I', LDU >= MAX(1,M); if JOBU = 'N', LDU >= 1. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             On entry, if JOBV = 'U', the leading N-by-MIN(M,N) part */
/*             of this array must contain a right transformation matrix */
/*             applied to the original matrix of the problem, and */
/*             on exit, the leading N-by-MIN(M,N) part of this array */
/*             contains the product of the input matrix V and the */
/*             right-hand Givens rotations. */
/*             On exit, if JOBV = 'I', then the leading N-by-MIN(M,N) */
/*             part of this array contains the matrix of accumulated */
/*             right-hand Givens rotations used. */
/*             If JOBV = 'N', the array V is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDV = 1 and */
/*             declare this array to be V(1,1) in the calling program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. If JOBV = 'U' or */
/*             JOBV = 'I', LDV >= MAX(1,N); if JOBV = 'N', LDV >= 1. */

/*     INUL    (input/output) LOGICAL array, dimension (MIN(M,N)) */
/*             On entry, the leading MIN(M,N) elements of this array must */
/*             be set to .FALSE. unless the i-th columns of U (if JOBU = */
/*             'U') and V (if JOBV = 'U') already contain a computed base */
/*             vector of the desired singular subspace of the original */
/*             matrix, in which case INUL(i) must be set to .TRUE. */
/*             for 1 <= i <= MIN(M,N). */
/*             On exit, the indices of the elements of this array with */
/*             value .TRUE. indicate the indices of the diagonal entries */
/*             of J which belong to those bidiagonal submatrices whose */
/*             singular values are all less than or equal to THETA. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             This parameter defines the multiplicity of singular values */
/*             by considering all singular values within an interval of */
/*             length TOL as coinciding. TOL is used in checking how many */
/*             singular values are less than or equal to THETA. Also in */
/*             computing an appropriate upper bound THETA by a bisection */
/*             method, TOL is used as a stopping criterion defining the */
/*             minimum (absolute) subinterval width. TOL is also taken */
/*             as an absolute tolerance for negligible elements in the */
/*             QR/QL iterations. If the user sets TOL to be less than or */
/*             equal to 0, then the tolerance is taken as */
/*             EPS * MAX(ABS(Q(i)), ABS(E(k))), where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH), */
/*             i = 1,2,...,MIN(M,N) and k = 1,2,...,MIN(M,N)-1. */

/*     RELTOL  DOUBLE PRECISION */
/*             This parameter specifies the minimum relative width of an */
/*             interval. When an interval is narrower than TOL, or than */
/*             RELTOL times the larger (in magnitude) endpoint, then it */
/*             is considered to be sufficiently small and bisection has */
/*             converged. If the user sets RELTOL to be less than */
/*             BASE * EPS, where BASE is machine radix and EPS is machine */
/*             precision (see LAPACK Library routine DLAMCH), then the */
/*             tolerance is taken as BASE * EPS. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,6*MIN(M,N)-5), if JOBU = 'I' or 'U', or */
/*                                               JOBV = 'I' or 'U'; */
/*             LDWORK >= MAX(1,4*MIN(M,N)-3), if JOBU = 'N' and */
/*                                               JOBV = 'N'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  if the rank of the bidiagonal matrix J (as specified */
/*                   by the user) has been lowered because a singular */
/*                   value of multiplicity larger than 1 was found. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; this includes values like RANK > MIN(M,N), or */
/*                   THETA < 0.0 and RANK < 0; */
/*             = 1:  if the maximum number of QR/QL iteration steps */
/*                   (30*MIN(M,N)) has been exceeded. */

/*     METHOD */

/*     If the upper bound THETA is not specified by the user, then it is */
/*     computed by the routine (using a bisection method) such that */
/*     precisely (MIN(M,N) - RANK) singular values of J are less than or */
/*     equal to THETA + TOL. */

/*     The method used by the routine (see [1]) then proceeds as follows. */

/*     The unreduced bidiagonal submatrices of J(j), where J(j) is the */
/*     transformed bidiagonal matrix after the j-th iteration step, are */
/*     classified into the following three classes: */

/*     - C1 contains the bidiagonal submatrices with all singular values */
/*       > THETA, */
/*     - C2 contains the bidiagonal submatrices with all singular values */
/*       <= THETA and */
/*     - C3 contains the bidiagonal submatrices with singular values */
/*       > THETA and also singular values <= THETA. */

/*     If C3 is empty, then the partial diagonalization is complete, and */
/*     RANK is the sum of the dimensions of the bidiagonal submatrices of */
/*     C1. */
/*     Otherwise, QR or QL iterations are performed on each bidiagonal */
/*     submatrix of C3, until this bidiagonal submatrix has been split */
/*     into two bidiagonal submatrices. These two submatrices are then */
/*     classified and the iterations are restarted. */
/*     If the upper left diagonal element of the bidiagonal submatrix is */
/*     larger than its lower right diagonal element, then QR iterations */
/*     are performed, else QL iterations are used. The shift is taken as */
/*     the smallest diagonal element of the bidiagonal submatrix (in */
/*     magnitude) unless its value exceeds THETA, in which case it is */
/*     taken as zero. */

/*     REFERENCES */

/*     [1] Van Huffel, S., Vandewalle, J. and Haegemans, A. */
/*         An efficient and reliable algorithm for computing the */
/*         singular subspace of a matrix associated with its smallest */
/*         singular values. */
/*         J. Comput. and Appl. Math., 19, pp. 313-330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     To avoid overflow, matrix J is scaled so that its largest element */
/*     is no greater than  overflow**(1/2) * underflow**(1/4) in absolute */
/*     value (and not much smaller than that, for maximal accuracy). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04QD by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     July 10, 1997. V. Sima. */
/*     November 25, 1997. V. Sima: Setting INUL(K) = .TRUE. when handling */
/*                                 2-by-2 submatrix. */

/*     KEYWORDS */

/*     Bidiagonal matrix, orthogonal transformation, singular values. */

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
    --q;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --inul;
    --dwork;

    /* Function Body */
    p = min(*m,*n);
    *info = 0;
    *iwarn = 0;
    ljobui = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
    ljobvi = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
    ljobua = ljobui || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
    ljobva = ljobvi || lsame_(jobv, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ljobua && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ljobva && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*rank > p) {
	*info = -5;
    } else if (*rank < 0 && *theta < 0.) {
	*info = -6;
    } else if (! ljobua && *ldu < 1 || ljobua && *ldu < max(1,*m)) {
	*info = -10;
    } else if (! ljobva && *ldv < 1 || ljobva && *ldv < max(1,*n)) {
	*info = -12;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = p * 6 - 5;
/* Computing MAX */
	i__3 = 1, i__4 = (p << 2) - 3;
	if ((ljobua || ljobva) && *ldwork < max(i__1,i__2) || ! (ljobua || 
		ljobva) && *ldwork < max(i__3,i__4)) {
	    *info = -17;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB04YD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (p == 0) {
	if (*rank >= 0) {
	    *theta = 0.;
	}
	*rank = 0;
	return 0;
    }

/*     Set tolerances and machine parameters. */

    tolabs = *tol;
    tolrel = *reltol;
    smax = (d__1 = q[p], abs(d__1));

    i__1 = p - 1;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__3 = smax, d__4 = (d__1 = q[j], abs(d__1)), d__3 = max(d__3,d__4), 
		d__4 = (d__2 = e[j], abs(d__2));
	smax = max(d__3,d__4);
/* L20: */
    }

    safemn = dlamch_("Safe minimum", (ftnlen)12);
    eps = dlamch_("Epsilon", (ftnlen)7);
    if (tolabs <= 0.) {
	tolabs = eps * smax;
    }
    x = dlamch_("Base", (ftnlen)4) * eps;
    if (tolrel <= x) {
	tolrel = x;
    }
/* Computing MAX */
/* Computing MIN */
    d__3 = 100., d__4 = pow_dd(&eps, &c_b13);
    d__1 = 10., d__2 = min(d__3,d__4);
    thresh = max(d__1,d__2) * eps;
    smlnum = safemn / eps;
    rmin = sqrt(smlnum);
/* Computing MIN */
    d__1 = 1. / rmin, d__2 = 1. / sqrt(sqrt(safemn));
    rmax = min(d__1,d__2);
    thetac = *theta;

/*     Scale the matrix to allowable range, if necessary, and set PIVMIN, */
/*     using the squares of Q and E (saved in DWORK). */

    iascl = 0;
    if (smax > 0. && smax < rmin) {
	iascl = 1;
	sigma = rmin / smax;
    } else if (smax > rmax) {
	iascl = 1;
	sigma = rmax / smax;
    }
    if (iascl == 1) {
	dscal_(&p, &sigma, &q[1], &c__1);
	i__1 = p - 1;
	dscal_(&i__1, &sigma, &e[1], &c__1);
	thetac = sigma * *theta;
	tolabs = sigma * tolabs;
    }

/* Computing 2nd power */
    d__1 = q[p];
    pivmin = d__1 * d__1;
    dwork[p] = pivmin;

    i__1 = p - 1;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = q[j];
	dwork[j] = d__1 * d__1;
/* Computing 2nd power */
	d__1 = e[j];
	dwork[p + j] = d__1 * d__1;
/* Computing MAX */
	d__1 = pivmin, d__2 = dwork[j], d__1 = max(d__1,d__2), d__2 = dwork[p 
		+ j];
	pivmin = max(d__1,d__2);
/* L40: */
    }

/* Computing MAX */
    d__1 = pivmin * safemn;
    pivmin = max(d__1,safemn);

/*     Initialize U and/or V to the identity matrix, if needed. */

    if (ljobui) {
	dlaset_("Full", m, &p, &c_b18, &c_b19, &u[u_offset], ldu, (ftnlen)4);
    }
    if (ljobvi) {
	dlaset_("Full", n, &p, &c_b18, &c_b19, &v[v_offset], ldv, (ftnlen)4);
    }

/*     Estimate THETA (if not fixed by the user), and set R. */

    if (*rank >= 0) {
	j = p - *rank;
	mb03md_(&p, &j, &thetac, &q[1], &e[1], &dwork[1], &dwork[p + 1], &
		pivmin, &tolabs, &tolrel, iwarn, &info1);
	*theta = thetac;
	if (iascl == 1) {
	    *theta /= sigma;
	}
	if (j <= 0) {
	    return 0;
	}
	r__ = p - j;
    } else {
	r__ = p - mb03nd_(&p, &thetac, &dwork[1], &dwork[p + 1], &pivmin, &
		info1);
    }

    *rank = p;

    i__1 = p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (inul[i__]) {
	    --(*rank);
	}
/* L60: */
    }

/*     From now on K is the smallest known index such that the elements */
/*     of the bidiagonal matrix J with indices larger than K belong to C1 */
/*     or C2. */
/*     RANK = P - SUM(dimensions of known bidiagonal matrices of C2). */

    k = p;
    oldi = -1;
    oldk = -1;
    iter = 0;
    maxit = p * 30;
/*     WHILE ( C3 NOT EMPTY ) DO */
L80:
    if (*rank > r__ && k > 0) {
/*        WHILE ( K.GT.0 .AND. INUL(K) ) DO */

/*        Search for the rightmost index of a bidiagonal submatrix, */
/*        not yet classified. */

L100:
	if (k > 0) {
	    if (inul[k]) {
		--k;
		goto L100;
	    }
	}
/*        END WHILE 100 */

	if (k == 0) {
	    return 0;
	}

	noc12 = TRUE_;
/*        WHILE ((ITER < MAXIT).AND.(No bidiagonal matrix of C1 or */
/*                C2 found)) DO */
L120:
	if (iter < maxit && noc12) {

/*           Search for negligible Q(I) or E(I-1) (for I > 1) and find */
/*           the shift. */

	    i__ = k;
	    x = (d__1 = q[i__], abs(d__1));
	    shift = x;
/*           WHILE ABS( Q(I) ) > TOLABS .AND. ABS( E(I-1) ) > TOLABS ) DO */
L140:
	    if (i__ > 1) {
		if (x > tolabs && (d__1 = e[i__ - 1], abs(d__1)) > tolabs) {
		    --i__;
		    x = (d__1 = q[i__], abs(d__1));
		    if (x < shift) {
			shift = x;
		    }
		    goto L140;
		}
	    }
/*           END WHILE 140 */

/*           Classify the bidiagonal submatrix (of order J) found. */

	    j = k - i__ + 1;
	    if (x <= tolabs || k == i__) {
		noc12 = FALSE_;
	    } else {
		numeig = mb03nd_(&j, &thetac, &dwork[i__], &dwork[p + i__], &
			pivmin, &info1);
		if (numeig >= j || numeig <= 0) {
		    noc12 = FALSE_;
		}
	    }
	    if (noc12) {
		if (j == 2) {

/*                 Handle separately the 2-by-2 submatrix. */

		    dlasv2_(&q[i__], &e[i__], &q[k], &sigmn, &sigmx, &sinr, &
			    cosr, &sinl, &cosl);
		    q[i__] = sigmx;
		    q[k] = sigmn;
		    e[i__] = 0.;
		    --(*rank);
		    inul[k] = TRUE_;
		    noc12 = FALSE_;

/*                 Update U and/or V, if needed. */

		    if (ljobua) {
			drot_(m, &u[i__ * u_dim1 + 1], &c__1, &u[k * u_dim1 + 
				1], &c__1, &cosl, &sinl);
		    }
		    if (ljobva) {
			drot_(n, &v[i__ * v_dim1 + 1], &c__1, &v[k * v_dim1 + 
				1], &c__1, &cosr, &sinr);
		    }
		} else {

/*                 If working on new submatrix, choose QR or */
/*                 QL iteration. */

		    if (i__ != oldi || k != oldk) {
			qrit = (d__1 = q[i__], abs(d__1)) >= (d__2 = q[k], 
				abs(d__2));
		    }
		    oldi = i__;
		    if (qrit) {
			if ((d__2 = e[k - 1], abs(d__2)) <= thresh * (d__1 = 
				q[k], abs(d__1))) {
			    e[k - 1] = 0.;
			}
		    } else {
			if ((d__2 = e[i__], abs(d__2)) <= thresh * (d__1 = q[
				i__], abs(d__1))) {
			    e[i__] = 0.;
			}
		    }

		    mb04yw_(&qrit, &ljobua, &ljobva, m, n, &i__, &k, &shift, &
			    q[1], &e[1], &u[u_offset], ldu, &v[v_offset], ldv,
			     &dwork[p * 2]);

		    if (qrit) {
			if ((d__1 = e[k - 1], abs(d__1)) <= tolabs) {
			    e[k - 1] = 0.;
			}
		    } else {
			if ((d__1 = e[i__], abs(d__1)) <= tolabs) {
			    e[i__] = 0.;
			}
		    }
/* Computing 2nd power */
		    d__1 = q[k];
		    dwork[k] = d__1 * d__1;

		    i__1 = k - 1;
		    for (i1 = i__; i1 <= i__1; ++i1) {
/* Computing 2nd power */
			d__1 = q[i1];
			dwork[i1] = d__1 * d__1;
/* Computing 2nd power */
			d__1 = e[i1];
			dwork[p + i1] = d__1 * d__1;
/* L160: */
		    }

		    ++iter;
		}
	    }
	    goto L120;
	}
/*        END WHILE 120 */

	if (iter >= maxit) {
	    *info = 1;
	    goto L200;
	}

	if (x <= tolabs) {

/*           Split at negligible diagonal element ABS( Q(I) ) <= TOLABS. */

	    mb02ny_(&ljobua, &ljobva, m, n, &i__, &k, &q[1], &e[1], &u[
		    u_offset], ldu, &v[v_offset], ldv, &dwork[p * 2]);
	    inul[i__] = TRUE_;
	    --(*rank);
	} else {

/*           A negligible superdiagonal element ABS( E(I-1) ) <= TOL */
/*           has been found, the corresponding bidiagonal submatrix */
/*           belongs to C1 or C2. Treat this bidiagonal submatrix. */

	    if (j >= 2) {
		if (numeig == j) {

		    i__1 = k;
		    for (i1 = i__; i1 <= i__1; ++i1) {
			inul[i1] = TRUE_;
/* L180: */
		    }

		    *rank -= j;
		    k -= j;
		} else {
		    k = i__ - 1;
		}
	    } else {
		if (x <= thetac + tolabs) {
		    inul[i__] = TRUE_;
		    --(*rank);
		}
		--k;
	    }
	    oldk = k;
	}
	goto L80;
    }
/*     END WHILE 80 */

/*     If matrix was scaled, then rescale Q and E appropriately. */

L200:
    if (iascl == 1) {
	d__1 = 1. / sigma;
	dscal_(&p, &d__1, &q[1], &c__1);
	i__1 = p - 1;
	d__1 = 1. / sigma;
	dscal_(&i__1, &d__1, &e[1], &c__1);
    }

    return 0;
/* *** Last line of MB04YD *** */
} /* mb04yd_ */

