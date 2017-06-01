/* SB10ZP.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b6 = 1.;
static integer c_n1 = -1;
static doublereal c_b35 = -1.;

/* Subroutine */ int sb10zp_(integer *discfl, integer *n, doublereal *a, 
	integer *lda, doublereal *b, doublereal *c__, doublereal *d__, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, iwa, imp, rep, imz, iwp, iwq, rez, idw1, idw2, idw3, 
	    ldw1, iwps, iwqs, info2;
    extern /* Subroutine */ int ab04md_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     ab07nd_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), td04ad_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);
    static doublereal scalb, scalc, scald;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), mc01pd_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), dgeev_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal rcond;
    static integer index[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer maxwrk;


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

/*     To transform a SISO (single-input single-output) system [A,B;C,D] */
/*     by mirroring its unstable poles and zeros in the boundary of the */
/*     stability domain, thus preserving the frequency response of the */
/*     system, but making it stable and minimum phase. Specifically, for */
/*     a continuous-time system, the positive real parts of its poles */
/*     and zeros are exchanged with their negatives. Discrete-time */
/*     systems are first converted to continuous-time systems using a */
/*     bilinear transformation, and finally converted back. */

/*     ARGUMENTS */

/*     Input/Output parameters */

/*     DISCFL  (input) INTEGER */
/*             Indicates the type of the system, as follows: */
/*             = 0: continuous-time system; */
/*             = 1: discrete-time system. */

/*     N       (input/output) INTEGER */
/*             On entry, the order of the original system.  N >= 0. */
/*             On exit, the order of the transformed, minimal system. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original system matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed matrix A, in an upper Hessenberg form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the original system */
/*             vector B. */
/*             On exit, this array contains the transformed vector B. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the original system */
/*             vector C. */
/*             On exit, this array contains the transformed vector C. */
/*             The first N-1 elements are zero (for the exit value of N). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (1) */
/*             On entry, this array must contain the original system */
/*             scalar D. */
/*             On exit, this array contains the transformed scalar D. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2,N+1) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max(N*N + 5*N, 6*N + 1 + min(1,N)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the discrete --> continuous transformation cannot */
/*                   be made; */
/*             = 2:  if the system poles cannot be found; */
/*             = 3:  if the inverse system cannot be found, i.e., D is */
/*                   (close to) zero; */
/*             = 4:  if the system zeros cannot be found; */
/*             = 5:  if the state-space representation of the new */
/*                   transfer function T(s) cannot be found; */
/*             = 6:  if the continuous --> discrete transformation cannot */
/*                   be made. */

/*     METHOD */

/*     First, if the system is discrete-time, it is transformed to */
/*     continuous-time using alpha = beta = 1 in the bilinear */
/*     transformation implemented in the SLICOT routine AB04MD. */
/*     Then the eigenvalues of A, i.e., the system poles, are found. */
/*     Then, the inverse of the original system is found and its poles, */
/*     i.e., the system zeros, are evaluated. */
/*     The obtained system poles Pi and zeros Zi are checked and if a */
/*     positive real part is detected, it is exchanged by -Pi or -Zi. */
/*     Then the polynomial coefficients of the transfer function */
/*     T(s) = Q(s)/P(s) are found. */
/*     The state-space representation of T(s) is then obtained. */
/*     The system matrices B, C, D are scaled so that the transformed */
/*     system has the same system gain as the original system. */
/*     If the original system is discrete-time, then the result (which is */
/*     continuous-time) is converted back to discrete-time. */

/*     CONTRIBUTORS */

/*     Asparuh Markovski, Technical University of Sofia, July 2003. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003. */

/*     KEYWORDS */

/*     Bilinear transformation, stability, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */

/*     Test input parameters and workspace. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --c__;
    --d__;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*discfl != 0 && *discfl != 1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = *n * *n + *n * 5, i__2 = *n * 6 + 1 + min(1,*n);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -10;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB10ZP", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Workspace usage 1. */

    rep = 1;
    imp = rep + *n;
    rez = imp + *n;
    imz = rez + *n;
    iwa = rez;
    idw1 = iwa + *n * *n;
    ldw1 = *ldwork - idw1 + 1;

/*     1. Discrete --> continuous transformation if needed. */

    if (*discfl == 1) {

/*        Workspace:  need    max(1,N); */
/*                            prefer  larger. */

	ab04md_("D", n, &c__1, &c__1, &c_b6, &c_b6, &a[a_offset], lda, &b[1], 
		lda, &c__[1], &c__1, &d__[1], &c__1, &iwork[1], &dwork[1], 
		ldwork, &info2, (ftnlen)1);
	if (info2 != 0) {
	    *info = 1;
	    return 0;
	}
	maxwrk = (integer) dwork[1];
    } else {
	maxwrk = 0;
    }

/*     2. Determine the factors for restoring system gain. */

    scald = d__[1];
    scalc = sqrt((abs(scald)));
    scalb = d_sign(&scalc, &scald);

/*     3. Find the system poles, i.e., the eigenvalues of A. */
/*        Workspace:  need    N*N + 2*N + 3*N; */
/*                            prefer  larger. */

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iwa], n, (ftnlen)4);

    dgeev_("N", "N", n, &dwork[iwa], n, &dwork[rep], &dwork[imp], &dwork[idw1]
	    , &c__1, &dwork[idw1], &c__1, &dwork[idw1], &ldw1, &info2, (
	    ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 2;
	return 0;
    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) (dwork[idw1] + idw1 - 1);
    maxwrk = max(i__1,i__2);

/*     4. Compute the inverse system [Ai, Bi; Ci, Di]. */
/*        Workspace:  need    N*N + 2*N + 4; */
/*                            prefer  larger. */

    ab07nd_(n, &c__1, &a[a_offset], lda, &b[1], lda, &c__[1], &c__1, &d__[1], 
	    &c__1, &rcond, &iwork[1], &dwork[idw1], &ldw1, &info2);
    if (info2 != 0) {
	*info = 3;
	return 0;
    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) (dwork[idw1] + idw1 - 1);
    maxwrk = max(i__1,i__2);

/*     5. Find the system zeros, i.e., the eigenvalues of Ai. */
/*        Workspace:  need    4*N + 3*N; */
/*                            prefer  larger. */

    idw1 = imz + *n;
    ldw1 = *ldwork - idw1 + 1;

    dgeev_("N", "N", n, &a[a_offset], lda, &dwork[rez], &dwork[imz], &dwork[
	    idw1], &c__1, &dwork[idw1], &c__1, &dwork[idw1], &ldw1, &info2, (
	    ftnlen)1, (ftnlen)1);
    if (info2 != 0) {
	*info = 4;
	return 0;
    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) (dwork[idw1] + idw1 - 1);
    maxwrk = max(i__1,i__2);

/*     6. Exchange the zeros and the poles with positive real parts with */
/*        their negatives. */

    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (dwork[rep + i__] > 0.) {
	    dwork[rep + i__] = -dwork[rep + i__];
	}
	if (dwork[rez + i__] > 0.) {
	    dwork[rez + i__] = -dwork[rez + i__];
	}
/* L10: */
    }

/*     Workspace usage 2. */

    iwp = idw1;
    idw2 = iwp + *n + 1;
    iwps = 1;

/*     7. Construct the nominator and the denominator */
/*        of the system transfer function T( s ) = Q( s )/P( s ). */
/*     8. Rearrange the coefficients in Q(s) and P(s) because */
/*        MC01PD subroutine produces them in increasing powers of s. */
/*        Workspace:  need    6*N + 2. */

    mc01pd_(n, &dwork[rep], &dwork[imp], &dwork[iwp], &dwork[idw2], &info2);
    i__1 = *n + 1;
    dcopy_(&i__1, &dwork[iwp], &c_n1, &dwork[iwps], &c__1);

/*     Workspace usage 3. */

    iwq = idw1;
    iwqs = iwps + *n + 1;
    idw3 = iwqs + *n + 1;

    mc01pd_(n, &dwork[rez], &dwork[imz], &dwork[iwq], &dwork[idw2], &info2);
    i__1 = *n + 1;
    dcopy_(&i__1, &dwork[iwq], &c_n1, &dwork[iwqs], &c__1);

/*     9. Make the conversion T(s) --> [A, B; C, D]. */
/*        Workspace:  need    2*N + 2 + N + max(N,3); */
/*                            prefer  larger. */

    index[0] = *n;
    i__1 = *ldwork - idw3 + 1;
    td04ad_("R", &c__1, &c__1, index, &dwork[iwps], &c__1, &dwork[iwqs], &
	    c__1, &c__1, n, &a[a_offset], lda, &b[1], lda, &c__[1], &c__1, &
	    d__[1], &c__1, &c_b35, &iwork[1], &dwork[idw3], &i__1, &info2, (
	    ftnlen)1);
    if (info2 != 0) {
	*info = 5;
	return 0;
    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) (dwork[idw3] + idw3 - 1);
    maxwrk = max(i__1,i__2);

/*    10. Scale the transformed system to the previous gain. */

    if (*n > 0) {
	dscal_(n, &scalb, &b[1], &c__1);
	c__[*n] = scalc * c__[*n];
    }

    d__[1] = scald;

/*     11. Continuous --> discrete transformation if needed. */

    if (*discfl == 1) {
	ab04md_("C", n, &c__1, &c__1, &c_b6, &c_b6, &a[a_offset], lda, &b[1], 
		lda, &c__[1], &c__1, &d__[1], &c__1, &iwork[1], &dwork[1], 
		ldwork, &info2, (ftnlen)1);
	if (info2 != 0) {
	    *info = 6;
	    return 0;
	}
    }

    dwork[1] = (doublereal) maxwrk;
    return 0;

/* *** Last line of SB10ZP *** */
} /* sb10zp_ */

