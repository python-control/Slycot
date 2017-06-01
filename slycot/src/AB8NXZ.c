/* AB8NXZ.f -- translated by f2c (version 20100827).
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

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;

/* Subroutine */ int ab8nxz_(integer *n, integer *m, integer *p, integer *ro, 
	integer *sigma, doublereal *svlmax, doublecomplex *abcd, integer *
	ldabcd, integer *ninfz, integer *infz, integer *kronl, integer *mu, 
	integer *nu, integer *nkrol, doublereal *tol, integer *iwork, 
	doublereal *dwork, doublecomplex *zwork, integer *lzwork, integer *
	info)
{
    /* System generated locals */
    integer abcd_dim1, abcd_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i1, nb, ik;
    static doublecomplex tc;
    static integer np, iz, mm1, ro1, mpm, tau, mnu, rank, itau;
    static doublereal sval[3];
    static integer irow, mntau, jwork;
    extern /* Subroutine */ int mb3oyz_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    integer *), mb3pyz_(integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zlaset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zlapmt_(logical *, integer *,
	     integer *, doublecomplex *, integer *, integer *);
    static logical lquery;
    extern /* Subroutine */ int zlatzm_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static integer wrkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmrq_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

/*     To extract from the (N+P)-by-(M+N) system */
/*                  ( B  A ) */
/*                  ( D  C ) */
/*     an (NU+MU)-by-(M+NU) "reduced" system */
/*                  ( B' A') */
/*                  ( D' C') */
/*     having the same transmission zeros but with D' of full row rank. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of state variables.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     RO      (input/output) INTEGER */
/*             On entry, */
/*             = P     for the original system; */
/*             = MAX(P-M, 0) for the pertransposed system. */
/*             On exit, RO contains the last computed rank. */

/*     SIGMA   (input/output) INTEGER */
/*             On entry, */
/*             = 0  for the original system; */
/*             = M  for the pertransposed system. */
/*             On exit, SIGMA contains the last computed value sigma in */
/*             the algorithm. */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             During each reduction step, the rank-revealing QR */
/*             factorization of a matrix stops when the estimated minimum */
/*             singular value is smaller than TOL * MAX(SVLMAX,EMSV), */
/*             where EMSV is the estimated maximum singular value. */
/*             SVLMAX >= 0. */

/*     ABCD    (input/output) COMPLEX*16 array, dimension (LDABCD,M+N) */
/*             On entry, the leading (N+P)-by-(M+N) part of this array */
/*             must contain the compound input matrix of the system. */
/*             On exit, the leading (NU+MU)-by-(M+NU) part of this array */
/*             contains the reduced compound input matrix of the system. */

/*     LDABCD  INTEGER */
/*             The leading dimension of array ABCD. */
/*             LDABCD >= MAX(1,N+P). */

/*     NINFZ   (input/output) INTEGER */
/*             On entry, the currently computed number of infinite zeros. */
/*             It should be initialized to zero on the first call. */
/*             NINFZ >= 0. */
/*             On exit, the number of infinite zeros. */

/*     INFZ    (input/output) INTEGER array, dimension (N) */
/*             On entry, INFZ(i) must contain the current number of */
/*             infinite zeros of degree i, where i = 1,2,...,N, found in */
/*             the previous call(s) of the routine. It should be */
/*             initialized to zero on the first call. */
/*             On exit, INFZ(i) contains the number of infinite zeros of */
/*             degree i, where i = 1,2,...,N. */

/*     KRONL   (input/output) INTEGER array, dimension (N+1) */
/*             On entry, this array must contain the currently computed */
/*             left Kronecker (row) indices found in the previous call(s) */
/*             of the routine. It should be initialized to zero on the */
/*             first call. */
/*             On exit, the leading NKROL elements of this array contain */
/*             the left Kronecker (row) indices. */

/*     MU      (output) INTEGER */
/*             The normal rank of the transfer function matrix of the */
/*             original system. */

/*     NU      (output) INTEGER */
/*             The dimension of the reduced system matrix and the number */
/*             of (finite) invariant zeros if D' is invertible. */

/*     NKROL   (output) INTEGER */
/*             The number of left Kronecker indices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             NOTE that when SVLMAX > 0, the estimated ranks could be */
/*             less than those defined above (see SVLMAX). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX(M,P)) */

/*     DWORK   DOUBLE PRECISION array, dimension (2*MAX(M,P)) */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ). */
/*             For optimum performance LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Svaricek, F. */
/*         Computation of the Structural Invariants of Linear */
/*         Multivariable Systems with an Extended Version of */
/*         the Program ZEROS. */
/*         System & Control Letters, 6, pp. 261-266, 1985. */

/*     [2] Emami-Naeini, A. and Van Dooren, P. */
/*         Computation of Zeros of Linear Multivariable Systems. */
/*         Automatica, 18, pp. 415-430, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008 with suggestions from P. Gahinet, */
/*     The MathWorks. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, unitary transformation, structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    abcd_dim1 = *ldabcd;
    abcd_offset = 1 + abcd_dim1;
    abcd -= abcd_offset;
    --infz;
    --kronl;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    np = *n + *p;
    mpm = min(*p,*m);
    *info = 0;
    lquery = *lzwork == -1;

/*     Test the input scalar arguments. */

    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = *p - *m;
	if (*ro != *p && *ro != max(i__1,0)) {
	    *info = -4;
	} else if (*sigma != 0 && *sigma != *m) {
	    *info = -5;
	} else if (*svlmax < 0.) {
	    *info = -6;
	} else if (*ldabcd < max(1,np)) {
	    *info = -8;
	} else if (*ninfz < 0) {
	    *info = -9;
	} else {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *m * 3 - 1;
/* Computing MAX */
	    i__4 = *p * 3 - 1, i__4 = max(i__4,np), i__5 = *n + *m;
	    i__1 = 1, i__2 = mpm + max(i__3,*n), i__1 = max(i__1,i__2), i__2 =
		     min(*p,*n) + max(i__4,i__5);
	    jwork = max(i__1,i__2);
	    if (lquery) {
		if (*m > 0) {
/* Computing MIN */
		    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", p, n, &
			    mpm, &c_n1, (ftnlen)6, (ftnlen)2);
		    nb = min(i__1,i__2);
/* Computing MAX */
		    i__1 = jwork, i__2 = mpm + max(1,*n) * nb;
		    wrkopt = max(i__1,i__2);
		} else {
		    wrkopt = jwork;
		}
/* Computing MIN */
		i__3 = min(*p,*n);
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", &np, n, &
			i__3, &c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
		i__1 = wrkopt, i__2 = min(*p,*n) + max(1,np) * nb;
		wrkopt = max(i__1,i__2);
/* Computing MIN */
		i__3 = *m + *n;
		i__4 = min(*p,*n);
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "LN", n, &i__3, &
			i__4, &c_n1, (ftnlen)6, (ftnlen)2);
		nb = min(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
		i__3 = 1, i__4 = *m + *n;
		i__1 = wrkopt, i__2 = min(*p,*n) + max(i__3,i__4) * nb;
		wrkopt = max(i__1,i__2);
	    } else if (*lzwork < jwork) {
		*info = -19;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB8NXZ", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
	return 0;
    }

    *mu = *p;
    *nu = *n;

    iz = 0;
    ik = 1;
    mm1 = *m + 1;
    itau = 1;
    *nkrol = 0;
    wrkopt = 1;

/*     Main reduction loop: */

/*            M   NU                  M     NU */
/*      NU  [ B   A ]           NU  [ B     A ] */
/*      MU  [ D   C ]  -->    SIGMA [ RD   C1 ]   (SIGMA = rank(D) = */
/*                             TAU  [ 0    C2 ]    row size of RD) */

/*                                    M   NU-RO  RO */
/*                            NU-RO [ B1   A11  A12 ] */
/*                     -->      RO  [ B2   A21  A22 ]  (RO = rank(C2) = */
/*                            SIGMA [ RD   C11  C12 ]   col size of LC) */
/*                             TAU  [ 0     0   LC  ] */

/*                                     M   NU-RO */
/*                            NU-RO [ B1   A11 ]     NU := NU - RO */
/*                                  [----------]     MU := RO + SIGMA */
/*                     -->      RO  [ B2   A21 ]      D := [B2;RD] */
/*                            SIGMA [ RD   C11 ]      C := [A21;C11] */

L20:
    if (*mu == 0) {
	goto L80;
    }

/*     (Note: Comments in the code beginning "xWorkspace:", where x is */
/*     I, D, or C, describe the minimal amount of integer, real and */
/*     complex workspace needed at that point in the code, respectively, */
/*     as well as the preferred amount for good performance.) */

    ro1 = *ro;
    mnu = *m + *nu;
    if (*m > 0) {
	if (*sigma != 0) {
	    irow = *nu + 1;

/*           Compress rows of D.  First exploit triangular shape. */
/*           CWorkspace: need   M+N-1. */

	    i__1 = *sigma;
	    for (i1 = 1; i1 <= i__1; ++i1) {
		i__2 = *ro + 1;
		zlarfg_(&i__2, &abcd[irow + i1 * abcd_dim1], &abcd[irow + 1 + 
			i1 * abcd_dim1], &c__1, &tc);
		i__2 = *ro + 1;
		i__3 = mnu - i1;
		d_cnjg(&z__1, &tc);
		zlatzm_("L", &i__2, &i__3, &abcd[irow + 1 + i1 * abcd_dim1], &
			c__1, &z__1, &abcd[irow + (i1 + 1) * abcd_dim1], &
			abcd[irow + 1 + (i1 + 1) * abcd_dim1], ldabcd, &zwork[
			1], (ftnlen)1);
		++irow;
/* L40: */
	    }
	    i__1 = *ro + *sigma - 1;
	    zlaset_("Lower", &i__1, sigma, &c_b1, &c_b1, &abcd[*nu + 2 + 
		    abcd_dim1], ldabcd, (ftnlen)5);
	}

/*        Continue with Householder with column pivoting. */

/*        The rank of D is the number of (estimated) singular values */
/*        that are greater than TOL * MAX(SVLMAX,EMSV). This number */
/*        includes the singular values of the first SIGMA columns. */
/*        IWorkspace: need   M; */
/*        RWorkspace: need   2*M; */
/*        CWorkspace: need   min(RO1,M) + 3*M - 1.  RO1 <= P. */

	if (*sigma < *m) {
	    jwork = itau + min(ro1,*m);
	    i1 = *sigma + 1;
	    irow = *nu + i1;
	    i__1 = *m - *sigma;
	    mb3oyz_(&ro1, &i__1, &abcd[irow + i1 * abcd_dim1], ldabcd, tol, 
		    svlmax, &rank, sval, &iwork[1], &zwork[itau], &dwork[1], &
		    zwork[jwork], info);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = jwork + *m * 3 - 2;
	    wrkopt = max(i__1,i__2);

/*           Apply the column permutations to matrices B and part of D. */

	    i__1 = *nu + *sigma;
	    i__2 = *m - *sigma;
	    zlapmt_(&c_true, &i__1, &i__2, &abcd[i1 * abcd_dim1 + 1], ldabcd, 
		    &iwork[1]);

	    if (rank > 0) {

/*              Apply the Householder transformations to the submatrix C. */
/*              CWorkspace:    need   min(RO1,M) + NU; */
/*                             prefer min(RO1,M) + NU*NB. */

		i__1 = *lzwork - jwork + 1;
		zunmqr_("Left", "Conjugate", &ro1, nu, &rank, &abcd[irow + i1 
			* abcd_dim1], ldabcd, &zwork[itau], &abcd[irow + mm1 *
			 abcd_dim1], ldabcd, &zwork[jwork], &i__1, info, (
			ftnlen)4, (ftnlen)9);
/* Computing MAX */
		i__3 = jwork;
		i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
		wrkopt = max(i__1,i__2);
		if (ro1 > 1) {
		    i__1 = ro1 - 1;
/* Computing MIN */
		    i__3 = ro1 - 1;
		    i__2 = min(i__3,rank);
		    zlaset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &abcd[irow + 
			    1 + i1 * abcd_dim1], ldabcd, (ftnlen)5);
		}
		ro1 -= rank;
	    }
	}
    }

    tau = ro1;
    *sigma = *mu - tau;

/*     Determination of the orders of the infinite zeros. */

    if (iz > 0) {
	infz[iz] = infz[iz] + *ro - tau;
	*ninfz += iz * (*ro - tau);
    }
    if (ro1 == 0) {
	goto L80;
    }
    ++iz;

    if (*nu <= 0) {
	*mu = *sigma;
	*nu = 0;
	*ro = 0;
    } else {

/*        Compress the columns of C2 using RQ factorization with row */
/*        pivoting, P * C2 = R * Q. */

	i1 = *nu + *sigma + 1;
	mntau = min(tau,*nu);
	jwork = itau + mntau;

/*        The rank of C2 is the number of (estimated) singular values */
/*        greater than TOL * MAX(SVLMAX,EMSV). */
/*        IWorkspace: need TAU; */
/*        RWorkspace: need 2*TAU; */
/*        CWorkspace: need min(TAU,NU) + 3*TAU - 1. */

	mb3pyz_(&tau, nu, &abcd[i1 + mm1 * abcd_dim1], ldabcd, tol, svlmax, &
		rank, sval, &iwork[1], &zwork[itau], &dwork[1], &zwork[jwork],
		 info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = jwork + tau * 3 - 1;
	wrkopt = max(i__1,i__2);
	if (rank > 0) {
	    irow = i1 + tau - rank;

/*           Apply Q' to the first NU columns of [A; C1] from the right. */
/*           CWorkspace: need   min(TAU,NU) + NU + SIGMA; SIGMA <= P; */
/*                       prefer min(TAU,NU) + (NU  + SIGMA)*NB. */

	    i__1 = i1 - 1;
	    i__2 = *lzwork - jwork + 1;
	    zunmrq_("Right", "ConjTranspose", &i__1, nu, &rank, &abcd[irow + 
		    mm1 * abcd_dim1], ldabcd, &zwork[mntau - rank + 1], &abcd[
		    mm1 * abcd_dim1 + 1], ldabcd, &zwork[jwork], &i__2, info, 
		    (ftnlen)5, (ftnlen)13);
/* Computing MAX */
	    i__3 = jwork;
	    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
	    wrkopt = max(i__1,i__2);

/*           Apply Q to the first NU rows and M + NU columns of [ B  A ] */
/*           from the left. */
/*           CWorkspace: need   min(TAU,NU) + M + NU; */
/*                       prefer min(TAU,NU) + (M + NU)*NB. */

	    i__1 = *lzwork - jwork + 1;
	    zunmrq_("Left", "NoTranspose", nu, &mnu, &rank, &abcd[irow + mm1 *
		     abcd_dim1], ldabcd, &zwork[mntau - rank + 1], &abcd[
		    abcd_offset], ldabcd, &zwork[jwork], &i__1, info, (ftnlen)
		    4, (ftnlen)11);
/* Computing MAX */
	    i__3 = jwork;
	    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
	    wrkopt = max(i__1,i__2);

	    i__1 = *nu - rank;
	    zlaset_("Full", &rank, &i__1, &c_b1, &c_b1, &abcd[irow + mm1 * 
		    abcd_dim1], ldabcd, (ftnlen)4);
	    if (rank > 1) {
		i__1 = rank - 1;
		i__2 = rank - 1;
		zlaset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &abcd[irow + 1 + 
			(mm1 + *nu - rank) * abcd_dim1], ldabcd, (ftnlen)5);
	    }
	}

	*ro = rank;
    }

/*     Determine the left Kronecker indices (row indices). */

    kronl[ik] = kronl[ik] + tau - *ro;
    *nkrol += kronl[ik];
    ++ik;

/*     C and D are updated to [A21 ; C11] and [B2 ; RD]. */

    *nu -= *ro;
    *mu = *sigma + *ro;
    if (*ro != 0) {
	goto L20;
    }

L80:
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
    return 0;
/* *** Last line of AB8NXZ *** */
} /* ab8nxz_ */

