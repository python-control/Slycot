/* MB03ZD.f -- translated by f2c (version 20100827).
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

static doublereal c_b26 = 1.;
static doublereal c_b34 = -1.;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b272 = 0.;

/* Subroutine */ int mb03zd_(char *which, char *meth, char *stab, char *
	balanc, char *ortbal, logical *select, integer *n, integer *mm, 
	integer *ilo, doublereal *scale, doublereal *s, integer *lds, 
	doublereal *t, integer *ldt, doublereal *g, integer *ldg, doublereal *
	u1, integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1, 
	integer *ldv1, doublereal *v2, integer *ldv2, integer *m, doublereal *
	wr, doublereal *wi, doublereal *us, integer *ldus, doublereal *uu, 
	integer *lduu, logical *lwork, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen which_len, ftnlen meth_len, 
	ftnlen stab_len, ftnlen balanc_len, ftnlen ortbal_len)
{
    /* System generated locals */
    integer g_dim1, g_offset, s_dim1, s_offset, t_dim1, t_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, us_dim1, us_offset, uu_dim1, 
	    uu_offset, v1_dim1, v1_offset, v2_dim1, v2_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, pw, pdw;
    static logical lus, luu, lbef, lbal, lall, pair;
    static integer ierr;
    static doublereal temp;
    static logical lext;
    extern /* Subroutine */ int mb04di_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb03td_(char *, char *, logical *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), mb03za_(char *, char *, char *, char *
	    , char *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ux_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), dgeqp3_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin, wrkopt;


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

/*     To compute the stable and unstable invariant subspaces for a */
/*     Hamiltonian matrix with no eigenvalues on the imaginary axis, */
/*     using the output of the SLICOT Library routine MB03XD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WHICH   CHARACTER*1 */
/*             Specifies the cluster of eigenvalues for which the */
/*             invariant subspaces are computed: */
/*             = 'A':  select all n eigenvalues; */
/*             = 'S':  select a cluster of eigenvalues specified by */
/*                     SELECT. */

/*     METH    CHARACTER*1 */
/*             If WHICH = 'A' this parameter specifies the method to be */
/*             used for computing bases of the invariant subspaces: */
/*             = 'S':  compute the n-dimensional basis from a set of */
/*                     n vectors; */
/*             = 'L':  compute the n-dimensional basis from a set of */
/*                     2*n vectors. */
/*             When in doubt, use METH = 'S'. In some cases, METH = 'L' */
/*             may result in more accurately computed invariant */
/*             subspaces, see [1]. */

/*     STAB    CHARACTER*1 */
/*             Specifies the type of invariant subspaces to be computed: */
/*             = 'S':  compute the stable invariant subspace, i.e., the */
/*                     invariant subspace belonging to those selected */
/*                     eigenvalues that have negative real part; */
/*             = 'U':  compute the unstable invariant subspace, i.e., */
/*                     the invariant subspace belonging to those */
/*                     selected eigenvalues that have positive real */
/*                     part; */
/*             = 'B':  compute both the stable and unstable invariant */
/*                     subspaces. */

/*     BALANC  CHARACTER*1 */
/*             Specifies the type of inverse balancing transformation */
/*             required: */
/*             = 'N':  do nothing; */
/*             = 'P':  do inverse transformation for permutation only; */
/*             = 'S':  do inverse transformation for scaling only; */
/*             = 'B':  do inverse transformations for both permutation */
/*                     and scaling. */
/*             BALANC must be the same as the argument BALANC supplied to */
/*             MB03XD. Note that if the data is further post-processed, */
/*             e.g., for solving an algebraic Riccati equation, it is */
/*             recommended to delay inverse balancing (in particular the */
/*             scaling part) and apply it to the final result only, */
/*             see [2]. */

/*     ORTBAL  CHARACTER*1 */
/*             If BALANC <> 'N', this option specifies how inverse */
/*             balancing is applied to the computed invariant subspaces: */
/*             = 'B':  apply inverse balancing before orthogonal bases */
/*                     for the invariant subspaces are computed; */
/*             = 'A':  apply inverse balancing after orthogonal bases */
/*                     for the invariant subspaces have been computed; */
/*                     this may yield non-orthogonal bases if */
/*                     BALANC = 'S' or BALANC = 'B'. */

/*     SELECT  (input) LOGICAL array, dimension (N) */
/*             If WHICH = 'S', SELECT specifies the eigenvalues */
/*             corresponding to the positive and negative square */
/*             roots of the eigenvalues of S*T in the selected cluster. */
/*             To select a real eigenvalue w(j), SELECT(j) must be set */
/*             to .TRUE.. To select a complex conjugate pair of */
/*             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2 */
/*             diagonal block, both SELECT(j) and SELECT(j+1) must be set */
/*             to .TRUE.; a complex conjugate pair of eigenvalues must be */
/*             either both included in the cluster or both excluded. */
/*             This array is not referenced if WHICH = 'A'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices S, T and G. N >= 0. */

/*     MM      (input) INTEGER */
/*             The number of columns in the arrays US and/or UU. */
/*             If WHICH = 'A' and METH = 'S',  MM >=   N; */
/*             if WHICH = 'A' and METH = 'L',  MM >= 2*N; */
/*             if WHICH = 'S',                 MM >=   M. */
/*             The minimal values above for MM give the numbers of */
/*             vectors to be used for computing a basis for the */
/*             invariant subspace(s). */

/*     ILO     (input) INTEGER */
/*             If BALANC <> 'N', then ILO is the integer returned by */
/*             MB03XD.  1 <= ILO <= N+1. */

/*     SCALE   (input) DOUBLE PRECISION array, dimension (N) */
/*             If BALANC <> 'N', the leading N elements of this array */
/*             must contain details of the permutation and scaling */
/*             factors, as returned by MB03XD. */
/*             This array is not referenced if BALANC = 'N'. */

/*     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix S in real Schur form. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S.  LDS >= max(1,N). */

/*     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix T. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= max(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, if METH = 'L', the leading N-by-N part of this */
/*             array must contain a general matrix G. */
/*             On exit, if METH = 'L', the leading N-by-N part of this */
/*             array is overwritten. */
/*             This array is not referenced if METH = 'S'. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= 1. */
/*             LDG >= max(1,N) if METH = 'L'. */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (1,1) block of an orthogonal symplectic */
/*             matrix U. */
/*             On exit, this array is overwritten. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= MAX(1,N). */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (2,1) block of an orthogonal symplectic */
/*             matrix U. */
/*             On exit, this array is overwritten. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= MAX(1,N). */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (1,1) block of an orthogonal symplectic */
/*             matrix V. */
/*             On exit, this array is overwritten. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= MAX(1,N). */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (2,1) block of an orthogonal symplectic */
/*             matrix V. */
/*             On exit, this array is overwritten. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= MAX(1,N). */

/*     M       (output) INTEGER */
/*             The number of selected eigenvalues. */

/*     WR      (output) DOUBLE PRECISION array, dimension (M) */
/*     WI      (output) DOUBLE PRECISION array, dimension (M) */
/*             On exit, the leading M elements of WR and WI contain the */
/*             real and imaginary parts, respectively, of the selected */
/*             eigenvalues that have nonpositive real part. Complex */
/*             conjugate pairs of eigenvalues with real part not equal */
/*             to zero will appear consecutively with the eigenvalue */
/*             having the positive imaginary part first. Note that, due */
/*             to roundoff errors, these numbers may differ from the */
/*             eigenvalues computed by MB03XD. */

/*     US      (output) DOUBLE PRECISION array, dimension (LDUS,MM) */
/*             On exit, if STAB = 'S' or STAB = 'B', the leading 2*N-by-M */
/*             part of this array contains a basis for the stable */
/*             invariant subspace belonging to the selected eigenvalues. */
/*             This basis is orthogonal unless ORTBAL = 'A'. */

/*     LDUS    INTEGER */
/*             The leading dimension of the array US.  LDUS >= 1. */
/*             If STAB = 'S' or STAB = 'B',  LDUS >= 2*N. */

/*     UU      (output) DOUBLE PRECISION array, dimension (LDUU,MM) */
/*             On exit, if STAB = 'U' or STAB = 'B', the leading 2*N-by-M */
/*             part of this array contains a basis for the unstable */
/*             invariant subspace belonging to the selected eigenvalues. */
/*             This basis is orthogonal unless ORTBAL = 'A'. */

/*     LDUU    INTEGER */
/*             The leading dimension of the array UU.  LDUU >= 1. */
/*             If STAB = 'U' or STAB = 'B',  LDUU >= 2*N. */

/*     Workspace */

/*     LWORK   LOGICAL array, dimension (2*N) */
/*             This array is only referenced if WHICH = 'A' and */
/*             METH = 'L'. */

/*     IWORK   INTEGER array, dimension (2*N), */
/*             This array is only referenced if WHICH = 'A' and */
/*             METH = 'L'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -35,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If WHICH = 'S' or METH = 'S': */
/*                LDWORK >= MAX( 1, 4*M*M + MAX( 8*M, 4*N ) ). */
/*             If WHICH = 'A' and METH = 'L' and */
/*                ( STAB = 'U' or STAB = 'S' ): */
/*                LDWORK >= MAX( 1, 2*N*N + 2*N, 8*N ). */
/*             If WHICH = 'A' and METH = 'L' and STAB = 'B': */
/*                LDWORK >= 8*N + 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  some of the selected eigenvalues are on or too close */
/*                   to the imaginary axis; */
/*             = 2:  reordering of the product S*T in routine MB03ZA */
/*                   failed because some eigenvalues are too close to */
/*                   separate; */
/*             = 3:  the QR algorithm failed to compute some Schur form */
/*                   in MB03ZA; */
/*             = 4:  reordering of the Hamiltonian Schur form in routine */
/*                   MB03TD failed because some eigenvalues are too close */
/*                   to separate. */

/*     METHOD */

/*     This is an implementation of Algorithm 1 in [1]. */

/*     NUMERICAL ASPECTS */

/*     The method is strongly backward stable for an embedded */
/*     (skew-)Hamiltonian matrix, see [1]. Although good results have */
/*     been reported if the eigenvalues are not too close to the */
/*     imaginary axis, the method is not backward stable for the original */
/*     Hamiltonian matrix itself. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix, J. Comput. Appl. Math., 86, */
/*         pp. 17-43, 1997. */

/*     [2] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHASUB). */

/*     KEYWORDS */

/*     Hamiltonian matrix, invariant subspace. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and check input parameters. */

    /* Parameter adjustments */
    --select;
    --scale;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    v1_dim1 = *ldv1;
    v1_offset = 1 + v1_dim1;
    v1 -= v1_offset;
    v2_dim1 = *ldv2;
    v2_offset = 1 + v2_dim1;
    v2 -= v2_offset;
    --wr;
    --wi;
    us_dim1 = *ldus;
    us_offset = 1 + us_dim1;
    us -= us_offset;
    uu_dim1 = *lduu;
    uu_offset = 1 + uu_dim1;
    uu -= uu_offset;
    --lwork;
    --iwork;
    --dwork;

    /* Function Body */
    lall = lsame_(which, "A", (ftnlen)1, (ftnlen)1);
    if (lall) {
	lext = lsame_(meth, "L", (ftnlen)1, (ftnlen)1);
    } else {
	lext = FALSE_;
    }
    lus = lsame_(stab, "S", (ftnlen)1, (ftnlen)1) || lsame_(stab, "B", (
	    ftnlen)1, (ftnlen)1);
    luu = lsame_(stab, "U", (ftnlen)1, (ftnlen)1) || lsame_(stab, "B", (
	    ftnlen)1, (ftnlen)1);
    lbal = lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1);
    lbef = FALSE_;
    if (lbal) {
	lbef = lsame_(ortbal, "B", (ftnlen)1, (ftnlen)1);
    }

    wrkmin = 1;
    wrkopt = wrkmin;

    *info = 0;

    if (! lall && ! lsame_(which, "S", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (lall && (! lext && ! lsame_(meth, "S", (ftnlen)1, (ftnlen)1))) 
	    {
	*info = -2;
    } else if (! lus && ! luu) {
	*info = -3;
    } else if (! lbal && ! lsame_(balanc, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -4;
    } else if (lbal && (! lbef && ! lsame_(ortbal, "A", (ftnlen)1, (ftnlen)1))
	    ) {
	*info = -5;
    } else {
	if (lall) {
	    *m = *n;
	} else {

/*           Set M to the dimension of the specified invariant subspace. */

	    *m = 0;
	    pair = FALSE_;
	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {
		if (pair) {
		    pair = FALSE_;
		} else {
		    if (k < *n) {
			if (s[k + 1 + k * s_dim1] == 0.) {
			    if (select[k]) {
				++(*m);
			    }
			} else {
			    pair = TRUE_;
			    if (select[k] || select[k + 1]) {
				*m += 2;
			    }
			}
		    } else {
			if (select[*n]) {
			    ++(*m);
			}
		    }
		}
/* L10: */
	    }
	}

/*        Compute workspace requirements. */

	if (! lext) {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *m << 3, i__4 = *n << 2;
	    i__1 = wrkopt, i__2 = (*m << 2) * *m + max(i__3,i__4);
	    wrkopt = max(i__1,i__2);
	} else {
	    if (lus && luu) {
/* Computing MAX */
		i__1 = wrkopt, i__2 = (*n << 3) + 1;
		wrkopt = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = wrkopt, i__2 = (*n << 1) * *n + (*n << 1), i__1 = max(
			i__1,i__2), i__2 = *n << 3;
		wrkopt = max(i__1,i__2);
	    }
	}

	if (*n < 0) {
	    *info = -7;
	} else if (*mm < *m || lext && *mm < *n << 1) {
	    *info = -8;
	} else if (lbal && (*ilo < 1 || *ilo > *n + 1)) {
	    *info = -9;
	} else if (*lds < max(1,*n)) {
	    *info = -12;
	} else if (*ldt < max(1,*n)) {
	    *info = -14;
	} else if (*ldg < 1 || lext && *ldg < *n) {
	    *info = -16;
	} else if (*ldu1 < max(1,*n)) {
	    *info = -18;
	} else if (*ldu2 < max(1,*n)) {
	    *info = -20;
	} else if (*ldv1 < max(1,*n)) {
	    *info = -22;
	} else if (*ldv2 < max(1,*n)) {
	    *info = -24;
	} else if (*ldus < 1 || lus && *ldus < *n << 1) {
	    *info = -29;
	} else if (*lduu < 1 || luu && *lduu < *n << 1) {
	    *info = -31;
	} else if (*ldwork < wrkmin) {
	    *info = -35;
	    dwork[1] = (doublereal) wrkmin;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03ZD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*m,*n) == 0) {
	dwork[1] = 1.;
	return 0;
    }
    wrkopt = wrkmin;

    if (! lext) {

/*        Workspace requirements: 4*M*M + MAX( 8*M, 4*N ). */

	pw = 1;
	pdw = pw + (*m << 2) * *m;
	i__1 = *m << 1;
	i__2 = *ldwork - pdw + 1;
	mb03za_("No Update", "Update", "Update", "Init", which, &select[1], n,
		 &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], ldg, &u1[
		u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[v1_offset], ldv1, 
		&v2[v2_offset], ldv2, &dwork[pw], &i__1, &wr[1], &wi[1], m, &
		dwork[pdw], &i__2, &ierr, (ftnlen)9, (ftnlen)6, (ftnlen)6, (
		ftnlen)4, (ftnlen)1);
	if (ierr != 0) {
	    goto L250;
	}

	pdw = pw + (*m << 1) * *m;
	i__1 = *m << 1;
	i__2 = *ldwork - pdw + 1;
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw], &
		i__1, &v1[v1_offset], ldv1, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	wrkopt = max(i__1,i__2);

	if (lus) {
	    dlacpy_("All", n, m, &v1[v1_offset], ldv1, &us[us_offset], ldus, (
		    ftnlen)3);
	}
	if (luu) {
	    dlacpy_("All", n, m, &v1[v1_offset], ldv1, &uu[uu_offset], lduu, (
		    ftnlen)3);
	}

	i__1 = *m << 1;
	i__2 = *ldwork - pdw + 1;
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw + *
		m], &i__1, &u1[u1_offset], ldu1, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

	if (lus) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &u1[j * u1_dim1 + 1], &c__1, &us[j * 
			us_dim1 + 1], &c__1);
/* L20: */
	    }
	}
	if (luu) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &u1[j * u1_dim1 + 1], &c__1, &uu[j * 
			uu_dim1 + 1], &c__1);
/* L30: */
	    }
	}

	i__1 = *m << 1;
	i__2 = *ldwork - pdw + 1;
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b34, &dwork[pw], &
		i__1, &v2[v2_offset], ldv2, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

	if (lus) {
	    dlacpy_("All", n, m, &v2[v2_offset], ldv2, &us[*n + 1 + us_dim1], 
		    ldus, (ftnlen)3);
	}
	if (luu) {
	    dlacpy_("All", n, m, &v2[v2_offset], ldv2, &uu[*n + 1 + uu_dim1], 
		    lduu, (ftnlen)3);
	}

	i__1 = *m << 1;
	i__2 = *ldwork - pdw + 1;
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw + *
		m], &i__1, &u2[u2_offset], ldu2, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

	if (lus) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &u2[j * u2_dim1 + 1], &c__1, &us[*n + 1 + j 
			* us_dim1], &c__1);
/* L40: */
	    }
	}
	if (luu) {
	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &u2[j * u2_dim1 + 1], &c__1, &uu[*n + 1 + j 
			* uu_dim1], &c__1);
/* L50: */
	    }
	}

/*        Orthonormalize obtained bases and apply inverse balancing */
/*        transformation. */

	if (lbal && lbef) {
	    if (lus) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	    if (luu) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	}

	if (lus) {
	    i__1 = *n << 1;
	    i__2 = *ldwork - *m;
	    dgeqrf_(&i__1, m, &us[us_offset], ldus, &dwork[1], &dwork[*m + 1],
		     &i__2, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
	    wrkopt = max(i__1,i__2);
	    i__1 = *n << 1;
	    i__2 = *ldwork - *m;
	    dorgqr_(&i__1, m, m, &us[us_offset], ldus, &dwork[1], &dwork[*m + 
		    1], &i__2, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
	    wrkopt = max(i__1,i__2);
	}
	if (luu) {
	    i__1 = *n << 1;
	    i__2 = *ldwork - *m;
	    dgeqrf_(&i__1, m, &uu[uu_offset], lduu, &dwork[1], &dwork[*m + 1],
		     &i__2, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
	    wrkopt = max(i__1,i__2);
	    i__1 = *n << 1;
	    i__2 = *ldwork - *m;
	    dorgqr_(&i__1, m, m, &uu[uu_offset], lduu, &dwork[1], &dwork[*m + 
		    1], &i__2, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
	    wrkopt = max(i__1,i__2);
	}

	if (lbal && ! lbef) {
	    if (lus) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	    if (luu) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	}

    } else {

	i__1 = *n << 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lwork[i__] = TRUE_;
/* L60: */
	}

	if (lus && ! luu) {

/*           Workspace requirements: MAX( 2*N*N + 2*N, 8*N ) */

	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &us[us_offset], 
		    ldus, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
	    if (ierr != 0) {
		goto L250;
	    }

	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &us[*n + 1 + (
		    *n + 1) * us_dim1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);

	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
/* L70: */
	    }
	    pdw = (*n << 1) * *n + 1;

/*           DW <- -[V1;V2]*W11 */

	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &dwork[1], &i__1, (
		    ftnlen)3);
	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &dwork[*n + 1], &i__1, 
		    (ftnlen)3);
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b34, &us[
		    us_offset], ldus, &dwork[1], &i__2, &dwork[pdw], &i__3, &
		    ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	    wrkopt = max(i__1,i__2);

/*           DW2 <- DW2 - U2*W21 */

	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &us[us_offset], ldus, (
		    ftnlen)3);
	    i__1 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + us_dim1], ldus, &us[us_offset], ldus, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &us[j * us_dim1 + 1], &c__1, &dwork[*n + (j 
			- 1 << 1) * *n + 1], &c__1);
/* L80: */
	    }

/*           US11 <- -U1*W21 - DW1 */

	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &us[us_offset], ldus, (
		    ftnlen)3);
	    i__1 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b34, &us[*n + 
		    1 + us_dim1], ldus, &us[us_offset], ldus, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &dwork[(j - 1 << 1) * *n + 1], &c__1, &us[j 
			* us_dim1 + 1], &c__1);
/* L90: */
	    }

/*           US21 <- DW2 */

	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &dwork[*n + 1], &i__1, &us[*n + 1 + us_dim1],
		     ldus, (ftnlen)3);

	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &us[(*n + 1) * us_dim1 
		    + 1], ldus, (ftnlen)3);
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &us[*n + 1 + (*n + 1) *
		     us_dim1], ldus, (ftnlen)3);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &u1[j * u1_dim1 + 1], &c__1, &us[(*n + j) * 
			us_dim1 + 1], &c__1);
/* L100: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &u2[j * u2_dim1 + 1], &c__1, &us[*n + 1 + (*
			n + j) * us_dim1], &c__1);
/* L110: */
	    }

	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &us[(*n + 1) * us_dim1 
		    + 1], ldus, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
	    if (ierr != 0) {
		*info = 4;
		return 0;
	    }
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, n, &us[*n + 1 
		    + (*n + 1) * us_dim1], ldus, &ierr, (ftnlen)7);

	} else if (! lus && luu) {

/*           Workspace requirements: MAX( 2*N*N + 2*N, 8*N ) */

	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &uu[uu_offset], 
		    lduu, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
	    if (ierr != 0) {
		goto L250;
	    }
	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &uu[*n + 1 + (
		    *n + 1) * uu_dim1], lduu, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
/* L120: */
	    }
	    pdw = (*n << 1) * *n + 1;

/*           DW <- -[V1;V2]*W11 */

	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &dwork[1], &i__1, (
		    ftnlen)3);
	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &dwork[*n + 1], &i__1, 
		    (ftnlen)3);
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b34, &uu[
		    uu_offset], lduu, &dwork[1], &i__2, &dwork[pdw], &i__3, &
		    ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
	    wrkopt = max(i__1,i__2);

/*           DW2 <- DW2 - U2*W21 */

	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &uu[uu_offset], lduu, (
		    ftnlen)3);
	    i__1 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b34, &uu[*n + 
		    1 + uu_dim1], lduu, &uu[uu_offset], lduu, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &uu[j * uu_dim1 + 1], &c__1, &dwork[*n + (j 
			- 1 << 1) * *n + 1], &c__1);
/* L130: */
	    }

/*           UU11 <- U1*W21 - DW1 */

	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &uu[uu_offset], lduu, (
		    ftnlen)3);
	    i__1 = *ldwork - pdw + 1;
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + uu_dim1], lduu, &uu[uu_offset], lduu, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b34, &dwork[(j - 1 << 1) * *n + 1], &c__1, &uu[j 
			* uu_dim1 + 1], &c__1);
/* L140: */
	    }

/*           UU21 <- DW2 */

	    i__1 = *n << 1;
	    dlacpy_("All", n, n, &dwork[*n + 1], &i__1, &uu[*n + 1 + uu_dim1],
		     lduu, (ftnlen)3);

	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + (*n + 1) * uu_dim1], lduu, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + (*n + 1) * uu_dim1], lduu, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, (ftnlen)3);
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &uu[*n + 1 + (*n + 1) *
		     uu_dim1], lduu, (ftnlen)3);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &u1[j * u1_dim1 + 1], &c__1, &uu[(*n + j) * 
			uu_dim1 + 1], &c__1);
/* L150: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(n, &c_b26, &u2[j * u2_dim1 + 1], &c__1, &uu[*n + 1 + (*
			n + j) * uu_dim1], &c__1);
/* L160: */
	    }

	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
	    if (ierr != 0) {
		*info = 4;
		return 0;
	    }
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, n, &uu[*n + 1 
		    + (*n + 1) * uu_dim1], lduu, &ierr, (ftnlen)7);
	} else {

/*           Workspace requirements: 8*N */

	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &us[us_offset], 
		    ldus, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
	    if (ierr != 0) {
		goto L250;
	    }
	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &us[*n + 1 + (
		    *n + 1) * us_dim1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
/* L170: */
	    }

/*           UU = [ V1 -V2; U1 -U2 ]*diag(W11,W21) */

	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &uu[uu_offset], lduu, (
		    ftnlen)3);
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &uu[*n + 1 + uu_dim1], 
		    lduu, (ftnlen)3);
	    i__1 = *n << 1;
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b26, &us[
		    us_offset], ldus, &uu[uu_offset], lduu, &dwork[1], ldwork,
		     &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);
	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, (ftnlen)3);
	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &uu[*n + 1 + (*n + 1) *
		     uu_dim1], lduu, (ftnlen)3);
	    i__1 = *n << 1;
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b26, &us[*
		    n + 1 + us_dim1], ldus, &uu[(*n + 1) * uu_dim1 + 1], lduu,
		     &dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)
		    12);
	    i__1 = *n << 1;
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, &i__1, &uu[*n 
		    + 1 + uu_dim1], lduu, &ierr, (ftnlen)7);

	    i__1 = *n << 1;
	    dlacpy_("All", &i__1, n, &uu[uu_offset], lduu, &us[us_offset], 
		    ldus, (ftnlen)3);
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n << 1;
		daxpy_(&i__2, &c_b34, &uu[(*n + j) * uu_dim1 + 1], &c__1, &us[
			j * us_dim1 + 1], &c__1);
/* L180: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n << 1;
		daxpy_(&i__2, &c_b26, &uu[(*n + j) * uu_dim1 + 1], &c__1, &uu[
			j * uu_dim1 + 1], &c__1);
/* L190: */
	    }

/*           V1 <- V1*W12-U1*W22 */
/*           U1 <- V1*W12+U1*W22 */
/*           V2 <- V2*W12-U2*W22 */
/*           U2 <- V2*W12+U2*W22 */

	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = v1[i__ + j * v1_dim1];
		    v1[i__ + j * v1_dim1] = temp - u1[i__ + j * u1_dim1];
		    u1[i__ + j * u1_dim1] = temp + u1[i__ + j * u1_dim1];
/* L200: */
		}
/* L210: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = v2[i__ + j * v2_dim1];
		    v2[i__ + j * v2_dim1] = temp - u2[i__ + j * u2_dim1];
		    u2[i__ + j * u2_dim1] = temp + u2[i__ + j * u2_dim1];
/* L220: */
		}
/* L230: */
	    }

	    i__1 = *n << 1;
	    dlaset_("All", &i__1, n, &c_b272, &c_b26, &us[(*n + 1) * us_dim1 
		    + 1], ldus, (ftnlen)3);
	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &us[(*n + 1) * us_dim1 
		    + 1], ldus, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
	    if (ierr != 0) {
		*info = 4;
		return 0;
	    }
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b26, &u1[
		    u1_offset], ldu1, &us[(*n + 1) * us_dim1 + 1], ldus, &
		    c_b272, &uu[(*n + 1) * uu_dim1 + 1], lduu, (ftnlen)12, (
		    ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u2[
		    u2_offset], ldu2, &us[*n + 1 + (*n + 1) * us_dim1], ldus, 
		    &c_b26, &uu[(*n + 1) * uu_dim1 + 1], lduu, (ftnlen)12, (
		    ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u1[
		    u1_offset], ldu1, &us[*n + 1 + (*n + 1) * us_dim1], ldus, 
		    &c_b272, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, (ftnlen)
		    12, (ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u2[
		    u2_offset], ldu2, &us[(*n + 1) * us_dim1 + 1], ldus, &
		    c_b26, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, (ftnlen)12,
		     (ftnlen)12);
	    dlacpy_("All", n, n, &us[(*n + 1) * us_dim1 + 1], ldus, &u1[
		    u1_offset], ldu1, (ftnlen)3);
	    dlacpy_("All", n, n, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &u2[
		    u2_offset], ldu2, (ftnlen)3);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b26, &v1[
		    v1_offset], ldv1, &u1[u1_offset], ldu1, &c_b272, &us[(*n 
		    + 1) * us_dim1 + 1], ldus, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v2[
		    v2_offset], ldv2, &u2[u2_offset], ldu2, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v1[
		    v1_offset], ldv1, &u2[u2_offset], ldu2, &c_b272, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, (ftnlen)12, (ftnlen)12);
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v2[
		    v2_offset], ldv2, &u1[u1_offset], ldu1, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, (ftnlen)12, (ftnlen)12);
	}

/*        Orthonormalize obtained bases and apply inverse balancing */
/*        transformation. */

	if (lbal && lbef) {
	    if (lus) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	    if (luu) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	}

/*        Workspace requirements: 8*N+1 */

	i__1 = *n << 1;
	for (j = 1; j <= i__1; ++j) {
	    iwork[j] = 0;
/* L240: */
	}
	if (lus) {
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - (*n << 1);
	    dgeqp3_(&i__1, &i__2, &us[us_offset], ldus, &iwork[1], &dwork[1], 
		    &dwork[(*n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
	    wrkopt = max(i__1,i__2);
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - (*n << 1);
	    dorgqr_(&i__1, &i__2, n, &us[us_offset], ldus, &dwork[1], &dwork[(
		    *n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
	    wrkopt = max(i__1,i__2);
	}
	if (luu) {
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - (*n << 1);
	    dgeqp3_(&i__1, &i__2, &uu[uu_offset], lduu, &iwork[1], &dwork[1], 
		    &dwork[(*n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
	    wrkopt = max(i__1,i__2);
	    i__1 = *n << 1;
	    i__2 = *n << 1;
	    i__3 = *ldwork - (*n << 1);
	    dorgqr_(&i__1, &i__2, n, &uu[uu_offset], lduu, &dwork[1], &dwork[(
		    *n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
	    wrkopt = max(i__1,i__2);
	}

	if (lbal && ! lbef) {
	    if (lus) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	    if (luu) {
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
	    }
	}
    }

    dscal_(m, &c_b34, &wr[1], &c__1);
    dwork[1] = (doublereal) wrkopt;

    return 0;
L250:
    if (ierr == 1) {
	*info = 2;
    } else if (ierr == 2 || ierr == 4) {
	*info = 1;
    } else if (ierr == 3) {
	*info = 3;
    }
    return 0;
/* *** Last line of MB03ZD *** */
} /* mb03zd_ */

