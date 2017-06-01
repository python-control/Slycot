/* SG03BU.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;
static doublereal c_b19 = -1.;
static doublereal c_b21 = 0.;
static doublereal c_b24 = 1.;
static integer c__4 = 4;
static doublereal c_b77 = .5;
static doublereal c_b78 = 2.;
static integer c__32 = 32;

/* Subroutine */ int sg03bu_(char *trans, integer *n, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *b, integer *ldb, 
	doublereal *scale, doublereal *dwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, m;
    static doublereal s, x, z__, m1[4]	/* was [2][2] */, m2[4]	/* was [2][2] 
	    */, m3[16]	/* was [4][4] */;
    static integer kb, kh, kl;
    static doublereal ui[4]	/* was [2][2] */;
    static integer iw[24];
    static doublereal tm[4]	/* was [2][2] */, rw[32], m3c[16]	/* 
	    was [4][4] */, eps;
    static integer wpt, ypt;
    static doublereal m3ew[4];
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ldws;
    static doublereal uflt;
    static integer info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sg03bw_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    sg03bx_(char *, char *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), drotg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer uiipt;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static doublereal scale1, delta1;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum, smlnum;
    extern /* Subroutine */ int dsyevx_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen, ftnlen);
    static logical notrns;


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

/*     To compute the Cholesky factor U of the matrix X, X = U**T * U or */
/*     X = U * U**T, which is the solution of the generalized d-stable */
/*     discrete-time Lyapunov equation */

/*         T            T                  2    T */
/*        A  * X * A - E  * X * E = - SCALE  * B  * B,                (1) */

/*     or the transposed equation */

/*                 T            T          2        T */
/*        A * X * A  - E * X * E  = - SCALE  * B * B ,                (2) */

/*     respectively, where A, E, B, and U are real N-by-N matrices. The */
/*     Cholesky factor U of the solution is computed without first */
/*     finding X. The pencil A - lambda * E must be in generalized Schur */
/*     form ( A upper quasitriangular, E upper triangular ). Moreover, it */
/*     must be d-stable, i.e. the moduli of its eigenvalues must be less */
/*     than one. B must be an upper triangular matrix with non-negative */
/*     entries on its main diagonal. */

/*     The resulting matrix U is upper triangular. The entries on its */
/*     main diagonal are non-negative. SCALE is an output scale factor */
/*     set to avoid overflow in U. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether equation (1) or equation (2) is to be */
/*             solved: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the quasitriangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the matrix B. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*             array contains the solution matrix U. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in U. */
/*             0 < SCALE <= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (6*N-6) */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the generalized Sylvester equation to be solved in */
/*                   step II (see METHOD) is (nearly) singular to working */
/*                   precision;  perturbed values were used to solve the */
/*                   equation (but the matrices A and E are unchanged); */
/*             = 2:  the generalized Schur form of the pencil */
/*                   A - lambda * E contains a 2-by-2 main diagonal block */
/*                   whose eigenvalues are not a pair of conjugate */
/*                   complex numbers; */
/*             = 3:  the pencil A - lambda * E is not d-stable, i.e. */
/*                   there are eigenvalues outside the open unit circle; */
/*             = 4:  the LAPACK routine DSYEVX utilized to factorize M3 */
/*                   failed to converge. This error is unlikely to occur. */

/*     METHOD */

/*     The method [2] used by the routine is an extension of Hammarling's */
/*     algorithm [1] to generalized Lyapunov equations. */

/*     We present the method for solving equation (1). Equation (2) can */
/*     be treated in a similar fashion. For simplicity, assume SCALE = 1. */

/*     The matrix A is an upper quasitriangular matrix, i.e. it is a */
/*     block triangular matrix with square blocks on the main diagonal */
/*     and the block order at most 2. We use the following partitioning */
/*     for the matrices A, E, B and the solution matrix U */

/*               ( A11   A12 )        ( E11   E12 ) */
/*           A = (           ),   E = (           ), */
/*               (   0   A22 )        (   0   E22 ) */

/*               ( B11   B12 )        ( U11   U12 ) */
/*           B = (           ),   U = (           ).                  (3) */
/*               (   0   B22 )        (   0   U22 ) */

/*     The size of the (1,1)-blocks is 1-by-1 (iff A(2,1) = 0.0) or */
/*     2-by-2. */

/*     We compute U11 and U12**T in three steps. */

/*     Step I: */

/*        From (1) and (3) we get the 1-by-1 or 2-by-2 equation */

/*                T      T                   T      T */
/*             A11  * U11  * U11 * A11  - E11  * U11  * U11 * E11 */

/*                    T */
/*             = - B11  * B11. */

/*        For brevity, details are omitted here. The technique for */
/*        computing U11 is similar to those applied to standard Lyapunov */
/*        equations in Hammarling's algorithm ([1], section 6). */

/*        Furthermore, the auxiliary matrices M1 and M2 defined as */
/*        follows */

/*                               -1      -1 */
/*           M1 = U11 * A11 * E11   * U11 */

/*                         -1      -1 */
/*           M2 = B11 * E11   * U11 */

/*        are computed in a numerically reliable way. */

/*     Step II: */

/*        We solve for U12**T the generalized Sylvester equation */

/*              T      T           T      T */
/*           A22  * U12  * M1 - E22  * U12 */

/*                  T           T      T      T      T */
/*           = - B12  * M2 + E12  * U11  - A12  * U11  * M1. */

/*     Step III: */

/*        One can show that */

/*              T      T                  T      T */
/*           A22  * U22  * U22 * A22 - E22  * U22  * U22 * E22  = */

/*                T              T */
/*           - B22  * B22 - y * y                                     (4) */

/*        holds, where y is defined as follows */

/*                  T      T      T      T */
/*           w = A12  * U11  + A22  * U12 */

/*                    T */
/*           y = ( B12   w ) * M3EV, */

/*        where M3EV is a matrix which fulfils */

/*                ( I-M2*M2**T   -M2*M1**T )              T */
/*           M3 = (                        ) = M3EV * M3EV . */
/*                (  -M1*M2**T  I-M1*M1**T ) */

/*        M3 is positive semidefinite and its rank is equal to the size */
/*        of U11. Therefore, a matrix M3EV can be found by solving the */
/*        symmetric eigenvalue problem for M3 such that y consists of */
/*        either 1 or 2 rows. */

/*        If B22_tilde is the square triangular matrix arising from the */
/*        QR-factorization */

/*               ( B22_tilde )     ( B22  ) */
/*           Q * (           )  =  (      ), */
/*               (     0     )     ( y**T ) */

/*        then */

/*                T              T                T */
/*           - B22  * B22 - y * y   =  - B22_tilde  * B22_tilde. */

/*        Replacing the right hand side in (4) by the term */
/*        - B22_tilde**T * B22_tilde leads to a generalized Lyapunov */
/*        equation of lower dimension compared to (1). */

/*     The solution U of the equation (1) can be obtained by recursive */
/*     application of the steps I to III. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-323, 1982. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The routine requires 2*N**3 flops. Note that we count a single */
/*     floating point arithmetic operation as one flop. */

/*     FURTHER COMMENTS */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if the pencil A - lambda * E has a pair of almost reciprocal */
/*     eigenvalues, then the Lyapunov equation will be ill-conditioned. */
/*     Perturbed values were used to solve the equation. */
/*     A condition estimate can be obtained from the routine SG03AD. */
/*     When setting the error indicator INFO, the routine does not test */
/*     for near instability in the equation but only for exact */
/*     instability. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

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

/*     Decode input parameter. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dwork;

    /* Function Body */
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*lde < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else {
	*info = 0;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SG03BU", &i__1, (ftnlen)6);
	return 0;
    }

    *scale = 1.;

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     Set constants to control overflow. */

    eps = dlamch_("P", (ftnlen)1);
    uflt = dlamch_("S", (ftnlen)1);
    smlnum = uflt / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

/*     Set work space pointers and leading dimension of matrices in */
/*     work space. */

    uiipt = 1;
    wpt = (*n << 1) - 1;
    ypt = (*n << 2) - 3;
    ldws = *n - 1;

    if (notrns) {

/*        Solve equation (1). */

/*        Main Loop. Compute block row U(KL:KH,KL:N). KB denotes the */
/*        number of rows in this block row. */

	kh = 0;
/*        WHILE ( KH .LT. N ) DO */
L20:
	if (kh < *n) {
	    kl = kh + 1;
	    if (kl == *n) {
		kh = *n;
		kb = 1;
	    } else {
		if (a[kl + 1 + kl * a_dim1] == 0.) {
		    kh = kl;
		    kb = 1;
		} else {
		    kh = kl + 1;
		    kb = 2;
		}
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

	    if (kb == 1) {
/* Computing 2nd power */
		d__1 = e[kl + kl * e_dim1];
/* Computing 2nd power */
		d__2 = a[kl + kl * a_dim1];
		delta1 = d__1 * d__1 - d__2 * d__2;
		if (delta1 <= 0.) {
		    *info = 3;
		    return 0;
		}
		delta1 = sqrt(delta1);
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
		if (z__ > delta1) {
		    scale1 = delta1 / z__;
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L40: */
		    }
		}
		ui[0] = b[kl + kl * b_dim1] / delta1;
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
		m2[0] = delta1 / e[kl + kl * e_dim1];
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

		sg03bx_("D", "N", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
		if (info1 != 0) {
		    *info = info1;
		    return 0;
		}
		if (scale1 != 1.) {
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L60: */
		    }
		}
	    }

	    if (kh < *n) {

/*              STEP II: Compute U(KL:KH,KH+1:N) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(KL:KH,KH+1:N) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

		i__1 = *n - kh;
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &b[kl + (kh + 1) * 
			b_dim1], ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
		i__1 = *n - kh;
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b24, &e[kl + (kh + 1) * 
			e_dim1], lde, ui, &c__2, &c_b24, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
		dgemm_("T", "N", &kb, &kb, &kb, &c_b24, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
		i__1 = *n - kh;
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &a[kl + (kh + 1) * 
			a_dim1], lda, tm, &c__2, &c_b24, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

		dlaset_("A", &kb, &kb, &c_b21, &c_b19, tm, &c__2, (ftnlen)1);
		i__1 = *n - kh;
		sg03bw_("N", &i__1, &kb, &a[kh + 1 + (kh + 1) * a_dim1], lda, 
			m1, &c__2, &e[kh + 1 + (kh + 1) * e_dim1], lde, tm, &
			c__2, &dwork[uiipt], &ldws, &scale1, &info1, (ftnlen)
			1);
		if (info1 != 0) {
		    *info = 1;
		}
		if (scale1 != 1.) {
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L80: */
		    }
		    dscal_(&c__4, &scale1, ui, &c__1);
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(KH+1:N,KH+1:N) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary matrices M3 and Y. The factorization */
/*              M3 = M3C * M3C**T is found by solving the symmetric */
/*              eigenvalue problem. */

		i__1 = kb << 1;
		i__2 = kb << 1;
		dlaset_("U", &i__1, &i__2, &c_b21, &c_b24, m3, &c__4, (ftnlen)
			1);
		dsyrk_("U", "N", &kb, &kb, &c_b19, m2, &c__2, &c_b24, m3, &
			c__4, (ftnlen)1, (ftnlen)1);
		dgemm_("N", "T", &kb, &kb, &kb, &c_b19, m2, &c__2, m1, &c__2, 
			&c_b21, &m3[(kb + 1 << 2) - 4], &c__4, (ftnlen)1, (
			ftnlen)1);
		dsyrk_("U", "N", &kb, &kb, &c_b19, m1, &c__2, &c_b24, &m3[kb 
			+ 1 + (kb + 1 << 2) - 5], &c__4, (ftnlen)1, (ftnlen)1)
			;
		i__1 = kb << 1;
		d__1 = uflt * 2.;
		dsyevx_("V", "V", "U", &i__1, m3, &c__4, &c_b77, &c_b78, &
			c__1, &c__4, &d__1, &m, m3ew, m3c, &c__4, rw, &c__32, 
			&iw[4], iw, &info1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
		if (info1 != 0) {
		    *info = 4;
		    return 0;
		}
		i__1 = *n - kh;
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b24, &b[kl + (kh + 1) * 
			b_dim1], ldb, m3c, &c__4, &c_b21, &dwork[ypt], &ldws, 
			(ftnlen)1, (ftnlen)1);
		i__1 = *n - kh;
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b24, &a[kl + (kh + 1) * 
			a_dim1], lda, ui, &c__2, &c_b21, &dwork[wpt], &ldws, (
			ftnlen)1, (ftnlen)1);
		i__1 = *n - kh;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
		    i__3 = i__ + 1, i__4 = *n - kh;
		    i__2 = min(i__3,i__4);
		    dgemv_("T", &i__2, &kb, &c_b24, &dwork[uiipt], &ldws, &a[
			    kh + 1 + (kh + i__) * a_dim1], &c__1, &c_b24, &
			    dwork[wpt + i__ - 1], &ldws, (ftnlen)1);
/* L100: */
		}
		i__1 = *n - kh;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &dwork[wpt], &ldws, 
			&m3c[kb], &c__4, &c_b24, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);

/*              Overwrite B(KH+1:N,KH+1:N) with the triangular matrix */
/*              from the QR-factorization of the (N-KH+KB)-by-(N-KH) */
/*              matrix */

/*                          (  B(KH+1:N,KH+1:N)  ) */
/*                          (                    ) */
/*                          (       Y**T         ) . */

		i__1 = kb;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n - kh;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			x = b[kh + i__ + (kh + i__) * b_dim1];
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
			drotg_(&x, &z__, &c__, &s);
			i__3 = *n - kh - i__ + 1;
			drot_(&i__3, &b[kh + i__ + (kh + i__) * b_dim1], ldb, 
				&dwork[ypt + i__ - 1 + (j - 1) * ldws], &c__1,
				 &c__, &s);
/* L120: */
		    }
/* L140: */
		}

/*              Make main diagonal elements of B(KH+1:N,KH+1:N) positive. */

		i__1 = *n;
		for (i__ = kh + 1; i__ <= i__1; ++i__) {
		    if (b[i__ + i__ * b_dim1] < 0.) {
			i__2 = *n - i__ + 1;
			dscal_(&i__2, &c_b19, &b[i__ + i__ * b_dim1], ldb);
		    }
/* L160: */
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

		i__1 = kh;
		for (j = kl; j <= i__1; ++j) {
		    i__2 = *n - kh;
		    dcopy_(&i__2, &dwork[uiipt + (j - kl) * ldws], &c__1, &b[
			    j + (kh + 1) * b_dim1], ldb);
/* L180: */
		}
	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

	    goto L20;
	}
/*        END WHILE 20 */

    } else {

/*        Solve equation (2). */

/*        Main Loop. Compute block column U(1:KH,KL:KH). KB denotes the */
/*        number of columns in this block column. */

	kl = *n + 1;
/*        WHILE ( KL .GT. 1 ) DO */
L200:
	if (kl > 1) {
	    kh = kl - 1;
	    if (kh == 1) {
		kl = 1;
		kb = 1;
	    } else {
		if (a[kh + (kh - 1) * a_dim1] == 0.) {
		    kl = kh;
		    kb = 1;
		} else {
		    kl = kh - 1;
		    kb = 2;
		}
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

	    if (kb == 1) {
/* Computing 2nd power */
		d__1 = e[kl + kl * e_dim1];
/* Computing 2nd power */
		d__2 = a[kl + kl * a_dim1];
		delta1 = d__1 * d__1 - d__2 * d__2;
		if (delta1 <= 0.) {
		    *info = 3;
		    return 0;
		}
		delta1 = sqrt(delta1);
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
		if (z__ > delta1) {
		    scale1 = delta1 / z__;
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L220: */
		    }
		}
		ui[0] = b[kl + kl * b_dim1] / delta1;
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
		m2[0] = delta1 / e[kl + kl * e_dim1];
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

		sg03bx_("D", "T", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
		if (info1 != 0) {
		    *info = info1;
		    return 0;
		}
		if (scale1 != 1.) {
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L240: */
		    }
		}
	    }

	    if (kl > 1) {

/*              STEP II: Compute U(1:KL-1,KL:KH) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(1:KL-1,KL:KH) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

		i__1 = kl - 1;
		dgemm_("N", "T", &i__1, &kb, &kb, &c_b19, &b[kl * b_dim1 + 1],
			 ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
		i__1 = kl - 1;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &e[kl * e_dim1 + 1],
			 lde, ui, &c__2, &c_b24, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
		dgemm_("N", "T", &kb, &kb, &kb, &c_b24, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
		i__1 = kl - 1;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b19, &a[kl * a_dim1 + 1],
			 lda, tm, &c__2, &c_b24, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

		dlaset_("A", &kb, &kb, &c_b21, &c_b19, tm, &c__2, (ftnlen)1);
		i__1 = kl - 1;
		sg03bw_("T", &i__1, &kb, &a[a_offset], lda, m1, &c__2, &e[
			e_offset], lde, tm, &c__2, &dwork[uiipt], &ldws, &
			scale1, &info1, (ftnlen)1);
		if (info1 != 0) {
		    *info = 1;
		}
		if (scale1 != 1.) {
		    *scale = scale1 * *scale;
		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
/* L260: */
		    }
		    dscal_(&c__4, &scale1, ui, &c__1);
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(1:KL-1,1:KL-1) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary matrices M3 and Y. The factorization */
/*              M3 = M3C * M3C**T is found by solving the symmetric */
/*              eigenvalue problem. */

		i__1 = kb << 1;
		i__2 = kb << 1;
		dlaset_("U", &i__1, &i__2, &c_b21, &c_b24, m3, &c__4, (ftnlen)
			1);
		dsyrk_("U", "T", &kb, &kb, &c_b19, m2, &c__2, &c_b24, m3, &
			c__4, (ftnlen)1, (ftnlen)1);
		dgemm_("T", "N", &kb, &kb, &kb, &c_b19, m2, &c__2, m1, &c__2, 
			&c_b21, &m3[(kb + 1 << 2) - 4], &c__4, (ftnlen)1, (
			ftnlen)1);
		dsyrk_("U", "T", &kb, &kb, &c_b19, m1, &c__2, &c_b24, &m3[kb 
			+ 1 + (kb + 1 << 2) - 5], &c__4, (ftnlen)1, (ftnlen)1)
			;
		i__1 = kb << 1;
		d__1 = uflt * 2.;
		dsyevx_("V", "V", "U", &i__1, m3, &c__4, &c_b77, &c_b78, &
			c__1, &c__4, &d__1, &m, m3ew, m3c, &c__4, rw, &c__32, 
			&iw[4], iw, &info1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
		if (info1 != 0) {
		    *info = 4;
		    return 0;
		}
		i__1 = kl - 1;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &b[kl * b_dim1 + 1],
			 ldb, m3c, &c__4, &c_b21, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);
		i__1 = kl - 1;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &a[kl * a_dim1 + 1],
			 lda, ui, &c__2, &c_b21, &dwork[wpt], &ldws, (ftnlen)
			1, (ftnlen)1);
		i__1 = kl - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
		    i__3 = kl - i__ + 1, i__4 = kl - 1;
		    i__2 = min(i__3,i__4);
/* Computing MAX */
		    i__5 = uiipt, i__6 = uiipt + i__ - 2;
/* Computing MAX */
		    i__7 = i__ - 1;
		    dgemv_("T", &i__2, &kb, &c_b24, &dwork[max(i__5,i__6)], &
			    ldws, &a[i__ + max(i__7,1) * a_dim1], lda, &c_b24,
			     &dwork[wpt + i__ - 1], &ldws, (ftnlen)1);
/* L280: */
		}
		i__1 = kl - 1;
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &dwork[wpt], &ldws, 
			&m3c[kb], &c__4, &c_b24, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);

/*              Overwrite B(1:KL-1,1:KL-1) with the triangular matrix */
/*              from the RQ-factorization of the (KL-1)-by-KH matrix */

/*                          (                        ) */
/*                          (  B(1:KL-1,1:KL-1)   Y  ) */
/*                          (                        ). */

		i__1 = kb;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = kl - 1; i__ >= 1; --i__) {
			x = b[i__ + i__ * b_dim1];
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
			drotg_(&x, &z__, &c__, &s);
			drot_(&i__, &b[i__ * b_dim1 + 1], &c__1, &dwork[ypt + 
				(j - 1) * ldws], &c__1, &c__, &s);
/* L300: */
		    }
/* L320: */
		}

/*              Make main diagonal elements of B(1:KL-1,1:KL-1) positive. */

		i__1 = kl - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    if (b[i__ + i__ * b_dim1] < 0.) {
			dscal_(&i__, &c_b19, &b[i__ * b_dim1 + 1], &c__1);
		    }
/* L340: */
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

		i__1 = kl - 1;
		dlacpy_("A", &i__1, &kb, &dwork[uiipt], &ldws, &b[kl * b_dim1 
			+ 1], ldb, (ftnlen)1);

	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

	    goto L200;
	}
/*        END WHILE 200 */

    }

    return 0;
/* *** Last line of SG03BU *** */
} /* sg03bu_ */

