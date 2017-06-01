/* IB01PX.f -- translated by f2c (version 20100827).
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

static doublereal c_b11 = 0.;
static doublereal c_b20 = 1.;
static integer c__1 = 1;

/* Subroutine */ int ib01px_(char *job, integer *nobr, integer *n, integer *m,
	 integer *l, doublereal *uf, integer *lduf, doublereal *un, integer *
	ldun, doublereal *ul, integer *ldul, doublereal *pgal, integer *
	ldpgal, doublereal *k, integer *ldk, doublereal *r__, integer *ldr, 
	doublereal *x, doublereal *b, integer *ldb, doublereal *d__, integer *
	ldd, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, d_dim1, d_offset, k_dim1, k_offset, pgal_dim1, 
	    pgal_offset, r_dim1, r_offset, uf_dim1, uf_offset, ul_dim1, 
	    ul_offset, un_dim1, un_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, lp1, np1, npl, rank, ierr;
    static doublereal toll;
    static integer ldun2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb01vd_(char *, char *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    static integer lnobr, mnobr;
    static logical withb, withd;
    static integer mkron, nkron, jwork;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dgelsy_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dtrcon_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer minwrk, maxwrk;


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

/*     To build and solve the least squares problem  T*X = Kv,  and */
/*     estimate the matrices B and D of a linear time-invariant (LTI) */
/*     state space model, using the solution  X,  and the singular */
/*     value decomposition information and other intermediate results, */
/*     provided by other routines. */

/*     The matrix  T  is computed as a sum of Kronecker products, */

/*        T = T + kron(Uf(:,(i-1)*m+1:i*m),N_i),  for i = 1 : s, */

/*     (with  T  initialized by zero), where  Uf  is the triangular */
/*     factor of the QR factorization of the future input part (see */
/*     SLICOT Library routine IB01ND),  N_i  is given by the i-th block */
/*     row of the matrix */

/*            [ Q_11  Q_12  ...  Q_1,s-2  Q_1,s-1  Q_1s ]   [ I_L  0  ] */
/*            [ Q_12  Q_13  ...  Q_1,s-1    Q_1s    0   ]   [         ] */
/*        N = [ Q_13  Q_14  ...    Q_1s      0      0   ] * [         ], */
/*            [  :     :            :        :      :   ]   [         ] */
/*            [ Q_1s   0    ...     0        0      0   ]   [  0  GaL ] */

/*     and where */

/*               [   -L_1|1    ]          [ M_i-1 - L_1|i ] */
/*        Q_11 = [             ],  Q_1i = [               ],  i = 2:s, */
/*               [ I_L - L_2|1 ]          [     -L_2|i    ] */

/*     are  (n+L)-by-L  matrices, and  GaL  is built from the first  n */
/*     relevant singular vectors,  GaL = Un(1:L(s-1),1:n),  computed */
/*     by IB01ND. */

/*     The vector  Kv  is vec(K), with the matrix  K  defined by */

/*        K = [ K_1  K_2  K_3  ...  K_s ], */

/*     where  K_i = K(:,(i-1)*m+1:i*m),  i = 1:s,  is  (n+L)-by-m. */
/*     The given matrices are  Uf,  GaL,  and */

/*            [ L_1|1  ...  L_1|s ] */
/*        L = [                   ],   (n+L)-by-L*s, */
/*            [ L_2|1  ...  L_2|s ] */

/*        M = [ M_1  ...  M_s-1 ],      n-by-L*(s-1),  and */
/*        K,                            (n+L)-by-m*s. */

/*     Matrix  M  is the pseudoinverse of the matrix  GaL,  computed by */
/*     SLICOT Library routine IB01PD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies which of the matrices B and D should be */
/*             computed, as follows: */
/*             = 'B':  compute the matrix B, but not the matrix D; */
/*             = 'D':  compute both matrices B and D. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             Hankel matrices processed by other routines.  NOBR > 1. */

/*     N       (input) INTEGER */
/*             The order of the system.  NOBR > N > 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     UF      (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDUF,M*NOBR ) */
/*             On entry, the leading  M*NOBR-by-M*NOBR  upper triangular */
/*             part of this array must contain the upper triangular */
/*             factor of the QR factorization of the future input part, */
/*             as computed by SLICOT Library routine IB01ND. */
/*             The strict lower triangle need not be set to zero. */
/*             On exit, the leading  M*NOBR-by-M*NOBR  upper triangular */
/*             part of this array is unchanged, and the strict lower */
/*             triangle is set to zero. */

/*     LDUF    INTEGER */
/*             The leading dimension of the array  UF. */
/*             LDUF >= MAX( 1, M*NOBR ). */

/*     UN      (input) DOUBLE PRECISION array, dimension ( LDUN,N ) */
/*             The leading  L*(NOBR-1)-by-N  part of this array must */
/*             contain the matrix  GaL,  i.e., the leading part of the */
/*             first  N  columns of the matrix  Un  of relevant singular */
/*             vectors. */

/*     LDUN    INTEGER */
/*             The leading dimension of the array  UN. */
/*             LDUN >= L*(NOBR-1). */

/*     UL      (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDUL,L*NOBR ) */
/*             On entry, the leading  (N+L)-by-L*NOBR  part of this array */
/*             must contain the given matrix  L. */
/*             On exit, if  M > 0,  the leading  (N+L)-by-L*NOBR  part of */
/*             this array is overwritten by the matrix */
/*             [ Q_11  Q_12  ...  Q_1,s-2  Q_1,s-1  Q_1s ]. */

/*     LDUL    INTEGER */
/*             The leading dimension of the array  UL.  LDUL >= N+L. */

/*     PGAL    (input) DOUBLE PRECISION array, dimension */
/*             ( LDPGAL,L*(NOBR-1) ) */
/*             The leading  N-by-L*(NOBR-1)  part of this array must */
/*             contain the pseudoinverse of the matrix  GaL,  computed by */
/*             SLICOT Library routine IB01PD. */

/*     LDPGAL  INTEGER */
/*             The leading dimension of the array  PGAL.  LDPGAL >= N. */

/*     K       (input) DOUBLE PRECISION array, dimension ( LDK,M*NOBR ) */
/*             The leading  (N+L)-by-M*NOBR  part of this array must */
/*             contain the given matrix  K. */

/*     LDK     INTEGER */
/*             The leading dimension of the array  K.  LDK >= N+L. */

/*     R       (output) DOUBLE PRECISION array, dimension ( LDR,M*(N+L) ) */
/*             The leading  (N+L)*M*NOBR-by-M*(N+L)  part of this array */
/*             contains details of the complete orthogonal factorization */
/*             of the coefficient matrix  T  of the least squares problem */
/*             which is solved for getting the system matrices B and D. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= MAX( 1, (N+L)*M*NOBR ). */

/*     X       (output) DOUBLE PRECISION array, dimension */
/*             ( (N+L)*M*NOBR ) */
/*             The leading  M*(N+L)  elements of this array contain the */
/*             least squares solution of the system  T*X = Kv. */
/*             The remaining elements are used as workspace (to store the */
/*             corresponding part of the vector Kv = vec(K)). */

/*     B       (output) DOUBLE PRECISION array, dimension ( LDB,M ) */
/*             The leading N-by-M part of this array contains the system */
/*             input matrix. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= N. */

/*     D       (output) DOUBLE PRECISION array, dimension ( LDD,M ) */
/*             If  JOB = 'D',  the leading L-by-M part of this array */
/*             contains the system input-output matrix. */
/*             If  JOB = 'B',  this array is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L, if  JOB = 'D'; */
/*             LDD >= 1, if  JOB = 'B'. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  an m-by-n matrix whose estimated */
/*             condition number is less than  1/TOL  is considered to */
/*             be of full rank.  If the user sets  TOL <= 0,  then an */
/*             implicitly computed, default tolerance, defined by */
/*             TOLDEF = m*n*EPS,  is used instead, where  EPS  is the */
/*             relative machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension ( M*(N+L) ) */

/*     DWORK   DOUBLE PRECISION array, dimension ( LDWORK ) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and, if  M > 0,  DWORK(2)  contains the */
/*             reciprocal condition number of the triangular factor of */
/*             the matrix  T. */
/*             On exit, if  INFO = -26,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( (N+L)*(N+L), 4*M*(N+L)+1 ). */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix  T  is computed, evaluating the sum of Kronecker */
/*     products, and then the linear system  T*X = Kv  is solved in a */
/*     least squares sense. The matrices  B  and  D  are then directly */
/*     obtained from the least squares solution. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error */
/*         state-space model identification class of algorithms. */
/*         Int. J. Control, 56, pp. 1187-1210, 1992. */

/*     [2] Van Overschee, P., and De Moor, B. */
/*         N4SID: Two Subspace Algorithms for the Identification */
/*         of Combined Deterministic-Stochastic Systems. */
/*         Automatica, Vol.30, No.1, pp. 75-93, 1994. */

/*     [3] Van Overschee, P. */
/*         Subspace Identification : Theory - Implementation - */
/*         Applications. */
/*         Ph. D. Thesis, Department of Electrical Engineering, */
/*         Katholieke Universiteit Leuven, Belgium, Feb. 1995. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Universiteit Leuven, Feb. 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Universiteit Leuven, Sep. 2001. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    uf_dim1 = *lduf;
    uf_offset = 1 + uf_dim1;
    uf -= uf_offset;
    un_dim1 = *ldun;
    un_offset = 1 + un_dim1;
    un -= un_offset;
    ul_dim1 = *ldul;
    ul_offset = 1 + ul_dim1;
    ul -= ul_offset;
    pgal_dim1 = *ldpgal;
    pgal_offset = 1 + pgal_dim1;
    pgal -= pgal_offset;
    k_dim1 = *ldk;
    k_offset = 1 + k_dim1;
    k -= k_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --iwork;
    --dwork;

    /* Function Body */
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1);
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
    mnobr = *m * *nobr;
    lnobr = *l * *nobr;
    ldun2 = lnobr - *l;
    lp1 = *l + 1;
    np1 = *n + 1;
    npl = *n + *l;
    *iwarn = 0;
    *info = 0;

/*     Check the scalar input parameters. */

    if (! withb) {
	*info = -1;
    } else if (*nobr <= 1) {
	*info = -2;
    } else if (*n >= *nobr || *n <= 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*l <= 0) {
	*info = -5;
    } else if (*lduf < max(1,mnobr)) {
	*info = -7;
    } else if (*ldun < ldun2) {
	*info = -9;
    } else if (*ldul < npl) {
	*info = -11;
    } else if (*ldpgal < *n) {
	*info = -13;
    } else if (*ldk < npl) {
	*info = -15;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = mnobr * npl;
	if (*ldr < max(i__1,i__2)) {
	    *info = -17;
	} else if (*ldb < *n) {
	    *info = -20;
	} else if (*ldd < 1 || withd && *ldd < *l) {
	    *info = -22;
	} else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance. */
/*         NB refers to the optimal block size for the immediately */
/*         following subroutine, as returned by ILAENV.) */

/* Computing MAX */
	    i__1 = npl * npl, i__2 = (*m << 2) * npl + 1;
	    minwrk = max(i__1,i__2);

	    if (*ldwork < minwrk) {
		*info = -26;
		dwork[1] = (doublereal) minwrk;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01PX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Construct the matrix  [ Q_11  Q_12  ...  Q_1,s-1  Q_1s ]  in  UL. */

    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {

	i__2 = npl;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ul[i__ + j * ul_dim1] = -ul[i__ + j * ul_dim1];
/* L10: */
	}

	ul[*n + j + j * ul_dim1] += 1.;
/* L20: */
    }

    i__1 = lnobr;
    for (j = lp1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ul[i__ + j * ul_dim1] = pgal[i__ + (j - *l) * pgal_dim1] - ul[i__ 
		    + j * ul_dim1];
/* L30: */
	}

	i__2 = npl;
	for (i__ = np1; i__ <= i__2; ++i__) {
	    ul[i__ + j * ul_dim1] = -ul[i__ + j * ul_dim1];
/* L40: */
	}

/* L50: */
    }

/*     Compute the coefficient matrix T using Kronecker products. */
/*     Workspace: (N+L)*(N+L). */
/*     In the same loop, vectorize K in X. */

    i__1 = mnobr * npl;
    i__2 = *m * npl;
    dlaset_("Full", &i__1, &i__2, &c_b11, &c_b11, &r__[r_offset], ldr, (
	    ftnlen)4);
    i__1 = mnobr - 1;
    i__2 = mnobr - 1;
    dlaset_("Lower", &i__1, &i__2, &c_b11, &c_b11, &uf[uf_dim1 + 2], lduf, (
	    ftnlen)5);
    jwork = npl * *l + 1;

    i__1 = *nobr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dlacpy_("Full", &npl, l, &ul[((i__ - 1) * *l + 1) * ul_dim1 + 1], 
		ldul, &dwork[1], &npl, (ftnlen)4);
	if (i__ < *nobr) {
	    i__2 = *l * (*nobr - i__);
	    dgemm_("NoTranspose", "NoTranspose", &npl, n, &i__2, &c_b20, &ul[(
		    i__ * *l + 1) * ul_dim1 + 1], ldul, &un[un_offset], ldun, 
		    &c_b11, &dwork[jwork], &npl, (ftnlen)11, (ftnlen)11);
	} else {
	    dlaset_("Full", &npl, n, &c_b11, &c_b11, &dwork[jwork], &npl, (
		    ftnlen)4);
	}
	mb01vd_("NoTranspose", "NoTranspose", &mnobr, m, &npl, &npl, &c_b20, &
		c_b20, &uf[((i__ - 1) * *m + 1) * uf_dim1 + 1], lduf, &dwork[
		1], &npl, &r__[r_offset], ldr, &mkron, &nkron, &ierr, (ftnlen)
		11, (ftnlen)11);
	dlacpy_("Full", &npl, m, &k[((i__ - 1) * *m + 1) * k_dim1 + 1], ldk, &
		x[(i__ - 1) * nkron + 1], &npl, (ftnlen)4);
/* L60: */
    }

/*     Compute the tolerance. */

    toll = *tol;
    if (toll <= 0.) {
	toll = mkron * nkron * dlamch_("Precision", (ftnlen)9);
    }

/*     Solve the least square problem  T*X = vec(K). */
/*     Workspace:  need   4*M*(N+L)+1; */
/*                 prefer 3*M*(N+L)+(M*(N+L)+1)*NB. */

    i__1 = nkron;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
/* L70: */
    }

    dgelsy_(&mkron, &nkron, &c__1, &r__[r_offset], ldr, &x[1], &mkron, &iwork[
	    1], &toll, &rank, &dwork[1], ldwork, &ierr);
    maxwrk = (integer) dwork[1];

/*     Compute the reciprocal of the condition number of the triangular */
/*     factor  R  of  T. */
/*     Workspace: need 3*M*(N+L). */

    dtrcon_("1-norm", "Upper", "NonUnit", &nkron, &r__[r_offset], ldr, &rcond,
	     &dwork[1], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (ftnlen)7);

    if (rank < nkron) {

/*        The least squares problem is rank-deficient. */

	*iwarn = 4;
    }

/*     Construct the matrix  D,  if needed. */

    if (withd) {
	dlacpy_("Full", l, m, &x[1], &npl, &d__[d_offset], ldd, (ftnlen)4);
    }

/*     Construct the matrix  B. */

    dlacpy_("Full", n, m, &x[lp1], &npl, &b[b_offset], ldb, (ftnlen)4);

/*     Return optimal workspace in DWORK(1) and reciprocal condition */
/*     number in  DWORK(2). */

    dwork[1] = (doublereal) max(minwrk,maxwrk);
    dwork[2] = rcond;

    return 0;

/* *** Last line of IB01PX *** */
} /* ib01px_ */

