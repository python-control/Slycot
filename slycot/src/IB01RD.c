/* IB01RD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b26 = 1.;
static integer c__0 = 0;
static doublereal c_b49 = -1.;
static doublereal c_b150 = .66666666666666663;
static doublereal c_b152 = 0.;

/* Subroutine */ int ib01rd_(char *job, integer *n, integer *m, integer *l, 
	integer *nsmp, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *u, integer *ldu, doublereal *y, integer *ldy, doublereal *
	x0, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, i__3, 
	    i__4, i__5;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, i2, ia, ic, ie, ig, nc, iq, nn, iu, ix, iy, ias, ldr;
    static doublereal dum[1];
    static integer isv, iut, ncp1, ldw1, ldw2, inih, lddw, irem, nrbl, rank;
    static logical ncyc;
    static integer ierr, inir, init, itau, irhs, nobs;
    static doublereal toll;
    static integer nrow;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb04od_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, ftnlen), 
	    mb01td_(integer *, doublereal *, integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical block;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical withd;
    static integer isize;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical first;
    static integer nsmpl, jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer iupnt;
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer iypnt;
    static logical power2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer inigam, icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static integer ncycle;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dgelss_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), dtrcon_(char *, char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical switch__;
    static integer iexpon, iutran, ixinit, minsmp, minwrk;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer maxwrk, minwls;


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

/*     To estimate the initial state of a linear time-invariant (LTI) */
/*     discrete-time system, given the system matrices  (A,B,C,D)  and */
/*     the input and output trajectories of the system. The model */
/*     structure is : */

/*           x(k+1) = Ax(k) + Bu(k),   k >= 0, */
/*           y(k)   = Cx(k) + Du(k), */

/*     where  x(k)  is the  n-dimensional state vector (at time k), */
/*            u(k)  is the  m-dimensional input vector, */
/*            y(k)  is the  l-dimensional output vector, */
/*     and  A, B, C, and D  are real matrices of appropriate dimensions. */
/*     Matrix  A  is assumed to be in a real Schur form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies whether or not the matrix D is zero, as follows: */
/*             = 'Z':  the matrix  D  is zero; */
/*             = 'N':  the matrix  D  is not zero. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples used,  t).  NSMP >= N. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix  A  in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix  B  (corresponding to the real Schur */
/*             form of  A). */
/*             If  N = 0  or  M = 0,  this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= N,  if  N > 0  and  M > 0; */
/*             LDB >= 1,  if  N = 0  or   M = 0. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array must contain the */
/*             system output matrix  C  (corresponding to the real Schur */
/*             form of  A). */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= L. */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading L-by-M part of this array must contain the */
/*             system input-output matrix. */
/*             If  M = 0  or  JOB = 'Z',  this array is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if  M > 0  and  JOB = 'N'; */
/*             LDD >= 1,  if  M = 0  or   JOB = 'Z'. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             If  M > 0,  the leading NSMP-by-M part of this array must */
/*             contain the t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
/*             If M = 0, this array is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,NSMP),  if M > 0; */
/*             LDU >= 1,            if M = 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             t-by-l output-data sequence matrix  Y, */
/*             Y = [y_1 y_2 ... y_l].  Column  j  of  Y  contains the */
/*             NSMP  values of the j-th output component for consecutive */
/*             time increments. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     X0      (output) DOUBLE PRECISION array, dimension (N) */
/*             The estimated initial state of the system,  x(0). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  a matrix whose estimated condition */
/*             number is less than  1/TOL  is considered to be of full */
/*             rank.  If the user sets  TOL <= 0,  then  EPS  is used */
/*             instead, where  EPS  is the relative machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL <= 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK and  DWORK(2)  contains the reciprocal condition */
/*             number of the triangular factor of the QR factorization of */
/*             the matrix  Gamma  (see METHOD). */
/*             On exit, if  INFO = -22,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 2, min( LDW1, LDW2 ) ),  where */
/*             LDW1 = t*L*(N + 1) + 2*N + max( 2*N*N, 4*N ), */
/*             LDW2 =   N*(N + 1) + 2*N + */
/*                      max( q*(N + 1) + 2*N*N + L*N, 4*N ), */
/*                q = N*L. */
/*             For good performance,  LDWORK  should be larger. */
/*             If  LDWORK >= LDW1,  then standard QR factorization of */
/*             the matrix  Gamma  (see METHOD) is used. Otherwise, the */
/*             QR factorization is computed sequentially by performing */
/*             NCYCLE  cycles, each cycle (except possibly the last one) */
/*             processing  s  samples, where  s  is chosen by equating */
/*             LDWORK  to  LDW2,  for  q  replaced by  s*L. */
/*             The computational effort may increase and the accuracy may */
/*             decrease with the decrease of  s.  Recommended value is */
/*             LDRWRK = LDW1,  assuming a large enough cache size, to */
/*             also accommodate  A, B, C, D, U,  and  Y. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     An extension and refinement of the method in [1] is used. */
/*     Specifically, the output y0(k) of the system for zero initial */
/*     state is computed for k = 0, 1, ...,  t-1 using the given model. */
/*     Then the following least squares problem is solved for x(0) */

/*                         (     C     )            (   y(0) - y0(0)   ) */
/*                         (    C*A    )            (   y(1) - y0(1)   ) */
/*        Gamma * x(0)  =  (     :     ) * x(0)  =  (        :         ). */
/*                         (     :     )            (        :         ) */
/*                         ( C*A^(t-1) )            ( y(t-1) - y0(t-1) ) */

/*     The coefficient matrix  Gamma  is evaluated using powers of A with */
/*     exponents 2^k. The QR decomposition of this matrix is computed. */
/*     If its triangular factor  R  is too ill conditioned, then singular */
/*     value decomposition of  R  is used. */

/*     If the coefficient matrix cannot be stored in the workspace (i.e., */
/*     LDWORK < LDW1),  the QR decomposition is computed sequentially. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Varga, A. */
/*         Some Experience with the MOESP Class of Subspace Model */
/*         Identification Methods in Identifying the BO105 Helicopter. */
/*         Report TR R165-94, DLR Oberpfaffenhofen, 1994. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     IBLOCK is a threshold value for switching to a block algorithm */
/*     for  U  (to avoid row by row passing through  U). */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --x0;
    --iwork;
    --dwork;

    /* Function Body */
    withd = lsame_(job, "N", (ftnlen)1, (ftnlen)1);
    *iwarn = 0;
    *info = 0;
    nn = *n * *n;
    minsmp = *n;

    if (! (lsame_(job, "Z", (ftnlen)1, (ftnlen)1) || withd)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*l <= 0) {
	*info = -4;
    } else if (*nsmp < minsmp) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < 1 || *ldb < *n && *m > 0) {
	*info = -9;
    } else if (*ldc < *l) {
	*info = -11;
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
	*info = -13;
    } else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
	*info = -15;
    } else if (*ldy < max(1,*nsmp)) {
	*info = -17;
    } else if (*tol > 1.) {
	*info = -19;
    }

/*     Compute workspace. */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

    nsmpl = *nsmp * *l;
    iq = minsmp * *l;
    ncp1 = *n + 1;
    isize = nsmpl * ncp1;
    ic = nn << 1;
    minwls = minsmp * ncp1;
    itau = ic + *l * *n;
/* Computing MAX */
    i__1 = ic, i__2 = *n << 2;
    ldw1 = isize + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
    i__1 = iq * ncp1 + itau, i__2 = *n << 2;
    ldw2 = minwls + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
    i__1 = min(ldw1,ldw2);
    minwrk = max(i__1,2);
    if (*info == 0 && *ldwork >= minwrk) {
/* Computing MAX */
	i__1 = *n * ilaenv_(&c__1, "DGEQRF", " ", &nsmpl, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1), i__2 = ilaenv_(&c__1, "DORMQR", "LT", &
		nsmpl, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)2);
	maxwrk = isize + (*n << 1) + max(i__1,i__2);
	maxwrk = max(maxwrk,minwrk);
    }

    if (*info == 0 && *ldwork < minwrk) {
	*info = -22;
	dwork[1] = (doublereal) minwrk;
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01RD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 2.;
	dwork[2] = 1.;
	return 0;
    }

/*     Set up the least squares problem, either directly, if enough */
/*     workspace, or sequentially, otherwise. */

    iypnt = 1;
    iupnt = 1;
    inir = 1;
    if (*ldwork >= ldw1) {

/*        Enough workspace for solving the problem directly. */

	ncycle = 1;
	nobs = *nsmp;
	lddw = nsmpl;
	inigam = 1;
    } else {

/*        NCYCLE > 1  cycles are needed for solving the problem */
/*        sequentially, taking  NOBS  samples in each cycle (or the */
/*        remaining samples in the last cycle). */

	jwork = *ldwork - minwls - (*n << 1) - itau;
	lddw = jwork / ncp1;
	nobs = lddw / *l;
	lddw = *l * nobs;
	ncycle = *nsmp / nobs;
	if (*nsmp % nobs != 0) {
	    ++ncycle;
	}
	inih = inir + nn;
	inigam = inih + *n;
    }

    ncyc = ncycle > 1;
    irhs = inigam + lddw * *n;
    ixinit = irhs + lddw;
    ic = ixinit + *n;
    if (ncyc) {
	ia = ic + *l * *n;
	ldr = *n;
	ie = inigam;
    } else {
	inih = irhs;
	ia = ic;
	ldr = lddw;
	ie = ixinit;
    }
    iutran = ia;
    ias = ia + nn;
    itau = ia;
    dum[0] = 0.;

/*     Set block parameters for passing through the array  U. */

    block = *m > 1 && *nsmp * *m >= 16384;
    if (block) {
	nrbl = (*ldwork - iutran + 1) / *m;
	nc = nobs / nrbl;
	if (nobs % nrbl != 0) {
	    ++nc;
	}
	init = (nc - 1) * nrbl;
	block = block && nrbl > 1;
    }

/*     Perform direct of sequential compression of the matrix  Gamma. */

    i__1 = ncycle;
    for (icycle = 1; icycle <= i__1; ++icycle) {
	first = icycle == 1;
	if (! first) {
	    if (icycle == ncycle) {
		nobs = *nsmp - (ncycle - 1) * nobs;
		lddw = *l * nobs;
		if (block) {
		    nc = nobs / nrbl;
		    if (nobs % nrbl != 0) {
			++nc;
		    }
		    init = (nc - 1) * nrbl;
		}
	    }
	}

/*        Compute the extended observability matrix  Gamma. */
/*        Workspace: need   s*L*(N + 1) + 2*N*N + 2*N + a + w, */
/*                   where  s = NOBS, */
/*                          a = 0,   w = 0,          if NCYCLE = 1, */
/*                          a = L*N, w = N*(N + 1),  if NCYCLE > 1; */
/*                   prefer as above, with  s = t,  a = w = 0. */

	jwork = ias + nn;
	iexpon = (integer) (log((doublereal) nobs) / log(2.));
	irem = *l * (nobs - pow_ii(&c__2, &iexpon));
	power2 = irem == 0;
	if (! power2) {
	    ++iexpon;
	}

	if (first) {
	    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[inigam], &lddw, 
		    (ftnlen)4);
	} else {
	    dlacpy_("Full", l, n, &dwork[ic], l, &dwork[inigam], &lddw, (
		    ftnlen)4);
	}
/*                                       p */
/*        Use powers of the matrix  A:  A ,  p = 2**(J-1). */

	dlacpy_("Upper", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)5);
	if (*n > 1) {
	    i__2 = *n - 1;
	    i__3 = *lda + 1;
	    i__4 = *n + 1;
	    dcopy_(&i__2, &a[a_dim1 + 2], &i__3, &dwork[ia + 1], &i__4);
	}
	i2 = *l;
	nrow = 0;

	i__2 = iexpon;
	for (j = 1; j <= i__2; ++j) {
	    ig = inigam;
	    if (j < iexpon || power2) {
		nrow = i2;
	    } else {
		nrow = irem;
	    }

	    dlacpy_("Full", &nrow, n, &dwork[ig], &lddw, &dwork[ig + i2], &
		    lddw, (ftnlen)4);
	    dtrmm_("Right", "Upper", "No Transpose", "Non Unit", &nrow, n, &
		    c_b26, &dwork[ia], n, &dwork[ig + i2], &lddw, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
/*                                                            p */
/*           Compute the contribution of the subdiagonal of  A   to the */
/*           product. */

	    i__3 = *n - 1;
	    for (ix = 1; ix <= i__3; ++ix) {
		daxpy_(&nrow, &dwork[ia + (ix - 1) * *n + ix], &dwork[ig + 
			lddw], &c__1, &dwork[ig + i2], &c__1);
		ig += lddw;
/* L10: */
	    }

	    if (j < iexpon) {
		dlacpy_("Upper", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)
			5);
		i__3 = *n - 1;
		i__4 = *n + 1;
		i__5 = *n + 1;
		dcopy_(&i__3, &dwork[ia + 1], &i__4, &dwork[ias + 1], &i__5);
		mb01td_(n, &dwork[ias], n, &dwork[ia], n, &dwork[jwork], &
			ierr);
		i2 <<= 1;
	    }
/* L20: */
	}

	if (ncyc) {
	    ig = inigam + i2 + nrow - *l;
	    dlacpy_("Full", l, n, &dwork[ig], &lddw, &dwork[ic], l, (ftnlen)4)
		    ;
	    dtrmm_("Right", "Upper", "No Transpose", "Non Unit", l, n, &c_b26,
		     &a[a_offset], lda, &dwork[ic], l, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);

/*           Compute the contribution of the subdiagonal of  A  to the */
/*           product. */

	    i__2 = *n - 1;
	    for (ix = 1; ix <= i__2; ++ix) {
		daxpy_(l, &a[ix + 1 + ix * a_dim1], &dwork[ig + lddw], &c__1, 
			&dwork[ic + (ix - 1) * *l], &c__1);
		ig += lddw;
/* L30: */
	    }

	}

/*        Setup (part of) the right hand side of the least squares */
/*        problem starting from  DWORK(IRHS);  use the estimated output */
/*        trajectory for zero initial state, or for the saved final state */
/*        value of the previous cycle. */
/*        A specialization of SLICOT Library routine TF01ND is used. */
/*        For large input sets  (NSMP*M >= IBLOCK),  chunks of  U  are */
/*        transposed, to reduce the number of row-wise passes. */
/*        Workspace: need   s*L*(N + 1) + N + w; */
/*                   prefer as above, with  s = t,  w = 0. */

	if (first) {
	    dcopy_(n, dum, &c__0, &dwork[ixinit], &c__1);
	}
	dcopy_(n, &dwork[ixinit], &c__1, &x0[1], &c__1);
	iy = irhs;

	i__2 = *l;
	for (j = 1; j <= i__2; ++j) {
	    dcopy_(&nobs, &y[iypnt + j * y_dim1], &c__1, &dwork[iy], l);
	    ++iy;
/* L40: */
	}

	iy = irhs;
	iu = iupnt;
	if (*m > 0) {
	    if (withd) {

		if (block) {
		    switch__ = TRUE_;
		    nrow = nrbl;

		    i__2 = nobs;
		    for (k = 1; k <= i__2; ++k) {
			if ((k - 1) % nrow == 0 && switch__) {
			    iut = iutran;
			    if (k > init) {
				nrow = nobs - init;
				switch__ = FALSE_;
			    }
			    ma02ad_("Full", &nrow, m, &u[iu + u_dim1], ldu, &
				    dwork[iut], m, (ftnlen)4);
			    iu += nrow;
			}
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
			dgemv_("No transpose", l, m, &c_b49, &d__[d_offset], 
				ldd, &dwork[iut], &c__1, &c_b26, &dwork[iy], &
				c__1, (ftnlen)12);
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

			i__3 = *n;
			for (ix = 2; ix <= i__3; ++ix) {
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
/* L50: */
			}

			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &dwork[iut], &c__1, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
			iy += *l;
			iut += *m;
/* L60: */
		    }

		} else {

		    i__2 = nobs;
		    for (k = 1; k <= i__2; ++k) {
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
			dgemv_("No transpose", l, m, &c_b49, &d__[d_offset], 
				ldd, &u[iu + u_dim1], ldu, &c_b26, &dwork[iy],
				 &c__1, (ftnlen)12);
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

			i__3 = *n;
			for (ix = 2; ix <= i__3; ++ix) {
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
/* L70: */
			}

			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &u[iu + u_dim1], ldu, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
			iy += *l;
			++iu;
/* L80: */
		    }

		}

	    } else {

		if (block) {
		    switch__ = TRUE_;
		    nrow = nrbl;

		    i__2 = nobs;
		    for (k = 1; k <= i__2; ++k) {
			if ((k - 1) % nrow == 0 && switch__) {
			    iut = iutran;
			    if (k > init) {
				nrow = nobs - init;
				switch__ = FALSE_;
			    }
			    ma02ad_("Full", &nrow, m, &u[iu + u_dim1], ldu, &
				    dwork[iut], m, (ftnlen)4);
			    iu += nrow;
			}
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

			i__3 = *n;
			for (ix = 2; ix <= i__3; ++ix) {
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
/* L90: */
			}

			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &dwork[iut], &c__1, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
			iy += *l;
			iut += *m;
/* L100: */
		    }

		} else {

		    i__2 = nobs;
		    for (k = 1; k <= i__2; ++k) {
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

			i__3 = *n;
			for (ix = 2; ix <= i__3; ++ix) {
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
/* L110: */
			}

			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &u[iu + u_dim1], ldu, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
			iy += *l;
			++iu;
/* L120: */
		    }

		}

	    }

	} else {

	    i__2 = nobs;
	    for (k = 1; k <= i__2; ++k) {
		dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], ldc, &x0[
			1], &c__1, &c_b26, &dwork[iy], &c__1, (ftnlen)12);
		dtrmv_("Upper", "No transpose", "Non-unit", n, &a[a_offset], 
			lda, &x0[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

		i__3 = *n;
		for (ix = 2; ix <= i__3; ++ix) {
		    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[ixinit + ix - 
			    2];
/* L130: */
		}

		dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
		iy += *l;
/* L140: */
	    }

	}

/*        Compress the data using (sequential) QR factorization. */
/*        Workspace: need   v + 2*N; */
/*                   where  v = s*L*(N + 1) + N + a + w. */

	jwork = itau + *n;
	if (first) {

/*           Compress the first data segment of  Gamma. */
/*           Workspace: need   v + 2*N, */
/*                      prefer v + N + N*NB. */

	    i__2 = *ldwork - jwork + 1;
	    dgeqrf_(&lddw, n, &dwork[inigam], &lddw, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);

/*           Apply the transformation to the right hand side part. */
/*           Workspace: need   v + N + 1, */
/*                      prefer v + N + NB. */

	    i__2 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &lddw, &c__1, n, &dwork[inigam], &
		    lddw, &dwork[itau], &dwork[irhs], &lddw, &dwork[jwork], &
		    i__2, &ierr, (ftnlen)4, (ftnlen)9);

	    if (ncyc) {

/*              Save the triangular factor of  Gamma  and the */
/*              corresponding right hand side. */

		dlacpy_("Upper", n, &ncp1, &dwork[inigam], &lddw, &dwork[inir]
			, &ldr, (ftnlen)5);
	    }
	} else {

/*           Compress the current (but not the first) data segment of */
/*           Gamma. */
/*           Workspace: need   v + N - 1. */

	    mb04od_("Full", n, &c__1, &lddw, &dwork[inir], &ldr, &dwork[
		    inigam], &lddw, &dwork[inih], &ldr, &dwork[irhs], &lddw, &
		    dwork[itau], &dwork[jwork], (ftnlen)4);
	}

	iupnt += nobs;
	iypnt += nobs;
/* L150: */
    }

/*     Estimate the reciprocal condition number of the triangular factor */
/*     of the QR decomposition. */
/*     Workspace: need  u + 3*N, where */
/*                      u = t*L*(N + 1), if NCYCLE = 1; */
/*                      u = w,           if NCYCLE > 1. */

    dtrcon_("1-norm", "Upper", "No Transpose", n, &dwork[inir], &ldr, &rcond, 
	    &dwork[ie], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (ftnlen)12);

    toll = *tol;
    if (toll <= 0.) {
	toll = dlamch_("Precision", (ftnlen)9);
    }
    if (rcond <= pow_dd(&toll, &c_b150)) {
	*iwarn = 4;

/*        The least squares problem is ill-conditioned. */
/*        Use SVD to solve it. */
/*        Workspace: need   u + 6*N; */
/*                   prefer larger. */

	i__1 = *n - 1;
	i__2 = *n - 1;
	dlaset_("Lower", &i__1, &i__2, &c_b152, &c_b152, &dwork[inir + 1], &
		ldr, (ftnlen)5);
	isv = ie;
	jwork = isv + *n;
	i__1 = *ldwork - jwork + 1;
	dgelss_(n, n, &c__1, &dwork[inir], &ldr, &dwork[inih], &ldr, &dwork[
		isv], &toll, &rank, &dwork[jwork], &i__1, &ierr);
	if (ierr > 0) {

/*           Return if SVD algorithm did not converge. */

	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] - jwork + 1;
	maxwrk = max(i__1,i__2);
    } else {

/*        Find the least squares solution using QR decomposition only. */

	dtrsv_("Upper", "No Transpose", "Non Unit", n, &dwork[inir], &ldr, &
		dwork[inih], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
    }

/*     Return the estimated initial state of the system  x0. */

    dcopy_(n, &dwork[inih], &c__1, &x0[1], &c__1);

    dwork[1] = (doublereal) maxwrk;
    dwork[2] = rcond;

    return 0;

/* *** End of IB01RD *** */
} /* ib01rd_ */

