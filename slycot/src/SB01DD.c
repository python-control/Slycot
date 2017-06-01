/* SB01DD.f -- translated by f2c (version 20100827).
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
static doublereal c_b24 = -1.;
static doublereal c_b26 = 1.;
static integer c__2 = 2;
static doublereal c_b99 = 0.;

/* Subroutine */ int sb01dd_(integer *n, integer *m, integer *indcon, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	nblk, doublereal *wr, doublereal *wi, doublereal *z__, integer *ldz, 
	doublereal *y, integer *count, doublereal *g, integer *ldg, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l;
    static doublereal p, q, r__, s;
    static integer m1, ia, nc, kk, mi, ni, ip, nj, mr, nr, lp1, mp1, np1, mr1,
	     nr1, kmr, rank;
    static doublereal sval[3];
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer iwrk, irmx;
    extern /* Subroutine */ int mb02qd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen),
	     dscal_(integer *, doublereal *, doublereal *, integer *), dlarf_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer indcn1, indcn2;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static integer nblkcr;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer indcrt;
    static logical complx;
    static integer maxwrk;
    static doublereal svlmax;


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

/*     To compute for a controllable matrix pair ( A, B ) a matrix G */
/*     such that the matrix A - B*G has the desired eigenstructure, */
/*     specified by desired eigenvalues and free eigenvector elements. */

/*     The pair ( A, B ) should be given in orthogonal canonical form */
/*     as returned by the SLICOT Library routine AB01ND. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and the number of rows of the */
/*             matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B.  M >= 0. */

/*     INDCON  (input) INTEGER */
/*             The controllability index of the pair ( A, B ). */
/*             0 <= INDCON <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N matrix A in orthogonal canonical form, */
/*             as returned by SLICOT Library routine AB01ND. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the real Schur form of the matrix A - B*G. */
/*             The elements below the real Schur form of A are set to */
/*             zero. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the N-by-M matrix B in orthogonal canonical form, */
/*             as returned by SLICOT Library routine AB01ND. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     NBLK    (input) INTEGER array, dimension (N) */
/*             The leading INDCON elements of this array must contain the */
/*             orders of the diagonal blocks in the orthogonal canonical */
/*             form of A, as returned by SLICOT Library routine AB01ND. */
/*             The values of these elements must satisfy the following */
/*             conditions: */
/*             NBLK(1) >= NBLK(2) >= ... >= NBLK(INDCON), */
/*             NBLK(1) + NBLK(2) + ... + NBLK(INDCON) = N. */

/*     WR      (input) DOUBLE PRECISION array, dimension (N) */
/*     WI      (input) DOUBLE PRECISION array, dimension (N) */
/*             These arrays must contain the real and imaginary parts, */
/*             respectively, of the desired poles of the closed-loop */
/*             system, i.e., the eigenvalues of A - B*G. The poles can be */
/*             unordered, except that complex conjugate pairs of poles */
/*             must appear consecutively. */
/*             The elements of WI for complex eigenvalues are modified */
/*             internally, but restored on exit. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the orthogonal matrix Z generated by SLICOT */
/*             Library routine AB01ND in the reduction of ( A, B ) to */
/*             orthogonal canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal transformation matrix which reduces A - B*G */
/*             to real Schur form. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= max(1,N). */

/*     Y       (input) DOUBLE PRECISION array, dimension (M*N) */
/*             Y contains elements which are used as free parameters */
/*             in the eigenstructure design. The values of these */
/*             parameters are often set by an external optimization */
/*             procedure. */

/*     COUNT   (output) INTEGER */
/*             The actual number of elements in Y used as free */
/*             eigenvector and feedback matrix elements in the */
/*             eigenstructure design. */

/*     G       (output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             The leading M-by-N part of this array contains the */
/*             feedback matrix which assigns the desired eigenstructure */
/*             of A - B*G. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= max(1,M). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where */
/*             EPS  is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(M*N,M*M+2*N+4*M+1). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the pair ( A, B ) is not controllable or the free */
/*                   parameters are not set appropriately. */

/*     METHOD */

/*     The routine implements the method proposed in [1], [2]. */

/*     REFERENCES */

/*     [1] Petkov, P.Hr., Konstantinov, M.M., Gu, D.W. and */
/*         Postlethwaite, I. */
/*         Optimal pole assignment design of linear multi-input systems. */
/*         Report 96-11, Department of Engineering, Leicester University, */
/*         1996. */

/*     [2] Petkov, P.Hr., Christov, N.D. and Konstantinov, M.M. */
/*         A computational algorithm for pole assignment of linear multi */
/*         input systems. IEEE Trans. Automatic Control, vol. AC-31, */
/*         pp. 1044-1047, 1986. */

/*     NUMERICAL ASPECTS */

/*     The method implemented is backward stable. */

/*     FURTHER COMMENTS */

/*     The eigenvalues of the real Schur form matrix As, returned in the */
/*     array A, are very close to the desired eigenvalues WR+WI*i. */
/*     However, the eigenvalues of the closed-loop matrix A - B*G, */
/*     computed by the QR algorithm using the matrices A and B, given on */
/*     entry, may be far from WR+WI*i, although the relative error */
/*        norm( Z'*(A - B*G)*Z - As )/norm( As ) */
/*     is close to machine accuracy. This may happen when the eigenvalue */
/*     problem for the matrix A - B*G is ill-conditioned. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, Technical University of Sofia, Oct. 1998. */
/*     V. Sima, Katholieke Universiteit Leuven, Jan. 1999, SLICOT Library */
/*     version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Closed loop spectrum, closed loop systems, eigenvalue assignment, */
/*     orthogonal canonical form, orthogonal transformation, pole */
/*     placement, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --nblk;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    nr = 0;
/* Computing MAX */
    i__1 = *m * *n, i__2 = *m * *m + (*n << 1) + (*m << 2) + 1;
    iwrk = max(i__1,i__2);
    i__1 = min(*indcon,*n);
    for (i__ = 1; i__ <= i__1; ++i__) {
	nr += nblk[i__];
	if (i__ > 1) {
	    if (nblk[i__ - 1] < nblk[i__]) {
		*info = -8;
	    }
	}
/* L10: */
    }
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*indcon < 0 || *indcon > *n) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (nr != *n) {
	*info = -8;
    } else if (*ldz < max(1,*n)) {
	*info = -12;
    } else if (*ldg < max(1,*m)) {
	*info = -16;
    } else if (*ldwork < iwrk) {
	*info = -20;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SB01DD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*m,*n);
    if (min(i__1,*indcon) == 0) {
	*count = 0;
	dwork[1] = 1.;
	return 0;
    }

    maxwrk = iwrk;
    toldef = *tol;
    if (toldef <= 0.) {

/*        Use the default tolerance, based on machine precision. */

	toldef = (doublereal) (*n * *n) * dlamch_("EPSILON", (ftnlen)7);
    }

    irmx = (*n << 1) + 1;
    iwrk = irmx + *m * *m;
    m1 = nblk[1];
    *count = 1;
    indcrt = *indcon;
    nblkcr = nblk[indcrt];

/*     Compute the Frobenius norm of [ B  A ] (used for rank estimation), */
/*     taking into account the structure. */

    nr = m1;
    nc = 1;
    svlmax = dlange_("Frobenius", &m1, m, &b[b_offset], ldb, &dwork[1], (
	    ftnlen)9);

    i__1 = indcrt - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nr += nblk[i__ + 1];
	d__1 = dlange_("Frobenius", &nr, &nblk[i__], &a[nc * a_dim1 + 1], lda,
		 &dwork[1], (ftnlen)9);
	svlmax = dlapy2_(&svlmax, &d__1);
	nc += nblk[i__];
/* L20: */
    }

    d__1 = dlange_("Frobenius", n, &nblkcr, &a[nc * a_dim1 + 1], lda, &dwork[
	    1], (ftnlen)9);
    svlmax = dlapy2_(&svlmax, &d__1);
    l = 1;
    mr = nblkcr;
    nr = *n - mr + 1;
L30:
/*     WHILE( INDCRT.GT.1 )LOOP */
    if (indcrt > 1) {

/*        Assign next eigenvalue/eigenvector. */

	lp1 = l + m1;
	indcn1 = indcrt - 1;
	mr1 = nblk[indcn1];
	nr1 = nr - mr1;
	complx = wi[l] != 0.;
	dcopy_(&mr, &y[*count], &c__1, &dwork[nr], &c__1);
	*count += mr;
	nc = 1;
	if (complx) {
	    dcopy_(&mr, &y[*count], &c__1, &dwork[*n + nr], &c__1);
	    *count += mr;
	    wi[l + 1] = wi[l] * wi[l + 1];
	    nc = 2;
	}

/*        Compute and transform eiegenvector. */

	i__1 = indcrt;
	for (ip = 1; ip <= i__1; ++ip) {
	    if (ip != indcrt) {
		dlacpy_("Full", &mr, &mr1, &a[nr + nr1 * a_dim1], lda, &dwork[
			irmx], m, (ftnlen)4);
		if (ip == 1) {
		    mp1 = mr;
		    np1 = nr + mp1;
		} else {
		    mp1 = mr + 1;
		    np1 = nr + mp1;
		    s = dasum_(&mp1, &dwork[nr], &c__1);
		    if (complx) {
			s += dasum_(&mp1, &dwork[*n + nr], &c__1);
		    }
		    if (s != 0.) {

/*                    Scale eigenvector elements. */

			d__1 = 1. / s;
			dscal_(&mp1, &d__1, &dwork[nr], &c__1);
			if (complx) {
			    d__1 = 1. / s;
			    dscal_(&mp1, &d__1, &dwork[*n + nr], &c__1);
			    if (np1 <= *n) {
				dwork[*n + np1] /= s;
			    }
			}
		    }
		}

/*              Compute the right-hand side of the eigenvector equations. */

		dcopy_(&mr, &dwork[nr], &c__1, &dwork[nr1], &c__1);
		dscal_(&mr, &wr[l], &dwork[nr1], &c__1);
		dgemv_("No transpose", &mr, &mp1, &c_b24, &a[nr + nr * a_dim1]
			, lda, &dwork[nr], &c__1, &c_b26, &dwork[nr1], &c__1, 
			(ftnlen)12);
		if (complx) {
		    daxpy_(&mr, &wi[l + 1], &dwork[*n + nr], &c__1, &dwork[
			    nr1], &c__1);
		    dcopy_(&mr, &dwork[nr], &c__1, &dwork[*n + nr1], &c__1);
		    daxpy_(&mr, &wr[l + 1], &dwork[*n + nr], &c__1, &dwork[*n 
			    + nr1], &c__1);
		    dgemv_("No transpose", &mr, &mp1, &c_b24, &a[nr + nr * 
			    a_dim1], lda, &dwork[*n + nr], &c__1, &c_b26, &
			    dwork[*n + nr1], &c__1, (ftnlen)12);
		    if (np1 <= *n) {
			d__1 = -dwork[*n + np1];
			daxpy_(&mr, &d__1, &a[nr + np1 * a_dim1], &c__1, &
				dwork[*n + nr1], &c__1);
		    }
		}

/*              Solve linear equations for eigenvector elements. */

		i__2 = *ldwork - iwrk + 1;
		mb02qd_("FreeElements", "NoPermuting", &mr, &mr1, &nc, &
			toldef, &svlmax, &dwork[irmx], m, &dwork[nr1], n, &y[*
			count], &iwork[1], &rank, sval, &dwork[iwrk], &i__2, 
			info, (ftnlen)12, (ftnlen)11);
/* Computing MAX */
		i__2 = maxwrk, i__3 = (integer) dwork[iwrk] + iwrk - 1;
		maxwrk = max(i__2,i__3);
		if (rank < mr) {
		    goto L80;
		}

		*count += (mr1 - mr) * nc;
		nj = nr1;
	    } else {
		nj = nr;
	    }
	    ni = nr + mr - 1;
	    if (ip == 1) {
		kmr = mr - 1;
	    } else {
		kmr = mr;
		if (ip == 2) {
		    ni += nblkcr;
		} else {
		    ni = ni + nblk[indcrt - ip + 2] + 1;
		    if (complx) {
/* Computing MIN */
			i__2 = ni + 1;
			ni = min(i__2,*n);
		    }
		}
	    }

	    i__2 = kmr;
	    for (kk = 1; kk <= i__2; ++kk) {
		k = nr + mr - kk;
		if (ip == 1) {
		    k = *n - kk;
		}
		dlartg_(&dwork[k], &dwork[k + 1], &p, &q, &r__);
		dwork[k] = r__;
		dwork[k + 1] = 0.;

/*              Transform  A. */

		i__3 = *n - nj + 1;
		drot_(&i__3, &a[k + nj * a_dim1], lda, &a[k + 1 + nj * a_dim1]
			, lda, &p, &q);
		drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1, &p, &q);

		if (k < lp1) {

/*                 Transform B. */

		    drot_(m, &b[k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb, &p,
			     &q);
		}

/*              Accumulate transformations. */

		drot_(n, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * z_dim1 + 
			1], &c__1, &p, &q);

		if (complx) {
		    drot_(&c__1, &dwork[*n + k], &c__1, &dwork[*n + k + 1], &
			    c__1, &p, &q);
		    ++k;
		    if (k < *n) {
			dlartg_(&dwork[*n + k], &dwork[*n + k + 1], &p, &q, &
				r__);
			dwork[*n + k] = r__;
			dwork[*n + k + 1] = 0.;

/*                    Transform  A. */

			i__3 = *n - nj + 1;
			drot_(&i__3, &a[k + nj * a_dim1], lda, &a[k + 1 + nj *
				 a_dim1], lda, &p, &q);
			drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
				a_dim1 + 1], &c__1, &p, &q);

			if (k <= lp1) {

/*                       Transform B. */

			    drot_(m, &b[k + b_dim1], ldb, &b[k + 1 + b_dim1], 
				    ldb, &p, &q);
			}

/*                    Accumulate transformations. */

			drot_(n, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * 
				z_dim1 + 1], &c__1, &p, &q);

		    }
		}
/* L40: */
	    }

	    if (ip != indcrt) {
		mr = mr1;
		nr = nr1;
		if (ip != indcn1) {
		    indcn2 = indcrt - ip - 1;
		    mr1 = nblk[indcn2];
		    nr1 -= mr1;
		}
	    }
/* L50: */
	}

	if (! complx) {

/*           Find one column of G. */

	    dlacpy_("Full", &m1, m, &b[l + 1 + b_dim1], ldb, &dwork[irmx], m, 
		    (ftnlen)4);
	    dcopy_(&m1, &a[l + 1 + l * a_dim1], &c__1, &g[l * g_dim1 + 1], &
		    c__1);
	} else {

/*           Find two columns of G. */

	    if (lp1 < *n) {
		++lp1;
		k = l + 2;
	    } else {
		k = l + 1;
	    }
	    dlacpy_("Full", &m1, m, &b[k + b_dim1], ldb, &dwork[irmx], m, (
		    ftnlen)4);
	    dlacpy_("Full", &m1, &c__2, &a[k + l * a_dim1], lda, &g[l * 
		    g_dim1 + 1], ldg, (ftnlen)4);
	    if (k == l + 1) {
		g[l * g_dim1 + 1] -= dwork[*n + l + 1] / dwork[l] * wi[l + 1];
		g[(l + 1) * g_dim1 + 1] = g[(l + 1) * g_dim1 + 1] - wr[l + 1] 
			+ dwork[*n + l] / dwork[l] * wi[l + 1];
	    }
	}

	i__1 = *ldwork - iwrk + 1;
	mb02qd_("FreeElements", "NoPermuting", &m1, m, &nc, &toldef, &svlmax, 
		&dwork[irmx], m, &g[l * g_dim1 + 1], ldg, &y[*count], &iwork[
		1], &rank, sval, &dwork[iwrk], &i__1, info, (ftnlen)12, (
		ftnlen)11);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,i__2);
	if (rank < m1) {
	    goto L80;
	}

	*count += (*m - m1) * nc;
	dgemm_("No transpose", "No transpose", &lp1, &nc, m, &c_b24, &b[
		b_offset], ldb, &g[l * g_dim1 + 1], ldg, &c_b26, &a[l * 
		a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
	++l;
	--nblkcr;
	if (nblkcr == 0) {
	    --indcrt;
	    nblkcr = nblk[indcrt];
	}
	if (complx) {
	    wi[l] = -wi[l - 1];
	    ++l;
	    --nblkcr;
	    if (nblkcr == 0) {
		--indcrt;
		if (indcrt > 0) {
		    nblkcr = nblk[indcrt];
		}
	    }
	}
	mr = nblkcr;
	nr = *n - mr + 1;
	goto L30;
    }
/*     END WHILE 30 */

    if (l <= *n) {

/*        Find the remaining columns of G. */

/*        QR decomposition of the free eigenvectors. */

	i__1 = mr - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia = l + i__ - 1;
	    mi = mr - i__ + 1;
	    dcopy_(&mi, &y[*count], &c__1, &dwork[1], &c__1);
	    *count += mi;
	    dlarfg_(&mi, &dwork[1], &dwork[2], &c__1, &r__);
	    dwork[1] = 1.;

/*           Transform A. */

	    dlarf_("Left", &mi, &mr, &dwork[1], &c__1, &r__, &a[ia + l * 
		    a_dim1], lda, &dwork[*n + 1], (ftnlen)4);
	    dlarf_("Right", n, &mi, &dwork[1], &c__1, &r__, &a[ia * a_dim1 + 
		    1], lda, &dwork[*n + 1], (ftnlen)5);

/*           Transform B. */

	    dlarf_("Left", &mi, m, &dwork[1], &c__1, &r__, &b[ia + b_dim1], 
		    ldb, &dwork[*n + 1], (ftnlen)4);

/*           Accumulate transformations. */

	    dlarf_("Right", n, &mi, &dwork[1], &c__1, &r__, &z__[ia * z_dim1 
		    + 1], ldz, &dwork[*n + 1], (ftnlen)5);
/* L60: */
	}

	i__ = 0;
/*        REPEAT */
L70:
	++i__;
	ia = l + i__ - 1;
	if (wi[ia] == 0.) {
	    dcopy_(&mr, &a[ia + l * a_dim1], lda, &g[i__ + l * g_dim1], ldg);
	    i__1 = mr - i__;
	    daxpy_(&i__1, &c_b24, &y[*count], &c__1, &g[i__ + (l + i__) * 
		    g_dim1], ldg);
	    *count = *count + mr - i__;
	    g[i__ + ia * g_dim1] -= wr[ia];
	} else {
	    dlacpy_("Full", &c__2, &mr, &a[ia + l * a_dim1], lda, &g[i__ + l *
		     g_dim1], ldg, (ftnlen)4);
	    i__1 = mr - i__ - 1;
	    daxpy_(&i__1, &c_b24, &y[*count], &c__2, &g[i__ + (l + i__ + 1) * 
		    g_dim1], ldg);
	    i__1 = mr - i__ - 1;
	    daxpy_(&i__1, &c_b24, &y[*count + 1], &c__2, &g[i__ + 1 + (l + 
		    i__ + 1) * g_dim1], ldg);
	    *count += mr - i__ - 1 << 1;
	    g[i__ + ia * g_dim1] -= wr[ia];
	    g[i__ + (ia + 1) * g_dim1] -= wi[ia];
	    g[i__ + 1 + ia * g_dim1] -= wi[ia + 1];
	    g[i__ + 1 + (ia + 1) * g_dim1] -= wr[ia + 1];
	    ++i__;
	}
	if (i__ < mr) {
	    goto L70;
	}
/*        UNTIL I.GE.MR */

	dlacpy_("Full", &mr, m, &b[l + b_dim1], ldb, &dwork[irmx], m, (ftnlen)
		4);
	i__1 = *ldwork - iwrk + 1;
	mb02qd_("FreeElements", "NoPermuting", &mr, m, &mr, &toldef, &svlmax, 
		&dwork[irmx], m, &g[l * g_dim1 + 1], ldg, &y[*count], &iwork[
		1], &rank, sval, &dwork[iwrk], &i__1, info, (ftnlen)12, (
		ftnlen)11);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,i__2);
	if (rank < mr) {
	    goto L80;
	}

	*count += (*m - mr) * mr;
	dgemm_("No transpose", "No transpose", n, &mr, m, &c_b24, &b[b_offset]
		, ldb, &g[l * g_dim1 + 1], ldg, &c_b26, &a[l * a_dim1 + 1], 
		lda, (ftnlen)12, (ftnlen)12);
    }

/*     Transform G: */
/*     G := G * Z'. */

    dgemm_("No transpose", "Transpose", m, n, n, &c_b26, &g[g_offset], ldg, &
	    z__[z_offset], ldz, &c_b99, &dwork[1], m, (ftnlen)12, (ftnlen)9);
    dlacpy_("Full", m, n, &dwork[1], m, &g[g_offset], ldg, (ftnlen)4);
    --(*count);

    if (*n > 2) {

/*        Set the elements of A below the Hessenberg part to zero. */

	i__1 = *n - 2;
	i__2 = *n - 2;
	dlaset_("Lower", &i__1, &i__2, &c_b99, &c_b99, &a[a_dim1 + 3], lda, (
		ftnlen)5);
    }
    dwork[1] = (doublereal) maxwrk;
    return 0;

/*     Exit with INFO = 1 if the pair ( A, B ) is not controllable or */
/*     the free parameters are not set appropriately. */

L80:
    *info = 1;
    return 0;
/* *** Last line of SB01DD *** */
} /* sb01dd_ */

