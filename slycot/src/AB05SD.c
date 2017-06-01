/* AB05SD.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static doublereal c_b11 = 1.;
static doublereal c_b14 = 0.;
static integer c__1 = 1;

/* Subroutine */ int ab05sd_(char *fbtype, char *jobd, integer *n, integer *m,
	 integer *p, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *f, integer *ldf, 
	doublereal *rcond, integer *iwork, doublereal *dwork, integer *ldwork,
	 integer *info, ftnlen fbtype_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, iw, ldwn, ldwp;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal enorm;
    static logical unitf;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal dummy[1];
    static logical outpf;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen), dgetrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);


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

/*     To construct for a given state space system (A,B,C,D) the closed- */
/*     loop system (Ac,Bc,Cc,Dc) corresponding to the output feedback */
/*     control law */

/*          u = alpha*F*y + v. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FBTYPE  CHARACTER*1 */
/*             Specifies the type of the feedback law as follows: */
/*             = 'I':  Unitary output feedback (F = I); */
/*             = 'O':  General output feedback. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of state variables, i.e. the order of the */
/*             matrix A, the number of rows of B and the number of */
/*             columns of C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of input variables, i.e. the number of columns */
/*             of matrices B and D, and the number of rows of F.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of output variables, i.e. the number of rows of */
/*             matrices C and D, and the number of columns of F.  P >= 0 */
/*             and P = M if FBTYPE = 'I'. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The coefficient alpha in the output feedback law. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state transition matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the state matrix Ac of the closed-loop system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the input matrix Bc of the closed-loop system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the output matrix Cc of the closed-loop system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P) if N > 0. */
/*             LDC >= 1 if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the system direct input/output transmission */
/*             matrix D. */
/*             On exit, if JOBD = 'D', the leading P-by-M part of this */
/*             array contains the direct input/output transmission */
/*             matrix Dc of the closed-loop system. */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P) if JOBD = 'D'. */
/*             LDD >= 1 if JOBD = 'Z'. */

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,P) */
/*             If FBTYPE = 'O', the leading M-by-P part of this array */
/*             must contain the output feedback matrix F. */
/*             If FBTYPE = 'I', then the feedback matrix is assumed to be */
/*             an M x M order identity matrix. */
/*             The array F is not referenced if FBTYPE = 'I' or */
/*             ALPHA = 0. */

/*     LDF     INTEGER */
/*             The leading dimension of array F. */
/*             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0. */
/*             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal condition number of the matrix */
/*             I - alpha*D*F. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,2*P) if JOBD = 'D'. */
/*             LIWORK >= 1 if JOBD = 'Z'. */
/*             IWORK is not referenced if JOBD = 'Z'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= wspace, where */
/*                       wspace = MAX( 1, M, P*P + 4*P ) if JOBD = 'D', */
/*                       wspace = MAX( 1, M ) if JOBD = 'Z'. */
/*             For best performance, LDWORK >= MAX( wspace, N*M, N*P ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix I - alpha*D*F is numerically singular. */

/*     METHOD */

/*     The matrices of the closed-loop system have the expressions: */

/*     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D, */
/*     Cc = E*C,                Dc = E*D, */

/*     where E = (I - alpha*D*F)**-1. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of computations basically depends on the conditioning */
/*     of the matrix I - alpha*D*F.  If RCOND is very small, it is likely */
/*     that the computed results are inaccurate. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Research Establishment, */
/*     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven, */
/*     Belgium, Nov. 1996. */

/*     REVISIONS */

/*     January 14, 1997. */
/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003. */

/*     KEYWORDS */

/*     Multivariable system, state-space model, state-space */
/*     representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External functions .. */
/*     .. External subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the input scalar arguments. */

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
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --iwork;
    --dwork;

    /* Function Body */
    unitf = lsame_(fbtype, "I", (ftnlen)1, (ftnlen)1);
    outpf = lsame_(fbtype, "O", (ftnlen)1, (ftnlen)1);
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
    ldwn = max(1,*n);
    ldwp = max(1,*p);

    *info = 0;

    if (! unitf && ! outpf) {
	*info = -1;
    } else if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0 || unitf && *p != *m) {
	*info = -5;
    } else if (*lda < ldwn) {
	*info = -7;
    } else if (*ldb < ldwn) {
	*info = -9;
    } else if (*n > 0 && *ldc < ldwp || *n == 0 && *ldc < 1) {
	*info = -11;
    } else if (ljobd && *ldd < ldwp || ! ljobd && *ldd < 1) {
	*info = -13;
    } else if (outpf && *alpha != 0. && *ldf < max(1,*m) || (unitf || *alpha 
	    == 0.) && *ldf < 1) {
	*info = -16;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m), i__2 = *p * *p + (*p << 2);
	if (ljobd && *ldwork < max(i__1,i__2) || ! ljobd && *ldwork < max(1,*
		m)) {
	    *info = -20;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB05SD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *rcond = 1.;
/* Computing MAX */
    i__1 = *n, i__2 = min(*m,*p);
    if (max(i__1,i__2) == 0 || *alpha == 0.) {
	return 0;
    }

    if (ljobd) {
	iw = *p * *p + 1;

/*        Compute I - alpha*D*F. */

	if (unitf) {
	    dlacpy_("F", p, p, &d__[d_offset], ldd, &dwork[1], &ldwp, (ftnlen)
		    1);
	    if (*alpha != -1.) {
		d__1 = -(*alpha);
		dlascl_("G", &c__0, &c__0, &c_b11, &d__1, p, p, &dwork[1], &
			ldwp, info, (ftnlen)1);
	    }
	} else {
	    d__1 = -(*alpha);
	    dgemm_("N", "N", p, p, m, &d__1, &d__[d_offset], ldd, &f[f_offset]
		    , ldf, &c_b14, &dwork[1], &ldwp, (ftnlen)1, (ftnlen)1);
	}

	dummy[0] = 1.;
	i__1 = *p + 1;
	daxpy_(p, &c_b11, dummy, &c__0, &dwork[1], &i__1);

/*        Compute Cc = E*C, Dc = E*D, where E = (I - alpha*D*F)**-1. */

	enorm = dlange_("1", p, p, &dwork[1], &ldwp, &dwork[iw], (ftnlen)1);
	dgetrf_(p, p, &dwork[1], &ldwp, &iwork[1], info);
	if (*info > 0) {

/*           Error return. */

	    *rcond = 0.;
	    *info = 1;
	    return 0;
	}
	dgecon_("1", p, &dwork[1], &ldwp, &enorm, rcond, &dwork[iw], &iwork[*
		p + 1], info, (ftnlen)1);
	if (*rcond <= dlamch_("E", (ftnlen)1)) {

/*           Error return. */

	    *info = 1;
	    return 0;
	}

	if (*n > 0) {
	    dgetrs_("N", p, n, &dwork[1], &ldwp, &iwork[1], &c__[c_offset], 
		    ldc, info, (ftnlen)1);
	}
	dgetrs_("N", p, m, &dwork[1], &ldwp, &iwork[1], &d__[d_offset], ldd, 
		info, (ftnlen)1);
    }

    if (*n == 0) {
	return 0;
    }

/*     Compute Ac = A + alpha*B*F*Cc and Bc = B + alpha*B*F*Dc. */

    if (unitf) {
	dgemm_("N", "N", n, n, m, alpha, &b[b_offset], ldb, &c__[c_offset], 
		ldc, &c_b11, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
	if (ljobd) {

	    if (*ldwork < *n * *m) {

/*              Not enough working space for using DGEMM. */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dcopy_(p, &b[i__ + b_dim1], ldb, &dwork[1], &c__1);
		    dgemv_("T", p, p, alpha, &d__[d_offset], ldd, &dwork[1], &
			    c__1, &c_b11, &b[i__ + b_dim1], ldb, (ftnlen)1);
/* L10: */
		}

	    } else {
		dlacpy_("F", n, m, &b[b_offset], ldb, &dwork[1], &ldwn, (
			ftnlen)1);
		dgemm_("N", "N", n, p, m, alpha, &dwork[1], &ldwn, &d__[
			d_offset], ldd, &c_b11, &b[b_offset], ldb, (ftnlen)1, 
			(ftnlen)1);
	    }
	}
    } else {

	if (*ldwork < *n * *p) {

/*           Not enough working space for using DGEMM. */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dgemv_("N", m, p, alpha, &f[f_offset], ldf, &c__[i__ * c_dim1 
			+ 1], &c__1, &c_b14, &dwork[1], &c__1, (ftnlen)1);
		dgemv_("N", n, m, &c_b11, &b[b_offset], ldb, &dwork[1], &c__1,
			 &c_b11, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)1);
/* L20: */
	    }

	    if (ljobd) {

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dgemv_("T", m, p, alpha, &f[f_offset], ldf, &b[i__ + 
			    b_dim1], ldb, &c_b14, &dwork[1], &c__1, (ftnlen)1)
			    ;
		    dgemv_("T", p, m, &c_b11, &d__[d_offset], ldd, &dwork[1], 
			    &c__1, &c_b11, &b[i__ + b_dim1], ldb, (ftnlen)1);
/* L30: */
		}

	    }
	} else {

	    dgemm_("N", "N", n, p, m, alpha, &b[b_offset], ldb, &f[f_offset], 
		    ldf, &c_b14, &dwork[1], &ldwn, (ftnlen)1, (ftnlen)1);
	    dgemm_("N", "N", n, n, p, &c_b11, &dwork[1], &ldwn, &c__[c_offset]
		    , ldc, &c_b11, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
	    if (ljobd) {
		dgemm_("N", "N", n, m, p, &c_b11, &dwork[1], &ldwn, &d__[
			d_offset], ldd, &c_b11, &b[b_offset], ldb, (ftnlen)1, 
			(ftnlen)1);
	    }
	}
    }

    return 0;
/* *** Last line of AB05SD *** */
} /* ab05sd_ */

