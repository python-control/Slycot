/* AB01MD.f -- translated by f2c (version 20100827).
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
static doublereal c_b10 = 0.;
static doublereal c_b18 = 1.;
static integer c__0 = 0;

/* Subroutine */ int ab01md_(char *jobz, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ncont, doublereal *z__, integer *ldz, 
	doublereal *tau, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal h__;
    static integer j;
    static doublereal b1, nblk[1];
    static integer itau;
    extern /* Subroutine */ int mb01pd_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), dlarf_(char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
    static logical ljobf, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm, bnorm;
    static logical ljobz;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlarfg_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal fanorm, fbnorm, thresh;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal wrkopt;


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

/*     To find a controllable realization for the linear time-invariant */
/*     single-input system */

/*             dX/dt = A * X + B * U, */

/*     where A is an N-by-N matrix and B is an N element vector which */
/*     are reduced by this routine to orthogonal canonical form using */
/*     (and optionally accumulating) orthogonal similarity */
/*     transformations. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal similarity transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form Z and do not store the orthogonal */
/*                     transformations; */
/*             = 'F':  Do not form Z, but store the orthogonal */
/*                     transformations in the factored form; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NCONT-by-NCONT upper Hessenberg */
/*             part of this array contains the canonical form of the */
/*             state dynamics matrix, given by Z' * A * Z, of a */
/*             controllable realization for the original system. The */
/*             elements below the first subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, the original input/state vector B. */
/*             On exit, the leading NCONT elements of this array contain */
/*             canonical form of the input/state vector, given by Z' * B, */
/*             with all elements but B(1) set to zero. */

/*     NCONT   (output) INTEGER */
/*             The order of the controllable state-space representation. */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If JOBZ = 'I', then the leading N-by-N part of this array */
/*             contains the matrix of accumulated orthogonal similarity */
/*             transformations which reduces the given system to */
/*             orthogonal canonical form. */
/*             If JOBZ = 'F', the elements below the diagonal, with the */
/*             array TAU, represent the orthogonal transformation matrix */
/*             as a product of elementary reflectors. The transformation */
/*             matrix can then be obtained by calling the LAPACK Library */
/*             routine DORGQR. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'I' or */
/*             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The elements of TAU contain the scalar factors of the */
/*             elementary reflectors used in the reduction of B and A. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the */
/*             controllability of (A,B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as an absolute tolerance; */
/*             elements with absolute value less than TOL are considered */
/*             neglijible. If the user sets TOL <= 0, then an implicitly */
/*             computed, default tolerance, defined by */
/*             TOLDEF = N*EPS*MAX( NORM(A), NORM(B) ) is used instead, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             routine DLAMCH). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. LDWORK >= MAX(1,N). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Householder matrix which reduces all but the first element */
/*     of vector B to zero is found and this orthogonal similarity */
/*     transformation is applied to the matrix A. The resulting A is then */
/*     reduced to upper Hessenberg form by a sequence of Householder */
/*     transformations. Finally, the order of the controllable state- */
/*     space representation (NCONT) is determined by finding the position */
/*     of the first sub-diagonal element of A which is below an */
/*     appropriate zero threshold, either TOL or TOLDEF (see parameter */
/*     TOL); if NORM(B) is smaller than this threshold, NCONT is set to */
/*     zero, and no computations for reducing the system to orthogonal */
/*     canonical form are performed. */

/*     REFERENCES */

/*     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D. */
/*         Orthogonal Invariants and Canonical Forms for Linear */
/*         Controllable Systems. */
/*         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981. */

/*     [2] Hammarling, S.J. */
/*         Notes on the use of orthogonal similarity transformations in */
/*         control. */
/*         NPL Report DITC 8/82, August 1982. */

/*     [3] Paige, C.C */
/*         Properties of numerical algorithms related to computing */
/*         controllability. */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1996. */
/*     Supersedes Release 2.0 routine AB01AD by T.W.C. Williams, */
/*     Kingston Polytechnic, United Kingdom, October 1982. */

/*     REVISIONS */

/*     V. Sima, February 16, 1998, October 19, 2001, February 2, 2005. */

/*     KEYWORDS */

/*     Controllability, minimal realization, orthogonal canonical form, */
/*     orthogonal transformation. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobf = lsame_(jobz, "F", (ftnlen)1, (ftnlen)1);
    ljobi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
    ljobz = ljobf || ljobi;

/*     Test the input scalar arguments. */

    if (! ljobz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (! ljobz && *ldz < 1 || ljobz && *ldz < max(1,*n)) {
	*info = -8;
    } else if (*ldwork < max(1,*n)) {
	*info = -12;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB01MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *ncont = 0;
    dwork[1] = 1.;
    if (*n == 0) {
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    wrkopt = 1.;

/*     Calculate the absolute norms of A and B (used for scaling). */

    anorm = dlange_("M", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
    bnorm = dlange_("M", n, &c__1, &b[1], n, &dwork[1], (ftnlen)1);

/*     Return if matrix B is zero. */

    if (bnorm == 0.) {
	if (ljobf) {
	    dlaset_("F", n, n, &c_b10, &c_b10, &z__[z_offset], ldz, (ftnlen)1)
		    ;
	    dlaset_("F", n, &c__1, &c_b10, &c_b10, &tau[1], n, (ftnlen)1);
	} else if (ljobi) {
	    dlaset_("F", n, n, &c_b10, &c_b18, &z__[z_offset], ldz, (ftnlen)1)
		    ;
	}
	return 0;
    }

/*     Scale (if needed) the matrices A and B. */

    mb01pd_("S", "G", n, n, &c__0, &c__0, &anorm, &c__0, nblk, &a[a_offset], 
	    lda, info, (ftnlen)1, (ftnlen)1);
    mb01pd_("S", "G", n, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[1], n, 
	    info, (ftnlen)1, (ftnlen)1);

/*     Calculate the Frobenius norm of A and the 1-norm of B (used for */
/*     controlability test). */

    fanorm = dlange_("F", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
    fbnorm = dlange_("1", n, &c__1, &b[1], n, &dwork[1], (ftnlen)1);

    toldef = *tol;
    if (toldef <= 0.) {

/*        Use the default tolerance in controllability determination. */

	thresh = (doublereal) (*n) * dlamch_("EPSILON", (ftnlen)7);
	toldef = thresh * max(fanorm,fbnorm);
    }

    itau = 1;
    if (fbnorm > toldef) {

/*        B is not negligible compared with A. */

	if (*n > 1) {

/*           Transform B by a Householder matrix Z1: store vector */
/*           describing this temporarily in B and in the local scalar H. */

	    dlarfg_(n, &b[1], &b[2], &c__1, &h__);

	    b1 = b[1];
	    b[1] = 1.;

/*           Form Z1 * A * Z1. */

	    dlarf_("R", n, n, &b[1], &c__1, &h__, &a[a_offset], lda, &dwork[1]
		    , (ftnlen)1);
	    dlarf_("L", n, n, &b[1], &c__1, &h__, &a[a_offset], lda, &dwork[1]
		    , (ftnlen)1);

	    b[1] = b1;
	    tau[1] = h__;
	    ++itau;
	} else {
	    b1 = b[1];
	}

/*        Reduce modified A to upper Hessenberg form by an orthogonal */
/*        similarity transformation with matrix Z2. */
/*        Workspace: need N;  prefer N*NB. */

	dgehrd_(n, &c__1, n, &a[a_offset], lda, &tau[itau], &dwork[1], ldwork,
		 info);
	wrkopt = dwork[1];

	if (ljobz) {

/*           Save the orthogonal transformations used, so that they could */
/*           be accumulated by calling DORGQR routine. */

	    if (*n > 1) {
		i__1 = *n - 1;
		i__2 = *n - 1;
		dlacpy_("F", &i__1, &c__1, &b[2], &i__2, &z__[z_dim1 + 2], 
			ldz, (ftnlen)1);
	    }
	    if (*n > 2) {
		i__1 = *n - 2;
		i__2 = *n - 2;
		dlacpy_("L", &i__1, &i__2, &a[a_dim1 + 3], lda, &z__[(z_dim1 
			<< 1) + 3], ldz, (ftnlen)1);
	    }
	    if (ljobi) {

/*              Form the orthogonal transformation matrix Z = Z1 * Z2. */
/*              Workspace: need N;  prefer N*NB. */

		dorgqr_(n, n, n, &z__[z_offset], ldz, &tau[1], &dwork[1], 
			ldwork, info);
		wrkopt = max(wrkopt,dwork[1]);
	    }
	}

/*        Annihilate the lower part of A and B. */

	if (*n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    dlaset_("L", &i__1, &i__2, &c_b10, &c_b10, &a[a_dim1 + 3], lda, (
		    ftnlen)1);
	}
	if (*n > 1) {
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("F", &i__1, &c__1, &c_b10, &c_b10, &b[2], &i__2, (ftnlen)
		    1);
	}

/*        Find NCONT by checking sizes of the sub-diagonal elements of */
/*        transformed A. */

	if (*tol <= 0.) {
/* Computing MAX */
	    d__1 = fanorm, d__2 = abs(b1);
	    toldef = thresh * max(d__1,d__2);
	}

	j = 1;

/*        WHILE ( J < N and ABS( A(J+1,J) ) > TOLDEF ) DO */

L10:
	if (j < *n) {
	    if ((d__1 = a[j + 1 + j * a_dim1], abs(d__1)) > toldef) {
		++j;
		goto L10;
	    }
	}

/*        END WHILE 10 */

/*        First negligible sub-diagonal element found, if any: set NCONT. */

	*ncont = j;
	if (j < *n) {
	    a[j + 1 + j * a_dim1] = 0.;
	}

/*        Undo scaling of A and B. */

	mb01pd_("U", "H", ncont, ncont, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		a_offset], lda, info, (ftnlen)1, (ftnlen)1);
	mb01pd_("U", "G", &c__1, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[
		1], n, info, (ftnlen)1, (ftnlen)1);
	if (*ncont < *n) {
	    i__1 = *n - *ncont;
	    mb01pd_("U", "G", n, &i__1, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		    (*ncont + 1) * a_dim1 + 1], lda, info, (ftnlen)1, (ftnlen)
		    1);
	}
    } else {

/*        B is negligible compared with A. No computations for reducing */
/*        the system to orthogonal canonical form have been performed, */
/*        except scaling (which is undoed). */

	if (ljobf) {
	    dlaset_("F", n, n, &c_b10, &c_b10, &z__[z_offset], ldz, (ftnlen)1)
		    ;
	    dlaset_("F", n, &c__1, &c_b10, &c_b10, &tau[1], n, (ftnlen)1);
	} else if (ljobi) {
	    dlaset_("F", n, n, &c_b10, &c_b18, &z__[z_offset], ldz, (ftnlen)1)
		    ;
	}
	mb01pd_("U", "G", n, n, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		a_offset], lda, info, (ftnlen)1, (ftnlen)1);
	mb01pd_("U", "G", n, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[1], 
		n, info, (ftnlen)1, (ftnlen)1);
    }

/*     Set optimal workspace dimension. */

    dwork[1] = wrkopt;

    return 0;
/* *** Last line of AB01MD *** */
} /* ab01md_ */

