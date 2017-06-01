/* MB04TS.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04ts_(char *trana, char *tranb, integer *n, integer *
	ilo, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *g, integer *ldg, doublereal *q, integer *ldq, doublereal *
	csl, doublereal *csr, doublereal *taul, doublereal *taur, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen trana_len, ftnlen 
	tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal s, nu;
    static logical ltra, ltrb;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal alpha;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen);


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

/*     To compute a symplectic URV (SURV) decomposition of a real */
/*     2N-by-2N matrix H: */

/*             [ op(A)   G   ]        T       [ op(R11)   R12   ]    T */
/*         H = [             ] = U R V  = U * [                 ] * V , */
/*             [  Q    op(B) ]                [   0     op(R22) ] */

/*     where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real */
/*     N-by-N upper triangular matrix, op(R22) is a real N-by-N lower */
/*     Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic */
/*     matrices. Unblocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op( A ) as follows: */
/*             = 'N': op( A ) = A; */
/*             = 'T': op( A ) = A'; */
/*             = 'C': op( A ) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op( B ) as follows: */
/*             = 'N': op( B ) = B; */
/*             = 'T': op( B ) = B'; */
/*             = 'C': op( B ) = B'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     ILO     (input) INTEGER */
/*             It is assumed that op(A) is already upper triangular, */
/*             op(B) is lower triangular and Q is zero in rows and */
/*             columns 1:ILO-1. ILO is normally set by a previous call */
/*             to MB04DD; otherwise it should be set to 1. */
/*             1 <= ILO <= N, if N > 0; ILO=1, if N=0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the triangular matrix R11, and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the Hessenberg matrix R22, and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix G. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix R12. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix Q. */
/*             On exit, the leading N-by-N part of this array contains */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDG >= MAX(1,N). */

/*     CSL     (output) DOUBLE PRECISION array, dimension (2N) */
/*             On exit, the first 2N elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the left-hand side used to compute the SURV */
/*             decomposition. */

/*     CSR     (output) DOUBLE PRECISION array, dimension (2N-2) */
/*             On exit, the first 2N-2 elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the right-hand side used to compute the SURV */
/*             decomposition. */

/*     TAUL    (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, the first N elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied from the left-hand side. */

/*     TAUR    (output) DOUBLE PRECISION array, dimension (N-1) */
/*             On exit, the first N-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied from the right-hand side. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrices U and V are represented as products of symplectic */
/*     reflectors and Givens rotators */

/*     U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) ) */
/*         diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) ) */
/*                              .... */
/*         diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ), */

/*     V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) ) */
/*         diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) ) */
/*                                   .... */
/*         diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ). */

/*     Each HU(i) has the form */

/*           HU(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in */
/*     Q(i+1:n,i), and tau in Q(i,i). */

/*     Each FU(i) has the form */

/*           FU(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i-1) = 0 and w(i) = 1; w(i+1:n) is stored on exit in */
/*     A(i+1:n,i), if op(A) = 'N', and in A(i,i+1:n), otherwise. The */
/*     scalar nu is stored in TAUL(i). */

/*     Each GU(i) is a Givens rotator acting on rows i and n+i, */
/*     where the cosine is stored in CSL(2*i-1) and the sine in */
/*     CSL(2*i). */

/*     Each HV(i) has the form */

/*           HV(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in */
/*     Q(i,i+2:n), and tau in Q(i,i+1). */

/*     Each FV(i) has the form */

/*           FV(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in */
/*     B(i,i+2:n), if op(B) = 'N', and in B(i+2:n,i), otherwise. */
/*     The scalar nu is stored in TAUR(i). */

/*     Each GV(i) is a Givens rotator acting on columns i+1 and n+i+1, */
/*     where the cosine is stored in CSR(2*i-1) and the sine in */
/*     CSR(2*i). */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 80/3 N**3 + 20 N**2 + O(N) floating point */
/*     operations and is numerically backward stable. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DGESUV). */

/*     KEYWORDS */

/*     Elementary matrix operations, Matrix decompositions, Hamiltonian */
/*     matrix */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --csl;
    --csr;
    --taul;
    --taur;
    --dwork;

    /* Function Body */
    *info = 0;
    ltra = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
    ltrb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
    if (! ltra && ! lsame_(trana, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ltrb && ! lsame_(tranb, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldg < max(1,*n)) {
	*info = -10;
    } else if (*ldq < max(1,*n)) {
	*info = -12;
    } else if (*ldwork < max(1,*n)) {
	dwork[1] = (doublereal) max(1,*n);
	*info = -18;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04TS", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

    i__1 = *n;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	alpha = q[i__ + i__ * q_dim1];
	if (i__ < *n) {

/*           Generate elementary reflector HU(i) to annihilate Q(i+1:n,i) */

	    i__2 = *n - i__ + 1;
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &nu);

/*           Apply HU(i) from the left. */

	    q[i__ + i__ * q_dim1] = 1.;
	    i__2 = *n - i__ + 1;
	    i__3 = *n - i__;
	    dlarf_("Left", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &nu, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[1], (ftnlen)4);
	    if (ltra) {
		i__2 = *n - i__ + 1;
		i__3 = *n - i__ + 1;
		dlarf_("Right", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &
			nu, &a[i__ + i__ * a_dim1], lda, &dwork[1], (ftnlen)5)
			;
	    } else {
		i__2 = *n - i__ + 1;
		i__3 = *n - i__ + 1;
		dlarf_("Left", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &
			nu, &a[i__ + i__ * a_dim1], lda, &dwork[1], (ftnlen)4)
			;
	    }
	    if (ltrb) {
		i__2 = *n - i__ + 1;
		dlarf_("Right", n, &i__2, &q[i__ + i__ * q_dim1], &c__1, &nu, 
			&b[i__ * b_dim1 + 1], ldb, &dwork[1], (ftnlen)5);
	    } else {
		i__2 = *n - i__ + 1;
		dlarf_("Left", &i__2, n, &q[i__ + i__ * q_dim1], &c__1, &nu, &
			b[i__ + b_dim1], ldb, &dwork[1], (ftnlen)4);
	    }
	    i__2 = *n - i__ + 1;
	    dlarf_("Left", &i__2, n, &q[i__ + i__ * q_dim1], &c__1, &nu, &g[
		    i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
	    q[i__ + i__ * q_dim1] = nu;
	} else {
	    q[i__ + i__ * q_dim1] = 0.;
	}

/*        Generate symplectic Givens rotator GU(i) to annihilate Q(i,i). */

	temp = a[i__ + i__ * a_dim1];
	dlartg_(&temp, &alpha, &c__, &s, &a[i__ + i__ * a_dim1]);

/*        Apply G(i) from the left. */

	if (ltra) {
	    i__2 = *n - i__;
	    drot_(&i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &q[i__ + (i__ + 1)
		     * q_dim1], ldq, &c__, &s);
	} else {
	    i__2 = *n - i__;
	    drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (i__ + 1)
		     * q_dim1], ldq, &c__, &s);
	}
	if (ltrb) {
	    drot_(n, &g[i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &c__1, &c__,
		     &s);
	} else {
	    drot_(n, &g[i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &c__, &s);
	}
	csl[(i__ << 1) - 1] = c__;
	csl[i__ * 2] = s;

	if (i__ < *n) {
	    if (ltra) {

/*              Generate elementary reflector FU(i) to annihilate */
/*              A(i,i+1:n). */

		i__2 = *n - i__ + 1;
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + (i__ + 1) * 
			a_dim1], lda, &taul[i__]);

/*              Apply FU(i) from the left. */

		temp = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__2 = *n - i__;
		i__3 = *n - i__ + 1;
		dlarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taul[i__], &a[i__ + 1 + i__ * a_dim1], lda, &dwork[1],
			 (ftnlen)5);
		i__2 = *n - i__ + 1;
		i__3 = *n - i__;
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taul[i__], &q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[
			1], (ftnlen)4);
		if (ltrb) {
		    i__2 = *n - i__ + 1;
		    dlarf_("Right", n, &i__2, &a[i__ + i__ * a_dim1], lda, &
			    taul[i__], &b[i__ * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *n - i__ + 1;
		    dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], lda, &
			    taul[i__], &b[i__ + b_dim1], ldb, &dwork[1], (
			    ftnlen)4);
		}
		i__2 = *n - i__ + 1;
		dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], lda, &taul[
			i__], &g[i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
		a[i__ + i__ * a_dim1] = temp;
	    } else {

/*              Generate elementary reflector FU(i) to annihilate */
/*              A(i+1:n,i). */

		i__2 = *n - i__ + 1;
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * 
			a_dim1], &c__1, &taul[i__]);

/*              Apply FU(i) from the left. */

		temp = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		i__2 = *n - i__ + 1;
		i__3 = *n - i__;
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			taul[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &dwork[
			1], (ftnlen)4);
		i__2 = *n - i__ + 1;
		i__3 = *n - i__;
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			taul[i__], &q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[
			1], (ftnlen)4);
		if (ltrb) {
		    i__2 = *n - i__ + 1;
		    dlarf_("Right", n, &i__2, &a[i__ + i__ * a_dim1], &c__1, &
			    taul[i__], &b[i__ * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)5);
		} else {
		    i__2 = *n - i__ + 1;
		    dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], &c__1, &
			    taul[i__], &b[i__ + b_dim1], ldb, &dwork[1], (
			    ftnlen)4);
		}
		i__2 = *n - i__ + 1;
		dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], &c__1, &taul[
			i__], &g[i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
		a[i__ + i__ * a_dim1] = temp;
	    }
	} else {
	    taul[i__] = 0.;
	}
	if (i__ < *n) {
	    alpha = q[i__ + (i__ + 1) * q_dim1];
	}
	if (i__ < *n - 1) {

/*           Generate elementary reflector HV(i) to annihilate Q(i,i+2:n) */

	    i__2 = *n - i__;
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &nu);

/*           Apply HV(i) from the right. */

	    q[i__ + (i__ + 1) * q_dim1] = 1.;
	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dlarf_("Right", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    nu, &q[i__ + 1 + (i__ + 1) * q_dim1], ldq, &dwork[1], (
		    ftnlen)5);
	    if (ltra) {
		i__2 = *n - i__;
		dlarf_("Left", &i__2, n, &q[i__ + (i__ + 1) * q_dim1], ldq, &
			nu, &a[i__ + 1 + a_dim1], lda, &dwork[1], (ftnlen)4);
	    } else {
		i__2 = *n - i__;
		dlarf_("Right", n, &i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &
			nu, &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (
			ftnlen)5);
	    }
	    if (ltrb) {
		i__2 = *n - i__;
		i__3 = *n - i__ + 1;
		dlarf_("Left", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], 
			ldq, &nu, &b[i__ + 1 + i__ * b_dim1], ldb, &dwork[1], 
			(ftnlen)4);
	    } else {
		i__2 = *n - i__ + 1;
		i__3 = *n - i__;
		dlarf_("Right", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], 
			ldq, &nu, &b[i__ + (i__ + 1) * b_dim1], ldb, &dwork[1]
			, (ftnlen)5);
	    }
	    i__2 = *n - i__;
	    dlarf_("Right", n, &i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &nu, 
		    &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1], (ftnlen)5);
	    q[i__ + (i__ + 1) * q_dim1] = nu;
	} else if (i__ < *n) {
	    q[i__ + (i__ + 1) * q_dim1] = 0.;
	}
	if (i__ < *n) {

/*           Generate symplectic Givens rotator GV(i) to annihilate */
/*           Q(i,i+1). */

	    if (ltrb) {
		temp = b[i__ + 1 + i__ * b_dim1];
		dlartg_(&temp, &alpha, &c__, &s, &b[i__ + 1 + i__ * b_dim1]);
		s = -s;
		i__2 = *n - i__;
		drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ 
			+ 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
	    } else {
		temp = b[i__ + (i__ + 1) * b_dim1];
		dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (i__ + 1) * b_dim1])
			;
		s = -s;
		i__2 = *n - i__;
		drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ 
			+ 1 + (i__ + 1) * b_dim1], &c__1, &c__, &s);
	    }
	    if (ltra) {
		drot_(n, &a[i__ + 1 + a_dim1], lda, &g[(i__ + 1) * g_dim1 + 1]
			, &c__1, &c__, &s);
	    } else {
		drot_(n, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(i__ + 1) * 
			g_dim1 + 1], &c__1, &c__, &s);
	    }
	    csr[(i__ << 1) - 1] = c__;
	    csr[i__ * 2] = s;
	}
	if (i__ < *n - 1) {
	    if (ltrb) {

/*              Generate elementary reflector FV(i) to annihilate */
/*              B(i+2:n,i). */

		i__2 = *n - i__;
		dlarfg_(&i__2, &b[i__ + 1 + i__ * b_dim1], &b[i__ + 2 + i__ * 
			b_dim1], &c__1, &taur[i__]);

/*              Apply FV(i) from the right. */

		temp = b[i__ + 1 + i__ * b_dim1];
		b[i__ + 1 + i__ * b_dim1] = 1.;
		i__2 = *n - i__;
		i__3 = *n - i__;
		dlarf_("Left", &i__2, &i__3, &b[i__ + 1 + i__ * b_dim1], &
			c__1, &taur[i__], &b[i__ + 1 + (i__ + 1) * b_dim1], 
			ldb, &dwork[1], (ftnlen)4);
		i__2 = *n - i__;
		i__3 = *n - i__;
		dlarf_("Right", &i__2, &i__3, &b[i__ + 1 + i__ * b_dim1], &
			c__1, &taur[i__], &q[i__ + 1 + (i__ + 1) * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
		if (ltra) {
		    i__2 = *n - i__;
		    dlarf_("Left", &i__2, n, &b[i__ + 1 + i__ * b_dim1], &
			    c__1, &taur[i__], &a[i__ + 1 + a_dim1], lda, &
			    dwork[1], (ftnlen)4);
		} else {
		    i__2 = *n - i__;
		    dlarf_("Right", n, &i__2, &b[i__ + 1 + i__ * b_dim1], &
			    c__1, &taur[i__], &a[(i__ + 1) * a_dim1 + 1], lda,
			     &dwork[1], (ftnlen)5);
		}
		i__2 = *n - i__;
		dlarf_("Right", n, &i__2, &b[i__ + 1 + i__ * b_dim1], &c__1, &
			taur[i__], &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1],
			 (ftnlen)5);
		b[i__ + 1 + i__ * b_dim1] = temp;
	    } else {

/*              Generate elementary reflector FV(i) to annihilate */
/*              B(i,i+2:n). */

		i__2 = *n - i__;
		dlarfg_(&i__2, &b[i__ + (i__ + 1) * b_dim1], &b[i__ + (i__ + 
			2) * b_dim1], ldb, &taur[i__]);

/*              Apply FV(i) from the right. */

		temp = b[i__ + (i__ + 1) * b_dim1];
		b[i__ + (i__ + 1) * b_dim1] = 1.;
		i__2 = *n - i__;
		i__3 = *n - i__;
		dlarf_("Right", &i__2, &i__3, &b[i__ + (i__ + 1) * b_dim1], 
			ldb, &taur[i__], &b[i__ + 1 + (i__ + 1) * b_dim1], 
			ldb, &dwork[1], (ftnlen)5);
		i__2 = *n - i__;
		i__3 = *n - i__;
		dlarf_("Right", &i__2, &i__3, &b[i__ + (i__ + 1) * b_dim1], 
			ldb, &taur[i__], &q[i__ + 1 + (i__ + 1) * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
		if (ltra) {
		    i__2 = *n - i__;
		    dlarf_("Left", &i__2, n, &b[i__ + (i__ + 1) * b_dim1], 
			    ldb, &taur[i__], &a[i__ + 1 + a_dim1], lda, &
			    dwork[1], (ftnlen)4);
		} else {
		    i__2 = *n - i__;
		    dlarf_("Right", n, &i__2, &b[i__ + (i__ + 1) * b_dim1], 
			    ldb, &taur[i__], &a[(i__ + 1) * a_dim1 + 1], lda, 
			    &dwork[1], (ftnlen)5);
		}
		i__2 = *n - i__;
		dlarf_("Right", n, &i__2, &b[i__ + (i__ + 1) * b_dim1], ldb, &
			taur[i__], &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1],
			 (ftnlen)5);
		b[i__ + (i__ + 1) * b_dim1] = temp;
	    }
	} else if (i__ < *n) {
	    taur[i__] = 0.;
	}
/* L10: */
    }
    dwork[1] = (doublereal) max(1,*n);
    return 0;
/* *** Last line of MB04TS *** */
} /* mb04ts_ */

