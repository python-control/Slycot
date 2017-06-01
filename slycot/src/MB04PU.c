/* MB04PU.f -- translated by f2c (version 20100827).
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
static doublereal c_b7 = 0.;
static doublereal c_b14 = -1.;

/* Subroutine */ int mb04pu_(integer *n, integer *ilo, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *cs, doublereal *tau, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal s, mu, nu;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr2_(char 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), daxpy_(integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen);


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

/*     To reduce a Hamiltonian matrix, */

/*                   [  A   G  ] */
/*              H =  [       T ] , */
/*                   [  Q  -A  ] */

/*     where A is an N-by-N matrix and G,Q are N-by-N symmetric matrices, */
/*     to Paige/Van Loan (PVL) form. That is, an orthogonal symplectic U */
/*     is computed so that */

/*               T       [  Aout   Gout  ] */
/*              U H U =  [             T ] , */
/*                       [  Qout  -Aout  ] */

/*     where Aout is upper Hessenberg and Qout is diagonal. */
/*     Unblocked version. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     ILO     (input) INTEGER */
/*             It is assumed that A is already upper triangular and Q is */
/*             zero in rows and columns 1:ILO-1. ILO is normally set by a */
/*             previous call to MB04DD; otherwise it should be set to 1. */
/*             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Aout and, in the zero part of Aout, */
/*             information about the elementary reflectors used to */
/*             compute the PVL factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain the lower triangular part of the matrix Q and */
/*             the upper triangular part of the matrix G. */
/*             On exit, the leading N-by-N+1 part of this array contains */
/*             the diagonal of the matrix Qout, the upper triangular part */
/*             of the matrix Gout and, in the zero parts of Qout, */
/*             information about the elementary reflectors used to */
/*             compute the PVL factorization. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     CS      (output) DOUBLE PRECISION array, dimension (2N-2) */
/*             On exit, the first 2N-2 elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations used */
/*             to compute the PVL factorization. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
/*             On exit, the first N-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -10,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N-1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix U is represented as a product of symplectic reflectors */
/*     and Givens rotators */

/*     U = diag( H(1),H(1) )     G(1)   diag( F(1),F(1) ) */
/*         diag( H(2),H(2) )     G(2)   diag( F(2),F(2) ) */
/*                                .... */
/*         diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ). */

/*     Each H(i) has the form */

/*           H(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in */
/*     QG(i+2:n,i), and tau in QG(i+1,i). */

/*     Each F(i) has the form */

/*           F(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in */
/*     A(i+2:n,i), and nu in TAU(i). */

/*     Each G(i) is a Givens rotator acting on rows i+1 and n+i+1, */
/*     where the cosine is stored in CS(2*i-1) and the sine in */
/*     CS(2*i). */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 40/3 N**3 + O(N) floating point operations */
/*     and is strongly backward stable. */

/*     REFERENCES */

/*     [1] C. F. VAN LOAN: */
/*         A symplectic method for approximating all the eigenvalues of */
/*         a Hamiltonian matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     CONTRIBUTORS */

/*     D. Kressner (Technical Univ. Berlin, Germany) and */
/*     P. Benner (Technical Univ. Chemnitz, Germany), December 2003. */

/*     REVISIONS */

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DHAPVL). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix. */

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
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    --cs;
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldqg < max(1,*n)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
	    i__1 = 1, i__2 = *n - 1;
	    dwork[1] = (doublereal) max(i__1,i__2);
	    *info = -10;
	}
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04PU", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n <= *ilo) {
	dwork[1] = 1.;
	return 0;
    }

    i__1 = *n - 1;
    for (i__ = *ilo; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate QG(i+2:n,i). */

	alpha = qg[i__ + 1 + i__ * qg_dim1];
	i__2 = *n - i__;
/* Computing MIN */
	i__3 = i__ + 2;
	dlarfg_(&i__2, &alpha, &qg[min(i__3,*n) + i__ * qg_dim1], &c__1, &nu);
	if (nu != 0.) {
	    qg[i__ + 1 + i__ * qg_dim1] = 1.;

/*           Apply H(i) from both sides to QG(i+1:n,i+1:n). */
/*           Compute  x := nu * QG(i+1:n,i+1:n) * v. */

	    i__2 = *n - i__;
	    dsymv_("Lower", &i__2, &nu, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &c_b7, &dwork[
		    1], &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * nu * (x'*v) * v. */

	    i__2 = *n - i__;
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &qg[i__ + 1 + i__ *
		     qg_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &mu, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &dwork[1],
		     &c__1);

/*           Apply the transformation as a rank-2 update: */
/*                QG := QG - v * w' - w * v'. */

	    i__2 = *n - i__;
	    dsyr2_("Lower", &i__2, &c_b14, &qg[i__ + 1 + i__ * qg_dim1], &
		    c__1, &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 1) * qg_dim1]
		    , ldqg, (ftnlen)5);

/*           Apply H(i) from the right hand side to QG(1:i,i+2:n+1). */

	    i__2 = *n - i__;
	    dlarf_("Right", &i__, &i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1, 
		    &nu, &qg[(i__ + 2) * qg_dim1 + 1], ldqg, &dwork[1], (
		    ftnlen)5);

/*           Apply H(i) from both sides to QG(i+1:n,i+2:n+1). */
/*           Compute  x := nu * QG(i+1:n,i+2:n+1) * v. */

	    i__2 = *n - i__;
	    dsymv_("Upper", &i__2, &nu, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &c_b7, &dwork[
		    1], &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * nu * (x'*v) * v. */

	    i__2 = *n - i__;
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &qg[i__ + 1 + i__ *
		     qg_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &mu, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &dwork[1],
		     &c__1);

/*           Apply the transformation as a rank-2 update: */
/*              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'. */

	    i__2 = *n - i__;
	    dsyr2_("Upper", &i__2, &c_b14, &qg[i__ + 1 + i__ * qg_dim1], &
		    c__1, &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 2) * qg_dim1]
		    , ldqg, (ftnlen)5);

/*           Apply H(i) from the left hand side to A(i+1:n,i:n). */

	    i__2 = *n - i__;
	    i__3 = *n - i__ + 1;
	    dlarf_("Left", &i__2, &i__3, &qg[i__ + 1 + i__ * qg_dim1], &c__1, 
		    &nu, &a[i__ + 1 + i__ * a_dim1], lda, &dwork[1], (ftnlen)
		    4);

/*           Apply H(i) from the right hand side to A(1:n,i+1:n). */

	    i__2 = *n - i__;
	    dlarf_("Right", n, &i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
		    nu, &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (ftnlen)5)
		    ;
	}
	qg[i__ + 1 + i__ * qg_dim1] = nu;

/*        Generate symplectic Givens rotation G(i) to annihilate */
/*        QG(i+1,i). */

	temp = a[i__ + 1 + i__ * a_dim1];
	dlartg_(&temp, &alpha, &c__, &s, &a[i__ + 1 + i__ * a_dim1]);

/*        Apply G(i) to [A(I+1,I+2:N); QG(I+2:N,I+1)']. */

	i__2 = *n - i__ - 1;
	drot_(&i__2, &a[i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[i__ + 2 + (
		i__ + 1) * qg_dim1], &c__1, &c__, &s);

/*        Apply G(i) to [A(1:I,I+1) QG(1:I,I+2)]. */

	drot_(&i__, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		qg_dim1 + 1], &c__1, &c__, &s);

/*        Apply G(i) to [A(I+2:N,I+1) QG(I+1, I+3:N+1)'] from the right. */

	i__2 = *n - i__ - 1;
	drot_(&i__2, &a[i__ + 2 + (i__ + 1) * a_dim1], &c__1, &qg[i__ + 1 + (
		i__ + 3) * qg_dim1], ldqg, &c__, &s);

/*        Fix the diagonal part. */

	temp = a[i__ + 1 + (i__ + 1) * a_dim1];
	ttemp = qg[i__ + 1 + (i__ + 2) * qg_dim1];
	a[i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[i__ + 1 + (i__ 
		+ 1) * qg_dim1];
	qg[i__ + 1 + (i__ + 2) * qg_dim1] = c__ * ttemp - s * temp;
	qg[i__ + 1 + (i__ + 1) * qg_dim1] = -s * temp + c__ * qg[i__ + 1 + (
		i__ + 1) * qg_dim1];
	ttemp = -s * ttemp - c__ * temp;
	temp = a[i__ + 1 + (i__ + 1) * a_dim1];
	qg[i__ + 1 + (i__ + 1) * qg_dim1] = c__ * qg[i__ + 1 + (i__ + 1) * 
		qg_dim1] + s * ttemp;
	a[i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[i__ + 1 + (i__ 
		+ 2) * qg_dim1];
	qg[i__ + 1 + (i__ + 2) * qg_dim1] = -s * temp + c__ * qg[i__ + 1 + (
		i__ + 2) * qg_dim1];
	cs[(i__ << 1) - 1] = c__;
	cs[i__ * 2] = s;

/*        Generate elementary reflector F(i) to annihilate A(i+2:n,i). */

	i__2 = *n - i__;
/* Computing MIN */
	i__3 = i__ + 2;
	dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ * 
		a_dim1], &c__1, &nu);
	if (nu != 0.) {
	    temp = a[i__ + 1 + i__ * a_dim1];
	    a[i__ + 1 + i__ * a_dim1] = 1.;

/*           Apply F(i) from the left hand side to A(i+1:n,i+1:n). */

	    i__2 = *n - i__;
	    i__3 = *n - i__;
	    dlarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &
		    nu, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &dwork[1], (
		    ftnlen)4);

/*           Apply G(i) from the right hand side to A(1:n,i+1:n). */

	    i__2 = *n - i__;
	    dlarf_("Right", n, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &nu, 
		    &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (ftnlen)5);

/*           Apply G(i) from both sides to QG(i+1:n,i+1:n). */
/*           Compute  x := nu * QG(i+1:n,i+1:n) * v. */

	    i__2 = *n - i__;
	    dsymv_("Lower", &i__2, &nu, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b7, &dwork[1],
		     &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * tau * (x'*v) * v. */

	    i__2 = *n - i__;
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &a[i__ + 1 + i__ * 
		    a_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &mu, &a[i__ + 1 + i__ * a_dim1], &c__1, &dwork[1], &
		    c__1);

/*           Apply the transformation as a rank-2 update: */
/*                QG := QG - v * w' - w * v'. */

	    i__2 = *n - i__;
	    dsyr2_("Lower", &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, 
		    &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, (ftnlen)5);

/*           Apply G(i) from the right hand side to QG(1:i,i+2:n+1). */

	    i__2 = *n - i__;
	    dlarf_("Right", &i__, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &
		    nu, &qg[(i__ + 2) * qg_dim1 + 1], ldqg, &dwork[1], (
		    ftnlen)5);

/*           Apply G(i) from both sides to QG(i+1:n,i+2:n+1). */
/*           Compute  x := nu * QG(i+1:n,i+2:n+1) * v. */

	    i__2 = *n - i__;
	    dsymv_("Upper", &i__2, &nu, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b7, &dwork[1],
		     &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * tau * (x'*v) * v. */

	    i__2 = *n - i__;
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &a[i__ + 1 + i__ * 
		    a_dim1], &c__1);
	    i__2 = *n - i__;
	    daxpy_(&i__2, &mu, &a[i__ + 1 + i__ * a_dim1], &c__1, &dwork[1], &
		    c__1);

/*           Apply the transformation as a rank-2 update: */
/*              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'. */

	    i__2 = *n - i__;
	    dsyr2_("Upper", &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, 
		    &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, (ftnlen)5);
	    a[i__ + 1 + i__ * a_dim1] = temp;
	}
	tau[i__] = nu;
/* L10: */
    }
/* Computing MAX */
    i__1 = 1, i__2 = *n - 1;
    dwork[1] = (doublereal) max(i__1,i__2);
    return 0;
/* *** Last line of MB04PU *** */
} /* mb04pu_ */

