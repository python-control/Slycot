/* AB04MD.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = -1.;
static doublereal c_b10 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int ab04md_(char *type__, integer *n, integer *m, integer *p,
	 doublereal *alpha, doublereal *beta, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen type_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, ip;
    static doublereal ab2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal pbeta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical ltype;
    static doublereal sqrab2;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal palpha;
    extern /* Subroutine */ int dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen), dgetri_(integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *), dgetrs_(char *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);


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

/*     To perform a transformation on the parameters (A,B,C,D) of a */
/*     system, which is equivalent to a bilinear transformation of the */
/*     corresponding transfer function matrix. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             Indicates the type of the original system and the */
/*             transformation to be performed as follows: */
/*             = 'D':  discrete-time   -> continuous-time; */
/*             = 'C':  continuous-time -> discrete-time. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     ALPHA,  (input) DOUBLE PRECISION */
/*     BETA    Parameters specifying the bilinear transformation. */
/*             Recommended values for stable systems: ALPHA = 1, */
/*             BETA = 1.  ALPHA <> 0, BETA <> 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state matrix A of the original system. */
/*             On exit, the leading N-by-N part of this array contains */
/*                              _ */
/*             the state matrix A of the transformed system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the original system. */
/*             On exit, the leading N-by-M part of this array contains */
/*                              _ */
/*             the input matrix B of the transformed system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C of the original system. */
/*             On exit, the leading P-by-N part of this array contains */
/*                               _ */
/*             the output matrix C of the transformed system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix D for the original system. */
/*             On exit, the leading P-by-M part of this array contains */
/*                                     _ */
/*             the input/output matrix D of the transformed system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */
/*             For optimum performance LDWORK >= MAX(1,N*NB), where NB */
/*             is the optimal blocksize. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix (ALPHA*I + A) is exactly singular; */
/*             = 2:  if the matrix  (BETA*I - A) is exactly singular. */

/*     METHOD */

/*     The parameters of the discrete-time system are transformed into */
/*     the parameters of the continuous-time system (TYPE = 'D'), or */
/*     vice-versa (TYPE = 'C') by the transformation: */

/*     1.  Discrete -> continuous */
/*         _                     -1 */
/*         A = beta*(alpha*I + A)  * (A - alpha*I) */
/*         _                                     -1 */
/*         B = sqrt(2*alpha*beta) * (alpha*I + A)  * B */
/*         _                                         -1 */
/*         C = sqrt(2*alpha*beta) * C * (alpha*I + A) */
/*         _                        -1 */
/*         D = D - C * (alpha*I + A)  * B */

/*     which is equivalent to the bilinear transformation */

/*                       z - alpha */
/*         z -> s = beta ---------  . */
/*                       z + alpha */

/*     of one transfer matrix onto the other. */

/*     2.  Continuous -> discrete */
/*         _                     -1 */
/*         A = alpha*(beta*I - A)  * (beta*I + A) */
/*         _                                    -1 */
/*         B = sqrt(2*alpha*beta) * (beta*I - A)  * B */
/*         _                                        -1 */
/*         C = sqrt(2*alpha*beta) * C * (beta*I - A) */
/*         _                       -1 */
/*         D = D + C * (beta*I - A)  * B */

/*     which is equivalent to the bilinear transformation */

/*                      beta + s */
/*       s -> z = alpha -------- . */
/*                      beta - s */

/*     of one transfer matrix onto the other. */

/*     REFERENCES */

/*     [1] Al-Saggaf, U.M. and Franklin, G.F. */
/*         Model reduction via balanced realizations: a extension and */
/*         frequency weighting techniques. */
/*         IEEE Trans. Autom. Contr., AC-33, pp. 687-692, 1988. */

/*     NUMERICAL ASPECTS */
/*                                                      3 */
/*     The time taken is approximately proportional to N . */
/*     The accuracy depends mainly on the condition number of the matrix */
/*     to be inverted. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and */
/*                  A. Varga, German Aerospace Research Establishment, */
/*                  Oberpfaffenhofen, Germany, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB04AD by W. van der Linden, and */
/*     A.J. Geurts, Technische Hogeschool Eindhoven, Holland. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bilinear transformation, continuous-time system, discrete-time */
/*     system, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

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
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    ltype = lsame_(type__, "D", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! ltype && ! lsame_(type__, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*alpha == 0.) {
	*info = -5;
    } else if (*beta == 0.) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldc < max(1,*p)) {
	*info = -12;
    } else if (*ldd < max(1,*p)) {
	*info = -14;
    } else if (*ldwork < max(1,*n)) {
	*info = -17;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB04MD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MAX */
    i__1 = max(*n,*m);
    if (max(i__1,*p) == 0) {
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    if (ltype) {

/*        Discrete-time to continuous-time with (ALPHA, BETA). */

	palpha = *alpha;
	pbeta = *beta;
    } else {

/*        Continuous-time to discrete-time with (ALPHA, BETA) is */
/*        equivalent with discrete-time to continuous-time with */
/*        (-BETA, -ALPHA), if B and C change the sign. */

	palpha = -(*beta);
	pbeta = -(*alpha);
    }

    ab2 = palpha * pbeta * 2.;
    d__1 = sqrt((abs(ab2)));
    sqrab2 = d_sign(&d__1, &palpha);
/*                          -1 */
/*     Compute (alpha*I + A)  . */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__ + i__ * a_dim1] += palpha;
/* L10: */
    }

    dgetrf_(n, n, &a[a_offset], lda, &iwork[1], info);

    if (*info != 0) {

/*        Error return. */

	if (ltype) {
	    *info = 1;
	} else {
	    *info = 2;
	}
	return 0;
    }
/*                         -1 */
/*     Compute  (alpha*I+A)  *B. */

    dgetrs_("No transpose", n, m, &a[a_offset], lda, &iwork[1], &b[b_offset], 
	    ldb, info, (ftnlen)12);
/*                               -1 */
/*     Compute  D - C*(alpha*I+A)  *B. */

    dgemm_("No transpose", "No transpose", p, m, n, &c_b9, &c__[c_offset], 
	    ldc, &b[b_offset], ldb, &c_b10, &d__[d_offset], ldd, (ftnlen)12, (
	    ftnlen)12);

/*     Scale B by  sqrt(2*alpha*beta). */

    dlascl_("General", &c__0, &c__0, &c_b10, &sqrab2, n, m, &b[b_offset], ldb,
	     info, (ftnlen)7);
/*                                                -1 */
/*     Compute  sqrt(2*alpha*beta)*C*(alpha*I + A)  . */

    dtrsm_("Right", "Upper", "No transpose", "Non-unit", p, n, &sqrab2, &a[
	    a_offset], lda, &c__[c_offset], ldc, (ftnlen)5, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

    dtrsm_("Right", "Lower", "No transpose", "Unit", p, n, &c_b10, &a[
	    a_offset], lda, &c__[c_offset], ldc, (ftnlen)5, (ftnlen)5, (
	    ftnlen)12, (ftnlen)4);

/*     Apply column interchanges to the solution matrix. */

    for (i__ = *n - 1; i__ >= 1; --i__) {
	ip = iwork[i__];
	if (ip != i__) {
	    dswap_(p, &c__[i__ * c_dim1 + 1], &c__1, &c__[ip * c_dim1 + 1], &
		    c__1);
	}
/* L20: */
    }
/*                               -1 */
/*     Compute beta*(alpha*I + A)  *(A - alpha*I) as */
/*                                        -1 */
/*     beta*I - 2*alpha*beta*(alpha*I + A)  . */

/*     Workspace: need N;  prefer N*NB. */

    dgetri_(n, &a[a_offset], lda, &iwork[1], &dwork[1], ldwork, info);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = -ab2;
	dscal_(n, &d__1, &a[i__ * a_dim1 + 1], &c__1);
	a[i__ + i__ * a_dim1] += pbeta;
/* L30: */
    }

    return 0;
/* *** Last line of AB04MD *** */
} /* ab04md_ */

