/* AB13AX.f -- translated by f2c (version 20100827).
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

static logical c_false = FALSE_;
static logical c_true = TRUE_;
static integer c__1 = 1;

doublereal ab13ax_(char *dico, integer *n, integer *m, integer *p, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *hsv, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j, kr, ks, ku, kw, ierr, ktau, mnmp;
    extern /* Subroutine */ int ma02dd_(char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), mb03ud_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dtpmv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

/*     To compute the Hankel-norm of the transfer-function matrix G of */
/*     a stable state-space system (A,B,C). The state dynamics matrix A */
/*     of the given system is an upper quasi-triangular matrix in */
/*     real Schur form. */

/*     FUNCTION VALUE */

/*     AB13AX  DOUBLE PRECISION */
/*             The Hankel-norm of G (if INFO = 0). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, i.e. the */
/*             order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state dynamics matrix A in a real Schur canonical form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, this array contains the Hankel singular */
/*             values of the given system ordered decreasingly. */
/*             HSV(1) is the Hankel norm of the given system. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(MAX(N,M,P)+5)+N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 2:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t)                               (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system, and let G be the corresponding */
/*     transfer-function matrix. The Hankel-norm of G is computed as the */
/*     the maximum Hankel singular value of the system (A,B,C). */
/*     The computation of the Hankel singular values is performed */
/*     by using the square-root method of [1]. */

/*     REFERENCES */

/*     [1] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The implemented method relies on a square-root technique. */
/*                                     3 */
/*     The algorithms require about 17N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine SHANRM. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Feb. 2000, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Multivariable system, state-space model, system norms. */

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
    --hsv;
    --dwork;

    /* Function Body */
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*p < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    } else if (*ldc < max(1,*p)) {
	*info = -10;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = max(*n,*m);
	i__1 = 1, i__2 = *n * (max(i__3,*p) + 5) + *n * (*n + 1) / 2;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -13;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB13AX", &i__1, (ftnlen)6);
	return ret_val;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0) {
	ret_val = 0.;
	dwork[1] = 1.;
	return ret_val;
    }

/*     Allocate N*MAX(N,M,P), N, and N*(N+1)/2 working storage for the */
/*     matrices S, TAU, and R, respectively. S shares the storage with U. */

    ku = 1;
    ks = 1;
/* Computing MAX */
    i__1 = max(*n,*m);
    mnmp = max(i__1,*p);
    ktau = ks + *n * mnmp;
    kr = ktau + *n;
    kw = kr;

/*     Copy C in U. */

    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], &mnmp, (ftnlen)4);

/*     If DISCR = .FALSE., solve for R the Lyapunov equation */
/*                                  2 */
/*     A'*(R'*R) + (R'*R)*A + scaleo  * C'*C = 0 . */

/*     If DISCR = .TRUE., solve for R the Lyapunov equation */
/*                         2 */
/*     A'*(R'*R)*A + scaleo  * C'*C = R'*R . */

/*     Workspace needed:      N*(MAX(N,M,P)+1); */
/*     Additional workspace:  need   4*N; */
/*                            prefer larger. */

    i__1 = *ldwork - kw + 1;
    sb03ou_(&discr, &c_false, n, p, &a[a_offset], lda, &dwork[ku], &mnmp, &
	    dwork[ktau], &dwork[ku], n, &scaleo, &dwork[kw], &i__1, &ierr);
    if (ierr != 0) {
	*info = 1;
	return ret_val;
    }

    wrkopt = dwork[kw] + (doublereal) (kw - 1);

/*     Pack the upper triangle of R in DWORK(KR). */
/*     Workspace needed:      N*(MAX(N,M,P) + 1) + N*(N+1)/2. */

    ma02dd_("Pack", "Upper", n, &dwork[ku], n, &dwork[kr], (ftnlen)4, (ftnlen)
	    5);

    kw = kr + *n * (*n + 1) / 2;

/*     Copy B in S (over U). */

    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ks], n, (ftnlen)4);

/*     If DISCR = .FALSE., solve for S the Lyapunov equation */
/*                                  2 */
/*     A*(S*S') + (S*S')*A' + scalec *B*B' = 0 . */

/*     If DISCR = .TRUE., solve for S the Lyapunov equation */
/*                         2 */
/*     A*(S*S')*A' + scalec *B*B' = S*S' . */

/*     Workspace needed:      N*(MAX(N,M,P) + 1) + N*(N+1)/2; */
/*     Additional workspace:  need   4*N; */
/*                            prefer larger. */

    i__1 = *ldwork - kw + 1;
    sb03ou_(&discr, &c_true, n, m, &a[a_offset], lda, &dwork[ks], n, &dwork[
	    ktau], &dwork[ks], n, &scalec, &dwork[kw], &i__1, &ierr);

/* Computing MAX */
    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
    wrkopt = max(d__1,d__2);

/*                             | x x | */
/*     Compute R*S in the form | 0 x | in S. Note that R is packed. */

    j = ks;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtpmv_("Upper", "NoTranspose", "NonUnit", &i__, &dwork[kr], &dwork[j],
		 &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
	j += *n;
/* L10: */
    }

/*     Compute the singular values of the upper triangular matrix R*S. */

/*     Workspace needed:      N*MAX(N,M,P); */
/*     Additional workspace:  need   MAX(1,5*N); */
/*                            prefer larger. */

    kw = ktau;
    i__1 = *ldwork - kw + 1;
    mb03ud_("NoVectors", "NoVectors", n, &dwork[ks], n, &dwork[1], &c__1, &
	    hsv[1], &dwork[kw], &i__1, &ierr, (ftnlen)9, (ftnlen)9);
    if (ierr != 0) {
	*info = 2;
	return ret_val;
    }

/*     Scale singular values. */

    d__1 = 1. / scalec / scaleo;
    dscal_(n, &d__1, &hsv[1], &c__1);
    ret_val = hsv[1];

/* Computing MAX */
    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
    dwork[1] = max(d__1,d__2);

    return ret_val;
/* *** Last line of AB13AX *** */
} /* ab13ax_ */

