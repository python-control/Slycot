/* AB09AD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ab09ad_(char *dico, char *job, char *equil, char *ordsel,
	 integer *n, integer *m, integer *p, integer *nr, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *hsv, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	dico_len, ftnlen job_len, ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer ki, nn, kr, kt, kw, kti, ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab09ax_(char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01wd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal maxred;
    static logical fixord;
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

/*     To compute a reduced order model (Ar,Br,Cr) for a stable original */
/*     state-space representation (A,B,C) by using either the square-root */
/*     or the balancing-free square-root Balance & Truncate (B & T) */
/*     model reduction method. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'N':  use the balancing-free square-root */
/*                     Balance & Truncate method. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR */
/*             is the desired order on entry and NMIN is the order of a */
/*             minimal realization of the given system; NMIN is */
/*             determined as the number of Hankel singular values greater */
/*             than N*EPS*HNORM(A,B,C), where EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH) and */
/*             HNORM(A,B,C) is the Hankel norm of the system (computed */
/*             in HSV(1)); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL,N*EPS*HNORM(A,B,C)). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the reduced */
/*             order system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, if INFO = 0, the leading NR-by-M part of this */
/*             array contains the input/state matrix Br of the reduced */
/*             order system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, if INFO = 0, the leading P-by-NR part of this */
/*             array contains the state/output matrix Cr of the reduced */
/*             order system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the Hankel singular values of */
/*             the original system ordered decreasingly. HSV(1) is the */
/*             Hankel norm of the system. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL = c*HNORM(A,B,C), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(A,B,C) is the */
/*             Hankel-norm of the given system (computed in HSV(1)). */
/*             For computing a minimal realization, the recommended */
/*             value is TOL = N*EPS*HNORM(A,B,C), where EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH). */
/*             This value is used by default if TOL <= 0 on entry. */
/*             If ORDSEL = 'F', the value of TOL is ignored. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = 0, if JOB = 'B'; */
/*             LIWORK = N, if JOB = 'N'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than the order of a minimal realization of the */
/*                   given system. In this case, the resulting NR is */
/*                   set automatically to a value corresponding to the */
/*                   order of a minimal realization of the system. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to the real Schur form failed; */
/*             = 2:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 3:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t)                               (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09AD determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t)                             (2) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     If JOB = 'B', the square-root Balance & Truncate method of [1] */
/*     is used and, for DICO = 'C', the resulting model is balanced. */
/*     By setting TOL <= 0, the routine can be used to compute balanced */
/*     minimal state-space realizations of stable systems. */

/*     If JOB = 'N', the balancing-free square-root version of the */
/*     Balance & Truncate method [2] is used. */
/*     By setting TOL <= 0, the routine can be used to compute minimal */
/*     state-space realizations of stable systems. */

/*     REFERENCES */

/*     [1] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [2] Varga A. */
/*         Efficient minimal realization procedure based on balancing. */
/*         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991, */
/*         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), */
/*         Vol. 2, pp. 42-46. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routines SRBT and SRBFT. */

/*     REVISIONS */

/*     May 2, 1998. */
/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */

/*     KEYWORDS */

/*     Balancing, minimal state-space representation, model reduction, */
/*     multivariable system, state-space model. */

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
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    *iwarn = 0;
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || lsame_(dico, "D", (
	    ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "N", (
	    ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*p < 0) {
	*info = -7;
    } else if (fixord && (*nr < 0 || *nr > *n)) {
	*info = -8;
    } else if (*lda < max(1,*n)) {
	*info = -10;
    } else if (*ldb < max(1,*n)) {
	*info = -12;
    } else if (*ldc < max(1,*p)) {
	*info = -14;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = max(*n,*m);
	i__1 = 1, i__2 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) /
		 2;
	if (*ldwork < max(i__1,i__2)) {
	    *info = -19;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB09AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*m);
    if (min(i__1,*p) == 0 || fixord && *nr == 0) {
	*nr = 0;
	dwork[1] = 1.;
	return 0;
    }

/*     Allocate working storage. */

    nn = *n * *n;
    kt = 1;
    kr = kt + nn;
    ki = kr + *n;
    kw = ki + *n;

    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */

	maxred = 100.;
	tb01id_("A", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb, &
		c__[c_offset], ldc, &dwork[1], info, (ftnlen)1);
    }

/*     Reduce A to the real Schur form using an orthogonal similarity */
/*     transformation A <- T'*A*T and apply the transformation to */
/*     B and C: B <- T'*B and C <- C*T. */

    i__1 = *ldwork - kw + 1;
    tb01wd_(n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], 
	    ldc, &dwork[kt], n, &dwork[kr], &dwork[ki], &dwork[kw], &i__1, &
	    ierr);
    if (ierr != 0) {
	*info = 1;
	return 0;
    }

    wrkopt = dwork[kw] + (doublereal) (kw - 1);
    kti = kt + nn;
    kw = kti + nn;

    i__1 = *ldwork - kw + 1;
    ab09ax_(dico, job, ordsel, n, m, p, nr, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, &hsv[1], &dwork[kt], n, &dwork[kti], n, 
	    tol, &iwork[1], &dwork[kw], &i__1, iwarn, &ierr, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);

    if (ierr != 0) {
	*info = ierr + 1;
	return 0;
    }

/* Computing MAX */
    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
    dwork[1] = max(d__1,d__2);

    return 0;
/* *** Last line of AB09AD *** */
} /* ab09ad_ */

