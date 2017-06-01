/* AB13DD.f -- translated by f2c (version 20100827).
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

static doublereal c_b20 = 1.;
static doublereal c_b21 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b163 = -1.;

/* Subroutine */ int ab13dd_(char *dico, char *jobe, char *equil, char *jobd, 
	integer *n, integer *m, integer *p, doublereal *fpeak, doublereal *a, 
	integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *gpeak, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, doublecomplex *cwork, integer *lcwork, integer *
	info, ftnlen dico_len, ftnlen jobe_len, ftnlen equil_len, ftnlen 
	jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), log(doublereal), 
	    sin(doublereal), cos(doublereal), atan(doublereal);

    /* Local variables */
    static integer i__, j, k, n2, ia, ib, ic, id, ie, ih, ii, im;
    static doublereal pi;
    static integer ir, is, nn, iu, iv, pm;
    static doublereal tm;
    static integer lw, ny, ih12, ihi, ipa, iar, ias, ibs, ibt, ipe, ibv, icu, 
	    ies, ilo, isb, isc, isl, nei;
    static doublereal eps, rat;
    static integer nws, n2pm, imin;
    static doublereal anrm;
    static char vect[1];
    static integer ierr, itau, iter;
    static doublereal enrm, temp[1];
    static integer iwrk;
    static doublereal wmax, rtol;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    tg01ad_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), tg01bd_(char *, 
	    char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal gamma;
    extern doublereal ab13dx_(char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal omega;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb01sd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, ftnlen), dggev_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), mb03xd_(char *, char *, char *, char *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond;
    static logical fulle;
    static doublereal bound, bnorm, cnorm;
    static logical withd;
    static integer minpm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical nodyn;
    static doublereal toler, wrmin;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebal_(
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), dggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal gammal, fpeaki;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal gammas;
    static logical ilascl;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static doublereal fpeaks;
    static logical ilescl;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmax, maxred, bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dhseqr_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    static integer maxcwk;
    static logical lequil;
    extern /* Subroutine */ int dormhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dorgqr_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *);
    static logical usepen;
    static integer mincwr;
    static doublereal anrmto, enrmto;
    static integer minwrk, maxwrk;
    static doublereal smlnum;


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

/*     To compute the L-infinity norm of a continuous-time or */
/*     discrete-time system, either standard or in the descriptor form, */

/*                                     -1 */
/*        G(lambda) = C*( lambda*E - A ) *B + D . */

/*     The norm is finite if and only if the matrix pair (A,E) has no */
/*     eigenvalue on the boundary of the stability domain, i.e., the */
/*     imaginary axis, or the unit circle, respectively. It is assumed */
/*     that the matrix E is nonsingular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system, as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether E is a general square or an identity */
/*             matrix, as follows: */
/*             = 'G':  E is a general square matrix; */
/*             = 'I':  E is the identity matrix. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the system (A,E,B,C) or (A,B,C), as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The row size of the matrix C.  P >= 0. */

/*     FPEAK   (input/output) DOUBLE PRECISION array, dimension (2) */
/*             On entry, this parameter must contain an estimate of the */
/*             frequency where the gain of the frequency response would */
/*             achieve its peak value. Setting FPEAK(2) = 0 indicates an */
/*             infinite frequency. An accurate estimate could reduce the */
/*             number of iterations of the iterative algorithm. If no */
/*             estimate is available, set FPEAK(1) = 0, and FPEAK(2) = 1. */
/*             FPEAK(1) >= 0, FPEAK(2) >= 0. */
/*             On exit, if INFO = 0, this array contains the frequency */
/*             OMEGA, where the gain of the frequency response achieves */
/*             its peak value GPEAK, i.e., */

/*                 || G ( j*OMEGA ) || = GPEAK ,  if DICO = 'C', or */

/*                         j*OMEGA */
/*                 || G ( e       ) || = GPEAK ,  if DICO = 'D', */

/*             where OMEGA = FPEAK(1), if FPEAK(2) > 0, and OMEGA is */
/*             infinite, if FPEAK(2) = 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state dynamics matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             If JOBE = 'G', the leading N-by-N part of this array must */
/*             contain the descriptor matrix E of the system. */
/*             If JOBE = 'I', then E is assumed to be the identity */
/*             matrix and is not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E. */
/*             LDE >= MAX(1,N), if JOBE = 'G'; */
/*             LDE >= 1,        if JOBE = 'I'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If JOBD = 'D', the leading P-by-M part of this array must */
/*             contain the direct transmission matrix D. */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     GPEAK   (output) DOUBLE PRECISION array, dimension (2) */
/*             The L-infinity norm of the system, i.e., the peak gain */
/*             of the frequency response (as measured by the largest */
/*             singular value in the MIMO case), coded in the same way */
/*             as FPEAK. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used to set the accuracy in determining the */
/*             norm.  0 <= TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= K, where K can be computed using the following */
/*             pseudo-code (or the Fortran code included in the routine) */

/*                d = 6*MIN(P,M); */
/*                c = MAX( 4*MIN(P,M) + MAX(P,M), d ); */
/*                if ( MIN(P,M) = 0 ) then */
/*                   K = 1; */
/*                else if( N = 0 or B = 0 or C = 0 ) then */
/*                   if( JOBD = 'D' ) then */
/*                      K = P*M + c; */
/*                   else */
/*                      K = 1; */
/*                   end */
/*                else */
/*                   if ( DICO = 'D' ) then */
/*                      b = 0;  e = d; */
/*                   else */
/*                      b = N*(N+M);  e = c; */
/*                      if ( JOBD = Z' ) then  b = b + P*M;  end */
/*                   end */
/*                   if ( JOBD = 'D' ) then */
/*                      r = P*M; */
/*                      if ( JOBE = 'I', DICO = 'C', */
/*                           N > 0, B <> 0, C <> 0 ) then */
/*                         K = P*P + M*M; */
/*                         r = r + N*(P+M); */
/*                      else */
/*                         K = 0; */
/*                      end */
/*                      K = K + r + c;  r = r + MIN(P,M); */
/*                   else */
/*                      r = 0;  K = 0; */
/*                   end */
/*                   r = r + N*(N+P+M); */
/*                   if ( JOBE = 'G' ) then */
/*                      r = r + N*N; */
/*                      if ( EQUIL = 'S' ) then */
/*                         K = MAX( K, r + 9*N ); */
/*                      end */
/*                      K = MAX( K, r + 4*N + MAX( M, 2*N*N, N+b+e ) ); */
/*                   else */
/*                      K = MAX( K, r + N + */
/*                                  MAX( M, P, N*N+2*N, 3*N+b+e ) ); */
/*                   end */
/*                   w = 0; */
/*                   if ( JOBE = 'I', DICO = 'C' ) then */
/*                      w = r + 4*N*N + 11*N; */
/*                      if ( JOBD = 'D' ) then */
/*                         w = w + MAX(M,P) + N*(P+M); */
/*                      end */
/*                   end */
/*                   if ( JOBE = 'E' or DICO = 'D' or JOBD = 'D' ) then */
/*                      w = MAX( w, r + 6*N + (2*N+P+M)*(2*N+P+M) + */
/*                               MAX( 2*(N+P+M), 8*N*N + 16*N ) ); */
/*                   end */
/*                   K = MAX( 1, K, w, r + 2*N + e ); */
/*                end */

/*             For good performance, LDWORK must generally be larger. */

/*             An easily computable upper bound is */

/*             K = MAX( 1, 15*N*N + P*P + M*M + (6*N+3)*(P+M) + 4*P*M + */
/*                         N*M + 22*N + 7*MIN(P,M) ). */

/*             The smallest workspace is obtained for DICO = 'C', */
/*             JOBE = 'I', and JOBD = 'Z', namely */

/*             K = MAX( 1, N*N + N*P + N*M + N + */
/*                         MAX( N*N + N*M + P*M + 3*N + c, */
/*                              4*N*N + 10*N ) ). */

/*             for which an upper bound is */

/*             K = MAX( 1, 6*N*N + N*P + 2*N*M + P*M + 11*N + MAX(P,M) + */
/*                         6*MIN(P,M) ). */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) contains the optimal */
/*             LCWORK. */

/*     LCWORK  INTEGER */
/*             The dimension of the array CWORK. */
/*             LCWORK >= 1,  if N = 0, or B = 0, or C = 0; */
/*             LCWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)), */
/*                           otherwise. */
/*             For good performance, LCWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the matrix E is (numerically) singular; */
/*             = 2:  the (periodic) QR (or QZ) algorithm for computing */
/*                   eigenvalues did not converge; */
/*             = 3:  the SVD algorithm for computing singular values did */
/*                   not converge; */
/*             = 4:  the tolerance is too small and the algorithm did */
/*                   not converge. */

/*     METHOD */

/*     The routine implements the method presented in [1], with */
/*     extensions and refinements for improving numerical robustness and */
/*     efficiency. Structure-exploiting eigenvalue computations for */
/*     Hamiltonian matrices are used if JOBE = 'I', DICO = 'C', and the */
/*     symmetric matrices to be implicitly inverted are not too ill- */
/*     conditioned. Otherwise, generalized eigenvalue computations are */
/*     used in the iterative algorithm of [1]. */

/*     REFERENCES */

/*     [1] Bruinsma, N.A. and Steinbuch, M. */
/*         A fast algorithm to compute the Hinfinity-norm of a transfer */
/*         function matrix. */
/*         Systems & Control Letters, vol. 14, pp. 287-293, 1990. */

/*     NUMERICAL ASPECTS */

/*     If the algorithm does not converge in MAXIT = 30 iterations */
/*     (INFO = 4), the tolerance must be increased. */

/*     FURTHER COMMENTS */

/*     If the matrix E is singular, other SLICOT Library routines */
/*     could be used before calling AB13DD, for removing the singular */
/*     part of the system. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */
/*     Partly based on SLICOT Library routine AB13CD by P.Hr. Petkov, */
/*     D.W. Gu and M.M. Konstantinov. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     May 2003, Aug. 2005, March 2008, May 2009, Sep. 2009. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, system norm. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
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

/*     Test the input scalar parameters. */

    /* Parameter adjustments */
    --fpeak;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --gpeak;
    --iwork;
    --dwork;
    --cwork;

    /* Function Body */
    n2 = *n << 1;
    nn = *n * *n;
    pm = *p + *m;
    n2pm = n2 + pm;
    minpm = min(*p,*m);
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    fulle = lsame_(jobe, "G", (ftnlen)1, (ftnlen)1);
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);

    if (! (discr || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (fulle || lsame_(jobe, "I", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lequil || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (! (withd || lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1))) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*p < 0) {
	*info = -7;
    } else if (min(fpeak[1],fpeak[2]) < 0.) {
	*info = -8;
    } else if (*lda < max(1,*n)) {
	*info = -10;
    } else if (*lde < 1 || fulle && *lde < *n) {
	*info = -12;
    } else if (*ldb < max(1,*n)) {
	*info = -14;
    } else if (*ldc < max(1,*p)) {
	*info = -16;
    } else if (*ldd < 1 || withd && *ldd < *p) {
	*info = -18;
    } else if (*tol < 0. || *tol >= 1.) {
	*info = -20;
    } else {
	bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)
		6);
	cnorm = dlange_("1-norm", p, n, &c__[c_offset], ldc, &dwork[1], (
		ftnlen)6);
	nodyn = *n == 0 || min(bnorm,cnorm) == 0.;
	usepen = fulle || discr;

/*        Compute workspace. */

	id = minpm * 6;
/* Computing MAX */
	i__1 = (minpm << 2) + max(*p,*m);
	ic = max(i__1,id);
	if (minpm == 0) {
	    minwrk = 1;
	} else if (nodyn) {
	    if (withd) {
		minwrk = *p * *m + ic;
	    } else {
		minwrk = 1;
	    }
	} else {
	    if (discr) {
		ib = 0;
		ie = id;
	    } else {
		ib = *n * (*n + *m);
		if (! withd) {
		    ib += *p * *m;
		}
		ie = ic;
	    }
	    if (withd) {
		ir = *p * *m;
		if (! usepen) {
		    minwrk = *p * *p + *m * *m;
		    ir += *n * pm;
		} else {
		    minwrk = 0;
		}
		minwrk = minwrk + ir + ic;
		ir += minpm;
	    } else {
		ir = 0;
		minwrk = 0;
	    }
	    ir += *n * (*n + pm);
	    if (fulle) {
		ir += nn;
		if (lequil) {
/* Computing MAX */
		    i__1 = minwrk, i__2 = ir + *n * 9;
		    minwrk = max(i__1,i__2);
		}
/* Computing MAX */
/* Computing MAX */
		i__3 = *m, i__4 = nn << 1, i__3 = max(i__3,i__4), i__4 = *n + 
			ib + ie;
		i__1 = minwrk, i__2 = ir + (*n << 2) + max(i__3,i__4);
		minwrk = max(i__1,i__2);
	    } else {
/* Computing MAX */
/* Computing MAX */
		i__3 = max(*m,*p), i__4 = nn + n2, i__3 = max(i__3,i__4), 
			i__4 = *n * 3 + ib + ie;
		i__1 = minwrk, i__2 = ir + *n + max(i__3,i__4);
		minwrk = max(i__1,i__2);
	    }
	    lw = 0;
	    if (! usepen) {
		lw = ir + (nn << 2) + *n * 11;
		if (withd) {
		    lw = lw + max(*m,*p) + *n * pm;
		}
	    }
	    if (usepen || withd) {
/* Computing MAX */
/* Computing MAX */
		i__3 = n2pm + pm, i__4 = nn + n2 << 3;
		i__1 = lw, i__2 = ir + *n * 6 + n2pm * n2pm + max(i__3,i__4);
		lw = max(i__1,i__2);
	    }
/* Computing MAX */
	    i__1 = max(1,minwrk), i__1 = max(i__1,lw), i__2 = ir + n2 + ie;
	    minwrk = max(i__1,i__2);
	}

	if (*ldwork < minwrk) {
	    *info = -23;
	} else {
	    if (nodyn) {
		mincwr = 1;
	    } else {
/* Computing MAX */
		i__1 = 1, i__2 = (*n + *m) * (*n + *p) + (minpm << 1) + max(*
			p,*m);
		mincwr = max(i__1,i__2);
	    }
	    if (*lcwork < mincwr) {
		*info = -25;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("AB13DD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *p == 0) {
	gpeak[1] = 0.;
	fpeak[1] = 0.;
	gpeak[2] = 1.;
	fpeak[2] = 1.;
	dwork[1] = 1.;
	cwork[1].r = 1., cwork[1].i = 0.;
	return 0;
    }

/*     Determine the maximum singular value of G(infinity) = D . */
/*     If JOBE = 'I' and DICO = 'C', the full SVD of D, D = U*S*V', is */
/*     computed and saved for later use. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    id = 1;
    if (withd) {
	is = id + *p * *m;
	if (usepen || nodyn) {
	    iu = is + minpm;
	    iv = iu;
	    iwrk = iv;
	    *(unsigned char *)vect = 'N';
	} else {
	    ibv = is + minpm;
	    icu = ibv + *n * *m;
	    iu = icu + *p * *n;
	    iv = iu + *p * *p;
	    iwrk = iv + *m * *m;
	    *(unsigned char *)vect = 'A';
	}

/*        Workspace: need   P*M + MIN(P,M) + V + */
/*                          MAX( 3*MIN(P,M) + MAX(P,M), 5*MIN(P,M) ), */
/*                          where V = N*(M+P) + P*P + M*M, */
/*                                        if JOBE = 'I' and DICO = 'C', */
/*                                        and N > 0, B <> 0, C <> 0, */
/*                                V = 0,  otherwise; */
/*                   prefer larger. */

	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)4);
	i__1 = *ldwork - iwrk + 1;
	dgesvd_(vect, vect, p, m, &dwork[id], p, &dwork[is], &dwork[iu], p, &
		dwork[iv], m, &dwork[iwrk], &i__1, &ierr, (ftnlen)1, (ftnlen)
		1);
	if (ierr > 0) {
	    *info = 3;
	    return 0;
	}
	gammal = dwork[is];
	maxwrk = (integer) dwork[iwrk] + iwrk - 1;

/*        Restore D for later calculations. */

	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)4);
    } else {
	iwrk = 1;
	gammal = 0.;
	maxwrk = 1;
    }

/*     Quick return if possible. */

    if (nodyn) {
	gpeak[1] = gammal;
	fpeak[1] = 0.;
	gpeak[2] = 1.;
	fpeak[2] = 1.;
	dwork[1] = (doublereal) maxwrk;
	cwork[1].r = 1., cwork[1].i = 0.;
	return 0;
    }

    if (! usepen && withd) {

/*        Standard continuous-time case, D <> 0: Compute B*V and C'*U . */

	dgemm_("No Transpose", "Transpose", n, m, m, &c_b20, &b[b_offset], 
		ldb, &dwork[iv], m, &c_b21, &dwork[ibv], n, (ftnlen)12, (
		ftnlen)9);
	dgemm_("Transpose", "No Transpose", n, p, p, &c_b20, &c__[c_offset], 
		ldc, &dwork[iu], p, &c_b21, &dwork[icu], n, (ftnlen)9, (
		ftnlen)12);

/*        U and V are no longer needed: free their memory space. */
/*        Total workspace here: need   P*M + MIN(P,M) + N*(M+P) */
/*        (JOBE = 'I', DICO = 'C', JOBD = 'D'). */

	iwrk = iu;
    }

/*     Get machine constants. */

    eps = dlamch_("Epsilon", (ftnlen)7);
    safmin = dlamch_("Safe minimum", (ftnlen)12);
    safmax = 1. / safmin;
    dlabad_(&safmin, &safmax);
    smlnum = sqrt(safmin) / dlamch_("Precision", (ftnlen)9);
    bignum = 1. / smlnum;
    toler = sqrt(eps);

/*     Initiate the transformation of the system to an equivalent one, */
/*     to be used for eigenvalue computations. */

/*     Additional workspace: need   N*N + N*M + P*N + 2*N, if JOBE = 'I'; */
/*                                2*N*N + N*M + P*N + 2*N, if JOBE = 'G'. */

    ia = iwrk;
    ie = ia + nn;
    if (fulle) {
	ib = ie + nn;
    } else {
	ib = ie;
    }
    ic = ib + *n * *m;
    ir = ic + *p * *n;
    ii = ir + *n;
    ibt = ii + *n;

    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ib], n, (ftnlen)4);
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ic], p, (ftnlen)4);

/*     Scale A if maximum element is outside the range [SMLNUM,BIGNUM]. */

    anrm = dlange_("Max", n, n, &dwork[ia], n, &dwork[1], (ftnlen)3);
    ilascl = FALSE_;
    if (anrm > 0. && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = TRUE_;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = TRUE_;
    }
    if (ilascl) {
	dlascl_("General", &c__0, &c__0, &anrm, &anrmto, n, n, &dwork[ia], n, 
		&ierr, (ftnlen)7);
    }

    if (fulle) {

/*        Descriptor system. */

/*        Additional workspace: need   N. */

	iwrk = ibt + *n;
	dlacpy_("Full", n, n, &e[e_offset], lde, &dwork[ie], n, (ftnlen)4);

/*        Scale E if maximum element is outside the range */
/*        [SMLNUM,BIGNUM]. */

	enrm = dlange_("Max", n, n, &dwork[ie], n, &dwork[1], (ftnlen)3);
	ilescl = FALSE_;
	if (enrm > 0. && enrm < smlnum) {
	    enrmto = smlnum;
	    ilescl = TRUE_;
	} else if (enrm > bignum) {
	    enrmto = bignum;
	    ilescl = TRUE_;
	} else if (enrm == 0.) {

/*           Error return: Matrix E is 0. */

	    *info = 1;
	    return 0;
	}
	if (ilescl) {
	    dlascl_("General", &c__0, &c__0, &enrm, &enrmto, n, n, &dwork[ie],
		     n, &ierr, (ftnlen)7);
	}

/*        Equilibrate the system, if required. */

/*        Additional workspace: need   6*N. */

	if (lequil) {
	    tg01ad_("All", n, n, m, p, &c_b21, &dwork[ia], n, &dwork[ie], n, &
		    dwork[ib], n, &dwork[ic], p, &dwork[ii], &dwork[ir], &
		    dwork[iwrk], &ierr, (ftnlen)3);
	}

/*        For efficiency of later calculations, the system (A,E,B,C) is */
/*        reduced to an equivalent one with the state matrix A in */
/*        Hessenberg form, and E upper triangular. */
/*        First, permute (A,E) to make it more nearly triangular. */

	dggbal_("Permute", n, &dwork[ia], n, &dwork[ie], n, &ilo, &ihi, &
		dwork[ii], &dwork[ir], &dwork[iwrk], &ierr, (ftnlen)7);

/*        Apply the permutations to (the copies of) B and C. */

	i__1 = ihi + 1;
	for (i__ = *n; i__ >= i__1; --i__) {
	    k = (integer) dwork[ii + i__ - 1];
	    if (k != i__) {
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
	    }
	    k = (integer) dwork[ir + i__ - 1];
	    if (k != i__) {
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
	    }
/* L10: */
	}

	i__1 = ilo - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = (integer) dwork[ii + i__ - 1];
	    if (k != i__) {
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
	    }
	    k = (integer) dwork[ir + i__ - 1];
	    if (k != i__) {
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
	    }
/* L20: */
	}

/*        Reduce (A,E) to generalized Hessenberg form and apply the */
/*        transformations to B and C. */
/*        Additional workspace: need   N + MAX(N,M); */
/*                              prefer N + MAX(N,M)*NB. */

	i__1 = *ldwork - iwrk + 1;
	tg01bd_("General", "No Q", "No Z", n, m, p, &ilo, &ihi, &dwork[ia], n,
		 &dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[1], &
		c__1, &dwork[1], &c__1, &dwork[iwrk], &i__1, &ierr, (ftnlen)7,
		 (ftnlen)4, (ftnlen)4);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*        Check whether matrix E is nonsingular. */
/*        Additional workspace: need   3*N. */

	dtrcon_("1-norm", "Upper", "Non Unit", n, &dwork[ie], n, &rcond, &
		dwork[iwrk], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (ftnlen)
		8);
	if (rcond <= (doublereal) (*n) * 10. * eps) {

/*           Error return: Matrix E is numerically singular. */

	    *info = 1;
	    return 0;
	}

/*        Perform QZ algorithm, computing eigenvalues. The generalized */
/*        Hessenberg form is saved for later use. */
/*        Additional workspace: need   2*N*N + N; */
/*                              prefer larger. */

	ias = iwrk;
	ies = ias + nn;
	iwrk = ies + nn;
	dlacpy_("Full", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)4);
	dlacpy_("Full", n, n, &dwork[ie], n, &dwork[ies], n, (ftnlen)4);
	i__1 = *ldwork - iwrk + 1;
	dhgeqz_("Eigenvalues", "No Vectors", "No Vectors", n, &ilo, &ihi, &
		dwork[ias], n, &dwork[ies], n, &dwork[ir], &dwork[ii], &dwork[
		ibt], &dwork[1], n, &dwork[1], n, &dwork[iwrk], &i__1, &ierr, 
		(ftnlen)11, (ftnlen)10, (ftnlen)10);
	if (ierr != 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*        Check if unscaling would cause over/underflow; if so, rescale */
/*        eigenvalues (DWORK( IR+I-1 ),DWORK( II+I-1 ),DWORK( IBT+I-1 )) */
/*        so DWORK( IBT+I-1 ) is on the order of E(I,I) and */
/*        DWORK( IR+I-1 ) and DWORK( II+I-1 ) are on the order of A(I,I). */

	if (ilascl) {

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (dwork[ii + i__ - 1] != 0.) {
		    if (dwork[ir + i__ - 1] / safmax > anrmto / anrm || 
			    safmin / dwork[ir + i__ - 1] > anrm / anrmto) {
			tm = (d__1 = dwork[ia + (i__ - 1) * *n + i__] / dwork[
				ir + i__ - 1], abs(d__1));
			dwork[ibt + i__ - 1] *= tm;
			dwork[ir + i__ - 1] *= tm;
			dwork[ii + i__ - 1] *= tm;
		    } else if (dwork[ii + i__ - 1] / safmax > anrmto / anrm ||
			     safmin / dwork[ii + i__ - 1] > anrm / anrmto) {
			tm = (d__1 = dwork[ia + i__ * *n + i__] / dwork[ii + 
				i__ - 1], abs(d__1));
			dwork[ibt + i__ - 1] *= tm;
			dwork[ir + i__ - 1] *= tm;
			dwork[ii + i__ - 1] *= tm;
		    }
		}
/* L30: */
	    }

	}

	if (ilescl) {

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (dwork[ii + i__ - 1] != 0.) {
		    if (dwork[ibt + i__ - 1] / safmax > enrmto / enrm || 
			    safmin / dwork[ibt + i__ - 1] > enrm / enrmto) {
			tm = (d__1 = dwork[ie + (i__ - 1) * *n + i__] / dwork[
				ibt + i__ - 1], abs(d__1));
			dwork[ibt + i__ - 1] *= tm;
			dwork[ir + i__ - 1] *= tm;
			dwork[ii + i__ - 1] *= tm;
		    }
		}
/* L40: */
	    }

	}

/*        Undo scaling. */

	if (ilascl) {
	    dlascl_("Hessenberg", &c__0, &c__0, &anrmto, &anrm, n, n, &dwork[
		    ia], n, &ierr, (ftnlen)10);
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ir], n, &ierr, (ftnlen)7);
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ii], n, &ierr, (ftnlen)7);
	}

	if (ilescl) {
	    dlascl_("Upper", &c__0, &c__0, &enrmto, &enrm, n, n, &dwork[ie], 
		    n, &ierr, (ftnlen)5);
	    dlascl_("General", &c__0, &c__0, &enrmto, &enrm, n, &c__1, &dwork[
		    ibt], n, &ierr, (ftnlen)7);
	}

    } else {

/*        Standard state-space system. */

	if (lequil) {

/*           Equilibrate the system. */

	    maxred = 100.;
	    tb01id_("All", n, m, p, &maxred, &dwork[ia], n, &dwork[ib], n, &
		    dwork[ic], p, &dwork[ii], &ierr, (ftnlen)3);
	}

/*        For efficiency of later calculations, the system (A,B,C) is */
/*        reduced to a similar one with the state matrix in Hessenberg */
/*        form. */

/*        First, permute the matrix A to make it more nearly triangular */
/*        and apply the permutations to B and C. */

	dgebal_("Permute", n, &dwork[ia], n, &ilo, &ihi, &dwork[ir], &ierr, (
		ftnlen)7);

	i__1 = ihi + 1;
	for (i__ = *n; i__ >= i__1; --i__) {
	    k = (integer) dwork[ir + i__ - 1];
	    if (k != i__) {
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
	    }
/* L50: */
	}

	i__1 = ilo - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = (integer) dwork[ir + i__ - 1];
	    if (k != i__) {
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
	    }
/* L60: */
	}

/*        Reduce A to upper Hessenberg form and apply the transformations */
/*        to B and C. */
/*        Additional workspace: need   N;   (from II) */
/*                              prefer N*NB. */

	itau = ir;
	iwrk = itau + *n;
	i__1 = *ldwork - iwrk + 1;
	dgehrd_(n, &ilo, &ihi, &dwork[ia], n, &dwork[itau], &dwork[iwrk], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*        Additional workspace: need   M; */
/*                              prefer M*NB. */

	i__1 = *ldwork - iwrk + 1;
	dormhr_("Left", "Transpose", n, m, &ilo, &ihi, &dwork[ia], n, &dwork[
		itau], &dwork[ib], n, &dwork[iwrk], &i__1, &ierr, (ftnlen)4, (
		ftnlen)9);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*        Additional workspace: need   P; */
/*                              prefer P*NB. */

	i__1 = *ldwork - iwrk + 1;
	dormhr_("Right", "NoTranspose", p, n, &ilo, &ihi, &dwork[ia], n, &
		dwork[itau], &dwork[ic], p, &dwork[iwrk], &i__1, &ierr, (
		ftnlen)5, (ftnlen)11);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*        Compute the eigenvalues. The Hessenberg form is saved for */
/*        later use. */
/*        Additional workspace:  need   N*N + N;   (from IBT) */
/*                               prefer larger. */

	ias = ibt;
	iwrk = ias + nn;
	dlacpy_("Full", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)4);
	i__1 = *ldwork - iwrk + 1;
	dhseqr_("Eigenvalues", "No Vectors", n, &ilo, &ihi, &dwork[ias], n, &
		dwork[ir], &dwork[ii], &dwork[1], n, &dwork[iwrk], &i__1, &
		ierr, (ftnlen)11, (ftnlen)10);
	if (ierr > 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

	if (ilascl) {

/*           Undo scaling for the Hessenberg form of A and eigenvalues. */

	    dlascl_("Hessenberg", &c__0, &c__0, &anrmto, &anrm, n, n, &dwork[
		    ia], n, &ierr, (ftnlen)10);
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ir], n, &ierr, (ftnlen)7);
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ii], n, &ierr, (ftnlen)7);
	}

    }

/*     Look for (generalized) eigenvalues on the boundary of the */
/*     stability domain. (Their existence implies an infinite norm.) */
/*     Additional workspace:  need   2*N.   (from IAS) */

    im = ias;
    iar = im + *n;
    imin = ii;
    wrmin = safmax;
    bound = eps * 1e3;

    if (discr) {
	gammal = 0.;

/*        For discrete-time case, compute the logarithm of the non-zero */
/*        eigenvalues and save their moduli and absolute real parts. */
/*        (The logarithms are overwritten on the eigenvalues.) */
/*        Also, find the minimum distance to the unit circle. */

	if (fulle) {

	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
			safmax * dwork[ibt + i__]) {
		    tm /= dwork[ibt + i__];
		} else {

/*                 The pencil has too large eigenvalues. SAFMAX is used. */

		    tm = safmax;
		}
		if (tm != 0.) {
		    dwork[ii + i__] = atan2(dwork[ii + i__], dwork[ir + i__]);
		    dwork[ir + i__] = log(tm);
		}
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		tm = (d__1 = 1. - tm, abs(d__1));
		if (tm < wrmin) {
		    imin = ii + i__;
		    wrmin = tm;
		}
		++im;
		dwork[iar + i__] = (d__1 = dwork[ir + i__], abs(d__1));
/* L70: */
	    }

	} else {

	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		if (tm != 0.) {
		    dwork[ii + i__] = atan2(dwork[ii + i__], dwork[ir + i__]);
		    dwork[ir + i__] = log(tm);
		}
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		tm = (d__1 = 1. - tm, abs(d__1));
		if (tm < wrmin) {
		    imin = ii + i__;
		    wrmin = tm;
		}
		++im;
		dwork[iar + i__] = (d__1 = dwork[ir + i__], abs(d__1));
/* L80: */
	    }

	}

    } else {

/*        For continuous-time case, save moduli of eigenvalues and */
/*        absolute real parts and find the maximum modulus and minimum */
/*        absolute real part. */

	wmax = 0.;

	if (fulle) {

	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		tm = (d__1 = dwork[ir + i__], abs(d__1));
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && dwork[
			im] < safmax * dwork[ibt + i__]) {
		    tm /= dwork[ibt + i__];
		    dwork[im] /= dwork[ibt + i__];
		} else {
		    if (tm < safmax * dwork[ibt + i__]) {
			tm /= dwork[ibt + i__];
		    } else {

/*                    The pencil has too large eigenvalues. */
/*                    SAFMAX is used. */

			tm = safmax;
		    }
		    dwork[im] = safmax;
		}
		if (tm < wrmin) {
		    imin = ii + i__;
		    wrmin = tm;
		}
		dwork[iar + i__] = tm;
		if (dwork[im] > wmax) {
		    wmax = dwork[im];
		}
		++im;
/* L90: */
	    }

	} else {

	    i__1 = *n - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
		tm = (d__1 = dwork[ir + i__], abs(d__1));
		if (tm < wrmin) {
		    imin = ii + i__;
		    wrmin = tm;
		}
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
		if (dwork[im] > wmax) {
		    wmax = dwork[im];
		}
		++im;
		dwork[iar + i__] = tm;
/* L100: */
	    }

	}

	bound += eps * wmax;

    }

    im -= *n;

    if (wrmin < bound) {

/*        The L-infinity norm was found as infinite. */

	gpeak[1] = 1.;
	gpeak[2] = 0.;
	tm = (d__1 = dwork[imin], abs(d__1));
	if (discr) {
	    tm = (d__1 = atan2(sin(tm), cos(tm)), abs(d__1));
	}
	fpeak[1] = tm;
	if (tm < safmax) {
	    fpeak[2] = 1.;
	} else {
	    fpeak[2] = 0.;
	}

	dwork[1] = (doublereal) maxwrk;
	cwork[1].r = 1., cwork[1].i = 0.;
	return 0;
    }

/*     Determine the maximum singular value of */
/*        G(lambda) = C*inv(lambda*E - A)*B + D, */
/*     over a selected set of frequencies. Besides the frequencies w = 0, */
/*     w = pi (if DICO = 'D'), and the given value FPEAK, this test set */
/*     contains the peak frequency for each mode (or an approximation */
/*     of it). The (generalized) Hessenberg form of the system is used. */

/*     First, determine the maximum singular value of G(0) and set FPEAK */
/*     accordingly. */
/*     Additional workspace: */
/*           complex: need   1, if DICO = 'C'; */
/*                           (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)), otherwise; */
/*                    prefer larger; */
/*           real:    need   LDW0+LDW1+LDW2, where */
/*                           LDW0 = N*N+N*M, if DICO = 'C'; */
/*                           LDW0 = 0,       if DICO = 'D'; */
/*                           LDW1 = P*M,     if DICO = 'C', JOBD = 'Z'; */
/*                           LDW1 = 0,       otherwise; */
/*                           LDW2 = MIN(P,M)+MAX(3*MIN(P,M)+MAX(P,M), */
/*                                               5*MIN(P,M)), */
/*                                              if DICO = 'C'; */
/*                           LDW2 = 6*MIN(P,M), otherwise. */
/*                    prefer larger. */

    if (discr) {
	ias = ia;
	ibs = ib;
	iwrk = iar + *n;
    } else {
	ias = iar + *n;
	ibs = ias + nn;
	iwrk = ibs + *n * *m;
	dlacpy_("Upper", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)5);
	i__1 = *n - 1;
	i__2 = *n + 1;
	i__3 = *n + 1;
	dcopy_(&i__1, &dwork[ia + 1], &i__2, &dwork[ias + 1], &i__3);
	dlacpy_("Full", n, m, &dwork[ib], n, &dwork[ibs], n, (ftnlen)4);
    }
    i__1 = *ldwork - iwrk + 1;
    gamma = ab13dx_(dico, jobe, jobd, n, m, p, &c_b21, &dwork[ias], n, &dwork[
	    ie], n, &dwork[ibs], n, &dwork[ic], p, &dwork[id], p, &iwork[1], &
	    dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (ftnlen)1, (ftnlen)
	    1, (ftnlen)1);
/* Computing MAX */
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
    maxwrk = max(i__1,maxwrk);
    if (ierr >= 1 && ierr <= *n) {
	gpeak[1] = 1.;
	fpeak[1] = 0.;
	gpeak[2] = 0.;
	fpeak[2] = 1.;
	goto L340;
    } else if (ierr == *n + 1) {
	*info = 3;
	return 0;
    }

    fpeaks = fpeak[1];
    fpeaki = fpeak[2];
    if (gammal < gamma) {
	gammal = gamma;
	fpeak[1] = 0.;
	fpeak[2] = 1.;
    } else if (! discr) {
	fpeak[1] = 1.;
	fpeak[2] = 0.;
    }

    maxcwk = (integer) cwork[1].r;

    if (discr) {

/*        Try the frequency w = pi. */

	pi = atan(1.) * 4.;
	i__1 = *ldwork - iwrk + 1;
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &pi, &dwork[ia], n, &dwork[
		ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &iwork[1]
		, &dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/* Computing MAX */
	i__1 = (integer) cwork[1].r;
	maxcwk = max(i__1,maxcwk);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);
	if (ierr >= 1 && ierr <= *n) {
	    gpeak[1] = 1.;
	    fpeak[1] = pi;
	    gpeak[2] = 0.;
	    fpeak[2] = 1.;
	    goto L340;
	} else if (ierr == *n + 1) {
	    *info = 3;
	    return 0;
	}

	if (gammal < gamma) {
	    gammal = gamma;
	    fpeak[1] = pi;
	    fpeak[2] = 1.;
	}

    } else {
	iwrk = ias;

/*        Restore D, if needed. */

	if (withd) {
	    dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)
		    4);
	}
    }

/*     Build the remaining set of frequencies. */
/*     Complex workspace:  need   (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)); */
/*                         prefer larger. */
/*     Real workspace:     need   LDW2, see above; */
/*                         prefer larger. */

    if (min(fpeaks,fpeaki) != 0.) {

/*        Compute also the norm at the given (finite) frequency. */

	i__1 = *ldwork - iwrk + 1;
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &fpeaks, &dwork[ia], n, &
		dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &
		iwork[1], &dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
	i__1 = (integer) cwork[1].r;
	maxcwk = max(i__1,maxcwk);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);
	if (discr) {
	    tm = (d__1 = atan2(sin(fpeaks), cos(fpeaks)), abs(d__1));
	} else {
	    tm = fpeaks;
	}
	if (ierr >= 1 && ierr <= *n) {
	    gpeak[1] = 1.;
	    fpeak[1] = tm;
	    gpeak[2] = 0.;
	    fpeak[2] = 1.;
	    goto L340;
	} else if (ierr == *n + 1) {
	    *info = 3;
	    return 0;
	}

	if (gammal < gamma) {
	    gammal = gamma;
	    fpeak[1] = tm;
	    fpeak[2] = 1.;
	}

    }

    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (dwork[ii + i__] >= 0. && dwork[im + i__] > 0.) {
	    if (dwork[im + i__] >= 1. || dwork[im + i__] < 1. && dwork[iar + 
		    i__] < safmax * dwork[im + i__]) {
		rat = dwork[iar + i__] / dwork[im + i__];
	    } else {
		rat = 1.;
	    }
/* Computing MAX */
/* Computing 2nd power */
	    d__3 = rat;
	    d__1 = .25, d__2 = 1. - d__3 * d__3 * 2.;
	    omega = dwork[im + i__] * sqrt((max(d__1,d__2)));

	    i__2 = *ldwork - iwrk + 1;
	    gamma = ab13dx_(dico, jobe, jobd, n, m, p, &omega, &dwork[ia], n, 
		    &dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], 
		    p, &iwork[1], &dwork[iwrk], &i__2, &cwork[1], lcwork, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
	    i__2 = (integer) cwork[1].r;
	    maxcwk = max(i__2,maxcwk);
/* Computing MAX */
	    i__2 = (integer) dwork[iwrk] + iwrk - 1;
	    maxwrk = max(i__2,maxwrk);
	    if (discr) {
		tm = (d__1 = atan2(sin(omega), cos(omega)), abs(d__1));
	    } else {
		tm = omega;
	    }
	    if (ierr >= 1 && ierr <= *n) {
		gpeak[1] = 1.;
		fpeak[1] = tm;
		gpeak[2] = 0.;
		fpeak[2] = 1.;
		goto L340;
	    } else if (ierr == *n + 1) {
		*info = 3;
		return 0;
	    }

	    if (gammal < gamma) {
		gammal = gamma;
		fpeak[1] = tm;
		fpeak[2] = 1.;
	    }

	}
/* L110: */
    }

/*     Return if the lower bound is zero. */

    if (gammal == 0.) {
	gpeak[1] = 0.;
	fpeak[1] = 0.;
	gpeak[2] = 1.;
	fpeak[2] = 1.;
	goto L340;
    }

/*     Start the modified gamma iteration for the Bruinsma-Steinbuch */
/*     algorithm. */

    if (! discr) {
	rtol = toler * 100.;
    }
    iter = 0;

/*     WHILE ( Iteration may continue ) DO */

L120:

    ++iter;
    gamma = (*tol + 1.) * gammal;
    usepen = fulle || discr;
    if (! usepen && withd) {

/*           Check whether one can use an explicit Hamiltonian matrix: */
/*           compute */
/*           min(rcond(GAMMA**2*Im - S'*S), rcond(GAMMA**2*Ip - S*S')). */
/*           If P = M = 1, then GAMMA**2 - S(1)**2 is used instead. */

	if (*m != *p) {
/* Computing 2nd power */
	    d__1 = dwork[is] / gamma;
	    rcond = 1. - d__1 * d__1;
	} else if (minpm > 1) {
/* Computing 2nd power */
	    d__1 = gamma;
/* Computing 2nd power */
	    d__2 = dwork[is];
/* Computing 2nd power */
	    d__3 = gamma;
/* Computing 2nd power */
	    d__4 = dwork[is + *p - 1];
	    rcond = (d__1 * d__1 - d__2 * d__2) / (d__3 * d__3 - d__4 * d__4);
	} else {
/* Computing 2nd power */
	    d__1 = gamma;
/* Computing 2nd power */
	    d__2 = dwork[is];
	    rcond = d__1 * d__1 - d__2 * d__2;
	}

	usepen = rcond < rtol;
    }

    if (usepen) {

/*           Use the QZ algorithm on a pencil. */
/*           Additional workspace here:  need   6*N.   (from IR) */

	ii = ir + n2;
	ibt = ii + n2;
	ih12 = ibt + n2;
	im = ih12;

/*           Set up the needed parts of the Hamiltonian pencil (H,J), */

/*                  ( H11  H12 ) */
/*              H = (          ) , */
/*                  ( H21  H22 ) */

/*           with */

/*                 ( A  0  )            ( 0  B )            ( E  0  ) */
/*           H11 = (       ),     H12 = (      )/nB,  J11 = (       ), */
/*                 ( 0 -A' )            ( C' 0 )            ( 0  E' ) */

/*                 ( C  0  )            ( Ip  D/g ) */
/*           H21 = (       )*nB,  H22 = (         ), */
/*                 ( 0 -B' )            ( D'/g Im ) */

/*           if DICO = 'C', and */

/*                 ( A  0  )            ( B  0  )            ( E  0 ) */
/*           H11 = (       ),     H12 = (       )/nB,  J11 = (      ), */
/*                 ( 0  E' )            ( 0  C' )            ( 0  A') */

/*                 ( 0  0  )            ( Im  D'/g )         ( 0  B') */
/*           H21 = (       )*nB,  H22 = (          ),  J21 = (      )*nB, */
/*                 ( C  0  )            ( D/g  Ip  )         ( 0  0 ) */

/*           if DICO = 'D', where g = GAMMA, and nB = norm(B,1). */
/*           First build [H12; H22]. */

	temp[0] = 0.;
	ih = ih12;

	if (discr) {

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = b[i__ + j * b_dim1] / bnorm;
		    ++ih;
/* L130: */
		}

		i__2 = *n + *m;
		dcopy_(&i__2, temp, &c__0, &dwork[ih], &c__1);
		dwork[ih + *n + j - 1] = 1.;
		ih = ih + *n + *m;

		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = d__[i__ + j * d_dim1] / gamma;
		    ++ih;
/* L140: */
		}

/* L150: */
	    }

	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
		ih += *n;

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = c__[j + i__ * c_dim1] / bnorm;
		    ++ih;
/* L160: */
		}

		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = d__[j + i__ * d_dim1] / gamma;
		    ++ih;
/* L170: */
		}

		dcopy_(p, temp, &c__0, &dwork[ih], &c__1);
		dwork[ih + j - 1] = 1.;
		ih += *p;
/* L180: */
	    }

	} else {

	    i__1 = *p;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
		ih += *n;

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = c__[j + i__ * c_dim1] / bnorm;
		    ++ih;
/* L190: */
		}

		dcopy_(p, temp, &c__0, &dwork[ih], &c__1);
		dwork[ih + j - 1] = 1.;
		ih += *p;

		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = d__[j + i__ * d_dim1] / gamma;
		    ++ih;
/* L200: */
		}

/* L210: */
	    }

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = b[i__ + j * b_dim1] / bnorm;
		    ++ih;
/* L220: */
		}

		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
		ih += *n;

		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[ih] = d__[i__ + j * d_dim1] / gamma;
		    ++ih;
/* L230: */
		}

		dcopy_(m, temp, &c__0, &dwork[ih], &c__1);
		dwork[ih + j - 1] = 1.;
		ih += *m;
/* L240: */
	    }

	}

/*           Compute the QR factorization of [H12; H22]. */
/*           For large P and M, it could be more efficient to exploit the */
/*           structure of [H12; H22] and use the factored form of Q. */
/*           Additional workspace: need   (2*N+P+M)*(2*N+P+M)+2*(P+M); */
/*                                 prefer (2*N+P+M)*(2*N+P+M)+P+M+ */
/*                                                           (P+M)*NB. */

	itau = ih12 + n2pm * n2pm;
	iwrk = itau + pm;
	i__1 = *ldwork - iwrk + 1;
	dgeqrf_(&n2pm, &pm, &dwork[ih12], &n2pm, &dwork[itau], &dwork[iwrk], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*           Apply part of the orthogonal transformation: */
/*           Q1 = Q(:,P+M+(1:2*N))' to the matrix [H11; H21/GAMMA]. */
/*           If DICO = 'C', apply Q(1:2*N,P+M+(1:2*N))' to the */
/*           matrix J11. */
/*           If DICO = 'D', apply Q1 to the matrix [J11; J21/GAMMA]. */
/*           H11, H21, J11, and J21 are not fully built. */
/*           First, build the (2*N+P+M)-by-(2*N+P+M) matrix Q. */
/*           Using Q will often provide better efficiency than the direct */
/*           use of the factored form of Q, especially when P+M < N. */
/*           Additional workspace: need   P+M+2*N+P+M; */
/*                                 prefer P+M+(2*N+P+M)*NB. */

	i__1 = *ldwork - iwrk + 1;
	dorgqr_(&n2pm, &n2pm, &pm, &dwork[ih12], &n2pm, &dwork[itau], &dwork[
		iwrk], &i__1, &ierr);
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

/*           Additional workspace: need   8*N*N. */

	ipa = itau;
	ipe = ipa + (nn << 2);
	iwrk = ipe + (nn << 2);
	dgemm_("Transpose", "No Transpose", &n2, n, n, &c_b20, &dwork[ih12 + 
		pm * n2pm], &n2pm, &a[a_offset], lda, &c_b21, &dwork[ipa], &
		n2, (ftnlen)9, (ftnlen)12);
	if (discr) {
	    d__1 = bnorm / gamma;
	    dgemm_("Transpose", "No Transpose", &n2, n, p, &d__1, &dwork[ih12 
		    + pm * n2pm + n2 + *m], &n2pm, &c__[c_offset], ldc, &
		    c_b20, &dwork[ipa], &n2, (ftnlen)9, (ftnlen)12);
	    if (fulle) {
		dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[
			ih12 + pm * n2pm + *n], &n2pm, &e[e_offset], lde, &
			c_b21, &dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (
			ftnlen)9);
	    } else {
		ma02ad_("Full", n, &n2, &dwork[ih12 + pm * n2pm + *n], &n2pm, 
			&dwork[ipa + (nn << 1)], &n2, (ftnlen)4);
		ny = *n;
	    }
	} else {
	    d__1 = bnorm / gamma;
	    dgemm_("Transpose", "No Transpose", &n2, n, p, &d__1, &dwork[ih12 
		    + pm * n2pm + n2], &n2pm, &c__[c_offset], ldc, &c_b20, &
		    dwork[ipa], &n2, (ftnlen)9, (ftnlen)12);
	    dgemm_("Transpose", "Transpose", &n2, n, n, &c_b163, &dwork[ih12 
		    + pm * n2pm + *n], &n2pm, &a[a_offset], lda, &c_b21, &
		    dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
	    d__1 = -bnorm / gamma;
	    dgemm_("Transpose", "Transpose", &n2, n, m, &d__1, &dwork[ih12 + 
		    pm * n2pm + n2 + *p], &n2pm, &b[b_offset], ldb, &c_b20, &
		    dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
	    ny = n2;
	}

	if (fulle) {
	    dgemm_("Transpose", "No Transpose", &n2, n, n, &c_b20, &dwork[
		    ih12 + pm * n2pm], &n2pm, &e[e_offset], lde, &c_b21, &
		    dwork[ipe], &n2, (ftnlen)9, (ftnlen)12);
	} else {
	    ma02ad_("Full", &ny, &n2, &dwork[ih12 + pm * n2pm], &n2pm, &dwork[
		    ipe], &n2, (ftnlen)4);
	}
	if (discr) {
	    dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[ih12 + 
		    pm * n2pm + *n], &n2pm, &a[a_offset], lda, &c_b21, &dwork[
		    ipe + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
	    d__1 = bnorm / gamma;
	    dgemm_("Transpose", "Transpose", &n2, n, m, &d__1, &dwork[ih12 + 
		    pm * n2pm + n2], &n2pm, &b[b_offset], ldb, &c_b20, &dwork[
		    ipe + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
	} else {
	    if (fulle) {
		dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[
			ih12 + pm * n2pm + *n], &n2pm, &e[e_offset], lde, &
			c_b21, &dwork[ipe + (nn << 1)], &n2, (ftnlen)9, (
			ftnlen)9);
	    }
	}

/*           Compute the eigenvalues of the Hamiltonian pencil. */
/*           Additional workspace: need   16*N; */
/*                                 prefer larger. */

	i__1 = *ldwork - iwrk + 1;
	dggev_("No Vectors", "No Vectors", &n2, &dwork[ipa], &n2, &dwork[ipe],
		 &n2, &dwork[ir], &dwork[ii], &dwork[ibt], &dwork[1], &n2, &
		dwork[1], &n2, &dwork[iwrk], &i__1, &ierr, (ftnlen)10, (
		ftnlen)10);
	if (ierr > 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__1,maxwrk);

    } else if (! withd) {

/*           Standard continuous-time case with D = 0. */
/*           Form the needed part of the Hamiltonian matrix explicitly: */
/*              H = H11 - H12*inv(H22)*H21/g. */
/*           Additional workspace: need   2*N*N+N.   (from IBT) */

	ih = ibt;
	ih12 = ih + nn;
	isl = ih12 + nn + *n;
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ih], n, (ftnlen)4);

/*           Compute triangles of -C'*C/GAMMA and B*B'/GAMMA. */

	d__1 = -1. / gamma;
	dsyrk_("Lower", "Transpose", n, p, &d__1, &c__[c_offset], ldc, &c_b21,
		 &dwork[ih12], n, (ftnlen)5, (ftnlen)9);
	d__1 = 1. / gamma;
	dsyrk_("Upper", "No Transpose", n, m, &d__1, &b[b_offset], ldb, &
		c_b21, &dwork[ih12 + *n], n, (ftnlen)5, (ftnlen)12);

    } else {

/*           Standard continuous-time case with D <> 0 and the SVD of D */
/*           can be used. Compute explicitly the needed part of the */
/*           Hamiltonian matrix: */

/*               (A+B1*S'*inv(g^2*Ip-S*S')*C1' g*B1*inv(g^2*Im-S'*S)*B1') */
/*           H = (                                                      ) */
/*               (  -g*C1*inv(g^2*Ip-S*S')*C1'            -H11'         ) */

/*           where g = GAMMA, B1 = B*V, C1 = C'*U, and H11 is the first */
/*           block of H. */
/*           Primary additional workspace: need   2*N*N+N   (from IBT) */
/*           (for building the relevant part of the Hamiltonian matrix). */

/*           Compute C1*sqrt(inv(g^2*Ip-S*S')) . */
/*           Additional workspace: need   MAX(M,P)+N*P. */

	ih = ibt;
	ih12 = ih + nn;
	isl = ih12 + nn + *n;

	i__1 = minpm - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = gamma;
/* Computing 2nd power */
	    d__2 = dwork[is + i__];
	    dwork[isl + i__] = 1. / sqrt(d__1 * d__1 - d__2 * d__2);
/* L250: */
	}

	if (*m < *p) {
	    dwork[isl + *m] = 1. / gamma;
	    i__1 = *p - *m - 1;
	    dcopy_(&i__1, &dwork[isl + *m], &c__0, &dwork[isl + *m + 1], &
		    c__1);
	}
	isc = isl + max(*m,*p);
	dlacpy_("Full", n, p, &dwork[icu], n, &dwork[isc], n, (ftnlen)4);
	mb01sd_("Column", n, p, &dwork[isc], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute B1*S' . */
/*           Additional workspace: need   N*M. */

	isb = isc + *p * *n;
	dlacpy_("Full", n, m, &dwork[ibv], n, &dwork[isb], n, (ftnlen)4);
	mb01sd_("Column", n, &minpm, &dwork[isb], n, &dwork[1], &dwork[is], (
		ftnlen)6);

/*           Compute B1*S'*sqrt(inv(g^2*Ip-S*S')) . */

	mb01sd_("Column", n, &minpm, &dwork[isb], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute H11 . */

	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ih], n, (ftnlen)4);
	dgemm_("No Transpose", "Transpose", n, n, &minpm, &c_b20, &dwork[isb],
		 n, &dwork[isc], n, &c_b20, &dwork[ih], n, (ftnlen)12, (
		ftnlen)9);

/*           Compute B1*sqrt(inv(g^2*Im-S'*S)) . */

	if (*p < *m) {
	    dwork[isl + *p] = 1. / gamma;
	    i__1 = *m - *p - 1;
	    dcopy_(&i__1, &dwork[isl + *p], &c__0, &dwork[isl + *p + 1], &
		    c__1);
	}
	dlacpy_("Full", n, m, &dwork[ibv], n, &dwork[isb], n, (ftnlen)4);
	mb01sd_("Column", n, m, &dwork[isb], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute the lower triangle of H21 and the upper triangle */
/*           of H12. */

	d__1 = -gamma;
	dsyrk_("Lower", "No Transpose", n, p, &d__1, &dwork[isc], n, &c_b21, &
		dwork[ih12], n, (ftnlen)5, (ftnlen)12);
	dsyrk_("Upper", "No Transpose", n, m, &gamma, &dwork[isb], n, &c_b21, 
		&dwork[ih12 + *n], n, (ftnlen)5, (ftnlen)12);
    }

    if (! usepen) {

/*           Compute the eigenvalues of the Hamiltonian matrix by the */
/*           symplectic URV and the periodic Schur decompositions. */
/*           Additional workspace: need   (2*N+8)*N; */
/*                                 prefer larger. */

	iwrk = isl + nn;
	i__1 = *ldwork - iwrk - *n + 1;
	mb03xd_("Both", "Eigenvalues", "No vectors", "No vectors", n, &dwork[
		ih], n, &dwork[ih12], n, &dwork[isl], n, temp, &c__1, temp, &
		c__1, temp, &c__1, temp, &c__1, &dwork[ir], &dwork[ii], &ilo, 
		&dwork[iwrk], &dwork[iwrk + *n], &i__1, &ierr, (ftnlen)4, (
		ftnlen)11, (ftnlen)10, (ftnlen)10);
	if (ierr > 0) {
	    *info = 2;
	    return 0;
	}
/* Computing MAX */
	i__1 = (integer) dwork[iwrk] + iwrk + *n - 1;
	maxwrk = max(i__1,maxwrk);
    }

/*        Detect eigenvalues on the boundary of the stability domain, */
/*        if any. The test is based on a round-off level of eps*rho(H) */
/*        (after balancing) resulting in worst-case perturbations of */
/*        order sqrt(eps*rho(H)), for continuous-time systems, on the */
/*        real part of poles of multiplicity two (typical as GAMMA */
/*        approaches the infinity norm). Similarly, in the discrete-time */
/*        case. Above, rho(H) is the maximum modulus of eigenvalues */
/*        (continuous-time case). */

/*        Compute maximum eigenvalue modulus and check the absolute real */
/*        parts (if DICO = 'C'), or moduli (if DICO = 'D'). */

    wmax = 0.;

    if (usepen) {

/*           Additional workspace: need   2*N, if DICO = 'D';   (from IM) */
/*                                        0,   if DICO = 'C'. */

	i__1 = n2 - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
	    if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
		    safmax * dwork[ibt + i__]) {
		tm /= dwork[ibt + i__];
	    } else {

/*                 The pencil has too large eigenvalues. SAFMAX is used. */

		tm = safmax;
	    }
	    wmax = max(wmax,tm);
	    if (discr) {
		dwork[im + i__] = tm;
	    }
/* L260: */
	}

    } else {

	i__1 = *n - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
	    wmax = max(wmax,tm);
/* L270: */
	}

    }

    nei = 0;

    if (usepen) {

	i__1 = n2 - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    if (discr) {
		tm = (d__1 = 1. - dwork[im + i__], abs(d__1));
	    } else {
		tm = (d__1 = dwork[ir + i__], abs(d__1));
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
			safmax * dwork[ibt + i__]) {
		    tm /= dwork[ibt + i__];
		} else {

/*                    The pencil has too large eigenvalues. */
/*                    SAFMAX is used. */

		    tm = safmax;
		}
	    }
	    if (tm <= toler * sqrt(wmax + 100.)) {
		dwork[ir + nei] = dwork[ir + i__] / dwork[ibt + i__];
		dwork[ii + nei] = dwork[ii + i__] / dwork[ibt + i__];
		++nei;
	    }
/* L280: */
	}

    } else {

	i__1 = *n - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    tm = (d__1 = dwork[ir + i__], abs(d__1));
	    if (tm <= toler * sqrt(wmax + 100.)) {
		dwork[ir + nei] = dwork[ir + i__];
		dwork[ii + nei] = dwork[ii + i__];
		++nei;
	    }
/* L290: */
	}

    }

    if (nei == 0) {

/*           There is no eigenvalue on the boundary of the stability */
/*           domain for G = ( ONE + TOL )*GAMMAL. The norm was found. */

	gpeak[1] = gammal;
	gpeak[2] = 1.;
	goto L340;
    }

/*        Compute the frequencies where the gain G is attained and */
/*        generate new test frequencies. */

    nws = 0;

    if (discr) {

	i__1 = nei - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    tm = atan2(dwork[ii + i__], dwork[ir + i__]);
	    dwork[ir + i__] = max(eps,tm);
	    ++nws;
/* L300: */
	}

    } else {

	j = 0;

	i__1 = nei - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    if (dwork[ii + i__] > eps) {
		dwork[ir + nws] = dwork[ii + i__];
		++nws;
	    } else if (dwork[ii + i__] == eps) {
		++j;
		if (j == 1) {
		    dwork[ir + nws] = eps;
		    ++nws;
		}
	    }
/* L310: */
	}

    }

    dlasrt_("Increasing", &nws, &dwork[ir], &ierr, (ftnlen)10);
    lw = 1;

    i__1 = nws - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (dwork[ir + lw - 1] != dwork[ir + i__]) {
	    dwork[ir + lw] = dwork[ir + i__];
	    ++lw;
	}
/* L320: */
    }

    if (lw == 1) {
	if (iter == 1 && nws >= 1) {

/*              Duplicate the frequency trying to force iteration. */

	    dwork[ir + 1] = dwork[ir];
	    ++lw;
	} else {

/*              The norm was found. */

	    gpeak[1] = gammal;
	    gpeak[2] = 1.;
	    goto L340;
	}
    }

/*        Form the vector of mid-points and compute the gain at new test */
/*        frequencies. Save the current lower bound. */

    iwrk = ir + lw;
    gammas = gammal;

    i__1 = lw - 2;
    for (i__ = 0; i__ <= i__1; ++i__) {
	if (discr) {
	    omega = (dwork[ir + i__] + dwork[ir + i__ + 1]) / 2.;
	} else {
	    omega = sqrt(dwork[ir + i__] * dwork[ir + i__ + 1]);
	}

/*           Additional workspace:  need   LDW2, see above; */
/*                                  prefer larger. */

	i__2 = *ldwork - iwrk + 1;
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &omega, &dwork[ia], n, &
		dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &
		iwork[1], &dwork[iwrk], &i__2, &cwork[1], lcwork, &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
	i__2 = (integer) cwork[1].r;
	maxcwk = max(i__2,maxcwk);
/* Computing MAX */
	i__2 = (integer) dwork[iwrk] + iwrk - 1;
	maxwrk = max(i__2,maxwrk);
	if (discr) {
	    tm = (d__1 = atan2(sin(omega), cos(omega)), abs(d__1));
	} else {
	    tm = omega;
	}
	if (ierr >= 1 && ierr <= *n) {
	    gpeak[1] = 1.;
	    fpeak[1] = tm;
	    gpeak[2] = 0.;
	    fpeak[2] = 1.;
	    goto L340;
	} else if (ierr == *n + 1) {
	    *info = 3;
	    return 0;
	}

	if (gammal < gamma) {
	    gammal = gamma;
	    fpeak[1] = tm;
	    fpeak[2] = 1.;
	}
/* L330: */
    }

/*        If the lower bound has not been improved, return. (This is a */
/*        safeguard against undetected modes of Hamiltonian matrix on the */
/*        boundary of the stability domain.) */

    if (gammal < gammas * (*tol / 10. + 1.)) {
	gpeak[1] = gammal;
	gpeak[2] = 1.;
	goto L340;
    }

/*     END WHILE */

    if (iter <= 30) {
	goto L120;
    } else {
	*info = 4;
	return 0;
    }

L340:
    dwork[1] = (doublereal) maxwrk;
    cwork[1].r = (doublereal) maxcwk, cwork[1].i = 0.;
    return 0;
/* *** Last line of AB13DD *** */
} /* ab13dd_ */

