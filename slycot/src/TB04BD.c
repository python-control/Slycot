/* TB04BD.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;

/* Subroutine */ int tb04bd_(char *jobd, char *order, char *equil, integer *n,
	 integer *m, integer *p, integer *md, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *ign, integer *ldign, integer *
	igd, integer *ldigd, doublereal *gn, doublereal *gd, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen jobd_len, ftnlen order_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal x, z__[1];
    static integer ia, ib, ic, jj, im, ip, iz, iac, icc;
    static doublereal dij;
    static integer ias, iip, irp, ipm1, ierr, itau;
    static doublereal epsn;
    static integer itau1;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mc01pd_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), tb01id_(char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb04bx_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    ), tb01zd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), mc01py_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *);
    static doublereal anorm;
    static logical dijnz, withd;
    static integer ncont;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork, jwork1;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical fndeig, ascend;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    extern /* Subroutine */ int dhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer wrkopt;


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

/*     To compute the transfer function matrix G of a state-space */
/*     representation (A,B,C,D) of a linear time-invariant multivariable */
/*     system, using the pole-zeros method. Each element of the transfer */
/*     function matrix is returned in a cancelled, minimal form, with */
/*     numerator and denominator polynomials stored either in increasing */
/*     or decreasing order of the powers of the indeterminate. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state-space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed to be a zero matrix. */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system (A,B,C,D).  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1. An */
/*             upper bound for MD is N+1.  MD >= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, if EQUIL = 'S', the leading N-by-N part of this */
/*             array contains the balanced matrix inv(S)*A*S, as returned */
/*             by SLICOT Library routine TB01ID. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the contents of B are destroyed: all elements but */
/*             those in the first row are set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, if EQUIL = 'S', the leading P-by-N part of this */
/*             array contains the balanced matrix C*S, as returned by */
/*             SLICOT Library routine TB01ID. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If JOBD = 'D', the leading P-by-M part of this array must */
/*             contain the matrix D. */
/*             If JOBD = 'Z', the array D is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     IGN     (output) INTEGER array, dimension (LDIGN,M) */
/*             The leading P-by-M part of this array contains the degrees */
/*             of the numerator polynomials in the transfer function */
/*             matrix G. Specifically, the (i,j) element of IGN contains */
/*             the degree of the numerator polynomial of the transfer */
/*             function G(i,j) from the j-th input to the i-th output. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (output) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array contains the degrees */
/*             of the denominator polynomials in the transfer function */
/*             matrix G. Specifically, the (i,j) element of IGD contains */
/*             the degree of the denominator polynomial of the transfer */
/*             function G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array contains the coefficients of the numerator */
/*             polynomials, Num(i,j), of the transfer function matrix G. */
/*             The polynomials are stored in a column-wise order, i.e., */
/*             Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), Num(2,2), */
/*             ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., Num(P,M); */
/*             MD memory locations are reserved for each polynomial, */
/*             hence, the (i,j) polynomial is stored starting from the */
/*             location ((j-1)*P+i-1)*MD+1. The coefficients appear in */
/*             increasing or decreasing order of the powers of the */
/*             indeterminate, according to ORDER. */

/*     GD      (output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array contains the coefficients of the denominator */
/*             polynomials, Den(i,j), of the transfer function matrix G. */
/*             The polynomials are stored in the same way as the */
/*             numerator polynomials. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the */
/*             controllability of a single-input system (A,b) or (A',c'), */
/*             where b and c' are columns in B and C' (C transposed). If */
/*             the user sets TOL > 0, then the given value of TOL is used */
/*             as an absolute tolerance; elements with absolute value */
/*             less than TOL are considered neglijible. If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by TOLDEF = N*EPS*MAX( NORM(A), NORM(bc) ) is used */
/*             instead, where EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH), and bc denotes the currently used */
/*             column in B or C' (see METHOD). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*(N+P) + */
/*                              MAX( N + MAX( N,P ), N*(2*N+5))) */
/*             If N >= P, N >= 1, the formula above can be written as */
/*             LDWORK >= N*(3*N + P + 5). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the QR algorithm failed to converge when trying to */
/*                   compute the zeros of a transfer function; */
/*             = 2:  the QR algorithm failed to converge when trying to */
/*                   compute the poles of a transfer function. */
/*                   The errors INFO = 1 or 2 are unlikely to appear. */

/*     METHOD */

/*     The routine implements the pole-zero method proposed in [1]. */
/*     This method is based on an algorithm for computing the transfer */
/*     function of a single-input single-output (SISO) system. */
/*     Let (A,b,c,d) be a SISO system. Its transfer function is computed */
/*     as follows: */

/*     1) Find a controllable realization (Ac,bc,cc) of (A,b,c). */
/*     2) Find an observable realization (Ao,bo,co) of (Ac,bc,cc). */
/*     3) Compute the r eigenvalues of Ao (the poles of (Ao,bo,co)). */
/*     4) Compute the zeros of (Ao,bo,co,d). */
/*     5) Compute the gain of (Ao,bo,co,d). */

/*     This algorithm can be implemented using only orthogonal */
/*     transformations [1]. However, for better efficiency, the */
/*     implementation in TB04BD uses one elementary transformation */
/*     in Step 4 and r elementary transformations in Step 5 (to reduce */
/*     an upper Hessenberg matrix to upper triangular form). These */
/*     special elementary transformations are numerically stable */
/*     in practice. */

/*     In the multi-input multi-output (MIMO) case, the algorithm */
/*     computes each element (i,j) of the transfer function matrix G, */
/*     for i = 1 : P, and for j = 1 : M. For efficiency reasons, Step 1 */
/*     is performed once for each value of j (each column of B). The */
/*     matrices Ac and Ao result in Hessenberg form. */

/*     REFERENCES */

/*     [1] Varga, A. and Sima, V. */
/*         Numerically Stable Algorithm for Transfer Function Matrix */
/*         Evaluation. */
/*         Int. J. Control, vol. 33, nr. 6, pp. 1123-1133, 1981. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable in practice and requires about */
/*     20*N**3 floating point operations at most, but usually much less. */

/*     FURTHER COMMENTS */

/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Partly based on the BIMASC Library routine TSMT1 by A. Varga. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue, state-space representation, transfer function, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

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
    ign_dim1 = *ldign;
    ign_offset = 1 + ign_dim1;
    ign -= ign_offset;
    igd_dim1 = *ldigd;
    igd_offset = 1 + igd_dim1;
    igd -= igd_offset;
    --gn;
    --gd;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
    ascend = lsame_(order, "I", (ftnlen)1, (ftnlen)1);
    if (! withd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ascend && ! lsame_(order, "D", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*md < 1) {
	*info = -7;
    } else if (*lda < max(1,*n)) {
	*info = -9;
    } else if (*ldb < max(1,*n)) {
	*info = -11;
    } else if (*ldc < max(1,*p)) {
	*info = -13;
    } else if (*ldd < 1 || withd && *ldd < *p) {
	*info = -15;
    } else if (*ldign < max(1,*p)) {
	*info = -17;
    } else if (*ldigd < max(1,*p)) {
	*info = -19;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n + max(*n,*p), i__4 = *n * ((*n << 1) + 5);
	i__1 = 1, i__2 = *n * (*n + *p) + max(i__3,i__4);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -25;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB04BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize GN and GD to zero. */

    z__[0] = 0.;
    i__1 = *p * *m * *md;
    dcopy_(&i__1, z__, &c__0, &gn[1], &c__1);
    i__1 = *p * *m * *md;
    dcopy_(&i__1, z__, &c__0, &gd[1], &c__1);

/*     Quick return if possible. */

/* Computing MIN */
    i__1 = min(*n,*p);
    if (min(i__1,*m) == 0) {
	if (min(*p,*m) > 0) {
	    k = 1;

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ign[i__ + j * ign_dim1] = 0;
		    igd[i__ + j * igd_dim1] = 0;
		    if (withd) {
			gn[k] = d__[i__ + j * d_dim1];
		    }
		    gd[k] = 1.;
		    k += *md;
/* L10: */
		}

/* L20: */
	    }

	}
	dwork[1] = 1.;
	return 0;
    }

/*     Prepare the computation of the default tolerance. */

    toldef = *tol;
    if (toldef <= 0.) {
	epsn = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	anorm = dlange_("Frobenius", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)9);
    }

/*     Initializations. */

    ia = 1;
    ic = ia + *n * *n;
    itau = ic + *p * *n;
    jwork = itau + *n;
    iac = itau;

    k = 1;
    dij = 0.;

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(S)*A*S,  B <- inv(S)*B and C <- C*S, where S is a */
/*        diagonal scaling matrix. */
/*        Workspace: need   N. */

	maxred = 100.;
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], &ierr, (ftnlen)3);
    }

/*     Compute the transfer function matrix of the system (A,B,C,D). */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {

/*        Save A and C. */
/*        Workspace: need   W1 = N*(N+P). */

	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ic], p, (ftnlen)4);

/*        Remove the uncontrollable part of the system (A,B(J),C). */
/*        Workspace: need   W1+N+MAX(N,P); */
/*                   prefer larger. */

	i__2 = *ldwork - jwork + 1;
	tb01zd_("No Z", n, p, &dwork[ia], n, &b[j * b_dim1 + 1], &dwork[ic], 
		p, &ncont, z__, &c__1, &dwork[itau], tol, &dwork[jwork], &
		i__2, &ierr, (ftnlen)4);
	if (j == 1) {
	    wrkopt = (integer) dwork[jwork] + jwork - 1;
	}

	ib = iac + ncont * ncont;
	icc = ib + ncont;
	itau1 = icc + ncont;
	irp = itau1;
	iip = irp + ncont;
	ias = iip + ncont;
	jwork1 = ias + ncont * ncont;

	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (withd) {
		dij = d__[i__ + j * d_dim1];
	    }
	    if (ncont > 0) {

/*              Form the matrices of the state-space representation of */
/*              the dual system for the controllable part. */
/*              Workspace: need   W2 = W1+N*(N+2). */

		ma02ad_("Full", &ncont, &ncont, &dwork[ia], n, &dwork[iac], &
			ncont, (ftnlen)4);
		dcopy_(&ncont, &b[j * b_dim1 + 1], &c__1, &dwork[ib], &c__1);
		dcopy_(&ncont, &dwork[ic + i__ - 1], p, &dwork[icc], &c__1);

/*              Remove the unobservable part of the system (A,B(J),C(I)). */
/*              Workspace: need   W2+2*N; */
/*                         prefer larger. */

		i__3 = *ldwork - iip + 1;
		tb01zd_("No Z", &ncont, &c__1, &dwork[iac], &ncont, &dwork[
			icc], &dwork[ib], &c__1, &ip, z__, &c__1, &dwork[
			itau1], tol, &dwork[iip], &i__3, &ierr, (ftnlen)4);
		if (i__ == 1) {
/* Computing MAX */
		    i__3 = wrkopt, i__4 = (integer) dwork[iip] + iip - 1;
		    wrkopt = max(i__3,i__4);
		}

		if (ip > 0) {

/*                 Save the state matrix of the minimal part. */
/*                 Workspace: need   W3 = W2+N*(N+2). */

		    dlacpy_("Full", &ip, &ip, &dwork[iac], &ncont, &dwork[ias]
			    , &ip, (ftnlen)4);

/*                 Compute the poles of the transfer function. */
/*                 Workspace: need   W3+N; */
/*                            prefer larger. */

		    i__3 = *ldwork - jwork1 + 1;
		    dhseqr_("Eigenvalues", "No vectors", &ip, &c__1, &ip, &
			    dwork[iac], &ncont, &dwork[irp], &dwork[iip], z__,
			     &c__1, &dwork[jwork1], &i__3, &ierr, (ftnlen)11, 
			    (ftnlen)10);
		    if (ierr != 0) {
			*info = 2;
			return 0;
		    }
/* Computing MAX */
		    i__3 = wrkopt, i__4 = (integer) dwork[jwork1] + jwork1 - 
			    1;
		    wrkopt = max(i__3,i__4);

/*                 Compute the zeros of the transfer function. */

		    ipm1 = ip - 1;
		    dijnz = withd && dij != 0.;
		    fndeig = dijnz || ipm1 > 0;
		    if (! fndeig) {
			iz = 0;
		    } else if (dijnz) {

/*                    Add the contribution due to D(i,j). */
/*                    Note that the matrix whose eigenvalues have to */
/*                    be computed remains in an upper Hessenberg form. */

			iz = ip;
			dlacpy_("Full", &iz, &iz, &dwork[ias], &ip, &dwork[
				iac], &ncont, (ftnlen)4);
			d__1 = -dwork[icc] / dij;
			daxpy_(&iz, &d__1, &dwork[ib], &c__1, &dwork[iac], &
				ncont);
		    } else {
			if (*tol <= 0.) {
/* Computing MAX */
			    d__1 = anorm, d__2 = dlange_("Frobenius", &ip, &
				    c__1, &dwork[ib], &c__1, &dwork[1], (
				    ftnlen)9);
			    toldef = epsn * max(d__1,d__2);
			}

			i__3 = ipm1;
			for (im = 1; im <= i__3; ++im) {
			    if ((d__1 = dwork[ib + im - 1], abs(d__1)) > 
				    toldef) {
				goto L40;
			    }
/* L30: */
			}

			iz = 0;
			goto L50;

L40:

/*                    Restore (part of) the saved state matrix. */

			iz = ip - im;
			dlacpy_("Full", &iz, &iz, &dwork[ias + im * (ip + 1)],
				 &ip, &dwork[iac], &ncont, (ftnlen)4);

/*                    Apply the output injection. */

			d__1 = -dwork[ias + im * (ip + 1) - ip] / dwork[ib + 
				im - 1];
			daxpy_(&iz, &d__1, &dwork[ib + im], &c__1, &dwork[iac]
				, &ncont);
		    }

		    if (fndeig) {

/*                    Find the zeros. */
/*                    Workspace: need   W3+N; */
/*                               prefer larger. */

			i__3 = *ldwork - jwork1 + 1;
			dhseqr_("Eigenvalues", "No vectors", &iz, &c__1, &iz, 
				&dwork[iac], &ncont, &gn[k], &gd[k], z__, &
				c__1, &dwork[jwork1], &i__3, &ierr, (ftnlen)
				11, (ftnlen)10);
			if (ierr != 0) {
			    *info = 1;
			    return 0;
			}
		    }

/*                 Compute the gain. */

L50:
		    if (dijnz) {
			x = dij;
		    } else {
			tb04bx_(&ip, &iz, &dwork[ias], &ip, &dwork[icc], &
				dwork[ib], &dij, &dwork[irp], &dwork[iip], &
				gn[k], &gd[k], &x, &iwork[1]);
		    }

/*                 Form the numerator coefficients in increasing or */
/*                 decreasing powers of the indeterminate. */
/*                 IAS is used here as pointer to the workspace. */

		    if (ascend) {
			mc01pd_(&iz, &gn[k], &gd[k], &dwork[ib], &dwork[ias], 
				&ierr);
		    } else {
			mc01py_(&iz, &gn[k], &gd[k], &dwork[ib], &dwork[ias], 
				&ierr);
		    }
		    jj = k;

		    i__3 = ib + iz;
		    for (l = ib; l <= i__3; ++l) {
			gn[jj] = dwork[l] * x;
			++jj;
/* L60: */
		    }

/*                 Form the denominator coefficients. */

		    if (ascend) {
			mc01pd_(&ip, &dwork[irp], &dwork[iip], &gd[k], &dwork[
				ias], &ierr);
		    } else {
			mc01py_(&ip, &dwork[irp], &dwork[iip], &gd[k], &dwork[
				ias], &ierr);
		    }
		    ign[i__ + j * ign_dim1] = iz;
		    igd[i__ + j * igd_dim1] = ip;
		} else {

/*                 Null element. */

		    ign[i__ + j * ign_dim1] = 0;
		    igd[i__ + j * igd_dim1] = 0;
		    gn[k] = dij;
		    gd[k] = 1.;
		}

	    } else {

/*              Null element. */

		ign[i__ + j * ign_dim1] = 0;
		igd[i__ + j * igd_dim1] = 0;
		gn[k] = dij;
		gd[k] = 1.;
	    }

	    k += *md;
/* L70: */
	}

/* L80: */
    }

    return 0;
/* *** Last line of TB04BD *** */
} /* tb04bd_ */

