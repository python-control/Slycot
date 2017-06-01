/* TB04CD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tb04cd_(char *jobd, char *equil, integer *n, integer *m, 
	integer *p, integer *npz, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, integer *nz, integer *ldnz, integer *np, integer *ldnp, 
	doublereal *zerosr, doublereal *zerosi, doublereal *polesr, 
	doublereal *polesi, doublereal *gains, integer *ldgain, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen jobd_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, gains_dim1, gains_offset, np_dim1, np_offset, nz_dim1, 
	    nz_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal z__[1];
    static integer ia, ib, ic, im, ip, iz, iac, icc;
    static doublereal dij;
    static integer ias, jwk, ipm1, ierr, itau;
    static doublereal epsn;
    static integer itau1;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    tb01id_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb04bx_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    ), tb01zd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal anorm;
    static logical dijnz, withd;
    static integer ncont;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork, jwork1;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical fndeig;
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
/*     system, using the pole-zeros method. The transfer function matrix */
/*     is returned in a minimal pole-zero-gain form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state-space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed to be a zero matrix. */

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

/*     NPZ     (input) INTEGER */
/*             The maximum number of poles or zeros of the single-input */
/*             single-output channels in the system. An upper bound */
/*             for NPZ is N.  NPZ >= 0. */

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

/*     NZ      (output) INTEGER array, dimension (LDNZ,M) */
/*             The leading P-by-M part of this array contains the numbers */
/*             of zeros of the elements of the transfer function */
/*             matrix G. Specifically, the (i,j) element of NZ contains */
/*             the number of zeros of the transfer function G(i,j) from */
/*             the j-th input to the i-th output. */

/*     LDNZ    INTEGER */
/*             The leading dimension of array NZ.  LDNZ >= max(1,P). */

/*     NP      (output) INTEGER array, dimension (LDNP,M) */
/*             The leading P-by-M part of this array contains the numbers */
/*             of poles of the elements of the transfer function */
/*             matrix G. Specifically, the (i,j) element of NP contains */
/*             the number of poles of the transfer function G(i,j). */

/*     LDNP    INTEGER */
/*             The leading dimension of array NP.  LDNP >= max(1,P). */

/*     ZEROSR  (output) DOUBLE PRECISION array, dimension (P*M*NPZ) */
/*             This array contains the real parts of the zeros of the */
/*             transfer function matrix G. The real parts of the zeros */
/*             are stored in a column-wise order, i.e., for the transfer */
/*             functions (1,1), (2,1), ..., (P,1), (1,2), (2,2), ..., */
/*             (P,2), ..., (1,M), (2,M), ..., (P,M); NPZ memory locations */
/*             are reserved for each transfer function, hence, the real */
/*             parts of the zeros for the (i,j) transfer function */
/*             are stored starting from the location ((j-1)*P+i-1)*NPZ+1. */
/*             Pairs of complex conjugate zeros are stored in consecutive */
/*             memory locations. Note that only the first NZ(i,j) entries */
/*             are initialized for the (i,j) transfer function. */

/*     ZEROSI  (output) DOUBLE PRECISION array, dimension (P*M*NPZ) */
/*             This array contains the imaginary parts of the zeros of */
/*             the transfer function matrix G, stored in a similar way */
/*             as the real parts of the zeros. */

/*     POLESR  (output) DOUBLE PRECISION array, dimension (P*M*NPZ) */
/*             This array contains the real parts of the poles of the */
/*             transfer function matrix G, stored in the same way as */
/*             the zeros. Note that only the first NP(i,j) entries are */
/*             initialized for the (i,j) transfer function. */

/*     POLESI  (output) DOUBLE PRECISION array, dimension (P*M*NPZ) */
/*             This array contains the imaginary parts of the poles of */
/*             the transfer function matrix G, stored in the same way as */
/*             the poles. */

/*     GAINS   (output) DOUBLE PRECISION array, dimension (LDGAIN,M) */
/*             The leading P-by-M part of this array contains the gains */
/*             of the transfer function matrix G. Specifically, */
/*             GAINS(i,j) contains the gain of the transfer function */
/*             G(i,j). */

/*     LDGAIN  INTEGER */
/*             The leading dimension of array GAINS.  LDGAIN >= max(1,P). */

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
/*                              MAX( N + MAX( N,P ), N*(2*N+3))) */
/*             If N >= P, N >= 1, the formula above can be written as */
/*             LDWORK >= N*(3*N + P + 3). */
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
/*     implementation in TB04CD uses one elementary transformation */
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

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */

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
    nz_dim1 = *ldnz;
    nz_offset = 1 + nz_dim1;
    nz -= nz_offset;
    np_dim1 = *ldnp;
    np_offset = 1 + np_dim1;
    np -= np_offset;
    --zerosr;
    --zerosi;
    --polesr;
    --polesi;
    gains_dim1 = *ldgain;
    gains_offset = 1 + gains_dim1;
    gains -= gains_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
    if (! withd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*npz < 0) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldc < max(1,*p)) {
	*info = -12;
    } else if (*ldd < 1 || withd && *ldd < *p) {
	*info = -14;
    } else if (*ldnz < max(1,*p)) {
	*info = -16;
    } else if (*ldnp < max(1,*p)) {
	*info = -18;
    } else if (*ldgain < max(1,*p)) {
	*info = -24;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n + max(*n,*p), i__4 = *n * ((*n << 1) + 3);
	i__1 = 1, i__2 = *n * (*n + *p) + max(i__3,i__4);
	if (*ldwork < max(i__1,i__2)) {
	    *info = -28;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TB04CD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    dij = 0.;
/* Computing MIN */
    i__1 = min(*n,*p);
    if (min(i__1,*m) == 0) {
	if (min(*p,*m) > 0) {

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {

		i__2 = *p;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    nz[i__ + j * nz_dim1] = 0;
		    np[i__ + j * np_dim1] = 0;
		    if (withd) {
			dij = d__[i__ + j * d_dim1];
		    }
		    gains[i__ + j * gains_dim1] = dij;
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

/*     Compute the transfer function matrix of the system (A,B,C,D), */
/*     in the pole-zero-gain form. */

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
	jwk = itau1 + ncont;
	ias = itau1;
	jwork1 = ias + ncont * ncont;

	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (ncont > 0) {
		if (withd) {
		    dij = d__[i__ + j * d_dim1];
		}

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

		i__3 = *ldwork - jwk + 1;
		tb01zd_("No Z", &ncont, &c__1, &dwork[iac], &ncont, &dwork[
			icc], &dwork[ib], &c__1, &ip, z__, &c__1, &dwork[
			itau1], tol, &dwork[jwk], &i__3, &ierr, (ftnlen)4);
		if (i__ == 1) {
/* Computing MAX */
		    i__3 = wrkopt, i__4 = (integer) dwork[jwk] + jwk - 1;
		    wrkopt = max(i__3,i__4);
		}

		if (ip > 0) {

/*                 Save the state matrix of the minimal part. */
/*                 Workspace: need   W3 = W2+N*N. */

		    dlacpy_("Full", &ip, &ip, &dwork[iac], &ncont, &dwork[ias]
			    , &ip, (ftnlen)4);

/*                 Compute the poles of the transfer function. */
/*                 Workspace: need   W3+N; */
/*                            prefer larger. */

		    i__3 = *ldwork - jwork1 + 1;
		    dhseqr_("Eigenvalues", "No vectors", &ip, &c__1, &ip, &
			    dwork[iac], &ncont, &polesr[k], &polesi[k], z__, &
			    c__1, &dwork[jwork1], &i__3, &ierr, (ftnlen)11, (
			    ftnlen)10);
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
				&dwork[iac], &ncont, &zerosr[k], &zerosi[k], 
				z__, &c__1, &dwork[jwork1], &i__3, &ierr, (
				ftnlen)11, (ftnlen)10);
			if (ierr != 0) {
			    *info = 1;
			    return 0;
			}
		    }

/*                 Compute the gain. */

L50:
		    if (dijnz) {
			gains[i__ + j * gains_dim1] = dij;
		    } else {
			tb04bx_(&ip, &iz, &dwork[ias], &ip, &dwork[icc], &
				dwork[ib], &dij, &polesr[k], &polesi[k], &
				zerosr[k], &zerosi[k], &gains[i__ + j * 
				gains_dim1], &iwork[1]);
		    }
		    nz[i__ + j * nz_dim1] = iz;
		    np[i__ + j * np_dim1] = ip;
		} else {

/*                 Null element. */

		    nz[i__ + j * nz_dim1] = 0;
		    np[i__ + j * np_dim1] = 0;
		}

	    } else {

/*              Null element. */

		nz[i__ + j * nz_dim1] = 0;
		np[i__ + j * np_dim1] = 0;
	    }

	    k += *npz;
/* L70: */
	}

/* L80: */
    }

    return 0;
/* *** Last line of TB04CD *** */
} /* tb04cd_ */

