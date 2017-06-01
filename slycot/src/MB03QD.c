/* MB03QD.f -- translated by f2c (version 20100827).
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

static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;

/* Subroutine */ int mb03qd_(char *dico, char *stdom, char *jobu, integer *n, 
	integer *nlow, integer *nsup, doublereal *alpha, doublereal *a, 
	integer *lda, doublereal *u, integer *ldu, integer *ndim, doublereal *
	dwork, integer *info, ftnlen dico_len, ftnlen stdom_len, ftnlen 
	jobu_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer l;
    static doublereal e1, e2;
    static integer ib, lm1, nup;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int mb03qy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal tlambd;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrexc_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, ftnlen);
    static logical lstdom;


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

/*     To reorder the diagonal blocks of a principal submatrix of an */
/*     upper quasi-triangular matrix A together with their eigenvalues by */
/*     constructing an orthogonal similarity transformation UT. */
/*     After reordering, the leading block of the selected submatrix of A */
/*     has eigenvalues in a suitably defined domain of interest, usually */
/*     related to stability/instability in a continuous- or discrete-time */
/*     sense. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the spectrum separation to be */
/*             performed as follows: */
/*             = 'C':  continuous-time sense; */
/*             = 'D':  discrete-time sense. */

/*     STDOM   CHARACTER*1 */
/*             Specifies whether the domain of interest is of stability */
/*             type (left part of complex plane or inside of a circle) */
/*             or of instability type (right part of complex plane or */
/*             outside of a circle) as follows: */
/*             = 'S':  stability type domain; */
/*             = 'U':  instability type domain. */

/*     JOBU    CHARACTER*1 */
/*             Indicates how the performed orthogonal transformations UT */
/*             are accumulated, as follows: */
/*             = 'I':  U is initialized to the unit matrix and the matrix */
/*                     UT is returned in U; */
/*             = 'U':  the given matrix U is updated and the matrix U*UT */
/*                     is returned in U. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and U.  N >= 1. */

/*     NLOW,   (input) INTEGER */
/*     NSUP    NLOW and NSUP specify the boundary indices for the rows */
/*             and columns of the principal submatrix of A whose diagonal */
/*             blocks are to be reordered.  1 <= NLOW <= NSUP <= N. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The boundary of the domain of interest for the eigenvalues */
/*             of A. If DICO = 'C', ALPHA is the boundary value for the */
/*             real parts of eigenvalues, while for DICO = 'D', */
/*             ALPHA >= 0 represents the boundary value for the moduli of */
/*             eigenvalues. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain a matrix in a real Schur form whose 1-by-1 and */
/*             2-by-2 diagonal blocks between positions NLOW and NSUP */
/*             are to be reordered. */
/*             On exit, the leading N-by-N part contains the ordered */
/*             real Schur matrix UT' * A * UT with the elements below the */
/*             first subdiagonal set to zero. */
/*             The leading NDIM-by-NDIM part of the principal submatrix */
/*             D = A(NLOW:NSUP,NLOW:NSUP) has eigenvalues in the domain */
/*             of interest and the trailing part of this submatrix has */
/*             eigenvalues outside the domain of interest. */
/*             The domain of interest for lambda(D), the eigenvalues of */
/*             D, is defined by the parameters ALPHA, DICO and STDOM as */
/*             follows: */
/*               For DICO = 'C': */
/*                  Real(lambda(D)) < ALPHA if STDOM = 'S'; */
/*                  Real(lambda(D)) > ALPHA if STDOM = 'U'. */
/*               For DICO = 'D': */
/*                  Abs(lambda(D)) < ALPHA if STDOM = 'S'; */
/*                  Abs(lambda(D)) > ALPHA if STDOM = 'U'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             On entry with JOBU = 'U', the leading N-by-N part of this */
/*             array must contain a transformation matrix (e.g. from a */
/*             previous call to this routine). */
/*             On exit, if JOBU = 'U', the leading N-by-N part of this */
/*             array contains the product of the input matrix U and the */
/*             orthogonal matrix UT used to reorder the diagonal blocks */
/*             of A. */
/*             On exit, if JOBU = 'I', the leading N-by-N part of this */
/*             array contains the matrix UT of the performed orthogonal */
/*             transformations. */
/*             Array U need not be set on entry if JOBU = 'I'. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= N. */

/*     NDIM    (output) INTEGER */
/*             The number of eigenvalues of the selected principal */
/*             submatrix lying inside the domain of interest. */
/*             If NLOW = 1, NDIM is also the dimension of the invariant */
/*             subspace corresponding to the eigenvalues of the leading */
/*             NDIM-by-NDIM submatrix. In this case, if U is the */
/*             orthogonal transformation matrix used to compute and */
/*             reorder the real Schur form of A, its first NDIM columns */
/*             form an orthonormal basis for the above invariant */
/*             subspace. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  A(NLOW,NLOW-1) is nonzero, i.e. A(NLOW,NLOW) is not */
/*                   the leading element of a 1-by-1 or 2-by-2 diagonal */
/*                   block of A, or A(NSUP+1,NSUP) is nonzero, i.e. */
/*                   A(NSUP,NSUP) is not the bottom element of a 1-by-1 */
/*                   or 2-by-2 diagonal block of A; */
/*             = 2:  two adjacent blocks are too close to swap (the */
/*                   problem is very ill-conditioned). */

/*     METHOD */

/*     Given an upper quasi-triangular matrix A with 1-by-1 or 2-by-2 */
/*     diagonal blocks, the routine reorders its diagonal blocks along */
/*     with its eigenvalues by performing an orthogonal similarity */
/*     transformation UT' * A * UT. The column transformation UT is also */
/*     performed on the given (initial) transformation U (resulted from */
/*     a possible previous step or initialized as the identity matrix). */
/*     After reordering, the eigenvalues inside the region specified by */
/*     the parameters ALPHA, DICO and STDOM appear at the top of */
/*     the selected diagonal block between positions NLOW and NSUP. */
/*     In other words, lambda(A(NLOW:NSUP,NLOW:NSUP)) are ordered such */
/*     that lambda(A(NLOW:NLOW+NDIM-1,NLOW:NLOW+NDIM-1)) are inside and */
/*     lambda(A(NLOW+NDIM:NSUP,NLOW+NDIM:NSUP)) are outside the domain */
/*     of interest. If NLOW = 1, the first NDIM columns of U*UT span the */
/*     corresponding invariant subspace of A. */

/*     REFERENCES */

/*     [1] Stewart, G.W. */
/*         HQR3 and EXCHQZ: FORTRAN subroutines for calculating and */
/*         ordering the eigenvalues of a real upper Hessenberg matrix. */
/*         ACM TOMS, 2, pp. 275-280, 1976. */

/*     NUMERICAL ASPECTS */
/*                                         3 */
/*     The algorithm requires less than 4*N  operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     April 1998. Based on the RASP routine SEOR1. */

/*     KEYWORDS */

/*     Eigenvalues, invariant subspace, orthogonal transformation, real */
/*     Schur form, similarity transformation. */

/*    ****************************************************************** */

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
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
    lstdom = lsame_(stdom, "S", (ftnlen)1, (ftnlen)1);

/*     Check input scalar arguments. */

    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
	*info = -1;
    } else if (! (lstdom || lsame_(stdom, "U", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (! (lsame_(jobu, "I", (ftnlen)1, (ftnlen)1) || lsame_(jobu, 
	    "U", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*n < 1) {
	*info = -4;
    } else if (*nlow < 1) {
	*info = -5;
    } else if (*nlow > *nsup || *nsup > *n) {
	*info = -6;
    } else if (discr && *alpha < 0.) {
	*info = -7;
    } else if (*lda < *n) {
	*info = -9;
    } else if (*ldu < *n) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB03QD", &i__1, (ftnlen)6);
	return 0;
    }

    if (*nlow > 1) {
	if (a[*nlow + (*nlow - 1) * a_dim1] != 0.) {
	    *info = 1;
	}
    }
    if (*nsup < *n) {
	if (a[*nsup + 1 + *nsup * a_dim1] != 0.) {
	    *info = 1;
	}
    }
    if (*info != 0) {
	return 0;
    }

/*     Initialize U with an identity matrix if necessary. */

    if (lsame_(jobu, "I", (ftnlen)1, (ftnlen)1)) {
	dlaset_("Full", n, n, &c_b11, &c_b12, &u[u_offset], ldu, (ftnlen)4);
    }

    *ndim = 0;
    l = *nsup;
    nup = *nsup;

/*     NUP is the minimal value such that the submatrix A(i,j) with */
/*     NUP+1 <= i,j <= NSUP contains no eigenvalues inside the domain of */
/*     interest. L is such that all the eigenvalues of the submatrix */
/*     A(i,j) with L+1 <= i,j <= NUP lie inside the domain of interest. */

/*     WHILE( L >= NLOW ) DO */

L10:
    if (l >= *nlow) {
	ib = 1;
	if (l > *nlow) {
	    lm1 = l - 1;
	    if (a[l + lm1 * a_dim1] != 0.) {
		mb03qy_(n, &lm1, &a[a_offset], lda, &u[u_offset], ldu, &e1, &
			e2, info);
		if (a[l + lm1 * a_dim1] != 0.) {
		    ib = 2;
		}
	    }
	}
	if (discr) {
	    if (ib == 1) {
		tlambd = (d__1 = a[l + l * a_dim1], abs(d__1));
	    } else {
		tlambd = dlapy2_(&e1, &e2);
	    }
	} else {
	    if (ib == 1) {
		tlambd = a[l + l * a_dim1];
	    } else {
		tlambd = e1;
	    }
	}
	if (lstdom && tlambd < *alpha || ! lstdom && tlambd > *alpha) {
	    *ndim += ib;
	    l -= ib;
	} else {
	    if (*ndim != 0) {
		dtrexc_("V", n, &a[a_offset], lda, &u[u_offset], ldu, &l, &
			nup, &dwork[1], info, (ftnlen)1);
		if (*info != 0) {
		    *info = 2;
		    return 0;
		}
		--nup;
		--l;
	    } else {
		nup -= ib;
		l -= ib;
	    }
	}
	goto L10;
    }

/*     END WHILE 10 */

    return 0;
/* *** Last line of MB03QD *** */
} /* mb03qd_ */

