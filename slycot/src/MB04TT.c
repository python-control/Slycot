/* MB04TT.f -- translated by f2c (version 20100827).
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
static logical c_false = FALSE_;

/* Subroutine */ int mb04tt_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *ifira, integer *ifica, integer *nca, doublereal *
	a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer *
	ldq, doublereal *z__, integer *ldz, integer *istair, integer *rank, 
	doublereal *tol, integer *iwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l, ii, kk, mj, ll, ip, nj;
    static doublereal sc, ss;
    static integer jc1, jc2, mk1;
    static doublereal bmx;
    static integer ist1, ist2, lsav;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jpvt;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer itype;
    static logical lzero;
    static integer ifica1, ifira1;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlapmt_(logical *, integer *, integer *, doublereal *, integer *, 
	    integer *);
    static integer mxrank;
    static doublereal eijpvt, bmxnrm;
    static integer istpvt;


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

/*     Let A and E be M-by-N matrices with E in column echelon form. */
/*     Let AA and EE be the following submatrices of A and E: */
/*       AA := A(IFIRA : M ; IFICA : N) */
/*       EE := E(IFIRA : M ; IFICA : N). */
/*     Let Aj and Ej be the following submatrices of AA and EE: */
/*       Aj := A(IFIRA : M ; IFICA : IFICA + NCA - 1) and */
/*       Ej := E(IFIRA : M ; IFICA + NCA : N). */

/*     To transform (AA,EE) such that Aj is row compressed while keeping */
/*     matrix Ej in column echelon form (which may be different from the */
/*     form on entry). */
/*     In fact the routine performs the j-th step of Algorithm 3.2.1 in */
/*     [1]. Furthermore, it determines the rank RANK of the submatrix Ej, */
/*     which is equal to the number of corner points in submatrix Ej. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPDATQ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the orthogonal row transformations, as follows: */
/*             = .FALSE.: Do not form Q; */
/*             = .TRUE.:  The given matrix Q is updated by the orthogonal */
/*                        row transformations used in the reduction. */

/*     UPDATZ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal column transformations, as */
/*             follows: */
/*             = .FALSE.: Do not form Z; */
/*             = .TRUE.:  The given matrix Z is updated by the orthogonal */
/*                        column transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             M is the number of rows of the matrices A, E and Q. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             N is the number of columns of the matrices A, E and Z. */
/*             N >= 0. */

/*     IFIRA   (input) INTEGER */
/*             IFIRA is the first row index of the submatrices Aj and Ej */
/*             in the matrices A and E, respectively. */

/*     IFICA   (input) INTEGER */
/*             IFICA and IFICA + NCA are the first column indices of the */
/*             submatrices Aj and Ej in the matrices A and E, */
/*             respectively. */

/*     NCA     (input) INTEGER */
/*             NCA is the number of columns of the submatrix Aj in A. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, A(IFIRA : M ; IFICA : IFICA + NCA - 1) contains */
/*             the matrix Aj. */
/*             On exit, it contains the matrix A with AA that has been */
/*             row compressed while keeping EE in column echelon form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, E(IFIRA : M ; IFICA + NCA : N) contains the */
/*             matrix Ej which is in column echelon form. */
/*             On exit, it contains the transformed matrix EE which is */
/*             kept in column echelon form. */

/*     LDE     INTEGER */
/*             The leading dimension of array E. LDE >= MAX(1,M). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             On entry, if UPDATQ = .TRUE., then the leading M-by-M */
/*             part of this array must contain a given matrix Q (e.g. */
/*             from a previous call to another SLICOT routine), and on */
/*             exit, the leading M-by-M part of this array contains the */
/*             product of the input matrix Q and the row transformation */
/*             matrix that has transformed the rows of the matrices A */
/*             and E. */
/*             If UPDATQ = .FALSE., the array Q is not referenced and */
/*             can be supplied as a dummy array (i.e. set parameter */
/*             LDQ = 1 and declare this array to be Q(1,1) in the calling */
/*             program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. If UPDATQ = .TRUE., */
/*             LDQ >= MAX(1,M); if UPDATQ = .FALSE., LDQ >= 1. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*) */
/*             On entry, if UPDATZ = .TRUE., then the leading N-by-N */
/*             part of this array must contain a given matrix Z (e.g. */
/*             from a previous call to another SLICOT routine), and on */
/*             exit, the leading N-by-N part of this array contains the */
/*             product of the input matrix Z and the column */
/*             transformation matrix that has transformed the columns of */
/*             the matrices A and E. */
/*             If UPDATZ = .FALSE., the array Z is not referenced and */
/*             can be supplied as a dummy array (i.e. set parameter */
/*             LDZ = 1 and declare this array to be Z(1,1) in the calling */
/*             program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If UPDATZ = .TRUE., */
/*             LDZ >= MAX(1,N); if UPDATZ = .FALSE., LDZ >= 1. */

/*     ISTAIR  (input/output) INTEGER array, dimension (M) */
/*             On entry, ISTAIR contains information on the column */
/*             echelon form of the input matrix E as follows: */
/*             ISTAIR(i) = +j: the boundary element E(i,j) is a corner */
/*                             point; */
/*                         -j: the boundary element E(i,j) is not a */
/*                             corner point (where i=1,...,M). */
/*             On exit, ISTAIR contains the same information for the */
/*             transformed matrix E. */

/*     RANK    (output) INTEGER */
/*             Numerical rank of the submatrix Aj in A (based on TOL). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance used when considering matrix elements */
/*             to be zero. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MB04FZ by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     June 13, 1997, V. Sima. */
/*     November 24, 1997, A. Varga: array starting point A(KK,LL) */
/*                                  correctly set when calling DLASET. */

/*     KEYWORDS */

/*     Echelon form, orthogonal transformation. */

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
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --istair;
    --iwork;

    /* Function Body */
    *rank = 0;
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

/*     Initialisation. */

/*     NJ = number of columns in submatrix Aj, */
/*     MJ = number of rows in submatrices Aj and Ej. */

    nj = *nca;
    mj = *m + 1 - *ifira;
    ifira1 = *ifira - 1;
    ifica1 = *ifica - 1;

    i__1 = nj;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = i__;
/* L20: */
    }

    k = 1;
    lzero = FALSE_;
    *rank = min(nj,mj);
    mxrank = *rank;

/*     WHILE ( K <= MXRANK ) and ( LZERO = FALSE ) DO */
L40:
    if (k <= mxrank && ! lzero) {

/*        Determine column in Aj with largest max-norm. */

	bmxnrm = 0.;
	lsav = k;
	kk = ifira1 + k;

	i__1 = nj;
	for (l = k; l <= i__1; ++l) {

/*           IDAMAX call gives the relative index in column L of Aj where */
/*           max element is found. */
/*           Note: the first element in column L is in row K of */
/*                 matrix Aj. */

	    ll = ifica1 + l;
	    i__2 = mj - k + 1;
	    bmx = (d__1 = a[idamax_(&i__2, &a[kk + ll * a_dim1], &c__1) + kk 
		    - 1 + ll * a_dim1], abs(d__1));
	    if (bmx > bmxnrm) {
		bmxnrm = bmx;
		lsav = l;
	    }
/* L60: */
	}

	ll = ifica1 + k;
	if (bmxnrm < *tol) {

/*           Set submatrix of Aj to zero. */

	    i__1 = mj - k + 1;
	    i__2 = nj - k + 1;
	    dlaset_("Full", &i__1, &i__2, &c_b7, &c_b7, &a[kk + ll * a_dim1], 
		    lda, (ftnlen)4);
	    lzero = TRUE_;
	    *rank = k - 1;
	} else {

/*           Check whether columns have to be interchanged. */

	    if (lsav != k) {

/*              Interchange the columns in A which correspond to the */
/*              columns lsav and k in Aj. Store the permutation in IWORK. */

		dswap_(m, &a[ll * a_dim1 + 1], &c__1, &a[(ifica1 + lsav) * 
			a_dim1 + 1], &c__1);
		ip = iwork[lsav];
		iwork[lsav] = iwork[k];
		iwork[k] = ip;
	    }

	    ++k;
	    mk1 = *n - ll + 1;

	    i__1 = k;
	    for (i__ = mj; i__ >= i__1; --i__) {

/*              II = absolute row number in A corresponding to row i in */
/*                   Aj. */

		ii = ifira1 + i__;

/*              Construct Givens transformation to annihilate Aj(i,k). */
/*              Apply the row transformation to whole matrix A */
/*              (NOT only to Aj). */
/*              Update row transformation matrix Q, if needed. */

		drotg_(&a[ii - 1 + ll * a_dim1], &a[ii + ll * a_dim1], &sc, &
			ss);
		i__2 = mk1 - 1;
		drot_(&i__2, &a[ii - 1 + (ll + 1) * a_dim1], lda, &a[ii + (ll 
			+ 1) * a_dim1], lda, &sc, &ss);
		a[ii + ll * a_dim1] = 0.;
		if (*updatq) {
		    drot_(m, &q[(ii - 1) * q_dim1 + 1], &c__1, &q[ii * q_dim1 
			    + 1], &c__1, &sc, &ss);
		}

/*              Determine boundary type of matrix E at rows II-1 and II. */

		ist1 = istair[ii - 1];
		ist2 = istair[ii];
		if (ist1 * ist2 > 0) {
		    if (ist1 > 0) {

/*                    boundary form = (* x) */
/*                                    (0 *) */

			itype = 1;
		    } else {

/*                    boundary form = (x x) */
/*                                    (x x) */

			itype = 3;
		    }
		} else {
		    if (ist1 < 0) {

/*                    boundary form = (x x) */
/*                                    (* x) */

			itype = 2;
		    } else {

/*                    boundary form = (* x) */
/*                                    (0 x) */

			itype = 4;
		    }
		}

/*              Apply row transformation also to matrix E. */

/*              JC1 = absolute number of the column in E in which stair */
/*                    element of row i-1 of Ej is present. */
/*              JC2 = absolute number of the column in E in which stair */
/*                    element of row i of Ej is present. */

/*              Note: JC1 < JC2   if ITYPE = 1. */
/*                    JC1 = JC2   if ITYPE = 2, 3 or 4. */

		jc1 = abs(ist1);
		jc2 = abs(ist2);
		jpvt = min(jc1,jc2);

		i__2 = *n - jpvt + 1;
		drot_(&i__2, &e[ii - 1 + jpvt * e_dim1], lde, &e[ii + jpvt * 
			e_dim1], lde, &sc, &ss);
		eijpvt = e[ii + jpvt * e_dim1];

		if (itype == 1) {

/*                 Construct column Givens transformation to annihilate */
/*                 E(ii,jpvt). */
/*                 Apply column Givens transformation to matrix E */
/*                 (NOT only to Ej). */

		    drotg_(&e[ii + (jpvt + 1) * e_dim1], &e[ii + jpvt * 
			    e_dim1], &sc, &ss);
		    i__2 = ii - 1;
		    drot_(&i__2, &e[(jpvt + 1) * e_dim1 + 1], &c__1, &e[jpvt *
			     e_dim1 + 1], &c__1, &sc, &ss);
		    e[ii + jpvt * e_dim1] = 0.;

/*                 Apply this transformation also to matrix A */
/*                 (NOT only to Aj). */
/*                 Update column transformation matrix Z, if needed. */

		    drot_(m, &a[(jpvt + 1) * a_dim1 + 1], &c__1, &a[jpvt * 
			    a_dim1 + 1], &c__1, &sc, &ss);
		    if (*updatz) {
			drot_(n, &z__[(jpvt + 1) * z_dim1 + 1], &c__1, &z__[
				jpvt * z_dim1 + 1], &c__1, &sc, &ss);
		    }

		} else if (itype == 2) {
		    if (abs(eijpvt) < *tol) {

/*                                                        (x x)    (* x) */
/*                    Boundary form has been changed from (* x) to (0 x). */

			istpvt = istair[ii];
			istair[ii - 1] = istpvt;
			istair[ii] = -(istpvt + 1);
			e[ii + jpvt * e_dim1] = 0.;
		    }

		} else if (itype == 4) {
		    if (abs(eijpvt) >= *tol) {

/*                                                        (* x)    (x x) */
/*                    Boundary form has been changed from (0 x) to (* x). */

			istpvt = istair[ii - 1];
			istair[ii - 1] = -istpvt;
			istair[ii] = istpvt;
		    }
		}
/* L80: */
	    }

	}
	goto L40;
    }
/*     END WHILE 40 */

/*     Permute columns of Aj to original order. */

    i__1 = ifira1 + *rank;
    dlapmt_(&c_false, &i__1, &nj, &a[*ifica * a_dim1 + 1], lda, &iwork[1]);

    return 0;
/* *** Last line of MB04TT *** */
} /* mb04tt_ */

