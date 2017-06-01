/* MB04TV.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04tv_(logical *updatz, integer *n, integer *nra, 
	integer *nca, integer *ifira, integer *ifica, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *z__, integer *ldz)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, z_dim1, z_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j;
    static doublereal sc, ss;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jpvt;
    extern /* Subroutine */ int drotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);
    static integer ifira1;


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

/*     To reduce a submatrix A(k) of A to upper triangular form by column */
/*     Givens rotations only. */
/*     Here A(k) = A(IFIRA:ma,IFICA:na) where ma = IFIRA - 1 + NRA, */
/*     na = IFICA - 1 + NCA. */
/*     Matrix A(k) is assumed to have full row rank on entry. Hence, no */
/*     pivoting is done during the reduction process. See Algorithm 2.3.1 */
/*     and Remark 2.3.4 in [1]. */
/*     The constructed column transformations are also applied to matrix */
/*     E(k) = E(1:IFIRA-1,IFICA:na). */
/*     Note that in E columns are transformed with the same column */
/*     indices as in A, but with row indices different from those in A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPDATZ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal column transformations, as */
/*             follows: */
/*             = .FALSE.: Do not form Z; */
/*             = .TRUE.:  The given matrix Z is updated by the orthogonal */
/*                        column transformations used in the reduction. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             Number of columns of A and E.  N >= 0. */

/*     NRA     (input) INTEGER */
/*             Number of rows in A to be transformed.  0 <= NRA <= LDA. */

/*     NCA     (input) INTEGER */
/*             Number of columns in A to be transformed.  0 <= NCA <= N. */

/*     IFIRA   (input) INTEGER */
/*             Index of the first row in A to be transformed. */

/*     IFICA   (input) INTEGER */
/*             Index of the first column in A to be transformed. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the elements of A(IFIRA:ma,IFICA:na) must */
/*             contain the submatrix A(k) of full row rank to be reduced */
/*             to upper triangular form. */
/*             On exit, it contains the transformed matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,NRA). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the elements of E(1:IFIRA-1,IFICA:na) must */
/*             contain the submatrix E(k). */
/*             On exit, it contains the transformed matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,IFIRA-1). */

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

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB04FV by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, orthogonal transformation, */
/*     staircase form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    if (*n <= 0 || *nra <= 0 || *nca <= 0) {
	return 0;
    }
    ifira1 = *ifira - 1;
    jpvt = *ifica + *nca;

    i__1 = *ifira;
    for (i__ = ifira1 + *nra; i__ >= i__1; --i__) {
	--jpvt;

	i__2 = *ifica;
	for (j = jpvt - 1; j >= i__2; --j) {

/*           Determine the Givens transformation on columns j and jpvt */
/*           to annihilate A(i,j). Apply the transformation to these */
/*           columns from rows 1 up to i. */
/*           Apply the transformation also to the E-matrix (from rows 1 */
/*           up to ifira1). */
/*           Update column transformation matrix Z, if needed. */

	    drotg_(&a[i__ + jpvt * a_dim1], &a[i__ + j * a_dim1], &sc, &ss);
	    i__3 = i__ - 1;
	    drot_(&i__3, &a[jpvt * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1, &sc, &ss);
	    a[i__ + j * a_dim1] = 0.;
	    drot_(&ifira1, &e[jpvt * e_dim1 + 1], &c__1, &e[j * e_dim1 + 1], &
		    c__1, &sc, &ss);
	    if (*updatz) {
		drot_(n, &z__[jpvt * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1, &sc, &ss);
	    }
/* L20: */
	}

/* L40: */
    }

    return 0;
/* *** Last line of MB04TV *** */
} /* mb04tv_ */

