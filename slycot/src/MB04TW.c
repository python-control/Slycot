/* MB04TW.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04tw_(logical *updatq, integer *m, integer *n, integer 
	*nre, integer *nce, integer *ifire, integer *ifice, integer *ifica, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, doublereal *
	q, integer *ldq)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j;
    static doublereal sc, ss;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ipvt;
    extern /* Subroutine */ int drotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);


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

/*     To reduce a submatrix E(k) of E to upper triangular form by row */
/*     Givens rotations only. */
/*     Here E(k) = E(IFIRE:me,IFICE:ne), where me = IFIRE - 1 + NRE, */
/*                                             ne = IFICE - 1 + NCE. */
/*     Matrix E(k) is assumed to have full column rank on entry. Hence, */
/*     no pivoting is done during the reduction process. See Algorithm */
/*     2.3.1 and Remark 2.3.4 in [1]. */
/*     The constructed row transformations are also applied to matrix */
/*     A(k) = A(IFIRE:me,IFICA:N). */
/*     Note that in A(k) rows are transformed with the same row indices */
/*     as in E but with column indices different from those in E. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPDATQ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the orthogonal row transformations, as follows: */
/*             = .FALSE.: Do not form Q; */
/*             = .TRUE.:  The given matrix Q is updated by the orthogonal */
/*                        row transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             Number of rows of A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             Number of columns of A and E.  N >= 0. */

/*     NRE     (input) INTEGER */
/*             Number of rows in E to be transformed.  0 <= NRE <= M. */

/*     NCE     (input) INTEGER */
/*             Number of columns in E to be transformed.  0 <= NCE <= N. */

/*     IFIRE   (input) INTEGER */
/*             Index of first row in E to be transformed. */

/*     IFICE   (input) INTEGER */
/*             Index of first column in E to be transformed. */

/*     IFICA   (input) INTEGER */
/*             Index of first column in A to be transformed. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, this array contains the submatrix A(k). */
/*             On exit, it contains the transformed matrix A(k). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, this array contains the submatrix E(k) of full */
/*             column rank to be reduced to upper triangular form. */
/*             On exit, it contains the transformed matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,M). */

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
/*     Supersedes Release 2.0 routine MB04FW by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     June 13, 1997. V. Sima. */
/*     December 30, 1997. A. Varga: Corrected column range to apply */
/*                                  transformations on the matrix E. */

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
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0 || *nre <= 0 || *nce <= 0) {
	return 0;
    }

    ipvt = *ifire - 1;

    i__1 = *ifice + *nce - 1;
    for (j = *ifice; j <= i__1; ++j) {
	++ipvt;

	i__2 = *ifire + *nre - 1;
	for (i__ = ipvt + 1; i__ <= i__2; ++i__) {

/*           Determine the Givens transformation on rows i and ipvt */
/*           to annihilate E(i,j). */
/*           Apply the transformation to these rows (in whole E-matrix) */
/*           from columns j up to n . */
/*           Apply the transformations also to the A-matrix */
/*           (from columns ifica up to n). */
/*           Update the row transformation matrix Q, if needed. */

	    drotg_(&e[ipvt + j * e_dim1], &e[i__ + j * e_dim1], &sc, &ss);
	    i__3 = *n - j;
	    drot_(&i__3, &e[ipvt + (j + 1) * e_dim1], lde, &e[i__ + (j + 1) * 
		    e_dim1], lde, &sc, &ss);
	    e[i__ + j * e_dim1] = 0.;
	    i__3 = *n - *ifica + 1;
	    drot_(&i__3, &a[ipvt + *ifica * a_dim1], lda, &a[i__ + *ifica * 
		    a_dim1], lda, &sc, &ss);
	    if (*updatq) {
		drot_(m, &q[ipvt * q_dim1 + 1], &c__1, &q[i__ * q_dim1 + 1], &
			c__1, &sc, &ss);
	    }
/* L20: */
	}

/* L40: */
    }

    return 0;
/* *** Last line of MB04TW *** */
} /* mb04tw_ */

