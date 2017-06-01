/* MB04TY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04ty_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *nblcks, integer *inuk, integer *imuk, doublereal 
	*a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer 
	*ldq, doublereal *z__, integer *ldz, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;

    /* Local variables */
    static integer k, muk, nuk, mukp1, ifica, ifice, ifire;
    extern /* Subroutine */ int mb04tv_(logical *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), mb04tw_(
	    logical *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer ismuk, isnuk1;


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

/*     To perform the triangularization of the submatrices having full */
/*     row and column rank in the pencil s*E(eps,inf)-A(eps,inf) below */

/*                    | s*E(eps,inf)-A(eps,inf) |     X       | */
/*          s*E - A = |-------------------------|-------------| , */
/*                    |            0            | s*E(r)-A(r) | */

/*     using Algorithm 3.3.1 in [1]. */
/*     On entry, it is assumed that the M-by-N matrices A and E have */
/*     been transformed to generalized Schur form by unitary */
/*     transformations (see Algorithm 3.2.1 in [1]), and that the pencil */
/*     s*E(eps,inf)-A(eps,inf) is in staircase form. */
/*     This pencil contains all Kronecker column indices and infinite */
/*     elementary divisors of the pencil s*E - A. */
/*     The pencil s*E(r)-A(r) contains all Kronecker row indices and */
/*     finite elementary divisors of s*E - A. */

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
/*             Number of rows in A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             Number of columns in A and E.  N >= 0. */

/*     NBLCKS  (input) INTEGER */
/*             Number of submatrices having full row rank (possibly zero) */
/*             in A(eps,inf). */

/*     INUK    (input) INTEGER array, dimension (NBLCKS) */
/*             The row dimensions nu(k) (k=1, 2, ..., NBLCKS) of the */
/*             submatrices having full row rank in the pencil */
/*             s*E(eps,inf)-A(eps,inf). */

/*     IMUK    (input) INTEGER array, dimension (NBLCKS) */
/*             The column dimensions mu(k) (k=1, 2, ..., NBLCKS) of the */
/*             submatrices having full column rank in the pencil. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, this array contains the matrix A to be reduced. */
/*             On exit, it contains the transformed matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, this array contains the matrix E to be reduced. */
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

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if incorrect dimensions of a full column rank */
/*                   submatrix; */
/*             = 2:  if incorrect dimensions of a full row rank */
/*                   submatrix. */

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
/*     Supersedes Release 2.0 routine MB04FY by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, orthogonal transformation, */
/*     staircase form. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --inuk;
    --imuk;
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

    /* Function Body */
    *info = 0;
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

/*     ISMUK  = sum(i=1,...,k) MU(i), */
/*     ISNUK1 = sum(i=1,...,k-1) NU(i). */

    ismuk = 0;
    isnuk1 = 0;

    i__1 = *nblcks;
    for (k = 1; k <= i__1; ++k) {
	ismuk += imuk[k];
	isnuk1 += inuk[k];
/* L20: */
    }

/*     Note:  ISNUK1 has not yet the correct value. */

    mukp1 = 0;

    for (k = *nblcks; k >= 1; --k) {
	muk = imuk[k];
	nuk = inuk[k];
	isnuk1 -= nuk;

/*        Determine left upper absolute co-ordinates of E(k) in E-matrix */
/*        and of A(k) in A-matrix. */

	ifire = isnuk1 + 1;
	ifice = ismuk + 1;
	ifica = ifice - muk;

/*        Reduce E(k) to upper triangular form using Givens */
/*        transformations on rows only. Apply the same transformations */
/*        to the rows of A(k). */

	if (mukp1 > nuk) {
	    *info = 1;
	    return 0;
	}

	mb04tw_(updatq, m, n, &nuk, &mukp1, &ifire, &ifice, &ifica, &a[
		a_offset], lda, &e[e_offset], lde, &q[q_offset], ldq);

/*        Reduce A(k) to upper triangular form using Givens */
/*        transformations on columns only. Apply the same transformations */
/*        to the columns in the E-matrix. */

	if (nuk > muk) {
	    *info = 2;
	    return 0;
	}

	mb04tv_(updatz, n, &nuk, &muk, &ifire, &ifica, &a[a_offset], lda, &e[
		e_offset], lde, &z__[z_offset], ldz);

	ismuk -= muk;
	mukp1 = muk;
/* L40: */
    }

    return 0;
/* *** Last line of MB04TY *** */
} /* mb04ty_ */

