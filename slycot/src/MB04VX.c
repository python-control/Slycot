/* MB04VX.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb04vx_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *nblcks, integer *inuk, integer *imuk, doublereal 
	*a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer 
	*ldq, doublereal *z__, integer *ldz, integer *mnei)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, ca, ra;
    static doublereal sc;
    static integer ip;
    static doublereal ss;
    static integer tp1, cja, cje, rje, muk, nuk, mup, nup, mup1, minf, sk1p1, 
	    tk1p1, meps, neps, mukp1;
    extern /* Subroutine */ int mb04tu_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), drotg_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer ismuk, isnuk;


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

/*     To separate the pencils s*E(eps)-A(eps) and s*E(inf)-A(inf) in */
/*     s*E(eps,inf)-A(eps,inf) using Algorithm 3.3.3 in [1]. */

/*     On entry, it is assumed that the M-by-N matrices A and E have */
/*     been obtained after applying the Algorithms 3.2.1 and 3.3.1 to */
/*     the pencil s*E - A as described in [1], i.e. */

/*                        | s*E(eps,inf)-A(eps,inf) |      X      | */
/*        Q'(s*E - A)Z  = |-------------------------|-------------| */
/*                        |             0           | s*E(r)-A(r) | */

/*     Here the pencil s*E(eps,inf)-A(eps,inf) is in staircase form. */
/*     This pencil contains all Kronecker column indices and infinite */
/*     elementary divisors of the pencil s*E - A. */
/*     The pencil s*E(r)-A(r) contains all Kronecker row indices and */
/*     finite elementary divisors of s*E - A. */
/*     Furthermore, the submatrices having full row and column rank in */
/*     the pencil s*E(eps,inf)-A(eps,inf) are assumed to be */
/*     triangularized. */

/*     On exit, the result then is */

/*                        Q'(s*E - A)Z = */

/*          | s*E(eps)-A(eps) |        X        |      X      | */
/*          |-----------------|-----------------|-------------| */
/*          |        0        | s*E(inf)-A(inf) |      X      | */
/*          |===================================|=============| */
/*          |                                   |             | */
/*          |                 0                 | s*E(r)-A(r) | */

/*     Note that the pencil s*E(r)-A(r) is not reduced further. */

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
/*             Number of rows of A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             Number of columns of A and E.  N >= 0. */

/*     NBLCKS  (input) INTEGER */
/*             The number of submatrices having full row rank (possibly */
/*             zero) in A(eps,inf). */

/*     INUK    (input/output) INTEGER array, dimension (NBLCKS) */
/*             On entry, this array contains the row dimensions nu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full row */
/*             rank in the pencil s*E(eps,inf)-A(eps,inf). */
/*             On exit, this array contains the row dimensions nu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full row */
/*             rank in the pencil s*E(eps)-A(eps). */

/*     IMUK    (input/output) INTEGER array, dimension (NBLCKS) */
/*             On entry, this array contains the column dimensions mu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full */
/*             column rank in the pencil s*E(eps,inf)-A(eps,inf). */
/*             On exit, this array contains the column dimensions mu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full */
/*             column rank in the pencil s*E(eps)-A(eps). */

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

/*     MNEI    (output) INTEGER array, dimension (3) */
/*             MNEI(1) = MEPS =    row dimension of sE(eps)-A(eps); */
/*             MNEI(2) = NEPS = column dimension of sE(eps)-A(eps); */
/*             MNEI(3) = MINF = order of the regular pencil */
/*                              sE(inf)-A(inf). */

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Based on Release 3.0 routine MB04TX modified by A. Varga, */
/*     German Aerospace Research Establishment, Oberpfaffenhofen, */
/*     Germany, Nov. 1997, as follows: */
/*     1) NBLCKS is only an input variable; */
/*     2) the significance of MNEI is changed. */

/*     REVISIONS */

/*     A. Varga, DLR Oberpfaffenhofen, March 2002. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, orthogonal */
/*     transformation, staircase form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
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
    --mnei;

    /* Function Body */
    mnei[1] = 0;
    mnei[2] = 0;
    mnei[3] = 0;
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

/*     Initialisation. */

    ismuk = 0;
    isnuk = 0;

    i__1 = *nblcks;
    for (k = 1; k <= i__1; ++k) {
	ismuk += imuk[k];
	isnuk += inuk[k];
/* L20: */
    }

/*     MEPS, NEPS are the dimensions of the pencil s*E(eps)-A(eps). */
/*     MEPS = Sum(k=1,...,nblcks) NU(k), */
/*     NEPS = Sum(k=1,...,nblcks) MU(k). */
/*     MINF is the order of the regular pencil s*E(inf)-A(inf). */

    meps = isnuk;
    neps = ismuk;
    minf = 0;

/*     MUKP1 = mu(k+1).  N.B. It is assumed that mu(NBLCKS + 1) = 0. */

    mukp1 = 0;

    for (k = *nblcks; k >= 1; --k) {
	nuk = inuk[k];
	muk = imuk[k];

/*        Reduce submatrix E(k,k+1) to square matrix. */
/*        NOTE that always NU(k) >= MU(k+1) >= 0. */

/*        WHILE ( NU(k) >  MU(k+1) ) DO */
L40:
	if (nuk > mukp1) {

/*           sk1p1 = sum(i=k+1,...,p-1) NU(i) */
/*           tk1p1 = sum(i=k+1,...,p-1) MU(i) */
/*           ismuk = sum(i=1,...,k) MU(i) */
/*           tp1   = sum(i=1,...,p-1) MU(i) = ismuk + tk1p1. */

	    sk1p1 = 0;
	    tk1p1 = 0;

	    i__1 = *nblcks;
	    for (ip = k + 1; ip <= i__1; ++ip) {

/*              Annihilate the elements originally present in the last */
/*              row of E(k,p+1) and A(k,p). */
/*              Start annihilating the first MU(p) - MU(p+1) elements by */
/*              applying column Givens rotations plus interchanging */
/*              elements. */
/*              Use original bottom diagonal element of A(k,k) as pivot. */
/*              Start position of pivot in A = (ra,ca). */

		tp1 = ismuk + tk1p1;
		ra = isnuk + sk1p1;
		ca = tp1;

		mup = imuk[ip];
		nup = inuk[ip];
		mup1 = nup;

		i__2 = ca + mup - nup - 1;
		for (cja = ca; cja <= i__2; ++cja) {

/*                 CJA = current column index of pivot in A. */

		    drotg_(&a[ra + cja * a_dim1], &a[ra + (cja + 1) * a_dim1],
			     &sc, &ss);

/*                 Apply transformations to A- and E-matrix. */
/*                 Interchange columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

		    i__3 = ra - 1;
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
		    a[ra + (cja + 1) * a_dim1] = a[ra + cja * a_dim1];
		    a[ra + cja * a_dim1] = 0.;
		    mb04tu_(&ra, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
		    if (*updatz) {
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
		    }
/* L60: */
		}

/*              Annihilate the remaining elements originally present in */
/*              the last row of E(k,p+1) and A(k,p) by alternatingly */
/*              applying row and column rotations plus interchanging */
/*              elements. */
/*              Use diagonal elements of E(p,p+1) and original bottom */
/*              diagonal element of A(k,k) as pivots, respectively. */
/*              (re,ce) and (ra,ca) are the starting positions of the */
/*              pivots in E and A. */

		cje = tp1 + mup;
		cja = cje - mup1 - 1;

		i__2 = ra + mup1;
		for (rje = ra + 1; rje <= i__2; ++rje) {

/*                 (RJE,CJE) = current position pivot in E. */

		    ++cje;
		    ++cja;

/*                 Determine the row transformations. */
/*                 Apply these transformations to E- and A-matrix. */
/*                 Interchange the rows simultaneously. */
/*                 Update row transformation matrix Q, if needed. */

		    drotg_(&e[rje + cje * e_dim1], &e[rje - 1 + cje * e_dim1],
			     &sc, &ss);
		    i__3 = *n - cje;
		    mb04tu_(&i__3, &e[rje + (cje + 1) * e_dim1], lde, &e[rje 
			    - 1 + (cje + 1) * e_dim1], lde, &sc, &ss);
		    e[rje - 1 + cje * e_dim1] = e[rje + cje * e_dim1];
		    e[rje + cje * e_dim1] = 0.;
		    i__3 = *n - cja + 1;
		    mb04tu_(&i__3, &a[rje + cja * a_dim1], lda, &a[rje - 1 + 
			    cja * a_dim1], lda, &sc, &ss);
		    if (*updatq) {
			mb04tu_(m, &q[rje * q_dim1 + 1], &c__1, &q[(rje - 1) *
				 q_dim1 + 1], &c__1, &sc, &ss);
		    }

/*                 Determine the column transformations. */
/*                 Apply these transformations to A- and E-matrix. */
/*                 Interchange the columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

		    drotg_(&a[rje + cja * a_dim1], &a[rje + (cja + 1) * 
			    a_dim1], &sc, &ss);
		    i__3 = rje - 1;
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
		    a[rje + (cja + 1) * a_dim1] = a[rje + cja * a_dim1];
		    a[rje + cja * a_dim1] = 0.;
		    mb04tu_(&rje, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
		    if (*updatz) {
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
		    }
/* L80: */
		}

		sk1p1 += nup;
		tk1p1 += mup;

/* L100: */
	    }

/*           Reduce A=A(eps,inf) and E=E(eps,inf) by ignoring their last */
/*           row and right most column. The row and column ignored */
/*           belong to the pencil s*E(inf)-A(inf). */
/*           Redefine blocks in new A and E. */

	    --muk;
	    --nuk;
	    --ismuk;
	    --isnuk;
	    --meps;
	    --neps;
	    ++minf;

	    goto L40;
	}
/*        END WHILE 40 */

	imuk[k] = muk;
	inuk[k] = nuk;

/*        Now submatrix E(k,k+1) is square. */

/*        Consider next submatrix (k:=k-1). */

	isnuk -= nuk;
	ismuk -= muk;
	mukp1 = muk;
/* L120: */
    }

/*     Store dimensions of the pencils s*E(eps)-A(eps) and */
/*     s*E(inf)-A(inf) in array MNEI. */

    mnei[1] = meps;
    mnei[2] = neps;
    mnei[3] = minf;

    return 0;
/* *** Last line of MB04VX *** */
} /* mb04vx_ */

