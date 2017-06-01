/* MB02NY.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02ny_(logical *updatu, logical *updatv, integer *m, 
	integer *n, integer *i__, integer *k, doublereal *q, doublereal *e, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *
	dwork)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal c__, f, g;
    static integer l;
    static doublereal r__, s;
    static integer i1, l1, irot, nrot;
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


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

/*     To separate a zero singular value of a bidiagonal submatrix of */
/*     order k, k <= p, of the bidiagonal matrix */

/*               |Q(1) E(1)  0    ...   0   | */
/*               | 0   Q(2) E(2)        .   | */
/*           J = | .                    .   | */
/*               | .                  E(p-1)| */
/*               | 0   ...  ...   ...  Q(p) | */

/*     with p = MIN(M,N), by annihilating one or two superdiagonal */
/*     elements E(i-1) (if i > 1) and/or E(i) (if i < k). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPDATU  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the left-hand Givens rotations S, as follows: */
/*             = .FALSE.:  Do not form U; */
/*             = .TRUE. :  The given matrix U is updated (postmultiplied) */
/*                         by the left-hand Givens rotations S. */

/*     UPDATV  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the right-hand Givens rotations T, as follows: */
/*             = .FALSE.:  Do not form V; */
/*             = .TRUE. :  The given matrix V is updated (postmultiplied) */
/*                         by the right-hand Givens rotations T. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix U.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix V.  N >= 0. */

/*     I       (input) INTEGER */
/*             The index of the negligible diagonal entry Q(I) of the */
/*             bidiagonal matrix J, I <= p. */

/*     K       (input) INTEGER */
/*             The index of the last diagonal entry of the considered */
/*             bidiagonal submatrix of J, i.e., E(K-1) is considered */
/*             negligible, K <= p. */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (p) */
/*             where p = MIN(M,N). */
/*             On entry, Q must contain the diagonal entries of the */
/*             bidiagonal matrix J. */
/*             On exit, Q contains the diagonal entries of the */
/*             transformed bidiagonal matrix S' J T. */

/*     E       (input/output) DOUBLE PRECISION array, dimension (p-1) */
/*             On entry, E must contain the superdiagonal entries of J. */
/*             On exit, E contains the superdiagonal entries of the */
/*             transformed bidiagonal matrix S' J T. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,p) */
/*             On entry, if UPDATU = .TRUE., U must contain the M-by-p */
/*             left transformation matrix. */
/*             On exit, if UPDATU = .TRUE., the Givens rotations S on the */
/*             left, annihilating E(i) if i < k, have been postmultiplied */
/*             into U. */
/*             U is not referenced if UPDATU = .FALSE.. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= max(1,M) if UPDATU = .TRUE.; */
/*             LDU >= 1        if UPDATU = .FALSE.. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,p) */
/*             On entry, if UPDATV = .TRUE., V must contain the N-by-p */
/*             right transformation matrix. */
/*             On exit, if UPDATV = .TRUE., the Givens rotations T on the */
/*             right, annihilating E(i-1) if i > 1,  have been */
/*             postmultiplied into V. */
/*             V is not referenced if UPDATV = .FALSE.. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= max(1,N) if UPDATV = .TRUE.; */
/*             LDV >= 1        if UPDATV = .FALSE.. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(1,LDWORK)) */
/*             LDWORK >= 2*MAX(K-I,I-1),  if UPDATV = UPDATU = .TRUE.; */
/*             LDWORK >= 2*(K-I), if UPDATU = .TRUE., UPDATV = .FALSE.; */
/*             LDWORK >= 2*(I-1), if UPDATV = .TRUE., UPDATU = .FALSE.; */
/*             LDWORK >= 1,       if UPDATU = UPDATV = .FALSE.. */

/*     METHOD */

/*     Let the considered bidiagonal submatrix be */

/*               |Q(1) E(1)  0                    ...   0   | */
/*               | 0   Q(2) E(2)                        .   | */
/*               | .                                    .   | */
/*               | .           Q(i-1) E(i-1)            .   | */
/*          Jk = | .                   Q(i) E(i)        .   |. */
/*               | .                       Q(i+1) .     .   | */
/*               | .                              ..    .   | */
/*               | .                                  E(k-1)| */
/*               | 0    ...                       ...  Q(k) | */

/*     A zero singular value of Jk manifests itself by a zero diagonal */
/*     entry Q(i) or in practice, a negligible value of Q(i). */
/*     When a negligible diagonal element Q(i) in Jk is present, the */
/*     bidiagonal submatrix Jk is split by the routine into 2 or 3 */
/*     unreduced bidiagonal submatrices by annihilating E(i) (if i < k) */
/*     using Givens rotations S on the left and by annihilating E(i-1) */
/*     (if i > 1) using Givens rotations T on the right until Jk is */
/*     reduced to the form: */

/*               |Q(1) E(1)  0                ...   0   | */
/*               | 0         .                ...   .   | */
/*               | .                          ...   .   | */
/*               | .       Q(i-1) 0                 .   | */
/*     S' Jk T = | .              0   0             .   |. */
/*               | .                 Q(i+1)   .     .   | */
/*               | .                          ..    .   | */
/*               | .                              E(k-1)| */
/*               | 0    ...                   ...  Q(k) | */

/*     For more details, see [1, pp.11.12-11.14]. */

/*     REFERENCES */

/*     [1] Dongarra, J.J., Bunch, J.R., Moler C.B. and Stewart, G.W. */
/*         LINPACK User's Guide. */
/*         SIAM, Philadelphia, 1979. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB02BZ by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bidiagonal matrix, orthogonal transformation, singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     For speed, no tests of the input scalar arguments are done. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    --q;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --dwork;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (*i__ <= min(*m,*n)) {
	q[*i__] = 0.;
    }

/*     Annihilate E(I) (if I < K). */

    if (*i__ < *k) {
	c__ = 0.;
	s = 1.;
	irot = 0;
	nrot = *k - *i__;

	i__1 = *k - 1;
	for (l = *i__; l <= i__1; ++l) {
	    g = e[l];
	    e[l] = c__ * g;
	    d__1 = s * g;
	    dlartg_(&q[l + 1], &d__1, &c__, &s, &r__);
	    q[l + 1] = r__;
	    if (*updatu) {
		++irot;
		dwork[irot] = c__;
		dwork[irot + nrot] = s;
	    }
/* L20: */
	}

	if (*updatu) {
	    i__1 = nrot + 1;
	    dlasr_("Right", "Top", "Forward", m, &i__1, &dwork[1], &dwork[
		    nrot + 1], &u[*i__ * u_dim1 + 1], ldu, (ftnlen)5, (ftnlen)
		    3, (ftnlen)7);
	}
    }

/*     Annihilate E(I-1) (if I > 1). */

    if (*i__ > 1) {
	i1 = *i__ - 1;
	f = e[i1];
	e[i1] = 0.;

	i__1 = i1 - 1;
	for (l1 = 1; l1 <= i__1; ++l1) {
	    l = *i__ - l1;
	    dlartg_(&q[l], &f, &c__, &s, &r__);
	    q[l] = r__;
	    if (*updatv) {
		dwork[l] = c__;
		dwork[l + i1] = s;
	    }
	    g = e[l - 1];
	    f = -s * g;
	    e[l - 1] = c__ * g;
/* L40: */
	}

	dlartg_(&q[1], &f, &c__, &s, &r__);
	q[1] = r__;
	if (*updatv) {
	    dwork[1] = c__;
	    dwork[*i__] = s;
	    dlasr_("Right", "Bottom", "Backward", n, i__, &dwork[1], &dwork[*
		    i__], &v[v_dim1 + 1], ldv, (ftnlen)5, (ftnlen)6, (ftnlen)
		    8);
	}
    }

    return 0;
/* *** Last line of MB02NY *** */
} /* mb02ny_ */

