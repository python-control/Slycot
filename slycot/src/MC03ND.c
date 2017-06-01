/* MC03ND.f -- translated by f2c (version 20100827).
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

static doublereal c_b13 = 0.;
static doublereal c_b19 = 1.;

/* Subroutine */ int mc03nd_(integer *mp, integer *np, integer *dp, 
	doublereal *p, integer *ldp1, integer *ldp2, integer *dk, integer *
	gam, doublereal *nullsp, integer *ldnull, doublereal *ker, integer *
	ldker1, integer *ldker2, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer ker_dim1, ker_dim2, ker_offset, nullsp_dim1, nullsp_offset, 
	    p_dim1, p_dim2, p_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer h__, i__, j, k, m, n, vc1, vr2, nca, nra, ncv, muk, nuk, 
	    gamj, mnei[3], ifir, tail, idiff;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb04ud_(char *, char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), mb04vd_(char *, char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    static integer ranke, sgamk;
    extern /* Subroutine */ int mc03nx_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *), mc03ny_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal toler;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    static integer nblcki, nblcks;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer jworka, jworke, jworkq, jworkv, jworkz;


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

/*     To compute the coefficients of a minimal polynomial basis */
/*                                                 DK */
/*         K(s) = K(0) + K(1) * s + ... + K(DK) * s */

/*     for the right nullspace of the MP-by-NP polynomial matrix of */
/*     degree DP, given by */
/*                                                 DP */
/*         P(s) = P(0) + P(1) * s + ... + P(DP) * s  , */

/*     which corresponds to solving the polynomial matrix equation */
/*     P(s) * K(s) = 0. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     MP      (input) INTEGER */
/*             The number of rows of the polynomial matrix P(s). */
/*             MP >= 0. */

/*     NP      (input) INTEGER */
/*             The number of columns of the polynomial matrix P(s). */
/*             NP >= 0. */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial matrix P(s).  DP >= 1. */

/*     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array must */
/*             contain the coefficients of the polynomial matrix P(s). */
/*             Specifically, P(i,j,k) must contain the (i,j)-th element */
/*             of P(k-1), which is the cofficient of s**(k-1) of P(s), */
/*             where i = 1,2,...,MP, j = 1,2,...,NP and k = 1,2,...,DP+1. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MAX(1,MP). */

/*     LDP2    INTEGER */
/*             The second dimension of array P.   LDP2 >= MAX(1,NP). */

/*     DK      (output) INTEGER */
/*             The degree of the minimal polynomial basis K(s) for the */
/*             right nullspace of P(s) unless DK = -1, in which case */
/*             there is no right nullspace. */

/*     GAM     (output) INTEGER array, dimension (DP*MP+1) */
/*             The leading (DK+1) elements of this array contain */
/*             information about the ordering of the right nullspace */
/*             vectors stored in array NULLSP. */

/*     NULLSP  (output) DOUBLE PRECISION array, dimension */
/*             (LDNULL,(DP*MP+1)*NP) */
/*             The leading NP-by-SUM(i*GAM(i)) part of this array */
/*             contains the right nullspace vectors of P(s) in condensed */
/*             form (as defined in METHOD), where i = 1,2,...,DK+1. */

/*     LDNULL  INTEGER */
/*             The leading dimension of array NULLSP. */
/*             LDNULL >= MAX(1,NP). */

/*     KER     (output) DOUBLE PRECISION array, dimension */
/*             (LDKER1,LDKER2,DP*MP+1) */
/*             The leading NP-by-nk-by-(DK+1) part of this array contains */
/*             the coefficients of the minimal polynomial basis K(s), */
/*             where nk = SUM(GAM(i)) and i = 1,2,...,DK+1. Specifically, */
/*             KER(i,j,m) contains the (i,j)-th element of K(m-1), which */
/*             is the coefficient of s**(m-1) of K(s), where i = 1,2,..., */
/*             NP, j = 1,2,...,nk and m = 1,2,...,DK+1. */

/*     LDKER1  INTEGER */
/*             The leading dimension of array KER.  LDKER1 >= MAX(1,NP). */

/*     LDKER2  INTEGER */
/*             The second dimension of array KER.   LDKER2 >= MAX(1,NP). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance below which matrix elements are considered */
/*             to be zero. If the user sets TOL to be less than */
/*             10 * EPS * MAX( ||A|| , ||E|| ), then the tolerance is */
/*                                  F       F */
/*             taken as 10 * EPS * MAX( ||A|| , ||E|| ), where EPS is the */
/*                                           F       F */
/*             machine precision (see LAPACK Library Routine DLAMCH) and */
/*             A and E are matrices (as defined in METHOD). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (m+2*MAX(n,m+1)+n), */
/*             where m = DP*MP and n = (DP-1)*MP + NP. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= m*n*n + 2*m*n + 2*n*n. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*             > 0:  if incorrect rank decisions were taken during the */
/*                   computations. This failure is not likely to occur. */
/*                   The possible values are: */
/*                     k, 1 <= k <= DK+1, the k-th diagonal submatrix had */
/*                           not a full row rank; */
/*                     DK+2, if incorrect dimensions of a full column */
/*                           rank submatrix; */
/*                     DK+3, if incorrect dimensions of a full row rank */
/*                           submatrix. */

/*     METHOD */

/*     The computation of the right nullspace of the MP-by-NP polynomial */
/*     matrix P(s) of degree DP given by */
/*                                                  DP-1            DP */
/*        P(s) = P(0) + P(1) * s + ... + P(DP-1) * s     + P(DP) * s */

/*     is performed via the pencil s*E - A, associated with P(s), where */

/*            | I              |           | 0         -P(DP) | */
/*            |   .            |           | I .          .   | */
/*        A = |     .          |  and  E = |   . .        .   |.      (1) */
/*            |       .        |           |     . 0      .   | */
/*            |         I      |           |       I 0 -P(2)  | */
/*            |           P(0) |           |         I -P(1)  | */

/*     The pencil s*E - A is transformed by unitary matrices Q and Z such */
/*     that */

/*                     | sE(eps)-A(eps) |        X       |      X     | */
/*                     |----------------|----------------|------------| */
/*                     |        0       | sE(inf)-A(inf) |      X     | */
/*        Q'(s*E-A)Z = |=================================|============|. */
/*                     |                                 |            | */
/*                     |                0                | sE(r)-A(r) | */

/*     Since s*E(inf)-A(inf) and s*E(r)-A(r) have full column rank, the */
/*     minimal polynomial basis for the right nullspace of Q'(s*E-A)Z */
/*     (and consequently the basis for the right nullspace of s*E - A) is */
/*     completely determined by s*E(eps)-A(eps). */

/*     Let Veps(s) be a minimal polynomial basis for the right nullspace */
/*     of s*E(eps)-A(eps). Then */

/*                   | Veps(s) | */
/*        V(s) = Z * |---------| */
/*                   |    0    | */

/*     is a minimal polynomial basis for the right nullspace of s*E - A. */
/*     From the structure of s*E - A it can be shown that if V(s) is */
/*     partitioned as */

/*               | Vo(s) | (DP-1)*MP */
/*        V(s) = |------ | */
/*               | Ve(s) | NP */

/*     then the columns of Ve(s) form a minimal polynomial basis for the */
/*     right nullspace of P(s). */

/*     The vectors of Ve(s) are computed and stored in array NULLSP in */
/*     the following condensed form: */

/*        ||      ||      |      ||      |      |      ||      |     | */
/*        || U1,0 || U2,0 | U2,1 || U3,0 | U3,1 | U3,2 || U4,0 | ... |, */
/*        ||      ||      |      ||      |      |      ||      |     | */

/*     where Ui,j is an NP-by-GAM(i) matrix which contains the i-th block */
/*     of columns of K(j), the j-th coefficient of the polynomial matrix */
/*     representation for the right nullspace */
/*                                                  DK */
/*        K(s) = K(0) + K(1) * s + . . . + K(DK) * s  . */

/*     The coefficients K(0), K(1), ..., K(DK) are NP-by-nk matrices */
/*     given by */

/*        K(0)  = | U1,0 | U2,0 | U3,0 | . . .          | U(DK+1,0) | */

/*        K(1)  = |  0   | U2,1 | U3,1 | . . .          | U(DK+1,1) | */

/*        K(2)  = |  0   |  0   | U3,2 | . . .          | U(DK+1,2) | */

/*          .     .     .     .     .     .     .     .     .     . */

/*        K(DK) = |  0   |  0   |  0   | . . .    |  0  | U(DK+1,DK)|. */

/*     Note that the degree of K(s) satisfies the inequality DK <= */
/*     DP * MIN(MP,NP) and that the dimension of K(s) satisfies the */
/*     inequality (NP-MP) <= nk <= NP. */

/*     REFERENCES */

/*     [1] Beelen, Th.G.J. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, 1987. */

/*     [2] Van Den Hurk, G.J.H.H. */
/*         New Algorithms for Solving Polynomial Matrix Problems. */
/*         Master's Thesis, Eindhoven University of Technology, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm used by the routine involves the construction of a */
/*     special block echelon form with pivots considered to be non-zero */
/*     when they are larger than TOL. These pivots are then inverted in */
/*     order to construct the columns of the kernel of the polynomial */
/*     matrix. If TOL is chosen to be too small then these inversions may */
/*     be sensitive whereas increasing TOL will make the inversions more */
/*     robust but will affect the block echelon form (and hence the */
/*     column degrees of the polynomial kernel). Furthermore, if the */
/*     elements of the computed polynomial kernel are large relative to */
/*     the polynomial matrix, then the user should consider trying */
/*     several values of TOL. */

/*     FURTHER COMMENTS */

/*     It also possible to compute a minimal polynomial basis for the */
/*     right nullspace of a pencil, since a pencil is a polynomial matrix */
/*     of degree 1. Thus for the pencil (s*E - A), the required input is */
/*     P(1)  = E and P(0) = -A. */

/*     The routine can also be used to compute a minimal polynomial */
/*     basis for the left nullspace of a polynomial matrix by simply */
/*     transposing P(s). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC03BD by A.J. Geurts and MC03BZ by */
/*     Th.G.J. Beelen, A.J. Geurts, and G.J.H.H. van den Hurk. */

/*     REVISIONS */

/*     Jan. 1998. */

/*     KEYWORDS */

/*     Echelon form, elementary polynomial operations, input output */
/*     description, polynomial matrix, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    p_dim1 = *ldp1;
    p_dim2 = *ldp2;
    p_offset = 1 + p_dim1 * (1 + p_dim2);
    p -= p_offset;
    --gam;
    nullsp_dim1 = *ldnull;
    nullsp_offset = 1 + nullsp_dim1;
    nullsp -= nullsp_offset;
    ker_dim1 = *ldker1;
    ker_dim2 = *ldker2;
    ker_offset = 1 + ker_dim1 * (1 + ker_dim2);
    ker -= ker_offset;
    --iwork;
    --dwork;

    /* Function Body */
    m = *dp * *mp;
    h__ = m - *mp;
    n = h__ + *np;
    *info = 0;
    if (*mp < 0) {
	*info = -1;
    } else if (*np < 0) {
	*info = -2;
    } else if (*dp <= 0) {
	*info = -3;
    } else if (*ldp1 < max(1,*mp)) {
	*info = -5;
    } else if (*ldp2 < max(1,*np)) {
	*info = -6;
    } else if (*ldnull < max(1,*np)) {
	*info = -10;
    } else if (*ldker1 < max(1,*np)) {
	*info = -12;
    } else if (*ldker2 < max(1,*np)) {
	*info = -13;
    } else if (*ldwork < n * (m * n + (m + n << 1))) {
	*info = -17;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC03ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*mp == 0 || *np == 0) {
	*dk = -1;
	return 0;
    }

    jworka = 1;
    jworke = jworka + m * n;
    jworkz = jworke + m * n;
    jworkv = jworkz + n * n;
    jworkq = jworka;

/*     Construct the matrices A and E in the pencil s*E-A in (1). */
/*     Workspace:  2*M*N. */

    mc03nx_(mp, np, dp, &p[p_offset], ldp1, ldp2, &dwork[jworka], &m, &dwork[
	    jworke], &m);

/*     Computation of the tolerance. */

/* Computing MAX */
    d__1 = dlange_("F", &m, np, &dwork[jworke + h__ * m], &m, &dwork[1], (
	    ftnlen)1), d__2 = dlange_("F", mp, np, &p[p_offset], ldp1, &dwork[
	    1], (ftnlen)1);
    toler = max(d__1,d__2);
    d__1 = sqrt((doublereal) h__);
    toler = dlamch_("Epsilon", (ftnlen)7) * 10. * dlapy2_(&toler, &d__1);
    if (toler <= *tol) {
	toler = *tol;
    }

/*     Reduction of E to column echelon form E0 = Q' x E x Z and */
/*     transformation of A, A0 = Q' x A x Z. */
/*     Workspace:  2*M*N + N*N + max(M,N). */

    mb04ud_("No Q", "Identity Z", &m, &n, &dwork[jworka], &m, &dwork[jworke], 
	    &m, &dwork[jworkq], &m, &dwork[jworkz], &n, &ranke, &iwork[1], &
	    toler, &dwork[jworkv], info, (ftnlen)4, (ftnlen)10);

/*     The contents of ISTAIR is transferred from MB04UD to MB04VD by */
/*     IWORK(i), i=1,...,M. */
/*     In the sequel the arrays IMUK and INUK are part of IWORK, namely: */
/*     IWORK(i), i = M+1,...,M+max(N,M+1), contains IMUK, */
/*     IWORK(i), i = M+max(N,M+1)+1,...,M+2*max(N,M+1), contains INUK. */
/*     IWORK(i), i = M+2*max(N,M+1)+1,...,M+2*max(N,M+1)+N, contains */
/*               IMUK0 (not needed), and is also used as workspace. */

    muk = m + 1;
/* Computing MAX */
    i__1 = n, i__2 = m + 1;
    nuk = muk + max(i__1,i__2);
/* Computing MAX */
    i__1 = n, i__2 = m + 1;
    tail = nuk + max(i__1,i__2);

    mb04vd_("Separation", "No Q", "Update Z", &m, &n, &ranke, &dwork[jworka], 
	    &m, &dwork[jworke], &m, &dwork[jworkq], &m, &dwork[jworkz], &n, &
	    iwork[1], &nblcks, &nblcki, &iwork[muk], &iwork[nuk], &iwork[tail]
	    , mnei, &toler, &iwork[tail], info, (ftnlen)10, (ftnlen)4, (
	    ftnlen)8);
    if (*info > 0) {

/*        Incorrect rank decisions. */

	*info += nblcks;
	return 0;
    }

/*     If NBLCKS < 1, or the column dimension of s*E(eps) - A(eps) is */
/*     zero, then there is no right nullspace. */

    if (nblcks < 1 || mnei[1] == 0) {
	*dk = -1;
	return 0;
    }

/*     Start of the computation of the minimal basis. */

    *dk = nblcks - 1;
    nra = mnei[0];
    nca = mnei[1];

/*     Determine a minimal basis VEPS(s) for the right nullspace of the */
/*     pencil s*E(eps)-A(eps) associated with the polynomial matrix P(s). */
/*     Workspace:  2*M*N + N*N + N*N*(M+1). */

    mc03ny_(&nblcks, &nra, &nca, &dwork[jworka], &m, &dwork[jworke], &m, &
	    iwork[muk], &iwork[nuk], &dwork[jworkv], &n, info);

    if (*info > 0) {
	return 0;
    }

    ncv = iwork[muk] - iwork[nuk];
    gam[1] = ncv;
    iwork[1] = 0;
    iwork[tail] = iwork[muk];

    i__1 = nblcks;
    for (i__ = 2; i__ <= i__1; ++i__) {
	idiff = iwork[muk + i__ - 1] - iwork[nuk + i__ - 1];
	gam[i__] = idiff;
	iwork[i__] = ncv;
	ncv += i__ * idiff;
	iwork[tail + i__ - 1] = iwork[tail + i__ - 2] + iwork[muk + i__ - 1];
/* L20: */
    }

/*     Determine a basis for the right nullspace of the polynomial */
/*     matrix P(s). This basis is stored in array NULLSP in condensed */
/*     form. */

    dlaset_("Full", np, &ncv, &c_b13, &c_b13, &nullsp[nullsp_offset], ldnull, 
	    (ftnlen)4);

/*                                                |VEPS(s)| */
/*     The last NP rows of the product matrix Z x |-------| contain the */
/*                                                |   0   | */
/*     polynomial basis for the right nullspace of the polynomial matrix */
/*     P(s) in condensed form. The multiplication is restricted to the */
/*     nonzero submatrices Vij,k of VEPS, the result is stored in the */
/*     array NULLSP. */

    vc1 = 1;

    i__1 = nblcks;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vr2 = iwork[tail + i__ - 1];

	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {

/*           Multiplication of Z(H+1:N,1:VR2) with V.i,j-1 stored in */
/*           VEPS(1:VR2,VC1:VC1+GAM(I)-1). */

	    dgemm_("No transpose", "No transpose", np, &gam[i__], &vr2, &
		    c_b19, &dwork[jworkz + h__], &n, &dwork[jworkv + (vc1 - 1)
		     * n], &n, &c_b13, &nullsp[vc1 * nullsp_dim1 + 1], ldnull,
		     (ftnlen)12, (ftnlen)12);
	    vc1 += gam[i__];
	    vr2 -= iwork[muk + i__ - j];
/* L40: */
	}

/* L60: */
    }

/*     Transfer of the columns of NULLSP to KER in order to obtain the */
/*     polynomial matrix representation of K(s), the right nullspace */
/*     of P(s). */

    sgamk = 1;

    i__1 = nblcks;
    for (k = 1; k <= i__1; ++k) {
	i__2 = sgamk - 1;
	dlaset_("Full", np, &i__2, &c_b13, &c_b13, &ker[(k * ker_dim2 + 1) * 
		ker_dim1 + 1], ldker1, (ftnlen)4);
	ifir = sgamk;

/*        Copy the appropriate columns of NULLSP into KER(k). */
/*        SGAMK = 1 + SUM(i=1,..,k-1) GAM(i), is the first nontrivial */
/*        column of KER(k), the first SGAMK - 1 columns of KER(k) are */
/*        zero. IFIR denotes the position of the first column in KER(k) */
/*        in the set of columns copied for a value of J. */
/*        VC1 is the first column of NULLSP to be copied. */

	i__2 = nblcks;
	for (j = k; j <= i__2; ++j) {
	    gamj = gam[j];
	    vc1 = iwork[j] + (k - 1) * gamj + 1;
	    dlacpy_("Full", np, &gamj, &nullsp[vc1 * nullsp_dim1 + 1], ldnull,
		     &ker[(ifir + k * ker_dim2) * ker_dim1 + 1], ldker1, (
		    ftnlen)4);
	    ifir += gamj;
/* L80: */
	}

	sgamk += gam[k];
/* L100: */
    }

    return 0;
/* *** Last line of MC03ND *** */
} /* mc03nd_ */

