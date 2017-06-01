/* MB04VD.f -- translated by f2c (version 20100827).
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
static doublereal c_b14 = 1.;

/* Subroutine */ int mb04vd_(char *mode, char *jobq, char *jobz, integer *m, 
	integer *n, integer *ranke, doublereal *a, integer *lda, doublereal *
	e, integer *lde, doublereal *q, integer *ldq, doublereal *z__, 
	integer *ldz, integer *istair, integer *nblcks, integer *nblcki, 
	integer *imuk, integer *inuk, integer *imuk0, integer *mnei, 
	doublereal *tol, integer *iwork, integer *info, ftnlen mode_len, 
	ftnlen jobq_len, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, jk, nca, nra, ifica, ifira, ranka;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04tt_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, doublereal *, integer *), 
	    mb04ty_(logical *, logical *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), mb04vx_(logical *, logical *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    static doublereal toler, dwork[1];
    static logical first;
    static integer ismuk, isnuk;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical lmodeb;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobqi, lmodes, lmodet, ljobzi, updatq, firsti, updatz;


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

/*     To compute orthogonal transformations Q and Z such that the */
/*     transformed pencil Q'(sE-A)Z is in upper block triangular form, */
/*     where E is an M-by-N matrix in column echelon form (see SLICOT */
/*     Library routine MB04UD) and A is an M-by-N matrix. */

/*     If MODE = 'B', then the matrices A and E are transformed into the */
/*     following generalized Schur form by unitary transformations Q1 */
/*     and Z1 : */

/*                      | sE(eps,inf)-A(eps,inf) |      X     | */
/*        Q1'(sE-A)Z1 = |------------------------|------------|.   (1) */
/*                      |            O           | sE(r)-A(r) | */

/*     The pencil sE(eps,inf)-A(eps,inf) is in staircase form, and it */
/*     contains all Kronecker column indices and infinite elementary */
/*     divisors of the pencil sE-A. The pencil sE(r)-A(r) contains all */
/*     Kronecker row indices and elementary divisors of sE-A. */
/*     Note: X is a pencil. */

/*     If MODE = 'T', then the submatrices having full row and column */
/*     rank in the pencil sE(eps,inf)-A(eps,inf) in (1) are */
/*     triangularized by applying unitary transformations Q2 and Z2 to */
/*     Q1'*(sE-A)*Z1. */

/*     If MODE = 'S', then the pencil sE(eps,inf)-A(eps,inf) in (1) is */
/*     separated into sE(eps)-A(eps) and sE(inf)-A(inf) by applying */
/*     unitary transformations Q3 and Z3 to Q2'*Q1'*(sE-A)*Z1*Z2. */

/*     This gives */

/*                | sE(eps)-A(eps) |        X       |      X     | */
/*                |----------------|----------------|------------| */
/*                |        O       | sE(inf)-A(inf) |      X     | */
/*     Q'(sE-A)Z =|=================================|============| (2) */
/*                |                                 |            | */
/*                |                O                | sE(r)-A(r) | */

/*     where Q = Q1*Q2*Q3 and Z = Z1*Z2*Z3. */
/*     Note: the pencil sE(r)-A(r) is not reduced further. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     MODE    CHARACTER*1 */
/*             Specifies the desired structure of the transformed */
/*             pencil Q'(sE-A)Z to be computed as follows: */
/*             = 'B':  Basic reduction given by (1); */
/*             = 'T':  Further reduction of (1) to triangular form; */
/*             = 'S':  Further separation of sE(eps,inf)-A(eps,inf) */
/*                     in (1) into the two pencils in (2). */

/*     JOBQ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the orthogonal row transformations, as follows: */
/*             = 'N':  Do not form Q; */
/*             = 'I':  Q is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Q is returned; */
/*             = 'U':  The given matrix Q is updated by the orthogonal */
/*                     row transformations used in the reduction. */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal column transformations, as */
/*             follows: */
/*             = 'N':  Do not form Z; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned; */
/*             = 'U':  The given matrix Z is updated by the orthogonal */
/*                     transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in the matrices A, E and the order of */
/*             the matrix Q.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in the matrices A, E and the order */
/*             of the matrix Z.  N >= 0. */

/*     RANKE   (input) INTEGER */
/*             The rank of the matrix E in column echelon form. */
/*             RANKE >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix to be row compressed. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix that has been row compressed while keeping */
/*             matrix E in column echelon form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix in column echelon form to be */
/*             transformed equivalent to matrix A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix that has been transformed equivalent to matrix */
/*             A. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,M). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             On entry, if JOBQ = 'U', then the leading M-by-M part of */
/*             this array must contain a given matrix Q (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading M-by-M part of this array contains the product of */
/*             the input matrix Q and the row transformation matrix used */
/*             to transform the rows of matrices A and E. */
/*             On exit, if JOBQ = 'I', then the leading M-by-M part of */
/*             this array contains the matrix of accumulated orthogonal */
/*             row transformations performed. */
/*             If JOBQ = 'N', the array Q is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDQ = 1 and */
/*             declare this array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. If JOBQ = 'U' or */
/*             JOBQ = 'I', LDQ >= MAX(1,M); if JOBQ = 'N', LDQ >= 1. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*) */
/*             On entry, if JOBZ = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix Z (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix Z and the column transformation matrix */
/*             used to transform the columns of matrices A and E. */
/*             On exit, if JOBZ = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated orthogonal */
/*             column transformations performed. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'U' or */
/*             JOBZ = 'I', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     ISTAIR  (input/output) INTEGER array, dimension (M) */
/*             On entry, this array must contain information on the */
/*             column echelon form of the unitary transformed matrix E. */
/*             Specifically, ISTAIR(i) must be set to +j if the first */
/*             non-zero element E(i,j) is a corner point and -j */
/*             otherwise, for i = 1,2,...,M. */
/*             On exit, this array contains no useful information. */

/*     NBLCKS  (output) INTEGER */
/*             The number of submatrices having full row rank greater */
/*             than or equal to 0 detected in matrix A in the pencil */
/*             sE(x)-A(x), */
/*                where  x = eps,inf  if MODE = 'B' or 'T', */
/*                or     x = eps      if MODE = 'S'. */

/*     NBLCKI  (output) INTEGER */
/*             If MODE = 'S', the number of diagonal submatrices in the */
/*             pencil sE(inf)-A(inf). If MODE = 'B' or 'T' then */
/*             NBLCKI = 0. */

/*     IMUK    (output) INTEGER array, dimension (MAX(N,M+1)) */
/*             The leading NBLCKS elements of this array contain the */
/*             column dimensions mu(1),...,mu(NBLCKS) of the submatrices */
/*             having full column rank in the pencil sE(x)-A(x), */
/*                where  x = eps,inf  if MODE = 'B' or 'T', */
/*                or     x = eps      if MODE = 'S'. */

/*     INUK    (output) INTEGER array, dimension (MAX(N,M+1)) */
/*             The leading NBLCKS elements of this array contain the */
/*             row dimensions nu(1),...,nu(NBLCKS) of the submatrices */
/*             having full row rank in the pencil sE(x)-A(x), */
/*                where  x = eps,inf  if MODE = 'B' or 'T', */
/*                or     x = eps      if MODE = 'S'. */

/*     IMUK0   (output) INTEGER array, dimension (limuk0), */
/*             where limuk0 = N if MODE = 'S' and 1, otherwise. */
/*             If MODE = 'S', then the leading NBLCKI elements of this */
/*             array contain the dimensions mu0(1),...,mu0(NBLCKI) */
/*             of the square diagonal submatrices in the pencil */
/*             sE(inf)-A(inf). */
/*             Otherwise, IMUK0 is not referenced and can be supplied */
/*             as a dummy array. */

/*     MNEI    (output) INTEGER array, dimension (3) */
/*             If MODE = 'B' or 'T' then */
/*             MNEI(1) contains the row dimension of */
/*                     sE(eps,inf)-A(eps,inf); */
/*             MNEI(2) contains the column dimension of */
/*                     sE(eps,inf)-A(eps,inf); */
/*             MNEI(3) = 0. */
/*             If MODE = 'S', then */
/*             MNEI(1) contains the row    dimension of sE(eps)-A(eps); */
/*             MNEI(2) contains the column dimension of sE(eps)-A(eps); */
/*             MNEI(3) contains the order of the regular pencil */
/*                     sE(inf)-A(inf). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance below which matrix elements are considered */
/*             to be zero. If the user sets TOL to be less than (or */
/*             equal to) zero then the tolerance is taken as */
/*             EPS * MAX( ABS(A(I,J)), ABS(E(I,J)) ), where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH), */
/*             I = 1,2,...,M and J = 1,2,...,N. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*             > 0:  if incorrect rank decisions were revealed during the */
/*                   triangularization phase. This failure is not likely */
/*                   to occur. The possible values are: */
/*             = 1:  if incorrect dimensions of a full column rank */
/*                   submatrix; */
/*             = 2:  if incorrect dimensions of a full row rank */
/*                   submatrix. */

/*     METHOD */

/*     Let sE - A be an arbitrary pencil. Prior to calling the routine, */
/*     this pencil must be transformed into a pencil with E in column */
/*     echelon form. This may be accomplished by calling the SLICOT */
/*     Library routine MB04UD. Depending on the value of MODE, */
/*     submatrices of A and E are then reduced to one of the forms */
/*     described above. Further details can be found in [1]. */

/*     REFERENCES */

/*     [1] Beelen, Th. and Van Dooren, P. */
/*         An improved algorithm for the computation of Kronecker's */
/*         canonical form of a singular pencil. */
/*         Linear Algebra and Applications, 105, pp. 9-65, 1988. */

/*     NUMERICAL ASPECTS */

/*     It is shown in [1] that the algorithm is numerically backward */
/*     stable. The operations count is proportional to (MAX(M,N))**3. */

/*     FURTHER COMMENTS */

/*     The difference mu(k)-nu(k), for k = 1,2,...,NBLCKS, is the number */
/*     of elementary Kronecker blocks of size k x (k+1). */

/*     If MODE = 'B' or 'T' on entry, then the difference nu(k)-mu(k+1), */
/*     for k = 1,2,...,NBLCKS, is the number of infinite elementary */
/*     divisors of degree k (with mu(NBLCKS+1) = 0). */

/*     If MODE = 'S' on entry, then the difference mu0(k)-mu0(k+1), */
/*     for k = 1,2,...,NBLCKI, is the number of infinite elementary */
/*     divisors of degree k (with mu0(NBLCKI+1) = 0). */
/*     In the pencil sE(r)-A(r), the pencils sE(f)-A(f) and */
/*     sE(eta)-A(eta) can be separated by pertransposing the pencil */
/*     sE(r)-A(r) and calling the routine with MODE set to 'B'. The */
/*     result has got to be pertransposed again. (For more details see */
/*     [1]). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Based on Release 3.0 routine MB04TD modified by A. Varga, */
/*     German Aerospace Research Establishment, Oberpfaffenhofen, */
/*     Germany, Nov. 1997, as follows: */
/*     1) NBLCKI is added; */
/*     2) the significance of IMUK0 and MNEI is changed; */
/*     3) INUK0 is removed. */

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
/*     .. Local Arrays .. */
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
    --imuk;
    --inuk;
    --imuk0;
    --mnei;
    --iwork;

    /* Function Body */
    *info = 0;
    lmodeb = lsame_(mode, "B", (ftnlen)1, (ftnlen)1);
    lmodet = lsame_(mode, "T", (ftnlen)1, (ftnlen)1);
    lmodes = lsame_(mode, "S", (ftnlen)1, (ftnlen)1);
    ljobqi = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
    updatq = ljobqi || lsame_(jobq, "U", (ftnlen)1, (ftnlen)1);
    ljobzi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
    updatz = ljobzi || lsame_(jobz, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

    if (! lmodeb && ! lmodet && ! lmodes) {
	*info = -1;
    } else if (! updatq && ! lsame_(jobq, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (! updatz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*ranke < 0) {
	*info = -6;
    } else if (*lda < max(1,*m)) {
	*info = -8;
    } else if (*lde < max(1,*m)) {
	*info = -10;
    } else if (! updatq && *ldq < 1 || updatq && *ldq < max(1,*m)) {
	*info = -12;
    } else if (! updatz && *ldz < 1 || updatz && *ldz < max(1,*n)) {
	*info = -14;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MB04VD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Initialize Q and Z to the identity matrices, if needed. */

    if (ljobqi) {
	dlaset_("Full", m, m, &c_b13, &c_b14, &q[q_offset], ldq, (ftnlen)4);
    }
    if (ljobzi) {
	dlaset_("Full", n, n, &c_b13, &c_b14, &z__[z_offset], ldz, (ftnlen)4);
    }

/*     Quick return if possible. */

    *nblcks = 0;
    *nblcki = 0;

    if (*n == 0) {
	mnei[1] = 0;
	mnei[2] = 0;
	mnei[3] = 0;
	return 0;
    }

    if (*m == 0) {
	*nblcks = *n;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    imuk[i__] = 1;
	    inuk[i__] = 0;
/* L10: */
	}
	mnei[1] = 0;
	mnei[2] = *n;
	mnei[3] = 0;
	return 0;
    }

    toler = *tol;
    if (toler <= 0.) {
/* Computing MAX */
	d__1 = dlange_("M", m, n, &a[a_offset], lda, dwork, (ftnlen)1), d__2 =
		 dlange_("M", m, n, &e[e_offset], lde, dwork, (ftnlen)1);
	toler = dlamch_("Epsilon", (ftnlen)7) * max(d__1,d__2);
    }

/*     A(k) is the submatrix in A that will be row compressed. */

/*     ISMUK = sum(i=1,..,k) MU(i), ISNUK = sum(i=1,...,k) NU(i), */
/*     IFIRA, IFICA: first row and first column index of A(k) in A. */
/*     NRA, NCA: number of rows and columns in A(k). */

    ifira = 1;
    ifica = 1;
    nra = *m;
    nca = *n - *ranke;
    isnuk = 0;
    ismuk = 0;
    k = 0;

/*     Initialization of the arrays INUK and IMUK. */

    i__1 = *m + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	inuk[i__] = -1;
/* L20: */
    }

/*     Note: it is necessary that array INUK has DIMENSION M+1 since it */
/*           is possible that M = 1 and NBLCKS = 2. */
/*           Example sE-A = (0 0 s -1). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	imuk[i__] = -1;
/* L40: */
    }

/*     Compress the rows of A while keeping E in column echelon form. */

/*     REPEAT */

L60:
    ++k;
    mb04tt_(&updatq, &updatz, m, n, &ifira, &ifica, &nca, &a[a_offset], lda, &
	    e[e_offset], lde, &q[q_offset], ldq, &z__[z_offset], ldz, &istair[
	    1], &ranka, &toler, &iwork[1]);
    imuk[k] = nca;
    ismuk += nca;

    inuk[k] = ranka;
    isnuk += ranka;
    ++(*nblcks);

/*        If the rank of A(k) is nra then A has full row rank; */
/*        JK = the first column index (in A) after the right most column */
/*        of matrix A(k+1). (In case A(k+1) is empty, then JK = N+1.) */

    ifira = isnuk + 1;
    ifica = ismuk + 1;
    if (ifira > *m) {
	jk = *n + 1;
    } else {
	jk = (i__1 = istair[ifira], abs(i__1));
    }
    nra = *m - isnuk;
    nca = jk - 1 - ismuk;

/*        If NCA > 0 then there can be done some more row compression */
/*        of matrix A while keeping matrix E in column echelon form. */

    if (nca > 0) {
	goto L60;
    }
/*     UNTIL NCA <= 0 */

/*     Matrix E(k+1) has full column rank since NCA = 0. */
/*     Reduce A and E by ignoring all rows and columns corresponding */
/*     to E(k+1). Ignoring these columns in E changes the ranks of the */
/*     submatrices E(i), (i=1,...,k-1). */

    mnei[1] = isnuk;
    mnei[2] = ismuk;
    mnei[3] = 0;

    if (lmodeb) {
	return 0;
    }

/*     Triangularization of the submatrices in A and E. */

    mb04ty_(&updatq, &updatz, m, n, nblcks, &inuk[1], &imuk[1], &a[a_offset], 
	    lda, &e[e_offset], lde, &q[q_offset], ldq, &z__[z_offset], ldz, 
	    info);

    if (*info > 0 || lmodet) {
	return 0;
    }

/*     Save the row dimensions of the diagonal submatrices in pencil */
/*     sE(eps,inf)-A(eps,inf). */

    i__1 = *nblcks;
    for (i__ = 1; i__ <= i__1; ++i__) {
	imuk0[i__] = inuk[i__];
/* L80: */
    }

/*     Reduction to square submatrices E(k)'s in E. */

    mb04vx_(&updatq, &updatz, m, n, nblcks, &inuk[1], &imuk[1], &a[a_offset], 
	    lda, &e[e_offset], lde, &q[q_offset], ldq, &z__[z_offset], ldz, &
	    mnei[1]);

/*     Determine the dimensions of the inf diagonal submatrices and */
/*     update block numbers if necessary. */

    first = TRUE_;
    firsti = TRUE_;
    *nblcki = *nblcks;
    k = *nblcks;

    for (i__ = k; i__ >= 1; --i__) {
	imuk0[i__] -= inuk[i__];
	if (firsti && imuk0[i__] == 0) {
	    --(*nblcki);
	} else {
	    firsti = FALSE_;
	}
	if (first && imuk[i__] == 0) {
	    --(*nblcks);
	} else {
	    first = FALSE_;
	}
/* L100: */
    }

    return 0;
/* *** Last line of MB04VD *** */
} /* mb04vd_ */

