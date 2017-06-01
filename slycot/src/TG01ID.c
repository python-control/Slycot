/* TG01ID.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tg01id_(char *jobobs, char *compq, char *compz, integer *
	n, integer *m, integer *p, doublereal *a, integer *lda, doublereal *e,
	 integer *lde, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *nobsv, integer *niuobs, integer *nlblck, integer *ctau, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *info, 
	ftnlen jobobs_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, nr, lba, lbe;
    static logical ilq;
    static doublereal dum[1];
    static logical ilz;
    static char jobq[1], jobz[1];
    extern /* Subroutine */ int ma02bd_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), ma02cd_(integer *, integer *, 
	    integer *, doublereal *, integer *), ab07md_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), tg01hx_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), dswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), xerbla_(char *,
	     integer *, ftnlen);
    static logical finobs, infobs;
    static integer icompq, icompz;


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

/*     To compute orthogonal transformation matrices Q and Z which */
/*     reduce the N-th order descriptor system (A-lambda*E,B,C) */
/*     to the form */

/*                ( Ano  * )             ( Eno  * )           ( Bno ) */
/*       Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (     ) , */
/*                ( 0   Ao )             ( 0   Eo )           ( Bo  ) */

/*          C*Z = ( 0   Co ) , */

/*     where the NOBSV-th order descriptor system (Ao-lambda*Eo,Bo,Co) */
/*     is a finite and/or infinite observable. The pencil */
/*     Ano - lambda*Eno is regular of order N-NOBSV and contains the */
/*     unobservable finite and/or infinite eigenvalues of the pencil */
/*     A-lambda*E. */

/*     For JOBOBS = 'O' or 'I', the pencil ( Eo-lambda*Ao ) has full */
/*                                         (      Co      ) */
/*     column rank NOBSV for all finite lambda and is in a staircase form */
/*     with */
/*                     _      _            _      _ */
/*                   ( Ek,k   Ek,k-1   ... Ek,2   Ek,1   ) */
/*                   ( _      _            _      _      ) */
/*       ( Eo ) =    ( Ek-1,k Ek-1,k-1 ... Ek-1,2 Ek-1,1 ) ,  (1) */
/*       ( Co )      (     ...         ... _      _      ) */
/*                   (  0       0      ... E1,2   E1,1   ) */
/*                   (                            _      ) */
/*                   (  0       0      ... 0      E0,1   ) */
/*                     _          _      _ */
/*                   ( Ak,k  ...  Ak,2   Ak,1 ) */
/*                   (       ...  _      _    ) */
/*         Ao      = (   0   ...  A2,2   A2,1 ) ,             (2) */
/*                   (                   _    ) */
/*                   (   0   ...    0    A1,1 ) */
/*           _ */
/*     where Ei-1,i is a CTAU(i-1)-by-CTAU(i) full column rank matrix */
/*                            _ */
/*     (with CTAU(0) = P) and Ai,i is a CTAU(i)-by-CTAU(i) */
/*     upper triangular matrix. */

/*     For JOBOBS = 'F', the pencil ( Ao-lambda*Eo ) has full */
/*                                  (      Co      ) */
/*     column rank NOBSV for all finite lambda and is in a staircase form */
/*     with */
/*                     _      _            _      _ */
/*                   ( Ak,k   Ak,k-1   ... Ak,2   Ak,1   ) */
/*                   ( _      _            _      _      ) */
/*       ( Ao ) =    ( Ak-1,k Ak-1,k-1 ... Ak-1,2 Ak-1,1 ) ,  (3) */
/*       ( Co )      (     ...         ... _      _      ) */
/*                   (  0       0      ... A1,2   A1,1   ) */
/*                   (                            _      ) */
/*                   (  0       0      ... 0      A0,1   ) */
/*                     _          _      _ */
/*                   ( Ek,k  ...  Ek,2   Ek,1 ) */
/*                   (       ...  _      _    ) */
/*         Eo      = (   0   ...  E2,2   E2,1 ) ,             (4) */
/*                   (                   _    ) */
/*                   (   0   ...    0    E1,1 ) */
/*           _ */
/*     where Ai-1,i is a CTAU(i-1)-by-CTAU(i) full column rank matrix */
/*                            _ */
/*     (with CTAU(0) = P) and Ei,i is a CTAU(i)-by-CTAU(i) */
/*     upper triangular matrix. */

/*     For JOBOBS = 'O', the (N-NOBSV)-by-(N-NOBSV) regular pencil */
/*     Ano - lambda*Eno has the form */

/*                         ( Afno - lambda*Efno         *          ) */
/*      Ano - lambda*Eno = (                                       ) , */
/*                         (        0           Aino - lambda*Eino ) */

/*     where: */
/*       1) the NIUOBS-by-NIUOBS regular pencil Aino - lambda*Eino, */
/*          with Aino upper triangular and nonsingular, contains the */
/*          unobservable infinite eigenvalues of A - lambda*E; */
/*       2) the (N-NOBSV-NIUOBS)-by-(N-NOBSV-NIUOBS) regular pencil */
/*          Afno - lambda*Efno, with Efno upper triangular and */
/*          nonsingular, contains the unobservable finite */
/*          eigenvalues of A - lambda*E. */

/*     Note: The significance of the two diagonal blocks can be */
/*           interchanged by calling the routine with the */
/*           arguments A and E interchanged. In this case, */
/*           Aino - lambda*Eino contains the unobservable zero */
/*           eigenvalues of A - lambda*E, while Afno - lambda*Efno */
/*           contains the unobservable nonzero finite and infinite */
/*           eigenvalues of A - lambda*E. */

/*     For JOBOBS = 'F', the pencil Ano - lambda*Eno has the form */

/*        Ano - lambda*Eno = Afno - lambda*Efno , */

/*     where the regular pencil Afno - lambda*Efno, with Efno */
/*     upper triangular and nonsingular, contains the unobservable */
/*     finite eigenvalues of A - lambda*E. */

/*     For JOBOBS = 'I', the pencil Ano - lambda*Eno has the form */

/*        Ano - lambda*Eno = Aino - lambda*Eino , */

/*     where the regular pencil Aino - lambda*Eino, with Aino */
/*     upper triangular and nonsingular, contains the unobservable */
/*     nonzero finite and infinite eigenvalues of A - lambda*E. */

/*     The left and/or right orthogonal transformations Q and Z */
/*     performed to reduce the system matrices can be optionally */
/*     accumulated. */

/*     The reduced order descriptor system (Ao-lambda*Eo,Bo,Co) has */
/*     the same transfer-function matrix as the original system */
/*     (A-lambda*E,B,C). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBOBS   CHARACTER*1 */
/*             = 'O':  separate both finite and infinite unobservable */
/*                     eigenvalues; */
/*             = 'F':  separate only finite unobservable eigenvalues; */
/*             = 'I':  separate only nonzero finite and infinite */
/*                     unobservable eigenvalues. */

/*     COMPQ   CHARACTER*1 */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     orthogonal matrix Q is returned; */
/*             = 'U':  Q must contain an orthogonal matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     orthogonal matrix Z is returned; */
/*             = 'U':  Z must contain an orthogonal matrix Z1 on entry, */
/*                     and the product Z1*Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the descriptor state vector; also the */
/*             order of square matrices A and E, the number of rows of */
/*             matrix B, and the number of columns of matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of descriptor system input vector; also the */
/*             number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of descriptor system output vector; also the */
/*             number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state matrix Q'*A*Z, */

/*                                ( Ano  *  ) */
/*                       Q'*A*Z = (         ) , */
/*                                ( 0    Ao ) */

/*             where Ao is NOBSV-by-NOBSV and Ano is */
/*             (N-NOBSV)-by-(N-NOBSV). */
/*             If JOBOBS = 'F', the matrix ( Ao ) is in the observability */
/*                                         ( Co ) */
/*             staircase form (3). */
/*             If JOBOBS = 'O' or 'I', the submatrix Ao is upper */
/*             triangular. */
/*             If JOBOBS = 'O', the submatrix Ano has the form */

/*                             ( Afno   *  ) */
/*                       Ano = (           ) , */
/*                             (  0   Aino ) */

/*             where the NIUOBS-by-NIUOBS matrix Aino is nonsingular and */
/*             upper triangular. */
/*             If JOBOBS = 'I', Ano is nonsingular and upper triangular. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N descriptor matrix E. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state matrix Q'*E*Z, */

/*                                ( Eno  *  ) */
/*                       Q'*E*Z = (         ) , */
/*                                ( 0    Eo ) */

/*             where Eo is NOBSV-by-NOBSV and Eno is */
/*             (N-NOBSV)-by-(N-NOBSV). */
/*             If JOBOBS = 'O' or 'I', the matrix ( Eo ) is in the */
/*                                                ( Co ) */
/*             observability staircase form (1). */
/*             If JOBOBS = 'F', the submatrix Eo is upper triangular. */
/*             If JOBOBS = 'O', the Eno matrix has the form */

/*                             ( Efno   *  ) */
/*                       Eno = (           ) , */
/*                             (  0   Eino ) */

/*             where the NIUOBS-by-NIUOBS matrix Eino is nilpotent */
/*             and the (N-NOBSV-NIUOBS)-by-(N-NOBSV-NIUOBS) matrix Efno */
/*             is nonsingular and upper triangular. */
/*             If JOBOBS = 'F', Eno is nonsingular and upper triangular. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the N-by-M input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix Q'*B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix */

/*                     C*Z = (  0   Co ) , */

/*             where Co is P-by-NOBSV. */
/*             If JOBOBS = 'O' or 'I', the matrix ( Eo ) is in the */
/*                                                ( Co ) */
/*             observability staircase form (1). */
/*             If JOBOBS = 'F', the matrix ( Ao ) is in the observability */
/*                                         ( Co ) */
/*             staircase form (3). */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If COMPQ = 'N': Q is not referenced. */
/*             If COMPQ = 'I': on entry, Q need not be set; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix Q, */
/*                             where Q' is the product of transformations */
/*                             which are applied to A, E, and B on */
/*                             the left. */
/*             If COMPQ = 'U': on entry, the leading N-by-N part of this */
/*                             array must contain an orthogonal matrix */
/*                             Qc; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix */
/*                             Qc*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,N), if COMPQ = 'U' or 'I'. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If COMPZ = 'N': Z is not referenced. */
/*             If COMPZ = 'I': on entry, Z need not be set; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix Z, */
/*                             which is the product of transformations */
/*                             applied to A, E, and C on the right. */
/*             If COMPZ = 'U': on entry, the leading N-by-N part of this */
/*                             array must contain an orthogonal matrix */
/*                             Zc; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix */
/*                             Zc*Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. */
/*             LDZ >= 1,        if COMPZ = 'N'; */
/*             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'. */

/*     NOBSV   (output) INTEGER */
/*             The order of the reduced matrices Ao and Eo, and the */
/*             number of columns of reduced matrix Co; also the order of */
/*             observable part of the pair (C, A-lambda*E). */

/*     NIUOBS  (output) INTEGER */
/*             For JOBOBS = 'O', the order of the reduced matrices */
/*             Aino and Eino; also the number of unobservable */
/*             infinite eigenvalues of the pencil A - lambda*E. */
/*             For JOBOBS = 'F' or 'I', NIUOBS has no significance */
/*             and is set to zero. */

/*     NLBLCK  (output) INTEGER */
/*             For JOBOBS = 'O' or 'I', the number k, of full column rank */
/*                    _ */
/*             blocks Ei-1,i in the staircase form of the pencil */
/*             (Eo-lambda*Ao) (see (1) and (2)). */
/*             (    Co      ) */
/*             For JOBOBS = 'F', the number k, of full column rank blocks */
/*             _ */
/*             Ai-1,i in the staircase form of the pencil (Ao-lambda*Eo) */
/*                                                        (     Co     ) */
/*             (see (3) and (4)). */

/*     CTAU    (output) INTEGER array, dimension (N) */
/*             CTAU(i), for i = 1, ..., NLBLCK, is the column dimension */
/*                                           _         _ */
/*             of the full column rank block Ei-1,i or Ai-1,i in the */
/*             staircase form (1) or (3) for JOBOBS = 'O' or 'I', or */
/*             for JOBOBS = 'F', respectively. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determinations when */
/*             transforming (A'-lambda*E',C')'. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             reciprocal condition numbers in rank determinations; a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS */
/*             is the machine precision (see LAPACK Library routine */
/*             DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (P) */

/*     DWORK   DOUBLE PRECISION array, dimension MAX(N,2*P) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the dual of the reduction */
/*     algorithms of [1]. */

/*     REFERENCES */

/*     [1] A. Varga */
/*         Computation of Irreducible Generalized State-Space */
/*         Realizations. */
/*         Kybernetika, vol. 26, pp. 89-106, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( N**3 )  floating point operations. */

/*     FURTHER COMMENTS */

/*     If the system matrices A, E and C are badly scaled, it is */
/*     generally recommendable to scale them with the SLICOT routine */
/*     TG01AD, before calling TG01ID. */

/*     CONTRIBUTOR */

/*     C. Oara, University "Politehnica" Bucharest. */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSCF. */

/*     REVISIONS */

/*     July 1999, V. Sima, Research Institute for Informatics, Bucharest. */
/*     May 2003, March 2004, V. Sima. */

/*     KEYWORDS */

/*     Observability, minimal realization, orthogonal canonical form, */
/*     orthogonal transformation. */

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

/*     Decode JOBOBS. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --ctau;
    --iwork;
    --dwork;

    /* Function Body */
    if (lsame_(jobobs, "O", (ftnlen)1, (ftnlen)1)) {
	finobs = TRUE_;
	infobs = TRUE_;
    } else if (lsame_(jobobs, "F", (ftnlen)1, (ftnlen)1)) {
	finobs = TRUE_;
	infobs = FALSE_;
    } else if (lsame_(jobobs, "I", (ftnlen)1, (ftnlen)1)) {
	finobs = FALSE_;
	infobs = TRUE_;
    } else {
	finobs = FALSE_;
	infobs = FALSE_;
    }

/*     Decode COMPQ. */

    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
	ilq = FALSE_;
	icompq = 1;
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
	ilq = TRUE_;
	icompq = 2;
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
	ilq = TRUE_;
	icompq = 3;
    } else {
	icompq = 0;
    }

/*     Decode COMPZ. */

    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
	ilz = FALSE_;
	icompz = 1;
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
	ilz = TRUE_;
	icompz = 2;
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
	ilz = TRUE_;
	icompz = 3;
    } else {
	icompz = 0;
    }

/*     Test the input scalar parameters. */

    *info = 0;
    if (! finobs && ! infobs) {
	*info = -1;
    } else if (icompq <= 0) {
	*info = -2;
    } else if (icompz <= 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*lda < max(1,*n)) {
	*info = -8;
    } else if (*lde < max(1,*n)) {
	*info = -10;
    } else if (*ldb < 1 || *m > 0 && *ldb < *n) {
	*info = -12;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*ldc < max(i__1,*p)) {
	    *info = -14;
	} else if (ilq && *ldq < *n || *ldq < 1) {
	    *info = -16;
	} else if (ilz && *ldz < *n || *ldz < 1) {
	    *info = -18;
	} else if (*tol >= 1.) {
	    *info = -23;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("TG01ID", &i__1, (ftnlen)6);
	return 0;
    }

    *(unsigned char *)jobq = *(unsigned char *)compq;
    *(unsigned char *)jobz = *(unsigned char *)compz;

/*     Build the dual system. */

    ab07md_("Z", n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset]
	    , ldc, dum, &c__1, info, (ftnlen)1);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	dswap_(&i__2, &e[i__ + e_dim1], lde, &e[i__ * e_dim1 + 1], &c__1);
/* L10: */
    }

    if (finobs) {

/*        Perform finite observability form reduction. */

/* Computing MAX */
	i__2 = 0, i__3 = *n - 1;
	i__1 = max(i__2,i__3);
	tg01hx_(jobz, jobq, n, n, p, m, n, &i__1, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], ldc, &z__[
		z_offset], ldz, &q[q_offset], ldq, &nr, nlblck, &ctau[1], tol,
		 &iwork[1], &dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (*nlblck > 1) {
	    lba = ctau[1] + ctau[2] - 1;
	} else if (*nlblck == 1) {
	    lba = ctau[1] - 1;
	} else {
	    lba = 0;
	}
	if (ilq) {
	    *(unsigned char *)jobq = 'U';
	}
	if (ilz) {
	    *(unsigned char *)jobz = 'U';
	}
	lbe = 0;
    } else {
	nr = *n;
/* Computing MAX */
	i__1 = 0, i__2 = *n - 1;
	lba = max(i__1,i__2);
	lbe = lba;
    }

    if (infobs) {

/*        Perform infinite observability form reduction. */

	tg01hx_(jobz, jobq, n, n, p, m, &nr, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, &z__[
		z_offset], ldz, &q[q_offset], ldq, nobsv, nlblck, &ctau[1], 
		tol, &iwork[1], &dwork[1], info, (ftnlen)1, (ftnlen)1);
	if (finobs) {
	    *niuobs = nr - *nobsv;
	} else {
	    *niuobs = 0;
	}
	if (*nlblck > 1) {
	    lbe = ctau[1] + ctau[2] - 1;
	} else if (*nlblck == 1) {
	    lbe = ctau[1] - 1;
	} else {
	    lbe = 0;
	}
	lba = 0;
    } else {
	*nobsv = nr;
	*niuobs = 0;
    }

/*     Compute the pertransposed dual system exploiting matrix shapes. */

/* Computing MAX */
    i__1 = lba, i__2 = *niuobs - 1, i__1 = max(i__1,i__2), i__2 = *n - *nobsv 
	    - *niuobs - 1;
    lba = max(i__1,i__2);
    if (*p == 0 || nr == 0) {
/* Computing MAX */
	i__1 = 0, i__2 = *n - 1;
	lbe = max(i__1,i__2);
    }
/* Computing MAX */
    i__2 = 0, i__3 = *n - 1;
    i__1 = max(i__2,i__3);
    tb01xd_("Z", n, p, m, &lba, &i__1, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, dum, &c__1, info, (ftnlen)1);
/* Computing MAX */
    i__2 = 0, i__3 = *n - 1;
    i__1 = max(i__2,i__3);
    ma02cd_(n, &lbe, &i__1, &e[e_offset], lde);
    if (ilz) {
	ma02bd_("Right", n, n, &z__[z_offset], ldz, (ftnlen)5);
    }
    if (ilq) {
	ma02bd_("Right", n, n, &q[q_offset], ldq, (ftnlen)5);
    }
    return 0;
/* *** Last line of TG01ID *** */
} /* tg01id_ */

