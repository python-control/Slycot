/* AB01OD.f -- translated by f2c (version 20100827).
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
static integer c__1 = 1;

/* Subroutine */ int ab01od_(char *stages, char *jobu, char *jobv, integer *n,
	 integer *m, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	 doublereal *u, integer *ldu, doublereal *v, integer *ldv, integer *
	ncont, integer *indcon, integer *kstair, doublereal *tol, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	stages_len, ftnlen jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, u_dim1, u_offset, v_dim1, 
	    v_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, i0, j0, mm, jini, itau, mcrt, ncrt;
    extern /* Subroutine */ int ab01nd_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dgerqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static logical lstagb, lstgab;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical ljobui, ljobvi;
    static integer ibstep;
    extern /* Subroutine */ int dorgrq_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dormrq_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer wrkopt;


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

/*     To reduce the matrices A and B using (and optionally accumulating) */
/*     state-space and input-space transformations U and V respectively, */
/*     such that the pair of matrices */

/*        Ac = U' * A * U,    Bc = U' * B * V */

/*     are in upper "staircase" form. Specifically, */

/*             [ Acont     *    ]         [ Bcont ] */
/*        Ac = [                ],   Bc = [       ], */
/*             [   0    Auncont ]         [   0   ] */

/*        and */

/*                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ] */
/*                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ] */
/*                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ] */
/*        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ], */
/*                [  .   .     . .    .     .  ]         [ .  ] */
/*                [  .   .       .    .     .  ]         [ .  ] */
/*                [  0   0   . . .  Ap,p-1 App ]         [ 0  ] */

/*     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and */
/*     p is the controllability index of the pair.  The size of the */
/*     block Auncont is equal to the dimension of the uncontrollable */
/*     subspace of the pair (A, B).  The first stage of the reduction, */
/*     the "forward" stage, accomplishes the reduction to the orthogonal */
/*     canonical form (see SLICOT library routine AB01ND). The blocks */
/*     B1, A21, ..., Ap,p-1 are further reduced in a second, "backward" */
/*     stage to upper triangular form using RQ factorization. Each of */
/*     these stages is optional. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STAGES  CHARACTER*1 */
/*             Specifies the reduction stages to be performed as follows: */
/*             = 'F':  Perform the forward stage only; */
/*             = 'B':  Perform the backward stage only; */
/*             = 'A':  Perform both (all) stages. */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the state-space transformations as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is internally initialized to the unit matrix (if */
/*                     STAGES <> 'B'), or updated (if STAGES = 'B'), and */
/*                     the orthogonal transformation matrix U is */
/*                     returned. */

/*     JOBV    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the input-space transformations as follows: */
/*             = 'N':  Do not form V; */
/*             = 'I':  V is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix V is returned. */
/*             JOBV is not referenced if STAGES = 'F'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the */
/*             matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The actual input dimension.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A to be transformed. */
/*             If STAGES = 'B', A should be in the orthogonal canonical */
/*             form, as returned by SLICOT library routine AB01ND. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state transition matrix U' * A * U. */
/*             The leading NCONT-by-NCONT part contains the upper block */
/*             Hessenberg state matrix Acont in Ac, given by U' * A * U, */
/*             of a controllable realization for the original system. */
/*             The elements below the first block-subdiagonal are set to */
/*             zero.  If STAGES <> 'F', the subdiagonal blocks of A are */
/*             triangularized by RQ factorization, and the annihilated */
/*             elements are explicitly zeroed. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B to be transformed. */
/*             If STAGES = 'B', B should be in the orthogonal canonical */
/*             form, as returned by SLICOT library routine AB01ND. */
/*             On exit with STAGES = 'F', the leading N-by-M part of */
/*             this array contains the transformed input matrix U' * B, */
/*             with all elements but the first block set to zero. */
/*             On exit with STAGES <> 'F', the leading N-by-M part of */
/*             this array contains the transformed input matrix */
/*             U' * B * V, with all elements but the first block set to */
/*             zero and the first block in upper triangular form. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             If STAGES <> 'B' or JOBU = 'N', then U need not be set */
/*             on entry. */
/*             If STAGES = 'B' and JOBU = 'I', then, on entry, the */
/*             leading N-by-N part of this array must contain the */
/*             transformation matrix U that reduced the pair to the */
/*             orthogonal canonical form. */
/*             On exit, if JOBU = 'I', the leading N-by-N part of this */
/*             array contains the transformation matrix U that performed */
/*             the specified reduction. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             If JOBU = 'I', LDU >= MAX(1,N);  if JOBU = 'N', LDU >= 1. */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,M) */
/*             If JOBV = 'I', then the leading M-by-M part of this array */
/*             contains the transformation matrix V. */
/*             If STAGES = 'F', or JOBV = 'N', the array V is not */
/*             referenced and can be supplied as a dummy array (i.e. set */
/*             parameter  LDV = 1 and declare this array to be V(1,1) in */
/*             the calling program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. */
/*             If STAGES <> 'F' and JOBV = 'I', LDV >= MAX(1,M); */
/*             if STAGES = 'F' or JOBV = 'N', LDV >= 1. */

/*     NCONT   (input/output) INTEGER */
/*             The order of the controllable state-space representation. */
/*             NCONT is input only if STAGES = 'B'. */

/*     INDCON  (input/output) INTEGER */
/*             The number of stairs in the staircase form (also, the */
/*             controllability index of the controllable part of the */
/*             system representation). */
/*             INDCON is input only if STAGES = 'B'. */

/*     KSTAIR  (input/output) INTEGER array, dimension (N) */
/*             The leading INDCON elements of this array contain the */
/*             dimensions of the stairs, or, also, the orders of the */
/*             diagonal blocks of Acont. */
/*             KSTAIR is input if STAGES = 'B', and output otherwise. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS */
/*             is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */
/*             TOL is not referenced if STAGES = 'B'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */
/*             IWORK is not referenced if STAGES = 'B'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If STAGES <> 'B', LDWORK >= MAX(1, N + MAX(N,3*M)); */
/*             If STAGES =  'B', LDWORK >= MAX(1, M + MAX(N,M)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Staircase reduction of the pencil [B|sI - A] is used. Orthogonal */
/*     transformations U and V are constructed such that */


/*                        |B |sI-A      *  . . .  *      *       | */
/*                        | 1|    11       .      .      .       | */
/*                        |  |  A    sI-A    .    .      .       | */
/*                        |  |   21      22    .  .      .       | */
/*                        |  |        .     .     *      *       | */
/*     [U'BV|sI - U'AU] = |0 |     0    .     .                  | */
/*                        |  |            A     sI-A     *       | */
/*                        |  |             p,p-1    pp           | */
/*                        |  |                                   | */
/*                        |0 |         0          0   sI-A       | */
/*                        |  |                            p+1,p+1| */


/*     where the i-th diagonal block of U'AU has dimension KSTAIR(i), */
/*     for i = 1,...,p. The value of p is returned in INDCON. The last */
/*     block contains the uncontrollable modes of the (A,B)-pair which */
/*     are also the generalized eigenvalues of the above pencil. */

/*     The complete reduction is performed in two stages. The first, */
/*     forward stage accomplishes the reduction to the orthogonal */
/*     canonical form. The second, backward stage consists in further */
/*     reduction to triangular form by applying left and right orthogonal */
/*     transformations. */

/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         The generalized eigenvalue problem in linear system theory. */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 111-129, 1981. */

/*     [2] Miminis, G. and Paige, C. */
/*         An algorithm for pole assignment of time-invariant multi-input */
/*         linear systems. */
/*         Proc. 21st IEEE CDC, Orlando, Florida, 1, pp. 62-67, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O((N + M) x N**2) operations and is */
/*     backward stable (see [1]). */

/*     FURTHER COMMENTS */

/*     If the system matrices A and B are badly scaled, it would be */
/*     useful to scale them with SLICOT routine TB01ID, before calling */
/*     the routine. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB01CD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     January 14, 1997, February 12, 1998, September 22, 2003. */

/*     KEYWORDS */

/*     Controllability, generalized eigenvalue problem, orthogonal */
/*     transformation, staircase form. */

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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --kstair;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobui = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);

    lstagb = lsame_(stages, "B", (ftnlen)1, (ftnlen)1);
    lstgab = lsame_(stages, "A", (ftnlen)1, (ftnlen)1) || lstagb;

    if (lstgab) {
	ljobvi = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
    }

/*     Test the input scalar arguments. */

    if (! lstgab && ! lsame_(stages, "F", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! ljobui && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*m < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (*ldu < 1 || ljobui && *ldu < *n) {
	*info = -11;
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
	i__3 = *n, i__4 = *m * 3;
	i__1 = 1, i__2 = *n + max(i__3,i__4);
/* Computing MAX */
	i__5 = 1, i__6 = *m + max(*n,*m);
	if (! lstagb && *ldwork < max(i__1,i__2) || lstagb && *ldwork < max(
		i__5,i__6)) {
	    *info = -20;
	} else if (lstagb && *ncont > *n) {
	    *info = -14;
	} else if (lstagb && *indcon > *n) {
	    *info = -15;
	} else if (lstgab) {
	    if (! ljobvi && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
		*info = -3;
	    } else if (*ldv < 1 || ljobvi && *ldv < *m) {
		*info = -13;
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("AB01OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*n,*m) == 0) {
	*ncont = 0;
	*indcon = 0;
	if (*n > 0 && ljobui) {
	    dlaset_("F", n, n, &c_b11, &c_b12, &u[u_offset], ldu, (ftnlen)1);
	}
	if (lstgab) {
	    if (*m > 0 && ljobvi) {
		dlaset_("F", m, m, &c_b11, &c_b12, &v[v_offset], ldv, (ftnlen)
			1);
	    }
	}
	dwork[1] = 1.;
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    itau = 1;
    wrkopt = 1;

    if (! lstagb) {

/*        Perform the forward stage computations of the staircase */
/*        algorithm on B and A: reduce the (A, B) pair to orthogonal */
/*        canonical form. */

/*        Workspace: N + MAX(N,3*M). */

	jwork = *n + 1;
	i__1 = *ldwork - jwork + 1;
	ab01nd_(jobu, n, m, &a[a_offset], lda, &b[b_offset], ldb, ncont, 
		indcon, &kstair[1], &u[u_offset], ldu, &dwork[itau], tol, &
		iwork[1], &dwork[jwork], &i__1, info, (ftnlen)1);

	wrkopt = (integer) dwork[jwork] + jwork - 1;
    }

/*     Exit if no further reduction to triangularize B1 and subdiagonal */
/*     blocks of A is required, or if the order of the controllable part */
/*     is 0. */

    if (! lstgab) {
	dwork[1] = (doublereal) wrkopt;
	return 0;
    } else if (*ncont == 0 || *indcon == 0) {
	if (ljobvi) {
	    dlaset_("F", m, m, &c_b11, &c_b12, &v[v_offset], ldv, (ftnlen)1);
	}
	dwork[1] = (doublereal) wrkopt;
	return 0;
    }

/*     Now perform the backward steps except the last one. */

    mcrt = kstair[*indcon];
    i0 = *ncont - mcrt + 1;
    jwork = *m + 1;

    for (ibstep = *indcon; ibstep >= 2; --ibstep) {
	ncrt = kstair[ibstep - 1];
	j0 = i0 - ncrt;
	mm = min(ncrt,mcrt);

/*        Compute the RQ factorization of the current subdiagonal block */
/*        of A, Ai,i-1 = R*Q (where i is IBSTEP), of dimension */
/*        MCRT-by-NCRT, starting in position (I0,J0). */
/*        The matrix Q' should postmultiply U, if required. */
/*        Workspace: need   M + MCRT; */
/*                   prefer M + MCRT*NB. */

	i__1 = *ldwork - jwork + 1;
	dgerqf_(&mcrt, &ncrt, &a[i0 + j0 * a_dim1], lda, &dwork[itau], &dwork[
		jwork], &i__1, info);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Set JINI to the first column number in A where the current */
/*        transformation Q is to be applied, taking the block Hessenberg */
/*        form into account. */

	if (ibstep > 2) {
	    jini = j0 - kstair[ibstep - 2];
	} else {
	    jini = 1;

/*           Premultiply the first block row (B1) of B by Q. */
/*           Workspace: need   2*M; */
/*                      prefer M + M*NB. */

	    i__1 = *ldwork - jwork + 1;
	    dormrq_("Left", "No transpose", &ncrt, m, &mm, &a[i0 + j0 * 
		    a_dim1], lda, &dwork[itau], &b[b_offset], ldb, &dwork[
		    jwork], &i__1, info, (ftnlen)4, (ftnlen)12);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__1,i__2);
	}

/*        Premultiply the appropriate block row of A by Q. */
/*        Workspace: need   M + N; */
/*                   prefer M + N*NB. */

	i__1 = *n - jini + 1;
	i__2 = *ldwork - jwork + 1;
	dormrq_("Left", "No transpose", &ncrt, &i__1, &mm, &a[i0 + j0 * 
		a_dim1], lda, &dwork[itau], &a[j0 + jini * a_dim1], lda, &
		dwork[jwork], &i__2, info, (ftnlen)4, (ftnlen)12);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

/*        Postmultiply the appropriate block column of A by Q'. */
/*        Workspace: need   M +  I0-1; */
/*                   prefer M + (I0-1)*NB. */

	i__1 = i0 - 1;
	i__2 = *ldwork - jwork + 1;
	dormrq_("Right", "Transpose", &i__1, &ncrt, &mm, &a[i0 + j0 * a_dim1],
		 lda, &dwork[itau], &a[j0 * a_dim1 + 1], lda, &dwork[jwork], &
		i__2, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);

	if (ljobui) {

/*           Update U, postmultiplying it by Q'. */
/*           Workspace: need   M + N; */
/*                      prefer M + N*NB. */

	    i__1 = *ldwork - jwork + 1;
	    dormrq_("Right", "Transpose", n, &ncrt, &mm, &a[i0 + j0 * a_dim1],
		     lda, &dwork[itau], &u[j0 * u_dim1 + 1], ldu, &dwork[
		    jwork], &i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	    wrkopt = max(i__1,i__2);
	}

/*        Zero the subdiagonal elements of the current subdiagonal block */
/*        of A. */

	i__1 = ncrt - mcrt;
	dlaset_("F", &mcrt, &i__1, &c_b11, &c_b11, &a[i0 + j0 * a_dim1], lda, 
		(ftnlen)1);
	if (i0 < *n) {
	    i__1 = mcrt - 1;
	    i__2 = mcrt - 1;
	    dlaset_("L", &i__1, &i__2, &c_b11, &c_b11, &a[i0 + 1 + (i0 - mcrt)
		     * a_dim1], lda, (ftnlen)1);
	}

	mcrt = ncrt;
	i0 = j0;

/* L10: */
    }

/*     Now perform the last backward step on B, V = Qb'. */

/*     Compute the RQ factorization of the first block of B, B1 = R*Qb. */
/*     Workspace: need   M + MCRT; */
/*                prefer M + MCRT*NB. */

    i__1 = *ldwork - jwork + 1;
    dgerqf_(&mcrt, m, &b[b_offset], ldb, &dwork[itau], &dwork[jwork], &i__1, 
	    info);
/* Computing MAX */
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
    wrkopt = max(i__1,i__2);

    if (ljobvi) {

/*        Accumulate the input-space transformations V. */
/*        Workspace: need 2*M;  prefer M + M*NB. */

	i__1 = *m - mcrt;
	dlacpy_("F", &mcrt, &i__1, &b[b_offset], ldb, &v[*m - mcrt + 1 + 
		v_dim1], ldv, (ftnlen)1);
	if (mcrt > 1) {
	    i__1 = mcrt - 1;
	    i__2 = mcrt - 1;
	    dlacpy_("L", &i__1, &i__2, &b[(*m - mcrt + 1) * b_dim1 + 2], ldb, 
		    &v[*m - mcrt + 2 + (*m - mcrt + 1) * v_dim1], ldv, (
		    ftnlen)1);
	}
	i__1 = *ldwork - jwork + 1;
	dorgrq_(m, m, &mcrt, &v[v_offset], ldv, &dwork[itau], &dwork[jwork], &
		i__1, info);

	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = i__ - 1;
	    dswap_(&i__2, &v[i__ + v_dim1], ldv, &v[i__ * v_dim1 + 1], &c__1);
/* L20: */
	}

/* Computing MAX */
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
	wrkopt = max(i__1,i__2);
    }

/*     Zero the subdiagonal elements of the submatrix B1. */

    i__1 = *m - mcrt;
    dlaset_("F", &mcrt, &i__1, &c_b11, &c_b11, &b[b_offset], ldb, (ftnlen)1);
    if (mcrt > 1) {
	i__1 = mcrt - 1;
	i__2 = mcrt - 1;
	dlaset_("L", &i__1, &i__2, &c_b11, &c_b11, &b[(*m - mcrt + 1) * 
		b_dim1 + 2], ldb, (ftnlen)1);
    }

/*     Set optimal workspace dimension. */

    dwork[1] = (doublereal) wrkopt;
    return 0;
/* *** Last line of AB01OD *** */
} /* ab01od_ */

