/* NF01BP.f -- translated by f2c (version 20100827).
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
static doublereal c_b53 = 1.;

/* Subroutine */ int nf01bp_(char *cond, integer *n, integer *ipar, integer *
	lipar, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag,
	 doublereal *qtb, doublereal *delta, doublereal *par, integer *ranks, 
	doublereal *x, doublereal *rx, doublereal *tol, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, n2, bn;
    static doublereal fp;
    static integer jw, st, bsm, bsn, lds;
    static doublereal sum, parc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ibsn, rank;
    static doublereal parl;
    static logical sing;
    static integer iter;
    static doublereal temp, paru;
    static integer nths;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static logical badrk;
    extern /* Subroutine */ int nf01bq_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd;
    extern /* Subroutine */ int nf01br_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    md03by_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, ftnlen);
    static char condl[1];
    static logical ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwarf;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal dmino;
    static logical ucond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal gnorm;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal dxnorm;


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

/*     To determine a value for the Levenberg-Marquardt parameter PAR */
/*     such that if x solves the system */

/*           J*x = b ,     sqrt(PAR)*D*x = 0 , */

/*     in the least squares sense, where J is an m-by-n matrix, D is an */
/*     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if */
/*     DELTA is a positive number, DXNORM is the Euclidean norm of D*x, */
/*     then either PAR is zero and */

/*           ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*     The matrix J is the current Jacobian matrix of a nonlinear least */
/*     squares problem, provided in a compressed form by SLICOT Library */
/*     routine NF01BD. It is assumed that a block QR factorization, with */
/*     column pivoting, of J is available, that is, J*P = Q*R, where P is */
/*     a permutation matrix, Q has orthogonal columns, and R is an upper */
/*     triangular matrix with diagonal elements of nonincreasing */
/*     magnitude for each block, as returned by SLICOT Library */
/*     routine NF01BS. The routine NF01BP needs the upper triangle of R */
/*     in compressed form, the permutation matrix P, and the first */
/*     n components of Q'*b (' denotes the transpose). On output, */
/*     NF01BP also provides a compressed representation of an upper */
/*     triangular matrix S, such that */

/*           P'*(J'*J + PAR*D*D)*P = S'*S . */

/*     Matrix S is used in the solution process. The matrix R has the */
/*     following structure */

/*         /   R_1    0    ..   0   |   L_1   \ */
/*         |    0    R_2   ..   0   |   L_2   | */
/*         |    :     :    ..   :   |    :    | , */
/*         |    0     0    ..  R_l  |   L_l   | */
/*         \    0     0    ..   0   |  R_l+1  / */

/*     where the submatrices R_k, k = 1:l, have the same order BSN, */
/*     and R_k, k = 1:l+1, are square and upper triangular. This matrix */
/*     is stored in the compressed form */

/*              /   R_1  |   L_1   \ */
/*              |   R_2  |   L_2   | */
/*       Rc =   |    :   |    :    | , */
/*              |   R_l  |   L_l   | */
/*              \    X   |  R_l+1  / */

/*     where the submatrix X is irrelevant. The matrix S has the same */
/*     structure as R, and its diagonal blocks are denoted by S_k, */
/*     k = 1:l+1. */

/*     If l <= 1, then the full upper triangle of the matrix R is stored. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the diagonal blocks R_k */
/*             and S_k of the matrices R and S should be estimated, */
/*             as follows: */
/*             = 'E' :  use incremental condition estimation for each */
/*                      diagonal block of R_k and S_k to find its */
/*                      numerical rank; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R_k and S_k for zero values; */
/*             = 'U' :  use the ranks already stored in RANKS (for R). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N = BN*BSN + ST >= 0. */
/*             (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix R, as follows: */
/*             IPAR(1) must contain ST, the number of columns of the */
/*                     submatrices L_k and the order of R_l+1.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, l, in the */
/*                     block diagonal part of R.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     R_k, k = 1:l.  BSM >= 0. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks R_k, k = 1:l.  BSN >= 0. */
/*             BSM is not used by this routine, but assumed equal to BSN. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             On entry, the leading N-by-NC part of this array must */
/*             contain the (compressed) representation (Rc) of the upper */
/*             triangular matrix R. If BN > 1, the submatrix X in Rc is */
/*             not referenced. The zero strict lower triangles of R_k, */
/*             k = 1:l+1, need not be set. If BN <= 1 or BSN = 0, then */
/*             the full upper triangle of R must be stored. */
/*             On exit, the full upper triangles of R_k, k = 1:l+1, and */
/*             L_k, k = 1:l, are unaltered, and the strict lower */
/*             triangles of R_k, k = 1:l+1, contain the corresponding */
/*             strict upper triangles (transposed) of the upper */
/*             triangular matrix S. */
/*             If BN <= 1 or BSN = 0, then the transpose of the strict */
/*             upper triangle of S is stored in the strict lower triangle */
/*             of R. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             An upper bound on the Euclidean norm of D*x.  DELTA > 0. */

/*     PAR     (input/output) DOUBLE PRECISION */
/*             On entry, PAR must contain an initial estimate of the */
/*             Levenberg-Marquardt parameter.  PAR >= 0. */
/*             On exit, it contains the final estimate of this parameter. */

/*     RANKS   (input or output) INTEGER array, dimension (r), where */
/*             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1; */
/*             r = BN,      if ST = 0 and BSN > 0; */
/*             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 ); */
/*             r = 0,       if ST = 0 and BSN = 0. */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical ranks of the submatrices R_k, k = 1:l(+1). */
/*             On exit, if N > 0, this array contains the numerical ranks */
/*             of the submatrices S_k, k = 1:l(+1). */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system J*x = b, sqrt(PAR)*D*x = 0. */

/*     RX      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the matrix-vector product -R*P'*x. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the submatrices R_k and S_k. If the user sets */
/*             TOL > 0, then the given value of TOL is used as a lower */
/*             bound for the reciprocal condition number;  a (sub)matrix */
/*             whose estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S. */
/*             If BN > 1 and BSN > 0, the elements N+1 : N+ST*(N-ST) */
/*             contain the submatrix (S(1:N-ST,N-ST+1:N))' of the */
/*             matrix S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*N,              if BN <= 1 or  BSN = 0 and */
/*                                                        COND <> 'E'; */
/*             LDWORK >= 4*N,              if BN <= 1 or  BSN = 0 and */
/*                                                        COND =  'E'; */
/*             LDWORK >= ST*(N-ST) + 2*N,  if BN >  1 and BSN > 0 and */
/*                                                        COND <> 'E'; */
/*             LDWORK >= ST*(N-ST) + 2*N + 2*MAX(BSN,ST), */
/*                                         if BN >  1 and BSN > 0 and */
/*                                                        COND =  'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm computes the Gauss-Newton direction. An approximate */
/*     basic least squares solution is found if the Jacobian is rank */
/*     deficient. The computations exploit the special structure and */
/*     storage scheme of the matrix R. If one or more of the submatrices */
/*     R_k or S_k, k = 1:l+1, is singular, then the computed result is */
/*     not the basic least squares solution for the whole problem, but a */
/*     concatenation of (least squares) solutions of the individual */
/*     subproblems involving R_k or S_k, k = 1:l+1 (with adapted right */
/*     hand sides). */

/*     If the Gauss-Newton direction is not acceptable, then an iterative */
/*     algorithm obtains improved lower and upper bounds for the */
/*     Levenberg-Marquardt parameter PAR. Only a few iterations are */
/*     generally needed for convergence of the algorithm. If, however, */
/*     the limit of ITMAX = 10 iterations is reached, then the output PAR */
/*     will contain the best value obtained so far. If the Gauss-Newton */
/*     step is acceptable, it is stored in x, and PAR is set to zero, */
/*     hence S = R. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N*(BSN+ST)) operations and is backward */
/*     stable, if R is nonsingular. */

/*     FURTHER COMMENTS */

/*     This routine is a structure-exploiting, LAPACK-based modification */
/*     of LMPAR from the MINPACK package [1], and with optional condition */
/*     estimation. The option COND = 'U' is useful when dealing with */
/*     several right-hand side vectors, but RANKS array should be reset. */
/*     If COND = 'E', but the matrix S is guaranteed to be nonsingular */
/*     and well conditioned relative to TOL, i.e., rank(R) = N, and */
/*     min(DIAG) > 0, then its condition is not estimated. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Feb. 2004. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --ipar;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --ipvt;
    --diag;
    --qtb;
    --ranks;
    --x;
    --rx;
    --dwork;

    /* Function Body */
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
    ucond = lsame_(cond, "U", (ftnlen)1, (ftnlen)1);
    *info = 0;
    n2 = *n << 1;
    if (! (econd || ncond || ucond)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lipar < 4) {
	*info = -4;
    } else if (*ldr < max(1,*n)) {
	*info = -6;
    } else if (*delta <= 0.) {
	*info = -10;
    } else if (*par < 0.) {
	*info = -11;
    } else {
	st = ipar[1];
	bn = ipar[2];
	bsm = ipar[3];
	bsn = ipar[4];
	nths = bn * bsn;
/* Computing MIN */
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
	if (min(i__1,bsn) < 0) {
	    *info = -3;
	} else if (*n != nths + st) {
	    *info = -2;
	} else {
	    if (*n > 0) {
		dmino = diag[1];
	    }
	    sing = FALSE_;

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (diag[j] < dmino) {
		    dmino = diag[j];
		}
		sing = sing || diag[j] == 0.;
/* L10: */
	    }

	    if (sing) {
		*info = -8;
	    } else if (ucond) {
		badrk = FALSE_;
		if (bn <= 1 || bsn == 0) {
		    if (*n > 0) {
			badrk = ranks[1] < 0 || ranks[1] > *n;
		    }
		} else {
		    rank = 0;

		    i__1 = bn;
		    for (k = 1; k <= i__1; ++k) {
			badrk = badrk || ranks[k] < 0 || ranks[k] > bsn;
			rank += ranks[k];
/* L20: */
		    }

		    if (st > 0) {
			badrk = badrk || ranks[bn + 1] < 0 || ranks[bn + 1] > 
				st;
			rank += ranks[bn + 1];
		    }
		}
		if (badrk) {
		    *info = -12;
		}
	    } else {
		jw = n2;
		if (bn <= 1 || bsn == 0) {
		    if (econd) {
			jw = *n << 2;
		    }
		} else {
		    jw = st * nths + jw;
		    if (econd) {
			jw = (max(bsn,st) << 1) + jw;
		    }
		}
		if (*ldwork < jw) {
		    *info = -17;
		}
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BP", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	*par = 0.;
	return 0;
    }

    if (bn <= 1 || bsn == 0) {

/*        Special case: R is just an upper triangular matrix. */
/*        Workspace: 4*N, if COND =  'E'; */
/*                   2*N, if COND <> 'E'. */

	md03by_(cond, n, &r__[r_offset], ldr, &ipvt[1], &diag[1], &qtb[1], 
		delta, par, &ranks[1], &x[1], &rx[1], tol, &dwork[1], ldwork, 
		info, (ftnlen)1);
	return 0;
    }

/*     General case: l > 1 and BSN > 0. */
/*     DWARF is the smallest positive magnitude. */

    dwarf = dlamch_("Underflow", (ftnlen)9);

/*     Compute and store in x the Gauss-Newton direction. If the */
/*     Jacobian is rank-deficient, obtain a least squares solution. */
/*     The array RX is used as workspace. */
/*     Workspace: 2*MAX(BSN,ST), if COND =  'E'; */
/*                0,             if COND <> 'E'. */

    dcopy_(n, &qtb[1], &c__1, &rx[1], &c__1);
    nf01br_(cond, "Upper", "No transpose", n, &ipar[1], lipar, &r__[r_offset],
	     ldr, &dwork[1], &dwork[1], &c__1, &rx[1], &ranks[1], tol, &dwork[
	    1], ldwork, info, (ftnlen)1, (ftnlen)5, (ftnlen)12);

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = rx[j];
/* L30: */
    }

/*     Initialize the iteration counter. */
/*     Evaluate the function at the origin, and test */
/*     for acceptance of the Gauss-Newton direction. */

    iter = 0;

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dwork[j] = diag[j] * x[j];
/* L40: */
    }

    dxnorm = dnrm2_(n, &dwork[1], &c__1);
    fp = dxnorm - *delta;
    if (fp > *delta * .1) {

/*        Set an appropriate option for estimating the condition of */
/*        the matrix S. */

	lds = max(1,st);
	jw = n2 + st * nths;
	if (ucond) {
	    if (*ldwork >= jw + (max(bsn,st) << 1)) {
		*(unsigned char *)condl = 'E';
		toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
	    } else {
		*(unsigned char *)condl = 'N';
		toldef = *tol;
	    }
	} else {
	    rank = 0;

	    i__1 = bn;
	    for (k = 1; k <= i__1; ++k) {
		rank += ranks[k];
/* L50: */
	    }

	    if (st > 0) {
		rank += ranks[bn + 1];
	    }
	    *(unsigned char *)condl = *(unsigned char *)cond;
	    toldef = *tol;
	}

/*        If the Jacobian is not rank deficient, the Newton */
/*        step provides a lower bound, PARL, for the zero of */
/*        the function. Otherwise set this bound to zero. */

	if (rank == *n) {

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		rx[j] = diag[l] * (dwork[l] / dxnorm);
/* L60: */
	    }

	    nf01br_("Use ranks", "Upper", "Transpose", n, &ipar[1], lipar, &
		    r__[r_offset], ldr, &dwork[1], &dwork[1], &c__1, &rx[1], &
		    ranks[1], tol, &dwork[1], ldwork, info, (ftnlen)9, (
		    ftnlen)5, (ftnlen)9);
	    temp = dnrm2_(n, &rx[1], &c__1);
	    parl = fp / *delta / temp / temp;

/*           For efficiency, use CONDL = 'U', if possible. */

	    if (! lsame_(condl, "U", (ftnlen)1, (ftnlen)1) && dmino > 0.) {
		*(unsigned char *)condl = 'U';
	    }
	} else {
	    parl = 0.;
	}

	ibsn = 0;
	k = 1;

/*        Calculate an upper bound, PARU, for the zero of the function. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ++ibsn;
	    if (j < nths) {
		sum = ddot_(&ibsn, &r__[k + ibsn * r_dim1], &c__1, &qtb[k], &
			c__1);
		if (ibsn == bsn) {
		    ibsn = 0;
		    k += bsn;
		}
	    } else if (j == nths) {
		sum = ddot_(&ibsn, &r__[k + ibsn * r_dim1], &c__1, &qtb[k], &
			c__1);
	    } else {
		sum = ddot_(&j, &r__[ibsn * r_dim1 + 1], &c__1, &qtb[1], &
			c__1);
	    }
	    l = ipvt[j];
	    rx[j] = sum / diag[l];
/* L70: */
	}

	gnorm = dnrm2_(n, &rx[1], &c__1);
	paru = gnorm / *delta;
	if (paru == 0.) {
	    paru = dwarf / min(*delta,.1) / .001;
	}

/*        If the input PAR lies outside of the interval (PARL,PARU), */
/*        set PAR to the closer endpoint. */

	*par = max(*par,parl);
	*par = min(*par,paru);
	if (*par == 0.) {
	    *par = gnorm / dxnorm;
	}

/*        Beginning of an iteration. */

L80:
	++iter;

/*           Evaluate the function at the current value of PAR. */

	if (*par == 0.) {
/* Computing MAX */
	    d__1 = dwarf, d__2 = paru * .001;
	    *par = max(d__1,d__2);
	}
	temp = sqrt(*par);

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    rx[j] = temp * diag[j];
/* L90: */
	}

/*           Solve the system J*x = b , sqrt(PAR)*D*x = 0 , in a least */
/*           square sense. */
/*           The first N elements of DWORK contain the diagonal elements */
/*           of the upper triangular matrix S, and the next N elements */
/*           contain the the vector z, so that x = P*z (see NF01BQ). */
/*           The vector z is not preserved, to reduce the workspace. */
/*           The elements 2*N+1 : 2*N+ST*(N-ST) contain the */
/*           submatrix (S(1:N-ST,N-ST+1:N))' of the matrix S. */
/*           Workspace: ST*(N-ST) + 2*N,                 if CONDL <> 'E'; */
/*                      ST*(N-ST) + 2*N + 2*MAX(BSN,ST), if CONDL =  'E'. */

	nf01bq_(condl, n, &ipar[1], lipar, &r__[r_offset], ldr, &ipvt[1], &rx[
		1], &qtb[1], &ranks[1], &x[1], &toldef, &dwork[1], ldwork, 
		info, (ftnlen)1);

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dwork[*n + j] = diag[j] * x[j];
/* L100: */
	}

	dxnorm = dnrm2_(n, &dwork[*n + 1], &c__1);
	temp = fp;
	fp = dxnorm - *delta;

/*           If the function is small enough, accept the current value */
/*           of PAR. Also test for the exceptional cases where PARL */
/*           is zero or the number of iterations has reached ITMAX. */

	if (abs(fp) > *delta * .1 && (parl != 0. || fp > temp || temp >= 0.) 
		&& iter < 10) {

/*              Compute the Newton correction. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		l = ipvt[j];
		rx[j] = diag[l] * (dwork[*n + l] / dxnorm);
/* L110: */
	    }

	    i__1 = *ldwork - jw;
	    nf01br_("Use ranks", "Lower", "Transpose", n, &ipar[1], lipar, &
		    r__[r_offset], ldr, &dwork[1], &dwork[n2 + 1], &lds, &rx[
		    1], &ranks[1], tol, &dwork[jw], &i__1, info, (ftnlen)9, (
		    ftnlen)5, (ftnlen)9);
	    temp = dnrm2_(n, &rx[1], &c__1);
	    parc = fp / *delta / temp / temp;

/*              Depending on the sign of the function, update PARL */
/*              or PARU. */

	    if (fp > 0.) {
		parl = max(parl,*par);
	    } else if (fp < 0.) {
		paru = min(paru,*par);
	    }

/*              Compute an improved estimate for PAR. */

/* Computing MAX */
	    d__1 = parl, d__2 = *par + parc;
	    *par = max(d__1,d__2);

/*              End of an iteration. */

	    goto L80;
	}
    }

/*     Compute -R*P'*x = -R*z. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	rx[j] = -x[l];
/* L120: */
    }

    i__1 = nths;
    i__2 = bsn;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dtrmv_("Upper", "NoTranspose", "NonUnit", &bsn, &r__[i__ + r_dim1], 
		ldr, &rx[i__], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
/* L130: */
    }

    if (st > 0) {
	dgemv_("NoTranspose", &nths, &st, &c_b53, &r__[(bsn + 1) * r_dim1 + 1]
		, ldr, &rx[nths + 1], &c__1, &c_b53, &rx[1], &c__1, (ftnlen)
		11);
	dtrmv_("Upper", "NoTranspose", "NonUnit", &st, &r__[nths + 1 + (bsn + 
		1) * r_dim1], ldr, &rx[nths + 1], &c__1, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
    }

/*     Termination. If PAR = 0, set S. */

    if (iter == 0) {
	*par = 0.;
	i__ = 1;

	i__2 = bn;
	for (k = 1; k <= i__2; ++k) {

	    i__1 = bsn;
	    for (j = 1; j <= i__1; ++j) {
		dwork[i__] = r__[i__ + j * r_dim1];
		i__3 = bsn - j + 1;
		dcopy_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
		++i__;
/* L140: */
	    }

/* L150: */
	}

	if (st > 0) {

	    i__2 = bsn + st;
	    for (j = bsn + 1; j <= i__2; ++j) {
		dcopy_(&nths, &r__[j * r_dim1 + 1], &c__1, &dwork[*n + j - 
			bsn], &st);
		dwork[i__] = r__[i__ + j * r_dim1];
		i__1 = bsn + st - j + 1;
		dcopy_(&i__1, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
		++i__;
/* L160: */
	    }

	}
    } else {

	i__2 = *n + st * nths;
	for (k = *n + 1; k <= i__2; ++k) {
	    dwork[k] = dwork[k + *n];
/* L170: */
	}

    }

    return 0;

/* *** Last line of NF01BP *** */
} /* nf01bp_ */

