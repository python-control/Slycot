/* IB01PD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static doublereal c_b43 = .66666666666666663;
static doublereal c_b63 = 1.;
static doublereal c_b66 = 0.;
static doublereal c_b173 = -1.;

/* Subroutine */ int ib01pd_(char *meth, char *job, char *jobcv, integer *
	nobr, integer *n, integer *m, integer *l, integer *nsmpl, doublereal *
	r__, integer *ldr, doublereal *a, integer *lda, doublereal *c__, 
	integer *ldc, doublereal *b, integer *ldb, doublereal *d__, integer *
	ldd, doublereal *q, integer *ldq, doublereal *ry, integer *ldry, 
	doublereal *s, integer *lds, doublereal *o, integer *ldo, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen meth_len, ftnlen job_len, ftnlen 
	jobcv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, o_dim1, o_offset, q_dim1, q_offset, r_dim1, r_offset, 
	    ry_dim1, ry_offset, s_dim1, s_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, n2, id, nn, iu, nr, ix, nr2, nr3, nr4, iaw, ldw;
    static doublereal eps;
    static integer npl, isv, iun2, igal;
    static char fact[1], jobp[1];
    static integer ncol, rank, ierr, itau;
    static doublereal sval[3], toll, rnrm;
    static integer nrow;
    static logical n4sid;
    static integer itau1, itau2, ldun2;
    static doublereal toll1;
    static integer nr4mn, nr4pl;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02ed_(char *, integer *, doublereal *, integer *, ftnlen), 
	    mb03od_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), mb02ud_(char *, char *, 
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static integer rank11;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ib01px_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), ib01py_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer rankm;
    extern /* Subroutine */ int mb02qy_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer lnobr, mnobr;
    static logical shift, withb;
    static integer ldunn;
    static logical withc;
    static char jobpy[1];
    static logical fullr, moesp;
    static integer ihous;
    static logical withd;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dtrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal rcond1, rcond2, rcond3, rcond4;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lmmnob;
    static logical withal;
    static integer lmnobr, lnobrn, mnobrn, iwarnl;
    static doublereal thresh;
    static integer lmmnol;
    static logical withco;
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dormqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal svlmax;
    extern /* Subroutine */ int dtrtrs_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


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

/*     To estimate the matrices A, C, B, and D of a linear time-invariant */
/*     (LTI) state space model, using the singular value decomposition */
/*     information provided by other routines. Optionally, the system and */
/*     noise covariance matrices, needed for the Kalman gain, are also */
/*     determined. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     JOB     CHARACTER*1 */
/*             Specifies which matrices should be computed, as follows: */
/*             = 'A':  compute all system matrices, A, B, C, and D; */
/*             = 'C':  compute the matrices A and C only; */
/*             = 'B':  compute the matrix B only; */
/*             = 'D':  compute the matrices B and D only. */

/*     JOBCV   CHARACTER*1 */
/*             Specifies whether or not the covariance matrices are to */
/*             be computed, as follows: */
/*             = 'C':  the covariance matrices should be computed; */
/*             = 'N':  the covariance matrices should not be computed. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             Hankel matrices processed by other routines.  NOBR > 1. */

/*     N       (input) INTEGER */
/*             The order of the system.  NOBR > N > 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMPL   (input) INTEGER */
/*             If JOBCV = 'C', the total number of samples used for */
/*             calculating the covariance matrices. */
/*             NSMPL >= 2*(M+L)*NOBR. */
/*             This parameter is not meaningful if  JOBCV = 'N'. */

/*     R       (input/workspace) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On entry, the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  part */
/*             of this array must contain the relevant data for the MOESP */
/*             or N4SID algorithms, as constructed by SLICOT Library */
/*             routines IB01AD or IB01ND. Let  R_ij,  i,j = 1:4,  be the */
/*             ij submatrix of  R  (denoted  S  in IB01AD and IB01ND), */
/*             partitioned by  M*NOBR,  L*NOBR,  M*NOBR,  and  L*NOBR */
/*             rows and columns. The submatrix  R_22  contains the matrix */
/*             of left singular vectors used. Also needed, for */
/*             METH = 'N'  or  JOBCV = 'C',  are the submatrices  R_11, */
/*             R_14 : R_44,  and, for  METH = 'M'  and  JOB <> 'C',  the */
/*             submatrices  R_31  and  R_12,  containing the processed */
/*             matrices  R_1c  and  R_2c,  respectively, as returned by */
/*             SLICOT Library routines IB01AD or IB01ND. */
/*             Moreover, if  METH = 'N'  and  JOB = 'A' or 'C',  the */
/*             block-row  R_41 : R_43  must contain the transpose of the */
/*             block-column  R_14 : R_34  as returned by SLICOT Library */
/*             routines IB01AD or IB01ND. */
/*             The remaining part of  R  is used as workspace. */
/*             On exit, part of this array is overwritten. Specifically, */
/*             if  METH = 'M',  R_22  and  R_31  are overwritten if */
/*                 JOB = 'B' or 'D',  and  R_12,  R_22,  R_14 : R_34, */
/*                 and possibly  R_11  are overwritten if  JOBCV = 'C'; */
/*             if  METH = 'N',  all needed submatrices are overwritten. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     A       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             On entry, if  METH = 'N'  and  JOB = 'B' or 'D',  the */
/*             leading N-by-N part of this array must contain the system */
/*             state matrix. */
/*             If  METH = 'M'  or  (METH = 'N'  and JOB = 'A' or 'C'), */
/*             this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  the */
/*             leading N-by-N part of this array contains the system */
/*             state matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= N,  if  JOB = 'A' or 'C',  or  METH = 'N'  and */
/*                            JOB = 'B' or 'D'; */
/*             LDA >= 1,  otherwise. */

/*     C       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDC,N) */
/*             On entry, if  METH = 'N'  and  JOB = 'B' or 'D',  the */
/*             leading L-by-N part of this array must contain the system */
/*             output matrix. */
/*             If  METH = 'M'  or  (METH = 'N'  and JOB = 'A' or 'C'), */
/*             this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  or */
/*             INFO = 3  (or  INFO >= 0,  for  METH = 'M'),  the leading */
/*             L-by-N part of this array contains the system output */
/*             matrix. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= L,  if  JOB = 'A' or 'C',  or  METH = 'N'  and */
/*                            JOB = 'B' or 'D'; */
/*             LDC >= 1,  otherwise. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If  M > 0,  JOB = 'A', 'B', or 'D'  and  INFO = 0,  the */
/*             leading N-by-M part of this array contains the system */
/*             input matrix. If  M = 0  or  JOB = 'C',  this array is */
/*             not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= N,  if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             LDB >= 1,  if M = 0 or  JOB = 'C'. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If  M > 0,  JOB = 'A' or 'D'  and  INFO = 0,  the leading */
/*             L-by-M part of this array contains the system input-output */
/*             matrix. If  M = 0  or  JOB = 'C' or 'B',  this array is */
/*             not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if M > 0 and JOB = 'A' or 'D'; */
/*             LDD >= 1,  if M = 0 or  JOB = 'C' or 'B'. */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If JOBCV = 'C', the leading N-by-N part of this array */
/*             contains the positive semidefinite state covariance matrix */
/*             to be used as state weighting matrix when computing the */
/*             Kalman gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= N,  if JOBCV = 'C'; */
/*             LDQ >= 1,  if JOBCV = 'N'. */

/*     RY      (output) DOUBLE PRECISION array, dimension (LDRY,L) */
/*             If JOBCV = 'C', the leading L-by-L part of this array */
/*             contains the positive (semi)definite output covariance */
/*             matrix to be used as output weighting matrix when */
/*             computing the Kalman gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDRY    INTEGER */
/*             The leading dimension of the array RY. */
/*             LDRY >= L,  if JOBCV = 'C'; */
/*             LDRY >= 1,  if JOBCV = 'N'. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,L) */
/*             If JOBCV = 'C', the leading N-by-L part of this array */
/*             contains the state-output cross-covariance matrix to be */
/*             used as cross-weighting matrix when computing the Kalman */
/*             gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= N,  if JOBCV = 'C'; */
/*             LDS >= 1,  if JOBCV = 'N'. */

/*     O       (output) DOUBLE PRECISION array, dimension ( LDO,N ) */
/*             If  METH = 'M'  and  JOBCV = 'C',  or  METH = 'N', */
/*             the leading  L*NOBR-by-N  part of this array contains */
/*             the estimated extended observability matrix, i.e., the */
/*             first  N  columns of the relevant singular vectors. */
/*             If  METH = 'M'  and  JOBCV = 'N',  this array is not */
/*             referenced. */

/*     LDO     INTEGER */
/*             The leading dimension of the array  O. */
/*             LDO >= L*NOBR,  if  JOBCV = 'C'  or  METH = 'N'; */
/*             LDO >= 1,       otherwise. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  an m-by-n matrix whose estimated */
/*             condition number is less than  1/TOL  is considered to */
/*             be of full rank.  If the user sets  TOL <= 0,  then an */
/*             implicitly computed, default tolerance, defined by */
/*             TOLDEF = m*n*EPS,  is used instead, where  EPS  is the */
/*             relative machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = N,                   if METH = 'M' and M = 0 */
/*                                        or JOB = 'C' and JOBCV = 'N'; */
/*             LIWORK = M*NOBR+N,            if METH = 'M', JOB = 'C', */
/*                                           and JOBCV = 'C'; */
/*             LIWORK = max(L*NOBR,M*NOBR),  if METH = 'M', JOB <> 'C', */
/*                                           and JOBCV = 'N'; */
/*             LIWORK = max(L*NOBR,M*NOBR+N),  if METH = 'M', JOB <> 'C', */
/*                                             and JOBCV = 'C'; */
/*             LIWORK = max(M*NOBR+N,M*(N+L)), if METH = 'N'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and  DWORK(2),  DWORK(3),  DWORK(4),  and */
/*             DWORK(5)  contain the reciprocal condition numbers of the */
/*             triangular factors of the matrices, defined in the code, */
/*             GaL  (GaL = Un(1:(s-1)*L,1:n)),  R_1c  (if  METH = 'M'), */
/*             M  (if  JOBCV = 'C'  or  METH = 'N'),  and  Q  or  T  (see */
/*             SLICOT Library routines IB01PY or IB01PX),  respectively. */
/*             If  METH = 'N',  DWORK(3)  is set to one without any */
/*             calculations. Similarly, if  METH = 'M'  and  JOBCV = 'N', */
/*             DWORK(4)  is set to one. If  M = 0  or  JOB = 'C', */
/*             DWORK(3)  and  DWORK(5)  are set to one. */
/*             On exit, if  INFO = -30,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( LDW1,LDW2 ), where, if METH = 'M', */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N ), */
/*                     if JOB = 'C' or JOB = 'A' and M = 0; */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+N*N+7*N, */
/*                          (L*NOBR-L)*N+N+6*M*NOBR, (L*NOBR-L)*N+N+ */
/*                          max( L+M*NOBR, L*NOBR + */
/*                                         max( 3*L*NOBR+1, M ) ) ) */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             LDW2 >= 0,                                 if JOBCV = 'N'; */
/*             LDW2 >= max( (L*NOBR-L)*N+Aw+2*N+max(5*N,(2*M+L)*NOBR+L), */
/*                          4*(M*NOBR+N)+1, M*NOBR+2*N+L ), */
/*                                                        if JOBCV = 'C', */
/*             where Aw = N+N*N, if M = 0 or JOB = 'C'; */
/*                   Aw = 0,     otherwise; */
/*             and, if METH = 'N', */
/*             LDW1 >= max( (L*NOBR-L)*N+2*N+(2*M+L)*NOBR+L, */
/*                          2*(L*NOBR-L)*N+N*N+8*N, N+4*(M*NOBR+N)+1, */
/*                          M*NOBR+3*N+L ); */
/*             LDW2 >= 0, if M = 0 or JOB = 'C'; */
/*             LDW2 >= M*NOBR*(N+L)*(M*(N+L)+1)+ */
/*                     max( (N+L)**2, 4*M*(N+L)+1 ), */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  a least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix; */
/*             = 5:  the computed covariance matrices are too small. */
/*                   The problem seems to be a deterministic one. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge; */
/*             = 3:  a singular upper triangular matrix was found. */

/*     METHOD */

/*     In the MOESP approach, the matrices  A  and  C  are first */
/*     computed from an estimated extended observability matrix [1], */
/*     and then, the matrices  B  and  D  are obtained by solving an */
/*     extended linear system in a least squares sense. */
/*     In the N4SID approach, besides the estimated extended */
/*     observability matrix, the solutions of two least squares problems */
/*     are used to build another least squares problem, whose solution */
/*     is needed to compute the system matrices  A,  C,  B,  and  D.  The */
/*     solutions of the two least squares problems are also optionally */
/*     used by both approaches to find the covariance matrices. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error state- */
/*         space model identification class of algorithms. */
/*         Int. J. Control, 56, pp. 1187-1210, 1992. */

/*     [2] Van Overschee, P., and De Moor, B. */
/*         N4SID: Two Subspace Algorithms for the Identification */
/*         of Combined Deterministic-Stochastic Systems. */
/*         Automatica, Vol.30, No.1, pp. 75-93, 1994. */

/*     [3] Van Overschee, P. */
/*         Subspace Identification : Theory - Implementation - */
/*         Applications. */
/*         Ph. D. Thesis, Department of Electrical Engineering, */
/*         Katholieke Universiteit Leuven, Belgium, Feb. 1995. */

/*     [4] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     FURTHER COMMENTS */

/*     In some applications, it is useful to compute the system matrices */
/*     using two calls to this routine, the first one with  JOB = 'C', */
/*     and the second one with  JOB = 'B' or 'D'.  This is slightly less */
/*     efficient than using a single call with  JOB = 'A',  because some */
/*     calculations are repeated. If  METH = 'N',  all the calculations */
/*     at the first call are performed again at the second call; */
/*     moreover, it is required to save the needed submatrices of  R */
/*     before the first call and restore them before the second call. */
/*     If the covariance matrices are desired,  JOBCV  should be set */
/*     to  'C'  at the second call. If  B  and  D  are both needed, they */
/*     should be computed at once. */
/*     It is possible to compute the matrices A and C using the MOESP */
/*     algorithm (METH = 'M'), and the matrices B and D using the N4SID */
/*     algorithm (METH = 'N'). This combination could be slightly more */
/*     efficient than N4SID algorithm alone and it could be more accurate */
/*     than MOESP algorithm. No saving/restoring is needed in such a */
/*     combination, provided  JOBCV  is set to  'N'  at the first call. */
/*     Recommended usage:  either one call with  JOB = 'A',  or */
/*        first  call with  METH = 'M',  JOB = 'C',  JOBCV = 'N', */
/*        second call with  METH = 'M',  JOB = 'D',  JOBCV = 'C',  or */
/*        first  call with  METH = 'M',  JOB = 'C',  JOBCV = 'N', */
/*        second call with  METH = 'N',  JOB = 'D',  JOBCV = 'C'. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 1999. */

/*     REVISIONS */

/*     March 2000, Feb. 2001, Sep. 2001, March 2005. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Array .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    ry_dim1 = *ldry;
    ry_offset = 1 + ry_dim1;
    ry -= ry_offset;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    o_dim1 = *ldo;
    o_offset = 1 + o_dim1;
    o -= o_offset;
    --iwork;
    --dwork;

    /* Function Body */
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
    withal = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
    withc = lsame_(job, "C", (ftnlen)1, (ftnlen)1) || withal;
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1) || withal;
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
    withco = lsame_(jobcv, "C", (ftnlen)1, (ftnlen)1);
    mnobr = *m * *nobr;
    lnobr = *l * *nobr;
    lmnobr = lnobr + mnobr;
    lmmnob = lnobr + (mnobr << 1);
    mnobrn = mnobr + *n;
    lnobrn = lnobr - *n;
    ldun2 = lnobr - *l;
    ldunn = ldun2 * *n;
    lmmnol = lmmnob + *l;
    nr = lmnobr + lmnobr;
    npl = *n + *l;
    n2 = *n + *n;
    nn = *n * *n;
    minwrk = 1;
    *iwarn = 0;
    *info = 0;

/*     Check the scalar input parameters. */

    if (! (moesp || n4sid)) {
	*info = -1;
    } else if (! (withb || withc)) {
	*info = -2;
    } else if (! (withco || lsame_(jobcv, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -3;
    } else if (*nobr <= 1) {
	*info = -4;
    } else if (*n <= 0 || *n >= *nobr) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*l <= 0) {
	*info = -7;
    } else if (withco && *nsmpl < nr) {
	*info = -8;
    } else if (*ldr < nr) {
	*info = -10;
    } else if (*lda < 1 || (withc || withb && n4sid) && *lda < *n) {
	*info = -12;
    } else if (*ldc < 1 || (withc || withb && n4sid) && *ldc < *l) {
	*info = -14;
    } else if (*ldb < 1 || withb && *ldb < *n && *m > 0) {
	*info = -16;
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
	*info = -18;
    } else if (*ldq < 1 || withco && *ldq < *n) {
	*info = -20;
    } else if (*ldry < 1 || withco && *ldry < *l) {
	*info = -22;
    } else if (*lds < 1 || withco && *lds < *n) {
	*info = -24;
    } else if (*ldo < 1 || (withco || n4sid) && *ldo < lnobr) {
	*info = -26;
    } else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance. */
/*         NB refers to the optimal block size for the immediately */
/*         following subroutine, as returned by ILAENV.) */

	iaw = 0;
	minwrk = ldunn + (*n << 2);
	maxwrk = ldunn + *n + *n * ilaenv_(&c__1, "DGEQRF", " ", &ldun2, n, &
		c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
	if (moesp) {
	    id = 0;
	    if (withc) {
/* Computing MAX */
		i__1 = minwrk, i__2 = (ldunn << 1) + n2, i__1 = max(i__1,i__2)
			, i__2 = ldunn + nn + *n * 7;
		minwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = (ldunn << 1) + *n + *n * ilaenv_(&c__1, 
			"DORMQR", "LT", &ldun2, n, n, &c_n1, (ftnlen)6, (
			ftnlen)2);
		maxwrk = max(i__1,i__2);
	    }
	} else {
	    id = *n;
	}

	if (*m > 0 && withb || n4sid) {
/* Computing MAX */
	    i__1 = minwrk, i__2 = (ldunn << 1) + nn + id + *n * 7;
	    minwrk = max(i__1,i__2);
	    if (moesp) {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
		i__5 = lnobr * 3 + 1;
		i__3 = *l + mnobr, i__4 = lnobr + max(i__5,*m);
		i__1 = minwrk, i__2 = ldunn + *n + mnobr * 6, i__1 = max(i__1,
			i__2), i__2 = ldunn + *n + max(i__3,i__4);
		minwrk = max(i__1,i__2);
	    }
	} else {
	    if (moesp) {
		iaw = *n + nn;
	    }
	}

	if (n4sid || withco) {
/* Computing MAX */
/* Computing MAX */
	    i__3 = *n * 5;
	    i__1 = minwrk, i__2 = ldunn + iaw + n2 + max(i__3,lmmnol), i__1 = 
		    max(i__1,i__2), i__2 = id + (mnobrn << 2) + 1, i__1 = max(
		    i__1,i__2), i__2 = id + mnobrn + npl;
	    minwrk = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
	    i__3 = *n * ilaenv_(&c__1, "DGEQRF", " ", &lnobr, n, &c_n1, &c_n1,
		     (ftnlen)6, (ftnlen)1), i__4 = lmmnob * ilaenv_(&c__1, 
		    "DORMQR", "LT", &lnobr, &lmmnob, n, &c_n1, (ftnlen)6, (
		    ftnlen)2), i__3 = max(i__3,i__4), i__4 = lmmnol * ilaenv_(
		    &c__1, "DORMQR", "LT", &ldun2, &lmmnol, n, &c_n1, (ftnlen)
		    6, (ftnlen)2);
	    i__1 = maxwrk, i__2 = ldunn + iaw + n2 + max(i__3,i__4), i__1 = 
		    max(i__1,i__2), i__2 = id + *n + *n * ilaenv_(&c__1, 
		    "DGEQRF", " ", &lmnobr, n, &c_n1, &c_n1, (ftnlen)6, (
		    ftnlen)1), i__1 = max(i__1,i__2), i__2 = id + *n + npl * 
		    ilaenv_(&c__1, "DORMQR", "LT", &lmnobr, &npl, n, &c_n1, (
		    ftnlen)6, (ftnlen)2);
	    maxwrk = max(i__1,i__2);
	    if (n4sid && (*m > 0 && withb)) {
/* Computing MAX */
/* Computing MAX */
/* Computing 2nd power */
		i__5 = npl;
		i__3 = i__5 * i__5, i__4 = (*m << 2) * npl + 1;
		i__1 = minwrk, i__2 = mnobr * npl * (*m * npl + 1) + max(i__3,
			i__4);
		minwrk = max(i__1,i__2);
	    }
	}
	maxwrk = max(minwrk,maxwrk);

	if (*ldwork < minwrk) {
	    *info = -30;
	    dwork[1] = (doublereal) minwrk;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01PD", &i__1, (ftnlen)6);
	return 0;
    }

    nr2 = mnobr + 1;
    nr3 = lmnobr + 1;
    nr4 = lmmnob + 1;

/*     Set the precision parameters. A threshold value  EPS**(2/3)  is */
/*     used for deciding to use pivoting or not, where  EPS  is the */
/*     relative machine precision (see LAPACK Library routine DLAMCH). */

    eps = dlamch_("Precision", (ftnlen)9);
    thresh = pow_dd(&eps, &c_b43);
    svlmax = 0.;
    rcond4 = 1.;

/*     Let  Un  be the matrix of left singular vectors (stored in  R_22). */
/*     Copy  un1 = GaL = Un(1:(s-1)*L,1:n)  in the workspace. */

    igal = 1;
    dlacpy_("Full", &ldun2, n, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[igal], &
	    ldun2, (ftnlen)4);

/*     Factor un1 = Q1*[r1'  0]' (' means transposition). */
/*     Workspace: need   L*(NOBR-1)*N+2*N, */
/*                prefer L*(NOBR-1)*N+N+N*NB. */

    itau1 = igal + ldunn;
    jwork = itau1 + *n;
    ldw = jwork;
    i__1 = *ldwork - jwork + 1;
    dgeqrf_(&ldun2, n, &dwork[igal], &ldun2, &dwork[itau1], &dwork[jwork], &
	    i__1, &ierr);

/*     Compute the reciprocal of the condition number of r1. */
/*     Workspace: need L*(NOBR-1)*N+4*N. */

    dtrcon_("1-norm", "Upper", "NonUnit", n, &dwork[igal], &ldun2, &rcond1, &
	    dwork[jwork], &iwork[1], info, (ftnlen)6, (ftnlen)5, (ftnlen)7);

    toll1 = *tol;
    if (toll1 <= 0.) {
	toll1 = nn * eps;
    }

    if (*m > 0 && withb || n4sid) {
	*(unsigned char *)jobp = 'P';
	if (withal) {
	    *(unsigned char *)jobpy = 'D';
	} else {
	    *(unsigned char *)jobpy = *(unsigned char *)job;
	}
    } else {
	*(unsigned char *)jobp = 'N';
    }

    if (moesp) {
	ncol = 0;
	iun2 = jwork;
	if (withc) {

/*           Set  C = Un(1:L,1:n)  and then compute the system matrix A. */

/*           Set  un2 = Un(L+1:L*s,1:n)  in  DWORK(IUN2). */
/*           Workspace: need   2*L*(NOBR-1)*N+N. */

	    dlacpy_("Full", l, n, &r__[nr2 + nr2 * r_dim1], ldr, &c__[
		    c_offset], ldc, (ftnlen)4);
	    dlacpy_("Full", &ldun2, n, &r__[nr2 + *l + nr2 * r_dim1], ldr, &
		    dwork[iun2], &ldun2, (ftnlen)4);

/*           Note that un1 has already been factored as */
/*           un1 = Q1*[r1'  0]'  and usually (generically, assuming */
/*           observability) has full column rank. */
/*           Update  un2 <-- Q1'*un2  in  DWORK(IUN2)  and save its */
/*           first  n  rows in  A. */
/*           Workspace: need   2*L*(NOBR-1)*N+2*N; */
/*                      prefer 2*L*(NOBR-1)*N+N+N*NB. */

	    jwork = iun2 + ldunn;
	    i__1 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &ldun2, n, n, &dwork[igal], &ldun2, &
		    dwork[itau1], &dwork[iun2], &ldun2, &dwork[jwork], &i__1, 
		    &ierr, (ftnlen)4, (ftnlen)9);
	    dlacpy_("Full", n, n, &dwork[iun2], &ldun2, &a[a_offset], lda, (
		    ftnlen)4);
	    ncol = *n;
	    jwork = iun2;
	}

	if (rcond1 > max(toll1,thresh)) {

/*           The triangular factor r1 is considered to be of full rank. */
/*           Solve for  A  (if requested),  r1*A = un2(1:n,:)  in  A. */

	    if (withc) {
		dtrtrs_("Upper", "NoTranspose", "NonUnit", n, n, &dwork[igal],
			 &ldun2, &a[a_offset], lda, &ierr, (ftnlen)5, (ftnlen)
			11, (ftnlen)7);
		if (ierr > 0) {
		    *info = 3;
		    return 0;
		}
	    }
	    rank = *n;
	} else {

/*           Rank-deficient triangular factor r1.  Use SVD of r1, */
/*           r1 = U*S*V',  also for computing  A  (if requested) from */
/*           r1*A = un2(1:n,:).  Matrix  U  is computed in  DWORK(IU), */
/*           and  V' overwrites  r1.  If  B  is requested, the */
/*           pseudoinverse of  r1  and then of  GaL  are also computed */
/*           in  R(NR3,NR2). */
/*           Workspace: need   c*L*(NOBR-1)*N+N*N+7*N, */
/*                             where  c = 1  if  B and D  are not needed, */
/*                                    c = 2  if  B and D  are needed; */
/*                      prefer larger. */

	    iu = iun2;
	    isv = iu + nn;
	    jwork = isv + *n;
	    if (*m > 0 && withb) {

/*              Save the elementary reflectors used for computing r1, */
/*              if  B, D  are needed. */
/*              Workspace: need   2*L*(NOBR-1)*N+2*N+N*N. */

		ihous = jwork;
		jwork = ihous + ldunn;
		dlacpy_("Lower", &ldun2, n, &dwork[igal], &ldun2, &dwork[
			ihous], &ldun2, (ftnlen)5);
	    } else {
		ihous = igal;
	    }

	    i__1 = *ldwork - jwork + 1;
	    mb02ud_("Not factored", "Left", "NoTranspose", jobp, n, &ncol, &
		    c_b63, &toll1, &rank, &dwork[igal], &ldun2, &dwork[iu], n,
		     &dwork[isv], &a[a_offset], lda, &r__[nr3 + nr2 * r_dim1],
		     ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)12, (ftnlen)4, 
		    (ftnlen)11, (ftnlen)1);
	    if (ierr != 0) {
		*info = 2;
		return 0;
	    }
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

	    if (rank == 0) {
		*(unsigned char *)jobp = 'N';
	    } else if (*m > 0 && withb) {

/*              Compute  pinv(GaL)  in  R(NR3,NR2)  if  B, D  are needed. */
/*              Workspace: need   2*L*(NOBR-1)*N+N*N+3*N; */
/*                         prefer 2*L*(NOBR-1)*N+N*N+2*N+N*NB. */

		i__1 = ldun2 - *n;
		dlaset_("Full", n, &i__1, &c_b66, &c_b66, &r__[nr3 + (nr2 + *
			n) * r_dim1], ldr, (ftnlen)4);
		i__1 = *ldwork - jwork + 1;
		dormqr_("Right", "Transpose", n, &ldun2, n, &dwork[ihous], &
			ldun2, &dwork[itau1], &r__[nr3 + nr2 * r_dim1], ldr, &
			dwork[jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
		maxwrk = max(i__1,i__2);
		if (withco) {

/*                 Save  pinv(GaL)  in  DWORK(IGAL). */

		    dlacpy_("Full", n, &ldun2, &r__[nr3 + nr2 * r_dim1], ldr, 
			    &dwork[igal], n, (ftnlen)4);
		}
		jwork = iun2;
	    }
	    ldw = jwork;
	}

	if (*m > 0 && withb) {

/*           Computation of  B  and  D. */

/*           Compute the reciprocal of the condition number of R_1c. */
/*           Workspace: need L*(NOBR-1)*N+N+3*M*NOBR. */

	    dtrcon_("1-norm", "Upper", "NonUnit", &mnobr, &r__[nr3 + r_dim1], 
		    ldr, &rcond2, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, 
		    (ftnlen)5, (ftnlen)7);

	    toll = *tol;
	    if (toll <= 0.) {
		toll = mnobr * mnobr * eps;
	    }

/*           Compute the right hand side and solve for  K  (in  R_23), */
/*              K*R_1c' = u2'*R_2c', */
/*           where  u2 = Un(:,n+1:L*s),  and  K  is  (Ls-n) x ms. */

	    dgemm_("Transpose", "Transpose", &lnobrn, &mnobr, &lnobr, &c_b63, 
		    &r__[nr2 + (nr2 + *n) * r_dim1], ldr, &r__[nr2 * r_dim1 + 
		    1], ldr, &c_b66, &r__[nr2 + nr3 * r_dim1], ldr, (ftnlen)9,
		     (ftnlen)9);

	    if (rcond2 > max(toll,thresh)) {

/*              The triangular factor R_1c is considered to be of full */
/*              rank. Solve for  K,  K*R_1c' = u2'*R_2c'. */

		dtrsm_("Right", "Upper", "Transpose", "Non-unit", &lnobrn, &
			mnobr, &c_b63, &r__[nr3 + r_dim1], ldr, &r__[nr2 + 
			nr3 * r_dim1], ldr, (ftnlen)5, (ftnlen)5, (ftnlen)9, (
			ftnlen)8);
	    } else {

/*              Rank-deficient triangular factor  R_1c.  Use SVD of  R_1c */
/*              for computing  K  from  K*R_1c' = u2'*R_2c',  where */
/*              R_1c = U1*S1*V1'.  Matrix  U1  is computed in  R_33, */
/*              and  V1'  overwrites  R_1c. */
/*              Workspace: need   L*(NOBR-1)*N+N+6*M*NOBR; */
/*                         prefer larger. */

		isv = ldw;
		jwork = isv + mnobr;
		i__1 = *ldwork - jwork + 1;
		mb02ud_("Not factored", "Right", "Transpose", "No pinv", &
			lnobrn, &mnobr, &c_b63, &toll, &rank11, &r__[nr3 + 
			r_dim1], ldr, &r__[nr3 + nr3 * r_dim1], ldr, &dwork[
			isv], &r__[nr2 + nr3 * r_dim1], ldr, &dwork[jwork], &
			c__1, &dwork[jwork], &i__1, &ierr, (ftnlen)12, (
			ftnlen)5, (ftnlen)9, (ftnlen)7);
		if (ierr != 0) {
		    *info = 2;
		    return 0;
		}
/* Computing MAX */
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
		maxwrk = max(i__1,i__2);
		jwork = ldw;
	    }

/*           Compute the triangular factor of the structured matrix  Q */
/*           and apply the transformations to the matrix  Kexpand,  where */
/*           Q  and  Kexpand  are defined in SLICOT Library routine */
/*           IB01PY.  Compute also the matrices  B,  D. */
/*           Workspace: need   L*(NOBR-1)*N+N+max(L+M*NOBR,L*NOBR+ */
/*                                                max(3*L*NOBR+1,M)); */
/*                      prefer larger. */

	    if (withco) {
		dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &o[
			o_offset], ldo, (ftnlen)4);
	    }
	    i__1 = *ldwork - jwork + 1;
	    ib01py_(meth, jobpy, nobr, n, m, l, &rank, &r__[nr2 + nr2 * 
		    r_dim1], ldr, &dwork[igal], &ldun2, &dwork[itau1], &r__[
		    nr3 + nr2 * r_dim1], ldr, &r__[nr2 + nr3 * r_dim1], ldr, &
		    r__[nr4 + nr2 * r_dim1], ldr, &r__[nr4 + nr3 * r_dim1], 
		    ldr, &b[b_offset], ldb, &d__[d_offset], ldd, tol, &iwork[
		    1], &dwork[jwork], &i__1, iwarn, info, (ftnlen)1, (ftnlen)
		    1);
	    if (*info != 0) {
		return 0;
	    }
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	    rcond4 = dwork[jwork + 1];
	    if (withco) {
		dlacpy_("Full", &lnobr, n, &o[o_offset], ldo, &r__[nr2 + 
			r_dim1], ldr, (ftnlen)4);
	    }

	} else {
	    rcond2 = 1.;
	}

	if (! withco) {
	    rcond3 = 1.;
	    goto L30;
	}
    } else {

/*        For N4SID, set  RCOND2  to one. */

	rcond2 = 1.;
    }

/*     If needed, save the first  n  columns, representing  Gam,  of the */
/*     matrix of left singular vectors,  Un,  in  R_21  and in  O. */

    if (n4sid || withc && ! withal) {
	if (*m > 0) {
	    dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &r__[
		    nr2 + r_dim1], ldr, (ftnlen)4);
	}
	dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &o[o_offset]
		, ldo, (ftnlen)4);
    }

/*     Computations for covariance matrices, and system matrices (N4SID). */
/*     Solve the least squares problems  Gam*Y = R4(1:L*s,1:(2*m+L)*s), */
/*                                       GaL*X = R4(L+1:L*s,:),  where */
/*     GaL = Gam(1:L*(s-1),:),  Gam  has full column rank, and */
/*     R4 = [ R_14' R_24' R_34' R_44L' ],  R_44L = R_44(1:L,:), as */
/*     returned by SLICOT Library routine  IB01ND. */
/*     First, find the  QR  factorization of  Gam,  Gam = Q*R. */
/*     Workspace: need   L*(NOBR-1)*N+Aw+3*N; */
/*                prefer L*(NOBR-1)*N+Aw+2*N+N*NB, where */
/*                Aw = N+N*N,  if  (M = 0  or  JOB = 'C'),  rank(r1) < N, */
/*                             and  METH = 'M'; */
/*                Aw = 0,      otherwise. */

    itau2 = ldw;
    jwork = itau2 + *n;
    i__1 = *ldwork - jwork + 1;
    dgeqrf_(&lnobr, n, &r__[nr2 + r_dim1], ldr, &dwork[itau2], &dwork[jwork], 
	    &i__1, &ierr);

/*     For METH = 'M' or when JOB = 'B' or 'D', transpose */
/*     [ R_14' R_24' R_34' ]'  in the last block-row of  R, obtaining  Z, */
/*     and for METH = 'N' and JOB = 'A' or 'C', use the matrix  Z */
/*     already available in the last block-row of  R,  and then apply */
/*     the transformations, Z <-- Q'*Z. */
/*     Workspace: need   L*(NOBR-1)*N+Aw+2*N+(2*M+L)*NOBR; */
/*                prefer L*(NOBR-1)*N+Aw+2*N+(2*M+L)*NOBR*NB. */

    if (moesp || withb && ! withal) {
	ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		nr4 + r_dim1], ldr, (ftnlen)4);
    }
    i__1 = *ldwork - jwork + 1;
    dormqr_("Left", "Transpose", &lnobr, &lmmnob, n, &r__[nr2 + r_dim1], ldr, 
	    &dwork[itau2], &r__[nr4 + r_dim1], ldr, &dwork[jwork], &i__1, &
	    ierr, (ftnlen)4, (ftnlen)9);

/*     Solve for  Y,  RY = Z  in  Z  and save the transpose of the */
/*     solution  Y  in the second block-column of  R. */

    dtrtrs_("Upper", "NoTranspose", "NonUnit", n, &lmmnob, &r__[nr2 + r_dim1],
	     ldr, &r__[nr4 + r_dim1], ldr, &ierr, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);
    if (ierr > 0) {
	*info = 3;
	return 0;
    }
    ma02ad_("Full", n, &lmmnob, &r__[nr4 + r_dim1], ldr, &r__[nr2 * r_dim1 + 
	    1], ldr, (ftnlen)4);
    nr4mn = nr4 - *n;
    nr4pl = nr4 + *l;
    nrow = lmmnol;

/*     SHIFT is .TRUE. if some columns of  R_14 : R_44L  should be */
/*     shifted to the right, to avoid overwriting useful information. */

    shift = *m == 0 && lnobr < n2;

    if (rcond1 > max(toll1,thresh)) {

/*        The triangular factor  r1  of  GaL  (GaL = Q1*r1)  is */
/*        considered to be of full rank. */

/*        Transpose  [ R_14' R_24' R_34' R_44L' ]'(:,L+1:L*s)  in the */
/*        last block-row of  R  (beginning with the  (L+1)-th  row), */
/*        obtaining  Z1,  and then apply the transformations, */
/*        Z1 <-- Q1'*Z1. */
/*        Workspace: need   L*(NOBR-1)*N+Aw+2*N+ (2*M+L)*NOBR + L; */
/*                   prefer L*(NOBR-1)*N+Aw+2*N+((2*M+L)*NOBR + L)*NB. */

	ma02ad_("Full", &lmmnol, &ldun2, &r__[nr4pl * r_dim1 + 1], ldr, &r__[
		nr4pl + r_dim1], ldr, (ftnlen)4);
	i__1 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &ldun2, &lmmnol, n, &dwork[igal], &ldun2,
		 &dwork[itau1], &r__[nr4pl + r_dim1], ldr, &dwork[jwork], &
		i__1, &ierr, (ftnlen)4, (ftnlen)9);

/*        Solve for  X,  r1*X = Z1  in  Z1,  and copy the transpose of  X */
/*        into the last part of the third block-column of  R. */

	dtrtrs_("Upper", "NoTranspose", "NonUnit", n, &lmmnol, &dwork[igal], &
		ldun2, &r__[nr4pl + r_dim1], ldr, &ierr, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
	if (ierr > 0) {
	    *info = 3;
	    return 0;
	}

	if (shift) {
	    nr4mn = nr4;

	    for (i__ = *l - 1; i__ >= 0; --i__) {
		dcopy_(&lmmnol, &r__[(nr4 + i__) * r_dim1 + 1], &c__1, &r__[(
			nr4 + *n + i__) * r_dim1 + 1], &c__1);
/* L10: */
	    }

	}
	ma02ad_("Full", n, &lmmnol, &r__[nr4pl + r_dim1], ldr, &r__[nr4mn * 
		r_dim1 + 1], ldr, (ftnlen)4);
	nrow = 0;
    }

    if (n4sid || nrow > 0) {

/*        METH = 'N'  or rank-deficient triangular factor r1. */
/*        For  METH = 'N',  use SVD of  r1,  r1 = U*S*V', for computing */
/*        X'  from  X'*GaL' = Z1',  if  rank(r1) < N.  Matrix  U  is */
/*        computed in  DWORK(IU)  and  V'  overwrites  r1.  Then, the */
/*        pseudoinverse of  GaL  is determined in  R(NR4+L,NR2). */
/*        For METH = 'M', the pseudoinverse of  GaL  is already available */
/*        if  M > 0  and  B  is requested;  otherwise, the SVD of  r1  is */
/*        available in  DWORK(IU),  DWORK(ISV),  and  DWORK(IGAL). */
/*        Workspace for N4SID: need   2*L*(NOBR-1)*N+N*N+8*N; */
/*                             prefer larger. */

	if (moesp) {
	    *(unsigned char *)fact = 'F';
	    if (*m > 0 && withb) {
		dlacpy_("Full", n, &ldun2, &dwork[igal], n, &r__[nr4pl + nr2 *
			 r_dim1], ldr, (ftnlen)4);
	    }
	} else {

/*           Save the elementary reflectors used for computing r1. */

	    ihous = jwork;
	    dlacpy_("Lower", &ldun2, n, &dwork[igal], &ldun2, &dwork[ihous], &
		    ldun2, (ftnlen)5);
	    *(unsigned char *)fact = 'N';
	    iu = ihous + ldunn;
	    isv = iu + nn;
	    jwork = isv + *n;
	}

	i__1 = *ldwork - jwork + 1;
	mb02ud_(fact, "Right", "Transpose", jobp, &nrow, n, &c_b63, &toll1, &
		rank, &dwork[igal], &ldun2, &dwork[iu], n, &dwork[isv], &r__[
		nr4pl * r_dim1 + 1], ldr, &r__[nr4pl + nr2 * r_dim1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)1, (ftnlen)5, (ftnlen)9, (
		ftnlen)1);
	if (nrow > 0) {
	    if (shift) {
		nr4mn = nr4;
		dlacpy_("Full", &lmmnol, l, &r__[nr4 * r_dim1 + 1], ldr, &r__[
			(nr4 - *l) * r_dim1 + 1], ldr, (ftnlen)4);
		dlacpy_("Full", &lmmnol, n, &r__[nr4pl * r_dim1 + 1], ldr, &
			r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)4);
		dlacpy_("Full", &lmmnol, l, &r__[(nr4 - *l) * r_dim1 + 1], 
			ldr, &r__[(nr4 + *n) * r_dim1 + 1], ldr, (ftnlen)4);
	    } else {
		dlacpy_("Full", &lmmnol, n, &r__[nr4pl * r_dim1 + 1], ldr, &
			r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)4);
	    }
	}

	if (n4sid) {
	    if (ierr != 0) {
		*info = 2;
		return 0;
	    }
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Compute  pinv(GaL)  in  R(NR4+L,NR2). */
/*           Workspace: need   2*L*(NOBR-1)*N+3*N; */
/*                      prefer 2*L*(NOBR-1)*N+2*N+N*NB. */

	    jwork = iu;
	    i__1 = ldun2 - *n;
	    dlaset_("Full", n, &i__1, &c_b66, &c_b66, &r__[nr4pl + (nr2 + *n) 
		    * r_dim1], ldr, (ftnlen)4);
	    i__1 = *ldwork - jwork + 1;
	    dormqr_("Right", "Transpose", n, &ldun2, n, &dwork[ihous], &ldun2,
		     &dwork[itau1], &r__[nr4pl + nr2 * r_dim1], ldr, &dwork[
		    jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	}
    }

/*     For METH = 'N', find part of the solution (corresponding to A */
/*     and C) and, optionally, for both  METH = 'M',  or  METH = 'N', */
/*     find the residual of the least squares problem that gives the */
/*     covariances,  M*V = N,  where */
/*         (     R_11 ) */
/*     M = (  Y'      ),  N = (  X'   R4'(:,1:L) ),  V = V(n+m*s, n+L), */
/*         (  0   0   ) */
/*     with  M((2*m+L)*s+L, n+m*s),  N((2*m+L)*s+L, n+L),  R4'  being */
/*     stored in the last block-column of  R.  The last  L  rows of  M */
/*     are not explicitly considered. Note that, for efficiency, the */
/*     last  m*s  columns of  M  are in the first positions of arrray  R. */
/*     This permutation does not affect the residual, only the */
/*     solution.  (The solution is not needed for METH = 'M'.) */
/*     Note that R_11 corresponds to the future outputs for both */
/*     METH = 'M', or METH = 'N' approaches.  (For  METH = 'N',  the */
/*     first two block-columns have been interchanged.) */
/*     For  METH = 'N',  A and C are obtained as follows: */
/*     [ A'  C' ] = V(m*s+1:m*s+n,:). */

/*     First, find the  QR  factorization of  Y'(m*s+1:(2*m+L)*s,:) */
/*     and apply the transformations to the corresponding part of N. */
/*     Compress the workspace for N4SID by moving the scalar reflectors */
/*     corresponding to  Q. */
/*     Workspace: need   d*N+2*N; */
/*                prefer d*N+N+N*NB; */
/*     where  d = 0,  for  MOESP,  and  d = 1,  for  N4SID. */

    if (moesp) {
	itau = 1;
    } else {
	dcopy_(n, &dwork[itau2], &c__1, &dwork[1], &c__1);
	itau = *n + 1;
    }

    jwork = itau + *n;
    i__1 = *ldwork - jwork + 1;
    dgeqrf_(&lmnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[itau], &dwork[
	    jwork], &i__1, &ierr);

/*     Workspace: need   d*N+N+(N+L); */
/*                prefer d*N+N+(N+L)*NB. */

    i__1 = *ldwork - jwork + 1;
    dormqr_("Left", "Transpose", &lmnobr, &npl, n, &r__[nr2 + nr2 * r_dim1], 
	    ldr, &dwork[itau], &r__[nr2 + nr4mn * r_dim1], ldr, &dwork[jwork],
	     &i__1, &ierr, (ftnlen)4, (ftnlen)9);

    i__1 = *l - 1;
    i__2 = *l - 1;
    dlaset_("Lower", &i__1, &i__2, &c_b66, &c_b66, &r__[nr4 + 1 + nr4 * 
	    r_dim1], ldr, (ftnlen)5);

/*     Now, matrix  M  with permuted block-columns has been */
/*     triangularized. */
/*     Compute the reciprocal of the condition number of its */
/*     triangular factor in  R(1:m*s+n,1:m*s+n). */
/*     Workspace: need d*N+3*(M*NOBR+N). */

    jwork = itau;
    dtrcon_("1-norm", "Upper", "NonUnit", &mnobrn, &r__[r_offset], ldr, &
	    rcond3, &dwork[jwork], &iwork[1], info, (ftnlen)6, (ftnlen)5, (
	    ftnlen)7);

    toll = *tol;
    if (toll <= 0.) {
	toll = mnobrn * mnobrn * eps;
    }
    if (rcond3 > max(toll,thresh)) {

/*        The triangular factor is considered to be of full rank. */
/*        Solve for  V(m*s+1:m*s+n,:),  giving  [ A'  C' ]. */

	fullr = TRUE_;
	rankm = mnobrn;
	if (n4sid) {
	    dtrsm_("Left", "Upper", "NoTranspose", "NonUnit", n, &npl, &c_b63,
		     &r__[nr2 + nr2 * r_dim1], ldr, &r__[nr2 + nr4mn * r_dim1]
		    , ldr, (ftnlen)4, (ftnlen)5, (ftnlen)11, (ftnlen)7);
	}
    } else {
	fullr = FALSE_;

/*        Use QR factorization (with pivoting). For METH = 'N', save */
/*        (and then restore) information about the QR factorization of */
/*        Gam,  for later use. Note that  R_11  could be modified by */
/*        MB03OD, but the corresponding part of  N  is also modified */
/*        accordingly. */
/*        Workspace: need   d*N+4*(M*NOBR+N)+1; */
/*                   prefer d*N+3*(M*NOBR+N)+(M*NOBR+N+1)*NB. */

	i__1 = mnobrn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iwork[i__] = 0;
/* L20: */
	}

	if (n4sid && (*m > 0 && withb)) {
	    dlacpy_("Full", &lnobr, n, &r__[nr2 + r_dim1], ldr, &r__[nr4 + 
		    r_dim1], ldr, (ftnlen)4);
	}
	jwork = itau + mnobrn;
	i__1 = mnobrn - 1;
	dlaset_("Lower", &i__1, &mnobrn, &c_b66, &c_b66, &r__[r_dim1 + 2], 
		ldr, (ftnlen)5);
	i__1 = *ldwork - jwork + 1;
	mb03od_("QR", &mnobrn, &mnobrn, &r__[r_offset], ldr, &iwork[1], &toll,
		 &svlmax, &dwork[itau], &rankm, sval, &dwork[jwork], &i__1, &
		ierr, (ftnlen)2);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__2);

/*        Workspace: need   d*N+M*NOBR+N+N+L; */
/*                   prefer d*N+M*NOBR+N+(N+L)*NB. */

	i__1 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &mnobrn, &npl, &mnobrn, &r__[r_offset], 
		ldr, &dwork[itau], &r__[nr4mn * r_dim1 + 1], ldr, &dwork[
		jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__2);
    }

    if (withco) {

/*        The residual (transposed) of the least squares solution */
/*        (multiplied by a matrix with orthogonal columns) is stored */
/*        in the rows  RANKM+1:(2*m+L)*s+L  of V,  and it should be */
/*        squared-up for getting the covariance matrices. (Generically, */
/*        RANKM = m*s+n.) */

	rnrm = 1. / (doublereal) (*nsmpl);
	if (moesp) {
	    i__1 = lmmnol - rankm;
	    dsyrk_("Upper", "Transpose", &npl, &i__1, &rnrm, &r__[rankm + 1 + 
		    nr4mn * r_dim1], ldr, &c_b66, &r__[r_offset], ldr, (
		    ftnlen)5, (ftnlen)9);
	    dlacpy_("Upper", n, n, &r__[r_offset], ldr, &q[q_offset], ldq, (
		    ftnlen)5);
	    dlacpy_("Full", n, l, &r__[(*n + 1) * r_dim1 + 1], ldr, &s[
		    s_offset], lds, (ftnlen)4);
	    dlacpy_("Upper", l, l, &r__[*n + 1 + (*n + 1) * r_dim1], ldr, &ry[
		    ry_offset], ldry, (ftnlen)5);
	} else {
	    i__1 = lmmnol - rankm;
	    dsyrk_("Upper", "Transpose", &npl, &i__1, &rnrm, &r__[rankm + 1 + 
		    nr4mn * r_dim1], ldr, &c_b66, &dwork[jwork], &npl, (
		    ftnlen)5, (ftnlen)9);
	    dlacpy_("Upper", n, n, &dwork[jwork], &npl, &q[q_offset], ldq, (
		    ftnlen)5);
	    dlacpy_("Full", n, l, &dwork[jwork + *n * npl], &npl, &s[s_offset]
		    , lds, (ftnlen)4);
	    dlacpy_("Upper", l, l, &dwork[jwork + *n * (npl + 1)], &npl, &ry[
		    ry_offset], ldry, (ftnlen)5);
	}
	ma02ed_("Upper", n, &q[q_offset], ldq, (ftnlen)5);
	ma02ed_("Upper", l, &ry[ry_offset], ldry, (ftnlen)5);

/*        Check the magnitude of the residual. */

	i__1 = lmmnol - rankm;
	rnrm = dlange_("1-norm", &i__1, &npl, &r__[rankm + 1 + nr4mn * r_dim1]
		, ldr, &dwork[jwork], (ftnlen)6);
	if (rnrm < thresh) {
	    *iwarn = 5;
	}
    }

    if (n4sid) {
	if (! fullr) {
	    *iwarn = 4;

/*           Compute part of the solution of the least squares problem, */
/*           M*V = N,  for the rank-deficient problem. */
/*           Remark: this computation should not be performed before the */
/*           symmetric updating operation above. */
/*           Workspace: need   M*NOBR+3*N+L; */
/*                      prefer larger. */

	    i__1 = *ldwork - jwork + 1;
	    mb03od_("No QR", n, n, &r__[nr2 + nr2 * r_dim1], ldr, &iwork[1], &
		    toll1, &svlmax, &dwork[itau], &rankm, sval, &dwork[jwork],
		     &i__1, &ierr, (ftnlen)5);
	    i__1 = *ldwork - jwork + 1;
	    mb02qy_(n, n, &npl, &rankm, &r__[nr2 + nr2 * r_dim1], ldr, &iwork[
		    1], &r__[nr2 + nr4mn * r_dim1], ldr, &dwork[itau + mnobr],
		     &dwork[jwork], &i__1, info);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	    jwork = itau;
	    if (*m > 0 && withb) {
		dlacpy_("Full", &lnobr, n, &r__[nr4 + r_dim1], ldr, &r__[nr2 
			+ r_dim1], ldr, (ftnlen)4);
	    }
	}

	if (withc) {

/*           Obtain  A  and  C,  noting that block-permutations have been */
/*           implicitly used. */

	    ma02ad_("Full", n, n, &r__[nr2 + nr4mn * r_dim1], ldr, &a[
		    a_offset], lda, (ftnlen)4);
	    ma02ad_("Full", n, l, &r__[nr2 + (nr4mn + *n) * r_dim1], ldr, &
		    c__[c_offset], ldc, (ftnlen)4);
	} else {

/*           Use the given  A  and  C. */

	    ma02ad_("Full", n, n, &a[a_offset], lda, &r__[nr2 + nr4mn * 
		    r_dim1], ldr, (ftnlen)4);
	    ma02ad_("Full", l, n, &c__[c_offset], ldc, &r__[nr2 + (nr4mn + *n)
		     * r_dim1], ldr, (ftnlen)4);
	}

	if (*m > 0 && withb) {

/*           Obtain  B  and  D. */
/*           First, compute the transpose of the matrix K as */
/*           N(1:m*s,:) - M(1:m*s,m*s+1:m*s+n)*[A'  C'],  in the first */
/*           m*s  rows of  R(1,NR4MN). */

	    dgemm_("NoTranspose", "NoTranspose", &mnobr, &npl, n, &c_b173, &
		    r__[nr2 * r_dim1 + 1], ldr, &r__[nr2 + nr4mn * r_dim1], 
		    ldr, &c_b63, &r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)11, (
		    ftnlen)11);

/*           Denote   M = pinv(GaL)  and construct */

/*                    [ [ A ]   -1   ]                      [ R ] */
/*           and  L = [ [   ]  R   0 ] Q',  where Gam = Q * [   ]. */
/*                    [ [ C ]        ]                      [ 0 ] */

/*           Then, solve the least squares problem. */

	    dlacpy_("Full", n, n, &a[a_offset], lda, &r__[nr2 + nr4 * r_dim1],
		     ldr, (ftnlen)4);
	    dlacpy_("Full", l, n, &c__[c_offset], ldc, &r__[nr2 + *n + nr4 * 
		    r_dim1], ldr, (ftnlen)4);
	    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", &npl, n, &
		    c_b63, &r__[nr2 + r_dim1], ldr, &r__[nr2 + nr4 * r_dim1], 
		    ldr, (ftnlen)5, (ftnlen)5, (ftnlen)11, (ftnlen)7);
	    dlaset_("Full", &npl, &lnobrn, &c_b66, &c_b66, &r__[nr2 + (nr4 + *
		    n) * r_dim1], ldr, (ftnlen)4);

/*           Workspace: need 2*N+L; prefer N + (N+L)*NB. */

	    i__1 = *ldwork - jwork + 1;
	    dormqr_("Right", "Transpose", &npl, &lnobr, n, &r__[nr2 + r_dim1],
		     ldr, &dwork[1], &r__[nr2 + nr4 * r_dim1], ldr, &dwork[
		    jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);

/*           Obtain the matrix  K  by transposition, and find  B  and  D. */
/*           Workspace: need   NOBR*(M*(N+L))**2+M*NOBR*(N+L)+ */
/*                             max((N+L)**2,4*M*(N+L)+1); */
/*                      prefer larger. */

	    ma02ad_("Full", &mnobr, &npl, &r__[nr4mn * r_dim1 + 1], ldr, &r__[
		    nr2 + nr3 * r_dim1], ldr, (ftnlen)4);
/* Computing 2nd power */
	    i__1 = npl;
	    ix = mnobr * (i__1 * i__1) * *m + 1;
	    jwork = ix + mnobr * npl;
	    i__1 = mnobr * npl;
	    i__2 = *ldwork - jwork + 1;
	    ib01px_(jobpy, nobr, n, m, l, &r__[r_offset], ldr, &o[o_offset], 
		    ldo, &r__[nr2 + nr4 * r_dim1], ldr, &r__[nr4pl + nr2 * 
		    r_dim1], ldr, &r__[nr2 + nr3 * r_dim1], ldr, &dwork[1], &
		    i__1, &dwork[ix], &b[b_offset], ldb, &d__[d_offset], ldd, 
		    tol, &iwork[1], &dwork[jwork], &i__2, &iwarnl, info, (
		    ftnlen)1);
	    if (*info != 0) {
		return 0;
	    }
	    *iwarn = max(*iwarn,iwarnl);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	    rcond4 = dwork[jwork + 1];

	}
    }

L30:

/*     Return optimal workspace in  DWORK(1)  and reciprocal condition */
/*     numbers in the next locations. */

    dwork[1] = (doublereal) maxwrk;
    dwork[2] = rcond1;
    dwork[3] = rcond2;
    dwork[4] = rcond3;
    dwork[5] = rcond4;
    return 0;

/* *** Last line of IB01PD *** */
} /* ib01pd_ */

