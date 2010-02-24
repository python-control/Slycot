      SUBROUTINE SG02AD( DICO, JOBB, FACT, UPLO, JOBL, SCAL, SORT, ACC,
     $                   N, M, P, A, LDA, E, LDE, B, LDB, Q, LDQ, R,
     $                   LDR, L, LDL, RCONDU, X, LDX, ALFAR, ALFAI,
     $                   BETA, S, LDS, T, LDT, U, LDU, TOL, IWORK,
     $                   DWORK, LDWORK, BWORK, IWARN, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To solve for X either the continuous-time algebraic Riccati
C     equation
C                                   -1
C        Q + A'XE + E'XA - (L+E'XB)R  (L+E'XB)' = 0 ,              (1)
C
C     or the discrete-time algebraic Riccati equation
C                                        -1
C        E'XE = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q ,         (2)
C
C     where A, E, B, Q, R, and L are N-by-N, N-by-N, N-by-M, N-by-N,
C     M-by-M and N-by-M matrices, respectively, such that Q = C'C,
C     R = D'D and L = C'D; X is an N-by-N symmetric matrix.
C     The routine also returns the computed values of the closed-loop
C     spectrum of the system, i.e., the stable eigenvalues
C     lambda(1),...,lambda(N) of the pencil (A - BF,E), where F is
C     the optimal gain matrix,
C             -1
C        F = R  (L+E'XB)' ,        for (1),
C
C     and
C                    -1
C        F = (R+B'XB)  (L+A'XB)' , for (2).
C                              -1
C     Optionally, matrix G = BR  B' may be given instead of B and R.
C     Other options include the case with Q and/or R given in a
C     factored form, Q = C'C, R = D'D, and with L a zero matrix.
C
C     The routine uses the method of deflating subspaces, based on
C     reordering the eigenvalues in a generalized Schur matrix pair.
C
C     It is assumed that E is nonsingular, but this condition is not
C     checked. Note that the definition (1) of the continuous-time
C     algebraic Riccati equation, and the formula for the corresponding
C     optimal gain matrix, require R to be nonsingular, but the
C     associated linear quadratic optimal problem could have a unique
C     solution even when matrix R is singular, under mild assumptions
C     (see METHOD). The routine SG02AD works accordingly in this case.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of Riccati equation to be solved as
C             follows:
C             = 'C':  Equation (1), continuous-time case;
C             = 'D':  Equation (2), discrete-time case.
C
C     JOBB    CHARACTER*1
C             Specifies whether or not the matrix G is given, instead
C             of the matrices B and R, as follows:
C             = 'B':  B and R are given;
C             = 'G':  G is given.
C
C     FACT    CHARACTER*1
C             Specifies whether or not the matrices Q and/or R (if
C             JOBB = 'B') are factored, as follows:
C             = 'N':  Not factored, Q and R are given;
C             = 'C':  C is given, and Q = C'C;
C             = 'D':  D is given, and R = D'D;
C             = 'B':  Both factors C and D are given, Q = C'C, R = D'D.
C
C     UPLO    CHARACTER*1
C             If JOBB = 'G', or FACT = 'N', specifies which triangle of
C             the matrices G, or Q and R, is stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
C
C     JOBL    CHARACTER*1
C             Specifies whether or not the matrix L is zero, as follows:
C             = 'Z':  L is zero;
C             = 'N':  L is nonzero.
C             JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed.
C             SLICOT Library routine SB02MT should be called just before
C             SG02AD, for obtaining the results when JOBB = 'G' and
C             JOBL = 'N'.
C
C     SCAL    CHARACTER*1
C             If JOBB = 'B', specifies whether or not a scaling strategy
C             should be used to scale Q, R, and L, as follows:
C             = 'G':  General scaling should be used;
C             = 'N':  No scaling should be used.
C             SCAL is not used if JOBB = 'G'.
C
C     SORT    CHARACTER*1
C             Specifies which eigenvalues should be obtained in the top
C             of the generalized Schur form, as follows:
C             = 'S':  Stable   eigenvalues come first;
C             = 'U':  Unstable eigenvalues come first.
C
C     ACC     CHARACTER*1
C             Specifies whether or not iterative refinement should be
C             used to solve the system of algebraic equations giving
C             the solution matrix X, as follows:
C             = 'R':  Use iterative refinement;
C             = 'N':  Do not use iterative refinement.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices A, E, Q, and X, and the number of rows of the
C             matrices B and L.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs. If JOBB = 'B', M is the
C             order of the matrix R, and the number of columns of the
C             matrix B.  M >= 0.
C             M is not used if JOBB = 'G'.
C
C     P       (input) INTEGER
C             The number of system outputs. If FACT = 'C' or 'D' or 'B',
C             P is the number of rows of the matrices C and/or D.
C             P >= 0.
C             Otherwise, P is not used.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the descriptor system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             The leading N-by-N part of this array must contain the
C             matrix E of the descriptor system.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,*)
C             If JOBB = 'B', the leading N-by-M part of this array must
C             contain the input matrix B of the system.
C             If JOBB = 'G', the leading N-by-N upper triangular part
C             (if UPLO = 'U') or lower triangular part (if UPLO = 'L')
C             of this array must contain the upper triangular part or
C             lower triangular part, respectively, of the matrix
C                   -1
C             G = BR  B'. The stricly lower triangular part (if
C             UPLO = 'U') or stricly upper triangular part (if
C             UPLO = 'L') is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N)
C             If FACT = 'N' or 'D', the leading N-by-N upper triangular
C             part (if UPLO = 'U') or lower triangular part (if UPLO =
C             'L') of this array must contain the upper triangular part
C             or lower triangular part, respectively, of the symmetric
C             state weighting matrix Q. The stricly lower triangular
C             part (if UPLO = 'U') or stricly upper triangular part (if
C             UPLO = 'L') is not referenced.
C             If FACT = 'C' or 'B', the leading P-by-N part of this
C             array must contain the output matrix C of the system.
C             If JOBB = 'B' and SCAL = 'G', then Q is modified
C             internally, but is restored on exit.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= MAX(1,N) if FACT = 'N' or 'D';
C             LDQ >= MAX(1,P) if FACT = 'C' or 'B'.
C
C     R       (input) DOUBLE PRECISION array, dimension (LDR,*)
C             If FACT = 'N' or 'C', the leading M-by-M upper triangular
C             part (if UPLO = 'U') or lower triangular part (if UPLO =
C             'L') of this array must contain the upper triangular part
C             or lower triangular part, respectively, of the symmetric
C             input weighting matrix R. The stricly lower triangular
C             part (if UPLO = 'U') or stricly upper triangular part (if
C             UPLO = 'L') is not referenced.
C             If FACT = 'D' or 'B', the leading P-by-M part of this
C             array must contain the direct transmission matrix D of the
C             system.
C             If JOBB = 'B' and SCAL = 'G', then R is modified
C             internally, but is restored on exit.
C             If JOBB = 'G', this array is not referenced.
C
C     LDR     INTEGER
C             The leading dimension of array R.
C             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C';
C             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B';
C             LDR >= 1        if JOBB = 'G'.
C
C     L       (input) DOUBLE PRECISION array, dimension (LDL,*)
C             If JOBL = 'N' and JOBB = 'B', the leading N-by-M part of
C             this array must contain the cross weighting matrix L.
C             If JOBB = 'B' and SCAL = 'G', then L is modified
C             internally, but is restored on exit.
C             If JOBL = 'Z' or JOBB = 'G', this array is not referenced.
C
C     LDL     INTEGER
C             The leading dimension of array L.
C             LDL >= MAX(1,N) if JOBL = 'N' and JOBB = 'B';
C             LDL >= 1        if JOBL = 'Z' or  JOBB = 'G'.
C
C     RCONDU  (output) DOUBLE PRECISION
C             If N > 0 and INFO = 0 or INFO = 7, an estimate of the
C             reciprocal of the condition number (in the 1-norm) of
C             the N-th order system of algebraic equations from which
C             the solution matrix X is obtained.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N)
C             If INFO = 0, the leading N-by-N part of this array
C             contains the solution matrix X of the problem.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= MAX(1,N).
C
C     ALFAR   (output) DOUBLE PRECISION array, dimension (2*N)
C     ALFAI   (output) DOUBLE PRECISION array, dimension (2*N)
C     BETA    (output) DOUBLE PRECISION array, dimension (2*N)
C             The generalized eigenvalues of the 2N-by-2N matrix pair,
C             ordered as specified by SORT (if INFO = 0, or INFO >= 5).
C             For instance, if SORT = 'S', the leading N elements of
C             these arrays contain the closed-loop spectrum of the
C             system. Specifically,
C                lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for
C             k = 1,2,...,N.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,*)
C             The leading 2N-by-2N part of this array contains the
C             ordered real Schur form S of the first matrix in the
C             reduced matrix pencil associated to the optimal problem,
C             corresponding to the scaled Q, R, and L, if JOBB = 'B'
C             and SCAL = 'G'. That is,
C
C                    (S   S  )
C                    ( 11  12)
C                S = (       ),
C                    (0   S  )
C                    (     22)
C
C             where S  , S   and S   are N-by-N matrices.
C                    11   12      22
C             Array S must have 2*N+M columns if JOBB = 'B', and 2*N
C             columns, otherwise.
C
C     LDS     INTEGER
C             The leading dimension of array S.
C             LDS >= MAX(1,2*N+M) if JOBB = 'B';
C             LDS >= MAX(1,2*N)   if JOBB = 'G'.
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,2*N)
C             The leading 2N-by-2N part of this array contains the
C             ordered upper triangular form T of the second matrix in
C             the reduced matrix pencil associated to the optimal
C             problem, corresponding to the scaled Q, R, and L, if
C             JOBB = 'B' and SCAL = 'G'. That is,
C
C                    (T   T  )
C                    ( 11  12)
C                T = (       ),
C                    (0   T  )
C                    (     22)
C
C             where T  , T   and T   are N-by-N matrices.
C                    11   12      22
C
C     LDT     INTEGER
C             The leading dimension of array T.
C             LDT >= MAX(1,2*N+M) if JOBB = 'B';
C             LDT >= MAX(1,2*N)   if JOBB = 'G'.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N)
C             The leading 2N-by-2N part of this array contains the right
C             transformation matrix U which reduces the 2N-by-2N matrix
C             pencil to the ordered generalized real Schur form (S,T).
C             That is,
C
C                    (U   U  )
C                    ( 11  12)
C                U = (       ),
C                    (U   U  )
C                    ( 21  22)
C
C             where U  , U  , U   and U   are N-by-N matrices.
C                    11   12   21      22
C             If JOBB = 'B' and SCAL = 'G', then U corresponds to the
C             scaled pencil. If a basis for the stable deflating
C             subspace of the original problem is needed, then the
C             submatrix U   must be multiplied by the scaling factor
C                        21
C             contained in DWORK(4).
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,2*N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the original matrix pencil, specifically of the triangular
C             M-by-M factor obtained during the reduction process. If
C             the user sets TOL > 0, then the given value of TOL is used
C             as a lower bound for the reciprocal condition number of
C             that matrix; a matrix whose estimated condition number is
C             less than 1/TOL is considered to be nonsingular. If the
C             user sets TOL <= 0, then a default tolerance, defined by
C             TOLDEF = EPS, is used instead, where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH).
C             This parameter is not referenced if JOBB = 'G'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK >= MAX(1,M,2*N) if JOBB = 'B';
C             LIWORK >= MAX(1,2*N)   if JOBB = 'G'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK. If JOBB = 'B' and N > 0, DWORK(2) returns the
C             reciprocal of the condition number of the M-by-M bottom
C             right lower triangular matrix obtained while compressing
C             the matrix pencil of order 2N+M to obtain a pencil of
C             order 2N. If ACC = 'R', and INFO = 0 or INFO = 7, DWORK(3)
C             returns the reciprocal pivot growth factor (see SLICOT
C             Library routine MB02PD) for the LU factorization of the
C             coefficient matrix of the system of algebraic equations
C             giving the solution matrix X; if DWORK(3) is much
C             less than 1, then the computed X and RCONDU could be
C             unreliable. If INFO = 0 or INFO = 7, DWORK(4) returns the
C             scaling factor used to scale Q, R, and L. DWORK(4) is set
C             to 1 if JOBB = 'G' or SCAL = 'N'.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G';
C             LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'.
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  the computed solution may be inaccurate due to poor
C                   scaling or eigenvalues too close to the boundary of
C                   the stability domain (the imaginary axis, if
C                   DICO = 'C', or the unit circle, if DICO = 'D').
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the computed extended matrix pencil is singular,
C                   possibly due to rounding errors;
C             = 2:  if the QZ algorithm failed;
C             = 3:  if reordering of the generalized eigenvalues failed;
C             = 4:  if after reordering, roundoff changed values of
C                   some complex eigenvalues so that leading eigenvalues
C                   in the generalized Schur form no longer satisfy the
C                   stability condition; this could also be caused due
C                   to scaling;
C             = 5:  if the computed dimension of the solution does not
C                   equal N;
C             = 6:  if the spectrum is too close to the boundary of
C                   the stability domain;
C             = 7:  if a singular matrix was encountered during the
C                   computation of the solution matrix X.
C
C     METHOD
C
C     The routine uses a variant of the method of deflating subspaces
C     proposed by van Dooren [1]. See also [2], [3], [4].
C     It is assumed that E is nonsingular, the triple (E,A,B) is
C     strongly stabilizable and detectable (see [3]); if, in addition,
C
C        -    [ Q   L ]
C        R := [       ] >= 0 ,
C             [ L'  R ]
C
C     then the pencils
C
C           discrete-time                   continuous-time
C
C     |A   0   B|     |E   0   0|    |A   0   B|     |E   0   0|
C     |Q  -E'  L| - z |0  -A'  0| ,  |Q   A'  L| - s |0  -E'  0| ,   (3)
C     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0|
C
C     are dichotomic, i.e., they have no eigenvalues on the boundary of
C     the stability domain. The above conditions are sufficient for
C     regularity of these pencils. A necessary condition is that
C     rank([ B'  L'  R']') = m.
C
C     Under these assumptions the algebraic Riccati equation is known to
C     have a unique non-negative definite solution.
C     The first step in the method of deflating subspaces is to form the
C     extended matrices in (3), of order 2N + M. Next, these pencils are
C     compressed to a form of order 2N (see [1])
C
C        lambda x A  - B .
C                  f    f
C
C     This generalized eigenvalue problem is then solved using the QZ
C     algorithm and the stable deflating subspace Ys is determined.
C     If [Y1'|Y2']' is a basis for Ys, then the required solution is
C                       -1
C            X = Y2 x Y1  .
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         A Generalized Eigenvalue Approach for Solving Riccati
C         Equations.
C         SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981.
C
C     [2] Arnold, III, W.F. and Laub, A.J.
C         Generalized Eigenproblem Algorithms and Software for
C         Algebraic Riccati Equations.
C         Proc. IEEE, 72, 1746-1754, 1984.
C
C     [3] Mehrmann, V.
C         The Autonomous Linear Quadratic Control Problem. Theory and
C         Numerical Solution.
C         Lect. Notes in Control and Information Sciences, vol. 163,
C         Springer-Verlag, Berlin, 1991.
C
C     [4] Sima, V.
C         Algorithms for Linear-Quadratic Optimization.
C         Pure and Applied Mathematics: A Series of Monographs and
C         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996.
C
C     NUMERICAL ASPECTS
C
C     This routine is particularly suited for systems where the matrix R
C     is ill-conditioned, or even singular.
C
C     FURTHER COMMENTS
C
C     To obtain a stabilizing solution of the algebraic Riccati
C     equations set SORT = 'S'.
C
C     The routine can also compute the anti-stabilizing solutions of
C     the algebraic Riccati equations, by specifying SORT = 'U'.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, June 2002.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, September 2002,
C     December 2002.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, discrete-time system, optimal regulator, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, P1, FOUR
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                    P1 = 0.1D0, FOUR = 4.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         ACC, DICO, FACT, JOBB, JOBL, SCAL, SORT, UPLO
      INTEGER           INFO, IWARN, LDA, LDB, LDE, LDL, LDQ, LDR, LDS,
     $                  LDT, LDU, LDWORK, LDX, M, N, P
      DOUBLE PRECISION  RCONDU, TOL
C     .. Array Arguments ..
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), ALFAI(*), ALFAR(*), B(LDB,*), BETA(*),
     $                  DWORK(*), E(LDE,*), L(LDL,*), Q(LDQ,*),
     $                  R(LDR,*), S(LDS,*), T(LDT,*), U(LDU,*), X(LDX,*)
C     .. Local Scalars ..
      CHARACTER         EQUED, QTYPE, RTYPE
      LOGICAL           COLEQU, DISCR, LFACB, LFACN, LFACQ, LFACR,
     $                  LJOBB, LJOBL, LJOBLN, LSCAL, LSORT, LUPLO,
     $                  REFINE, ROWEQU
      INTEGER           I, INFO1, IW, IWB, IWC, IWF, IWR, J, LDW, MP,
     $                  NDIM, NN, NNM, NP, NP1, WRKOPT
      DOUBLE PRECISION  ASYM, EPS, PIVOTU, RCONDL, RNORM, SCALE, SEPS,
     $                  U12M, UNORM
C     .. External Functions ..
      LOGICAL           LSAME, SB02OU, SB02OV, SB02OW
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANSY
      EXTERNAL          DLAMCH, DLANGE, DLANSY, LSAME, SB02OU, SB02OV,
     $                  SB02OW
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGECON, DGEMM, DGEQRF, DGGES,
     $                  DLACPY, DLASCL, DLASET, DORGQR, DSCAL, DSWAP,
     $                  MB01SD, MB02PD, MB02VD, SB02OY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX, SQRT
C     .. Executable Statements ..
C
      IWARN  = 0
      INFO   = 0
      DISCR  = LSAME( DICO, 'D' )
      LJOBB  = LSAME( JOBB, 'B' )
      LFACN  = LSAME( FACT, 'N' )
      LFACQ  = LSAME( FACT, 'C' )
      LFACR  = LSAME( FACT, 'D' )
      LFACB  = LSAME( FACT, 'B' )
      LUPLO  = LSAME( UPLO, 'U' )
      LSORT  = LSAME( SORT, 'S' )
      REFINE = LSAME( ACC,  'R' )
      NN = 2*N
      IF ( LJOBB ) THEN
         LJOBL  = LSAME( JOBL, 'Z' )
         LJOBLN = LSAME( JOBL, 'N' )
         LSCAL  = LSAME( SCAL, 'G' )
         NNM = NN + M
         LDW = MAX( NNM, 3*M )
      ELSE
         LSCAL = .FALSE.
         NNM = NN
         LDW = 1
      END IF
      NP1 = N + 1
C
C     Test the input scalar arguments.
C
      IF( .NOT.DISCR .AND. .NOT.LSAME( DICO, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LJOBB .AND. .NOT.LSAME( JOBB, 'G' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LFACQ .AND. .NOT.LFACR .AND. .NOT.LFACB
     $                                     .AND. .NOT.LFACN ) THEN
         INFO = -3
      ELSE IF( .NOT.LJOBB .OR. LFACN ) THEN
         IF( .NOT.LUPLO .AND. .NOT.LSAME( UPLO, 'L' ) )
     $      INFO = -4
      END IF
      IF( INFO.EQ.0 .AND. LJOBB ) THEN
         IF( .NOT.LJOBL .AND. .NOT.LJOBLN ) THEN
            INFO = -5
         ELSE IF( .NOT.LSCAL .AND. .NOT. LSAME( SCAL, 'N' ) ) THEN
            INFO = -6
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.LSORT .AND. .NOT.LSAME( SORT, 'U' ) ) THEN
            INFO = -7
         ELSE IF( .NOT.REFINE .AND. .NOT.LSAME( ACC, 'N' ) ) THEN
            INFO = -8
         ELSE IF( N.LT.0 ) THEN
            INFO = -9
         ELSE IF( LJOBB ) THEN
            IF( M.LT.0 )
     $         INFO = -10
         END IF
      END IF
      IF( INFO.EQ.0 .AND. .NOT.LFACN ) THEN
         IF( P.LT.0 )
     $      INFO = -11
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -13
         ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
            INFO = -15
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -17
         ELSE IF( ( ( LFACN.OR.LFACR ) .AND. LDQ.LT.MAX( 1, N ) ) .OR.
     $            ( ( LFACQ.OR.LFACB ) .AND. LDQ.LT.MAX( 1, P ) ) ) THEN
            INFO = -19
         ELSE IF( LJOBB ) THEN
            IF ( ( LFACN.OR.LFACQ ) .AND. LDR.LT.MAX( 1, M ) .OR.
     $           ( LFACR.OR.LFACB ) .AND. LDR.LT.MAX( 1, P ) ) THEN
               INFO = -21
            ELSE IF( ( LJOBLN .AND. LDL.LT.MAX( 1, N ) ) .OR.
     $               ( LJOBL  .AND. LDL.LT.1 ) ) THEN
               INFO = -23
            END IF
         ELSE
            IF( LDR.LT.1 ) THEN
               INFO = -21
            ELSE IF( LDL.LT.1 ) THEN
               INFO = -23
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDX.LT.MAX( 1, N ) ) THEN
            INFO = -26
         ELSE IF( LDS.LT.MAX( 1, NNM ) ) THEN
            INFO = -31
         ELSE IF( LDT.LT.MAX( 1, NNM ) ) THEN
            INFO = -33
         ELSE IF( LDU.LT.MAX( 1, NN ) ) THEN
            INFO = -35
         ELSE IF( LDWORK.LT.MAX( 14*N + 23, 16*N, LDW ) ) THEN
            INFO = -39
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SG02AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = FOUR
         DWORK(4) = ONE
         RETURN
      END IF
C
C     Start computations.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      LSCAL = LSCAL .AND. LJOBB
      IF ( LSCAL ) THEN
C
C        Scale the matrices Q, R (or G), and L so that
C           norm(Q) + norm(R) + norm(L) = 1,
C        using the 1-norm. If Q and/or R are factored, the norms of
C        the factors are used.
C        Workspace: need   max(N,M), if FACT = 'N';
C                          N,        if FACT = 'D';
C                          M,        if FACT = 'C'.
C
         IF ( LFACN .OR. LFACR ) THEN
            SCALE = DLANSY( '1-norm', UPLO, N, Q, LDQ, DWORK )
            QTYPE = UPLO
            NP = N
         ELSE
            SCALE = DLANGE( '1-norm', P, N, Q, LDQ, DWORK )
            QTYPE = 'G'
            NP = P
         END IF
C
         IF ( LFACN .OR. LFACQ ) THEN
            RNORM = DLANSY( '1-norm', UPLO, M, R, LDR, DWORK )
            RTYPE = UPLO
            MP = M
         ELSE
            RNORM = DLANGE( '1-norm', P, M, R, LDR, DWORK )
            RTYPE = 'G'
            MP = P
         END IF
         SCALE = SCALE + RNORM
C
         IF ( LJOBLN )
     $      SCALE = SCALE + DLANGE( '1-norm', N, M, L, LDL, DWORK )
         IF ( SCALE.EQ.ZERO )
     $      SCALE = ONE
C
         CALL DLASCL( QTYPE, 0, 0, SCALE, ONE, NP, N, Q, LDQ, INFO1 )
         CALL DLASCL( RTYPE, 0, 0, SCALE, ONE, MP, M, R, LDR, INFO1 )
         IF ( LJOBLN )
     $      CALL DLASCL( 'G', 0, 0, SCALE, ONE, N, M, L, LDL, INFO1 )
      ELSE
         SCALE = ONE
      END IF
C
C     Construct the extended matrix pair.
C     Workspace: need   1,                if JOBB = 'G',
C                       max(1,2*N+M,3*M), if JOBB = 'B';
C                prefer larger.
C
      CALL SB02OY( 'Optimal control', DICO, JOBB, FACT, UPLO, JOBL,
     $             'Not identity E', N, M, P, A, LDA, B, LDB, Q, LDQ, R,
     $             LDR, L, LDL, E, LDE, S, LDS, T, LDT, TOL, IWORK,
     $             DWORK, LDWORK, INFO )
C
      IF ( LSCAL ) THEN
C
C        Undo scaling of the data arrays.
C
         CALL DLASCL( QTYPE, 0, 0, ONE, SCALE, NP, N, Q, LDQ, INFO1 )
         CALL DLASCL( RTYPE, 0, 0, ONE, SCALE, MP, M, R, LDR, INFO1 )
         IF ( LJOBLN )
     $      CALL DLASCL( 'G', 0, 0, ONE, SCALE, N, M, L, LDL, INFO1 )
      END IF
C
      IF ( INFO.NE.0 )
     $   RETURN
      WRKOPT = DWORK(1)
      IF ( LJOBB )
     $   RCONDL = DWORK(2)
C
C     Workspace: need   max(7*(2*N+1)+16,16*N);
C                prefer larger.
C
      IF ( DISCR ) THEN
         IF ( LSORT ) THEN
C
C           The natural tendency of the QZ algorithm to get the largest
C           eigenvalues in the leading part of the matrix pair is
C           exploited, by computing the unstable eigenvalues of the
C           permuted matrix pair.
C
            CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OV, NN, T,
     $                  LDT, S, LDS, NDIM, ALFAR, ALFAI, BETA, U, LDU,
     $                  U, LDU, DWORK, LDWORK, BWORK, INFO1 )
            CALL DSWAP( N, ALFAR(NP1), 1, ALFAR, 1 )
            CALL DSWAP( N, ALFAI(NP1), 1, ALFAI, 1 )
            CALL DSWAP( N, BETA (NP1), 1, BETA , 1 )
         ELSE
            CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OV, NN, S,
     $                  LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U, LDU,
     $                  U, LDU, DWORK, LDWORK, BWORK, INFO1 )
         END IF
      ELSE
         IF ( LSORT ) THEN
            CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OW, NN, S,
     $                  LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U, LDU,
     $                  U, LDU, DWORK, LDWORK, BWORK, INFO1 )
         ELSE
            CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OU, NN, S,
     $                  LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U, LDU,
     $                  U, LDU, DWORK, LDWORK, BWORK, INFO1 )
         END IF
      END IF
      IF ( INFO1.GT.0 .AND. INFO1.LE.NN+1 ) THEN
         INFO = 2
      ELSE IF ( INFO1.EQ.NN+2 ) THEN
         INFO = 4
      ELSE IF ( INFO1.EQ.NN+3 ) THEN
         INFO = 3
      ELSE IF ( NDIM.NE.N ) THEN
         INFO = 5
      END IF
      IF ( INFO.NE.0 )
     $   RETURN
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
C
C     Take the non-identity matrix E into account and orthogonalize the
C     basis. Use the array X as workspace.
C     Workspace: need   N;
C                prefer N*NB.
C
      CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE, E, LDE,
     $            U, LDU, ZERO, X, LDX )
      CALL DLACPY( 'Full', N, N, X, LDX, U, LDU )
      CALL DGEQRF( NN, N, U, LDU, X, DWORK, LDWORK, INFO1 )
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      CALL DORGQR( NN, N, N, U, LDU, X, DWORK, LDWORK, INFO1 )
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
C
C     Check for the symmetry of the solution. The array X is again used
C     as workspace.
C
      CALL DGEMM( 'Transpose', 'No transpose', N, N, N, ONE, U, LDU,
     $            U(NP1,1), LDU, ZERO, X, LDX )
      U12M = ZERO
      ASYM = ZERO
C
      DO 20 J = 1, N
C
         DO 10 I = 1, N
            U12M = MAX( U12M, ABS( X(I,J) ) )
            ASYM = MAX( ASYM, ABS( X(I,J) - X(J,I) ) )
   10    CONTINUE
C
   20 CONTINUE
C
      EPS  = DLAMCH( 'Epsilon' )
      SEPS = SQRT( EPS )
      ASYM = ASYM - SEPS
      IF ( ASYM.GT.P1*U12M ) THEN
         INFO = 6
         RETURN
      ELSE IF ( ASYM.GT.SEPS ) THEN
         IWARN = 1
      END IF
C
C     Compute the solution of X*U(1,1) = U(2,1). Use the (2,1) block
C     of S as a workspace for factoring U(1,1).
C
      IF ( REFINE ) THEN
C
C        Use LU factorization and iterative refinement for finding X.
C        Workspace:  need   8*N.
C
C        First transpose U(2,1) in-situ.
C
         DO 30 I = 1, N - 1
            CALL DSWAP( N-I, U(N+I,I+1), LDU, U(N+I+1,I), 1 )
   30    CONTINUE
C
         IWR = 1
         IWC = IWR + N
         IWF = IWC + N
         IWB = IWF + N
         IW  = IWB + N
C
         CALL MB02PD( 'Equilibrate', 'Transpose', N, N, U, LDU,
     $                S(NP1,1), LDS, IWORK, EQUED, DWORK(IWR),
     $                DWORK(IWC), U(NP1,1), LDU, X, LDX, RCONDU,
     $                DWORK(IWF), DWORK(IWB), IWORK(NP1), DWORK(IW),
     $                INFO1 )
C
C        Transpose U(2,1) back in-situ.
C
         DO 40 I = 1, N - 1
            CALL DSWAP( N-I, U(N+I,I+1), LDU, U(N+I+1,I), 1 )
   40    CONTINUE
C
         IF( .NOT.LSAME( EQUED, 'N' ) ) THEN
C
C           Undo the equilibration of U(1,1) and U(2,1).
C
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
C
            IF( ROWEQU ) THEN
C
               DO 50 I = 0, N - 1
                  DWORK(IWR+I) = ONE / DWORK(IWR+I)
   50          CONTINUE
C
               CALL MB01SD( 'Row scaling', N, N, U, LDU, DWORK(IWR),
     $                      DWORK(IWC) )
            END IF
C
            IF( COLEQU ) THEN
C
               DO 60 I = 0, N - 1
                  DWORK(IWC+I) = ONE / DWORK(IWC+I)
   60          CONTINUE
C
               CALL MB01SD( 'Column scaling', NN, N, U, LDU, DWORK(IWR),
     $                      DWORK(IWC) )
             END IF
         END IF
C
         PIVOTU = DWORK(IW)
C
         IF ( INFO1.GT.0 ) THEN
C
C           Singular matrix. Set INFO and DWORK for error return.
C
            INFO = 7
            GO TO 80
         END IF
C
      ELSE
C
C        Use LU factorization and a standard solution algorithm.
C
         CALL DLACPY( 'Full', N, N, U, LDU, S(NP1,1), LDS )
         CALL DLACPY( 'Full', N, N, U(NP1,1), LDU, X, LDX )
C
C        Solve the system X*U(1,1) = U(2,1).
C
         CALL MB02VD( 'No Transpose', N, N, S(NP1,1), LDS, IWORK, X,
     $                LDX, INFO1 )
C
         IF ( INFO1.NE.0 ) THEN
            INFO   = 7
            RCONDU = ZERO
            GO TO 80
         ELSE
C
C           Compute the norm of U(1,1).
C
            UNORM = DLANGE( '1-norm', N, N, U, LDU, DWORK )
C
C           Estimate the reciprocal condition of U(1,1).
C           Workspace: need 4*N.
C
            CALL DGECON( '1-norm', N, S(NP1,1), LDS, UNORM, RCONDU,
     $                   DWORK, IWORK(NP1), INFO )
C
            IF ( RCONDU.LT.EPS ) THEN
C
C              Nearly singular matrix. Set IWARN for warning indication.
C
               IWARN = 1
            END IF
            WRKOPT = MAX( WRKOPT, 4*N )
         END IF
      END IF
C
C     Set S(2,1) to zero.
C
      CALL DLASET( 'Full', N, N, ZERO, ZERO, S(NP1,1), LDS )
C
C     Make sure the solution matrix X is symmetric.
C
      DO 70 I = 1, N - 1
         CALL DAXPY( N-I, ONE, X(I,I+1), LDX, X(I+1,I), 1 )
         CALL DSCAL( N-I, HALF, X(I+1,I), 1 )
         CALL DCOPY( N-I, X(I+1,I), 1, X(I,I+1), LDX )
   70 CONTINUE
C
      IF ( LSCAL ) THEN
C
C        Undo scaling for the solution X.
C
         CALL DLASCL( 'G', 0, 0, ONE, SCALE, N, N, X, LDX, INFO1 )
      END IF
C
      DWORK(1) = WRKOPT
C
   80 CONTINUE
      IF ( LJOBB )
     $   DWORK(2) = RCONDL
      IF ( REFINE )
     $   DWORK(3) = PIVOTU
      DWORK(4) = SCALE
C
      RETURN
C *** Last line of SG02AD ***
      END
