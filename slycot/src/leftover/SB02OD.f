      SUBROUTINE SB02OD( DICO, JOBB, FACT, UPLO, JOBL, SORT, N, M, P, A,
     $                   LDA, B, LDB, Q, LDQ, R, LDR, L, LDL, RCOND, X,
     $                   LDX, ALFAR, ALFAI, BETA, S, LDS, T, LDT, U,
     $                   LDU, TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
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
C                              -1
C        Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1)
C
C     or the discrete-time algebraic Riccati equation
C                                     -1
C        X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2)
C
C     where A, B, Q, R, and L are N-by-N, N-by-M, N-by-N, M-by-M and
C     N-by-M matrices, respectively, such that Q = C'C, R = D'D and
C     L = C'D; X is an N-by-N symmetric matrix.
C     The routine also returns the computed values of the closed-loop
C     spectrum of the system, i.e., the stable eigenvalues lambda(1),
C     ..., lambda(N) of the corresponding Hamiltonian or symplectic
C     pencil, in the continuous-time case or discrete-time case,
C     respectively.
C                              -1
C     Optionally, matrix G = BR  B' may be given instead of B and R.
C     Other options include the case with Q and/or R given in a
C     factored form, Q = C'C, R = D'D, and with L a zero matrix.
C
C     The routine uses the method of deflating subspaces, based on
C     reordering the eigenvalues in a generalized Schur matrix pair.
C     A standard eigenproblem is solved in the continuous-time case
C     if G is given.
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
C             the matrices G and Q (if FACT = 'N'), or Q and R (if
C             JOBB = 'B'), is stored, as follows:
C             = 'U':  Upper triangle is stored;
C             = 'L':  Lower triangle is stored.
C
C     JOBL    CHARACTER*1
C             Specifies whether or not the matrix L is zero, as follows:
C             = 'Z':  L is zero;
C             = 'N':  L is nonzero.
C             JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed.
C             SLICOT Library routine SB02MT should be called just before
C             SB02OD, for obtaining the results when JOBB = 'G' and
C             JOBL = 'N'.
C
C     SORT    CHARACTER*1
C             Specifies which eigenvalues should be obtained in the top
C             of the generalized Schur form, as follows:
C             = 'S':  Stable   eigenvalues come first;
C             = 'U':  Unstable eigenvalues come first.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e. the order of the matrices
C             A, Q, and X, and the number of rows of the matrices B
C             and L.  N >= 0.
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
C             state matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
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
C             If JOBB = 'B', the triangular part of this array defined
C             by UPLO is modified internally, but is restored on exit.
C             If FACT = 'C' or 'B', the leading P-by-N part of this
C             array must contain the output matrix C of the system.
C             If JOBB = 'B', this part is modified internally, but is
C             restored on exit.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= MAX(1,N) if FACT = 'N' or 'D',
C             LDQ >= MAX(1,P) if FACT = 'C' or 'B'.
C
C     R       (input) DOUBLE PRECISION array, dimension (LDR,M)
C             If FACT = 'N' or 'C', the leading M-by-M upper triangular
C             part (if UPLO = 'U') or lower triangular part (if UPLO =
C             'L') of this array must contain the upper triangular part
C             or lower triangular part, respectively, of the symmetric
C             input weighting matrix R. The stricly lower triangular
C             part (if UPLO = 'U') or stricly upper triangular part (if
C             UPLO = 'L') is not referenced.
C             The triangular part of this array defined by UPLO is
C             modified internally, but is restored on exit.
C             If FACT = 'D' or 'B', the leading P-by-M part of this
C             array must contain the direct transmission matrix D of the
C             system. This part is modified internally, but is restored
C             on exit.
C             If JOBB = 'G', this array is not referenced.
C
C     LDR     INTEGER
C             The leading dimension of array R.
C             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C';
C             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B';
C             LDR >= 1        if JOBB = 'G'.
C
C     L       (input) DOUBLE PRECISION array, dimension (LDL,M)
C             If JOBL = 'N' (and JOBB = 'B'), the leading N-by-M part of
C             this array must contain the cross weighting matrix L.
C             This part is modified internally, but is restored on exit.
C             If JOBL = 'Z' or JOBB = 'G', this array is not referenced.
C
C     LDL     INTEGER
C             The leading dimension of array L.
C             LDL >= MAX(1,N) if JOBL = 'N' and JOBB = 'B';
C             LDL >= 1        if JOBL = 'Z' or  JOBB = 'G'.
C
C     RCOND   (output) DOUBLE PRECISION
C             An estimate of the reciprocal of the condition number (in
C             the 1-norm) of the N-th order system of algebraic
C             equations from which the solution matrix X is obtained.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N)
C             The leading N-by-N part of this array contains the
C             solution matrix X of the problem.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= MAX(1,N).
C
C     ALFAR   (output) DOUBLE PRECISION array, dimension (2*N)
C     ALFAI   (output) DOUBLE PRECISION array, dimension (2*N)
C     BETA    (output) DOUBLE PRECISION array, dimension (2*N)
C             The generalized eigenvalues of the 2N-by-2N matrix pair,
C             ordered as specified by SORT (if INFO = 0). For instance,
C             if SORT = 'S', the leading N elements of these arrays
C             contain the closed-loop spectrum of the system matrix
C             A - BF, where F is the optimal feedback matrix computed
C             based on the solution matrix X. Specifically,
C                lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for
C             k = 1,2,...,N.
C             If DICO = 'C' and JOBB = 'G', the elements of BETA are
C             set to 1.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,*)
C             The leading 2N-by-2N part of this array contains the
C             ordered real Schur form S of the first matrix in the
C             reduced matrix pencil associated to the optimal problem,
C             or of the corresponding Hamiltonian matrix, if DICO = 'C'
C             and JOBB = 'G'. That is,
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
C             LDS >= MAX(1,2*N+M) if JOBB = 'B',
C             LDS >= MAX(1,2*N)   if JOBB = 'G'.
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,2*N)
C             If DICO = 'D' or JOBB = 'B', the leading 2N-by-2N part of
C             this array contains the ordered upper triangular form T of
C             the second matrix in the reduced matrix pencil associated
C             to the optimal problem. That is,
C
C                    (T   T  )
C                    ( 11  12)
C                T = (       ),
C                    (0   T  )
C                    (     22)
C
C             where T  , T   and T   are N-by-N matrices.
C                    11   12      22
C             If DICO = 'C' and JOBB = 'G' this array is not referenced.
C
C     LDT     INTEGER
C             The leading dimension of array T.
C             LDT >= MAX(1,2*N+M) if JOBB = 'B',
C             LDT >= MAX(1,2*N)   if JOBB = 'G' and DICO = 'D',
C             LDT >= 1            if JOBB = 'G' and DICO = 'C'.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N)
C             The leading 2N-by-2N part of this array contains the right
C             transformation matrix U which reduces the 2N-by-2N matrix
C             pencil to the ordered generalized real Schur form (S,T),
C             or the Hamiltonian matrix to the ordered real Schur
C             form S, if DICO = 'C' and JOBB = 'G'. That is,
C
C                    (U   U  )
C                    ( 11  12)
C                U = (       ),
C                    (U   U  )
C                    ( 21  22)
C
C             where U  , U  , U   and U   are N-by-N matrices.
C                    11   12   21      22
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,2*N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used to test for near singularity of
C             the original matrix pencil, specifically of the triangular
C             factor obtained during the reduction process. If the user
C             sets TOL > 0, then the given value of TOL is used as a
C             lower bound for the reciprocal condition number of that
C             matrix; a matrix whose estimated condition number is less
C             than 1/TOL is considered to be nonsingular. If the user
C             sets TOL <= 0, then a default tolerance, defined by
C             TOLDEF = EPS, is used instead, where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH).
C             This parameter is not referenced if JOBB = 'G'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK >= MAX(1,M,2*N) if JOBB = 'B',
C             LIWORK >= MAX(1,2*N)   if JOBB = 'G'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK. If JOBB = 'B' and N > 0, DWORK(2) returns the
C             reciprocal of the condition number of the M-by-M lower
C             triangular matrix obtained after compressing the matrix
C             pencil of order 2N+M to obtain a pencil of order 2N.
C             If INFO = 0 or INFO = 6, DWORK(3) returns the scaling
C             factor used internally, which should multiply the
C             submatrix Y2 to recover X from the first N columns of U
C             (see METHOD).
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(3,6*N),                       if JOBB = 'G',
C                                                            DICO = 'C';
C             LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G',
C                                                            DICO = 'D';
C             LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'.
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the computed extended matrix pencil is singular,
C                   possibly due to rounding errors;
C             = 2:  if the QZ (or QR) algorithm failed;
C             = 3:  if reordering of the (generalized) eigenvalues
C                   failed;
C             = 4:  if after reordering, roundoff changed values of
C                   some complex eigenvalues so that leading eigenvalues
C                   in the (generalized) Schur form no longer satisfy
C                   the stability condition; this could also be caused
C                   due to scaling;
C             = 5:  if the computed dimension of the solution does not
C                   equal N;
C             = 6:  if a singular matrix was encountered during the
C                   computation of the solution matrix X.
C
C     METHOD
C
C     The routine uses a variant of the method of deflating subspaces
C     proposed by van Dooren [1]. See also [2], [3].
C     It is assumed that (A,B) is stabilizable and (C,A) is detectable.
C     Under these assumptions the algebraic Riccati equation is known to
C     have a unique non-negative definite solution.
C     The first step in the method of deflating subspaces is to form the
C     extended Hamiltonian matrices, dimension 2N + M given by
C
C           discrete-time                   continuous-time
C
C     |A   0   B|     |I   0   0|    |A   0   B|     |I   0   0|
C     |Q  -I   L| - z |0  -A'  0|,   |Q   A'  L| - s |0  -I   0|.
C     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0|
C
C     Next, these pencils are compressed to a form (see [1])
C
C        lambda x A  - B .
C                  f    f
C
C     This generalized eigenvalue problem is then solved using the QZ
C     algorithm and the stable deflating subspace Ys is determined.
C     If [Y1'|Y2']' is a basis for Ys, then the required solution is
C                       -1
C            X = Y2 x Y1  .
C     A standard eigenvalue problem is solved using the QR algorithm in
C     the continuous-time case when G is given (DICO = 'C', JOBB = 'G').
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         A Generalized Eigenvalue Approach for Solving Riccati
C         Equations.
C         SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981.
C
C     [2] Mehrmann, V.
C         The Autonomous Linear Quadratic Control Problem. Theory and
C         Numerical Solution.
C         Lect. Notes in Control and Information Sciences, vol. 163,
C         Springer-Verlag, Berlin, 1991.
C
C     [3] Sima, V.
C         Algorithms for Linear-Quadratic Optimization.
C         Pure and Applied Mathematics: A Series of Monographs and
C         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996.
C
C     NUMERICAL ASPECTS
C
C     This routine is particularly suited for systems where the matrix R
C     is ill-conditioned. Internal scaling is used.
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
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997.
C     Supersedes Release 2.0 routine SB02CD by T.G.J. Beelen, Philips,
C     Eindhoven, Holland.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, June 2002,
C     December 2002, January 2005.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, closed loop system, continuous-time
C     system, discrete-time system, optimal regulator, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, THREE
      PARAMETER         ( ZERO  = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                    THREE = 3.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOBB, JOBL, SORT, UPLO
      INTEGER           INFO, LDA, LDB, LDL, LDQ, LDR, LDS, LDT, LDU,
     $                  LDWORK, LDX, M, N, P
      DOUBLE PRECISION  RCOND, TOL
C     .. Array Arguments ..
      LOGICAL           BWORK(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), ALFAI(*), ALFAR(*), B(LDB,*), BETA(*),
     $                  DWORK(*), L(LDL,*), Q(LDQ,*), R(LDR,*),
     $                  S(LDS,*), T(LDT,*), U(LDU,*), X(LDX,*)
C     .. Local Scalars ..
      CHARACTER         QTYPE, RTYPE
      LOGICAL           DISCR, LFACB, LFACN, LFACQ, LFACR, LJOBB, LJOBL,
     $                  LJOBLN, LSCAL, LSCL, LSORT, LUPLO
      INTEGER           I, INFO1, J, LDW, MP, NDIM, NN, NNM, NP, NP1,
     $                  WRKOPT
      DOUBLE PRECISION  QSCAL, RCONDL, RNORM, RSCAL, SCALE, UNORM
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      LOGICAL           LSAME, SB02MR, SB02MV, SB02OU, SB02OV, SB02OW
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANSY
      EXTERNAL          DLAMCH, DLANGE, DLANSY, LSAME, SB02MR, SB02MV,
     $                  SB02OU, SB02OV, SB02OW
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGECON, DGEES, DGETRF, DGETRS,
     $                  DGGES, DLACPY, DLASCL, DLASET, DSCAL, DSWAP,
     $                  SB02OY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO  = 0
      DISCR = LSAME( DICO, 'D' )
      LJOBB = LSAME( JOBB, 'B' )
      LFACN = LSAME( FACT, 'N' )
      LFACQ = LSAME( FACT, 'C' )
      LFACR = LSAME( FACT, 'D' )
      LFACB = LSAME( FACT, 'B' )
      LUPLO = LSAME( UPLO, 'U' )
      LSORT = LSAME( SORT, 'S' )
C
      NN = 2*N
      IF ( LJOBB ) THEN
         LJOBL  = LSAME( JOBL, 'Z' )
         LJOBLN = LSAME( JOBL, 'N' )
         NNM = NN + M
         LDW = MAX( NNM, 3*M )
      ELSE
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
         IF( .NOT.LJOBL .AND. .NOT.LJOBLN )
     $      INFO = -5
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( .NOT.LSORT .AND. .NOT.LSAME( SORT, 'U' ) ) THEN
            INFO = -6
         ELSE IF( N.LT.0 ) THEN
            INFO = -7
         ELSE IF( LJOBB ) THEN
            IF( M.LT.0 )
     $         INFO = -8
         END IF
      END IF
      IF( INFO.EQ.0 .AND. .NOT.LFACN ) THEN
         IF( P.LT.0 )
     $      INFO = -9
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -11
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -13
         ELSE IF( ( ( LFACN.OR.LFACR ) .AND. LDQ.LT.MAX( 1, N ) ) .OR.
     $            ( ( LFACQ.OR.LFACB ) .AND. LDQ.LT.MAX( 1, P ) ) ) THEN
            INFO = -15
         ELSE IF( LDR.LT.1 ) THEN
            INFO = -17
         ELSE IF( LDL.LT.1 ) THEN
            INFO = -19
         ELSE IF( LJOBB ) THEN
            IF ( ( LFACN.OR.LFACQ ) .AND. LDR.LT.M .OR.
     $           ( LFACR.OR.LFACB ) .AND. LDR.LT.P ) THEN
               INFO = -17
            ELSE IF( LJOBLN .AND. LDL.LT.N ) THEN
               INFO = -19
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDX.LT.MAX( 1, N ) ) THEN
            INFO = -22
         ELSE IF( LDS.LT.MAX( 1, NNM ) ) THEN
            INFO = -27
         ELSE IF( LDT.LT.1 ) THEN
            INFO = -29
         ELSE IF( LDU.LT.MAX( 1, NN ) ) THEN
            INFO = -31
         ELSE IF( LDWORK.LT.MAX( 3, 6*N ) ) THEN
            INFO = -35
         ELSE IF( DISCR .OR. LJOBB ) THEN
            IF( LDT.LT.NNM ) THEN
               INFO = -29
            ELSE IF( LDWORK.LT.MAX( 14*N + 23, 16*N, LDW ) ) THEN
               INFO = -35
            END IF
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB02OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         RCOND = ONE
         DWORK(1) = THREE
         DWORK(3) = ONE
         RETURN
      END IF
C
C     Always scale the matrix pencil.
C
      LSCAL = .TRUE.
C
C     Start computations.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( LSCAL .AND. LJOBB ) THEN
C
C        Scale the matrices Q, R, and L so that
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
         IF ( LFACN .OR. LFACR ) THEN
            QSCAL = SCALE
         ELSE
            QSCAL = SQRT( SCALE )
         END IF
C
         IF ( LFACN .OR. LFACQ ) THEN
            RSCAL = SCALE
         ELSE
            RSCAL = SQRT( SCALE )
         END IF
C
         CALL DLASCL( QTYPE, 0, 0, QSCAL, ONE, NP, N, Q, LDQ, INFO1 )
         CALL DLASCL( RTYPE, 0, 0, RSCAL, ONE, MP, M, R, LDR, INFO1 )
         IF ( LJOBLN )
     $      CALL DLASCL( 'G', 0, 0, SCALE, ONE, N, M, L, LDL, INFO1 )
      END IF
C
C     Construct the extended matrix pair.
C
C     Workspace: need   1,                if JOBB = 'G',
C                       max(1,2*N+M,3*M), if JOBB = 'B';
C                prefer larger.
C
      CALL SB02OY( 'Optimal control', DICO, JOBB, FACT, UPLO, JOBL,
     $             'Identity E', N, M, P, A, LDA, B, LDB, Q, LDQ, R,
     $             LDR, L, LDL, U, 1, S, LDS, T, LDT, TOL, IWORK, DWORK,
     $             LDWORK, INFO )
C
      IF ( LSCAL .AND. LJOBB ) THEN
C
C        Undo scaling of the data arrays.
C
         CALL DLASCL( QTYPE, 0, 0, ONE, QSCAL, NP, N, Q, LDQ, INFO1 )
         CALL DLASCL( RTYPE, 0, 0, ONE, RSCAL, MP, M, R, LDR, INFO1 )
         IF ( LJOBLN )
     $      CALL DLASCL( 'G', 0, 0, ONE, SCALE, N, M, L, LDL, INFO1 )
      END IF
C
      IF ( INFO.NE.0 )
     $   RETURN
      WRKOPT = DWORK(1)
      IF ( LJOBB ) RCONDL = DWORK(2)
C
      IF ( LSCAL .AND. .NOT.LJOBB ) THEN
C
C        This part of the code is used when G is given (JOBB = 'G').
C        A standard eigenproblem is solved in the continuous-time case.
C        Scale the Hamiltonian matrix S, if DICO = 'C', or the
C        symplectic pencil (S,T), if DICO = 'D', using the square roots
C        of the norms of the matrices Q and G.
C        Workspace: need   N.
C
         IF ( LFACN .OR. LFACR ) THEN
            SCALE = SQRT( DLANSY( '1-norm', UPLO, N, Q, LDQ, DWORK ) )
         ELSE
            SCALE = DLANGE( '1-norm', P, N, Q, LDQ, DWORK )
         END IF
         RNORM = SQRT( DLANSY( '1-norm', UPLO, N, B, LDB, DWORK ) )
C
         LSCL = MIN( SCALE, RNORM ).GT.ZERO .AND. SCALE.NE.RNORM
C
         IF( LSCL ) THEN
            IF( DISCR ) THEN
               CALL DLASCL( 'G', 0, 0, SCALE, RNORM, N, N, S(NP1,1),
     $                      LDS, INFO1 )
               CALL DLASCL( 'G', 0, 0, RNORM, SCALE, N, N, T(1,NP1),
     $                      LDT, INFO1 )
            ELSE
               CALL DLASCL( 'G', 0, 0, SCALE, -RNORM, N, N, S(NP1,1),
     $                      LDS, INFO1 )
               CALL DLASCL( 'G', 0, 0, RNORM, SCALE, N, N, S(1,NP1),
     $                      LDS, INFO1 )
               CALL DLASCL( 'G', 0, 0, ONE, -ONE, N, N, S(NP1,NP1),
     $                      LDS, INFO1 )
            END IF
         ELSE
            IF( .NOT.DISCR ) THEN
               CALL DLASCL( 'G', 0, 0, ONE, -ONE, N, NN, S(NP1,1), LDS,
     $                      INFO1 )
            END IF
         END IF
      ELSE
         LSCL = .FALSE.
      END IF
C
C     Workspace: need   max(7*(2*N+1)+16,16*N),
C                                          if JOBB = 'B' or  DICO = 'D';
C                       6*N,               if JOBB = 'G' and DICO = 'C';
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
            CALL DSWAP( N, BETA (NP1), 1, BETA,  1 )
         ELSE
            CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OV, NN, S,
     $                  LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U, LDU,
     $                  U, LDU, DWORK, LDWORK, BWORK, INFO1 )
         END IF
      ELSE
         IF ( LJOBB ) THEN
            IF ( LSORT ) THEN
               CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OW, NN,
     $                     S, LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U,
     $                     LDU, U, LDU, DWORK, LDWORK, BWORK, INFO1 )
            ELSE
               CALL DGGES( 'No vectors', 'Vectors', 'Sort', SB02OU, NN,
     $                     S, LDS, T, LDT, NDIM, ALFAR, ALFAI, BETA, U,
     $                     LDU, U, LDU, DWORK, LDWORK, BWORK, INFO1 )
            END IF
         ELSE
            IF ( LSORT ) THEN
               CALL DGEES( 'Vectors', 'Sort', SB02MV, NN, S, LDS, NDIM,
     $                     ALFAR, ALFAI, U, LDU, DWORK, LDWORK, BWORK,
     $                     INFO1 )
            ELSE
               CALL DGEES( 'Vectors', 'Sort', SB02MR, NN, S, LDS, NDIM,
     $                     ALFAR, ALFAI, U, LDU, DWORK, LDWORK, BWORK,
     $                     INFO1 )
            END IF
            DUM(1) = ONE
            CALL DCOPY( NN, DUM, 0, BETA, 1 )
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
C     Select submatrices U1 and U2 out of the array U which define the
C     solution X = U2 x inv(U1).
C     Since X = X' we may obtain X as the solution of the system of
C     linear equations U1' x X = U2', where
C        U1 = U(1:n, 1:n),
C        U2 = U(n+1:2n, 1:n).
C     Use the (2,1) block of S as a workspace for factoring U1.
C
      DO 20 J = 1, N
         CALL DCOPY( N, U(NP1,J), 1, X(J,1), LDX )
   20 CONTINUE
C
      CALL DLACPY( 'Full', N, N, U, LDU, S(NP1,1), LDS )
C
C     Check if U1 is singular.
C
      UNORM = DLANGE( '1-norm', N, N, S(NP1,1), LDS, DWORK )
C
C     Solve the system U1' x X = U2'.
C
      CALL DGETRF( N, N, S(NP1,1), LDS, IWORK, INFO1 )
      IF ( INFO1.NE.0 ) THEN
         INFO = 6
         DWORK(3) = ONE
         IF ( LSCAL ) THEN
            IF ( LJOBB ) THEN
               DWORK(3) = SCALE
            ELSE IF ( LSCL ) THEN
               DWORK(3) = SCALE / RNORM
            END IF
         END IF
         RETURN
      ELSE
C
C        Estimate the reciprocal condition of U1.
C        Workspace: need 3*N.
C
         CALL DGECON( '1-norm', N, S(NP1,1), LDS, UNORM, RCOND, DWORK,
     $                IWORK(NP1), INFO )
C
         IF ( RCOND.LT.DLAMCH( 'Epsilon' ) ) THEN
C
C           Nearly singular matrix.  Set INFO for error return.
C
            INFO = 6
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, 3*N )
         CALL DGETRS( 'Transpose', N, N, S(NP1,1), LDS, IWORK, X, LDX,
     $                INFO1 )
C
C        Set S(2,1) to zero.
C
         CALL DLASET( 'Full', N, N, ZERO, ZERO, S(NP1,1), LDS )
C
         IF ( LSCAL ) THEN
C
C           Prepare to undo scaling for the solution X.
C
            IF ( .NOT.LJOBB ) THEN
               IF ( LSCL ) THEN
                  SCALE = SCALE / RNORM
               ELSE
                  SCALE = ONE
               END IF
            END IF
            DWORK(3) = SCALE
            SCALE = HALF*SCALE
         ELSE
            DWORK(3) = ONE
            SCALE = HALF
         END IF
C
C        Make sure the solution matrix X is symmetric.
C
         DO 40 I = 1, N
            CALL DAXPY( N-I+1, ONE, X(I,I), LDX, X(I,I), 1 )
            CALL DSCAL( N-I+1, SCALE, X(I,I), 1 )
            CALL DCOPY( N-I+1, X(I,I), 1, X(I,I), LDX )
   40    CONTINUE
      END IF
C
      DWORK(1) = WRKOPT
      IF ( LJOBB ) DWORK(2) = RCONDL
C
      RETURN
C *** Last line of SB02OD ***
      END
