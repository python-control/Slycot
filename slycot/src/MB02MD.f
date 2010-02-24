      SUBROUTINE MB02MD( JOB, M, N, L, RANK, C, LDC, S, X, LDX, TOL,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     To solve the Total Least Squares (TLS) problem using a Singular
C     Value Decomposition (SVD) approach.
C     The TLS problem assumes an overdetermined set of linear equations
C     AX = B, where both the data matrix A as well as the observation
C     matrix B are inaccurate. The routine also solves determined and
C     underdetermined sets of equations by computing the minimum norm
C     solution.
C     It is assumed that all preprocessing measures (scaling, coordinate
C     transformations, whitening, ... ) of the data have been performed
C     in advance.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Determines whether the values of the parameters RANK and
C             TOL are to be specified by the user or computed by the
C             routine as follows:
C             = 'R':  Compute RANK only;
C             = 'T':  Compute TOL only;
C             = 'B':  Compute both RANK and TOL;
C             = 'N':  Compute neither RANK nor TOL.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in the data matrix A and the
C             observation matrix B.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns in the data matrix A.  N >= 0.
C
C     L       (input) INTEGER
C             The number of columns in the observation matrix B.
C             L >= 0.
C
C     RANK    (input/output) INTEGER
C             On entry, if JOB = 'T' or JOB = 'N', then RANK must
C             specify r, the rank of the TLS approximation [A+DA|B+DB].
C             RANK <= min(M,N).
C             Otherwise, r is computed by the routine.
C             On exit, if JOB = 'R' or JOB = 'B', and INFO = 0, then
C             RANK contains the computed (effective) rank of the TLS
C             approximation [A+DA|B+DB].
C             Otherwise, the user-supplied value of RANK may be
C             changed by the routine on exit if the RANK-th and the
C             (RANK+1)-th singular values of C = [A|B] are considered
C             to be equal, or if the upper triangular matrix F (as
C             defined in METHOD) is (numerically) singular.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N+L)
C             On entry, the leading M-by-(N+L) part of this array must
C             contain the matrices A and B. Specifically, the first N
C             columns must contain the data matrix A and the last L
C             columns the observation matrix B (right-hand sides).
C             On exit, the leading (N+L)-by-(N+L) part of this array
C             contains the (transformed) right singular vectors,
C             including null space vectors, if any, of C = [A|B].
C             Specifically, the leading (N+L)-by-RANK part of this array
C             always contains the first RANK right singular vectors,
C             corresponding to the largest singular values of C. If
C             L = 0, or if RANK = 0 and IWARN <> 2, the remaining
C             (N+L)-by-(N+L-RANK) top-right part of this array contains
C             the remaining N+L-RANK right singular vectors. Otherwise,
C             this part contains the matrix V2 transformed as described
C             in Step 3 of the TLS algorithm (see METHOD).
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= max(1,M,N+L).
C
C     S       (output) DOUBLE PRECISION array, dimension (min(M,N+L))
C             If INFO = 0, the singular values of matrix C, ordered
C             such that S(1) >= S(2) >= ... >= S(p-1) >= S(p) >= 0,
C             where p = min(M,N+L).
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,L)
C             If INFO = 0, the leading N-by-L part of this array
C             contains the solution X to the TLS problem specified
C             by A and B.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= max(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used to determine the rank of the TLS
C             approximation [A+DA|B+DB] and to check the multiplicity
C             of the singular values of matrix C. Specifically, S(i)
C             and S(j) (i < j) are considered to be equal if
C             SQRT(S(i)**2 - S(j)**2) <= TOL, and the TLS approximation
C             [A+DA|B+DB] has rank r if S(i) > TOL*S(1) (or S(i) > TOL,
C             if TOL specifies sdev (see below)), for i = 1,2,...,r.
C             TOL is also used to check the singularity of the upper
C             triangular matrix F (as defined in METHOD).
C             If JOB = 'R' or JOB = 'N', then TOL must specify the
C             desired tolerance. If the user sets TOL to be less than or
C             equal to 0, the tolerance is taken as EPS, where EPS is
C             the machine precision (see LAPACK Library routine DLAMCH).
C             Otherwise, the tolerance is computed by the routine and
C             the user must supply the non-negative value sdev, i.e. the
C             estimated standard deviation of the error on each element
C             of the matrix C, as input value of TOL.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (L)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, and DWORK(2) returns the reciprocal of the
C             condition number of the matrix F.
C             If INFO > 0, DWORK(1:min(M,N+L)-1) contain the unconverged
C             non-diagonal elements of the bidiagonal matrix whose
C             diagonal is in S (see LAPACK Library routine DGESVD).
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK = max(2, 3*(N+L) + M, 5*(N+L)),       if M >= N+L;
C             LDWORK = max(2, M*(N+L) + max( 3M+N+L, 5*M), 3*L),
C                                                          if M <  N+L.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warnings;
C             = 1:  if the rank of matrix C has been lowered because a
C                   singular value of multiplicity greater than 1 was
C                   found;
C             = 2:  if the rank of matrix C has been lowered because the
C                   upper triangular matrix F is (numerically) singular.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if the SVD algorithm (in LAPACK Library routine
C                   DBDSQR) has failed to converge. In this case, S(1),
C                   S(2), ..., S(INFO) may not have been found
C                   correctly and the remaining singular values may
C                   not be the smallest. This failure is not likely
C                   to occur.
C
C     METHOD
C
C     The method used is an extension (see [3,4,5]) of the classical
C     TLS algorithm proposed by Golub and Van Loan [1].
C
C     Let [A|B] denote the matrix formed by adjoining the columns of B
C     to the columns of A on the right.
C
C     Total Least Squares (TLS) definition:
C     -------------------------------------
C
C       Given matrices A and B, find a matrix X satisfying
C
C            (A + DA) X = B + DB,
C
C       where A and DA are M-by-N matrices, B and DB are M-by-L matrices
C       and X is an N-by-L matrix.
C       The solution X must be such that the Frobenius norm of [DA|DB]
C       is a minimum and each column of B + DB is in the range of
C       A + DA. Whenever the solution is not unique, the routine singles
C       out the minimum norm solution X.
C
C     Define matrix C = [A|B] and s(i) as its i-th singular value for
C     i = 1,2,...,min(M,NL), where NL = N + L. If M < NL, then s(j) = 0
C     for j = M+1,...,NL.
C
C     The Classical TLS algorithm proceeds as follows (see [3,4,5]):
C
C     Step 1: Compute part of the singular value decomposition (SVD)
C             USV' of C = [A|B], namely compute S and V'. (An initial
C             QR factorization of C is used when M is larger enough
C             than NL.)
C
C     Step 2: If not fixed by the user, compute the rank r0 of the data
C             [A|B] based on TOL as follows: if JOB = 'R' or JOB = 'N',
C
C                s(1) >= ... >= s(r0) > TOL*s(1) >= ... >= s(NL).
C
C             Otherwise, using [2], TOL can be computed from the
C             standard deviation sdev of the errors on [A|B]:
C
C                TOL = SQRT(2 * max(M,NL)) * sdev,
C
C             and the rank r0 is determined (if JOB = 'R' or 'B') using
C
C                s(1) >= ... >= s(r0) > TOL >= ... >= s(NL).
C
C             The rank r of the approximation [A+DA|B+DB] is then equal
C             to the minimum of N and r0.
C
C     Step 3: Let V2 be the matrix of the columns of V corresponding to
C             the (NL - r) smallest singular values of C, i.e. the last
C             (NL - r) columns of V.
C             Compute with Householder transformations the orthogonal
C             matrix Q such that:
C
C                       |VH   Y|
C              V2 x Q = |      |
C                       |0    F|
C
C             where VH is an N-by-(N - r) matrix, Y is an N-by-L matrix
C             and F is an L-by-L upper triangular matrix.
C             If F is singular, then lower the rank r with the
C             multiplicity of s(r) and repeat this step.
C
C     Step 4: If F is nonsingular then the solution X is obtained by
C             solving the following equations by forward elimination:
C
C                X F = -Y.
C
C     Notes :
C     The TLS solution is unique if r = N, F is nonsingular and
C     s(N) > s(N+1).
C     If F is singular, however, then the computed solution is infinite
C     and hence does not satisfy the second TLS criterion (see TLS
C     definition). For these cases, Golub and Van Loan [1] claim that
C     the TLS problem has no solution. The properties of these so-called
C     nongeneric problems are described in [4] and the TLS computations
C     are generalized in order to solve them. As proven in [4], the
C     proposed generalization satisfies the TLS criteria for any
C     number L of observation vectors in B provided that, in addition,
C     the solution | X| is constrained to be orthogonal to all vectors
C                  |-I|
C     of the form |w| which belong to the space generated by the columns
C                 |0|
C     of the submatrix |Y|.
C                      |F|
C
C     REFERENCES
C
C     [1] Golub, G.H. and Van Loan, C.F.
C         An Analysis of the Total Least-Squares Problem.
C         SIAM J. Numer. Anal., 17, pp. 883-893, 1980.
C
C     [2] Staar, J., Vandewalle, J. and Wemans, M.
C         Realization of Truncated Impulse Response Sequences with
C         Prescribed Uncertainty.
C         Proc. 8th IFAC World Congress, Kyoto, I, pp. 7-12, 1981.
C
C     [3] Van Huffel, S.
C         Analysis of the Total Least Squares Problem and its Use in
C         Parameter Estimation.
C         Doctoral dissertation, Dept. of Electr. Eng., Katholieke
C         Universiteit Leuven, Belgium, June 1987.
C
C     [4] Van Huffel, S. and Vandewalle, J.
C         Analysis and Solution of the Nongeneric Total Least Squares
C         Problem.
C         SIAM J. Matr. Anal. and Appl., 9, pp. 360-372, 1988.
C
C     [5] Van Huffel, S. and Vandewalle, J.
C         The Total Least Squares Problem: Computational Aspects and
C         Analysis.
C         Series "Frontiers in Applied Mathematics", Vol. 9,
C         SIAM, Philadelphia, 1991.
C
C     NUMERICAL ASPECTS
C
C     The algorithm consists in (backward) stable steps.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB02AD by S. Van Huffel, Katholieke
C     University, Leuven, Belgium.
C
C     REVISIONS
C
C     June 24, 1997, Feb. 27, 2000, Oct. 19, 2003, Feb. 21, 2004.
C
C     KEYWORDS
C
C     Least-squares approximation, singular subspace, singular value
C     decomposition, singular values, total least-squares.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOB
      INTEGER           INFO, IWARN, L, LDC, LDWORK, LDX, M, N, RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  C(LDC,*), DWORK(*), S(*), X(LDX,*)
C     .. Local Scalars ..
      LOGICAL           CRANK, CTOL, LJOBN, LJOBR, LJOBT
      INTEGER           ITAU, J, JWORK, LDW, K, MINMNL, N1, NL, P, R1,
     $                  WRKOPT
      DOUBLE PRECISION  FNORM, RCOND, SMAX, TOLTMP
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE, DLANTR
      EXTERNAL          DLAMCH, DLANGE, DLANTR, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGERQF, DGESVD, DLACPY, DLASET, DORMRQ, DSWAP,
     $                  DTRCON, DTRSM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      IWARN = 0
      INFO = 0
      NL = N + L
      K = MAX( M, NL )
      P = MIN( M, N )
      MINMNL = MIN( M, NL )
      LDW = MAX( 3*MINMNL + K, 5*MINMNL )
      LJOBR = LSAME( JOB, 'R' )
      LJOBT = LSAME( JOB, 'T' )
      LJOBN = LSAME( JOB, 'N' )
C
C     Determine whether RANK or/and TOL is/are to be computed.
C
      CRANK = .NOT.LJOBT .AND. .NOT.LJOBN
      CTOL  = .NOT.LJOBR .AND. .NOT.LJOBN
C
C     Test the input scalar arguments.
C
      IF( CTOL .AND. CRANK .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( L.LT.0 ) THEN
         INFO = -4
      ELSE IF( .NOT.CRANK .AND. RANK.GT.P ) THEN
         INFO = -5
      ELSE IF( LDC.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( CTOL .AND. TOL.LT.ZERO ) THEN
         INFO = -11
      ELSE IF( ( M.GE.NL .AND. LDWORK.LT.MAX( 2, LDW ) ).OR.
     $         ( M.LT.NL .AND. LDWORK.LT.MAX( 2, M*NL + LDW, 3*L ) ) )
     $       THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB02MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( CRANK )
     $   RANK = P
      IF ( MIN( M, NL ).EQ.0 ) THEN
         IF ( M.EQ.0 ) THEN
            CALL DLASET( 'Full', NL, NL, ZERO, ONE, C, LDC )
            CALL DLASET( 'Full', N, L, ZERO, ZERO, X, LDX )
         END IF
         DWORK(1) = TWO
         DWORK(2) = ONE
         RETURN
      END IF
C
C     Subroutine MB02MD solves a set of linear equations by a Total
C     Least Squares Approximation.
C
C     Step 1: Compute part of the singular value decomposition (SVD)
C             USV' of C = [A   |B   ], namely compute S and V'.
C                           M,N  M,L
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( M.GE.NL ) THEN
C
C        M >= N + L:  Overwrite V' on C.
C        Workspace: need max(3*min(M,N+L) + max(M,N+L), 5*min(M,N+L)).
C
         JWORK = 1
         CALL DGESVD( 'No left vectors', 'Overwritten on C', M, NL, C,
     $                LDC, S, DWORK, 1, DWORK, 1, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
      ELSE
C
C        M < N + L:  Save C in the workspace and compute V' in C.
C        Note that the previous DGESVD call cannot be used in this case.
C        Workspace: need M*(N+L) + max(3*min(M,N+L) + max(M,N+L),
C                                      5*min(M,N+L)).
C
         CALL DLACPY( 'Full', M, NL, C, LDC, DWORK, M )
         JWORK = M*NL + 1
         CALL DGESVD( 'No left vectors', 'All right vectors', M, NL,
     $                DWORK, M, S, DWORK, 1, C, LDC, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
      END IF
C
      IF ( INFO.GT.0 ) THEN
C
C        Save the unconverged non-diagonal elements of the bidiagonal
C        matrix and exit.
C
         DO 10 J = 1, MINMNL - 1
            DWORK(J) = DWORK(JWORK+J)
   10    CONTINUE
C
         RETURN
      END IF
      WRKOPT = MAX( 2, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C     Transpose V' in-situ (in C).
C
      DO 20 J = 2, NL
         CALL DSWAP( J-1, C(J,1), LDC, C(1,J), 1 )
   20 CONTINUE
C
C     Step 2: Compute the rank of the approximation [A+DA|B+DB].
C
      IF ( CTOL ) THEN
         TOLTMP = SQRT( TWO*DBLE( K ) )*TOL
         SMAX = TOLTMP
      ELSE
         TOLTMP = TOL
         IF ( TOLTMP.LE.ZERO ) TOLTMP = DLAMCH( 'Precision' )
         SMAX = MAX( TOLTMP*S(1), DLAMCH( 'Safe minimum' ) )
      END IF
C
      IF ( CRANK ) THEN
C        WHILE ( RANK .GT. 0 ) .AND. ( S(RANK) .LE. SMAX ) DO
   40    IF ( RANK.GT.0 ) THEN
            IF ( S(RANK).LE.SMAX ) THEN
               RANK = RANK - 1
               GO TO 40
            END IF
         END IF
C        END WHILE 40
      END IF
C
      IF ( L.EQ.0 ) THEN
         DWORK(1) = WRKOPT
         DWORK(2) = ONE
         RETURN
      END IF
C
      N1 = N + 1
      ITAU  = 1
      JWORK = ITAU + L
C
C     Step 3: Compute the orthogonal matrix Q and matrices F and Y
C     such that F is nonsingular.
C
C     REPEAT
C
C        Adjust the rank if S(RANK) has multiplicity greater than 1.
C
   60    CONTINUE
         R1 = RANK + 1
         IF ( RANK.LT.MINMNL ) THEN
C           WHILE RANK.GT.0 .AND. S(RANK)**2 - S(R1)**2.LE.TOL**2 DO
   80       IF ( RANK.GT.0 ) THEN
               IF ( ONE - ( S(R1)/S(RANK) )**2.LE.( TOLTMP/S(RANK) )**2
     $            ) THEN
                  RANK = RANK - 1
                  IWARN = 1
                  GO TO 80
               END IF
            END IF
C           END WHILE 80
         END IF
C
         IF ( RANK.EQ.0 ) THEN
C
C           Return zero solution.
C
            CALL DLASET( 'Full', N, L, ZERO, ZERO, X, LDX )
            DWORK(1) = WRKOPT
            DWORK(2) = ONE
            RETURN
         END IF
C
C        Compute the orthogonal matrix Q (in factorized form) and the
C        matrices F and Y using RQ factorization. It is assumed that,
C        generically, the last L rows of V2 matrix have full rank.
C        The code could not be the most efficient one when RANK has been
C        lowered, because the already created zero pattern of the last
C        L rows of V2 matrix is not exploited.
C        Workspace: need 2*L;  prefer L + L*NB.
C
         R1 = RANK + 1
         CALL DGERQF( L, NL-RANK, C(N1,R1), LDC, DWORK(ITAU),
     $                DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Workspace: need N+L;  prefer L + N*NB.
C
         CALL DORMRQ( 'Right', 'Transpose', N, NL-RANK, L, C(N1,R1),
     $                LDC, DWORK(ITAU), C(1,R1), LDC, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
         CALL DLASET( 'Full', L, N-RANK, ZERO, ZERO, C(N1,R1), LDC )
         IF ( L.GT.1 )
     $      CALL DLASET( 'Lower', L-1, L-1, ZERO, ZERO, C(N1+1,N1),
     $                   LDC )
C
C        Estimate the reciprocal condition number of the matrix F,
C        and lower the rank if F can be considered as singular.
C        Workspace: need 3*L.
C
         CALL DTRCON( '1-norm', 'Upper', 'Non-unit', L, C(N1,N1), LDC,
     $                RCOND, DWORK, IWORK, INFO )
         WRKOPT = MAX( WRKOPT, 3*L )
C
         FNORM = DLANTR( '1-norm', 'Upper', 'Non-unit', L, L, C(N1,N1),
     $                   LDC, DWORK )
         IF ( RCOND.LE.TOLTMP*FNORM ) THEN
            RANK = RANK - 1
            IWARN = 2
            GO TO 60
         ELSE IF ( FNORM.LE.TOLTMP*DLANGE( '1-norm', N, L, C(1,N1), LDC,
     $                                     DWORK ) ) THEN
            RANK = RANK - L
            IWARN = 2
            GO TO 60
         END IF
C     UNTIL ( F nonsingular, i.e., RCOND.GT.TOL*FNORM or
C                                  FNORM.GT.TOL*norm(Y) )
C
C     Step 4: Solve X F = -Y by forward elimination,
C             (F is upper triangular).
C
      CALL DLACPY( 'Full', N, L, C(1,N1), LDC, X, LDX )
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', N, L,
     $            -ONE, C(N1,N1), LDC, X, LDX )
C
C     Set the optimal workspace and reciprocal condition number of F.
C
      DWORK(1) = WRKOPT
      DWORK(2) = RCOND
C
      RETURN
C *** Last line of MB02MD ***
      END
