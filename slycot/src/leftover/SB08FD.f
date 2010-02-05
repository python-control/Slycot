      SUBROUTINE SB08FD( DICO, N, M, P, ALPHA, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, NQ, NR, CR, LDCR, DR, LDDR, TOL, DWORK,
     $                   LDWORK, IWARN, INFO )
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
C     To construct, for a given system G = (A,B,C,D), a feedback
C     matrix F and an orthogonal transformation matrix Z, such that
C     the systems
C
C          Q = (Z'*(A+B*F)*Z, Z'*B, (C+D*F)*Z, D)
C     and
C          R = (Z'*(A+B*F)*Z, Z'*B, F*Z, I)
C
C     provide a stable right coprime factorization of G in the form
C                       -1
C              G = Q * R  ,
C
C     where G, Q and R are the corresponding transfer-function matrices.
C     The resulting state dynamics matrix of the systems Q and R has
C     eigenvalues lying inside a given stability domain.
C     The Z matrix is not explicitly computed.
C
C     Note: If the given state-space representation is not stabilizable,
C     the unstabilizable part of the original system is automatically
C     deflated and the order of the systems Q and R is accordingly
C     reduced.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the state vector, i.e. the order of the
C             matrix A, and also the number of rows of the matrix B and
C             the number of columns of the matrices C and CR.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector, i.e. the number of columns
C             of the matrices B, D and DR and the number of rows of the
C             matrices CR and DR.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of output vector, i.e. the number of rows
C             of the matrices C and D.  P >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION array, dimension (2)
C             ALPHA(1) contains the desired stability degree to be
C             assigned for the eigenvalues of A+B*F, and ALPHA(2)
C             the stability margin. The eigenvalues outside the
C             ALPHA(2)-stability region will be assigned to have the
C             real parts equal to ALPHA(1) < 0 and unmodified
C             imaginary parts for a continuous-time system
C             (DICO = 'C'), or moduli equal to 0 <= ALPHA(2) < 1
C             for a discrete-time system (DICO = 'D').
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading NQ-by-NQ part of this array contains
C             the leading NQ-by-NQ part of the matrix Z'*(A+B*F)*Z, the
C             state dynamics matrix of the numerator factor Q, in a
C             real Schur form. The trailing NR-by-NR part of this matrix
C             represents the state dynamics matrix of a minimal
C             realization of the denominator factor R.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix.
C             On exit, the leading NQ-by-M part of this array contains
C             the leading NQ-by-M part of the matrix Z'*B, the
C             input/state matrix of the numerator factor Q. The last
C             NR rows of this matrix form the input/state matrix of
C             a minimal realization of the denominator factor R.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, the leading P-by-NQ part of this array contains
C             the leading P-by-NQ part of the matrix (C+D*F)*Z,
C             the state/output matrix of the numerator factor Q.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             input/output matrix. D represents also the input/output
C             matrix of the numerator factor Q.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     NQ      (output) INTEGER
C             The order of the resulting factors Q and R.
C             Generally, NQ = N - NS, where NS is the number of
C             uncontrollable eigenvalues outside the stability region.
C
C     NR      (output) INTEGER
C             The order of the minimal realization of the factor R.
C             Generally, NR is the number of controllable eigenvalues
C             of A outside the stability region (the number of modified
C             eigenvalues).
C
C     CR      (output) DOUBLE PRECISION array, dimension (LDCR,N)
C             The leading M-by-NQ part of this array contains the
C             leading M-by-NQ part of the feedback matrix F*Z, which
C             moves the eigenvalues of A lying outside the ALPHA-stable
C             region to values which are on the ALPHA-stability
C             boundary.  The last NR columns of this matrix form the
C             state/output matrix of a minimal realization of the
C             denominator factor R.
C
C     LDCR    INTEGER
C             The leading dimension of array CR.  LDCR >= MAX(1,M).
C
C     DR      (output) DOUBLE PRECISION array, dimension (LDDR,M)
C             The leading M-by-M part of this array contains an
C             identity matrix representing the input/output matrix
C             of the denominator factor R.
C
C     LDDR    INTEGER
C             The leading dimension of array DR.  LDDR >= MAX(1,M).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of
C             B are considered zero (used for controllability tests).
C             If the user sets TOL <= 0, then an implicitly computed,
C             default tolerance, defined by  TOLDEF = N*EPS*NORM(B),
C             is used instead, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH) and NORM(B) denotes
C             the 1-norm of B.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.
C             LWORK >= MAX( 1, N*(N+5), 5*M, 4*P ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = K:  K violations of the numerical stability condition
C                   NORM(F) <= 10*NORM(A)/NORM(B) occured during the
C                   assignment of eigenvalues.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A to a real Schur form failed;
C             = 2:  a failure was detected during the ordering of the
C                   real Schur form of A, or in the iterative process
C                   for reordering the eigenvalues of Z'*(A + B*F)*Z
C                   along the diagonal.
C
C     METHOD
C
C     The subroutine is based on the factorization algorithm of [1].
C
C     REFERENCES
C
C     [1] Varga A.
C         Coprime factors model reduction method based on
C         square-root balancing-free techniques.
C         System Analysis, Modelling and Simulation,
C         vol. 11, pp. 303-311, 1993.
C
C     NUMERICAL ASPECTS
C                                            3
C     The algorithm requires no more than 14N  floating point
C     operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine RCFS.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     Mar. 2003, May 2003, A. Varga, German Aerospace Center.
C     May 2003, V. Sima, Research Institute for Informatics, Bucharest.
C     Sep. 2005, A. Varga, German Aerospace Center.
C
C     KEYWORDS
C
C     Coprime factorization, eigenvalue, eigenvalue assignment,
C     feedback control, pole placement, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TEN, ZERO
      PARAMETER         ( ONE = 1.0D0, TEN = 1.0D1, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDCR, LDD, LDDR,
     $                  LDWORK, M, N, NQ, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ALPHA(*), B(LDB,*), C(LDC,*),
     $                  CR(LDCR,*), D(LDD,*), DR(LDDR,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           DISCR
      INTEGER           I, IB, IB1, J, K, KFI, KG, KW, KWI, KWR, KZ, L,
     $                  L1, NB, NCUR, NCUR1, NFP, NLOW, NMOVES, NSUP
      DOUBLE PRECISION  BNORM, CS, PR, RMAX, SM, SN, TOLER, WRKOPT, X, Y
C     .. Local Arrays ..
      DOUBLE PRECISION  A2(2,2), Z(4,4)
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLANGE, DLAPY2
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, DLANGE, DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLACPY, DLAEXC, DLANV2, DLASET, DROT,
     $                  SB01BY, TB01LD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SIGN, SQRT
C
C     .. Executable Statements ..
C
      DISCR = LSAME( DICO, 'D' )
      IWARN = 0
      INFO  = 0
C
C     Check the scalar input parameters.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( ( DISCR .AND. ( ALPHA(1).LT.ZERO .OR. ALPHA(1).GE.ONE
     $                  .OR.   ALPHA(2).LT.ZERO .OR. ALPHA(2).GE.ONE ) )
     $    .OR.
     $    ( .NOT.DISCR .AND. ( ALPHA(1).GE.ZERO .OR. ALPHA(2).GE.ZERO )
     $    ) ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDCR.LT.MAX( 1, M ) ) THEN
         INFO = -17
      ELSE IF( LDDR.LT.MAX( 1, M ) ) THEN
         INFO = -19
      ELSE IF( LDWORK.LT.MAX( 1, N*(N+5), 5*M, 4*P ) ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 )THEN
C
C        Error return.
C
         CALL XERBLA( 'SB08FD', -INFO )
         RETURN
      END IF
C
C     Set DR = I and quick return if possible.
C
      NR = 0
      CALL DLASET( 'Full', M, M, ZERO, ONE, DR, LDDR )
      IF( MIN( N, M ).EQ.0 ) THEN
         NQ = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Set F = 0 in the array CR.
C
      CALL DLASET( 'Full', M, N, ZERO, ZERO, CR, LDCR )
C
C     Compute the norm of B and set the default tolerance if necessary.
C
      BNORM = DLANGE( '1-norm', N, M, B, LDB, DWORK )
      TOLER = TOL
      IF( TOLER.LE.ZERO )
     $   TOLER = DBLE( N ) * BNORM * DLAMCH( 'Epsilon' )
      IF( BNORM.LE.TOLER ) THEN
         NQ = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Compute the bound for the numerical stability condition.
C
      RMAX = TEN * DLANGE( '1-norm', N, N, A, LDA, DWORK ) / BNORM
C
C     Allocate working storage.
C
      KZ  = 1
      KWR = KZ  + N*N
      KWI = KWR + N
      KW  = KWI + N
C
C     Reduce A to an ordered real Schur form using an orthogonal
C     similarity transformation A <- Z'*A*Z and accumulate the
C     transformations in Z.  The separation of spectrum of A is
C     performed such that the leading NFP-by-NFP submatrix of A
C     corresponds to the "stable" eigenvalues which will be not
C     modified. The bottom (N-NFP)-by-(N-NFP) diagonal block of A
C     corresponds to the "unstable" eigenvalues to be modified.
C     Apply the transformation to B and C: B <- Z'*B and C <- C*Z.
C
C     Workspace needed:      N*(N+2);
C     Additional workspace:  need   3*N;
C                            prefer larger.
C
      CALL TB01LD( DICO, 'Stable', 'General', N, M, P, ALPHA(2), A, LDA,
     $             B, LDB, C, LDC, NFP, DWORK(KZ), N, DWORK(KWR),
     $             DWORK(KWI), DWORK(KW), LDWORK-KW+1, INFO )
      IF( INFO.NE.0 )
     $    RETURN
C
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
C     Perform the pole assignment if there exist "unstable" eigenvalues.
C
      NQ = N
      IF( NFP.LT.N ) THEN
         KG  = 1
         KFI = KG  + 2*M
         KW  = KFI + 2*M
C
C        Set the limits for the bottom diagonal block.
C
         NLOW = NFP + 1
         NSUP = N
C
C        WHILE (NLOW <= NSUP) DO
   10    IF( NLOW.LE.NSUP ) THEN
C
C           Main loop for assigning one or two poles.
C
C           Determine the dimension of the last block.
C
            IB = 1
            IF( NLOW.LT.NSUP ) THEN
               IF( A(NSUP,NSUP-1).NE.ZERO ) IB = 2
            END IF
            L = NSUP - IB + 1
C
C           Save the last IB rows of B in G.
C
            CALL DLACPY( 'Full', IB, M, B(L,1), LDB, DWORK(KG), IB )
C
C           Check the controllability of the last block.
C
            IF( DLANGE( '1-norm', IB, M, DWORK(KG), IB, DWORK(KW) )
     $            .LE.TOLER )THEN
C
C              Deflate the uncontrollable block and resume the
C              main loop.
C
               NSUP = NSUP - IB
            ELSE
C
C              Form the IBxIB matrix A2 from the last diagonal block and
C              set the pole(s) to be assigned.
C
               A2(1,1) = A(L,L)
               IF( IB.EQ.1 ) THEN
                  SM = ALPHA(1)
                  IF( DISCR ) SM = SIGN( ALPHA(1), A2(1,1) )
                  PR = ALPHA(1)
               ELSE
                  A2(1,2) = A(L,NSUP)
                  A2(2,1) = A(NSUP,L)
                  A2(2,2) = A(NSUP,NSUP)
                  SM = ALPHA(1) + ALPHA(1)
                  PR = ALPHA(1)*ALPHA(1)
                  IF( DISCR ) THEN
                     X  = A2(1,1)
                     Y  = SQRT( ABS( A2(1,2)*A2(2,1) ) )
                     SM = SM * X / DLAPY2( X, Y )
                  ELSE
                     PR = PR - A2(1,2)*A2(2,1)
                  END IF
               END IF
C
C              Determine the M-by-IB feedback matrix FI which assigns
C              the selected IB poles for the pair (A2,G).
C
C              Workspace needed: 5*M.
C
               CALL SB01BY( IB, M, SM, PR, A2, DWORK(KG), DWORK(KFI),
     $                      TOLER, DWORK(KW), INFO )
               IF( INFO.NE.0 ) THEN
C
C                 Uncontrollable 2x2 block with double real eigenvalues
C                 which due to roundoff appear as a pair of complex
C                 conjugated eigenvalues.
C                 One of them can be elliminated using the information
C                 in DWORK(KFI) and DWORK(KFI+M).
C
                  CS = DWORK(KFI)
                  SN = -DWORK(KFI+M)
C
C                 Apply the Givens transformation to A, B, C and F.
C
                  L1 = L + 1
                  CALL DROT( NSUP-L+1, A(L1,L), LDA, A(L,L),
     $                       LDA, CS, SN )
                  CALL DROT( L1, A(1,L1), 1, A(1,L), 1, CS, SN )
                  CALL DROT( M, B(L1,1), LDB, B(L,1), LDB, CS, SN )
                  IF( P.GT.0 )
     $               CALL DROT( P, C(1,L1), 1, C(1,L), 1, CS, SN )
                  CALL DROT( M, CR(1,L1), 1, CR(1,L), 1, CS, SN )
C
C                 Deflate the uncontrollable block and resume the
C                 main loop.
C
                  A(L1,L) = ZERO
                  NSUP = NSUP - 1
                  INFO = 0
                  GO TO 10
               END IF
C
C              Check for possible numerical instability.
C
               IF( DLANGE( '1-norm', M, IB, DWORK(KFI), M, DWORK(KW) )
     $             .GT.RMAX ) IWARN = IWARN + 1
C
C              Update the feedback matrix F <-- F + [0 FI] in CR.
C
               K = KFI
               DO 30 J = L, L + IB - 1
                  DO 20 I = 1, M
                     CR(I,J) = CR(I,J) + DWORK(K)
                     K = K + 1
   20             CONTINUE
   30          CONTINUE
C
C              Update the state matrix A <-- A + B*[0 FI].
C
               CALL DGEMM( 'NoTranspose', 'NoTranspose', NSUP, IB, M,
     $                     ONE, B, LDB, DWORK(KFI), M, ONE, A(1,L),
     $                     LDA )
               IF( IB.EQ.2 ) THEN
C
C                 Try to split the 2x2 block and standardize it.
C
                  L1 = L + 1
                  CALL DLANV2( A(L,L), A(L,L1), A(L1,L), A(L1,L1),
     $                         X, Y, PR, SM, CS, SN )
C
C                 Apply the transformation to A, B, C and F.
C
                  IF( L1.LT.NSUP )
     $               CALL DROT( NSUP-L1, A(L,L1+1), LDA, A(L1,L1+1),
     $                          LDA, CS, SN )
                  CALL DROT( L-1, A(1,L), 1, A(1,L1), 1, CS, SN )
                  CALL DROT( M, B(L,1), LDB, B(L1,1), LDB, CS, SN )
                  IF( P.GT.0 )
     $               CALL DROT( P, C(1,L), 1, C(1,L1), 1, CS, SN )
                  CALL DROT( M, CR(1,L), 1, CR(1,L1), 1, CS, SN )
               END IF
               IF( NLOW+IB.LE.NSUP ) THEN
C
C                 Move the last block(s) to the leading position(s) of
C                 the bottom block.
C
C                 Workspace:     need MAX(4*N, 4*M, 4*P).
C
                  NCUR1 = NSUP - IB
                  NMOVES = 1
                  IF( IB.EQ.2 .AND. A(NSUP,NSUP-1).EQ.ZERO ) THEN
                     IB = 1
                     NMOVES = 2
                  END IF
C
C                 WHILE (NMOVES > 0) DO
   40             IF( NMOVES.GT.0 ) THEN
                     NCUR = NCUR1
C
C                    WHILE (NCUR >= NLOW) DO
   50                IF( NCUR.GE.NLOW ) THEN
C
C                       Loop for positioning of the last block.
C
C                       Determine the dimension of the current block.
C
                        IB1 = 1
                        IF( NCUR.GT.NLOW ) THEN
                           IF( A(NCUR,NCUR-1).NE.ZERO ) IB1 = 2
                        END IF
                        NB = IB1 + IB
C
C                       Initialize the local transformation matrix Z.
C
                        CALL DLASET( 'Full', NB, NB, ZERO, ONE, Z, 4 )
                        L = NCUR - IB1 + 1
C
C                       Exchange two adjacent blocks and accumulate the
C                       transformations in Z.
C
                        CALL DLAEXC( .TRUE., NB, A(L,L), LDA, Z, 4, 1,
     $                               IB1, IB, DWORK, INFO )
                        IF( INFO.NE.0 ) THEN
                           INFO = 2
                           RETURN
                        END IF
C
C                       Apply the transformation to the rest of A.
C
                        L1 = L + NB
                        IF( L1.LE.NSUP ) THEN
                           CALL DGEMM( 'Transpose', 'NoTranspose', NB,
     $                                 NSUP-L1+1, NB, ONE, Z, 4,
     $                                 A(L,L1), LDA, ZERO, DWORK, NB )
                           CALL DLACPY( 'Full', NB, NSUP-L1+1, DWORK,
     $                                  NB, A(L,L1), LDA )
                        END IF
                        CALL DGEMM( 'NoTranspose', 'NoTranspose', L-1,
     $                              NB, NB, ONE, A(1,L), LDA, Z, 4,
     $                              ZERO, DWORK, N )
                        CALL DLACPY( 'Full', L-1, NB, DWORK, N, A(1,L),
     $                               LDA )
C
C                       Apply the transformation to B, C and F.
C
                        CALL DGEMM( 'Transpose', 'NoTranspose', NB, M,
     $                              NB, ONE, Z, 4, B(L,1), LDB, ZERO,
     $                              DWORK, NB )
                        CALL DLACPY( 'Full', NB, M, DWORK, NB, B(L,1),
     $                               LDB )
C
                        IF( P.GT.0 ) THEN
                           CALL DGEMM( 'NoTranspose', 'NoTranspose', P,
     $                                 NB, NB, ONE, C(1,L), LDC, Z, 4,
     $                                 ZERO, DWORK, P )
                           CALL DLACPY( 'Full', P, NB, DWORK, P,
     $                                  C(1,L), LDC )
                        END IF
C
                        CALL DGEMM( 'NoTranspose', 'NoTranspose', M, NB,
     $                              NB, ONE, CR(1,L), LDCR, Z, 4, ZERO,
     $                              DWORK, M )
                        CALL DLACPY( 'Full', M, NB, DWORK, M, CR(1,L),
     $                               LDCR )
C
                        NCUR = NCUR - IB1
                        GO TO 50
                     END IF
C                    END WHILE 50
C
                     NMOVES = NMOVES - 1
                     NCUR1  = NCUR1 + 1
                     NLOW   = NLOW + IB
                     GO TO 40
                  END IF
C                 END WHILE 40
C
               ELSE
                  NLOW = NLOW + IB
               END IF
            END IF
            GO TO 10
         END IF
C        END WHILE 10
C
         NQ = NSUP
         NR = NSUP - NFP
C
C        Annihilate the elements below the first subdiagonal of A.
C
         IF( NQ.GT.2 )
     $      CALL DLASET( 'Lower', NQ-2, NQ-2, ZERO, ZERO, A(3,1), LDA )
      END IF
C
C     Compute C <-- CQ = C + D*F.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NQ, M, ONE, D, LDD,
     $            CR, LDCR, ONE, C, LDC )
C
      DWORK(1) = MAX( WRKOPT, DBLE( MAX( 5*M, 4*P ) ) )
C
      RETURN
C *** Last line of SB08FD ***
      END
