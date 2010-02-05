      SUBROUTINE SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI,
     $                   NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK,
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
C     To determine the state feedback matrix F for a given system (A,B)
C     such that the closed-loop state matrix A+B*F has specified
C     eigenvalues.
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
C             the number of columns of the matrix F.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector, i.e. the number of columns
C             of the matrix B and the number of rows of the matrix F.
C             M >= 0.
C
C     NP      (input) INTEGER
C             The number of given eigenvalues. At most N eigenvalues
C             can be assigned.  0 <= NP.
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the maximum admissible value, either for real
C             parts, if DICO = 'C', or for moduli, if DICO = 'D',
C             of the eigenvalues of A which will not be modified by
C             the eigenvalue assignment algorithm.
C             ALPHA >= 0 if DICO = 'D'.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Z'*(A+B*F)*Z in a real Schur form.
C             The leading NFP-by-NFP diagonal block of A corresponds
C             to the fixed (unmodified) eigenvalues having real parts
C             less than ALPHA, if DICO = 'C', or moduli less than ALPHA,
C             if DICO = 'D'. The trailing NUP-by-NUP diagonal block of A
C             corresponds to the uncontrollable eigenvalues detected by
C             the eigenvalue assignment algorithm. The elements under
C             the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     WR,WI   (input/output) DOUBLE PRECISION array, dimension (NP)
C             On entry, these arrays must contain the real and imaginary
C             parts, respectively, of the desired eigenvalues of the
C             closed-loop system state-matrix A+B*F. The eigenvalues
C             can be unordered, except that complex conjugate pairs
C             must appear consecutively in these arrays.
C             On exit, if INFO = 0, the leading NAP elements of these
C             arrays contain the real and imaginary parts, respectively,
C             of the assigned eigenvalues. The trailing NP-NAP elements
C             contain the unassigned eigenvalues.
C
C     NFP     (output) INTEGER
C             The number of eigenvalues of A having real parts less than
C             ALPHA, if DICO = 'C', or moduli less than ALPHA, if
C             DICO = 'D'. These eigenvalues are not modified by the
C             eigenvalue assignment algorithm.
C
C     NAP     (output) INTEGER
C             The number of assigned eigenvalues. If INFO = 0 on exit,
C             then NAP = N-NFP-NUP.
C
C     NUP     (output) INTEGER
C             The number of uncontrollable eigenvalues detected by the
C             eigenvalue assignment algorithm (see METHOD).
C
C     F       (output) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading M-by-N part of this array contains the state
C             feedback F, which assigns NAP closed-loop eigenvalues and
C             keeps unaltered N-NAP open-loop eigenvalues.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N)
C             The leading N-by-N part of this array contains the
C             orthogonal matrix Z which reduces the closed-loop
C             system state matrix A + B*F to upper real Schur form.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of A
C             or B are considered zero (used for controllability tests).
C             If the user sets TOL <= 0, then the default tolerance
C             TOL = N * EPS * max(NORM(A),NORM(B)) is used, where EPS is
C             the machine precision (see LAPACK Library routine DLAMCH)
C             and NORM(A) denotes the 1-norm of A.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of working array DWORK.
C             LDWORK >= MAX( 1,5*M,5*N,2*N+4*M ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = K:  K violations of the numerical stability condition
C                   NORM(F) <= 100*NORM(A)/NORM(B) occured during the
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
C             = 3:  the number of eigenvalues to be assigned is less
C                   than the number of possibly assignable eigenvalues;
C                   NAP eigenvalues have been properly assigned,
C                   but some assignable eigenvalues remain unmodified.
C             = 4:  an attempt is made to place a complex conjugate
C                   pair on the location of a real eigenvalue. This
C                   situation can only appear when N-NFP is odd,
C                   NP > N-NFP-NUP is even, and for the last real
C                   eigenvalue to be modified there exists no available
C                   real eigenvalue to be assigned. However, NAP
C                   eigenvalues have been already properly assigned.
C
C     METHOD
C
C     SB01BD is based on the factorization algorithm of [1].
C     Given the matrices A and B of dimensions N-by-N and N-by-M,
C     respectively, this subroutine constructs an M-by-N matrix F such
C     that A + BF has eigenvalues as follows.
C     Let NFP eigenvalues of A have real parts less than ALPHA, if
C     DICO = 'C', or moduli less then ALPHA, if DICO = 'D'. Then:
C     1) If the pair (A,B) is controllable, then A + B*F has
C        NAP = MIN(NP,N-NFP) eigenvalues assigned from those specified
C        by WR + j*WI and N-NAP unmodified eigenvalues;
C     2) If the pair (A,B) is uncontrollable, then the number of
C        assigned eigenvalues NAP satifies generally the condition
C        NAP <= MIN(NP,N-NFP).
C
C     At the beginning of the algorithm, F = 0 and the matrix A is
C     reduced to an ordered real Schur form by separating its spectrum
C     in two parts. The leading NFP-by-NFP part of the Schur form of
C     A corresponds to the eigenvalues which will not be modified.
C     These eigenvalues have real parts less than ALPHA, if
C     DICO = 'C', or moduli less than ALPHA, if DICO = 'D'.
C     The performed orthogonal transformations are accumulated in Z.
C     After this preliminary reduction, the algorithm proceeds
C     recursively.
C
C     Let F be the feedback matrix at the beginning of a typical step i.
C     At each step of the algorithm one real eigenvalue or two complex
C     conjugate eigenvalues are placed by a feedback Fi of rank 1 or
C     rank 2, respectively. Since the feedback Fi affects only the
C     last 1 or 2 columns of Z'*(A+B*F)*Z, the matrix Z'*(A+B*F+B*Fi)*Z
C     therefore remains in real Schur form. The assigned eigenvalue(s)
C     is (are) then moved to another diagonal position of the real
C     Schur form using reordering techniques and a new block is
C     transfered in the last diagonal position. The feedback matrix F
C     is updated as F <-- F + Fi. The eigenvalue(s) to be assigned at
C     each step is (are) chosen such that the norm of each Fi is
C     minimized.
C
C     If uncontrollable eigenvalues are encountered in the last diagonal
C     position of the real Schur matrix Z'*(A+B*F)*Z, the algorithm
C     deflates them at the bottom of the real Schur form and redefines
C     accordingly the position of the "last" block.
C
C     Note: Not all uncontrollable eigenvalues of the pair (A,B) are
C     necessarily detected by the eigenvalue assignment algorithm.
C     Undetected uncontrollable eigenvalues may exist if NFP > 0 and/or
C     NP < N-NFP.
C
C     REFERENCES
C
C     [1] Varga A.
C         A Schur method for pole assignment.
C         IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519, 1981.
C
C     NUMERICAL ASPECTS
C                                            3
C     The algorithm requires no more than 14N  floating point
C     operations. Although no proof of numerical stability is known,
C     the algorithm has always been observed to yield reliable
C     numerical results.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     February 1999. Based on the RASP routine SB01BD.
C
C     REVISIONS
C
C     March 30, 1999, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     April 4, 1999. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C     May 18, 2003. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C     Feb. 15, 2004, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     May 12, 2005. A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen.
C
C     KEYWORDS
C
C     Eigenvalues, eigenvalue assignment, feedback control,
C     pole placement, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION HUNDR, ONE, TWO, ZERO
      PARAMETER        ( HUNDR = 1.0D2, ONE = 1.0D0, TWO = 2.0D0,
     $                   ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO
      INTEGER          INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N,
     $                 NAP, NFP, NP, NUP
      DOUBLE PRECISION ALPHA, TOL
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*),
     $                 WI(*), WR(*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL          CEIG, DISCR, SIMPLB
      INTEGER          I, IB, IB1, IERR, IPC, J, K, KFI, KG, KW, KWI,
     $                 KWR, NCUR, NCUR1, NL, NLOW, NMOVES, NPC, NPR,
     $                 NSUP, WRKOPT
      DOUBLE PRECISION ANORM, BNORM, C, P, RMAX, S, X, Y, TOLER, TOLERB
C     .. Local Arrays ..
      LOGICAL          BWORK(1)
      DOUBLE PRECISION A2(2,2)
C     .. External Functions ..
      LOGICAL          LSAME, SELECT
      DOUBLE PRECISION DLAMCH, DLANGE
      EXTERNAL         DLAMCH, DLANGE, LSAME, SELECT
C     .. External Subroutines ..
      EXTERNAL         DGEES, DGEMM, DLAEXC, DLASET, DROT, DSWAP,
     $                 MB03QD, MB03QY, SB01BX, SB01BY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        DBLE, INT, MAX
C     ..
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
      ELSE IF( NP.LT.0 ) THEN
         INFO = -4
      ELSE IF( DISCR .AND. ( ALPHA.LT.ZERO ) ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -16
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -18
      ELSE IF( LDWORK.LT.MAX( 1, 5*M, 5*N, 2*N + 4*M ) ) THEN
         INFO = -21
      END IF
      IF( INFO.NE.0 )THEN
C
C        Error return.
C
         CALL XERBLA( 'SB01BD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         NFP = 0
         NAP = 0
         NUP = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Compute the norms of A and B, and set default tolerances
C     if necessary.
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK )
      BNORM = DLANGE( '1-norm', N, M, B, LDB, DWORK )
      IF( TOL.LE.ZERO ) THEN
         X = DLAMCH( 'Epsilon' )
         TOLER  = DBLE( N ) * MAX( ANORM, BNORM ) * X
         TOLERB = DBLE( N ) * BNORM * X
      ELSE
         TOLER  = TOL
         TOLERB = TOL
      END IF
C
C     Allocate working storage.
C
      KWR = 1
      KWI = KWR + N
      KW  = KWI + N
C
C     Reduce A to real Schur form using an orthogonal similarity
C     transformation A <- Z'*A*Z and accumulate the transformation in Z.
C
C     Workspace:  need   5*N;
C                 prefer larger.
C
      CALL DGEES( 'Vectors', 'No ordering', SELECT, N, A, LDA, NCUR,
     $            DWORK(KWR), DWORK(KWI), Z, LDZ, DWORK(KW),
     $            LDWORK-KW+1, BWORK, INFO )
      WRKOPT = KW - 1 + INT( DWORK( KW ) )
      IF( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
C
C     Reduce A to an ordered real Schur form using an orthogonal
C     similarity transformation A <- Z'*A*Z and accumulate the
C     transformations in Z. The separation of the spectrum of A is
C     performed such that the leading NFP-by-NFP submatrix of A
C     corresponds to the "good" eigenvalues which will not be
C     modified. The bottom (N-NFP)-by-(N-NFP) diagonal block of A
C     corresponds to the "bad" eigenvalues to be modified.
C
C     Workspace needed:  N.
C
      CALL MB03QD( DICO, 'Stable', 'Update', N, 1, N, ALPHA,
     $             A, LDA, Z, LDZ, NFP, DWORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
C
C     Set F = 0.
C
      CALL DLASET( 'Full', M, N, ZERO, ZERO, F, LDF )
C
C     Return if B is negligible (uncontrollable system).
C
      IF( BNORM.LE.TOLERB ) THEN
         NAP = 0
         NUP = N
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Compute the bound for the numerical stability condition.
C
      RMAX = HUNDR * ANORM / BNORM
C
C     Perform eigenvalue assignment if there exist "bad" eigenvalues.
C
      NAP = 0
      NUP = 0
      IF( NFP .LT. N ) THEN
         KG  = 1
         KFI = KG  + 2*M
         KW  = KFI + 2*M
C
C        Set the limits for the bottom diagonal block.
C
         NLOW = NFP + 1
         NSUP = N
C
C        Separate and count real and complex eigenvalues to be assigned.
C
         NPR = 0
         DO 10 I = 1, NP
            IF( WI(I) .EQ. ZERO ) THEN
               NPR = NPR + 1
               K = I - NPR
               IF( K .GT. 0 ) THEN
                  S = WR(I)
                  DO 5 J = NPR + K - 1, NPR, -1
                     WR(J+1) = WR(J)
                     WI(J+1) = WI(J)
    5             CONTINUE
                  WR(NPR) = S
                  WI(NPR) = ZERO
               END IF
            END IF
   10    CONTINUE
         NPC = NP - NPR
C
C        The first NPR elements of WR and WI contain the real
C        eigenvalues, the last NPC elements contain the complex
C        eigenvalues. Set the pointer to complex eigenvalues.
C
         IPC = NPR + 1
C
C        Main loop for assigning one or two eigenvalues.
C
C        Terminate if all eigenvalues were assigned, or if there
C        are no more eigenvalues to be assigned, or if a non-fatal
C        error condition was set.
C
C        WHILE (NLOW <= NSUP and INFO = 0) DO
C
   20    IF( NLOW.LE.NSUP .AND. INFO.EQ.0 ) THEN
C
C           Determine the dimension of the last block.
C
            IB = 1
            IF( NLOW.LT.NSUP ) THEN
               IF( A(NSUP,NSUP-1).NE.ZERO ) IB = 2
            END IF
C
C           Compute G, the current last IB rows of Z'*B.
C
            NL = NSUP - IB + 1
            CALL DGEMM( 'Transpose', 'NoTranspose', IB, M, N, ONE,
     $                  Z(1,NL), LDZ, B, LDB, ZERO, DWORK(KG), IB )
C
C           Check the controllability for a simple block.
C
            IF( DLANGE( '1', IB, M, DWORK(KG), IB, DWORK(KW) )
     $          .LE. TOLERB ) THEN
C
C              Deflate the uncontrollable block and resume the
C              main loop.
C
               NSUP = NSUP - IB
               NUP = NUP + IB
               GO TO 20
            END IF
C
C           Test for termination with INFO = 3.
C
            IF( NAP.EQ.NP) THEN
               INFO = 3
C
C              Test for compatibility. Terminate if an attempt occurs
C              to place a complex conjugate pair on a 1x1 block.
C
            ELSE IF( IB.EQ.1 .AND. NPR.EQ.0 .AND. NLOW.EQ.NSUP ) THEN
               INFO = 4
            ELSE
C
C              Set the simple block flag.
C
               SIMPLB = .TRUE.
C
C              Form a 2-by-2 block if necessary from two 1-by-1 blocks.
C              Consider special case IB = 1, NPR = 1 and
C              NPR+NPC > NSUP-NLOW+1 to avoid incompatibility.
C
               IF( ( IB.EQ.1 .AND. NPR.EQ.0 ) .OR.
     $             ( IB.EQ.1 .AND. NPR.EQ.1 .AND. NSUP.GT.NLOW .AND.
     $               NPR+NPC.GT.NSUP-NLOW+1 ) ) THEN
                  IF( NSUP.GT.2 ) THEN
                     IF( A(NSUP-1,NSUP-2) .NE. ZERO ) THEN
C
C                       Interchange with the adjacent 2x2 block.
C
C                       Workspace needed: N.
C
                        CALL DLAEXC( .TRUE., N, A, LDA, Z, LDZ, NSUP-2,
     $                               2, 1, DWORK(KW), INFO )
                        IF( INFO .NE. 0 ) THEN
                           INFO = 2
                           RETURN
                        END IF
                     ELSE
C
C                       Form a non-simple block by extending the last
C                       block with a 1x1 block.
C
                        SIMPLB = .FALSE.
                     END IF
                  ELSE
                     SIMPLB = .FALSE.
                  END IF
                  IB = 2
               END IF
               NL = NSUP - IB + 1
C
C              Compute G, the current last IB rows of Z'*B.
C
               CALL DGEMM( 'Transpose', 'NoTranspose', IB, M, N, ONE,
     $                      Z(1,NL), LDZ, B, LDB, ZERO, DWORK(KG), IB )
C
C              Check the controllability for the current block.
C
               IF( DLANGE( '1', IB, M, DWORK(KG), IB, DWORK(KW) )
     $            .LE. TOLERB ) THEN
C
C                 Deflate the uncontrollable block and resume the
C                 main loop.
C
                  NSUP = NSUP - IB
                  NUP = NUP + IB
                  GO TO 20
               END IF
C
               IF( NAP+IB .GT. NP ) THEN
C
C                 No sufficient eigenvalues to be assigned.
C
                  INFO = 3
               ELSE
                  IF( IB .EQ. 1 ) THEN
C
C                    A 1-by-1 block.
C
C                    Assign the real eigenvalue nearest to A(NSUP,NSUP).
C
                     X = A(NSUP,NSUP)
                     CALL SB01BX( .TRUE., NPR, X, X, WR, X, S, P )
                     NPR  = NPR - 1
                     CEIG = .FALSE.
                  ELSE
C
C                    A 2-by-2 block.
C
                     IF( SIMPLB ) THEN
C
C                       Simple 2-by-2 block with complex eigenvalues.
C                       Compute the eigenvalues of the last block.
C
                        CALL MB03QY( N, NL, A, LDA, Z, LDZ, X, Y, INFO )
                        IF( NPC .GT. 1 ) THEN
                           CALL SB01BX( .FALSE., NPC, X, Y,
     $                                  WR(IPC), WI(IPC), S, P )
                           NPC  = NPC - 2
                           CEIG = .TRUE.
                        ELSE
C
C                          Choose the nearest two real eigenvalues.
C
                           CALL SB01BX( .TRUE., NPR, X, X, WR, X, S, P )
                           CALL SB01BX( .TRUE., NPR-1, X, X, WR, X,
     $                                  Y, P )
                           P = S * Y
                           S = S + Y
                           NPR = NPR - 2
                           CEIG = .FALSE.
                        END IF
                     ELSE
C
C                       Non-simple 2x2 block with real eigenvalues.
C                       Choose the nearest pair of complex eigenvalues.
C
                        X = ( A(NL,NL) + A(NSUP,NSUP) )/TWO
                        CALL SB01BX( .FALSE., NPC, X, ZERO, WR(IPC),
     $                               WI(IPC), S, P )
                        NPC = NPC - 2
                     END IF
                  END IF
C
C                 Form the IBxIB matrix A2 from the current diagonal
C                 block.
C
                  A2(1,1) = A(NL,NL)
                  IF( IB .GT. 1 ) THEN
                     A2(1,2) = A(NL,NSUP)
                     A2(2,1) = A(NSUP,NL)
                     A2(2,2) = A(NSUP,NSUP)
                  END IF
C
C                 Determine the M-by-IB feedback matrix FI which
C                 assigns the chosen IB eigenvalues for the pair (A2,G).
C
C                 Workspace needed: 5*M.
C
                  CALL SB01BY( IB, M, S, P, A2, DWORK(KG), DWORK(KFI),
     $                         TOLER, DWORK(KW), IERR )
                  IF( IERR .NE. 0 ) THEN
                     IF( IB.EQ.1 .OR. SIMPLB ) THEN
C
C                       The simple 1x1 block is uncontrollable.
C
                        NSUP = NSUP - IB
                        IF( CEIG ) THEN
                           NPC = NPC + IB
                        ELSE
                           NPR = NPR + IB
                        END IF
                        NUP  = NUP + IB
                     ELSE
C
C                       The non-simple 2x2 block is uncontrollable.
C                       Eliminate its uncontrollable part by using
C                       the information in elements FI(1,1) and F(1,2).
C
                        C = DWORK(KFI)
                        S = DWORK(KFI+IB)
C
C                       Apply the transformation to A and accumulate it
C                       in Z.
C
                        CALL DROT( N-NL+1, A(NL,NL), LDA,
     $                             A(NSUP,NL), LDA, C, S )
                        CALL DROT( N, A(1,NL), 1, A(1,NSUP), 1, C, S )
                        CALL DROT( N, Z(1,NL), 1, Z(1,NSUP), 1, C, S )
C
C                       Annihilate the subdiagonal element of the last
C                       block, redefine the upper limit for the bottom
C                       block and resume the main loop.
C
                        A(NSUP,NL) = ZERO
                        NSUP = NL
                        NUP  = NUP + 1
                        NPC  = NPC + 2
                     END IF
                  ELSE
C
C                    Successful assignment of IB eigenvalues.
C
C                    Update the feedback matrix F <-- F + [0 FI]*Z'.
C
                     CALL DGEMM( 'NoTranspose', 'Transpose', M, N,
     $                           IB, ONE, DWORK(KFI), M, Z(1,NL),
     $                           LDZ, ONE, F, LDF )
C
C                    Check for possible numerical instability.
C
                     IF( DLANGE( '1', M, IB, DWORK(KFI), M, DWORK(KW) )
     $                           .GT. RMAX ) IWARN = IWARN + 1
C
C                    Update the state matrix A <-- A + Z'*B*[0 FI].
C                    Workspace needed: 2*N+4*M.
C
                     CALL DGEMM( 'NoTranspose', 'NoTranspose', N, IB,
     $                           M, ONE, B, LDB, DWORK(KFI), M, ZERO,
     $                           DWORK(KW), N )
                     CALL DGEMM( 'Transpose', 'NoTranspose', NSUP,
     $                           IB, N, ONE, Z, LDZ, DWORK(KW), N,
     $                           ONE, A(1,NL), LDA )
C
C                    Try to split the 2x2 block.
C
                     IF( IB .EQ. 2 )
     $                  CALL MB03QY( N, NL, A, LDA, Z, LDZ, X, Y,
     $                               INFO )
                     NAP = NAP + IB
                     IF( NLOW+IB.LE.NSUP ) THEN
C
C                       Move the last block(s) to the leading
C                       position(s) of the bottom block.
C
                        NCUR1 = NSUP - IB
                        NMOVES = 1
                        IF( IB.EQ.2 .AND. A(NSUP,NSUP-1).EQ.ZERO ) THEN
                           IB = 1
                           NMOVES = 2
                        END IF
C
C                       WHILE (NMOVES > 0) DO
   30                   IF( NMOVES .GT. 0 ) THEN
                           NCUR = NCUR1
C
C                          WHILE (NCUR >= NLOW) DO
   40                      IF( NCUR .GE. NLOW ) THEN
C
C                             Loop for the last block positioning.
C
                              IB1 = 1
                              IF( NCUR.GT.NLOW ) THEN
                                 IF( A(NCUR,NCUR-1).NE.ZERO ) IB1 = 2
                              END IF
                              CALL DLAEXC( .TRUE., N, A, LDA, Z, LDZ,
     $                                     NCUR-IB1+1, IB1, IB,
     $                                     DWORK(KW), INFO )
                              IF( INFO .NE. 0 ) THEN
                                 INFO = 2
                                 RETURN
                              END IF
                              NCUR = NCUR - IB1
                              GO TO 40
                           END IF
C
C                          END WHILE 40
C
                           NMOVES = NMOVES - 1
                           NCUR1 = NCUR1 + 1
                           NLOW = NLOW + IB
                           GO TO 30
                        END IF
C
C                       END WHILE 30
C
                     ELSE
                        NLOW = NLOW + IB
                     END IF
                  END IF
               END IF
            END IF
            IF( INFO.EQ.0 ) GO TO 20
C
C        END WHILE 20
C
         END IF
C
         WRKOPT = MAX( WRKOPT, 5*M, 2*N + 4*M )
      END IF
C
C     Annihilate the elements below the first subdiagonal of A.
C
      IF( N .GT. 2)
     $   CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
      IF( NAP .GT. 0 ) THEN
C
C        Move the assigned eigenvalues in the first NAP positions of
C        WR and WI.
C
         K = IPC - NPR - 1
         IF( K .GT. 0 ) CALL DSWAP( K, WR(NPR+1), 1, WR, 1 )
         J = NAP - K
         IF( J .GT. 0 ) THEN
            CALL DSWAP( J, WR(IPC+NPC), 1, WR(K+1), 1 )
            CALL DSWAP( J, WI(IPC+NPC), 1, WI(K+1), 1 )
         END IF
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB01BD ***
      END
