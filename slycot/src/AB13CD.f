      DOUBLE PRECISION FUNCTION AB13CD( N, M, NP, A, LDA, B, LDB, C,
     $                                  LDC, D, LDD, TOL, IWORK, DWORK,
     $                                  LDWORK, CWORK, LCWORK, BWORK,
     $                                  INFO )
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
C     To compute the H-infinity norm of the continuous-time stable
C     system
C
C                          | A | B |
C                   G(s) = |---|---| .
C                          | C | D |
C
C     FUNCTION VALUE
C
C     AB13CD  DOUBLE PRECISION
C             If INFO = 0, the H-infinity norm of the system, HNORM,
C             i.e., the peak gain of the frequency response (as measured
C             by the largest singular value in the MIMO case).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     NP      (input) INTEGER
C             The row size of the matrix C.  NP >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading NP-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used to set the accuracy in determining the
C             norm.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension N
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK, and DWORK(2) contains the frequency where the
C             gain of the frequency response achieves its peak value
C             HNORM.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= max(2,4*N*N+2*M*M+3*M*N+M*NP+2*(N+NP)*NP+10*N+
C                             6*max(M,NP)).
C             For good performance, LDWORK must generally be larger.
C
C     CWORK   COMPLEX*16 array, dimension (LCWORK)
C             On exit, if INFO = 0, CWORK(1) contains the optimal value
C             of LCWORK.
C
C     LCWORK  INTEGER
C             The dimension of the array CWORK.
C             LCWORK >= max(1,(N+M)*(N+NP)+3*max(M,NP)).
C             For good performance, LCWORK must generally be larger.
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the system is unstable;
C             = 2:  the tolerance is too small (the algorithm for
C                   computing the H-infinity norm did not converge);
C             = 3:  errors in computing the eigenvalues of A or of the
C                   Hamiltonian matrix (the QR algorithm did not
C                   converge);
C             = 4:  errors in computing singular values.
C
C     METHOD
C
C     The routine implements the method presented in [1].
C
C     REFERENCES
C
C     [1] Bruinsma, N.A. and Steinbuch, M.
C         A fast algorithm to compute the Hinfinity-norm of a transfer
C         function matrix.
C         Systems & Control Letters, vol. 14, pp. 287-293, 1990.
C
C     NUMERICAL ASPECTS
C
C     If the algorithm does not converge (INFO = 2), the tolerance must
C     be increased.
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999,
C     Oct. 2000.
C     P.Hr. Petkov, October 2000.
C     A. Varga, October 2000.
C     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     H-infinity optimal control, robust control, system norm.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 10 )
      COMPLEX*16         CONE, JIMAG
      PARAMETER          ( CONE  = ( 1.0D0, 0.0D0 ),
     $                     JIMAG = ( 0.0D0, 1.0D0 ) )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION   HUGE
      PARAMETER          ( HUGE = 10.0D+0**30 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LDC, LCWORK, LDD, LDWORK, M, N,
     $                   NP
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      COMPLEX*16         CWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), DWORK( * )
      LOGICAL            BWORK( * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, ICW2, ICW3, ICW4, ICWRK, INFO2, ITER, IW10,
     $                   IW11, IW12, IW2, IW3, IW4, IW5, IW6, IW7, IW8,
     $                   IW9, IWRK, J, K, L, LCWAMX, LWAMAX, MINCWR,
     $                   MINWRK, SDIM
      DOUBLE PRECISION   DEN, FPEAK, GAMMA, GAMMAL, GAMMAU, OMEGA, RAT,
     $                   RATMAX, TEMP, WIMAX, WRMIN
      LOGICAL            COMPLX
C
C     .. External Functions ..
      DOUBLE PRECISION   DLAPY2
      LOGICAL            SB02MV, SB02CX
      EXTERNAL           DLAPY2, SB02MV, SB02CX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEES, DGEMM, DGESV, DGESVD, DLACPY, DPOSV,
     $                   DPOTRF, DPOTRS, DSYRK, MA02ED, MB01RX, XERBLA,
     $                   ZGEMM, ZGESV, ZGESVD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( NP.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, NP ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, NP ) ) THEN
         INFO = -11
      END IF
C
C     Compute workspace.
C
      MINWRK =  MAX( 2, 4*N*N + 2*M*M + 3*M*N + M*NP + 2*( N + NP )*NP +
     $                  10*N + 6*MAX( M, NP ) )
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -15
      END IF
      MINCWR = MAX( 1, ( N + M )*( N + NP ) + 3*MAX( M, NP ) )
      IF( LCWORK.LT.MINCWR ) THEN
         INFO = -17
      END IF
      IF( INFO.NE.0 ) THEN
         AB13CD   = ZERO
         CALL XERBLA( 'AB13CD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( M.EQ.0 .OR. NP.EQ.0 ) THEN
         AB13CD   = ZERO
         RETURN
      END IF
C
C     Workspace usage.
C
      IW2  = N
      IW3  = IW2 + N
      IW4  = IW3 + N*N
      IW5  = IW4 + N*M
      IW6  = IW5 + NP*M
      IWRK = IW6 + MIN( NP, M )
C
C     Determine the maximum singular value of G(infinity) = D .
C
      CALL DLACPY( 'Full', NP, M, D, LDD, DWORK( IW5+1 ), NP )
      CALL DGESVD( 'N', 'N', NP, M, DWORK( IW5+1 ), NP, DWORK( IW6+1 ),
     $             DWORK, NP, DWORK, M, DWORK( IWRK+1 ), LDWORK-IWRK,
     $             INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 4
         AB13CD   = ZERO
         RETURN
      END IF
      GAMMAL = DWORK( IW6+1 )
      FPEAK  = HUGE
      LWAMAX = INT( DWORK( IWRK+1 ) ) + IWRK
C
C     Quick return if N = 0 .
C
      IF( N.EQ.0 ) THEN
         AB13CD   = GAMMAL
         DWORK(1) = TWO
         DWORK(2) = ZERO
         CWORK(1) = ONE
         RETURN
      END IF
C
C     Stability check.
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IW3+1 ), N )
      CALL DGEES( 'N', 'S', SB02MV, N, DWORK( IW3+1 ), N, SDIM, DWORK,
     $            DWORK( IW2+1 ), DWORK, N, DWORK( IWRK+1 ),
     $            LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 3
         RETURN
      END IF
      IF( SDIM.LT.N ) THEN
         INFO = 1
         RETURN
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK+1 ) ) + IWRK, LWAMAX )
C
C     Determine the maximum singular value of G(0) = -C*inv(A)*B + D .
C
      CALL DLACPY( 'Full', N,  N, A, LDA, DWORK( IW3+1 ), N )
      CALL DLACPY( 'Full', N,  M, B, LDB, DWORK( IW4+1 ), N )
      CALL DLACPY( 'Full', NP, M, D, LDD, DWORK( IW5+1 ), NP )
      CALL DGESV( N, M, DWORK( IW3+1 ), N, IWORK, DWORK( IW4+1 ), N,
     $            INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      CALL DGEMM( 'N', 'N', NP, M, N, -ONE, C, LDC, DWORK( IW4+1 ), N,
     $            ONE, DWORK( IW5+1 ), NP )
      CALL DGESVD( 'N', 'N', NP, M, DWORK( IW5+1 ), NP, DWORK( IW6+1 ),
     $             DWORK, NP, DWORK, M, DWORK( IWRK+1 ), LDWORK-IWRK,
     $             INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 4
         RETURN
      END IF
      IF( GAMMAL.LT.DWORK( IW6+1 ) ) THEN
         GAMMAL = DWORK( IW6+1 )
         FPEAK  = ZERO
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK+1 ) ) + IWRK, LWAMAX )
C
C     Find a frequency which is close to the peak frequency.
C
      COMPLX = .FALSE.
      DO 10 I = 1, N
         IF( DWORK( IW2+I ).NE.ZERO ) COMPLX = .TRUE.
   10 CONTINUE
      IF( .NOT.COMPLX ) THEN
         WRMIN = ABS( DWORK( 1 ) )
         DO 20 I = 2, N
            IF( WRMIN.GT.ABS( DWORK( I ) ) ) WRMIN = ABS( DWORK( I ) )
   20    CONTINUE
         OMEGA = WRMIN
      ELSE
         RATMAX = ZERO
         DO 30 I = 1, N
            DEN = DLAPY2( DWORK( I ), DWORK( IW2+I ) )
            RAT = ABS( ( DWORK( IW2+I )/DWORK( I ) )/DEN )
            IF( RATMAX.LT.RAT ) THEN
               RATMAX = RAT
               WIMAX  = DEN
            END IF
   30    CONTINUE
         OMEGA = WIMAX
      END IF
C
C     Workspace usage.
C
      ICW2  = N*N
      ICW3  = ICW2 + N*M
      ICW4  = ICW3 + NP*N
      ICWRK = ICW4 + NP*M
C
C     Determine the maximum singular value of
C     G(omega) = C*inv(j*omega*In - A)*B + D .
C
      DO 50 J = 1, N
         DO 40 I = 1, N
            CWORK( I+(J-1)*N ) = -A( I, J )
   40    CONTINUE
         CWORK( J+(J-1)*N ) = JIMAG*OMEGA - A( J, J )
   50 CONTINUE
      DO 70 J = 1, M
         DO 60 I = 1, N
            CWORK( ICW2+I+(J-1)*N ) = B( I, J )
   60    CONTINUE
   70 CONTINUE
      DO 90 J = 1, N
         DO 80 I = 1, NP
            CWORK( ICW3+I+(J-1)*NP ) = C( I, J )
   80    CONTINUE
   90 CONTINUE
      DO 110 J = 1, M
         DO 100 I = 1, NP
            CWORK( ICW4+I+(J-1)*NP ) = D( I, J )
  100    CONTINUE
  110 CONTINUE
      CALL ZGESV( N, M, CWORK, N, IWORK, CWORK( ICW2+1 ), N, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
      CALL ZGEMM( 'N', 'N', NP, M, N, CONE, CWORK( ICW3+1 ), NP,
     $            CWORK( ICW2+1 ), N, CONE, CWORK( ICW4+1 ), NP )
      CALL ZGESVD( 'N', 'N', NP, M, CWORK( ICW4+1 ), NP, DWORK( IW6+1 ),
     $             CWORK, NP, CWORK, M, CWORK( ICWRK+1 ), LCWORK-ICWRK,
     $             DWORK( IWRK+1 ), INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 4
         RETURN
      END IF
      IF( GAMMAL.LT.DWORK( IW6+1 ) ) THEN
         GAMMAL = DWORK( IW6+1 )
         FPEAK  = OMEGA
      END IF
      LCWAMX = INT( CWORK( ICWRK+1 ) ) + ICWRK
C
C     Workspace usage.
C
      IW2  = M*N
      IW3  = IW2  + M*M
      IW4  = IW3  + NP*NP
      IW5  = IW4  + M*M
      IW6  = IW5  + M*N
      IW7  = IW6  + M*N
      IW8  = IW7  + NP*NP
      IW9  = IW8  + NP*N
      IW10 = IW9  + 4*N*N
      IW11 = IW10 + 2*N
      IW12 = IW11 + 2*N
      IWRK = IW12 + MIN( NP, M )
C
C     Compute D'*C .
C
      CALL DGEMM( 'T', 'N', M, N, NP, ONE, D, LDD, C, LDC, ZERO,
     $            DWORK, M )
C
C     Compute D'*D .
C
      CALL DSYRK( 'U', 'T', M, NP, ONE, D, LDD, ZERO, DWORK( IW2+1 ),
     $            M )
C
C     Compute D*D' .
C
      CALL DSYRK( 'U', 'N', NP, M, ONE, D, LDD, ZERO, DWORK( IW3+1 ),
     $            NP )
C
C     Main iteration loop for gamma.
C
      ITER = 0
  120 ITER = ITER + 1
      IF( ITER.GT.MAXIT ) THEN
         INFO = 2
         RETURN
      END IF
      GAMMA = ( ONE + TWO*TOL )*GAMMAL
C
C     Compute R = GAMMA^2*Im - D'*D .
C
      DO 140 J = 1, M
         DO 130 I = 1, J
            DWORK( IW4+I+(J-1)*M ) = -DWORK( IW2+I+(J-1)*M )
  130    CONTINUE
         DWORK( IW4+J+(J-1)*M ) = GAMMA**2 - DWORK( IW2+J+(J-1)*M )
  140 CONTINUE
C
C     Compute inv(R)*D'*C .
C
      CALL DLACPY( 'Full', M, N, DWORK, M, DWORK( IW5+1 ), M )
      CALL DPOTRF( 'U', M, DWORK( IW4+1 ), M, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
      CALL DPOTRS( 'U', M, N, DWORK( IW4+1 ), M, DWORK( IW5+1 ), M,
     $             INFO2 )
C
C     Compute inv(R)*B' .
C
      DO 160 J = 1, N
         DO 150 I = 1, M
            DWORK( IW6+I+(J-1)*M ) = B( J, I )
  150    CONTINUE
  160 CONTINUE
      CALL DPOTRS( 'U', M, N, DWORK( IW4+1 ), M, DWORK( IW6+1 ), M,
     $             INFO2 )
C
C     Compute S = GAMMA^2*Ip - D*D' .
C
      DO 180 J = 1, NP
         DO 170 I = 1, J
            DWORK( IW7+I+(J-1)*NP ) = -DWORK( IW3+I+(J-1)*NP )
  170    CONTINUE
         DWORK( IW7+J+(J-1)*NP ) = GAMMA**2 - DWORK( IW3+J+(J-1)*NP )
  180 CONTINUE
C
C     Compute inv(S)*C .
C
      CALL DLACPY( 'Full', NP, N, C, LDC, DWORK( IW8+1 ), NP )
      CALL DPOSV( 'U', NP, N, DWORK( IW7+1 ), NP, DWORK( IW8+1 ), NP,
     $            INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
C
C     Construct the Hamiltonian matrix .
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IW9+1 ), 2*N )
      CALL DGEMM( 'N', 'N', N, N, M, ONE, B, LDB, DWORK( IW5+1 ), M,
     $            ONE, DWORK( IW9+1 ), 2*N )
      CALL MB01RX( 'Left', 'Upper', 'Transpose', N, NP, ZERO, -GAMMA,
     $             DWORK( IW9+N+1 ), 2*N, C, LDC, DWORK( IW8+1 ), NP,
     $             INFO2 )
      CALL MA02ED( 'Upper', N, DWORK( IW9+N+1 ), 2*N )
      CALL MB01RX( 'Left', 'Upper', 'NoTranspose', N, M, ZERO, GAMMA,
     $             DWORK( IW9+2*N*N+1 ), 2*N, B, LDB, DWORK( IW6+1 ), M,
     $             INFO2 )
      CALL MA02ED( 'Upper', N, DWORK( IW9+2*N*N+1 ), 2*N )
      DO 200 J = 1, N
         DO 190 I = 1, N
            DWORK( IW9+2*N*N+N+I+(J-1)*2*N ) = -DWORK( IW9+J+(I-1)*2*N )
  190    CONTINUE
  200 CONTINUE
C
C     Compute the eigenvalues of the Hamiltonian matrix.
C
      CALL DGEES( 'N', 'S', SB02CX, 2*N, DWORK( IW9+1 ), 2*N, SDIM,
     $            DWORK( IW10+1 ), DWORK( IW11+1 ), DWORK, 2*N,
     $            DWORK( IWRK+1 ), LDWORK-IWRK, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 3
         RETURN
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK+1 ) ) + IWRK, LWAMAX )
C
      IF( SDIM.EQ.0 ) THEN
         GAMMAU = GAMMA
         GO TO 330
      END IF
C
C     Store the positive imaginary parts.
C
      J = 0
      DO 210 I = 1, SDIM-1, 2
         J = J + 1
         DWORK( IW10+J ) = DWORK( IW11+I )
  210 CONTINUE
      K = J
C
      IF( K.GE.2 ) THEN
C
C        Reorder the imaginary parts.
C
         DO 230 J = 1, K-1
            DO 220 L = J+1, K
               IF( DWORK( IW10+J ).LE. DWORK( IW10+L ) ) GO TO 220
                  TEMP = DWORK( IW10+J )
                  DWORK( IW10+J ) = DWORK( IW10+L )
                  DWORK( IW10+L ) = TEMP
  220       CONTINUE
  230    CONTINUE
C
C        Determine the next frequency.
C
         DO 320 L = 1, K - 1
            OMEGA = ( DWORK( IW10+L ) + DWORK( IW10+L+1 ) )/TWO
            DO 250 J = 1, N
               DO 240 I = 1, N
                  CWORK( I+(J-1)*N ) = -A( I, J )
  240          CONTINUE
               CWORK( J+(J-1)*N ) = JIMAG*OMEGA - A( J, J )
  250       CONTINUE
            DO 270 J = 1, M
               DO 260 I = 1, N
                  CWORK( ICW2+I+(J-1)*N ) = B( I, J )
  260          CONTINUE
  270       CONTINUE
            DO 290 J = 1, N
               DO 280 I = 1, NP
                  CWORK( ICW3+I+(J-1)*NP ) = C( I, J )
  280          CONTINUE
  290       CONTINUE
            DO 310 J = 1, M
               DO 300 I = 1, NP
                  CWORK( ICW4+I+(J-1)*NP ) = D( I, J )
  300          CONTINUE
  310       CONTINUE
            CALL ZGESV( N, M, CWORK, N, IWORK, CWORK( ICW2+1 ), N,
     $                  INFO2 )
            IF( INFO2.GT.0 ) THEN
               INFO = 1
               RETURN
            END IF
            CALL ZGEMM( 'N', 'N', NP, M, N, CONE, CWORK( ICW3+1 ), NP,
     $                   CWORK( ICW2+1 ), N, CONE, CWORK( ICW4+1 ), NP )
            CALL ZGESVD( 'N', 'N', NP, M, CWORK( ICW4+1 ), NP,
     $                   DWORK( IW6+1 ), CWORK, NP, CWORK, M,
     $                   CWORK( ICWRK+1 ), LCWORK-ICWRK,
     $                   DWORK( IWRK+1 ), INFO2 )
            IF( INFO2.GT.0 ) THEN
               INFO = 4
               RETURN
            END IF
            IF( GAMMAL.LT.DWORK( IW6+1 ) ) THEN
               GAMMAL = DWORK( IW6+1 )
               FPEAK  = OMEGA
            END IF
            LCWAMX = MAX( INT( CWORK( ICWRK+1 ) ) + ICWRK, LCWAMX )
  320    CONTINUE
      END IF
      GO TO 120
  330 AB13CD = ( GAMMAL + GAMMAU )/TWO
C
      DWORK( 1 ) = LWAMAX
      DWORK( 2 ) = FPEAK
      CWORK( 1 ) = LCWAMX
      RETURN
C *** End of AB13CD ***
      END
