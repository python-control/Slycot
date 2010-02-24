      SUBROUTINE TB01VY( APPLY, N, M, L, THETA, LTHETA, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, X0, DWORK, LDWORK, INFO )
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
C     To convert the linear discrete-time system given as its output
C     normal form [1], with parameter vector THETA, into the state-space
C     representation (A, B, C, D), with the initial state x0.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     APPLY   CHARACTER*1
C             Specifies whether or not the parameter vector should be
C             transformed using a bijective mapping, as follows:
C             = 'A' : apply the bijective mapping to the N vectors in
C                     THETA corresponding to the matrices A and C;
C             = 'N' : do not apply the bijective mapping.
C             The transformation performed when APPLY = 'A' allows
C             to get rid of the constraints norm(THETAi) < 1, i = 1:N.
C             A call of the SLICOT Library routine TB01VD associated to
C             a call of TB01VY must use the same value of APPLY.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     L       (input) INTEGER
C             The number of system outputs.  L >= 0.
C
C     THETA   (input) DOUBLE PRECISION array, dimension (LTHETA)
C             The leading N*(L+M+1)+L*M part of this array must contain
C             the parameter vector that defines a system (A, B, C, D),
C             with the initial state x0. The parameters are:
C
C             THETA(1:N*L)                      : parameters for A, C;
C             THETA(N*L+1:N*(L+M))              : parameters for B;
C             THETA(N*(L+M)+1:N*(L+M)+L*M)      : parameters for D;
C             THETA(N*(L+M)+L*M+1:N*(L+M+1)+L*M): parameters for x0.
C
C     LTHETA  INTEGER
C             The length of array THETA.  LTHETA >= N*(L+M+1)+L*M.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the system
C             state matrix corresponding to the output normal form with
C             parameter vector THETA.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array contains the system
C             input matrix corresponding to the output normal form with
C             parameter vector THETA.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading L-by-N part of this array contains the system
C             output matrix corresponding to the output normal form with
C             parameter vector THETA.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,L).
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading L-by-M part of this array contains the system
C             input/output matrix corresponding to the output normal
C             form with parameter vector THETA.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,L).
C
C     X0      (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the initial state of the system, x0,
C             corresponding to the output normal form with parameter
C             vector THETA.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= N*(N+L+1).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The parameters characterizing A and C are used to build N
C     orthogonal transformations, which are then applied to recover
C     these matrices.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Feb. 2002, Feb. 2004.
C
C     KEYWORDS
C
C     Asymptotically stable, output normal form, parameter estimation,
C     similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )
C     .. Scalar Arguments ..
      CHARACTER         APPLY
      INTEGER           INFO, L, LDA, LDB, LDC, LDD, LDWORK, LTHETA, M,
     $                  N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), THETA(*), X0(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR, RI, TI, TOBYPI
      INTEGER           CA, JWORK, I, IN, J, K, LDCA
      LOGICAL           LAPPLY
C     .. External Functions ..
      EXTERNAL          DNRM2, LSAME
      DOUBLE PRECISION  DNRM2
      LOGICAL           LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DGER, DLACPY, DSCAL,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      LAPPLY = LSAME( APPLY, 'A' )
C
      INFO = 0
      IF ( .NOT.( LAPPLY .OR. LSAME( APPLY, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( N.LT.0 ) THEN
         INFO = -2
      ELSEIF ( M.LT.0 ) THEN
         INFO = -3
      ELSEIF ( L.LT.0 ) THEN
         INFO = -4
      ELSEIF ( LTHETA.LT.( N*( L + M + 1 ) + L*M ) ) THEN
         INFO = -6
      ELSEIF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSEIF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSEIF ( LDC.LT.MAX( 1, L ) ) THEN
         INFO = -12
      ELSEIF ( LDD.LT.MAX( 1, L ) ) THEN
         INFO = -14
      ELSEIF ( LDWORK.LT.N*( N + L + 1 ) ) THEN
         INFO = -17
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TB01VY', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, L ).EQ.0 )
     $   RETURN
C
      IF ( M.GT.0 ) THEN
C
C        Copy the matrix B from THETA.
C
         CALL DLACPY( 'Full', N, M, THETA(N*L+1), N, B, LDB )
C
C        Copy the matrix D.
C
         CALL DLACPY( 'Full', L, M, THETA(N*(L+M)+1), L, D, LDD )
      ENDIF
C
      IF ( N.EQ.0 ) THEN
         RETURN
      ELSE IF ( L.EQ.0 ) THEN
         CALL DCOPY( N, THETA(N*M+1), 1, X0, 1 )
         RETURN
      END IF
C
C     Initialize the indices in the workspace.
C
      LDCA = N + L
C
      CA = 1
C
      JWORK  = CA + N*LDCA
      TOBYPI = HALF/ATAN( ONE )
C
C     Generate the matrices C and A from their parameters.
C     Start with the block matrix [0; I], where 0 is a block of zeros
C     of size L-by-N, and I is the identity matrix of order N.
C
      DWORK(CA) = ZERO
      CALL DCOPY( N*(L+N), DWORK(CA), 0, DWORK(CA), 1 )
      DWORK(CA+L) = ONE
      CALL DCOPY( N, DWORK(CA+L), 0, DWORK(CA+L), LDCA+1 )
C
C     Now, read out THETA(1 : N*L) and perform the transformations
C     defined by the parameters in THETA.
C
      DO 30 I = N, 1, -1
C
C        Save THETAi in the first column of C and use the copy for
C        further processing.
C
         CALL DCOPY( L, THETA((I-1)*L+1), 1, C, 1 )
         TI = DNRM2( L, C, 1 )
         IF ( LAPPLY .AND. TI.NE.ZERO ) THEN
C
C           Apply the bijective mapping which guarantees that TI < 1.
C
            FACTOR = TOBYPI*ATAN( TI )/TI
C
C           Scale THETAi and apply the same scaling on TI.
C
            CALL DSCAL( L, FACTOR, C, 1 )
            TI = TI*FACTOR
         END IF
C
C        RI = sqrt( 1 - TI**2 ).
C
         RI = SQRT( ( ONE - TI )*( ONE + TI ) )
C
C        Multiply a certain part of DWORK(CA) with Ui' from the left,
C        where Ui = [ -THETAi, Si; RI, THETAi' ] is (L+1)-by-(L+1), but
C        Ui is not stored.
C
         CALL DGEMV( 'Transpose', L, N, -ONE, DWORK(CA+N-I), LDCA, C, 1,
     $               ZERO, DWORK(JWORK), 1 )
C
         IF ( TI.GT.ZERO ) THEN
            CALL DGER( L, N, (ONE-RI)/TI/TI, C, 1, DWORK(JWORK), 1,
     $                 DWORK(CA+N-I), LDCA )
         ELSE
C
C           The call below is for the limiting case.
C
            CALL DGER( L, N, HALF, C, 1, DWORK(JWORK), 1,
     $                 DWORK(CA+N-I), LDCA )
         ENDIF
C
         CALL DGER( L, N, ONE, C, 1, DWORK(CA+N-I+L), LDCA,
     $              DWORK(CA+N-I), LDCA )
         CALL DAXPY( N, RI, DWORK(CA+N-I+L), LDCA, DWORK(JWORK), 1 )
C
C        Move these results to their appropriate locations.
C
         DO 20 J = 1, N
            IN = CA + N - I + ( J - 1 )*LDCA
            DO 10 K = IN + L, IN + 1, -1
               DWORK(K) = DWORK(K-1)
   10       CONTINUE
            DWORK(IN) = DWORK(JWORK+J-1)
   20    CONTINUE
C
   30 CONTINUE
C
C     Now, DWORK(CA) = [C; A]. Copy to C and A.
C
      DO 40 I = 1, N
         CALL DCOPY( L, DWORK(CA+(I-1)*LDCA),   1, C(1,I), 1 )
         CALL DCOPY( N, DWORK(CA+L+(I-1)*LDCA), 1, A(1,I), 1 )
   40 CONTINUE
C
C     Copy the initial state x0.
C
      CALL DCOPY( N, THETA(N*(L+M)+L*M+1), 1, X0, 1 )
C
      RETURN
C
C *** Last line of TB01VY ***
      END
