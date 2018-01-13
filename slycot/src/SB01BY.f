      SUBROUTINE SB01BY( N, M, S, P, A, B, F, TOL, DWORK, INFO )
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
C     To solve an N-by-N pole placement problem for the simple cases
C     N = 1 or N = 2: given the N-by-N matrix A and N-by-M matrix B,
C     construct an M-by-N matrix F such that A + B*F has prescribed
C     eigenvalues. These eigenvalues are specified by their sum S and
C     product P (if N = 2). The resulting F has minimum Frobenius norm.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A and also the number of rows of
C             the matrix B and the number of columns of the matrix F.
C             N is either 1, if a single real eigenvalue is prescribed
C             or 2, if a complex conjugate pair or a set of two real
C             eigenvalues are prescribed.
C
C     M       (input) INTEGER
C             The number of columns of the matrix B and also the number
C             of rows of the matrix F.  M >= 1.
C
C     S       (input) DOUBLE PRECISION
C             The sum of the prescribed eigenvalues if N = 2 or the
C             value of prescribed eigenvalue if N = 1.
C
C     P       (input) DOUBLE PRECISION
C             The product of the prescribed eigenvalues if N = 2.
C             Not referenced if N = 1.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N,N)
C             On entry, this array must contain the N-by-N state
C             dynamics matrix whose eigenvalues have to be moved to
C             prescribed locations.
C             On exit, this array contains no useful information.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N,M)
C             On entry, this array must contain the N-by-M input/state
C             matrix B.
C             On exit, this array contains no useful information.
C
C     F       (output) DOUBLE PRECISION array, dimension (M,N)
C             The state feedback matrix F which assigns one pole or two
C             poles of the closed-loop matrix A + B*F.
C             If N = 2 and the pair (A,B) is not controllable
C             (INFO = 1), then F(1,1) and F(1,2) contain the elements of
C             an orthogonal rotation which can be used to remove the
C             uncontrollable part of the pair (A,B).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The absolute tolerance level below which the elements of A
C             and B are considered zero (used for controllability test).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (M)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if uncontrollability of the pair (A,B) is detected.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine SB01BY.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     May  2003, A. Varga, German Aerospace Center.
C
C     KEYWORDS
C
C     Eigenvalue, eigenvalue assignment, feedback control, pole
C     placement, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  FOUR, ONE, THREE, TWO, ZERO
      PARAMETER         ( FOUR = 4.0D0,  ONE = 1.0D0, THREE = 3.0D0,
     $                    TWO  = 2.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, M, N
      DOUBLE PRECISION  P, S, TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,*), B(N,*), DWORK(*), F(M,*)
C     .. Local Scalars ..
      INTEGER           IR, J
      DOUBLE PRECISION  ABSR, B1, B2, B21, C, C0, C1, C11, C12, C21,
     $                  C22, C3, C4, CS, CU, CV, DC0, DC2, DC3, DIFFR,
     $                  R, RN, S12, S21, SIG, SN, SU, SV, TAU1, TAU2,
     $                  WI, WI1, WR, WR1, X, Y, Z
C     .. External Functions ..
      DOUBLE PRECISION  DLAMC3, DLAMCH
      EXTERNAL          DLAMC3, DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DLANV2, DLARFG, DLASET, DLASV2, SLCT_DLATZM, DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
C
C     For efficiency reasons, the parameters are not checked.
C
      INFO = 0
      IF( N.EQ.1 ) THEN
C
C        The case N = 1.
C
         IF( M.GT.1 )
     $      CALL DLARFG( M, B(1,1), B(1,2), N, TAU1 )
         B1 = B(1,1)
         IF( ABS( B1 ).LE.TOL ) THEN
C
C           The pair (A,B) is uncontrollable.
C
            INFO = 1
            RETURN
         END IF
C
         F(1,1) = ( S - A(1,1) )/B1
         IF( M.GT.1 ) THEN
            CALL DLASET( 'Full', M-1, 1, ZERO, ZERO, F(2,1), M )
            CALL SLCT_DLATZM( 'Left', M, N, B(1,2), N, TAU1,
     $                   F(1,1), F(2,1),
     $                   M, DWORK )
         END IF
         RETURN
      END IF
C
C     In the sequel N = 2.
C
C     Compute the singular value decomposition of B in the form
C
C                    ( V  0 )                ( B1 0  )
C     B = U*( G1 0 )*(      )*H2*H1 ,   G1 = (       ),
C                    ( 0  I )                ( 0  B2 )
C
C               ( CU   SU )          ( CV   SV )
C     where U = (         )  and V = (         )  are orthogonal
C               (-SU   CU )          (-SV   CV )
C
C     rotations and H1 and H2 are elementary Householder reflectors.
C     ABS(B1) and ABS(B2) are the singular values of matrix B,
C     with ABS(B1) >= ABS(B2).
C
C     Reduce first B to the lower bidiagonal form  ( B1  0  ... 0 ).
C                                                  ( B21 B2 ... 0 )
      IF( M.EQ.1 ) THEN
C
C        Initialization for the case M = 1; no reduction required.
C
         B1  = B(1,1)
         B21 = B(2,1)
         B2  = ZERO
      ELSE
C
C        Postmultiply B with elementary Householder reflectors H1
C        and H2.
C
         CALL DLARFG( M, B(1,1), B(1,2), N, TAU1 )
         CALL SLCT_DLATZM( 'Right', N-1, M, B(1,2), N, TAU1,
     $                B(2,1), B(2,2),
     $                N, DWORK )
         B1  = B(1,1)
         B21 = B(2,1)
         IF( M.GT.2 )
     $      CALL DLARFG( M-1, B(2,2), B(2,3), N, TAU2 )
         B2  = B(2,2)
      END IF
C
C     Reduce B to a diagonal form by premultiplying and postmultiplying
C     it with orthogonal rotations U and V, respectively, and order the
C     diagonal elements to have decreasing magnitudes.
C     Note: B2 has been set to zero if M = 1. Thus in the following
C     computations the case M = 1 need not to be distinguished.
C     Note also that LAPACK routine DLASV2 assumes an upper triangular
C     matrix, so the results should be adapted.
C
      CALL DLASV2( B1, B21, B2, X, Y, SU, CU, SV, CV )
      SU = -SU
      B1 =  Y
      B2 =  X
C
C     Compute  A1 = U'*A*U.
C
      CALL DROT( 2, A(2,1), 2, A(1,1), 2, CU, SU )
      CALL DROT( 2, A(1,2), 1, A(1,1), 1, CU, SU )
C
C     Compute the rank of B and check the controllability of the
C     pair (A,B).
C
      IR = 0
      IF( ABS( B2 ).GT.TOL ) IR = IR + 1
      IF( ABS( B1 ).GT.TOL ) IR = IR + 1
      IF( IR.EQ.0 .OR. ( IR.EQ.1 .AND. ABS( A(2,1) ).LE.TOL ) ) THEN
         F(1,1) =  CU
         F(1,2) = -SU
C
C        The pair (A,B) is uncontrollable.
C
         INFO = 1
         RETURN
      END IF
C
C     Compute F1 which assigns N poles for the reduced pair (A1,G1).
C
      X = DLAMC3( B1, B2 )
      IF( X.EQ.B1 ) THEN
C
C        Rank one G1.
C
         F(1,1) = ( S - ( A(1,1) + A(2,2) ) )/B1
         F(1,2) = -( A(2,2)*( A(2,2) - S ) + A(2,1)*A(1,2) + P )/
     $             A(2,1)/B1
         IF( M.GT.1 ) THEN
            F(2,1) = ZERO
            F(2,2) = ZERO
         END IF
      ELSE
C
C        Rank two G1.
C
         Z = ( S - ( A(1,1) + A(2,2) ) )/( B1*B1 + B2*B2 )
         F(1,1) = B1*Z
         F(2,2) = B2*Z
C
C        Compute an approximation for the minimum norm parameter
C        selection.
C
         X = A(1,1) + B1*F(1,1)
         C = X*( S - X ) - P
         IF( C.GE.ZERO ) THEN
            SIG =  ONE
         ELSE
            SIG = -ONE
         END IF
         S12 = B1/B2
         S21 = B2/B1
         C11 = ZERO
         C12 = ONE
         C21 = SIG*S12*C
         C22 = A(1,2) - SIG*S12*A(2,1)
         CALL DLANV2( C11, C12, C21, C22, WR, WI, WR1, WI1, CS, SN )
         IF( ABS( WR - A(1,2) ).GT.ABS( WR1 - A(1,2) ) ) THEN
            R = WR1
         ELSE
            R = WR
         END IF
C
C        Perform Newton iteration to solve the equation for minimum.
C
         C0 = -C*C
         C1 =  C*A(2,1)
         C4 =  S21*S21
         C3 = -C4*A(1,2)
         DC0 = C1
         DC2 = THREE*C3
         DC3 = FOUR*C4
C
         DO 10 J = 1, 10
            X  = C0 + R*( C1 + R*R*( C3 + R*C4 ) )
            Y  = DC0 + R*R*( DC2 + R*DC3 )
            IF( Y.EQ.ZERO ) GO TO 20
            RN = R - X/Y
            ABSR  = ABS( R )
            DIFFR = ABS( R - RN )
            Z = DLAMC3( ABSR, DIFFR )
            IF( Z.EQ.ABSR )
     $         GO TO 20
            R = RN
   10    CONTINUE
C
   20    CONTINUE
         IF( R.EQ.ZERO ) R = DLAMCH( 'Epsilon' )
         F(1,2) = (  R  - A(1,2) )/B1
         F(2,1) = ( C/R - A(2,1) )/B2
      END IF
C
C     Back-transform F1. Compute first F1*U'.
C
      CALL DROT( MIN( M, 2 ), F(1,1), 1, F(1,2), 1, CU, SU )
      IF( M.EQ.1 )
     $   RETURN
C
C     Compute V'*F1.
C
      CALL DROT( 2, F(2,1), M, F(1,1), M, CV, SV )
C
C               ( F1 )
C     Form  F = (    ) .
C               ( 0  )
C
      IF( M.GT.N )
     $   CALL DLASET( 'Full', M-N, N, ZERO, ZERO, F(N+1,1), M )
C
C     Compute H1*H2*F.
C
      IF( M.GT.2 )
     $     CALL SLCT_DLATZM( 'Left', M-1, N, B(2,3), N, TAU2,
     $                F(2,1), F(3,1), M, DWORK )
      CALL SLCT_DLATZM( 'Left', M, N, B(1,2), N, TAU1,
     $                  F(1,1), F(2,1), M, DWORK )
C
      RETURN
C *** Last line of SB01BY ***
      END
