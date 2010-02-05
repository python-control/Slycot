      SUBROUTINE TB04AY( N, MWORK, PWORK, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, NCONT, INDEXD, DCOEFF, LDDCOE, UCOEFF,
     $                   LDUCO1, LDUCO2, AT, N1, TAU, TOL1, TOL2, IWORK,
     $                   DWORK, LDWORK, INFO )
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
C     Calculates the (PWORK x MWORK) transfer matrix T(s), in the form
C     of polynomial row vectors over monic least common denominator
C     polynomials, of a given state-space representation (ssr).  Each
C     such row of T(s) is simply a single-output relatively left prime
C     polynomial matrix representation (pmr), so can be calculated by
C     applying a simplified version of the Orthogonal Structure
C     Theorem to a minimal ssr for the corresponding row of the given
C     system: such an ssr is obtained by using the Orthogonal Canon-
C     ical Form to first separate out a completely controllable one
C     for the overall system and then, for each row in turn, applying
C     it again to the resulting dual SIMO system.  The Orthogonal
C     Structure Theorem produces non-monic denominator and V:I(s)
C     polynomials: this is avoided here by first scaling AT (the
C     transpose of the controllable part of A, found in this routine)
C     by suitable products of its sub-diagonal elements (these are then
C     no longer needed, so freeing the entire lower triangle for
C     storing the coefficients of V(s) apart from the leading 1's,
C     which are treated implicitly).  These polynomials are calculated
C     in reverse order (IW = NMINL - 1,...,1), the monic denominator
C     D:I(s) found exactly as if it were V:0(s), and finally the
C     numerator vector U:I(s) obtained from the Orthogonal Structure
C     Theorem relation.
C
C     ******************************************************************
C
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDDCOE, LDUCO1,
     $                  LDUCO2, LDWORK, MWORK, N, N1, NCONT, PWORK
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           INDEXD(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AT(N1,*), B(LDB,*), C(LDC,*),
     $                  D(LDD,*), DCOEFF(LDDCOE,*), DWORK(*),
     $                  UCOEFF(LDUCO1,LDUCO2,*), TAU(*)
C     .. Local Scalars ..
      INTEGER           I, IB, IBI, IC, INDCON, IS, IV, IVMIN1, IWPLUS,
     $                  IZ, J, JWORK, K, L, LWORK, MAXM, NMINL, NPLUS,
     $                  WRKOPT
      DOUBLE PRECISION  TEMP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, TB01UD, TB01ZD
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C     .. Executable Statements ..
C
C     Separate out controllable subsystem (of order NCONT).
C
C     Workspace: MAX(N, 3*MWORK, PWORK).
C
      CALL TB01UD( 'No Z', N, MWORK, PWORK, A, LDA, B, LDB, C, LDC,
     $             NCONT, INDCON, IWORK, AT, 1, TAU, TOL2, IWORK(N+1),
     $             DWORK, LDWORK, INFO )
      WRKOPT = INT( DWORK(1) )
C
      IS = 1
      IC = IS + NCONT
      IZ = IC
      IB = IC + NCONT
      LWORK = IB + MWORK*NCONT
      MAXM  = MAX( 1, MWORK )
C
C     Calculate each row of T(s) in turn.
C
      DO 140 I = 1, PWORK
C
C        Form the dual of I-th NCONT-order MISO subsystem ...
C
         CALL DCOPY( NCONT, C(I,1), LDC, DWORK(IC), 1 )
C
         DO 10 J = 1, NCONT
            CALL DCOPY( NCONT, A(J,1), LDA, AT(1,J), 1 )
            CALL DCOPY( MWORK, B(J,1), LDB, DWORK((J-1)*MAXM+IB), 1 )
   10    CONTINUE
C
C        and separate out its controllable part, giving minimal
C        state-space realization for row I.
C
C        Workspace: MWORK*NCONT + 2*NCONT + MAX(NCONT,MWORK).
C
         CALL TB01ZD( 'No Z', NCONT, MWORK, AT, N1, DWORK(IC),
     $                DWORK(IB), MAXM, NMINL, DWORK(IZ), 1, TAU, TOL1,
     $                DWORK(LWORK), LDWORK-LWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(LWORK) )+LWORK-1 )
C
C        Store degree of (monic) denominator, and leading coefficient
C        vector of numerator.
C
         INDEXD(I) = NMINL
         DCOEFF(I,1) = ONE
         CALL DCOPY( MWORK, D(I,1), LDD, UCOEFF(I,1,1), LDUCO1 )
C
         IF ( NMINL.EQ.1 ) THEN
C
C           Finish off numerator, denominator for simple case NMINL=1.
C
            TEMP = -AT(1,1)
            DCOEFF(I,2) = TEMP
            CALL DCOPY( MWORK, D(I,1), LDD, UCOEFF(I,1,2), LDUCO1 )
            CALL DSCAL( MWORK, TEMP, UCOEFF(I,1,2), LDUCO1 )
            CALL DAXPY( MWORK, DWORK(IC), DWORK(IB), 1, UCOEFF(I,1,2),
     $                  LDUCO1 )
         ELSE IF ( NMINL.GT.1 ) THEN
C
C           Set up factors for scaling upper triangle of AT ...
C
            CALL DCOPY( NMINL-1, AT(2,1), N1+1, DWORK(IC+1), 1 )
            NPLUS = NMINL + 1
C
            DO 20 L = IS, IS + NMINL - 1
               DWORK(L) = ONE
   20       CONTINUE
C
C           and scale it, row by row, starting with row NMINL.
C
            DO 40 JWORK = NMINL, 1, -1
C
               DO 30 J = JWORK, NMINL
                  AT(JWORK,J) = DWORK(IS+J-1)*AT(JWORK,J)
   30          CONTINUE
C
C              Update scale factors for next row.
C
               CALL DSCAL( NMINL-JWORK+1, DWORK(IC+JWORK-1),
     $                     DWORK(IS+JWORK-1), 1 )
   40       CONTINUE
C
C           Calculate each monic polynomial V:JWORK(s) in turn:
C           K-th coefficient stored as AT(IV,K-1).
C
            DO 70 IV = 2, NMINL
               JWORK = NPLUS - IV
               IWPLUS = JWORK + 1
               IVMIN1 = IV - 1
C
C              Set up coefficients due to leading 1's of existing
C              V:I(s)'s.
C
               DO 50 K = 1, IVMIN1
                  AT(IV,K) = -AT(IWPLUS,JWORK+K)
   50          CONTINUE
C
               IF ( IV.NE.2 ) THEN
C
C                 Then add contribution from s * V:JWORK+1(s) term.
C
                  CALL DAXPY( IV-2, ONE, AT(IVMIN1,1), N1, AT(IV,1),
     $                        N1 )
C
C                 Finally, add effect of lower coefficients of existing
C                 V:I(s)'s.
C
                  DO 60 K = 2, IVMIN1
                     AT(IV,K) = AT(IV,K) - DDOT( K-1,
     $                                           AT(IWPLUS,JWORK+1), N1,
     $                                           AT(IV-K+1,1), -(N1+1) )
   60             CONTINUE
C
               END IF
   70       CONTINUE
C
C           Determine denominator polynomial D(s) as if it were V:0(s).
C
            DO 80 K = 2, NPLUS
               DCOEFF(I,K) = -AT(1,K-1)
   80       CONTINUE
C
            CALL DAXPY( NMINL-1, ONE, AT(NMINL,1), N1, DCOEFF(I,2),
     $                  LDDCOE )
C
            DO 90 K = 3, NPLUS
               DCOEFF(I,K) = DCOEFF(I,K) - DDOT( K-2, AT, N1,
     $                       AT(NMINL-K+3,1), -(N1+1) )
   90       CONTINUE
C
C           Scale (B' * Z), stored in DWORK(IB).
C
            IBI = IB
C
            DO 100 L = 1, NMINL
               CALL DSCAL( MWORK, DWORK(IS+L-1), DWORK(IBI), 1 )
               IBI = IBI + MAXM
  100       CONTINUE
C
C           Evaluate numerator polynomial vector (V(s) * B) + (D(s)
C           * D:I): first set up coefficients due to D:I and leading
C           1's of V(s).
C
            IBI = IB
C
            DO 110 K = 2, NPLUS
               CALL DCOPY( MWORK, DWORK(IBI), 1, UCOEFF(I,1,K), LDUCO1 )
               CALL DAXPY( MWORK, DCOEFF(I,K), D(I,1), LDD,
     $                     UCOEFF(I,1,K), LDUCO1 )
               IBI = IBI + MAXM
  110       CONTINUE
C
C           Add contribution from lower coefficients of V(s).
C
            DO 130 K = 3, NPLUS
C
               DO 120 J = 1, MWORK
                  UCOEFF(I,J,K) = UCOEFF(I,J,K) + DDOT( K-2,
     $                                  AT(NMINL-K+3,1), -(N1+1),
     $                                  DWORK(IB+J-1), MAXM )
  120          CONTINUE
C
  130       CONTINUE
C
         END IF
  140 CONTINUE
C
C     Set optimal workspace dimension.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TB04AY ***
      END
