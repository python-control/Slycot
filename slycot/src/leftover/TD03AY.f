      SUBROUTINE TD03AY( MWORK, PWORK, INDEX, DCOEFF, LDDCOE, UCOEFF,
     $                   LDUCO1, LDUCO2, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, INFO )
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
C     Calculates a state-space representation for a (PWORK x MWORK)
C     transfer matrix given in the form of polynomial row vectors over
C     common denominators (not necessarily lcd's).  Such a description
C     is simply the polynomial matrix representation
C
C          T(s) = inv(D(s)) * U(s),
C
C     where D(s) is diagonal with (I,I)-th element D:I(s) of degree
C     INDEX(I); applying Wolovich's Observable Structure Theorem to
C     this left matrix fraction then yields an equivalent state-space
C     representation in observable companion form, of order
C     N = sum(INDEX(I)).  As D(s) is diagonal, the PWORK ordered
C     'non-trivial' columns of C and A are very simply calculated, these
C     submatrices being diagonal and (INDEX(I) x 1) - block diagonal,
C     respectively: finding B and D is also somewhat simpler than for
C     general P(s) as dealt with in TC04AD. Finally, the state-space
C     representation obtained here is not necessarily controllable
C     (as D(s) and U(s) are not necessarily relatively left prime), but
C     it is theoretically completely observable: however, its
C     observability matrix may be poorly conditioned, so it is safer
C     not to assume observability either.
C
C     REVISIONS
C
C     May 13, 1998.
C
C     KEYWORDS
C
C     Coprime matrix fraction, elementary polynomial operations,
C     polynomial matrix, state-space representation, transfer matrix.
C
C     ******************************************************************
C
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDDCOE, LDUCO1,
     $                  LDUCO2, MWORK, N, PWORK
C     .. Array Arguments ..
      INTEGER           INDEX(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DCOEFF(LDDCOE,*), UCOEFF(LDUCO1,LDUCO2,*)
C     .. Local Scalars ..
      INTEGER           I, IA, IBIAS, INDCUR, JA, JMAX1, K
      DOUBLE PRECISION  ABSDIA, ABSDMX, BIGNUM, DIAG, SMLNUM, UMAX1,
     $                  TEMP
C     .. External Functions ..
      INTEGER           IDAMAX
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DLASET, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      INFO = 0
C
C     Initialize A and C to be zero, apart from 1's on the subdiagonal
C     of A.
C
      CALL DLASET( 'Upper', N, N, ZERO, ZERO, A, LDA )
      IF ( N.GT.1 ) CALL DLASET( 'Lower', N-1, N-1, ZERO, ONE, A(2,1),
     $                           LDA )
C
      CALL DLASET( 'Full', PWORK, N, ZERO, ZERO, C, LDC )
C
C     Calculate B and D, as well as 'non-trivial' elements of A and C.
C     Check if any leading coefficient of D(s) nearly zero: if so, exit.
C     Caution is taken to avoid overflow.
C
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
C
      IBIAS = 2
      JA = 0
C
      DO 20 I = 1, PWORK
         ABSDIA = ABS( DCOEFF(I,1) )
         JMAX1  = IDAMAX( MWORK, UCOEFF(I,1,1), LDUCO1 )
         UMAX1  = ABS( UCOEFF(I,JMAX1,1) )
         IF ( ( ABSDIA.LT.SMLNUM ) .OR.
     $        ( ABSDIA.LT.ONE .AND. UMAX1.GT.ABSDIA*BIGNUM ) ) THEN
C
C           Error return.
C
            INFO = I
            RETURN
         END IF
         DIAG   = ONE/DCOEFF(I,1)
         INDCUR = INDEX(I)
         IF ( INDCUR.NE.0 ) THEN
            IBIAS = IBIAS + INDCUR
            JA = JA + INDCUR
            IF ( INDCUR.GE.1 ) THEN
               JMAX1  = IDAMAX( INDCUR, DCOEFF(I,2), LDDCOE )
               ABSDMX = ABS( DCOEFF(I,JMAX1) )
               IF ( ABSDIA.GE.ONE ) THEN
                  IF ( UMAX1.GT.ONE ) THEN
                     IF ( ( ABSDMX/ABSDIA ).GT.( BIGNUM/UMAX1 ) ) THEN
C
C                       Error return.
C
                        INFO = I
                        RETURN
                     END IF
                  END IF
               ELSE
                  IF ( UMAX1.GT.ONE ) THEN
                     IF ( ABSDMX.GT.( BIGNUM*ABSDIA )/UMAX1 ) THEN
C
C                       Error return.
C
                        INFO = I
                        RETURN
                     END IF
                  END IF
               END IF
            END IF
C
C           I-th 'non-trivial' sub-vector of A given from coefficients
C           of D:I(s), while I-th row block of B given from this and
C           row I of U(s).
C
            DO 10 K = 2, INDCUR + 1
               IA = IBIAS - K
               TEMP = -DIAG*DCOEFF(I,K)
               A(IA,JA) = TEMP
C
               CALL DCOPY( MWORK, UCOEFF(I,1,K), LDUCO1, B(IA,1), LDB )
               CALL DAXPY( MWORK, TEMP, UCOEFF(I,1,1), LDUCO1, B(IA,1),
     $                     LDB )
   10       CONTINUE
C
            IF ( JA.LT.N ) A(JA+1,JA) = ZERO
C
C           Finally, I-th 'non-trivial' entry of C and row of D obtained
C           also.
C
            C(I,JA) = DIAG
         END IF
C
         CALL DCOPY( MWORK, UCOEFF(I,1,1), LDUCO1, D(I,1), LDD )
         CALL DSCAL( MWORK, DIAG, D(I,1), LDD )
   20 CONTINUE
C
      RETURN
C *** Last line of TD03AY ***
      END
