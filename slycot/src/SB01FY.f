      SUBROUTINE SB01FY( DISCR, N, M, A, LDA, B, LDB, F, LDF, V, LDV,
     $                   INFO )
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
C     To compute the inner denominator of a right-coprime factorization
C     of a system of order N, where N is either 1 or 2. Specifically,
C     given the N-by-N unstable system state matrix A and the N-by-M
C     system input matrix B, an M-by-N state-feedback matrix F and
C     an M-by-M matrix V are constructed, such that the system
C     (A + B*F, B*V, F, V) is inner.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DISCR   LOGICAL
C             Specifies the type of system as follows:
C             = .FALSE.:  continuous-time system;
C             = .TRUE. :  discrete-time system.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A and also the number of rows of
C             the matrix B and the number of columns of the matrix F.
C             N is either 1 or 2.
C
C     M       (input) INTEGER
C             The number of columns of the matrices B and V, and also
C             the number of rows of the matrix F.  M >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A whose eigenvalues must have positive
C             real parts if DISCR = .FALSE. or moduli greater than unity
C             if DISCR = .TRUE..
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= N.
C
C     F       (output) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading M-by-N part of this array contains the state-
C             feedback matrix F which assigns one eigenvalue (if N = 1)
C             or two eigenvalues (if N = 2) of the matrix A + B*F in
C             symmetric positions with respect to the imaginary axis
C             (if DISCR = .FALSE.) or the unit circle (if
C             DISCR = .TRUE.).
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,M)
C             The leading M-by-M upper triangular part of this array
C             contains the input/output matrix V of the resulting inner
C             system in upper triangular form.
C             If DISCR = .FALSE., the resulting V is an identity matrix.
C
C     LDV     INTEGER
C             The leading dimension of array V.  LDF >= MAX(1,M).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if uncontrollability of the pair (A,B) is detected;
C             = 2:  if A is stable or at the stability limit;
C             = 3:  if N = 2 and A has a pair of real eigenvalues.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine RCFID2.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     Feb. 1999, A. Varga, DLR Oberpfaffenhofen.
C
C     KEYWORDS
C
C     Coprime factorization, eigenvalue, eigenvalue assignment,
C     feedback control, pole placement, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           DISCR
      INTEGER           INFO, LDA, LDB, LDF, LDV, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), F(LDF,*), V(LDV,*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  CS, R11, R12, R22, SCALE, SN, TEMP
C     .. Local Arrays ..
      DOUBLE PRECISION  AT(2,2), DUMMY(2,2), U(2,2)
C     .. External Functions ..
      DOUBLE PRECISION  DLAPY2, DLAPY3
      EXTERNAL          DLAPY2, DLAPY3
C     .. External Subroutines ..
      EXTERNAL          DLARFG, DLASET, SLCT_DLATZM, DROTG, DTRTRI,
     $                  MA02AD, MB04OX, SB03OY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     For efficiency reasons, the parameters are not checked.
C
      INFO = 0
C
C     Compute an N-by-N upper triangular R such that R'*R = B*B' and
C     find an upper triangular matrix U in the equation
C
C     A'*U'*U + U'*U*A = R'*R if DISCR = .FALSE. or
C     A'*U'*U*A - U'*U = R'*R if DISCR = .TRUE. .
C
      CALL MA02AD( 'Full', N, M, B, LDB, F, LDF )
C
      IF( N.EQ.1 ) THEN
C
C        The N = 1 case.
C
         IF( M.GT.1 )
     $      CALL DLARFG( M, F(1,1), F(2,1), 1, TEMP )
         R11 = ABS( F(1,1) )
C
C        Make sure A is unstable or divergent and find U.
C
         IF( DISCR ) THEN
            TEMP = ABS( A(1,1) )
            IF( TEMP.LE.ONE ) THEN
               INFO = 2
               RETURN
            ELSE
               TEMP = R11 / SQRT( ( TEMP - ONE )*( TEMP + ONE ) )
            END IF
         ELSE
            IF( A(1,1).LE.ZERO ) THEN
               INFO = 2
               RETURN
            ELSE
               TEMP = R11 / SQRT( ABS( TWO*A(1,1) ) )
            END IF
         END IF
         U(1,1) = TEMP
         SCALE  = ONE
      ELSE
C
C        The N = 2 case.
C
         IF( M.GT.1 ) THEN
            CALL DLARFG( M, F(1,1), F(2,1), 1, TEMP )
            CALL SLCT_DLATZM( 'Left', M, N-1, F(2,1), 1, TEMP, F(1,2),
     $                   F(2,2), LDF, V )
         END IF
         R11 = F(1,1)
         R12 = F(1,2)
         IF( M.GT.2 )
     $      CALL DLARFG( M-1, F(2,2), F(3,2), 1, TEMP )
         IF( M.EQ.1 ) THEN
            R22 = ZERO
         ELSE
            R22 = F(2,2)
         END IF
         AT(1,1) = A(1,1)
         AT(1,2) = A(2,1)
         AT(2,1) = A(1,2)
         AT(2,2) = A(2,2)
         U(1,1)  = R11
         U(1,2)  = R12
         U(2,2)  = R22
         CALL SB03OY( DISCR, .FALSE., -1, AT, 2, U, 2, DUMMY, 2,
     $                SCALE, INFO )
         IF( INFO.NE.0 ) THEN
            IF( INFO.NE.4 ) THEN
               INFO = 2
            ELSE
               INFO = 3
            END IF
            RETURN
         END IF
      END IF
C
C     Check the controllability of the pair (A,B).
C
C     Warning. Only an exact controllability check is performed.
C              If the pair (A,B) is nearly uncontrollable, then
C              the computed results may be inaccurate.
C
      DO 10 I = 1, N
         IF( U(I,I).EQ.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
   10 CONTINUE
C
C     Set V = I.
C
      CALL DLASET( 'Upper', M, M, ZERO, ONE, V, LDV )
C
      IF( DISCR ) THEN
C
C        Compute an upper triangular matrix V such that
C                                 -1
C        V*V' = (I+B'*inv(U'*U)*B)  .
C
C        First compute F = B'*inv(U) and the Cholesky factorization
C        of I + F*F'.
C
         DO 20 I = 1, M
            F(I,1) = B(1,I)/U(1,1)*SCALE
   20    CONTINUE
         IF( N.EQ.2 ) THEN
            DO 30 I = 1, M
               F(I,2) = ( B(2,I) - F(I,1)*U(1,2) )/U(2,2)*SCALE
   30       CONTINUE
            CALL MB04OX( M, V, LDV, F(1,2), 1 )
         END IF
         CALL MB04OX( M, V, LDV, F(1,1), 1 )
         CALL DTRTRI( 'Upper', 'NonUnit', M, V, LDV, INFO )
      END IF
C
C     Compute the feedback matrix F as:
C
C     1)   If DISCR = .FALSE.
C
C             F = -B'*inv(U'*U);
C
C     2)   If DISCR = .TRUE.
C                                -1
C             F = -B'*(U'*U+B*B')  *A.
C
      IF( N.EQ.1 ) THEN
         IF( DISCR ) THEN
            TEMP = -A(1,1)
            R11  = DLAPY2( U(1,1), R11 )
            DO 40 I = 1, M
               F(I,1) = ( ( B(1,I)/R11 )/R11 )*TEMP
   40       CONTINUE
         ELSE
            R11 = U(1,1)
            DO 50 I = 1, M
               F(I,1) = -( ( B(1,I)/R11 )/R11 )
   50       CONTINUE
         END IF
      ELSE
C
C        Set R = U  if DISCR = .FALSE. or compute the Cholesky
C        factorization of R'*R = U'*U+B*B' if DISCR = .TRUE..
C
         IF( DISCR ) THEN
            TEMP = U(1,1)
            CALL DROTG( R11, TEMP, CS, SN )
            TEMP = -SN*R12 + CS*U(1,2)
            R12  =  CS*R12 + SN*U(1,2)
            R22  = DLAPY3( R22, TEMP, U(2,2) )
         ELSE
            R11 = U(1,1)
            R12 = U(1,2)
            R22 = U(2,2)
         END IF
C
C        Compute F = -B'*inv(R'*R).
C
         DO 60 I = 1, M
            F(I,1) = -B(1,I)/R11
            F(I,2) = -( B(2,I) + F(I,1)*R12 )/R22
            F(I,2) =  F(I,2)/R22
            F(I,1) =  ( F(I,1) - F(I,2)*R12 )/R11
   60    CONTINUE
         IF( DISCR ) THEN
C
C           Compute F <-- F*A.
C
            DO 70 I = 1, M
               TEMP   = F(I,1)*A(1,1) + F(I,2)*A(2,1)
               F(I,2) = F(I,1)*A(1,2) + F(I,2)*A(2,2)
               F(I,1) = TEMP
   70       CONTINUE
         END IF
      END IF
C
      RETURN
C *** Last line of SB01FY ***
      END
