      SUBROUTINE TC01OD( LERI, M, P, INDLIM, PCOEFF, LDPCO1, LDPCO2,
     $                   QCOEFF, LDQCO1, LDQCO2, INFO )
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
C     To find the dual right (left) polynomial matrix representation of
C     a given left (right) polynomial matrix representation, where the
C     right and left polynomial matrix representations are of the form
C     Q(s)*inv(P(s)) and inv(P(s))*Q(s) respectively.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LERI    CHARACTER*1
C             Indicates whether a left or right matrix fraction is input
C             as follows:
C             = 'L':  A left matrix fraction is input;
C             = 'R':  A right matrix fraction is input.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     INDLIM  (input) INTEGER
C             The highest value of K for which PCOEFF(.,.,K) and
C             QCOEFF(.,.,K) are to be transposed.
C             K = kpcoef + 1, where kpcoef is the maximum degree of the
C             polynomials in P(s).  INDLIM >= 1.
C
C     PCOEFF  (input/output) DOUBLE PRECISION array, dimension
C             (LDPCO1,LDPCO2,INDLIM)
C             If LERI = 'L' then porm = P, otherwise porm = M.
C             On entry, the leading porm-by-porm-by-INDLIM part of this
C             array must contain the coefficients of the denominator
C             matrix P(s).
C             PCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of
C             polynomial (I,J) of P(s), where K = 1,2,...,INDLIM.
C             On exit, the leading porm-by-porm-by-INDLIM part of this
C             array contains the coefficients of the denominator matrix
C             P'(s) of the dual system.
C
C     LDPCO1  INTEGER
C             The leading dimension of array PCOEFF.
C             LDPCO1 >= MAX(1,P) if LERI = 'L',
C             LDPCO1 >= MAX(1,M) if LERI = 'R'.
C
C     LDPCO2  INTEGER
C             The second dimension of array PCOEFF.
C             LDPCO2 >= MAX(1,P) if LERI = 'L',
C             LDPCO2 >= MAX(1,M) if LERI = 'R'.
C
C     QCOEFF  (input/output) DOUBLE PRECISION array, dimension
C             (LDQCO1,LDQCO2,INDLIM)
C             On entry, the leading P-by-M-by-INDLIM part of this array
C             must contain the coefficients of the numerator matrix
C             Q(s).
C             QCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of
C             polynomial (I,J) of Q(s), where K = 1,2,...,INDLIM.
C             On exit, the leading M-by-P-by-INDLIM part of the array
C             contains the coefficients of the numerator matrix Q'(s)
C             of the dual system.
C
C     LDQCO1  INTEGER
C             The leading dimension of array QCOEFF.
C             LDQCO1 >= MAX(1,M,P).
C
C     LDQCO2  INTEGER
C             The second dimension of array QCOEFF.
C             LDQCO2 >= MAX(1,M,P).
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
C     If the given M-input/P-output left (right) polynomial matrix
C     representation has numerator matrix Q(s) and denominator matrix
C     P(s), its dual P-input/M-output right (left) polynomial matrix
C     representation simply has numerator matrix Q'(s) and denominator
C     matrix P'(s).
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TC01CD by T.W.C.Williams, Kingston
C     Polytechnic, United Kingdom, March 1982.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Coprime matrix fraction, elementary polynomial operations,
C     polynomial matrix, state-space representation, transfer matrix.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER         LERI
      INTEGER           INFO, INDLIM, LDPCO1, LDPCO2, LDQCO1, LDQCO2, M,
     $                  P
C     .. Array Arguments ..
      DOUBLE PRECISION  PCOEFF(LDPCO1,LDPCO2,*), QCOEFF(LDQCO1,LDQCO2,*)
C     .. Local Scalars ..
      LOGICAL           LLERI
      INTEGER           J, K, MINMP, MPLIM, PORM
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LLERI = LSAME( LERI, 'L' )
      MPLIM = MAX( M, P )
      MINMP = MIN( M, P )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LLERI .AND. .NOT.LSAME( LERI, 'R' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( INDLIM.LT.1 ) THEN
         INFO = -4
      ELSE IF( ( LLERI .AND. LDPCO1.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO1.LT.MAX( 1, M ) ) ) THEN
         INFO = -6
      ELSE IF( ( LLERI .AND. LDPCO2.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO2.LT.MAX( 1, M ) ) ) THEN
         INFO = -7
      ELSE IF( LDQCO1.LT.MAX( 1, MPLIM ) ) THEN
         INFO = -9
      ELSE IF( LDQCO2.LT.MAX( 1, MPLIM ) ) THEN
         INFO = -10
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TC01OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( M.EQ.0 .OR. P.EQ.0 )
     $   RETURN
C
      IF ( MPLIM.NE.1 ) THEN
C
C        Non-scalar system: transpose numerator matrix Q(s).
C
         DO 20 K = 1, INDLIM
C
            DO 10 J = 1, MPLIM
               IF ( J.LT.MINMP ) THEN
                  CALL DSWAP( MINMP-J, QCOEFF(J+1,J,K), 1,
     $                        QCOEFF(J,J+1,K), LDQCO1 )
               ELSE IF ( J.GT.P ) THEN
                  CALL DCOPY( P, QCOEFF(1,J,K), 1, QCOEFF(J,1,K),
     $                        LDQCO1 )
               ELSE IF ( J.GT.M ) THEN
                  CALL DCOPY( M, QCOEFF(J,1,K), LDQCO1, QCOEFF(1,J,K),
     $                        1 )
               END IF
   10       CONTINUE
C
   20    CONTINUE
C
C        Find dimension of denominator matrix P(s): M (P) for
C        right (left) polynomial matrix representation.
C
         PORM = M
         IF ( LLERI ) PORM = P
         IF ( PORM.NE.1 ) THEN
C
C           Non-scalar P(s): transpose it.
C
            DO 40 K = 1, INDLIM
C
               DO 30 J = 1, PORM - 1
                  CALL DSWAP( PORM-J, PCOEFF(J+1,J,K), 1,
     $                        PCOEFF(J,J+1,K), LDPCO1 )
   30          CONTINUE
C
   40       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of TC01OD ***
      END
