      SUBROUTINE MC03MD( RP1, CP1, CP2, DP1, DP2, DP3, ALPHA, P1,
     $                   LDP11, LDP12, P2, LDP21, LDP22, P3, LDP31,
     $                   LDP32, DWORK, INFO )
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
C     To compute the coefficients of the real polynomial matrix
C
C        P(x) = P1(x) * P2(x) + alpha * P3(x),
C
C     where P1(x), P2(x) and P3(x) are given real polynomial matrices
C     and alpha is a real scalar.
C
C     Each of the polynomial matrices P1(x), P2(x) and P3(x) may be the
C     zero matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     RP1     (input) INTEGER
C             The number of rows of the matrices P1(x) and P3(x).
C             RP1 >= 0.
C
C     CP1     (input) INTEGER
C             The number of columns of matrix P1(x) and the number of
C             rows of matrix P2(x).  CP1 >= 0.
C
C     CP2     (input) INTEGER
C             The number of columns of the matrices P2(x) and P3(x).
C             CP2 >= 0.
C
C     DP1     (input) INTEGER
C             The degree of the polynomial matrix P1(x).  DP1 >= -1.
C
C     DP2     (input) INTEGER
C             The degree of the polynomial matrix P2(x).  DP2 >= -1.
C
C     DP3     (input/output) INTEGER
C             On entry, the degree of the polynomial matrix P3(x).
C             DP3 >= -1.
C             On exit, the degree of the polynomial matrix P(x).
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar value alpha of the problem.
C
C     P1      (input) DOUBLE PRECISION array, dimension (LDP11,LDP12,*)
C             If DP1 >= 0, then the leading RP1-by-CP1-by-(DP1+1) part
C             of this array must contain the coefficients of the
C             polynomial matrix P1(x). Specifically, P1(i,j,k) must
C             contain the coefficient of x**(k-1) of the polynomial
C             which is the (i,j)-th element of P1(x), where i = 1,2,...,
C             RP1, j = 1,2,...,CP1 and k = 1,2,...,DP1+1.
C             If DP1 = -1, then P1(x) is taken to be the zero polynomial
C             matrix, P1 is not referenced and can be supplied as a
C             dummy array (i.e. set the parameters LDP11 = LDP12 = 1 and
C             declare this array to be P1(1,1,1) in the calling
C             program).
C
C     LDP11   INTEGER
C             The leading dimension of array P1.
C             LDP11 >= MAX(1,RP1) if DP1 >= 0,
C             LDP11 >= 1          if DP1 = -1.
C
C     LDP12   INTEGER
C             The second dimension of array P1.
C             LDP12 >= MAX(1,CP1) if DP1 >= 0,
C             LDP12 >= 1          if DP1 = -1.
C
C     P2      (input) DOUBLE PRECISION array, dimension (LDP21,LDP22,*)
C             If DP2 >= 0, then the leading CP1-by-CP2-by-(DP2+1) part
C             of this array must contain the coefficients of the
C             polynomial matrix P2(x). Specifically, P2(i,j,k) must
C             contain the coefficient of x**(k-1) of the polynomial
C             which is the (i,j)-th element of P2(x), where i = 1,2,...,
C             CP1, j = 1,2,...,CP2 and k = 1,2,...,DP2+1.
C             If DP2 = -1, then P2(x) is taken to be the zero polynomial
C             matrix, P2 is not referenced and can be supplied as a
C             dummy array (i.e. set the parameters LDP21 = LDP22 = 1 and
C             declare this array to be P2(1,1,1) in the calling
C             program).
C
C     LDP21   INTEGER
C             The leading dimension of array P2.
C             LDP21 >= MAX(1,CP1) if DP2 >= 0,
C             LDP21 >= 1          if DP2 = -1.
C
C     LDP22   INTEGER
C             The second dimension of array P2.
C             LDP22 >= MAX(1,CP2) if DP2 >= 0,
C             LDP22 >= 1          if DP2 = -1.
C
C     P3      (input/output) DOUBLE PRECISION array, dimension
C             (LDP31,LDP32,n), where n = MAX(DP1+DP2,DP3,0)+1.
C             On entry, if DP3 >= 0, then the leading
C             RP1-by-CP2-by-(DP3+1) part of this array must contain the
C             coefficients of the polynomial matrix P3(x). Specifically,
C             P3(i,j,k) must contain the coefficient of x**(k-1) of the
C             polynomial which is the (i,j)-th element of P3(x), where
C             i = 1,2,...,RP1, j = 1,2,...,CP2 and k = 1,2,...,DP3+1.
C             If DP3 = -1, then P3(x) is taken to be the zero polynomial
C             matrix.
C             On exit, if DP3 >= 0 on exit (ALPHA <> 0.0 and DP3 <> -1,
C             on entry, or DP1 <> -1 and DP2 <> -1), then the leading
C             RP1-by-CP2-by-(DP3+1) part of this array contains the
C             coefficients of P(x). Specifically, P3(i,j,k) contains the
C             coefficient of x**(k-1) of the polynomial which is the
C             (i,j)-th element of P(x), where i = 1,2,...,RP1, j = 1,2,
C             ...,CP2 and k = 1,2,...,DP3+1.
C             If DP3 = -1 on exit, then the coefficients of P(x) (the
C             zero polynomial matrix) are not stored in the array.
C
C     LDP31   INTEGER
C             The leading dimension of array P3.  LDP31 >= MAX(1,RP1).
C
C     LDP32   INTEGER
C             The second dimension of array P3.   LDP32 >= MAX(1,CP2).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (CP1)
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
C     Given real polynomial matrices
C
C                DP1            i
C        P1(x) = SUM (A(i+1) * x ),
C                i=0
C
C                DP2            i
C        P2(x) = SUM (B(i+1) * x ),
C                i=0
C
C                DP3            i
C        P3(x) = SUM (C(i+1) * x )
C                i=0
C
C     and a real scalar alpha, the routine computes the coefficients
C     d ,d ,..., of the polynomial matrix
C      1  2
C
C        P(x) = P1(x) * P2(x) + alpha * P3(x)
C
C     from the formula
C
C                 s
C        d    =  SUM (A(k+1) * B(i-k+1)) + alpha * C(i+1),
C         i+1    k=r
C
C     where i = 0,1,...,DP1+DP2 and r and s depend on the value of i
C     (e.g. if i <= DP1 and i <= DP2, then r = 0 and s = i).
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     FURTHER COMMENTS
C
C     Other elementary operations involving polynomial matrices can
C     easily be obtained by calling the appropriate BLAS routine(s).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC03AD by A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, input output description,
C     polynomial matrix, polynomial operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           CP1, CP2, DP1, DP2, DP3, INFO, LDP11, LDP12,
     $                  LDP21, LDP22, LDP31, LDP32, RP1
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(*), P1(LDP11,LDP12,*), P2(LDP21,LDP22,*),
     $                  P3(LDP31,LDP32,*)
C     .. Local Scalars ..
      LOGICAL           CFZERO
      INTEGER           DPOL3, E, H, I, J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLASET, DSCAL, XERBLA
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( RP1.LT.0 ) THEN
         INFO = -1
      ELSE IF( CP1.LT.0 ) THEN
         INFO = -2
      ELSE IF( CP2.LT.0 ) THEN
         INFO = -3
      ELSE IF( DP1.LT.-1 ) THEN
         INFO = -4
      ELSE IF( DP2.LT.-1 ) THEN
         INFO = -5
      ELSE IF( DP3.LT.-1 ) THEN
         INFO = -6
      ELSE IF( ( DP1.EQ.-1 .AND. LDP11.LT.1               ) .OR.
     $         ( DP1.GE. 0 .AND. LDP11.LT.MAX( 1, RP1 ) ) ) THEN
         INFO = -9
      ELSE IF( ( DP1.EQ.-1 .AND. LDP12.LT.1               ) .OR.
     $         ( DP1.GE. 0 .AND. LDP12.LT.MAX( 1, CP1 ) ) ) THEN
         INFO = -10
      ELSE IF( ( DP2.EQ.-1 .AND. LDP21.LT.1               ) .OR.
     $         ( DP2.GE. 0 .AND. LDP21.LT.MAX( 1, CP1 ) ) ) THEN
         INFO = -12
      ELSE IF( ( DP2.EQ.-1 .AND. LDP22.LT.1               ) .OR.
     $         ( DP2.GE. 0 .AND. LDP22.LT.MAX( 1, CP2 ) ) ) THEN
         INFO = -13
      ELSE IF( LDP31.LT.MAX( 1, RP1 ) ) THEN
         INFO = -15
      ELSE IF( LDP32.LT.MAX( 1, CP2 ) ) THEN
         INFO = -16
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC03MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( RP1.EQ.0 .OR. CP2.EQ.0 )
     $   RETURN
C
      IF ( ALPHA.EQ.ZERO )
     $   DP3 = -1
C
      IF ( DP3.GE.0 ) THEN
C
C        P3(x) := ALPHA * P3(x).
C
         DO 40 K = 1, DP3 + 1
C
            DO 20 J = 1, CP2
               CALL DSCAL( RP1, ALPHA, P3(1,J,K), 1 )
   20       CONTINUE
C
   40    CONTINUE
      END IF
C
      IF ( ( DP1.EQ.-1 ) .OR. ( DP2.EQ.-1 ) .OR. ( CP1.EQ.0 ) )
     $   RETURN
C
C     Neither of P1(x) and P2(x) is the zero polynomial.
C
      DPOL3 = DP1 + DP2
      IF ( DPOL3.GT.DP3 ) THEN
C
C        Initialize the additional part of P3(x) to zero.
C
         DO 80 K = DP3 + 2, DPOL3 + 1
            CALL DLASET( 'Full', RP1, CP2, ZERO, ZERO, P3(1,1,K),
     $                   LDP31 )
   80    CONTINUE
C
         DP3 = DPOL3
      END IF
C                                                              k-1
C     The inner product of the j-th row of the coefficient of x    of P1
C                                                i-1
C     and the h-th column of the coefficient of x    of P2(x) contribute
C                                                 k+i-2
C     the (j,h)-th element of the coefficient of x      of P3(x).
C
      DO 160 K = 1, DP1 + 1
C
         DO 140 J = 1, RP1
            CALL DCOPY( CP1, P1(J,1,K), LDP11, DWORK, 1 )
C
            DO 120 I = 1, DP2 + 1
               E = K + I - 1
C
               DO 100 H = 1, CP2
                  P3(J,H,E) = DDOT( CP1, DWORK, 1, P2(1,H,I), 1 ) +
     $                        P3(J,H,E)
  100          CONTINUE
C
  120       CONTINUE
C
  140    CONTINUE
C
  160 CONTINUE
C
C     Computation of the exact degree of P3(x).
C
      CFZERO = .TRUE.
C     WHILE ( DP3 >= 0 and CFZERO ) DO
  180 IF ( ( DP3.GE.0 ) .AND. CFZERO ) THEN
         DPOL3 = DP3 + 1
C
         DO 220 J = 1, CP2
C
            DO 200 I = 1, RP1
               IF ( P3(I,J,DPOL3 ).NE.ZERO ) CFZERO = .FALSE.
  200       CONTINUE
C
  220    CONTINUE
C
         IF ( CFZERO ) DP3 = DP3 - 1
         GO TO 180
      END IF
C     END WHILE 180
C
      RETURN
C *** Last line of MC03MD ***
      END
