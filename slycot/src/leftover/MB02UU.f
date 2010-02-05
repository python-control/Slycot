      SUBROUTINE MB02UU( N, A, LDA, RHS, IPIV, JPIV, SCALE )
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
C     To solve for x in A * x = scale * RHS, using the LU factorization
C     of the N-by-N matrix A computed by SLICOT Library routine MB02UV.
C     The factorization has the form A = P * L * U * Q, where P and Q
C     are permutation matrices, L is unit lower triangular and U is
C     upper triangular.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA, N)
C             The leading N-by-N part of this array must contain
C             the LU part of the factorization of the matrix A computed
C             by SLICOT Library routine MB02UV:  A = P * L * U * Q.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1, N).
C
C     RHS     (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the right hand side
C             of the system.
C             On exit, this array contains the solution of the system.
C
C     IPIV    (input) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= i <= N, row i of the
C             matrix has been interchanged with row IPIV(i).
C
C     JPIV    (input) INTEGER array, dimension (N)
C             The pivot indices; for 1 <= j <= N, column j of the
C             matrix has been interchanged with column JPIV(j).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scale factor, chosen 0 < SCALE <= 1 to prevent
C             overflow in the solution.
C
C     FURTHER COMMENTS
C
C     In the interest of speed, this routine does not check the input
C     for errors. It should only be used if the order of the matrix A
C     is very small.
C
C     CONTRIBUTOR
C
C     Bo Kagstrom and P. Poromaa, Univ. of Umea, Sweden, Nov. 1993.
C
C     REVISIONS
C
C     April 1998 (T. Penzl).
C     Sep. 1998 (V. Sima).
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
C     .. Scalar Arguments ..
      INTEGER            LDA, N
      DOUBLE PRECISION   SCALE
C     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), RHS( * )
C     .. Local Scalars ..
      INTEGER            I, IP, J
      DOUBLE PRECISION   BIGNUM, EPS, FACTOR, SMLNUM, TEMP
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, IDAMAX
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DLABAD, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX
C     .. Executable Statements ..
C
C     Set constants to control owerflow.
C
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = DLAMCH( 'Safe minimum' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
C
C     Apply permutations IPIV to RHS.
C
      DO 20 I = 1, N - 1
         IP = IPIV(I)
         IF ( IP.NE.I ) THEN
            TEMP    = RHS(I)
            RHS(I)  = RHS(IP)
            RHS(IP) = TEMP
         ENDIF
   20 CONTINUE
C
C     Solve for L part.
C
      DO 40 I = 1, N - 1
         CALL DAXPY( N-I, -RHS(I), A(I+1, I), 1, RHS(I+1), 1 )
  40  CONTINUE
C
C     Solve for U part.
C
C     Check for scaling.
C
      FACTOR = TWO * DBLE( N )
      I = 1
   60 CONTINUE
      IF ( ( FACTOR * SMLNUM ) * ABS( RHS(I) ) .LE. ABS( A(I, I) ) )
     $       THEN
         I = I + 1
         IF ( I .LE. N ) GO TO 60
         SCALE = ONE
      ELSE
         SCALE = ( ONE / FACTOR ) / ABS( RHS( IDAMAX( N, RHS, 1 ) ) )
         CALL DSCAL( N, SCALE, RHS, 1 )
      END IF
C
      DO 100 I = N, 1, -1
         TEMP = ONE / A(I, I)
         RHS(I) = RHS(I) * TEMP
         DO 80 J = I + 1, N
            RHS(I) = RHS(I) - RHS(J) * ( A(I, J) * TEMP )
   80    CONTINUE
  100 CONTINUE
C
C     Apply permutations JPIV to the solution (RHS).
C
      DO 120 I = N - 1, 1, -1
         IP = JPIV(I)
         IF ( IP.NE.I ) THEN
            TEMP    = RHS(I)
            RHS(I)  = RHS(IP)
            RHS(IP) = TEMP
         ENDIF
  120 CONTINUE
C
      RETURN
C *** Last line of MB02UU ***
      END
