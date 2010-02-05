      SUBROUTINE MB04PY( SIDE, M, N, V, TAU, C, LDC, DWORK )
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
C     To apply a real elementary reflector H to a real m-by-n matrix
C     C, from either the left or the right. H is represented in the form
C                                        ( 1 )
C           H = I - tau * u *u',    u  = (   ),
C                                        ( v )
C     where tau is a real scalar and v is a real vector.
C
C     If tau = 0, then H is taken to be the unit matrix.
C
C     In-line code is used if H has order < 11.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Indicates whether the elementary reflector should be
C             applied from the left or from the right, as follows:
C             = 'L':  Compute H * C;
C             = 'R':  Compute C * H.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix C.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix C.  N >= 0.
C
C     V       (input) DOUBLE PRECISION array, dimension
C             (M-1), if SIDE = 'L', or
C             (N-1), if SIDE = 'R'.
C             The vector v in the representation of H.
C
C     TAU     (input) DOUBLE PRECISION
C             The scalar factor of the elementary reflector H.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix C.
C             On exit, the leading M-by-N part of this array contains
C             the matrix H * C, if SIDE = 'L', or C * H, if SIDE = 'R'.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N), if SIDE = 'L', or
C                                               (M), if SIDE = 'R'.
C             DWORK is not referenced if H has order less than 11.
C
C     METHOD
C
C     The routine applies the elementary reflector H, taking its special
C     structure into account. The multiplications by the first component
C     of u (which is 1) are avoided, to increase the efficiency.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999.
C     This is a modification of LAPACK Library routine DLARFX.
*
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary matrix operations, elementary reflector, orthogonal
C     transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      DOUBLE PRECISION   TAU
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), DWORK( * ), V( * )
C     ..
C     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   SUM, T1, T2, T3, T4, T5, T6, T7, T8, T9,
     $                   V1, V2, V3, V4, V5, V6, V7, V8, V9
C     ..
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGER
C     ..
C     .. Executable Statements ..
C
      IF( TAU.EQ.ZERO )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
C
C        Form  H * C, where H has order m.
C
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150,
     $           170, 190 ) M
C
C        Code for general M.
C
C        w := C'*u.
C
         CALL DCOPY( N, C, LDC, DWORK, 1 )
         CALL DGEMV( 'Transpose', M-1, N, ONE, C( 2, 1 ), LDC, V, 1,
     $               ONE, DWORK, 1 )
C
C        C := C - tau * u * w'.
C
         CALL DAXPY( N, -TAU, DWORK, 1, C, LDC )
         CALL DGER( M-1, N, -TAU, V, 1, DWORK, 1, C( 2, 1 ), LDC )
         GO TO 410
   10    CONTINUE
C
C        Special code for 1 x 1 Householder.
C
         T1 = ONE - TAU
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
C
C        Special code for 2 x 2 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         DO 40 J = 1, N
            SUM = C( 1, J ) + V1*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
   40    CONTINUE
         GO TO 410
   50    CONTINUE
C
C        Special code for 3 x 3 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 60 J = 1, N
            SUM = C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
   60    CONTINUE
         GO TO 410
   70    CONTINUE
C
C        Special code for 4 x 4 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 80 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
   80    CONTINUE
         GO TO 410
   90    CONTINUE
C
C        Special code for 5 x 5 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 100 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
            C( 5, J ) = C( 5, J ) - SUM*T4
  100    CONTINUE
         GO TO 410
  110    CONTINUE
C
C        Special code for 6 x 6 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 120 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J ) + V5*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
            C( 5, J ) = C( 5, J ) - SUM*T4
            C( 6, J ) = C( 6, J ) - SUM*T5
  120    CONTINUE
         GO TO 410
  130    CONTINUE
C
C        Special code for 7 x 7 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 140 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J ) + V5*C( 6, J ) +
     $            V6*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
            C( 5, J ) = C( 5, J ) - SUM*T4
            C( 6, J ) = C( 6, J ) - SUM*T5
            C( 7, J ) = C( 7, J ) - SUM*T6
  140    CONTINUE
         GO TO 410
  150    CONTINUE
C
C        Special code for 8 x 8 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 160 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J ) + V5*C( 6, J ) +
     $            V6*C( 7, J ) + V7*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
            C( 5, J ) = C( 5, J ) - SUM*T4
            C( 6, J ) = C( 6, J ) - SUM*T5
            C( 7, J ) = C( 7, J ) - SUM*T6
            C( 8, J ) = C( 8, J ) - SUM*T7
  160    CONTINUE
         GO TO 410
  170    CONTINUE
C
C        Special code for 9 x 9 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 180 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J ) + V5*C( 6, J ) +
     $            V6*C( 7, J ) + V7*C( 8, J ) + V8*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*TAU
            C( 2, J ) = C( 2, J ) - SUM*T1
            C( 3, J ) = C( 3, J ) - SUM*T2
            C( 4, J ) = C( 4, J ) - SUM*T3
            C( 5, J ) = C( 5, J ) - SUM*T4
            C( 6, J ) = C( 6, J ) - SUM*T5
            C( 7, J ) = C( 7, J ) - SUM*T6
            C( 8, J ) = C( 8, J ) - SUM*T7
            C( 9, J ) = C( 9, J ) - SUM*T8
  180    CONTINUE
         GO TO 410
  190    CONTINUE
C
C        Special code for 10 x 10 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 200 J = 1, N
            SUM =    C( 1, J ) + V1*C( 2, J ) + V2*C( 3, J ) +
     $            V3*C( 4, J ) + V4*C( 5, J ) + V5*C( 6, J ) +
     $            V6*C( 7, J ) + V7*C( 8, J ) + V8*C( 9, J ) +
     $            V9*C( 10, J )
            C(  1, J ) = C(  1, J ) - SUM*TAU
            C(  2, J ) = C(  2, J ) - SUM*T1
            C(  3, J ) = C(  3, J ) - SUM*T2
            C(  4, J ) = C(  4, J ) - SUM*T3
            C(  5, J ) = C(  5, J ) - SUM*T4
            C(  6, J ) = C(  6, J ) - SUM*T5
            C(  7, J ) = C(  7, J ) - SUM*T6
            C(  8, J ) = C(  8, J ) - SUM*T7
            C(  9, J ) = C(  9, J ) - SUM*T8
            C( 10, J ) = C( 10, J ) - SUM*T9
  200    CONTINUE
         GO TO 410
      ELSE
C
C        Form  C * H, where H has order n.
C
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350,
     $           370, 390 ) N
C
C        Code for general N.
C
C        w := C * u.
C
         CALL DCOPY( M, C, 1, DWORK, 1 )
         CALL DGEMV( 'No transpose', M, N-1, ONE, C( 1, 2 ), LDC, V, 1,
     $               ONE, DWORK, 1 )
C
C        C := C - tau * w * u'.
C
         CALL DAXPY( M, -TAU, DWORK, 1, C, 1 )
         CALL DGER( M, N-1, -TAU, DWORK, 1, V, 1, C( 1, 2 ), LDC )
         GO TO 410
  210    CONTINUE
C
C        Special code for 1 x 1 Householder.
C
         T1 = ONE - TAU
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
C
C        Special code for 2 x 2 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         DO 240 J = 1, M
            SUM = C( J, 1 ) + V1*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
  240    CONTINUE
         GO TO 410
  250    CONTINUE
C
C        Special code for 3 x 3 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 260 J = 1, M
            SUM = C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
  260    CONTINUE
         GO TO 410
  270    CONTINUE
C
C        Special code for 4 x 4 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 280 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
  280    CONTINUE
         GO TO 410
  290    CONTINUE
C
C        Special code for 5 x 5 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 300 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
            C( J, 5 ) = C( J, 5 ) - SUM*T4
  300    CONTINUE
         GO TO 410
  310    CONTINUE
C
C        Special code for 6 x 6 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 320 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 ) + V5*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
            C( J, 5 ) = C( J, 5 ) - SUM*T4
            C( J, 6 ) = C( J, 6 ) - SUM*T5
  320    CONTINUE
         GO TO 410
  330    CONTINUE
C
C        Special code for 7 x 7 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 340 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 ) + V5*C( J, 6 ) +
     $            V6*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
            C( J, 5 ) = C( J, 5 ) - SUM*T4
            C( J, 6 ) = C( J, 6 ) - SUM*T5
            C( J, 7 ) = C( J, 7 ) - SUM*T6
  340    CONTINUE
         GO TO 410
  350    CONTINUE
C
C        Special code for 8 x 8 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 360 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 ) + V5*C( J, 6 ) +
     $            V6*C( J, 7 ) + V7*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
            C( J, 5 ) = C( J, 5 ) - SUM*T4
            C( J, 6 ) = C( J, 6 ) - SUM*T5
            C( J, 7 ) = C( J, 7 ) - SUM*T6
            C( J, 8 ) = C( J, 8 ) - SUM*T7
  360    CONTINUE
         GO TO 410
  370    CONTINUE
C
C        Special code for 9 x 9 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 380 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 ) + V5*C( J, 6 ) +
     $            V6*C( J, 7 ) + V7*C( J, 8 ) + V8*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*TAU
            C( J, 2 ) = C( J, 2 ) - SUM*T1
            C( J, 3 ) = C( J, 3 ) - SUM*T2
            C( J, 4 ) = C( J, 4 ) - SUM*T3
            C( J, 5 ) = C( J, 5 ) - SUM*T4
            C( J, 6 ) = C( J, 6 ) - SUM*T5
            C( J, 7 ) = C( J, 7 ) - SUM*T6
            C( J, 8 ) = C( J, 8 ) - SUM*T7
            C( J, 9 ) = C( J, 9 ) - SUM*T8
  380    CONTINUE
         GO TO 410
  390    CONTINUE
C
C        Special code for 10 x 10 Householder.
C
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 400 J = 1, M
            SUM =    C( J, 1 ) + V1*C( J, 2 ) + V2*C( J, 3 ) +
     $            V3*C( J, 4 ) + V4*C( J, 5 ) + V5*C( J, 6 ) +
     $            V6*C( J, 7 ) + V7*C( J, 8 ) + V8*C( J, 9 ) +
     $            V9*C( J, 10 )
            C( J,  1 ) = C( J,  1 ) - SUM*TAU
            C( J,  2 ) = C( J,  2 ) - SUM*T1
            C( J,  3 ) = C( J,  3 ) - SUM*T2
            C( J,  4 ) = C( J,  4 ) - SUM*T3
            C( J,  5 ) = C( J,  5 ) - SUM*T4
            C( J,  6 ) = C( J,  6 ) - SUM*T5
            C( J,  7 ) = C( J,  7 ) - SUM*T6
            C( J,  8 ) = C( J,  8 ) - SUM*T7
            C( J,  9 ) = C( J,  9 ) - SUM*T8
            C( J, 10 ) = C( J, 10 ) - SUM*T9
  400    CONTINUE
         GO TO 410
      END IF
  410 CONTINUE
      RETURN
C
C *** Last line of MB04PY ***
      END
