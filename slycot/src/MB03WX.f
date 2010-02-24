      SUBROUTINE MB03WX( N, P, T, LDT1, LDT2, WR, WI, INFO )
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
C     To compute the eigenvalues of a product of matrices,
C     T = T_1*T_2*...*T_p, where T_1 is an upper quasi-triangular
C     matrix and T_2, ..., T_p are upper triangular matrices.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix T.  N >= 0.
C
C     P       (input) INTEGER
C             The number of matrices in the product T_1*T_2*...*T_p.
C             P >= 1.
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT1,LDT2,P)
C             The leading N-by-N part of T(*,*,1) must contain the upper
C             quasi-triangular matrix T_1 and the leading N-by-N part of
C             T(*,*,j) for j > 1 must contain the upper-triangular
C             matrix T_j, j = 2, ..., p.
C             The elements below the subdiagonal of T(*,*,1) and below
C             the diagonal of T(*,*,j), j = 2, ..., p, are not
C             referenced.
C
C     LDT1    INTEGER
C             The first leading dimension of the array T.
C             LDT1 >= max(1,N).
C
C     LDT2    INTEGER
C             The second leading dimension of the array T.
C             LDT2 >= max(1,N).
C
C     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N)
C             The real and imaginary parts, respectively, of the
C             eigenvalues of T. The eigenvalues are stored in the same
C             order as on the diagonal of T_1. If T(i:i+1,i:i+1,1) is a
C             2-by-2 diagonal block with complex conjugated eigenvalues
C             then WI(i) > 0 and WI(i+1) = -WI(i).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalue, eigenvalue decomposition, periodic systems,
C     real Schur form, triangular form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      INTEGER          INFO, LDT1, LDT2, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION T( LDT1, LDT2, * ), WI( * ), WR( * )
C     .. Local Scalars ..
      INTEGER          I, I1, INEXT, J
      DOUBLE PRECISION A11, A12, A21, A22, CS, SN, T11, T12, T22
C     .. External Subroutines ..
      EXTERNAL         DLANV2, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        MAX
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( P.LT.1 ) THEN
         INFO = -2
      ELSE IF( LDT1.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDT2.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03WX', -INFO )
         RETURN
      END IF
C
      INEXT = 1
      DO 30 I = 1, N
         IF( I.LT.INEXT )
     $      GO TO 30
         IF( I.NE.N ) THEN
            IF( T( I+1, I, 1 ).NE.ZERO ) THEN
C
C              A pair of eigenvalues. First compute the corresponding
C              elements of T(I:I+1,I:I+1).
C
               INEXT = I + 2
               I1  = I + 1
               T11 = ONE
               T12 = ZERO
               T22 = ONE
C
               DO 10 J = 2, P
                  T22 = T22*T( I1, I1, J )
                  T12 = T11*T( I,  I1, J ) + T12*T( I1, I1, J )
                  T11 = T11*T( I,  I,  J )
   10          CONTINUE
C
               A11 = T( I,  I, 1 )*T11
               A12 = T( I,  I, 1 )*T12 + T( I,  I1, 1 )*T22
               A21 = T( I1, I, 1 )*T11
               A22 = T( I1, I, 1 )*T12 + T( I1, I1, 1 )*T22
C
               CALL DLANV2( A11, A12, A21, A22, WR( I ), WI( I ),
     $                      WR( I1 ), WI( I1 ), CS, SN )
               GO TO 30
            END IF
         END IF
C
C        Simple eigenvalue. Compute the corresponding element of T(I,I).
C
         INEXT = I + 1
         T11 = ONE
C
         DO 20 J = 1, P
            T11 = T11*T( I, I, J )
   20    CONTINUE
C
         WR( I ) = T11
         WI( I ) = ZERO
   30 CONTINUE
C
      RETURN
C *** Last line of MB03WX ***
      END
