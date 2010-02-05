      SUBROUTINE MB03QX( N, T, LDT, WR, WI, INFO )
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
C     To compute the eigenvalues of an upper quasi-triangular matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix T.  N >= 0.
C
C     T       (input) DOUBLE PRECISION array, dimension(LDT,N)
C             The upper quasi-triangular matrix T.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= max(1,N).
C
C     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N)
C             The real and imaginary parts, respectively, of the
C             eigenvalues of T. The eigenvalues are stored in the same
C             order as on the diagonal of T. If T(i:i+1,i:i+1) is a
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
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     March 1998. Based on the RASP routine SEIG.
C
C     ******************************************************************
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER        ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER          INFO, LDT, N
C     .. Array Arguments ..
      DOUBLE PRECISION T(LDT, *), WI(*), WR(*)
C     .. Local Scalars ..
      INTEGER          I, I1, INEXT
      DOUBLE PRECISION A11, A12, A21, A22, CS, SN
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
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03QX', -INFO )
         RETURN
      END IF
C
      INEXT = 1
      DO 10 I = 1, N
         IF( I.LT.INEXT )
     $      GO TO 10
         IF( I.NE.N ) THEN
            IF( T(I+1,I).NE.ZERO ) THEN
C
C              A pair of eigenvalues.
C
               INEXT = I + 2
               I1 = I + 1
               A11 = T(I,I)
               A12 = T(I,I1)
               A21 = T(I1,I)
               A22 = T(I1,I1)
               CALL DLANV2( A11, A12, A21, A22, WR(I), WI(I), WR(I1),
     $                      WI(I1), CS, SN )
               GO TO 10
            END IF
         END IF
C
C        Simple eigenvalue.
C
         INEXT = I + 1
         WR(I) = T(I,I)
         WI(I) = ZERO
   10 CONTINUE
C
      RETURN
C *** Last line of MB03QX ***
      END
