      SUBROUTINE MB04OX( N, A, LDA, X, INCX )
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
C     To perform the QR factorization
C
C        (U ) = Q*(R),
C        (x')     (0)
C
C     where U and R are n-by-n upper triangular matrices, x is an
C     n element vector and Q is an (n+1)-by-(n+1) orthogonal matrix.
C
C     U must be supplied in the n-by-n upper triangular part of the
C     array A and this is overwritten by R.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N      (input) INTEGER
C            The number of elements of X and the order of the square
C            matrix A.  N >= 0.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C            On entry, the leading N-by-N upper triangular part of this
C            array must contain the upper triangular matrix U.
C            On exit, the leading N-by-N upper triangular part of this
C            array contains the upper triangular matrix R.
C            The strict lower triangle of A is not referenced.
C
C     LDA    INTEGER
C            The leading dimension of the array A.  LDA >= max(1,N).
C
C     X      (input/output) DOUBLE PRECISION array, dimension
C            (1+(N-1)*INCX)
C            On entry, the incremented array X must contain the
C            vector x. On exit, the content of X is changed.
C
C     INCX   (input) INTEGER.
C            Specifies the increment for the elements of X.  INCX > 0.
C
C     METHOD
C
C     The matrix Q is formed as a sequence of plane rotations in planes
C     (1, n+1), (2, n+1), ..., (n, n+1), the rotation in the (j, n+1)th
C     plane, Q(j), being chosen to annihilate the jth element of x.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine DUTUPD.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION   CI, SI, TEMP
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           DLARTG, DROT
C
C     .. Executable Statements ..
C
C     For efficiency reasons, the parameters are not checked.
C
      IX = 1
C
      DO 20 I = 1, N - 1
         CALL DLARTG( A(I,I), X(IX), CI, SI, TEMP )
         A(I,I) = TEMP
         IX = IX + INCX
         CALL DROT( N-I, A(I,I+1), LDA, X(IX), INCX, CI, SI )
   20 CONTINUE
C
      CALL DLARTG( A(N,N), X(IX), CI, SI, TEMP )
      A(N,N) = TEMP
C
      RETURN
C *** Last line of MB04OX ***
      END
