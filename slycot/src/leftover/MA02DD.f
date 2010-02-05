      SUBROUTINE MA02DD( JOB, UPLO, N, A, LDA, AP )
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
C     To pack/unpack the upper or lower triangle of a symmetric matrix.
C     The packed matrix is stored column-wise in the one-dimensional
C     array AP.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies whether the matrix should be packed or unpacked,
C             as follows:
C             = 'P':  The matrix should be packed;
C             = 'U':  The matrix should be unpacked.
C
C     UPLO    CHARACTER*1
C             Specifies the part of the matrix to be packed/unpacked,
C             as follows:
C             = 'U':  Upper triangular part;
C             = 'L':  Lower triangular part.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input or output) DOUBLE PRECISION array, dimension
C             (LDA,N)
C             This array is an input parameter if JOB = 'P', and an
C             output parameter if JOB = 'U'.
C             On entry, if JOB = 'P', the leading N-by-N upper
C             triangular part (if UPLO = 'U'), or lower triangular part
C             (if UPLO = 'L'), of this array must contain the
C             corresponding upper or lower triangle of the symmetric
C             matrix A, and the other strictly triangular part is not
C             referenced.
C             On exit, if JOB = 'U', the leading N-by-N upper triangular
C             part (if UPLO = 'U'), or lower triangular part (if
C             UPLO = 'L'), of this array contains the corresponding
C             upper or lower triangle of the symmetric matrix A; the
C             other strictly triangular part is not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     AP      (output or input) DOUBLE PRECISION array, dimension
C             (N*(N+1)/2)
C             This array is an output parameter if JOB = 'P', and an
C             input parameter if JOB = 'U'.
C             On entry, if JOB = 'U', the leading N*(N+1)/2 elements of
C             this array must contain the upper (if UPLO = 'U') or lower
C             (if UPLO = 'L') triangle of the symmetric matrix A, packed
C             column-wise. That is, the elements are stored in the order
C             11, 12, 22, ..., 1n, 2n, 3n, ..., nn,      if UPLO = 'U';
C             11, 21, 31, ..., n1, 22, 32, ..., n2, ..., if UPLO = 'L'.
C             On exit, if JOB = 'P', the leading N*(N+1)/2 elements of
C             this array contain the upper (if UPLO = 'U') or lower
C             (if UPLO = 'L') triangle of the symmetric matrix A, packed
C             column-wise, as described above.
C
C     CONTRIBUTOR
C
C     V. Sima, Research Institute for Informatics, Bucharest, Romania,
C     Oct. 1998.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER          JOB, UPLO
      INTEGER            LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), AP(*)
C     .. Local Scalars ..
      LOGICAL            LUPLO
      INTEGER            IJ, J
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DCOPY
C
C     .. Executable Statements ..
C
C     For efficiency reasons, the parameters are not checked for errors.
C
      LUPLO = LSAME( UPLO, 'L' )
      IJ = 1
      IF( LSAME( JOB, 'P' ) ) THEN
         IF( LUPLO ) THEN
C
C           Pack the lower triangle of A.
C
            DO 20 J = 1, N
               CALL DCOPY( N-J+1, A(J,J), 1, AP(IJ), 1 )
               IJ = IJ + N - J + 1
   20       CONTINUE
C
         ELSE
C
C           Pack the upper triangle of A.
C
            DO 40 J = 1, N
               CALL DCOPY( J, A(1,J), 1, AP(IJ), 1 )
               IJ = IJ + J
   40       CONTINUE
C
         END IF
      ELSE
         IF( LUPLO ) THEN
C
C           Unpack the lower triangle of A.
C
            DO 60 J = 1, N
               CALL DCOPY( N-J+1, AP(IJ), 1, A(J,J), 1 )
               IJ = IJ + N - J + 1
   60       CONTINUE
C
         ELSE
C
C           Unpack the upper triangle of A.
C
            DO 80 J = 1, N
               CALL DCOPY( J, AP(IJ), 1, A(1,J), 1 )
               IJ = IJ + J
   80       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of MA02DD ***
      END
