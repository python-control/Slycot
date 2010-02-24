      SUBROUTINE MB04OW( M, N, P, A, LDA, T, LDT, X, INCX, B, LDB,
     $                   C, LDC, D, INCD )
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
C        ( U  ) = Q*( R ),  where  U = ( U1  U2 ),  R = ( R1  R2 ),
C        ( x' )     ( 0 )              ( 0   T  )       ( 0   R3 )
C
C     where U and R are (m+n)-by-(m+n) upper triangular matrices, x is
C     an m+n element vector, U1 is m-by-m, T is n-by-n, stored
C     separately, and Q is an (m+n+1)-by-(m+n+1) orthogonal matrix.
C
C     The matrix ( U1 U2 ) must be supplied in the m-by-(m+n) upper
C     trapezoidal part of the array A and this is overwritten by the
C     corresponding part ( R1 R2 ) of R. The remaining upper triangular
C     part of R, R3, is overwritten on the array T.
C
C     The transformations performed are also applied to the (m+n+1)-by-p
C     matrix ( B' C' d )' (' denotes transposition), where B, C, and d'
C     are m-by-p, n-by-p, and 1-by-p matrices, respectively.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M      (input) INTEGER
C            The number of rows of the matrix ( U1  U2 ).  M >= 0.
C
C     N      (input) INTEGER
C            The order of the matrix T.  N >= 0.
C
C     P      (input) INTEGER
C            The number of columns of the matrices B and C.  P >= 0.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C            On entry, the leading M-by-(M+N) upper trapezoidal part of
C            this array must contain the upper trapezoidal matrix
C            ( U1 U2 ).
C            On exit, the leading M-by-(M+N) upper trapezoidal part of
C            this array contains the upper trapezoidal matrix ( R1 R2 ).
C            The strict lower triangle of A is not referenced.
C
C     LDA    INTEGER
C            The leading dimension of the array A.  LDA >= max(1,M).
C
C     T      (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C            On entry, the leading N-by-N upper triangular part of this
C            array must contain the upper triangular matrix T.
C            On exit, the leading N-by-N upper triangular part of this
C            array contains the upper triangular matrix R3.
C            The strict lower triangle of T is not referenced.
C
C     LDT    INTEGER
C            The leading dimension of the array T.  LDT >= max(1,N).
C
C     X      (input/output) DOUBLE PRECISION array, dimension
C            (1+(M+N-1)*INCX), if M+N > 0, or dimension (0), if M+N = 0.
C            On entry, the incremented array X must contain the
C            vector x. On exit, the content of X is changed.
C
C     INCX   (input) INTEGER
C            Specifies the increment for the elements of X.  INCX > 0.
C
C     B      (input/output) DOUBLE PRECISION array, dimension (LDB,P)
C            On entry, the leading M-by-P part of this array must
C            contain the matrix B.
C            On exit, the leading M-by-P part of this array contains
C            the transformed matrix B.
C            If M = 0 or P = 0, this array is not referenced.
C
C     LDB    INTEGER
C            The leading dimension of the array B.
C            LDB >= max(1,M), if P > 0;
C            LDB >= 1,        if P = 0.
C
C     C      (input/output) DOUBLE PRECISION array, dimension (LDC,P)
C            On entry, the leading N-by-P part of this array must
C            contain the matrix C.
C            On exit, the leading N-by-P part of this array contains
C            the transformed matrix C.
C            If N = 0 or P = 0, this array is not referenced.
C
C     LDC    INTEGER
C            The leading dimension of the array C.
C            LDC >= max(1,N), if P > 0;
C            LDC >= 1,        if P = 0.
C
C     D      (input/output) DOUBLE PRECISION array, dimension
C            (1+(P-1)*INCD), if P > 0, or dimension (0), if P = 0.
C            On entry, the incremented array D must contain the
C            vector d.
C            On exit, this incremented array contains the transformed
C            vector d.
C            If P = 0, this array is not referenced.
C
C     INCD   (input) INTEGER
C            Specifies the increment for the elements of D.  INCD > 0.
C
C     METHOD
C
C     Let q = m+n. The matrix Q is formed as a sequence of plane
C     rotations in planes (1, q+1), (2, q+1), ..., (q, q+1), the
C     rotation in the (j, q+1)th plane, Q(j), being chosen to
C     annihilate the jth element of x.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0((M+N)*(M+N+P)) operations and is backward
C     stable.
C
C     FURTHER COMMENTS
C
C     For P = 0, this routine produces the same result as SLICOT Library
C     routine MB04OX, but matrix T may not be stored in the array A.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Matrix operations, plane rotations.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER            INCD, INCX, LDA, LDB, LDC, LDT, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), C(LDC,*), D(*), T(LDT,*),
     $                   X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION   CI, SI, TEMP
      INTEGER            I, IX, MN
C     .. External Subroutines ..
      EXTERNAL           DLARTG, DROT
C
C     .. Executable Statements ..
C
C     For efficiency reasons, the parameters are not checked.
C
      MN = M + N
      IF ( INCX.GT.1 ) THEN
C
C        Code for increment INCX > 1.
C
         IX = 1
         IF ( M.GT.0 ) THEN
C
            DO 10 I = 1, M - 1
               CALL DLARTG( A(I,I), X(IX), CI, SI, TEMP )
               A(I,I) = TEMP
               IX = IX + INCX
               CALL DROT( MN-I, A(I,I+1), LDA, X(IX), INCX, CI, SI )
               IF ( P.GT.0 )
     $            CALL DROT( P, B(I,1), LDB, D, INCD, CI, SI )
   10       CONTINUE
C
            CALL DLARTG( A(M,M), X(IX), CI, SI, TEMP )
            A(M,M) = TEMP
            IX = IX + INCX
            IF ( N.GT.0 )
     $         CALL DROT( N, A(M,M+1), LDA, X(IX), INCX, CI, SI )
            IF ( P.GT.0 )
     $         CALL DROT( P, B(M,1), LDB, D, INCD, CI, SI )
         END IF
C
         IF ( N.GT.0 ) THEN
C
            DO 20 I = 1, N - 1
               CALL DLARTG( T(I,I), X(IX), CI, SI, TEMP )
               T(I,I) = TEMP
               IX = IX + INCX
               CALL DROT( N-I, T(I,I+1), LDT, X(IX), INCX, CI, SI )
               IF ( P.GT.0 )
     $            CALL DROT( P, C(I,1), LDC, D, INCD, CI, SI )
   20       CONTINUE
C
            CALL DLARTG( T(N,N), X(IX), CI, SI, TEMP )
            T(N,N) = TEMP
            IF ( P.GT.0 )
     $         CALL DROT( P, C(N,1), LDC, D, INCD, CI, SI )
         END IF
C
      ELSEIF ( INCX.EQ.1 ) THEN
C
C        Code for increment INCX = 1.
C
         IF ( M.GT.0 ) THEN
C
            DO 30 I = 1, M - 1
               CALL DLARTG( A(I,I), X(I), CI, SI, TEMP )
               A(I,I) = TEMP
               CALL DROT( MN-I, A(I,I+1), LDA, X(I+1), 1, CI, SI )
               IF ( P.GT.0 )
     $            CALL DROT( P, B(I,1), LDB, D, INCD, CI, SI )
   30       CONTINUE
C
            CALL DLARTG( A(M,M), X(M), CI, SI, TEMP )
            A(M,M) = TEMP
            IF ( N.GT.0 )
     $         CALL DROT( N, A(M,M+1), LDA, X(M+1), 1, CI, SI )
            IF ( P.GT.0 )
     $         CALL DROT( P, B(M,1), LDB, D, INCD, CI, SI )
         END IF
C
         IF ( N.GT.0 ) THEN
            IX = M + 1
C
            DO 40 I = 1, N - 1
               CALL DLARTG( T(I,I), X(IX), CI, SI, TEMP )
               T(I,I) = TEMP
               IX = IX + 1
               CALL DROT( N-I, T(I,I+1), LDT, X(IX), 1, CI, SI )
               IF ( P.GT.0 )
     $            CALL DROT( P, C(I,1), LDC, D, INCD, CI, SI )
   40       CONTINUE
C
            CALL DLARTG( T(N,N), X(IX), CI, SI, TEMP )
            T(N,N) = TEMP
            IF ( P.GT.0 )
     $         CALL DROT( P, C(N,1), LDC, D, INCD, CI, SI )
         END IF
      END IF
C
      RETURN
C *** Last line of MB04OW ***
      END
