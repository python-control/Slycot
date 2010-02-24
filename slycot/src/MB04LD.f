      SUBROUTINE MB04LD( UPLO, N, M, P, L, LDL, A, LDA, B, LDB, C, LDC,
     $                   TAU, DWORK )
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
C     To calculate an LQ factorization of the first block row and apply
C     the orthogonal transformations (from the right) also to the second
C     block row of a structured matrix, as follows
C                        _
C        [ L   A ]     [ L   0 ]
C        [       ]*Q = [       ]
C        [ 0   B ]     [ C   D ]
C                 _
C     where L and L are lower triangular. The matrix A can be full or
C     lower trapezoidal/triangular. The problem structure is exploited.
C     This computation is useful, for instance, in combined measurement
C     and time update of one iteration of the Kalman filter (square
C     root covariance filter).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPLO    CHARACTER*1
C             Indicates if the matrix A is or not triangular as follows:
C             = 'L':  Matrix A is lower trapezoidal/triangular;
C             = 'F':  Matrix A is full.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER                 _
C             The order of the matrices L and L.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrices A, B and D.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of the matrices B, C and D.  P >= 0.
C
C     L       (input/output) DOUBLE PRECISION array, dimension (LDL,N)
C             On entry, the leading N-by-N lower triangular part of this
C             array must contain the lower triangular matrix L.
C             On exit, the leading N-by-N lower triangular part of this
C                                                        _
C             array contains the lower triangular matrix L.
C             The strict upper triangular part of this array is not
C             referenced.
C
C     LDL     INTEGER
C             The leading dimension of array L.  LDL >= MAX(1,N).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
C             On entry, if UPLO = 'F', the leading N-by-M part of this
C             array must contain the matrix A. If UPLO = 'L', the
C             leading N-by-MIN(N,M) part of this array must contain the
C             lower trapezoidal (lower triangular if N <= M) matrix A,
C             and the elements above the diagonal are not referenced.
C             On exit, the leading N-by-M part (lower trapezoidal or
C             triangular, if UPLO = 'L') of this array contains the
C             trailing components (the vectors v, see Method) of the
C             elementary reflectors used in the factorization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading P-by-M part of this array must
C             contain the matrix B.
C             On exit, the leading P-by-M part of this array contains
C             the computed matrix D.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,P).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array contains the
C             computed matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N)
C             The scalar factors of the elementary reflectors used.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     METHOD
C
C     The routine uses N Householder transformations exploiting the zero
C     pattern of the block matrix.  A Householder matrix has the form
C
C                                     ( 1 ),
C        H  = I - tau *u *u',    u  = ( v )
C         i          i  i  i      i   (  i)
C
C     where v  is an M-vector, if UPLO = 'F', or an min(i,M)-vector, if
C            i
C     UPLO = 'L'.  The components of v  are stored in the i-th row of A,
C                                     i
C     and tau  is stored in TAU(i).
C            i
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary reflector, LQ factorization, orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         UPLO
      INTEGER           LDA, LDB, LDC, LDL, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*),
     $                  L(LDL,*), TAU(*)
C     .. Local Scalars ..
      LOGICAL           LUPLO
      INTEGER           I, IM
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DGER, DLARFG, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      IF( MIN( M, N ).EQ.0 )
     $   RETURN
C
      LUPLO = LSAME( UPLO, 'L' )
      IM = M
C
      DO 10 I = 1, N
C
C        Annihilate the I-th row of A and apply the transformations to
C        the entire block matrix, exploiting its structure.
C
         IF( LUPLO ) IM = MIN( I, M )
         CALL DLARFG( IM+1, L(I,I), A(I,1), LDA, TAU(I) )
         IF( TAU(I).NE.ZERO ) THEN
C
C           [    w   ]    [ L(I+1:N,I) A(I+1:N,1:IM) ]   [ 1 ]
C           [        ] := [                          ] * [   ]
C           [ C(:,I) ]    [      0        B(:,1:IM)  ]   [ v ]
C
            IF( I.LT.N ) THEN
               CALL DCOPY( N-I, L(I+1,I), 1, DWORK, 1 )
               CALL DGEMV( 'No transpose', N-I, IM, ONE, A(I+1,1), LDA,
     $                     A(I,1), LDA, ONE, DWORK, 1 )
            END IF
            CALL DGEMV( 'No transpose', P, IM, ONE, B, LDB, A(I,1),
     $                  LDA, ZERO, C(1,I), 1 )
C
C           [ L(I+1:N,I) A(I+1:N,1:IM) ]    [ L(I+1:N,I) A(I+1:N,1:IM) ]
C           [                          ] := [                          ]
C           [   C(:,I)      D(:,1:IM)  ]    [      0        B(:,1:IM)  ]
C
C                                                 [    w   ]
C                                         - tau * [        ] * [ 1 , v']
C                                                 [ C(:,I) ]
C
            IF( I.LT.N ) THEN
               CALL DAXPY( N-I, -TAU(I), DWORK, 1, L(I+1,I), 1 )
               CALL DGER( N-I, IM, -TAU(I), DWORK, 1, A(I,1), LDA,
     $                    A(I+1,1), LDA )
            END IF
            CALL DSCAL( P, -TAU(I), C(1,I), 1 )
            CALL DGER( P, IM, ONE, C(1,I), 1, A(I,1), LDA, B, LDB )
         END IF
   10 CONTINUE
C
      RETURN
C *** Last line of MB04LD ***
      END
