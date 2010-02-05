      SUBROUTINE MB04GD( M, N, A, LDA, JPVT, TAU, DWORK, INFO )
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
C     To compute an RQ factorization with row pivoting of a
C     real m-by-n matrix A: P*A = R*Q.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the m-by-n matrix A.
C             On exit,
C             if m <= n, the upper triangle of the subarray
C             A(1:m,n-m+1:n) contains the m-by-m upper triangular
C             matrix R;
C             if m >= n, the elements on and above the (m-n)-th
C             subdiagonal contain the m-by-n upper trapezoidal matrix R;
C             the remaining elements, with the array TAU, represent the
C             orthogonal matrix Q as a product of min(m,n) elementary
C             reflectors (see METHOD).
C
C     LDA     INTEGER
C             The leading dimension of the array A. LDA >= max(1,M).
C
C     JPVT    (input/output) INTEGER array, dimension (M)
C             On entry, if JPVT(i) .ne. 0, the i-th row of A is permuted
C             to the bottom of P*A (a trailing row); if JPVT(i) = 0,
C             the i-th row of A is a free row.
C             On exit, if JPVT(i) = k, then the i-th row of P*A
C             was the k-th row of A.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
C             The scalar factors of the elementary reflectors.
C
C     Workspace
C
C     DWORK    DOUBLE PRECISION array, dimension (3*M)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The matrix Q is represented as a product of elementary reflectors
C
C        Q = H(1) H(2) . . . H(k), where k = min(m,n).
C
C     Each H(i) has the form
C
C        H = I - tau * v * v'
C
C     where tau is a real scalar, and v is a real vector with
C     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit
C     in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
C
C     The matrix P is represented in jpvt as follows: If
C        jpvt(j) = i
C     then the jth row of P is the ith canonical unit vector.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997.
C     Based on LAPACK Library routines DGEQPF and DGERQ2.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Factorization, matrix algebra, matrix operations, orthogonal
C     transformation, triangular form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, P05
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, P05 = 0.05D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
C     ..
C     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), DWORK( * ), TAU( * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, ITEMP, J, K, MA, MKI, NFREE, NKI, PVT
      DOUBLE PRECISION   AII, TEMP, TEMP2
C     ..
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DNRM2
      EXTERNAL           DNRM2, IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGERQ2, DLARF, DLARFG, DORMR2, DSWAP, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04GD', -INFO )
         RETURN
      END IF
C
      K = MIN( M, N )
C
C     Move non-free rows bottom.
C
      ITEMP = M
      DO 10 I = M, 1, -1
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL DSWAP( N, A( I, 1 ), LDA, A( ITEMP, 1 ), LDA )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP - 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      NFREE = M - ITEMP
C
C     Compute the RQ factorization and update remaining rows.
C
      IF( NFREE.GT.0 ) THEN
         MA = MIN( NFREE, N )
         CALL DGERQ2( MA, N, A(M-MA+1,1), LDA, TAU(K-MA+1), DWORK,
     $                INFO )
         CALL DORMR2( 'Right', 'Transpose', M-MA, N, MA, A(M-MA+1,1),
     $                LDA, TAU(K-MA+1), A, LDA, DWORK, INFO )
      END IF
C
      IF( NFREE.LT.K ) THEN
C
C        Initialize partial row norms. The first ITEMP elements of
C        DWORK store the exact row norms. (Here, ITEMP is the number of
C        free rows, which have been permuted to be the first ones.)
C
         DO 20 I = 1, ITEMP
            DWORK( I ) = DNRM2( N-NFREE, A( I, 1 ), LDA )
            DWORK( M+I ) = DWORK( I )
   20    CONTINUE
C
C        Compute factorization.
C
         DO 40 I = K-NFREE, 1, -1
C
C           Determine ith pivot row and swap if necessary.
C
            MKI = M - K + I
            NKI = N - K + I
            PVT = IDAMAX( MKI, DWORK, 1 )
C
            IF( PVT.NE.MKI ) THEN
               CALL DSWAP( N, A( PVT, 1 ), LDA, A( MKI, 1 ), LDA )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( MKI )
               JPVT( MKI ) = ITEMP
               DWORK( PVT )   = DWORK( MKI )
               DWORK( M+PVT ) = DWORK( M+MKI )
            END IF
C
C           Generate elementary reflector H(i) to annihilate
C           A(m-k+i,1:n-k+i-1), k = min(m,n).
C
            CALL DLARFG( NKI, A( MKI, NKI ), A( MKI, 1 ), LDA, TAU( I )
     $                 )
C
C           Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right.
C
            AII = A( MKI, NKI )
            A( MKI, NKI ) = ONE
            CALL DLARF( 'Right', MKI-1, NKI, A( MKI, 1 ), LDA,
     $                  TAU( I ), A, LDA, DWORK( 2*M+1 ) )
            A( MKI, NKI ) = AII
C
C           Update partial row norms.
C
            DO 30 J = 1, MKI - 1
               IF( DWORK( J ).NE.ZERO ) THEN
                  TEMP = ONE - ( ABS( A( J, NKI ) ) / DWORK( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = ONE + P05*TEMP*
     $                    ( DWORK( J ) / DWORK( M+J ) )**2
                  IF( TEMP2.EQ.ONE ) THEN
                     DWORK( J ) = DNRM2( NKI-1, A( J, 1 ), LDA )
                     DWORK( M+J ) = DWORK( J )
                  ELSE
                     DWORK( J ) = DWORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
C
   40    CONTINUE
      END IF
C
      RETURN
C *** Last line of MB04GD ***
      END
