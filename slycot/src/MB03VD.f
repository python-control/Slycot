      SUBROUTINE MB03VD( N, P, ILO, IHI, A, LDA1, LDA2, TAU, LDTAU,
     $                   DWORK, INFO )
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
C     To reduce a product of p real general matrices A = A_1*A_2*...*A_p
C     to upper Hessenberg form, H = H_1*H_2*...*H_p, where H_1 is
C     upper Hessenberg, and H_2, ..., H_p are upper triangular, by using
C     orthogonal similarity transformations on A,
C
C             Q_1' * A_1 * Q_2 = H_1,
C             Q_2' * A_2 * Q_3 = H_2,
C                    ...
C             Q_p' * A_p * Q_1 = H_p.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the square matrices A_1, A_2, ..., A_p.
C             N >= 0.
C
C     P       (input) INTEGER
C             The number of matrices in the product A_1*A_2*...*A_p.
C             P >= 1.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that all matrices A_j, j = 2, ..., p, are
C             already upper triangular in rows and columns 1:ILO-1 and
C             IHI+1:N, and A_1 is upper Hessenberg in rows and columns
C             1:ILO-1 and IHI+1:N, with A_1(ILO,ILO-1) = 0 (unless
C             ILO = 1), and A_1(IHI+1,IHI) = 0 (unless IHI = N).
C             If this is not the case, ILO and IHI should be set to 1
C             and N, respectively.
C             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C             (LDA1,LDA2,P)
C             On entry, the leading N-by-N-by-P part of this array must
C             contain the matrices of factors to be reduced;
C             specifically, A(*,*,j) must contain A_j, j = 1, ..., p.
C             On exit, the leading N-by-N upper triangle and the first
C             subdiagonal of A(*,*,1) contain the upper Hessenberg
C             matrix H_1, and the elements below the first subdiagonal,
C             with the first column of the array TAU represent the
C             orthogonal matrix Q_1 as a product of elementary
C             reflectors. See FURTHER COMMENTS.
C             For j > 1, the leading N-by-N upper triangle of A(*,*,j)
C             contains the upper triangular matrix H_j, and the elements
C             below the diagonal, with the j-th column of the array TAU
C             represent the orthogonal matrix Q_j as a product of
C             elementary reflectors. See FURTHER COMMENTS.
C
C     LDA1    INTEGER
C             The first leading dimension of the array A.
C             LDA1 >= max(1,N).
C
C     LDA2    INTEGER
C             The second leading dimension of the array A.
C             LDA2 >= max(1,N).
C
C     TAU     (output) DOUBLE PRECISION array, dimension (LDTAU,P)
C             The leading N-1 elements in the j-th column contain the
C             scalar factors of the elementary reflectors used to form
C             the matrix Q_j, j = 1, ..., P. See FURTHER COMMENTS.
C
C     LDTAU   INTEGER
C             The leading dimension of the array TAU.
C             LDTAU >= max(1,N-1).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
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
C     The algorithm consists in ihi-ilo major steps. In each such
C     step i, ilo <= i <= ihi-1, the subdiagonal elements in the i-th
C     column of A_j are annihilated using a Householder transformation
C     from the left, which is also applied to A_(j-1) from the right,
C     for j = p:-1:2. Then, the elements below the subdiagonal of the
C     i-th column of A_1 are annihilated, and the Householder
C     transformation is also applied to A_p from the right.
C     See FURTHER COMMENTS.
C
C     REFERENCES
C
C     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P.
C         The periodic Schur decomposition: algorithms and applications.
C         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42,
C         1992.
C
C     [2] Sreedhar, J. and Van Dooren, P.
C         Periodic Schur form and some matrix equations.
C         Proc. of the Symposium on the Mathematical Theory of Networks
C         and Systems (MTNS'93), Regensburg, Germany (U. Helmke,
C         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable.
C
C     FURTHER COMMENTS
C
C     Each matrix Q_j is represented as a product of (ihi-ilo)
C     elementary reflectors,
C
C        Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1).
C
C     Each H_j(i), i = ilo, ..., ihi-1, has the form
C
C        H_j(i) = I - tau_j * v_j * v_j',
C
C     where tau_j is a real scalar, and v_j is a real vector with
C     v_j(1:i) = 0, v_j(i+1) = 1 and v_j(ihi+1:n) = 0; v_j(i+2:ihi)
C     is stored on exit in A_j(i+2:ihi,i), and tau_j in TAU(i,j).
C
C     The contents of A_1 are illustrated by the following example
C     for n = 7, ilo = 2, and ihi = 6:
C
C     on entry                         on exit
C
C     ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a )
C     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   a )
C     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h )
C     ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a )
C
C     where a denotes an element of the original matrix A_1, h denotes
C     a modified element of the upper Hessenberg matrix H_1, and vi
C     denotes an element of the vector defining H_1(i).
C
C     The contents of A_j, j > 1, are illustrated by the following
C     example for n = 7, ilo = 2, and ihi = 6:
C
C     on entry                         on exit
C
C     ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a )
C     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h )
C     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  v5  h   h )
C     ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a )
C
C     where a denotes an element of the original matrix A_j, h denotes
C     a modified element of the upper triangular matrix H_j, and vi
C     denotes an element of the vector defining H_j(i). (The element
C     (1,2) in A_p is also unchanged for this example.)
C
C     Note that for P = 1, the LAPACK Library routine DGEHRD could be
C     more efficient on some computer architectures than this routine
C     (a BLAS 2 version).
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga,
C     German Aerospace Center, DLR Oberpfaffenhofen, February 1999.
C     Partly based on the routine PSHESS by A. Varga
C     (DLR Oberpfaffenhofen), November 26, 1995.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Hessenberg form, orthogonal transformation, periodic systems,
C     similarity transformation, triangular form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER           IHI, ILO, INFO, LDA1, LDA2, LDTAU, N, P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA1, LDA2, * ), DWORK( * ), TAU( LDTAU, * )
C     ..
C     .. Local Scalars ..
      INTEGER           I, I1, I2, J, NH
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY( 1 )
C     ..
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLARFG, MB04PY, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( P.LT.1 ) THEN
         INFO = -2
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -4
      ELSE IF( LDA1.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDA2.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDTAU.LT.MAX( 1, N-1 ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03VD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      NH = IHI - ILO + 1
      IF ( NH.LE.1 )
     $   RETURN
C
      DUMMY( 1 ) = ZERO
C
      DO 20 I = ILO, IHI - 1
         I1 = I + 1
         I2 = MIN( I+2, N )
C
         DO 10 J = P, 2, -1
C
C           Set the elements 1:ILO-1 and IHI:N-1 of TAU(*,J) to zero.
C
            CALL DCOPY( ILO-1, DUMMY, 0, TAU( 1, J ), 1 )
            IF ( IHI.LT.N )
     $         CALL DCOPY( N-IHI, DUMMY, 0, TAU( IHI, J ), 1 )
C
C           Compute elementary reflector H_j(i) to annihilate
C           A_j(i+1:ihi,i).
C
            CALL DLARFG( IHI-I+1, A( I, I, J ), A( I1, I, J ), 1,
     $                   TAU( I, J ) )
C
C           Apply H_j(i) to A_(j-1)(1:ihi,i:ihi) from the right.
C
            CALL MB04PY( 'Right', IHI, IHI-I+1, A( I1, I, J ),
     $                   TAU( I, J ), A( 1, I, J-1 ), LDA1, DWORK )
C
C           Apply H_j(i) to A_j(i:ihi,i+1:n) from the left.
C
            CALL MB04PY( 'Left', IHI-I+1, N-I, A( I1, I, J ),
     $                   TAU( I, J ), A( I, I1, J ), LDA1, DWORK )
   10    CONTINUE
C
C        Compute elementary reflector H_1(i) to annihilate
C        A_1(i+2:ihi,i).
C
         CALL DLARFG( IHI-I, A( I1, I, 1 ), A( I2, I, 1 ), 1,
     $                TAU( I, 1 ) )
C
C        Apply H_1(i) to A_p(1:ihi,i+1:ihi) from the right.
C
         CALL MB04PY( 'Right', IHI, IHI-I, A( I2, I, 1 ), TAU( I, 1 ),
     $                A( 1, I1, P ), LDA1, DWORK )
C
C        Apply H_1(i) to A_1(i+1:ihi,i+1:n) from the left.
C
         CALL MB04PY( 'Left', IHI-I, N-I, A( I2, I, 1 ), TAU( I, 1 ),
     $                A( I1, I1, 1 ), LDA1, DWORK )
   20 CONTINUE
C
      RETURN
C
C *** Last line of MB03VD ***
      END
