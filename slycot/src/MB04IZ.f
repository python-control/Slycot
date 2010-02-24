      SUBROUTINE MB04IZ( N, M, P, L, A, LDA, B, LDB, TAU, ZWORK, LZWORK,
     $                   INFO )
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
C     To compute a QR factorization of an n-by-m matrix A (A = Q * R),
C     having a p-by-min(p,m) zero triangle in the lower left-hand side
C     corner, as shown below, for n = 8, m = 7, and p = 2:
C
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C            [ x x x x x x x ]
C        A = [ x x x x x x x ],
C            [ x x x x x x x ]
C            [ 0 x x x x x x ]
C            [ 0 0 x x x x x ]
C
C     and optionally apply the transformations to an n-by-l matrix B
C     (from the left). The problem structure is exploited. This
C     computation is useful, for instance, in combined measurement and
C     time update of one iteration of the time-invariant Kalman filter
C     (square root information filter).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of rows of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrix A.  M >= 0.
C
C     P       (input) INTEGER
C             The order of the zero triagle.  P >= 0.
C
C     L       (input) INTEGER
C             The number of columns of the matrix B.  L >= 0.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA,M)
C             On entry, the leading N-by-M part of this array must
C             contain the matrix A. The elements corresponding to the
C             zero P-by-MIN(P,M) lower trapezoidal/triangular part
C             (if P > 0) are not referenced.
C             On exit, the elements on and above the diagonal of this
C             array contain the MIN(N,M)-by-M upper trapezoidal matrix
C             R (R is upper triangular, if N >= M) of the QR
C             factorization, and the relevant elements below the
C             diagonal contain the trailing components (the vectors v,
C             see Method) of the elementary reflectors used in the
C             factorization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) COMPLEX*16 array, dimension (LDB,L)
C             On entry, the leading N-by-L part of this array must
C             contain the matrix B.
C             On exit, the leading N-by-L part of this array contains
C             the updated matrix B.
C             If L = 0, this array is not referenced.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,N) if L > 0;
C             LDB >= 1        if L = 0.
C
C     TAU     (output) COMPLEX*16 array, dimension MIN(N,M)
C             The scalar factors of the elementary reflectors used.
C
C     Workspace
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  The length of the array ZWORK.
C             LZWORK >= MAX(1,M-1,M-P,L).
C             For optimum performance LZWORK should be larger.
C
C             If LZWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             ZWORK array, returns this value as the first entry of
C             the ZWORK array, and no error message related to LZWORK
C             is issued by XERBLA.
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
C     The routine uses min(N,M) Householder transformations exploiting
C     the zero pattern of the matrix.  A Householder matrix has the form
C
C                                     ( 1 ),
C        H  = I - tau *u *u',    u  = ( v )
C         i          i  i  i      i   (  i)
C
C     where v  is an (N-P+I-2)-vector.  The components of v  are stored
C            i                                             i
C     in the i-th column of A, beginning from the location i+1, and
C     tau  is stored in TAU(i).
C        i
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTORS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009.
C
C     KEYWORDS
C
C     Elementary reflector, QR factorization, unitary transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      PARAMETER         ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                    ONE  = ( 1.0D+0, 0.0D+0 ) )
C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, LDB, LZWORK, M, N, P
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), TAU(*), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           LQUERY
      INTEGER           I, NB, WRKOPT
      COMPLEX*16        FIRST
C     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV
C     .. External Subroutines ..
      EXTERNAL          XERBLA, ZGEQRF, ZLARF, ZLARFG, ZUNMQR
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, INT, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO = 0
      LQUERY = ( LZWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( L.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.1 .OR. ( L.GT.0 .AND. LDB.LT.N ) ) THEN
         INFO = -8
      ELSE
         I = MAX( 1, M - 1, M - P, L )
         IF( LQUERY ) THEN
            IF( M.GT.P ) THEN
               NB = ILAENV( 1, 'ZGEQRF', ' ', N-P, M-P, -1, -1 )
               WRKOPT = MAX( I, ( M - P )*NB )
               IF ( L.GT.0 ) THEN
                  NB = MIN( 64, ILAENV( 1, 'ZUNMQR', 'LC', N-P, L,
     $                                  MIN(N,M)-P, -1 ) )
                  WRKOPT = MAX( WRKOPT, MAX( 1, L )*NB )
               END IF
            END IF
         ELSE IF( LZWORK.LT.I ) THEN
            INFO = -11
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB04IZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         ZWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( M, N ).EQ.0 ) THEN
         ZWORK(1) = ONE
         RETURN
      ELSE IF( N.LE.P+1 ) THEN
         DO 5 I = 1, MIN( N, M )
            TAU(I) = ZERO
    5    CONTINUE
         ZWORK(1) = ONE
         RETURN
      END IF
C
C     Annihilate the subdiagonal elements of A and apply the
C     transformations to B, if L > 0.
C     Workspace: need MAX(M-1,L).
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of complex workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      DO 10 I = 1, MIN( P, M )
C
C        Exploit the structure of the I-th column of A.
C
         CALL ZLARFG( N-P, A(I,I), A(I+1,I), 1, TAU(I) )
         IF( TAU(I).NE.ZERO ) THEN
C
            FIRST = A(I,I)
            A(I,I) = ONE
C
            IF ( I.LT.M ) CALL ZLARF( 'Left', N-P, M-I, A(I,I), 1,
     $                                DCONJG( TAU(I) ), A(I,I+1), LDA,
     $                                ZWORK )
            IF ( L.GT.0 ) CALL ZLARF( 'Left', N-P, L, A(I,I), 1,
     $                                DCONJG( TAU(I) ), B(I,1), LDB,
     $                                ZWORK )
C
            A(I,I) = FIRST
         END IF
   10 CONTINUE
C
      WRKOPT = MAX( 1, M - 1, L )
C
C     Fast QR factorization of the remaining right submatrix, if any.
C     Workspace: need M-P;  prefer (M-P)*NB.
C
      IF( M.GT.P ) THEN
         CALL ZGEQRF( N-P, M-P, A(P+1,P+1), LDA, TAU(P+1), ZWORK,
     $                LZWORK, INFO )
         WRKOPT = MAX( WRKOPT, INT( ZWORK(1) ) )
C
         IF ( L.GT.0 ) THEN
C
C           Apply the transformations to B.
C           Workspace: need L;  prefer L*NB.
C
            CALL ZUNMQR( 'Left', 'Conjugate', N-P, L, MIN(N,M)-P,
     $                   A(P+1,P+1), LDA, TAU(P+1), B(P+1,1), LDB,
     $                   ZWORK, LZWORK, INFO )
            WRKOPT = MAX( WRKOPT, INT( ZWORK(1) ) )
         END IF
      END IF
C
      ZWORK(1) = WRKOPT
      RETURN
C *** Last line of MB04IZ ***
      END
