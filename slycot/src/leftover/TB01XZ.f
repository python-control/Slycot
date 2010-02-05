      SUBROUTINE TB01XZ( JOBD, N, M, P, KL, KU, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, INFO )
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
C     To apply a special transformation to a system given as a triple
C     (A,B,C),
C
C        A <-- P * A' * P,  B <-- P * C',  C <-- B' * P,
C
C     where P is a matrix with 1 on the secondary diagonal, and with 0
C     in the other entries. Matrix A can be specified as a band matrix.
C     Optionally, matrix D of the system can be transposed. This
C     transformation is actually a special similarity transformation of
C     the dual system.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears in
C             the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A, the number of rows of matrix B
C             and the number of columns of matrix C.
C             N represents the dimension of the state vector.  N >= 0.
C
C     M       (input) INTEGER.
C             The number of columns of matrix B.
C             M represents the dimension of input vector.  M >= 0.
C
C     P       (input) INTEGER.
C             The number of rows of matrix C.
C             P represents the dimension of output vector.  P >= 0.
C
C     KL      (input) INTEGER
C             The number of subdiagonals of A to be transformed.
C             MAX( 0, N-1 ) >= KL >= 0.
C
C     KU      (input) INTEGER
C             The number of superdiagonals of A to be transformed.
C             MAX( 0, N-1 ) >= KU >= 0.
C
C     A       (input/output) COMPLEX*16 array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the system state matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the transformed (pertransposed) matrix P*A'*P.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) COMPLEX*16 array, dimension (LDB,MAX(M,P))
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, the leading N-by-P part of this array contains
C             the dual input/state matrix P*C'.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,N) if M > 0 or  P > 0.
C             LDB >= 1        if M = 0 and P = 0.
C
C     C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, the leading M-by-N part of this array contains
C             the dual state/output matrix B'*P.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,M,P) if N > 0.
C             LDC >= 1          if N = 0.
C
C     D       (input/output) COMPLEX*16 array, dimension (LDD,MAX(M,P))
C             On entry, if JOBD = 'D', the leading P-by-M part of this
C             array must contain the original direct transmission
C             matrix D.
C             On exit, if JOBD = 'D', the leading M-by-P part of this
C             array contains the transposed direct transmission matrix
C             D'. The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,M,P) if JOBD = 'D'.
C             LDD >= 1          if JOBD = 'Z'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit.
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The rows and/or columns of the matrices of the triplet (A,B,C)
C     and, optionally, of the matrix D are swapped in a special way.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Nov. 2008.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Matrix algebra, matrix operations, similarity transformation.
C
C  *********************************************************************
C
C     ..
C     .. Scalar Arguments ..
      CHARACTER          JOBD
      INTEGER            INFO, KL, KU, LDA, LDB, LDC, LDD, M, N, P
C     ..
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * )
C     ..
C     .. Local Scalars ..
      LOGICAL            LJOBD
      INTEGER            J, J1, LDA1, MAXMP, MINMP, NM1
C     ..
C     .. External functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input arguments.
C
      INFO  = 0
      LJOBD = LSAME( JOBD, 'D' )
      MAXMP = MAX( M, P )
      MINMP = MIN( M, P )
      NM1   = N - 1
C
      IF( .NOT.LJOBD .AND. .NOT.LSAME( JOBD, 'Z' )  ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( KL.LT.0 .OR. KL.GT.MAX( 0, NM1 ) ) THEN
         INFO = -5
      ELSE IF( KU.LT.0 .OR. KU.GT.MAX( 0, NM1 ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( ( MAXMP.GT.0 .AND. LDB.LT.MAX( 1, N ) ) .OR.
     $         ( MINMP.EQ.0 .AND. LDB.LT.1 ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.1 .OR. ( N.GT.0 .AND. LDC.LT.MAXMP ) ) THEN
         INFO = -12
      ELSE IF( LDD.LT.1 .OR. ( LJOBD  .AND. LDD.LT.MAXMP ) ) THEN
         INFO = -14
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01XZ', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( LJOBD ) THEN
C
C        Replace D by D', if non-scalar.
C
         DO 5 J = 1, MAXMP
            IF ( J.LT.MINMP ) THEN
               CALL ZSWAP( MINMP-J, D(J+1,J), 1, D(J,J+1), LDD )
            ELSE IF ( J.GT.P ) THEN
               CALL ZCOPY( P, D(1,J), 1, D(J,1), LDD )
            ELSE IF ( J.GT.M ) THEN
               CALL ZCOPY( M, D(J,1), LDD, D(1,J), 1 )
            END IF
    5    CONTINUE
C
      END IF
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Replace matrix A by P*A'*P.
C
      IF ( KL.EQ.NM1 .AND. KU.EQ.NM1 ) THEN
C
C        Full matrix A.
C
         DO 10 J = 1, NM1
            CALL ZSWAP( N-J, A( 1, J ), 1, A( N-J+1, J+1 ), -LDA )
   10    CONTINUE
C
      ELSE
C
C        Band matrix A.
C
         LDA1 = LDA + 1
C
C        Pertranspose the KL subdiagonals.
C
         DO 20 J = 1, MIN( KL, N-2 )
            J1 = ( N - J )/2
            CALL ZSWAP( J1, A(J+1,1), LDA1, A(N-J1+1,N-J1+1-J), -LDA1 )
   20    CONTINUE
C
C        Pertranspose the KU superdiagonals.
C
         DO 30 J = 1, MIN( KU, N-2 )
            J1 = ( N - J )/2
            CALL ZSWAP( J1, A(1,J+1), LDA1, A(N-J1+1-J,N-J1+1), -LDA1 )
   30    CONTINUE
C
C        Pertranspose the diagonal.
C
         J1 = N/2
         CALL ZSWAP( J1, A(1,1), LDA1, A(N-J1+1,N-J1+1), -LDA1 )
C
      END IF
C
C     Replace matrix B by P*C' and matrix C by B'*P.
C
      DO 40 J = 1, MAXMP
         IF ( J.LE.MINMP ) THEN
            CALL ZSWAP( N, B(1,J), 1, C(J,1), -LDC )
         ELSE IF ( J.GT.P ) THEN
            CALL ZCOPY( N, B(1,J), 1, C(J,1), -LDC )
         ELSE
            CALL ZCOPY( N, C(J,1), -LDC, B(1,J), 1 )
         END IF
   40 CONTINUE
C
      RETURN
C *** Last line of TB01XZ ***
      END
