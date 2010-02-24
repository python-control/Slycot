      SUBROUTINE AB07MD( JOBD, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
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
C     To find the dual of a given state-space representation.
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
C             The order of the state-space representation.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the dual state dynamics matrix A'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,MAX(M,P))
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, the leading N-by-P part of this array contains
C             the dual input/state matrix C'.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, the leading M-by-N part of this array contains
C             the dual state/output matrix B'.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,M,P) if N > 0.
C             LDC >= 1 if N = 0.
C
C     D       (input/output) DOUBLE PRECISION array, dimension
C             (LDD,MAX(M,P))
C             On entry, if JOBD = 'D', the leading P-by-M part of this
C             array must contain the original direct transmission
C             matrix D.
C             On exit, if JOBD = 'D', the leading M-by-P part of this
C             array contains the dual direct transmission matrix D'.
C             The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,M,P) if JOBD = 'D'.
C             LDD >= 1 if JOBD = 'Z'.
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
C     If the given state-space representation is the M-input/P-output
C     (A,B,C,D), its dual is simply the P-input/M-output (A',C',B',D').
C
C     REFERENCES
C
C     None
C
C     NUMERICAL ASPECTS
C
C     None
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine AB07AD by T.W.C.Williams, Kingston
C     Polytechnic, United Kingdom, March 1982.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
C
C     KEYWORDS
C
C     Dual system, state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER         JOBD
      INTEGER           INFO, LDA, LDB, LDC, LDD, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*)
C     .. Local Scalars ..
      LOGICAL           LJOBD
      INTEGER           J, MINMP, MPLIM
C     .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External subroutines ..
      EXTERNAL          DCOPY, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LJOBD = LSAME( JOBD, 'D' )
      MPLIM = MAX( M, P )
      MINMP = MIN( M, P )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBD .AND. .NOT.LSAME( JOBD, 'Z' )  ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( ( N.GT.0 .AND. LDC.LT.MAX( 1, MPLIM ) ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -10
      ELSE IF( ( LJOBD .AND. LDD.LT.MAX( 1, MPLIM ) ) .OR.
     $    ( .NOT.LJOBD .AND. LDD.LT.1 ) ) THEN
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB07MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, MINMP ).EQ.0 )
     $   RETURN
C
      IF ( N.GT.0 ) THEN
C
C        Transpose A, if non-scalar.
C
         DO 10 J = 1, N - 1
            CALL DSWAP( N-J, A(J+1,J), 1, A(J,J+1), LDA )
   10    CONTINUE
C
C        Replace B by C' and C by B'.
C
         DO 20 J = 1, MPLIM
            IF ( J.LE.MINMP ) THEN
               CALL DSWAP( N, B(1,J), 1, C(J,1), LDC )
            ELSE IF ( J.GT.P ) THEN
               CALL DCOPY( N, B(1,J), 1, C(J,1), LDC )
            ELSE
               CALL DCOPY( N, C(J,1), LDC, B(1,J), 1 )
            END IF
   20    CONTINUE
C
      END IF
C
      IF ( LJOBD .AND. MINMP.GT.0 ) THEN
C
C        Transpose D, if non-scalar.
C
         DO 30 J = 1, MPLIM
            IF ( J.LT.MINMP ) THEN
               CALL DSWAP( MINMP-J, D(J+1,J), 1, D(J,J+1), LDD )
            ELSE IF ( J.GT.P ) THEN
               CALL DCOPY( P, D(1,J), 1, D(J,1), LDD )
            ELSE IF ( J.GT.M ) THEN
               CALL DCOPY( M, D(J,1), LDD, D(1,J), 1 )
            END IF
   30    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of AB07MD ***
      END
