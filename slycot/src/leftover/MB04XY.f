      SUBROUTINE MB04XY( JOBU, JOBV, M, N, X, LDX, TAUP, TAUQ, U,
     $                   LDU, V, LDV, INUL, INFO )
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
C     To apply the Householder transformations Pj stored in factored
C     form into the columns of the array X, to the desired columns of
C     the matrix U by premultiplication, and/or the Householder
C     transformations Qj stored in factored form into the rows of the
C     array X, to the desired columns of the matrix V by
C     premultiplication. The Householder transformations Pj and Qj
C     are stored as produced by LAPACK Library routine DGEBRD.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBU    CHARACTER*1
C             Specifies whether to transform the columns in U as
C             follows:
C             = 'N':  Do not transform the columns in U;
C             = 'A':  Transform the columns in U (U has M columns);
C             = 'S':  Transform the columns in U (U has min(M,N)
C                     columns).
C
C     JOBV    CHARACTER*1
C             Specifies whether to transform the columns in V as
C             follows:
C             = 'N':  Do not transform the columns in V;
C             = 'A':  Transform the columns in V (V has N columns);
C             = 'S':  Transform the columns in V (V has min(M,N)
C                     columns).
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix X.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix X.  N >= 0.
C
C     X       (input) DOUBLE PRECISION array, dimension (LDX,N)
C             The leading M-by-N part contains in the columns of its
C             lower triangle the Householder transformations Pj, and
C             in the rows of its upper triangle the Householder
C             transformations Qj in factored form.
C             X is modified by the routine but restored on exit.
C
C     LDX     INTEGER
C             The leading dimension of the array X.   LDX >= MAX(1,M).
C
C     TAUP    (input) DOUBLE PRECISION array, dimension (MIN(M,N))
C             The scalar factors of the Householder transformations Pj.
C
C     TAUQ    (input) DOUBLE PRECISION array, dimension (MIN(M,N))
C             The scalar factors of the Householder transformations Qj.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*)
C             On entry, U contains the M-by-M (if JOBU = 'A') or
C             M-by-min(M,N) (if JOBU = 'S') matrix U.
C             On exit, the Householder transformations Pj have been
C             applied to each column i of U corresponding to a parameter
C             INUL(i) = .TRUE.
C             NOTE that U is not referenced if JOBU = 'N'.
C
C     LDU     INTEGER
C             The leading dimension of the array U.
C             LDU >= MAX(1,M), if JOBU = 'A' or JOBU = 'S';
C             LDU >= 1,        if JOBU = 'N'.
C
C     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*)
C             On entry, V contains the N-by-N (if JOBV = 'A') or
C             N-by-min(M,N) (if JOBV = 'S') matrix V.
C             On exit, the Householder transformations Qj have been
C             applied to each column i of V corresponding to a parameter
C             INUL(i) = .TRUE.
C             NOTE that V is not referenced if JOBV = 'N'.
C
C     LDV     INTEGER
C             The leading dimension of the array V.
C             LDV >= MAX(1,M), if JOBV = 'A' or JOBV = 'S';
C             LDV >= 1,        if JOBV = 'N'.
C
C     INUL    (input) LOGICAL array, dimension (MAX(M,N))
C             INUL(i) = .TRUE. if the i-th column of U and/or V is to be
C             transformed, and INUL(i) = .FALSE., otherwise.
C             (1 <= i <= MAX(M,N)).
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
C     The Householder transformations Pj or Qj are applied to the
C     columns of U or V indexed by I for which INUL(I) = .TRUE..
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB04PZ by S. Van Huffel, Katholieke
C     University Leuven, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Bidiagonalization, orthogonal transformation, singular subspace,
C     singular value decomposition.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBU, JOBV
      INTEGER           INFO, LDU, LDV, LDX, M, N
C     .. Array Arguments ..
      LOGICAL           INUL(*)
      DOUBLE PRECISION  TAUP(*), TAUQ(*), U(LDU,*), V(LDV,*),
     $                  X(LDX,*)
C     .. Local Scalars ..
      LOGICAL           LJOBUA, LJOBUS, LJOBVA, LJOBVS, WANTU, WANTV
      INTEGER           I, IM, IOFF, L, NCOL, P
      DOUBLE PRECISION  FIRST
C     .. Local Arrays ..
      DOUBLE PRECISION  DWORK(1)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLARF, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MAX
C     .. Executable Statements ..
C
      INFO = 0
      LJOBUA = LSAME( JOBU, 'A' )
      LJOBUS = LSAME( JOBU, 'S' )
      LJOBVA = LSAME( JOBV, 'A' )
      LJOBVS = LSAME( JOBV, 'S' )
      WANTU  = LJOBUA.OR.LJOBUS
      WANTV  = LJOBVA.OR.LJOBVS
C
C     Test the input scalar arguments.
C
      IF( .NOT.WANTU .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.WANTV .AND. .NOT.LSAME( JOBV, 'N' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDX.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( ( WANTU.AND.LDU.LT.MAX( 1, M ) ) .OR.
     $    ( .NOT.WANTU.AND.LDU.LT.1 ) ) THEN
         INFO = -10
      ELSE IF( ( WANTV.AND.LDV.LT.MAX( 1, N ) ) .OR.
     $    ( .NOT.WANTV.AND.LDV.LT.1 ) ) THEN
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return
C
         CALL XERBLA( 'MB04XY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      P = MIN( M, N )
      IF ( P.EQ.0 )
     $   RETURN
C
      IF ( M.LT.N ) THEN
         IOFF = 1
      ELSE
         IOFF = 0
      END IF
C
C     Apply the Householder transformations Pj onto the desired
C     columns of U.
C
      IM = MIN( M-1, N )
      IF ( WANTU .AND. ( IM.GT.0 ) ) THEN
         IF ( LJOBUA ) THEN
            NCOL = M
         ELSE
            NCOL = P
         END IF
C
         DO 40 I = 1, NCOL
            IF ( INUL(I) ) THEN
C
               DO 20 L = IM, 1, -1
                  IF ( TAUP(L).NE.ZERO ) THEN
                     FIRST = X(L+IOFF,L)
                     X(L+IOFF,L) = ONE
                     CALL DLARF( 'Left', M-L+1-IOFF, 1, X(L+IOFF,L), 1,
     $                           TAUP(L), U(L+IOFF,I), LDU, DWORK )
                     X(L+IOFF,L) = FIRST
                  END IF
   20          CONTINUE
C
            END IF
   40    CONTINUE
C
      END IF
C
C     Apply the Householder transformations Qj onto the desired columns
C     of V.
C
      IM = MIN( N-1, M )
      IF ( WANTV .AND. ( IM.GT.0 ) ) THEN
         IF ( LJOBVA ) THEN
            NCOL = N
         ELSE
            NCOL = P
         END IF
C
         DO 80 I = 1, NCOL
            IF ( INUL(I) ) THEN
C
               DO 60 L = IM, 1, -1
                  IF ( TAUQ(L).NE.ZERO ) THEN
                     FIRST = X(L,L+1-IOFF)
                     X(L,L+1-IOFF) = ONE
                     CALL DLARF( 'Left', N-L+IOFF, 1, X(L,L+1-IOFF),
     $                           LDX, TAUQ(L), V(L+1-IOFF,I), LDV,
     $                           DWORK )
                     X(L,L+1-IOFF) = FIRST
                  END IF
   60          CONTINUE
C
            END IF
   80    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of MB04XY ***
      END
