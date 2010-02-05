      SUBROUTINE TG01AD( JOB, L, N, M, P, THRESH, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, LSCALE, RSCALE, DWORK, INFO )
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
C     To balance the matrices of the system pencil
C
C             S =  ( A  B ) - lambda ( E  0 ) :=  Q - lambda Z,
C                  ( C  0 )          ( 0  0 )
C
C     corresponding to the descriptor triple (A-lambda E,B,C),
C     by balancing. This involves diagonal similarity transformations
C     (Dl*A*Dr - lambda Dl*E*Dr, Dl*B, C*Dr) applied to the system
C     (A-lambda E,B,C) to make the rows and columns of system pencil
C     matrices
C
C                  diag(Dl,I) * S * diag(Dr,I)
C
C     as close in norm as possible. Balancing may reduce the 1-norms
C     of the matrices of the system pencil S.
C
C     The balancing can be performed optionally on the following
C     particular system pencils
C
C              S = A-lambda E,
C
C              S = ( A-lambda E  B ),    or
C
C              S = ( A-lambda E ).
C                  (     C      )
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates which matrices are involved in balancing, as
C             follows:
C             = 'A':  All matrices are involved in balancing;
C             = 'B':  B, A and E matrices are involved in balancing;
C             = 'C':  C, A and E matrices are involved in balancing;
C             = 'N':  B and C matrices are not involved in balancing.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of rows of matrices A, B, and E.  L >= 0.
C
C     N       (input) INTEGER
C             The number of columns of matrices A, E, and C.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of matrix C.  P >= 0.
C
C     THRESH  (input) DOUBLE PRECISION
C             Threshold value for magnitude of elements:
C             elements with magnitude less than or equal to
C             THRESH are ignored for balancing. THRESH >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading L-by-N part of this array contains
C             the balanced matrix Dl*A*Dr.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the descriptor matrix E.
C             On exit, the leading L-by-N part of this array contains
C             the balanced matrix Dl*E*Dr.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,L).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading L-by-M part of this array must
C             contain the input/state matrix B.
C             On exit, if M > 0, the leading L-by-M part of this array
C             contains the balanced matrix Dl*B.
C             The array B is not referenced if M = 0.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, if P > 0, the leading P-by-N part of this array
C             contains the balanced matrix C*Dr.
C             The array C is not referenced if P = 0.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     LSCALE  (output) DOUBLE PRECISION array, dimension (L)
C             The scaling factors applied to S from left.  If Dl(j) is
C             the scaling factor applied to row j, then
C             SCALE(j) = Dl(j), for j = 1,...,L.
C
C     RSCALE  (output) DOUBLE PRECISION array, dimension (N)
C             The scaling factors applied to S from right.  If Dr(j) is
C             the scaling factor applied to column j, then
C             SCALE(j) = Dr(j), for j = 1,...,N.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (3*(L+N))
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
C     Balancing consists of applying a diagonal similarity
C     transformation
C                            -1
C                  diag(Dl,I)  * S * diag(Dr,I)
C
C     to make the 1-norms of each row of the first L rows of S and its
C     corresponding N columns nearly equal.
C
C     Information about the diagonal matrices Dl and Dr are returned in
C     the vectors LSCALE and RSCALE, respectively.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     [2] R.C. Ward, R. C.
C         Balancing the generalized eigenvalue problem.
C         SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the LAPACK routine DGGBAL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 1999,
C     May 2003, March 2004, Jan. 2009.
C
C     KEYWORDS
C
C     Balancing, eigenvalue, matrix algebra, matrix operations,
C     similarity transformation.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   HALF, ONE, ZERO
      PARAMETER          ( HALF = 0.5D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      DOUBLE PRECISION   SCLFAC, THREE
      PARAMETER          ( SCLFAC = 1.0D+1, THREE = 3.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, L, LDA, LDB, LDC, LDE, M, N, P
      DOUBLE PRECISION   THRESH
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK( * ), E( LDE, * ), LSCALE( * ),
     $                   RSCALE( * )
C     .. Local Scalars ..
      LOGICAL            WITHB, WITHC
      INTEGER            I, ICAB, IR, IRAB, IT, J, JC, KOUNT, KW1, KW2,
     $                   KW3, KW4, KW5, LCAB, LRAB, LSFMAX, LSFMIN,
     $                   NRP2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2,
     $                   COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX,
     $                   SFMIN, SUM, T, TA, TB, TC, TE
C     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           DDOT, DLAMCH, IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG10, MAX, MIN, SIGN
C
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO  = 0
      WITHB = LSAME( JOB, 'A' ) .OR. LSAME( JOB, 'B' )
      WITHC = LSAME( JOB, 'A' ) .OR. LSAME( JOB, 'C' )
C
      IF( .NOT.WITHB .AND. .NOT.WITHC .AND. .NOT.LSAME( JOB, 'N' ) )
     $   THEN
         INFO = -1
      ELSE IF( L.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -8
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.L ) ) THEN
         INFO = -12
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( L.EQ.0 .OR. N.EQ.0 ) THEN
         DUM( 1 ) = ONE
         IF( L.GT.0 ) THEN
            CALL DCOPY( L, DUM, 0, LSCALE, 1 )
         ELSE IF( N.GT.0 ) THEN
            CALL DCOPY( N, DUM, 0, RSCALE, 1 )
         END IF
         RETURN
      END IF
C
C     Initialize balancing and allocate work storage.
C
      KW1 = N
      KW2 = KW1 + L
      KW3 = KW2 + L
      KW4 = KW3 + N
      KW5 = KW4 + L
      DUM( 1 ) = ZERO
      CALL DCOPY( L, DUM, 0, LSCALE, 1 )
      CALL DCOPY( N, DUM, 0, RSCALE, 1 )
      CALL DCOPY( 3*(L+N), DUM, 0, DWORK, 1 )
C
C     Compute right side vector in resulting linear equations.
C
      BASL = LOG10( SCLFAC )
      DO 20 I = 1, L
         DO 10 J = 1, N
            TE = ABS( E( I, J ) )
            TA = ABS( A( I, J ) )
            IF( TA.GT.THRESH ) THEN
               TA = LOG10( TA ) / BASL
            ELSE
               TA = ZERO
            END IF
            IF( TE.GT.THRESH ) THEN
               TE = LOG10( TE ) / BASL
            ELSE
               TE = ZERO
            END IF
            DWORK( I+KW4 ) = DWORK( I+KW4 ) - TA - TE
            DWORK( J+KW5 ) = DWORK( J+KW5 ) - TA - TE
  10     CONTINUE
  20  CONTINUE
C
      IF( M.EQ.0 ) THEN
         WITHB = .FALSE.
         TB    = ZERO
      END IF
      IF( P.EQ.0 ) THEN
         WITHC = .FALSE.
         TC    = ZERO
      END IF
C
      IF( WITHB ) THEN
         DO 30 I = 1, L
            J = IDAMAX( M, B( I, 1 ), LDB )
            TB = ABS( B( I, J ) )
            IF( TB.GT.THRESH ) THEN
               TB = LOG10( TB ) / BASL
               DWORK( I+KW4 ) = DWORK( I+KW4 ) - TB
            END IF
  30     CONTINUE
      END IF
C
      IF( WITHC ) THEN
         DO 40 J = 1, N
            I = IDAMAX( P, C( 1, J ), 1 )
            TC = ABS( C( I, J ) )
            IF( TC.GT.THRESH ) THEN
               TC = LOG10( TC ) / BASL
               DWORK( J+KW5 ) = DWORK( J+KW5 ) - TC
            END IF
  40     CONTINUE
      END IF
C
      COEF  = ONE / DBLE( L+N )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = MAX( L, N ) + 2
      BETA = ZERO
      IT = 1
C
C     Start generalized conjugate gradient iteration.
C
  50  CONTINUE
C
      GAMMA = DDOT( L, DWORK( 1+KW4 ), 1, DWORK( 1+KW4 ), 1 ) +
     $        DDOT( N, DWORK( 1+KW5 ), 1, DWORK( 1+KW5 ), 1 )
C
      EW = ZERO
      DO 60 I = 1, L
         EW = EW + DWORK( I+KW4 )
  60  CONTINUE
C
      EWC = ZERO
      DO 70 I = 1, N
         EWC = EWC + DWORK( I+KW5 )
  70  CONTINUE
C
      GAMMA = COEF*GAMMA - COEF2*( EW**2 + EWC**2 ) -
     $                     COEF5*( EW - EWC )**2
      IF( GAMMA.EQ.ZERO )
     $   GO TO 160
      IF( IT.NE.1 )
     $   BETA = GAMMA / PGAMMA
      T  = COEF5*( EWC - THREE*EW  )
      TC = COEF5*( EW  - THREE*EWC )
C
      CALL DSCAL( N+L, BETA, DWORK, 1 )
C
      CALL DAXPY( L, COEF, DWORK( 1+KW4 ), 1, DWORK( 1+KW1 ), 1 )
      CALL DAXPY( N, COEF, DWORK( 1+KW5 ), 1, DWORK, 1 )
C
      DO 80 J = 1, N
         DWORK( J ) = DWORK( J ) + TC
   80 CONTINUE
C
      DO 90 I = 1, L
         DWORK( I+KW1 ) = DWORK( I+KW1 ) + T
   90 CONTINUE
C
C     Apply matrix to vector.
C
      DO 110 I = 1, L
         KOUNT = 0
         SUM = ZERO
         DO 100 J = 1, N
            IF( ABS( A( I, J ) ).GT.THRESH ) THEN
               KOUNT = KOUNT + 1
               SUM = SUM + DWORK( J )
            END IF
            IF( ABS( E( I, J ) ).GT.THRESH ) THEN
               KOUNT = KOUNT + 1
               SUM = SUM + DWORK( J )
            END IF
  100    CONTINUE
         IF( WITHB ) THEN
            J = IDAMAX( M, B( I, 1 ), LDB )
            IF( ABS( B( I, J ) ).GT.THRESH ) KOUNT = KOUNT + 1
         END IF
         DWORK( I+KW2 ) = DBLE( KOUNT )*DWORK( I+KW1 ) + SUM
  110 CONTINUE
C
      DO 130 J = 1, N
         KOUNT = 0
         SUM = ZERO
         DO 120 I = 1, L
            IF( ABS( A( I, J ) ).GT.THRESH ) THEN
               KOUNT = KOUNT + 1
               SUM = SUM + DWORK( I+KW1 )
            END IF
            IF( ABS( E( I, J ) ).GT.THRESH ) THEN
               KOUNT = KOUNT + 1
               SUM = SUM + DWORK( I+KW1 )
            END IF
  120    CONTINUE
         IF( WITHC ) THEN
            I = IDAMAX( P, C( 1, J ), 1 )
            IF( ABS( C( I, J ) ).GT.THRESH ) KOUNT = KOUNT + 1
         END IF
         DWORK( J+KW3 ) = DBLE( KOUNT )*DWORK( J ) + SUM
  130 CONTINUE
C
      SUM = DDOT( L, DWORK( 1+KW1 ), 1, DWORK( 1+KW2 ), 1 ) +
     $      DDOT( N, DWORK, 1, DWORK( 1+KW3 ), 1 )
      ALPHA = GAMMA / SUM
C
C     Determine correction to current iteration.
C
      CMAX = ZERO
      DO 140 I = 1, L
         COR = ALPHA*DWORK( I+KW1 )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
  140 CONTINUE
C
      DO 150 J = 1, N
         COR = ALPHA*DWORK( J )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         RSCALE( J ) = RSCALE( J ) + COR
  150 CONTINUE
      IF( CMAX.LT.HALF )
     $   GO TO 160
C
      CALL DAXPY( L, -ALPHA, DWORK( 1+KW2 ), 1, DWORK( 1+KW4 ), 1 )
      CALL DAXPY( N, -ALPHA, DWORK( 1+KW3 ), 1, DWORK( 1+KW5 ), 1 )
C
      PGAMMA = GAMMA
      IT = IT + 1
      IF( IT.LE.NRP2 )
     $   GO TO 50
C
C     End generalized conjugate gradient iteration.
C
  160 CONTINUE
      SFMIN = DLAMCH( 'Safe minimum' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL + ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
C
C     Compute left diagonal scaling matrix.
C
      DO 170 I = 1, L
         IRAB = IDAMAX( N, A( I, 1 ), LDA )
         RAB  = ABS( A( I, IRAB ) )
         IRAB = IDAMAX( N, E( I, 1 ), LDE )
         RAB  = MAX( RAB, ABS( E( I, IRAB ) ) )
         IF( WITHB ) THEN
            IRAB = IDAMAX( M, B( I, 1 ), LDB )
            RAB  = MAX( RAB, ABS( B( I, IRAB ) ) )
         END IF
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL + ONE )
         IR = LSCALE( I ) + SIGN( HALF, LSCALE( I ) )
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
  170 CONTINUE
C
C     Compute right diagonal scaling matrix.
C
      DO 180 J = 1, N
         ICAB = IDAMAX( L, A( 1, J ), 1 )
         CAB  = ABS( A( ICAB, J ) )
         ICAB = IDAMAX( L, E( 1, J ), 1 )
         CAB  = MAX( CAB, ABS( E( ICAB, J ) ) )
         IF( WITHC ) THEN
            ICAB = IDAMAX( P, C( 1, J ), 1 )
            CAB  = MAX( CAB, ABS( C( ICAB, J ) ) )
         END IF
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL + ONE )
         JC = RSCALE( J ) + SIGN( HALF, RSCALE( J ) )
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( J ) = SCLFAC**JC
  180 CONTINUE
C
C     Row scaling of matrices A, E and B.
C
      DO 190 I = 1, L
         CALL DSCAL( N, LSCALE( I ), A( I, 1 ), LDA )
         CALL DSCAL( N, LSCALE( I ), E( I, 1 ), LDE )
         IF( WITHB )
     $      CALL DSCAL( M, LSCALE( I ), B( I, 1 ), LDB )
  190 CONTINUE
C
C     Column scaling of matrices A, E and C.
C
      DO 200 J = 1, N
         CALL DSCAL( L, RSCALE( J ), A( 1, J ), 1 )
         CALL DSCAL( L, RSCALE( J ), E( 1, J ), 1 )
         IF( WITHC )
     $      CALL DSCAL( P, RSCALE( J ), C( 1, J ), 1 )
  200 CONTINUE
C
      RETURN
C *** Last line of TG01AD ***
      END
