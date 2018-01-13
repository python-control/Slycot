      SUBROUTINE SB06ND( N, M, KMAX, A, LDA, B, LDB, KSTAIR, U, LDU, F,
     $                   LDF, DWORK, INFO )
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
C     To construct the minimum norm feedback matrix F to perform
C     "deadbeat control" on a (A,B)-pair of a state-space model (which
C     must be preliminarily reduced to upper "staircase" form using
C     SLICOT Library routine AB01OD) such that the matrix R = A + BFU'
C     is nilpotent.
C     (The transformation matrix U reduces R to upper Schur form with
C     zero blocks on its diagonal (of dimension KSTAIR(i)) and
C     therefore contains bases for the i-th controllable subspaces,
C     where i = 1,...,KMAX).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The actual state dimension, i.e. the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The actual input dimension.  M >= 0.
C
C     KMAX    (input) INTEGER
C             The number of "stairs" in the staircase form as produced
C             by SLICOT Library routine AB01OD.  0 <= KMAX <= N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the transformed state-space matrix of the
C             (A,B)-pair with triangular stairs, as produced by SLICOT
C             Library routine AB01OD (with option STAGES = 'A').
C             On exit, the leading N-by-N part of this array contains
C             the matrix U'AU + U'BF.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the transformed triangular input matrix of the
C             (A,B)-pair as produced by SLICOT Library routine AB01OD
C             (with option STAGES = 'A').
C             On exit, the leading N-by-M part of this array contains
C             the matrix U'B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     KSTAIR  (input) INTEGER array, dimension (KMAX)
C             The leading KMAX elements of this array must contain the
C             dimensions of each "stair" as produced by SLICOT Library
C             routine AB01OD.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             On entry, the leading N-by-N part of this array must
C             contain either a transformation matrix (e.g. from a
C             previous call to other SLICOT routine) or be initialised
C             as the identity matrix.
C             On exit, the leading N-by-N part of this array contains
C             the product of the input matrix U and the state-space
C             transformation matrix which reduces A + BFU' to real
C             Schur form.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,N).
C
C     F       (output) DOUBLE PRECISION array, dimension (LDF,N)
C             The leading M-by-N part of this array contains the
C             deadbeat feedback matrix F.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (2*N)
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
C     Starting from the (A,B)-pair in "staircase form" with "triangular"
C     stairs, dimensions KSTAIR(i+1) x KSTAIR(i), (described by the
C     vector KSTAIR):
C
C                    | B | A      *  . . .  *  |
C                    |  1|  11       .      .  |
C                    |   | A     A     .    .  |
C                    |   |  21    22     .  .  |
C                    |   |    .      .     .   |
C      [ B | A ]  =  |   |      .      .    *  |
C                    |   |        .      .     |
C                    | 0 |   0                 |
C                    |   |          A      A   |
C                    |   |           r,r-1  rr |
C
C     where the i-th diagonal block of A has dimension KSTAIR(i), for
C     i = 1,2,...,r, the feedback matrix F is constructed recursively in
C     r steps (where the number of "stairs" r is given by KMAX). In each
C     step a unitary state-space transformation U and a part of F are
C     updated in order to achieve the final form:
C
C                       | 0   A      *   . . .  *   |
C                       |      12      .        .   |
C                       |                .      .   |
C                       |     0    A       .    .   |
C                       |           23       .  .   |
C                       |         .      .          |
C     [ U'AU + U'BF ] = |           .      .    *   | .
C                       |             .      .      |
C                       |                           |
C                       |                     A     |
C                       |                      r-1,r|
C                       |                           |
C                       |                       0   |
C
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         Deadbeat control: a special inverse eigenvalue problem.
C         BIT, 24, pp. 681-699, 1984.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O((N + M) * N**2) operations and is mixed
C     numerical stable (see [1]).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997.
C     Supersedes Release 2.0 routine SB06BD by M. Vanbegin, and
C     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium.
C
C     REVISIONS
C
C     1997, December 10; 2003, September 27.
C
C     KEYWORDS
C
C     Canonical form, deadbeat control, eigenvalue assignment, feedback
C     control, orthogonal transformation, real Schur form, staircase
C     form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, KMAX, LDA, LDB, LDF, LDU, M, N
C     .. Array Arguments ..
      INTEGER           KSTAIR(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DWORK(*), F(LDF,*), U(LDU,*)
C     .. Local Scalars ..
      INTEGER           J, J0, JCUR, JKCUR, JMKCUR, KCUR, KK, KMIN,
     $                  KSTEP, MKCUR, NCONT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DLACPY, DLARFG, DLASET,
     $                  SLCT_DLATZM, DTRSM, XERBLA
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( KMAX.LT.0 .OR. KMAX.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDU.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -12
      ELSE
         NCONT = 0
C
         DO 10 KK = 1, KMAX
            NCONT = NCONT + KSTAIR(KK)
   10    CONTINUE
C
         IF( NCONT.GT.N )
     $      INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB06ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
C
      DO 120 KMIN = 1, KMAX
         JCUR = NCONT
         KSTEP = KMAX - KMIN
C
C        Triangularize bottom part of A (if KSTEP > 0).
C
         DO 40 KK = KMAX, KMAX - KSTEP + 1, -1
            KCUR = KSTAIR(KK)
C
C           Construct Ukk and store in Fkk.
C
            DO 20 J = 1, KCUR
               JMKCUR = JCUR - KCUR
               CALL DCOPY( KCUR, A(JCUR,JMKCUR), LDA, F(1,JCUR), 1 )
               CALL DLARFG( KCUR+1, A(JCUR,JCUR), F(1,JCUR), 1,
     $                      DWORK(JCUR) )
               CALL DLASET( 'Full', 1, KCUR, ZERO, ZERO, A(JCUR,JMKCUR),
     $                      LDA )
C
C              Backmultiply A and U with Ukk.
C
               CALL SLCT_DLATZM( 'Right', JCUR-1, KCUR+1, F(1,JCUR), 1,
     $                      DWORK(JCUR), A(1,JCUR), A(1,JMKCUR), LDA,
     $                      DWORK )
C
               CALL SLCT_DLATZM( 'Right', N, KCUR+1, F(1,JCUR), 1,
     $                      DWORK(JCUR), U(1,JCUR), U(1,JMKCUR), LDU,
     $                      DWORK(N+1) )
               JCUR = JCUR - 1
   20       CONTINUE
C
   40    CONTINUE
C
C        Eliminate diagonal block Aii by feedback Fi.
C
         KCUR = KSTAIR(KMIN)
         J0 = JCUR - KCUR + 1
         MKCUR = M - KCUR + 1
C
C        Solve for Fi and add B x Fi to A.
C
         CALL DLACPY( 'Full', KCUR, KCUR, A(J0,J0), LDA, F(MKCUR,J0),
     $                LDF )
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', KCUR,
     $               KCUR, -ONE, B(J0,MKCUR), LDB, F(MKCUR,J0), LDF )
         IF ( J0.GT.1 )
     $      CALL DGEMM( 'No transpose', 'No transpose', J0-1, KCUR,
     $                  KCUR, ONE, B(1,MKCUR), LDB, F(MKCUR,J0), LDF,
     $                  ONE, A(1,J0), LDA )
         CALL DLASET( 'Full',   KCUR, KCUR, ZERO, ZERO, A(J0,J0), LDA )
         CALL DLASET( 'Full', M-KCUR, KCUR, ZERO, ZERO, F(1,J0), LDF )
C
         IF ( KSTEP.NE.0 ) THEN
            JKCUR = NCONT
C
C           Premultiply A with Ukk.
C
            DO 80 KK = KMAX, KMAX - KSTEP + 1, -1
               KCUR = KSTAIR(KK)
               JCUR = JKCUR - KCUR
C
               DO 60 J = 1, KCUR
                  CALL SLCT_DLATZM( 'Left', KCUR+1, N-JCUR+1,
     $                         F(1,JKCUR), 1,
     $                         DWORK(JKCUR), A(JKCUR,JCUR),
     $                         A(JCUR,JCUR), LDA, DWORK(N+1) )
                  JCUR = JCUR - 1
                  JKCUR = JKCUR - 1
   60          CONTINUE
C
   80       CONTINUE
C
C           Premultiply B with Ukk.
C
            JCUR  = JCUR + KCUR
            JKCUR = JCUR + KCUR
C
            DO 100 J = M, M - KCUR + 1, -1
               CALL SLCT_DLATZM( 'Left', KCUR+1, M-J+1, F(1,JKCUR), 1,
     $                      DWORK(JKCUR), B(JKCUR,J), B(JCUR,J), LDB,
     $                      DWORK(N+1) )
               JCUR = JCUR - 1
               JKCUR = JKCUR - 1
  100       CONTINUE
C
         END IF
  120 CONTINUE
C
      IF ( NCONT.NE.N )
     $   CALL DLASET( 'Full', M, N-NCONT, ZERO, ZERO, F(1,NCONT+1),
     $                LDF )
C
      RETURN
C *** Last line of SB06ND ***
      END
