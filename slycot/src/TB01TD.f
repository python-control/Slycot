      SUBROUTINE TB01TD( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD, LOW,
     $                   IGH, SCSTAT, SCIN, SCOUT, DWORK, INFO )
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
C     To reduce a given state-space representation (A,B,C,D) to
C     balanced form by means of state permutations and state, input and
C     output scalings.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state-space representation, i.e. the
C             order of the original state dynamics matrix A.  N >= 0.
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
C             the balanced state dynamics matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the balanced input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the balanced state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the original direct transmission matrix D.
C             On exit, the leading P-by-M part of this array contains
C             the scaled direct transmission matrix D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     LOW     (output) INTEGER
C             The index of the lower end of the balanced submatrix of A.
C
C     IGH     (output) INTEGER
C             The index of the upper end of the balanced submatrix of A.
C
C     SCSTAT  (output) DOUBLE PRECISION array, dimension (N)
C             This array contains the information defining the
C             similarity transformations used to permute and balance
C             the state dynamics matrix A, as returned from the LAPACK
C             library routine DGEBAL.
C
C     SCIN    (output) DOUBLE PRECISION array, dimension (M)
C             Contains the scalars used to scale the system inputs so
C             that the columns of the final matrix B have norms roughly
C             equal to the column sums of the balanced matrix A
C             (see FURTHER COMMENTS).
C             The j-th input of the balanced state-space representation
C             is SCIN(j)*(j-th column of the permuted and balanced
C             input/state matrix B).
C
C     SCOUT   (output) DOUBLE PRECISION array, dimension (P)
C             Contains the scalars used to scale the system outputs so
C             that the rows of the final matrix C have norms roughly
C             equal to the row sum of the balanced matrix A.
C             The i-th output of the balanced state-space representation
C             is SCOUT(i)*(i-th row of the permuted and balanced
C             state/ouput matrix C).
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
C     Similarity transformations are used to permute the system states
C     and balance the corresponding row and column sum norms of a
C     submatrix of the state dynamics matrix A. These operations are
C     also applied to the input/state matrix B and the system inputs
C     are then scaled (see parameter SCIN) so that the columns of the
C     final matrix B have norms roughly equal to the column sum norm of
C     the balanced matrix A (see FURTHER COMMENTS).
C     The above operations are also applied to the matrix C, and the
C     system outputs are then scaled (see parameter SCOUT) so that the
C     rows of the final matrix C have norms roughly equal to the row sum
C     norm of the balanced matrix A (see FURTHER COMMENTS).
C     Finally, the (I,J)-th element of the direct transmission matrix D
C     is scaled as
C          D(I,J) = D(I,J)*(1.0/SCIN(J))*SCOUT(I), where I = 1,2,...,P
C     and J = 1,2,...,M.
C
C     Scaling performed to balance the row/column sum norms is by
C     integer powers of the machine base so as to avoid introducing
C     rounding errors.
C
C     REFERENCES
C
C     [1] Wilkinson, J.H. and Reinsch, C.
C         Handbook for Automatic Computation, (Vol II, Linear Algebra).
C         Springer-Verlag, 1971, (contribution II/11).
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     FURTHER COMMENTS
C
C     The columns (rows) of the final matrix B (matrix C) have norms
C     'roughly' equal to the column (row) sum norm of the balanced
C     matrix A, i.e.
C        size/BASE < abssum <= size
C     where
C        BASE   = the base of the arithmetic used on the computer, which
C                 can be obtained from the LAPACK Library routine
C                 DLAMCH;
C
C        size   = column or row sum norm of the balanced matrix A;
C        abssum = column sum norm of the balanced matrix B or row sum
C                 norm of the balanced matrix C.
C
C     The routine is BASE dependent.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TB01HD by T.W.C.Williams, Kingston
C     Polytechnic, United Kingdom, October 1982.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Balanced form, orthogonal transformation, similarity
C     transformation, state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           IGH, INFO, LDA, LDB, LDC, LDD, LOW, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), SCIN(*), SCOUT(*), SCSTAT(*)
C     .. Local Scalars ..
      INTEGER           I, J, K, KNEW, KOLD
      DOUBLE PRECISION  ACNORM, ARNORM, SCALE
C     .. External Functions ..
      DOUBLE PRECISION  DLANGE
      EXTERNAL          DLANGE
C     .. External Subroutines ..
      EXTERNAL          DGEBAL, DSCAL, DSWAP, TB01TY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
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
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -9
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01TD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, P ).EQ.0 ) THEN
         LOW = 1
         IGH = N
         RETURN
      END IF
C
C     Permute states, and balance a submatrix of A.
C
      CALL DGEBAL( 'Both', N, A, LDA, LOW, IGH, SCSTAT, INFO )
C
C     Use the information in SCSTAT on state scalings and reorderings
C     to transform B and C.
C
      DO 10 K = 1, N
         KOLD = K
         IF ( ( LOW.GT.KOLD ) .OR. ( KOLD.GT.IGH ) ) THEN
            IF ( KOLD.LT.LOW ) KOLD = LOW - KOLD
            KNEW = INT( SCSTAT(KOLD) )
            IF ( KNEW.NE.KOLD ) THEN
C
C              Exchange rows KOLD and KNEW of B.
C
               CALL DSWAP( M, B(KOLD,1), LDB, B(KNEW,1), LDB )
C
C              Exchange columns KOLD and KNEW of C.
C
               CALL DSWAP( P, C(1,KOLD), 1, C(1,KNEW), 1 )
            END IF
         END IF
   10 CONTINUE
C
      IF ( IGH.NE.LOW ) THEN
C
         DO 20 K = LOW, IGH
            SCALE = SCSTAT(K)
C
C           Scale the K-th row of permuted B.
C
            CALL DSCAL( M, ONE/SCALE, B(K,1), LDB )
C
C           Scale the K-th column of permuted C.
C
            CALL DSCAL( P, SCALE, C(1,K), 1 )
   20    CONTINUE
C
      END IF
C
C     Calculate the column and row sum norms of A.
C
      ACNORM = DLANGE( '1-norm', N, N, A, LDA, DWORK )
      ARNORM = DLANGE( 'I-norm', N, N, A, LDA, DWORK )
C
C     Scale the columns of B (i.e. inputs) to have norms roughly ACNORM.
C
      CALL TB01TY( 1, 0, 0, N, M, ACNORM, B, LDB, SCIN )
C
C     Scale the rows of C (i.e. outputs) to have norms roughly ARNORM.
C
      CALL TB01TY( 0, 0, 0, P, N, ARNORM, C, LDC, SCOUT )
C
C     Finally, apply these input and output scalings to D and set SCIN.
C
      DO 40 J = 1, M
         SCALE = SCIN(J)
C
         DO 30 I = 1, P
            D(I,J) = D(I,J)*( SCALE*SCOUT(I) )
   30    CONTINUE
C
         SCIN(J) = ONE/SCALE
   40 CONTINUE
C
      RETURN
C *** Last line of TB01TD ***
      END
