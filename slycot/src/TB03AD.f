      SUBROUTINE TB03AD( LERI, EQUIL, N, M, P, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, NR, INDEX, PCOEFF, LDPCO1, LDPCO2,
     $                   QCOEFF, LDQCO1, LDQCO2, VCOEFF, LDVCO1, LDVCO2,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To find a relatively prime left polynomial matrix representation
C     inv(P(s))*Q(s) or right polynomial matrix representation
C     Q(s)*inv(P(s)) with the same transfer matrix T(s) as that of a
C     given state-space representation, i.e.
C
C        inv(P(s))*Q(s) = Q(s)*inv(P(s)) = T(s) = C*inv(s*I-A)*B + D.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LERI    CHARACTER*1
C             Indicates whether the left polynomial matrix
C             representation or the right polynomial matrix
C             representation is required as follows:
C             = 'L':  A left matrix fraction is required;
C             = 'R':  A right matrix fraction is required.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to balance the triplet
C             (A,B,C), before computing a minimal state-space
C             representation, as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
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
C             On exit, the leading NR-by-NR part of this array contains
C             the upper block Hessenberg state dynamics matrix Amin of a
C             minimal realization for the original system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,MAX(M,P))
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B; the remainder
C             of the leading N-by-MAX(M,P) part is used as internal
C             workspace.
C             On exit, the leading NR-by-M part of this array contains
C             the transformed input/state matrix Bmin.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C; the remainder
C             of the leading MAX(M,P)-by-N part is used as internal
C             workspace.
C             On exit, the leading P-by-NR part of this array contains
C             the transformed state/output matrix Cmin.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,MAX(M,P))
C             The leading P-by-M part of this array must contain the
C             original direct transmission matrix D; the remainder of
C             the leading MAX(M,P)-by-MAX(M,P) part is used as internal
C             workspace.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,M,P).
C
C     NR      (output) INTEGER
C             The order of the minimal state-space representation
C             (Amin,Bmin,Cmin).
C
C     INDEX   (output) INTEGER array, dimension (P), if LERI = 'L', or
C                                     dimension (M), if LERI = 'R'.
C             If LERI = 'L', INDEX(I), I = 1,2,...,P, contains the
C             maximum degree of the polynomials in the I-th row of the
C             denominator matrix P(s) of the left polynomial matrix
C             representation.
C             These elements are ordered so that
C             INDEX(1) >= INDEX(2) >= ... >= INDEX(P).
C             If LERI = 'R', INDEX(I), I = 1,2,...,M, contains the
C             maximum degree of the polynomials in the I-th column of
C             the denominator matrix P(s) of the right polynomial
C             matrix representation.
C             These elements are ordered so that
C             INDEX(1) >= INDEX(2) >= ... >= INDEX(M).
C
C     PCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDPCO1,LDPCO2,N+1)
C             If LERI = 'L' then porm = P, otherwise porm = M.
C             The leading porm-by-porm-by-kpcoef part of this array
C             contains the coefficients of the denominator matrix P(s),
C             where kpcoef = MAX(INDEX(I)) + 1.
C             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1)
C             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if
C             LERI = 'L' then iorj = I, otherwise iorj = J.
C             Thus for LERI = 'L', P(s) =
C             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...).
C
C     LDPCO1  INTEGER
C             The leading dimension of array PCOEFF.
C             LDPCO1 >= MAX(1,P), if LERI = 'L';
C             LDPCO1 >= MAX(1,M), if LERI = 'R'.
C
C     LDPCO2  INTEGER
C             The second dimension of array PCOEFF.
C             LDPCO2 >= MAX(1,P), if LERI = 'L';
C             LDPCO2 >= MAX(1,M), if LERI = 'R'.
C
C     QCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDQCO1,LDQCO2,N+1)
C             If LERI = 'L' then porp = M, otherwise porp = P.
C             If LERI = 'L', the leading porm-by-porp-by-kpcoef part
C             of this array contains the coefficients of the numerator
C             matrix Q(s).
C             If LERI = 'R', the leading porp-by-porm-by-kpcoef part
C             of this array contains the coefficients of the numerator
C             matrix Q(s).
C             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K).
C
C     LDQCO1  INTEGER
C             The leading dimension of array QCOEFF.
C             LDQCO1 >= MAX(1,P),   if LERI = 'L';
C             LDQCO1 >= MAX(1,M,P), if LERI = 'R'.
C
C     LDQCO2  INTEGER
C             The second dimension of array QCOEFF.
C             LDQCO2 >= MAX(1,M),   if LERI = 'L';
C             LDQCO2 >= MAX(1,M,P), if LERI = 'R'.
C
C     VCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDVCO1,LDVCO2,N+1)
C             The leading porm-by-NR-by-kpcoef part of this array
C             contains the coefficients of the intermediate matrix V(s).
C             VCOEFF(I,J,K) is defined as for PCOEFF(I,J,K).
C
C     LDVCO1  INTEGER
C             The leading dimension of array VCOEFF.
C             LDVCO1 >= MAX(1,P), if LERI = 'L';
C             LDVCO1 >= MAX(1,M), if LERI = 'R'.
C
C     LDVCO2  INTEGER
C             The second dimension of array VCOEFF.  LDVCO2 >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determination when
C             transforming (A, B, C). If the user sets TOL > 0, then
C             the given value of TOL is used as a lower bound for the
C             reciprocal condition number (see the description of the
C             argument RCOND in the SLICOT routine MB03OD);  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance
C             (determined by the SLICOT routine TB01UD) is used instead.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N+MAX(M,P))
C             On exit, if INFO = 0, the first nonzero elements of
C             IWORK(1:N) return the orders of the diagonal blocks of A.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P), PM*(PM + 2))
C             where  PM = P, if LERI = 'L';
C                    PM = M, if LERI = 'R'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if a singular matrix was encountered during the
C                   computation of V(s);
C             = 2:  if a singular matrix was encountered during the
C                   computation of P(s).
C
C     METHOD
C
C     The method for a left matrix fraction will be described here:
C     right matrix fractions are dealt with by constructing a left
C     fraction for the dual of the original system. The first step is to
C     obtain, by means of orthogonal similarity transformations, a
C     minimal state-space representation (Amin,Bmin,Cmin,D) for the
C     original system (A,B,C,D), where Amin is lower block Hessenberg
C     with all its superdiagonal blocks upper triangular and Cmin has
C     all but its first rank(C) columns zero.  The number and dimensions
C     of the blocks of Amin now immediately yield the row degrees of
C     P(s) with P(s) row proper: furthermore, the P-by-NR polynomial
C     matrix V(s) (playing a similar role to S(s) in Wolovich's
C     Structure Theorem) can be calculated a column block at a time, in
C     reverse order, from Amin. P(s) is then found as if it were the
C     O-th column block of V(s) (using Cmin as well as Amin), while
C     Q(s) = (V(s) * Bmin) + (P(s) * D). Finally, a special similarity
C     transformation is used to put Amin in an upper block Hessenberg
C     form.
C
C     REFERENCES
C
C     [1] Williams, T.W.C.
C         An Orthogonal Structure Theorem for Linear Systems.
C         Kingston Polytechnic Control Systems Research Group,
C         Internal Report 82/2, July 1982.
C
C     [2] Patel, R.V.
C         On Computing Matrix Fraction Descriptions and Canonical
C         Forms of Linear Time-Invariant Systems.
C         UMIST Control Systems Centre Report 489, 1980.
C         (Algorithms 1 and 2, extensively modified).
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998.
C     Supersedes Release 3.0 routine TB01SD.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000.
C
C     KEYWORDS
C
C     Canonical form, coprime matrix fraction, dual system, elementary
C     polynomial operations, Hessenberg form, minimal realization,
C     orthogonal transformation, polynomial matrix, state-space
C     representation, transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL, LERI
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDPCO1, LDPCO2,
     $                  LDQCO1, LDQCO2, LDVCO1, LDVCO2, LDWORK, M, N,
     $                  NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), PCOEFF(LDPCO1,LDPCO2,*),
     $                  QCOEFF(LDQCO1,LDQCO2,*), VCOEFF(LDVCO1,LDVCO2,*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LLERIL, LLERIR
      INTEGER           I, IC, IFIRST, INDBLK, INPLUS, IOFF, IRANKC,
     $                  ISTART, ISTOP, ITAU, IZ, JOFF, JWORK, K, KMAX,
     $                  KPCOEF, KPLUS, KWORK, LDWRIC, MAXMP, MPLIM,
     $                  MWORK, NCOL, NCONT, NREFLC, NROW, PWORK, WRKOPT
      DOUBLE PRECISION  MAXRED
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07MD, DGEMM, DGEQRF, DGETRF, DLACPY, DLASET,
     $                  DORMQR, DTRSM, MA02GD, TB01ID, TB01UD, TB01YD,
     $                  TB03AY, TC01OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C     .. Executable Statements ..
C
      INFO   = 0
      LLERIL = LSAME( LERI,  'L' )
      LLERIR = LSAME( LERI,  'R' )
      LEQUIL = LSAME( EQUIL, 'S' )
      MAXMP  = MAX( M, P )
      MPLIM  = MAX( 1, MAXMP )
      IF ( LLERIR ) THEN
C
C        Initialization for right matrix fraction.
C
         PWORK = M
         MWORK = P
      ELSE
C
C        Initialization for left matrix fraction.
C
         PWORK = P
         MWORK = M
      END IF
C
C     Test the input scalar arguments.
C
      IF( .NOT.LLERIL .AND. .NOT.LLERIR ) THEN
         INFO = -1
      ELSE IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MPLIM ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MPLIM ) THEN
         INFO = -13
      ELSE IF( LDPCO1.LT.MAX( 1, PWORK ) ) THEN
         INFO = -17
      ELSE IF( LDPCO2.LT.MAX( 1, PWORK ) ) THEN
         INFO = -18
      ELSE IF( LDQCO1.LT.MAX( 1, PWORK ) .OR. LLERIR .AND.
     $         LDQCO1.LT.MPLIM ) THEN
         INFO = -20
      ELSE IF( LDQCO2.LT.MAX( 1, MWORK ) .OR. LLERIR .AND.
     $         LDQCO2.LT.MPLIM ) THEN
         INFO = -21
      ELSE IF( LDVCO1.LT.MAX( 1, PWORK ) ) THEN
         INFO = -23
      ELSE IF( LDVCO2.LT.MAX( 1, N ) ) THEN
         INFO = -24
      ELSE IF( LDWORK.LT.MAX( 1, N + MAX( N, 3*MAXMP ),
     $                           PWORK*( PWORK + 2 ) ) ) THEN
         INFO = -28
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB03AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, P ).EQ.0 ) THEN
         NR = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( LLERIR ) THEN
C
C        For right matrix fraction, obtain dual system.
C
         CALL AB07MD( 'D', N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                INFO )
      END IF
C
C     Obtain minimal realization, in canonical form, for this system.
C     Part of the code in SLICOT routine TB01PD is included in-line
C     here. (TB01PD cannot be directly used.)
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
C     If required, balance the triplet (A,B,C) (default MAXRED).
C     Workspace: need N.
C
      IF ( LEQUIL ) THEN
         MAXRED = ZERO
         CALL TB01ID( 'A', N, MWORK, PWORK, MAXRED, A, LDA, B, LDB, C,
     $                LDC, DWORK, INFO )
      END IF
C
      IZ    = 1
      ITAU  = 1
      JWORK = ITAU + N
C
C     Separate out controllable subsystem (of order NCONT):
C     A <-- Z'*A*Z,  B <-- Z'*B,  C <-- C*Z.
C
C     Workspace: need   N + MAX(N, 3*MWORK, PWORK).
C                prefer larger.
C
      CALL TB01UD( 'No Z', N, MWORK, PWORK, A, LDA, B, LDB, C, LDC,
     $             NCONT, INDBLK, IWORK, DWORK(IZ), 1, DWORK(ITAU), TOL,
     $             IWORK(N+1), DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
      WRKOPT = INT( DWORK(JWORK) ) + JWORK - 1
C
C     Separate out the observable subsystem (of order NR):
C     Form the dual of the subsystem of order NCONT (which is
C     controllable), leaving rest as it is.
C
      CALL AB07MD( 'Z', NCONT, MWORK, PWORK, A, LDA, B, LDB, C, LDC,
     $             DWORK, 1, INFO )
C
C     And separate out the controllable part of this dual subsystem.
C
C     Workspace: need   NCONT + MAX(NCONT, 3*PWORK, MWORK).
C                prefer larger.
C
      CALL TB01UD( 'No Z', NCONT, PWORK, MWORK, A, LDA, B, LDB, C, LDC,
     $             NR, INDBLK, IWORK, DWORK(IZ), 1, DWORK(ITAU), TOL,
     $             IWORK(N+1), DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     Retranspose, giving controllable and observable (i.e. minimal)
C     part of original system.
C
      CALL AB07MD( 'Z', NR, PWORK, MWORK, A, LDA, B, LDB, C, LDC, DWORK,
     $             1, INFO )
C
C     Annihilate the trailing components of IWORK(1:N).
C
      DO 10 I = INDBLK + 1, N
         IWORK(I) = 0
   10 CONTINUE
C
C     Initialize polynomial matrices P(s), Q(s) and V(s) to zero.
C
      DO 20 K = 1, N + 1
         CALL DLASET( 'Full', PWORK, PWORK, ZERO, ZERO, PCOEFF(1,1,K),
     $                LDPCO1 )
         CALL DLASET( 'Full', PWORK, MWORK, ZERO, ZERO, QCOEFF(1,1,K),
     $                LDQCO1 )
         CALL DLASET( 'Full', PWORK, NR, ZERO, ZERO, VCOEFF(1,1,K),
     $                LDVCO1 )
   20 CONTINUE
C
C     Finish initializing V(s), and set up row degrees of P(s).
C
      INPLUS = INDBLK + 1
      ISTART = 1
      JOFF = NR
C
      DO 40 K = 1, INDBLK
         KWORK = INPLUS - K
         KPLUS = KWORK + 1
         ISTOP = IWORK(KWORK)
         JOFF  = JOFF - ISTOP
C
         DO 30 I = ISTART, ISTOP
            INDEX(I) = KWORK
            VCOEFF(I,JOFF+I,KPLUS) = ONE
   30    CONTINUE
C
         ISTART = ISTOP + 1
   40 CONTINUE
C
C     ISTART = IWORK(1)+1 now: if .LE. PWORK, set up final rows of P(s).
C
      DO 50 I = ISTART, PWORK
         INDEX(I) = 0
         PCOEFF(I,I,1) = ONE
   50 CONTINUE
C
C     Triangularize the superdiagonal blocks of Amin.
C
      NROW = IWORK(INDBLK)
      IOFF = NR - NROW
      KMAX = INDBLK - 1
      ITAU = 1
      IFIRST = 0
      IF ( INDBLK.GT.2 ) IFIRST = IOFF - IWORK(KMAX)
C
C     QR decomposition of each superdiagonal block of A in turn
C     (done in reverse order to preserve upper triangular blocks in A).
C
      DO 60 K = 1, KMAX
C
C        Calculate dimensions of new block & its position in A.
C
         KWORK = INDBLK - K
         NCOL = NROW
         NROW = IWORK(KWORK)
         JOFF = IOFF
         IOFF = IOFF - NROW
         NREFLC = MIN( NROW, NCOL )
         JWORK  = ITAU + NREFLC
         IF ( KWORK.GE.2 ) IFIRST = IFIRST - IWORK(KWORK-1)
C
C        Find QR decomposition of this (full rank) block:
C        block = QR.  No pivoting is needed.
C
C        Workspace: need   MIN(NROW,NCOL) + NCOL;
C                   prefer MIN(NROW,NCOL) + NCOL*NB.
C
         CALL DGEQRF( NROW, NCOL, A(IOFF+1,JOFF+1), LDA, DWORK(ITAU),
     $                DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Premultiply appropriate row block of A by Q'.
C
C        Workspace: need   MIN(NROW,NCOL) + JOFF;
C                   prefer MIN(NROW,NCOL) + JOFF*NB.
C
         CALL DORMQR( 'Left', 'Transpose', NROW, JOFF, NREFLC,
     $                A(IOFF+1,JOFF+1), LDA, DWORK(ITAU), A(IOFF+1,1),
     $                LDA, DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Premultiply appropriate row block of B by Q' also.
C
C        Workspace: need   MIN(NROW,NCOL) + MWORK;
C                   prefer MIN(NROW,NCOL) + MWORK*NB.
C
         CALL DORMQR( 'Left', 'Transpose', NROW, MWORK, NREFLC,
     $                A(IOFF+1,JOFF+1), LDA, DWORK(ITAU), B(IOFF+1,1),
     $                LDB, DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        And postmultiply the non-zero part of appropriate column
C        block of A by Q.
C
C        Workspace: need   MIN(NROW,NCOL) + NR;
C                   prefer MIN(NROW,NCOL) + NR*NB.
C
         CALL DORMQR( 'Right', 'No Transpose', NR-IFIRST, NROW, NREFLC,
     $                A(IOFF+1,JOFF+1), LDA, DWORK(ITAU),
     $                A(IFIRST+1,IOFF+1), LDA, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Annihilate the lower triangular part of the block in A.
C
         IF ( K.NE.KMAX .AND. NROW.GT.1 )
     $      CALL DLASET( 'Lower', NROW-1, NCOL, ZERO, ZERO,
     $                   A(IOFF+2,JOFF+1), LDA )
C
   60 CONTINUE
C
C     Finally: postmultiply non-zero columns of C by Q (K = KMAX).
C
C     Workspace: need   MIN(NROW,NCOL) + PWORK;
C                prefer MIN(NROW,NCOL) + PWORK*NB.
C
      CALL DORMQR( 'Right', 'No Transpose', PWORK, NROW, NREFLC,
     $             A(IOFF+1,JOFF+1), LDA, DWORK(ITAU), C, LDC,
     $             DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C     Annihilate the lower triangular part of the block in A.
C
      IF ( NROW.GT.1 )
     $   CALL DLASET( 'Lower', NROW-1, NCOL, ZERO, ZERO,
     $                A(IOFF+2,JOFF+1), LDA )
C
C     Calculate the (PWORK x NR) polynomial matrix V(s) ...
C
      CALL TB03AY( NR, A, LDA, INDBLK, IWORK, VCOEFF, LDVCO1, LDVCO2,
     $             PCOEFF, LDPCO1, LDPCO2, INFO)
C
      IF ( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      ELSE
C
C        And then use this matrix to calculate P(s): first store
C        C1 from C.
C
         IC = 1
         IRANKC = IWORK(1)
         LDWRIC = MAX( 1, PWORK )
         CALL DLACPY( 'Full', PWORK, IRANKC, C, LDC, DWORK(IC), LDWRIC )
C
         IF ( IRANKC.LT.PWORK ) THEN
C
C           rank(C) .LT. PWORK: obtain QR decomposition of C1,
C           giving R and Q.
C
C           Workspace: need   PWORK*IRANKC + 2*IRANKC;
C                      prefer PWORK*IRANKC +   IRANKC + IRANKC*NB.
C
            ITAU  = IC + LDWRIC*IRANKC
            JWORK = ITAU + IRANKC
C
            CALL DGEQRF( PWORK, IRANKC, DWORK(IC), LDWRIC, DWORK(ITAU),
     $                   DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
            WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C           First IRANKC rows of Pbar(s) are given by Wbar(s) * inv(R).
C           Check for zero diagonal elements of R.
C
            DO 70 I = 1, IRANKC
               IF ( DWORK(IC+(I-1)*LDWRIC+I-1).EQ.ZERO ) THEN
C
C                 Error return.
C
                  INFO = 2
                  RETURN
               END IF
   70       CONTINUE
C
            NROW = IRANKC
C
            DO 80 K = 1, INPLUS
               CALL DTRSM( 'Right', 'Upper', 'No Transpose', 'Non-unit',
     $                     NROW, IRANKC, ONE, DWORK(IC), LDWRIC,
     $                     PCOEFF(1,1,K), LDPCO1 )
               NROW = IWORK(K)
   80       CONTINUE
C
C           P(s) itself is now given by Pbar(s) * Q'.
C
            NROW = PWORK
C
            DO 90 K = 1, INPLUS
C
C              Workspace: need   PWORK*IRANKC + IRANKC + NROW;
C                         prefer PWORK*IRANKC + IRANKC + NROW*NB.
C
               CALL DORMQR( 'Right', 'Transpose', NROW, PWORK, IRANKC,
     $                      DWORK(IC), LDWRIC, DWORK(ITAU),
     $                      PCOEFF(1,1,K), LDPCO1, DWORK(JWORK),
     $                      LDWORK-JWORK+1, INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
               NROW = IWORK(K)
   90       CONTINUE
C
         ELSE
C
C           Special case rank(C) = PWORK, full:
C           no QR decomposition (P(s)=Wbar(s)*inv(C1)).
C
            CALL DGETRF( PWORK, PWORK, DWORK(IC), LDWRIC, IWORK(N+1),
     $                   INFO )
C
            IF ( INFO.NE.0 ) THEN
C
C              Error return.
C
               INFO = 2
               RETURN
            ELSE
C
               NROW = IRANKC
C
C              Workspace: need   PWORK*IRANKC + N.
C
               DO 100 K = 1, INPLUS
                  CALL DTRSM( 'Right', 'Upper', 'No Transpose',
     $                        'Non-unit', NROW, PWORK, ONE, DWORK(IC),
     $                        LDWRIC, PCOEFF(1,1,K), LDPCO1 )
                  CALL DTRSM( 'Right', 'Lower', 'No Transpose', 'Unit',
     $                        NROW, PWORK, ONE, DWORK(IC), LDWRIC,
     $                        PCOEFF(1,1,K), LDPCO1 )
                  CALL MA02GD( NROW, PCOEFF(1,1,K), LDPCO1, 1, PWORK,
     $                         IWORK(N+1), -1 )
                  NROW = IWORK(K)
  100          CONTINUE
            END IF
         END IF
C
C        Finally, Q(s) = V(s) * B + P(s) * D can now be evaluated.
C
         NROW = PWORK
C
         DO 110 K = 1, INPLUS
            CALL DGEMM( 'No transpose', 'No transpose', NROW, MWORK,
     $                  NR, ONE, VCOEFF(1,1,K), LDVCO1, B, LDB, ZERO,
     $                  QCOEFF(1,1,K), LDQCO1 )
            CALL DGEMM( 'No transpose', 'No transpose', NROW, MWORK,
     $                  PWORK, ONE, PCOEFF(1,1,K), LDPCO1, D, LDD, ONE,
     $                  QCOEFF(1,1,K), LDQCO1 )
            NROW = IWORK(K)
  110    CONTINUE
C
      END IF
C
      IF ( LLERIR ) THEN
C
C        For right matrix fraction, return to original (dual of dual)
C        system.
C
         CALL AB07MD( 'Z', NR, MWORK, PWORK, A, LDA, B, LDB, C, LDC,
     $                DWORK, 1, INFO )
C
C        Also, obtain the dual of the polynomial matrix representation.
C
         KPCOEF = 0
C
         DO 120 I = 1, PWORK
            KPCOEF = MAX( KPCOEF, INDEX(I) )
  120    CONTINUE
C
         KPCOEF = KPCOEF + 1
         CALL TC01OD( 'L', MWORK, PWORK, KPCOEF, PCOEFF, LDPCO1,
     $                LDPCO2, QCOEFF, LDQCO1, LDQCO2, INFO )
      ELSE
C
C        Reorder the rows and columns of the system, to get an upper
C        block Hessenberg matrix A of the minimal system.
C
         CALL TB01YD( NR, M, P, A, LDA, B, LDB, C, LDC, INFO )
      END IF
C
C     Set optimal workspace dimension.
C
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of TB03AD ***
      END
