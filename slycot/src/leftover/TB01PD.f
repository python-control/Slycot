      SUBROUTINE TB01PD( JOB, EQUIL, N, M, P, A, LDA, B, LDB, C, LDC,
     $                   NR, TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To find a reduced (controllable, observable, or minimal) state-
C     space representation (Ar,Br,Cr) for any original state-space
C     representation (A,B,C). The matrix Ar is in upper block
C     Hessenberg form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates whether the user wishes to remove the
C             uncontrollable and/or unobservable parts as follows:
C             = 'M':  Remove both the uncontrollable and unobservable
C                     parts to get a minimal state-space representation;
C             = 'C':  Remove the uncontrollable part only to get a
C                     controllable state-space representation;
C             = 'O':  Remove the unobservable part only to get an
C                     observable state-space representation.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily balance
C             the triplet (A,B,C) as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation, i.e.
C             the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.   P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, the leading NR-by-NR part of this array contains
C             the upper block Hessenberg state dynamics matrix Ar of a
C             minimal, controllable, or observable realization for the
C             original system, depending on the value of JOB, JOB = 'M',
C             JOB = 'C', or JOB = 'O', respectively.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M),
C             if JOB = 'C', or (LDB,MAX(M,P)), otherwise.
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B; if JOB = 'M',
C             or JOB = 'O', the remainder of the leading N-by-MAX(M,P)
C             part is used as internal workspace.
C             On exit, the leading NR-by-M part of this array contains
C             the transformed input/state matrix Br of a minimal,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'M',
C             JOB = 'C', or JOB = 'O', respectively.
C             If JOB = 'C', only the first IWORK(1) rows of B are
C             nonzero.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C; if JOB = 'M',
C             or JOB = 'O', the remainder of the leading MAX(M,P)-by-N
C             part is used as internal workspace.
C             On exit, the leading P-by-NR part of this array contains
C             the transformed state/output matrix Cr of a minimal,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'M',
C             JOB = 'C', or JOB = 'O', respectively.
C             If JOB = 'M', or JOB = 'O', only the last IWORK(1) columns
C             (in the first NR columns) of C are nonzero.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,M,P) if N > 0.
C             LDC >= 1          if N = 0.
C
C     NR      (output) INTEGER
C             The order of the reduced state-space representation
C             (Ar,Br,Cr) of a minimal, controllable, or observable
C             realization for the original system, depending on
C             JOB = 'M', JOB = 'C', or JOB = 'O'.
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
C             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P)).
C             For optimum performance LDWORK should be larger.
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
C     If JOB = 'M', the matrices A and B are operated on by orthogonal
C     similarity transformations (made up of products of Householder
C     transformations) so as to produce an upper block Hessenberg matrix
C     A1 and a matrix B1 with all but its first rank(B) rows zero; this
C     separates out the controllable part of the original system.
C     Applying the same algorithm to the dual of this subsystem,
C     therefore separates out the controllable and observable (i.e.
C     minimal) part of the original system representation, with the
C     final Ar upper block Hessenberg (after using pertransposition).
C     If JOB = 'C', or JOB = 'O', only the corresponding part of the
C     above procedure is applied.
C
C     REFERENCES
C
C     [1] Van Dooren, P.
C         The Generalized Eigenstructure Problem in Linear System
C         Theory. (Algorithm 1)
C         IEEE Trans. Auto. Contr., AC-26, pp. 111-129, 1981.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998.
C
C     REVISIONS
C
C     A. Varga, DLR Oberpfaffenhofen, July 1998.
C     A. Varga, DLR Oberpfaffenhofen, April 28, 1999.
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004.
C
C     KEYWORDS
C
C     Hessenberg form, minimal realization, multivariable system,
C     orthogonal transformation, state-space model, state-space
C     representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDIZ
      PARAMETER         ( LDIZ = 1 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL, JOB
      INTEGER           INFO, LDA, LDB, LDC, LDWORK, M, N, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LNJOBC, LNJOBO
      INTEGER           I, INDCON, ITAU, IZ, JWORK, KL, MAXMP, NCONT,
     $                  WRKOPT
      DOUBLE PRECISION  MAXRED
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07MD, TB01ID, TB01UD, TB01XD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      MAXMP = MAX( M, P )
      LNJOBC = .NOT.LSAME( JOB,   'C' )
      LNJOBO = .NOT.LSAME( JOB,   'O' )
      LEQUIL =      LSAME( EQUIL, 'S' )
C
C     Test the input scalar arguments.
C
      IF( LNJOBC .AND. LNJOBO .AND. .NOT.LSAME( JOB, 'M' ) ) THEN
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
      ELSE IF( LDC.LT.1 .OR. ( N.GT.0 .AND. LDC.LT.MAXMP ) ) THEN
         INFO = -11
      ELSE IF( LDWORK.LT.MAX( 1, N + MAX( N, 3*MAXMP ) ) ) THEN
         INFO = -16
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB01PD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. ( LNJOBC .AND. MIN( N, P ).EQ.0 ) .OR.
     $                 ( LNJOBO .AND. MIN( N, M ).EQ.0 ) ) THEN
         NR = 0
C
         DO 5 I = 1, N
            IWORK(I) = 0
    5    CONTINUE
C
         DWORK(1) = ONE
         RETURN
      END IF
C
C     If required, balance the triplet (A,B,C) (default MAXRED).
C     Workspace: need N.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the code,
C     as well as the preferred amount for good performance.)
C
      IF ( LEQUIL ) THEN
         MAXRED = ZERO
         CALL TB01ID( 'A', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
         WRKOPT = N
      ELSE
         WRKOPT = 1
      END IF
C
      IZ    = 1
      ITAU  = 1
      JWORK = ITAU + N
      IF ( LNJOBO ) THEN
C
C        Separate out controllable subsystem (of order NCONT):
C        A <-- Z'*A*Z,  B <-- Z'*B,  C <-- C*Z.
C
C        Workspace: need   N + MAX(N, 3*M, P).
C                   prefer larger.
C
         CALL TB01UD( 'No Z', N, M, P, A, LDA, B, LDB, C, LDC, NCONT,
     $                INDCON, IWORK, DWORK(IZ), LDIZ, DWORK(ITAU), TOL,
     $                IWORK(N+1), DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = INT( DWORK(JWORK) ) + JWORK - 1
      ELSE
         NCONT = N
      END IF
C
      IF ( LNJOBC ) THEN
C
C        Separate out the observable subsystem (of order NR):
C        Form the dual of the subsystem of order NCONT (which is
C        controllable, if JOB = 'M'), leaving rest as it is.
C
         CALL AB07MD( 'Z', NCONT, M, P, A, LDA, B, LDB, C, LDC, DWORK,
     $                1, INFO )
C
C        And separate out the controllable part of this dual subsystem.
C
C        Workspace: need   NCONT + MAX(NCONT, 3*P, M).
C                   prefer larger.
C
         CALL TB01UD( 'No Z', NCONT, P, M, A, LDA, B, LDB, C, LDC, NR,
     $                INDCON, IWORK, DWORK(IZ), LDIZ, DWORK(ITAU), TOL,
     $                IWORK(N+1), DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Transpose and reorder (to get a block upper Hessenberg
C        matrix A), giving, for JOB = 'M', the controllable and
C        observable (i.e., minimal) part of original system.
C
         IF( INDCON.GT.0 ) THEN
            KL = IWORK(1) - 1
            IF ( INDCON.GE.2 )
     $         KL = KL + IWORK(2)
         ELSE
            KL = 0
         END IF
         CALL TB01XD( 'Zero D', NR, P, M, KL, MAX( 0, NR-1 ), A, LDA,
     $                B, LDB, C, LDC, DWORK, 1, INFO )
      ELSE
         NR = NCONT
      END IF
C
C     Annihilate the trailing components of IWORK(1:N).
C
      DO 10 I = INDCON + 1, N
         IWORK(I) = 0
   10 CONTINUE
C
C     Set optimal workspace dimension.
C
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of TB01PD ***
      END
