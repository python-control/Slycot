      SUBROUTINE TC04AD( LERI, M, P, INDEX, PCOEFF, LDPCO1, LDPCO2,
     $                   QCOEFF, LDQCO1, LDQCO2, N, RCOND, A, LDA, B,
     $                   LDB, C, LDC, D, LDD, IWORK, DWORK, LDWORK,
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
C     To find a state-space representation (A,B,C,D) with the same
C     transfer matrix T(s) as that of a given left or right polynomial
C     matrix representation, i.e.
C
C        C*inv(sI-A)*B + D = T(s) = inv(P(s))*Q(s) = Q(s)*inv(P(s)).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LERI    CHARACTER*1
C             Indicates whether a left polynomial matrix representation
C             or a right polynomial matrix representation is input as
C             follows:
C             = 'L':  A left matrix fraction is input;
C             = 'R':  A right matrix fraction is input.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     INDEX   (input) INTEGER array, dimension (MAX(M,P))
C             If LERI = 'L', INDEX(I), I = 1,2,...,P, must contain the
C             maximum degree of the polynomials in the I-th row of the
C             denominator matrix P(s) of the given left polynomial
C             matrix representation.
C             If LERI = 'R', INDEX(I), I = 1,2,...,M, must contain the
C             maximum degree of the polynomials in the I-th column of
C             the denominator matrix P(s) of the given right polynomial
C             matrix representation.
C
C     PCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDPCO1,LDPCO2,kpcoef), where kpcoef = MAX(INDEX(I)) + 1.
C             If LERI = 'L' then porm = P, otherwise porm = M.
C             The leading porm-by-porm-by-kpcoef part of this array must
C             contain the coefficients of the denominator matrix P(s).
C             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1)
C             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if
C             LERI = 'L' then iorj = I, otherwise iorj = J.
C             Thus for LERI = 'L', P(s) =
C             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...).
C             If LERI = 'R', PCOEFF is modified by the routine but
C             restored on exit.
C
C     LDPCO1  INTEGER
C             The leading dimension of array PCOEFF.
C             LDPCO1 >= MAX(1,P) if LERI = 'L',
C             LDPCO1 >= MAX(1,M) if LERI = 'R'.
C
C     LDPCO2  INTEGER
C             The second dimension of array PCOEFF.
C             LDPCO2 >= MAX(1,P) if LERI = 'L',
C             LDPCO2 >= MAX(1,M) if LERI = 'R'.
C
C     QCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDQCO1,LDQCO2,kpcoef)
C             If LERI = 'L' then porp = M, otherwise porp = P.
C             The leading porm-by-porp-by-kpcoef part of this array must
C             contain the coefficients of the numerator matrix Q(s).
C             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K).
C             If LERI = 'R', QCOEFF is modified by the routine but
C             restored on exit.
C
C     LDQCO1  INTEGER
C             The leading dimension of array QCOEFF.
C             LDQCO1 >= MAX(1,P)   if LERI = 'L',
C             LDQCO1 >= MAX(1,M,P) if LERI = 'R'.
C
C     LDQCO2  INTEGER
C             The second dimension of array QCOEFF.
C             LDQCO2 >= MAX(1,M)   if LERI = 'L',
C             LDQCO2 >= MAX(1,M,P) if LERI = 'R'.
C
C     N       (output) INTEGER
C             The order of the resulting state-space representation.
C                          porm
C             That is, N = SUM INDEX(I).
C                          I=1
C
C     RCOND   (output) DOUBLE PRECISION
C             The estimated reciprocal of the condition number of the
C             leading row (if LERI = 'L') or the leading column (if
C             LERI = 'R') coefficient matrix of P(s).
C             If RCOND is nearly zero, P(s) is nearly row or column
C             non-proper.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the state
C             dynamics matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P))
C             The leading N-by-M part of this array contains the
C             input/state matrix B; the remainder of the leading
C             N-by-MAX(M,P) part is used as internal workspace.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array contains the
C             state/output matrix C; the remainder of the leading
C             MAX(M,P)-by-N part is used as internal workspace.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M,P).
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,MAX(M,P))
C             The leading P-by-M part of this array contains the direct
C             transmission matrix D; the remainder of the leading
C             MAX(M,P)-by-MAX(M,P) part is used as internal workspace.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,M,P).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*MAX(M,P))
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,MAX(M,P)*(MAX(M,P)+4)).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if P(s) is not row (if LERI = 'L') or column
C                   (if LERI = 'R') proper. Consequently, no state-space
C                   representation is calculated.
C
C     METHOD
C
C     The method for a left matrix fraction will be described here;
C     right matrix fractions are dealt with by obtaining the dual left
C     polynomial matrix representation and constructing an equivalent
C     state-space representation for this. The first step is to check
C     if the denominator matrix P(s) is row proper; if it is not then
C     the routine returns with the Error Indicator (INFO) set to 1.
C     Otherwise, Wolovich's Observable  Structure Theorem is used to
C     construct a state-space representation (A,B,C,D) in observable
C     companion form. The sizes of the blocks of matrix A and matrix C
C     here are precisely the row degrees of P(s), while their
C     'non-trivial' columns are given easily from its coefficients.
C     Similarly, the matrix D is obtained from the leading coefficients
C     of P(s) and of the numerator matrix Q(s), while matrix B is given
C     by the relation Sbar(s)B = Q(s) - P(s)D, where Sbar(s) is a
C     polynomial matrix whose (j,k)(th) element is given by
C
C                  j-u(k-1)-1
C               ( s           , j = u(k-1)+1,u(k-1)+2,....,u(k)
C     Sbar    = (
C        j,k    (           0 , otherwise
C
C             k
C     u(k) = SUM d , k = 1,2,...,M and d ,d ,...,d  are the
C            i=1  i                     1  2      M
C     controllability indices. For convenience in solving this, C' and B
C     are initially set up to contain the coefficients of P(s) and Q(s),
C     respectively, stored by rows.
C
C     REFERENCES
C
C     [1] Wolovich, W.A.
C         Linear Multivariate Systems, (Theorem 4.3.3).
C         Springer-Verlag, 1974.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TC01BD by T.W.C.Williams, Kingston
C     Polytechnic, United Kingdom, March 1982.
C
C     REVISIONS
C
C     February 22, 1998 (changed the name of TC01ND).
C     May 12, 1998.
C
C     KEYWORDS
C
C     Coprime matrix fraction, elementary polynomial operations,
C     polynomial matrix, state-space representation, transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         LERI
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDPCO1, LDPCO2,
     $                  LDQCO1, LDQCO2, LDWORK, M, N, P
      DOUBLE PRECISION  RCOND
C     .. Array Arguments ..
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), PCOEFF(LDPCO1,LDPCO2,*),
     $                  QCOEFF(LDQCO1,LDQCO2,*)
C     .. Local Scalars ..
      LOGICAL           LLERI
      INTEGER           I, IA, IBIAS, J, JA, JC, JW, JWORK, LDW, K,
     $                  KPCOEF, KSTOP, MAXIND, MINDEX, MWORK, PWORK,
     $                  WRKOPT
      DOUBLE PRECISION  DWNORM
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          LSAME, DLAMCH, DLANGE
C     .. External Subroutines ..
      EXTERNAL          AB07MD, DCOPY, DGECON, DGEMM, DGETRF, DGETRI,
     $                  DGETRS, DLACPY, DLASET, TC01OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX
C     .. Executable Statements ..
C
      INFO = 0
      LLERI = LSAME( LERI, 'L' )
      MINDEX = MAX( M, P )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LLERI .AND. .NOT.LSAME( LERI, 'R' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( LLERI .AND. LDPCO1.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO1.LT.MAX( 1, M ) ) ) THEN
         INFO = -6
      ELSE IF( ( LLERI .AND. LDPCO2.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO2.LT.MAX( 1, M ) ) ) THEN
         INFO = -7
      ELSE IF( ( LLERI .AND. LDQCO1.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDQCO1.LT.MAX( 1, MINDEX ) ) ) THEN
         INFO = -9
      ELSE IF( ( LLERI .AND. LDQCO2.LT.MAX( 1, M ) ) .OR.
     $    ( .NOT.LLERI .AND. LDQCO2.LT.MAX( 1, MINDEX ) ) ) THEN
         INFO = -10
      END IF
C
      N = 0
      IF ( INFO.EQ.0 ) THEN
         IF ( LLERI ) THEN
            PWORK = P
            MWORK = M
         ELSE
            PWORK = M
            MWORK = P
         END IF
C
         MAXIND = 0
         DO 10 I = 1, PWORK
            N = N + INDEX(I)
            IF ( INDEX(I).GT.MAXIND ) MAXIND = INDEX(I)
   10    CONTINUE
         KPCOEF = MAXIND + 1
      END IF
C
      IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF( LDC.LT.MAX( 1, MINDEX ) ) THEN
         INFO = -18
      ELSE IF( LDD.LT.MAX( 1, MINDEX ) ) THEN
         INFO = -20
      ELSE IF( LDWORK.LT.MAX( 1, MINDEX*( MINDEX + 4 ) ) ) THEN
         INFO = -23
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TC04AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( M.EQ.0 .OR. P.EQ.0 ) THEN
         N = 0
         RCOND = ONE
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( .NOT.LLERI ) THEN
C
C        Initialization for right matrix fraction: obtain the dual
C        system.
C
         CALL TC01OD( 'R', M, P, KPCOEF, PCOEFF, LDPCO1, LDPCO2,
     $                QCOEFF, LDQCO1, LDQCO2, INFO )
      END IF
C
C     Store leading row coefficient matrix of P(s).
C
      LDW = MAX( 1, PWORK )
      CALL DLACPY( 'Full', PWORK, PWORK, PCOEFF, LDPCO1, DWORK, LDW )
C
C     Check if P(s) is row proper: if not, exit.
C
      DWNORM = DLANGE( '1-norm', PWORK, PWORK, DWORK, LDW, DWORK )
C
      CALL DGETRF( PWORK, PWORK, DWORK, LDW, IWORK, INFO )
C
C     Workspace: need  PWORK*(PWORK + 4).
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      JWORK = LDW*PWORK + 1
C
      CALL DGECON( '1-norm', PWORK, DWORK, LDW, DWNORM, RCOND,
     $             DWORK(JWORK), IWORK(PWORK+1), INFO )
C
      WRKOPT = MAX( 1, PWORK*(PWORK + 4) )
C
      IF ( RCOND.LE.DLAMCH( 'Epsilon' ) ) THEN
C
C        Error return: P(s) is not row proper.
C
         INFO = 1
         RETURN
      ELSE
C
C        Calculate the order of equivalent state-space representation,
C        and initialize A.
C
         CALL DLASET( 'Full', N, N, ZERO, ZERO, A, LDA )
C
         DWORK(JWORK) = ONE
         IF ( N.GT.1 ) CALL DCOPY( N-1, DWORK(JWORK), 0, A(2,1), LDA+1 )
C
C        Find the PWORK ordered 'non-trivial' columns row by row,
C        in PWORK row blocks, the I-th having INDEX(I) rows.
C
         IBIAS = 2
C
         DO 50 I = 1, PWORK
            KSTOP = INDEX(I) + 1
            IF ( KSTOP.NE.1 ) THEN
               IBIAS = IBIAS + INDEX(I)
C
C              These rows given from the lower coefficients of row I
C              of P(s).
C
               DO 40 K = 2, KSTOP
                  IA = IBIAS - K
C
                  DO 20 J = 1, PWORK
                     DWORK(JWORK+J-1) = -PCOEFF(I,J,K)
   20             CONTINUE
C
                  CALL DGETRS( 'Transpose', PWORK, 1, DWORK, LDW,
     $                         IWORK, DWORK(JWORK), LDW, INFO )
C
                  JA = 0
C
                  DO 30 J = 1, PWORK
                     IF ( INDEX(J).NE.0 ) THEN
                        JA = JA + INDEX(J)
                        A(IA,JA) = DWORK(JWORK+J-1)
                     END IF
   30             CONTINUE
C
C                 Also, set up B and C (temporarily) for use when
C                 finding B.
C
                  CALL DCOPY( MWORK, QCOEFF(I,1,K), LDQCO1, B(IA,1),
     $                        LDB )
                  CALL DCOPY( PWORK, PCOEFF(I,1,K), LDPCO1, C(1,IA), 1 )
   40          CONTINUE
C
            END IF
   50    CONTINUE
C
C        Calculate D from the leading coefficients of P and Q.
C
         CALL DLACPY( 'Full', PWORK, MWORK, QCOEFF, LDQCO1, D, LDD )
C
         CALL DGETRS( 'No transpose', PWORK, MWORK, DWORK, LDW, IWORK,
     $                D, LDD, INFO )
C
C        For B and C as set up above, desired B = B - (C' * D).
C
         CALL DGEMM( 'Transpose', 'No transpose', N, MWORK, PWORK, -ONE,
     $               C, LDC, D, LDD, ONE, B, LDB )
C
C        Finally, calculate C: zero, apart from ...
C
         CALL DLASET( 'Full', PWORK, N, ZERO, ZERO, C, LDC )
C
C        PWORK ordered 'non-trivial' columns, equal to those
C        of inv(DWORK).
C
C        Workspace: need   PWORK*(PWORK + 1);
C                   prefer PWORK*PWORK + PWORK*NB.
C
         CALL DGETRI( PWORK, DWORK, LDW, IWORK, DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
         JC = 0
         JW = 1
C
         DO 60 J = 1, PWORK
            IF ( INDEX(J).NE.0 ) THEN
               JC = JC + INDEX(J)
               CALL DCOPY( PWORK, DWORK(JW), 1, C(1,JC), 1 )
            END IF
            JW = JW + LDW
   60    CONTINUE
C
      END IF
C
C     For right matrix fraction, return to original (dual of dual)
C     system.
C
      IF ( .NOT.LLERI ) THEN
         CALL TC01OD( 'L', MWORK, PWORK, KPCOEF, PCOEFF, LDPCO1,
     $                LDPCO2, QCOEFF, LDQCO1, LDQCO2, INFO )
C
C        Also, obtain dual of state-space representation.
C
         CALL AB07MD( 'D', N, MWORK, PWORK, A, LDA, B, LDB, C, LDC, D,
     $                LDD, INFO )
      END IF
C
C     Set optimal workspace dimension.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TC04AD ***
      END
