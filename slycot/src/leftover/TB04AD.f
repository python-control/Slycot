      SUBROUTINE TB04AD( ROWCOL, N, M, P, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, NR, INDEX, DCOEFF, LDDCOE, UCOEFF, LDUCO1,
     $                   LDUCO2, TOL1, TOL2, IWORK, DWORK, LDWORK,
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
C     To find the transfer matrix T(s) of a given state-space
C     representation (A,B,C,D). T(s) is expressed as either row or
C     column polynomial vectors over monic least common denominator
C     polynomials.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ROWCOL  CHARACTER*1
C             Indicates whether the transfer matrix T(s) is required
C             as rows or columns over common denominators as follows:
C             = 'R':  T(s) is required as rows over common denominators;
C             = 'C':  T(s) is required as columns over common
C                     denominators.
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
C             the upper block Hessenberg state dynamics matrix A of a
C             transformed representation for the original system: this
C             is completely controllable if ROWCOL = 'R', or completely
C             observable if ROWCOL = 'C'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M),
C             if ROWCOL = 'R', and (LDB,MAX(M,P)) if ROWCOL = 'C'.
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B; if
C             ROWCOL = 'C', the remainder of the leading N-by-MAX(M,P)
C             part is used as internal workspace.
C             On exit, the leading NR-by-M part of this array contains
C             the transformed input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C; if
C             ROWCOL = 'C', the remainder of the leading MAX(M,P)-by-N
C             part is used as internal workspace.
C             On exit, the leading P-by-NR part of this array contains
C             the transformed state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P)   if ROWCOL = 'R';
C             LDC >= MAX(1,M,P) if ROWCOL = 'C'.
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M),
C             if ROWCOL = 'R', and (LDD,MAX(M,P)) if ROWCOL = 'C'.
C             The leading P-by-M part of this array must contain the
C             original direct transmission matrix D; if ROWCOL = 'C',
C             this array is modified internally, but restored on exit,
C             and the remainder of the leading MAX(M,P)-by-MAX(M,P)
C             part is used as internal workspace.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P)   if ROWCOL = 'R';
C             LDD >= MAX(1,M,P) if ROWCOL = 'C'.
C
C     NR      (output) INTEGER
C             The order of the transformed state-space representation.
C
C     INDEX   (output) INTEGER array, dimension (porm), where porm = P,
C             if ROWCOL = 'R', and porm = M, if ROWCOL = 'C'.
C             The degrees of the denominator polynomials.
C
C     DCOEFF  (output) DOUBLE PRECISION array, dimension (LDDCOE,N+1)
C             The leading porm-by-kdcoef part of this array contains
C             the coefficients of each denominator polynomial, where
C             kdcoef = MAX(INDEX(I)) + 1.
C             DCOEFF(I,K) is the coefficient in s**(INDEX(I)-K+1) of
C             the I-th denominator polynomial, where K = 1,2,...,kdcoef.
C
C     LDDCOE  INTEGER
C             The leading dimension of array DCOEFF.
C             LDDCOE >= MAX(1,P) if ROWCOL = 'R';
C             LDDCOE >= MAX(1,M) if ROWCOL = 'C'.
C
C     UCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDUCO1,LDUCO2,N+1)
C             If ROWCOL = 'R' then porp = M, otherwise porp = P.
C             The leading porm-by-porp-by-kdcoef part of this array
C             contains the coefficients of the numerator matrix U(s).
C             UCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1)
C             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef;
C             if ROWCOL = 'R' then iorj = I, otherwise iorj = J.
C             Thus for ROWCOL = 'R', U(s) =
C             diag(s**INDEX(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...).
C
C     LDUCO1  INTEGER
C             The leading dimension of array UCOEFF.
C             LDUCO1 >= MAX(1,P) if ROWCOL = 'R';
C             LDUCO1 >= MAX(1,M) if ROWCOL = 'C'.
C
C     LDUCO2  INTEGER
C             The second dimension of array UCOEFF.
C             LDUCO2 >= MAX(1,M) if ROWCOL = 'R';
C             LDUCO2 >= MAX(1,P) if ROWCOL = 'C'.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             The tolerance to be used in determining the i-th row of
C             T(s), where i = 1,2,...,porm. If the user sets TOL1 > 0,
C             then the given value of TOL1 is used as an absolute
C             tolerance; elements with absolute value less than TOL1 are
C             considered neglijible. If the user sets TOL1 <= 0, then
C             an implicitly computed, default tolerance, defined in
C             the SLICOT Library routine TB01ZD, is used instead.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance to be used to separate out a controllable
C             subsystem of (A,B,C). If the user sets TOL2 > 0, then
C             the given value of TOL2 is used as a lower bound for the
C             reciprocal condition number (see the description of the
C             argument RCOND in the SLICOT routine MB03OD);  a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL2 is considered to be of full rank.  If the user sets
C             TOL2 <= 0, then an implicitly computed, default tolerance,
C             defined in the SLICOT Library routine TB01UD, is used
C             instead.
C
C     Workspace
C
C     IWORK   DOUBLE PRECISION array, dimension (N+MAX(M,P))
C             On exit, if INFO = 0, the first nonzero elements of
C             IWORK(1:N) return the orders of the diagonal blocks of A.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N*(N + 1) + MAX(N*MP + 2*N + MAX(N,MP),
C                                       3*MP, PM)),
C             where MP = M, PM = P, if ROWCOL = 'R';
C                   MP = P, PM = M, if ROWCOL = 'C'.
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
C     The method for transfer matrices factorized by rows will be
C     described here: T(s) factorized by columns is dealt with by
C     operating on the dual of the original system.  Each row of
C     T(s) is simply a single-output relatively left prime polynomial
C     matrix representation, so can be calculated by applying a
C     simplified version of the Orthogonal Structure Theorem to a
C     minimal state-space representation for the corresponding row of
C     the given system. A minimal state-space representation is obtained
C     using the Orthogonal Canonical Form to first separate out a
C     completely controllable one for the overall system and then, for
C     each row in turn, applying it again to the resulting dual SIMO
C     (single-input multi-output) system. Note that the elements of the
C     transformed matrix A so calculated are individually scaled in a
C     way which guarantees a monic denominator polynomial.
C
C     REFERENCES
C
C     [1] Williams, T.W.C.
C         An Orthogonal Structure Theorem for Linear Systems.
C         Control Systems Research Group, Kingston Polytechnic,
C         Internal Report 82/2, 1982.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998.
C     Supersedes Release 3.0 routine TB01QD.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Controllability, dual system, minimal realization, orthogonal
C     canonical form, orthogonal transformation, polynomial matrix,
C     transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         ROWCOL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDDCOE, LDUCO1,
     $                  LDUCO2, LDWORK, M, N, NR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DCOEFF(LDDCOE,*), DWORK(*),
     $                  UCOEFF(LDUCO1,LDUCO2,*)
C     .. Local Scalars ..
      LOGICAL           LROCOC, LROCOR
      CHARACTER*1       JOBD
      INTEGER           I, IA, ITAU, J, JWORK, K, KDCOEF, MAXMP, MAXMPN,
     $                  MPLIM, MWORK, N1, PWORK
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07MD, DLASET, DSWAP, TB01XD, TB04AY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C     .. Executable Statements ..
C
      INFO = 0
      LROCOR = LSAME( ROWCOL, 'R' )
      LROCOC = LSAME( ROWCOL, 'C' )
      MAXMP  = MAX( M, P )
      MPLIM  = MAX( 1, MAXMP )
      MAXMPN = MAX( MAXMP, N )
      N1 = MAX( 1, N )
      IF ( LROCOR ) THEN
C
C        T(s) given as rows over common denominators.
C
         PWORK = P
         MWORK = M
      ELSE
C
C        T(s) given as columns over common denominators.
C
         PWORK = M
         MWORK = P
      END IF
C
C     Test the input scalar arguments.
C
      IF( .NOT.LROCOR .AND. .NOT.LROCOC ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -8
      ELSE IF( ( LROCOC .AND. LDC.LT.MPLIM )
     $                   .OR. LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( ( LROCOC .AND. LDD.LT.MPLIM )
     $                   .OR. LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDDCOE.LT.MAX( 1, PWORK ) ) THEN
         INFO = -16
      ELSE IF( LDUCO1.LT.MAX( 1, PWORK ) ) THEN
         INFO = -18
      ELSE IF( LDUCO2.LT.MAX( 1, MWORK ) ) THEN
         INFO = -19
      ELSE IF( LDWORK.LT.MAX( 1, N*( N + 1 ) +
     $                           MAX( N*MWORK + 2*N + MAX( N, MWORK ),
     $                                3*MWORK, PWORK ) ) ) THEN
         INFO = -24
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB04AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAXMPN.EQ.0 )
     $   RETURN
C
      JOBD = 'D'
      IA = 1
      ITAU  = IA + N*N
      JWORK = ITAU + N
C
      IF ( LROCOC ) THEN
C
C        Initialization for T(s) given as columns over common
C        denominators.
C
         CALL AB07MD( JOBD, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                INFO )
      END IF
C
C     Initialize polynomial matrix U(s) to zero.
C
      DO 10 K = 1, N + 1
         CALL DLASET( 'Full', PWORK, MWORK, ZERO, ZERO, UCOEFF(1,1,K),
     $                LDUCO1 )
   10 CONTINUE
C
C     Calculate T(s) by applying the Orthogonal Structure Theorem to
C     each of the PWORK MISO subsystems (A,B,C:I,D:I) in turn.
C
      CALL TB04AY( N, MWORK, PWORK, A, LDA, B, LDB, C, LDC, D, LDD,
     $             NR, INDEX, DCOEFF, LDDCOE, UCOEFF, LDUCO1, LDUCO2,
     $             DWORK(IA), N1, DWORK(ITAU), TOL1, TOL2, IWORK,
     $             DWORK(JWORK), LDWORK-JWORK+1, INFO )
      DWORK(1) = DWORK(JWORK) + DBLE( JWORK-1 )
C
      IF ( LROCOC ) THEN
C
C        For T(s) factorized by columns, return to original (dual of
C        dual) system, and reorder the rows and columns to get an upper
C        block Hessenberg state dynamics matrix.
C
         CALL TB01XD( JOBD, N, MWORK, PWORK, IWORK(1)+IWORK(2)-1, N-1,
     $                A, LDA, B, LDB, C, LDC, D, LDD, INFO )
C
         IF ( MPLIM.NE.1 ) THEN
C
C           Also, transpose U(s) (not 1-by-1).
C
            KDCOEF = 0
C
            DO 20 I = 1, PWORK
               KDCOEF = MAX( KDCOEF, INDEX(I) )
   20       CONTINUE
C
            KDCOEF = KDCOEF + 1
C
            DO 50 K = 1, KDCOEF
C
               DO 40 J = 1, MPLIM - 1
                  CALL DSWAP( MPLIM-J, UCOEFF(J+1,J,K), 1,
     $                        UCOEFF(J,J+1,K), LDUCO1 )
   40          CONTINUE
C
   50       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of TB04AD ***
      END
