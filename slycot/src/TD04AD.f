      SUBROUTINE TD04AD( ROWCOL, M, P, INDEX, DCOEFF, LDDCOE, UCOEFF,
     $                   LDUCO1, LDUCO2, NR, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To find a minimal state-space representation (A,B,C,D) for a
C     proper transfer matrix T(s) given as either row or column
C     polynomial vectors over denominator polynomials, possibly with
C     uncancelled common terms.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ROWCOL  CHARACTER*1
C             Indicates whether the transfer matrix T(s) is given as
C             rows or columns over common denominators as follows:
C             = 'R':  T(s) is given as rows over common denominators;
C             = 'C':  T(s) is given as columns over common denominators.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     INDEX   (input) INTEGER array, dimension (porm), where porm = P,
C             if ROWCOL = 'R', and porm = M, if ROWCOL = 'C'.
C             This array must contain the degrees of the denominator
C             polynomials in D(s).
C
C     DCOEFF  (input) DOUBLE PRECISION array, dimension (LDDCOE,kdcoef),
C             where kdcoef = MAX(INDEX(I)) + 1.
C             The leading porm-by-kdcoef part of this array must contain
C             the coefficients of each denominator polynomial.
C             DCOEFF(I,K) is the coefficient in s**(INDEX(I)-K+1) of the
C             I-th denominator polynomial in D(s), where
C             K = 1,2,...,kdcoef.
C
C     LDDCOE  INTEGER
C             The leading dimension of array DCOEFF.
C             LDDCOE >= MAX(1,P) if ROWCOL = 'R';
C             LDDCOE >= MAX(1,M) if ROWCOL = 'C'.
C
C     UCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDUCO1,LDUCO2,kdcoef)
C             The leading P-by-M-by-kdcoef part of this array must
C             contain the numerator matrix U(s); if ROWCOL = 'C', this
C             array is modified internally but restored on exit, and the
C             remainder of the leading MAX(M,P)-by-MAX(M,P)-by-kdcoef
C             part is used as internal workspace.
C             UCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1)
C             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef;
C             if ROWCOL = 'R' then iorj = I, otherwise iorj = J.
C             Thus for ROWCOL = 'R', U(s) =
C             diag(s**INDEX(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...).
C
C     LDUCO1  INTEGER
C             The leading dimension of array UCOEFF.
C             LDUCO1 >= MAX(1,P)   if ROWCOL = 'R';
C             LDUCO1 >= MAX(1,M,P) if ROWCOL = 'C'.
C
C     LDUCO2  INTEGER
C             The second dimension of array UCOEFF.
C             LDUCO2 >= MAX(1,M)   if ROWCOL = 'R';
C             LDUCO2 >= MAX(1,M,P) if ROWCOL = 'C'.
C
C     NR      (output) INTEGER
C             The order of the resulting minimal realization, i.e. the
C             order of the state dynamics matrix A.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N),
C                       porm
C             where N = SUM INDEX(I).
C                       I=1
C             The leading NR-by-NR part of this array contains the upper
C             block Hessenberg state dynamics matrix A of a minimal
C             realization.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P))
C             The leading NR-by-M part of this array contains the
C             input/state matrix B of a minimal realization; the
C             remainder of the leading N-by-MAX(M,P) part is used as
C             internal workspace.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-NR part of this array contains the
C             state/output matrix C of a minimal realization; the
C             remainder of the leading MAX(M,P)-by-N part is used as
C             internal workspace.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M,P).
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M),
C             if ROWCOL = 'R', and (LDD,MAX(M,P)) if ROWCOL = 'C'.
C             The leading P-by-M part of this array contains the direct
C             transmission matrix D; if ROWCOL = 'C', the remainder of
C             the leading MAX(M,P)-by-MAX(M,P) part is used as internal
C             workspace.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P)   if ROWCOL = 'R';
C             LDD >= MAX(1,M,P) if ROWCOL = 'C'.
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
C                   value;
C             > 0:  if INFO = i, then i is the first integer for which
C                   ABS( DCOEFF(I,1) ) is so small that the calculations
C                   would overflow (see SLICOT Library routine TD03AY);
C                   that is, the leading coefficient of a polynomial is
C                   nearly zero; no state-space representation is
C                   calculated.
C
C     METHOD
C
C     The method for transfer matrices factorized by rows will be
C     described here: T(s) factorized by columns is dealt with by
C     operating on the dual T'(s). This description for T(s) is
C     actually the left polynomial matrix representation
C
C          T(s) = inv(D(s))*U(s),
C
C     where D(s) is diagonal with its (I,I)-th polynomial element of
C     degree INDEX(I). The first step is to check whether the leading
C     coefficient of any polynomial element of D(s) is approximately
C     zero; if so the routine returns with INFO > 0. Otherwise,
C     Wolovich's Observable Structure Theorem is used to construct a
C     state-space representation in observable companion form which
C     is equivalent to the above polynomial matrix representation.
C     The method is particularly easy here due to the diagonal form
C     of D(s). This state-space representation is not necessarily
C     controllable (as D(s) and U(s) are not necessarily relatively
C     left prime), but it is in theory completely observable; however,
C     its observability matrix may be poorly conditioned, so it is
C     treated as a general state-space representation and SLICOT
C     Library routine TB01PD is then called to separate out a minimal
C     realization from this general state-space representation by means
C     of orthogonal similarity transformations.
C
C     REFERENCES
C
C     [1] Patel, R.V.
C         Computation of Minimal-Order State-Space Realizations and
C         Observability Indices using Orthogonal Transformations.
C         Int. J. Control, 33, pp. 227-246, 1981.
C
C     [2] Wolovich, W.A.
C         Linear Multivariable Systems, (Theorem 4.3.3).
C         Springer-Verlag, 1974.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998.
C     Supersedes Release 3.0 routine TD01OD.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Controllability, elementary polynomial operations, minimal
C     realization, polynomial matrix, state-space representation,
C     transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         ROWCOL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDDCOE, LDUCO1,
     $                  LDUCO2, LDWORK, M, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DCOEFF(LDDCOE,*), DWORK(*),
     $                  UCOEFF(LDUCO1,LDUCO2,*)
C     .. Local Scalars ..
      LOGICAL           LROCOC, LROCOR
      INTEGER           I, J, JSTOP, K, KDCOEF, MPLIM, MWORK, N, PWORK,
     $                  KU
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLASET, DSWAP, TB01PD, TB01XD, TD03AY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO = 0
      LROCOR = LSAME( ROWCOL, 'R' )
      LROCOC = LSAME( ROWCOL, 'C' )
      MPLIM = MAX( 1, M, P )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LROCOR .AND. .NOT.LROCOC ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( LROCOR .AND. LDDCOE.LT.MAX( 1, P ) ) .OR.
     $         ( LROCOC .AND. LDDCOE.LT.MAX( 1, M ) ) ) THEN
         INFO = -6
      ELSE IF( ( LROCOR .AND. LDUCO1.LT.MAX( 1, P ) ) .OR.
     $         ( LROCOC .AND. LDUCO1.LT.MPLIM ) ) THEN
         INFO = -8
      ELSE IF( ( LROCOR .AND. LDUCO2.LT.MAX( 1, M ) ) .OR.
     $         ( LROCOC .AND. LDUCO2.LT.MPLIM ) ) THEN
         INFO = -9
      END IF
C
      N = 0
      IF ( INFO.EQ.0 ) THEN
         IF ( LROCOR ) THEN
C
C           Initialization for T(s) given as rows over common
C           denominators.
C
            PWORK = P
            MWORK = M
         ELSE
C
C           Initialization for T(s) given as columns over common
C           denominators.
C
            PWORK = M
            MWORK = P
         END IF
C
C        Calculate N, the order of the resulting state-space
C        representation.
C
         KDCOEF = 0
C
         DO 10 I = 1, PWORK
            KDCOEF = MAX( KDCOEF, INDEX(I) )
            N = N + INDEX(I)
   10    CONTINUE
C
         KDCOEF = KDCOEF + 1
C
         IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -12
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDC.LT.MPLIM ) THEN
            INFO = -16
         ELSE IF( ( LROCOR .AND. LDD.LT.MAX( 1, P ) ) .OR.
     $         ( LROCOC .AND. LDD.LT.MPLIM ) ) THEN
            INFO = -18
         ELSE IF( LDWORK.LT.MAX( 1, N + MAX( N, 3*M, 3*P ) ) ) THEN
            INFO = -22
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TD04AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, P ).EQ.0 ) THEN
         NR  = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( LROCOC ) THEN
C
C        Initialize the remainder of the leading
C        MPLIM-by-MPLIM-by-KDCOEF part of U(s) to zero.
C
         IF ( P.LT.M ) THEN
C
            DO 20 K = 1, KDCOEF
               CALL DLASET( 'Full', M-P, MPLIM, ZERO, ZERO,
     $                      UCOEFF(P+1,1,K), LDUCO1 )
   20       CONTINUE
C
         ELSE IF ( P.GT.M ) THEN
C
            DO 30 K = 1, KDCOEF
               CALL DLASET( 'Full', MPLIM, P-M, ZERO, ZERO,
     $                      UCOEFF(1,M+1,K), LDUCO1 )
   30       CONTINUE
C
         END IF
C
         IF ( MPLIM.NE.1 ) THEN
C
C           Non-scalar T(s) factorized by columns: transpose it (i.e.
C           U(s)).
C
            JSTOP = MPLIM - 1
C
            DO 50 K = 1, KDCOEF
C
               DO 40 J = 1, JSTOP
                  CALL DSWAP( MPLIM-J, UCOEFF(J+1,J,K), 1,
     $                        UCOEFF(J,J+1,K), LDUCO1 )
   40          CONTINUE
C
   50       CONTINUE
C
         END IF
      END IF
C
C     Construct non-minimal state-space representation (by Wolovich's
C     Structure Theorem) which has transfer matrix T(s) or T'(s) as
C     appropriate ...
C
      CALL TD03AY( MWORK, PWORK, INDEX, DCOEFF, LDDCOE, UCOEFF, LDUCO1,
     $             LDUCO2, N, A, LDA, B, LDB, C, LDC, D, LDD, INFO )
      IF ( INFO.GT.0 )
     $   RETURN
C
C     and then separate out a minimal realization from this.
C
C     Workspace: need  N + MAX(N, 3*MWORK, 3*PWORK).
C
      CALL TB01PD( 'Minimal', 'Scale', N, MWORK, PWORK, A, LDA, B, LDB,
     $             C, LDC, NR, TOL, IWORK, DWORK, LDWORK, INFO )
C
      IF ( LROCOC ) THEN
C
C        If T(s) originally factorized by columns, find dual of minimal
C        state-space representation, and reorder the rows and columns
C        to get an upper block Hessenberg state dynamics matrix.
C
C        IWORK contains the orders of the diagnonal blocks
C        RvP, In TB01PD, IWORK is zeroed from INDCON to N, beyond N it may
C        contain nonsense?         
         K = -1 
         DO 55 I = 1, N
            K = K + IWORK(I)
 55      CONTINUE
C
C        RvP 180615 Try to protect against re-working an empty [] A
C        matrix, failed with K < 0
C                  
         CALL TB01XD( 'D', NR, MWORK, PWORK, MAX(0, K), MAX(0,NR-1),
     $                A, LDA, B, LDB, C, LDC, D, LDD, INFO )
         IF ( MPLIM.NE.1 ) THEN
C
C           Also, retranspose U(s) if this is non-scalar.
C
            DO 70 K = 1, KDCOEF
C 
               DO 60 J = 1, JSTOP
                  CALL DSWAP( MPLIM-J, UCOEFF(J+1,J,K), 1,
     $                           UCOEFF(J,J+1,K), LDUCO1 )
   60          CONTINUE
C
   70       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of TD04AD ***
      END
