      SUBROUTINE TD03AD( ROWCOL, LERI, EQUIL, M, P, INDEXD, DCOEFF,
     $                   LDDCOE, UCOEFF, LDUCO1, LDUCO2, NR, A, LDA, B,
     $                   LDB, C, LDC, D, LDD, INDEXP, PCOEFF, LDPCO1,
     $                   LDPCO2, QCOEFF, LDQCO1, LDQCO2, VCOEFF, LDVCO1,
     $                   LDVCO2, TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To find a relatively prime left or right polynomial matrix
C     representation for a proper transfer matrix T(s) given as either
C     row or column polynomial vectors over common denominator
C     polynomials, possibly with uncancelled common terms.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ROWCOL  CHARACTER*1
C             Indicates whether T(s) is to be factorized by rows or by
C             columns as follows:
C             = 'R':  T(s) is factorized by rows;
C             = 'C':  T(s) is factorized by columns.
C
C     LERI    CHARACTER*1
C             Indicates whether a left or a right polynomial matrix
C             representation is required as follows:
C             = 'L':  A left polynomial matrix representation
C                     inv(P(s))*Q(s) is required;
C             = 'R':  A right polynomial matrix representation
C                     Q(s)*inv(P(s)) is required.
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
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     INDEXD  (input) INTEGER array, dimension (P), if ROWCOL = 'R', or
C                                    dimension (M), if ROWCOL = 'C'.
C             The leading pormd elements of this array must contain the
C             row degrees of the denominator polynomials in D(s).
C             pormd = P if the transfer matrix T(s) is given as row
C             polynomial vectors over denominator polynomials;
C             pormd = M if the transfer matrix T(s) is given as column
C             polynomial vectors over denominator polynomials.
C
C     DCOEFF  (input) DOUBLE PRECISION array, dimension (LDDCOE,kdcoef),
C             where kdcoef = MAX(INDEXD(I)) + 1.
C             The leading pormd-by-kdcoef part of this array must
C             contain the coefficients of each denominator polynomial.
C             DCOEFF(I,K) is the coefficient in s**(INDEXD(I)-K+1) of
C             the I-th denominator polynomial in D(s), where K = 1,2,
C             ...,kdcoef.
C
C     LDDCOE  INTEGER
C             The leading dimension of array DCOEFF.
C             LDDCOE >= MAX(1,P), if ROWCOL = 'R';
C             LDDCOE >= MAX(1,M), if ROWCOL = 'C'.
C
C     UCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDUCO1,LDUCO2,kdcoef)
C             The leading P-by-M-by-kdcoef part of this array must
C             contain the coefficients of the numerator matrix U(s);
C             if ROWCOL = 'C', this array is modified internally but
C             restored on exit, and the remainder of the leading
C             MAX(M,P)-by-MAX(M,P)-by-kdcoef part is used as internal
C             workspace.
C             UCOEFF(I,J,K) is the coefficient in s**(INDEXD(iorj)-K+1)
C             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef;
C             iorj = I if T(s) is given as row polynomial vectors over
C             denominator polynomials; iorj = J if T(s) is given as
C             column polynomial vectors over denominator polynomials.
C             Thus for ROWCOL = 'R', U(s) =
C             diag(s**INDEXD(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...).
C
C     LDUCO1  INTEGER
C             The leading dimension of array UCOEFF.
C             LDUCO1 >= MAX(1,P),   if ROWCOL = 'R';
C             LDUCO1 >= MAX(1,M,P), if ROWCOL = 'C'.
C
C     LDUCO2  INTEGER
C             The second dimension of array UCOEFF.
C             LDUCO2 >= MAX(1,M),   if ROWCOL = 'R';
C             LDUCO2 >= MAX(1,M,P), if ROWCOL = 'C'.
C
C     NR      (output) INTEGER
C             The order of the resulting minimal realization, i.e. the
C             order of the state dynamics matrix A.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N),
C                      pormd
C             where N = SUM INDEXD(I)
C                       I=1
C             The leading NR-by-NR part of this array contains the upper
C             block Hessenberg state dynamics matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P))
C             The leading NR-by-M part of this array contains the
C             input/state matrix B; the remainder of the leading
C             N-by-MAX(M,P) part is used as internal workspace.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-NR part of this array contains the
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
C     INDEXP  (output) INTEGER array, dimension (P), if ROWCOL = 'R', or
C                                     dimension (M), if ROWCOL = 'C'.
C             The leading pormp elements of this array contain the
C             row (column if ROWCOL = 'C') degrees of the denominator
C             matrix P(s).
C             pormp = P if a left polynomial matrix representation
C             is requested; pormp = M if a right polynomial matrix
C             representation is requested.
C             These elements are ordered so that
C             INDEXP(1) >= INDEXP(2) >= ... >= INDEXP(pormp).
C
C     PCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDPCO1,LDPCO2,N+1)
C             The leading pormp-by-pormp-by-kpcoef part of this array
C             contains the coefficients of the denominator matrix P(s),
C             where kpcoef = MAX(INDEXP(I)) + 1.
C             PCOEFF(I,J,K) is the coefficient in s**(INDEXP(iorj)-K+1)
C             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef;
C             iorj = I if a left polynomial matrix representation is
C             requested; iorj = J if a right polynomial matrix
C             representation is requested.
C             Thus for a left polynomial matrix representation, P(s) =
C             diag(s**INDEXP(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...).
C
C     LDPCO1  INTEGER
C             The leading dimension of array PCOEFF.
C             LDPCO1 >= MAX(1,P), if ROWCOL = 'R';
C             LDPCO1 >= MAX(1,M), if ROWCOL = 'C'.
C
C     LDPCO2  INTEGER
C             The second dimension of array PCOEFF.
C             LDPCO2 >= MAX(1,P), if ROWCOL = 'R';
C             LDPCO2 >= MAX(1,M), if ROWCOL = 'C'.
C
C     QCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDQCO1,LDQCO2,N+1)
C             The leading pormp-by-pormd-by-kpcoef part of this array
C             contains the coefficients of the numerator matrix Q(s).
C             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K).
C
C     LDQCO1  INTEGER
C             The leading dimension of array QCOEFF.
C             If LERI = 'L', LDQCO1 >= MAX(1,PM),
C                                      where PM = P, if ROWCOL = 'R';
C                                            PM = M, if ROWCOL = 'C'.
C             If LERI = 'R', LDQCO1 >= MAX(1,M,P).
C
C     LDQCO2  INTEGER
C             The second dimension of array QCOEFF.
C             If LERI = 'L', LDQCO2 >= MAX(1,MP),
C                                      where MP = M, if ROWCOL = 'R';
C                                            MP = P, if ROWCOL = 'C'.
C             If LERI = 'R', LDQCO2 >= MAX(1,M,P).
C
C     VCOEFF  (output) DOUBLE PRECISION array, dimension
C             (LDVCO1,LDVCO2,N+1)
C             The leading pormp-by-NR-by-kpcoef part of this array
C             contains the coefficients of the intermediate matrix
C             V(s) as produced by SLICOT Library routine TB03AD.
C
C     LDVCO1  INTEGER
C             The leading dimension of array VCOEFF.
C             LDVCO1 >= MAX(1,P), if ROWCOL = 'R';
C             LDVCO1 >= MAX(1,M), if ROWCOL = 'C'.
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
C             where  PM = P, if ROWCOL = 'R';
C                    PM = M, if ROWCOL = 'C'.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i (i <= k = pormd), then i is the first
C                   integer I for which ABS( DCOEFF(I,1) ) is so small
C                   that the calculations would overflow (see SLICOT
C                   Library routine TD03AY); that is, the leading
C                   coefficient of a polynomial is nearly zero; no
C                   state-space representation or polynomial matrix
C                   representation is calculated;
C             = k+1:  if a singular matrix was encountered during the
C                   computation of V(s);
C             = k+2:  if a singular matrix was encountered during the
C                   computation of P(s).
C
C     METHOD
C
C     The method for transfer matrices factorized by rows will be
C     described here; T(s) factorized by columns is dealt with by
C     operating on the dual T'(s). The description for T(s) is actually
C     the left polynomial matrix representation
C
C          T(s) = inv(D(s))*U(s),
C
C     where D(s) is diagonal with its (I,I)-th polynomial element of
C     degree INDEXD(I). The first step is to check whether the leading
C     coefficient of any polynomial element of D(s) is approximately
C     zero, if so the routine returns with INFO > 0. Otherwise,
C     Wolovich's Observable Structure Theorem is used to construct a
C     state-space representation in observable companion form which is
C     equivalent to the above polynomial matrix representation. The
C     method is particularly easy here due to the diagonal form of D(s).
C     This state-space representation is not necessarily controllable
C     (as D(s) and U(s) are not necessarily relatively left prime), but
C     it is in theory completely observable; however, its observability
C     matrix may be poorly conditioned, so it is treated as a general
C     state-space representation and SLICOT Library routine TB03AD is
C     used to separate out a minimal realization for T(s) from it by
C     means of orthogonal similarity transformations and then to
C     calculate a relatively prime (left or right) polynomial matrix
C     representation which is equivalent to this.
C
C     REFERENCES
C
C     [1] Patel, R.V.
C         On Computing Matrix Fraction Descriptions and Canonical
C         Forms of Linear Time-Invariant Systems.
C         UMIST Control Systems Centre Report 489, 1980.
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
C     V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1998.
C     Supersedes Release 3.0 routine TD01ND.
C
C     REVISIONS
C
C     -
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
      CHARACTER         EQUIL, LERI, ROWCOL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDDCOE, LDPCO1,
     $                  LDPCO2, LDQCO1, LDQCO2, LDUCO1, LDUCO2, LDVCO1,
     $                  LDVCO2, LDWORK, M, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INDEXD(*), INDEXP(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DCOEFF(LDDCOE,*), DWORK(*),
     $                  PCOEFF(LDPCO1,LDPCO2,*),
     $                  QCOEFF(LDQCO1,LDQCO2,*),
     $                  UCOEFF(LDUCO1,LDUCO2,*), VCOEFF(LDVCO1,LDVCO2,*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LLERI, LROWCO
      INTEGER           I, IDUAL, ITEMP, J, JSTOP, K, KDCOEF, KPCOEF,
     $                  MAXMP, MPLIM, MWORK, N, PWORK
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB07MD, DLACPY, DSWAP, TB01XD, TB03AD, TC01OD,
     $                  TD03AY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO = 0
      LROWCO = LSAME( ROWCOL, 'R' )
      LLERI  = LSAME( LERI,   'L' )
      LEQUIL = LSAME( EQUIL,  'S' )
C
C     Test the input scalar arguments.
C
      MAXMP  = MAX( M, P )
      MPLIM  = MAX( 1, MAXMP )
      IF ( LROWCO ) THEN
C
C        Initialization for T(s) given as rows over common denominators.
C
         PWORK = P
         MWORK = M
      ELSE
C
C        Initialization for T(s) given as columns over common
C        denominators.
C
         PWORK = M
         MWORK = P
      END IF
C
      IF( .NOT.LROWCO .AND. .NOT.LSAME( ROWCOL, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LLERI .AND. .NOT.LSAME( LERI, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDDCOE.LT.MAX( 1, PWORK ) ) THEN
         INFO = -8
      ELSE IF( LDUCO1.LT.MAX( 1, PWORK ) .OR. ( .NOT.LROWCO .AND.
     $         LDUCO1.LT.MPLIM ) ) THEN
         INFO = -10
      ELSE IF( LDUCO2.LT.MAX( 1, MWORK ) .OR. ( .NOT.LROWCO .AND.
     $         LDUCO2.LT.MPLIM ) ) THEN
         INFO = -11
      END IF
C
      N = 0
      IF ( INFO.EQ.0 ) THEN
C
C        Calculate N, the order of the resulting state-space
C        representation, and the index kdcoef.
C
         KDCOEF = 0
C
         DO 10 I = 1, PWORK
            KDCOEF = MAX( KDCOEF, INDEXD(I) )
            N = N + INDEXD(I)
   10    CONTINUE
C
         KDCOEF = KDCOEF + 1
C
         IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -16
         ELSE IF( LDC.LT.MPLIM ) THEN
            INFO = -18
         ELSE IF( LDD.LT.MPLIM ) THEN
            INFO = -20
         ELSE IF( LDPCO1.LT.PWORK ) THEN
            INFO = -23
         ELSE IF( LDPCO2.LT.PWORK ) THEN
            INFO = -24
         ELSE IF( LDQCO1.LT.MAX( 1, PWORK ) .OR. ( .NOT.LLERI .AND.
     $            LDQCO1.LT.MPLIM ) ) THEN
            INFO = -26
         ELSE IF( LDQCO2.LT.MAX( 1, MWORK ) .OR. ( .NOT.LLERI .AND.
     $            LDQCO2.LT.MPLIM ) ) THEN
            INFO = -27
         ELSE IF( LDVCO1.LT.MAX( 1, PWORK ) ) THEN
            INFO = -29
         ELSE IF( LDVCO2.LT.MAX( 1, N ) ) THEN
            INFO = -30
C
         ELSE IF( LDWORK.LT.MAX( 1, N + MAX( N, 3*MAXMP ),
     $                              PWORK*( PWORK + 2 ) ) ) THEN
            INFO = -34
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TD03AD', -INFO )
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
C     IDUAL = 1 iff precisely ROWCOL = 'C' or (exclusively) LERI = 'R',
C     i.e. iff AB07MD call is required before TB03AD.
C
      IDUAL = 0
      IF ( .NOT.LROWCO ) IDUAL = 1
      IF ( .NOT.LLERI  ) IDUAL = IDUAL + 1
C
      IF ( .NOT.LROWCO ) THEN
C
C        Initialize the remainder of the leading
C        MPLIM-by-MPLIM-by-KDCOEF part of U(s) to zero.
C
         IF ( P.LT.M ) THEN
C
            DO 20 K = 1, KDCOEF
               CALL DLACPY( 'Full', M-P, MPLIM, ZERO, ZERO,
     $                      UCOEFF(P+1,1,K), LDUCO1 )
   20       CONTINUE
C
         ELSE IF ( P.GT.M ) THEN
C
            DO 30 K = 1, KDCOEF
               CALL DLACPY( 'Full', MPLIM, P-M, ZERO, ZERO,
     $                      UCOEFF(1,M+1,K), LDUCO1 )
   30       CONTINUE
C
         END IF
C
         IF ( MPLIM.NE.1 ) THEN
C
C           Non-scalar T(s) factorized by columns: transpose it
C           (i.e. U(s)).
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
C     appropriate,
C
      CALL TD03AY( MWORK, PWORK, INDEXD, DCOEFF, LDDCOE, UCOEFF, LDUCO1,
     $             LDUCO2, N, A, LDA, B, LDB, C, LDC, D, LDD, INFO )
      IF ( INFO.GT.0 )
     $   RETURN
C
      IF ( IDUAL.EQ.1 ) THEN
C
C        and then obtain (MWORK x PWORK) dual of this system if
C        appropriate.
C
         CALL AB07MD( 'D', N, MWORK, PWORK, A, LDA, B, LDB, C, LDC, D,
     $                LDD, INFO )
         ITEMP = PWORK
         PWORK = MWORK
         MWORK = ITEMP
      END IF
C
C     Find left polynomial matrix representation (and minimal
C     state-space representation en route) for the relevant state-space
C     representation ...
C
      CALL TB03AD( 'Left', EQUIL, N, MWORK, PWORK, A, LDA, B, LDB, C,
     $             LDC, D, LDD, NR, INDEXP, PCOEFF, LDPCO1, LDPCO2,
     $             QCOEFF, LDQCO1, LDQCO2, VCOEFF, LDVCO1, LDVCO2, TOL,
     $             IWORK, DWORK, LDWORK, INFO )
C
      IF ( INFO.GT.0 ) THEN
         INFO = PWORK + INFO
         RETURN
      END IF
C
      IF ( .NOT.LLERI ) THEN
C
C        and, if a right polynomial matrix representation is required,
C        transpose and reorder (to get a block upper Hessenberg
C        matrix A).
C
         K = IWORK(1) - 1
         IF ( N.GE.2 )
     $      K = K + IWORK(2)
         CALL TB01XD( 'D', NR, MWORK, PWORK, K, NR-1, A, LDA, B, LDB, C,
     $                LDC, D, LDD, INFO )
C
         KPCOEF = 0
C
         DO 60 I = 1, PWORK
            KPCOEF = MAX( KPCOEF, INDEXP(I) )
   60    CONTINUE
C
         KPCOEF = KPCOEF + 1
         CALL TC01OD( 'L', MWORK, PWORK, KPCOEF, PCOEFF, LDPCO1, LDPCO2,
     $                QCOEFF, LDQCO1, LDQCO2, INFO )
      END IF
C
      IF ( ( .NOT.LROWCO ) .AND. ( MPLIM.NE.1 ) ) THEN
C
C        If non-scalar T(s) originally given by columns,
C        retranspose U(s).
C
         DO 80 K = 1, KDCOEF
C
            DO 70 J = 1, JSTOP
               CALL DSWAP( MPLIM-J, UCOEFF(J+1,J,K), 1, UCOEFF(J,J+1,K),
     $                     LDUCO1 )
   70       CONTINUE
C
   80    CONTINUE
C
      END IF
      RETURN
C *** Last line of TD03AD ***
      END
