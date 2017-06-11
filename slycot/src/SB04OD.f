      SUBROUTINE SB04OD( REDUCE, TRANS, JOBD, M, N, A, LDA, B, LDB, C,
     $                   LDC, D, LDD, E, LDE, F, LDF, SCALE, DIF, P,
     $                   LDP, Q, LDQ, U, LDU, V, LDV, IWORK, DWORK,
     $                   LDWORK, INFO )
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
C     To solve for R and L one of the generalized Sylvester equations
C
C        A * R - L * B = scale * C )
C                                  )                                 (1)
C        D * R - L * E = scale * F )
C
C     or
C
C        A' * R + D' * L = scale * C    )
C                                       )                            (2)
C        R * B' + L * E' = scale * (-F) )
C
C     where A and D are M-by-M matrices, B and E are N-by-N matrices and
C     C, F, R and L are M-by-N matrices.
C
C     The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an
C     output scaling factor chosen to avoid overflow.
C
C     The routine also optionally computes a Dif estimate, which
C     measures the separation of the spectrum of the matrix pair (A,D)
C     from the spectrum of the matrix pair (B,E), Dif[(A,D),(B,E)].
C
C     ARGUMENTS
C
C     MODE PARAMETERS
C
C     REDUCE  CHARACTER*1
C             Indicates whether the matrix pairs (A,D) and/or (B,E) are
C             to be reduced to generalized Schur form as follows:
C             = 'R':  The matrix pairs (A,D) and (B,E) are to be reduced
C                     to generalized (real) Schur canonical form;
C             = 'A':  The matrix pair (A,D) only is to be reduced
C                     to generalized (real) Schur canonical form,
C                     and the matrix pair (B,E) already is in this form;
C             = 'B':  The matrix pair (B,E) only is to be reduced
C                     to generalized (real) Schur canonical form,
C                     and the matrix pair (A,D) already is in this form;
C             = 'N':  The matrix pairs (A,D) and (B,E) are already in
C                     generalized (real) Schur canonical form, as
C                     produced by LAPACK routine DGEES.
C
C     TRANS   CHARACTER*1
C             Indicates which of the equations, (1) or (2), is to be
C             solved as follows:
C             = 'N':  The generalized Sylvester equation (1) is to be
C                     solved;
C             = 'T':  The "transposed" generalized Sylvester equation
C                     (2) is to be solved.
C
C     JOBD    CHARACTER*1
C             Indicates whether the Dif estimator is to be computed as
C             follows:
C             = '1':  Only the one-norm-based Dif estimate is computed
C                     and stored in DIF;
C             = '2':  Only the Frobenius norm-based Dif estimate is
C                     computed and stored in DIF;
C             = 'D':  The equation (1) is solved and the one-norm-based
C                     Dif estimate is computed and stored in DIF;
C             = 'F':  The equation (1) is solved and the Frobenius norm-
C                     based Dif estimate is computed and stored in DIF;
C             = 'N':  The Dif estimator is not required and hence DIF is
C                     not referenced. (Solve either (1) or (2) only.)
C             JOBD is not referenced if TRANS = 'T'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of the matrices A and D and the number of rows
C             of the matrices C, F, R and L.  M >= 0.
C
C     N       (input) INTEGER
C             The order of the matrices B and E and the number of
C             columns of the matrices C, F, R and L.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
C             On entry, the leading M-by-M part of this array must
C             contain the coefficient matrix A of the equation; A must
C             be in upper quasi-triangular form if REDUCE = 'B' or 'N'.
C             On exit, the leading M-by-M part of this array contains
C             the upper quasi-triangular form of A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the coefficient matrix B of the equation; B must
C             be in upper quasi-triangular form if REDUCE = 'A' or 'N'.
C             On exit, the leading N-by-N part of this array contains
C             the upper quasi-triangular form of B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand side matrix C of the first equation
C             in (1) or (2).
C             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N
C             part of this array contains the solution matrix R of the
C             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading
C             M-by-N part of this array contains the solution matrix R
C             achieved during the computation of the Dif estimate.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading M-by-M part of this array must
C             contain the coefficient matrix D of the equation; D must
C             be in upper triangular form if REDUCE = 'B' or 'N'.
C             On exit, the leading M-by-M part of this array contains
C             the upper triangular form of D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the coefficient matrix E of the equation; E must
C             be in upper triangular form if REDUCE = 'A' or 'N'.
C             On exit, the leading N-by-N part of this array contains
C             the upper triangular form of E.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand side matrix F of the second
C             equation in (1) or (2).
C             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N
C             part of this array contains the solution matrix L of the
C             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading
C             M-by-N part of this array contains the solution matrix L
C             achieved during the computation of the Dif estimate.
C
C     LDF     INTEGER
C             The leading dimension of array F.  LDF >= MAX(1,M).
C
C     SCALE   (output) DOUBLE PRECISION
C             The scaling factor in (1) or (2). If 0 < SCALE < 1, C and
C             F hold the solutions R and L, respectively, to a slightly
C             perturbed system (but the input or computed generalized
C             (real) Schur canonical form matrices A, B, D, and E
C             have not been changed). If SCALE = 0, C and F hold the
C             solutions R and L, respectively, to the homogeneous system
C             with C = F = 0. Normally, SCALE = 1.
C
C     DIF     (output) DOUBLE PRECISION
C             If TRANS = 'N' and JOBD <> 'N', then DIF contains the
C             value of the Dif estimator, which is an upper bound of
C                                                    -1
C             Dif[(A,D),(B,E)] = sigma_min(Z) = 1/||Z  ||, in either the
C             one-norm, or Frobenius norm, respectively (see METHOD).
C             Otherwise, DIF is not referenced.
C
C     P       (output) DOUBLE PRECISION array, dimension (LDP,*)
C             If REDUCE = 'R' or 'A', then the leading M-by-M part of
C             this array contains the (left) transformation matrix used
C             to reduce (A,D) to generalized Schur form.
C             Otherwise, P is not referenced and can be supplied as a
C             dummy array (i.e. set parameter LDP = 1 and declare this
C             array to be P(1,1) in the calling program).
C
C     LDP     INTEGER
C             The leading dimension of array P.
C             LDP >= MAX(1,M) if REDUCE = 'R' or 'A',
C             LDP >= 1        if REDUCE = 'B' or 'N'.
C
C     Q       (output) DOUBLE PRECISION array, dimension (LDQ,*)
C             If REDUCE = 'R' or 'A', then the leading M-by-M part of
C             this array contains the (right) transformation matrix used
C             to reduce (A,D) to generalized Schur form.
C             Otherwise, Q is not referenced and can be supplied as a
C             dummy array (i.e. set parameter LDQ = 1 and declare this
C             array to be Q(1,1) in the calling program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= MAX(1,M) if REDUCE = 'R' or 'A',
C             LDQ >= 1        if REDUCE = 'B' or 'N'.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,*)
C             If REDUCE = 'R' or 'B', then the leading N-by-N part of
C             this array contains the (left) transformation matrix used
C             to reduce (B,E) to generalized Schur form.
C             Otherwise, U is not referenced and can be supplied as a
C             dummy array (i.e. set parameter LDU = 1 and declare this
C             array to be U(1,1) in the calling program).
C
C     LDU     INTEGER
C             The leading dimension of array U.
C             LDU >= MAX(1,N) if REDUCE = 'R' or 'B',
C             LDU >= 1        if REDUCE = 'A' or 'N'.
C
C     V       (output) DOUBLE PRECISION array, dimension (LDV,*)
C             If REDUCE = 'R' or 'B', then the leading N-by-N part of
C             this array contains the (right) transformation matrix used
C             to reduce (B,E) to generalized Schur form.
C             Otherwise, V is not referenced and can be supplied as a
C             dummy array (i.e. set parameter LDV = 1 and declare this
C             array to be V(1,1) in the calling program).
C
C     LDV     INTEGER
C             The leading dimension of array V.
C             LDV >= MAX(1,N) if REDUCE = 'R' or 'B',
C             LDV >= 1        if REDUCE = 'A' or 'N'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M+N+6)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If TRANS = 'N' and JOBD = 'D' or 'F', then
C                LDWORK = MAX(1,7*M,7*N,2*M*N) if REDUCE = 'R';
C                LDWORK = MAX(1,7*M,2*M*N)     if REDUCE = 'A';
C                LDWORK = MAX(1,7*N,2*M*N)     if REDUCE = 'B';
C                LDWORK = MAX(1,2*M*N)         if REDUCE = 'N'.
C             Otherwise, the term 2*M*N above should be omitted.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if REDUCE <> 'N' and either (A,D) and/or (B,E)
C                   cannot be reduced to generalized Schur form;
C             = 2:  if REDUCE = 'N' and either A or B is not in
C                   upper quasi-triangular form;
C             = 3:  if a singular matrix was encountered during the
C                   computation of the solution matrices R and L, that
C                   is (A,D) and (B,E) have common or close eigenvalues.
C
C     METHOD
C
C     For the case TRANS = 'N', and REDUCE = 'R' or 'N', the algorithm
C     used by the routine consists of four steps (see [1] and [2]) as
C     follows:
C
C        (a) if REDUCE = 'R', then the matrix pairs (A,D) and (B,E) are
C            transformed to generalized Schur form, i.e. orthogonal
C            matrices P, Q, U and V are computed such that P' * A * Q
C            and U' * B * V are in upper quasi-triangular form and
C            P' * D * Q and U' * E * V are in upper triangular form;
C        (b) if REDUCE = 'R', then the matrices C and F are transformed
C            to give P' * C * V and P' * F * V respectively;
C        (c) if REDUCE = 'R', then the transformed system
C
C            P' * A * Q * R1 - L1 * U' * B * V = scale * P' * C * V
C            P' * D * Q * R1 - L1 * U' * E * V = scale * P' * F * V
C
C            is solved to give R1 and L1; otherwise, equation (1) is
C            solved to give R and L directly. The Dif estimator
C            is also computed if JOBD <> 'N'.
C        (d) if REDUCE = 'R', then the solution is transformed back
C            to give R = Q * R1 * V' and L = P * L1 * U'.
C
C     By using Kronecker products, equation (1) can also be written as
C     the system of linear equations Z * x = scale*y (see [1]), where
C
C            | I*A    I*D  |
C        Z = |             |.
C            |-B'*I  -E'*I |
C
C                                              -1
C     If JOBD <> 'N', then a lower bound on ||Z  ||, in either the one-
C     norm or Frobenius norm, is computed, which in most cases is
C     a reliable estimate of the true value. Notice that since Z is a
C     matrix of order 2 * M * N, the exact value of Dif (i.e., in the
C     Frobenius norm case, the smallest singular value of Z) may be very
C     expensive to compute.
C
C     The case TRANS = 'N', and REDUCE = 'A' or 'B', is similar, but
C     only one of the matrix pairs should be reduced and the
C     calculations simplify.
C
C     For the case TRANS = 'T', and REDUCE = 'R' or 'N', the algorithm
C     is similar, but the steps (b), (c), and (d) are as follows:
C
C        (b) if REDUCE = 'R', then the matrices C and F are transformed
C            to give Q' * C * V and P' * F * U respectively;
C        (c) if REDUCE = 'R', then the transformed system
C
C            Q' * A' * P * R1 + Q' * D' * P * L1 =  scale * Q' * C * V
C            R1 * V' * B' * U + L1 * V' * E' * U = -scale * P' * F * U
C
C            is solved to give R1 and L1; otherwise, equation (2) is
C            solved to give R and L directly.
C        (d) if REDUCE = 'R', then the solution is transformed back
C            to give R = P * R1 * V' and L = P * L1 * V'.
C
C     REFERENCES
C
C     [1] Kagstrom, B. and Westin, L.
C         Generalized Schur Methods with Condition Estimators for
C         Solving the Generalized Sylvester Equation.
C         IEEE Trans. Auto. Contr., 34, pp. 745-751, 1989.
C     [2] Kagstrom, B. and Westin, L.
C         GSYLV - Fortran Routines for the Generalized Schur Method with
C         Dif Estimators for Solving the Generalized Sylvester
C         Equation.
C         Report UMINF-132.86, Institute of Information Processing,
C         Univ. of Umea, Sweden, July 1987.
C     [3] Golub, G.H., Nash, S. and Van Loan, C.F.
C         A Hessenberg-Schur Method for the Problem AX + XB = C.
C         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979.
C     [4] Kagstrom, B. and Van Dooren, P.
C         Additive Decomposition of a Transfer Function with respect to
C         a Specified Region.
C         In: "Signal Processing, Scattering and Operator Theory, and
C         Numerical Methods" (Eds. M.A. Kaashoek et al.).
C         Proceedings of MTNS-89, Vol. 3, pp. 469-477, Birkhauser Boston
C         Inc., 1990.
C     [5] Kagstrom, B. and Van Dooren, P.
C         A Generalized State-space Approach for the Additive
C         Decomposition of a Transfer Matrix.
C         Report UMINF-91.12, Institute of Information Processing, Univ.
C         of Umea, Sweden, April 1991.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable. A reliable estimate for the
C     condition number of Z in the Frobenius norm, is (see [1])
C
C        K(Z) = SQRT(  ||A||**2 + ||B||**2 + ||C||**2 + ||D||**2 )/DIF.
C
C     If mu is an upper bound on the relative error of the elements of
C     the matrices A, B, C, D, E and F, then the relative error in the
C     actual solution is approximately mu * K(Z).
C
C     The relative error in the computed solution (due to rounding
C     errors) is approximately EPS * K(Z), where EPS is the machine
C     precision (see LAPACK Library routine DLAMCH).
C
C     FURTHER COMMENTS
C
C     For applications of the generalized Sylvester equation in control
C     theory, see [4] and [5].
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB04CD by Bo Kagstrom and Lars
C     Westin.
C
C     REVISIONS
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, Dec. 1999,
C     May 2009.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, orthogonal transformation, real
C     Schur form, Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBD, REDUCE, TRANS
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDF, LDP, LDQ,
     $                  LDU, LDV, LDWORK, M, N
      DOUBLE PRECISION  DIF, SCALE
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), E(LDE,*), F(LDF,*), P(LDP,*),
     $                  Q(LDQ,*), U(LDU,*), V(LDV,*)
C     .. Local Scalars ..
      LOGICAL           ILASCL, ILBSCL, ILDSCL, ILESCL, LJOB1, LJOB2,
     $                  LJOBD, LJOBDF, LJOBF, LREDRA, LREDRB, LREDUA,
     $                  LREDUB, LREDUC, LREDUR, LTRANN, SUFWRK
      INTEGER           I, IERR, IJOB, MINWRK, MN, WRKOPT, SDIM
      DOUBLE PRECISION  ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, DNRM,
     $                  DNRMTO, ENRM, ENRMTO, SAFMAX, SAFMIN, SMLNUM
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGGES, DGEMM, DGEMV, DLABAD, DLACPY,
     $                  DLASCL, DTGSYL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      MN   = MAX( M, N )
      LREDUR = LSAME( REDUCE, 'R' )
      LREDUA = LSAME( REDUCE, 'A' )
      LREDUB = LSAME( REDUCE, 'B' )
      LREDRA = LREDUR.OR.LREDUA
      LREDRB = LREDUR.OR.LREDUB
      LREDUC = LREDRA.OR.LREDUB
      IF ( LREDUR ) THEN
         MINWRK = MAX( 1, 7*MN )
      ELSE IF ( LREDUA ) THEN
         MINWRK = MAX( 1, 7*M )
      ELSE IF ( LREDUB ) THEN
         MINWRK = MAX( 1, 7*N )
      ELSE
         MINWRK = 1
      END IF
      LTRANN = LSAME( TRANS,  'N' )
      IF ( LTRANN ) THEN
         LJOB1  = LSAME( JOBD, '1' )
         LJOB2  = LSAME( JOBD, '2' )
         LJOBD  = LSAME( JOBD, 'D' )
         LJOBF  = LSAME( JOBD, 'F' )
         LJOBDF = LJOB1.OR.LJOB2.OR.LJOBD.OR.LJOBF
         IF ( LJOBD.OR.LJOBF ) MINWRK = MAX( MINWRK, 2*M*N )
      END IF
C
C     Test the input scalar arguments.
C
      IF( .NOT.LREDUC .AND. .NOT.LSAME( REDUCE, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LTRANN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( LTRANN ) THEN
         IF( .NOT.LJOBDF .AND. .NOT.LSAME( JOBD, 'N' ) )
     $      INFO = -3
      END IF
      IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -17
      ELSE IF( ( .NOT.LREDRA .AND. LDP.LT.1 )             .OR.
     $         (      LREDRA .AND. LDP.LT.MAX( 1, M ) ) ) THEN
         INFO = -21
      ELSE IF( ( .NOT.LREDRA .AND. LDQ.LT.1 )             .OR.
     $         (      LREDRA .AND. LDQ.LT.MAX( 1, M ) ) ) THEN
         INFO = -23
      ELSE IF( ( .NOT.LREDRB .AND. LDU.LT.1 )             .OR.
     $         (      LREDRB .AND. LDU.LT.MAX( 1, N ) ) ) THEN
         INFO = -25
      ELSE IF( ( .NOT.LREDRB .AND. LDV.LT.1 )             .OR.
     $         (      LREDRB .AND. LDV.LT.MAX( 1, N ) ) ) THEN
         INFO = -27
      ELSE IF( LDWORK.LT.MINWRK ) THEN
         INFO = -30
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB04OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. M.EQ.0 ) THEN
         SCALE = ONE
         DWORK(1) = ONE
         IF ( LTRANN ) THEN
            IF ( LJOBDF ) DIF = ONE
         END IF
         RETURN
      END IF
      WRKOPT = 1
      SUFWRK = LDWORK.GE.M*N
C
C     STEP 1: Reduce (A,D) and/or (B,E) to generalized Schur form.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( LREDUC ) THEN
C
C        Get machine constants.
C
         SAFMIN = DLAMCH( 'Safe minimum' )
         SAFMAX = ONE / SAFMIN
         CALL DLABAD( SAFMIN, SAFMAX )
         SMLNUM = SQRT( SAFMIN ) / DLAMCH( 'Precision' )
         BIGNUM = ONE / SMLNUM
C
         IF ( .NOT.LREDUB ) THEN
C
C           Scale A if max element outside range [SMLNUM,BIGNUM].
C
            ANRM = DLANGE( 'M', M, M, A, LDA, DWORK )
            ILASCL = .FALSE.
            IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
               ANRMTO = SMLNUM
               ILASCL = .TRUE.
            ELSE IF( ANRM.GT.BIGNUM ) THEN
               ANRMTO = BIGNUM
               ILASCL = .TRUE.
            END IF
            IF( ILASCL )
     $         CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, M, M, A, LDA,
     $                      IERR )
C
C           Scale D if max element outside range [SMLNUM,BIGNUM]
C
            DNRM = DLANGE( 'M', M, M, D, LDD, DWORK )
            ILDSCL = .FALSE.
            IF( DNRM.GT.ZERO .AND. DNRM.LT.SMLNUM ) THEN
               DNRMTO = SMLNUM
               ILDSCL = .TRUE.
            ELSE IF( DNRM.GT.BIGNUM ) THEN
               DNRMTO = BIGNUM
               ILDSCL = .TRUE.
            END IF
            IF( ILDSCL )
     $         CALL DLASCL( 'G', 0, 0, DNRM, DNRMTO, M, M, D, LDD,
     $                      IERR )
C
C           Reduce (A,D) to generalized Schur form.
C           Workspace:  need   7*M;
C                       prefer 5*M + M*(NB+1).
C
C            CALL DGEGS( 'Vectors left', 'Vectors right', M, A, LDA, D,
C     $                  LDD, DWORK, DWORK(M+1), DWORK(2*M+1), P, LDP, Q,
C     $                  LDQ, DWORK(3*M+1), LDWORK-3*M, INFO )
         CALL DGGES( 'Vectors left', 'Vectors right', 'N', 0, N, A, LDA,
     $               D, LDD, SDIM, DWORK, DWORK(M+1), DWORK(2*M+1), P, LDP, Q,
     $               LDQ, DWORK(3*M+1), LDWORK-3*M, 0, INFO )

C
C           Undo scaling
C
            IF( ILASCL )
     $         CALL DLASCL( 'H', 0, 0, ANRMTO, ANRM, M, M, A, LDA,
     $                      IERR )
C
            IF( ILDSCL )
     $         CALL DLASCL( 'U', 0, 0, DNRMTO, DNRM, M, M, D, LDD,
     $                      IERR )
C
            IF ( INFO.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(3*M+1) ) + 3*M )
         END IF
         IF ( .NOT.LREDUA ) THEN
C
C           Scale B if max element outside range [SMLNUM,BIGNUM]
C
            BNRM = DLANGE( 'M', N, N, B, LDB, DWORK )
            ILBSCL = .FALSE.
            IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
               BNRMTO = SMLNUM
               ILBSCL = .TRUE.
            ELSE IF( BNRM.GT.BIGNUM ) THEN
               BNRMTO = BIGNUM
               ILBSCL = .TRUE.
            END IF
            IF( ILBSCL )
     $         CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB,
     $                      IERR )
C
C           Scale E if max element outside range [SMLNUM,BIGNUM]
C
            ENRM = DLANGE( 'M', N, N, E, LDE, DWORK )
            ILESCL = .FALSE.
            IF( ENRM.GT.ZERO .AND. ENRM.LT.SMLNUM ) THEN
               ENRMTO = SMLNUM
               ILESCL = .TRUE.
            ELSE IF( ENRM.GT.BIGNUM ) THEN
               ENRMTO = BIGNUM
               ILESCL = .TRUE.
            END IF
            IF( ILESCL )
     $         CALL DLASCL( 'G', 0, 0, ENRM, ENRMTO, N, N, E, LDE,
     $                      IERR )
C
C           Reduce (B,E) to generalized Schur form.
C           Workspace:  need   7*N;
C                       prefer 5*N + N*(NB+1).
C
C            CALL DGEGS( 'Vectors left', 'Vectors right', N, B, LDB, E,
C     $                  LDE, DWORK, DWORK(N+1), DWORK(2*N+1), U, LDU, V,
C     $                  LDV, DWORK(3*N+1), LDWORK-3*N, INFO )
            CALL DGGES( 'Vectors left', 'Vectors right', 'N', 0, N, B, LDB, E,
     $                  LDE, SDIM, DWORK, DWORK(N+1), DWORK(2*N+1), U, LDU, V,
     $                  LDV, DWORK(3*N+1), LDWORK-3*N, 0, INFO )
C
C           Undo scaling
C
            IF( ILBSCL )
     $         CALL DLASCL( 'H', 0, 0, BNRMTO, BNRM, N, N, B, LDB,
     $                      IERR )
C
            IF( ILESCL )
     $         CALL DLASCL( 'U', 0, 0, ENRMTO, ENRM, N, N, E, LDE,
     $                      IERR )
C
            IF ( INFO.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            WRKOPT = MAX( WRKOPT, INT( DWORK(3*N+1) ) + 3*N )
         END IF
      END IF
C
      IF (.NOT.LREDUR ) THEN
C
C        Set INFO = 2 if A and/or B are/is not in quasi-triangular form.
C
         IF (.NOT.LREDUA ) THEN
            I = 1
C
   20       CONTINUE
            IF ( I.LE.M-2 ) THEN
               IF ( A(I+1,I).NE.ZERO ) THEN
                  IF ( A(I+2,I+1).NE.ZERO ) THEN
                     INFO = 2
                     RETURN
                  ELSE
                     I = I + 1
                  END IF
               END IF
               I = I + 1
               GO TO 20
            END IF
         END IF
C
         IF (.NOT.LREDUB ) THEN
            I = 1
C
   40       CONTINUE
            IF ( I.LE.N-2 ) THEN
               IF ( B(I+1,I).NE.ZERO ) THEN
                  IF ( B(I+2,I+1).NE.ZERO ) THEN
                     INFO = 2
                     RETURN
                  ELSE
                     I = I + 1
                  END IF
               END IF
               I = I + 1
               GO TO 40
            END IF
         END IF
      END IF
C
C     STEP 2: Modify right hand sides (C,F).
C
      IF ( LREDUC ) THEN
         WRKOPT = MAX( WRKOPT, M*N )
         IF ( SUFWRK ) THEN
C
C           Enough workspace for a BLAS 3 calculation.
C
            IF ( LTRANN ) THEN
C
C              Equation (1).
C
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'Transpose', 'No transpose', M, N, M, ONE,
     $                        P, LDP, C, LDC, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, C, LDC, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, N,
     $                        ONE, DWORK, M, V, LDV, ZERO, C, LDC )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, C, LDC )
               END IF
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'Transpose', 'No transpose', M, N, M, ONE,
     $                        P, LDP, F, LDF, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, F, LDF, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, N,
     $                        ONE, DWORK, M, V, LDV, ZERO, F, LDF )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, F, LDF )
               END IF
            ELSE
C
C              Equation (2).
C
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'Transpose', 'No transpose', M, N, M, ONE,
     $                        Q, LDQ, C, LDC, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, C, LDC, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, N,
     $                        ONE, DWORK, M, V, LDV, ZERO, C, LDC )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, C, LDC )
               END IF
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'Transpose', 'No transpose', M, N, M, ONE,
     $                        P, LDP, F, LDF, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, F, LDF, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, N,
     $                        ONE, DWORK, M, U, LDU, ZERO, F, LDF )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, F, LDF )
               END IF
            END IF
         ELSE
C
C           Use a BLAS 2 calculation.
C
            IF ( LTRANN ) THEN
C
C              Equation (1).
C
               IF ( .NOT.LREDUB ) THEN
C
                  DO 60 I = 1, N
                     CALL DGEMV( 'Transpose', M, M, ONE, P, LDP, C(1,I),
     $                           1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, C(1,I), 1 )
   60             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 80 I = 1, M
                     CALL DGEMV( 'Transpose', N, N, ONE, V, LDV, C(I,1),
     $                           LDC, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, C(I,1), LDC )
   80             CONTINUE
C
               END IF
               IF ( .NOT.LREDUB ) THEN
C
                  DO 100 I = 1, N
                     CALL DGEMV( 'Transpose', M, M, ONE, P, LDP, F(1,I),
     $                           1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, F(1,I), 1 )
  100             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 120 I = 1, M
                     CALL DGEMV( 'Transpose', N, N, ONE, V, LDV, F(I,1),
     $                           LDF, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, F(I,1), LDF )
  120             CONTINUE
C
               END IF
            ELSE
C
C              Equation (2).
C
               IF ( .NOT.LREDUB ) THEN
C
                  DO 140 I = 1, N
                     CALL DGEMV( 'Transpose', M, M, ONE, Q, LDQ, C(1,I),
     $                           1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, C(1,I), 1 )
  140             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 160 I = 1, M
                     CALL DGEMV( 'Transpose', N, N, ONE, V, LDV, C(I,1),
     $                           LDC, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, C(I,1), LDC )
  160             CONTINUE
C
               END IF
               IF ( .NOT.LREDUB ) THEN
C
                  DO 180 I = 1, N
                     CALL DGEMV( 'Transpose', M, M, ONE, P, LDP, F(1,I),
     $                           1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, F(1,I), 1 )
  180             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 200 I = 1, M
                     CALL DGEMV( 'Transpose', N, N, ONE, U, LDU, F(I,1),
     $                           LDF, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, F(I,1), LDF )
  200             CONTINUE
C
               END IF
            END IF
         END IF
      END IF
C
C     STEP 3: Solve the transformed system and compute the Dif
C             estimator.
C
      IF ( LTRANN ) THEN
         IF ( LJOBD ) THEN
            IJOB = 1
         ELSE IF ( LJOBF ) THEN
            IJOB = 2
         ELSE IF ( LJOB1 ) THEN
            IJOB = 3
         ELSE IF ( LJOB2 ) THEN
            IJOB = 4
         ELSE
            IJOB = 0
         END IF
      ELSE
         IJOB = 0
      END IF
C
C     Workspace:  need 2*M*N if TRANS = 'N' and JOBD = 'D' or 'F';
C                      1, otherwise.
C
      CALL DTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
     $             E, LDE, F, LDF, SCALE, DIF, DWORK, LDWORK, IWORK,
     $             INFO )
      IF ( INFO.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      IF ( LTRANN ) THEN
         IF ( LJOBD.OR.LJOBF )
     $      WRKOPT = MAX( WRKOPT, 2*M*N )
      END IF
C
C     STEP 4: Back transformation of the solution.
C
      IF ( LREDUC ) THEN
         IF (SUFWRK ) THEN
C
C           Enough workspace for a BLAS 3 calculation.
C
            IF ( LTRANN ) THEN
C
C              Equation (1).
C
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, M,
     $                        ONE, Q, LDQ, C, LDC, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, C, LDC, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'Transpose', M, N, N, ONE,
     $                        DWORK, M, V, LDV, ZERO, C, LDC )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, C, LDC )
               END IF
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, M,
     $                        ONE, P, LDP, F, LDF, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, F, LDF, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'Transpose', M, N, N, ONE,
     $                        DWORK, M, U, LDU, ZERO, F, LDF )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, F, LDF )
               END IF
            ELSE
C
C              Equation (2).
C
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, M,
     $                        ONE, P, LDP, C, LDC, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, C, LDC, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'Transpose', M, N, N,
     $                        ONE, DWORK, M, V, LDV, ZERO, C, LDC )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, C, LDC )
               END IF
               IF ( .NOT.LREDUB ) THEN
                  CALL DGEMM( 'No transpose', 'No transpose', M, N, M,
     $                        ONE, P, LDP, F, LDF, ZERO, DWORK, M )
               ELSE
                  CALL DLACPY( 'Full', M, N, F, LDF, DWORK, M )
               END IF
               IF ( .NOT.LREDUA ) THEN
                  CALL DGEMM( 'No transpose', 'Transpose', M, N, N,
     $                        ONE, DWORK, M, V, LDV, ZERO, F, LDF )
               ELSE
                  CALL DLACPY( 'Full', M, N, DWORK, M, F, LDF )
               END IF
            END IF
         ELSE
C
C           Use a BLAS 2 calculation.
C
            IF ( LTRANN ) THEN
C
C              Equation (1).
C
               IF ( .NOT.LREDUB ) THEN
C
                  DO 220 I = 1, N
                     CALL DGEMV( 'No transpose', M, M, ONE, Q, LDQ,
     $                           C(1,I), 1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, C(1,I), 1 )
  220             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 240 I = 1, M
                     CALL DGEMV( 'No transpose', N, N, ONE, V, LDV,
     $                           C(I,1), LDC, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, C(I,1), LDC )
  240             CONTINUE
C
               END IF
               IF ( .NOT.LREDUB ) THEN
C
                  DO 260 I = 1, N
                     CALL DGEMV( 'No transpose', M, M, ONE, P, LDP,
     $                           F(1,I), 1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, F(1,I), 1 )
  260             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 280 I = 1, M
                     CALL DGEMV( 'No transpose', N, N, ONE, U, LDU,
     $                           F(I,1), LDF, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, F(I,1), LDF )
  280             CONTINUE
C
               END IF
            ELSE
C
C              Equation (2).
C
               IF ( .NOT.LREDUB ) THEN
C
                  DO 300 I = 1, N
                     CALL DGEMV( 'No transpose', M, M, ONE, P, LDP,
     $                           C(1,I), 1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, C(1,I), 1 )
  300             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 320 I = 1, M
                     CALL DGEMV( 'No transpose', N, N, ONE, V, LDV,
     $                           C(I,1), LDC, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, C(I,1), LDC )
  320             CONTINUE
C
               END IF
               IF ( .NOT.LREDUB ) THEN
C
                  DO 340 I = 1, N
                     CALL DGEMV( 'No transpose', M, M, ONE, P, LDP,
     $                           F(1,I), 1, ZERO, DWORK, 1 )
                     CALL DCOPY( M, DWORK, 1, F(1,I), 1 )
  340             CONTINUE
C
               END IF
               IF ( .NOT.LREDUA ) THEN
C
                  DO 360 I = 1, M
                     CALL DGEMV( 'No transpose', N, N, ONE, V, LDV,
     $                           F(I,1), LDF, ZERO, DWORK, 1 )
                     CALL DCOPY( N, DWORK, 1, F(I,1), LDF )
  360             CONTINUE
C
               END IF
            END IF
         END IF
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of SB04OD ***
      END
