      SUBROUTINE AB09HX( DICO, JOB, ORDSEL, N, M, P, NR, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, HSV, T, LDT, TI, LDTI, TOL1,
     $                   TOL2, IWORK, DWORK, LDWORK, BWORK, IWARN,
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for an original
C     stable state-space representation (A,B,C,D) by using the
C     stochastic balancing approach in conjunction with the square-root
C     or the balancing-free square-root Balance & Truncate (B&T) or
C     Singular Perturbation Approximation (SPA) model reduction methods.
C     The state dynamics matrix A of the original system is an upper
C     quasi-triangular matrix in real Schur canonical form and D must be
C     full row rank.
C
C     For the B&T approach, the matrices of the reduced order system
C     are computed using the truncation formulas:
C
C          Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T .     (1)
C
C     For the SPA approach, the matrices of a minimal realization
C     (Am,Bm,Cm) are computed using the truncation formulas:
C
C          Am = TI * A * T ,  Bm = TI * B ,  Cm = C * T .     (2)
C
C     Am, Bm, Cm and D serve further for computing the SPA of the given
C     system.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOB     CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root Balance & Truncate method;
C             = 'F':  use the balancing-free square-root
C                     Balance & Truncate method;
C             = 'S':  use the square-root Singular Perturbation
C                     Approximation method;
C             = 'P':  use the balancing-free square-root
C                     Singular Perturbation Approximation method.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation,
C             i.e., the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  M >= P >= 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of
C             the resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. NR is set as follows:
C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
C             is the desired order on entry and NMIN is the order of a
C             minimal realization of the given system; NMIN is
C             determined as the number of Hankel singular values greater
C             than N*EPS, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH);
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A in a real Schur
C             canonical form.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the
C             reduced order system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, if INFO = 0, the leading NR-by-M part of this
C             array contains the input/state matrix Br of the reduced
C             order system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, if INFO = 0, the leading P-by-NR part of this
C             array contains the state/output matrix Cr of the reduced
C             order system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the original input/output matrix D.
C             On exit, if INFO = 0, the leading P-by-M part of this
C             array contains the input/output matrix Dr of the reduced
C             order system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the Hankel singular values,
C             ordered decreasingly, of the phase system. All singular
C             values are less than or equal to 1.
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,N)
C             If INFO = 0 and NR > 0, the leading N-by-NR part of this
C             array contains the right truncation matrix T in (1), for
C             the B&T approach, or in (2), for the SPA approach.
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= MAX(1,N).
C
C     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N)
C             If INFO = 0 and NR > 0, the leading NR-by-N part of this
C             array contains the left truncation matrix TI in (1), for
C             the B&T approach, or in (2), for the SPA approach.
C
C     LDTI    INTEGER
C             The leading dimension of array TI.  LDTI >= MAX(1,N).
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value lies in the
C             interval [0.00001,0.001].
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = N*EPS, where EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the phase system (see METHOD) corresponding
C             to the given system.
C             The recommended value is TOL2 = N*EPS.
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(1,2*N)
C             On exit with INFO = 0, IWORK(1) contains the order of the
C             minimal realization of the system.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) contains RCOND, the reciprocal
C             condition number of the U11 matrix from the expression
C             used to compute the solution X = U21*inv(U11) of the
C             Riccati equation for spectral factorization.
C             A small value RCOND indicates possible ill-conditioning
C             of the respective Riccati equation.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 2, N*(MAX(N,M,P)+5),
C                            2*N*P+MAX(P*(M+2),10*N*(N+1) ) ).
C             For optimum performance LDWORK should be larger.
C
C     BWORK   LOGICAL array, dimension 2*N
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than the order of a minimal realization of the
C                   given system. In this case, the resulting NR is
C                   set automatically to a value corresponding to the
C                   order of a minimal realization of the system.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the state matrix A is not stable (if DICO = 'C')
C                   or not convergent (if DICO = 'D'), or it is not in
C                   a real Schur form;
C             = 2:  the reduction of Hamiltonian matrix to real
C                   Schur form failed;
C             = 3:  the reordering of the real Schur form of the
C                   Hamiltonian matrix failed;
C             = 4:  the Hamiltonian matrix has less than N stable
C                   eigenvalues;
C             = 5:  the coefficient matrix U11 in the linear system
C                   X*U11 = U21, used to determine X, is singular to
C                   working precision;
C             = 6:  the feedthrough matrix D has not a full row rank P;
C             = 7:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t),                             (3)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09HX determines for
C     the given system (3), the matrices of a reduced NR-rder system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                         (4)
C
C     such that
C
C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     If JOB = 'B', the square-root stochastic Balance & Truncate
C     method of [1] is used and the resulting model is balanced.
C
C     If JOB = 'F', the balancing-free square-root version of the
C     stochastic Balance & Truncate method [1] is used.
C
C     If JOB = 'S', the stochastic balancing method, in conjunction
C     with the square-root version of the Singular Perturbation
C     Approximation method [2,3] is used.
C
C     If JOB = 'P', the stochastic balancing method, in conjunction
C     with the balancing-free square-root version of the Singular
C     Perturbation Approximation method [2,3] is used.
C
C     By setting TOL1 = TOL2, the routine can be also used to compute
C     Balance & Truncate approximations.
C
C     REFERENCES
C
C     [1] Varga A. and Fasol K.H.
C         A new square-root balancing-free stochastic truncation
C         model reduction algorithm.
C         Proc. of 12th IFAC World Congress, Sydney, 1993.
C
C     [2] Liu Y. and Anderson B.D.O.
C         Singular Perturbation Approximation of balanced systems.
C         Int. J. Control, Vol. 50, pp. 1379-1405, 1989.
C
C     [3] Varga A.
C         Balancing-free square-root algorithm for computing singular
C         perturbation approximations.
C         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
C         Vol. 2, pp. 1062-1065.
C
C     NUMERICAL ASPECTS
C
C     The implemented method relies on accuracy enhancing square-root
C     or balancing-free square-root methods. The effectiveness of the
C     accuracy enhancing technique depends on the accuracy of the
C     solution of a Riccati equation. Ill-conditioned Riccati solution
C     typically results when D is nearly rank deficient.
C                                      3
C     The algorithm requires about 100N  floating point operations.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000.
C     D. Sima, University of Bucharest, May 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2000.
C     Partly based on the RASP routine SRBFS1, by A. Varga, 1992.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001.
C
C     KEYWORDS
C
C     Balance and truncate, minimal state-space representation,
C     model reduction, multivariable system,
C     singular perturbation approximation, state-space model,
C     stochastic balancing.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDT, LDTI,
     $                  LDWORK, M, N, NR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*), T(LDT,*), TI(LDTI,*)
      LOGICAL           BWORK(*)
C     .. Local Scalars ..
      LOGICAL           BAL, BTA, DISCR, FIXORD, SPA
      INTEGER           IERR, IJ, J, K, KTAU, KU, KV, KW, LDW, LW,
     $                  NMINR, NR1, NS, WRKOPT
      DOUBLE PRECISION  ATOL, RCOND, RICOND, SCALEC, SCALEO, TEMP,
     $                  TOLDEF
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB04MD, AB09DD, AB09HY, DGEMM, DGEMV, DGEQRF,
     $                  DGETRF, DGETRS, DLACPY, DORGQR, DSCAL, DTRMM,
     $                  DTRMV, MA02AD, MB03UD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      BTA    = LSAME( JOB,    'B' ) .OR. LSAME( JOB, 'F' )
      SPA    = LSAME( JOB,    'S' ) .OR. LSAME( JOB, 'P' )
      BAL    = LSAME( JOB,    'B' ) .OR. LSAME( JOB, 'S' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LW = MAX( 2, N*(MAX( N, M, P )+5),
     $          2*N*P+MAX( P*(M+2), 10*N*(N+1) ) )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 .OR. P.GT.M ) THEN
         INFO = -6
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -18
      ELSE IF( LDTI.LT.MAX( 1, N ) ) THEN
         INFO = -20
      ELSE IF( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 ) THEN
         INFO = -22
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -25
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09HX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         IWORK(1) = 0
         DWORK(1) = TWO
         DWORK(2) = ONE
         RETURN
      END IF
C
C     For discrete-time case, apply the discrete-to-continuous bilinear
C     transformation.
C
      IF( DISCR ) THEN
C
C        Real workspace:    need  N, prefer larger;
C        Integer workspace: need  N.
C
         CALL AB04MD( 'Discrete', N, M, P, ONE, ONE, A, LDA, B, LDB,
     $                C, LDC, D, LDD, IWORK, DWORK, LDWORK, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         WRKOPT = MAX( N, INT( DWORK(1) ) )
      ELSE
         WRKOPT = 0
      END IF
C
C     Compute in TI and T the Cholesky factors Su and Ru of the
C     controllability and observability Grammians, respectively.
C     Real workspace:    need  MAX( 2, N*(MAX(N,M,P)+5),
C                                   2*N*P+MAX(P*(M+2),10*N*(N+1) ) );
C                        prefer larger.
C     Integer workspace: need  2*N.
C
      CALL AB09HY( N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $             SCALEC, SCALEO, TI, LDTI, T, LDT, IWORK,
     $             DWORK, LDWORK, BWORK, INFO )
      IF( INFO.NE.0)
     $   RETURN
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      RICOND = DWORK(2)
C
C     Save Su in V.
C
      KU = 1
      KV = KU + N*N
      KW = KV + N*N
      CALL DLACPY( 'Upper', N, N, TI, LDTI, DWORK(KV), N )
C                               | x x |
C     Compute Ru*Su in the form | 0 x | in TI.
C
      DO 10 J = 1, N
         CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', J, T, LDT,
     $               TI(1,J), 1 )
   10 CONTINUE
C
C     Compute the singular value decomposition Ru*Su = V*S*UT
C     of the upper triangular matrix Ru*Su, with UT in TI and V in U.
C
C     Workspace:  need   2*N*N + 5*N;
C                 prefer larger.
C
      CALL MB03UD( 'Vectors', 'Vectors', N, TI, LDTI, DWORK(KU), N, HSV,
     $             DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 7
         RETURN
      ENDIF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Scale the singular values.
C
      CALL DSCAL( N, ONE / SCALEC / SCALEO, HSV, 1 )
C
C     Partition S, U and V conformally as:
C
C     S = diag(S1,S2,S3),  U = [U1,U2,U3] (U' in TI) and V = [V1,V2,V3]
C     (in U).
C
C     Compute the order NR of reduced system, as the order of S1.
C
      TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
      ATOL   = TOLDEF
      IF( FIXORD ) THEN
         IF( NR.GT.0 ) THEN
            IF( HSV(NR).LE.ATOL ) THEN
               NR = 0
               IWARN  = 1
               FIXORD = .FALSE.
            ENDIF
         ENDIF
      ELSE
         ATOL = MAX( TOL1, ATOL )
         NR = 0
      ENDIF
      IF( .NOT.FIXORD ) THEN
         DO 20 J = 1, N
            IF( HSV(J).LE.ATOL ) GO TO 30
            NR = NR + 1
   20    CONTINUE
   30    CONTINUE
      ENDIF
C
C     Compute the order of minimal realization as the order of [S1 S2].
C
      NR1   = NR + 1
      NMINR = NR
      IF( NR.LT.N ) THEN
         IF( SPA ) ATOL = MAX( TOL2, TOLDEF )
         DO 40 J = NR1, N
            IF( HSV(J).LE.ATOL ) GO TO 50
            NMINR = NMINR + 1
   40    CONTINUE
   50    CONTINUE
      END IF
C
C     Finish if the order is zero.
C
      IF( NR.EQ.0 ) THEN
         IF( SPA ) THEN
             CALL AB09DD( 'Continuous', N, M, P, NR, A, LDA, B, LDB,
     $                    C, LDC, D, LDD, RCOND, IWORK, DWORK, IERR )
            IWORK(1) = NMINR
         ELSE
            IWORK(1) = 0
         END IF
         DWORK(1) = WRKOPT
         DWORK(2) = RICOND
         RETURN
      END IF
C
C     Compute NS, the order of S2.
C     Note: For BTA, NS is always zero, because NMINR = NR.
C
      NS = NMINR - NR
C
C     Compute the truncation matrices.
C
C     Compute TI' = | TI1' TI2' | = Ru'*| V1 V2 | in U.
C
      CALL DTRMM( 'Left', 'Upper', 'Transpose', 'NonUnit', N, NMINR,
     $            ONE, T, LDT, DWORK(KU), N )
C
C     Compute  T = | T1 T2 | = Su*| U1 U2 | .
C
      CALL MA02AD( 'Full', NMINR, N, TI, LDTI, T, LDT )
      CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N,
     $            NMINR, ONE, DWORK(KV), N, T, LDT )
      KTAU = KV
C
      IF( BAL ) THEN
         IJ = KU
C
C        Square-Root B&T/SPA method.
C
C        Compute the truncation matrices for balancing
C                    -1/2            -1/2
C               T1*S1     and TI1'*S1    .
C
         DO 70 J = 1, NR
            TEMP = ONE/SQRT( HSV(J) )
            CALL DSCAL( N, TEMP, T(1,J), 1 )
            CALL DSCAL( N, TEMP, DWORK(IJ), 1 )
            IJ = IJ + N
   70    CONTINUE
      ELSE
C
C        Balancing-Free B&T/SPA method.
C
C        Compute orthogonal bases for the images of matrices T1 and
C        TI1'.
C
C        Workspace:  need   N*MAX(N,M,P) + 2*NR;
C                    prefer N*MAX(N,M,P) + NR*(NB+1)
C                           (NB determined by ILAENV for DGEQRF).
C
         KW  = KTAU + NR
         LDW = LDWORK - KW + 1
         CALL DGEQRF( N, NR, T, LDT, DWORK(KTAU), DWORK(KW), LDW, IERR )
         CALL DORGQR( N, NR, NR, T, LDT, DWORK(KTAU), DWORK(KW), LDW,
     $                IERR )
         CALL DGEQRF( N, NR, DWORK(KU), N, DWORK(KTAU), DWORK(KW), LDW,
     $                IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         CALL DORGQR( N, NR, NR, DWORK(KU), N, DWORK(KTAU), DWORK(KW),
     $                LDW, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      ENDIF
      IF( NS.GT.0 ) THEN
C
C        Compute orthogonal bases for the images of matrices T2 and
C        TI2'.
C
C        Workspace:  need   N*MAX(N,M,P) + 2*NS;
C                    prefer N*MAX(N,M,P) + NS*(NB+1)
C                           (NB determined by ILAENV for DGEQRF).
         KW  = KTAU + NS
         LDW = LDWORK - KW + 1
         CALL DGEQRF( N, NS, T(1,NR1), LDT, DWORK(KTAU), DWORK(KW), LDW,
     $                IERR )
         CALL DORGQR( N, NS, NS, T(1,NR1), LDT, DWORK(KTAU), DWORK(KW),
     $                LDW, IERR )
         CALL DGEQRF( N, NS, DWORK(KU+N*NR), N, DWORK(KTAU), DWORK(KW),
     $                LDW, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
         CALL DORGQR( N, NS, NS, DWORK(KU+N*NR), N, DWORK(KTAU),
     $                DWORK(KW), LDW, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      ENDIF
C
C     Transpose TI' in TI.
C
      CALL MA02AD( 'Full', N, NMINR, DWORK(KU), N, TI, LDTI )
C
      IF( .NOT.BAL ) THEN
C                        -1
C        Compute (TI1*T1)  *TI1 in TI.
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, NR, N, ONE, TI,
     $               LDTI, T, LDT, ZERO, DWORK(KU), N )
         CALL DGETRF( NR, NR, DWORK(KU), N, IWORK, IERR )
         CALL DGETRS( 'NoTranspose', NR, N, DWORK(KU), N, IWORK, TI,
     $                LDTI, IERR )
C
         IF( NS.GT.0 ) THEN
C                           -1
C           Compute (TI2*T2)  *TI2 in TI2.
C
            CALL DGEMM( 'NoTranspose', 'NoTranspose', NS, NS, N, ONE,
     $                  TI(NR1,1), LDTI, T(1,NR1), LDT, ZERO, DWORK(KU),
     $                  N )
            CALL DGETRF( NS, NS, DWORK(KU), N, IWORK, IERR )
            CALL DGETRS( 'NoTranspose', NS, N, DWORK(KU), N, IWORK,
     $                   TI(NR1,1), LDTI, IERR )
         END IF
      END IF
C
C     Compute TI*A*T (A is in RSF).
C
      IJ = KU
      DO 80 J = 1, N
         K = MIN( J+1, N )
         CALL DGEMV( 'NoTranspose', NMINR, K, ONE, TI, LDTI, A(1,J), 1,
     $               ZERO, DWORK(IJ), 1 )
         IJ = IJ + N
   80 CONTINUE
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NMINR, NMINR, N, ONE,
     $            DWORK(KU), N, T, LDT, ZERO, A, LDA )
C
C     Compute TI*B and C*T.
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), N )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NMINR, M, N, ONE, TI,
     $            LDTI, DWORK(KU), N, ZERO, B, LDB )
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), P )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NMINR, N, ONE,
     $            DWORK(KU), P, T, LDT, ZERO, C, LDC )
C
C     Compute the singular perturbation approximation if possible.
C     Note that IERR = 1 on exit from AB09DD cannot appear here.
C
C     Workspace:  need real    4*(NMINR-NR);
C                 need integer 2*(NMINR-NR).
C
      CALL AB09DD( 'Continuous', NMINR, M, P, NR, A, LDA, B, LDB,
     $             C, LDC, D, LDD, RCOND, IWORK, DWORK, IERR )
C
C     For discrete-time case, apply the continuous-to-discrete
C     bilinear transformation.
C
      IF( DISCR ) THEN
         CALL AB04MD( 'Continuous', NR, M, P, ONE, ONE, A, LDA, B, LDB,
     $                C, LDC, D, LDD, IWORK, DWORK, LDWORK, IERR )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      END IF
      IWORK(1) = NMINR
      DWORK(1) = WRKOPT
      DWORK(2) = RICOND
C
      RETURN
C *** Last line of AB09HX ***
      END
