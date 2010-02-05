      SUBROUTINE AB09IX( DICO, JOB, FACT, ORDSEL, N, M, P, NR,
     $                   SCALEC, SCALEO, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   TI, LDTI, T, LDT, NMINR, HSV, TOL1, TOL2,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     state-space representation (A,B,C,D) by using the square-root or
C     balancing-free square-root Balance & Truncate (B&T) or
C     Singular Perturbation Approximation (SPA) model reduction methods.
C     The computation of truncation matrices TI and T is based on
C     the Cholesky factor S of a controllability Grammian P = S*S'
C     and the Cholesky factor R of an observability Grammian Q = R'*R,
C     where S and R are given upper triangular matrices.
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
C             = 'B':  use the square-root B&T method;
C             = 'F':  use the balancing-free square-root B&T method;
C             = 'S':  use the square-root SPA method;
C             = 'P':  use the balancing-free square-root SPA method.
C
C     FACT    CHARACTER*1
C             Specifies whether or not, on entry, the matrix A is in a
C             real Schur form, as follows:
C             = 'S':  A is in a real Schur form;
C             = 'N':  A is a general dense square matrix.
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
C             The number of system outputs.  P >= 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of
C             the resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. NR is set as follows:
C             if ORDSEL = 'F', NR is equal to MIN(NR,NMINR), where NR
C             is the desired order on entry and NMINR is the number of
C             the Hankel singular values greater than N*EPS*S1, where
C             EPS is the machine precision (see LAPACK Library Routine
C             DLAMCH) and S1 is the largest Hankel singular value
C             (computed in HSV(1));
C             NR can be further reduced to ensure HSV(NR) > HSV(NR+1);
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS*S1).
C
C     SCALEC  (input) DOUBLE PRECISION
C             Scaling factor for the Cholesky factor S of the
C             controllability Grammian, i.e., S/SCALEC is used to
C             compute the Hankel singular values.  SCALEC > 0.
C
C     SCALEO  (input) DOUBLE PRECISION
C             Scaling factor for the Cholesky factor R of the
C             observability Grammian, i.e., R/SCALEO is used to
C             compute the Hankel singular values.  SCALEO > 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A. If FACT = 'S',
C             A is in a real Schur form.
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
C             On entry, if JOB = 'S' or JOB = 'P', the leading P-by-M
C             part of this array must contain the original input/output
C             matrix D.
C             On exit, if INFO = 0 and JOB = 'S' or JOB = 'P', the
C             leading P-by-M part of this array contains the
C             input/output matrix Dr of the reduced order system.
C             If JOB = 'B' or JOB = 'F', this array is not referenced.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= 1,        if JOB = 'B' or JOB = 'F';
C             LDD >= MAX(1,P), if JOB = 'S' or JOB = 'P'.
C
C     TI      (input/output) DOUBLE PRECISION array, dimension (LDTI,N)
C             On entry, the leading N-by-N upper triangular part of
C             this array must contain the Cholesky factor S of a
C             controllability Grammian P = S*S'.
C             On exit, if INFO = 0, and NR > 0, the leading NMINR-by-N
C             part of this array contains the left truncation matrix
C             TI in (1), for the B&T approach, or in (2), for the
C             SPA approach.
C
C     LDTI    INTEGER
C             The leading dimension of array TI.  LDTI >= MAX(1,N).
C
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C             On entry, the leading N-by-N upper triangular part of
C             this array must contain the Cholesky factor R of an
C             observability Grammian Q = R'*R.
C             On exit, if INFO = 0, and NR > 0, the leading N-by-NMINR
C             part of this array contains the right truncation matrix
C             T in (1), for the B&T approach, or in (2), for the
C             SPA approach.
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= MAX(1,N).
C
C     NMINR   (output) INTEGER
C             The number of Hankel singular values greater than
C             MAX(TOL2,N*EPS*S1).
C             Note: If S and R are the Cholesky factors of the
C             controllability and observability Grammians of the
C             original system (A,B,C,D), respectively, then NMINR is
C             the order of a minimal realization of the original system.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the Hankel singular values,
C             ordered decreasingly. The Hankel singular values are
C             singular values of the product R*S.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of the reduced system.
C             For model reduction, the recommended value lies in the
C             interval [0.00001,0.001].
C             If TOL1 <= 0 on entry, the used default value is
C             TOL1 = N*EPS*S1, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH) and S1 is the largest
C             Hankel singular value (computed in HSV(1)).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the system.
C             The recommended value is TOL2 = N*EPS*S1.
C             This value is used by default if TOL2 <= 0 on entry.
C             If TOL2 > 0, and ORDSEL = 'A', then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension LIWORK, where
C             LIWORK = 0,   if JOB = 'B';
C             LIWORK = N,   if JOB = 'F';
C             LIWORK = 2*N, if JOB = 'S' or 'P'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, 2*N*N + 5*N, N*MAX(M,P) ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than NMINR, the order of a minimal realization of
C                   the given system; in this case, the resulting NR is
C                   set automatically to NMINR;
C             = 2:  with ORDSEL = 'F', the selected order NR corresponds
C                   to repeated singular values, which are neither all
C                   included nor all excluded from the reduced model;
C                   in this case, the resulting NR is set automatically
C                   to the largest value such that HSV(NR) > HSV(NR+1).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t),                             (3)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09IX determines for
C     the given system (3), the matrices of a reduced NR order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t),                         (4)
C
C     by using the square-root or balancing-free square-root
C     Balance & Truncate (B&T) or Singular Perturbation Approximation
C     (SPA) model reduction methods.
C
C     The projection matrices TI and T are determined using the
C     Cholesky factors S and R of a controllability Grammian P and an
C     observability Grammian Q.
C     The Hankel singular values HSV(1), ...., HSV(N) are computed as
C     singular values of the product R*S.
C
C     If JOB = 'B', the square-root Balance & Truncate technique
C     of [1] is used.
C
C     If JOB = 'F', the balancing-free square-root version of the
C     Balance & Truncate technique [2] is used.
C
C     If JOB = 'S', the square-root version of the Singular Perturbation
C     Approximation method [3,4] is used.
C
C     If JOB = 'P', the balancing-free square-root version of the
C     Singular Perturbation Approximation method [3,4] is used.
C
C     REFERENCES
C
C     [1] Tombs M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [2] Varga A.
C         Efficient minimal realization procedure based on balancing.
C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
C         A. El Moudni, P. Borne, S. G. Tzafestas (Eds.),
C         Vol. 2, pp. 42-46.
C
C     [3] Liu Y. and Anderson B.D.O.
C         Singular Perturbation Approximation of balanced systems.
C         Int. J. Control, Vol. 50, pp. 1379-1405, 1989.
C
C     [4] Varga A.
C         Balancing-free square-root algorithm for computing singular
C         perturbation approximations.
C         Proc. 30-th CDC, Brighton, Dec. 11-13, 1991,
C         Vol. 2, pp. 1062-1065.
C
C     NUMERICAL ASPECTS
C
C     The implemented method relies on accuracy enhancing square-root
C     or balancing-free square-root methods.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000.
C     D. Sima, University of Bucharest, August 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000,
C              Sep. 2001.
C
C     KEYWORDS
C
C     Balance and truncate, minimal state-space representation,
C     model reduction, multivariable system,
C     singular perturbation approximation, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, FACT, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDT, LDTI,
     $                  LDWORK, M, N, NMINR, NR, P
      DOUBLE PRECISION  SCALEC, SCALEO, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*), T(LDT,*), TI(LDTI,*)
C     .. Local Scalars ..
      LOGICAL           BAL, BTA, DISCR, FIXORD, RSF, SPA
      INTEGER           IERR, IJ, J, K, KTAU, KU, KV, KW, LDW, LW,
     $                  NRED, NR1, NS, WRKOPT
      DOUBLE PRECISION  ATOL, RCOND, SKP, TEMP, TOLDEF
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09DD, DGEMM,  DGEMV, DGEQRF, DGETRF, DGETRS,
     $                  DLACPY, DORGQR, DSCAL, DTRMM,  DTRMV,  MA02AD,
     $                  MB03UD, XERBLA
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
      RSF    = LSAME( FACT,   'S' )
      FIXORD = LSAME( ORDSEL, 'F' )
C
      LW = MAX( 1, 2*N*N + 5*N, N*MAX( M, P ) )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( BTA .OR. SPA ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( RSF .OR. LSAME( FACT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -8
      ELSE IF( SCALEC.LE.ZERO ) THEN
         INFO = -9
      ELSE IF( SCALEO.LE.ZERO ) THEN
         INFO = -10
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDD.LT.1 .OR. ( SPA .AND. LDD.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDTI.LT.MAX( 1, N ) ) THEN
         INFO = -20
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -22
      ELSE IF( TOL2.GT.ZERO .AND. .NOT.FIXORD .AND. TOL2.GT.TOL1 ) THEN
         INFO = -26
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -29
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09IX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         NMINR = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Save S in DWORK(KV).
C
      KV = 1
      KU = KV + N*N
      KW = KU + N*N
      CALL DLACPY( 'Upper', N, N, TI, LDTI, DWORK(KV), N )
C                             | x x |
C     Compute R*S in the form | 0 x | in TI.
C
      DO 10 J = 1, N
         CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', J, T, LDT,
     $               TI(1,J), 1 )
   10 CONTINUE
C
C     Compute the singular value decomposition R*S = V*Sigma*UT of the
C     upper triangular matrix R*S, with UT in TI and V in DWORK(KU).
C
C     Workspace:  need   2*N*N + 5*N;
C                 prefer larger.
C
      CALL MB03UD( 'Vectors', 'Vectors', N, TI, LDTI, DWORK(KU), N, HSV,
     $             DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 1
         RETURN
      ENDIF
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Scale the singular values.
C
      CALL DSCAL( N, ONE / SCALEC / SCALEO, HSV, 1 )
C
C     Partition Sigma, U and V conformally as:
C
C     Sigma = diag(Sigma1,Sigma2,Sigma3),  U = [U1,U2,U3] (U' in TI) and
C     V = [V1,V2,V3] (in DWORK(KU)).
C
C     Compute NMINR, the order of a minimal realization, as the order
C     of [Sigma1 Sigma2].
C
      TOLDEF = DBLE( N )*DLAMCH( 'Epsilon' )
      ATOL   = MAX( TOL2, TOLDEF*HSV(1) )
      NMINR  = N
   20 IF( NMINR.GT.0 ) THEN
         IF( HSV(NMINR).LE.ATOL ) THEN
            NMINR = NMINR - 1
            GO TO 20
         END IF
      END IF
C
C     Compute the order NR of reduced system, as the order of Sigma1.
C
      IF( FIXORD ) THEN
C
C        Check if the desired order is less than the order of a minimal
C        realization.
C
         IF( NR.GT.NMINR ) THEN
C
C           Reduce the order to NMINR.
C
            NR = NMINR
            IWARN = 1
         END IF
C
C        Check for singular value multiplicity at cut-off point.
C
         IF( NR.GT.0 .AND. NR.LT.NMINR ) THEN
            SKP = HSV(NR)
            IF( SKP-HSV(NR+1).LE.TOLDEF*SKP ) THEN
               IWARN = 2
C
C              Reduce the order such that HSV(NR) > HSV(NR+1).
C
   30          NR = NR - 1
               IF( NR.GT.0 ) THEN
                  IF( HSV(NR)-SKP.LE.TOLDEF*SKP ) GO TO 30
               END IF
            END IF
         END IF
      ELSE
C
C        The order is given as the number of singular values
C        exceeding MAX( TOL1, N*EPS*HSV(1) ).
C
         ATOL = MAX( TOL1, ATOL )
         NR   = 0
         DO 40 J = 1, NMINR
            IF( HSV(J).LE.ATOL ) GO TO 50
            NR = NR + 1
   40    CONTINUE
   50    CONTINUE
      ENDIF
C
C     Finish if the order is zero.
C
      IF( NR.EQ.0 ) THEN
         IF( SPA )
     $      CALL AB09DD( DICO, N, M, P, NR, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, RCOND, IWORK, DWORK, IERR )
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Compute NS, the order of Sigma2. For BTA, NS = 0.
C
      IF( SPA ) THEN
         NRED = NMINR
      ELSE
         NRED = NR
      END IF
      NS = NRED - NR
C
C     Compute the truncation matrices.
C
C     Compute TI' = | TI1' TI2' | = R'*| V1 V2 | in DWORK(KU).
C
      CALL DTRMM( 'Left', 'Upper', 'Transpose', 'NonUnit', N, NRED,
     $            ONE, T, LDT, DWORK(KU), N )
C
C     Compute  T = | T1 T2 | = S*| U1 U2 | .
C
      CALL MA02AD( 'Full', NRED, N, TI, LDTI, T, LDT )
      CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N,
     $            NRED, ONE, DWORK(KV), N, T, LDT )
C
      KTAU = KW
      IF( BAL ) THEN
         IJ = KU
C
C        Square-Root B&T/SPA method.
C
C        Compute the truncation matrices for balancing
C                        -1/2                -1/2
C               T1*Sigma1     and TI1'*Sigma1    .
C
         DO 60 J = 1, NR
            TEMP = ONE/SQRT( HSV(J) )
            CALL DSCAL( N, TEMP, T(1,J), 1 )
            CALL DSCAL( N, TEMP, DWORK(IJ), 1 )
            IJ = IJ + N
   60    CONTINUE
C
      ELSE
C
C        Balancing-Free B&T/SPA method.
C
C        Compute orthogonal bases for the images of matrices T1 and
C        TI1'.
C
C        Workspace:  need   2*N*N + 2*N;
C                    prefer larger.
C
         KW   = KTAU + NR
         LDW  = LDWORK - KW + 1
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
C
      IF( NS.GT.0 ) THEN
C
C        Compute orthogonal bases for the images of matrices T2 and
C        TI2'.
C
C        Workspace:  need   2*N*N + 2*N;
C                    prefer larger.
C
         NR1 = NR + 1
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
      CALL MA02AD( 'Full', N, NRED, DWORK(KU), N, TI, LDTI )
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
C     Compute TI*A*T. Exploit RSF of A if possible.
C     Workspace:  need   N*N.
C
      IF( RSF ) THEN
         IJ = 1
         DO 80 J = 1, N
            K = MIN( J+1, N )
            CALL DGEMV( 'NoTranspose', NRED, K, ONE, TI, LDTI,
     $                  A(1,J), 1, ZERO, DWORK(IJ), 1 )
            IJ = IJ + N
   80    CONTINUE
      ELSE
         CALL DGEMM( 'NoTranspose', 'NoTranspose', NRED, N, N, ONE,
     $               TI, LDTI, A, LDA, ZERO, DWORK, N )
      END IF
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NRED, NRED, N, ONE,
     $            DWORK, N, T, LDT, ZERO, A, LDA )
C
C     Compute TI*B and C*T.
C     Workspace:  need   N*MAX(M,P).
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NRED, M, N, ONE, TI,
     $            LDTI, DWORK, N, ZERO, B, LDB )
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK, P )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NRED, N, ONE,
     $            DWORK, P, T, LDT, ZERO, C, LDC )
C
C     Compute the singular perturbation approximation if possible.
C     Note that IERR = 1 on exit from AB09DD cannot appear here.
C
C     Workspace:  need real    4*(NMINR-NR);
C                 need integer 2*(NMINR-NR).
C
      IF( SPA) THEN
         CALL AB09DD( DICO, NRED, M, P, NR, A, LDA, B, LDB,
     $                C, LDC, D, LDD, RCOND, IWORK, DWORK, IERR )
      ELSE
         NMINR = NR
      END IF
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09IX ***
      END
