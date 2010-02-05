      SUBROUTINE AB09BX( DICO, JOB, ORDSEL, N, M, P, NR, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, HSV, T, LDT, TI, LDTI, TOL1,
     $                   TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
C     original state-space representation (A,B,C,D) by using either the
C     square-root or the balancing-free square-root
C     Singular Perturbation Approximation (SPA) model reduction method.
C     The state dynamics matrix A of the original system is an upper
C     quasi-triangular matrix in real Schur canonical form. The matrices
C     of a minimal realization are computed using the truncation
C     formulas:
C
C          Am = TI * A * T ,  Bm = TI * B ,  Cm = C * T .      (1)
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
C             = 'B':  use the square-root SPA method;
C             = 'N':  use the balancing-free square-root SPA method.
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
C             The order of the original state-space representation, i.e.
C             the order of the matrix A.  N >= 0.
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
C             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR
C             is the desired order on entry and NMIN is the order of a
C             minimal realization of the given system; NMIN is
C             determined as the number of Hankel singular values greater
C             than N*EPS*HNORM(A,B,C), where EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH) and
C             HNORM(A,B,C) is the Hankel norm of the system (computed
C             in HSV(1));
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)).
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
C             If INFO = 0, it contains the Hankel singular values of
C             the original system ordered decreasingly. HSV(1) is the
C             Hankel norm of the system.
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,N)
C             If INFO = 0 and NR > 0, the leading N-by-NR part of this
C             array contains the right truncation matrix T in (1).
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= MAX(1,N).
C
C     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N)
C             If INFO = 0 and NR > 0, the leading NR-by-N part of this
C             array contains the left truncation matrix TI in (1).
C
C     LDTI    INTEGER
C             The leading dimension of array TI.  LDTI >= MAX(1,N).
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(A,B,C), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(A,B,C) is the
C             Hankel-norm of the given system (computed in HSV(1)).
C             For computing a minimal realization, the recommended
C             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH).
C             This value is used by default if TOL1 <= 0 on entry.
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the given system. The recommended value is
C             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default
C             if TOL2 <= 0 on entry.
C             If TOL2 > 0, then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension MAX(1,2*N)
C             On exit with INFO = 0, IWORK(1) contains the order of the
C             minimal realization of the system.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*(MAX(N,M,P)+5) + N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
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
C                   or not convergent (if DICO = 'D');
C             = 2:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t)                              (2)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09BX determines for
C     the given system (1), the matrices of a reduced NR order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t)                          (3)
C
C     such that
C
C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     If JOB = 'B', the balancing-based square-root SPA method of [1]
C     is used and the resulting model is balanced.
C
C     If JOB = 'N', the balancing-free square-root SPA method of [2]
C     is used.
C     By setting TOL1 = TOL2, the routine can be also used to compute
C     Balance & Truncate approximations.
C
C     REFERENCES
C
C     [1] Liu Y. and Anderson B.D.O.
C         Singular Perturbation Approximation of Balanced Systems,
C         Int. J. Control, Vol. 50, pp. 1379-1405, 1989.
C
C     [2] Varga A.
C         Balancing-free square-root algorithm for computing singular
C         perturbation approximations.
C         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991,
C         Vol. 2, pp. 1062-1065.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root or
C     balancing-free square-root techniques.
C                                         3
C     The algorithms require less than 30N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine SRBFP1.
C
C     REVISIONS
C
C     May 2, 1998.
C     November 11, 1998, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     December 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C     February 14, 1999, A. Varga, German Aerospace Center.
C     February 22, 1999, V. Sima, Research Institute for Informatics.
C     February 27, 2000, V. Sima, Research Institute for Informatics.
C     May 26, 2000, A. Varga, German Aerospace Center.
C
C     KEYWORDS
C
C     Balancing, minimal state-space representation, model reduction,
C     multivariable system, singular perturbation approximation,
C     state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDT, LDTI,
     $                  LDWORK, M, N, NR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*), T(LDT,*), TI(LDTI,*)
C     .. Local Scalars ..
      LOGICAL           BAL, DISCR, FIXORD, PACKED
      INTEGER           IERR, IJ, J, K, KTAU, KU, KV, KW, LDW, NMINR,
     $                  NR1, NS, WRKOPT
      DOUBLE PRECISION  ATOL, RCOND, RTOL, SCALEC, SCALEO, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09DD, DGEMM, DGEMV, DGEQRF, DGETRF, DGETRS,
     $                  DLACPY, DORGQR, DSCAL, DTPMV, DTRMM, DTRMV,
     $                  MA02AD, MA02DD, MB03UD, SB03OU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      BAL    = LSAME( JOB,    'B' )
      FIXORD = LSAME( ORDSEL, 'F' )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( BAL .OR. LSAME( JOB, 'N') ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
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
      ELSE IF( TOL2.GT.ZERO .AND. TOL2.GT.TOL1 ) THEN
         INFO = -22
      ELSE IF( LDWORK.LT.MAX( 1, N*( MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -25
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09BX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      RTOL = DBLE( N )*DLAMCH( 'Epsilon' )
C
C     Allocate N*MAX(N,M,P) and N working storage for the matrices U
C     and TAU, respectively.
C
      KU   = 1
      KTAU = KU + N*MAX( N, M, P )
      KW   = KTAU + N
      LDW  = LDWORK - KW + 1
C
C     Copy B in U.
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), N )
C
C     If DISCR = .FALSE., solve for Su the Lyapunov equation
C                                      2
C     A*(Su*Su') + (Su*Su')*A' + scalec *B*B' = 0 .
C
C     If DISCR = .TRUE., solve for Su the Lyapunov equation
C                           2
C     A*(Su*Su')*A' + scalec *B*B' = Su*Su' .
C
C     Workspace:  need   N*(MAX(N,M,P) + 5);
C                 prefer larger.
C
      CALL SB03OU( DISCR, .TRUE., N, M, A, LDA, DWORK(KU), N,
     $             DWORK(KTAU), TI, LDTI, SCALEC, DWORK(KW), LDW, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 1
         RETURN
      ENDIF
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Copy C in U.
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), P )
C
C     If DISCR = .FALSE., solve for Ru the Lyapunov equation
C                                      2
C     A'*(Ru'*Ru) + (Ru'*Ru)*A + scaleo  * C'*C = 0 .
C
C     If DISCR = .TRUE., solve for Ru the Lyapunov equation
C                           2
C     A'*(Ru'*Ru)*A + scaleo  * C'*C = Ru'*Ru .
C
C     Workspace:  need   N*(MAX(N,M,P) + 5);
C                 prefer larger.
C
      CALL SB03OU( DISCR, .FALSE., N, P, A, LDA, DWORK(KU), P,
     $             DWORK(KTAU), T, LDT, SCALEO, DWORK(KW), LDW, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Allocate N*(N+1)/2 (or, if possible, N*N) working storage for the
C     matrix V, a packed (or unpacked) copy of Su, and save Su in V.
C     (The locations for TAU are reused here.)
C
      KV = KTAU
      IF ( LDWORK-KV+1.LT.N*( N + 5 ) ) THEN
         PACKED = .TRUE.
         CALL MA02DD( 'Pack', 'Upper', N, TI, LDTI, DWORK(KV) )
         KW = KV + ( N*( N + 1 ) )/2
      ELSE
         PACKED = .FALSE.
         CALL DLACPY( 'Upper', N, N, TI, LDTI, DWORK(KV), N )
         KW = KV + N*N
      END IF
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
C     Workspace:  need   N*MAX(N,M,P) + N*(N+1)/2 + 5*N;
C                 prefer larger.
C
      CALL MB03UD( 'Vectors', 'Vectors', N, TI, LDTI, DWORK(KU), N, HSV,
     $             DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 2
         RETURN
      ENDIF
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C     Scale singular values.
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
      ATOL = RTOL*HSV(1)
      IF( FIXORD ) THEN
         IF( NR.GT.0 ) THEN
            IF( HSV(NR).LE.ATOL ) THEN
               NR = 0
               IWARN = 1
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
C     Finish if the order of the reduced model is zero.
C
      IF( NR.EQ.0 ) THEN
C
C       Compute only Dr using singular perturbation formulas.
C       Workspace:  need real    4*N;
C                   need integer 2*N.
C
         CALL AB09DD( DICO, N, M, P, NR, A, LDA, B, LDB, C, LDC, D,
     $                LDD, RCOND, IWORK, DWORK, IERR )
         IWORK(1) = 0
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Compute the order of minimal realization as the order of [S1 S2].
C
      NR1 = NR + 1
      NMINR = NR
      IF( NR.LT.N ) THEN
         ATOL = MAX( TOL2, RTOL*HSV(1) )
         DO 40 J = NR1, N
            IF( HSV(J).LE.ATOL ) GO TO 50
            NMINR = NMINR + 1
   40    CONTINUE
   50    CONTINUE
      END IF
C
C     Compute the order of S2.
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
C     Compute  T = | T1 T2 | = Su*| U1 U2 |
C     (with Su packed, if not enough workspace).
C
      CALL MA02AD( 'Full', NMINR, N, TI, LDTI, T, LDT )
      IF ( PACKED ) THEN
         DO 60 J = 1, NMINR
            CALL DTPMV( 'Upper', 'NoTranspose', 'NonUnit', N, DWORK(KV),
     $                  T(1,J), 1 )
   60    CONTINUE
      ELSE
         CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N,
     $               NMINR, ONE, DWORK(KV), N, T, LDT )
      END IF
C
      IF( BAL ) THEN
         IJ = KU
C
C        Square-Root SPA method.
C
C        Compute the truncation matrices for balancing
C                    -1/2            -1/2
C               T1*S1     and TI1'*S1
C
         DO 70 J = 1, NR
            TEMP = ONE/SQRT( HSV(J) )
            CALL DSCAL( N, TEMP, T(1,J), 1 )
            CALL DSCAL( N, TEMP, DWORK(IJ), 1 )
            IJ = IJ + N
   70    CONTINUE
      ELSE
C
C        Balancing-Free SPA method.
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
      CALL AB09DD( DICO, NMINR, M, P, NR, A, LDA, B, LDB, C, LDC, D,
     $             LDD, RCOND, IWORK, DWORK, IERR )
C
      IWORK(1) = NMINR
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09BX ***
      END
