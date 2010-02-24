      SUBROUTINE AB09AX( DICO, JOB, ORDSEL, N, M, P, NR, A, LDA, B, LDB,
     $                   C, LDC, HSV, T, LDT, TI, LDTI, TOL, IWORK,
     $                   DWORK, LDWORK, IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr) for a stable original
C     state-space representation (A,B,C) by using either the square-root
C     or the balancing-free square-root Balance & Truncate model
C     reduction method. The state dynamics matrix A of the original
C     system is an upper quasi-triangular matrix in real Schur canonical
C     form. The matrices of the reduced order system are computed using
C     the truncation formulas:
C
C          Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T .
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
C             = 'N':  use the balancing-free square-root
C                     Balance & Truncate method.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL.
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
C             On entry with ORDSEL = 'F', NR is the desired order of the
C             resulting reduced order system.  0 <= NR <= N.
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
C             singular values greater than MAX(TOL,N*EPS*HNORM(A,B,C)).
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
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the Hankel singular values of
C             the original system ordered decreasingly. HSV(1) is the
C             Hankel norm of the system.
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,N)
C             If INFO = 0 and NR > 0, the leading N-by-NR part of this
C             array contains the right truncation matrix T.
C
C     LDT     INTEGER
C             The leading dimension of array T.  LDT >= MAX(1,N).
C
C     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N)
C             If INFO = 0 and NR > 0, the leading NR-by-N part of this
C             array contains the left truncation matrix TI.
C
C     LDTI    INTEGER
C             The leading dimension of array TI.  LDTI >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If ORDSEL = 'A', TOL contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL = c*HNORM(A,B,C), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(A,B,C) is the
C             Hankel-norm of the given system (computed in HSV(1)).
C             For computing a minimal realization, the recommended
C             value is TOL = N*EPS*HNORM(A,B,C), where EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH).
C             This value is used by default if TOL <= 0 on entry.
C             If ORDSEL = 'F', the value of TOL is ignored.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = 0, if JOB = 'B', or
C             LIWORK = N, if JOB = 'N'.
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
C          y(t)    = Cx(t)                               (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09AX determines for
C     the given system (1), the matrices of a reduced NR order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t)                             (2)
C
C     such that
C
C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     If JOB = 'B', the square-root Balance & Truncate method of [1]
C     is used and, for DICO = 'C', the resulting model is balanced.
C     By setting TOL <= 0, the routine can be used to compute balanced
C     minimal state-space realizations of stable systems.
C
C     If JOB = 'N', the balancing-free square-root version of the
C     Balance & Truncate method [2] is used.
C     By setting TOL <= 0, the routine can be used to compute minimal
C     state-space realizations of stable systems.
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
C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.),
C         Vol. 2, pp. 42-46.
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
C     Based on the RASP routines SRBT1 and SRBFT1.
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
C
C     KEYWORDS
C
C     Balancing, minimal state-space representation, model reduction,
C     multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDT, LDTI, LDWORK,
     $                  M, N, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), HSV(*),
     $                  T(LDT,*), TI(LDTI,*)
C     .. Local Scalars ..
      LOGICAL           BAL, DISCR, FIXORD, PACKED
      INTEGER           IERR, IJ, J, K, KTAU, KU, KV, KW, LDW, WRKOPT
      DOUBLE PRECISION  ATOL, RTOL, SCALEC, SCALEO, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGEMV, DGEQRF, DGETRF, DGETRS, DLACPY,
     $                  DORGQR, DSCAL, DTPMV, DTRMM, DTRMV, MA02AD,
     $                  MA02DD, MB03UD, SB03OU, XERBLA
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
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -16
      ELSE IF( LDTI.LT.MAX( 1, N ) ) THEN
         INFO = -18
      ELSE IF( LDWORK.LT.MAX( 1, N*( MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -22
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09AX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 .OR. ( FIXORD .AND. NR.EQ.0 ) ) THEN
         NR = 0
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
C     S = diag(S1,S2),  U = [U1,U2] (U' in TI) and V = [V1,V2] (in U).
C
C     Compute the order of reduced system, as the order of S1.
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
         ATOL = MAX( TOL, ATOL )
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
      IF( NR.EQ.0 ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Compute the truncation matrices.
C
C     Compute TI' =  Ru'*V1 in U.
C
      CALL DTRMM( 'Left', 'Upper', 'Transpose', 'NonUnit', N, NR, ONE,
     $             T, LDT, DWORK(KU), N )
C
C     Compute T = Su*U1 (with Su packed, if not enough workspace).
C
      CALL MA02AD( 'Full', NR, N, TI, LDTI, T, LDT )
      IF ( PACKED ) THEN
         DO 40 J = 1, NR
            CALL DTPMV( 'Upper', 'NoTranspose', 'NonUnit', N, DWORK(KV),
     $                  T(1,J), 1 )
   40    CONTINUE
      ELSE
         CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N, NR,
     $               ONE, DWORK(KV), N, T, LDT )
      END IF
C
      IF( BAL ) THEN
         IJ = KU
C
C        Square-Root B & T method.
C
C        Compute the truncation matrices for balancing
C                    -1/2           -1/2
C                T*S1     and TI'*S1
C
         DO 50 J = 1, NR
            TEMP = ONE/SQRT( HSV(J) )
            CALL DSCAL( N, TEMP, T(1,J), 1 )
            CALL DSCAL( N, TEMP, DWORK(IJ), 1 )
            IJ = IJ + N
   50    CONTINUE
      ELSE
C
C        Balancing-Free B & T method.
C
C        Compute orthogonal bases for the images of matrices T and TI'.
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
      END IF
C
C     Transpose TI' to obtain TI.
C
      CALL MA02AD( 'Full', N, NR, DWORK(KU), N, TI, LDTI )
C
      IF( .NOT.BAL ) THEN
C                      -1
C        Compute (TI*T)  *TI in TI.
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, NR, N, ONE, TI,
     $               LDTI, T, LDT, ZERO, DWORK(KU), N )
         CALL DGETRF( NR, NR, DWORK(KU), N, IWORK, IERR )
         CALL DGETRS( 'NoTranspose', NR, N, DWORK(KU), N, IWORK, TI,
     $                LDTI, IERR )
      END IF
C
C     Compute TI*A*T (A is in RSF).
C
      IJ = KU
      DO 60 J = 1, N
         K = MIN( J+1, N )
         CALL DGEMV( 'NoTranspose', NR, K, ONE, TI, LDTI, A(1,J), 1,
     $               ZERO, DWORK(IJ), 1 )
         IJ = IJ + N
   60 CONTINUE
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, NR, N, ONE,
     $            DWORK(KU), N, T, LDT, ZERO, A, LDA )
C
C     Compute TI*B and C*T.
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), N )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, M, N, ONE, TI, LDTI,
     $            DWORK(KU), N, ZERO, B, LDB )
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), P )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NR, N, ONE,
     $            DWORK(KU), P, T, LDT, ZERO, C, LDC )
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09AX ***
      END
