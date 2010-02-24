      SUBROUTINE AB09FD( DICO, JOBCF, FACT, JOBMR, EQUIL, ORDSEL, N, M,
     $                   P, NR, ALPHA, A, LDA, B, LDB, C, LDC, NQ, HSV,
     $                   TOL1, TOL2, IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr) for an original
C     state-space representation (A,B,C) by using either the square-root
C     or the balancing-free square-root Balance & Truncate (B & T)
C     model reduction method in conjunction with stable coprime
C     factorization techniques.
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
C     JOBCF   CHARACTER*1
C             Specifies whether left or right coprime factorization is
C             to be used as follows:
C             = 'L':  use left coprime factorization;
C             = 'R':  use right coprime factorization.
C
C     FACT    CHARACTER*1
C             Specifies the type of coprime factorization to be computed
C             as follows:
C             = 'S':  compute a coprime factorization with prescribed
C                     stability degree ALPHA;
C             = 'I':  compute a coprime factorization with inner
C                     denominator.
C
C     JOBMR   CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root Balance & Truncate method;
C             = 'N':  use the balancing-free square-root
C                     Balance & Truncate method.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
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
C             On entry with ORDSEL = 'F', NR is the desired order of the
C             resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. NR is set as follows:
C             if ORDSEL = 'F', NR is equal to MIN(NR,NQ,NMIN), where NR
C             is the desired order on entry, NQ is the order of the
C             computed coprime factorization of the given system, and
C             NMIN is the order of a minimal realization of the extended
C             system (see METHOD); NMIN is determined as the number of
C             Hankel singular values greater than NQ*EPS*HNORM(Ge),
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH) and HNORM(Ge) is the Hankel norm of the
C             extended system (computed in HSV(1));
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,NQ*EPS*HNORM(Ge)).
C
C     ALPHA   (input) DOUBLE PRECISION
C             If FACT = 'S', the desired stability degree for the
C             factors of the coprime factorization (see SLICOT Library
C             routines SB08ED/SB08FD).
C             ALPHA < 0 for a continuous-time system (DICO = 'C'), and
C             0 <= ALPHA < 1 for a discrete-time system (DICO = 'D').
C             If FACT = 'I', ALPHA is not used.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the reduced
C             order system.
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
C     NQ      (output) INTEGER
C             The order of the computed extended system Ge (see METHOD).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the NQ Hankel singular values of
C             the extended system Ge ordered decreasingly (see METHOD).
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced extended system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(Ge), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(Ge) is the
C             Hankel-norm of the extended system (computed in HSV(1)).
C             The value TOL1 = NQ*EPS*HNORM(Ge) is used by default if
C             TOL1 <= 0 on entry, where EPS is the machine precision
C             (see LAPACK Library Routine DLAMCH).
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The absolute tolerance level below which the elements of
C             B or C are considered zero (used for controllability or
C             observability tests).
C             If the user sets TOL2 <= 0, then an implicitly computed,
C             default tolerance TOLDEF is used:
C             TOLDEF = N*EPS*NORM(C'), if JOBCF = 'L', or
C             TOLDEF = N*EPS*NORM(B),  if JOBCF = 'R',
C             where EPS is the machine precision, and NORM(.) denotes
C             the 1-norm of a matrix.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = PM,        if JOBMR = 'B',
C             LIWORK = MAX(N,PM), if JOBMR = 'N', where
C             PM = P, if JOBCF = 'L',
C             PM = M, if JOBCF = 'R'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,LW1) if JOBCF = 'L' and FACT = 'S',
C             LDWORK >= MAX(1,LW2) if JOBCF = 'L' and FACT = 'I',
C             LDWORK >= MAX(1,LW3) if JOBCF = 'R' and FACT = 'S',
C             LDWORK >= MAX(1,LW4) if JOBCF = 'R' and FACT = 'I', where
C             LW1 = N*(2*MAX(M,P) + P) + MAX(M,P)*(MAX(M,P) + P) +
C                   MAX( N*P+MAX(N*(N+5), 5*P, 4*M), LWR ),
C             LW2 = N*(2*MAX(M,P) + P) + MAX(M,P)*(MAX(M,P) + P) +
C                   MAX( N*P+MAX(N*(N+5), P*(P+2), 4*P, 4*M), LWR ),
C             LW3 = (N+M)*(M+P) + MAX( 5*M, 4*P, LWR ),
C             LW4 = (N+M)*(M+P) + MAX( M*(M+2), 4*M, 4*P, LWR ), and
C             LWR = 2*N*N + N*(MAX(N,M+P)+5) + N*(N+1)/2.
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 10*K+I:
C               I = 1:  with ORDSEL = 'F', the selected order NR is
C                       greater than the order of the computed coprime
C                       factorization of the given system. In this case,
C                       the resulting NR is set automatically to a value
C                       corresponding to the order of a minimal
C                       realization of the system;
C               K > 0:  K violations of the numerical stability
C                       condition occured when computing the coprime
C                       factorization using pole assignment (see SLICOT
C                       Library routines SB08CD/SB08ED, SB08DD/SB08FD).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction of A to a real Schur form failed;
C             = 2:  a failure was detected during the ordering of the
C                   real Schur form of A, or in the iterative process
C                   for reordering the eigenvalues of Z'*(A + H*C)*Z
C                   (or Z'*(A + B*F)*Z) along the diagonal; see SLICOT
C                   Library routines SB08CD/SB08ED (or SB08DD/SB08FD);
C             = 3:  the matrix A has an observable or controllable
C                   eigenvalue on the imaginary axis if DICO = 'C' or
C                   on the unit circle if DICO = 'D';
C             = 4:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t)                               (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system, and let G be the corresponding
C     transfer-function matrix. The subroutine AB09FD determines
C     the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t)                             (2)
C
C     with the transfer-function matrix Gr, by using the
C     balanced-truncation model reduction method in conjunction with
C     a left coprime factorization (LCF) or a right coprime
C     factorization (RCF) technique:
C
C     1. Compute the appropriate stable coprime factorization of G:
C                     -1                   -1
C                G = R  *Q (LCF) or G = Q*R   (RCF).
C
C     2. Perform the model reduction algorithm on the extended system
C                                           ( Q )
C                Ge = ( Q R ) (LCF) or Ge = ( R )  (RCF)
C
C        to obtain a reduced extended system with reduced factors
C                                               ( Qr )
C                Ger = ( Qr Rr ) (LCF) or Ger = ( Rr )  (RCF).
C
C     3. Recover the reduced system from the reduced factors as
C                       -1                       -1
C                Gr = Rr  *Qr (LCF) or Gr = Qr*Rr   (RCF).
C
C     The approximation error for the extended system satisfies
C
C        HSV(NR) <= INFNORM(Ge-Ger) <= 2*[HSV(NR+1) + ... + HSV(NQ)],
C
C     where INFNORM(G) is the infinity-norm of G.
C
C     If JOBMR = 'B', the square-root Balance & Truncate method of [1]
C     is used for model reduction.
C     If JOBMR = 'N', the balancing-free square-root version of the
C     Balance & Truncate method [2] is used for model reduction.
C
C     If FACT = 'S', the stable coprime factorization with prescribed
C     stability degree ALPHA is computed by using the algorithm of [3].
C     If FACT = 'I', the stable coprime factorization with inner
C     denominator is computed by using the algorithm of [4].
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
C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), Vol. 2,
C         pp. 42-46, 1991.
C
C     [3] Varga A.
C         Coprime factors model reduction method based on square-root
C         balancing-free techniques.
C         System Analysis, Modelling and Simulation, Vol. 11,
C         pp. 303-311, 1993.
C
C     [4] Varga A.
C         A Schur method for computing coprime factorizations with
C         inner denominators and applications in model reduction.
C         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993.
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
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, August 1998.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven.
C
C     KEYWORDS
C
C     Balancing, coprime factorization, minimal realization,
C     model reduction, multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  C100, ONE, ZERO
      PARAMETER         ( C100 = 100.0D0, ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, FACT, JOBCF, JOBMR, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDWORK, M, N, NQ,
     $                  NR, P
      DOUBLE PRECISION  ALPHA, TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           DISCR, FIXORD, LEFT, STABD
      INTEGER           IERR, IWARNK, KB, KBR, KBT, KC, KCR, KD, KDR,
     $                  KDT, KT, KTI, KW, LW1, LW2, LW3, LW4, LWR,
     $                  MAXMP, MP, NDR, PM, WRKOPT
      DOUBLE PRECISION  MAXRED
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09AX, DLACPY, DLASET, SB08CD, SB08DD, SB08ED,
     $                  SB08FD, SB08GD, SB08HD, TB01ID, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
      LEFT   = LSAME( JOBCF,  'L' )
      STABD  = LSAME( FACT,   'S' )
      MAXMP  = MAX( M, P )
C
      LWR = 2*N*N + N*( MAX( N, M + P ) + 5 ) + ( N*( N + 1 ) )/2
      LW1 = N*( 2*MAXMP + P ) + MAXMP*( MAXMP + P )
      LW2 = LW1 +
     $      MAX( N*P + MAX( N*( N + 5 ), P*( P+2 ), 4*P, 4*M ), LWR )
      LW1 = LW1 + MAX( N*P + MAX( N*( N + 5 ), 5*P, 4*M ), LWR )
      LW3 = ( N + M )*( M + P ) + MAX( 5*M, 4*P, LWR )
      LW4 = ( N + M )*( M + P ) + MAX( M*( M + 2 ), 4*M, 4*P, LWR )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LEFT .OR. LSAME( JOBCF, 'R' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( STABD .OR. LSAME( FACT, 'I' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LSAME( JOBMR, 'B' ) .OR.
     $                 LSAME( JOBMR, 'N' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -5
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -6
      ELSE IF( STABD .AND. ( ( .NOT.DISCR .AND. ALPHA.GE.ZERO ) .OR.
     $       ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GE.ONE ) ) ) )
     $      THEN
         INFO = -7
      ELSE IF( N.LT.0 ) THEN
         INFO = -8
      ELSE IF( M.LT.0 ) THEN
         INFO = -9
      ELSE IF( P.LT.0 ) THEN
         INFO = -10
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -11
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -17
      ELSE IF( ( LDWORK.LT.1 ) .OR.
     $   (      STABD .AND.      LEFT .AND. LDWORK.LT.LW1 ) .OR.
     $   ( .NOT.STABD .AND.      LEFT .AND. LDWORK.LT.LW2 ) .OR.
     $   (      STABD .AND. .NOT.LEFT .AND. LDWORK.LT.LW3 ) .OR.
     $   ( .NOT.STABD .AND. .NOT.LEFT .AND. LDWORK.LT.LW4 ) ) THEN
         INFO = -24
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09FD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 .OR. ( FIXORD .AND. NR.EQ.0 ) ) THEN
         NR = 0
         NQ = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C
         MAXRED = C100
         CALL TB01ID( 'A', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
      END IF
C
C     Perform the coprime factor model reduction procedure.
C
      KD = 1
      IF( LEFT ) THEN
C                           -1
C        Compute a LCF G = R  *Q.
C
         MP  = M + P
         KDR = KD  + MAXMP*MAXMP
         KC  = KDR + MAXMP*P
         KB  = KC  + MAXMP*N
         KBR = KB  + N*MAXMP
         KW  = KBR + N*P
         LWR = LDWORK - KW + 1
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KB), N )
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KC), MAXMP )
         CALL DLASET( 'Full', P, M, ZERO, ZERO, DWORK(KD), MAXMP )
C
         IF( STABD ) THEN
C
C           Compute a LCF with prescribed stability degree.
C
C           Workspace needed:      N*(2*MAX(M,P)+P) +
C                                  MAX(M,P)*(MAX(M,P)+P);
C           Additional workspace:  need   N*P+MAX(N*(N+5),5*P,4*M);
C                                  prefer larger.
C
            CALL SB08ED( DICO, N, M, P, ALPHA, A, LDA, DWORK(KB), N,
     $                   DWORK(KC), MAXMP, DWORK(KD), MAXMP, NQ, NDR,
     $                   DWORK(KBR), N, DWORK(KDR), MAXMP, TOL2,
     $                   DWORK(KW), LWR, IWARN, INFO )
         ELSE
C
C           Compute a LCF with inner denominator.
C
C           Workspace needed:      N*(2*MAX(M,P)+P) +
C                                  MAX(M,P)*(MAX(M,P)+P);
C           Additional workspace:  need   N*P +
C                                         MAX(N*(N+5),P*(P+2),4*P,4*M).
C                                  prefer larger;
C
            CALL SB08CD( DICO, N, M, P, A, LDA, DWORK(KB), N,
     $                   DWORK(KC), MAXMP, DWORK(KD), MAXMP, NQ, NDR,
     $                   DWORK(KBR), N, DWORK(KDR), MAXMP, TOL2,
     $                   DWORK(KW), LWR, IWARN, INFO )
         END IF
C
         IWARN = 10*IWARN
         IF( INFO.NE.0 )
     $      RETURN
C
         WRKOPT = INT( DWORK(KW) ) + KW - 1
C
         IF( NQ.EQ.0 ) THEN
            NR = 0
            DWORK(1) = WRKOPT
            RETURN
         END IF
C
         IF( MAXMP.GT.M ) THEN
C
C           Form the matrices ( BQ, BR ) and ( DQ, DR ) in consecutive
C           columns (see SLICOT Library routines SB08CD/SB08ED).
C
            KBT = KBR
            KBR = KB + N*M
            KDT = KDR
            KDR = KD + MAXMP*M
            CALL DLACPY( 'Full', NQ, P, DWORK(KBT), N, DWORK(KBR), N )
            CALL DLACPY( 'Full', P, P, DWORK(KDT), MAXMP, DWORK(KDR),
     $                   MAXMP )
         END IF
C
C        Perform model reduction on ( Q, R ) to determine ( Qr, Rr ).
C
C        Workspace needed:      N*(2*MAX(M,P)+P) +
C                               MAX(M,P)*(MAX(M,P)+P) + 2*N*N;
C        Additional workspace:  need   N*(MAX(N,M+P)+5) + N*(N+1)/2;
C                               prefer larger.
C
         KT  = KW
         KTI = KT  + NQ*NQ
         KW  = KTI + NQ*NQ
         CALL AB09AX( DICO, JOBMR, ORDSEL, NQ, MP, P, NR, A, LDA,
     $                DWORK(KB), N, DWORK(KC), MAXMP, HSV, DWORK(KT),
     $                N, DWORK(KTI), N, TOL1, IWORK, DWORK(KW),
     $                LDWORK-KW+1, IWARNK, IERR )
C
         IWARN = IWARN + IWARNK
         IF( IERR.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C                                                -1
C        Compute the reduced order system Gr = Rr  *Qr.
C
C        Workspace needed:      N*(2*MAX(M,P)+P) +
C                               MAX(M,P)*(MAX(M,P)+P);
C        Additional workspace:  need   4*P.
C
         KW = KT
         CALL SB08GD( NR, M, P, A, LDA, DWORK(KB), N, DWORK(KC), MAXMP,
     $                DWORK(KD), MAXMP, DWORK(KBR), N, DWORK(KDR),
     $                MAXMP, IWORK, DWORK(KW), INFO )
C
C        Copy the reduced system matrices Br and Cr to B and C.
C
         CALL DLACPY( 'Full', NR, M, DWORK(KB), N, B, LDB )
         CALL DLACPY( 'Full', P, NR, DWORK(KC), MAXMP, C, LDC )
C
      ELSE
C                             -1
C        Compute a RCF G = Q*R  .
C
         PM  = P  + M
         KDR = KD + P
         KC  = KD + PM*M
         KCR = KC + P
         KW  = KC + PM*N
         LWR = LDWORK - KW + 1
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KC), PM )
         CALL DLASET( 'Full', P, M, ZERO, ZERO, DWORK(KD), PM )
C
         IF( STABD ) THEN
C
C           Compute a RCF with prescribed stability degree.
C
C           Workspace needed:      (N+M)*(M+P);
C           Additional workspace:  need   MAX( N*(N+5), 5*M, 4*P );
C                                  prefer larger.
C
            CALL SB08FD( DICO, N, M, P, ALPHA, A, LDA, B, LDB,
     $                   DWORK(KC), PM, DWORK(KD), PM, NQ, NDR,
     $                   DWORK(KCR), PM, DWORK(KDR), PM, TOL2,
     $                   DWORK(KW), LWR, IWARN, INFO )
         ELSE
C
C           Compute a RCF with inner denominator.
C
C           Workspace needed:      (N+M)*(M+P);
C           Additional workspace:  need   MAX(N*(N+5),M*(M+2),4*M,4*P);
C                                  prefer larger.
C
            CALL SB08DD( DICO, N, M, P, A, LDA, B, LDB,
     $                   DWORK(KC), PM, DWORK(KD), PM, NQ, NDR,
     $                   DWORK(KCR), PM, DWORK(KDR), PM, TOL2,
     $                   DWORK(KW), LWR, IWARN, INFO )
         END IF
C
         IWARN = 10*IWARN
         IF( INFO.NE.0 )
     $      RETURN
C
         WRKOPT = INT( DWORK(KW) ) + KW - 1
C
         IF( NQ.EQ.0 ) THEN
            NR = 0
            DWORK(1) = WRKOPT
            RETURN
         END IF
C                                   ( Q )              ( Qr )
C        Perform model reduction on ( R ) to determine ( Rr ).
C
C        Workspace needed:      (N+M)*(M+P) + 2*N*N;
C        Additional workspace:  need   N*(MAX(N,M+P)+5) + N*(N+1)/2;
C                               prefer larger.
C
         KT  = KW
         KTI = KT  + NQ*NQ
         KW  = KTI + NQ*NQ
         CALL AB09AX( DICO, JOBMR, ORDSEL, NQ, M, PM, NR, A, LDA, B,
     $                LDB, DWORK(KC), PM, HSV, DWORK(KT), N, DWORK(KTI),
     $                N, TOL1, IWORK, DWORK(KW), LDWORK-KW+1, IWARNK,
     $                IERR )
C
         IWARN = IWARN + IWARNK
         IF( IERR.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C                                                   -1
C        Compute the reduced order system Gr = Qr*Rr  .
C
C        Workspace needed:      (N+M)*(M+P);
C        Additional workspace:  need 4*M.
C
         KW = KT
         CALL SB08HD( NR, M, P, A, LDA, B, LDB, DWORK(KC), PM,
     $                DWORK(KD), PM, DWORK(KCR), PM, DWORK(KDR), PM,
     $                IWORK, DWORK(KW), INFO )
C
C        Copy the reduced system matrix Cr to C.
C
         CALL DLACPY( 'Full', P, NR, DWORK(KC), PM, C, LDC )
      END IF
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09FD ***
      END
