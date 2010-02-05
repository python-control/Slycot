      SUBROUTINE AB09JV( JOB, DICO, JOBEV, STBCHK, N, M, P, NV, PV,
     $                   A, LDA, B, LDB, C, LDC, D, LDD, AV, LDAV,
     $                   EV, LDEV, BV, LDBV, CV, LDCV, DV, LDDV, IWORK,
     $                   DWORK, LDWORK, INFO )
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
C     To construct a state-space representation (A,BS,CS,DS) of the
C     projection of V*G or conj(V)*G containing the poles of G, from the
C     state-space representations (A,B,C,D) and (AV-lambda*EV,BV,CV,DV),
C     of the transfer-function matrices G and V, respectively.
C     G is assumed to be a stable transfer-function matrix and
C     the state matrix A must be in a real Schur form.
C     When computing the stable projection of V*G, it is assumed
C     that G and V have completely distinct poles.
C     When computing the stable projection of conj(V)*G, it is assumed
C     that G and conj(V) have completely distinct poles.
C
C     Note: For a transfer-function matrix G, conj(G) denotes the
C     conjugate of G given by G'(-s) for a continuous-time system or
C     G'(1/z) for a discrete-time system.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies the projection to be computed as follows:
C             = 'V':  compute the projection of V*G containing
C                     the poles of G;
C             = 'C':  compute the projection of conj(V)*G containing
C                     the poles of G.
C
C     DICO    CHARACTER*1
C             Specifies the type of the systems as follows:
C             = 'C':  G and V are continuous-time systems;
C             = 'D':  G and V are discrete-time systems.
C
C     JOBEV   CHARACTER*1
C             Specifies whether EV is a general square or an identity
C             matrix as follows:
C             = 'G':  EV is a general square matrix;
C             = 'I':  EV is the identity matrix.
C
C     STBCHK  CHARACTER*1
C             Specifies whether stability/antistability of V is to be
C             checked as follows:
C             = 'C':  check stability if JOB = 'C' or antistability if
C                     JOB = 'V';
C             = 'N':  do not check stability or antistability.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the state vector of the system with
C             the transfer-function matrix G.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of the input vector of the system with
C             the transfer-function matrix G.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of the output vector of the system with the
C             transfer-function matrix G, and also the dimension of
C             the input vector if JOB = 'V', or of the output vector
C             if JOB = 'C', of the system with the transfer-function
C             matrix V.  P >= 0.
C
C     NV      (input) INTEGER
C             The dimension of the state vector of the system with
C             the transfer-function matrix V.  NV >= 0.
C
C     PV      (input) INTEGER
C             The dimension of the output vector, if JOB = 'V', or
C             of the input vector, if JOB = 'C', of the system with
C             the transfer-function matrix V.  PV >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the system with the transfer-function
C             matrix G in a real Schur form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain
C             the input/state matrix B of the system with the
C             transfer-function matrix G. The matrix BS is equal to B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading PV-by-N part of this
C             array contains the output matrix CS of the projection of
C             V*G, if JOB = 'V', or of conj(V)*G, if JOB = 'C'.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,P,PV).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the feedthrough matrix D of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading PV-by-M part of
C             this array contains the feedthrough matrix DS of the
C             projection of V*G, if JOB = 'V', or of conj(V)*G,
C             if JOB = 'C'.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,P,PV).
C
C     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV)
C             On entry, the leading NV-by-NV part of this array must
C             contain the state matrix AV of the system with the
C             transfer-function matrix V.
C             On exit, if INFO = 0, the leading NV-by-NV part of this
C             array contains a condensed matrix as follows:
C             if JOBEV = 'I', it contains the real Schur form of AV;
C             if JOBEV = 'G' and JOB = 'V', it contains a quasi-upper
C             triangular matrix representing the real Schur matrix
C             in the real generalized Schur form of the pair (AV,EV);
C             if JOBEV = 'G', JOB = 'C' and DICO = 'C', it contains a
C             quasi-upper triangular matrix corresponding to the
C             generalized real Schur form of the pair (AV',EV');
C             if JOBEV = 'G', JOB = 'C' and DICO = 'D', it contains an
C             upper triangular matrix corresponding to the generalized
C             real Schur form of the pair (EV',AV').
C
C     LDAV    INTEGER
C             The leading dimension of the array AV.  LDAV >= MAX(1,NV).
C
C     EV      (input/output) DOUBLE PRECISION array, dimension (LDEV,NV)
C             On entry, if JOBEV = 'G', the leading NV-by-NV part of
C             this array must contain the descriptor matrix EV of the
C             system with the transfer-function matrix V.
C             If JOBEV = 'I', EV is assumed to be an identity matrix
C             and is not referenced.
C             On exit, if INFO = 0 and JOBEV = 'G', the leading NV-by-NV
C             part of this array contains a condensed matrix as follows:
C             if JOB = 'V', it contains an upper triangular matrix
C             corresponding to the real generalized Schur form of the
C             pair (AV,EV);
C             if JOB = 'C' and DICO = 'C', it contains an upper
C             triangular matrix corresponding to the generalized real
C             Schur form of the pair (AV',EV');
C             if JOB = 'C' and DICO = 'D', it contains a quasi-upper
C             triangular matrix corresponding to the generalized
C             real Schur form of the pair (EV',AV').
C
C     LDEV    INTEGER
C             The leading dimension of the array EV.
C             LDEV >= MAX(1,NV), if JOBEV = 'G';
C             LDEV >= 1,         if JOBEV = 'I'.
C
C     BV      (input/output) DOUBLE PRECISION array,
C             dimension (LDBV,MBV), where MBV = P, if JOB = 'V', and
C             MBV = PV, if JOB = 'C'.
C             On entry, the leading NV-by-MBV part of this array must
C             contain the input matrix BV of the system with the
C             transfer-function matrix V.
C             On exit, if INFO = 0, the leading NV-by-MBV part of this
C             array contains Q'*BV, where Q is the orthogonal matrix
C             that reduces AV to the real Schur form or the left
C             orthogonal matrix used to reduce the pair (AV,EV),
C             (AV',EV') or (EV',AV') to the generalized real Schur form.
C
C     LDBV    INTEGER
C             The leading dimension of the array BV.  LDBV >= MAX(1,NV).
C
C     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV)
C             On entry, the leading PCV-by-NV part of this array must
C             contain the output matrix CV of the system with the
C             transfer-function matrix V, where PCV = PV, if JOB = 'V',
C             or PCV = P, if JOB = 'C'.
C             On exit, if INFO = 0, the leading PCV-by-NV part of this
C             array contains CV*Q, where Q is the orthogonal matrix that
C             reduces AV to the real Schur form, or CV*Z, where Z is the
C             right orthogonal matrix used to reduce the pair (AV,EV),
C             (AV',EV') or (EV',AV') to the generalized real Schur form.
C
C     LDCV    INTEGER
C             The leading dimension of the array CV.
C             LDCV >= MAX(1,PV) if JOB = 'V';
C             LDCV >= MAX(1,P)  if JOB = 'C'.
C
C     DV      (input) DOUBLE PRECISION array,
C             dimension (LDDV,MBV), where MBV = P, if JOB = 'V', and
C             MBV = PV, if JOB = 'C'.
C             The leading PCV-by-MBV part of this array must contain
C             the feedthrough matrix DV of the system with the
C             transfer-function matrix V, where PCV = PV, if JOB = 'V',
C             or PCV = P, if JOB = 'C'.
C
C     LDDV    INTEGER
C             The leading dimension of the array DV.
C             LDDV >= MAX(1,PV) if JOB = 'V';
C             LDDV >= MAX(1,P)  if JOB = 'C'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK =   0,    if JOBEV = 'I';
C             LIWORK = NV+N+6, if JOBEV = 'G'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= LW1, if JOBEV = 'I',
C             LDWORK >= LW2, if JOBEV = 'G', where
C               LW1 = MAX( 1, NV*(NV+5), NV*N + MAX( a, PV*N, PV*M ) )
C                     a = 0,    if DICO = 'C' or  JOB = 'V',
C                     a = 2*NV, if DICO = 'D' and JOB = 'C';
C               LW2 = MAX( 2*NV*NV + MAX( 11*NV+16, P*NV, PV*NV ),
C                          NV*N + MAX( NV*N+N*N, PV*N, PV*M ) ).
C             For good performance, LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  the reduction of the pair (AV,EV) to the real
C                    generalized Schur form failed (JOBEV = 'G'),
C                    or the reduction of the matrix AV to the real
C                    Schur form failed (JOBEV = 'I);
C             =  2:  the solution of the Sylvester equation failed
C                    because the matrix A and the pencil AV-lambda*EV
C                    have common eigenvalues (if JOB = 'V'), or the
C                    pencil -AV-lambda*EV and A have common eigenvalues
C                    (if JOB = 'C' and DICO = 'C'), or the pencil
C                    AV-lambda*EV has an eigenvalue which is the
C                    reciprocal of one of eigenvalues of A
C                    (if JOB = 'C' and DICO = 'D');
C             =  3:  the solution of the Sylvester equation failed
C                    because the matrices A and AV have common
C                    eigenvalues (if JOB = 'V'), or the matrices A
C                    and -AV have common eigenvalues (if JOB = 'C' and
C                    DICO = 'C'), or the matrix A has an eigenvalue
C                    which is the reciprocal of one of eigenvalues of AV
C                    (if JOB = 'C' and DICO = 'D');
C             =  4:  JOB = 'V' and the pair (AV,EV) has not completely
C                    unstable generalized eigenvalues, or JOB = 'C' and
C                    the pair (AV,EV) has not completely stable
C                    generalized eigenvalues.
C
C     METHOD
C
C     If JOB = 'V', the matrices of the stable projection of V*G are
C     computed as
C
C       BS = B,  CS = CV*X + DV*C,  DS = DV*D,
C
C     where X satisfies the generalized Sylvester equation
C
C       AV*X - EV*X*A + BV*C = 0.
C
C     If JOB = 'C', the matrices of the stable projection of conj(V)*G
C     are computed using the following formulas:
C
C     - for a continuous-time system, the matrices BS, CS and DS of
C       the stable projection are computed as
C
C         BS = B,  CS = BV'*X + DV'*C,  DS = DV'*D,
C
C       where X satisfies the generalized Sylvester equation
C
C         AV'*X + EV'*X*A + CV'*C = 0.
C
C     - for a discrete-time system, the matrices BS, CS and DS of
C       the stable projection are computed as
C
C         BS = B,  CS = BV'*X*A + DV'*C,  DS = DV'*D + BV'*X*B,
C
C       where X satisfies the generalized Sylvester equation
C
C         EV'*X - AV'*X*A = CV'*C.
C
C     REFERENCES
C
C     [1] Varga, A.
C         Efficient and numerically reliable implementation of the
C         frequency-weighted Hankel-norm approximation model reduction
C         approach.
C         Proc. 2001 ECC, Porto, Portugal, 2001.
C
C     [2] Zhou, K.
C         Frequency-weighted H-infinity norm and optimal Hankel norm
C         model reduction.
C         IEEE Trans. Autom. Control, vol. 40, pp. 1687-1699, 1995.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on numerically stable algorithms.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000.
C     D. Sima, University of Bucharest, March 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     REVISIONS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001.
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, Nov. 2003.
C
C     KEYWORDS
C
C     Frequency weighting, model reduction, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, JOB, JOBEV, STBCHK
      INTEGER           INFO, LDA, LDAV, LDB, LDBV, LDC, LDCV,
     $                  LDD, LDDV, LDEV, LDWORK, M, N, NV, P, PV
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AV(LDAV,*), B(LDB,*), BV(LDBV,*),
     $                  C(LDC,*), CV(LDCV,*), D(LDD,*), DV(LDDV,*),
     $                  DWORK(*), EV(LDEV,*)
C     .. Local Scalars ..
      CHARACTER*1       EVTYPE, STDOM
      LOGICAL           CONJS, DISCR, STABCK, UNITEV
      DOUBLE PRECISION  ALPHA, DIF, SCALE, TOLINF, WORK
      INTEGER           I, IA, IERR, KAI, KAR, KB, KC, KE, KF, KQ, KW,
     $                  KZ, LDW, LDWN, LW, SDIM
C     .. Local Arrays ..
      LOGICAL           BWORK(1)
C     .. External Functions ..
      LOGICAL           DELCTG, LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DELCTG, DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09JX, DGEMM, DGGES, DLACPY, DLASET, DSWAP,
     $                  DTGSYL, DTRSYL, SB04PY, TB01WD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, SQRT
C
C     .. Executable Statements ..
C
      CONJS  = LSAME( JOB,    'C' )
      DISCR  = LSAME( DICO,   'D' )
      UNITEV = LSAME( JOBEV,  'I' )
      STABCK = LSAME( STBCHK, 'C' )
C
      INFO = 0
      IF( UNITEV ) THEN
         IF ( DISCR .AND. CONJS ) THEN
            IA = 2*NV
         ELSE
            IA = 0
         END IF
         LW = MAX( 1, NV*( NV + 5 ), NV*N + MAX( IA, PV*N, PV*M ) )
      ELSE
         LW = MAX( 2*NV*NV + MAX( 11*NV+16, P*NV, PV*NV ),
     $             NV*N + MAX( NV*N + N*N, PV*N, PV*M ) )
      END IF
C
C     Test the input scalar arguments.
C
      LDWN = MAX( 1, N )
      LDW  = MAX( 1, NV )
      IF( .NOT. ( LSAME( JOB, 'V' ) .OR. CONJS ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( DICO,  'C' ) .OR. DISCR  ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( JOBEV, 'G' ) .OR. UNITEV ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LSAME( STBCHK, 'N' ) .OR. STABCK ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( NV.LT.0 ) THEN
         INFO = -8
      ELSE IF( PV.LT.0 ) THEN
         INFO = -9
      ELSE IF( LDA.LT.LDWN ) THEN
         INFO = -11
      ELSE IF( LDB.LT.LDWN ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P, PV ) ) THEN
         INFO = -15
      ELSE IF( LDD.LT.MAX( 1, P, PV ) ) THEN
         INFO = -17
      ELSE IF( LDAV.LT.LDW ) THEN
         INFO = -19
      ELSE IF( LDEV.LT.1 .OR. ( .NOT.UNITEV .AND. LDEV.LT.NV ) ) THEN
         INFO = -21
      ELSE IF( LDBV.LT.LDW ) THEN
         INFO = -23
      ELSE IF( ( .NOT.CONJS .AND. LDCV.LT.MAX( 1, PV ) ) .OR.
     $         (      CONJS .AND. LDCV.LT.MAX( 1, P  ) ) ) THEN
         INFO = -25
      ELSE IF( ( .NOT.CONJS .AND. LDDV.LT.MAX( 1, PV ) ) .OR.
     $         (      CONJS .AND. LDDV.LT.MAX( 1, P  ) ) ) THEN
         INFO = -27
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -30
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09JV', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( P.EQ.0 .OR. PV.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Set options for stability/antistability checking.
C
      IF( DISCR ) THEN
         ALPHA = ONE
      ELSE
         ALPHA = ZERO
      END IF
C
      WORK = ONE
      TOLINF = DLAMCH( 'Epsilon' )
C
      IF( UNITEV ) THEN
C
C        EV is the identity matrix.
C
         IF( NV.GT.0 ) THEN
C
C           Reduce AV to the real Schur form using an orthogonal
C           similarity transformation AV <- Q'*AV*Q and apply the
C           transformation to BV and CV: BV <- Q'*BV and CV <- CV*Q.
C
C           Workspace needed:  NV*(NV+5);
C                              prefer larger.
C
            KW = NV*( NV + 2 ) + 1
            IF( CONJS ) THEN
               STDOM = 'S'
               ALPHA = ALPHA + SQRT( TOLINF )
               CALL TB01WD( NV, PV, P, AV, LDAV, BV, LDBV, CV, LDCV,
     $                      DWORK(2*NV+1), NV, DWORK, DWORK(NV+1),
     $                      DWORK(KW), LDWORK-KW+1, IERR )
            ELSE
               STDOM = 'U'
               ALPHA = ALPHA - SQRT( TOLINF )
               CALL TB01WD( NV, P, PV, AV, LDAV, BV, LDBV, CV, LDCV,
     $                      DWORK(2*NV+1), NV, DWORK, DWORK(NV+1),
     $                      DWORK(KW), LDWORK-KW+1, IERR )
            END IF
            IF( IERR.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            IF( STABCK ) THEN
C
C              Check stability/antistability of eigenvalues of AV.
C
               CALL AB09JX( DICO, STDOM, 'S', NV, ALPHA, DWORK,
     $                      DWORK(NV+1), DWORK, TOLINF, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 4
                  RETURN
               END IF
            END IF
C
            WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
         END IF
C
         KW = NV*N + 1
         IF( CONJS ) THEN
C
C           Compute the projection of conj(V)*G.
C
C           Total workspace needed:  NV*N + MAX( a, PV*N, PV*M ), where
C                                    a = 0,    if DICO = 'C',
C                                    a = 2*NV, if DICO = 'D'.
C
C           Compute -CV'*C.
C           Workspace needed: NV*N.
C
            CALL DGEMM( 'T', 'N', NV, N, P, -ONE, CV, LDCV, C, LDC,
     $                  ZERO, DWORK, LDW )
C
            IF( DISCR ) THEN
C
C              Compute X and SCALE satisfying
C
C              AV'*X*A - X = -SCALE*CV'*C.
C
C              Additional workspace needed: 2*NV.
C
               CALL SB04PY( 'T', 'N', -1, NV, N, AV, LDAV, A, LDA,
     $                      DWORK, LDW, SCALE, DWORK(KW), IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 3
                  RETURN
               END IF
C
C              Construct CS = DV'*C + BV'*X*A/SCALE,
C                        DS = DV'*D + BV'*X*B/SCALE.
C
C              Additional workspace needed: MAX( PV*N, PV*M ).
C
C              C <- DV'*C.
C
               CALL DGEMM( 'T', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
C
C              D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
C
C              C <- C + BV'*X*A/SCALE.
C
               CALL DGEMM( 'T', 'N', PV, N, NV, ONE / SCALE, BV, LDBV,
     $                     DWORK, LDW, ZERO, DWORK(KW), PV )
               CALL DGEMM( 'N', 'N', PV, N, N, ONE, DWORK(KW), PV,
     $                     A, LDA, ONE, C, LDC )
C
C              D <- D + BV'*X*B/SCALE.
C
               CALL DGEMM( 'N', 'N', PV, M, N, ONE, DWORK(KW), PV,
     $                     B, LDB, ONE, D, LDD )
            ELSE
C
C              Compute X and SCALE satisfying
C
C              AV'*X + X*A + SCALE*CV'*C = 0.
C
               IF( N.GT.0 ) THEN
                  CALL DTRSYL( 'T', 'N', 1, NV, N, AV, LDAV, A, LDA,
     $                         DWORK, LDW, SCALE, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 3
                     RETURN
                  END IF
               END IF
C
C              Construct CS = DV'*C + BV'*X/SCALE,
C                        DS = DV'*D.
C              Additional workspace needed: MAX( PV*N, PV*M ).
C
C              Construct C <- DV'*C + BV'*X/SCALE.
C
               CALL DGEMM( 'T', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
               CALL DGEMM( 'T', 'N', PV, N, NV, ONE / SCALE, BV, LDBV,
     $                     DWORK, LDW, ONE, C, LDC )
C
C              Construct D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of V*G.
C
C           Total workspace needed:  NV*N + MAX( PV*N, PV*M ).
C
C           Compute -BV*C.
C           Workspace needed: NV*N.
C
            CALL DGEMM( 'N', 'N', NV, N, P, -ONE, BV, LDBV, C, LDC,
     $                  ZERO, DWORK, LDW )
C
C           Compute X and SCALE satisfying
C
C           AV*X - X*A + SCALE*BV*C = 0.
C
            IF( N.GT.0 ) THEN
               CALL DTRSYL( 'N', 'N', -1, NV, N, AV, LDAV, A, LDA,
     $                      DWORK, LDW, SCALE, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 3
                  RETURN
               END IF
            END IF
C
C           Construct CS = DV*C + CV*X/SCALE,
C                     DS = DV*D.
C           Additional workspace needed: MAX( PV*N, PV*M ).
C
C           Construct C <- DV*C + CV*X/SCALE.
C
            CALL DGEMM( 'N', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                  ZERO, DWORK(KW), PV )
            CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
            CALL DGEMM( 'N', 'N', PV, N, NV, ONE / SCALE, CV, LDCV,
     $                  DWORK, LDW, ONE, C, LDC )
C
C           Construct D <- DV*D.
C
            CALL DGEMM( 'N', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                  ZERO, DWORK(KW), PV )
            CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
         END IF
      ELSE
C
C        EV is a general matrix.
C
         IF( NV.GT.0 ) THEN
            TOLINF = TOLINF * DLANGE( '1', NV, NV, EV, LDEV, DWORK )
C
C           Reduce (AV,EV), or (AV',EV') or (EV',AV') to a generalized
C           real Schur form using an orthogonal equivalence
C           transformation and apply the orthogonal transformation
C           appropriately to BV and CV, or CV' and BV'.
C
C           Workspace needed:  2*NV*NV + MAX( 11*NV+16, NV*P, NV*PV );
C                              prefer larger.
C
            KQ  = 1
            KZ  = KQ  + NV*NV
            KAR = KZ  + NV*NV
            KAI = KAR + NV
            KB  = KAI + NV
            KW  = KB  + NV
C
            IF( CONJS ) THEN
               STDOM = 'S'
               ALPHA = ALPHA + SQRT( TOLINF )
C
C              Transpose AV and EV, if non-scalar.
C
               DO 10 I = 1, NV - 1
                  CALL DSWAP( NV-I, AV(I+1,I), 1, AV(I,I+1), LDAV )
                  CALL DSWAP( NV-I, EV(I+1,I), 1, EV(I,I+1), LDEV )
   10          CONTINUE
C
               IF( DISCR ) THEN
C
C                 Reduce (EV',AV') to a generalized real Schur form
C                 using orthogonal transformation matrices Q and Z
C                 such that Q'*EV'*Z results in a quasi-triangular form
C                 and Q'*AV'*Z results upper triangular.
C                 Total workspace needed: 2*NV*NV + 11*NV + 16.
C
                  EVTYPE = 'R'
                  CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                        DELCTG, NV, EV, LDEV, AV, LDAV, SDIM,
     $                        DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                        DWORK(KQ), LDW, DWORK(KZ), LDW,
     $                        DWORK(KW), LDWORK-KW+1, BWORK, IERR )
               ELSE
C
C                 Reduce (AV',EV') to a generalized real Schur form
C                 using orthogonal transformation matrices Q and Z
C                 such that Q'*AV'*Z results in a quasi-triangular form
C                 and Q'*EV'*Z results upper triangular.
C                 Total workspace needed: 2*NV*NV + 11*NV + 16.
C
                  EVTYPE = 'G'
                  CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                        DELCTG, NV, AV, LDAV, EV, LDEV, SDIM,
     $                        DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                        DWORK(KQ), LDW, DWORK(KZ), LDW,
     $                        DWORK(KW), LDWORK-KW+1, BWORK, IERR )
               END IF
               IF( IERR.NE.0 ) THEN
                  INFO = 1
                  RETURN
               END IF
               IF( STABCK ) THEN
C
C                 Check stability/antistability of generalized
C                 eigenvalues of the pair (AV,EV).
C
                  CALL AB09JX( DICO, STDOM, EVTYPE, NV, ALPHA,
     $                         DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                         TOLINF, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 4
                     RETURN
                  END IF
               END IF
               WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C              Compute Z'*BV and CV*Q.
C              Total workspace needed: 2*NV*NV + NV*MAX(P,PV).
C
               KW = KAR
               CALL DLACPY( 'Full', NV, PV, BV, LDBV, DWORK(KW), LDW )
               CALL DGEMM( 'T', 'N', NV, PV, NV, ONE, DWORK(KZ), LDW,
     $                     DWORK(KW), LDW, ZERO, BV, LDBV )
               CALL DLACPY( 'Full', P, NV, CV, LDCV, DWORK(KW), P )
               CALL DGEMM( 'N', 'N', P, NV, NV, ONE, DWORK(KW), P,
     $                     DWORK(KQ), LDW, ZERO, CV, LDCV )
            ELSE
C
C              Reduce (AV,EV) to a generalized real Schur form
C              using orthogonal transformation matrices Q and Z
C              such that Q'*AV*Z results in a quasi-triangular form
C              and Q'*EV*Z results upper triangular.
C              Total workspace needed: 2*NV*NV + 11*NV + 16.
C
               STDOM  = 'U'
               EVTYPE = 'G'
               ALPHA  = ALPHA - SQRT( TOLINF )
               CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                     DELCTG, NV, AV, LDAV, EV, LDEV, SDIM,
     $                     DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                     DWORK(KQ), LDW, DWORK(KZ), LDW,
     $                     DWORK(KW), LDWORK-KW+1, BWORK, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 1
                  RETURN
               END IF
               IF( STABCK ) THEN
C
C                 Check stability/antistability of generalized
C                 eigenvalues of the pair (AV,EV).
C
                  CALL AB09JX( DICO, STDOM, EVTYPE, NV, ALPHA,
     $                         DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                         TOLINF, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 4
                     RETURN
                  END IF
               END IF
               WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C              Compute Q'*BV and CV*Z.
C              Total workspace needed: 2*NV*NV + NV*MAX(P,PV).
C
               KW = KAR
               CALL DLACPY( 'Full', NV, P, BV, LDBV, DWORK(KW), LDW )
               CALL DGEMM( 'T', 'N', NV, P, NV, ONE, DWORK(KQ), LDW,
     $                     DWORK(KW), LDW, ZERO, BV, LDBV )
               CALL DLACPY( 'Full', PV, NV, CV, LDCV, DWORK(KW), PV )
               CALL DGEMM( 'N', 'N', PV, NV, NV, ONE, DWORK(KW), PV,
     $                     DWORK(KZ), LDW, ZERO, CV, LDCV )
            END IF
            WORK = MAX( WORK, DBLE( 2*NV*NV + NV*MAX( P, PV ) ) )
C
         END IF
C
         KC = 1
         KF = KC + NV*N
         KE = KF + NV*N
         KW = KE + N*N
         CALL DLASET( 'Full', NV, N, ZERO, ZERO, DWORK(KF), LDW )
C
         IF( CONJS ) THEN
C
C           Compute the projection of conj(V)*G.
C
C           Total workspace needed: NV*N + MAX( NV*N+N*N, PV*N, PV*M )
C
C           Compute CV'*C.
C           Workspace needed: NV*N.
C
            CALL DGEMM( 'T', 'N', NV, N, P, ONE, CV, LDCV, C, LDC,
     $                  ZERO, DWORK(KC), LDW )
C
            IF( DISCR ) THEN
C
C              Compute X and SCALE satisfying
C
C              EV'*X - AV'*X*A = SCALE*CV'*C by solving equivalently
C
C              EV'*X - Y*A = SCALE*CV'*C,
C              AV'*X - Y   = 0.
C
C              Additional workspace needed:
C              real    NV*N + N*N;
C              integer NV+N+6.
C
               IF( N.GT.0 ) THEN
                  CALL DLASET( 'Full', N, N, ZERO, ONE, DWORK(KE), LDWN
     $                       )
                  CALL DTGSYL( 'N', 0, NV, N, EV, LDEV, A, LDA,
     $                         DWORK(KC), LDW, AV, LDAV, DWORK(KE),
     $                         LDWN, DWORK(KF), LDW, SCALE, DIF,
     $                         DWORK(KW), LDWORK-KW+1, IWORK, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 2
                     RETURN
                  END IF
               END IF
C
C              Construct C <- DV'*C + BV'*X*A/SCALE,
C                        D <- DV'*D + BV'*X*B/SCALE.
C
C              Additional workspace needed: MAX( PV*N, PV*M ).
C
C              C <- DV'*C.
C
               KW = KF
               CALL DGEMM( 'T', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
C
C              D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
C
C              C <- C + BV'*X*A/SCALE.
C
               CALL DGEMM( 'T', 'N', PV, N, NV, ONE / SCALE, BV, LDBV,
     $                     DWORK(KC), LDW, ZERO, DWORK(KW), PV )
               CALL DGEMM( 'N', 'N', PV, N, N, ONE, DWORK(KW), PV,
     $                     A, LDA, ONE, C, LDC )
C
C              D <- D + BV'*X*B/SCALE.
C
               CALL DGEMM( 'N', 'N', PV, M, N, ONE, DWORK(KW), PV,
     $                     B, LDB, ONE, D, LDD )
            ELSE
C
C              Compute X and SCALE satisfying
C
C              AV'*X + EV'*X*A + SCALE*CV'*C = 0 by solving equivalently
C
C              AV'*X - Y*A    = -SCALE*CV'*C,
C              EV'*X - Y*(-I) = 0.
C
C              Additional workspace needed:
C              real    NV*N+N*N;
C              integer NV+N+6.
C
               IF( N.GT.0 ) THEN
                  CALL DLASET( 'Full', N, N, ZERO, -ONE, DWORK(KE), LDWN
     $                       )
                  CALL DTGSYL( 'N', 0, NV, N, AV, LDAV, A, LDA,
     $                         DWORK(KC), LDW, EV, LDEV, DWORK(KE),
     $                         LDWN, DWORK(KF), LDW, SCALE, DIF,
     $                         DWORK(KW), LDWORK-KW+1, IWORK, IERR )
C
C                 Note that the computed solution in DWORK(KC) is -X.
C
                  IF( IERR.NE.0 ) THEN
                     INFO = 2
                     RETURN
                  END IF
               END IF
C
C              Construct C <- DV'*C + BV'*X/SCALE.
C
               KW = KF
               CALL DGEMM( 'T', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
               CALL DGEMM( 'T', 'N', PV, N, NV, -ONE / SCALE, BV, LDBV,
     $                     DWORK(KC), LDW, ONE, C, LDC )
C
C              Construct D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                     ZERO, DWORK(KW), PV )
               CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of V*G.
C
C           Total workspace needed: NV*N + MAX( NV*N+N*N, PV*N, PV*M )
C
C           Compute -BV*C.
C           Workspace needed: NV*N.
C
            CALL DGEMM( 'N', 'N', NV, N, P, -ONE, BV, LDBV, C, LDC,
     $                  ZERO, DWORK, LDW )
C
C           Compute X and SCALE satisfying
C
C           AV*X - EV*X*A + SCALE*BV*C = 0 by solving equivalently
C
C           AV*X - Y*A = -SCALE*BV*C,
C           EV*X - Y   = 0.
C
C           Additional workspace needed:
C           real    NV*N + N*N;
C           integer NV+N+6.
C
            IF( N.GT.0 ) THEN
               CALL DLASET( 'Full', N, N, ZERO, ONE, DWORK(KE), LDWN )
               CALL DTGSYL( 'N', 0, NV, N, AV, LDAV, A, LDA,
     $                      DWORK(KC), LDW, EV, LDEV, DWORK(KE), LDWN,
     $                      DWORK(KF), LDW, SCALE, DIF, DWORK(KW),
     $                      LDWORK-KW+1, IWORK, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 2
                  RETURN
               END IF
            END IF
C
C           Construct C <- DV*C + CV*X/SCALE.
C
            KW = KF
            CALL DGEMM( 'N', 'N', PV, N, P, ONE, DV, LDDV, C, LDC,
     $                  ZERO, DWORK(KW), PV )
            CALL DLACPY( 'Full', PV, N, DWORK(KW), PV, C, LDC )
            CALL DGEMM( 'N', 'N', PV, N, NV, ONE / SCALE, CV, LDCV,
     $                  DWORK, LDW, ONE, C, LDC )
C
C           Construct D <- DV*D.
C
            CALL DGEMM( 'N', 'N', PV, M, P, ONE, DV, LDDV, D, LDD,
     $                  ZERO, DWORK(KW), PV )
            CALL DLACPY( 'Full', PV, M, DWORK(KW), PV, D, LDD )
         END IF
      END IF
C
      DWORK(1) = MAX( WORK, DBLE( LW ) )
C
      RETURN
C *** Last line of AB09JV ***
      END
