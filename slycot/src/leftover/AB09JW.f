      SUBROUTINE AB09JW( JOB, DICO, JOBEW, STBCHK, N, M, P, NW, MW,
     $                   A, LDA, B, LDB, C, LDC, D, LDD, AW, LDAW,
     $                   EW, LDEW, BW, LDBW, CW, LDCW, DW, LDDW, IWORK,
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
C     projection of G*W or G*conj(W) containing the poles of G, from the
C     state-space representations (A,B,C,D) and (AW-lambda*EW,BW,CW,DW),
C     of the transfer-function matrices G and W, respectively.
C     G is assumed to be a stable transfer-function matrix and
C     the state matrix A must be in a real Schur form.
C     When computing the stable projection of G*W, it is assumed
C     that G and W have completely distinct poles.
C     When computing the stable projection of G*conj(W), it is assumed
C     that G and conj(W) have completely distinct poles.
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
C             = 'W':  compute the projection of G*W containing
C                     the poles of G;
C             = 'C':  compute the projection of G*conj(W) containing
C                     the poles of G.
C
C     DICO    CHARACTER*1
C             Specifies the type of the systems as follows:
C             = 'C':  G and W are continuous-time systems;
C             = 'D':  G and W are discrete-time systems.
C
C     JOBEW   CHARACTER*1
C             Specifies whether EW is a general square or an identity
C             matrix as follows:
C             = 'G':  EW is a general square matrix;
C             = 'I':  EW is the identity matrix.
C
C     STBCHK  CHARACTER*1
C             Specifies whether stability/antistability of W is to be
C             checked as follows:
C             = 'C':  check stability if JOB = 'C' or antistability if
C                     JOB = 'W';
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
C             the transfer-function matrix G, and also the dimension
C             of the output vector if JOB = 'W', or of the input vector
C             if JOB = 'C', of the system with the transfer-function
C             matrix W.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of the output vector of the system with the
C             transfer-function matrix G.  P >= 0.
C
C     NW      (input) INTEGER
C             The dimension of the state vector of the system with the
C             transfer-function matrix W.  NW >= 0.
C
C     MW      (input) INTEGER
C             The dimension of the input vector, if JOB = 'W', or of
C             the output vector, if JOB = 'C', of the system with the
C             transfer-function matrix W.  MW >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the system with the transfer-function
C             matrix G in a real Schur form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array,
C             dimension (LDB,MAX(M,MW))
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading N-by-MW part of this
C             array contains the input matrix BS of the projection of
C             G*W, if JOB = 'W', or of G*conj(W), if JOB = 'C'.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain
C             the output/state matrix C of the system with the
C             transfer-function matrix G. The matrix CS is equal to C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array,
C             dimension (LDB,MAX(M,MW))
C             On entry, the leading P-by-M part of this array must
C             contain the feedthrough matrix D of the system with
C             the transfer-function matrix G.
C             On exit, if INFO = 0, the leading P-by-MW part of
C             this array contains the feedthrough matrix DS of the
C             projection of G*W, if JOB = 'W', or of G*conj(W),
C             if JOB = 'C'.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,P).
C
C     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW)
C             On entry, the leading NW-by-NW part of this array must
C             contain the state matrix AW of the system with the
C             transfer-function matrix W.
C             On exit, if INFO = 0, the leading NW-by-NW part of this
C             array contains a condensed matrix as follows:
C             if JOBEW = 'I', it contains the real Schur form of AW;
C             if JOBEW = 'G' and JOB = 'W', it contains a quasi-upper
C             triangular matrix representing the real Schur matrix
C             in the real generalized Schur form of the pair (AW,EW);
C             if JOBEW = 'G', JOB = 'C' and DICO = 'C', it contains a
C             quasi-upper triangular matrix corresponding to the
C             generalized real Schur form of the pair (AW',EW');
C             if JOBEW = 'G', JOB = 'C' and DICO = 'D', it contains an
C             upper triangular matrix corresponding to the generalized
C             real Schur form of the pair (EW',AW').
C
C     LDAW    INTEGER
C             The leading dimension of the array AW.  LDAW >= MAX(1,NW).
C
C     EW      (input/output) DOUBLE PRECISION array, dimension (LDEW,NW)
C             On entry, if JOBEW = 'G', the leading NW-by-NW part of
C             this array must contain the descriptor matrix EW of the
C             system with the transfer-function matrix W.
C             If JOBEW = 'I', EW is assumed to be an identity matrix
C             and is not referenced.
C             On exit, if INFO = 0 and JOBEW = 'G', the leading NW-by-NW
C             part of this array contains a condensed matrix as follows:
C             if JOB = 'W', it contains an upper triangular matrix
C             corresponding to the real generalized Schur form of the
C             pair (AW,EW);
C             if JOB = 'C' and DICO = 'C', it contains an upper
C             triangular matrix corresponding to the generalized real
C             Schur form of the pair (AW',EW');
C             if JOB = 'C' and DICO = 'D', it contains a quasi-upper
C             triangular matrix corresponding to the generalized
C             real Schur form of the pair (EW',AW').
C
C     LDEW    INTEGER
C             The leading dimension of the array EW.
C             LDEW >= MAX(1,NW), if JOBEW = 'G';
C             LDEW >= 1,         if JOBEW = 'I'.
C
C     BW      (input/output) DOUBLE PRECISION array,
C             dimension (LDBW,MBW), where MBW = MW, if JOB = 'W', and
C             MBW = M, if JOB = 'C'.
C             On entry, the leading NW-by-MBW part of this array must
C             contain the input matrix BW of the system with the
C             transfer-function matrix W.
C             On exit, if INFO = 0, the leading NW-by-MBW part of this
C             array contains Q'*BW, where Q is the orthogonal matrix
C             that reduces AW to the real Schur form or the left
C             orthogonal matrix used to reduce the pair (AW,EW),
C             (AW',EW') or (EW',AW') to the generalized real Schur form.
C
C     LDBW    INTEGER
C             The leading dimension of the array BW.  LDBW >= MAX(1,NW).
C
C     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW)
C             On entry, the leading PCW-by-NW part of this array must
C             contain the output matrix CW of the system with the
C             transfer-function matrix W, where PCW = M if JOB = 'W' or
C             PCW = MW if JOB = 'C'.
C             On exit, if INFO = 0, the leading PCW-by-NW part of this
C             array contains CW*Q, where Q is the orthogonal matrix that
C             reduces AW to the real Schur form, or CW*Z, where Z is the
C             right orthogonal matrix used to reduce the pair (AW,EW),
C             (AW',EW') or (EW',AW') to the generalized real Schur form.
C
C     LDCW    INTEGER
C             The leading dimension of the array CW.
C             LDCW >= MAX(1,PCW), where PCW = M if JOB = 'W', or
C             PCW = MW if JOB = 'C'.
C
C     DW      (input) DOUBLE PRECISION array,
C             dimension (LDDW,MBW), where MBW = MW if JOB = 'W', and
C             MBW = M if JOB = 'C'.
C             The leading PCW-by-MBW part of this array must contain
C             the feedthrough matrix DW of the system with the
C             transfer-function matrix W, where PCW = M if JOB = 'W',
C             or PCW = MW if JOB = 'C'.
C
C     LDDW    INTEGER
C             LDDW >= MAX(1,PCW), where PCW = M if JOB = 'W', or
C             PCW = MW if JOB = 'C'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK =   0,    if JOBEW = 'I';
C             LIWORK = NW+N+6, if JOBEW = 'G'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= LW1, if JOBEW = 'I',
C             LDWORK >= LW2, if JOBEW = 'G', where
C               LW1 = MAX( 1, NW*(NW+5), NW*N + MAX( a, N*MW, P*MW ) )
C                     a = 0,    if DICO = 'C' or  JOB = 'W',
C                     a = 2*NW, if DICO = 'D' and JOB = 'C';
C               LW2 = MAX( 2*NW*NW + MAX( 11*NW+16, NW*M, MW*NW ),
C                          NW*N + MAX( NW*N+N*N, MW*N, P*MW ) ).
C             For good performance, LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  the reduction of the pair (AW,EW) to the real
C                    generalized Schur form failed (JOBEW = 'G'),
C                    or the reduction of the matrix AW to the real
C                    Schur form failed (JOBEW = 'I);
C             =  2:  the solution of the Sylvester equation failed
C                    because the matrix A and the pencil AW-lambda*EW
C                    have common eigenvalues (if JOB = 'W'), or the
C                    pencil -AW-lambda*EW and A have common eigenvalues
C                    (if JOB = 'C' and DICO = 'C'), or the pencil
C                    AW-lambda*EW has an eigenvalue which is the
C                    reciprocal of one of eigenvalues of A
C                    (if JOB = 'C' and DICO = 'D');
C             =  3:  the solution of the Sylvester equation failed
C                    because the matrices A and AW have common
C                    eigenvalues (if JOB = 'W'), or the matrices A
C                    and -AW have common eigenvalues (if JOB = 'C' and
C                    DICO = 'C'), or the matrix A has an eigenvalue
C                    which is the reciprocal of one of eigenvalues of AW
C                    (if JOB = 'C' and DICO = 'D');
C             =  4:  JOB = 'W' and the pair (AW,EW) has not completely
C                    unstable generalized eigenvalues, or JOB = 'C' and
C                    the pair (AW,EW) has not completely stable
C                    generalized eigenvalues.
C
C     METHOD
C
C     If JOB = 'W', the matrices of the stable projection of G*W are
C     computed as
C
C       BS = B*DW + Y*BW,  CS = C,  DS = D*DW,
C
C     where Y satisfies the generalized Sylvester equation
C
C       -A*Y*EW + Y*AW + B*CW = 0.
C
C     If JOB = 'C', the matrices of the stable projection of G*conj(W)
C     are computed using the following formulas:
C
C     - for a continuous-time system, the matrices BS, CS and DS of
C       the stable projection are computed as
C
C         BS = B*DW' + Y*CW',  CS = C,  DS = D*DW',
C
C       where Y satisfies the generalized Sylvester equation
C
C         A*Y*EW' + Y*AW' + B*BW' = 0.
C
C     - for a discrete-time system, the matrices BS, CS and DS of
C       the stable projection are computed as
C
C         BS = B*DW' + A*Y*CW',  CS = C,  DS = D*DW' + C*Y*CW',
C
C       where Y satisfies the generalized Sylvester equation
C
C         Y*EW' - A*Y*AW' = B*BW'.
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
      CHARACTER         DICO, JOB, JOBEW, STBCHK
      INTEGER           INFO, LDA, LDAW, LDB, LDBW, LDC, LDCW,
     $                  LDD, LDDW, LDEW, LDWORK, M, MW, N, NW, P
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), AW(LDAW,*), B(LDB,*), BW(LDBW,*),
     $                  C(LDC,*), CW(LDCW,*), D(LDD,*), DW(LDDW,*),
     $                  DWORK(*), EW(LDEW,*)
C     .. Local Scalars ..
      CHARACTER*1       EVTYPE, STDOM
      LOGICAL           CONJS, DISCR, STABCK, UNITEW
      DOUBLE PRECISION  ALPHA, DIF, SCALE, TOLINF, WORK
      INTEGER           I, IA, IERR, KAI, KAR, KB, KC, KE, KF, KQ, KW,
     $                  KZ, LDW, LDWM, LDWN, LDWP, LW, SDIM
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
      UNITEW = LSAME( JOBEW,  'I' )
      STABCK = LSAME( STBCHK, 'C' )
C
      INFO = 0
      IF( UNITEW ) THEN
         IF ( DISCR .AND. CONJS ) THEN
            IA = 2*NW
         ELSE
            IA = 0
         END IF
         LW = MAX( 1, NW*( NW + 5 ), NW*N + MAX( IA, N*MW, P*MW ) )
      ELSE
         LW = MAX( 2*NW*NW + MAX( 11*NW+16, NW*M, MW*NW ),
     $             NW*N + MAX( NW*N + N*N, MW*N, P*MW ) )
      END IF
C
C     Test the input scalar arguments.
C
      LDW  = MAX( 1, NW )
      LDWM = MAX( 1, MW )
      LDWN = MAX( 1, N )
      LDWP = MAX( 1, P )
      IF( .NOT. ( LSAME( JOB, 'W' ) .OR. CONJS ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( DICO,   'C' ) .OR. DISCR  ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( JOBEW,  'G' ) .OR. UNITEW ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( LSAME( STBCHK, 'N' ) .OR. STABCK ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( NW.LT.0 ) THEN
         INFO = -8
      ELSE IF( MW.LT.0 ) THEN
         INFO = -9
      ELSE IF( LDA.LT.LDWN ) THEN
         INFO = -11
      ELSE IF( LDB.LT.LDWN ) THEN
         INFO = -13
      ELSE IF( LDC.LT.LDWP ) THEN
         INFO = -15
      ELSE IF( LDD.LT.LDWP ) THEN
         INFO = -17
      ELSE IF( LDAW.LT.LDW ) THEN
         INFO = -19
      ELSE IF( LDEW.LT.1 .OR. ( .NOT.UNITEW .AND. LDEW.LT.NW ) ) THEN
         INFO = -21
      ELSE IF( LDBW.LT.LDW ) THEN
         INFO = -23
      ELSE IF( ( .NOT.CONJS .AND. LDCW.LT.MAX( 1, M  ) ) .OR.
     $         (      CONJS .AND. LDCW.LT.LDWM ) ) THEN
         INFO = -25
      ELSE IF( ( .NOT.CONJS .AND. LDDW.LT.MAX( 1, M  ) ) .OR.
     $         (      CONJS .AND. LDDW.LT.LDWM ) ) THEN
         INFO = -27
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -30
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09JW', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( M.EQ.0 ) THEN
         CALL DLASET( 'Full', N, MW, ZERO, ZERO, B, LDB )
         CALL DLASET( 'Full', P, MW, ZERO, ZERO, D, LDD )
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
      WORK   = ONE
      TOLINF = DLAMCH( 'Epsilon' )
C
      IF( UNITEW ) THEN
C
C        EW is the identity matrix.
C
         IF( NW.GT.0 ) THEN
C
C           Reduce AW to the real Schur form using an orthogonal
C           similarity transformation AW <- Q'*AW*Q and apply the
C           transformation to BW and CW: BW <- Q'*BW and CW <- CW*Q.
C
C           Workspace needed:  NW*(NW+5);
C                              prefer larger.
C
            KW = NW*( NW + 2 ) + 1
            IF( CONJS ) THEN
               STDOM = 'S'
               ALPHA = ALPHA + SQRT( TOLINF )
               CALL TB01WD( NW, M, MW, AW, LDAW, BW, LDBW, CW, LDCW,
     $                      DWORK(2*NW+1), NW, DWORK, DWORK(NW+1),
     $                      DWORK(KW), LDWORK-KW+1, IERR )
            ELSE
               STDOM = 'U'
               ALPHA = ALPHA - SQRT( TOLINF )
               CALL TB01WD( NW, MW, M, AW, LDAW, BW, LDBW, CW, LDCW,
     $                      DWORK(2*NW+1), NW, DWORK, DWORK(NW+1),
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
               CALL AB09JX( DICO, STDOM, 'S', NW, ALPHA, DWORK,
     $                      DWORK(NW+1), DWORK, TOLINF, IERR )
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
        KW = NW*N + 1
         IF( CONJS ) THEN
C
C           Compute the projection of G*conj(W).
C
C           Total workspace needed:  NW*N + MAX( a, N*MW, P*MW ), where
C                                    a = 0,    if DICO = 'C',
C                                    a = 2*NW, if DICO = 'D'.
C
C           Compute -BW*B'.
C           Workspace needed: NW*N.
C
            CALL DGEMM( 'N', 'T', NW, N, M, -ONE, BW, LDBW, B, LDB,
     $                  ZERO, DWORK, LDW )
C
            IF( DISCR ) THEN
C
C              Compute Y' and SCALE satisfying
C
C              AW*Y'*A' - Y' = -SCALE*BW*B'.
C
C              Additional workspace needed: 2*NW.
C
               CALL SB04PY( 'N', 'T', -1, NW, N, AW, LDAW, A, LDA,
     $                      DWORK, LDW, SCALE, DWORK(KW), IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 3
                  RETURN
               END IF
C
C              Construct BS = B*DW' + A*Y*CW'/SCALE,
C                        DS = D*DW' + C*Y*CW'/SCALE.
C
C              Additional workspace needed: MAX( N*MW, P*MW ).
C
C              B <- B*DW'.
C
               CALL DGEMM( 'N', 'T', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWP )
               CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
C
C              B <- B + A*Y*CW'/SCALE.
C
               CALL DGEMM( 'T', 'T', N, MW, NW, ONE / SCALE, DWORK, LDW,
     $                     CW, LDCW, ZERO, DWORK(KW), LDWN )
               CALL DGEMM( 'N', 'N', N, MW, N, ONE, A, LDA,
     $                     DWORK(KW), LDWN, ONE, B, LDB )
C
C              D <- D + C*Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'N', P, MW, N, ONE, C, LDC,
     $                     DWORK(KW), LDWN, ONE, D, LDD )
            ELSE
C
C              Compute Y' and SCALE satisfying
C
C              AW*Y' + Y'*A' + SCALE*BW*B' = 0.
C
               IF( N.GT.0 ) THEN
                  CALL DTRSYL( 'N', 'T', 1, NW, N, AW, LDAW, A, LDA,
     $                         DWORK, LDW, SCALE, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 3
                     RETURN
                  END IF
               END IF
C
C              Construct BS = B*DW' + Y*CW'/SCALE,
C                        DS = D*DW'.
C
C              Additional workspace needed: MAX( N*MW, P*MW ).
C
C              Construct B <- B*DW' + Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'T', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
               CALL DGEMM( 'T', 'T', N, MW, NW, ONE / SCALE, DWORK, LDW,
     $                     CW, LDCW, ONE, B, LDB)
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWP )
               CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of G*W.
C
C           Total workspace needed:  NW*N + MAX( N*MW, P*MW ).
C
C           Compute B*CW.
C           Workspace needed: N*NW.
C
            CALL DGEMM( 'N', 'N', N, NW, M, ONE, B, LDB, CW, LDCW,
     $                  ZERO, DWORK, LDWN )
C
C           Compute Y and SCALE satisfying
C
C           A*Y - Y*AW - SCALE*B*CW = 0.
C
            IF( N.GT.0 ) THEN
               CALL DTRSYL( 'N', 'N', -1, N, NW, A, LDA, AW, LDAW,
     $                      DWORK, LDWN, SCALE, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 3
                  RETURN
               END IF
            END IF
C
C           Construct BS = B*DW + Y*BW/SCALE,
C                     DS = D*DW.
C
C           Additional workspace needed: MAX( N*MW, P*MW ).
C           Construct B <- B*DW + Y*BW/SCALE.
C
            CALL DGEMM( 'N', 'N', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                  ZERO, DWORK(KW), LDWN )
            CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
            CALL DGEMM( 'N', 'N', N, MW, NW, ONE / SCALE, DWORK, LDWN,
     $                  BW, LDBW, ONE, B, LDB)
C
C           D <- D*DW.
C
            CALL DGEMM( 'N', 'N', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                  ZERO, DWORK(KW), LDWP )
            CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
         END IF
      ELSE
C
C        EW is a general matrix.
C
         IF( NW.GT.0 ) THEN
            TOLINF = TOLINF * DLANGE( '1', NW, NW, EW, LDEW, DWORK )
C
C           Reduce (AW,EW), or (AW',EW') or (EW',AW') to a generalized
C           real Schur form using an orthogonal equivalence
C           transformation and apply the orthogonal transformation
C           appropriately to BW and CW, or CW' and BW'.
C
C           Workspace needed:  2*NW*NW + MAX( 11*NW+16, NW*M, MW*NW );
C                              prefer larger.
C
            KQ = 1
            KZ  = KQ  + NW*NW
            KAR = KZ  + NW*NW
            KAI = KAR + NW
            KB  = KAI + NW
            KW  = KB  + NW
C
            IF( CONJS ) THEN
               STDOM = 'S'
               ALPHA = ALPHA + SQRT( TOLINF )
C
C              Transpose AW and EW, if non-scalar.
C
               DO 10 I = 1, NW - 1
                  CALL DSWAP( NW-I, AW(I+1,I), 1, AW(I,I+1), LDAW )
                  CALL DSWAP( NW-I, EW(I+1,I), 1, EW(I,I+1), LDEW )
   10          CONTINUE
C
               IF( DISCR ) THEN
C
C                 Reduce (EW',AW') to a generalized real Schur form
C                 using orthogonal transformation matrices Q and Z
C                 such that Q'*EW'*Z results in a quasi-triangular form
C                 and Q'*AW'*Z results upper triangular.
C                 Total workspace needed: 2*NW*NW + 11*NW + 16.
C
                  EVTYPE = 'R'
                  CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                        DELCTG, NW, EW, LDEW, AW, LDAW, SDIM,
     $                        DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                        DWORK(KQ), LDW, DWORK(KZ), LDW,
     $                        DWORK(KW), LDWORK-KW+1, BWORK, IERR )
               ELSE
C
C                 Reduce (AW',EW') to a generalized real Schur form
C                 using orthogonal transformation matrices Q and Z
C                 such that Q'*AW'*Z results in a quasi-triangular form
C                 and Q'*EW'*Z results upper triangular.
C                 Total workspace needed: 2*NW*NW + 11*NW + 16.
C
                  EVTYPE = 'G'
                  CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                        DELCTG, NW, AW, LDAW, EW, LDEW, SDIM,
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
                  CALL AB09JX( DICO, STDOM, EVTYPE, NW, ALPHA,
     $                         DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                         TOLINF, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 4
                     RETURN
                  END IF
               END IF
               WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C              Compute Z'*BW and CW*Q.
C              Total workspace needed: 2*NW*NW + NW*MAX(M,MW).
C
               KW = KAR
               CALL DLACPY( 'Full', NW, M, BW, LDBW, DWORK(KW), LDW )
               CALL DGEMM( 'T', 'N', NW, M, NW, ONE, DWORK(KZ), LDW,
     $                     DWORK(KW), LDW, ZERO, BW, LDBW )
               CALL DLACPY( 'Full', MW, NW, CW, LDCW, DWORK(KW), LDWM )
               CALL DGEMM( 'N', 'N', MW, NW, NW, ONE, DWORK(KW), LDWM,
     $                     DWORK(KQ), LDW, ZERO, CW, LDCW )
            ELSE
C
C              Reduce (AW,EW) to a generalized real Schur form
C              using orthogonal transformation matrices Q and Z
C              such that Q'*AW*Z results in a quasi-triangular form
C              and Q'*EW*Z results upper triangular.
C              Total workspace needed: 2*NW*NW + 11*NW + 16.
C
               STDOM  = 'U'
               EVTYPE = 'G'
               ALPHA  = ALPHA - SQRT( TOLINF )
               CALL DGGES( 'Vectors', 'Vectors', 'Not ordered',
     $                     DELCTG, NW, AW, LDAW, EW, LDEW, SDIM,
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
                  CALL AB09JX( DICO, STDOM, EVTYPE, NW, ALPHA,
     $                         DWORK(KAR), DWORK(KAI), DWORK(KB),
     $                         TOLINF, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 4
                     RETURN
                  END IF
               END IF
               WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C              Compute Q'*BW and CW*Z.
C              Total workspace needed: 2*NW*NW + NW*MAX(M,MW).
C
               KW = KAR
               CALL DLACPY( 'Full', NW, MW, BW, LDBW, DWORK(KW), LDW )
               CALL DGEMM( 'T', 'N', NW, MW, NW, ONE, DWORK(KQ), LDW,
     $                     DWORK(KW), LDW, ZERO, BW, LDBW )
               CALL DLACPY( 'Full', M, NW, CW, LDCW, DWORK(KW), M )
               CALL DGEMM( 'N', 'N', M, NW, NW, ONE, DWORK(KW), M,
     $                     DWORK(KZ), LDW, ZERO, CW, LDCW )
            END IF
            WORK = MAX( WORK, DBLE( 2*NW*NW + NW*MAX( M, MW ) ) )
C
         END IF
C
         KC = 1
         KF = KC + NW*N
         KE = KF + NW*N
         KW = KE + N*N
         CALL DLASET( 'Full', N, NW, ZERO, ZERO, DWORK(KF), LDWN )
C
         IF( CONJS ) THEN
C
C           Compute the projection of G*conj(W).
C
C           Total workspace needed: NW*N + MAX( NW*N+N*N, MW*N, P*MW )
C
C           Compute B*BW'.
C           Workspace needed: N*NW.
C
            CALL DGEMM( 'N', 'T', N, NW, M, ONE, B, LDB, BW, LDBW,
     $                  ZERO, DWORK(KC), LDWN )
C
            IF( DISCR ) THEN
C
C              Compute Y and SCALE satisfying
C
C              Y*EW' - A*Y*AW' = SCALE*B*BW' by solving equivalently
C
C              A*X - Y*EW' = -SCALE*B*BW',
C              X   - Y*AW' = 0.
C
C              Additional workspace needed:
C              real    N*NW + N*N;
C              integer NW+N+6.
C
C
               IF( N.GT.0 ) THEN
                  CALL DLASET( 'Full', N, N, ZERO, ONE, DWORK(KE), LDWN
     $                       )
                  CALL DTGSYL( 'N', 0, N, NW, A, LDA, EW, LDEW,
     $                         DWORK(KC), LDWN, DWORK(KE), LDWN, AW,
     $                         LDAW, DWORK(KF), LDWN, SCALE, DIF,
     $                         DWORK(KW), LDWORK-KW+1, IWORK, IERR )
C
C                 Note that the computed solution in DWORK(KC) is -Y.
C
                  IF( IERR.NE.0 ) THEN
                     INFO = 2
                     RETURN
                  END IF
               END IF
C
C              Construct BS = B*DW' + A*Y*CW'/SCALE,
C                        DS = D*DW' + C*Y*CW'/SCALE.
C
C              Additional workspace needed: MAX( N*MW, P*MW ).
C
C              B <- B*DW'.
C
               CALL DGEMM( 'N', 'T', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWP )
               CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
C
C              B <- B + A*Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'T', N, MW, NW, -ONE / SCALE,
     $                     DWORK(KF), LDWN, CW, LDCW, ZERO,
     $                     DWORK(KW), LDWN )
               CALL DGEMM( 'N', 'N', N, MW, N, ONE, A, LDA,
     $                     DWORK(KW), LDWN, ONE, B, LDB )
C
C              D <- D + C*Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'N', P, MW, N, ONE, C, LDC,
     $                     DWORK(KW), LDWN, ONE, D, LDD )
            ELSE
C
C              Compute Y and SCALE satisfying
C
C              A*Y*EW' + Y*AW' + SCALE*B*BW' = 0 by solving equivalently
C
C              A*X    - Y*AW' = SCALE*B*BW',
C              (-I)*X - Y*EW' = 0.
C
C              Additional workspace needed:
C              real    N*NW+N*N;
C              integer NW+N+6.
C
               IF( N.GT.0 ) THEN
                  CALL DLASET( 'Full', N, N, ZERO, -ONE, DWORK(KE), LDWN
     $                       )
                  CALL DTGSYL( 'N', 0, N, NW, A, LDA, AW, LDAW,
     $                         DWORK(KC), LDWN, DWORK(KE), LDWN, EW,
     $                         LDEW, DWORK(KF), LDWN, SCALE, DIF,
     $                         DWORK(KW), LDWORK-KW+1, IWORK, IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = 2
                     RETURN
                  END IF
               END IF
C
C              Construct BS = B*DW' + Y*CW'/SCALE,
C                        DS = D*DW'.
C
C              Additional workspace needed: MAX( N*MW, P*MW ).
C
C              Construct B <- B*DW' + Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'T', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
               CALL DGEMM( 'N', 'T', N, MW, NW, ONE / SCALE,
     $                     DWORK(KF), LDWN, CW, LDCW, ONE, B, LDB )
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWP )
               CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of G*W.
C
C           Total workspace needed: NW*N + MAX( NW*N+N*N, MW*N, P*MW )
C
C           Compute B*CW.
C           Workspace needed: N*NW.
C
            CALL DGEMM( 'N', 'N', N, NW, M, ONE, B, LDB, CW, LDCW,
     $                  ZERO, DWORK(KC), LDWN )
C
C           Compute Y and SCALE satisfying
C
C           -A*Y*EW + Y*AW + B*CW = 0 by solving equivalently
C
C           A*X - Y*AW = SCALE*B*CW,
C           X   - Y*EW = 0.
C
C           Additional workspace needed:
C           real    N*NW + N*N;
C           integer NW+N+6.
C
            IF( N.GT.0 ) THEN
               CALL DLASET( 'Full', N, N, ZERO, ONE, DWORK(KE), LDWN )
               CALL DTGSYL( 'N', 0, N, NW, A, LDA, AW, LDAW,
     $                      DWORK(KC), LDWN, DWORK(KE), LDWN, EW, LDEW,
     $                      DWORK(KF), LDWN, SCALE, DIF, DWORK(KW),
     $                      LDWORK-KW+1, IWORK, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 2
                  RETURN
               END IF
            END IF
C
C           Construct BS = B*DW + Y*BW/SCALE,
C                     DS = D*DW.
C
C           Additional workspace needed: MAX( N*MW, P*MW ).
C           Construct B <- B*DW + Y*BW/SCALE.
C
            CALL DGEMM( 'N', 'N', N, MW, M, ONE, B, LDB, DW, LDDW,
     $                  ZERO, DWORK(KW), LDWN )
            CALL DLACPY( 'Full', N, MW, DWORK(KW), LDWN, B, LDB )
            CALL DGEMM( 'N', 'N', N, MW, NW, ONE / SCALE,
     $                  DWORK(KF), LDWN, BW, LDBW, ONE, B, LDB)
C
C           D <- D*DW.
C
            CALL DGEMM( 'N', 'N', P, MW, M, ONE, D, LDD, DW, LDDW,
     $                  ZERO, DWORK(KW), LDWP )
            CALL DLACPY( 'Full', P, MW, DWORK(KW), LDWP, D, LDD )
         END IF
      END IF
C
      DWORK(1) = MAX( WORK, DBLE( LW ) )
C
      RETURN
C *** Last line of AB09JW ***
      END
