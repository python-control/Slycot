      SUBROUTINE AB09KX( JOB, DICO, WEIGHT, N, NV, NW, M, P,
     $                   A,  LDA,  B,  LDB,  C,  LDC,  D,  LDD,
     $                   AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                   AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
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
C     To construct a state-space representation (A,BS,CS,DS) of the
C     stable projection of V*G*W or conj(V)*G*conj(W) from the
C     state-space representations (A,B,C,D), (AV,BV,CV,DV), and
C     (AW,BW,CW,DW) of the transfer-function matrices G, V and W,
C     respectively. G is assumed to be a stable transfer-function
C     matrix and the state matrix A must be in a real Schur form.
C     When computing the stable projection of V*G*W, V and W are assumed
C     to be completely unstable transfer-function matrices.
C     When computing the stable projection of conj(V)*G*conj(W),
C     V and W are assumed to be stable transfer-function matrices.
C
C     For a transfer-function matrix G, conj(G) denotes the conjugate
C     of G given by G'(-s) for a continuous-time system or G'(1/z)
C     for a discrete-time system.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Specifies which projection to be computed as follows:
C             = 'N':  compute the stable projection of V*G*W;
C             = 'C':  compute the stable projection of
C                     conj(V)*G*conj(W).
C
C     DICO    CHARACTER*1
C             Specifies the type of the systems as follows:
C             = 'C':  G, V and W are continuous-time systems;
C             = 'D':  G, V and W are discrete-time systems.
C
C     WEIGHT  CHARACTER*1
C             Specifies the type of frequency weighting, as follows:
C             = 'N':  no weightings are used (V = I, W = I);
C             = 'L':  only left weighting V is used (W = I);
C             = 'R':  only right weighting W is used (V = I);
C             = 'B':  both left and right weightings V and W are used.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. Also the number of rows of
C             the matrix B and the number of columns of the matrix C.
C             N represents the dimension of the state vector of the
C             system with the transfer-function matrix G.  N >= 0.
C
C     NV      (input) INTEGER
C             The order of the matrix AV. Also the number of rows of
C             the matrix BV and the number of columns of the matrix CV.
C             NV represents the dimension of the state vector of the
C             system with the transfer-function matrix V.  NV >= 0.
C
C     NW      (input) INTEGER
C             The order of the matrix AW. Also the number of rows of
C             the matrix BW and the number of columns of the matrix CW.
C             NW represents the dimension of the state vector of the
C             system with the transfer-function matrix W.  NW >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrices B, D, BW and DW
C             and number of rows of the matrices CW and DW.  M >= 0.
C             M represents the dimension of input vectors of the
C             systems with the transfer-function matrices G and W and
C             also the dimension of the output vector of the system
C             with the transfer-function matrix W.
C
C     P       (input) INTEGER
C             The number of rows of the matrices C, D, CV and DV and the
C             number of columns of the matrices BV and DV.  P >= 0.
C             P represents the dimension of output vectors of the
C             systems with the transfer-function matrices G and V and
C             also the dimension of the input vector of the system
C             with the transfer-function matrix V.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must
C             contain the state matrix A of the system with the
C             transfer-function matrix G in a real Schur form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading N-by-M part of this
C             array contains the input matrix BS of the stable
C             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W)
C             if JOB = 'C'.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading P-by-N part of this
C             array contains the output matrix CS of the stable
C             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W)
C             if JOB = 'C'.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the feedthrough matrix D of the system with the
C             transfer-function matrix G.
C             On exit, if INFO = 0, the leading P-by-M part of this
C             array contains the feedthrough matrix DS of the stable
C             projection of V*G*W if JOB = 'N', and of conj(V)*G*conj(W)
C             if JOB = 'C'.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,P).
C
C     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV)
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV
C             part of this array must contain the state matrix AV of
C             the system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             NV-by-NV part of this array contains a real Schur form
C             of AV.
C             AV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDAV    INTEGER
C             The leading dimension of the array AV.
C             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDAV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P)
C             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part
C             of this array must contain the input matrix BV of the
C             system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             NV-by-P part of this array contains the transformed input
C             matrix BV.
C             BV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDBV    INTEGER
C             The leading dimension of the array BV.
C             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDBV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV)
C             On entry, if WEIGHT = 'L' or 'B', the leading P-by-NV part
C             of this array must contain the output matrix CV of the
C             system with the transfer-function matrix V.
C             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading
C             P-by-NV part of this array contains the transformed output
C             matrix CV.
C             CV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDCV    INTEGER
C             The leading dimension of the array CV.
C             LDCV >= MAX(1,P), if WEIGHT = 'L' or 'B';
C             LDCV >= 1,        if WEIGHT = 'R' or 'N'.
C
C     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P)
C             If WEIGHT = 'L' or 'B', the leading P-by-P part of this
C             array must contain the feedthrough matrix DV of the system
C             with the transfer-function matrix V.
C             DV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDDV    INTEGER
C             The leading dimension of the array DV.
C             LDDV >= MAX(1,P), if WEIGHT = 'L' or 'B';
C             LDDV >= 1,        if WEIGHT = 'R' or 'N'.
C
C     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW
C             part of this array must contain the state matrix AW of
C             the system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             NW-by-NW part of this array contains a real Schur form
C             of AW.
C             AW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDAW    INTEGER
C             The leading dimension of the array AW.
C             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDAW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M)
C             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-M part
C             of this array must contain the input matrix BW of the
C             system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             NW-by-M part of this array contains the transformed input
C             matrix BW.
C             BW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDBW    INTEGER
C             The leading dimension of the array BW.
C             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDBW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW)
C             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part
C             of this array must contain the output matrix CW of the
C             system with the transfer-function matrix W.
C             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading
C             M-by-NW part of this array contains the transformed output
C             matrix CW.
C             CW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDCW    INTEGER
C             The leading dimension of the array CW.
C             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDCW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     DW      (input) DOUBLE PRECISION array, dimension (LDDW,M)
C             If WEIGHT = 'R' or 'B', the leading M-by-M part of this
C             array must contain the feedthrough matrix DW of the system
C             with the transfer-function matrix W.
C             DW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDDW    INTEGER
C             The leading dimension of the array DW.
C             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDDW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, LDW1, LDW2 ), where
C               LDW1 = 0 if WEIGHT = 'R' or 'N' and
C               LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) )
C                      if WEIGHT = 'L' or WEIGHT = 'B',
C               LDW2 = 0 if WEIGHT = 'L' or 'N' and
C               LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) )
C                      if WEIGHT = 'R' or WEIGHT = 'B',
C               a = 0,    b = 0,     if DICO = 'C' or  JOB = 'N',
C               a = 2*NV, b = 2*NW,  if DICO = 'D' and JOB = 'C'.
C             For good performance, LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             =  0:  no warning;
C             =  1:  JOB = 'N' and AV is not completely unstable, or
C                    JOB = 'C' and AV is not stable;
C             =  2:  JOB = 'N' and AW is not completely unstable, or
C                    JOB = 'C' and AW is not stable;
C             =  3:  both above conditions appear.
C
C     Error Indicator
C
C     INFO    INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  the reduction of AV to a real Schur form failed;
C             =  2:  the reduction of AW to a real Schur form failed;
C             =  3:  the solution of the Sylvester equation failed
C                    because the matrices A and AV have common
C                    eigenvalues (if JOB = 'N'), or -AV and A have
C                    common eigenvalues (if JOB = 'C' and DICO = 'C'),
C                    or AV has an eigenvalue which is the reciprocal of
C                    one of the eigenvalues of A (if JOB = 'C' and
C                    DICO = 'D');
C             =  4:  the solution of the Sylvester equation failed
C                    because the matrices A and AW have common
C                    eigenvalues (if JOB = 'N'), or -AW and A have
C                    common eigenvalues (if JOB = 'C' and DICO = 'C'),
C                    or AW has an eigenvalue which is the reciprocal of
C                    one of the eigenvalues of A (if JOB = 'C' and
C                    DICO = 'D').
C
C     METHOD
C
C     The matrices of the stable projection of V*G*W are computed as
C
C       BS = B*DW + Y*BW,  CS = CV*X + DV*C,  DS = DV*D*DW,
C
C     where X and Y satisfy the continuous-time Sylvester equations
C
C       AV*X - X*A  + BV*C = 0,
C       -A*Y + Y*AW + B*CW = 0.
C
C     The matrices of the stable projection of conj(V)*G*conj(W) are
C     computed using the explicit formulas established in [1].
C
C     For a continuous-time system, the matrices BS, CS and DS of
C     the stable projection are computed as
C
C       BS = B*DW' + Y*CW',  CS = BV'*X + DV'*C,  DS = DV'*D*DW',
C
C     where X and Y satisfy the continuous-time Sylvester equations
C
C       AV'*X + X*A   + CV'*C = 0,
C         A*Y + Y*AW' + B*BW' = 0.
C
C     For a discrete-time system, the matrices BS, CS and DS of
C     the stable projection are computed as
C
C       BS = B*DW' + A*Y*CW',  CS = BV'*X*A + DV'*C,
C       DS = DV'*D*DW' + BV'*X*B*DW' + DV'*C*Y*CW' + BV'*X*A*Y*CW',
C
C     where X and Y satisfy the discrete-time Sylvester equations
C
C       AV'*X*A + CV'*C = X,
C       A*Y*AW' + B*BW' = Y.
C
C     REFERENCES
C
C     [1] Varga A.
C         Explicit formulas for an efficient implementation
C         of the frequency-weighting model reduction approach.
C         Proc. 1993 European Control Conference, Groningen, NL,
C         pp. 693-696, 1993.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on numerically stable algorithms.
C
C     FURTHER COMMENTS
C
C     The matrix A must be stable, but its stability is not checked by
C     this routine.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, April 2000.
C     D. Sima, University of Bucharest, May 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2000.
C     Based on the RASP routines SFRLW, SFRLW1, SFRRW and SFRRW1,
C     by A. Varga, 1992.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Frequency weighting, model reduction, multivariable system,
C     state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO, JOB, WEIGHT
      INTEGER          INFO, IWARN, LDA, LDAV, LDAW, LDB, LDBV, LDBW,
     $                 LDC, LDCV, LDCW, LDD, LDDV, LDDW, LDWORK, M, N,
     $                 NV, NW, P
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),   B(LDB,*),   C(LDC,*),   D(LDD,*),
     $                 AV(LDAV,*), BV(LDBV,*), CV(LDCV,*), DV(LDDV,*),
     $                 AW(LDAW,*), BW(LDBW,*), CW(LDCW,*), DW(LDDW,*),
     $                 DWORK(*)
C     .. Local Scalars
      LOGICAL          CONJS, DISCR, FRWGHT, LEFTW, RIGHTW
      DOUBLE PRECISION SCALE, WORK
      INTEGER          I, IA, IB, IERR, KW, LDW, LDWN, LW
C     .. External Functions ..
      LOGICAL          LSAME
      DOUBLE PRECISION DLAPY2
      EXTERNAL         DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL         DGEMM, DLACPY, DTRSYL, SB04PY, TB01WD, XERBLA
C     .. Executable Statements ..
C
      CONJS  = LSAME( JOB,    'C' )
      DISCR  = LSAME( DICO,   'D' )
      LEFTW  = LSAME( WEIGHT, 'L' ) .OR. LSAME( WEIGHT, 'B' )
      RIGHTW = LSAME( WEIGHT, 'R' ) .OR. LSAME( WEIGHT, 'B' )
      FRWGHT = LEFTW .OR. RIGHTW
C
      IWARN = 0
      INFO  = 0
      IF ( DISCR .AND. CONJS ) THEN
         IA = 2*NV
         IB = 2*NW
      ELSE
         IA = 0
         IB = 0
      END IF
      LW = 1
      IF( LEFTW )
     $   LW = MAX( LW, NV*( NV + 5 ), NV*N + MAX( IA, P*N, P*M ) )
      IF( RIGHTW )
     $   LW = MAX( LW, NW*( NW + 5 ), NW*N + MAX( IB, M*N, P*M ) )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( JOB, 'N' ) .OR. CONJS ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NV.LT.0 ) THEN
         INFO = -5
      ELSE IF( NW.LT.0 ) THEN
         INFO = -6
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( P.LT.0 ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDAV.LT.1 .OR. ( LEFTW  .AND. LDAV.LT.NV ) ) THEN
         INFO = -18
      ELSE IF( LDBV.LT.1 .OR. ( LEFTW  .AND. LDBV.LT.NV ) ) THEN
         INFO = -20
      ELSE IF( LDCV.LT.1 .OR. ( LEFTW  .AND. LDCV.LT.P  ) ) THEN
         INFO = -22
      ELSE IF( LDDV.LT.1 .OR. ( LEFTW  .AND. LDDV.LT.P  ) ) THEN
         INFO = -24
      ELSE IF( LDAW.LT.1 .OR. ( RIGHTW .AND. LDAW.LT.NW ) ) THEN
         INFO = -26
      ELSE IF( LDBW.LT.1 .OR. ( RIGHTW .AND. LDBW.LT.NW ) ) THEN
         INFO = -28
      ELSE IF( LDCW.LT.1 .OR. ( RIGHTW .AND. LDCW.LT.M  ) ) THEN
         INFO = -30
      ELSE IF( LDDW.LT.1 .OR. ( RIGHTW .AND. LDDW.LT.M  ) ) THEN
         INFO = -32
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -34
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09KX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( .NOT.FRWGHT .OR. MIN( M, P ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      WORK = ONE
      IF( LEFTW .AND. NV.GT.0 ) THEN
C
C        Reduce AV to a real Schur form using an orthogonal similarity
C        transformation AV <- Q'*AV*Q and apply the transformation to
C        BV and CV: BV <- Q'*BV and CV <- CV*Q.
C
C        Workspace needed:  NV*(NV+5);
C                           prefer larger.
C
         KW = NV*( NV + 2 ) + 1
         CALL TB01WD( NV, P, P, AV, LDAV, BV, LDBV, CV, LDCV,
     $                DWORK(2*NV+1), NV, DWORK, DWORK(NV+1), DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
         IF( CONJS ) THEN
C
C           Check the stability of the eigenvalues of AV.
C
            IF ( DISCR ) THEN
               DO 10 I = 1, NV
                  IF( DLAPY2( DWORK(I), DWORK(NV+I) ).GE.ONE) THEN
                     IWARN = 1
                     GO TO 50
                  END IF
   10          CONTINUE
            ELSE
               DO 20 I = 1, NV
                  IF( DWORK(I).GE.ZERO ) THEN
                     IWARN = 1
                     GO TO 50
                  END IF
   20          CONTINUE
            END IF
         ELSE
C
C           Check the anti-stability of the eigenvalues of AV.
C
            IF ( DISCR ) THEN
               DO 30 I = 1, NV
                  IF( DLAPY2( DWORK(I), DWORK(NV+I) ).LE.ONE) THEN
                     IWARN = 1
                     GO TO 50
                  END IF
   30          CONTINUE
            ELSE
               DO 40 I = 1, NV
                  IF( DWORK(I).LE.ZERO ) THEN
                     IWARN = 1
                     GO TO 50
                  END IF
   40          CONTINUE
            END IF
         END IF
   50    CONTINUE
C
      END IF
C
      IF( RIGHTW .AND. NW.GT.0 ) THEN
C
C        Reduce AW to a real Schur form using an orthogonal similarity
C        transformation AW <- T'*AW*T and apply the transformation to
C        BW and CW: BW <- T'*BW and CW <- CW*T.
C
C        Workspace needed:  NW*(NW+5);
C                           prefer larger.
C
         KW = NW*( NW + 2 ) + 1
         CALL TB01WD( NW, M, M, AW, LDAW, BW, LDBW, CW, LDCW,
     $                DWORK(2*NW+1), NW, DWORK, DWORK(NW+1), DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
         IF( CONJS ) THEN
C
C           Check the stability of the eigenvalues of AW.
C
            IF ( DISCR ) THEN
               DO 60 I = 1, NW
                  IF( DLAPY2( DWORK(I), DWORK(NW+I) ).GE.ONE) THEN
                     IWARN = IWARN + 2
                     GO TO 100
                  END IF
   60          CONTINUE
            ELSE
               DO 70 I = 1, NW
                  IF( DWORK(I).GE.ZERO ) THEN
                     IWARN = IWARN + 2
                     GO TO 100
                  END IF
   70          CONTINUE
            END IF
         ELSE
C
C           Check the anti-stability of the eigenvalues of AW.
C
            IF ( DISCR ) THEN
               DO 80 I = 1, NW
                  IF( DLAPY2( DWORK(I), DWORK(NW+I) ).LE.ONE) THEN
                     IWARN = IWARN + 2
                     GO TO 100
                  END IF
   80          CONTINUE
            ELSE
               DO 90 I = 1, NW
                  IF( DWORK(I).LE.ZERO ) THEN
                     IWARN = IWARN + 2
                     GO TO 100
                  END IF
   90          CONTINUE
            END IF
         END IF
  100    CONTINUE
      END IF
C
      IF( LEFTW ) THEN
         LDW = MAX( NV, 1 )
         KW  = NV*N + 1
         IF( CONJS ) THEN
C
C           Compute the projection of conj(V)*G.
C
C           Total workspace needed:  NV*N + MAX( a, P*N, P*M ), where
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
C              Construct C <- DV'*C + BV'*X*A/SCALE,
C                        D <- DV'*D + BV'*X*B/SCALE.
C
C              Additional workspace needed: MAX( P*N, P*M ).
C
C              C <- DV'*C.
C
               CALL DGEMM( 'T', 'N', P, N, P, ONE, DV, LDDV, C, LDC,
     $                      ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, N, DWORK(KW), P, C, LDC )
C
C              D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', P, M, P, ONE, DV, LDDV, D, LDD,
     $                      ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
C
C              C <- C + BV'*X*A/SCALE.
C
               CALL DGEMM( 'T', 'N', P, N, NV, ONE / SCALE, BV, LDBV,
     $                     DWORK, LDW, ZERO, DWORK(KW), P )
               CALL DGEMM( 'N', 'N', P, N, N, ONE, DWORK(KW), P, A, LDA,
     $                     ONE, C, LDC )
C
C              D <- D + BV'*X*B/SCALE.
C
               CALL DGEMM( 'N', 'N', P, M, N, ONE, DWORK(KW), P, B, LDB,
     $                     ONE, D, LDD )
            ELSE
C
C              Compute X and SCALE satisfying
C
C              AV'*X + X*A + SCALE*CV'*C = 0.
C
               CALL DTRSYL( 'T', 'N', 1, NV, N, AV, LDAV, A, LDA,
     $                      DWORK, LDW, SCALE, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 3
                  RETURN
               END IF
C
C              Construct C and D.
C              Additional workspace needed: MAX( P*N, P*M ).
C
C              Construct C <- BV'*X/SCALE + DV'*C.
C
               CALL DGEMM( 'T', 'N', P, N, P, ONE, DV, LDDV, C, LDC,
     $                     ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, N, DWORK(KW), P, C, LDC )
               CALL DGEMM( 'T', 'N', P, N, NV, ONE / SCALE, BV, LDBV,
     $                     DWORK, LDW, ONE, C, LDC )
C
C              Construct D <- DV'*D.
C
               CALL DGEMM( 'T', 'N', P, M, P, ONE, DV, LDDV, D, LDD,
     $                     ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of V*G.
C
C           Total workspace needed:  NV*N + MAX( P*N, P*M ).
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
            CALL DTRSYL( 'N', 'N', -1, NV, N, AV, LDAV, A, LDA,
     $                   DWORK, LDW, SCALE, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 3
               RETURN
            END IF
C
C           Construct C <- CV*X/SCALE + DV*C.
C
            CALL DGEMM( 'N', 'N', P, N, P, ONE, DV, LDDV, C, LDC,
     $                  ZERO, DWORK(KW), P )
            CALL DLACPY( 'F', P, N, DWORK(KW), P, C, LDC )
            CALL DGEMM( 'N', 'N', P, N, NV, ONE / SCALE, CV, LDCV,
     $                  DWORK, LDW, ONE, C, LDC )
C
C           Construct D <- DV*D.
C
            CALL DGEMM( 'N', 'N', P, M, P, ONE, DV, LDDV, D, LDD,
     $                  ZERO, DWORK(KW), P )
            CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
         END IF
      END IF
C
      IF( RIGHTW ) THEN
         LDWN = MAX( N, 1 )
         KW   = N*NW + 1
         IF( CONJS ) THEN
C
C           Compute the projection of G*conj(W) or of conj(V)*G*conj(W).
C
C           Total workspace needed:  NW*N + MAX( b, M*N, P*M ), where
C                                    b = 0,    if DICO = 'C',
C                                    b = 2*NW, if DICO = 'D'.
C
C           Compute -BW*B'.
C           Workspace needed: N*NW.
C
            LDW = MAX( NW, 1 )
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
                  INFO = 4
                  RETURN
               END IF
C
C              Construct B <- B*DW' + A*Y*CW'/SCALE,
C                        D <- D*DW' + C*Y*CW'/SCALE.
C
C              Additional workspace needed: MAX( N*M, P*M ).
C
C              B <- B*DW'.
C
               CALL DGEMM( 'N', 'T', N, M, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'F', N, M, DWORK(KW), LDWN, B, LDB )
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, M, M, ONE, D, LDD, DW, LDDW,
     $                    ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
C
C              B <- B + A*Y*CW'/SCALE.
C
               CALL DGEMM( 'T', 'T', N, M, NW, ONE / SCALE, DWORK, LDW,
     $                     CW, LDCW, ZERO, DWORK(KW), LDWN )
               CALL DGEMM( 'N', 'N', N, M, N, ONE, A, LDA,
     $                     DWORK(KW), LDWN, ONE, B, LDB )
C
C              D <- D + C*Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'N', P, M, N, ONE, C, LDC,
     $                     DWORK(KW), LDWN, ONE, D, LDD )
            ELSE
C
C              Compute Y' and SCALE satisfying
C
C              AW*Y' + Y'*A' + SCALE*BW*B' = 0.
C
               CALL DTRSYL( 'N', 'T', 1, NW, N, AW, LDAW, A, LDA,
     $                      DWORK, LDW, SCALE, IERR )
               IF( IERR.NE.0 ) THEN
                  INFO = 4
                  RETURN
               END IF
C
C              Construct B and D.
C              Additional workspace needed: MAX( N*M, P*M ).
C
C              Construct B <- B*DW' + Y*CW'/SCALE.
C
               CALL DGEMM( 'N', 'T', N, M, M, ONE, B, LDB, DW, LDDW,
     $                     ZERO, DWORK(KW), LDWN )
               CALL DLACPY( 'F', N, M, DWORK(KW), LDWN, B, LDB )
               CALL DGEMM( 'T', 'T', N, M, NW, ONE / SCALE, DWORK, LDW,
     $                     CW, LDCW, ONE, B, LDB)
C
C              D <- D*DW'.
C
               CALL DGEMM( 'N', 'T', P, M, M, ONE, D, LDD, DW, LDDW,
     $                     ZERO, DWORK(KW), P )
               CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
            END IF
         ELSE
C
C           Compute the projection of G*W or of V*G*W.
C
C           Total workspace needed:  NW*N + MAX( M*N, P*M ).
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
            CALL DTRSYL( 'N', 'N', -1, N, NW, A, LDA, AW, LDAW,
     $                   DWORK, LDWN, SCALE, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 4
               RETURN
            END IF
C
C           Construct B and D.
C           Additional workspace needed: MAX( N*M, P*M ).
C           Construct B <- B*DW + Y*BW/SCALE.
C
            CALL DGEMM( 'N', 'N', N, M, M, ONE, B, LDB, DW, LDDW,
     $                  ZERO, DWORK(KW), LDWN )
            CALL DLACPY( 'F', N, M, DWORK(KW), LDWN, B, LDB )
            CALL DGEMM( 'N', 'N', N, M, NW, ONE / SCALE, DWORK, LDWN,
     $                  BW, LDBW, ONE, B, LDB)
C
C           D <- D*DW.
C
            CALL DGEMM( 'N', 'N', P, M, M, ONE, D, LDD, DW, LDDW,
     $                  ZERO, DWORK(KW), P )
            CALL DLACPY( 'F', P, M, DWORK(KW), P, D, LDD )
         END IF
      END IF
C
      DWORK(1) = MAX( WORK, DBLE( LW ) )
C
      RETURN
C *** Last line of AB09KX ***
      END
