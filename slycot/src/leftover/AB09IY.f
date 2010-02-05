      SUBROUTINE AB09IY( DICO, JOBC, JOBO, WEIGHT, N, M, P, NV, PV,
     $                   NW, MW, ALPHAC, ALPHAO, A, LDA, B, LDB, C, LDC,
     $                   AV, LDAV, BV, LDBV, CV, LDCV, DV, LDDV,
     $                   AW, LDAW, BW, LDBW, CW, LDCW, DW, LDDW,
     $                   SCALEC, SCALEO, S, LDS, R, LDR,
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
C     To compute for given state-space representations
C     (A,B,C,0), (AV,BV,CV,DV), and (AW,BW,CW,DW) of the
C     transfer-function matrices G, V and W, respectively,
C     the Cholesky factors of the frequency-weighted
C     controllability and observability Grammians corresponding
C     to a frequency-weighted model reduction problem.
C     G, V and W must be stable transfer-function matrices with
C     the state matrices A, AV, and AW in real Schur form.
C     It is assumed that the state space realizations (AV,BV,CV,DV)
C     and (AW,BW,CW,DW) are minimal. In case of possible pole-zero
C     cancellations in forming V*G and/or G*W, the parameters for the
C     choice of frequency-weighted Grammians ALPHAO and/or ALPHAC,
C     respectively, must be different from 1.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the systems as follows:
C             = 'C':  G, V and W are continuous-time systems;
C             = 'D':  G, V and W are discrete-time systems.
C
C     JOBC    CHARACTER*1
C             Specifies the choice of frequency-weighted controllability
C             Grammian as follows:
C             = 'S': choice corresponding to a combination method [4]
C                    of the approaches of Enns [1] and Lin-Chiu [2,3];
C             = 'E': choice corresponding to the stability enhanced
C                    modified combination method of [4].
C
C     JOBO    CHARACTER*1
C             Specifies the choice of frequency-weighted observability
C             Grammian as follows:
C             = 'S': choice corresponding to a combination method [4]
C                    of the approaches of Enns [1] and Lin-Chiu [2,3];
C             = 'E': choice corresponding to the stability enhanced
C                    modified combination method of [4].
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
C             The order of the state-space representation of G, i.e.,
C             the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrix B and
C             the number of rows of the matrices CW and DW.  M >= 0.
C             M represents the dimension of the input vector of the
C             system with the transfer-function matrix G and
C             also the dimension of the output vector of the system
C             with the transfer-function matrix W.
C
C     P       (input) INTEGER
C             The number of rows of the matrix C and the
C             number of columns of the matrices BV and DV.  P >= 0.
C             P represents the dimension of the output vector of the
C             system with the transfer-function matrix G and
C             also the dimension of the input vector of the system
C             with the transfer-function matrix V.
C
C     NV      (input) INTEGER
C             The order of the matrix AV. Also the number of rows of
C             the matrix BV and the number of columns of the matrix CV.
C             NV represents the dimension of the state vector of the
C             system with the transfer-function matrix V.  NV >= 0.
C
C     PV      (input) INTEGER
C             The number of rows of the matrices CV and DV.  PV >= 0.
C             PV represents the dimension of the output vector of the
C             system with the transfer-function matrix V.
C
C     NW      (input) INTEGER
C             The order of the matrix AW. Also the number of rows of
C             the matrix BW and the number of columns of the matrix CW.
C             NW represents the dimension of the state vector of the
C             system with the transfer-function matrix W.  NW >= 0.
C
C     MW      (input) INTEGER
C             The number of columns of the matrices BW and DW.  MW >= 0.
C             MW represents the dimension of the input vector of the
C             system with the transfer-function matrix W.
C
C     ALPHAC  (input) DOUBLE PRECISION
C             Combination method parameter for defining the
C             frequency-weighted controllability Grammian (see METHOD);
C             ABS(ALPHAC) <= 1.
C
C     ALPHAO  (input) DOUBLE PRECISION
C             Combination method parameter for defining the
C             frequency-weighted observability Grammian (see METHOD);
C             ABS(ALPHAO) <= 1.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must
C             contain the state matrix A (of the system with the
C             transfer-function matrix G) in a real Schur form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     AV      (input) DOUBLE PRECISION array, dimension (LDAV,NV)
C             If WEIGHT = 'L' or 'B', the leading NV-by-NV part of this
C             array must contain the state matrix AV (of the system with
C             the transfer-function matrix V) in a real Schur form.
C             AV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDAV    INTEGER
C             The leading dimension of array AV.
C             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDAV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     BV      (input) DOUBLE PRECISION array, dimension (LDBV,P)
C             If WEIGHT = 'L' or 'B', the leading NV-by-P part of this
C             array must contain the input matrix BV of the system with
C             the transfer-function matrix V.
C             BV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDBV    INTEGER
C             The leading dimension of array BV.
C             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B';
C             LDBV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     CV      (input) DOUBLE PRECISION array, dimension (LDCV,NV)
C             If WEIGHT = 'L' or 'B', the leading PV-by-NV part of this
C             array must contain the output matrix CV of the system with
C             the transfer-function matrix V.
C             CV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDCV    INTEGER
C             The leading dimension of array CV.
C             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B';
C             LDCV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P)
C             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this
C             array must contain the feedthrough matrix DV of the system
C             with the transfer-function matrix V.
C             DV is not referenced if WEIGHT = 'R' or 'N'.
C
C     LDDV    INTEGER
C             The leading dimension of array DV.
C             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B';
C             LDDV >= 1,         if WEIGHT = 'R' or 'N'.
C
C     AW      (input) DOUBLE PRECISION array, dimension (LDAW,NW)
C             If WEIGHT = 'R' or 'B', the leading NW-by-NW part of this
C             array must contain the state matrix AW (of the system with
C             the transfer-function matrix W) in a real Schur form.
C             AW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDAW    INTEGER
C             The leading dimension of array AW.
C             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDAW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     BW      (input) DOUBLE PRECISION array, dimension (LDBW,MW)
C             If WEIGHT = 'R' or 'B', the leading NW-by-MW part of this
C             array must contain the input matrix BW of the system with
C             the transfer-function matrix W.
C             BW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDBW    INTEGER
C             The leading dimension of array BW.
C             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B';
C             LDBW >= 1,         if WEIGHT = 'L' or 'N'.
C
C     CW      (input) DOUBLE PRECISION array, dimension (LDCW,NW)
C             If WEIGHT = 'R' or 'B', the leading M-by-NW part of this
C             array must contain the output matrix CW of the system with
C             the transfer-function matrix W.
C             CW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDCW    INTEGER
C             The leading dimension of array CW.
C             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDCW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW)
C             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this
C             array must contain the feedthrough matrix DW of the system
C             with the transfer-function matrix W.
C             DW is not referenced if WEIGHT = 'L' or 'N'.
C
C     LDDW    INTEGER
C             The leading dimension of array DW.
C             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B';
C             LDDW >= 1,        if WEIGHT = 'L' or 'N'.
C
C     SCALEC  (output) DOUBLE PRECISION
C             Scaling factor for the controllability Grammian in (1)
C             or (3). See METHOD.
C
C     SCALEO  (output) DOUBLE PRECISION
C             Scaling factor for the observability Grammian in (2)
C             or (4). See METHOD.
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor S of the frequency-weighted
C             cotrollability Grammian P = S*S'. See METHOD.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= MAX(1,N).
C
C     R       (output) DOUBLE PRECISION array, dimension (LDR,N)
C             The leading N-by-N upper triangular part of this array
C             contains the Cholesky factor R of the frequency-weighted
C             observability Grammian Q = R'*R. See METHOD.
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, LLEFT, LRIGHT ),
C             where
C             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5)
C                              if WEIGHT = 'L' or 'B' and PV > 0;
C             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0;
C             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5)
C                              if WEIGHT = 'R' or 'B' and MW > 0;
C             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0.
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the state matrices A and/or AV are not stable or
C                   not in a real Schur form;
C             = 2:  if the state matrices A and/or AW are not stable or
C                   not in a real Schur form;
C             = 3:  eigenvalues computation failure.
C
C     METHOD
C
C     Let Pi = Si*Si' and Qo = Ro'*Ro be the Cholesky factored
C     controllability and observability Grammians satisfying
C     in the continuous-time case
C
C            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0,       (1)
C
C            Ao'*Qo + Qo*Ao +  scaleo^2*Co'*Co = 0,       (2)
C
C     and in the discrete-time case
C
C            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0,       (3)
C
C            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0,       (4)
C
C     where
C
C           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) ,
C                ( 0   Aw  )          (  Bw  )
C
C           Ao = (  A   0  ) ,   Co = ( Dv*C  Cv ) .
C                ( Bv*C Av )
C
C     Consider the partitioned Grammians
C
C           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) ,
C                ( P12' P22 )               ( Q12' Q22 )
C
C     where P11 and Q11 are the leading N-by-N parts of Pi and Qo,
C     respectively, and let P0 and Q0 be non-negative definite matrices
C     defined in the combination method [4]
C                                        -1
C            P0 = P11 - ALPHAC**2*P12*P22 *P21 ,
C                                        -1
C            Q0 = Q11 - ALPHAO**2*Q12*Q22 *Q21.
C
C     The frequency-weighted controllability and observability
C     Grammians, P and Q, respectively, are defined as follows:
C     P = P0 if JOBC = 'S' (standard combination method [4]);
C     P = P1 >= P0 if JOBC = 'E', where P1 is the controllability
C     Grammian defined to enforce stability for a modified combination
C     method of [4];
C     Q = Q0 if JOBO = 'S' (standard combination method [4]);
C     Q = Q1 >= Q0 if JOBO = 'E', where Q1 is the observability
C     Grammian defined to enforce stability for a modified combination
C     method of [4].
C
C     If JOBC = JOBO = 'S' and ALPHAC = ALPHAO = 0, the choice of
C     Grammians corresponds to the method of Enns [1], while if
C     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds to the
C     method of Lin and Chiu [2,3].
C
C     The routine computes directly the Cholesky factors S and R
C     such that P = S*S' and Q = R'*R according to formulas
C     developed in [4]. No matrix inversions are involved.
C
C     REFERENCES
C
C     [1] Enns, D.
C         Model reduction with balanced realizations: An error bound
C         and a frequency weighted generalization.
C         Proc. CDC, Las Vegas, pp. 127-132, 1984.
C
C     [2] Lin, C.-A. and Chiu, T.-Y.
C         Model reduction via frequency-weighted balanced realization.
C         Control Theory and Advanced Technology, vol. 8,
C         pp. 341-351, 1992.
C
C     [3] Sreeram, V., Anderson, B.D.O and Madievski, A.G.
C         New results on frequency weighted balanced reduction
C         technique.
C         Proc. ACC, Seattle, Washington, pp. 4004-4009, 1995.
C
C     [4] Varga, A. and Anderson, B.D.O.
C         Square-root balancing-free methods for the frequency-weighted
C         balancing related model reduction.
C         (report in preparation)
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000.
C     D. Sima, University of Bucharest, August 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000.
C
C     REVISIONS
C
C     A. Varga, Australian National University, Canberra, November 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000.
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2001.
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
      CHARACTER        DICO, JOBC, JOBO, WEIGHT
      INTEGER          INFO, LDA, LDAV, LDAW, LDB, LDBV, LDBW,
     $                 LDC, LDCV, LDCW, LDDV, LDDW, LDR, LDS, LDWORK,
     $                 M, MW, N, NV, NW, P, PV
      DOUBLE PRECISION ALPHAC, ALPHAO, SCALEC, SCALEO
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), AV(LDAV,*), AW(LDAW,*),
     $                 B(LDB,*), BV(LDBV,*), BW(LDBW,*),
     $                 C(LDC,*), CV(LDCV,*), CW(LDCW,*),
     $                           DV(LDDV,*), DW(LDDW,*),
     $                 DWORK(*), R(LDR,*),   S(LDS,*)
C     .. Local Scalars ..
      LOGICAL          DISCR, FRWGHT, LEFTW, RIGHTW
      INTEGER          I, IERR, J, KAW, KTAU, KU, KW, LDU, LW, MBBAR,
     $                 NNV, NNW, PCBAR
      DOUBLE PRECISION T, TOL, WORK
C     .. Local Arrays ..
      DOUBLE PRECISION DUM(1)
C     .. External Functions ..
      LOGICAL          LSAME
      DOUBLE PRECISION DLAMCH
      EXTERNAL         DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL         DCOPY, DGEMM, DLACPY, DLASET, DSCAL, DSYEV,
     $                 MB01WD, MB04ND, MB04OD, SB03OU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        ABS, DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      DISCR  = LSAME( DICO,   'D' )
      LEFTW  = LSAME( WEIGHT, 'L' ) .OR. LSAME( WEIGHT, 'B' )
      RIGHTW = LSAME( WEIGHT, 'R' ) .OR. LSAME( WEIGHT, 'B' )
      FRWGHT = LEFTW .OR. RIGHTW
C
      INFO = 0
      LW   = 1
      NNV  = N + NV
      NNW  = N + NW
      IF( LEFTW .AND. PV.GT.0 ) THEN
         LW = MAX( LW, NNV*( NNV + MAX( NNV, PV ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( P + 5 ) )
      END IF
      IF( RIGHTW .AND. MW.GT.0 ) THEN
         LW = MAX( LW, NNW*( NNW + MAX( NNW, MW ) + 5 ) )
      ELSE
         LW = MAX( LW, N*( M + 5 ) )
      END IF
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( JOBC, 'S' ) .OR. LSAME( JOBC, 'E' ) ) )
     $     THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( JOBO, 'S' ) .OR. LSAME( JOBO, 'E' ) ) )
     $     THEN
         INFO = -3
      ELSE IF( .NOT.( FRWGHT .OR. LSAME( WEIGHT, 'N' ) ) ) THEN
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
      ELSE IF( NW.LT.0 ) THEN
         INFO = -10
      ELSE IF( MW.LT.0 ) THEN
         INFO = -11
      ELSE IF( ABS( ALPHAC ).GT.ONE  ) THEN
         INFO = -12
      ELSE IF( ABS( ALPHAO ).GT.ONE  ) THEN
         INFO = -13
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -19
      ELSE IF( LDAV.LT.1 .OR. ( LEFTW  .AND. LDAV.LT.NV ) ) THEN
         INFO = -21
      ELSE IF( LDBV.LT.1 .OR. ( LEFTW  .AND. LDBV.LT.NV ) ) THEN
         INFO = -23
      ELSE IF( LDCV.LT.1 .OR. ( LEFTW  .AND. LDCV.LT.PV ) ) THEN
         INFO = -25
      ELSE IF( LDDV.LT.1 .OR. ( LEFTW  .AND. LDDV.LT.PV ) ) THEN
         INFO = -27
      ELSE IF( LDAW.LT.1 .OR. ( RIGHTW .AND. LDAW.LT.NW ) ) THEN
         INFO = -29
      ELSE IF( LDBW.LT.1 .OR. ( RIGHTW .AND. LDBW.LT.NW ) ) THEN
         INFO = -31
      ELSE IF( LDCW.LT.1 .OR. ( RIGHTW .AND. LDCW.LT.M  ) ) THEN
         INFO = -33
      ELSE IF( LDDW.LT.1 .OR. ( RIGHTW .AND. LDDW.LT.M  ) ) THEN
         INFO = -35
      ELSE IF( LDS.LT.MAX( 1, N ) ) THEN
         INFO = -39
      ELSE IF( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -41
      ELSE IF( LDWORK.LT.LW ) THEN
         INFO = -43
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09IY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      SCALEC = ONE
      SCALEO = ONE
      IF( MIN( N, M, P ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      WORK = 1
      IF( LEFTW .AND. PV.GT.0 ) THEN
C
C        Build the extended permuted matrices
C
C           Ao = ( Av  Bv*C ) ,   Co = ( Cv Dv*C ) .
C                ( 0     A  )
C
         KAW = 1
         KU  = KAW + NNV*NNV
         LDU = MAX( NNV, PV )
         CALL DLACPY( 'Full', NV, NV, AV, LDAV, DWORK(KAW), NNV )
         CALL DLASET( 'Full', N, NV, ZERO, ZERO, DWORK(KAW+NV), NNV )
         CALL DGEMM( 'No-transpose', 'No-transpose', NV, N, P, ONE,
     $               BV, LDBV, C, LDC, ZERO, DWORK(KAW+NNV*NV), NNV )
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(KAW+NNV*NV+NV), NNV )
C
         CALL DLACPY( 'Full', PV, NV, CV, LDCV, DWORK(KU), LDU  )
         CALL DGEMM( 'No-transpose', 'No-transpose', PV, N, P, ONE,
     $               DV, LDDV, C, LDC, ZERO, DWORK(KU+LDU*NV), LDU )
C
C        Solve for the Cholesky factor Ro of Qo, Qo = Ro'*Ro,
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            Ao'*Qo + Qo*Ao  +  scaleo^2*Co'*Co = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0.
C
C        Workspace:  need   (N+NV)*(N+NV+MAX(N+NV,PV)+5);
C                           prefer larger.
C
         KTAU = KU + LDU*NNV
         KW   = KTAU + NNV
C
         CALL SB03OU( DISCR, .FALSE., NNV, PV, DWORK(KAW), NNV,
     $                DWORK(KU), LDU, DWORK(KTAU), DWORK(KU), LDU,
     $                SCALEO, DWORK(KW), LDWORK-KW+1, IERR )
C
         IF( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C        Partition Ro as Ro = ( R11 R12 ) and compute R such that
C                             (  0  R22 )
C
C        R'*R = R22'*R22 + (1-ALPHAO**2)*R12'*R12.
C
         KW = KU + LDU*NV + NV
         CALL DLACPY( 'Upper', N, N, DWORK(KW), LDU, R, LDR )
         IF( ALPHAO.NE.ZERO ) THEN
            T = SQRT( ONE - ALPHAO*ALPHAO )
            DO 10 J = KU + LDU*NV, KU + LDU*(NNV-1), LDU
               CALL DSCAL( NV, T, DWORK(J), 1 )
   10       CONTINUE
         END IF
         IF( ALPHAO.LT.ONE .AND. NV.GT.0 ) THEN
            KTAU = 1
            CALL MB04OD( 'Full', N, 0, NV, R, LDR, DWORK(KU+LDU*NV),
     $                   LDU, DUM, 1, DUM, 1, DWORK(KTAU), DWORK(KW) )
C
            DO 30 J = 1, N
               DWORK(J) = R(J,J)
               DO 20 I = 1, J
                  IF ( DWORK(I).LT.ZERO ) R(I,J) = -R(I,J)
   20          CONTINUE
   30       CONTINUE
C
         END IF
C
         IF( LSAME( JOBO, 'E' ) .AND. ALPHAO.LT.ONE ) THEN
C
C           Form Y = -A'*(R'*R)-(R'*R)*A if DICO = 'C', or
C                Y = -A'*(R'*R)*A+(R'*R) if DICO = 'D'.
C
            CALL DLACPY( 'Upper', N, N, R, LDR, DWORK(KU), N )
            CALL MB01WD( DICO, 'Upper', 'No-transpose', 'Hessenberg', N,
     $                   -ONE, ZERO, R, LDR, DWORK(KAW+NNV*NV+NV), NNV,
     $                   DWORK(KU), N, IERR )
C
C           Compute the eigendecomposition of Y as Y = Z*Sigma*Z'.
C
            KU = N + 1
            CALL DSYEV( 'Vectors', 'Upper', N, R, LDR, DWORK, DWORK(KU),
     $                  LDWORK-N, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 3
               RETURN
            END IF
            WORK = MAX( WORK, DWORK(KU) + DBLE( N ) )
C
C           Partition Sigma = (Sigma1,Sigma2), such that
C           Sigma1 <= 0, Sigma2 > 0.
C           Partition correspondingly Z = [Z1 Z2].
C
            TOL = MAX( ABS( DWORK(1) ), ABS( DWORK(N) ) )
     $            * DLAMCH( 'Epsilon')
C                _
C           Form C = [ sqrt(Sigma2)*Z2' ]
C
            PCBAR = 0
            DO 40 J = 1, N
               IF( DWORK(J).GT.TOL ) THEN
                  CALL DSCAL( N, SQRT( DWORK(J) ), R(1,J), 1 )
                  CALL DCOPY( N, R(1,J), 1, DWORK(KU+PCBAR), N )
                  PCBAR = PCBAR + 1
               END IF
   40       CONTINUE
C
C           Solve for the Cholesky factor R of Q, Q = R'*R,
C           the continuous-time Lyapunov equation (if DICO = 'C')
C                                      _  _
C                   A'*Q + Q*A  +  t^2*C'*C = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C                                      _  _
C                   A'*Q*A - Q  +  t^2*C'*C = 0.
C
C           Workspace:  need   N*(N + 6);
C                              prefer larger.
C
            KTAU = KU + N*N
            KW   = KTAU + N
C
            CALL SB03OU( DISCR, .FALSE., N, PCBAR, A, LDA, DWORK(KU), N,
     $                   DWORK(KTAU), R, LDR, T, DWORK(KW), LDWORK-KW+1,
     $                   IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            SCALEO = SCALEO*T
            WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
         END IF
C
      ELSE
C
C        Solve for the Cholesky factor R of Q, Q = R'*R,
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            A'*Q + Q*A  +  scaleo^2*C'*C = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            A'*Q*A - Q +  scaleo^2*C'*C = 0.
C
C        Workspace:  need   N*(P + 5);
C                           prefer larger.
C
         KU   = 1
         KTAU = KU + P*N
         KW   = KTAU + N
C
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), P )
         CALL SB03OU( DISCR, .FALSE., N, P, A, LDA, DWORK(KU), P,
     $                DWORK(KTAU), R, LDR, SCALEO, DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
      END IF
C
      IF( RIGHTW .AND. MW.GT.0 ) THEN
C
C        Build the extended matrices
C
C           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) .
C                ( 0   Aw  )          (  Bw  )
C
         KAW = 1
         KU  = KAW + NNW*NNW
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(KAW), NNW )
         CALL DLASET( 'Full', NW, N, ZERO, ZERO, DWORK(KAW+N), NNW )
         CALL DGEMM( 'No-transpose', 'No-transpose', N, NW, M, ONE,
     $               B, LDB, CW, LDCW, ZERO, DWORK(KAW+NNW*N), NNW )
         CALL DLACPY( 'Full', NW, NW, AW, LDAW,
     $                DWORK(KAW+NNW*N+N), NNW )
C
         CALL DGEMM( 'No-transpose', 'No-transpose', N, MW, M, ONE,
     $               B, LDB, DW, LDDW, ZERO, DWORK(KU), NNW  )
         CALL DLACPY( 'Full', NW, MW, BW, LDBW, DWORK(KU+N), NNW )
C
C        Solve for the Cholesky factor Si of Pi, Pi = Si*Si',
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0.
C
C        Workspace:  need   (N+NW)*(N+NW+MAX(N+NW,MW)+5);
C                           prefer larger.
C
         KTAU = KU + NNW*MAX( NNW, MW )
         KW   = KTAU + NNW
C
         CALL SB03OU( DISCR, .TRUE., NNW, MW, DWORK(KAW), NNW,
     $                DWORK(KU), NNW, DWORK(KTAU), DWORK(KU), NNW,
     $                SCALEC, DWORK(KW), LDWORK-KW+1, IERR )
C
         IF( IERR.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
C
C        Partition Si as Si = ( S11 S12 ) and compute S such that
C                             (  0  S22 )
C
C        S*S' = S11*S11' + (1-ALPHAC**2)*S12*S12'.
C
         CALL DLACPY( 'Upper', N, N, DWORK(KU), NNW, S, LDS )
         IF( ALPHAC.NE.ZERO ) THEN
            T = SQRT( ONE - ALPHAC*ALPHAC )
            DO 50 J = KU + NNW*N, KU + NNW*(NNW-1), NNW
               CALL DSCAL( N, T, DWORK(J), 1 )
   50       CONTINUE
         END IF
         IF( ALPHAC.LT.ONE .AND. NW.GT.0 ) THEN
            KTAU = N*NNW + 1
            KW   = KTAU  + N
            CALL MB04ND( 'Full', N, 0, NW, S, LDS, DWORK(KU+NNW*N), NNW,
     $                   DUM, 1, DUM, 1, DWORK(KTAU), DWORK(KW) )
C
            DO 70 J = 1, N
               IF ( S(J,J).LT.ZERO ) THEN
                  DO 60 I = 1, J
                     S(I,J) = -S(I,J)
   60             CONTINUE
               END IF
   70       CONTINUE
         END IF
C
         IF( LSAME( JOBC, 'E' ) .AND. ALPHAC.LT.ONE ) THEN
C
C           Form X = -A*(S*S')-(S*S')*A' if DICO = 'C', or
C                X = -A*(S*S')*A'+(S*S') if DICO = 'D'.
C
            CALL DLACPY( 'Upper', N, N, S, LDS, DWORK(KU), N )
            CALL MB01WD( DICO, 'Upper', 'Transpose', 'Hessenberg', N,
     $                   -ONE, ZERO, S, LDS, DWORK(KAW), NNW, DWORK(KU),
     $                   N, IERR )
C
C           Compute the eigendecomposition of X as X = Z*Sigma*Z'.
C
            KU = N + 1
            CALL DSYEV( 'Vectors', 'Upper', N, S, LDS, DWORK, DWORK(KU),
     $                  LDWORK-N, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 3
               RETURN
            END IF
            WORK = MAX( WORK, DWORK(KU) + DBLE( N ) )
C
C           Partition Sigma = (Sigma1,Sigma2), such that
C           Sigma1 =< 0, Sigma2 > 0.
C           Partition correspondingly Z = [Z1 Z2].
C
            TOL = MAX( ABS( DWORK(1) ), ABS( DWORK(N) ) )
     $            * DLAMCH( 'Epsilon')
C                _
C           Form B = [ Z2*sqrt(Sigma2) ]
C
            MBBAR = 0
            I = KU
            DO 80 J = 1, N
               IF( DWORK(J).GT.TOL ) THEN
                  MBBAR = MBBAR + 1
                  CALL DSCAL( N, SQRT( DWORK(J) ), S(1,J), 1 )
                  CALL DCOPY( N, S(1,J), 1, DWORK(I), 1 )
                  I = I + N
               END IF
   80       CONTINUE
C
C           Solve for the Cholesky factor S of P, P = S*S',
C           the continuous-time Lyapunov equation (if DICO = 'C')
C                                      _ _
C                   A*P + P*A'  +  t^2*B*B' = 0,
C
C           or the discrete-time Lyapunov equation (if DICO = 'D')
C                                      _ _
C                   A*P*A' - P  +  t^2*B*B' = 0.
C
C           Workspace:  need   maximum N*(N + 6);
C                              prefer larger.
C
            KTAU = KU + MBBAR*N
            KW   = KTAU + N
C
            CALL SB03OU( DISCR, .TRUE., N, MBBAR, A, LDA, DWORK(KU), N,
     $                   DWORK(KTAU), S, LDS, T, DWORK(KW), LDWORK-KW+1,
     $                   IERR )
            IF( IERR.NE.0 ) THEN
               INFO = 2
               RETURN
            END IF
            SCALEC = SCALEC*T
            WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
         END IF
C
      ELSE
C
C        Solve for the Cholesky factor S of P, P = S*S',
C        the continuous-time Lyapunov equation (if DICO = 'C')
C
C            A*P + P*A' +  scalec^2*B*B' = 0,
C
C        or the discrete-time Lyapunov equation (if DICO = 'D')
C
C            A*P*A' - P +  scalec^2*B*B' = 0.
C
C        Workspace:  need   N*(M+5);
C                           prefer larger.
C
         KU   = 1
         KTAU = KU + N*M
         KW   = KTAU + N
C
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KU), N )
         CALL SB03OU( DISCR, .TRUE., N, M, A, LDA, DWORK(KU), N,
     $                DWORK(KTAU), S, LDS, SCALEC, DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
         WORK = MAX( WORK, DWORK(KW) + DBLE( KW - 1 ) )
      END IF
C
C     Save optimal workspace.
C
      DWORK(1) = WORK
C
      RETURN
C *** Last line of AB09IY ***
      END
