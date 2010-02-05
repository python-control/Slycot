      SUBROUTINE SB10MD( NC, MP, LENDAT, F, ORD, MNB, NBLOCK, ITYPE,
     $                   QUTOL, A, LDA, B, LDB, C, LDC, D, LDD, OMEGA,
     $                   TOTORD, AD, LDAD, BD, LDBD, CD, LDCD, DD, LDDD,
     $                   MJU, IWORK, LIWORK, DWORK, LDWORK, ZWORK,
     $                   LZWORK, INFO )
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
C     To perform the D-step in the D-K iteration. It handles
C     continuous-time case.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NC      (input) INTEGER
C             The order of the matrix A.  NC >= 0.
C
C     MP      (input) INTEGER
C             The order of the matrix D.  MP >= 0.
C
C     LENDAT  (input) INTEGER
C             The length of the vector OMEGA.  LENDAT >= 2.
C
C     F       (input) INTEGER
C             The number of the measurements and controls, i.e.,
C             the size of the block I_f in the D-scaling system.
C             F >= 0.
C
C     ORD     (input/output) INTEGER
C             The MAX order of EACH block in the fitting procedure.
C             ORD <= LENDAT-1.
C             On exit, if ORD < 1 then ORD = 1.
C
C     MNB     (input) INTEGER
C             The number of diagonal blocks in the block structure of
C             the uncertainty, and the length of the vectors NBLOCK
C             and ITYPE.  1 <= MNB <= MP.
C
C     NBLOCK  (input) INTEGER array, dimension (MNB)
C             The vector of length MNB containing the block structure
C             of the uncertainty. NBLOCK(I), I = 1:MNB, is the size of
C             each block.
C
C     ITYPE   (input) INTEGER array, dimension (MNB)
C             The vector of length MNB indicating the type of each
C             block.
C             For I = 1 : MNB,
C             ITYPE(I) = 1 indicates that the corresponding block is a
C             real block. IN THIS CASE ONLY MJU(JW) WILL BE ESTIMATED
C             CORRECTLY, BUT NOT D(S)!
C             ITYPE(I) = 2 indicates that the corresponding block is a
C             complex block. THIS IS THE ONLY ALLOWED VALUE NOW!
C             NBLOCK(I) must be equal to 1 if ITYPE(I) is equal to 1.
C
C     QUTOL   (input) DOUBLE PRECISION
C             The acceptable mean relative error between the D(jw) and
C             the frequency responce of the estimated block
C             [ADi,BDi;CDi,DDi]. When it is reached, the result is
C             taken as good enough.
C             A good value is QUTOL = 2.0.
C             If QUTOL < 0 then only mju(jw) is being estimated,
C             not D(s).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,NC)
C             On entry, the leading NC-by-NC part of this array must
C             contain the A matrix of the closed-loop system.
C             On exit, if MP > 0, the leading NC-by-NC part of this
C             array contains an upper Hessenberg matrix similar to A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,NC).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,MP)
C             On entry, the leading NC-by-MP part of this array must
C             contain the B matrix of the closed-loop system.
C             On exit, the leading NC-by-MP part of this array contains
C             the transformed B matrix corresponding to the Hessenberg
C             form of A.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,NC).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,NC)
C             On entry, the leading MP-by-NC part of this array must
C             contain the C matrix of the closed-loop system.
C             On exit, the leading MP-by-NC part of this array contains
C             the transformed C matrix corresponding to the Hessenberg
C             form of A.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,MP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,MP)
C             The leading MP-by-MP part of this array must contain the
C             D matrix of the closed-loop system.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,MP).
C
C     OMEGA   (input) DOUBLE PRECISION array, dimension (LENDAT)
C             The vector with the frequencies.
C
C     TOTORD  (output) INTEGER
C             The TOTAL order of the D-scaling system.
C             TOTORD is set to zero, if QUTOL < 0.
C
C     AD      (output) DOUBLE PRECISION array, dimension (LDAD,MP*ORD)
C             The leading TOTORD-by-TOTORD part of this array contains
C             the A matrix of the D-scaling system.
C             Not referenced if QUTOL < 0.
C
C     LDAD    INTEGER
C             The leading dimension of the array AD.
C             LDAD >= MAX(1,MP*ORD), if QUTOL >= 0;
C             LDAD >= 1,             if QUTOL <  0.
C
C     BD      (output) DOUBLE PRECISION array, dimension (LDBD,MP+F)
C             The leading TOTORD-by-(MP+F) part of this array contains
C             the B matrix of the D-scaling system.
C             Not referenced if QUTOL < 0.
C
C     LDBD    INTEGER
C             The leading dimension of the array BD.
C             LDBD >= MAX(1,MP*ORD), if QUTOL >= 0;
C             LDBD >= 1,             if QUTOL <  0.
C
C     CD      (output) DOUBLE PRECISION array, dimension (LDCD,MP*ORD)
C             The leading (MP+F)-by-TOTORD part of this array contains
C             the C matrix of the D-scaling system.
C             Not referenced if QUTOL < 0.
C
C     LDCD    INTEGER
C             The leading dimension of the array CD.
C             LDCD >= MAX(1,MP+F), if QUTOL >= 0;
C             LDCD >= 1,           if QUTOL <  0.
C
C     DD      (output) DOUBLE PRECISION array, dimension (LDDD,MP+F)
C             The leading (MP+F)-by-(MP+F) part of this array contains
C             the D matrix of the D-scaling system.
C             Not referenced if QUTOL < 0.
C
C     LDDD    INTEGER
C             The leading dimension of the array DD.
C             LDDD >= MAX(1,MP+F), if QUTOL >= 0;
C             LDDD >= 1,           if QUTOL <  0.
C
C     MJU     (output) DOUBLE PRECISION array, dimension (LENDAT)
C             The vector with the upper bound of the structured
C             singular value (mju) for each frequency in OMEGA.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C
C     LIWORK  INTEGER
C             The length of the array IWORK.
C             LIWORK >= MAX( NC, 4*MNB-2, MP, 2*ORD+1 ), if QUTOL >= 0;
C             LIWORK >= MAX( NC, 4*MNB-2, MP ),          if QUTOL <  0.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK, DWORK(2) returns the optimal value of LZWORK,
C             and DWORK(3) returns an estimate of the minimum reciprocal
C             of the condition numbers (with respect to inversion) of
C             the generated Hessenberg matrices.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 3, LWM, LWD ), where
C             LWM = LWA + MAX( NC + MAX( NC, MP-1 ),
C                              2*MP*MP*MNB - MP*MP + 9*MNB*MNB +
C                              MP*MNB + 11*MP + 33*MNB - 11 );
C             LWD = LWB + MAX( 2, LW1, LW2, LW3, LW4, 2*ORD ),
C                              if QUTOL >= 0;
C             LWD = 0,         if QUTOL <  0;
C             LWA = MP*LENDAT + 2*MNB + MP - 1;
C             LWB = LENDAT*(MP + 2) + ORD*(ORD + 2) + 1;
C             LW1 = 2*LENDAT + 4*HNPTS;  HNPTS = 2048;
C             LW2 =   LENDAT + 6*HNPTS;  MN  = MIN( 2*LENDAT, 2*ORD+1 );
C             LW3 = 2*LENDAT*(2*ORD + 1) + MAX( 2*LENDAT, 2*ORD + 1 ) +
C                   MAX( MN + 6*ORD + 4, 2*MN + 1 );
C             LW4 = MAX( ORD*ORD + 5*ORD, 6*ORD + 1 + MIN( 1, ORD ) ).
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK >= MAX( LZM, LZD ), where
C             LZM = MAX( MP*MP + NC*MP + NC*NC + 2*NC,
C                        6*MP*MP*MNB + 13*MP*MP + 6*MNB + 6*MP - 3 );
C             LZD = MAX( LENDAT*(2*ORD + 3), ORD*ORD + 3*ORD + 1 ),
C                              if QUTOL >= 0;
C             LZD = 0,         if QUTOL <  0.
C
C     Error indicator
C
C     INFO    (output) INTEGER
C             =  0:  successful exit;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  if one or more values w in OMEGA are (close to
C                    some) poles of the closed-loop system, i.e., the
C                    matrix jw*I - A is (numerically) singular;
C             =  2:  the block sizes must be positive integers;
C             =  3:  the sum of block sizes must be equal to MP;
C             =  4:  the size of a real block must be equal to 1;
C             =  5:  the block type must be either 1 or 2;
C             =  6:  errors in solving linear equations or in matrix
C                    inversion;
C             =  7:  errors in computing eigenvalues or singular values.
C             = 1i:  INFO on exit from SB10YD is i. (1i means 10 + i.)
C
C     METHOD
C
C     I.   First, W(jw) for the given closed-loop system is being
C          estimated.
C     II.  Now, AB13MD SLICOT subroutine can obtain the D(jw) scaling
C          system with respect to NBLOCK and ITYPE, and colaterally,
C          mju(jw).
C          If QUTOL < 0 then the estimations stop and the routine exits.
C     III. Now that we have D(jw), SB10YD subroutine can do block-by-
C          block fit. For each block it tries with an increasing order
C          of the fit, starting with 1 until the
C          (mean quadratic error + max quadratic error)/2
C          between the Dii(jw) and the estimated frequency responce
C          of the block becomes less than or equal to the routine
C          argument QUTOL, or the order becomes equal to ORD.
C     IV.  Arrange the obtained blocks in the AD, BD, CD and DD
C          matrices and estimate the total order of D(s), TOTORD.
C     V.   Add the system I_f to the system obtained in IV.
C
C     REFERENCES
C
C     [1] Balas, G., Doyle, J., Glover, K., Packard, A. and Smith, R.
C         Mu-analysis and Synthesis toolbox - User's Guide,
C         The Mathworks Inc., Natick, MA, USA, 1998.
C
C     CONTRIBUTORS
C
C     Asparuh Markovski, Technical University of Sofia, July 2003.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003.
C     A. Markovski, V. Sima, October 2003.
C
C     KEYWORDS
C
C     Frequency response, H-infinity optimal control, robust control,
C     structured singular value.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO  = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     THREE = 3.0D+0 )
      INTEGER            HNPTS
      PARAMETER          ( HNPTS = 2048 )
C     ..
C     .. Scalar Arguments ..
      INTEGER           F, INFO, LDA, LDAD, LDB, LDBD, LDC, LDCD, LDD,
     $                  LDDD, LDWORK, LENDAT, LIWORK, LZWORK, MNB, MP,
     $                  NC, ORD, TOTORD
      DOUBLE PRECISION  QUTOL
C     ..
C     .. Array Arguments ..
      INTEGER           ITYPE(*), IWORK(*), NBLOCK(*)
      DOUBLE PRECISION  A(LDA, *), AD(LDAD, *), B(LDB, *), BD(LDBD, *),
     $                  C(LDC, *), CD(LDCD, *), D(LDD, *), DD(LDDD, *),
     $                  DWORK(*), MJU(*), OMEGA(*)
      COMPLEX*16        ZWORK(*)
C     ..
C     .. Local Scalars ..
      CHARACTER         BALEIG, INITA
      INTEGER           CLWMAX, CORD, DLWMAX, I, IC, ICWRK, IDWRK, II,
     $                  INFO2, IWAD, IWB, IWBD, IWCD, IWDD, IWGJOM,
     $                  IWIFRD, IWRFRD, IWX, K, LCSIZE, LDSIZE, LORD,
     $                  LW1, LW2, LW3, LW4, LWA, LWB, MAXCWR, MAXWRK,
     $                  MN, W
      DOUBLE PRECISION  MAQE, MEQE, MOD1, MOD2, RCND, RCOND, RQE, TOL,
     $                  TOLER
      COMPLEX*16        FREQ
C     ..
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL          AB13MD, DCOPY, DLACPY, DLASET, DSCAL, SB10YD,
     $                  TB05AD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX, INT, MAX, MIN, SQRT
C
C     Decode and test input parameters.
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     Workspace usage 1.
C
C     real
C
      IWX    = 1 + MP*LENDAT
      IWGJOM = IWX + 2*MNB - 1
      IDWRK  = IWGJOM + MP
      LDSIZE = LDWORK - IDWRK + 1
C
C     complex
C
      IWB    = MP*MP + 1
      ICWRK  = IWB + NC*MP
      LCSIZE = LZWORK - ICWRK + 1
C
      INFO = 0
      IF ( NC.LT.0 ) THEN
         INFO = -1
      ELSE IF( MP.LT.0 ) THEN
         INFO = -2
      ELSE IF( LENDAT.LT.2 ) THEN
         INFO = -3
      ELSE IF( F.LT.0 ) THEN
         INFO = -4
      ELSE IF( ORD.GT.LENDAT - 1 ) THEN
         INFO = -5
      ELSE IF( MNB.LT.1 .OR. MNB.GT.MP ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, NC ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, NC ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, MP ) ) THEN
         INFO = -15
      ELSE IF( LDD.LT.MAX( 1, MP ) ) THEN
         INFO = -17
      ELSE IF( LDAD.LT.1 .OR. ( QUTOL.GE.ZERO .AND. LDAD.LT.MP*ORD ) )
     $      THEN
         INFO = -21
      ELSE IF( LDBD.LT.1 .OR. ( QUTOL.GE.ZERO .AND. LDBD.LT.MP*ORD ) )
     $      THEN
         INFO = -23
      ELSE IF( LDCD.LT.1 .OR. ( QUTOL.GE.ZERO .AND. LDCD.LT.MP + F ) )
     $      THEN
         INFO = -25
      ELSE IF( LDDD.LT.1 .OR. ( QUTOL.GE.ZERO .AND. LDDD.LT.MP + F ) )
     $      THEN
         INFO = -27
      ELSE
C
C        Compute workspace.
C
         II  = MAX( NC, 4*MNB - 2, MP )
         MN  = MIN( 2*LENDAT, 2*ORD + 1 )
         LWA = IDWRK - 1
         LWB = LENDAT*( MP + 2 ) + ORD*( ORD + 2 ) + 1
         LW1 = 2*LENDAT + 4*HNPTS
         LW2 =   LENDAT + 6*HNPTS
         LW3 = 2*LENDAT*( 2*ORD + 1 ) + MAX( 2*LENDAT, 2*ORD + 1 ) +
     $                                  MAX( MN + 6*ORD + 4, 2*MN + 1 )
         LW4 = MAX( ORD*ORD + 5*ORD, 6*ORD + 1 + MIN( 1, ORD ) )
C
         DLWMAX = LWA + MAX( NC + MAX( NC, MP - 1 ),
     $                       2*MP*MP*MNB - MP*MP + 9*MNB*MNB + MP*MNB +
     $                       11*MP + 33*MNB - 11 )
C
         CLWMAX = MAX( ICWRK - 1 + NC*NC + 2*NC,
     $                 6*MP*MP*MNB + 13*MP*MP + 6*MNB + 6*MP - 3 )
C
         IF ( QUTOL.GE.ZERO ) THEN
            II     = MAX( II, 2*ORD + 1 )
            DLWMAX = MAX( DLWMAX,
     $                    LWB + MAX( 2, LW1, LW2, LW3, LW4, 2*ORD ) )
            CLWMAX = MAX( CLWMAX, LENDAT*( 2*ORD + 3 ),
     $                    ORD*( ORD + 3 ) + 1 )
         END IF
         IF ( LIWORK.LT.II ) THEN
            INFO = -30
         ELSE IF ( LDWORK.LT.MAX( 3, DLWMAX ) ) THEN
            INFO = -32
         ELSE IF ( LZWORK.LT.CLWMAX ) THEN
            INFO = -34
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB10MD', -INFO )
         RETURN
      END IF
C
      ORD    = MAX( 1, ORD )
      TOTORD = 0
C
C     Quick return if possible.
C
      IF( NC.EQ.0 .OR. MP.EQ.0 ) THEN
         DWORK(1) = THREE
         DWORK(2) = ZERO
         DWORK(3) = ONE
         RETURN
      END IF
C
      TOLER = SQRT( DLAMCH( 'Epsilon' ) )
C
      BALEIG = 'C'
      RCOND  = ONE
      MAXCWR = CLWMAX
C
C     @@@ 1. Estimate W(jw) for the closed-loop system, @@@
C     @@@      D(jw) and mju(jw) for each frequency.    @@@
C
      DO 30 W = 1, LENDAT
         FREQ = DCMPLX( ZERO, OMEGA(W) )
         IF ( W.EQ.1 ) THEN
            INITA = 'G'
         ELSE
            INITA = 'H'
         END IF
C
C        Compute C*inv(jw*I-A)*B.
C        Integer workspace: need   NC.
C        Real workspace:    need   LWA + NC + MAX(NC,MP-1);
C                           prefer larger,
C                           where  LWA = MP*LENDAT + 2*MNB + MP - 1.
C        Complex workspace: need   MP*MP + NC*MP + NC*NC + 2*NC.
C
         CALL TB05AD( BALEIG, INITA, NC, MP, MP, FREQ, A, LDA, B, LDB,
     $                C, LDC, RCND, ZWORK, MP, DWORK, DWORK, ZWORK(IWB),
     $                NC, IWORK, DWORK(IDWRK), LDSIZE, ZWORK(ICWRK),
     $                LCSIZE, INFO2 )
C
         IF ( INFO2.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
C
         RCOND = MIN( RCOND, RCND )
         IF ( W.EQ.1 )
     $      MAXWRK = INT( DWORK(IDWRK) + IDWRK - 1 )
         IC = 0
C
C        D + C*inv(jw*I-A)*B
C
         DO 20 K = 1, MP
            DO 10 I = 1, MP
               IC = IC + 1
               ZWORK(IC) = ZWORK(IC) + DCMPLX ( D(I,K), ZERO )
   10       CONTINUE
   20    CONTINUE
C
C        Estimate D(jw) and mju(jw).
C        Integer workspace: need   MAX(4*MNB-2,MP).
C        Real workspace:    need   LWA + 2*MP*MP*MNB - MP*MP + 9*MNB*MNB
C                                  + MP*MNB + 11*MP + 33*MNB - 11;
C                           prefer larger.
C        Complex workspace: need   6*MP*MP*MNB + 13*MP*MP + 6*MNB +
C                                  6*MP - 3.
C
         CALL AB13MD( 'N', MP, ZWORK, MP, MNB, NBLOCK, ITYPE,
     $                DWORK(IWX), MJU(W), DWORK((W-1)*MP+1),
     $                DWORK(IWGJOM), IWORK, DWORK(IDWRK), LDSIZE,
     $                ZWORK(IWB), LZWORK-IWB+1, INFO2 )
C
         IF ( INFO2.NE.0 ) THEN
            INFO = INFO2 + 1
            RETURN
         END IF
C
         IF ( W.EQ.1 ) THEN
            MAXWRK = MAX( MAXWRK, INT( DWORK(IDWRK) ) + IDWRK - 1 )
            MAXCWR = MAX( MAXCWR, INT( ZWORK(IWB) ) + IWB - 1 )
         END IF
C
C        Normalize D(jw) through it's last entry.
C
         IF ( DWORK(W*MP).NE.ZERO )
     $      CALL DSCAL( MP, ONE/DWORK(W*MP), DWORK((W-1)*MP+1), 1 )
C
   30 CONTINUE
C
C     Quick return if needed.
C
      IF ( QUTOL.LT.ZERO ) THEN
         DWORK(1) = MAXWRK
         DWORK(2) = MAXCWR
         DWORK(3) = RCOND
         RETURN
      END IF
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     Workspace usage 2.
C
C     real
C
      IWRFRD = IWX
      IWIFRD = IWRFRD + LENDAT
      IWAD   = IWIFRD + LENDAT
      IWBD   = IWAD + ORD*ORD
      IWCD   = IWBD + ORD
      IWDD   = IWCD + ORD
      IDWRK  = IWDD + 1
      LDSIZE = LDWORK - IDWRK + 1
C
C     complex
C
      ICWRK  = ORD + 2
      LCSIZE = LZWORK - ICWRK + 1
      INITA  = 'H'
C
C     Use default tolerance for SB10YD.
C
      TOL = -ONE
C
C     @@@ 2. Clear imag parts of D(jw) for SB10YD. @@@
C
      DO 40 I = 1, LENDAT
         DWORK(IWIFRD+I-1) = ZERO
   40 CONTINUE
C
C     @@@ 3. Clear AD, BD, CD and initialize DD with I_(mp+f). @@@
C
      CALL DLASET( 'Full', MP*ORD, MP*ORD, ZERO, ZERO, AD, LDAD )
      CALL DLASET( 'Full', MP*ORD, MP+F,   ZERO, ZERO, BD, LDBD )
      CALL DLASET( 'Full', MP+F,   MP*ORD, ZERO, ZERO, CD, LDCD )
      CALL DLASET( 'Full', MP+F,   MP+F,   ZERO, ONE,  DD, LDDD )
C
C     @@@ 4. Block by block frequency identification. @@@
C
      DO 80 II = 1, MP
C
         CALL DCOPY( LENDAT, DWORK(II), MP, DWORK(IWRFRD), 1 )
C
C        Increase CORD from 1 to ORD for every block, if needed.
C
         CORD = 1
C
   50    CONTINUE
            LORD = CORD
C
C           Now, LORD is the desired order.
C           Integer workspace: need   2*N+1, where N = LORD.
C           Real workspace:    need   LWB + MAX( 2, LW1, LW2, LW3, LW4),
C                                     where
C                                     LWB = LENDAT*(MP+2) +
C                                           ORD*(ORD+2) + 1,
C                                     HNPTS = 2048, and
C                                     LW1 = 2*LENDAT + 4*HNPTS;
C                                     LW2 =   LENDAT + 6*HNPTS;
C                                     MN  = min( 2*LENDAT, 2*N+1 )
C                                     LW3 = 2*LENDAT*(2*N+1) +
C                                           max( 2*LENDAT, 2*N+1 ) +
C                                           max( MN + 6*N + 4, 2*MN+1 );
C                                     LW4 = max( N*N + 5*N,
C                                                6*N + 1 + min( 1,N ) );
C                              prefer larger.
C           Complex workspace: need   LENDAT*(2*N+3).
C
            CALL SB10YD( 0, 1, LENDAT, DWORK(IWRFRD), DWORK(IWIFRD),
     $                   OMEGA, LORD, DWORK(IWAD), ORD, DWORK(IWBD),
     $                   DWORK(IWCD), DWORK(IWDD), TOL, IWORK,
     $                   DWORK(IDWRK), LDSIZE, ZWORK, LZWORK, INFO2 )
C
C           At this point, LORD is the actual order reached by SB10YD,
C           0 <= LORD <= CORD.
C           [ADi,BDi; CDi,DDi] is a minimal realization with ADi in
C           upper Hessenberg form.
C           The leading LORD-by-LORD part of ORD-by-ORD DWORK(IWAD)
C           contains ADi, the leading LORD-by-1 part of ORD-by-1
C           DWORK(IWBD) contains BDi, the leading 1-by-LORD part of
C           1-by-ORD DWORK(IWCD) contains CDi, DWORK(IWDD) contains DDi.
C
            IF ( INFO2.NE.0 ) THEN
               INFO = 10 + INFO2
               RETURN
            END IF
C
C          Compare the original D(jw) with the fitted one.
C
            MEQE = ZERO
            MAQE = ZERO
C
            DO 60 W = 1, LENDAT
               FREQ = DCMPLX( ZERO, OMEGA(W) )
C
C              Compute CD*inv(jw*I-AD)*BD.
C              Integer workspace: need   LORD.
C              Real workspace:    need   LWB + 2*LORD;
C                                 prefer larger.
C              Complex workspace: need   1 + ORD + LORD*LORD + 2*LORD.
C
               CALL TB05AD( BALEIG, INITA, LORD, 1, 1, FREQ,
     $                      DWORK(IWAD), ORD, DWORK(IWBD), ORD,
     $                      DWORK(IWCD), 1, RCND, ZWORK, 1,
     $                      DWORK(IDWRK), DWORK(IDWRK), ZWORK(2), ORD,
     $                      IWORK, DWORK(IDWRK), LDSIZE, ZWORK(ICWRK),
     $                      LCSIZE, INFO2 )
C
               IF ( INFO2.GT.0 ) THEN
                  INFO = 1
                  RETURN
               END IF
C
               RCOND = MIN( RCOND, RCND )
               IF ( W.EQ.1 )
     $            MAXWRK = MAX( MAXWRK, INT( DWORK(IDWRK) ) + IDWRK - 1)
C
C              DD + CD*inv(jw*I-AD)*BD
C
               ZWORK(1) = ZWORK(1) + DCMPLX( DWORK(IWDD), ZERO )
C
               MOD1 = ABS( DWORK(IWRFRD+W-1) )
               MOD2 = ABS( ZWORK(1) )
               RQE  = ABS( ( MOD1 - MOD2 )/( MOD1 + TOLER ) )
               MEQE = MEQE + RQE
               MAQE = MAX( MAQE, RQE )
C
   60       CONTINUE
C
            MEQE = MEQE/LENDAT
C
            IF ( ( ( MEQE + MAQE )/TWO.LE.QUTOL ) .OR.
     $           ( CORD.EQ.ORD ) ) THEN
               GOTO 70
            END IF
C
            CORD = CORD + 1
         GOTO 50
C
   70    TOTORD = TOTORD + LORD
C
C        Copy ad(ii), bd(ii) and cd(ii) to AD, BD and CD, respectively.
C
         CALL DLACPY( 'Full', LORD, LORD, DWORK(IWAD), ORD,
     $                AD(TOTORD-LORD+1,TOTORD-LORD+1), LDAD )
         CALL DCOPY(  LORD, DWORK(IWBD), 1, BD(TOTORD-LORD+1,II), 1 )
         CALL DCOPY(  LORD, DWORK(IWCD), 1, CD(II,TOTORD-LORD+1), LDCD )
C
C        Copy dd(ii) to DD.
C
         DD(II,II) = DWORK(IWDD)
C
   80 CONTINUE
C
      DWORK(1) = MAXWRK
      DWORK(2) = MAXCWR
      DWORK(3) = RCOND
      RETURN
C
C *** Last line of SB10MD ***
      END
