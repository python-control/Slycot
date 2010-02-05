      SUBROUTINE SB10YD( DISCFL, FLAG, LENDAT, RFRDAT, IFRDAT, OMEGA, N,
     $                   A, LDA, B, C, D, TOL, IWORK, DWORK, LDWORK,
     $                   ZWORK, LZWORK, INFO )
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
C     To fit a supplied frequency response data with a stable, minimum
C     phase SISO (single-input single-output) system represented by its
C     matrices A, B, C, D. It handles both discrete- and continuous-time
C     cases.
C
C     ARGUMENTS
C
C     Input/Output parameters
C
C     DISCFL  (input) INTEGER
C             Indicates the type of the system, as follows:
C             = 0: continuous-time system;
C             = 1: discrete-time system.
C
C     FLAG    (input) INTEGER
C             If FLAG = 0, then the system zeros and poles are not
C             constrained.
C             If FLAG = 1, then the system zeros and poles will have
C             negative real parts in the continuous-time case, or moduli
C             less than 1 in the discrete-time case. Consequently, FLAG
C             must be equal to 1 in mu-synthesis routines.
C
C     LENDAT  (input) INTEGER
C             The length of the vectors RFRDAT, IFRDAT and OMEGA.
C             LENDAT >= 2.
C
C     RFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT)
C             The real part of the frequency data to be fitted.
C
C     IFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT)
C             The imaginary part of the frequency data to be fitted.
C
C     OMEGA   (input) DOUBLE PRECISION array, dimension (LENDAT)
C             The frequencies corresponding to RFRDAT and IFRDAT.
C             These values must be nonnegative and monotonically
C             increasing. Additionally, for discrete-time systems
C             they must be between 0 and PI.
C
C     N       (input/output) INTEGER
C             On entry, the desired order of the system to be fitted.
C             N <= LENDAT-1.
C             On exit, the order of the obtained system. The value of N
C             could only be modified if N > 0 and FLAG = 1.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the
C             matrix A. If FLAG = 1, then A is in an upper Hessenberg
C             form, and corresponds to a minimal realization.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (output) DOUBLE PRECISION array, dimension (N)
C             The computed vector B.
C
C     C       (output) DOUBLE PRECISION array, dimension (N)
C             The computed vector C. If FLAG = 1, the first N-1 elements
C             are zero (for the exit value of N).
C
C     D       (output) DOUBLE PRECISION array, dimension (1)
C             The computed scalar D.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used for determining the effective
C             rank of matrices. If the user sets TOL > 0, then the given
C             value of TOL is used as a lower bound for the reciprocal
C             condition number;  a (sub)matrix whose estimated condition
C             number is less than 1/TOL is considered to be of full
C             rank.  If the user sets TOL <= 0, then an implicitly
C             computed, default tolerance, defined by TOLDEF = SIZE*EPS,
C             is used instead, where SIZE is the product of the matrix
C             dimensions, and EPS is the machine precision (see LAPACK
C             Library routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(2,2*N+1)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK and DWORK(2) contains the optimal value of
C             LZWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK = max( 2, LW1, LW2, LW3, LW4 ), where
C             LW1 = 2*LENDAT + 4*HNPTS;  HNPTS = 2048;
C             LW2 =   LENDAT + 6*HNPTS;
C             MN  = min( 2*LENDAT, 2*N+1 )
C             LW3 = 2*LENDAT*(2*N+1) + max( 2*LENDAT, 2*N+1 ) +
C                   max( MN + 6*N + 4, 2*MN + 1 ), if N > 0;
C             LW3 = 4*LENDAT + 5                 , if N = 0;
C             LW4 = max( N*N + 5*N, 6*N + 1 + min( 1,N ) ), if FLAG = 1;
C             LW4 = 0,                                      if FLAG = 0.
C             For optimum performance LDWORK should be larger.
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK = LENDAT*(2*N+3), if N > 0;
C             LZWORK = LENDAT,         if N = 0.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the discrete --> continuous transformation cannot
C                   be made;
C             = 2:  if the system poles cannot be found;
C             = 3:  if the inverse system cannot be found, i.e., D is
C                   (close to) zero;
C             = 4:  if the system zeros cannot be found;
C             = 5:  if the state-space representation of the new
C                   transfer function T(s) cannot be found;
C             = 6:  if the continuous --> discrete transformation cannot
C                   be made.
C
C     METHOD
C
C     First, if the given frequency data are corresponding to a
C     continuous-time system, they are changed to a discrete-time
C     system using a bilinear transformation with a scaled alpha.
C     Then, the magnitude is obtained from the supplied data.
C     Then, the frequency data are linearly interpolated around
C     the unit-disc.
C     Then, Oppenheim and Schafer complex cepstrum method is applied
C     to get frequency data corresponding to a stable, minimum-
C     phase system. This is done in the following steps:
C     - Obtain LOG (magnitude)
C     - Obtain IFFT of the result (DG01MD SLICOT subroutine);
C     - halve the data at 0;
C     - Obtain FFT of the halved data (DG01MD SLICOT subroutine);
C     - Obtain EXP of the result.
C     Then, the new frequency data are interpolated back to the
C     original frequency.
C     Then, based on these newly obtained data, the system matrices
C     A, B, C, D are constructed; the very identification is
C     performed by Least Squares Method using DGELSY LAPACK subroutine.
C     If needed, a discrete-to-continuous time transformation is
C     applied on the system matrices by AB04MD SLICOT subroutine.
C     Finally, if requested, the poles and zeros of the system are
C     checked. If some of them have positive real parts in the
C     continuous-time case (or are not inside the unit disk in the
C     complex plane in the discrete-time case), they are exchanged with
C     their negatives (or reciprocals, respectively), to preserve the
C     frequency response, while getting a minimum phase and stable
C     system. This is done by SB10ZP SLICOT subroutine.
C
C     REFERENCES
C
C     [1] Oppenheim, A.V. and Schafer, R.W.
C         Discrete-Time Signal Processing.
C         Prentice-Hall Signal Processing Series, 1989.
C
C     [2] Balas, G., Doyle, J., Glover, K., Packard, A., and Smith, R.
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
C     A. Markovski, Technical University of Sofia, October 2003.
C
C     KEYWORDS
C
C     Bilinear transformation, frequency response, least-squares
C     approximation, stability.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16         ZZERO, ZONE
      PARAMETER          ( ZZERO = ( 0.0D+0, 0.0D+0 ),
     $                     ZONE  = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR, TEN
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     FOUR = 4.0D+0, TEN = 1.0D+1 )
      INTEGER            HNPTS
      PARAMETER          ( HNPTS = 2048 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            DISCFL, FLAG, INFO, LDA, LDWORK, LENDAT,
     $                   LZWORK, N
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK(*)
      DOUBLE PRECISION   A(LDA, *), B(*), C(*), D(*), DWORK(*),
     $                   IFRDAT(*), OMEGA(*), RFRDAT(*)
      COMPLEX*16         ZWORK(*)
C     ..
C     .. Local Scalars ..
      INTEGER            CLWMAX, DLWMAX, I, II, INFO2, IP1, IP2, ISTART,
     $                   ISTOP, IWA0, IWAB, IWBMAT, IWBP, IWBX, IWDME,
     $                   IWDOMO, IWMAG, IWS, IWVAR, IWXI, IWXR, IWYMAG,
     $                   K, LW1, LW2, LW3, LW4, MN, N1, N2, P, RANK
      DOUBLE PRECISION   P1, P2, PI, PW, RAT, TOLB, TOLL
      COMPLEX*16         XHAT(HNPTS/2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
C     ..
C     .. External Subroutines ..
      EXTERNAL           AB04MD, DCOPY, DG01MD, DGELSY, DLASET, DSCAL,
     $                   SB10ZP, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ACOS, ATAN, COS, DBLE, DCMPLX, DIMAG, EXP, LOG,
     $                   MAX, MIN, SIN, SQRT
C
C     Test input parameters and workspace.
C
      PI = FOUR*ATAN( ONE )
      PW = OMEGA(1)
      N1 = N + 1
      N2 = N + N1
C
      INFO = 0
      IF( DISCFL.NE.0 .AND. DISCFL.NE.1 ) THEN
         INFO = -1
      ELSE IF( FLAG.NE.0 .AND. FLAG.NE.1 ) THEN
         INFO = -2
      ELSE IF ( LENDAT.LT.2 ) THEN
         INFO = -3
      ELSE IF ( PW.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( N.GT.LENDAT - 1 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
C
         DO 10 K = 2, LENDAT
            IF ( OMEGA(K).LT.PW )
     $         INFO = -6
            PW = OMEGA(K)
   10    CONTINUE
C
         IF ( DISCFL.EQ.1 .AND. OMEGA(LENDAT).GT.PI )
     $      INFO = -6
      END IF
C
      IF ( INFO.EQ.0 ) THEN
C
C        Workspace.
C
         LW1 = 2*LENDAT + 4*HNPTS
         LW2 =   LENDAT + 6*HNPTS
         MN  = MIN( 2*LENDAT, N2 )
C
         IF ( N.GT.0 ) THEN
            LW3 = 2*LENDAT*N2 + MAX( 2*LENDAT, N2 ) +
     $                          MAX( MN + 6*N + 4, 2*MN + 1 )
         ELSE
            LW3 = 4*LENDAT + 5
         END IF
C
         IF ( FLAG.EQ.0 ) THEN
            LW4 = 0
         ELSE
            LW4 = MAX( N*N + 5*N, 6*N + 1 + MIN ( 1, N ) )
         END IF
C
         DLWMAX = MAX( 2, LW1, LW2, LW3, LW4 )
C
         IF ( N.GT.0 ) THEN
            CLWMAX = LENDAT*( N2 + 2 )
         ELSE
            CLWMAX = LENDAT
         END IF
C
         IF ( LDWORK.LT.DLWMAX ) THEN
            INFO = -16
         ELSE IF ( LZWORK.LT.CLWMAX ) THEN
            INFO = -18
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB10YD', -INFO )
         RETURN
      END IF
C
C     Set tolerances.
C
      TOLB = DLAMCH( 'Epsilon' )
      TOLL = TOL
      IF ( TOLL.LE.ZERO )
     $   TOLL = FOUR*DBLE( LENDAT*N )*TOLB
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Workspace usage 1.
C     Workspace:  need  2*LENDAT + 4*HNPTS.
C
      IWDOMO = 1
      IWDME  = IWDOMO + LENDAT
      IWYMAG = IWDME  + 2*HNPTS
      IWMAG  = IWYMAG + 2*HNPTS
C
C     Bilinear transformation.
C
      IF ( DISCFL.EQ.0 ) THEN
         PW = SQRT( OMEGA(1)*OMEGA(LENDAT) + SQRT( TOLB ) )
C
         DO 20 K = 1, LENDAT
            DWORK(IWDME+K-1)  = ( OMEGA(K)/PW )**2
            DWORK(IWDOMO+K-1) =
     $         ACOS( ( ONE - DWORK(IWDME+K-1) )/
     $               ( ONE + DWORK(IWDME+K-1) ) )
   20    CONTINUE
C
      ELSE
         CALL DCOPY( LENDAT, OMEGA, 1, DWORK(IWDOMO), 1 )
      END IF
C
C     Linear interpolation.
C
      DO 30 K = 1, LENDAT
         DWORK(IWMAG+K-1) = DLAPY2( RFRDAT(K), IFRDAT(K) )
         DWORK(IWMAG+K-1) = ( ONE/LOG( TEN ) ) * LOG( DWORK(IWMAG+K-1) )
   30 CONTINUE
C
      DO 40 K = 1, HNPTS
         DWORK(IWDME+K-1)  = ( K - 1 )*PI/HNPTS
         DWORK(IWYMAG+K-1) = ZERO
C
         IF ( DWORK(IWDME+K-1).LT.DWORK(IWDOMO) ) THEN
            DWORK(IWYMAG+K-1) = DWORK(IWMAG)
         ELSE IF ( DWORK(IWDME+K-1).GE.DWORK(IWDOMO+LENDAT-1) ) THEN
            DWORK(IWYMAG+K-1) = DWORK(IWMAG+LENDAT-1)
         END IF
C
   40 CONTINUE
C
      DO 60 I = 2, LENDAT
         P1 = HNPTS*DWORK(IWDOMO+I-2)/PI + ONE
C
         IP1 = INT( P1 )
         IF ( DBLE( IP1 ).NE.P1 )
     $      IP1 = IP1 + 1
C
         P2 = HNPTS*DWORK(IWDOMO+I-1)/PI + ONE
C
         IP2 = INT( P2 )
         IF ( DBLE( IP2 ).NE.P2 )
     $      IP2 = IP2 + 1
C
         DO 50 P = IP1, IP2 - 1
            RAT = DWORK(IWDME+P-1) - DWORK(IWDOMO+I-2)
            RAT = RAT/( DWORK(IWDOMO+I-1) - DWORK(IWDOMO+I-2) )
            DWORK(IWYMAG+P-1) = ( ONE - RAT )*DWORK(IWMAG+I-2) +
     $                          RAT*DWORK(IWMAG+I-1)
   50    CONTINUE
C
   60 CONTINUE
C
      DO 70 K = 1, HNPTS
         DWORK(IWYMAG+K-1) = EXP( LOG( TEN )*DWORK(IWYMAG+K-1) )
   70 CONTINUE
C
C     Duplicate data around disc.
C
      DO 80 K = 1, HNPTS
         DWORK(IWDME+HNPTS+K-1)  = TWO*PI - DWORK(IWDME+HNPTS-K)
         DWORK(IWYMAG+HNPTS+K-1) = DWORK(IWYMAG+HNPTS-K)
   80 CONTINUE
C
C     Complex cepstrum to get min phase:
C     LOG (Magnitude)
C
      DO 90 K = 1, 2*HNPTS
         DWORK(IWYMAG+K-1) = TWO*LOG( DWORK(IWYMAG+K-1) )
   90 CONTINUE
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Workspace usage 2.
C     Workspace:  need  LENDAT + 6*HNPTS.
C
      IWXR = IWYMAG
      IWXI = IWMAG
C
      DO 100 K = 1, 2*HNPTS
         DWORK(IWXI+K-1) = ZERO
  100 CONTINUE
C
C     IFFT
C
      CALL DG01MD( 'I', 2*HNPTS, DWORK(IWXR), DWORK(IWXI), INFO2 )
C
C     Rescale, because DG01MD doesn't do it.
C
      CALL DSCAL( HNPTS, ONE/( TWO*HNPTS ), DWORK(IWXR), 1 )
      CALL DSCAL( HNPTS, ONE/( TWO*HNPTS ), DWORK(IWXI), 1 )
C
C     Halve the result at 0.
C
      DWORK(IWXR) = DWORK(IWXR)/TWO
      DWORK(IWXI) = DWORK(IWXI)/TWO
C
C     FFT
C
      CALL DG01MD( 'D', HNPTS, DWORK(IWXR), DWORK(IWXI), INFO2 )
C
C     Get the EXP of the result.
C
      DO 110 K = 1, HNPTS/2
         XHAT(K) = EXP( DWORK(IWXR+K-1) )*
     $         DCMPLX ( COS( DWORK(IWXI+K-1)), SIN( DWORK(IWXI+K-1) ) )
         DWORK(IWDME+K-1) = DWORK(IWDME+2*K-2)
  110 CONTINUE
C
C     Interpolate back to original frequency data.
C
      ISTART = 1
      ISTOP  = LENDAT
C
      DO 120 I = 1, LENDAT
         ZWORK(I) = ZZERO
         IF ( DWORK(IWDOMO+I-1).LE.DWORK(IWDME) ) THEN
            ZWORK(I) = XHAT(1)
            ISTART = I + 1
         ELSE IF ( DWORK(IWDOMO+I-1).GE.DWORK(IWDME+HNPTS/2-1) )
     $         THEN
            ZWORK(I) = XHAT(HNPTS/2)
            ISTOP = ISTOP - 1
         END IF
  120 CONTINUE
C
      DO 140 I = ISTART, ISTOP
         II = HNPTS/2
  130    CONTINUE
            IF ( DWORK(IWDME+II-1).GE.DWORK(IWDOMO+I-1) )
     $         P = II
            II = II - 1
         IF ( II.GT.0 )
     $      GOTO 130
         RAT = ( DWORK(IWDOMO+I-1) - DWORK(IWDME+P-2) )/
     $         ( DWORK(IWDME+P-1)  - DWORK(IWDME+P-2) )
         ZWORK(I) = RAT*XHAT(P) + ( ONE - RAT )*XHAT(P-1)
  140 CONTINUE
C
C     CASE N > 0.
C     This is the only allowed case in mu-synthesis subroutines.
C
      IF ( N.GT.0 ) THEN
C
C        Preparation for frequency identification.
C
C        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C        Complex workspace usage 1.
C        Complex workspace:  need  2*LENDAT + LENDAT*(N+1).
C
         IWA0  = 1 + LENDAT
         IWVAR = IWA0 + LENDAT*N1
C
         DO 150 K = 1, LENDAT
            IF ( DISCFL.EQ.0 ) THEN
               ZWORK(IWVAR+K-1) = DCMPLX( COS( DWORK(IWDOMO+K-1) ),
     $                                    SIN( DWORK(IWDOMO+K-1) ) )
            ELSE
               ZWORK(IWVAR+K-1) = DCMPLX( COS( OMEGA(K) ),
     $                                    SIN( OMEGA(K) ) )
            END IF
  150    CONTINUE
C
C        Array for DGELSY.
C
         DO 160 K = 1, N2
            IWORK(K) = 0
  160    CONTINUE
C
C        Constructing A0.
C
         DO 170 K = 1, LENDAT
            ZWORK(IWA0+N*LENDAT+K-1) = ZONE
  170    CONTINUE
C
         DO 190 I = 1, N
            DO 180 K = 1, LENDAT
               ZWORK(IWA0+(N-I)*LENDAT+K-1) =
     $            ZWORK(IWA0+(N1-I)*LENDAT+K-1)*ZWORK(IWVAR+K-1)
  180       CONTINUE
  190    CONTINUE
C
C        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C        Complex workspace usage 2.
C        Complex workspace:  need  2*LENDAT + LENDAT*(2*N+1).
C
         IWBP = IWVAR
         IWAB = IWBP + LENDAT
C
C        Constructing BP.
C
         DO 200 K = 1, LENDAT
            ZWORK(IWBP+K-1) = ZWORK(IWA0+K-1)*ZWORK(K)
  200    CONTINUE
C
C        Constructing AB.
C
         DO 220 I = 1, N
            DO 210 K = 1, LENDAT
               ZWORK(IWAB+(I-1)*LENDAT+K-1) = -ZWORK(K)*
     $             ZWORK(IWA0+I*LENDAT+K-1)
  210       CONTINUE
  220    CONTINUE
C
C        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C        Workspace usage 3.
C        Workspace:  need  LW3 = 2*LENDAT*(2*N+1) + max(2*LENDAT,2*N+1).
C
         IWBX = 1 + 2*LENDAT*N2
         IWS  = IWBX + MAX( 2*LENDAT, N2 )
C
C        Constructing AX.
C
         DO 240 I = 1, N1
            DO 230 K = 1, LENDAT
               DWORK(2*(I-1)*LENDAT+K) =
     $            DBLE( ZWORK(IWA0+(I-1)*LENDAT+K-1) )
               DWORK((2*I-1)*LENDAT+K) =
     $            DIMAG( ZWORK(IWA0+(I-1)*LENDAT+K-1) )
  230       CONTINUE
  240    CONTINUE
C
         DO 260 I = 1, N
            DO 250 K = 1, LENDAT
               DWORK(2*N1*LENDAT+2*(I-1)*LENDAT+K) =
     $            DBLE( ZWORK(IWAB+(I-1)*LENDAT+K-1) )
               DWORK(2*N1*LENDAT+(2*I-1)*LENDAT+K) =
     $            DIMAG( ZWORK(IWAB+(I-1)*LENDAT+K-1) )
  250       CONTINUE
  260    CONTINUE
C
C        Constructing BX.
C
         DO 270 K = 1, LENDAT
            DWORK(IWBX+K-1) = DBLE( ZWORK(IWBP+K-1) )
            DWORK(IWBX+LENDAT+K-1) = DIMAG( ZWORK(IWBP+K-1) )
  270    CONTINUE
C
C        Estimating X.
C        Workspace:  need    LW3 + max( MN+3*(2*N+1)+1, 2*MN+1 ),
C                            where MN = min( 2*LENDAT, 2*N+1 );
C                            prefer  larger.
C
         CALL DGELSY( 2*LENDAT, N2, 1, DWORK, 2*LENDAT, DWORK(IWBX),
     $                MAX( 2*LENDAT, N2 ), IWORK, TOLL, RANK,
     $                DWORK(IWS), LDWORK-IWS+1, INFO2 )
         DLWMAX = MAX( DLWMAX, INT( DWORK(IWS) + IWS - 1 ) )
C
C        Constructing A matrix.
C
         DO 280 K = 1, N
            A(K,1) = -DWORK(IWBX+N1+K-1)
  280    CONTINUE
C
         IF ( N.GT.1 )
     $      CALL DLASET( 'Full', N, N-1, ZERO, ONE, A(1,2), LDA )
C
C        Constructing B matrix.
C
         DO 290 K = 1, N
            B(K) = DWORK(IWBX+N1+K-1)*DWORK(IWBX) - DWORK(IWBX+K)
  290    CONTINUE
C
C        Constructing C matrix.
C
         C(1) = -ONE
C
         DO 300 K = 2, N
            C(K) = ZERO
  300    CONTINUE
C
C        Constructing D matrix.
C
         D(1) = DWORK(IWBX)
C
C        Transform to continuous-time case, if needed.
C        Workspace:  need    max(1,N);
C                            prefer  larger.
C
         IF ( DISCFL.EQ.0 ) THEN
            CALL AB04MD( 'D', N, 1, 1, ONE, PW, A, LDA, B, LDA, C, 1,
     $                   D, 1, IWORK, DWORK, LDWORK, INFO2 )
            IF ( INFO2.NE.0 ) THEN
               INFO = 1
               RETURN
            END IF
            DLWMAX = MAX( DLWMAX, INT( DWORK(1) ) )
         END IF
C
C        Make all the real parts of the poles and the zeros negative.
C
         IF ( FLAG.EQ.1 ) THEN
C
C           Workspace:  need    max(N*N + 5*N, 6*N + 1 + min(1,N));
C                               prefer  larger.
            CALL SB10ZP( DISCFL, N, A, LDA, B, C, D, IWORK, DWORK,
     $                   LDWORK, INFO )
            IF ( INFO.NE.0 )
     $         RETURN
            DLWMAX = MAX( DLWMAX, INT( DWORK(1) ) )
         END IF
C
       ELSE
C
C        CASE N = 0.
C
C        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C        Workspace usage 4.
C        Workspace:  need  4*LENDAT.
C
         IWBMAT = 1 + 2*LENDAT
         IWS    = IWBMAT + 2*LENDAT
C
C        Constructing AMAT and BMAT.
C
         DO 310 K = 1, LENDAT
            DWORK(K) = ONE
            DWORK(K+LENDAT) = ZERO
            DWORK(IWBMAT+K-1) = DBLE( ZWORK(K) )
            DWORK(IWBMAT+LENDAT+K-1) = DIMAG( ZWORK(K) )
  310    CONTINUE
C
C        Estimating D matrix.
C        Workspace:  need    4*LENDAT + 5;
C                            prefer  larger.
C
         IWORK(1) = 0
         CALL DGELSY( 2*LENDAT, 1, 1, DWORK, 2*LENDAT, DWORK(IWBMAT),
     $                2*LENDAT, IWORK, TOLL, RANK, DWORK(IWS),
     $                LDWORK-IWS+1, INFO2 )
         DLWMAX = MAX( DLWMAX, INT( DWORK(IWS) + IWS - 1 ) )
C
         D(1) = DWORK(IWBMAT)
C
      END IF
C
C     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      DWORK(1) = DLWMAX
      DWORK(2) = CLWMAX
      RETURN
C
C *** Last line of SB10YD ***
      END
