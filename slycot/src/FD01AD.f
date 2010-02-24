      SUBROUTINE FD01AD( JP, L, LAMBDA, XIN, YIN, EFOR, XF, EPSBCK,
     $                   CTETA, STETA, YQ, EPOS, EOUT, SALPH, IWARN,
     $                   INFO )
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
C     To solve the least-squares filtering problem recursively in time.
C     Each subroutine call implements one time update of the solution.
C     The algorithm uses a fast QR-decomposition based approach.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JP      CHARACTER*1
C             Indicates whether the user wishes to apply both prediction
C             and filtering parts, as follows:
C             = 'B':  Both prediction and filtering parts are to be
C                     applied;
C             = 'P':  Only the prediction section is to be applied.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The length of the impulse response of the equivalent
C             transversal filter model.  L >= 1.
C
C     LAMBDA  (input) DOUBLE PRECISION
C             Square root of the forgetting factor.
C             For tracking capabilities and exponentially stable error
C             propagation, LAMBDA < 1.0 (strict inequality) should
C             be used.  0.0 < LAMBDA <= 1.0.
C
C     XIN     (input) DOUBLE PRECISION
C             The input sample at instant n.
C             (The situation just before and just after the call of
C             the routine are denoted by instant (n-1) and instant n,
C             respectively.)
C
C     YIN     (input) DOUBLE PRECISION
C             If JP = 'B', then YIN must contain the reference sample
C             at instant n.
C             Otherwise, YIN is not referenced.
C
C     EFOR    (input/output) DOUBLE PRECISION
C             On entry, this parameter must contain the square root of
C             exponentially weighted forward prediction error energy
C             at instant (n-1).  EFOR >= 0.0.
C             On exit, this parameter contains the square root of the
C             exponentially weighted forward prediction error energy
C             at instant n.
C
C     XF      (input/output) DOUBLE PRECISION array, dimension (L)
C             On entry, this array must contain the transformed forward
C             prediction variables at instant (n-1).
C             On exit, this array contains the transformed forward
C             prediction variables at instant n.
C
C     EPSBCK  (input/output) DOUBLE PRECISION array, dimension (L+1)
C             On entry, the leading L elements of this array must
C             contain the normalized a posteriori backward prediction
C             error residuals of orders zero through L-1, respectively,
C             at instant (n-1), and EPSBCK(L+1) must contain the
C             square-root of the so-called "conversion factor" at
C             instant (n-1).
C             On exit, this array contains the normalized a posteriori
C             backward prediction error residuals, plus the square root
C             of the conversion factor at instant n.
C
C     CTETA   (input/output) DOUBLE PRECISION array, dimension (L)
C             On entry, this array must contain the cosines of the
C             rotation angles used in time updates, at instant (n-1).
C             On exit, this array contains the cosines of the rotation
C             angles at instant n.
C
C     STETA   (input/output) DOUBLE PRECISION array, dimension (L)
C             On entry, this array must contain the sines of the
C             rotation angles used in time updates, at instant (n-1).
C             On exit, this array contains the sines of the rotation
C             angles at instant n.
C
C     YQ      (input/output) DOUBLE PRECISION array, dimension (L)
C             On entry, if JP = 'B', then this array must contain the
C             orthogonally transformed reference vector at instant
C             (n-1). These elements are also the tap multipliers of an
C             equivalent normalized lattice least-squares filter.
C             Otherwise, YQ is not referenced and can be supplied as
C             a dummy array (i.e., declare this array to be YQ(1) in
C             the calling program).
C             On exit, if JP = 'B', then this array contains the
C             orthogonally transformed reference vector at instant n.
C
C     EPOS    (output) DOUBLE PRECISION
C             The a posteriori forward prediction error residual.
C
C     EOUT    (output) DOUBLE PRECISION
C             If JP = 'B', then EOUT contains the a posteriori output
C             error residual from the least-squares filter at instant n.
C
C     SALPH   (output) DOUBLE PRECISION array, dimension (L)
C             The element SALPH(i), i=1,...,L, contains the opposite of
C             the i-(th) reflection coefficient for the least-squares
C             normalized lattice predictor (whose value is -SALPH(i)).
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  an element to be annihilated by a rotation is less
C                   than the machine precision (see LAPACK Library
C                   routine DLAMCH).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The output error EOUT at instant n, denoted by EOUT(n), is the
C     reference sample minus a linear combination of L successive input
C     samples:
C
C                           L-1
C        EOUT(n) = YIN(n) - SUM h_i * XIN(n-i),
C                           i=0
C
C     where YIN(n) and XIN(n) are the scalar samples at instant n.
C     A least-squares filter uses those h_0,...,h_{L-1} which minimize
C     an exponentially weighted sum of successive output errors squared:
C
C         n
C        SUM [LAMBDA**(2(n-k)) * EOUT(k)**2].
C        k=1
C
C     Each subroutine call performs a time update of the least-squares
C     filter using a fast least-squares algorithm derived from a
C     QR decomposition, as described in references [1] and [2] (the
C     notation from [2] is followed in the naming of the arrays).
C     The algorithm does not compute the parameters h_0,...,h_{L-1} from
C     the above formula, but instead furnishes the parameters of an
C     equivalent normalized least-squares lattice filter, which are
C     available from the arrays SALPH (reflection coefficients) and YQ
C     (tap multipliers), as well as the exponentially weighted input
C     signal energy
C
C         n                                              L
C        SUM [LAMBDA**(2(n-k)) * XIN(k)**2] = EFOR**2 + SUM XF(i)**2.
C        k=1                                            i=1
C
C     For more details on reflection coefficients and tap multipliers,
C     references [2] and [4] are recommended.
C
C     REFERENCES
C
C     [1]  Proudler, I. K., McWhirter, J. G., and Shepherd, T. J.
C          Fast QRD based algorithms for least-squares linear
C          prediction.
C          Proceedings IMA Conf. Mathematics in Signal Processing
C          Warwick, UK, December 1988.
C
C     [2]  Regalia, P. A., and Bellanger, M. G.
C          On the duality between QR methods and lattice methods in
C          least-squares adaptive filtering.
C          IEEE Trans. Signal Processing, SP-39, pp. 879-891,
C          April 1991.
C
C     [3]  Regalia, P. A.
C          Numerical stability properties of a QR-based fast
C          least-squares algorithm.
C          IEEE Trans. Signal Processing, SP-41, June 1993.
C
C     [4]  Lev-Ari, H., Kailath, T., and Cioffi, J.
C          Least-squares adaptive lattice and transversal filters:
C          A unified geometric theory.
C          IEEE Trans. Information Theory, IT-30, pp. 222-236,
C          March 1984.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(L) operations for each subroutine call.
C     It is backward consistent for all input sequences XIN, and
C     backward stable for persistently exciting input sequences,
C     assuming LAMBDA < 1.0 (see [3]).
C     If the condition of the signal is very poor (IWARN = 1), then the
C     results are not guaranteed to be reliable.
C
C     FURTHER COMMENTS
C
C     1.  For tracking capabilities and exponentially stable error
C         propagation, LAMBDA < 1.0 should be used.  LAMBDA is typically
C         chosen slightly less than 1.0 so that "past" data are
C         exponentially forgotten.
C     2.  Prior to the first subroutine call, the variables must be
C         initialized. The following initial values are recommended:
C
C         XF(i) = 0.0,        i=1,...,L
C         EPSBCK(i) = 0.0     i=1,...,L
C         EPSBCK(L+1) = 1.0
C         CTETA(i) = 1.0      i=1,...,L
C         STETA(i) = 0.0      i=1,...,L
C         YQ(i) = 0.0         i=1,...,L
C
C         EFOR = 0.0          (exact start)
C         EFOR = "small positive constant" (soft start).
C
C         Soft starts are numerically more reliable, but result in a
C         biased least-squares solution during the first few iterations.
C         This bias decays exponentially fast provided LAMBDA < 1.0.
C         If sigma is the standard deviation of the input sequence
C         XIN, then initializing EFOR = sigma*1.0E-02 usually works
C         well.
C
C     CONTRIBUTOR
C
C     P. A. Regalia (October 1994).
C     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Kalman filtering, least-squares estimator, optimal filtering,
C     orthogonal transformation, recursive estimation, QR decomposition.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JP
      INTEGER           INFO, IWARN, L
      DOUBLE PRECISION  EFOR, EOUT, EPOS, LAMBDA, XIN, YIN
C     .. Array Arguments ..
      DOUBLE PRECISION  CTETA(*), EPSBCK(*), SALPH(*), STETA(*), XF(*),
     $                  YQ(*)
C     .. Local Scalars ..
      LOGICAL           BOTH
      INTEGER           I
      DOUBLE PRECISION  CTEMP, EPS, FNODE, NORM, TEMP, XFI, YQI
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLAPY2, DNRM2
      EXTERNAL          DLAMCH, DLAPY2, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DLARTG, XERBLA
C     .. Intrinsic Functions
      INTRINSIC         ABS, SQRT
C     .. Executable statements ..
C
C     Test the input scalar arguments.
C
      BOTH  = LSAME( JP, 'B' )
      IWARN = 0
      INFO  = 0
C
      IF( .NOT.BOTH .AND. .NOT.LSAME( JP, 'P' ) ) THEN
         INFO = -1
      ELSE IF( L.LT.1 ) THEN
         INFO = -2
      ELSE IF( ( LAMBDA.LE.ZERO ) .OR. ( LAMBDA.GT.ONE  ) ) THEN
         INFO = -3
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'FD01AD', -INFO )
         RETURN
      END IF
C
C     Computation of the machine precision EPS.
C
      EPS = DLAMCH( 'Epsilon' )
C
C     Forward prediction rotations.
C
      FNODE = XIN
C
      DO 10  I = 1, L
         XFI   = XF(I) * LAMBDA
         XF(I) = STETA(I) * FNODE + CTETA(I) * XFI
         FNODE = CTETA(I) * FNODE - STETA(I) * XFI
   10 CONTINUE
C
      EPOS = FNODE * EPSBCK(L+1)
C
C     Update the square root of the prediction energy.
C
      EFOR = EFOR * LAMBDA
      TEMP = DLAPY2( FNODE, EFOR )
      IF ( TEMP.LT.EPS ) THEN
         FNODE = ZERO
         IWARN = 1
      ELSE
         FNODE = FNODE * EPSBCK(L+1)/TEMP
      END IF
      EFOR = TEMP
C
C     Calculate the reflection coefficients and the backward prediction
C     errors.
C
      DO 20 I = L, 1, -1
         IF ( ABS( XF(I) ).LT.EPS )
     $      IWARN = 1
         CALL DLARTG( TEMP, XF(I), CTEMP, SALPH(I), NORM )
         EPSBCK(I+1) = CTEMP * EPSBCK(I) - SALPH(I) * FNODE
         FNODE = CTEMP * FNODE + SALPH(I) * EPSBCK(I)
         TEMP  = NORM
   20 CONTINUE
C
      EPSBCK(1) = FNODE
C
C     Update to new rotation angles.
C
      NORM = DNRM2( L, EPSBCK, 1 )
      TEMP = SQRT( ( ONE + NORM )*( ONE - NORM ) )
      EPSBCK(L+1) = TEMP
C
      DO 30 I = L, 1, -1
         IF ( ABS( EPSBCK(I) ).LT.EPS )
     $      IWARN = 1
         CALL DLARTG( TEMP, EPSBCK(I), CTETA(I), STETA(I), NORM )
         TEMP = NORM
   30 CONTINUE
C
C     Joint process section.
C
      IF ( BOTH) THEN
         FNODE = YIN
C
         DO 40  I = 1, L
            YQI   = YQ(I) * LAMBDA
            YQ(I) = STETA(I) * FNODE + CTETA(I) * YQI
            FNODE = CTETA(I) * FNODE - STETA(I) * YQI
   40    CONTINUE
C
         EOUT = FNODE * EPSBCK(L+1)
      END IF
C
      RETURN
C *** Last line of FD01AD ***
      END
