      SUBROUTINE MB04YD( JOBU, JOBV, M, N, RANK, THETA, Q, E, U, LDU, V,
     $                   LDV, INUL, TOL, RELTOL, DWORK, LDWORK, IWARN,
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
C     To partially diagonalize the bidiagonal matrix
C
C               |q(1) e(1)  0    ...       0      |
C               | 0   q(2) e(2)            .      |
C           J = | .                        .      |                  (1)
C               | .                  e(MIN(M,N)-1)|
C               | 0   ...        ...  q(MIN(M,N)) |
C
C     using QR or QL iterations in such a way that J is split into
C     unreduced bidiagonal submatrices whose singular values are either
C     all larger than a given bound or are all smaller than (or equal
C     to) this bound. The left- and right-hand Givens rotations
C     performed on J (corresponding to each QR or QL iteration step) may
C     be optionally accumulated in the arrays U and V.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBU    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix U the left-hand Givens rotations, as follows:
C             = 'N':  Do not form U;
C             = 'I':  U is initialized to the M-by-MIN(M,N) submatrix of
C                     the unit matrix and the left-hand Givens rotations
C                     are accumulated in U;
C             = 'U':  The given matrix U is updated by the left-hand
C                     Givens rotations used in the calculation.
C
C     JOBV    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix V the right-hand Givens rotations, as follows:
C             = 'N':  Do not form V;
C             = 'I':  V is initialized to the N-by-MIN(M,N) submatrix of
C                     the unit matrix and the right-hand Givens
C                     rotations are accumulated in V;
C             = 'U':  The given matrix V is updated by the right-hand
C                     Givens rotations used in the calculation.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in matrix U.  M >= 0.
C
C     N       (input) INTEGER
C             The number of rows in matrix V.  N >= 0.
C
C     RANK    (input/output) INTEGER
C             On entry, if RANK < 0, then the rank of matrix J is
C             computed by the routine as the number of singular values
C             larger than THETA.
C             Otherwise, RANK must specify the rank of matrix J.
C             RANK <= MIN(M,N).
C             On exit, if RANK < 0 on entry, then RANK contains the
C             computed rank of J. That is, the number of singular
C             values of J larger than THETA.
C             Otherwise, the user-supplied value of RANK may be
C             changed by the routine on exit if the RANK-th and the
C             (RANK+1)-th singular values of J are considered to be
C             equal. See also the parameter TOL.
C
C     THETA   (input/output) DOUBLE PRECISION
C             On entry, if RANK < 0, then THETA must specify an upper
C             bound on the smallest singular values of J. THETA >= 0.0.
C             Otherwise, THETA must specify an initial estimate (t say)
C             for computing an upper bound such that precisely RANK
C             singular values are greater than this bound.
C             If THETA < 0.0, then t is computed by the routine.
C             On exit, if RANK >= 0 on entry, then THETA contains the
C             computed upper bound such that precisely RANK singular
C             values of J are greater than THETA + TOL.
C             Otherwise, THETA is unchanged.
C
C     Q       (input/output) DOUBLE PRECISION array, dimension
C             (MIN(M,N))
C             On entry, this array must contain the diagonal elements
C             q(1),q(2),...,q(MIN(M,N)) of the bidiagonal matrix J. That
C             is, Q(i) = J(i,i) for i = 1,2,...,MIN(M,N).
C             On exit, this array contains the leading diagonal of the
C             transformed bidiagonal matrix J.
C
C     E       (input/output) DOUBLE PRECISION array, dimension
C             (MIN(M,N)-1)
C             On entry, this array must contain the superdiagonal
C             elements e(1),e(2),...,e(MIN(M,N)-1) of the bidiagonal
C             matrix J. That is, E(k) = J(k,k+1) for k = 1,2,...,
C             MIN(M,N)-1.
C             On exit, this array contains the superdiagonal of the
C             transformed bidiagonal matrix J.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*)
C             On entry, if JOBU = 'U', the leading M-by-MIN(M,N) part
C             of this array must contain a left transformation matrix
C             applied to the original matrix of the problem, and
C             on exit, the leading M-by-MIN(M,N) part of this array
C             contains the product of the input matrix U and the
C             left-hand Givens rotations.
C             On exit, if JOBU = 'I', then the leading M-by-MIN(M,N)
C             part of this array contains the matrix of accumulated
C             left-hand Givens rotations used.
C             If JOBU = 'N', the array U is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDU = 1 and
C             declare this array to be U(1,1) in the calling program).
C
C     LDU     INTEGER
C             The leading dimension of array U. If JOBU = 'U' or
C             JOBU = 'I', LDU >= MAX(1,M); if JOBU = 'N', LDU >= 1.
C
C     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*)
C             On entry, if JOBV = 'U', the leading N-by-MIN(M,N) part
C             of this array must contain a right transformation matrix
C             applied to the original matrix of the problem, and
C             on exit, the leading N-by-MIN(M,N) part of this array
C             contains the product of the input matrix V and the
C             right-hand Givens rotations.
C             On exit, if JOBV = 'I', then the leading N-by-MIN(M,N)
C             part of this array contains the matrix of accumulated
C             right-hand Givens rotations used.
C             If JOBV = 'N', the array V is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDV = 1 and
C             declare this array to be V(1,1) in the calling program).
C
C     LDV     INTEGER
C             The leading dimension of array V. If JOBV = 'U' or
C             JOBV = 'I', LDV >= MAX(1,N); if JOBV = 'N', LDV >= 1.
C
C     INUL    (input/output) LOGICAL array, dimension (MIN(M,N))
C             On entry, the leading MIN(M,N) elements of this array must
C             be set to .FALSE. unless the i-th columns of U (if JOBU =
C             'U') and V (if JOBV = 'U') already contain a computed base
C             vector of the desired singular subspace of the original
C             matrix, in which case INUL(i) must be set to .TRUE.
C             for 1 <= i <= MIN(M,N).
C             On exit, the indices of the elements of this array with
C             value .TRUE. indicate the indices of the diagonal entries
C             of J which belong to those bidiagonal submatrices whose
C             singular values are all less than or equal to THETA.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             This parameter defines the multiplicity of singular values
C             by considering all singular values within an interval of
C             length TOL as coinciding. TOL is used in checking how many
C             singular values are less than or equal to THETA. Also in
C             computing an appropriate upper bound THETA by a bisection
C             method, TOL is used as a stopping criterion defining the
C             minimum (absolute) subinterval width. TOL is also taken
C             as an absolute tolerance for negligible elements in the
C             QR/QL iterations. If the user sets TOL to be less than or
C             equal to 0, then the tolerance is taken as
C             EPS * MAX(ABS(Q(i)), ABS(E(k))), where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH),
C             i = 1,2,...,MIN(M,N) and k = 1,2,...,MIN(M,N)-1.
C
C     RELTOL  DOUBLE PRECISION
C             This parameter specifies the minimum relative width of an
C             interval. When an interval is narrower than TOL, or than
C             RELTOL times the larger (in magnitude) endpoint, then it
C             is considered to be sufficiently small and bisection has
C             converged. If the user sets RELTOL to be less than
C             BASE * EPS, where BASE is machine radix and EPS is machine
C             precision (see LAPACK Library routine DLAMCH), then the
C             tolerance is taken as BASE * EPS.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,6*MIN(M,N)-5), if JOBU = 'I' or 'U', or
C                                               JOBV = 'I' or 'U';
C             LDWORK >= MAX(1,4*MIN(M,N)-3), if JOBU = 'N' and
C                                               JOBV = 'N'.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  if the rank of the bidiagonal matrix J (as specified
C                   by the user) has been lowered because a singular
C                   value of multiplicity larger than 1 was found.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value; this includes values like RANK > MIN(M,N), or
C                   THETA < 0.0 and RANK < 0;
C             = 1:  if the maximum number of QR/QL iteration steps
C                   (30*MIN(M,N)) has been exceeded.
C
C     METHOD
C
C     If the upper bound THETA is not specified by the user, then it is
C     computed by the routine (using a bisection method) such that
C     precisely (MIN(M,N) - RANK) singular values of J are less than or
C     equal to THETA + TOL.
C
C     The method used by the routine (see [1]) then proceeds as follows.
C
C     The unreduced bidiagonal submatrices of J(j), where J(j) is the
C     transformed bidiagonal matrix after the j-th iteration step, are
C     classified into the following three classes:
C
C     - C1 contains the bidiagonal submatrices with all singular values
C       > THETA,
C     - C2 contains the bidiagonal submatrices with all singular values
C       <= THETA and
C     - C3 contains the bidiagonal submatrices with singular values
C       > THETA and also singular values <= THETA.
C
C     If C3 is empty, then the partial diagonalization is complete, and
C     RANK is the sum of the dimensions of the bidiagonal submatrices of
C     C1.
C     Otherwise, QR or QL iterations are performed on each bidiagonal
C     submatrix of C3, until this bidiagonal submatrix has been split
C     into two bidiagonal submatrices. These two submatrices are then
C     classified and the iterations are restarted.
C     If the upper left diagonal element of the bidiagonal submatrix is
C     larger than its lower right diagonal element, then QR iterations
C     are performed, else QL iterations are used. The shift is taken as
C     the smallest diagonal element of the bidiagonal submatrix (in
C     magnitude) unless its value exceeds THETA, in which case it is
C     taken as zero.
C
C     REFERENCES
C
C     [1] Van Huffel, S., Vandewalle, J. and Haegemans, A.
C         An efficient and reliable algorithm for computing the
C         singular subspace of a matrix associated with its smallest
C         singular values.
C         J. Comput. and Appl. Math., 19, pp. 313-330, 1987.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     To avoid overflow, matrix J is scaled so that its largest element
C     is no greater than  overflow**(1/2) * underflow**(1/4) in absolute
C     value (and not much smaller than that, for maximal accuracy).
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB04QD by S. Van Huffel, Katholieke
C     University Leuven, Belgium.
C
C     REVISIONS
C
C     July 10, 1997. V. Sima.
C     November 25, 1997. V. Sima: Setting INUL(K) = .TRUE. when handling
C                                 2-by-2 submatrix.
C
C     KEYWORDS
C
C     Bidiagonal matrix, orthogonal transformation, singular values.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TEN, HNDRD
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0,
     $                    HNDRD = 100.0D0 )
      DOUBLE PRECISION  MEIGTH
      PARAMETER         ( MEIGTH = -0.125D0 )
      INTEGER           MAXITR
      PARAMETER         ( MAXITR = 30 )
C     .. Scalar Arguments ..
      CHARACTER         JOBU, JOBV
      INTEGER           INFO, IWARN, LDU, LDV, LDWORK, M, N, RANK
      DOUBLE PRECISION  RELTOL, THETA, TOL
C     .. Array Arguments ..
      LOGICAL           INUL(*)
      DOUBLE PRECISION  DWORK(*), E(*), Q(*), U(LDU,*), V(LDV,*)
C     .. Local Scalars ..
      LOGICAL           LJOBUA, LJOBUI, LJOBVA, LJOBVI, NOC12, QRIT
      INTEGER           I, I1, IASCL, INFO1, ITER, J, K, MAXIT, NUMEIG,
     $                  OLDI, OLDK, P, R
      DOUBLE PRECISION  COSL, COSR, EPS, PIVMIN, RMAX, RMIN, SAFEMN,
     $                  SHIFT, SIGMA, SIGMN, SIGMX, SINL, SINR, SMAX,
     $                  SMLNUM, THETAC, THRESH, TOLABS, TOLREL, X
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           MB03ND
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME, MB03ND
C     .. External Subroutines ..
      EXTERNAL          DLASET, DLASV2, DROT, DSCAL, MB02NY, MB03MD,
     $                  MB04YW, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      P = MIN( M, N )
      INFO = 0
      IWARN = 0
      LJOBUI = LSAME( JOBU, 'I' )
      LJOBVI = LSAME( JOBV, 'I' )
      LJOBUA = LJOBUI.OR.LSAME( JOBU, 'U' )
      LJOBVA = LJOBVI.OR.LSAME( JOBV, 'U' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBUA .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LJOBVA .AND. .NOT.LSAME( JOBV, 'N' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( RANK.GT.P ) THEN
         INFO = -5
      ELSE IF( RANK.LT.0 .AND. THETA.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( .NOT.LJOBUA .AND. LDU.LT.1 .OR.
     $              LJOBUA .AND. LDU.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( .NOT.LJOBVA .AND. LDV.LT.1 .OR.
     $              LJOBVA .AND. LDV.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( ( ( LJOBUA.OR.LJOBVA ) .AND. LDWORK.LT.MAX( 1, 6*P-5 ) )
     $ .OR.(.NOT.( LJOBUA.OR.LJOBVA ) .AND. LDWORK.LT.MAX( 1, 4*P-3 ) )
     $       ) THEN
         INFO = -17
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB04YD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( P.EQ.0 ) THEN
         IF ( RANK.GE.0 )
     $      THETA = ZERO
         RANK = 0
         RETURN
      END IF
C
C     Set tolerances and machine parameters.
C
      TOLABS = TOL
      TOLREL = RELTOL
      SMAX = ABS( Q(P) )
C
      DO 20 J = 1, P - 1
         SMAX = MAX( SMAX, ABS( Q(J) ), ABS( E(J) ) )
   20 CONTINUE
C
      SAFEMN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Epsilon' )
      IF ( TOLABS.LE.ZERO ) TOLABS = EPS*SMAX
      X = DLAMCH( 'Base' )*EPS
      IF ( TOLREL.LE.X ) TOLREL = X
      THRESH = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )*EPS
      SMLNUM = SAFEMN / EPS
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( ONE / RMIN, ONE / SQRT( SQRT( SAFEMN ) ) )
      THETAC = THETA
C
C     Scale the matrix to allowable range, if necessary, and set PIVMIN,
C     using the squares of Q and E (saved in DWORK).
C
      IASCL = 0
      IF( SMAX.GT.ZERO .AND. SMAX.LT.RMIN ) THEN
         IASCL = 1
         SIGMA = RMIN / SMAX
      ELSE IF( SMAX.GT.RMAX ) THEN
         IASCL = 1
         SIGMA = RMAX / SMAX
      END IF
      IF( IASCL.EQ.1 ) THEN
         CALL DSCAL( P, SIGMA, Q, 1 )
         CALL DSCAL( P-1, SIGMA, E, 1 )
         THETAC = SIGMA*THETA
         TOLABS = SIGMA*TOLABS
      END IF
C
      PIVMIN = Q(P)**2
      DWORK(P) = PIVMIN
C
      DO 40 J = 1, P - 1
         DWORK(J)   = Q(J)**2
         DWORK(P+J) = E(J)**2
         PIVMIN = MAX( PIVMIN, DWORK(J), DWORK(P+J) )
   40 CONTINUE
C
      PIVMIN = MAX( PIVMIN*SAFEMN, SAFEMN )
C
C     Initialize U and/or V to the identity matrix, if needed.
C
      IF ( LJOBUI )
     $   CALL DLASET( 'Full', M, P, ZERO, ONE, U, LDU )
      IF ( LJOBVI )
     $   CALL DLASET( 'Full', N, P, ZERO, ONE, V, LDV )
C
C     Estimate THETA (if not fixed by the user), and set R.
C
      IF ( RANK.GE.0 ) THEN
         J = P - RANK
         CALL MB03MD( P, J, THETAC, Q, E, DWORK(1), DWORK(P+1), PIVMIN,
     $                TOLABS, TOLREL, IWARN, INFO1 )
         THETA = THETAC
         IF ( IASCL.EQ.1 ) THETA = THETA / SIGMA
         IF ( J.LE.0 )
     $      RETURN
         R = P - J
      ELSE
         R = P - MB03ND( P, THETAC, DWORK, DWORK(P+1), PIVMIN, INFO1 )
      END IF
C
      RANK = P
C
      DO 60 I = 1, P
         IF ( INUL(I) ) RANK = RANK - 1
   60 CONTINUE
C
C     From now on K is the smallest known index such that the elements
C     of the bidiagonal matrix J with indices larger than K belong to C1
C     or C2.
C     RANK = P - SUM(dimensions of known bidiagonal matrices of C2).
C
      K = P
      OLDI = -1
      OLDK = -1
      ITER = 0
      MAXIT = MAXITR*P
C     WHILE ( C3 NOT EMPTY ) DO
   80 IF ( RANK.GT.R .AND. K.GT.0 ) THEN
C        WHILE ( K.GT.0 .AND. INUL(K) ) DO
C
C        Search for the rightmost index of a bidiagonal submatrix,
C        not yet classified.
C
  100    IF ( K.GT.0 ) THEN
            IF ( INUL(K) ) THEN
               K = K - 1
               GO TO 100
            END IF
         END IF
C        END WHILE 100
C
         IF ( K.EQ.0 )
     $      RETURN
C
         NOC12 = .TRUE.
C        WHILE ((ITER < MAXIT).AND.(No bidiagonal matrix of C1 or
C                C2 found)) DO
  120    IF ( ( ITER.LT.MAXIT ) .AND. NOC12 ) THEN
C
C           Search for negligible Q(I) or E(I-1) (for I > 1) and find
C           the shift.
C
            I = K
            X = ABS( Q(I) )
            SHIFT = X
C           WHILE ABS( Q(I) ) > TOLABS .AND. ABS( E(I-1) ) > TOLABS ) DO
  140       IF ( I.GT.1 ) THEN
               IF ( ( X.GT.TOLABS ).AND.( ABS( E(I-1) ).GT.TOLABS ) )
     $               THEN
                  I = I - 1
                  X = ABS( Q(I) )
                  IF ( X.LT.SHIFT ) SHIFT = X
                  GO TO 140
               END IF
            END IF
C           END WHILE 140
C
C           Classify the bidiagonal submatrix (of order J) found.
C
            J = K - I + 1
            IF ( ( X.LE.TOLABS ) .OR. ( K.EQ.I ) ) THEN
               NOC12 = .FALSE.
            ELSE
               NUMEIG = MB03ND( J, THETAC, DWORK(I), DWORK(P+I), PIVMIN,
     $                          INFO1 )
               IF ( NUMEIG.GE.J .OR. NUMEIG.LE.0 ) NOC12 = .FALSE.
            END IF
            IF ( NOC12 ) THEN
               IF ( J.EQ.2 ) THEN
C
C                 Handle separately the 2-by-2 submatrix.
C
                  CALL DLASV2( Q(I), E(I), Q(K), SIGMN, SIGMX, SINR,
     $                         COSR, SINL, COSL )
                  Q(I) = SIGMX
                  Q(K) = SIGMN
                  E(I) = ZERO
                  RANK = RANK - 1
                  INUL(K) = .TRUE.
                  NOC12 = .FALSE.
C
C                 Update U and/or V, if needed.
C
                  IF( LJOBUA )
     $               CALL DROT( M, U(1,I), 1, U(1,K), 1, COSL, SINL )
                  IF( LJOBVA )
     $               CALL DROT( N, V(1,I), 1, V(1,K), 1, COSR, SINR )
               ELSE
C
C                 If working on new submatrix, choose QR or
C                 QL iteration.
C
                  IF ( I.NE.OLDI .OR. K.NE.OLDK )
     $               QRIT = ABS( Q(I) ).GE.ABS( Q(K) )
                  OLDI = I
                  IF ( QRIT ) THEN
                     IF ( ABS( E(K-1) ).LE.THRESH*ABS( Q(K) ) )
     $                         E(K-1) = ZERO
                  ELSE
                     IF ( ABS( E(I) ).LE.THRESH*ABS( Q(I) ) )
     $                         E(I) = ZERO
                  END IF
C
                  CALL MB04YW( QRIT, LJOBUA, LJOBVA, M, N, I, K, SHIFT,
     $                         Q, E, U, LDU, V, LDV, DWORK(2*P) )
C
                  IF ( QRIT ) THEN
                     IF ( ABS( E(K-1) ).LE.TOLABS ) E(K-1) = ZERO
                  ELSE
                     IF ( ABS( E(I) ).LE.TOLABS ) E(I) = ZERO
                  END IF
                  DWORK(K) = Q(K)**2
C
                  DO 160 I1 = I, K - 1
                     DWORK(I1)   = Q(I1)**2
                     DWORK(P+I1) = E(I1)**2
  160             CONTINUE
C
                  ITER = ITER + 1
               END IF
            END IF
            GO TO 120
         END IF
C        END WHILE 120
C
         IF ( ITER.GE.MAXIT ) THEN
            INFO = 1
            GO TO 200
         END IF
C
         IF ( X.LE.TOLABS ) THEN
C
C           Split at negligible diagonal element ABS( Q(I) ) <= TOLABS.
C
            CALL MB02NY( LJOBUA, LJOBVA, M, N, I, K, Q, E, U, LDU, V,
     $                   LDV, DWORK(2*P) )
            INUL(I) = .TRUE.
            RANK = RANK - 1
         ELSE
C
C           A negligible superdiagonal element ABS( E(I-1) ) <= TOL
C           has been found, the corresponding bidiagonal submatrix
C           belongs to C1 or C2. Treat this bidiagonal submatrix.
C
            IF ( J.GE.2 ) THEN
               IF ( NUMEIG.EQ.J ) THEN
C
                  DO 180 I1 = I, K
                     INUL(I1) = .TRUE.
  180             CONTINUE
C
                  RANK = RANK - J
                  K = K - J
               ELSE
                  K = I - 1
               END IF
            ELSE
               IF ( X.LE.( THETAC + TOLABS ) ) THEN
                  INUL(I) = .TRUE.
                  RANK = RANK - 1
               END IF
               K = K - 1
            END IF
            OLDK = K
         END IF
         GO TO 80
      END IF
C     END WHILE 80
C
C     If matrix was scaled, then rescale Q and E appropriately.
C
  200 CONTINUE
      IF( IASCL.EQ.1 ) THEN
         CALL DSCAL( P,   ONE / SIGMA, Q, 1 )
         CALL DSCAL( P-1, ONE / SIGMA, E, 1 )
      END IF
C
      RETURN
C *** Last line of MB04YD ***
      END
