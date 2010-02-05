      SUBROUTINE MB05OD( BALANC, N, NDIAG, DELTA, A, LDA, MDIG, IDIG,
     $                   IWORK, DWORK, LDWORK, IWARN, INFO )
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
C     To compute exp(A*delta) where A is a real N-by-N matrix and delta
C     is a scalar value. The routine also returns the minimal number of
C     accurate digits in the 1-norm of exp(A*delta) and the number of
C     accurate digits in the 1-norm of exp(A*delta) at 95% confidence
C     level.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALANC  CHARACTER*1
C             Specifies whether or not a balancing transformation (done
C             by SLICOT Library routine MB04MD) is required, as follows:
C             = 'N', do not use balancing;
C             = 'S', use balancing (scaling).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     NDIAG   (input) INTEGER
C             The specified order of the diagonal Pade approximant.
C             In the absence of further information NDIAG should
C             be set to 9.  NDIAG should not exceed 15.  NDIAG >= 1.
C
C     DELTA   (input) DOUBLE PRECISION
C             The scalar value delta of the problem.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On input, the leading N-by-N part of this array must
C             contain the matrix A of the problem. (This is not needed
C             if DELTA = 0.)
C             On exit, if INFO = 0, the leading N-by-N part of this
C             array contains the solution matrix exp(A*delta).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     MDIG    (output) INTEGER
C             The minimal number of accurate digits in the 1-norm of
C             exp(A*delta).
C
C     IDIG    (output) INTEGER
C             The number of accurate digits in the 1-norm of
C             exp(A*delta) at 95% confidence level.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= N*(2*N+NDIAG+1)+NDIAG, if N >  1.
C             LDWORK >= 1,                     if N <= 1.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  if MDIG = 0 and IDIG > 0, warning for possible
C                   inaccuracy (the exponential has been computed);
C             = 2:  if MDIG = 0 and IDIG = 0, warning for severe
C                   inaccuracy (the exponential has been computed);
C             = 3:  if balancing has been requested, but it failed to
C                   reduce the matrix norm and was not actually used.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the norm of matrix A*delta (after a possible
C                   balancing) is too large to obtain an accurate
C                   result;
C             = 2:  if the coefficient matrix (the denominator of the
C                   Pade approximant) is exactly singular; try a
C                   different value of NDIAG;
C             = 3:  if the solution exponential would overflow, possibly
C                   due to a too large value DELTA; the calculations
C                   stopped prematurely. This error is not likely to
C                   appear.
C
C     METHOD
C
C     The exponential of the matrix A is evaluated from a diagonal Pade
C     approximant. This routine is a modification of the subroutine
C     PADE, described in reference [1]. The routine implements an
C     algorithm which exploits the identity
C
C         (exp[(2**-m)*A]) ** (2**m) = exp(A),
C
C     where m is an integer determined by the algorithm, to improve the
C     accuracy for matrices with large norms.
C
C     REFERENCES
C
C     [1] Ward, R.C.
C         Numerical computation of the matrix exponential with accuracy
C         estimate.
C         SIAM J. Numer. Anal., 14, pp. 600-610, 1977.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05CD by T.W.C. Williams, Kingston
C     Polytechnic, March 1982.
C
C     REVISIONS
C
C     June 14, 1997, April 25, 2003, December 12, 2004.
C
C     KEYWORDS
C
C     Continuous-time system, matrix algebra, matrix exponential,
C     matrix operations, Pade approximation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO, FOUR, EIGHT, TEN, TWELVE,
     $                  NINTEN, TWO4, FOUR7, TWOHND
      PARAMETER         ( ZERO = 0.0D0,    HALF = 0.5D0,   ONE = 1.0D0,
     $                    TWO  = 2.0D0,    FOUR = 4.0D0, EIGHT = 8.0D0,
     $                    TEN  = 10.0D0, TWELVE = 12.0D0,
     $                    NINTEN = 19.0D0, TWO4 = 24.0D0,
     $                    FOUR7  = 47.0D0, TWOHND = 200.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         BALANC
      INTEGER           IDIG, INFO, IWARN, LDA, LDWORK, MDIG, N,
     $                  NDIAG
      DOUBLE PRECISION  DELTA
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LBALS
      CHARACTER         ACTBAL
      INTEGER           BASE, I, IFAIL, IJ, IK, IM1, J, JWORA1, JWORA2,
     $                  JWORA3, JWORV1, JWORV2, K, M, MPOWER, NDAGM1,
     $                  NDAGM2, NDEC, NDECM1
      DOUBLE PRECISION  ANORM, AVGEV, BD, BIG, EABS, EAVGEV, EMNORM,
     $                  EPS, FACTOR, FN, GN, MAXRED, OVRTH2, OVRTHR, P,
     $                  RERL, RERR, S, SD2, SIZE, SMALL, SS, SUM2D,
     $                  TEMP, TMP1, TR, U, UNDERF, VAR, VAREPS, XN
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DASUM, DLAMCH, DLANGE, DNRM2
      EXTERNAL          DASUM, DLAMCH, DLANGE, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DGEMV, DGETRF, DGETRS, DLACPY,
     $                  DLASCL, DLASET, DSCAL, MB04MD, MB05OY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, EXP, INT, LOG, LOG10, MAX, MIN, MOD, SQRT
C     .. Executable Statements ..
C
      IWARN = 0
      INFO  = 0
      LBALS = LSAME( BALANC, 'S' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LBALS ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NDIAG.LT.1 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDWORK.LT.1 .OR.
     $       ( LDWORK.LT.N*( 2*N + NDIAG + 1 ) + NDIAG .AND. N.GT.1 )
     $       ) THEN
         INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB05OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      EPS  = DLAMCH( 'Epsilon' )
      NDEC = INT( LOG10( ONE/EPS ) + ONE )
C
      IF ( N.EQ.0 ) THEN
         MDIG = NDEC
         IDIG = NDEC
         RETURN
      END IF
C
C     Set some machine parameters.
C
      BASE = DLAMCH( 'Base' )
      NDECM1 = NDEC - 1
      UNDERF = DLAMCH( 'Underflow' )
      OVRTHR = DLAMCH( 'Overflow' )
      OVRTH2 = SQRT( OVRTHR )
C
      IF ( DELTA.EQ.ZERO ) THEN
C
C        The DELTA = 0 case.
C
         CALL DLASET( 'Full', N, N, ZERO, ONE, A, LDA )
         MDIG = NDECM1
         IDIG = NDECM1
         RETURN
      END IF
C
      IF ( N.EQ.1 ) THEN
C
C        The 1-by-1 case.
C
         A(1,1) = EXP( A(1,1)*DELTA )
         MDIG = NDECM1
         IDIG = NDECM1
         RETURN
      END IF
C
C     Set pointers for the workspace.
C
      JWORA1 = 1
      JWORA2 = JWORA1 + N*N
      JWORA3 = JWORA2 + N*NDIAG
      JWORV1 = JWORA3 + N*N
      JWORV2 = JWORV1 + N
C
C     Compute Pade coefficients in DWORK(JWORV2:JWORV2+NDIAG-1).
C
      DWORK(JWORV2) = HALF
C
      DO 20 I = 2, NDIAG
         IM1 = I - 1
         DWORK(JWORV2+IM1) = DWORK(JWORV2+I-2)*DBLE( NDIAG - IM1 )/
     $                                   DBLE( I*( 2*NDIAG - IM1 ) )
   20 CONTINUE
C
      VAREPS = EPS**2*( ( DBLE( BASE )**2 - ONE )/
     $        ( TWO4*LOG( DBLE( BASE ) ) ) )
      XN = DBLE( N )
      TR = ZERO
C
C     Apply a translation with the mean of the eigenvalues of A*DELTA.
C
      DO 40 I = 1, N
         CALL DSCAL( N, DELTA, A(1,I), 1 )
         TR = TR + A(I,I)
   40 CONTINUE
C
      AVGEV = TR/XN
      IF ( AVGEV.GT.LOG( OVRTHR ) .OR. AVGEV.LT.LOG( UNDERF ) )
     $   AVGEV = ZERO
      IF ( AVGEV.NE.ZERO ) THEN
         ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
C
         DO 60 I = 1, N
            A(I,I) = A(I,I) - AVGEV
   60    CONTINUE
C
         TEMP = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
         IF ( TEMP.GT.HALF*ANORM ) THEN
C
            DO 80 I = 1, N
               A(I,I) = A(I,I) + AVGEV
   80       CONTINUE
C
            AVGEV = ZERO
         END IF
      END IF
      ACTBAL = BALANC
      IF ( LBALS ) THEN
C
C        Balancing (scaling) has been requested.  First, save A.
C
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(JWORA1), N )
         MAXRED = TWOHND
         CALL MB04MD( N, MAXRED, A, LDA, DWORK(JWORV1), INFO )
         IF ( MAXRED.LT.ONE ) THEN
C
C           Recover the matrix and reset DWORK(JWORV1,...,JWORV1+N-1)
C           to 1, as no reduction of the norm occured (unlikely event).
C
            CALL DLACPY( 'Full', N, N, DWORK(JWORA1), N, A, LDA )
            ACTBAL = 'N'
            DWORK(JWORV1) = ONE
            CALL DCOPY( N-1, DWORK(JWORV1), 0, DWORK(JWORV1+1), 1 )
            IWARN = 3
         END IF
      END IF
C
C     Scale the matrix by 2**(-M), where M is the minimum integer
C     so that the resulted matrix has the 1-norm less than 0.5.
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
      M = 0
      IF ( ANORM.GE.HALF ) THEN
         MPOWER = INT( LOG( OVRTHR )/LOG( TWO ) )
         M = INT( LOG( ANORM )/LOG( TWO ) ) + 1
         IF ( M.GT.MPOWER ) THEN
C
C           Error return: The norm of A*DELTA is too large.
C
            INFO = 1
            RETURN
         END IF
         FACTOR = TWO**M
         IF ( M+1.LT.MPOWER ) THEN
            M = M + 1
            FACTOR = FACTOR*TWO
         END IF
C
         DO 120 I = 1, N
            CALL DSCAL( N, ONE/FACTOR, A(1,I), 1 )
  120    CONTINUE
C
      END IF
      NDAGM1 = NDIAG - 1
      NDAGM2 = NDAGM1 - 1
      IJ = 0
C
C     Compute the factors of the diagonal Pade approximant.
C     The loop 200 takes the accuracy requirements into account:
C     Pade coefficients decrease with K, so the calculations should
C     be performed in backward order, one column at a time.
C     (A BLAS 3 implementation in forward order, using DGEMM, could
C     possibly be less accurate.)
C
      DO 200 J = 1, N
         CALL DGEMV( 'No transpose', N, N, ONE, A, LDA, A(1,J), 1, ZERO,
     $               DWORK(JWORA2), 1 )
         IK = 0
C
         DO 140 K = 1, NDAGM2
            CALL DGEMV( 'No transpose', N, N, ONE, A, LDA,
     $                  DWORK(JWORA2+IK), 1, ZERO, DWORK(JWORA2+IK+N),
     $                  1 )
            IK = IK + N
  140    CONTINUE
C
         DO 180 I = 1, N
            S = ZERO
            U = ZERO
            IK = NDAGM2*N + I - 1
C
            DO 160 K = NDAGM1, 1, -1
               P = DWORK(JWORV2+K)*DWORK(JWORA2+IK)
               IK = IK - N
               S = S + P
               IF ( MOD( K+1, 2 ).EQ.0 ) THEN
                  U = U + P
               ELSE
                  U = U - P
               END IF
  160       CONTINUE
C
            P = DWORK(JWORV2)*A(I,J)
            S = S + P
            U = U - P
            IF ( I.EQ.J ) THEN
               S = S + ONE
               U = U + ONE
            END IF
            DWORK(JWORA3+IJ) = S
            DWORK(JWORA1+IJ) = U
            IJ = IJ + 1
  180    CONTINUE
C
  200 CONTINUE
C
C     Compute the exponential of the scaled matrix, using diagonal Pade
C     approximants.  As, in theory [1], the denominator of the Pade
C     approximant should be very well conditioned, no condition estimate
C     is computed.
C
      CALL DGETRF( N, N, DWORK(JWORA1), N, IWORK, IFAIL )
      IF ( IFAIL.GT.0 ) THEN
C
C        Error return: The matrix is exactly singular.
C
         INFO = 2
         RETURN
      END IF
C
      CALL DLACPY( 'Full', N, N, DWORK(JWORA3), N, A, LDA )
      CALL DGETRS( 'No transpose', N, N, DWORK(JWORA1), N, IWORK, A,
     $             LDA, IFAIL )
C
C     Prepare for the calculation of the accuracy estimates.
C     Note that ANORM here is in the range [1, e].
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
      IF ( ANORM.GE.ONE ) THEN
         EABS = ( NINTEN*XN + FOUR7 )*( EPS*ANORM )
      ELSE
         EABS = ( ( NINTEN*XN + FOUR7 )*EPS )*ANORM
      END IF
      IF ( M.NE.0 ) THEN
         VAR = XN*VAREPS
         FN = ( FOUR*XN )/( ( XN + TWO )*( XN + ONE ) )
         GN = ( ( TWO*XN + TEN )*XN - FOUR )/( ( ( XN + TWO )**2 )
     $                                        *( ( XN + ONE )**2 ) )
C
C        Square-up the computed exponential matrix M times, with caution
C        for avoiding overflows.
C
         DO 220 K = 1, M
            IF ( ANORM.GT.OVRTH2 ) THEN
C
C              The solution could overflow.
C
               CALL DGEMM( 'No transpose', 'No transpose', N, N, N,
     $                     ONE/ANORM, A, LDA, A, LDA, ZERO,
     $                     DWORK(JWORA1), N )
               S = DLANGE( '1-norm', N, N, DWORK(JWORA1), N,
     $                     DWORK(JWORA1) )
               IF ( ANORM.LE.OVRTHR/S ) THEN
                  CALL DLASCL( 'General', N, N, ONE, ANORM, N, N,
     $                         DWORK(JWORA1), N, INFO )
                  TEMP = OVRTHR
               ELSE
C
C                 Error return: The solution would overflow.
C                 This will not happen on most machines, due to the
C                 selection of M.
C
                  INFO = 3
                  RETURN
               END IF
            ELSE
               CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $                     A, LDA, A, LDA, ZERO, DWORK(JWORA1), N )
               TEMP = ANORM**2
            END IF
            IF ( EABS.LT.ONE ) THEN
               EABS = ( TWO*ANORM + EABS )*EABS  + XN*( EPS*TEMP )
            ELSE IF ( EABS.LT.SQRT( ONE - XN*EPS + OVRTHR/TEMP )*ANORM -
     $                        ANORM ) THEN
               EABS = XN*( EPS*TEMP ) + TWO*( ANORM*EABS ) + EABS**2
            ELSE
               EABS = OVRTHR
            END IF
C
            TMP1 = FN*VAR + GN*( TEMP*VAREPS )
            IF ( TMP1.GT.OVRTHR/TEMP ) THEN
               VAR = OVRTHR
            ELSE
               VAR = TMP1*TEMP
            END IF
C
            CALL DLACPY( 'Full', N, N, DWORK(JWORA1), N, A, LDA )
            ANORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
  220    CONTINUE
C
      ELSE
         VAR = ( TWELVE*XN )*VAREPS
      END IF
C
C     Apply back transformations, if balancing was effectively used.
C
      CALL MB05OY( ACTBAL, N, 1, N, A, LDA, DWORK(JWORV1), INFO )
      EAVGEV = EXP( AVGEV )
      EMNORM = DLANGE( '1-norm', N, N, A, LDA, DWORK(JWORA1) )
C
C     Compute auxiliary quantities needed for the accuracy estimates.
C
      BIG = ONE
      SMALL = ONE
      IF ( LBALS ) THEN
C
C        Compute norms of the diagonal scaling matrix and its inverse.
C
         DO 240 I = 1, N
            U = DWORK(JWORV1+I-1)
            IF (   BIG.LT.U )   BIG = U
            IF ( SMALL.GT.U ) SMALL = U
  240    CONTINUE
C
         SUM2D = DNRM2( N, DWORK(JWORV1), 1 )
      ELSE
         SUM2D = SQRT( XN )
      END IF
C
C     Update the exponential for the initial translation, and update the
C     auxiliary quantities needed for the accuracy estimates.
C
      SD2 = SQRT( EIGHT*XN*VAREPS )*ANORM
      BD  = SQRT( VAR )
      SS  = MAX( BD, SD2 )
      BD  = MIN( BD, SD2 )
      SD2 = SS*SQRT( ONE + ( BD/SS )**2 )
      IF ( SD2.LE.ONE ) THEN
         SD2 = ( TWO/XN )*SUM2D*SD2
      ELSE IF ( SUM2D/XN.LT.OVRTHR/TWO/SD2 ) THEN
         SD2 = ( TWO/XN )*SUM2D*SD2
      ELSE
         SD2 = OVRTHR
      END IF
      IF ( LBALS ) THEN
         SIZE = ZERO
      ELSE
         IF ( SD2.LT.OVRTHR - EMNORM ) THEN
            SIZE = EMNORM + SD2
         ELSE
            SIZE = OVRTHR
         END IF
      END IF
C
      DO 260 J = 1, N
         SS = DASUM( N, A(1,J), 1 )
         CALL DSCAL( N, EAVGEV, A(1,J), 1 )
         IF ( LBALS ) THEN
            BD = DWORK(JWORV1+J-1)
            SIZE = MAX( SIZE, SS + SD2/BD )
         END IF
  260 CONTINUE
C
C     Set the accuracy estimates and warning errors, if any.
C
      RERR = LOG10( BIG ) + LOG10( EABS ) - LOG10( SMALL ) -
     $       LOG10( EMNORM ) - LOG10( EPS )
      IF ( SIZE.GT.EMNORM ) THEN
         RERL = LOG10( ( SIZE/EMNORM - ONE )/EPS )
      ELSE
         RERL = ZERO
      END IF
      MDIG = MIN( NDEC - INT( RERR + HALF ), NDECM1 )
      IDIG = MIN( NDEC - INT( RERL + HALF ), NDECM1 )
C
      IF ( MDIG.LE.0 ) THEN
         MDIG = 0
         IWARN = 1
      END IF
      IF ( IDIG.LE.0 ) THEN
         IDIG = 0
         IWARN = 2
      END IF
C
      RETURN
C *** Last line of MB05OD ***
      END
