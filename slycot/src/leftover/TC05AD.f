      SUBROUTINE TC05AD( LERI, M, P, SVAL, INDEX, PCOEFF, LDPCO1,
     $                   LDPCO2, QCOEFF, LDQCO1, LDQCO2, RCOND, CFREQR,
     $                   LDCFRE, IWORK, DWORK, ZWORK, INFO )
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
C     To evaluate the transfer matrix T(s) of a left polynomial matrix
C     representation [T(s) = inv(P(s))*Q(s)] or a right polynomial
C     matrix representation [T(s) = Q(s)*inv(P(s))] at any specified
C     complex frequency s = SVAL.
C
C     This routine will calculate the standard frequency response
C     matrix at frequency omega if SVAL is supplied as (0.0,omega).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     LERI    CHARACTER*1
C             Indicates whether a left polynomial matrix representation
C             or a right polynomial matrix representation is to be used
C             to evaluate the transfer matrix as follows:
C             = 'L':  A left matrix fraction is input;
C             = 'R':  A right matrix fraction is input.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     SVAL    (input) COMPLEX*16
C             The frequency at which the transfer matrix or the
C             frequency respose matrix is to be evaluated.
C             For a standard frequency response set the real part
C             of SVAL to zero.
C
C     INDEX   (input) INTEGER array, dimension (MAX(M,P))
C             If LERI = 'L', INDEX(I), I = 1,2,...,P, must contain the
C             maximum degree of the polynomials in the I-th row of the
C             denominator matrix P(s) of the given left polynomial
C             matrix representation.
C             If LERI = 'R', INDEX(I), I = 1,2,...,M, must contain the
C             maximum degree of the polynomials in the I-th column of
C             the denominator matrix P(s) of the given right polynomial
C             matrix representation.
C
C     PCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDPCO1,LDPCO2,kpcoef), where kpcoef = MAX(INDEX(I)) + 1.
C             If LERI = 'L' then porm = P, otherwise porm = M.
C             The leading porm-by-porm-by-kpcoef part of this array must
C             contain the coefficients of the denominator matrix P(s).
C             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1)
C             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if
C             LERI = 'L' then iorj = I, otherwise iorj = J.
C             Thus for LERI = 'L', P(s) =
C             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...).
C             If LERI = 'R', PCOEFF is modified by the routine but
C             restored on exit.
C
C     LDPCO1  INTEGER
C             The leading dimension of array PCOEFF.
C             LDPCO1 >= MAX(1,P) if LERI = 'L',
C             LDPCO1 >= MAX(1,M) if LERI = 'R'.
C
C     LDPCO2  INTEGER
C             The second dimension of array PCOEFF.
C             LDPCO2 >= MAX(1,P) if LERI = 'L',
C             LDPCO2 >= MAX(1,M) if LERI = 'R'.
C
C     QCOEFF  (input) DOUBLE PRECISION array, dimension
C             (LDQCO1,LDQCO2,kpcoef)
C             If LERI = 'L' then porp = M, otherwise porp = P.
C             The leading porm-by-porp-by-kpcoef part of this array must
C             contain the coefficients of the numerator matrix Q(s).
C             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K).
C             If LERI = 'R', QCOEFF is modified by the routine but
C             restored on exit.
C
C     LDQCO1  INTEGER
C             The leading dimension of array QCOEFF.
C             LDQCO1 >= MAX(1,P)   if LERI = 'L',
C             LDQCO1 >= MAX(1,M,P) if LERI = 'R'.
C
C     LDQCO2  INTEGER
C             The second dimension of array QCOEFF.
C             LDQCO2 >= MAX(1,M)   if LERI = 'L',
C             LDQCO2 >= MAX(1,M,P) if LERI = 'R'.
C
C     RCOND   (output) DOUBLE PRECISION
C             The estimated reciprocal of the condition number of the
C             denominator matrix P(SVAL).
C             If RCOND is nearly zero, SVAL is approximately a system
C             pole.
C
C     CFREQR  (output) COMPLEX*16 array, dimension (LDCFRE,MAX(M,P))
C             The leading porm-by-porp part of this array contains the
C             frequency response matrix T(SVAL).
C
C     LDCFRE  INTEGER
C             The leading dimension of array CFREQR.
C             LDCFRE >= MAX(1,P)   if LERI = 'L',
C             LDCFRE >= MAX(1,M,P) if LERI = 'R'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (liwork)
C             where liwork = P, if LERI = 'L',
C                   liwork = M, if LERI = 'R'.
C
C     DWORK   DOUBLE PRECISION array, dimension (ldwork)
C             where ldwork = 2*P, if LERI = 'L',
C                   ldwork = 2*M, if LERI = 'R'.
C
C     ZWORK   COMPLEX*16 array, dimension (lzwork),
C             where lzwork = P*(P+2), if LERI = 'L',
C                   lzwork = M*(M+2), if LERI = 'R'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if P(SVAL) is exactly or nearly singular;
C                   no frequency response is calculated.
C
C     METHOD
C
C     The method for a left matrix fraction will be described here;
C     right matrix fractions are dealt with by obtaining the dual left
C     fraction and calculating its frequency response (see SLICOT
C     Library routine TC01OD). The first step is to calculate the
C     complex value P(SVAL) of the denominator matrix P(s) at the
C     desired frequency SVAL. If P(SVAL) is approximately singular,
C     SVAL is approximately a pole of this system and so the frequency
C     response matrix T(SVAL) is not calculated; in this case, the
C     routine returns with the Error Indicator (INFO) set to 1.
C     Otherwise, the complex value Q(SVAL) of the numerator matrix Q(s)
C     at frequency SVAL is calculated in a similar way to P(SVAL), and
C     the desired response matrix T(SVAL) = inv(P(SVAL))*Q(SVAL) is
C     found by solving the corresponding system of complex linear
C     equations.
C
C     REFERENCES
C
C     None
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TC01AD by T.W.C.Williams, Kingston
C     Polytechnic, United Kingdom, March 1982.
C
C     REVISIONS
C
C     February 22, 1998 (changed the name of TC01MD).
C
C     KEYWORDS
C
C     Coprime matrix fraction, elementary polynomial operations,
C     polynomial matrix, state-space representation, transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         LERI
      INTEGER           INFO, LDCFRE, LDPCO1, LDPCO2, LDQCO1, LDQCO2, M,
     $                  P
      DOUBLE PRECISION  RCOND
      COMPLEX*16        SVAL
C     .. Array Arguments ..
      INTEGER           INDEX(*), IWORK(*)
      DOUBLE PRECISION  DWORK(*), PCOEFF(LDPCO1,LDPCO2,*),
     $                  QCOEFF(LDQCO1,LDQCO2,*)
      COMPLEX*16        CFREQR(LDCFRE,*), ZWORK(*)
C     .. Local Scalars ..
      LOGICAL           LLERI
      INTEGER           I, IZWORK, IJ, INFO1, J, K, KPCOEF, LDZWOR,
     $                  MAXIND, MINMP, MPLIM, MWORK, PWORK
      DOUBLE PRECISION  CNORM
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, ZLANGE
      EXTERNAL          DLAMCH, LSAME, ZLANGE
C     .. External Subroutines ..
      EXTERNAL          TC01OD, XERBLA, ZCOPY, ZGECON, ZGETRF, ZGETRS,
     $                  ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LLERI = LSAME( LERI, 'L' )
      MPLIM = MAX( M, P )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LLERI .AND. .NOT.LSAME( LERI, 'R' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( LLERI .AND. LDPCO1.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO1.LT.MAX( 1, M ) ) ) THEN
         INFO = -7
      ELSE IF( ( LLERI .AND. LDPCO2.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDPCO2.LT.MAX( 1, M ) ) ) THEN
         INFO = -8
      ELSE IF( ( LLERI .AND. LDQCO1.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDQCO1.LT.MAX( 1, M, P ) ) ) THEN
         INFO = -10
      ELSE IF( ( LLERI .AND. LDQCO2.LT.MAX( 1, M ) ) .OR.
     $    ( .NOT.LLERI .AND. LDQCO2.LT.MAX( 1, MPLIM ) ) ) THEN
         INFO = -11
      ELSE IF( ( LLERI .AND. LDCFRE.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LLERI .AND. LDCFRE.LT.MAX( 1, MPLIM ) ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TC05AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( M.EQ.0 .OR. P.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
C
      IF ( LLERI ) THEN
C
C        Initialization for left matrix fraction.
C
         PWORK = P
         MWORK = M
      ELSE
C
C        Initialization for right matrix fraction: obtain dual system.
C
         PWORK = M
         MWORK = P
         IF ( MPLIM.GT.1 )
     $      CALL TC01OD( 'R', M, P, KPCOEF, PCOEFF, LDPCO1, LDPCO2,
     $                   QCOEFF, LDQCO1, LDQCO2, INFO )
      END IF
C
      LDZWOR = PWORK
      IZWORK = LDZWOR*LDZWOR + 1
      MAXIND = 0
C
      DO 10 I = 1, PWORK
         IF ( INDEX(I).GT.MAXIND ) MAXIND = INDEX(I)
   10 CONTINUE
C
      KPCOEF = MAXIND + 1
C
C     Calculate the complex denominator matrix P(SVAL), row by row.
C
      DO 50 I = 1, PWORK
         IJ = I
C
         DO 20 J = 1, PWORK
            ZWORK(IJ) = DCMPLX( PCOEFF(I,J,1), ZERO )
            IJ = IJ + PWORK
   20    CONTINUE
C
C        Possibly non-constant row: finish evaluating it.
C
         DO 40 K = 2, INDEX(I) + 1
C
            IJ = I
C
            DO 30 J = 1, PWORK
               ZWORK(IJ) = ( SVAL*ZWORK(IJ) ) +
     $                       DCMPLX( PCOEFF(I,J,K), ZERO )
               IJ = IJ + PWORK
   30       CONTINUE
C
   40    CONTINUE
C
   50 CONTINUE
C
C     Check if this P(SVAL) is singular: if so, don't compute T(SVAL).
C     Note that DWORK is not actually referenced in ZLANGE routine.
C
      CNORM = ZLANGE( '1-norm', PWORK, PWORK, ZWORK, LDZWOR, DWORK )
C
      CALL ZGETRF( PWORK, PWORK, ZWORK, LDZWOR, IWORK, INFO )
C
      IF ( INFO.GT.0 ) THEN
C
C        Singular matrix.  Set INFO and RCOND for error return.
C
         INFO  = 1
         RCOND = ZERO
      ELSE
C
C        Estimate the reciprocal condition of P(SVAL).
C        Workspace: ZWORK: PWORK*PWORK + 2*PWORK, DWORK: 2*PWORK.
C
         CALL ZGECON( '1-norm', PWORK, ZWORK, LDZWOR, CNORM, RCOND,
     $                ZWORK(IZWORK), DWORK, INFO )
C
         IF ( RCOND.LE.DLAMCH( 'Epsilon' ) ) THEN
C
C           Nearly singular matrix.  Set INFO for error return.
C
            INFO  = 1
         ELSE
C
C           Calculate the complex numerator matrix Q(SVAL), row by row.
C
            DO 90 I = 1, PWORK
C
               DO 60 J = 1, MWORK
                  CFREQR(I,J) = DCMPLX( QCOEFF(I,J,1), ZERO )
   60          CONTINUE
C
C              Possibly non-constant row: finish evaluating it.
C
               DO 80 K = 2, INDEX(I) + 1
C
                  DO 70 J = 1, MWORK
                     CFREQR(I,J) = ( SVAL*CFREQR(I,J) ) +
     $                             DCMPLX( QCOEFF(I,J,K), ZERO )
   70             CONTINUE
C
   80          CONTINUE
C
   90       CONTINUE
C
C           Now calculate frequency response T(SVAL).
C
            CALL ZGETRS( 'No transpose', PWORK, MWORK, ZWORK, LDZWOR,
     $                   IWORK, CFREQR, LDCFRE, INFO )
         END IF
      END IF
C
C     For right matrix fraction, return to original (dual of the dual)
C     system.
C
      IF ( ( .NOT.LLERI ) .AND. ( MPLIM.NE.1 ) ) THEN
         CALL TC01OD( 'L', MWORK, PWORK, KPCOEF, PCOEFF, LDPCO1,
     $                LDPCO2, QCOEFF, LDQCO1, LDQCO2, INFO1 )
C
         IF ( INFO.EQ.0 ) THEN
C
C           Also, transpose T(SVAL) here if this was successfully
C           calculated.
C
            MINMP = MIN( M, P )
C
            DO 100 J = 1, MPLIM
               IF ( J.LT.MINMP ) THEN
                  CALL ZSWAP( MINMP-J, CFREQR(J+1,J), 1, CFREQR(J,J+1),
     $                        LDCFRE )
               ELSE IF ( J.GT.P ) THEN
                  CALL ZCOPY( P, CFREQR(1,J), 1, CFREQR(J,1), LDCFRE )
               ELSE IF ( J.GT.M ) THEN
                  CALL ZCOPY( M, CFREQR(J,1), LDCFRE, CFREQR(1,J), 1 )
               END IF
  100       CONTINUE
C
         END IF
      END IF
C
      RETURN
C *** Last line of TC05AD ***
      END
