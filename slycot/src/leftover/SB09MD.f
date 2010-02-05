      SUBROUTINE SB09MD( N, NC, NB, H1, LDH1, H2, LDH2, SS, LDSS, SE,
     $                   LDSE, PRE, LDPRE, TOL, INFO )
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
C     To compare two multivariable sequences M1(k) and M2(k) for
C     k = 1,2,...,N, and evaluate their closeness. Each of the
C     parameters M1(k) and M2(k) is an NC by NB matrix.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of parameters.  N >= 0.
C
C     NC      (input) INTEGER
C             The number of rows in M1(k) and M2(k).  NC >= 0.
C
C     NB      (input) INTEGER
C             The number of columns in M1(k) and M2(k).  NB >= 0.
C
C     H1      (input) DOUBLE PRECISION array, dimension (LDH1,N*NB)
C             The leading NC-by-N*NB part of this array must contain
C             the multivariable sequence M1(k), where k = 1,2,...,N.
C             Each parameter M1(k) is an NC-by-NB matrix, whose
C             (i,j)-th element must be stored in H1(i,(k-1)*NB+j) for
C             i = 1,2,...,NC and j = 1,2,...,NB.
C
C     LDH1    INTEGER
C             The leading dimension of array H1.  LDH1 >= MAX(1,NC).
C
C     H2      (input) DOUBLE PRECISION array, dimension (LDH2,N*NB)
C             The leading NC-by-N*NB part of this array must contain
C             the multivariable sequence M2(k), where k = 1,2,...,N.
C             Each parameter M2(k) is an NC-by-NB matrix, whose
C             (i,j)-th element must be stored in H2(i,(k-1)*NB+j) for
C             i = 1,2,...,NC and j = 1,2,...,NB.
C
C     LDH2    INTEGER
C             The leading dimension of array H2.  LDH2 >= MAX(1,NC).
C
C     SS      (output) DOUBLE PRECISION array, dimension (LDSS,NB)
C             The leading NC-by-NB part of this array contains the
C             matrix SS.
C
C     LDSS    INTEGER
C             The leading dimension of array SS.  LDSS >= MAX(1,NC).
C
C     SE      (output) DOUBLE PRECISION array, dimension (LDSE,NB)
C             The leading NC-by-NB part of this array contains the
C             quadratic error matrix SE.
C
C     LDSE    INTEGER
C             The leading dimension of array SE.  LDSE >= MAX(1,NC).
C
C     PRE     (output) DOUBLE PRECISION array, dimension (LDPRE,NB)
C             The leading NC-by-NB part of this array contains the
C             percentage relative error matrix PRE.
C
C     LDPRE   INTEGER
C             The leading dimension of array PRE.  LDPRE >= MAX(1,NC).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in the computation of the error
C             matrices SE and PRE. If the user sets TOL to be less than
C             EPS then the tolerance is taken as EPS, where EPS is the
C             machine precision (see LAPACK Library routine DLAMCH).
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
C     The (i,j)-th element of the matrix SS is defined by:
C                        N          2
C               SS    = SUM  M1  (k) .                            (1)
C                 ij    k=1    ij
C
C     The (i,j)-th element of the quadratic error matrix SE is defined
C     by:
C                        N                      2
C               SE    = SUM  (M1  (k) - M2  (k)) .                (2)
C                 ij    k=1     ij        ij
C
C     The (i,j)-th element of the percentage relative error matrix PRE
C     is defined by:
C
C               PRE   = 100 x SQRT( SE  / SS  ).                  (3)
C                  ij                 ij    ij
C
C     The following precautions are taken by the routine to guard
C     against underflow and overflow:
C
C     (i) if ABS( M1  (k) ) > 1/TOL or ABS( M1  (k) - M2  (k) ) > 1/TOL,
C                   ij                        ij        ij
C
C         then SE   and SS   are set to 1/TOL and PRE   is set to 1; and
C                ij       ij                         ij
C
C     (ii) if ABS( SS  ) <= TOL, then PRE   is set to 100.
C                    ij                  ij
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately
C        2xNBxNCx(N+1) multiplications/divisions,
C        4xNBxNCxN     additions/subtractions and
C          NBxNC       square roots.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB09AD by S. Van Huffel, Katholieke
C     University Leuven, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Closeness multivariable sequences, elementary matrix operations,
C     real signals, system response.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HUNDRD
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, HUNDRD = 100.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH1, LDH2, LDPRE, LDSE, LDSS, N, NB, NC
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      DOUBLE PRECISION  H1(LDH1,*), H2(LDH2,*), PRE(LDPRE,*),
     $                  SE(LDSE,*), SS(LDSS,*)
C     .. Local Scalars ..
      LOGICAL           NOFLOW
      INTEGER           I, J, K
      DOUBLE PRECISION  EPSO, SSE, SSS, TOLER, VAR, VARE
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NC.LT.0 ) THEN
         INFO = -2
      ELSE IF( NB.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDH1.LT.MAX( 1, NC ) ) THEN
         INFO = -5
      ELSE IF( LDH2.LT.MAX( 1, NC ) ) THEN
         INFO = -7
      ELSE IF( LDSS.LT.MAX( 1, NC ) ) THEN
         INFO = -9
      ELSE IF( LDSE.LT.MAX( 1, NC ) ) THEN
         INFO = -11
      ELSE IF( LDPRE.LT.MAX( 1, NC ) ) THEN
         INFO = -13
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB09MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. NC.EQ.0 .OR. NB.EQ.0 )
     $   RETURN
C
      TOLER = MAX( TOL, DLAMCH( 'Epsilon' ) )
      EPSO = ONE/TOLER
C
      DO 60 J = 1, NB
C
         DO 40 I = 1, NC
            SSE = ZERO
            SSS = ZERO
            NOFLOW = .TRUE.
            K  = 0
C
C           WHILE ( ( NOFLOW .AND. ( K .LT. N*NB ) ) DO
   20       IF ( ( NOFLOW ) .AND. ( K.LT.N*NB ) ) THEN
               VAR  = H1(I,K+J)
               VARE = H2(I,K+J) - VAR
               IF ( ABS( VAR ).GT.EPSO .OR. ABS( VARE ).GT.EPSO )
     $            THEN
                  SE(I,J) = EPSO
                  SS(I,J) = EPSO
                  PRE(I,J) = ONE
                  NOFLOW = .FALSE.
               ELSE
                  IF ( ABS( VARE ).GT.TOLER ) SSE = SSE + VARE*VARE
                  IF ( ABS( VAR  ).GT.TOLER ) SSS = SSS + VAR*VAR
                  K  = K + NB
               END IF
               GO TO 20
            END IF
C           END WHILE 20
C
            IF ( NOFLOW ) THEN
               SE(I,J) = SSE
               SS(I,J) = SSS
               PRE(I,J) = HUNDRD
               IF ( SSS.GT.TOLER ) PRE(I,J) = SQRT( SSE/SSS )*HUNDRD
            END IF
   40    CONTINUE
C
   60 CONTINUE
C
      RETURN
C *** Last line of SB09MD ***
      END
