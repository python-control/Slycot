      SUBROUTINE MB05ND( N, DELTA, A, LDA, EX, LDEX, EXINT, LDEXIN,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To compute
C
C     (a)    F(delta) =  exp(A*delta) and
C
C     (b)    H(delta) =  Int[F(s) ds] from s = 0 to s = delta,
C
C     where A is a real N-by-N matrix and delta is a scalar value.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     DELTA   (input) DOUBLE PRECISION
C             The scalar value delta of the problem.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             matrix A of the problem. (Array A need not be set if
C             DELTA = 0.)
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,N).
C
C     EX      (output) DOUBLE PRECISION array, dimension (LDEX,N)
C             The leading N-by-N part of this array contains an
C             approximation to F(delta).
C
C     LDEX    INTEGER
C             The leading dimension of array EX.  LDEX >= MAX(1,N).
C
C     EXINT   (output) DOUBLE PRECISION array, dimension (LDEXIN,N)
C             The leading N-by-N part of this array contains an
C             approximation to H(delta).
C
C     LDEXIN  INTEGER
C             The leading dimension of array EXINT.  LDEXIN >= MAX(1,N).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the order of the
C             Pade approximation to H(t), where t is a scale factor
C             determined by the routine. A reasonable value for TOL may
C             be SQRT(EPS), where EPS is the machine precision (see
C             LAPACK Library routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK. LDWORK >= MAX(1,N*(N+1)).
C             For optimum performance LDWORK should be larger (2*N*N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, the (i,i) element of the denominator of
C                   the Pade approximation is zero, so the denominator
C                   is exactly singular;
C             = N+1:  if DELTA = (delta * frobenius norm of matrix A) is
C                   probably too large to permit meaningful computation.
C                   That is, DELTA > SQRT(BIG), where BIG is a
C                   representable number near the overflow threshold of
C                   the machine (see LAPACK Library Routine DLAMCH).
C
C     METHOD
C
C     This routine uses a Pade approximation to H(t) for some small
C     value of t (where 0 < t <= delta) and then calculates F(t) from
C     H(t). Finally, the results are re-scaled to give F(delta) and
C     H(delta). For a detailed description of the implementation of this
C     algorithm see [1].
C
C     REFERENCES
C
C     [1] Benson, C.J.
C         The numerical evaluation of the matrix exponential and its
C         integral.
C         Report 82/03, Control Systems Research Group,
C         School of Electronic Engineering and Computer
C         Science, Kingston Polytechnic, January 1982.
C
C     [2] Ward, R.C.
C         Numerical computation of the matrix exponential with accuracy
C         estimate.
C         SIAM J. Numer. Anal., 14, pp. 600-610, 1977.
C
C     [3] Moler, C.B. and Van Loan, C.F.
C         Nineteen Dubious Ways to Compute the Exponential of a Matrix.
C         SIAM Rev., 20, pp. 801-836, 1978.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997.
C     Supersedes Release 2.0 routine MB05BD by C.J. Benson, Kingston
C     Polytechnic, January 1982.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Continuous-time system, matrix algebra, matrix exponential,
C     matrix operations, Pade approximation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, ONE64, THREE, FOUR8
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                    ONE64 = 1.64D0, THREE = 3.0D0, FOUR8 = 4.8D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDEX, LDEXIN, LDWORK, N
      DOUBLE PRECISION  DELTA, TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), EX(LDEX,*), EXINT(LDEXIN,*)
C     .. Local Scalars ..
      INTEGER           I, I2IQ1, IJ, IQ, J, JSCAL, KK, L, NN
      DOUBLE PRECISION  COEFFD, COEFFN, DELSC, EPS, ERR, F2IQ1,
     $                  FNORM, FNORM2, QMAX, SMALL
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DGEMV, DGESV, DLACPY,
     $                  DLASET, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, EXP, MAX, MOD, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      NN = N*N
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDEX.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDEXIN.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDWORK.LT.MAX( 1, NN + N ) ) THEN
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB05ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      DWORK(1) = ONE
      IF ( N.EQ.0 )
     $   RETURN
C
      CALL DLASET( 'Full', N, N, ZERO, ZERO, EX, LDEX )
      CALL DLASET( 'Full', N, N, ZERO, ZERO, EXINT, LDEXIN )
C
      IF ( DELTA.EQ.ZERO ) THEN
         CALL DLASET( 'Upper', N, N, ZERO, ONE, EX, LDEX )
         RETURN
      END IF
C
      IF ( N.EQ.1 ) THEN
         EX(1,1) = EXP( DELTA*A(1,1) )
         IF ( A(1,1).EQ.ZERO ) THEN
            EXINT(1,1) = DELTA
         ELSE
            EXINT(1,1) = ( ( ONE/A(1,1) )*EX(1,1) ) - ( ONE/A(1,1) )
         END IF
         RETURN
      END IF
C
C     Set some machine parameters.
C
      EPS = DLAMCH( 'Epsilon' )
      SMALL = DLAMCH( 'Safe minimum' )/EPS
C
C     First calculate the Frobenius norm of A, and the scaling factor.
C
      FNORM = DELTA*DLANGE( 'Frobenius', N, N, A, LDA, DWORK )
C
      IF ( FNORM.GT.SQRT( ONE/SMALL ) ) THEN
         INFO = N + 1
         RETURN
      END IF
C
      JSCAL = 0
      DELSC = DELTA
C     WHILE ( FNORM >= HALF ) DO
   20 CONTINUE
      IF ( FNORM.GE.HALF ) THEN
         JSCAL = JSCAL + 1
         DELSC = DELSC*HALF
         FNORM = FNORM*HALF
         GO TO 20
      END IF
C     END WHILE 20
C
C     Calculate the order of the Pade approximation needed to satisfy
C     the requested relative error  TOL.
C
      FNORM2 = FNORM**2
      IQ = 1
      QMAX = FNORM/THREE
      ERR  = DELTA/DELSC*FNORM2**2/FOUR8
C     WHILE ( ERR > TOL*( 2*IQ + 3 - FNORM )/1.64 and QMAX >= EPS ) DO
   40 CONTINUE
      IF ( ERR.GT.TOL*( DBLE( 2*IQ + 3 ) - FNORM )/ONE64 ) THEN
         IQ = IQ + 1
         QMAX = QMAX*DBLE( IQ + 1 )*FNORM/DBLE( 2*IQ*( 2*IQ + 1 ) )
         IF ( QMAX.GE.EPS ) THEN
            ERR = ERR*FNORM2*DBLE( 2*IQ + 5 )/DBLE( ( 2*IQ + 3 )**2
     $                          *( 2*IQ + 4 ) )
            GO TO 40
         END IF
      END IF
C     END WHILE 40
C
C     Initialise DWORK (to contain succesive powers of A),
C                EXINT (to contain the numerator) and
C                EX    (to contain the denominator).
C
      I2IQ1 = 2*IQ + 1
      F2IQ1 = DBLE( I2IQ1 )
      COEFFD = -DBLE( IQ )/F2IQ1
      COEFFN = HALF/F2IQ1
      IJ = 1
C
      DO 80 J = 1, N
C
         DO 60 I = 1, N
            DWORK(IJ)  = DELSC*A(I,J)
            EXINT(I,J) = COEFFN*DWORK(IJ)
            EX(I,J)    = COEFFD*DWORK(IJ)
            IJ = IJ + 1
   60    CONTINUE
C
         EXINT(J,J) = EXINT(J,J) + ONE
         EX(J,J) = EX(J,J) + ONE
   80 CONTINUE
C
      DO 140 KK = 2, IQ
C
C        Calculate the next power of  A*DELSC,  and update the numerator
C        and denominator.
C
         COEFFD = -COEFFD*DBLE( IQ+1-KK )/DBLE( KK*( I2IQ1+1-KK ) )
         IF ( MOD( KK, 2 ).EQ.0 ) THEN
            COEFFN = COEFFD/DBLE( KK + 1 )
         ELSE
            COEFFN = -COEFFD/DBLE( I2IQ1 - KK )
         END IF
         IJ = 1
C
         IF ( LDWORK.GE.2*NN ) THEN
C
C           Enough space for a BLAS 3 calculation.
C
            CALL DGEMM( 'No transpose', 'No transpose', N, N, N, DELSC,
     $                  A, LDA, DWORK, N, ZERO, DWORK(NN+1), N )
            CALL DCOPY( NN, DWORK(NN+1), 1, DWORK, 1 )
C
            DO 100 J = 1, N
               CALL DAXPY( N, COEFFN, DWORK(IJ), 1, EXINT(1,J), 1 )
               CALL DAXPY( N, COEFFD, DWORK(IJ), 1, EX(1,J), 1 )
               IJ = IJ + N
  100       CONTINUE
C
         ELSE
C
C           Not enough space for a BLAS 3 calculation. Use BLAS 2.
C
            DO 120 J = 1, N
               CALL DGEMV( 'No transpose', N, N, ONE, A, LDA, DWORK(IJ),
     $                     1, ZERO, DWORK(NN+1), 1 )
               CALL DCOPY( N, DWORK(NN+1), 1, DWORK(IJ), 1 )
               CALL DSCAL( N, DELSC, DWORK(IJ), 1 )
               CALL DAXPY( N, COEFFN, DWORK(IJ), 1, EXINT(1,J), 1 )
               CALL DAXPY( N, COEFFD, DWORK(IJ), 1, EX(1,J), 1 )
               IJ = IJ + N
  120       CONTINUE
C
         END IF
  140 CONTINUE
C
C     We now have numerator in EXINT, denominator in EX.
C
C     Solve the set of N systems of linear equations for the columns of
C     EXINT  using the LU factorization of EX.
C
      CALL DGESV( N, N, EX, LDEX, IWORK, EXINT, LDEXIN, INFO )
      IF ( INFO.NE.0 )
     $   RETURN
C
C     Now we can form EX from EXINT using the formula:
C     EX = EXINT * A  +  I
C
      DO 160 J = 1, N
         CALL DSCAL( N, DELSC, EXINT(1,J), 1 )
  160 CONTINUE
C
      CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE, EXINT,
     $            LDEXIN, A, LDA, ZERO, EX, LDEX  )
C
      DO 180 J = 1, N
         EX(J,J) = EX(J,J) + ONE
  180 CONTINUE
C
C     EX  and  EXINT  have been evaluated at  DELSC,  so the results
C     must be re-scaled to give the function values at  DELTA.
C
C     EXINT(2t) = EXINT(t) * ^ EX(t) + I [
C     EX(2t) = EX(t) * EX(t)
C
C     DWORK  is used to accumulate products.
C
      DO 200 L = 1, JSCAL
         CALL DLACPY( 'Full', N, N, EXINT, LDEXIN, DWORK, N )
         CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $               DWORK, N, EX, LDEX, ONE, EXINT, LDEXIN )
         CALL DLACPY( 'Full', N, N, EX, LDEX, DWORK, N )
         CALL DGEMM( 'No transpose', 'No transpose', N, N, N, ONE,
     $               DWORK, N, DWORK, N, ZERO, EX, LDEX )
  200 CONTINUE
C
      DWORK(1) = 2*NN
      RETURN
C *** Last line of MB05ND ***
      END
