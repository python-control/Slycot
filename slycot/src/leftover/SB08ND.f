      SUBROUTINE SB08ND( ACONA, DA, A, RES, E, DWORK, LDWORK, INFO )
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
C     To compute a real polynomial E(z) such that
C
C        (a)  E(1/z) * E(z) = A(1/z) * A(z) and
C        (b)  E(z) is stable - that is, E(z) has no zeros with modulus
C             greater than 1,
C
C     which corresponds to computing the spectral factorization of the
C     real polynomial A(z) arising from discrete optimality problems.
C
C     The input polynomial may be supplied either in the form
C
C     A(z) = a(0) + a(1) * z + ... + a(DA) * z**DA
C
C     or as
C
C     B(z) = A(1/z) * A(z)
C          = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA)
C                                                                    (1)
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ACONA   CHARACTER*1
C             Indicates whether the coefficients of A(z) or B(z) =
C             A(1/z) * A(z) are to be supplied as follows:
C             = 'A':  The coefficients of A(z) are to be supplied;
C             = 'B':  The coefficients of B(z) are to be supplied.
C
C     Input/Output Parameters
C
C     DA      (input) INTEGER
C             The degree of the polynomials A(z) and E(z).  DA >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (DA+1)
C             On entry, if ACONA = 'A', this array must contain the
C             coefficients of the polynomial A(z) in increasing powers
C             of z, and if ACONA = 'B', this array must contain the
C             coefficients b ,b ,...,b   of the polynomial B(z) in
C                           0  1      DA
C             equation (1). That is, A(i) = b    for i = 1,2,...,DA+1.
C                                            i-1
C             On exit, this array contains the coefficients of the
C             polynomial B(z) in eqation (1). Specifically, A(i)
C             contains b   ,  for i = 1,2,...DA+1.
C                       i-1
C
C     RES     (output) DOUBLE PRECISION
C             An estimate of the accuracy with which the coefficients of
C             the polynomial E(z) have been computed (see also METHOD
C             and NUMERICAL ASPECTS).
C
C     E       (output) DOUBLE PRECISION array, dimension (DA+1)
C             The coefficients of the spectral factor E(z) in increasing
C             powers of z.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= 5*DA+5.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 2:  if on entry, ACONA = 'B' but the supplied
C                   coefficients of the polynomial B(z) are not the
C                   coefficients of A(1/z) * A(z) for some real A(z);
C                   in this case, RES and E are unassigned;
C             = 3:  if the iterative process (see METHOD) has failed to
C                   converge in 30 iterations;
C             = 4:  if the last computed iterate (see METHOD) is
C                   unstable. If ACONA = 'B', then the supplied
C                   coefficients of the polynomial B(z) may not be the
C                   coefficients of A(1/z) * A(z) for some real A(z).
C
C     METHOD
C         _                                               _
C     Let A(z) be the conjugate polynomial of A(z), i.e., A(z) = A(1/z).
C
C     The method used by the routine is based on applying the
C     Newton-Raphson iteration to the function
C               _       _
C        F(e) = A * A - e * e,
C
C     which leads to the iteration formulae (see [1] and [2])
C
C        _(i)   (i)  _(i)   (i)     _      )
C        q   * x   + x   * q    = 2 A * A  )
C                                          )   for i = 0, 1, 2,...
C         (i+1)    (i)   (i)               )
C        q     = (q   + x   )/2            )
C
C     The iteration starts from
C
C         (0)                                        DA
C        q   (z) = (b(0) + b(1) * z + ... + b(DA) * z  ) / SQRT( b(0))
C
C     which is a Hurwitz polynomial that has no zeros in the closed unit
C                                            (i)
C     circle (see [2], Theorem 3). Then lim q   = e, the convergence is
C     uniform and e is a Hurwitz polynomial.
C
C     The iterates satisfy the following conditions:
C              (i)
C        (a)  q    has no zeros in the closed unit circle,
C              (i)     (i-1)
C        (b)  q    <= q     and
C              0       0
C              DA   (i) 2    DA     2
C        (c)  SUM (q   )  - SUM (A )  >= 0.
C             k=0   k       k=0   k
C                                     (i)
C     The iterative process stops if q    violates (a), (b) or (c),
C     or if the condition
C                       _(i) (i)  _
C        (d)  RES  = ||(q   q   - A A)|| < tol,
C
C     is satisfied, where || . || denotes the largest coefficient of
C                     _(i) (i)  _
C     the polynomial (q   q   - A A) and tol is an estimate of the
C                                                    _(i)  (i)
C     rounding error in the computed coefficients of q    q   . If
C                                            (i-1)
C     condition (a) or (b) is violated then q      is taken otherwise
C      (i)
C     q    is used. Thus the computed reciprocal polynomial E(z) = z**DA
C     * q(1/z) is stable. If there is no convergence after 30 iterations
C     then the routine returns with the Error Indicator (INFO) set to 3,
C     and the value of RES may indicate whether or not the last computed
C     iterate is close to the solution.
C                                               (0)
C     If ACONA = 'B', then it is possible that q    is not a Hurwitz
C     polynomial, in which case the equation e(1/z) * e(z) = B(z) has no
C     real solution (see [2], Theorem 3).
C
C     REFERENCES
C
C     [1] Kucera, V.
C         Discrete Linear Control, The polynomial Approach.
C         John Wiley & Sons, Chichester, 1979.
C
C     [2] Vostry, Z.
C         New Algorithm for Polynomial Spectral Factorization with
C         Quadratic Convergence I.
C         Kybernetika, 11, pp. 415-422, 1975.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB08BD by F. Delebecque and
C     A.J. Geurts.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Factorization, Laplace transform, optimal control, optimal
C     filtering, polynomial operations, spectral factorization, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                    TWO  = 2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         ACONA
      INTEGER           DA, INFO, LDWORK
      DOUBLE PRECISION  RES
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), DWORK(*), E(*)
C     .. Local Scalars ..
      LOGICAL           CONV, HURWTZ, LACONA
      INTEGER           I, J, K, LALPHA, LAMBDA, LETA, LQ, LRO, NC, NCK
      DOUBLE PRECISION  A0, RES0, S, SA0, TOLQ, W
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, DSWAP, SB08NY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      LACONA = LSAME( ACONA, 'A' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LACONA .AND. .NOT.LSAME( ACONA, 'B' ) ) THEN
         INFO = -1
      ELSE IF( DA.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDWORK.LT.5*DA + 5 ) THEN
         INFO = -7
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB08ND', -INFO )
         RETURN
      END IF
C
      NC = DA + 1
      IF ( .NOT.LACONA ) THEN
         IF ( A(1).LE.ZERO ) THEN
            INFO = 2
            RETURN
         END IF
         CALL DCOPY( NC, A, 1, E, 1 )
      ELSE
         CALL SB08NY( DA, A, E, W )
      END IF
C
C     Initialization.
C
      LALPHA = 1
      LRO = LALPHA + NC
      LETA = LRO + NC
      LAMBDA = LETA + NC
      LQ = LAMBDA + NC
C
      A0 = E(1)
      SA0 = SQRT( A0 )
      S = ZERO
C
      DO 20 J = 1, NC
         W = E(J)
         A(J) = W
         W = W/SA0
         E(J) = W
         DWORK(LQ-1+J) = W
         S = S + W**2
   20 CONTINUE
C
      RES0 = S - A0
C
C     The contents of the arrays is, cf [1], Section 7.6,
C
C     E : the last computed Hurwitz polynomial q   ;
C                                               i-1
C     DWORK(LALPHA,..,LALPHA+DA-K)  : alpha(k,0),...alpha(k,n-k);
C          (LRO,...,LRO+DA-K)       : alpha(k,n-k),...,alpha(k);
C          (LETA,...,LETA+DA)       : eta(0),...,eta(n);
C          (LAMBDA,...,LAMBDA+DA-1) : lambda(0),...,lambda(n-1)
C
C     DWORK(LQ,...,LQ+DA) : the last computed polynomial q .
C                                                         i
      I = 0
      CONV = .FALSE.
      HURWTZ = .TRUE.
C
C     WHILE ( I < 30 and CONV = FALSE and HURWTZ = TRUE ) DO
   40 IF ( I.LT.30 .AND. .NOT.CONV .AND. HURWTZ ) THEN
         I = I + 1
         CALL DCOPY( NC, A, 1, DWORK(LETA), 1 )
         CALL DSCAL( NC, TWO, DWORK(LETA), 1 )
         CALL DCOPY( NC, DWORK(LQ), 1, DWORK(LALPHA), 1 )
C
C        Computation of lambda(k) and eta(k).
C
         K = 1
C
C        WHILE ( K <= DA and HURWTZ = TRUE ) DO
   60    IF ( ( K.LE.DA ) .AND. HURWTZ ) THEN
            NCK = NC - K
            CALL DCOPY( NCK+1, DWORK(LALPHA), -1, DWORK(LRO), 1 )
            W = DWORK(LALPHA+NCK)/DWORK(LRO+NCK)
            IF ( ABS( W ).GE.ONE ) HURWTZ = .FALSE.
            IF ( HURWTZ ) THEN
               DWORK(LAMBDA+K-1) = W
               CALL DAXPY( NCK, -W, DWORK(LRO), 1, DWORK(LALPHA), 1 )
               W = DWORK(LETA+NCK)/DWORK(LALPHA)
               DWORK(LETA+NCK) = W
               CALL DAXPY( NCK-1, -W, DWORK(LALPHA+1), -1,
     $                     DWORK(LETA+1), 1 )
               K = K + 1
            END IF
            GO TO 60
         END IF
C        END WHILE 60
C
C        HURWTZ = The polynomial q    is a Hurwitz polynomial.
C                                 i-1
         IF ( HURWTZ ) THEN
            CALL DCOPY( NC, DWORK(LQ), 1, E, 1 )
C
C           Accuracy test.
C
            CALL SB08NY( DA, E, DWORK(LQ), TOLQ )
            CALL DAXPY( NC, -ONE, A, 1, DWORK(LQ), 1 )
            RES = ABS( DWORK( IDAMAX( NC, DWORK(LQ), 1 ) + LQ - 1 ) )
            CONV = ( RES.LT.TOLQ ) .OR. ( RES0.LT.ZERO )
C
            IF ( .NOT.CONV ) THEN
               DWORK(LETA) = HALF*DWORK(LETA)/DWORK(LALPHA)
C
C              Computation of x  and q .
C                              i      i
C              DWORK(LETA,...,LETA+DA)   : eta(k,0),...,eta(k,n)
C                   (LRO,...,LRO+DA-K+1) : eta(k,n-k+1),...,eta(k,0)
C
               DO 80 K = DA, 1, -1
                  NCK = NC - K + 1
                  CALL DCOPY( NCK, DWORK(LETA), -1, DWORK(LRO), 1 )
                  W = DWORK(LAMBDA+K-1)
                  CALL DAXPY( NCK, -W, DWORK(LRO), 1, DWORK(LETA), 1 )
   80          CONTINUE
C
               S = ZERO
C
               DO 100 J = 0, DA
                  W = HALF*( DWORK(LETA+J) + E(J+1) )
                  DWORK(LQ+J) = W
                  S = S + W**2
  100          CONTINUE
C
               RES0 = S - A0
C
C              Test on the monotonicity of q .
C                                           0
               CONV = DWORK(LQ).GT.E(1)
               GO TO 40
            END IF
         END IF
      END IF
C     END WHILE 40
C
C     Reverse the order of the coefficients in the array E.
C
      CALL DSWAP( NC, E, 1, DWORK, -1 )
      CALL DSWAP( NC, DWORK, 1, E, 1 )
C
      IF ( .NOT.CONV ) THEN
         IF ( HURWTZ ) THEN
            INFO = 3
         ELSE IF ( I.EQ.1 ) THEN
            INFO = 2
         ELSE
            INFO = 4
         END IF
      END IF
C
      RETURN
C *** Last line of SB08ND ***
      END
