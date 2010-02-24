      SUBROUTINE TB01VD( APPLY, N, M, L, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   X0, THETA, LTHETA, DWORK, LDWORK, INFO )
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
C     To convert the linear discrete-time system given as (A, B, C, D),
C     with initial state x0, into the output normal form [1], with
C     parameter vector THETA. The matrix A is assumed to be stable.
C     The matrices A, B, C, D and the vector x0 are converted, so that
C     on exit they correspond to the system defined by THETA.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     APPLY   CHARACTER*1
C             Specifies whether or not the parameter vector should be
C             transformed using a bijective mapping, as follows:
C             = 'A' : apply the bijective mapping to the N vectors in
C                     THETA corresponding to the matrices A and C;
C             = 'N' : do not apply the bijective mapping.
C             The transformation performed when APPLY = 'A' allows
C             to get rid of the constraints norm(THETAi) < 1, i = 1:N.
C             A call of the SLICOT Library routine TB01VY associated to
C             a call of TB01VD must use the same value of APPLY.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     L       (input) INTEGER
C             The number of system outputs.  L >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the system state matrix A, assumed to be stable.
C             On exit, the leading N-by-N part of this array contains
C             the transformed system state matrix corresponding to the
C             output normal form with parameter vector THETA.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the system input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed system input matrix corresponding to the
C             output normal form with parameter vector THETA.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading L-by-N part of this array must
C             contain the system output matrix C.
C             On exit, the leading L-by-N part of this array contains
C             the transformed system output matrix corresponding to the
C             output normal form with parameter vector THETA.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,L).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading L-by-M part of this array must contain the
C             system input/output matrix D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,L).
C
C     X0      (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the initial state of the
C             system, x0.
C             On exit, this array contains the transformed initial state
C             of the system, corresponding to the output normal form
C             with parameter vector THETA.
C
C     THETA   (output) DOUBLE PRECISION array, dimension (LTHETA)
C             The leading N*(L+M+1)+L*M part of this array contains the
C             parameter vector that defines a system (A, B, C, D, x0)
C             which is equivalent up to a similarity transformation to
C             the system given on entry. The parameters are:
C
C             THETA(1:N*L)                      : parameters for A, C;
C             THETA(N*L+1:N*(L+M))              : parameters for B;
C             THETA(N*(L+M)+1:N*(L+M)+L*M)      : parameters for D;
C             THETA(N*(L+M)+L*M+1:N*(L+M+1)+L*M): parameters for x0.
C
C     LTHETA  INTEGER
C             The length of array THETA.  LTHETA >= N*(L+M+1)+L*M.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N*N*L + N*L + N,
C                           N*N + MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L),
C                                     N*M)).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the Lyapunov equation A'*Q*A - Q = -scale^2*C'*C
C                   could only be solved with scale = 0;
C             = 2:  if matrix A is not discrete-time stable;
C             = 3:  if the QR algorithm failed to converge for
C                   matrix A.
C
C     METHOD
C
C     The matrices A and C are converted to output normal form.
C     First, the Lyapunov equation
C
C        A'*Q*A - Q = -scale^2*C'*C,
C
C     is solved in the Cholesky factor T, T'*T = Q, and then T is used
C     to get the transformation matrix.
C
C     The matrix B and the initial state x0 are transformed accordingly.
C
C     Then, the QR factorization of the transposed observability matrix
C     is computed, and the matrix Q is used to further transform the
C     system matrices. The parameters characterizing A and C are finally
C     obtained by applying a set of N orthogonal transformations.
C
C     REFERENCES
C
C     [1] Peeters, R.L.M., Hanzon, B., and Olivi, M.
C         Balanced realizations of discrete-time stable all-pass
C         systems and the tangential Schur algorithm.
C         Proceedings of the European Control Conference,
C         31 August - 3 September 1999, Karlsruhe, Germany.
C         Session CP-6, Discrete-time Systems, 1999.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Oct. 2000, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Feb. 2002, Feb. 2004.
C
C     KEYWORDS
C
C     Asymptotically stable, Lyapunov equation, output normal form,
C     parameter estimation, similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, HALF
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                    HALF = 0.5D0 )
C     .. Scalar Arguments ..
      CHARACTER         APPLY
      INTEGER           INFO, L, LDA, LDB, LDC, LDD, LDWORK, LTHETA, M,
     $                  N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), THETA(*), X0(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  PIBY2, RI, SCALE, TI
      INTEGER           CA, I, IA, IN, IQ, IR, IT, ITAU, IU, IWI, IWR,
     $                  J, JWORK, K, LDCA, LDT, WRKOPT
      LOGICAL           LAPPLY
C     .. External Functions ..
      EXTERNAL          DNRM2, LSAME
      DOUBLE PRECISION  DNRM2
      LOGICAL           LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DGEMV, DGEQRF, DGER,
     $                  DLACPY, DLASET, DORMQR, DSCAL, DTRMM, DTRMV,
     $                  DTRSM, MA02AD, SB03OD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, INT, MAX, MIN, SQRT, TAN
C     ..
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      LAPPLY = LSAME( APPLY, 'A' )
C
      INFO = 0
      IF ( .NOT.( LAPPLY .OR. LSAME( APPLY, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( N.LT.0 ) THEN
         INFO = -2
      ELSEIF ( M.LT.0 ) THEN
         INFO = -3
      ELSEIF ( L.LT.0 ) THEN
         INFO = -4
      ELSEIF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSEIF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSEIF ( LDC.LT.MAX( 1, L ) ) THEN
         INFO = -10
      ELSEIF ( LDD.LT.MAX( 1, L ) ) THEN
         INFO = -12
      ELSEIF ( LTHETA.LT.( N*( M + L + 1 ) + L*M ) ) THEN
         INFO = -15
      ELSEIF ( LDWORK.LT.MAX( 1, N*N*L + N*L + N, N*N +
     $                        MAX( N*( N + MAX( N, L ) + 6 ) +
     $                             MIN( N, L ), N*M ) ) ) THEN
         INFO = -17
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TB01VD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, L ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      ELSE IF ( N.EQ.0 ) THEN
         CALL DLACPY( 'Full', L, M, D, LDD, THETA, MAX( 1, L ) )
         DWORK(1) = ONE
         RETURN
      ELSE IF ( L.EQ.0 ) THEN
         CALL DLACPY( 'Full', N, M, B, LDB, THETA, N )
         CALL DCOPY( N, X0, 1, THETA(N*M+1), 1 )
         DWORK(1) = ONE
         RETURN
      ENDIF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      WRKOPT = 1
      PIBY2  = TWO*ATAN( ONE )
C
C     Convert A and C to output normal form.
C     First, solve the Lyapunov equation
C        A'*Q*A - Q = -scale^2*C'*C,
C     in the Cholesky factor T, T'*T = Q, and use T to get the
C     transformation matrix. Copy A and C, to preserve them.
C
C     Workspace: need   N*(2*N + MAX(N,L) + 6) + MIN(N,L).
C                prefer larger.
C
C     Initialize the indices in the workspace.
C
      LDT = MAX( N, L )
      CA  = 1
      IA  = 1
      IT  = IA  + N*N
      IU  = IT  + LDT*N
      IWR = IU  + N*N
      IWI = IWR + N
C
      JWORK = IWI + N
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK(IA), N )
      CALL DLACPY( 'Full', L, N, C, LDC, DWORK(IT), LDT )
C
      CALL SB03OD( 'Discrete', 'NotFactored', 'NoTranspose', N, L,
     $             DWORK(IA), N, DWORK(IU), N, DWORK(IT), LDT, SCALE,
     $             DWORK(IWR), DWORK(IWI), DWORK(JWORK), LDWORK-JWORK+1,
     $             INFO )
      IF ( INFO.NE.0 ) THEN
         IF ( INFO.EQ.6 ) THEN
            INFO = 3
         ELSE
            INFO = 2
         ENDIF
         RETURN
      ENDIF
      WRKOPT = INT( DWORK(JWORK) ) + JWORK - 1
C
      IF ( SCALE.EQ.ZERO ) THEN
         INFO = 1
         RETURN
      ENDIF
C
C     Compute A = T*A*T^(-1).
C
      CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N, N, ONE,
     $            DWORK(IT), LDT, A, LDA )
C
      CALL DTRSM( 'Right', 'Upper', 'NoTranspose', 'NonUnit', N, N, ONE,
     $            DWORK(IT), LDT, A, LDA )
      IF ( M.GT.0 ) THEN
C
C        Compute B = (1/scale)*T*B.
C
         CALL DTRMM( 'Left', 'Upper', 'NoTranspose', 'NonUnit', N, M,
     $               ONE/SCALE, DWORK(IT), LDT, B, LDB )
      ENDIF
C
C     Compute x0 = (1/scale)*T*x0.
C
      CALL DTRMV( 'Upper', 'NoTranspose', 'NonUnit', N, DWORK(IT), LDT,
     $            X0, 1 )
      CALL DSCAL( N, ONE/SCALE, X0, 1 )
C
C     Compute C = scale*C*T^(-1).
C
      CALL DTRSM( 'Right', 'Upper', 'NoTranspose', 'NonUnit', L, N,
     $            SCALE, DWORK(IT), LDT, C, LDC )
C
C     Now, the system has been transformed to the output normal form.
C     Build the transposed observability matrix in DWORK(CA) and compute
C     its QR factorization.
C
      CALL MA02AD( 'Full', L, N, C, LDC, DWORK(CA), N )
C
      DO 10 I = 1, N - 1
         CALL DGEMM( 'Transpose', 'NoTranspose', N, L, N, ONE, A, LDA,
     $               DWORK(CA+(I-1)*N*L), N, ZERO, DWORK(CA+I*N*L), N )
   10  CONTINUE
C
C     Compute the QR factorization.
C
C     Workspace: need   N*N*L + N + L*N.
C                prefer N*N*L + N + NB*L*N.
C
      ITAU  = CA   + N*N*L
      JWORK = ITAU + N
      CALL DGEQRF( N, L*N, DWORK(CA), N, DWORK(ITAU), DWORK(JWORK),
     $             LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C     Compute Q such that R has all diagonal elements nonnegative.
C     Only the first N*N part of R is needed. Move the details
C     of the QR factorization process, to gain memory and efficiency.
C
C     Workspace: need   2*N*N + 2*N.
C                prefer 2*N*N + N + NB*N.
C
      IR = N*N + 1
      IF ( L.NE.2 )
     $   CALL DCOPY( N, DWORK(ITAU), 1, DWORK(IR+N*N), 1 )
      CALL DLACPY( 'Lower', N, N, DWORK(CA), N, DWORK(IR), N )
      ITAU  = IR + N*N
      JWORK = ITAU + N
C
      IQ = 1
      CALL DLASET( 'Full', N, N, ZERO, ONE, DWORK(IQ), N )
C
      DO 20 I = 1, N
         IF ( DWORK(IR+(I-1)*(N+1)).LT.ZERO )
     $        DWORK(IQ+(I-1)*(N+1))= -ONE
   20 CONTINUE
C
      CALL DORMQR( 'Left', 'NoTranspose', N, N, N, DWORK(IR), N,
     $             DWORK(ITAU), DWORK(IQ), N, DWORK(JWORK),
     $             LDWORK-JWORK+1, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
      JWORK  = IR
C
C     Now, the transformation matrix Q is in DWORK(IQ).
C
C     Compute A = Q'*A*Q.
C
      CALL DGEMM( 'Transpose', 'NoTranspose', N, N, N, ONE, DWORK(IQ),
     $            N, A, LDA, ZERO, DWORK(JWORK), N )
      CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, N, ONE,
     $            DWORK(JWORK), N, DWORK(IQ), N, ZERO, A, LDA )
C
      IF ( M.GT.0 ) THEN
C
C        Compute B = Q'*B.
C        Workspace: need   N*N + N*M.
C
         CALL DLACPY( 'Full', N, M, B, LDB, DWORK(JWORK), N )
         CALL DGEMM(  'Transpose', 'NoTranspose', N, M, N, ONE,
     $                DWORK(IQ), N, DWORK(JWORK), N, ZERO, B, LDB )
      ENDIF
C
C     Compute C = C*Q.
C     Workspace: need   N*N + N*L.
C
      CALL DLACPY( 'Full', L, N, C, LDC, DWORK(JWORK), L )
      CALL DGEMM(  'NoTranspose', 'NoTranspose', L, N, N, ONE,
     $             DWORK(JWORK), L, DWORK(IQ), N, ZERO, C, LDC )
C
C     Compute x0 = Q'*x0.
C
      CALL DCOPY( N, X0, 1, DWORK(JWORK), 1 )
      CALL DGEMV( 'Transpose', N, N, ONE, DWORK(IQ), N, DWORK(JWORK),
     $            1, ZERO, X0, 1 )
C
C     Now, copy C and A into the workspace to make it easier to read out
C     the corresponding part of THETA, and to apply the transformations.
C
      LDCA = N + L
C
      DO 30 I = 1, N
         CALL DCOPY( L, C(1,I), 1, DWORK(CA+(I-1)*LDCA), 1 )
         CALL DCOPY( N, A(1,I), 1, DWORK(CA+L+(I-1)*LDCA), 1 )
   30 CONTINUE
C
      JWORK = CA + LDCA*N
C
C     The parameters characterizing A and C are extracted in this loop.
C     Workspace: need   N*(N + L + 1).
C
      DO 60 I = 1, N
         CALL DCOPY( L, DWORK(CA+1+(N-I)*(LDCA+1)), 1, THETA((I-1)*L+1),
     $               1 )
         RI = DWORK(CA+(N-I)*(LDCA+1))
         TI = DNRM2( L, THETA((I-1)*L+1), 1 )
C
C        Multiply the part of [C; A] which will be currently transformed
C        with Ui = [ -THETAi, Si; RI, THETAi' ] from the left, without
C        storing Ui. Ui has the size (L+1)-by-(L+1).
C
         CALL DGEMV( 'Transpose', L, N, ONE, DWORK(CA+N-I+1), LDCA,
     $               THETA((I-1)*L+1), 1, ZERO, DWORK(JWORK), 1 )
C
         IF ( TI.GT.ZERO ) THEN
            CALL DGER( L, N, (RI-ONE)/TI/TI, THETA((I-1)*L+1), 1,
     $                 DWORK(JWORK), 1, DWORK(CA+N-I+1), LDCA )
         ELSE
C
C           The call below is for the limiting case.
C
            CALL DGER( L, N, -HALF, THETA((I-1)*L+1), 1,
     $                 DWORK(JWORK), 1, DWORK(CA+N-I+1), LDCA )
         ENDIF
C
         CALL DGER( L, N, -ONE, THETA((I-1)*L+1), 1, DWORK(CA+N-I),
     $              LDCA, DWORK(CA+N-I+1), LDCA )
         CALL DAXPY( N, RI, DWORK(CA+N-I), LDCA, DWORK(JWORK), 1 )
C
C        Move these results to their appropriate locations.
C
         DO 50 J = 1, N
            IN = CA + N - I + ( J - 1 )*LDCA
            DO 40 K = IN + 1, IN + L
               DWORK(K-1) = DWORK(K)
   40       CONTINUE
            DWORK(IN+L) = DWORK(JWORK+J-1)
   50    CONTINUE
C
C        Now, apply the bijective mapping, which allows to get rid
C        of the constraint norm(THETAi) < 1.
C
         IF ( LAPPLY .AND. TI.NE.ZERO )
     $     CALL DSCAL( L, TAN( TI*PIBY2 )/TI, THETA((I-1)*L+1), 1 )
C
   60 CONTINUE
C
      IF ( M.GT.0 ) THEN
C
C        The next part of THETA is B.
C
         CALL DLACPY( 'Full', N, M, B, LDB, THETA(N*L+1), N )
C
C        Copy the matrix D.
C
         CALL DLACPY( 'Full', L, M, D, LDD, THETA(N*(L+M)+1), L )
      ENDIF
C
C     Copy the initial state x0.
C
      CALL DCOPY( N, X0, 1, THETA(N*(L+M)+L*M+1), 1 )
C
      DWORK(1) = WRKOPT
      RETURN
C
C *** Last line of TB01VD ***
      END
