      SUBROUTINE NF01BD( CJTE, NSMP, M, L, IPAR, LIPAR, X, LX, U, LDU,
     $                   E, J, LDJ, JTE, DWORK, LDWORK, INFO )
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
C     To calculate the Jacobian dy/dX of the Wiener system
C
C        x(t+1) = A*x(t) + B*u(t)
C        z(t)   = C*x(t) + D*u(t),
C
C        y(t,i) = sum( ws(k, i)*f(w(k, i)*z(t) + b(k,i)) ) + b(k+1,i),
C
C     where t = 1, 2, ...,  NSMP,
C           i = 1, 2, ...,  L,
C           k = 1, 2, ...,  NN.
C
C     NN is arbitrary eligible and has to be provided in IPAR(2), and
C     X = ( wb(1), ..., wb(L), theta ) is described below.
C
C     Denoting y(j) = y(1:NSMP,j), the Jacobian J has the block form
C
C       dy(1)/dwb(1)       0         .....       0         dy(1)/dtheta
C            0        dy(2)/dwb(2)   .....       0         dy(2)/dtheta
C          .....         .....       .....     .....          .....
C            0           .....         0    dy(L)/dwb(L)   dy(L)/dtheta
C
C     but it will be returned without the zero blocks, in the form
C
C     dy(1)/dwb(1)    dy(1)/dtheta
C                  ...
C     dy(L)/dwb(L)    dy(L)/dtheta.
C
C     dy(i)/dwb(i) depends on f and is calculated by the routine NF01BY;
C     dy(i)/dtheta is computed by a forward-difference approximation.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     CJTE    CHARACTER*1
C             Specifies whether the matrix-vector product J'*e should be
C             computed or not, as follows:
C             = 'C' :  compute J'*e;
C             = 'N' :  do not compute J'*e.
C
C     Input/Output Parameters
C
C     NSMP    (input) INTEGER
C             The number of training samples.  NSMP >= 0.
C
C     M       (input) INTEGER
C             The length of each input sample.  M >= 0.
C
C     L       (input) INTEGER
C             The length of each output sample.  L >= 0.
C
C     IPAR    (input/output) INTEGER array, dimension (LIPAR)
C             On entry, the first entries of this array must contain
C             the integer parameters needed; specifically,
C             IPAR(1)  must contain the order of the linear part, N;
C                      actually, N = abs(IPAR(1)), since setting
C                      IPAR(1) < 0 has a special meaning (see below);
C             IPAR(2)  must contain the number of neurons for the
C                      nonlinear part, NN, NN >= 0.
C             On exit, if IPAR(1) < 0 on entry, then no computations are
C             performed, except the needed tests on input parameters,
C             but the following values are returned:
C             IPAR(1) contains the length of the array J, LJ;
C             LDJ     contains the leading dimension of array J.
C             Otherwise, IPAR(1) and LDJ are unchanged on exit.
C
C     LIPAR   (input) INTEGER
C             The length of the array IPAR.  LIPAR >= 2.
C
C     X       (input) DOUBLE PRECISION array, dimension (LX)
C             The leading LPAR entries of this array must contain the
C             set of system parameters, where
C                LPAR = (NN*(L + 2) + 1)*L + N*(M + L + 1) + L*M.
C             X has the form (wb(1), ..., wb(L), theta), where the
C             vectors wb(i) have the structure
C              (w(1,1), ..., w(1,L), ..., w(NN,1), ..., w(NN,L),
C                ws(1), ..., ws(NN), b(1), ..., b(NN+1) ),
C             and the vector theta represents the matrices A, B, C, D
C             and x(1), and it can be retrieved from these matrices
C             by SLICOT Library routine TB01VD and retranslated by
C             TB01VY.
C
C     LX      (input) INTEGER
C             The length of X.
C             LX >= (NN*(L + 2) + 1)*L + N*(M + L + 1) + L*M.
C
C     U       (input) DOUBLE PRECISION array, dimension (LDU, M)
C             The leading NSMP-by-M part of this array must contain the
C             set of input samples,
C             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ).
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= MAX(1,NSMP).
C
C     E       (input) DOUBLE PRECISION array, dimension (NSMP*L)
C             If CJTE = 'C', this array must contain a vector e, which
C             will be premultiplied with J', e = vec( Y - y ), where
C             Y is set of output samples, and vec denotes the
C             concatenation of the columns of a matrix.
C             If CJTE = 'N', this array is not referenced.
C
C     J       (output) DOUBLE PRECISION array, dimension (LDJ, *)
C             The leading NSMP*L-by-NCOLJ part of this array contains
C             the Jacobian of the error function stored in a compressed
C             form, as described above, where
C             NCOLJ = NN*(L + 2) + 1 + N*(M + L + 1) + L*M.
C
C     LDJ     INTEGER
C             The leading dimension of array J.  LDJ >= MAX(1,NSMP*L).
C             Note that LDJ is an input parameter, except for
C             IPAR(1) < 0 on entry, when it is an output parameter.
C
C     JTE     (output) DOUBLE PRECISION array, dimension (LPAR)
C             If CJTE = 'C', this array contains the matrix-vector
C             product J'*e.
C             If CJTE = 'N', this array is not referenced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= 2*NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N +
C                                       MAX( N*(N + L), N + M + L ) )
C                                                              if M > 0;
C             LDWORK >= 2*NSMP*L + MAX( 2*NN, (N + L)*N + 2*N +
C                                       MAX( N*(N + L), L ) ), if M = 0.
C             A larger value of LDWORK could improve the efficiency.
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
C     BLAS routines are used for the matrix-vector multiplications, and
C     the SLICOT Library routine TB01VY is called for the conversion of
C     the output normal form parameters to an LTI-system; the routine
C     NF01AD is then used for the simulation of the system with given
C     parameters, and the routine NF01BY is called for the (analytically
C     performed) calculation of the parts referring to the parameters
C     of the static nonlinearity.
C
C     CONTRIBUTORS
C
C     A. Riedel, R. Schneider, Chemnitz University of Technology,
C     Mar. 2001, during a stay at University of Twente, NL.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001,
C     Dec. 2001.
C
C     KEYWORDS
C
C     Jacobian matrix, nonlinear system, output normal form, simulation,
C     state-space representation, Wiener system.
C
C     ******************************************************************
C
C     .. Parameters ..
C     .. EPSFCN is related to the error in computing the functions ..
C     .. For EPSFCN = 0.0D0, the square root of the machine precision
C     .. is used for finite difference approximation of the derivatives.
      DOUBLE PRECISION  ZERO, ONE, EPSFCN
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, EPSFCN = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         CJTE
      INTEGER           INFO, L, LDJ, LDU, LDWORK, LX, LIPAR, M, NSMP
C     .. Array Arguments ..
      INTEGER           IPAR(*)
      DOUBLE PRECISION  DWORK(*), E(*), J(LDJ, *), JTE(*), U(LDU,*),
     $                  X(*)
C     .. Local Scalars ..
      LOGICAL           WJTE
      DOUBLE PRECISION  EPS, H, PARSAV
      INTEGER           AC, BD, BSN, I, IX, IY, JW, K, KCOL, LDAC, LPAR,
     $                  LTHS, N, NN, NSML, NTHS, Z
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      LOGICAL           LSAME
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, NF01AD, NF01AY, NF01BY, TB01VY,
     $                  TF01MX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
      N    = IPAR(1)
      NN   = IPAR(2)
      BSN  = NN*( L + 2 ) + 1
      NSML = NSMP*L
      NTHS = BSN*L
      LTHS = N*( M + L + 1 ) + L*M
      LPAR = NTHS + LTHS
      WJTE = LSAME( CJTE, 'C' )
C
C     Check the scalar input parameters.
C
      INFO = 0
      IF( .NOT.( WJTE .OR. LSAME( CJTE, 'N' ) ) ) THEN
         INFO = -1
      ELSEIF ( NSMP.LT.0 ) THEN
         INFO = -2
      ELSEIF ( M.LT.0 ) THEN
         INFO = -3
      ELSEIF ( L.LT.0 ) THEN
         INFO = -4
      ELSEIF ( NN.LT.0 ) THEN
         INFO = -5
      ELSEIF ( LIPAR.LT.2 ) THEN
         INFO = -6
      ELSEIF ( IPAR(1).LT.0 ) THEN
         IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'NF01BD', -INFO )
         ELSE
            IPAR(1) = NSML*( ABS( N )*( M + L + 1 ) + L*M + BSN )
            LDJ     = MAX( 1, NSML )
         ENDIF
         RETURN
      ELSEIF ( LX.LT.LPAR ) THEN
         INFO = -8
      ELSEIF ( LDU.LT.MAX( 1, NSMP ) ) THEN
         INFO = -10
      ELSEIF ( LDJ.LT.MAX( 1, NSML ) ) THEN
         INFO = -13
      ELSE
         LDAC = N + L
         IF ( M.GT.0 ) THEN
            JW = MAX( N*LDAC, N + M + L )
         ELSE
            JW = MAX( N*LDAC, L )
         END IF
         IF ( LDWORK.LT.2*NSML + MAX( 2*NN, LDAC*( N + M ) + 2*N + JW ))
     $      INFO = -16
      ENDIF
C
C     Return if there are illegal arguments.
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'NF01BD', -INFO )
         RETURN
      ENDIF
C
C     Quick return if possible.
C
      IF ( MIN( NSMP, L ).EQ.0 ) THEN
         IF ( WJTE .AND. LPAR.GE.1 ) THEN
            JTE(1) = ZERO
            CALL DCOPY( LPAR, JTE(1), 0, JTE(1), 1 )
         END IF
         RETURN
      END IF
C
C     Compute the output of the linear part.
C     Workspace: need  2*NSMP*L + (N + L)*(N + M) + N + N*(N + L + 1).
C     (2*NSMP*L locations are reserved for computing two times the
C     output of the linear part.)
C
      IY = 1
      Z  = IY + NSML
      AC = Z  + NSML
      BD = AC + LDAC*N
      IX = BD + LDAC*M
      JW = IX + N
C
      CALL TB01VY( 'Apply', N, M, L, X(NTHS+1), LTHS, DWORK(AC), LDAC,
     $             DWORK(BD), LDAC, DWORK(AC+N), LDAC, DWORK(BD+N),
     $             LDAC, DWORK(IX), DWORK(JW), LDWORK-JW+1, INFO )
C
C     Workspace: need   2*NSMP*L + (N + L)*(N + M) + 3*N + M + L,
C                                                             if M > 0;
C                       2*NSMP*L + (N + L)*N + 2*N + L,       if M = 0;
C                prefer larger.
C
      CALL TF01MX( N, M, L, NSMP, DWORK(AC), LDAC, U, LDU, DWORK(IX),
     $             DWORK(Z), NSMP, DWORK(JW), LDWORK-JW+1, INFO )
C
C     Fill the blocks dy(i)/dwb(i) and the corresponding parts of JTE,
C     if needed.
C
      JW = AC
      IF ( WJTE ) THEN
C
         DO 10 I = 0, L - 1
            CALL NF01BY( CJTE, NSMP, L, 1, IPAR(2), LIPAR-1, X(I*BSN+1),
     $                   BSN, DWORK(Z), NSMP, E(I*NSMP+1),
     $                   J(I*NSMP+1,1), LDJ, JTE(I*BSN+1), DWORK(JW),
     $                   LDWORK-JW+1, INFO )
   10    CONTINUE
C
      ELSE
C
         DO 20 I = 0, L - 1
            CALL NF01BY( CJTE, NSMP, L, 1, IPAR(2), LIPAR-1, X(I*BSN+1),
     $                   BSN, DWORK(Z), NSMP, DWORK, J(I*NSMP+1,1), LDJ,
     $                   DWORK, DWORK(JW), LDWORK-JW+1, INFO )
   20    CONTINUE
C
      END IF
C
C     Compute the output of the system with unchanged parameters.
C     Workspace: need   2*NSMP*L + 2*NN;
C                prefer larger.
C
      CALL NF01AY( NSMP, L, L, IPAR(2), LIPAR-1, X, NTHS, DWORK(Z),
     $             NSMP, DWORK(IY), NSMP, DWORK(JW), LDWORK-JW+1,
     $             INFO )
C
C     Compute dy/dtheta numerically by forward-difference approximation.
C     Workspace: need   2*NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N +
C                                       MAX( N*(N + L), N + M + L ) ),
C                                                              if M > 0;
C                       2*NSMP*L + MAX( 2*NN, (N + L)*N + 2*N +
C                                       MAX( N*(N + L), L ) ), if M = 0;
C                prefer larger.
C
      JW  = Z
      EPS = SQRT( MAX( EPSFCN, DLAMCH( 'Epsilon' ) ) )
C
      DO 40 K = NTHS + 1, LPAR
         KCOL   = K - NTHS + BSN
         PARSAV = X(K)
         IF ( PARSAV.EQ.ZERO ) THEN
            H = EPS
         ELSE
            H = EPS*ABS( PARSAV )
         END IF
         X(K) = X(K) + H
         CALL NF01AD( NSMP, M, L, IPAR, LIPAR, X, LPAR, U, LDU,
     $                J(1,KCOL), NSMP, DWORK(JW), LDWORK-JW+1,
     $                INFO )
         X(K) = PARSAV
C
         DO 30 I = 1, NSML
            J(I,KCOL) = ( J(I,KCOL) - DWORK(I) ) / H
   30    CONTINUE
C
   40 CONTINUE
C
      IF ( WJTE ) THEN
C
C        Compute the last part of J'e in JTE.
C
         CALL DGEMV( 'Transpose', NSML, LTHS, ONE, J(1,BSN+1), LDJ, E,
     $               1, ZERO, JTE(NTHS+1), 1 )
      END IF
C
      RETURN
C
C *** Last line of NF01BD ***
      END
