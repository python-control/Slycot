      DOUBLE PRECISION FUNCTION AB13AX( DICO, N, M, P, A, LDA, B, LDB,
     $                                  C, LDC, HSV, DWORK, LDWORK,
     $                                  INFO )
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
C     To compute the Hankel-norm of the transfer-function matrix G of
C     a stable state-space system (A,B,C). The state dynamics matrix A
C     of the given system is an upper quasi-triangular matrix in
C     real Schur form.
C
C     FUNCTION VALUE
C
C     AB13AX  DOUBLE PRECISION
C             The Hankel-norm of G (if INFO = 0).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state-space representation, i.e. the
C             order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A in a real Schur canonical form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, this array contains the Hankel singular
C             values of the given system ordered decreasingly.
C             HSV(1) is the Hankel norm of the given system.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*(MAX(N,M,P)+5)+N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the state matrix A is not stable (if DICO = 'C')
C                   or not convergent (if DICO = 'D');
C             = 2:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t)                               (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system, and let G be the corresponding
C     transfer-function matrix. The Hankel-norm of G is computed as the
C     the maximum Hankel singular value of the system (A,B,C).
C     The computation of the Hankel singular values is performed
C     by using the square-root method of [1].
C
C     REFERENCES
C
C     [1] Tombs M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     NUMERICAL ASPECTS
C
C     The implemented method relies on a square-root technique.
C                                     3
C     The algorithms require about 17N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, July 1998.
C     Based on the RASP routine SHANRM.
C
C     REVISIONS
C
C     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest.
C     Feb. 2000, V. Sima, Research Institute for Informatics, Bucharest.
C     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     Multivariable system, state-space model, system norms.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO
      INTEGER           INFO, LDA, LDB, LDC, LDWORK, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           DISCR
      INTEGER           I, IERR, J, KR, KS, KTAU, KU, KW, MNMP
      DOUBLE PRECISION  SCALEC, SCALEO, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DSCAL, DTPMV, MA02DD, MB03UD, SB03OU,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C     .. Executable Statements ..
C
      INFO  = 0
      DISCR = LSAME( DICO, 'D' )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDWORK.LT.MAX( 1, N*( MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -13
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         AB13AX = ZERO
         CALL XERBLA( 'AB13AX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         AB13AX = ZERO
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Allocate N*MAX(N,M,P), N, and N*(N+1)/2 working storage for the
C     matrices S, TAU, and R, respectively. S shares the storage with U.
C
      KU   = 1
      KS   = 1
      MNMP = MAX( N, M, P )
      KTAU = KS + N*MNMP
      KR   = KTAU + N
      KW   = KR
C
C     Copy C in U.
C
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK(KU), MNMP )
C
C     If DISCR = .FALSE., solve for R the Lyapunov equation
C                                  2
C     A'*(R'*R) + (R'*R)*A + scaleo  * C'*C = 0 .
C
C     If DISCR = .TRUE., solve for R the Lyapunov equation
C                         2
C     A'*(R'*R)*A + scaleo  * C'*C = R'*R .
C
C     Workspace needed:      N*(MAX(N,M,P)+1);
C     Additional workspace:  need   4*N;
C                            prefer larger.
C
      CALL SB03OU( DISCR, .FALSE., N, P, A, LDA, DWORK(KU), MNMP,
     $             DWORK(KTAU), DWORK(KU), N, SCALEO, DWORK(KW),
     $             LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 1
         RETURN
      ENDIF
C
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
C     Pack the upper triangle of R in DWORK(KR).
C     Workspace needed:      N*(MAX(N,M,P) + 1) + N*(N+1)/2.
C
      CALL MA02DD( 'Pack', 'Upper', N, DWORK(KU), N, DWORK(KR) )
C
      KW = KR + ( N*( N + 1 ) )/2
C
C     Copy B in S (over U).
C
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK(KS), N )
C
C     If DISCR = .FALSE., solve for S the Lyapunov equation
C                                  2
C     A*(S*S') + (S*S')*A' + scalec *B*B' = 0 .
C
C     If DISCR = .TRUE., solve for S the Lyapunov equation
C                         2
C     A*(S*S')*A' + scalec *B*B' = S*S' .
C
C     Workspace needed:      N*(MAX(N,M,P) + 1) + N*(N+1)/2;
C     Additional workspace:  need   4*N;
C                            prefer larger.
C
      CALL SB03OU( DISCR, .TRUE., N, M, A, LDA, DWORK(KS), N,
     $             DWORK(KTAU), DWORK(KS), N, SCALEC, DWORK(KW),
     $             LDWORK-KW+1, IERR )
C
      WRKOPT = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
C                             | x x |
C     Compute R*S in the form | 0 x | in S. Note that R is packed.
C
      J = KS
      DO 10 I = 1, N
         CALL DTPMV( 'Upper', 'NoTranspose', 'NonUnit', I, DWORK(KR),
     $               DWORK(J), 1 )
         J = J + N
   10 CONTINUE
C
C     Compute the singular values of the upper triangular matrix R*S.
C
C     Workspace needed:      N*MAX(N,M,P);
C     Additional workspace:  need   MAX(1,5*N);
C                            prefer larger.
C
      KW = KTAU
      CALL MB03UD( 'NoVectors', 'NoVectors', N, DWORK(KS), N, DWORK, 1,
     $             HSV, DWORK(KW), LDWORK-KW+1, IERR )
      IF( IERR.NE.0 ) THEN
         INFO = 2
         RETURN
      ENDIF
C
C     Scale singular values.
C
      CALL DSCAL( N, ONE / SCALEC / SCALEO, HSV, 1 )
      AB13AX = HSV(1)
C
      DWORK(1) = MAX( WRKOPT, DWORK(KW) + DBLE( KW-1 ) )
C
      RETURN
C *** Last line of AB13AX ***
      END
