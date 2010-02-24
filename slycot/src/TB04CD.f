      SUBROUTINE TB04CD( JOBD, EQUIL, N, M, P, NPZ, A, LDA, B, LDB, C,
     $                   LDC, D, LDD, NZ, LDNZ, NP, LDNP, ZEROSR,
     $                   ZEROSI, POLESR, POLESI, GAINS, LDGAIN, TOL,
     $                   IWORK, DWORK, LDWORK, INFO )
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
C     To compute the transfer function matrix G of a state-space
C     representation (A,B,C,D) of a linear time-invariant multivariable
C     system, using the pole-zeros method. The transfer function matrix
C     is returned in a minimal pole-zero-gain form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears in
C             the given state-space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed to be a zero matrix.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system (A,B,C,D).  N >= 0.
C
C     M       (input) INTEGER
C             The number of the system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of the system outputs.  P >= 0.
C
C     NPZ     (input) INTEGER
C             The maximum number of poles or zeros of the single-input
C             single-output channels in the system. An upper bound
C             for NPZ is N.  NPZ >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state dynamics matrix A.
C             On exit, if EQUIL = 'S', the leading N-by-N part of this
C             array contains the balanced matrix inv(S)*A*S, as returned
C             by SLICOT Library routine TB01ID.
C             If EQUIL = 'N', this array is unchanged on exit.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B.
C             On exit, the contents of B are destroyed: all elements but
C             those in the first row are set to zero.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C.
C             On exit, if EQUIL = 'S', the leading P-by-N part of this
C             array contains the balanced matrix C*S, as returned by
C             SLICOT Library routine TB01ID.
C             If EQUIL = 'N', this array is unchanged on exit.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             If JOBD = 'D', the leading P-by-M part of this array must
C             contain the matrix D.
C             If JOBD = 'Z', the array D is not referenced.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P), if JOBD = 'D';
C             LDD >= 1,        if JOBD = 'Z'.
C
C     NZ      (output) INTEGER array, dimension (LDNZ,M)
C             The leading P-by-M part of this array contains the numbers
C             of zeros of the elements of the transfer function
C             matrix G. Specifically, the (i,j) element of NZ contains
C             the number of zeros of the transfer function G(i,j) from
C             the j-th input to the i-th output.
C
C     LDNZ    INTEGER
C             The leading dimension of array NZ.  LDNZ >= max(1,P).
C
C     NP      (output) INTEGER array, dimension (LDNP,M)
C             The leading P-by-M part of this array contains the numbers
C             of poles of the elements of the transfer function
C             matrix G. Specifically, the (i,j) element of NP contains
C             the number of poles of the transfer function G(i,j).
C
C     LDNP    INTEGER
C             The leading dimension of array NP.  LDNP >= max(1,P).
C
C     ZEROSR  (output) DOUBLE PRECISION array, dimension (P*M*NPZ)
C             This array contains the real parts of the zeros of the
C             transfer function matrix G. The real parts of the zeros
C             are stored in a column-wise order, i.e., for the transfer
C             functions (1,1), (2,1), ..., (P,1), (1,2), (2,2), ...,
C             (P,2), ..., (1,M), (2,M), ..., (P,M); NPZ memory locations
C             are reserved for each transfer function, hence, the real
C             parts of the zeros for the (i,j) transfer function
C             are stored starting from the location ((j-1)*P+i-1)*NPZ+1.
C             Pairs of complex conjugate zeros are stored in consecutive
C             memory locations. Note that only the first NZ(i,j) entries
C             are initialized for the (i,j) transfer function.
C
C     ZEROSI  (output) DOUBLE PRECISION array, dimension (P*M*NPZ)
C             This array contains the imaginary parts of the zeros of
C             the transfer function matrix G, stored in a similar way
C             as the real parts of the zeros.
C
C     POLESR  (output) DOUBLE PRECISION array, dimension (P*M*NPZ)
C             This array contains the real parts of the poles of the
C             transfer function matrix G, stored in the same way as
C             the zeros. Note that only the first NP(i,j) entries are
C             initialized for the (i,j) transfer function.
C
C     POLESI  (output) DOUBLE PRECISION array, dimension (P*M*NPZ)
C             This array contains the imaginary parts of the poles of
C             the transfer function matrix G, stored in the same way as
C             the poles.
C
C     GAINS   (output) DOUBLE PRECISION array, dimension (LDGAIN,M)
C             The leading P-by-M part of this array contains the gains
C             of the transfer function matrix G. Specifically,
C             GAINS(i,j) contains the gain of the transfer function
C             G(i,j).
C
C     LDGAIN  INTEGER
C             The leading dimension of array GAINS.  LDGAIN >= max(1,P).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the
C             controllability of a single-input system (A,b) or (A',c'),
C             where b and c' are columns in B and C' (C transposed). If
C             the user sets TOL > 0, then the given value of TOL is used
C             as an absolute tolerance; elements with absolute value
C             less than TOL are considered neglijible. If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by TOLDEF = N*EPS*MAX( NORM(A), NORM(bc) ) is used
C             instead, where EPS is the machine precision (see LAPACK
C             Library routine DLAMCH), and bc denotes the currently used
C             column in B or C' (see METHOD).
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
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N*(N+P) +
C                              MAX( N + MAX( N,P ), N*(2*N+3)))
C             If N >= P, N >= 1, the formula above can be written as
C             LDWORK >= N*(3*N + P + 3).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the QR algorithm failed to converge when trying to
C                   compute the zeros of a transfer function;
C             = 2:  the QR algorithm failed to converge when trying to
C                   compute the poles of a transfer function.
C                   The errors INFO = 1 or 2 are unlikely to appear.
C
C     METHOD
C
C     The routine implements the pole-zero method proposed in [1].
C     This method is based on an algorithm for computing the transfer
C     function of a single-input single-output (SISO) system.
C     Let (A,b,c,d) be a SISO system. Its transfer function is computed
C     as follows:
C
C     1) Find a controllable realization (Ac,bc,cc) of (A,b,c).
C     2) Find an observable realization (Ao,bo,co) of (Ac,bc,cc).
C     3) Compute the r eigenvalues of Ao (the poles of (Ao,bo,co)).
C     4) Compute the zeros of (Ao,bo,co,d).
C     5) Compute the gain of (Ao,bo,co,d).
C
C     This algorithm can be implemented using only orthogonal
C     transformations [1]. However, for better efficiency, the
C     implementation in TB04CD uses one elementary transformation
C     in Step 4 and r elementary transformations in Step 5 (to reduce
C     an upper Hessenberg matrix to upper triangular form). These
C     special elementary transformations are numerically stable
C     in practice.
C
C     In the multi-input multi-output (MIMO) case, the algorithm
C     computes each element (i,j) of the transfer function matrix G,
C     for i = 1 : P, and for j = 1 : M. For efficiency reasons, Step 1
C     is performed once for each value of j (each column of B). The
C     matrices Ac and Ao result in Hessenberg form.
C
C     REFERENCES
C
C     [1] Varga, A. and Sima, V.
C         Numerically Stable Algorithm for Transfer Function Matrix
C         Evaluation.
C         Int. J. Control, vol. 33, nr. 6, pp. 1123-1133, 1981.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable in practice and requires about
C     20*N**3 floating point operations at most, but usually much less.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2002.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalue, state-space representation, transfer function, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, C100
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, C100 = 100.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          EQUIL, JOBD
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDGAIN, LDNP, LDNZ,
     $                   LDWORK, M, N, NPZ, P
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                   DWORK(*), GAINS(LDGAIN,*), POLESI(*),
     $                   POLESR(*), ZEROSI(*), ZEROSR(*)
      INTEGER            IWORK(*), NP(LDNP,*), NZ(LDNZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION   ANORM, DIJ, EPSN, MAXRED, TOLDEF
      INTEGER            I, IA, IAC, IAS, IB, IC, ICC, IERR, IM, IP,
     $                   IPM1, ITAU, ITAU1, IZ, J, JWK, JWORK, JWORK1,
     $                   K, NCONT, WRKOPT
      LOGICAL            DIJNZ, FNDEIG, WITHD
C     .. Local Arrays ..
      DOUBLE PRECISION   Z(1)
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DHSEQR, DLACPY, MA02AD, TB01ID,
     $                   TB01ZD, TB04BX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      INFO  = 0
      WITHD = LSAME( JOBD, 'D' )
      IF( .NOT.WITHD .AND. .NOT.LSAME( JOBD, 'Z' ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( NPZ.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -14
      ELSE IF( LDNZ.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDNP.LT.MAX( 1, P ) ) THEN
         INFO = -18
      ELSE IF( LDGAIN.LT.MAX( 1, P ) ) THEN
         INFO = -24
      ELSE IF( LDWORK.LT.MAX( 1, N*( N + P ) +
     $                           MAX( N + MAX( N, P ), N*( 2*N + 3 ) ) )
     $       ) THEN
         INFO = -28
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB04CD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      DIJ = ZERO
      IF( MIN( N, P, M ).EQ.0 ) THEN
         IF( MIN( P, M ).GT.0 ) THEN
C
            DO 20 J = 1, M
C
               DO 10 I = 1, P
                  NZ(I,J) = 0
                  NP(I,J) = 0
                  IF ( WITHD )
     $               DIJ = D(I,J)
                  GAINS(I,J) = DIJ
   10          CONTINUE
C
   20       CONTINUE
C
         END IF
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Prepare the computation of the default tolerance.
C
      TOLDEF = TOL
      IF( TOLDEF.LE.ZERO ) THEN
         EPSN  = DBLE( N )*DLAMCH( 'Epsilon' )
         ANORM = DLANGE( 'Frobenius', N, N, A, LDA, DWORK )
      END IF
C
C     Initializations.
C
      IA    = 1
      IC    = IA + N*N
      ITAU  = IC + P*N
      JWORK = ITAU + N
      IAC   = ITAU
C
      K = 1
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(S)*A*S,  B <- inv(S)*B and C <- C*S, where S is a
C        diagonal scaling matrix.
C        Workspace: need   N.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, IERR )
      END IF
C
C     Compute the transfer function matrix of the system (A,B,C,D),
C     in the pole-zero-gain form.
C
      DO 80 J = 1, M
C
C        Save A and C.
C        Workspace: need   W1 = N*(N+P).
C
         CALL DLACPY( 'Full', N, N, A, LDA, DWORK(IA), N )
         CALL DLACPY( 'Full', P, N, C, LDC, DWORK(IC), P )
C
C        Remove the uncontrollable part of the system (A,B(J),C).
C        Workspace: need   W1+N+MAX(N,P);
C                   prefer larger.
C
         CALL TB01ZD( 'No Z', N, P, DWORK(IA), N, B(1,J), DWORK(IC), P,
     $                NCONT, Z, 1, DWORK(ITAU), TOL, DWORK(JWORK),
     $                LDWORK-JWORK+1, IERR )
         IF ( J.EQ.1 )
     $      WRKOPT = INT( DWORK(JWORK) ) + JWORK - 1
C
         IB     = IAC   + NCONT*NCONT
         ICC    = IB    + NCONT
         ITAU1  = ICC   + NCONT
         JWK    = ITAU1 + NCONT
         IAS    = ITAU1
         JWORK1 = IAS   + NCONT*NCONT
C
         DO 70 I = 1, P
            IF ( NCONT.GT.0 ) THEN
               IF ( WITHD )
     $            DIJ = D(I,J)
C
C              Form the matrices of the state-space representation of
C              the dual system for the controllable part.
C              Workspace: need   W2 = W1+N*(N+2).
C
               CALL MA02AD( 'Full', NCONT, NCONT, DWORK(IA), N,
     $                      DWORK(IAC), NCONT )
               CALL DCOPY( NCONT, B(1,J), 1, DWORK(IB), 1 )
               CALL DCOPY( NCONT, DWORK(IC+I-1), P, DWORK(ICC), 1 )
C
C              Remove the unobservable part of the system (A,B(J),C(I)).
C              Workspace: need   W2+2*N;
C                         prefer larger.
C
               CALL TB01ZD( 'No Z', NCONT, 1, DWORK(IAC), NCONT,
     $                      DWORK(ICC), DWORK(IB), 1, IP, Z, 1,
     $                      DWORK(ITAU1), TOL, DWORK(JWK), LDWORK-JWK+1,
     $                      IERR )
               IF ( I.EQ.1 )
     $            WRKOPT = MAX( WRKOPT, INT( DWORK(JWK) ) + JWK - 1 )
C
               IF ( IP.GT.0 ) THEN
C
C                 Save the state matrix of the minimal part.
C                 Workspace: need   W3 = W2+N*N.
C
                  CALL DLACPY( 'Full', IP, IP, DWORK(IAC), NCONT,
     $                         DWORK(IAS), IP )
C
C                 Compute the poles of the transfer function.
C                 Workspace: need   W3+N;
C                            prefer larger.
C
                  CALL DHSEQR( 'Eigenvalues', 'No vectors', IP, 1, IP,
     $                         DWORK(IAC), NCONT, POLESR(K), POLESI(K),
     $                         Z, 1, DWORK(JWORK1), LDWORK-JWORK1+1,
     $                         IERR )
                  IF ( IERR.NE.0 ) THEN
                     INFO = 2
                     RETURN
                  END IF
                  WRKOPT = MAX( WRKOPT,
     $                          INT( DWORK(JWORK1) ) + JWORK1 - 1 )
C
C                 Compute the zeros of the transfer function.
C
                  IPM1   = IP - 1
                  DIJNZ  = WITHD .AND. DIJ.NE.ZERO
                  FNDEIG = DIJNZ .OR. IPM1.GT.0
                  IF ( .NOT.FNDEIG ) THEN
                     IZ = 0
                  ELSE IF ( DIJNZ ) THEN
C
C                    Add the contribution due to D(i,j).
C                    Note that the matrix whose eigenvalues have to
C                    be computed remains in an upper Hessenberg form.
C
                     IZ = IP
                     CALL DLACPY( 'Full', IZ, IZ, DWORK(IAS), IP,
     $                            DWORK(IAC), NCONT )
                     CALL DAXPY( IZ, -DWORK(ICC)/DIJ, DWORK(IB), 1,
     $                           DWORK(IAC), NCONT )
                  ELSE
                     IF( TOL.LE.ZERO )
     $                  TOLDEF = EPSN*MAX( ANORM,
     $                                     DLANGE( 'Frobenius', IP, 1,
     $                                             DWORK(IB), 1, DWORK )
     $                                           )
C
                     DO 30 IM = 1, IPM1
                        IF ( ABS( DWORK(IB+IM-1) ).GT.TOLDEF ) GO TO 40
   30                CONTINUE
C
                     IZ = 0
                     GO TO 50
C
   40                CONTINUE
C
C                    Restore (part of) the saved state matrix.
C
                     IZ = IP - IM
                     CALL DLACPY( 'Full', IZ, IZ, DWORK(IAS+IM*(IP+1)),
     $                             IP, DWORK(IAC), NCONT )
C
C                    Apply the output injection.
C
                     CALL DAXPY( IZ, -DWORK(IAS+IM*(IP+1)-IP)/
     $                           DWORK(IB+IM-1), DWORK(IB+IM), 1,
     $                           DWORK(IAC), NCONT )
                  END IF
C
                  IF ( FNDEIG ) THEN
C
C                    Find the zeros.
C                    Workspace: need   W3+N;
C                               prefer larger.
C
                     CALL DHSEQR( 'Eigenvalues', 'No vectors', IZ, 1,
     $                            IZ, DWORK(IAC), NCONT, ZEROSR(K),
     $                            ZEROSI(K), Z, 1, DWORK(JWORK1),
     $                            LDWORK-JWORK1+1, IERR )
                     IF ( IERR.NE.0 ) THEN
                        INFO = 1
                        RETURN
                     END IF
                  END IF
C
C                 Compute the gain.
C
   50             CONTINUE
                  IF ( DIJNZ ) THEN
                     GAINS(I,J) = DIJ
                  ELSE
                     CALL TB04BX( IP, IZ, DWORK(IAS), IP, DWORK(ICC),
     $                            DWORK(IB), DIJ, POLESR(K), POLESI(K),
     $                            ZEROSR(K), ZEROSI(K), GAINS(I,J),
     $                            IWORK )
                  END IF
                  NZ(I,J) = IZ
                  NP(I,J) = IP
               ELSE
C
C                 Null element.
C
                  NZ(I,J) = 0
                  NP(I,J) = 0
               END IF
C
            ELSE
C
C              Null element.
C
               NZ(I,J) = 0
               NP(I,J) = 0
            END IF
C
            K = K + NPZ
   70    CONTINUE
C
   80 CONTINUE
C
      RETURN
C *** Last line of TB04CD ***
      END
