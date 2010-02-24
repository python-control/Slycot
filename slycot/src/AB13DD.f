      SUBROUTINE AB13DD( DICO, JOBE, EQUIL, JOBD, N, M, P, FPEAK,
     $                   A, LDA, E, LDE, B, LDB, C, LDC, D, LDD, GPEAK,
     $                   TOL, IWORK, DWORK, LDWORK, CWORK, LCWORK,
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
C     To compute the L-infinity norm of a continuous-time or
C     discrete-time system, either standard or in the descriptor form,
C
C                                     -1
C        G(lambda) = C*( lambda*E - A ) *B + D .
C
C     The norm is finite if and only if the matrix pair (A,E) has no
C     eigenvalue on the boundary of the stability domain, i.e., the
C     imaginary axis, or the unit circle, respectively. It is assumed
C     that the matrix E is nonsingular.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system, as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOBE    CHARACTER*1
C             Specifies whether E is a general square or an identity
C             matrix, as follows:
C             = 'G':  E is a general square matrix;
C             = 'I':  E is the identity matrix.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the system (A,E,B,C) or (A,B,C), as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears in
C             the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The row size of the matrix C.  P >= 0.
C
C     FPEAK   (input/output) DOUBLE PRECISION array, dimension (2)
C             On entry, this parameter must contain an estimate of the
C             frequency where the gain of the frequency response would
C             achieve its peak value. Setting FPEAK(2) = 0 indicates an
C             infinite frequency. An accurate estimate could reduce the
C             number of iterations of the iterative algorithm. If no
C             estimate is available, set FPEAK(1) = 0, and FPEAK(2) = 1.
C             FPEAK(1) >= 0, FPEAK(2) >= 0.
C             On exit, if INFO = 0, this array contains the frequency
C             OMEGA, where the gain of the frequency response achieves
C             its peak value GPEAK, i.e.,
C
C                 || G ( j*OMEGA ) || = GPEAK ,  if DICO = 'C', or
C
C                         j*OMEGA
C                 || G ( e       ) || = GPEAK ,  if DICO = 'D',
C
C             where OMEGA = FPEAK(1), if FPEAK(2) > 0, and OMEGA is
C             infinite, if FPEAK(2) = 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             If JOBE = 'G', the leading N-by-N part of this array must
C             contain the descriptor matrix E of the system.
C             If JOBE = 'I', then E is assumed to be the identity
C             matrix and is not referenced.
C
C     LDE     INTEGER
C             The leading dimension of the array E.
C             LDE >= MAX(1,N), if JOBE = 'G';
C             LDE >= 1,        if JOBE = 'I'.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             If JOBD = 'D', the leading P-by-M part of this array must
C             contain the direct transmission matrix D.
C             The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P), if JOBD = 'D';
C             LDD >= 1,        if JOBD = 'Z'.
C
C     GPEAK   (output) DOUBLE PRECISION array, dimension (2)
C             The L-infinity norm of the system, i.e., the peak gain
C             of the frequency response (as measured by the largest
C             singular value in the MIMO case), coded in the same way
C             as FPEAK.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used to set the accuracy in determining the
C             norm.  0 <= TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= K, where K can be computed using the following
C             pseudo-code (or the Fortran code included in the routine)
C
C                d = 6*MIN(P,M);
C                c = MAX( 4*MIN(P,M) + MAX(P,M), d );
C                if ( MIN(P,M) = 0 ) then
C                   K = 1;
C                else if( N = 0 or B = 0 or C = 0 ) then
C                   if( JOBD = 'D' ) then
C                      K = P*M + c;
C                   else
C                      K = 1;
C                   end
C                else
C                   if ( DICO = 'D' ) then
C                      b = 0;  e = d;
C                   else
C                      b = N*(N+M);  e = c;
C                      if ( JOBD = Z' ) then  b = b + P*M;  end
C                   end
C                   if ( JOBD = 'D' ) then
C                      r = P*M;
C                      if ( JOBE = 'I', DICO = 'C',
C                           N > 0, B <> 0, C <> 0 ) then
C                         K = P*P + M*M;
C                         r = r + N*(P+M);
C                      else
C                         K = 0;
C                      end
C                      K = K + r + c;  r = r + MIN(P,M);
C                   else
C                      r = 0;  K = 0;
C                   end
C                   r = r + N*(N+P+M);
C                   if ( JOBE = 'G' ) then
C                      r = r + N*N;
C                      if ( EQUIL = 'S' ) then
C                         K = MAX( K, r + 9*N );
C                      end
C                      K = MAX( K, r + 4*N + MAX( M, 2*N*N, N+b+e ) );
C                   else
C                      K = MAX( K, r + N +
C                                  MAX( M, P, N*N+2*N, 3*N+b+e ) );
C                   end
C                   w = 0;
C                   if ( JOBE = 'I', DICO = 'C' ) then
C                      w = r + 4*N*N + 11*N;
C                      if ( JOBD = 'D' ) then
C                         w = w + MAX(M,P) + N*(P+M);
C                      end
C                   end
C                   if ( JOBE = 'E' or DICO = 'D' or JOBD = 'D' ) then
C                      w = MAX( w, r + 6*N + (2*N+P+M)*(2*N+P+M) +
C                               MAX( 2*(N+P+M), 8*N*N + 16*N ) );
C                   end
C                   K = MAX( 1, K, w, r + 2*N + e );
C                end
C
C             For good performance, LDWORK must generally be larger.
C
C             An easily computable upper bound is
C
C             K = MAX( 1, 15*N*N + P*P + M*M + (6*N+3)*(P+M) + 4*P*M +
C                         N*M + 22*N + 7*MIN(P,M) ).
C
C             The smallest workspace is obtained for DICO = 'C',
C             JOBE = 'I', and JOBD = 'Z', namely
C
C             K = MAX( 1, N*N + N*P + N*M + N +
C                         MAX( N*N + N*M + P*M + 3*N + c,
C                              4*N*N + 10*N ) ).
C
C             for which an upper bound is
C
C             K = MAX( 1, 6*N*N + N*P + 2*N*M + P*M + 11*N + MAX(P,M) +
C                         6*MIN(P,M) ).
C
C     CWORK   COMPLEX*16 array, dimension (LCWORK)
C             On exit, if INFO = 0, CWORK(1) contains the optimal
C             LCWORK.
C
C     LCWORK  INTEGER
C             The dimension of the array CWORK.
C             LCWORK >= 1,  if N = 0, or B = 0, or C = 0;
C             LCWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)),
C                           otherwise.
C             For good performance, LCWORK must generally be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the matrix E is (numerically) singular;
C             = 2:  the (periodic) QR (or QZ) algorithm for computing
C                   eigenvalues did not converge;
C             = 3:  the SVD algorithm for computing singular values did
C                   not converge;
C             = 4:  the tolerance is too small and the algorithm did
C                   not converge.
C
C     METHOD
C
C     The routine implements the method presented in [1], with
C     extensions and refinements for improving numerical robustness and
C     efficiency. Structure-exploiting eigenvalue computations for
C     Hamiltonian matrices are used if JOBE = 'I', DICO = 'C', and the
C     symmetric matrices to be implicitly inverted are not too ill-
C     conditioned. Otherwise, generalized eigenvalue computations are
C     used in the iterative algorithm of [1].
C
C     REFERENCES
C
C     [1] Bruinsma, N.A. and Steinbuch, M.
C         A fast algorithm to compute the Hinfinity-norm of a transfer
C         function matrix.
C         Systems & Control Letters, vol. 14, pp. 287-293, 1990.
C
C     NUMERICAL ASPECTS
C
C     If the algorithm does not converge in MAXIT = 30 iterations
C     (INFO = 4), the tolerance must be increased.
C
C     FURTHER COMMENTS
C
C     If the matrix E is singular, other SLICOT Library routines
C     could be used before calling AB13DD, for removing the singular
C     part of the system.
C
C     CONTRIBUTORS
C
C     D. Sima, University of Bucharest, May 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2001.
C     Partly based on SLICOT Library routine AB13CD by P.Hr. Petkov,
C     D.W. Gu and M.M. Konstantinov.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001,
C     May 2003, Aug. 2005, March 2008, May 2009, Sep. 2009.
C
C     KEYWORDS
C
C     H-infinity optimal control, robust control, system norm.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR, P25
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                     FOUR = 4.0D+0, P25 = 0.25D+0 )
      DOUBLE PRECISION   TEN, HUNDRD, THOUSD
      PARAMETER          ( TEN    = 1.0D+1, HUNDRD = 1.0D+2,
     $                     THOUSD = 1.0D+3 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          DICO, EQUIL, JOBD, JOBE
      INTEGER            INFO, LCWORK, LDA, LDB, LDC, LDD, LDE, LDWORK,
     $                   M, N, P
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      COMPLEX*16         CWORK(  * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), DWORK(  * ), E( LDE, * ),
     $                   FPEAK(  2 ), GPEAK(  2 )
      INTEGER            IWORK(  * )
C     ..
C     .. Local Scalars ..
      CHARACTER          VECT
      LOGICAL            DISCR, FULLE, ILASCL, ILESCL, LEQUIL, NODYN,
     $                   USEPEN, WITHD
      INTEGER            I, IA, IAR, IAS, IB, IBS, IBT, IBV, IC, ICU,
     $                   ID, IE, IERR, IES, IH, IH12, IHI, II, ILO, IM,
     $                   IMIN, IPA, IPE, IR, IS, ISB, ISC, ISL, ITAU,
     $                   ITER, IU, IV, IWRK, J, K, LW, MAXCWK, MAXWRK,
     $                   MINCWR, MINPM, MINWRK, N2, N2PM, NEI, NN, NWS,
     $                   NY, PM
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNORM, BOUND, CNORM,
     $                   ENRM, ENRMTO, EPS, FPEAKI, FPEAKS, GAMMA,
     $                   GAMMAL, GAMMAS, MAXRED, OMEGA, PI, RAT, RCOND,
     $                   RTOL, SAFMAX, SAFMIN, SMLNUM, TM, TOLER, WMAX,
     $                   WRMIN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION   TEMP( 1 )
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   AB13DX, DLAMCH, DLANGE, DLAPY2
      LOGICAL            LSAME
      EXTERNAL           AB13DX, DLAMCH, DLANGE, DLAPY2, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEBAL, DGEHRD, DGEMM, DGEQRF, DGESVD,
     $                   DGGBAL, DGGEV, DHGEQZ, DHSEQR, DLABAD, DLACPY,
     $                   DLASCL, DLASRT, DORGQR, DORMHR, DSWAP, DSYRK,
     $                   DTRCON, MA02AD, MB01SD, MB03XD, TB01ID, TG01AD,
     $                   TG01BD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, ATAN, ATAN2, COS, DBLE, INT, LOG, MAX,
     $                   MIN, SIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      N2     = 2*N
      NN     = N*N
      PM     = P  + M
      N2PM   = N2 + PM
      MINPM  = MIN( P, M )
      INFO   = 0
      DISCR  = LSAME( DICO,  'D' )
      FULLE  = LSAME( JOBE,  'G' )
      LEQUIL = LSAME( EQUIL, 'S' )
      WITHD  = LSAME( JOBD,  'D' )
C
      IF( .NOT. ( DISCR .OR. LSAME( DICO, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( FULLE  .OR. LSAME( JOBE,  'I' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LEQUIL .OR. LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( WITHD  .OR. LSAME( JOBD,  'Z' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( MIN( FPEAK( 1 ), FPEAK( 2 ) ).LT.ZERO ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDE.LT.1 .OR. ( FULLE .AND. LDE.LT.N ) ) THEN
         INFO = -12
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -18
      ELSE IF( TOL.LT.ZERO .OR. TOL.GE.ONE ) THEN
         INFO = -20
      ELSE
         BNORM  = DLANGE( '1-norm', N, M, B, LDB, DWORK )
         CNORM  = DLANGE( '1-norm', P, N, C, LDC, DWORK )
         NODYN  = N.EQ.0 .OR. MIN( BNORM, CNORM ).EQ.ZERO
         USEPEN = FULLE .OR. DISCR
C
C        Compute workspace.
C
         ID = 6*MINPM
         IC = MAX( 4*MINPM + MAX( P, M ), ID )
         IF( MINPM.EQ.0 ) THEN
            MINWRK = 1
         ELSE IF( NODYN ) THEN
            IF( WITHD ) THEN
               MINWRK = P*M + IC
            ELSE
               MINWRK = 1
            END IF
         ELSE
            IF ( DISCR ) THEN
               IB = 0
               IE = ID
            ELSE
               IB = N*( N + M )
               IF ( .NOT.WITHD )
     $            IB = IB + P*M
               IE = IC
            END IF
            IF ( WITHD ) THEN
               IR = P*M
               IF ( .NOT.USEPEN ) THEN
                  MINWRK = P*P + M*M
                  IR = IR + N*PM
               ELSE
                  MINWRK = 0
               END IF
               MINWRK = MINWRK + IR + IC
               IR = IR + MINPM
            ELSE
               IR = 0
               MINWRK = 0
            END IF
            IR = IR + N*( N + PM )
            IF ( FULLE ) THEN
               IR = IR + NN
               IF ( LEQUIL )
     $            MINWRK = MAX( MINWRK, IR + 9*N )
               MINWRK = MAX( MINWRK, IR + 4*N + MAX( M, 2*NN,
     $                                               N + IB + IE ) )
            ELSE
               MINWRK = MAX( MINWRK, IR + N + MAX( M, P, NN + N2,
     $                                             3*N + IB + IE ) )
            END IF
            LW = 0
            IF ( .NOT.USEPEN ) THEN
               LW = IR + 4*NN + 11*N
               IF ( WITHD )
     $            LW = LW + MAX( M, P ) + N*PM
            END IF
            IF ( USEPEN .OR. WITHD )
     $         LW = MAX( LW, IR + 6*N + N2PM*N2PM +
     $                       MAX( N2PM + PM, 8*( NN + N2 ) ) )
            MINWRK = MAX( 1, MINWRK, LW, IR + N2 + IE )
         END IF
C
         IF( LDWORK.LT.MINWRK ) THEN
            INFO = -23
         ELSE
            IF ( NODYN ) THEN
               MINCWR = 1
            ELSE
               MINCWR = MAX( 1, ( N + M )*( N + P ) +
     $                          2*MINPM + MAX( P, M ) )
            END IF
            IF( LCWORK.LT.MINCWR )
     $         INFO = -25
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'AB13DD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( M.EQ.0 .OR. P.EQ.0 ) THEN
         GPEAK( 1 ) = ZERO
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         DWORK( 1 ) = ONE
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine the maximum singular value of G(infinity) = D .
C     If JOBE = 'I' and DICO = 'C', the full SVD of D, D = U*S*V', is
C     computed and saved for later use.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      ID = 1
      IF ( WITHD ) THEN
         IS = ID + P*M
         IF ( USEPEN .OR. NODYN ) THEN
            IU   = IS + MINPM
            IV   = IU
            IWRK = IV
            VECT = 'N'
         ELSE
            IBV  = IS  + MINPM
            ICU  = IBV + N*M
            IU   = ICU + P*N
            IV   = IU  + P*P
            IWRK = IV  + M*M
            VECT = 'A'
         END IF
C
C        Workspace: need   P*M + MIN(P,M) + V +
C                          MAX( 3*MIN(P,M) + MAX(P,M), 5*MIN(P,M) ),
C                          where V = N*(M+P) + P*P + M*M,
C                                        if JOBE = 'I' and DICO = 'C',
C                                        and N > 0, B <> 0, C <> 0,
C                                V = 0,  otherwise;
C                   prefer larger.
C
         CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
         CALL DGESVD( VECT, VECT, P, M, DWORK( ID ), P, DWORK( IS ),
     $                DWORK( IU ), P, DWORK( IV ), M, DWORK( IWRK ),
     $                LDWORK-IWRK+1, IERR )
         IF( IERR.GT.0 ) THEN
            INFO = 3
            RETURN
         END IF
         GAMMAL = DWORK( IS )
         MAXWRK = INT( DWORK( IWRK ) ) + IWRK - 1
C
C        Restore D for later calculations.
C
         CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
      ELSE
         IWRK   = 1
         GAMMAL = ZERO
         MAXWRK = 1
      END IF
C
C     Quick return if possible.
C
      IF( NODYN ) THEN
         GPEAK( 1 ) = GAMMAL
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         DWORK( 1 ) = MAXWRK
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
      IF ( .NOT.USEPEN .AND. WITHD ) THEN
C
C        Standard continuous-time case, D <> 0: Compute B*V and C'*U .
C
         CALL DGEMM( 'No Transpose', 'Transpose', N, M, M, ONE, B, LDB,
     $               DWORK( IV ), M, ZERO, DWORK( IBV ), N )
         CALL DGEMM( 'Transpose', 'No Transpose', N, P, P, ONE, C,
     $               LDC, DWORK( IU ), P, ZERO, DWORK( ICU ), N )
C
C        U and V are no longer needed: free their memory space.
C        Total workspace here: need   P*M + MIN(P,M) + N*(M+P)
C        (JOBE = 'I', DICO = 'C', JOBD = 'D').
C
         IWRK = IU
      END IF
C
C     Get machine constants.
C
      EPS    = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SQRT( SAFMIN ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      TOLER  = SQRT( EPS )
C
C     Initiate the transformation of the system to an equivalent one,
C     to be used for eigenvalue computations.
C
C     Additional workspace: need   N*N + N*M + P*N + 2*N, if JOBE = 'I';
C                                2*N*N + N*M + P*N + 2*N, if JOBE = 'G'.
C
      IA = IWRK
      IE = IA + NN
      IF ( FULLE ) THEN
         IB = IE + NN
      ELSE
         IB = IE
      END IF
      IC  = IB + N*M
      IR  = IC + P*N
      II  = IR + N
      IBT = II + N
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IA ), N )
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK( IB ), N )
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK( IC ), P )
C
C     Scale A if maximum element is outside the range [SMLNUM,BIGNUM].
C
      ANRM = DLANGE( 'Max', N, N, DWORK( IA ), N, DWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL DLASCL( 'General', 0, 0, ANRM, ANRMTO, N, N, DWORK( IA ),
     $                N, IERR )
C
      IF ( FULLE ) THEN
C
C        Descriptor system.
C
C        Additional workspace: need   N.
C
         IWRK = IBT + N
         CALL DLACPY( 'Full', N, N, E, LDE, DWORK( IE ), N )
C
C        Scale E if maximum element is outside the range
C        [SMLNUM,BIGNUM].
C
         ENRM = DLANGE( 'Max', N, N, DWORK( IE ), N, DWORK )
         ILESCL = .FALSE.
         IF( ENRM.GT.ZERO .AND. ENRM.LT.SMLNUM ) THEN
            ENRMTO = SMLNUM
            ILESCL = .TRUE.
         ELSE IF( ENRM.GT.BIGNUM ) THEN
            ENRMTO = BIGNUM
            ILESCL = .TRUE.
         ELSE IF( ENRM.EQ.ZERO ) THEN
C
C           Error return: Matrix E is 0.
C
            INFO = 1
            RETURN
         END IF
         IF( ILESCL )
     $      CALL DLASCL( 'General', 0, 0, ENRM, ENRMTO, N, N,
     $                   DWORK( IE ), N, IERR )
C
C        Equilibrate the system, if required.
C
C        Additional workspace: need   6*N.
C
         IF( LEQUIL )
     $      CALL TG01AD( 'All', N, N, M, P, ZERO, DWORK( IA ), N,
     $                   DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ), P,
     $                   DWORK( II ), DWORK( IR ), DWORK( IWRK ),
     $                   IERR )
C
C        For efficiency of later calculations, the system (A,E,B,C) is
C        reduced to an equivalent one with the state matrix A in
C        Hessenberg form, and E upper triangular.
C        First, permute (A,E) to make it more nearly triangular.
C
         CALL DGGBAL( 'Permute', N, DWORK( IA ), N, DWORK( IE ), N, ILO,
     $                IHI, DWORK( II ), DWORK( IR ), DWORK( IWRK ),
     $                IERR )
C
C        Apply the permutations to (the copies of) B and C.
C
         DO 10 I = N, IHI + 1, -1
            K = DWORK( II+I-1 )
            IF( K.NE.I )
     $         CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
            K = DWORK( IR+I-1 )
            IF( K.NE.I )
     $         CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
   10    CONTINUE
C
         DO 20 I = 1, ILO - 1
            K = DWORK( II+I-1 )
            IF( K.NE.I )
     $         CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
            K = DWORK( IR+I-1 )
            IF( K.NE.I )
     $         CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
   20    CONTINUE
C
C        Reduce (A,E) to generalized Hessenberg form and apply the
C        transformations to B and C.
C        Additional workspace: need   N + MAX(N,M);
C                              prefer N + MAX(N,M)*NB.
C
         CALL TG01BD( 'General', 'No Q', 'No Z', N, M, P, ILO, IHI,
     $                DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ), N,
     $                DWORK( IC ), P, DWORK, 1, DWORK, 1, DWORK( IWRK ),
     $                LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Check whether matrix E is nonsingular.
C        Additional workspace: need   3*N.
C
         CALL DTRCON( '1-norm', 'Upper', 'Non Unit', N, DWORK( IE ), N,
     $                RCOND, DWORK( IWRK ), IWORK, IERR )
         IF( RCOND.LE.TEN*DBLE( N )*EPS ) THEN
C
C           Error return: Matrix E is numerically singular.
C
            INFO = 1
            RETURN
         END IF
C
C        Perform QZ algorithm, computing eigenvalues. The generalized
C        Hessenberg form is saved for later use.
C        Additional workspace: need   2*N*N + N;
C                              prefer larger.
C
         IAS  = IWRK
         IES  = IAS + NN
         IWRK = IES + NN
         CALL DLACPY( 'Full', N, N, DWORK( IA ), N, DWORK( IAS ), N )
         CALL DLACPY( 'Full', N, N, DWORK( IE ), N, DWORK( IES ), N )
         CALL DHGEQZ( 'Eigenvalues', 'No Vectors', 'No Vectors', N, ILO,
     $                IHI, DWORK( IAS ), N, DWORK( IES ), N,
     $                DWORK( IR ), DWORK( II ), DWORK( IBT ), DWORK, N,
     $                DWORK, N, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Check if unscaling would cause over/underflow; if so, rescale
C        eigenvalues (DWORK( IR+I-1 ),DWORK( II+I-1 ),DWORK( IBT+I-1 ))
C        so DWORK( IBT+I-1 ) is on the order of E(I,I) and
C        DWORK( IR+I-1 ) and DWORK( II+I-1 ) are on the order of A(I,I).
C
         IF( ILASCL ) THEN
C
            DO 30 I = 1, N
               IF( DWORK( II+I-1 ).NE.ZERO ) THEN
                  IF( ( DWORK( IR+I-1 ) / SAFMAX ).GT.( ANRMTO / ANRM )
     $                                                              .OR.
     $                ( SAFMIN / DWORK( IR+I-1 ) ).GT.( ANRM / ANRMTO )
     $              ) THEN
                     TM = ABS( DWORK( IA+(I-1)*N+I ) / DWORK( IR+I-1 ) )
                     DWORK( IBT+I-1 ) = DWORK( IBT+I-1 )*TM
                     DWORK(  IR+I-1 ) = DWORK(  IR+I-1 )*TM
                     DWORK(  II+I-1 ) = DWORK(  II+I-1 )*TM
                  ELSE IF( ( DWORK( II+I-1 ) / SAFMAX ).GT.
     $                     ( ANRMTO / ANRM ) .OR.
     $               ( SAFMIN / DWORK( II+I-1 ) ).GT.( ANRM / ANRMTO ) )
     $                     THEN
                     TM = ABS( DWORK( IA+I*N+I ) / DWORK( II+I-1 ) )
                     DWORK( IBT+I-1 ) = DWORK( IBT+I-1 )*TM
                     DWORK(  IR+I-1 ) = DWORK(  IR+I-1 )*TM
                     DWORK(  II+I-1 ) = DWORK(  II+I-1 )*TM
                  END IF
               END IF
   30       CONTINUE
C
         END IF
C
         IF( ILESCL ) THEN
C
            DO 40 I = 1, N
               IF( DWORK( II+I-1 ).NE.ZERO ) THEN
                  IF( ( DWORK( IBT+I-1 ) / SAFMAX ).GT.( ENRMTO / ENRM )
     $                                                              .OR.
     $                ( SAFMIN / DWORK( IBT+I-1 ) ).GT.( ENRM / ENRMTO )
     $              ) THEN
                     TM = ABS( DWORK( IE+(I-1)*N+I ) / DWORK( IBT+I-1 ))
                     DWORK( IBT+I-1 ) = DWORK( IBT+I-1 )*TM
                     DWORK(  IR+I-1 ) = DWORK(  IR+I-1 )*TM
                     DWORK(  II+I-1 ) = DWORK(  II+I-1 )*TM
                  END IF
               END IF
   40       CONTINUE
C
         END IF
C
C        Undo scaling.
C
         IF( ILASCL ) THEN
            CALL DLASCL( 'Hessenberg', 0, 0, ANRMTO, ANRM, N, N,
     $                   DWORK( IA ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( IR ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( II ), N, IERR )
         END IF
C
         IF( ILESCL ) THEN
            CALL DLASCL( 'Upper', 0, 0, ENRMTO, ENRM, N, N,
     $                   DWORK( IE ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ENRMTO, ENRM, N, 1,
     $                   DWORK( IBT ), N, IERR )
         END IF
C
      ELSE
C
C        Standard state-space system.
C
         IF( LEQUIL ) THEN
C
C           Equilibrate the system.
C
            MAXRED = HUNDRD
            CALL TB01ID( 'All', N, M, P, MAXRED, DWORK( IA ), N,
     $                   DWORK( IB ), N,  DWORK( IC ), P, DWORK( II ),
     $                   IERR )
         END IF
C
C        For efficiency of later calculations, the system (A,B,C) is
C        reduced to a similar one with the state matrix in Hessenberg
C        form.
C
C        First, permute the matrix A to make it more nearly triangular
C        and apply the permutations to B and C.
C
         CALL DGEBAL( 'Permute', N, DWORK( IA ), N, ILO, IHI,
     $                DWORK( IR ), IERR )
C
         DO 50 I = N, IHI + 1, -1
            K = DWORK( IR+I-1 )
            IF( K.NE.I ) THEN
               CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
               CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
            END IF
   50    CONTINUE
C
         DO 60 I = 1, ILO - 1
            K = DWORK( IR+I-1 )
            IF( K.NE.I ) THEN
               CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
               CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
            END IF
   60    CONTINUE
C
C        Reduce A to upper Hessenberg form and apply the transformations
C        to B and C.
C        Additional workspace: need   N;   (from II)
C                              prefer N*NB.
C
         ITAU = IR
         IWRK = ITAU + N
         CALL DGEHRD( N, ILO, IHI, DWORK( IA ), N, DWORK( ITAU ),
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Additional workspace: need   M;
C                              prefer M*NB.
C
         CALL DORMHR( 'Left', 'Transpose', N, M, ILO, IHI, DWORK( IA ),
     $                N, DWORK( ITAU ), DWORK( IB ), N, DWORK( IWRK ),
     $                LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Additional workspace: need   P;
C                              prefer P*NB.
C
         CALL DORMHR( 'Right', 'NoTranspose', P, N, ILO, IHI,
     $                DWORK( IA ), N, DWORK( ITAU ), DWORK( IC ), P,
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Compute the eigenvalues. The Hessenberg form is saved for
C        later use.
C        Additional workspace:  need   N*N + N;   (from IBT)
C                               prefer larger.
C
         IAS  = IBT
         IWRK = IAS + NN
         CALL DLACPY( 'Full', N, N, DWORK( IA ), N, DWORK( IAS ), N )
         CALL DHSEQR( 'Eigenvalues', 'No Vectors', N, ILO, IHI,
     $                DWORK( IAS ), N, DWORK( IR ), DWORK( II ), DWORK,
     $                N, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         IF( IERR.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
         IF( ILASCL ) THEN
C
C           Undo scaling for the Hessenberg form of A and eigenvalues.
C
            CALL DLASCL( 'Hessenberg', 0, 0, ANRMTO, ANRM, N, N,
     $                   DWORK( IA ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( IR ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( II ), N, IERR )
         END IF
C
      END IF
C
C     Look for (generalized) eigenvalues on the boundary of the
C     stability domain. (Their existence implies an infinite norm.)
C     Additional workspace:  need   2*N.   (from IAS)
C
      IM    = IAS
      IAR   = IM + N
      IMIN  = II
      WRMIN = SAFMAX
      BOUND = EPS*THOUSD
C
      IF ( DISCR ) THEN
         GAMMAL = ZERO
C
C        For discrete-time case, compute the logarithm of the non-zero
C        eigenvalues and save their moduli and absolute real parts.
C        (The logarithms are overwritten on the eigenvalues.)
C        Also, find the minimum distance to the unit circle.
C
         IF ( FULLE ) THEN
C
            DO 70 I = 0, N - 1
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( ( DWORK( IBT+I ).GE.ONE ) .OR.
     $              ( DWORK( IBT+I ).LT.ONE  .AND.
     $                TM.LT.SAFMAX*DWORK( IBT+I ) ) ) THEN
                  TM = TM / DWORK( IBT+I )
               ELSE
C
C                 The pencil has too large eigenvalues. SAFMAX is used.
C
                  TM = SAFMAX
               END IF
               IF ( TM.NE.ZERO ) THEN
                  DWORK( II+I ) = ATAN2( DWORK( II+I ), DWORK( IR+I ) )
                  DWORK( IR+I ) = LOG( TM )
               END IF
               DWORK( IM ) = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               TM = ABS( ONE - TM )
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + I
                  WRMIN = TM
               END IF
               IM = IM + 1
               DWORK( IAR+I ) = ABS( DWORK( IR+I ) )
   70       CONTINUE
C
         ELSE
C
            DO 80 I = 0, N - 1
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( TM.NE.ZERO ) THEN
                  DWORK( II+I ) = ATAN2( DWORK( II+I ), DWORK( IR+I ) )
                  DWORK( IR+I ) = LOG( TM )
               END IF
               DWORK( IM ) = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               TM = ABS( ONE - TM )
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + I
                  WRMIN = TM
               END IF
               IM = IM + 1
               DWORK( IAR+I ) = ABS( DWORK( IR+I ) )
   80       CONTINUE
C
         END IF
C
      ELSE
C
C        For continuous-time case, save moduli of eigenvalues and
C        absolute real parts and find the maximum modulus and minimum
C        absolute real part.
C
         WMAX = ZERO
C
         IF ( FULLE ) THEN
C
            DO 90 I = 0, N - 1
               TM = ABS( DWORK( IR+I ) )
               DWORK( IM ) = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( ( DWORK( IBT+I ).GE.ONE ) .OR.
     $              ( DWORK( IBT+I ).LT.ONE  .AND.
     $                DWORK( IM ).LT.SAFMAX*DWORK( IBT+I ) ) )
     $                   THEN
                  TM = TM / DWORK( IBT+I )
                  DWORK( IM ) = DWORK( IM ) / DWORK( IBT+I )
               ELSE
                  IF ( TM.LT.SAFMAX*DWORK( IBT+I ) ) THEN
                     TM = TM / DWORK( IBT+I )
                  ELSE
C
C                    The pencil has too large eigenvalues.
C                    SAFMAX is used.
C
                     TM = SAFMAX
                  END IF
                  DWORK( IM ) = SAFMAX
               END IF
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + I
                  WRMIN = TM
               END IF
               DWORK( IAR+I ) = TM
               IF( DWORK( IM ).GT.WMAX )
     $            WMAX = DWORK( IM )
               IM = IM + 1
   90       CONTINUE
C
         ELSE
C
            DO 100 I = 0, N - 1
               TM = ABS( DWORK( IR+I ) )
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + I
                  WRMIN = TM
               END IF
               DWORK( IM ) = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF( DWORK( IM ).GT.WMAX )
     $            WMAX = DWORK( IM )
               IM = IM + 1
               DWORK( IAR+I ) = TM
  100       CONTINUE
C
         END IF
C
         BOUND = BOUND + EPS*WMAX
C
      END IF
C
      IM = IM - N
C
      IF( WRMIN.LT.BOUND ) THEN
C
C        The L-infinity norm was found as infinite.
C
         GPEAK( 1 ) = ONE
         GPEAK( 2 ) = ZERO
         TM = ABS( DWORK( IMIN ) )
         IF ( DISCR )
     $      TM = ABS( ATAN2( SIN( TM ), COS( TM ) ) )
         FPEAK( 1 ) = TM
         IF ( TM.LT.SAFMAX ) THEN
            FPEAK( 2 ) = ONE
         ELSE
            FPEAK( 2 ) = ZERO
         END IF
C
         DWORK( 1 ) = MAXWRK
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine the maximum singular value of
C        G(lambda) = C*inv(lambda*E - A)*B + D,
C     over a selected set of frequencies. Besides the frequencies w = 0,
C     w = pi (if DICO = 'D'), and the given value FPEAK, this test set
C     contains the peak frequency for each mode (or an approximation
C     of it). The (generalized) Hessenberg form of the system is used.
C
C     First, determine the maximum singular value of G(0) and set FPEAK
C     accordingly.
C     Additional workspace:
C           complex: need   1, if DICO = 'C';
C                           (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)), otherwise;
C                    prefer larger;
C           real:    need   LDW0+LDW1+LDW2, where
C                           LDW0 = N*N+N*M, if DICO = 'C';
C                           LDW0 = 0,       if DICO = 'D';
C                           LDW1 = P*M,     if DICO = 'C', JOBD = 'Z';
C                           LDW1 = 0,       otherwise;
C                           LDW2 = MIN(P,M)+MAX(3*MIN(P,M)+MAX(P,M),
C                                               5*MIN(P,M)),
C                                              if DICO = 'C';
C                           LDW2 = 6*MIN(P,M), otherwise.
C                    prefer larger.
C
      IF ( DISCR ) THEN
         IAS  = IA
         IBS  = IB
         IWRK = IAR + N
      ELSE
         IAS  = IAR + N
         IBS  = IAS + NN
         IWRK = IBS + N*M
         CALL DLACPY( 'Upper', N, N, DWORK( IA ), N, DWORK( IAS ), N )
         CALL DCOPY(  N-1, DWORK( IA+1 ), N+1, DWORK( IAS+1 ), N+1 )
         CALL DLACPY( 'Full', N, M, DWORK( IB ), N, DWORK( IBS ), N )
      END IF
      GAMMA = AB13DX( DICO, JOBE, JOBD, N, M, P, ZERO, DWORK( IAS ), N,
     $                DWORK( IE ), N, DWORK( IBS ), N, DWORK( IC ), P,
     $                DWORK( ID ), P, IWORK, DWORK( IWRK ),
     $                LDWORK-IWRK+1, CWORK, LCWORK, IERR )
      MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
      IF( IERR.GE.1 .AND. IERR.LE.N ) THEN
         GPEAK( 1 ) = ONE
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ZERO
         FPEAK( 2 ) = ONE
         GO TO 340
      ELSE IF( IERR.EQ.N+1 ) THEN
         INFO = 3
         RETURN
      END IF
C
      FPEAKS = FPEAK( 1 )
      FPEAKI = FPEAK( 2 )
      IF( GAMMAL.LT.GAMMA ) THEN
         GAMMAL     = GAMMA
         FPEAK( 1 ) = ZERO
         FPEAK( 2 ) = ONE
      ELSE IF( .NOT.DISCR ) THEN
         FPEAK( 1 ) = ONE
         FPEAK( 2 ) = ZERO
      END IF
C
      MAXCWK = INT( CWORK( 1 ) )
C
      IF( DISCR ) THEN
C
C        Try the frequency w = pi.
C
         PI    = FOUR*ATAN( ONE )
         GAMMA = AB13DX( DICO, JOBE, JOBD, N, M, P, PI, DWORK( IA ),
     $                   N, DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ),
     $                   P, DWORK( ID ), P, IWORK, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, CWORK, LCWORK, IERR )
         MAXCWK = MAX( INT( CWORK( 1 ) ), MAXCWK )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         IF( IERR.GE.1 .AND. IERR.LE.N ) THEN
            GPEAK( 1 ) = ONE
            FPEAK( 1 ) = PI
            GPEAK( 2 ) = ZERO
            FPEAK( 2 ) = ONE
            GO TO 340
         ELSE IF( IERR.EQ.N+1 ) THEN
            INFO = 3
            RETURN
         END IF
C
         IF( GAMMAL.LT.GAMMA ) THEN
            GAMMAL     = GAMMA
            FPEAK( 1 ) = PI
            FPEAK( 2 ) = ONE
         END IF
C
      ELSE
         IWRK = IAS
C
C        Restore D, if needed.
C
         IF ( WITHD )
     $      CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
      END IF
C
C     Build the remaining set of frequencies.
C     Complex workspace:  need   (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M));
C                         prefer larger.
C     Real workspace:     need   LDW2, see above;
C                         prefer larger.
C
      IF ( MIN( FPEAKS, FPEAKI ).NE.ZERO ) THEN
C
C        Compute also the norm at the given (finite) frequency.
C
         GAMMA = AB13DX( DICO, JOBE, JOBD, N, M, P, FPEAKS, DWORK( IA ),
     $                   N, DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ),
     $                   P, DWORK( ID ), P, IWORK, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, CWORK, LCWORK, IERR )
         MAXCWK = MAX( INT( CWORK( 1 ) ), MAXCWK )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         IF ( DISCR ) THEN
            TM = ABS( ATAN2( SIN( FPEAKS ), COS( FPEAKS ) ) )
         ELSE
            TM = FPEAKS
         END IF
         IF( IERR.GE.1 .AND. IERR.LE.N ) THEN
            GPEAK( 1 ) = ONE
            FPEAK( 1 ) = TM
            GPEAK( 2 ) = ZERO
            FPEAK( 2 ) = ONE
            GO TO 340
         ELSE IF( IERR.EQ.N+1 ) THEN
            INFO = 3
            RETURN
         END IF
C
         IF( GAMMAL.LT.GAMMA ) THEN
            GAMMAL     = GAMMA
            FPEAK( 1 ) = TM
            FPEAK( 2 ) = ONE
         END IF
C
      END IF
C
      DO 110 I = 0, N - 1
         IF( DWORK( II+I ).GE.ZERO .AND. DWORK( IM+I ).GT.ZERO ) THEN
            IF ( ( DWORK( IM+I ).GE.ONE ) .OR. ( DWORK( IM+I ).LT.ONE
     $            .AND. DWORK( IAR+I ).LT.SAFMAX*DWORK( IM+I ) ) ) THEN
               RAT = DWORK( IAR+I ) / DWORK( IM+I )
            ELSE
               RAT = ONE
            END IF
            OMEGA = DWORK( IM+I )*SQRT( MAX( P25, ONE - TWO*RAT**2 ) )
C
            GAMMA = AB13DX( DICO, JOBE, JOBD, N, M, P, OMEGA,
     $                      DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ),
     $                      N, DWORK( IC ), P, DWORK( ID ), P, IWORK,
     $                      DWORK( IWRK ), LDWORK-IWRK+1, CWORK, LCWORK,
     $                      IERR )
            MAXCWK = MAX( INT( CWORK( 1 ) ), MAXCWK )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            IF ( DISCR ) THEN
               TM = ABS( ATAN2( SIN( OMEGA ), COS( OMEGA ) ) )
            ELSE
               TM = OMEGA
            END IF
            IF( IERR.GE.1 .AND. IERR.LE.N ) THEN
               GPEAK( 1 ) = ONE
               FPEAK( 1 ) = TM
               GPEAK( 2 ) = ZERO
               FPEAK( 2 ) = ONE
               GO TO 340
            ELSE IF( IERR.EQ.N+1 ) THEN
               INFO = 3
               RETURN
            END IF
C
            IF( GAMMAL.LT.GAMMA ) THEN
               GAMMAL     = GAMMA
               FPEAK( 1 ) = TM
               FPEAK( 2 ) = ONE
            END IF
C
         END IF
  110 CONTINUE
C
C     Return if the lower bound is zero.
C
      IF( GAMMAL.EQ.ZERO ) THEN
         GPEAK( 1 ) = ZERO
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         GO TO 340
      END IF
C
C     Start the modified gamma iteration for the Bruinsma-Steinbuch
C     algorithm.
C
      IF ( .NOT.DISCR )
     $   RTOL = HUNDRD*TOLER
      ITER = 0
C
C     WHILE ( Iteration may continue ) DO
C
  120 CONTINUE
C
         ITER   = ITER + 1
         GAMMA  = ( ONE + TOL )*GAMMAL
         USEPEN = FULLE .OR. DISCR
         IF ( .NOT.USEPEN .AND. WITHD ) THEN
C
C           Check whether one can use an explicit Hamiltonian matrix:
C           compute
C           min(rcond(GAMMA**2*Im - S'*S), rcond(GAMMA**2*Ip - S*S')).
C           If P = M = 1, then GAMMA**2 - S(1)**2 is used instead.
C
            IF ( M.NE.P ) THEN
               RCOND = ONE - ( DWORK( IS ) / GAMMA )**2
            ELSE IF ( MINPM.GT.1 ) THEN
               RCOND = ( GAMMA**2 - DWORK( IS )**2 ) /
     $                 ( GAMMA**2 - DWORK( IS+P-1 )**2 )
            ELSE
               RCOND = GAMMA**2 - DWORK( IS )**2
            END IF
C
            USEPEN = RCOND.LT.RTOL
         END IF
C
         IF ( USEPEN ) THEN
C
C           Use the QZ algorithm on a pencil.
C           Additional workspace here:  need   6*N.   (from IR)
C
            II   = IR  + N2
            IBT  = II  + N2
            IH12 = IBT + N2
            IM   = IH12
C
C           Set up the needed parts of the Hamiltonian pencil (H,J),
C
C                  ( H11  H12 )
C              H = (          ) ,
C                  ( H21  H22 )
C
C           with
C
C                 ( A  0  )            ( 0  B )            ( E  0  )
C           H11 = (       ),     H12 = (      )/nB,  J11 = (       ),
C                 ( 0 -A' )            ( C' 0 )            ( 0  E' )
C
C                 ( C  0  )            ( Ip  D/g )
C           H21 = (       )*nB,  H22 = (         ),
C                 ( 0 -B' )            ( D'/g Im )
C
C           if DICO = 'C', and
C
C                 ( A  0  )            ( B  0  )            ( E  0 )
C           H11 = (       ),     H12 = (       )/nB,  J11 = (      ),
C                 ( 0  E' )            ( 0  C' )            ( 0  A')
C
C                 ( 0  0  )            ( Im  D'/g )         ( 0  B')
C           H21 = (       )*nB,  H22 = (          ),  J21 = (      )*nB,
C                 ( C  0  )            ( D/g  Ip  )         ( 0  0 )
C
C           if DICO = 'D', where g = GAMMA, and nB = norm(B,1).
C           First build [H12; H22].
C
            TEMP( 1 ) = ZERO
            IH = IH12
C
            IF ( DISCR ) THEN
C
               DO 150 J = 1, M
C
                  DO 130 I = 1, N
                     DWORK( IH ) = B( I, J ) / BNORM
                     IH = IH + 1
  130             CONTINUE
C
                  CALL DCOPY( N+M, TEMP, 0, DWORK( IH ), 1 )
                  DWORK( IH+N+J-1 ) = ONE
                  IH = IH + N + M
C
                  DO 140 I = 1, P
                     DWORK( IH ) = D( I, J ) / GAMMA
                     IH = IH + 1
  140             CONTINUE
C
  150          CONTINUE
C
               DO 180 J = 1, P
                  CALL DCOPY( N, TEMP, 0, DWORK( IH ), 1 )
                  IH = IH + N
C
                  DO 160 I = 1, N
                     DWORK( IH ) = C( J, I ) / BNORM
                     IH = IH + 1
  160             CONTINUE
C
                  DO 170 I = 1, M
                     DWORK( IH ) = D( J, I ) / GAMMA
                     IH = IH + 1
  170             CONTINUE
C
                  CALL DCOPY( P, TEMP, 0, DWORK( IH ), 1 )
                  DWORK( IH+J-1 ) = ONE
                  IH = IH + P
  180          CONTINUE
C
            ELSE
C
               DO 210 J = 1, P
                  CALL DCOPY( N, TEMP, 0, DWORK( IH ), 1 )
                  IH = IH + N
C
                  DO 190 I = 1, N
                     DWORK( IH ) = C( J, I ) / BNORM
                     IH = IH + 1
  190             CONTINUE
C
                  CALL DCOPY( P, TEMP, 0, DWORK( IH ), 1 )
                  DWORK( IH+J-1 ) = ONE
                  IH = IH + P
C
                  DO 200 I = 1, M
                     DWORK( IH ) = D( J, I ) / GAMMA
                     IH = IH + 1
  200             CONTINUE
C
  210          CONTINUE
C
               DO 240 J = 1, M
C
                  DO 220 I = 1, N
                     DWORK( IH ) = B( I, J ) / BNORM
                     IH = IH + 1
  220             CONTINUE
C
                  CALL DCOPY( N, TEMP, 0, DWORK( IH ), 1 )
                  IH = IH + N
C
                  DO 230 I = 1, P
                     DWORK( IH ) = D( I, J ) / GAMMA
                     IH = IH + 1
  230             CONTINUE
C
                  CALL DCOPY( M, TEMP, 0, DWORK( IH ), 1 )
                  DWORK( IH+J-1 ) = ONE
                  IH = IH + M
  240          CONTINUE
C
            END IF
C
C           Compute the QR factorization of [H12; H22].
C           For large P and M, it could be more efficient to exploit the
C           structure of [H12; H22] and use the factored form of Q.
C           Additional workspace: need   (2*N+P+M)*(2*N+P+M)+2*(P+M);
C                                 prefer (2*N+P+M)*(2*N+P+M)+P+M+
C                                                           (P+M)*NB.
C
            ITAU = IH12 + N2PM*N2PM
            IWRK = ITAU + PM
            CALL DGEQRF( N2PM, PM, DWORK( IH12 ), N2PM, DWORK( ITAU ),
     $                   DWORK( IWRK ), LDWORK-IWRK+1, IERR )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C           Apply part of the orthogonal transformation:
C           Q1 = Q(:,P+M+(1:2*N))' to the matrix [H11; H21/GAMMA].
C           If DICO = 'C', apply Q(1:2*N,P+M+(1:2*N))' to the
C           matrix J11.
C           If DICO = 'D', apply Q1 to the matrix [J11; J21/GAMMA].
C           H11, H21, J11, and J21 are not fully built.
C           First, build the (2*N+P+M)-by-(2*N+P+M) matrix Q.
C           Using Q will often provide better efficiency than the direct
C           use of the factored form of Q, especially when P+M < N.
C           Additional workspace: need   P+M+2*N+P+M;
C                                 prefer P+M+(2*N+P+M)*NB.
C
            CALL DORGQR( N2PM, N2PM, PM, DWORK( IH12 ), N2PM,
     $                   DWORK( ITAU ), DWORK( IWRK ), LDWORK-IWRK+1,
     $                   IERR )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C           Additional workspace: need   8*N*N.
C
            IPA  = ITAU
            IPE  = IPA + 4*NN
            IWRK = IPE + 4*NN
            CALL DGEMM( 'Transpose', 'No Transpose', N2, N, N, ONE,
     $                  DWORK( IH12+PM*N2PM ), N2PM, A, LDA, ZERO,
     $                  DWORK( IPA ), N2 )
            IF ( DISCR ) THEN
               CALL DGEMM( 'Transpose', 'No Transpose', N2, N, P,
     $                     BNORM/GAMMA, DWORK( IH12+PM*N2PM+N2+M), N2PM,
     $                     C, LDC, ONE, DWORK( IPA ), N2 )
               IF ( FULLE ) THEN
                  CALL DGEMM( 'Transpose', 'Transpose', N2, N, N, ONE,
     $                        DWORK( IH12+PM*N2PM+N ), N2PM, E, LDE,
     $                        ZERO, DWORK( IPA+2*NN ), N2 )
               ELSE
                  CALL MA02AD( 'Full', N, N2, DWORK( IH12+PM*N2PM+N ),
     $                         N2PM, DWORK( IPA+2*NN ), N2 )
                  NY = N
               END IF
            ELSE
               CALL DGEMM( 'Transpose', 'No Transpose', N2, N, P,
     $                     BNORM/GAMMA, DWORK( IH12+PM*N2PM+N2), N2PM,
     $                     C, LDC, ONE, DWORK( IPA ), N2 )
               CALL DGEMM( 'Transpose', 'Transpose', N2, N, N, -ONE,
     $                     DWORK( IH12+PM*N2PM+N ), N2PM, A, LDA, ZERO,
     $                     DWORK( IPA+2*NN ), N2 )
               CALL DGEMM( 'Transpose', 'Transpose', N2, N, M,
     $                     -BNORM/GAMMA, DWORK( IH12+PM*N2PM+N2+P),
     $                     N2PM, B, LDB, ONE, DWORK( IPA+2*NN ), N2 )
               NY = N2
            END IF
C
            IF ( FULLE ) THEN
               CALL DGEMM( 'Transpose', 'No Transpose', N2, N, N, ONE,
     $                     DWORK( IH12+PM*N2PM ), N2PM, E, LDE, ZERO,
     $                     DWORK( IPE ), N2 )
            ELSE
               CALL MA02AD( 'Full', NY, N2, DWORK( IH12+PM*N2PM ),
     $                      N2PM, DWORK( IPE ), N2 )
            END IF
            IF ( DISCR ) THEN
               CALL DGEMM( 'Transpose', 'Transpose', N2, N, N, ONE,
     $                     DWORK( IH12+PM*N2PM+N ), N2PM, A, LDA,
     $                     ZERO, DWORK( IPE+2*NN ), N2 )
               CALL DGEMM( 'Transpose', 'Transpose', N2, N, M,
     $                     BNORM/GAMMA, DWORK( IH12+PM*N2PM+N2 ), N2PM,
     $                     B, LDB, ONE, DWORK( IPE+2*NN ), N2 )
            ELSE
               IF ( FULLE )
     $            CALL DGEMM( 'Transpose', 'Transpose', N2, N, N, ONE,
     $                        DWORK( IH12+PM*N2PM+N ), N2PM, E, LDE,
     $                        ZERO, DWORK( IPE+2*NN ), N2 )
            END IF
C
C           Compute the eigenvalues of the Hamiltonian pencil.
C           Additional workspace: need   16*N;
C                                 prefer larger.
C
            CALL DGGEV( 'No Vectors', 'No Vectors', N2, DWORK( IPA ),
     $                  N2, DWORK( IPE ), N2, DWORK( IR ), DWORK( II ),
     $                  DWORK( IBT ), DWORK, N2, DWORK, N2,
     $                  DWORK( IWRK ), LDWORK-IWRK+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 2
               RETURN
            END IF
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
         ELSE IF ( .NOT.WITHD ) THEN
C
C           Standard continuous-time case with D = 0.
C           Form the needed part of the Hamiltonian matrix explicitly:
C              H = H11 - H12*inv(H22)*H21/g.
C           Additional workspace: need   2*N*N+N.   (from IBT)
C
            IH   = IBT
            IH12 = IH   + NN
            ISL  = IH12 + NN + N
            CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IH ), N )
C
C           Compute triangles of -C'*C/GAMMA and B*B'/GAMMA.
C
            CALL DSYRK( 'Lower', 'Transpose', N, P, -ONE/GAMMA, C, LDC,
     $                  ZERO, DWORK( IH12 ), N )
            CALL DSYRK( 'Upper', 'No Transpose', N, M, ONE/GAMMA, B,
     $                  LDB, ZERO, DWORK( IH12+N ), N )
C
         ELSE
C
C           Standard continuous-time case with D <> 0 and the SVD of D
C           can be used. Compute explicitly the needed part of the
C           Hamiltonian matrix:
C
C               (A+B1*S'*inv(g^2*Ip-S*S')*C1' g*B1*inv(g^2*Im-S'*S)*B1')
C           H = (                                                      )
C               (  -g*C1*inv(g^2*Ip-S*S')*C1'            -H11'         )
C
C           where g = GAMMA, B1 = B*V, C1 = C'*U, and H11 is the first
C           block of H.
C           Primary additional workspace: need   2*N*N+N   (from IBT)
C           (for building the relevant part of the Hamiltonian matrix).
C
C           Compute C1*sqrt(inv(g^2*Ip-S*S')) .
C           Additional workspace: need   MAX(M,P)+N*P.
C
            IH   = IBT
            IH12 = IH   + NN
            ISL  = IH12 + NN + N
C
            DO 250 I = 0, MINPM - 1
               DWORK( ISL+I ) = ONE/SQRT( GAMMA**2 - DWORK( IS+I )**2 )
  250       CONTINUE
C
            IF ( M.LT.P ) THEN
               DWORK( ISL+M ) = ONE / GAMMA
               CALL DCOPY( P-M-1, DWORK( ISL+M ), 0, DWORK( ISL+M+1 ),
     $                     1 )
            END IF
            ISC = ISL + MAX( M, P )
            CALL DLACPY( 'Full', N, P, DWORK( ICU ), N, DWORK( ISC ),
     $                   N )
            CALL MB01SD( 'Column', N, P, DWORK( ISC ), N, DWORK,
     $                   DWORK( ISL ) )
C
C           Compute B1*S' .
C           Additional workspace: need   N*M.
C
            ISB = ISC + P*N
            CALL DLACPY( 'Full', N, M, DWORK( IBV ), N, DWORK( ISB ),
     $                   N )
            CALL MB01SD( 'Column', N, MINPM, DWORK( ISB ), N, DWORK,
     $                   DWORK( IS ) )
C
C           Compute B1*S'*sqrt(inv(g^2*Ip-S*S')) .
C
            CALL MB01SD( 'Column', N, MINPM, DWORK( ISB ), N, DWORK,
     $                   DWORK( ISL ) )
C
C           Compute H11 .
C
            CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IH ), N )
            CALL DGEMM( 'No Transpose', 'Transpose', N, N, MINPM, ONE,
     $                  DWORK( ISB ), N, DWORK( ISC ), N, ONE,
     $                  DWORK( IH ), N )
C
C           Compute B1*sqrt(inv(g^2*Im-S'*S)) .
C
            IF ( P.LT.M ) THEN
               DWORK( ISL+P ) = ONE / GAMMA
               CALL DCOPY( M-P-1, DWORK( ISL+P ), 0, DWORK( ISL+P+1 ),
     $                     1 )
            END IF
            CALL DLACPY( 'Full', N, M, DWORK( IBV ), N, DWORK( ISB ),
     $                   N )
            CALL MB01SD( 'Column', N, M, DWORK( ISB ), N, DWORK,
     $                   DWORK( ISL ) )
C
C           Compute the lower triangle of H21 and the upper triangle
C           of H12.
C
            CALL DSYRK( 'Lower', 'No Transpose', N, P, -GAMMA,
     $                  DWORK( ISC ), N, ZERO, DWORK( IH12 ), N )
            CALL DSYRK( 'Upper', 'No Transpose', N, M, GAMMA,
     $                  DWORK( ISB ), N, ZERO, DWORK( IH12+N ), N )
         END IF
C
         IF ( .NOT.USEPEN ) THEN
C
C           Compute the eigenvalues of the Hamiltonian matrix by the
C           symplectic URV and the periodic Schur decompositions.
C           Additional workspace: need   (2*N+8)*N;
C                                 prefer larger.
C
            IWRK = ISL + NN
            CALL MB03XD( 'Both',  'Eigenvalues', 'No vectors',
     $                   'No vectors', N, DWORK( IH ), N, DWORK( IH12 ),
     $                   N, DWORK( ISL ), N, TEMP, 1, TEMP, 1, TEMP, 1,
     $                   TEMP, 1, DWORK( IR ), DWORK( II ), ILO,
     $                   DWORK( IWRK ), DWORK( IWRK+N ),
     $                   LDWORK-IWRK-N+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 2
               RETURN
            END IF
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK + N - 1, MAXWRK )
         END IF
C
C        Detect eigenvalues on the boundary of the stability domain,
C        if any. The test is based on a round-off level of eps*rho(H)
C        (after balancing) resulting in worst-case perturbations of
C        order sqrt(eps*rho(H)), for continuous-time systems, on the
C        real part of poles of multiplicity two (typical as GAMMA
C        approaches the infinity norm). Similarly, in the discrete-time
C        case. Above, rho(H) is the maximum modulus of eigenvalues
C        (continuous-time case).
C
C        Compute maximum eigenvalue modulus and check the absolute real
C        parts (if DICO = 'C'), or moduli (if DICO = 'D').
C
         WMAX = ZERO
C
         IF ( USEPEN ) THEN
C
C           Additional workspace: need   2*N, if DICO = 'D';   (from IM)
C                                        0,   if DICO = 'C'.
C
            DO 260 I = 0, N2 - 1
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( ( DWORK( IBT+I ).GE.ONE ) .OR.
     $              ( DWORK( IBT+I ).LT.ONE  .AND.
     $               TM.LT.SAFMAX*DWORK( IBT+I ) ) ) THEN
                  TM = TM / DWORK( IBT+I )
               ELSE
C
C                 The pencil has too large eigenvalues. SAFMAX is used.
C
                  TM = SAFMAX
               END IF
               WMAX = MAX( WMAX, TM )
               IF ( DISCR )
     $            DWORK( IM+I ) = TM
  260       CONTINUE
C
         ELSE
C
            DO 270 I = 0, N - 1
               TM   = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               WMAX = MAX( WMAX, TM )
  270       CONTINUE
C
         END IF
C
         NEI = 0
C
         IF ( USEPEN ) THEN
C
            DO 280 I = 0, N2 - 1
               IF ( DISCR ) THEN
                  TM = ABS( ONE - DWORK( IM+I ) )
               ELSE
                  TM = ABS( DWORK( IR+I ) )
                  IF ( ( DWORK( IBT+I ).GE.ONE ) .OR.
     $                 ( DWORK( IBT+I ).LT.ONE  .AND.
     $                  TM.LT.SAFMAX*DWORK( IBT+I ) ) ) THEN
                     TM = TM / DWORK( IBT+I )
                  ELSE
C
C                    The pencil has too large eigenvalues.
C                    SAFMAX is used.
C
                     TM = SAFMAX
                  END IF
               END IF
               IF ( TM.LE.TOLER*SQRT( HUNDRD + WMAX ) ) THEN
                  DWORK( IR+NEI ) = DWORK( IR+I ) / DWORK( IBT+I )
                  DWORK( II+NEI ) = DWORK( II+I ) / DWORK( IBT+I )
                  NEI = NEI + 1
               END IF
  280       CONTINUE
C
         ELSE
C
            DO 290 I = 0, N - 1
               TM = ABS( DWORK( IR+I ) )
               IF ( TM.LE.TOLER*SQRT( HUNDRD + WMAX ) ) THEN
                  DWORK( IR+NEI ) = DWORK( IR+I )
                  DWORK( II+NEI ) = DWORK( II+I )
                  NEI = NEI + 1
               END IF
  290       CONTINUE
C
         END IF
C
         IF( NEI.EQ.0 ) THEN
C
C           There is no eigenvalue on the boundary of the stability
C           domain for G = ( ONE + TOL )*GAMMAL. The norm was found.
C
            GPEAK( 1 ) = GAMMAL
            GPEAK( 2 ) = ONE
            GO TO 340
         END IF
C
C        Compute the frequencies where the gain G is attained and
C        generate new test frequencies.
C
         NWS = 0
C
         IF ( DISCR ) THEN
C
            DO 300 I = 0, NEI - 1
               TM = ATAN2( DWORK( II+I ), DWORK( IR+I ) )
               DWORK( IR+I ) = MAX( EPS, TM )
               NWS = NWS + 1
  300       CONTINUE
C
         ELSE
C
            J = 0
C
            DO 310 I = 0, NEI - 1
               IF ( DWORK( II+I ).GT.EPS ) THEN
                  DWORK( IR+NWS ) = DWORK( II+I )
                  NWS = NWS + 1
               ELSE IF ( DWORK( II+I ).EQ.EPS ) THEN
                  J = J + 1
                  IF ( J.EQ.1 ) THEN
                     DWORK( IR+NWS ) = EPS
                     NWS = NWS + 1
                  END IF
               END IF
  310       CONTINUE
C
         END IF
C
         CALL DLASRT( 'Increasing', NWS, DWORK( IR ), IERR )
         LW = 1
C
         DO 320 I = 0, NWS - 1
            IF ( DWORK( IR+LW-1 ).NE.DWORK( IR+I ) ) THEN
               DWORK( IR+LW ) = DWORK( IR+I )
               LW = LW + 1
            END IF
  320    CONTINUE
C
         IF ( LW.EQ.1 ) THEN
            IF ( ITER.EQ.1 .AND. NWS.GE.1 ) THEN
C
C              Duplicate the frequency trying to force iteration.
C
               DWORK( IR+1 ) = DWORK( IR )
               LW = LW + 1
            ELSE
C
C              The norm was found.
C
               GPEAK( 1 ) = GAMMAL
               GPEAK( 2 ) = ONE
               GO TO 340
            END IF
         END IF
C
C        Form the vector of mid-points and compute the gain at new test
C        frequencies. Save the current lower bound.
C
         IWRK   = IR + LW
         GAMMAS = GAMMAL
C
         DO 330 I = 0, LW - 2
            IF ( DISCR ) THEN
               OMEGA = ( DWORK( IR+I ) + DWORK( IR+I+1 ) ) / TWO
            ELSE
               OMEGA = SQRT( DWORK( IR+I )*DWORK( IR+I+1 ) )
            END IF
C
C           Additional workspace:  need   LDW2, see above;
C                                  prefer larger.
C
            GAMMA = AB13DX( DICO, JOBE, JOBD, N, M, P, OMEGA,
     $                      DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ),
     $                      N, DWORK( IC ), P, DWORK( ID ), P, IWORK,
     $                      DWORK( IWRK ), LDWORK-IWRK+1, CWORK, LCWORK,
     $                      IERR )
            MAXCWK = MAX( INT( CWORK( 1 ) ), MAXCWK )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            IF ( DISCR ) THEN
               TM = ABS( ATAN2( SIN( OMEGA ), COS( OMEGA ) ) )
            ELSE
               TM = OMEGA
            END IF
            IF( IERR.GE.1 .AND. IERR.LE.N ) THEN
               GPEAK( 1 ) = ONE
               FPEAK( 1 ) = TM
               GPEAK( 2 ) = ZERO
               FPEAK( 2 ) = ONE
               GO TO 340
            ELSE IF( IERR.EQ.N+1 ) THEN
               INFO = 3
               RETURN
            END IF
C
            IF( GAMMAL.LT.GAMMA ) THEN
               GAMMAL     = GAMMA
               FPEAK( 1 ) = TM
               FPEAK( 2 ) = ONE
            END IF
  330    CONTINUE
C
C        If the lower bound has not been improved, return. (This is a
C        safeguard against undetected modes of Hamiltonian matrix on the
C        boundary of the stability domain.)
C
         IF ( GAMMAL.LT.GAMMAS*( ONE + TOL/TEN ) ) THEN
            GPEAK( 1 ) = GAMMAL
            GPEAK( 2 ) = ONE
            GO TO 340
         END IF
C
C     END WHILE
C
      IF ( ITER.LE.MAXIT ) THEN
         GO TO 120
      ELSE
         INFO = 4
         RETURN
      END IF
C
  340 CONTINUE
      DWORK( 1 ) = MAXWRK
      CWORK( 1 ) = MAXCWK
      RETURN
C *** Last line of AB13DD ***
      END
