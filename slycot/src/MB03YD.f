      SUBROUTINE MB03YD( WANTT, WANTQ, WANTZ, N, ILO, IHI, ILOQ, IHIQ,
     $                   A, LDA, B, LDB, Q, LDQ, Z, LDZ, ALPHAR, ALPHAI,
     $                   BETA, DWORK, LDWORK, INFO )
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
C     To deal with small subtasks of the product eigenvalue problem.
C
C     MB03YD is an auxiliary routine called by SLICOT Library routine
C     MB03XP.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     WANTT   LOGICAL
C             Indicates whether the user wishes to compute the full
C             Schur form or the eigenvalues only, as follows:
C             = .TRUE. :  Compute the full Schur form;
C             = .FALSE.:  compute the eigenvalues only.
C
C     WANTQ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Q as follows:
C             = .TRUE. :  The matrix Q is updated;
C             = .FALSE.:  the matrix Q is not required.
C
C     WANTZ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Z as follows:
C             = .TRUE. :  The matrix Z is updated;
C             = .FALSE.:  the matrix Z is not required.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and B. N >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that the matrices A and B are already
C             (quasi) upper triangular in rows and columns 1:ILO-1 and
C             IHI+1:N. The routine works primarily with the submatrices
C             in rows and columns ILO to IHI, but applies the
C             transformations to all the rows and columns of the
C             matrices A and B, if WANTT = .TRUE..
C             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N.
C
C     ILOQ    (input) INTEGER
C     IHIQ    (input) INTEGER
C             Specify the rows of Q and Z to which transformations
C             must be applied if WANTQ = .TRUE. and WANTZ = .TRUE.,
C             respectively.
C             1 <= ILOQ <= ILO; IHI <= IHIQ <= N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper Hessenberg matrix A.
C             On exit, if WANTT = .TRUE., the leading N-by-N part of
C             this array is upper quasi-triangular in rows and columns
C             ILO:IHI.
C             If WANTT = .FALSE., the diagonal elements and 2-by-2
C             diagonal blocks of A will be correct, but the remaining
C             parts of A are unspecified on exit.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular matrix B.
C             On exit, if WANTT = .TRUE., the leading N-by-N part of
C             this array contains the transformed upper triangular
C             matrix. 2-by-2 blocks in B corresponding to 2-by-2 blocks
C             in A will be reduced to positive diagonal form. (I.e., if
C             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j)
C             and B(j+1,j+1) will be positive.)
C             If WANTT = .FALSE., the elements corresponding to diagonal
C             elements and 2-by-2 diagonal blocks in A will be correct,
C             but the remaining parts of B are unspecified on exit.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, if WANTQ = .TRUE., then the leading N-by-N part
C             of this array must contain the current matrix Q of
C             transformations accumulated by MB03XP.
C             On exit, if WANTQ = .TRUE., then the leading N-by-N part
C             of this array contains the matrix Q updated in the
C             submatrix Q(ILOQ:IHIQ,ILO:IHI).
C             If WANTQ = .FALSE., Q is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= 1.
C             If WANTQ = .TRUE., LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, if WANTZ = .TRUE., then the leading N-by-N part
C             of this array must contain the current matrix Z of
C             transformations accumulated by MB03XP.
C             On exit, if WANTZ = .TRUE., then the leading N-by-N part
C             of this array contains the matrix Z updated in the
C             submatrix Z(ILOQ:IHIQ,ILO:IHI).
C             If WANTZ = .FALSE., Z is not referenced.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= 1.
C             If WANTZ = .TRUE., LDZ >= MAX(1,N).
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C             The i-th (ILO <= i <= IHI) computed eigenvalue is given
C             by BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two
C             eigenvalues are computed as a complex conjugate pair,
C             they are stored in consecutive elements of ALPHAR, ALPHAI
C             and BETA. If WANTT = .TRUE., the eigenvalues are stored in
C             the same order as on the diagonals of the Schur forms of
C             A and B.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if  INFO = -19,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, then MB03YD failed to compute the Schur
C                   form in a total of 30*(IHI-ILO+1) iterations;
C                   elements i+1:n of ALPHAR, ALPHAI and BETA contain
C                   successfully computed eigenvalues.
C
C     METHOD
C
C     The implemented algorithm is a double-shift version of the
C     periodic QR algorithm described in [1,3] with some minor
C     modifications [2]. The eigenvalues are computed via an implicit
C     complex single shift algorithm.
C
C     REFERENCES
C
C     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P.
C         The periodic Schur decomposition: Algorithms and applications.
C         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42,
C         1992.
C
C     [2] Kressner, D.
C         An efficient and reliable implementation of the periodic QZ
C         algorithm. Proc. of the IFAC Workshop on Periodic Control
C         Systems, pp. 187-192, 2001.
C
C     [3] Van Loan, C.
C         Generalized Singular Values with Algorithms and Applications.
C         Ph. D. Thesis, University of Michigan, 1973.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(N**3) floating point operations and is
C     backward stable.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPQR).
C
C     KEYWORDS
C
C     Eigenvalue, eigenvalue decomposition, Hessenberg form, orthogonal
C     transformation, (periodic) Schur form
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTT, WANTZ
      INTEGER            IHI, IHIQ, ILO, ILOQ, INFO, LDA, LDB, LDQ,
     $                   LDWORK, LDZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), ALPHAI(*), ALPHAR(*), B(LDB,*),
     $                   BETA(*), DWORK(*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, K, KK, L, NH, NQ, NR
      DOUBLE PRECISION   ALPHA, BETAX, CS1, CS2, CS3, DELTA, GAMMA,
     $                   OVFL, SMLNUM, SN1, SN2, SN3, TAUV, TAUW,
     $                   TEMP, TST, ULP, UNFL
C     .. Local Arrays ..
      INTEGER            ISEED(4)
      DOUBLE PRECISION   V(3), W(3)
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANHS
      EXTERNAL           DLAMCH, DLANHS
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DLABAD, DLARFG, DLARFX, DLARNV, DLARTG,
     $                   DROT, MB03YA, MB03YT, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      NH = IHI - ILO + 1
      NQ = IHIQ - ILOQ + 1
      IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -6
      ELSE IF ( ILOQ.LT.1 .OR. ILOQ.GT.ILO ) THEN
         INFO = -7
      ELSE IF ( IHIQ.LT.IHI .OR. IHIQ.GT.N ) THEN
         INFO = -8
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF ( LDQ.LT.1 .OR. WANTQ .AND. LDQ.LT.N ) THEN
         INFO = -14
      ELSE IF ( LDZ.LT.1 .OR. WANTZ .AND. LDZ.LT.N ) THEN
         INFO = -16
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = DBLE( MAX( 1, N ) )
         INFO = -21
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03YD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
C     Set machine-dependent constants for the stopping criterion.
C
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
C
C     I1 and I2 are the indices of the first rows and last columns of
C     A and B to which transformations must be applied.
C
      I1 = 1
      I2 = N
      ISEED(1) = 1
      ISEED(2) = 0
      ISEED(3) = 0
      ISEED(4) = 1
C
C     ITN is the maximal number of QR iterations.
C
      ITN = 30*NH
C
C     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO
C     or A(L,L-1) is negligible.
C
      I = IHI
   10 CONTINUE
      L = ILO
      IF ( I.LT.ILO )
     $   GO TO 120
C
C     Perform periodic QR iteration on rows and columns ILO to I of A
C     and B until a submatrix of order 1 or 2 splits off at the bottom.
C
      DO 70  ITS = 0, ITN
C
C        Look for deflations in A.
C
         DO 20  K = I, L + 1, -1
            TST = ABS( A(K-1,K-1) ) + ABS( A(K,K) )
            IF ( TST.EQ.ZERO )
     $         TST = DLANHS( '1', I-L+1, A(L,L), LDA, DWORK )
            IF ( ABS( A(K,K-1) ).LE.MAX( ULP*TST, SMLNUM ) )
     $         GO TO 30
   20    CONTINUE
   30    CONTINUE
C
C        Look for deflation in B if problem size is greater than 1.
C
         IF ( I-K.GE.1 ) THEN
            DO 40  KK = I, K, -1
               IF ( KK.EQ.I ) THEN
                  TST = ABS( B(KK-1,KK) )
               ELSE IF ( KK.EQ.K ) THEN
                  TST = ABS( B(KK,KK+1) )
               ELSE
                  TST = ABS( B(KK-1,KK) ) + ABS( B(KK,KK+1) )
               END IF
               IF ( TST.EQ.ZERO )
     $            TST = DLANHS( '1', I-K+1, B(K,K), LDB, DWORK )
               IF ( ABS( B(KK,KK) ).LE.MAX( ULP*TST, SMLNUM ) )
     $            GO TO 50
   40       CONTINUE
         ELSE
            KK = K-1
         END IF
   50    CONTINUE
         IF ( KK.GE.K ) THEN
C
C           B has an element close to zero at position (KK,KK).
C
            B(KK,KK) = ZERO
            CALL MB03YA( WANTT, WANTQ, WANTZ, N, K, I, ILOQ, IHIQ, KK,
     $                   A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO )
            K = KK+1
         END IF
         L = K
         IF( L.GT.ILO ) THEN
C
C           A(L,L-1) is negligible.
C
            A(L,L-1) = ZERO
         END IF
C
C        Exit from loop if a submatrix of order 1 or 2 has split off.
C
         IF ( L.GE.I-1 )
     $      GO TO 80
C
C        The active submatrices are now in rows and columns L:I.
C
         IF ( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
         IF ( ITS.EQ.10.OR.ITS.EQ.20 ) THEN
C
C           Exceptional shift. The first column of the shift polynomial
C           is a pseudo-random vector.
C
            CALL DLARNV( 3, ISEED, 3, V )
         ELSE
C
C           The implicit double shift is constructed via a partial
C           product QR factorization [2].
C
            CALL DLARTG( B(L,L), B(I,I), CS2, SN2, TEMP )
            CALL DLARTG( TEMP, B(I-1,I), CS1, SN1, ALPHA )
C
            ALPHA = A(L,L)*CS2 - A(I,I)*SN2
            BETAX = CS1*( CS2*A(L+1,L) )
            GAMMA = CS1*( SN2*A(I-1,I) ) + SN1*A(I-1,I-1)
            ALPHA = ALPHA*CS1 - A(I,I-1)*SN1
            CALL DLARTG( ALPHA, BETAX, CS1, SN1, TEMP )
C
            CALL DLARTG( TEMP, GAMMA, CS2, SN2, ALPHA )
            ALPHA = CS2
            GAMMA = ( A(I-1,I-1)*CS1 )*CS2 + A(I,I-1)*SN2
            DELTA = ( A(I-1,I-1)*SN1 )*CS2
            CALL DLARTG( GAMMA, DELTA, CS3, SN3, TEMP )
            CALL DLARTG( ALPHA, TEMP, CS2, SN2, ALPHA )
C
            ALPHA = ( B(L,L)*CS1 + B(L,L+1)*SN1 )*CS2
            BETAX = ( B(L+1,L+1)*SN1 )*CS2
            GAMMA = B(I-1,I-1)*SN2
            CALL DLARTG( ALPHA, BETAX, CS1, SN1, TEMP )
            CALL DLARTG( TEMP, GAMMA, CS2, SN2, ALPHA )
C
            ALPHA = CS1*A(L,L) + SN1*A(L,L+1)
            BETAX = CS1*A(L+1,L) + SN1*A(L+1,L+1)
            GAMMA = SN1*A(L+2,L+1)
C
            V(1) = CS2*ALPHA - SN2*CS3
            V(2) = CS2*BETAX - SN2*SN3
            V(3) = GAMMA*CS2
         END IF
C
C        Double-shift QR step
C
         DO 60  K = L, I-1
C
            NR = MIN( 3,I-K+1 )
            IF ( K.GT.L )
     $         CALL DCOPY( NR, A(K,K-1), 1, V, 1 )
            CALL DLARFG( NR, V(1), V(2), 1, TAUV )
            IF ( K.GT.L ) THEN
               A(K,K-1) = V(1)
               A(K+1,K-1) = ZERO
               IF ( K.LT.I-1 )
     $            A(K+2,K-1) = ZERO
            END IF
C
C           Apply reflector V from the right to B in rows I1:min(K+2,I).
C
            V(1) = ONE
            CALL DLARFX( 'Right', MIN(K+2,I)-I1+1, NR, V, TAUV, B(I1,K),
     $                   LDB, DWORK )
C
C           Annihilate the introduced nonzeros in the K-th column.
C
            CALL DCOPY( NR, B(K,K), 1, W, 1 )
            CALL DLARFG( NR, W(1), W(2), 1, TAUW )
            B(K,K) = W(1)
            B(K+1,K) = ZERO
            IF ( K.LT.I-1 )
     $         B(K+2,K) = ZERO
C
C           Apply reflector W from the left to transform the rows of the
C           matrix B in columns K+1:I2.
C
            W(1) = ONE
            CALL DLARFX( 'Left', NR, I2-K, W, TAUW, B(K,K+1), LDB,
     $                   DWORK )
C
C           Apply reflector V from the left to transform the rows of the
C           matrix A in columns K:I2.
C
            CALL DLARFX( 'Left', NR, I2-K+1, V, TAUV, A(K,K), LDA,
     $                   DWORK )
C
C           Apply reflector W from the right to transform the columns of
C           the matrix A in rows I1:min(K+3,I).
C
            CALL DLARFX( 'Right', MIN(K+3,I)-I1+1, NR, W, TAUW, A(I1,K),
     $                   LDA, DWORK )
C
C           Accumulate transformations in the matrices Q and Z.
C
            IF ( WANTQ )
     $         CALL DLARFX( 'Right', NQ, NR, V, TAUV, Q(ILOQ,K), LDQ,
     $                      DWORK )
            IF ( WANTZ )
     $         CALL DLARFX( 'Right', NQ, NR, W, TAUW, Z(ILOQ,K), LDZ,
     $                      DWORK )
   60    CONTINUE
   70 CONTINUE
C
C     Failure to converge.
C
      INFO = I
      RETURN
C
   80 CONTINUE
C
C     Compute 1-by-1 or 2-by-2 subproblem.
C
      IF ( L.EQ.I ) THEN
C
C        Standardize B, set ALPHAR, ALPHAI and BETA.
C
         IF ( B(I,I).LT.ZERO ) THEN
            IF ( WANTT ) THEN
               DO 90  K = I1, I
                  B(K,I) = -B(K,I)
   90          CONTINUE
               DO 100  K = I, I2
                  A(I,K) = -A(I,K)
  100          CONTINUE
            ELSE
               B(I,I) = -B(I,I)
               A(I,I) = -A(I,I)
            END IF
            IF ( WANTQ ) THEN
               DO 110  K = ILOQ, IHIQ
                  Q(K,I) = -Q(K,I)
  110          CONTINUE
            END IF
         END IF
         ALPHAR(I) = A(I,I)
         ALPHAI(I) = ZERO
         BETA(I) = B(I,I)
      ELSE IF( L.EQ.I-1 ) THEN
C
C        A double block has converged.
C        Compute eigenvalues and standardize double block.
C
         CALL MB03YT( A(I-1,I-1), LDA, B(I-1,I-1), LDB, ALPHAR(I-1),
     $                ALPHAI(I-1), BETA(I-1), CS1, SN1, CS2, SN2 )
C
C        Apply transformation to rest of A and B.
C
         IF ( I2.GT.I )
     $      CALL DROT( I2-I, A(I-1,I+1), LDA, A(I,I+1), LDA, CS1, SN1 )
         CALL DROT( I-I1-1, A(I1,I-1), 1, A(I1,I), 1, CS2, SN2 )
         IF ( I2.GT.I )
     $      CALL DROT( I2-I, B(I-1,I+1), LDB, B(I,I+1), LDB, CS2, SN2 )
         CALL DROT( I-I1-1, B(I1,I-1), 1, B(I1,I), 1, CS1, SN1 )
C
C        Apply transformation to rest of Q and Z if desired.
C
         IF ( WANTQ )
     $      CALL DROT( NQ, Q(ILOQ,I-1), 1, Q(ILOQ,I), 1, CS1, SN1 )
         IF ( WANTZ )
     $      CALL DROT( NQ, Z(ILOQ,I-1), 1, Z(ILOQ,I), 1, CS2, SN2 )
      END IF
C
C     Decrement number of remaining iterations, and return to start of
C     the main loop with new value of I.
C
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
C
  120 CONTINUE
      DWORK(1) = DBLE( MAX( 1, N ) )
      RETURN
C *** Last line of MB03YD ***
      END
