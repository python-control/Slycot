      SUBROUTINE MB03XP( JOB, COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,
     $                   Q, LDQ, Z, LDZ, ALPHAR, ALPHAI, BETA, DWORK,
     $                   LDWORK, INFO )
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
C     To compute the periodic Schur decomposition and the eigenvalues of
C     a product of matrices, H = A*B, with A upper Hessenberg and B
C     upper triangular without evaluating any part of the product.
C     Specifically, the matrices Q and Z are computed, so that
C
C          Q' * A * Z = S,    Z' * B * Q = T
C
C     where S is in real Schur form, and T is upper triangular.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates whether the user wishes to compute the full
C             Schur form or the eigenvalues only, as follows:
C             = 'E':  Compute the eigenvalues only;
C             = 'S':  compute the factors S and T of the full
C                     Schur form.
C
C     COMPQ   CHARACTER*1
C             Indicates whether or not the user wishes to accumulate
C             the matrix Q as follows:
C             = 'N':  The matrix Q is not required;
C             = 'I':  Q is initialized to the unit matrix and the
C                     orthogonal transformation matrix Q is returned;
C             = 'V':  Q must contain an orthogonal matrix U on entry,
C                     and the product U*Q is returned.
C
C     COMPZ   CHARACTER*1
C             Indicates whether or not the user wishes to accumulate
C             the matrix Z as follows:
C             = 'N':  The matrix Z is not required;
C             = 'I':  Z is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z is returned;
C             = 'V':  Z must contain an orthogonal matrix U on entry,
C                     and the product U*Z is returned.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and B. N >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that the matrices A and B are already upper
C             triangular in rows and columns 1:ILO-1 and IHI+1:N.
C             The routine works primarily with the submatrices in rows
C             and columns ILO to IHI, but applies the transformations to
C             all the rows and columns of the matrices A and B, if
C             JOB = 'S'.
C             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array A must
C             contain the upper Hessenberg matrix A.
C             On exit, if JOB = 'S', the leading N-by-N part of this
C             array is upper quasi-triangular with any 2-by-2 diagonal
C             blocks corresponding to a pair of complex conjugated
C             eigenvalues.
C             If JOB = 'E', the diagonal elements and 2-by-2 diagonal
C             blocks of A will be correct, but the remaining parts of A
C             are unspecified on exit.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array B must
C             contain the upper triangular matrix B.
C             On exit, if JOB = 'S', the leading N-by-N part of this
C             array contains the transformed upper triangular matrix.
C             2-by-2 blocks in B corresponding to 2-by-2 blocks in A
C             will be reduced to positive diagonal form. (I.e., if
C             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j)
C             and B(j+1,j+1) will be positive.)
C             If JOB = 'E', the elements corresponding to diagonal
C             elements and 2-by-2 diagonal blocks in A will be correct,
C             but the remaining parts of B are unspecified on exit.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, if COMPQ = 'V', then the leading N-by-N part of
C             this array must contain a matrix Q which is assumed to be
C             equal to the unit matrix except for the submatrix
C             Q(ILO:IHI,ILO:IHI).
C             If COMPQ = 'I', Q need not be set on entry.
C             On exit, if COMPQ = 'V' or COMPQ = 'I' the leading N-by-N
C             part of this array contains the transformation matrix
C             which produced the Schur form.
C             If COMPQ = 'N', Q is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= 1.
C             If COMPQ <> 'N', LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, if COMPZ = 'V', then the leading N-by-N part of
C             this array must contain a matrix Z which is assumed to be
C             equal to the unit matrix except for the submatrix
C             Z(ILO:IHI,ILO:IHI).
C             If COMPZ = 'I', Z need not be set on entry.
C             On exit, if COMPZ = 'V' or COMPZ = 'I' the leading N-by-N
C             part of this array contains the transformation matrix
C             which produced the Schur form.
C             If COMPZ = 'N', Z is not referenced.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= 1.
C             If COMPZ <> 'N', LDZ >= MAX(1,N).
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
C     BETA    (output) DOUBLE PRECISION array, dimension (N)
C             The i-th (1 <= i <= N) computed eigenvalue is given by
C             BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two
C             eigenvalues are computed as a complex conjugate pair,
C             they are stored in consecutive elements of ALPHAR, ALPHAI
C             and BETA. If JOB = 'S', the eigenvalues are stored in the
C             same order as on the diagonales of the Schur forms of A
C             and B.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
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
C             > 0:  if INFO = i, then MB03XP failed to compute the Schur
C                   form in a total of 30*(IHI-ILO+1) iterations;
C                   elements 1:ilo-1 and i+1:n of ALPHAR, ALPHAI and
C                   BETA contain successfully computed eigenvalues.
C
C     METHOD
C
C     The implemented algorithm is a multi-shift version of the periodic
C     QR algorithm described in [1,3] with some minor modifications
C     proposed in [2].
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
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHGPQR).
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
      INTEGER            NSMAX, LDAS, LDBS
      PARAMETER          ( NSMAX = 15, LDAS = NSMAX, LDBS = NSMAX )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDWORK, LDZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), ALPHAI(*), ALPHAR(*), B(LDB,*),
     $                   BETA(*), DWORK(*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL            INITQ, INITZ, WANTQ, WANTT, WANTZ
      INTEGER            DUM, I, I1, I2, IERR, ITEMP, ITN, ITS, J, K,
     $                   KK, L, MAXB, NH, NR, NS, NV, PV2, PV3
      DOUBLE PRECISION   OVFL, SMLNUM, TAUV, TAUW, TEMP, TST, ULP, UNFL
C     .. Local Arrays ..
      INTEGER            ISEED(4)
      DOUBLE PRECISION   AS(LDAS,LDAS), BS(LDBS,LDBS), V(3*NSMAX+6)
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, UE01MD
      DOUBLE PRECISION   DLAMCH, DLANHS
      EXTERNAL           DLAMCH, DLANHS, IDAMAX, LSAME, UE01MD
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DLABAD, DLACPY, DLARFG,
     $                   DLARFX, DLARNV, DLASET, DSCAL, DTRMV, MB03YA,
     $                   MB03YD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      WANTT = LSAME( JOB, 'S' )
      INITQ = LSAME( COMPQ, 'I' )
      WANTQ = INITQ.OR.LSAME( COMPQ, 'V' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ.OR.LSAME( COMPZ, 'V' )
C
C     Check the scalar input parameters.
C
      INFO = 0
      IF ( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF ( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF ( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -3
      ELSE IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF ( IHI.LT.MIN( ILO,N ).OR.IHI.GT.N ) THEN
         INFO = -6
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF ( LDQ.LT.1 .OR. WANTQ .AND. LDQ.LT.N ) THEN
         INFO = -12
      ELSE IF ( LDZ.LT.1 .OR. WANTZ .AND. LDZ.LT.N ) THEN
         INFO = -14
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = DBLE( MAX( 1, N ) )
         INFO = -19
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03XP', -INFO )
         RETURN
      END IF
C
C     Initialize Q and Z, if necessary.
C
      IF ( INITQ )
     $   CALL DLASET( 'All', N, N, ZERO, ONE, Q, LDQ )
      IF ( INITZ )
     $   CALL DLASET( 'All', N, N, ZERO, ONE, Z, LDZ )
C
C     Store isolated eigenvalues and standardize B.
C
C     FOR I = [1:ILO-1, IHI+1:N]
      I = 1
   10 CONTINUE
      IF ( I.EQ.ILO ) THEN
         I = IHI+1
      END IF
      IF ( I.LE.N ) THEN
         IF ( B(I,I).LT.ZERO ) THEN
            IF ( WANTT ) THEN
               DO 20  K = ILO, I
                  B(K,I) = -B(K,I)
   20          CONTINUE
               DO 30  K = I, IHI
                  A(I,K) = -A(I,K)
   30          CONTINUE
            ELSE
               B(I,I) = -B(I,I)
               A(I,I) = -A(I,I)
            END IF
            IF ( WANTQ ) THEN
               DO 40  K = ILO, IHI
                  Q(K,I) = -Q(K,I)
   40          CONTINUE
            END IF
         END IF
         ALPHAR(I) = A(I,I)
         ALPHAI(I) = ZERO
         BETA(I) = B(I,I)
         I = I + 1
C        END FOR
         GO TO 10
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. ILO.EQ.IHI+1 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Set rows and coloms ILO to IHI of B (A) to zero below the first
C     (sub)diagonal.
C
      DO  60 J = ILO, IHI - 2
         DO 50  I = J + 2, N
            A(I,J) = ZERO
   50    CONTINUE
   60 CONTINUE
      DO  80 J = ILO, IHI - 1
         DO 70  I = J + 1, N
            B(I,J) = ZERO
   70    CONTINUE
   80 CONTINUE
      NH = IHI - ILO + 1
C
C     Suboptimal choice of the number of shifts.
C
      IF ( WANTQ ) THEN
         NS   = UE01MD( 4, 'MB03XP', JOB // COMPQ, N, ILO, IHI )
         MAXB = UE01MD( 8, 'MB03XP', JOB // COMPQ, N, ILO, IHI )
      ELSE
         NS   = UE01MD( 4, 'MB03XP', JOB // COMPZ, N, ILO, IHI )
         MAXB = UE01MD( 8, 'MB03XP', JOB // COMPZ, N, ILO, IHI )
      END IF
C
      IF ( NS.LE.2 .OR. NS.GT.NH .OR. MAXB.GE.NH ) THEN
C
C        Standard double-shift product QR.
C
         CALL MB03YD( WANTT, WANTQ, WANTZ, N, ILO, IHI, ILO, IHI, A,
     $                LDA, B, LDB, Q, LDQ, Z, LDZ, ALPHAR, ALPHAI, BETA,
     $                DWORK, LDWORK, INFO )
         RETURN
      END IF
      MAXB = MAX( 3, MAXB )
      NS = MIN( NS, MAXB, NSMAX )
C
C     Set machine-dependent constants for the stopping criterion.
C     If max(norm(A),norm(B)) <= sqrt(OVFL), then overflow should not
C     occur.
C
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( DBLE( NH ) / ULP )
C
C     I1 and I2 are the indices of the first rows and last columns of
C     A and B to which transformations must be applied.
C
      IF ( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
      ISEED(1) = 1
      ISEED(2) = 0
      ISEED(3) = 0
      ISEED(4) = 1
C
C     ITN is the maximal number of QR iterations.
C
      ITN = 30*NH
      DUM = 0
C
C     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO
C     or A(L,L-1) is negligible.
C
      I = IHI
   90 CONTINUE
      L = ILO
      IF ( I.LT.ILO )
     $   GO TO 210
C
      DO 190  ITS = 0, ITN
         DUM = DUM + (IHI-ILO)*(IHI-ILO)
C
C        Look for deflations in A.
C
         DO 100  K = I, L + 1, -1
            TST = ABS( A(K-1,K-1) ) + ABS( A(K,K) )
            IF ( TST.EQ.ZERO )
     $         TST = DLANHS( '1', I-L+1, A(L,L), LDA, DWORK )
            IF ( ABS( A(K,K-1) ).LE.MAX( ULP*TST, SMLNUM ) )
     $         GO TO 110
  100    CONTINUE
  110    CONTINUE
C
C        Look for deflation in B if problem size is greater than 1.
C
         IF ( I-K.GE.1 ) THEN
            DO 120  KK = I, K, -1
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
     $            GO TO 130
  120       CONTINUE
         ELSE
            KK = K-1
         END IF
  130    CONTINUE
         IF ( KK.GE.K ) THEN
C
C           B has an element close to zero at position (KK,KK).
C
            B(KK,KK) = ZERO
            CALL MB03YA( WANTT, WANTQ, WANTZ, N, K, I, ILO, IHI, KK,
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
C        Exit from loop if a submatrix of order <= MAXB has split off.
C
         IF ( L.GE.I-MAXB+1 )
     $      GO TO 200
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
            CALL DLARNV( 3, ISEED, NS+1, V )
         ELSE
C
C           Use eigenvalues of trailing submatrix as shifts.
C
            CALL DLACPY( 'Full', NS, NS, A(I-NS+1,I-NS+1), LDA, AS,
     $                   LDAS )
            CALL DLACPY( 'Full', NS, NS, B(I-NS+1,I-NS+1), LDB, BS,
     $                   LDBS )
            CALL MB03YD( .FALSE., .FALSE., .FALSE., NS, 1, NS, 1, NS,
     $                   AS, LDAS, BS, LDBS, Q, LDQ, Z, LDZ,
     $                   ALPHAR(I-NS+1), ALPHAI(I-NS+1), BETA(I-NS+1),
     $                   DWORK, LDWORK, IERR )
         END IF
C
C        Compute the nonzero elements of the first column of
C        (A*B-w(1)) (A*B-w(2)) .. (A*B-w(ns)).
C
         V(1) = ONE
         NV = 1
C        WHILE NV <= NS
  140    CONTINUE
         IF ( NV.LE.NS ) THEN
            IF ( NV.EQ.NS .OR. AS(NV+1,NV).EQ.ZERO ) THEN
C
C              Real shift.
C
               V(NV+1) = ZERO
               PV2 = NV+2
               CALL DCOPY( NV, V, 1, V(PV2), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'No unit diagonal',
     $                     NV, B(L,L), LDB, V(PV2), 1 )
               CALL DSCAL( NV, BS(NV,NV), V, 1 )
               ITEMP = IDAMAX( 2*NV+1, V, 1 )
               TEMP = ONE / MAX( ABS( V(ITEMP) ), SMLNUM )
               CALL DSCAL( 2*NV+1, TEMP, V, 1 )
               CALL DGEMV( 'No transpose', NV+1, NV, ONE, A(L,L), LDA,
     $                     V(PV2), 1, -AS(NV,NV), V, 1 )
               NV = NV + 1
            ELSE
C
C              Double shift using a product formulation of the shift
C              polynomial [2].
C
               V(NV+1) = ZERO
               V(NV+2) = ZERO
               PV2 = NV+3
               PV3 = 2*NV+5
               CALL DCOPY( NV+2, V, 1, V(PV2), 1 )
               CALL DCOPY( NV+1, V, 1, V(PV3), 1 )
               CALL DSCAL( NV, BS(NV+1,NV+1), V(PV2), 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'No unit diagonal',
     $                     NV, B(L,L), LDB, V(PV3), 1 )
               ITEMP = IDAMAX( 2*NV+3, V(PV2), 1 )
               TEMP = ONE / MAX( ABS( V(PV2+ITEMP-1) ), SMLNUM )
               CALL DSCAL( 2*NV+3, TEMP, V(PV2), 1 )
C
               CALL DCOPY( NV, V(PV2), 1, V, 1 )
               CALL DGEMV( 'No transpose', NV+1, NV, -ONE, A(L,L), LDA,
     $                     V(PV3), 1, AS(NV+1,NV+1), V(PV2), 1 )
               CALL DSCAL( NV, AS(NV,NV+1), V, 1 )
               ITEMP = IDAMAX( 2*NV+3, V, 1 )
               TEMP = ONE / MAX( ABS( V(ITEMP) ), SMLNUM )
               CALL DSCAL( 2*NV+3, TEMP, V, 1 )
C
               CALL DSCAL( NV, -AS(NV+1,NV), V, 1 )
               CALL DAXPY( NV+1, AS(NV,NV), V(PV2), 1, V, 1)
               ITEMP = IDAMAX( 2*NV+3, V, 1 )
               TEMP = ONE / MAX( ABS( V(ITEMP) ), SMLNUM )
               CALL DSCAL( 2*NV+3, TEMP, V, 1 )
C
               CALL DSCAL( NV+1, BS(NV,NV), V, 1 )
               CALL DTRMV( 'Upper', 'No transpose', 'No unit diagonal',
     $                     NV+1, B(L,L), LDB, V(PV2), 1 )
               ITEMP = IDAMAX( 2*NV+3, V, 1 )
               TEMP = ONE / MAX( ABS( V(ITEMP) ), SMLNUM )
               CALL DSCAL( 2*NV+3, TEMP, V, 1 )
C
               CALL DGEMV( 'No transpose', NV+2, NV+1, -ONE, A(L,L),
     $                     LDA, V(PV2), 1, ONE, V, 1 )
               NV = NV + 2
            END IF
            ITEMP = IDAMAX( NV, V, 1 )
            TEMP = ABS( V(ITEMP) )
            IF ( TEMP.EQ.ZERO ) THEN
               V(1) = ONE
               DO 150  K = 2, NV
                  V(K) = ZERO
  150          CONTINUE
            ELSE
               TEMP = MAX( TEMP, SMLNUM )
               CALL DSCAL( NV, ONE/TEMP, V, 1 )
            END IF
            GO TO 140
C        END WHILE
         END IF
C
C        Multi-shift product QR step.
C
         PV2 = NS+2
         DO 180  K = L,I-1
            NR = MIN( NS+1,I-K+1 )
            IF ( K.GT.L )
     $         CALL DCOPY( NR, A(K,K-1), 1, V, 1 )
            CALL DLARFG( NR, V(1), V(2), 1, TAUV )
            IF ( K.GT.L ) THEN
               A(K,K-1) = V(1)
               DO 160  KK = K+1,I
                  A(KK,K-1) = ZERO
  160          CONTINUE
            END IF
C
C           Apply reflector V from the right to B in rows
C           I1:min(K+NS,I).
C
            V(1) = ONE
            CALL DLARFX( 'Right', MIN(K+NS,I)-I1+1, NR, V, TAUV,
     $                   B(I1,K), LDB, DWORK )
C
C           Annihilate the introduced nonzeros in the K-th column.
C
            CALL DCOPY( NR, B(K,K), 1, V(PV2), 1 )
            CALL DLARFG( NR, V(PV2), V(PV2+1), 1, TAUW )
            B(K,K) = V(PV2)
            DO 170  KK = K+1,I
               B(KK,K) = ZERO
  170       CONTINUE
            V(PV2) = ONE
C
C           Apply reflector W from the left to transform the rows of the
C           matrix B in columns K+1:I2.
C
            CALL DLARFX( 'Left', NR, I2-K, V(PV2), TAUW, B(K,K+1), LDB,
     $                   DWORK )
C
C           Apply reflector V from the left to transform the rows of the
C           matrix A in columns K:I2.
C
            CALL DLARFX( 'Left', NR, I2-K+1, V, TAUV, A(K,K), LDA,
     $                   DWORK )
C
C           Apply reflector W from the right to transform the columns of
C           the matrix A in rows I1:min(K+NS,I).
C
            CALL DLARFX( 'Right', MIN(K+NS+1,I)-I1+1, NR, V(PV2), TAUW,
     $                   A(I1,K), LDA, DWORK )
C
C           Accumulate transformations in the matrices Q and Z.
C
            IF ( WANTQ )
     $         CALL DLARFX( 'Right', NH, NR, V, TAUV, Q(ILO,K), LDQ,
     $                      DWORK )
            IF ( WANTZ )
     $         CALL DLARFX( 'Right', NH, NR, V(PV2), TAUW, Z(ILO,K),
     $                      LDZ, DWORK )
  180    CONTINUE
  190 CONTINUE
C
C     Failure to converge.
C
      INFO = I
      RETURN
  200 CONTINUE
C
C     Submatrix of order <= MAXB has split off. Use double-shift
C     periodic QR algorithm.
C
      CALL MB03YD( WANTT, WANTQ, WANTZ, N, L, I, ILO, IHI, A, LDA, B,
     $             LDB, Q, LDQ, Z, LDZ, ALPHAR, ALPHAI, BETA, DWORK,
     $             LDWORK, INFO )
      IF ( INFO.GT.0 )
     $   RETURN
      ITN = ITN - ITS
      I = L - 1
      GO TO 90
C
  210 CONTINUE
      DWORK(1) = DBLE( MAX( 1,N ) )
      RETURN
C *** Last line of MB03XP ***
      END
