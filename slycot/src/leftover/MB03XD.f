      SUBROUTINE MB03XD( BALANC, JOB, JOBU, JOBV, N, A, LDA, QG, LDQG,
     $                   T, LDT, U1, LDU1, U2, LDU2, V1, LDV1, V2, LDV2,
     $                   WR, WI, ILO, SCALE, DWORK, LDWORK, INFO )
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
C     To compute the eigenvalues of a Hamiltonian matrix,
C
C                   [  A   G  ]         T        T
C             H  =  [       T ],   G = G,   Q = Q,                  (1)
C                   [  Q  -A  ]
C
C     where A, G and Q are real n-by-n matrices.
C
C     Due to the structure of H all eigenvalues appear in pairs
C     (lambda,-lambda). This routine computes the eigenvalues of H
C     using an algorithm based on the symplectic URV and the periodic
C     Schur decompositions as described in [1],
C
C           T       [  T   G  ]
C          U H V =  [       T ],                                    (2)
C                   [  0  -S  ]
C
C     where U and V are 2n-by-2n orthogonal symplectic matrices,
C     S is in real Schur form and T is upper triangular.
C
C     The algorithm is backward stable and preserves the eigenvalue
C     pairings in finite precision arithmetic.
C
C     Optionally, a symplectic balancing transformation to improve the
C     conditioning of eigenvalues is computed (see MB04DD). In this
C     case, the matrix H in decomposition (2) must be replaced by the
C     balanced matrix.
C
C     The SLICOT Library routine MB03ZD can be used to compute invariant
C     subspaces of H from the output of this routine.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALANC  CHARACTER*1
C             Indicates how H should be diagonally scaled and/or
C             permuted to reduce its norm.
C             = 'N': Do not diagonally scale or permute;
C             = 'P': Perform symplectic permutations to make the matrix
C                    closer to Hamiltonian Schur form. Do not diagonally
C                    scale;
C             = 'S': Diagonally scale the matrix, i.e., replace A, G and
C                    Q by D*A*D**(-1), D*G*D and D**(-1)*Q*D**(-1) where
C                    D is a diagonal matrix chosen to make the rows and
C                    columns of H more equal in norm. Do not permute;
C             = 'B': Both diagonally scale and permute A, G and Q.
C             Permuting does not change the norm of H, but scaling does.
C
C     JOB     CHARACTER*1
C             Indicates whether the user wishes to compute the full
C             decomposition (2) or the eigenvalues only, as follows:
C             = 'E': compute the eigenvalues only;
C             = 'S': compute matrices T and S of (2);
C             = 'G': compute matrices T, S and G of (2).
C
C     JOBU    CHARACTER*1
C             Indicates whether or not the user wishes to compute the
C             orthogonal symplectic matrix U of (2) as follows:
C             = 'N': the matrix U is not computed;
C             = 'U': the matrix U is computed.
C
C     JOBV    CHARACTER*1
C             Indicates whether or not the user wishes to compute the
C             orthogonal symplectic matrix V of (2) as follows:
C             = 'N': the matrix V is not computed;
C             = 'V': the matrix V is computed.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C             On exit, this array is overwritten. If JOB = 'S' or
C             JOB = 'G', the leading N-by-N part of this array contains
C             the matrix S in real Schur form of decomposition (2).
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C                            (LDQG,N+1)
C             On entry, the leading N-by-N+1 part of this array must
C             contain in columns 1:N the lower triangular part of the
C             matrix Q and in columns 2:N+1 the upper triangular part
C             of the matrix G.
C             On exit, this array is overwritten. If JOB = 'G', the
C             leading N-by-N+1 part of this array contains in columns
C             2:N+1 the matrix G of decomposition (2).
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= max(1,N).
C
C     T       (output) DOUBLE PRECISION array, dimension (LDT,N)
C             On exit, if JOB = 'S' or JOB = 'G', the leading N-by-N
C             part of this array contains the upper triangular matrix T
C             of the decomposition (2). Otherwise, this array is used as
C             workspace.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= MAX(1,N).
C
C     U1      (output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On exit, if JOBU = 'U', the leading N-by-N part of this
C             array contains the (1,1) block of the orthogonal
C             symplectic matrix U of decomposition (2).
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.  LDU1 >= 1.
C             LDU1 >= N,    if JOBU = 'U'.
C
C     U2      (output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On exit, if JOBU = 'U', the leading N-by-N part of this
C             array contains the (2,1) block of the orthogonal
C             symplectic matrix U of decomposition (2).
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.  LDU2 >= 1.
C             LDU2 >= N,    if JOBU = 'U'.
C
C     V1      (output) DOUBLE PRECISION array, dimension (LDV1,N)
C             On exit, if JOBV = 'V', the leading N-by-N part of this
C             array contains the (1,1) block of the orthogonal
C             symplectic matrix V of decomposition (2).
C
C     LDV1    INTEGER
C             The leading dimension of the array V1.  LDV1 >= 1.
C             LDV1 >= N,    if JOBV = 'V'.
C
C     V2      (output) DOUBLE PRECISION array, dimension (LDV2,N)
C             On exit, if JOBV = 'V', the leading N-by-N part of this
C             array contains the (2,1) block of the orthogonal
C             symplectic matrix V of decomposition (2).
C
C     LDV2    INTEGER
C             The leading dimension of the array V2.  LDV2 >= 1.
C             LDV2 >= N,    if JOBV = 'V'.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             On exit, the leading N elements of WR and WI contain the
C             real and imaginary parts, respectively, of N eigenvalues
C             that have nonpositive real part. Complex conjugate pairs
C             of eigenvalues with real part not equal to zero will
C             appear consecutively with the eigenvalue having the
C             positive imaginary part first. For complex conjugate pairs
C             of eigenvalues on the imaginary axis only the eigenvalue
C             having nonnegative imaginary part will be returned.
C
C     ILO     (output) INTEGER
C             ILO is an integer value determined when H was balanced.
C             The balanced A(i,j) = 0 if I > J and J = 1,...,ILO-1.
C             The balanced Q(i,j) = 0 if J = 1,...,ILO-1 or
C             I = 1,...,ILO-1.
C
C     SCALE   (output) DOUBLE PRECISION array, dimension (N)
C             On exit, if SCALE = 'S', the leading N elements of this
C             array contain details of the permutation and scaling
C             factors applied when balancing H, see MB04DD.
C             This array is not referenced if BALANC = 'N'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -25,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  (input) INTEGER
C             The dimension of the array DWORK. LDWORK >= max( 1, 8*N ).
C             Moreover:
C             If JOB = 'E' or 'S' and JOBU = 'N' and JOBV = 'N',
C                LDWORK >= 7*N+N*N.
C             If JOB = 'G' and JOBU = 'N' and JOBV = 'N',
C                LDWORK >= max( 7*N+N*N, 2*N+3*N*N ).
C             If JOB = 'G' and JOBU = 'U' and JOBV = 'N',
C                LDWORK >= 7*N+2*N*N.
C             If JOB = 'G' and JOBU = 'N' and JOBV = 'V',
C                LDWORK >= 7*N+2*N*N.
C             If JOB = 'G' and JOBU = 'U' and JOBV = 'V',
C                LDWORK >= 7*N+N*N.
C             For good performance, LDWORK must generally be larger.
C
C     Error Indicator
C
C     INFO     (output) INTEGER
C              = 0:  successful exit;
C              < 0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C              > 0:  if INFO = i, the periodic QR algorithm failed to
C                    compute all the eigenvalues, elements i+1:N of WR
C                    and WI contain eigenvalues which have converged.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A numerically stable, structure preserving method for
C         computing the eigenvalues of real Hamiltonian or symplectic
C         pencils.
C         Numer. Math., Vol. 78(3), pp. 329-358, 1998.
C
C     [2] Benner, P., Mehrmann, V., and Xu, H.
C         A new method for computing the stable invariant subspace of a
C         real Hamiltonian matrix,  J. Comput. Appl. Math., vol. 86,
C         pp. 17-43, 1997.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAESU).
C
C     KEYWORDS
C
C     Eigenvalues, invariant subspace, Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          BALANC, JOB, JOBU, JOBV
      INTEGER            ILO, INFO, LDA, LDQG, LDT, LDU1, LDU2, LDV1,
     $                   LDV2, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), DWORK(*), QG(LDQG,*), SCALE(*),
     $                   T(LDT,*), U1(LDU1,*), U2(LDU2,*), V1(LDV1,*),
     $                   V2(LDV2,*), WI(*), WR(*)
C     .. Local Scalars ..
      CHARACTER          UCHAR, VCHAR
      LOGICAL            LPERM, LSCAL, SCALEH, WANTG, WANTS, WANTU,
     $                   WANTV
      INTEGER            I, IERR, ILO1, J, K, L, PBETA, PCSL, PCSR, PDW,
     $                   PQ, PTAUL, PTAUR, PZ, WRKMIN, WRKOPT
      DOUBLE PRECISION   BIGNUM, CSCALE, EPS, HNRM, SMLNUM, TEMP, TEMPI,
     $                   TEMPR
C     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, MA02ID
      EXTERNAL           DLAMCH, LSAME, MA02ID
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLABAD, DLACPY, DLASCL, DLASET,
     $                   DSCAL, MA01AD, MB03XP, MB04DD, MB04QB, MB04TB,
     $                   XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO = 0
      LPERM = LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' )
      LSCAL = LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' )
      WANTS = LSAME( JOB,    'S' ) .OR. LSAME( JOB,    'G' )
      WANTG = LSAME( JOB,  'G' )
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
C
      IF ( WANTG ) THEN
         IF ( WANTU ) THEN
            IF ( WANTV ) THEN
               WRKMIN = MAX( 1, 7*N+N*N )
            ELSE
               WRKMIN = MAX( 1, 7*N+2*N*N )
            END IF
         ELSE
            IF ( WANTV ) THEN
               WRKMIN = MAX( 1, 7*N+2*N*N )
            ELSE
               WRKMIN = MAX( 1, 7*N+N*N, 2*N+3*N*N )
            END IF
         END IF
      ELSE
         IF ( WANTU ) THEN
            IF ( WANTV ) THEN
               WRKMIN = MAX( 1, 8*N )
            ELSE
               WRKMIN = MAX( 1, 8*N )
            END IF
         ELSE
            IF ( WANTV ) THEN
               WRKMIN = MAX( 1, 8*N )
            ELSE
               WRKMIN = MAX( 1, 7*N+N*N )
            END IF
         END IF
      END IF
C
      WRKOPT = WRKMIN
C
C     Test the scalar input parameters.
C
      IF ( .NOT.LPERM .AND. .NOT.LSCAL
     $     .AND. .NOT.LSAME( BALANC, 'N' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.WANTS .AND. .NOT.LSAME( JOB, 'E' ) ) THEN
         INFO = -2
      ELSE IF ( .NOT.WANTU .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -3
      ELSE IF ( .NOT.WANTV .AND. .NOT.LSAME( JOBV, 'N' ) ) THEN
         INFO = -4
      ELSE IF ( N.LT.0 ) THEN
         INFO = -5
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF ( LDQG.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF ( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF ( LDU1.LT.1 .OR. ( WANTU .AND. LDU1.LT.N ) ) THEN
         INFO = -13
      ELSE IF ( LDU2.LT.1 .OR. ( WANTU .AND. LDU2.LT.N ) ) THEN
         INFO = -15
      ELSE IF ( LDV1.LT.1 .OR. ( WANTV .AND. LDV1.LT.N ) ) THEN
         INFO = -17
      ELSE IF ( LDV2.LT.1 .OR. ( WANTV .AND. LDV2.LT.N ) ) THEN
         INFO = -19
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         INFO = -25
         DWORK(1) = DBLE( WRKMIN )
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03XD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      ILO = 0
      IF( N.EQ.0 )
     $   RETURN
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
C
C     Scale H if maximal element is outside range [SMLNUM,BIGNUM].
C
      HNRM = MA02ID( 'Hamiltonian', 'MaxElement', N, A, LDA, QG, LDQG,
     $               DWORK )
      SCALEH = .FALSE.
      IF ( HNRM.GT.ZERO .AND. HNRM.LT.SMLNUM ) THEN
         SCALEH = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( HNRM.GT.BIGNUM ) THEN
         SCALEH = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF ( SCALEH ) THEN
        CALL DLASCL( 'General', 0, 0, HNRM, CSCALE, N, N, A, LDA, IERR )
        CALL DLASCL( 'General', 0, 0, HNRM, CSCALE, N, N+1, QG, LDQG,
     $               IERR )
      END IF
C
C     Balance the matrix.
C
      CALL MB04DD( BALANC, N, A, LDA, QG, LDQG, ILO, SCALE, IERR )
C
C     Copy A to T and multiply A by -1.
C
      CALL DLACPY( 'All', N, N, A, LDA, T, LDT )
      CALL DLASCL( 'General', 0, 0, ONE, -ONE, N, N, A, LDA, IERR )
C
C     ---------------------------------------------
C     Step 1: Compute symplectic URV decomposition.
C     ---------------------------------------------
C
      PCSL  = 1
      PCSR  = PCSL  + 2*N
      PTAUL = PCSR  + 2*N
      PTAUR = PTAUL + N
      PDW   = PTAUR + N

      IF ( .NOT.WANTU .AND. .NOT.WANTV  ) THEN
C
C         Copy Q and Q' to workspace.
C
         PQ   = PDW
         PDW  = PDW + N*N
         DO 20  J = 1, N
            K = PQ + (N+1)*(J-1)
            L = K
            DWORK(K) = QG(J,J)
            DO 10  I = J+1, N
               K = K + 1
               L = L + N
               TEMP = QG(I,J)
               DWORK(K) = TEMP
               DWORK(L) = TEMP
   10       CONTINUE
   20    CONTINUE
      ELSE IF ( WANTU ) THEN
C
C         Copy Q and Q' to U2.
C
         DO 40  J = 1, N
            U2(J,J) = QG(J,J)
            DO 30 I = J+1, N
               TEMP = QG(I,J)
               U2(I,J) = TEMP
               U2(J,I) = TEMP
   30       CONTINUE
   40    CONTINUE
      ELSE
C
C         Copy Q and Q' to V2.
C
         DO 60  J = 1, N
            V2(J,J) = QG(J,J)
            DO 50 I = J+1, N
               TEMP = QG(I,J)
               V2(I,J) = TEMP
               V2(J,I) = TEMP
   50       CONTINUE
   60    CONTINUE
      END IF
C
C     Transpose G.
C
      DO 80 J = 1, N
         DO 70 I = J+1, N
            QG(I,J+1) = QG(J,I+1)
   70    CONTINUE
   80 CONTINUE
C
      IF ( .NOT.WANTU .AND. .NOT.WANTV  ) THEN
         CALL MB04TB( 'Not Transposed', 'Transposed', N, ILO, T, LDT, A,
     $                LDA, QG(1,2), LDQG, DWORK(PQ), N, DWORK(PCSL),
     $                DWORK(PCSR), DWORK(PTAUL), DWORK(PTAUR),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
      ELSE IF ( WANTU ) THEN
         CALL MB04TB( 'Not Transposed', 'Transposed', N, ILO, T, LDT, A,
     $                LDA, QG(1,2), LDQG, U2, LDU2, DWORK(PCSL),
     $                DWORK(PCSR), DWORK(PTAUL), DWORK(PTAUR),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
      ELSE
         CALL MB04TB( 'Not Transposed', 'Transposed', N, ILO, T, LDT, A,
     $                LDA, QG(1,2), LDQG, V2, LDV2, DWORK(PCSL),
     $                DWORK(PCSR), DWORK(PTAUL), DWORK(PTAUR),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
      END IF
      WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
      IF ( WANTU .AND. .NOT.WANTV .AND. .NOT.WANTG ) THEN
         IF ( N.GT.1 )
     $     CALL DLACPY( 'Lower', N-1, N-1, T(2,1), LDT, QG(2,1), LDQG )
      ELSE IF ( .NOT.WANTU .AND. WANTV .AND. .NOT.WANTG ) THEN
         IF ( N.GT.1 ) THEN
            CALL DLACPY( 'Lower', N-1, N-1, A(2,1), LDA, QG(2,1), LDQG )
            CALL DLACPY( 'Upper', N-1, N-1, V2(1,2), LDV2, QG(1,2),
     $                   LDQG )
         END IF
      ELSE IF ( WANTU .AND. WANTV .AND. .NOT.WANTG ) THEN
         IF ( N.GT.1 ) THEN
            CALL DLACPY( 'Lower', N-1, N-1, T(2,1), LDT, V2(2,1), LDV2 )
            CALL DLACPY( 'Lower', N-1, N-1, A(2,1), LDA, QG(2,1), LDQG )
         END IF
      ELSE IF ( WANTU .AND. .NOT.WANTV .AND. WANTG ) THEN
         IF ( N.GT.1 )
     $     CALL DLACPY( 'Lower', N-1, N-1, T(2,1), LDT,
     $                  DWORK(PDW+N*N+N), N-1 )
      ELSE IF ( .NOT.WANTU .AND. WANTV .AND. WANTG ) THEN
         IF ( N.GT.2 )
     $     CALL DLACPY( 'Lower', N-2, N-2, A(3,1), LDA,
     $                  DWORK(PDW+N*N+N), N-2 )
      ELSE IF ( WANTU .AND. WANTV .AND. WANTG ) THEN
         IF ( N.GT.1 )
     $     CALL DLACPY( 'Lower', N-1, N-1, T(2,1), LDT,
     $                  DWORK(PDW+N), N-1 )
         IF ( N.GT.2 )
     $     CALL DLACPY( 'Lower', N-2, N-2, A(3,1), LDA, V2(3,1), LDV2 )
      END IF
C
C     ----------------------------------------------
C     Step 2:  Compute periodic Schur decomposition.
C     ----------------------------------------------
C
      IF ( N.GT.2 )
     $   CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
      IF ( N.GT.1 )
     $   CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, T(2,1), LDT )
      IF ( .NOT.WANTU .AND. .NOT.WANTV ) THEN
         PBETA = 1
      ELSE
         PBETA = PDW
      END IF
C
      IF ( .NOT.WANTG ) THEN
C
C        Workspace requirements: 2*N (8*N with U or V).
C
         PDW = PBETA + N
         IF ( WANTU ) THEN
            UCHAR = 'I'
         ELSE
            UCHAR = 'N'
         END IF
         IF ( WANTV ) THEN
            VCHAR = 'I'
         ELSE
            VCHAR = 'N'
         END IF
         CALL MB03XP( JOB, VCHAR, UCHAR, N, ILO, N, A, LDA, T, LDT, V1,
     $                LDV1, U1, LDU1, WR, WI, DWORK(PBETA), DWORK(PDW),
     $                LDWORK-PDW+1, INFO )
         IF ( INFO.NE.0 )
     $      GO TO 90
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
      ELSE IF ( .NOT.WANTU .AND. .NOT.WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 3*N*N + 2*N.
C
         PQ    = PBETA + N
         PZ    = PQ + N*N
         PDW   = PZ + N*N
         CALL MB03XP( 'Schur', 'Init', 'Init', N, ILO, N, A, LDA, T,
     $                LDT, DWORK(PQ), N, DWORK(PZ), N, WR, WI,
     $                DWORK(PBETA), DWORK(PDW), LDWORK-PDW+1, INFO )
         IF ( INFO.NE.0 )
     $      GO TO 90
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PZ), N, QG(1,2), LDQG, ZERO, DWORK(PDW), N )
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PDW), N, DWORK(PQ), N, ZERO, QG(1,2), LDQG )
      ELSE IF ( WANTU .AND. .NOT.WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 2*N*N + 7*N.
C
         PQ    = PBETA + N
         PDW   = PQ + N*N
         CALL MB03XP( 'Schur', 'Init', 'Init', N, ILO, N, A, LDA, T,
     $                LDT, DWORK(PQ), N, U1, LDU1, WR, WI, DWORK(PBETA),
     $                DWORK(PDW+(N-1)*(N-1)), LDWORK-PDW-(N-1)*(N-1)+1,
     $                INFO )
         IF ( INFO.NE.0 )
     $      GO TO 90
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+(N-1)*(N-1)) ) + PDW
     $                         + (N-1)*(N-1) - 1 )
         IF ( N.GT.1 )
     $      CALL DLACPY( 'Lower', N-1, N-1, DWORK(PDW), N-1, T(2,1),
     $                   LDT )
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE,
     $                U1, LDU1, QG(1,2), LDQG, ZERO, DWORK(PDW), N )
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PDW), N, DWORK(PQ), N, ZERO, QG(1,2), LDQG )
C
      ELSE IF ( .NOT.WANTU .AND. WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 2*N*N + 7*N
C
         PZ    = PBETA + N
         PDW   = PZ + N*N
         CALL MB03XP( 'Schur', 'Init', 'Init', N, ILO, N, A, LDA, T,
     $                LDT, V1, LDV1, DWORK(PZ), N, WR, WI, DWORK(PBETA),
     $                DWORK(PDW+(N-1)*(N-1)), LDWORK-PDW-(N-1)*(N-1)+1,
     $                INFO )
         IF ( INFO.NE.0 )
     $      GO TO 90
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+(N-1)*(N-1)) ) + PDW
     $                         + (N-1)*(N-1) - 1 )
         IF ( N.GT.2 )
     $      CALL DLACPY( 'Lower', N-2, N-2, DWORK(PDW), N-2, A(3,1),
     $                   LDA )
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PZ), N, QG(1,2), LDQG, ZERO, DWORK(PDW), N )
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PDW), N, V1, LDV1, ZERO, QG(1,2), LDQG )
C
      ELSE IF ( WANTU .AND. WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: N*N + 7*N.
C
         PDW = PBETA + N
         CALL MB03XP( 'Schur', 'Init', 'Init', N, ILO, N, A, LDA, T,
     $                LDT, V1, LDV1, U1, LDU1, WR, WI, DWORK(PBETA),
     $                DWORK(PDW+(N-1)*(N-1)), LDWORK-PDW-(N-1)*(N-1)+1,
     $                INFO )
         IF ( INFO.NE.0 )
     $      GO TO 90
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+(N-1)*(N-1)) ) + PDW
     $                         + (N-1)*(N-1) - 1 )
         IF ( N.GT.1 )
     $      CALL DLACPY( 'Lower', N-1, N-1, DWORK(PDW), N-1, T(2,1),
     $                   LDT )
         IF ( N.GT.2 )
     $      CALL DLACPY( 'Lower', N-2, N-2, V2(3,1), LDV2, A(3,1), LDA )
         CALL DGEMM( 'Transpose', 'No Transpose', N, N, N, ONE,
     $                   U1, LDU1, QG(1,2), LDQG, ZERO, DWORK(PDW), N )
         CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $               DWORK(PDW), N, V1, LDV1, ZERO, QG(1,2), LDQG )
      END IF
C
   90 CONTINUE
C
C     Compute square roots of eigenvalues and rescale.
C
      DO 100 I = INFO + 1, N
         TEMPR = WR(I)
         TEMPI = WI(I)
         TEMP  = DWORK(PBETA + I - 1)
         IF ( TEMP.GT.ZERO )
     $      TEMPR = -TEMPR
         TEMP  = ABS( TEMP )
         IF ( TEMPI.EQ.ZERO ) THEN
            IF ( TEMPR.LT.ZERO ) THEN
               WR(I) = ZERO
               WI(I) =  SQRT( TEMP ) * SQRT( -TEMPR )
            ELSE
               WR(I) = -SQRT( TEMP ) * SQRT(  TEMPR )
               WI(I) = ZERO
            END IF
         ELSE
            CALL MA01AD( TEMPR, TEMPI, WR(I), WI(I) )
            WR(I) = -WR(I) * SQRT( TEMP )
            IF ( TEMP.GT.0 ) THEN
               WI(I) = WI(I) * SQRT( TEMP )
            ELSE
               WI(I) = ZERO
            END IF
         END IF
  100 CONTINUE
C
      IF ( SCALEH ) THEN
C
C        Undo scaling.
C
         CALL DLASCL( 'Hessenberg', 0, 0, CSCALE, HNRM, N, N, A, LDA,
     $                IERR )
         CALL DLASCL( 'Upper', 0, 0, CSCALE, HNRM, N, N, T, LDT, IERR )
         If ( WANTG )
     $      CALL DLASCL( 'General', 0, 0, CSCALE, HNRM, N, N, QG(1,2),
     $                   LDQG, IERR )
         CALL DLASCL( 'General', 0, 0, CSCALE, HNRM, N, 1, WR, N, IERR )
         CALL DLASCL( 'General', 0, 0, CSCALE, HNRM, N, 1, WI, N, IERR )
      END IF
C
      IF ( INFO.NE.0 )
     $   RETURN
C
C     -----------------------------------------------
C     Step 3:  Compute orthogonal symplectic factors.
C     -----------------------------------------------
C
C     Fix CSL and CSR for MB04QB.
C
      IF ( WANTU )
     $   CALL DSCAL( N, -ONE, DWORK(PCSL+1), 2 )
      IF ( WANTV )
     $   CALL DSCAL( N-1, -ONE, DWORK(PCSR+1), 2 )
      ILO1 = MIN( N, ILO + 1 )
C
      IF ( WANTU .AND. .NOT.WANTV .AND. .NOT.WANTG ) THEN
C
C        Workspace requirements: 7*N.
C
         PDW = PTAUR
         CALL DCOPY( N, T(1,1), LDT+1, DWORK(PDW), 1 )
         CALL DLACPY( 'Lower', N, N, U2, LDU2, T, LDT )
         CALL DLASET( 'All', N, N, ZERO, ZERO, U2, LDU2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Columnwise', N-ILO+1, N, N-ILO+1,
     $                QG(ILO,ILO), LDQG, T(ILO,ILO), LDT, U1(ILO,1),
     $                LDU1, U2(ILO,1), LDU2, DWORK(PCSL+2*ILO-2),
     $                DWORK(PTAUL+ILO-1), DWORK(PDW+N), LDWORK-PDW-N+1,
     $                IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+N) ) + PDW + N - 1 )
         CALL DCOPY( N, DWORK(PDW), 1, T(1,1), LDT+1 )
         IF ( N.GT.1 )
     $      CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, T(2,1), LDT )
C
      ELSE IF ( .NOT.WANTU .AND. WANTV .AND. .NOT.WANTG ) THEN
C
C        Workspace requirements: 7*N.
C
         PDW = PTAUR + N
         CALL DLASET( 'All', N, N, ZERO, ZERO, V2, LDV2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Rowwise', MAX(0,N-ILO), N,
     $                MAX(0,N-ILO), QG(ILO1,ILO), LDQG, QG(ILO,ILO1),
     $                LDQG, V1(ILO1,1), LDV1, V2(ILO1,1), LDV2,
     $                DWORK(PCSR+2*ILO-2), DWORK(PTAUR+ILO-1),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
      ELSE IF ( WANTU .AND. WANTV .AND. .NOT.WANTG ) THEN
C
C        Workspace requirements: 8*N.
C
         PDW = PTAUR + N
         CALL DCOPY( N, T(1,1), LDT+1, DWORK(PDW), 1 )
         CALL DLACPY( 'Lower', N, N, V2, LDV2, T, LDT )
         CALL DLASET( 'All', N, N, ZERO, ZERO, V2, LDV2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Rowwise', MAX(0,N-ILO), N,
     $                MAX(0,N-ILO), QG(ILO1,ILO), LDQG, U2(ILO,ILO1),
     $                LDU2, V1(ILO1,1), LDV1, V2(ILO1,1), LDV2,
     $                DWORK(PCSR+2*ILO-2), DWORK(PTAUR+ILO-1),
     $                DWORK(PDW+N), LDWORK-PDW-N+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+N) ) + PDW + N - 1 )
C
         CALL DLACPY( 'Lower', N, N, U2, LDU2, QG, LDQG )
         CALL DLASET( 'All', N, N, ZERO, ZERO, U2, LDU2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Columnwise', N-ILO+1, N, N-ILO+1,
     $                T(ILO,ILO), LDT, QG(ILO,ILO), LDQG, U1(ILO,1),
     $                LDU1, U2(ILO,1), LDU2, DWORK(PCSL+2*ILO-2),
     $                DWORK(PTAUL+ILO-1), DWORK(PDW+N), LDWORK-PDW-N+1,
     $                IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+N) ) + PDW + N - 1 )
         CALL DCOPY( N, DWORK(PDW), 1, T(1,1), LDT+1 )
         IF ( N.GT.1 )
     $      CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, T(2,1), LDT )
C
      ELSE IF ( WANTU .AND. .NOT.WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 6*N + N*N.
C
         PQ  = PTAUR
         PDW = PQ + N*N
         CALL DLACPY( 'Lower', N, N, U2, LDU2, DWORK(PQ), N )
         CALL DLASET( 'All', N, N, ZERO, ZERO, U2, LDU2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Columnwise', N-ILO+1, N, N-ILO+1,
     $                T(ILO,ILO), LDT, DWORK(PQ+(ILO-1)*(N+1)), N,
     $                U1(ILO,1), LDU1, U2(ILO,1), LDU2,
     $                DWORK(PCSL+2*ILO-2), DWORK(PTAUL+ILO-1),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         IF ( N.GT.1 )
     $      CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, T(2,1), LDT )
C
      ELSE IF ( .NOT.WANTU .AND. WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 7*N + N*N.
C
         PQ  = PTAUR+N
         PDW = PQ + N*N
         CALL DLACPY( 'Upper', N, N, V2, LDV2, DWORK(PQ), N )
         CALL DLASET( 'All', N, N, ZERO, ZERO, V2, LDV2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Rowwise', MAX(0,N-ILO), N,
     $                MAX(0,N-ILO), A(ILO1,ILO), LDA,
     $                DWORK(PQ+ILO*N+ILO-1), N, V1(ILO1,1), LDV1,
     $                V2(ILO1,1), LDV2, DWORK(PCSR+2*ILO-2),
     $                DWORK(PTAUR+ILO-1), DWORK(PDW+N),
     $                LDWORK-PDW-N+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW+N) ) + PDW + N - 1 )
         IF ( N.GT.2 )
     $      CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
C
      ELSE IF ( WANTU .AND. WANTV .AND. WANTG ) THEN
C
C        Workspace requirements: 6*N + N*N.
C
         PDW = PTAUR + N
         CALL DLASET( 'All', N, N, ZERO, ZERO, V2, LDV2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Rowwise', MAX(0,N-ILO), N,
     $                MAX(0,N-ILO), A(ILO1,ILO), LDA, U2(ILO,ILO1),
     $                LDU2, V1(ILO1,1), LDV1, V2(ILO1,1), LDV2,
     $                DWORK(PCSR+2*ILO-2), DWORK(PTAUR+ILO-1),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
         PQ  = PTAUR
         PDW = PQ + N*N
         CALL DLACPY( 'Lower', N, N, U2, LDU2, DWORK(PQ), N )
         CALL DLASET( 'All', N, N, ZERO, ZERO, U2, LDU2 )
         CALL MB04QB( 'No Transpose', 'No Transpose', 'No Transpose',
     $                'Columnwise', 'Columnwise', N-ILO+1, N, N-ILO+1,
     $                T(ILO,ILO), LDT, DWORK(PQ+(ILO-1)*(N+1)), N,
     $                U1(ILO,1), LDU1, U2(ILO,1), LDU2,
     $                DWORK(PCSL+2*ILO-2), DWORK(PTAUL+ILO-1),
     $                DWORK(PDW), LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
         IF ( N.GT.2 )
     $      CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, A(3,1), LDA )
         IF ( N.GT.1 )
     $      CALL DLASET( 'Lower', N-1, N-1, ZERO, ZERO, T(2,1), LDT )
      END IF
C
      DWORK(1) = DBLE( WRKOPT )
      RETURN
C *** Last line of MB03XD ***
      END
