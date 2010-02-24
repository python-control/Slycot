      SUBROUTINE MB03WD( JOB, COMPZ, N, P, ILO, IHI, ILOZ, IHIZ, H,
     $                   LDH1, LDH2, Z, LDZ1, LDZ2, WR, WI, DWORK,
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
C     To compute the Schur decomposition and the eigenvalues of a
C     product of matrices, H = H_1*H_2*...*H_p, with H_1 an upper
C     Hessenberg matrix and H_2, ..., H_p upper triangular matrices,
C     without evaluating the product. Specifically, the matrices Z_i
C     are computed, such that
C
C             Z_1' * H_1 * Z_2 = T_1,
C             Z_2' * H_2 * Z_3 = T_2,
C                    ...
C             Z_p' * H_p * Z_1 = T_p,
C
C     where T_1 is in real Schur form, and T_2, ..., T_p are upper
C     triangular.
C
C     The routine works primarily with the Hessenberg and triangular
C     submatrices in rows and columns ILO to IHI, but optionally applies
C     the transformations to all the rows and columns of the matrices
C     H_i, i = 1,...,p. The transformations can be optionally
C     accumulated.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates whether the user wishes to compute the full
C             Schur form or the eigenvalues only, as follows:
C             = 'E':  Compute the eigenvalues only;
C             = 'S':  Compute the factors T_1, ..., T_p of the full
C                     Schur form, T = T_1*T_2*...*T_p.
C
C     COMPZ   CHARACTER*1
C             Indicates whether or not the user wishes to accumulate
C             the matrices Z_1, ..., Z_p, as follows:
C             = 'N':  The matrices Z_1, ..., Z_p are not required;
C             = 'I':  Z_i is initialized to the unit matrix and the
C                     orthogonal transformation matrix Z_i is returned,
C                     i = 1, ..., p;
C             = 'V':  Z_i must contain an orthogonal matrix Q_i on
C                     entry, and the product Q_i*Z_i is returned,
C                     i = 1, ..., p.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix H.  N >= 0.
C
C     P       (input) INTEGER
C             The number of matrices in the product H_1*H_2*...*H_p.
C             P >= 1.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that all matrices H_j, j = 2, ..., p, are
C             already upper triangular in rows and columns 1:ILO-1 and
C             IHI+1:N, and H_1 is upper quasi-triangular in rows and
C             columns 1:ILO-1 and IHI+1:N, with H_1(ILO,ILO-1) = 0
C             (unless ILO = 1), and H_1(IHI+1,IHI) = 0 (unless IHI = N).
C             The routine works primarily with the Hessenberg submatrix
C             in rows and columns ILO to IHI, but applies the
C             transformations to all the rows and columns of the
C             matrices H_i, i = 1,...,p, if JOB = 'S'.
C             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N.
C
C     ILOZ    (input) INTEGER
C     IHIZ    (input) INTEGER
C             Specify the rows of Z to which the transformations must be
C             applied if COMPZ = 'I' or COMPZ = 'V'.
C             1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
C
C     H       (input/output) DOUBLE PRECISION array, dimension
C             (LDH1,LDH2,P)
C             On entry, the leading N-by-N part of H(*,*,1) must contain
C             the upper Hessenberg matrix H_1 and the leading N-by-N
C             part of H(*,*,j) for j > 1 must contain the upper
C             triangular matrix H_j, j = 2, ..., p.
C             On exit, if JOB = 'S', the leading N-by-N part of H(*,*,1)
C             is upper quasi-triangular in rows and columns ILO:IHI,
C             with any 2-by-2 diagonal blocks corresponding to a pair of
C             complex conjugated eigenvalues, and the leading N-by-N
C             part of H(*,*,j) for j > 1 contains the resulting upper
C             triangular matrix T_j.
C             If JOB = 'E', the contents of H are unspecified on exit.
C
C     LDH1    INTEGER
C             The first leading dimension of the array H.
C             LDH1 >= max(1,N).
C
C     LDH2    INTEGER
C             The second leading dimension of the array H.
C             LDH2 >= max(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension
C             (LDZ1,LDZ2,P)
C             On entry, if COMPZ = 'V', the leading N-by-N-by-P part of
C             this array must contain the current matrix Q of
C             transformations accumulated by SLICOT Library routine
C             MB03VY.
C             If COMPZ = 'I', Z need not be set on entry.
C             On exit, if COMPZ = 'V', or COMPZ = 'I', the leading
C             N-by-N-by-P part of this array contains the transformation
C             matrices which produced the Schur form; the
C             transformations are applied only to the submatrices
C             Z_j(ILOZ:IHIZ,ILO:IHI), j = 1, ..., P.
C             If COMPZ = 'N', Z is not referenced.
C
C     LDZ1    INTEGER
C             The first leading dimension of the array Z.
C             LDZ1 >= 1,        if COMPZ = 'N';
C             LDZ1 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'.
C
C     LDZ2    INTEGER
C             The second leading dimension of the array Z.
C             LDZ2 >= 1,        if COMPZ = 'N';
C             LDZ2 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             The real and imaginary parts, respectively, of the
C             computed eigenvalues ILO to IHI are stored in the
C             corresponding elements of WR and WI. If two eigenvalues
C             are computed as a complex conjugate pair, they are stored
C             in consecutive elements of WR and WI, say the i-th and
C             (i+1)th, with WI(i) > 0 and WI(i+1) < 0. If JOB = 'S', the
C             eigenvalues are stored in the same order as on the
C             diagonal of the Schur form returned in H.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION work array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= IHI-ILO+P-1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, ILO <= i <= IHI, the QR algorithm
C                   failed to compute all the eigenvalues ILO to IHI
C                   in a total of 30*(IHI-ILO+1) iterations;
C                   the elements i+1:IHI of WR and WI contain those
C                   eigenvalues which have been successfully computed.
C
C     METHOD
C
C     A refined version of the QR algorithm proposed in [1] and [2] is
C     used. The elements of the subdiagonal, diagonal, and the first
C     supradiagonal of current principal submatrix of H are computed
C     in the process.
C
C     REFERENCES
C
C     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P.
C         The periodic Schur decomposition: algorithms and applications.
C         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42,
C         1992.
C
C     [2] Sreedhar, J. and Van Dooren, P.
C         Periodic Schur form and some matrix equations.
C         Proc. of the Symposium on the Mathematical Theory of Networks
C         and Systems (MTNS'93), Regensburg, Germany (U. Helmke,
C         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable.
C
C     FURTHER COMMENTS
C
C     Note that for P = 1, the LAPACK Library routine DHSEQR could be
C     more efficient on some computer architectures than this routine,
C     because DHSEQR uses a block multishift QR algorithm.
C     When P is large and JOB = 'S', it could be more efficient to
C     compute the product matrix H, and use the LAPACK Library routines.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga,
C     German Aerospace Center, DLR Oberpfaffenhofen, February 1999.
C     Partly based on the routine PSHQR by A. Varga
C     (DLR Oberpfaffenhofen), January 22, 1996.
C
C     REVISIONS
C
C     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest.
C
C     KEYWORDS
C
C     Eigenvalue, eigenvalue decomposition, Hessenberg form,
C     orthogonal transformation, periodic systems, (periodic) Schur
C     form, real Schur form, similarity transformation, triangular form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
      DOUBLE PRECISION  DAT1, DAT2
      PARAMETER         ( DAT1 = 0.75D+0, DAT2 = -0.4375D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER         COMPZ, JOB
      INTEGER           IHI, IHIZ, ILO, ILOZ, INFO, LDH1, LDH2, LDWORK,
     $                  LDZ1, LDZ2, N, P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK( * ), H( LDH1, LDH2, * ), WI( * ),
     $                  WR( * ), Z( LDZ1, LDZ2, * )
C     ..
C     .. Local Scalars ..
      LOGICAL           INITZ, WANTT, WANTZ
      INTEGER           I, I1, I2, ITN, ITS, J, JMAX, JMIN, K, L, M,
     $                  NH, NR, NROW, NZ
      DOUBLE PRECISION  AVE, CS, DISC, H11, H12, H21, H22, H33, H33S,
     $                  H43H34, H44, H44S, HH10, HH11, HH12, HH21, HH22,
     $                  HP00, HP01, HP02, HP11, HP12, HP22, OVFL, S,
     $                  SMLNUM, SN, TAU, TST1, ULP, UNFL, V1, V2, V3
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION  V( 3 )
C     ..
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANHS, DLANTR
      EXTERNAL          DLAMCH, DLANHS, DLANTR, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DLABAD, DLANV2, DLARFG, DLARFX, DLARTG,
     $                  DLASET, DROT, MB04PY, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SIGN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      WANTT = LSAME( JOB,   'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = LSAME( COMPZ, 'V' ) .OR. INITZ
      INFO = 0
      IF( .NOT. ( WANTT .OR. LSAME( JOB, 'E' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( WANTZ .OR. LSAME( COMPZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.1 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -6
      ELSE IF( ILOZ.LT.1 .OR. ILOZ.GT.ILO ) THEN
         INFO = -7
      ELSE IF( IHIZ.LT.IHI .OR. IHIZ.GT.N ) THEN
         INFO = -8
      ELSE IF( LDH1.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDH2.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDZ1.LT.1 .OR. ( WANTZ .AND. LDZ1.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LDZ2.LT.1 .OR. ( WANTZ .AND. LDZ2.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDWORK.LT.IHI - ILO + P - 1 ) THEN
         INFO = -18
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( ILO.GT.1 ) THEN
            IF( H( ILO, ILO-1, 1 ).NE.ZERO )
     $         INFO = -5
         ELSE IF( IHI.LT.N ) THEN
            IF( H( IHI+1, IHI, 1 ).NE.ZERO )
     $         INFO = -6
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03WD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Initialize Z, if necessary.
C
      IF( INITZ ) THEN
C
         DO 10 J = 1, P
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z( 1, 1, J ), LDZ1 )
   10    CONTINUE
C
      END IF
C
      NH = IHI - ILO + 1
C
      IF( NH.EQ.1 ) THEN
         HP00 = ONE
C
         DO 20 J = 1, P
            HP00 = HP00 * H( ILO, ILO, J )
   20    CONTINUE
C
         WR( ILO ) = HP00
         WI( ILO ) = ZERO
         RETURN
      END IF
C
C     Set machine-dependent constants for the stopping criterion.
C     If norm(H) <= sqrt(OVFL), overflow should not occur.
C
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( DBLE( NH ) / ULP )
C
C     Set the elements in rows and columns ILO to IHI to zero below the
C     first subdiagonal in H(*,*,1) and below the first diagonal in
C     H(*,*,j), j >= 2. In the same loop, compute and store in
C     DWORK(NH:NH+P-2) the 1-norms of the matrices H_2, ..., H_p, to be
C     used later.
C
      I = NH
      S = ULP * DBLE( N )
      IF( NH.GT.2 )
     $   CALL DLASET( 'Lower', NH-2, NH-2, ZERO, ZERO,
     $                H( ILO+2, ILO, 1 ), LDH1 )
C
      DO 30 J = 2, P
         CALL DLASET( 'Lower', NH-1, NH-1, ZERO, ZERO,
     $                H( ILO+1, ILO, J ), LDH1 )
         DWORK( I ) = S * DLANTR( '1-norm', 'Upper', 'NonUnit', NH, NH,
     $                            H( ILO, ILO, J ), LDH1, DWORK )
         I = I + 1
   30 CONTINUE
C
C     I1 and I2 are the indices of the first row and last column of H
C     to which transformations must be applied. If eigenvalues only are
C     being computed, I1 and I2 are set inside the main loop.
C
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
C
      IF( WANTZ )
     $   NZ = IHIZ - ILOZ + 1
C
C     ITN is the total number of QR iterations allowed.
C
      ITN = 30*NH
C
C     The main loop begins here. I is the loop index and decreases from
C     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
C     with the active submatrix in rows and columns L to I.
C     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
C     H(L,L-1) is negligible so that the matrix splits.
C
      I = IHI
C
   40 CONTINUE
      L = ILO
C
C     Perform QR iterations on rows and columns ILO to I until a
C     submatrix of order 1 or 2 splits off at the bottom because a
C     subdiagonal element has become negligible.
C
C     Let T = H_2*...*H_p, and H = H_1*T. Part of the currently
C     free locations of WR and WI are temporarily used as workspace.
C
C     WR(L:I):      the current diagonal elements of h = H(L:I,L:I);
C     WI(L+1:I):    the current elements of the first subdiagonal of h;
C     DWORK(NH-I+L:NH-1): the current elements of the first
C                   supradiagonal of h.
C
      DO 160 ITS = 0, ITN
C
C        Initialization: compute H(I,I) (and H(I,I-1) if I > L).
C
         HP22 = ONE
         IF( I.GT.L ) THEN
            HP12 = ZERO
            HP11 = ONE
C
            DO 50 J = 2, P
               HP22 = HP22*H( I,   I,   J )
               HP12 = HP11*H( I-1, I,   J ) + HP12*H( I, I, J )
               HP11 = HP11*H( I-1, I-1, J )
   50       CONTINUE
C
            HH21 = H( I, I-1, 1 )*HP11
            HH22 = H( I, I-1, 1 )*HP12 + H( I, I, 1 )*HP22
C
            WR( I ) = HH22
            WI( I ) = HH21
         ELSE
C
            DO 60 J = 1, P
               HP22 = HP22*H( I, I, J )
   60       CONTINUE
C
            WR( I ) = HP22
         END IF
C
C        Look for a single small subdiagonal element.
C        The loop also computes the needed current elements of the
C        diagonal and the first two supradiagonals of T, as well as
C        the current elements of the central tridiagonal of H.
C
         DO 80 K = I, L + 1, -1
C
C           Evaluate H(K-1,K-1), H(K-1,K) (and H(K-1,K-2) if K > L+1).
C
            HP00 = ONE
            HP01 = ZERO
            IF( K.GT.L+1 ) THEN
               HP02 = ZERO
C
               DO 70 J = 2, P
                  HP02 = HP00*H( K-2, K,   J ) + HP01*H( K-1, K,   J )
     $                                         + HP02*H( K,   K,   J )
                  HP01 = HP00*H( K-2, K-1, J ) + HP01*H( K-1, K-1, J )
                  HP00 = HP00*H( K-2, K-2, J )
   70          CONTINUE
C
               HH10 = H( K-1, K-2, 1 )*HP00
               HH11 = H( K-1, K-2, 1 )*HP01 + H( K-1, K-1, 1 )*HP11
               HH12 = H( K-1, K-2, 1 )*HP02 + H( K-1, K-1, 1 )*HP12
     $                                      + H( K-1, K,   1 )*HP22
               WI( K-1 ) = HH10
            ELSE
               HH10 = ZERO
               HH11 = H( K-1, K-1, 1 )*HP11
               HH12 = H( K-1, K-1, 1 )*HP12 + H( K-1, K, 1 )*HP22
            END IF
            WR( K-1 ) = HH11
            DWORK( NH-I+K-1) = HH12
C
C           Test for a negligible subdiagonal element.
C
            TST1 = ABS( HH11 ) + ABS( HH22 )
            IF( TST1.EQ.ZERO )
     $         TST1 = DLANHS( '1-norm', I-L+1, H( L, L, 1 ), LDH1,
     $                        DWORK )
            IF( ABS( HH21 ).LE.MAX( ULP*TST1, SMLNUM ) )
     $         GO TO 90
C
C           Update the values for the next cycle.
C
            HP22 = HP11
            HP11 = HP00
            HP12 = HP01
            HH22 = HH11
            HH21 = HH10
   80    CONTINUE
C
   90    CONTINUE
         L = K
C
         IF( L.GT.ILO ) THEN
C
C           H(L,L-1) is negligible.
C
            IF( WANTT ) THEN
C
C              If H(L,L-1,1) is also negligible, set it to 0; otherwise,
C              annihilate the subdiagonal elements bottom-up, and
C              restore the triangular form of H(*,*,j). Since H(L,L-1)
C              is negligible, the second case can only appear when the
C              product of H(L-1,L-1,j), j >= 2, is negligible.
C
               TST1 = ABS( H( L-1, L-1, 1 ) ) + ABS( H( L, L, 1 ) )
               IF( TST1.EQ.ZERO )
     $            TST1 = DLANHS( '1-norm', I-L+1, H( L, L, 1 ), LDH1,
     $                           DWORK )
               IF( ABS( H( L, L-1, 1 ) ).GT.MAX( ULP*TST1, SMLNUM ) )
     $            THEN
C
                  DO 110 K = I, L, -1
C
                     DO 100 J = 1, P - 1
C
C                       Compute G to annihilate from the right the
C                       (K,K-1) element of the matrix H_j.
C
                        V( 1 ) = H( K, K-1, J )
                        CALL DLARFG( 2, H( K, K, J ), V, 1, TAU )
                        H( K, K-1, J ) = ZERO
                        V( 2 ) = ONE
C
C                       Apply G from the right to transform the columns
C                       of the matrix H_j in rows I1 to K-1.
C
                        CALL DLARFX( 'Right', K-I1, 2, V, TAU,
     $                               H( I1, K-1, J ), LDH1, DWORK )
C
C                       Apply G from the left to transform the rows of
C                       the matrix H_(j+1) in columns K-1 to I2.
C
                        CALL DLARFX( 'Left', 2, I2-K+2, V, TAU,
     $                               H( K-1, K-1, J+1 ), LDH1, DWORK )
C
                        IF( WANTZ ) THEN
C
C                          Accumulate transformations in the matrix
C                          Z_(j+1).
C
                           CALL DLARFX( 'Right', NZ, 2, V, TAU,
     $                                  Z( ILOZ, K-1, J+1 ), LDZ1,
     $                                  DWORK )
                        END IF
  100                CONTINUE
C
                     IF( K.LT.I ) THEN
C
C                       Compute G to annihilate from the right the
C                       (K+1,K) element of the matrix H_p.
C
                        V( 1 ) = H( K+1, K, P )
                        CALL DLARFG( 2, H( K+1, K+1, P ), V, 1, TAU )
                        H( K+1, K, P ) = ZERO
                        V( 2 ) = ONE
C
C                       Apply G from the right to transform the columns
C                       of the matrix H_p in rows I1 to K.
C
                        CALL DLARFX( 'Right', K-I1+1, 2, V, TAU,
     $                               H( I1, K, P ), LDH1, DWORK )
C
C                       Apply G from the left to transform the rows of
C                       the matrix H_1 in columns K to I2.
C
                        CALL DLARFX( 'Left', 2, I2-K+1, V, TAU,
     $                               H( K, K, 1 ), LDH1, DWORK )
C
                        IF( WANTZ ) THEN
C
C                          Accumulate transformations in the matrix Z_1.
C
                           CALL DLARFX( 'Right', NZ, 2, V, TAU,
     $                                  Z( ILOZ, K, 1 ), LDZ1, DWORK )
                        END IF
                     END IF
  110             CONTINUE
C
                  H( L, L-1, P ) = ZERO
               END IF
               H( L, L-1, 1 ) = ZERO
            END IF
         END IF
C
C        Exit from loop if a submatrix of order 1 or 2 has split off.
C
         IF( L.GE.I-1 )
     $      GO TO 170
C
C        Now the active submatrix is in rows and columns L to I. If
C        eigenvalues only are being computed, only the active submatrix
C        need be transformed.
C
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
C
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
C
C           Exceptional shift.
C
            S   = ABS( WI( I ) ) + ABS( WI( I-1 ) )
            H44 = DAT1*S + WR( I )
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
C
C           Prepare to use Francis' double shift (i.e., second degree
C           generalized Rayleigh quotient).
C
            H44 = WR( I )
            H33 = WR( I-1 )
            H43H34 = WI( I )*DWORK( NH-1 )
            DISC = ( H33 - H44 )*HALF
            DISC = DISC*DISC + H43H34
            IF( DISC.GT.ZERO ) THEN
C
C              Real roots: use Wilkinson's shift twice.
C
               DISC = SQRT( DISC )
               AVE  = HALF*( H33 + H44 )
               IF( ABS( H33 )-ABS( H44 ).GT.ZERO ) THEN
                  H33 = H33*H44 - H43H34
                  H44 = H33 / ( SIGN( DISC, AVE ) + AVE )
               ELSE
                  H44 = SIGN( DISC, AVE ) + AVE
               END IF
               H33 = H44
               H43H34 = ZERO
            END IF
         END IF
C
C        Look for two consecutive small subdiagonal elements.
C
         DO 120 M = I - 2, L, -1
C
C           Determine the effect of starting the double-shift QR
C           iteration at row M, and see if this would make H(M,M-1)
C           negligible.
C
            H11  = WR( M )
            H12  = DWORK( NH-I+M )
            H21  = WI( M+1 )
            H22  = WR( M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1   = ( H33S*H44S - H43H34 ) / H21 + H12
            V2   = H22 - H11 - H33S - H44S
            V3   = WI( M+2 )
            S    = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
            V1   = V1 / S
            V2   = V2 / S
            V3   = V3 / S
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            IF( M.EQ.L )
     $         GO TO 130
            TST1 = ABS( V1 )*( ABS( WR( M-1 ) ) +
     $                         ABS( H11 ) + ABS( H22 ) )
            IF( ABS( WI( M ) )*( ABS( V2 ) + ABS( V3 ) ).LE.ULP*TST1 )
     $         GO TO 130
  120    CONTINUE
C
  130    CONTINUE
C
C        Double-shift QR step.
C
         DO 150 K = M, I - 1
C
C           The first iteration of this loop determines a reflection G
C           from the vector V and applies it from left and right to H,
C           thus creating a nonzero bulge below the subdiagonal.
C
C           Each subsequent iteration determines a reflection G to
C           restore the Hessenberg form in the (K-1)th column, and thus
C           chases the bulge one step toward the bottom of the active
C           submatrix. NR is the order of G.
C
            NR   = MIN( 3, I-K+1 )
            NROW = MIN( K+NR, I ) - I1 + 1
            IF( K.GT.M )
     $         CALL DCOPY( NR, H( K, K-1, 1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, TAU )
            IF( K.GT.M ) THEN
               H( K,   K-1, 1 ) = V( 1 )
               H( K+1, K-1, 1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1, 1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
               H( K, K-1, 1 ) = -H( K, K-1, 1 )
            END IF
C
C           Apply G from the left to transform the rows of the matrix
C           H_1 in columns K to I2.
C
            CALL MB04PY( 'Left', NR, I2-K+1, V( 2 ), TAU, H( K, K, 1 ),
     $                   LDH1, DWORK )
C
C           Apply G from the right to transform the columns of the
C           matrix H_p in rows I1 to min(K+NR,I).
C
            CALL MB04PY( 'Right', NROW, NR, V( 2 ), TAU, H( I1, K, P ),
     $                   LDH1, DWORK )
C
            IF( WANTZ ) THEN
C
C              Accumulate transformations in the matrix Z_1.
C
               CALL MB04PY( 'Right', NZ, NR, V( 2 ), TAU,
     $                      Z( ILOZ, K, 1 ), LDZ1, DWORK )
            END IF
C
            DO 140 J = P, 2, -1
C
C              Apply G1 (and G2, if NR = 3) from the left to transform
C              the NR-by-NR submatrix of H_j in position (K,K) to upper
C              triangular form.
C
C              Compute G1.
C
               CALL DCOPY( NR-1, H( K+1, K, J ), 1, V, 1 )
               CALL DLARFG( NR, H( K, K, J ), V, 1, TAU )
               H( K+1, K, J ) = ZERO
               IF( NR.EQ.3 )
     $            H( K+2, K, J ) = ZERO
C
C              Apply G1 from the left to transform the rows of the
C              matrix H_j in columns K+1 to I2.
C
               CALL MB04PY( 'Left', NR, I2-K, V, TAU, H( K, K+1, J ),
     $                      LDH1, DWORK )
C
C              Apply G1 from the right to transform the columns of the
C              matrix H_(j-1) in rows I1 to min(K+NR,I).
C
               CALL MB04PY( 'Right', NROW, NR, V, TAU, H( I1, K, J-1 ),
     $                      LDH1, DWORK )
C
               IF( WANTZ ) THEN
C
C                 Accumulate transformations in the matrix Z_j.
C
                  CALL MB04PY( 'Right', NZ, NR, V, TAU, Z( ILOZ, K, J ),
     $                         LDZ1, DWORK )
               END IF
C
               IF( NR.EQ.3 ) THEN
C
C                 Compute G2.
C
                  V( 1 ) = H( K+2, K+1, J )
                  CALL DLARFG( 2, H( K+1, K+1, J ), V, 1, TAU )
                  H( K+2, K+1, J ) = ZERO
C
C                 Apply G2 from the left to transform the rows of the
C                 matrix H_j in columns K+2 to I2.
C
                  CALL MB04PY( 'Left', 2, I2-K-1, V, TAU,
     $                         H( K+1, K+2, J ), LDH1, DWORK )
C
C                 Apply G2 from the right to transform the columns of
C                 the matrix H_(j-1) in rows I1 to min(K+3,I).
C
                  CALL MB04PY( 'Right', NROW, 2, V, TAU,
     $                         H( I1, K+1, J-1 ), LDH1, DWORK )
C
                  IF( WANTZ ) THEN
C
C                    Accumulate transformations in the matrix Z_j.
C
                     CALL MB04PY( 'Right', NZ, 2, V, TAU,
     $                            Z( ILOZ, K+1, J ), LDZ1, DWORK )
                  END IF
               END IF
  140       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
C     Failure to converge in remaining number of iterations.
C
      INFO = I
      RETURN
C
  170 CONTINUE
C
      IF( L.EQ.I ) THEN
C
C        H(I,I-1,1) is negligible: one eigenvalue has converged.
C        Note that WR(I) has already been set.
C
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
C
C        H(I-1,I-2,1) is negligible: a pair of eigenvalues have
C        converged.
C
C        Transform the 2-by-2 submatrix of H_1*H_2*...*H_p in position
C        (I-1,I-1) to standard Schur form, and compute and store its
C        eigenvalues. If the Schur form is not required, then the
C        previously stored values of a similar submatrix are used.
C        For real eigenvalues, a Givens transformation is used to
C        triangularize the submatrix.
C
         IF( WANTT ) THEN
            HP22 = ONE
            HP12 = ZERO
            HP11 = ONE
C
            DO 180 J = 2, P
               HP22 = HP22*H( I,   I,   J )
               HP12 = HP11*H( I-1, I,   J ) + HP12*H( I, I, J )
               HP11 = HP11*H( I-1, I-1, J )
  180       CONTINUE
C
            HH21 = H( I,   I-1, 1 )*HP11
            HH22 = H( I,   I-1, 1 )*HP12 + H( I,   I, 1 )*HP22
            HH11 = H( I-1, I-1, 1 )*HP11
            HH12 = H( I-1, I-1, 1 )*HP12 + H( I-1, I, 1 )*HP22
         ELSE
            HH11 = WR( I-1 )
            HH12 = DWORK( NH-1 )
            HH21 = WI( I )
            HH22 = WR( I )
         END IF
C
         CALL DLANV2( HH11, HH12, HH21, HH22, WR( I-1 ), WI( I-1 ),
     $                WR( I ), WI( I ), CS, SN )
C
         IF( WANTT ) THEN
C
C           Detect negligible diagonal elements in positions (I-1,I-1)
C           and (I,I) in H_j, J > 1.
C
            JMIN = 0
            JMAX = 0
C
            DO 190 J = 2, P
               IF( JMIN.EQ.0 ) THEN
                  IF( ABS( H( I-1, I-1, J ) ).LE.DWORK( NH+J-2 ) )
     $                JMIN = J
               END IF
               IF( ABS( H( I, I, J ) ).LE.DWORK( NH+J-2 ) ) JMAX = J
  190       CONTINUE
C
            IF( JMIN.NE.0 .AND. JMAX.NE.0 ) THEN
C
C              Choose the shorter path if zero elements in both
C              (I-1,I-1) and (I,I) positions are present.
C
               IF( JMIN-1.LE.P-JMAX+1 ) THEN
                  JMAX = 0
               ELSE
                  JMIN = 0
               END IF
            END IF
C
            IF( JMIN.NE.0 ) THEN
C
               DO 200 J = 1, JMIN - 1
C
C                 Compute G to annihilate from the right the (I,I-1)
C                 element of the matrix H_j.
C
                  V( 1 ) = H( I, I-1, J )
                  CALL DLARFG( 2, H( I, I, J ), V, 1, TAU )
                  H( I, I-1, J ) = ZERO
                  V( 2 ) = ONE
C
C                 Apply G from the right to transform the columns of the
C                 matrix H_j in rows I1 to I-1.
C
                  CALL DLARFX( 'Right', I-I1, 2, V, TAU,
     $                         H( I1, I-1, J ), LDH1, DWORK )
C
C                 Apply G from the left to transform the rows of the
C                 matrix H_(j+1) in columns I-1 to I2.
C
                  CALL DLARFX( 'Left', 2, I2-I+2, V, TAU,
     $                         H( I-1, I-1, J+1 ), LDH1, DWORK )
C
                  IF( WANTZ ) THEN
C
C                    Accumulate transformations in the matrix Z_(j+1).
C
                     CALL DLARFX( 'Right', NZ, 2, V, TAU,
     $                            Z( ILOZ, I-1, J+1 ), LDZ1, DWORK )
                  END IF
  200          CONTINUE
C
               H( I, I-1, JMIN ) = ZERO
C
            ELSE
               IF( JMAX.GT.0 .AND. WI( I-1 ).EQ.ZERO )
     $            CALL DLARTG( H( I-1, I-1, 1 ), H( I, I-1, 1 ), CS, SN,
     $                         TAU )
C
C              Apply the transformation to H.
C
               CALL DROT( I2-I+2, H( I-1, I-1, 1 ), LDH1,
     $                    H( I, I-1, 1 ), LDH1, CS, SN )
               CALL DROT( I-I1+1, H( I1, I-1, P ), 1, H( I1, I, P ), 1,
     $                    CS, SN )
               IF( WANTZ ) THEN
C
C                 Apply transformation to Z_1.
C
                  CALL DROT( NZ, Z( ILOZ, I-1, 1 ), 1, Z( ILOZ, I, 1 ),
     $                       1, CS, SN )
               END IF
C
               DO 210 J = P, MAX( 2, JMAX+1 ), -1
C
C                 Compute G1 to annihilate from the left the (I,I-1)
C                 element of the matrix H_j.
C
                  V( 1 ) = H( I, I-1, J )
                  CALL DLARFG( 2, H( I-1, I-1, J ), V, 1, TAU )
                  H( I, I-1, J ) = ZERO
C
C                 Apply G1 from the left to transform the rows of the
C                 matrix H_j in columns I to I2.
C
                  CALL MB04PY( 'Left', 2, I2-I+1, V, TAU,
     $                         H( I-1, I, J ), LDH1, DWORK )
C
C                 Apply G1 from the right to transform the columns of
C                 the matrix H_(j-1) in rows I1 to I.
C
                  CALL MB04PY( 'Right', I-I1+1, 2, V, TAU,
     $                         H( I1, I-1, J-1 ), LDH1, DWORK )
C
                  IF( WANTZ ) THEN
C
C                    Apply G1 to Z_j.
C
                     CALL MB04PY( 'Right', NZ, 2, V, TAU,
     $                            Z( ILOZ, I-1, J ), LDZ1, DWORK )
                  END IF
  210          CONTINUE
C
               IF( JMAX.GT.0 ) THEN
                  H( I, I-1, 1 )    = ZERO
                  H( I, I-1, JMAX ) = ZERO
               ELSE
                  IF( HH21.EQ.ZERO )
     $               H( I, I-1, 1 ) = ZERO
               END IF
            END IF
         END IF
      END IF
C
C     Decrement number of remaining iterations, and return to start of
C     the main loop with new value of I.
C
      ITN = ITN - ITS
      I = L - 1
      IF( I.GE.ILO )
     $   GO TO 40
C
      RETURN
C
C *** Last line of MB03WD ***
      END
