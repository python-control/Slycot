      SUBROUTINE MB03ZD( WHICH, METH, STAB, BALANC, ORTBAL, SELECT, N,
     $                   MM, ILO, SCALE, S, LDS, T, LDT, G, LDG, U1,
     $                   LDU1, U2, LDU2, V1, LDV1, V2, LDV2, M, WR, WI,
     $                   US, LDUS, UU, LDUU, LWORK, IWORK, DWORK,
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
C     To compute the stable and unstable invariant subspaces for a
C     Hamiltonian matrix with no eigenvalues on the imaginary axis,
C     using the output of the SLICOT Library routine MB03XD.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     WHICH   CHARACTER*1
C             Specifies the cluster of eigenvalues for which the
C             invariant subspaces are computed:
C             = 'A':  select all n eigenvalues;
C             = 'S':  select a cluster of eigenvalues specified by
C                     SELECT.
C
C     METH    CHARACTER*1
C             If WHICH = 'A' this parameter specifies the method to be
C             used for computing bases of the invariant subspaces:
C             = 'S':  compute the n-dimensional basis from a set of
C                     n vectors;
C             = 'L':  compute the n-dimensional basis from a set of
C                     2*n vectors.
C             When in doubt, use METH = 'S'. In some cases, METH = 'L'
C             may result in more accurately computed invariant
C             subspaces, see [1].
C
C     STAB    CHARACTER*1
C             Specifies the type of invariant subspaces to be computed:
C             = 'S':  compute the stable invariant subspace, i.e., the
C                     invariant subspace belonging to those selected
C                     eigenvalues that have negative real part;
C             = 'U':  compute the unstable invariant subspace, i.e.,
C                     the invariant subspace belonging to those
C                     selected eigenvalues that have positive real
C                     part;
C             = 'B':  compute both the stable and unstable invariant
C                     subspaces.
C
C     BALANC  CHARACTER*1
C             Specifies the type of inverse balancing transformation
C             required:
C             = 'N':  do nothing;
C             = 'P':  do inverse transformation for permutation only;
C             = 'S':  do inverse transformation for scaling only;
C             = 'B':  do inverse transformations for both permutation
C                     and scaling.
C             BALANC must be the same as the argument BALANC supplied to
C             MB03XD. Note that if the data is further post-processed,
C             e.g., for solving an algebraic Riccati equation, it is
C             recommended to delay inverse balancing (in particular the
C             scaling part) and apply it to the final result only,
C             see [2].
C
C     ORTBAL  CHARACTER*1
C             If BALANC <> 'N', this option specifies how inverse
C             balancing is applied to the computed invariant subspaces:
C             = 'B':  apply inverse balancing before orthogonal bases
C                     for the invariant subspaces are computed;
C             = 'A':  apply inverse balancing after orthogonal bases
C                     for the invariant subspaces have been computed;
C                     this may yield non-orthogonal bases if
C                     BALANC = 'S' or BALANC = 'B'.
C
C     SELECT  (input) LOGICAL array, dimension (N)
C             If WHICH = 'S', SELECT specifies the eigenvalues
C             corresponding to the positive and negative square
C             roots of the eigenvalues of S*T in the selected cluster.
C             To select a real eigenvalue w(j), SELECT(j) must be set
C             to .TRUE.. To select a complex conjugate pair of
C             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2
C             diagonal block, both SELECT(j) and SELECT(j+1) must be set
C             to .TRUE.; a complex conjugate pair of eigenvalues must be
C             either both included in the cluster or both excluded.
C             This array is not referenced if WHICH = 'A'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices S, T and G. N >= 0.
C
C     MM      (input) INTEGER
C             The number of columns in the arrays US and/or UU.
C             If WHICH = 'A' and METH = 'S',  MM >=   N;
C             if WHICH = 'A' and METH = 'L',  MM >= 2*N;
C             if WHICH = 'S',                 MM >=   M.
C             The minimal values above for MM give the numbers of
C             vectors to be used for computing a basis for the
C             invariant subspace(s).
C
C     ILO     (input) INTEGER
C             If BALANC <> 'N', then ILO is the integer returned by
C             MB03XD.  1 <= ILO <= N+1.
C
C     SCALE   (input) DOUBLE PRECISION array, dimension (N)
C             If BALANC <> 'N', the leading N elements of this array
C             must contain details of the permutation and scaling
C             factors, as returned by MB03XD.
C             This array is not referenced if BALANC = 'N'.
C
C     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix S in real Schur form.
C             On exit, the leading N-by-N part of this array is
C             overwritten.
C
C     LDS     INTEGER
C             The leading dimension of the array S.  LDS >= max(1,N).
C
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular matrix T.
C             On exit, the leading N-by-N part of this array is
C             overwritten.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= max(1,N).
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N)
C             On entry, if METH = 'L', the leading N-by-N part of this
C             array must contain a general matrix G.
C             On exit, if METH = 'L', the leading N-by-N part of this
C             array is overwritten.
C             This array is not referenced if METH = 'S'.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= 1.
C             LDG >= max(1,N) if METH = 'L'.
C
C     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On entry, the leading N-by-N part of this array must
C             contain the (1,1) block of an orthogonal symplectic
C             matrix U.
C             On exit, this array is overwritten.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.  LDU1 >= MAX(1,N).
C
C     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On entry, the leading N-by-N part of this array must
C             contain the (2,1) block of an orthogonal symplectic
C             matrix U.
C             On exit, this array is overwritten.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.  LDU2 >= MAX(1,N).
C
C     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N)
C             On entry, the leading N-by-N part of this array must
C             contain the (1,1) block of an orthogonal symplectic
C             matrix V.
C             On exit, this array is overwritten.
C
C     LDV1    INTEGER
C             The leading dimension of the array V1.  LDV1 >= MAX(1,N).
C
C     V2      (input/output) DOUBLE PRECISION array, dimension (LDV1,N)
C             On entry, the leading N-by-N part of this array must
C             contain the (2,1) block of an orthogonal symplectic
C             matrix V.
C             On exit, this array is overwritten.
C
C     LDV2    INTEGER
C             The leading dimension of the array V2.  LDV2 >= MAX(1,N).
C
C     M       (output) INTEGER
C             The number of selected eigenvalues.
C
C     WR      (output) DOUBLE PRECISION array, dimension (M)
C     WI      (output) DOUBLE PRECISION array, dimension (M)
C             On exit, the leading M elements of WR and WI contain the
C             real and imaginary parts, respectively, of the selected
C             eigenvalues that have nonpositive real part. Complex
C             conjugate pairs of eigenvalues with real part not equal
C             to zero will appear consecutively with the eigenvalue
C             having the positive imaginary part first. Note that, due
C             to roundoff errors, these numbers may differ from the
C             eigenvalues computed by MB03XD.
C
C     US      (output) DOUBLE PRECISION array, dimension (LDUS,MM)
C             On exit, if STAB = 'S' or STAB = 'B', the leading 2*N-by-M
C             part of this array contains a basis for the stable
C             invariant subspace belonging to the selected eigenvalues.
C             This basis is orthogonal unless ORTBAL = 'A'.
C
C     LDUS    INTEGER
C             The leading dimension of the array US.  LDUS >= 1.
C             If STAB = 'S' or STAB = 'B',  LDUS >= 2*N.
C
C     UU      (output) DOUBLE PRECISION array, dimension (LDUU,MM)
C             On exit, if STAB = 'U' or STAB = 'B', the leading 2*N-by-M
C             part of this array contains a basis for the unstable
C             invariant subspace belonging to the selected eigenvalues.
C             This basis is orthogonal unless ORTBAL = 'A'.
C
C     LDUU    INTEGER
C             The leading dimension of the array UU.  LDUU >= 1.
C             If STAB = 'U' or STAB = 'B',  LDUU >= 2*N.
C
C     Workspace
C
C     LWORK   LOGICAL array, dimension (2*N)
C             This array is only referenced if WHICH = 'A' and
C             METH = 'L'.
C
C     IWORK   INTEGER array, dimension (2*N),
C             This array is only referenced if WHICH = 'A' and
C             METH = 'L'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -35,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             If WHICH = 'S' or METH = 'S':
C                LDWORK >= MAX( 1, 4*M*M + MAX( 8*M, 4*N ) ).
C             If WHICH = 'A' and METH = 'L' and
C                ( STAB = 'U' or STAB = 'S' ):
C                LDWORK >= MAX( 1, 2*N*N + 2*N, 8*N ).
C             If WHICH = 'A' and METH = 'L' and STAB = 'B':
C                LDWORK >= 8*N + 1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  some of the selected eigenvalues are on or too close
C                   to the imaginary axis;
C             = 2:  reordering of the product S*T in routine MB03ZA
C                   failed because some eigenvalues are too close to
C                   separate;
C             = 3:  the QR algorithm failed to compute some Schur form
C                   in MB03ZA;
C             = 4:  reordering of the Hamiltonian Schur form in routine
C                   MB03TD failed because some eigenvalues are too close
C                   to separate.
C
C     METHOD
C
C     This is an implementation of Algorithm 1 in [1].
C
C     NUMERICAL ASPECTS
C
C     The method is strongly backward stable for an embedded
C     (skew-)Hamiltonian matrix, see [1]. Although good results have
C     been reported if the eigenvalues are not too close to the
C     imaginary axis, the method is not backward stable for the original
C     Hamiltonian matrix itself.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A new method for computing the stable invariant subspace of a
C         real Hamiltonian matrix, J. Comput. Appl. Math., 86,
C         pp. 17-43, 1997.
C
C     [2] Benner, P.
C         Symplectic balancing of Hamiltonian matrices.
C         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHASUB).
C
C     KEYWORDS
C
C     Hamiltonian matrix, invariant subspace.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         BALANC, METH, ORTBAL, STAB, WHICH
      INTEGER           ILO, INFO, LDG, LDS, LDT, LDU1, LDU2, LDUS,
     $                  LDUU, LDV1, LDV2, LDWORK, M, MM, N
C     .. Array Arguments ..
      LOGICAL           LWORK(*), SELECT(*)
      INTEGER           IWORK(*)
      DOUBLE PRECISION  DWORK(*), G(LDG,*), S(LDS,*), SCALE(*),
     $                  T(LDT,*), U1(LDU1,*), U2(LDU2,*), US(LDUS,*),
     $                  UU(LDUU,*), V1(LDV1,*), V2(LDV2,*), WI(*),
     $                  WR(*)
C     .. Local Scalars ..
      LOGICAL           LALL, LBAL, LBEF, LEXT, LUS, LUU, PAIR
      INTEGER           I, IERR, J, K, PDW, PW, WRKMIN, WRKOPT
      DOUBLE PRECISION  TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMM, DGEQP3, DGEQRF, DLACPY, DLASCL,
     $                  DLASET, DORGQR, DSCAL, MB01UX, MB03TD, MB03ZA,
     $                  MB04DI, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode and check input parameters.
C
      LALL = LSAME( WHICH, 'A' )
      IF ( LALL ) THEN
         LEXT = LSAME( METH, 'L' )
      ELSE
         LEXT = .FALSE.
      END IF
      LUS  = LSAME( STAB,   'S' ) .OR. LSAME( STAB,   'B' )
      LUU  = LSAME( STAB,   'U' ) .OR. LSAME( STAB,   'B' )
      LBAL = LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'S' ) .OR.
     $       LSAME( BALANC, 'B' )
      LBEF = .FALSE.
      IF ( LBAL )
     $   LBEF = LSAME( ORTBAL, 'B' )
C
      WRKMIN = 1
      WRKOPT = WRKMIN
C
      INFO = 0
C
      IF ( .NOT.LALL .AND. .NOT.LSAME( WHICH, 'S' ) ) THEN
         INFO = -1
      ELSE IF ( LALL .AND. ( .NOT.LEXT .AND.
     $                       .NOT.LSAME( METH, 'S' ) ) ) THEN
         INFO = -2
      ELSE IF ( .NOT.LUS .AND. .NOT.LUU ) THEN
         INFO = -3
      ELSE IF ( .NOT.LBAL .AND. .NOT.LSAME( BALANC, 'N' ) ) THEN
         INFO = -4
      ELSE IF ( LBAL .AND. ( .NOT.LBEF .AND.
     $                       .NOT.LSAME( ORTBAL, 'A' ) ) ) THEN
         INFO = -5
      ELSE
         IF ( LALL ) THEN
            M = N
         ELSE
C
C           Set M to the dimension of the specified invariant subspace.
C
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               IF ( PAIR ) THEN
                  PAIR = .FALSE.
               ELSE
                  IF ( K.LT.N ) THEN
                     IF ( S(K+1,K).EQ.ZERO ) THEN
                        IF ( SELECT(K) )
     $                     M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF ( SELECT(K) .OR. SELECT(K+1) )
     $                     M = M + 2
                     END IF
                  ELSE
                     IF ( SELECT(N) )
     $                  M = M + 1
                  END IF
               END IF
   10       CONTINUE
         END IF
C
C        Compute workspace requirements.
C
         IF ( .NOT.LEXT ) THEN
            WRKOPT = MAX( WRKOPT, 4*M*M + MAX( 8*M, 4*N ) )
         ELSE
            IF ( LUS.AND.LUU ) THEN
               WRKOPT = MAX( WRKOPT, 8*N + 1 )
            ELSE
               WRKOPT = MAX( WRKOPT, 2*N*N + 2*N, 8*N )
            END IF
         END IF
C
         IF ( N.LT.0 ) THEN
            INFO = -7
         ELSE IF ( MM.LT.M .OR. ( LEXT .AND. MM.LT.2*N ) ) THEN
            INFO = -8
         ELSE IF ( LBAL .AND. ( ILO.LT.1 .OR. ILO.GT.N+1 ) ) THEN
            INFO = -9
         ELSE IF ( LDS.LT.MAX( 1, N ) ) THEN
            INFO = -12
         ELSE IF ( LDT.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF ( LDG.LT.1 .OR. ( LEXT .AND. LDG.LT.N ) ) THEN
            INFO = -16
         ELSE IF ( LDU1.LT.MAX( 1, N ) ) THEN
            INFO = -18
         ELSE IF ( LDU2.LT.MAX( 1, N ) ) THEN
            INFO = -20
         ELSE IF ( LDV1.LT.MAX( 1, N ) ) THEN
            INFO = -22
         ELSE IF ( LDV2.LT.MAX( 1, N ) ) THEN
            INFO = -24
         ELSE IF ( LDUS.LT.1 .OR. ( LUS .AND. LDUS.LT.2*N ) ) THEN
            INFO = -29
         ELSE IF ( LDUU.LT.1 .OR. ( LUU .AND. LDUU.LT.2*N ) ) THEN
            INFO = -31
         ELSE IF ( LDWORK.LT.WRKMIN ) THEN
            INFO = -35
            DWORK(1) = DBLE( WRKMIN )
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03ZD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( M, N ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
      WRKOPT = WRKMIN
C
      IF ( .NOT.LEXT ) THEN
C
C        Workspace requirements: 4*M*M + MAX( 8*M, 4*N ).
C
         PW  = 1
         PDW = PW + 4*M*M
         CALL MB03ZA( 'No Update', 'Update', 'Update', 'Init', WHICH,
     $                SELECT, N, S, LDS, T, LDT, G, LDG, U1, LDU1, U2,
     $                LDU2, V1, LDV1, V2, LDV2, DWORK(PW), 2*M, WR, WI,
     $                M, DWORK(PDW), LDWORK-PDW+1, IERR )
         IF ( IERR.NE.0 )
     $      GO TO 250
C
         PDW = PW + 2*M*M
         CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, M, ONE,
     $                DWORK(PW), 2*M, V1, LDV1, DWORK(PDW),
     $                LDWORK-PDW+1, IERR )
         WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
         IF ( LUS )
     $      CALL DLACPY( 'All', N, M, V1, LDV1, US, LDUS )
         IF ( LUU )
     $      CALL DLACPY( 'All', N, M, V1, LDV1, UU, LDUU )
C
         CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, M, ONE,
     $                DWORK(PW+M), 2*M, U1, LDU1, DWORK(PDW),
     $                LDWORK-PDW+1, IERR )
C
         IF ( LUS ) THEN
            DO 20 J = 1, M
               CALL DAXPY( N, -ONE, U1(1,J), 1, US(1,J), 1 )
   20       CONTINUE
         END IF
         IF ( LUU ) THEN
            DO 30 J = 1, M
               CALL DAXPY( N, ONE, U1(1,J), 1, UU(1,J), 1 )
   30       CONTINUE
         END IF
C
         CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, M, -ONE,
     $                DWORK(PW), 2*M, V2, LDV2, DWORK(PDW),
     $                LDWORK-PDW+1, IERR )
C
         IF ( LUS )
     $      CALL DLACPY( 'All', N, M, V2, LDV2, US(N+1,1), LDUS )
         IF ( LUU )
     $      CALL DLACPY( 'All', N, M, V2, LDV2, UU(N+1,1), LDUU )
C
         CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, M, ONE,
     $                DWORK(PW+M), 2*M, U2, LDU2, DWORK(PDW),
     $                LDWORK-PDW+1, IERR )
C
         IF ( LUS ) THEN
            DO 40 J = 1, M
               CALL DAXPY( N, ONE, U2(1,J), 1, US(N+1,J), 1 )
   40       CONTINUE
         END IF
         IF ( LUU ) THEN
            DO 50 J = 1, M
               CALL DAXPY( N, -ONE, U2(1,J), 1, UU(N+1,J), 1 )
   50       CONTINUE
         END IF
C
C        Orthonormalize obtained bases and apply inverse balancing
C        transformation.
C
         IF ( LBAL .AND. LBEF ) THEN
            IF ( LUS )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, M, US,
     $                      LDUS, US(N+1,1), LDUS, IERR )
            IF ( LUU )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, M, UU,
     $                      LDUU, UU(N+1,1), LDUU, IERR )
         END IF
C
         IF ( LUS ) THEN
            CALL DGEQRF( 2*N, M, US, LDUS, DWORK(1), DWORK(M+1),
     $                   LDWORK-M, IERR )
               WRKOPT = MAX( WRKOPT, INT( DWORK(M+1) ) + M )
            CALL DORGQR( 2*N, M, M, US, LDUS, DWORK(1), DWORK(M+1),
     $                   LDWORK-M, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(M+1) ) + M )
         END IF
         IF ( LUU ) THEN
            CALL DGEQRF( 2*N, M, UU, LDUU, DWORK(1), DWORK(M+1),
     $                   LDWORK-M, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(M+1) ) + M )
            CALL DORGQR( 2*N, M, M, UU, LDUU, DWORK(1), DWORK(M+1),
     $                   LDWORK-M, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(M+1) ) + M )
         END IF
C
         IF ( LBAL .AND. .NOT.LBEF ) THEN
            IF ( LUS )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, M, US,
     $                      LDUS, US(N+1,1), LDUS, IERR )
            IF ( LUU )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, M, UU,
     $                      LDUU, UU(N+1,1), LDUU, IERR )
         END IF
C
      ELSE
C
         DO 60 I = 1, 2*N
            LWORK(I) = .TRUE.
   60    CONTINUE
C
         IF ( LUS .AND.( .NOT.LUU ) ) THEN
C
C           Workspace requirements: MAX( 2*N*N + 2*N, 8*N )
C
            CALL MB03ZA( 'Update', 'Update', 'Update', 'Init', WHICH,
     $                   SELECT, N, S, LDS, T, LDT, G, LDG, U1, LDU1,
     $                   U2, LDU2, V1, LDV1, V2, LDV2, US, LDUS, WR,
     $                   WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 )
     $         GO TO 250
C
            CALL MB01UX( 'Left', 'Lower', 'Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
C
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
C
            DO 70 J = 1, N
               CALL DAXPY( J, ONE, G(J,1), LDG, G(1,J), 1 )
   70       CONTINUE
            PDW = 2*N*N+1
C
C           DW <- -[V1;V2]*W11
C
            CALL DLACPY( 'All', N, N, V1, LDV1, DWORK, 2*N )
            CALL DLACPY( 'All', N, N, V2, LDV2, DWORK(N+1), 2*N )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', 2*N, N, -ONE,
     $                   US, LDUS, DWORK, 2*N, DWORK(PDW), LDWORK-PDW+1,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
C           DW2 <- DW2 - U2*W21
C
            CALL DLACPY( 'All', N, N, U2, LDU2, US, LDUS )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, N, ONE,
     $                   US(N+1,1), LDUS, US, LDUS, DWORK(PDW),
     $                   LDWORK-PDW+1, IERR )
            DO 80 J = 1, N
               CALL DAXPY( N, ONE, US(1,J), 1, DWORK(N+2*(J-1)*N+1), 1 )
   80       CONTINUE
C
C           US11 <- -U1*W21 - DW1
C
            CALL DLACPY( 'All', N, N, U1, LDU1, US, LDUS )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, N, -ONE,
     $                   US(N+1,1), LDUS, US, LDUS, DWORK(PDW),
     $                   LDWORK-PDW+1, IERR )
            DO 90 J = 1, N
               CALL DAXPY( N, -ONE, DWORK(2*(J-1)*N+1), 1, US(1,J), 1 )
   90       CONTINUE
C
C           US21 <- DW2
C
            CALL DLACPY( 'All', N, N, DWORK(N+1), 2*N, US(N+1,1), LDUS )
C
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, V1, LDV1, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, V2, LDV2, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, U1, LDU1, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, U2, LDU2, DWORK, LDWORK,
     $                   IERR )
            CALL DLACPY( 'All', N, N, V1, LDV1, US(1,N+1), LDUS )
            CALL DLACPY( 'All', N, N, V2, LDV2, US(N+1,N+1), LDUS )
            DO 100 J = 1, N
               CALL DAXPY( N, -ONE, U1(1,J), 1, US(1,N+J), 1 )
  100       CONTINUE
            DO 110 J = 1, N
               CALL DAXPY( N, -ONE, U2(1,J), 1, US(N+1,N+J), 1 )
  110       CONTINUE
C
            CALL MB03TD( 'Hamiltonian', 'Update', LWORK, LWORK(N+1), N,
     $                   S, LDS, G, LDG, US(1,N+1), LDUS, US(N+1,N+1),
     $                   LDUS, WR, WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 ) THEN
               INFO = 4
               RETURN
            END IF
            CALL DLASCL( 'General', 0, 0, ONE, -ONE, N, N, US(N+1,N+1),
     $                   LDUS, IERR )
C
         ELSE IF ( ( .NOT.LUS ).AND.LUU ) THEN
C
C           Workspace requirements: MAX( 2*N*N + 2*N, 8*N )
C
            CALL MB03ZA( 'Update', 'Update', 'Update', 'Init', WHICH,
     $                   SELECT, N, S, LDS, T, LDT, G, LDG, U1, LDU1,
     $                   U2, LDU2, V1, LDV1, V2, LDV2, UU, LDUU, WR,
     $                   WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 )
     $         GO TO 250
            CALL MB01UX( 'Left', 'Lower', 'Transpose', N, N, ONE,
     $                   UU(N+1,N+1), LDUU, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   UU(1,N+1), LDUU, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            DO 120 J = 1, N
               CALL DAXPY( J, ONE, G(J,1), LDG, G(1,J), 1 )
  120       CONTINUE
            PDW = 2*N*N+1
C
C           DW <- -[V1;V2]*W11
C
            CALL DLACPY( 'All', N, N, V1, LDV1, DWORK, 2*N )
            CALL DLACPY( 'All', N, N, V2, LDV2, DWORK(N+1), 2*N )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', 2*N, N, -ONE,
     $                   UU, LDUU, DWORK, 2*N, DWORK(PDW), LDWORK-PDW+1,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(PDW) ) + PDW - 1 )
C
C           DW2 <- DW2 - U2*W21
C
            CALL DLACPY( 'All', N, N, U2, LDU2, UU, LDUU )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, N, -ONE,
     $                   UU(N+1,1), LDUU, UU, LDUU, DWORK(PDW),
     $                   LDWORK-PDW+1, IERR )
            DO 130 J = 1, N
               CALL DAXPY( N, ONE, UU(1,J), 1, DWORK(N+2*(J-1)*N+1), 1 )
  130       CONTINUE
C
C           UU11 <- U1*W21 - DW1
C
            CALL DLACPY( 'All', N, N, U1, LDU1, UU, LDUU )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', N, N, ONE,
     $                   UU(N+1,1), LDUU, UU, LDUU, DWORK(PDW),
     $                   LDWORK-PDW+1, IERR )
            DO 140 J = 1, N
               CALL DAXPY( N, -ONE, DWORK(2*(J-1)*N+1), 1, UU(1,J), 1 )
  140       CONTINUE
C
C           UU21 <- DW2
C
            CALL DLACPY( 'All', N, N, DWORK(N+1), 2*N, UU(N+1,1), LDUU )
C
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   UU(1,N+1), LDUU, V1, LDV1, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   UU(1,N+1), LDUU, V2, LDV2, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   UU(N+1,N+1), LDUU, U1, LDU1, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   UU(N+1,N+1), LDUU, U2, LDU2, DWORK, LDWORK,
     $                   IERR )
            CALL DLACPY( 'All', N, N, V1, LDV1, UU(1,N+1), LDUU )
            CALL DLACPY( 'All', N, N, V2, LDV2, UU(N+1,N+1), LDUU )
            DO 150 J = 1, N
               CALL DAXPY( N, ONE, U1(1,J), 1, UU(1,N+J), 1 )
  150       CONTINUE
            DO 160 J = 1, N
               CALL DAXPY( N, ONE, U2(1,J), 1, UU(N+1,N+J), 1 )
  160       CONTINUE
C
            CALL MB03TD( 'Hamiltonian', 'Update', LWORK, LWORK(N+1), N,
     $                   S, LDS, G, LDG, UU(1,N+1), LDUU, UU(N+1,N+1),
     $                   LDUU, WR, WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 ) THEN
               INFO = 4
               RETURN
            END IF
            CALL DLASCL( 'General', 0, 0, ONE, -ONE, N, N, UU(N+1,N+1),
     $                   LDUU, IERR )
         ELSE
C
C           Workspace requirements: 8*N
C
            CALL MB03ZA( 'Update', 'Update', 'Update', 'Init', WHICH,
     $                   SELECT, N, S, LDS, T, LDT, G, LDG, U1, LDU1,
     $                   U2, LDU2, V1, LDV1, V2, LDV2, US, LDUS, WR,
     $                   WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 )
     $         GO TO 250
            CALL MB01UX( 'Left', 'Lower', 'Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, G, LDG, DWORK, LDWORK,
     $                   IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            DO 170 J = 1, N
               CALL DAXPY( J, ONE, G(J,1), LDG, G(1,J), 1 )
  170       CONTINUE
C
C           UU = [ V1 -V2; U1 -U2 ]*diag(W11,W21)
C
            CALL DLACPY( 'All', N, N, V1, LDV1, UU, LDUU )
            CALL DLACPY( 'All', N, N, V2, LDV2, UU(N+1,1), LDUU )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', 2*N, N, ONE,
     $                   US, LDUS, UU, LDUU, DWORK, LDWORK, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            CALL DLACPY( 'All', N, N, U1, LDU1, UU(1,N+1), LDUU )
            CALL DLACPY( 'All', N, N, U2, LDU2, UU(N+1,N+1), LDUU )
            CALL MB01UX( 'Right', 'Upper', 'No Transpose', 2*N, N, ONE,
     $                   US(N+1,1), LDUS, UU(1,N+1), LDUU, DWORK,
     $                   LDWORK, IERR )
            CALL DLASCL( 'General', 0, 0, ONE, -ONE, N, 2*N, UU(N+1,1),
     $                   LDUU, IERR )
C
            CALL DLACPY( 'All', 2*N, N, UU, LDUU, US, LDUS )
            DO 180 J = 1, N
               CALL DAXPY( 2*N, -ONE, UU(1,N+J), 1, US(1,J), 1 )
  180       CONTINUE
            DO 190 J = 1, N
               CALL DAXPY( 2*N, ONE, UU(1,N+J), 1, UU(1,J), 1 )
  190       CONTINUE
C
C           V1 <- V1*W12-U1*W22
C           U1 <- V1*W12+U1*W22
C           V2 <- V2*W12-U2*W22
C           U2 <- V2*W12+U2*W22
C
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, V1, LDV1, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(1,N+1), LDUS, V2, LDV2, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, U1, LDU1, DWORK, LDWORK,
     $                   IERR )
            CALL MB01UX( 'Right', 'Lower', 'No Transpose', N, N, ONE,
     $                   US(N+1,N+1), LDUS, U2, LDU2, DWORK, LDWORK,
     $                   IERR )
            DO 210 J = 1, N
               DO 200 I = 1, N
                  TEMP = V1(I,J)
                  V1(I,J) = TEMP - U1(I,J)
                  U1(I,J) = TEMP + U1(I,J)
  200          CONTINUE
  210       CONTINUE
            DO 230 J = 1, N
               DO 220 I = 1, N
                  TEMP = V2(I,J)
                  V2(I,J) = TEMP - U2(I,J)
                  U2(I,J) = TEMP + U2(I,J)
  220          CONTINUE
  230       CONTINUE
C
            CALL DLASET( 'All', 2*N, N, ZERO, ONE, US(1,N+1), LDUS )
            CALL MB03TD( 'Hamiltonian', 'Update', LWORK, LWORK(N+1), N,
     $                   S, LDS, G, LDG, US(1,N+1), LDUS, US(N+1,N+1),
     $                   LDUS, WR, WI, M, DWORK, LDWORK, IERR )
            IF ( IERR.NE.0 ) THEN
               INFO = 4
               RETURN
            END IF
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $                  U1, LDU1, US(1,N+1), LDUS, ZERO, UU(1,N+1),
     $                  LDUU )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  U2, LDU2, US(N+1,N+1), LDUS, ONE, UU(1,N+1),
     $                  LDUU )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  U1, LDU1, US(N+1,N+1), LDUS, ZERO, UU(N+1,N+1),
     $                  LDUU )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  U2, LDU2, US(1,N+1), LDUS, ONE, UU(N+1,N+1),
     $                  LDUU )
            CALL DLACPY( 'All', N, N, US(1,N+1), LDUS, U1, LDU1 )
            CALL DLACPY( 'All', N, N, US(N+1,N+1), LDUS, U2, LDU2 )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, ONE,
     $                  V1, LDV1, U1, LDU1, ZERO, US(1,N+1), LDUS )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  V2, LDV2, U2, LDU2, ONE, US(1,N+1), LDUS )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  V1, LDV1, U2, LDU2, ZERO, US(N+1,N+1), LDUS )
            CALL DGEMM( 'No Transpose', 'No Transpose', N, N, N, -ONE,
     $                  V2, LDV2, U1, LDU1, ONE, US(N+1,N+1), LDUS )
         END IF
C
C        Orthonormalize obtained bases and apply inverse balancing
C        transformation.
C
         IF ( LBAL .AND. LBEF ) THEN
            IF ( LUS )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, N, US,
     $                      LDUS, US(N+1,1), LDUS, IERR )
            IF ( LUU )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, N, UU,
     $                      LDUU, UU(N+1,1), LDUU, IERR )
         END IF
C
C        Workspace requirements: 8*N+1
C
         DO 240  J = 1, 2*N
            IWORK(J) = 0
  240    CONTINUE
         IF ( LUS ) THEN
            CALL DGEQP3( 2*N, 2*N, US, LDUS, IWORK, DWORK, DWORK(2*N+1),
     $                   LDWORK-2*N, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(2*N+1) ) + 2*N )
            CALL DORGQR( 2*N, 2*N, N, US, LDUS, DWORK, DWORK(2*N+1),
     $                   LDWORK-2*N, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(2*N+1) ) + 2*N )
         END IF
         IF ( LUU ) THEN
            CALL DGEQP3( 2*N, 2*N, UU, LDUU, IWORK, DWORK, DWORK(2*N+1),
     $                   LDWORK-2*N, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(2*N+1) ) + 2*N )
            CALL DORGQR( 2*N, 2*N, N, UU, LDUU, DWORK, DWORK(2*N+1),
     $                   LDWORK-2*N, IERR )
            WRKOPT = MAX( WRKOPT, INT( DWORK(2*N+1) ) + 2*N )
         END IF
C
         IF ( LBAL .AND. .NOT.LBEF ) THEN
            IF ( LUS )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, N, US,
     $                      LDUS, US(N+1,1), LDUS, IERR )
            IF ( LUU )
     $         CALL MB04DI( BALANC, 'Positive', N, ILO, SCALE, N, UU,
     $                      LDUU, UU(N+1,1), LDUU, IERR )
         END IF
      END IF
C
      CALL DSCAL( M, -ONE, WR, 1 )
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
  250 CONTINUE
      IF ( IERR.EQ.1 ) THEN
         INFO = 2
      ELSE IF ( IERR.EQ.2 .OR. IERR.EQ.4 ) THEN
         INFO = 1
      ELSE IF ( IERR.EQ.3 ) THEN
         INFO = 3
      END IF
      RETURN
C *** Last line of MB03ZD ***
      END
