      SUBROUTINE MB03ZA( COMPC, COMPU, COMPV, COMPW, WHICH, SELECT, N,
     $                   A, LDA, B, LDB, C, LDC, U1, LDU1, U2, LDU2, V1,
     $                   LDV1, V2, LDV2, W, LDW, WR, WI, M, DWORK,
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
C     1. To compute, for a given matrix pair (A,B) in periodic Schur
C        form, orthogonal matrices Ur and Vr so that
C
C            T           [ A11  A12 ]     T           [ B11  B12 ]
C          Vr * A * Ur = [          ],  Ur * B * Vr = [          ], (1)
C                        [  0   A22 ]                 [  0   B22 ]
C
C        is in periodic Schur form, and the eigenvalues of A11*B11
C        form a selected cluster of eigenvalues.
C
C     2. To compute an orthogonal matrix W so that
C
C                   T  [  0  -A11 ]       [  R11   R12 ]
C                  W * [          ] * W = [            ],           (2)
C                      [ B11   0  ]       [   0    R22 ]
C
C        where the eigenvalues of R11 and -R22 coincide and have
C        positive real part.
C
C     Optionally, the matrix C is overwritten by Ur'*C*Vr.
C
C     All eigenvalues of A11*B11 must either be complex or real and
C     negative.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     COMPC   CHARACTER*1
C             = 'U':  update the matrix C;
C             = 'N':  do not update C.
C
C     COMPU   CHARACTER*1
C             = 'U':  update the matrices U1 and U2;
C             = 'N':  do not update U1 and U2.
C             See the description of U1 and U2.
C
C     COMPV   CHARACTER*1
C             = 'U':  update the matrices V1 and V2;
C             = 'N':  do not update V1 and V2.
C             See the description of V1 and V2.
C
C     COMPW   CHARACTER*1
C             Indicates whether or not the user wishes to accumulate
C             the matrix W as follows:
C             = 'N':  the matrix W is not required;
C             = 'I':  W is initialized to the unit matrix and the
C                     orthogonal transformation matrix W is returned;
C             = 'V':  W must contain an orthogonal matrix Q on entry,
C                     and the product Q*W is returned.
C
C     WHICH   CHARACTER*1
C             = 'A':  select all eigenvalues, this effectively means
C                     that Ur and Vr are identity matrices and A11 = A,
C                     B11 = B;
C             = 'S':  select a cluster of eigenvalues specified by
C                     SELECT.
C
C     SELECT  LOGICAL array, dimension (N)
C             If WHICH = 'S', then SELECT specifies the eigenvalues of
C             A*B in the selected cluster. To select a real eigenvalue
C             w(j), SELECT(j) must be set to .TRUE.. To select a complex
C             conjugate pair of eigenvalues w(j) and w(j+1),
C             corresponding to a 2-by-2 diagonal block in A, both
C             SELECT(j) and SELECT(j+1) must be set to .TRUE.; a complex
C             conjugate pair of eigenvalues must be either both included
C             in the cluster or both excluded.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi-triangular matrix A of the matrix
C             pair (A,B) in periodic Schur form.
C             On exit, the leading M-by-M part of this array contains
C             the matrix R22 in (2).
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular matrix B of the matrix pair
C             (A,B) in periodic Schur form.
C             On exit, the leading N-by-N part of this array is
C             overwritten.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, if COMPC = 'U', the leading N-by-N part of this
C             array must contain a general matrix C.
C             On exit, if COMPC = 'U', the leading N-by-N part of this
C             array contains the updated matrix Ur'*C*Vr.
C             If COMPC = 'N' or WHICH = 'A', this array is not
C             referenced.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= 1.
C             LDC >= N,  if COMPC = 'U' and WHICH = 'S'.
C
C     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On entry, if COMPU = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array must contain U1, the (1,1)
C             block of an orthogonal symplectic matrix
C             U = [ U1, U2; -U2, U1 ].
C             On exit, if COMPU = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array contains U1*Ur.
C             If COMPU = 'N' or WHICH = 'A', this array is not
C             referenced.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.  LDU1 >= 1.
C             LDU1 >= N,  if COMPU = 'U' and WHICH = 'S'.
C
C     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On entry, if COMPU = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array must contain U2, the (1,2)
C             block of an orthogonal symplectic matrix
C             U = [ U1, U2; -U2, U1 ].
C             On exit, if COMPU = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array contains U2*Ur.
C             If COMPU = 'N' or WHICH = 'A', this array is not
C             referenced.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.  LDU2 >= 1.
C             LDU2 >= N,  if COMPU = 'U' and WHICH = 'S'.
C
C     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N)
C             On entry, if COMPV = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array must contain V1, the (1,1)
C             block of an orthogonal symplectic matrix
C             V = [ V1, V2; -V2, V1 ].
C             On exit, if COMPV = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array contains V1*Vr.
C             If COMPV = 'N' or WHICH = 'A', this array is not
C             referenced.
C
C     LDV1    INTEGER
C             The leading dimension of the array V1.  LDV1 >= 1.
C             LDV1 >= N,  if COMPV = 'U' and WHICH = 'S'.
C
C     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,N)
C             On entry, if COMPV = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array must contain V2, the (1,2)
C             block of an orthogonal symplectic matrix
C             V = [ V1, V2; -V2, V1 ].
C             On exit, if COMPV = 'U' and WHICH = 'S', the leading
C             N-by-N part of this array contains V2*Vr.
C             If COMPV = 'N' or WHICH = 'A', this array is not
C             referenced.
C
C     LDV2    INTEGER
C             The leading dimension of the array V2.  LDV2 >= 1.
C             LDV2 >= N,  if COMPV = 'U' and WHICH = 'S'.
C
C     W       (input/output) DOUBLE PRECISION array, dimension (LDW,2*M)
C             On entry, if COMPW = 'V', then the leading 2*M-by-2*M part
C             of this array must contain a matrix W.
C             If COMPW = 'I', then W need not be set on entry, W is set
C             to the identity matrix.
C             On exit, if COMPW = 'I' or 'V' the leading 2*M-by-2*M part
C             of this array is post-multiplied by the transformation
C             matrix that produced (2).
C             If COMPW = 'N', this array is not referenced.
C
C     LDW     INTEGER
C             The leading dimension of the array W.  LDW >= 1.
C             LDW >= 2*M,  if COMPW = 'I' or COMPW = 'V'.
C
C     WR      (output) DOUBLE PRECISION array, dimension (M)
C     WI      (output) DOUBLE PRECISION array, dimension (M)
C             The real and imaginary parts, respectively, of the
C             eigenvalues of R22. The eigenvalues are stored in the same
C             order as on the diagonal of R22, with
C             WR(i) = R22(i,i) and, if R22(i:i+1,i:i+1) is a 2-by-2
C             diagonal block, WI(i) > 0 and WI(i+1) = -WI(i).
C             In exact arithmetic, these eigenvalue are the positive
C             square roots of the selected eigenvalues of the product
C             A*B. However, if an eigenvalue is sufficiently
C             ill-conditioned, then its value may differ significantly.
C
C     M       (output) INTEGER
C             The number of selected eigenvalues.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if  INFO = -28,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, 4*N, 8*M ).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  reordering of the product A*B in Step 1 failed
C                   because some eigenvalues are too close to separate;
C             = 2:  reordering of some submatrix in Step 2 failed
C                   because some eigenvalues are too close to separate;
C             = 3:  the QR algorithm failed to compute the Schur form
C                   of some submatrix in Step 2;
C             = 4:  the condition that all eigenvalues of A11*B11 must
C                   either be complex or real and negative is
C                   numerically violated.
C
C     METHOD
C
C     Step 1 is performed using a reordering technique analogous to the
C     LAPACK routine DTGSEN for reordering matrix pencils [1,2]. Step 2
C     is an implementation of Algorithm 2 in [3]. It requires O(M*N*N)
C     floating point operations.
C
C     REFERENCES
C
C     [1] Kagstrom, B.
C         A direct method for reordering eigenvalues in the generalized
C         real Schur form of a regular matrix pair (A,B), in M.S. Moonen
C         et al (eds), Linear Algebra for Large Scale and Real-Time
C         Applications, Kluwer Academic Publ., 1993, pp. 195-218.
C
C     [2] Kagstrom, B. and Poromaa P.:
C         Computing eigenspaces with specified eigenvalues of a regular
C         matrix pair (A, B) and condition estimation: Theory,
C         algorithms and software, Numer. Algorithms, 1996, vol. 12,
C         pp. 369-407.
C
C     [3] Benner, P., Mehrmann, V., and Xu, H.
C         A new method for computing the stable invariant subspace of a
C         real Hamiltonian matrix,  J. Comput. Appl. Math., 86,
C         pp. 17-43, 1997.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLABMX).
C
C     KEYWORDS
C
C     Hamiltonian matrix, invariant subspace.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDQZ
      PARAMETER         ( LDQZ = 4 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COMPC, COMPU, COMPV, COMPW, WHICH
      INTEGER           INFO, LDA, LDB, LDC, LDU1, LDU2, LDV1, LDV2,
     $                  LDW, LDWORK, M, N
C     .. Array Arguments ..
      LOGICAL           SELECT(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*),
     $                  U1(LDU1,*), U2(LDU2,*), V1(LDV1,*), V2(LDV2,*),
     $                  W(LDW,*), WI(*), WR(*)
C     .. Local Scalars ..
      LOGICAL           CMPALL, INITW, PAIR, SWAP, WANTC, WANTU, WANTV,
     $                  WANTW
      INTEGER           HERE, I, IERR, IFST, ILST, K, KS, L, LEN, MM,
     $                  NB, NBF, NBL, NBNEXT, POS, PW, PWC, PWCK, PWD,
     $                  PWDL, WRKMIN
      DOUBLE PRECISION  TEMP
C     .. Local Arrays ..
      LOGICAL           LDUM(1), SELNEW(4)
      DOUBLE PRECISION  DW12(12), Q(LDQZ,LDQZ), T(LDQZ,LDQZ), WINEW(4),
     $                  WRNEW(4), Z(LDQZ,LDQZ)
      INTEGER           IDUM(1)
C     .. External Functions ..
      LOGICAL           LFDUM, LSAME
      EXTERNAL          LFDUM, LSAME
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEES, DGEMM, DLACPY, DLASET, DSCAL,
     $                  DTRSEN, MB03WA, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Decode and check input parameters
C
      WANTC  = LSAME( COMPC, 'U' )
      WANTU  = LSAME( COMPU, 'U' )
      WANTV  = LSAME( COMPV, 'U' )
      INITW  = LSAME( COMPW, 'I' )
      WANTW  = INITW .OR. LSAME( COMPW, 'V' )
      CMPALL = LSAME( WHICH, 'A' )
      WRKMIN = MAX( 1, 4*N )
C
      INFO = 0
      IF ( .NOT.WANTC .AND. .NOT.LSAME( COMPC, 'N' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.WANTU .AND. .NOT.LSAME( COMPU, 'N' ) ) THEN
         INFO = -2
      ELSE IF ( .NOT.WANTV .AND. .NOT.LSAME( COMPV, 'N' ) ) THEN
         INFO = -3
      ELSE IF ( .NOT.WANTW .AND. .NOT.LSAME( COMPW, 'N' ) ) THEN
         INFO = -4
      ELSE IF ( .NOT.CMPALL .AND. .NOT.LSAME( WHICH, 'S' ) ) THEN
         INFO = -5
      ELSE
         IF ( CMPALL ) THEN
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
                     IF ( A(K+1,K).EQ.ZERO ) THEN
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
         WRKMIN = MAX( WRKMIN, 8*M )
C
         IF ( N.LT.0 ) THEN
            INFO = -7
         ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -9
         ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -11
         ELSE IF ( LDC.LT.1  .OR. ( WANTC .AND. .NOT.CMPALL .AND.
     $                              LDC.LT.N ) ) THEN
            INFO = -13
         ELSE IF ( LDU1.LT.1 .OR. ( WANTU .AND. .NOT.CMPALL .AND.
     $                              LDU1.LT.N ) ) THEN
            INFO = -15
         ELSE IF ( LDU2.LT.1 .OR. ( WANTU .AND. .NOT.CMPALL .AND.
     $                              LDU2.LT.N ) ) THEN
            INFO = -17
         ELSE IF ( LDV1.LT.1 .OR. ( WANTV .AND. .NOT.CMPALL .AND.
     $                              LDV1.LT.N ) ) THEN
            INFO = -19
         ELSE IF ( LDV2.LT.1 .OR. ( WANTV .AND. .NOT.CMPALL .AND.
     $                              LDV2.LT.N ) ) THEN
            INFO = -21
         ELSE IF ( LDW.LT.1 .OR. ( WANTW .AND. LDW.LT.2*M ) ) THEN
            INFO = -23
         ELSE IF ( LDWORK.LT.WRKMIN ) THEN
            INFO = -28
            DWORK(1) = DBLE( WRKMIN )
         END IF
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03ZA', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Jump immediately to Step 2, if all eigenvalues are requested.
C
      IF ( CMPALL )
     $   GO TO 50
C
C     Step 1: Collect the selected blocks at the top-left corner of A*B.
C
      KS = 0
      PAIR = .FALSE.
      DO 40 K = 1, N
         IF ( PAIR ) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT(K)
            IF  ( K.LT.N ) THEN
               IF ( A(K+1,K).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT(K+1)
               END IF
            END IF
C
            IF ( PAIR ) THEN
               NBF = 2
            ELSE
               NBF = 1
            END IF
C
            IF ( SWAP ) THEN
               KS = KS + 1
               IFST = K
C
C              Swap the K-th block to position KS.
C
               ILST = KS
               NBL = 1
               IF ( ILST.GT.1 ) THEN
                  IF ( A(ILST,ILST-1).NE.ZERO ) THEN
                     ILST = ILST - 1
                     NBL = 2
                  END IF
               END IF
C
               IF ( ILST.EQ.IFST )
     $            GO TO 30
C
               HERE = IFST
   20          CONTINUE
C
C              Swap block with next one above.
C
               IF ( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
C
C                 Current block either 1-by-1 or 2-by-2.
C
                  NBNEXT = 1
                  IF ( HERE.GE.3 ) THEN
                     IF ( A(HERE-1,HERE-2).NE.ZERO )
     $                  NBNEXT = 2
                  END IF
                  POS = HERE - NBNEXT
                  NB  = NBNEXT + NBF
                  CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                  CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                  CALL MB03WA( .TRUE., .TRUE., NBNEXT, NBF, A(POS,POS),
     $                         LDA, B(POS,POS), LDB, Q, LDQZ, Z, LDQZ,
     $                         IERR )
C
                  IF ( IERR.NE.0 ) THEN
                     DWORK(1) = DBLE( WRKMIN )
                     INFO = 1
                     RETURN
                  END IF
C
C                 Update rest of A.
C
                  IF ( POS.GT.1 ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', POS-1,
     $                           NB, NB, ONE, A(1,POS), LDA, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', POS-1, NB, DWORK, N, A(1,POS),
     $                            LDA )
                  END IF
                  IF ( POS+NB.LE.N ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                           A(POS,POS+NB), LDA, ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                            A(POS,POS+NB), LDA )
                  END IF
C
C                 Update rest of B.
C
                  IF ( POS.GT.1 ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', POS-1,
     $                           NB, NB, ONE, B(1,POS), LDB, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', POS-1, NB, DWORK, N, B(1,POS),
     $                            LDB )
                  END IF
                  IF ( POS+NB.LE.N ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                           B(POS,POS+NB), LDB, ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                            B(POS,POS+NB), LDB )
                  END IF
C
C                 Update C.
C
                  IF ( WANTC ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, C(1,POS), LDC, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, C(1,POS),
     $                            LDC )
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N, NB, ONE, Z, LDQZ, C(POS,1), LDC,
     $                           ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N, DWORK, NB, C(POS,1),
     $                            LDC )
                  END IF
C
C                 Update U.
C
                  IF ( WANTU ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, U1(1,POS), LDU1, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, U1(1,POS),
     $                            LDU1 )
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, U2(1,POS), LDU2, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, U2(1,POS),
     $                            LDU2 )
                  END IF
C
C                 Update V.
C
                  IF ( WANTV ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, V1(1,POS), LDV1, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, V1(1,POS),
     $                            LDV1 )
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, V2(1,POS), LDV2, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, V2(1,POS),
     $                            LDV2 )
                  END IF
C
                  HERE = HERE - NBNEXT
C
C                 Test if 2-by-2 block breaks into two 1-by-1 blocks.
C
                  IF ( NBF.EQ.2 ) THEN
                     IF ( A(HERE+1,HERE).EQ.ZERO )
     $                  NBF = 3
                  END IF
C
               ELSE
C
C                 Current block consists of two 1 by 1 blocks each of
C                 which must be swapped individually.
C
                  NBNEXT = 1
                  IF ( HERE.GE.3 ) THEN
                     IF ( A(HERE-1,HERE-2).NE.ZERO )
     $                  NBNEXT = 2
                  END IF
                  POS = HERE - NBNEXT
                  NB  = NBNEXT + 1
                  CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                  CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                  CALL MB03WA( .TRUE., .TRUE., NBNEXT, 1, A(POS,POS),
     $                         LDA, B(POS,POS), LDB, Q, LDQZ, Z, LDQZ,
     $                         IERR )
C
                  IF ( IERR.NE.0 ) THEN
                     DWORK(1) = DBLE( WRKMIN )
                     INFO = 1
                     RETURN
                  END IF
C
C                 Update rest of A.
C
                  IF ( POS.GT.1 ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', POS-1,
     $                           NB, NB, ONE, A(1,POS), LDA, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', POS-1, NB, DWORK, N, A(1,POS),
     $                            LDA )
                  END IF
                  IF ( POS+NB.LE.N ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                           A(POS,POS+NB), LDA, ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                            A(POS,POS+NB), LDA )
                  END IF
C
C                 Update rest of B.
C
                  IF ( POS.GT.1 ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', POS-1,
     $                           NB, NB, ONE, B(1,POS), LDB, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', POS-1, NB, DWORK, N, B(1,POS),
     $                            LDB )
                  END IF
                  IF ( POS+NB.LE.N ) THEN
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                           B(POS,POS+NB), LDB, ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                            B(POS,POS+NB), LDB )
                  END IF
C
C                 Update C.
C
                  IF ( WANTC ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, C(1,POS), LDC, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, C(1,POS),
     $                            LDC )
                     CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                           N, NB, ONE, Z, LDQZ, C(POS,1), LDC,
     $                           ZERO, DWORK, NB )
                     CALL DLACPY( 'All', NB, N, DWORK, NB, C(POS,1),
     $                            LDC )
                  END IF
C
C                 Update U.
C
                  IF ( WANTU ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, U1(1,POS), LDU1, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, U1(1,POS),
     $                            LDU1 )
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, U2(1,POS), LDU2, Z, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, U2(1,POS),
     $                            LDU2 )
                  END IF
C
C                 Update V.
C
                  IF ( WANTV ) THEN
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, V1(1,POS), LDV1, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, V1(1,POS),
     $                            LDV1 )
                     CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                           NB, NB, ONE, V2(1,POS), LDV2, Q, LDQZ,
     $                           ZERO, DWORK, N )
                     CALL DLACPY( 'All', N, NB, DWORK, N, V2(1,POS),
     $                            LDV2 )
                  END IF
C
                  IF ( NBNEXT.EQ.1 ) THEN
C
C                    Swap two 1-by-1 blocks.
C
                     POS = HERE
                     NB  = NBNEXT + 1
                     CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                     CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                     CALL MB03WA( .TRUE., .TRUE., NBNEXT, 1, A(POS,POS),
     $                            LDA, B(POS,POS), LDB, Q, LDQZ, Z,
     $                            LDQZ, IERR )
C
                     IF ( IERR.NE.0 ) THEN
                        DWORK(1) = DBLE( WRKMIN )
                        INFO = 1
                        RETURN
                     END IF
C
C                    Update rest of A.
C
                     IF ( POS.GT.1 ) THEN
                        CALL DGEMM( 'No Transpose', 'No Transpose',
     $                              POS-1, NB, NB, ONE, A(1,POS), LDA,
     $                              Z, LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                               A(1,POS), LDA )
                     END IF
                     IF ( POS+NB.LE.N ) THEN
                        CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                              N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                              A(POS,POS+NB), LDA, ZERO, DWORK,
     $                              NB )
                        CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                               A(POS,POS+NB), LDA )
                     END IF
C
C                    Update rest of B.
C
                     IF ( POS.GT.1 ) THEN
                        CALL DGEMM( 'No Transpose', 'No Transpose',
     $                              POS-1, NB, NB, ONE, B(1,POS), LDB,
     $                              Q, LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                               B(1,POS), LDB )
                     END IF
                     IF ( POS+NB.LE.N ) THEN
                        CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                              N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                              B(POS,POS+NB), LDB, ZERO, DWORK,
     $                              NB )
                        CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK, NB,
     $                               B(POS,POS+NB), LDB )
                     END IF
C
C                    Update C.
C
                     IF ( WANTC ) THEN
                        CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                              NB, NB, ONE, C(1,POS), LDC, Q, LDQZ,
     $                              ZERO, DWORK, N )
                        CALL DLACPY( 'All', N, NB, DWORK, N, C(1,POS),
     $                               LDC )
                        CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                              N, NB, ONE, Z, LDQZ, C(POS,1), LDC,
     $                              ZERO, DWORK, NB )
                        CALL DLACPY( 'All', NB, N, DWORK, NB, C(POS,1),
     $                               LDC )
                     END IF
C
C                    Update U.
C
                     IF ( WANTU ) THEN
                        CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                              NB, NB, ONE, U1(1,POS), LDU1, Z,
     $                              LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', N, NB, DWORK, N, U1(1,POS),
     $                               LDU1 )
                        CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                              NB, NB, ONE, U2(1,POS), LDU2, Z,
     $                              LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', N, NB, DWORK, N, U2(1,POS),
     $                               LDU2 )
                     END IF
C
C                    Update V.
C
                     IF ( WANTV ) THEN
                        CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                              NB, NB, ONE, V1(1,POS), LDV1, Q,
     $                              LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', N, NB, DWORK, N, V1(1,POS),
     $                               LDV1 )
                        CALL DGEMM( 'No Transpose', 'No Transpose', N,
     $                              NB, NB, ONE, V2(1,POS), LDV2, Q,
     $                              LDQZ, ZERO, DWORK, N )
                        CALL DLACPY( 'All', N, NB, DWORK, N, V2(1,POS),
     $                               LDV2 )
                     END IF
C
                     HERE = HERE - 1
                  ELSE
C
C                    Recompute NBNEXT in case 2-by-2 split.
C
                     IF ( A(HERE,HERE-1).EQ.ZERO )
     $                  NBNEXT = 1
C
                     IF ( NBNEXT.EQ.2 ) THEN
C
C                       2-by-2 block did not split.
C
                        POS = HERE - 1
                        NB  = 3
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                        CALL MB03WA( .TRUE., .TRUE., 2, 1, A(POS,POS),
     $                               LDA, B(POS,POS), LDB, Q, LDQZ, Z,
     $                               LDQZ, IERR )
C
                        IF ( IERR.NE.0 ) THEN
                           DWORK(1) = DBLE( WRKMIN )
                           INFO = 1
                           RETURN
                        END IF
C
C                       Update rest of A.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, A(1,POS),
     $                                 LDA, Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  A(1,POS), LDA )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                                 A(POS,POS+NB), LDA, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, A(POS,POS+NB), LDA )
                        END IF
C
C                       Update rest of B.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, B(1,POS),
     $                                 LDB, Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  B(1,POS), LDB )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                                 B(POS,POS+NB), LDB, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, B(POS,POS+NB), LDB )
                        END IF
C
C                       Update C.
C
                        IF ( WANTC ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, C(1,POS), LDC, Q,
     $                                 LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  C(1,POS), LDC )
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N, NB, ONE, Z, LDQZ, C(POS,1),
     $                                 LDC, ZERO, DWORK, NB )
                           CALL DLACPY( 'All', NB, N, DWORK, NB,
     $                                  C(POS,1), LDC )
                        END IF
C
C                       Update U.
C
                        IF ( WANTU ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U1(1,POS), LDU1,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                 U1(1,POS), LDU1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U2(1,POS), LDU2,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  U2(1,POS), LDU2 )
                        END IF
C
C                       Update V.
C
                        IF ( WANTV ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V1(1,POS), LDV1,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V1(1,POS), LDV1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V2(1,POS), LDV2,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V2(1,POS), LDV2 )
                        END IF
C
                        HERE = HERE - 2
                     ELSE
C
C                       2-by-2 block did split.
C
                        POS = HERE
                        NB = 2
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                        CALL MB03WA( .TRUE., .TRUE., 2, 1, A(POS,POS),
     $                               LDA, B(POS,POS), LDB, Q, LDQZ, Z,
     $                               LDQZ, IERR )
C
                        IF ( IERR.NE.0 ) THEN
                           DWORK(1) = DBLE( WRKMIN )
                           INFO = 1
                           RETURN
                        END IF
C
C                       Update rest of A.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, A(1,POS),
     $                                 LDA, Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  A(1,POS), LDA )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                                 A(POS,POS+NB), LDA, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, A(POS,POS+NB), LDA )
                        END IF
C
C                       Update rest of B.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, B(1,POS),
     $                                 LDB, Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  B(1,POS), LDB )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                                 B(POS,POS+NB), LDB, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, B(POS,POS+NB), LDB )
                        END IF
C
C                       Update C.
C
                        IF ( WANTC ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, C(1,POS), LDC, Q,
     $                                 LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  C(1,POS), LDC )
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N, NB, ONE, Z, LDQZ, C(POS,1),
     $                                 LDC, ZERO, DWORK, NB )
                           CALL DLACPY( 'All', NB, N, DWORK, NB,
     $                                  C(POS,1), LDC )
                        END IF
C
C                       Update U.
C
                        IF ( WANTU ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U1(1,POS), LDU1,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                 U1(1,POS), LDU1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U2(1,POS), LDU2,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  U2(1,POS), LDU2 )
                        END IF
C
C                       Update V.
C
                        IF ( WANTV ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V1(1,POS), LDV1,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V1(1,POS), LDV1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V2(1,POS), LDV2,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V2(1,POS), LDV2 )
                        END IF
C
                        POS = HERE - 1
                        NB = 2
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Q, LDQZ )
                        CALL DLASET( 'All', NB, NB, ZERO, ONE, Z, LDQZ )
C
                        CALL MB03WA( .TRUE., .TRUE., 2, 1, A(POS,POS),
     $                               LDA, B(POS,POS), LDB, Q, LDQZ, Z,
     $                               LDQZ, IERR )
C
                        IF ( IERR.NE.0 ) THEN
                           DWORK(1) = DBLE( WRKMIN )
                           INFO = 1
                           RETURN
                        END IF
C
C                       Update rest of A.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, A(1,POS),
     $                                 LDA, Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  A(1,POS), LDA )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Q, LDQZ,
     $                                 A(POS,POS+NB), LDA, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, A(POS,POS+NB), LDA )
                        END IF
C
C                       Update rest of B.
C
                        IF ( POS.GT.1 ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 POS-1, NB, NB, ONE, B(1,POS),
     $                                 LDB, Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', POS-1, NB, DWORK, N,
     $                                  B(1,POS), LDB )
                        END IF
                        IF ( POS+NB.LE.N ) THEN
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N-POS-NB+1, NB, ONE, Z, LDQZ,
     $                                 B(POS,POS+NB), LDB, ZERO, DWORK,
     $                                 NB )
                           CALL DLACPY( 'All', NB, N-POS-NB+1, DWORK,
     $                                  NB, B(POS,POS+NB), LDB )
                        END IF
C
C                       Update C.
C
                        IF ( WANTC ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, C(1,POS), LDC, Q,
     $                                 LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  C(1,POS), LDC )
                           CALL DGEMM( 'Transpose', 'No Transpose', NB,
     $                                 N, NB, ONE, Z, LDQZ, C(POS,1),
     $                                 LDC, ZERO, DWORK, NB )
                           CALL DLACPY( 'All', NB, N, DWORK, NB,
     $                                  C(POS,1), LDC )
                        END IF
C
C                       Update U.
C
                        IF ( WANTU ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U1(1,POS), LDU1,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                 U1(1,POS), LDU1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, U2(1,POS), LDU2,
     $                                 Z, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  U2(1,POS), LDU2 )
                        END IF
C
C                       Update V.
C
                        IF ( WANTV ) THEN
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V1(1,POS), LDV1,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V1(1,POS), LDV1 )
                           CALL DGEMM( 'No Transpose', 'No Transpose',
     $                                 N, NB, NB, ONE, V2(1,POS), LDV2,
     $                                 Q, LDQZ, ZERO, DWORK, N )
                           CALL DLACPY( 'All', N, NB, DWORK, N,
     $                                  V2(1,POS), LDV2 )
                        END IF
C
                        HERE = HERE - 2
                     END IF
                  END IF
               END IF
C
               IF ( HERE.GT.ILST )
     $            GO TO 20
C
   30          CONTINUE
               IF ( PAIR )
     $            KS = KS + 1
            END IF
         END IF
   40 CONTINUE
C
   50 CONTINUE
C
C     Step 2: Compute an ordered Schur decomposition of
C             [ 0, -A11; B11, 0 ].
C
      IF ( INITW )
     $   CALL DLASET( 'All', 2*M, 2*M, ZERO, ONE, W, LDW )
      PWC = 1
      PWD = PWC + 2*M
      PW  = PWD + 2*M
      PAIR = .FALSE.
      NB = 1
C
      DO 80 K = 1, M
         IF ( PAIR ) THEN
            PAIR = .FALSE.
            NB = 1
         ELSE
            IF ( K.LT.N ) THEN
               IF ( A(K+1,K).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  NB = 2
               END IF
            END IF
            PWCK = PWC + 2*( K - 1 )
            PWDL = PWD + 2*( K - 1 )
            CALL DLASET( 'All', NB, M-K+1, ZERO, ZERO, DWORK(PWCK), 2 )
            CALL DLACPY( 'All', NB, M-K+1, A(K,K), LDA, DWORK(PWDL), 2 )
            CALL DLASET( 'All', NB, M-K+1, ZERO, ZERO, A(K,K), LDA )
C
            L = K
C
C           WHILE L >= 1 DO
C
   60       CONTINUE
C
               IF ( K.EQ.L ) THEN
C
C                 Annihilate B(k,k).
C
                  NBL = NB
                  CALL DLASET( 'All', NB+NBL, NB+NBL, ZERO, ZERO, T,
     $                         LDQZ )
                  CALL DLACPY( 'Upper', NBL, NBL, B(L,L), LDB,
     $                         T(NB+1,1), LDQZ )
                  IF ( NB.EQ.1 ) THEN
                     DWORK(PWDL) = -DWORK(PWDL)
                  ELSE
                     CALL DSCAL( 2*NB, -ONE, DWORK(PWDL), 1 )
                  END IF
                  CALL DLACPY( 'All', NB, NB, DWORK(PWDL), 2, T(1,NB+1),
     $                         LDQZ )
               ELSE
C
C                 Annihilate B(l,k).
C
                  CALL DLASET( 'All', NBL+NB, NBL+NB, ZERO, ZERO, T,
     $                         LDQZ )
                  CALL DLACPY( 'All', NBL, NBL, A(L,L), LDA, T, LDQZ )
                  CALL DLACPY( 'All', NBL, NB, B(L,K), LDB, T(1,NBL+1),
     $                         LDQZ )
                  CALL DLACPY( 'All', NB, NB, DWORK(PWCK), 2,
     $                         T(NBL+1,NBL+1), LDQZ )
                  PWDL = PWD + 2*( L - 1 )
               END IF
C
               CALL DGEES( 'V', 'Not Sorted', LFDUM, NB+NBL, T, LDQZ,
     $                     MM, WRNEW, WINEW, Q, LDQZ, DW12, 12, LDUM,
     $                     IERR )
               IF ( IERR.NE.0 ) THEN
                  DWORK(1) = DBLE( WRKMIN )
                  INFO = 3
                  RETURN
               END IF
C
C              Reorder Schur form.
C
               MM = 0
               DO 70  I = 1, NB+NBL
                  IF ( WRNEW(I).GT.0 ) THEN
                     MM = MM + 1
                     SELNEW(I) = .TRUE.
                  ELSE
                     SELNEW(I) = .FALSE.
                  END IF
   70          CONTINUE
               IF ( MM.LT.NB ) THEN
                  DWORK(1) = DBLE( WRKMIN )
                  INFO = 4
                  RETURN
               END IF
               CALL DTRSEN( 'None', 'V', SELNEW, NB+NBL, T, LDQZ, Q,
     $                      LDQZ, WRNEW, WINEW, MM, TEMP, TEMP, DW12,
     $                      4, IDUM, 1, IERR )
               IF ( IERR.NE.0 ) THEN
                  DWORK(1) = DBLE( WRKMIN )
                  INFO = 2
                  RETURN
               END IF
C
C              Permute Q if necessary.
C
               IF ( K.NE.L ) THEN
                  CALL DLACPY( 'All', NBL, NB+NBL, Q, LDQZ, Z(NB+1,1),
     $                         LDQZ )
                  CALL DLACPY( 'All', NB, NB+NBL, Q(NBL+1,1), LDQZ,
     $                         Z, LDQZ )
                  CALL DLACPY( 'All', NB+NBL, NB+NBL, Z, LDQZ, Q, LDQZ )
               END IF
C
C              Update "diagonal" blocks.
C
               CALL DLACPY( 'All', NB, NB, T, LDQZ, DWORK(PWCK), 2 )
               CALL DLACPY( 'All', NB, NBL, T(1,NB+1), LDQZ,
     $                      DWORK(PWDL), 2 )
               IF ( NB.EQ.1 ) THEN
                  CALL DSCAL( NBL, -ONE, DWORK(PWDL), 2 )
               ELSE
                  CALL DSCAL( 2*NBL, -ONE, DWORK(PWDL), 1 )
               END IF
               CALL DLACPY( 'All', NBL, NBL, T(NB+1,NB+1), LDQZ,
     $                      A(L,L), LDA )
C
C              Update block columns of A and B.
C
               LEN = L - 1
               IF ( LEN.GT.0 ) THEN
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NB, ONE, B(1,K), LDB, Q, LDQZ, ZERO,
     $                        DWORK(PW), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NB, ONE, B(1,K), LDB, Q(1,NB+1), LDQZ,
     $                        ZERO, DWORK(PW+2*M), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NBL, ONE, A(1,L), LDA, Q(NB+1,1), LDQZ,
     $                        ONE, DWORK(PW), M )
                  CALL DLACPY( 'All', LEN, NB, DWORK(PW), M, B(1,K),
     $                         LDB )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NBL, ONE, A(1,L), LDA, Q(NB+1,NB+1),
     $                        LDQZ, ONE, DWORK(PW+2*M), M )
                  CALL DLACPY( 'All', LEN, NBL, DWORK(PW+2*M), M,
     $                         A(1,L), LDA )
               END IF
C
C              Update block column of A.
C
               LEN = M - L - NBL + 1
               IF ( LEN.GT.0 ) THEN
                  CALL DGEMM( 'Transpose', 'No Transpose', NB, LEN, NB,
     $                        ONE, Q, LDQZ, DWORK(PWDL+2*NBL), 2, ZERO,
     $                        DWORK(PW), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NBL, LEN, NB,
     $                        -ONE, Q(1,NB+1), LDQZ, DWORK(PWDL+2*NBL),
     $                        2, ZERO, DWORK(PW+2*M), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NB, LEN, NBL,
     $                        -ONE, Q(NB+1,1), LDQZ, A(L,L+NBL), LDA,
     $                        ONE, DWORK(PW), 2 )
                  CALL DLACPY( 'All', NB, LEN, DWORK(PW), 2,
     $                         DWORK(PWDL+2*NBL), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NBL, LEN,
     $                        NBL, ONE, Q(NB+1,NB+1), LDQZ, A(L,L+NBL),
     $                        LDA, ONE, DWORK(PW+2*M), 2 )
                  CALL DLACPY( 'All', NBL, LEN, DWORK(PW+2*M), 2,
     $                         A(L,L+NBL), LDA )
               END IF
C
C              Update block row of B.
C
               LEN = M - K - NB + 1
               IF ( LEN.GT.0 ) THEN
                  CALL DGEMM( 'Transpose', 'No Transpose', NB, LEN, NB,
     $                        ONE, Q, LDQZ, DWORK(PWCK+2*NB), 2, ZERO,
     $                        DWORK(PW), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NBL, LEN, NB,
     $                        ONE, Q(1,NB+1), LDQZ, DWORK(PWCK+2*NB), 2,
     $                        ZERO, DWORK(PW+2*M), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NB, LEN, NBL,
     $                        ONE, Q(NB+1,1), LDQZ, B(L,K+NB), LDB, ONE,
     $                        DWORK(PW), 2 )
                  CALL DLACPY( 'All', NB, LEN, DWORK(PW), 2,
     $                         DWORK(PWCK+2*NB), 2 )
                  CALL DGEMM( 'Transpose', 'No Transpose', NBL, LEN,
     $                        NBL, ONE, Q(NB+1,NB+1), LDQZ, B(L,K+NB),
     $                        LDB, ONE, DWORK(PW+2*M), 2 )
                  CALL DLACPY( 'All', NBL, LEN, DWORK(PW+2*M), 2,
     $                         B(L,K+NB), LDB )
               END IF
C
C              Update W.
C
               IF ( WANTW ) THEN
                  IF ( INITW ) THEN
                     POS = L
                     LEN = K + NB - L
                  ELSE
                     POS = 1
                     LEN = M
                  END IF
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NB, ONE, W(POS,K), LDW, Q, LDQZ, ZERO,
     $                        DWORK(PW), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NB, ONE, W(POS,K), LDW, Q(1,NB+1), LDQZ,
     $                        ZERO, DWORK(PW+2*M), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NBL, ONE, W(POS,M+L), LDW, Q(NB+1,1),
     $                        LDQZ, ONE, DWORK(PW), M )
                  CALL DLACPY( 'All', LEN, NB, DWORK(PW), M, W(POS,K),
     $                         LDW )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NBL, ONE, W(POS,M+L), LDW, Q(NB+1,NB+1),
     $                        LDQZ, ONE, DWORK(PW+2*M), M )
                  CALL DLACPY( 'All', LEN, NBL, DWORK(PW+2*M), M,
     $                         W(POS,M+L), LDW )
C
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NB, ONE, W(M+POS,K), LDW, Q, LDQZ, ZERO,
     $                        DWORK(PW), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NB, ONE, W(M+POS,K), LDW, Q(1,NB+1), LDQZ,
     $                        ZERO, DWORK(PW+2*M), M )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NB,
     $                        NBL, ONE, W(M+POS,M+L), LDW, Q(NB+1,1),
     $                        LDQZ, ONE, DWORK(PW), M )
                  CALL DLACPY( 'All', LEN, NB, DWORK(PW), M, W(M+POS,K),
     $                         LDW )
                  CALL DGEMM( 'No Transpose', 'No Transpose', LEN, NBL,
     $                        NBL, ONE, W(M+POS,M+L), LDW, Q(NB+1,NB+1),
     $                        LDQZ, ONE, DWORK(PW+2*M), M )
                  CALL DLACPY( 'All', LEN, NBL, DWORK(PW+2*M), M,
     $                         W(M+POS,M+L), LDW )
               END IF
C
               L = L - 1
               NBL = 1
               IF ( L.GT.1 ) THEN
                  IF ( A(L,L-1).NE.ZERO ) THEN
                     NBL = 2
                     L = L - 1
                  END IF
               END IF
C
C           END WHILE L >= 1 DO
C
            IF ( L.GE.1 )
     $         GO TO 60
C
C           Copy recomputed eigenvalues.
C
            CALL DCOPY( NB, WRNEW, 1, WR(K), 1 )
            CALL DCOPY( NB, WINEW, 1, WI(K), 1 )
         END IF
   80 CONTINUE
      DWORK(1) = DBLE( WRKMIN )
      RETURN
C *** Last line of MB03ZA ***
      END
C
      LOGICAL FUNCTION LFDUM( X, Y )
C
C     Void logical function for DGEES.
C
      DOUBLE PRECISION X, Y
      LFDUM = .FALSE.
      RETURN
C *** Last line of LFDUM ***
      END
