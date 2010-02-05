      SUBROUTINE MB03TD( TYP, COMPU, SELECT, LOWER, N, A, LDA, G, LDG,
     $                   U1, LDU1, U2, LDU2, WR, WI, M, DWORK, LDWORK,
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
C     To reorder a matrix X in skew-Hamiltonian Schur form:
C
C                   [  A   G  ]          T
C             X  =  [       T ],   G = -G,
C                   [  0   A  ]
C
C     or in Hamiltonian Schur form:
C
C                   [  A   G  ]          T
C             X  =  [       T ],   G =  G,
C                   [  0  -A  ]
C
C     where A is in upper quasi-triangular form, so that a selected
C     cluster of eigenvalues appears in the leading diagonal blocks
C     of the matrix A (in X) and the leading columns of [ U1; -U2 ] form
C     an orthonormal basis for the corresponding right invariant
C     subspace.
C
C     If X is skew-Hamiltonian, then each eigenvalue appears twice; one
C     copy corresponds to the j-th diagonal element and the other to the
C     (n+j)-th diagonal element of X. The logical array LOWER controls
C     which copy is to be reordered to the leading part of A.
C
C     If X is Hamiltonian then the eigenvalues appear in pairs
C     (lambda,-lambda); lambda corresponds to the j-th diagonal
C     element and -lambda to the (n+j)-th diagonal element of X.
C     The logical array LOWER controls whether lambda or -lambda is to
C     be reordered to the leading part of A.
C
C     The matrix A must be in Schur canonical form (as returned by the
C     LAPACK routine DHSEQR), that is, block upper triangular with
C     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
C     its diagonal elements equal and its off-diagonal elements of
C     opposite sign.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYP     CHARACTER*1
C             Specifies the type of the input matrix X:
C             = 'S': X is skew-Hamiltonian;
C             = 'H': X is Hamiltonian.
C
C     COMPU   CHARACTER*1
C             = 'U': update the matrices U1 and U2 containing the
C                    Schur vectors;
C             = 'N': do not update U1 and U2.
C
C     SELECT  (input/output) LOGICAL array, dimension (N)
C             SELECT specifies the eigenvalues in the selected cluster.
C             To select a real eigenvalue w(j), SELECT(j) must be set
C             to .TRUE.. To select a complex conjugate pair of
C             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2
C             diagonal block, both SELECT(j) and SELECT(j+1) must be set
C             to .TRUE.; a complex conjugate pair of eigenvalues must be
C             either both included in the cluster or both excluded.
C
C     LOWER   (input/output) LOGICAL array, dimension (N)
C             LOWER controls which copy of a selected eigenvalue is
C             included in the cluster. If SELECT(j) is set to .TRUE.
C             for a real eigenvalue w(j); then LOWER(j) must be set to
C             .TRUE. if the eigenvalue corresponding to the (n+j)-th
C             diagonal element of X is to be reordered to the leading
C             part; and LOWER(j) must be set to .FALSE. if the
C             eigenvalue corresponding to the j-th diagonal element of
C             X is to be reordered to the leading part. Similarly, for
C             a complex conjugate pair of eigenvalues w(j) and w(j+1),
C             both LOWER(j) and LOWER(j+1) must be set to .TRUE. if the
C             eigenvalues corresponding to the (n+j:n+j+1,n+j:n+j+1)
C             diagonal block of X are to be reordered to the leading
C             part.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi-triangular matrix A in Schur
C             canonical form.
C             On exit, the leading N-by-N part of this array contains
C             the reordered matrix A, again in Schur canonical form,
C             with the selected eigenvalues in the diagonal blocks.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N)
C             On entry, if TYP = 'S', the leading N-by-N part of this
C             array must contain the strictly upper triangular part of
C             the skew-symmetric matrix G. The rest of this array is not
C             referenced.
C             On entry, if TYP = 'H', the leading N-by-N part of this
C             array must contain the upper triangular part of the
C             symmetric matrix G. The rest of this array is not
C             referenced.
C             On exit, if TYP = 'S', the leading N-by-N part of this
C             array contains the strictly upper triangular part of the
C             skew-symmetric matrix G, updated by the orthogonal
C             symplectic transformation which reorders X.
C             On exit, if TYP = 'H', the leading N-by-N part of this
C             array contains the upper triangular part of the symmetric
C             matrix G, updated by the orthogonal symplectic
C             transformation which reorders X.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= MAX(1,N).
C
C     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On entry, if COMPU = 'U', the leading N-by-N part of this
C             array must contain U1, the (1,1) block of an orthogonal
C             symplectic matrix U = [ U1, U2; -U2, U1 ].
C             On exit, if COMPU = 'U', the leading N-by-N part of this
C             array contains the (1,1) block of the matrix U,
C             postmultiplied by the orthogonal symplectic transformation
C             which reorders X. The leading M columns of U form an
C             orthonormal basis for the specified invariant subspace.
C             If COMPU = 'N', this array is not referenced.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.
C             LDU1 >= MAX(1,N),  if COMPU = 'U';
C             LDU1 >= 1,         otherwise.
C
C     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On entry, if COMPU = 'U', the leading N-by-N part of this
C             array must contain U2, the (1,2) block of an orthogonal
C             symplectic matrix U = [ U1, U2; -U2, U1 ].
C             On exit, if COMPU = 'U', the leading N-by-N part of this
C             array contains the (1,2) block of the matrix U,
C             postmultiplied by the orthogonal symplectic transformation
C             which reorders X.
C             If COMPU = 'N', this array is not referenced.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.
C             LDU2 >= MAX(1,N),  if COMPU = 'U';
C             LDU2 >= 1,         otherwise.
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             The real and imaginary parts, respectively, of the
C             reordered eigenvalues of A. The eigenvalues are stored
C             in the same order as on the diagonal of A, with
C             WR(i) = A(i,i) and, if A(i:i+1,i:i+1) is a 2-by-2 diagonal
C             block, WI(i) > 0 and WI(i+1) = -WI(i). Note that if an
C             eigenvalue is sufficiently ill-conditioned, then its value
C             may differ significantly from its value before reordering.
C
C     M       (output) INTEGER
C             The dimension of the specified invariant subspace.
C             0 <= M <= N.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -18,  DWORK(1)  returns the minimum
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
C                   value.
C             = 1:  reordering of X failed because some eigenvalue pairs
C                   are too close to separate (the problem is very
C                   ill-conditioned); X may have been partially
C                   reordered, and WR and WI contain the eigenvalues in
C                   the same order as in X.
C
C     REFERENCES
C
C     [1] Bai, Z. and Demmel, J.W.
C         On Swapping Diagonal Blocks in Real Schur Form.
C         Linear Algebra Appl., 186, pp. 73-95, 1993.
C
C     [2] Benner, P., Kressner, D., and Mehrmann, V.
C         Skew-Hamiltonian and Hamiltonian Eigenvalue Problems: Theory,
C         Algorithms and Applications. Techn. Report, TU Berlin, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAORD).
C
C     KEYWORDS
C
C     Hamiltonian matrix, skew-Hamiltonian matrix, invariant subspace.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         COMPU, TYP
      INTEGER           INFO, LDA, LDG, LDU1, LDU2, LDWORK, M, N
C     .. Array Arguments ..
      LOGICAL           LOWER(*), SELECT(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), G(LDG,*), U1(LDU1,*),
     $                  U2(LDU2,*), WI(*), WR(*)
C     .. Local Scalars ..
      LOGICAL           FLOW, ISHAM, PAIR, SWAP, WANTU
      INTEGER           HERE, IERR, IFST, ILST, K, KS, NBF, NBL, NBNEXT,
     $                  WRKMIN
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          MB03TS, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Decode and check input parameters.
C
      ISHAM = LSAME( TYP, 'H' )
      WANTU = LSAME( COMPU, 'U' )
      WRKMIN = MAX( 1, N )
      INFO = 0
      IF ( .NOT.ISHAM .AND. .NOT.LSAME( TYP, 'S' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.WANTU .AND. .NOT.LSAME( COMPU, 'N' ) ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -5
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF ( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF ( LDU1.LT.1 .OR. ( WANTU .AND. LDU1.LT.N ) ) THEN
         INFO = -11
      ELSE IF ( LDU2.LT.1 .OR. ( WANTU .AND. LDU2.LT.N ) ) THEN
         INFO = -13
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         INFO = -18
         DWORK(1) = DBLE( WRKMIN )
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03TD', -INFO )
         RETURN
      END IF
C
C     Set M to the dimension of the specified invariant subspace.
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
     $               M = M + 1
               ELSE
                  PAIR = .TRUE.
                  IF ( SELECT(K) .OR. SELECT(K+1) )
     $               M = M + 2
               END IF
            ELSE
               IF ( SELECT(N) )
     $            M = M + 1
            END IF
         END IF
   10 CONTINUE
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Collect the selected blocks at the top-left corner of X.
C
      KS = 0
      PAIR = .FALSE.
      DO 60 K = 1, N
         IF ( PAIR ) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT(K)
            FLOW = LOWER(K)
            IF  ( K.LT.N ) THEN
               IF ( A(K+1,K).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP.OR.SELECT(K+1)
                  FLOW = FLOW.OR.LOWER(K+1)
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
               IF ( FLOW ) THEN
C
C                 Step 1: Swap the K-th block to position N.
C
                  IFST = K
                  ILST = N
                  NBL = 1
                  IF ( ILST.GT.1 ) THEN
                     IF ( A(ILST,ILST-1).NE.ZERO ) THEN
                        ILST = ILST - 1
                        NBL = 2
                     END IF
                  END IF
C
C                 Update ILST.
C
                  IF ( NBF.EQ.2 .AND. NBL.EQ.1 )
     $               ILST = ILST - 1
                  IF ( NBF.EQ.1 .AND. NBL.EQ.2 )
     $               ILST = ILST + 1
C
                  IF ( ILST.EQ.IFST )
     $               GO TO 30
C
                  HERE = IFST
C
   20             CONTINUE
C
C                 Swap block with next one below.
C
                  IF ( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
C
C                    Current block is either 1-by-1 or 2-by-2.
C
                     NBNEXT = 1
                     IF ( HERE+NBF+1.LE.N ) THEN
                        IF ( A(HERE+NBF+1,HERE+NBF).NE.ZERO )
     $                     NBNEXT = 2
                     END IF
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, HERE, NBF, NBNEXT,
     $                            DWORK, IERR )
                     IF ( IERR.NE.0 ) THEN
                        INFO = 1
                        GO TO 70
                     END IF
                     HERE = HERE + NBNEXT
C
C                    Test if 2-by-2 block breaks into two 1-by-1 blocks.
C
                     IF ( NBF.EQ.2 ) THEN
                        IF ( A(HERE+1,HERE).EQ.ZERO )
     $                     NBF = 3
                     END IF
C
                  ELSE
C
C                    Current block consists of two 1-by-1 blocks each of
C                    which must be swapped individually.
C
                     NBNEXT = 1
                     IF ( HERE+3.LE.N ) THEN
                        IF ( A(HERE+3,HERE+2).NE.ZERO )
     $                     NBNEXT = 2
                     END IF
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, HERE+1, 1, NBNEXT,
     $                            DWORK, IERR )
                     IF ( IERR.NE.0 ) THEN
                        INFO = 1
                        GO TO 70
                     END IF
                     IF ( NBNEXT.EQ.1 ) THEN
C
C                       Swap two 1-by-1 blocks, no problems possible.
C
                        CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                               U1, LDU1, U2, LDU2, HERE, 1,
     $                               NBNEXT, DWORK, IERR )
                        HERE = HERE + 1
                     ELSE
C
C                       Recompute NBNEXT in case 2 by 2 split.
C
                        IF ( A(HERE+2,HERE+1).EQ.ZERO )
     $                     NBNEXT = 1
                        IF ( NBNEXT.EQ.2 ) THEN
C
C                          2-by-2 block did not split
C
                           CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                                  U1, LDU1, U2, LDU2, HERE, 1,
     $                                  NBNEXT, DWORK, IERR )
                           IF ( IERR.NE.0 ) THEN
                              INFO = 1
                              GO TO 70
                           END IF
                           HERE = HERE + 2
                        ELSE
C
C                          2-by-2 block did split
C
                           CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                                  U1, LDU1, U2, LDU2, HERE, 1, 1,
     $                                  DWORK, IERR )
                           CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                                  U1, LDU1, U2, LDU2, HERE+1, 1,
     $                                  1, DWORK, IERR )
                           HERE = HERE + 2
                        END IF
                     END IF
                  END IF
                  IF ( HERE.LT.ILST )
     $               GO TO 20
C
   30             CONTINUE
C
C                 Step 2: Apply an orthogonal symplectic transformation
C                         to swap the last blocks in A and -A' (or A').
C
                  IF ( NBF.EQ.1 ) THEN
C
C                    Exchange columns/rows N <-> 2*N. No problems
C                    possible.
C
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                            U1, LDU1, U2, LDU2, N, 1, 1,
     $                            DWORK, IERR )
C
                  ELSE IF ( NBF.EQ.2 ) THEN
C
C                    Swap last block with its equivalent by an
C                    orthogonal symplectic transformation.
C
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                            U1, LDU1, U2, LDU2, N-1, 2, 2,
     $                            DWORK, IERR )
                     IF ( IERR.NE.0 ) THEN
                        INFO = 1
                        GO TO 70
                     END IF
C
C                    Test if 2-by-2 block breaks into two 1-by-1 blocks.
C
                     IF ( A(N-1,N).EQ.ZERO )
     $                  NBF = 3
                  ELSE
C
C                    Block did split. Swap (N-1)-th and N-th elements
C                    consecutively by symplectic generalized
C                    permutations and one rotation.
C
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, N, 1, 1, DWORK, IERR )
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, N-1, 1, 1, DWORK,
     $                            IERR )
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, N, 1, 1, DWORK, IERR )
                  END IF
                  IFST = N
                  IF ( PAIR )
     $               IFST = N-1
               ELSE
                  IFST = K
               END IF
C
C              Step 3: Swap the K-th / N-th block to position KS.
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
     $            GO TO 50
C
               HERE = IFST
   40          CONTINUE
C
C              Swap block with next one above.
C
               IF ( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
C
C                 Current block either 1 by 1 or 2 by 2.
C
                  NBNEXT = 1
                  IF ( HERE.GE.3 ) THEN
                     IF ( A(HERE-1,HERE-2).NE.ZERO )
     $                  NBNEXT = 2
                  END IF
                  CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                         LDU1, U2, LDU2, HERE-NBNEXT, NBNEXT,
     $                         NBF, DWORK, IERR )
                  IF ( IERR.NE.0 ) THEN
                     INFO = 1
                     GO TO 70
                  END IF
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
                  CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                         LDU1, U2, LDU2, HERE-NBNEXT, NBNEXT,
     $                         1, DWORK, IERR )
                  IF ( IERR.NE.0 ) THEN
                     INFO = 1
                     GO TO 70
                  END IF
                  IF ( NBNEXT.EQ.1 ) THEN
C
C                    Swap two 1-by-1 blocks, no problems possible.
C
                     CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1,
     $                            LDU1, U2, LDU2, HERE, NBNEXT, 1,
     $                            DWORK, IERR )

                     HERE = HERE - 1
                  ELSE
C
C                    Recompute NBNEXT in case 2-by-2 split.
C
                     IF ( A(HERE,HERE-1).EQ.ZERO )
     $                  NBNEXT = 1
                     IF ( NBNEXT.EQ.2 ) THEN
C
C                       2-by-2 block did not split
C
                        CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                               U1, LDU1, U2, LDU2, HERE-1, 2, 1,
     $                               DWORK, IERR )
                        IF ( IERR.NE.0 ) THEN
                           INFO = 1
                           GO TO 70
                        END IF
                        HERE = HERE - 2
                     ELSE
C
C                       2-by-2 block did split
C
                        CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                               U1, LDU1, U2, LDU2, HERE, 1, 1,
     $                               DWORK, IERR )
                        CALL MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG,
     $                               U1, LDU1, U2, LDU2, HERE-1, 1, 1,
     $                               DWORK, IERR )
                        HERE = HERE - 2
                     END IF
                  END IF
               END IF
C
               IF ( HERE.GT.ILST )
     $            GO TO 40
C
   50          CONTINUE
               IF ( PAIR )
     $            KS = KS + 1
            END IF
         END IF
   60 CONTINUE
C
   70 CONTINUE
C
C     Store eigenvalues.
C
      DO 80 K = 1, N
         WR(K) = A(K,K)
         WI(K) = ZERO
   80 CONTINUE
      DO 90 K = 1, N - 1
         IF ( A(K+1,K).NE.ZERO ) THEN
            WI(K) = SQRT( ABS( A(K,K+1) ) )*
     $              SQRT( ABS( A(K+1,K) ) )
            WI(K+1) = -WI(K)
         END IF
   90 CONTINUE
C
      DWORK(1) = DBLE( WRKMIN )
C
      RETURN
C *** Last line of MB03TD ***
      END
