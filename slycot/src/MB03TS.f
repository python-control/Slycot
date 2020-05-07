      SUBROUTINE MB03TS( ISHAM, WANTU, N, A, LDA, G, LDG, U1, LDU1, U2,
     $                   LDU2, J1, N1, N2, DWORK, INFO )
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
C     To swap diagonal blocks A11 and A22 of order 1 or 2 in the upper
C     quasi-triangular matrix A contained in a skew-Hamiltonian matrix
C
C                   [  A   G  ]          T
C             X  =  [       T ],   G = -G,
C                   [  0   A  ]
C
C     or in a Hamiltonian matrix
C
C                   [  A   G  ]          T
C             X  =  [       T ],   G =  G.
C                   [  0  -A  ]
C
C     This routine is a modified version of the LAPACK subroutine
C     DLAEX2.
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
C     ISHAM   LOGIGAL
C             Specifies the type of X:
C             = .TRUE.:   X is a Hamiltonian matrix;
C             = .FALSE.:  X is a skew-Hamiltonian matrix.
C
C     WANTU   LOGIGAL
C             = .TRUE.:   update the matrices U1 and U2 containing the
C                         Schur vectors;
C             = .FALSE.:  do not update U1 and U2.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi-triangular matrix A, in Schur
C             canonical form.
C             On exit, the leading N-by-N part of this array contains
C             the reordered matrix A, again in Schur canonical form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular part of the symmetric
C             matrix G, if ISHAM = .TRUE., or the strictly upper
C             triangular part of the skew-symmetric matrix G, otherwise.
C             The rest of this array is not referenced.
C             On exit, the leading N-by-N part of this array contains
C             the upper or strictly upper triangular part of the
C             symmetric or skew-symmetric matrix G, respectively,
C             updated by the orthogonal transformation which reorders A.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= MAX(1,N).
C
C     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N)
C             On entry, if WANTU = .TRUE., the leading N-by-N part of
C             this array must contain the matrix U1.
C             On exit, if WANTU = .TRUE., the leading N-by-N part of
C             this array contains U1, postmultiplied by the orthogonal
C             transformation which reorders A. See the description in
C             the SLICOT subroutine MB03TD for further details.
C             If WANTU = .FALSE., this array is not referenced.
C
C     LDU1    INTEGER
C             The leading dimension of the array U1.
C             LDU1 >= MAX(1,N),  if WANTU = .TRUE.;
C             LDU1 >= 1,         otherwise.
C
C     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N)
C             On entry, if WANTU = .TRUE., the leading N-by-N part of
C             this array must contain the matrix U2.
C             On exit, if WANTU = .TRUE., the leading N-by-N part of
C             this array contains U2, postmultiplied by the orthogonal
C             transformation which reorders A.
C             If WANTU = .FALSE., this array is not referenced.
C
C     LDU2    INTEGER
C             The leading dimension of the array U2.
C             LDU2 >= MAX(1,N),  if WANTU = .TRUE.;
C             LDU2 >= 1,         otherwise.
C
C     J1      (input) INTEGER
C             The index of the first row of the first block A11.
C             If J1+N1 < N, then A11 is swapped with the block starting
C             at (J1+N1+1)-th diagonal element.
C             If J1+N1 > N, then A11 is the last block in A and swapped
C             with -A11', if ISHAM = .TRUE.,
C             or    A11', if ISHAM = .FALSE..
C
C     N1      (input) INTEGER
C             The order of the first block A11. N1 = 0, 1 or 2.
C
C     N2      (input) INTEGER
C             The order of the second block A22. N2 = 0, 1 or 2.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  the transformed matrix A would be too far from Schur
C                   form; the blocks are not swapped and A, G, U1 and
C                   U2 are unchanged.
C
C     REFERENCES
C
C     [1] Bai, Z., and Demmel, J.W.
C        On swapping diagonal blocks in real Schur form.
C        Linear Algebra Appl., 186, pp. 73-95, 1993.
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
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAEX2).
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, THIRTY, FORTY
      PARAMETER          ( ZERO  = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
     $                     TWO   = 2.0D+0, THIRTY = 3.0D+1,
     $                     FORTY = 4.0D+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
C     .. Scalar Arguments ..
      LOGICAL            ISHAM, WANTU
      INTEGER            INFO, J1, LDA, LDG, LDU1, LDU2, N, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), DWORK(*), G(LDG,*), U1(LDU1,*),
     $                   U2(LDU2,*)
C     .. Local Scalars ..
      LOGICAL            LBLK
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   A11, A22, A33, CS, DNORM, EPS, SCALE, SMLNUM,
     $                   SN, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
     $                   WR1, WR2, XNORM
C     .. Local Arrays ..
      DOUBLE PRECISION   D(LDD,4), V(3), V1(3), V2(3), X(LDX,2)
C     .. External Functions ..
      DOUBLE PRECISION   DDOT, DLAMCH, DLANGE
      EXTERNAL           DDOT, DLAMCH, DLANGE
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DLACPY, DLANV2, DLARFG, DLARFX, DLARTG,
     $                   DLASET, DLASY2, DROT, DSCAL, DSWAP, DSYMV,
     $                   DSYR2, MB01MD, MB01ND
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible.
C
      IF ( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
     $   RETURN
      LBLK = ( J1+N1.GT.N )
C
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
C
      IF ( LBLK .AND. N1.EQ.1 ) THEN
C
         IF ( ISHAM ) THEN
            A11 = A(N,N)
            CALL DLARTG( G(N,N), -TWO*A11, CS, SN, TEMP )
            CALL DROT( N-1, A(1,N), 1, G(1,N), 1, CS, SN )
            A(N,N) = -A11
            IF ( WANTU )
     $         CALL DROT( N, U1(1,N), 1, U2(1,N), 1, CS, SN )
         ELSE
            CALL DSWAP( N-1, A(1,N), 1, G(1,N), 1 )
            CALL DSCAL( N-1, -ONE, A(1,N), 1 )
            IF ( WANTU ) THEN
               CALL DSWAP( N, U1(1,N), 1, U2(1,N), 1 )
               CALL DSCAL( N, -ONE, U1(1,N), 1 )
            END IF
         END IF
C
      ELSE IF ( LBLK .AND. N1.EQ.2 ) THEN
C
         IF ( ISHAM ) THEN
C
C           Reorder Hamiltonian matrix:
C
C                      [ A11  G11  ]
C                      [         T ].
C                      [  0  -A11  ]
C
            ND = 4
            CALL DLACPY( 'Full', 2, 2, A(N-1,N-1), LDA, D, LDD )
            CALL DLASET( 'All', 2, 2, ZERO, ZERO, D(3,1), LDD )
            CALL DLACPY( 'Upper', 2, 2, G(N-1,N-1), LDG, D(1,3), LDD )
            D(2,3) = D(1,4)
            D(3,3) = -D(1,1)
            D(4,3) = -D(1,2)
            D(3,4) = -D(2,1)
            D(4,4) = -D(2,2)
            DNORM = DLANGE( 'Max', ND, ND, D, LDD, DWORK )
C
C           Compute machine-dependent threshold for test for accepting
C           swap.
C
            EPS = DLAMCH( 'P' )
            SMLNUM = DLAMCH( 'S' ) / EPS
            THRESH = MAX( FORTY*EPS*DNORM, SMLNUM )
C
C           Solve A11*X + X*A11' = scale*G11 for X.
C
            CALL DLASY2( .FALSE., .FALSE., -1, 2, 2, D, LDD, D(3,3),
     $                   LDD, D(1,3), LDD, SCALE, X, LDX, XNORM, IERR )
C
C           Compute symplectic QR decomposition of
C
C                  (  -X11  -X12 )
C                  (  -X21  -X22 ).
C                  ( scale    0  )
C                  (    0  scale )
C
            TEMP = -X(1,1)
            CALL DLARTG( TEMP, SCALE, V1(1), V2(1), X(1,1) )
            CALL DLARTG( X(1,1), -X(2,1), V1(2), V2(2), TEMP )
            X(1,2) = -X(1,2)
            X(2,2) = -X(2,2)
            X(1,1) = ZERO
            X(2,1) = SCALE
            CALL DROT( 1, X(1,2), 1, X(1,1), 1, V1(1), V2(1) )
            CALL DROT( 1, X(1,2), 1, X(2,2), 1, V1(2), V2(2) )
            CALL DROT( 1, X(1,1), 1, X(2,1), 1, V1(2), V2(2) )
            CALL DLARTG( X(2,2), X(2,1), V1(3), V2(3), TEMP )
C
C           Perform swap provisionally on D.
C
            CALL DROT( 4, D(1,1), LDD, D(3,1), LDD, V1(1), V2(1) )
            CALL DROT( 4, D(1,1), LDD, D(2,1), LDD, V1(2), V2(2) )
            CALL DROT( 4, D(3,1), LDD, D(4,1), LDD, V1(2), V2(2) )
            CALL DROT( 4, D(2,1), LDD, D(4,1), LDD, V1(3), V2(3) )
            CALL DROT( 4, D(1,1), 1, D(1,3), 1, V1(1), V2(1) )
            CALL DROT( 4, D(1,1), 1, D(1,2), 1, V1(2), V2(2) )
            CALL DROT( 4, D(1,3), 1, D(1,4), 1, V1(2), V2(2) )
            CALL DROT( 4, D(1,2), 1, D(1,4), 1, V1(3), V2(3) )
C
C           Test whether to reject swap.
C
            IF ( MAX( ABS( D(3,1) ), ABS( D(3,2) ), ABS( D(4,1) ),
     $                ABS( D(4,2) ) ).GT.THRESH ) GO TO 50
C
            CALL DLACPY( 'All', 2, 2, D(1,1), LDD, A(N-1,N-1), LDA )
            CALL DLACPY( 'Upper', 2, 2, D(1,3), LDD, G(N-1,N-1), LDG )
C
            IF ( N.GT.2 ) THEN
               CALL DROT( N-2, A(1,N-1), 1, G(1,N-1), 1, V1(1), V2(1) )
               CALL DROT( N-2, A(1,N-1), 1, A(1,N), 1, V1(2), V2(2) )
               CALL DROT( N-2, G(1,N-1), 1, G(1,N), 1, V1(2), V2(2) )
               CALL DROT( N-2, A(1,N), 1, G(1,N), 1, V1(3), V2(3) )
            END IF
C
            IF ( WANTU ) THEN
               CALL DROT( N, U1(1,N-1), 1, U2(1,N-1), 1, V1(1), V2(1) )
               CALL DROT( N, U1(1,N-1), 1, U1(1,N), 1, V1(2), V2(2) )
               CALL DROT( N, U2(1,N-1), 1, U2(1,N), 1, V1(2), V2(2) )
               CALL DROT( N, U1(1,N), 1, U2(1,N), 1, V1(3), V2(3) )
            END IF
C
         ELSE
C
            IF ( ABS( A(N-1,N) ).GT.ABS( A(N,N-1) ) ) THEN
               TEMP = G(N-1,N)
               CALL DLARTG( TEMP, A(N-1,N), CS, SN, G(N-1,N) )
               SN = -SN
               CALL DROT(N-2, A(1,N), 1, G(1,N), 1, CS, SN )
C
               A(N-1,N) = -SN*A(N,N-1)
               TEMP = -CS*A(N,N-1)
               A(N,N-1) = G(N-1,N)
               G(N-1,N) = TEMP
               IF ( WANTU )
     $            CALL DROT( N, U1(1,N), 1, U2(1,N), 1, CS, SN )
               CALL DSWAP( N-2, A(1,N-1), 1, G(1,N-1), 1 )
               CALL DSCAL( N-2, -ONE, A(1,N-1), 1 )
               IF ( WANTU ) THEN
                  CALL DSWAP( N, U1(1,N-1), 1, U2(1,N-1), 1 )
                  CALL DSCAL( N, -ONE, U1(1,N-1), 1 )
               END IF
            ELSE
               TEMP = G(N-1,N)
               CALL DLARTG( TEMP, A(N,N-1), CS, SN, G(N-1,N) )
               CALL DROT( N-2, A(1,N-1), 1, G(1,N-1), 1, CS, SN )
               A(N,N-1) = -SN*A(N-1,N)
               A(N-1,N) = CS*A(N-1,N)
               IF ( WANTU )
     $            CALL DROT( N, U1(1,N-1), 1, U2(1,N-1), 1, CS, SN )
               CALL DSWAP( N-1, A(1,N), 1, G(1,N), 1 )
               CALL DSCAL( N-1, -ONE, A(1,N), 1 )
               IF ( WANTU ) THEN
                  CALL DSWAP( N, U1(1,N), 1, U2(1,N), 1 )
                  CALL DSCAL( N, -ONE, U1(1,N), 1 )
               END IF
            END IF
         END IF
C
C        Standardize new 2-by-2 block.
C
         CALL DLANV2( A(N-1,N-1), A(N-1,N), A(N,N-1),
     $                A(N,N), WR1, WI1, WR2, WI2, CS, SN )
         CALL DROT( N-2, A(1,N-1), 1, A(1,N), 1, CS, SN )
         IF ( ISHAM ) THEN
            TEMP = G(N-1,N)
            CALL DROT( N-1, G(1,N-1), 1, G(1,N), 1, CS, SN )
            TAU = CS*TEMP + SN*G(N,N)
            G(N,N) = CS*G(N,N) - SN*TEMP
            G(N-1,N-1) = CS*G(N-1,N-1) + SN*TAU
            CALL DROT( 1, G(N-1,N), LDG, G(N,N), LDG, CS, SN )
         ELSE
            CALL DROT( N-2, G(1,N-1), 1, G(1,N), 1, CS, SN )
         END IF
         IF ( WANTU ) THEN
            CALL DROT( N, U1(1,N-1), 1, U1(1,N), 1, CS, SN )
            CALL DROT( N, U2(1,N-1), 1, U2(1,N), 1, CS, SN )
         END IF
C
      ELSE IF ( N1.EQ.1 .AND. N2.EQ.1 ) THEN
C
C        Swap two 1-by-1 blocks.
C
         A11 = A(J1,J1)
         A22 = A(J2,J2)
C
C        Determine the transformation to perform the interchange.
C
         CALL DLARTG( A(J1,J2), A22-A11, CS, SN, TEMP )
C
C        Apply transformation to the matrix A.
C
         IF ( J3.LE.N )
     $      CALL DROT( N-J1-1, A(J1,J3), LDA, A(J2,J3), LDA, CS, SN )
         CALL DROT( J1-1, A(1,J1), 1, A(1,J2), 1, CS, SN )
C
         A(J1,J1) = A22
         A(J2,J2) = A11
C
C        Apply transformation to the matrix G.
C
         IF ( ISHAM ) THEN
            TEMP = G(J1,J2)
            CALL DROT( J1, G(1,J1), 1, G(1,J2), 1, CS, SN )
            TAU = CS*TEMP + SN*G(J2,J2)
            G(J2,J2) = CS*G(J2,J2) - SN*TEMP
            G(J1,J1) = CS*G(J1,J1) + SN*TAU
            CALL DROT( N-J1, G(J1,J2), LDG, G(J2,J2), LDG, CS, SN )
         ELSE
            IF ( N.GT.J1+1 )
     $         CALL DROT( N-J1-1, G(J1,J1+2), LDG, G(J2,J1+2), LDG, CS,
     $                    SN )
            CALL DROT( J1-1, G(1,J1), 1, G(1,J2), 1, CS, SN )
         END IF
         IF ( WANTU ) THEN
C
C           Accumulate transformation in the matrices U1 and U2.
C
            CALL DROT( N, U1(1,J1), 1, U1(1,J2), 1, CS, SN )
            CALL DROT( N, U2(1,J1), 1, U2(1,J2), 1, CS, SN )
         END IF
C
      ELSE
C
C        Swapping involves at least one 2-by-2 block.
C
C        Copy the diagonal block of order N1+N2 to the local array D
C        and compute its norm.
C
         ND = N1 + N2
         CALL DLACPY( 'Full', ND, ND, A(J1,J1), LDA, D, LDD )
         DNORM = DLANGE( 'Max', ND, ND, D, LDD, DWORK )
C
C        Compute machine-dependent threshold for test for accepting
C        swap.
C
         EPS = DLAMCH( 'P' )
         SMLNUM = DLAMCH( 'S' ) / EPS
         THRESH = MAX( THIRTY*EPS*DNORM, SMLNUM )
C
C        Solve A11*X - X*A22 = scale*A12 for X.
C
         CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
     $                D(N1+1,N1+1), LDD, D(1,N1+1), LDD, SCALE, X, LDX,
     $                XNORM, IERR )
C
C        Swap the adjacent diagonal blocks.
C
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
C
   10    CONTINUE
C
C        N1 = 1, N2 = 2: generate elementary reflector H so that:
C
C        ( scale, X11, X12 ) H = ( 0, 0, * ).
C
         V(1) = SCALE
         V(2) = X(1,1)
         V(3) = X(1,2)
         CALL DLARFG( 3, V(3), V, 1, TAU )
         V(3) = ONE
         A11 = A(J1,J1)
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL DLARFX( 'Left', 3, 3, V, TAU, D, LDD, DWORK )
         CALL DLARFX( 'Right', 3, 3, V, TAU, D, LDD, DWORK )
C
C        Test whether to reject swap.
C
         IF ( MAX( ABS( D(3,1) ), ABS( D(3,2) ), ABS( D(3,3)-A11 ) )
     $        .GT.THRESH ) GO TO 50
C
C        Accept swap: apply transformation to the entire matrix A.
C
         CALL DLARFX( 'Left', 3, N-J1+1, V, TAU, A(J1,J1), LDA, DWORK )
         CALL DLARFX( 'Right', J2, 3, V, TAU, A(1,J1), LDA, DWORK )
C
         A(J3,J1) = ZERO
         A(J3,J2) = ZERO
         A(J3,J3) = A11
C
C        Apply transformation to G.
C
         IF ( ISHAM ) THEN
            CALL DLARFX( 'Right', J1-1, 3, V, TAU, G(1,J1), LDG, DWORK )
            CALL DSYMV( 'Upper', 3, TAU, G(J1,J1), LDG, V, 1, ZERO,
     $                  DWORK, 1 )
            TEMP = -HALF*TAU*DDOT( 3, DWORK, 1, V, 1 )
            CALL DAXPY( 3, TEMP, V, 1, DWORK, 1 )
            CALL DSYR2( 'Upper', 3, -ONE, V, 1, DWORK, 1,
     $                  G(J1,J1), LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V, TAU, G(J1,J1+3), LDG,
     $                      DWORK )
         ELSE
            CALL DLARFX( 'Right', J1-1, 3, V, TAU, G(1,J1), LDG, DWORK )
            CALL MB01MD( 'Upper', 3, TAU, G(J1,J1), LDG, V, 1, ZERO,
     $                   DWORK, 1 )
            CALL MB01ND( 'Upper', 3, ONE, V, 1, DWORK, 1, G(J1,J1),
     $                   LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V, TAU, G(J1,J1+3), LDG,
     $                      DWORK )
         END IF
C
         IF ( WANTU ) THEN
C
C           Accumulate transformation in the matrices U1 and U2.
C
            CALL DLARFX( 'R', N, 3, V, TAU, U1(1,J1), LDU1, DWORK )
            CALL DLARFX( 'R', N, 3, V, TAU, U2(1,J1), LDU2, DWORK )
         END IF
         GO TO 40
C
   20    CONTINUE
C
C        N1 = 2, N2 = 1: generate elementary reflector H so that:
C
C        H (  -X11 ) = ( * )
C          (  -X21 ) = ( 0 ).
C          ( scale ) = ( 0 )
C
         V(1) = -X(1,1)
         V(2) = -X(2,1)
         V(3) = SCALE
         CALL DLARFG( 3, V(1), V(2), 1, TAU )
         V(1) = ONE
         A33 = A(J3,J3)
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL DLARFX( 'L', 3, 3, V, TAU, D, LDD, DWORK )
         CALL DLARFX( 'R', 3, 3, V, TAU, D, LDD, DWORK )
C
C        Test whether to reject swap.
C
         IF ( MAX( ABS( D(2,1) ), ABS( D(3,1) ), ABS( D(1,1)-A33 ) )
     $       .GT. THRESH ) GO TO 50
C
C        Accept swap: apply transformation to the entire matrix A.
C
         CALL DLARFX( 'Right', J3, 3, V, TAU, A(1,J1), LDA, DWORK )
         CALL DLARFX( 'Left', 3, N-J1, V, TAU, A(J1,J2), LDA, DWORK )
C
         A(J1,J1) = A33
         A(J2,J1) = ZERO
         A(J3,J1) = ZERO
C
C        Apply transformation to G.
C
         IF ( ISHAM ) THEN
            CALL DLARFX( 'Right', J1-1, 3, V, TAU, G(1,J1), LDG, DWORK )
            CALL DSYMV( 'Upper', 3, TAU, G(J1,J1), LDG, V, 1, ZERO,
     $                  DWORK, 1 )
            TEMP = -HALF*TAU*DDOT( 3, DWORK, 1, V, 1 )
            CALL DAXPY( 3, TEMP, V, 1, DWORK, 1 )
            CALL DSYR2( 'Upper', 3, -ONE, V, 1, DWORK, 1, G(J1,J1),
     $                  LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V, TAU, G(J1,J1+3), LDG,
     $                      DWORK )
         ELSE
            CALL DLARFX( 'Right', J1-1, 3, V, TAU, G(1,J1), LDG, DWORK )
            CALL MB01MD( 'Upper', 3, TAU, G(J1,J1), LDG, V, 1, ZERO,
     $                   DWORK, 1 )
            CALL MB01ND( 'Upper', 3, ONE, V, 1, DWORK, 1, G(J1,J1),
     $                   LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V, TAU, G(J1,J1+3), LDG,
     $                      DWORK )
         END IF
C
         IF ( WANTU ) THEN
C
C           Accumulate transformation in the matrices U1 and U2.
C
            CALL DLARFX( 'R', N, 3, V, TAU, U1(1,J1), LDU1, DWORK )
            CALL DLARFX( 'R', N, 3, V, TAU, U2(1,J1), LDU2, DWORK )
         END IF
         GO TO 40
C
   30    CONTINUE
C
C        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
C        that:
C
C        H(2) H(1) (  -X11  -X12 ) = (  *  * )
C                  (  -X21  -X22 )   (  0  * ).
C                  ( scale    0  )   (  0  0 )
C                  (    0  scale )   (  0  0 )
C
         V1(1) = -X(1,1)
         V1(2) = -X(2,1)
         V1(3) = SCALE
         CALL DLARFG( 3, V1(1), V1(2), 1, TAU1 )
         V1(1) = ONE
C
         TEMP  = -TAU1*( X(1,2)+V1(2)*X(2,2) )
         V2(1) = -TEMP*V1(2) - X(2,2)
         V2(2) = -TEMP*V1(3)
         V2(3) = SCALE
         CALL DLARFG( 3, V2(1), V2(2), 1, TAU2 )
         V2(1) = ONE
C
C        Perform swap provisionally on diagonal block in D.
C
         CALL DLARFX( 'L', 3, 4, V1, TAU1, D, LDD, DWORK )
         CALL DLARFX( 'R', 4, 3, V1, TAU1, D, LDD, DWORK )
         CALL DLARFX( 'L', 3, 4, V2, TAU2, D(2,1), LDD, DWORK )
         CALL DLARFX( 'R', 4, 3, V2, TAU2, D(1,2), LDD, DWORK )
C
C        Test whether to reject swap.
C
         IF ( MAX( ABS( D(3,1) ), ABS( D(3,2) ), ABS( D(4,1) ),
     $             ABS( D(4,2) ) ).GT.THRESH ) GO TO 50
C
C        Accept swap: apply transformation to the entire matrix A.
C
         CALL DLARFX( 'L', 3, N-J1+1, V1, TAU1, A(J1,J1), LDA, DWORK )
         CALL DLARFX( 'R', J4, 3, V1, TAU1, A(1,J1), LDA, DWORK )
         CALL DLARFX( 'L', 3, N-J1+1, V2, TAU2, A(J2,J1), LDA, DWORK )
         CALL DLARFX( 'R', J4, 3, V2, TAU2, A(1,J2), LDA, DWORK )
C
         A(J3,J1) = ZERO
         A(J3,J2) = ZERO
         A(J4,J1) = ZERO
         A(J4,J2) = ZERO
C
C        Apply transformation to G.
C
         IF ( ISHAM ) THEN
            CALL DLARFX( 'Right', J1-1, 3, V1, TAU1, G(1,J1), LDG,
     $                   DWORK )
            CALL DSYMV( 'Upper', 3, TAU1, G(J1,J1), LDG, V1, 1, ZERO,
     $                  DWORK, 1 )
            TEMP = -HALF*TAU1*DDOT( 3, DWORK, 1, V1, 1 )
            CALL DAXPY( 3, TEMP, V1, 1, DWORK, 1 )
            CALL DSYR2( 'Upper', 3, -ONE, V1, 1, DWORK, 1,
     $                  G(J1,J1), LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V1, TAU1, G(J1,J1+3),
     $                      LDG, DWORK )
C
            CALL DLARFX( 'Right', J2-1, 3, V2, TAU2, G(1,J2), LDG,
     $                   DWORK )
            CALL DSYMV( 'Upper', 3, TAU2, G(J2,J2), LDG, V2, 1, ZERO,
     $                  DWORK, 1 )
            TEMP = -HALF*TAU2*DDOT( 3, DWORK, 1, V2, 1 )
            CALL DAXPY( 3, TEMP, V2, 1, DWORK, 1 )
            CALL DSYR2( 'Upper', 3, -ONE, V2, 1, DWORK, 1, G(J2,J2),
     $                  LDG )
            IF ( N.GT.J2+2 )
     $         CALL DLARFX( 'Left', 3, N-J2-2, V2, TAU2, G(J2,J2+3),
     $                      LDG, DWORK )
         ELSE
            CALL DLARFX( 'Right', J1-1, 3, V1, TAU1, G(1,J1), LDG,
     $                   DWORK )
            CALL MB01MD( 'Upper', 3, TAU1, G(J1,J1), LDG, V1, 1, ZERO,
     $                   DWORK, 1 )
            CALL MB01ND( 'Upper', 3, ONE, V1, 1, DWORK, 1, G(J1,J1),
     $                   LDG )
            IF ( N.GT.J1+2 )
     $         CALL DLARFX( 'Left', 3, N-J1-2, V1, TAU1, G(J1,J1+3),
     $                      LDG, DWORK )
            CALL DLARFX( 'Right', J2-1, 3, V2, TAU2, G(1,J2), LDG,
     $                   DWORK )
            CALL MB01MD( 'Upper', 3, TAU2, G(J2,J2), LDG, V2, 1, ZERO,
     $                   DWORK, 1 )
            CALL MB01ND( 'Upper', 3, ONE, V2, 1, DWORK, 1, G(J2,J2),
     $                   LDG )
            IF ( N.GT.J2+2 )
     $         CALL DLARFX( 'Left', 3, N-J2-2, V2, TAU2, G(J2,J2+3),
     $                      LDG, DWORK )
         END IF
C
         IF ( WANTU ) THEN
C
C           Accumulate transformation in the matrices U1 and U2.
C
            CALL DLARFX( 'R', N, 3, V1, TAU1, U1(1,J1), LDU1, DWORK )
            CALL DLARFX( 'R', N, 3, V2, TAU2, U1(1,J2), LDU1, DWORK )
            CALL DLARFX( 'R', N, 3, V1, TAU1, U2(1,J1), LDU2, DWORK )
            CALL DLARFX( 'R', N, 3, V2, TAU2, U2(1,J2), LDU2, DWORK )
         END IF
C
   40    CONTINUE
C
         IF ( N2.EQ.2 ) THEN
C
C           Standardize new 2-by-2 block A11.
C
            CALL DLANV2( A(J1,J1), A(J1,J2), A(J2,J1), A(J2,J2), WR1,
     $                   WI1, WR2, WI2, CS, SN )
            CALL DROT( N-J1-1, A(J1,J1+2), LDA, A(J2,J1+2), LDA, CS,
     $                 SN )
            CALL DROT( J1-1, A(1,J1), 1, A(1,J2), 1, CS, SN )
            IF ( ISHAM ) THEN
               TEMP = G(J1,J2)
               CALL DROT( J1, G(1,J1), 1, G(1,J2), 1, CS, SN )
               TAU = CS*TEMP + SN*G(J2,J2)
               G(J2,J2) = CS*G(J2,J2) - SN*TEMP
               G(J1,J1) = CS*G(J1,J1) + SN*TAU
               CALL DROT( N-J1, G(J1,J2), LDG, G(J2,J2), LDG, CS, SN )
            ELSE
               IF ( N.GT.J1+1 )
     $            CALL DROT( N-J1-1, G(J1,J1+2), LDG, G(J2,J1+2), LDG,
     $                       CS, SN )
               CALL DROT( J1-1, G(1,J1), 1, G(1,J2), 1, CS, SN )
            END IF
            IF ( WANTU ) THEN
               CALL DROT( N, U1(1,J1), 1, U1(1,J2), 1, CS, SN )
               CALL DROT( N, U2(1,J1), 1, U2(1,J2), 1, CS, SN )
            END IF
         END IF
C
         IF ( N1.EQ.2 ) THEN
C
C           Standardize new 2-by-2 block A22.
C
            J3 = J1 + N2
            J4 = J3 + 1
            CALL DLANV2( A(J3,J3), A(J3,J4), A(J4,J3), A(J4,J4), WR1,
     $                   WI1, WR2, WI2, CS, SN )
            IF ( J3+2.LE.N )
     $         CALL DROT( N-J3-1, A(J3,J3+2), LDA, A(J4,J3+2), LDA, CS,
     $                    SN )
            CALL DROT( J3-1, A(1,J3), 1, A(1,J4), 1, CS, SN )
            IF ( ISHAM ) THEN
               TEMP = G(J3,J4)
               CALL DROT( J3, G(1,J3), 1, G(1,J4), 1, CS, SN )
               TAU = CS*TEMP + SN*G(J4,J4)
               G(J4,J4) = CS*G(J4,J4) - SN*TEMP
               G(J3,J3) = CS*G(J3,J3) + SN*TAU
               CALL DROT( N-J3, G(J3,J4), LDG, G(J4,J4), LDG, CS, SN )
            ELSE
               IF ( N.GT.J3+1 )
     $            CALL DROT( N-J3-1, G(J3,J3+2), LDG, G(J4,J3+2), LDG,
     $                       CS, SN )
               CALL DROT( J3-1, G(1,J3), 1, G(1,J4), 1, CS, SN )
            END IF
            IF ( WANTU ) THEN
               CALL DROT( N, U1(1,J3), 1, U1(1,J4), 1, CS, SN )
               CALL DROT( N, U2(1,J3), 1, U2(1,J4), 1, CS, SN )
            END IF
         END IF
C
      END IF
      RETURN
C
C     Exit with INFO = 1 if swap was rejected.
C
   50 CONTINUE
      INFO = 1
      RETURN
C *** Last line of MB03TS ***
      END
