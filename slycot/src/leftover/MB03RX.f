      SUBROUTINE MB03RX( JOBV, N, KL, KU, A, LDA, X, LDX, WR, WI,
     $                   DWORK )
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
C     To reorder the diagonal blocks of the principal submatrix between
C     the indices KL and KU (KU >= KL) of a real Schur form matrix A
C     together with their eigenvalues, using orthogonal similarity
C     transformations, such that the block specified by KU is moved in
C     the position KL. The transformations are optionally postmultiplied
C     in a given matrix X.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBV    CHARACTER*1
C             Specifies whether or not the transformations are
C             accumulated, as follows:
C             = 'N':  The transformations are not accumulated;
C             = 'V':  The transformations are accumulated in X (the
C                     given matrix X is updated).
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and X.  N >= 0.
C
C     KL      (input) INTEGER
C             The lower boundary index for the rows and columns of the
C             principal submatrix of A whose diagonal blocks are to be
C             reordered, and also the target position for the block to
C             be moved.  1 <= KL <= KU <= N.
C
C     KU      (input/output) INTEGER
C             On entry, KU specifies the upper boundary index for the
C             rows and columns of the principal submatrix of A whose
C             diagonal blocks are to be reordered, and also the original
C             position for the block to be moved.  1 <= KL <= KU <= N.
C             On exit, KU specifies the upper boundary index for the
C             rows and columns of the principal submatrix of A whose
C             diagonal blocks have been reordered. The given value will
C             be increased by 1 if the moved block was 2-by-2 and it has
C             been replaced by two 1-by-1 blocks. Otherwise, its input
C             value is preserved.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A in real Schur canonical form.
C             On exit, the leading N-by-N part of this array contains
C             the ordered real Schur canonical form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N)
C             On entry, if JOBV = 'V', the leading N-by-N part of this
C             array must contain a given matrix X.
C             On exit, if JOBV = 'V', the leading N-by-N part of this
C             array contains the product of the given matrix X and the
C             transformation matrix that performed the reordering of A.
C             If JOBV = 'N', this array is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of array X.
C             LDX >= 1,        if JOBV = 'N';
C             LDX >= MAX(1,N), if JOBV = 'V'.
C
C     WR,     (input/output) DOUBLE PRECISION arrays, dimension (N)
C     WI      On entry, these arrays must contain the real and imaginary
C             parts, respectively, of the eigenvalues of the matrix A.
C             On exit, these arrays contain the real and imaginary
C             parts, respectively, of the eigenvalues of the matrix A,
C             possibly reordered.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     METHOD
C
C     An attempt is made to move the block in the position (KU,KU) to
C     the position (KL,KL) by a sequence of orthogonal similarity
C     transformations, each swapping two consecutive blocks. The
C     standard algorithm [1], [2] usually succeeds to perform this
C     reordering. A failure of this algorithm means that two consecutive
C     blocks (one of them being the desired block possibly moved) are
C     too close to swap. In such a case, the leading block of the two
C     is tried to be moved in the position (KL,KL) and the procedure is
C     repeated.
C
C     REFERENCES
C
C     [1] Stewart, G.W.
C         HQR3 and EXCHQZ: FORTRAN subroutines for calculating and
C         ordering the eigenvalues of a real upper Hessenberg matrix.
C         ACM TOMS, 2, pp. 275-280, 1976.
C
C     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable. If some eigenvalues are
C     ill-conditioned, their returned values could differ much from
C     their input values.
C
C     CONTRIBUTOR
C
C     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalue, orthogonal transformation, real Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBV
      INTEGER           KL, KU, LDA, LDX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), WI(*), WR(*), X(LDX,*)
C     .. Local Scalars ..
      INTEGER           IERR, IFST, ILST, L
C     .. External Subroutines ..
      EXTERNAL          DTREXC
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
      IF ( KU.GT.KL ) THEN
C
C        Try to move the block in position (KU,KU) to position (KL,KL).
C
         IFST = KU
C        REPEAT
   10    CONTINUE
         ILST = KL
         CALL DTREXC( JOBV, N, A, LDA, X, LDX, IFST, ILST, DWORK, IERR )
         IF ( IERR.NE.0 ) THEN
C
C           During calculations, two adjacent blocks were too close
C           to swap; the desired block cannot be moved further, but the
C           block above it is suitable and is tried for moving. The
C           number of repeat cycles is usually 1, and at most the number
C           of blocks between the current position and the position KL.
C
            IFST = ILST - 1
            IF ( IFST.GT.1 ) THEN
               IF ( A(IFST,IFST-1).NE.ZERO )
     $            IFST = ILST - 2
            END IF
            IF ( ILST.GT.KL )
     $         GO TO 10
         END IF
C        UNTIL ( ILST.EQ.KL on output from DTREXC )
C
C        Recompute the eigenvalues for the modified part of A.
C        Note that KU must be incremented if the moved block was 2-by-2
C        and it has been replaced by two 1-by-1 blocks.
C
         IF ( WI(KU).NE.ZERO ) THEN
            IF ( A(KU+1,KU).EQ.ZERO )
     $         KU = KU + 1
         END IF
C
         L = KL
C        WHILE ( L.LT.KU .OR. ( L.EQ.KU .AND. L.LT.N ) ) DO
   20    IF ( L.LT.KU .OR. ( L.EQ.KU .AND. L.LT.N ) ) THEN
            IF ( A(L+1,L).NE.ZERO ) THEN
C
C              A 2x2 block.
C
               WR(L)   =  A(L,L)
               WR(L+1) =  WR(L)
               WI(L)   =  SQRT( ABS( A(L,L+1) ) )*
     $                    SQRT( ABS( A(L+1,L) ) )
               WI(L+1) = -WI(L)
               L = L + 2
            ELSE
C
C              An 1x1 block.
C
               WR(L) = A(L,L)
               WI(L) = ZERO
               L = L + 1
            END IF
            GO TO 20
         ELSE IF ( L.EQ.N ) THEN
            WR(L) = A(L,L)
            WI(L) = ZERO
         END IF
C        END WHILE 20
      END IF
C
      RETURN
C *** Last line of MB03RX ***
      END
