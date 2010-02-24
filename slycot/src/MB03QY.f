      SUBROUTINE MB03QY( N, L, A, LDA, U, LDU, E1, E2, INFO )
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
C     To compute the eigenvalues of a selected 2-by-2 diagonal block
C     of an upper quasi-triangular matrix, to reduce the selected block
C     to the standard form and to split the block in the case of real
C     eigenvalues by constructing an orthogonal transformation UT.
C     This transformation is applied to A (by similarity) and to
C     another matrix U from the right.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and UT.  N >= 2.
C
C     L       (input) INTEGER
C             Specifies the position of the block.  1 <= L < N.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi-triangular matrix A whose
C             selected 2-by-2 diagonal block is to be processed.
C             On exit, the leading N-by-N part of this array contains
C             the upper quasi-triangular matrix A after its selected
C             block has been splitt and/or put in the LAPACK standard
C             form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             On entry, the leading N-by-N part of this array must
C             contain a transformation matrix U.
C             On exit, the leading N-by-N part of this array contains
C             U*UT, where UT is the transformation matrix used to
C             split and/or standardize the selected block.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= N.
C
C     E1, E2  (output) DOUBLE PRECISION
C             E1 and E2 contain either the real eigenvalues or the real
C             and positive imaginary parts, respectively, of the complex
C             eigenvalues of the selected 2-by-2 diagonal block of A.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Let A1 = ( A(L,L)    A(L,L+1)   )
C              ( A(L+1,L)  A(L+1,L+1) )
C     be the specified 2-by-2 diagonal block of matrix A.
C     If the eigenvalues of A1 are complex, then they are computed and
C     stored in E1 and E2, where the real part is stored in E1 and the
C     positive imaginary part in E2. The 2-by-2 block is reduced if
C     necessary to the standard form, such that A(L,L) = A(L+1,L+1), and
C     A(L,L+1) and A(L+1,L) have oposite signs. If the eigenvalues are
C     real, the 2-by-2 block is reduced to an upper triangular form such
C     that ABS(A(L,L)) >= ABS(A(L+1,L+1)).
C     In both cases, an orthogonal rotation U1' is constructed such that
C     U1'*A1*U1 has the appropriate form. Let UT be an extension of U1
C     to an N-by-N orthogonal matrix, using identity submatrices. Then A
C     is replaced by UT'*A*UT and the contents of array U is U * UT.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     March 1998. Based on the RASP routine SPLITB.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalues, orthogonal transformation, real Schur form,
C     similarity transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER        ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER          INFO, L, LDA, LDU, N
      DOUBLE PRECISION E1, E2
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), U(LDU,*)
C     .. Local Scalars ..
      INTEGER          L1
      DOUBLE PRECISION EW1, EW2, CS, SN
C     .. External Subroutines ..
      EXTERNAL         DLANV2, DROT, XERBLA
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.2 ) THEN
         INFO = -1
      ELSE IF( L.LT.1 .OR. L.GE.N ) THEN
         INFO = -2
      ELSE IF( LDA.LT.N ) THEN
         INFO = -4
      ELSE IF( LDU.LT.N ) THEN
         INFO = -6
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03QY', -INFO )
         RETURN
      END IF
C
C     Compute the eigenvalues and the elements of the Givens
C     transformation.
C
      L1 = L + 1
      CALL DLANV2( A(L,L), A(L,L1), A(L1,L), A(L1,L1), E1, E2,
     $             EW1, EW2, CS, SN )
      IF( E2.EQ.ZERO ) E2 = EW1
C
C     Apply the transformation to A.
C
      IF( L1.LT.N )
     $   CALL DROT( N-L1, A(L,L1+1), LDA, A(L1,L1+1), LDA, CS, SN )
      CALL DROT( L-1, A(1,L), 1, A(1,L1), 1, CS, SN )
C
C     Accumulate the transformation in U.
C
      CALL DROT( N, U(1,L), 1, U(1,L1), 1, CS, SN )
C
      RETURN
C *** Last line of MB03QY ***
      END
