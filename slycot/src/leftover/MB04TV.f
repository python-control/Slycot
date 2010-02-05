      SUBROUTINE MB04TV( UPDATZ, N, NRA, NCA, IFIRA, IFICA, A, LDA, E,
     $                   LDE, Z, LDZ )
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
C     To reduce a submatrix A(k) of A to upper triangular form by column
C     Givens rotations only.
C     Here A(k) = A(IFIRA:ma,IFICA:na) where ma = IFIRA - 1 + NRA,
C     na = IFICA - 1 + NCA.
C     Matrix A(k) is assumed to have full row rank on entry. Hence, no
C     pivoting is done during the reduction process. See Algorithm 2.3.1
C     and Remark 2.3.4 in [1].
C     The constructed column transformations are also applied to matrix
C     E(k) = E(1:IFIRA-1,IFICA:na).
C     Note that in E columns are transformed with the same column
C     indices as in A, but with row indices different from those in A.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPDATZ  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the orthogonal column transformations, as
C             follows:
C             = .FALSE.: Do not form Z;
C             = .TRUE.:  The given matrix Z is updated by the orthogonal
C                        column transformations used in the reduction.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             Number of columns of A and E.  N >= 0.
C
C     NRA     (input) INTEGER
C             Number of rows in A to be transformed.  0 <= NRA <= LDA.
C
C     NCA     (input) INTEGER
C             Number of columns in A to be transformed.  0 <= NCA <= N.
C
C     IFIRA   (input) INTEGER
C             Index of the first row in A to be transformed.
C
C     IFICA   (input) INTEGER
C             Index of the first column in A to be transformed.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the elements of A(IFIRA:ma,IFICA:na) must
C             contain the submatrix A(k) of full row rank to be reduced
C             to upper triangular form.
C             On exit, it contains the transformed matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,NRA).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the elements of E(1:IFIRA-1,IFICA:na) must
C             contain the submatrix E(k).
C             On exit, it contains the transformed matrix E.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,IFIRA-1).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*)
C             On entry, if UPDATZ = .TRUE., then the leading N-by-N
C             part of this array must contain a given matrix Z (e.g.
C             from a previous call to another SLICOT routine), and on
C             exit, the leading N-by-N part of this array contains the
C             product of the input matrix Z and the column
C             transformation matrix that has transformed the columns of
C             the matrices A and E.
C             If UPDATZ = .FALSE., the array Z is not referenced and
C             can be supplied as a dummy array (i.e. set parameter
C             LDZ = 1 and declare this array to be Z(1,1) in the calling
C             program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If UPDATZ = .TRUE.,
C             LDZ >= MAX(1,N); if UPDATZ = .FALSE., LDZ >= 1.
C
C     REFERENCES
C
C     [1] Beelen, Th.
C         New Algorithms for Computing the Kronecker structure of a
C         Pencil with Applications to Systems and Control Theory.
C         Ph.D.Thesis, Eindhoven University of Technology,
C         The Netherlands, 1987.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB04FV by Th.G.J. Beelen,
C     Philips Glass Eindhoven, Holland.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, orthogonal transformation,
C     staircase form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           UPDATZ
      INTEGER           IFICA, IFIRA, LDA, LDE, LDZ, N, NCA, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER           I, IFIRA1, J, JPVT
      DOUBLE PRECISION  SC, SS
C     .. External Subroutines ..
      EXTERNAL          DROT, DROTG
C     .. Executable Statements ..
C
      IF ( N.LE.0 .OR. NRA.LE.0 .OR. NCA.LE.0 )
     $   RETURN
      IFIRA1 = IFIRA - 1
      JPVT = IFICA + NCA
C
      DO 40 I = IFIRA1 + NRA, IFIRA, -1
         JPVT = JPVT - 1
C
         DO 20 J = JPVT - 1, IFICA, -1
C
C           Determine the Givens transformation on columns j and jpvt
C           to annihilate A(i,j). Apply the transformation to these
C           columns from rows 1 up to i.
C           Apply the transformation also to the E-matrix (from rows 1
C           up to ifira1).
C           Update column transformation matrix Z, if needed.
C
            CALL DROTG( A(I,JPVT), A(I,J), SC, SS )
            CALL DROT( I-1, A(1,JPVT), 1, A(1,J), 1, SC, SS )
            A(I,J) = ZERO
            CALL DROT( IFIRA1, E(1,JPVT), 1, E(1,J), 1, SC, SS )
            IF( UPDATZ ) CALL DROT( N, Z(1,JPVT), 1, Z(1,J), 1, SC, SS )
   20    CONTINUE
C
   40 CONTINUE
C
      RETURN
C *** Last line of MB04TV ***
      END
