      SUBROUTINE MB04TW( UPDATQ, M, N, NRE, NCE, IFIRE, IFICE, IFICA, A,
     $                   LDA, E, LDE, Q, LDQ )
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
C     To reduce a submatrix E(k) of E to upper triangular form by row
C     Givens rotations only.
C     Here E(k) = E(IFIRE:me,IFICE:ne), where me = IFIRE - 1 + NRE,
C                                             ne = IFICE - 1 + NCE.
C     Matrix E(k) is assumed to have full column rank on entry. Hence,
C     no pivoting is done during the reduction process. See Algorithm
C     2.3.1 and Remark 2.3.4 in [1].
C     The constructed row transformations are also applied to matrix
C     A(k) = A(IFIRE:me,IFICA:N).
C     Note that in A(k) rows are transformed with the same row indices
C     as in E but with column indices different from those in E.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     UPDATQ  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix Q the orthogonal row transformations, as follows:
C             = .FALSE.: Do not form Q;
C             = .TRUE.:  The given matrix Q is updated by the orthogonal
C                        row transformations used in the reduction.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             Number of rows of A and E.  M >= 0.
C
C     N       (input) INTEGER
C             Number of columns of A and E.  N >= 0.
C
C     NRE     (input) INTEGER
C             Number of rows in E to be transformed.  0 <= NRE <= M.
C
C     NCE     (input) INTEGER
C             Number of columns in E to be transformed.  0 <= NCE <= N.
C
C     IFIRE   (input) INTEGER
C             Index of first row in E to be transformed.
C
C     IFICE   (input) INTEGER
C             Index of first column in E to be transformed.
C
C     IFICA   (input) INTEGER
C             Index of first column in A to be transformed.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, this array contains the submatrix A(k).
C             On exit, it contains the transformed matrix A(k).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, this array contains the submatrix E(k) of full
C             column rank to be reduced to upper triangular form.
C             On exit, it contains the transformed matrix E.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,M).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*)
C             On entry, if UPDATQ = .TRUE., then the leading M-by-M
C             part of this array must contain a given matrix Q (e.g.
C             from a previous call to another SLICOT routine), and on
C             exit, the leading M-by-M part of this array contains the
C             product of the input matrix Q and the row transformation
C             matrix that has transformed the rows of the matrices A
C             and E.
C             If UPDATQ = .FALSE., the array Q is not referenced and
C             can be supplied as a dummy array (i.e. set parameter
C             LDQ = 1 and declare this array to be Q(1,1) in the calling
C             program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q. If UPDATQ = .TRUE.,
C             LDQ >= MAX(1,M); if UPDATQ = .FALSE., LDQ >= 1.
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
C     Supersedes Release 2.0 routine MB04FW by Th.G.J. Beelen,
C     Philips Glass Eindhoven, Holland.
C
C     REVISIONS
C
C     June 13, 1997. V. Sima.
C     December 30, 1997. A. Varga: Corrected column range to apply
C                                  transformations on the matrix E.
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
      LOGICAL           UPDATQ
      INTEGER           IFICA, IFICE, IFIRE, LDA, LDE, LDQ, M, N, NCE,
     $                  NRE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Q(LDQ,*)
C     .. Local Scalars ..
      INTEGER           I, IPVT, J
      DOUBLE PRECISION  SC, SS
C     .. External Subroutines ..
      EXTERNAL          DROT, DROTG
C     .. Executable Statements ..
C
      IF ( M.LE.0 .OR. N.LE.0 .OR. NRE.LE.0 .OR. NCE.LE.0 )
     $   RETURN
C
      IPVT = IFIRE - 1
C
      DO 40 J = IFICE, IFICE + NCE - 1
         IPVT = IPVT + 1
C
         DO 20 I = IPVT + 1, IFIRE + NRE - 1
C
C           Determine the Givens transformation on rows i and ipvt
C           to annihilate E(i,j).
C           Apply the transformation to these rows (in whole E-matrix)
C           from columns j up to n .
C           Apply the transformations also to the A-matrix
C           (from columns ifica up to n).
C           Update the row transformation matrix Q, if needed.
C
            CALL DROTG( E(IPVT,J), E(I,J), SC, SS )
            CALL DROT( N-J, E(IPVT,J+1), LDE, E(I,J+1), LDE, SC, SS )
            E(I,J) = ZERO
            CALL DROT( N-IFICA+1, A(IPVT,IFICA), LDA, A(I,IFICA), LDA,
     $                 SC, SS )
            IF( UPDATQ )
     $         CALL DROT( M, Q(1,IPVT), 1, Q(1,I), 1, SC, SS )
   20    CONTINUE
C
   40 CONTINUE
C
      RETURN
C *** Last line of MB04TW ***
      END
