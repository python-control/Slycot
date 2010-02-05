      SUBROUTINE MB04TY( UPDATQ, UPDATZ, M, N, NBLCKS, INUK, IMUK, A,
     $                   LDA, E, LDE, Q, LDQ, Z, LDZ, INFO )
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
C     To perform the triangularization of the submatrices having full
C     row and column rank in the pencil s*E(eps,inf)-A(eps,inf) below
C
C                    | s*E(eps,inf)-A(eps,inf) |     X       |
C          s*E - A = |-------------------------|-------------| ,
C                    |            0            | s*E(r)-A(r) |
C
C     using Algorithm 3.3.1 in [1].
C     On entry, it is assumed that the M-by-N matrices A and E have
C     been transformed to generalized Schur form by unitary
C     transformations (see Algorithm 3.2.1 in [1]), and that the pencil
C     s*E(eps,inf)-A(eps,inf) is in staircase form.
C     This pencil contains all Kronecker column indices and infinite
C     elementary divisors of the pencil s*E - A.
C     The pencil s*E(r)-A(r) contains all Kronecker row indices and
C     finite elementary divisors of s*E - A.
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
C     M       (input) INTEGER
C             Number of rows in A and E.  M >= 0.
C
C     N       (input) INTEGER
C             Number of columns in A and E.  N >= 0.
C
C     NBLCKS  (input) INTEGER
C             Number of submatrices having full row rank (possibly zero)
C             in A(eps,inf).
C
C     INUK    (input) INTEGER array, dimension (NBLCKS)
C             The row dimensions nu(k) (k=1, 2, ..., NBLCKS) of the
C             submatrices having full row rank in the pencil
C             s*E(eps,inf)-A(eps,inf).
C
C     IMUK    (input) INTEGER array, dimension (NBLCKS)
C             The column dimensions mu(k) (k=1, 2, ..., NBLCKS) of the
C             submatrices having full column rank in the pencil.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, this array contains the matrix A to be reduced.
C             On exit, it contains the transformed matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, this array contains the matrix E to be reduced.
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
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  if incorrect dimensions of a full column rank
C                   submatrix;
C             = 2:  if incorrect dimensions of a full row rank
C                   submatrix.
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
C     Supersedes Release 2.0 routine MB04FY by Th.G.J. Beelen,
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
C     .. Scalar Arguments ..
      LOGICAL           UPDATQ, UPDATZ
      INTEGER           INFO, LDA, LDE, LDQ, LDZ, M, N, NBLCKS
C     .. Array Arguments ..
      INTEGER           IMUK(*), INUK(*)
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER           IFICA, IFICE, IFIRE, ISMUK, ISNUK1, K, MUK,
     $                  MUKP1, NUK
C     .. External Subroutines ..
      EXTERNAL          MB04TV, MB04TW
C     .. Executable Statements ..
C
      INFO = 0
      IF ( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
C     ISMUK  = sum(i=1,...,k) MU(i),
C     ISNUK1 = sum(i=1,...,k-1) NU(i).
C
      ISMUK = 0
      ISNUK1 = 0
C
      DO 20 K = 1, NBLCKS
         ISMUK  = ISMUK  + IMUK(K)
         ISNUK1 = ISNUK1 + INUK(K)
   20 CONTINUE
C
C     Note:  ISNUK1 has not yet the correct value.
C
      MUKP1 = 0
C
      DO 40 K = NBLCKS, 1, -1
         MUK = IMUK(K)
         NUK = INUK(K)
         ISNUK1 = ISNUK1 - NUK
C
C        Determine left upper absolute co-ordinates of E(k) in E-matrix
C        and of A(k) in A-matrix.
C
         IFIRE = 1 + ISNUK1
         IFICE = 1 + ISMUK
         IFICA = IFICE - MUK
C
C        Reduce E(k) to upper triangular form using Givens
C        transformations on rows only. Apply the same transformations
C        to the rows of A(k).
C
         IF ( MUKP1.GT.NUK ) THEN
            INFO = 1
            RETURN
         END IF
C
         CALL MB04TW( UPDATQ, M, N, NUK, MUKP1, IFIRE, IFICE, IFICA, A,
     $                LDA, E, LDE, Q, LDQ )
C
C        Reduce A(k) to upper triangular form using Givens
C        transformations on columns only. Apply the same transformations
C        to the columns in the E-matrix.
C
         IF ( NUK.GT.MUK ) THEN
            INFO = 2
            RETURN
         END IF
C
         CALL MB04TV( UPDATZ, N, NUK, MUK, IFIRE, IFICA, A, LDA, E, LDE,
     $                Z, LDZ )
C
         ISMUK = ISMUK - MUK
         MUKP1 = MUK
   40 CONTINUE
C
      RETURN
C *** Last line of MB04TY ***
      END
