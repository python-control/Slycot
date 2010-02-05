      SUBROUTINE MB04TX( UPDATQ, UPDATZ, M, N, NBLCKS, INUK, IMUK, A,
     $                   LDA, E, LDE, Q, LDQ, Z, LDZ, MNEI )
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
C     To separate the pencils s*E(eps)-A(eps) and s*E(inf)-A(inf) in
C     s*E(eps,inf)-A(eps,inf) using Algorithm 3.3.3 in [1].
C
C     On entry, it is assumed that the M-by-N matrices A and E have
C     been obtained after applying the Algorithms 3.2.1 and 3.3.1 to
C     the pencil s*E - A as described in [1], i.e.
C
C                        | s*E(eps,inf)-A(eps,inf) |      X      |
C        Q'(s*E - A)Z  = |-------------------------|-------------|
C                        |             0           | s*E(r)-A(r) |
C
C     Here the pencil s*E(eps,inf)-A(eps,inf) is in staircase form.
C     This pencil contains all Kronecker column indices and infinite
C     elementary divisors of the pencil s*E - A.
C     The pencil s*E(r)-A(r) contains all Kronecker row indices and
C     finite elementary divisors of s*E - A.
C     Furthermore, the submatrices having full row and column rank in
C     the pencil s*E(eps,inf)-A(eps,inf) are assumed to be
C     triangularized.
C
C     On exit, the result then is
C
C                        Q'(s*E - A)Z =
C
C          | s*E(eps)-A(eps) |        X        |      X      |
C          |-----------------|-----------------|-------------|
C          |        0        | s*E(inf)-A(inf) |      X      |
C          |===================================|=============|
C          |                                   |             |
C          |                 0                 | s*E(r)-A(r) |
C
C     Note that the pencil s*E(r)-A(r) is not reduced further.
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
C             Number of rows of A and E.  M >= 0.
C
C     N       (input) INTEGER
C             Number of columns of A and E.  N >= 0.
C
C     NBLCKS  (input/output) INTEGER
C             On entry, the number of submatrices having full row rank
C             (possibly zero) in A(eps,inf).
C             On exit, the input value has been reduced by one, if the
C             last submatrix is a 0-by-0 (empty) matrix.
C
C     INUK    (input/output) INTEGER array, dimension (NBLCKS)
C             On entry, this array contains the row dimensions nu(k),
C             (k=1, 2, ..., NBLCKS) of the submatrices having full row
C             rank in the pencil s*E(eps,inf)-A(eps,inf).
C             On exit, this array contains the row dimensions nu(k),
C             (k=1, 2, ..., NBLCKS) of the submatrices having full row
C             rank in the pencil s*E(eps)-A(eps).
C
C     IMUK    (input/output) INTEGER array, dimension (NBLCKS)
C             On entry, this array contains the column dimensions mu(k),
C             (k=1, 2, ..., NBLCKS) of the submatrices having full
C             column rank in the pencil s*E(eps,inf)-A(eps,inf).
C             On exit, this array contains the column dimensions mu(k),
C             (k=1, 2, ..., NBLCKS) of the submatrices having full
C             column rank in the pencil s*E(eps)-A(eps).
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
C     MNEI    (output) INTEGER array, dimension (4)
C             MNEI(1) = MEPS = row    dimension of s*E(eps)-A(eps),
C             MNEI(2) = NEPS = column dimension of s*E(eps)-A(eps),
C             MNEI(3) = MINF = row    dimension of s*E(inf)-A(inf),
C             MNEI(4) = NINF = column dimension of s*E(inf)-A(inf).
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
C     Supersedes Release 2.0 routine MB04FX by Th.G.J. Beelen,
C     Philips Glass Eindhoven, Holland.
C
C     REVISIONS
C
C     June 13, 1997, V. Sima.
C     November 24, 1997, A. Varga: initialization of MNEI to 0, instead
C                                  of ZERO.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, orthogonal
C     transformation, staircase form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           UPDATQ, UPDATZ
      INTEGER           LDA, LDE, LDQ, LDZ, M, N, NBLCKS
C     .. Array Arguments ..
      INTEGER           IMUK(*), INUK(*), MNEI(4)
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER           CA, CE, CJA, CJE, IP, ISMUK, ISNUK, K, MEPS,
     $                  MINF, MUK, MUKP1, MUP, MUP1, NEPS, NINF, NUK,
     $                  NUP, RA, RJE, SK1P1, TK1P1, TP1
      DOUBLE PRECISION  SC, SS
C     .. External Subroutines ..
      EXTERNAL          DROTG, MB04TU
C     .. Executable Statements ..
C
      MNEI(1) = 0
      MNEI(2) = 0
      MNEI(3) = 0
      MNEI(4) = 0
      IF ( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
C     Initialisation.
C
      ISMUK = 0
      ISNUK = 0
C
      DO 20 K = 1, NBLCKS
         ISMUK = ISMUK + IMUK(K)
         ISNUK = ISNUK + INUK(K)
   20 CONTINUE
C
C     MEPS, NEPS are the dimensions of the pencil s*E(eps)-A(eps).
C     MEPS = Sum(k=1,...,nblcks) NU(k),
C     NEPS = Sum(k=1,...,nblcks) MU(k).
C     MINF, NINF are the dimensions of the pencil s*E(inf)-A(inf).
C
      MEPS = ISNUK
      NEPS = ISMUK
      MINF = 0
      NINF = 0
C
C     MUKP1 = mu(k+1).  N.B. It is assumed that mu(NBLCKS + 1) = 0.
C
      MUKP1 = 0
C
      DO 120 K = NBLCKS, 1, -1
         NUK = INUK(K)
         MUK = IMUK(K)
C
C        Reduce submatrix E(k,k+1) to square matrix.
C        NOTE that always NU(k) >= MU(k+1) >= 0.
C
C        WHILE ( NU(k) >  MU(k+1) ) DO
   40    IF ( NUK.GT.MUKP1 ) THEN
C
C           sk1p1 = sum(i=k+1,...,p-1) NU(i)
C           tk1p1 = sum(i=k+1,...,p-1) MU(i)
C           ismuk = sum(i=1,...,k) MU(i)
C           tp1   = sum(i=1,...,p-1) MU(i) = ismuk + tk1p1.
C
            SK1P1 = 0
            TK1P1 = 0
C
            DO 100 IP = K + 1, NBLCKS
C
C              Annihilate the elements originally present in the last
C              row of E(k,p+1) and A(k,p).
C              Start annihilating the first MU(p) - MU(p+1) elements by
C              applying column Givens rotations plus interchanging
C              elements.
C              Use original bottom diagonal element of A(k,k) as pivot.
C              Start position of pivot in A = (ra,ca).
C
               TP1 = ISMUK + TK1P1
               RA  = ISNUK + SK1P1
               CA  = TP1
C
               MUP  = IMUK(IP)
               NUP  = INUK(IP)
               MUP1 = NUP
C
               DO 60 CJA = CA, CA + MUP - NUP - 1
C
C                 CJA = current column index of pivot in A.
C
                  CALL DROTG( A(RA,CJA), A(RA,CJA+1), SC, SS )
C
C                 Apply transformations to A- and E-matrix.
C                 Interchange columns simultaneously.
C                 Update column transformation matrix Z, if needed.
C
                  CALL MB04TU( RA-1, A(1,CJA), 1, A(1,CJA+1), 1, SC,
     $                         SS )
                  A(RA,CJA+1) = A(RA,CJA)
                  A(RA,CJA) = ZERO
                  CALL MB04TU( RA, E(1,CJA), 1, E(1,CJA+1), 1, SC, SS )
                  IF( UPDATZ ) CALL MB04TU( N, Z(1,CJA), 1, Z(1,CJA+1),
     $                                      1, SC, SS )
   60          CONTINUE
C
C              Annihilate the remaining elements originally present in
C              the last row of E(k,p+1) and A(k,p) by alternatingly
C              applying row and column rotations plus interchanging
C              elements.
C              Use diagonal elements of E(p,p+1) and original bottom
C              diagonal element of A(k,k) as pivots, respectively.
C              (re,ce) and (ra,ca) are the starting positions of the
C              pivots in E and A.
C
               CE = TP1 + MUP
               CA = CE - MUP1 - 1
C
               DO 80 RJE = RA + 1, RA + MUP1
C
C                 (RJE,CJE) = current position pivot in E.
C
                  CJE = CE + 1
                  CJA = CA + 1
C
C                 Determine the row transformations.
C                 Apply these transformations to E- and A-matrix.
C                 Interchange the rows simultaneously.
C                 Update row transformation matrix Q, if needed.
C
                  CALL DROTG( E(RJE,CJE), E(RJE-1,CJE), SC, SS )
                  CALL MB04TU( N-CJE, E(RJE,CJE+1), LDE, E(RJE-1,CJE+1),
     $                         LDE, SC, SS )
                  E(RJE-1,CJE) = E(RJE,CJE)
                  E(RJE,CJE) = ZERO
                  CALL MB04TU( N-CJA+1, A(RJE,CJA), LDA, A(RJE-1,CJA),
     $                         LDA, SC, SS )
                  IF( UPDATQ ) CALL MB04TU( M, Q(1,RJE), 1,
     $                                      Q(1,RJE-1), 1, SC, SS )
C
C                 Determine the column transformations.
C                 Apply these transformations to A- and E-matrix.
C                 Interchange the columns simultaneously.
C                 Update column transformation matrix Z, if needed.
C
                  CALL DROTG( A(RJE,CJA), A(RJE,CJA+1), SC, SS )
                  CALL MB04TU( RJE-1, A(1,CJA), 1, A(1,CJA+1), 1, SC,
     $                         SS )
                  A(RJE,CJA+1) = A(RJE,CJA)
                  A(RJE,CJA) = ZERO
                  CALL MB04TU( RJE, E(1,CJA), 1, E(1,CJA+1), 1, SC, SS )
                  IF( UPDATZ ) CALL MB04TU( N, Z(1,CJA), 1, Z(1,CJA+1),
     $                                      1, SC, SS )
   80          CONTINUE
C
               SK1P1 = SK1P1 + NUP
               TK1P1 = TK1P1 + MUP
C
  100       CONTINUE
C
C           Reduce A=A(eps,inf) and E=E(eps,inf) by ignoring their last
C           row and right most column. The row and column ignored
C           belong to the pencil s*E(inf)-A(inf).
C           Redefine blocks in new A and E.
C
            MUK = MUK - 1
            NUK = NUK - 1
            ISMUK = ISMUK - 1
            ISNUK = ISNUK - 1
            MEPS = MEPS - 1
            NEPS = NEPS - 1
            MINF = MINF + 1
            NINF = NINF + 1
C
            GO TO 40
         END IF
C        END WHILE 40
C
         IMUK(K) = MUK
         INUK(K) = NUK
C
C        Now submatrix E(k,k+1) is square.
C
C        Consider next submatrix (k:=k-1).
C
         ISNUK = ISNUK - NUK
         ISMUK = ISMUK - MUK
         MUKP1 = MUK
  120 CONTINUE
C
C     If mu(NBLCKS) = 0, then the last submatrix counted in NBLCKS is
C     a 0-by-0 (empty) matrix. This "matrix" must be removed.
C
      IF ( IMUK(NBLCKS).EQ.0 ) NBLCKS = NBLCKS - 1
C
C     Store dimensions of the pencils s*E(eps)-A(eps) and
C     s*E(inf)-A(inf) in array MNEI.
C
      MNEI(1) = MEPS
      MNEI(2) = NEPS
      MNEI(3) = MINF
      MNEI(4) = NINF
C
      RETURN
C *** Last line of MB04TX ***
      END
