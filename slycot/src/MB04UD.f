      SUBROUTINE MB04UD( JOBQ, JOBZ, M, N, A, LDA, E, LDE, Q, LDQ,
     $                   Z, LDZ, RANKE, ISTAIR, TOL, DWORK, INFO )
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
C     To compute orthogonal transformations Q and Z such that the
C     transformed pencil Q'(sE-A)Z has the E matrix in column echelon
C     form, where E and A are M-by-N matrices.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBQ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Q the unitary row permutations, as follows:
C             = 'N':  Do not form Q;
C             = 'I':  Q is initialized to the unit matrix and the
C                     unitary row permutation matrix Q is returned;
C             = 'U':  The given matrix Q is updated by the unitary
C                     row permutations used in the reduction.
C
C     JOBZ    CHARACTER*1
C             Indicates whether the user wishes to accumulate in a
C             matrix Z the unitary column transformations, as follows:
C             = 'N':  Do not form Z;
C             = 'I':  Z is initialized to the unit matrix and the
C                     unitary transformation matrix Z is returned;
C             = 'U':  The given matrix Z is updated by the unitary
C                     transformations used in the reduction.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows in the matrices A, E and the order of
C             the matrix Q.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns in the matrices A, E and the order
C             of the matrix Z.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the A matrix of the pencil sE-A.
C             On exit, the leading M-by-N part of this array contains
C             the unitary transformed matrix Q' * A * Z.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading M-by-N part of this array must
C             contain the E matrix of the pencil sE-A, to be reduced to
C             column echelon form.
C             On exit, the leading M-by-N part of this array contains
C             the unitary transformed matrix Q' * E * Z, which is in
C             column echelon form.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,M).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*)
C             On entry, if JOBQ = 'U', then the leading M-by-M part of
C             this array must contain a given matrix Q (e.g. from a
C             previous call to another SLICOT routine), and on exit, the
C             leading M-by-M part of this array contains the product of
C             the input matrix Q and the row permutation matrix used to
C             transform the rows of matrix E.
C             On exit, if JOBQ = 'I', then the leading M-by-M part of
C             this array contains the matrix of accumulated unitary
C             row transformations performed.
C             If JOBQ = 'N', the array Q is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDQ = 1 and
C             declare this array to be Q(1,1) in the calling program).
C
C     LDQ     INTEGER
C             The leading dimension of array Q. If JOBQ = 'U' or
C             JOBQ = 'I', LDQ >= MAX(1,M); if JOBQ = 'N', LDQ >= 1.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*)
C             On entry, if JOBZ = 'U', then the leading N-by-N part of
C             this array must contain a given matrix Z (e.g. from a
C             previous call to another SLICOT routine), and on exit, the
C             leading N-by-N part of this array contains the product of
C             the input matrix Z and the column transformation matrix
C             used to transform the columns of matrix E.
C             On exit, if JOBZ = 'I', then the leading N-by-N part of
C             this array contains the matrix of accumulated unitary
C             column transformations performed.
C             If JOBZ = 'N', the array Z is not referenced and can be
C             supplied as a dummy array (i.e. set parameter LDZ = 1 and
C             declare this array to be Z(1,1) in the calling program).
C
C     LDZ     INTEGER
C             The leading dimension of array Z. If JOBZ = 'U' or
C             JOBZ = 'I', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1.
C
C     RANKE   (output) INTEGER
C             The computed rank of the unitary transformed matrix E.
C
C     ISTAIR  (output) INTEGER array, dimension (M)
C             This array contains information on the column echelon form
C             of the unitary transformed matrix E. Specifically,
C             ISTAIR(i) = +j if the first non-zero element E(i,j)
C             is a corner point and -j otherwise, for i = 1,2,...,M.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance below which matrix elements are considered
C             to be zero. If the user sets TOL to be less than (or
C             equal to) zero then the tolerance is taken as
C             EPS * MAX(ABS(E(I,J))), where EPS is the machine
C             precision (see LAPACK Library routine DLAMCH),
C             I = 1,2,...,M and J = 1,2,...,N.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension MAX(M,N)
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
C     Given an M-by-N matrix pencil sE-A with E not necessarily regular,
C     the routine computes a unitary transformed pencil Q'(sE-A)Z such
C     that the matrix Q' * E * Z is in column echelon form (trapezoidal
C     form).  Further details can be found in [1].
C
C     [An M-by-N matrix E with rank(E) = r is said to be in column
C     echelon form if the following conditions are satisfied:
C     (a) the first (N - r) columns contain only zero elements; and
C     (b) if E(i(k),k) is the last nonzero element in column k for
C         k = N-r+1,...,N, i.e. E(i(k),k) <> 0 and E(j,k) = 0 for
C         j > i(k), then 1 <= i(N-r+1) < i(N-r+2) < ... < i(N) <= M.]
C
C     REFERENCES
C
C     [1] Beelen, Th. and Van Dooren, P.
C         An improved algorithm for the computation of Kronecker's
C         canonical form of a singular pencil.
C         Linear Algebra and Applications, 105, pp. 9-65, 1988.
C
C     NUMERICAL ASPECTS
C
C     It is shown in [1] that the algorithm is numerically backward
C     stable. The operations count is proportional to (MAX(M,N))**3.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998.
C     Based on Release 3.0 routine MB04SD modified by A. Varga,
C     German Aerospace Research Establishment, Oberpfaffenhofen,
C     Germany, Dec. 1997, to transform also the matrix A.
C
C     REVISIONS
C
C     A. Varga, DLR Oberpfaffenhofen, June 2005.
C
C     KEYWORDS
C
C     Echelon form, orthogonal transformation, staircase form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         JOBQ, JOBZ
      INTEGER           INFO, LDA, LDE, LDQ, LDZ, M, N, RANKE
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           ISTAIR(*)
      DOUBLE PRECISION  A(LDA,*), DWORK(*), E(LDE,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           LJOBQI, LJOBZI, LZERO, UPDATQ, UPDATZ
      INTEGER           I, K, KM1, L, LK, MNK, NR1
      DOUBLE PRECISION  EMX, EMXNRM, TAU, TOLER
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           IDAMAX
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL          DLARF, DLARFG, DLASET, DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LJOBQI = LSAME( JOBQ, 'I' )
      UPDATQ = LJOBQI.OR.LSAME( JOBQ, 'U' )
      LJOBZI = LSAME( JOBZ, 'I' )
      UPDATZ = LJOBZI.OR.LSAME( JOBZ, 'U' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.UPDATQ .AND. .NOT.LSAME( JOBQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPDATZ .AND. .NOT.LSAME( JOBZ, 'N' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDE.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( .NOT.UPDATQ .AND. LDQ.LT.1 .OR.
     $              UPDATQ .AND. LDQ.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( .NOT.UPDATZ .AND. LDZ.LT.1 .OR.
     $              UPDATZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB04UD', -INFO )
         RETURN
      END IF
C
C     Initialize Q and Z to the identity matrices, if needed.
C
      IF ( LJOBQI )
     $   CALL DLASET( 'Full', M, M, ZERO, ONE, Q, LDQ )
      IF ( LJOBZI )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Quick return if possible.
C
      RANKE = MIN( M, N )
C
      IF ( RANKE.EQ.0 )
     $   RETURN
C
      TOLER = TOL
      IF ( TOLER.LE.ZERO )
     $   TOLER = DLAMCH( 'Epsilon' )*DLANGE( 'M', M, N, E, LDE, DWORK )
C
      K = N
      LZERO = .FALSE.
C
C     WHILE ( ( K > 0 ) AND ( NOT a zero submatrix encountered ) ) DO
   20 IF ( ( K.GT.0 ) .AND. ( .NOT. LZERO ) ) THEN
C
C         Intermediate form of E
C
C                     <--k--><--n-k->
C                l=1 |x....x|       |
C                    |      |       |
C                    |  Ek  |   X   |
C                    |      |       |
C            l=m-n+k |x....x|       |
C                    ----------------
C                    |      |x ... x|  }
C                    |  O   |  x x x|  }
C                    |      |    x x|  } n-k
C                    |      |      x|  }
C
C        where submatrix Ek = E[1:m-n+k;1:k].
C
C        Determine row LK in submatrix Ek with largest max-norm
C        (starting with row m-n+k).
C
         MNK = M - N + K
         EMXNRM = ZERO
         LK = MNK
C
         DO 40 L = MNK, 1, -1
            EMX = ABS( E(L,IDAMAX( K, E(L,1), LDE )) )
            IF ( EMX.GT.EMXNRM ) THEN
               EMXNRM = EMX
               LK = L
            END IF
   40    CONTINUE
C
         IF ( EMXNRM.LE.TOLER ) THEN
C
C           Set submatrix Ek to zero.
C
            CALL DLASET( 'Full', MNK, K, ZERO, ZERO, E, LDE )
            LZERO = .TRUE.
            RANKE = N - K
         ELSE
C
C           Submatrix Ek is not considered to be identically zero.
C           Check whether rows have to be interchanged.
C
            IF ( LK.NE.MNK ) THEN
C
C              Interchange rows lk and m-n+k in whole A- and E-matrix
C              and update the row transformation matrix Q, if needed.
C              (For Q, the number of elements involved is m.)
C
               CALL DSWAP( N, E(LK,1), LDE, E(MNK,1), LDE )
               CALL DSWAP( N, A(LK,1), LDA, A(MNK,1), LDA )
               IF( UPDATQ ) CALL DSWAP( M, Q(1,LK), 1, Q(1,MNK), 1 )
            END IF
C
            KM1 = K - 1
C
C           Determine a Householder transformation to annihilate
C           E(m-n+k,1:k-1) using E(m-n+k,k) as pivot.
C           Apply the transformation to the columns of A and Ek
C           (number of elements involved is m for A and m-n+k for Ek).
C           Update the column transformation matrix Z, if needed
C           (number of elements involved is n).
C
            CALL DLARFG( K, E(MNK,K), E(MNK,1), LDE, TAU )
            EMX = E(MNK,K)
            E(MNK,K) = ONE
            CALL DLARF( 'Right', MNK-1, K, E(MNK,1), LDE, TAU, E, LDE,
     $                   DWORK )
            CALL DLARF( 'Right', M, K, E(MNK,1), LDE, TAU, A, LDA,
     $                   DWORK )
            IF( UPDATZ ) CALL DLARF( 'Right', N, K, E(MNK,1), LDE, TAU,
     $                               Z, LDZ, DWORK )
            E(MNK,K) = EMX
            CALL DLASET( 'Full', 1, KM1, ZERO, ZERO, E(MNK,1), LDE )
C
            K = KM1
         END IF
         GO TO 20
      END IF
C     END WHILE 20
C
C     Initialise administration staircase form, i.e.
C     ISTAIR(i) =  j  if E(i,j) is a nonzero corner point
C               = -j  if E(i,j) is on the boundary but is no corner
C                     point.
C     Thus,
C     ISTAIR(m-k) =   n-k           for k=0,...,rank(E)-1
C                 = -(n-rank(E)+1)  for k=rank(E),...,m-1.
C
      DO 60 I = 0, RANKE - 1
         ISTAIR(M-I) = N - I
   60 CONTINUE
C
      NR1 = -(N - RANKE + 1)
C
      DO 80 I = 1, M - RANKE
         ISTAIR(I) = NR1
   80 CONTINUE
C
      RETURN
C *** Last line of MB04UD ***
      END
