      SUBROUTINE MB04TT( UPDATQ, UPDATZ, M, N, IFIRA, IFICA, NCA, A,
     $                   LDA, E, LDE, Q, LDQ, Z, LDZ, ISTAIR, RANK, TOL,
     $                   IWORK )
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
C     Let A and E be M-by-N matrices with E in column echelon form.
C     Let AA and EE be the following submatrices of A and E:
C       AA := A(IFIRA : M ; IFICA : N)
C       EE := E(IFIRA : M ; IFICA : N).
C     Let Aj and Ej be the following submatrices of AA and EE:
C       Aj := A(IFIRA : M ; IFICA : IFICA + NCA - 1) and
C       Ej := E(IFIRA : M ; IFICA + NCA : N).
C
C     To transform (AA,EE) such that Aj is row compressed while keeping
C     matrix Ej in column echelon form (which may be different from the
C     form on entry).
C     In fact the routine performs the j-th step of Algorithm 3.2.1 in
C     [1]. Furthermore, it determines the rank RANK of the submatrix Ej,
C     which is equal to the number of corner points in submatrix Ej.
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
C             M is the number of rows of the matrices A, E and Q.
C             M >= 0.
C
C     N       (input) INTEGER
C             N is the number of columns of the matrices A, E and Z.
C             N >= 0.
C
C     IFIRA   (input) INTEGER
C             IFIRA is the first row index of the submatrices Aj and Ej
C             in the matrices A and E, respectively.
C
C     IFICA   (input) INTEGER
C             IFICA and IFICA + NCA are the first column indices of the
C             submatrices Aj and Ej in the matrices A and E,
C             respectively.
C
C     NCA     (input) INTEGER
C             NCA is the number of columns of the submatrix Aj in A.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, A(IFIRA : M ; IFICA : IFICA + NCA - 1) contains
C             the matrix Aj.
C             On exit, it contains the matrix A with AA that has been
C             row compressed while keeping EE in column echelon form.
C
C     LDA     INTEGER
C             The leading dimension of array A. LDA >= MAX(1,M).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, E(IFIRA : M ; IFICA + NCA : N) contains the
C             matrix Ej which is in column echelon form.
C             On exit, it contains the transformed matrix EE which is
C             kept in column echelon form.
C
C     LDE     INTEGER
C             The leading dimension of array E. LDE >= MAX(1,M).
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
C     ISTAIR  (input/output) INTEGER array, dimension (M)
C             On entry, ISTAIR contains information on the column
C             echelon form of the input matrix E as follows:
C             ISTAIR(i) = +j: the boundary element E(i,j) is a corner
C                             point;
C                         -j: the boundary element E(i,j) is not a
C                             corner point (where i=1,...,M).
C             On exit, ISTAIR contains the same information for the
C             transformed matrix E.
C
C     RANK    (output) INTEGER
C             Numerical rank of the submatrix Aj in A (based on TOL).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance used when considering matrix elements
C             to be zero.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
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
C                               3
C     The algorithm requires 0(N ) operations and is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MB04FZ by Th.G.J. Beelen,
C     Philips Glass Eindhoven, Holland.
C
C     REVISIONS
C
C     June 13, 1997, V. Sima.
C     November 24, 1997, A. Varga: array starting point A(KK,LL)
C                                  correctly set when calling DLASET.
C
C     KEYWORDS
C
C     Echelon form, orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           UPDATQ, UPDATZ
      INTEGER           IFICA, IFIRA, LDA, LDE, LDQ, LDZ, M, N, NCA,
     $                  RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           ISTAIR(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           LZERO
      INTEGER           I, IFICA1, IFIRA1, II, IP, IST1, IST2, ISTPVT,
     $                  ITYPE, JC1, JC2, JPVT, K, KK, L, LL, LSAV, MJ,
     $                  MK1, MXRANK, NJ
      DOUBLE PRECISION  BMX, BMXNRM, EIJPVT, SC, SS
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DLAPMT, DLASET, DROT, DROTG, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
C
      RANK = 0
      IF ( M.LE.0 .OR. N.LE.0 )
     $   RETURN
C
C     Initialisation.
C
C     NJ = number of columns in submatrix Aj,
C     MJ = number of rows in submatrices Aj and Ej.
C
      NJ = NCA
      MJ = M + 1 - IFIRA
      IFIRA1 = IFIRA - 1
      IFICA1 = IFICA - 1
C
      DO 20 I = 1, NJ
         IWORK(I) = I
   20 CONTINUE
C
      K = 1
      LZERO = .FALSE.
      RANK = MIN( NJ, MJ )
      MXRANK = RANK
C
C     WHILE ( K <= MXRANK ) and ( LZERO = FALSE ) DO
   40 IF ( ( K.LE.MXRANK ) .AND. ( .NOT.LZERO ) ) THEN
C
C        Determine column in Aj with largest max-norm.
C
         BMXNRM = ZERO
         LSAV = K
         KK = IFIRA1 + K
C
         DO 60 L = K, NJ
C
C           IDAMAX call gives the relative index in column L of Aj where
C           max element is found.
C           Note: the first element in column L is in row K of
C                 matrix Aj.
C
            LL = IFICA1 + L
            BMX = ABS( A(IDAMAX( MJ-K+1, A(KK,LL), 1 )+KK-1,LL) )
            IF ( BMX.GT.BMXNRM ) THEN
               BMXNRM = BMX
               LSAV = L
            END IF
   60    CONTINUE
C
         LL = IFICA1 + K
         IF ( BMXNRM.LT.TOL ) THEN
C
C           Set submatrix of Aj to zero.
C
            CALL DLASET( 'Full', MJ-K+1, NJ-K+1, ZERO, ZERO, A(KK,LL),
     $                   LDA )
            LZERO = .TRUE.
            RANK = K - 1
         ELSE
C
C           Check whether columns have to be interchanged.
C
            IF ( LSAV.NE.K ) THEN
C
C              Interchange the columns in A which correspond to the
C              columns lsav and k in Aj. Store the permutation in IWORK.
C
               CALL DSWAP( M, A(1,LL), 1, A(1,IFICA1+LSAV), 1 )
               IP = IWORK(LSAV)
               IWORK(LSAV) = IWORK(K)
               IWORK(K) = IP
            END IF
C
            K = K + 1
            MK1 = N - LL + 1
C
            DO 80 I = MJ, K, -1
C
C              II = absolute row number in A corresponding to row i in
C                   Aj.
C
               II = IFIRA1 + I
C
C              Construct Givens transformation to annihilate Aj(i,k).
C              Apply the row transformation to whole matrix A
C              (NOT only to Aj).
C              Update row transformation matrix Q, if needed.
C
               CALL DROTG( A(II-1,LL), A(II,LL), SC, SS )
               CALL DROT( MK1-1, A(II-1,LL+1), LDA, A(II,LL+1), LDA, SC,
     $                    SS )
               A(II,LL) = ZERO
               IF ( UPDATQ )
     $            CALL DROT( M, Q(1,II-1), 1, Q(1,II), 1, SC, SS )
C
C              Determine boundary type of matrix E at rows II-1 and II.
C
               IST1 = ISTAIR(II-1)
               IST2 = ISTAIR(II)
               IF ( ( IST1*IST2 ).GT.0 ) THEN
                  IF ( IST1.GT.0 ) THEN
C
C                    boundary form = (* x)
C                                    (0 *)
C
                     ITYPE = 1
                  ELSE
C
C                    boundary form = (x x)
C                                    (x x)
C
                     ITYPE = 3
                  END IF
               ELSE
                  IF ( IST1.LT.0 ) THEN
C
C                    boundary form = (x x)
C                                    (* x)
C
                     ITYPE = 2
                  ELSE
C
C                    boundary form = (* x)
C                                    (0 x)
C
                     ITYPE = 4
                  END IF
               END IF
C
C              Apply row transformation also to matrix E.
C
C              JC1 = absolute number of the column in E in which stair
C                    element of row i-1 of Ej is present.
C              JC2 = absolute number of the column in E in which stair
C                    element of row i of Ej is present.
C
C              Note: JC1 < JC2   if ITYPE = 1.
C                    JC1 = JC2   if ITYPE = 2, 3 or 4.
C
               JC1 = ABS( IST1 )
               JC2 = ABS( IST2 )
               JPVT = MIN( JC1, JC2 )
C
               CALL DROT( N-JPVT+1, E(II-1,JPVT), LDE, E(II,JPVT), LDE,
     $                    SC, SS )
               EIJPVT = E(II,JPVT)
C
               IF ( ITYPE.EQ.1 ) THEN
C
C                 Construct column Givens transformation to annihilate
C                 E(ii,jpvt).
C                 Apply column Givens transformation to matrix E
C                 (NOT only to Ej).
C
                  CALL DROTG( E(II,JPVT+1), E(II,JPVT), SC, SS )
                  CALL DROT( II-1, E(1,JPVT+1), 1, E(1,JPVT), 1, SC,
     $                       SS )
                  E(II,JPVT) = ZERO
C
C                 Apply this transformation also to matrix A
C                 (NOT only to Aj).
C                 Update column transformation matrix Z, if needed.
C
                  CALL DROT( M, A(1,JPVT+1), 1, A(1,JPVT), 1, SC, SS )
                  IF ( UPDATZ ) CALL DROT( N, Z(1,JPVT+1), 1, Z(1,JPVT),
     $                                     1, SC, SS )
C
               ELSE IF ( ITYPE.EQ.2 ) THEN
                  IF ( ABS( EIJPVT ).LT.TOL ) THEN
C
C                                                        (x x)    (* x)
C                    Boundary form has been changed from (* x) to (0 x).
C
                     ISTPVT = ISTAIR(II)
                     ISTAIR(II-1) = ISTPVT
                     ISTAIR(II) = -(ISTPVT+1 )
                     E(II,JPVT) = ZERO
                  END IF
C
               ELSE IF ( ITYPE.EQ.4 ) THEN
                  IF ( ABS( EIJPVT ).GE.TOL ) THEN
C
C                                                        (* x)    (x x)
C                    Boundary form has been changed from (0 x) to (* x).
C
                     ISTPVT = ISTAIR(II-1)
                     ISTAIR(II-1) = -ISTPVT
                     ISTAIR(II) = ISTPVT
                  END IF
               END IF
   80       CONTINUE
C
         END IF
         GO TO 40
      END IF
C     END WHILE 40
C
C     Permute columns of Aj to original order.
C
      CALL DLAPMT( .FALSE., IFIRA1+RANK, NJ, A(1,IFICA), LDA, IWORK )
C
      RETURN
C *** Last line of MB04TT ***
      END
