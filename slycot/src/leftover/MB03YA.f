      SUBROUTINE MB03YA( WANTT, WANTQ, WANTZ, N, ILO, IHI, ILOQ, IHIQ,
     $                   POS, A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO )
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
C     To annihilate one or two entries on the subdiagonal of the
C     Hessenberg matrix A for dealing with zero elements on the diagonal
C     of the triangular matrix B.
C
C     MB03YA is an auxiliary routine called by SLICOT Library routines
C     MB03XP and MB03YD.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     WANTT   LOGICAL
C             Indicates whether the user wishes to compute the full
C             Schur form or the eigenvalues only, as follows:
C             = .TRUE. :  Compute the full Schur form;
C             = .FALSE.:  compute the eigenvalues only.
C
C     WANTQ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Q as follows:
C             = .TRUE. :  The matrix Q is updated;
C             = .FALSE.:  the matrix Q is not required.
C
C     WANTZ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Z as follows:
C             = .TRUE. :  The matrix Z is updated;
C             = .FALSE.:  the matrix Z is not required.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and B. N >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that the matrices A and B are already
C             (quasi) upper triangular in rows and columns 1:ILO-1 and
C             IHI+1:N. The routine works primarily with the submatrices
C             in rows and columns ILO to IHI, but applies the
C             transformations to all the rows and columns of the
C             matrices A and B, if WANTT = .TRUE..
C             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N.
C
C     ILOQ    (input) INTEGER
C     IHIQ    (input) INTEGER
C             Specify the rows of Q and Z to which transformations
C             must be applied if WANTQ = .TRUE. and WANTZ = .TRUE.,
C             respectively.
C             1 <= ILOQ <= ILO; IHI <= IHIQ <= N.
C
C     POS     (input) INTEGER
C             The position of the zero element on the diagonal of B.
C             ILO <= POS <= IHI.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper Hessenberg matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the updated matrix A where A(POS,POS-1) = 0, if POS > ILO,
C             and A(POS+1,POS) = 0, if POS < IHI.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain an upper triangular matrix B with B(POS,POS) = 0.
C             On exit, the leading N-by-N part of this array contains
C             the updated upper triangular matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, if WANTQ = .TRUE., then the leading N-by-N part
C             of this array must contain the current matrix Q of
C             transformations accumulated by MB03XP.
C             On exit, if WANTQ = .TRUE., then the leading N-by-N part
C             of this array contains the matrix Q updated in the
C             submatrix Q(ILOQ:IHIQ,ILO:IHI).
C             If WANTQ = .FALSE., Q is not referenced.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= 1.
C             If WANTQ = .TRUE., LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, if WANTZ = .TRUE., then the leading N-by-N part
C             of this array must contain the current matrix Z of
C             transformations accumulated by MB03XP.
C             On exit, if WANTZ = .TRUE., then the leading N-by-N part
C             of this array contains the matrix Z updated in the
C             submatrix Z(ILOQ:IHIQ,ILO:IHI).
C             If WANTZ = .FALSE., Z is not referenced.
C
C     LDZ     INTEGER
C             The leading dimension of the array Z.  LDZ >= 1.
C             If WANTZ = .TRUE., LDZ >= MAX(1,N).
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
C     The method is illustrated by Wilkinson diagrams for N = 5,
C     POS = 3:
C
C           [ x x x x x ]       [ x x x x x ]
C           [ x x x x x ]       [ o x x x x ]
C       A = [ o x x x x ],  B = [ o o o x x ].
C           [ o o x x x ]       [ o o o x x ]
C           [ o o o x x ]       [ o o o o x ]
C
C     First, a QR factorization is applied to A(1:3,1:3) and the
C     resulting nonzero in the updated matrix B is immediately
C     annihilated by a Givens rotation acting on columns 1 and 2:
C
C           [ x x x x x ]       [ x x x x x ]
C           [ x x x x x ]       [ o x x x x ]
C       A = [ o o x x x ],  B = [ o o o x x ].
C           [ o o x x x ]       [ o o o x x ]
C           [ o o o x x ]       [ o o o o x ]
C
C     Secondly, an RQ factorization is applied to A(4:5,4:5) and the
C     resulting nonzero in the updated matrix B is immediately
C     annihilated by a Givens rotation acting on rows 4 and 5:
C
C           [ x x x x x ]       [ x x x x x ]
C           [ x x x x x ]       [ o x x x x ]
C       A = [ o o x x x ],  B = [ o o o x x ].
C           [ o o o x x ]       [ o o o x x ]
C           [ o o o x x ]       [ o o o o x ]
C
C     REFERENCES
C
C     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P.
C         The periodic Schur decomposition: Algorithms and applications.
C         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42,
C         1992.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(N**2) floating point operations and is
C     backward stable.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLADFB).
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTT, WANTZ
      INTEGER            IHI, IHIQ, ILO, ILOQ, INFO, LDA, LDB, LDQ, LDZ,
     $                   N, POS
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      INTEGER            I1, I2, J, NQ
      DOUBLE PRECISION   CS, SN, TEMP
C     .. External Subroutines ..
      EXTERNAL           DLARTG, DROT, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      NQ = IHIQ - ILOQ + 1
      IF ( N.LT.0 ) THEN
         INFO = -4
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF ( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -6
      ELSE IF ( ILOQ.LT.1 .OR. ILOQ.GT.ILO ) THEN
         INFO = -7
      ELSE IF ( IHIQ.LT.IHI .OR. IHIQ.GT.N ) THEN
         INFO = -8
      ELSE IF ( POS.LT.ILO .OR. POS.GT.IHI ) THEN
         INFO = -9
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF ( LDQ.LT.1 .OR. WANTQ .AND. LDQ.LT.N ) THEN
         INFO = -15
      ELSE IF ( LDZ.LT.1 .OR. WANTZ .AND. LDZ.LT.N ) THEN
         INFO = -17
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB03YA', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 )
     $   RETURN
C
      IF ( WANTT ) THEN
         I1 = 1
         I2 = N
      ELSE
         I1 = ILO
         I2 = IHI
      END IF
C
C     Apply a zero-shifted QR step.
C
      DO 10  J = ILO, POS-1
         TEMP = A(J,J)
         CALL DLARTG( TEMP, A(J+1,J), CS, SN, A(J,J) )
         A(J+1,J) = ZERO
         CALL DROT( I2-J, A(J,J+1), LDA, A(J+1,J+1), LDA, CS, SN )
         CALL DROT( MIN(J,POS-2)-I1+2, B(I1,J), 1, B(I1,J+1), 1, CS,
     $              SN )
         IF ( WANTQ )
     $      CALL DROT( NQ, Q(ILOQ,J), 1, Q(ILOQ,J+1), 1, CS, SN )
   10 CONTINUE
      DO 20  J = ILO, POS-2
         TEMP = B(J,J)
         CALL DLARTG( TEMP, B(J+1,J), CS, SN, B(J,J) )
         B(J+1,J) = ZERO
         CALL DROT( I2-J, B(J,J+1), LDB, B(J+1,J+1), LDB, CS, SN )
         CALL DROT( J-I1+2, A(I1,J), 1, A(I1,J+1), 1, CS, SN )
         IF ( WANTZ )
     $      CALL DROT( NQ, Z(ILOQ,J), 1, Z(ILOQ,J+1), 1, CS, SN )
   20 CONTINUE
C
C     Apply a zero-shifted RQ step.
C
      DO 30  J = IHI, POS+1, -1
         TEMP = A(J,J)
         CALL DLARTG( TEMP, A(J,J-1), CS, SN, A(J,J) )
         A(J,J-1) = ZERO
         SN = -SN
         CALL DROT( J-I1, A(I1,J-1), 1, A(I1,J), 1, CS, SN )
         CALL DROT( I2 - MAX( J-1,POS+1 ) + 1, B(J-1,MAX( J-1,POS+1 )),
     $              LDB, B(J,MAX(J-1,POS+1)), LDB, CS, SN )
         IF ( WANTZ )
     $      CALL DROT( NQ, Z(ILOQ,J-1), 1, Z(ILOQ,J), 1, CS, SN )
   30 CONTINUE
      DO 40  J = IHI, POS+2, -1
         TEMP = B(J,J)
         CALL DLARTG( TEMP, B(J,J-1), CS, SN, B(J,J) )
         B(J,J-1) = ZERO
         SN = -SN
         CALL DROT( J-I1, B(I1,J-1), 1, B(I1,J), 1, CS, SN )
         CALL DROT( I2-J+2, A(J-1,J-1), LDA, A(J,J-1), LDA, CS, SN )
         IF ( WANTQ )
     $      CALL DROT( NQ, Q(ILOQ,J-1), 1, Q(ILOQ,J), 1, CS, SN )
   40 CONTINUE
      RETURN
C *** Last line of MB03YA ***
      END
