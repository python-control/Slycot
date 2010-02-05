      SUBROUTINE AB09DD( DICO, N, M, P, NR, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, RCOND, IWORK, DWORK, INFO )
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
C     To compute a reduced order model by using singular perturbation
C     approximation formulas.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the state vector, i.e. the order of the
C             matrix A; also the number of rows of matrix B and the
C             number of columns of the matrix C.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector, i.e. the number of columns
C             of matrices B and D.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of output vector, i.e. the number of rows of
C             matrices C and D.  P >= 0.
C
C     NR      (input) INTEGER
C             The order of the reduced order system.  N >= NR >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix of the original system.
C             On exit, the leading NR-by-NR part of this array contains
C             the state dynamics matrix Ar of the reduced order system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix of the original system.
C             On exit, the leading NR-by-M part of this array contains
C             the input/state matrix Br of the reduced order system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix of the original system.
C             On exit, the leading P-by-NR part of this array contains
C             the state/output matrix Cr of the reduced order system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the input/output matrix of the original system.
C             On exit, the leading P-by-M part of this array contains
C             the input/output matrix Dr of the reduced order system.
C             If NR = 0 and the given system is stable, then D contains
C             the steady state gain of the system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal condition number of the matrix A22-g*I
C             (see METHOD).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension 2*(N-NR)
C
C     DWORK   DOUBLE PRECISION array, dimension 4*(N-NR)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1: if the matrix A22-g*I (see METHOD) is numerically
C                  singular.
C
C     METHOD
C
C     Given the system (A,B,C,D), partition the system matrices as
C
C            ( A11 A12 )        ( B1 )
C        A = (         ) ,  B = (    ) ,  C = ( C1  C2 ),
C            ( A21 A22 )        ( B2 )
C
C     where A11 is NR-by-NR, B1 is NR-by-M, C1 is P-by-NR, and the other
C     submatrices have appropriate dimensions.
C
C     The matrices of the reduced order system (Ar,Br,Cr,Dr) are
C     computed according to the following residualization formulas:
C                                -1                               -1
C        Ar = A11 + A12*(g*I-A22)  *A21 ,  Br = B1 + A12*(g*I-A22)  *B2
C                              -1                               -1
C        Cr = C1 + C2*(g*I-A22)  *A21   ,  Dr = D + C2*(g*I-A22)  *B2
C
C     where g = 0 if DICO = 'C' and g = 1 if DICO = 'D'.
C
C     CONTRIBUTOR
C
C     C. Oara and A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, March 1998.
C     Based on the RASP routine SRESID.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Model reduction, multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO
      INTEGER           INFO, LDA, LDB, LDC, LDD, M, N, NR, P
      DOUBLE PRECISION  RCOND
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), DWORK(*)
      INTEGER           IWORK(*)
C     .. Local Scalars
      LOGICAL           DISCR
      INTEGER           I, J, K, NS
      DOUBLE PRECISION  A22NRM
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGECON, DGEMM, DGETRF, DGETRS, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Check the input scalar arguments.
C
      INFO = 0
      DISCR = LSAME( DICO, 'D' )
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( NR.LT.0 .OR. NR.GT.N ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -13
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09DD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( NR.EQ.N ) THEN
         RCOND = ONE
         RETURN
      END IF
C
      K  = NR + 1
      NS = N - NR
C
C     Compute: T = -A22   if  DICO = 'C' and
C              T = -A22+I if  DICO = 'D'.
C
      DO 20 J = K, N
         DO 10 I = K, N
            A(I,J) = -A(I,J)
   10    CONTINUE
         IF( DISCR ) A(J,J) = A(J,J) + ONE
   20 CONTINUE
C
C     Compute the LU decomposition of T.
C
      A22NRM = DLANGE( '1-norm', NS, NS, A(K,K), LDA, DWORK )
      CALL DGETRF( NS, NS, A(K,K), LDA, IWORK, INFO )
      IF( INFO.GT.0 ) THEN
C
C        Error return.
C
         RCOND = ZERO
         INFO = 1
         RETURN
      END IF
      CALL DGECON( '1-norm', NS, A(K,K), LDA, A22NRM, RCOND, DWORK,
     $             IWORK(NS+1), INFO )
      IF( RCOND.LE.DLAMCH('E') ) THEN
C
C        Error return.
C
         INFO = 1
         RETURN
      END IF
C
C     Compute A21 <- INV(T)*A21.
C
      CALL DGETRS( 'NoTranspose', NS, NR, A(K,K), LDA, IWORK, A(K,1),
     $             LDA, INFO )
C
C     Compute B2 <- INV(T)*B2.
C
      CALL DGETRS( 'NoTranspose', NS, M, A(K,K), LDA, IWORK, B(K,1),
     $             LDB, INFO )
C
C     Compute the residualized systems matrices.
C     Ar = A11 + A12*INV(T)*A21.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, NR, NS, ONE, A(1,K),
     $            LDA, A(K,1), LDA, ONE, A, LDA )
C
C     Br = B1 + A12*INV(T)*B2.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', NR, M, NS, ONE, A(1,K),
     $            LDA, B(K,1), LDB, ONE, B, LDB )
C
C     Cr = C1 + C2*INV(T)*A21.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, NR, NS, ONE, C(1,K),
     $            LDC, A(K,1), LDA, ONE, C, LDC )
C
C     Dr = D + C2*INV(T)*B2.
C
      CALL DGEMM( 'NoTranspose', 'NoTranspose', P, M, NS, ONE, C(1,K),
     $            LDC, B(K,1), LDB, ONE, D, LDD )
C
      RETURN
C *** Last line of AB09DD ***
      END
