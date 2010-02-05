      SUBROUTINE MB04TB( TRANA, TRANB, N, ILO, A, LDA, B, LDB, G, LDG,
     $                   Q, LDQ, CSL, CSR, TAUL, TAUR, DWORK, LDWORK,
     $                   INFO )
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
C     To compute a symplectic URV (SURV) decomposition of a real
C     2N-by-2N matrix H,
C
C            [ op(A)   G   ]                 [ op(R11)   R12   ]
C        H = [             ] = U R V'  = U * [                 ] * V' ,
C            [  Q    op(B) ]                 [   0     op(R22) ]
C
C     where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real
C     N-by-N upper triangular matrix, op(R22) is a real N-by-N lower
C     Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic
C     matrices. Blocked version.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANA   CHARACTER*1
C             Specifies the form of op( A ) as follows:
C             = 'N':  op( A ) = A;
C             = 'T':  op( A ) = A';
C             = 'C':  op( A ) = A'.
C
C     TRANB   CHARACTER*1
C             Specifies the form of op( B ) as follows:
C             = 'N':  op( B ) = B;
C             = 'T':  op( B ) = B';
C             = 'C':  op( B ) = B'.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A. N >= 0.
C
C     ILO     (input) INTEGER
C             It is assumed that op(A) is already upper triangular,
C             op(B) is lower triangular and Q is zero in rows and
C             columns 1:ILO-1. ILO is normally set by a previous call
C             to MB04DD; otherwise it should be set to 1.
C             1 <= ILO <= N, if N > 0; ILO=1, if N=0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the triangular matrix R11, and in the zero part
C             information about the elementary reflectors used to
C             compute the SURV decomposition.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix B.
C             On exit, the leading N-by-N part of this array contains
C             the Hessenberg matrix R22, and in the zero part
C             information about the elementary reflectors used to
C             compute the SURV decomposition.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix G.
C             On exit, the leading N-by-N part of this array contains
C             the matrix R12.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= MAX(1,N).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix Q.
C             On exit, the leading N-by-N part of this array contains
C             information about the elementary reflectors used to
C             compute the SURV decomposition.
C
C     LDQ     INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1,N).
C
C     CSL     (output) DOUBLE PRECISION array, dimension (2N)
C             On exit, the first 2N elements of this array contain the
C             cosines and sines of the symplectic Givens rotations
C             applied from the left-hand side used to compute the SURV
C             decomposition.
C
C     CSR     (output) DOUBLE PRECISION array, dimension (2N-2)
C             On exit, the first 2N-2 elements of this array contain the
C             cosines and sines of the symplectic Givens rotations
C             applied from the right-hand side used to compute the SURV
C             decomposition.
C
C     TAUL    (output) DOUBLE PRECISION array, dimension (N)
C             On exit, the first N elements of this array contain the
C             scalar factors of some of the elementary reflectors
C             applied form the left-hand side.
C
C     TAUR    (output) DOUBLE PRECISION array, dimension (N-1)
C             On exit, the first N-1 elements of this array contain the
C             scalar factors of some of the elementary reflectors
C             applied form the right-hand side.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK, (16*N + 5)*NB, where NB is the optimal
C             block size determined by the function UE01MD.
C             On exit, if  INFO = -16,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
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
C     The matrices U and V are represented as products of symplectic
C     reflectors and Givens rotators
C
C     U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) )
C         diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) )
C                              ....
C         diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ),
C
C     V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) )
C         diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) )
C                                   ....
C         diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ).
C
C     Each HU(i) has the form
C
C           HU(i) = I - tau * v * v'
C
C     where tau is a real scalar, and v is a real vector with
C     v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in
C     Q(i+1:n,i), and tau in Q(i,i).
C
C     Each FU(i) has the form
C
C           FU(i) = I - nu * w * w'
C
C     where nu is a real scalar, and w is a real vector with
C     w(1:i-1) = 0 and w(i) = 1; w(i+1:n) is stored on exit in
C     A(i+1:n,i), if op(A) = 'N', and in A(i,i+1:n), otherwise. The
C     scalar nu is stored in TAUL(i).
C
C     Each GU(i) is a Givens rotator acting on rows i and n+i,
C     where the cosine is stored in CSL(2*i-1) and the sine in
C     CSL(2*i).
C
C     Each HV(i) has the form
C
C           HV(i) = I - tau * v * v'
C
C     where tau is a real scalar, and v is a real vector with
C     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
C     Q(i,i+2:n), and tau in Q(i,i+1).
C
C     Each FV(i) has the form
C
C           FV(i) = I - nu * w * w'
C
C     where nu is a real scalar, and w is a real vector with
C     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in
C     B(i,i+2:n), if op(B) = 'N', and in B(i+2:n,i), otherwise.
C     The scalar nu is stored in TAUR(i).
C
C     Each GV(i) is a Givens rotator acting on columns i+1 and n+i+1,
C     where the cosine is stored in CSR(2*i-1) and the sine in
C     CSR(2*i).
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 80/3*N**3 + ( 64*NB + 77 )*N**2 +
C     ( -16*NB + 48 )*NB*N + O(N) floating point operations, where
C     NB is the used block size, and is numerically backward stable.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A numerically stable, structure preserving method for
C         computing the eigenvalues of real Hamiltonian or symplectic
C         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998.
C
C     [2] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DGESUB).
C
C     KEYWORDS
C
C     Elementary matrix operations, Matrix decompositions, Hamiltonian
C     matrix
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TRANA, TRANB
      INTEGER           ILO, INFO, LDA, LDB, LDG, LDQ, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), CSL(*), CSR(*), DWORK(*),
     $                  G(LDG,*), Q(LDQ,*), TAUL(*), TAUR(*)
C     .. Local Scalars ..
      LOGICAL           LTRA, LTRB
      INTEGER           I, IB, IERR, NB, NBMIN, NH, NIB, NNB, NX, PDW,
     $                  PXA, PXB, PXG, PXQ, PYA, PYB, PYG, PYQ, WRKOPT
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           UE01MD
      EXTERNAL          LSAME, UE01MD
C     .. External Subroutines ..
      EXTERNAL          DGEMM, MB03XU, MB04TS, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      LTRA = LSAME( TRANA, 'T' ) .OR. LSAME( TRANA, 'C' )
      LTRB = LSAME( TRANB, 'T' ) .OR. LSAME( TRANB, 'C' )
      IF ( .NOT.LTRA .AND. .NOT.LSAME( TRANA, 'N' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.LTRB .AND. .NOT.LSAME( TRANB, 'N' ) ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF ( LDG.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF ( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = DBLE( MAX( 1, N ) )
         INFO = -18
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04TB', -INFO )
         RETURN
      END IF
C
C     Set elements 1:ILO-1 of CSL, CSR, TAUL and TAUR to their default
C     values.
C
      DO 10  I = 1, ILO - 1
         CSL(2*I-1) = ONE
         CSL(2*I) = ZERO
         CSR(2*I-1) = ONE
         CSR(2*I) = ZERO
         TAUL(I) = ZERO
         TAUR(I) = ZERO
   10 CONTINUE
C
C     Quick return if possible.
C
      NH = N - ILO + 1
      IF ( NH.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Determine the block size.
C
      NB = UE01MD( 1, 'MB04TB', TRANA // TRANB, N, ILO, -1 )
      NBMIN = 2
      WRKOPT = N
      IF ( NB.GT.1 .AND. NB.LT.NH ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( NB, UE01MD( 3, 'MB04TB', TRANA // TRANB, N, ILO, -1 )
     $           )
         IF ( NX.LT.NH ) THEN
C
C           Check whether workspace is large enough for blocked code.
C
            WRKOPT = 16*N*NB + 5*NB
            IF ( LDWORK.LT.WRKOPT ) THEN
C
C              Not enough workspace available. Determine minimum value
C              of NB, and reduce NB.
C
               NBMIN = MAX( 2, UE01MD( 2, 'MB04TB', TRANA // TRANB, N,
     $                                 ILO, -1 ) )
               NB = LDWORK / ( 16*N + 5 )
            END IF
         END IF
      END IF
C
      NNB = N*NB
      PYB = 1
      PYQ = PYB + 2*NNB
      PYA = PYQ + 2*NNB
      PYG = PYA + 2*NNB
      PXQ = PYG + 2*NNB
      PXA = PXQ + 2*NNB
      PXG = PXA + 2*NNB
      PXB = PXG + 2*NNB
      PDW = PXB + 2*NNB
C
      IF ( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
C
C        Use unblocked code.
C
         I = ILO
C
      ELSE IF ( LTRA .AND. LTRB ) THEN
         DO 20  I = ILO, N-NX-1, NB
            IB = MIN( NB, N-I )
            NIB = N*IB
C
C           Reduce rows and columns i:i+nb-1 to symplectic URV form and
C           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which
C           are needed to update the unreduced parts of the matrices.
C
            CALL MB03XU( LTRA, LTRB, N-I+1, I-1, IB, A(I,1), LDA,
     $                   B(1,I), LDB, G, LDG, Q(I,I), LDQ, DWORK(PXA),
     $                   N, DWORK(PXB), N, DWORK(PXG), N, DWORK(PXQ), N,
     $                   DWORK(PYA), N, DWORK(PYB), N, DWORK(PYG), N,
     $                   DWORK(PYQ), N, CSL(2*I-1), CSR(2*I-1), TAUL(I),
     $                   TAUR(I), DWORK(PDW) )
C
C           Update the submatrix A(i+1+ib:n,1:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, DWORK(PXA+NB+1), N, Q(I+IB,I), LDQ,
     $                  ONE, A(I+IB+1,I+IB), LDA )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB,
     $                  N-I-IB+1, IB, ONE, DWORK(PXA+NIB+NB+1), N,
     $                  A(I,I+IB), LDA, ONE, A(I+IB+1,I+IB), LDA )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB, N, IB,
     $                  ONE, Q(I,I+IB+1), LDQ, DWORK(PYA), N, ONE,
     $                  A(I+IB+1,1), LDA )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N, IB,
     $                  ONE, B(I+IB+1,I), LDB, DWORK(PYA+NIB), N, ONE,
     $                  A(I+IB+1,1), LDA )
C
C           Update the submatrix Q(i+ib:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXQ+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I,I+IB), LDA, DWORK(PXQ+NIB+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NB), N,
     $                  Q(I,I+IB+1), LDQ, ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NIB+NB), N,
     $                  B(I+IB+1,I), LDB, ONE, Q(I+IB,I+IB+1), LDQ )
C
C           Update the matrix G.
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, Q(I+IB,I), LDQ, DWORK(PXG), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I,I+IB), LDA, DWORK(PXG+NIB), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG), N, Q(I,I+IB+1), LDQ, ONE,
     $                  G(1,I+IB+1), LDG )
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG+NIB), N, B(I+IB+1,I), LDB, ONE,
     $                  G(1,I+IB+1), LDG )
C
C           Update the submatrix B(1:n,i+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB+1,
     $                  IB, ONE, DWORK(PXB), N, Q(I+IB,I), LDQ,
     $                  ONE, B(1,I+IB), LDB )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB+1, IB,
     $                  ONE, DWORK(PXB+NIB), N, A(I,I+IB), LDA, ONE,
     $                  B(1,I+IB), LDB )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, Q(I,I+IB+1), LDQ, DWORK(PYB+NB), N,
     $                  ONE, B(I+IB+1,I+IB), LDB )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, B(I+IB+1,I), LDB, DWORK(PYB+NIB+NB), N,
     $                  ONE, B(I+IB+1,I+IB), LDB )
   20    CONTINUE
C
      ELSE IF ( LTRA ) THEN
         DO 30  I = ILO, N-NX-1, NB
            IB = MIN( NB, N-I )
            NIB = N*IB
C
C           Reduce rows and columns i:i+nb-1 to symplectic URV form and
C           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which
C           are needed to update the unreduced parts of the matrices.
C
            CALL MB03XU( LTRA, LTRB, N-I+1, I-1, IB, A(I,1), LDA,
     $                   B(I,1), LDB, G, LDG, Q(I,I), LDQ, DWORK(PXA),
     $                   N, DWORK(PXB), N, DWORK(PXG), N, DWORK(PXQ), N,
     $                   DWORK(PYA), N, DWORK(PYB), N, DWORK(PYG), N,
     $                   DWORK(PYQ), N, CSL(2*I-1), CSR(2*I-1), TAUL(I),
     $                   TAUR(I), DWORK(PDW) )
C
C           Update the submatrix A(i+1+ib:n,1:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, DWORK(PXA+NB+1), N, Q(I+IB,I), LDQ,
     $                  ONE, A(I+IB+1,I+IB), LDA )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB,
     $                  N-I-IB+1, IB, ONE, DWORK(PXA+NIB+NB+1), N,
     $                  A(I,I+IB), LDA, ONE, A(I+IB+1,I+IB), LDA )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB, N, IB,
     $                  ONE, Q(I,I+IB+1), LDQ, DWORK(PYA), N, ONE,
     $                  A(I+IB+1,1), LDA )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB, N, IB,
     $                  ONE, B(I,I+IB+1), LDB, DWORK(PYA+NIB), N, ONE,
     $                  A(I+IB+1,1), LDA )
C
C           Update the submatrix Q(i+ib:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXQ+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I,I+IB), LDA, DWORK(PXQ+NIB+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NB), N,
     $                  Q(I,I+IB+1), LDQ, ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NIB+NB), N,
     $                  B(I,I+IB+1), LDB, ONE, Q(I+IB,I+IB+1), LDQ )
C
C           Update the matrix G.
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, Q(I+IB,I), LDQ, DWORK(PXG), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I,I+IB), LDA, DWORK(PXG+NIB), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG), N, Q(I,I+IB+1), LDQ, ONE,
     $                  G(1,I+IB+1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG+NIB), N, B(I,I+IB+1), LDB, ONE,
     $                  G(1,I+IB+1), LDG )
C
C           Update the submatrix B(i+ib:n,1:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXB), N,
     $                  ONE, B(I+IB,1), LDB )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I,I+IB), LDA, DWORK(PXB+NIB), N, ONE,
     $                  B(I+IB,1), LDB )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYB+NB), N, Q(I,I+IB+1),
     $                  LDQ, ONE, B(I+IB,I+IB+1), LDB )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYB+NIB+NB), N,
     $                  B(I,I+IB+1), LDB, ONE, B(I+IB,I+IB+1), LDB )
   30    CONTINUE
C
      ELSE IF ( LTRB ) THEN
         DO 40  I = ILO, N-NX-1, NB
            IB = MIN( NB, N-I )
            NIB = N*IB
C
C           Reduce rows and columns i:i+nb-1 to symplectic URV form and
C           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which
C           are needed to update the unreduced parts of the matrices.
C
            CALL MB03XU( LTRA, LTRB, N-I+1, I-1, IB, A(1,I), LDA,
     $                   B(1,I), LDB, G, LDG, Q(I,I), LDQ, DWORK(PXA),
     $                   N, DWORK(PXB), N, DWORK(PXG), N, DWORK(PXQ), N,
     $                   DWORK(PYA), N, DWORK(PYB), N, DWORK(PYG), N,
     $                   DWORK(PYQ), N, CSL(2*I-1), CSR(2*I-1), TAUL(I),
     $                   TAUR(I), DWORK(PDW) )
C
C           Update the submatrix A(1:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXA+NB+1), N,
     $                  ONE, A(I+IB,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I+IB,I), LDA, DWORK(PXA+NIB+NB+1), N,
     $                  ONE, A(I+IB,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYA), N, Q(I,I+IB+1), LDQ, ONE,
     $                  A(1,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYA+NIB), N, B(I+IB+1,I), LDB, ONE,
     $                  A(1,I+IB+1), LDA )
C
C           Update the submatrix Q(i+ib:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXQ+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I+IB,I), LDA, DWORK(PXQ+NIB+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NB), N,
     $                  Q(I,I+IB+1), LDQ, ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No Transpose', 'Transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NIB+NB), N,
     $                  B(I+IB+1,I), LDB, ONE, Q(I+IB,I+IB+1), LDQ )
C
C           Update the matrix G.
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, Q(I+IB,I), LDQ, DWORK(PXG), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I+IB,I), LDA, DWORK(PXG+NIB), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG), N, Q(I,I+IB+1), LDQ, ONE,
     $                  G(1,I+IB+1), LDG )
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG+NIB), N, B(I+IB+1,I), LDB, ONE,
     $                  G(1,I+IB+1), LDG )
C
C           Update the submatrix B(1:n,i+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB+1,
     $                  IB, ONE, DWORK(PXB), N, Q(I+IB,I), LDQ,
     $                  ONE, B(1,I+IB), LDB )
            CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB+1, IB,
     $                  ONE, DWORK(PXB+NIB), N, A(I+IB,I), LDA, ONE,
     $                  B(1,I+IB), LDB )
            CALL DGEMM( 'Transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, Q(I,I+IB+1), LDQ, DWORK(PYB+NB), N,
     $                  ONE, B(I+IB+1,I+IB), LDB )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB+1,
     $                  IB, ONE, B(I+IB+1,I), LDB, DWORK(PYB+NIB+NB), N,
     $                  ONE, B(I+IB+1,I+IB), LDB )
   40    CONTINUE
C
      ELSE
         DO 50  I = ILO, N-NX-1, NB
            IB = MIN( NB, N-I )
            NIB = N*IB
C
C           Reduce rows and columns i:i+nb-1 to symplectic URV form and
C           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which
C           are needed to update the unreduced parts of the matrices.
C
            CALL MB03XU( LTRA, LTRB, N-I+1, I-1, IB, A(1,I), LDA,
     $                   B(I,1), LDB, G, LDG, Q(I,I), LDQ, DWORK(PXA),
     $                   N, DWORK(PXB), N, DWORK(PXG), N, DWORK(PXQ), N,
     $                   DWORK(PYA), N, DWORK(PYB), N, DWORK(PYG), N,
     $                   DWORK(PYQ), N, CSL(2*I-1), CSR(2*I-1), TAUL(I),
     $                   TAUR(I), DWORK(PDW) )
C
C           Update the submatrix A(1:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXA+NB+1), N,
     $                  ONE, A(I+IB,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I+IB,I), LDA, DWORK(PXA+NIB+NB+1), N,
     $                  ONE, A(I+IB,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYA), N, Q(I,I+IB+1), LDQ, ONE,
     $                  A(1,I+IB+1), LDA )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYA+NIB), N, B(I,I+IB+1), LDB, ONE,
     $                  A(1,I+IB+1), LDA )
C
C           Update the submatrix Q(i+ib:n,i+1+ib:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXQ+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N-I-IB,
     $                  IB, ONE, A(I+IB,I), LDA, DWORK(PXQ+NIB+NB+1), N,
     $                  ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NB), N,
     $                  Q(I,I+IB+1), LDQ, ONE, Q(I+IB,I+IB+1), LDQ )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYQ+NIB+NB), N,
     $                  B(I,I+IB+1), LDB, ONE, Q(I+IB,I+IB+1), LDQ )
C
C           Update the matrix G.
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, Q(I+IB,I), LDQ, DWORK(PXG), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I+IB,I), LDA, DWORK(PXG+NIB), N, ONE,
     $                  G(I+IB,1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG), N, Q(I,I+IB+1), LDQ, ONE,
     $                  G(1,I+IB+1), LDG )
            CALL DGEMM( 'No transpose', 'No transpose', N, N-I-IB, IB,
     $                  ONE, DWORK(PYG+NIB), N, B(I,I+IB+1), LDB, ONE,
     $                  G(1,I+IB+1), LDG )
C
C           Update the submatrix B(i+ib:n,1:n).
C
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N,
     $                  IB, ONE, Q(I+IB,I), LDQ, DWORK(PXB), N,
     $                  ONE, B(I+IB,1), LDB )
            CALL DGEMM( 'No transpose', 'Transpose', N-I-IB+1, N, IB,
     $                  ONE, A(I+IB,I), LDA, DWORK(PXB+NIB), N, ONE,
     $                  B(I+IB,1), LDB )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYB+NB), N, Q(I,I+IB+1),
     $                  LDQ, ONE, B(I+IB,I+IB+1), LDB )
            CALL DGEMM( 'No transpose', 'No transpose', N-I-IB+1,
     $                  N-I-IB, IB, ONE, DWORK(PYB+NIB+NB), N,
     $                  B(I,I+IB+1), LDB, ONE, B(I+IB,I+IB+1), LDB )
   50    CONTINUE
      END IF
C
C     Unblocked code to reduce the rest of the matrices.
C
      CALL MB04TS( TRANA, TRANB, N, I, A, LDA, B, LDB, G, LDG, Q, LDQ,
     $             CSL, CSR, TAUL, TAUR, DWORK, LDWORK, IERR )
C
      DWORK(1) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of MB04TB ***
      END
