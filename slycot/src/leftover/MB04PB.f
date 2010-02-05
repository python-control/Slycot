      SUBROUTINE MB04PB( N, ILO, A, LDA, QG, LDQG, CS, TAU, DWORK,
     $                   LDWORK, INFO )
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
C     To reduce a Hamiltonian matrix,
C
C                   [  A   G  ]
C              H =  [       T ] ,
C                   [  Q  -A  ]
C
C     where A is an N-by-N matrix and G,Q are N-by-N symmetric matrices,
C     to Paige/Van Loan (PVL) form. That is, an orthogonal symplectic U
C     is computed so that
C
C               T       [  Aout   Gout  ]
C              U H U =  [             T ] ,
C                       [  Qout  -Aout  ]
C
C     where Aout is upper Hessenberg and Qout is diagonal.
C     Blocked version.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     ILO     (input) INTEGER
C             It is assumed that A is already upper triangular and Q is
C             zero in rows and columns 1:ILO-1. ILO is normally set by a
C             previous call to MB04DD; otherwise it should be set to 1.
C             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Aout and, in the zero part of Aout,
C             information about the elementary reflectors used to
C             compute the PVL factorization.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     QG      (input/output) DOUBLE PRECISION array, dimension
C                            (LDQG,N+1)
C             On entry, the leading N-by-N+1 part of this array must
C             contain the lower triangular part of the matrix Q and
C             the upper triangular part of the matrix G.
C             On exit, the leading N-by-N+1 part of this array contains
C             the diagonal of the matrix Qout, the upper triangular part
C             of the matrix Gout and, in the zero parts of Qout,
C             information about the elementary reflectors used to
C             compute the PVL factorization.
C
C     LDQG    INTEGER
C             The leading dimension of the array QG.  LDQG >= MAX(1,N).
C
C     CS      (output) DOUBLE PRECISION array, dimension (2N-2)
C             On exit, the first 2N-2 elements of this array contain the
C             cosines and sines of the symplectic Givens rotations used
C             to compute the PVL factorization.
C
C     TAU     (output) DOUBLE PRECISION array, dimension (N-1)
C             On exit, the first N-1 elements of this array contain the
C             scalar factors of some of the elementary reflectors.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK, 8*N*NB + 3*NB, where NB is the optimal
C             block size determined by the function UE01MD.
C             On exit, if  INFO = -10,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N-1).
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
C     The matrix U is represented as a product of symplectic reflectors
C     and Givens rotators
C
C     U = diag( H(1),H(1) )     G(1)   diag( F(1),F(1) )
C         diag( H(2),H(2) )     G(2)   diag( F(2),F(2) )
C                                ....
C         diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ).
C
C     Each H(i) has the form
C
C           H(i) = I - tau * v * v'
C
C     where tau is a real scalar, and v is a real vector with
C     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
C     QG(i+2:n,i), and tau in QG(i+1,i).
C
C     Each F(i) has the form
C
C           F(i) = I - nu * w * w'
C
C     where nu is a real scalar, and w is a real vector with
C     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in
C     A(i+2:n,i), and nu in TAU(i).
C
C     Each G(i) is a Givens rotator acting on rows i+1 and n+i+1,
C     where the cosine is stored in CS(2*i-1) and the sine in
C     CS(2*i).
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires O(N**3) floating point operations and is
C     strongly backward stable.
C
C     REFERENCES
C
C     [1] C. F. VAN LOAN:
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications, 61, pp. 233-251, 1984.
C
C     [2] D. KRESSNER:
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     REVISIONS
C
C     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DHAPVB).
C
C     KEYWORDS
C
C     Elementary matrix operations, Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           ILO, INFO, LDA, LDQG, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), CS(*), DWORK(*), QG(LDQG,*), TAU(*)
C     .. Local Scalars ..
      INTEGER           I, IB, IERR, NB, NBMIN, NH, NIB, NNB, NX, PDW,
     $                  PXA, PXG, PXQ, PYA, WRKOPT
C     .. External Functions ..
      INTEGER           UE01MD
      EXTERNAL          UE01MD
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DSYR2K, MB04PA, MB04PU, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDQG.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDWORK.LT.MAX( 1, N-1 ) ) THEN
         DWORK(1) = DBLE( MAX( 1, N-1 ) )
         INFO = -10
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04PB', -INFO )
         RETURN
      END IF
C
C     Set elements 1:ILO-1 of TAU and CS.
C
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
         CS(2*I-1) = ONE
         CS(2*I) = ZERO
   10 CONTINUE
C
C     Quick return if possible.
C
      IF ( N.LE.ILO ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Determine the block size.
C
      NH = N - ILO + 1
      NB = UE01MD( 1, 'MB04PB', ' ', N, ILO, -1 )
      NBMIN = 2
      WRKOPT = N-1
      IF ( NB.GT.1 .AND. NB.LT.NH ) THEN
C
C        Determine when to cross over from blocked to unblocked code.
C
         NX = MAX( NB, UE01MD( 3, 'MB04PB', ' ', N, ILO, -1 ) )
         IF ( NX.LT.NH ) THEN
C
C           Check whether workspace is large enough for blocked code.
C
            WRKOPT = 8*N*NB + 3*NB
            IF ( LDWORK.LT.WRKOPT ) THEN
C
C              Not enough workspace available. Determine minimum value
C              of NB, and reduce NB.
C
               NBMIN = MAX( 2, UE01MD( 2, 'MB04PB', ' ', N, ILO, -1 ) )
               NB = LDWORK / ( 8*N + 3 )
            END IF
         END IF
      END IF
C
      NNB = N*NB
      PXA = 1
      PYA = PXA + 2*NNB
      PXQ = PYA + 2*NNB
      PXG = PXQ + 2*NNB
      PDW = PXG + 2*NNB
C
      IF ( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
C
C        Use unblocked code.
C
         I = ILO
C
      ELSE
         DO 20  I = ILO, N-NX-1, NB
            IB = MIN( NB, N-I )
            NIB = N*IB
C
C           Reduce rows and columns i:i+nb-1 to PVL form and return the
C           matrices XA, XG, XQ, and YA which are needed to update the
C           unreduced parts of the matrices.
C
            CALL MB04PA( .TRUE., N-I+1, I-1, IB, A(1,I), LDA, QG(1,I),
     $                   LDQG, DWORK(PXA), N, DWORK(PXG), N,
     $                   DWORK(PXQ), N, DWORK(PYA), N, CS(2*I-1),
     $                   TAU(I), DWORK(PDW) )
            IF ( N.GT.I+IB ) THEN
C
C              Update the submatrix A(1:n,i+ib+1:n).
C
               CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB,
     $                     IB, ONE, QG(I+IB+1,I), LDQG, DWORK(PXA+IB+1),
     $                     N, ONE, A(I+IB+1,I+IB+1), LDA )
               CALL DGEMM( 'No transpose', 'Transpose', N-I-IB, N-I-IB,
     $                     IB, ONE, A(I+IB+1,I), LDA,
     $                     DWORK(PXA+NIB+IB+1), N, ONE,
     $                     A(I+IB+1,I+IB+1), LDA )
               CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB, IB,
     $                     ONE, DWORK(PYA), N, QG(I+IB+1,I), LDQG, ONE,
     $                     A(1,I+IB+1), LDA )
               CALL DGEMM( 'No transpose', 'Transpose', N, N-I-IB, IB,
     $                     ONE, DWORK(PYA+NIB), N, A(I+IB+1,I), LDA,
     $                     ONE, A(1,I+IB+1), LDA )
C
C              Update the submatrix Q(i+ib+1:n,i+ib+1:n).
C
               CALL DSYR2K( 'Lower', 'No Transpose', N-I-IB, IB, ONE,
     $                      DWORK(PXQ+IB+1), N, QG(I+IB+1,I), LDQG, ONE,
     $                      QG(I+IB+1,I+IB+1), LDQG )
               CALL DSYR2K( 'Lower', 'No Transpose', N-I-IB, IB, ONE,
     $                      DWORK(PXQ+NIB+IB+1), N, A(I+IB+1,I), LDA,
     $                      ONE, QG(I+IB+1,I+IB+1), LDQG )
C
C              Update the submatrix G(1:n,1:n).
C
               CALL DGEMM( 'No transpose', 'Transpose', I+IB, N-I-IB,
     $                     IB, ONE, DWORK(PXG), N, QG(I+IB+1,I), LDQG,
     $                     ONE, QG(1,I+IB+2), LDQG )
               CALL DGEMM( 'No transpose', 'Transpose', I+IB, N-I-IB,
     $                     IB, ONE, DWORK(PXG+NIB), N, A(I+IB+1,I), LDA,
     $                     ONE, QG(1,I+IB+2), LDQG )
               CALL DSYR2K( 'Upper', 'No Transpose', N-I-IB, IB, ONE,
     $                      DWORK(PXG+IB+I), N, QG(I+IB+1,I), LDQG, ONE,
     $                      QG(I+IB+1,I+IB+2), LDQG )
               CALL DSYR2K( 'Upper', 'No Transpose', N-I-IB, IB, ONE,
     $                      DWORK(PXG+NIB+IB+I), N, A(I+IB+1,I), LDA,
     $                      ONE, QG(I+IB+1,I+IB+2), LDQG )
            END IF
   20    CONTINUE
      END IF
C
C     Unblocked code to reduce the rest of the matrices.
C
      CALL MB04PU( N, I, A, LDA, QG, LDQG, CS, TAU, DWORK, LDWORK,
     $             IERR )
C
      DWORK( 1 ) = DBLE( WRKOPT )
C
      RETURN
C *** Last line of MB04PB ***
      END
