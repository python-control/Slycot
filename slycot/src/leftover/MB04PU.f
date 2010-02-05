      SUBROUTINE MB04PU( N, ILO, A, LDA, QG, LDQG, CS, TAU, DWORK,
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
C     Unblocked version.
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
C             value of LDWORK.
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
C     The algorithm requires 40/3 N**3 + O(N) floating point operations
C     and is strongly backward stable.
C
C     REFERENCES
C
C     [1] C. F. VAN LOAN:
C         A symplectic method for approximating all the eigenvalues of
C         a Hamiltonian matrix.
C         Linear Algebra and its Applications, 61, pp. 233-251, 1984.
C
C     CONTRIBUTORS
C
C     D. Kressner (Technical Univ. Berlin, Germany) and
C     P. Benner (Technical Univ. Chemnitz, Germany), December 2003.
C
C     REVISIONS
C
C     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DHAPVL).
C
C     KEYWORDS
C
C     Elementary matrix operations, Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           ILO, INFO, LDA, LDQG, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), CS(*), DWORK(*), QG(LDQG,*), TAU(*)
C     .. Local Scalars ..
      INTEGER           I
      DOUBLE PRECISION  ALPHA, C, MU, NU, S, TEMP, TTEMP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLARF, DLARFG, DLARTG, DROT, DSYMV,
     $                  DSYR2, XERBLA
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
         CALL XERBLA( 'MB04PU', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.LE.ILO ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      DO 10 I = ILO, N-1
C
C        Generate elementary reflector H(i) to annihilate QG(i+2:n,i).
C
         ALPHA = QG(I+1,I)
         CALL DLARFG( N-I, ALPHA, QG(MIN( I+2,N ),I), 1, NU )
         IF ( NU.NE.ZERO ) THEN
            QG(I+1,I) = ONE
C
C           Apply H(i) from both sides to QG(i+1:n,i+1:n).
C           Compute  x := nu * QG(i+1:n,i+1:n) * v.
C
            CALL DSYMV( 'Lower', N-I, NU, QG(I+1,I+1), LDQG, QG(I+1,I),
     $                  1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( N-I, DWORK, 1, QG(I+1,I), 1 )
            CALL DAXPY( N-I, MU, QG(I+1,I), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C                QG := QG - v * w' - w * v'.
C
            CALL DSYR2( 'Lower', N-I, -ONE, QG(I+1,I), 1, DWORK, 1,
     $                  QG(I+1,I+1), LDQG )
C
C           Apply H(i) from the right hand side to QG(1:i,i+2:n+1).
C
            CALL DLARF( 'Right', I, N-I, QG(I+1,I), 1, NU, QG(1,I+2),
     $                  LDQG, DWORK )
C
C           Apply H(i) from both sides to QG(i+1:n,i+2:n+1).
C           Compute  x := nu * QG(i+1:n,i+2:n+1) * v.
C
            CALL DSYMV( 'Upper', N-I, NU, QG(I+1,I+2), LDQG, QG(I+1,I),
     $                  1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * nu * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( N-I, DWORK, 1, QG(I+1,I), 1 )
            CALL DAXPY( N-I, MU, QG(I+1,I), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'.
C
            CALL DSYR2( 'Upper', N-I, -ONE, QG(I+1,I), 1, DWORK, 1,
     $                  QG(I+1,I+2), LDQG )
C
C           Apply H(i) from the left hand side to A(i+1:n,i:n).
C
            CALL DLARF( 'Left', N-I, N-I+1, QG(I+1,I), 1, NU,
     $                  A(I+1,I), LDA, DWORK )
C
C           Apply H(i) from the right hand side to A(1:n,i+1:n).
C
            CALL DLARF( 'Right', N, N-I,  QG(I+1,I), 1, NU,
     $                  A(1,I+1), LDA, DWORK )
         END IF
         QG(I+1,I) = NU
C
C        Generate symplectic Givens rotation G(i) to annihilate
C        QG(i+1,i).
C
         TEMP = A(I+1,I)
         CALL DLARTG( TEMP, ALPHA, C, S, A(I+1,I) )
C
C        Apply G(i) to [A(I+1,I+2:N); QG(I+2:N,I+1)'].
C
         CALL DROT( N-I-1, A(I+1,I+2), LDA, QG(I+2,I+1), 1, C, S )
C
C        Apply G(i) to [A(1:I,I+1) QG(1:I,I+2)].
C
         CALL DROT(I, A(1,I+1), 1, QG(1,I+2), 1, C, S )
C
C        Apply G(i) to [A(I+2:N,I+1) QG(I+1, I+3:N+1)'] from the right.
C
         CALL DROT(N-I-1, A(I+2,I+1), 1, QG(I+1,I+3), LDQG, C, S )
C
C        Fix the diagonal part.
C
         TEMP  = A(I+1,I+1)
         TTEMP = QG(I+1,I+2)
         A(I+1,I+1)  = C*TEMP  + S*QG(I+1,I+1)
         QG(I+1,I+2) = C*TTEMP - S * TEMP
         QG(I+1,I+1) = -S*TEMP + C*QG(I+1,I+1)
         TTEMP = -S*TTEMP - C*TEMP
         TEMP = A(I+1,I+1)
         QG(I+1,I+1) =  C*QG(I+1,I+1) + S*TTEMP
         A(I+1,I+1)  =  C*TEMP + S*QG(I+1,I+2)
         QG(I+1,I+2) = -S*TEMP + C*QG(I+1,I+2)
         CS(2*I-1) = C
         CS(2*I)   = S
C
C        Generate elementary reflector F(i) to annihilate A(i+2:n,i).
C
         CALL DLARFG( N-I, A(I+1,I), A(MIN( I+2,N ),I), 1, NU )
         IF ( NU.NE.ZERO ) THEN
            TEMP = A(I+1,I)
            A(I+1,I) = ONE
C
C           Apply F(i) from the left hand side to A(i+1:n,i+1:n).
C
            CALL DLARF( 'Left', N-I, N-I, A(I+1,I), 1, NU, A(I+1,I+1),
     $                  LDA, DWORK )
C
C           Apply G(i) from the right hand side to A(1:n,i+1:n).
C
            CALL DLARF( 'Right', N, N-I,  A(I+1,I), 1, NU,
     $                  A(1,I+1), LDA, DWORK )
C
C           Apply G(i) from both sides to QG(i+1:n,i+1:n).
C           Compute  x := nu * QG(i+1:n,i+1:n) * v.
C
            CALL DSYMV( 'Lower', N-I, NU, QG(I+1,I+1), LDQG, A(I+1,I),
     $                  1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * tau * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( N-I, DWORK, 1, A(I+1,I), 1 )
            CALL DAXPY( N-I, MU, A(I+1,I), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C                QG := QG - v * w' - w * v'.
C
            CALL DSYR2( 'Lower', N-I, -ONE, A(I+1,I), 1, DWORK, 1,
     $                  QG(I+1,I+1), LDQG )
C
C           Apply G(i) from the right hand side to QG(1:i,i+2:n+1).
C
            CALL DLARF( 'Right', I, N-I, A(I+1,I), 1, NU, QG(1,I+2),
     $                  LDQG, DWORK )
C
C           Apply G(i) from both sides to QG(i+1:n,i+2:n+1).
C           Compute  x := nu * QG(i+1:n,i+2:n+1) * v.
C
            CALL DSYMV( 'Upper', N-I, NU, QG(I+1,I+2), LDQG, A(I+1,I),
     $                  1, ZERO, DWORK, 1 )
C
C           Compute  w := x - 1/2 * tau * (x'*v) * v.
C
            MU = -HALF*NU*DDOT( N-I, DWORK, 1, A(I+1,I), 1 )
            CALL DAXPY( N-I, MU, A(I+1,I), 1, DWORK, 1 )
C
C           Apply the transformation as a rank-2 update:
C              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'.
C
            CALL DSYR2( 'Upper', N-I, -ONE, A(I+1,I), 1, DWORK, 1,
     $                  QG(I+1,I+2), LDQG )
            A(I+1,I) = TEMP
         END IF
         TAU(I) = NU
   10 CONTINUE
      DWORK(1) = DBLE( MAX( 1, N-1 ) )
      RETURN
C *** Last line of MB04PU ***
      END
