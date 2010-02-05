      SUBROUTINE MB04TS( TRANA, TRANB, N, ILO, A, LDA, B, LDB, G, LDG,
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
C     2N-by-2N matrix H:
C
C             [ op(A)   G   ]        T       [ op(R11)   R12   ]    T
C         H = [             ] = U R V  = U * [                 ] * V ,
C             [  Q    op(B) ]                [   0     op(R22) ]
C
C     where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real
C     N-by-N upper triangular matrix, op(R22) is a real N-by-N lower
C     Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic
C     matrices. Unblocked version.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TRANA   CHARACTER*1
C             Specifies the form of op( A ) as follows:
C             = 'N': op( A ) = A;
C             = 'T': op( A ) = A';
C             = 'C': op( A ) = A'.
C
C     TRANB   CHARACTER*1
C             Specifies the form of op( B ) as follows:
C             = 'N': op( B ) = B;
C             = 'T': op( B ) = B';
C             = 'C': op( B ) = B'.
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
C             The leading dimension of the array Q.  LDG >= MAX(1,N).
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
C             applied from the left-hand side.
C
C     TAUR    (output) DOUBLE PRECISION array, dimension (N-1)
C             On exit, the first N-1 elements of this array contain the
C             scalar factors of some of the elementary reflectors
C             applied from the right-hand side.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
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
C     The algorithm requires 80/3 N**3 + 20 N**2 + O(N) floating point
C     operations and is numerically backward stable.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A numerically stable, structure preserving method for
C         computing the eigenvalues of real Hamiltonian or symplectic
C         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DGESUV).
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
      INTEGER           I
      DOUBLE PRECISION  ALPHA, C, NU, S, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLARF, DLARFG, DLARTG, DROT, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO  = 0
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
         CALL XERBLA( 'MB04TS', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      DO 10 I = ILO, N
         ALPHA = Q(I,I)
         IF ( I.LT.N ) THEN
C
C           Generate elementary reflector HU(i) to annihilate Q(i+1:n,i)
C
            CALL DLARFG( N-I+1, ALPHA, Q(I+1,I), 1, NU )
C
C           Apply HU(i) from the left.
C
            Q(I,I) = ONE
            CALL DLARF( 'Left', N-I+1, N-I, Q(I,I), 1, NU, Q(I,I+1),
     $                  LDQ, DWORK )
            IF ( LTRA ) THEN
               CALL DLARF( 'Right', N-I+1, N-I+1, Q(I,I), 1, NU, A(I,I),
     $                     LDA, DWORK )
            ELSE
               CALL DLARF( 'Left', N-I+1, N-I+1, Q(I,I), 1, NU, A(I,I),
     $                     LDA, DWORK )
            END IF
            IF ( LTRB ) THEN
               CALL DLARF( 'Right', N, N-I+1, Q(I,I), 1, NU, B(1,I),
     $                     LDB, DWORK )
            ELSE
               CALL DLARF( 'Left', N-I+1, N, Q(I,I), 1, NU, B(I,1), LDB,
     $                     DWORK )
            END IF
            CALL DLARF( 'Left', N-I+1, N, Q(I,I), 1, NU, G(I,1), LDG,
     $                  DWORK )
            Q(I,I) = NU
         ELSE
            Q(I,I) = ZERO
         END IF
C
C        Generate symplectic Givens rotator GU(i) to annihilate Q(i,i).
C
         TEMP = A(I,I)
         CALL DLARTG( TEMP, ALPHA, C, S, A(I,I) )
C
C        Apply G(i) from the left.
C
         IF ( LTRA ) THEN
            CALL DROT( N-I, A(I+1,I), 1, Q(I,I+1), LDQ, C, S )
         ELSE
            CALL DROT( N-I, A(I,I+1), LDA, Q(I,I+1), LDQ, C, S )
         END IF
         IF ( LTRB ) THEN
            CALL DROT( N, G(I,1), LDG, B(1,I), 1, C, S )
         ELSE
            CALL DROT( N, G(I,1), LDG, B(I,1), LDB, C, S )
         END IF
         CSL(2*I-1) = C
         CSL(2*I)   = S
C
         IF ( I.LT.N ) THEN
            IF ( LTRA ) THEN
C
C              Generate elementary reflector FU(i) to annihilate
C              A(i,i+1:n).
C
               CALL DLARFG( N-I+1, A(I,I), A(I,I+1), LDA, TAUL(I) )
C
C              Apply FU(i) from the left.
C
               TEMP = A(I,I)
               A(I,I) = ONE
               CALL DLARF( 'Right', N-I, N-I+1, A(I,I), LDA, TAUL(I),
     $                     A(I+1,I), LDA, DWORK )
               CALL DLARF( 'Left', N-I+1, N-I, A(I,I), LDA, TAUL(I),
     $                     Q(I,I+1), LDQ, DWORK )
               IF ( LTRB ) THEN
                  CALL DLARF( 'Right', N, N-I+1, A(I,I), LDA, TAUL(I),
     $                        B(1,I), LDB, DWORK )
               ELSE
                  CALL DLARF( 'Left', N-I+1, N, A(I,I), LDA, TAUL(I),
     $                        B(I,1), LDB, DWORK )
               END IF
               CALL DLARF( 'Left', N-I+1, N, A(I,I), LDA, TAUL(I),
     $                     G(I,1), LDG, DWORK )
               A(I,I) = TEMP
            ELSE
C
C              Generate elementary reflector FU(i) to annihilate
C              A(i+1:n,i).
C
               CALL DLARFG( N-I+1, A(I,I), A(I+1,I), 1, TAUL(I) )
C
C              Apply FU(i) from the left.
C
               TEMP = A(I,I)
               A(I,I) = ONE
               CALL DLARF( 'Left', N-I+1, N-I, A(I,I), 1, TAUL(I),
     $                     A(I,I+1), LDA, DWORK )
               CALL DLARF( 'Left', N-I+1, N-I, A(I,I), 1, TAUL(I),
     $                     Q(I,I+1), LDQ, DWORK )
               IF ( LTRB ) THEN
                  CALL DLARF( 'Right', N, N-I+1, A(I,I), 1, TAUL(I),
     $                        B(1,I), LDB, DWORK )
               ELSE
                  CALL DLARF( 'Left', N-I+1, N, A(I,I), 1, TAUL(I),
     $                        B(I,1), LDB, DWORK )
               END IF
               CALL DLARF( 'Left', N-I+1, N, A(I,I), 1, TAUL(I), G(I,1),
     $                     LDG, DWORK )
               A(I,I) = TEMP
            END IF
         ELSE
            TAUL(I) = ZERO
         END IF
         IF ( I.LT.N )
     $      ALPHA = Q(I,I+1)
         IF ( I.LT.N-1 ) THEN
C
C           Generate elementary reflector HV(i) to annihilate Q(i,i+2:n)
C
            CALL DLARFG( N-I, ALPHA, Q(I,I+2), LDQ, NU )
C
C           Apply HV(i) from the right.
C
            Q(I,I+1) = ONE
            CALL DLARF( 'Right', N-I, N-I, Q(I,I+1), LDQ, NU,
     $                  Q(I+1,I+1), LDQ, DWORK )
            IF ( LTRA ) THEN
               CALL DLARF( 'Left', N-I, N, Q(I,I+1), LDQ, NU,
     $                     A(I+1,1), LDA, DWORK )
            ELSE
               CALL DLARF( 'Right', N, N-I, Q(I,I+1), LDQ, NU,
     $                     A(1,I+1), LDA, DWORK )
            END IF
            IF ( LTRB ) THEN
               CALL DLARF( 'Left', N-I, N-I+1, Q(I,I+1), LDQ, NU,
     $                     B(I+1,I), LDB, DWORK )
            ELSE
               CALL DLARF( 'Right', N-I+1, N-I, Q(I,I+1), LDQ, NU,
     $                     B(I,I+1), LDB, DWORK )
            END IF
            CALL DLARF( 'Right', N, N-I, Q(I,I+1), LDQ, NU,
     $                  G(1,I+1), LDG, DWORK )
            Q(I,I+1) = NU
         ELSE IF ( I.LT.N ) THEN
            Q(I,I+1) = ZERO
         END IF
         IF ( I.LT.N ) THEN
C
C           Generate symplectic Givens rotator GV(i) to annihilate
C           Q(i,i+1).
C
            IF ( LTRB ) THEN
               TEMP = B(I+1,I)
               CALL DLARTG( TEMP, ALPHA, C, S, B(I+1,I) )
               S = -S
               CALL DROT( N-I, Q(I+1,I+1), 1, B(I+1,I+1), LDB, C, S )
            ELSE
               TEMP = B(I,I+1)
               CALL DLARTG( TEMP, ALPHA, C, S, B(I,I+1) )
               S = -S
               CALL DROT( N-I, Q(I+1,I+1), 1, B(I+1,I+1), 1, C, S )
            END IF
            IF ( LTRA ) THEN
               CALL DROT( N, A(I+1,1), LDA, G(1,I+1), 1, C, S )
            ELSE
               CALL DROT( N, A(1,I+1), 1, G(1,I+1), 1, C, S )
            END IF
            CSR(2*I-1) = C
            CSR(2*I)   = S
         END IF
         IF ( I.LT.N-1 ) THEN
            IF ( LTRB ) THEN
C
C              Generate elementary reflector FV(i) to annihilate
C              B(i+2:n,i).
C
               CALL DLARFG( N-I, B(I+1,I), B(I+2,I), 1, TAUR(I) )
C
C              Apply FV(i) from the right.
C
               TEMP = B(I+1,I)
               B(I+1,I) = ONE
               CALL DLARF( 'Left', N-I, N-I, B(I+1,I), 1, TAUR(I),
     $                     B(I+1,I+1), LDB, DWORK )
               CALL DLARF( 'Right', N-I, N-I, B(I+1,I), 1, TAUR(I),
     $                     Q(I+1,I+1), LDQ, DWORK )
               IF ( LTRA ) THEN
                  CALL DLARF( 'Left', N-I, N, B(I+1,I), 1,
     $                        TAUR(I), A(I+1,1), LDA, DWORK )
               ELSE
                  CALL DLARF( 'Right', N, N-I, B(I+1,I), 1,
     $                        TAUR(I), A(1,I+1), LDA, DWORK )
               END IF
               CALL DLARF( 'Right', N, N-I, B(I+1,I), 1, TAUR(I),
     $                     G(1,I+1), LDG, DWORK )
               B(I+1,I) = TEMP
            ELSE
C
C              Generate elementary reflector FV(i) to annihilate
C              B(i,i+2:n).
C
               CALL DLARFG( N-I, B(I,I+1), B(I,I+2), LDB, TAUR(I) )
C
C              Apply FV(i) from the right.
C
               TEMP = B(I,I+1)
               B(I,I+1) = ONE
               CALL DLARF( 'Right', N-I, N-I, B(I,I+1), LDB, TAUR(I),
     $                     B(I+1,I+1), LDB, DWORK )
               CALL DLARF( 'Right', N-I, N-I, B(I,I+1), LDB, TAUR(I),
     $                     Q(I+1,I+1), LDQ, DWORK )
               IF ( LTRA ) THEN
                  CALL DLARF( 'Left', N-I, N, B(I,I+1), LDB, TAUR(I),
     $                        A(I+1,1), LDA, DWORK )
               ELSE
                  CALL DLARF( 'Right', N, N-I, B(I,I+1), LDB,
     $                        TAUR(I), A(1,I+1), LDA, DWORK )
               END IF
               CALL DLARF( 'Right', N, N-I, B(I,I+1), LDB, TAUR(I),
     $                     G(1,I+1), LDG, DWORK )
               B(I,I+1) = TEMP
            END IF
         ELSE IF ( I.LT.N ) THEN
            TAUR(I) = ZERO
         END IF
   10 CONTINUE
      DWORK(1) = DBLE( MAX( 1, N ) )
      RETURN
C *** Last line of MB04TS ***
      END
