      SUBROUTINE SB04OW( M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE,
     $                   F, LDF, SCALE, IWORK, INFO )
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
C     To solve a periodic Sylvester equation
C
C              A * R - L * B = scale * C                           (1)
C              D * L - R * E = scale * F,
C
C     using Level 1 and 2 BLAS, where R and L are unknown M-by-N
C     matrices, (A, D), (B, E) and (C, F) are given matrix pairs of
C     size M-by-M, N-by-N and M-by-N, respectively, with real entries.
C     (A, D) and (B, E) must be in periodic Schur form, i.e. A, B are
C     upper quasi triangular and D, E are upper triangular. The solution
C     (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output scaling
C     factor chosen to avoid overflow.
C
C     This routine is largely based on the LAPACK routine DTGSY2
C     developed by Bo Kagstrom and Peter Poromaa.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The order of A and D, and the row dimension of C, F, R
C             and L.  M >= 0.
C
C     N       (input) INTEGER
C             The order of B and E, and the column dimension of C, F, R
C             and L.  N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,M)
C             On entry, the leading M-by-M part of this array must
C             contain the upper quasi triangular matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,M).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper quasi triangular matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand-side of the first matrix equation
C             in (1).
C             On exit, the leading M-by-N part of this array contains
C             the solution R.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,M).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading M-by-M part of this array must
C             contain the upper triangular matrix D.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,M).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular matrix E.
C
C     LDE     INTEGER
C             The leading dimension of the array E.  LDE >= MAX(1,N).
C
C     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N)
C             On entry, the leading M-by-N part of this array must
C             contain the right-hand-side of the second matrix equation
C             in (1).
C             On exit, the leading M-by-N part of this array contains
C             the solution L.
C
C     LDF     INTEGER
C             The leading dimension of the array F.  LDF >= MAX(1,M).
C
C     SCALE   (output) DOUBLE PRECISION
C             On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the arrays
C             C and F will hold the solutions R and L, respectively, to
C             a slightly perturbed system but the input matrices A, B, D
C             and E have not been changed. If SCALE = 0, C and F will
C             hold solutions to the homogeneous system with C = F = 0.
C             Normally, SCALE = 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (M+N+2)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  the matrix products A*D and B*E have common or very
C                   close eigenvalues.
C
C     METHOD
C
C     In matrix notation solving equation (1) corresponds to solving
C     Z*x = scale*b, where Z is defined as
C
C         Z = [  kron(In, A)  -kron(B', Im) ]            (2)
C             [ -kron(E', Im)  kron(In, D)  ],
C
C     Ik is the identity matrix of size k and X' is the transpose of X.
C     kron(X, Y) is the Kronecker product between the matrices X and Y.
C     In the process of solving (1), we solve a number of such systems
C     where Dim(Im), Dim(In) = 1 or 2.
C
C     REFERENCES
C
C     [1] Kagstrom, B.
C         A Direct Method for Reordering Eigenvalues in the Generalized
C         Real Schur Form of a Regular Matrix Pair (A,B). M.S. Moonen
C         et al (eds.), Linear Algebra for Large Scale and Real-Time
C         Applications, Kluwer Academic Publ., pp. 195-218, 1993.
C
C     [2] Sreedhar, J. and Van Dooren, P.
C         A Schur approach for solving some periodic matrix equations.
C         U. Helmke et al (eds.), Systems and Networks: Mathematical
C         Theory and Applications, Akademie Verlag, Berlin, vol. 77,
C         pp. 339-362, 1994.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DTGPY2).
C
C     KEYWORDS
C
C     Matrix equation, periodic Sylvester equation.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDZ
      PARAMETER         ( LDZ = 8 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N
      DOUBLE PRECISION  SCALE
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  E(LDE,*), F(LDF,*)
C     .. Local Scalars ..
      INTEGER           I, IE, IERR, II, IS, ISP1, J, JE, JJ, JS, JSP1,
     $                  K, MB, NB, P, Q, ZDIM
      DOUBLE PRECISION  SCALOC
C     .. Local Arrays ..
      INTEGER           IPIV(LDZ), JPIV(LDZ)
      DOUBLE PRECISION  RHS(LDZ), Z(LDZ,LDZ)
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMM, DGEMV, DGER, DGESC2,
     $                  DGETC2, DLASET, DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO = 0
      IERR = 0
      IF ( M.LE.0 ) THEN
         INFO = -1
      ELSE IF ( N.LE.0 ) THEN
         INFO = -2
      ELSE IF ( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF ( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF ( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF ( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF ( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -14
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB04OW', -INFO )
         RETURN
      END IF
C
C     Determine block structure of A.
C
      P = 0
      I = 1
   10 CONTINUE
      IF ( I.GT.M )
     $   GO TO 20
      P = P + 1
      IWORK(P) = I
      IF( I.EQ.M )
     $   GO TO 20
      IF ( A(I+1,I).NE.ZERO ) THEN
         I = I + 2
      ELSE
         I = I + 1
      END IF
      GO TO 10
   20 CONTINUE
      IWORK(P+1) = M + 1
C
C     Determine block structure of B.
C
      Q = P + 1
      J = 1
   30 CONTINUE
      IF ( J.GT.N )
     $   GO TO 40
      Q = Q + 1
      IWORK(Q) = J
      IF( J.EQ.N )
     $   GO TO 40
      IF ( B(J+1,J).NE.ZERO ) THEN
         J = J + 2
      ELSE
         J = J + 1
      END IF
      GO TO 30
   40 CONTINUE
      IWORK(Q+1) = N + 1
C
C     Solve (I, J) - subsystem
C       A(I,I) * R(I,J) - L(I,J) * B(J,J) = C(I,J)
C       D(I,I) * L(I,J) - R(I,J) * E(J,J) = F(I,J)
C     for I = P, P - 1, ..., 1; J = 1, 2, ..., Q.
C
      SCALE  = ONE
      SCALOC = ONE
      DO 120 J = P + 2, Q
         JS = IWORK(J)
         JSP1 = JS + 1
         JE = IWORK(J+1) - 1
         NB = JE - JS + 1
         DO 110 I = P, 1, -1
C
            IS = IWORK(I)
            ISP1 = IS + 1
            IE = IWORK(I+1) - 1
            MB = IE - IS + 1
            ZDIM = MB*NB*2
C
            IF ( ( MB.EQ.1 ).AND.( NB.EQ.1 ) ) THEN
C
C              Build a 2-by-2 system Z * x = RHS.
C
               Z(1,1) =  A(IS,IS)
               Z(2,1) = -E(JS,JS)
               Z(1,2) = -B(JS,JS)
               Z(2,2) =  D(IS,IS)
C
C              Set up right hand side(s).
C
               RHS(1) = C(IS,JS)
               RHS(2) = F(IS,JS)
C
C              Solve Z * x = RHS.
C
               CALL DGETC2( ZDIM, Z, LDZ, IPIV, JPIV, IERR )
               IF ( IERR.GT.0 )
     $            INFO = IERR
C
               CALL DGESC2( ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF ( SCALOC.NE.ONE ) THEN
                  DO 50 K = 1, N
                     CALL DSCAL( M, SCALOC, C(1,K), 1 )
                     CALL DSCAL( M, SCALOC, F(1,K), 1 )
   50             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
C
C              Unpack solution vector(s).
C
               C(IS,JS) = RHS(1)
               F(IS,JS) = RHS(2)
C
C              Substitute R(I,J) and L(I,J) into remaining equation.
C
               IF ( I.GT.1 ) THEN
                  CALL DAXPY( IS-1, -RHS(1), A(1,IS), 1, C(1,JS), 1 )
                  CALL DAXPY( IS-1, -RHS(2), D(1,IS), 1, F(1,JS), 1 )
               END IF
               IF ( J.LT.Q ) THEN
                  CALL DAXPY( N-JE, RHS(2), B(JS,JE+1), LDB, C(IS,JE+1),
     $                        LDC )
                  CALL DAXPY( N-JE, RHS(1), E(JS,JE+1), LDE, F(IS,JE+1),
     $                        LDF )
               END IF
C
            ELSE IF ( ( MB.EQ.1 ).AND.( NB.EQ.2 ) ) THEN
C
C              Build a 4-by-4 system Z * x = RHS.
C
               Z(1,1) =  A(IS,IS)
               Z(2,1) =  ZERO
               Z(3,1) = -E(JS,JS)
               Z(4,1) = -E(JS,JSP1)
C
               Z(1,2) =  ZERO
               Z(2,2) =  A(IS,IS)
               Z(3,2) =  ZERO
               Z(4,2) = -E(JSP1,JSP1)
C
               Z(1,3) = -B(JS,JS)
               Z(2,3) = -B(JS,JSP1)
               Z(3,3) =  D(IS,IS)
               Z(4,3) =  ZERO
C
               Z(1,4) = -B(JSP1,JS)
               Z(2,4) = -B(JSP1,JSP1)
               Z(3,4) =  ZERO
               Z(4,4) =  D(IS,IS)
C
C              Set up right hand side(s).
C
               RHS(1) = C(IS,JS)
               RHS(2) = C(IS,JSP1)
               RHS(3) = F(IS,JS)
               RHS(4) = F(IS,JSP1)
C
C              Solve Z * x = RHS.
C
               CALL DGETC2( ZDIM, Z, LDZ, IPIV, JPIV, IERR )
               IF ( IERR.GT.0 )
     $            INFO = IERR
C
               CALL DGESC2( ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF ( SCALOC.NE.ONE ) THEN
                  DO 60 K = 1, N
                     CALL DSCAL( M, SCALOC, C(1,K), 1 )
                     CALL DSCAL( M, SCALOC, F(1,K), 1 )
   60             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
C
C              Unpack solution vector(s).
C
               C(IS,JS)   = RHS(1)
               C(IS,JSP1) = RHS(2)
               F(IS,JS)   = RHS(3)
               F(IS,JSP1) = RHS(4)
C
C              Substitute R(I,J) and L(I,J) into remaining equation.
C
               IF ( I.GT.1 ) THEN
                  CALL DGER( IS-1, NB, -ONE, A(1,IS), 1, RHS(1), 1,
     $                       C(1,JS), LDC )
                  CALL DGER( IS-1, NB, -ONE, D(1,IS), 1, RHS(3), 1,
     $                       F(1,JS), LDF )
               END IF
               IF ( J.LT.Q ) THEN
                  CALL DAXPY( N-JE, RHS(3), B(JS,JE+1), LDB, C(IS,JE+1),
     $                        LDC )
                  CALL DAXPY( N-JE, RHS(1), E(JS,JE+1), LDE, F(IS,JE+1),
     $                        LDF )
                  CALL DAXPY( N-JE, RHS(4), B(JSP1,JE+1), LDB,
     $                        C(IS,JE+1), LDC )
                  CALL DAXPY( N-JE, RHS(2), E(JSP1,JE+1), LDE,
     $                        F(IS,JE+1), LDF )
               END IF
C
            ELSE IF( ( MB.EQ.2 ) .AND. ( NB.EQ.1 ) ) THEN
C
C              Build a 4-by-4 system Z * x = RHS.
C
               Z(1,1) =  A(IS,IS)
               Z(2,1) =  A(ISP1,IS)
               Z(3,1) = -E(JS,JS)
               Z(4,1) =  ZERO
C
               Z(1,2) =  A(IS,ISP1)
               Z(2,2) =  A(ISP1,ISP1)
               Z(3,2) =  ZERO
               Z(4,2) = -E(JS,JS)
C
               Z(1,3) = -B(JS,JS)
               Z(2,3) =  ZERO
               Z(3,3) =  D(IS,IS)
               Z(4,3) =  ZERO
C
               Z(1,4) =  ZERO
               Z(2,4) = -B(JS,JS)
               Z(3,4) =  D(IS,ISP1)
               Z(4,4) =  D(ISP1,ISP1)
C
C              Set up right hand side(s).
C
               RHS(1) = C(IS,JS)
               RHS(2) = C(ISP1,JS)
               RHS(3) = F(IS,JS)
               RHS(4) = F(ISP1,JS)
C
C              Solve Z * x = RHS.
C
               CALL DGETC2( ZDIM, Z, LDZ, IPIV, JPIV, IERR )
               IF ( IERR.GT.0 )
     $            INFO = IERR
C
               CALL DGESC2( ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF ( SCALOC.NE.ONE ) THEN
                  DO 70 K = 1, N
                     CALL DSCAL( M, SCALOC, C(1,K), 1 )
                     CALL DSCAL( M, SCALOC, F(1,K), 1 )
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
C
C              Unpack solution vector(s).
C
               C(IS,JS)   = RHS(1)
               C(ISP1,JS) = RHS(2)
               F(IS,JS)   = RHS(3)
               F(ISP1,JS) = RHS(4)
C
C              Substitute R(I,J) and L(I,J) into remaining equation.
C
               IF ( I.GT.1 ) THEN
                  CALL DGEMV( 'N', IS-1, MB, -ONE, A(1,IS), LDA, RHS(1),
     $                        1, ONE, C(1,JS), 1 )
                  CALL DGEMV( 'N', IS-1, MB, -ONE, D(1,IS), LDD, RHS(3),
     $                        1, ONE, F(1,JS), 1 )
               END IF
               IF ( J.LT.Q ) THEN
                  CALL DGER( MB, N-JE, ONE, RHS(3), 1, B(JS,JE+1), LDB,
     $                       C(IS,JE+1), LDC )
                  CALL DGER( MB, N-JE, ONE, RHS(1), 1, E(JS,JE+1), LDE,
     $                       F(IS,JE+1), LDF )
               END IF
C
            ELSE IF ( ( MB.EQ.2 ).AND.( NB.EQ.2 ) ) THEN
C
C              Build an 8-by-8 system Z * x = RHS.
C
               CALL DLASET( 'All', LDZ, LDZ, ZERO, ZERO, Z, LDZ )
C
               Z(1,1) =  A(IS,IS)
               Z(2,1) =  A(ISP1,IS)
               Z(5,1) = -E(JS,JS)
               Z(7,1) = -E(JS,JSP1)
C
               Z(1,2) =  A(IS,ISP1)
               Z(2,2) =  A(ISP1,ISP1)
               Z(6,2) = -E(JS,JS)
               Z(8,2) = -E(JS,JSP1)
C
               Z(3,3) =  A(IS,IS)
               Z(4,3) =  A(ISP1,IS)
               Z(7,3) = -E(JSP1,JSP1)
C
               Z(3,4) =  A(IS,ISP1)
               Z(4,4) =  A(ISP1,ISP1)
               Z(8,4) = -E(JSP1,JSP1)
C
               Z(1,5) = -B(JS,JS)
               Z(3,5) = -B(JS,JSP1)
               Z(5,5) =  D(IS,IS)
C
               Z(2,6) = -B(JS,JS)
               Z(4,6) = -B(JS,JSP1)
               Z(5,6) =  D(IS,ISP1)
               Z(6,6) =  D(ISP1,ISP1)
C
               Z(1,7) = -B(JSP1,JS)
               Z(3,7) = -B(JSP1,JSP1)
               Z(7,7) =  D(IS,IS)
C
               Z(2,8) = -B(JSP1,JS)
               Z(4,8) = -B(JSP1,JSP1)
C
               Z(7,8) = D(IS,ISP1)
               Z(8,8) = D(ISP1,ISP1)
C
C              Set up right hand side(s).
C
               K  = 1
               II = MB*NB + 1
               DO 80 JJ = 0, NB - 1
                  CALL DCOPY( MB, C(IS,JS+JJ), 1, RHS(K),  1 )
                  CALL DCOPY( MB, F(IS,JS+JJ), 1, RHS(II), 1 )
                  K  = K  + MB
                  II = II + MB
   80          CONTINUE
C
C              Solve Z * x = RHS.
C
               CALL DGETC2( ZDIM, Z, LDZ, IPIV, JPIV, IERR )
               IF ( IERR.GT.0 )
     $            INFO = IERR
C
               CALL DGESC2( ZDIM, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF ( SCALOC.NE.ONE ) THEN
                  DO 90 K = 1, N
                     CALL DSCAL( M, SCALOC, C(1,K), 1 )
                     CALL DSCAL( M, SCALOC, F(1,K), 1 )
   90             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
C
C              Unpack solution vector(s).
C
               K  = 1
               II = MB*NB + 1
               DO 100 JJ = 0, NB - 1
                  CALL DCOPY( MB, RHS(K),  1, C(IS,JS+JJ), 1 )
                  CALL DCOPY( MB, RHS(II), 1, F(IS,JS+JJ), 1 )
                  K  = K  + MB
                  II = II + MB
  100          CONTINUE
C
C              Substitute R(I,J) and L(I,J) into remaining equation.
C
               K = MB*NB + 1
               IF ( I.GT.1 ) THEN
                  CALL DGEMM( 'N', 'N', IS-1, NB, MB, -ONE, A(1,IS),
     $                        LDA, RHS(1), MB, ONE, C(1,JS), LDC )
                  CALL DGEMM( 'N', 'N', IS-1, NB, MB, -ONE, D(1,IS),
     $                        LDD, RHS(K), MB, ONE, F(1,JS), LDF )
               END IF
               IF ( J.LT.Q ) THEN
                  CALL DGEMM( 'N', 'N', MB, N-JE, NB, ONE, RHS(K), MB,
     $                        B(JS,JE+1), LDB, ONE, C(IS,JE+1), LDC )
                  CALL DGEMM( 'N', 'N', MB, N-JE, NB, ONE, RHS(1), MB,
     $                        E(JS,JE+1), LDE, ONE, F(IS,JE+1), LDF )
               END IF
C
            END IF
C
  110    CONTINUE
  120 CONTINUE
      RETURN
C *** Last line of SB04OW ***
      END
