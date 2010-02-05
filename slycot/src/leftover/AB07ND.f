      SUBROUTINE AB07ND( N, M, A, LDA, B, LDB, C, LDC, D, LDD, RCOND,
     $                   IWORK, DWORK, LDWORK, INFO )
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
C     To compute the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D).
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs and outputs.  M >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state matrix A of the original system.
C             On exit, the leading N-by-N part of this array contains
C             the state matrix Ai of the inverse system.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the original system.
C             On exit, the leading N-by-M part of this array contains
C             the input matrix Bi of the inverse system.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading M-by-N part of this array must
C             contain the output matrix C of the original system.
C             On exit, the leading M-by-N part of this array contains
C             the output matrix Ci of the inverse system.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,M).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading M-by-M part of this array must
C             contain the feedthrough matrix D of the original system.
C             On exit, the leading M-by-M part of this array contains
C             the feedthrough matrix Di of the inverse system.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,M).
C
C     RCOND   (output) DOUBLE PRECISION
C             The estimated reciprocal condition number of the
C             feedthrough matrix D of the original system.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*M)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 or M+1, DWORK(1) returns the optimal
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,4*M).
C             For good performance, LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = i:  the matrix D is exactly singular; the (i,i) diagonal
C                   element is zero, i <= M; RCOND was set to zero;
C             = M+1:  the matrix D is numerically singular, i.e., RCOND
C                   is less than the relative machine precision, EPS
C                   (see LAPACK Library routine DLAMCH). The
C                   calculations have been completed, but the results
C                   could be very inaccurate.
C
C     METHOD
C
C     The matrices of the inverse system are computed with the formulas:
C                   -1              -1         -1           -1
C       Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  .
C
C     NUMERICAL ASPECTS
C
C     The accuracy depends mainly on the condition number of the matrix
C     D to be inverted. The estimated reciprocal condition number is
C     returned in RCOND.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000.
C     D. Sima, University of Bucharest, April 2000.
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
C     Based on the routine SYSINV, A. Varga, 1992.
C
C     REVISIONS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000.
C
C     KEYWORDS
C
C     Inverse system, state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   RCOND
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDWORK, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                   DWORK(*)
      INTEGER            IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION   DNORM
      INTEGER            BL, CHUNK, I, IERR, J, MAXWRK
      LOGICAL            BLAS3, BLOCK
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      INTEGER            ILAENV
      EXTERNAL           DLAMCH, DLANGE, ILAENV
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGECON, DGEMM, DGEMV, DGETRF, DGETRI,
     $                   DLACPY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDWORK.LT.MAX( 1, 4*M ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB07ND', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( M.EQ.0 ) THEN
         RCOND    = ONE
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Factorize D.
C
      CALL DGETRF( M, M, D, LDD, IWORK, INFO )
      IF ( INFO.NE.0 ) THEN
         RCOND = ZERO
         RETURN
      END IF
C
C     Compute the reciprocal condition number of the matrix D.
C     Workspace: need   4*M.
C     (Note: Comments in the code beginning "Workspace:" describe the
C      minimal amount of workspace needed at that point in the code,
C      as well as the preferred amount for good performance.
C      NB refers to the optimal block size for the immediately
C      following subroutine, as returned by ILAENV.)
C
      DNORM = DLANGE( '1-norm', M, M, D, LDD, DWORK )
      CALL DGECON( '1-norm', M, D, LDD, DNORM, RCOND, DWORK, IWORK(M+1),
     $             IERR )
      IF ( RCOND.LT.DLAMCH( 'Epsilon' ) )
     $   INFO = M + 1
C                   -1
C     Compute Di = D  .
C     Workspace: need   M;
C                prefer M*NB.
C
      MAXWRK = MAX( 4*M, M*ILAENV( 1, 'DGETRI', ' ', M, -1, -1, -1 ) )
      CALL DGETRI( M, D, LDD, IWORK, DWORK, LDWORK, IERR )
      IF ( N.GT.0 ) THEN
         CHUNK = LDWORK / M
         BLAS3 = CHUNK.GE.N .AND. M.GT.1
         BLOCK = MIN( CHUNK, M ).GT.1
C                          -1
C        Compute  Bi = -B*D  .
C
         IF ( BLAS3 ) THEN
C
C           Enough workspace for a fast BLAS 3 algorithm.
C
            CALL DLACPY( 'Full', N, M, B, LDB, DWORK, N )
            CALL DGEMM( 'NoTranspose', 'NoTranspose', N, M, M, -ONE,
     $                  DWORK, N, D, LDD, ZERO, B, LDB )
C
         ELSE IF( BLOCK ) THEN
C
C           Use as many rows of B as possible.
C
            DO 10 I = 1, N, CHUNK
               BL = MIN( N-I+1, CHUNK )
               CALL DLACPY( 'Full', BL, M, B(I,1), LDB, DWORK, BL )
               CALL DGEMM( 'NoTranspose', 'NoTranspose', BL, M, M, -ONE,
     $                     DWORK, BL, D, LDD, ZERO, B(I,1), LDB )
   10       CONTINUE
C
         ELSE
C
C           Use a BLAS 2 algorithm.
C
            DO 20 I = 1, N
               CALL DCOPY( M, B(I,1), LDB, DWORK, 1 )
               CALL DGEMV( 'Transpose', M, M, -ONE, D, LDD, DWORK, 1,
     $                     ZERO, B(I,1), LDB )
   20       CONTINUE
C
         END IF
C
C        Compute  Ai = A + Bi*C.
C
         CALL DGEMM( 'NoTranspose', 'NoTranspose', N, N, M, ONE, B, LDB,
     $               C, LDC, ONE, A, LDA )
C                        -1
C        Compute  C <-- D  *C.
C
         IF ( BLAS3 ) THEN
C
C           Enough workspace for a fast BLAS 3 algorithm.
C
            CALL DLACPY( 'Full', M, N, C, LDC, DWORK, M )
            CALL DGEMM( 'NoTranspose', 'NoTranspose', M, N, M, ONE,
     $                  D, LDD, DWORK, M, ZERO, C, LDC )
C
         ELSE IF( BLOCK ) THEN
C
C           Use as many columns of C as possible.
C
            DO 30 J = 1, N, CHUNK
               BL = MIN( N-J+1, CHUNK )
               CALL DLACPY( 'Full', M, BL, C(1,J), LDC, DWORK, M )
               CALL DGEMM( 'NoTranspose', 'NoTranspose', M, BL, M, ONE,
     $                     D, LDD, DWORK, M, ZERO, C(1,J), LDC )
   30       CONTINUE
C
         ELSE
C
C           Use a BLAS 2 algorithm.
C
            DO 40 J = 1, N
               CALL DCOPY( M, C(1,J), 1, DWORK, 1 )
               CALL DGEMV( 'NoTranspose', M, M, ONE, D, LDD, DWORK, 1,
     $                     ZERO, C(1,J), 1 )
   40       CONTINUE
C
         END IF
      END IF
C
C     Return optimal workspace in DWORK(1).
C
      DWORK(1) = DBLE( MAX( MAXWRK, N*M ) )
      RETURN
C
C *** Last line of AB07ND ***
      END
