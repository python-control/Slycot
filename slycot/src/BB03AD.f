      SUBROUTINE BB03AD(DEF, NR, DPAR, IPAR, VEC, N, M, E, LDE, A, LDA,
     1                  Y, LDY, B, LDB, X, LDX, U, LDU, NOTE, DWORK,
     2                  LDWORK, INFO)
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
C     To generate benchmark examples of (generalized) continuous-time
C     Lyapunov equations
C
C        T           T
C       A  X  E  +  E  X A  =  Y .
C
C     In some examples, the right hand side has the form
C
C                T
C       Y  =  - B  B
C
C     and the solution can be represented as a product of Cholesky
C     factors
C
C              T
C       X  =  U  U .
C
C     E, A, Y, X, and U are real N-by-N matrices, and B is M-by-N. Note
C     that E can be the identity matrix. For some examples, B, X, or U
C     are not provided.
C
C     This routine is an implementation of the benchmark library
C     CTLEX (Version 1.0) described in [1].
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DEF     CHARACTER*1
C             Specifies the kind of values used as parameters when
C             generating parameter-dependent and scalable examples
C             (i.e., examples with NR(1) = 2, 3, or 4):
C             DEF = 'D' or 'd': Default values are used.
C             DEF = 'N' or 'n': Values set in DPAR and IPAR are used.
C             This parameter is not referenced if NR(1) = 1.
C             Note that the scaling parameter of examples with
C             NR(1) = 3 or 4 is considered as a regular parameter in
C             this context.
C
C     Input/Output Parameters
C
C     NR      (input) INTEGER array, dimension 2
C             Specifies the index of the desired example according
C             to [1].
C             NR(1) defines the group:
C                   1 : parameter-free problems of fixed size
C                   2 : parameter-dependent problems of fixed size
C                   3 : parameter-free problems of scalable size
C                   4 : parameter-dependent problems of scalable size
C             NR(2) defines the number of the benchmark example
C             within a certain group according to [1].
C
C     DPAR    (input/output) DOUBLE PRECISION array, dimension 2
C             On entry, if DEF = 'N' or 'n' and the desired example
C             depends on real parameters, then the array DPAR must
C             contain the values for these parameters.
C             For an explanation of the parameters see [1].
C             For Example 4.1, DPAR(1) and DPAR(2) define 'r' and 's',
C             respectively.
C             For Example 4.2, DPAR(1) and DPAR(2) define 'lambda' and
C             's', respectively.
C             For Examples 4.3 and 4.4, DPAR(1) defines the parameter
C             't'.
C             On exit, if DEF = 'D' or 'd' and the desired example
C             depends on real parameters, then the array DPAR is
C             overwritten by the default values given in [1].
C
C     IPAR    (input/output) INTEGER array of DIMENSION at least 1
C             On entry, if DEF = 'N' or 'n' and the desired example
C             depends on integer parameters, then the array IPAR must
C             contain the values for these parameters.
C             For an explanation of the parameters see [1].
C             For Examples 4.1, 4.2, and 4.3, IPAR(1) defines 'n'.
C             For Example 4.4, IPAR(1) defines 'q'.
C             On exit, if DEF = 'D' or 'd' and the desired example
C             depends on integer parameters, then the array IPAR is
C             overwritten by the default values given in [1].
C
C     VEC     (output) LOGICAL array, dimension 8
C             Flag vector which displays the availability of the output
C             data:
C             VEC(1) and VEC(2) refer to N and M, respectively, and are
C             always .TRUE.
C             VEC(3) is .TRUE. iff E is NOT the identity matrix.
C             VEC(4) and VEC(5) refer to A and Y, respectively, and are
C             always .TRUE.
C             VEC(6) is .TRUE. iff B is provided.
C             VEC(7) is .TRUE. iff the solution matrix X is provided.
C             VEC(8) is .TRUE. iff the Cholesky factor U is provided.
C
C     N       (output) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices E and A.
C
C     M       (output) INTEGER
C             The number of rows in the matrix B. If B is not provided
C             for the desired example, M = 0 is returned.
C
C     E       (output) DOUBLE PRECISION array, dimension (LDE,N)
C             The leading N-by-N part of this array contains the
C             matrix E.
C             NOTE that this array is overwritten (by the identity
C             matrix), if VEC(3) = .FALSE.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= N.
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the
C             matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     Y       (output) DOUBLE PRECISION array, dimension (LDY,N)
C             The leading N-by-N part of this array contains the
C             matrix Y.
C
C     LDY     INTEGER
C             The leading dimension of array Y.  LDY >= N.
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,N)
C             The leading M-by-N part of this array contains the
C             matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= M.
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,N)
C             The leading N-by-N part of this array contains the
C             matrix X.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= N.
C
C     U       (output) DOUBLE PRECISION array, dimension (LDU,N)
C             The leading N-by-N part of this array contains the
C             matrix U.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= N.
C
C     NOTE    (output) CHARACTER*70
C             String containing short information about the chosen
C             example.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             For Examples 4.1 and 4.2., LDWORK >= 2*IPAR(1) is
C             required.
C             For the other examples, no workspace is needed, i.e.,
C             LDWORK >= 1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value; in particular, INFO = -3 or -4 indicates
C                   that at least one of the parameters in DPAR or
C                   IPAR, respectively, has an illegal value.
C
C     REFERENCES
C
C     [1]  D. Kressner, V. Mehrmann, and T. Penzl.
C          CTLEX - a Collection of Benchmark Examples for Continuous-
C          Time Lyapunov Equations.
C          SLICOT Working Note 1999-6, 1999.
C
C     NUMERICAL ASPECTS
C
C     None
C
C     CONTRIBUTOR
C
C     D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz)
C
C     For questions concerning the collection or for the submission of
C     test examples, please contact Volker Mehrmann
C     (Email: volker.mehrmann@mathematik.tu-chemnitz.de).
C
C     REVISIONS
C
C     June 1999, V. Sima.
C
C     KEYWORDS
C
C     continuous-time Lyapunov equations
C
C     ********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, THREE, FOUR
      PARAMETER         (ZERO = .0D0, ONE = .1D1, TWO = .2D1,
     1                   THREE = .3D1, FOUR = .4D1)
C     .. Scalar Arguments ..
      CHARACTER         DEF
      CHARACTER*70      NOTE
      INTEGER           INFO, LDA, LDB, LDE, LDU, LDWORK, LDX, LDY, M, N
C     .. Array Arguments ..
      LOGICAL           VEC(8)
      INTEGER           IPAR(*), NR(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), DPAR(*), DWORK(LDWORK),
     1                  E(LDE,*), U(LDU,*), X(LDX,*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           I, J, K
      DOUBLE PRECISION  TEMP, TTM1, TTP1, TWOBYN
C     .. Local Arrays ..
      LOGICAL           VECDEF(8)
C     .. External Functions ..
C     . BLAS .
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     . LAPACK .
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
C     . BLAS .
      EXTERNAL          DGEMV, DGER, DAXPY
C     . LAPACK .
      EXTERNAL          DLASET
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MIN, MOD
C     .. Data Statements ..
C     . default values for availabilities .
      DATA VECDEF /.TRUE., .TRUE., .FALSE., .TRUE.,
     1             .TRUE., .FALSE., .FALSE., .FALSE./
C
C     .. Executable Statements ..
C
      INFO = 0
      DO 10  I = 1, 8
        VEC(I) = VECDEF(I)
  10  CONTINUE
C
      IF (NR(1) .EQ. 4) THEN
        IF (.NOT. (LSAME(DEF,'D') .OR. LSAME(DEF,'N'))) THEN
          INFO = -1
          RETURN
        END IF
C
        IF (NR(2) .EQ. 1) THEN
          NOTE = 'CTLEX: Example 4.1'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 10
            DPAR(1) = .15D1
            DPAR(2) = .15D1
          END IF
          IF ((DPAR(1) .LE. ONE) .OR. (DPAR(2) .LE. ONE)) INFO = -3
          IF (IPAR(1) .LT. 2) INFO = -4
          N = IPAR(1)
          M = 1
          IF (LDE .LT. N) INFO = -9
          IF (LDA .LT. N) INFO = -11
          IF (LDY .LT. N) INFO = -13
          IF (LDB .LT. M) INFO = -15
          IF (LDX .LT. N) INFO = -17
          IF (LDWORK .LT. N*2) INFO = -22
          IF (INFO .NE. 0) RETURN
C
          VEC(6) = .TRUE.
          VEC(7) = .TRUE.
          TWOBYN = TWO / DBLE( N )
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N, N, ZERO, ZERO, Y, LDY)
          CALL DLASET('A', M, N, ZERO, ZERO, B, LDB)
          CALL DLASET('A', N, N, ZERO, ZERO, X, LDX)
          DO 30  J = 1, N
            TEMP = DPAR(1) ** (J-1)
            A(J,J) = -TEMP
            DWORK(J) = ONE
            DO 20  I = 1, N
              X(I,J) = DBLE( I*J ) / (TEMP + DPAR(1)**(I-1))
  20        CONTINUE
  30      CONTINUE
C         H1 * A
          CALL DGEMV('T', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, A, LDA)
C         A * H1
          CALL DGEMV('N', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, A, LDA)
C         H1 * X
          CALL DGEMV('T', N,N, ONE, X, LDX, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, X, LDX)
C         X * H1
          CALL DGEMV('N', N,N, ONE, X, LDX, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, X, LDX)
C         S A INV(S), INV(S) X INV(S), B INV(S)
          DO 50  J = 1, N
            B(1,J) = DBLE( J-N-1 ) / (DPAR(2)**(J-1))
            DO 40  I = 1, N
              X(I,J) = X(I,J) / (DPAR(2)**(I+J-2))
              A(I,J) = A(I,J) * (DPAR(2)**(I-J))
  40        CONTINUE
            DWORK(J) = ONE - TWO * MOD(J,2)
  50      CONTINUE
C         H2 * A
          CALL DGEMV('T', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, A, LDA)
C         A * H2
          CALL DGEMV('N', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, A, LDA)
C         H2 * X
          CALL DGEMV('T', N,N, ONE, X, LDX, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, X, LDX)
C         X * H2
          CALL DGEMV('N', N,N, ONE, X, LDX, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, X, LDX)
C         B * H2
          CALL DAXPY(N, -TWOBYN * DDOT(N, B, LDB, DWORK, 1), DWORK, 1,
     1               B, LDB)
C         Y = -B' * B
          CALL DGER(N ,N, -ONE, B, LDB, B, LDB, Y, LDY)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'CTLEX: Example 4.2'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 10
            DPAR(1) = -.5D0
            DPAR(2) = .15D1
          END IF
          IF ((DPAR(1) .GE. ZERO) .OR. (DPAR(2) .LE. ONE)) INFO = -3
          IF (IPAR(1) .LT. 2) INFO = -4
          N = IPAR(1)
          M = 1
          IF (LDE .LT. N) INFO = -9
          IF (LDA .LT. N) INFO = -11
          IF (LDY .LT. N) INFO = -13
          IF (LDB .LT. M) INFO = -15
          IF (LDWORK .LT. N*2) INFO = -22
          IF (INFO .NE. 0) RETURN
C
          VEC(6) = .TRUE.
          TWOBYN = TWO / DBLE( N )
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, DPAR(1), A, LDA)
          CALL DLASET('A', N, N, ZERO, ZERO, Y, LDY)
          CALL DLASET('A', M, N, -TWOBYN, ONE - TWOBYN, B, LDB)
          DO 60  I = 1, N-1
            DWORK(I) = ONE
            A(I,I+1) = ONE
  60      CONTINUE
          DWORK(N) = ONE
C         H1 * A
          CALL DGEMV('T', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, A, LDA)
C         A * H1
          CALL DGEMV('N', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, A, LDA)
C         S A INV(S), B INV(S)
          DO 80  J = 1, N
            B(1,J) = B(1,J) / (DPAR(2)**(J-1))
            DO 70  I = 1, N
              A(I,J) = A(I,J) * (DPAR(2)**(I-J))
  70        CONTINUE
            DWORK(J) = ONE - TWO * MOD(J,2)
  80      CONTINUE
C         H2 * A
          CALL DGEMV('T', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK, 1, DWORK(N+1), 1, A, LDA)
C         A * H2
          CALL DGEMV('N', N,N, ONE, A, LDA, DWORK,1, ZERO, DWORK(N+1),1)
          CALL DGER(N, N, -TWOBYN, DWORK(N+1), 1, DWORK, 1, A, LDA)
C         B * H2
          CALL DAXPY(N, -TWOBYN * DDOT(N, B, LDB, DWORK, 1), DWORK, 1,
     1               B, LDB)
C         Y = -B' * B
          CALL DGER(N ,N, -ONE, B, LDB, B, LDB, Y, LDY)
C
        ELSE IF (NR(2) .EQ. 3) THEN
          NOTE = 'CTLEX: Example 4.3'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 10
            DPAR(1) = .1D2
          END IF
          IF (DPAR(1) .LT. ZERO) INFO = -3
          IF (IPAR(1) .LT. 2) INFO = -4
          N = IPAR(1)
          M = 0
          IF (LDE .LT. N) INFO = -9
          IF (LDA .LT. N) INFO = -11
          IF (LDY .LT. N) INFO = -13
          IF (LDX .LT. N) INFO = -17
          IF (INFO .NE. 0) RETURN
C
          VEC(3) = .TRUE.
          VEC(7) = .TRUE.
          TEMP = TWO ** (-DPAR(1))
          CALL DLASET('U', N, N, ZERO, ZERO, E, LDE)
          CALL DLASET('L', N, N, TEMP, ONE, E, LDE)
          CALL DLASET('L', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('U', N, N, ONE, ZERO, A, LDA)
          CALL DLASET('A', N, N, ONE, ONE, X, LDX)
          DO 90  I = 1, N
            A(I,I) = DBLE( I - 1 ) + TEMP
  90      CONTINUE
          Y(1,1) = TWO * TEMP + TWO * DBLE( N-1 ) * TEMP**2
          TTP1 = TWO * DBLE( N+1 ) * TEMP + TWO - TEMP**2
          TTM1 = TWO * DBLE( N-1 ) * TEMP + TWO - TEMP**2
          DO 100  I = 2, N
            Y(I,1) = Y(1,1) + DBLE( I-1 ) * TTM1
 100      CONTINUE
          DO 120  J = 2, N
            DO 110  I = 1, N
              Y(I,J) = Y(I,1) + DBLE( J-1 ) * (TTP1 - FOUR * I * TEMP)
 110        CONTINUE
 120     CONTINUE
C
        ELSE IF (NR(2) .EQ. 4) THEN
          NOTE = 'CTLEX: Example 4.4'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 10
            DPAR(1) = .15D1
          END IF
          IF (DPAR(1) .LT. ONE) INFO = -3
          IF (IPAR(1) .LT. 1) INFO = -4
          N = IPAR(1) * 3
          M = 1
          IF (LDE .LT. N) INFO = -9
          IF (LDA .LT. N) INFO = -11
          IF (LDY .LT. N) INFO = -13
          IF (LDB .LT. M) INFO = -15
          IF (INFO .NE. 0) RETURN
C
          VEC(3) = .TRUE.
          VEC(6) = .TRUE.
          CALL DLASET('A', N, N, ZERO, ZERO, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          DO 150  I = 1, IPAR(1)
            TEMP  = -DPAR(1)**I
            DO 140  J = 1, I - 1
              DO 130  K = 0, 2
                A(N - I*3+3, J*3-K) = TEMP
                A(N - I*3+2, J*3-K) = TWO * TEMP
 130          CONTINUE
 140        CONTINUE
            A(N - I*3+3, I*3-2) = TEMP
            A(N - I*3+2, I*3-2) = TWO * TEMP
            A(N - I*3+2, I*3-1) = TWO * TEMP
            A(N - I*3+2, I*3  ) = TEMP
            A(N - I*3+1, I*3  ) = TEMP
 150      CONTINUE
          DO 170  J = 1, N
            IF (J .GT. 1) CALL DAXPY(N, ONE, A(J-1,1), LDA, A(J,1), LDA)
            B(1, J) = DBLE( J )
            DO 160  I = 1, N
              E(I,N-J+1) = DBLE( MIN( I, J ) )
              Y(I,J) = -DBLE( I*J )
 160        CONTINUE
 170      CONTINUE
C
        ELSE
          INFO = -2
        END IF
      ELSE
        INFO = -2
      END IF
C
      RETURN
C *** Last Line of BB03AD ***
      END
