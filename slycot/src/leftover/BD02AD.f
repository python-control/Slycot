      SUBROUTINE BD02AD( DEF, NR, DPAR, IPAR, VEC, N, M, P, E, LDE, A,
     1                   LDA, B, LDB, C, LDC, D, LDD, NOTE, DWORK,
     2                   LDWORK, INFO )
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
C     To generate benchmark examples for time-invariant,
C     discrete-time dynamical systems
C
C       E x_k+1 = A x_k + B u_k
C
C           y_k = C x_k + D u_k
C
C     E, A are real N-by-N matrices, B is N-by-M, C is P-by-N, and
C     D is P-by-M. In many examples, E is the identity matrix and D is
C     the zero matrix.
C
C     This routine is an implementation of the benchmark library
C     DTDSX (Version 1.0) described in [1].
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DEF     CHARACTER*1
C             Specifies the kind of values used as parameters when
C             generating parameter-dependent and scalable examples
C             (i.e., examples with NR(1) = 2, 3, or 4):
C             = 'D':  Default values defined in [1] are used;
C             = 'N':  Values set in DPAR and IPAR are used.
C             This parameter is not referenced if NR(1) = 1.
C             Note that the scaling parameter of examples with
C             NR(1) = 3 or 4 is considered as a regular parameter in
C             this context.
C
C     Input/Output Parameters
C
C     NR      (input) INTEGER array, dimension (2)
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
C     DPAR    (input/output) DOUBLE PRECISION array, dimension (7)
C             On entry, if DEF = 'N' and the desired example depends on
C             real parameters, then the array DPAR must contain the
C             values for these parameters.
C             For an explanation of the parameters see [1].
C             For Example 2.1, DPAR(1), ..., DPAR(3) define the
C             parameters 'tau', 'delta', 'K', respectively.
C             On exit, if DEF = 'D' and the desired example depends on
C             real parameters, then the array DPAR is overwritten by the
C             default values given in [1].
C
C     IPAR    (input/output) INTEGER array, dimension (1)
C             On entry, if DEF = 'N' and the desired example depends on
C             integer parameters, then the array IPAR must contain the
C             values for these parameters.
C             For an explanation of the parameters see [1].
C             For Example 3.1, IPAR(1) defines the parameter 'n'.
C             On exit, if DEF = 'D' and the desired example depends on
C             integer parameters, then the array IPAR is overwritten by
C             the default values given in [1].
C
C     VEC     (output) LOGICAL array, dimension (8)
C             Flag vector which displays the availabilty of the output
C             data:
C             VEC(1), ..., VEC(3) refer to N, M, and P, respectively,
C             and are always .TRUE..
C             VEC(4) is .TRUE. iff E is NOT the identity matrix.
C             VEC(5), ..., VEC(7) refer to A, B, and C, respectively,
C             and are always .TRUE..
C             VEC(8) is .TRUE. iff D is NOT the zero matrix.
C
C     N       (output) INTEGER
C             The actual state dimension, i.e., the order of the
C             matrices E and A.
C
C     M       (output) INTEGER
C             The number of columns in the matrices B and D.
C
C     P       (output) INTEGER
C             The number of rows in the matrices C and D.
C
C     E       (output) DOUBLE PRECISION array, dimension (LDE,N)
C             The leading N-by-N part of this array contains the
C             matrix E.
C             NOTE that this array is overwritten (by the identity
C             matrix), if VEC(4) = .FALSE..
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
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array contains the
C             matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= N.
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array contains the
C             matrix C.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= P.
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array contains the
C             matrix D.
C             NOTE that this array is overwritten (by the zero
C             matrix), if VEC(8) = .FALSE..
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= P.
C
C     NOTE    (output) CHARACTER*70
C             String containing short information about the chosen
C             example.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             NOTE that DWORK is not used in the current version
C             of BD02AD.
C
C     LDWORK  INTEGER
C             LDWORK >= 1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value; in particular, INFO = -3 or -4 indicates
C                   that at least one of the parameters in DPAR or
C                   IPAR, respectively, has an illegal value;
C             = 1:  data file can not be opened or has wrong format.
C
C     REFERENCES
C
C     [1]  Kressner, D., Mehrmann, V. and Penzl, T.
C          DTDSX - a Collection of Benchmark Examples for State-Space
C          Realizations of Discrete-Time Dynamical Systems.
C          SLICOT Working Note 1998-10. 1998.
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
C     discrete-time dynamical systems
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, THREE, FOUR, PI
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     1                    THREE = 3.0D0, FOUR = 4.0D0,
     2                    PI = .3141592653589793D1 )
C     .. Scalar Arguments ..
      CHARACTER         DEF
      CHARACTER*70      NOTE
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDWORK, M, N, P
C     .. Array Arguments ..
      LOGICAL           VEC(8)
      INTEGER           IPAR(*), NR(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), DPAR(*),
     1                  DWORK(*), E(LDE,*)
C     .. Local Scalars ..
      CHARACTER*12      DATAF
      INTEGER           I, J, STATUS
      DOUBLE PRECISION  TEMP
C     .. Local Arrays ..
      LOGICAL           VECDEF(8)
C     .. External Functions ..
C     . LAPACK .
      LOGICAL          LSAME
      EXTERNAL         LSAME
C     .. External Subroutines ..
C     . LAPACK .
      EXTERNAL         DLASET
C     .. Data Statements ..
C     . default values for availabities .
      DATA VECDEF /.TRUE., .TRUE., .TRUE., .FALSE.,
     1             .TRUE., .TRUE., .TRUE., .FALSE./
C
C     .. Executable Statements ..
C
      INFO = 0
      DO 10  I = 1, 8
        VEC(I) = VECDEF(I)
10    CONTINUE
C
      IF (NR(1) .EQ. 1) THEN
C
        IF (NR(2) .EQ. 1) THEN
          NOTE = 'Laub 1979, Ex. 2: uncontrollable-unobservable data'
          N = 2
          M = 1
          P = 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          A(1,1) = FOUR
          A(2,1) = -.45D1
          A(1,2) = THREE
          A(2,2) = -.35D1
          CALL DLASET('A', N, M, -ONE, ONE, B, LDB)
          C(1,1) = 3.0D0
          C(1,2) = 2.0D0
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'Laub 1979, Ex. 3'
          N = 2
          M = 2
          P = 2
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,1) = .9512D0
          A(2,2) = .9048D0
          B(1,1) = .4877D1
          B(1,2) = .4877D1
          B(2,1) = -.11895D1
          B(2,2) = .3569D1
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 3) THEN
          NOTE = 'Van Dooren 1981, Ex. II'
          N = 2
          M = 1
          P = 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          A(1,1) =  TWO
          A(2,1) =  ONE
          A(1,2) = -ONE
          A(2,2) = ZERO
          CALL DLASET('A', N, M, ZERO, ONE, B, LDB)
          CALL DLASET('A', P, N, ONE, ZERO, C, LDC)
          D(1,1) = ZERO
C
        ELSE IF (NR(2) .EQ. 4) THEN
          NOTE = 'Ionescu/Weiss 1992'
          N = 2
          M = 2
          P = 2
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,2) =  ONE
          A(2,2) = -ONE
          CALL DLASET('A', N, M, ZERO, ONE, B, LDB)
          B(2,1) =  TWO
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 5) THEN
          NOTE = 'Jonckheere 1981'
          N = 2
          M = 1
          P = 2
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,2) = ONE
          CALL DLASET('A', N, M, ONE, ZERO, B, LDB)
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 6) THEN
          NOTE = 'Ackerson/Fu 1970: satellite control problem'
          N = 4
          M = 2
          P = 4
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 7) THEN
          NOTE = 'Litkouhi 1983: system with slow and fast modes'
          N = 4
          M = 2
          P = 4
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 8) THEN
          NOTE = 'Lu/Lin 1993, Ex. 4.3'
          N = 4
          M = 4
          P = 4
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('U', P, N, ONE, ONE, C, LDC)
          C(1,3) = TWO
          C(1,4) = FOUR
          C(2,4) = TWO
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 9) THEN
          NOTE = 'Gajic/Shen 1993, Section 2.7.4: chemical plant'
          N = 5
          M = 2
          P = 5
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 10) THEN
          NOTE = 'Davison/Wang 1974'
          N = 6
          M = 2
          P = 2
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
          VEC(8) = .TRUE.
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,2) = ONE
          A(2,3) = ONE
          A(4,5) = ONE
          A(5,6) = ONE
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(3,1) = ONE
          B(6,2) = ONE
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,1) = ONE
          C(1,2) = ONE
          C(2,4) = ONE
          C(2,5) = -ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
          D(1,1) = ONE
          D(2,1) = ONE
C
        ELSE IF (NR(2) .EQ. 11) THEN
          NOTE = 'Patnaik et al. 1980: tubular ammonia reactor'
          N = 9
          M = 3
          P = 2
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,1) = ONE
          C(2,5) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 12) THEN
          NOTE = 'Smith 1969: two-stand cold rolling mill'
          N = 10
          M = 3
          P = 5
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
          VEC(8) = .TRUE.
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N, N, ZERO, ONE, A(2,1), LDA)
          A(1,10) = .112D0
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(1,1) =  .276D1
          B(1,2) = -.135D1
          B(1,3) = -.46D0
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,1)  = ONE
          C(2,10) =  .894D0
          C(3,10) = -.1693D2
          C(4,10) =  .7D-1
          C(5,10) =  .398D0
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = 'BD02112.dat')
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 110  I = 1, P
              READ (1, FMT = *, IOSTAT = STATUS) (D(I,J), J = 1, M)
              IF (STATUS .NE. 0)  INFO = 1
110         CONTINUE
          END IF
          CLOSE(1)
C
        ELSE
          INFO = -2
        END IF
C
        IF (((NR(2) .GE. 6) .AND. (NR(2) .LE. 9)) .OR.
     1       (NR(2) .EQ. 11))  THEN
C         .. loading data files
          WRITE (DATAF(1:11), '(A,I2.2,A)') 'BD021', NR(2), '.dat'
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = DATAF(1:11))
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 120  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (A(I,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
120         CONTINUE
            DO 130  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (B(I,J), J = 1, M)
              IF (STATUS .NE. 0)  INFO = 1
130         CONTINUE
          END IF
          CLOSE(1)
        END IF
C
      ELSE IF (NR(1) .EQ. 2) THEN
        IF (.NOT. (LSAME(DEF,'D') .OR. LSAME(DEF,'N'))) THEN
          INFO = -1
          RETURN
        END IF
C
        IF (NR(2) .EQ. 1) THEN
          NOTE = 'Pappas et al. 1980: process control of paper machine'
          IF (LSAME(DEF,'D')) THEN
            DPAR(1) = .1D9
            DPAR(2) = ONE
            DPAR(3) = ONE
          END IF
          IF (DPAR(1) .EQ. ZERO) INFO = -3
          N = 4
          M = 1
          P = 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          TEMP = DPAR(2) / DPAR(1)
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N-1, N-1, ZERO, ONE, A(2,1), LDA)
          A(1,1) = ONE - TEMP
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(1,1) = DPAR(3) * TEMP
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,4) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE
          INFO = -2
        END IF
C
      ELSE IF (NR(1) .EQ. 3) THEN
        IF (.NOT. (LSAME(DEF,'D') .OR. LSAME(DEF,'N'))) THEN
          INFO = -1
          RETURN
        END IF
C
        IF (NR(2) .EQ. 1) THEN
          NOTE = 'Pappas et al. 1980, Ex. 3'
          IF (LSAME(DEF,'D')) IPAR(1) = 100
          IF (IPAR(1) .LT. 2) INFO = -4
          N = IPAR(1)
          M = 1
          P = N
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N-1, N-1, ZERO, ONE, A(1,2), LDA)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(N,1) = ONE
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE
          INFO = -2
        END IF
C
      ELSE
        INFO = -2
      END IF
C
      RETURN
C *** Last Line of BD02AD ***
      END
