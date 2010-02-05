      SUBROUTINE BD01AD( DEF, NR, DPAR, IPAR, VEC, N, M, P, E, LDE, A,
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
C     continuous-time dynamical systems
C
C         .
C       E x(t) = A x(t) + B u(t)
C
C         y(t) = C x(t) + D u(t)
C
C     E, A are real N-by-N matrices, B is N-by-M, C is P-by-N, and
C     D is P-by-M. In many examples, E is the identity matrix and D is
C     the zero matrix.
C
C     This routine is an implementation of the benchmark library
C     CTDSX (Version 1.0) described in [1].
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
C             For Examples 2.1 and 2.2, DPAR(1) defines the parameter
C             'epsilon'.
C             For Example 2.4, DPAR(1), ..., DPAR(7) define 'b', 'mu',
C             'r', 'r_c', 'k_l', 'sigma', 'a', respectively.
C             For Example 2.7, DPAR(1) and DPAR(2) define 'mu' and 'nu',
C             respectively.
C             For Example 4.1, DPAR(1), ..., DPAR(7) define 'a', 'b',
C             'c', 'beta_1', 'beta_2', 'gamma_1', 'gamma_2',
C             respectively.
C             For Example 4.2, DPAR(1), ..., DPAR(3) define 'mu',
C             'delta', 'kappa', respectively.
C             On exit, if DEF = 'D' and the desired example depends on
C             real parameters, then the array DPAR is overwritten by the
C             default values given in [1].
C
C     IPAR    (input/output) INTEGER array, dimension (1)
C             On entry, if DEF = 'N' and the desired example depends on
C             integer parameters, then the array IPAR must contain the
C             values for these parameters.
C             For an explanation of the parameters see [1].
C             For Examples 2.3, 2.5, and 2.6, IPAR(1) defines the
C             parameter 's'.
C             For Example 3.1, IPAR(1) defines 'q'.
C             For Examples 3.2 and 3.3, IPAR(1) defines 'n'.
C             For Example 3.4, IPAR(1) defines 'l'.
C             For Example 4.1, IPAR(1) defines 'n'.
C             For Example 4.2, IPAR(1) defines 'l'.
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
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             For Example 3.4, LDWORK >= 4*IPAR(1) is required.
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
C                   IPAR, respectively, has an illegal value;
C             = 1:  data file can not be opened or has wrong format.
C
C
C     REFERENCES
C
C     [1]  Kressner, D., Mehrmann, V. and Penzl, T.
C          CTDSX - a Collection of Benchmark Examples for State-Space
C          Realizations of Continuous-Time Dynamical Systems.
C          SLICOT Working Note 1998-9. 1998.
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
C     continuous-time dynamical systems
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
      INTEGER           I, J, L, STATUS
      DOUBLE PRECISION  APPIND, B1, B2, C1, C2, TEMP, TTEMP
C     .. Local Arrays ..
      LOGICAL           VECDEF(8)
C     .. External Functions ..
C     . LAPACK .
      LOGICAL          LSAME
      EXTERNAL         LSAME
C     .. External Subroutines ..
C     . BLAS .
      EXTERNAL         DSCAL
C     . LAPACK .
      EXTERNAL         DLASET
C     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN, MOD
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
          NOTE = 'Laub 1979, Ex.1'
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
          B(1,1) = ZERO
          B(2,1) = ONE
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'Laub 1979, Ex.2: uncontrollable-unobservable data'
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
          A(1,2) = .3D1
          A(2,2) = -.35D1
          B(1,1) = ONE
          B(2,1) = -ONE
          C(1,1) = THREE
          C(1,2) = TWO
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 3) THEN
          NOTE = 'Beale/Shafai 1989: model of L-1011 aircraft'
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
        ELSE IF (NR(2) .EQ. 4) THEN
          NOTE = 'Bhattacharyya et al. 1983: binary distillation column'
          N = 8
          M = 2
          P = 8
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
        ELSE IF (NR(2) .EQ. 5) THEN
          NOTE = 'Patnaik et al. 1980: tubular ammonia reactor'
          N = 9
          M = 3
          P = 9
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
        ELSE IF (NR(2) .EQ. 6) THEN
          NOTE = 'Davison/Gesing 1978: J-100 jet engine'
          N = 30
          M = 3
          P = 5
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 7) THEN
          NOTE = 'Davison 1967: binary distillation column'
          N = 11
          M = 3
          P = 3
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(2,1) = ONE
          C(1,10) = ONE
          C(3,11) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)

        ELSE IF (NR(2) .EQ. 8) THEN
          NOTE = 'Chien/Ergin/Ling/Lee 1958: drum boiler'
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
          C(1,6) = ONE
          C(2,9) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 9) THEN
          NOTE = 'Ly, Gangsaas 1981: B-767 airplane'
          N = 55
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
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 10) THEN
          NOTE = 'control surface servo for an underwater vehicle'
          N = 8
          M = 2
          P = 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,7) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
        ELSE
          INFO = -2
        END IF
C
        IF ((NR(2) .GE. 3) .AND. (NR(2) .LE. 10)) THEN
C         .. loading data files
          WRITE (DATAF(1:11), '(A,I2.2,A)') 'BD011', NR(2), '.dat'
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = DATAF(1:11))
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 110  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (A(I,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
110         CONTINUE
            DO 120  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (B(I,J), J = 1, M)
              IF (STATUS .NE. 0)  INFO = 1
120         CONTINUE
            IF ((NR(2) .EQ. 6) .OR. (NR(2) .EQ. 9)) THEN
              DO 130  I = 1, P
                READ (1, FMT = *, IOSTAT = STATUS) (C(I,J), J = 1, N)
                IF (STATUS .NE. 0)  INFO = 1
130           CONTINUE
            END IF
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
          NOTE = 'Chow/Kokotovic 1976: magnetic tape control system'
          IF (LSAME(DEF,'D')) DPAR(1) = 1D-6
          IF (DPAR(1) .EQ. ZERO) INFO = -3
          N = 4
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
          A(1,2) =  .400D0
          A(2,3) =  .345D0
          A(3,2) = -.524D0/DPAR(1)
          A(3,3) = -.465D0/DPAR(1)
          A(3,4) =  .262D0/DPAR(1)
          A(4,4) = -ONE/DPAR(1)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(4,1) = ONE/DPAR(1)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,1) = ONE
          C(2,3) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'Arnold/Laub 1984'
          IF (LSAME(DEF,'D')) DPAR(1) = 1D-6
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
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, DPAR(1), A, LDA)
          A(1,1) = -DPAR(1)
          A(2,1) = -ONE
          A(1,2) =  ONE
          A(2,2) = -DPAR(1)
          A(4,3) = -ONE
          A(3,4) =  ONE
          CALL DLASET('A', N, M, ONE, ONE, B, LDB)
          CALL DLASET('A', P, N, ONE, ONE, C, LDC)
          D(1,1) = ZERO
C
        ELSE IF (NR(2) .EQ. 3) THEN
          NOTE = 'Vertical acceleration of a rigid guided missile'
          IF (LSAME(DEF,'D')) IPAR(1) = 1
          IF ((IPAR(1) .LT. 1) .OR. (IPAR(1) .GT. 10)) INFO = -4
          N = 3
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
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(2,1) = ONE
          A(3,3) = -.19D3
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(3,1) = .19D3
          D(1,1) = ZERO
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = 'BD01203.dat')
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 210  I = 1, IPAR(1)
              READ (1, FMT = *, IOSTAT = STATUS) (A(1,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
              READ (1, FMT = *, IOSTAT = STATUS) (A(2,J), J = 2, N)
              IF (STATUS .NE. 0)  INFO = 1
              READ (1, FMT = *, IOSTAT = STATUS) (C(1,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
210         CONTINUE
          END IF
          CLOSE(1)
C
        ELSE IF (NR(2) .EQ. 4) THEN
          NOTE = 'Senning 1980: hydraulic positioning system'
          IF (LSAME(DEF,'D')) THEN
            DPAR(1) = .14D5
            DPAR(2) = .1287D0
            DPAR(3) = .15D0
            DPAR(4) = .1D-1
            DPAR(5) = .2D-2
            DPAR(6) = .24D0
            DPAR(7) = .1075D2
          END IF
          IF (((DPAR(1) .LE. .9D4)    .OR. (DPAR(1) .GE. .16D5))   .OR.
     1        ((DPAR(2) .LE. .5D-1)   .OR. (DPAR(2) .GE. .3D0))    .OR.
     2        ((DPAR(3) .LE. .5D-1)   .OR. (DPAR(3) .GE. .5D1))    .OR.
     3        ((DPAR(4) .LE. ZERO)    .OR. (DPAR(4) .GE. .5D-1))   .OR.
     4        ((DPAR(5) .LE. .103D-3) .OR. (DPAR(5) .GE. .35D-2))  .OR.
     5        ((DPAR(6) .LE. .1D-2)   .OR. (DPAR(6) .GE. .15D2))   .OR.
     6        ((DPAR(7) .LE. .105D2)  .OR. (DPAR(7) .GE. .111D2))) THEN
            INFO = -3
          END IF
          N = 3
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
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,2) = ONE
          A(2,2) = -(DPAR(3) + FOUR*DPAR(4)/PI) / DPAR(2)
          A(2,3) = DPAR(7) / DPAR(2)
          A(3,2) = -FOUR * DPAR(7) * DPAR(1) / .874D3
          A(3,3) = -FOUR * DPAR(1) * (DPAR(6) + DPAR(5)) / .874D3
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(3,1) = -FOUR * DPAR(1) / .874D3
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          D(1,1) = 0
C
        ELSE IF (NR(2) .EQ. 5) THEN
          NOTE = 'Kwakernaak/Westdyk 1985: cascade of inverted pendula'
          IF (LSAME(DEF,'D')) IPAR(1) = 1
          IF ((IPAR(1) .LT. 1) .OR. (IPAR(1) .GT. 7)) INFO = -4
          IF (IPAR(1) .LE. 6) THEN
            M = IPAR(1)
          ELSE
            M = 10
          END IF
          N = 2 * M
          P = M
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          WRITE (DATAF(1:12), '(A,I1,A)') 'BD01205', IPAR(1), '.dat'
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = DATAF(1:12))
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 220  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (A(I,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
220         CONTINUE
            DO 230  I = 1, N
              READ (1, FMT = *, IOSTAT = STATUS) (B(I,J), J = 1, M)
              IF (STATUS .NE. 0)  INFO = 1
230         CONTINUE
            DO 240  I = 1, P
              READ (1, FMT = *, IOSTAT = STATUS) (C(I,J), J = 1, N)
              IF (STATUS .NE. 0)  INFO = 1
240         CONTINUE
          END IF
          CLOSE(1)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 6) THEN
          NOTE = 'Kallstrom/Astrom 1981: regulation of a ship heading'
          IF (LSAME(DEF,'D')) IPAR(1) = 1
          IF ((IPAR(1) .LT. 1) .OR. (IPAR(1) .GT. 5)) INFO = -4
          N = 3
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
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(3,2) = ONE
          B(3,1) = ZERO
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,3) = ONE
          D(1,1) = ZERO
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = 'BD01206.dat')
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 250  I = 1, IPAR(1)
              READ (1, FMT = *, IOSTAT = STATUS) (A(1,J), J = 1, 2)
              IF (STATUS .NE. 0)  INFO = 1
              READ (1, FMT = *, IOSTAT = STATUS) (A(2,J), J = 1, 2)
              IF (STATUS .NE. 0)  INFO = 1
              READ (1, FMT = *, IOSTAT = STATUS) (B(J,1), J = 1, 2)
              IF (STATUS .NE. 0)  INFO = 1
250         CONTINUE
          END IF
          CLOSE(1)
C
        ELSE IF (NR(2) .EQ. 7) THEN
          NOTE = 'Ackermann 1989: track-guided bus'
          IF (LSAME(DEF,'D')) THEN
            DPAR(1) = .15D2
            DPAR(2) = .1D2
          END IF
          IF ((DPAR(1) .LT. .995D1) .OR. (DPAR(1) .GT. .16D2)) INFO = -3
          IF ((DPAR(1) .LT. .1D1) .OR. (DPAR(1) .GT. .2D2)) INFO = -3
          N = 5
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
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          A(1,1) = -.668D3 / (DPAR(1)*DPAR(2))
          A(1,2) = -ONE + .1804D3 / (DPAR(1)*DPAR(2)**2)
          A(2,1) = .1804D3 / (.1086D2*DPAR(1))
          A(2,2) = -.44175452D4 / (.1086D2*DPAR(1)*DPAR(2))
          A(1,5) = 198 / (DPAR(1)*DPAR(2))
          A(2,5) = .72666D3 / (.1086D2*DPAR(1))
          A(3,1) = DPAR(2)
          A(3,4) = DPAR(2)
          A(4,2) = ONE
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(5,1) = ONE
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,3) = ONE
          C(1,4) = .612D1
          D(1,1) = 0
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
          NOTE = 'Laub 1979, Ex.4: string of high speed vehicles'
          IF (LSAME(DEF,'D')) IPAR(1) = 20
          IF (IPAR(1) .LT. 2) INFO = -4
          N = 2*IPAR(1) - 1
          M = IPAR(1)
          P = IPAR(1) - 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          DO 310  I = 1, N
            IF (MOD(I,2) .EQ. 1) THEN
              A(I,I)       = -ONE
              B(I,(I+1)/2) =  ONE
            ELSE
              A(I,I-1) =  ONE
              A(I,I+1) = -ONE
              C(I/2,I) =  ONE
            END IF
310       CONTINUE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'Hodel et al. 1996: heat flow in a thin rod'
          IF (LSAME(DEF,'D')) IPAR(1) = 100
          IF (IPAR(1) .LT. 1) INFO = -4
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
          TEMP = DBLE(N + 1)
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          CALL DLASET('A', N, N, ZERO, -TWO * TEMP, A, LDA)
          A(1,1) = -TEMP
          DO 320  I = 1, N - 1
            A(I,I+1) = TEMP
            A(I+1,I) = TEMP
320       CONTINUE
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(N,1) = TEMP
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 3) THEN
          NOTE = 'Laub 1979, Ex.6'
          IF (LSAME(DEF,'D')) IPAR(1) = 21
          IF (IPAR(1) .LT. 1) INFO = -4
          N = IPAR(1)
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
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          CALL DLASET('A', N-1, N-1, ZERO, ONE, A(1,2), LDA)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(N,1) = ONE
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          C(1,1) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 4) THEN
          NOTE = 'Lang/Penzl 1994: rotating axle'
          IF (LSAME(DEF,'D')) IPAR(1) = 211
          IF ((IPAR(1) .LT. 1) .OR. (IPAR(1) .GT. 211)) INFO = -4
          N = 2*IPAR(1) - 1
          M = IPAR(1)
          P = IPAR(1)
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (LDWORK .LT. M*4) INFO = -21
          IF (INFO .NE. 0) RETURN
C
          OPEN(1, IOSTAT = STATUS, STATUS = 'OLD', FILE = 'BD01304.dat')
          IF (STATUS .NE. 0) THEN
            INFO = 1
          ELSE
            DO 330  I = 1, M*4
              READ (1, FMT = *, IOSTAT = STATUS) DWORK(I)
              IF (STATUS .NE. 0)  INFO = 1
330         CONTINUE
          END IF
          CLOSE(1)
          IF (INFO .NE. 0) RETURN
          CALL DLASET('A', N, N, ZERO, ONE, E, LDE)
          E(1,1) = DWORK(1)
          DO 340  I = 2, M
            E(I,I-1) = DWORK((I-2) * 4 + 1)
            E(I,I) = -DWORK((I-1) * 4 + 1)
340       CONTINUE
          E(M,M) = -E(M,M)
          DO 350  I = M-1, 1, -1
            DO 345  J = I, M
             IF (I .EQ. 1) THEN
               E(J,I) = E(J,I) - E(J,I+1)
             ELSE
               E(J,I) = E(J,I+1) - E(J,I)
             END IF
345         CONTINUE
350       CONTINUE
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          DO 360  I = 2, M
            A(I-1,I) = DWORK((I-2) * 4 + 3)
            A(I,I) = -TWO * DWORK((I-2) * 4 + 3) - DWORK((I-1) * 4 + 2)
            A(I,1) = DWORK((I-1) * 4 + 2) - DWORK((I-2) * 4 + 2)
            A(I-1,M+I-1) = DWORK((I-1) * 4)
            A(I,M+I-1) = -TWO * DWORK((I-1) * 4)
            IF (I .LT. M) THEN
              A(I+1,I) = DWORK((I-2) * 4 + 3)
              DO 355  J = I+1, M
                A(J,I) = A(J,I) + DWORK((J-2) * 4 + 2)
     1                          - DWORK((J-1) * 4 + 2)
355           CONTINUE
              A(I+1,M+I-1) = DWORK((I-1) * 4)
            END IF
360       CONTINUE
          A(1,1) = -DWORK(2)
          A(1,2) = -DWORK(3)
          A(1,M+1) = -A(1,M+1)
          CALL DLASET('A', M-1, M-1, ZERO, ONE, A(M+1,2), LDA)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          CALL DLASET('A', P, N, ZERO, ZERO, C, LDC)
          DO 370 I = 2, M
            B(I,I) = -ONE
            B(I,I-1) = ONE
            C(I,I) = DWORK((I-2) * 4 + 3)
            C(I,M+I-1) = DWORK((I-1) * 4)
370       CONTINUE
          B(1,1) = ONE
          C(1,1) = ONE
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE
          INFO = -2
        END IF
C
      ELSE IF (NR(1) .EQ. 4) THEN
        IF (.NOT. (LSAME(DEF,'D') .OR. LSAME(DEF,'N'))) THEN
          INFO = -1
          RETURN
        END IF
C
        IF (NR(2) .EQ. 1) THEN
          NOTE = 'Rosen/Wang 1995: control of 1-dim. heat flow'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 100
            DPAR(1) = .1D-1
            DPAR(2) = ONE
            DPAR(3) = ONE
            DPAR(4) = .2D0
            DPAR(5) = .3D0
            DPAR(6) = .2D0
            DPAR(7) = .3D0
          END IF
          IF (IPAR(1) .LT. 2) INFO = -4
          N = IPAR(1)
          M = 1
          P = 1
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          VEC(4) = .TRUE.
          APPIND = DBLE(N + 1)
          TTEMP = -DPAR(1) * APPIND
          TEMP = 1 / (.6D1 * APPIND)
          CALL DLASET('A', N, N, ZERO, FOUR*TEMP, E, LDE)
          CALL DLASET('A', N, N, ZERO, TWO*TTEMP, A, LDA)
          DO 410  I = 1, N - 1
            A(I+1,I) = -TTEMP
            A(I,I+1) = -TTEMP
            E(I+1,I) = TEMP
            E(I,I+1) = TEMP
410       CONTINUE
          DO 420  I = 1, N
            B1 = MAX(DBLE(I-1)/APPIND, DPAR(4))
            B2 = MIN(DBLE(I+1)/APPIND, DPAR(5))
            C1 = MAX(DBLE(I-1)/APPIND, DPAR(6))
            C2 = MIN(DBLE(I+1)/APPIND, DPAR(7))
            IF (B1 .GE. B2) THEN
              B(I,1) = ZERO
            ELSE
              B(I,1) = B2 - B1
              TEMP   = MIN(B2, DBLE(I)/APPIND)
              IF (B1 .LT. TEMP) THEN
                B(I,1) = B(I,1) + APPIND*(TEMP**2 - B1**2)/TWO
                B(I,1) = B(I,1) + DBLE(I)*(B1 - TEMP)
              END IF
              TEMP = MAX(B1, DBLE(I)/APPIND)
              IF (TEMP .LT. B2) THEN
                B(I,1) = B(I,1) - APPIND*(B2**2 - TEMP**2)/TWO
                B(I,1) = B(I,1) - DBLE(I)*(TEMP - B2)
              END IF
            END IF
            IF (C1 .GE. C2) THEN
              C(1,I) = ZERO
            ELSE
              C(1,I) = C2 - C1
              TEMP   = MIN(C2, DBLE(I)/APPIND)
              IF (C1 .LT. TEMP) THEN
                C(1,I) = C(1,I) + APPIND*(TEMP**2 - C1**2)/TWO
                C(1,I) = C(1,I) + DBLE(I)*(C1 - TEMP)
              END IF
              TEMP = MAX(C1, DBLE(I)/APPIND)
              IF (TEMP .LT. C2) THEN
                C(1,I) = C(1,I) - APPIND*(C2**2 - TEMP**2)/TWO
                C(1,I) = C(1,I) - DBLE(I)*(TEMP - C2)
              END IF
            END IF
420       CONTINUE
          CALL DSCAL(N, DPAR(2), B(1,1), 1)
          CALL DSCAL(N, DPAR(3), C(1,1), LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          NOTE = 'Hench et al. 1995: coupled springs, dashpots, masses'
          IF (LSAME(DEF,'D')) THEN
            IPAR(1) = 30
            DPAR(1) = FOUR
            DPAR(2) = FOUR
            DPAR(3) = ONE
          END IF
          IF (IPAR(1) .LT. 2) INFO = -4
          L = IPAR(1)
          N = 2*L
          M = 2
          P = 2*L
          IF (LDE .LT. N) INFO = -10
          IF (LDA .LT. N) INFO = -12
          IF (LDB .LT. N) INFO = -14
          IF (LDC .LT. P) INFO = -16
          IF (LDD .LT. P) INFO = -18
          IF (INFO .NE. 0) RETURN
C
          VEC(4) = .TRUE.
          CALL DLASET('A', N, N, ZERO, DPAR(1), E, LDE)
          CALL DLASET('A', N, N, ZERO, ZERO, A, LDA)
          TEMP = -TWO * DPAR(3)
          DO 430  I = 1, L
            E(I,I) = ONE
            A(I,I+L) = ONE
            A(I+L,I+L) = -DPAR(2)
            IF (I .LT. L) THEN
              A(I+L,I+1) = DPAR(3)
              A(I+L+1,I) = DPAR(3)
              IF (I .GT. 1) THEN
                A(I+L,I) = TEMP
              END IF
            END IF
 430      CONTINUE
          A(L+1,1) = -DPAR(3)
          A(N,L) = -DPAR(3)
          CALL DLASET('A', N, M, ZERO, ZERO, B, LDB)
          B(L+1,1) = ONE
          B(N,2) = -ONE
          CALL DLASET('A', P, N, ZERO, ONE, C, LDC)
          CALL DLASET('A', P, M, ZERO, ZERO, D, LDD)
C
        ELSE
          INFO = -2
        END IF
      ELSE
        INFO = -2
      END IF
C
      RETURN
C *** Last Line of BD01AD ***
      END
