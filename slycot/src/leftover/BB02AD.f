      SUBROUTINE BB02AD(DEF, NR, DPAR, IPAR, BPAR, CHPAR, VEC, N, M, P,
     1                  A, LDA, B, LDB, C, LDC, Q, LDQ, R, LDR, S, LDS,
     2                  X, LDX, DWORK, LDWORK, INFO)
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
C     To generate the benchmark examples for the numerical solution of
C     discrete-time algebraic Riccati equations (DAREs) of the form
C
C            T                T               T    -1  T       T
C     0  =  A X A  -  X  -  (A X B + S) (R + B X B)  (B X A + S )  +  Q
C
C     as presented in [1]. Here, A,Q,X are real N-by-N matrices, B,S are
C     N-by-M, and R is M-by-M. The matrices Q and R are symmetric and Q
C     may be given in factored form
C
C                   T
C     (I)    Q  =  C Q0 C .
C
C     Here, C is P-by-N and Q0 is P-by-P. If R is nonsingular and S = 0,
C     the DARE can be rewritten equivalently as
C
C                  T             -1
C     0  =  X  -  A X (I_n + G X)  A  -  Q,
C
C     where I_n is the N-by-N identity matrix and
C
C                   -1  T
C     (II)   G = B R   B .
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DEF     CHARACTER
C             This parameter specifies if the default parameters are
C             to be used or not.
C             = 'N' or 'n' : The parameters given in the input vectors
C                            xPAR (x = 'D', 'I', 'B', 'CH') are used.
C             = 'D' or 'd' : The default parameters for the example
C                            are used.
C             This parameter is not meaningful if NR(1) = 1.
C
C     Input/Output Parameters
C
C     NR      (input) INTEGER array, dimension (2)
C             This array determines the example for which DAREX returns
C             data. NR(1) is the group of examples.
C             NR(1) = 1 : parameter-free problems of fixed size.
C             NR(1) = 2 : parameter-dependent problems of fixed size.
C             NR(1) = 3 : parameter-free problems of scalable size.
C             NR(1) = 4 : parameter-dependent problems of scalable size.
C             NR(2) is the number of the example in group NR(1).
C             Let NEXi be the number of examples in group i. Currently,
C             NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1.
C             1 <= NR(1) <= 4;
C             0 <= NR(2) <= NEXi, where i = NR(1).
C
C     DPAR    (input/output) DOUBLE PRECISION array, dimension (4)
C             Double precision parameter vector. For explanation of the
C             parameters see [1].
C             DPAR(1) defines the parameter 'epsilon' for
C             examples NR = 2.2,2.3,2.4, the parameter 'tau'
C             for NR = 2.5, and the 1-by-1 matrix R for NR = 2.1,4.1.
C             For Example 2.5, DPAR(2) - DPAR(4) define in
C             consecutive order 'D', 'K', and 'r'.
C             NOTE that DPAR is overwritten with default values
C             if DEF = 'D' or 'd'.
C
C     IPAR    (input/output) INTEGER array, dimension (3)
C             On input, IPAR(1) determines the actual state dimension,
C             i.e., the order of the matrix A as follows:
C             NR(1) = 1, NR(1) = 2   : IPAR(1) is ignored.
C             NR = NR(1).NR(2) = 4.1 : IPAR(1) determines the order of
C                                      the output matrix A.
C             NOTE that IPAR(1) is overwritten for Examples 1.1-2.3. For
C             the other examples, IPAR(1) is overwritten if the default
C             parameters are to be used.
C             On output, IPAR(1) contains the order of the matrix A.
C
C             On input, IPAR(2) is the number of colums in the matrix B
C             and the order of the matrix R (in control problems, the
C             number of inputs of the system). Currently, IPAR(2) is
C             fixed for all examples and thus is not referenced on
C             input.
C             On output, IPAR(2) is the number of columns of the
C             matrix B from (I).
C
C             On input, IPAR(3) is the number of rows in the matrix C
C             (in control problems, the number of outputs of the
C             system). Currently, IPAR(3) is fixed for all examples
C             and thus is not referenced on input.
C             On output, IPAR(3) is the number of rows of the matrix C
C             from (I).
C
C             NOTE that IPAR(2) and IPAR(3) are overwritten and
C             IPAR(2) <= IPAR(1) and IPAR(3) <= IPAR(1) for all
C             examples.
C
C     BPAR    (input) LOGICAL array, dimension (7)
C             This array defines the form of the output of the examples
C             and the storage mode of the matrices Q, G or R.
C             BPAR(1) = .TRUE.  : Q is returned.
C             BPAR(1) = .FALSE. : Q is returned in factored form, i.e.,
C                                 Q0 and C from (I) are returned.
C             BPAR(2) = .TRUE.  : The matrix returned in array Q (i.e.,
C                                 Q if BPAR(1) = .TRUE. and Q0 if
C                                 BPAR(1) = .FALSE.) is stored as full
C                                 matrix.
C             BPAR(2) = .FALSE. : The matrix returned in array Q is
C                                 provided in packed storage mode.
C             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix
C                                 returned in array Q is stored in upper
C                                 packed mode, i.e., the upper triangle
C                                 of a symmetric n-by-n matrix is stored
C                                 by columns, e.g., the matrix entry
C                                 Q(i,j) is stored in the array entry
C                                 Q(i+j*(j-1)/2) for i <= j.
C                                 Otherwise, this entry is ignored.
C             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix
C                                 returned in array Q is stored in lower
C                                 packed mode, i.e., the lower triangle
C                                 of a symmetric n-by-n matrix is stored
C                                 by columns, e.g., the matrix entry
C                                 Q(i,j) is stored in the array entry
C                                 Q(i+(2*n-j)*(j-1)/2) for j <= i.
C                                 Otherwise, this entry is ignored.
C             BPAR(4) = .TRUE.  : The product G in (II) is returned.
C             BPAR(4) = .FALSE. : G is returned in factored form, i.e.,
C                                 B and R from (II) are returned.
C             BPAR(5) = .TRUE.  : The matrix returned in array R (i.e.,
C                                 G if BPAR(4) = .TRUE. and R if
C                                 BPAR(4) = .FALSE.) is stored as full
C                                 matrix.
C             BPAR(5) = .FALSE. : The matrix returned in array R is
C                                 provided in packed storage mode.
C             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix
C                                 returned in array R is stored in upper
C                                 packed mode (see above).
C                                 Otherwise, this entry is ignored.
C             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix
C                                 returned in array R is stored in lower
C                                 packed mode (see above).
C                                 Otherwise, this entry is ignored.
C             BPAR(7) = .TRUE.  : The coefficient matrix S of the DARE
C                                 is returned in array S.
C             BPAR(7) = .FALSE. : The coefficient matrix S of the DARE
C                                 is not returned.
C             NOTE that there are no default values for BPAR.  If all
C             entries are declared to be .TRUE., then matrices Q, G or R
C             are returned in conventional storage mode, i.e., as
C             N-by-N or M-by-M arrays where the array element Z(I,J)
C             contains the matrix entry Z_{i,j}.
C
C     CHPAR   (output) CHARACTER*255
C             On output, this string contains short information about
C             the chosen example.
C
C     VEC     (output) LOGICAL array, dimension (10)
C             Flag vector which displays the availability of the output
C             data:
C             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and
C             are always .TRUE.
C             VEC(4) refers to A and is always .TRUE.
C             VEC(5) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors B
C             and R from (II) are returned.
C             VEC(6) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors C
C             and Q0 from (I) are returned.
C             VEC(7) refers to Q and is always .TRUE.
C             VEC(8) refers to R and is always .TRUE.
C             VEC(9) is .TRUE. if BPAR(7) = .TRUE., i.e., the matrix S
C             is returned.
C             VEC(10) refers to X and is .TRUE. if the exact solution
C             matrix is available.
C             NOTE that VEC(i) = .FALSE. for i = 1 to 10 if on exit
C             INFO .NE. 0.
C
C     N       (output) INTEGER
C             The order of the matrices A, X, G if BPAR(4) = .TRUE., and
C             Q if BPAR(1) = .TRUE.
C
C     M       (output) INTEGER
C             The number of columns in the matrix B (or the dimension of
C             the control input space of the underlying dynamical
C             system).
C
C     P       (output) INTEGER
C             The number of rows in the matrix C (or the dimension of
C             the output space of the underlying dynamical system).
C
C     A       (output) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array contains the
C             coefficient matrix A of the DARE.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
C             If (BPAR(4) = .FALSE.), then the leading N-by-M part
C             of this array contains the coefficient matrix B of
C             the DARE.  Otherwise, B is used as workspace.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= N.
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             If (BPAR(1) = .FALSE.), then the leading P-by-N part
C             of this array contains the matrix C of the factored
C             form (I) of Q.  Otherwise, C is used as workspace.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= P.
C
C     Q       (output) DOUBLE PRECISION array, dimension (NQ)
C             If (BPAR(1) = .TRUE.) and (BPAR(2) = .TRUE.), then
C             NQ = LDQ*N.
C             IF (BPAR(1) = .TRUE.) and (BPAR(2) = .FALSE.), then
C             NQ = N*(N+1)/2.
C             If (BPAR(1) = .FALSE.) and (BPAR(2) = .TRUE.), then
C             NQ = LDQ*P.
C             IF (BPAR(1) = .FALSE.) and (BPAR(2) = .FALSE.), then
C             NQ = P*(P+1)/2.
C             The symmetric matrix contained in array Q is stored
C             according to BPAR(2) and BPAR(3).
C
C     LDQ     INTEGER
C             If conventional storage mode is used for Q, i.e.,
C             BPAR(2) = .TRUE., then Q is stored like a 2-dimensional
C             array with leading dimension LDQ. If packed symmetric
C             storage mode is used, then LDQ is irrelevant.
C             LDQ >= N if BPAR(1) = .TRUE.;
C             LDQ >= P if BPAR(1) = .FALSE..
C
C     R       (output) DOUBLE PRECISION array, dimension (MR)
C             If (BPAR(4) = .TRUE.) and (BPAR(5) = .TRUE.), then
C             MR = LDR*N.
C             IF (BPAR(4) = .TRUE.) and (BPAR(5) = .FALSE.), then
C             MR = N*(N+1)/2.
C             If (BPAR(4) = .FALSE.) and (BPAR(5) = .TRUE.), then
C             MR = LDR*M.
C             IF (BPAR(4) = .FALSE.) and (BPAR(5) = .FALSE.), then
C             MR = M*(M+1)/2.
C             The symmetric matrix contained in array R is stored
C             according to BPAR(5) and BPAR(6).
C
C     LDR     INTEGER
C             If conventional storage mode is used for R, i.e.,
C             BPAR(5) = .TRUE., then R is stored like a 2-dimensional
C             array with leading dimension LDR. If packed symmetric
C             storage mode is used, then LDR is irrelevant.
C             LDR >= N  if BPAR(4) =  .TRUE.;
C             LDR >= M  if BPAR(4) = .FALSE..
C
C     S       (output) DOUBLE PRECISION array, dimension (LDS,M)
C             If (BPAR(7) = .TRUE.), then the leading N-by-M part of
C             this array contains the coefficient matrix S of the DARE.
C
C     LDS     INTEGER
C             The leading dimension of array S.  LDS >= 1, and
C             LDS >= N if BPAR(7) = .TRUE..
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,NX)
C             If an exact solution is available (NR = 1.1,1.3,1.4,2.1,
C             2.3,2.4,2.5,4.1), then NX = N and the leading N-by-N part
C             of this array contains the solution matrix X.
C             Otherwise, X is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= 1, and
C             LDX >= N if an exact solution is available.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= N*N.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0 : successful exit;
C             < 0 : if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1 : data file could not be opened or had wrong format;
C             = 2 : division by zero;
C             = 3 : G can not be computed as in (II) due to a singular R
C                   matrix. This error can only occur if
C                   BPAR(4) = .TRUE..
C
C     REFERENCES
C
C     [1] Abels, J. and Benner, P.
C         DAREX - A Collection of Benchmark Examples for Discrete-Time
C         Algebraic Riccati Equations (Version 2.0).
C         SLICOT Working Note 1999-16, November 1999. Available from
C         http://www.win.tue.nl/niconet/NIC2/reports.html.
C
C     This is an updated and extended version of
C
C     [2] Benner, P., Laub, A.J., and Mehrmann, V.
C         A Collection of Benchmark Examples for the Numerical Solution
C         of Algebraic Riccati Equations II: Discrete-Time Case.
C         Technical Report SPC 95_23, Fak. f. Mathematik,
C         TU Chemnitz-Zwickau (Germany), December 1995.
C
C     FURTHER COMMENTS
C
C     Some benchmark examples read data from the data files provided
C     with the collection.
C
C     CONTRIBUTOR
C
C     Peter Benner (Universitaet Bremen), November 25, 1999.
C
C     For questions concerning the collection or for the submission of
C     test examples, please send e-mail to benner@math.uni-bremen.de.
C
C     REVISIONS
C
C     1999, December 23 (V. Sima).
C
C     KEYWORDS
C
C     Discrete-time algebraic Riccati equation.
C
C     ******************************************************************
C
C     .. Parameters ..
C     . # of examples available , # of examples with fixed size. .
      INTEGER          NEX1, NEX2, NEX3, NEX4, NMAX
      PARAMETER        ( NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1 )
      PARAMETER        ( NMAX = 13 )
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, FIVE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     1                   THREE = 3.0D0, FOUR = 4.0D0, FIVE = 5.0D0 )
C
C     .. Scalar Arguments ..
      INTEGER          INFO, LDA, LDB, LDC, LDQ, LDR, LDS, LDWORK, LDX,
     $                 M, N, P
      CHARACTER        DEF
C
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), DPAR(*), DWORK(*),
     1                 Q(*), R(*), S(LDS,*), X(LDX,*)
      INTEGER          IPAR(3), NR(2)
      CHARACTER        CHPAR*255
      LOGICAL          BPAR(7), VEC(10)
C
C     .. Local Scalars ..
      INTEGER          I, IOS, ISYMM, J, MSYMM, NSYMM, PSYMM, QDIMM,
     1                 RDIMM
      DOUBLE PRECISION ALPHA, BETA, TEMP
C
C     ..Local Arrays ..
      INTEGER          MDEF(2,NMAX), NDEF(4,NMAX), NEX(4), PDEF(2,NMAX)
      CHARACTER        IDENT*4
      CHARACTER*255    NOTES(4,NMAX)
C
C     .. External Functions ..
C     . LAPACK .
      LOGICAL          LSAME
      EXTERNAL         LSAME
C
C     .. External Subroutines ..
C     . BLAS .
      EXTERNAL         DCOPY, DGEMV, DSPMV, DSPR, DSYMM, DSYRK
C     . LAPACK .
      EXTERNAL         DLASET, DPPTRF, DPPTRI, DRSCL, XERBLA
C     . SLICOT .
      EXTERNAL         MA02DD, MA02ED
C
C     .. Intrinsic Functions ..
      INTRINSIC        SQRT
C
C     .. Data Statements ..
C     . default values for dimensions .
      DATA NEX /NEX1, NEX2, NEX3, NEX4/
      DATA (NDEF(1,I), I = 1, NEX1) /2, 2, 2, 3, 4, 4, 4, 5, 6, 9,
     1                               11, 13, 26/
      DATA (NDEF(2,I), I = 1, NEX2) /2, 2, 2, 3, 4/
      DATA (NDEF(4,I), I = 1, NEX4) /100/
      DATA (MDEF(1,I), I = 1, NEX1) /1, 2, 1, 2, 2, 2, 4, 2, 2, 3,
     1                               2, 2, 6/
      DATA (MDEF(2,I), I = 1, NEX2) /1, 2, 1, 3, 1/
      DATA (PDEF(1,I), I = 1, NEX1) /1, 2, 2, 3, 4, 4, 4, 5, 2, 2,
     1                               4, 4, 12/
      DATA (PDEF(2,I), I = 1, NEX2) /2, 2, 2, 3, 1/
C     . comments on examples .
      DATA (NOTES(1,I), I = 1, 10) /
     1'Van Dooren 1981, Ex. II: singular R matrix', 'Ionescu/Weiss 1992
     2: singular R matrix, nonzero S matrix', 'Jonckheere 1981: (A,B) co
     3ntrollable, no solution X <= 0', 'Sun 1998: R singular, Q non-defi
     4nite', 'Ackerson/Fu 1970 : satellite control problem', 'Litkouhi 1
     5983 : system with slow and fast modes', 'Lu/Lin 1993, Ex. 4.3', 'G
     6ajic/Shen 1993, Section 2.7.4: chemical plant', 'Davison/Wang 1974
     7: nonzero S matrix', 'Patnaik et al. 1980: tubular ammonia reactor
     8'/
      DATA (NOTES(1,I), I = 11, NEX1) /
     1'Sima 1996, Sec. 1.2.2: paper machine model error integrators', 'S
     2ima 1996, Ex. 2.6: paper machine model with with disturbances', 'P
     3ower plant model, Katayama et al., 1985'/
      DATA (NOTES(2,I), I = 1, NEX2) /
     1'Laub 1979, Ex. 2: uncontrollable-unobservable data', 'Laub 1979,
     2Ex. 3: increasingly ill-conditioned R-matrix', 'increasingly bad s
     3caled system as eps -> oo','Petkov et. al. 1989 : increasingly bad
     4 scaling as eps -> oo', 'Pappas et al. 1980: process control of pa
     5per machine'/
      DATA (NOTES(4,I), I = 1, NEX4) /'Pappas et al. 1980, Ex. 3'/
C
C     .. Executable Statements ..
C
      INFO = 0
      DO 1 I = 1, 10
        VEC(I) = .FALSE.
    1 CONTINUE
C
      IF (NR(1) .GE. 3) THEN
        IF (LSAME(DEF, 'D'))  IPAR(1) = NDEF(NR(1),NR(2))
        IPAR(2) = 1
        IPAR(3) = IPAR(1)
      ELSE
        IPAR(1) = NDEF(NR(1),NR(2))
        IPAR(2) = MDEF(NR(1),NR(2))
        IPAR(3) = PDEF(NR(1),NR(2))
      END IF
C
      IF ((NR(1) .GE. 2) .AND. .NOT. ((LSAME(DEF,'D')) .OR.
     $                                (LSAME(DEF,'N')))) THEN
        INFO = -1
      ELSE IF ((NR(1) .LT. 1) .OR. (NR(1) .GT. 4) .OR. (NR(2) .LT. 0)
     1    .OR. (NR(2) .GT. NEX(NR(1)))) THEN
        INFO = -2
      ELSE IF (IPAR(1) .LT. 1) THEN
        INFO = -4
      ELSE IF (IPAR(1) .GT. LDA) THEN
        INFO = -12
      ELSE IF (IPAR(1) .GT. LDB) THEN
        INFO = -14
      ELSE IF (IPAR(3) .GT. LDC) THEN
        INFO = -16
      ELSE IF (BPAR(2) .AND. (((.NOT. BPAR(1)) .AND.
     1       (IPAR(3) .GT. LDQ)) .OR. (BPAR(1) .AND.
     2       (IPAR(1) .GT. LDQ)))) THEN
        INFO = -18
      ELSE IF (BPAR(5) .AND. ((BPAR(4) .AND. (IPAR(1) .GT. LDR)) .OR.
     1         ((.NOT. BPAR(4)) .AND. (IPAR(2) .GT. LDR)))) THEN
        INFO = -20
      ELSE IF (LDS .LT. 1 .OR. (BPAR(7) .AND. (IPAR(1) .GT. LDS))) THEN
        INFO = -22
      ELSE IF (LDX .LT. 1) THEN
        INFO = -24
      ELSE IF (((NR(1) .EQ. 1) .AND. ((NR(2) .EQ. 1)   .OR.
     1          (NR(2) .EQ. 3) .OR.   (NR(2) .EQ. 4))) .OR.
     2         ((NR(1) .EQ. 2) .AND. ((NR(2). EQ. 1)   .OR.
     3          (NR(2) .GE. 3))) .OR. (NR(1) .EQ. 4)) THEN
C    .. solution X available ..
        IF (IPAR(1) .GT. LDX) THEN
          INFO = -24
        ELSE
          CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ZERO, X, LDX)
        END IF
      ELSE IF (LDWORK .LT. N*N) THEN
        INFO = -26
      END IF
      IF (INFO .NE. 0)  THEN
        CALL XERBLA( 'BB02AD', -INFO )
        RETURN
      END IF
C
      NSYMM = (IPAR(1)*(IPAR(1)+1))/2
      MSYMM = (IPAR(2)*(IPAR(2)+1))/2
      PSYMM = (IPAR(3)*(IPAR(3)+1))/2
C
      CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ZERO, A, LDA)
      CALL DLASET('A', IPAR(1), IPAR(2), ZERO, ZERO, B, LDB)
      CALL DLASET('A', IPAR(3), IPAR(1), ZERO, ZERO, C, LDC)
      CALL DLASET('L', PSYMM, 1, ZERO, ZERO, Q, 1)
      CALL DLASET('L', MSYMM, 1, ZERO, ZERO, R, 1)
      IF (BPAR(7))  CALL DLASET('A', IPAR(1), IPAR(2), ZERO, ZERO,
     1                          S, LDS)
C
      IF(NR(1) .EQ. 1) THEN
C
        IF (NR(2) .EQ. 1) THEN
          A(1,1) =  TWO
          A(2,1) =  ONE
          A(1,2) = -ONE
          B(1,1) =  ONE
          Q(1)   = ONE
          C(1,2) = ONE
          R(1)   = ZERO
          CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ONE, X, LDX)
          IDENT  = '0000'
C
        ELSE IF (NR(2) .EQ. 2) THEN
          A(1,2) =  ONE
          A(2,2) = -ONE
          B(1,1) =  ONE
          B(2,1) =  TWO
          B(2,2) =  ONE
          R(1)   = 9.0D0
          R(2)   = THREE
          R(3)   = ONE
          CALL DLASET('A', PSYMM, 1, -FOUR, -FOUR, Q, PSYMM)
          Q(3)   = 7.0D0
          CALL DRSCL(MSYMM, 11.0D0, Q, 1)
          IF (BPAR(7)) THEN
            S(1,1) =  THREE
            S(2,1) = -ONE
            S(1,2) =  ONE
            S(2,2) =  7.0D0
          END IF
          IDENT  = '0100'
C
        ELSE IF (NR(2) .EQ. 3) THEN
          A(1,2) = ONE
          B(2,1) = ONE
          Q(1)   = ONE
          Q(2)   = TWO
          Q(3)   = FOUR
          X(1,1) = ONE
          X(2,1) = TWO
          X(1,2) = TWO
          X(2,2) = TWO + SQRT(FIVE)
          IDENT  = '0101'
C
        ELSE IF (NR(2) .EQ. 4) THEN
          A(1,2) = .1000D+00
          A(2,3) = .0100D+00
          B(1,1) = ONE
          B(3,2) = ONE
          R(3) = ONE
          Q(1) = .1D+06
          Q(4) = .1D+04
          Q(6) = -.1D+02
          X(1,1) = .1D+06
          X(2,2) = .1D+04
          IDENT = '0100'
C
        ELSE IF (((NR(2) .GE. 5) .AND. (NR(2) .LE. 8)) .OR.
     1           (NR(2) .EQ. 10) .OR. (NR(2) .EQ. 11) .OR.
     2           (NR(2) .EQ. 13)) THEN
          IF (NR(2) .LT. 10) THEN
            WRITE (CHPAR(1:11), '(A,I1,A,I1,A)')
     1                              'BB02', NR(1), '0', NR(2), '.dat'
            OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:11))
          ELSE
            WRITE (CHPAR(1:11), '(A,I1,I2,A)')
     1                                   'BB02', NR(1), NR(2), '.dat'
            OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:11))
          END IF
          IF (IOS .NE. 0) THEN
            INFO = 1
          ELSE
            IF (.NOT. (NR(2) .EQ. 13)) THEN
              DO 10  I = 1, IPAR(1)
                READ (1, FMT = *, IOSTAT = IOS) (A(I,J), J = 1, IPAR(1))
                IF (IOS .NE. 0)  INFO = 1
   10         CONTINUE
              DO 20  I = 1, IPAR(1)
                READ (1, FMT = *, IOSTAT = IOS) (B(I,J), J = 1, IPAR(2))
                IF (IOS .NE. 0)  INFO = 1
   20         CONTINUE
            END IF
            IF (NR(2) .EQ. 5) THEN
              Q(1)  =  .187D1
              Q(4)  = -.244D0
              Q(5)  =  .744D0
              Q(6)  =  .205D0
              Q(8)  =  .589D0
              Q(10) =  .1048D1
            ELSE IF (NR(2) .EQ. 6) THEN
              Q(1)  = .1D-1
              Q(5)  = .1D-1
              Q(8)  = .1D-1
              Q(10) = .1D-1
            ELSE IF (NR(2) .EQ. 7) THEN
              CALL DLASET('U', IPAR(3), IPAR(1), ONE, ONE, C, LDC)
              C(1,3) =  TWO
              C(1,4) =  FOUR
              C(2,4) =  TWO
              Q(1)   =  TWO
              Q(2)   = -ONE
              Q(5)   =  TWO
              Q(6)   = -ONE
              Q(8)   =  TWO
            ELSE IF (NR(2) .EQ. 10) THEN
              C(1,1) = ONE
              C(2,5) = ONE
              Q(1)   = 50.0D0
              Q(3)   = 50.0D0
            ELSE IF (NR(2) .EQ. 11) THEN
              A(10,10) = ONE
              A(11,11) = ONE
              C(1,6) = 15.0D0
              C(2,7) =  7.0D0
              C(2,8) = -.5357D+01
              C(2,9) = -.3943D+01
              C(3,10) = ONE
              C(4,11) = ONE
              Q(1) = 0.5D0
              Q(5) = 5.0D0
              Q(8) = 0.5D0
              Q(10) = 5.0D0
              R(1) = 400.0D0
              R(3) = 700.0D0
              IDENT = '0000'
C
            ELSE IF (NR(2) .EQ. 13) THEN
              DO 24 I = 1, IPAR(1)-6
                READ (1, FMT = *, IOSTAT = IOS)
     1                            (A(I,J), J = 1, IPAR(1)-6)
                IF (IOS .NE. 0)  INFO = 1
   24         CONTINUE
              DO 25 I = 1, IPAR(1)-6
                READ (1, FMT = *, IOSTAT = IOS)
     1                            (B(I,J), J = 1, IPAR(2))
                IF (IOS .NE. 0)  INFO = 1
   25         CONTINUE
              DO 26 I = 1, IPAR(2)
                READ (1, FMT = *, IOSTAT = IOS)
     1                             (C(I,J), J = 1, IPAR(1)-6)
                IF (IOS .NE. 0)  INFO = 1
   26         CONTINUE
              DO 27 I = 1, 6
                A(20+I,20+I) = ONE
                C(6+I,20+I)  = ONE
   27         CONTINUE
              J = 58
              DO 28 I = 7, 12
                READ (1, FMT = *, IOSTAT = IOS) Q(J)
                IF (IOS .NE. 0)  INFO = 1
                J = J + (13 - I)
   28         CONTINUE
              J = 1
              DO 29 I = 1, 6
                READ (1, FMT = *, IOSTAT = IOS) R(J)
                IF (IOS .NE. 0)  INFO = 1
                J = J + (7 - I)
   29         CONTINUE
              DO 31 I = 1, 6
                DO 30 J = 1, 20
                  A(I+20,J) = -C(I,J)
   30           CONTINUE
   31         CONTINUE
              IDENT = '0000'
            END IF
          END IF
          CLOSE(1)
          IF ((NR(2) .EQ. 5) .OR. (NR(2) .EQ. 6)) THEN
            IDENT = '0101'
          ELSE IF ((NR(2) .EQ. 7) .OR. (NR(2) .EQ. 10)) THEN
            IDENT = '0001'
          ELSE IF (NR(2) .EQ. 8) THEN
            IDENT = '0111'
          END IF
C
        ELSE IF (NR(2). EQ. 9) THEN
          A(1,2) = ONE
          A(2,3) = ONE
          A(4,5) = ONE
          A(5,6) = ONE
          B(3,1) = ONE
          B(6,2) = ONE
          C(1,1) = ONE
          C(1,2) = ONE
          C(2,4) = ONE
          C(2,5) = -ONE
          R(1)   = THREE
          R(3)   = ONE
          IF (BPAR(7)) THEN
            S(1,1) = ONE
            S(2,1) = ONE
            S(4,1) = ONE
            S(5,1) = -ONE
          END IF
          IDENT  = '0010'
        ELSE IF (NR(2) .EQ. 12) THEN
          DO 32 I = 1, 10
            A(I,I+1) = ONE
   32     CONTINUE
          A(6,7) = ZERO
          A(8,9) = ZERO
          A(12,12) = ONE
          A(13,13) = ONE
          A(12,1) = -.3318D+01
          A(13,1) = -.15484D+01
          A(6,6) = .7788D+00
          A(8,7) = -.4724D+00
          A(13,7) = .3981D+00
          A(8,8) = .13746D+01
          A(13,8) = .5113D+00
          A(13,9) = .57865D+01
          A(11,11) = .8071D+00
          B(6,1) = ONE
          B(8,2) = ONE
          C(1,1) = .3318D+01
          C(2,1) = .15484D+01
          C(2,7) = -.3981D+00
          C(2,8) = -.5113D+00
          C(2,9) = -.57865D+01
          C(3,12) = ONE
          C(4,13) = ONE
          Q(1) = 0.5D0
          Q(5) = 5.0D0
          Q(8) = 0.5D0
          Q(10) = 5.0D0
          R(1) = 400.0D0
          R(3) = 700.0D0
          IDENT = '0000'
        END IF
C
      ELSE IF (NR(1) .EQ. 2) THEN
        IF (NR(2) .EQ. 1) THEN
          IF (LSAME(DEF, 'D')) DPAR(1) = .1D+07
          A(1,1) = FOUR
          A(2,1) = -.45D1
          A(1,2) = THREE
          A(2,2) = -.35D1
          CALL DLASET('A', IPAR(1), IPAR(2), -ONE, ONE, B, LDB)
          R(1) = DPAR(1)
          Q(1)  = 9.0D0
          Q(2)  = 6.0D0
          Q(3)  = FOUR
          TEMP  = (ONE + SQRT(ONE+FOUR*DPAR(1))) / TWO
          X(1,1) = TEMP*Q(1)
          X(2,1) = TEMP*Q(2)
          X(1,2) = X(2,1)
          X(2,2) = TEMP*Q(3)
          IDENT = '0100'
C
        ELSE IF (NR(2) .EQ. 2) THEN
          IF (LSAME(DEF, 'D')) DPAR(1) = .1D+07
          IF (DPAR(1) .EQ. ZERO) THEN
            INFO = 2
          ELSE
            A(1,1) = .9512D0
            A(2,2) = .9048D0
            CALL DLASET('A', 1, IPAR(2), .4877D1, .4877D1, B, LDB)
            B(2,1) = -.11895D1
            B(2,2) = .3569D1
            R(1)   = ONE / (THREE*DPAR(1))
            R(3)   = THREE*DPAR(1)
            Q(1)   = .5D-2
            Q(3)   = .2D-1
            IDENT  = '0100'
          END IF
C
        ELSE IF (NR(2) .EQ. 3) THEN
          IF (LSAME(DEF,'D'))  DPAR(1) = .1D7
          A(1,2) = DPAR(1)
          B(2,1) = ONE
          X(1,1) = ONE
          X(2,2) = ONE + DPAR(1)*DPAR(1)
          IDENT  = '0111'
C
        ELSE IF (NR(2) .EQ. 4) THEN
          IF (LSAME(DEF,'D'))  DPAR(1) = .1D7
          A(2,2) = ONE
          A(3,3) = THREE
          R(1)   = DPAR(1)
          R(4)   = DPAR(1)
          R(6)   = DPAR(1)
C     .. set C = V ..
          TEMP   = TWO/THREE
          CALL DLASET('A', IPAR(3), IPAR(1), -TEMP, ONE - TEMP, C, LDC)
C     .. and compute A <- C' A C
          CALL DSYMM('L', 'L', IPAR(1), IPAR(1), ONE, C, LDC, A, LDA,
     1               ZERO, DWORK, IPAR(1))
          CALL DSYMM('R', 'L', IPAR(1), IPAR(1), ONE, C, LDC, DWORK,
     1               IPAR(1), ZERO, A, LDA)
          Q(1)   = DPAR(1)
          Q(4)   = DPAR(1)
          Q(6)   = DPAR(1)
          X(1,1) = DPAR(1)
          X(2,2) = DPAR(1) * (ONE + SQRT(FIVE)) / TWO
          X(3,3) = DPAR(1) * (9.0D0 + SQRT(85.0D0)) / TWO
          CALL DSYMM('L', 'L', IPAR(1), IPAR(1), ONE, C, LDC, X, LDX,
     1               ZERO, DWORK, IPAR(1))
          CALL DSYMM('R', 'L', IPAR(1), IPAR(1), ONE, C, LDC, DWORK,
     1               IPAR(1), ZERO, X, LDX)
          IDENT  = '1000'
C
        ELSE IF (NR(2) .EQ. 5) THEN
          IF (LSAME(DEF, 'D')) THEN
            DPAR(4) = .25D0
            DPAR(3) = ONE
            DPAR(2) = ONE
            DPAR(1) = .1D9
          END IF
          IF (DPAR(1) .EQ. ZERO) THEN
            INFO = 2
          ELSE
            TEMP  = DPAR(2) / DPAR(1)
            BETA  = DPAR(3) * TEMP
            ALPHA = ONE - TEMP
            A(1,1) = ALPHA
            CALL DLASET('A', IPAR(1)-1, IPAR(1)-1, ZERO, ONE, A(2,1),
     1                  LDA)
            B(1,1) = BETA
            C(1,4) = ONE
            R(1)  = DPAR(4)
            IF (BETA .EQ. ZERO) THEN
              INFO = 2
            ELSE
              CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ONE, X, LDX)
              BETA   = BETA * BETA
              TEMP   = DPAR(4) * (ALPHA + ONE) * (ALPHA - ONE) + BETA
              X(1,1) = (TEMP + SQRT(TEMP*TEMP + FOUR*BETA*DPAR(4)))
              X(1,1) = X(1,1) / TWO / BETA
            END IF
            IDENT = '0010'
          END IF
        END IF
C
      ELSE IF (NR(1) .EQ. 4) THEN
        IF (NR(2) .EQ. 1) THEN
          IF (LSAME(DEF,'D'))  DPAR(1) = ONE
          CALL DLASET('A', IPAR(1)-1, IPAR(1)-1, ZERO, ONE, A(1,2), LDA)
          B(IPAR(1),1) = ONE
          R(1) = DPAR(1)
          DO 40  I = 1, IPAR(1)
            X(I,I) = DBLE(I)
   40     CONTINUE
          IDENT  = '0110'
        END IF
      END IF
C
      IF (INFO .NE. 0)  GOTO 2001
C     .. set up data in required format ..
C
      IF (BPAR(4)) THEN
C     .. G is to be returned in product form ..
        RDIMM = IPAR(1)
        IF (IDENT(4:4) .EQ. '0') THEN
C       .. invert R using Cholesky factorization, ..
          CALL DPPTRF('L', IPAR(2), R, INFO)
          IF (INFO .EQ. 0) THEN
            CALL DPPTRI('L', IPAR(2), R, INFO)
            IF (IDENT(1:1) .EQ. '0') THEN
C           .. B is not identity matrix ..
              DO 100  I = 1, IPAR(1)
                CALL DSPMV('L', IPAR(2), ONE, R, B(I,1), LDB, ZERO,
     1                     DWORK((I-1)*IPAR(1)+1), 1)
  100         CONTINUE
              CALL DGEMV('T', IPAR(2), IPAR(1), ONE, DWORK, IPAR(1),
     1                   B(1,1), LDB, ZERO, R, 1)
              ISYMM = IPAR(1) + 1
              DO 110  I = 2, IPAR(1)
                CALL DGEMV('T', IPAR(2), IPAR(1), ONE, DWORK, IPAR(1),
     1                     B(I,1), LDB, ZERO, B(1,1), LDB)
                CALL DCOPY(IPAR(1) - I + 1, B(1,I), LDB, R(ISYMM), 1)
                ISYMM = ISYMM + (IPAR(1) - I + 1)
  110         CONTINUE
            END IF
          ELSE
            IF (INFO .GT. 0)  THEN
              INFO = 3
              GOTO 2001
            END IF
          END IF
        ELSE
C       .. R = identity ..
          IF (IDENT(1:1) .EQ. '0') THEN
C         .. B not identity matrix ..
            IF (IPAR(2) .EQ. 1) THEN
              CALL DLASET('L', NSYMM, 1, ZERO, ZERO, R, 1)
              CALL DSPR('L', IPAR(1), ONE, B, 1, R)
            ELSE
              CALL DSYRK('L', 'N', IPAR(1), IPAR(2), ONE, B, LDB, ZERO,
     1                   DWORK, IPAR(1))
              CALL MA02DD('Pack', 'Lower', IPAR(1), DWORK, IPAR(1), R)
            END IF
          ELSE
C         .. B = R = identity ..
            ISYMM = 1
            DO 120  I = IPAR(1), 1, -1
              R(ISYMM) = ONE
              ISYMM = ISYMM + I
  120       CONTINUE
          END IF
        END IF
      ELSE
        RDIMM = IPAR(2)
        IF (IDENT(1:1) .EQ. '1')
     1    CALL DLASET('A', IPAR(1), IPAR(2), ZERO, ONE, B, LDB)
        IF (IDENT(4:4) .EQ. '1') THEN
          ISYMM = 1
          DO 130  I = IPAR(2), 1, -1
            R(ISYMM) = ONE
            ISYMM = ISYMM + I
  130     CONTINUE
        END IF
      END IF
C
      IF (BPAR(1)) THEN
C     .. Q is to be returned in product form ..
        QDIMM = IPAR(1)
        IF (IDENT(3:3) .EQ. '0') THEN
          IF (IDENT(2:2) .EQ. '0') THEN
C         .. C is not identity matrix ..
            DO 140  I = 1, IPAR(1)
              CALL DSPMV('L', IPAR(3), ONE, Q, C(1,I), 1, ZERO,
     1                   DWORK((I-1)*IPAR(1)+1), 1)
  140       CONTINUE
C         .. use Q(1:IPAR(1)) as workspace and compute the first column
C            of Q at the end ..
            ISYMM = IPAR(1) + 1
            DO 150  I = 2, IPAR(1)
              CALL DGEMV('T', IPAR(3), IPAR(1), ONE, DWORK, IPAR(1),
     1                   C(1,I), 1, ZERO, Q(1), 1)
              CALL DCOPY(IPAR(1) - I + 1, Q(I), 1, Q(ISYMM), 1)
              ISYMM = ISYMM + (IPAR(1) - I + 1)
  150       CONTINUE
            CALL DGEMV('T', IPAR(3), IPAR(1), ONE, DWORK, IPAR(1),
     1                 C(1,1), 1, ZERO, Q, 1)
          END IF
        ELSE
C       .. Q = identity ..
          IF (IDENT(2:2) .EQ. '0') THEN
C         .. C is not identity matrix ..
            IF (IPAR(3) .EQ. 1) THEN
              CALL DLASET('L', NSYMM, 1, ZERO, ZERO, Q, 1)
              CALL DSPR('L', IPAR(1), ONE, C, LDC, Q)
            ELSE
              CALL DSYRK('L', 'T', IPAR(1), IPAR(3), ONE, C, LDC, ZERO,
     1                   DWORK, IPAR(1))
              CALL MA02DD('Pack', 'Lower', IPAR(1), DWORK, IPAR(1), Q)
            END IF
          ELSE
C         .. C = Q = identity ..
            ISYMM = 1
            DO 160  I = IPAR(1), 1, -1
              Q(ISYMM) = ONE
              ISYMM    = ISYMM + I
  160       CONTINUE
          END IF
        END IF
      ELSE
        QDIMM = IPAR(3)
        IF (IDENT(2:2) .EQ. '1')
     1    CALL DLASET('A', IPAR(3), IPAR(1), ZERO, ONE, C, LDC)
        IF (IDENT(3:3) .EQ. '1') THEN
          ISYMM = 1
          DO 170  I = IPAR(3), 1, -1
            Q(ISYMM) = ONE
            ISYMM    = ISYMM + I
  170     CONTINUE
        END IF
      END IF
C
C     .. unpack symmetric matrices if required ..
      IF (BPAR(2)) THEN
        ISYMM = (QDIMM * (QDIMM + 1)) / 2
        CALL DCOPY(ISYMM, Q, 1, DWORK, 1)
        CALL MA02DD('Unpack', 'Lower', QDIMM, Q, LDQ, DWORK)
        CALL MA02ED('Lower', QDIMM, Q, LDQ)
      ELSE IF (BPAR(3)) THEN
        CALL MA02DD('Unpack', 'Lower', QDIMM, DWORK, QDIMM, Q)
        CALL MA02ED('Lower', QDIMM, DWORK, QDIMM)
        CALL MA02DD('Pack', 'Upper', QDIMM, DWORK, QDIMM, Q)
      END IF
      IF (BPAR(5)) THEN
        ISYMM = (RDIMM * (RDIMM + 1)) / 2
        CALL DCOPY(ISYMM, R, 1, DWORK, 1)
        CALL MA02DD('Unpack', 'Lower', RDIMM, R, LDR, DWORK)
        CALL MA02ED('Lower', RDIMM, R, LDR)
      ELSE IF (BPAR(6)) THEN
        CALL MA02DD('Unpack', 'Lower', RDIMM, DWORK, RDIMM, R)
        CALL MA02ED('Lower', RDIMM, DWORK, RDIMM)
        CALL MA02DD('Pack', 'Upper', RDIMM, DWORK, RDIMM, R)
      END IF
C
C     ...set VEC...
      VEC(1) = .TRUE.
      VEC(2) = .TRUE.
      VEC(3) = .TRUE.
      VEC(4) = .TRUE.
      VEC(5) = .NOT. BPAR(4)
      VEC(6) = .NOT. BPAR(1)
      VEC(7) = .TRUE.
      VEC(8) = .TRUE.
      VEC(9) = BPAR(7)
      IF (((NR(1) .EQ. 1) .AND. ((NR(2) .EQ. 1)   .OR.
     1     (NR(2) .EQ. 3) .OR.   (NR(2) .EQ. 4))) .OR.
     2    ((NR(1) .EQ. 2) .AND. ((NR(2). EQ. 1)   .OR.
     3     (NR(2) .GE. 3))) .OR. (NR(1) .EQ. 4)) THEN
        VEC(10) = .TRUE.
      END IF
      CHPAR = NOTES(NR(1),NR(2))
      N = IPAR(1)
      M = IPAR(2)
      P = IPAR(3)
C
 2001 CONTINUE
      RETURN
C *** Last line of BB02AD ***
      END
