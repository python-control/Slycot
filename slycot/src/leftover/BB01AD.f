      SUBROUTINE BB01AD(DEF, NR, DPAR, IPAR, BPAR, CHPAR, VEC, N, M, P,
     1                  A, LDA, B, LDB, C, LDC, G, LDG, Q, LDQ, X, LDX,
     2                  DWORK, LDWORK, INFO)
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
C     continuous-time algebraic Riccati equations (CAREs) of the form
C
C       0 = Q + A'X + XA - XGX
C
C     corresponding to the Hamiltonian matrix
C
C            (  A   G  )
C        H = (       T ).
C            (  Q  -A  )
C
C     A,G,Q,X are real N-by-N matrices, Q and G are symmetric and may
C     be given in factored form
C
C                   -1 T                         T
C      (I)   G = B R  B  ,           (II)   Q = C W C .
C
C     Here, C is P-by-N, W P-by-P, B N-by-M, and R M-by-M, where W
C     and R are symmetric. In linear-quadratic optimal control problems,
C     usually W is positive semidefinite and R positive definite.  The
C     factorized form can be used if the CARE is solved using the
C     deflating subspaces of the extended Hamiltonian pencil
C
C                  (  A   0   B  )       (  I   0   0  )
C                  (       T     )       (             )
C        H - s K = (  Q   A   0  )  -  s (  0  -I   0  ) ,
C                  (       T     )       (             )
C                  (  0   B   R  )       (  0   0   0  )
C
C     where I and 0 denote the identity and zero matrix, respectively,
C     of appropriate dimensions.
C
C     NOTE: the formulation of the CARE and the related matrix (pencils)
C           used here does not include CAREs as they arise in robust
C           control (H_infinity optimization).
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
C             This array determines the example for which CAREX returns
C             data. NR(1) is the group of examples.
C             NR(1) = 1 : parameter-free problems of fixed size.
C             NR(1) = 2 : parameter-dependent problems of fixed size.
C             NR(1) = 3 : parameter-free problems of scalable size.
C             NR(1) = 4 : parameter-dependent problems of scalable size.
C             NR(2) is the number of the example in group NR(1).
C             Let NEXi be the number of examples in group i. Currently,
C             NEX1 = 6, NEX2 = 9, NEX3 = 2, NEX4 = 4.
C             1 <= NR(1) <= 4;
C             1 <= NR(2) <= NEXi , where i = NR(1).
C
C     DPAR    (input/output) DOUBLE PRECISION array, dimension (7)
C             Double precision parameter vector. For explanation of the
C             parameters see [1].
C             DPAR(1)           : defines the parameters
C                                 'delta' for NR(1) = 3,
C                                 'q' for NR(1).NR(2) = 4.1,
C                                 'a' for NR(1).NR(2) = 4.2, and
C                                 'mu' for NR(1).NR(2) = 4.3.
C             DPAR(2)           : defines parameters
C                                 'r' for NR(1).NR(2) = 4.1,
C                                 'b' for NR(1).NR(2) = 4.2, and
C                                 'delta' for NR(1).NR(2) = 4.3.
C             DPAR(3)           : defines parameters
C                                 'c' for NR(1).NR(2) = 4.2 and
C                                 'kappa' for NR(1).NR(2) = 4.3.
C             DPAR(j), j=4,5,6,7: These arguments are only used to
C                                 generate Example 4.2 and define in
C                                 consecutive order the intervals
C                                 ['beta_1', 'beta_2'],
C                                 ['gamma_1', 'gamma_2'].
C             NOTE that if DEF = 'D' or 'd', the values of DPAR entries
C             on input are ignored and, on output, they are overwritten
C             with the default parameters.
C
C     IPAR    (input/output) INTEGER array, dimension (3)
C             On input, IPAR(1) determines the actual state dimension,
C             i.e., the order of the matrix A as follows, where
C             NO = NR(1).NR(2).
C             NR(1) = 1 or 2.1-2.8: IPAR(1) is ignored.
C             NO = 2.9            : IPAR(1) = 1 generates the CARE for
C                                   optimal state feedback (default);
C                                   IPAR(1) = 2 generates the Kalman
C                                   filter CARE.
C             NO = 3.1            : IPAR(1) is the number of vehicles
C                                   (parameter 'l' in the description
C                                    in [1]).
C             NO = 3.2, 4.1 or 4.2: IPAR(1) is the order of the matrix
C                                   A.
C             NO = 4.3 or 4.4     : IPAR(1) determines the dimension of
C                                   the second-order system, i.e., the
C                                   order of the stiffness matrix for
C                                   Examples 4.3 and 4.4 (parameter 'l'
C                                   in the description in [1]).
C
C             The order of the output matrix A is N = 2*IPAR(1) for
C             Example 4.3 and N = 2*IPAR(1)-1 for Examples 3.1 and 4.4.
C             NOTE that IPAR(1) is overwritten for Examples 1.1-2.8. For
C             the other examples, IPAR(1) is overwritten if the default
C             parameters are to be used.
C             On output, IPAR(1) contains the order of the matrix A.
C
C             On input, IPAR(2) is the number of colums in the matrix B
C             in (I) (in control problems, the number of inputs of the
C             system). Currently, IPAR(2) is fixed or determined by
C             IPAR(1) for all examples and thus is not referenced on
C             input.
C             On output, IPAR(2) is the number of columns of the
C             matrix B from (I).
C             NOTE that currently IPAR(2) is overwritten and that
C             rank(G) <= IPAR(2).
C
C             On input, IPAR(3) is the number of rows in the matrix C
C             in (II) (in control problems, the number of outputs of the
C             system). Currently, IPAR(3) is fixed or determined by
C             IPAR(1) for all examples and thus is not referenced on
C             input.
C             On output, IPAR(3) contains the number of rows of the
C             matrix C in (II).
C             NOTE that currently IPAR(3) is overwritten and that
C             rank(Q) <= IPAR(3).
C
C     BPAR    (input) BOOLEAN array, dimension (6)
C             This array defines the form of the output of the examples
C             and the storage mode of the matrices G and Q.
C             BPAR(1) = .TRUE.  : G is returned.
C             BPAR(1) = .FALSE. : G is returned in factored form, i.e.,
C                                 B and R from (I) are returned.
C             BPAR(2) = .TRUE.  : The matrix returned in array G (i.e.,
C                                 G if BPAR(1) = .TRUE. and R if
C                                 BPAR(1) = .FALSE.) is stored as full
C                                 matrix.
C             BPAR(2) = .FALSE. : The matrix returned in array G is
C                                 provided in packed storage mode.
C             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix
C                                 returned in array G is stored in upper
C                                 packed mode, i.e., the upper triangle
C                                 of a symmetric n-by-n matrix is stored
C                                 by columns, e.g., the matrix entry
C                                 G(i,j) is stored in the array entry
C                                 G(i+j*(j-1)/2) for i <= j.
C                                 Otherwise, this entry is ignored.
C             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix
C                                 returned in array G is stored in lower
C                                 packed mode, i.e., the lower triangle
C                                 of a symmetric n-by-n matrix is stored
C                                 by columns, e.g., the matrix entry
C                                 G(i,j) is stored in the array entry
C                                 G(i+(2*n-j)*(j-1)/2) for j <= i.
C                                 Otherwise, this entry is ignored.
C             BPAR(4) = .TRUE.  : Q is returned.
C             BPAR(4) = .FALSE. : Q is returned in factored form, i.e.,
C                                 C and W from (II) are returned.
C             BPAR(5) = .TRUE.  : The matrix returned in array Q (i.e.,
C                                 Q if BPAR(4) = .TRUE. and W if
C                                 BPAR(4) = .FALSE.) is stored as full
C                                 matrix.
C             BPAR(5) = .FALSE. : The matrix returned in array Q is
C                                 provided in packed storage mode.
C             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix
C                                 returned in array Q is stored in upper
C                                 packed mode (see above).
C                                 Otherwise, this entry is ignored.
C             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix
C                                 returned in array Q is stored in lower
C                                 packed mode (see above).
C                                 Otherwise, this entry is ignored.
C             NOTE that there are no default values for BPAR.  If all
C             entries are declared to be .TRUE., then matrices G and Q
C             are returned in conventional storage mode, i.e., as
C             N-by-N arrays where the array element Z(I,J) contains the
C             matrix entry Z_{i,j}.
C
C     CHPAR   (input/output) CHARACTER*255
C             On input, this is the name of a data file supplied by the
C             user.
C             In the current version, only Example 4.4 allows a
C             user-defined data file. This file must contain
C             consecutively DOUBLE PRECISION vectors mu, delta, gamma,
C             and kappa. The length of these vectors is determined by
C             the input value for IPAR(1).
C             If on entry, IPAR(1) = L, then mu and delta must each
C             contain L DOUBLE PRECISION values, and gamma and kappa
C             must each contain L-1 DOUBLE PRECISION values.
C             On output, this string contains short information about
C             the chosen example.
C
C     VEC     (output) LOGICAL array, dimension (9)
C             Flag vector which displays the availability of the output
C             data:
C             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and
C             are always .TRUE.
C             VEC(4) refers to A and is always .TRUE.
C             VEC(5) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors B
C             and R from (I) are returned.
C             VEC(6) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors C
C             and W from (II) are returned.
C             VEC(7) refers to G and is always .TRUE.
C             VEC(8) refers to Q and is always .TRUE.
C             VEC(9) refers to X and is .TRUE. if the exact solution
C             matrix is available.
C             NOTE that VEC(i) = .FALSE. for i = 1 to 9 if on exit
C             INFO .NE. 0.
C
C     N       (output) INTEGER
C             The order of the matrices A, X, G if BPAR(1) = .TRUE., and
C             Q if BPAR(4) = .TRUE.
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
C             coefficient matrix A of the CARE.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     B       (output) DOUBLE PRECISION array, dimension (LDB,M)
C             If (BPAR(1) = .FALSE.), then the leading N-by-M part of
C             this array contains the matrix B of the factored form (I)
C             of G. Otherwise, B is used as workspace.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= N.
C
C     C       (output) DOUBLE PRECISION array, dimension (LDC,N)
C             If (BPAR(4) = .FALSE.), then the leading P-by-N part of
C             this array contains the matrix C of the factored form (II)
C             of Q. Otherwise, C is used as workspace.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= P, where P is the number of rows of the matrix C,
C             i.e., the output value of IPAR(3). (For all examples,
C             P <= N, where N equals the output value of the argument
C             IPAR(1), i.e., LDC >= LDA is always safe.)
C
C     G       (output) DOUBLE PRECISION array, dimension (NG)
C             If (BPAR(2) = .TRUE.)  then NG = LDG*N.
C             If (BPAR(2) = .FALSE.) then NG = N*(N+1)/2.
C             If (BPAR(1) = .TRUE.), then array G contains the
C             coefficient matrix G of the CARE.
C             If (BPAR(1) = .FALSE.), then array G contains the 'control
C             weighting matrix' R of G's factored form as in (I). (For
C             all examples, M <= N.) The symmetric matrix contained in
C             array G is stored according to BPAR(2) and BPAR(3).
C
C     LDG     INTEGER
C             If conventional storage mode is used for G, i.e.,
C             BPAR(2) = .TRUE., then G is stored like a 2-dimensional
C             array with leading dimension LDG. If packed symmetric
C             storage mode is used, then LDG is not referenced.
C             LDG >= N if BPAR(2) = .TRUE..
C
C     Q       (output) DOUBLE PRECISION array, dimension (NQ)
C             If (BPAR(5) = .TRUE.)  then NQ = LDQ*N.
C             If (BPAR(5) = .FALSE.) then NQ = N*(N+1)/2.
C             If (BPAR(4) = .TRUE.), then array Q contains the
C             coefficient matrix Q of the CARE.
C             If (BPAR(4) = .FALSE.), then array Q contains the 'output
C             weighting matrix' W of Q's factored form as in (II).
C             The symmetric matrix contained in array Q is stored
C             according to BPAR(5) and BPAR(6).
C
C     LDQ     INTEGER
C             If conventional storage mode is used for Q, i.e.,
C             BPAR(5) = .TRUE., then Q is stored like a 2-dimensional
C             array with leading dimension LDQ. If packed symmetric
C             storage mode is used, then LDQ is not referenced.
C             LDQ >= N if BPAR(5) = .TRUE..
C
C     X       (output) DOUBLE PRECISION array, dimension (LDX,IPAR(1))
C             If an exact solution is available (NR = 1.1, 1.2, 2.1,
C             2.3-2.6, 3.2), then the leading N-by-N part of this array
C             contains the solution matrix X in conventional storage
C             mode. Otherwise, X is not referenced.
C
C     LDX     INTEGER
C             The leading dimension of array X.  LDX >= 1, and
C             LDX >= N if NR = 1.1, 1.2, 2.1, 2.3-2.6, 3.2.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= N*MAX(4,N).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0 : successful exit;
C             < 0 : if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1 : data file could not be opened or had wrong format;
C             = 2 : division by zero;
C             = 3 : G can not be computed as in (I) due to a singular R
C                   matrix.
C
C     REFERENCES
C
C     [1] Abels, J. and Benner, P.
C         CAREX - A Collection of Benchmark Examples for Continuous-Time
C         Algebraic Riccati Equations (Version 2.0).
C         SLICOT Working Note 1999-14, November 1999. Available from
C         http://www.win.tue.nl/niconet/NIC2/reports.html.
C
C     This is an updated and extended version of
C
C     [2] Benner, P., Laub, A.J., and Mehrmann, V.
C         A Collection of Benchmark Examples for the Numerical Solution
C         of Algebraic Riccati Equations I: Continuous-Time Case.
C         Technical Report SPC 95_22, Fak. f. Mathematik,
C         TU Chemnitz-Zwickau (Germany), October 1995.
C
C     NUMERICAL ASPECTS
C
C     If the original data as taken from the literature is given via
C     matrices G and Q, but factored forms are requested as output, then
C     these factors are obtained from Cholesky or LDL' decompositions of
C     G and Q, i.e., the output data will be corrupted by roundoff
C     errors.
C
C     FURTHER COMMENTS
C
C     Some benchmark examples read data from the data files provided
C     with the collection.
C
C     CONTRIBUTOR
C
C     Peter Benner (Universitaet Bremen), November 15, 1999.
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
C     Algebraic Riccati equation, Hamiltonian matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
C     . # of examples available , # of examples with fixed size. .
      INTEGER          NEX1, NEX2, NEX3, NEX4, NMAX
      PARAMETER        ( NMAX = 9, NEX1 = 6, NEX2 = 9, NEX3 = 2,
     1                   NEX4 = 4 )
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, PI
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     1                   THREE = 3.0D0, FOUR = 4.0D0,
     2                   PI = .3141592653589793D1 )
C
C     .. Scalar Arguments ..
      INTEGER          INFO, LDA, LDB, LDC, LDG, LDQ, LDWORK, LDX, M, N,
     $                 P
      CHARACTER        DEF
C
C     .. Array Arguments ..
      INTEGER          IPAR(3), NR(2)
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), DPAR(*), DWORK(*),
     1                 G(*), Q(*), X(LDX,*)
      CHARACTER        CHPAR*255
      LOGICAL          BPAR(6), VEC(9)
C
C     .. Local Scalars ..
      INTEGER          GDIMM, I, IOS, ISYMM, J, K, L, MSYMM, NSYMM, POS,
     1                 PSYMM, QDIMM
      DOUBLE PRECISION APPIND, B1, B2, C1, C2, SUM, TEMP, TTEMP
C
C     ..Local Arrays ..
      INTEGER          MDEF(2,NMAX), NDEF(4,NMAX), NEX(4), PDEF(2,NMAX)
      DOUBLE PRECISION PARDEF(4,NMAX)
      CHARACTER        IDENT*4
      CHARACTER*255    NOTES(4,NMAX)
C
C     .. External Functions ..
C     . BLAS .
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C     . LAPACK .
      LOGICAL          LSAME
      DOUBLE PRECISION DLAPY2
      EXTERNAL         LSAME, DLAPY2
C
C     .. External Subroutines ..
C     . BLAS .
      EXTERNAL         DCOPY, DGEMV, DSCAL, DSPMV, DSPR, DSYMM, DSYRK
C     . LAPACK .
      EXTERNAL         DLASET, DPPTRF, DPPTRI, DPTTRF, DPTTRS, XERBLA
C     . SLICOT .
      EXTERNAL         MA02DD, MA02ED
C
C     .. Intrinsic Functions ..
      INTRINSIC        COS, MAX, MIN, MOD, SQRT
C
C     .. Data Statements ..
C     . default values for dimensions .
      DATA (NEX(I), I = 1, 4) /NEX1, NEX2, NEX3, NEX4/
      DATA (NDEF(1,I), I = 1, NEX1) /2, 2, 4, 8, 9, 30/
      DATA (NDEF(2,I), I = 1, NEX2) /2, 2, 2, 2, 2, 3, 4, 4, 55/
      DATA (NDEF(3,I), I = 1, NEX3) /20, 64/
      DATA (NDEF(4,I), I = 1, NEX4) /21, 100, 30, 211/
      DATA (MDEF(1,I), I = 1, NEX1) /1, 1, 2, 2, 3, 3/
      DATA (MDEF(2,I), I = 1, NEX2) /1, 2, 1, 2, 1, 3, 1, 1, 2/
      DATA (PDEF(1,I), I = 1, NEX1) /2, 2, 4, 8, 9, 5/
      DATA (PDEF(2,I), I = 1, NEX2) /1, 1, 2, 2, 2, 3, 2, 1, 10/
C     . default values for parameters .
      DATA (PARDEF(1,I), I = 1, NEX1) /ZERO, ZERO, ZERO, ZERO, ZERO,
     1                                 ZERO/
      DATA (PARDEF(2,I), I = 1, NEX2) /.1D-5, .1D-7, .1D7, .1D-6, ZERO,
     1                                 .1D7, .1D-5, .1D-5, .1D1/
      DATA (PARDEF(3,I), I = 1, NEX3) /ZERO, ZERO/
      DATA (PARDEF(4,I), I = 1, NEX4) /ONE, .1D-1, FOUR, ZERO/
C     . comments on examples .
      DATA (NOTES(1,I), I = 1, NEX1) /
     1'Laub 1979, Ex.1', 'Laub 1979, Ex.2: uncontrollable-unobservable d
     2ata', 'Beale/Shafai 1989: model of L-1011 aircraft', 'Bhattacharyy
     3a et al. 1983: binary distillation column', 'Patnaik et al. 1980:
     4tubular ammonia reactor', 'Davison/Gesing 1978: J-100 jet engine'/
      DATA (NOTES(2,I), I = 1, NEX2) /
     1'Arnold/Laub 1984, Ex.1: (A,B) unstabilizable as EPS -> 0', 'Arnol
     2d/Laub 1984, Ex.3: control weighting matrix singular as EPS -> 0',
     3'Kenney/Laub/Wette 1989, Ex.2: ARE ill conditioned for EPS -> oo',
     4'Bai/Qian 1994: ill-conditioned Hamiltonian for EPS -> 0', 'Laub 1
     5992: H-infinity problem, eigenvalues  +/- EPS +/- i', 'Petkov et a
     6l. 1987: increasingly badly scaled Hamiltonian as EPS -> oo', 'Cho
     7w/Kokotovic 1976: magnetic tape control system', 'Arnold/Laub 1984
     8, Ex.2: poor sep. of closed-loop spectrum as EPS -> 0', 'IFAC Benc
     9hmark Problem #90-06: LQG design for modified Boing B-767 at flutt
     1er condition'/
      DATA (NOTES(3,I), I = 1, NEX3) /
     1'Laub 1979, Ex.4: string of high speed vehicles', 'Laub 1979, Ex.5
     2: circulant matrices'/
      DATA (NOTES(4,I), I = 1, NEX4) /
     1'Laub 1979, Ex.6: ill-conditioned Riccati equation', 'Rosen/Wang 1
     2992: lq control of 1-dimensional heat flow','Hench et al. 1995: co
     3upled springs, dashpots and masses','Lang/Penzl 1994: rotating axl
     4e' /
C
C     .. Executable Statements ..
C
      INFO = 0
      DO 5 I = 1, 9
        VEC(I) = .FALSE.
    5 CONTINUE
C
      IF ((NR(1) .NE. 1) .AND. (.NOT. (LSAME(DEF,'N')
     1    .OR. LSAME(DEF,'D')))) THEN
        INFO = -1
      ELSE IF ((NR(1) .LT. 1) .OR. (NR(2) .LT. 1) .OR.
     1  (NR(1) .GT. 4) .OR. (NR(2) .GT. NEX(NR(1)))) THEN
        INFO = -2
      ELSE IF (NR(1) .GT. 2) THEN
        IF (.NOT. LSAME(DEF,'N')) IPAR(1) = NDEF(NR(1),NR(2))
        IF (NR(1) .EQ. 3) THEN
          IF (NR(2) .EQ. 1) THEN
            IPAR(2) = IPAR(1)
            IPAR(3) = IPAR(1) - 1
            IPAR(1) = 2*IPAR(1) - 1
          ELSE IF (NR(2) .EQ. 2) THEN
            IPAR(2) = IPAR(1)
            IPAR(3) = IPAR(1)
          ELSE
            IPAR(2) = 1
            IPAR(3) = 1
          END IF
        ELSE IF (NR(1) .EQ. 4) THEN
          IF (NR(2) .EQ. 3) THEN
            L = IPAR(1)
            IPAR(2) = 2
            IPAR(3) = 2*L
            IPAR(1) = 2*L
          ELSE IF (NR(2) .EQ. 4) THEN
            L = IPAR(1)
            IPAR(2) = L
            IPAR(3) = L
            IPAR(1) = 2*L-1
          ELSE
            IPAR(2) = 1
            IPAR(3) = 1
          END IF
        END IF
      ELSE IF ((NR(1) .EQ. 2) .AND. (NR(2) .EQ. 9) .AND.
     1         (IPAR(1) . EQ. 2)) THEN
        IPAR(1) = NDEF(NR(1),NR(2))
        IPAR(2) = MDEF(NR(1),NR(2))
        IPAR(3) = 3
      ELSE
        IPAR(1) = NDEF(NR(1),NR(2))
        IPAR(2) = MDEF(NR(1),NR(2))
        IPAR(3) = PDEF(NR(1),NR(2))
      END IF
      IF (INFO .NE. 0)  GOTO 7
C
      IF (IPAR(1) .LT. 1) THEN
        INFO = -4
      ELSE IF (IPAR(1) .GT. LDA) THEN
        INFO = -12
      ELSE IF (IPAR(1) .GT. LDB) THEN
        INFO = -14
      ELSE IF (IPAR(3) .GT. LDC) THEN
        INFO = -16
      ELSE IF (BPAR(2) .AND. (IPAR(1).GT. LDG)) THEN
        INFO = -18
      ELSE IF (BPAR(5) .AND. (IPAR(1).GT. LDQ)) THEN
        INFO = -20
      ELSE IF (LDX.LT.1) THEN
        INFO = -22
      ELSE IF ((NR(1) .EQ. 1) .AND.
     $        ((NR(2) .EQ. 1) .OR. (NR(2) .EQ.2))) THEN
        IF (IPAR(1) .GT. LDX) INFO = -22
      ELSE IF ((NR(1) .EQ. 2) .AND. (NR(2) .EQ. 1)) THEN
        IF (IPAR(1) .GT. LDX) INFO = -22
      ELSE IF ((NR(1) .EQ. 2) .AND. ((NR(2) .GE. 3) .AND.
     1         (NR(2) .LE. 6))) THEN
        IF (IPAR(1) .GT. LDX) INFO = -22
      ELSE IF ((NR(1) .EQ. 3) .AND. (NR(2) .EQ. 2)) THEN
        IF (IPAR(1) .GT. LDX) INFO = -22
      ELSE IF (LDWORK .LT. N*(MAX(4,N))) THEN
        INFO = -24
      END IF
C
    7 CONTINUE
      IF (INFO .NE. 0)  THEN
        CALL XERBLA( 'BB01AD', -INFO )
        RETURN
      END IF
C
      NSYMM = (IPAR(1)*(IPAR(1)+1))/2
      MSYMM = (IPAR(2)*(IPAR(2)+1))/2
      PSYMM = (IPAR(3)*(IPAR(3)+1))/2
      IF (.NOT. LSAME(DEF,'N')) DPAR(1) = PARDEF(NR(1),NR(2))
C
      CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ZERO, A, LDA)
      CALL DLASET('A', IPAR(1), IPAR(2), ZERO, ZERO, B, LDB)
      CALL DLASET('A', IPAR(3), IPAR(1), ZERO, ZERO, C, LDC)
      CALL DLASET('L', MSYMM, 1, ZERO, ZERO, G, 1)
      CALL DLASET('L', PSYMM, 1, ZERO, ZERO, Q, 1)
C
      IF (NR(1) .EQ. 1) THEN
        IF (NR(2) .EQ. 1) THEN
          A(1,2) = ONE
          B(2,1) = ONE
          Q(1)   = ONE
          Q(3)   = TWO
          IDENT  = '0101'
          CALL DLASET('A', IPAR(1), IPAR(1), ONE, TWO, X, LDX)
C
        ELSE IF (NR(2) .EQ. 2) THEN
          A(1,1) = FOUR
          A(2,1) = -.45D1
          A(1,2) = THREE
          A(2,2) = -.35D1
          CALL DLASET('A', IPAR(1), IPAR(2), -ONE, ONE, B, LDB)
          Q(1)  = 9.0D0
          Q(2)  = 6.0D0
          Q(3)  = FOUR
          IDENT = '0101'
          TEMP  = ONE + SQRT(TWO)
          CALL DLASET('A', IPAR(1), IPAR(1), 6.0D0*TEMP, FOUR*TEMP, X,
     1                LDX)
          X(1,1) = 9.0D0*TEMP
C
        ELSE IF ((NR(2) .GE. 3) .AND. (NR(2) .LE. 6)) THEN
          WRITE (CHPAR(1:11), '(A,I1,A,I1,A)') 'BB01', NR(1), '0',
     1                                          NR(2) , '.dat'
          IF ((NR(2) .EQ. 3) .OR. (NR(2) .EQ. 4)) THEN
            IDENT = '0101'
          ELSE IF (NR(2) .EQ. 5) THEN
            IDENT = '0111'
          ELSE IF (NR(2) .EQ. 6) THEN
            IDENT = '0011'
          END IF
          OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:11))
          IF (IOS .NE. 0) THEN
            INFO = 1
          ELSE IF (NR(2) .LE. 6) THEN
            DO 10  I = 1, IPAR(1)
              READ (1, FMT = *, IOSTAT = IOS)
     1                                 (A(I,J), J = 1, IPAR(1))
              IF (IOS .NE. 0) INFO = 1
   10       CONTINUE
            DO 20  I = 1, IPAR(1)
              READ (1, FMT = *, IOSTAT = IOS)
     1                                 (B(I,J), J = 1, IPAR(2))
              IF (IOS .NE. 0) INFO = 1
   20       CONTINUE
            IF (NR(2) .LE. 4) THEN
              DO 30  I = 1, IPAR(1)
                POS = (I-1)*IPAR(1)
                READ (1, FMT = *, IOSTAT = IOS) (DWORK(POS+J),
     1                                          J = 1,IPAR(1))
   30         CONTINUE
              IF (IOS .NE. 0) THEN
                INFO = 1
              ELSE
                CALL MA02DD('Pack', 'Lower', IPAR(1), DWORK, IPAR(1), Q)
              END IF
            ELSE IF (NR(2) .EQ. 6) THEN
              DO 35  I = 1, IPAR(3)
                READ (1, FMT = *, IOSTAT = IOS)
     1                                   (C(I,J), J = 1, IPAR(1))
                IF (IOS .NE. 0) INFO = 1
   35         CONTINUE
            END IF
            CLOSE(1)
          END IF
        END IF
C
      ELSE IF (NR(1) .EQ. 2) THEN
        IF (NR(2) .EQ. 1) THEN
          A(1,1) =  ONE
          A(2,2) = -TWO
          B(1,1) = DPAR(1)
          CALL DLASET('U', IPAR(3), IPAR(1), ONE, ONE, C, LDC)
          IDENT  = '0011'
          IF (DPAR(1) .NE. ZERO) THEN
            TEMP   = DLAPY2(ONE, DPAR(1))
            X(1,1) = (ONE + TEMP)/DPAR(1)/DPAR(1)
            X(2,1) = ONE/(TWO + TEMP)
            X(1,2) = X(2,1)
            TTEMP  = DPAR(1)*X(1,2)
            TEMP   = (ONE - TTEMP) * (ONE + TTEMP)
            X(2,2) = TEMP / FOUR
          ELSE
            INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 2) THEN
          A(1,1) = -.1D0
          A(2,2) = -.2D-1
          B(1,1) =  .1D0
          B(2,1) =  .1D-2
          B(2,2) =  .1D-1
          CALL DLASET('L', MSYMM, 1, ONE, ONE, G, MSYMM)
          G(1)   = G(1) + DPAR(1)
          C(1,1) = .1D2
          C(1,2) = .1D3
          IDENT  = '0010'
C
        ELSE IF (NR(2) .EQ. 3) THEN
          A(1,2) = DPAR(1)
          B(2,1) = ONE
          IDENT  = '0111'
          IF (DPAR(1) .NE. ZERO) THEN
            TEMP   = SQRT(ONE + TWO*DPAR(1))
            CALL DLASET('A', IPAR(1), IPAR(1), ONE, TEMP, X, LDX)
            X(1,1) = X(1,1)/DPAR(1)
          ELSE
            INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 4) THEN
          TEMP = DPAR(1) + ONE
          CALL DLASET('A', IPAR(1), IPAR(1), ONE, TEMP, A, LDA)
          Q(1) = DPAR(1)**2
          Q(3) = Q(1)
          IDENT = '1101'
          X(1,1) = TWO*TEMP + SQRT(TWO)*(SQRT(TEMP**2 + ONE) + DPAR(1))
          X(1,1) = X(1,1)/TWO
          X(2,2) = X(1,1)
          TTEMP  = X(1,1) - TEMP
          IF (TTEMP .NE. ZERO) THEN
             X(2,1) = X(1,1) / TTEMP
             X(1,2) = X(2,1)
          ELSE
             INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 5) THEN
          A(1,1) = THREE - DPAR(1)
          A(2,1) = FOUR
          A(1,2) = ONE
          A(2,2) = TWO - DPAR(1)
          CALL DLASET('L', IPAR(1), IPAR(2), ONE, ONE, B, LDB)
          Q(1)   = FOUR*DPAR(1) - 11.0D0
          Q(2)   = TWO*DPAR(1)  - 5.0D0
          Q(3)   = TWO*DPAR(1)  - TWO
          IDENT  = '0101'
          CALL DLASET('A', IPAR(1), IPAR(1), ONE, ONE, X, LDX)
          X(1,1) = TWO
C
        ELSE IF (NR(2) .EQ. 6) THEN
          IF (DPAR(1) .NE. ZERO) THEN
            A(1,1) = DPAR(1)
            A(2,2) = DPAR(1)*TWO
            A(3,3) = DPAR(1)*THREE
C     .. set C = V ..
            TEMP   = TWO/THREE
            CALL DLASET('A', IPAR(3), IPAR(1), -TEMP, ONE - TEMP,
     1                  C, LDC)
            CALL DSYMM('L', 'L', IPAR(1), IPAR(1), ONE, C, LDC, A, LDA,
     1                 ZERO, DWORK, IPAR(1))
            CALL DSYMM('R', 'L', IPAR(1), IPAR(1), ONE, C, LDC, DWORK,
     1                 IPAR(1), ZERO, A, LDA)
C     .. G = R ! ..
            G(1) = DPAR(1)
            G(4) = DPAR(1)
            G(6) = DPAR(1)
            Q(1) = ONE/DPAR(1)
            Q(4) = ONE
            Q(6) = DPAR(1)
            IDENT = '1000'
            CALL DLASET('A', IPAR(1), IPAR(1), ZERO, ZERO, X, LDX)
            TEMP   = DPAR(1)**2
            X(1,1) = TEMP + SQRT(TEMP**2 + ONE)
            X(2,2) = TEMP*TWO + SQRT(FOUR*TEMP**2 + DPAR(1))
            X(3,3) = TEMP*THREE + DPAR(1)*SQRT(9.0D0*TEMP + ONE)
            CALL DSYMM('L', 'L', IPAR(1), IPAR(1), ONE, C, LDC, X, LDX,
     1                  ZERO, DWORK, IPAR(1))
            CALL DSYMM('R', 'L', IPAR(1), IPAR(1), ONE, C, LDC, DWORK,
     1                 IPAR(1), ZERO, X, LDX)
          ELSE
            INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 7) THEN
          IF (DPAR(1) .NE. ZERO) THEN
            A(1,2) =  .400D0
            A(2,3) =  .345D0
            A(3,2) = -.524D0/DPAR(1)
            A(3,3) = -.465D0/DPAR(1)
            A(3,4) =  .262D0/DPAR(1)
            A(4,4) = -ONE/DPAR(1)
            B(4,1) =  ONE/DPAR(1)
            C(1,1) =  ONE
            C(2,3) =  ONE
            IDENT  = '0011'
          ELSE
            INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 8) THEN
          A(1,1) = -DPAR(1)
          A(2,1) = -ONE
          A(1,2) =  ONE
          A(2,2) = -DPAR(1)
          A(3,3) =  DPAR(1)
          A(4,3) = -ONE
          A(3,4) =  ONE
          A(4,4) =  DPAR(1)
          CALL DLASET('L', IPAR(1), IPAR(2), ONE, ONE, B, LDB)
          CALL DLASET('U', IPAR(3), IPAR(1), ONE, ONE, C, LDC)
          IDENT = '0011'
C
        ELSE IF (NR(2) .EQ. 9) THEN
          IF (IPAR(3) .EQ. 10) THEN
C     .. read LQR CARE ...
            WRITE (CHPAR(1:12), '(A,I1,A,I1,A)') 'BB01', NR(1), '0',
     1                                           NR(2), '1.dat'
            OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:12))
            IF (IOS .NE. 0) THEN
              INFO = 1
            ELSE
              DO 36 I = 1, 27, 2
                READ (1, FMT = *, IOSTAT = IOS)
     1               ((A(I+J,I+K), K = 0, 1), J = 0, 1)
                IF (IOS .NE. 0) INFO = 1
   36         CONTINUE
              DO 37 I = 30, 44, 2
                READ (1, FMT = *, IOSTAT = IOS)
     1               ((A(I+J,I+K), K = 0, 1), J = 0, 1)
                IF (IOS .NE. 0) INFO = 1
   37         CONTINUE
              DO 38 I = 1, IPAR(1)
                READ (1, FMT = *, IOSTAT = IOS)
     1               (A(I,J), J = 46, IPAR(1))
                IF (IOS .NE. 0) INFO = 1
   38         CONTINUE
              A(29,29) = -.5301D1
              B(48,1) = .8D06
              B(51,2) = .8D06
              G(1) = .3647D03
              G(3) = .1459D02
              DO 39 I = 1,6
                READ (1, FMT = *, IOSTAT = IOS)
     1               (C(I,J), J = 1,45)
                IF (IOS .NE. 0) INFO = 1
   39         CONTINUE
              C(7,47)  = ONE
              C(8,46)  = ONE
              C(9,50)  = ONE
              C(10,49) = ONE
              Q(11) = .376D-13
              Q(20) = .120D-12
              Q(41) = .245D-11
            END IF
          ELSE
C     .. read Kalman filter CARE ..
            WRITE (CHPAR(1:12), '(A,I1,A,I1,A)') 'BB01', NR(1), '0',
     1                                           NR(2), '2.dat'
            OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:12))
            IF (IOS .NE. 0) THEN
              INFO = 1
            ELSE
              DO 40 I = 1, 27, 2
                READ (1, FMT = *, IOSTAT = IOS)
     1               ((A(I+K,I+J), K = 0, 1), J = 0, 1)
                IF (IOS .NE. 0) INFO = 1
   40         CONTINUE
              DO 41 I = 30, 44, 2
                READ (1, FMT = *, IOSTAT = IOS)
     1               ((A(I+K,I+J), K = 0, 1), J = 0, 1)
                IF (IOS .NE. 0) INFO = 1
   41         CONTINUE
              DO 42 I = 1, IPAR(1)
                READ (1, FMT = *, IOSTAT = IOS)
     1               (A(J,I), J = 46, IPAR(1))
                IF (IOS .NE. 0) INFO = 1
   42         CONTINUE
              A(29,29) = -.5301D1
              DO 43 J = 1, IPAR(2)
                READ (1, FMT = *, IOSTAT = IOS)
     1                                 (B(I,J), I = 1, IPAR(1))
                IF (IOS .NE. 0) INFO = 1
   43         CONTINUE
              G(1) = .685D-5
              G(3) = .373D3
              C(1,52) = .3713
              C(1,53) = .1245D1
              C(2,48) = .8D6
              C(2,54) = ONE
              C(3,51) = .8D6
              C(3,55) = ONE
              Q(1) = .28224D5
              Q(4) = .2742D-4
              Q(6) = .6854D-3
            END IF
          END IF
          CLOSE(1)
          IDENT = '0000'
        END IF
C
      ELSE IF (NR(1) .EQ. 3) THEN
        IF (NR(2) .EQ. 1) THEN
          DO 45  I = 1, IPAR(1)
            IF (MOD(I,2) .EQ. 1) THEN
              A(I,I)       = -ONE
              B(I,(I+1)/2) =  ONE
            ELSE
              A(I,I-1) =  ONE
              A(I,I+1) = -ONE
              C(I/2,I) =  ONE
            END IF
   45     CONTINUE
          ISYMM = 1
          DO 50  I = IPAR(3), 1, -1
            Q(ISYMM) = 10.0D0
            ISYMM    = ISYMM + I
   50     CONTINUE
          IDENT = '0001'
C
        ELSE IF (NR(2) .EQ. 2) THEN
          DO 60  I = 1, IPAR(1)
            A(I,I) = -TWO
            IF (I .LT. IPAR(1)) THEN
              A(I,I+1) = ONE
              A(I+1,I) = ONE
            END IF
   60	    CONTINUE
          A(1,IPAR(1)) = ONE
          A(IPAR(1),1) = ONE
          IDENT = '1111'
          TEMP = TWO * PI / DBLE(IPAR(1))
          DO 70  I = 1, IPAR(1)
            DWORK(I)   = COS(TEMP*DBLE(I-1))
            DWORK(IPAR(1)+I) = -TWO + TWO*DWORK(I) +
     1                 SQRT(5.0D0 + FOUR*DWORK(I)*(DWORK(I) - TWO))
   70     CONTINUE
          DO 90  J = 1, IPAR(1)
            DO 80  I = 1, IPAR(1)
               DWORK(2*IPAR(1)+I) = COS(TEMP*DBLE(I-1)*DBLE(J-1))
   80       CONTINUE
            X(J,1) = DDOT(IPAR(1), DWORK(IPAR(1)+1), 1,
     1                    DWORK(2*IPAR(1)+1), 1)/DBLE(IPAR(1))
   90     CONTINUE
C         .. set up circulant solution matrix ..
          DO 100  I = 2, IPAR(1)
            CALL DCOPY(IPAR(1)-I+1, X(1,1),   1, X(I,I), 1)
            CALL DCOPY(I-1, X(IPAR(1)-I+2,1), 1, X(1,I), 1)
  100     CONTINUE
        END IF
C
      ELSE IF (NR(1) .EQ. 4) THEN
        IF (NR(2) .EQ. 1) THEN
C       .. set up remaining parameter ..
          IF (.NOT. LSAME(DEF,'N')) THEN
            DPAR(1) = ONE
            DPAR(2) = ONE
          END IF
          CALL DLASET('A', IPAR(1)-1, IPAR(1)-1, ZERO, ONE, A(1,2), LDA)
          B(IPAR(1),1) = ONE
          C(1,1) = ONE
          Q(1)   = DPAR(1)
          G(1)   = DPAR(2)
          IDENT  = '0000'
C
        ELSE IF (NR(2) .EQ. 2) THEN
C         .. set up remaining parameters ..
          APPIND = DBLE(IPAR(1) + 1)
          IF (.NOT. LSAME(DEF,'N')) THEN
            DPAR(1) = PARDEF(NR(1), NR(2))
            DPAR(2) = ONE
            DPAR(3) = ONE
            DPAR(4) = .2D0
            DPAR(5) = .3D0
            DPAR(6) = .2D0
            DPAR(7) = .3D0
          END IF
C         .. set up stiffness matrix ..
          TEMP = -DPAR(1)*APPIND
          CALL DLASET('A', IPAR(1), IPAR(1), ZERO, TWO*TEMP, A, LDA)
          DO 110  I = 1, IPAR(1) - 1
            A(I+1,I) = -TEMP
            A(I,I+1) = -TEMP
  110     CONTINUE
C         .. set up Gramian, stored by diagonals ..
          TEMP = ONE/(6.0D0*APPIND)
          CALL DLASET('L', IPAR(1), 1, FOUR*TEMP, FOUR*TEMP, DWORK,
     1                IPAR(1))
          CALL DLASET('L', IPAR(1)-1, 1, TEMP, TEMP, DWORK(IPAR(1)+1),
     1                IPAR(1))
          CALL DPTTRF(IPAR(1), DWORK(1), DWORK(IPAR(1)+1), INFO)
C         .. A = (inverse of Gramian) * (stiffness matrix) ..
          CALL DPTTRS(IPAR(1), IPAR(1), DWORK(1), DWORK(IPAR(1)+1),
     1                A, LDA, INFO)
C         .. compute B, C ..
          DO 120  I = 1, IPAR(1)
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
  120     CONTINUE
          CALL DSCAL(IPAR(1), DPAR(2), B(1,1), 1)
          CALL DSCAL(IPAR(1), DPAR(3), C(1,1), LDC)
          CALL DPTTRS(IPAR(1), 1, DWORK(1), DWORK(IPAR(1)+1), B, LDB,
     1                INFO)
          IDENT = '0011'
C
        ELSE IF (NR(2) .EQ. 3) THEN
C         .. set up remaining parameters ..
          IF (.NOT. LSAME(DEF,'N')) THEN
            DPAR(1) = PARDEF(NR(1),NR(2))
            DPAR(2) = FOUR
            DPAR(3) = ONE
          END IF
          IF (DPAR(1) . NE. 0) THEN
            CALL DLASET('A', L, L, ZERO, ONE, A(1,L+1), LDA)
            TEMP  = DPAR(3) / DPAR(1)
            A(L+1,1) = -TEMP
            A(L+1,2) =  TEMP
            A(IPAR(1),L-1) =  TEMP
            A(IPAR(1),L)   = -TEMP
            TTEMP = TWO*TEMP
            DO 130  I = 2, L-1
              A(L+I,I)   = -TTEMP
              A(L+I,I+1) =  TEMP
              A(L+I,I-1) =  TEMP
  130       CONTINUE
            CALL DLASET('A', L, L, ZERO, -DPAR(2)/DPAR(1), A(L+1,L+1),
     1                  LDA)
            B(L+1,1) =  ONE / DPAR(1)
            B(IPAR(1),IPAR(2)) = -ONE / DPAR(1)
            IDENT = '0111'
          ELSE
            INFO = 2
          END IF
C
        ELSE IF (NR(2) .EQ. 4) THEN
          IF (.NOT. LSAME(DEF,'N')) WRITE (CHPAR(1:11), '(A,I1,A,I1,A)')
     1                               'BB01', NR(1), '0', NR(2), '.dat'
          OPEN(1, IOSTAT = IOS, STATUS = 'OLD', FILE = CHPAR(1:11))
          IF (IOS .NE. 0) THEN
            INFO = 1
          ELSE
            READ (1, FMT = *, IOSTAT = IOS) (DWORK(I), I = 1, 4*L-2)
            IF (IOS .NE. 0)  INFO = 1
          END IF
          CLOSE(1)
          IF (INFO .EQ. 0) THEN
            CALL DLASET('A', L-1, L-1, ZERO, ONE, A(L+1,2), LDA)
            POS    = 2*L + 1
            A(1,2) = - DWORK(POS) / DWORK(1)
            DO 140  I = 2, L
              TEMP  = DWORK(POS) / DWORK(I-1)
              TTEMP = DWORK(POS) / DWORK(I)
              IF (I .GT. 2)  A(I-1,I) = TEMP
              A(I,I) = -(TEMP + TTEMP)
              IF (I .LT. L)  A(I+1,I) = TTEMP
              POS = POS + 1
  140       CONTINUE
            POS    = L
            TEMP   = DWORK(POS+1) / DWORK(1)
            A(1,1) = -TEMP
            DO 160  I = 2, L
              TTEMP  = TEMP
              TEMP   = DWORK(POS+I) / DWORK(I)
              SUM    = TTEMP - TEMP
              A(I,1) = -SUM
              A(I,I) = A(I,I) - TEMP
              DO 150  J = 2, I-2
                A(I,J) = SUM
  150         CONTINUE
              IF (I .GT. 2)  A(I,I-1) = A(I,I-1) + SUM
  160       CONTINUE
            POS      = 3*L
            A(1,L+1) = -DWORK(3*L)/DWORK(1)
            DO 170  I = 2, L
              TEMP  = DWORK(POS) / DWORK(I-1)
              TTEMP = DWORK(POS) / DWORK(I)
              IF (I .GT. 2)  A(I-1,L+I-1) = TEMP
              A(I,L+I-1) = -(TEMP + TTEMP)
              IF (I .LT. L)  A(I+1,L+I-1) = TTEMP
              POS = POS + 1
  170       CONTINUE
            B(1,1) = ONE/DWORK(1)
            DO 180  I = 1, L
              TEMP = ONE/DWORK(I)
              IF (I .GT. 1)  B(I,I)   = -TEMP
              IF (I .LT. L)  B(I+1,I) =  TEMP
  180       CONTINUE
            C(1,1) = ONE
            Q(1)   = ONE
            POS    = 2*L - 1
            ISYMM  = L + 1
            DO 190  I = 2, L
              TEMP       = DWORK(POS+I)
              TTEMP      = DWORK(POS+L+I-1)
              C(I,I)     = TEMP
              C(I,L+I-1) = TTEMP
              Q(ISYMM)   = ONE / (TEMP*TEMP + TTEMP*TTEMP)
              ISYMM      = ISYMM + L - I + 1
  190       CONTINUE
            IDENT = '0001'
          END IF
        END IF
      END IF
C
      IF (INFO .NE. 0)  GOTO 2001
C     .. set up data in required format ..
C
      IF (BPAR(1)) THEN
C     .. G is to be returned in product form ..
        GDIMM = IPAR(1)
        IF (IDENT(4:4) .EQ. '0') THEN
C       .. invert R using Cholesky factorization, store in G ..
          CALL DPPTRF('L', IPAR(2), G, INFO)
          IF (INFO .EQ. 0) THEN
            CALL DPPTRI('L', IPAR(2), G, INFO)
            IF (IDENT(1:1) .EQ. '0') THEN
C         .. B is not identity matrix ..
              DO 200  I = 1, IPAR(1)
                CALL DSPMV('L', IPAR(2), ONE, G, B(I,1), LDB, ZERO,
     1                     DWORK((I-1)*IPAR(1)+1), 1)
  200         CONTINUE
              CALL DGEMV('T', IPAR(2), IPAR(1), ONE, DWORK, IPAR(1),
     1                   B(1,1), LDB, ZERO, G, 1)
              ISYMM = IPAR(1) + 1
              DO 210  I = 2, IPAR(1)
                CALL DGEMV('T', IPAR(2), IPAR(1), ONE, DWORK, IPAR(1),
     1                     B(I,1), LDB, ZERO, B(1,1), LDB)
                CALL DCOPY(IPAR(1) - I + 1, B(1,I), LDB, G(ISYMM), 1)
                ISYMM = ISYMM + (IPAR(1) - I + 1)
  210         CONTINUE
            END IF
          ELSE
            IF (INFO .GT. 0) THEN
              INFO = 3
              GOTO 2001
            END IF
          END IF
        ELSE
C       .. R = identity ..
          IF (IDENT(1:1) .EQ. '0') THEN
C         .. B is not identity matrix ..
            IF (IPAR(2) .EQ. 1) THEN
              CALL DLASET('L', NSYMM, 1, ZERO, ZERO, G, 1)
              CALL DSPR('L', IPAR(1), ONE, B, 1, G)
            ELSE
              CALL DSYRK('L', 'N', IPAR(1), IPAR(2), ONE,
     1                    B, LDB, ZERO, DWORK, IPAR(1))
              CALL MA02DD('Pack', 'Lower', IPAR(1), DWORK, IPAR(1), G)
            END IF
          ELSE
C         .. B = R = identity ..
            ISYMM = 1
            DO 220  I = IPAR(1), 1, -1
              G(ISYMM) = ONE
              ISYMM = ISYMM + I
  220       CONTINUE
          END IF
        END IF
      ELSE
        GDIMM = IPAR(2)
        IF (IDENT(1:1) .EQ. '1')
     1    CALL DLASET('A', IPAR(1), IPAR(2), ZERO, ONE, B, LDB)
        IF (IDENT(4:4) .EQ. '1') THEN
          ISYMM = 1
          DO 230  I = IPAR(2), 1, -1
            G(ISYMM) = ONE
            ISYMM = ISYMM + I
  230     CONTINUE
        END IF
      END IF
C
      IF (BPAR(4)) THEN
C     .. Q is to be returned in product form ..
        QDIMM = IPAR(1)
        IF (IDENT(3:3) .EQ. '0') THEN
          IF (IDENT(2:2) .EQ. '0') THEN
C         .. C is not identity matrix ..
            DO 240  I = 1, IPAR(1)
              CALL DSPMV('L', IPAR(3), ONE, Q, C(1,I), 1, ZERO,
     1                   DWORK((I-1)*IPAR(1)+1), 1)
  240       CONTINUE
C         .. use Q(1:IPAR(1)) as workspace and compute the first column
C            of Q in the end ..
            ISYMM = IPAR(1) + 1
            DO 250  I = 2, IPAR(1)
              CALL DGEMV('T', IPAR(3), IPAR(1), ONE, DWORK, IPAR(1),
     1                   C(1,I), 1, ZERO, Q(1), 1)
              CALL DCOPY(IPAR(1) - I + 1, Q(I), 1, Q(ISYMM), 1)
              ISYMM = ISYMM + (IPAR(1) - I + 1)
  250       CONTINUE
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
              CALL DSYRK('L', 'T', IPAR(1), IPAR(3), ONE, C, LDC,
     1                   ZERO, DWORK, IPAR(1))
              CALL MA02DD('Pack', 'Lower', IPAR(1), DWORK, IPAR(1), Q)
            END IF
          ELSE
C         .. C = Q = identity ..
            ISYMM = 1
            DO 260  I = IPAR(1), 1, -1
              Q(ISYMM) = ONE
              ISYMM    = ISYMM + I
  260       CONTINUE
          END IF
        END IF
      ELSE
        QDIMM = IPAR(3)
        IF (IDENT(2:2) .EQ. '1')
     1    CALL DLASET('A', IPAR(3), IPAR(1), ZERO, ONE, C, LDC)
        IF (IDENT(3:3) .EQ. '1') THEN
          ISYMM = 1
          DO 270  I = IPAR(3), 1, -1
            Q(ISYMM) = ONE
            ISYMM    = ISYMM + I
  270     CONTINUE
        END IF
      END IF
C
C     .. unpack symmetric matrices if desired ..
      IF (BPAR(2)) THEN
        ISYMM = (GDIMM * (GDIMM + 1)) / 2
        CALL DCOPY(ISYMM, G, 1, DWORK, 1)
        CALL MA02DD('Unpack', 'Lower', GDIMM, G, LDG, DWORK)
        CALL MA02ED('Lower', GDIMM, G, LDG)
      ELSE IF (BPAR(3)) THEN
        CALL MA02DD('Unpack', 'Lower', GDIMM, DWORK, GDIMM, G)
        CALL MA02ED('Lower', GDIMM, DWORK, GDIMM)
        CALL MA02DD('Pack', 'Upper', GDIMM, DWORK, GDIMM, G)
      END IF
      IF (BPAR(5)) THEN
        ISYMM = (QDIMM * (QDIMM + 1)) / 2
        CALL DCOPY(ISYMM, Q, 1, DWORK, 1)
        CALL MA02DD('Unpack', 'Lower', QDIMM, Q, LDQ, DWORK)
        CALL MA02ED('Lower', QDIMM, Q, LDQ)
      ELSE IF (BPAR(6)) THEN
        CALL MA02DD('Unpack', 'Lower', QDIMM, DWORK, QDIMM, Q)
        CALL MA02ED('Lower', QDIMM, DWORK, QDIMM)
        CALL MA02DD('Pack', 'Upper', QDIMM, DWORK, QDIMM, Q)
      END IF
C
C     ...set VEC...
      VEC(1) = .TRUE.
      VEC(2) = .TRUE.
      VEC(3) = .TRUE.
      VEC(4) = .TRUE.
      VEC(5) = .NOT. BPAR(1)
      VEC(6) = .NOT. BPAR(4)
      VEC(7) = .TRUE.
      VEC(8) = .TRUE.
      IF (NR(1) .EQ. 1) THEN
        IF ((NR(2) .EQ. 1) .OR. (NR(2) .EQ. 2)) VEC(9) = .TRUE.
      ELSE IF (NR(1) .EQ. 2) THEN
        IF ((NR(2) .EQ. 1) .OR. ((NR(2) .GE. 3) .AND. (NR(2) .LE. 6)))
     1     VEC(9) = .TRUE.
      ELSE IF (NR(1) .EQ. 3) THEN
        IF (NR(2) .EQ. 2) VEC(9) = .TRUE.
      END IF
      CHPAR  = NOTES(NR(1),NR(2))
      N = IPAR(1)
      M = IPAR(2)
      P = IPAR(3)
 2001 CONTINUE
      RETURN
C *** Last line of BB01AD ***
      END
