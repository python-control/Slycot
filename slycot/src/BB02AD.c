/* BB02AD.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b7 = 0.;
static integer c__1 = 1;
static doublereal c_b34 = 1.;
static doublereal c_b38 = -4.;
static doublereal c_b40 = 11.;
static integer c__5 = 5;
static doublereal c_b112 = -1.;
static doublereal c_b118 = 4.877;

/* Subroutine */ int bb02ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *bpar, char *chpar, logical *vec, integer *n, 
	integer *m, integer *p, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *q, integer *
	ldq, doublereal *r__, integer *ldr, doublereal *s, integer *lds, 
	doublereal *x, integer *ldx, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen def_len, ftnlen chpar_len)
{
    /* Initialized data */

    static integer nex[4] = { 13,5,0,1 };
    static integer ndef[52]	/* was [4][13] */ = { 2,2,0,100,2,2,0,0,2,2,0,
	    0,3,3,0,0,4,4,0,0,4,0,0,0,4,0,0,0,5,0,0,0,6,0,0,0,9,0,0,0,11,0,0,
	    0,13,0,0,0,26 };
    static integer mdef[26]	/* was [2][13] */ = { 1,1,2,2,1,1,2,3,2,1,2,0,
	    4,0,2,0,2,0,3,0,2,0,2,0,6 };
    static integer pdef[26]	/* was [2][13] */ = { 1,2,2,2,2,2,3,3,4,1,4,0,
	    4,0,5,0,2,0,2,0,4,0,4,0,12 };
    static struct {
	char e_1[510];
	char fill_2[255];
	char e_3[765];
	char fill_4[510];
	char e_5[510];
	char fill_6[510];
	char e_7[510];
	char fill_8[510];
	char e_9[510];
	char fill_10[510];
	char e_11[255];
	char fill_12[765];
	char e_13[255];
	char fill_14[765];
	char e_15[255];
	char fill_16[765];
	char e_17[255];
	char fill_18[765];
	char e_19[255];
	char fill_20[765];
	char e_21[255];
	char fill_22[765];
	char e_23[255];
	char fill_24[765];
	char e_25[255];
	char fill_26[765];
	} equiv_26 = { "Van Dooren 1981, Ex. II: singular R matrix          "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Laub 1979, Ex. 2: uncontrollable-unob"
		"servable data                                               "
		"                                                            "
		"                                                            "
		"                                      ", {0}, "Pappas et al."
		" 1980, Ex. 3                                                "
		"                                                            "
		"                                                            "
		"                                                            "
		"  Ionescu/Weiss 1992 : singular R matrix, nonzero S matrix  "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 Laub 1979, Ex. 3: increasingly ill-conditio"
		"ned R-matrix                                                "
		"                                                            "
		"                                                            "
		"                                ", {0}, "Jonckheere 1981: (A"
		",B) controllable, no solution X <= 0                        "
		"                                                            "
		"                                                            "
		"                                                        incr"
		"easingly bad scaled system as eps -> oo                     "
		"                                                            "
		"                                                            "
		"                                                            "
		"           ", {0}, "Sun 1998: R singular, Q non-definite    "
		"                                                            "
		"                                                            "
		"                                                            "
		"                                   Petkov et. al. 1989 : inc"
		"reasingly bad scaling as eps -> oo                          "
		"                                                            "
		"                                                            "
		"                                                  ", {0}, 
		"Ackerson/Fu 1970 : satellite control problem               "
		"                                                            "
		"                                                            "
		"                                                            "
		"                Pappas et al. 1980: process control of paper"
		" machine                                                    "
		"                                                            "
		"                                                            "
		"                               ", {0}, "Litkouhi 1983 : syst"
		"em with slow and fast modes                                 "
		"                                                            "
		"                                                            "
		"                                                       ", {0},
		 "Lu/Lin 1993, Ex. 4.3                                      "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 ", {0}, "Gajic/Shen 1993, Section 2.7.4: ch"
		"emical plant                                                "
		"                                                            "
		"                                                            "
		"                                         ", {0}, "Davison/Wa"
		"ng 1974: nonzero S matrix                                   "
		"                                                            "
		"                                                            "
		"                                                            "
		"     ", {0}, "Patnaik et al. 1980: tubular ammonia reactor  "
		"                                                            "
		"                                                            "
		"                                                            "
		"                             ", {0}, "Sima 1996, Sec. 1.2.2:"
		" paper machine model error integrators                      "
		"                                                            "
		"                                                            "
		"                                                     ", {0}, 
		"Sima 1996, Ex. 2.6: paper machine model with with disturban"
		"ces                                                         "
		"                                                            "
		"                                                            "
		"                ", {0}, "Power plant model, Katayama et al.,"
		" 1985                                                       "
		"                                                            "
		"                                                            "
		"                                        " };

#define notes ((char *)&equiv_26)


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, s_dim1, 
	    s_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, ios;
    static doublereal beta, temp;
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), ma02dd_(char *, 
	    char *, integer *, doublereal *, integer *, doublereal *, ftnlen, 
	    ftnlen), ma02ed_(char *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static char ident[4];
    static integer qdimm, rdimm;
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), dspmv_(char *, integer *, doublereal *, doublereal 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dsymm_(char *, char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen, ftnlen), dsyrk_(char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer isymm, msymm, nsymm, psymm;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpptrf_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpptri_(char *, integer *, 
	    doublereal *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___12 = { 1, 1, 1, 0, 0 };
    static cilist io___14 = { 1, 1, 1, 0, 0 };
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___16 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___18 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };



/*     SLICOT RELEASE 5.0. */

/*     Copyright (c) 2002-2009 NICONET e.V. */

/*     This program is free software: you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation, either version 2 of */
/*     the License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see */
/*     <http://www.gnu.org/licenses/>. */

/*     PURPOSE */

/*     To generate the benchmark examples for the numerical solution of */
/*     discrete-time algebraic Riccati equations (DAREs) of the form */

/*            T                T               T    -1  T       T */
/*     0  =  A X A  -  X  -  (A X B + S) (R + B X B)  (B X A + S )  +  Q */

/*     as presented in [1]. Here, A,Q,X are real N-by-N matrices, B,S are */
/*     N-by-M, and R is M-by-M. The matrices Q and R are symmetric and Q */
/*     may be given in factored form */

/*                   T */
/*     (I)    Q  =  C Q0 C . */

/*     Here, C is P-by-N and Q0 is P-by-P. If R is nonsingular and S = 0, */
/*     the DARE can be rewritten equivalently as */

/*                  T             -1 */
/*     0  =  X  -  A X (I_n + G X)  A  -  Q, */

/*     where I_n is the N-by-N identity matrix and */

/*                   -1  T */
/*     (II)   G = B R   B . */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER */
/*             This parameter specifies if the default parameters are */
/*             to be used or not. */
/*             = 'N' or 'n' : The parameters given in the input vectors */
/*                            xPAR (x = 'D', 'I', 'B', 'CH') are used. */
/*             = 'D' or 'd' : The default parameters for the example */
/*                            are used. */
/*             This parameter is not meaningful if NR(1) = 1. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension (2) */
/*             This array determines the example for which DAREX returns */
/*             data. NR(1) is the group of examples. */
/*             NR(1) = 1 : parameter-free problems of fixed size. */
/*             NR(1) = 2 : parameter-dependent problems of fixed size. */
/*             NR(1) = 3 : parameter-free problems of scalable size. */
/*             NR(1) = 4 : parameter-dependent problems of scalable size. */
/*             NR(2) is the number of the example in group NR(1). */
/*             Let NEXi be the number of examples in group i. Currently, */
/*             NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1. */
/*             1 <= NR(1) <= 4; */
/*             0 <= NR(2) <= NEXi, where i = NR(1). */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (4) */
/*             Double precision parameter vector. For explanation of the */
/*             parameters see [1]. */
/*             DPAR(1) defines the parameter 'epsilon' for */
/*             examples NR = 2.2,2.3,2.4, the parameter 'tau' */
/*             for NR = 2.5, and the 1-by-1 matrix R for NR = 2.1,4.1. */
/*             For Example 2.5, DPAR(2) - DPAR(4) define in */
/*             consecutive order 'D', 'K', and 'r'. */
/*             NOTE that DPAR is overwritten with default values */
/*             if DEF = 'D' or 'd'. */

/*     IPAR    (input/output) INTEGER array, dimension (3) */
/*             On input, IPAR(1) determines the actual state dimension, */
/*             i.e., the order of the matrix A as follows: */
/*             NR(1) = 1, NR(1) = 2   : IPAR(1) is ignored. */
/*             NR = NR(1).NR(2) = 4.1 : IPAR(1) determines the order of */
/*                                      the output matrix A. */
/*             NOTE that IPAR(1) is overwritten for Examples 1.1-2.3. For */
/*             the other examples, IPAR(1) is overwritten if the default */
/*             parameters are to be used. */
/*             On output, IPAR(1) contains the order of the matrix A. */

/*             On input, IPAR(2) is the number of colums in the matrix B */
/*             and the order of the matrix R (in control problems, the */
/*             number of inputs of the system). Currently, IPAR(2) is */
/*             fixed for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(2) is the number of columns of the */
/*             matrix B from (I). */

/*             On input, IPAR(3) is the number of rows in the matrix C */
/*             (in control problems, the number of outputs of the */
/*             system). Currently, IPAR(3) is fixed for all examples */
/*             and thus is not referenced on input. */
/*             On output, IPAR(3) is the number of rows of the matrix C */
/*             from (I). */

/*             NOTE that IPAR(2) and IPAR(3) are overwritten and */
/*             IPAR(2) <= IPAR(1) and IPAR(3) <= IPAR(1) for all */
/*             examples. */

/*     BPAR    (input) LOGICAL array, dimension (7) */
/*             This array defines the form of the output of the examples */
/*             and the storage mode of the matrices Q, G or R. */
/*             BPAR(1) = .TRUE.  : Q is returned. */
/*             BPAR(1) = .FALSE. : Q is returned in factored form, i.e., */
/*                                 Q0 and C from (I) are returned. */
/*             BPAR(2) = .TRUE.  : The matrix returned in array Q (i.e., */
/*                                 Q if BPAR(1) = .TRUE. and Q0 if */
/*                                 BPAR(1) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(2) = .FALSE. : The matrix returned in array Q is */
/*                                 provided in packed storage mode. */
/*             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array Q is stored in upper */
/*                                 packed mode, i.e., the upper triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 Q(i,j) is stored in the array entry */
/*                                 Q(i+j*(j-1)/2) for i <= j. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array Q is stored in lower */
/*                                 packed mode, i.e., the lower triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 Q(i,j) is stored in the array entry */
/*                                 Q(i+(2*n-j)*(j-1)/2) for j <= i. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(4) = .TRUE.  : The product G in (II) is returned. */
/*             BPAR(4) = .FALSE. : G is returned in factored form, i.e., */
/*                                 B and R from (II) are returned. */
/*             BPAR(5) = .TRUE.  : The matrix returned in array R (i.e., */
/*                                 G if BPAR(4) = .TRUE. and R if */
/*                                 BPAR(4) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(5) = .FALSE. : The matrix returned in array R is */
/*                                 provided in packed storage mode. */
/*             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array R is stored in upper */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array R is stored in lower */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(7) = .TRUE.  : The coefficient matrix S of the DARE */
/*                                 is returned in array S. */
/*             BPAR(7) = .FALSE. : The coefficient matrix S of the DARE */
/*                                 is not returned. */
/*             NOTE that there are no default values for BPAR.  If all */
/*             entries are declared to be .TRUE., then matrices Q, G or R */
/*             are returned in conventional storage mode, i.e., as */
/*             N-by-N or M-by-M arrays where the array element Z(I,J) */
/*             contains the matrix entry Z_{i,j}. */

/*     CHPAR   (output) CHARACTER*255 */
/*             On output, this string contains short information about */
/*             the chosen example. */

/*     VEC     (output) LOGICAL array, dimension (10) */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and */
/*             are always .TRUE. */
/*             VEC(4) refers to A and is always .TRUE. */
/*             VEC(5) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors B */
/*             and R from (II) are returned. */
/*             VEC(6) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors C */
/*             and Q0 from (I) are returned. */
/*             VEC(7) refers to Q and is always .TRUE. */
/*             VEC(8) refers to R and is always .TRUE. */
/*             VEC(9) is .TRUE. if BPAR(7) = .TRUE., i.e., the matrix S */
/*             is returned. */
/*             VEC(10) refers to X and is .TRUE. if the exact solution */
/*             matrix is available. */
/*             NOTE that VEC(i) = .FALSE. for i = 1 to 10 if on exit */
/*             INFO .NE. 0. */

/*     N       (output) INTEGER */
/*             The order of the matrices A, X, G if BPAR(4) = .TRUE., and */
/*             Q if BPAR(1) = .TRUE. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrix B (or the dimension of */
/*             the control input space of the underlying dynamical */
/*             system). */

/*     P       (output) INTEGER */
/*             The number of rows in the matrix C (or the dimension of */
/*             the output space of the underlying dynamical system). */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             coefficient matrix A of the DARE. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If (BPAR(4) = .FALSE.), then the leading N-by-M part */
/*             of this array contains the coefficient matrix B of */
/*             the DARE.  Otherwise, B is used as workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If (BPAR(1) = .FALSE.), then the leading P-by-N part */
/*             of this array contains the matrix C of the factored */
/*             form (I) of Q.  Otherwise, C is used as workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= P. */

/*     Q       (output) DOUBLE PRECISION array, dimension (NQ) */
/*             If (BPAR(1) = .TRUE.) and (BPAR(2) = .TRUE.), then */
/*             NQ = LDQ*N. */
/*             IF (BPAR(1) = .TRUE.) and (BPAR(2) = .FALSE.), then */
/*             NQ = N*(N+1)/2. */
/*             If (BPAR(1) = .FALSE.) and (BPAR(2) = .TRUE.), then */
/*             NQ = LDQ*P. */
/*             IF (BPAR(1) = .FALSE.) and (BPAR(2) = .FALSE.), then */
/*             NQ = P*(P+1)/2. */
/*             The symmetric matrix contained in array Q is stored */
/*             according to BPAR(2) and BPAR(3). */

/*     LDQ     INTEGER */
/*             If conventional storage mode is used for Q, i.e., */
/*             BPAR(2) = .TRUE., then Q is stored like a 2-dimensional */
/*             array with leading dimension LDQ. If packed symmetric */
/*             storage mode is used, then LDQ is irrelevant. */
/*             LDQ >= N if BPAR(1) = .TRUE.; */
/*             LDQ >= P if BPAR(1) = .FALSE.. */

/*     R       (output) DOUBLE PRECISION array, dimension (MR) */
/*             If (BPAR(4) = .TRUE.) and (BPAR(5) = .TRUE.), then */
/*             MR = LDR*N. */
/*             IF (BPAR(4) = .TRUE.) and (BPAR(5) = .FALSE.), then */
/*             MR = N*(N+1)/2. */
/*             If (BPAR(4) = .FALSE.) and (BPAR(5) = .TRUE.), then */
/*             MR = LDR*M. */
/*             IF (BPAR(4) = .FALSE.) and (BPAR(5) = .FALSE.), then */
/*             MR = M*(M+1)/2. */
/*             The symmetric matrix contained in array R is stored */
/*             according to BPAR(5) and BPAR(6). */

/*     LDR     INTEGER */
/*             If conventional storage mode is used for R, i.e., */
/*             BPAR(5) = .TRUE., then R is stored like a 2-dimensional */
/*             array with leading dimension LDR. If packed symmetric */
/*             storage mode is used, then LDR is irrelevant. */
/*             LDR >= N  if BPAR(4) =  .TRUE.; */
/*             LDR >= M  if BPAR(4) = .FALSE.. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,M) */
/*             If (BPAR(7) = .TRUE.), then the leading N-by-M part of */
/*             this array contains the coefficient matrix S of the DARE. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= 1, and */
/*             LDS >= N if BPAR(7) = .TRUE.. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,NX) */
/*             If an exact solution is available (NR = 1.1,1.3,1.4,2.1, */
/*             2.3,2.4,2.5,4.1), then NX = N and the leading N-by-N part */
/*             of this array contains the solution matrix X. */
/*             Otherwise, X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 1, and */
/*             LDX >= N if an exact solution is available. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= N*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0 : successful exit; */
/*             < 0 : if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : data file could not be opened or had wrong format; */
/*             = 2 : division by zero; */
/*             = 3 : G can not be computed as in (II) due to a singular R */
/*                   matrix. This error can only occur if */
/*                   BPAR(4) = .TRUE.. */

/*     REFERENCES */

/*     [1] Abels, J. and Benner, P. */
/*         DAREX - A Collection of Benchmark Examples for Discrete-Time */
/*         Algebraic Riccati Equations (Version 2.0). */
/*         SLICOT Working Note 1999-16, November 1999. Available from */
/*         http://www.win.tue.nl/niconet/NIC2/reports.html. */

/*     This is an updated and extended version of */

/*     [2] Benner, P., Laub, A.J., and Mehrmann, V. */
/*         A Collection of Benchmark Examples for the Numerical Solution */
/*         of Algebraic Riccati Equations II: Discrete-Time Case. */
/*         Technical Report SPC 95_23, Fak. f. Mathematik, */
/*         TU Chemnitz-Zwickau (Germany), December 1995. */

/*     FURTHER COMMENTS */

/*     Some benchmark examples read data from the data files provided */
/*     with the collection. */

/*     CONTRIBUTOR */

/*     Peter Benner (Universitaet Bremen), November 25, 1999. */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please send e-mail to benner@math.uni-bremen.de. */

/*     REVISIONS */

/*     1999, December 23 (V. Sima). */

/*     KEYWORDS */

/*     Discrete-time algebraic Riccati equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     . # of examples available , # of examples with fixed size. . */

/*     .. Scalar Arguments .. */

/*     .. Array Arguments .. */

/*     .. Local Scalars .. */

/*     ..Local Arrays .. */

/*     .. External Functions .. */
/*     . LAPACK . */

/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     . SLICOT . */

/*     .. Intrinsic Functions .. */

/*     .. Data Statements .. */
/*     . default values for dimensions . */
    /* Parameter adjustments */
    --nr;
    --dpar;
    --ipar;
    --bpar;
    --vec;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --q;
    --r__;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --dwork;

    /* Function Body */
/*     . comments on examples . */

/*     .. Executable Statements .. */

    *info = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	vec[i__] = FALSE_;
/* L1: */
    }

    if (nr[1] >= 3) {
	if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
	    ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
	}
	ipar[2] = 1;
	ipar[3] = ipar[1];
    } else {
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
	ipar[3] = pdef[nr[1] + (nr[2] << 1) - 3];
    }

    if (nr[1] >= 2 && ! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def,
	     "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (nr[1] < 1 || nr[1] > 4 || nr[2] < 0 || nr[2] > nex[nr[1] - 1]) 
	    {
	*info = -2;
    } else if (ipar[1] < 1) {
	*info = -4;
    } else if (ipar[1] > *lda) {
	*info = -12;
    } else if (ipar[1] > *ldb) {
	*info = -14;
    } else if (ipar[3] > *ldc) {
	*info = -16;
    } else if (bpar[2] && (! bpar[1] && ipar[3] > *ldq || bpar[1] && ipar[1] 
	    > *ldq)) {
	*info = -18;
    } else if (bpar[5] && (bpar[4] && ipar[1] > *ldr || ! bpar[4] && ipar[2] 
	    > *ldr)) {
	*info = -20;
    } else if (*lds < 1 || bpar[7] && ipar[1] > *lds) {
	*info = -22;
    } else if (*ldx < 1) {
	*info = -24;
    } else if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 3 || nr[2] == 4) || nr[1]
	     == 2 && (nr[2] == 1 || nr[2] >= 3) || nr[1] == 4) {
/*    .. solution X available .. */
	if (ipar[1] > *ldx) {
	    *info = -24;
	} else {
	    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b7, &x[x_offset], ldx, 
		    (ftnlen)1);
	}
    } else if (*ldwork < *n * *n) {
	*info = -26;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("BB02AD", &i__1, (ftnlen)6);
	return 0;
    }

    nsymm = ipar[1] * (ipar[1] + 1) / 2;
    msymm = ipar[2] * (ipar[2] + 1) / 2;
    psymm = ipar[3] * (ipar[3] + 1) / 2;

    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b7, &a[a_offset], lda, (ftnlen)
	    1);
    dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b7, &b[b_offset], ldb, (ftnlen)
	    1);
    dlaset_("A", &ipar[3], &ipar[1], &c_b7, &c_b7, &c__[c_offset], ldc, (
	    ftnlen)1);
    dlaset_("L", &psymm, &c__1, &c_b7, &c_b7, &q[1], &c__1, (ftnlen)1);
    dlaset_("L", &msymm, &c__1, &c_b7, &c_b7, &r__[1], &c__1, (ftnlen)1);
    if (bpar[7]) {
	dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b7, &s[s_offset], lds, (
		ftnlen)1);
    }

    if (nr[1] == 1) {

	if (nr[2] == 1) {
	    a[a_dim1 + 1] = 2.;
	    a[a_dim1 + 2] = 1.;
	    a[(a_dim1 << 1) + 1] = -1.;
	    b[b_dim1 + 1] = 1.;
	    q[1] = 1.;
	    c__[(c_dim1 << 1) + 1] = 1.;
	    r__[1] = 0.;
	    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b34, &x[x_offset], ldx,
		     (ftnlen)1);
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 2) {
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = -1.;
	    b[b_dim1 + 1] = 1.;
	    b[b_dim1 + 2] = 2.;
	    b[(b_dim1 << 1) + 2] = 1.;
	    r__[1] = 9.;
	    r__[2] = 3.;
	    r__[3] = 1.;
	    dlaset_("A", &psymm, &c__1, &c_b38, &c_b38, &q[1], &psymm, (
		    ftnlen)1);
	    q[3] = 7.;
	    drscl_(&msymm, &c_b40, &q[1], &c__1);
	    if (bpar[7]) {
		s[s_dim1 + 1] = 3.;
		s[s_dim1 + 2] = -1.;
		s[(s_dim1 << 1) + 1] = 1.;
		s[(s_dim1 << 1) + 2] = 7.;
	    }
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 3) {
	    a[(a_dim1 << 1) + 1] = 1.;
	    b[b_dim1 + 2] = 1.;
	    q[1] = 1.;
	    q[2] = 2.;
	    q[3] = 4.;
	    x[x_dim1 + 1] = 1.;
	    x[x_dim1 + 2] = 2.;
	    x[(x_dim1 << 1) + 1] = 2.;
	    x[(x_dim1 << 1) + 2] = sqrt(5.) + 2.;
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 4) {
	    a[(a_dim1 << 1) + 1] = .1;
	    a[a_dim1 * 3 + 2] = .01;
	    b[b_dim1 + 1] = 1.;
	    b[(b_dim1 << 1) + 3] = 1.;
	    r__[3] = 1.;
	    q[1] = 1e5;
	    q[4] = 1e3;
	    q[6] = -10.;
	    x[x_dim1 + 1] = 1e5;
	    x[(x_dim1 << 1) + 2] = 1e3;
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] >= 5 && nr[2] <= 8 || nr[2] == 10 || nr[2] == 11 || 
		nr[2] == 13) {
	    if (nr[2] < 10) {
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 11;
		ici__1.iciunit = chpar;
		ici__1.icifmt = "(A,I1,A,I1,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "BB02", (ftnlen)4);
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
		do_fio(&c__1, "0", (ftnlen)1);
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
		do_fio(&c__1, ".dat", (ftnlen)4);
		e_wsfi();
		o__1.oerr = 1;
		o__1.ounit = 1;
		o__1.ofnmlen = 11;
		o__1.ofnm = chpar;
		o__1.orl = 0;
		o__1.osta = "OLD";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		ios = f_open(&o__1);
	    } else {
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 11;
		ici__1.iciunit = chpar;
		ici__1.icifmt = "(A,I1,I2,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "BB02", (ftnlen)4);
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
		do_fio(&c__1, ".dat", (ftnlen)4);
		e_wsfi();
		o__1.oerr = 1;
		o__1.ounit = 1;
		o__1.ofnmlen = 11;
		o__1.ofnm = chpar;
		o__1.orl = 0;
		o__1.osta = "OLD";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		ios = f_open(&o__1);
	    }
	    if (ios != 0) {
		*info = 1;
	    } else {
		if (! (nr[2] == 13)) {
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___12);
			if (ios != 0) {
			    goto L100001;
			}
			i__2 = ipar[1];
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100001;
			    }
			}
			ios = e_rsle();
L100001:
			if (ios != 0) {
			    *info = 1;
			}
/* L10: */
		    }
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___14);
			if (ios != 0) {
			    goto L100002;
			}
			i__2 = ipar[2];
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100002;
			    }
			}
			ios = e_rsle();
L100002:
			if (ios != 0) {
			    *info = 1;
			}
/* L20: */
		    }
		}
		if (nr[2] == 5) {
		    q[1] = 1.87;
		    q[4] = -.244;
		    q[5] = .744;
		    q[6] = .205;
		    q[8] = .589;
		    q[10] = 1.048;
		} else if (nr[2] == 6) {
		    q[1] = .01;
		    q[5] = .01;
		    q[8] = .01;
		    q[10] = .01;
		} else if (nr[2] == 7) {
		    dlaset_("U", &ipar[3], &ipar[1], &c_b34, &c_b34, &c__[
			    c_offset], ldc, (ftnlen)1);
		    c__[c_dim1 * 3 + 1] = 2.;
		    c__[(c_dim1 << 2) + 1] = 4.;
		    c__[(c_dim1 << 2) + 2] = 2.;
		    q[1] = 2.;
		    q[2] = -1.;
		    q[5] = 2.;
		    q[6] = -1.;
		    q[8] = 2.;
		} else if (nr[2] == 10) {
		    c__[c_dim1 + 1] = 1.;
		    c__[c_dim1 * 5 + 2] = 1.;
		    q[1] = 50.;
		    q[3] = 50.;
		} else if (nr[2] == 11) {
		    a[a_dim1 * 10 + 10] = 1.;
		    a[a_dim1 * 11 + 11] = 1.;
		    c__[c_dim1 * 6 + 1] = 15.;
		    c__[c_dim1 * 7 + 2] = 7.;
		    c__[(c_dim1 << 3) + 2] = -5.357;
		    c__[c_dim1 * 9 + 2] = -3.943;
		    c__[c_dim1 * 10 + 3] = 1.;
		    c__[c_dim1 * 11 + 4] = 1.;
		    q[1] = .5;
		    q[5] = 5.;
		    q[8] = .5;
		    q[10] = 5.;
		    r__[1] = 400.;
		    r__[3] = 700.;
		    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

		} else if (nr[2] == 13) {
		    i__1 = ipar[1] - 6;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___15);
			if (ios != 0) {
			    goto L100003;
			}
			i__2 = ipar[1] - 6;
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100003;
			    }
			}
			ios = e_rsle();
L100003:
			if (ios != 0) {
			    *info = 1;
			}
/* L24: */
		    }
		    i__1 = ipar[1] - 6;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___16);
			if (ios != 0) {
			    goto L100004;
			}
			i__2 = ipar[2];
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100004;
			    }
			}
			ios = e_rsle();
L100004:
			if (ios != 0) {
			    *info = 1;
			}
/* L25: */
		    }
		    i__1 = ipar[2];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___17);
			if (ios != 0) {
			    goto L100005;
			}
			i__2 = ipar[1] - 6;
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100005;
			    }
			}
			ios = e_rsle();
L100005:
			if (ios != 0) {
			    *info = 1;
			}
/* L26: */
		    }
		    for (i__ = 1; i__ <= 6; ++i__) {
			a[i__ + 20 + (i__ + 20) * a_dim1] = 1.;
			c__[i__ + 6 + (i__ + 20) * c_dim1] = 1.;
/* L27: */
		    }
		    j = 58;
		    for (i__ = 7; i__ <= 12; ++i__) {
			ios = s_rsle(&io___18);
			if (ios != 0) {
			    goto L100006;
			}
			ios = do_lio(&c__5, &c__1, (char *)&q[j], (ftnlen)
				sizeof(doublereal));
			if (ios != 0) {
			    goto L100006;
			}
			ios = e_rsle();
L100006:
			if (ios != 0) {
			    *info = 1;
			}
			j += 13 - i__;
/* L28: */
		    }
		    j = 1;
		    for (i__ = 1; i__ <= 6; ++i__) {
			ios = s_rsle(&io___19);
			if (ios != 0) {
			    goto L100007;
			}
			ios = do_lio(&c__5, &c__1, (char *)&r__[j], (ftnlen)
				sizeof(doublereal));
			if (ios != 0) {
			    goto L100007;
			}
			ios = e_rsle();
L100007:
			if (ios != 0) {
			    *info = 1;
			}
			j += 7 - i__;
/* L29: */
		    }
		    for (i__ = 1; i__ <= 6; ++i__) {
			for (j = 1; j <= 20; ++j) {
			    a[i__ + 20 + j * a_dim1] = -c__[i__ + j * c_dim1];
/* L30: */
			}
/* L31: */
		    }
		    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    if (nr[2] == 5 || nr[2] == 6) {
		s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
	    } else if (nr[2] == 7 || nr[2] == 10) {
		s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);
	    } else if (nr[2] == 8) {
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
	    }

	} else if (nr[2] == 9) {
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[a_dim1 * 3 + 2] = 1.;
	    a[a_dim1 * 5 + 4] = 1.;
	    a[a_dim1 * 6 + 5] = 1.;
	    b[b_dim1 + 3] = 1.;
	    b[(b_dim1 << 1) + 6] = 1.;
	    c__[c_dim1 + 1] = 1.;
	    c__[(c_dim1 << 1) + 1] = 1.;
	    c__[(c_dim1 << 2) + 2] = 1.;
	    c__[c_dim1 * 5 + 2] = -1.;
	    r__[1] = 3.;
	    r__[3] = 1.;
	    if (bpar[7]) {
		s[s_dim1 + 1] = 1.;
		s[s_dim1 + 2] = 1.;
		s[s_dim1 + 4] = 1.;
		s[s_dim1 + 5] = -1.;
	    }
	    s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);
	} else if (nr[2] == 12) {
	    for (i__ = 1; i__ <= 10; ++i__) {
		a[i__ + (i__ + 1) * a_dim1] = 1.;
/* L32: */
	    }
	    a[a_dim1 * 7 + 6] = 0.;
	    a[a_dim1 * 9 + 8] = 0.;
	    a[a_dim1 * 12 + 12] = 1.;
	    a[a_dim1 * 13 + 13] = 1.;
	    a[a_dim1 + 12] = -3.318;
	    a[a_dim1 + 13] = -1.5484;
	    a[a_dim1 * 6 + 6] = .7788;
	    a[a_dim1 * 7 + 8] = -.4724;
	    a[a_dim1 * 7 + 13] = .3981;
	    a[(a_dim1 << 3) + 8] = 1.3746;
	    a[(a_dim1 << 3) + 13] = .5113;
	    a[a_dim1 * 9 + 13] = 5.7865;
	    a[a_dim1 * 11 + 11] = .8071;
	    b[b_dim1 + 6] = 1.;
	    b[(b_dim1 << 1) + 8] = 1.;
	    c__[c_dim1 + 1] = 3.318;
	    c__[c_dim1 + 2] = 1.5484;
	    c__[c_dim1 * 7 + 2] = -.3981;
	    c__[(c_dim1 << 3) + 2] = -.5113;
	    c__[c_dim1 * 9 + 2] = -5.7865;
	    c__[c_dim1 * 12 + 3] = 1.;
	    c__[c_dim1 * 13 + 4] = 1.;
	    q[1] = .5;
	    q[5] = 5.;
	    q[8] = .5;
	    q[10] = 5.;
	    r__[1] = 400.;
	    r__[3] = 700.;
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
	}

    } else if (nr[1] == 2) {
	if (nr[2] == 1) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e6;
	    }
	    a[a_dim1 + 1] = 4.;
	    a[a_dim1 + 2] = -4.5;
	    a[(a_dim1 << 1) + 1] = 3.;
	    a[(a_dim1 << 1) + 2] = -3.5;
	    dlaset_("A", &ipar[1], &ipar[2], &c_b112, &c_b34, &b[b_offset], 
		    ldb, (ftnlen)1);
	    r__[1] = dpar[1];
	    q[1] = 9.;
	    q[2] = 6.;
	    q[3] = 4.;
	    temp = (sqrt(dpar[1] * 4. + 1.) + 1.) / 2.;
	    x[x_dim1 + 1] = temp * q[1];
	    x[x_dim1 + 2] = temp * q[2];
	    x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
	    x[(x_dim1 << 1) + 2] = temp * q[3];
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 2) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e6;
	    }
	    if (dpar[1] == 0.) {
		*info = 2;
	    } else {
		a[a_dim1 + 1] = .9512;
		a[(a_dim1 << 1) + 2] = .9048;
		dlaset_("A", &c__1, &ipar[2], &c_b118, &c_b118, &b[b_offset], 
			ldb, (ftnlen)1);
		b[b_dim1 + 2] = -1.1895;
		b[(b_dim1 << 1) + 2] = 3.569;
		r__[1] = 1. / (dpar[1] * 3.);
		r__[3] = dpar[1] * 3.;
		q[1] = .005;
		q[3] = .02;
		s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);
	    }

	} else if (nr[2] == 3) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e6;
	    }
	    a[(a_dim1 << 1) + 1] = dpar[1];
	    b[b_dim1 + 2] = 1.;
	    x[x_dim1 + 1] = 1.;
	    x[(x_dim1 << 1) + 2] = dpar[1] * dpar[1] + 1.;
	    s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 4) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e6;
	    }
	    a[(a_dim1 << 1) + 2] = 1.;
	    a[a_dim1 * 3 + 3] = 3.;
	    r__[1] = dpar[1];
	    r__[4] = dpar[1];
	    r__[6] = dpar[1];
/*     .. set C = V .. */
	    temp = .66666666666666663;
	    d__1 = -temp;
	    d__2 = 1. - temp;
	    dlaset_("A", &ipar[3], &ipar[1], &d__1, &d__2, &c__[c_offset], 
		    ldc, (ftnlen)1);
/*     .. and compute A <- C' A C */
	    dsymm_("L", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &a[a_offset], lda, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, 
		    (ftnlen)1);
	    dsymm_("R", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &dwork[1], &ipar[1], &c_b7, &a[a_offset], lda, (ftnlen)1, 
		    (ftnlen)1);
	    q[1] = dpar[1];
	    q[4] = dpar[1];
	    q[6] = dpar[1];
	    x[x_dim1 + 1] = dpar[1];
	    x[(x_dim1 << 1) + 2] = dpar[1] * (sqrt(5.) + 1.) / 2.;
	    x[x_dim1 * 3 + 3] = dpar[1] * (sqrt(85.) + 9.) / 2.;
	    dsymm_("L", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &x[x_offset], ldx, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, 
		    (ftnlen)1);
	    dsymm_("R", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &dwork[1], &ipar[1], &c_b7, &x[x_offset], ldx, (ftnlen)1, 
		    (ftnlen)1);
	    s_copy(ident, "1000", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 5) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[4] = .25;
		dpar[3] = 1.;
		dpar[2] = 1.;
		dpar[1] = 1e8;
	    }
	    if (dpar[1] == 0.) {
		*info = 2;
	    } else {
		temp = dpar[2] / dpar[1];
		beta = dpar[3] * temp;
		alpha = 1. - temp;
		a[a_dim1 + 1] = alpha;
		i__1 = ipar[1] - 1;
		i__2 = ipar[1] - 1;
		dlaset_("A", &i__1, &i__2, &c_b7, &c_b34, &a[a_dim1 + 2], lda,
			 (ftnlen)1);
		b[b_dim1 + 1] = beta;
		c__[(c_dim1 << 2) + 1] = 1.;
		r__[1] = dpar[4];
		if (beta == 0.) {
		    *info = 2;
		} else {
		    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b34, &x[
			    x_offset], ldx, (ftnlen)1);
		    beta *= beta;
		    temp = dpar[4] * (alpha + 1.) * (alpha - 1.) + beta;
		    x[x_dim1 + 1] = temp + sqrt(temp * temp + beta * 4. * 
			    dpar[4]);
		    x[x_dim1 + 1] = x[x_dim1 + 1] / 2. / beta;
		}
		s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);
	    }
	}

    } else if (nr[1] == 4) {
	if (nr[2] == 1) {
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1.;
	    }
	    i__1 = ipar[1] - 1;
	    i__2 = ipar[1] - 1;
	    dlaset_("A", &i__1, &i__2, &c_b7, &c_b34, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
	    b[ipar[1] + b_dim1] = 1.;
	    r__[1] = dpar[1];
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__ + i__ * x_dim1] = (doublereal) i__;
/* L40: */
	    }
	    s_copy(ident, "0110", (ftnlen)4, (ftnlen)4);
	}
    }

    if (*info != 0) {
	goto L2001;
    }
/*     .. set up data in required format .. */

    if (bpar[4]) {
/*     .. G is to be returned in product form .. */
	rdimm = ipar[1];
	if (*(unsigned char *)&ident[3] == '0') {
/*       .. invert R using Cholesky factorization, .. */
	    dpptrf_("L", &ipar[2], &r__[1], info, (ftnlen)1);
	    if (*info == 0) {
		dpptri_("L", &ipar[2], &r__[1], info, (ftnlen)1);
		if (*(unsigned char *)ident == '0') {
/*           .. B is not identity matrix .. */
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dspmv_("L", &ipar[2], &c_b34, &r__[1], &b[i__ + 
				b_dim1], ldb, &c_b7, &dwork[(i__ - 1) * ipar[
				1] + 1], &c__1, (ftnlen)1);
/* L100: */
		    }
		    dgemv_("T", &ipar[2], &ipar[1], &c_b34, &dwork[1], &ipar[
			    1], &b[b_dim1 + 1], ldb, &c_b7, &r__[1], &c__1, (
			    ftnlen)1);
		    isymm = ipar[1] + 1;
		    i__1 = ipar[1];
		    for (i__ = 2; i__ <= i__1; ++i__) {
			dgemv_("T", &ipar[2], &ipar[1], &c_b34, &dwork[1], &
				ipar[1], &b[i__ + b_dim1], ldb, &c_b7, &b[
				b_dim1 + 1], ldb, (ftnlen)1);
			i__2 = ipar[1] - i__ + 1;
			dcopy_(&i__2, &b[i__ * b_dim1 + 1], ldb, &r__[isymm], 
				&c__1);
			isymm += ipar[1] - i__ + 1;
/* L110: */
		    }
		}
	    } else {
		if (*info > 0) {
		    *info = 3;
		    goto L2001;
		}
	    }
	} else {
/*       .. R = identity .. */
	    if (*(unsigned char *)ident == '0') {
/*         .. B not identity matrix .. */
		if (ipar[2] == 1) {
		    dlaset_("L", &nsymm, &c__1, &c_b7, &c_b7, &r__[1], &c__1, 
			    (ftnlen)1);
		    dspr_("L", &ipar[1], &c_b34, &b[b_offset], &c__1, &r__[1],
			     (ftnlen)1);
		} else {
		    dsyrk_("L", "N", &ipar[1], &ipar[2], &c_b34, &b[b_offset],
			     ldb, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, (
			    ftnlen)1);
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    r__[1], (ftnlen)4, (ftnlen)5);
		}
	    } else {
/*         .. B = R = identity .. */
		isymm = 1;
		for (i__ = ipar[1]; i__ >= 1; --i__) {
		    r__[isymm] = 1.;
		    isymm += i__;
/* L120: */
		}
	    }
	}
    } else {
	rdimm = ipar[2];
	if (*(unsigned char *)ident == '1') {
	    dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b34, &b[b_offset], ldb,
		     (ftnlen)1);
	}
	if (*(unsigned char *)&ident[3] == '1') {
	    isymm = 1;
	    for (i__ = ipar[2]; i__ >= 1; --i__) {
		r__[isymm] = 1.;
		isymm += i__;
/* L130: */
	    }
	}
    }

    if (bpar[1]) {
/*     .. Q is to be returned in product form .. */
	qdimm = ipar[1];
	if (*(unsigned char *)&ident[2] == '0') {
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
		i__1 = ipar[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dspmv_("L", &ipar[3], &c_b34, &q[1], &c__[i__ * c_dim1 + 
			    1], &c__1, &c_b7, &dwork[(i__ - 1) * ipar[1] + 1],
			     &c__1, (ftnlen)1);
/* L140: */
		}
/*         .. use Q(1:IPAR(1)) as workspace and compute the first column */
/*            of Q at the end .. */
		isymm = ipar[1] + 1;
		i__1 = ipar[1];
		for (i__ = 2; i__ <= i__1; ++i__) {
		    dgemv_("T", &ipar[3], &ipar[1], &c_b34, &dwork[1], &ipar[
			    1], &c__[i__ * c_dim1 + 1], &c__1, &c_b7, &q[1], &
			    c__1, (ftnlen)1);
		    i__2 = ipar[1] - i__ + 1;
		    dcopy_(&i__2, &q[i__], &c__1, &q[isymm], &c__1);
		    isymm += ipar[1] - i__ + 1;
/* L150: */
		}
		dgemv_("T", &ipar[3], &ipar[1], &c_b34, &dwork[1], &ipar[1], &
			c__[c_dim1 + 1], &c__1, &c_b7, &q[1], &c__1, (ftnlen)
			1);
	    }
	} else {
/*       .. Q = identity .. */
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
		if (ipar[3] == 1) {
		    dlaset_("L", &nsymm, &c__1, &c_b7, &c_b7, &q[1], &c__1, (
			    ftnlen)1);
		    dspr_("L", &ipar[1], &c_b34, &c__[c_offset], ldc, &q[1], (
			    ftnlen)1);
		} else {
		    dsyrk_("L", "T", &ipar[1], &ipar[3], &c_b34, &c__[
			    c_offset], ldc, &c_b7, &dwork[1], &ipar[1], (
			    ftnlen)1, (ftnlen)1);
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    q[1], (ftnlen)4, (ftnlen)5);
		}
	    } else {
/*         .. C = Q = identity .. */
		isymm = 1;
		for (i__ = ipar[1]; i__ >= 1; --i__) {
		    q[isymm] = 1.;
		    isymm += i__;
/* L160: */
		}
	    }
	}
    } else {
	qdimm = ipar[3];
	if (*(unsigned char *)&ident[1] == '1') {
	    dlaset_("A", &ipar[3], &ipar[1], &c_b7, &c_b34, &c__[c_offset], 
		    ldc, (ftnlen)1);
	}
	if (*(unsigned char *)&ident[2] == '1') {
	    isymm = 1;
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
		q[isymm] = 1.;
		isymm += i__;
/* L170: */
	    }
	}
    }

/*     .. unpack symmetric matrices if required .. */
    if (bpar[2]) {
	isymm = qdimm * (qdimm + 1) / 2;
	dcopy_(&isymm, &q[1], &c__1, &dwork[1], &c__1);
	ma02dd_("Unpack", "Lower", &qdimm, &q[1], ldq, &dwork[1], (ftnlen)6, (
		ftnlen)5);
	ma02ed_("Lower", &qdimm, &q[1], ldq, (ftnlen)5);
    } else if (bpar[3]) {
	ma02dd_("Unpack", "Lower", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)
		6, (ftnlen)5);
	ma02ed_("Lower", &qdimm, &dwork[1], &qdimm, (ftnlen)5);
	ma02dd_("Pack", "Upper", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)4, 
		(ftnlen)5);
    }
    if (bpar[5]) {
	isymm = rdimm * (rdimm + 1) / 2;
	dcopy_(&isymm, &r__[1], &c__1, &dwork[1], &c__1);
	ma02dd_("Unpack", "Lower", &rdimm, &r__[1], ldr, &dwork[1], (ftnlen)6,
		 (ftnlen)5);
	ma02ed_("Lower", &rdimm, &r__[1], ldr, (ftnlen)5);
    } else if (bpar[6]) {
	ma02dd_("Unpack", "Lower", &rdimm, &dwork[1], &rdimm, &r__[1], (
		ftnlen)6, (ftnlen)5);
	ma02ed_("Lower", &rdimm, &dwork[1], &rdimm, (ftnlen)5);
	ma02dd_("Pack", "Upper", &rdimm, &dwork[1], &rdimm, &r__[1], (ftnlen)
		4, (ftnlen)5);
    }

/*     ...set VEC... */
    vec[1] = TRUE_;
    vec[2] = TRUE_;
    vec[3] = TRUE_;
    vec[4] = TRUE_;
    vec[5] = ! bpar[4];
    vec[6] = ! bpar[1];
    vec[7] = TRUE_;
    vec[8] = TRUE_;
    vec[9] = bpar[7];
    if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 3 || nr[2] == 4) || nr[1] == 2 
	    && (nr[2] == 1 || nr[2] >= 3) || nr[1] == 4) {
	vec[10] = TRUE_;
    }
    s_copy(chpar, notes + (nr[1] + (nr[2] << 2) - 5) * 255, (ftnlen)255, (
	    ftnlen)255);
    *n = ipar[1];
    *m = ipar[2];
    *p = ipar[3];

L2001:
    return 0;
/* *** Last line of BB02AD *** */
} /* bb02ad_ */

#undef notes


