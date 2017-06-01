/* BB01AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static integer c__1 = 1;
static doublereal c_b30 = 1.;
static doublereal c_b31 = 2.;
static doublereal c_b33 = -1.;
static integer c__5 = 5;

/* Subroutine */ int bb01ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *bpar, char *chpar, logical *vec, integer *n, 
	integer *m, integer *p, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *g, integer *
	ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen def_len, 
	ftnlen chpar_len)
{
    /* Initialized data */

    static integer nex[4] = { 6,9,2,4 };
    static integer ndef[36]	/* was [4][9] */ = { 2,2,20,21,2,2,64,100,4,2,
	    0,30,8,2,0,211,9,2,0,0,30,3,0,0,0,4,0,0,0,4,0,0,0,55 };
    static integer mdef[18]	/* was [2][9] */ = { 1,1,1,2,2,1,2,2,3,1,3,3,
	    0,1,0,1,0,2 };
    static integer pdef[18]	/* was [2][9] */ = { 2,1,2,1,4,2,8,2,9,2,5,3,
	    0,2,0,1,0,10 };
    static doublereal pardef[36]	/* was [4][9] */ = { 0.,1e-6,0.,1.,0.,
	    1e-8,0.,.01,0.,1e6,0.0,4.,0.,1e-7,0.0,0.,0.,0.,0.0,0.0,0.,1e6,0.0,
	    0.0,0.0,1e-6,0.0,0.0,0.0,1e-6,0.0,0.0,0.0,1. };
    static struct {
	char e_1[2550];
	char fill_2[255];
	char e_3[765];
	char fill_4[255];
	char e_5[765];
	char fill_6[510];
	char e_7[510];
	char fill_8[765];
	char e_9[255];
	char fill_10[765];
	char e_11[255];
	char fill_12[765];
	char e_13[255];
	char fill_14[510];
	} equiv_41 = { "Laub 1979, Ex.1                                     "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Arnold/Laub 1984, Ex.1: (A,B) unstabi"
		"lizable as EPS -> 0                                         "
		"                                                            "
		"                                                            "
		"                                      Laub 1979, Ex.4: strin"
		"g of high speed vehicles                                    "
		"                                                            "
		"                                                            "
		"                                                     Laub 19"
		"79, Ex.6: ill-conditioned Riccati equation                  "
		"                                                            "
		"                                                            "
		"                                                            "
		"        Laub 1979, Ex.2: uncontrollable-unobservable data   "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Arnold/Laub 1984, Ex.3: control weigh"
		"ting matrix singular as EPS -> 0                            "
		"                                                            "
		"                                                            "
		"                                      Laub 1979, Ex.5: circu"
		"lant matrices                                               "
		"                                                            "
		"                                                            "
		"                                                     Rosen/W"
		"ang 1992: lq control of 1-dimensional heat flow             "
		"                                                            "
		"                                                            "
		"                                                            "
		"        Beale/Shafai 1989: model of L-1011 aircraft         "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Kenney/Laub/Wette 1989, Ex.2: ARE ill"
		" conditioned for EPS -> oo                                  "
		"                                                            "
		"                                                            "
		"                                      ", {0}, "Hench et al. "
		"1995: coupled springs, dashpots and masses                  "
		"                                                            "
		"                                                            "
		"                                                            "
		"  Bhattacharyya et al. 1983: binary distillation column     "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 Bai/Qian 1994: ill-conditioned Hamiltonian "
		"for EPS -> 0                                                "
		"                                                            "
		"                                                            "
		"                                ", {0}, "Lang/Penzl 1994: ro"
		"tating axle                                                 "
		"                                                            "
		"                                                            "
		"                                                        Patn"
		"aik et al. 1980: tubular ammonia reactor                    "
		"                                                            "
		"                                                            "
		"                                                            "
		"           Laub 1992: H-infinity problem, eigenvalues  +/- E"
		"PS +/- i                                                    "
		"                                                            "
		"                                                            "
		"                          ", {0}, "Davison/Gesing 1978: J-10"
		"0 jet engine                                                "
		"                                                            "
		"                                                            "
		"                                                  Petkov et "
		"al. 1987: increasingly badly scaled Hamiltonian as EPS -> oo"
		"                                                            "
		"                                                            "
		"                                                            "
		"     ", {0}, "Chow/Kokotovic 1976: magnetic tape control sys"
		"tem                                                         "
		"                                                            "
		"                                                            "
		"                             ", {0}, "Arnold/Laub 1984, Ex.2"
		": poor sep. of closed-loop spectrum as EPS -> 0             "
		"                                                            "
		"                                                            "
		"                                                     ", {0}, 
		"IFAC Benchmark Problem #90-06: LQG design for modified Boin"
		"g B-767 at flutter condition                                "
		"                                                            "
		"                                                            "
		"                " };

#define notes ((char *)&equiv_41)


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, x_dim1, 
	    x_offset, i__1, i__2;
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
    double cos(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal b1, b2, c1, c2;
    static integer ios, pos;
    static doublereal sum;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), ma02dd_(char *, 
	    char *, integer *, doublereal *, integer *, doublereal *, ftnlen, 
	    ftnlen), ma02ed_(char *, integer *, doublereal *, integer *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    static integer gdimm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static char ident[4];
    static integer qdimm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dspmv_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen), dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer isymm, msymm, nsymm, psymm;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal appind;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpptrf_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpptri_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpttrf_(integer *, doublereal *,
	     doublereal *, integer *), dpttrs_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };
    static cilist io___20 = { 1, 1, 1, 0, 0 };
    static cilist io___22 = { 1, 1, 1, 0, 0 };
    static cilist io___24 = { 1, 1, 1, 0, 0 };
    static cilist io___25 = { 1, 1, 1, 0, 0 };
    static cilist io___26 = { 1, 1, 1, 0, 0 };
    static cilist io___27 = { 1, 1, 1, 0, 0 };
    static cilist io___28 = { 1, 1, 1, 0, 0 };
    static cilist io___29 = { 1, 1, 1, 0, 0 };
    static cilist io___30 = { 1, 1, 1, 0, 0 };
    static cilist io___37 = { 1, 1, 1, 0, 0 };



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
/*     continuous-time algebraic Riccati equations (CAREs) of the form */

/*       0 = Q + A'X + XA - XGX */

/*     corresponding to the Hamiltonian matrix */

/*            (  A   G  ) */
/*        H = (       T ). */
/*            (  Q  -A  ) */

/*     A,G,Q,X are real N-by-N matrices, Q and G are symmetric and may */
/*     be given in factored form */

/*                   -1 T                         T */
/*      (I)   G = B R  B  ,           (II)   Q = C W C . */

/*     Here, C is P-by-N, W P-by-P, B N-by-M, and R M-by-M, where W */
/*     and R are symmetric. In linear-quadratic optimal control problems, */
/*     usually W is positive semidefinite and R positive definite.  The */
/*     factorized form can be used if the CARE is solved using the */
/*     deflating subspaces of the extended Hamiltonian pencil */

/*                  (  A   0   B  )       (  I   0   0  ) */
/*                  (       T     )       (             ) */
/*        H - s K = (  Q   A   0  )  -  s (  0  -I   0  ) , */
/*                  (       T     )       (             ) */
/*                  (  0   B   R  )       (  0   0   0  ) */

/*     where I and 0 denote the identity and zero matrix, respectively, */
/*     of appropriate dimensions. */

/*     NOTE: the formulation of the CARE and the related matrix (pencils) */
/*           used here does not include CAREs as they arise in robust */
/*           control (H_infinity optimization). */

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
/*             This array determines the example for which CAREX returns */
/*             data. NR(1) is the group of examples. */
/*             NR(1) = 1 : parameter-free problems of fixed size. */
/*             NR(1) = 2 : parameter-dependent problems of fixed size. */
/*             NR(1) = 3 : parameter-free problems of scalable size. */
/*             NR(1) = 4 : parameter-dependent problems of scalable size. */
/*             NR(2) is the number of the example in group NR(1). */
/*             Let NEXi be the number of examples in group i. Currently, */
/*             NEX1 = 6, NEX2 = 9, NEX3 = 2, NEX4 = 4. */
/*             1 <= NR(1) <= 4; */
/*             1 <= NR(2) <= NEXi , where i = NR(1). */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (7) */
/*             Double precision parameter vector. For explanation of the */
/*             parameters see [1]. */
/*             DPAR(1)           : defines the parameters */
/*                                 'delta' for NR(1) = 3, */
/*                                 'q' for NR(1).NR(2) = 4.1, */
/*                                 'a' for NR(1).NR(2) = 4.2, and */
/*                                 'mu' for NR(1).NR(2) = 4.3. */
/*             DPAR(2)           : defines parameters */
/*                                 'r' for NR(1).NR(2) = 4.1, */
/*                                 'b' for NR(1).NR(2) = 4.2, and */
/*                                 'delta' for NR(1).NR(2) = 4.3. */
/*             DPAR(3)           : defines parameters */
/*                                 'c' for NR(1).NR(2) = 4.2 and */
/*                                 'kappa' for NR(1).NR(2) = 4.3. */
/*             DPAR(j), j=4,5,6,7: These arguments are only used to */
/*                                 generate Example 4.2 and define in */
/*                                 consecutive order the intervals */
/*                                 ['beta_1', 'beta_2'], */
/*                                 ['gamma_1', 'gamma_2']. */
/*             NOTE that if DEF = 'D' or 'd', the values of DPAR entries */
/*             on input are ignored and, on output, they are overwritten */
/*             with the default parameters. */

/*     IPAR    (input/output) INTEGER array, dimension (3) */
/*             On input, IPAR(1) determines the actual state dimension, */
/*             i.e., the order of the matrix A as follows, where */
/*             NO = NR(1).NR(2). */
/*             NR(1) = 1 or 2.1-2.8: IPAR(1) is ignored. */
/*             NO = 2.9            : IPAR(1) = 1 generates the CARE for */
/*                                   optimal state feedback (default); */
/*                                   IPAR(1) = 2 generates the Kalman */
/*                                   filter CARE. */
/*             NO = 3.1            : IPAR(1) is the number of vehicles */
/*                                   (parameter 'l' in the description */
/*                                    in [1]). */
/*             NO = 3.2, 4.1 or 4.2: IPAR(1) is the order of the matrix */
/*                                   A. */
/*             NO = 4.3 or 4.4     : IPAR(1) determines the dimension of */
/*                                   the second-order system, i.e., the */
/*                                   order of the stiffness matrix for */
/*                                   Examples 4.3 and 4.4 (parameter 'l' */
/*                                   in the description in [1]). */

/*             The order of the output matrix A is N = 2*IPAR(1) for */
/*             Example 4.3 and N = 2*IPAR(1)-1 for Examples 3.1 and 4.4. */
/*             NOTE that IPAR(1) is overwritten for Examples 1.1-2.8. For */
/*             the other examples, IPAR(1) is overwritten if the default */
/*             parameters are to be used. */
/*             On output, IPAR(1) contains the order of the matrix A. */

/*             On input, IPAR(2) is the number of colums in the matrix B */
/*             in (I) (in control problems, the number of inputs of the */
/*             system). Currently, IPAR(2) is fixed or determined by */
/*             IPAR(1) for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(2) is the number of columns of the */
/*             matrix B from (I). */
/*             NOTE that currently IPAR(2) is overwritten and that */
/*             rank(G) <= IPAR(2). */

/*             On input, IPAR(3) is the number of rows in the matrix C */
/*             in (II) (in control problems, the number of outputs of the */
/*             system). Currently, IPAR(3) is fixed or determined by */
/*             IPAR(1) for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(3) contains the number of rows of the */
/*             matrix C in (II). */
/*             NOTE that currently IPAR(3) is overwritten and that */
/*             rank(Q) <= IPAR(3). */

/*     BPAR    (input) BOOLEAN array, dimension (6) */
/*             This array defines the form of the output of the examples */
/*             and the storage mode of the matrices G and Q. */
/*             BPAR(1) = .TRUE.  : G is returned. */
/*             BPAR(1) = .FALSE. : G is returned in factored form, i.e., */
/*                                 B and R from (I) are returned. */
/*             BPAR(2) = .TRUE.  : The matrix returned in array G (i.e., */
/*                                 G if BPAR(1) = .TRUE. and R if */
/*                                 BPAR(1) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(2) = .FALSE. : The matrix returned in array G is */
/*                                 provided in packed storage mode. */
/*             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array G is stored in upper */
/*                                 packed mode, i.e., the upper triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 G(i,j) is stored in the array entry */
/*                                 G(i+j*(j-1)/2) for i <= j. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array G is stored in lower */
/*                                 packed mode, i.e., the lower triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 G(i,j) is stored in the array entry */
/*                                 G(i+(2*n-j)*(j-1)/2) for j <= i. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(4) = .TRUE.  : Q is returned. */
/*             BPAR(4) = .FALSE. : Q is returned in factored form, i.e., */
/*                                 C and W from (II) are returned. */
/*             BPAR(5) = .TRUE.  : The matrix returned in array Q (i.e., */
/*                                 Q if BPAR(4) = .TRUE. and W if */
/*                                 BPAR(4) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(5) = .FALSE. : The matrix returned in array Q is */
/*                                 provided in packed storage mode. */
/*             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array Q is stored in upper */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array Q is stored in lower */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             NOTE that there are no default values for BPAR.  If all */
/*             entries are declared to be .TRUE., then matrices G and Q */
/*             are returned in conventional storage mode, i.e., as */
/*             N-by-N arrays where the array element Z(I,J) contains the */
/*             matrix entry Z_{i,j}. */

/*     CHPAR   (input/output) CHARACTER*255 */
/*             On input, this is the name of a data file supplied by the */
/*             user. */
/*             In the current version, only Example 4.4 allows a */
/*             user-defined data file. This file must contain */
/*             consecutively DOUBLE PRECISION vectors mu, delta, gamma, */
/*             and kappa. The length of these vectors is determined by */
/*             the input value for IPAR(1). */
/*             If on entry, IPAR(1) = L, then mu and delta must each */
/*             contain L DOUBLE PRECISION values, and gamma and kappa */
/*             must each contain L-1 DOUBLE PRECISION values. */
/*             On output, this string contains short information about */
/*             the chosen example. */

/*     VEC     (output) LOGICAL array, dimension (9) */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and */
/*             are always .TRUE. */
/*             VEC(4) refers to A and is always .TRUE. */
/*             VEC(5) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors B */
/*             and R from (I) are returned. */
/*             VEC(6) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors C */
/*             and W from (II) are returned. */
/*             VEC(7) refers to G and is always .TRUE. */
/*             VEC(8) refers to Q and is always .TRUE. */
/*             VEC(9) refers to X and is .TRUE. if the exact solution */
/*             matrix is available. */
/*             NOTE that VEC(i) = .FALSE. for i = 1 to 9 if on exit */
/*             INFO .NE. 0. */

/*     N       (output) INTEGER */
/*             The order of the matrices A, X, G if BPAR(1) = .TRUE., and */
/*             Q if BPAR(4) = .TRUE. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrix B (or the dimension of */
/*             the control input space of the underlying dynamical */
/*             system). */

/*     P       (output) INTEGER */
/*             The number of rows in the matrix C (or the dimension of */
/*             the output space of the underlying dynamical system). */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             coefficient matrix A of the CARE. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If (BPAR(1) = .FALSE.), then the leading N-by-M part of */
/*             this array contains the matrix B of the factored form (I) */
/*             of G. Otherwise, B is used as workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If (BPAR(4) = .FALSE.), then the leading P-by-N part of */
/*             this array contains the matrix C of the factored form (II) */
/*             of Q. Otherwise, C is used as workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= P, where P is the number of rows of the matrix C, */
/*             i.e., the output value of IPAR(3). (For all examples, */
/*             P <= N, where N equals the output value of the argument */
/*             IPAR(1), i.e., LDC >= LDA is always safe.) */

/*     G       (output) DOUBLE PRECISION array, dimension (NG) */
/*             If (BPAR(2) = .TRUE.)  then NG = LDG*N. */
/*             If (BPAR(2) = .FALSE.) then NG = N*(N+1)/2. */
/*             If (BPAR(1) = .TRUE.), then array G contains the */
/*             coefficient matrix G of the CARE. */
/*             If (BPAR(1) = .FALSE.), then array G contains the 'control */
/*             weighting matrix' R of G's factored form as in (I). (For */
/*             all examples, M <= N.) The symmetric matrix contained in */
/*             array G is stored according to BPAR(2) and BPAR(3). */

/*     LDG     INTEGER */
/*             If conventional storage mode is used for G, i.e., */
/*             BPAR(2) = .TRUE., then G is stored like a 2-dimensional */
/*             array with leading dimension LDG. If packed symmetric */
/*             storage mode is used, then LDG is not referenced. */
/*             LDG >= N if BPAR(2) = .TRUE.. */

/*     Q       (output) DOUBLE PRECISION array, dimension (NQ) */
/*             If (BPAR(5) = .TRUE.)  then NQ = LDQ*N. */
/*             If (BPAR(5) = .FALSE.) then NQ = N*(N+1)/2. */
/*             If (BPAR(4) = .TRUE.), then array Q contains the */
/*             coefficient matrix Q of the CARE. */
/*             If (BPAR(4) = .FALSE.), then array Q contains the 'output */
/*             weighting matrix' W of Q's factored form as in (II). */
/*             The symmetric matrix contained in array Q is stored */
/*             according to BPAR(5) and BPAR(6). */

/*     LDQ     INTEGER */
/*             If conventional storage mode is used for Q, i.e., */
/*             BPAR(5) = .TRUE., then Q is stored like a 2-dimensional */
/*             array with leading dimension LDQ. If packed symmetric */
/*             storage mode is used, then LDQ is not referenced. */
/*             LDQ >= N if BPAR(5) = .TRUE.. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,IPAR(1)) */
/*             If an exact solution is available (NR = 1.1, 1.2, 2.1, */
/*             2.3-2.6, 3.2), then the leading N-by-N part of this array */
/*             contains the solution matrix X in conventional storage */
/*             mode. Otherwise, X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 1, and */
/*             LDX >= N if NR = 1.1, 1.2, 2.1, 2.3-2.6, 3.2. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= N*MAX(4,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0 : successful exit; */
/*             < 0 : if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : data file could not be opened or had wrong format; */
/*             = 2 : division by zero; */
/*             = 3 : G can not be computed as in (I) due to a singular R */
/*                   matrix. */

/*     REFERENCES */

/*     [1] Abels, J. and Benner, P. */
/*         CAREX - A Collection of Benchmark Examples for Continuous-Time */
/*         Algebraic Riccati Equations (Version 2.0). */
/*         SLICOT Working Note 1999-14, November 1999. Available from */
/*         http://www.win.tue.nl/niconet/NIC2/reports.html. */

/*     This is an updated and extended version of */

/*     [2] Benner, P., Laub, A.J., and Mehrmann, V. */
/*         A Collection of Benchmark Examples for the Numerical Solution */
/*         of Algebraic Riccati Equations I: Continuous-Time Case. */
/*         Technical Report SPC 95_22, Fak. f. Mathematik, */
/*         TU Chemnitz-Zwickau (Germany), October 1995. */

/*     NUMERICAL ASPECTS */

/*     If the original data as taken from the literature is given via */
/*     matrices G and Q, but factored forms are requested as output, then */
/*     these factors are obtained from Cholesky or LDL' decompositions of */
/*     G and Q, i.e., the output data will be corrupted by roundoff */
/*     errors. */

/*     FURTHER COMMENTS */

/*     Some benchmark examples read data from the data files provided */
/*     with the collection. */

/*     CONTRIBUTOR */

/*     Peter Benner (Universitaet Bremen), November 15, 1999. */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please send e-mail to benner@math.uni-bremen.de. */

/*     REVISIONS */

/*     1999, December 23 (V. Sima). */

/*     KEYWORDS */

/*     Algebraic Riccati equation, Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     . # of examples available , # of examples with fixed size. . */

/*     .. Scalar Arguments .. */

/*     .. Array Arguments .. */

/*     .. Local Scalars .. */

/*     ..Local Arrays .. */

/*     .. External Functions .. */
/*     . BLAS . */
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
    --g;
    --q;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --dwork;

    /* Function Body */
/*     . default values for parameters . */
/*     . comments on examples . */

/*     .. Executable Statements .. */

    *info = 0;
    for (i__ = 1; i__ <= 9; ++i__) {
	vec[i__] = FALSE_;
/* L5: */
    }

    if (nr[1] != 1 && ! (lsame_(def, "N", (ftnlen)1, (ftnlen)1) || lsame_(def,
	     "D", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (nr[1] < 1 || nr[2] < 1 || nr[1] > 4 || nr[2] > nex[nr[1] - 1]) 
	    {
	*info = -2;
    } else if (nr[1] > 2) {
	if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
	    ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
	}
	if (nr[1] == 3) {
	    if (nr[2] == 1) {
		ipar[2] = ipar[1];
		ipar[3] = ipar[1] - 1;
		ipar[1] = (ipar[1] << 1) - 1;
	    } else if (nr[2] == 2) {
		ipar[2] = ipar[1];
		ipar[3] = ipar[1];
	    } else {
		ipar[2] = 1;
		ipar[3] = 1;
	    }
	} else if (nr[1] == 4) {
	    if (nr[2] == 3) {
		l = ipar[1];
		ipar[2] = 2;
		ipar[3] = l << 1;
		ipar[1] = l << 1;
	    } else if (nr[2] == 4) {
		l = ipar[1];
		ipar[2] = l;
		ipar[3] = l;
		ipar[1] = (l << 1) - 1;
	    } else {
		ipar[2] = 1;
		ipar[3] = 1;
	    }
	}
    } else if (nr[1] == 2 && nr[2] == 9 && ipar[1] == 2) {
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
	ipar[3] = 3;
    } else {
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
	ipar[3] = pdef[nr[1] + (nr[2] << 1) - 3];
    }
    if (*info != 0) {
	goto L7;
    }

    if (ipar[1] < 1) {
	*info = -4;
    } else if (ipar[1] > *lda) {
	*info = -12;
    } else if (ipar[1] > *ldb) {
	*info = -14;
    } else if (ipar[3] > *ldc) {
	*info = -16;
    } else if (bpar[2] && ipar[1] > *ldg) {
	*info = -18;
    } else if (bpar[5] && ipar[1] > *ldq) {
	*info = -20;
    } else if (*ldx < 1) {
	*info = -22;
    } else if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 2)) {
	if (ipar[1] > *ldx) {
	    *info = -22;
	}
    } else if (nr[1] == 2 && nr[2] == 1) {
	if (ipar[1] > *ldx) {
	    *info = -22;
	}
    } else if (nr[1] == 2 && (nr[2] >= 3 && nr[2] <= 6)) {
	if (ipar[1] > *ldx) {
	    *info = -22;
	}
    } else if (nr[1] == 3 && nr[2] == 2) {
	if (ipar[1] > *ldx) {
	    *info = -22;
	}
    } else if (*ldwork < *n * max(4,*n)) {
	*info = -24;
    }

L7:
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("BB01AD", &i__1, (ftnlen)6);
	return 0;
    }

    nsymm = ipar[1] * (ipar[1] + 1) / 2;
    msymm = ipar[2] * (ipar[2] + 1) / 2;
    psymm = ipar[3] * (ipar[3] + 1) / 2;
    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
	dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
    }

    dlaset_("A", &ipar[1], &ipar[1], &c_b10, &c_b10, &a[a_offset], lda, (
	    ftnlen)1);
    dlaset_("A", &ipar[1], &ipar[2], &c_b10, &c_b10, &b[b_offset], ldb, (
	    ftnlen)1);
    dlaset_("A", &ipar[3], &ipar[1], &c_b10, &c_b10, &c__[c_offset], ldc, (
	    ftnlen)1);
    dlaset_("L", &msymm, &c__1, &c_b10, &c_b10, &g[1], &c__1, (ftnlen)1);
    dlaset_("L", &psymm, &c__1, &c_b10, &c_b10, &q[1], &c__1, (ftnlen)1);

    if (nr[1] == 1) {
	if (nr[2] == 1) {
	    a[(a_dim1 << 1) + 1] = 1.;
	    b[b_dim1 + 2] = 1.;
	    q[1] = 1.;
	    q[3] = 2.;
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &c_b31, &x[x_offset], 
		    ldx, (ftnlen)1);

	} else if (nr[2] == 2) {
	    a[a_dim1 + 1] = 4.;
	    a[a_dim1 + 2] = -4.5;
	    a[(a_dim1 << 1) + 1] = 3.;
	    a[(a_dim1 << 1) + 2] = -3.5;
	    dlaset_("A", &ipar[1], &ipar[2], &c_b33, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
	    q[1] = 9.;
	    q[2] = 6.;
	    q[3] = 4.;
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
	    temp = sqrt(2.) + 1.;
	    d__1 = temp * 6.;
	    d__2 = temp * 4.;
	    dlaset_("A", &ipar[1], &ipar[1], &d__1, &d__2, &x[x_offset], ldx, 
		    (ftnlen)1);
	    x[x_dim1 + 1] = temp * 9.;

	} else if (nr[2] >= 3 && nr[2] <= 6) {
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 11;
	    ici__1.iciunit = chpar;
	    ici__1.icifmt = "(A,I1,A,I1,A)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, "BB01", (ftnlen)4);
	    do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, "0", (ftnlen)1);
	    do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".dat", (ftnlen)4);
	    e_wsfi();
	    if (nr[2] == 3 || nr[2] == 4) {
		s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
	    } else if (nr[2] == 5) {
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
	    } else if (nr[2] == 6) {
		s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
	    }
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
	    if (ios != 0) {
		*info = 1;
	    } else if (nr[2] <= 6) {
		i__1 = ipar[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ios = s_rsle(&io___15);
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
		    ios = s_rsle(&io___17);
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
		if (nr[2] <= 4) {
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			pos = (i__ - 1) * ipar[1];
			ios = s_rsle(&io___19);
			if (ios != 0) {
			    goto L100003;
			}
			i__2 = ipar[1];
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&dwork[pos + j]
				    , (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100003;
			    }
			}
			ios = e_rsle();
L100003:
/* L30: */
			;
		    }
		    if (ios != 0) {
			*info = 1;
		    } else {
			ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1]
				, &q[1], (ftnlen)4, (ftnlen)5);
		    }
		} else if (nr[2] == 6) {
		    i__1 = ipar[3];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___20);
			if (ios != 0) {
			    goto L100004;
			}
			i__2 = ipar[1];
			for (j = 1; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100004;
			    }
			}
			ios = e_rsle();
L100004:
			if (ios != 0) {
			    *info = 1;
			}
/* L35: */
		    }
		}
		cl__1.cerr = 0;
		cl__1.cunit = 1;
		cl__1.csta = 0;
		f_clos(&cl__1);
	    }
	}

    } else if (nr[1] == 2) {
	if (nr[2] == 1) {
	    a[a_dim1 + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = -2.;
	    b[b_dim1 + 1] = dpar[1];
	    dlaset_("U", &ipar[3], &ipar[1], &c_b30, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
	    if (dpar[1] != 0.) {
		temp = dlapy2_(&c_b30, &dpar[1]);
		x[x_dim1 + 1] = (temp + 1.) / dpar[1] / dpar[1];
		x[x_dim1 + 2] = 1. / (temp + 2.);
		x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
		ttemp = dpar[1] * x[(x_dim1 << 1) + 1];
		temp = (1. - ttemp) * (ttemp + 1.);
		x[(x_dim1 << 1) + 2] = temp / 4.;
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 2) {
	    a[a_dim1 + 1] = -.1;
	    a[(a_dim1 << 1) + 2] = -.02;
	    b[b_dim1 + 1] = .1;
	    b[b_dim1 + 2] = .001;
	    b[(b_dim1 << 1) + 2] = .01;
	    dlaset_("L", &msymm, &c__1, &c_b30, &c_b30, &g[1], &msymm, (
		    ftnlen)1);
	    g[1] += dpar[1];
	    c__[c_dim1 + 1] = 10.;
	    c__[(c_dim1 << 1) + 1] = 100.;
	    s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 3) {
	    a[(a_dim1 << 1) + 1] = dpar[1];
	    b[b_dim1 + 2] = 1.;
	    s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
	    if (dpar[1] != 0.) {
		temp = sqrt(dpar[1] * 2. + 1.);
		dlaset_("A", &ipar[1], &ipar[1], &c_b30, &temp, &x[x_offset], 
			ldx, (ftnlen)1);
		x[x_dim1 + 1] /= dpar[1];
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 4) {
	    temp = dpar[1] + 1.;
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &temp, &a[a_offset], lda,
		     (ftnlen)1);
/* Computing 2nd power */
	    d__1 = dpar[1];
	    q[1] = d__1 * d__1;
	    q[3] = q[1];
	    s_copy(ident, "1101", (ftnlen)4, (ftnlen)4);
/* Computing 2nd power */
	    d__1 = temp;
	    x[x_dim1 + 1] = temp * 2. + sqrt(2.) * (sqrt(d__1 * d__1 + 1.) + 
		    dpar[1]);
	    x[x_dim1 + 1] /= 2.;
	    x[(x_dim1 << 1) + 2] = x[x_dim1 + 1];
	    ttemp = x[x_dim1 + 1] - temp;
	    if (ttemp != 0.) {
		x[x_dim1 + 2] = x[x_dim1 + 1] / ttemp;
		x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 5) {
	    a[a_dim1 + 1] = 3. - dpar[1];
	    a[a_dim1 + 2] = 4.;
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = 2. - dpar[1];
	    dlaset_("L", &ipar[1], &ipar[2], &c_b30, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
	    q[1] = dpar[1] * 4. - 11.;
	    q[2] = dpar[1] * 2. - 5.;
	    q[3] = dpar[1] * 2. - 2.;
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &c_b30, &x[x_offset], 
		    ldx, (ftnlen)1);
	    x[x_dim1 + 1] = 2.;

	} else if (nr[2] == 6) {
	    if (dpar[1] != 0.) {
		a[a_dim1 + 1] = dpar[1];
		a[(a_dim1 << 1) + 2] = dpar[1] * 2.;
		a[a_dim1 * 3 + 3] = dpar[1] * 3.;
/*     .. set C = V .. */
		temp = .66666666666666663;
		d__1 = -temp;
		d__2 = 1. - temp;
		dlaset_("A", &ipar[3], &ipar[1], &d__1, &d__2, &c__[c_offset],
			 ldc, (ftnlen)1);
		dsymm_("L", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &a[a_offset], lda, &c_b10, &dwork[1], &ipar[1], (
			ftnlen)1, (ftnlen)1);
		dsymm_("R", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &dwork[1], &ipar[1], &c_b10, &a[a_offset], lda, (
			ftnlen)1, (ftnlen)1);
/*     .. G = R ! .. */
		g[1] = dpar[1];
		g[4] = dpar[1];
		g[6] = dpar[1];
		q[1] = 1. / dpar[1];
		q[4] = 1.;
		q[6] = dpar[1];
		s_copy(ident, "1000", (ftnlen)4, (ftnlen)4);
		dlaset_("A", &ipar[1], &ipar[1], &c_b10, &c_b10, &x[x_offset],
			 ldx, (ftnlen)1);
/* Computing 2nd power */
		d__1 = dpar[1];
		temp = d__1 * d__1;
/* Computing 2nd power */
		d__1 = temp;
		x[x_dim1 + 1] = temp + sqrt(d__1 * d__1 + 1.);
/* Computing 2nd power */
		d__1 = temp;
		x[(x_dim1 << 1) + 2] = temp * 2. + sqrt(d__1 * d__1 * 4. + 
			dpar[1]);
		x[x_dim1 * 3 + 3] = temp * 3. + dpar[1] * sqrt(temp * 9. + 1.)
			;
		dsymm_("L", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &x[x_offset], ldx, &c_b10, &dwork[1], &ipar[1], (
			ftnlen)1, (ftnlen)1);
		dsymm_("R", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &dwork[1], &ipar[1], &c_b10, &x[x_offset], ldx, (
			ftnlen)1, (ftnlen)1);
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 7) {
	    if (dpar[1] != 0.) {
		a[(a_dim1 << 1) + 1] = .4;
		a[a_dim1 * 3 + 2] = .345;
		a[(a_dim1 << 1) + 3] = -.524 / dpar[1];
		a[a_dim1 * 3 + 3] = -.465 / dpar[1];
		a[(a_dim1 << 2) + 3] = .262 / dpar[1];
		a[(a_dim1 << 2) + 4] = -1. / dpar[1];
		b[b_dim1 + 4] = 1. / dpar[1];
		c__[c_dim1 + 1] = 1.;
		c__[c_dim1 * 3 + 2] = 1.;
		s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 8) {
	    a[a_dim1 + 1] = -dpar[1];
	    a[a_dim1 + 2] = -1.;
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = -dpar[1];
	    a[a_dim1 * 3 + 3] = dpar[1];
	    a[a_dim1 * 3 + 4] = -1.;
	    a[(a_dim1 << 2) + 3] = 1.;
	    a[(a_dim1 << 2) + 4] = dpar[1];
	    dlaset_("L", &ipar[1], &ipar[2], &c_b30, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
	    dlaset_("U", &ipar[3], &ipar[1], &c_b30, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 9) {
	    if (ipar[3] == 10) {
/*     .. read LQR CARE ... */
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 12;
		ici__1.iciunit = chpar;
		ici__1.icifmt = "(A,I1,A,I1,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "BB01", (ftnlen)4);
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
		do_fio(&c__1, "0", (ftnlen)1);
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
		do_fio(&c__1, "1.dat", (ftnlen)5);
		e_wsfi();
		o__1.oerr = 1;
		o__1.ounit = 1;
		o__1.ofnmlen = 12;
		o__1.ofnm = chpar;
		o__1.orl = 0;
		o__1.osta = "OLD";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		ios = f_open(&o__1);
		if (ios != 0) {
		    *info = 1;
		} else {
		    for (i__ = 1; i__ <= 27; i__ += 2) {
			ios = s_rsle(&io___22);
			if (ios != 0) {
			    goto L100005;
			}
			for (j = 0; j <= 1; ++j) {
			    for (k = 0; k <= 1; ++k) {
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j 
					+ (i__ + k) * a_dim1], (ftnlen)sizeof(
					doublereal));
				if (ios != 0) {
				    goto L100005;
				}
			    }
			}
			ios = e_rsle();
L100005:
			if (ios != 0) {
			    *info = 1;
			}
/* L36: */
		    }
		    for (i__ = 30; i__ <= 44; i__ += 2) {
			ios = s_rsle(&io___24);
			if (ios != 0) {
			    goto L100006;
			}
			for (j = 0; j <= 1; ++j) {
			    for (k = 0; k <= 1; ++k) {
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j 
					+ (i__ + k) * a_dim1], (ftnlen)sizeof(
					doublereal));
				if (ios != 0) {
				    goto L100006;
				}
			    }
			}
			ios = e_rsle();
L100006:
			if (ios != 0) {
			    *info = 1;
			}
/* L37: */
		    }
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___25);
			if (ios != 0) {
			    goto L100007;
			}
			i__2 = ipar[1];
			for (j = 46; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100007;
			    }
			}
			ios = e_rsle();
L100007:
			if (ios != 0) {
			    *info = 1;
			}
/* L38: */
		    }
		    a[a_dim1 * 29 + 29] = -5.301;
		    b[b_dim1 + 48] = 8e5;
		    b[(b_dim1 << 1) + 51] = 8e5;
		    g[1] = 364.7;
		    g[3] = 14.59;
		    for (i__ = 1; i__ <= 6; ++i__) {
			ios = s_rsle(&io___26);
			if (ios != 0) {
			    goto L100008;
			}
			for (j = 1; j <= 45; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100008;
			    }
			}
			ios = e_rsle();
L100008:
			if (ios != 0) {
			    *info = 1;
			}
/* L39: */
		    }
		    c__[c_dim1 * 47 + 7] = 1.;
		    c__[c_dim1 * 46 + 8] = 1.;
		    c__[c_dim1 * 50 + 9] = 1.;
		    c__[c_dim1 * 49 + 10] = 1.;
		    q[11] = 3.76e-14;
		    q[20] = 1.2e-13;
		    q[41] = 2.45e-12;
		}
	    } else {
/*     .. read Kalman filter CARE .. */
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 12;
		ici__1.iciunit = chpar;
		ici__1.icifmt = "(A,I1,A,I1,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "BB01", (ftnlen)4);
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
		do_fio(&c__1, "0", (ftnlen)1);
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
		do_fio(&c__1, "2.dat", (ftnlen)5);
		e_wsfi();
		o__1.oerr = 1;
		o__1.ounit = 1;
		o__1.ofnmlen = 12;
		o__1.ofnm = chpar;
		o__1.orl = 0;
		o__1.osta = "OLD";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		ios = f_open(&o__1);
		if (ios != 0) {
		    *info = 1;
		} else {
		    for (i__ = 1; i__ <= 27; i__ += 2) {
			ios = s_rsle(&io___27);
			if (ios != 0) {
			    goto L100009;
			}
			for (j = 0; j <= 1; ++j) {
			    for (k = 0; k <= 1; ++k) {
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + k 
					+ (i__ + j) * a_dim1], (ftnlen)sizeof(
					doublereal));
				if (ios != 0) {
				    goto L100009;
				}
			    }
			}
			ios = e_rsle();
L100009:
			if (ios != 0) {
			    *info = 1;
			}
/* L40: */
		    }
		    for (i__ = 30; i__ <= 44; i__ += 2) {
			ios = s_rsle(&io___28);
			if (ios != 0) {
			    goto L100010;
			}
			for (j = 0; j <= 1; ++j) {
			    for (k = 0; k <= 1; ++k) {
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + k 
					+ (i__ + j) * a_dim1], (ftnlen)sizeof(
					doublereal));
				if (ios != 0) {
				    goto L100010;
				}
			    }
			}
			ios = e_rsle();
L100010:
			if (ios != 0) {
			    *info = 1;
			}
/* L41: */
		    }
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			ios = s_rsle(&io___29);
			if (ios != 0) {
			    goto L100011;
			}
			i__2 = ipar[1];
			for (j = 46; j <= i__2; ++j) {
			    ios = do_lio(&c__5, &c__1, (char *)&a[j + i__ * 
				    a_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100011;
			    }
			}
			ios = e_rsle();
L100011:
			if (ios != 0) {
			    *info = 1;
			}
/* L42: */
		    }
		    a[a_dim1 * 29 + 29] = -5.301;
		    i__1 = ipar[2];
		    for (j = 1; j <= i__1; ++j) {
			ios = s_rsle(&io___30);
			if (ios != 0) {
			    goto L100012;
			}
			i__2 = ipar[1];
			for (i__ = 1; i__ <= i__2; ++i__) {
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
			    if (ios != 0) {
				goto L100012;
			    }
			}
			ios = e_rsle();
L100012:
			if (ios != 0) {
			    *info = 1;
			}
/* L43: */
		    }
		    g[1] = 6.85e-6;
		    g[3] = 373.;
		    c__[c_dim1 * 52 + 1] = .3713f;
		    c__[c_dim1 * 53 + 1] = 1.245;
		    c__[c_dim1 * 48 + 2] = 8e5;
		    c__[c_dim1 * 54 + 2] = 1.;
		    c__[c_dim1 * 51 + 3] = 8e5;
		    c__[c_dim1 * 55 + 3] = 1.;
		    q[1] = 28224.;
		    q[4] = 2.742e-5;
		    q[6] = 6.854e-4;
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
	}

    } else if (nr[1] == 3) {
	if (nr[2] == 1) {
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (i__ % 2 == 1) {
		    a[i__ + i__ * a_dim1] = -1.;
		    b[i__ + (i__ + 1) / 2 * b_dim1] = 1.;
		} else {
		    a[i__ + (i__ - 1) * a_dim1] = 1.;
		    a[i__ + (i__ + 1) * a_dim1] = -1.;
		    c__[i__ / 2 + i__ * c_dim1] = 1.;
		}
/* L45: */
	    }
	    isymm = 1;
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
		q[isymm] = 10.;
		isymm += i__;
/* L50: */
	    }
	    s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 2) {
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + i__ * a_dim1] = -2.;
		if (i__ < ipar[1]) {
		    a[i__ + (i__ + 1) * a_dim1] = 1.;
		    a[i__ + 1 + i__ * a_dim1] = 1.;
		}
/* L60: */
	    }
	    a[ipar[1] * a_dim1 + 1] = 1.;
	    a[ipar[1] + a_dim1] = 1.;
	    s_copy(ident, "1111", (ftnlen)4, (ftnlen)4);
	    temp = 6.2831853071795862 / (doublereal) ipar[1];
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dwork[i__] = cos(temp * (doublereal) (i__ - 1));
		dwork[ipar[1] + i__] = dwork[i__] * 2. - 2. + sqrt(dwork[i__] 
			* 4. * (dwork[i__] - 2.) + 5.);
/* L70: */
	    }
	    i__1 = ipar[1];
	    for (j = 1; j <= i__1; ++j) {
		i__2 = ipar[1];
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dwork[(ipar[1] << 1) + i__] = cos(temp * (doublereal) (
			    i__ - 1) * (doublereal) (j - 1));
/* L80: */
		}
		x[j + x_dim1] = ddot_(&ipar[1], &dwork[ipar[1] + 1], &c__1, &
			dwork[(ipar[1] << 1) + 1], &c__1) / (doublereal) ipar[
			1];
/* L90: */
	    }
/*         .. set up circulant solution matrix .. */
	    i__1 = ipar[1];
	    for (i__ = 2; i__ <= i__1; ++i__) {
		i__2 = ipar[1] - i__ + 1;
		dcopy_(&i__2, &x[x_dim1 + 1], &c__1, &x[i__ + i__ * x_dim1], &
			c__1);
		i__2 = i__ - 1;
		dcopy_(&i__2, &x[ipar[1] - i__ + 2 + x_dim1], &c__1, &x[i__ * 
			x_dim1 + 1], &c__1);
/* L100: */
	    }
	}

    } else if (nr[1] == 4) {
	if (nr[2] == 1) {
/*       .. set up remaining parameter .. */
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1.;
		dpar[2] = 1.;
	    }
	    i__1 = ipar[1] - 1;
	    i__2 = ipar[1] - 1;
	    dlaset_("A", &i__1, &i__2, &c_b10, &c_b30, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
	    b[ipar[1] + b_dim1] = 1.;
	    c__[c_dim1 + 1] = 1.;
	    q[1] = dpar[1];
	    g[1] = dpar[2];
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 2) {
/*         .. set up remaining parameters .. */
	    appind = (doublereal) (ipar[1] + 1);
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
		dpar[2] = 1.;
		dpar[3] = 1.;
		dpar[4] = .2;
		dpar[5] = .3;
		dpar[6] = .2;
		dpar[7] = .3;
	    }
/*         .. set up stiffness matrix .. */
	    temp = -dpar[1] * appind;
	    d__1 = temp * 2.;
	    dlaset_("A", &ipar[1], &ipar[1], &c_b10, &d__1, &a[a_offset], lda,
		     (ftnlen)1);
	    i__1 = ipar[1] - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + 1 + i__ * a_dim1] = -temp;
		a[i__ + (i__ + 1) * a_dim1] = -temp;
/* L110: */
	    }
/*         .. set up Gramian, stored by diagonals .. */
	    temp = 1. / (appind * 6.);
	    d__1 = temp * 4.;
	    d__2 = temp * 4.;
	    dlaset_("L", &ipar[1], &c__1, &d__1, &d__2, &dwork[1], &ipar[1], (
		    ftnlen)1);
	    i__1 = ipar[1] - 1;
	    dlaset_("L", &i__1, &c__1, &temp, &temp, &dwork[ipar[1] + 1], &
		    ipar[1], (ftnlen)1);
	    dpttrf_(&ipar[1], &dwork[1], &dwork[ipar[1] + 1], info);
/*         .. A = (inverse of Gramian) * (stiffness matrix) .. */
	    dpttrs_(&ipar[1], &ipar[1], &dwork[1], &dwork[ipar[1] + 1], &a[
		    a_offset], lda, info);
/*         .. compute B, C .. */
	    i__1 = ipar[1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__1 = (doublereal) (i__ - 1) / appind;
		b1 = max(d__1,dpar[4]);
/* Computing MIN */
		d__1 = (doublereal) (i__ + 1) / appind;
		b2 = min(d__1,dpar[5]);
/* Computing MAX */
		d__1 = (doublereal) (i__ - 1) / appind;
		c1 = max(d__1,dpar[6]);
/* Computing MIN */
		d__1 = (doublereal) (i__ + 1) / appind;
		c2 = min(d__1,dpar[7]);
		if (b1 >= b2) {
		    b[i__ + b_dim1] = 0.;
		} else {
		    b[i__ + b_dim1] = b2 - b1;
/* Computing MIN */
		    d__1 = b2, d__2 = (doublereal) i__ / appind;
		    temp = min(d__1,d__2);
		    if (b1 < temp) {
/* Computing 2nd power */
			d__1 = temp;
/* Computing 2nd power */
			d__2 = b1;
			b[i__ + b_dim1] += appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
			b[i__ + b_dim1] += (doublereal) i__ * (b1 - temp);
		    }
/* Computing MAX */
		    d__1 = b1, d__2 = (doublereal) i__ / appind;
		    temp = max(d__1,d__2);
		    if (temp < b2) {
/* Computing 2nd power */
			d__1 = b2;
/* Computing 2nd power */
			d__2 = temp;
			b[i__ + b_dim1] -= appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
			b[i__ + b_dim1] -= (doublereal) i__ * (temp - b2);
		    }
		}
		if (c1 >= c2) {
		    c__[i__ * c_dim1 + 1] = 0.;
		} else {
		    c__[i__ * c_dim1 + 1] = c2 - c1;
/* Computing MIN */
		    d__1 = c2, d__2 = (doublereal) i__ / appind;
		    temp = min(d__1,d__2);
		    if (c1 < temp) {
/* Computing 2nd power */
			d__1 = temp;
/* Computing 2nd power */
			d__2 = c1;
			c__[i__ * c_dim1 + 1] += appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
			c__[i__ * c_dim1 + 1] += (doublereal) i__ * (c1 - 
				temp);
		    }
/* Computing MAX */
		    d__1 = c1, d__2 = (doublereal) i__ / appind;
		    temp = max(d__1,d__2);
		    if (temp < c2) {
/* Computing 2nd power */
			d__1 = c2;
/* Computing 2nd power */
			d__2 = temp;
			c__[i__ * c_dim1 + 1] -= appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
			c__[i__ * c_dim1 + 1] -= (doublereal) i__ * (temp - 
				c2);
		    }
		}
/* L120: */
	    }
	    dscal_(&ipar[1], &dpar[2], &b[b_dim1 + 1], &c__1);
	    dscal_(&ipar[1], &dpar[3], &c__[c_dim1 + 1], ldc);
	    dpttrs_(&ipar[1], &c__1, &dwork[1], &dwork[ipar[1] + 1], &b[
		    b_offset], ldb, info);
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);

	} else if (nr[2] == 3) {
/*         .. set up remaining parameters .. */
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
		dpar[2] = 4.;
		dpar[3] = 1.;
	    }
	    if (dpar[1] != 0.) {
		dlaset_("A", &l, &l, &c_b10, &c_b30, &a[(l + 1) * a_dim1 + 1],
			 lda, (ftnlen)1);
		temp = dpar[3] / dpar[1];
		a[l + 1 + a_dim1] = -temp;
		a[l + 1 + (a_dim1 << 1)] = temp;
		a[ipar[1] + (l - 1) * a_dim1] = temp;
		a[ipar[1] + l * a_dim1] = -temp;
		ttemp = temp * 2.;
		i__1 = l - 1;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    a[l + i__ + i__ * a_dim1] = -ttemp;
		    a[l + i__ + (i__ + 1) * a_dim1] = temp;
		    a[l + i__ + (i__ - 1) * a_dim1] = temp;
/* L130: */
		}
		d__1 = -dpar[2] / dpar[1];
		dlaset_("A", &l, &l, &c_b10, &d__1, &a[l + 1 + (l + 1) * 
			a_dim1], lda, (ftnlen)1);
		b[l + 1 + b_dim1] = 1. / dpar[1];
		b[ipar[1] + ipar[2] * b_dim1] = -1. / dpar[1];
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
	    } else {
		*info = 2;
	    }

	} else if (nr[2] == 4) {
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
		ici__1.icierr = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 11;
		ici__1.iciunit = chpar;
		ici__1.icifmt = "(A,I1,A,I1,A)";
		s_wsfi(&ici__1);
		do_fio(&c__1, "BB01", (ftnlen)4);
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
		do_fio(&c__1, "0", (ftnlen)1);
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
		do_fio(&c__1, ".dat", (ftnlen)4);
		e_wsfi();
	    }
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
	    if (ios != 0) {
		*info = 1;
	    } else {
		ios = s_rsle(&io___37);
		if (ios != 0) {
		    goto L100013;
		}
		i__1 = (l << 2) - 2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ios = do_lio(&c__5, &c__1, (char *)&dwork[i__], (ftnlen)
			    sizeof(doublereal));
		    if (ios != 0) {
			goto L100013;
		    }
		}
		ios = e_rsle();
L100013:
		if (ios != 0) {
		    *info = 1;
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    if (*info == 0) {
		i__1 = l - 1;
		i__2 = l - 1;
		dlaset_("A", &i__1, &i__2, &c_b10, &c_b30, &a[l + 1 + (a_dim1 
			<< 1)], lda, (ftnlen)1);
		pos = (l << 1) + 1;
		a[(a_dim1 << 1) + 1] = -dwork[pos] / dwork[1];
		i__1 = l;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    temp = dwork[pos] / dwork[i__ - 1];
		    ttemp = dwork[pos] / dwork[i__];
		    if (i__ > 2) {
			a[i__ - 1 + i__ * a_dim1] = temp;
		    }
		    a[i__ + i__ * a_dim1] = -(temp + ttemp);
		    if (i__ < l) {
			a[i__ + 1 + i__ * a_dim1] = ttemp;
		    }
		    ++pos;
/* L140: */
		}
		pos = l;
		temp = dwork[pos + 1] / dwork[1];
		a[a_dim1 + 1] = -temp;
		i__1 = l;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    ttemp = temp;
		    temp = dwork[pos + i__] / dwork[i__];
		    sum = ttemp - temp;
		    a[i__ + a_dim1] = -sum;
		    a[i__ + i__ * a_dim1] -= temp;
		    i__2 = i__ - 2;
		    for (j = 2; j <= i__2; ++j) {
			a[i__ + j * a_dim1] = sum;
/* L150: */
		    }
		    if (i__ > 2) {
			a[i__ + (i__ - 1) * a_dim1] += sum;
		    }
/* L160: */
		}
		pos = l * 3;
		a[(l + 1) * a_dim1 + 1] = -dwork[l * 3] / dwork[1];
		i__1 = l;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    temp = dwork[pos] / dwork[i__ - 1];
		    ttemp = dwork[pos] / dwork[i__];
		    if (i__ > 2) {
			a[i__ - 1 + (l + i__ - 1) * a_dim1] = temp;
		    }
		    a[i__ + (l + i__ - 1) * a_dim1] = -(temp + ttemp);
		    if (i__ < l) {
			a[i__ + 1 + (l + i__ - 1) * a_dim1] = ttemp;
		    }
		    ++pos;
/* L170: */
		}
		b[b_dim1 + 1] = 1. / dwork[1];
		i__1 = l;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    temp = 1. / dwork[i__];
		    if (i__ > 1) {
			b[i__ + i__ * b_dim1] = -temp;
		    }
		    if (i__ < l) {
			b[i__ + 1 + i__ * b_dim1] = temp;
		    }
/* L180: */
		}
		c__[c_dim1 + 1] = 1.;
		q[1] = 1.;
		pos = (l << 1) - 1;
		isymm = l + 1;
		i__1 = l;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    temp = dwork[pos + i__];
		    ttemp = dwork[pos + l + i__ - 1];
		    c__[i__ + i__ * c_dim1] = temp;
		    c__[i__ + (l + i__ - 1) * c_dim1] = ttemp;
		    q[isymm] = 1. / (temp * temp + ttemp * ttemp);
		    isymm = isymm + l - i__ + 1;
/* L190: */
		}
		s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);
	    }
	}
    }

    if (*info != 0) {
	goto L2001;
    }
/*     .. set up data in required format .. */

    if (bpar[1]) {
/*     .. G is to be returned in product form .. */
	gdimm = ipar[1];
	if (*(unsigned char *)&ident[3] == '0') {
/*       .. invert R using Cholesky factorization, store in G .. */
	    dpptrf_("L", &ipar[2], &g[1], info, (ftnlen)1);
	    if (*info == 0) {
		dpptri_("L", &ipar[2], &g[1], info, (ftnlen)1);
		if (*(unsigned char *)ident == '0') {
/*         .. B is not identity matrix .. */
		    i__1 = ipar[1];
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dspmv_("L", &ipar[2], &c_b30, &g[1], &b[i__ + b_dim1],
				 ldb, &c_b10, &dwork[(i__ - 1) * ipar[1] + 1],
				 &c__1, (ftnlen)1);
/* L200: */
		    }
		    dgemv_("T", &ipar[2], &ipar[1], &c_b30, &dwork[1], &ipar[
			    1], &b[b_dim1 + 1], ldb, &c_b10, &g[1], &c__1, (
			    ftnlen)1);
		    isymm = ipar[1] + 1;
		    i__1 = ipar[1];
		    for (i__ = 2; i__ <= i__1; ++i__) {
			dgemv_("T", &ipar[2], &ipar[1], &c_b30, &dwork[1], &
				ipar[1], &b[i__ + b_dim1], ldb, &c_b10, &b[
				b_dim1 + 1], ldb, (ftnlen)1);
			i__2 = ipar[1] - i__ + 1;
			dcopy_(&i__2, &b[i__ * b_dim1 + 1], ldb, &g[isymm], &
				c__1);
			isymm += ipar[1] - i__ + 1;
/* L210: */
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
/*         .. B is not identity matrix .. */
		if (ipar[2] == 1) {
		    dlaset_("L", &nsymm, &c__1, &c_b10, &c_b10, &g[1], &c__1, 
			    (ftnlen)1);
		    dspr_("L", &ipar[1], &c_b30, &b[b_offset], &c__1, &g[1], (
			    ftnlen)1);
		} else {
		    dsyrk_("L", "N", &ipar[1], &ipar[2], &c_b30, &b[b_offset],
			     ldb, &c_b10, &dwork[1], &ipar[1], (ftnlen)1, (
			    ftnlen)1);
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    g[1], (ftnlen)4, (ftnlen)5);
		}
	    } else {
/*         .. B = R = identity .. */
		isymm = 1;
		for (i__ = ipar[1]; i__ >= 1; --i__) {
		    g[isymm] = 1.;
		    isymm += i__;
/* L220: */
		}
	    }
	}
    } else {
	gdimm = ipar[2];
	if (*(unsigned char *)ident == '1') {
	    dlaset_("A", &ipar[1], &ipar[2], &c_b10, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
	}
	if (*(unsigned char *)&ident[3] == '1') {
	    isymm = 1;
	    for (i__ = ipar[2]; i__ >= 1; --i__) {
		g[isymm] = 1.;
		isymm += i__;
/* L230: */
	    }
	}
    }

    if (bpar[4]) {
/*     .. Q is to be returned in product form .. */
	qdimm = ipar[1];
	if (*(unsigned char *)&ident[2] == '0') {
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
		i__1 = ipar[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dspmv_("L", &ipar[3], &c_b30, &q[1], &c__[i__ * c_dim1 + 
			    1], &c__1, &c_b10, &dwork[(i__ - 1) * ipar[1] + 1]
			    , &c__1, (ftnlen)1);
/* L240: */
		}
/*         .. use Q(1:IPAR(1)) as workspace and compute the first column */
/*            of Q in the end .. */
		isymm = ipar[1] + 1;
		i__1 = ipar[1];
		for (i__ = 2; i__ <= i__1; ++i__) {
		    dgemv_("T", &ipar[3], &ipar[1], &c_b30, &dwork[1], &ipar[
			    1], &c__[i__ * c_dim1 + 1], &c__1, &c_b10, &q[1], 
			    &c__1, (ftnlen)1);
		    i__2 = ipar[1] - i__ + 1;
		    dcopy_(&i__2, &q[i__], &c__1, &q[isymm], &c__1);
		    isymm += ipar[1] - i__ + 1;
/* L250: */
		}
		dgemv_("T", &ipar[3], &ipar[1], &c_b30, &dwork[1], &ipar[1], &
			c__[c_dim1 + 1], &c__1, &c_b10, &q[1], &c__1, (ftnlen)
			1);
	    }
	} else {
/*       .. Q = identity .. */
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
		if (ipar[3] == 1) {
		    dlaset_("L", &nsymm, &c__1, &c_b10, &c_b10, &q[1], &c__1, 
			    (ftnlen)1);
		    dspr_("L", &ipar[1], &c_b30, &c__[c_offset], ldc, &q[1], (
			    ftnlen)1);
		} else {
		    dsyrk_("L", "T", &ipar[1], &ipar[3], &c_b30, &c__[
			    c_offset], ldc, &c_b10, &dwork[1], &ipar[1], (
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
/* L260: */
		}
	    }
	}
    } else {
	qdimm = ipar[3];
	if (*(unsigned char *)&ident[1] == '1') {
	    dlaset_("A", &ipar[3], &ipar[1], &c_b10, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
	}
	if (*(unsigned char *)&ident[2] == '1') {
	    isymm = 1;
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
		q[isymm] = 1.;
		isymm += i__;
/* L270: */
	    }
	}
    }

/*     .. unpack symmetric matrices if desired .. */
    if (bpar[2]) {
	isymm = gdimm * (gdimm + 1) / 2;
	dcopy_(&isymm, &g[1], &c__1, &dwork[1], &c__1);
	ma02dd_("Unpack", "Lower", &gdimm, &g[1], ldg, &dwork[1], (ftnlen)6, (
		ftnlen)5);
	ma02ed_("Lower", &gdimm, &g[1], ldg, (ftnlen)5);
    } else if (bpar[3]) {
	ma02dd_("Unpack", "Lower", &gdimm, &dwork[1], &gdimm, &g[1], (ftnlen)
		6, (ftnlen)5);
	ma02ed_("Lower", &gdimm, &dwork[1], &gdimm, (ftnlen)5);
	ma02dd_("Pack", "Upper", &gdimm, &dwork[1], &gdimm, &g[1], (ftnlen)4, 
		(ftnlen)5);
    }
    if (bpar[5]) {
	isymm = qdimm * (qdimm + 1) / 2;
	dcopy_(&isymm, &q[1], &c__1, &dwork[1], &c__1);
	ma02dd_("Unpack", "Lower", &qdimm, &q[1], ldq, &dwork[1], (ftnlen)6, (
		ftnlen)5);
	ma02ed_("Lower", &qdimm, &q[1], ldq, (ftnlen)5);
    } else if (bpar[6]) {
	ma02dd_("Unpack", "Lower", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)
		6, (ftnlen)5);
	ma02ed_("Lower", &qdimm, &dwork[1], &qdimm, (ftnlen)5);
	ma02dd_("Pack", "Upper", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)4, 
		(ftnlen)5);
    }

/*     ...set VEC... */
    vec[1] = TRUE_;
    vec[2] = TRUE_;
    vec[3] = TRUE_;
    vec[4] = TRUE_;
    vec[5] = ! bpar[1];
    vec[6] = ! bpar[4];
    vec[7] = TRUE_;
    vec[8] = TRUE_;
    if (nr[1] == 1) {
	if (nr[2] == 1 || nr[2] == 2) {
	    vec[9] = TRUE_;
	}
    } else if (nr[1] == 2) {
	if (nr[2] == 1 || nr[2] >= 3 && nr[2] <= 6) {
	    vec[9] = TRUE_;
	}
    } else if (nr[1] == 3) {
	if (nr[2] == 2) {
	    vec[9] = TRUE_;
	}
    }
    s_copy(chpar, notes + (nr[1] + (nr[2] << 2) - 5) * 255, (ftnlen)255, (
	    ftnlen)255);
    *n = ipar[1];
    *m = ipar[2];
    *p = ipar[3];
L2001:
    return 0;
/* *** Last line of BB01AD *** */
} /* bb01ad_ */

#undef notes


