/* BD01AD.f -- translated by f2c (version 20100827).
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

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static integer c__5 = 5;

/* Subroutine */ int bd01ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *vec, integer *n, integer *m, integer *p, 
	doublereal *e, integer *lde, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, 
	integer *ldd, char *note, doublereal *dwork, integer *ldwork, integer 
	*info, ftnlen def_len, ftnlen note_len)
{
    /* Initialized data */

    static logical vecdef[8] = { TRUE_,TRUE_,TRUE_,FALSE_,TRUE_,TRUE_,TRUE_,
	    FALSE_ };

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, l;
    static doublereal b1, b2, c1, c2, temp;
    static char dataf[12];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal ttemp, appind;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static integer status;

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, dataf, 0, "(A,I2.2,A)", 11, 1 };
    static cilist io___6 = { 1, 1, 1, 0, 0 };
    static cilist io___8 = { 1, 1, 1, 0, 0 };
    static cilist io___9 = { 1, 1, 1, 0, 0 };
    static cilist io___10 = { 1, 1, 1, 0, 0 };
    static cilist io___11 = { 1, 1, 1, 0, 0 };
    static cilist io___12 = { 1, 1, 1, 0, 0 };
    static icilist io___13 = { 0, dataf, 0, "(A,I1,A)", 12, 1 };
    static cilist io___14 = { 1, 1, 1, 0, 0 };
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___16 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___18 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };
    static cilist io___21 = { 1, 1, 1, 0, 0 };



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

/*     To generate benchmark examples for time-invariant, */
/*     continuous-time dynamical systems */

/*         . */
/*       E x(t) = A x(t) + B u(t) */

/*         y(t) = C x(t) + D u(t) */

/*     E, A are real N-by-N matrices, B is N-by-M, C is P-by-N, and */
/*     D is P-by-M. In many examples, E is the identity matrix and D is */
/*     the zero matrix. */

/*     This routine is an implementation of the benchmark library */
/*     CTDSX (Version 1.0) described in [1]. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER*1 */
/*             Specifies the kind of values used as parameters when */
/*             generating parameter-dependent and scalable examples */
/*             (i.e., examples with NR(1) = 2, 3, or 4): */
/*             = 'D':  Default values defined in [1] are used; */
/*             = 'N':  Values set in DPAR and IPAR are used. */
/*             This parameter is not referenced if NR(1) = 1. */
/*             Note that the scaling parameter of examples with */
/*             NR(1) = 3 or 4 is considered as a regular parameter in */
/*             this context. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension (2) */
/*             Specifies the index of the desired example according */
/*             to [1]. */
/*             NR(1) defines the group: */
/*                   1 : parameter-free problems of fixed size */
/*                   2 : parameter-dependent problems of fixed size */
/*                   3 : parameter-free problems of scalable size */
/*                   4 : parameter-dependent problems of scalable size */
/*             NR(2) defines the number of the benchmark example */
/*             within a certain group according to [1]. */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (7) */
/*             On entry, if DEF = 'N' and the desired example depends on */
/*             real parameters, then the array DPAR must contain the */
/*             values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 2.1 and 2.2, DPAR(1) defines the parameter */
/*             'epsilon'. */
/*             For Example 2.4, DPAR(1), ..., DPAR(7) define 'b', 'mu', */
/*             'r', 'r_c', 'k_l', 'sigma', 'a', respectively. */
/*             For Example 2.7, DPAR(1) and DPAR(2) define 'mu' and 'nu', */
/*             respectively. */
/*             For Example 4.1, DPAR(1), ..., DPAR(7) define 'a', 'b', */
/*             'c', 'beta_1', 'beta_2', 'gamma_1', 'gamma_2', */
/*             respectively. */
/*             For Example 4.2, DPAR(1), ..., DPAR(3) define 'mu', */
/*             'delta', 'kappa', respectively. */
/*             On exit, if DEF = 'D' and the desired example depends on */
/*             real parameters, then the array DPAR is overwritten by the */
/*             default values given in [1]. */

/*     IPAR    (input/output) INTEGER array, dimension (1) */
/*             On entry, if DEF = 'N' and the desired example depends on */
/*             integer parameters, then the array IPAR must contain the */
/*             values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 2.3, 2.5, and 2.6, IPAR(1) defines the */
/*             parameter 's'. */
/*             For Example 3.1, IPAR(1) defines 'q'. */
/*             For Examples 3.2 and 3.3, IPAR(1) defines 'n'. */
/*             For Example 3.4, IPAR(1) defines 'l'. */
/*             For Example 4.1, IPAR(1) defines 'n'. */
/*             For Example 4.2, IPAR(1) defines 'l'. */
/*             On exit, if DEF = 'D' and the desired example depends on */
/*             integer parameters, then the array IPAR is overwritten by */
/*             the default values given in [1]. */

/*     VEC     (output) LOGICAL array, dimension (8) */
/*             Flag vector which displays the availabilty of the output */
/*             data: */
/*             VEC(1), ..., VEC(3) refer to N, M, and P, respectively, */
/*             and are always .TRUE.. */
/*             VEC(4) is .TRUE. iff E is NOT the identity matrix. */
/*             VEC(5), ..., VEC(7) refer to A, B, and C, respectively, */
/*             and are always .TRUE.. */
/*             VEC(8) is .TRUE. iff D is NOT the zero matrix. */

/*     N       (output) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices E and A. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrices B and D. */

/*     P       (output) INTEGER */
/*             The number of rows in the matrices C and D. */

/*     E       (output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix E. */
/*             NOTE that this array is overwritten (by the identity */
/*             matrix), if VEC(4) = .FALSE.. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= N. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array contains the */
/*             matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array contains the */
/*             matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= P. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array contains the */
/*             matrix D. */
/*             NOTE that this array is overwritten (by the zero */
/*             matrix), if VEC(8) = .FALSE.. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= P. */

/*     NOTE    (output) CHARACTER*70 */
/*             String containing short information about the chosen */
/*             example. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             For Example 3.4, LDWORK >= 4*IPAR(1) is required. */
/*             For the other examples, no workspace is needed, i.e., */
/*             LDWORK >= 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; in particular, INFO = -3 or -4 indicates */
/*                   that at least one of the parameters in DPAR or */
/*                   IPAR, respectively, has an illegal value; */
/*             = 1:  data file can not be opened or has wrong format. */


/*     REFERENCES */

/*     [1]  Kressner, D., Mehrmann, V. and Penzl, T. */
/*          CTDSX - a Collection of Benchmark Examples for State-Space */
/*          Realizations of Continuous-Time Dynamical Systems. */
/*          SLICOT Working Note 1998-9. 1998. */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz) */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please contact Volker Mehrmann */
/*     (Email: volker.mehrmann@mathematik.tu-chemnitz.de). */

/*     REVISIONS */

/*     June 1999, V. Sima. */

/*     KEYWORDS */

/*     continuous-time dynamical systems */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     . LAPACK . */
/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. Intrinsic Functions .. */
/*     .. Data Statements .. */
/*     . default values for availabities . */
    /* Parameter adjustments */
    --nr;
    --dpar;
    --ipar;
    --vec;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --dwork;

    /* Function Body */

/*     .. Executable Statements .. */

    *info = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	vec[i__] = vecdef[i__ - 1];
/* L10: */
    }

    if (nr[1] == 1) {

	if (nr[2] == 1) {
	    s_copy(note, "Laub 1979, Ex.1", (ftnlen)70, (ftnlen)15);
	    *n = 2;
	    *m = 1;
	    *p = 2;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[(a_dim1 << 1) + 1] = 1.;
	    b[b_dim1 + 1] = 0.;
	    b[b_dim1 + 2] = 1.;
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 2) {
	    s_copy(note, "Laub 1979, Ex.2: uncontrollable-unobservable data", 
		    (ftnlen)70, (ftnlen)49);
	    *n = 2;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    a[a_dim1 + 1] = 4.;
	    a[a_dim1 + 2] = -4.5;
	    a[(a_dim1 << 1) + 1] = 3.;
	    a[(a_dim1 << 1) + 2] = -3.5;
	    b[b_dim1 + 1] = 1.;
	    b[b_dim1 + 2] = -1.;
	    c__[c_dim1 + 1] = 3.;
	    c__[(c_dim1 << 1) + 1] = 2.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 3) {
	    s_copy(note, "Beale/Shafai 1989: model of L-1011 aircraft", (
		    ftnlen)70, (ftnlen)43);
	    *n = 4;
	    *m = 2;
	    *p = 4;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 4) {
	    s_copy(note, "Bhattacharyya et al. 1983: binary distillation col"
		    "umn", (ftnlen)70, (ftnlen)53);
	    *n = 8;
	    *m = 2;
	    *p = 8;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 5) {
	    s_copy(note, "Patnaik et al. 1980: tubular ammonia reactor", (
		    ftnlen)70, (ftnlen)44);
	    *n = 9;
	    *m = 3;
	    *p = 9;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 6) {
	    s_copy(note, "Davison/Gesing 1978: J-100 jet engine", (ftnlen)70, 
		    (ftnlen)37);
	    *n = 30;
	    *m = 3;
	    *p = 5;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 7) {
	    s_copy(note, "Davison 1967: binary distillation column", (ftnlen)
		    70, (ftnlen)40);
	    *n = 11;
	    *m = 3;
	    *p = 3;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 + 2] = 1.;
	    c__[c_dim1 * 10 + 1] = 1.;
	    c__[c_dim1 * 11 + 3] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);
	} else if (nr[2] == 8) {
	    s_copy(note, "Chien/Ergin/Ling/Lee 1958: drum boiler", (ftnlen)70,
		     (ftnlen)38);
	    *n = 9;
	    *m = 3;
	    *p = 2;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 * 6 + 1] = 1.;
	    c__[c_dim1 * 9 + 2] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 9) {
	    s_copy(note, "Ly, Gangsaas 1981: B-767 airplane", (ftnlen)70, (
		    ftnlen)33);
	    *n = 55;
	    *m = 2;
	    *p = 2;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 10) {
	    s_copy(note, "control surface servo for an underwater vehicle", (
		    ftnlen)70, (ftnlen)47);
	    *n = 8;
	    *m = 2;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 * 7 + 1] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);
	} else {
	    *info = -2;
	}

	if (nr[2] >= 3 && nr[2] <= 10) {
/*         .. loading data files */
	    s_wsfi(&io___4);
	    do_fio(&c__1, "BD011", (ftnlen)5);
	    do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".dat", (ftnlen)4);
	    e_wsfi();
	    o__1.oerr = 1;
	    o__1.ounit = 1;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = dataf;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    status = f_open(&o__1);
	    if (status != 0) {
		*info = 1;
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___6);
		    if (status != 0) {
			goto L100001;
		    }
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100001;
			}
		    }
		    status = e_rsle();
L100001:
		    if (status != 0) {
			*info = 1;
		    }
/* L110: */
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___8);
		    if (status != 0) {
			goto L100002;
		    }
		    i__2 = *m;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100002;
			}
		    }
		    status = e_rsle();
L100002:
		    if (status != 0) {
			*info = 1;
		    }
/* L120: */
		}
		if (nr[2] == 6 || nr[2] == 9) {
		    i__1 = *p;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			status = s_rsle(&io___9);
			if (status != 0) {
			    goto L100003;
			}
			i__2 = *n;
			for (j = 1; j <= i__2; ++j) {
			    status = do_lio(&c__5, &c__1, (char *)&c__[i__ + 
				    j * c_dim1], (ftnlen)sizeof(doublereal));
			    if (status != 0) {
				goto L100003;
			    }
			}
			status = e_rsle();
L100003:
			if (status != 0) {
			    *info = 1;
			}
/* L130: */
		    }
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

    } else if (nr[1] == 2) {
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
	    *info = -1;
	    return 0;
	}

	if (nr[2] == 1) {
	    s_copy(note, "Chow/Kokotovic 1976: magnetic tape control system", 
		    (ftnlen)70, (ftnlen)49);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e-6;
	    }
	    if (dpar[1] == 0.) {
		*info = -3;
	    }
	    *n = 4;
	    *m = 1;
	    *p = 2;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[(a_dim1 << 1) + 1] = .4;
	    a[a_dim1 * 3 + 2] = .345;
	    a[(a_dim1 << 1) + 3] = -.524 / dpar[1];
	    a[a_dim1 * 3 + 3] = -.465 / dpar[1];
	    a[(a_dim1 << 2) + 3] = .262 / dpar[1];
	    a[(a_dim1 << 2) + 4] = -1. / dpar[1];
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[b_dim1 + 4] = 1. / dpar[1];
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 + 1] = 1.;
	    c__[c_dim1 * 3 + 2] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 2) {
	    s_copy(note, "Arnold/Laub 1984", (ftnlen)70, (ftnlen)16);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1e-6;
	    }
	    *n = 4;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &dpar[1], &a[a_offset], lda, (ftnlen)1);
	    a[a_dim1 + 1] = -dpar[1];
	    a[a_dim1 + 2] = -1.;
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = -dpar[1];
	    a[a_dim1 * 3 + 4] = -1.;
	    a[(a_dim1 << 2) + 3] = 1.;
	    dlaset_("A", n, m, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
	    dlaset_("A", p, n, &c_b6, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    d__[d_dim1 + 1] = 0.;

	} else if (nr[2] == 3) {
	    s_copy(note, "Vertical acceleration of a rigid guided missile", (
		    ftnlen)70, (ftnlen)47);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 1;
	    }
	    if (ipar[1] < 1 || ipar[1] > 10) {
		*info = -4;
	    }
	    *n = 3;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[a_dim1 + 2] = 1.;
	    a[a_dim1 * 3 + 3] = -190.;
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[b_dim1 + 3] = 190.;
	    d__[d_dim1 + 1] = 0.;
	    o__1.oerr = 1;
	    o__1.ounit = 1;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = "BD01203.dat";
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    status = f_open(&o__1);
	    if (status != 0) {
		*info = 1;
	    } else {
		i__1 = ipar[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___10);
		    if (status != 0) {
			goto L100004;
		    }
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100004;
			}
		    }
		    status = e_rsle();
L100004:
		    if (status != 0) {
			*info = 1;
		    }
		    status = s_rsle(&io___11);
		    if (status != 0) {
			goto L100005;
		    }
		    i__2 = *n;
		    for (j = 2; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				2], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100005;
			}
		    }
		    status = e_rsle();
L100005:
		    if (status != 0) {
			*info = 1;
		    }
		    status = s_rsle(&io___12);
		    if (status != 0) {
			goto L100006;
		    }
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&c__[j * c_dim1 
				+ 1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100006;
			}
		    }
		    status = e_rsle();
L100006:
		    if (status != 0) {
			*info = 1;
		    }
/* L210: */
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);

	} else if (nr[2] == 4) {
	    s_copy(note, "Senning 1980: hydraulic positioning system", (
		    ftnlen)70, (ftnlen)42);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 1.4e4;
		dpar[2] = .1287;
		dpar[3] = .15;
		dpar[4] = .01;
		dpar[5] = .002;
		dpar[6] = .24;
		dpar[7] = 10.75;
	    }
	    if (dpar[1] <= 9e3 || dpar[1] >= 1.6e4 || (dpar[2] <= .05 || dpar[
		    2] >= .3) || (dpar[3] <= .05 || dpar[3] >= 5.) || (dpar[4]
		     <= 0. || dpar[4] >= .05) || (dpar[5] <= 1.03e-4 || dpar[
		    5] >= .0035) || (dpar[6] <= .001 || dpar[6] >= 15.) || (
		    dpar[7] <= 10.5 || dpar[7] >= 11.1)) {
		*info = -3;
	    }
	    *n = 3;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[(a_dim1 << 1) + 1] = 1.;
	    a[(a_dim1 << 1) + 2] = -(dpar[3] + dpar[4] * 4. / 
		    3.141592653589793) / dpar[2];
	    a[a_dim1 * 3 + 2] = dpar[7] / dpar[2];
	    a[(a_dim1 << 1) + 3] = dpar[7] * -4. * dpar[1] / 874.;
	    a[a_dim1 * 3 + 3] = dpar[1] * -4. * (dpar[6] + dpar[5]) / 874.;
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[b_dim1 + 3] = dpar[1] * -4. / 874.;
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    d__[d_dim1 + 1] = 0.;

	} else if (nr[2] == 5) {
	    s_copy(note, "Kwakernaak/Westdyk 1985: cascade of inverted pendu"
		    "la", (ftnlen)70, (ftnlen)52);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 1;
	    }
	    if (ipar[1] < 1 || ipar[1] > 7) {
		*info = -4;
	    }
	    if (ipar[1] <= 6) {
		*m = ipar[1];
	    } else {
		*m = 10;
	    }
	    *n = *m << 1;
	    *p = *m;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    s_wsfi(&io___13);
	    do_fio(&c__1, "BD01205", (ftnlen)7);
	    do_fio(&c__1, (char *)&ipar[1], (ftnlen)sizeof(integer));
	    do_fio(&c__1, ".dat", (ftnlen)4);
	    e_wsfi();
	    o__1.oerr = 1;
	    o__1.ounit = 1;
	    o__1.ofnmlen = 12;
	    o__1.ofnm = dataf;
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    status = f_open(&o__1);
	    if (status != 0) {
		*info = 1;
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___14);
		    if (status != 0) {
			goto L100007;
		    }
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100007;
			}
		    }
		    status = e_rsle();
L100007:
		    if (status != 0) {
			*info = 1;
		    }
/* L220: */
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___15);
		    if (status != 0) {
			goto L100008;
		    }
		    i__2 = *m;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100008;
			}
		    }
		    status = e_rsle();
L100008:
		    if (status != 0) {
			*info = 1;
		    }
/* L230: */
		}
		i__1 = *p;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___16);
		    if (status != 0) {
			goto L100009;
		    }
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				c_dim1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100009;
			}
		    }
		    status = e_rsle();
L100009:
		    if (status != 0) {
			*info = 1;
		    }
/* L240: */
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 6) {
	    s_copy(note, "Kallstrom/Astrom 1981: regulation of a ship heading"
		    , (ftnlen)70, (ftnlen)51);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 1;
	    }
	    if (ipar[1] < 1 || ipar[1] > 5) {
		*info = -4;
	    }
	    *n = 3;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[(a_dim1 << 1) + 3] = 1.;
	    b[b_dim1 + 3] = 0.;
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 * 3 + 1] = 1.;
	    d__[d_dim1 + 1] = 0.;
	    o__1.oerr = 1;
	    o__1.ounit = 1;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = "BD01206.dat";
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    status = f_open(&o__1);
	    if (status != 0) {
		*info = 1;
	    } else {
		i__1 = ipar[1];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___17);
		    if (status != 0) {
			goto L100010;
		    }
		    for (j = 1; j <= 2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				1], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100010;
			}
		    }
		    status = e_rsle();
L100010:
		    if (status != 0) {
			*info = 1;
		    }
		    status = s_rsle(&io___18);
		    if (status != 0) {
			goto L100011;
		    }
		    for (j = 1; j <= 2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				2], (ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100011;
			}
		    }
		    status = e_rsle();
L100011:
		    if (status != 0) {
			*info = 1;
		    }
		    status = s_rsle(&io___19);
		    if (status != 0) {
			goto L100012;
		    }
		    for (j = 1; j <= 2; ++j) {
			status = do_lio(&c__5, &c__1, (char *)&b[j + b_dim1], 
				(ftnlen)sizeof(doublereal));
			if (status != 0) {
			    goto L100012;
			}
		    }
		    status = e_rsle();
L100012:
		    if (status != 0) {
			*info = 1;
		    }
/* L250: */
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);

	} else if (nr[2] == 7) {
	    s_copy(note, "Ackermann 1989: track-guided bus", (ftnlen)70, (
		    ftnlen)32);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		dpar[1] = 15.;
		dpar[2] = 10.;
	    }
	    if (dpar[1] < 9.95 || dpar[1] > 16.) {
		*info = -3;
	    }
	    if (dpar[1] < 1. || dpar[1] > 20.) {
		*info = -3;
	    }
	    *n = 5;
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    a[a_dim1 + 1] = -668. / (dpar[1] * dpar[2]);
/* Computing 2nd power */
	    d__1 = dpar[2];
	    a[(a_dim1 << 1) + 1] = 180.4 / (dpar[1] * (d__1 * d__1)) - 1.;
	    a[a_dim1 + 2] = 180.4 / (dpar[1] * 10.86);
	    a[(a_dim1 << 1) + 2] = -4417.5452 / (dpar[1] * 10.86 * dpar[2]);
	    a[a_dim1 * 5 + 1] = 198 / (dpar[1] * dpar[2]);
	    a[a_dim1 * 5 + 2] = 726.66 / (dpar[1] * 10.86);
	    a[a_dim1 + 3] = dpar[2];
	    a[(a_dim1 << 2) + 3] = dpar[2];
	    a[(a_dim1 << 1) + 4] = 1.;
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[b_dim1 + 5] = 1.;
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 * 3 + 1] = 1.;
	    c__[(c_dim1 << 2) + 1] = 6.12;
	    d__[d_dim1 + 1] = 0.;

	} else {
	    *info = -2;
	}

    } else if (nr[1] == 3) {
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
	    *info = -1;
	    return 0;
	}

	if (nr[2] == 1) {
	    s_copy(note, "Laub 1979, Ex.4: string of high speed vehicles", (
		    ftnlen)70, (ftnlen)46);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 20;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    *n = (ipar[1] << 1) - 1;
	    *m = ipar[1];
	    *p = ipar[1] - 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (i__ % 2 == 1) {
		    a[i__ + i__ * a_dim1] = -1.;
		    b[i__ + (i__ + 1) / 2 * b_dim1] = 1.;
		} else {
		    a[i__ + (i__ - 1) * a_dim1] = 1.;
		    a[i__ + (i__ + 1) * a_dim1] = -1.;
		    c__[i__ / 2 + i__ * c_dim1] = 1.;
		}
/* L310: */
	    }
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 2) {
	    s_copy(note, "Hodel et al. 1996: heat flow in a thin rod", (
		    ftnlen)70, (ftnlen)42);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 100;
	    }
	    if (ipar[1] < 1) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 1;
	    *p = *n;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    temp = (doublereal) (*n + 1);
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    d__1 = temp * -2.;
	    dlaset_("A", n, n, &c_b5, &d__1, &a[a_offset], lda, (ftnlen)1);
	    a[a_dim1 + 1] = -temp;
	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + (i__ + 1) * a_dim1] = temp;
		a[i__ + 1 + i__ * a_dim1] = temp;
/* L320: */
	    }
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[*n + b_dim1] = temp;
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 3) {
	    s_copy(note, "Laub 1979, Ex.6", (ftnlen)70, (ftnlen)15);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 21;
	    }
	    if (ipar[1] < 1) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[*n + b_dim1] = 1.;
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    c__[c_dim1 + 1] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 4) {
	    s_copy(note, "Lang/Penzl 1994: rotating axle", (ftnlen)70, (
		    ftnlen)30);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 211;
	    }
	    if (ipar[1] < 1 || ipar[1] > 211) {
		*info = -4;
	    }
	    *n = (ipar[1] << 1) - 1;
	    *m = ipar[1];
	    *p = ipar[1];
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*ldwork < *m << 2) {
		*info = -21;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    o__1.oerr = 1;
	    o__1.ounit = 1;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = "BD01304.dat";
	    o__1.orl = 0;
	    o__1.osta = "OLD";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    status = f_open(&o__1);
	    if (status != 0) {
		*info = 1;
	    } else {
		i__1 = *m << 2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    status = s_rsle(&io___21);
		    if (status != 0) {
			goto L100013;
		    }
		    status = do_lio(&c__5, &c__1, (char *)&dwork[i__], (
			    ftnlen)sizeof(doublereal));
		    if (status != 0) {
			goto L100013;
		    }
		    status = e_rsle();
L100013:
		    if (status != 0) {
			*info = 1;
		    }
/* L330: */
		}
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 1;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    if (*info != 0) {
		return 0;
	    }
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
	    e[e_dim1 + 1] = dwork[1];
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		e[i__ + (i__ - 1) * e_dim1] = dwork[(i__ - 2 << 2) + 1];
		e[i__ + i__ * e_dim1] = -dwork[(i__ - 1 << 2) + 1];
/* L340: */
	    }
	    e[*m + *m * e_dim1] = -e[*m + *m * e_dim1];
	    for (i__ = *m - 1; i__ >= 1; --i__) {
		i__1 = *m;
		for (j = i__; j <= i__1; ++j) {
		    if (i__ == 1) {
			e[j + i__ * e_dim1] -= e[j + (i__ + 1) * e_dim1];
		    } else {
			e[j + i__ * e_dim1] = e[j + (i__ + 1) * e_dim1] - e[j 
				+ i__ * e_dim1];
		    }
/* L345: */
		}
/* L350: */
	    }
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		a[i__ - 1 + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3];
		a[i__ + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3] * -2. - 
			dwork[(i__ - 1 << 2) + 2];
		a[i__ + a_dim1] = dwork[(i__ - 1 << 2) + 2] - dwork[(i__ - 2 
			<< 2) + 2];
		a[i__ - 1 + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 4];
		a[i__ + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 4] * -2.;
		if (i__ < *m) {
		    a[i__ + 1 + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3];
		    i__2 = *m;
		    for (j = i__ + 1; j <= i__2; ++j) {
			a[j + i__ * a_dim1] = a[j + i__ * a_dim1] + dwork[(j 
				- 2 << 2) + 2] - dwork[(j - 1 << 2) + 2];
/* L355: */
		    }
		    a[i__ + 1 + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 
			    4];
		}
/* L360: */
	    }
	    a[a_dim1 + 1] = -dwork[2];
	    a[(a_dim1 << 1) + 1] = -dwork[3];
	    a[(*m + 1) * a_dim1 + 1] = -a[(*m + 1) * a_dim1 + 1];
	    i__1 = *m - 1;
	    i__2 = *m - 1;
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[*m + 1 + (a_dim1 << 1)
		    ], lda, (ftnlen)1);
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		b[i__ + i__ * b_dim1] = -1.;
		b[i__ + (i__ - 1) * b_dim1] = 1.;
		c__[i__ + i__ * c_dim1] = dwork[(i__ - 2 << 2) + 3];
		c__[i__ + (*m + i__ - 1) * c_dim1] = dwork[(i__ - 1) * 4];
/* L370: */
	    }
	    b[b_dim1 + 1] = 1.;
	    c__[c_dim1 + 1] = 1.;
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else {
	    *info = -2;
	}

    } else if (nr[1] == 4) {
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
	    *info = -1;
	    return 0;
	}

	if (nr[2] == 1) {
	    s_copy(note, "Rosen/Wang 1995: control of 1-dim. heat flow", (
		    ftnlen)70, (ftnlen)44);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 100;
		dpar[1] = .01;
		dpar[2] = 1.;
		dpar[3] = 1.;
		dpar[4] = .2;
		dpar[5] = .3;
		dpar[6] = .2;
		dpar[7] = .3;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    *n = ipar[1];
	    *m = 1;
	    *p = 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[4] = TRUE_;
	    appind = (doublereal) (*n + 1);
	    ttemp = -dpar[1] * appind;
	    temp = 1 / (appind * 6.);
	    d__1 = temp * 4.;
	    dlaset_("A", n, n, &c_b5, &d__1, &e[e_offset], lde, (ftnlen)1);
	    d__1 = ttemp * 2.;
	    dlaset_("A", n, n, &c_b5, &d__1, &a[a_offset], lda, (ftnlen)1);
	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__ + 1 + i__ * a_dim1] = -ttemp;
		a[i__ + (i__ + 1) * a_dim1] = -ttemp;
		e[i__ + 1 + i__ * e_dim1] = temp;
		e[i__ + (i__ + 1) * e_dim1] = temp;
/* L410: */
	    }
	    i__1 = *n;
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
/* L420: */
	    }
	    dscal_(n, &dpar[2], &b[b_dim1 + 1], &c__1);
	    dscal_(n, &dpar[3], &c__[c_dim1 + 1], ldc);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else if (nr[2] == 2) {
	    s_copy(note, "Hench et al. 1995: coupled springs, dashpots, mass"
		    "es", (ftnlen)70, (ftnlen)52);
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
		ipar[1] = 30;
		dpar[1] = 4.;
		dpar[2] = 4.;
		dpar[3] = 1.;
	    }
	    if (ipar[1] < 2) {
		*info = -4;
	    }
	    l = ipar[1];
	    *n = l << 1;
	    *m = 2;
	    *p = l << 1;
	    if (*lde < *n) {
		*info = -10;
	    }
	    if (*lda < *n) {
		*info = -12;
	    }
	    if (*ldb < *n) {
		*info = -14;
	    }
	    if (*ldc < *p) {
		*info = -16;
	    }
	    if (*ldd < *p) {
		*info = -18;
	    }
	    if (*info != 0) {
		return 0;
	    }

	    vec[4] = TRUE_;
	    dlaset_("A", n, n, &c_b5, &dpar[1], &e[e_offset], lde, (ftnlen)1);
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
	    temp = dpar[3] * -2.;
	    i__1 = l;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		e[i__ + i__ * e_dim1] = 1.;
		a[i__ + (i__ + l) * a_dim1] = 1.;
		a[i__ + l + (i__ + l) * a_dim1] = -dpar[2];
		if (i__ < l) {
		    a[i__ + l + (i__ + 1) * a_dim1] = dpar[3];
		    a[i__ + l + 1 + i__ * a_dim1] = dpar[3];
		    if (i__ > 1) {
			a[i__ + l + i__ * a_dim1] = temp;
		    }
		}
/* L430: */
	    }
	    a[l + 1 + a_dim1] = -dpar[3];
	    a[*n + l * a_dim1] = -dpar[3];
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
	    b[l + 1 + b_dim1] = 1.;
	    b[*n + (b_dim1 << 1)] = -1.;
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

	} else {
	    *info = -2;
	}
    } else {
	*info = -2;
    }

    return 0;
/* *** Last Line of BD01AD *** */
} /* bd01ad_ */

