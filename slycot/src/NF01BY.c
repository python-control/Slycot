/* NF01BY.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b15 = -2.;
static doublereal c_b22 = 1.;
static doublereal c_b24 = 0.;

/* Subroutine */ int nf01by_(char *cjte, integer *nsmp, integer *nz, integer *
	l, integer *ipar, integer *lipar, doublereal *wb, integer *lwb, 
	doublereal *z__, integer *ldz, doublereal *e, doublereal *j, integer *
	ldj, doublereal *jte, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen cjte_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, k, m, ib, di, nn, ws, bp1, nwb;
    static doublereal tmp;
    static logical wjte;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


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

/*     To compute the Jacobian of the error function for a neural network */
/*     of the structure */

/*             - tanh(w1*z+b1) - */
/*           /      :            \ */
/*         z ---    :          --- sum(ws(i)*...)+ b(n+1)  --- y, */
/*           \      :            / */
/*             - tanh(wn*z+bn) - */

/*     for the single-output case. The Jacobian has the form */

/*                d e(1)  / d WB(1)   ...    d e(1)  / d WB(NWB) */
/*         J =            :                          :           , */
/*              d e(NSMP) / d WB(1)   ...  d e(NSMP) / d WB(NWB) */

/*     where e(z) is the error function, WB is the set of weights and */
/*     biases of the network (for the considered output), and NWB is */
/*     the number of elements of this set, NWB = IPAR(1)*(NZ+2)+1 */
/*     (see below). */

/*     In the multi-output case, this routine should be called for each */
/*     output. */

/*     NOTE: this routine must have the same arguments as SLICOT Library */
/*     routine NF01BD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CJTE    CHARACTER*1 */
/*             Specifies whether the matrix-vector product J'*e should be */
/*             computed or not, as follows: */
/*             = 'C' :  compute J'*e; */
/*             = 'N' :  do not compute J'*e. */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     NZ      (input) INTEGER */
/*             The length of each input sample.  NZ >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample. */
/*             Currently, L must be 1. */

/*     IPAR    (input/output) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             On entry, the first element of this array must contain */
/*             a value related to the number of neurons, n; specifically, */
/*             n = abs(IPAR(1)), since setting IPAR(1) < 0 has a special */
/*             meaning (see below). */
/*             On exit, if IPAR(1) < 0 on entry, then no computations are */
/*             performed, except the needed tests on input parameters, */
/*             but the following values are returned: */
/*             IPAR(1) contains the length of the array J, LJ; */
/*             LDJ     contains the leading dimension of array J. */
/*             Otherwise, IPAR(1) and LDJ are unchanged on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the vector IPAR.  LIPAR >= 1. */

/*     WB      (input) DOUBLE PRECISION array, dimension (LWB) */
/*             The leading NWB = IPAR(1)*(NZ+2)+1 part of this array */
/*             must contain the weights and biases of the network, */
/*             WB = ( w(1,1), ..., w(1,NZ), ..., w(n,1), ...,  w(n,NZ), */
/*                    ws(1), ..., ws(n), b(1), ..., b(n+1) ), */
/*             where w(i,j) are the weights of the hidden layer, */
/*             ws(i) are the weights of the linear output layer and */
/*             b(i) are the biases. */

/*     LWB     (input) INTEGER */
/*             The length of array WB.  LWB >= NWB. */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ) */
/*             The leading NSMP-by-NZ part of this array must contain the */
/*             set of input samples, */
/*             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,NSMP). */

/*     E       (input) DOUBLE PRECISION array, dimension (NSMP) */
/*             If CJTE = 'C', this array must contain the error vector e. */
/*             If CJTE = 'N', this array is not referenced. */

/*     J       (output) DOUBLE PRECISION array, dimension (LDJ, NWB) */
/*             The leading NSMP-by-NWB part of this array contains the */
/*             Jacobian of the error function. */

/*     LDJ     INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NSMP). */
/*             Note that LDJ is an input parameter, except for */
/*             IPAR(1) < 0 on entry, when it is an output parameter. */

/*     JTE     (output) DOUBLE PRECISION array, dimension (NWB) */
/*             If CJTE = 'C', this array contains the matrix-vector */
/*             product J'*e. */
/*             If CJTE = 'N', this array is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             This argument is included for combatibility with SLICOT */
/*             Library routine NF01BD. */

/*     LDWORK  INTEGER */
/*             Normally, the length of the array DWORK.  LDWORK >= 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Jacobian is computed analytically. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Input output description, neural network, nonlinear system, */
/*     optimization, system response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipar;
    --wb;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --e;
    j_dim1 = *ldj;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --jte;
    --dwork;

    /* Function Body */
    wjte = lsame_(cjte, "C", (ftnlen)1, (ftnlen)1);
    *info = 0;
    nn = ipar[1];
    nwb = nn * (*nz + 2) + 1;
    if (! (wjte || lsame_(cjte, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*nsmp < 0) {
	*info = -2;
    } else if (*nz < 0) {
	*info = -3;
    } else if (*l != 1) {
	*info = -4;
    } else if (*lipar < 1) {
	*info = -6;
    } else if (ipar[1] < 0) {
	if (*info != 0) {
	    i__1 = -(*info);
	    xerbla_("NF01BY", &i__1, (ftnlen)6);
	} else {
	    ipar[1] = *nsmp * (abs(nn) * (*nz + 2) + 1);
	    *ldj = *nsmp;
	}
	return 0;
    } else if (*lwb < nwb) {
	*info = -8;
    } else if (*ldz < max(1,*nsmp)) {
	*info = -10;
    } else if (*ldj < max(1,*nsmp)) {
	*info = -13;
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*nsmp,*nz) == 0) {
	return 0;
    }

/*     Set parameters to avoid overflows and increase accuracy for */
/*     extreme values. */

    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = log(smlnum);
    bignum = log(bignum);

    ws = *nz * nn + 1;
    ib = ws + nn;
    bp1 = ib + nn;

    j[bp1 * j_dim1 + 1] = 1.;
    dcopy_(nsmp, &j[bp1 * j_dim1 + 1], &c__0, &j[bp1 * j_dim1 + 1], &c__1);

    i__1 = nn - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
	dcopy_(nsmp, &wb[ib + i__], &c__0, &j[(ws + i__) * j_dim1 + 1], &c__1)
		;
/* L10: */
    }

    dgemm_("NoTranspose", "NoTranspose", nsmp, &nn, nz, &c_b15, &z__[z_offset]
	    , ldz, &wb[1], nz, &c_b15, &j[ws * j_dim1 + 1], ldj, (ftnlen)11, (
	    ftnlen)11);
    di = 1;

    i__1 = nn - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {

	i__2 = *nsmp;
	for (k = 1; k <= i__2; ++k) {
	    tmp = j[k + (ws + i__) * j_dim1];
	    if (abs(tmp) >= bignum) {
		if (tmp > 0.) {
		    j[k + (ws + i__) * j_dim1] = -1.;
		} else {
		    j[k + (ws + i__) * j_dim1] = 1.;
		}
	    } else if (abs(tmp) <= smlnum) {
		j[k + (ws + i__) * j_dim1] = 0.;
	    } else {
		j[k + (ws + i__) * j_dim1] = 2. / (exp(tmp) + 1.) - 1.;
	    }
/* Computing 2nd power */
	    d__1 = j[k + (ws + i__) * j_dim1];
	    j[k + (ib + i__) * j_dim1] = wb[ws + i__] * (1. - d__1 * d__1);
/* L20: */
	}

	i__2 = *nz - 1;
	for (k = 0; k <= i__2; ++k) {

	    i__3 = *nsmp;
	    for (m = 1; m <= i__3; ++m) {
		j[m + (di + k) * j_dim1] = j[m + (ib + i__) * j_dim1] * z__[m 
			+ (k + 1) * z_dim1];
/* L30: */
	    }

/* L40: */
	}

	di += *nz;
/* L50: */
    }

    if (wjte) {

/*        Compute J'e. */

	dgemv_("Transpose", nsmp, &nwb, &c_b22, &j[j_offset], ldj, &e[1], &
		c__1, &c_b24, &jte[1], &c__1, (ftnlen)9);
    }

    return 0;

/* *** Last line of NF01BY *** */
} /* nf01by_ */

