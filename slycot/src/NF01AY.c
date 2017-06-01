/* NF01AY.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = -2.;
static doublereal c_b11 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b17 = 1.;

/* Subroutine */ int nf01ay_(integer *nsmp, integer *nz, integer *l, integer *
	ipar, integer *lipar, doublereal *wb, integer *lwb, doublereal *z__, 
	integer *ldz, doublereal *y, integer *ldy, doublereal *dwork, integer 
	*ldwork, integer *info)
{
    /* System generated locals */
    integer y_dim1, y_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal df;
    static integer ib, mf, lj, lk, nn, nv, ws;
    static doublereal tmp;
    static integer ldwb;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical last;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlabad_(doublereal *, doublereal *);
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

/*     To calculate the output of a set of neural networks with the */
/*     structure */

/*             - tanh(w1'*z+b1) - */
/*           /      :             \ */
/*         z ---    :           --- sum(ws(i)*...)+ b(n+1)  --- y, */
/*           \      :             / */
/*             - tanh(wn'*z+bn) - */

/*     given the input z and the parameter vectors wi, ws, and b, */
/*     where z, w1, ..., wn are vectors of length NZ, ws is a vector */
/*     of length n, b(1), ..., b(n+1) are scalars, and n is called the */
/*     number of neurons in the hidden layer, or just number of neurons. */
/*     Such a network is used for each L output variables. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     NZ      (input) INTEGER */
/*             The length of each input sample.  NZ >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample.  L >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             IPAR(1) must contain the number of neurons, n, per output */
/*             variable, denoted NN in the sequel.  NN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the vector IPAR.  LIPAR >= 1. */

/*     WB      (input) DOUBLE PRECISION array, dimension (LWB) */
/*             The leading (NN*(NZ+2)+1)*L part of this array must */
/*             contain the weights and biases of the network. This vector */
/*             is partitioned into L vectors of length NN*(NZ+2)+1, */
/*             WB = [ wb(1), ..., wb(L) ]. Each wb(k), k = 1, ..., L, */
/*             corresponds to one output variable, and has the structure */
/*             wb(k) = [ w1(1), ..., w1(NZ), ..., wn(1), ..., wn(NZ), */
/*                       ws(1), ..., ws(n), b(1), ..., b(n+1) ], */
/*             where wi(j) are the weights of the hidden layer, */
/*             ws(i) are the weights of the linear output layer, and */
/*             b(i) are the biases, as in the scheme above. */

/*     LWB     (input) INTEGER */
/*             The length of the array WB. */
/*             LWB >= ( NN*(NZ + 2) + 1 )*L. */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ) */
/*             The leading NSMP-by-NZ part of this array must contain the */
/*             set of input samples, */
/*             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ). */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= MAX(1,NSMP). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array contains the set */
/*             of output samples, */
/*             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ). */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 2*NN. */
/*             For better performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     BLAS routines are used to compute the matrix-vector products. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Input output description, neural network, nonlinear system, */
/*     simulation, system response. */

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
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    nn = ipar[1];
    ldwb = nn * (*nz + 2) + 1;
    if (*nsmp < 0) {
	*info = -1;
    } else if (*nz < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (nn < 0) {
	*info = -4;
    } else if (*lipar < 1) {
	*info = -5;
    } else if (*lwb < ldwb * *l) {
	*info = -7;
    } else if (*ldz < max(1,*nsmp)) {
	*info = -9;
    } else if (*ldy < max(1,*nsmp)) {
	*info = -11;
    } else if (*ldwork < nn << 1) {
	*info = -13;
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01AY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*nsmp,*l) == 0) {
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
    ib = ws + nn - 1;
    lk = 0;
    if (min(*nz,nn) == 0) {
	nv = 2;
    } else {
	nv = (*ldwork - nn) / nn;
    }

    if (nv > 2) {
	mf = *nsmp / nv * nv;
	last = *nsmp % nv != 0;

/*        Some BLAS 3 calculations can be used. */

	i__1 = *l - 1;
	for (k = 0; k <= i__1; ++k) {
	    tmp = wb[ib + nn + 1 + lk];

	    i__2 = nn;
	    for (j = 1; j <= i__2; ++j) {
		dwork[j] = wb[ib + j + lk] * 2.;
/* L10: */
	    }

	    i__2 = mf;
	    i__3 = nv;
	    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {

/*              Compute -2*[w1 w2 ... wn]'*Z', where */
/*              Z = [z(i)';...; z(i+NV-1)']. */

		dgemm_("Transpose", "Transpose", &nn, &nv, nz, &c_b10, &wb[lk 
			+ 1], nz, &z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn 
			+ 1], &nn, (ftnlen)9, (ftnlen)9);
		lj = nn;

		i__4 = nv;
		for (m = 1; m <= i__4; ++m) {
		    i__5 = nn;
		    for (j = 1; j <= i__5; ++j) {

/*                    Compute tanh(wj'*z(i) + bj), j = 1:n. */

			++lj;
			df = dwork[lj] - dwork[j];
			if (abs(df) >= bignum) {
			    if (df > 0.) {
				dwork[lj] = -1.;
			    } else {
				dwork[lj] = 1.;
			    }
			} else if (abs(df) <= smlnum) {
			    dwork[lj] = 0.;
			} else {
			    dwork[lj] = 2. / (exp(df) + 1.) - 1.;
			}
/* L20: */
		    }

/* L30: */
		}

		y[i__ + (k + 1) * y_dim1] = tmp;
		i__4 = nv - 1;
		dcopy_(&i__4, &y[i__ + (k + 1) * y_dim1], &c__0, &y[i__ + 1 + 
			(k + 1) * y_dim1], &c__1);
		dgemv_("Transpose", &nn, &nv, &c_b17, &dwork[nn + 1], &nn, &
			wb[ws + lk], &c__1, &c_b17, &y[i__ + (k + 1) * y_dim1]
			, &c__1, (ftnlen)9);
/* L40: */
	    }

	    if (last) {

/*              Process the last samples. */

		nv = *nsmp - mf;
		i__ = mf + 1;

/*              Compute -2*[w1 w2 ... wn]'*Z', where */
/*              Z = [z(i)';...; z(NSMP)']. */

		dgemm_("Transpose", "Transpose", &nn, &nv, nz, &c_b10, &wb[lk 
			+ 1], nz, &z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn 
			+ 1], &nn, (ftnlen)9, (ftnlen)9);
		lj = nn;

		i__3 = nv;
		for (m = 1; m <= i__3; ++m) {
		    i__2 = nn;
		    for (j = 1; j <= i__2; ++j) {

/*                    Compute tanh(wj'*z(i) + bj), j = 1:n. */

			++lj;
			df = dwork[lj] - dwork[j];
			if (abs(df) >= bignum) {
			    if (df > 0.) {
				dwork[lj] = -1.;
			    } else {
				dwork[lj] = 1.;
			    }
			} else if (abs(df) <= smlnum) {
			    dwork[lj] = 0.;
			} else {
			    dwork[lj] = 2. / (exp(df) + 1.) - 1.;
			}
/* L50: */
		    }

/* L60: */
		}

		y[i__ + (k + 1) * y_dim1] = tmp;
		if (nv > 1) {
		    i__3 = nv - 1;
		    dcopy_(&i__3, &y[i__ + (k + 1) * y_dim1], &c__0, &y[i__ + 
			    1 + (k + 1) * y_dim1], &c__1);
		}
		dgemv_("Transpose", &nn, &nv, &c_b17, &dwork[nn + 1], &nn, &
			wb[ws + lk], &c__1, &c_b17, &y[i__ + (k + 1) * y_dim1]
			, &c__1, (ftnlen)9);
	    }

	    lk += ldwb;
/* L70: */
	}

    } else {

/*        BLAS 2 calculations only can be used. */

	i__1 = *l - 1;
	for (k = 0; k <= i__1; ++k) {
	    tmp = wb[ib + nn + 1 + lk];

	    i__3 = nn;
	    for (j = 1; j <= i__3; ++j) {
		dwork[j] = wb[ib + j + lk] * 2.;
/* L80: */
	    }

	    i__3 = *nsmp;
	    for (i__ = 1; i__ <= i__3; ++i__) {

/*              Compute -2*[w1 w2 ... wn]'*z(i). */

		if (*nz == 0) {
		    dwork[nn + 1] = 0.;
		    dcopy_(&nn, &dwork[nn + 1], &c__0, &dwork[nn + 1], &c__1);
		} else {
		    dgemv_("Transpose", nz, &nn, &c_b10, &wb[lk + 1], nz, &
			    z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn + 1], &
			    c__1, (ftnlen)9);
		}

		i__2 = nn << 1;
		for (j = nn + 1; j <= i__2; ++j) {

/*                 Compute tanh(wj'*z(i) + bj), j = 1:n. */

		    df = dwork[j] - dwork[j - nn];
		    if (abs(df) >= bignum) {
			if (df > 0.) {
			    dwork[j] = -1.;
			} else {
			    dwork[j] = 1.;
			}
		    } else if (abs(df) <= smlnum) {
			dwork[j] = 0.;
		    } else {
			dwork[j] = 2. / (exp(df) + 1.) - 1.;
		    }
/* L90: */
		}

		y[i__ + (k + 1) * y_dim1] = ddot_(&nn, &wb[ws + lk], &c__1, &
			dwork[nn + 1], &c__1) + tmp;
/* L100: */
	    }

	    lk += ldwb;
/* L110: */
	}

    }
    return 0;

/* *** Last line of NF01AY *** */
} /* nf01ay_ */

