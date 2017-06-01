/* SB10YD.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static integer c__4096 = 4096;
static integer c__2048 = 2048;
static doublereal c_b21 = 2.44140625e-4;
static doublereal c_b48 = 0.;
static doublereal c_b49 = 1.;

/* Subroutine */ int sb10yd_(integer *discfl, integer *flag__, integer *
	lendat, doublereal *rfrdat, doublereal *ifrdat, doublereal *omega, 
	integer *n, doublereal *a, integer *lda, doublereal *b, doublereal *
	c__, doublereal *d__, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, doublecomplex *zwork, integer *lzwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), acos(doublereal), log(
	    doublereal), exp(doublereal), cos(doublereal), sin(doublereal), 
	    d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, k, p, n1, n2;
    static doublereal p1, p2;
    static integer ii;
    static doublereal pi;
    static integer mn;
    static doublereal pw;
    static integer ip1, ip2, lw1, lw2, lw3, lw4;
    static doublereal rat;
    static integer iws, iwa0, iwab, rank;
    static doublereal tolb;
    static integer iwbp;
    static doublecomplex xhat[1024];
    static integer iwbx;
    static doublereal toll;
    static integer iwxi, iwxr, info2;
    extern /* Subroutine */ int ab04md_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dg01md_(char *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *)
	    ;
    static integer iwmag, iwdme;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb10zp_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer iwvar, istop;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer iwbmat;
    extern /* Subroutine */ int dgelsy_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    static integer clwmax, dlwmax, iwymag, iwdomo, istart;


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

/*     To fit a supplied frequency response data with a stable, minimum */
/*     phase SISO (single-input single-output) system represented by its */
/*     matrices A, B, C, D. It handles both discrete- and continuous-time */
/*     cases. */

/*     ARGUMENTS */

/*     Input/Output parameters */

/*     DISCFL  (input) INTEGER */
/*             Indicates the type of the system, as follows: */
/*             = 0: continuous-time system; */
/*             = 1: discrete-time system. */

/*     FLAG    (input) INTEGER */
/*             If FLAG = 0, then the system zeros and poles are not */
/*             constrained. */
/*             If FLAG = 1, then the system zeros and poles will have */
/*             negative real parts in the continuous-time case, or moduli */
/*             less than 1 in the discrete-time case. Consequently, FLAG */
/*             must be equal to 1 in mu-synthesis routines. */

/*     LENDAT  (input) INTEGER */
/*             The length of the vectors RFRDAT, IFRDAT and OMEGA. */
/*             LENDAT >= 2. */

/*     RFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The real part of the frequency data to be fitted. */

/*     IFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The imaginary part of the frequency data to be fitted. */

/*     OMEGA   (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The frequencies corresponding to RFRDAT and IFRDAT. */
/*             These values must be nonnegative and monotonically */
/*             increasing. Additionally, for discrete-time systems */
/*             they must be between 0 and PI. */

/*     N       (input/output) INTEGER */
/*             On entry, the desired order of the system to be fitted. */
/*             N <= LENDAT-1. */
/*             On exit, the order of the obtained system. The value of N */
/*             could only be modified if N > 0 and FLAG = 1. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. If FLAG = 1, then A is in an upper Hessenberg */
/*             form, and corresponds to a minimal realization. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (N) */
/*             The computed vector B. */

/*     C       (output) DOUBLE PRECISION array, dimension (N) */
/*             The computed vector C. If FLAG = 1, the first N-1 elements */
/*             are zero (for the exit value of N). */

/*     D       (output) DOUBLE PRECISION array, dimension (1) */
/*             The computed scalar D. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for determining the effective */
/*             rank of matrices. If the user sets TOL > 0, then the given */
/*             value of TOL is used as a lower bound for the reciprocal */
/*             condition number;  a (sub)matrix whose estimated condition */
/*             number is less than 1/TOL is considered to be of full */
/*             rank.  If the user sets TOL <= 0, then an implicitly */
/*             computed, default tolerance, defined by TOLDEF = SIZE*EPS, */
/*             is used instead, where SIZE is the product of the matrix */
/*             dimensions, and EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2,2*N+1) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) contains the optimal value of */
/*             LZWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max( 2, LW1, LW2, LW3, LW4 ), where */
/*             LW1 = 2*LENDAT + 4*HNPTS;  HNPTS = 2048; */
/*             LW2 =   LENDAT + 6*HNPTS; */
/*             MN  = min( 2*LENDAT, 2*N+1 ) */
/*             LW3 = 2*LENDAT*(2*N+1) + max( 2*LENDAT, 2*N+1 ) + */
/*                   max( MN + 6*N + 4, 2*MN + 1 ), if N > 0; */
/*             LW3 = 4*LENDAT + 5                 , if N = 0; */
/*             LW4 = max( N*N + 5*N, 6*N + 1 + min( 1,N ) ), if FLAG = 1; */
/*             LW4 = 0,                                      if FLAG = 0. */
/*             For optimum performance LDWORK should be larger. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK = LENDAT*(2*N+3), if N > 0; */
/*             LZWORK = LENDAT,         if N = 0. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the discrete --> continuous transformation cannot */
/*                   be made; */
/*             = 2:  if the system poles cannot be found; */
/*             = 3:  if the inverse system cannot be found, i.e., D is */
/*                   (close to) zero; */
/*             = 4:  if the system zeros cannot be found; */
/*             = 5:  if the state-space representation of the new */
/*                   transfer function T(s) cannot be found; */
/*             = 6:  if the continuous --> discrete transformation cannot */
/*                   be made. */

/*     METHOD */

/*     First, if the given frequency data are corresponding to a */
/*     continuous-time system, they are changed to a discrete-time */
/*     system using a bilinear transformation with a scaled alpha. */
/*     Then, the magnitude is obtained from the supplied data. */
/*     Then, the frequency data are linearly interpolated around */
/*     the unit-disc. */
/*     Then, Oppenheim and Schafer complex cepstrum method is applied */
/*     to get frequency data corresponding to a stable, minimum- */
/*     phase system. This is done in the following steps: */
/*     - Obtain LOG (magnitude) */
/*     - Obtain IFFT of the result (DG01MD SLICOT subroutine); */
/*     - halve the data at 0; */
/*     - Obtain FFT of the halved data (DG01MD SLICOT subroutine); */
/*     - Obtain EXP of the result. */
/*     Then, the new frequency data are interpolated back to the */
/*     original frequency. */
/*     Then, based on these newly obtained data, the system matrices */
/*     A, B, C, D are constructed; the very identification is */
/*     performed by Least Squares Method using DGELSY LAPACK subroutine. */
/*     If needed, a discrete-to-continuous time transformation is */
/*     applied on the system matrices by AB04MD SLICOT subroutine. */
/*     Finally, if requested, the poles and zeros of the system are */
/*     checked. If some of them have positive real parts in the */
/*     continuous-time case (or are not inside the unit disk in the */
/*     complex plane in the discrete-time case), they are exchanged with */
/*     their negatives (or reciprocals, respectively), to preserve the */
/*     frequency response, while getting a minimum phase and stable */
/*     system. This is done by SB10ZP SLICOT subroutine. */

/*     REFERENCES */

/*     [1] Oppenheim, A.V. and Schafer, R.W. */
/*         Discrete-Time Signal Processing. */
/*         Prentice-Hall Signal Processing Series, 1989. */

/*     [2] Balas, G., Doyle, J., Glover, K., Packard, A., and Smith, R. */
/*         Mu-analysis and Synthesis toolbox - User's Guide, */
/*         The Mathworks Inc., Natick, MA, USA, 1998. */

/*     CONTRIBUTORS */

/*     Asparuh Markovski, Technical University of Sofia, July 2003. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003. */
/*     A. Markovski, Technical University of Sofia, October 2003. */

/*     KEYWORDS */

/*     Bilinear transformation, frequency response, least-squares */
/*     approximation, stability. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */

/*     Test input parameters and workspace. */

    /* Parameter adjustments */
    --rfrdat;
    --ifrdat;
    --omega;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;
    --c__;
    --d__;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    pi = atan(1.) * 4.;
    pw = omega[1];
    n1 = *n + 1;
    n2 = *n + n1;

    *info = 0;
    if (*discfl != 0 && *discfl != 1) {
	*info = -1;
    } else if (*flag__ != 0 && *flag__ != 1) {
	*info = -2;
    } else if (*lendat < 2) {
	*info = -3;
    } else if (pw < 0.) {
	*info = -6;
    } else if (*n > *lendat - 1) {
	*info = -7;
    } else if (*lda < max(1,*n)) {
	*info = -9;
    } else {

	i__1 = *lendat;
	for (k = 2; k <= i__1; ++k) {
	    if (omega[k] < pw) {
		*info = -6;
	    }
	    pw = omega[k];
/* L10: */
	}

	if (*discfl == 1 && omega[*lendat] > pi) {
	    *info = -6;
	}
    }

    if (*info == 0) {

/*        Workspace. */

	lw1 = (*lendat << 1) + 8192;
	lw2 = *lendat + 12288;
/* Computing MIN */
	i__1 = *lendat << 1;
	mn = min(i__1,n2);

	if (*n > 0) {
/* Computing MAX */
	    i__1 = *lendat << 1;
/* Computing MAX */
	    i__2 = mn + *n * 6 + 4, i__3 = (mn << 1) + 1;
	    lw3 = (*lendat << 1) * n2 + max(i__1,n2) + max(i__2,i__3);
	} else {
	    lw3 = (*lendat << 2) + 5;
	}

	if (*flag__ == 0) {
	    lw4 = 0;
	} else {
/* Computing MAX */
	    i__1 = *n * *n + *n * 5, i__2 = *n * 6 + 1 + min(1,*n);
	    lw4 = max(i__1,i__2);
	}

/* Computing MAX */
	i__1 = max(2,lw1), i__1 = max(i__1,lw2), i__1 = max(i__1,lw3);
	dlwmax = max(i__1,lw4);

	if (*n > 0) {
	    clwmax = *lendat * (n2 + 2);
	} else {
	    clwmax = *lendat;
	}

	if (*ldwork < dlwmax) {
	    *info = -16;
	} else if (*lzwork < clwmax) {
	    *info = -18;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB10YD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Set tolerances. */

    tolb = dlamch_("Epsilon", (ftnlen)7);
    toll = *tol;
    if (toll <= 0.) {
	toll = (doublereal) (*lendat * *n) * 4. * tolb;
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     Workspace usage 1. */
/*     Workspace:  need  2*LENDAT + 4*HNPTS. */

    iwdomo = 1;
    iwdme = iwdomo + *lendat;
    iwymag = iwdme + 4096;
    iwmag = iwymag + 4096;

/*     Bilinear transformation. */

    if (*discfl == 0) {
	pw = sqrt(omega[1] * omega[*lendat] + sqrt(tolb));

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
	    d__1 = omega[k] / pw;
	    dwork[iwdme + k - 1] = d__1 * d__1;
	    dwork[iwdomo + k - 1] = acos((1. - dwork[iwdme + k - 1]) / (dwork[
		    iwdme + k - 1] + 1.));
/* L20: */
	}

    } else {
	dcopy_(lendat, &omega[1], &c__1, &dwork[iwdomo], &c__1);
    }

/*     Linear interpolation. */

    i__1 = *lendat;
    for (k = 1; k <= i__1; ++k) {
	dwork[iwmag + k - 1] = dlapy2_(&rfrdat[k], &ifrdat[k]);
	dwork[iwmag + k - 1] = 1. / log(10.) * log(dwork[iwmag + k - 1]);
/* L30: */
    }

    for (k = 1; k <= 2048; ++k) {
	dwork[iwdme + k - 1] = (k - 1) * pi / 2048;
	dwork[iwymag + k - 1] = 0.;

	if (dwork[iwdme + k - 1] < dwork[iwdomo]) {
	    dwork[iwymag + k - 1] = dwork[iwmag];
	} else if (dwork[iwdme + k - 1] >= dwork[iwdomo + *lendat - 1]) {
	    dwork[iwymag + k - 1] = dwork[iwmag + *lendat - 1];
	}

/* L40: */
    }

    i__1 = *lendat;
    for (i__ = 2; i__ <= i__1; ++i__) {
	p1 = dwork[iwdomo + i__ - 2] * 2048 / pi + 1.;

	ip1 = (integer) p1;
	if ((doublereal) ip1 != p1) {
	    ++ip1;
	}

	p2 = dwork[iwdomo + i__ - 1] * 2048 / pi + 1.;

	ip2 = (integer) p2;
	if ((doublereal) ip2 != p2) {
	    ++ip2;
	}

	i__2 = ip2 - 1;
	for (p = ip1; p <= i__2; ++p) {
	    rat = dwork[iwdme + p - 1] - dwork[iwdomo + i__ - 2];
	    rat /= dwork[iwdomo + i__ - 1] - dwork[iwdomo + i__ - 2];
	    dwork[iwymag + p - 1] = (1. - rat) * dwork[iwmag + i__ - 2] + rat 
		    * dwork[iwmag + i__ - 1];
/* L50: */
	}

/* L60: */
    }

    for (k = 1; k <= 2048; ++k) {
	dwork[iwymag + k - 1] = exp(log(10.) * dwork[iwymag + k - 1]);
/* L70: */
    }

/*     Duplicate data around disc. */

    for (k = 1; k <= 2048; ++k) {
	dwork[iwdme + 2048 + k - 1] = pi * 2. - dwork[iwdme + 2048 - k];
	dwork[iwymag + 2048 + k - 1] = dwork[iwymag + 2048 - k];
/* L80: */
    }

/*     Complex cepstrum to get min phase: */
/*     LOG (Magnitude) */

    for (k = 1; k <= 4096; ++k) {
	dwork[iwymag + k - 1] = log(dwork[iwymag + k - 1]) * 2.;
/* L90: */
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     Workspace usage 2. */
/*     Workspace:  need  LENDAT + 6*HNPTS. */

    iwxr = iwymag;
    iwxi = iwmag;

    for (k = 1; k <= 4096; ++k) {
	dwork[iwxi + k - 1] = 0.;
/* L100: */
    }

/*     IFFT */

    dg01md_("I", &c__4096, &dwork[iwxr], &dwork[iwxi], &info2, (ftnlen)1);

/*     Rescale, because DG01MD doesn't do it. */

    dscal_(&c__2048, &c_b21, &dwork[iwxr], &c__1);
    dscal_(&c__2048, &c_b21, &dwork[iwxi], &c__1);

/*     Halve the result at 0. */

    dwork[iwxr] /= 2.;
    dwork[iwxi] /= 2.;

/*     FFT */

    dg01md_("D", &c__2048, &dwork[iwxr], &dwork[iwxi], &info2, (ftnlen)1);

/*     Get the EXP of the result. */

    for (k = 1; k <= 1024; ++k) {
	i__1 = k - 1;
	d__1 = exp(dwork[iwxr + k - 1]);
	d__2 = cos(dwork[iwxi + k - 1]);
	d__3 = sin(dwork[iwxi + k - 1]);
	z__2.r = d__2, z__2.i = d__3;
	z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
	xhat[i__1].r = z__1.r, xhat[i__1].i = z__1.i;
	dwork[iwdme + k - 1] = dwork[iwdme + (k << 1) - 2];
/* L110: */
    }

/*     Interpolate back to original frequency data. */

    istart = 1;
    istop = *lendat;

    i__1 = *lendat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	zwork[i__2].r = 0., zwork[i__2].i = 0.;
	if (dwork[iwdomo + i__ - 1] <= dwork[iwdme]) {
	    i__2 = i__;
	    zwork[i__2].r = xhat[0].r, zwork[i__2].i = xhat[0].i;
	    istart = i__ + 1;
	} else if (dwork[iwdomo + i__ - 1] >= dwork[iwdme + 1023]) {
	    i__2 = i__;
	    zwork[i__2].r = xhat[1023].r, zwork[i__2].i = xhat[1023].i;
	    --istop;
	}
/* L120: */
    }

    i__1 = istop;
    for (i__ = istart; i__ <= i__1; ++i__) {
	ii = 1024;
L130:
	if (dwork[iwdme + ii - 1] >= dwork[iwdomo + i__ - 1]) {
	    p = ii;
	}
	--ii;
	if (ii > 0) {
	    goto L130;
	}
	rat = (dwork[iwdomo + i__ - 1] - dwork[iwdme + p - 2]) / (dwork[iwdme 
		+ p - 1] - dwork[iwdme + p - 2]);
	i__2 = i__;
	i__3 = p - 1;
	z__2.r = rat * xhat[i__3].r, z__2.i = rat * xhat[i__3].i;
	d__1 = 1. - rat;
	i__4 = p - 2;
	z__3.r = d__1 * xhat[i__4].r, z__3.i = d__1 * xhat[i__4].i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L140: */
    }

/*     CASE N > 0. */
/*     This is the only allowed case in mu-synthesis subroutines. */

    if (*n > 0) {

/*        Preparation for frequency identification. */

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Complex workspace usage 1. */
/*        Complex workspace:  need  2*LENDAT + LENDAT*(N+1). */

	iwa0 = *lendat + 1;
	iwvar = iwa0 + *lendat * n1;

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
	    if (*discfl == 0) {
		i__2 = iwvar + k - 1;
		d__1 = cos(dwork[iwdomo + k - 1]);
		d__2 = sin(dwork[iwdomo + k - 1]);
		z__1.r = d__1, z__1.i = d__2;
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
	    } else {
		i__2 = iwvar + k - 1;
		d__1 = cos(omega[k]);
		d__2 = sin(omega[k]);
		z__1.r = d__1, z__1.i = d__2;
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
	    }
/* L150: */
	}

/*        Array for DGELSY. */

	i__1 = n2;
	for (k = 1; k <= i__1; ++k) {
	    iwork[k] = 0;
/* L160: */
	}

/*        Constructing A0. */

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = iwa0 + *n * *lendat + k - 1;
	    zwork[i__2].r = 1., zwork[i__2].i = 0.;
/* L170: */
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *lendat;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = iwa0 + (*n - i__) * *lendat + k - 1;
		i__4 = iwa0 + (n1 - i__) * *lendat + k - 1;
		i__5 = iwvar + k - 1;
		z__1.r = zwork[i__4].r * zwork[i__5].r - zwork[i__4].i * 
			zwork[i__5].i, z__1.i = zwork[i__4].r * zwork[i__5].i 
			+ zwork[i__4].i * zwork[i__5].r;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L180: */
	    }
/* L190: */
	}

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Complex workspace usage 2. */
/*        Complex workspace:  need  2*LENDAT + LENDAT*(2*N+1). */

	iwbp = iwvar;
	iwab = iwbp + *lendat;

/*        Constructing BP. */

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = iwbp + k - 1;
	    i__3 = iwa0 + k - 1;
	    i__4 = k;
	    z__1.r = zwork[i__3].r * zwork[i__4].r - zwork[i__3].i * zwork[
		    i__4].i, z__1.i = zwork[i__3].r * zwork[i__4].i + zwork[
		    i__3].i * zwork[i__4].r;
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L200: */
	}

/*        Constructing AB. */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *lendat;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = iwab + (i__ - 1) * *lendat + k - 1;
		i__4 = k;
		z__2.r = -zwork[i__4].r, z__2.i = -zwork[i__4].i;
		i__5 = iwa0 + i__ * *lendat + k - 1;
		z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, 
			z__1.i = z__2.r * zwork[i__5].i + z__2.i * zwork[i__5]
			.r;
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
/* L210: */
	    }
/* L220: */
	}

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Workspace usage 3. */
/*        Workspace:  need  LW3 = 2*LENDAT*(2*N+1) + max(2*LENDAT,2*N+1). */

	iwbx = (*lendat << 1) * n2 + 1;
/* Computing MAX */
	i__1 = *lendat << 1;
	iws = iwbx + max(i__1,n2);

/*        Constructing AX. */

	i__1 = n1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *lendat;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = iwa0 + (i__ - 1) * *lendat + k - 1;
		dwork[(i__ - 1 << 1) * *lendat + k] = zwork[i__3].r;
		dwork[((i__ << 1) - 1) * *lendat + k] = d_imag(&zwork[iwa0 + (
			i__ - 1) * *lendat + k - 1]);
/* L230: */
	    }
/* L240: */
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *lendat;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = iwab + (i__ - 1) * *lendat + k - 1;
		dwork[(n1 << 1) * *lendat + (i__ - 1 << 1) * *lendat + k] = 
			zwork[i__3].r;
		dwork[(n1 << 1) * *lendat + ((i__ << 1) - 1) * *lendat + k] = 
			d_imag(&zwork[iwab + (i__ - 1) * *lendat + k - 1]);
/* L250: */
	    }
/* L260: */
	}

/*        Constructing BX. */

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = iwbp + k - 1;
	    dwork[iwbx + k - 1] = zwork[i__2].r;
	    dwork[iwbx + *lendat + k - 1] = d_imag(&zwork[iwbp + k - 1]);
/* L270: */
	}

/*        Estimating X. */
/*        Workspace:  need    LW3 + max( MN+3*(2*N+1)+1, 2*MN+1 ), */
/*                            where MN = min( 2*LENDAT, 2*N+1 ); */
/*                            prefer  larger. */

	i__1 = *lendat << 1;
	i__2 = *lendat << 1;
/* Computing MAX */
	i__4 = *lendat << 1;
	i__3 = max(i__4,n2);
	i__5 = *ldwork - iws + 1;
	dgelsy_(&i__1, &n2, &c__1, &dwork[1], &i__2, &dwork[iwbx], &i__3, &
		iwork[1], &toll, &rank, &dwork[iws], &i__5, &info2);
/* Computing MAX */
	i__1 = dlwmax, i__2 = (integer) (dwork[iws] + iws - 1);
	dlwmax = max(i__1,i__2);

/*        Constructing A matrix. */

	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    a[k + a_dim1] = -dwork[iwbx + n1 + k - 1];
/* L280: */
	}

	if (*n > 1) {
	    i__1 = *n - 1;
	    dlaset_("Full", n, &i__1, &c_b48, &c_b49, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)4);
	}

/*        Constructing B matrix. */

	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    b[k] = dwork[iwbx + n1 + k - 1] * dwork[iwbx] - dwork[iwbx + k];
/* L290: */
	}

/*        Constructing C matrix. */

	c__[1] = -1.;

	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    c__[k] = 0.;
/* L300: */
	}

/*        Constructing D matrix. */

	d__[1] = dwork[iwbx];

/*        Transform to continuous-time case, if needed. */
/*        Workspace:  need    max(1,N); */
/*                            prefer  larger. */

	if (*discfl == 0) {
	    ab04md_("D", n, &c__1, &c__1, &c_b49, &pw, &a[a_offset], lda, &b[
		    1], lda, &c__[1], &c__1, &d__[1], &c__1, &iwork[1], &
		    dwork[1], ldwork, &info2, (ftnlen)1);
	    if (info2 != 0) {
		*info = 1;
		return 0;
	    }
/* Computing MAX */
	    i__1 = dlwmax, i__2 = (integer) dwork[1];
	    dlwmax = max(i__1,i__2);
	}

/*        Make all the real parts of the poles and the zeros negative. */

	if (*flag__ == 1) {

/*           Workspace:  need    max(N*N + 5*N, 6*N + 1 + min(1,N)); */
/*                               prefer  larger. */
	    sb10zp_(discfl, n, &a[a_offset], lda, &b[1], &c__[1], &d__[1], &
		    iwork[1], &dwork[1], ldwork, info);
	    if (*info != 0) {
		return 0;
	    }
/* Computing MAX */
	    i__1 = dlwmax, i__2 = (integer) dwork[1];
	    dlwmax = max(i__1,i__2);
	}

    } else {

/*        CASE N = 0. */

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Workspace usage 4. */
/*        Workspace:  need  4*LENDAT. */

	iwbmat = (*lendat << 1) + 1;
	iws = iwbmat + (*lendat << 1);

/*        Constructing AMAT and BMAT. */

	i__1 = *lendat;
	for (k = 1; k <= i__1; ++k) {
	    dwork[k] = 1.;
	    dwork[k + *lendat] = 0.;
	    i__2 = k;
	    dwork[iwbmat + k - 1] = zwork[i__2].r;
	    dwork[iwbmat + *lendat + k - 1] = d_imag(&zwork[k]);
/* L310: */
	}

/*        Estimating D matrix. */
/*        Workspace:  need    4*LENDAT + 5; */
/*                            prefer  larger. */

	iwork[1] = 0;
	i__1 = *lendat << 1;
	i__2 = *lendat << 1;
	i__3 = *lendat << 1;
	i__4 = *ldwork - iws + 1;
	dgelsy_(&i__1, &c__1, &c__1, &dwork[1], &i__2, &dwork[iwbmat], &i__3, 
		&iwork[1], &toll, &rank, &dwork[iws], &i__4, &info2);
/* Computing MAX */
	i__1 = dlwmax, i__2 = (integer) (dwork[iws] + iws - 1);
	dlwmax = max(i__1,i__2);

	d__[1] = dwork[iwbmat];

    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

    dwork[1] = (doublereal) dlwmax;
    dwork[2] = (doublereal) clwmax;
    return 0;

/* *** Last line of SB10YD *** */
} /* sb10yd_ */

