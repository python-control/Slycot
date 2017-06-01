/* IB01OD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ib01od_(char *ctrl, integer *nobr, integer *l, 
	doublereal *sv, integer *n, doublereal *tol, integer *iwarn, integer *
	info, ftnlen ctrl_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_lg10(doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal gap;
    static integer ierr;
    static doublereal toll, rnrm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ib01oy_(integer *, integer *, integer *, 
	    doublereal *, integer *);
    static integer lnobr;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical contrl;


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

/*     To estimate the system order, based on the singular values of the */
/*     relevant part of the triangular factor of the concatenated block */
/*     Hankel matrices. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CTRL    CHARACTER*1 */
/*             Specifies whether or not the user's confirmation of the */
/*             system order estimate is desired, as follows: */
/*             = 'C':  user's confirmation; */
/*             = 'N':  no confirmation. */
/*             If  CTRL = 'C',  a reverse communication routine,  IB01OY, */
/*             is called, and, after inspecting the singular values and */
/*             system order estimate,  n,  the user may accept  n  or set */
/*             a new value. */
/*             IB01OY  is not called by the routine if CTRL = 'N'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the processed input and */
/*             output block Hankel matrices.  NOBR > 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     SV      (input) DOUBLE PRECISION array, dimension ( L*NOBR ) */
/*             The singular values of the relevant part of the triangular */
/*             factor from the QR factorization of the concatenated block */
/*             Hankel matrices. */

/*     N       (output) INTEGER */
/*             The estimated order of the system. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Absolute tolerance used for determining an estimate of */
/*             the system order. If  TOL >= 0,  the estimate is */
/*             indicated by the index of the last singular value greater */
/*             than or equal to  TOL.  (Singular values less than  TOL */
/*             are considered as zero.) When  TOL = 0,  an internally */
/*             computed default value,  TOL = NOBR*EPS*SV(1),  is used, */
/*             where  SV(1)  is the maximal singular value, and  EPS  is */
/*             the relative machine precision (see LAPACK Library routine */
/*             DLAMCH). When  TOL < 0,  the estimate is indicated by the */
/*             index of the singular value that has the largest */
/*             logarithmic gap to its successor. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 3:  all singular values were exactly zero, hence  N = 0. */
/*                   (Both input and output were identically zero.) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The singular values are compared to the given, or default TOL, and */
/*     the estimated order  n  is returned, possibly after user's */
/*     confirmation. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     August 2000. */

/*     KEYWORDS */

/*     Identification methods, multivariable systems, singular value */
/*     decomposition. */

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

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --sv;

    /* Function Body */
    contrl = lsame_(ctrl, "C", (ftnlen)1, (ftnlen)1);
    lnobr = *l * *nobr;
    *iwarn = 0;
    *info = 0;
    if (! (contrl || lsame_(ctrl, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*nobr <= 0) {
	*info = -2;
    } else if (*l <= 0) {
	*info = -3;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01OD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Set  TOL  if necessay. */

    toll = *tol;
    if (toll == 0.) {
	toll = dlamch_("Precision", (ftnlen)9) * sv[1] * (doublereal) (*nobr);
    }

/*     Obtain the system order. */

    *n = 0;
    if (sv[1] != 0.) {
	*n = *nobr;
	if (toll >= 0.) {

/*           Estimate  n  based on the tolerance  TOLL. */

	    i__1 = *nobr - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (sv[i__ + 1] < toll) {
		    *n = i__;
		    goto L30;
		}
/* L10: */
	    }
	} else {

/*           Estimate  n  based on the largest logarithmic gap between */
/*           two consecutive singular values. */

	    gap = 0.;
	    i__1 = *nobr - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		rnrm = sv[i__ + 1];
		if (rnrm != 0.) {
		    rnrm = d_lg10(&sv[i__]) - d_lg10(&rnrm);
		    if (rnrm > gap) {
			gap = rnrm;
			*n = i__;
		    }
		} else {
		    if (gap == 0.) {
			*n = i__;
		    }
		    goto L30;
		}
/* L20: */
	    }
	}
    }

L30:
    if (*n == 0) {

/*        Return with  N = 0  if all singular values are zero. */

	*iwarn = 3;
	return 0;
    }

    if (contrl) {

/*        Ask confirmation of the system order. */

	i__1 = *nobr - 1;
	ib01oy_(&lnobr, &i__1, n, &sv[1], &ierr);
    }
    return 0;

/* *** Last line of IB01OD *** */
} /* ib01od_ */

