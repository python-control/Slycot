/* MA02HD.f -- translated by f2c (version 20100827).
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

logical ma02hd_(char *job, integer *m, integer *n, doublereal *diag, 
	doublereal *a, integer *lda, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    logical ret_val;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


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

/*     To check if A = DIAG*I, where I is an M-by-N matrix with ones on */
/*     the diagonal and zeros elsewhere. */

/*     FUNCTION VALUE */

/*     MA02HD  LOGICAL */
/*             The function value is set to .TRUE. if A = DIAG*I, and to */
/*             .FALSE., otherwise. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the part of the matrix A to be checked out, */
/*             as follows: */
/*             = 'U': Upper triangular/trapezoidal part; */
/*             = 'L': Lower triangular/trapezoidal part. */
/*             Otherwise:  All of the matrix A. */

/*     Input/Output Parameters */

/*     M      (input) INTEGER */
/*            The number of rows of the matrix A.  M >= 0. */

/*     N      (input) INTEGER */
/*            The number of columns of the matrix A.  N >= 0. */

/*     DIAG   (input) DOUBLE PRECISION */
/*            The scalar DIAG. */

/*     A      (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*            The leading M-by-N part of this array must contain the */
/*            matrix A.  If JOB = 'U', only the upper triangle or */
/*            trapezoid is accessed; if JOB = 'L', only the lower */
/*            triangle or trapezoid is accessed. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,M). */

/*     METHOD */

/*     The routine returns immediately after detecting a diagonal element */
/*     which differs from DIAG, or a nonzero off-diagonal element in the */
/*     searched part of A. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2003. */

/*     KEYWORDS */

/*     Elementary operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Do not check parameters, for efficiency. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (lsame_(job, "U", (ftnlen)1, (ftnlen)1)) {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = min(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (a[i__ + j * a_dim1] != 0.) {
		    ret_val = FALSE_;
		    return ret_val;
		}
/* L10: */
	    }

	    if (j <= *m) {
		if (a[j + j * a_dim1] != *diag) {
		    ret_val = FALSE_;
		    return ret_val;
		}
	    }
/* L20: */
	}

    } else if (lsame_(job, "L", (ftnlen)1, (ftnlen)1)) {

	i__1 = min(*m,*n);
	for (j = 1; j <= i__1; ++j) {
	    if (a[j + j * a_dim1] != *diag) {
		ret_val = FALSE_;
		return ret_val;
	    }

	    if (j != *m) {

/* Computing MIN */
		i__2 = j + 1;
		i__3 = *m;
		for (i__ = min(i__2,*m); i__ <= i__3; ++i__) {
		    if (a[i__ + j * a_dim1] != 0.) {
			ret_val = FALSE_;
			return ret_val;
		    }
/* L30: */
		}

	    }
/* L40: */
	}

    } else {

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
	    i__2 = j - 1;
	    i__3 = min(i__2,*m);
	    for (i__ = 1; i__ <= i__3; ++i__) {
		if (a[i__ + j * a_dim1] != 0.) {
		    ret_val = FALSE_;
		    return ret_val;
		}
/* L50: */
	    }

	    if (j <= *m) {
		if (a[j + j * a_dim1] != *diag) {
		    ret_val = FALSE_;
		    return ret_val;
		}
	    }

	    if (j < *m) {

/* Computing MIN */
		i__3 = j + 1;
		i__2 = *m;
		for (i__ = min(i__3,*m); i__ <= i__2; ++i__) {
		    if (a[i__ + j * a_dim1] != 0.) {
			ret_val = FALSE_;
			return ret_val;
		    }
/* L60: */
		}

	    }
/* L70: */
	}

    }

    ret_val = TRUE_;

    return ret_val;
/* *** Last line of MA02HD *** */
} /* ma02hd_ */

