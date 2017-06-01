/* MB02TZ.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mb02tz_(char *norm, integer *n, doublereal *hnorm, 
	doublecomplex *h__, integer *ldh, integer *ipiv, doublereal *rcond, 
	doublereal *dwork, doublecomplex *zwork, integer *info, ftnlen 
	norm_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer j;
    static doublecomplex t;
    static integer jp, ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlacon_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical onenrm;
    static doublereal hinvnm;
    extern /* Subroutine */ int zdrscl_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    static char normin[1];
    static doublereal smlnum;
    extern /* Subroutine */ int zlatrs_(char *, char *, char *, char *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);


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

/*     To estimate the reciprocal of the condition number of a complex */
/*     upper Hessenberg matrix H, in either the 1-norm or the */
/*     infinity-norm, using the LU factorization computed by MB02SZ. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     NORM    CHARACTER*1 */
/*             Specifies whether the 1-norm condition number or the */
/*             infinity-norm condition number is required: */
/*             = '1' or 'O':  1-norm; */
/*             = 'I':         Infinity-norm. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     HNORM   (input) DOUBLE PRECISION */
/*             If NORM = '1' or 'O', the 1-norm of the original matrix H. */
/*             If NORM = 'I', the infinity-norm of the original matrix H. */

/*     H       (input) COMPLEX*16 array, dimension (LDH,N) */
/*             The factors L and U from the factorization H = P*L*U */
/*             as computed by MB02SZ. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the matrix */
/*             was interchanged with row IPIV(i). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal of the condition number of the matrix H, */
/*             computed as RCOND = 1/(norm(H) * norm(inv(H))). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     ZWORK   COMPLEX*16 array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     An estimate is obtained for norm(inv(H)), and the reciprocal of */
/*     the condition number is computed as */
/*        RCOND = 1 / ( norm(H) * norm(inv(H)) ). */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N ) complex operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01FY by A.J. Laub, University of */
/*     Southern California, United States of America, May 1980. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2005. */

/*     KEYWORDS */

/*     Frequency response, Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */

/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipiv;
    --dwork;
    --zwork;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*hnorm < 0.) {
	*info = -3;
    } else if (*ldh < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB02TZ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*hnorm == 0.) {
	return 0;
    }

    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(H). */

    hinvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
L10:
    zlacon_(n, &zwork[*n + 1], &zwork[1], &hinvnm, &kase);
    if (kase != 0) {
	if (kase == kase1) {

/*           Multiply by inv(L). */

	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		jp = ipiv[j];
		i__2 = jp;
		t.r = zwork[i__2].r, t.i = zwork[i__2].i;
		if (jp != j) {
		    i__2 = jp;
		    i__3 = j;
		    zwork[i__2].r = zwork[i__3].r, zwork[i__2].i = zwork[i__3]
			    .i;
		    i__2 = j;
		    zwork[i__2].r = t.r, zwork[i__2].i = t.i;
		}
		i__2 = j + 1;
		i__3 = j + 1;
		i__4 = j + 1 + j * h_dim1;
		z__2.r = t.r * h__[i__4].r - t.i * h__[i__4].i, z__2.i = t.r *
			 h__[i__4].i + t.i * h__[i__4].r;
		z__1.r = zwork[i__3].r - z__2.r, z__1.i = zwork[i__3].i - 
			z__2.i;
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
/* L20: */
	    }

/*           Multiply by inv(U). */

	    zlatrs_("Upper", "No transpose", "Non-unit", normin, n, &h__[
		    h_offset], ldh, &zwork[1], &scale, &dwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {

/*           Multiply by inv(U'). */

	    zlatrs_("Upper", "Conjugate transpose", "Non-unit", normin, n, &
		    h__[h_offset], ldh, &zwork[1], &scale, &dwork[1], info, (
		    ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L'). */

	    for (j = *n - 1; j >= 1; --j) {
		i__1 = j;
		i__2 = j;
		d_cnjg(&z__3, &h__[j + 1 + j * h_dim1]);
		i__3 = j + 1;
		z__2.r = z__3.r * zwork[i__3].r - z__3.i * zwork[i__3].i, 
			z__2.i = z__3.r * zwork[i__3].i + z__3.i * zwork[i__3]
			.r;
		z__1.r = zwork[i__2].r - z__2.r, z__1.i = zwork[i__2].i - 
			z__2.i;
		zwork[i__1].r = z__1.r, zwork[i__1].i = z__1.i;
		jp = ipiv[j];
		if (jp != j) {
		    i__1 = jp;
		    t.r = zwork[i__1].r, t.i = zwork[i__1].i;
		    i__1 = jp;
		    i__2 = j;
		    zwork[i__1].r = zwork[i__2].r, zwork[i__1].i = zwork[i__2]
			    .i;
		    i__1 = j;
		    zwork[i__1].r = t.r, zwork[i__1].i = t.i;
		}
/* L30: */
	    }
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = izamax_(n, &zwork[1], &c__1);
	    i__1 = ix;
	    if (scale < ((d__1 = zwork[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    zwork[ix]), abs(d__2))) * smlnum || scale == 0.) {
		goto L40;
	    }
	    zdrscl_(n, &scale, &zwork[1], &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (hinvnm != 0.) {
	*rcond = 1. / hinvnm / *hnorm;
    }

L40:
    return 0;
/* *** Last line of MB02TZ *** */
} /* mb02tz_ */

