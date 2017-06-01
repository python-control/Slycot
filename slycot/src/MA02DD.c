/* MA02DD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ma02dd_(char *job, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *ap, ftnlen job_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer j, ij;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical luplo;


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

/*     To pack/unpack the upper or lower triangle of a symmetric matrix. */
/*     The packed matrix is stored column-wise in the one-dimensional */
/*     array AP. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies whether the matrix should be packed or unpacked, */
/*             as follows: */
/*             = 'P':  The matrix should be packed; */
/*             = 'U':  The matrix should be unpacked. */

/*     UPLO    CHARACTER*1 */
/*             Specifies the part of the matrix to be packed/unpacked, */
/*             as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             This array is an input parameter if JOB = 'P', and an */
/*             output parameter if JOB = 'U'. */
/*             On entry, if JOB = 'P', the leading N-by-N upper */
/*             triangular part (if UPLO = 'U'), or lower triangular part */
/*             (if UPLO = 'L'), of this array must contain the */
/*             corresponding upper or lower triangle of the symmetric */
/*             matrix A, and the other strictly triangular part is not */
/*             referenced. */
/*             On exit, if JOB = 'U', the leading N-by-N upper triangular */
/*             part (if UPLO = 'U'), or lower triangular part (if */
/*             UPLO = 'L'), of this array contains the corresponding */
/*             upper or lower triangle of the symmetric matrix A; the */
/*             other strictly triangular part is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     AP      (output or input) DOUBLE PRECISION array, dimension */
/*             (N*(N+1)/2) */
/*             This array is an output parameter if JOB = 'P', and an */
/*             input parameter if JOB = 'U'. */
/*             On entry, if JOB = 'U', the leading N*(N+1)/2 elements of */
/*             this array must contain the upper (if UPLO = 'U') or lower */
/*             (if UPLO = 'L') triangle of the symmetric matrix A, packed */
/*             column-wise. That is, the elements are stored in the order */
/*             11, 12, 22, ..., 1n, 2n, 3n, ..., nn,      if UPLO = 'U'; */
/*             11, 21, 31, ..., n1, 22, 32, ..., n2, ..., if UPLO = 'L'. */
/*             On exit, if JOB = 'P', the leading N*(N+1)/2 elements of */
/*             this array contain the upper (if UPLO = 'U') or lower */
/*             (if UPLO = 'L') triangle of the symmetric matrix A, packed */
/*             column-wise, as described above. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked for errors. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ap;

    /* Function Body */
    luplo = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
    ij = 1;
    if (lsame_(job, "P", (ftnlen)1, (ftnlen)1)) {
	if (luplo) {

/*           Pack the lower triangle of A. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		dcopy_(&i__2, &a[j + j * a_dim1], &c__1, &ap[ij], &c__1);
		ij = ij + *n - j + 1;
/* L20: */
	    }

	} else {

/*           Pack the upper triangle of A. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(&j, &a[j * a_dim1 + 1], &c__1, &ap[ij], &c__1);
		ij += j;
/* L40: */
	    }

	}
    } else {
	if (luplo) {

/*           Unpack the lower triangle of A. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		dcopy_(&i__2, &ap[ij], &c__1, &a[j + j * a_dim1], &c__1);
		ij = ij + *n - j + 1;
/* L60: */
	    }

	} else {

/*           Unpack the upper triangle of A. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		dcopy_(&j, &ap[ij], &c__1, &a[j * a_dim1 + 1], &c__1);
		ij += j;
/* L80: */
	    }

	}
    }

    return 0;
/* *** Last line of MA02DD *** */
} /* ma02dd_ */

