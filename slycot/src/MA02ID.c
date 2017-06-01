/* MA02ID.f -- translated by f2c (version 20100827).
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

doublereal ma02id_(char *typ, char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, doublereal *dwork, ftnlen typ_len,
	 ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static logical lsh;
    static doublereal sum, dscl, temp, dsum, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlange_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


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

/*     To compute the value of the one norm, or the Frobenius norm, or */
/*     the infinity norm, or the element of largest absolute value */
/*     of a real skew-Hamiltonian matrix */

/*                   [  A   G  ]          T         T */
/*             X  =  [       T ],   G = -G,   Q = -Q, */
/*                   [  Q   A  ] */

/*     or of a real Hamiltonian matrix */

/*                   [  A   G  ]          T         T */
/*             X  =  [       T ],   G =  G,   Q =  Q, */
/*                   [  Q  -A  ] */

/*     where A, G and Q are real n-by-n matrices. */

/*     Note that for this kind of matrices the infinity norm is equal */
/*     to the one norm. */

/*     FUNCTION VALUE */

/*     MA02ID  DOUBLE PRECISION */
/*             The computed norm. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYP     CHARACTER*1 */
/*             Specifies the type of the input matrix X: */
/*             = 'S':         X is skew-Hamiltonian; */
/*             = 'H':         X is Hamiltonian. */

/*     NORM    CHARACTER*1 */
/*             Specifies the value to be returned in MA02ID: */
/*             = '1' or 'O':  one norm of X; */
/*             = 'F' or 'E':  Frobenius norm of X; */
/*             = 'I':         infinity norm of X; */
/*             = 'M':         max(abs(X(i,j)). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the lower triangular part of the */
/*             matrix Q and in columns 2:N+1 the upper triangular part */
/*             of the matrix G. If TYP = 'S', the parts containing the */
/*             diagonal and the first supdiagonal of this array are not */
/*             referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     Workspace */


/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             where LDWORK >= 2*N when NORM = '1', NORM = 'I' or */
/*             NORM = 'O'; otherwise, DWORK is not referenced. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLANHA). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix, skew-Hamiltonian */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    --dwork;

    /* Function Body */
    lsh = lsame_(typ, "S", (ftnlen)1, (ftnlen)1);

    if (*n == 0) {
	value = 0.;

    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1) && lsh) {

/*        Find max(abs(A(i,j))). */

	value = dlange_("MaxElement", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)10);
	if (*n > 1) {
	    i__1 = *n + 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j - 2;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = qg[i__ + j * qg_dim1], abs(
			    d__1));
		    value = max(d__2,d__3);
/* L10: */
		}
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = qg[i__ + j * qg_dim1], abs(
			    d__1));
		    value = max(d__2,d__3);
/* L20: */
		}
/* L30: */
	    }
	}

    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max( abs( A(i,j) ), abs( QG(i,j) ) ). */

/* Computing MAX */
	i__1 = *n + 1;
	d__1 = dlange_("MaxElement", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)10), d__2 = dlange_("MaxElement", n, &i__1, &qg[
		qg_offset], ldqg, &dwork[1], (ftnlen)10);
	value = max(d__1,d__2);

    } else if ((lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) && lsh) {

/*        Find the column and row sums of A (in one pass). */

	value = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[i__] = 0.;
/* L40: */
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		sum += temp;
		dwork[i__] += temp;
/* L50: */
	    }
	    dwork[*n + j] = sum;
/* L60: */
	}

/*        Compute the maximal absolute column sum. */

	i__1 = *n + 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 2;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
		dwork[i__] += temp;
		dwork[j - 1] += temp;
/* L70: */
	    }
	    if (j < *n + 1) {
		sum = dwork[*n + j];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
		    sum += temp;
		    dwork[*n + i__] += temp;
/* L80: */
		}
		value = max(value,sum);
	    }
/* L90: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = value, d__2 = dwork[i__];
	    value = max(d__1,d__2);
/* L100: */
	}

    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find the column and row sums of A (in one pass). */

	value = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dwork[i__] = 0.;
/* L110: */
	}

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		sum += temp;
		dwork[i__] += temp;
/* L120: */
	    }
	    dwork[*n + j] = sum;
/* L130: */
	}

/*        Compute the maximal absolute column sum. */

	i__1 = *n + 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j - 2;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
		dwork[i__] += temp;
		dwork[j - 1] += temp;
/* L140: */
	    }
	    if (j > 1) {
		dwork[j - 1] += (d__1 = qg[j - 1 + j * qg_dim1], abs(d__1));
	    }
	    if (j < *n + 1) {
		sum = dwork[*n + j] + (d__1 = qg[j + j * qg_dim1], abs(d__1));
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
		    sum += temp;
		    dwork[*n + i__] += temp;
/* L150: */
		}
		value = max(value,sum);
	    }
/* L160: */
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = value, d__2 = dwork[i__];
	    value = max(d__1,d__2);
/* L170: */
	}

    } else if ((lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) && lsh) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dlassq_(n, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
/* L180: */
	}

/*        Add normF(G) and normF(Q). */

	i__1 = *n + 1;
	for (j = 1; j <= i__1; ++j) {
	    if (j > 2) {
		i__2 = j - 2;
		dlassq_(&i__2, &qg[j * qg_dim1 + 1], &c__1, &scale, &sum);
	    }
	    if (j < *n) {
		i__2 = *n - j;
		dlassq_(&i__2, &qg[j + 1 + j * qg_dim1], &c__1, &scale, &sum);
	    }
/* L190: */
	}
	value = sqrt(2.) * scale * sqrt(sum);
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {
	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dlassq_(n, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
/* L200: */
	}
	dscl = 0.;
	dsum = 1.;
	i__1 = *n + 1;
	for (j = 1; j <= i__1; ++j) {
	    if (j > 1) {
		i__2 = j - 2;
		dlassq_(&i__2, &qg[j * qg_dim1 + 1], &c__1, &scale, &sum);
		dlassq_(&c__1, &qg[j - 1 + j * qg_dim1], &c__1, &dscl, &dsum);
	    }
	    if (j < *n + 1) {
		dlassq_(&c__1, &qg[j + j * qg_dim1], &c__1, &dscl, &dsum);
		i__2 = *n - j;
		dlassq_(&i__2, &qg[j + 1 + j * qg_dim1], &c__1, &scale, &sum);
	    }
/* L210: */
	}
	d__1 = sqrt(2.) * scale * sqrt(sum);
	d__2 = dscl * sqrt(dsum);
	value = dlapy2_(&d__1, &d__2);
    }

    ret_val = value;
    return ret_val;
/* *** Last line of MA02ID *** */
} /* ma02id_ */

