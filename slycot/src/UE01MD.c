/* UE01MD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__8 = 8;

integer ue01md_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, ftnlen name_len, ftnlen opts_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char c1[1], c2[2], c3[1];
    static integer ic, nb, iz, nx;
    static logical cname, sname;
    static integer nbmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char subnam[6];


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

/*     To provide an extension of the LAPACK routine ILAENV to */
/*     machine-specific parameters for SLICOT routines. */

/*     The default values in this version aim to give good performance on */
/*     a wide range of computers. For optimal performance, however, the */
/*     user is advised to modify this routine. Note that an optimized */
/*     BLAS is a crucial prerequisite for any speed gains. For further */
/*     details, see ILAENV. */

/*     FUNCTION VALUE */

/*     UE01MD  INTEGER */
/*             The function value set according to ISPEC. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     ISPEC   (input) INTEGER */
/*             Specifies the parameter to be returned as the value of */
/*             UE01MD, as follows: */
/*             = 1: the optimal blocksize; if the returned value is 1, an */
/*                  unblocked algorithm will give the best performance; */
/*             = 2: the minimum block size for which the block routine */
/*                  should be used; if the usable block size is less than */
/*                  this value, an unblocked routine should be used; */
/*             = 3: the crossover point (in a block routine, for N less */
/*                  than this value, an unblocked routine should be used) */
/*             = 4: the number of shifts, used in the product eigenvalue */
/*                  routine; */
/*             = 8: the crossover point for the multishift QR method for */
/*                  product eigenvalue problems. */

/*     NAME    (input) CHARACTER*(*) */
/*             The name of the calling subroutine, in either upper case */
/*             or lower case. */

/*     OPTS    (input) CHARACTER*(*) */
/*             The character options to the subroutine NAME, concatenated */
/*             into a single character string. */

/*     N1      (input) INTEGER */
/*     N2      (input) INTEGER */
/*     N3      (input) INTEGER */
/*             Problem dimensions for the subroutine NAME; these may not */
/*             all be required. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine ILAHAP). */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */

/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

    if (*ispec == 1 || *ispec == 2 || *ispec == 3) {

/*        Convert NAME to upper case if the first character is lower */
/*        case. */

	ret_val = 1;
	s_copy(subnam, name__, (ftnlen)6, name_len);
	ic = *(unsigned char *)subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122) {

/*           ASCII character set. */

	    if (ic >= 97 && ic <= 122) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 97 && ic <= 122) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
/* L10: */
		}
	    }

	} else if (iz == 233 || iz == 169) {

/*           EBCDIC character set. */

	    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
		    && ic <= 169) {
		*(unsigned char *)subnam = (char) (ic + 64);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
			    ic >= 162 && ic <= 169) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		    }
/* L20: */
		}
	    }

	} else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128. */

	    if (ic >= 225 && ic <= 250) {
		*(unsigned char *)subnam = (char) (ic - 32);
		for (i__ = 2; i__ <= 6; ++i__) {
		    ic = *(unsigned char *)&subnam[i__ - 1];
		    if (ic >= 225 && ic <= 250) {
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		    }
/* L30: */
		}
	    }
	}

	*(unsigned char *)c1 = *(unsigned char *)subnam;
	sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
	cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
	if (! (cname || sname)) {
	    return ret_val;
	}
	s_copy(c2, subnam + 3, (ftnlen)2, (ftnlen)2);
	*(unsigned char *)c3 = *(unsigned char *)&subnam[5];

	if (*ispec == 1) {

/*           Block size. */

	    nb = 1;
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
		    nb = ilaenv_(&c__1, "DGEQRF", " ", n1, n2, &c_n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		} else if (*(unsigned char *)c3 == 'T') {
		    nb = ilaenv_(&c__1, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 4;
		}
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
		    nb = ilaenv_(&c__1, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		}
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'D') {
		    nb = ilaenv_(&c__1, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		} else if (*(unsigned char *)c3 == 'B') {
		    nb = ilaenv_(&c__1, "DORMQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NB = ILAENV( 1, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2 */
/* *             END IF */
	    }
	    ret_val = nb;
	} else if (*ispec == 2) {

/*           Minimum block size. */

	    nbmin = 2;
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", n1, n2, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1) / 2;
		    nbmin = max(i__1,i__2);
		} else if (*(unsigned char *)c3 == 'T') {
/* Computing MAX */
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n1, n2, n1,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 4;
		    nbmin = max(i__1,i__2);
		}
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n1, n2, n1,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 4;
		    nbmin = max(i__1,i__2);
		}
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'D') {
/* Computing MAX */
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", n1, n2, n3,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 2;
		    nbmin = max(i__1,i__2);
		} else if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQR", " ", n1, n2, n3,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 2;
		    nbmin = max(i__1,i__2);
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N1, N2, N1, */
/* *   $                                    -1 ) / 4 ) */
/* *             END IF */
	    }
	    ret_val = nbmin;
	} else if (*ispec == 3) {

/*           Crossover point. */

	    nx = 0;
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
		    nx = ilaenv_(&c__3, "DGEQRF", " ", n1, n2, &c_n1, &c_n1, (
			    ftnlen)6, (ftnlen)1);
		} else if (*(unsigned char *)c3 == 'T') {
		    nx = ilaenv_(&c__3, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		}
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'B') {
		    nx = ilaenv_(&c__3, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
		}
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
		if (*(unsigned char *)c3 == 'D') {
		    nx = ilaenv_(&c__3, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1);
		} else if (*(unsigned char *)c3 == 'B') {
		    nx = ilaenv_(&c__3, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1);
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NX = ILAENV( 3, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2 */
/* *             END IF */
	    }
	    ret_val = nx;
	}
    } else if (*ispec == 4) {

/*        Number of shifts (used by MB03XP). */

	ret_val = ilaenv_(&c__4, "DHSEQR", opts, n1, n2, n3, &c_n1, (ftnlen)6,
		 opts_len);
    } else if (*ispec == 8) {

/*        Crossover point for multishift (used by MB03XP). */

	ret_val = ilaenv_(&c__8, "DHSEQR", opts, n1, n2, n3, &c_n1, (ftnlen)6,
		 opts_len);
    } else {

/*        Invalid value for ISPEC. */

	ret_val = -1;
    }
    return ret_val;
/* *** Last line of UE01MD *** */
} /* ue01md_ */

