/* IB01OY.f -- translated by f2c (version 20100827).
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
static integer c__3 = 3;

/* Subroutine */ int ib01oy_(integer *ns, integer *nmax, integer *n, 
	doublereal *sv, integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void), s_rsle(cilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_rsle(void);

    /* Local variables */
    static integer i__;
    static char ans[1];
    static logical yes;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, "(/' Singular values (in descending or"
	    "der) used',                  ' to estimate the system order:', /"
	    "/                               (5D15.8) )", 0 };
    static cilist io___3 = { 0, 6, 0, "(/' Estimated order of the system,  n"
	    " = ', I5 )", 0 };
    static cilist io___4 = { 0, 6, 0, "(/' Do you want this value of  n  to "
	    "be used',                    ' to determine the system matrices?"
	    "' )", 0 };
    static cilist io___5 = { 0, 6, 0, "(/'  Type \"yes\" or \"no\":  ' )", 0 }
	    ;
    static cilist io___6 = { 0, 5, 0, "( A )", 0 };
    static cilist io___9 = { 0, 6, 0, "(/' n  should be less than or equal',"
	    "                             ' to ', I5 )", 0 };
    static cilist io___10 = { 0, 6, 0, "( ' (It may be useful to restart',  "
	    "                              ' with a larger tolerance.)' )", 0 }
	    ;
    static cilist io___11 = { 0, 6, 0, "(/' Enter the desired value of n (n "
	    "<= ', I5,                     ');  n = ' )", 0 };
    static cilist io___12 = { 0, 5, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, "(/' n  should be larger than zero.' )"
	    , 0 };
    static cilist io___14 = { 0, 6, 0, "(/' n  should be less than or equal "
	    "to ',                    I5 )", 0 };



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

/*     To ask for user's confirmation of the system order found by */
/*     SLICOT Library routine IB01OD. This routine may be modified, */
/*     but its interface must be preserved. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NS      (input) INTEGER */
/*             The number of singular values.  NS > 0. */

/*     NMAX    (input) INTEGER */
/*             The maximum value of the system order.  0 <= NMAX <= NS. */

/*     N       (input/output) INTEGER */
/*             On entry, the estimate of the system order computed by */
/*             IB01OD routine.  0 <= N <= NS. */
/*             On exit, the user's estimate of the system order, which */
/*             could be identical with the input value of  N. */
/*             Note that the output value of  N  should be less than */
/*             or equal to  NMAX. */

/*     SV      (input) DOUBLE PRECISION array, dimension ( NS ) */
/*             The singular values, in descending order, used for */
/*             determining the system order. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Identification, parameter estimation, singular values, structure */
/*     identification. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*        INTRMN is the unit number for the (terminal) input device. */
/*        OUTRMN is the unit number for the (terminal) output device. */
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

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    --sv;

    /* Function Body */
    *info = 0;
    if (*ns <= 0) {
	*info = -1;
    } else if (*nmax < 0 || *nmax > *ns) {
	*info = -2;
    } else if (*n < 0 || *n > *ns) {
	*info = -3;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01OY", &i__1, (ftnlen)6);
	return 0;
    }

    s_wsfe(&io___1);
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&sv[i__], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();
    s_wsfe(&io___3);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___4);
    e_wsfe();

L10:
    s_wsfe(&io___5);
    e_wsfe();
    s_rsfe(&io___6);
    do_fio(&c__1, ans, (ftnlen)1);
    e_rsfe();
    yes = lsame_(ans, "Y", (ftnlen)1, (ftnlen)1);
    if (yes) {
	if (*n <= *nmax) {

/*              The value of n is adequate and has been confirmed. */

	    return 0;
	} else {

/*              The estimated value of n is not acceptable. */

	    s_wsfe(&io___9);
	    do_fio(&c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___10);
	    e_wsfe();
	    goto L20;
	}

    } else if (lsame_(ans, "N", (ftnlen)1, (ftnlen)1)) {
	goto L20;
    } else {

/*           Wrong answer should be re-entered. */

	goto L10;
    }

/*     Enter the desired value of n. */

L20:
    s_wsfe(&io___11);
    do_fio(&c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
    e_wsfe();
    s_rsle(&io___12);
    do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_rsle();
    if (*n < 0) {

/*           The specified value of n is not acceptable. */

	s_wsfe(&io___13);
	e_wsfe();
	goto L20;
    } else if (*n > *nmax) {

/*           The specified value of n is not acceptable. */

	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
	e_wsfe();
	goto L20;
    }

    return 0;

/* *** Last line of IB01OY *** */
} /* ib01oy_ */

