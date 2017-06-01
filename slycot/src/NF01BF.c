/* NF01BF.f -- translated by f2c (version 20100827).
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

static doublereal c_b3 = -1.;
static integer c__1 = 1;

/* Subroutine */ int nf01bf_(integer *iflag, integer *nfun, integer *lx, 
	integer *ipar, integer *lipar, doublereal *u, integer *ldu, 
	doublereal *y, integer *ldy, doublereal *x, integer *nfevl, 
	doublereal *e, doublereal *j, integer *ldj, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, l, m, n, nn, st, bsn;
    static doublereal err;
    static logical full;
    static integer nsmp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int nf01ad_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), nf01bd_(char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer jwork;

    /* Fortran I/O blocks */
    static cilist io___12 = { 0, 6, 0, "(' Norm of current error = ', D15.6)",
	     0 };



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

/*     This is the FCN routine for optimizing all parameters of a Wiener */
/*     system using SLICOT Library routine MD03BD. See the argument FCN */
/*     in the routine MD03BD for the description of parameters. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. CJTE is initialized to avoid the calculation of J'*e .. */
/*     .. NOUT is the unit number for printing intermediate results .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipar;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --x;
    --e;
    j_dim1 = *ldj;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --dwork;

    /* Function Body */
    l = ipar[2];
    m = ipar[5];
    n = ipar[6];
    if (l == 0) {
	nsmp = *nfun;
    } else {
	nsmp = *nfun / l;
    }

    *info = 0;
    if (*iflag == 1) {

/*        Call NF01AD to compute the output y of the Wiener system (in E) */
/*        and then the error functions (also in E). The array U must */
/*        contain the input to the linear part of the Wiener system, and */
/*        Y must contain the original output Y of the Wiener system. */
/*        IPAR(6) must contain the number of states of the linear part, n. */
/*        Workspace: need:    NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                        MAX( N*(N + L), N + M + L ) ), */
/*                                                               if M>0, */
/*                            NFUN + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                        MAX( N*(N + L), L ) ), if M=0, */
/*                            where NN = IPAR(7) (number of neurons); */
/*                   prefer:  larger. */

	i__1 = *lipar - 2;
	nf01ad_(&nsmp, &m, &l, &ipar[6], &i__1, &x[1], lx, &u[u_offset], ldu, 
		&e[1], &nsmp, &dwork[1], ldwork, info);

	i__1 = l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    daxpy_(&nsmp, &c_b3, &y[i__ * y_dim1 + 1], &c__1, &e[(i__ - 1) * 
		    nsmp + 1], &c__1);
/* L10: */
	}

/* Computing MAX */
/* Computing MAX */
	i__3 = n * (n + l), i__4 = n + m + l;
	i__1 = ipar[7] << 1, i__2 = (n + l) * (n + m) + (n << 1) + max(i__3,
		i__4);
	dwork[1] = (doublereal) (*nfun + max(i__1,i__2));

    } else if (*iflag == 2) {

/*        Call NF01BD to compute the Jacobian in a compressed form. */
/*        Workspace: need:    2*NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                          MAX( N*(N + L), N + M + L )), */
/*                                                              if M > 0, */
/*                            2*NFUN + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                          MAX( N*(N + L), L ) ), */
/*                                                              if M > 0; */
/*                   prefer:  larger. */

	i__1 = *lipar - 2;
	nf01bd_("N", &nsmp, &m, &l, &ipar[6], &i__1, &x[1], lx, &u[u_offset], 
		ldu, &e[1], &j[j_offset], ldj, &dwork[1], &dwork[1], ldwork, 
		info, (ftnlen)1);
	*nfevl = ipar[6] * (m + l + 1) + l * m;
/* Computing MAX */
/* Computing MAX */
	i__3 = n * (n + l), i__4 = n + m + l;
	i__1 = ipar[7] << 1, i__2 = (n + l) * (n + m) + (n << 1) + max(i__3,
		i__4);
	dwork[1] = (doublereal) ((*nfun << 1) + max(i__1,i__2));

    } else if (*iflag == 3) {

/*        Set the parameter LDJ, the length of the array J, and the sizes */
/*        of the workspace for FCN (IFLAG = 1 or 2), QRFACT and LMPARM. */
/*        Condition estimation (COND = 'E') is assumed in these routines. */

	st = ipar[1];
	bsn = ipar[4];
	nn = ipar[7];
	full = l <= 1 || bsn == 0;

	*ldj = *nfun;
	ipar[1] = *ldj * (bsn + st);
	if (m > 0) {
/* Computing MAX */
	    i__1 = n * (n + l), i__2 = n + m + l;
	    jwork = max(i__1,i__2);
	} else {
/* Computing MAX */
	    i__1 = n * (n + l);
	    jwork = max(i__1,l);
	}
/* Computing MAX */
	i__1 = (n + l) * (n + m) + (n << 1) + jwork, i__2 = nn << 1;
	ipar[2] = *ldj + max(i__1,i__2);
	ipar[3] = *ldj + ipar[2];
	jwork = 1;
	if (full) {
	    jwork = (*lx << 2) + 1;
	} else if (bsn > 0) {
/* Computing MAX */
	    i__1 = bsn * 3 + 1;
	    jwork = bsn + max(i__1,st);
	    if (nsmp > bsn) {
/* Computing MAX */
		i__1 = jwork, i__2 = (st << 2) + 1;
		jwork = max(i__1,i__2);
		if (nsmp < bsn << 1) {
/* Computing MAX */
		    i__1 = jwork, i__2 = (nsmp - bsn) * (l - 1);
		    jwork = max(i__1,i__2);
		}
	    }
	}
	ipar[4] = jwork;
	if (full) {
	    jwork = *lx << 2;
	} else {
	    jwork = st * (*lx - st) + (*lx << 1) + (max(bsn,st) << 1);
	}
	ipar[5] = jwork;

    } else if (*iflag == 0) {

/*        Special call for printing intermediate results. */

	err = dnrm2_(nfun, &e[1], &c__1);
	s_wsfe(&io___12);
	do_fio(&c__1, (char *)&err, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;

/* *** Last line of NF01BF *** */
} /* nf01bf_ */

