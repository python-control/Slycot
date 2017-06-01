/* NF01AD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int nf01ad_(integer *nsmp, integer *m, integer *l, integer *
	ipar, integer *lipar, doublereal *x, integer *lx, doublereal *u, 
	integer *ldu, doublereal *y, integer *ldy, doublereal *dwork, integer 
	*ldwork, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer n, z__, ac, bd, nn, ix, jw, ldac, lths, nths;
    extern /* Subroutine */ int nf01ay_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), tf01mx_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    tb01vy_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);


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

/*     To calculate the output y of the Wiener system */

/*        x(t+1) = A*x(t) + B*u(t) */
/*        z(t)   = C*x(t) + D*u(t), */

/*        y(t)   = f(z(t),wb(1:L)), */

/*     where t = 1, 2, ..., NSMP, and f is a nonlinear function, */
/*     evaluated by the SLICOT Library routine NF01AY. The parameter */
/*     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ), */
/*     where wb(i), i = 1:L, correspond to the nonlinear part, theta */
/*     corresponds to the linear part, and the notation is fully */
/*     described below. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     M       (input) INTEGER */
/*             The length of each input sample.  M >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample.  L >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             IPAR(1)  must contain the order of the linear part, */
/*                      referred to as N below.  N >= 0. */
/*             IPAR(2)  must contain the number of neurons for the */
/*                      nonlinear part, referred to as NN below. */
/*                      NN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of IPAR.  LIPAR >= 2. */

/*     X       (input) DOUBLE PRECISION array, dimension (LX) */
/*             The parameter vector, partitioned as */
/*             X = (wb(1), ..., wb(L), theta), where the vectors */
/*             wb(i), of length NN*(L+2)+1, are parameters for the */
/*             static nonlinearity, which is simulated by the */
/*             SLICOT Library routine NF01AY. See the documentation of */
/*             NF01AY for further details. The vector theta, of length */
/*             N*(M + L + 1) + L*M, represents the matrices A, B, C, */
/*             D and x(1), and it can be retrieved from these matrices */
/*             by SLICOT Library routine TB01VD and retranslated by */
/*             TB01VY. */

/*     LX      (input) INTEGER */
/*             The length of the array X. */
/*             LX >= ( NN*(L+2)+1 )*L + N*(M + L + 1) + L*M. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU, M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             set of input samples, */
/*             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ). */

/*     LDU     INTEGER */
/*             The leading dimension of the array U.  LDU >= MAX(1,NSMP). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array contains the */
/*             simulated output. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                     MAX( N*(N + L), N + M + L ) ) */
/*                                                              if M > 0; */
/*             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                     MAX( N*(N + L), L ) ),   if M = 0. */
/*             A larger value of LDWORK could improve the efficiency. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*     METHOD */

/*     BLAS routines are used for the matrix-vector multiplications and */
/*     the routine NF01AY is called for the calculation of the nonlinear */
/*     function. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Mar. 2001, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Dec. 2001. */

/*     KEYWORDS */

/*     Nonlinear system, output normal form, simulation, state-space */
/*     representation, Wiener system. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipar;
    --x;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*nsmp < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*l < 0) {
	*info = -3;
    } else if (*lipar < 2) {
	*info = -5;
    } else {

	n = ipar[1];
	nn = ipar[2];
	ldac = n + *l;
	nths = (nn * (*l + 2) + 1) * *l;
	lths = n * (*m + *l + 1) + *l * *m;

	if (n < 0 || nn < 0) {
	    *info = -4;
	} else if (*lx < nths + lths) {
	    *info = -7;
	} else if (*ldu < max(1,*nsmp)) {
	    *info = -9;
	} else if (*ldy < max(1,*nsmp)) {
	    *info = -11;
	} else {
	    if (*m > 0) {
/* Computing MAX */
		i__1 = n * ldac, i__2 = n + *m + *l;
		jw = max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = n * ldac;
		jw = max(i__1,*l);
	    }
/* Computing MAX */
	    i__1 = nn << 1, i__2 = ldac * (n + *m) + (n << 1) + jw;
	    if (*ldwork < *nsmp * *l + max(i__1,i__2)) {
		*info = -13;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*nsmp,*l) == 0) {
	return 0;
    }

/*     Compute the output of the linear part. */
/*     Workspace: need   NSMP*L + (N + L)*(N + M) + N + N*(N + L + 1). */
/*     (NSMP*L locations are reserved for the output of the linear part.) */

    z__ = 1;
    ac = z__ + *nsmp * *l;
    bd = ac + ldac * n;
    ix = bd + ldac * *m;
    jw = ix + n;

    i__1 = *ldwork - jw + 1;
    tb01vy_("Apply", &n, m, l, &x[nths + 1], &lths, &dwork[ac], &ldac, &dwork[
	    bd], &ldac, &dwork[ac + n], &ldac, &dwork[bd + n], &ldac, &dwork[
	    ix], &dwork[jw], &i__1, info, (ftnlen)5);

/*     Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, if M>0; */
/*                       NSMP*L + (N + L)*N + 2*N + L,           if M=0; */
/*                prefer larger. */

    i__1 = *ldwork - jw + 1;
    tf01mx_(&n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &dwork[ix], 
	    &dwork[z__], nsmp, &dwork[jw], &i__1, info);

/*     Simulate the static nonlinearity. */
/*     Workspace: need   NSMP*L + 2*NN; */
/*                prefer larger. */

    jw = ac;
    i__1 = *lipar - 1;
    i__2 = *ldwork - jw + 1;
    nf01ay_(nsmp, l, l, &ipar[2], &i__1, &x[1], &nths, &dwork[z__], nsmp, &y[
	    y_offset], ldy, &dwork[jw], &i__2, info);

    return 0;

/* *** Last line of NF01AD *** */
} /* nf01ad_ */

