/* NF01BD.f -- translated by f2c (version 20100827).
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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b17 = 1.;
static doublereal c_b19 = 0.;

/* Subroutine */ int nf01bd_(char *cjte, integer *nsmp, integer *m, integer *
	l, integer *ipar, integer *lipar, doublereal *x, integer *lx, 
	doublereal *u, integer *ldu, doublereal *e, doublereal *j, integer *
	ldj, doublereal *jte, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen cjte_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, u_dim1, u_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal h__;
    static integer i__, k, n, z__, ac, bd, nn, ix, iy, jw, bsn;
    static doublereal eps;
    static integer ldac, kcol, lpar;
    static logical wjte;
    static integer lths, nsml, nths;
    extern /* Subroutine */ int nf01ad_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), nf01ay_(integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    nf01by_(char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), tf01mx_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *), tb01vy_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal parsav;


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

/*     To calculate the Jacobian dy/dX of the Wiener system */

/*        x(t+1) = A*x(t) + B*u(t) */
/*        z(t)   = C*x(t) + D*u(t), */

/*        y(t,i) = sum( ws(k, i)*f(w(k, i)*z(t) + b(k,i)) ) + b(k+1,i), */

/*     where t = 1, 2, ...,  NSMP, */
/*           i = 1, 2, ...,  L, */
/*           k = 1, 2, ...,  NN. */

/*     NN is arbitrary eligible and has to be provided in IPAR(2), and */
/*     X = ( wb(1), ..., wb(L), theta ) is described below. */

/*     Denoting y(j) = y(1:NSMP,j), the Jacobian J has the block form */

/*       dy(1)/dwb(1)       0         .....       0         dy(1)/dtheta */
/*            0        dy(2)/dwb(2)   .....       0         dy(2)/dtheta */
/*          .....         .....       .....     .....          ..... */
/*            0           .....         0    dy(L)/dwb(L)   dy(L)/dtheta */

/*     but it will be returned without the zero blocks, in the form */

/*     dy(1)/dwb(1)    dy(1)/dtheta */
/*                  ... */
/*     dy(L)/dwb(L)    dy(L)/dtheta. */

/*     dy(i)/dwb(i) depends on f and is calculated by the routine NF01BY; */
/*     dy(i)/dtheta is computed by a forward-difference approximation. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CJTE    CHARACTER*1 */
/*             Specifies whether the matrix-vector product J'*e should be */
/*             computed or not, as follows: */
/*             = 'C' :  compute J'*e; */
/*             = 'N' :  do not compute J'*e. */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     M       (input) INTEGER */
/*             The length of each input sample.  M >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample.  L >= 0. */

/*     IPAR    (input/output) INTEGER array, dimension (LIPAR) */
/*             On entry, the first entries of this array must contain */
/*             the integer parameters needed; specifically, */
/*             IPAR(1)  must contain the order of the linear part, N; */
/*                      actually, N = abs(IPAR(1)), since setting */
/*                      IPAR(1) < 0 has a special meaning (see below); */
/*             IPAR(2)  must contain the number of neurons for the */
/*                      nonlinear part, NN, NN >= 0. */
/*             On exit, if IPAR(1) < 0 on entry, then no computations are */
/*             performed, except the needed tests on input parameters, */
/*             but the following values are returned: */
/*             IPAR(1) contains the length of the array J, LJ; */
/*             LDJ     contains the leading dimension of array J. */
/*             Otherwise, IPAR(1) and LDJ are unchanged on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 2. */

/*     X       (input) DOUBLE PRECISION array, dimension (LX) */
/*             The leading LPAR entries of this array must contain the */
/*             set of system parameters, where */
/*                LPAR = (NN*(L + 2) + 1)*L + N*(M + L + 1) + L*M. */
/*             X has the form (wb(1), ..., wb(L), theta), where the */
/*             vectors wb(i) have the structure */
/*              (w(1,1), ..., w(1,L), ..., w(NN,1), ..., w(NN,L), */
/*                ws(1), ..., ws(NN), b(1), ..., b(NN+1) ), */
/*             and the vector theta represents the matrices A, B, C, D */
/*             and x(1), and it can be retrieved from these matrices */
/*             by SLICOT Library routine TB01VD and retranslated by */
/*             TB01VY. */

/*     LX      (input) INTEGER */
/*             The length of X. */
/*             LX >= (NN*(L + 2) + 1)*L + N*(M + L + 1) + L*M. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU, M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             set of input samples, */
/*             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NSMP). */

/*     E       (input) DOUBLE PRECISION array, dimension (NSMP*L) */
/*             If CJTE = 'C', this array must contain a vector e, which */
/*             will be premultiplied with J', e = vec( Y - y ), where */
/*             Y is set of output samples, and vec denotes the */
/*             concatenation of the columns of a matrix. */
/*             If CJTE = 'N', this array is not referenced. */

/*     J       (output) DOUBLE PRECISION array, dimension (LDJ, *) */
/*             The leading NSMP*L-by-NCOLJ part of this array contains */
/*             the Jacobian of the error function stored in a compressed */
/*             form, as described above, where */
/*             NCOLJ = NN*(L + 2) + 1 + N*(M + L + 1) + L*M. */

/*     LDJ     INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NSMP*L). */
/*             Note that LDJ is an input parameter, except for */
/*             IPAR(1) < 0 on entry, when it is an output parameter. */

/*     JTE     (output) DOUBLE PRECISION array, dimension (LPAR) */
/*             If CJTE = 'C', this array contains the matrix-vector */
/*             product J'*e. */
/*             If CJTE = 'N', this array is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                       MAX( N*(N + L), N + M + L ) ) */
/*                                                              if M > 0; */
/*             LDWORK >= 2*NSMP*L + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                       MAX( N*(N + L), L ) ), if M = 0. */
/*             A larger value of LDWORK could improve the efficiency. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     BLAS routines are used for the matrix-vector multiplications, and */
/*     the SLICOT Library routine TB01VY is called for the conversion of */
/*     the output normal form parameters to an LTI-system; the routine */
/*     NF01AD is then used for the simulation of the system with given */
/*     parameters, and the routine NF01BY is called for the (analytically */
/*     performed) calculation of the parts referring to the parameters */
/*     of the static nonlinearity. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Mar. 2001, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Dec. 2001. */

/*     KEYWORDS */

/*     Jacobian matrix, nonlinear system, output normal form, simulation, */
/*     state-space representation, Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. EPSFCN is related to the error in computing the functions .. */
/*     .. For EPSFCN = 0.0D0, the square root of the machine precision */
/*     .. is used for finite difference approximation of the derivatives. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
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
    --e;
    j_dim1 = *ldj;
    j_offset = 1 + j_dim1;
    j -= j_offset;
    --jte;
    --dwork;

    /* Function Body */
    n = ipar[1];
    nn = ipar[2];
    bsn = nn * (*l + 2) + 1;
    nsml = *nsmp * *l;
    nths = bsn * *l;
    lths = n * (*m + *l + 1) + *l * *m;
    lpar = nths + lths;
    wjte = lsame_(cjte, "C", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

    *info = 0;
    if (! (wjte || lsame_(cjte, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (*nsmp < 0) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (nn < 0) {
	*info = -5;
    } else if (*lipar < 2) {
	*info = -6;
    } else if (ipar[1] < 0) {
	if (*info != 0) {
	    i__1 = -(*info);
	    xerbla_("NF01BD", &i__1, (ftnlen)6);
	} else {
	    ipar[1] = nsml * (abs(n) * (*m + *l + 1) + *l * *m + bsn);
	    *ldj = max(1,nsml);
	}
	return 0;
    } else if (*lx < lpar) {
	*info = -8;
    } else if (*ldu < max(1,*nsmp)) {
	*info = -10;
    } else if (*ldj < max(1,nsml)) {
	*info = -13;
    } else {
	ldac = n + *l;
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
	if (*ldwork < (nsml << 1) + max(i__1,i__2)) {
	    *info = -16;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("NF01BD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*nsmp,*l) == 0) {
	if (wjte && lpar >= 1) {
	    jte[1] = 0.;
	    dcopy_(&lpar, &jte[1], &c__0, &jte[1], &c__1);
	}
	return 0;
    }

/*     Compute the output of the linear part. */
/*     Workspace: need  2*NSMP*L + (N + L)*(N + M) + N + N*(N + L + 1). */
/*     (2*NSMP*L locations are reserved for computing two times the */
/*     output of the linear part.) */

    iy = 1;
    z__ = iy + nsml;
    ac = z__ + nsml;
    bd = ac + ldac * n;
    ix = bd + ldac * *m;
    jw = ix + n;

    i__1 = *ldwork - jw + 1;
    tb01vy_("Apply", &n, m, l, &x[nths + 1], &lths, &dwork[ac], &ldac, &dwork[
	    bd], &ldac, &dwork[ac + n], &ldac, &dwork[bd + n], &ldac, &dwork[
	    ix], &dwork[jw], &i__1, info, (ftnlen)5);

/*     Workspace: need   2*NSMP*L + (N + L)*(N + M) + 3*N + M + L, */
/*                                                             if M > 0; */
/*                       2*NSMP*L + (N + L)*N + 2*N + L,       if M = 0; */
/*                prefer larger. */

    i__1 = *ldwork - jw + 1;
    tf01mx_(&n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &dwork[ix], 
	    &dwork[z__], nsmp, &dwork[jw], &i__1, info);

/*     Fill the blocks dy(i)/dwb(i) and the corresponding parts of JTE, */
/*     if needed. */

    jw = ac;
    if (wjte) {

	i__1 = *l - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__2 = *lipar - 1;
	    i__3 = *ldwork - jw + 1;
	    nf01by_(cjte, nsmp, l, &c__1, &ipar[2], &i__2, &x[i__ * bsn + 1], 
		    &bsn, &dwork[z__], nsmp, &e[i__ * *nsmp + 1], &j[i__ * *
		    nsmp + 1 + j_dim1], ldj, &jte[i__ * bsn + 1], &dwork[jw], 
		    &i__3, info, (ftnlen)1);
/* L10: */
	}

    } else {

	i__1 = *l - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__2 = *lipar - 1;
	    i__3 = *ldwork - jw + 1;
	    nf01by_(cjte, nsmp, l, &c__1, &ipar[2], &i__2, &x[i__ * bsn + 1], 
		    &bsn, &dwork[z__], nsmp, &dwork[1], &j[i__ * *nsmp + 1 + 
		    j_dim1], ldj, &dwork[1], &dwork[jw], &i__3, info, (ftnlen)
		    1);
/* L20: */
	}

    }

/*     Compute the output of the system with unchanged parameters. */
/*     Workspace: need   2*NSMP*L + 2*NN; */
/*                prefer larger. */

    i__1 = *lipar - 1;
    i__2 = *ldwork - jw + 1;
    nf01ay_(nsmp, l, l, &ipar[2], &i__1, &x[1], &nths, &dwork[z__], nsmp, &
	    dwork[iy], nsmp, &dwork[jw], &i__2, info);

/*     Compute dy/dtheta numerically by forward-difference approximation. */
/*     Workspace: need   2*NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                       MAX( N*(N + L), N + M + L ) ), */
/*                                                              if M > 0; */
/*                       2*NSMP*L + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                       MAX( N*(N + L), L ) ), if M = 0; */
/*                prefer larger. */

    jw = z__;
/* Computing MAX */
    d__1 = 0., d__2 = dlamch_("Epsilon", (ftnlen)7);
    eps = sqrt((max(d__1,d__2)));

    i__1 = lpar;
    for (k = nths + 1; k <= i__1; ++k) {
	kcol = k - nths + bsn;
	parsav = x[k];
	if (parsav == 0.) {
	    h__ = eps;
	} else {
	    h__ = eps * abs(parsav);
	}
	x[k] += h__;
	i__2 = *ldwork - jw + 1;
	nf01ad_(nsmp, m, l, &ipar[1], lipar, &x[1], &lpar, &u[u_offset], ldu, 
		&j[kcol * j_dim1 + 1], nsmp, &dwork[jw], &i__2, info);
	x[k] = parsav;

	i__2 = nsml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    j[i__ + kcol * j_dim1] = (j[i__ + kcol * j_dim1] - dwork[i__]) / 
		    h__;
/* L30: */
	}

/* L40: */
    }

    if (wjte) {

/*        Compute the last part of J'e in JTE. */

	dgemv_("Transpose", &nsml, &lths, &c_b17, &j[(bsn + 1) * j_dim1 + 1], 
		ldj, &e[1], &c__1, &c_b19, &jte[nths + 1], &c__1, (ftnlen)9);
    }

    return 0;

/* *** Last line of NF01BD *** */
} /* nf01bd_ */

