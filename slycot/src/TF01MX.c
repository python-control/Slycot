/* TF01MX.f -- translated by f2c (version 20100827).
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

static doublereal c_b4 = 0.;
static doublereal c_b8 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int tf01mx_(integer *n, integer *m, integer *p, integer *ny, 
	doublereal *s, integer *lds, doublereal *u, integer *ldu, doublereal *
	x, doublereal *y, integer *ldy, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer s_dim1, s_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, ic, nb, nf, nm, iu, np, iw, jw, iy, ns, n2m, 
	    n2p;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the output sequence of a linear time-invariant */
/*     open-loop system given by its discrete-time state-space model */
/*     with an (N+P)-by-(N+M) general system matrix S, */

/*            ( A  B ) */
/*        S = (      ) . */
/*            ( C  D ) */

/*     The initial state vector x(1) must be supplied by the user. */

/*     The input and output trajectories are stored as in the SLICOT */
/*     Library routine TF01MY. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NY      (input) INTEGER */
/*             The number of output vectors y(k) to be computed. */
/*             NY >= 0. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N+M) */
/*             The leading (N+P)-by-(N+M) part of this array must contain */
/*             the system matrix S. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N+P). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             The leading NY-by-M part of this array must contain the */
/*             input vector sequence u(k), for k = 1,2,...,NY. */
/*             Specifically, the k-th row of U must contain u(k)'. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NY). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state vector */
/*             x(1) which consists of the N initial states of the system. */
/*             On exit, this array contains the final state vector */
/*             x(NY+1) of the N states of the system at instant NY+1. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,P) */
/*             The leading NY-by-P part of this array contains the output */
/*             vector sequence y(1),y(2),...,y(NY) such that the k-th */
/*             row of Y contains y(k)' (the outputs at instant k), */
/*             for k = 1,2,...,NY. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,NY). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 0,        if MIN(N,P,NY) = 0;  otherwise, */
/*             LDWORK >= N+P,      if M = 0; */
/*             LDWORK >= 2*N+M+P,  if M > 0. */
/*             For better performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an initial state vector x(1), the output vector sequence */
/*     y(1), y(2),..., y(NY) is obtained via the formulae */

/*        ( x(k+1) )     ( x(k) ) */
/*        (        ) = S (      ) , */
/*        (  y(k)  )     ( u(k) ) */

/*     where each element y(k) is a vector of length P containing the */
/*     outputs at instant k, and k = 1,2,...,NY. */

/*     REFERENCES */

/*     [1] Luenberger, D.G. */
/*         Introduction to Dynamic Systems: Theory, Models and */
/*         Applications. */
/*         John Wiley & Sons, New York, 1979. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately (N + M) x (N + P) x NY */
/*     multiplications and additions. */

/*     FURTHER COMMENTS */

/*     The implementation exploits data locality as much as possible, */
/*     given the workspace length. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 2002. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, multivariable system, state-space model, */
/*     state-space representation, time response. */

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
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --x;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --dwork;

    /* Function Body */
    *info = 0;

/*     Test the input scalar arguments. */

    np = *n + *p;
    nm = *n + *m;
    iw = nm + np;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*ny < 0) {
	*info = -4;
    } else if (*lds < max(1,np)) {
	*info = -6;
    } else if (*ldu < max(1,*ny)) {
	*info = -8;
    } else if (*ldy < max(1,*ny)) {
	*info = -11;
    } else {
/* Computing MIN */
	i__1 = min(*n,*p);
	if (min(i__1,*ny) == 0) {
	    jw = 0;
	} else if (*m == 0) {
	    jw = np;
	} else {
	    jw = iw;
	}
	if (*ldwork < jw) {
	    *info = -13;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TF01MX", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*ny,*p) == 0) {
	return 0;
    } else if (*n == 0) {

/*        Non-dynamic system: compute the output vectors. */

	if (*m == 0) {
	    dlaset_("Full", ny, p, &c_b4, &c_b4, &y[y_offset], ldy, (ftnlen)4)
		    ;
	} else {
	    dgemm_("No transpose", "Transpose", ny, p, m, &c_b8, &u[u_offset],
		     ldu, &s[s_offset], lds, &c_b4, &y[y_offset], ldy, (
		    ftnlen)12, (ftnlen)9);
	}
	return 0;
    }

/*     Determine the block size (taken as for LAPACK routine DGETRF). */

    i__1 = max(*m,*p);
    nb = ilaenv_(&c__1, "DGETRF", " ", ny, &i__1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     Find the number of state vectors, extended with inputs (if M > 0) */
/*     and outputs, that can be accommodated in the provided workspace. */

/* Computing MIN */
    i__1 = *ldwork / jw, i__2 = nb * nb / jw, i__1 = min(i__1,i__2);
    ns = min(i__1,*ny);
    n2p = *n + np;

    if (*m == 0) {

/*        System with no inputs. */
/*        Workspace: need   N + P; */
/*                   prefer larger. */

	if (ns <= 1 || *ny * *p <= nb * nb) {
	    iy = *n + 1;

/*           LDWORK < 2*(N+P), or small problem. */
/*           One row of array Y is computed for each loop index value. */

	    i__1 = *ny;
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Compute */

/*              /x(i+1)\    /A\ */
/*              |      | =  | | * x(i). */
/*              \ y(i) /    \C/ */

		dgemv_("NoTranspose", &np, n, &c_b8, &s[s_offset], lds, &x[1],
			 &c__1, &c_b4, &dwork[1], &c__1, (ftnlen)11);
		dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);
		dcopy_(p, &dwork[iy], &c__1, &y[i__ + y_dim1], ldy);
/* L10: */
	    }

	} else {

/*           LDWORK >= 2*(N+P), and large problem. */
/*           NS rows of array Y are computed before being saved. */

	    nf = *ny / ns * ns;
	    dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

	    i__1 = nf;
	    i__2 = ns;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*              Compute the current NS extended state vectors in the */
/*              workspace: */

/*              /x(i+1)\    /A\ */
/*              |      | =  | | * x(i),  i = 1 : ns - 1. */
/*              \ y(i) /    \C/ */

		i__3 = (ns - 1) * np;
		i__4 = np;
		for (ic = 1; i__4 < 0 ? ic >= i__3 : ic <= i__3; ic += i__4) {
		    dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			    dwork[ic], &c__1, &c_b4, &dwork[ic + np], &c__1, (
			    ftnlen)12);
/* L20: */
		}

/*              Prepare the next iteration. */

		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[(ns - 1) * np + 1], &c__1, &c_b4, &dwork[1], &
			c__1, (ftnlen)12);

/*              Transpose the NS output vectors in the corresponding part */
/*              of Y (column-wise). */

		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    i__3 = ns - 1;
		    dcopy_(&i__3, &dwork[n2p + j], &np, &y[i__ + j * y_dim1], 
			    &c__1);
		    y[i__ + ns - 1 + j * y_dim1] = dwork[*n + j];
/* L30: */
		}

/* L40: */
	    }

	    ns = *ny - nf;

	    if (ns > 1) {

/*              Compute similarly the last NS output vectors. */

		i__2 = (ns - 1) * np;
		i__1 = np;
		for (ic = 1; i__1 < 0 ? ic >= i__2 : ic <= i__2; ic += i__1) {
		    dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			    dwork[ic], &c__1, &c_b4, &dwork[ic + np], &c__1, (
			    ftnlen)12);
/* L50: */
		}

		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[(ns - 1) * np + 1], &c__1, &c_b4, &dwork[1], &
			c__1, (ftnlen)12);

		i__1 = *p;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = ns - 1;
		    dcopy_(&i__2, &dwork[n2p + j], &np, &y[nf + 1 + j * 
			    y_dim1], &c__1);
		    y[nf + ns + j * y_dim1] = dwork[*n + j];
/* L60: */
		}

	    } else if (ns == 1) {

/*              Compute similarly the last NS = 1 output vectors. */

		dcopy_(n, &dwork[1], &c__1, &dwork[np + 1], &c__1);
		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[np + 1], &c__1, &c_b4, &dwork[1], &c__1, (
			ftnlen)12);
		dcopy_(p, &dwork[*n + 1], &c__1, &y[nf + 1 + y_dim1], ldy);

	    }

/*           Set the final state vector. */

	    dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);

	}

    } else {

/*        General case. */
/*        Workspace: need   2*N + M + P; */
/*                   prefer larger. */

	dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

	if (ns <= 1 || *ny * (*m + *p) <= nb * nb) {
	    iu = *n + 1;
	    jw = iu + *m;
	    iy = jw + *n;

/*           LDWORK < 2*(2*N+M+P), or small problem. */
/*           One row of array Y is computed for each loop index value. */

	    i__1 = *ny;
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Compute */

/*              /x(i+1)\    /A, B\   /x(i)\ */
/*              |      | =  |    | * |    | . */
/*              \ y(i) /    \C, D/   \u(i)/ */

		dcopy_(m, &u[i__ + u_dim1], ldu, &dwork[iu], &c__1);
		dgemv_("NoTranspose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[1], &c__1, &c_b4, &dwork[jw], &c__1, (ftnlen)11)
			;
		dcopy_(n, &dwork[jw], &c__1, &dwork[1], &c__1);
		dcopy_(p, &dwork[iy], &c__1, &y[i__ + y_dim1], ldy);
/* L70: */
	    }

	} else {

/*           LDWORK >= 2*(2*N+M+P), and large problem. */
/*           NS rows of array Y are computed before being saved. */

	    nf = *ny / ns * ns;
	    n2m = *n + nm;

	    i__1 = nf;
	    i__2 = ns;
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
		jw = 1;

/*              Compute the current NS extended state vectors in the */
/*              workspace: */

/*              /x(i+1)\    /A, B\   /x(i)\ */
/*              |      | =  |    | * |    | ,  i = 1 : ns - 1. */
/*              \ y(i) /    \C, D/   \u(i)/ */

		i__4 = *m;
		for (j = 1; j <= i__4; ++j) {
		    dcopy_(&ns, &u[i__ + j * u_dim1], &c__1, &dwork[*n + j], &
			    iw);
/* L80: */
		}

		i__4 = ns - 1;
		for (k = 1; k <= i__4; ++k) {
		    dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds,
			     &dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1,
			     (ftnlen)12);
		    jw += nm;
		    dcopy_(n, &dwork[jw], &c__1, &dwork[jw + np], &c__1);
		    jw += np;
/* L90: */
		}

/*              Prepare the next iteration. */

		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1, (
			ftnlen)12);
		dcopy_(n, &dwork[jw + nm], &c__1, &dwork[1], &c__1);

/*              Transpose the NS output vectors in the corresponding part */
/*              of Y (column-wise). */

		i__4 = *p;
		for (j = 1; j <= i__4; ++j) {
		    dcopy_(&ns, &dwork[n2m + j], &iw, &y[i__ + j * y_dim1], &
			    c__1);
/* L100: */
		}

/* L110: */
	    }

	    ns = *ny - nf;

	    if (ns > 1) {
		jw = 1;

/*              Compute similarly the last NS output vectors. */

		i__2 = *m;
		for (j = 1; j <= i__2; ++j) {
		    dcopy_(&ns, &u[nf + 1 + j * u_dim1], &c__1, &dwork[*n + j]
			    , &iw);
/* L120: */
		}

		i__2 = ns - 1;
		for (k = 1; k <= i__2; ++k) {
		    dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds,
			     &dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1,
			     (ftnlen)12);
		    jw += nm;
		    dcopy_(n, &dwork[jw], &c__1, &dwork[jw + np], &c__1);
		    jw += np;
/* L130: */
		}

		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1, (
			ftnlen)12);
		dcopy_(n, &dwork[jw + nm], &c__1, &dwork[1], &c__1);

		i__2 = *p;
		for (j = 1; j <= i__2; ++j) {
		    dcopy_(&ns, &dwork[n2m + j], &iw, &y[nf + 1 + j * y_dim1],
			     &c__1);
/* L140: */
		}

	    } else if (ns == 1) {

/*              Compute similarly the last NS = 1 output vectors. */

		dcopy_(n, &dwork[1], &c__1, &dwork[np + 1], &c__1);
		dcopy_(m, &u[nf + 1 + u_dim1], ldu, &dwork[n2p + 1], &c__1);
		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[np + 1], &c__1, &c_b4, &dwork[1], &c__1, (
			ftnlen)12);
		dcopy_(p, &dwork[*n + 1], &c__1, &y[nf + 1 + y_dim1], ldy);

	    }

	}

/*        Set the final state vector. */

	dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);

    }

    return 0;
/* *** Last line of TF01MX *** */
} /* tf01mx_ */

