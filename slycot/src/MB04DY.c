/* MB04DY.f -- translated by f2c (version 20100827).
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
static doublereal c_b19 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb04dy_(char *jobscl, integer *n, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *d__, doublereal *
	dwork, integer *info, ftnlen jobscl_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal y;
    static integer ihi;
    static doublereal ofl;
    static integer ilo;
    static doublereal ufl, eps, rho, tau, base, anrm;
    static logical none;
    static integer ierr;
    static doublereal gnrm;
    static logical norm;
    static doublereal qnrm;
    static logical symp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sfmin, sfmax;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebal_(
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);


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

/*     To perform a symplectic scaling on the Hamiltonian matrix */

/*              ( A    G  ) */
/*          H = (       T ),                                          (1) */
/*              ( Q   -A  ) */

/*     i.e., perform either the symplectic scaling transformation */

/*                                   -1 */
/*                 ( A'   G'  )   ( D   0 ) ( A   G  ) ( D  0   ) */
/*          H' <-- (        T ) = (       ) (      T ) (     -1 ),    (2) */
/*                 ( Q'  -A'  )   ( 0   D ) ( Q  -A  ) ( 0  D   ) */

/*     where D is a diagonal scaling matrix, or the symplectic norm */
/*     scaling transformation */

/*                  ( A''   G''  )    1  (   A   G/tau ) */
/*          H'' <-- (          T ) = --- (           T ),             (3) */
/*                  ( Q''  -A''  )   tau ( tau Q   -A  ) */

/*     where tau is a real scalar.  Note that if tau is not equal to 1, */
/*     then (3) is NOT a similarity transformation.  The eigenvalues */
/*     of H are then tau times the eigenvalues of H''. */

/*     For symplectic scaling (2), D is chosen to give the rows and */
/*     columns of A' approximately equal 1-norms and to give Q' and G' */
/*     approximately equal norms.  (See METHOD below for details.) For */
/*     norm scaling, tau = MAX(1, ||A||, ||G||, ||Q||) where ||.|| */
/*     denotes the 1-norm (column sum norm). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBSCL  CHARACTER*1 */
/*             Indicates which scaling strategy is used, as follows: */
/*             = 'S'       :  do the symplectic scaling (2); */
/*             = '1' or 'O':  do the 1-norm scaling (3); */
/*             = 'N'       :  do nothing; set INFO and return. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On input, if JOBSCL <> 'N', the leading N-by-N part of */
/*             this array must contain the upper left block A of the */
/*             Hamiltonian matrix H in (1). */
/*             On output, if JOBSCL <> 'N', the leading N-by-N part of */
/*             this array contains the leading N-by-N part of the scaled */
/*             Hamiltonian matrix H' in (2) or H'' in (3), depending on */
/*             the setting of JOBSCL. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if JOBSCL <> 'N'; */
/*             LDA >= 1,        if JOBSCL =  'N'. */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQG,N+1) */
/*             On input, if JOBSCL <> 'N', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangle of the lower left symmetric block Q of the */
/*             Hamiltonian matrix H in (1), and the N-by-N upper */
/*             triangular part of the submatrix in the columns 2 to N+1 */
/*             of this array must contain the upper triangle of the upper */
/*             right symmetric block G of H in (1). */
/*             So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j) */
/*             and G(i,j) = G(j,i) is stored in QG(j,i+1). */
/*             On output, if JOBSCL <> 'N', the leading N-by-N lower */
/*             triangular part of this array contains the lower triangle */
/*             of the lower left symmetric block Q' or Q'', and the */
/*             N-by-N upper triangular part of the submatrix in the */
/*             columns 2 to N+1 of this array contains the upper triangle */
/*             of the upper right symmetric block G' or G'' of the scaled */
/*             Hamiltonian matrix H' in (2) or H'' in (3), depending on */
/*             the setting of JOBSCL. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG. */
/*             LDQG >= MAX(1,N), if JOBSCL <> 'N'; */
/*             LDQG >= 1,        if JOBSCL =  'N'. */

/*     D       (output) DOUBLE PRECISION array, dimension (nd) */
/*             If JOBSCL = 'S', then nd = N and D contains the diagonal */
/*             elements of the diagonal scaling matrix in (2). */
/*             If JOBSCL = '1' or 'O', then nd = 1 and D(1) is set to tau */
/*             from (3). In this case, no other elements of D are */
/*             referenced. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, then the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     1. Symplectic scaling (JOBSCL = 'S'): */

/*     First, LAPACK subroutine DGEBAL is used to equilibrate the 1-norms */
/*     of the rows and columns of A using a diagonal scaling matrix D_A. */
/*     Then, H is similarily transformed by the symplectic diagonal */
/*     matrix D1 = diag(D_A,D_A**(-1)).  Next, the off-diagonal blocks of */
/*     the resulting Hamiltonian matrix are equilibrated in the 1-norm */
/*     using the symplectic diagonal matrix D2 of the form */

/*                 ( I/rho    0   ) */
/*            D2 = (              ) */
/*                 (   0    rho*I ) */

/*     where rho is a real scalar. Thus, in (2), D = D1*D2. */

/*     2. Norm scaling (JOBSCL = '1' or 'O'): */

/*     The norm of the matrices A and G of (1) is reduced by setting */
/*     A := A/tau  and  G := G/(tau**2) where tau is the power of the */
/*     base of the arithmetic closest to MAX(1, ||A||, ||G||, ||Q||) and */
/*     ||.|| denotes the 1-norm. */

/*     REFERENCES */

/*     [1] Benner, P., Byers, R., and Barth, E. */
/*         Fortran 77 Subroutines for Computing the Eigenvalues of */
/*         Hamiltonian Matrices. I: The Square-Reduced Method. */
/*         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000. */

/*     NUMERICAL ASPECTS */

/*     For symplectic scaling, the complexity of the used algorithms is */
/*     hard to estimate and depends upon how well the rows and columns of */
/*     A in (1) are equilibrated.  In one sweep, each row/column of A is */
/*     scaled once, i.e., the cost of one sweep is N**2 multiplications. */
/*     Usually, 3-6 sweeps are enough to equilibrate the norms of the */
/*     rows and columns of a matrix.  Roundoff errors are possible as */
/*     LAPACK routine DGEBAL does NOT use powers of the machine base for */
/*     scaling. The second stage (equilibrating ||G|| and ||Q||) requires */
/*     N**2 multiplications. */
/*     For norm scaling, 3*N**2 + O(N) multiplications are required and */
/*     NO rounding errors occur as all multiplications are performed with */
/*     powers of the machine base. */

/*     CONTRIBUTOR */

/*     P. Benner, Universitaet Bremen, Germany, and */
/*     R. Byers, University of Kansas, Lawrence, USA. */
/*     Aug. 1998, routine DHABL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998, SLICOT Library version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2009. */

/*     KEYWORDS */

/*     Balancing, Hamiltonian matrix, norms, symplectic similarity */
/*     transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */

/*     .. Scalar Arguments .. */
/*    .. */
/*    .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    qg_dim1 = *ldqg;
    qg_offset = 1 + qg_dim1;
    qg -= qg_offset;
    --d__;
    --dwork;

    /* Function Body */
    *info = 0;
    symp = lsame_(jobscl, "S", (ftnlen)1, (ftnlen)1);
    norm = lsame_(jobscl, "1", (ftnlen)1, (ftnlen)1) || lsame_(jobscl, "O", (
	    ftnlen)1, (ftnlen)1);
    none = lsame_(jobscl, "N", (ftnlen)1, (ftnlen)1);

    if (! symp && ! norm && ! none) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < 1 || ! none && *lda < *n) {
	*info = -4;
    } else if (*ldqg < 1 || ! none && *ldqg < *n) {
	*info = -6;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB04DY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || none) {
	return 0;
    }

/*     Set some machine dependant constants. */

    base = dlamch_("Base", (ftnlen)4);
    eps = dlamch_("Precision", (ftnlen)9);
    ufl = dlamch_("Safe minimum", (ftnlen)12);
    ofl = 1. / ufl;
    dlabad_(&ufl, &ofl);
    sfmax = eps / base / ufl;
    sfmin = 1. / sfmax;

    if (norm) {

/*        Compute norms. */

	anrm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);
	gnrm = dlansy_("1-norm", "Upper", n, &qg[(qg_dim1 << 1) + 1], ldqg, &
		dwork[1], (ftnlen)6, (ftnlen)5);
	qnrm = dlansy_("1-norm", "Lower", n, &qg[qg_offset], ldqg, &dwork[1], 
		(ftnlen)6, (ftnlen)5);
/* Computing MAX */
	d__1 = max(1.,anrm), d__1 = max(d__1,gnrm);
	y = max(d__1,qnrm);
	tau = 1.;

/*        WHILE ( TAU < Y ) DO */
L10:
	if (tau < y && tau < sqrt(sfmax)) {
	    tau *= base;
	    goto L10;
	}
/*        END WHILE 10 */
	if (tau > 1.) {
	    if ((d__1 = tau / base - y, abs(d__1)) < (d__2 = tau - y, abs(
		    d__2))) {
		tau /= base;
	    }
	    dlascl_("General", &c__0, &c__0, &tau, &c_b19, n, n, &a[a_offset],
		     lda, &ierr, (ftnlen)7);
	    dlascl_("Upper", &c__0, &c__0, &tau, &c_b19, n, n, &qg[(qg_dim1 <<
		     1) + 1], ldqg, &ierr, (ftnlen)5);
	    dlascl_("Upper", &c__0, &c__0, &tau, &c_b19, n, n, &qg[(qg_dim1 <<
		     1) + 1], ldqg, &ierr, (ftnlen)5);
	}

	d__[1] = tau;

    } else {
	dgebal_("Scale", n, &a[a_offset], lda, &ilo, &ihi, &d__[1], &ierr, (
		ftnlen)5);

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

	    i__2 = *n;
	    for (i__ = j; i__ <= i__2; ++i__) {
		qg[i__ + j * qg_dim1] = qg[i__ + j * qg_dim1] * d__[j] * d__[
			i__];
/* L20: */
	    }

/* L30: */
	}

	i__1 = *n + 1;
	for (j = 2; j <= i__1; ++j) {

	    i__2 = j - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		qg[i__ + j * qg_dim1] = qg[i__ + j * qg_dim1] / d__[j - 1] / 
			d__[i__];
/* L40: */
	    }

/* L50: */
	}

	gnrm = dlansy_("1-norm", "Upper", n, &qg[(qg_dim1 << 1) + 1], ldqg, &
		dwork[1], (ftnlen)6, (ftnlen)5);
	qnrm = dlansy_("1-norm", "Lower", n, &qg[qg_offset], ldqg, &dwork[1], 
		(ftnlen)6, (ftnlen)5);
	if (gnrm == 0.) {
	    if (qnrm == 0.) {
		rho = 1.;
	    } else {
		rho = sfmax;
	    }
	} else if (qnrm == 0.) {
	    rho = sfmin;
	} else {
	    rho = sqrt(qnrm) / sqrt(gnrm);
	}

	dlascl_("Lower", &c__0, &c__0, &rho, &c_b19, n, n, &qg[qg_offset], 
		ldqg, &ierr, (ftnlen)5);
	dlascl_("Upper", &c__0, &c__0, &c_b19, &rho, n, n, &qg[(qg_dim1 << 1) 
		+ 1], ldqg, &ierr, (ftnlen)5);
	d__1 = sqrt(rho);
	drscl_(n, &d__1, &d__[1], &c__1);
    }

    return 0;
/*     *** Last line of MB04DY *** */
} /* mb04dy_ */

