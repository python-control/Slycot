/* MB03WD.f -- translated by f2c (version 20100827).
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

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03wd_(char *job, char *compz, integer *n, integer *p, 
	integer *ilo, integer *ihi, integer *iloz, integer *ihiz, doublereal *
	h__, integer *ldh1, integer *ldh2, doublereal *z__, integer *ldz1, 
	integer *ldz2, doublereal *wr, doublereal *wi, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen job_len, ftnlen compz_len)
{
    /* System generated locals */
    integer h_dim1, h_dim2, h_offset, z_dim1, z_dim2, z_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s, v[3];
    static integer i1, i2;
    static doublereal v1, v2, v3, h11, h12, h21, h22, h33, h44;
    static integer nh;
    static doublereal cs;
    static integer nr;
    static doublereal sn;
    static integer nz;
    static doublereal hh10, hh11, hh12, hh21, hh22, hp00, hp01, ave, hp02, 
	    hp11, hp12, hp22, h33s, h44s, tau;
    static integer itn, its;
    static doublereal ulp, tst1, h43h34, disc;
    static integer jmin, jmax;
    static doublereal unfl, ovfl;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nrow;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04py_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *);
    static logical initz, wantt, wantz;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen), dlarfx_(char *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;


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

/*     To compute the Schur decomposition and the eigenvalues of a */
/*     product of matrices, H = H_1*H_2*...*H_p, with H_1 an upper */
/*     Hessenberg matrix and H_2, ..., H_p upper triangular matrices, */
/*     without evaluating the product. Specifically, the matrices Z_i */
/*     are computed, such that */

/*             Z_1' * H_1 * Z_2 = T_1, */
/*             Z_2' * H_2 * Z_3 = T_2, */
/*                    ... */
/*             Z_p' * H_p * Z_1 = T_p, */

/*     where T_1 is in real Schur form, and T_2, ..., T_p are upper */
/*     triangular. */

/*     The routine works primarily with the Hessenberg and triangular */
/*     submatrices in rows and columns ILO to IHI, but optionally applies */
/*     the transformations to all the rows and columns of the matrices */
/*     H_i, i = 1,...,p. The transformations can be optionally */
/*     accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = 'E':  Compute the eigenvalues only; */
/*             = 'S':  Compute the factors T_1, ..., T_p of the full */
/*                     Schur form, T = T_1*T_2*...*T_p. */

/*     COMPZ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrices Z_1, ..., Z_p, as follows: */
/*             = 'N':  The matrices Z_1, ..., Z_p are not required; */
/*             = 'I':  Z_i is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z_i is returned, */
/*                     i = 1, ..., p; */
/*             = 'V':  Z_i must contain an orthogonal matrix Q_i on */
/*                     entry, and the product Q_i*Z_i is returned, */
/*                     i = 1, ..., p. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number of matrices in the product H_1*H_2*...*H_p. */
/*             P >= 1. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that all matrices H_j, j = 2, ..., p, are */
/*             already upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N, and H_1 is upper quasi-triangular in rows and */
/*             columns 1:ILO-1 and IHI+1:N, with H_1(ILO,ILO-1) = 0 */
/*             (unless ILO = 1), and H_1(IHI+1,IHI) = 0 (unless IHI = N). */
/*             The routine works primarily with the Hessenberg submatrix */
/*             in rows and columns ILO to IHI, but applies the */
/*             transformations to all the rows and columns of the */
/*             matrices H_i, i = 1,...,p, if JOB = 'S'. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     ILOZ    (input) INTEGER */
/*     IHIZ    (input) INTEGER */
/*             Specify the rows of Z to which the transformations must be */
/*             applied if COMPZ = 'I' or COMPZ = 'V'. */
/*             1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */

/*     H       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDH1,LDH2,P) */
/*             On entry, the leading N-by-N part of H(*,*,1) must contain */
/*             the upper Hessenberg matrix H_1 and the leading N-by-N */
/*             part of H(*,*,j) for j > 1 must contain the upper */
/*             triangular matrix H_j, j = 2, ..., p. */
/*             On exit, if JOB = 'S', the leading N-by-N part of H(*,*,1) */
/*             is upper quasi-triangular in rows and columns ILO:IHI, */
/*             with any 2-by-2 diagonal blocks corresponding to a pair of */
/*             complex conjugated eigenvalues, and the leading N-by-N */
/*             part of H(*,*,j) for j > 1 contains the resulting upper */
/*             triangular matrix T_j. */
/*             If JOB = 'E', the contents of H are unspecified on exit. */

/*     LDH1    INTEGER */
/*             The first leading dimension of the array H. */
/*             LDH1 >= max(1,N). */

/*     LDH2    INTEGER */
/*             The second leading dimension of the array H. */
/*             LDH2 >= max(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDZ1,LDZ2,P) */
/*             On entry, if COMPZ = 'V', the leading N-by-N-by-P part of */
/*             this array must contain the current matrix Q of */
/*             transformations accumulated by SLICOT Library routine */
/*             MB03VY. */
/*             If COMPZ = 'I', Z need not be set on entry. */
/*             On exit, if COMPZ = 'V', or COMPZ = 'I', the leading */
/*             N-by-N-by-P part of this array contains the transformation */
/*             matrices which produced the Schur form; the */
/*             transformations are applied only to the submatrices */
/*             Z_j(ILOZ:IHIZ,ILO:IHI), j = 1, ..., P. */
/*             If COMPZ = 'N', Z is not referenced. */

/*     LDZ1    INTEGER */
/*             The first leading dimension of the array Z. */
/*             LDZ1 >= 1,        if COMPZ = 'N'; */
/*             LDZ1 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'. */

/*     LDZ2    INTEGER */
/*             The second leading dimension of the array Z. */
/*             LDZ2 >= 1,        if COMPZ = 'N'; */
/*             LDZ2 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             computed eigenvalues ILO to IHI are stored in the */
/*             corresponding elements of WR and WI. If two eigenvalues */
/*             are computed as a complex conjugate pair, they are stored */
/*             in consecutive elements of WR and WI, say the i-th and */
/*             (i+1)th, with WI(i) > 0 and WI(i+1) < 0. If JOB = 'S', the */
/*             eigenvalues are stored in the same order as on the */
/*             diagonal of the Schur form returned in H. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION work array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= IHI-ILO+P-1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, ILO <= i <= IHI, the QR algorithm */
/*                   failed to compute all the eigenvalues ILO to IHI */
/*                   in a total of 30*(IHI-ILO+1) iterations; */
/*                   the elements i+1:IHI of WR and WI contain those */
/*                   eigenvalues which have been successfully computed. */

/*     METHOD */

/*     A refined version of the QR algorithm proposed in [1] and [2] is */
/*     used. The elements of the subdiagonal, diagonal, and the first */
/*     supradiagonal of current principal submatrix of H are computed */
/*     in the process. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P. */
/*         The periodic Schur decomposition: algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         Periodic Schur form and some matrix equations. */
/*         Proc. of the Symposium on the Mathematical Theory of Networks */
/*         and Systems (MTNS'93), Regensburg, Germany (U. Helmke, */
/*         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     FURTHER COMMENTS */

/*     Note that for P = 1, the LAPACK Library routine DHSEQR could be */
/*     more efficient on some computer architectures than this routine, */
/*     because DHSEQR uses a block multishift QR algorithm. */
/*     When P is large and JOB = 'S', it could be more efficient to */
/*     compute the product matrix H, and use the LAPACK Library routines. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga, */
/*     German Aerospace Center, DLR Oberpfaffenhofen, February 1999. */
/*     Partly based on the routine PSHQR by A. Varga */
/*     (DLR Oberpfaffenhofen), January 22, 1996. */

/*     REVISIONS */

/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, Hessenberg form, */
/*     orthogonal transformation, periodic systems, (periodic) Schur */
/*     form, real Schur form, similarity transformation, triangular form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    h_dim1 = *ldh1;
    h_dim2 = *ldh2;
    h_offset = 1 + h_dim1 * (1 + h_dim2);
    h__ -= h_offset;
    z_dim1 = *ldz1;
    z_dim2 = *ldz2;
    z_offset = 1 + z_dim1 * (1 + z_dim2);
    z__ -= z_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    wantt = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
    wantz = lsame_(compz, "V", (ftnlen)1, (ftnlen)1) || initz;
    *info = 0;
    if (! (wantt || lsame_(job, "E", (ftnlen)1, (ftnlen)1))) {
	*info = -1;
    } else if (! (wantz || lsame_(compz, "N", (ftnlen)1, (ftnlen)1))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*p < 1) {
	*info = -4;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -5;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -6;
    } else if (*iloz < 1 || *iloz > *ilo) {
	*info = -7;
    } else if (*ihiz < *ihi || *ihiz > *n) {
	*info = -8;
    } else if (*ldh1 < max(1,*n)) {
	*info = -10;
    } else if (*ldh2 < max(1,*n)) {
	*info = -11;
    } else if (*ldz1 < 1 || wantz && *ldz1 < *n) {
	*info = -13;
    } else if (*ldz2 < 1 || wantz && *ldz2 < *n) {
	*info = -14;
    } else if (*ldwork < *ihi - *ilo + *p - 1) {
	*info = -18;
    }
    if (*info == 0) {
	if (*ilo > 1) {
	    if (h__[*ilo + (*ilo - 1 + h_dim2) * h_dim1] != 0.) {
		*info = -5;
	    }
	} else if (*ihi < *n) {
	    if (h__[*ihi + 1 + (*ihi + h_dim2) * h_dim1] != 0.) {
		*info = -6;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03WD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     Initialize Z, if necessary. */

    if (initz) {

	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    dlaset_("Full", n, n, &c_b10, &c_b11, &z__[(j * z_dim2 + 1) * 
		    z_dim1 + 1], ldz1, (ftnlen)4);
/* L10: */
	}

    }

    nh = *ihi - *ilo + 1;

    if (nh == 1) {
	hp00 = 1.;

	i__1 = *p;
	for (j = 1; j <= i__1; ++j) {
	    hp00 *= h__[*ilo + (*ilo + j * h_dim2) * h_dim1];
/* L20: */
	}

	wr[*ilo] = hp00;
	wi[*ilo] = 0.;
	return 0;
    }

/*     Set machine-dependent constants for the stopping criterion. */
/*     If norm(H) <= sqrt(OVFL), overflow should not occur. */

    unfl = dlamch_("Safe minimum", (ftnlen)12);
    ovfl = 1. / unfl;
    dlabad_(&unfl, &ovfl);
    ulp = dlamch_("Precision", (ftnlen)9);
    smlnum = unfl * ((doublereal) nh / ulp);

/*     Set the elements in rows and columns ILO to IHI to zero below the */
/*     first subdiagonal in H(*,*,1) and below the first diagonal in */
/*     H(*,*,j), j >= 2. In the same loop, compute and store in */
/*     DWORK(NH:NH+P-2) the 1-norms of the matrices H_2, ..., H_p, to be */
/*     used later. */

    i__ = nh;
    s = ulp * (doublereal) (*n);
    if (nh > 2) {
	i__1 = nh - 2;
	i__2 = nh - 2;
	dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &h__[*ilo + 2 + (*ilo 
		+ h_dim2) * h_dim1], ldh1, (ftnlen)5);
    }

    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = nh - 1;
	i__3 = nh - 1;
	dlaset_("Lower", &i__2, &i__3, &c_b10, &c_b10, &h__[*ilo + 1 + (*ilo 
		+ j * h_dim2) * h_dim1], ldh1, (ftnlen)5);
	dwork[i__] = s * dlantr_("1-norm", "Upper", "NonUnit", &nh, &nh, &h__[
		*ilo + (*ilo + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
		ftnlen)6, (ftnlen)5, (ftnlen)7);
	++i__;
/* L30: */
    }

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

    if (wantt) {
	i1 = 1;
	i2 = *n;
    }

    if (wantz) {
	nz = *ihiz - *iloz + 1;
    }

/*     ITN is the total number of QR iterations allowed. */

    itn = nh * 30;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

    i__ = *ihi;

L40:
    l = *ilo;

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

/*     Let T = H_2*...*H_p, and H = H_1*T. Part of the currently */
/*     free locations of WR and WI are temporarily used as workspace. */

/*     WR(L:I):      the current diagonal elements of h = H(L:I,L:I); */
/*     WI(L+1:I):    the current elements of the first subdiagonal of h; */
/*     DWORK(NH-I+L:NH-1): the current elements of the first */
/*                   supradiagonal of h. */

    i__1 = itn;
    for (its = 0; its <= i__1; ++its) {

/*        Initialization: compute H(I,I) (and H(I,I-1) if I > L). */

	hp22 = 1.;
	if (i__ > l) {
	    hp12 = 0.;
	    hp11 = 1.;

	    i__2 = *p;
	    for (j = 2; j <= i__2; ++j) {
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
		hp12 = hp11 * h__[i__ - 1 + (i__ + j * h_dim2) * h_dim1] + 
			hp12 * h__[i__ + (i__ + j * h_dim2) * h_dim1];
		hp11 *= h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1];
/* L50: */
	    }

	    hh21 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp11;
	    hh22 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[i__ + (
		    i__ + h_dim2) * h_dim1] * hp22;

	    wr[i__] = hh22;
	    wi[i__] = hh21;
	} else {

	    i__2 = *p;
	    for (j = 1; j <= i__2; ++j) {
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
/* L60: */
	    }

	    wr[i__] = hp22;
	}

/*        Look for a single small subdiagonal element. */
/*        The loop also computes the needed current elements of the */
/*        diagonal and the first two supradiagonals of T, as well as */
/*        the current elements of the central tridiagonal of H. */

	i__2 = l + 1;
	for (k = i__; k >= i__2; --k) {

/*           Evaluate H(K-1,K-1), H(K-1,K) (and H(K-1,K-2) if K > L+1). */

	    hp00 = 1.;
	    hp01 = 0.;
	    if (k > l + 1) {
		hp02 = 0.;

		i__3 = *p;
		for (j = 2; j <= i__3; ++j) {
		    hp02 = hp00 * h__[k - 2 + (k + j * h_dim2) * h_dim1] + 
			    hp01 * h__[k - 1 + (k + j * h_dim2) * h_dim1] + 
			    hp02 * h__[k + (k + j * h_dim2) * h_dim1];
		    hp01 = hp00 * h__[k - 2 + (k - 1 + j * h_dim2) * h_dim1] 
			    + hp01 * h__[k - 1 + (k - 1 + j * h_dim2) * 
			    h_dim1];
		    hp00 *= h__[k - 2 + (k - 2 + j * h_dim2) * h_dim1];
/* L70: */
		}

		hh10 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp00;
		hh11 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp01 + h__[k 
			- 1 + (k - 1 + h_dim2) * h_dim1] * hp11;
		hh12 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp02 + h__[k 
			- 1 + (k - 1 + h_dim2) * h_dim1] * hp12 + h__[k - 1 + 
			(k + h_dim2) * h_dim1] * hp22;
		wi[k - 1] = hh10;
	    } else {
		hh10 = 0.;
		hh11 = h__[k - 1 + (k - 1 + h_dim2) * h_dim1] * hp11;
		hh12 = h__[k - 1 + (k - 1 + h_dim2) * h_dim1] * hp12 + h__[k 
			- 1 + (k + h_dim2) * h_dim1] * hp22;
	    }
	    wr[k - 1] = hh11;
	    dwork[nh - i__ + k - 1] = hh12;

/*           Test for a negligible subdiagonal element. */

	    tst1 = abs(hh11) + abs(hh22);
	    if (tst1 == 0.) {
		i__3 = i__ - l + 1;
		tst1 = dlanhs_("1-norm", &i__3, &h__[l + (l + h_dim2) * 
			h_dim1], ldh1, &dwork[1], (ftnlen)6);
	    }
/* Computing MAX */
	    d__1 = ulp * tst1;
	    if (abs(hh21) <= max(d__1,smlnum)) {
		goto L90;
	    }

/*           Update the values for the next cycle. */

	    hp22 = hp11;
	    hp11 = hp00;
	    hp12 = hp01;
	    hh22 = hh11;
	    hh21 = hh10;
/* L80: */
	}

L90:
	l = k;

	if (l > *ilo) {

/*           H(L,L-1) is negligible. */

	    if (wantt) {

/*              If H(L,L-1,1) is also negligible, set it to 0; otherwise, */
/*              annihilate the subdiagonal elements bottom-up, and */
/*              restore the triangular form of H(*,*,j). Since H(L,L-1) */
/*              is negligible, the second case can only appear when the */
/*              product of H(L-1,L-1,j), j >= 2, is negligible. */

		tst1 = (d__1 = h__[l - 1 + (l - 1 + h_dim2) * h_dim1], abs(
			d__1)) + (d__2 = h__[l + (l + h_dim2) * h_dim1], abs(
			d__2));
		if (tst1 == 0.) {
		    i__2 = i__ - l + 1;
		    tst1 = dlanhs_("1-norm", &i__2, &h__[l + (l + h_dim2) * 
			    h_dim1], ldh1, &dwork[1], (ftnlen)6);
		}
/* Computing MAX */
		d__2 = ulp * tst1;
		if ((d__1 = h__[l + (l - 1 + h_dim2) * h_dim1], abs(d__1)) > 
			max(d__2,smlnum)) {

		    i__2 = l;
		    for (k = i__; k >= i__2; --k) {

			i__3 = *p - 1;
			for (j = 1; j <= i__3; ++j) {

/*                       Compute G to annihilate from the right the */
/*                       (K,K-1) element of the matrix H_j. */

			    v[0] = h__[k + (k - 1 + j * h_dim2) * h_dim1];
			    dlarfg_(&c__2, &h__[k + (k + j * h_dim2) * h_dim1]
				    , v, &c__1, &tau);
			    h__[k + (k - 1 + j * h_dim2) * h_dim1] = 0.;
			    v[1] = 1.;

/*                       Apply G from the right to transform the columns */
/*                       of the matrix H_j in rows I1 to K-1. */

			    i__4 = k - i1;
			    dlarfx_("Right", &i__4, &c__2, v, &tau, &h__[i1 + 
				    (k - 1 + j * h_dim2) * h_dim1], ldh1, &
				    dwork[1], (ftnlen)5);

/*                       Apply G from the left to transform the rows of */
/*                       the matrix H_(j+1) in columns K-1 to I2. */

			    i__4 = i2 - k + 2;
			    dlarfx_("Left", &c__2, &i__4, v, &tau, &h__[k - 1 
				    + (k - 1 + (j + 1) * h_dim2) * h_dim1], 
				    ldh1, &dwork[1], (ftnlen)4);

			    if (wantz) {

/*                          Accumulate transformations in the matrix */
/*                          Z_(j+1). */

				dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*
					iloz + (k - 1 + (j + 1) * z_dim2) * 
					z_dim1], ldz1, &dwork[1], (ftnlen)5);
			    }
/* L100: */
			}

			if (k < i__) {

/*                       Compute G to annihilate from the right the */
/*                       (K+1,K) element of the matrix H_p. */

			    v[0] = h__[k + 1 + (k + *p * h_dim2) * h_dim1];
			    dlarfg_(&c__2, &h__[k + 1 + (k + 1 + *p * h_dim2) 
				    * h_dim1], v, &c__1, &tau);
			    h__[k + 1 + (k + *p * h_dim2) * h_dim1] = 0.;
			    v[1] = 1.;

/*                       Apply G from the right to transform the columns */
/*                       of the matrix H_p in rows I1 to K. */

			    i__3 = k - i1 + 1;
			    dlarfx_("Right", &i__3, &c__2, v, &tau, &h__[i1 + 
				    (k + *p * h_dim2) * h_dim1], ldh1, &dwork[
				    1], (ftnlen)5);

/*                       Apply G from the left to transform the rows of */
/*                       the matrix H_1 in columns K to I2. */

			    i__3 = i2 - k + 1;
			    dlarfx_("Left", &c__2, &i__3, v, &tau, &h__[k + (
				    k + h_dim2) * h_dim1], ldh1, &dwork[1], (
				    ftnlen)4);

			    if (wantz) {

/*                          Accumulate transformations in the matrix Z_1. */

				dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*
					iloz + (k + z_dim2) * z_dim1], ldz1, &
					dwork[1], (ftnlen)5);
			    }
			}
/* L110: */
		    }

		    h__[l + (l - 1 + *p * h_dim2) * h_dim1] = 0.;
		}
		h__[l + (l - 1 + h_dim2) * h_dim1] = 0.;
	    }
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

	if (l >= i__ - 1) {
	    goto L170;
	}

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

	if (! wantt) {
	    i1 = l;
	    i2 = i__;
	}

	if (its == 10 || its == 20) {

/*           Exceptional shift. */

	    s = (d__1 = wi[i__], abs(d__1)) + (d__2 = wi[i__ - 1], abs(d__2));
	    h44 = s * .75 + wr[i__];
	    h33 = h44;
	    h43h34 = s * -.4375 * s;
	} else {

/*           Prepare to use Francis' double shift (i.e., second degree */
/*           generalized Rayleigh quotient). */

	    h44 = wr[i__];
	    h33 = wr[i__ - 1];
	    h43h34 = wi[i__] * dwork[nh - 1];
	    disc = (h33 - h44) * .5;
	    disc = disc * disc + h43h34;
	    if (disc > 0.) {

/*              Real roots: use Wilkinson's shift twice. */

		disc = sqrt(disc);
		ave = (h33 + h44) * .5;
		if (abs(h33) - abs(h44) > 0.) {
		    h33 = h33 * h44 - h43h34;
		    h44 = h33 / (d_sign(&disc, &ave) + ave);
		} else {
		    h44 = d_sign(&disc, &ave) + ave;
		}
		h33 = h44;
		h43h34 = 0.;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l;
	for (m = i__ - 2; m >= i__2; --m) {

/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible. */

	    h11 = wr[m];
	    h12 = dwork[nh - i__ + m];
	    h21 = wi[m + 1];
	    h22 = wr[m + 1];
	    h44s = h44 - h11;
	    h33s = h33 - h11;
	    v1 = (h33s * h44s - h43h34) / h21 + h12;
	    v2 = h22 - h11 - h33s - h44s;
	    v3 = wi[m + 2];
	    s = abs(v1) + abs(v2) + abs(v3);
	    v1 /= s;
	    v2 /= s;
	    v3 /= s;
	    v[0] = v1;
	    v[1] = v2;
	    v[2] = v3;
	    if (m == l) {
		goto L130;
	    }
	    tst1 = abs(v1) * ((d__1 = wr[m - 1], abs(d__1)) + abs(h11) + abs(
		    h22));
	    if ((d__1 = wi[m], abs(d__1)) * (abs(v2) + abs(v3)) <= ulp * tst1)
		     {
		goto L130;
	    }
/* L120: */
	}

L130:

/*        Double-shift QR step. */

	i__2 = i__ - 1;
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
	    i__3 = 3, i__4 = i__ - k + 1;
	    nr = min(i__3,i__4);
/* Computing MIN */
	    i__3 = k + nr;
	    nrow = min(i__3,i__) - i1 + 1;
	    if (k > m) {
		dcopy_(&nr, &h__[k + (k - 1 + h_dim2) * h_dim1], &c__1, v, &
			c__1);
	    }
	    dlarfg_(&nr, v, &v[1], &c__1, &tau);
	    if (k > m) {
		h__[k + (k - 1 + h_dim2) * h_dim1] = v[0];
		h__[k + 1 + (k - 1 + h_dim2) * h_dim1] = 0.;
		if (k < i__ - 1) {
		    h__[k + 2 + (k - 1 + h_dim2) * h_dim1] = 0.;
		}
	    } else if (m > l) {
		h__[k + (k - 1 + h_dim2) * h_dim1] = -h__[k + (k - 1 + h_dim2)
			 * h_dim1];
	    }

/*           Apply G from the left to transform the rows of the matrix */
/*           H_1 in columns K to I2. */

	    i__3 = i2 - k + 1;
	    mb04py_("Left", &nr, &i__3, &v[1], &tau, &h__[k + (k + h_dim2) * 
		    h_dim1], ldh1, &dwork[1], (ftnlen)4);

/*           Apply G from the right to transform the columns of the */
/*           matrix H_p in rows I1 to min(K+NR,I). */

	    mb04py_("Right", &nrow, &nr, &v[1], &tau, &h__[i1 + (k + *p * 
		    h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)5);

	    if (wantz) {

/*              Accumulate transformations in the matrix Z_1. */

		mb04py_("Right", &nz, &nr, &v[1], &tau, &z__[*iloz + (k + 
			z_dim2) * z_dim1], ldz1, &dwork[1], (ftnlen)5);
	    }

	    for (j = *p; j >= 2; --j) {

/*              Apply G1 (and G2, if NR = 3) from the left to transform */
/*              the NR-by-NR submatrix of H_j in position (K,K) to upper */
/*              triangular form. */

/*              Compute G1. */

		i__3 = nr - 1;
		dcopy_(&i__3, &h__[k + 1 + (k + j * h_dim2) * h_dim1], &c__1, 
			v, &c__1);
		dlarfg_(&nr, &h__[k + (k + j * h_dim2) * h_dim1], v, &c__1, &
			tau);
		h__[k + 1 + (k + j * h_dim2) * h_dim1] = 0.;
		if (nr == 3) {
		    h__[k + 2 + (k + j * h_dim2) * h_dim1] = 0.;
		}

/*              Apply G1 from the left to transform the rows of the */
/*              matrix H_j in columns K+1 to I2. */

		i__3 = i2 - k;
		mb04py_("Left", &nr, &i__3, v, &tau, &h__[k + (k + 1 + j * 
			h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)4);

/*              Apply G1 from the right to transform the columns of the */
/*              matrix H_(j-1) in rows I1 to min(K+NR,I). */

		mb04py_("Right", &nrow, &nr, v, &tau, &h__[i1 + (k + (j - 1) *
			 h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)5);

		if (wantz) {

/*                 Accumulate transformations in the matrix Z_j. */

		    mb04py_("Right", &nz, &nr, v, &tau, &z__[*iloz + (k + j * 
			    z_dim2) * z_dim1], ldz1, &dwork[1], (ftnlen)5);
		}

		if (nr == 3) {

/*                 Compute G2. */

		    v[0] = h__[k + 2 + (k + 1 + j * h_dim2) * h_dim1];
		    dlarfg_(&c__2, &h__[k + 1 + (k + 1 + j * h_dim2) * h_dim1]
			    , v, &c__1, &tau);
		    h__[k + 2 + (k + 1 + j * h_dim2) * h_dim1] = 0.;

/*                 Apply G2 from the left to transform the rows of the */
/*                 matrix H_j in columns K+2 to I2. */

		    i__3 = i2 - k - 1;
		    mb04py_("Left", &c__2, &i__3, v, &tau, &h__[k + 1 + (k + 
			    2 + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)4);

/*                 Apply G2 from the right to transform the columns of */
/*                 the matrix H_(j-1) in rows I1 to min(K+3,I). */

		    mb04py_("Right", &nrow, &c__2, v, &tau, &h__[i1 + (k + 1 
			    + (j - 1) * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)5);

		    if (wantz) {

/*                    Accumulate transformations in the matrix Z_j. */

			mb04py_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (k 
				+ 1 + j * z_dim2) * z_dim1], ldz1, &dwork[1], 
				(ftnlen)5);
		    }
		}
/* L140: */
	    }

/* L150: */
	}

/* L160: */
    }

/*     Failure to converge in remaining number of iterations. */

    *info = i__;
    return 0;

L170:

    if (l == i__) {

/*        H(I,I-1,1) is negligible: one eigenvalue has converged. */
/*        Note that WR(I) has already been set. */

	wi[i__] = 0.;
    } else if (l == i__ - 1) {

/*        H(I-1,I-2,1) is negligible: a pair of eigenvalues have */
/*        converged. */

/*        Transform the 2-by-2 submatrix of H_1*H_2*...*H_p in position */
/*        (I-1,I-1) to standard Schur form, and compute and store its */
/*        eigenvalues. If the Schur form is not required, then the */
/*        previously stored values of a similar submatrix are used. */
/*        For real eigenvalues, a Givens transformation is used to */
/*        triangularize the submatrix. */

	if (wantt) {
	    hp22 = 1.;
	    hp12 = 0.;
	    hp11 = 1.;

	    i__1 = *p;
	    for (j = 2; j <= i__1; ++j) {
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
		hp12 = hp11 * h__[i__ - 1 + (i__ + j * h_dim2) * h_dim1] + 
			hp12 * h__[i__ + (i__ + j * h_dim2) * h_dim1];
		hp11 *= h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1];
/* L180: */
	    }

	    hh21 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp11;
	    hh22 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[i__ + (
		    i__ + h_dim2) * h_dim1] * hp22;
	    hh11 = h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1] * hp11;
	    hh12 = h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[
		    i__ - 1 + (i__ + h_dim2) * h_dim1] * hp22;
	} else {
	    hh11 = wr[i__ - 1];
	    hh12 = dwork[nh - 1];
	    hh21 = wi[i__];
	    hh22 = wr[i__];
	}

	dlanv2_(&hh11, &hh12, &hh21, &hh22, &wr[i__ - 1], &wi[i__ - 1], &wr[
		i__], &wi[i__], &cs, &sn);

	if (wantt) {

/*           Detect negligible diagonal elements in positions (I-1,I-1) */
/*           and (I,I) in H_j, J > 1. */

	    jmin = 0;
	    jmax = 0;

	    i__1 = *p;
	    for (j = 2; j <= i__1; ++j) {
		if (jmin == 0) {
		    if ((d__1 = h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1]
			    , abs(d__1)) <= dwork[nh + j - 2]) {
			jmin = j;
		    }
		}
		if ((d__1 = h__[i__ + (i__ + j * h_dim2) * h_dim1], abs(d__1))
			 <= dwork[nh + j - 2]) {
		    jmax = j;
		}
/* L190: */
	    }

	    if (jmin != 0 && jmax != 0) {

/*              Choose the shorter path if zero elements in both */
/*              (I-1,I-1) and (I,I) positions are present. */

		if (jmin - 1 <= *p - jmax + 1) {
		    jmax = 0;
		} else {
		    jmin = 0;
		}
	    }

	    if (jmin != 0) {

		i__1 = jmin - 1;
		for (j = 1; j <= i__1; ++j) {

/*                 Compute G to annihilate from the right the (I,I-1) */
/*                 element of the matrix H_j. */

		    v[0] = h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1];
		    dlarfg_(&c__2, &h__[i__ + (i__ + j * h_dim2) * h_dim1], v,
			     &c__1, &tau);
		    h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1] = 0.;
		    v[1] = 1.;

/*                 Apply G from the right to transform the columns of the */
/*                 matrix H_j in rows I1 to I-1. */

		    i__2 = i__ - i1;
		    dlarfx_("Right", &i__2, &c__2, v, &tau, &h__[i1 + (i__ - 
			    1 + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)5);

/*                 Apply G from the left to transform the rows of the */
/*                 matrix H_(j+1) in columns I-1 to I2. */

		    i__2 = i2 - i__ + 2;
		    dlarfx_("Left", &c__2, &i__2, v, &tau, &h__[i__ - 1 + (
			    i__ - 1 + (j + 1) * h_dim2) * h_dim1], ldh1, &
			    dwork[1], (ftnlen)4);

		    if (wantz) {

/*                    Accumulate transformations in the matrix Z_(j+1). */

			dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (
				i__ - 1 + (j + 1) * z_dim2) * z_dim1], ldz1, &
				dwork[1], (ftnlen)5);
		    }
/* L200: */
		}

		h__[i__ + (i__ - 1 + jmin * h_dim2) * h_dim1] = 0.;

	    } else {
		if (jmax > 0 && wi[i__ - 1] == 0.) {
		    dlartg_(&h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1], &h__[
			    i__ + (i__ - 1 + h_dim2) * h_dim1], &cs, &sn, &
			    tau);
		}

/*              Apply the transformation to H. */

		i__1 = i2 - i__ + 2;
		drot_(&i__1, &h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1], 
			ldh1, &h__[i__ + (i__ - 1 + h_dim2) * h_dim1], ldh1, &
			cs, &sn);
		i__1 = i__ - i1 + 1;
		drot_(&i__1, &h__[i1 + (i__ - 1 + *p * h_dim2) * h_dim1], &
			c__1, &h__[i1 + (i__ + *p * h_dim2) * h_dim1], &c__1, 
			&cs, &sn);
		if (wantz) {

/*                 Apply transformation to Z_1. */

		    drot_(&nz, &z__[*iloz + (i__ - 1 + z_dim2) * z_dim1], &
			    c__1, &z__[*iloz + (i__ + z_dim2) * z_dim1], &
			    c__1, &cs, &sn);
		}

/* Computing MAX */
		i__2 = 2, i__3 = jmax + 1;
		i__1 = max(i__2,i__3);
		for (j = *p; j >= i__1; --j) {

/*                 Compute G1 to annihilate from the left the (I,I-1) */
/*                 element of the matrix H_j. */

		    v[0] = h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1];
		    dlarfg_(&c__2, &h__[i__ - 1 + (i__ - 1 + j * h_dim2) * 
			    h_dim1], v, &c__1, &tau);
		    h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1] = 0.;

/*                 Apply G1 from the left to transform the rows of the */
/*                 matrix H_j in columns I to I2. */

		    i__2 = i2 - i__ + 1;
		    mb04py_("Left", &c__2, &i__2, v, &tau, &h__[i__ - 1 + (
			    i__ + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)4);

/*                 Apply G1 from the right to transform the columns of */
/*                 the matrix H_(j-1) in rows I1 to I. */

		    i__2 = i__ - i1 + 1;
		    mb04py_("Right", &i__2, &c__2, v, &tau, &h__[i1 + (i__ - 
			    1 + (j - 1) * h_dim2) * h_dim1], ldh1, &dwork[1], 
			    (ftnlen)5);

		    if (wantz) {

/*                    Apply G1 to Z_j. */

			mb04py_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (
				i__ - 1 + j * z_dim2) * z_dim1], ldz1, &dwork[
				1], (ftnlen)5);
		    }
/* L210: */
		}

		if (jmax > 0) {
		    h__[i__ + (i__ - 1 + h_dim2) * h_dim1] = 0.;
		    h__[i__ + (i__ - 1 + jmax * h_dim2) * h_dim1] = 0.;
		} else {
		    if (hh21 == 0.) {
			h__[i__ + (i__ - 1 + h_dim2) * h_dim1] = 0.;
		    }
		}
	    }
	}
    }

/*     Decrement number of remaining iterations, and return to start of */
/*     the main loop with new value of I. */

    itn -= its;
    i__ = l - 1;
    if (i__ >= *ilo) {
	goto L40;
    }

    return 0;

/* *** Last line of MB03WD *** */
} /* mb03wd_ */

