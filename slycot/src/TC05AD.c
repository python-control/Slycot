/* TC05AD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int tc05ad_(char *leri, integer *m, integer *p, 
	doublecomplex *sval, integer *index, doublereal *pcoeff, integer *
	ldpco1, integer *ldpco2, doublereal *qcoeff, integer *ldqco1, integer 
	*ldqco2, doublereal *rcond, doublecomplex *cfreqr, integer *ldcfre, 
	integer *iwork, doublereal *dwork, doublecomplex *zwork, integer *
	info, ftnlen leri_len)
{
    /* System generated locals */
    integer pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, qcoeff_dim2,
	     qcoeff_offset, cfreqr_dim1, cfreqr_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, j, k, ij, info1;
    extern /* Subroutine */ int tc01od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lleri;
    static integer mplim;
    static doublereal cnorm;
    static integer minmp, mwork, pwork;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static integer kpcoef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer maxind;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int zgecon_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), zgetrf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *), zgetrs_(char *,
	     integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer ldzwor, izwork;


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

/*     To evaluate the transfer matrix T(s) of a left polynomial matrix */
/*     representation [T(s) = inv(P(s))*Q(s)] or a right polynomial */
/*     matrix representation [T(s) = Q(s)*inv(P(s))] at any specified */
/*     complex frequency s = SVAL. */

/*     This routine will calculate the standard frequency response */
/*     matrix at frequency omega if SVAL is supplied as (0.0,omega). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether a left polynomial matrix representation */
/*             or a right polynomial matrix representation is to be used */
/*             to evaluate the transfer matrix as follows: */
/*             = 'L':  A left matrix fraction is input; */
/*             = 'R':  A right matrix fraction is input. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     SVAL    (input) COMPLEX*16 */
/*             The frequency at which the transfer matrix or the */
/*             frequency respose matrix is to be evaluated. */
/*             For a standard frequency response set the real part */
/*             of SVAL to zero. */

/*     INDEX   (input) INTEGER array, dimension (MAX(M,P)) */
/*             If LERI = 'L', INDEX(I), I = 1,2,...,P, must contain the */
/*             maximum degree of the polynomials in the I-th row of the */
/*             denominator matrix P(s) of the given left polynomial */
/*             matrix representation. */
/*             If LERI = 'R', INDEX(I), I = 1,2,...,M, must contain the */
/*             maximum degree of the polynomials in the I-th column of */
/*             the denominator matrix P(s) of the given right polynomial */
/*             matrix representation. */

/*     PCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,kpcoef), where kpcoef = MAX(INDEX(I)) + 1. */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             The leading porm-by-porm-by-kpcoef part of this array must */
/*             contain the coefficients of the denominator matrix P(s). */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if */
/*             LERI = 'L' then iorj = I, otherwise iorj = J. */
/*             Thus for LERI = 'L', P(s) = */
/*             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...). */
/*             If LERI = 'R', PCOEFF is modified by the routine but */
/*             restored on exit. */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO1 >= MAX(1,M) if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO2 >= MAX(1,M) if LERI = 'R'. */

/*     QCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,kpcoef) */
/*             If LERI = 'L' then porp = M, otherwise porp = P. */
/*             The leading porm-by-porp-by-kpcoef part of this array must */
/*             contain the coefficients of the numerator matrix Q(s). */
/*             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */
/*             If LERI = 'R', QCOEFF is modified by the routine but */
/*             restored on exit. */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,P)   if LERI = 'L', */
/*             LDQCO1 >= MAX(1,M,P) if LERI = 'R'. */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M)   if LERI = 'L', */
/*             LDQCO2 >= MAX(1,M,P) if LERI = 'R'. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The estimated reciprocal of the condition number of the */
/*             denominator matrix P(SVAL). */
/*             If RCOND is nearly zero, SVAL is approximately a system */
/*             pole. */

/*     CFREQR  (output) COMPLEX*16 array, dimension (LDCFRE,MAX(M,P)) */
/*             The leading porm-by-porp part of this array contains the */
/*             frequency response matrix T(SVAL). */

/*     LDCFRE  INTEGER */
/*             The leading dimension of array CFREQR. */
/*             LDCFRE >= MAX(1,P)   if LERI = 'L', */
/*             LDCFRE >= MAX(1,M,P) if LERI = 'R'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (liwork) */
/*             where liwork = P, if LERI = 'L', */
/*                   liwork = M, if LERI = 'R'. */

/*     DWORK   DOUBLE PRECISION array, dimension (ldwork) */
/*             where ldwork = 2*P, if LERI = 'L', */
/*                   ldwork = 2*M, if LERI = 'R'. */

/*     ZWORK   COMPLEX*16 array, dimension (lzwork), */
/*             where lzwork = P*(P+2), if LERI = 'L', */
/*                   lzwork = M*(M+2), if LERI = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if P(SVAL) is exactly or nearly singular; */
/*                   no frequency response is calculated. */

/*     METHOD */

/*     The method for a left matrix fraction will be described here; */
/*     right matrix fractions are dealt with by obtaining the dual left */
/*     fraction and calculating its frequency response (see SLICOT */
/*     Library routine TC01OD). The first step is to calculate the */
/*     complex value P(SVAL) of the denominator matrix P(s) at the */
/*     desired frequency SVAL. If P(SVAL) is approximately singular, */
/*     SVAL is approximately a pole of this system and so the frequency */
/*     response matrix T(SVAL) is not calculated; in this case, the */
/*     routine returns with the Error Indicator (INFO) set to 1. */
/*     Otherwise, the complex value Q(SVAL) of the numerator matrix Q(s) */
/*     at frequency SVAL is calculated in a similar way to P(SVAL), and */
/*     the desired response matrix T(SVAL) = inv(P(SVAL))*Q(SVAL) is */
/*     found by solving the corresponding system of complex linear */
/*     equations. */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TC01AD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     February 22, 1998 (changed the name of TC01MD). */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

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
    --index;
    pcoeff_dim1 = *ldpco1;
    pcoeff_dim2 = *ldpco2;
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
    pcoeff -= pcoeff_offset;
    qcoeff_dim1 = *ldqco1;
    qcoeff_dim2 = *ldqco2;
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
    qcoeff -= qcoeff_offset;
    cfreqr_dim1 = *ldcfre;
    cfreqr_offset = 1 + cfreqr_dim1;
    cfreqr -= cfreqr_offset;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    *info = 0;
    lleri = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
    mplim = max(*m,*p);

/*     Test the input scalar arguments. */

    if (! lleri && ! lsame_(leri, "R", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (lleri && *ldpco1 < max(1,*p) || ! lleri && *ldpco1 < max(1,*m))
	     {
	*info = -7;
    } else if (lleri && *ldpco2 < max(1,*p) || ! lleri && *ldpco2 < max(1,*m))
	     {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (lleri && *ldqco1 < max(1,*p) || ! lleri && *ldqco1 < max(i__1,*p))
		 {
	    *info = -10;
	} else if (lleri && *ldqco2 < max(1,*m) || ! lleri && *ldqco2 < max(1,
		mplim)) {
	    *info = -11;
	} else if (lleri && *ldcfre < max(1,*p) || ! lleri && *ldcfre < max(1,
		mplim)) {
	    *info = -14;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("TC05AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *p == 0) {
	*rcond = 1.;
	return 0;
    }

    if (lleri) {

/*        Initialization for left matrix fraction. */

	pwork = *p;
	mwork = *m;
    } else {

/*        Initialization for right matrix fraction: obtain dual system. */

	pwork = *m;
	mwork = *p;
	if (mplim > 1) {
	    tc01od_("R", m, p, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		    ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (
		    ftnlen)1);
	}
    }

    ldzwor = pwork;
    izwork = ldzwor * ldzwor + 1;
    maxind = 0;

    i__1 = pwork;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (index[i__] > maxind) {
	    maxind = index[i__];
	}
/* L10: */
    }

    kpcoef = maxind + 1;

/*     Calculate the complex denominator matrix P(SVAL), row by row. */

    i__1 = pwork;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ij = i__;

	i__2 = pwork;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ij;
	    i__4 = i__ + (j + pcoeff_dim2) * pcoeff_dim1;
	    z__1.r = pcoeff[i__4], z__1.i = 0.;
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
	    ij += pwork;
/* L20: */
	}

/*        Possibly non-constant row: finish evaluating it. */

	i__2 = index[i__] + 1;
	for (k = 2; k <= i__2; ++k) {

	    ij = i__;

	    i__3 = pwork;
	    for (j = 1; j <= i__3; ++j) {
		i__4 = ij;
		i__5 = ij;
		z__2.r = sval->r * zwork[i__5].r - sval->i * zwork[i__5].i, 
			z__2.i = sval->r * zwork[i__5].i + sval->i * zwork[
			i__5].r;
		i__6 = i__ + (j + k * pcoeff_dim2) * pcoeff_dim1;
		z__3.r = pcoeff[i__6], z__3.i = 0.;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
		ij += pwork;
/* L30: */
	    }

/* L40: */
	}

/* L50: */
    }

/*     Check if this P(SVAL) is singular: if so, don't compute T(SVAL). */
/*     Note that DWORK is not actually referenced in ZLANGE routine. */

    cnorm = zlange_("1-norm", &pwork, &pwork, &zwork[1], &ldzwor, &dwork[1], (
	    ftnlen)6);

    zgetrf_(&pwork, &pwork, &zwork[1], &ldzwor, &iwork[1], info);

    if (*info > 0) {

/*        Singular matrix.  Set INFO and RCOND for error return. */

	*info = 1;
	*rcond = 0.;
    } else {

/*        Estimate the reciprocal condition of P(SVAL). */
/*        Workspace: ZWORK: PWORK*PWORK + 2*PWORK, DWORK: 2*PWORK. */

	zgecon_("1-norm", &pwork, &zwork[1], &ldzwor, &cnorm, rcond, &zwork[
		izwork], &dwork[1], info, (ftnlen)6);

	if (*rcond <= dlamch_("Epsilon", (ftnlen)7)) {

/*           Nearly singular matrix.  Set INFO for error return. */

	    *info = 1;
	} else {

/*           Calculate the complex numerator matrix Q(SVAL), row by row. */

	    i__1 = pwork;
	    for (i__ = 1; i__ <= i__1; ++i__) {

		i__2 = mwork;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = i__ + j * cfreqr_dim1;
		    i__4 = i__ + (j + qcoeff_dim2) * qcoeff_dim1;
		    z__1.r = qcoeff[i__4], z__1.i = 0.;
		    cfreqr[i__3].r = z__1.r, cfreqr[i__3].i = z__1.i;
/* L60: */
		}

/*              Possibly non-constant row: finish evaluating it. */

		i__2 = index[i__] + 1;
		for (k = 2; k <= i__2; ++k) {

		    i__3 = mwork;
		    for (j = 1; j <= i__3; ++j) {
			i__4 = i__ + j * cfreqr_dim1;
			i__5 = i__ + j * cfreqr_dim1;
			z__2.r = sval->r * cfreqr[i__5].r - sval->i * cfreqr[
				i__5].i, z__2.i = sval->r * cfreqr[i__5].i + 
				sval->i * cfreqr[i__5].r;
			i__6 = i__ + (j + k * qcoeff_dim2) * qcoeff_dim1;
			z__3.r = qcoeff[i__6], z__3.i = 0.;
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
			cfreqr[i__4].r = z__1.r, cfreqr[i__4].i = z__1.i;
/* L70: */
		    }

/* L80: */
		}

/* L90: */
	    }

/*           Now calculate frequency response T(SVAL). */

	    zgetrs_("No transpose", &pwork, &mwork, &zwork[1], &ldzwor, &
		    iwork[1], &cfreqr[cfreqr_offset], ldcfre, info, (ftnlen)
		    12);
	}
    }

/*     For right matrix fraction, return to original (dual of the dual) */
/*     system. */

    if (! lleri && mplim != 1) {
	tc01od_("L", &mwork, &pwork, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, &info1, (
		ftnlen)1);

	if (*info == 0) {

/*           Also, transpose T(SVAL) here if this was successfully */
/*           calculated. */

	    minmp = min(*m,*p);

	    i__1 = mplim;
	    for (j = 1; j <= i__1; ++j) {
		if (j < minmp) {
		    i__2 = minmp - j;
		    zswap_(&i__2, &cfreqr[j + 1 + j * cfreqr_dim1], &c__1, &
			    cfreqr[j + (j + 1) * cfreqr_dim1], ldcfre);
		} else if (j > *p) {
		    zcopy_(p, &cfreqr[j * cfreqr_dim1 + 1], &c__1, &cfreqr[j 
			    + cfreqr_dim1], ldcfre);
		} else if (j > *m) {
		    zcopy_(m, &cfreqr[j + cfreqr_dim1], ldcfre, &cfreqr[j * 
			    cfreqr_dim1 + 1], &c__1);
		}
/* L100: */
	    }

	}
    }

    return 0;
/* *** Last line of TC05AD *** */
} /* tc05ad_ */

