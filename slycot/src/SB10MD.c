/* SB10MD.f -- translated by f2c (version 20100827).
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
static doublereal c_b15 = 0.;
static doublereal c_b25 = 1.;
static integer c__0 = 0;

/* Subroutine */ int sb10md_(integer *nc, integer *mp, integer *lendat, 
	integer *f, integer *ord, integer *mnb, integer *nblock, integer *
	itype, doublereal *qutol, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *omega, integer *totord, doublereal *ad, integer *
	ldad, doublereal *bd, integer *ldbd, doublereal *cd, integer *ldcd, 
	doublereal *dd, integer *lddd, doublereal *mju, integer *iwork, 
	integer *liwork, doublereal *dwork, integer *ldwork, doublecomplex *
	zwork, integer *lzwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ad_dim1, ad_offset, b_dim1, b_offset, bd_dim1, 
	    bd_offset, c_dim1, c_offset, cd_dim1, cd_offset, d_dim1, d_offset,
	     dd_dim1, dd_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, k, w, ic, ii, mn, lw1, lw2, lw3, lw4, iwb, lwa, lwb;
    static doublereal rqe, tol;
    static integer iwx;
    static doublereal mod1, mod2, maqe;
    static integer iwad, iwbd, iwcd, cord, iwdd;
    static doublereal meqe, rcnd;
    static doublecomplex freq;
    static integer lord, info2;
    extern /* Subroutine */ int ab13md_(char *, integer *, doublecomplex *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen), 
	    tb05ad_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, doublecomplex *, integer *
	    , doublereal *, doublereal *, doublecomplex *, integer *, integer 
	    *, doublereal *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *), sb10yd_(integer *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublecomplex *, integer *, integer *);
    static char inita[1];
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer icwrk, idwrk;
    static doublereal toler;
    static char baleig[1];
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer iwifrd, lcsize, ldsize, clwmax, dlwmax, iwgjom, iwrfrd, 
	    maxcwr, maxwrk;


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

/*     To perform the D-step in the D-K iteration. It handles */
/*     continuous-time case. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NC      (input) INTEGER */
/*             The order of the matrix A.  NC >= 0. */

/*     MP      (input) INTEGER */
/*             The order of the matrix D.  MP >= 0. */

/*     LENDAT  (input) INTEGER */
/*             The length of the vector OMEGA.  LENDAT >= 2. */

/*     F       (input) INTEGER */
/*             The number of the measurements and controls, i.e., */
/*             the size of the block I_f in the D-scaling system. */
/*             F >= 0. */

/*     ORD     (input/output) INTEGER */
/*             The MAX order of EACH block in the fitting procedure. */
/*             ORD <= LENDAT-1. */
/*             On exit, if ORD < 1 then ORD = 1. */

/*     MNB     (input) INTEGER */
/*             The number of diagonal blocks in the block structure of */
/*             the uncertainty, and the length of the vectors NBLOCK */
/*             and ITYPE.  1 <= MNB <= MP. */

/*     NBLOCK  (input) INTEGER array, dimension (MNB) */
/*             The vector of length MNB containing the block structure */
/*             of the uncertainty. NBLOCK(I), I = 1:MNB, is the size of */
/*             each block. */

/*     ITYPE   (input) INTEGER array, dimension (MNB) */
/*             The vector of length MNB indicating the type of each */
/*             block. */
/*             For I = 1 : MNB, */
/*             ITYPE(I) = 1 indicates that the corresponding block is a */
/*             real block. IN THIS CASE ONLY MJU(JW) WILL BE ESTIMATED */
/*             CORRECTLY, BUT NOT D(S)! */
/*             ITYPE(I) = 2 indicates that the corresponding block is a */
/*             complex block. THIS IS THE ONLY ALLOWED VALUE NOW! */
/*             NBLOCK(I) must be equal to 1 if ITYPE(I) is equal to 1. */

/*     QUTOL   (input) DOUBLE PRECISION */
/*             The acceptable mean relative error between the D(jw) and */
/*             the frequency responce of the estimated block */
/*             [ADi,BDi;CDi,DDi]. When it is reached, the result is */
/*             taken as good enough. */
/*             A good value is QUTOL = 2.0. */
/*             If QUTOL < 0 then only mju(jw) is being estimated, */
/*             not D(s). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,NC) */
/*             On entry, the leading NC-by-NC part of this array must */
/*             contain the A matrix of the closed-loop system. */
/*             On exit, if MP > 0, the leading NC-by-NC part of this */
/*             array contains an upper Hessenberg matrix similar to A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,NC). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,MP) */
/*             On entry, the leading NC-by-MP part of this array must */
/*             contain the B matrix of the closed-loop system. */
/*             On exit, the leading NC-by-MP part of this array contains */
/*             the transformed B matrix corresponding to the Hessenberg */
/*             form of A. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,NC). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,NC) */
/*             On entry, the leading MP-by-NC part of this array must */
/*             contain the C matrix of the closed-loop system. */
/*             On exit, the leading MP-by-NC part of this array contains */
/*             the transformed C matrix corresponding to the Hessenberg */
/*             form of A. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,MP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,MP) */
/*             The leading MP-by-MP part of this array must contain the */
/*             D matrix of the closed-loop system. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,MP). */

/*     OMEGA   (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The vector with the frequencies. */

/*     TOTORD  (output) INTEGER */
/*             The TOTAL order of the D-scaling system. */
/*             TOTORD is set to zero, if QUTOL < 0. */

/*     AD      (output) DOUBLE PRECISION array, dimension (LDAD,MP*ORD) */
/*             The leading TOTORD-by-TOTORD part of this array contains */
/*             the A matrix of the D-scaling system. */
/*             Not referenced if QUTOL < 0. */

/*     LDAD    INTEGER */
/*             The leading dimension of the array AD. */
/*             LDAD >= MAX(1,MP*ORD), if QUTOL >= 0; */
/*             LDAD >= 1,             if QUTOL <  0. */

/*     BD      (output) DOUBLE PRECISION array, dimension (LDBD,MP+F) */
/*             The leading TOTORD-by-(MP+F) part of this array contains */
/*             the B matrix of the D-scaling system. */
/*             Not referenced if QUTOL < 0. */

/*     LDBD    INTEGER */
/*             The leading dimension of the array BD. */
/*             LDBD >= MAX(1,MP*ORD), if QUTOL >= 0; */
/*             LDBD >= 1,             if QUTOL <  0. */

/*     CD      (output) DOUBLE PRECISION array, dimension (LDCD,MP*ORD) */
/*             The leading (MP+F)-by-TOTORD part of this array contains */
/*             the C matrix of the D-scaling system. */
/*             Not referenced if QUTOL < 0. */

/*     LDCD    INTEGER */
/*             The leading dimension of the array CD. */
/*             LDCD >= MAX(1,MP+F), if QUTOL >= 0; */
/*             LDCD >= 1,           if QUTOL <  0. */

/*     DD      (output) DOUBLE PRECISION array, dimension (LDDD,MP+F) */
/*             The leading (MP+F)-by-(MP+F) part of this array contains */
/*             the D matrix of the D-scaling system. */
/*             Not referenced if QUTOL < 0. */

/*     LDDD    INTEGER */
/*             The leading dimension of the array DD. */
/*             LDDD >= MAX(1,MP+F), if QUTOL >= 0; */
/*             LDDD >= 1,           if QUTOL <  0. */

/*     MJU     (output) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The vector with the upper bound of the structured */
/*             singular value (mju) for each frequency in OMEGA. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */

/*     LIWORK  INTEGER */
/*             The length of the array IWORK. */
/*             LIWORK >= MAX( NC, 4*MNB-2, MP, 2*ORD+1 ), if QUTOL >= 0; */
/*             LIWORK >= MAX( NC, 4*MNB-2, MP ),          if QUTOL <  0. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the optimal value of LZWORK, */
/*             and DWORK(3) returns an estimate of the minimum reciprocal */
/*             of the condition numbers (with respect to inversion) of */
/*             the generated Hessenberg matrices. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 3, LWM, LWD ), where */
/*             LWM = LWA + MAX( NC + MAX( NC, MP-1 ), */
/*                              2*MP*MP*MNB - MP*MP + 9*MNB*MNB + */
/*                              MP*MNB + 11*MP + 33*MNB - 11 ); */
/*             LWD = LWB + MAX( 2, LW1, LW2, LW3, LW4, 2*ORD ), */
/*                              if QUTOL >= 0; */
/*             LWD = 0,         if QUTOL <  0; */
/*             LWA = MP*LENDAT + 2*MNB + MP - 1; */
/*             LWB = LENDAT*(MP + 2) + ORD*(ORD + 2) + 1; */
/*             LW1 = 2*LENDAT + 4*HNPTS;  HNPTS = 2048; */
/*             LW2 =   LENDAT + 6*HNPTS;  MN  = MIN( 2*LENDAT, 2*ORD+1 ); */
/*             LW3 = 2*LENDAT*(2*ORD + 1) + MAX( 2*LENDAT, 2*ORD + 1 ) + */
/*                   MAX( MN + 6*ORD + 4, 2*MN + 1 ); */
/*             LW4 = MAX( ORD*ORD + 5*ORD, 6*ORD + 1 + MIN( 1, ORD ) ). */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX( LZM, LZD ), where */
/*             LZM = MAX( MP*MP + NC*MP + NC*NC + 2*NC, */
/*                        6*MP*MP*MNB + 13*MP*MP + 6*MNB + 6*MP - 3 ); */
/*             LZD = MAX( LENDAT*(2*ORD + 3), ORD*ORD + 3*ORD + 1 ), */
/*                              if QUTOL >= 0; */
/*             LZD = 0,         if QUTOL <  0. */

/*     Error indicator */

/*     INFO    (output) INTEGER */
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  if one or more values w in OMEGA are (close to */
/*                    some) poles of the closed-loop system, i.e., the */
/*                    matrix jw*I - A is (numerically) singular; */
/*             =  2:  the block sizes must be positive integers; */
/*             =  3:  the sum of block sizes must be equal to MP; */
/*             =  4:  the size of a real block must be equal to 1; */
/*             =  5:  the block type must be either 1 or 2; */
/*             =  6:  errors in solving linear equations or in matrix */
/*                    inversion; */
/*             =  7:  errors in computing eigenvalues or singular values. */
/*             = 1i:  INFO on exit from SB10YD is i. (1i means 10 + i.) */

/*     METHOD */

/*     I.   First, W(jw) for the given closed-loop system is being */
/*          estimated. */
/*     II.  Now, AB13MD SLICOT subroutine can obtain the D(jw) scaling */
/*          system with respect to NBLOCK and ITYPE, and colaterally, */
/*          mju(jw). */
/*          If QUTOL < 0 then the estimations stop and the routine exits. */
/*     III. Now that we have D(jw), SB10YD subroutine can do block-by- */
/*          block fit. For each block it tries with an increasing order */
/*          of the fit, starting with 1 until the */
/*          (mean quadratic error + max quadratic error)/2 */
/*          between the Dii(jw) and the estimated frequency responce */
/*          of the block becomes less than or equal to the routine */
/*          argument QUTOL, or the order becomes equal to ORD. */
/*     IV.  Arrange the obtained blocks in the AD, BD, CD and DD */
/*          matrices and estimate the total order of D(s), TOTORD. */
/*     V.   Add the system I_f to the system obtained in IV. */

/*     REFERENCES */

/*     [1] Balas, G., Doyle, J., Glover, K., Packard, A. and Smith, R. */
/*         Mu-analysis and Synthesis toolbox - User's Guide, */
/*         The Mathworks Inc., Natick, MA, USA, 1998. */

/*     CONTRIBUTORS */

/*     Asparuh Markovski, Technical University of Sofia, July 2003. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003. */
/*     A. Markovski, V. Sima, October 2003. */

/*     KEYWORDS */

/*     Frequency response, H-infinity optimal control, robust control, */
/*     structured singular value. */

/*     ****************************************************************** */

/*     .. Parameters .. */
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
/*     .. */
/*     .. Intrinsic Functions .. */

/*     Decode and test input parameters. */

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*     Workspace usage 1. */

/*     real */

    /* Parameter adjustments */
    --nblock;
    --itype;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    d_dim1 = *ldd;
    d_offset = 1 + d_dim1;
    d__ -= d_offset;
    --omega;
    ad_dim1 = *ldad;
    ad_offset = 1 + ad_dim1;
    ad -= ad_offset;
    bd_dim1 = *ldbd;
    bd_offset = 1 + bd_dim1;
    bd -= bd_offset;
    cd_dim1 = *ldcd;
    cd_offset = 1 + cd_dim1;
    cd -= cd_offset;
    dd_dim1 = *lddd;
    dd_offset = 1 + dd_dim1;
    dd -= dd_offset;
    --mju;
    --iwork;
    --dwork;
    --zwork;

    /* Function Body */
    iwx = *mp * *lendat + 1;
    iwgjom = iwx + (*mnb << 1) - 1;
    idwrk = iwgjom + *mp;
    ldsize = *ldwork - idwrk + 1;

/*     complex */

    iwb = *mp * *mp + 1;
    icwrk = iwb + *nc * *mp;
    lcsize = *lzwork - icwrk + 1;

    *info = 0;
    if (*nc < 0) {
	*info = -1;
    } else if (*mp < 0) {
	*info = -2;
    } else if (*lendat < 2) {
	*info = -3;
    } else if (*f < 0) {
	*info = -4;
    } else if (*ord > *lendat - 1) {
	*info = -5;
    } else if (*mnb < 1 || *mnb > *mp) {
	*info = -6;
    } else if (*lda < max(1,*nc)) {
	*info = -11;
    } else if (*ldb < max(1,*nc)) {
	*info = -13;
    } else if (*ldc < max(1,*mp)) {
	*info = -15;
    } else if (*ldd < max(1,*mp)) {
	*info = -17;
    } else if (*ldad < 1 || *qutol >= 0. && *ldad < *mp * *ord) {
	*info = -21;
    } else if (*ldbd < 1 || *qutol >= 0. && *ldbd < *mp * *ord) {
	*info = -23;
    } else if (*ldcd < 1 || *qutol >= 0. && *ldcd < *mp + *f) {
	*info = -25;
    } else if (*lddd < 1 || *qutol >= 0. && *lddd < *mp + *f) {
	*info = -27;
    } else {

/*        Compute workspace. */

/* Computing MAX */
	i__1 = *nc, i__2 = (*mnb << 2) - 2, i__1 = max(i__1,i__2);
	ii = max(i__1,*mp);
/* Computing MIN */
	i__1 = *lendat << 1, i__2 = (*ord << 1) + 1;
	mn = min(i__1,i__2);
	lwa = idwrk - 1;
	lwb = *lendat * (*mp + 2) + *ord * (*ord + 2) + 1;
	lw1 = (*lendat << 1) + 8192;
	lw2 = *lendat + 12288;
/* Computing MAX */
	i__1 = *lendat << 1, i__2 = (*ord << 1) + 1;
/* Computing MAX */
	i__3 = mn + *ord * 6 + 4, i__4 = (mn << 1) + 1;
	lw3 = (*lendat << 1) * ((*ord << 1) + 1) + max(i__1,i__2) + max(i__3,
		i__4);
/* Computing MAX */
	i__1 = *ord * *ord + *ord * 5, i__2 = *ord * 6 + 1 + min(1,*ord);
	lw4 = max(i__1,i__2);

/* Computing MAX */
/* Computing MAX */
	i__3 = *nc, i__4 = *mp - 1;
	i__1 = *nc + max(i__3,i__4), i__2 = (*mp << 1) * *mp * *mnb - *mp * *
		mp + *mnb * 9 * *mnb + *mp * *mnb + *mp * 11 + *mnb * 33 - 11;
	dlwmax = lwa + max(i__1,i__2);

/* Computing MAX */
	i__1 = icwrk - 1 + *nc * *nc + (*nc << 1), i__2 = *mp * 6 * *mp * *
		mnb + *mp * 13 * *mp + *mnb * 6 + *mp * 6 - 3;
	clwmax = max(i__1,i__2);

	if (*qutol >= 0.) {
/* Computing MAX */
	    i__1 = ii, i__2 = (*ord << 1) + 1;
	    ii = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
	    i__3 = max(2,lw1), i__3 = max(i__3,lw2), i__3 = max(i__3,lw3), 
		    i__3 = max(i__3,lw4), i__4 = *ord << 1;
	    i__1 = dlwmax, i__2 = lwb + max(i__3,i__4);
	    dlwmax = max(i__1,i__2);
/* Computing MAX */
	    i__1 = clwmax, i__2 = *lendat * ((*ord << 1) + 3), i__1 = max(
		    i__1,i__2), i__2 = *ord * (*ord + 3) + 1;
	    clwmax = max(i__1,i__2);
	}
	if (*liwork < ii) {
	    *info = -30;
	} else if (*ldwork < max(3,dlwmax)) {
	    *info = -32;
	} else if (*lzwork < clwmax) {
	    *info = -34;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("SB10MD", &i__1, (ftnlen)6);
	return 0;
    }

    *ord = max(1,*ord);
    *totord = 0;

/*     Quick return if possible. */

    if (*nc == 0 || *mp == 0) {
	dwork[1] = 3.;
	dwork[2] = 0.;
	dwork[3] = 1.;
	return 0;
    }

    toler = sqrt(dlamch_("Epsilon", (ftnlen)7));

    *(unsigned char *)baleig = 'C';
    rcond = 1.;
    maxcwr = clwmax;

/*     @@@ 1. Estimate W(jw) for the closed-loop system, @@@ */
/*     @@@      D(jw) and mju(jw) for each frequency.    @@@ */

    i__1 = *lendat;
    for (w = 1; w <= i__1; ++w) {
	i__2 = w;
	z__1.r = 0., z__1.i = omega[i__2];
	freq.r = z__1.r, freq.i = z__1.i;
	if (w == 1) {
	    *(unsigned char *)inita = 'G';
	} else {
	    *(unsigned char *)inita = 'H';
	}

/*        Compute C*inv(jw*I-A)*B. */
/*        Integer workspace: need   NC. */
/*        Real workspace:    need   LWA + NC + MAX(NC,MP-1); */
/*                           prefer larger, */
/*                           where  LWA = MP*LENDAT + 2*MNB + MP - 1. */
/*        Complex workspace: need   MP*MP + NC*MP + NC*NC + 2*NC. */

	tb05ad_(baleig, inita, nc, mp, mp, &freq, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &rcnd, &zwork[1], mp, &
		dwork[1], &dwork[1], &zwork[iwb], nc, &iwork[1], &dwork[idwrk]
		, &ldsize, &zwork[icwrk], &lcsize, &info2, (ftnlen)1, (ftnlen)
		1);

	if (info2 > 0) {
	    *info = 1;
	    return 0;
	}

	rcond = min(rcond,rcnd);
	if (w == 1) {
	    maxwrk = (integer) (dwork[idwrk] + idwrk - 1);
	}
	ic = 0;

/*        D + C*inv(jw*I-A)*B */

	i__2 = *mp;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		++ic;
		i__4 = ic;
		i__5 = ic;
		i__6 = i__ + k * d_dim1;
		z__2.r = d__[i__6], z__2.i = 0.;
		z__1.r = zwork[i__5].r + z__2.r, z__1.i = zwork[i__5].i + 
			z__2.i;
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
/* L10: */
	    }
/* L20: */
	}

/*        Estimate D(jw) and mju(jw). */
/*        Integer workspace: need   MAX(4*MNB-2,MP). */
/*        Real workspace:    need   LWA + 2*MP*MP*MNB - MP*MP + 9*MNB*MNB */
/*                                  + MP*MNB + 11*MP + 33*MNB - 11; */
/*                           prefer larger. */
/*        Complex workspace: need   6*MP*MP*MNB + 13*MP*MP + 6*MNB + */
/*                                  6*MP - 3. */

	i__2 = *lzwork - iwb + 1;
	ab13md_("N", mp, &zwork[1], mp, mnb, &nblock[1], &itype[1], &dwork[
		iwx], &mju[w], &dwork[(w - 1) * *mp + 1], &dwork[iwgjom], &
		iwork[1], &dwork[idwrk], &ldsize, &zwork[iwb], &i__2, &info2, 
		(ftnlen)1);

	if (info2 != 0) {
	    *info = info2 + 1;
	    return 0;
	}

	if (w == 1) {
/* Computing MAX */
	    i__2 = maxwrk, i__3 = (integer) dwork[idwrk] + idwrk - 1;
	    maxwrk = max(i__2,i__3);
/* Computing MAX */
	    i__4 = iwb;
	    i__2 = maxcwr, i__3 = (integer) zwork[i__4].r + iwb - 1;
	    maxcwr = max(i__2,i__3);
	}

/*        Normalize D(jw) through it's last entry. */

	if (dwork[w * *mp] != 0.) {
	    d__1 = 1. / dwork[w * *mp];
	    dscal_(mp, &d__1, &dwork[(w - 1) * *mp + 1], &c__1);
	}

/* L30: */
    }

/*     Quick return if needed. */

    if (*qutol < 0.) {
	dwork[1] = (doublereal) maxwrk;
	dwork[2] = (doublereal) maxcwr;
	dwork[3] = rcond;
	return 0;
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*     Workspace usage 2. */

/*     real */

    iwrfrd = iwx;
    iwifrd = iwrfrd + *lendat;
    iwad = iwifrd + *lendat;
    iwbd = iwad + *ord * *ord;
    iwcd = iwbd + *ord;
    iwdd = iwcd + *ord;
    idwrk = iwdd + 1;
    ldsize = *ldwork - idwrk + 1;

/*     complex */

    icwrk = *ord + 2;
    lcsize = *lzwork - icwrk + 1;
    *(unsigned char *)inita = 'H';

/*     Use default tolerance for SB10YD. */

    tol = -1.;

/*     @@@ 2. Clear imag parts of D(jw) for SB10YD. @@@ */

    i__1 = *lendat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dwork[iwifrd + i__ - 1] = 0.;
/* L40: */
    }

/*     @@@ 3. Clear AD, BD, CD and initialize DD with I_(mp+f). @@@ */

    i__1 = *mp * *ord;
    i__2 = *mp * *ord;
    dlaset_("Full", &i__1, &i__2, &c_b15, &c_b15, &ad[ad_offset], ldad, (
	    ftnlen)4);
    i__1 = *mp * *ord;
    i__2 = *mp + *f;
    dlaset_("Full", &i__1, &i__2, &c_b15, &c_b15, &bd[bd_offset], ldbd, (
	    ftnlen)4);
    i__1 = *mp + *f;
    i__2 = *mp * *ord;
    dlaset_("Full", &i__1, &i__2, &c_b15, &c_b15, &cd[cd_offset], ldcd, (
	    ftnlen)4);
    i__1 = *mp + *f;
    i__2 = *mp + *f;
    dlaset_("Full", &i__1, &i__2, &c_b15, &c_b25, &dd[dd_offset], lddd, (
	    ftnlen)4);

/*     @@@ 4. Block by block frequency identification. @@@ */

    i__1 = *mp;
    for (ii = 1; ii <= i__1; ++ii) {

	dcopy_(lendat, &dwork[ii], mp, &dwork[iwrfrd], &c__1);

/*        Increase CORD from 1 to ORD for every block, if needed. */

	cord = 1;

L50:
	lord = cord;

/*           Now, LORD is the desired order. */
/*           Integer workspace: need   2*N+1, where N = LORD. */
/*           Real workspace:    need   LWB + MAX( 2, LW1, LW2, LW3, LW4), */
/*                                     where */
/*                                     LWB = LENDAT*(MP+2) + */
/*                                           ORD*(ORD+2) + 1, */
/*                                     HNPTS = 2048, and */
/*                                     LW1 = 2*LENDAT + 4*HNPTS; */
/*                                     LW2 =   LENDAT + 6*HNPTS; */
/*                                     MN  = min( 2*LENDAT, 2*N+1 ) */
/*                                     LW3 = 2*LENDAT*(2*N+1) + */
/*                                           max( 2*LENDAT, 2*N+1 ) + */
/*                                           max( MN + 6*N + 4, 2*MN+1 ); */
/*                                     LW4 = max( N*N + 5*N, */
/*                                                6*N + 1 + min( 1,N ) ); */
/*                              prefer larger. */
/*           Complex workspace: need   LENDAT*(2*N+3). */

	sb10yd_(&c__0, &c__1, lendat, &dwork[iwrfrd], &dwork[iwifrd], &omega[
		1], &lord, &dwork[iwad], ord, &dwork[iwbd], &dwork[iwcd], &
		dwork[iwdd], &tol, &iwork[1], &dwork[idwrk], &ldsize, &zwork[
		1], lzwork, &info2);

/*           At this point, LORD is the actual order reached by SB10YD, */
/*           0 <= LORD <= CORD. */
/*           [ADi,BDi; CDi,DDi] is a minimal realization with ADi in */
/*           upper Hessenberg form. */
/*           The leading LORD-by-LORD part of ORD-by-ORD DWORK(IWAD) */
/*           contains ADi, the leading LORD-by-1 part of ORD-by-1 */
/*           DWORK(IWBD) contains BDi, the leading 1-by-LORD part of */
/*           1-by-ORD DWORK(IWCD) contains CDi, DWORK(IWDD) contains DDi. */

	if (info2 != 0) {
	    *info = info2 + 10;
	    return 0;
	}

/*          Compare the original D(jw) with the fitted one. */

	meqe = 0.;
	maqe = 0.;

	i__2 = *lendat;
	for (w = 1; w <= i__2; ++w) {
	    i__3 = w;
	    z__1.r = 0., z__1.i = omega[i__3];
	    freq.r = z__1.r, freq.i = z__1.i;

/*              Compute CD*inv(jw*I-AD)*BD. */
/*              Integer workspace: need   LORD. */
/*              Real workspace:    need   LWB + 2*LORD; */
/*                                 prefer larger. */
/*              Complex workspace: need   1 + ORD + LORD*LORD + 2*LORD. */

	    tb05ad_(baleig, inita, &lord, &c__1, &c__1, &freq, &dwork[iwad], 
		    ord, &dwork[iwbd], ord, &dwork[iwcd], &c__1, &rcnd, &
		    zwork[1], &c__1, &dwork[idwrk], &dwork[idwrk], &zwork[2], 
		    ord, &iwork[1], &dwork[idwrk], &ldsize, &zwork[icwrk], &
		    lcsize, &info2, (ftnlen)1, (ftnlen)1);

	    if (info2 > 0) {
		*info = 1;
		return 0;
	    }

	    rcond = min(rcond,rcnd);
	    if (w == 1) {
/* Computing MAX */
		i__3 = maxwrk, i__4 = (integer) dwork[idwrk] + idwrk - 1;
		maxwrk = max(i__3,i__4);
	    }

/*              DD + CD*inv(jw*I-AD)*BD */

	    i__3 = iwdd;
	    z__2.r = dwork[i__3], z__2.i = 0.;
	    z__1.r = zwork[1].r + z__2.r, z__1.i = zwork[1].i + z__2.i;
	    zwork[1].r = z__1.r, zwork[1].i = z__1.i;

	    mod1 = (d__1 = dwork[iwrfrd + w - 1], abs(d__1));
	    mod2 = z_abs(&zwork[1]);
	    rqe = (d__1 = (mod1 - mod2) / (mod1 + toler), abs(d__1));
	    meqe += rqe;
	    maqe = max(maqe,rqe);

/* L60: */
	}

	meqe /= *lendat;

	if ((meqe + maqe) / 2. <= *qutol || cord == *ord) {
	    goto L70;
	}

	++cord;
	goto L50;

L70:
	*totord += lord;

/*        Copy ad(ii), bd(ii) and cd(ii) to AD, BD and CD, respectively. */

	dlacpy_("Full", &lord, &lord, &dwork[iwad], ord, &ad[*totord - lord + 
		1 + (*totord - lord + 1) * ad_dim1], ldad, (ftnlen)4);
	dcopy_(&lord, &dwork[iwbd], &c__1, &bd[*totord - lord + 1 + ii * 
		bd_dim1], &c__1);
	dcopy_(&lord, &dwork[iwcd], &c__1, &cd[ii + (*totord - lord + 1) * 
		cd_dim1], ldcd);

/*        Copy dd(ii) to DD. */

	dd[ii + ii * dd_dim1] = dwork[iwdd];

/* L80: */
    }

    dwork[1] = (doublereal) maxwrk;
    dwork[2] = (doublereal) maxcwr;
    dwork[3] = rcond;
    return 0;

/* *** Last line of SB10MD *** */
} /* sb10md_ */

