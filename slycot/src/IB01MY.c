/* IB01MY.f -- translated by f2c (version 20100827).
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
static doublereal c_b18 = 1.;
static integer c__0 = 0;
static doublereal c_b111 = 0.;

/* Subroutine */ int ib01my_(char *meth, char *batch, char *conct, integer *
	nobr, integer *m, integer *l, integer *nsmp, doublereal *u, integer *
	ldu, doublereal *y, integer *ldy, doublereal *r__, integer *ldr, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, 
	integer *info, ftnlen meth_len, ftnlen batch_len, ftnlen conct_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, jd;
    static doublereal cs;
    static integer nr;
    static doublereal sn;
    static integer ns, icj, ing, ipg, jds;
    static doublereal dum[1];
    static integer nrg;
    static doublereal upd, tau;
    static integer nsm, ipy;
    static doublereal beta;
    static integer ingc, ipgc, icol, imax, ingp, ierr, itau, lnrg, mnrg, irev;
    static logical last, n4sid;
    static integer nobr2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), ma02fd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), mb04id_(integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), mb04od_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *),
	     dlarf_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen), 
	    dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nobr21, iconn, lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp, first;
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static logical onebch;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static logical connec;
    static integer icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer llnobr, mmnobr;
    static logical interm;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer ldrwrk, minwrk, maxwrk, nsmpsm;


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

/*     To construct an upper triangular factor  R  of the concatenated */
/*     block Hankel matrices using input-output data, via a fast QR */
/*     algorithm based on displacement rank.  The input-output data can, */
/*     optionally, be processed sequentially. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     BATCH   CHARACTER*1 */
/*             Specifies whether or not sequential data processing is to */
/*             be used, and, for sequential processing, whether or not */
/*             the current data block is the first block, an intermediate */
/*             block, or the last block, as follows: */
/*             = 'F':  the first block in sequential data processing; */
/*             = 'I':  an intermediate block in sequential data */
/*                     processing; */
/*             = 'L':  the last block in sequential data processing; */
/*             = 'O':  one block only (non-sequential data processing). */
/*             NOTE that when  100  cycles of sequential data processing */
/*                  are completed for  BATCH = 'I',  a warning is */
/*                  issued, to prevent for an infinite loop. */

/*     CONCT   CHARACTER*1 */
/*             Specifies whether or not the successive data blocks in */
/*             sequential data processing belong to a single experiment, */
/*             as follows: */
/*             = 'C':  the current data block is a continuation of the */
/*                     previous data block and/or it will be continued */
/*                     by the next data block; */
/*             = 'N':  there is no connection between the current data */
/*                     block and the previous and/or the next ones. */
/*             This parameter is not used if BATCH = 'O'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             block Hankel matrices to be processed.  NOBR > 0. */
/*             (In the MOESP theory,  NOBR  should be larger than  n, the */
/*             estimated dimension of state vector.) */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */
/*             When M = 0, no system inputs are processed. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples,  t). (When sequential data processing is used, */
/*             NSMP  is the number of samples of the current data */
/*             block.) */
/*             NSMP >= 2*(M+L+1)*NOBR - 1,  for non-sequential */
/*                                          processing; */
/*             NSMP >= 2*NOBR,  for sequential processing. */
/*             The total number of samples when calling the routine with */
/*             BATCH = 'L'  should be at least  2*(M+L+1)*NOBR - 1. */
/*             The  NSMP  argument may vary from a cycle to another in */
/*             sequential data processing, but  NOBR, M,  and  L  should */
/*             be kept constant. For efficiency, it is advisable to use */
/*             NSMP  as large as possible. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
/*             If M = 0, this array is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= NSMP, if M > 0; */
/*             LDU >= 1,    if M = 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             t-by-l output-data sequence matrix  Y, */
/*             Y = [y_1 y_2 ... y_l].  Column  j  of  Y  contains the */
/*             NSMP  values of the j-th output component for consecutive */
/*             time increments. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= NSMP. */

/*     R       (output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             If INFO = 0 and BATCH = 'L' or 'O', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the upper triangular factor R from the */
/*             QR factorization of the concatenated block Hankel */
/*             matrices. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first (M+L)*2*NOBR*(M+L+c) elements of  DWORK  should */
/*             be preserved during successive calls of the routine */
/*             with  BATCH = 'F'  or  'I',  till the final call with */
/*             BATCH = 'L',  where */
/*             c = 1,  if the successive data blocks do not belong to a */
/*                     single experiment  (CONCT = 'N'); */
/*             c = 2,  if the successive data blocks belong to a single */
/*                     experiment  (CONCT = 'C'). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+3), */
/*                              if BATCH <> 'O' and CONCT = 'C'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+1), */
/*                              if BATCH = 'F' or 'I' and CONCT = 'N'; */
/*             LDWORK >= (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR, */
/*                              if BATCH = 'L' and CONCT = 'N', */
/*                              or BATCH = 'O'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the number of 100 cycles in sequential data */
/*                   processing has been exhausted without signaling */
/*                   that the last block of data was get. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the fast QR factorization algorithm failed. The */
/*                   matrix H'*H is not (numerically) positive definite. */

/*     METHOD */

/*     Consider the  t x 2(m+l)s  matrix H of concatenated block Hankel */
/*     matrices */

/*          H = [ Uf'         Up'      Y'      ],  for METH = 'M', */
/*                  s+1,2s,t    1,s,t   1,2s,t */

/*          H = [ U'       Y'      ],              for METH = 'N', */
/*                 1,2s,t   1,2s,t */

/*     where  Up     , Uf        , U      , and  Y        are block */
/*              1,s,t    s+1,2s,t   1,2s,t        1,2s,t */
/*     Hankel matrices defined in terms of the input and output data [3]. */
/*     The fast QR algorithm uses a factorization of H'*H which exploits */
/*     the block-Hankel structure, via a displacement rank technique [5]. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error */
/*         state-space model identification class of algorithms. */
/*         Int. J. Control, 56, pp. 1187-1210, 1992. */

/*     [2] Verhaegen M. */
/*         Subspace Model Identification. Part 3: Analysis of the */
/*         ordinary output-error state-space model identification */
/*         algorithm. */
/*         Int. J. Control, 58, pp. 555-586, 1993. */

/*     [3] Verhaegen M. */
/*         Identification of the deterministic part of MIMO state space */
/*         models given in innovations form from input-output data. */
/*         Automatica, Vol.30, No.1, pp.61-74, 1994. */

/*     [4] Van Overschee, P., and De Moor, B. */
/*         N4SID: Subspace Algorithms for the Identification of */
/*         Combined Deterministic-Stochastic Systems. */
/*         Automatica, Vol.30, No.1, pp. 75-93, 1994. */

/*     [5] Kressner, D., Mastronardi, N., Sima, V., Van Dooren, P., and */
/*         Van Huffel, S. */
/*         A Fast Algorithm for Subspace State-space System */
/*         Identification via Exploitation of the Displacement Structure. */
/*         J. Comput. Appl. Math., Vol.132, No.1, pp. 71-81, 2001. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is reliable and efficient. Numerical */
/*     difficulties are possible when the matrix H'*H is nearly rank */
/*     defficient. The method cannot be used if the matrix H'*H is not */
/*     numerically positive definite. */
/*                                     2           3 2 */
/*     The algorithm requires 0(2t(m+l) s)+0(4(m+l) s ) floating point */
/*     operations. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Universiteit Leuven, June 2000. */
/*     Partly based on Matlab codes developed by N. Mastronardi, */
/*     Katholieke Universiteit Leuven, February 2000. */

/*     REVISIONS */

/*     V. Sima, July 2000, August 2000, Feb. 2004, May 2009. */

/*     KEYWORDS */

/*     Displacement rank, Hankel matrix, Householder transformation, */
/*     identification methods, multivariable systems. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save Statement .. */
/*        ICYCLE  is used to count the cycles for  BATCH = 'I'. */
/*        MAXWRK  is used to store the optimal workspace. */
/*        NSMPSM  is used to sum up the  NSMP  values for  BATCH <> 'O'. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --iwork;
    --dwork;

    /* Function Body */
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
    onebch = lsame_(batch, "O", (ftnlen)1, (ftnlen)1);
    first = lsame_(batch, "F", (ftnlen)1, (ftnlen)1) || onebch;
    interm = lsame_(batch, "I", (ftnlen)1, (ftnlen)1);
    last = lsame_(batch, "L", (ftnlen)1, (ftnlen)1) || onebch;
    if (! onebch) {
	connec = lsame_(conct, "C", (ftnlen)1, (ftnlen)1);
    } else {
	connec = FALSE_;
    }
    mnobr = *m * *nobr;
    lnobr = *l * *nobr;
    mmnobr = mnobr + mnobr;
    llnobr = lnobr + lnobr;
    nobr2 = *nobr << 1;
    nobr21 = nobr2 - 1;
    *iwarn = 0;
    *info = 0;
    if (first) {
	icycle = 1;
	maxwrk = 1;
	nsmpsm = 0;
    }
    nsmpsm += *nsmp;
    nr = mmnobr + llnobr;

/*     Check the scalar input parameters. */

    if (! (moesp || n4sid)) {
	*info = -1;
    } else if (! (first || interm || last)) {
	*info = -2;
    } else if (! onebch) {
	if (! (connec || lsame_(conct, "N", (ftnlen)1, (ftnlen)1))) {
	    *info = -3;
	}
    }
    if (*info == 0) {
	if (*nobr <= 0) {
	    *info = -4;
	} else if (*m < 0) {
	    *info = -5;
	} else if (*l <= 0) {
	    *info = -6;
	} else if (*nsmp < nobr2 || last && nsmpsm < nr + nobr21) {
	    *info = -7;
	} else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
	    *info = -9;
	} else if (*ldy < *nsmp) {
	    *info = -11;
	} else if (*ldr < nr) {
	    *info = -13;
	} else {

/*           Compute workspace. */
/*           NRG is the number of positive (or negative) generators. */

	    nrg = *m + *l + 1;
	    if (! onebch && connec) {
		minwrk = nr * (nrg + 2);
	    } else if (first || interm) {
		minwrk = nr * nrg;
	    } else {
		minwrk = (nr << 1) * nrg + nr;
	    }
	    maxwrk = max(minwrk,maxwrk);

	    if (*ldwork < minwrk) {
		*info = -16;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	nsmpsm = 0;
	if (*info == -16) {
	    dwork[1] = (doublereal) minwrk;
	}
	i__1 = -(*info);
	xerbla_("IB01MY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compute the  R  factor from a fast QR factorization of the */
/*     matrix  H,  a concatenation of two block Hankel matrices. */
/*     Specifically, a displacement rank technique is applied to */
/*     the block Toeplitz matrix,  G = (P*H)'*(P*H),  where  P  is a */
/*     2-by-2 block diagonal matrix, having as diagonal blocks identity */
/*     matrices with columns taken in the reverse order. */
/*     The technique builds and processes the generators of  G.  The */
/*     matrices  G  and  G1 = H'*H  have the same  R  factor. */

/*     Set the parameters for constructing the correlations of the */
/*     current block. */
/*     NSM is the number of processed samples in U and Y, t - 2s. */
/*     IPG and ING are pointers to the "positive" and "negative" */
/*     generators, stored row-wise in the workspace. All "positive" */
/*     generators are stored before any "negative" generators. */
/*     If BATCH <> 'O' and CONCT = 'C', the "connection" elements of */
/*     two successive batches are stored in the same workspace as the */
/*     "negative" generators (which will be computed later on). */
/*     IPY is a pointer to the Y part of the "positive" generators. */
/*     LDRWRK is used as a leading dimension for the workspace part used */
/*     to store the "connection" elements. */

    ns = *nsmp - nobr21;
    nsm = ns - 1;
    mnrg = *m * nrg;
    lnrg = *l * nrg;

    ldrwrk = nobr2 << 1;
    if (first) {
	upd = 0.;
    } else {
	upd = 1.;
    }
    dum[0] = 0.;

    ipg = 1;
    ipy = ipg + *m;
    ing = ipg + nrg * nr;
    iconn = ing;

    if (! first && connec) {

/*        Restore the saved (M+L)*2*NOBR "connection" elements of */
/*        U  and  Y  into their appropriate position in sequential */
/*        processing. The process is performed column-wise, in */
/*        reverse order, first for  Y  and then for  U. */
/*        ICONN is a pointer to the first saved "connection" element. */
/*        Workspace: need   (M+L)*2*NOBR*(M+L+3). */

	irev = iconn + nr;
	icol = iconn + (nr << 1);

	i__1 = *m + *l;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    irev -= nobr2;
	    icol -= ldrwrk;
	    dcopy_(&nobr2, &dwork[irev], &c__1, &dwork[icol], &c__1);
/* L10: */
	}

	if (*m > 0) {
	    dlacpy_("Full", &nobr2, m, &u[u_offset], ldu, &dwork[iconn + 
		    nobr2], &ldrwrk, (ftnlen)4);
	}
	dlacpy_("Full", &nobr2, l, &y[y_offset], ldy, &dwork[iconn + ldrwrk * 
		*m + nobr2], &ldrwrk, (ftnlen)4);
    }

    if (*m > 0) {

/*        Let  Guu(i,j) = Guu0(i,j) + u_i*u_j' + u_(i+1)*u_(j+1)' + */
/*                              ... + u_(i+NSM-1)*u_(j+NSM-1)', */
/*        where  u_i'  is the i-th row of  U,  j = 1 : 2s,  i = 1 : j, */
/*        NSM = NSMP - 2s,  and  Guu0(i,j)  is a zero matrix for */
/*        BATCH = 'O' or 'F', and it is the matrix Guu(i,j) computed */
/*        till the current block for BATCH = 'I' or 'L'. The matrix */
/*        Guu(i,j)  is  m-by-m,  and  Guu(j,j)  is symmetric. The */
/*        submatrices of the first block-row, Guu(1,j), are needed only. */

/*        Compute/update  Guu(1,1). */

	if (! first && connec) {
	    dsyrk_("Upper", "Transpose", m, &nobr2, &c_b18, &dwork[iconn], &
		    ldrwrk, &upd, &dwork[ipg], &nrg, (ftnlen)5, (ftnlen)9);
	}
	dsyrk_("Upper", "Transpose", m, &nsm, &c_b18, &u[u_offset], ldu, &upd,
		 &dwork[ipg], &nrg, (ftnlen)5, (ftnlen)9);
	ma02ed_("Upper", m, &dwork[ipg], &nrg, (ftnlen)5);

	jd = 1;

	if (first || ! connec) {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *m;

/*              Compute/update  Guu(1,j). */

		dgemm_("Transpose", "NoTranspose", m, m, &nsm, &c_b18, &u[
			u_offset], ldu, &u[j + u_dim1], ldu, &upd, &dwork[ipg 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L20: */
	    }

	} else {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *m;

/*              Compute/update  Guu(1,j)  for sequential processing */
/*              with connected blocks. */

		dgemm_("Transpose", "NoTranspose", m, m, &nobr2, &c_b18, &
			dwork[iconn], &ldrwrk, &dwork[iconn + j - 1], &ldrwrk,
			 &upd, &dwork[ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, 
			(ftnlen)11);
		dgemm_("Transpose", "NoTranspose", m, m, &nsm, &c_b18, &u[
			u_offset], ldu, &u[j + u_dim1], ldu, &c_b18, &dwork[
			ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L30: */
	    }

	}

/*        Let  Guy(i,j) = Guy0(i,j) + u_i*y_j' + u_(i+1)*y_(j+1)' + */
/*                              ... + u_(i+NSM-1)*y_(j+NSM-1)', */
/*        where  u_i'  is the i-th row of  U,  y_j'  is the j-th row */
/*        of  Y,  j = 1 : 2s,  i = 1 : 2s,  NSM = NSMP - 2s,  and */
/*        Guy0(i,j)  is a zero matrix for  BATCH = 'O' or 'F', and it */
/*        is the matrix Guy(i,j) computed till the current block for */
/*        BATCH = 'I' or 'L'.  Guy(i,j) is m-by-L. The submatrices */
/*        of the first block-row, Guy(1,j), as well as the transposes */
/*        of the submatrices of the first block-column, i.e., Gyu(1,j), */
/*        are needed only. */

	jd = mmnobr + 1;

	if (first || ! connec) {

	    i__1 = nobr2;
	    for (j = 1; j <= i__1; ++j) {

/*              Compute/update  Guy(1,j). */

		dgemm_("Transpose", "NoTranspose", m, l, &nsm, &c_b18, &u[
			u_offset], ldu, &y[j + y_dim1], ldy, &upd, &dwork[ipg 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
		jd += *l;
/* L40: */
	    }

	} else {

	    i__1 = nobr2;
	    for (j = 1; j <= i__1; ++j) {

/*              Compute/update  Guy(1,j)  for sequential processing */
/*              with connected blocks. */

		dgemm_("Transpose", "NoTranspose", m, l, &nobr2, &c_b18, &
			dwork[iconn], &ldrwrk, &dwork[iconn + ldrwrk * *m + j 
			- 1], &ldrwrk, &upd, &dwork[ipg + (jd - 1) * nrg], &
			nrg, (ftnlen)9, (ftnlen)11);
		dgemm_("Transpose", "NoTranspose", m, l, &nsm, &c_b18, &u[
			u_offset], ldu, &y[j + y_dim1], ldy, &c_b18, &dwork[
			ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
		jd += *l;
/* L50: */
	    }

	}

/*        Now, the first M "positive" generators have been built. */
/*        Transpose  Guy(1,1)  in the first block of the  Y  part of the */
/*        "positive" generators. */

	i__1 = *l;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(m, &dwork[ipg + (mmnobr + j - 1) * nrg], &c__1, &dwork[ipy 
		    + j - 1], &nrg);
/* L60: */
	}

	jd = 1;

	if (first || ! connec) {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *m;

/*              Compute/update  Gyu(1,j). */

		dgemm_("Transpose", "NoTranspose", l, m, &nsm, &c_b18, &y[
			y_offset], ldy, &u[j + u_dim1], ldu, &upd, &dwork[ipy 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L70: */
	    }

	} else {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *m;

/*              Compute/update  Gyu(1,j)  for sequential processing */
/*              with connected blocks. */

		dgemm_("Transpose", "NoTranspose", l, m, &nobr2, &c_b18, &
			dwork[iconn + ldrwrk * *m], &ldrwrk, &dwork[iconn + j 
			- 1], &ldrwrk, &upd, &dwork[ipy + (jd - 1) * nrg], &
			nrg, (ftnlen)9, (ftnlen)11);
		dgemm_("Transpose", "NoTranspose", l, m, &nsm, &c_b18, &y[
			y_offset], ldy, &u[j + u_dim1], ldu, &c_b18, &dwork[
			ipy + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L80: */
	    }

	}

    }

/*     Let  Gyy(i,j) = Gyy0(i,j) + y_i*y_i' + y_(i+1)*y_(i+1)' + ... + */
/*                                 y_(i+NSM-1)*y_(i+NSM-1)', */
/*     where  y_i'  is the i-th row of  Y,  j = 1 : 2s,  i = 1 : j, */
/*     NSM = NSMP - 2s,  and  Gyy0(i,j)  is a zero matrix for */
/*     BATCH = 'O' or 'F', and it is the matrix Gyy(i,j) computed till */
/*     the current block for BATCH = 'I' or 'L'.  Gyy(i,j) is L-by-L, */
/*     and  Gyy(j,j)  is symmetric. The submatrices of the first */
/*     block-row, Gyy(1,j), are needed only. */

    jd = mmnobr + 1;

/*     Compute/update  Gyy(1,1). */

    if (! first && connec) {
	dsyrk_("Upper", "Transpose", l, &nobr2, &c_b18, &dwork[iconn + ldrwrk 
		* *m], &ldrwrk, &upd, &dwork[ipy + mmnobr * nrg], &nrg, (
		ftnlen)5, (ftnlen)9);
    }
    dsyrk_("Upper", "Transpose", l, &nsm, &c_b18, &y[y_offset], ldy, &upd, &
	    dwork[ipy + mmnobr * nrg], &nrg, (ftnlen)5, (ftnlen)9);
    ma02ed_("Upper", l, &dwork[ipy + mmnobr * nrg], &nrg, (ftnlen)5);

    if (first || ! connec) {

	i__1 = nobr2;
	for (j = 2; j <= i__1; ++j) {
	    jd += *l;

/*           Compute/update  Gyy(1,j). */

	    dgemm_("Transpose", "NoTranspose", l, l, &nsm, &c_b18, &y[
		    y_offset], ldy, &y[j + y_dim1], ldy, &upd, &dwork[ipy + (
		    jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L90: */
	}

    } else {

	i__1 = nobr2;
	for (j = 2; j <= i__1; ++j) {
	    jd += *l;

/*           Compute/update  Gyy(1,j)  for sequential processing with */
/*           connected blocks. */

	    dgemm_("Transpose", "NoTranspose", l, l, &nobr2, &c_b18, &dwork[
		    iconn + ldrwrk * *m], &ldrwrk, &dwork[iconn + ldrwrk * *m 
		    + j - 1], &ldrwrk, &upd, &dwork[ipy + (jd - 1) * nrg], &
		    nrg, (ftnlen)9, (ftnlen)11);
	    dgemm_("Transpose", "NoTranspose", l, l, &nsm, &c_b18, &y[
		    y_offset], ldy, &y[j + y_dim1], ldy, &c_b18, &dwork[ipy + 
		    (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
/* L100: */
	}

    }

    if (! last) {
	if (first) {

/*           For sequential processing, save the first 2*NOBR-1 rows of */
/*           the first block of  U  and  Y  in the appropriate */
/*           (M+L)*(2*NOBR-1) locations of  DWORK  starting at (1+M)*NRG. */
/*           These will be used to construct the last negative generator. */

	    jd = nrg;
	    if (*m > 0) {
		dcopy_(m, dum, &c__0, &dwork[jd], &nrg);

		i__1 = nobr21;
		for (j = 1; j <= i__1; ++j) {
		    jd += mnrg;
		    dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
/* L110: */
		}

		jd += mnrg;
	    }
	    dcopy_(l, dum, &c__0, &dwork[jd], &nrg);

	    i__1 = nobr21;
	    for (j = 1; j <= i__1; ++j) {
		jd += lnrg;
		dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
/* L120: */
	    }

	}

	if (connec) {

/*           For sequential processing with connected data blocks, */
/*           save the remaining ("connection") elements of  U  and  Y */
/*           in (M+L)*2*NOBR locations of  DWORK  starting at ICONN. */

	    if (*m > 0) {
		dlacpy_("Full", &nobr2, m, &u[ns + u_dim1], ldu, &dwork[iconn]
			, &nobr2, (ftnlen)4);
	    }
	    dlacpy_("Full", &nobr2, l, &y[ns + y_dim1], ldy, &dwork[iconn + 
		    mmnobr], &nobr2, (ftnlen)4);
	}

/*        Return to get new data. */

	++icycle;
	if (icycle > 100) {
	    *iwarn = 1;
	}
	return 0;
    }

    if (last) {

/*        Try to compute the R factor. */

/*        Scale the first M+L positive generators and set the first */
/*        M+L negative generators. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+M+L. */

	jwork = (nrg << 1) * nr + 1;
	i__1 = nrg + 1;
	dcopy_(m, &dwork[ipg], &i__1, &dwork[jwork], &c__1);
	i__1 = nrg + 1;
	dcopy_(l, &dwork[ipy + mmnobr * nrg], &i__1, &dwork[jwork + *m], &
		c__1);

	i__1 = *m + *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *m + *l;
	    iwork[i__] = idamax_(&i__2, &dwork[jwork], &c__1);
	    dwork[jwork + iwork[i__] - 1] = 0.;
/* L130: */
	}

	i__1 = *m + *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    imax = iwork[i__];
	    if (imax <= *m) {
		icol = imax;
	    } else {
		icol = mmnobr - *m + imax;
	    }
	    beta = sqrt((d__1 = dwork[ipg + imax - 1 + (icol - 1) * nrg], abs(
		    d__1)));
	    if (beta == 0.) {

/*              Error exit. */

		*info = 1;
		return 0;
	    }
	    d__1 = 1. / beta;
	    dscal_(&nr, &d__1, &dwork[ipg + imax - 1], &nrg);
	    dcopy_(&nr, &dwork[ipg + imax - 1], &nrg, &dwork[ing + imax - 1], 
		    &nrg);
	    dwork[ipg + imax - 1 + (icol - 1) * nrg] = beta;
	    dwork[ing + imax - 1 + (icol - 1) * nrg] = 0.;

	    i__2 = *m + *l;
	    for (j = i__ + 1; j <= i__2; ++j) {
		dwork[ipg + iwork[j] - 1 + (icol - 1) * nrg] = 0.;
/* L140: */
	    }

/* L150: */
	}

/*        Compute the last two generators. */

	if (! first) {

/*           For sequential processing, move the stored last negative */
/*           generator. */

	    dcopy_(&nr, &dwork[nrg], &nrg, &dwork[ing + nrg - 1], &nrg);
	}

	jd = nrg;
	if (*m > 0) {

	    i__1 = *nsmp;
	    for (j = ns; j <= i__1; ++j) {
		dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
		jd += mnrg;
/* L160: */
	    }

	}

	i__1 = *nsmp;
	for (j = ns; j <= i__1; ++j) {
	    dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
	    jd += lnrg;
/* L170: */
	}

	if (first) {
	    if (*m > 0) {
		dcopy_(m, dum, &c__0, &dwork[jd], &nrg);

		i__1 = nobr21;
		for (j = 1; j <= i__1; ++j) {
		    jd += mnrg;
		    dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
/* L180: */
		}

		jd += mnrg;
	    }
	    dcopy_(l, dum, &c__0, &dwork[jd], &nrg);

	    i__1 = nobr21;
	    for (j = 1; j <= i__1; ++j) {
		jd += lnrg;
		dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
/* L190: */
	    }

	}

	itau = jwork;
	ipgc = ipg + mmnobr * nrg;

	if (*m > 0) {

/*           Process the input part of the generators. */

	    jwork = itau + *m;

/*           Reduce the first M columns of the matrix G1 of positive */
/*           generators to an upper triangular form. */
/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*M; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+M+M*NB. */

	    ingc = ing;
	    i__1 = *ldwork - jwork + 1;
	    dgeqrf_(&nrg, m, &dwork[ipg], &nrg, &dwork[itau], &dwork[jwork], &
		    i__1, &ierr);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR; */
/*                      prefer (M+L)*4*NOBR*(M+L+1)+M+ */
/*                                                 ((M+L)*2*NOBR-M)*NB. */

	    i__1 = nr - *m;
	    i__2 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &nrg, &i__1, m, &dwork[ipg], &nrg, &
		    dwork[itau], &dwork[ipg + mnrg], &nrg, &dwork[jwork], &
		    i__2, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Annihilate, column by column, the first M columns of the */
/*           matrix G2 of negative generators, using Householder */
/*           transformations and modified hyperbolic plane rotations. */
/*           In the DLARF calls, ITAU is a pointer to the workspace */
/*           array. */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
		beta = dwork[ingc];
		dwork[ingc] = 1.;
		ingp = ingc + nrg;
		i__2 = nr - j;
		dlarf_("Left", &nrg, &i__2, &dwork[ingc], &c__1, &tau, &dwork[
			ingp], &nrg, &dwork[itau], (ftnlen)4);
		dwork[ingc] = beta;

/*              Compute the coefficients of the modified hyperbolic */
/*              rotation. */

		ma02fd_(&dwork[ipg + (j - 1) * (nrg + 1)], &dwork[ingc], &cs, 
			&sn, &ierr);
		if (ierr != 0) {

/*                 Error return: the matrix H'*H is not (numerically) */
/*                 positive definite. */

		    *info = 1;
		    return 0;
		}

		i__2 = (nr - 1) * nrg;
		i__3 = nrg;
		for (i__ = j * nrg; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ 
			+= i__3) {
		    dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn 
			    * dwork[ing + i__]) / cs;
		    dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			    dwork[ing + i__];
/* L200: */
		}

		ingc = ingp;
/* L210: */
	    }

/*           Save one block row of R, and shift the generators for the */
/*           calculation of the following row. */

	    dlacpy_("Upper", m, &nr, &dwork[ipg], &nrg, &r__[r_offset], ldr, (
		    ftnlen)5);

	    i__1 = mnrg;
	    i__3 = -mnrg;
	    for (i__ = (mmnobr - *m) * nrg; i__3 < 0 ? i__ >= i__1 : i__ <= 
		    i__1; i__ += i__3) {
		dlacpy_("Full", m, m, &dwork[ipg + i__ - mnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
/* L220: */
	    }

	    i__3 = (mmnobr + *l) * nrg;
	    i__1 = -lnrg;
	    for (i__ = (nr - *l) * nrg; i__1 < 0 ? i__ >= i__3 : i__ <= i__3; 
		    i__ += i__1) {
		dlacpy_("Full", m, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
/* L230: */
	    }

	    dlaset_("Full", m, l, &c_b111, &c_b111, &dwork[ipgc], &nrg, (
		    ftnlen)4);

/*           Update the input part of generators using Schur algorithm. */
/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*NOBR*(M+L)-M. */

	    jds = mnrg;
	    icol = *m;

	    i__1 = nobr2;
	    for (k = 2; k <= i__1; ++k) {
		i__3 = nr - icol - *m;
		i__2 = *l + 1;
		mb04od_("Full", m, &i__3, &i__2, &dwork[ipg + jds], &nrg, &
			dwork[ipy + jds], &nrg, &dwork[ipg + jds + mnrg], &
			nrg, &dwork[ipy + jds + mnrg], &nrg, &dwork[itau], &
			dwork[jwork], (ftnlen)4);

		i__3 = *m;
		for (j = 1; j <= i__3; ++j) {
		    icj = icol + j;
		    dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau)
			    ;
		    beta = dwork[ingc];
		    dwork[ingc] = 1.;
		    ingp = ingc + nrg;
		    i__2 = nr - icj;
		    dlarf_("Left", &nrg, &i__2, &dwork[ingc], &c__1, &tau, &
			    dwork[ingp], &nrg, &dwork[itau], (ftnlen)4);
		    dwork[ingc] = beta;

/*                 Compute the coefficients of the modified hyperbolic */
/*                 rotation. */

		    ma02fd_(&dwork[ipg + j - 1 + (icj - 1) * nrg], &dwork[
			    ingc], &cs, &sn, &ierr);
		    if (ierr != 0) {

/*                    Error return: the matrix H'*H is not (numerically) */
/*                    positive definite. */

			*info = 1;
			return 0;
		    }

		    i__2 = (nr - 1) * nrg;
		    i__4 = nrg;
		    for (i__ = icj * nrg; i__4 < 0 ? i__ >= i__2 : i__ <= 
			    i__2; i__ += i__4) {
			dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] 
				- sn * dwork[ing + i__]) / cs;
			dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + 
				cs * dwork[ing + i__];
/* L240: */
		    }

		    ingc = ingp;
/* L250: */
		}

/*              Save one block row of R, and shift the generators for the */
/*              calculation of the following row. */

		i__3 = nr - icol;
		dlacpy_("Upper", m, &i__3, &dwork[ipg + jds], &nrg, &r__[icol 
			+ 1 + (icol + 1) * r_dim1], ldr, (ftnlen)5);
		icol += *m;

		i__3 = icol * nrg;
		i__4 = -mnrg;
		for (i__ = (mmnobr - *m) * nrg; i__4 < 0 ? i__ >= i__3 : i__ 
			<= i__3; i__ += i__4) {
		    dlacpy_("Full", m, m, &dwork[ipg + i__ - mnrg], &nrg, &
			    dwork[ipg + i__], &nrg, (ftnlen)4);
/* L260: */
		}

		i__4 = (mmnobr + *l) * nrg;
		i__3 = -lnrg;
		for (i__ = (nr - *l) * nrg; i__3 < 0 ? i__ >= i__4 : i__ <= 
			i__4; i__ += i__3) {
		    dlacpy_("Full", m, l, &dwork[ipg + i__ - lnrg], &nrg, &
			    dwork[ipg + i__], &nrg, (ftnlen)4);
/* L270: */
		}

		dlaset_("Full", m, l, &c_b111, &c_b111, &dwork[ipgc], &nrg, (
			ftnlen)4);
		jds += mnrg;
/* L280: */
	    }

	}

/*        Process the output part of the generators. */

	jwork = itau + *l;

/*        Reduce the first L columns of the submatrix */
/*        G1(1:M+L+1,2*M*NOBR+1:2*(M+L)*NOBR) to upper triangular form. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*L; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+L+L*NB. */

	ingc = ing + mmnobr * nrg;
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(&nrg, l, &dwork[ipgc], &nrg, &dwork[itau], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
	i__1 = maxwrk, i__3 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__3);

/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+L*2*NOBR; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+L+(L*2*NOBR-L)*NB. */

	i__1 = llnobr - *l;
	i__3 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &nrg, &i__1, l, &dwork[ipgc], &nrg, &
		dwork[itau], &dwork[ipgc + lnrg], &nrg, &dwork[jwork], &i__3, 
		&ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
	i__1 = maxwrk, i__3 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__3);

/*        Annihilate, column by column, the first L columns of the */
/*        output part of the matrix G2 of negative generators, using */
/*        Householder transformations and modified hyperbolic rotations. */

	i__1 = *l;
	for (j = 1; j <= i__1; ++j) {
	    dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
	    beta = dwork[ingc];
	    dwork[ingc] = 1.;
	    ingp = ingc + nrg;
	    i__3 = llnobr - j;
	    dlarf_("Left", &nrg, &i__3, &dwork[ingc], &c__1, &tau, &dwork[
		    ingp], &nrg, &dwork[itau], (ftnlen)4);
	    dwork[ingc] = beta;

/*           Compute the coefficients of the modified hyperbolic */
/*           rotation. */

	    ma02fd_(&dwork[ipgc + (j - 1) * (nrg + 1)], &dwork[ingc], &cs, &
		    sn, &ierr);
	    if (ierr != 0) {

/*              Error return: the matrix H'*H is not (numerically) */
/*              positive definite. */

		*info = 1;
		return 0;
	    }

	    i__3 = (nr - 1) * nrg;
	    i__4 = nrg;
	    for (i__ = (j + mmnobr) * nrg; i__4 < 0 ? i__ >= i__3 : i__ <= 
		    i__3; i__ += i__4) {
		dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn * 
			dwork[ing + i__]) / cs;
		dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			dwork[ing + i__];
/* L290: */
	    }

	    ingc = ingp;
/* L300: */
	}

/*        Save one block row of R, and shift the generators for the */
/*        calculation of the following row. */

	dlacpy_("Upper", l, &llnobr, &dwork[ipgc], &nrg, &r__[mmnobr + 1 + (
		mmnobr + 1) * r_dim1], ldr, (ftnlen)5);

	i__1 = (mmnobr + *l) * nrg;
	i__4 = -lnrg;
	for (i__ = (nr - *l) * nrg; i__4 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__4) {
	    dlacpy_("Full", l, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[ipg 
		    + i__], &nrg, (ftnlen)4);
/* L310: */
	}

/*        Update the output part of generators using the Schur algorithm. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*NOBR*L-L. */

	jds = lnrg;
	icol = *l;

	i__4 = nobr2;
	for (k = 2; k <= i__4; ++k) {
	    i__1 = llnobr - icol - *l;
	    i__3 = *m + 1;
	    mb04od_("Full", l, &i__1, &i__3, &dwork[ipgc + jds], &nrg, &dwork[
		    ipgc + *l + jds], &nrg, &dwork[ipgc + jds + lnrg], &nrg, &
		    dwork[ipgc + *l + jds + lnrg], &nrg, &dwork[itau], &dwork[
		    jwork], (ftnlen)4);

	    i__1 = *l;
	    for (j = 1; j <= i__1; ++j) {
		icj = icol + j;
		dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
		beta = dwork[ingc];
		dwork[ingc] = 1.;
		ingp = ingc + nrg;
		i__3 = llnobr - icj;
		dlarf_("Left", &nrg, &i__3, &dwork[ingc], &c__1, &tau, &dwork[
			ingp], &nrg, &dwork[itau], (ftnlen)4);
		dwork[ingc] = beta;

/*              Compute the coefficients of the modified hyperbolic */
/*              rotation. */

		ma02fd_(&dwork[ipgc + j - 1 + (icj - 1) * nrg], &dwork[ingc], 
			&cs, &sn, &ierr);
		if (ierr != 0) {

/*                 Error return: the matrix H'*H is not (numerically) */
/*                 positive definite. */

		    *info = 1;
		    return 0;
		}

		i__3 = (nr - 1) * nrg;
		i__2 = nrg;
		for (i__ = (icj + mmnobr) * nrg; i__2 < 0 ? i__ >= i__3 : i__ 
			<= i__3; i__ += i__2) {
		    dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn 
			    * dwork[ing + i__]) / cs;
		    dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			    dwork[ing + i__];
/* L320: */
		}

		ingc = ingp;
/* L330: */
	    }

/*           Save one block row of R, and shift the generators for the */
/*           calculation of the following row. */

	    i__1 = llnobr - icol;
	    dlacpy_("Upper", l, &i__1, &dwork[ipgc + jds], &nrg, &r__[mmnobr 
		    + icol + 1 + (mmnobr + icol + 1) * r_dim1], ldr, (ftnlen)
		    5);

	    i__1 = (mmnobr + icol) * nrg;
	    i__2 = -lnrg;
	    for (i__ = (nr - *l) * nrg; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; 
		    i__ += i__2) {
		dlacpy_("Full", l, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
/* L340: */
	    }

	    icol += *l;
	    jds += lnrg;
/* L350: */
	}

	if (moesp && *m > 0) {

/*           For the MOESP algorithm, interchange the past and future */
/*           input parts of the R factor, and compute the new R factor */
/*           using a specialized QR factorization.  A tailored fast */
/*           QR factorization for the MOESP algorithm could be slightly */
/*           more efficient. */

	    i__4 = mnobr;
	    for (j = 1; j <= i__4; ++j) {
		dswap_(&j, &r__[j * r_dim1 + 1], &c__1, &r__[(mnobr + j) * 
			r_dim1 + 1], &c__1);
		dcopy_(&mnobr, &r__[j + 1 + (mnobr + j) * r_dim1], &c__1, &
			r__[j + 1 + j * r_dim1], &c__1);
		i__2 = mmnobr - j;
		dcopy_(&i__2, dum, &c__0, &r__[j + 1 + (mnobr + j) * r_dim1], 
			&c__1);
/* L360: */
	    }

/*           Triangularize the first two block columns (using structure), */
/*           and apply the transformation to the corresponding part of */
/*           the remaining block columns. */
/*           Workspace: need 2*(M+L)*NOBR. */

	    itau = 1;
	    jwork = itau + mmnobr;
	    i__4 = mnobr - 1;
	    i__2 = *ldwork - jwork + 1;
	    mb04id_(&mmnobr, &mmnobr, &i__4, &llnobr, &r__[r_offset], ldr, &
		    r__[(mmnobr + 1) * r_dim1 + 1], ldr, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);
/* Computing MAX */
	    i__4 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__4,i__2);
	}
    }

    nsmpsm = 0;
    icycle = 1;

/*     Return optimal workspace in  DWORK(1). */

    dwork[1] = (doublereal) maxwrk;
    maxwrk = 1;
    return 0;

/* *** Last line of IB01MY *** */
} /* ib01my_ */

