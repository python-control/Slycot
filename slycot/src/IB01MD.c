/* IB01MD.f -- translated by f2c (version 20100827).
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
static integer c_n1 = -1;
static doublereal c_b24 = -1.;
static doublereal c_b29 = 1.;
static doublereal c_b207 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ib01md_(char *meth, char *alg, char *batch, char *conct, 
	integer *nobr, integer *m, integer *l, integer *nsmp, doublereal *u, 
	integer *ldu, doublereal *y, integer *ldy, doublereal *r__, integer *
	ldr, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen meth_len, ftnlen alg_len, ftnlen 
	batch_len, ftnlen conct_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, id, jd, ii, nr, ns;
    static doublereal dum[1];
    static integer nsf;
    static doublereal upd;
    static integer inu, nsl, iny;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer icol, ierr, itau, init;
    static logical last;
    static doublereal temp;
    static integer irev;
    static logical linr, n4sid;
    static integer nobr2;
    static logical chalg;
    extern /* Subroutine */ int mb04od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ib01my_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer nobr21;
    static logical qralg;
    static integer initi, lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp;
    static integer lldrw, mldrw;
    static logical first;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static integer nobrm1, ishft2;
    static logical onebch, connec;
    static integer icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static logical fqralg;
    static integer ncycle, inicyc;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer nicycl;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer lmnobr, mmnobr;
    static logical interm;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer ishftu, nslast, ldrwrk, ishfty, minwrk, maxwrk, ldrwmx, 
	    nsmpsm;


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
/*     block Hankel matrices using input-output data.  The input-output */
/*     data can, optionally, be processed sequentially. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     ALG     CHARACTER*1 */
/*             Specifies the algorithm for computing the triangular */
/*             factor R, as follows: */
/*             = 'C':  Cholesky algorithm applied to the correlation */
/*                     matrix of the input-output data; */
/*             = 'F':  Fast QR algorithm; */
/*             = 'Q':  QR algorithm applied to the concatenated block */
/*                     Hankel matrices. */

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
/*             (In the MOESP theory,  NOBR  should be larger than  n, */
/*             the estimated dimension of state vector.) */

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

/*     R       (output or input/output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On exit, if INFO = 0 and ALG = 'Q', or (ALG = 'C' or 'F', */
/*             and BATCH = 'L' or 'O'), the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of */
/*             this array contains the (current) upper triangular factor */
/*             R from the QR factorization of the concatenated block */
/*             Hankel matrices. The diagonal elements of R are positive */
/*             when the Cholesky algorithm was successfully used. */
/*             On exit, if ALG = 'C' and BATCH = 'F' or 'I', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the current upper triangular part of the */
/*             correlation matrix in sequential data processing. */
/*             If ALG = 'F' and BATCH = 'F' or 'I', the array R is not */
/*             referenced. */
/*             On entry, if ALG = 'C', or ALG = 'Q', and BATCH = 'I' or */
/*             'L', the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  upper */
/*             triangular part of this array must contain the upper */
/*             triangular matrix R computed at the previous call of this */
/*             routine in sequential data processing. The array R need */
/*             not be set on entry if ALG = 'F' or if BATCH = 'F' or 'O'. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= M+L, if ALG = 'F'; */
/*             LIWORK >= 0,   if ALG = 'C' or 'Q'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             Let */
/*             k = 0,               if CONCT = 'N' and ALG = 'C' or 'Q'; */
/*             k = 2*NOBR-1,        if CONCT = 'C' and ALG = 'C' or 'Q'; */
/*             k = 2*NOBR*(M+L+1),  if CONCT = 'N' and ALG = 'F'; */
/*             k = 2*NOBR*(M+L+2),  if CONCT = 'C' and ALG = 'F'. */
/*             The first (M+L)*k elements of  DWORK  should be preserved */
/*             during successive calls of the routine with  BATCH = 'F' */
/*             or  'I',  till the final call with  BATCH = 'L'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= (4*NOBR-2)*(M+L), if ALG = 'C', BATCH <> 'O' and */
/*                                     CONCT = 'C'; */
/*             LDWORK >= 1,            if ALG = 'C', BATCH = 'O' or */
/*                                     CONCT = 'N'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+3), if ALG = 'F', */
/*                                     BATCH <> 'O' and CONCT = 'C'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+1), if ALG = 'F', */
/*                                     BATCH = 'F', 'I' and CONCT = 'N'; */
/*             LDWORK >= (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR, if ALG = 'F', */
/*                                     BATCH = 'L' and CONCT = 'N', or */
/*                                     BATCH = 'O'; */
/*             LDWORK >= 4*(M+L)*NOBR, if ALG = 'Q', BATCH = 'F' or 'O', */
/*                                     and LDR >= NS = NSMP - 2*NOBR + 1; */
/*             LDWORK >= 6*(M+L)*NOBR, if ALG = 'Q', BATCH = 'F' or 'O', */
/*                                     and LDR < NS, or BATCH = 'I' or */
/*                                     'L' and CONCT = 'N'; */
/*             LDWORK >= 4*(NOBR+1)*(M+L)*NOBR, if ALG = 'Q', BATCH = 'I' */
/*                                     or 'L' and CONCT = 'C'. */
/*             The workspace used for ALG = 'Q' is */
/*                       LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*             where LDRWRK = LDWORK/(2*(M+L)*NOBR) - 2; recommended */
/*             value LDRWRK = NS, assuming a large enough cache size. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the number of 100 cycles in sequential data */
/*                   processing has been exhausted without signaling */
/*                   that the last block of data was get; the cycle */
/*                   counter was reinitialized; */
/*             = 2:  a fast algorithm was requested (ALG = 'C' or 'F'), */
/*                   but it failed, and the QR algorithm was then used */
/*                   (non-sequential data processing). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  a fast algorithm was requested (ALG = 'C', or 'F') */
/*                   in sequential data processing, but it failed. The */
/*                   routine can be repeatedly called again using the */
/*                   standard QR algorithm. */

/*     METHOD */

/*     1) For non-sequential data processing using QR algorithm, a */
/*     t x 2(m+l)s  matrix H is constructed, where */

/*          H = [ Uf'         Up'      Y'      ],  for METH = 'M', */
/*                  s+1,2s,t    1,s,t   1,2s,t */

/*          H = [ U'       Y'      ],              for METH = 'N', */
/*                 1,2s,t   1,2s,t */

/*     and  Up     , Uf        , U      , and  Y        are block Hankel */
/*            1,s,t    s+1,2s,t   1,2s,t        1,2s,t */
/*     matrices defined in terms of the input and output data [3]. */
/*     A QR factorization is used to compress the data. */
/*     The fast QR algorithm uses a QR factorization which exploits */
/*     the block-Hankel structure. Actually, the Cholesky factor of H'*H */
/*     is computed. */

/*     2) For sequential data processing using QR algorithm, the QR */
/*     decomposition is done sequentially, by updating the upper */
/*     triangular factor  R.  This is also performed internally if the */
/*     workspace is not large enough to accommodate an entire batch. */

/*     3) For non-sequential or sequential data processing using */
/*     Cholesky algorithm, the correlation matrix of input-output data is */
/*     computed (sequentially, if requested), taking advantage of the */
/*     block Hankel structure [7].  Then, the Cholesky factor of the */
/*     correlation matrix is found, if possible. */

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

/*     [5] Peternell, K., Scherrer, W. and Deistler, M. */
/*         Statistical Analysis of Novel Subspace Identification Methods. */
/*         Signal Processing, 52, pp. 161-177, 1996. */

/*     [6] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     [7] Sima, V. */
/*         Cholesky or QR Factorization for Data Compression in */
/*         Subspace-based Identification ? */
/*         Proceedings of the Second NICONET Workshop on ``Numerical */
/*         Control Software: SLICOT, a Useful Tool in Industry'', */
/*         December 3, 1999, INRIA Rocquencourt, France, pp. 75-80, 1999. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable (when QR algorithm is */
/*     used), reliable and efficient. The fast Cholesky or QR algorithms */
/*     are more efficient, but the accuracy could diminish by forming the */
/*     correlation matrix. */
/*                                        2 */
/*     The QR algorithm needs 0(t(2(m+l)s) ) floating point operations. */
/*                                           2              3 */
/*     The Cholesky algorithm needs 0(2t(m+l) s)+0((2(m+l)s) ) floating */
/*     point operations. */
/*                                          2           3 2 */
/*     The fast QR algorithm needs 0(2t(m+l) s)+0(4(m+l) s ) floating */
/*     point operations. */

/*     FURTHER COMMENTS */

/*     For ALG = 'Q', BATCH = 'O' and LDR < NS, or BATCH <> 'O', the */
/*     calculations could be rather inefficient if only minimal workspace */
/*     (see argument LDWORK) is provided. It is advisable to provide as */
/*     much workspace as possible. Almost optimal efficiency can be */
/*     obtained for  LDWORK = (NS+2)*(2*(M+L)*NOBR),  assuming that the */
/*     cache size is large enough to accommodate R, U, Y, and DWORK. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     Feb. 2000, Aug. 2000, Feb. 2004. */

/*     KEYWORDS */

/*     Cholesky decomposition, Hankel matrix, identification methods, */
/*     multivariable systems, QR decomposition. */

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
/*        ICYCLE  is used to count the cycles for  BATCH = 'I'. It is */
/*                reinitialized at each MAXCYC cycles. */
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
    fqralg = lsame_(alg, "F", (ftnlen)1, (ftnlen)1);
    qralg = lsame_(alg, "Q", (ftnlen)1, (ftnlen)1);
    chalg = lsame_(alg, "C", (ftnlen)1, (ftnlen)1);
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
    lmnobr = lnobr + mnobr;
    mmnobr = mnobr + mnobr;
    nobrm1 = *nobr - 1;
    nobr21 = *nobr + nobrm1;
    nobr2 = nobr21 + 1;
    *iwarn = 0;
    *info = 0;
    ierr = 0;
    if (first) {
	icycle = 1;
	maxwrk = 1;
	nsmpsm = 0;
    }
    nsmpsm += *nsmp;
    nr = lmnobr + lmnobr;

/*     Check the scalar input parameters. */

    if (! (moesp || n4sid)) {
	*info = -1;
    } else if (! (fqralg || qralg || chalg)) {
	*info = -2;
    } else if (! (first || interm || last)) {
	*info = -3;
    } else if (! onebch) {
	if (! (connec || lsame_(conct, "N", (ftnlen)1, (ftnlen)1))) {
	    *info = -4;
	}
    }
    if (*info == 0) {
	if (*nobr <= 0) {
	    *info = -5;
	} else if (*m < 0) {
	    *info = -6;
	} else if (*l <= 0) {
	    *info = -7;
	} else if (*nsmp < nobr2 || last && nsmpsm < nr + nobr21) {
	    *info = -8;
	} else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
	    *info = -10;
	} else if (*ldy < *nsmp) {
	    *info = -12;
	} else if (*ldr < nr) {
	    *info = -14;
	} else {

/*           Compute workspace. */
/*           (Note: Comments in the code beginning "Workspace:" describe */
/*           the minimal amount of workspace needed at that point in the */
/*           code, as well as the preferred amount for good performance. */
/*           NB refers to the optimal block size for the immediately */
/*           following subroutine, as returned by ILAENV.) */

	    ns = *nsmp - nobr21;
	    if (chalg) {
		if (! onebch && connec) {
		    minwrk = nr - *m - *l << 1;
		} else {
		    minwrk = 1;
		}
	    } else if (fqralg) {
		if (! onebch && connec) {
		    minwrk = nr * (*m + *l + 3);
		} else if (first || interm) {
		    minwrk = nr * (*m + *l + 1);
		} else {
		    minwrk = (nr << 1) * (*m + *l + 1) + nr;
		}
	    } else {
		minwrk = nr << 1;
		maxwrk = nr + nr * ilaenv_(&c__1, "DGEQRF", " ", &ns, &nr, &
			c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
		if (first) {
		    if (*ldr < ns) {
			minwrk += nr;
			maxwrk = ns * nr + maxwrk;
		    }
		} else {
		    if (connec) {
			minwrk *= *nobr + 1;
		    } else {
			minwrk += nr;
		    }
		    maxwrk = ns * nr + maxwrk;
		}
	    }
	    maxwrk = max(minwrk,maxwrk);

	    if (*ldwork < minwrk) {
		*info = -17;
		dwork[1] = (doublereal) minwrk;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01MD", &i__1, (ftnlen)6);
	return 0;
    }

    if (chalg) {

/*        Compute the  R  factor from a Cholesky factorization of the */
/*        input-output data correlation matrix. */

/*        Set the parameters for constructing the correlations of the */
/*        current block. */

	ldrwrk = (nobr2 << 1) - 2;
	if (first) {
	    upd = 0.;
	} else {
	    upd = 1.;
	}

	if (! first && connec) {

/*           Restore the saved (M+L)*(2*NOBR-1) "connection" elements of */
/*           U  and  Y  into their appropriate position in sequential */
/*           processing. The process is performed column-wise, in */
/*           reverse order, first for  Y  and then for  U. */
/*           Workspace: need   (4*NOBR-2)*(M+L). */

	    irev = nr - *m - *l - nobr21 + 1;
	    icol = (nr - *m - *l << 1) - ldrwrk + 1;

	    i__1 = *m + *l;
	    for (j = 2; j <= i__1; ++j) {
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
		    dwork[icol + i__] = dwork[irev + i__];
/* L5: */
		}
		irev -= nobr21;
		icol -= ldrwrk;
/* L10: */
	    }

	    if (*m > 0) {
		dlacpy_("Full", &nobr21, m, &u[u_offset], ldu, &dwork[nobr2], 
			&ldrwrk, (ftnlen)4);
	    }
	    dlacpy_("Full", &nobr21, l, &y[y_offset], ldy, &dwork[ldrwrk * *m 
		    + nobr2], &ldrwrk, (ftnlen)4);
	}

	if (*m > 0) {

/*           Let  Guu(i,j) = Guu0(i,j) + u_i*u_j' + u_(i+1)*u_(j+1)' + */
/*                                 ... + u_(i+NS-1)*u_(j+NS-1)', */
/*           where  u_i'  is the i-th row of  U,  j = 1 : 2s,  i = 1 : j, */
/*           NS = NSMP - 2s + 1,  and  Guu0(i,j)  is a zero matrix for */
/*           BATCH = 'O' or 'F', and it is the matrix Guu(i,j) computed */
/*           till the current block for BATCH = 'I' or 'L'. The matrix */
/*           Guu(i,j)  is  m-by-m,  and  Guu(j,j)  is symmetric. The */
/*           upper triangle of the U-U correlations,  Guu,  is computed */
/*           (or updated) column-wise in the array  R,  that is, in the */
/*           order  Guu(1,1),  Guu(1,2),  Guu(2,2),  ...,  Guu(2s,2s). */
/*           Only the submatrices of the first block-row are fully */
/*           computed (or updated). The remaining ones are determined */
/*           exploiting the block-Hankel structure, using the updating */
/*           formula */

/*           Guu(i+1,j+1) = Guu0(i+1,j+1) - Guu0(i,j) + Guu(i,j) + */
/*                                 u_(i+NS)*u_(j+NS)' - u_i*u_j'. */

	    if (! first) {

/*              Subtract the contribution of the previous block of data */
/*              in sequential processing. The columns must be processed */
/*              in backward order. */

		for (i__ = nobr21 * *m; i__ >= 1; --i__) {
		    daxpy_(&i__, &c_b24, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
			    m + 1 + (*m + i__) * r_dim1], &c__1);
/* L20: */
		}

	    }

/*           Compute/update  Guu(1,1). */

	    if (! first && connec) {
		dsyrk_("Upper", "Transpose", m, &nobr21, &c_b29, &dwork[1], &
			ldrwrk, &upd, &r__[r_offset], ldr, (ftnlen)5, (ftnlen)
			9);
	    }
	    dsyrk_("Upper", "Transpose", m, &ns, &c_b29, &u[u_offset], ldu, &
		    upd, &r__[r_offset], ldr, (ftnlen)5, (ftnlen)9);

	    jd = 1;

	    if (first || ! connec) {

		i__1 = nobr2;
		for (j = 2; j <= i__1; ++j) {
		    jd += *m;
		    id = *m + 1;

/*                 Compute/update  Guu(1,j). */

		    dgemm_("Transpose", "NoTranspose", m, m, &ns, &c_b29, &u[
			    u_offset], ldu, &u[j + u_dim1], ldu, &upd, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guu(2:j,j), exploiting the */
/*                 block-Hankel structure. */

		    if (first) {

			i__2 = jd - 1;
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
			    dcopy_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*m + i__) * r_dim1], &c__1);
/* L30: */
			}

		    } else {

			i__2 = jd - 1;
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
			    daxpy_(&i__, &c_b29, &r__[i__ * r_dim1 + 1], &
				    c__1, &r__[*m + 1 + (*m + i__) * r_dim1], 
				    &c__1);
/* L40: */
			}

		    }

		    i__2 = j - 1;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			dger_(m, m, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				u[ns + j - 1 + u_dim1], ldu, &r__[id + jd * 
				r_dim1], ldr);
			dger_(m, m, &c_b24, &u[i__ - 1 + u_dim1], ldu, &u[j - 
				1 + u_dim1], ldu, &r__[id + jd * r_dim1], ldr)
				;
			id += *m;
/* L50: */
		    }

		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			daxpy_(&i__, &u[ns + j - 1 + i__ * u_dim1], &u[ns + j 
				- 1 + u_dim1], ldu, &r__[jd + (jd + i__ - 1) *
				 r_dim1], &c__1);
			d__1 = -u[j - 1 + i__ * u_dim1];
			daxpy_(&i__, &d__1, &u[j - 1 + u_dim1], ldu, &r__[jd 
				+ (jd + i__ - 1) * r_dim1], &c__1);
/* L60: */
		    }

/* L70: */
		}

	    } else {

		i__1 = nobr2;
		for (j = 2; j <= i__1; ++j) {
		    jd += *m;
		    id = *m + 1;

/*                 Compute/update  Guu(1,j)  for sequential processing */
/*                 with connected blocks. */

		    dgemm_("Transpose", "NoTranspose", m, m, &nobr21, &c_b29, 
			    &dwork[1], &ldrwrk, &dwork[j], &ldrwrk, &upd, &
			    r__[jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);
		    dgemm_("Transpose", "NoTranspose", m, m, &ns, &c_b29, &u[
			    u_offset], ldu, &u[j + u_dim1], ldu, &c_b29, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guu(2:j,j)  for sequential processing */
/*                 with connected blocks, exploiting the block-Hankel */
/*                 structure. */

		    if (first) {

			i__2 = jd - 1;
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
			    dcopy_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*m + i__) * r_dim1], &c__1);
/* L80: */
			}

		    } else {

			i__2 = jd - 1;
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
			    daxpy_(&i__, &c_b29, &r__[i__ * r_dim1 + 1], &
				    c__1, &r__[*m + 1 + (*m + i__) * r_dim1], 
				    &c__1);
/* L90: */
			}

		    }

		    i__2 = j - 1;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			dger_(m, m, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				u[ns + j - 1 + u_dim1], ldu, &r__[id + jd * 
				r_dim1], ldr);
			dger_(m, m, &c_b24, &dwork[i__ - 1], &ldrwrk, &dwork[
				j - 1], &ldrwrk, &r__[id + jd * r_dim1], ldr);
			id += *m;
/* L100: */
		    }

		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			daxpy_(&i__, &u[ns + j - 1 + i__ * u_dim1], &u[ns + j 
				- 1 + u_dim1], ldu, &r__[jd + (jd + i__ - 1) *
				 r_dim1], &c__1);
			d__1 = -dwork[(i__ - 1) * ldrwrk + j - 1];
			daxpy_(&i__, &d__1, &dwork[j - 1], &ldrwrk, &r__[jd + 
				(jd + i__ - 1) * r_dim1], &c__1);
/* L110: */
		    }

/* L120: */
		}

	    }

	    if (last && moesp) {

/*              Interchange past and future parts for MOESP algorithm. */
/*              (Only the upper triangular parts are interchanged, and */
/*              the (1,2) part is transposed in-situ.) */

		temp = r__[r_dim1 + 1];
		r__[r_dim1 + 1] = r__[mnobr + 1 + (mnobr + 1) * r_dim1];
		r__[mnobr + 1 + (mnobr + 1) * r_dim1] = temp;

		i__1 = mnobr;
		for (j = 2; j <= i__1; ++j) {
		    dswap_(&j, &r__[j * r_dim1 + 1], &c__1, &r__[mnobr + 1 + (
			    mnobr + j) * r_dim1], &c__1);
		    i__2 = j - 1;
		    dswap_(&i__2, &r__[(mnobr + j) * r_dim1 + 1], &c__1, &r__[
			    j + (mnobr + 1) * r_dim1], ldr);
/* L130: */
		}

	    }

/*           Let  Guy(i,j) = Guy0(i,j) + u_i*y_j' + u_(i+1)*y_(j+1)' + */
/*                                 ... + u_(i+NS-1)*y_(j+NS-1)', */
/*           where  u_i'  is the i-th row of  U,  y_j'  is the j-th row */
/*           of  Y,  j = 1 : 2s,  i = 1 : 2s,  NS = NSMP - 2s + 1,  and */
/*           Guy0(i,j)  is a zero matrix for  BATCH = 'O' or 'F', and it */
/*           is the matrix Guy(i,j) computed till the current block for */
/*           BATCH = 'I' or 'L'.  Guy(i,j) is m-by-L. The U-Y */
/*           correlations,  Guy,  are computed (or updated) column-wise */
/*           in the array  R. Only the submatrices of the first block- */
/*           column and block-row are fully computed (or updated). The */
/*           remaining ones are determined exploiting the block-Hankel */
/*           structure, using the updating formula */

/*           Guy(i+1,j+1) = Guy0(i+1,j+1) - Guy0(i,j) + Guy(i,j) + */
/*                                 u_(i+NS)*y(j+NS)' - u_i*y_j'. */

	    ii = mmnobr - *m;
	    if (! first) {

/*              Subtract the contribution of the previous block of data */
/*              in sequential processing. The columns must be processed */
/*              in backward order. */

		i__1 = mmnobr + 1;
		for (i__ = nr - *l; i__ >= i__1; --i__) {
		    daxpy_(&ii, &c_b24, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
			    m + 1 + (*l + i__) * r_dim1], &c__1);
/* L140: */
		}

	    }

/*           Compute/update the first block-column of  Guy,  Guy(i,1). */

	    if (first || ! connec) {

		i__1 = nobr2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    i__ + u_dim1], ldu, &y[y_offset], ldy, &upd, &r__[
			    (i__ - 1) * *m + 1 + (mmnobr + 1) * r_dim1], ldr, 
			    (ftnlen)9, (ftnlen)11);
/* L150: */
		}

	    } else {

		i__1 = nobr2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dgemm_("Transpose", "NoTranspose", m, l, &nobr21, &c_b29, 
			    &dwork[i__], &ldrwrk, &dwork[ldrwrk * *m + 1], &
			    ldrwrk, &upd, &r__[(i__ - 1) * *m + 1 + (mmnobr + 
			    1) * r_dim1], ldr, (ftnlen)9, (ftnlen)11);
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    i__ + u_dim1], ldu, &y[y_offset], ldy, &c_b29, &
			    r__[(i__ - 1) * *m + 1 + (mmnobr + 1) * r_dim1], 
			    ldr, (ftnlen)9, (ftnlen)11);
/* L160: */
		}

	    }

	    jd = mmnobr + 1;

	    if (first || ! connec) {

		i__1 = nobr2;
		for (j = 2; j <= i__1; ++j) {
		    jd += *l;
		    id = *m + 1;

/*                 Compute/update  Guy(1,j). */

		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    u_offset], ldu, &y[j + y_dim1], ldy, &upd, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guy(2:2*s,j), exploiting the */
/*                 block-Hankel structure. */

		    if (first) {

			i__2 = jd - 1;
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
			    dcopy_(&ii, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*l + i__) * r_dim1], &c__1);
/* L170: */
			}

		    } else {

			i__2 = jd - 1;
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
			    daxpy_(&ii, &c_b29, &r__[i__ * r_dim1 + 1], &c__1,
				     &r__[*m + 1 + (*l + i__) * r_dim1], &
				    c__1);
/* L180: */
			}

		    }

		    i__2 = nobr2;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			dger_(m, l, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				y[ns + j - 1 + y_dim1], ldy, &r__[id + jd * 
				r_dim1], ldr);
			dger_(m, l, &c_b24, &u[i__ - 1 + u_dim1], ldu, &y[j - 
				1 + y_dim1], ldy, &r__[id + jd * r_dim1], ldr)
				;
			id += *m;
/* L190: */
		    }

/* L200: */
		}

	    } else {

		i__1 = nobr2;
		for (j = 2; j <= i__1; ++j) {
		    jd += *l;
		    id = *m + 1;

/*                 Compute/update  Guy(1,j)  for sequential processing */
/*                 with connected blocks. */

		    dgemm_("Transpose", "NoTranspose", m, l, &nobr21, &c_b29, 
			    &dwork[1], &ldrwrk, &dwork[ldrwrk * *m + j], &
			    ldrwrk, &upd, &r__[jd * r_dim1 + 1], ldr, (ftnlen)
			    9, (ftnlen)11);
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    u_offset], ldu, &y[j + y_dim1], ldy, &c_b29, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guy(2:2*s,j)  for sequential */
/*                 processing with connected blocks, exploiting the */
/*                 block-Hankel structure. */

		    if (first) {

			i__2 = jd - 1;
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
			    dcopy_(&ii, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*l + i__) * r_dim1], &c__1);
/* L210: */
			}

		    } else {

			i__2 = jd - 1;
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
			    daxpy_(&ii, &c_b29, &r__[i__ * r_dim1 + 1], &c__1,
				     &r__[*m + 1 + (*l + i__) * r_dim1], &
				    c__1);
/* L220: */
			}

		    }

		    i__2 = nobr2;
		    for (i__ = 2; i__ <= i__2; ++i__) {
			dger_(m, l, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				y[ns + j - 1 + y_dim1], ldy, &r__[id + jd * 
				r_dim1], ldr);
			dger_(m, l, &c_b24, &dwork[i__ - 1], &ldrwrk, &dwork[
				ldrwrk * *m + j - 1], &ldrwrk, &r__[id + jd * 
				r_dim1], ldr);
			id += *m;
/* L230: */
		    }

/* L240: */
		}

	    }

	    if (last && moesp) {

/*              Interchange past and future parts of U-Y correlations */
/*              for MOESP algorithm. */

		i__1 = nr;
		for (j = mmnobr + 1; j <= i__1; ++j) {
		    dswap_(&mnobr, &r__[j * r_dim1 + 1], &c__1, &r__[mnobr + 
			    1 + j * r_dim1], &c__1);
/* L250: */
		}

	    }
	}

/*        Let  Gyy(i,j) = Gyy0(i,j) + y_i*y_i' + y_(i+1)*y_(i+1)' + ... + */
/*                                    y_(i+NS-1)*y_(i+NS-1)', */
/*        where  y_i'  is the i-th row of  Y,  j = 1 : 2s,  i = 1 : j, */
/*        NS = NSMP - 2s + 1,  and  Gyy0(i,j)  is a zero matrix for */
/*        BATCH = 'O' or 'F', and it is the matrix Gyy(i,j) computed till */
/*        the current block for BATCH = 'I' or 'L'.  Gyy(i,j) is L-by-L, */
/*        and  Gyy(j,j)  is symmetric. The upper triangle of the Y-Y */
/*        correlations,  Gyy,  is computed (or updated) column-wise in */
/*        the corresponding part of the array  R,  that is, in the order */
/*        Gyy(1,1),  Gyy(1,2),  Gyy(2,2),  ...,  Gyy(2s,2s).  Only the */
/*        submatrices of the first block-row are fully computed (or */
/*        updated). The remaining ones are determined exploiting the */
/*        block-Hankel structure, using the updating formula */

/*        Gyy(i+1,j+1) = Gyy0(i+1,j+1) - Gyy0(i,j) + Gyy(i,j) + */
/*                              y_(i+NS)*y_(j+NS)' - y_i*y_j'. */

	jd = mmnobr + 1;

	if (! first) {

/*           Subtract the contribution of the previous block of data */
/*           in sequential processing. The columns must be processed in */
/*           backward order. */

	    i__1 = mmnobr + 1;
	    for (i__ = nr - *l; i__ >= i__1; --i__) {
		i__2 = i__ - mmnobr;
		daxpy_(&i__2, &c_b24, &r__[jd + i__ * r_dim1], &c__1, &r__[jd 
			+ *l + (*l + i__) * r_dim1], &c__1);
/* L260: */
	    }

	}

/*        Compute/update  Gyy(1,1). */

	if (! first && connec) {
	    dsyrk_("Upper", "Transpose", l, &nobr21, &c_b29, &dwork[ldrwrk * *
		    m + 1], &ldrwrk, &upd, &r__[jd + jd * r_dim1], ldr, (
		    ftnlen)5, (ftnlen)9);
	}
	dsyrk_("Upper", "Transpose", l, &ns, &c_b29, &y[y_offset], ldy, &upd, 
		&r__[jd + jd * r_dim1], ldr, (ftnlen)5, (ftnlen)9);

	if (first || ! connec) {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *l;
		id = mmnobr + *l + 1;

/*              Compute/update  Gyy(1,j). */

		dgemm_("Transpose", "NoTranspose", l, l, &ns, &c_b29, &y[
			y_offset], ldy, &y[j + y_dim1], ldy, &upd, &r__[
			mmnobr + 1 + jd * r_dim1], ldr, (ftnlen)9, (ftnlen)11)
			;

/*              Compute/update  Gyy(2:j,j), exploiting the block-Hankel */
/*              structure. */

		if (first) {

		    i__2 = jd - 1;
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
			i__3 = i__ - mmnobr;
			dcopy_(&i__3, &r__[mmnobr + 1 + i__ * r_dim1], &c__1, 
				&r__[mmnobr + *l + 1 + (*l + i__) * r_dim1], &
				c__1);
/* L270: */
		    }

		} else {

		    i__2 = jd - 1;
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
			i__3 = i__ - mmnobr;
			daxpy_(&i__3, &c_b29, &r__[mmnobr + 1 + i__ * r_dim1],
				 &c__1, &r__[mmnobr + *l + 1 + (*l + i__) * 
				r_dim1], &c__1);
/* L280: */
		    }

		}

		i__2 = j - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    dger_(l, l, &c_b29, &y[ns + i__ - 1 + y_dim1], ldy, &y[ns 
			    + j - 1 + y_dim1], ldy, &r__[id + jd * r_dim1], 
			    ldr);
		    dger_(l, l, &c_b24, &y[i__ - 1 + y_dim1], ldy, &y[j - 1 + 
			    y_dim1], ldy, &r__[id + jd * r_dim1], ldr);
		    id += *l;
/* L290: */
		}

		i__2 = *l;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    daxpy_(&i__, &y[ns + j - 1 + i__ * y_dim1], &y[ns + j - 1 
			    + y_dim1], ldy, &r__[jd + (jd + i__ - 1) * r_dim1]
			    , &c__1);
		    d__1 = -y[j - 1 + i__ * y_dim1];
		    daxpy_(&i__, &d__1, &y[j - 1 + y_dim1], ldy, &r__[jd + (
			    jd + i__ - 1) * r_dim1], &c__1);
/* L300: */
		}

/* L310: */
	    }

	} else {

	    i__1 = nobr2;
	    for (j = 2; j <= i__1; ++j) {
		jd += *l;
		id = mmnobr + *l + 1;

/*              Compute/update  Gyy(1,j)  for sequential processing with */
/*              connected blocks. */

		dgemm_("Transpose", "NoTranspose", l, l, &nobr21, &c_b29, &
			dwork[ldrwrk * *m + 1], &ldrwrk, &dwork[ldrwrk * *m + 
			j], &ldrwrk, &upd, &r__[mmnobr + 1 + jd * r_dim1], 
			ldr, (ftnlen)9, (ftnlen)11);
		dgemm_("Transpose", "NoTranspose", l, l, &ns, &c_b29, &y[
			y_offset], ldy, &y[j + y_dim1], ldy, &c_b29, &r__[
			mmnobr + 1 + jd * r_dim1], ldr, (ftnlen)9, (ftnlen)11)
			;

/*              Compute/update  Gyy(2:j,j)  for sequential processing */
/*              with connected blocks, exploiting the block-Hankel */
/*              structure. */

		if (first) {

		    i__2 = jd - 1;
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
			i__3 = i__ - mmnobr;
			dcopy_(&i__3, &r__[mmnobr + 1 + i__ * r_dim1], &c__1, 
				&r__[mmnobr + *l + 1 + (*l + i__) * r_dim1], &
				c__1);
/* L320: */
		    }

		} else {

		    i__2 = jd - 1;
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
			i__3 = i__ - mmnobr;
			daxpy_(&i__3, &c_b29, &r__[mmnobr + 1 + i__ * r_dim1],
				 &c__1, &r__[mmnobr + *l + 1 + (*l + i__) * 
				r_dim1], &c__1);
/* L330: */
		    }

		}

		i__2 = j - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    dger_(l, l, &c_b29, &y[ns + i__ - 1 + y_dim1], ldy, &y[ns 
			    + j - 1 + y_dim1], ldy, &r__[id + jd * r_dim1], 
			    ldr);
		    dger_(l, l, &c_b24, &dwork[ldrwrk * *m + i__ - 1], &
			    ldrwrk, &dwork[ldrwrk * *m + j - 1], &ldrwrk, &
			    r__[id + jd * r_dim1], ldr);
		    id += *l;
/* L340: */
		}

		i__2 = *l;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    daxpy_(&i__, &y[ns + j - 1 + i__ * y_dim1], &y[ns + j - 1 
			    + y_dim1], ldy, &r__[jd + (jd + i__ - 1) * r_dim1]
			    , &c__1);
		    d__1 = -dwork[ldrwrk * (*m + i__ - 1) + j - 1];
		    daxpy_(&i__, &d__1, &dwork[ldrwrk * *m + j - 1], &ldrwrk, 
			    &r__[jd + (jd + i__ - 1) * r_dim1], &c__1);
/* L350: */
		}

/* L360: */
	    }

	}

	if (! last) {
	    if (connec) {

/*              For sequential processing with connected data blocks, */
/*              save the remaining ("connection") elements of  U  and  Y */
/*              in the first  (M+L)*(2*NOBR-1)  locations of  DWORK. */

		if (*m > 0) {
		    dlacpy_("Full", &nobr21, m, &u[ns + 1 + u_dim1], ldu, &
			    dwork[1], &nobr21, (ftnlen)4);
		}
		dlacpy_("Full", &nobr21, l, &y[ns + 1 + y_dim1], ldy, &dwork[
			mmnobr - *m + 1], &nobr21, (ftnlen)4);
	    }

/*           Return to get new data. */

	    ++icycle;
	    if (icycle > 100) {
		*iwarn = 1;
	    }
	    return 0;

	} else {

/*           Try to compute the Cholesky factor of the correlation */
/*           matrix. */

	    dpotrf_("Upper", &nr, &r__[r_offset], ldr, &ierr, (ftnlen)5);
	    goto L370;
	}
    } else if (fqralg) {

/*        Compute the  R  factor from a fast QR factorization of the */
/*        input-output data correlation matrix. */

	ib01my_(meth, batch, conct, nobr, m, l, nsmp, &u[u_offset], ldu, &y[
		y_offset], ldy, &r__[r_offset], ldr, &iwork[1], &dwork[1], 
		ldwork, iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
	if (! last) {
	    return 0;
	}
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[1];
	maxwrk = max(i__1,i__2);
    }

L370:

    if (ierr != 0) {

/*        Error return from a fast factorization algorithm of the */
/*        input-output data correlation matrix. */

	if (onebch) {
	    qralg = TRUE_;
	    *iwarn = 2;
	    minwrk = nr << 1;
	    maxwrk = nr + nr * ilaenv_(&c__1, "DGEQRF", " ", &ns, &nr, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
	    if (*ldr < ns) {
		minwrk += nr;
		maxwrk = ns * nr + maxwrk;
	    }
	    maxwrk = max(minwrk,maxwrk);

	    if (*ldwork < minwrk) {
		*info = -17;

/*              Return: Not enough workspace. */

		dwork[1] = (doublereal) minwrk;
		i__1 = -(*info);
		xerbla_("IB01MD", &i__1, (ftnlen)6);
		return 0;
	    }
	} else {
	    *info = 1;
	    return 0;
	}
    }

    if (qralg) {

/*        Compute the  R  factor from a QR factorization of the matrix  H */
/*        of concatenated block Hankel matrices. */

/*        Construct the matrix  H. */

/*        Set the parameters for constructing the current segment of the */
/*        Hankel matrix, taking the available memory space into account. */
/*        INITI+1 points to the beginning rows of  U  and  Y  from which */
/*                data are taken when NCYCLE > 1 inner cycles are needed, */
/*                or for sequential processing with connected blocks. */
/*        LDRWMX is the number of rows that can fit in the working space. */
/*        LDRWRK is the actual number of rows processed in this space. */
/*        NSLAST is the number of samples to be processed at the last */
/*               inner cycle. */

	initi = 0;
	ldrwmx = *ldwork / nr - 2;
	ncycle = 1;
	nslast = *nsmp;
	linr = FALSE_;
	if (first) {
	    linr = *ldr >= ns;
	    ldrwrk = ns;
	} else if (connec) {
	    ldrwrk = *nsmp;
	} else {
	    ldrwrk = ns;
	}
	inicyc = 1;

	if (! linr) {
	    if (ldrwmx < ldrwrk) {

/*              Not enough working space for doing a single inner cycle. */
/*              NCYCLE inner cycles are to be performed for the current */
/*              data block using the working space. */

		ncycle = ldrwrk / ldrwmx;
		nslast = ldrwrk % ldrwmx;
		if (nslast != 0) {
		    ++ncycle;
		} else {
		    nslast = ldrwmx;
		}
		ldrwrk = ldrwmx;
		ns = ldrwrk;
		if (first) {
		    inicyc = 2;
		}
	    }
	    mldrw = *m * ldrwrk;
	    lldrw = *l * ldrwrk;
	    inu = mldrw * *nobr + 1;
	    iny = mldrw * nobr2 + 1;
	}

/*        Process the data given at the current call. */

	if (! first && connec) {

/*           Restore the saved (M+L)*(2*NOBR-1) "connection" elements of */
/*           U  and  Y  into their appropriate position in sequential */
/*           processing. The process is performed column-wise, in */
/*           reverse order, first for  Y  and then for  U. */

	    irev = nr - *m - *l - nobr21 + 1;
	    icol = iny + lldrw - ldrwrk;

	    i__1 = *l;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
		    dwork[icol + i__] = dwork[irev + i__];
/* L375: */
		}
		irev -= nobr21;
		icol -= ldrwrk;
/* L380: */
	    }

	    if (moesp) {
		icol = inu + mldrw - ldrwrk;
	    } else {
		icol = mldrw - ldrwrk + 1;
	    }

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
		    dwork[icol + i__] = dwork[irev + i__];
/* L385: */
		}
		irev -= nobr21;
		icol -= ldrwrk;
/* L390: */
	    }

	    if (moesp) {
		dlacpy_("Full", &nobrm1, m, &dwork[inu + *nobr], &ldrwrk, &
			dwork[1], &ldrwrk, (ftnlen)4);
	    }
	}

/*        Data compression using QR factorization. */

	if (first) {

/*           Non-sequential data processing or first block in */
/*           sequential data processing: */
/*           Use the general QR factorization algorithm. */

	    if (linr) {

/*              Put the input-output data in the array  R. */

		if (*m > 0) {
		    if (moesp) {

			i__1 = *nobr;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[*nobr + i__ + u_dim1], 
				    ldu, &r__[(*m * (i__ - 1) + 1) * r_dim1 + 
				    1], ldr, (ftnlen)4);
/* L400: */
			}

			i__1 = *nobr;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    r__[(mnobr + *m * (i__ - 1) + 1) * r_dim1 
				    + 1], ldr, (ftnlen)4);
/* L410: */
			}

		    } else {

			i__1 = nobr2;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    r__[(*m * (i__ - 1) + 1) * r_dim1 + 1], 
				    ldr, (ftnlen)4);
/* L420: */
			}

		    }
		}

		i__1 = nobr2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dlacpy_("Full", &ns, l, &y[i__ + y_dim1], ldy, &r__[(
			    mmnobr + *l * (i__ - 1) + 1) * r_dim1 + 1], ldr, (
			    ftnlen)4);
/* L430: */
		}

/*              Workspace: need   4*(M+L)*NOBR, */
/*                         prefer 2*(M+L)*NOBR+2*(M+L)*NOBR*NB. */

		itau = 1;
		jwork = itau + nr;
		i__1 = *ldwork - jwork + 1;
		dgeqrf_(&ns, &nr, &r__[r_offset], ldr, &dwork[itau], &dwork[
			jwork], &i__1, &ierr);
	    } else {

/*              Put the input-output data in the array  DWORK. */

		if (*m > 0) {
		    ishftu = 1;
		    if (moesp) {
			ishft2 = inu;

			i__1 = *nobr;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[*nobr + i__ + u_dim1], 
				    ldu, &dwork[ishftu], &ldrwrk, (ftnlen)4);
			    ishftu += mldrw;
/* L440: */
			}

			i__1 = *nobr;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    dwork[ishft2], &ldrwrk, (ftnlen)4);
			    ishft2 += mldrw;
/* L450: */
			}

		    } else {

			i__1 = nobr2;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    dwork[ishftu], &ldrwrk, (ftnlen)4);
			    ishftu += mldrw;
/* L460: */
			}

		    }
		}

		ishfty = iny;

		i__1 = nobr2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    dlacpy_("Full", &ns, l, &y[i__ + y_dim1], ldy, &dwork[
			    ishfty], &ldrwrk, (ftnlen)4);
		    ishfty += lldrw;
/* L470: */
		}

/*              Workspace: need   2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*                         prefer NS*2*(M+L)*NOBR + 2*(M+L)*NOBR */
/*                                                + 2*(M+L)*NOBR*NB, */
/*                         used LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*                         where  NS = NSMP - 2*NOBR + 1, */
/*                            LDRWRK = min(NS, LDWORK/(2*(M+L)*NOBR)-2). */

		itau = ldrwrk * nr + 1;
		jwork = itau + nr;
		i__1 = *ldwork - jwork + 1;
		dgeqrf_(&ns, &nr, &dwork[1], &ldrwrk, &dwork[itau], &dwork[
			jwork], &i__1, &ierr);
		i__1 = min(ns,nr);
		dlacpy_("Upper ", &i__1, &nr, &dwork[1], &ldrwrk, &r__[
			r_offset], ldr, (ftnlen)6);
	    }

	    if (ns < nr) {
		i__1 = nr - ns;
		i__2 = nr - ns;
		dlaset_("Upper ", &i__1, &i__2, &c_b207, &c_b207, &r__[ns + 1 
			+ (ns + 1) * r_dim1], ldr, (ftnlen)6);
	    }
	    initi += ns;
	}

	if (ncycle > 1 || ! first) {

/*           Remaining segments of the first data block or */
/*           remaining segments/blocks in sequential data processing: */
/*           Use a structure-exploiting QR factorization algorithm. */

	    nsl = ldrwrk;
	    if (! connec) {
		nsl = ns;
	    }
	    itau = ldrwrk * nr + 1;
	    jwork = itau + nr;

	    i__1 = ncycle;
	    for (nicycl = inicyc; nicycl <= i__1; ++nicycl) {

/*              INIT  denotes the beginning row where new data are put. */

		if (connec && nicycl == 1) {
		    init = nobr2;
		} else {
		    init = 1;
		}
		if (ncycle > 1 && nicycl == ncycle) {

/*                 Last samples in the last data segment of a block. */

		    ns = nslast;
		    nsl = nslast;
		}

/*              Put the input-output data in the array  DWORK. */

		nsf = ns;
		if (init > 1 && ncycle > 1) {
		    nsf -= nobr21;
		}
		if (*m > 0) {
		    ishftu = init;

		    if (moesp) {
			ishft2 = init + inu - 1;

			i__2 = *nobr;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    dlacpy_("Full", &nsf, m, &u[initi + *nobr + i__ + 
				    u_dim1], ldu, &dwork[ishftu], &ldrwrk, (
				    ftnlen)4);
			    ishftu += mldrw;
/* L480: */
			}

			i__2 = *nobr;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    dlacpy_("Full", &nsf, m, &u[initi + i__ + u_dim1],
				     ldu, &dwork[ishft2], &ldrwrk, (ftnlen)4);
			    ishft2 += mldrw;
/* L490: */
			}

		    } else {

			i__2 = nobr2;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    dlacpy_("Full", &nsf, m, &u[initi + i__ + u_dim1],
				     ldu, &dwork[ishftu], &ldrwrk, (ftnlen)4);
			    ishftu += mldrw;
/* L500: */
			}

		    }
		}

		ishfty = init + iny - 1;

		i__2 = nobr2;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    dlacpy_("Full", &nsf, l, &y[initi + i__ + y_dim1], ldy, &
			    dwork[ishfty], &ldrwrk, (ftnlen)4);
		    ishfty += lldrw;
/* L510: */
		}

		if (init > 1) {

/*                 Prepare the connection to the previous block of data */
/*                 in sequential processing. */

		    if (moesp && *m > 0) {
			dlacpy_("Full", nobr, m, &u[u_offset], ldu, &dwork[*
				nobr], &ldrwrk, (ftnlen)4);
		    }

/*                 Shift the elements from the connection to the previous */
/*                 block of data in sequential processing. */

		    if (*m > 0) {
			ishftu = mldrw + 1;

			if (moesp) {
			    ishft2 = mldrw + inu;

			    i__2 = nobrm1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				dlacpy_("Full", &nobr21, m, &dwork[ishftu - 
					mldrw + 1], &ldrwrk, &dwork[ishftu], &
					ldrwrk, (ftnlen)4);
				ishftu += mldrw;
/* L520: */
			    }

			    i__2 = nobrm1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				dlacpy_("Full", &nobr21, m, &dwork[ishft2 - 
					mldrw + 1], &ldrwrk, &dwork[ishft2], &
					ldrwrk, (ftnlen)4);
				ishft2 += mldrw;
/* L530: */
			    }

			} else {

			    i__2 = nobr21;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				dlacpy_("Full", &nobr21, m, &dwork[ishftu - 
					mldrw + 1], &ldrwrk, &dwork[ishftu], &
					ldrwrk, (ftnlen)4);
				ishftu += mldrw;
/* L540: */
			    }

			}
		    }

		    ishfty = lldrw + iny;

		    i__2 = nobr21;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			dlacpy_("Full", &nobr21, l, &dwork[ishfty - lldrw + 1]
				, &ldrwrk, &dwork[ishfty], &ldrwrk, (ftnlen)4)
				;
			ishfty += lldrw;
/* L550: */
		    }

		}

/*              Workspace: need LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR. */

		mb04od_("Full", &nr, &c__0, &nsl, &r__[r_offset], ldr, &dwork[
			1], &ldrwrk, dum, &nr, dum, &nr, &dwork[itau], &dwork[
			jwork], (ftnlen)4);
		initi += nsf;
/* L560: */
	    }

	}

	if (! last) {
	    if (connec) {

/*              For sequential processing with connected data blocks, */
/*              save the remaining ("connection") elements of  U  and  Y */
/*              in the first  (M+L)*(2*NOBR-1)  locations of  DWORK. */

		if (*m > 0) {
		    dlacpy_("Full", &nobr21, m, &u[initi + 1 + u_dim1], ldu, &
			    dwork[1], &nobr21, (ftnlen)4);
		}
		dlacpy_("Full", &nobr21, l, &y[initi + 1 + y_dim1], ldy, &
			dwork[mmnobr - *m + 1], &nobr21, (ftnlen)4);
	    }

/*           Return to get new data. */

	    ++icycle;
	    if (icycle <= 100) {
		return 0;
	    }
	    *iwarn = 1;
	    icycle = 1;

	}

    }

/*     Return optimal workspace in  DWORK(1). */

    dwork[1] = (doublereal) maxwrk;
    if (last) {
	icycle = 1;
	maxwrk = 1;
	nsmpsm = 0;
    }
    return 0;

/* *** Last line of IB01MD *** */
} /* ib01md_ */

