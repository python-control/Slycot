/* IB01ND.f -- translated by f2c (version 20100827).
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
static doublereal c_b37 = .66666666666666663;
static integer c__0 = 0;
static doublereal c_b50 = 0.;

/* Subroutine */ int ib01nd_(char *meth, char *jobd, integer *nobr, integer *
	m, integer *l, doublereal *r__, integer *ldr, doublereal *sv, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen meth_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, nr, nr2, nr3, nr4;
    static doublereal dum[1], eps;
    static integer rank, ierr, itau;
    static doublereal sval[3], toll;
    static integer rank1;
    static logical n4sid;
    static integer itau2, itau3;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb04id_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), mb04od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), mb03ud_(char *, char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen);
    static logical jobdm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04iy_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp;
    static integer jwork;
    static doublereal rcond1, rcond2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer llmnob, lmmnob, llnobr, lmnobr, mmnobr;
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal thresh;
    static integer nrsave;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal svlmax;


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

/*     To find the singular value decomposition (SVD) giving the system */
/*     order, using the triangular factor of the concatenated block */
/*     Hankel matrices. Related preliminary calculations needed for */
/*     computing the system matrices are also performed. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not the matrices B and D should later */
/*             be computed using the MOESP approach, as follows: */
/*             = 'M':  the matrices B and D should later be computed */
/*                     using the MOESP approach; */
/*             = 'N':  the matrices B and D should not be computed using */
/*                     the MOESP approach. */
/*             This parameter is not relevant for METH = 'N'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             block Hankel matrices.  NOBR > 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On entry, the leading 2*(M+L)*NOBR-by-2*(M+L)*NOBR upper */
/*             triangular part of this array must contain the upper */
/*             triangular factor R from the QR factorization of the */
/*             concatenated block Hankel matrices. Denote  R_ij, */
/*             i,j = 1:4,  the ij submatrix of  R,  partitioned by */
/*             M*NOBR,  M*NOBR,  L*NOBR,  and  L*NOBR  rows and columns. */
/*             On exit, if INFO = 0, the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the matrix S, the processed upper */
/*             triangular factor R, as required by other subroutines. */
/*             Specifically, let  S_ij, i,j = 1:4,  be the ij submatrix */
/*             of  S,  partitioned by  M*NOBR,  L*NOBR,  M*NOBR,  and */
/*             L*NOBR  rows and columns. The submatrix  S_22  contains */
/*             the matrix of left singular vectors needed subsequently. */
/*             Useful information is stored in  S_11  and in the */
/*             block-column  S_14 : S_44.  For METH = 'M' and JOBD = 'M', */
/*             the upper triangular part of  S_31  contains the upper */
/*             triangular factor in the QR factorization of the matrix */
/*             R_1c = [ R_12'  R_22'  R_11' ]',  and  S_12  contains the */
/*             corresponding leading part of the transformed matrix */
/*             R_2c = [ R_13'  R_23'  R_14' ]'.  For  METH = 'N',  the */
/*             subarray  S_41 : S_43  contains the transpose of the */
/*             matrix contained in  S_14 : S_34. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= MAX( 2*(M+L)*NOBR, 3*M*NOBR ), */
/*                                  for METH = 'M' and JOBD = 'M'; */
/*             LDR >= 2*(M+L)*NOBR, for METH = 'M' and JOBD = 'N' or */
/*                                  for METH = 'N'. */

/*     SV      (output) DOUBLE PRECISION array, dimension ( L*NOBR ) */
/*             The singular values of the relevant part of the triangular */
/*             factor from the QR factorization of the concatenated block */
/*             Hankel matrices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  an m-by-n matrix whose estimated */
/*             condition number is less than  1/TOL  is considered to */
/*             be of full rank.  If the user sets  TOL <= 0,  then an */
/*             implicitly computed, default tolerance, defined by */
/*             TOLDEF = m*n*EPS,  is used instead, where  EPS  is the */
/*             relative machine precision (see LAPACK Library routine */
/*             DLAMCH). */
/*             This parameter is not used for  METH = 'M'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension ((M+L)*NOBR) */
/*             This parameter is not referenced for METH = 'M'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and, for  METH = 'N',  DWORK(2)  and  DWORK(3) */
/*             contain the reciprocal condition numbers of the */
/*             triangular factors of the matrices  U_f  and  r_1  [6]. */
/*             On exit, if  INFO = -12,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( (2*M-1)*NOBR, (M+L)*NOBR, 5*L*NOBR ), */
/*                                         if METH = 'M' and JOBD = 'M'; */
/*             LDWORK >=  5*L*NOBR,        if METH = 'M' and JOBD = 'N'; */
/*             LDWORK >=  5*(M+L)*NOBR+1,  if METH = 'N'. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problems with coefficient matrix */
/*                   U_f,  used for computing the weighted oblique */
/*                   projection (for METH = 'N'), have a rank-deficient */
/*                   coefficient matrix; */
/*             = 5:  the least squares problem with coefficient matrix */
/*                   r_1  [6], used for computing the weighted oblique */
/*                   projection (for METH = 'N'), has a rank-deficient */
/*                   coefficient matrix. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     A singular value decomposition (SVD) of a certain matrix is */
/*     computed, which reveals the order  n  of the system as the number */
/*     of "non-zero" singular values. For the MOESP approach, this matrix */
/*     is  [ R_24'  R_34' ]' := R(ms+1:(2m+l)s,(2m+l)s+1:2(m+l)s), */
/*     where  R  is the upper triangular factor  R  constructed by SLICOT */
/*     Library routine  IB01MD.  For the N4SID approach, a weighted */
/*     oblique projection is computed from the upper triangular factor  R */
/*     and its SVD is then found. */

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

/*     [5] Van Overschee, P., and De Moor, B. */
/*         Subspace Identification for Linear Systems: Theory - */
/*         Implementation - Applications. */
/*         Kluwer Academic Publishers, Boston/London/Dordrecht, 1996. */

/*     [6] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */
/*                                      3 */
/*     The algorithm requires 0(((m+l)s) ) floating point operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     Feb. 2000, Feb. 2001, Feb. 2004, March 2005. */

/*     KEYWORDS */

/*     Identification methods, multivariable systems, QR decomposition, */
/*     singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --sv;
    --iwork;
    --dwork;

    /* Function Body */
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
    jobdm = lsame_(jobd, "M", (ftnlen)1, (ftnlen)1);
    mnobr = *m * *nobr;
    lnobr = *l * *nobr;
    llnobr = lnobr + lnobr;
    lmnobr = lnobr + mnobr;
    mmnobr = mnobr + mnobr;
    lmmnob = mmnobr + lnobr;
    nr = lmnobr + lmnobr;
    *iwarn = 0;
    *info = 0;

/*     Check the scalar input parameters. */

    if (! (moesp || n4sid)) {
	*info = -1;
    } else if (moesp && ! (jobdm || lsame_(jobd, "N", (ftnlen)1, (ftnlen)1))) 
	    {
	*info = -2;
    } else if (*nobr <= 0) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*l <= 0) {
	*info = -5;
    } else if (*ldr < nr || moesp && jobdm && *ldr < mnobr * 3) {
	*info = -7;
    } else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance. */
/*         NB refers to the optimal block size for the immediately */
/*         following subroutine, as returned by ILAENV.) */

	minwrk = 1;
	if (*ldwork >= 1) {
	    if (moesp) {
		minwrk = lnobr * 5;
		if (jobdm) {
/* Computing MAX */
		    i__1 = mmnobr - *nobr, i__1 = max(i__1,lmnobr);
		    minwrk = max(i__1,minwrk);
		}
		maxwrk = lnobr + lnobr * ilaenv_(&c__1, "DGEQRF", " ", &
			lmnobr, &lnobr, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
	    } else {

/* Computing MAX */
		i__1 = minwrk, i__2 = lmnobr * 5 + 1;
		minwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = mnobr + mnobr * ilaenv_(&c__1, "DGEQRF", " ", &mmnobr, 
			&mnobr, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1), i__2 = 
			mnobr + llnobr * ilaenv_(&c__1, "DORMQR", "LT", &
			mmnobr, &llnobr, &mnobr, &c_n1, (ftnlen)6, (ftnlen)2);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = mnobr + lnobr * ilaenv_(&c__1, "DORMQR",
			 "LN", &mmnobr, &lnobr, &mnobr, &c_n1, (ftnlen)6, (
			ftnlen)2);
		maxwrk = max(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = lnobr + lnobr * ilaenv_(&c__1, "DGEQRF",
			 " ", &lmmnob, &lnobr, &c_n1, &c_n1, (ftnlen)6, (
			ftnlen)1);
		maxwrk = max(i__1,i__2);
	    }
	    maxwrk = max(minwrk,maxwrk);
	}

	if (*ldwork < minwrk) {
	    *info = -12;
	    dwork[1] = (doublereal) minwrk;
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01ND", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compute pointers to the needed blocks of  R. */

    nr2 = mnobr + 1;
    nr3 = mmnobr + 1;
    nr4 = lmmnob + 1;
    itau = 1;
    jwork = itau + mnobr;

    if (moesp) {

/*        MOESP approach. */

	if (*m > 0 && jobdm) {

/*           Rearrange the blocks of  R: */
/*           Copy the (1,1) block into the position (3,2) and */
/*           copy the (1,4) block into (3,3). */

	    dlacpy_("Upper", &mnobr, &mnobr, &r__[r_offset], ldr, &r__[nr3 + 
		    nr2 * r_dim1], ldr, (ftnlen)5);
	    dlacpy_("Full", &mnobr, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		    nr3 + nr3 * r_dim1], ldr, (ftnlen)4);

/*           Using structure, triangularize the matrix */
/*              R_1c = [ R_12'  R_22'  R_11' ]' */
/*           and then apply the transformations to the matrix */
/*              R_2c = [ R_13'  R_23'  R_14' ]'. */
/*           Workspace: need M*NOBR + MAX(M-1,L)*NOBR. */

	    mb04od_("Upper", &mnobr, &lnobr, &mnobr, &r__[nr2 + nr2 * r_dim1],
		     ldr, &r__[nr3 + nr2 * r_dim1], ldr, &r__[nr2 + nr3 * 
		    r_dim1], ldr, &r__[nr3 + nr3 * r_dim1], ldr, &dwork[itau],
		     &dwork[jwork], (ftnlen)5);
	    i__1 = mnobr - 1;
	    i__2 = *ldwork - jwork + 1;
	    mb04id_(&mmnobr, &mnobr, &i__1, &lnobr, &r__[nr2 * r_dim1 + 1], 
		    ldr, &r__[nr3 * r_dim1 + 1], ldr, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Copy the leading  M*NOBR x M*NOBR  and  M*NOBR x L*NOBR */
/*           submatrices of  R_1c  and  R_2c,  respectively, into their */
/*           final positions, required by SLICOT Library routine  IB01PD. */

	    dlacpy_("Upper", &mnobr, &mnobr, &r__[nr2 * r_dim1 + 1], ldr, &
		    r__[lmnobr + 1 + r_dim1], ldr, (ftnlen)5);
	    dlacpy_("Full", &mnobr, &lnobr, &r__[nr3 * r_dim1 + 1], ldr, &r__[
		    nr2 * r_dim1 + 1], ldr, (ftnlen)4);
	}

/*        Copy [ R_24'  R_34' ]'  in  [ R_22'  R_32' ]'. */

	dlacpy_("Full", &lmnobr, &lnobr, &r__[nr2 + nr4 * r_dim1], ldr, &r__[
		nr2 + nr2 * r_dim1], ldr, (ftnlen)4);

/*        Triangularize the matrix in  [ R_22'  R_32' ]'. */
/*        Workspace: need 2*L*NOBR; prefer L*NOBR + L*NOBR*NB. */

	jwork = itau + lnobr;
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(&lmnobr, &lnobr, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[itau], 
		&dwork[jwork], &i__1, &ierr);

    } else {

/*        N4SID approach. */

	dum[0] = 0.;
	llmnob = llnobr + mnobr;

/*        Set the precision parameters. A threshold value  EPS**(2/3)  is */
/*        used for deciding to use pivoting or not, where  EPS  is the */
/*        relative machine precision (see LAPACK Library routine DLAMCH). */

	toll = *tol;
	eps = dlamch_("Precision", (ftnlen)9);
	thresh = pow_dd(&eps, &c_b37);

	if (*m > 0) {

/*           For efficiency of later calculations, interchange the first */
/*           two block-columns. The corresponding submatrices are */
/*           redefined according to their new position. */

	    i__1 = mnobr;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dswap_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[(mnobr + i__)
			 * r_dim1 + 1], &c__1);
		dcopy_(&mnobr, &r__[i__ + 1 + (mnobr + i__) * r_dim1], &c__1, 
			&r__[i__ + 1 + i__ * r_dim1], &c__1);
		i__2 = mmnobr - i__;
		dcopy_(&i__2, dum, &c__0, &r__[i__ + 1 + (mnobr + i__) * 
			r_dim1], &c__1);
/* L10: */
	    }

/*           Now, */

/*           U_f = [ R_11'  R_21'    0      0   ]', */
/*           U_p = [ R_12'    0      0      0   ]', */
/*           Y_p = [ R_13'  R_23'  R_33'    0   ]',  and */
/*           Y_f = [ R_14'  R_24'  R_34'  R_44' ]', */

/*           where  R_21,  R_12,  R_33,  and  R_44  are upper triangular. */
/*           Define  W_p := [ U_p  Y_p ]. */

/*           Prepare the computation of residuals of the two least */
/*           squares problems giving the weighted oblique projection P: */

/*           r_1 = W_p - U_f X_1,   X_1 = arg min || U_f X - W_p ||, */
/*           r_2 = Y_f - U_f X_2,   X_2 = arg min || U_f X - Y_f ||, */

/*           P = (arg min || r_1 X - r_2 ||)' r_1'.                   (1) */

/*           Alternately,  P'  is given by the projection */
/*              P' = Q_1 (Q_1)' r_2, */
/*           where  Q_1  contains the first  k  columns of the orthogonal */
/*           matrix in the  QR  factorization of  r_1,  k := rank(r_1). */

/*           Triangularize the matrix  U_f = q r  (using structure), and */
/*           apply the transformation  q'  to the corresponding part of */
/*           the matrices  W_p,  and  Y_f. */
/*           Workspace: need 2*(M+L)*NOBR. */

	    i__1 = mnobr - 1;
	    i__2 = *ldwork - jwork + 1;
	    mb04id_(&mmnobr, &mnobr, &i__1, &llmnob, &r__[r_offset], ldr, &
		    r__[nr2 * r_dim1 + 1], ldr, &dwork[itau], &dwork[jwork], &
		    i__2, &ierr);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Save updated  Y_f  (transposed) in the last block-row of  R. */

	    ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &
		    r__[nr4 + r_dim1], ldr, (ftnlen)4);

/*           Check the condition of the triangular factor  r  and decide */
/*           to use pivoting or not. */
/*           Workspace: need 4*M*NOBR. */

	    dtrcon_("1-norm", "Upper", "NonUnit", &mnobr, &r__[r_offset], ldr,
		     &rcond1, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, (
		    ftnlen)5, (ftnlen)7);

	    if (toll <= 0.) {
		toll = mnobr * mnobr * eps;
	    }
	    if (rcond1 > max(toll,thresh)) {

/*              U_f is considered full rank and no pivoting is used. */

		dlaset_("Full", &mnobr, &llmnob, &c_b50, &c_b50, &r__[nr2 * 
			r_dim1 + 1], ldr, (ftnlen)4);
	    } else {

/*              Save information about  q  in the (2,1) block of  R. */
/*              Use QR factorization with column pivoting,  r P = Q R. */
/*              Information on  Q  is stored in the strict lower triangle */
/*              of R_11  and in  DWORK(ITAU2). */

		i__1 = mnobr - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = nr2;
		    for (j = mmnobr; j >= i__2; --j) {
			r__[j + i__ * r_dim1] = r__[j - mnobr + i__ + i__ * 
				r_dim1];
/* L15: */
		    }
		    i__2 = mnobr - i__;
		    dcopy_(&i__2, dum, &c__0, &r__[i__ + 1 + i__ * r_dim1], &
			    c__1);
		    iwork[i__] = 0;
/* L20: */
		}

		iwork[mnobr] = 0;

/*              Workspace: need   5*M*NOBR+1. */
/*                         prefer 4*M*NOBR + (M*NOBR+1)*NB. */

		itau2 = jwork;
		jwork = itau2 + mnobr;
		svlmax = 0.;
		i__1 = *ldwork - jwork + 1;
		mb03od_("QR", &mnobr, &mnobr, &r__[r_offset], ldr, &iwork[1], 
			&toll, &svlmax, &dwork[itau2], &rank, sval, &dwork[
			jwork], &i__1, &ierr, (ftnlen)2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
		maxwrk = max(i__1,i__2);

/*              Workspace: need   2*M*NOBR + (M+2*L)*NOBR; */
/*                         prefer 2*M*NOBR + (M+2*L)*NOBR*NB. */

		i__1 = *ldwork - jwork + 1;
		dormqr_("Left", "Transpose", &mnobr, &llmnob, &mnobr, &r__[
			r_offset], ldr, &dwork[itau2], &r__[nr2 * r_dim1 + 1],
			 ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)
			9);
/* Computing MAX */
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
		maxwrk = max(i__1,i__2);
		if (rank < mnobr) {

/*                 The least squares problem is rank-deficient. */

		    *iwarn = 4;
		}

/*              Determine residuals r_1 and r_2: premultiply by  Q  and */
/*              then by  q. */
/*              Workspace: need   2*M*NOBR + (M+2*L)*NOBR); */
/*                         prefer 2*M*NOBR + (M+2*L)*NOBR*NB. */

		dlaset_("Full", &rank, &llmnob, &c_b50, &c_b50, &r__[nr2 * 
			r_dim1 + 1], ldr, (ftnlen)4);
		i__1 = *ldwork - jwork + 1;
		dormqr_("Left", "NoTranspose", &mnobr, &llmnob, &mnobr, &r__[
			r_offset], ldr, &dwork[itau2], &r__[nr2 * r_dim1 + 1],
			 ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)
			11);
/* Computing MAX */
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
		maxwrk = max(i__1,i__2);
		jwork = itau2;

/*              Restore the transformation  q. */

		i__1 = mnobr - 1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = mmnobr;
		    for (j = nr2; j <= i__2; ++j) {
			r__[j - mnobr + i__ + i__ * r_dim1] = r__[j + i__ * 
				r_dim1];
/* L25: */
		    }
/* L30: */
		}

	    }

/*           Premultiply by the transformation  q  (apply transformations */
/*           in backward order). */
/*           Workspace: need   M*NOBR + (M+2*L)*NOBR; */
/*                      prefer larger. */

	    i__1 = mnobr - 1;
	    i__2 = *ldwork - jwork + 1;
	    mb04iy_("Left", "NoTranspose", &mmnobr, &llmnob, &mnobr, &i__1, &
		    r__[r_offset], ldr, &dwork[itau], &r__[nr2 * r_dim1 + 1], 
		    ldr, &dwork[jwork], &i__2, &ierr, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

	} else {

/*           Save  Y_f  (transposed) in the last block-row of  R. */

	    ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &
		    r__[nr4 + r_dim1], ldr, (ftnlen)4);
	    rcond1 = 1.;
	}

/*        Triangularize the matrix  r_1  for determining the oblique */
/*        projection  P  in least squares problem in (1).  Exploit the */
/*        fact that the third block-row of r_1  has the structure */
/*        [ 0  T ],  where  T  is an upper triangular matrix.  Then apply */
/*        the corresponding transformations  Q'  to the matrix  r_2. */
/*        Workspace: need   2*M*NOBR; */
/*                   prefer   M*NOBR + M*NOBR*NB. */

	i__1 = *ldwork - jwork + 1;
	dgeqrf_(&mmnobr, &mnobr, &r__[nr2 * r_dim1 + 1], ldr, &dwork[itau], &
		dwork[jwork], &i__1, &ierr);

/*        Workspace: need   M*NOBR + 2*L*NOBR; */
/*                   prefer M*NOBR + 2*L*NOBR*NB. */

	i__1 = *ldwork - jwork + 1;
	dormqr_("Left", "Transpose", &mmnobr, &llnobr, &mnobr, &r__[nr2 * 
		r_dim1 + 1], ldr, &dwork[itau], &r__[nr3 * r_dim1 + 1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
	nrsave = nr2;

	itau2 = jwork;
	jwork = itau2 + lnobr;
	i__1 = lnobr - 1;
	i__2 = *ldwork - jwork + 1;
	mb04id_(&lmnobr, &lnobr, &i__1, &lnobr, &r__[nr2 + nr3 * r_dim1], ldr,
		 &r__[nr2 + nr4 * r_dim1], ldr, &dwork[itau2], &dwork[jwork], 
		&i__2, &ierr);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__2);

/*        Check the condition of the triangular matrix of order  (m+l)*s */
/*        just determined, and decide to use pivoting or not. */
/*        Workspace: need 4*(M+L)*NOBR. */

	dtrcon_("1-norm", "Upper", "NonUnit", &lmnobr, &r__[nr2 * r_dim1 + 1],
		 ldr, &rcond2, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, (
		ftnlen)5, (ftnlen)7);

	if (*tol <= 0.) {
	    toll = lmnobr * lmnobr * eps;
	}
	if (rcond2 <= max(toll,thresh)) {
	    if (*m > 0) {

/*              Save information about  Q  in  R_11  (in the strict lower */
/*              triangle),  R_21  and  R_31  (transposed information). */

		i__1 = mmnobr - 1;
		dlacpy_("Lower", &i__1, &mnobr, &r__[nr2 * r_dim1 + 2], ldr, &
			r__[r_dim1 + 2], ldr, (ftnlen)5);
		nrsave = 1;

		i__1 = lmnobr;
		for (i__ = nr2; i__ <= i__1; ++i__) {
		    dcopy_(&mnobr, &r__[i__ + 1 + (mnobr + i__) * r_dim1], &
			    c__1, &r__[mnobr + i__ + r_dim1], ldr);
/* L40: */
		}

	    }

	    i__1 = lmnobr - 1;
	    i__2 = lmnobr - 1;
	    dlaset_("Lower", &i__1, &i__2, &c_b50, &c_b50, &r__[nr2 * r_dim1 
		    + 2], ldr, (ftnlen)5);

/*           Use QR factorization with column pivoting. */
/*           Workspace: need   5*(M+L)*NOBR+1. */
/*                      prefer 4*(M+L)*NOBR + ((M+L)*NOBR+1)*NB. */

	    i__1 = lmnobr;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = 0;
/* L50: */
	    }

	    itau3 = jwork;
	    jwork = itau3 + lmnobr;
	    svlmax = 0.;
	    i__1 = *ldwork - jwork + 1;
	    mb03od_("QR", &lmnobr, &lmnobr, &r__[nr2 * r_dim1 + 1], ldr, &
		    iwork[1], &toll, &svlmax, &dwork[itau3], &rank1, sval, &
		    dwork[jwork], &i__1, &ierr, (ftnlen)2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);

/*           Workspace: need   2*(M+L)*NOBR + L*NOBR; */
/*                      prefer 2*(M+L)*NOBR + L*NOBR*NB. */

	    i__1 = *ldwork - jwork + 1;
	    dormqr_("Left", "Transpose", &lmnobr, &lnobr, &lmnobr, &r__[nr2 * 
		    r_dim1 + 1], ldr, &dwork[itau3], &r__[nr4 * r_dim1 + 1], 
		    ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	    if (rank1 < lmnobr) {

/*              The least squares problem is rank-deficient. */

		*iwarn = 5;
	    }

/*           Apply the orthogonal transformations, in backward order, to */
/*           [r_2(1:rank(r_1),:)' 0]',  to obtain  P'. */
/*           Workspace: need   2*(M+L)*NOBR + L*NOBR; */
/*                      prefer 2*(M+L)*NOBR + L*NOBR*NB. */

	    i__1 = lmnobr - rank1;
	    dlaset_("Full", &i__1, &lnobr, &c_b50, &c_b50, &r__[rank1 + 1 + 
		    nr4 * r_dim1], ldr, (ftnlen)4);
	    i__1 = *ldwork - jwork + 1;
	    dormqr_("Left", "NoTranspose", &lmnobr, &lnobr, &lmnobr, &r__[nr2 
		    * r_dim1 + 1], ldr, &dwork[itau3], &r__[nr4 * r_dim1 + 1],
		     ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	    maxwrk = max(i__1,i__2);
	    jwork = itau3;

	    if (*m > 0) {

/*              Restore the saved transpose matrix from  R_31. */

		i__1 = lmnobr;
		for (i__ = nr2; i__ <= i__1; ++i__) {
		    dcopy_(&mnobr, &r__[mnobr + i__ + r_dim1], ldr, &r__[i__ 
			    + 1 + (mnobr + i__) * r_dim1], &c__1);
/* L60: */
		}

	    }

	}

/*        Workspace: need   M*NOBR + L*NOBR; */
/*                   prefer larger. */

	i__1 = lnobr - 1;
	i__2 = *ldwork - jwork + 1;
	mb04iy_("Left", "NoTranspose", &lmnobr, &lnobr, &lnobr, &i__1, &r__[
		nr2 + nr3 * r_dim1], ldr, &dwork[itau2], &r__[nr2 + nr4 * 
		r_dim1], ldr, &dwork[jwork], &i__2, &ierr, (ftnlen)4, (ftnlen)
		11);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
	maxwrk = max(i__1,i__2);

/*        Workspace: need   M*NOBR + L*NOBR; */
/*                   prefer M*NOBR + L*NOBR*NB. */

	jwork = itau2;
	i__1 = *ldwork - jwork + 1;
	dormqr_("Left", "NoTranspose", &mmnobr, &lnobr, &mnobr, &r__[nrsave * 
		r_dim1 + 1], ldr, &dwork[itau], &r__[nr4 * r_dim1 + 1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)11);

/*        Now, the matrix  P'  is available in  R_14 : R_34. */
/*        Triangularize the matrix  P'. */
/*        Workspace: need   2*L*NOBR; */
/*                   prefer   L*NOBR + L*NOBR*NB. */

	jwork = itau + lnobr;
	i__1 = *ldwork - jwork + 1;
	dgeqrf_(&lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &dwork[itau], &
		dwork[jwork], &i__1, &ierr);

/*        Copy the triangular factor to its final position,  R_22. */

	dlacpy_("Upper", &lnobr, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		nr2 + nr2 * r_dim1], ldr, (ftnlen)5);

/*        Restore  Y_f. */

	ma02ad_("Full", &lnobr, &lmmnob, &r__[nr4 + r_dim1], ldr, &r__[nr4 * 
		r_dim1 + 1], ldr, (ftnlen)4);
    }

/*     Find the singular value decomposition of  R_22. */
/*     Workspace: need 5*L*NOBR. */

    mb03ud_("NoVectors", "Vectors", &lnobr, &r__[nr2 + nr2 * r_dim1], ldr, 
	    dum, &c__1, &sv[1], &dwork[1], ldwork, &ierr, (ftnlen)9, (ftnlen)
	    7);
    if (ierr != 0) {
	*info = 2;
	return 0;
    }
/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) dwork[1];
    maxwrk = max(i__1,i__2);

/*     Transpose  R(m*s+1:(m+L)*s,m*s+1:(m+L)*s)  in-situ; its */
/*     columns will then be the singular vectors needed subsequently. */

    i__1 = lmnobr;
    for (i__ = nr2 + 1; i__ <= i__1; ++i__) {
	i__2 = lmnobr - i__ + 1;
	dswap_(&i__2, &r__[i__ + (i__ - 1) * r_dim1], &c__1, &r__[i__ - 1 + 
		i__ * r_dim1], ldr);
/* L70: */
    }

/*     Return optimal workspace in  DWORK(1)  and reciprocal condition */
/*     numbers, if  METH = 'N'. */

    dwork[1] = (doublereal) maxwrk;
    if (n4sid) {
	dwork[2] = rcond1;
	dwork[3] = rcond2;
    }
    return 0;

/* *** Last line of IB01ND *** */
} /* ib01nd_ */

