/* IB01AD.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int ib01ad_(char *meth, char *alg, char *jobd, char *batch, 
	char *conct, char *ctrl, integer *nobr, integer *m, integer *l, 
	integer *nsmp, doublereal *u, integer *ldu, doublereal *y, integer *
	ldy, integer *n, doublereal *r__, integer *ldr, doublereal *sv, 
	doublereal *rcond, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, integer *iwarn, integer *info, ftnlen meth_len, 
	ftnlen alg_len, ftnlen jobd_len, ftnlen batch_len, ftnlen conct_len, 
	ftnlen ctrl_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer nr, ns;
    static logical last, n4sid;
    extern /* Subroutine */ int ib01md_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), ib01nd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen);
    static logical chalg;
    extern /* Subroutine */ int ib01od_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    static logical jobdm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nobr21;
    static logical qralg;
    static integer lnobr, mnobr;
    static logical moesp, first, onebch, connec, fqralg;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lmnobr, iwarnl;
    static logical interm, contrl;
    static integer minwrk, maxwrk, nsmpsm;


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

/*     To preprocess the input-output data for estimating the matrices */
/*     of a linear time-invariant dynamical system and to find an */
/*     estimate of the system order. The input-output data can, */
/*     optionally, be processed sequentially. */

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

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not the matrices B and D should later */
/*             be computed using the MOESP approach, as follows: */
/*             = 'M':  the matrices B and D should later be computed */
/*                     using the MOESP approach; */
/*             = 'N':  the matrices B and D should not be computed using */
/*                     the MOESP approach. */
/*             This parameter is not relevant for METH = 'N'. */

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

/*     CTRL    CHARACTER*1 */
/*             Specifies whether or not the user's confirmation of the */
/*             system order estimate is desired, as follows: */
/*             = 'C':  user's confirmation; */
/*             = 'N':  no confirmation. */
/*             If  CTRL = 'C',  a reverse communication routine,  IB01OY, */
/*             is indirectly called (by SLICOT Library routine IB01OD), */
/*             and, after inspecting the singular values and system order */
/*             estimate,  n,  the user may accept  n  or set a new value. */
/*             IB01OY  is not called if CTRL = 'N'. */

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

/*     N       (output) INTEGER */
/*             The estimated order of the system. */
/*             If  CTRL = 'C',  the estimated order has been reset to a */
/*             value specified by the user. */

/*     R       (output or input/output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On exit, if ALG = 'C' and BATCH = 'F' or 'I', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the current upper triangular part of the */
/*             correlation matrix in sequential data processing. */
/*             If ALG = 'F' and BATCH = 'F' or 'I', the array R is not */
/*             referenced. */
/*             On exit, if INFO = 0, ALG = 'Q', and BATCH = 'F' or 'I', */
/*             the leading 2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular */
/*             part of this array contains the current upper triangular */
/*             factor R from the QR factorization of the concatenated */
/*             block Hankel matrices. Denote  R_ij, i,j = 1:4,  the */
/*             ij submatrix of  R,  partitioned by M*NOBR,  M*NOBR, */
/*             L*NOBR,  and  L*NOBR  rows and columns. */
/*             On exit, if INFO = 0 and BATCH = 'L' or 'O', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of */
/*             this array contains the matrix S, the processed upper */
/*             triangular factor R from the QR factorization of the */
/*             concatenated block Hankel matrices, as required by other */
/*             subroutines. Specifically, let  S_ij, i,j = 1:4,  be the */
/*             ij submatrix of  S,  partitioned by M*NOBR,  L*NOBR, */
/*             M*NOBR,  and  L*NOBR  rows and columns. The submatrix */
/*             S_22  contains the matrix of left singular vectors needed */
/*             subsequently. Useful information is stored in  S_11  and */
/*             in the block-column  S_14 : S_44.  For METH = 'M' and */
/*             JOBD = 'M', the upper triangular part of  S_31  contains */
/*             the upper triangular factor in the QR factorization of the */
/*             matrix  R_1c = [ R_12'  R_22'  R_11' ]',  and  S_12 */
/*             contains the corresponding leading part of the transformed */
/*             matrix  R_2c = [ R_13'  R_23'  R_14' ]'.  For  METH = 'N', */
/*             the subarray  S_41 : S_43  contains the transpose of the */
/*             matrix contained in  S_14 : S_34. */
/*             The details of the contents of R need not be known if this */
/*             routine is followed by SLICOT Library routine IB01BD. */
/*             On entry, if ALG = 'C', or ALG = 'Q', and BATCH = 'I' or */
/*             'L', the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  upper */
/*             triangular part of this array must contain the upper */
/*             triangular matrix R computed at the previous call of this */
/*             routine in sequential data processing. The array R need */
/*             not be set on entry if ALG = 'F' or if BATCH = 'F' or 'O'. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= MAX( 2*(M+L)*NOBR, 3*M*NOBR ), */
/*                                  for METH = 'M' and JOBD = 'M'; */
/*             LDR >= 2*(M+L)*NOBR, for METH = 'M' and JOBD = 'N' or */
/*                                  for METH = 'N'. */

/*     SV      (output) DOUBLE PRECISION array, dimension ( L*NOBR ) */
/*             The singular values used to estimate the system order. */

/*     Tolerances */

/*     RCOND   DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  RCOND > 0,  the given value */
/*             of  RCOND  is used as a lower bound for the reciprocal */
/*             condition number;  an m-by-n matrix whose estimated */
/*             condition number is less than  1/RCOND  is considered to */
/*             be of full rank.  If the user sets  RCOND <= 0,  then an */
/*             implicitly computed, default tolerance, defined by */
/*             RCONDEF = m*n*EPS,  is used instead, where  EPS  is the */
/*             relative machine precision (see LAPACK Library routine */
/*             DLAMCH). */
/*             This parameter is not used for  METH = 'M'. */

/*     TOL     DOUBLE PRECISION */
/*             Absolute tolerance used for determining an estimate of */
/*             the system order. If  TOL >= 0,  the estimate is */
/*             indicated by the index of the last singular value greater */
/*             than or equal to  TOL.  (Singular values less than  TOL */
/*             are considered as zero.) When  TOL = 0,  an internally */
/*             computed default value,  TOL = NOBR*EPS*SV(1),  is used, */
/*             where  SV(1)  is the maximal singular value, and  EPS  is */
/*             the relative machine precision (see LAPACK Library routine */
/*             DLAMCH). When  TOL < 0,  the estimate is indicated by the */
/*             index of the singular value that has the largest */
/*             logarithmic gap to its successor. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= (M+L)*NOBR, if METH = 'N'; */
/*             LIWORK >= M+L, if METH = 'M' and ALG = 'F'; */
/*             LIWORK >= 0,   if METH = 'M' and ALG = 'C' or 'Q'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and, for  METH = 'N',  and  BATCH = 'L'  or */
/*             'O',  DWORK(2)  and  DWORK(3)  contain the reciprocal */
/*             condition numbers of the triangular factors of the */
/*             matrices  U_f  and  r_1  [6]. */
/*             On exit, if  INFO = -23,  DWORK(1)  returns the minimum */
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
/*             LDWORK >= (4*NOBR-2)*(M+L), if ALG = 'C', BATCH = 'F' or */
/*                             'I' and CONCT = 'C'; */
/*             LDWORK >= 1, if ALG = 'C', BATCH = 'F' or 'I' and */
/*                             CONCT = 'N'; */
/*             LDWORK >= max((4*NOBR-2)*(M+L), 5*L*NOBR), if METH = 'M', */
/*                             ALG = 'C', BATCH = 'L' and CONCT = 'C'; */
/*             LDWORK >= max((2*M-1)*NOBR, (M+L)*NOBR, 5*L*NOBR), */
/*                             if METH = 'M', JOBD = 'M', ALG = 'C', */
/*                              BATCH = 'O', or */
/*                             (BATCH = 'L' and CONCT = 'N'); */
/*             LDWORK >= 5*L*NOBR, if METH = 'M', JOBD = 'N', ALG = 'C', */
/*                              BATCH = 'O', or */
/*                             (BATCH = 'L' and CONCT = 'N'); */
/*             LDWORK >= 5*(M+L)*NOBR+1, if METH = 'N', ALG = 'C', and */
/*                             BATCH = 'L' or 'O'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+3), if ALG = 'F', */
/*                             BATCH <> 'O' and CONCT = 'C'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+1), if ALG = 'F', */
/*                             BATCH = 'F', 'I' and CONCT = 'N'; */
/*             LDWORK >= (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR, if ALG = 'F', */
/*                             BATCH = 'L' and CONCT = 'N', or */
/*                             BATCH = 'O'; */
/*             LDWORK >= 4*(M+L)*NOBR, if ALG = 'Q', BATCH = 'F', and */
/*                             LDR >= NS = NSMP - 2*NOBR + 1; */
/*             LDWORK >= max(4*(M+L)*NOBR, 5*L*NOBR), if METH = 'M', */
/*                             ALG = 'Q', BATCH = 'O', and LDR >= NS; */
/*             LDWORK >= 5*(M+L)*NOBR+1, if METH = 'N', ALG = 'Q', */
/*                             BATCH = 'O', and LDR >= NS; */
/*             LDWORK >= 6*(M+L)*NOBR, if ALG = 'Q', (BATCH = 'F' or 'O', */
/*                             and LDR < NS), or (BATCH = 'I' or */
/*                             'L' and CONCT = 'N'); */
/*             LDWORK >= 4*(NOBR+1)*(M+L)*NOBR, if ALG = 'Q', BATCH = 'I' */
/*                             or 'L' and CONCT = 'C'. */
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
/*                   (non-sequential data processing); */
/*             = 3:  all singular values were exactly zero, hence  N = 0 */
/*                   (both input and output were identically zero); */
/*             = 4:  the least squares problems with coefficient matrix */
/*                   U_f,  used for computing the weighted oblique */
/*                   projection (for METH = 'N'), have a rank-deficient */
/*                   coefficient matrix; */
/*             = 5:  the least squares problem with coefficient matrix */
/*                   r_1  [6], used for computing the weighted oblique */
/*                   projection (for METH = 'N'), has a rank-deficient */
/*                   coefficient matrix. */
/*             NOTE: the values 4 and 5 of IWARN have no significance */
/*                   for the identification problem. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  a fast algorithm was requested (ALG = 'C', or 'F') */
/*                   in sequential data processing, but it failed; the */
/*                   routine can be repeatedly called again using the */
/*                   standard QR algorithm; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     The procedure consists in three main steps, the first step being */
/*     performed by one of the three algorithms included. */

/*     1.a) For non-sequential data processing using QR algorithm, a */
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

/*     1.b) For sequential data processing using QR algorithm, the QR */
/*     decomposition is done sequentially, by updating the upper */
/*     triangular factor  R.  This is also performed internally if the */
/*     workspace is not large enough to accommodate an entire batch. */

/*     1.c) For non-sequential or sequential data processing using */
/*     Cholesky algorithm, the correlation matrix of input-output data is */
/*     computed (sequentially, if requested), taking advantage of the */
/*     block Hankel structure [7].  Then, the Cholesky factor of the */
/*     correlation matrix is found, if possible. */

/*     2) A singular value decomposition (SVD) of a certain matrix is */
/*     then computed, which reveals the order  n  of the system as the */
/*     number of "non-zero" singular values. For the MOESP approach, this */
/*     matrix is  [ R_24'  R_34' ]' := R(ms+1:(2m+l)s,(2m+l)s+1:2(m+l)s), */
/*     where  R  is the upper triangular factor  R  constructed by SLICOT */
/*     Library routine  IB01MD.  For the N4SID approach, a weighted */
/*     oblique projection is computed from the upper triangular factor  R */
/*     and its SVD is then found. */

/*     3) The singular values are compared to the given, or default TOL, */
/*     and the estimated order  n  is returned, possibly after user's */
/*     confirmation. */

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
/*     The most time-consuming computational step is step 1: */
/*                                        2 */
/*     The QR algorithm needs 0(t(2(m+l)s) ) floating point operations. */
/*                                           2              3 */
/*     The Cholesky algorithm needs 0(2t(m+l) s)+0((2(m+l)s) ) floating */
/*     point operations. */
/*                                          2           3 2 */
/*     The fast QR algorithm needs 0(2t(m+l) s)+0(4(m+l) s ) floating */
/*     point operations. */
/*                                                3 */
/*     Step 2 of the algorithm requires 0(((m+l)s) ) floating point */
/*     operations. */

/*     FURTHER COMMENTS */

/*     For ALG = 'Q', BATCH = 'O' and LDR < NS, or BATCH <> 'O', the */
/*     calculations could be rather inefficient if only minimal workspace */
/*     (see argument LDWORK) is provided. It is advisable to provide as */
/*     much workspace as possible. Almost optimal efficiency can be */
/*     obtained for  LDWORK = (NS+2)*(2*(M+L)*NOBR),  assuming that the */
/*     cache size is large enough to accommodate R, U, Y, and DWORK. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Universiteit Leuven, Feb. 2000. */

/*     REVISIONS */

/*     August 2000, March 2005. */

/*     KEYWORDS */

/*     Cholesky decomposition, Hankel matrix, identification methods, */
/*     multivariable systems, QR decomposition, singular value */
/*     decomposition. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save Statement .. */
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
    --sv;
    --iwork;
    --dwork;

    /* Function Body */
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
    fqralg = lsame_(alg, "F", (ftnlen)1, (ftnlen)1);
    qralg = lsame_(alg, "Q", (ftnlen)1, (ftnlen)1);
    chalg = lsame_(alg, "C", (ftnlen)1, (ftnlen)1);
    jobdm = lsame_(jobd, "M", (ftnlen)1, (ftnlen)1);
    onebch = lsame_(batch, "O", (ftnlen)1, (ftnlen)1);
    first = lsame_(batch, "F", (ftnlen)1, (ftnlen)1) || onebch;
    interm = lsame_(batch, "I", (ftnlen)1, (ftnlen)1);
    last = lsame_(batch, "L", (ftnlen)1, (ftnlen)1) || onebch;
    contrl = lsame_(ctrl, "C", (ftnlen)1, (ftnlen)1);

    if (! onebch) {
	connec = lsame_(conct, "C", (ftnlen)1, (ftnlen)1);
    } else {
	connec = FALSE_;
    }

    mnobr = *m * *nobr;
    lnobr = *l * *nobr;
    lmnobr = lnobr + mnobr;
    nr = lmnobr + lmnobr;
    nobr21 = (*nobr << 1) - 1;
    *iwarn = 0;
    *info = 0;
    if (first) {
	maxwrk = 1;
	nsmpsm = 0;
    }
    nsmpsm += *nsmp;

/*     Check the scalar input parameters. */

    if (! (moesp || n4sid)) {
	*info = -1;
    } else if (! (fqralg || qralg || chalg)) {
	*info = -2;
    } else if (moesp && ! (jobdm || lsame_(jobd, "N", (ftnlen)1, (ftnlen)1))) 
	    {
	*info = -3;
    } else if (! (first || interm || last)) {
	*info = -4;
    } else if (! onebch) {
	if (! (connec || lsame_(conct, "N", (ftnlen)1, (ftnlen)1))) {
	    *info = -5;
	}
    }
    if (*info == 0) {
	if (! (contrl || lsame_(ctrl, "N", (ftnlen)1, (ftnlen)1))) {
	    *info = -6;
	} else if (*nobr <= 0) {
	    *info = -7;
	} else if (*m < 0) {
	    *info = -8;
	} else if (*l <= 0) {
	    *info = -9;
	} else if (*nsmp < *nobr << 1 || last && nsmpsm < nr + nobr21) {
	    *info = -10;
	} else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
	    *info = -12;
	} else if (*ldy < *nsmp) {
	    *info = -14;
	} else if (*ldr < nr || moesp && jobdm && *ldr < mnobr * 3) {
	    *info = -17;
	} else {

/*           Compute workspace. */
/*           (Note: Comments in the code beginning "Workspace:" describe */
/*           the minimal amount of workspace needed at that point in the */
/*           code, as well as the preferred amount for good performance.) */

	    ns = *nsmp - nobr21;
	    if (chalg) {
		if (! last) {
		    if (connec) {
			minwrk = nr - *m - *l << 1;
		    } else {
			minwrk = 1;
		    }
		} else if (moesp) {
		    if (connec && ! onebch) {
/* Computing MAX */
			i__1 = nr - *m - *l << 1, i__2 = lnobr * 5;
			minwrk = max(i__1,i__2);
		    } else {
			minwrk = lnobr * 5;
			if (jobdm) {
/* Computing MAX */
			    i__1 = (mnobr << 1) - *nobr, i__1 = max(i__1,
				    lmnobr);
			    minwrk = max(i__1,minwrk);
			}
		    }
		} else {
		    minwrk = lmnobr * 5 + 1;
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
		if (onebch && *ldr >= ns) {
		    if (moesp) {
/* Computing MAX */
			i__1 = minwrk, i__2 = lnobr * 5;
			minwrk = max(i__1,i__2);
		    } else {
			minwrk = lmnobr * 5 + 1;
		    }
		}
		if (first) {
		    if (*ldr < ns) {
			minwrk += nr;
		    }
		} else {
		    if (connec) {
			minwrk *= *nobr + 1;
		    } else {
			minwrk += nr;
		    }
		}
	    }

	    maxwrk = minwrk;

	    if (*ldwork < minwrk) {
		*info = -23;
		dwork[1] = (doublereal) minwrk;
	    }
	}
    }

/*     Return if there are illegal arguments. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("IB01AD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compress the input-output data. */
/*     Workspace: need   c*(M+L)*NOBR, where c is a constant depending */
/*                       on the algorithm and the options used */
/*                       (see SLICOT Library routine IB01MD); */
/*                prefer larger. */

    ib01md_(meth, alg, batch, conct, nobr, m, l, nsmp, &u[u_offset], ldu, &y[
	    y_offset], ldy, &r__[r_offset], ldr, &iwork[1], &dwork[1], ldwork,
	     iwarn, info, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

    if (*info == 1) {

/*        Error return: A fast algorithm was requested (ALG = 'C', 'F') */
/*        in sequential data processing, but it failed. */

	return 0;
    }

/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) dwork[1];
    maxwrk = max(i__1,i__2);

    if (! last) {

/*        Return to get new data. */

	return 0;
    }

/*     Find the singular value decomposition (SVD) giving the system */
/*     order, and perform related preliminary calculations needed for */
/*     computing the system matrices. */
/*     Workspace: need   max( (2*M-1)*NOBR, (M+L)*NOBR, 5*L*NOBR ), */
/*                                            if METH = 'M'; */
/*                            5*(M+L)*NOBR+1, if METH = 'N'; */
/*                prefer larger. */

    ib01nd_(meth, jobd, nobr, m, l, &r__[r_offset], ldr, &sv[1], rcond, &
	    iwork[1], &dwork[1], ldwork, &iwarnl, info, (ftnlen)1, (ftnlen)1);
    *iwarn = max(*iwarn,iwarnl);

    if (*info == 2) {

/*        Error return: the singular value decomposition (SVD) algorithm */
/*        did not converge. */

	return 0;
    }

/*     Estimate the system order. */

    ib01od_(ctrl, nobr, l, &sv[1], n, tol, &iwarnl, info, (ftnlen)1);
    *iwarn = max(*iwarn,iwarnl);

/*     Return optimal workspace in  DWORK(1). */

/* Computing MAX */
    i__1 = maxwrk, i__2 = (integer) dwork[1];
    dwork[1] = (doublereal) max(i__1,i__2);
    return 0;

/* *** Last line of IB01AD *** */
} /* ib01ad_ */

