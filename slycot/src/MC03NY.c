/* MC03NY.f -- translated by f2c (version 20100827).
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

static doublereal c_b9 = -1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;
static integer c__0 = 0;
static doublereal c_b26 = 1.;

/* Subroutine */ int mc03ny_(integer *nblcks, integer *nra, integer *nca, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, integer *
	imuk, integer *inuk, doublereal *veps, integer *ldveps, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, veps_dim1, veps_offset, i__1, 
	    i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, ac1, ac2, ec1, ar1, er1, vc1, vc2, wc1, vr1, 
	    vr2, wr1, dif, ari, ark, ncv, mui, nui, nrv, smui, smui1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal dummy[1];
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrtrs_(char *, char *, char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);


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

/*     To determine a minimal basis of the right nullspace of the */
/*     subpencil s*E(eps)-A(eps) using the method given in [1] (see */
/*     Eqs.(4.6.8), (4.6.9)). */
/*     This pencil only contains Kronecker column indices, and it must be */
/*     in staircase form as supplied by SLICOT Library Routine MB04VD. */
/*     The basis vectors are represented by matrix V(s) having the form */

/*                | V11(s) V12(s) V13(s)   . .   V1n(s) | */
/*                |        V22(s) V23(s)         V2n(s) | */
/*                |               V33(s)           .    | */
/*         V(s) = |                  .             .    | */
/*                |                      .         .    | */
/*                |                          .     .    | */
/*                |                              Vnn(s) | */

/*     where n is the number of full row rank blocks in matrix A(eps) and */

/*                                               k               j-i */
/*         Vij(s) = Vij,0 + Vij,1*s +...+ Vij,k*s +...+ Vij,j-i*s   . (1) */

/*     In other words, Vij,k is the coefficient corresponding to degree k */
/*     in the matrix polynomial Vij(s). */
/*     Vij,k has dimensions mu(i)-by-(mu(j)-nu(j)). */
/*     The coefficients Vij,k are stored in the matrix VEPS as follows */
/*     (for the case n = 3): */

/*         sizes      m1-n1    m2-n2   m2-n2    m3-n3   m3-n3   m3-n3 */

/*             m1 { | V11,0 || V12,0 | V12,1 || V13,0 | V13,1 | V13,2 || */
/*                  |       ||       |       ||       |       |       || */
/*      VEPS = m2 { |       || V22,0 |       || V23,0 | V23,1 |       || */
/*                  |       ||       |       ||       |       |       || */
/*             m3 { |       ||       |       || V33,0 |       |       || */

/*     where mi = mu(i), ni = nu(i). */
/*     Matrix VEPS has dimensions nrv-by-ncv where */
/*       nrv = Sum(i=1,...,n) mu(i) */
/*       ncv = Sum(i=1,...,n) i*(mu(i)-nu(i)) */

/*     ================================================================== */
/*     REMARK: This routine is intended to be called only from the SLICOT */
/*             routine MC03ND. */
/*     ================================================================== */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NBLCKS  (input) INTEGER */
/*             Number of full row rank blocks in subpencil */
/*             s*E(eps)-A(eps) that contains all Kronecker column indices */
/*             of s*E-A.  NBLCKS >= 0. */

/*     NRA     (input) INTEGER */
/*             Number of rows of the subpencil s*E(eps)-A(eps) in s*E-A. */
/*             NRA = nu(1) + nu(2) + ... + nu(NBLCKS).  NRA >= 0. */

/*     NCA     (input) INTEGER */
/*             Number of columns of the subpencil s*E(eps)-A(eps) in */
/*             s*E-A. */
/*             NCA = mu(1) + mu(2) + ... + mu(NBLCKS).  NCA >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,NCA) */
/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,NCA) */
/*             On entry, the leading NRA-by-NCA part of these arrays must */
/*             contain the matrices A and E, where s*E-A is the */
/*             transformed pencil s*E0-A0 which is the pencil associated */
/*             with P(s) as described in [1] Section 4.6. The pencil */
/*             s*E-A is assumed to be in generalized Schur form. */
/*             On exit, these arrays contain no useful information. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,NRA). */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,NRA). */

/*     IMUK    (input) INTEGER array, dimension (NBLCKS) */
/*             This array must contain the column dimensions mu(k) of the */
/*             full column rank blocks in the subpencil s*E(eps)-A(eps) */
/*             of s*E-A. The content of IMUK is modified by the routine */
/*             but restored on exit. */

/*     INUK    (input) INTEGER array, dimension (NBLCKS) */
/*             This array must contain the row dimensions nu(k) of the */
/*             full row rank blocks in the subpencil s*E(eps)-A(eps) of */
/*             s*E-A. */

/*     VEPS    (output) DOUBLE PRECISION array, dimension (LDVEPS,ncv) */
/*             Let nrv = Sum(i=1,...,NBLCKS) mu(i) = NCA, */
/*                 ncv = Sum(i=1,...,NBLCKS) i*(mu(i)-nu(i)). */
/*             The leading nrv-by-ncv part of this array contains the */
/*             column vectors of a minimal polynomial basis for the right */
/*             nullspace of the subpencil s*E(eps)-A(eps). (See [1] */
/*             Section 4.6.4.) An upper bound for ncv is (NRA+1)*NCA. */

/*     LDVEPS  INTEGER */
/*             The leading dimension of array VEPS. */
/*             LDVEPS >= MAX(1,NCA). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = k, the k-th diagonal block of A had not a */
/*                   full row rank. */

/*     REFERENCES */

/*     [1] Th.G.J. Beelen, New Algorithms for Computing the Kronecker */
/*         structure of a Pencil with Applications to Systems and */
/*         Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, 1987. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC03BY by Th.G.J. Beelen, */
/*     A.J. Geurts, and G.J.H.H. van den Hurk. */

/*     REVISIONS */

/*     Dec. 1997. */

/*     KEYWORDS */

/*     Elementary polynomial operations, Kronecker form, polynomial */
/*     matrix, polynomial operations, staircase form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;
    --imuk;
    --inuk;
    veps_dim1 = *ldveps;
    veps_offset = 1 + veps_dim1;
    veps -= veps_offset;

    /* Function Body */
    *info = 0;
    if (*nblcks < 0) {
	*info = -1;
    } else if (*nra < 0) {
	*info = -2;
    } else if (*nca < 0) {
	*info = -3;
    } else if (*lda < max(1,*nra)) {
	*info = -5;
    } else if (*lde < max(1,*nra)) {
	*info = -7;
    } else if (*ldveps < max(1,*nca)) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla_("MC03NY", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible. */

    if (*nblcks == 0 || *nra == 0 || *nca == 0) {
	return 0;
    }

/*     Computation of the nonzero parts of W1 and W2: */

/*          | AH11 AH12 ... AH1n |       | EH11 EH12 ... EH1n | */
/*          |      AH22     AH2n |       |      EH22     EH2n | */
/*     W1 = |         .       .  |, W2 = |         .       .  | */
/*          |           .     .  |       |           .     .  | */
/*          |               AHnn |       |               EHnn | */

/*     with AHij = -pinv(Aii) * Aij, EHij = pinv(Aii) * Eij and EHii = 0, */
/*     AHij and EHij have dimensions mu(i)-by-mu(j), Aii = [ Oi | Ri ], */
/*     and */
/*       Ri is a regular nu(i)-by-nu(i) upper triangular matrix; */
/*       Oi is a not necessarily square null matrix. */
/*     Note that the first mu(i)-nu(i) rows in AHij and EHij are zero. */
/*     For memory savings, the nonzero parts of W1 and W2 are constructed */
/*     over A and E, respectively. */

/*     (AR1,AC1) denotes the position of the first element of the */
/*     submatrix Ri in matrix Aii. */
/*     EC1 is the index of the first column of Ai,i+1/Ei,i+1. */

    ec1 = 1;
    ar1 = 1;

    i__1 = *nblcks - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nui = inuk[i__];
	if (nui == 0) {
	    goto L60;
	}
	mui = imuk[i__];
	ec1 += mui;
	ac1 = ec1 - nui;
	i__2 = *nca - ec1 + 1;
	dtrtrs_("Upper", "No transpose", "Non-unit", &nui, &i__2, &a[ar1 + 
		ac1 * a_dim1], lda, &e[ar1 + ec1 * e_dim1], lde, info, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	if (*info > 0) {
	    *info = i__;
	    return 0;
	}

	i__2 = nui;
	for (j = 1; j <= i__2; ++j) {
	    dscal_(&j, &c_b9, &a[ar1 + (ac1 + j - 1) * a_dim1], &c__1);
/* L20: */
	}

	i__2 = *nca - ec1 + 1;
	dtrtrs_("Upper", "No transpose", "Non-unit", &nui, &i__2, &a[ar1 + 
		ac1 * a_dim1], lda, &a[ar1 + ec1 * a_dim1], lda, info, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
	ar1 += nui;
/* L40: */
    }

L60:

/*     The contents of the array IMUK is changed for temporary use in */
/*     this routine as follows: */

/*        IMUK(i) = Sum(j=1,...,i) mu(j). */

/*     On return, the original contents of IMUK is restored. */
/*     In the same loop the actual number of columns of VEPS is computed. */
/*     The number of rows of VEPS is NCA. */

/*        NRV = Sum(i=1,...,NBLCKS) mu(i) = NCA, */
/*        NCV = Sum(i=1,...,NBLCKS) i*(mu(i)-nu(i)). */

    smui = 0;
    ncv = 0;

    i__1 = *nblcks;
    for (i__ = 1; i__ <= i__1; ++i__) {
	mui = imuk[i__];
	smui += mui;
	imuk[i__] = smui;
	ncv += i__ * (mui - inuk[i__]);
/* L80: */
    }

    nrv = *nca;

/*     Computation of the matrix VEPS. */

/*     Initialisation of VEPS to zero. */

    dlaset_("Full", &nrv, &ncv, &c_b16, &c_b16, &veps[veps_offset], ldveps, (
	    ftnlen)4);
/*                                                           | I | */
/*     Set Vii,0 = Kii in VEPS , i=1,...,NBLCKS, where Kii = |---| */
/*                                                           | O | */
/*     and I is an identity matrix of size mu(i)-nu(i), */
/*         O is a null matrix, dimensions nu(i)-by-(mu(i)-nu(i)). */

/*     WR1 := Sum(j=1,...,i-1) mu(j) + 1 */
/*            is the index of the first row in Vii,0 in VEPS. */
/*     WC1 := Sum(j=1,...,i-1) j*(mu(j)-nu(j)) + 1 */
/*            is the index of the first column in Vii,0 in VEPS. */

    dummy[0] = 1.;
    nui = imuk[1] - inuk[1];
    i__1 = *ldveps + 1;
    dcopy_(&nui, dummy, &c__0, &veps[veps_offset], &i__1);
    wr1 = imuk[1] + 1;
    wc1 = nui + 1;

    i__1 = *nblcks;
    for (i__ = 2; i__ <= i__1; ++i__) {
	nui = imuk[i__] - imuk[i__ - 1] - inuk[i__];
	i__2 = *ldveps + 1;
	dcopy_(&nui, dummy, &c__0, &veps[wr1 + wc1 * veps_dim1], &i__2);
	wr1 = imuk[i__] + 1;
	wc1 += i__ * nui;
/* L100: */
    }

/*     Determination of the remaining nontrivial matrices in Vij,k */
/*     block column by block column with decreasing block row index. */

/*     The computation starts with the second block column since V11,0 */
/*     has already been determined. */
/*     The coefficients Vij,k satisfy the recurrence relation: */

/*        Vij,k = Sum(r=i+1,...,j-k)   AHir*Vrj,k + */
/*              + Sum(r=i+1,...,j-k+1) EHir*Vrj,k-1,   i + k < j, */

/*              = EHi,i+1 * Vi+1,j,k-1                 i + k = j. */

/*     This recurrence relation can be derived from [1], (4.6.8) */
/*     and formula (1) in Section PURPOSE. */

    vc1 = imuk[1] - inuk[1] + 1;
    ari = 1;

    i__1 = *nblcks;
    for (j = 2; j <= i__1; ++j) {
	dif = imuk[j] - imuk[j - 1] - inuk[j];
	ari += inuk[j - 1];
	ark = ari;

/*        Computation of the matrices Vij,k where i + k < j. */
/*        Each matrix Vij,k has dimension mu(i)-by-(mu(j) - nu(j)). */

	i__2 = j - 2;
	for (k = 0; k <= i__2; ++k) {

/*           VC1, VC2 are the first and last column index of Vij,k. */

	    vc2 = vc1 + dif - 1;
	    ac2 = imuk[j - k];
	    ar1 = ark;
	    ark -= inuk[j - k - 1];

	    for (i__ = j - k - 1; i__ >= 1; --i__) {

/*              Compute the first part of Vij,k in decreasing order: */
/*              Vij,k := Vij,k + Sum(r=i+1,..,j-k) AHir*Vrj,k. */
/*              The non-zero parts of AHir are stored in */
/*              A(AR1:AR1+nu(i)-1,AC1:AC2) and Vrj,k are stored in */
/*              VEPS(AC1:AC2,VC1:VC2). */
/*              The non-zero part of the result is stored in */
/*              VEPS(VR1:VR2,VC1:VC2). */

		vr2 = imuk[i__];
		ac1 = vr2 + 1;
		vr1 = ac1 - inuk[i__];
		ar1 -= inuk[i__];
		i__3 = ac2 - vr2;
		dgemm_("No transpose", "No transpose", &inuk[i__], &dif, &
			i__3, &c_b26, &a[ar1 + ac1 * a_dim1], lda, &veps[ac1 
			+ vc1 * veps_dim1], ldveps, &c_b26, &veps[vr1 + vc1 * 
			veps_dim1], ldveps, (ftnlen)12, (ftnlen)12);
/* L120: */
	    }

	    er1 = 1;

	    i__3 = j - k - 1;
	    for (i__ = 1; i__ <= i__3; ++i__) {

/*              Compute the second part of Vij,k+1 in normal order: */
/*              Vij,k+1 := Sum(r=i+1,..,j-k) EHir*Vrj,k. */
/*              The non-zero parts of EHir are stored in */
/*              E(ER1:ER1+nu(i)-1,EC1:AC2) and Vrj,k are stored in */
/*              VEPS(EC1:AC2,VC1:VC2). */
/*              The non-zero part of the result is stored in */
/*              VEPS(VR1:VR2,VC2+1:VC2+DIF), where */
/*              DIF = VC2 - VC1 + 1 = mu(j) - nu(j). */
/*              This code portion also computes Vij,k+1 for i + k = j. */

		vr2 = imuk[i__];
		ec1 = vr2 + 1;
		vr1 = ec1 - inuk[i__];
		i__4 = ac2 - vr2;
		dgemm_("No transpose", "No transpose", &inuk[i__], &dif, &
			i__4, &c_b26, &e[er1 + ec1 * e_dim1], lde, &veps[ec1 
			+ vc1 * veps_dim1], ldveps, &c_b16, &veps[vr1 + (vc2 
			+ 1) * veps_dim1], ldveps, (ftnlen)12, (ftnlen)12);
		er1 += inuk[i__];
/* L140: */
	    }

	    vc1 = vc2 + 1;
/* L160: */
	}

	vc1 += dif;
/* L180: */
    }

/*     Restore original contents of the array IMUK. */

/*     Since, at the moment: */
/*       IMUK(i) = Sum(j=1,...,i) mu(j),   (i=1,...,NBLCKS), */
/*     the original values are: */
/*       mu(i) = IMUK(i) - IMUK(i-1)  with IMUK(0 ) = 0. */

    smui1 = 0;

    i__1 = *nblcks;
    for (i__ = 1; i__ <= i__1; ++i__) {
	smui = imuk[i__];
	imuk[i__] = smui - smui1;
	smui1 = smui;
/* L200: */
    }

    return 0;
/* *** Last line of MC03NY *** */
} /* mc03ny_ */

