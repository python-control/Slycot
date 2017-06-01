/* MB03TD.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;

/* Subroutine */ int mb03td_(char *typ, char *compu, logical *select, logical 
	*lower, integer *n, doublereal *a, integer *lda, doublereal *g, 
	integer *ldg, doublereal *u1, integer *ldu1, doublereal *u2, integer *
	ldu2, doublereal *wr, doublereal *wi, integer *m, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen typ_len, ftnlen compu_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, u1_dim1, u1_offset, u2_dim1, 
	    u2_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, ks, nbf, nbl, here;
    static logical pair;
    static integer ierr, ifst;
    static logical flow, swap;
    static integer ilst;
    static logical isham;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03ts_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *);
    static logical wantu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer nbnext, wrkmin;


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

/*     To reorder a matrix X in skew-Hamiltonian Schur form: */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G = -G, */
/*                   [  0   A  ] */

/*     or in Hamiltonian Schur form: */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G =  G, */
/*                   [  0  -A  ] */

/*     where A is in upper quasi-triangular form, so that a selected */
/*     cluster of eigenvalues appears in the leading diagonal blocks */
/*     of the matrix A (in X) and the leading columns of [ U1; -U2 ] form */
/*     an orthonormal basis for the corresponding right invariant */
/*     subspace. */

/*     If X is skew-Hamiltonian, then each eigenvalue appears twice; one */
/*     copy corresponds to the j-th diagonal element and the other to the */
/*     (n+j)-th diagonal element of X. The logical array LOWER controls */
/*     which copy is to be reordered to the leading part of A. */

/*     If X is Hamiltonian then the eigenvalues appear in pairs */
/*     (lambda,-lambda); lambda corresponds to the j-th diagonal */
/*     element and -lambda to the (n+j)-th diagonal element of X. */
/*     The logical array LOWER controls whether lambda or -lambda is to */
/*     be reordered to the leading part of A. */

/*     The matrix A must be in Schur canonical form (as returned by the */
/*     LAPACK routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYP     CHARACTER*1 */
/*             Specifies the type of the input matrix X: */
/*             = 'S': X is skew-Hamiltonian; */
/*             = 'H': X is Hamiltonian. */

/*     COMPU   CHARACTER*1 */
/*             = 'U': update the matrices U1 and U2 containing the */
/*                    Schur vectors; */
/*             = 'N': do not update U1 and U2. */

/*     SELECT  (input/output) LOGICAL array, dimension (N) */
/*             SELECT specifies the eigenvalues in the selected cluster. */
/*             To select a real eigenvalue w(j), SELECT(j) must be set */
/*             to .TRUE.. To select a complex conjugate pair of */
/*             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2 */
/*             diagonal block, both SELECT(j) and SELECT(j+1) must be set */
/*             to .TRUE.; a complex conjugate pair of eigenvalues must be */
/*             either both included in the cluster or both excluded. */

/*     LOWER   (input/output) LOGICAL array, dimension (N) */
/*             LOWER controls which copy of a selected eigenvalue is */
/*             included in the cluster. If SELECT(j) is set to .TRUE. */
/*             for a real eigenvalue w(j); then LOWER(j) must be set to */
/*             .TRUE. if the eigenvalue corresponding to the (n+j)-th */
/*             diagonal element of X is to be reordered to the leading */
/*             part; and LOWER(j) must be set to .FALSE. if the */
/*             eigenvalue corresponding to the j-th diagonal element of */
/*             X is to be reordered to the leading part. Similarly, for */
/*             a complex conjugate pair of eigenvalues w(j) and w(j+1), */
/*             both LOWER(j) and LOWER(j+1) must be set to .TRUE. if the */
/*             eigenvalues corresponding to the (n+j:n+j+1,n+j:n+j+1) */
/*             diagonal block of X are to be reordered to the leading */
/*             part. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A in Schur */
/*             canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the reordered matrix A, again in Schur canonical form, */
/*             with the selected eigenvalues in the diagonal blocks. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, if TYP = 'S', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the skew-symmetric matrix G. The rest of this array is not */
/*             referenced. */
/*             On entry, if TYP = 'H', the leading N-by-N part of this */
/*             array must contain the upper triangular part of the */
/*             symmetric matrix G. The rest of this array is not */
/*             referenced. */
/*             On exit, if TYP = 'S', the leading N-by-N part of this */
/*             array contains the strictly upper triangular part of the */
/*             skew-symmetric matrix G, updated by the orthogonal */
/*             symplectic transformation which reorders X. */
/*             On exit, if TYP = 'H', the leading N-by-N part of this */
/*             array contains the upper triangular part of the symmetric */
/*             matrix G, updated by the orthogonal symplectic */
/*             transformation which reorders X. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if COMPU = 'U', the leading N-by-N part of this */
/*             array must contain U1, the (1,1) block of an orthogonal */
/*             symplectic matrix U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U', the leading N-by-N part of this */
/*             array contains the (1,1) block of the matrix U, */
/*             postmultiplied by the orthogonal symplectic transformation */
/*             which reorders X. The leading M columns of U form an */
/*             orthonormal basis for the specified invariant subspace. */
/*             If COMPU = 'N', this array is not referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1. */
/*             LDU1 >= MAX(1,N),  if COMPU = 'U'; */
/*             LDU1 >= 1,         otherwise. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if COMPU = 'U', the leading N-by-N part of this */
/*             array must contain U2, the (1,2) block of an orthogonal */
/*             symplectic matrix U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U', the leading N-by-N part of this */
/*             array contains the (1,2) block of the matrix U, */
/*             postmultiplied by the orthogonal symplectic transformation */
/*             which reorders X. */
/*             If COMPU = 'N', this array is not referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2. */
/*             LDU2 >= MAX(1,N),  if COMPU = 'U'; */
/*             LDU2 >= 1,         otherwise. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             reordered eigenvalues of A. The eigenvalues are stored */
/*             in the same order as on the diagonal of A, with */
/*             WR(i) = A(i,i) and, if A(i:i+1,i:i+1) is a 2-by-2 diagonal */
/*             block, WI(i) > 0 and WI(i+1) = -WI(i). Note that if an */
/*             eigenvalue is sufficiently ill-conditioned, then its value */
/*             may differ significantly from its value before reordering. */

/*     M       (output) INTEGER */
/*             The dimension of the specified invariant subspace. */
/*             0 <= M <= N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -18,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*             = 1:  reordering of X failed because some eigenvalue pairs */
/*                   are too close to separate (the problem is very */
/*                   ill-conditioned); X may have been partially */
/*                   reordered, and WR and WI contain the eigenvalues in */
/*                   the same order as in X. */

/*     REFERENCES */

/*     [1] Bai, Z. and Demmel, J.W. */
/*         On Swapping Diagonal Blocks in Real Schur Form. */
/*         Linear Algebra Appl., 186, pp. 73-95, 1993. */

/*     [2] Benner, P., Kressner, D., and Mehrmann, V. */
/*         Skew-Hamiltonian and Hamiltonian Eigenvalue Problems: Theory, */
/*         Algorithms and Applications. Techn. Report, TU Berlin, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAORD). */

/*     KEYWORDS */

/*     Hamiltonian matrix, skew-Hamiltonian matrix, invariant subspace. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and check input parameters. */

    /* Parameter adjustments */
    --select;
    --lower;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    --wr;
    --wi;
    --dwork;

    /* Function Body */
    isham = lsame_(typ, "H", (ftnlen)1, (ftnlen)1);
    wantu = lsame_(compu, "U", (ftnlen)1, (ftnlen)1);
    wrkmin = max(1,*n);
    *info = 0;
    if (! isham && ! lsame_(typ, "S", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (! wantu && ! lsame_(compu, "N", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldg < max(1,*n)) {
	*info = -9;
    } else if (*ldu1 < 1 || wantu && *ldu1 < *n) {
	*info = -11;
    } else if (*ldu2 < 1 || wantu && *ldu2 < *n) {
	*info = -13;
    } else if (*ldwork < wrkmin) {
	*info = -18;
	dwork[1] = (doublereal) wrkmin;
    }

/*     Return if there were illegal values. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("MB03TD", &i__1, (ftnlen)6);
	return 0;
    }

/*     Set M to the dimension of the specified invariant subspace. */

    *m = 0;
    pair = FALSE_;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (pair) {
	    pair = FALSE_;
	} else {
	    if (k < *n) {
		if (a[k + 1 + k * a_dim1] == 0.) {
		    if (select[k]) {
			++(*m);
		    }
		} else {
		    pair = TRUE_;
		    if (select[k] || select[k + 1]) {
			*m += 2;
		    }
		}
	    } else {
		if (select[*n]) {
		    ++(*m);
		}
	    }
	}
/* L10: */
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Collect the selected blocks at the top-left corner of X. */

    ks = 0;
    pair = FALSE_;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (pair) {
	    pair = FALSE_;
	} else {
	    swap = select[k];
	    flow = lower[k];
	    if (k < *n) {
		if (a[k + 1 + k * a_dim1] != 0.) {
		    pair = TRUE_;
		    swap = swap || select[k + 1];
		    flow = flow || lower[k + 1];
		}
	    }

	    if (pair) {
		nbf = 2;
	    } else {
		nbf = 1;
	    }

	    if (swap) {
		++ks;
		if (flow) {

/*                 Step 1: Swap the K-th block to position N. */

		    ifst = k;
		    ilst = *n;
		    nbl = 1;
		    if (ilst > 1) {
			if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
			    --ilst;
			    nbl = 2;
			}
		    }

/*                 Update ILST. */

		    if (nbf == 2 && nbl == 1) {
			--ilst;
		    }
		    if (nbf == 1 && nbl == 2) {
			++ilst;
		    }

		    if (ilst == ifst) {
			goto L30;
		    }

		    here = ifst;

L20:

/*                 Swap block with next one below. */

		    if (nbf == 1 || nbf == 2) {

/*                    Current block is either 1-by-1 or 2-by-2. */

			nbnext = 1;
			if (here + nbf + 1 <= *n) {
			    if (a[here + nbf + 1 + (here + nbf) * a_dim1] != 
				    0.) {
				nbnext = 2;
			    }
			}
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &here, &nbf, &nbnext, &
				dwork[1], &ierr);
			if (ierr != 0) {
			    *info = 1;
			    goto L70;
			}
			here += nbnext;

/*                    Test if 2-by-2 block breaks into two 1-by-1 blocks. */

			if (nbf == 2) {
			    if (a[here + 1 + here * a_dim1] == 0.) {
				nbf = 3;
			    }
			}

		    } else {

/*                    Current block consists of two 1-by-1 blocks each of */
/*                    which must be swapped individually. */

			nbnext = 1;
			if (here + 3 <= *n) {
			    if (a[here + 3 + (here + 2) * a_dim1] != 0.) {
				nbnext = 2;
			    }
			}
			i__2 = here + 1;
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__1, &nbnext, &
				dwork[1], &ierr);
			if (ierr != 0) {
			    *info = 1;
			    goto L70;
			}
			if (nbnext == 1) {

/*                       Swap two 1-by-1 blocks, no problems possible. */

			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &here, &c__1, &nbnext, &
				    dwork[1], &ierr);
			    ++here;
			} else {

/*                       Recompute NBNEXT in case 2 by 2 split. */

			    if (a[here + 2 + (here + 1) * a_dim1] == 0.) {
				nbnext = 1;
			    }
			    if (nbnext == 2) {

/*                          2-by-2 block did not split */

				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &here, &
					c__1, &nbnext, &dwork[1], &ierr);
				if (ierr != 0) {
				    *info = 1;
				    goto L70;
				}
				here += 2;
			    } else {

/*                          2-by-2 block did split */

				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &here, &
					c__1, &c__1, &dwork[1], &ierr);
				i__2 = here + 1;
				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &i__2, &
					c__1, &c__1, &dwork[1], &ierr);
				here += 2;
			    }
			}
		    }
		    if (here < ilst) {
			goto L20;
		    }

L30:

/*                 Step 2: Apply an orthogonal symplectic transformation */
/*                         to swap the last blocks in A and -A' (or A'). */

		    if (nbf == 1) {

/*                    Exchange columns/rows N <-> 2*N. No problems */
/*                    possible. */

			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);

		    } else if (nbf == 2) {

/*                    Swap last block with its equivalent by an */
/*                    orthogonal symplectic transformation. */

			i__2 = *n - 1;
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__2, &c__2, &dwork[
				1], &ierr);
			if (ierr != 0) {
			    *info = 1;
			    goto L70;
			}

/*                    Test if 2-by-2 block breaks into two 1-by-1 blocks. */

			if (a[*n - 1 + *n * a_dim1] == 0.) {
			    nbf = 3;
			}
		    } else {

/*                    Block did split. Swap (N-1)-th and N-th elements */
/*                    consecutively by symplectic generalized */
/*                    permutations and one rotation. */

			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);
			i__2 = *n - 1;
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__1, &c__1, &dwork[
				1], &ierr);
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);
		    }
		    ifst = *n;
		    if (pair) {
			ifst = *n - 1;
		    }
		} else {
		    ifst = k;
		}

/*              Step 3: Swap the K-th / N-th block to position KS. */

		ilst = ks;
		nbl = 1;
		if (ilst > 1) {
		    if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
			--ilst;
			nbl = 2;
		    }
		}

		if (ilst == ifst) {
		    goto L50;
		}

		here = ifst;
L40:

/*              Swap block with next one above. */

		if (nbf == 1 || nbf == 2) {

/*                 Current block either 1 by 1 or 2 by 2. */

		    nbnext = 1;
		    if (here >= 3) {
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
			    nbnext = 2;
			}
		    }
		    i__2 = here - nbnext;
		    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[g_offset]
			    , ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2,
			     &i__2, &nbnext, &nbf, &dwork[1], &ierr);
		    if (ierr != 0) {
			*info = 1;
			goto L70;
		    }
		    here -= nbnext;

/*                 Test if 2-by-2 block breaks into two 1-by-1 blocks. */

		    if (nbf == 2) {
			if (a[here + 1 + here * a_dim1] == 0.) {
			    nbf = 3;
			}
		    }

		} else {

/*                 Current block consists of two 1 by 1 blocks each of */
/*                 which must be swapped individually. */

		    nbnext = 1;
		    if (here >= 3) {
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
			    nbnext = 2;
			}
		    }
		    i__2 = here - nbnext;
		    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[g_offset]
			    , ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2,
			     &i__2, &nbnext, &c__1, &dwork[1], &ierr);
		    if (ierr != 0) {
			*info = 1;
			goto L70;
		    }
		    if (nbnext == 1) {

/*                    Swap two 1-by-1 blocks, no problems possible. */

			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &here, &nbnext, &c__1, &
				dwork[1], &ierr);
			--here;
		    } else {

/*                    Recompute NBNEXT in case 2-by-2 split. */

			if (a[here + (here - 1) * a_dim1] == 0.) {
			    nbnext = 1;
			}
			if (nbnext == 2) {

/*                       2-by-2 block did not split */

			    i__2 = here - 1;
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &i__2, &c__2, &c__1, &
				    dwork[1], &ierr);
			    if (ierr != 0) {
				*info = 1;
				goto L70;
			    }
			    here += -2;
			} else {

/*                       2-by-2 block did split */

			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &here, &c__1, &c__1, &
				    dwork[1], &ierr);
			    i__2 = here - 1;
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &i__2, &c__1, &c__1, &
				    dwork[1], &ierr);
			    here += -2;
			}
		    }
		}

		if (here > ilst) {
		    goto L40;
		}

L50:
		if (pair) {
		    ++ks;
		}
	    }
	}
/* L60: */
    }

L70:

/*     Store eigenvalues. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	wr[k] = a[k + k * a_dim1];
	wi[k] = 0.;
/* L80: */
    }
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	if (a[k + 1 + k * a_dim1] != 0.) {
	    wi[k] = sqrt((d__1 = a[k + (k + 1) * a_dim1], abs(d__1))) * sqrt((
		    d__2 = a[k + 1 + k * a_dim1], abs(d__2)));
	    wi[k + 1] = -wi[k];
	}
/* L90: */
    }

    dwork[1] = (doublereal) wrkmin;

    return 0;
/* *** Last line of MB03TD *** */
} /* mb03td_ */

