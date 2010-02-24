      SUBROUTINE MB05MY( BALANC, N, A, LDA, WR, WI, R, LDR, Q, LDQ,
     $                   DWORK, LDWORK, INFO )
C
C     SLICOT RELEASE 5.0.
C
C     Copyright (c) 2002-2009 NICONET e.V.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 2 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To compute, for an N-by-N real nonsymmetric matrix A, the
C     orthogonal matrix Q reducing it to real Schur form T, the
C     eigenvalues, and the right eigenvectors of T.
C
C     The right eigenvector r(j) of T satisfies
C                      T * r(j) = lambda(j) * r(j)
C     where lambda(j) is its eigenvalue.
C
C     The matrix of right eigenvectors R is upper triangular, by
C     construction.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALANC  CHARACTER*1
C             Indicates how the input matrix should be diagonally scaled
C             to improve the conditioning of its eigenvalues as follows:
C             = 'N':  Do not diagonally scale;
C             = 'S':  Diagonally scale the matrix, i.e. replace A by
C                     D*A*D**(-1), where D is a diagonal matrix chosen
C                     to make the rows and columns of A more equal in
C                     norm. Do not permute.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the given matrix A.
C             On exit, the leading N-by-N upper quasi-triangular part of
C             this array contains the real Schur canonical form of A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,N).
C
C     WR      (output) DOUBLE PRECISION array, dimension (N)
C     WI      (output) DOUBLE PRECISION array, dimension (N)
C             WR and WI contain the real and imaginary parts,
C             respectively, of the computed eigenvalues. Complex
C             conjugate pairs of eigenvalues appear consecutively
C             with the eigenvalue having the positive imaginary part
C             first.
C
C     R       (output) DOUBLE PRECISION array, dimension (LDR,N)
C             The leading N-by-N upper triangular part of this array
C             contains the matrix of right eigenvectors R, in the same
C             order as their eigenvalues. The real and imaginary parts
C             of a complex eigenvector corresponding to an eigenvalue
C             with positive imaginary part are stored in consecutive
C             columns. (The corresponding conjugate eigenvector is not
C             stored.) The eigenvectors are not backward transformed
C             for balancing (when BALANC = 'S').
C
C     LDR     INTEGER
C             The leading dimension of array R.  LDR >= max(1,N).
C
C     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
C             The leading N-by-N part of this array contains the
C             orthogonal matrix Q which has reduced A to real Schur
C             form.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.  LDQ >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK.
C             If BALANC = 'S', DWORK(2),...,DWORK(N+1) return the
C             scaling factors used for balancing.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= max(1,4*N).
C             For good performance, LDWORK must generally be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, the QR algorithm failed to compute all
C                   the eigenvalues, and no eigenvectors have been
C                   computed; elements i+1:N of WR and WI contain
C                   eigenvalues which have converged.
C
C     METHOD
C
C     This routine uses the QR algorithm to obtain the real Schur form
C     T of matrix A. Then, the right eigenvectors of T are computed,
C     but they are not backtransformed into the eigenvectors of A.
C     MB05MY is a modification of the LAPACK driver routine DGEEV.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997.
C     Supersedes Release 2.0 routine MB05AY.
C
C     REVISIONS
C
C     V. Sima, April 25, 2003, Feb. 15, 2004.
C
C     KEYWORDS
C
C     Eigenvalue, eigenvector decomposition, real Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         BALANC
      INTEGER           INFO, LDA, LDQ, LDR, LDWORK, N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), DWORK( * ), Q( LDQ, * ),
     $                  R( LDR, * ), WI( * ), WR( * )
C     ..
C     .. Local Scalars ..
      LOGICAL           SCALE, SCALEA
      INTEGER           HSDWOR, IBAL, IERR, IHI, ILO, ITAU, JWORK, K,
     $                  MAXB, MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION  ANRM, BIGNUM, CSCALE, EPS, SMLNUM
C     ..
C     .. Local Arrays ..
      LOGICAL           SELECT( 1 )
      DOUBLE PRECISION  DUM( 1 )
C     ..
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, ILAENV, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL          DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY, DLASCL,
     $                  DORGHR, DTREVC, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar arguments.
C
      INFO  = 0
      SCALE = LSAME( BALANC, 'S' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. SCALE ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDR.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
C
C     Compute workspace.
C      (Note: Comments in the code beginning "Workspace:" describe the
C       minimal amount of workspace needed at that point in the code,
C       as well as the preferred amount for good performance.
C       NB refers to the optimal block size for the immediately
C       following subroutine, as returned by ILAENV.
C       HSDWOR refers to the workspace preferred by DHSEQR, as
C       calculated below. HSDWOR is computed assuming ILO=1 and IHI=N,
C       the worst case.)
C
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LDWORK.GE.1 ) THEN
         MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
         MINWRK = MAX( 1, 4*N )
         MAXWRK = MAX( MAXWRK, 2*N+( N-1 )*
     $            ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )
         MAXB = MAX( ILAENV( 8, 'DHSEQR', 'SV', N, 1, N, -1 ), 2 )
         K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'SV', N, 1,
     $       N, -1 ) ) )
         HSDWOR = MAX( K*( K+2 ), 2*N )
         MAXWRK = MAX( MAXWRK, N+1, N+HSDWOR )
         MAXWRK = MAX( MAXWRK, 4*N )
         DWORK( 1 ) = MAXWRK
      END IF
      IF( LDWORK.LT.MINWRK ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB05MY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
C     Get machine constants.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
C
C     Scale A if max element outside range [SMLNUM,BIGNUM].
C
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
C
C     Balance the matrix, if requested. (Permutation is not possible.)
C     (Workspace: need N)
C
      IBAL = 1
      CALL DGEBAL( BALANC, N, A, LDA, ILO, IHI, DWORK( IBAL ), IERR )
C
C     Reduce to upper Hessenberg form.
C     (Workspace: need 3*N, prefer 2*N+N*NB)
C
      ITAU = IBAL + N
      JWORK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, DWORK( ITAU ), DWORK( JWORK ),
     $             LDWORK-JWORK+1, IERR )
C
C     Compute right eigenvectors of T.
C     Copy Householder vectors to Q.
C
      CALL DLACPY( 'Lower', N, N, A, LDA, Q, LDQ )
C
C     Generate orthogonal matrix in Q.
C     (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
C
      CALL DORGHR( N, ILO, IHI, Q, LDQ, DWORK( ITAU ), DWORK( JWORK ),
     $             LDWORK-JWORK+1, IERR )
C
C     Perform QR iteration, accumulating Schur vectors in Q.
C     (Workspace: need N+1, prefer N+HSDWOR (see comments) )
C
      JWORK = ITAU
      CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, Q, LDQ,
     $             DWORK( JWORK ), LDWORK-JWORK+1, INFO )
C
C     If INFO > 0 from DHSEQR, then quit.
C
      IF( INFO.GT.0 )
     $   GO TO 10
C
C     Compute right eigenvectors of T in R.
C     (Workspace: need 4*N)
C
      CALL DTREVC( 'Right', 'All', SELECT, N, A, LDA, DUM, 1, R, LDR, N,
     $             NOUT, DWORK( JWORK ), IERR )
C
C     Undo scaling if necessary.
C
   10 CONTINUE
      IF( SCALEA ) THEN
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N,
     $                   IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N,
     $                   IERR )
         END IF
      END IF
C
      IF ( SCALE ) THEN
         DO 20 K = N, 1, -1
            DWORK( K+1 ) = DWORK( K )
   20    CONTINUE
      END IF
      DWORK( 1 ) = MAXWRK
C
      RETURN
C *** Last line of MB05MY ***
      END
