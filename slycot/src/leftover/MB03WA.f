      SUBROUTINE MB03WA( WANTQ, WANTZ, N1, N2, A, LDA, B, LDB, Q, LDQ,
     $                   Z, LDZ, INFO )
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
C     To swap adjacent diagonal blocks A11*B11 and A22*B22 of size
C     1-by-1 or 2-by-2 in an upper (quasi) triangular matrix product
C     A*B by an orthogonal equivalence transformation.
C
C     (A, B) must be in periodic real Schur canonical form (as returned
C     by SLICOT Library routine MB03XP), i.e., A is block upper
C     triangular with 1-by-1 and 2-by-2 diagonal blocks, and B is upper
C     triangular.
C
C     Optionally, the matrices Q and Z of generalized Schur vectors are
C     updated.
C
C         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)',
C         Z(in) * B(in) * Q(in)' = Z(out) * B(out) * Q(out)'.
C
C     This routine is largely based on the LAPACK routine DTGEX2
C     developed by Bo Kagstrom and Peter Poromaa.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     WANTQ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Q as follows:
C             = .TRUE. :  The matrix Q is updated;
C             = .FALSE.:  the matrix Q is not required.
C
C     WANTZ   LOGICAL
C             Indicates whether or not the user wishes to accumulate
C             the matrix Z as follows:
C             = .TRUE. :  The matrix Z is updated;
C             = .FALSE.:  the matrix Z is not required.
C
C     Input/Output Parameters
C
C     N1      (input) INTEGER
C             The order of the first block A11*B11. N1 = 0, 1 or 2.
C
C     N2      (input) INTEGER
C             The order of the second block A22*B22. N2 = 0, 1 or 2.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C             (LDA,N1+N2)
C             On entry, the leading (N1+N2)-by-(N1+N2) part of this
C             array must contain the matrix A.
C             On exit, the leading (N1+N2)-by-(N1+N2) part of this array
C             contains the matrix A of the reordered pair.
C
C     LDA     INTEGER
C             The leading dimension of the array A. LDA >= MAX(1,N1+N2).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,N1+N2)
C             On entry, the leading (N1+N2)-by-(N1+N2) part of this
C             array must contain the matrix B.
C             On exit, the leading (N1+N2)-by-(N1+N2) part of this array
C             contains the matrix B of the reordered pair.
C
C     LDB     INTEGER
C             The leading dimension of the array B. LDB >= MAX(1,N1+N2).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension
C             (LDQ,N1+N2)
C             On entry, if WANTQ = .TRUE., the leading
C             (N1+N2)-by-(N1+N2) part of this array must contain the
C             orthogonal matrix Q.
C             On exit, the leading (N1+N2)-by-(N1+N2) part of this array
C             contains the updated matrix Q. Q will be a rotation
C             matrix for N1=N2=1.
C             This array is not referenced if WANTQ = .FALSE..
C
C     LDQ     INTEGER
C             The leading dimension of the array Q. LDQ >= 1.
C             If WANTQ = .TRUE., LDQ >= N1+N2.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension
C             (LDZ,N1+N2)
C             On entry, if WANTZ = .TRUE., the leading
C             (N1+N2)-by-(N1+N2) part of this array must contain the
C             orthogonal matrix Z.
C             On exit, the leading (N1+N2)-by-(N1+N2) part of this array
C             contains the updated matrix Z. Z will be a rotation
C             matrix for N1=N2=1.
C             This array is not referenced if WANTZ = .FALSE..
C
C     LDZ     INTEGER
C             The leading dimension of the array Z. LDZ >= 1.
C             If WANTZ = .TRUE., LDZ >= N1+N2.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             = 1:  the transformed matrix (A, B) would be
C                   too far from periodic Schur form; the blocks are
C                   not swapped and (A,B) and (Q,Z) are unchanged.
C
C     METHOD
C
C     In the current code both weak and strong stability tests are
C     performed. The user can omit the strong stability test by changing
C     the internal logical parameter WANDS to .FALSE.. See ref. [2] for
C     details.
C
C     REFERENCES
C
C     [1] Kagstrom, B.
C         A direct method for reordering eigenvalues in the generalized
C         real Schur form of a regular matrix pair (A,B), in M.S. Moonen
C         et al (eds.), Linear Algebra for Large Scale and Real-Time
C         Applications, Kluwer Academic Publ., 1993, pp. 195-218.
C
C     [2] Kagstrom, B., and Poromaa, P.
C         Computing eigenspaces with specified eigenvalues of a regular
C         matrix pair (A, B) and condition estimation: Theory,
C         algorithms and software, Numer. Algorithms, 1996, vol. 12,
C         pp. 369-407.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DTGPX2).
C
C     KEYWORDS
C
C     Eigenvalue, periodic Schur form, reordering
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+01 )
      INTEGER            LDST
      PARAMETER          ( LDST = 4 )
      LOGICAL            WANDS
      PARAMETER          ( WANDS = .TRUE. )
C     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, LDA, LDB, LDQ, LDZ, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(LDB,*), Q(LDQ,*), Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL            DTRONG, WEAK
      INTEGER            I, LINFO, M
      DOUBLE PRECISION   BQRA21, BRQA21, DDUM, DNORM, DSCALE, DSUM, EPS,
     $                   F, G, SA, SB, SCALE, SMLNUM, SS, THRESH, WS
C     .. Local Arrays ..
      INTEGER            IWORK( LDST )
      DOUBLE PRECISION   AI(2), AR(2), BE(2), DWORK(32), IR(LDST,LDST),
     $                   IRCOP(LDST,LDST), LI(LDST,LDST),
     $                   LICOP(LDST,LDST), S(LDST,LDST),
     $                   SCPY(LDST,LDST), T(LDST,LDST), TAUL(LDST),
     $                   TAUR(LDST), TCPY(LDST,LDST)
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEQR2, DGERQ2, DLACPY, DLARTG, DLASET,
     $                   DLASSQ, DORG2R, DORGR2, DORM2R, DORMR2, DROT,
     $                   DSCAL, MB03YT, SB04OW
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
C
C     .. Executable Statements ..
C
      INFO = 0
C
C     Quick return if possible.
C     For efficiency, the arguments are not checked.
C
      IF ( N1.LE.0 .OR. N2.LE.0 )
     $   RETURN
      M = N1 + N2
C
      WEAK = .FALSE.
      DTRONG = .FALSE.
C
C     Make a local copy of selected block.
C
      CALL DLASET( 'All', LDST, LDST, ZERO, ZERO, LI, LDST )
      CALL DLASET( 'All', LDST, LDST, ZERO, ZERO, IR, LDST )
      CALL DLACPY( 'Full', M, M, A, LDA, S, LDST )
      CALL DLACPY( 'Full', M, M, B, LDB, T, LDST )
C
C     Compute threshold for testing acceptance of swapping.
C
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      DSCALE = ZERO
      DSUM = ONE
      CALL DLACPY( 'Full', M, M, S, LDST, DWORK, M )
      CALL DLASSQ( M*M, DWORK, 1, DSCALE, DSUM )
      CALL DLACPY( 'Full', M, M, T, LDST, DWORK, M )
      CALL DLASSQ( M*M, DWORK, 1, DSCALE, DSUM )
      DNORM = DSCALE*SQRT( DSUM )
      THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
C
      IF ( M.EQ.2 ) THEN
C
C        CASE 1: Swap 1-by-1 and 1-by-1 blocks.
C
C        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks
C        using Givens rotations and perform the swap tentatively.
C
         F =  S(2,2)*T(2,2) - T(1,1)*S(1,1)
         G = -S(2,2)*T(1,2) - T(1,1)*S(1,2)
         SB = ABS( T(1,1) )
         SA = ABS( S(2,2) )
         CALL DLARTG( F, G, IR(1,2), IR(1,1), DDUM )
         IR(2,1) = -IR(1,2)
         IR(2,2) =  IR(1,1)
         CALL DROT( 2, S(1,1), 1, S(1,2), 1, IR(1,1), IR(2,1) )
         CALL DROT( 2, T(1,1), LDST, T(2,1), LDST, IR(1,1), IR(2,1) )
         IF( SA.GE.SB ) THEN
            CALL DLARTG( S(1,1), S(2,1), LI(1,1), LI(2,1), DDUM )
         ELSE
            CALL DLARTG( T(2,2), T(2,1), LI(1,1), LI(2,1), DDUM )
            LI(2,1) = -LI(2,1)
         END IF
         CALL DROT( 2, S(1,1), LDST, S(2,1), LDST, LI(1,1), LI(2,1) )
         CALL DROT( 2, T(1,1), 1, T(1,2), 1, LI(1,1), LI(2,1) )
         LI(2,2) =  LI(1,1)
         LI(1,2) = -LI(2,1)
C
C        Weak stability test:
C           |S21| + |T21| <= O(EPS * F-norm((S, T))).
C
         WS = ABS( S(2,1) ) + ABS( T(2,1) )
         WEAK = WS.LE.THRESH
         IF ( .NOT.WEAK )
     $      GO TO 50
C
         IF ( WANDS ) THEN
C
C           Strong stability test:
C             F-norm((A-QL'*S*QR, B-QR'*T*QL)) <= O(EPS*F-norm((A,B))).
C
            CALL DLACPY( 'Full', M, M, A, LDA, DWORK(M*M+1), M )
            CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE,
     $                  LI, LDST, S, LDST, ZERO, DWORK, M )
            CALL DGEMM( 'No Transpose', 'Transpose', M, M, M, -ONE,
     $                  DWORK, M, IR, LDST, ONE, DWORK(M*M+1), M )
            DSCALE = ZERO
            DSUM = ONE
            CALL DLASSQ( M*M, DWORK(M*M+1), 1, DSCALE, DSUM )
C
            CALL DLACPY( 'Full', M, M, B, LDB, DWORK(M*M+1), M )
            CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE,
     $                  IR, LDST, T, LDST, ZERO, DWORK, M )
            CALL DGEMM( 'No Transpose', 'Transpose', M, M, M, -ONE,
     $                  DWORK, M, LI, LDST, ONE, DWORK(M*M+1), M )
            CALL DLASSQ( M*M, DWORK(M*M+1), 1, DSCALE, DSUM )
            SS = DSCALE*SQRT( DSUM )
            DTRONG = SS.LE.THRESH
            IF( .NOT.DTRONG )
     $         GO TO 50
         END IF
C
C        Update A and B.
C
         CALL DLACPY( 'All', M, M, S, LDST, A, LDA )
         CALL DLACPY( 'All', M, M, T, LDST, B, LDB )
C
C        Set  N1-by-N2 (2,1) - blocks to ZERO.
C
         A(2,1) = ZERO
         B(2,1) = ZERO
C
C        Accumulate transformations into Q and Z if requested.
C
         IF ( WANTQ )
     $      CALL DROT( 2, Q(1,1), 1, Q(1,2), 1, LI(1,1), LI(2,1) )
         IF ( WANTZ )
     $      CALL DROT( 2, Z(1,1), 1, Z(1,2), 1, IR(1,1), IR(2,1) )
C
C        Exit with INFO = 0 if swap was successfully performed.
C
         RETURN
C
      ELSE
C
C        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2
C                and 2-by-2 blocks.
C
C        Solve the periodic Sylvester equation
C                 S11 * R - L * S22 = SCALE * S12
C                 T11 * L - R * T22 = SCALE * T12
C        for R and L. Solutions in IR and LI.
C
         CALL DLACPY( 'Full', N1, N2, T(1,N1+1), LDST, LI, LDST )
         CALL DLACPY( 'Full', N1, N2, S(1,N1+1), LDST, IR(N2+1,N1+1),
     $                LDST )
         CALL SB04OW( N1, N2, S, LDST, S(N1+1,N1+1), LDST,
     $                IR(N2+1,N1+1), LDST, T, LDST, T(N1+1,N1+1), LDST,
     $                LI, LDST, SCALE, IWORK, LINFO )
         IF ( LINFO.NE.0 )
     $      GO TO 50
C
C        Compute orthogonal matrix QL:
C
C                    QL' * LI = [ TL ]
C                               [ 0  ]
C        where
C                    LI =  [      -L              ].
C                          [ SCALE * identity(N2) ]
C
         DO 10 I = 1, N2
            CALL DSCAL( N1, -ONE, LI(1,I), 1 )
            LI(N1+I,I) = SCALE
   10    CONTINUE
         CALL DGEQR2( M, N2, LI, LDST, TAUL, DWORK, LINFO )
         CALL DORG2R( M, M, N2, LI, LDST, TAUL, DWORK, LINFO )
C
C        Compute orthogonal matrix RQ:
C
C                    IR * RQ' =   [ 0  TR],
C
C         where IR = [ SCALE * identity(N1), R ].
C
         DO 20 I = 1, N1
            IR(N2+I,I) = SCALE
   20    CONTINUE
         CALL DGERQ2( N1, M, IR(N2+1,1), LDST, TAUR, DWORK, LINFO )
         CALL DORGR2( M, M, N1, IR, LDST, TAUR, DWORK, LINFO )
C
C        Perform the swapping tentatively:
C
         CALL DGEMM( 'Transpose', 'No Transpose', M, M, M, ONE, LI,
     $               LDST, S, LDST, ZERO, DWORK, M )
         CALL DGEMM( 'No Transpose', 'Transpose', M, M, M, ONE, DWORK,
     $               M, IR, LDST, ZERO, S, LDST )
         CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE, IR,
     $               LDST, T, LDST, ZERO, DWORK, M )
         CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE,
     $               DWORK, M, LI, LDST, ZERO, T, LDST )
         CALL DLACPY( 'All', M, M, S, LDST, SCPY, LDST )
         CALL DLACPY( 'All', M, M, T, LDST, TCPY, LDST )
         CALL DLACPY( 'All', M, M, IR, LDST, IRCOP, LDST )
         CALL DLACPY( 'All', M, M, LI, LDST, LICOP, LDST )
C
C        Triangularize the B-part by a QR factorization.
C        Apply transformation (from left) to A-part, giving S.
C
         CALL DGEQR2( M, M, T, LDST, TAUR, DWORK, LINFO )
         CALL DORM2R( 'Right', 'No Transpose', M, M, M, T, LDST, TAUR,
     $                S, LDST, DWORK, LINFO )
         CALL DORM2R( 'Left', 'Transpose', M, M, M, T, LDST, TAUR,
     $                IR, LDST, DWORK, LINFO )
C
C        Compute F-norm(S21) in BRQA21. (T21 is 0.)
C
         DSCALE = ZERO
         DSUM = ONE
         DO 30 I = 1, N2
            CALL DLASSQ( N1, S(N2+1,I), 1, DSCALE, DSUM )
   30    CONTINUE
         BRQA21 = DSCALE*SQRT( DSUM )
C
C        Triangularize the B-part by an RQ factorization.
C        Apply transformation (from right) to A-part, giving S.
C
         CALL DGERQ2( M, M, TCPY, LDST, TAUL, DWORK, LINFO )
         CALL DORMR2( 'Left', 'No Transpose', M, M, M, TCPY, LDST,
     $                TAUL, SCPY, LDST, DWORK, LINFO )
         CALL DORMR2( 'Right', 'Transpose', M, M, M, TCPY, LDST,
     $                TAUL, LICOP, LDST, DWORK, LINFO )
C
C        Compute F-norm(S21) in BQRA21. (T21 is 0.)
C
         DSCALE = ZERO
         DSUM = ONE
         DO 40 I = 1, N2
            CALL DLASSQ( N1, SCPY(N2+1,I), 1, DSCALE, DSUM )
   40    CONTINUE
         BQRA21 = DSCALE*SQRT( DSUM )
C
C        Decide which method to use.
C          Weak stability test:
C             F-norm(S21) <= O(EPS * F-norm((S, T)))
C
         IF ( BQRA21.LE.BRQA21 .AND. BQRA21.LE.THRESH ) THEN
            CALL DLACPY( 'All', M, M, SCPY, LDST, S, LDST )
            CALL DLACPY( 'All', M, M, TCPY, LDST, T, LDST )
            CALL DLACPY( 'All', M, M, IRCOP, LDST, IR, LDST )
            CALL DLACPY( 'All', M, M, LICOP, LDST, LI, LDST )
         ELSE IF ( BRQA21.GE.THRESH ) THEN
            GO TO 50
         END IF
C
C        Set lower triangle of B-part to zero
C
         CALL DLASET( 'Lower', M-1, M-1, ZERO, ZERO, T(2,1), LDST )
C
         IF ( WANDS ) THEN
C
C           Strong stability test:
C              F-norm((A-QL*S*QR', B-QR*T*QL')) <= O(EPS*F-norm((A,B)))
C
            CALL DLACPY( 'All', M, M, A, LDA, DWORK(M*M+1), M )
            CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE,
     $                  LI, LDST, S, LDST, ZERO, DWORK, M )
            CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, -ONE,
     $                  DWORK, M, IR, LDST, ONE, DWORK(M*M+1), M )
            DSCALE = ZERO
            DSUM = ONE
            CALL DLASSQ( M*M, DWORK(M*M+1), 1, DSCALE, DSUM )
C
            CALL DLACPY( 'All', M, M, B, LDB, DWORK(M*M+1), M )
            CALL DGEMM( 'Transpose', 'No Transpose', M, M, M, ONE,
     $                  IR, LDST, T, LDST, ZERO, DWORK, M )
            CALL DGEMM( 'No Transpose', 'Transpose', M, M, M, -ONE,
     $                  DWORK, M, LI, LDST, ONE, DWORK(M*M+1), M )
            CALL DLASSQ( M*M, DWORK(M*M+1), 1, DSCALE, DSUM )
            SS = DSCALE*SQRT( DSUM )
            DTRONG = ( SS.LE.THRESH )
            IF( .NOT.DTRONG )
     $         GO TO 50
C
         END IF
C
C        If the swap is accepted ("weakly" and "strongly"), apply the
C        transformations and set N1-by-N2 (2,1)-block to zero.
C
         CALL DLASET( 'All', N1, N2, ZERO, ZERO, S(N2+1,1), LDST )
C
C        Copy (S,T) to (A,B).
C
         CALL DLACPY( 'All', M, M, S, LDST, A, LDA )
         CALL DLACPY( 'All', M, M, T, LDST, B, LDB )
         CALL DLASET( 'All', LDST, LDST, ZERO, ZERO, T, LDST )
C
C        Standardize existing 2-by-2 blocks.
C
         CALL DLASET( 'All', M, M, ZERO, ZERO, DWORK, M )
         DWORK(1) = ONE
         T(1,1) = ONE
         IF ( N2.GT.1 ) THEN
            CALL MB03YT( A, LDA, B, LDB, AR, AI, BE, DWORK(1), DWORK(2),
     $                   T(1,1), T(2,1) )
            DWORK(M+1) = -DWORK(2)
            DWORK(M+2) = DWORK(1)
            T(N2,N2) = T(1,1)
            T(1,2) = -T(2,1)
         END IF
         DWORK(M*M) = ONE
         T(M,M) = ONE
C
         IF ( N1.GT.1 ) THEN
            CALL MB03YT( A(N2+1,N2+1), LDA, B(N2+1,N2+1), LDB, TAUR,
     $                   TAUL, DWORK(M*M+1), DWORK(N2*M+N2+1),
     $                   DWORK(N2*M+N2+2), T(N2+1,N2+1), T(M,M-1) )
            DWORK(M*M) = DWORK(N2*M+N2+1)
            DWORK(M*M-1 ) = -DWORK(N2*M+N2+2)
            T(M,M) = T(N2+1,N2+1)
            T(M-1,M) = -T(M,M-1)
         END IF
C
         CALL DGEMM( 'Transpose', 'No Transpose', N2, N1, N2, ONE,
     $               DWORK, M, A(1,N2+1), LDA, ZERO, DWORK(M*M+1), N2 )
         CALL DLACPY( 'All', N2, N1, DWORK(M*M+1), N2, A(1,N2+1), LDA )
         CALL DGEMM( 'Transpose', 'No Transpose', N2, N1, N2, ONE,
     $               T(1,1), LDST, B(1,N2+1), LDB, ZERO,
     $               DWORK(M*M+1), N2 )
         CALL DLACPY( 'All', N2, N1, DWORK(M*M+1), N2, B(1,N2+1), LDB )
         CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE, LI,
     $               LDST, DWORK, M, ZERO, DWORK(M*M+1), M )
         CALL DLACPY( 'All', M, M, DWORK(M*M+1), M, LI, LDST )
         CALL DGEMM( 'No Transpose', 'No Transpose', N2, N1, N1, ONE,
     $               A(1,N2+1), LDA, T(N2+1,N2+1), LDST, ZERO,
     $               DWORK(M*M+1), M )
         CALL DLACPY( 'All', N2, N1, DWORK(M*M+1), M, A(1,N2+1), LDA )
         CALL DGEMM( 'No Transpose', 'No Transpose', N2, N1, N1, ONE,
     $               B(1,N2+1), LDB, DWORK(N2*M+N2+1), M, ZERO,
     $               DWORK(M*M+1), M )
         CALL DLACPY( 'All', N2, N1, DWORK(M*M+1), M, B(1,N2+1), LDB )
         CALL DGEMM( 'Transpose', 'No Transpose', M, M, M, ONE, T,
     $               LDST, IR, LDST, ZERO, DWORK, M )
         CALL DLACPY( 'All', M, M, DWORK, M, IR, LDST )
C
C        Accumulate transformations into Q and Z if requested.
C
         IF( WANTQ ) THEN
            CALL DGEMM( 'No Transpose', 'No Transpose', M, M, M, ONE, Q,
     $                  LDQ, LI, LDST, ZERO, DWORK, M )
            CALL DLACPY( 'All', M, M, DWORK, M, Q, LDQ )
         END IF
C
         IF( WANTZ ) THEN
            CALL DGEMM( 'No Transpose', 'Transpose', M, M, M, ONE, Z,
     $                  LDZ, IR, LDST, ZERO, DWORK, M )
            CALL DLACPY( 'Full', M, M, DWORK, M, Z, LDZ )
C
         END IF
C
C        Exit with INFO = 0 if swap was successfully performed.
C
         RETURN
C
      END IF
C
C     Exit with INFO = 1 if swap was rejected.
C
   50 CONTINUE
C
      INFO = 1
      RETURN
C *** Last line of MB03WA ***
      END
