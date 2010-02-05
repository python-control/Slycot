      SUBROUTINE SB10ZP( DISCFL, N, A, LDA, B, C, D, IWORK, DWORK,
     $                   LDWORK, INFO )
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
C     To transform a SISO (single-input single-output) system [A,B;C,D]
C     by mirroring its unstable poles and zeros in the boundary of the
C     stability domain, thus preserving the frequency response of the
C     system, but making it stable and minimum phase. Specifically, for
C     a continuous-time system, the positive real parts of its poles
C     and zeros are exchanged with their negatives. Discrete-time
C     systems are first converted to continuous-time systems using a
C     bilinear transformation, and finally converted back.
C
C     ARGUMENTS
C
C     Input/Output parameters
C
C     DISCFL  (input) INTEGER
C             Indicates the type of the system, as follows:
C             = 0: continuous-time system;
C             = 1: discrete-time system.
C
C     N       (input/output) INTEGER
C             On entry, the order of the original system.  N >= 0.
C             On exit, the order of the transformed, minimal system.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original system matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the transformed matrix A, in an upper Hessenberg form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the original system
C             vector B.
C             On exit, this array contains the transformed vector B.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the original system
C             vector C.
C             On exit, this array contains the transformed vector C.
C             The first N-1 elements are zero (for the exit value of N).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (1)
C             On entry, this array must contain the original system
C             scalar D.
C             On exit, this array contains the transformed scalar D.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(2,N+1)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= max(N*N + 5*N, 6*N + 1 + min(1,N)).
C             For optimum performance LDWORK should be larger.
C
C     Error indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the discrete --> continuous transformation cannot
C                   be made;
C             = 2:  if the system poles cannot be found;
C             = 3:  if the inverse system cannot be found, i.e., D is
C                   (close to) zero;
C             = 4:  if the system zeros cannot be found;
C             = 5:  if the state-space representation of the new
C                   transfer function T(s) cannot be found;
C             = 6:  if the continuous --> discrete transformation cannot
C                   be made.
C
C     METHOD
C
C     First, if the system is discrete-time, it is transformed to
C     continuous-time using alpha = beta = 1 in the bilinear
C     transformation implemented in the SLICOT routine AB04MD.
C     Then the eigenvalues of A, i.e., the system poles, are found.
C     Then, the inverse of the original system is found and its poles,
C     i.e., the system zeros, are evaluated.
C     The obtained system poles Pi and zeros Zi are checked and if a
C     positive real part is detected, it is exchanged by -Pi or -Zi.
C     Then the polynomial coefficients of the transfer function
C     T(s) = Q(s)/P(s) are found.
C     The state-space representation of T(s) is then obtained.
C     The system matrices B, C, D are scaled so that the transformed
C     system has the same system gain as the original system.
C     If the original system is discrete-time, then the result (which is
C     continuous-time) is converted back to discrete-time.
C
C     CONTRIBUTORS
C
C     Asparuh Markovski, Technical University of Sofia, July 2003.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003.
C
C     KEYWORDS
C
C     Bilinear transformation, stability, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            DISCFL, INFO, LDA, LDWORK, N
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( * ), C( * ), D( * ), DWORK( * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, IDW1, IDW2, IDW3, IMP, IMZ, INFO2, IWA, IWP,
     $                   IWPS, IWQ, IWQS, LDW1, MAXWRK, REP, REZ
      DOUBLE PRECISION   RCOND, SCALB, SCALC, SCALD
C     ..
C     .. Local Arrays ..
      INTEGER            INDEX(1)
C     ..
C     .. External Subroutines ..
      EXTERNAL           AB04MD, AB07ND, DCOPY, DGEEV, DLACPY, DSCAL,
     $                   MC01PD, TD04AD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX, MIN, SIGN, SQRT
C
C     Test input parameters and workspace.
C
      INFO = 0
      IF ( DISCFL.NE.0 .AND. DISCFL.NE.1 ) THEN
         INFO = -1
      ELSE IF ( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDWORK.LT.MAX( N*N + 5*N, 6*N + 1 + MIN( 1, N ) ) ) THEN
         INFO = -10
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'SB10ZP', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Workspace usage 1.
C
      REP  = 1
      IMP  = REP + N
      REZ  = IMP + N
      IMZ  = REZ + N
      IWA  = REZ
      IDW1 = IWA + N*N
      LDW1 = LDWORK - IDW1 + 1
C
C     1. Discrete --> continuous transformation if needed.
C
      IF ( DISCFL.EQ.1 ) THEN
C
C        Workspace:  need    max(1,N);
C                            prefer  larger.
C
         CALL AB04MD( 'D', N, 1, 1, ONE, ONE, A, LDA, B, LDA, C, 1,
     $                D, 1, IWORK, DWORK, LDWORK, INFO2 )
         IF ( INFO2.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
         MAXWRK = INT( DWORK(1) )
      ELSE
         MAXWRK = 0
      END IF
C
C     2. Determine the factors for restoring system gain.
C
      SCALD = D(1)
      SCALC = SQRT( ABS( SCALD ) )
      SCALB = SIGN( SCALC, SCALD )
C
C     3. Find the system poles, i.e., the eigenvalues of A.
C        Workspace:  need    N*N + 2*N + 3*N;
C                            prefer  larger.
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK(IWA), N )
C
      CALL DGEEV( 'N', 'N', N, DWORK(IWA), N, DWORK(REP), DWORK(IMP),
     $            DWORK(IDW1), 1, DWORK(IDW1), 1, DWORK(IDW1), LDW1,
     $            INFO2 )
      IF ( INFO2.NE.0 ) THEN
         INFO = 2
         RETURN
      END IF
      MAXWRK = MAX( MAXWRK, INT( DWORK(IDW1) + IDW1 - 1 ) )
C
C     4. Compute the inverse system [Ai, Bi; Ci, Di].
C        Workspace:  need    N*N + 2*N + 4;
C                            prefer  larger.
C
      CALL AB07ND( N, 1, A, LDA, B, LDA, C, 1, D, 1, RCOND, IWORK,
     $             DWORK(IDW1), LDW1, INFO2 )
      IF ( INFO2.NE.0 ) THEN
         INFO = 3
         RETURN
      END IF
      MAXWRK = MAX( MAXWRK, INT( DWORK(IDW1) + IDW1 - 1 ) )
C
C     5. Find the system zeros, i.e., the eigenvalues of Ai.
C        Workspace:  need    4*N + 3*N;
C                            prefer  larger.
C
      IDW1 = IMZ + N
      LDW1 = LDWORK - IDW1 + 1
C
      CALL DGEEV( 'N', 'N', N, A, LDA, DWORK(REZ), DWORK(IMZ),
     $            DWORK(IDW1), 1, DWORK(IDW1), 1, DWORK(IDW1), LDW1,
     $            INFO2 )
      IF ( INFO2.NE.0 ) THEN
         INFO = 4
         RETURN
      END IF
      MAXWRK = MAX( MAXWRK, INT( DWORK(IDW1) + IDW1 - 1 ) )
C
C     6. Exchange the zeros and the poles with positive real parts with
C        their negatives.
C
      DO 10 I = 0, N - 1
         IF ( DWORK(REP+I).GT.ZERO )
     $      DWORK(REP+I) = -DWORK(REP+I)
         IF ( DWORK(REZ+I).GT.ZERO )
     $      DWORK(REZ+I) = -DWORK(REZ+I)
   10 CONTINUE
C
C     Workspace usage 2.
C
      IWP  = IDW1
      IDW2 = IWP + N + 1
      IWPS = 1
C
C     7. Construct the nominator and the denominator
C        of the system transfer function T( s ) = Q( s )/P( s ).
C     8. Rearrange the coefficients in Q(s) and P(s) because
C        MC01PD subroutine produces them in increasing powers of s.
C        Workspace:  need    6*N + 2.
C
      CALL MC01PD( N, DWORK(REP), DWORK(IMP), DWORK(IWP), DWORK(IDW2),
     $             INFO2 )
      CALL DCOPY( N+1, DWORK(IWP), -1, DWORK(IWPS), 1 )
C
C     Workspace usage 3.
C
      IWQ  = IDW1
      IWQS = IWPS + N + 1
      IDW3 = IWQS + N + 1
C
      CALL MC01PD( N, DWORK(REZ), DWORK(IMZ), DWORK(IWQ), DWORK(IDW2),
     $             INFO2 )
      CALL DCOPY( N+1, DWORK(IWQ), -1, DWORK(IWQS), 1 )
C
C     9. Make the conversion T(s) --> [A, B; C, D].
C        Workspace:  need    2*N + 2 + N + max(N,3);
C                            prefer  larger.
C
      INDEX(1) = N
      CALL TD04AD( 'R', 1, 1, INDEX, DWORK(IWPS), 1, DWORK(IWQS), 1, 1,
     $             N, A, LDA, B, LDA, C, 1, D, 1, -ONE, IWORK,
     $             DWORK(IDW3), LDWORK-IDW3+1, INFO2 )
      IF ( INFO2.NE.0 ) THEN
         INFO = 5
         RETURN
      END IF
      MAXWRK = MAX( MAXWRK, INT( DWORK(IDW3) + IDW3 - 1 ) )
C
C    10. Scale the transformed system to the previous gain.
C
      IF ( N.GT.0 ) THEN
         CALL DSCAL( N, SCALB, B, 1 )
         C(N) = SCALC*C(N)
      END IF
C
      D(1) = SCALD
C
C     11. Continuous --> discrete transformation if needed.
C
      IF ( DISCFL.EQ.1 ) THEN
         CALL AB04MD( 'C', N, 1, 1, ONE, ONE, A, LDA, B, LDA, C, 1,
     $                D, 1, IWORK, DWORK, LDWORK, INFO2 )

         IF ( INFO2.NE.0 ) THEN
            INFO = 6
            RETURN
         END IF
      END IF
C
      DWORK(1) = MAXWRK
      RETURN
C
C *** Last line of SB10ZP ***
      END
