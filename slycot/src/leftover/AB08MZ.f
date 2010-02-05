      SUBROUTINE AB08MZ( EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   RANK, TOL, IWORK, DWORK, ZWORK, LZWORK, INFO )
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
C     To compute the normal rank of the transfer-function matrix of a
C     state-space model (A,B,C,D).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to balance the compound
C             matrix (see METHOD) as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables, i.e., the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input) COMPLEX*16 array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) COMPLEX*16 array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B of the system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) COMPLEX*16 array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C of the system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) COMPLEX*16 array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             direct transmission matrix D of the system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     RANK    (output) INTEGER
C             The normal rank of the transfer-function matrix.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (2*N+MAX(M,P)+1)
C
C     DWORK   DOUBLE PRECISION array, dimension (2*MAX(M,P))
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) returns the optimal value
C             of LZWORK.
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK >= (N+P)*(N+M) + MAX(MIN(P,M) + MAX(3*M-1,N), 1,
C                                         MIN(P,N) + MAX(3*P-1,N+P,N+M))
C             For optimum performance LZWORK should be larger.
C
C             If LZWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             ZWORK array, returns this value as the first entry of
C             the ZWORK array, and no error message related to LZWORK
C             is issued by XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The routine reduces the (N+P)-by-(M+N) compound matrix (B  A)
C                                                            (D  C)
C
C     to one with the same invariant zeros and with D of full row rank.
C     The normal rank of the transfer-function matrix is the rank of D.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable (see [2] and [1]).
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C     Complex version: V. Sima, Research Institute for Informatics,
C     Bucharest, Dec. 2008.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2009,
C     Apr. 2009.
C
C     KEYWORDS
C
C     Multivariable system, unitary transformation,
C     structural invariant.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL
      INTEGER           INFO, LDA, LDB, LDC, LDD, LZWORK, M, N, P, RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), ZWORK(*)
      DOUBLE PRECISION  DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LQUERY
      INTEGER           I, KW, MU, NB, NINFZ, NKROL, NM, NP, NU, RO,
     $                  SIGMA, WRKOPT
      DOUBLE PRECISION  MAXRED, SVLMAX, THRESH, TOLER
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, ZLANGE
      EXTERNAL          DLAMCH, LSAME, ZLANGE
C     .. External Subroutines ..
      EXTERNAL          AB8NXZ, TB01IZ, XERBLA, ZLACPY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      NP = N + P
      NM = N + M
      INFO = 0
      LEQUIL = LSAME( EQUIL, 'S' )
      LQUERY = ( LZWORK.EQ.-1 )
      WRKOPT = NP*NM
C
C     Test the input scalar arguments.
C
      IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE
         KW = WRKOPT + MAX( MIN( P, M ) + MAX( 3*M-1, N ), 1,
     $                      MIN( P, N ) + MAX( 3*P-1, NP, NM ) )
         IF( LQUERY ) THEN
            SVLMAX = ZERO
            NINFZ  = 0
            CALL AB8NXZ( N, M, P, P, 0, SVLMAX, ZWORK, MAX( 1, NP ),
     $                   NINFZ, IWORK, IWORK, MU, NU, NKROL, TOL, IWORK,
     $                   DWORK, ZWORK, -1, INFO )
            WRKOPT = MAX( KW, WRKOPT + INT( ZWORK(1) ) )
         ELSE IF( LZWORK.LT.KW ) THEN
            INFO = -17
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB08MZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         ZWORK(1) = WRKOPT
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( M, P ).EQ.0 ) THEN
         RANK = 0
         ZWORK(1) = ONE
         RETURN
      END IF
C
      DO 10 I = 1, 2*N+1
         IWORK(I) = 0
   10 CONTINUE
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of workspace needed at that point in the code,
C     as well as the preferred amount for good performance.)
C
C     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N).
C                                    ( D  C )
C     Complex workspace: need   (N+P)*(N+M).
C
      CALL ZLACPY( 'Full', N, M, B, LDB, ZWORK, NP )
      CALL ZLACPY( 'Full', P, M, D, LDD, ZWORK(N+1), NP )
      CALL ZLACPY( 'Full', N, N, A, LDA, ZWORK(NP*M+1), NP )
      CALL ZLACPY( 'Full', P, N, C, LDC, ZWORK(NP*M+N+1), NP )
C
C     If required, balance the compound matrix (default MAXRED).
C     Real Workspace: need   N.
C
      KW = WRKOPT + 1
      IF ( LEQUIL ) THEN
         MAXRED = ZERO
         CALL TB01IZ( 'A', N, M, P, MAXRED, ZWORK(NP*M+1), NP, ZWORK,
     $                NP, ZWORK(NP*M+N+1), NP, DWORK, INFO )
      END IF
C
C     If required, set tolerance.
C
      THRESH = SQRT( DBLE( NP*NM ) )*DLAMCH( 'Precision' )
      TOLER = TOL
      IF ( TOLER.LT.THRESH ) TOLER = THRESH
      SVLMAX = ZLANGE( 'Frobenius', NP, NM, ZWORK, NP, DWORK )
C
C     Reduce this system to one with the same invariant zeros and with
C     D full row rank MU (the normal rank of the original system).
C     Real workspace:    need   2*MAX(M,P);
C     Complex workspace: need   (N+P)*(N+M) +
C                               MAX( 1, MIN(P,M) + MAX(3*M-1,N),
C                                       MIN(P,N) + MAX(3*P-1,N+P,N+M) );
C                        prefer larger.
C     Integer workspace: 2*N+MAX(M,P)+1.
C
      RO = P
      SIGMA = 0
      NINFZ = 0
      CALL AB8NXZ( N, M, P, RO, SIGMA, SVLMAX, ZWORK, NP, NINFZ, IWORK,
     $             IWORK(N+1), MU, NU, NKROL, TOLER, IWORK(2*N+2),
     $             DWORK, ZWORK(KW), LZWORK-KW+1, INFO )
      RANK = MU
C
      ZWORK(1) = MAX( WRKOPT, INT( ZWORK(KW) ) + KW - 1 )
      RETURN
C *** Last line of AB08MZ ***
      END
