      SUBROUTINE SB10ED( N, M, NP, NCON, NMEAS, A, LDA, B, LDB, C, LDC,
     $                   D, LDD, AK, LDAK, BK, LDBK, CK, LDCK, DK, LDDK,
     $                   RCOND, TOL, IWORK, DWORK, LDWORK, BWORK, INFO )
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
C     To compute the matrices of the H2 optimal n-state controller
C
C                           | AK | BK |
C                       K = |----|----|
C                           | CK | DK |
C
C     for the discrete-time system
C
C                   | A  | B1  B2  |   | A | B |
C               P = |----|---------| = |---|---| ,
C                   | C1 |  0  D12 |   | C | D |
C                   | C2 | D21 D22 |
C
C     where B2 has as column size the number of control inputs (NCON)
C     and C2 has as row size the number of measurements (NMEAS) being
C     provided to the controller.
C
C     It is assumed that
C
C     (A1) (A,B2) is stabilizable and (C2,A) is detectable,
C
C     (A2) D12 is full column rank and D21 is full row rank,
C
C               j*Theta
C     (A3) | A-e       *I  B2  | has full column rank for all
C          |    C1         D12 |
C
C          0 <= Theta < 2*Pi ,
C
C
C               j*Theta
C     (A4) | A-e       *I  B1  |  has full row rank for all
C          |    C2         D21 |
C
C          0 <= Theta < 2*Pi .
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     NP      (input) INTEGER
C             The row size of the matrix C.  NP >= 0.
C
C     NCON    (input) INTEGER
C             The number of control inputs (M2).  M >= NCON >= 0,
C             NP-NMEAS >= NCON.
C
C     NMEAS   (input) INTEGER
C             The number of measurements (NP2).  NP >= NMEAS >= 0,
C             M-NCON >= NMEAS.
C
C     A       (input/worksp.) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             system state matrix A.
C             This array is modified internally, but it is restored on
C             exit.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading NP-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,NP).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading NP-by-M part of this array must contain the
C             system input/output matrix D.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= max(1,NP).
C
C     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N)
C             The leading N-by-N part of this array contains the
C             controller state matrix AK.
C
C     LDAK    INTEGER
C             The leading dimension of the array AK.  LDAK >= max(1,N).
C
C     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS)
C             The leading N-by-NMEAS part of this array contains the
C             controller input matrix BK.
C
C     LDBK    INTEGER
C             The leading dimension of the array BK.  LDBK >= max(1,N).
C
C     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N)
C             The leading NCON-by-N part of this array contains the
C             controller output matrix CK.
C
C     LDCK    INTEGER
C             The leading dimension of the array CK.
C             LDCK >= max(1,NCON).
C
C     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS)
C             The leading NCON-by-NMEAS part of this array contains the
C             controller input/output matrix DK.
C
C     LDDK    INTEGER
C             The leading dimension of the array DK.
C             LDDK >= max(1,NCON).
C
C     RCOND   (output) DOUBLE PRECISION array, dimension (7)
C             RCOND contains estimates the reciprocal condition
C             numbers of the matrices which are to be inverted and the
C             reciprocal condition numbers of the Riccati equations
C             which have to be solved during the computation of the
C             controller. (See the description of the algorithm in [2].)
C             RCOND(1) contains the reciprocal condition number of the
C                      control transformation matrix TU;
C             RCOND(2) contains the reciprocal condition number of the
C                      measurement transformation matrix TY;
C             RCOND(3) contains the reciprocal condition number of the
C                      matrix Im2 + B2'*X2*B2;
C             RCOND(4) contains the reciprocal condition number of the
C                      matrix Ip2 + C2*Y2*C2';
C             RCOND(5) contains the reciprocal condition number of the
C                      X-Riccati equation;
C             RCOND(6) contains the reciprocal condition number of the
C                      Y-Riccati equation;
C             RCOND(7) contains the reciprocal condition number of the
C                      matrix Im2 + DKHAT*D22 .
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             Tolerance used for controlling the accuracy of the
C             transformations applied for diagonalizing D12 and D21,
C             and for checking the nonsingularity of the matrices to be
C             inverted. If TOL <= 0, then a default value equal to
C             sqrt(EPS) is used, where EPS is the relative machine
C             precision.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension max(2*M2,2*N,N*N,NP2)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal
C             LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 +
C                       max(1,LW1,LW2,LW3,LW4,LW5,LW6), where
C             LW1 = (N+NP1+1)*(N+M2) + max(3*(N+M2)+N+NP1,5*(N+M2)),
C             LW2 = (N+NP2)*(N+M1+1) + max(3*(N+NP2)+N+M1,5*(N+NP2)),
C             LW3 = M2 + NP1*NP1 + max(NP1*max(N,M1),3*M2+NP1,5*M2),
C             LW4 = NP2 + M1*M1 + max(max(N,NP1)*M1,3*NP2+M1,5*NP2),
C             LW5 = 2*N*N+max(1,14*N*N+6*N+max(14*N+23,16*N),M2*(N+M2+
C                             max(3,M1)),NP2*(N+NP2+3)),
C             LW6 = max(N*M2,N*NP2,M2*NP2,M2*M2+4*M2),
C             with M1 = M - M2 and NP1 = NP - NP2.
C             For good performance, LDWORK must generally be larger.
C             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is
C             2*Q*(3*Q+2*N)+max(1,(N+Q)*(N+Q+6),Q*(Q+max(N,Q,5)+1),
C                     2*N*N+max(1,14*N*N+6*N+max(14*N+23,16*N),
C                               Q*(N+Q+max(Q,3)))).
C
C     BWORK   LOGICAL array, dimension (2*N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C                                      j*Theta
C             = 1:  if the matrix | A-e       *I  B2  | had not full
C                                 |      C1       D12 |
C                   column rank in respect to the tolerance EPS;
C                                      j*Theta
C             = 2:  if the matrix | A-e       *I  B1  |  had not full
C                                 |      C2       D21 |
C                   row rank in respect to the tolerance EPS;
C             = 3:  if the matrix D12 had not full column rank in
C                   respect to the tolerance TOL;
C             = 4:  if the matrix D21 had not full row rank in respect
C                   to the tolerance TOL;
C             = 5:  if the singular value decomposition (SVD) algorithm
C                   did not converge (when computing the SVD of one of
C                   the matrices |A-I  B2 |, |A-I  B1 |, D12 or D21).
C                                |C1   D12|  |C2   D21|
C             = 6:  if the X-Riccati equation was not solved
C                   successfully;
C             = 7:  if the matrix Im2 + B2'*X2*B2 is not positive
C                   definite, or it is numerically singular (with
C                   respect to the tolerance TOL);
C             = 8:  if the Y-Riccati equation was not solved
C                   successfully;
C             = 9:  if the matrix Ip2 + C2*Y2*C2' is not positive
C                   definite, or it is numerically singular (with
C                   respect to the tolerance TOL);
C             =10:  if the matrix Im2 + DKHAT*D22 is singular, or its
C                   estimated condition number is larger than or equal
C                   to 1/TOL.
C
C     METHOD
C
C     The routine implements the formulas given in [1].
C
C     REFERENCES
C
C     [1] Zhou, K., Doyle, J.C., and Glover, K.
C         Robust and Optimal Control.
C         Prentice-Hall, Upper Saddle River, NJ, 1996.
C
C     [2] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M.
C         Fortran 77 routines for Hinf and H2 design of linear
C         discrete-time control systems.
C         Report 99-8, Department of Engineering, Leicester University,
C         April 1999.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of the result depends on the condition numbers of the
C     matrices which are to be inverted and on the condition numbers of
C     the matrix Riccati equations which are to be solved in the
C     computation of the controller. (The corresponding reciprocal
C     condition numbers are given in the output array RCOND.)
C
C     CONTRIBUTORS
C
C     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 1999,
C     Sept. 1999, Feb. 2000, Nov. 2005.
C
C     KEYWORDS
C
C     Algebraic Riccati equation, H2 optimal control, optimal regulator,
C     robust control.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDAK, LDB, LDBK, LDC, LDCK, LDD,
     $                   LDDK, LDWORK, M, N, NCON, NMEAS, NP
      DOUBLE PRECISION   TOL
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AK( LDAK, * ), B( LDB, * ),
     $                   BK( LDBK, * ), C( LDC, * ), CK( LDCK, * ),
     $                   D( LDD, * ), DK( LDDK, * ), DWORK( * ),
     $                   RCOND( * )
      LOGICAL            BWORK( * )
C     ..
C     .. Local Scalars ..
      INTEGER            I, INFO2, IWC, IWD, IWRK, IWTU, IWTY, IWX, IWY,
     $                   LW1, LW2, LW3, LW4, LW5, LW6, LWAMAX, M1, M2,
     $                   M2L, MINWRK, NL, NLP, NP1, NP2, NPL
      DOUBLE PRECISION   TOLL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLACPY, SB10PD, SB10SD, SB10TD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
C     ..
C     .. Executable Statements ..
C
C     Decode and Test input parameters.
C
      M1  = M - NCON
      M2  = NCON
      NP1 = NP - NMEAS
      NP2 = NMEAS
      NL  = MAX( 1, N )
      NPL = MAX( 1, NP )
      M2L = MAX( 1, M2 )
      NLP = MAX( 1, NP2 )
C
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( NP.LT.0 ) THEN
         INFO = -3
      ELSE IF( NCON.LT.0 .OR. M1.LT.0 .OR. M2.GT.NP1 ) THEN
         INFO = -4
      ELSE IF( NMEAS.LT.0 .OR. NP1.LT.0 .OR. NP2.GT.M1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.NL ) THEN
         INFO = -7
      ELSE IF( LDB.LT.NL ) THEN
         INFO = -9
      ELSE IF( LDC.LT.NPL ) THEN
         INFO = -11
      ELSE IF( LDD.LT.NPL ) THEN
         INFO = -13
      ELSE IF( LDAK.LT.NL ) THEN
         INFO = -15
      ELSE IF( LDBK.LT.NL ) THEN
         INFO = -17
      ELSE IF( LDCK.LT.M2L ) THEN
         INFO = -19
      ELSE IF( LDDK.LT.M2L ) THEN
         INFO = -21
      ELSE
C
C        Compute workspace.
C
         LW1 = ( N + NP1 + 1 )*( N + M2 ) + MAX( 3*( N + M2 ) + N + NP1,
     $                                           5*( N + M2 ) )
         LW2 = ( N + NP2 )*( N + M1 + 1 ) + MAX( 3*( N + NP2 ) + N +
     $                                           M1, 5*( N + NP2 ) )
         LW3 = M2 + NP1*NP1 + MAX( NP1*MAX( N, M1 ), 3*M2 + NP1, 5*M2 )
         LW4 = NP2 + M1*M1  + MAX( MAX( N, NP1 )*M1, 3*NP2 + M1, 5*NP2 )
         LW5 = 2*N*N + MAX( 1, 14*N*N +
     $                         6*N + MAX( 14*N + 23, 16*N ),
     $                         M2*( N + M2 +  MAX( 3, M1 ) ),
     $                         NP2*( N + NP2 + 3 ) )
         LW6 = MAX( N*M2, N*NP2, M2*NP2, M2*M2 + 4*M2 )
         MINWRK = N*M + NP*( N + M ) + M2*M2 + NP2*NP2 +
     $            MAX( 1, LW1, LW2, LW3, LW4, LW5, LW6 )
         IF( LDWORK.LT.MINWRK )
     $      INFO = -26
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SB10ED', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 .AND. MAX( M2, NP2 ).EQ.0 ) THEN
          RCOND( 1 ) = ONE
          RCOND( 2 ) = ONE
          RCOND( 3 ) = ONE
          RCOND( 4 ) = ONE
          RCOND( 5 ) = ONE
          RCOND( 6 ) = ONE
          RCOND( 7 ) = ONE
          DWORK( 1 ) = ONE
          RETURN
      END IF
C
      TOLL = TOL
      IF( TOLL.LE.ZERO ) THEN
C
C        Set the default value of the tolerance for rank tests.
C
         TOLL = SQRT( DLAMCH( 'Epsilon' ) )
      END IF
C
C     Workspace usage.
C
      IWC  = N*M  + 1
      IWD  = IWC  + NP*N
      IWTU = IWD  + NP*M
      IWTY = IWTU + M2*M2
      IWRK = IWTY + NP2*NP2
C
      CALL DLACPY( 'Full', N,  M, B, LDB, DWORK, NL )
      CALL DLACPY( 'Full', NP, N, C, LDC, DWORK( IWC ), NPL )
      CALL DLACPY( 'Full', NP, M, D, LDD, DWORK( IWD ), NPL )
C
C     Transform the system so that D12 and D21 satisfy the formulas
C     in the computation of the H2 optimal controller.
C     Since SLICOT Library routine SB10PD performs the tests
C     corresponding to the continuous-time counterparts of the
C     assumptions (A3) and (A4), for the frequency w = 0, the
C     next SB10PD routine call uses A - I.
C
      DO 10 I = 1, N
         A(I,I) = A(I,I) - ONE
   10 CONTINUE
C
      CALL SB10PD( N, M, NP, NCON, NMEAS, A, LDA, DWORK, NL,
     $             DWORK( IWC ), NPL, DWORK( IWD ), NPL, DWORK( IWTU ),
     $             M2L, DWORK( IWTY ), NLP, RCOND, TOLL, DWORK( IWRK ),
     $             LDWORK-IWRK+1, INFO2 )
C
      DO 20 I = 1, N
         A(I,I) = A(I,I) + ONE
   20 CONTINUE
C
      IF( INFO2.GT.0 ) THEN
         INFO = INFO2
         RETURN
      END IF
      LWAMAX = INT( DWORK( IWRK ) ) + IWRK - 1
C
      IWX  = IWRK
      IWY  = IWX + N*N
      IWRK = IWY + N*N
C
C     Compute the optimal H2 controller for the normalized system.
C
      CALL SB10SD( N, M, NP, NCON, NMEAS, A, LDA, DWORK, NL,
     $             DWORK( IWC ), NPL, DWORK( IWD ), NPL, AK, LDAK, BK,
     $             LDBK, CK, LDCK, DK, LDDK, DWORK( IWX ), NL,
     $             DWORK( IWY ), NL, RCOND( 3 ), TOLL, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, BWORK, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = INFO2 + 5
         RETURN
      END IF
      LWAMAX = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, LWAMAX )
C
      IWRK = IWX
C
C     Compute the H2 optimal controller for the original system.
C
      CALL SB10TD( N, M, NP, NCON, NMEAS, DWORK( IWD ), NPL,
     $             DWORK( IWTU ), M2L, DWORK( IWTY ), NLP, AK, LDAK, BK,
     $             LDBK, CK, LDCK, DK, LDDK, RCOND( 7 ), TOLL, IWORK,
     $             DWORK( IWRK ), LDWORK-IWRK+1, INFO2 )
      IF( INFO2.GT.0 ) THEN
         INFO = 10
         RETURN
      END IF
C
      DWORK( 1 ) = DBLE( LWAMAX )
      RETURN
C *** Last line of SB10ED ***
      END
