      SUBROUTINE AB05RD( FBTYPE, JOBD, N, M, P, MV, PZ, ALPHA, BETA, A,
     $                   LDA, B, LDB, C, LDC, D, LDD, F, LDF, K, LDK,
     $                   G, LDG, H, LDH, RCOND, BC, LDBC, CC, LDCC,
     $                   DC, LDDC, IWORK, DWORK, LDWORK, INFO )
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
C     To construct for a given state space system (A,B,C,D) the closed-
C     loop system (Ac,Bc,Cc,Dc) corresponding to the mixed output and
C     state feedback control law
C
C          u = alpha*F*y + beta*K*x + G*v
C          z = H*y.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     FBTYPE  CHARACTER*1
C             Specifies the type of the feedback law as follows:
C             = 'I':  Unitary output feedback (F = I);
C             = 'O':  General output feedback.
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears
C             in the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of state vector x, i.e. the order of the
C             matrix A, the number of rows of B and the number of
C             columns of C.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of input vector u, i.e. the number of
C             columns of matrices B and D, and the number of rows of F.
C             M >= 0.
C
C     P       (input) INTEGER
C             The dimension of output vector y, i.e. the number of rows
C             of matrices C and D, and the number of columns of F.
C             P >= 0 and P = M if FBTYPE = 'I'.
C
C     MV      (input) INTEGER
C             The dimension of the new input vector v, i.e. the number
C             of columns of matrix G.  MV >= 0.
C
C     PZ      (input) INTEGER.
C             The dimension of the new output vector z, i.e. the number
C             of rows of matrix H.  PZ >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The coefficient alpha in the output feedback law.
C
C     BETA    (input) DOUBLE PRECISION.
C             The coefficient beta in the state feedback law.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the system state transition matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the state matrix Ac of the closed-loop system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the system input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the intermediary input matrix B1 (see METHOD).
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the system output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the intermediary output matrix C1+BETA*D1*K (see METHOD).
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P) if N > 0.
C             LDC >= 1 if N = 0.
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, if JOBD = 'D', the leading P-by-M part of this
C             array must contain the system direct input/output
C             transmission matrix D.
C             On exit, the leading P-by-M part of this array contains
C             the intermediary direct input/output transmission matrix
C             D1 (see METHOD).
C             The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P) if JOBD = 'D'.
C             LDD >= 1 if JOBD = 'Z'.
C
C     F       (input) DOUBLE PRECISION array, dimension (LDF,P)
C             If FBTYPE = 'O', the leading M-by-P part of this array
C             must contain the output feedback matrix F.
C             If FBTYPE = 'I', then the feedback matrix is assumed to be
C             an M x M order identity matrix.
C             The array F is not referenced if FBTYPE = 'I'  or
C             ALPHA = 0.
C
C     LDF     INTEGER
C             The leading dimension of array F.
C             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0.
C             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0.
C
C     K       (input) DOUBLE PRECISION array, dimension (LDK,N)
C             The leading M-by-N part of this array must contain the
C             state feedback matrix K.
C             The array K is not referenced if BETA = 0.
C
C     LDK     INTEGER
C             The leading dimension of the array K.
C             LDK >= MAX(1,M) if BETA <> 0.
C             LDK >= 1 if BETA = 0.
C
C     G       (input) DOUBLE PRECISION array, dimension (LDG,MV)
C             The leading M-by-MV part of this array must contain the
C             system input scaling matrix G.
C
C     LDG     INTEGER
C             The leading dimension of the array G.  LDG >= MAX(1,M).
C
C     H       (input) DOUBLE PRECISION array, dimension (LDH,P)
C             The leading PZ-by-P part of this array must contain the
C             system output scaling matrix H.
C
C     LDH     INTEGER
C             The leading dimension of the array H.  LDH >= MAX(1,PZ).
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal condition number of the matrix
C             I - alpha*D*F.
C
C     BC      (output) DOUBLE PRECISION array, dimension (LDBC,MV)
C             The leading N-by-MV part of this array contains the input
C             matrix Bc of the closed-loop system.
C
C     LDBC    INTEGER
C             The leading dimension of array BC.  LDBC >= MAX(1,N).
C
C     CC      (output) DOUBLE PRECISION array, dimension (LDCC,N)
C             The leading PZ-by-N part of this array contains the
C             system output matrix Cc of the closed-loop system.
C
C     LDCC    INTEGER
C             The leading dimension of array CC.
C             LDCC >= MAX(1,PZ) if N > 0.
C             LDCC >= 1 if N = 0.
C
C     DC      (output) DOUBLE PRECISION array, dimension (LDDC,MV)
C             If JOBD = 'D', the leading PZ-by-MV part of this array
C             contains the direct input/output transmission matrix Dc
C             of the closed-loop system.
C             The array DC is not referenced if JOBD = 'Z'.
C
C     LDDC    INTEGER
C             The leading dimension of array DC.
C             LDDC >= MAX(1,PZ) if JOBD = 'D'.
C             LDDC >= 1 if JOBD = 'Z'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK >= MAX(1,2*P) if JOBD = 'D'.
C             LIWORK >= 1 if JOBD = 'Z'.
C             IWORK is not referenced if JOBD = 'Z'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= wspace, where
C                   wspace = MAX( 1, M, P*MV, P*P + 4*P ) if JOBD = 'D',
C                   wspace = MAX( 1, M ) if JOBD = 'Z'.
C             For best performance, LDWORK >= MAX( wspace, N*M, N*P ).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix I - alpha*D*F is numerically singular.
C
C     METHOD
C
C     The matrices of the closed-loop system have the expressions:
C
C     Ac = A1 + beta*B1*K,      Bc = B1*G,
C     Cc = H*(C1 + beta*D1*K),  Dc = H*D1*G,
C
C     where
C
C     A1 = A + alpha*B*F*E*C,   B1 = B + alpha*B*F*E*D,
C     C1 = E*C,                 D1 = E*D,
C
C     with E = (I - alpha*D*F)**-1.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of computations basically depends on the conditioning
C     of the matrix I - alpha*D*F. If RCOND is very small, it is likely
C     that the computed results are inaccurate.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Research Establishment,
C     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven,
C     Belgium, Nov. 1996.
C
C     REVISIONS
C
C     January 14, 1997, February 18, 1998.
C     V. Sima, Research Institute for Informatics, Bucharest, July 2003,
C     Jan. 2005.
C
C     KEYWORDS
C
C     Multivariable system, state-space model, state-space
C     representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         FBTYPE, JOBD
      INTEGER           INFO, LDA, LDB, LDBC, LDC, LDCC, LDD, LDDC,
     $                  LDF, LDG, LDH, LDK, LDWORK, M, MV, N, P, PZ
      DOUBLE PRECISION  ALPHA, BETA, RCOND
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), BC(LDBC,*), C(LDC,*),
     $                  CC(LDCC,*), D(LDD,*), DC(LDDC,*), DWORK(*),
     $                  F(LDF,*), G(LDG,*), H(LDH,*), K(LDK,*)
C     .. Local Scalars ..
      LOGICAL           LJOBD, OUTPF, UNITF
      INTEGER           LDWP
C     .. External functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External subroutines ..
      EXTERNAL          AB05SD, DGEMM, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C
C     .. Executable Statements ..
C
C     Check the input scalar arguments.
C
      UNITF = LSAME( FBTYPE, 'I' )
      OUTPF = LSAME( FBTYPE, 'O' )
      LJOBD = LSAME( JOBD, 'D' )
C
      INFO = 0
C
      IF( .NOT.UNITF .AND. .NOT.OUTPF ) THEN
         INFO = -1
      ELSE IF( .NOT.LJOBD .AND. .NOT.LSAME( JOBD, 'Z' )  ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 .OR. UNITF.AND.P.NE.M ) THEN
         INFO = -5
      ELSE IF( MV.LT.0 ) THEN
         INFO = -6
      ELSE IF( PZ.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( ( N.GT.0 .AND. LDC.LT.MAX( 1, P ) ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -15
      ELSE IF( ( LJOBD .AND. LDD.LT.MAX( 1, P ) ) .OR.
     $    ( .NOT.LJOBD .AND. LDD.LT.1 ) ) THEN
         INFO = -17
      ELSE IF( ( OUTPF .AND. ALPHA.NE.ZERO .AND. LDF.LT.MAX( 1, M ) )
     $  .OR. ( ( UNITF .OR.  ALPHA.EQ.ZERO ) .AND. LDF.LT.1 ) ) THEN
         INFO = -19
      ELSE IF( ( BETA.NE.ZERO .AND. LDK.LT.MAX( 1, M ) ) .OR.
     $         ( BETA.EQ.ZERO .AND. LDK.LT.1 ) ) THEN
         INFO = -21
      ELSE IF( LDG.LT.MAX( 1, M ) ) THEN
         INFO = -23
      ELSE IF( LDH.LT.MAX( 1, PZ ) ) THEN
         INFO = -25
      ELSE IF( LDBC.LT.MAX( 1, N ) ) THEN
         INFO = -28
      ELSE IF( ( N.GT.0 .AND. LDCC.LT.MAX( 1, PZ ) ) .OR.
     $         ( N.EQ.0 .AND. LDCC.LT.1 ) ) THEN
         INFO = -30
      ELSE IF( ( ( LJOBD .AND. LDDC.LT.MAX( 1, PZ ) ) .OR.
     $      ( .NOT.LJOBD .AND. LDDC.LT.1 ) ) ) THEN
         INFO = -32
      ELSE IF( ( LJOBD .AND. LDWORK.LT.MAX( 1, M, P*MV, P*P + 4*P ) )
     $   .OR. (  .NOT.LJOBD .AND. LDWORK.LT.MAX( 1, M ) ) ) THEN
         INFO = -35
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB05RD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, MIN( M, P ), MIN( MV, PZ ) ).EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
C
C     Apply the partial output feedback u = alpha*F*y + v1
C
      CALL AB05SD( FBTYPE, JOBD, N, M, P, ALPHA, A, LDA, B, LDB, C,
     $             LDC, D, LDD, F, LDF, RCOND, IWORK, DWORK, LDWORK,
     $             INFO )
      IF ( INFO.NE.0 ) RETURN
C
C     Apply the partial state feedback v1 = beta*K*x + v2.
C
C     Compute Ac = A1 + beta*B1*K and C1 <- C1 + beta*D1*K.
C
      IF( BETA.NE.ZERO .AND. N.GT.0 ) THEN
         CALL DGEMM( 'N', 'N', N, N, M, BETA, B, LDB, K, LDK, ONE, A,
     $               LDA )
         IF( LJOBD )
     $       CALL DGEMM( 'N', 'N', P, N, M, BETA, D, LDD, K, LDK, ONE,
     $                   C, LDC )
      END IF
C
C     Apply the input and output conversions v2 = G*v, z = H*y.
C
C     Compute Bc = B1*G.
C
      CALL DGEMM( 'N', 'N', N, MV, M, ONE, B, LDB, G, LDG, ZERO, BC,
     $            LDBC )
C
C     Compute Cc = H*C1.
C
      IF( N.GT.0 )
     $   CALL DGEMM( 'N', 'N', PZ, N, P, ONE, H, LDH, C, LDC, ZERO, CC,
     $               LDCC )
C
C     Compute Dc = H*D1*G.
C
      IF( LJOBD ) THEN
         LDWP = MAX( 1, P )
         CALL DGEMM( 'N', 'N', P, MV, M, ONE, D, LDD, G, LDG, ZERO,
     $               DWORK, LDWP )
         CALL DGEMM( 'N', 'N', PZ, MV, P, ONE, H, LDH, DWORK, LDWP,
     $               ZERO, DC, LDDC )
      END IF
C
      RETURN
C *** Last line of AB05RD ***
      END
