      SUBROUTINE AB05SD( FBTYPE, JOBD, N, M, P, ALPHA, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, F, LDF, RCOND, IWORK, DWORK,
     $                   LDWORK, INFO)
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
C     loop system (Ac,Bc,Cc,Dc) corresponding to the output feedback
C     control law
C
C          u = alpha*F*y + v.
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
C             Specifies whether or not a non-zero matrix D appears in
C             the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables, i.e. the order of the
C             matrix A, the number of rows of B and the number of
C             columns of C.  N >= 0.
C
C     M       (input) INTEGER
C             The number of input variables, i.e. the number of columns
C             of matrices B and D, and the number of rows of F.  M >= 0.
C
C     P       (input) INTEGER
C             The number of output variables, i.e. the number of rows of
C             matrices C and D, and the number of columns of F.  P >= 0
C             and P = M if FBTYPE = 'I'.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The coefficient alpha in the output feedback law.
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
C             the input matrix Bc of the closed-loop system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the system output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the output matrix Cc of the closed-loop system.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,P) if N > 0.
C             LDC >= 1 if N = 0.
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the system direct input/output transmission
C             matrix D.
C             On exit, if JOBD = 'D', the leading P-by-M part of this
C             array contains the direct input/output transmission
C             matrix Dc of the closed-loop system.
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
C             The array F is not referenced if FBTYPE = 'I' or
C             ALPHA = 0.
C
C     LDF     INTEGER
C             The leading dimension of array F.
C             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0.
C             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0.
C
C     RCOND   (output) DOUBLE PRECISION
C             The reciprocal condition number of the matrix
C             I - alpha*D*F.
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
C                       wspace = MAX( 1, M, P*P + 4*P ) if JOBD = 'D',
C                       wspace = MAX( 1, M ) if JOBD = 'Z'.
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
C     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D,
C     Cc = E*C,                Dc = E*D,
C
C     where E = (I - alpha*D*F)**-1.
C
C     NUMERICAL ASPECTS
C
C     The accuracy of computations basically depends on the conditioning
C     of the matrix I - alpha*D*F.  If RCOND is very small, it is likely
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
C     January 14, 1997.
C     V. Sima, Research Institute for Informatics, Bucharest, July 2003.
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
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDF, LDWORK, M, N, P
      DOUBLE PRECISION  ALPHA, RCOND
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), F(LDF,*)
C     .. Local Scalars ..
      LOGICAL           LJOBD, OUTPF, UNITF
      INTEGER           I, IW, LDWN, LDWP
      DOUBLE PRECISION  ENORM
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME
C     .. External subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGECON, DGEMM, DGEMV, DGETRF,
     $                  DGETRS, DLACPY, DLASCL, XERBLA
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
      LDWN = MAX( 1, N )
      LDWP = MAX( 1, P )
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
      ELSE IF( LDA.LT.LDWN ) THEN
         INFO = -7
      ELSE IF( LDB.LT.LDWN ) THEN
         INFO = -9
      ELSE IF( ( N.GT.0 .AND. LDC.LT.LDWP ) .OR.
     $         ( N.EQ.0 .AND. LDC.LT.1 ) ) THEN
         INFO = -11
      ELSE IF( ( LJOBD .AND. LDD.LT.LDWP ) .OR.
     $    ( .NOT.LJOBD .AND. LDD.LT.1 ) ) THEN
         INFO = -13
      ELSE IF( ( OUTPF .AND. ALPHA.NE.ZERO .AND. LDF.LT.MAX( 1, M ) )
     $  .OR. ( ( UNITF .OR.  ALPHA.EQ.ZERO ) .AND. LDF.LT.1 ) ) THEN
         INFO = -16
      ELSE IF( ( LJOBD .AND. LDWORK.LT.MAX( 1, M, P*P + 4*P ) ) .OR.
     $    ( .NOT.LJOBD .AND. LDWORK.LT.MAX( 1, M ) ) ) THEN
         INFO = -20
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB05SD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      RCOND = ONE
      IF ( MAX( N, MIN( M, P ) ).EQ.0 .OR. ALPHA.EQ.ZERO )
     $   RETURN
C
      IF (LJOBD) THEN
         IW = P*P + 1
C
C        Compute I - alpha*D*F.
C
         IF( UNITF) THEN
            CALL DLACPY( 'F', P, P, D, LDD, DWORK, LDWP )
            IF ( ALPHA.NE.-ONE )
     $         CALL DLASCL( 'G', 0, 0, ONE, -ALPHA, P, P, DWORK, LDWP,
     $                      INFO )
         ELSE
            CALL DGEMM( 'N', 'N', P, P, M, -ALPHA, D, LDD, F, LDF, ZERO,
     $                  DWORK, LDWP )
         END IF
C
         DUMMY(1) = ONE
         CALL DAXPY( P, ONE, DUMMY, 0, DWORK, P+1 )
C
C        Compute Cc = E*C, Dc = E*D, where E = (I - alpha*D*F)**-1.
C
         ENORM = DLANGE( '1', P, P, DWORK, LDWP, DWORK(IW) )
         CALL DGETRF( P, P, DWORK, LDWP, IWORK, INFO )
         IF( INFO.GT.0 ) THEN
C
C           Error return.
C
            RCOND = ZERO
            INFO = 1
            RETURN
         END IF
         CALL DGECON( '1', P, DWORK, LDWP, ENORM, RCOND, DWORK(IW),
     $                IWORK(P+1), INFO )
         IF( RCOND.LE.DLAMCH('E') ) THEN
C
C           Error return.
C
            INFO = 1
            RETURN
         END IF
C
         IF( N.GT.0 )
     $      CALL DGETRS( 'N', P, N, DWORK, LDWP, IWORK, C, LDC, INFO )
         CALL DGETRS( 'N', P, M, DWORK, LDWP, IWORK, D, LDD, INFO )
      END IF
C
      IF ( N.EQ.0 )
     $   RETURN
C
C     Compute Ac = A + alpha*B*F*Cc and Bc = B + alpha*B*F*Dc.
C
      IF( UNITF ) THEN
         CALL DGEMM( 'N', 'N', N, N, M, ALPHA, B, LDB, C, LDC, ONE, A,
     $               LDA )
         IF( LJOBD ) THEN
C
            IF( LDWORK.LT.N*M ) THEN
C
C              Not enough working space for using DGEMM.
C
               DO 10 I = 1, N
                  CALL DCOPY( P, B(I,1), LDB, DWORK, 1 )
                  CALL DGEMV( 'T', P, P, ALPHA, D, LDD, DWORK, 1, ONE,
     $                        B(I,1), LDB )
   10          CONTINUE
C
            ELSE
               CALL DLACPY( 'F', N, M, B, LDB, DWORK, LDWN )
               CALL DGEMM( 'N', 'N', N, P, M, ALPHA, DWORK, LDWN, D,
     $                     LDD, ONE, B, LDB )
            END IF
         END IF
      ELSE
C
         IF( LDWORK.LT.N*P ) THEN
C
C           Not enough working space for using DGEMM.
C
            DO 20 I = 1, N
               CALL DGEMV( 'N', M, P, ALPHA, F, LDF, C(1,I), 1, ZERO,
     $                     DWORK, 1 )
               CALL DGEMV( 'N', N, M, ONE, B, LDB, DWORK, 1, ONE,
     $                     A(1,I), 1 )
   20       CONTINUE
C
            IF( LJOBD ) THEN
C
               DO 30 I = 1, N
                  CALL DGEMV( 'T', M, P, ALPHA, F, LDF, B(I,1), LDB,
     $                        ZERO, DWORK, 1 )
                  CALL DGEMV( 'T', P, M, ONE, D, LDD, DWORK, 1, ONE,
     $                        B(I,1), LDB )
   30          CONTINUE
C
            END IF
         ELSE
C
            CALL DGEMM( 'N', 'N', N, P, M, ALPHA, B, LDB, F, LDF,
     $                  ZERO, DWORK, LDWN )
            CALL DGEMM( 'N', 'N', N, N, P, ONE, DWORK, LDWN, C, LDC,
     $                  ONE, A, LDA )
            IF( LJOBD )
     $         CALL DGEMM( 'N', 'N', N, M, P, ONE, DWORK, LDWN, D, LDD,
     $                     ONE, B, LDB )
         END IF
      END IF
C
      RETURN
C *** Last line of AB05SD ***
      END
