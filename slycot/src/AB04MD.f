      SUBROUTINE AB04MD( TYPE, N, M, P, ALPHA, BETA, A, LDA, B, LDB, C,
     $                   LDC, D, LDD, IWORK, DWORK, LDWORK, INFO )
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
C     To perform a transformation on the parameters (A,B,C,D) of a
C     system, which is equivalent to a bilinear transformation of the
C     corresponding transfer function matrix.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPE    CHARACTER*1
C             Indicates the type of the original system and the
C             transformation to be performed as follows:
C             = 'D':  discrete-time   -> continuous-time;
C             = 'C':  continuous-time -> discrete-time.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the state matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     ALPHA,  (input) DOUBLE PRECISION
C     BETA    Parameters specifying the bilinear transformation.
C             Recommended values for stable systems: ALPHA = 1,
C             BETA = 1.  ALPHA <> 0, BETA <> 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state matrix A of the original system.
C             On exit, the leading N-by-N part of this array contains
C                              _
C             the state matrix A of the transformed system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input matrix B of the original system.
C             On exit, the leading N-by-M part of this array contains
C                              _
C             the input matrix B of the transformed system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the output matrix C of the original system.
C             On exit, the leading P-by-N part of this array contains
C                               _
C             the output matrix C of the transformed system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the input/output matrix D for the original system.
C             On exit, the leading P-by-M part of this array contains
C                                     _
C             the input/output matrix D of the transformed system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
C             For optimum performance LDWORK >= MAX(1,N*NB), where NB
C             is the optimal blocksize.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the matrix (ALPHA*I + A) is exactly singular;
C             = 2:  if the matrix  (BETA*I - A) is exactly singular.
C
C     METHOD
C
C     The parameters of the discrete-time system are transformed into
C     the parameters of the continuous-time system (TYPE = 'D'), or
C     vice-versa (TYPE = 'C') by the transformation:
C
C     1.  Discrete -> continuous
C         _                     -1
C         A = beta*(alpha*I + A)  * (A - alpha*I)
C         _                                     -1
C         B = sqrt(2*alpha*beta) * (alpha*I + A)  * B
C         _                                         -1
C         C = sqrt(2*alpha*beta) * C * (alpha*I + A)
C         _                        -1
C         D = D - C * (alpha*I + A)  * B
C
C     which is equivalent to the bilinear transformation
C
C                       z - alpha
C         z -> s = beta ---------  .
C                       z + alpha
C
C     of one transfer matrix onto the other.
C
C     2.  Continuous -> discrete
C         _                     -1
C         A = alpha*(beta*I - A)  * (beta*I + A)
C         _                                    -1
C         B = sqrt(2*alpha*beta) * (beta*I - A)  * B
C         _                                        -1
C         C = sqrt(2*alpha*beta) * C * (beta*I - A)
C         _                       -1
C         D = D + C * (beta*I - A)  * B
C
C     which is equivalent to the bilinear transformation
C
C                      beta + s
C       s -> z = alpha -------- .
C                      beta - s
C
C     of one transfer matrix onto the other.
C
C     REFERENCES
C
C     [1] Al-Saggaf, U.M. and Franklin, G.F.
C         Model reduction via balanced realizations: a extension and
C         frequency weighting techniques.
C         IEEE Trans. Autom. Contr., AC-33, pp. 687-692, 1988.
C
C     NUMERICAL ASPECTS
C                                                      3
C     The time taken is approximately proportional to N .
C     The accuracy depends mainly on the condition number of the matrix
C     to be inverted.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and
C                  A. Varga, German Aerospace Research Establishment,
C                  Oberpfaffenhofen, Germany, Nov. 1996.
C     Supersedes Release 2.0 routine AB04AD by W. van der Linden, and
C     A.J. Geurts, Technische Hogeschool Eindhoven, Holland.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Bilinear transformation, continuous-time system, discrete-time
C     system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TYPE
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDWORK, M, N, P
      DOUBLE PRECISION  ALPHA, BETA
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LTYPE
      INTEGER           I, IP
      DOUBLE PRECISION  AB2, PALPHA, PBETA, SQRAB2
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGETRF, DGETRS, DGETRI, DLASCL, DSCAL,
     $                  DSWAP, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      LTYPE = LSAME( TYPE, 'D' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LTYPE .AND. .NOT.LSAME( TYPE, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( ALPHA.EQ.ZERO ) THEN
         INFO = -5
      ELSE IF( BETA.EQ.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( LDWORK.LT.MAX( 1, N ) ) THEN
         INFO = -17
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB04MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( N, M, P ).EQ.0 )
     $   RETURN
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF (LTYPE) THEN
C
C        Discrete-time to continuous-time with (ALPHA, BETA).
C
         PALPHA = ALPHA
         PBETA = BETA
      ELSE
C
C        Continuous-time to discrete-time with (ALPHA, BETA) is
C        equivalent with discrete-time to continuous-time with
C        (-BETA, -ALPHA), if B and C change the sign.
C
         PALPHA = -BETA
         PBETA = -ALPHA
      END IF
C
      AB2 = PALPHA*PBETA*TWO
      SQRAB2 = SIGN( SQRT( ABS( AB2 ) ), PALPHA )
C                          -1
C     Compute (alpha*I + A)  .
C
      DO 10 I = 1, N
         A(I,I)  =  A(I,I) + PALPHA
   10 CONTINUE
C
      CALL DGETRF( N, N, A, LDA, IWORK, INFO )
C
      IF (INFO.NE.0) THEN
C
C        Error return.
C
         IF (LTYPE) THEN
            INFO = 1
         ELSE
            INFO = 2
         END IF
         RETURN
      END IF
C                         -1
C     Compute  (alpha*I+A)  *B.
C
      CALL DGETRS( 'No transpose', N, M, A, LDA, IWORK, B, LDB, INFO )
C                               -1
C     Compute  D - C*(alpha*I+A)  *B.
C
      CALL DGEMM( 'No transpose', 'No transpose', P, M, N, -ONE, C,
     $            LDC, B, LDB, ONE, D, LDD )
C
C     Scale B by  sqrt(2*alpha*beta).
C
      CALL DLASCL( 'General', 0, 0, ONE, SQRAB2, N, M, B, LDB, INFO )
C                                                -1
C     Compute  sqrt(2*alpha*beta)*C*(alpha*I + A)  .
C
      CALL DTRSM( 'Right', 'Upper', 'No transpose', 'Non-unit', P, N,
     $            SQRAB2, A, LDA, C, LDC )
C
      CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', P, N, ONE,
     $            A, LDA, C, LDC )
C
C     Apply column interchanges to the solution matrix.
C
      DO 20 I = N-1, 1, -1
         IP = IWORK(I)
         IF ( IP.NE.I )
     $      CALL DSWAP( P, C(1,I), 1, C(1,IP), 1 )
  20  CONTINUE
C                               -1
C     Compute beta*(alpha*I + A)  *(A - alpha*I) as
C                                        -1
C     beta*I - 2*alpha*beta*(alpha*I + A)  .
C
C     Workspace: need N;  prefer N*NB.
C
      CALL DGETRI( N, A, LDA, IWORK, DWORK, LDWORK, INFO )
C
      DO 30 I = 1, N
         CALL DSCAL(N, -AB2, A(1,I), 1)
         A(I,I) = A(I,I) + PBETA
   30 CONTINUE
C
      RETURN
C *** Last line of AB04MD ***
      END
