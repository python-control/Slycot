      SUBROUTINE TF01RD( NA, NB, NC, N, A, LDA, B, LDB, C, LDC, H, LDH,
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
C     To compute N Markov parameters M(1), M(2),..., M(N) from the
C     parameters (A,B,C) of a linear time-invariant system, where each
C     M(k) is an NC-by-NB matrix and k = 1,2,...,N.
C
C     All matrices are treated as dense, and hence TF01RD is not
C     intended for large sparse problems.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NA      (input) INTEGER
C             The order of the matrix A.  NA >= 0.
C
C     NB      (input) INTEGER
C             The number of system inputs.  NB >= 0.
C
C     NC      (input) INTEGER
C             The number of system outputs.  NC >= 0.
C
C     N       (input) INTEGER
C             The number of Markov parameters M(k) to be computed.
C             N >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
C             The leading NA-by-NA part of this array must contain the
C             state matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,NA).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,NB)
C             The leading NA-by-NB part of this array must contain the
C             input matrix B of the system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,NA).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,NA)
C             The leading NC-by-NA part of this array must contain the
C             output matrix C of the system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,NC).
C
C     H       (output) DOUBLE PRECISION array, dimension (LDH,N*NB)
C             The leading NC-by-N*NB part of this array contains the
C             multivariable parameters M(k), where each parameter M(k)
C             is an NC-by-NB matrix and k = 1,2,...,N. The Markov
C             parameters are stored such that H(i,(k-1)xNB+j) contains
C             the (i,j)-th element of M(k) for i = 1,2,...,NC and
C             j = 1,2,...,NB.
C
C     LDH     INTEGER
C             The leading dimension of array H.  LDH >= MAX(1,NC).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, 2*NA*NC).
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
C     For the linear time-invariant discrete-time system
C
C            x(k+1) = A x(k) + B u(k)
C             y(k)  = C x(k) + D u(k),
C
C     the transfer function matrix G(z) is given by
C                            -1
C              G(z) = C(zI-A)  B + D
C                             -1        -2     2   -3
C                   = D + CB z   + CAB z   + CA B z   + ...          (1)
C
C     Using Markov parameters, G(z) can also be written as
C                                 -1        -2        -3
C              G(z) = M(0) + M(1)z   + M(2)z   + M(3)z   + ...       (2)
C
C                                                               k-1
C     Equating (1) and (2), we find that M(0) = D and M(k) = C A    B
C     for k > 0, from which the Markov parameters M(1),M(2)...,M(N) are
C     computed.
C
C     REFERENCES
C
C     [1] Chen, C.T.
C         Introduction to Linear System Theory.
C         H.R.W. Series in Electrical Engineering, Electronics and
C         Systems, Holt, Rinehart and Winston Inc., London, 1970.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires approximately (NA + NB) x NA x NC x N
C     multiplications and additions.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TF01FD by S. Van Huffel, Katholieke
C     Univ. Leuven, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Markov parameters, multivariable system, time-invariant system,
C     transfer function, transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, LDC, LDH, LDWORK, N, NA, NB, NC
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), H(LDH,*)
C     .. Local Scalars ..
      INTEGER           I, JWORK, K, LDW
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DLACPY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( NA.LT.0 ) THEN
         INFO = -1
      ELSE IF( NB.LT.0 ) THEN
         INFO = -2
      ELSE IF( NC.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, NA ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, NA ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, NC ) ) THEN
         INFO = -10
      ELSE IF( LDH.LT.MAX( 1, NC ) ) THEN
         INFO = -12
      ELSE IF( LDWORK.LT.MAX( 1, 2*NA*NC ) ) THEN
         INFO = -14
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TF01RD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( NA, NB, NC, N ).EQ.0 )
     $   RETURN
C
      JWORK = 1 + NC*NA
      LDW   = MAX( 1, NC )
      I = 1
C
C     Copy C in the workspace beginning from the position JWORK.
C     This workspace will contain the product C*A**(K-1), K = 1,2,...,N.
C
      CALL DLACPY( 'Full', NC, NA, C, LDC, DWORK(JWORK), LDW )
C
C     Form M(1), M(2), ..., M(N).
C
      DO 10 K = 1, N
         CALL DLACPY( 'Full', NC, NA, DWORK(JWORK), LDW, DWORK, LDW )
C
C        Form (C * A**(K-1)) * B = M(K).
C
         CALL DGEMM( 'No transpose', 'No transpose', NC, NB, NA, ONE,
     $               DWORK, LDW, B, LDB, ZERO, H(1,I), LDH )
C
         IF ( K.NE.N ) THEN
C
C           Form C * A**K.
C
            CALL DGEMM( 'No transpose', 'No transpose', NC, NA, NA, ONE,
     $                  DWORK, LDW, A, LDA, ZERO, DWORK(JWORK), LDW )
C
            I = I + NB
         END IF
   10 CONTINUE
C
      RETURN
C *** Last line of TF01RD ***
      END
