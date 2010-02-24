      SUBROUTINE MB02CX( TYPET, P, Q, K, A, LDA, B, LDB, CS, LCS,
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
C     To bring the first blocks of a generator in proper form.
C     The columns / rows of the positive and negative generators
C     are contained in the arrays A and B, respectively.
C     Transformation information will be stored and can be applied
C     via SLICOT Library routine MB02CY.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPET   CHARACTER*1
C             Specifies the type of the generator, as follows:
C             = 'R':  A and B are the first blocks of the rows of the
C                     positive and negative generators;
C             = 'C':  A and B are the first blocks of the columns of the
C                     positive and negative generators.
C             Note:   in the sequel, the notation x / y means that
C                     x corresponds to TYPET = 'R' and y corresponds to
C                     TYPET = 'C'.
C
C     Input/Output Parameters
C
C     P       (input)  INTEGER
C             The number of rows / columns in A containing the positive
C             generators.  P >= 0.
C
C     Q       (input)  INTEGER
C             The number of rows / columns in B containing the negative
C             generators.  Q >= 0.
C
C     K       (input)  INTEGER
C             The number of columns / rows in A and B to be processed.
C             Normally, the size of the first block.  P >= K >= 0.
C
C     A       (input/output)  DOUBLE PRECISION array, dimension
C             (LDA, K) / (LDA, P)
C             On entry, the leading P-by-K upper / K-by-P lower
C             triangular part of this array must contain the rows /
C             columns of the positive part in the first block of the
C             generator.
C             On exit, the leading P-by-K upper / K-by-P lower
C             triangular part of this array contains the rows / columns
C             of the positive part in the first block of the proper
C             generator.
C             The lower / upper trapezoidal part is not referenced.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,P),    if TYPET = 'R';
C             LDA >= MAX(1,K),    if TYPET = 'C'.
C
C     B       (input/output)  DOUBLE PRECISION array, dimension
C             (LDB, K) / (LDB, Q)
C             On entry, the leading Q-by-K / K-by-Q part of this array
C             must contain the rows / columns of the negative part in
C             the first block of the generator.
C             On exit, the leading Q-by-K / K-by-Q part of this array
C             contains part of the necessary information for the
C             Householder transformations.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,Q),    if TYPET = 'R';
C             LDB >= MAX(1,K),    if TYPET = 'C'.
C
C     CS      (output)  DOUBLE PRECISION array, dimension (LCS)
C             On exit, the leading 2*K + MIN(K,Q) part of this array
C             contains necessary information for the SLICOT Library
C             routine MB02CY (modified hyperbolic rotation parameters
C             and scalar factors of the Householder transformations).
C
C     LCS     INTEGER
C             The length of the array CS.  LCS >= 2*K + MIN(K,Q).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if  INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -12,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,K).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  succesful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the reduction algorithm failed. The matrix
C                   associated with the generator is not (numerically)
C                   positive definite.
C
C     METHOD
C
C     If  TYPET = 'R',  a QR decomposition of B is first computed.
C     Then, the elements below the first row of each column i of B
C     are annihilated by a Householder transformation modifying the
C     first element in that column. This first element, in turn, is
C     then annihilated by a modified hyperbolic rotation, acting also
C     on the i-th row of A.
C
C     If  TYPET = 'C',  an LQ decomposition of B is first computed.
C     Then, the elements on the right of the first column of each row i
C     of B are annihilated by a Householder transformation modifying the
C     first element in that row. This first element, in turn, is
C     then annihilated by a modified hyperbolic rotation, acting also
C     on the i-th column of A.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2000,
C     February 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, Householder transformation, matrix
C     operations, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TYPET
      INTEGER           INFO, K, LDA, LDB, LCS, LDWORK, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA, *), B(LDB, *), CS(*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           ISROW
      INTEGER           I, IERR
      DOUBLE PRECISION  ALPHA, BETA, C, MAXWRK, S, TAU
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGELQF, DGEQRF, DLARF, DLARFG, DSCAL,
     $                  MA02FD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      ISROW = LSAME( TYPET, 'R' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( ISROW .OR. LSAME( TYPET, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF ( P.LT.0 ) THEN
         INFO = -2
      ELSE IF ( Q.LT.0 ) THEN
         INFO = -3
      ELSE IF ( K.LT.0 .OR. K.GT.P ) THEN
         INFO = -4
      ELSE IF ( LDA.LT.1 .OR. ( ISROW .AND. LDA.LT.P ) .OR.
     $                   ( .NOT.ISROW .AND. LDA.LT.K ) ) THEN
         INFO = -6
      ELSE IF ( LDB.LT.1 .OR. ( ISROW .AND. LDB.LT.Q ) .OR.
     $                   ( .NOT.ISROW .AND. LDB.LT.K ) ) THEN
         INFO = -8
      ELSE IF ( LCS.LT.2*K + MIN( K, Q ) ) THEN
         INFO = -10
      ELSE IF ( LDWORK.LT.MAX( 1, K ) ) THEN
         DWORK(1) = MAX( 1, K )
         INFO = -12
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02CX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( Q, K ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( ISROW ) THEN
C
C        The generator is row wise stored.
C
C        Step 0: Do QR decomposition of B.
C
         CALL DGEQRF ( Q, K, B, LDB, CS(2*K+1), DWORK(1), LDWORK, IERR )
         MAXWRK = DWORK(1)
C
         DO 10  I = 1, K
C
C           Step 1: annihilate the i-th column of B.
C
            IF ( Q.GT.1 ) THEN
               CALL DLARFG( MIN( I, Q ), B(1,I), B(2,I), 1, TAU )
               ALPHA  = B(1,I)
               B(1,I) = ONE
               IF ( K.GT.I )
     $            CALL DLARF( 'Left', MIN( I, Q ), K-I, B(1,I), 1, TAU,
     $                        B(1,I+1), LDB, DWORK )
               B(1,I) = ALPHA
            ELSE
               ALPHA = B(1,I)
               TAU   = ZERO
            END IF
C
C           Step 2: annihilate the top entry of the column.
C
            BETA = A(I,I)
            CALL MA02FD( BETA, ALPHA, C, S, IERR )
            IF ( IERR.NE.0 ) THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            CS(I*2-1) = C
            CS(I*2)   = S
            CALL DSCAL( K-I+1, ONE/C, A(I,I), LDA )
            CALL DAXPY( K-I+1,  -S/C, B(1,I), LDB, A(I,I), LDA )
            CALL DSCAL( K-I+1,     C, B(1,I), LDB )
            CALL DAXPY( K-I+1,    -S, A(I,I), LDA, B(1,I), LDB )
            B(1,I) = TAU
   10    CONTINUE
C
      ELSE
C
C        The generator is column wise stored.
C
C        Step 0: Do LQ decomposition of B.
C
         CALL DGELQF ( K, Q, B, LDB, CS(2*K+1), DWORK(1), LDWORK, IERR )
         MAXWRK = DWORK(1)
C
         DO 20  I = 1, K
C
C           Step 1: annihilate the i-th row of B.
C
            IF ( Q.GT.1 ) THEN
               CALL DLARFG( MIN( I, Q ), B(I,1), B(I,2), LDB, TAU )
               ALPHA  = B(I,1)
               B(I,1) = ONE
               IF ( K.GT.I )
     $            CALL DLARF( 'Right', K-I, MIN( I, Q ), B(I,1), LDB,
     $                        TAU, B(I+1,1), LDB, DWORK )
               B(I,1) = ALPHA
            ELSE
               ALPHA = B(I,1)
               TAU   = ZERO
            END IF
C
C           Step 2: annihilate the left entry of the row.
C
            BETA = A(I,I)
            CALL MA02FD( BETA, ALPHA, C, S, IERR )
            IF ( IERR.NE.0 ) THEN
C
C              Error return:  The matrix is not positive definite.
C
               INFO = 1
               RETURN
            END IF
C
            CS(I*2-1) = C
            CS(I*2)   = S
            CALL DSCAL( K-I+1, ONE/C, A(I,I), 1 )
            CALL DAXPY( K-I+1,  -S/C, B(I,1), 1, A(I,I), 1 )
            CALL DSCAL( K-I+1,     C, B(I,1), 1 )
            CALL DAXPY( K-I+1,    -S, A(I,I), 1, B(I,1), 1 )
            B(I,1) = TAU
   20    CONTINUE
C
      END IF
C
      DWORK(1) = MAXWRK
C
      RETURN
C
C *** Last line of MB02CX ***
      END
