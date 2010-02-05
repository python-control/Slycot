      SUBROUTINE MB02CY( TYPET, STRUCG, P, Q, N, K, A, LDA, B, LDB, H,
     $                   LDH, CS, LCS, DWORK, LDWORK, INFO )
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
C     To apply the transformations created by the SLICOT Library
C     routine MB02CX on other columns / rows of the generator,
C     contained in the arrays A and B of positive and negative
C     generators, respectively.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPET   CHARACTER*1
C             Specifies the type of the generator, as follows:
C             = 'R':  A and B are additional columns of the generator;
C             = 'C':  A and B are additional rows of the generator.
C             Note:   in the sequel, the notation x / y means that
C                     x corresponds to TYPET = 'R' and y corresponds to
C                     TYPET = 'C'.
C
C     STRUCG  CHARACTER*1
C             Information about the structure of the two generators,
C             as follows:
C             = 'T':  the trailing block of the positive generator
C                     is lower / upper triangular, and the trailing
C                     block of the negative generator is zero;
C             = 'N':  no special structure to mention.
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
C     N       (input)  INTEGER
C             The number of columns / rows in A and B to be processed.
C             N >= 0.
C
C     K       (input)  INTEGER
C             The number of columns / rows in H.  P >= K >= 0.
C
C     A       (input/output)  DOUBLE PRECISION array, dimension
C             (LDA, N) / (LDA, P)
C             On entry, the leading P-by-N / N-by-P part of this array
C             must contain the positive part of the generator.
C             On exit, the leading P-by-N / N-by-P part of this array
C             contains the transformed positive part of the generator.
C
C     LDA     INTEGER
C             The leading dimension of the array A.
C             LDA >= MAX(1,P),    if TYPET = 'R';
C             LDA >= MAX(1,N),    if TYPET = 'C'.
C
C     B       (input/output)  DOUBLE PRECISION array, dimension
C             (LDB, N) / (LDB, Q)
C             On entry, the leading Q-by-N / N-by-Q part of this array
C             must contain the negative part of the generator.
C             On exit, the leading Q-by-N / N-by-Q part of this array
C             contains the transformed negative part of the generator.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             LDB >= MAX(1,Q),    if TYPET = 'R';
C             LDB >= MAX(1,N),    if TYPET = 'C'.
C
C     H       (input)  DOUBLE PRECISION array, dimension
C             (LDH, K) / (LDH, Q)
C             The leading Q-by-K / K-by-Q part of this array must
C             contain part of the necessary information for the
C             Householder transformations computed by SLICOT Library
C             routine MB02CX.
C
C     LDH     INTEGER
C             The leading dimension of the array H.
C             LDH >= MAX(1,Q),    if TYPET = 'R';
C             LDH >= MAX(1,K),    if TYPET = 'C'.
C
C     CS      (input)  DOUBLE PRECISION array, dimension (LCS)
C             The leading 2*K + MIN(K,Q) part of this array must
C             contain the necessary information for modified hyperbolic
C             rotations and the scalar factors of the Householder
C             transformations computed by SLICOT Library routine MB02CX.
C
C     LCS     INTEGER
C             The length of the array CS.  LCS >= 2*K + MIN(K,Q).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if  INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -16,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= MAX(1,N).
C             For optimum performance LDWORK should be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  succesful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The Householder transformations and modified hyperbolic rotations
C     computed by SLICOT Library routine MB02CX are applied to the
C     corresponding parts of the matrices A and B.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, July 2000,
C     February 2004, March 2007.
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
      INTEGER           INFO, K, LDA, LDB, LCS, LDH, LDWORK, N, P, Q
      CHARACTER         STRUCG, TYPET
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA, *), B(LDB, *), CS(*), DWORK(*), H(LDH,*)
C     .. Local Scalars ..
      LOGICAL           ISLWR, ISROW
      INTEGER           I, IERR, CI, MAXWRK
      DOUBLE PRECISION  C, S, TAU
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLARF, DLASET, DORMLQ, DORMQR, DSCAL,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      ISROW = LSAME( TYPET,  'R' )
      ISLWR = LSAME( STRUCG, 'T' )
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( ISROW .OR. LSAME( TYPET, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( ISLWR .OR. LSAME( STRUCG, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF ( P.LT.0 ) THEN
         INFO = -3
      ELSE IF ( Q.LT.0 ) THEN
         INFO = -4
      ELSE IF ( N.LT.0 ) THEN
         INFO = -5
      ELSE IF ( K.LT.0 .OR. K.GT.P ) THEN
         INFO = -6
      ELSE IF ( LDA.LT.1 .OR. ( ISROW .AND. LDA.LT.P ) .OR.
     $                   ( .NOT.ISROW .AND. LDA.LT.N ) ) THEN
         INFO = -8
      ELSE IF ( LDB.LT.1 .OR. ( ISROW .AND. LDB.LT.Q ) .OR.
     $                   ( .NOT.ISROW .AND. LDB.LT.N ) ) THEN
         INFO = -10
      ELSE IF ( LDH.LT.1 .OR. ( ISROW .AND. LDH.LT.Q ) .OR.
     $                   ( .NOT.ISROW .AND. LDH.LT.K ) ) THEN
         INFO = -12
      ELSE IF ( LCS.LT.2*K + MIN( K, Q ) ) THEN
         INFO = -14
      ELSE IF ( LDWORK.LT.MAX( 1, N ) ) THEN
         DWORK(1) = MAX( 1, N )
         INFO = -16
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02CY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MIN( N, K, Q ).EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
C     Applying the transformations.
C
      IF ( ISROW ) THEN
C
C        The generator is row wise stored.
C
         IF ( ISLWR ) THEN
C
            DO 10  I = 1, K
C
C              Apply Householder transformation avoiding touching of
C              zero blocks.
C
               CI  = N - K + I - 1
               TAU = H(1,I)
               H(1,I) = ONE
               CALL DLARF( 'Left', MIN( I, Q ), CI, H(1,I), 1, TAU, B,
     $                     LDB, DWORK )
               H(1,I) = TAU
C
C              Now apply the hyperbolic rotation under the assumption
C              that A(I, N-K+I+1:N) and B(1, N-K+I:N) are zero.
C
               C = CS(I*2-1)
               S = CS(I*2)
C
               CALL DSCAL( CI, ONE/C, A(I,1), LDA )
               CALL DAXPY( CI,  -S/C, B(1,1), LDB, A(I,1), LDA )
               CALL DSCAL( CI,     C, B(1,1), LDB )
               CALL DAXPY( CI,    -S, A(I,1), LDA, B(1,1), LDB )
C
               B(1,N-K+I) =  -S/C * A(I,N-K+I)
               A(I,N-K+I) = ONE/C * A(I,N-K+I)
C
C              All below B(1,N-K+I) should be zero.
C
               IF( Q.GT.1 )
     $            CALL DLASET( 'All', Q-1, 1, ZERO, ZERO, B(2,N-K+I),
     $                         1 )
   10       CONTINUE
C
         ELSE
C
C           Apply the QR reduction on B.
C
            CALL DORMQR( 'Left', 'Transpose', Q, N, MIN( K, Q ), H,
     $                   LDH, CS(2*K+1), B, LDB, DWORK, LDWORK, IERR )
            MAXWRK = DWORK(1)
C
            DO 20  I = 1, K
C
C              Apply Householder transformation.
C
               TAU = H(1,I)
               H(1,I) = ONE
               CALL DLARF( 'Left', MIN( I, Q ), N, H(1,I), 1, TAU, B,
     $                     LDB, DWORK )
               H(1,I) = TAU
C
C              Apply Hyperbolic Rotation.
C
               C = CS(I*2-1)
               S = CS(I*2)
C
               CALL DSCAL( N, ONE/C, A(I,1), LDA )
               CALL DAXPY( N,  -S/C, B(1,1), LDB, A(I,1), LDA )
               CALL DSCAL( N,     C, B(1,1), LDB )
               CALL DAXPY( N,    -S, A(I,1), LDA, B(1,1), LDB )
   20       CONTINUE
C
         END IF
C
      ELSE
C
C        The generator is column wise stored.
C
         IF ( ISLWR ) THEN
C
            DO 30  I = 1, K
C
C              Apply Householder transformation avoiding touching zeros.
C
               CI  = N - K + I - 1
               TAU = H(I,1)
               H(I,1) = ONE
               CALL DLARF( 'Right', CI, MIN( I, Q ), H(I,1), LDH, TAU,
     $                     B, LDB, DWORK )
               H(I,1) = TAU
C
C              Apply Hyperbolic Rotation.
C
               C = CS(I*2-1)
               S = CS(I*2)
C
               CALL DSCAL( CI, ONE/C, A(1,I), 1 )
               CALL DAXPY( CI,  -S/C, B(1,1), 1, A(1,I), 1 )
               CALL DSCAL( CI,     C, B(1,1), 1 )
               CALL DAXPY( CI,    -S, A(1,I), 1, B(1,1), 1 )
C
               B(N-K+I,1) =  -S/C * A(N-K+I,I)
               A(N-K+I,I) = ONE/C * A(N-K+I,I)
C
C              All elements right behind B(N-K+I,1) should be zero.
C
               IF( Q.GT.1 )
     $            CALL DLASET( 'All', 1, Q-1, ZERO, ZERO, B(N-K+I,2),
     $                         LDB )
   30       CONTINUE
C
         ELSE
C
C           Apply the LQ reduction on B.
C
            CALL DORMLQ( 'Right', 'Transpose', N, Q, MIN( K, Q ), H,
     $                   LDH, CS(2*K+1), B, LDB, DWORK, LDWORK, IERR )
            MAXWRK = DWORK(1)
C
            DO 40  I = 1, K
C
C              Apply Householder transformation.
C
               TAU = H(I,1)
               H(I,1) = ONE
               CALL DLARF( 'Right', N, MIN( I, Q ), H(I,1), LDH, TAU, B,
     $                     LDB, DWORK )
               H(I,1) = TAU
C
C              Apply Hyperbolic Rotation.
C
               C = CS(I*2-1)
               S = CS(I*2)
C
               CALL DSCAL( N, ONE/C, A(1,I), 1 )
               CALL DAXPY( N,  -S/C, B(1,1), 1, A(1,I), 1 )
               CALL DSCAL( N,     C, B(1,1), 1 )
               CALL DAXPY( N,    -S, A(1,I), 1, B(1,1), 1 )
   40       CONTINUE
C
         END IF
C
      END IF
C
      DWORK(1) = MAX( MAXWRK, N )
C
      RETURN
C
C *** Last line of MB02CY ***
      END
