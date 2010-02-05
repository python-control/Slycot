      SUBROUTINE MB01UX( SIDE, UPLO, TRANS, M, N, ALPHA, T, LDT, A, LDA,
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
C     To compute one of the matrix products
C
C       A : = alpha*op( T ) * A, or A : = alpha*A * op( T ),
C
C     where alpha is a scalar, A is an m-by-n matrix, T is a quasi-
C     triangular matrix, and op( T ) is one of
C
C        op( T ) = T   or   op( T ) = T',  the transpose of T.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SIDE    CHARACTER*1
C             Specifies whether the upper quasi-triangular matrix H
C             appears on the left or right in the matrix product as
C             follows:
C             = 'L':  A := alpha*op( T ) * A;
C             = 'R':  A := alpha*A * op( T ).
C
C     UPLO    CHARACTER*1.
C             Specifies whether the matrix T is an upper or lower
C             quasi-triangular matrix as follows:
C             = 'U':  T is an upper quasi-triangular matrix;
C             = 'L':  T is a lower quasi-triangular matrix.
C
C     TRANS   CHARACTER*1
C             Specifies the form of op( T ) to be used in the matrix
C             multiplication as follows:
C             = 'N':  op( T ) = T;
C             = 'T':  op( T ) = T';
C             = 'C':  op( T ) = T'.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix A.  M >= 0.
C
C     N       (input) INTEGER
C             The number of columns of the matrix A.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The scalar alpha. When alpha is zero then T is not
C             referenced and A need not be set before entry.
C
C     T       (input) DOUBLE PRECISION array, dimension (LDT,k)
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C             On entry with UPLO = 'U', the leading k-by-k upper
C             Hessenberg part of this array must contain the upper
C             quasi-triangular matrix T. The elements below the
C             subdiagonal are not referenced.
C             On entry with UPLO = 'L', the leading k-by-k lower
C             Hessenberg part of this array must contain the lower
C             quasi-triangular matrix T. The elements above the
C             supdiagonal are not referenced.
C
C     LDT     INTEGER
C             The leading dimension of the array T.  LDT >= max(1,k),
C             where k is M when SIDE = 'L' and is N when SIDE = 'R'.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading M-by-N part of this array must
C             contain the matrix A.
C             On exit, the leading M-by-N part of this array contains
C             the computed product.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,M).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0 and ALPHA<>0,  DWORK(1)  returns the
C             optimal value of LDWORK.
C             On exit, if  INFO = -12,  DWORK(1)  returns the minimum
C             value of LDWORK.
C             This array is not referenced when alpha = 0.
C
C     LDWORK  The length of the array DWORK.
C             LDWORK >= 1,       if alpha =  0 or MIN(M,N) = 0;
C             LDWORK >= 2*(M-1), if SIDE  = 'L';
C             LDWORK >= 2*(N-1), if SIDE  = 'R'.
C             For maximal efficiency LDWORK should be at least
C             NOFF*N + M - 1,    if SIDE  = 'L';
C             NOFF*M + N - 1,    if SIDE  = 'R';
C             where NOFF is the number of nonzero elements on the
C             subdiagonal (if UPLO = 'U') or supdiagonal (if UPLO = 'L')
C             of T.
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
C     The technique used in this routine is similiar to the technique
C     used in the SLICOT [1] subroutine MB01UW developed by Vasile Sima.
C     The required matrix product is computed in two steps. In the first
C     step, the triangle of T specified by UPLO is used; in the second
C     step, the contribution of the sub-/supdiagonal is added. If the
C     workspace can accommodate parts of A, a fast BLAS 3 DTRMM
C     operation is used in the first step.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., Sima, V., Van Huffel, S., and
C         Varga, A.
C         SLICOT - A subroutine library in systems and control theory.
C         In: Applied and computational control, signals, and circuits,
C         Vol. 1, pp. 499-539, Birkhauser, Boston, 1999.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, May 2008 (SLICOT version of the HAPACK routine DTRQML).
C
C     KEYWORDS
C
C     Elementary matrix operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS, UPLO
      INTEGER           INFO, LDA, LDT, LDWORK, M, N
      DOUBLE PRECISION  ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), DWORK(*), T(LDT,*)
C     .. Local Scalars ..
      LOGICAL           LSIDE, LTRAN, LUP
      CHARACTER         ATRAN
      INTEGER           I, IERR, J, K, NOFF, PDW, PSAV, WRKMIN, WRKOPT,
     $                  XDIF
      DOUBLE PRECISION  TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DLASCL, DLASET, DTRMM, DTRMV,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C
C     .. Executable Statements ..
C
C     Decode and test the input scalar arguments.
C
      INFO  = 0
      LSIDE = LSAME( SIDE,  'L' )
      LUP   = LSAME( UPLO,  'U' )
      LTRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      IF ( LSIDE ) THEN
         K = M
      ELSE
         K = N
      END IF
      WRKMIN = 2*( K - 1 )
C
      IF ( ( .NOT.LSIDE ).AND.( .NOT.LSAME( SIDE, 'R' ) ) ) THEN
         INFO = -1
      ELSE IF ( ( .NOT.LUP ).AND.( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSE IF ( ( .NOT.LTRAN ).AND.( .NOT.LSAME( TRANS, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF ( M.LT.0 ) THEN
         INFO = -4
      ELSE IF ( N.LT.0 ) THEN
         INFO = -5
      ELSE IF ( LDT.LT.MAX( 1, K ) ) THEN
         INFO = -8
      ELSE IF ( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF ( LDWORK.LT.0 .OR.
     $          ( ALPHA.NE.ZERO .AND. MIN( M, N ).GT.0 .AND.
     $            LDWORK.LT.WRKMIN ) ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB01UX', -INFO )
         RETURN
      END IF
C
C     Quick return, if possible.
C
      IF ( MIN( M, N ).EQ.0 )
     $   RETURN
C
      IF ( ALPHA.EQ.ZERO ) THEN
C
C        Set A to zero and return.
C
         CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
         RETURN
      END IF
C
C     Save and count off-diagonal entries of T.
C
      IF ( LUP ) THEN
         CALL DCOPY( K-1, T(2,1), LDT+1, DWORK, 1 )
      ELSE
         CALL DCOPY( K-1, T(1,2), LDT+1, DWORK, 1 )
      END IF
      NOFF = 0
      DO 5 I = 1, K-1
         IF ( DWORK(I).NE.ZERO )
     $      NOFF = NOFF + 1
    5 CONTINUE
C
C     Compute optimal workspace.
C
      IF ( LSIDE ) THEN
         WRKOPT = NOFF*N + M - 1
      ELSE
         WRKOPT = NOFF*M + N - 1
      END IF
      PSAV = K
      IF ( .NOT.LTRAN ) THEN
         XDIF = 0
      ELSE
         XDIF = 1
      END IF
      IF ( .NOT.LUP )
     $   XDIF = 1 - XDIF
      IF ( .NOT.LSIDE )
     $   XDIF = 1 - XDIF
C
      IF ( LDWORK.GE.WRKOPT ) THEN
C
C        Enough workspace for a fast BLAS 3 calculation.
C        Save relevant parts of A in the workspace and compute one of
C        the matrix products
C          A : = alpha*op( triu( T ) ) * A, or
C          A : = alpha*A * op( triu( T ) ),
C        involving the upper/lower triangle of T.
C
         PDW = PSAV
         IF ( LSIDE ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, M-1
                  IF ( DWORK(I).NE.ZERO ) THEN
                     DWORK(PDW) = A(I+XDIF,J)
                     PDW = PDW + 1
                  END IF
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 30 J = 1, N-1
               IF ( DWORK(J).NE.ZERO ) THEN
                  CALL DCOPY( M, A(1,J+XDIF), 1, DWORK(PDW), 1 )
                  PDW = PDW + M
               END IF
   30       CONTINUE
         END IF
         CALL DTRMM( SIDE, UPLO, TRANS, 'Non-unit', M, N, ALPHA, T,
     $               LDT, A, LDA )
C
C        Add the contribution of the offdiagonal of T.
C
         PDW = PSAV
         XDIF = 1 - XDIF
         IF( LSIDE ) THEN
            DO 50 J = 1, N
               DO 40 I = 1, M-1
                  TEMP = DWORK(I)
                  IF ( TEMP.NE.ZERO ) THEN
                     A(I+XDIF,J) = A(I+XDIF,J) + ALPHA * TEMP *
     $                                           DWORK(PDW)
                     PDW = PDW + 1
                  END IF
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 60 J = 1, N-1
               TEMP = DWORK(J)*ALPHA
               IF ( TEMP.NE.ZERO ) THEN
                  CALL DAXPY( M, TEMP, DWORK(PDW), 1, A(1,J+XDIF), 1 )
                  PDW = PDW + M
               END IF
   60       CONTINUE
         END IF
      ELSE
C
C        Use a BLAS 2 calculation.
C
         IF ( LSIDE ) THEN
            DO 80 J = 1, N
C
C              Compute the contribution of the offdiagonal of T to
C              the j-th column of the product.
C
               DO 70 I = 1, M - 1
                  DWORK(PSAV+I-1) = DWORK(I)*A(I+XDIF,J)
   70          CONTINUE
C
C              Multiply the triangle of T by the j-th column of A,
C              and add to the above result.
C
               CALL DTRMV( UPLO, TRANS, 'Non-unit', M, T, LDT, A(1,J),
     $                     1 )
               CALL DAXPY( M-1, ONE, DWORK(PSAV), 1, A(2-XDIF,J), 1 )
   80       CONTINUE
         ELSE
            IF ( LTRAN ) THEN
               ATRAN = 'N'
            ELSE
               ATRAN = 'T'
            END IF
            DO 100 I = 1, M
C
C              Compute the contribution of the offdiagonal of T to
C              the i-th row of the product.
C
               DO 90 J = 1, N - 1
                  DWORK(PSAV+J-1) = A(I,J+XDIF)*DWORK(J)
   90          CONTINUE
C
C              Multiply the i-th row of A by the triangle of T,
C              and add to the above result.
C
               CALL DTRMV( UPLO, ATRAN, 'Non-unit', N, T, LDT, A(I,1),
     $                     LDA )
               CALL DAXPY( N-1, ONE, DWORK(PSAV), 1, A(I,2-XDIF), LDA )
  100       CONTINUE
         END IF
C
C        Scale the result by alpha.
C
         IF ( ALPHA.NE.ONE )
     $      CALL DLASCL( 'General', 0, 0, ONE, ALPHA, M, N, A, LDA,
     $                   IERR )
      END IF
      DWORK(1) = DBLE( MAX( WRKMIN, WRKOPT ) )
      RETURN
C *** Last line of MB01UX ***
      END
