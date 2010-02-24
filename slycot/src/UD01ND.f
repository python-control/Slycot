      SUBROUTINE UD01ND( MP, NP, DP, L, NOUT, P, LDP1, LDP2, TEXT,
     $                   INFO )
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
C     To print the MP-by-NP coefficient matrices of a matrix polynomial
C                                                    dp-1           dp
C        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  .
C
C     The elements of the matrices are output to 7 significant figures.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     MP      (input) INTEGER
C             The number of rows of the matrix polynomial P(s).
C             MP >= 1.
C
C     NP      (input) INTEGER
C             The number of columns of the matrix polynomial P(s).
C             NP >= 1.
C
C     DP      (input) INTEGER
C             The degree of the matrix polynomial P(s).  DP >= 0.
C
C     L       (input) INTEGER
C             The number of elements of the coefficient matrices to be
C             printed per line.  1 <= L <= 5.
C
C     NOUT    (input) INTEGER
C             The output channel to which the results are sent.
C             NOUT >= 0.
C
C     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1)
C             The leading MP-by-NP-by-(DP+1) part of this array must
C             contain the coefficients of the matrix polynomial P(s).
C             Specifically, P(i,j,k) must contain the coefficient of
C             s**(k-1) of the polynomial which is the (i,j)-th element
C             of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and
C             k = 1,2,...,DP+1.
C
C     LDP1    INTEGER
C             The leading dimension of array P.  LDP1 >= MP.
C
C     LDP2    INTEGER
C             The second dimension of array P.  LDP2 >= NP.
C
C     TEXT    (input) CHARACTER*72
C             Title caption of the coefficient matrices to be printed.
C             TEXT is followed by the degree of the coefficient matrix,
C             within brackets. If TEXT = ' ', then the coefficient
C             matrices are separated by an empty line.
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
C     For i = 1, 2, ..., DP + 1 the routine first prints the contents of
C     TEXT followed by (i-1) as a title, followed by the elements of the
C     MP-by-NP coefficient matrix P(i) such that
C     (i)  if NP < L, then the leading MP-by-NP part is printed;
C     (ii) if NP = k*L + p (where k, p > 0), then k MP-by-L blocks of
C          consecutive columns of P(i) are printed one after another
C          followed by one MP-by-p block containing the last p columns
C          of P(i).
C     Row numbers are printed on the left of each row and a column
C     number on top of each column.
C
C     REFERENCES
C
C     None.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998.
C     Based on routine PRMAPO by A.J. Geurts, Eindhoven University of
C     Technology, Holland.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           DP, INFO, L, LDP1, LDP2, MP, NP, NOUT
      CHARACTER*(*)     TEXT
C     .. Array Arguments ..
      DOUBLE PRECISION  P(LDP1,LDP2,*)
C     .. Local Scalars ..
      INTEGER           I, J, J1, J2, JJ, K, LENTXT, LTEXT, N1
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MIN
C
C     .. Executable Statements ..
C
      INFO = 0
C
C     Check the input scalar arguments.
C
      IF( MP.LT.1 ) THEN
         INFO = -1
      ELSE IF( NP.LT.1 ) THEN
         INFO = -2
      ELSE IF( DP.LT.0 ) THEN
         INFO = -3
      ELSE IF( L.LT.1 .OR. L.GT.5 ) THEN
         INFO = -4
      ELSE IF( NOUT.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDP1.LT.MP ) THEN
         INFO = -7
      ELSE IF( LDP2.LT.NP ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'UD01ND', -INFO )
         RETURN
      END IF
C
      LENTXT = LEN( TEXT )
      LTEXT  = MIN( 72, LENTXT )
C     WHILE ( TEXT(LTEXT:LTEXT) =  ' ' ) DO
   10 IF ( TEXT(LTEXT:LTEXT).EQ.' ' ) THEN
         LTEXT = LTEXT - 1
         GO TO 10
      END IF
C     END WHILE 10
C
      DO 50 K = 1, DP + 1
         IF ( LTEXT.EQ.0 ) THEN
            WRITE ( NOUT, FMT = 99999 )
         ELSE
            WRITE ( NOUT, FMT = 99998 ) TEXT(1:LTEXT), K - 1, MP, NP
         END IF
         N1 = ( NP - 1 )/L
         J1 = 1
         J2 = L
C
         DO 30 J = 1, N1
            WRITE ( NOUT, FMT = 99997 ) ( JJ, JJ = J1, J2 )
C
            DO 20 I = 1, MP
               WRITE ( NOUT, FMT = 99996 ) I, ( P(I,JJ,K), JJ = J1, J2 )
   20       CONTINUE
C
            J1 = J1 + L
            J2 = J2 + L
   30    CONTINUE
C
         WRITE ( NOUT, FMT = 99997 ) ( J, J = J1, NP )
C
         DO 40 I = 1, MP
            WRITE ( NOUT, FMT = 99996 ) I, ( P(I,JJ,K), JJ = J1, NP )
   40    CONTINUE
C
   50 CONTINUE
C
      WRITE ( NOUT, FMT = 99999 )
C
      RETURN
C
99999 FORMAT (' ')
99998 FORMAT (/, 1X, A, '(', I2, ')', ' (', I2, 'X', I2, ')')
99997 FORMAT (5X, 5(6X, I2, 7X))
99996 FORMAT (1X, I2, 2X, 5D15.7)
C
C *** Last line of UD01ND ***
      END
