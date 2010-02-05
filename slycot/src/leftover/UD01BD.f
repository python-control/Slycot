      SUBROUTINE UD01BD( MP, NP, DP, NIN, P, LDP1, LDP2, INFO )
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
C     To read the coefficients of a matrix polynomial
C                                                    dp-1           dp
C        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  .
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
C     NIN     (input) INTEGER
C             The input channel from which the elements of P(s) are
C             read.  NIN >= 0.
C
C     P       (output) DOUBLE PRECISION array, dimension
C             (LDP1,LDP2,DP+1)
C             The leading MP-by-NP-by-(DP+1) part of this array contains
C             the coefficients of the matrix polynomial P(s).
C             Specifically, P(i,j,k) contains the coefficient of
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
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The coefficients P(i), i = 0, ..., DP, which are MP-by-NP
C     matrices, are read from the input file NIN row by row. Each P(i)
C     must be preceded by a text line. This text line can be used to
C     indicate the coefficient matrices.
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
C     Based on routine RDMAPO by A.J. Geurts, Eindhoven University of
C     Technology, Holland.
C
C     REVISIONS
C
C     -
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DP, INFO, LDP1, LDP2, MP, NP, NIN
C     .. Array Arguments ..
      DOUBLE PRECISION  P(LDP1,LDP2,*)
C     .. Local Scalars ..
      INTEGER           I, J, K
C     .. External Subroutines ..
      EXTERNAL          XERBLA
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
      ELSE IF( NIN.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDP1.LT.MP ) THEN
         INFO = -6
      ELSE IF( LDP2.LT.NP ) THEN
         INFO = -7
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'UD01BD', -INFO )
         RETURN
      END IF
C
C     Skip the text line preceding P(i) and read P(i), i = 0, ..., DP,
C     row after row.
C
      DO 20 K = 1, DP + 1
         READ ( NIN, FMT = '()' )
C
         DO 10 I = 1, MP
            READ ( NIN, FMT = * ) ( P(I,J,K), J = 1, NP )
   10    CONTINUE
C
   20 CONTINUE
C
      RETURN
C *** Last line of UD01BD ***
      END
