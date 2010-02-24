      SUBROUTINE UD01CD( MP, NP, DP, NIN, P, LDP1, LDP2, INFO )
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
C     To read the elements of a sparse matrix polynomial
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
C             The not assigned elements are set to zero.
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
C                   value;
C             = 1 : if a row index i is read with i < 1 or i > MP or
C                   a column index j is read with j < 1 or j > NP or
C                   a coefficient degree d is read with d < 0 or
C                   d > DP + 1. This is a warning.
C
C     METHOD
C
C     First, the elements P(i,j,k) with 1 <= i <= MP, 1 <= j <= NP and
C     1 <= k <= DP + 1 are set to zero. Next the nonzero (polynomial)
C     elements are read from the input file NIN. Each nonzero element is
C     given by the values i, j, d, P(i,j,k), k = 1, ..., d+1, where d is
C     the degree and P(i,j,k) is the coefficient of s**(k-1) in the
C     (i,j)-th element of P(s), i.e., let
C                                                              d
C         P   (s) = P   (0) + P   (1) * s + . . . + P   (d) * s
C          i,j       i,j       i,j                   i,j
C
C     be the nonzero (i,j)-th element of the matrix polynomial P(s).
C
C     Then P(i,j,k) corresponds to coefficient P   (k-1), k = 1,...,d+1.
C                                               i,j
C     For each nonzero element, the values i, j, and d are read as one
C     record of the file NIN, and the values P(i,j,k), k = 1,...,d+1,
C     are read as the following record.
C     The routine terminates after the last line has been read.
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
C     Based on routine RDSPOM by A.J. Geurts, Eindhoven University of
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
      INTEGER           D, I, J, K
C     .. External Subroutines ..
      EXTERNAL          DLASET, XERBLA
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
         CALL XERBLA( 'UD01CD', -INFO )
         RETURN
      END IF
C
      DO 10 K = 1, DP+1
         CALL DLASET( 'Full', MP, NP, ZERO, ZERO, P(1,1,K), LDP1 )
   10 CONTINUE
C
C     Read (i, j, d, P(i,j,k), k=1,...,d+1) of the nonzero elements one
C     by one.
C
   20 READ( NIN, FMT = *, END = 30 ) I, J, D
      IF ( I.LT.1 .OR. I.GT.MP .OR. J.LT.1 .OR. J.GT.NP .OR.
     $     D.LT.0 .OR. D.GT.(DP+1) ) THEN
         INFO = 1
         READ ( NIN, FMT = * )
      ELSE
         READ ( NIN, FMT = * ) ( P(I,J,K), K = 1, D+1 )
      END IF
      GO TO 20
C
   30 CONTINUE
      RETURN
C *** Last line of UD01CD ***
      END
