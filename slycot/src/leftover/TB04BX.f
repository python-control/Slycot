      SUBROUTINE TB04BX( IP, IZ, A, LDA, B, C, D, PR, PI, ZR, ZI, GAIN,
     $                   IWORK )
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
C     To compute the gain of a single-input single-output linear system,
C     given its state-space representation (A,b,c,d), and its poles and
C     zeros. The matrix A is assumed to be in an upper Hessenberg form.
C     The gain is computed using the formula
C
C                          -1         IP              IZ
C        g = (c*( S0*I - A ) *b + d)*Prod( S0 - Pi )/Prod( S0 - Zi ) ,
C                                     i=1             i=1            (1)
C
C     where Pi, i = 1 : IP, and Zj, j = 1 : IZ, are the poles and zeros,
C     respectively, and S0 is a real scalar different from all poles and
C     zeros.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     IP      (input) INTEGER
C             The number of the system poles.  IP >= 0.
C
C     IZ      (input) INTEGER
C             The number of the system zeros.  IZ >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,IP)
C             On entry, the leading IP-by-IP part of this array must
C             contain the state dynamics matrix A in an upper Hessenberg
C             form. The elements below the second diagonal are not
C             referenced.
C             On exit, the leading IP-by-IP upper Hessenberg part of
C             this array contains the LU factorization of the matrix
C             A - S0*I, as computed by SLICOT Library routine MB02SD.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= max(1,IP).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (IP)
C             On entry, this array must contain the system input
C             vector b.
C             On exit, this array contains the solution of the linear
C             system ( A - S0*I )x = b .
C
C     C       (input) DOUBLE PRECISION array, dimension (IP)
C             This array must contain the system output vector c.
C
C     D       (input) DOUBLE PRECISION
C             The variable must contain the system feedthrough scalar d.
C
C     PR      (input) DOUBLE PRECISION array, dimension (IP)
C             This array must contain the real parts of the system
C             poles. Pairs of complex conjugate poles must be stored in
C             consecutive memory locations.
C
C     PI      (input) DOUBLE PRECISION array, dimension (IP)
C             This array must contain the imaginary parts of the system
C             poles.
C
C     ZR      (input) DOUBLE PRECISION array, dimension (IZ)
C             This array must contain the real parts of the system
C             zeros. Pairs of complex conjugate zeros must be stored in
C             consecutive memory locations.
C
C     ZI      (input) DOUBLE PRECISION array, dimension (IZ)
C             This array must contain the imaginary parts of the system
C             zeros.
C
C     GAIN    (output) DOUBLE PRECISION
C             The gain of the linear system (A,b,c,d), given by (1).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (IP)
C             On exit, it contains the pivot indices; for 1 <= i <= IP,
C             row i of the matrix A - S0*I was interchanged with
C             row IWORK(i).
C
C     METHOD
C
C     The routine implements the method presented in [1]. A suitable
C     value of S0 is chosen based on the system poles and zeros.
C     Then, the LU factorization of the upper Hessenberg, nonsingular
C     matrix A - S0*I is computed and used to solve the linear system
C     in (1).
C
C     REFERENCES
C
C     [1] Varga, A. and Sima, V.
C         Numerically Stable Algorithm for Transfer Function Matrix
C         Evaluation.
C         Int. J. Control, vol. 33, nr. 6, pp. 1123-1133, 1981.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable in practice and requires
C     O(IP*IP) floating point operations.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2002.
C     Partly based on the BIMASC Library routine GAIN by A. Varga.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalue, state-space representation, transfer function, zeros.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, P1, ONEP1
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                     P1 = 0.1D0, ONEP1 = 1.1D0 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   D, GAIN
      INTEGER            IP, IZ, LDA
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), B(*), C(*), PI(*), PR(*), ZI(*),
     $                   ZR(*)
      INTEGER            IWORK(*)
C     .. Local Scalars ..
      INTEGER            I, INFO
      DOUBLE PRECISION   S0, S
C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
C     .. External Subroutines ..
      EXTERNAL           MB02RD, MB02SD
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     ..
C     .. Executable Statements ..
C
C     For efficiency, the input scalar parameters are not checked.
C
C     Quick return if possible.
C
      IF( IP.EQ.0 ) THEN
         GAIN = ZERO
         RETURN
      END IF
C
C     Compute a suitable value for S0 .
C
      S0 = ZERO
C
      DO 10 I = 1, IP
         S = ABS( PR(I) )
         IF ( PI(I).NE.ZERO )
     $      S = S + ABS( PI(I) )
         S0 = MAX( S0, S )
   10 CONTINUE
C
      DO 20 I = 1, IZ
         S = ABS( ZR(I) )
         IF ( ZI(I).NE.ZERO )
     $      S = S + ABS( ZI(I) )
         S0 = MAX( S0, S )
   20 CONTINUE
C
      S0 = TWO*S0 + P1
      IF ( S0.LE.ONE )
     $   S0 = ONEP1
C
C     Form A - S0*I .
C
      DO 30 I = 1, IP
         A(I,I) = A(I,I) - S0
   30 CONTINUE
C
C     Compute the LU factorization of the matrix A - S0*I
C     (guaranteed to be nonsingular).
C
      CALL MB02SD( IP, A, LDA, IWORK, INFO )
C
C     Solve the linear system (A - S0*I)*x = b .
C
      CALL MB02RD( 'No Transpose', IP, 1, A, LDA, IWORK, B, IP, INFO )
C                        -1
C     Compute c*(S0*I - A) *b + d .
C
      GAIN = D - DDOT( IP, C, 1, B, 1 )
C
C     Multiply by the products in terms of poles and zeros in (1).
C
      I = 1
C
C     WHILE ( I <= IP ) DO
C
   40 IF ( I.LE.IP ) THEN
         IF ( PI(I).EQ.ZERO ) THEN
            GAIN = GAIN*( S0 - PR(I) )
            I = I + 1
         ELSE
            GAIN = GAIN*( S0*( S0  - TWO*PR(I) ) + PR(I)**2 + PI(I)**2 )
            I = I + 2
         END IF
         GO TO 40
      END IF
C
C     END WHILE 40
C
      I = 1
C
C     WHILE ( I <= IZ ) DO
C
   50 IF ( I.LE.IZ ) THEN
         IF ( ZI(I).EQ.ZERO ) THEN
            GAIN = GAIN/( S0 - ZR(I) )
            I = I + 1
         ELSE
            GAIN = GAIN/( S0*( S0  - TWO*ZR(I) ) + ZR(I)**2 + ZI(I)**2 )
            I = I + 2
         END IF
         GO TO 50
      END IF
C
C     END WHILE 50
C
      RETURN
C *** Last line of TB04BX ***
      END
