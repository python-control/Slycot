      SUBROUTINE TF01QD( NC, NB, N, IORD, AR, MA, H, LDH, INFO )
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
C     To compute N Markov parameters M(1), M(2),..., M(N) from a
C     multivariable system whose transfer function matrix G(z) is given.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NC      (input) INTEGER
C             The number of system outputs, i.e. the number of rows in
C             the transfer function matrix G(z).  NC >= 0.
C
C     NB      (input) INTEGER
C             The number of system inputs, i.e. the number of columns in
C             the transfer function matrix G(z).  NB >= 0.
C
C     N       (input) INTEGER
C             The number of Markov parameters M(k) to be computed.
C             N >= 0.
C
C     IORD    (input) INTEGER array, dimension (NC*NB)
C             This array must contain the order r of the elements of the
C             transfer function matrix G(z), stored row by row.
C             For example, the order of the (i,j)-th element of G(z) is
C             given by IORD((i-1)xNB+j).
C
C     AR      (input) DOUBLE PRECISION array, dimension (NA), where
C             NA = IORD(1) + IORD(2) + ... + IORD(NC*NB).
C             The leading NA elements of this array must contain the
C             denominator coefficients AR(1),...,AR(r) in equation (1)
C             of the (i,j)-th element of the transfer function matrix
C             G(z), stored row by row, i.e. in the order
C             (1,1),(1,2),...,(1,NB), (2,1),(2,2),...,(2,NB), ...,
C             (NC,1),(NC,2),...,(NC,NB). The coefficients must be given
C             in decreasing order of powers of z; the coefficient of the
C             highest order term is assumed to be equal to 1.
C
C     MA      (input) DOUBLE PRECISION array, dimension (NA)
C             The leading NA elements of this array must contain the
C             numerator coefficients MA(1),...,MA(r) in equation (1)
C             of the (i,j)-th element of the transfer function matrix
C             G(z), stored row by row, i.e. in the order
C             (1,1),(1,2),...,(1,NB), (2,1),(2,2),...,(2,NB), ...,
C             (NC,1),(NC,2),...,(NC,NB). The coefficients must be given
C             in decreasing order of powers of z.
C
C     H       (output) DOUBLE PRECISION array, dimension (LDH,N*NB)
C             The leading NC-by-N*NB part of this array contains the
C             multivariable Markov parameter sequence M(k), where each
C             parameter M(k) is an NC-by-NB matrix and k = 1,2,...,N.
C             The Markov parameters are stored such that H(i,(k-1)xNB+j)
C             contains the (i,j)-th element of M(k) for i = 1,2,...,NC
C             and j = 1,2,...,NB.
C
C     LDH     INTEGER
C             The leading dimension of array H.  LDH >= MAX(1,NC).
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
C     The (i,j)-th element of G(z), defining the particular I/O transfer
C     between output i and input j, has the following form:
C
C                          -1         -2               -r
C                    MA(1)z   + MA(2)z   + ... + MA(r)z
C         G  (z) = ----------------------------------------.         (1)
C          ij                -1         -2               -r
C                  1 + AR(1)z   + AR(2)z   + ... + AR(r)z
C
C     The (i,j)-th element of G(z) is defined by its order r, its r
C     moving average coefficients (= numerator) MA(1),...,MA(r) and its
C     r autoregressive coefficients (= denominator) AR(1),...,AR(r). The
C     coefficient of the constant term in the denominator is assumed to
C     be equal to 1.
C
C     The relationship between the (i,j)-th element of the Markov
C     parameters M(1),M(2),...,M(N) and the corresponding element of the
C     transfer function matrix G(z) is given by:
C
C                               -1          -2                -k
C      G  (z) = M  (0) + M  (1)z   + M  (2)z   + ... + M  (k)z  + ...(2)
C       ij       ij       ij          ij                ij
C
C     Equating (1) and (2), we find that the relationship between the
C     (i,j)-th element of the Markov parameters M(k) and the ARMA
C     parameters AR(1),...,AR(r) and MA(1),...,MA(r) of the (i,j)-th
C     element of the transfer function matrix G(z) is as follows:
C
C        M  (1)   = MA(1),
C         ij
C                           k-1
C        M  (k)   = MA(k) - SUM AR(p) x M  (k-p) for 1 < k <= r and
C         ij                p=1          ij
C                      r
C        M  (k+r) = - SUM AR(p) x M  (k+r-p) for k > 0.
C         ij          p=1          ij
C
C     From these expressions the Markov parameters M(k) are computed
C     element by element.
C
C     REFERENCES
C
C     [1] Luenberger, D.G.
C         Introduction to Dynamic Systems: Theory, Models and
C         Applications.
C         John Wiley & Sons, New York, 1979.
C
C     NUMERICAL ASPECTS
C
C     The computation of the (i,j)-th element of M(k) requires:
C        (k-1) multiplications and k additions if k <= r;
C          r   multiplications and r additions if k > r.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TF01ED by S. Van Huffel, Katholieke
C     Univ. Leuven, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Markov parameters, multivariable system, transfer function,
C     transfer matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, N, NB, NC
C     .. Array Arguments ..
      INTEGER           IORD(*)
      DOUBLE PRECISION  AR(*), H(LDH,*), MA(*)
C     .. Local Scalars ..
      INTEGER           I, J, JJ, JK, K, KI, LDHNB, NL, NORD
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( NC.LT.0 ) THEN
         INFO = -1
      ELSE IF( NB.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDH.LT.MAX( 1, NC ) ) THEN
         INFO = -8
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TF01QD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( MAX( NC, NB, N ).EQ.0 )
     $   RETURN
C
      LDHNB = LDH*NB
      NL = 1
      K = 1
C
      DO 60 I = 1, NC
C
         DO 50 J = 1, NB
            NORD = IORD(K)
            H(I,J) = MA(NL)
            JK = J
C
            DO 20 KI = 1, NORD - 1
               JK = JK + NB
               H(I,JK) = MA(NL+KI) - DDOT( KI, AR(NL), 1, H(I,J),
     $                                     -LDHNB )
   20       CONTINUE
C
            DO 40 JJ = J, J + (N - NORD - 1)*NB, NB
               JK = JK + NB
               H(I,JK) = -DDOT( NORD, AR(NL), 1, H(I,JJ), -LDHNB )
   40       CONTINUE
C
            NL = NL + NORD
            K = K + 1
   50    CONTINUE
C
   60 CONTINUE
C
      RETURN
C *** Last line of TF01QD ***
      END
