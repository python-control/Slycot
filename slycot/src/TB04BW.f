      SUBROUTINE TB04BW( ORDER, P, M, MD, IGN, LDIGN, IGD, LDIGD, GN,
     $                   GD, D, LDD, INFO )
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
C     To compute the sum of an P-by-M rational matrix G and a real
C     P-by-M matrix D.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ORDER   CHARACTER*1
C             Specifies the order in which the polynomial coefficients
C             of the rational matrix are stored, as follows:
C             = 'I':  Increasing order of powers of the indeterminate;
C             = 'D':  Decreasing order of powers of the indeterminate.
C
C     Input/Output Parameters
C
C     P       (input) INTEGER
C             The number of the system outputs.  P >= 0.
C
C     M       (input) INTEGER
C             The number of the system inputs.  M >= 0.
C
C     MD      (input) INTEGER
C             The maximum degree of the polynomials in G, plus 1, i.e.,
C             MD = MAX(IGN(I,J),IGD(I,J)) + 1.
C                  I,J
C
C     IGN     (input/output) INTEGER array, dimension (LDIGN,M)
C             On entry, the leading P-by-M part of this array must
C             contain the degrees of the numerator polynomials in G:
C             the (i,j) element of IGN must contain the degree of the
C             numerator polynomial of the polynomial ratio G(i,j).
C             On exit, the leading P-by-M part of this array contains
C             the degrees of the numerator polynomials in G + D.
C
C     LDIGN   INTEGER
C             The leading dimension of array IGN.  LDIGN >= max(1,P).
C
C     IGD     (input) INTEGER array, dimension (LDIGD,M)
C             The leading P-by-M part of this array must contain the
C             degrees of the denominator polynomials in G (and G + D):
C             the (i,j) element of IGD contains the degree of the
C             denominator polynomial of the polynomial ratio G(i,j).
C
C     LDIGD   INTEGER
C             The leading dimension of array IGD.  LDIGD >= max(1,P).
C
C     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD)
C             On entry, this array must contain the coefficients of the
C             numerator polynomials, Num(i,j), of the rational matrix G.
C             The polynomials are stored in a column-wise order, i.e.,
C             Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), Num(2,2),
C             ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., Num(P,M);
C             MD memory locations are reserved for each polynomial,
C             hence, the (i,j) polynomial is stored starting from the
C             location ((j-1)*P+i-1)*MD+1. The coefficients appear in
C             increasing or decreasing order of the powers of the
C             indeterminate, according to ORDER.
C             On exit, this array contains the coefficients of the
C             numerator polynomials of the rational matrix G + D,
C             stored similarly.
C
C     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD)
C             This array must contain the coefficients of the
C             denominator polynomials, Den(i,j), of the rational
C             matrix G. The polynomials are stored as for the
C             numerator polynomials.
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             matrix D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= max(1,P).
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
C     The (i,j) entry of the real matrix D is added to the (i,j) entry
C     of the matrix G, g(i,j), which is a ratio of two polynomials,
C     for i = 1 : P, and for j = 1 : M. If g(i,j) = 0, it is assumed
C     that its denominator is 1.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically stable.
C
C     FURTHER COMMENTS
C
C     Often, the rational matrix G is found from a state-space
C     representation (A,B,C), and D corresponds to the direct
C     feedthrough matrix of the system. The sum G + D gives the
C     transfer function matrix of the system (A,B,C,D).
C     For maximum efficiency of index calculations, GN and GD are
C     implemented as one-dimensional arrays.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2002.
C     Based on the BIMASC Library routine TMCADD by A. Varga.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004.
C
C     KEYWORDS
C
C     State-space representation, transfer function.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          ORDER
      INTEGER            INFO, LDD, LDIGD, LDIGN, M, MD, P
C     .. Array Arguments ..
      DOUBLE PRECISION   D(LDD,*), GD(*), GN(*)
      INTEGER            IGD(LDIGD,*), IGN(LDIGN,*)
C     .. Local Scalars ..
      LOGICAL            ASCEND
      INTEGER            I, II, J, K, KK, KM, ND, NN
      DOUBLE PRECISION   DIJ
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      INFO   = 0
      ASCEND = LSAME( ORDER, 'I' )
      IF( .NOT.ASCEND .AND. .NOT.LSAME( ORDER, 'D' ) ) THEN
         INFO = -1
      ELSE IF( P.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( MD.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDIGN.LT.MAX( 1, P ) ) THEN
         INFO = -6
      ELSE IF( LDIGD.LT.MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TB04BW', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( P, M ).EQ.0 )
     $   RETURN
C
      K = 1
C
      IF ( ASCEND ) THEN
C
C        Polynomial coefficients are stored in increasing order.
C
         DO 30 J = 1, M
C
            DO 20 I = 1, P
               DIJ = D(I,J)
               IF ( DIJ.NE.ZERO ) THEN
                  NN = IGN(I,J)
                  ND = IGD(I,J)
                  IF ( NN.EQ.0 .AND. ND.EQ.0 ) THEN
                     IF ( GN(K).EQ.ZERO ) THEN
                        GN(K) = DIJ
                     ELSE
                        GN(K) = GN(K) + DIJ*GD(K)
                     ENDIF
                  ELSE
                     KM = MIN( NN, ND ) + 1
                     CALL DAXPY( KM, DIJ, GD(K), 1, GN(K), 1 )
                     IF ( NN.LT.ND ) THEN
C
                        DO 10 II = K + KM, K + ND
                           GN(II) = DIJ*GD(II)
   10                   CONTINUE
C
                        IGN(I,J) = ND
                      ENDIF
                  ENDIF
               ENDIF
               K = K + MD
   20       CONTINUE
C
   30    CONTINUE
C
      ELSE
C
C        Polynomial coefficients are stored in decreasing order.
C
         DO 60 J = 1, M
C
            DO 50 I = 1, P
               DIJ = D(I,J)
               IF ( DIJ.NE.ZERO ) THEN
                  NN = IGN(I,J)
                  ND = IGD(I,J)
                  IF ( NN.EQ.0 .AND. ND.EQ.0 ) THEN
                     IF ( GN(K).EQ.ZERO ) THEN
                        GN(K) = DIJ
                     ELSE
                        GN(K) = GN(K) + DIJ*GD(K)
                     ENDIF
                  ELSE
                     KM = MIN( NN, ND ) + 1
                     IF ( NN.LT.ND ) THEN
                        KK = K + ND - NN
C
                        DO 35 II = K + NN, K, -1
                           GN(II+ND-NN) = GN(II)
   35                   CONTINUE
C
                        DO 40 II = K, KK - 1
                           GN(II) = DIJ*GD(II)
   40                   CONTINUE
C
                        IGN(I,J) = ND
                        CALL DAXPY( KM, DIJ, GD(KK), 1, GN(KK), 1 )
                     ELSE
                        KK = K + NN - ND
                        CALL DAXPY( KM, DIJ, GD(K), 1, GN(KK), 1 )
                     ENDIF
                  ENDIF
               ENDIF
               K = K + MD
   50       CONTINUE
C
   60    CONTINUE
C
      ENDIF
C
      RETURN
C *** Last line of TB04BW ***
      END
