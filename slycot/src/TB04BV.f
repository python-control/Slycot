      SUBROUTINE TB04BV( ORDER, P, M, MD, IGN, LDIGN, IGD, LDIGD, GN,
     $                   GD, D, LDD, TOL, INFO )
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
C     To separate the strictly proper part G0 from the constant part D
C     of an P-by-M proper transfer function matrix G.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     ORDER   CHARACTER*1
C             Specifies the order in which the polynomial coefficients
C             of the transfer function matrix are stored, as follows:
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
C             MD = MAX(IGD(I,J)) + 1.
C                  I,J
C
C     IGN     (input/output) INTEGER array, dimension (LDIGN,M)
C             On entry, the leading P-by-M part of this array must
C             contain the degrees of the numerator polynomials in G:
C             the (i,j) element of IGN must contain the degree of the
C             numerator polynomial of the polynomial ratio G(i,j).
C             On exit, the leading P-by-M part of this array contains
C             the degrees of the numerator polynomials in G0.
C
C     LDIGN   INTEGER
C             The leading dimension of array IGN.  LDIGN >= max(1,P).
C
C     IGD     (input) INTEGER array, dimension (LDIGD,M)
C             The leading P-by-M part of this array must contain the
C             degrees of the denominator polynomials in G (and G0):
C             the (i,j) element of IGD contains the degree of the
C             denominator polynomial of the polynomial ratio G(i,j).
C
C     LDIGD   INTEGER
C             The leading dimension of array IGD.  LDIGD >= max(1,P).
C
C     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD)
C             On entry, this array must contain the coefficients of the
C             numerator polynomials, Num(i,j), of the transfer function
C             matrix G. The polynomials are stored in a column-wise
C             order, i.e., Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2),
C             Num(2,2), ..., Num(P,2), ..., Num(1,M), Num(2,M), ...,
C             Num(P,M); MD memory locations are reserved for each
C             polynomial, hence, the (i,j) polynomial is stored starting
C             from the location ((j-1)*P+i-1)*MD+1. The coefficients
C             appear in increasing or decreasing order of the powers
C             of the indeterminate, according to ORDER.
C             On exit, this array contains the coefficients of the
C             numerator polynomials of the strictly proper part G0 of
C             the transfer function matrix G, stored similarly.
C
C     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD)
C             This array must contain the coefficients of the
C             denominator polynomials, Den(i,j), of the transfer
C             function matrix G. The polynomials are stored as for the
C             numerator polynomials.
C
C     D       (output) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array contains the
C             matrix D.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= max(1,P).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in determining the degrees of
C             the numerators Num0(i,j) of the strictly proper part of
C             the transfer function matrix G. If the user sets TOL > 0,
C             then the given value of TOL is used as an absolute
C             tolerance; the leading coefficients with absolute value
C             less than TOL are considered neglijible. If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by TOLDEF = IGN(i,j)*EPS*NORM( Num(i,j) ) is used
C             instead, where EPS is the machine precision (see LAPACK
C             Library routine DLAMCH), and NORM denotes the infinity
C             norm (the maximum coefficient in absolute value).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if the transfer function matrix is not proper;
C             = 2:  if a denominator polynomial is null.
C
C     METHOD
C
C     The (i,j) entry of the real matrix D is zero, if the degree of
C     Num(i,j), IGN(i,j), is less than the degree of Den(i,j), IGD(i,j),
C     and it is given by the ratio of the leading coefficients of
C     Num(i,j) and Den(i,j), if IGN(i,j) is equal to IGD(i,j),
C     for i = 1 : P, and for j = 1 : M.
C
C     FURTHER COMMENTS
C
C     For maximum efficiency of index calculations, GN and GD are
C     implemented as one-dimensional arrays.
C
C     CONTRIBUTORS
C
C     V. Sima, Research Institute for Informatics, Bucharest, May 2002.
C     Based on the BIMASC Library routine TMPRP by A. Varga.
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
      DOUBLE PRECISION   TOL
      INTEGER            INFO, LDD, LDIGD, LDIGN, M, MD, P
C     .. Array Arguments ..
      DOUBLE PRECISION   D(LDD,*), GD(*), GN(*)
      INTEGER            IGD(LDIGD,*), IGN(LDIGN,*)
C     .. Local Scalars ..
      LOGICAL            ASCEND
      INTEGER            I, II, J, K, KK, KM, ND, NN
      DOUBLE PRECISION   DIJ, EPS, TOLDEF
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, IDAMAX, LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
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
         CALL XERBLA( 'TB04BV', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( P, M ).EQ.0 )
     $   RETURN
C
C     Prepare the computation of the default tolerance.
C
      TOLDEF = TOL
      IF( TOLDEF.LE.ZERO )
     $   EPS = DLAMCH( 'Epsilon' )
C
      K = 1
C
      IF ( ASCEND ) THEN
C
C        Polynomial coefficients are stored in increasing order.
C
         DO 40 J = 1, M
C
            DO 30 I = 1, P
               NN = IGN(I,J)
               ND = IGD(I,J)
               IF ( NN.GT.ND ) THEN
C
C                 Error return: the transfer function matrix is
C                               not proper.
C
                  INFO = 1
                  RETURN
               ELSE IF ( NN.LT.ND .OR. ( ND.EQ.0 .AND. GN(K).EQ.ZERO ) )
     $               THEN
                  D(I,J) = ZERO
               ELSE
C
C                 Here NN = ND.
C
                  KK = K + NN
C
                  IF ( GD(KK).EQ.ZERO ) THEN
C
C                    Error return: the denominator is null.
C
                     INFO = 2
                     RETURN
                  ENDIF
C
                  DIJ    = GN(KK) / GD(KK)
                  D(I,J) = DIJ
                  GN(KK) = ZERO
                  IF ( NN.GT.0 ) THEN
                     CALL DAXPY( NN, -DIJ, GD(K), 1, GN(K), 1 )
                     IF ( TOL.LE.ZERO )
     $                  TOLDEF = DBLE( NN )*EPS*
     $                           ABS( GN(IDAMAX( NN, GN(K), 1 ) ) )
                     KM = NN
                     DO 10 II = 1, KM
                        KK = KK - 1
                        NN = NN - 1
                        IF ( ABS( GN(KK) ).GT.TOLDEF )
     $                     GO TO 20
   10                CONTINUE
C
   20                CONTINUE
C
                     IGN(I,J) = NN
                  ENDIF
               ENDIF
               K = K + MD
   30       CONTINUE
C
   40    CONTINUE
C
      ELSE
C
C        Polynomial coefficients are stored in decreasing order.
C
         DO 90 J = 1, M
C
            DO 80 I = 1, P
               NN = IGN(I,J)
               ND = IGD(I,J)
               IF ( NN.GT.ND ) THEN
C
C                 Error return: the transfer function matrix is
C                               not proper.
C
                  INFO = 1
                  RETURN
               ELSE IF ( NN.LT.ND .OR. ( ND.EQ.0 .AND. GN(K).EQ.ZERO ) )
     $               THEN
                  D(I,J) = ZERO
               ELSE
C
C                 Here NN = ND.
C
                  KK = K
C
                  IF ( GD(KK).EQ.ZERO ) THEN
C
C                    Error return: the denominator is null.
C
                     INFO = 2
                     RETURN
                  ENDIF
C
                  DIJ    = GN(KK) / GD(KK)
                  D(I,J) = DIJ
                  GN(KK) = ZERO
                  IF ( NN.GT.0 ) THEN
                     CALL DAXPY( NN, -DIJ, GD(K+1), 1, GN(K+1), 1 )
                     IF ( TOL.LE.ZERO )
     $                  TOLDEF = DBLE( NN )*EPS*
     $                           ABS( GN(IDAMAX( NN, GN(K+1), 1 ) ) )
                     KM = NN
                     DO 50 II = 1, KM
                        KK = KK + 1
                        NN = NN - 1
                        IF ( ABS( GN(KK) ).GT.TOLDEF )
     $                     GO TO 60
   50                CONTINUE
C
   60                CONTINUE
C
                     IGN(I,J) = NN
                     DO 70 II = 0, NN
                        GN(K+II) = GN(KK+II)
   70                CONTINUE
C
                  ENDIF
               ENDIF
               K = K + MD
   80       CONTINUE
C
   90    CONTINUE
C
      ENDIF
C
      RETURN
C *** Last line of TB04BV ***
      END
