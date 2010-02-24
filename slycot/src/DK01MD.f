      SUBROUTINE DK01MD( TYPE, N, A, INFO )
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
C     To apply an anti-aliasing window to a real signal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPE    CHARACTER*1
C             Indicates the type of window to be applied to the signal
C             as follows:
C             = 'M':  Hamming window;
C             = 'N':  Hann window;
C             = 'Q':  Quadratic window.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of samples.  N >= 1.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the signal to be
C             processed.
C             On exit, this array contains the windowing function.
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
C     If TYPE = 'M', then a Hamming window is applied to A(1),...,A(N),
C     which yields
C       _
C       A(i) = (0.54 + 0.46*cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N.
C
C     If TYPE = 'N', then a Hann window is applied to A(1),...,A(N),
C     which yields
C       _
C       A(i) = 0.5*(1 + cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N.
C
C     If TYPE = 'Q', then a quadratic window is applied to A(1),...,
C     A(N), which yields
C       _
C       A(i) = (1 - 2*((i-1)/(N-1))**2)*(1 - (i-1)/(N-1))*A(i),
C                                             i = 1,2,...,(N-1)/2+1;
C       _
C       A(i) = 2*(1 - ((i-1)/(N-1))**3)*A(i), i = (N-1)/2+2,...,N.
C
C     REFERENCES
C
C     [1] Rabiner, L.R. and Rader, C.M.
C         Digital Signal Processing.
C         IEEE Press, 1972.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0( N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine DK01AD by R. Dekeyser, State
C     University of Gent, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Digital signal processing, Hamming window, Hann window, real
C     signals, windowing.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  PT46, HALF, PT54, ONE, TWO, FOUR
      PARAMETER         ( PT46=0.46D0, HALF=0.5D0, PT54=0.54D0,
     $                    ONE = 1.0D0, TWO=2.0D0, FOUR=4.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         TYPE
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*)
C     .. Local Scalars ..
      LOGICAL           MTYPE, MNTYPE, NTYPE
      INTEGER           I, N1
      DOUBLE PRECISION  BUF, FN, TEMP
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, COS, DBLE
C     .. Executable Statements ..
C
      INFO = 0
      MTYPE = LSAME( TYPE, 'M' )
      NTYPE = LSAME( TYPE, 'N' )
      MNTYPE = MTYPE.OR.NTYPE
C
C     Test the input scalar arguments.
C
      IF( .NOT.MNTYPE .AND. .NOT.LSAME( TYPE, 'Q' ) )
     $   THEN
         INFO = -1
      ELSE IF( N.LE.0 ) THEN
         INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DK01MD', -INFO )
         RETURN
      END IF
C
      FN = DBLE( N-1 )
      IF( MNTYPE ) TEMP = FOUR*ATAN( ONE )/FN
C
      IF ( MTYPE ) THEN
C
C        Hamming window.
C
         DO 10 I = 1, N
            A(I) = A(I)*( PT54 + PT46*COS( TEMP*DBLE( I-1 ) ) )
   10    CONTINUE
C
      ELSE IF ( NTYPE ) THEN
C
C        Hann window.
C
         DO 20 I = 1, N
            A(I) = A(I)*HALF*( ONE + COS( TEMP*DBLE( I-1 ) ) )
   20    CONTINUE
C
      ELSE
C
C        Quadratic window.
C
         N1 = ( N-1 )/2 + 1
C
         DO 30 I = 1, N
            BUF = DBLE( I-1 )/FN
            TEMP = BUF**2
            IF ( I.LE.N1 ) THEN
               A(I) = A(I)*( ONE - TWO*TEMP )*( ONE - BUF )
            ELSE
               A(I) = A(I)*TWO*( ONE - BUF*TEMP )
            END IF
   30    CONTINUE
C
      END IF
C
      RETURN
C *** Last line of DK01MD ***
      END
