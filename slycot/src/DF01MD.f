      SUBROUTINE DF01MD( SICO, N, DT, A, DWORK, INFO )
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
C     To compute the sine transform or cosine transform of a real
C     signal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SICO    CHARACTER*1
C             Indicates whether the sine transform or cosine transform
C             is to be computed as follows:
C             = 'S':  The sine transform is computed;
C             = 'C':  The cosine transform is computed.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of samples.  N must be a power of 2 plus 1.
C             N >= 5.
C
C     DT      (input) DOUBLE PRECISION
C             The sampling time of the signal.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry, this array must contain the signal to be
C             processed.
C             On exit, this array contains either the sine transform, if
C             SICO = 'S', or the cosine transform, if SICO = 'C', of the
C             given signal.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N+1)
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
C     Let A(1), A(2),..., A(N) be a real signal of N samples.
C
C     If SICO = 'S', the routine computes the sine transform of A as
C     follows. First, transform A(i), i = 1,2,...,N, into the complex
C     signal B(i), i = 1,2,...,(N+1)/2, where
C
C        B(1) = -2*A(2),
C        B(i) = {A(2i-2) - A(2i)} - j*A(2i-1) for i = 2,3,...,(N-1)/2,
C        B((N+1)/2) = 2*A(N-1) and j**2 = -1.
C
C     Next, perform a discrete inverse Fourier transform on B(i) by
C     calling SLICOT Library Routine DG01ND, to give the complex signal
C     Z(i), i = 1,2,...,(N-1)/2, from which the real signal C(i) may be
C     obtained as follows:
C
C        C(2i-1) = Re(Z(i)),  C(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2.
C
C     Finally, compute the sine transform coefficients S ,S ,...,S
C                                                       1  2      N
C     given by
C
C        S  = 0,
C         1
C                {                     [C(k) + C(N+1-k)]     }
C        S  = DT*{[C(k) - C(N+1-k)] - -----------------------},
C         k      {                    [2*sin(pi*(k-1)/(N-1))]}
C
C           for k = 2,3,...,N-1, and
C
C        S = 0.
C         N
C
C     If SICO = 'C', the routine computes the cosine transform of A as
C     follows. First, transform A(i), i = 1,2,...,N, into the complex
C     signal B(i), i = 1,2,...,(N+1)/2, where
C
C        B(1) = 2*A(1),
C        B(i) = 2*A(2i-1) + 2*j*{[A(2i-2) - A(2i)]}
C        for i = 2,3,...,(N-1)/2 and B((N+1)/2) = 2*A(N).
C
C     Next, perform a discrete inverse Fourier transform on B(i) by
C     calling SLICOT Library Routine DG01ND, to give the complex signal
C     Z(i), i = 1,2,...,(N-1)/2, from which the real signal D(i) may be
C     obtained as follows:
C
C        D(2i-1) = Re(Z(i)),  D(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2.
C
C     Finally, compute the cosine transform coefficients S ,S ,...,S
C                                                         1  2      N
C     given by
C
C        S  = 2*DT*[D(1) + A0],
C         1
C                {                     [D(k) - D(N+1-k)]     }
C        S  = DT*{[D(k) + D(N+1-k)] - -----------------------},
C         k      {                    [2*sin(pi*(k-1)/(N-1))]}
C
C
C           for k = 2,3,...,N-1, and
C
C        S  = 2*DT*[D(1) - A0],
C         N
C                 (N-1)/2
C     where A0 = 2*SUM   A(2i).
C                  i=1
C
C     REFERENCES
C
C     [1] Rabiner, L.R. and Rader, C.M.
C         Digital Signal Processing.
C         IEEE Press, 1972.
C
C     [2] Oppenheim, A.V. and Schafer, R.W.
C         Discrete-Time Signal Processing.
C         Prentice-Hall Signal Processing Series, 1989.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0( N*log(N) ) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine DF01AD by F. Dumortier, and
C     R.M.C. Dekeyser, State University of Gent, Belgium.
C
C     REVISIONS
C
C     V. Sima, Jan. 2003.
C
C     KEYWORDS
C
C     Digital signal processing, fast Fourier transform, complex
C     signals.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, FOUR
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                    FOUR = 4.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         SICO
      INTEGER           INFO, N
      DOUBLE PRECISION  DT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LSICO, LSIG
      INTEGER           I, I2, IND1, IND2, M, MD2
      DOUBLE PRECISION  A0, PIBYM, W1, W2, W3
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DG01ND, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, DBLE, MOD, SIN
C     .. Executable Statements ..
C
      INFO = 0
      LSICO = LSAME( SICO, 'S' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LSICO .AND. .NOT.LSAME( SICO, 'C' ) ) THEN
         INFO = -1
      ELSE
         M = 0
         IF( N.GT.4 ) THEN
            M = N - 1
C           WHILE ( MOD( M, 2 ).EQ.0 ) DO
   10       CONTINUE
            IF ( MOD( M, 2 ).EQ.0 ) THEN
               M = M/2
               GO TO 10
            END IF
C           END WHILE 10
         END IF
         IF ( M.NE.1 ) INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DF01MD', -INFO )
         RETURN
      END IF
C
C     Initialisation.
C
      M = N - 1
      MD2 = ( N + 1 )/2
      PIBYM = FOUR*ATAN( ONE )/DBLE( M )
      I2 = 1
      DWORK(MD2+1) = ZERO
      DWORK(2*MD2) = ZERO
C
      IF ( LSICO ) THEN
C
C        Sine transform.
C
         LSIG = .TRUE.
         DWORK(1)  = -TWO*A(2)
         DWORK(MD2) = TWO*A(M)
C
         DO 20 I = 4, M, 2
            I2 = I2 + 1
            DWORK(I2) = A(I-2) - A(I)
            DWORK(MD2+I2) = -A(I-1)
   20    CONTINUE
C
      ELSE
C
C        Cosine transform.
C
         LSIG = .FALSE.
         DWORK(1)   = TWO*A(1)
         DWORK(MD2) = TWO*A(N)
         A0 = A(2)
C
         DO 30 I = 4, M, 2
            I2 = I2 + 1
            DWORK(I2) = TWO*A(I-1)
            DWORK(MD2+I2) = TWO*( A(I-2) - A(I) )
            A0 = A0 + A(I)
   30    CONTINUE
C
         A0 = TWO*A0
      END IF
C
C     Inverse Fourier transform.
C
      CALL DG01ND( 'Inverse', MD2-1, DWORK(1), DWORK(MD2+1), INFO )
C
C     Sine or cosine coefficients.
C
      IF ( LSICO ) THEN
         A(1) = ZERO
         A(N) = ZERO
      ELSE
         A(1) = TWO*DT*( DWORK(1) + A0 )
         A(N) = TWO*DT*( DWORK(1) - A0 )
      END IF
C
      IND1 = MD2 + 1
      IND2 = N
C
      DO 40 I = 1, M - 1, 2
         W1 = DWORK(IND1)
         W2 = DWORK(IND2)
         IF ( LSIG ) W2 = -W2
         W3 = TWO*SIN( PIBYM*DBLE( I ) )
         A(I+1) = DT*( W1 + W2 - ( W1 - W2 )/W3 )
         IND1 = IND1 + 1
         IND2 = IND2 - 1
   40 CONTINUE
C
      IND1 = 2
      IND2 = MD2 - 1
C
      DO 50 I = 2, M - 2, 2
         W1 = DWORK(IND1)
         W2 = DWORK(IND2)
         IF ( LSIG ) W2 = -W2
         W3 = TWO*SIN( PIBYM*DBLE( I ) )
         A(I+1) = DT*( W1 + W2 - ( W1 - W2 )/W3 )
         IND1 = IND1 + 1
         IND2 = IND2 - 1
   50 CONTINUE
C
      RETURN
C *** Last line of DF01MD ***
      END
