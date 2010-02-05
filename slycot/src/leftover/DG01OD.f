      SUBROUTINE DG01OD( SCR, WGHT, N, A, W, INFO )
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
C     To compute the (scrambled) discrete Hartley transform of
C     a real signal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     SCR     CHARACTER*1
C             Indicates whether the signal is scrambled on input or
C             on output as follows:
C             = 'N':  the signal is not scrambled at all;
C             = 'I':  the input signal is bit-reversed;
C             = 'O':  the output transform is bit-reversed.
C
C     WGHT    CHARACTER*1
C             Indicates whether the precomputed weights are available
C             or not, as follows:
C             = 'A':  available;
C             = 'N':  not available.
C             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is
C             set to 'A' on exit.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             Number of real samples. N must be a power of 2.
C             N >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (N)
C             On entry with SCR = 'N' or SCR = 'O', this array must
C             contain the input signal.
C             On entry with SCR = 'I', this array must contain the
C             bit-reversed input signal.
C             On exit with SCR = 'N' or SCR = 'I', this array contains
C             the Hartley transform of the input signal.
C             On exit with SCR = 'O', this array contains the
C             bit-reversed Hartley transform.
C
C     W       (input/output) DOUBLE PRECISION array,
C                            dimension (N - LOG2(N))
C             On entry with WGHT = 'A', this array must contain the long
C             weight vector computed by a previous call of this routine
C             with the same value of N. If WGHT = 'N', the contents of
C             this array on entry is ignored.
C             On exit, this array contains the long weight vector.
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
C     This routine uses a Hartley butterfly algorithm as described
C     in [1].
C
C     REFERENCES
C
C     [1] Van Loan, Charles.
C         Computational frameworks for the fast Fourier transform.
C         SIAM, 1992.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable and requires O(N log(N))
C     floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, April 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000.
C
C     KEYWORDS
C
C     Digital signal processing, fast Hartley transform, real signals.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TWO, FOUR
      PARAMETER         ( ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         SCR, WGHT
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), W(*)
C     .. Local Scalars ..
      INTEGER           I, J, L, LEN, M, P1, P2, Q1, Q2, R1, R2, S1, S2,
     $                  WPOS
      LOGICAL           LFWD, LSCR, LWGHT
      DOUBLE PRECISION  CF, SF, T1, T2, TH
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN, COS, DBLE, MOD, SIN
C     .. Executable Statements ..
C
      INFO  = 0
      LFWD  = LSAME( SCR, 'N' ) .OR. LSAME( SCR, 'I' )
      LSCR  = LSAME( SCR, 'I' ) .OR. LSAME( SCR, 'O' )
      LWGHT = LSAME( WGHT, 'A' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.( LFWD .OR. LSCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LWGHT .AND. .NOT.LSAME( WGHT, 'N' ) ) THEN
         INFO = -2
      ELSE
         M = 0
         J = 0
         IF( N.GE.1 ) THEN
            J = N
C           WHILE ( MOD( J, 2 ).EQ.0 ) DO
   10       CONTINUE
            IF ( MOD( J, 2 ).EQ.0 ) THEN
               J = J/2
               M = M + 1
               GO TO 10
            END IF
C           END WHILE 10
            IF ( J.NE.1 ) INFO = -3
         ELSE IF ( N.LT.0 ) THEN
            INFO = -3
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DG01OD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.LE.1 )
     $   RETURN
C
      IF ( .NOT. LWGHT ) THEN
C
C        Compute the long weight vector via subvector scaling.
C
         R1  = 1
         LEN = 1
         TH  = FOUR*ATAN( ONE ) / DBLE( N )
C
         DO 30 L = 1, M - 2
            LEN = 2*LEN
            TH  = TWO*TH
            CF  = COS(TH)
            SF  = SIN(TH)
            W(R1)   = CF
            W(R1+1) = SF
            R1 = R1 + 2
C
            DO 20 I = 1, LEN - 2, 2
               W(R1)   = CF*W(I) - SF*W(I+1)
               W(R1+1) = SF*W(I) + CF*W(I+1)
               R1 = R1 + 2
   20       CONTINUE
C
   30    CONTINUE
C
         P1 = 3
         Q1 = R1 - 2
C
         DO 50 L = M - 2, 1, -1
C
            DO 40 I = P1, Q1, 4
               W(R1)   = W(I)
               W(R1+1) = W(I+1)
               R1 = R1 + 2
   40       CONTINUE
C
            P1 = Q1 + 4
            Q1 = R1 - 2
   50    CONTINUE
C
         WGHT = 'A'
C
      END IF
C
      IF ( LFWD .AND. .NOT.LSCR ) THEN
C
C        Inplace shuffling of data.
C
         J = 1
C
         DO 70 I = 1, N
            IF ( J.GT.I ) THEN
               T1   = A(I)
               A(I) = A(J)
               A(J) = T1
            END IF
            L = N/2
C           REPEAT
   60       IF ( J.GT.L ) THEN
               J = J - L
               L = L/2
               IF ( L.GE.2 ) GO TO 60
            END IF
C           UNTIL ( L.LT.2 )
            J = J + L
   70    CONTINUE
C
      END IF
C
      IF ( LFWD ) THEN
C
C        Compute Hartley transform with butterfly operators.
C
         DO 110 J = 2, N, 2
            T1 = A(J)
            A(J)   = A(J-1) - T1
            A(J-1) = A(J-1) + T1
  110    CONTINUE
C
         LEN  = 1
         WPOS = N - 2*M + 1
C
         DO 140 L = 1, M - 1
            LEN = 2*LEN
            P2  = 1
            Q2  = LEN + 1
            R2  = LEN / 2 + 1
            S2  = R2 + Q2 - 1
C
            DO 130 I = 0, N/( 2*LEN ) - 1
               T1    = A(Q2)
               A(Q2) = A(P2) - T1
               A(P2) = A(P2) + T1
               T1    = A(S2)
               A(S2) = A(R2) - T1
               A(R2) = A(R2) + T1
C
               P1 = P2 + 1
               Q1 = P1 + LEN
               R1 = Q1 - 2
               S1 = R1 + LEN
C
               DO 120 J = WPOS, WPOS + LEN - 3, 2
                  CF = W(J)
                  SF = W(J+1)
                  T1 =  CF*A(Q1) + SF*A(S1)
                  T2 = -CF*A(S1) + SF*A(Q1)
                  A(Q1) = A(P1) - T1
                  A(P1) = A(P1) + T1
                  A(S1) = A(R1) - T2
                  A(R1) = A(R1) + T2
                  P1 = P1 + 1
                  Q1 = Q1 + 1
                  R1 = R1 - 1
                  S1 = S1 - 1
  120          CONTINUE
C
               P2 = P2 + 2*LEN
               Q2 = Q2 + 2*LEN
               R2 = R2 + 2*LEN
               S2 = S2 + 2*LEN
  130       CONTINUE
C
            WPOS = WPOS - 2*LEN + 2
  140    CONTINUE
C
      ELSE
C
C        Compute Hartley transform with transposed butterfly operators.
C
         WPOS = 1
         LEN  = N
C
         DO 230 L = M - 1, 1, -1
            LEN = LEN / 2
            P2  = 1
            Q2  = LEN + 1
            R2  = LEN / 2 + 1
            S2  = R2 + Q2 - 1
C
            DO 220 I = 0, N/( 2*LEN ) - 1
               T1    = A(Q2)
               A(Q2) = A(P2) - T1
               A(P2) = A(P2) + T1
               T1    = A(S2)
               A(S2) = A(R2) - T1
               A(R2) = A(R2) + T1
C
               P1 = P2 + 1
               Q1 = P1 + LEN
               R1 = Q1 - 2
               S1 = R1 + LEN
C
               DO 210 J = WPOS, WPOS + LEN - 3, 2
                  CF = W(J)
                  SF = W(J+1)
                  T1 = A(P1) - A(Q1)
                  T2 = A(R1) - A(S1)
                  A(P1) = A(P1) + A(Q1)
                  A(R1) = A(R1) + A(S1)
                  A(Q1) =  CF*T1 + SF*T2
                  A(S1) = -CF*T2 + SF*T1
                  P1 = P1 + 1
                  Q1 = Q1 + 1
                  R1 = R1 - 1
                  S1 = S1 - 1
  210          CONTINUE
C
               P2 = P2 + 2*LEN
               Q2 = Q2 + 2*LEN
               R2 = R2 + 2*LEN
               S2 = S2 + 2*LEN
  220       CONTINUE
C
            WPOS = WPOS + LEN - 2
  230    CONTINUE
C
         DO 240 J = 2, N, 2
            T1 = A(J)
            A(J)   = A(J-1) - T1
            A(J-1) = A(J-1) + T1
  240    CONTINUE
C
      END IF
      RETURN
C *** Last line of DG01OD ***
      END
