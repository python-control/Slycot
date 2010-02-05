      SUBROUTINE DG01ND( INDI, N, XR, XI, INFO )
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
C     To compute the discrete Fourier transform, or inverse Fourier
C     transform, of a real signal.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     INDI    CHARACTER*1
C             Indicates whether a Fourier transform or inverse Fourier
C             transform is to be performed as follows:
C             = 'D':  (Direct) Fourier transform;
C             = 'I':  Inverse Fourier transform.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             Half the number of real samples.  N must be a power of 2.
C             N >= 2.
C
C     XR      (input/output) DOUBLE PRECISION array, dimension (N+1)
C             On entry with INDI = 'D', the first N elements of this
C             array must contain the odd part of the input signal; for
C             example, XR(I) = A(2*I-1) for I = 1,2,...,N.
C             On entry with INDI = 'I', the first N+1 elements of this
C             array must contain the the real part of the input discrete
C             Fourier transform (computed, for instance, by a previous
C             call of the routine).
C             On exit with INDI = 'D', the first N+1 elements of this
C             array contain the real part of the output signal, that is
C             of the computed discrete Fourier transform.
C             On exit with INDI = 'I', the first N elements of this
C             array contain the odd part of the output signal, that is
C             of the computed inverse discrete Fourier transform.
C
C     XI      (input/output) DOUBLE PRECISION array, dimension (N+1)
C             On entry with INDI = 'D', the first N elements of this
C             array must contain the even part of the input signal; for
C             example, XI(I) = A(2*I) for I = 1,2,...,N.
C             On entry with INDI = 'I', the first N+1 elements of this
C             array must contain the the imaginary part of the input
C             discrete Fourier transform (computed, for instance, by a
C             previous call of the routine).
C             On exit with INDI = 'D', the first N+1 elements of this
C             array contain the imaginary part of the output signal,
C             that is of the computed discrete Fourier transform.
C             On exit with INDI = 'I', the first N elements of this
C             array contain the even part of the output signal, that is
C             of the computed inverse discrete Fourier transform.
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
C     Let A(1),....,A(2*N) be a real signal of 2*N samples. Then the
C     first N+1 samples of the discrete Fourier transform of this signal
C     are given by the formula:
C
C                  2*N           ((m-1)*(i-1))
C          FA(m) = SUM ( A(i) * W              ),
C                  i=1
C                                                  2
C     where m = 1,2,...,N+1, W = exp(-pi*j/N) and j = -1.
C
C     This transform can be computed as follows. First, transform A(i),
C     i = 1,2,...,2*N, into the complex signal Z(i) = (X(i),Y(i)),
C     i = 1,2,...,N. That is, X(i) = A(2*i-1) and Y(i) = A(2*i). Next,
C     perform a discrete Fourier transform on Z(i) by calling SLICOT
C     Library routine DG01MD. This gives a new complex signal FZ(k),
C     such that
C
C                   N            ((k-1)*(i-1))
C          FZ(k) = SUM ( Z(i) * V              ),
C                  i=1
C
C     where k = 1,2,...,N, V = exp(-2*pi*j/N).  Using the values of
C     FZ(k), the components of the discrete Fourier transform FA can be
C     computed by simple linear relations, implemented in the DG01NY
C     subroutine.
C
C     Finally, let
C
C          XR(k) = Re(FZ(k)), XI(k) = Im(FZ(k)),   k = 1,2,...,N,
C
C     be the contents of the arrays XR and XI on entry to DG01NY with
C     INDI = 'D', then on exit XR and XI contain the real and imaginary
C     parts of the Fourier transform of the original real signal A.
C     That is,
C
C          XR(m) = Re(FA(m)),  XI(m) = Im(FA(m)),
C
C     where m = 1,2,...,N+1.
C
C     If INDI = 'I', then the routine evaluates the inverse Fourier
C     transform of a complex signal which may itself be the discrete
C     Fourier transform of a real signal.
C
C     Let FA(m), m = 1,2,...,2*N, denote the full discrete Fourier
C     transform of a real signal A(i), i=1,2,...,2*N. The relationship
C     between FA and A is given by the formula:
C
C                 2*N            ((m-1)*(i-1))
C          A(i) = SUM ( FA(m) * W              ),
C                 m=1
C
C     where W = exp(pi*j/N).
C
C     Let
C
C          XR(m) = Re(FA(m)) and XI(m) = Im(FA(m)) for m = 1,2,...,N+1,
C
C     be the contents of the arrays XR and XI on entry to the routine
C     DG01NY with INDI = 'I', then on exit the first N samples of the
C     complex signal FZ are returned in XR and XI such that
C
C          XR(k) = Re(FZ(k)), XI(k) = Im(FZ(k)) and k = 1,2,...,N.
C
C     Next, an inverse Fourier transform is performed on FZ (e.g. by
C     calling SLICOT Library routine DG01MD), to give the complex signal
C     Z, whose i-th component is given by the formula:
C
C                  N             ((k-1)*(i-1))
C          Z(i) = SUM ( FZ(k) * V              ),
C                 k=1
C
C     where i = 1,2,...,N and V = exp(2*pi*j/N).
C
C     Finally, the 2*N samples of the real signal A can then be obtained
C     directly from Z. That is,
C
C          A(2*i-1) = Re(Z(i)) and A(2*i) = Im(Z(i)), for i = 1,2,...N.
C
C     Note that a discrete Fourier transform, followed by an inverse
C     transform will result in a signal which is a factor 2*N larger
C     than the original input signal.
C
C     REFERENCES
C
C     [1] Rabiner, L.R. and Rader, C.M.
C         Digital Signal Processing.
C         IEEE Press, 1972.
C
C     NUMERICAL ASPECTS
C
C     The algorithm requires 0( N*log(N) ) operations.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997.
C     Supersedes Release 2.0 routine DG01BD by R. Dekeyser, and
C     F. Dumortier, State University of Gent, Belgium.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Complex signals, digital signal processing, fast Fourier
C     transform, real signals.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER         INDI
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  XI(*), XR(*)
C     .. Local Scalars ..
      INTEGER           J
      LOGICAL           LINDI
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DG01MD, DG01NY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      INFO = 0
      LINDI = LSAME( INDI, 'D' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LINDI .AND. .NOT.LSAME( INDI, 'I' ) ) THEN
         INFO = -1
      ELSE
         J = 0
         IF( N.GE.2 ) THEN
            J = N
C           WHILE ( MOD( J, 2 ).EQ.0 ) DO
   10       CONTINUE
            IF ( MOD( J, 2 ).EQ.0 ) THEN
               J = J/2
               GO TO 10
            END IF
C           END WHILE 10
         END IF
         IF ( J.NE.1 ) INFO = -2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'DG01ND', -INFO )
         RETURN
      END IF
C
C     Compute the Fourier transform of Z = (XR,XI).
C
      IF ( .NOT.LINDI ) CALL DG01NY( INDI, N, XR, XI )
C
      CALL DG01MD( INDI, N, XR, XI, INFO )
C
      IF ( LINDI ) CALL DG01NY( INDI, N, XR, XI )
C
      RETURN
C *** Last line of DG01ND ***
      END
