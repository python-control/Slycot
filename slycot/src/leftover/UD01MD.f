      SUBROUTINE UD01MD( M, N, L, NOUT, A, LDA, TEXT, INFO )
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
C     To print an M-by-N real matrix A row by row. The elements of A
C     are output to 7 significant figures.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of matrix A to be printed.  M >= 1.
C
C     N       (input) INTEGER
C             The number of columns of matrix A to be printed.  N >= 1.
C
C     L       (input) INTEGER
C             The number of elements of matrix A to be printed per line.
C             1 <= L <= 5.
C
C     NOUT    (input) INTEGER
C             The output channel to which the results are sent.
C             NOUT >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading M-by-N part of this array must contain the
C             matrix to be printed.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= M.
C
C     TEXT    (input) CHARACTER*72.
C             Title caption of the matrix to be printed (up to a
C             maximum of 72 characters). For example, TEXT = 'Matrix A'.
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
C     The routine first prints the contents of TEXT as a title, followed
C     by the elements of the matrix A such that
C
C     (i)  if N <= L, the leading M-by-N part is printed;
C     (ii) if N = k*L + p (where k,p > 0), then k M-by-L blocks of
C          consecutive columns of A are printed one after another
C          followed by one M-by-p block containing the last p columns
C          of A.
C
C     Row numbers are printed on the left of each row and a column
C     number appears on top of each column.
C     The routine uses 2 + (k + 1)*(m + 1) lines and 8 + 15*c positions
C     per line where c is the actual number of columns, (i.e. c = L
C     or c = p).
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
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1997.
C     Supersedes Release 2.0 routine UD01AD by H. Willemsen, Eindhoven
C     University of Technology, Holland.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2009.
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           INFO, L, LDA, M, N, NOUT
      CHARACTER*(*)     TEXT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           I, J, J1, J2, JJ, LENTXT, LTEXT, N1
C     .. External Subroutines ..
      EXTERNAL          XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( M.LT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.1 ) THEN
         INFO = -2
      ELSE IF( L.LT.1 .OR. L.GT.5 ) THEN
         INFO = -3
      ELSE IF( NOUT.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.M ) THEN
         INFO = -6
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'UD01MD', -INFO )
         RETURN
      END IF
C
      LENTXT = LEN( TEXT )
C
      DO 20 LTEXT = MIN( 72, LENTXT ), 2, -1
         IF ( TEXT(LTEXT:LTEXT).NE.' ' ) GO TO 40
   20 CONTINUE
C
   40 CONTINUE
      WRITE ( NOUT, FMT=99996 ) TEXT(1:LTEXT), M, N
      N1 = ( N-1 )/L
      J1 = 1
      J2 = L
C
      DO 80 J = 1, N1
         WRITE ( NOUT, FMT=99999 ) ( JJ, JJ=J1, J2 )
C
         DO 60 I = 1, M
            WRITE ( NOUT, FMT=99997 ) I, ( A(I,JJ), JJ=J1, J2 )
   60    CONTINUE
C
         WRITE ( NOUT, FMT=99998 )
         J1 = J1 + L
         J2 = J2 + L
   80 CONTINUE
C
      WRITE ( NOUT, FMT=99999 ) ( J, J=J1, N )
C
      DO 100 I = 1, M
         WRITE ( NOUT, FMT=99997 ) I, ( A(I,JJ), JJ=J1, N )
  100 CONTINUE
C
      WRITE ( NOUT, FMT=99998 )
C
      RETURN
C
99999 FORMAT (8X,5(5X,I5,5X) )
99998 FORMAT (' ' )
99997 FORMAT (1X,I5,2X,5D15.7 )
99996 FORMAT (1X,A,' (',I5,'X',I5,')',/ )
C *** Last line of UD01MD ***
      END
