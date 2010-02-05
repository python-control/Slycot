      SUBROUTINE MB04WR( JOB, TRANS, N, ILO, Q1, LDQ1, Q2, LDQ2, CS,
     $                   TAU, DWORK, LDWORK, INFO )
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
C     To generate orthogonal symplectic matrices U or V, defined as
C     products of symplectic reflectors and Givens rotators
C
C     U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) )
C         diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) )
C                              ....
C         diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ),
C
C     V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) )
C         diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) )
C                                   ....
C         diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ),
C
C     as returned by the SLICOT Library routines MB04TS or MB04TB. The
C     matrices U and V are returned in terms of their first N/2 rows:
C
C                 [  U1   U2 ]           [  V1   V2 ]
C             U = [          ],      V = [          ].
C                 [ -U2   U1 ]           [ -V2   V1 ]
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     JOB     CHARACTER*1
C             Specifies whether the matrix U or the matrix V is
C             required:
C             = 'U':  generate U;
C             = 'V':  generate V.
C
C     TRANS   CHARACTER*1
C             If  JOB = 'U'  then TRANS must have the same value as
C             the argument TRANA in the previous call of MB04TS or
C             MB04TB.
C             If  JOB = 'V'  then TRANS must have the same value as
C             the argument TRANB in the previous call of MB04TS or
C             MB04TB.
C
C     N       (input) INTEGER
C             The order of the matrices Q1 and Q2. N >= 0.
C
C     ILO     (input) INTEGER
C             ILO must have the same value as in the previous call of
C             MB04TS or MB04TB. U and V are equal to the unit matrix
C             except in the submatrices
C             U([ilo:n n+ilo:2*n], [ilo:n n+ilo:2*n]) and
C             V([ilo+1:n n+ilo+1:2*n], [ilo+1:n n+ilo+1:2*n]),
C             respectively.
C             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0.
C
C     Q1      (input/output) DOUBLE PRECISION array, dimension (LDQ1,N)
C             On entry, if  JOB = 'U'  and  TRANS = 'N'  then the
C             leading N-by-N part of this array must contain in its i-th
C             column the vector which defines the elementary reflector
C             FU(i).
C             If  JOB = 'U'  and  TRANS = 'T'  or  TRANS = 'C' then the
C             leading N-by-N part of this array must contain in its i-th
C             row the vector which defines the elementary reflector
C             FU(i).
C             If  JOB = 'V'  and  TRANS = 'N'  then the leading N-by-N
C             part of this array must contain in its i-th row the vector
C             which defines the elementary reflector FV(i).
C             If  JOB = 'V'  and  TRANS = 'T'  or  TRANS = 'C' then the
C             leading N-by-N part of this array must contain in its i-th
C             column the vector which defines the elementary reflector
C             FV(i).
C             On exit, if  JOB = 'U'  and  TRANS = 'N'  then the leading
C             N-by-N part of this array contains the matrix U1.
C             If  JOB = 'U'  and  TRANS = 'T'  or  TRANS = 'C' then the
C             leading N-by-N part of this array contains the matrix
C             U1**T.
C             If  JOB = 'V'  and  TRANS = 'N'  then the leading N-by-N
C             part of this array contains the matrix V1**T.
C             If  JOB = 'V'  and  TRANS = 'T'  or  TRANS = 'C' then the
C             leading N-by-N part of this array contains the matrix V1.
C
C     LDQ1    INTEGER
C             The leading dimension of the array Q1.  LDQ1 >= MAX(1,N).
C
C     Q2      (input/output) DOUBLE PRECISION array, dimension (LDQ2,N)
C             On entry, if  JOB = 'U'  then the leading N-by-N part of
C             this array must contain in its i-th column the vector
C             which defines the elementary reflector HU(i).
C             If  JOB = 'V'  then the leading N-by-N part of this array
C             must contain in its i-th row the vector which defines the
C             elementary reflector HV(i).
C             On exit, if  JOB = 'U'  then the leading N-by-N part of
C             this array contains the matrix U2.
C             If  JOB = 'V'  then the leading N-by-N part of this array
C             contains the matrix V2**T.
C
C     LDQ2    INTEGER
C             The leading dimension of the array Q2.  LDQ2 >= MAX(1,N).
C
C     CS      (input) DOUBLE PRECISION array, dimension (2N)
C             On entry, if  JOB = 'U'  then the first 2N elements of
C             this array must contain the cosines and sines of the
C             symplectic Givens rotators GU(i).
C             If  JOB = 'V'  then the first 2N-2 elements of this array
C             must contain the cosines and sines of the symplectic
C             Givens rotators GV(i).
C
C     TAU     (input) DOUBLE PRECISION array, dimension (N)
C             On entry, if  JOB = 'U'  then the first N elements of
C             this array must contain the scalar factors of the
C             elementary reflectors FU(i).
C             If  JOB = 'V'  then the first N-1 elements of this array
C             must contain the scalar factors of the elementary
C             reflectors FV(i).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0,  DWORK(1)  returns the optimal
C             value of LDWORK.
C             On exit, if  INFO = -12,  DWORK(1)  returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,2*(N-ILO+1)).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     REFERENCES
C
C     [1] Benner, P., Mehrmann, V., and Xu, H.
C         A numerically stable, structure preserving method for
C         computing the eigenvalues of real Hamiltonian or symplectic
C         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998.
C
C     [2] Kressner, D.
C         Block algorithms for orthogonal symplectic factorizations.
C         BIT, 43 (4), pp. 775-790, 2003.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSU).
C
C     KEYWORDS
C
C     Elementary matrix operations, Hamiltonian matrix, orthogonal
C     symplectic matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER         JOB, TRANS
      INTEGER           ILO, INFO, LDQ1, LDQ2, LDWORK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CS(*), DWORK(*), Q1(LDQ1,*), Q2(LDQ2,*), TAU(*)
C     .. Local Scalars ..
      LOGICAL           COMPU, LTRAN
      INTEGER           I, IERR, J, NH
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLASET, MB04WD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX
C
C     .. Executable Statements ..
C
C     Check the scalar input parameters.
C
      INFO  = 0
      LTRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      COMPU = LSAME( JOB, 'U' )
      IF ( .NOT.COMPU .AND. .NOT.LSAME( JOB, 'V' ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.LTRAN .AND. .NOT.LSAME( TRANS, 'N' ) ) THEN
         INFO = -2
      ELSE IF ( N.LT.0 ) THEN
         INFO = -3
      ELSE IF ( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LDQ1.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF ( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF ( LDWORK.LT.MAX( 1, 2*( N-ILO+1 ) ) ) THEN
         DWORK(1) = DBLE( MAX( 1, 2*( N-ILO+1 ) ) )
         INFO = -12
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04WR', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF ( COMPU ) THEN
         CALL DLASET( 'All', N, ILO-1, ZERO, ONE, Q1, LDQ1 )
         CALL DLASET( 'All', ILO-1, N-ILO+1, ZERO, ZERO, Q1(1,ILO),
     $                LDQ1 )
         CALL DLASET( 'All', N, ILO-1, ZERO, ZERO, Q2, LDQ2 )
         CALL DLASET( 'All', ILO-1, N-ILO+1, ZERO, ZERO, Q2(1,ILO),
     $                LDQ2 )
         NH = N - ILO + 1
      END IF
      IF ( COMPU .AND. .NOT.LTRAN ) THEN
C
C        Generate U1 and U2.
C
         IF ( NH.GT.0 ) THEN
            CALL MB04WD( 'No Transpose', 'No Transpose', NH, NH, NH,
     $                   Q1(ILO,ILO), LDQ1, Q2(ILO,ILO), LDQ2, CS(ILO),
     $                   TAU(ILO), DWORK, LDWORK, IERR )
         END IF
      ELSE IF ( COMPU.AND.LTRAN ) THEN
C
C        Generate U1**T and U2.
C
         IF ( NH.GT.0 ) THEN
            CALL MB04WD( 'Transpose', 'No Transpose', NH, NH, NH,
     $                   Q1(ILO,ILO), LDQ1, Q2(ILO,ILO), LDQ2, CS(ILO),
     $                   TAU(ILO), DWORK, LDWORK, IERR )
         END IF
      ELSE IF ( .NOT.COMPU .AND. .NOT.LTRAN ) THEN
C
C        Generate V1**T and V2**T.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the bottom, and set the first ilo rows and
C        columns to those of the unit matrix.
C
         DO 40  I = 1, N
            DO 10  J = N, MAX( I, ILO )+1, -1
               Q1(J,I) = ZERO
   10       CONTINUE
            DO 20  J = MAX( I, ILO ), ILO+1, -1
               Q1(J,I) = Q1(J-1,I)
   20       CONTINUE
            DO 30  J = ILO, 1, -1
               Q1(J,I) = ZERO
   30       CONTINUE
            IF ( I.LE.ILO )  Q1(I,I) = ONE
   40    CONTINUE
         DO 80  I = 1, N
            DO 50  J = N, MAX( I, ILO )+1, -1
               Q2(J,I) = ZERO
   50       CONTINUE
            DO 60  J = MAX( I, ILO ), ILO+1, -1
               Q2(J,I) = Q2(J-1,I)
   60       CONTINUE
            DO 70  J = ILO, 1, -1
               Q2(J,I) = ZERO
   70       CONTINUE
   80    CONTINUE
C
         NH = N - ILO
         IF ( NH.GT.0 ) THEN
            CALL MB04WD( 'Transpose', 'Transpose', NH, NH, NH,
     $                   Q1(ILO+1,ILO+1), LDQ1, Q2(ILO+1,ILO+1), LDQ2,
     $                   CS(ILO), TAU(ILO), DWORK, LDWORK, IERR )
         END IF
      ELSE IF ( .NOT.COMPU .AND. LTRAN ) THEN
C
C        Generate V1 and V2**T.
C
C        Shift the vectors which define the elementary reflectors one
C        column to the right/bottom, and set the first ilo rows and
C        columns to those of the unit matrix.
C
         DO 110  J = N, ILO + 1, -1
            DO 90  I = 1, J-1
               Q1(I,J) = ZERO
   90       CONTINUE
            DO 100  I = J+1, N
               Q1(I,J) = Q1(I,J-1)
  100       CONTINUE
  110    CONTINUE
         CALL DLASET( 'All', N, ILO, ZERO, ONE, Q1, LDQ1 )
         DO 150 I = 1, N
            DO 120  J = N, MAX( I, ILO )+1, -1
               Q2(J,I) = ZERO
  120       CONTINUE
            DO 130  J = MAX( I, ILO ), ILO+1, -1
               Q2(J,I) = Q2(J-1,I)
  130       CONTINUE
            DO 140  J = ILO, 1, -1
               Q2(J,I) = ZERO
  140       CONTINUE
  150    CONTINUE
         NH = N - ILO
C
         IF ( NH.GT.0 ) THEN
            CALL MB04WD( 'No Transpose', 'Transpose', NH, NH, NH,
     $                   Q1(ILO+1,ILO+1), LDQ1, Q2(ILO+1,ILO+1), LDQ2,
     $                   CS(ILO), TAU(ILO), DWORK, LDWORK, IERR )
         END IF
      END IF
      RETURN
C *** Last line of MB04WR ***
      END
