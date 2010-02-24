      SUBROUTINE MB04MD( N, MAXRED, A, LDA, SCALE, INFO )
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
C     To reduce the 1-norm of a general real matrix A by balancing.
C     This involves diagonal similarity transformations applied
C     iteratively to A to make the rows and columns as close in norm as
C     possible.
C
C     This routine can be used instead LAPACK Library routine DGEBAL,
C     when no reduction of the 1-norm of the matrix is possible with
C     DGEBAL, as for upper triangular matrices. LAPACK Library routine
C     DGEBAK, with parameters ILO = 1, IHI = N, and JOB = 'S', should
C     be used to apply the backward transformation.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrix A.  N >= 0.
C
C     MAXRED  (input/output) DOUBLE PRECISION
C             On entry, the maximum allowed reduction in the 1-norm of
C             A (in an iteration) if zero rows or columns are
C             encountered.
C             If MAXRED > 0.0, MAXRED must be larger than one (to enable
C             the norm reduction).
C             If MAXRED <= 0.0, then the value 10.0 for MAXRED is
C             used.
C             On exit, if the 1-norm of the given matrix A is non-zero,
C             the ratio between the 1-norm of the given matrix and the
C             1-norm of the balanced matrix. Usually, this ratio will be
C             larger than one, but it can sometimes be one, or even less
C             than one (for instance, for some companion matrices).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the input matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the balanced matrix.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     SCALE   (output) DOUBLE PRECISION array, dimension (N)
C             The scaling factors applied to A.  If D(j) is the scaling
C             factor applied to row and column j, then SCALE(j) = D(j),
C             for j = 1,...,N.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit.
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     Balancing consists of applying a diagonal similarity
C     transformation inv(D) * A * D to make the 1-norms of each row
C     of A and its corresponding column nearly equal.
C
C     Information about the diagonal matrix D is returned in the vector
C     SCALE.
C
C     REFERENCES
C
C     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routine MB04AD by T.W.C. Williams,
C     Kingston Polytechnic, United Kingdom, October 1984.
C     This subroutine is based on LAPACK routine DGEBAL, and routine
C     BALABC (A. Varga, German Aerospace Research Establishment, DLR).
C
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Balancing, eigenvalue, matrix algebra, matrix operations,
C     similarity transformation.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 1.0D+1 )
      DOUBLE PRECISION   FACTOR, MAXR
      PARAMETER          ( FACTOR = 0.95D+0, MAXR = 10.0D+0 )
C     ..
C     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   MAXRED
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
C     ..
C     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IRA, J
      DOUBLE PRECISION   ANORM, C, CA, F, G, MAXNRM, R, RA, S, SFMAX1,
     $                   SFMAX2, SFMIN1, SFMIN2, SRED
C     ..
C     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE, IDAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL           DSCAL, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
C
C     Test the scalar input arguments.
C
      INFO  = 0
C
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( MAXRED.GT.ZERO .AND. MAXRED.LT.ONE ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB04MD', -INFO )
         RETURN
      END IF
C
      IF( N.EQ.0 )
     $   RETURN
C
      DO 10 I = 1, N
         SCALE( I ) = ONE
   10 CONTINUE
C
C     Compute the 1-norm of matrix A and exit if it is zero.
C
      ANORM = DLANGE( '1-norm', N, N, A, LDA, SCALE )
      IF( ANORM.EQ.ZERO )
     $   RETURN
C
C     Set some machine parameters and the maximum reduction in the
C     1-norm of A if zero rows or columns are encountered.
C
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
C
      SRED = MAXRED
      IF( SRED.LE.ZERO ) SRED = MAXR
C
      MAXNRM = MAX( ANORM/SRED, SFMIN1 )
C
C     Balance the matrix.
C
C     Iterative loop for norm reduction.
C
   20 CONTINUE
      NOCONV = .FALSE.
C
      DO 80 I = 1, N
         C = ZERO
         R = ZERO
C
         DO 30 J = 1, N
            IF( J.EQ.I )
     $         GO TO 30
            C = C + ABS( A( J, I ) )
            R = R + ABS( A( I, J ) )
   30    CONTINUE
         ICA = IDAMAX( N, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IDAMAX( N, A( I, 1 ), LDA )
         RA = ABS( A( I, IRA ) )
C
C        Special case of zero C and/or R.
C
         IF( C.EQ.ZERO .AND. R.EQ.ZERO )
     $      GO TO 80
         IF( C.EQ.ZERO ) THEN
            IF( R.LE.MAXNRM)
     $         GO TO 80
            C = MAXNRM
         END IF
         IF( R.EQ.ZERO ) THEN
            IF( C.LE.MAXNRM )
     $         GO TO 80
            R = MAXNRM
         END IF
C
C        Guard against zero C or R due to underflow.
C
         G = R / SCLFAC
         F = ONE
         S = C + R
   40    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR.
     $       MIN( R, G, RA ).LE.SFMIN2 )GO TO 50
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 40
C
   50    CONTINUE
         G = C / SCLFAC
   60    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR.
     $       MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 70
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 60
C
C        Now balance.
C
   70    CONTINUE
         IF( ( C+R ).GE.FACTOR*S )
     $      GO TO 80
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 )
     $         GO TO 80
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F )
     $         GO TO 80
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
C
         CALL DSCAL( N, G, A( I, 1 ), LDA )
         CALL DSCAL( N, F, A( 1, I ), 1 )
C
   80 CONTINUE
C
      IF( NOCONV )
     $   GO TO 20
C
C     Set the norm reduction parameter.
C
      MAXRED = ANORM/DLANGE( '1-norm', N, N, A, LDA, SCALE )
C
      RETURN
C *** End of MB04MD ***
      END
