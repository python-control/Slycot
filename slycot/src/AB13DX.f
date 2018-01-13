      DOUBLE PRECISION FUNCTION AB13DX( DICO, JOBE, JOBD, N, M, P,
     $                                  OMEGA, A, LDA, E, LDE, B, LDB,
     $                                  C, LDC, D, LDD, IWORK, DWORK,
     $                                  LDWORK, CWORK, LCWORK, INFO )
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
C     To compute the maximum singular value of a given continuous-time
C     or discrete-time transfer-function matrix, either standard or in
C     the descriptor form,
C
C                                     -1
C        G(lambda) = C*( lambda*E - A ) *B + D ,
C
C     for a given complex value lambda, where lambda = j*omega, in the
C     continuous-time case, and lambda = exp(j*omega), in the
C     discrete-time case. The matrices A, E, B, C, and D are real
C     matrices of appropriate dimensions. Matrix A must be in an upper
C     Hessenberg form, and if JOBE ='G', the matrix E must be upper
C     triangular. The matrices B and C must correspond to the system
C     in (generalized) Hessenberg form.
C
C     FUNCTION VALUE
C
C     AB13DX   DOUBLE PRECISION
C              The maximum singular value of G(lambda).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system, as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOBE    CHARACTER*1
C             Specifies whether E is an upper triangular or an identity
C             matrix, as follows:
C             = 'G':  E is a general upper triangular matrix;
C             = 'I':  E is the identity matrix.
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears in
C             the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The row size of the matrix C.  P >= 0.
C
C     OMEGA   (input) DOUBLE PRECISION
C             The frequency value for which the calculations should be
C             done.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N upper Hessenberg part of this
C             array must contain the state dynamics matrix A in upper
C             Hessenberg form. The elements below the subdiagonal are
C             not referenced.
C             On exit, if M > 0, P > 0, OMEGA = 0, DICO = 'C', B <> 0,
C             and C <> 0, the leading N-by-N upper Hessenberg part of
C             this array contains the factors L and U from the LU
C             factorization of A (A = P*L*U); the unit diagonal elements
C             of L are not stored, L is lower bidiagonal, and P is
C             stored in IWORK (see SLICOT Library routine MB02SD).
C             Otherwise, this array is unchanged on exit.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             If JOBE = 'G', the leading N-by-N upper triangular part of
C             this array must contain the upper triangular descriptor
C             matrix E of the system. The elements of the strict lower
C             triangular part of this array are not referenced.
C             If JOBE = 'I', then E is assumed to be the identity
C             matrix and is not referenced.
C
C     LDE     INTEGER
C             The leading dimension of the array E.
C             LDE >= MAX(1,N), if JOBE = 'G';
C             LDE >= 1,        if JOBE = 'I'.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the system input matrix B.
C             On exit, if M > 0, P > 0, OMEGA = 0, DICO = 'C', B <> 0,
C             C <> 0, and INFO = 0 or N+1, the leading N-by-M part of
C             this array contains the solution of the system A*X = B.
C             Otherwise, this array is unchanged on exit.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, if JOBD = 'D', the leading P-by-M part of this
C             array must contain the direct transmission matrix D.
C             On exit, if (N = 0, or B = 0, or C = 0) and JOBD = 'D',
C             or (OMEGA = 0, DICO = 'C', JOBD = 'D', and INFO = 0 or
C             N+1), the contents of this array is destroyed.
C             Otherwise, this array is unchanged on exit.
C             This array is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P), if JOBD = 'D';
C             LDD >= 1,        if JOBD = 'Z'.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK), where
C             LIWORK = N, if N > 0, M > 0, P > 0, B <> 0, and C <> 0;
C             LIWORK = 0, otherwise.
C             This array contains the pivot indices in the LU
C             factorization of the matrix lambda*E - A; for 1 <= i <= N,
C             row i of the matrix was interchanged with row IWORK(i).
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK, and DWORK(2), ..., DWORK(MIN(P,M)) contain the
C             singular values of G(lambda), except for the first one,
C             which is returned in the function value AB13DX.
C             If (N = 0, or B = 0, or C = 0) and JOBD = 'Z', the last
C             MIN(P,M)-1 zero singular values of G(lambda) are not
C             stored in DWORK(2), ..., DWORK(MIN(P,M)).
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= MAX(1, LDW1 + LDW2 ),
C             LDW1 = P*M, if N > 0, B <> 0, C <> 0, OMEGA = 0,
C                            DICO = 'C', and JOBD = 'Z';
C             LDW1 = 0,   otherwise;
C             LDW2 = MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M), 5*MIN(P,M)),
C                         if (N = 0, or B = 0, or C = 0) and JOBD = 'D',
C                         or (N > 0, B <> 0, C <> 0, OMEGA = 0, and
C                             DICO = 'C');
C             LDW2 = 0,   if (N = 0, or B = 0, or C = 0) and JOBD = 'Z',
C                         or MIN(P,M) = 0;
C             LDW2 = 6*MIN(P,M), otherwise.
C             For good performance, LDWORK must generally be larger.
C
C     CWORK   COMPLEX*16 array, dimension (LCWORK)
C             On exit, if INFO = 0, CWORK(1) contains the optimal
C             LCWORK.
C
C     LCWORK  INTEGER
C             The dimension of the array CWORK.
C             LCWORK >= 1, if N = 0, or B = 0, or C = 0, or (OMEGA = 0
C                             and DICO = 'C') or MIN(P,M) = 0;
C             LCWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)),
C                          otherwise.
C             For good performance, LCWORK must generally be larger.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = i, U(i,i) is exactly zero; the LU
C                   factorization of the matrix lambda*E - A has been
C                   completed, but the factor U is exactly singular,
C                   i.e., the matrix lambda*E - A is exactly singular;
C             = N+1:  the SVD algorithm for computing singular values
C                   did not converge.
C
C     METHOD
C
C     The routine implements standard linear algebra calculations,
C     taking problem structure into account. LAPACK Library routines
C     DGESVD and ZGESVD are used for finding the singular values.
C
C     CONTRIBUTORS
C
C     D. Sima, University of Bucharest, May 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Sep. 2005.
C
C     KEYWORDS
C
C     H-infinity optimal control, robust control, system norm.
C
C     ******************************************************************
C
C     .. Parameters ..
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          DICO, JOBD, JOBE
      INTEGER            INFO, LCWORK, LDA, LDB, LDC, LDD, LDE, LDWORK,
     $                   M, N, P
      DOUBLE PRECISION   OMEGA
C     ..
C     .. Array Arguments ..
      COMPLEX*16         CWORK(  * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), DWORK(  * ), E( LDE, * )
      INTEGER            IWORK(  * )
C     ..
C     .. Local Scalars ..
      LOGICAL            DISCR, FULLE, NODYN, SPECL, WITHD
      INTEGER            I, ICB, ICC, ICD, ICWK, ID, IERR, IS, IWRK, J,
     $                   MAXWRK, MINCWR, MINPM, MINWRK
      DOUBLE PRECISION   BNORM, CNORM, LAMBDI, LAMBDR, UPD
C
C     .. External Functions ..
      DOUBLE PRECISION   DLANGE
      LOGICAL            LSAME
      EXTERNAL           DLANGE, LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DGEMM, DGESVD, MB02RD, MB02RZ, MB02SD, MB02SZ,
     $                   XERBLA, ZGEMM, ZGESVD, ZLACP2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          COS, DCMPLX, INT, MAX, MIN, SIN
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      INFO  = 0
      DISCR = LSAME( DICO, 'D' )
      FULLE = LSAME( JOBE, 'G' )
      WITHD = LSAME( JOBD, 'D' )
C
      IF( .NOT. ( DISCR .OR. LSAME( DICO, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( FULLE .OR. LSAME( JOBE, 'I' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( WITHD .OR. LSAME( JOBD, 'Z' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDE.LT.1 .OR. ( FULLE .AND. LDE.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -17
      ELSE
         BNORM = DLANGE( '1-norm', N, M, B, LDB, DWORK )
         CNORM = DLANGE( '1-norm', P, N, C, LDC, DWORK )
         NODYN = N.EQ.0 .OR. MIN( BNORM, CNORM ).EQ.ZERO
         SPECL = .NOT.NODYN .AND. OMEGA.EQ.ZERO .AND. .NOT.DISCR
         MINPM = MIN( P, M )
C
C        Compute workspace.
C
         IF( MINPM.EQ.0 ) THEN
            MINWRK = 0
         ELSE IF( SPECL .OR. ( NODYN .AND. WITHD ) ) THEN
            MINWRK = MINPM + MAX( 3*MINPM + MAX( P, M ), 5*MINPM )
            IF ( SPECL .AND. .NOT.WITHD )
     $         MINWRK = MINWRK + P*M
         ELSE IF ( NODYN .AND. .NOT.WITHD ) THEN
            MINWRK = 0
         ELSE
            MINWRK = 6*MINPM
         END IF
         MINWRK = MAX( 1, MINWRK )
C
         IF( LDWORK.LT.MINWRK ) THEN
            INFO = -20
         ELSE
            IF ( NODYN .OR. ( OMEGA.EQ.ZERO .AND. .NOT.DISCR ) .OR.
     $           MINPM.EQ.0 ) THEN
               MINCWR = 1
            ELSE
               MINCWR = MAX( 1, ( N + M )*( N + P ) +
     $                          2*MINPM + MAX( P, M ) )
            END IF
            IF( LCWORK.LT.MINCWR )
     $         INFO = -22
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
         AB13DX = ZERO
         CALL XERBLA( 'AB13DX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MINPM.EQ.0 ) THEN
         AB13DX = ZERO
C
         DWORK( 1 ) = ONE
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      IS   = 1
      IWRK = IS + MINPM
C
      IF( NODYN ) THEN
C
C        No dynamics: Determine the maximum singular value of G = D .
C
         IF ( WITHD ) THEN
C
C           Workspace: need   MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M),
C                                            5*MIN(P,M));
C                      prefer larger.
C
            CALL DGESVD( 'No Vectors', 'No Vectors', P, M, D, LDD,
     $                   DWORK( IS ), DWORK, P, DWORK, M, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = N + 1
               AB13DX = ZERO
               RETURN
            END IF
            AB13DX = DWORK( IS )
            MAXWRK = INT( DWORK( IWRK ) ) + IWRK - 1
         ELSE
            AB13DX = ZERO
            MAXWRK = 1
         END IF
C
         DWORK( 1 ) = MAXWRK
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine the maximum singular value of
C        G(lambda) = C*inv(lambda*E - A)*B + D.
C     The (generalized) Hessenberg form of the system is used.
C
      IF ( SPECL ) THEN
C
C        Special continuous-time case:
C        Determine the maximum singular value of the real matrix G(0).
C        Workspace: need   MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M),
C                                         5*MIN(P,M));
C                   prefer larger.
C
         CALL MB02SD( N, A, LDA, IWORK, IERR )
         IF( IERR.GT.0 ) THEN
            INFO = IERR
            DWORK( 1 ) = ONE
            CWORK( 1 ) = ONE
            AB13DX = ZERO
            RETURN
         END IF
         CALL MB02RD( 'No Transpose', N, M, A, LDA, IWORK, B, LDB,
     $                IERR )
         IF ( WITHD ) THEN
            CALL DGEMM(  'No Transpose', 'No Transpose', P, M, N, -ONE,
     $                   C, LDC, B, LDB, ONE, D, LDD )
            CALL DGESVD( 'No Vectors', 'No Vectors', P, M, D, LDD,
     $                   DWORK( IS ), DWORK, P, DWORK, M, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, IERR )
         ELSE
C
C           Additional workspace: need   P*M.
C
            ID   = IWRK
            IWRK = ID + P*M
            CALL DGEMM(  'No Transpose', 'No Transpose', P, M, N, -ONE,
     $                   C, LDC, B, LDB, ZERO, DWORK( ID ), P )
            CALL DGESVD( 'No Vectors', 'No Vectors', P, M, DWORK( ID ),
     $                   P, DWORK( IS ), DWORK, P, DWORK, M,
     $                   DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         END IF
         IF( IERR.GT.0 ) THEN
            INFO = N + 1
            AB13DX = ZERO
            RETURN
         END IF
C
         AB13DX = DWORK( IS )
         DWORK( 1 ) = INT( DWORK( IWRK ) ) + IWRK - 1
         CWORK( 1 ) = ONE
         RETURN
      END IF
C
C     General case: Determine the maximum singular value of G(lambda).
C     Complex workspace:  need   N*N + N*M + P*N + P*M.
C
      ICB  = 1   + N*N
      ICC  = ICB + N*M
      ICD  = ICC + P*N
      ICWK = ICD + P*M
C
      IF ( WITHD ) THEN
         UPD = ONE
      ELSE
         UPD = ZERO
      END IF
C
      IF ( DISCR ) THEN
         LAMBDR = COS( OMEGA )
         LAMBDI = SIN( OMEGA )
C
C        Build lambda*E - A .
C
         IF ( FULLE ) THEN
C
            DO 20 J = 1, N
C
               DO 10 I = 1, J
                  CWORK( I+(J-1)*N ) =
     $               DCMPLX( LAMBDR*E( I, J ) - A( I, J ),
     $                       LAMBDI*E( I, J ) )
   10          CONTINUE
C
               IF( J.LT.N )
     $            CWORK( J+1+(J-1)*N ) = DCMPLX( -A( J+1, J ), ZERO )
   20      CONTINUE
C
         ELSE
C
            DO 40 J = 1, N
C
               DO 30 I = 1, MIN( J+1, N )
                  CWORK( I+(J-1)*N ) = -A( I, J )
   30          CONTINUE
C
               CWORK( J+(J-1)*N ) = DCMPLX( LAMBDR - A( J, J ), LAMBDI )
   40      CONTINUE
C
         END IF
C
      ELSE
C
C        Build j*omega*E - A.
C
         IF ( FULLE ) THEN
C
            DO 60 J = 1, N
C
               DO 50 I = 1, J
                  CWORK( I+(J-1)*N ) =
     $               DCMPLX( -A( I, J ), OMEGA*E( I, J ) )
   50         CONTINUE
C
               IF( J.LT.N )
     $            CWORK( J+1+(J-1)*N ) = DCMPLX( -A( J+1, J ), ZERO )
   60      CONTINUE
C
         ELSE
C
            DO 80 J = 1, N
C
               DO 70 I = 1, MIN( J+1, N )
                  CWORK( I+(J-1)*N ) = -A( I, J )
   70          CONTINUE
C
               CWORK( J+(J-1)*N ) = DCMPLX( -A( J, J ), OMEGA )
   80      CONTINUE
C
         END IF
C
      END IF
C
C     Build G(lambda) .
C
      CALL ZLACP2( 'Full', N, M, B, LDB, CWORK( ICB ), N )
      CALL ZLACP2( 'Full', P, N, C, LDC, CWORK( ICC ), P )
      IF ( WITHD )
     $   CALL ZLACP2( 'Full', P, M, D, LDD, CWORK( ICD ), P )
C
      CALL MB02SZ( N, CWORK, N, IWORK, IERR )
      IF( IERR.GT.0 ) THEN
         INFO = IERR
         DWORK( 1 ) = ONE
         CWORK( 1 ) = ICWK - 1
         AB13DX = ZERO
         RETURN
      END IF
      CALL MB02RZ( 'No Transpose', N, M, CWORK, N, IWORK,
     $             CWORK( ICB ), N, IERR )
      CALL ZGEMM(  'No Transpose', 'No Transpose', P, M, N, CONE,
     $             CWORK( ICC ), P, CWORK( ICB ), N,
     $             DCMPLX( UPD, ZERO ), CWORK( ICD ), P )
C
C     Additional workspace, complex: need   2*MIN(P,M) + MAX(P,M);
C                                    prefer larger;
C                           real:    need   5*MIN(P,M).
C
      CALL ZGESVD( 'No Vectors', 'No Vectors', P, M, CWORK( ICD ), P,
     $             DWORK( IS ), CWORK, P, CWORK, M, CWORK( ICWK ),
     $             LCWORK-ICWK+1, DWORK( IWRK ), IERR )
      IF( IERR.GT.0 ) THEN
         INFO = N + 1
         RETURN
      END IF
      AB13DX = DWORK( IS )
C
      DWORK( 1 ) = 6*MINPM
      CWORK( 1 ) = INT( CWORK( ICWK ) ) + ICWK - 1
C
      RETURN
C *** Last line of AB13DX ***
      END
