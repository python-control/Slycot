      SUBROUTINE TG01BD( JOBE, COMPQ, COMPZ, N, M, P, ILO, IHI, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, DWORK,
     $                   LDWORK, INFO )
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
C     To reduce the matrices A and E of the system pencil
C
C             S =  ( A  B ) - lambda ( E  0 ) ,
C                  ( C  0 )          ( 0  0 )
C
C     corresponding to the descriptor triple (A-lambda E,B,C),
C     to generalized upper Hessenberg form using orthogonal
C     transformations,
C
C          Q' * A * Z = H,   Q' * E * Z = T,
C
C     where H is upper Hessenberg, T is upper triangular, Q and Z
C     are orthogonal, and ' means transpose. The corresponding
C     transformations, written compactly as diag(Q',I) * S * diag(Z,I),
C     are also applied to B and C, getting Q' * B and C * Z.
C
C     The orthogonal matrices Q and Z are determined as products of
C     Givens rotations. They may either be formed explicitly, or they
C     may be postmultiplied into input matrices Q1 and Z1, so that
C
C          Q1 * A * Z1' = (Q1*Q) * H * (Z1*Z)'
C          Q1 * E * Z1' = (Q1*Q) * T * (Z1*Z)'.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBE    CHARACTER*1
C             Specifies whether E is a general square or an upper
C             triangular matrix, as follows:
C             = 'G':  E is a general square matrix;
C             = 'U':  E is an upper triangular matrix.
C
C     COMPQ   CHARACTER*1
C             Indicates what should be done with matrix Q, as follows:
C             = 'N':  do not compute Q;
C             = 'I':  Q is initialized to the unit matrix, and the
C                     orthogonal matrix Q is returned;
C             = 'V':  Q must contain an orthogonal matrix Q1 on entry,
C                     and the product Q1*Q is returned.
C
C     COMPZ   CHARACTER*1
C             Indicates what should be done with matrix Z, as follows:
C             = 'N':  do not compute Z;
C             = 'I':  Z is initialized to the unit matrix, and the
C                     orthogonal matrix Z is returned;
C             = 'V':  Z must contain an orthogonal matrix Z1 on entry,
C                     and the product Z1*Z is returned.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A, E, and the number of rows of
C             the matrix B.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of the matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of the matrix C.  P >= 0.
C
C     ILO     (input) INTEGER
C     IHI     (input) INTEGER
C             It is assumed that A and E are already upper triangular in
C             rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI could
C             normally be set by a previous call to LAPACK Library
C             routine DGGBAL; otherwise they should be set to 1 and N,
C             respectively.
C             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
C             If JOBE = 'U', the matrix E is assumed upper triangular.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the upper Hessenberg matrix H = Q' * A * Z. The elements
C             below the first subdiagonal are set to zero.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the descriptor matrix E. If JOBE = 'U', this
C             matrix is assumed upper triangular.
C             On exit, the leading N-by-N part of this array contains
C             the upper triangular matrix T = Q' * E * Z. The elements
C             below the diagonal are set to zero.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix B.
C             On exit, if M > 0, the leading N-by-M part of this array
C             contains the transformed matrix Q' * B.
C             The array B is not referenced if M = 0.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,N) if M > 0;  LDB >= 1 if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, if P > 0, the leading P-by-N part of this array
C             contains the transformed matrix C * Z.
C             The array C is not referenced if P = 0.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             If COMPQ = 'N':  Q is not referenced;
C             If COMPQ = 'I':  on entry, Q need not be set, and on exit
C                              it contains the orthogonal matrix Q,
C                              where Q' is the product of the Givens
C                              transformations which are applied to A,
C                              E, and B on the left;
C             If COMPQ = 'V':  on entry, Q must contain an orthogonal
C                              matrix Q1, and on exit this is
C                              overwritten by Q1*Q.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= 1,        if COMPQ = 'N';
C             LDQ >= MAX(1,N), if COMPQ = 'I' or 'V'.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If COMPZ = 'N':  Z is not referenced;
C             If COMPZ = 'I':  on entry, Z need not be set, and on exit
C                              it contains the orthogonal matrix Z,
C                              which is the product of the Givens
C                              transformations applied to A, E, and C
C                              on the right;
C             If COMPZ = 'V':  on entry, Z must contain an orthogonal
C                              matrix Z1, and on exit this is
C                              overwritten by Z1*Z.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.
C             LDZ >= 1,        if COMPZ = 'N';
C             LDZ >= MAX(1,N), if COMPZ = 'I' or 'V'.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 1,                          if JOBE = 'U';
C             LDWORK >= MAX(1,IHI+1-ILO+MAX(NI,M)), if JOBE = 'G', where
C             NI = N+1-ILO, if COMPQ = 'N', and NI = N, otherwise.
C             For good performance, if JOBE = 'G', LDWORK must generally
C             be larger, LDWORK >= MAX(1,IHI+1-ILO+MAX(NI,M)*NB), where
C             NB is the optimal block size.
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
C     First, this routine computes the QR factorization of E and applies
C     the transformations to A, B, and possibly Q. Then, the routine
C     reduces A to upper Hessenberg form, preserving E triangular, by
C     an unblocked reduction [1], using two sequences of plane rotations
C     applied alternately from the left and from the right. The
C     corresponding transformations may be accumulated and/or applied
C     to the matrices B and C. If JOBE = 'U', the initial reduction of E
C     to upper triangular form is skipped.
C
C     This routine is a modification and extension of the LAPACK Library
C     routine DGGHRD [2].
C
C     REFERENCES
C
C     [1] Golub, G.H. and van Loan, C.F.
C         Matrix Computations. Third Edition.
C         M. D. Johns Hopkins University Press, Baltimore, 1996.
C
C     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J.,
C         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A.,
C         Ostrouchov, S., and Sorensen, D.
C         LAPACK Users' Guide: Second Edition.
C         SIAM, Philadelphia, 1995.
C
C     CONTRIBUTOR
C
C     D. Sima, University of Bucharest, May 2001.
C     V. Sima, Research Institute for Informatics, Bucharest, May 2001.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Eigenvalue, matrix algebra, matrix operations, similarity
C     transformation.
C
C  *********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOBE
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDC, LDE, LDQ,
     $                   LDWORK, LDZ, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   DWORK(  * ), E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      LOGICAL            ILQ, ILZ, INQ, INZ, UPPER, WITHB, WITHC
      INTEGER            IERR, ITAU, IWRK, JCOL, JROW, MAXWRK, MINWRK
      DOUBLE PRECISION   CS, S, TEMP
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLARTG, DLASET, DORMQR, DROT, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX
C
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      UPPER = LSAME( JOBE,  'U' )
      INQ   = LSAME( COMPQ, 'I' )
      ILQ   = LSAME( COMPQ, 'V' ) .OR. INQ
      INZ   = LSAME( COMPZ, 'I' )
      ILZ   = LSAME( COMPZ, 'V' ) .OR. INZ
      WITHB = M.GT.0
      WITHC = P.GT.0
C
      INFO = 0
      IF( .NOT.( UPPER .OR. LSAME( JOBE, 'G' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ILQ .OR. LSAME( COMPQ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( ILZ .OR. LSAME( COMPZ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -7
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -8
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( ( WITHB .AND. LDB.LT.N ) .OR. LDB.LT.1 ) THEN
         INFO = -14
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -16
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -18
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -20
      ELSE
         JROW = IHI + 1 - ILO
         JCOL = N + 1 - ILO
         IF( UPPER ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            IF( ILQ ) THEN
               MINWRK = N
            ELSE
               MINWRK = JCOL
            END IF
            MINWRK = MAX( 1, JROW + MAX( MINWRK, M ) )
         END IF
         IF( LDWORK.LT.MINWRK )
     $      INFO = -22
      END IF
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01BD', -INFO )
         RETURN
      END IF
C
C     Initialize Q and Z if desired.
C
      IF( INQ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( INZ )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
C
C     Quick return if possible.
C
      IF( N.LE.1 ) THEN
         DWORK( 1 ) = ONE
         RETURN
      END IF
C
      IF( .NOT.UPPER ) THEN
C
C        Reduce E to triangular form (QR decomposition of E).
C
C        (Note: Comments in the code beginning "Workspace:" describe the
C        minimal amount of real workspace needed at that point in the
C        code, as well as the preferred amount for good performance.
C        NB refers to the optimal block size for the immediately
C        following subroutine, as returned by ILAENV.)
C
C        Workspace: need   IHI+1-ILO+N+1-ILO;
C                   prefer IHI+1-ILO+(N+1-ILO)*NB.
C
         ITAU = 1
         IWRK = ITAU + JROW
         CALL DGEQRF( JROW, JCOL, E( ILO, ILO ), LDE, DWORK( ITAU ),
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MINWRK )
C
C        Apply the orthogonal transformation to matrices A, B, and Q.
C        Workspace: need   IHI+1-ILO+N+1-ILO;
C                   prefer IHI+1-ILO+(N+1-ILO)*NB.
C
         CALL DORMQR( 'Left', 'Transpose', JROW, JCOL, JROW,
     $                E( ILO, ILO ), LDE, DWORK( ITAU ), A( ILO, ILO ),
     $                LDA, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
         IF ( WITHB ) THEN
C
C           Workspace: need   IHI+1-ILO+M;
C                      prefer IHI+1-ILO+M*NB.
C
            CALL DORMQR( 'Left', 'Transpose', JROW, M, JROW,
     $                   E( ILO, ILO ), LDE, DWORK( ITAU ), B( ILO, 1 ),
     $                   LDB, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         END IF
C
         IF( ILQ ) THEN
C
C           Workspace: need   IHI+1-ILO+N;
C                      prefer IHI+1-ILO+N*NB.
C
            CALL DORMQR( 'Right', 'No Transpose', N, JROW, JROW,
     $                   E( ILO, ILO ), LDE, DWORK( ITAU ), Q( 1, ILO ),
     $                   LDQ, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         END IF
      END IF
C
C     Zero out lower triangle of E.
C
      IF( JROW.GT.1 )
     $   CALL DLASET( 'Lower', JROW-1, JROW-1, ZERO, ZERO,
     $                E( ILO+1, ILO ), LDE )
C
C     Reduce A and E and apply the transformations to B, C, Q and Z.
C
      DO 20 JCOL = ILO, IHI - 2
C
         DO 10 JROW = IHI, JCOL + 2, -1
C
C           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL).
C
            TEMP = A( JROW-1, JCOL )
            CALL DLARTG( TEMP, A( JROW, JCOL ), CS, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = ZERO
            CALL DROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, CS, S )
            CALL DROT( N+2-JROW, E( JROW-1, JROW-1 ), LDE,
     $                 E( JROW, JROW-1 ), LDE, CS, S )
            IF( WITHB )
     $         CALL DROT( M, B( JROW-1, 1 ), LDB, B( JROW, 1 ), LDB,
     $                    CS, S )
            IF( ILQ )
     $         CALL DROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, CS, S )
C
C           Step 2: rotate columns JROW, JROW-1 to kill E(JROW,JROW-1).
C
            TEMP = E( JROW, JROW )
            CALL DLARTG( TEMP, E( JROW, JROW-1 ), CS, S,
     $                   E( JROW, JROW ) )
            E( JROW, JROW-1 ) = ZERO
            CALL DROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, CS, S )
            CALL DROT( JROW-1, E( 1, JROW ), 1, E( 1, JROW-1 ), 1, CS,
     $                 S )
            IF( WITHC )
     $         CALL DROT( P, C( 1, JROW ), 1, C( 1, JROW-1 ), 1, CS, S )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, CS, S )
   10    CONTINUE
C
   20 CONTINUE
C
      DWORK( 1 ) = MAXWRK
      RETURN
C *** Last line of TG01BD ***
      END
