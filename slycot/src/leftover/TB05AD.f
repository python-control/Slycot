      SUBROUTINE TB05AD( BALEIG, INITA, N, M, P, FREQ, A, LDA, B, LDB,
     $                   C, LDC, RCOND, G, LDG, EVRE, EVIM, HINVB,
     $                   LDHINV, IWORK, DWORK, LDWORK, ZWORK, LZWORK,
     $                   INFO )
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
C     To find the complex frequency response matrix (transfer matrix)
C     G(freq) of the state-space representation (A,B,C) given by
C                                   -1
C        G(freq) = C * ((freq*I - A)  ) * B
C
C     where A, B and C are real N-by-N, N-by-M and P-by-N matrices
C     respectively and freq is a complex scalar.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     BALEIG  CHARACTER*1
C             Determines whether the user wishes to balance matrix A
C             and/or compute its eigenvalues and/or estimate the
C             condition number of the problem as follows:
C             = 'N':  The matrix A should not be balanced and neither
C                     the eigenvalues of A nor the condition number
C                     estimate of the problem are to be calculated;
C             = 'C':  The matrix A should not be balanced and only an
C                     estimate of the condition number of the problem
C                     is to be calculated;
C             = 'B' or 'E' and INITA = 'G':  The matrix A is to be
C                     balanced and its eigenvalues calculated;
C             = 'A' and INITA = 'G':  The matrix A is to be balanced,
C                     and its eigenvalues and an estimate of the
C                     condition number of the problem are to be
C                     calculated.
C
C     INITA   CHARACTER*1
C             Specifies whether or not the matrix A is already in upper
C             Hessenberg form as follows:
C             = 'G':  The matrix A is a general matrix;
C             = 'H':  The matrix A is in upper Hessenberg form and
C                     neither balancing nor the eigenvalues of A are
C                     required.
C             INITA must be set to 'G' for the first call to the
C             routine, unless the matrix A is already in upper
C             Hessenberg form and neither balancing nor the eigenvalues
C             of A are required. Thereafter, it must be set to 'H' for
C             all subsequent calls.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of states, i.e. the order of the state
C             transition matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of inputs, i.e. the number of columns in the
C             matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of outputs, i.e. the number of rows in the
C             matrix C.  P >= 0.
C
C     FREQ    (input) COMPLEX*16
C             The frequency freq at which the frequency response matrix
C             (transfer matrix) is to be evaluated.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state transition matrix A.
C             If INITA = 'G', then, on exit, the leading N-by-N part of
C             this array contains an upper Hessenberg matrix similar to
C             (via an orthogonal matrix consisting of a sequence of
C             Householder transformations) the original state transition
C             matrix A.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the input/state matrix B.
C             If INITA = 'G', then, on exit, the leading N-by-M part of
C             this array contains the product of the transpose of the
C             orthogonal transformation matrix used to reduce A to upper
C             Hessenberg form and the original input/state matrix B.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             If INITA = 'G', then, on exit, the leading P-by-N part of
C             this array contains the product of the original output/
C             state matrix C and the orthogonal transformation matrix
C             used to reduce A to upper Hessenberg form.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     RCOND   (output) DOUBLE PRECISION
C             If BALEIG = 'C' or BALEIG = 'A', then RCOND contains an
C             estimate of the reciprocal of the condition number of
C             matrix H with respect to inversion (see METHOD).
C
C     G       (output) COMPLEX*16 array, dimension (LDG,M)
C             The leading P-by-M part of this array contains the
C             frequency response matrix G(freq).
C
C     LDG     INTEGER
C             The leading dimension of array G.  LDG >= MAX(1,P).
C
C     EVRE,   (output) DOUBLE PRECISION arrays, dimension (N)
C     EVIM    If INITA = 'G' and BALEIG = 'B' or 'E' or BALEIG = 'A',
C             then these arrays contain the real and imaginary parts,
C             respectively, of the eigenvalues of the matrix A.
C             Otherwise, these arrays are not referenced.
C
C     HINVB   (output) COMPLEX*16 array, dimension (LDHINV,M)
C             The leading N-by-M part of this array contains the
C                      -1
C             product H  B.
C
C     LDHINV  INTEGER
C             The leading dimension of array HINVB.  LDHINV >= MAX(1,N).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (N)
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1, N - 1 + MAX(N,M,P)),
C                       if INITA = 'G' and BALEIG = 'N', or 'B', or 'E';
C             LDWORK >= MAX(1, N + MAX(N,M-1,P-1)),
C                       if INITA = 'G' and BALEIG = 'C', or 'A';
C             LDWORK >= MAX(1, 2*N),
C                       if INITA = 'H' and BALEIG = 'C', or 'A';
C             LDWORK >= 1, otherwise.
C             For optimum performance when INITA = 'G' LDWORK should be
C             larger.
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C
C     LZWORK  INTEGER
C             The length of the array ZWORK.
C             LZWORK >= MAX(1,N*N+2*N), if BALEIG = 'C', or 'A';
C             LZWORK >= MAX(1,N*N),     otherwise.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if more than 30*N iterations are required to
C                   isolate all the eigenvalues of the matrix A; the
C                   computations are continued;
C             = 2:  if either FREQ is too near to an eigenvalue of the
C                   matrix A, or RCOND is less than EPS, where EPS is
C                   the machine  precision (see LAPACK Library routine
C                   DLAMCH).
C
C     METHOD
C
C     The matrix A is first balanced (if BALEIG = 'B' or 'E', or
C     BALEIG = 'A') and then reduced to upper Hessenberg form; the same
C     transformations are applied to the matrix B and the matrix C.
C     The complex Hessenberg matrix  H = (freq*I - A) is then used
C                       -1
C     to solve for C * H  * B.
C
C     Depending on the input values of parameters BALEIG and INITA,
C     the eigenvalues of matrix A and the condition number of
C     matrix H with respect to inversion are also calculated.
C
C     REFERENCES
C
C     [1] Laub, A.J.
C         Efficient Calculation of Frequency Response Matrices from
C         State-Space Models.
C         ACM TOMS, 12, pp. 26-33, 1986.
C
C     NUMERICAL ASPECTS
C                               3
C     The algorithm requires 0(N ) operations.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996.
C     Supersedes Release 2.0 routine TB01FD by A.J.Laub, University of
C     Southern California, Los Angeles, CA 90089, United States of
C     America, June 1982.
C
C     REVISIONS
C
C     V. Sima, February 22, 1998 (changed the name of TB01RD).
C     V. Sima, February 12, 1999, August 7, 2003.
C     A. Markovski, Technical University of Sofia, September 30, 2003.
C     V. Sima, October 1, 2003.
C
C     KEYWORDS
C
C     Frequency response, Hessenberg form, matrix algebra, input output
C     description, multivariable system, orthogonal transformation,
C     similarity transformation, state-space representation, transfer
C     matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16        CZERO
      PARAMETER         ( CZERO = ( 0.0D0, 0.0D0 ) )
C     .. Scalar Arguments ..
      CHARACTER         BALEIG, INITA
      INTEGER           INFO, LDA, LDB, LDC, LDG, LDHINV, LDWORK,
     $                  LZWORK, M, N, P
      DOUBLE PRECISION  RCOND
      COMPLEX*16        FREQ
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), EVIM(*),
     $                  EVRE(*)
      COMPLEX*16        ZWORK(*), G(LDG,*), HINVB(LDHINV,*)
C     .. Local Scalars ..
      CHARACTER         BALANC
      LOGICAL           LBALBA, LBALEA, LBALEB, LBALEC, LINITA
      INTEGER           I, IGH, IJ, ITAU, J, JJ, JP, JWORK, K, LOW,
     $                  WRKOPT
      DOUBLE PRECISION  HNORM, T
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DASUM, DLAMCH
      EXTERNAL          DASUM, DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          DGEBAL, DGEHRD, DHSEQR, DORMHR, DSCAL, DSWAP,
     $                  MB02RZ, MB02SZ, MB02TZ, XERBLA, ZLASET
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LBALEC = LSAME( BALEIG, 'C' )
      LBALEB = LSAME( BALEIG, 'B' ) .OR. LSAME( BALEIG, 'E' )
      LBALEA = LSAME( BALEIG, 'A' )
      LBALBA = LBALEB.OR.LBALEA
      LINITA = LSAME( INITA,  'G' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LBALEC .AND. .NOT.LBALBA .AND.
     $    .NOT.LSAME( BALEIG, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LINITA .AND. .NOT.LSAME( INITA, 'H' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDG.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDHINV.LT.MAX( 1, N ) ) THEN
         INFO = -19
      ELSE IF( ( LINITA .AND. .NOT.LBALEC .AND. .NOT.LBALEA .AND.
     $           LDWORK.LT.N - 1 + MAX( N, M, P ) ) .OR.
     $         ( LINITA .AND. ( LBALEC .OR. LBALEA ) .AND.
     $           LDWORK.LT.N + MAX( N, M-1, P-1 ) ) .OR.
     $         ( .NOT.LINITA .AND. ( LBALEC .OR. LBALEA ) .AND.
     $           LDWORK.LT.2*N ) .OR. ( LDWORK.LT.1 ) ) THEN
         INFO = -22
      ELSE IF( ( ( LBALEC .OR. LBALEA ) .AND. LZWORK.LT.N*( N + 2 ) )
     $      .OR. ( LZWORK.LT.MAX( 1, N*N ) ) ) THEN
         INFO = -24
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return
C
         CALL XERBLA( 'TB05AD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         IF ( MIN( M, P ).GT.0 )
     $      CALL ZLASET( 'Full', P, M, CZERO, CZERO, G, LDG )
         RCOND = ONE
         DWORK(1) = ONE
         RETURN
      END IF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      WRKOPT = 1
C
      IF ( LINITA ) THEN
         BALANC = 'N'
         IF ( LBALBA ) BALANC = 'B'
C
C        Workspace: need N.
C
         CALL DGEBAL( BALANC, N, A, LDA, LOW, IGH, DWORK, INFO )
         IF ( LBALBA ) THEN
C
C           Adjust B and C matrices based on information in the
C           vector DWORK which describes the balancing of A and is
C           defined in the subroutine DGEBAL.
C
            DO 10 J = 1, N
               JJ = J
               IF ( JJ.LT.LOW .OR. JJ.GT.IGH ) THEN
                  IF ( JJ.LT.LOW ) JJ = LOW - JJ
                  JP = DWORK(JJ)
                  IF ( JP.NE.JJ ) THEN
C
C                    Permute rows of B.
C
                     IF ( M.GT.0 )
     $                  CALL DSWAP( M, B(JJ,1), LDB, B(JP,1), LDB )
C
C                    Permute columns of C.
C
                     IF ( P.GT.0 )
     $                  CALL DSWAP( P, C(1,JJ), 1, C(1,JP), 1 )
                  END IF
               END IF
   10       CONTINUE
C
            IF ( IGH.NE.LOW ) THEN
C
               DO 20 J = LOW, IGH
                  T = DWORK(J)
C
C                 Scale rows of permuted B.
C
                  IF ( M.GT.0 )
     $               CALL DSCAL( M, ONE/T, B(J,1), LDB )
C
C                 Scale columns of permuted C.
C
                  IF ( P.GT.0 )
     $               CALL DSCAL( P, T, C(1,J), 1 )
   20          CONTINUE
C
            END IF
         END IF
C
C        Reduce A to Hessenberg form by orthogonal similarities and
C        accumulate the orthogonal transformations into B and C.
C        Workspace: need 2*N - 1;  prefer N - 1 + N*NB.
C
         ITAU = 1
         JWORK = ITAU + N - 1
         CALL DGEHRD( N, LOW, IGH, A, LDA, DWORK(ITAU), DWORK(JWORK),
     $                LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Workspace: need N - 1 + M;  prefer N - 1 + M*NB.
C
         CALL DORMHR( 'Left', 'Transpose', N, M, LOW, IGH, A, LDA,
     $                DWORK(ITAU), B, LDB, DWORK(JWORK), LDWORK-JWORK+1,
     $                INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
C
C        Workspace: need N - 1 + P;  prefer N - 1 + P*NB.
C
         CALL DORMHR( 'Right', 'No transpose', P, N, LOW, IGH, A, LDA,
     $                DWORK(ITAU), C, LDC, DWORK(JWORK), LDWORK-JWORK+1,
     $                INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) )+JWORK-1 )
         IF ( LBALBA ) THEN
C
C           Temporarily store Hessenberg form of A in array ZWORK.
C
            IJ = 0
            DO 40 J = 1, N
C
               DO 30 I = 1, N
                  IJ = IJ + 1
                  ZWORK(IJ) = DCMPLX( A(I,J), ZERO )
   30          CONTINUE
C
   40       CONTINUE
C
C           Compute the eigenvalues of A if that option is requested.
C           Workspace: need N.
C
            CALL DHSEQR( 'Eigenvalues', 'No Schur', N, LOW, IGH, A, LDA,
     $                   EVRE, EVIM, DWORK, 1, DWORK, LDWORK, INFO )
C
C           Restore upper Hessenberg form of A.
C
            IJ = 0
            DO 60 J = 1, N
C
               DO 50 I = 1, N
                  IJ = IJ + 1
                  A(I,J) = DBLE( ZWORK(IJ) )
   50          CONTINUE
C
   60       CONTINUE
C
            IF ( INFO.GT.0 ) THEN
C
C              DHSEQR could not evaluate the eigenvalues of A.
C
               INFO = 1
            END IF
         END IF
      END IF
C
C     Update  H := (FREQ * I) - A   with appropriate value of FREQ.
C
      IJ = 0
      JJ = 1
      DO 80 J = 1, N
C
         DO 70 I = 1, N
            IJ = IJ + 1
            ZWORK(IJ) = -DCMPLX( A(I,J), ZERO )
   70    CONTINUE
C
         ZWORK(JJ) = FREQ + ZWORK(JJ)
         JJ = JJ + N + 1
   80 CONTINUE
C
      IF ( LBALEC .OR. LBALEA ) THEN
C
C        Efficiently compute the 1-norm of the matrix for condition
C        estimation.
C
         HNORM = ZERO
         JJ = 1
C
         DO 90 J = 1, N
            T = ABS( ZWORK(JJ) ) + DASUM( J-1, A(1,J), 1 )
            IF ( J.LT.N ) T = T + ABS( A(J+1,J) )
            HNORM = MAX( HNORM, T )
            JJ = JJ + N + 1
   90    CONTINUE
C
      END IF
C
C     Factor the complex Hessenberg matrix.
C
      CALL MB02SZ( N, ZWORK, N, IWORK, INFO )
      IF ( INFO.NE.0 ) INFO = 2
C
      IF ( LBALEC .OR. LBALEA ) THEN
C
C        Estimate the condition of the matrix.
C
C        Workspace: need 2*N.
C
         CALL MB02TZ( '1-norm', N, HNORM, ZWORK, N, IWORK, RCOND, DWORK,
     $                ZWORK(N*N+1), INFO )
         WRKOPT = MAX( WRKOPT, 2*N )
         IF ( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = 2
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return: Linear system is numerically or exactly singular.
C
         RETURN
      END IF
C
C     Compute  (H-INVERSE)*B.
C
      DO 110 J = 1, M
C
         DO 100 I = 1, N
            HINVB(I,J) = DCMPLX( B(I,J), ZERO )
  100    CONTINUE
C
  110 CONTINUE
C
      CALL MB02RZ( 'No transpose', N, M, ZWORK, N, IWORK, HINVB, LDHINV,
     $             INFO )
C
C     Compute  C*(H-INVERSE)*B.
C
      DO 150 J = 1, M
C
         DO 120 I = 1, P
            G(I,J) = CZERO
  120    CONTINUE
C
         DO 140 K = 1, N
C
            DO 130 I = 1, P
               G(I,J) = G(I,J) + DCMPLX( C(I,K), ZERO )*HINVB(K,J)
  130       CONTINUE
C
  140    CONTINUE
C
  150 CONTINUE
C
C     G now contains the desired frequency response matrix.
C     Set the optimal workspace.
C
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of TB05AD ***
      END
