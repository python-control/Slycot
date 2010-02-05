      SUBROUTINE MB03YT( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL,
     $                   CSR, SNR )
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
C     To compute the periodic Schur factorization of a real 2-by-2
C     matrix pair (A,B) where B is upper triangular. This routine
C     computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
C     SNR such that
C
C     1) if the pair (A,B) has two real eigenvalues, then
C
C        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
C        [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
C
C        [ b11 b12 ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ]
C        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ],
C
C     2) if the pair (A,B) has a pair of complex conjugate eigenvalues,
C        then
C
C        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
C        [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
C
C        [ b11  0  ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ]
C        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ].
C
C     This is a modified version of the LAPACK routine DLAGV2 for
C     computing the real, generalized Schur decomposition of a
C     two-by-two matrix pencil.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,2)
C             On entry, the leading 2-by-2 part of this array must
C             contain the matrix A.
C             On exit, the leading 2-by-2 part of this array contains
C             the matrix A of the pair in periodic Schur form.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= 2.
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,2)
C             On entry, the leading 2-by-2 part of this array must
C             contain the upper triangular matrix B.
C             On exit, the leading 2-by-2 part of this array contains
C             the matrix B of the pair in periodic Schur form.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= 2.
C
C     ALPHAR  (output) DOUBLE PRECISION array, dimension (2)
C     ALPHAI  (output) DOUBLE PRECISION array, dimension (2)
C     BETA    (output) DOUBLE PRECISION array, dimension (2)
C             (ALPHAR(k)+i*ALPHAI(k))*BETA(k) are the eigenvalues of the
C             pair (A,B), k=1,2, i = sqrt(-1). ALPHAI(1) >= 0.
C
C     CSL     (output) DOUBLE PRECISION
C             The cosine of the first rotation matrix.
C
C     SNL     (output) DOUBLE PRECISION
C             The sine of the first rotation matrix.
C
C     CSR     (output) DOUBLE PRECISION
C             The cosine of the second rotation matrix.
C
C     SNR     (output) DOUBLE PRECISION
C             The sine of the second rotation matrix.
C
C     REFERENCES
C
C     [1] Van Loan, C.
C         Generalized Singular Values with Algorithms and Applications.
C         Ph. D. Thesis, University of Michigan, 1973.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPV2).
C     V. Sima, July 2008, May 2009.
C
C     KEYWORDS
C
C     Eigenvalue, periodic Schur form
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           LDA, LDB
      DOUBLE PRECISION  CSL, CSR, SNL, SNR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ALPHAI(2), ALPHAR(2), B(LDB,*),
     $                  BETA(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORM, BNORM, H1, H2, H3, QQ, R, RR, SAFMIN,
     $                  SCALE1, SCALE2, T, ULP, WI, WR1, WR2
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH, DLAPY2
      EXTERNAL          DLAMCH, DLAPY2
C     .. External Subroutines ..
      EXTERNAL          DLAG2, DLARTG, DLASV2, DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C
C     .. Executable Statements ..
C
      SAFMIN = DLAMCH( 'S' )
      ULP    = DLAMCH( 'P' )
C
C     Scale A.
C
      ANORM = MAX( ABS( A(1,1) ) + ABS( A(2,1) ),
     $             ABS( A(1,2) ) + ABS( A(2,2) ), SAFMIN )
      A(1,1) = A(1,1) / ANORM
      A(1,2) = A(1,2) / ANORM
      A(2,1) = A(2,1) / ANORM
      A(2,2) = A(2,2) / ANORM
C
C     Scale B.
C
      BNORM = MAX( ABS( B(1,1) ), ABS( B(1,2) ) + ABS( B(2,2) ), SAFMIN)
      B(1,1) = B(1,1) / BNORM
      B(1,2) = B(1,2) / BNORM
      B(2,2) = B(2,2) / BNORM
C
C     Check if A can be deflated.
C
      IF ( ABS( A(2,1) ).LE.ULP ) THEN
         CSL = ONE
         SNL = ZERO
         CSR = ONE
         SNR = ZERO
         WI  = ZERO
         A(2,1) = ZERO
         B(2,1) = ZERO
C
C     Check if B is singular.
C
      ELSE IF ( ABS( B(1,1) ).LE.ULP ) THEN
         CALL DLARTG( A(2,2), A(2,1), CSR, SNR, T )
         SNR = -SNR
         CALL DROT( 2, A(1,1), 1, A(1,2), 1, CSR, SNR )
         CALL DROT( 2, B(1,1), LDB, B(2,1), LDB, CSR, SNR )
         CSL = ONE
         SNL = ZERO
         WI  = ZERO
         A(2,1) = ZERO
         B(1,1) = ZERO
         B(2,1) = ZERO
      ELSE IF( ABS( B(2,2) ).LE.ULP ) THEN
         CALL DLARTG( A(1,1), A(2,1), CSL, SNL, R )
         CSR = ONE
         SNR = ZERO
         WI  = ZERO
         CALL DROT( 2, A(1,1), LDA, A(2,1), LDA, CSL, SNL )
         CALL DROT( 2, B(1,1), 1, B(1,2), 1, CSL, SNL )
         A(2,1) = ZERO
         B(2,1) = ZERO
         B(2,2) = ZERO
      ELSE
C
C        B is nonsingular, first compute the eigenvalues of A / adj(B).
C
         R = B(1,1)
         B(1,1) = B(2,2)
         B(2,2) = R
         B(1,2) = -B(1,2)
         CALL DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2,
     $               WI )
C
         IF( WI.EQ.ZERO ) THEN
C
C           Two real eigenvalues, compute s*A-w*B.
C
            H1 = SCALE1*A(1,1) - WR1*B(1,1)
            H2 = SCALE1*A(1,2) - WR1*B(1,2)
            H3 = SCALE1*A(2,2) - WR1*B(2,2)
C
            RR = DLAPY2( H1, H2 )
            QQ = DLAPY2( SCALE1*A(2,1), H3 )
C
            IF ( RR.GT.QQ ) THEN
C
C              Find right rotation matrix to zero 1,1 element of
C              (sA - wB).
C
               CALL DLARTG( H2, H1, CSR, SNR, T )
C
            ELSE
C
C              Find right rotation matrix to zero 2,1 element of
C              (sA - wB).
C
               CALL DLARTG( H3, SCALE1*A(2,1), CSR, SNR, T )
C
            END IF
C
            SNR = -SNR
            CALL DROT( 2, A(1,1), 1, A(1,2), 1, CSR, SNR )
            CALL DROT( 2, B(1,1), 1, B(1,2), 1, CSR, SNR )
C
C           Compute inf norms of A and B.
C
            H1 = MAX( ABS( A(1,1) ) + ABS( A(1,2) ),
     $                ABS( A(2,1) ) + ABS( A(2,2) ) )
            H2 = MAX( ABS( B(1,1) ) + ABS( B(1,2) ),
     $                ABS( B(2,1) ) + ABS( B(2,2) ) )
C
            IF( ( SCALE1*H1 ).GE.ABS( WR1 )*H2 ) THEN
C
C              Find left rotation matrix Q to zero out B(2,1).
C
               CALL DLARTG( B(1,1), B(2,1), CSL, SNL, R )
C
            ELSE
C
C              Find left rotation matrix Q to zero out A(2,1).
C
               CALL DLARTG( A(1,1), A(2,1), CSL, SNL, R )
C
            END IF
C
            CALL DROT( 2, A(1,1), LDA, A(2,1), LDA, CSL, SNL )
            CALL DROT( 2, B(1,1), LDB, B(2,1), LDB, CSL, SNL )
C
            A(2,1) = ZERO
            B(2,1) = ZERO
C
C           Re-adjoint B.
C
            R = B(1,1)
            B(1,1) = B(2,2)
            B(2,2) = R
            B(1,2) = -B(1,2)
C
         ELSE
C
C           A pair of complex conjugate eigenvalues:
C           first compute the SVD of the matrix adj(B).
C
            R = B(1,1)
            B(1,1) = B(2,2)
            B(2,2) = R
            B(1,2) = -B(1,2)
            CALL DLASV2( B(1,1), B(1,2), B(2,2), R, T, SNL, CSL,
     $                   SNR, CSR )
C
C           Form (A,B) := Q(A,adj(B))Z' where Q is left rotation matrix
C           and Z is right rotation matrix computed from DLASV2.
C
            CALL DROT( 2, A(1,1), LDA, A(2,1), LDA, CSL, SNL )
            CALL DROT( 2, B(1,1), LDB, B(2,1), LDB, CSR, SNR )
            CALL DROT( 2, A(1,1), 1, A(1,2), 1, CSR, SNR )
            CALL DROT( 2, B(1,1), 1, B(1,2), 1, CSL, SNL )
C
            B(2,1) = ZERO
            B(1,2) = ZERO
         END IF
C
      END IF
C
C     Unscaling
C
      R = B(1,1)
      T = B(2,2)
      A(1,1) = ANORM*A(1,1)
      A(2,1) = ANORM*A(2,1)
      A(1,2) = ANORM*A(1,2)
      A(2,2) = ANORM*A(2,2)
      B(1,1) = BNORM*B(1,1)
      B(2,1) = BNORM*B(2,1)
      B(1,2) = BNORM*B(1,2)
      B(2,2) = BNORM*B(2,2)
C
      IF( WI.EQ.ZERO ) THEN
         ALPHAR(1) = A(1,1)
         ALPHAR(2) = A(2,2)
         ALPHAI(1) = ZERO
         ALPHAI(2) = ZERO
         BETA(1) = B(1,1)
         BETA(2) = B(2,2)
      ELSE
         WR1 = ANORM*WR1
         WI  = ANORM*WI
         IF ( ABS( WR1 ).GT.ONE .OR. WI.GT.ONE ) THEN
            WR1 = WR1*R
            WI = WI*R
            R = ONE
         END IF
         IF ( ABS( WR1 ).GT.ONE .OR. ABS( WI ).GT.ONE ) THEN
            WR1 = WR1*T
            WI = WI*T
            T = ONE
         END IF
         ALPHAR(1) = ( WR1 / SCALE1 )*R*T
         ALPHAI(1) = ABS( ( WI  / SCALE1 )*R*T )
         ALPHAR(2) =  ALPHAR(1)
         ALPHAI(2) = -ALPHAI(1)
         BETA(1) = BNORM
         BETA(2) = BNORM
      END IF
      RETURN
C *** Last line of MB03YT ***
      END
