      SUBROUTINE MB04YW( QRIT, UPDATU, UPDATV, M, N, L, K, SHIFT, D, E,
     $                   U, LDU, V, LDV, DWORK )
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
C     To perform either one QR or QL iteration step onto the unreduced
C     bidiagonal submatrix Jk:
C
C              |D(l) E(l)    0  ...    0   |
C              | 0   D(l+1) E(l+1)     .   |
C         Jk = | .                     .   |
C              | .                     .   |
C              | .                   E(k-1)|
C              | 0   ...        ...   D(k) |
C
C     with k <= p and l >= 1, p = MIN(M,N), of the bidiagonal matrix J:
C
C              |D(1) E(1)  0    ...   0   |
C              | 0   D(2) E(2)        .   |
C          J = | .                    .   |.
C              | .                    .   |
C              | .                  E(p-1)|
C              | 0   ...        ...  D(p) |
C
C     Hereby, Jk is transformed to  S' Jk T with S and T products of
C     Givens rotations. These Givens rotations S (respectively, T) are
C     postmultiplied into U (respectively, V), if UPDATU (respectively,
C     UPDATV) is .TRUE..
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     QRIT    LOGICAL
C             Indicates whether a QR or QL iteration step is to be
C             taken (from larger end diagonal element towards smaller),
C             as follows:
C             = .TRUE. :  QR iteration step (chase bulge from top to
C                         bottom);
C             = .FALSE.:  QL iteration step (chase bulge from bottom to
C                         top).
C
C     UPDATU  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix U the left-hand Givens rotations S, as follows:
C             = .FALSE.:  Do not form U;
C             = .TRUE. :  The given matrix U is updated (postmultiplied)
C                         by the left-hand Givens rotations S.
C
C     UPDATV  LOGICAL
C             Indicates whether the user wishes to accumulate in a
C             matrix V the right-hand Givens rotations S, as follows:
C             = .FALSE.:  Do not form V;
C             = .TRUE. :  The given matrix V is updated (postmultiplied)
C                         by the right-hand Givens rotations T.
C
C     Input/Output Parameters
C
C     M       (input) INTEGER
C             The number of rows of the matrix U.  M >= 0.
C
C     N       (input) INTEGER
C             The number of rows of the matrix V.  N >= 0.
C
C     L       (input) INTEGER
C             The index of the first diagonal entry of the considered
C             unreduced bidiagonal submatrix Jk of J.
C
C     K       (input) INTEGER
C             The index of the last diagonal entry of the considered
C             unreduced bidiagonal submatrix Jk of J.
C
C     SHIFT   (input) DOUBLE PRECISION
C             Value of the shift used in the QR or QL iteration step.
C
C     D       (input/output) DOUBLE PRECISION array, dimension (p)
C             where p = MIN(M,N)
C             On entry, D must contain the diagonal entries of the
C             bidiagonal matrix J.
C             On exit, D contains the diagonal entries of the
C             transformed bidiagonal matrix S' J T.
C
C     E       (input/output) DOUBLE PRECISION array, dimension (p-1)
C             On entry, E must contain the superdiagonal entries of J.
C             On exit, E contains the superdiagonal entries of the
C             transformed matrix S' J T.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,p)
C             On entry, if UPDATU = .TRUE., U must contain the M-by-p
C             left transformation matrix.
C             On exit, if UPDATU = .TRUE., the Givens rotations S on the
C             left have been postmultiplied into U, i.e., U * S is
C             returned.
C             U is not referenced if UPDATU = .FALSE..
C
C     LDU     INTEGER
C             The leading dimension of the array U.
C             LDU >= max(1,M) if UPDATU = .TRUE.;
C             LDU >= 1        if UPDATU = .FALSE..
C
C     V       (input/output) DOUBLE PRECISION array, dimension (LDV,p)
C             On entry, if UPDATV = .TRUE., V must contain the N-by-p
C             right transformation matrix.
C             On exit, if UPDATV = .TRUE., the Givens rotations T on the
C             right have been postmultiplied into V, i.e., V * T is
C             returned.
C             V is not referenced if UPDATV = .FALSE..
C
C     LDV     INTEGER
C             The leading dimension of the array V.
C             LDV >= max(1,N) if UPDATV = .TRUE.;
C             LDV >= 1        if UPDATV = .FALSE..
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (MAX(1,LDWORK))
C             LDWORK >= 4*MIN(M,N)-4, if UPDATU = UPDATV = .TRUE.;
C             LDWORK >= 2*MIN(M,N)-2, if
C                             UPDATU = .TRUE. and UPDATV = .FALSE. or
C                             UPDATV = .TRUE. and UPDATU = .FALSE.;
C             LDWORK >= 1, if UPDATU = UPDATV = .FALSE..
C
C     METHOD
C
C     QR iterations diagonalize the bidiagonal matrix by zeroing the
C     super-diagonal elements of Jk from bottom to top.
C     QL iterations diagonalize the bidiagonal matrix by zeroing the
C     super-diagonal elements of Jk from top to bottom.
C     The routine overwrites Jk with the bidiagonal matrix S' Jk T,
C     where S and T are products of Givens rotations.
C     T is essentially the orthogonal matrix that would be obtained by
C     applying one implicit symmetric shift QR (QL) step onto the matrix
C     Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a
C     product of an orthogonal matrix T and a upper (lower) triangular
C     matrix. See [1,Sec.8.2-8.3] and [2] for more details.
C
C     REFERENCES
C
C     [1] Golub, G.H. and Van Loan, C.F.
C         Matrix Computations.
C         The Johns Hopkins University Press, Baltimore, Maryland, 1983.
C
C     [2] Bowdler, H., Martin, R.S. and Wilkinson, J.H.
C         The QR and QL algorithms for symmetric matrices.
C         Numer. Math., 11, pp. 293-306, 1968.
C
C     [3] Demmel, J. and Kahan, W.
C         Computing small singular values of bidiagonal matrices with
C         guaranteed high relative accuracy.
C         SIAM J. Sci. Statist. Comput., 11, pp. 873-912, 1990.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997.
C     Supersedes Release 2.0 routines MB04QY and MB04QZ by S. Van
C     Huffel, Katholieke University Leuven, Belgium.
C     This subroutine is based on the QR/QL step implemented in LAPACK
C     routine DBDSQR.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Bidiagonal matrix, orthogonal transformation, singular values.
C
C     ******************************************************************
C
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      LOGICAL           QRIT, UPDATU, UPDATV
      INTEGER           K, L, LDU, LDV, M, N
      DOUBLE PRECISION  SHIFT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  D( * ), DWORK( * ), E( * ), U( LDU, * ),
     $                  V( LDV, * )
C     ..
C     .. Local Scalars ..
      INTEGER           I, IROT, NCV, NM1, NM12, NM13
      DOUBLE PRECISION  COSL, COSR, CS, F, G, H, OLDCS, OLDSN, R, SINL,
     $                  SINR, SN
C     ..
C     .. External Subroutines ..
      EXTERNAL          DLARTG, DLASR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, SIGN
C     ..
C     .. Executable Statements ..
C
C     For speed, no tests of the input scalar arguments are done.
C
C     Quick return if possible.
C
      NCV = MIN( M, N )
      IF ( NCV.LE.1 .OR. L.EQ.K )
     $   RETURN
C
      NM1 = NCV - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IF ( .NOT.UPDATV ) THEN
         NM12 = 0
         NM13 = NM1
      END IF
C
C     If SHIFT = 0, do simplified QR iteration.
C
      IF( SHIFT.EQ.ZERO ) THEN
         IF( QRIT ) THEN
C
C           Chase bulge from top to bottom.
C           Save cosines and sines for later U and/or V updates,
C           if needed.
C
            CS = ONE
            OLDCS = ONE
            CALL DLARTG( D( L )*CS, E( L ), CS, SN, R )
            CALL DLARTG( OLDCS*R, D( L+1 )*SN, OLDCS, OLDSN, D( L ) )
            IF ( UPDATV ) THEN
               DWORK( 1 ) = CS
               DWORK( 1+NM1 ) = SN
            END IF
            IF ( UPDATU ) THEN
               DWORK( 1+NM12 ) = OLDCS
               DWORK( 1+NM13 ) = OLDSN
            END IF
            IROT = 1
C
            DO 110 I = L + 1, K - 1
               CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
               E( I-1 ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
               IROT = IROT + 1
               IF ( UPDATV ) THEN
                  DWORK( IROT ) = CS
                  DWORK( IROT+NM1 ) = SN
               END IF
               IF ( UPDATU ) THEN
                  DWORK( IROT+NM12 ) = OLDCS
                  DWORK( IROT+NM13 ) = OLDSN
               END IF
  110       CONTINUE
C
            H = D( K )*CS
            D( K ) = H*OLDCS
            E( K-1 ) = H*OLDSN
C
C           Update U and/or V.
C
            IF( UPDATV )
     $         CALL DLASR( 'R', 'V', 'F', N, K-L+1, DWORK( 1 ),
     $                     DWORK( NCV ), V( 1, L ), LDV )
            IF( UPDATU )
     $         CALL DLASR( 'R', 'V', 'F', M, K-L+1, DWORK( NM12+1 ),
     $                     DWORK( NM13+1 ), U( 1, L ), LDU )
C
         ELSE
C
C           Chase bulge from bottom to top.
C           Save cosines and sines for later U and/or V updates,
C           if needed.
C
            CS = ONE
            OLDCS = ONE
            CALL DLARTG( D( K )*CS, E( K-1 ), CS, SN, R )
            CALL DLARTG( OLDCS*R, D( K-1 )*SN, OLDCS, OLDSN, D( K ) )
            IF ( UPDATV ) THEN
               DWORK( K-L ) = OLDCS
               DWORK( K-L+NM1 ) = -OLDSN
            END IF
            IF ( UPDATU ) THEN
               DWORK( K-L+NM12 ) = CS
               DWORK( K-L+NM13 ) = -SN
            END IF
            IROT = K - L
C
            DO 120 I = K - 1, L + 1, -1
               CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
               E( I ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
               IROT = IROT - 1
               IF ( UPDATV ) THEN
                  DWORK( IROT ) = OLDCS
                  DWORK( IROT+NM1 ) = -OLDSN
               END IF
               IF ( UPDATU ) THEN
                  DWORK( IROT+NM12 ) = CS
                  DWORK( IROT+NM13 ) = -SN
               END IF
  120       CONTINUE
C
            H = D( L )*CS
            D( L ) = H*OLDCS
            E( L ) = H*OLDSN
C
C           Update U and/or V.
C
            IF( UPDATV )
     $         CALL DLASR( 'R', 'V', 'B', N, K-L+1, DWORK( 1 ),
     $                     DWORK( NCV ), V( 1, L ), LDV )
            IF( UPDATU )
     $         CALL DLASR( 'R', 'V', 'B', M, K-L+1, DWORK( NM12+1 ),
     $                     DWORK( NM13+1 ), U( 1, L ), LDU )
         END IF
      ELSE
C
C        Use nonzero shift.
C
         IF( QRIT ) THEN
C
C           Chase bulge from top to bottom.
C           Save cosines and sines for later U and/or V updates,
C           if needed.
C
            F = ( ABS( D( L ) ) - SHIFT )*
     $          ( SIGN( ONE, D( L ) ) + SHIFT / D( L ) )
            G = E( L )
            CALL DLARTG( F, G, COSR, SINR, R )
            F = COSR*D( L ) + SINR*E( L )
            E( L ) = COSR*E( L ) - SINR*D( L )
            G = SINR*D( L+1 )
            D( L+1 ) = COSR*D( L+1 )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( L ) = R
            F = COSL*E( L ) + SINL*D( L+1 )
            D( L+1 ) = COSL*D( L+1 ) - SINL*E( L )
            G = SINL*E( L+1 )
            E( L+1 ) = COSL*E( L+1 )
            IF ( UPDATV ) THEN
               DWORK( 1 ) = COSR
               DWORK( 1+NM1 ) = SINR
            END IF
            IF ( UPDATU ) THEN
               DWORK( 1+NM12 ) = COSL
               DWORK( 1+NM13 ) = SINL
            END IF
            IROT = 1
C
            DO 130 I = L + 1, K - 2
               CALL DLARTG( F, G, COSR, SINR, R )
               E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               G = SINL*E( I+1 )
               E( I+1 ) = COSL*E( I+1 )
               IROT = IROT + 1
               IF ( UPDATV ) THEN
                  DWORK( IROT ) = COSR
                  DWORK( IROT+NM1 ) = SINR
               END IF
               IF ( UPDATU ) THEN
                  DWORK( IROT+NM12 ) = COSL
                  DWORK( IROT+NM13 ) = SINL
               END IF
  130       CONTINUE
C
            IF ( L.LT.K-1 ) THEN
               CALL DLARTG( F, G, COSR, SINR, R )
               E( K-2 ) = R
               F = COSR*D( K-1 ) + SINR*E( K-1 )
               E( K-1 ) = COSR*E( K-1 ) - SINR*D( K-1 )
               G = SINR*D( K )
               D( K ) = COSR*D( K )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( K-1 ) = R
               F = COSL*E( K-1 ) + SINL*D( K )
               D( K ) = COSL*D( K ) - SINL*E( K-1 )
               IROT = IROT + 1
               IF ( UPDATV ) THEN
                  DWORK( IROT ) = COSR
                  DWORK( IROT+NM1 ) = SINR
               END IF
               IF ( UPDATU ) THEN
                  DWORK( IROT+NM12 ) = COSL
                  DWORK( IROT+NM13 ) = SINL
               END IF
            END IF
            E( K-1 ) = F
C
C           Update U and/or V.
C
            IF( UPDATV )
     $         CALL DLASR( 'R', 'V', 'F', N, K-L+1, DWORK( 1 ),
     $                     DWORK( NCV ), V( 1, L ), LDV )
            IF( UPDATU )
     $         CALL DLASR( 'R', 'V', 'F', M, K-L+1, DWORK( NM12+1 ),
     $                     DWORK( NM13+1 ), U( 1, L ), LDU )
C
         ELSE
C
C           Chase bulge from bottom to top.
C           Save cosines and sines for later U and/or V updates,
C           if needed.
C
            F = ( ABS( D( K ) ) - SHIFT )*
     $          ( SIGN( ONE, D( K ) ) + SHIFT / D( K ) )
            G = E( K-1 )
            IF ( L.LT.K-1 ) THEN
               CALL DLARTG( F, G, COSR, SINR, R )
               F = COSR*D( K ) + SINR*E( K-1 )
               E( K-1 ) = COSR*E( K-1 ) - SINR*D( K )
               G = SINR*D( K-1 )
               D( K-1 ) = COSR*D( K-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( K ) = R
               F = COSL*E( K-1 ) + SINL*D( K-1 )
               D( K-1 ) = COSL*D( K-1 ) - SINL*E( K-1 )
               G = SINL*E( K-2 )
               E( K-2 ) = COSL*E( K-2 )
               IF ( UPDATV ) THEN
                  DWORK( K-L ) = COSL
                  DWORK( K-L+NM1 ) = -SINL
               END IF
               IF ( UPDATU ) THEN
                  DWORK( K-L+NM12 ) = COSR
                  DWORK( K-L+NM13 ) = -SINR
               END IF
               IROT = K - L
            ELSE
               IROT = K - L + 1
            END IF
C
            DO 140 I = K - 1, L + 2, -1
               CALL DLARTG( F, G, COSR, SINR, R )
               E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               G = SINL*E( I-2 )
               E( I-2 ) = COSL*E( I-2 )
               IROT = IROT - 1
               IF ( UPDATV ) THEN
                  DWORK( IROT ) = COSL
                  DWORK( IROT+NM1 ) = -SINL
               END IF
               IF ( UPDATU ) THEN
                  DWORK( IROT+NM12 ) = COSR
                  DWORK( IROT+NM13 ) = -SINR
               END IF
  140       CONTINUE
C
            CALL DLARTG( F, G, COSR, SINR, R )
            E( L+1 ) = R
            F = COSR*D( L+1 ) + SINR*E( L )
            E( L ) = COSR*E( L ) - SINR*D( L+1 )
            G = SINR*D( L )
            D( L ) = COSR*D( L )
            CALL DLARTG( F, G, COSL, SINL, R )
            D( L+1 ) = R
            F = COSL*E( L ) + SINL*D( L )
            D( L ) = COSL*D( L ) - SINL*E( L )
            IROT = IROT - 1
            IF ( UPDATV ) THEN
               DWORK( IROT ) = COSL
               DWORK( IROT+NM1 ) = -SINL
            END IF
            IF ( UPDATU ) THEN
               DWORK( IROT+NM12 ) = COSR
               DWORK( IROT+NM13 ) = -SINR
            END IF
            E( L ) = F
C
C           Update U and/or V if desired.
C
            IF( UPDATV )
     $         CALL DLASR( 'R', 'V', 'B', N, K-L+1, DWORK( 1 ),
     $                     DWORK( NCV ), V( 1, L ), LDV )
            IF( UPDATU )
     $         CALL DLASR( 'R', 'V', 'B', M, K-L+1, DWORK( NM12+1 ),
     $                     DWORK( NM13+1 ), U( 1, L ), LDU )
         END IF
      END IF
C
      RETURN
C *** Last line of MB04YW ***
      END
