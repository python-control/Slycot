      SUBROUTINE TB03AY( NR, A, LDA, INDBLK, NBLK, VCOEFF, LDVCO1,
     $                   LDVCO2, PCOEFF, LDPCO1, LDPCO2, INFO )
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
C     To calculate the (PWORK-by-NR) polynomial matrix V(s) one
C     (PWORK-by-NBLK(L-1)) block V:L-1(s) at a time, in reverse order
C     (L = INDBLK,...,1).  At each stage, the (NBLK(L)-by-NBLK(L)) poly-
C     nomial matrix W(s) = V2(s) * A2 is formed, where V2(s) is that
C     part of V(s) already computed and A2 is the subdiagonal (incl.)
C     part of the L-th column block of A; W(s) is temporarily stored in
C     the top left part of P(s), as is subsequently the further matrix
C     Wbar(s) = s * V:L(s) - W(s).  Then, except for the final stage
C     L = 1 (when the next step is to calculate P(s) itself, not here),
C     the top left part of V:L-1(s) is given by Wbar(s) * inv(R), where
C     R is the upper triangular part of the L-th superdiagonal block of
C     A.  Finally, note that the coefficient matrices W(.,.,K) can only
C     be non-zero for K = L + 1,...,INPLUS, with each of these matrices
C     having only its first NBLK(L-1) rows non-trivial.  Similarly,
C     Wbar(.,.,K) (and so clearly V:L-1(.,.,K) ) can only be non-zero
C     for K = L,...,INPLUS, with each of these having only its first
C     NBLK(K-1) rows non-trivial except for K = L, which has NBLK(L)
C     such rows.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Coprime matrix fraction, elementary polynomial operations,
C     polynomial matrix, state-space representation, transfer matrix.
C
C     NOTE: In the interests of speed, this routine does not check the
C           inputs for errors.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INDBLK, INFO, LDA, LDPCO1, LDPCO2, LDVCO1,
     $                  LDVCO2, NR
C     .. Array Arguments ..
      INTEGER           NBLK(*)
      DOUBLE PRECISION  A(LDA,*), PCOEFF(LDPCO1,LDPCO2,*),
     $                  VCOEFF(LDVCO1,LDVCO2,*)
C     .. Local Scalars ..
      INTEGER           I, INPLUS, IOFF, J, JOFF, K, KPLUS, L, LSTART,
     $                  LSTOP, LWORK, NCOL, NROW
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMM, DLACPY, DSCAL, DTRSM
C     .. Executable Statements ..
C
      INFO = 0
      INPLUS = INDBLK + 1
      JOFF = NR
C
C     Calculate each column block V:LWORK-1(s) of V(s) in turn.
C
      DO 70 L = 1, INDBLK
         LWORK = INPLUS - L
C
C        Determine number of columns of V:LWORK(s) & its position in V.
C
         NCOL = NBLK(LWORK)
         JOFF = JOFF - NCOL
C
C        Find limits for V2(s) * A2 calculation: skips zero rows
C        in V(s).
C
         LSTART = JOFF + 1
         LSTOP = JOFF
C
C        Calculate W(s) and store (temporarily) in top left part
C        of P(s).
C
         DO 10 K = LWORK + 1, INPLUS
            NROW = NBLK(K-1)
            LSTOP = LSTOP + NROW
            CALL DGEMM( 'No transpose', 'No transpose', NROW, NCOL,
     $                  LSTOP-LSTART+1, ONE, VCOEFF(1,LSTART,K), LDVCO1,
     $                  A(LSTART,JOFF+1), LDA, ZERO, PCOEFF(1,1,K),
     $                  LDPCO1 )
   10    CONTINUE
C
C        Replace W(s) by Wbar(s) = s * V:L(s) - W(s).
C
         NROW = NCOL
C
         DO 30 K = LWORK, INDBLK
            KPLUS = K + 1
C
            DO 20 J = 1, NCOL
               CALL DSCAL( NROW, -ONE, PCOEFF(1,J,K), 1 )
               CALL DAXPY( NROW, ONE, VCOEFF(1,JOFF+J,KPLUS), 1,
     $                     PCOEFF(1,J,K), 1 )
   20       CONTINUE
C
            NROW = NBLK(K)
   30    CONTINUE
C
         DO 40 J = 1, NCOL
            CALL DSCAL( NROW, -ONE, PCOEFF(1,J,INPLUS), 1 )
   40    CONTINUE
C
         IF ( LWORK.NE.1 ) THEN
C
C           If not final stage, use the upper triangular R (from A)
C           to calculate V:L-1(s), finally storing this new block.
C
            IOFF = JOFF - NBLK(LWORK-1)
C
            DO 50 I = 1, NCOL
               IF ( A(IOFF+I,JOFF+I).EQ.ZERO ) THEN
C
C                 Error return.
C
                  INFO = I
                  RETURN
               END IF
   50       CONTINUE
C
            NROW = NBLK(LWORK)
C
            DO 60 K = LWORK, INPLUS
               CALL DLACPY( 'Full', NROW, NCOL, PCOEFF(1,1,K), LDPCO1,
     $                      VCOEFF(1,IOFF+1,K), LDVCO1 )
               CALL DTRSM( 'Right', 'Upper', 'No Transpose', 'Non-unit',
     $                     NROW, NCOL, ONE, A(IOFF+1,JOFF+1), LDA,
     $                     VCOEFF(1,IOFF+1,K), LDVCO1 )
               NROW = NBLK(K)
   60       CONTINUE
C
         END IF
   70 CONTINUE
C
      RETURN
C *** Last line of TB03AY ***
      END
