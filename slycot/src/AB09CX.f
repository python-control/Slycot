      SUBROUTINE AB09CX( DICO, ORDSEL, N, M, P, NR, A, LDA, B, LDB,
     $                   C, LDC, D, LDD, HSV, TOL1, TOL2, IWORK,
     $                   DWORK, LDWORK, IWARN, INFO )
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
C     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable
C     original state-space representation (A,B,C,D) by using the optimal
C     Hankel-norm approximation method in conjunction with square-root
C     balancing. The state dynamics matrix A of the original system is
C     an upper quasi-triangular matrix in real Schur canonical form.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the original system as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL1.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the original state-space representation, i.e.
C             the order of the matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     NR      (input/output) INTEGER
C             On entry with ORDSEL = 'F', NR is the desired order of
C             the resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. NR is set as follows:
C             if ORDSEL = 'F', NR is equal to MIN(MAX(0,NR-KR+1),NMIN),
C             where KR is the multiplicity of the Hankel singular value
C             HSV(NR+1), NR is the desired order on entry, and NMIN is
C             the order of a minimal realization of the given system;
C             NMIN is determined as the number of Hankel singular values
C             greater than N*EPS*HNORM(A,B,C), where EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH) and
C             HNORM(A,B,C) is the Hankel norm of the system (computed
C             in HSV(1));
C             if ORDSEL = 'A', NR is equal to the number of Hankel
C             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)).
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A in a real Schur
C             canonical form.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the
C             reduced order system in a real Schur form.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the original input/state matrix B.
C             On exit, if INFO = 0, the leading NR-by-M part of this
C             array contains the input/state matrix Br of the reduced
C             order system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original state/output matrix C.
C             On exit, if INFO = 0, the leading P-by-NR part of this
C             array contains the state/output matrix Cr of the reduced
C             order system.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M)
C             On entry, the leading P-by-M part of this array must
C             contain the original input/output matrix D.
C             On exit, if INFO = 0, the leading P-by-M part of this
C             array contains the input/output matrix Dr of the reduced
C             order system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, it contains the Hankel singular values of
C             the original system ordered decreasingly. HSV(1) is the
C             Hankel norm of the system.
C
C     Tolerances
C
C     TOL1    DOUBLE PRECISION
C             If ORDSEL = 'A', TOL1 contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL1 = c*HNORM(A,B,C), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(A,B,C) is the
C             Hankel-norm of the given system (computed in HSV(1)).
C             For computing a minimal realization, the recommended
C             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the
C             machine precision (see LAPACK Library Routine DLAMCH).
C             This value is used by default if TOL1 <= 0 on entry.
C             If ORDSEL = 'F', the value of TOL1 is ignored.
C
C     TOL2    DOUBLE PRECISION
C             The tolerance for determining the order of a minimal
C             realization of the given system. The recommended value is
C             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default
C             if TOL2 <= 0 on entry.
C             If TOL2 > 0, then TOL2 <= TOL1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = MAX(1,M),   if DICO = 'C';
C             LIWORK = MAX(1,N,M), if DICO = 'D'.
C             On exit, if INFO = 0, IWORK(1) contains NMIN, the order of
C             the computed minimal realization.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( LDW1,LDW2 ), where
C             LDW1 = N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2,
C             LDW2 = N*(M+P+2) + 2*M*P + MIN(N,M) +
C                    MAX( 3*M+1, MIN(N,M)+P ).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than the order of a minimal realization of the
C                   given system. In this case, the resulting NR is set
C                   automatically to a value corresponding to the order
C                   of a minimal realization of the system.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the state matrix A is not stable (if DICO = 'C')
C                   or not convergent (if DICO = 'D');
C             = 2:  the computation of Hankel singular values failed;
C             = 3:  the computation of stable projection failed;
C             = 4:  the order of computed stable projection differs
C                   from the order of Hankel-norm approximation.
C
C     METHOD
C
C     Let be the stable linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t) + Du(t)                           (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09CX determines for
C     the given system (1), the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t) + Dr*u(t)                       (2)
C
C     such that
C
C           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     The optimal Hankel-norm approximation method of [1], based on the
C     square-root balancing projection formulas of [2], is employed.
C
C     REFERENCES
C
C     [1] Glover, K.
C         All optimal Hankel norm approximation of linear
C         multivariable systems and their L-infinity error bounds.
C         Int. J. Control, Vol. 36, pp. 1145-1193, 1984.
C
C     [2] Tombs M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on an accuracy enhancing square-root
C     technique.
C                                         3
C     The algorithms require less than 30N  floating point operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center,
C     DLR Oberpfaffenhofen, April 1998.
C     Based on the RASP routine OHNAP1.
C
C     REVISIONS
C
C     November 11, 1998, V. Sima, Research Institute for Informatics,
C     Bucharest.
C     April 24, 2000, A. Varga, DLR Oberpfaffenhofen.
C     April  8, 2001, A. Varga, DLR Oberpfaffenhofen.
C     March 26, 2005, V. Sima, Research Institute for Informatics.
C
C     KEYWORDS
C
C     Balancing, Hankel-norm approximation, model reduction,
C     multivariable system, state-space model.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDD, LDWORK,
     $                  M, N, NR, P
      DOUBLE PRECISION  TOL1, TOL2
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), HSV(*)
C     .. Local Scalars
      LOGICAL           DISCR, FIXORD
      INTEGER           I, I1, IERR, IRANK, J, KB1, KB2, KC1, KC2T,
     $                  KHSVP, KHSVP2, KR, KT, KTI, KU, KW, KW1, KW2,
     $                  LDB1, LDB2, LDC1, LDC2T, NA, NDIM, NKR1, NMINR,
     $                  NR1, NU, WRKOPT
      DOUBLE PRECISION  ATOL, RTOL, SKP, SKP2, SRRTOL
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB04MD, AB09AX, DAXPY, DCOPY, DGELSY, DGEMM,
     $                  DLACPY, DSWAP, MA02AD, MB01SD, TB01KD, TB01WD,
     $                  XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
C
C     Check the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( TOL2.GT.ZERO .AND. TOL2.GT.TOL1 ) THEN
         INFO = -17
      ELSE IF( LDWORK.LT.MAX( N*( 2*N + MAX( N, M, P ) + 5 ) +
     $                        ( N*( N + 1 ) )/2,
     $                        N*( M + P + 2 ) + 2*M*P + MIN( N, M ) +
     $                        MAX ( 3*M + 1, MIN( N, M ) + P ) ) ) THEN
         INFO = -20
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09CX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      RTOL   = DBLE( N )*DLAMCH( 'Epsilon' )
      SRRTOL = SQRT( RTOL )
C
C     Allocate working storage.
C
      KT  = 1
      KTI = KT  + N*N
      KW  = KTI + N*N
C
C     Compute a minimal order balanced realization of the given system.
C     Workspace: need   N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2;
C                prefer larger.
C
      CALL AB09AX( DICO, 'Balanced', 'Automatic', N, M, P, NMINR, A,
     $             LDA, B, LDB, C, LDC, HSV, DWORK(KT), N, DWORK(KTI),
     $             N, TOL2, IWORK, DWORK(KW), LDWORK-KW+1, IWARN, INFO )
C
      IF( INFO.NE.0 )
     $   RETURN
      WRKOPT = INT( DWORK(KW) ) + KW - 1
C
C     Compute the order of reduced system.
C
      ATOL = RTOL*HSV(1)
      IF( FIXORD ) THEN
         IF( NR.GT.0 ) THEN
            IF( NR.GT.NMINR ) THEN
               NR = NMINR
               IWARN = 1
            ENDIF
         ENDIF
      ELSE
         ATOL = MAX( TOL1, ATOL )
         NR = 0
         DO 10 I = 1, NMINR
            IF( HSV(I).LE.ATOL ) GO TO 20
            NR = NR + 1
   10    CONTINUE
   20    CONTINUE
      ENDIF
C
      IF( NR.EQ.NMINR ) THEN
         IWORK(1) = NMINR
         DWORK(1) = WRKOPT
         KW = N*(N+2)+1
C
C        Reduce Ar to a real Schur form.
C
         CALL TB01WD( NMINR, M, P, A, LDA, B, LDB, C, LDC,
     $                DWORK(2*N+1), N, DWORK, DWORK(N+1), DWORK(KW),
     $                LDWORK-KW+1, IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
         RETURN
      END IF
      SKP = HSV(NR+1)
C
C     If necessary, reduce the order such that HSV(NR) > HSV(NR+1).
C
   30 IF( NR.GT.0 ) THEN
         IF( ABS( HSV(NR)-SKP ).LE.SRRTOL*SKP ) THEN
            NR = NR - 1
            GO TO 30
         END IF
      END IF
C
C     Determine KR, the multiplicity of HSV(NR+1).
C
      KR = 1
      DO 40 I = NR+2, NMINR
         IF( ABS( HSV(I)-SKP ).GT.SRRTOL*SKP ) GO TO 50
         KR = KR + 1
   40 CONTINUE
   50 CONTINUE
C
C     For discrete-time case, apply the discrete-to-continuous bilinear
C     transformation.
C
      IF( DISCR ) THEN
C
C        Workspace: need   N;
C                   prefer larger.
C
         CALL AB04MD( 'Discrete', NMINR, M, P, ONE, ONE, A, LDA, B, LDB,
     $                C, LDC, D, LDD, IWORK, DWORK, LDWORK, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      END IF
C
C     Define leading dimensions and offsets for temporary data.
C
      NU     = NMINR - NR - KR
      NA     = NR + NU
      LDB1   = NA
      LDC1   = P
      LDB2   = KR
      LDC2T  = MAX( KR, M )
      NR1    = NR  + 1
      NKR1   = MIN( NMINR, NR1 + KR )
C
      KHSVP  = 1
      KHSVP2 = KHSVP  + NA
      KU     = KHSVP2 + NA
      KB1    = KU     + P*M
      KB2    = KB1    + LDB1*M
      KC1    = KB2    + LDB2*M
      KC2T   = KC1    + LDC1*NA
      KW     = KC2T   + LDC2T*P
C
C     Save B2 and C2'.
C
      CALL DLACPY( 'Full', KR, M, B(NR1,1), LDB, DWORK(KB2), LDB2 )
      CALL MA02AD( 'Full', P, KR, C(1,NR1), LDC, DWORK(KC2T), LDC2T )
      IF( NR.GT.0 ) THEN
C
C        Permute the elements of HSV and of matrices A, B, C.
C
         CALL DCOPY( NR, HSV(1), 1, DWORK(KHSVP), 1 )
         CALL DCOPY( NU, HSV(NKR1), 1, DWORK(KHSVP+NR), 1 )
         CALL DLACPY( 'Full', NMINR, NU, A(1,NKR1), LDA, A(1,NR1), LDA )
         CALL DLACPY( 'Full', NU, NA, A(NKR1,1), LDA, A(NR1,1), LDA )
         CALL DLACPY( 'Full', NU, M, B(NKR1,1), LDB, B(NR1,1), LDB )
         CALL DLACPY( 'Full', P, NU, C(1,NKR1), LDC, C(1,NR1), LDC )
C
C        Save B1 and C1.
C
         CALL DLACPY( 'Full', NA, M, B, LDB, DWORK(KB1), LDB1 )
         CALL DLACPY( 'Full', P, NA, C, LDC, DWORK(KC1), LDC1 )
      END IF
C
C     Compute U = C2*pinv(B2').
C     Workspace: need   N*(M+P+2) + 2*M*P +
C                       max(min(KR,M)+3*M+1,2*min(KR,M)+P);
C                prefer N*(M+P+2) + 2*M*P +
C                       max(min(KR,M)+2*M+(M+1)*NB,2*min(KR,M)+P*NB),
C                where  NB  is the maximum of the block sizes for
C                DGEQP3, DTZRZF, DTZRQF, DORMQR, and DORMRZ.
C
      DO 55 J = 1, M
         IWORK(J) = 0
   55 CONTINUE
      CALL DGELSY( KR, M, P, DWORK(KB2), LDB2, DWORK(KC2T), LDC2T,
     $             IWORK, RTOL, IRANK, DWORK(KW), LDWORK-KW+1, IERR )
      WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
      CALL MA02AD( 'Full', M, P, DWORK(KC2T), LDC2T, DWORK(KU), P )
C
C     Compute D <- D + HSV(NR+1)*U.
C
      I = KU
      DO 60 J = 1, M
         CALL DAXPY( P, SKP, DWORK(I), 1, D(1,J), 1 )
         I = I + P
   60 CONTINUE
C
      IF( NR.GT.0 ) THEN
         SKP2 = SKP*SKP
C
C        Compute G = inv(S1*S1-skp*skp*I), where S1 is the diagonal
C        matrix of relevant singular values (of order NMINR - KR).
C
         I1 = KHSVP2
         DO 70 I = KHSVP, KHSVP+NA-1
            DWORK(I1) = ONE / ( DWORK(I)*DWORK(I) - SKP2 )
            I1 = I1 + 1
   70    CONTINUE
C
C        Compute C <- C1*S1-skp*U*B1'.
C
         CALL MB01SD( 'Column', P, NA, C, LDC, DWORK, DWORK(KHSVP) )
         CALL DGEMM( 'NoTranspose', 'Transpose', P, NA, M, -SKP,
     $               DWORK(KU), P, DWORK(KB1), LDB1, ONE, C, LDC )
C
C        Compute B <- G*(S1*B1-skp*C1'*U).
C
         CALL MB01SD( 'Row', NA, M, B, LDB, DWORK(KHSVP), DWORK )
         CALL DGEMM( 'Transpose', 'NoTranspose', NA, M, P, -SKP,
     $               DWORK(KC1), LDC1, DWORK(KU), P, ONE, B, LDB )
         CALL MB01SD( 'Row', NA, M, B, LDB, DWORK(KHSVP2), DWORK )
C
C        Compute A <- -A1' - B*B1'.
C
         DO 80 J = 2, NA
            CALL DSWAP( J-1, A(1,J), 1, A(J,1), LDA )
   80    CONTINUE
         CALL DGEMM( 'NoTranspose', 'Transpose', NA, NA, M, -ONE, B,
     $               LDB, DWORK(KB1), LDB1, -ONE, A, LDA )
C
C        Extract stable part.
C        Workspace:  need   N*N+5*N;
C                    prefer larger.
C
         KW1 = NA*NA + 1
         KW2 = KW1 + NA
         KW  = KW2 + NA
         CALL TB01KD( 'Continuous', 'Stability', 'General', NA, M, P,
     $                ZERO, A, LDA, B, LDB, C, LDC, NDIM, DWORK, NA,
     $                DWORK(KW1), DWORK(KW2), DWORK(KW), LDWORK-KW+1,
     $                IERR )
         IF( IERR.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
C
         IF( NDIM.NE.NR ) THEN
            INFO = 4
            RETURN
         END IF
         WRKOPT = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
C        For discrete-time case, apply the continuous-to-discrete
C        bilinear transformation.
C
         IF( DISCR )
     $      CALL AB04MD( 'Continuous', NR, M, P, ONE, ONE, A, LDA, B,
     $                   LDB, C, LDC, D, LDD, IWORK, DWORK, LDWORK,
     $                   INFO )
      END IF
      IWORK(1) = NMINR
      DWORK(1) = WRKOPT
C
      RETURN
C *** Last line of AB09CX ***
      END
