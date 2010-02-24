      SUBROUTINE AB09MD( DICO, JOB, EQUIL, ORDSEL, N, M, P, NR, ALPHA,
     $                   A, LDA, B, LDB, C, LDC, NS, HSV, TOL, IWORK,
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
C     To compute a reduced order model (Ar,Br,Cr) for an original
C     state-space representation (A,B,C) by using either the square-root
C     or the balancing-free square-root Balance & Truncate (B & T)
C     model reduction method for the ALPHA-stable part of the system.
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
C     JOB     CHARACTER*1
C             Specifies the model reduction approach to be used
C             as follows:
C             = 'B':  use the square-root Balance & Truncate method;
C             = 'N':  use the balancing-free square-root
C                     Balance & Truncate method.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the triplet (A,B,C) as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     ORDSEL  CHARACTER*1
C             Specifies the order selection method as follows:
C             = 'F':  the resulting order NR is fixed;
C             = 'A':  the resulting order NR is automatically determined
C                     on basis of the given tolerance TOL.
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
C             On entry with ORDSEL = 'F', NR is the desired order of the
C             resulting reduced order system.  0 <= NR <= N.
C             On exit, if INFO = 0, NR is the order of the resulting
C             reduced order model. For a system with NU ALPHA-unstable
C             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N),
C             NR is set as follows: if ORDSEL = 'F', NR is equal to
C             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order
C             on entry, and NMIN is the order of a minimal realization
C             of the ALPHA-stable part of the given system; NMIN is
C             determined as the number of Hankel singular values greater
C             than NS*EPS*HNORM(As,Bs,Cs), where EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH) and
C             HNORM(As,Bs,Cs) is the Hankel norm of the ALPHA-stable
C             part of the given system (computed in HSV(1));
C             if ORDSEL = 'A', NR is the sum of NU and the number of
C             Hankel singular values greater than
C             MAX(TOL,NS*EPS*HNORM(As,Bs,Cs)).
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the ALPHA-stability boundary for the eigenvalues
C             of the state dynamics matrix A. For a continuous-time
C             system (DICO = 'C'), ALPHA <= 0 is the boundary value for
C             the real parts of eigenvalues, while for a discrete-time
C             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the
C             boundary value for the moduli of eigenvalues.
C             The ALPHA-stability domain does not include the boundary.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the state dynamics matrix A.
C             On exit, if INFO = 0, the leading NR-by-NR part of this
C             array contains the state dynamics matrix Ar of the reduced
C             order system.
C             The resulting A has a block-diagonal form with two blocks.
C             For a system with NU ALPHA-unstable eigenvalues and
C             NS ALPHA-stable eigenvalues (NU+NS = N), the leading
C             NU-by-NU block contains the unreduced part of A
C             corresponding to ALPHA-unstable eigenvalues in an
C             upper real Schur form.
C             The trailing (NR+NS-N)-by-(NR+NS-N) block contains
C             the reduced part of A corresponding to ALPHA-stable
C             eigenvalues.
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
C     NS      (output) INTEGER
C             The dimension of the ALPHA-stable subsystem.
C
C     HSV     (output) DOUBLE PRECISION array, dimension (N)
C             If INFO = 0, the leading NS elements of HSV contain the
C             Hankel singular values of the ALPHA-stable part of the
C             original system ordered decreasingly.
C             HSV(1) is the Hankel norm of the ALPHA-stable subsystem.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If ORDSEL = 'A', TOL contains the tolerance for
C             determining the order of reduced system.
C             For model reduction, the recommended value is
C             TOL = c*HNORM(As,Bs,Cs), where c is a constant in the
C             interval [0.00001,0.001], and HNORM(As,Bs,Cs) is the
C             Hankel-norm of the ALPHA-stable part of the given system
C             (computed in HSV(1)).
C             If TOL <= 0 on entry, the used default value is
C             TOL = NS*EPS*HNORM(As,Bs,Cs), where NS is the number of
C             ALPHA-stable eigenvalues of A and EPS is the machine
C             precision (see LAPACK Library Routine DLAMCH).
C             This value is appropriate to compute a minimal realization
C             of the ALPHA-stable part.
C             If ORDSEL = 'F', the value of TOL is ignored.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK = 0, if JOB = 'B';
C             LIWORK = N, if JOB = 'N'.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2).
C             For optimum performance LDWORK should be larger.
C
C     Warning Indicator
C
C     IWARN   INTEGER
C             = 0:  no warning;
C             = 1:  with ORDSEL = 'F', the selected order NR is greater
C                   than NSMIN, the sum of the order of the
C                   ALPHA-unstable part and the order of a minimal
C                   realization of the ALPHA-stable part of the given
C                   system. In this case, the resulting NR is set equal
C                   to NSMIN.
C             = 2:  with ORDSEL = 'F', the selected order NR is less
C                   than the order of the ALPHA-unstable part of the
C                   given system. In this case NR is set equal to the
C                   order of the ALPHA-unstable part.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  the computation of the ordered real Schur form of A
C                   failed;
C             = 2:  the separation of the ALPHA-stable/unstable diagonal
C                   blocks failed because of very close eigenvalues;
C             = 3:  the computation of Hankel singular values failed.
C
C     METHOD
C
C     Let be the following linear system
C
C          d[x(t)] = Ax(t) + Bu(t)
C          y(t)    = Cx(t)                               (1)
C
C     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1)
C     for a discrete-time system. The subroutine AB09MD determines for
C     the given system (1), the matrices of a reduced order system
C
C          d[z(t)] = Ar*z(t) + Br*u(t)
C          yr(t)   = Cr*z(t)                             (2)
C
C     such that
C
C     HSV(NR+NS-N) <= INFNORM(G-Gr) <= 2*[HSV(NR+NS-N+1)+...+HSV(NS)],
C
C     where G and Gr are transfer-function matrices of the systems
C     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the
C     infinity-norm of G.
C
C     The following procedure is used to reduce a given G:
C
C     1) Decompose additively G as
C
C          G = G1 + G2
C
C        such that G1 = (As,Bs,Cs) has only ALPHA-stable poles and
C        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles.
C
C     2) Determine G1r, a reduced order approximation of the
C        ALPHA-stable part G1.
C
C     3) Assemble the reduced model Gr as
C
C           Gr = G1r + G2.
C
C     To reduce the ALPHA-stable part G1, if JOB = 'B', the square-root
C     Balance & Truncate method of [1] is used, and for an ALPHA-stable
C     continuous-time system (DICO = 'C'), the resulting reduced model
C     is balanced. For ALPHA-stable systems, setting TOL < 0, the
C     routine can be used to compute balanced minimal state-space
C     realizations.
C
C     If JOB = 'N', the balancing-free square-root version of the
C     Balance & Truncate method [2] is used to reduce the ALPHA-stable
C     part G1.
C
C     REFERENCES
C
C     [1] Tombs M.S. and Postlethwaite I.
C         Truncated balanced realization of stable, non-minimal
C         state-space systems.
C         Int. J. Control, Vol. 46, pp. 1319-1330, 1987.
C
C     [2] Varga A.
C         Efficient minimal realization procedure based on balancing.
C         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991,
C         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.),
C         Vol. 2, pp. 42-46.
C
C     NUMERICAL ASPECTS
C
C     The implemented methods rely on accuracy enhancing square-root or
C     balancing-free square-root techniques.
C                                         3
C     The algorithms require less than 30N  floating point operations.
C
C     CONTRIBUTOR
C
C     C. Oara, University "Politehnica" Bucharest.
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     February 1999. Based on the RASP routines SADSDC, SRBT and SRBFT.
C
C     REVISIONS
C
C     Mar. 1999, V. Sima, Research Institute for Informatics, Bucharest.
C     Nov. 2000, A. Varga, DLR Oberpfaffenhofen.
C
C     KEYWORDS
C
C     Balancing, minimal realization, model reduction, multivariable
C     system, state-space model, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, C100
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, C100 = 100.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         DICO, EQUIL, JOB, ORDSEL
      INTEGER           INFO, IWARN, LDA, LDB, LDC, LDWORK, M, N, NR,
     $                  NS, P
      DOUBLE PRECISION  ALPHA, TOL
C     .. Array Arguments ..
      INTEGER           IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), HSV(*)
C     .. Local Scalars ..
      LOGICAL           DISCR, FIXORD
      INTEGER           IERR, IWARNL, KT, KTI, KU, KW, KWI, KWR, LWR,
     $                  NN, NRA, NU, NU1, WRKOPT
      DOUBLE PRECISION  ALPWRK, MAXRED
C     .. External Functions ..
      LOGICAL           LSAME
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB09AX, TB01ID, TB01KD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO   = 0
      IWARN  = 0
      DISCR  = LSAME( DICO,   'D' )
      FIXORD = LSAME( ORDSEL, 'F' )
C
C     Test the input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSAME( JOB, 'B' ) .OR. LSAME( JOB, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( EQUIL, 'S' ) .OR.
     $                 LSAME( EQUIL, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( FIXORD .OR. LSAME( ORDSEL, 'A' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( P.LT.0 ) THEN
         INFO = -7
      ELSE IF( FIXORD .AND. ( NR.LT.0 .OR. NR.GT.N ) ) THEN
         INFO = -8
      ELSE IF( ( DISCR .AND. ( ALPHA.LT.ZERO .OR. ALPHA.GT.ONE ) ) .OR.
     $    ( .NOT.DISCR .AND.   ALPHA.GT.ZERO ) ) THEN
         INFO = -9
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( LDWORK.LT.MAX( 1, N*( 2*N + MAX( N, M, P ) + 5 ) +
     $                         ( N*( N + 1 ) )/2 ) ) THEN
         INFO = -21
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MIN( N, M, P ).EQ.0 ) THEN
         NR = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
      IF( LSAME( EQUIL, 'S' ) ) THEN
C
C        Scale simultaneously the matrices A, B and C:
C        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a
C        diagonal matrix.
C        Workspace: N.
C
         MAXRED = C100
         CALL TB01ID( 'All', N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                DWORK, INFO )
      END IF
C
C     Correct the value of ALPHA to ensure stability.
C
      ALPWRK = ALPHA
      IF( DISCR ) THEN
         IF( ALPHA.EQ.ONE ) ALPWRK = ONE - SQRT( DLAMCH( 'E' ) )
      ELSE
         IF( ALPHA.EQ.ZERO ) ALPWRK = -SQRT( DLAMCH( 'E' ) )
      END IF
C
C     Allocate working storage.
C
      NN  = N*N
      KU  = 1
      KWR = KU + NN
      KWI = KWR + N
      KW  = KWI + N
      LWR = LDWORK - KW + 1
C
C     Reduce A to a block-diagonal real Schur form, with the
C     ALPHA-unstable part in the leading diagonal position, using a
C     non-orthogonal similarity transformation A <- inv(T)*A*T and
C     apply the transformation to B and C: B <- inv(T)*B and C <- C*T.
C
C     Workspace needed:      N*(N+2);
C     Additional workspace:  need   3*N;
C                            prefer larger.
C
      CALL TB01KD( DICO, 'Unstable', 'General', N, M, P, ALPWRK, A, LDA,
     $             B, LDB, C, LDC, NU, DWORK(KU), N, DWORK(KWR),
     $             DWORK(KWI), DWORK(KW), LWR, IERR )
C
      IF( IERR.NE.0 ) THEN
         IF( IERR.NE.3 ) THEN
            INFO = 1
         ELSE
            INFO = 2
         END IF
         RETURN
      END IF
C
      WRKOPT = DWORK(KW) + DBLE( KW-1 )
C
      IWARNL = 0
      NS = N - NU
      IF( FIXORD ) THEN
         NRA = MAX( 0, NR-NU )
         IF( NR.LT.NU )
     $      IWARNL = 2
      ELSE
         NRA = 0
      END IF
C
C     Finish if only unstable part is present.
C
      IF( NS.EQ.0 ) THEN
         NR = NU
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      NU1 = NU + 1
C
C     Allocate working storage.
C
      KT  = 1
      KTI = KT  + NN
      KW  = KTI + NN
C
C     Compute a B & T approximation of the stable part.
C     Workspace: need   N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2;
C                prefer larger.
C
      CALL AB09AX( DICO, JOB, ORDSEL, NS, M, P, NRA, A(NU1,NU1), LDA,
     $             B(NU1,1), LDB, C(1,NU1), LDC, HSV, DWORK(KT), N,
     $             DWORK(KTI), N, TOL, IWORK, DWORK(KW), LDWORK-KW+1,
     $             IWARN, IERR )
      IWARN = MAX( IWARN, IWARNL )
C
      IF( IERR.NE.0 ) THEN
         INFO = IERR + 1
         RETURN
      END IF
C
      NR = NRA + NU
C
      DWORK(1) = MAX( WRKOPT, INT( DWORK(KW) ) + KW - 1 )
C
      RETURN
C *** Last line of AB09MD ***
      END
