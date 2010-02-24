      SUBROUTINE AG08BD( EQUIL, L, N, M, P, A, LDA, E, LDE, B, LDB,
     $                   C, LDC, D, LDD, NFZ, NRANK, NIZ, DINFZ, NKROR,
     $                   NINFE, NKROL, INFZ, KRONR, INFE, KRONL,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
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
C     To extract from the system pencil
C
C                       ( A-lambda*E B )
C           S(lambda) = (              )
C                       (      C     D )
C
C     a regular pencil Af-lambda*Ef which has the finite Smith zeros of
C     S(lambda) as generalized eigenvalues. The routine also computes
C     the orders of the infinite Smith zeros and determines the singular
C     and infinite Kronecker structure of system pencil, i.e., the right
C     and left Kronecker indices, and the multiplicities of infinite
C     eigenvalues.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to balance the system
C             matrix as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
C
C     Input/Output Parameters
C
C     L       (input) INTEGER
C             The number of rows of matrices A, B, and E.  L >= 0.
C
C     N       (input) INTEGER
C             The number of columns of matrices A, E, and C.  N >= 0.
C
C     M       (input) INTEGER
C             The number of columns of matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The number of rows of matrix C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading L-by-N part of this array must
C             contain the state dynamics matrix A of the system.
C             On exit, the leading NFZ-by-NFZ part of this array
C             contains the matrix Af of the reduced pencil.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,L).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading L-by-N part of this array must
C             contain the descriptor matrix E of the system.
C             On exit, the leading NFZ-by-NFZ part of this array
C             contains the matrix Ef of the reduced pencil.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,L).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading L-by-M part of this array must
C             contain the input/state matrix B of the system.
C             On exit, this matrix does not contain useful information.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,L) if M > 0;
C             LDB >= 1        if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C of the system.
C             On exit, this matrix does not contain useful information.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading P-by-M part of this array must contain the
C             direct transmission matrix D of the system.
C
C     LDD     INTEGER
C             The leading dimension of array D.  LDD >= MAX(1,P).
C
C     NFZ     (output) INTEGER
C             The number of finite zeros.
C
C     NRANK   (output) INTEGER
C             The normal rank of the system pencil.
C
C     NIZ     (output) INTEGER
C             The number of infinite zeros.
C
C     DINFZ   (output) INTEGER
C             The maximal multiplicity of infinite Smith zeros.
C
C     NKROR   (output) INTEGER
C             The number of right Kronecker indices.
C
C     NINFE   (output) INTEGER
C             The number of elementary infinite blocks.
C
C     NKROL   (output) INTEGER
C             The number of left Kronecker indices.
C
C     INFZ    (output) INTEGER array, dimension (N+1)
C             The leading DINFZ elements of INFZ contain information
C             on the infinite elementary divisors as follows:
C             the system has INFZ(i) infinite elementary divisors of
C             degree i in the Smith form, where i = 1,2,...,DINFZ.
C
C     KRONR   (output) INTEGER array, dimension (N+M+1)
C             The leading NKROR elements of this array contain the
C             right Kronecker (column) indices.
C
C     INFE    (output) INTEGER array, dimension (1+MIN(L+P,N+M))
C             The leading NINFE elements of INFE contain the
C             multiplicities of infinite eigenvalues.
C
C     KRONL   (output) INTEGER array, dimension (L+P+1)
C             The leading NKROL elements of this array contain the
C             left Kronecker (row) indices.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             If the user sets TOL <= 0, then default tolerances are
C             used instead, as follows: TOLDEF = L*N*EPS in TG01FD
C             (to determine the rank of E) and TOLDEF = (L+P)*(N+M)*EPS
C             in the rest, where EPS is the machine precision
C             (see LAPACK Library routine DLAMCH).  TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension N+max(1,M)
C             On output, IWORK(1) contains the normal rank of the
C             transfer function matrix.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= max( 4*(L+N), LDW ), if EQUIL = 'S',
C             LDWORK >= LDW,                 if EQUIL = 'N', where
C             LDW = max(L+P,M+N)*(M+N) + max(1,5*max(L+P,M+N)).
C             For optimum performance LDWORK should be larger.
C
C             If LDWORK = -1, then a workspace query is assumed;
C             the routine only calculates the optimal size of the
C             DWORK array, returns this value as the first entry of
C             the DWORK array, and no error message related to LDWORK
C             is issued by XERBLA.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value.
C
C     METHOD
C
C     The routine extracts from the system matrix of a descriptor
C     system (A-lambda*E,B,C,D) a regular pencil Af-lambda*Ef which
C     has the finite zeros of the system as generalized eigenvalues.
C     The procedure has the following main computational steps:
C
C        (a) construct the (L+P)-by-(N+M) system pencil
C
C             S(lambda) = ( B  A )-lambda*( 0  E );
C                         ( D  C )        ( 0  0 )
C
C        (b) reduce S(lambda) to S1(lambda) with the same finite
C            zeros and right Kronecker structure but with E
C            upper triangular and nonsingular;
C
C        (c) reduce S1(lambda) to S2(lambda) with the same finite
C            zeros and right Kronecker structure but with D of
C            full row rank;
C
C        (d) reduce S2(lambda) to S3(lambda) with the same finite zeros
C            and with D square invertible;
C
C        (e) perform a unitary transformation on the columns of
C
C            S3(lambda) = (A-lambda*E   B) in order to reduce it to
C                         (     C       D)
C
C            (Af-lambda*Ef   X), with Y and Ef square invertible;
C            (     0         Y)
C
C        (f) compute the right and left Kronecker indices of the system
C            matrix, which together with the multiplicities of the
C            finite and infinite eigenvalues constitute the
C            complete set of structural invariants under strict
C            equivalence transformations of a linear system.
C
C     REFERENCES
C
C     [1] P. Misra, P. Van Dooren and A. Varga.
C         Computation of structural invariants of generalized
C         state-space systems.
C         Automatica, 30, pp. 1921-1936, 1994.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable (see [1]).
C
C     FURTHER COMMENTS
C
C     In order to compute the finite Smith zeros of the system
C     explicitly, a call to this routine may be followed by a
C     call to the LAPACK Library routines DGEGV or DGGEV.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     May 1999.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Sep. 1999,
C     Jan. 2009, Mar. 2009, Apr. 2009.
C     A. Varga, DLR Oberpfaffenhofen, Nov. 1999, Feb. 2002, Mar. 2002.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, multivariable
C     system, orthogonal transformation, structural invariant.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL
      INTEGER           DINFZ, INFO, L, LDA, LDB, LDC, LDD, LDE, LDWORK,
     $                  M, N, NFZ, NINFE, NIZ, NKROL, NKROR, NRANK, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INFE(*), INFZ(*), IWORK(*), KRONL(*), KRONR(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), D(LDD,*),
     $                  DWORK(*), E(LDE,*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LQUERY
      INTEGER           I, I0, I1, II, IPD, ITAU, J, JWORK, KABCD,
     $                  LABCD2, LDABCD, LDW, MM, MU, N2, NB, NN, NSINFE,
     $                  NU, NUMU, PP, WRKOPT
      DOUBLE PRECISION  SVLMAX, TOLER
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          AG08BY, DLACPY, DLASET, DORMRZ, DTZRZF, MA02BD,
     $                  MA02CD, TB01XD, TG01AD, TG01FD, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
      LDABCD = MAX( L+P, N+M )
      LABCD2 = LDABCD*( N+M )
      LEQUIL = LSAME( EQUIL, 'S' )
      LQUERY = ( LDWORK.EQ.-1 )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -1
      ELSE IF( L.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, L ) ) THEN
         INFO = -7
      ELSE IF( LDE.LT.MAX( 1, L ) ) THEN
         INFO = -9
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.L ) ) THEN
         INFO = -11
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -15
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -27
      ELSE
         I0  = MIN( L+P, M+N )
         I1  = MIN( L, N )
         II  = MIN( M, P )
         LDW = LABCD2 + MAX( 1, 5*LDABCD )
         IF( LEQUIL )
     $      LDW = MAX( 4*( L + N ), LDW )
         IF( LQUERY ) THEN
            CALL TG01FD( 'N', 'N', 'N', L, N, M, P, A, LDA, E, LDE, B,
     $                   LDB, C, LDC, DUM, 1, DUM, 1, NN, N2, TOL,
     $                   IWORK, DWORK, -1, INFO )
            WRKOPT = MAX( LDW, INT( DWORK(1) ) )
            SVLMAX = ZERO
            CALL AG08BY( .TRUE., I1, M+N, P+L, SVLMAX, DWORK, LDABCD+I1,
     $                   E, LDE, NU, MU, NIZ, DINFZ, NKROL, INFZ, KRONL,
     $                   TOL, IWORK, DWORK, -1, INFO )
            WRKOPT = MAX( WRKOPT, LABCD2 + INT( DWORK(1) ) )
            CALL AG08BY( .FALSE., I1, II, M+N, SVLMAX, DWORK, LDABCD+I1,
     $                   E, LDE, NU, MU, NIZ, DINFZ, NKROL, INFZ, KRONL,
     $                   TOL, IWORK, DWORK, -1, INFO )
            WRKOPT = MAX( WRKOPT, LABCD2 + INT( DWORK(1) ) )
            NB = ILAENV( 1, 'ZGERQF', ' ', II, I1+II, -1, -1 )
            WRKOPT = MAX( WRKOPT, LABCD2 + II + II*NB )
            NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RT', I1, I1+II, II,
     $                            -1 ) )
            WRKOPT = MAX( WRKOPT, LABCD2 + II + I1*NB )
         ELSE IF( LDWORK.LT.LDW ) THEN
            INFO = -30
         END IF
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AG08BD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      NIZ = 0
      NKROL = 0
      NKROR = 0
C
C     Quick return if possible.
C
      IF( MAX( L, N, M, P ).EQ.0 ) THEN
         NFZ = 0
         DINFZ = 0
         NINFE = 0
         NRANK = 0
         IWORK(1) = 0
         DWORK(1) = ONE
         RETURN
      END IF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      WRKOPT = 1
      KABCD  = 1
      JWORK  = KABCD + LABCD2
C
C     If required, balance the system pencil.
C     Workspace: need   4*(L+N).
C
      IF( LEQUIL ) THEN
         CALL TG01AD( 'A', L, N, M, P, ZERO, A, LDA, E, LDE, B, LDB,
     $                C, LDC, DWORK, DWORK(L+1), DWORK(L+N+1), INFO )
         WRKOPT = 4*(L+N)
      END IF
      SVLMAX = DLANGE( 'Frobenius', L, N, E, LDE, DWORK )
C
C     Reduce the system matrix to QR form,
C
C          ( A11-lambda*E11 A12 B1 )
C          (     A21        A22 B2 ) ,
C          (     C1         C2  D  )
C
C     with E11 invertible and upper triangular.
C     Real workspace: need   max( 1, N+P, min(L,N)+max(3*N-1,M,L) );
C                     prefer larger.
C     Integer workspace: N.
C
      CALL TG01FD( 'N', 'N', 'N', L, N, M, P, A, LDA, E, LDE, B, LDB,
     $             C, LDC, DUM, 1, DUM, 1, NN, N2, TOL, IWORK, DWORK,
     $             LDWORK, INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
C
C     Construct the system pencil
C
C                          MM         NN
C                      ( B1 A12 A11-lambda*E11 ) NN
C        S1(lambda) =  ( B2 A22      A21       ) L-NN
C                      ( D  C2       C1        ) P
C
C     of dimension (L+P)-by-(M+N).
C     Workspace: need  LABCD2 = max( L+P, N+M )*( N+M ).
C
      N2 = N - NN
      MM = M + N2
      PP = P + ( L - NN )
      CALL DLACPY( 'Full', L, M, B, LDB, DWORK(KABCD), LDABCD )
      CALL DLACPY( 'Full', P, M, D, LDD, DWORK(KABCD+L), LDABCD )
      CALL DLACPY( 'Full', L, N2, A(1,NN+1), LDA,
     $              DWORK(KABCD+LDABCD*M), LDABCD )
      CALL DLACPY( 'Full', P, N2, C(1,NN+1), LDC,
     $              DWORK(KABCD+LDABCD*M+L), LDABCD )
      CALL DLACPY( 'Full', L, NN, A, LDA,
     $              DWORK(KABCD+LDABCD*MM), LDABCD )
      CALL DLACPY( 'Full', P, NN, C, LDC,
     $              DWORK(KABCD+LDABCD*MM+L), LDABCD )
C
C     If required, set tolerance.
C
      TOLER = TOL
      IF( TOLER.LE.ZERO ) THEN
         TOLER = DBLE( ( L + P )*( M + N ) ) * DLAMCH( 'Precision' )
      END IF
      SVLMAX = MAX( SVLMAX,
     $              DLANGE( 'Frobenius', NN+PP, NN+MM, DWORK(KABCD),
     $                      LDABCD, DWORK(JWORK) ) )
C
C     Extract the reduced pencil S2(lambda)
C
C             ( Bc  Ac-lambda*Ec )
C             ( Dc      Cc       )
C
C     having the same finite Smith zeros as the system pencil
C     S(lambda) but with Dc, a MU-by-MM full row rank
C     left upper trapezoidal matrix, and Ec, an NU-by-NU
C     upper triangular nonsingular matrix.
C
C     Real workspace: need   max( min(P+L,M+N)+max(min(L,N),3*(M+N)-1),
C                                  5*(P+L), 1 ) + LABCD2;
C                     prefer larger.
C     Integer workspace: MM, MM <= M+N; PP <= P+L.
C
      CALL AG08BY( .TRUE., NN, MM, PP, SVLMAX, DWORK(KABCD), LDABCD,
     $             E, LDE, NU, MU, NIZ, DINFZ, NKROL, INFZ, KRONL,
     $             TOLER, IWORK, DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
      WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C     Set the number of simple (nondynamic) infinite eigenvalues
C     and the normal rank of the system pencil.
C
      NSINFE = MU
      NRANK  = NN + MU
C
C     Pertranspose the system.
C
      CALL TB01XD( 'D', NU, MM, MM, MAX( 0, NU-1 ), MAX( 0, NU-1 ),
     $              DWORK(KABCD+LDABCD*MM), LDABCD,
     $              DWORK(KABCD), LDABCD,
     $              DWORK(KABCD+LDABCD*MM+NU), LDABCD,
     $              DWORK(KABCD+NU), LDABCD, INFO )
      CALL MA02BD( 'Right', NU+MM, MM, DWORK(KABCD), LDABCD )
      CALL MA02BD( 'Left',  MM, NU+MM, DWORK(KABCD+NU), LDABCD )
      CALL MA02CD( NU, 0, MAX( 0, NU-1 ), E, LDE )
C
      IF( MU.NE.MM ) THEN
         NN = NU
         PP = MM
         MM = MU
         KABCD = KABCD + ( PP - MM )*LDABCD
C
C        Extract the reduced pencil S3(lambda),
C
C             ( Br  Ar-lambda*Er ) ,
C             ( Dr      Cr       )
C
C        having the same finite Smith zeros as the pencil S(lambda),
C        but with Dr, an MU-by-MU invertible upper triangular matrix,
C        and Er, an NU-by-NU upper triangular nonsingular matrix.
C
C        Workspace: need   max( 1, 5*(M+N) ) + LABCD2.
C                   prefer larger.
C        No integer workspace necessary.
C
         CALL AG08BY( .FALSE., NN, MM, PP, SVLMAX, DWORK(KABCD), LDABCD,
     $                E, LDE, NU, MU, I0, I1, NKROR, IWORK, KRONR,
     $                TOLER, IWORK, DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
      END IF
C
      IF( NU.NE.0 ) THEN
C
C        Perform a unitary transformation on the columns of
C                     ( Br Ar-lambda*Er )
C                     ( Dr     Cr       )
C        in order to reduce it to
C                     ( *  Af-lambda*Ef )
C                     ( Y       0       )
C        with Y and Ef square invertible.
C
C        Compute Af by reducing  ( Br Ar ) to  ( *  Af ) .
C                                ( Dr Cr )     ( Y   0 )
C
         NUMU  = NU + MU
         IPD   = KABCD + NU
         ITAU  = JWORK
         JWORK = ITAU + MU
C
C        Workspace: need   LABCD2 + 2*min(M,P);
C                   prefer LABCD2 + min(M,P) + min(M,P)*NB.
C
         CALL DTZRZF( MU, NUMU, DWORK(IPD), LDABCD, DWORK(ITAU),
     $                DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C        Workspace: need   LABCD2 + min(M,P) + min(L,N);
C                   prefer LABCD2 + min(M,P) + min(L,N)*NB.
C
         CALL DORMRZ( 'Right', 'Transpose', NU, NUMU, MU, NU,
     $                DWORK(IPD), LDABCD, DWORK(ITAU), DWORK(KABCD),
     $                LDABCD, DWORK(JWORK), LDWORK-JWORK+1, INFO )
         WRKOPT = MAX( WRKOPT, INT( DWORK(JWORK) ) + JWORK - 1 )
C
C        Save Af.
C
         CALL DLACPY( 'Full', NU, NU, DWORK(KABCD+LDABCD*MU), LDABCD, A,
     $                LDA )
C
C        Compute Ef by applying the saved transformations from previous
C        reduction to ( 0  Er ) .
C
         CALL DLASET( 'Full', NU, MU, ZERO, ZERO, DWORK(KABCD), LDABCD )
         CALL DLACPY( 'Full', NU, NU, E, LDE, DWORK(KABCD+LDABCD*MU),
     $                LDABCD )
C
         CALL DORMRZ( 'Right', 'Transpose', NU, NUMU, MU, NU,
     $                DWORK(IPD), LDABCD, DWORK(ITAU), DWORK(KABCD),
     $                LDABCD, DWORK(JWORK), LDWORK-JWORK+1, INFO )
C
C        Save Ef.
C
         CALL DLACPY( 'Full', NU, NU, DWORK(KABCD+LDABCD*MU), LDABCD, E,
     $                LDE )
      END IF
C
      NFZ = NU
C
C     Set right Kronecker indices (column indices).
C
      DO 10 I = 1, NKROR
         IWORK(I) = KRONR(I)
   10 CONTINUE
C
      J = 0
      DO 30 I = 1, NKROR
         DO 20 II = J + 1, J + IWORK(I)
            KRONR(II) = I - 1
   20    CONTINUE
         J = J + IWORK(I)
   30 CONTINUE
C
      NKROR = J
C
C     Set left Kronecker indices (row indices).
C
      DO 40 I = 1, NKROL
         IWORK(I) = KRONL(I)
   40 CONTINUE
C
      J = 0
      DO 60 I = 1, NKROL
         DO 50 II = J + 1, J + IWORK(I)
            KRONL(II) = I - 1
   50    CONTINUE
         J = J + IWORK(I)
   60 CONTINUE
C
      NKROL = J
C
C     Determine the number of simple infinite blocks
C     as the difference between the number of infinite blocks
C     of order greater than one and the order of Dr.
C
      NINFE = 0
      DO 70 I = 1, DINFZ
         NINFE = NINFE + INFZ(I)
   70 CONTINUE
      NINFE = NSINFE - NINFE
      DO 80 I = 1, NINFE
         INFE(I) = 1
   80 CONTINUE
C
C     Set the structure of infinite eigenvalues.
C
      DO 100 I = 1, DINFZ
         DO 90 II = NINFE + 1, NINFE + INFZ(I)
            INFE(II) = I + 1
   90    CONTINUE
         NINFE = NINFE + INFZ(I)
  100 CONTINUE
C
      IWORK(1) = NSINFE
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of AG08BD ***
      END
