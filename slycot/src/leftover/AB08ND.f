      SUBROUTINE AB08ND( EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                   NU, RANK, DINFZ, NKROR, NKROL, INFZ, KRONR,
     $                   KRONL, AF, LDAF, BF, LDBF, TOL, IWORK, DWORK,
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
C     To construct for a linear multivariable system described by a
C     state-space model (A,B,C,D) a regular pencil (A - lambda*B ) which
C                                                    f          f
C     has the invariant zeros of the system as generalized eigenvalues.
C     The routine also computes the orders of the infinite zeros and the
C     right and left Kronecker indices of the system (A,B,C,D).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to balance the compound
C             matrix (see METHOD) as follows:
C             = 'S':  Perform balancing (scaling);
C             = 'N':  Do not perform balancing.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The number of state variables, i.e., the order of the
C             matrix A.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs.  M >= 0.
C
C     P       (input) INTEGER
C             The number of system outputs.  P >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A of the system.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input/state matrix B of the system.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             state/output matrix C of the system.
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
C     NU      (output) INTEGER
C             The number of (finite) invariant zeros.
C
C     RANK    (output) INTEGER
C             The normal rank of the transfer function matrix.
C
C     DINFZ   (output) INTEGER
C             The maximum degree of infinite elementary divisors.
C
C     NKROR   (output) INTEGER
C             The number of right Kronecker indices.
C
C     NKROL   (output) INTEGER
C             The number of left Kronecker indices.
C
C     INFZ    (output) INTEGER array, dimension (N)
C             The leading DINFZ elements of INFZ contain information
C             on the infinite elementary divisors as follows:
C             the system has INFZ(i) infinite elementary divisors
C             of degree i, where i = 1,2,...,DINFZ.
C
C     KRONR   (output) INTEGER array, dimension (MAX(N,M)+1)
C             The leading NKROR elements of this array contain the
C             right Kronecker (column) indices.
C
C     KRONL   (output) INTEGER array, dimension (MAX(N,P)+1)
C             The leading NKROL elements of this array contain the
C             left Kronecker (row) indices.
C
C     AF      (output) DOUBLE PRECISION array, dimension
C             (LDAF,N+MIN(P,M))
C             The leading NU-by-NU part of this array contains the
C             coefficient matrix A  of the reduced pencil. The remainder
C                                 f
C             of the leading (N+M)-by-(N+MIN(P,M)) part is used as
C             internal workspace.
C
C     LDAF    INTEGER
C             The leading dimension of array AF.  LDAF >= MAX(1,N+M).
C
C     BF      (output) DOUBLE PRECISION array, dimension (LDBF,N+M)
C             The leading NU-by-NU part of this array contains the
C             coefficient matrix B  of the reduced pencil. The
C                                 f
C             remainder of the leading (N+P)-by-(N+M) part is used as
C             internal workspace.
C
C     LDBF    INTEGER
C             The leading dimension of array BF.  LDBF >= MAX(1,N+P).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             A tolerance used in rank decisions to determine the
C             effective rank, which is defined as the order of the
C             largest leading (or trailing) triangular submatrix in the
C             QR (or RQ) factorization with column (or row) pivoting
C             whose estimated condition number is less than 1/TOL.
C             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS
C             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS,
C             where EPS is the machine precision (see LAPACK Library
C             Routine DLAMCH).
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (MAX(M,P))
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N),
C                               MIN(P,N) + MAX(3*P-1,N+P,N+M),
C                               MIN(M,N) + MAX(3*M-1,N+M) ).
C             An upper bound is MAX(s,N) + MAX(3*s-1,N+s), with
C             s = MAX(M,P).
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
C     The routine extracts from the system matrix of a state-space
C     system (A,B,C,D) a regular pencil A - lambda*B  which has the
C                                        f          f
C     invariant zeros of the system as generalized eigenvalues as
C     follows:
C
C        (a) construct the (N+P)-by-(N+M) compound matrix (B  A);
C                                                         (D  C)
C
C        (b) reduce the above system to one with the same invariant
C            zeros and with D of full row rank;
C
C        (c) pertranspose the system;
C
C        (d) reduce the system to one with the same invariant zeros and
C            with D square invertible;
C
C        (e) perform a unitary transformation on the columns of
C            (A - lambda*I  B) in order to reduce it to
C            (      C       D)
C
C            (A  - lambda*B   X)
C            ( f           f   ), with Y and B  square invertible;
C            (      0         Y)              f
C
C        (f) compute the right and left Kronecker indices of the system
C            (A,B,C,D), which together with the orders of the infinite
C            zeros (determined by steps (a) - (e)) constitute the
C            complete set of structural invariants under strict
C            equivalence transformations of a linear system.
C
C     REFERENCES
C
C     [1] Svaricek, F.
C         Computation of the Structural Invariants of Linear
C         Multivariable Systems with an Extended Version of
C         the Program ZEROS.
C         System & Control Letters, 6, pp. 261-266, 1985.
C
C     [2] Emami-Naeini, A. and Van Dooren, P.
C         Computation of Zeros of Linear Multivariable Systems.
C         Automatica, 18, pp. 415-430, 1982.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is backward stable (see [2] and [1]).
C
C     FURTHER COMMENTS
C
C     In order to compute the invariant zeros of the system explicitly,
C     a call to this routine may be followed by a call to the LAPACK
C     Library routine DGGEV with A = A , B = B  and N = NU.
C                                     f       f
C     If RANK = 0, the routine DGEEV can be used (since B = I).
C                                                        f
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
C     Supersedes Release 2.0 routine AB08BD by F. Svaricek.
C
C     REVISIONS
C
C     Oct. 1997, Feb. 1998, Dec. 2003, March 2004, Jan. 2009, Mar. 2009,
C     Apr. 2009.
C
C     KEYWORDS
C
C     Generalized eigenvalue problem, Kronecker indices, multivariable
C     system, orthogonal transformation, structural invariant.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL
      INTEGER           DINFZ, INFO, LDA, LDAF, LDB, LDBF, LDC, LDD,
     $                  LDWORK, M, N, NKROL, NKROR, NU, P, RANK
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INFZ(*), IWORK(*), KRONL(*), KRONR(*)
      DOUBLE PRECISION  A(LDA,*), AF(LDAF,*), B(LDB,*), BF(LDBF,*),
     $                  C(LDC,*), D(LDD,*), DWORK(*)
C     .. Local Scalars ..
      LOGICAL           LEQUIL, LQUERY
      INTEGER           I, I1, II, J, MM, MNU, MU, NB, NINFZ, NN, NU1,
     $                  NUMU, NUMU1, PP, RO, SIGMA, WRKOPT
      DOUBLE PRECISION  MAXRED, SVLMAX, THRESH, TOLER
C     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, ILAENV, LSAME
C     .. External Subroutines ..
      EXTERNAL          AB08NX, DCOPY, DLACPY, DLASET, DORMRZ, DTZRZF,
     $                  TB01ID, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO = 0
      LEQUIL = LSAME( EQUIL, 'S' )
      LQUERY = ( LDWORK.EQ.-1 )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( P.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDAF.LT.MAX( 1, N + M ) ) THEN
         INFO = -22
      ELSE IF( LDBF.LT.MAX( 1, N + P ) ) THEN
         INFO = -24
      ELSE
         II = MIN( P, M )
         I  = MAX(          II + MAX( 3*M - 1, N ),
     $             MIN( P, N ) + MAX( 3*P - 1, N+P, N+M ),
     $             MIN( M, N ) + MAX( 3*M - 1, N+M ), 1 )
         IF( LQUERY ) THEN
            SVLMAX = ZERO
            NINFZ  = 0
            CALL AB08NX( N, M, P, P, 0, SVLMAX, BF, LDBF, NINFZ, INFZ,
     $                   KRONL, MU, NU, NKROL, TOL, IWORK, DWORK, -1,
     $                   INFO )
            WRKOPT = MAX( I, INT( DWORK(1) ) )
            CALL AB08NX( N, II, M, M-II, II, SVLMAX, AF, LDAF, NINFZ,
     $                   INFZ, KRONL, MU, NU, NKROL, TOL, IWORK, DWORK,
     $                   -1, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
            NB = ILAENV( 1, 'DGERQF', ' ', II, N+II, -1, -1 )
            WRKOPT = MAX( WRKOPT, II + II*NB )
            NB = MIN( 64, ILAENV( 1, 'DORMRQ', 'RT', N, N+II, II, -1 ) )
            WRKOPT = MAX( WRKOPT, II + N*NB )
         ELSE IF( LDWORK.LT.I ) THEN
            INFO = -28
         END IF
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB08ND', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         DWORK(1) = WRKOPT
         RETURN
      END IF
C
      DINFZ = 0
      NKROL = 0
      NKROR = 0
C
C     Quick return if possible.
C
      IF ( N.EQ.0 ) THEN
         IF ( MIN( M, P ).EQ.0 ) THEN
            NU   = 0
            RANK = 0
            DWORK(1) = ONE
            RETURN
         END IF
      END IF
C
      MM = M
      NN = N
      PP = P
C
      DO 20 I = 1, N
         INFZ(I) = 0
   20 CONTINUE
C
      IF ( M.GT.0 ) THEN
         DO 40 I = 1, N + 1
            KRONR(I) = 0
   40    CONTINUE
      END IF
C
      IF ( P.GT.0 ) THEN
         DO 60 I = 1, N + 1
            KRONL(I) = 0
   60    CONTINUE
      END IF
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.)
C
      WRKOPT = 1
C
C     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N).
C                                    ( D  C )
C
      CALL DLACPY( 'Full', NN, MM, B, LDB, BF, LDBF )
      IF ( PP.GT.0 )
     $   CALL DLACPY( 'Full', PP, MM, D, LDD, BF(1+NN,1), LDBF )
      IF ( NN.GT.0 ) THEN
         CALL DLACPY( 'Full', NN, NN, A, LDA, BF(1,1+MM), LDBF )
         IF ( PP.GT.0 )
     $      CALL DLACPY( 'Full', PP, NN, C, LDC, BF(1+NN,1+MM), LDBF )
      END IF
C
C     If required, balance the compound matrix (default MAXRED).
C     Workspace: need   N.
C
      IF ( LEQUIL .AND. NN.GT.0 .AND. PP.GT.0 ) THEN
         MAXRED = ZERO
         CALL TB01ID( 'A', NN, MM, PP, MAXRED, BF(1,1+MM), LDBF, BF,
     $                LDBF, BF(1+NN,1+MM), LDBF, DWORK, INFO )
         WRKOPT = N
      END IF
C
C     If required, set tolerance.
C
      THRESH = SQRT( DBLE( (N + P)*(N + M) ) )*DLAMCH( 'Precision' )
      TOLER = TOL
      IF ( TOLER.LT.THRESH ) TOLER = THRESH
      SVLMAX = DLANGE( 'Frobenius', NN+PP, NN+MM, BF, LDBF, DWORK )
C
C     Reduce this system to one with the same invariant zeros and with
C     D upper triangular of full row rank MU (the normal rank of the
C     original system).
C     Workspace: need   MAX( 1, MIN(P,M) + MAX(3*M-1,N),
C                               MIN(P,N) + MAX(3*P-1,N+P,N+M) );
C                prefer larger.
C
      RO = PP
      SIGMA = 0
      NINFZ = 0
      CALL AB08NX( NN, MM, PP, RO, SIGMA, SVLMAX, BF, LDBF, NINFZ, INFZ,
     $             KRONL, MU, NU, NKROL, TOLER, IWORK, DWORK, LDWORK,
     $             INFO )
      WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
      RANK = MU
C
C     Pertranspose the system.
C
      NUMU = NU + MU
      IF ( NUMU.NE.0 ) THEN
         MNU = MM + NU
         NUMU1 = NUMU + 1
C
         DO 80 I = 1, NUMU
            CALL DCOPY( MNU, BF(I,1), LDBF, AF(1,NUMU1-I), -1 )
   80    CONTINUE
C
         IF ( MU.NE.MM ) THEN
C
C           Here MU < MM and MM > 0 (since MM = 0 implies MU = 0 = MM).
C
            PP = MM
            NN = NU
            MM = MU
C
C           Reduce the system to one with the same invariant zeros and
C           with D square invertible.
C           Workspace: need  MAX( 1, MU + MAX(3*MU-1,N),
C                                 MIN(M,N) + MAX(3*M-1,N+M) );
C                prefer larger. Note that MU <= MIN(P,M).
C
            RO = PP - MM
            SIGMA = MM
            CALL AB08NX( NN, MM, PP, RO, SIGMA, SVLMAX, AF, LDAF, NINFZ,
     $                   INFZ, KRONR, MU, NU, NKROR, TOLER, IWORK,
     $                   DWORK, LDWORK, INFO )
            WRKOPT = MAX( WRKOPT, INT( DWORK(1) ) )
         END IF
C
         IF ( NU.NE.0 ) THEN
C
C           Perform a unitary transformation on the columns of
C                     ( B  A-lambda*I )
C                     ( D       C     )
C           in order to reduce it to
C                     ( X  AF-lambda*BF )
C                     ( Y       0       )
C           with Y and BF square invertible.
C
            CALL DLASET( 'Full', NU, MU, ZERO, ZERO, BF, LDBF )
            CALL DLASET( 'Full', NU, NU, ZERO, ONE,  BF(1,MU+1), LDBF )
C
            IF ( RANK.NE.0 ) THEN
               NU1 = NU + 1
               I1  = NU + MU
C
C              Workspace: need   2*MIN(M,P);
C                         prefer MIN(M,P) + MIN(M,P)*NB.
C
               CALL DTZRZF( MU, I1, AF(NU1,1), LDAF, DWORK, DWORK(MU+1),
     $                      LDWORK-MU, INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(MU+1) ) + MU )
C
C              Workspace: need   MIN(M,P) + N;
C                         prefer MIN(M,P) + N*NB.
C
               CALL DORMRZ( 'Right', 'Transpose', NU, I1, MU, NU,
     $                      AF(NU1,1), LDAF, DWORK, AF, LDAF,
     $                      DWORK(MU+1), LDWORK-MU, INFO )
               WRKOPT = MAX( WRKOPT, INT( DWORK(MU+1) ) + MU )
C
               CALL DORMRZ( 'Right', 'Transpose', NU, I1, MU, NU,
     $                      AF(NU1,1), LDAF, DWORK, BF, LDBF,
     $                      DWORK(MU+1), LDWORK-MU, INFO )
C
            END IF
C
C           Move AF and BF in the first columns. This assumes that
C           DLACPY moves column by column.
C
            CALL DLACPY( 'Full', NU, NU, AF(1,MU+1), LDAF, AF, LDAF )
            IF ( RANK.NE.0 )
     $         CALL DLACPY( 'Full', NU, NU, BF(1,MU+1), LDBF, BF, LDBF )
C
         END IF
      END IF
C
C     Set right Kronecker indices (column indices).
C
      IF ( NKROR.GT.0 ) THEN
         J = 1
C
         DO 120 I = 1, N + 1
C
            DO 100 II = J, J + KRONR(I) - 1
               IWORK(II) = I - 1
  100       CONTINUE
C
            J = J + KRONR(I)
            KRONR(I) = 0
  120    CONTINUE
C
         NKROR = J - 1
C
         DO 140 I = 1, NKROR
            KRONR(I) = IWORK(I)
  140    CONTINUE
C
      END IF
C
C     Set left Kronecker indices (row indices).
C
      IF ( NKROL.GT.0 ) THEN
         J = 1
C
         DO 180 I = 1, N + 1
C
            DO 160 II = J, J + KRONL(I) - 1
               IWORK(II) = I - 1
  160       CONTINUE
C
            J = J + KRONL(I)
            KRONL(I) = 0
  180    CONTINUE
C
         NKROL = J - 1
C
         DO 200 I = 1, NKROL
            KRONL(I) = IWORK(I)
  200    CONTINUE
C
      END IF
C
      IF ( N.GT.0 ) THEN
         DINFZ = N
C
  220    CONTINUE
         IF ( INFZ(DINFZ).EQ.0 ) THEN
            DINFZ = DINFZ - 1
            IF ( DINFZ.GT.0 )
     $         GO TO 220
         END IF
      END IF
C
      DWORK(1) = WRKOPT
      RETURN
C *** Last line of AB08ND ***
      END
