      SUBROUTINE TG01JD( JOB, SYSTYP, EQUIL, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, NR, INFRED, TOL, IWORK, DWORK,
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
C     To find a reduced (controllable, observable, or irreducible)
C     descriptor representation (Ar-lambda*Er,Br,Cr) for an original
C     descriptor representation (A-lambda*E,B,C).
C     The pencil Ar-lambda*Er is in an upper block Hessenberg form, with
C     either Ar or Er upper triangular.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOB     CHARACTER*1
C             Indicates whether the user wishes to remove the
C             uncontrollable and/or unobservable parts as follows:
C             = 'I':  Remove both the uncontrollable and unobservable
C                     parts to get an irreducible descriptor
C                     representation;
C             = 'C':  Remove the uncontrollable part only to get a
C                     controllable descriptor representation;
C             = 'O':  Remove the unobservable part only to get an
C                     observable descriptor representation.
C
C     SYSTYP  CHARACTER*1
C             Indicates the type of descriptor system algorithm
C             to be applied according to the assumed
C             transfer-function matrix as follows:
C             = 'R':  Rational transfer-function matrix;
C             = 'S':  Proper (standard) transfer-function matrix;
C             = 'P':  Polynomial transfer-function matrix.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily scale
C             the system (A-lambda*E,B,C) as follows:
C             = 'S':  Perform scaling;
C             = 'N':  Do not perform scaling.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of the descriptor state vector; also the
C             order of square matrices A and E, the number of rows of
C             matrix B, and the number of columns of matrix C.  N >= 0.
C
C     M       (input) INTEGER
C             The dimension of descriptor system input vector; also the
C             number of columns of matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The dimension of descriptor system output vector; also the
C             number of rows of matrix C.  P >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original state matrix A.
C             On exit, the leading NR-by-NR part of this array contains
C             the reduced order state matrix Ar of an irreducible,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'I',
C             JOB = 'C', or JOB = 'O', respectively.
C             The matrix Ar is upper triangular if SYSTYP = 'R' or 'P'.
C             If SYSTYP = 'S' and JOB = 'C', the matrix [Br Ar]
C             is in a controllable staircase form (see TG01HD).
C             If SYSTYP = 'S' and JOB = 'I' or 'O', the matrix ( Ar )
C                                                              ( Cr )
C             is in an observable staircase form (see TG01HD).
C             The block structure of staircase forms is contained
C             in the leading INFRED(7) elements of IWORK.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the original descriptor matrix E.
C             On exit, the leading NR-by-NR part of this array contains
C             the reduced order descriptor matrix Er of an irreducible,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'I',
C             JOB = 'C', or JOB = 'O', respectively.
C             The resulting Er has INFRED(6) nonzero sub-diagonals.
C             If at least for one k = 1,...,4, INFRED(k) >= 0, then the
C             resulting Er is structured being either upper triangular
C             or block Hessenberg, in accordance to the last
C             performed order reduction phase (see METHOD).
C             The block structure of staircase forms is contained
C             in the leading INFRED(7) elements of IWORK.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M),
C             if JOB = 'C', or (LDB,MAX(M,P)), otherwise.
C             On entry, the leading N-by-M part of this array must
C             contain the original input matrix B; if JOB = 'I',
C             or JOB = 'O', the remainder of the leading N-by-MAX(M,P)
C             part is used as internal workspace.
C             On exit, the leading NR-by-M part of this array contains
C             the reduced input matrix Br of an irreducible,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'I',
C             JOB = 'C', or JOB = 'O', respectively.
C             If JOB = 'C', only the first IWORK(1) rows of B are
C             nonzero.
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the original output matrix C; if JOB = 'I',
C             or JOB = 'O', the remainder of the leading MAX(M,P)-by-N
C             part is used as internal workspace.
C             On exit, the leading P-by-NR part of this array contains
C             the transformed state/output matrix Cr of an irreducible,
C             controllable, or observable realization for the original
C             system, depending on the value of JOB, JOB = 'I',
C             JOB = 'C', or JOB = 'O', respectively.
C             If JOB = 'I', or JOB = 'O', only the last IWORK(1) columns
C             (in the first NR columns) of C are nonzero.
C
C     LDC     INTEGER
C             The leading dimension of array C.
C             LDC >= MAX(1,M,P) if N > 0.
C             LDC >= 1          if N = 0.
C
C     NR      (output) INTEGER
C             The order of the reduced descriptor representation
C             (Ar-lambda*Er,Br,Cr) of an irreducible, controllable,
C             or observable realization for the original system,
C             depending on JOB = 'I', JOB = 'C', or JOB = 'O',
C             respectively.
C
C     INFRED  (output) INTEGER array, dimension 7
C             This array contains information on performed reduction
C             and on structure of resulting system matrices as follows:
C             INFRED(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction
C                            (see METHOD) has been performed. In this
C                            case, INFRED(k) is the achieved order
C                            reduction in Phase k.
C             INFRED(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not
C                            performed.
C             INFRED(5)  -   the number of nonzero sub-diagonals of A.
C             INFRED(6)  -   the number of nonzero sub-diagonals of E.
C             INFRED(7)  -   the number of blocks in the resulting
C                            staircase form at last performed reduction
C                            phase. The block dimensions are contained
C                            in the first INFRED(7) elements of IWORK.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determinations when
C             transforming (A-lambda*E,B,C). If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             reciprocal condition numbers in rank determinations; a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*N*EPS,  is used instead, where
C             EPS is the machine precision (see LAPACK Library routine
C             DLAMCH).  TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension N+MAX(M,P)
C             On exit, if INFO = 0, the leading INFRED(7) elements of
C             IWORK contain the orders of the diagonal blocks of
C             Ar-lambda*Er.
C
C     DWORK   DOUBLE PRECISION array, dimension LDWORK
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(8*N,2*M,2*P), if EQUIL = 'S';
C             LDWORK >= MAX(N,2*M,2*P),   if EQUIL = 'N'.
C             If LDWORK >= MAX(2*N*N+N*M+N*P)+MAX(N,2*M,2*P) then more
C             accurate results are to be expected by performing only
C             those reductions phases (see METHOD), where effective
C             order reduction occurs. This is achieved by saving the
C             system matrices before each phase and restoring them if no
C             order reduction took place.
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
C     The subroutine is based on the reduction algorithms of [1].
C     The order reduction is performed in 4 phases:
C     Phase 1: Eliminate all finite uncontrolable eigenvalues.
C              The resulting matrix ( Br Ar ) is in a controllable
C              staircase form (see SLICOT Library routine TG01HD), and
C              Er is upper triangular.
C              This phase is performed if JOB = 'I' or 'C' and
C              SYSTYP = 'R' or 'S'.
C     Phase 2: Eliminate all infinite and finite nonzero uncontrollable
C              eigenvalues. The resulting matrix ( Br Er ) is in a
C              controllable staircase form (see TG01HD), and Ar is
C              upper triangular.
C              This phase is performed if JOB = 'I' or 'C' and
C              SYSTYP = 'R' or 'P'.
C     Phase 3: Eliminate all finite unobservable eigenvalues.
C              The resulting matrix ( Ar ) is in an observable
C                                   ( Cr )
C              staircase form (see SLICOT Library routine TG01ID), and
C              Er is upper triangular.
C              This phase is performed if JOB = 'I' or 'O' and
C              SYSTYP = 'R' or 'S'.
C     Phase 4: Eliminate all infinite and finite nonzero unobservable
C              eigenvalues. The resulting matrix ( Er ) is in an
C                                                ( Cr )
C              observable staircase form (see TG01ID), and Ar is
C              upper triangular.
C              This phase is performed if JOB = 'I' or 'O' and
C              SYSTYP = 'R' or 'P'.
C
C     REFERENCES
C
C     [1] A. Varga
C         Computation of Irreducible Generalized State-Space
C         Realizations.
C         Kybernetika, vol. 26, pp. 89-106, 1990.
C
C     NUMERICAL ASPECTS
C
C     The algorithm is numerically backward stable and requires
C     0( N**3 )  floating point operations.
C
C     FURTHER COMMENTS
C
C     If the pencil (A-lambda*E) has no zero eigenvalues, then an
C     irreducible realization can be computed skipping Phases 1 and 3
C     by using the setting: JOB = 'I' and SYSTYP = 'P'.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     April 1999. Based on the RASP routine RPDSIR.
C
C     REVISIONS
C
C     July 1999, V. Sima, Research Institute for Informatics, Bucharest.
C     May 2003, A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     May 2003, March 2004, V. Sima.
C
C     KEYWORDS
C
C     Controllability, irreducible realization, observability,
C     orthogonal canonical form, orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER         EQUIL, JOB, SYSTYP
      INTEGER           INFO, LDA, LDB, LDC, LDE, LDWORK, M, N, NR, P
      DOUBLE PRECISION  TOL
C     .. Array Arguments ..
      INTEGER           INFRED(*), IWORK(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), DWORK(*), E(LDE,*)
C     .. Local Scalars ..
      CHARACTER         JOBQ, JOBZ
      LOGICAL           FINCON, FINOBS, INFCON, INFOBS, LEQUIL, LJOBC,
     $                  LJOBIR, LJOBO, LSPACE, LSYSP, LSYSR, LSYSS
      INTEGER           KWA, KWB, KWC, KWE, LBA, LBE, LDM, LDP, LDQ,
     $                  LDZ, M1, MAXMP, N1, NBLCK, NC, P1
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
C     .. External Subroutines ..
      EXTERNAL          DLACPY, MA02CD, TB01XD, TG01AD, TG01HX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO  = 0
      MAXMP = MAX( M, P )
      N1    = MAX( 1, N )
C
C     Decode JOB.
C
      LJOBIR = LSAME( JOB, 'I' )
      LJOBC  = LJOBIR .OR. LSAME( JOB, 'C' )
      LJOBO  = LJOBIR .OR. LSAME( JOB, 'O' )
C
C     Decode SYSTYP.
C
      LSYSR  = LSAME( SYSTYP, 'R' )
      LSYSS  = LSYSR .OR. LSAME( SYSTYP, 'S' )
      LSYSP  = LSYSR .OR. LSAME( SYSTYP, 'P' )
C
      LEQUIL =  LSAME( EQUIL, 'S' )
C
C     Test the input scalar arguments.
C
      IF( .NOT.LJOBC .AND. .NOT.LJOBO ) THEN
         INFO = -1
      ELSE IF( .NOT.LSYSS .AND. .NOT.LSYSP ) THEN
         INFO = -2
      ELSE IF( .NOT.LEQUIL .AND. .NOT.LSAME( EQUIL, 'N' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.N1 ) THEN
         INFO = -8
      ELSE IF( LDE.LT.N1 ) THEN
         INFO = -10
      ELSE IF( LDB.LT.N1 ) THEN
         INFO = -12
      ELSE IF( LDC.LT.1 .OR. ( N.GT.0 .AND. LDC.LT.MAXMP ) ) THEN
         INFO = -14
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -17
      ELSE IF( ( .NOT.LEQUIL .AND. LDWORK.LT.MAX( N, 2*MAXMP ) ) .OR.
     $         ( LEQUIL .AND. LDWORK.LT.MAX( 8*N, 2*MAXMP ) ) ) THEN
         INFO = -20
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'TG01JD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      INFRED(1) = -1
      INFRED(2) = -1
      INFRED(3) = -1
      INFRED(4) = -1
      INFRED(5) =  0
      INFRED(6) =  0
      INFRED(7) =  0
C
      IF( MAX( N, MAXMP ).EQ.0 ) THEN
         NR = 0
         RETURN
      END IF
C
      M1  = MAX( 1, M )
      P1  = MAX( 1, P )
      LDM = MAX( LDC, M )
      LDP = MAX( LDC, P )
C
C     Set controllability/observability determination options.
C
      FINCON = LJOBC .AND. LSYSS
      INFCON = LJOBC .AND. LSYSP
      FINOBS = LJOBO .AND. LSYSS
      INFOBS = LJOBO .AND. LSYSP
C
C     Set large workspace option and determine offsets.
C
      LSPACE = LDWORK.GE.N*( 2*N + M + P ) + MAX( N, 2*MAXMP )
      KWA = MAX( N, 2*MAXMP ) + 1
      KWE = KWA + N*N
      KWB = KWE + N*N
      KWC = KWB + N*M
C
C     If required, scale the system (A-lambda*E,B,C).
C     Workspace: need 8*N.
C
      IF( LEQUIL ) THEN
         CALL TG01AD( 'All', N, N, M, P, ZERO, A, LDA, E, LDE, B, LDB,
     $                C, LDP, DWORK(1), DWORK(N+1), DWORK(2*N+1), INFO )
      END IF
C
      JOBQ = 'N'
      JOBZ = 'N'
      LDQ = 1
      LDZ = 1
      LBA = MAX( 0, N-1 )
      LBE = LBA
      NC = N
      NR = N
C
      IF( FINCON ) THEN
C
C        Phase 1: Eliminate all finite uncontrolable eigenvalues.
C
         IF( LSPACE) THEN
C
C           Save system matrices.
C
            CALL DLACPY( 'Full', NC, NC, A, LDA, DWORK(KWA), N1 )
            CALL DLACPY( 'Full', NC, NC, E, LDE, DWORK(KWE), N1 )
            CALL DLACPY( 'Full', NC, M,  B, LDB, DWORK(KWB), N1 )
            CALL DLACPY( 'Full', P,  NC, C, LDC, DWORK(KWC), P1 )
         END IF
C
C        Perform finite controllability form reduction.
C        Workspace: need   MAX(N,2*M).
C
         CALL TG01HX( JOBQ, JOBZ, NC, NC, M, P, NC, LBE, A, LDA,
     $                E, LDE, B, LDB, C, LDP, DUM, LDQ, DUM, LDZ, NR,
     $                NBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         IF( NR.LT.NC .OR. .NOT.LSPACE ) THEN
            IF( NBLCK.GT.1 ) THEN
               LBA = IWORK(1) + IWORK(2) - 1
            ELSE IF( NBLCK.EQ.1 ) THEN
               LBA = IWORK(1) - 1
            ELSE
               LBA = 0
            END IF
            LBE = 0
            INFRED(1) = NC - NR
            INFRED(7) = NBLCK
            NC = NR
         ELSE
C
C           Restore system matrices.
C
            CALL DLACPY( 'Full', NC, NC, DWORK(KWA), N1, A, LDA )
            CALL DLACPY( 'Full', NC, NC, DWORK(KWE), N1, E, LDE )
            CALL DLACPY( 'Full', NC, M,  DWORK(KWB), N1, B, LDB )
            CALL DLACPY( 'Full', P,  NC, DWORK(KWC), P1, C, LDC )
         END IF
      END IF
C
      IF( INFCON ) THEN
C
C        Phase 2: Eliminate all infinite and all finite nonzero
C                 uncontrolable eigenvalues.
C
         IF( LSPACE ) THEN
C
C           Save system matrices.
C
            CALL DLACPY( 'Full', NC, NC, A, LDA, DWORK(KWA), N1 )
            CALL DLACPY( 'Full', NC, NC, E, LDE, DWORK(KWE), N1 )
            CALL DLACPY( 'Full', NC, M,  B, LDB, DWORK(KWB), N1 )
            CALL DLACPY( 'Full', P,  NC, C, LDC, DWORK(KWC), P1 )
         END IF
C
C        Perform infinite controllability form reduction.
C        Workspace: need   MAX(N,2*M).
C
         CALL TG01HX( JOBQ, JOBZ, NC, NC, M, P, NC, LBA, E, LDE,
     $                A, LDA, B, LDB, C, LDP, DUM, LDQ, DUM, LDZ, NR,
     $                NBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         IF( NR.LT.NC .OR. .NOT.LSPACE ) THEN
            IF( NBLCK.GT.1 ) THEN
               LBE = IWORK(1) + IWORK(2) - 1
            ELSE IF( NBLCK.EQ.1 ) THEN
               LBE = IWORK(1) - 1
            ELSE
               LBE = 0
            END IF
            LBA = 0
            INFRED(2) = NC - NR
            INFRED(7) = NBLCK
            NC = NR
         ELSE
C
C           Restore system matrices.
C
            CALL DLACPY( 'Full', NC, NC, DWORK(KWA), N1, A, LDA )
            CALL DLACPY( 'Full', NC, NC, DWORK(KWE), N1, E, LDE )
            CALL DLACPY( 'Full', NC, M,  DWORK(KWB), N1, B, LDB )
            CALL DLACPY( 'Full', P,  NC, DWORK(KWC), P1, C, LDC )
         END IF
      END IF
C
      IF( FINOBS .OR. INFOBS) THEN
C
C        Compute the pertransposed dual system exploiting matrix shapes.
C
         CALL TB01XD( 'Z', NC, M, P, LBA, MAX( 0, NC-1 ), A, LDA,
     $                B, LDB, C, LDC, DUM, 1, INFO )
         CALL MA02CD( NC, LBE, MAX( 0, NC-1 ), E, LDE )
      END IF
C
      IF( FINOBS ) THEN
C
C        Phase 3: Eliminate all finite unobservable eigenvalues.
C
         IF( LSPACE ) THEN
C
C           Save system matrices.
C
            CALL DLACPY( 'Full', NC, NC, A, LDA, DWORK(KWA), N1 )
            CALL DLACPY( 'Full', NC, NC, E, LDE, DWORK(KWE), N1 )
            CALL DLACPY( 'Full', NC, P,  B, LDB, DWORK(KWC), N1 )
            CALL DLACPY( 'Full', M,  NC, C, LDC, DWORK(KWB), M1 )
         END IF
C
C        Perform finite observability form reduction.
C        Workspace: need   MAX(N,2*P).
C
         CALL TG01HX( JOBZ, JOBQ, NC, NC, P, M, NC, LBE, A, LDA,
     $                E, LDE, B, LDB, C, LDM, DUM, LDZ, DUM, LDQ, NR,
     $                NBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         IF( NR.LT.NC .OR. .NOT.LSPACE ) THEN
            IF( NBLCK.GT.1 ) THEN
               LBA = IWORK(1) + IWORK(2) - 1
            ELSE IF( NBLCK.EQ.1 ) THEN
               LBA = IWORK(1) - 1
            ELSE
               LBA = 0
            END IF
            LBE = 0
            INFRED(3) = NC - NR
            INFRED(7) = NBLCK
            NC = NR
         ELSE
C
C           Restore system matrices.
C
            CALL DLACPY( 'Full', NC, NC, DWORK(KWA), N1, A, LDA )
            CALL DLACPY( 'Full', NC, NC, DWORK(KWE), N1, E, LDE )
            CALL DLACPY( 'Full', NC, P,  DWORK(KWC), N1, B, LDB )
            CALL DLACPY( 'Full', M,  NC, DWORK(KWB), M1, C, LDC )
         END IF
      END IF
C
      IF( INFOBS ) THEN
C
C        Phase 4: Eliminate all infinite and all finite nonzero
C                 unobservable eigenvalues.
C
         IF( LSPACE) THEN
C
C           Save system matrices.
C
            CALL DLACPY( 'Full', NC, NC, A, LDA, DWORK(KWA), N1 )
            CALL DLACPY( 'Full', NC, NC, E, LDE, DWORK(KWE), N1 )
            CALL DLACPY( 'Full', NC, P,  B, LDB, DWORK(KWC), N1 )
            CALL DLACPY( 'Full', M,  NC, C, LDC, DWORK(KWB), M1 )
         END IF
C
C        Perform infinite observability form reduction.
C        Workspace: need   MAX(N,2*P).
C
         CALL TG01HX( JOBZ, JOBQ, NC, NC, P, M, NC, LBA, E, LDE,
     $                A, LDA, B, LDB, C, LDM, DUM, LDZ, DUM, LDQ, NR,
     $                NBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         IF( NR.LT.NC .OR. .NOT.LSPACE ) THEN
            IF( NBLCK.GT.1 ) THEN
               LBE = IWORK(1) + IWORK(2) - 1
            ELSE IF( NBLCK.EQ.1 ) THEN
               LBE = IWORK(1) - 1
            ELSE
               LBE = 0
            END IF
            LBA = 0
            INFRED(4) = NC - NR
            INFRED(7) = NBLCK
            NC = NR
         ELSE
C
C           Restore system matrices.
C
            CALL DLACPY( 'Full', NC, NC, DWORK(KWA), N1, A, LDA )
            CALL DLACPY( 'Full', NC, NC, DWORK(KWE), N1, E, LDE )
            CALL DLACPY( 'Full', NC, P,  DWORK(KWC), N1, B, LDB )
            CALL DLACPY( 'Full', M,  NC, DWORK(KWB), M1, C, LDC )
         END IF
      END IF
C
      IF( FINOBS .OR. INFOBS ) THEN
C
C        Compute the pertransposed dual system exploiting matrix shapes.
C
         CALL TB01XD( 'Z', NC, P, M, LBA, MAX( 0, NC-1 ), A, LDA,
     $                B, LDB, C, LDC, DUM, 1, INFO )
         CALL MA02CD( NC, LBE, MAX( 0, NC-1 ), E, LDE )
      END IF
C
C     Set structural information on A and E.
C
      INFRED(5) = LBA
      INFRED(6) = LBE
C
      RETURN
C *** Last line of TG01JD ***
      END
