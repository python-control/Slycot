      SUBROUTINE TG01ID( JOBOBS, COMPQ, COMPZ, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, NOBSV, NIUOBS,
     $                   NLBLCK, CTAU, TOL, IWORK, DWORK, INFO )
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
C     To compute orthogonal transformation matrices Q and Z which
C     reduce the N-th order descriptor system (A-lambda*E,B,C)
C     to the form
C
C                ( Ano  * )             ( Eno  * )           ( Bno )
C       Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (     ) ,
C                ( 0   Ao )             ( 0   Eo )           ( Bo  )
C
C          C*Z = ( 0   Co ) ,
C
C     where the NOBSV-th order descriptor system (Ao-lambda*Eo,Bo,Co)
C     is a finite and/or infinite observable. The pencil
C     Ano - lambda*Eno is regular of order N-NOBSV and contains the
C     unobservable finite and/or infinite eigenvalues of the pencil
C     A-lambda*E.
C
C     For JOBOBS = 'O' or 'I', the pencil ( Eo-lambda*Ao ) has full
C                                         (      Co      )
C     column rank NOBSV for all finite lambda and is in a staircase form
C     with
C                     _      _            _      _
C                   ( Ek,k   Ek,k-1   ... Ek,2   Ek,1   )
C                   ( _      _            _      _      )
C       ( Eo ) =    ( Ek-1,k Ek-1,k-1 ... Ek-1,2 Ek-1,1 ) ,  (1)
C       ( Co )      (     ...         ... _      _      )
C                   (  0       0      ... E1,2   E1,1   )
C                   (                            _      )
C                   (  0       0      ... 0      E0,1   )
C                     _          _      _
C                   ( Ak,k  ...  Ak,2   Ak,1 )
C                   (       ...  _      _    )
C         Ao      = (   0   ...  A2,2   A2,1 ) ,             (2)
C                   (                   _    )
C                   (   0   ...    0    A1,1 )
C           _
C     where Ei-1,i is a CTAU(i-1)-by-CTAU(i) full column rank matrix
C                            _
C     (with CTAU(0) = P) and Ai,i is a CTAU(i)-by-CTAU(i)
C     upper triangular matrix.
C
C     For JOBOBS = 'F', the pencil ( Ao-lambda*Eo ) has full
C                                  (      Co      )
C     column rank NOBSV for all finite lambda and is in a staircase form
C     with
C                     _      _            _      _
C                   ( Ak,k   Ak,k-1   ... Ak,2   Ak,1   )
C                   ( _      _            _      _      )
C       ( Ao ) =    ( Ak-1,k Ak-1,k-1 ... Ak-1,2 Ak-1,1 ) ,  (3)
C       ( Co )      (     ...         ... _      _      )
C                   (  0       0      ... A1,2   A1,1   )
C                   (                            _      )
C                   (  0       0      ... 0      A0,1   )
C                     _          _      _
C                   ( Ek,k  ...  Ek,2   Ek,1 )
C                   (       ...  _      _    )
C         Eo      = (   0   ...  E2,2   E2,1 ) ,             (4)
C                   (                   _    )
C                   (   0   ...    0    E1,1 )
C           _
C     where Ai-1,i is a CTAU(i-1)-by-CTAU(i) full column rank matrix
C                            _
C     (with CTAU(0) = P) and Ei,i is a CTAU(i)-by-CTAU(i)
C     upper triangular matrix.
C
C     For JOBOBS = 'O', the (N-NOBSV)-by-(N-NOBSV) regular pencil
C     Ano - lambda*Eno has the form
C
C                         ( Afno - lambda*Efno         *          )
C      Ano - lambda*Eno = (                                       ) ,
C                         (        0           Aino - lambda*Eino )
C
C     where:
C       1) the NIUOBS-by-NIUOBS regular pencil Aino - lambda*Eino,
C          with Aino upper triangular and nonsingular, contains the
C          unobservable infinite eigenvalues of A - lambda*E;
C       2) the (N-NOBSV-NIUOBS)-by-(N-NOBSV-NIUOBS) regular pencil
C          Afno - lambda*Efno, with Efno upper triangular and
C          nonsingular, contains the unobservable finite
C          eigenvalues of A - lambda*E.
C
C     Note: The significance of the two diagonal blocks can be
C           interchanged by calling the routine with the
C           arguments A and E interchanged. In this case,
C           Aino - lambda*Eino contains the unobservable zero
C           eigenvalues of A - lambda*E, while Afno - lambda*Efno
C           contains the unobservable nonzero finite and infinite
C           eigenvalues of A - lambda*E.
C
C     For JOBOBS = 'F', the pencil Ano - lambda*Eno has the form
C
C        Ano - lambda*Eno = Afno - lambda*Efno ,
C
C     where the regular pencil Afno - lambda*Efno, with Efno
C     upper triangular and nonsingular, contains the unobservable
C     finite eigenvalues of A - lambda*E.
C
C     For JOBOBS = 'I', the pencil Ano - lambda*Eno has the form
C
C        Ano - lambda*Eno = Aino - lambda*Eino ,
C
C     where the regular pencil Aino - lambda*Eino, with Aino
C     upper triangular and nonsingular, contains the unobservable
C     nonzero finite and infinite eigenvalues of A - lambda*E.
C
C     The left and/or right orthogonal transformations Q and Z
C     performed to reduce the system matrices can be optionally
C     accumulated.
C
C     The reduced order descriptor system (Ao-lambda*Eo,Bo,Co) has
C     the same transfer-function matrix as the original system
C     (A-lambda*E,B,C).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBOBS   CHARACTER*1
C             = 'O':  separate both finite and infinite unobservable
C                     eigenvalues;
C             = 'F':  separate only finite unobservable eigenvalues;
C             = 'I':  separate only nonzero finite and infinite
C                     unobservable eigenvalues.
C
C     COMPQ   CHARACTER*1
C             = 'N':  do not compute Q;
C             = 'I':  Q is initialized to the unit matrix, and the
C                     orthogonal matrix Q is returned;
C             = 'U':  Q must contain an orthogonal matrix Q1 on entry,
C                     and the product Q1*Q is returned.
C
C     COMPZ   CHARACTER*1
C             = 'N':  do not compute Z;
C             = 'I':  Z is initialized to the unit matrix, and the
C                     orthogonal matrix Z is returned;
C             = 'U':  Z must contain an orthogonal matrix Z1 on entry,
C                     and the product Z1*Z is returned.
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
C             contain the N-by-N state matrix A.
C             On exit, the leading N-by-N part of this array contains
C             the transformed state matrix Q'*A*Z,
C
C                                ( Ano  *  )
C                       Q'*A*Z = (         ) ,
C                                ( 0    Ao )
C
C             where Ao is NOBSV-by-NOBSV and Ano is
C             (N-NOBSV)-by-(N-NOBSV).
C             If JOBOBS = 'F', the matrix ( Ao ) is in the observability
C                                         ( Co )
C             staircase form (3).
C             If JOBOBS = 'O' or 'I', the submatrix Ao is upper
C             triangular.
C             If JOBOBS = 'O', the submatrix Ano has the form
C
C                             ( Afno   *  )
C                       Ano = (           ) ,
C                             (  0   Aino )
C
C             where the NIUOBS-by-NIUOBS matrix Aino is nonsingular and
C             upper triangular.
C             If JOBOBS = 'I', Ano is nonsingular and upper triangular.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the N-by-N descriptor matrix E.
C             On exit, the leading N-by-N part of this array contains
C             the transformed state matrix Q'*E*Z,
C
C                                ( Eno  *  )
C                       Q'*E*Z = (         ) ,
C                                ( 0    Eo )
C
C             where Eo is NOBSV-by-NOBSV and Eno is
C             (N-NOBSV)-by-(N-NOBSV).
C             If JOBOBS = 'O' or 'I', the matrix ( Eo ) is in the
C                                                ( Co )
C             observability staircase form (1).
C             If JOBOBS = 'F', the submatrix Eo is upper triangular.
C             If JOBOBS = 'O', the Eno matrix has the form
C
C                             ( Efno   *  )
C                       Eno = (           ) ,
C                             (  0   Eino )
C
C             where the NIUOBS-by-NIUOBS matrix Eino is nilpotent
C             and the (N-NOBSV-NIUOBS)-by-(N-NOBSV-NIUOBS) matrix Efno
C             is nonsingular and upper triangular.
C             If JOBOBS = 'F', Eno is nonsingular and upper triangular.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension
C             (LDB,MAX(M,P))
C             On entry, the leading N-by-M part of this array must
C             contain the N-by-M input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed input matrix Q'*B.
C
C     LDB     INTEGER
C             The leading dimension of array B.
C             LDB >= MAX(1,N) if M > 0 or LDB >= 1 if M = 0.
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed matrix
C
C                     C*Z = (  0   Co ) ,
C
C             where Co is P-by-NOBSV.
C             If JOBOBS = 'O' or 'I', the matrix ( Eo ) is in the
C                                                ( Co )
C             observability staircase form (1).
C             If JOBOBS = 'F', the matrix ( Ao ) is in the observability
C                                         ( Co )
C             staircase form (3).
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,M,P).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             If COMPQ = 'N': Q is not referenced.
C             If COMPQ = 'I': on entry, Q need not be set;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix Q,
C                             where Q' is the product of transformations
C                             which are applied to A, E, and B on
C                             the left.
C             If COMPQ = 'U': on entry, the leading N-by-N part of this
C                             array must contain an orthogonal matrix
C                             Qc;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix
C                             Qc*Q.
C
C     LDQ     INTEGER
C             The leading dimension of array Q.
C             LDQ >= 1,        if COMPQ = 'N';
C             LDQ >= MAX(1,N), if COMPQ = 'U' or 'I'.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             If COMPZ = 'N': Z is not referenced.
C             If COMPZ = 'I': on entry, Z need not be set;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix Z,
C                             which is the product of transformations
C                             applied to A, E, and C on the right.
C             If COMPZ = 'U': on entry, the leading N-by-N part of this
C                             array must contain an orthogonal matrix
C                             Zc;
C                             on exit, the leading N-by-N part of this
C                             array contains the orthogonal matrix
C                             Zc*Z.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.
C             LDZ >= 1,        if COMPZ = 'N';
C             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'.
C
C     NOBSV   (output) INTEGER
C             The order of the reduced matrices Ao and Eo, and the
C             number of columns of reduced matrix Co; also the order of
C             observable part of the pair (C, A-lambda*E).
C
C     NIUOBS  (output) INTEGER
C             For JOBOBS = 'O', the order of the reduced matrices
C             Aino and Eino; also the number of unobservable
C             infinite eigenvalues of the pencil A - lambda*E.
C             For JOBOBS = 'F' or 'I', NIUOBS has no significance
C             and is set to zero.
C
C     NLBLCK  (output) INTEGER
C             For JOBOBS = 'O' or 'I', the number k, of full column rank
C                    _
C             blocks Ei-1,i in the staircase form of the pencil
C             (Eo-lambda*Ao) (see (1) and (2)).
C             (    Co      )
C             For JOBOBS = 'F', the number k, of full column rank blocks
C             _
C             Ai-1,i in the staircase form of the pencil (Ao-lambda*Eo)
C                                                        (     Co     )
C             (see (3) and (4)).
C
C     CTAU    (output) INTEGER array, dimension (N)
C             CTAU(i), for i = 1, ..., NLBLCK, is the column dimension
C                                           _         _
C             of the full column rank block Ei-1,i or Ai-1,i in the
C             staircase form (1) or (3) for JOBOBS = 'O' or 'I', or
C             for JOBOBS = 'F', respectively.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determinations when
C             transforming (A'-lambda*E',C')'. If the user sets TOL > 0,
C             then the given value of TOL is used as a lower bound for
C             reciprocal condition numbers in rank determinations; a
C             (sub)matrix whose estimated condition number is less than
C             1/TOL is considered to be of full rank.  If the user sets
C             TOL <= 0, then an implicitly computed, default tolerance,
C             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS
C             is the machine precision (see LAPACK Library routine
C             DLAMCH).  TOL < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (P)
C
C     DWORK   DOUBLE PRECISION array, dimension MAX(N,2*P)
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
C     The subroutine is based on the dual of the reduction
C     algorithms of [1].
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
C     If the system matrices A, E and C are badly scaled, it is
C     generally recommendable to scale them with the SLICOT routine
C     TG01AD, before calling TG01ID.
C
C     CONTRIBUTOR
C
C     C. Oara, University "Politehnica" Bucharest.
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen.
C     March 1999. Based on the RASP routine RPDSCF.
C
C     REVISIONS
C
C     July 1999, V. Sima, Research Institute for Informatics, Bucharest.
C     May 2003, March 2004, V. Sima.
C
C     KEYWORDS
C
C     Observability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOBOBS
      INTEGER            INFO, LDA, LDB, LDC, LDE, LDQ, LDZ,
     $                   M, N, NIUOBS, NLBLCK, NOBSV, P
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      INTEGER            CTAU( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, *  ),
     $                   DWORK( * ), E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      CHARACTER          JOBQ, JOBZ
      LOGICAL            FINOBS, ILQ, ILZ, INFOBS
      INTEGER            I, ICOMPQ, ICOMPZ, LBA, LBE, NR
C     .. Local Arrays ..
      DOUBLE PRECISION   DUM(1)
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           AB07MD, DSWAP, MA02BD, MA02CD, TB01XD,
     $                   TG01HX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C
C     .. Executable Statements ..
C
C     Decode JOBOBS.
C
      IF( LSAME( JOBOBS, 'O') ) THEN
         FINOBS = .TRUE.
         INFOBS = .TRUE.
      ELSE IF( LSAME( JOBOBS, 'F') ) THEN
         FINOBS = .TRUE.
         INFOBS = .FALSE.
      ELSE IF( LSAME( JOBOBS, 'I') ) THEN
         FINOBS = .FALSE.
         INFOBS = .TRUE.
      ELSE
         FINOBS = .FALSE.
         INFOBS = .FALSE.
      END IF
C
C     Decode COMPQ.
C
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'U' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
C
C     Decode COMPZ.
C
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'U' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
C
C     Test the input scalar parameters.
C
      INFO = 0
      IF( .NOT.FINOBS .AND. .NOT.INFOBS ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.LE.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.1 .OR. ( M.GT.0 .AND. LDB.LT.N ) ) THEN
         INFO = -12
      ELSE IF( LDC.LT.MAX( 1, M, P ) ) THEN
         INFO = -14
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -16
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -18
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -23
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01ID', -INFO )
         RETURN
      END IF
C
      JOBQ = COMPQ
      JOBZ = COMPZ
C
C     Build the dual system.
C
      CALL AB07MD( 'Z', N, M, P, A, LDA, B, LDB, C, LDC, DUM, 1,
     $             INFO )
      DO 10 I = 2, N
         CALL DSWAP( I-1, E(I,1), LDE, E(1,I), 1 )
   10 CONTINUE
C
      IF( FINOBS ) THEN
C
C        Perform finite observability form reduction.
C
         CALL TG01HX( JOBZ, JOBQ, N, N, P, M, N, MAX( 0, N-1 ), A, LDA,
     $                E, LDE, B, LDB, C, LDC, Z, LDZ, Q, LDQ, NR,
     $                NLBLCK, CTAU, TOL, IWORK, DWORK, INFO )
         IF( NLBLCK.GT.1 ) THEN
            LBA = CTAU(1) + CTAU(2) - 1
         ELSE IF( NLBLCK.EQ.1 ) THEN
            LBA = CTAU(1) - 1
         ELSE
            LBA = 0
         END IF
         IF( ILQ ) JOBQ = 'U'
         IF( ILZ ) JOBZ = 'U'
         LBE = 0
      ELSE
         NR = N
         LBA = MAX( 0, N-1 )
         LBE = LBA
      END IF
C
      IF( INFOBS ) THEN
C
C        Perform infinite observability form reduction.
C
         CALL TG01HX( JOBZ, JOBQ, N, N, P, M, NR, LBA, E, LDE,
     $                A, LDA, B, LDB, C, LDC, Z, LDZ, Q, LDQ, NOBSV,
     $                NLBLCK, CTAU, TOL, IWORK, DWORK, INFO )
         IF( FINOBS ) THEN
            NIUOBS = NR - NOBSV
         ELSE
            NIUOBS = 0
         END IF
         IF( NLBLCK.GT.1 ) THEN
            LBE = CTAU(1) + CTAU(2) - 1
         ELSE IF( NLBLCK.EQ.1 ) THEN
            LBE = CTAU(1) - 1
         ELSE
            LBE = 0
         END IF
         LBA = 0
      ELSE
         NOBSV = NR
         NIUOBS = 0
      END IF
C
C     Compute the pertransposed dual system exploiting matrix shapes.
C
      LBA = MAX( LBA, NIUOBS-1, N-NOBSV-NIUOBS-1 )
      IF ( P.EQ.0 .OR. NR.EQ.0 )
     $   LBE = MAX( 0, N - 1 )
      CALL TB01XD( 'Z', N, P, M, LBA, MAX( 0, N-1 ), A, LDA, B, LDB,
     $             C, LDC, DUM, 1, INFO )
      CALL MA02CD( N, LBE, MAX( 0, N-1 ), E, LDE )
      IF( ILZ ) CALL MA02BD( 'Right', N, N, Z, LDZ )
      IF( ILQ ) CALL MA02BD( 'Right', N, N, Q, LDQ )
      RETURN
C *** Last line of TG01ID ***
      END
