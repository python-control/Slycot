      SUBROUTINE TG01HD( JOBCON, COMPQ, COMPZ, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, NCONT, NIUCON,
     $                   NRBLCK, RTAU, TOL, IWORK, DWORK, INFO )
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
C                ( Ac  *  )             ( Ec  *  )           ( Bc )
C       Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (    ) ,
C                ( 0  Anc )             ( 0  Enc )           ( 0  )
C
C          C*Z = ( Cc Cnc ) ,
C
C     where the NCONT-th order descriptor system (Ac-lambda*Ec,Bc,Cc)
C     is a finite and/or infinite controllable. The pencil
C     Anc - lambda*Enc is regular of order N-NCONT and contains the
C     uncontrollable finite and/or infinite eigenvalues of the pencil
C     A-lambda*E.
C
C     For JOBCON = 'C' or 'I', the pencil ( Bc Ec-lambda*Ac ) has full
C     row rank NCONT for all finite lambda and is in a staircase form
C     with
C                     _      _          _        _
C                   ( E1,0   E1,1  ...  E1,k-1   E1,k  )
C                   (        _          _        _     )
C       ( Bc Ec ) = (  0     E2,1  ...  E2,k-1   E2,k  ) ,  (1)
C                   (              ...  _        _     )
C                   (  0       0   ...  Ek,k-1   Ek,k  )
C
C                     _          _        _
C                   ( A1,1  ...  A1,k-1   A1,k  )
C                   (            _        _     )
C         Ac      = (   0   ...  A2,k-1   A2,k  ) ,         (2)
C                   (       ...           _     )
C                   (   0   ...    0      Ak,k  )
C           _
C     where Ei,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
C                            _
C     (with rtau(0) = M) and Ai,i is an rtau(i)-by-rtau(i)
C     upper triangular matrix.
C
C     For JOBCON = 'F', the pencil ( Bc Ac-lambda*Ec ) has full
C     row rank NCONT for all finite lambda and is in a staircase form
C     with
C                     _     _          _        _
C                   ( A1,0  A1,1  ...  A1,k-1   A1,k  )
C                   (       _          _        _     )
C       ( Bc Ac ) = (  0    A2,1  ...  A2,k-1   A2,k  ) ,   (3)
C                   (             ...  _        _     )
C                   (  0      0   ...  Ak,k-1   Ak,k  )
C
C                     _          _        _
C                   ( E1,1  ...  E1,k-1   E1,k  )
C                   (            _        _     )
C         Ec      = (   0   ...  E2,k-1   E2,k  ) ,         (4)
C                   (       ...           _     )
C                   (   0   ...    0      Ek,k  )
C           _
C     where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix
C                            _
C     (with rtau(0) = M) and Ei,i is an rtau(i)-by-rtau(i)
C     upper triangular matrix.
C
C     For JOBCON = 'C', the (N-NCONT)-by-(N-NCONT) regular pencil
C     Anc - lambda*Enc has the form
C
C                         ( Ainc - lambda*Einc         *          )
C      Anc - lambda*Enc = (                                       ) ,
C                         (        0           Afnc - lambda*Efnc )
C
C     where:
C       1) the NIUCON-by-NIUCON regular pencil Ainc - lambda*Einc,
C          with Ainc upper triangular and nonsingular, contains the
C          uncontrollable infinite eigenvalues of A - lambda*E;
C       2) the (N-NCONT-NIUCON)-by-(N-NCONT-NIUCON) regular pencil
C          Afnc - lambda*Efnc, with Efnc upper triangular and
C          nonsingular, contains the uncontrollable finite
C          eigenvalues of A - lambda*E.
C
C     Note: The significance of the two diagonal blocks can be
C           interchanged by calling the routine with the
C           arguments A and E interchanged. In this case,
C           Ainc - lambda*Einc contains the uncontrollable zero
C           eigenvalues of A - lambda*E, while Afnc - lambda*Efnc
C           contains the uncontrollable nonzero finite and infinite
C           eigenvalues of A - lambda*E.
C
C     For JOBCON = 'F', the pencil Anc - lambda*Enc has the form
C
C        Anc - lambda*Enc = Afnc - lambda*Efnc ,
C
C     where the regular pencil Afnc - lambda*Efnc, with Efnc
C     upper triangular and nonsingular, contains the uncontrollable
C     finite eigenvalues of A - lambda*E.
C
C     For JOBCON = 'I', the pencil Anc - lambda*Enc has the form
C
C        Anc - lambda*Enc = Ainc - lambda*Einc ,
C
C     where the regular pencil Ainc - lambda*Einc, with Ainc
C     upper triangular and nonsingular, contains the uncontrollable
C     nonzero finite and infinite eigenvalues of A - lambda*E.
C
C     The left and/or right orthogonal transformations Q and Z
C     performed to reduce the system matrices can be optionally
C     accumulated.
C
C     The reduced order descriptor system (Ac-lambda*Ec,Bc,Cc) has
C     the same transfer-function matrix as the original system
C     (A-lambda*E,B,C).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBCON  CHARACTER*1
C             = 'C':  separate both finite and infinite uncontrollable
C                     eigenvalues;
C             = 'F':  separate only finite uncontrollable eigenvalues:
C             = 'I':  separate only nonzero finite and infinite
C                     uncontrollable eigenvalues.
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
C                                ( Ac   *  )
C                       Q'*A*Z = (         ) ,
C                                ( 0   Anc )
C
C             where Ac is NCONT-by-NCONT and Anc is
C             (N-NCONT)-by-(N-NCONT).
C             If JOBCON = 'F', the matrix ( Bc Ac ) is in the
C             controllability staircase form (3).
C             If JOBCON = 'C' or 'I', the submatrix Ac is upper
C             triangular.
C             If JOBCON = 'C', the Anc matrix has the form
C
C                             ( Ainc   *  )
C                       Anc = (           ) ,
C                             (  0   Afnc )
C
C             where the NIUCON-by-NIUCON matrix Ainc is nonsingular and
C             upper triangular.
C             If JOBCON = 'I', Anc is nonsingular and upper triangular.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,N).
C
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N)
C             On entry, the leading N-by-N part of this array must
C             contain the N-by-N descriptor matrix E.
C             On exit, the leading N-by-N part of this array contains
C             the transformed descriptor matrix Q'*E*Z,
C
C                                ( Ec   *  )
C                       Q'*E*Z = (         ) ,
C                                ( 0   Enc )
C
C             where Ec is NCONT-by-NCONT and Enc is
C             (N-NCONT)-by-(N-NCONT).
C             If JOBCON = 'C' or 'I', the matrix ( Bc Ec ) is in the
C             controllability staircase form (1).
C             If JOBCON = 'F', the submatrix Ec is upper triangular.
C             If JOBCON = 'C', the Enc matrix has the form
C
C                             ( Einc   *  )
C                       Enc = (           ) ,
C                             (  0   Efnc )
C
C             where the NIUCON-by-NIUCON matrix Einc is nilpotent
C             and the (N-NCONT-NIUCON)-by-(N-NCONT-NIUCON) matrix Efnc
C             is nonsingular and upper triangular.
C             If JOBCON = 'F', Enc is nonsingular and upper triangular.
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,N).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M)
C             On entry, the leading N-by-M part of this array must
C             contain the N-by-M input matrix B.
C             On exit, the leading N-by-M part of this array contains
C             the transformed input matrix
C
C                              ( Bc )
C                       Q'*B = (    ) ,
C                              ( 0  )
C
C              where Bc is NCONT-by-M.
C              For JOBCON = 'C' or 'I', the matrix ( Bc Ec ) is in the
C              controllability staircase form (1).
C              For JOBCON = 'F', the matrix ( Bc Ac ) is in the
C              controllability staircase form (3).
C
C     LDB     INTEGER
C             The leading dimension of array B.  LDB >= MAX(1,N).
C
C     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C             On entry, the leading P-by-N part of this array must
C             contain the state/output matrix C.
C             On exit, the leading P-by-N part of this array contains
C             the transformed matrix C*Z.
C
C     LDC     INTEGER
C             The leading dimension of array C.  LDC >= MAX(1,P).
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
C     NCONT   (output) INTEGER
C             The order of the reduced matrices Ac and Ec, and the
C             number of rows of reduced matrix Bc; also the order of
C             the controllable part of the pair (A-lambda*E,B).
C
C     NIUCON  (output) INTEGER
C             For JOBCON = 'C', the order of the reduced matrices
C             Ainc and Einc; also the number of uncontrollable
C             infinite eigenvalues of the pencil A - lambda*E.
C             For JOBCON = 'F' or 'I', NIUCON has no significance
C             and is set to zero.
C
C     NRBLCK  (output) INTEGER
C             For JOBCON = 'C' or 'I', the number k, of full row rank
C                    _
C             blocks Ei,i in the staircase form of the pencil
C             (Bc Ec-lambda*Ac) (see (1) and (2)).
C             For JOBCON = 'F', the number k, of full row rank blocks
C             _
C             Ai,i in the staircase form of the pencil (Bc Ac-lambda*Ec)
C             (see (3) and (4)).
C
C     RTAU    (output) INTEGER array, dimension (N)
C             RTAU(i), for i = 1, ..., NRBLCK, is the row dimension of
C                                     _         _
C             the full row rank block Ei,i-1 or Ai,i-1 in the staircase
C             form (1) or (3) for JOBCON = 'C' or 'I', or
C             for JOBCON = 'F', respectively.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             The tolerance to be used in rank determinations when
C             transforming (A-lambda*E, B). If the user sets TOL > 0,
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
C     IWORK   INTEGER array, dimension (M)
C
C     DWORK   DOUBLE PRECISION array, dimension MAX(N,2*M)
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
C     If the system matrices A, E and B are badly scaled, it is
C     generally recommendable to scale them with the SLICOT routine
C     TG01AD, before calling TG01HD.
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
C
C     KEYWORDS
C
C     Controllability, minimal realization, orthogonal canonical form,
C     orthogonal transformation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOBCON
      INTEGER            INFO, LDA, LDB, LDC, LDE, LDQ, LDZ,
     $                   M, N, NCONT, NIUCON, NRBLCK, P
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      INTEGER            IWORK( * ), RTAU( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, *  ),
     $                   DWORK( * ),  E( LDE, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
C     .. Local Scalars ..
      CHARACTER          JOBQ, JOBZ
      LOGICAL            FINCON, ILQ, ILZ, INFCON
      INTEGER            ICOMPQ, ICOMPZ, LBA, NR
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           TG01HX, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C
C     .. Executable Statements ..
C
C     Decode JOBCON.
C
      IF( LSAME( JOBCON, 'C' ) ) THEN
         FINCON = .TRUE.
         INFCON = .TRUE.
      ELSE IF( LSAME( JOBCON, 'F' ) ) THEN
         FINCON = .TRUE.
         INFCON = .FALSE.
      ELSE IF( LSAME( JOBCON, 'I' ) ) THEN
         FINCON = .FALSE.
         INFCON = .TRUE.
      ELSE
         FINCON = .FALSE.
         INFCON = .FALSE.
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
      IF( .NOT.FINCON .AND. .NOT.INFCON ) THEN
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
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -14
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -16
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -18
      ELSE IF( TOL.GE.ONE ) THEN
         INFO = -23
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'TG01HD', -INFO )
         RETURN
      END IF
C
      JOBQ = COMPQ
      JOBZ = COMPZ
C
      IF( FINCON ) THEN
C
C        Perform finite controllability form reduction.
C
         CALL TG01HX( JOBQ, JOBZ, N, N, M, P, N, MAX( 0, N-1 ), A, LDA,
     $                E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, NR,
     $                NRBLCK, RTAU, TOL, IWORK, DWORK, INFO )
         IF( NRBLCK.GT.1 ) THEN
            LBA = RTAU(1) + RTAU(2) - 1
         ELSE IF( NRBLCK.EQ.1 ) THEN
            LBA = RTAU(1) - 1
         ELSE
            LBA = 0
         END IF
         IF( ILQ ) JOBQ = 'U'
         IF( ILZ ) JOBZ = 'U'
      ELSE
         NR = N
         LBA = MAX( 0, N-1 )
      END IF
C
      IF( INFCON ) THEN
C
C        Perform infinite controllability form reduction.
C
         CALL TG01HX( JOBQ, JOBZ, N, N, M, P, NR, LBA, E, LDE,
     $                A, LDA, B, LDB, C, LDC, Q, LDQ, Z, LDZ, NCONT,
     $                NRBLCK, RTAU, TOL, IWORK, DWORK, INFO )
         IF( FINCON ) THEN
            NIUCON = NR - NCONT
         ELSE
            NIUCON = 0
         END IF
      ELSE
         NCONT  = NR
         NIUCON = 0
      END IF
C
      RETURN
C
C *** Last line of TG01HD ***
      END
