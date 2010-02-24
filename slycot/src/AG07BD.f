      SUBROUTINE AG07BD( JOBE, N, M, A, LDA, E, LDE, B, LDB, C, LDC,
     $                   D, LDD, AI, LDAI, EI, LDEI, BI, LDBI, CI, LDCI,
     $                   DI, LDDI, INFO )
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
C     To compute the inverse (Ai-lambda*Ei,Bi,Ci,Di) of a given
C     descriptor system (A-lambda*E,B,C,D).
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     JOBE    CHARACTER*1
C             Specifies whether E is a general square or an identity
C             matrix as follows:
C             = 'G':  E is a general square matrix;
C             = 'I':  E is the identity matrix.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the square matrices A and E;
C             also the number of rows of matrix B and the number of
C             columns of matrix C.  N >= 0.
C
C     M       (input) INTEGER
C             The number of system inputs and outputs, i.e., the number
C             of columns of matrices B and D and the number of rows of
C             matrices C and D.  M >= 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state matrix A of the original system.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= MAX(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,N)
C             If JOBE = 'G', the leading N-by-N part of this array must
C             contain the descriptor matrix E of the original system.
C             If JOBE = 'I', then E is assumed to be the identity
C             matrix and is not referenced.
C
C     LDE     INTEGER
C             The leading dimension of the array E.
C             LDE >= MAX(1,N), if JOBE = 'G';
C             LDE >= 1,        if JOBE = 'I'.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             input matrix B of the original system.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= MAX(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading M-by-N part of this array must contain the
C             output matrix C of the original system.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= MAX(1,M).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             The leading M-by-M part of this array must contain the
C             feedthrough matrix D of the original system.
C
C     LDD     INTEGER
C             The leading dimension of the array D.  LDD >= MAX(1,M).
C
C     AI      (output) DOUBLE PRECISION array, dimension (LDAI,N+M)
C             The leading (N+M)-by-(N+M) part of this array contains
C             the state matrix Ai of the inverse system.
C             If LDAI = LDA >= N+M, then AI and A can share the same
C             storage locations.
C
C     LDAI    INTEGER
C             The leading dimension of the array AI.
C             LDAI >= MAX(1,N+M).
C
C     EI      (output) DOUBLE PRECISION array, dimension (LDEI,N+M)
C             The leading (N+M)-by-(N+M) part of this array contains
C             the descriptor matrix Ei of the inverse system.
C             If LDEI = LDE >= N+M, then EI and E can share the same
C             storage locations.
C
C     LDEI    INTEGER
C             The leading dimension of the array EI.
C             LDEI >= MAX(1,N+M).
C
C     BI      (output) DOUBLE PRECISION array, dimension (LDBI,M)
C             The leading (N+M)-by-M part of this array contains
C             the input matrix Bi of the inverse system.
C             If LDBI = LDB >= N+M, then BI and B can share the same
C             storage locations.
C
C     LDBI    INTEGER
C             The leading dimension of the array BI.
C             LDBI >= MAX(1,N+M).
C
C     CI      (output) DOUBLE PRECISION array, dimension (LDCI,N+M)
C             The leading M-by-(N+M) part of this array contains
C             the output matrix Ci of the inverse system.
C             If LDCI = LDC, CI and C can share the same storage
C             locations.
C
C     LDCI    INTEGER
C             The leading dimension of the array CI.  LDCI >= MAX(1,M).
C
C     DI      (output) DOUBLE PRECISION array, dimension (LDDI,M)
C             The leading M-by-M part of this array contains
C             the feedthrough matrix Di = 0 of the inverse system.
C             DI and D can share the same storage locations.
C
C     LDDI    INTEGER
C             The leading dimension of the array DI.  LDDI >= MAX(1,M).
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
C     The matrices of the inverse system are computed with the formulas
C
C                ( E  0 )        ( A  B )         (  0 )
C           Ei = (      ) , Ai = (      ) ,  Bi = (    ),
C                ( 0  0 )        ( C  D )         ( -I )
C
C           Ci = ( 0  I ),  Di = 0.
C
C     FURTHER COMMENTS
C
C     The routine does not perform an invertibility test. This check can
C     be performed by using the SLICOT routines AB08NX or AG08BY.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001.
C
C     KEYWORDS
C
C     Descriptor system, inverse system, state-space representation.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER          JOBE
      INTEGER            INFO, LDA, LDAI, LDB, LDBI, LDC, LDCI,
     $                   LDD, LDDI, LDE, LDEI, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*), AI(LDAI,*), B(LDB,*), BI(LDBI,*),
     $                   C(LDC,*), CI(LDCI,*), D(LDD,*), DI(LDDI,*),
     $                   E(LDE,*), EI(LDEI,*)
C     .. Local Scalars ..
      LOGICAL            UNITE
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           DLACPY, DLASET, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      UNITE = LSAME( JOBE, 'I' )
      IF( .NOT. ( LSAME( JOBE, 'G' ) .OR. UNITE ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDE.LT.1 .OR. ( .NOT.UNITE .AND. LDE.LT.N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -13
      ELSE IF( LDAI.LT.MAX( 1, N+M ) ) THEN
         INFO = -15
      ELSE IF( LDEI.LT.MAX( 1, N+M ) ) THEN
         INFO = -17
      ELSE IF( LDBI.LT.MAX( 1, N+M ) ) THEN
         INFO = -19
      ELSE IF( LDCI.LT.MAX( 1, M ) ) THEN
         INFO = -21
      ELSE IF( LDDI.LT.MAX( 1, M ) ) THEN
         INFO = -23
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AG07BD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( M.EQ.0 )
     $   RETURN
C
C     Form Ai.
C
      CALL DLACPY( 'Full', N, N, A, LDA, AI, LDAI )
      CALL DLACPY( 'Full', M, N, C, LDC, AI(N+1,1), LDAI )
      CALL DLACPY( 'Full', N, M, B, LDB, AI(1,N+1), LDAI )
      CALL DLACPY( 'Full', M, M, D, LDD, AI(N+1,N+1), LDAI )
C
C     Form Ei.
C
      IF( UNITE ) THEN
         CALL DLASET( 'Full', N+M, N, ZERO, ONE, EI, LDEI )
      ELSE
         CALL DLACPY( 'Full', N, N, E, LDE, EI, LDEI )
         CALL DLASET( 'Full', M, N, ZERO, ZERO, EI(N+1,1), LDEI )
      END IF
      CALL DLASET( 'Full', N+M, M, ZERO, ZERO, EI(1,N+1), LDEI )
C
C     Form Bi.
C
      CALL DLASET( 'Full', N, M, ZERO, ZERO, BI, LDBI )
      CALL DLASET( 'Full', M, M, ZERO, -ONE, BI(N+1,1), LDBI )
C
C     Form Ci.
C
      CALL DLASET( 'Full', M, N, ZERO, ZERO, CI, LDCI )
      CALL DLASET( 'Full', M, M, ZERO, ONE, CI(1,N+1), LDCI )
C
C     Set Di.
C
      CALL DLASET( 'Full', M, M, ZERO, ZERO, DI, LDDI )
C
      RETURN
C *** Last line of AG07BD ***
      END
