      SUBROUTINE MB03QD( DICO, STDOM, JOBU, N, NLOW, NSUP, ALPHA,
     $                   A, LDA, U, LDU, NDIM, DWORK, INFO )
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
C     To reorder the diagonal blocks of a principal submatrix of an
C     upper quasi-triangular matrix A together with their eigenvalues by
C     constructing an orthogonal similarity transformation UT.
C     After reordering, the leading block of the selected submatrix of A
C     has eigenvalues in a suitably defined domain of interest, usually
C     related to stability/instability in a continuous- or discrete-time
C     sense.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the spectrum separation to be
C             performed as follows:
C             = 'C':  continuous-time sense;
C             = 'D':  discrete-time sense.
C
C     STDOM   CHARACTER*1
C             Specifies whether the domain of interest is of stability
C             type (left part of complex plane or inside of a circle)
C             or of instability type (right part of complex plane or
C             outside of a circle) as follows:
C             = 'S':  stability type domain;
C             = 'U':  instability type domain.
C
C     JOBU    CHARACTER*1
C             Indicates how the performed orthogonal transformations UT
C             are accumulated, as follows:
C             = 'I':  U is initialized to the unit matrix and the matrix
C                     UT is returned in U;
C             = 'U':  the given matrix U is updated and the matrix U*UT
C                     is returned in U.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices A and U.  N >= 1.
C
C     NLOW,   (input) INTEGER
C     NSUP    NLOW and NSUP specify the boundary indices for the rows
C             and columns of the principal submatrix of A whose diagonal
C             blocks are to be reordered.  1 <= NLOW <= NSUP <= N.
C
C     ALPHA   (input) DOUBLE PRECISION
C             The boundary of the domain of interest for the eigenvalues
C             of A. If DICO = 'C', ALPHA is the boundary value for the
C             real parts of eigenvalues, while for DICO = 'D',
C             ALPHA >= 0 represents the boundary value for the moduli of
C             eigenvalues.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C             On entry, the leading N-by-N part of this array must
C             contain a matrix in a real Schur form whose 1-by-1 and
C             2-by-2 diagonal blocks between positions NLOW and NSUP
C             are to be reordered.
C             On exit, the leading N-by-N part contains the ordered
C             real Schur matrix UT' * A * UT with the elements below the
C             first subdiagonal set to zero.
C             The leading NDIM-by-NDIM part of the principal submatrix
C             D = A(NLOW:NSUP,NLOW:NSUP) has eigenvalues in the domain
C             of interest and the trailing part of this submatrix has
C             eigenvalues outside the domain of interest.
C             The domain of interest for lambda(D), the eigenvalues of
C             D, is defined by the parameters ALPHA, DICO and STDOM as
C             follows:
C               For DICO = 'C':
C                  Real(lambda(D)) < ALPHA if STDOM = 'S';
C                  Real(lambda(D)) > ALPHA if STDOM = 'U'.
C               For DICO = 'D':
C                  Abs(lambda(D)) < ALPHA if STDOM = 'S';
C                  Abs(lambda(D)) > ALPHA if STDOM = 'U'.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= N.
C
C     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N)
C             On entry with JOBU = 'U', the leading N-by-N part of this
C             array must contain a transformation matrix (e.g. from a
C             previous call to this routine).
C             On exit, if JOBU = 'U', the leading N-by-N part of this
C             array contains the product of the input matrix U and the
C             orthogonal matrix UT used to reorder the diagonal blocks
C             of A.
C             On exit, if JOBU = 'I', the leading N-by-N part of this
C             array contains the matrix UT of the performed orthogonal
C             transformations.
C             Array U need not be set on entry if JOBU = 'I'.
C
C     LDU     INTEGER
C             The leading dimension of array U.  LDU >= N.
C
C     NDIM    (output) INTEGER
C             The number of eigenvalues of the selected principal
C             submatrix lying inside the domain of interest.
C             If NLOW = 1, NDIM is also the dimension of the invariant
C             subspace corresponding to the eigenvalues of the leading
C             NDIM-by-NDIM submatrix. In this case, if U is the
C             orthogonal transformation matrix used to compute and
C             reorder the real Schur form of A, its first NDIM columns
C             form an orthonormal basis for the above invariant
C             subspace.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (N)
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  A(NLOW,NLOW-1) is nonzero, i.e. A(NLOW,NLOW) is not
C                   the leading element of a 1-by-1 or 2-by-2 diagonal
C                   block of A, or A(NSUP+1,NSUP) is nonzero, i.e.
C                   A(NSUP,NSUP) is not the bottom element of a 1-by-1
C                   or 2-by-2 diagonal block of A;
C             = 2:  two adjacent blocks are too close to swap (the
C                   problem is very ill-conditioned).
C
C     METHOD
C
C     Given an upper quasi-triangular matrix A with 1-by-1 or 2-by-2
C     diagonal blocks, the routine reorders its diagonal blocks along
C     with its eigenvalues by performing an orthogonal similarity
C     transformation UT' * A * UT. The column transformation UT is also
C     performed on the given (initial) transformation U (resulted from
C     a possible previous step or initialized as the identity matrix).
C     After reordering, the eigenvalues inside the region specified by
C     the parameters ALPHA, DICO and STDOM appear at the top of
C     the selected diagonal block between positions NLOW and NSUP.
C     In other words, lambda(A(NLOW:NSUP,NLOW:NSUP)) are ordered such
C     that lambda(A(NLOW:NLOW+NDIM-1,NLOW:NLOW+NDIM-1)) are inside and
C     lambda(A(NLOW+NDIM:NSUP,NLOW+NDIM:NSUP)) are outside the domain
C     of interest. If NLOW = 1, the first NDIM columns of U*UT span the
C     corresponding invariant subspace of A.
C
C     REFERENCES
C
C     [1] Stewart, G.W.
C         HQR3 and EXCHQZ: FORTRAN subroutines for calculating and
C         ordering the eigenvalues of a real upper Hessenberg matrix.
C         ACM TOMS, 2, pp. 275-280, 1976.
C
C     NUMERICAL ASPECTS
C                                         3
C     The algorithm requires less than 4*N  operations.
C
C     CONTRIBUTOR
C
C     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen,
C     April 1998. Based on the RASP routine SEOR1.
C
C     KEYWORDS
C
C     Eigenvalues, invariant subspace, orthogonal transformation, real
C     Schur form, similarity transformation.
C
C    ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        ( ONE = 1.0D0, ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO, JOBU, STDOM
      INTEGER          INFO, LDA, LDU, N, NDIM, NLOW, NSUP
      DOUBLE PRECISION ALPHA
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), DWORK(*), U(LDU,*)
C     .. Local Scalars ..
      LOGICAL          DISCR, LSTDOM
      INTEGER          IB, L, LM1, NUP
      DOUBLE PRECISION E1, E2, TLAMBD
C     .. External Functions ..
      LOGICAL          LSAME
      DOUBLE PRECISION DLAPY2
      EXTERNAL         DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL         DLASET, DTREXC, MB03QY, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        ABS
C     .. Executable Statements ..
C
      INFO = 0
      DISCR = LSAME( DICO, 'D' )
      LSTDOM = LSAME( STDOM, 'S' )
C
C     Check input scalar arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( LSTDOM .OR. LSAME( STDOM, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LSAME( JOBU, 'I' ) .OR.
     $                 LSAME( JOBU, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.1 ) THEN
         INFO = -4
      ELSE IF( NLOW.LT.1 ) THEN
         INFO = -5
      ELSE IF( NLOW.GT.NSUP .OR. NSUP.GT.N ) THEN
         INFO = -6
      ELSE IF( DISCR .AND. ALPHA.LT.ZERO ) THEN
         INFO = -7
      ELSE IF( LDA.LT.N ) THEN
         INFO = -9
      ELSE IF( LDU.LT.N ) THEN
         INFO = -11
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MB03QD', -INFO )
         RETURN
      END IF
C
      IF( NLOW.GT.1 ) THEN
         IF( A(NLOW,NLOW-1).NE.ZERO ) INFO = 1
      END IF
      IF( NSUP.LT.N ) THEN
         IF( A(NSUP+1,NSUP).NE.ZERO ) INFO = 1
      END IF
      IF( INFO.NE.0 )
     $   RETURN
C
C     Initialize U with an identity matrix if necessary.
C
      IF( LSAME( JOBU, 'I' ) )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, U, LDU )
C
      NDIM = 0
      L = NSUP
      NUP = NSUP
C
C     NUP is the minimal value such that the submatrix A(i,j) with
C     NUP+1 <= i,j <= NSUP contains no eigenvalues inside the domain of
C     interest. L is such that all the eigenvalues of the submatrix
C     A(i,j) with L+1 <= i,j <= NUP lie inside the domain of interest.
C
C     WHILE( L >= NLOW ) DO
C
   10 IF( L.GE.NLOW ) THEN
         IB = 1
         IF( L.GT.NLOW ) THEN
            LM1 = L - 1
            IF( A(L,LM1).NE.ZERO ) THEN
               CALL MB03QY( N, LM1, A, LDA, U, LDU, E1, E2, INFO )
               IF( A(L,LM1).NE.ZERO ) IB = 2
            END IF
         END IF
         IF( DISCR ) THEN
            IF( IB.EQ.1 ) THEN
               TLAMBD = ABS( A(L,L) )
            ELSE
               TLAMBD = DLAPY2( E1, E2 )
            END IF
         ELSE
            IF( IB.EQ.1 ) THEN
               TLAMBD = A(L,L)
            ELSE
               TLAMBD = E1
            END IF
         END IF
         IF( (      LSTDOM .AND. TLAMBD.LT.ALPHA ) .OR.
     $       ( .NOT.LSTDOM .AND. TLAMBD.GT.ALPHA ) ) THEN
            NDIM = NDIM + IB
            L = L - IB
         ELSE
            IF( NDIM.NE.0 ) THEN
               CALL DTREXC( 'V', N, A, LDA, U, LDU, L, NUP, DWORK,
     $                      INFO )
               IF( INFO.NE.0 ) THEN
                  INFO = 2
                  RETURN
               END IF
               NUP = NUP - 1
               L = L - 1
            ELSE
               NUP = NUP - IB
               L = L - IB
            END IF
         END IF
         GO TO 10
      END IF
C
C     END WHILE 10
C
      RETURN
C *** Last line of MB03QD ***
      END
