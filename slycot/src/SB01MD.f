      SUBROUTINE SB01MD( NCONT, N, A, LDA, B, WR, WI, Z, LDZ, G, DWORK,
     $                   INFO )
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
C     To determine the one-dimensional state feedback matrix G of the
C     linear time-invariant single-input system
C
C           dX/dt = A * X + B * U,
C
C     where A is an NCONT-by-NCONT matrix and B is an NCONT element
C     vector such that the closed-loop system
C
C           dX/dt = (A - B * G) * X
C
C     has desired poles. The system must be preliminarily reduced
C     to orthogonal canonical form using the SLICOT Library routine
C     AB01MD.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NCONT   (input) INTEGER
C             The order of the matrix A as produced by SLICOT Library
C             routine AB01MD.  NCONT >= 0.
C
C     N       (input) INTEGER
C             The order of the matrix Z.  N >= NCONT.
C
C     A       (input/output) DOUBLE PRECISION array, dimension
C             (LDA,NCONT)
C             On entry, the leading NCONT-by-NCONT part of this array
C             must contain the canonical form of the state dynamics
C             matrix A as produced by SLICOT Library routine AB01MD.
C             On exit, the leading NCONT-by-NCONT part of this array
C             contains the upper quasi-triangular form S of the closed-
C             loop system matrix (A - B * G), that is triangular except
C             for possible 2-by-2 diagonal blocks.
C             (To reconstruct the closed-loop system matrix see
C             FURTHER COMMENTS below.)
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,NCONT).
C
C     B       (input/output) DOUBLE PRECISION array, dimension (NCONT)
C             On entry, this array must contain the canonical form of
C             the input/state vector B as produced by SLICOT Library
C             routine AB01MD.
C             On exit, this array contains the transformed vector Z * B
C             of the closed-loop system.
C
C     WR      (input) DOUBLE PRECISION array, dimension (NCONT)
C     WI      (input) DOUBLE PRECISION array, dimension (NCONT)
C             These arrays must contain the real and imaginary parts,
C             respectively, of the desired poles of the closed-loop
C             system. The poles can be unordered, except that complex
C             conjugate pairs of poles must appear consecutively.
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, the leading N-by-N part of this array must
C             contain the orthogonal transformation matrix as produced
C             by SLICOT Library routine AB01MD, which reduces the system
C             to canonical form.
C             On exit, the leading NCONT-by-NCONT part of this array
C             contains the orthogonal matrix Z which reduces the closed-
C             loop system matrix (A - B * G) to upper quasi-triangular
C             form.
C
C     LDZ     INTEGER
C             The leading dimension of array Z.  LDZ >= MAX(1,N).
C
C     G       (output) DOUBLE PRECISION array, dimension (NCONT)
C             This array contains the one-dimensional state feedback
C             matrix G of the original system.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (3*NCONT)
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
C     The method is based on the orthogonal reduction of the closed-loop
C     system matrix (A - B * G) to upper quasi-triangular form S whose
C     1-by-1 and 2-by-2 diagonal blocks correspond to the desired poles.
C     That is, S = Z'*(A - B * G)*Z, where Z is an orthogonal matrix.
C
C     REFERENCES
C
C     [1] Petkov, P. Hr.
C         A Computational Algorithm for Pole Assignment of Linear
C         Single Input Systems.
C         Internal Report 81/2, Control Systems Research Group, School
C         of Electronic Engineering and Computer Science, Kingston
C         Polytechnic, 1981.
C
C     NUMERICAL ASPECTS
C                                   3
C     The algorithm requires 0(NCONT ) operations and is backward
C     stable.
C
C     FURTHER COMMENTS
C
C     If required, the closed-loop system matrix (A - B * G) can be
C     formed from the matrix product Z * S * Z' (where S and Z are the
C     matrices output in arrays A and Z respectively).
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997.
C     Supersedes Release 2.0 routine SB01AD by Control Systems Research
C     Group, Kingston Polytechnic, United Kingdom, May 1981.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Closed loop spectrum, closed loop systems, eigenvalue assignment,
C     orthogonal canonical form, orthogonal transformation, pole
C     placement, Schur form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDZ, N, NCONT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(*), DWORK(*), G(*), WI(*), WR(*),
     $                  Z(LDZ,*)
C     .. Local Scalars ..
      LOGICAL           COMPL
      INTEGER           I, IM1, K, L, LL, LP1, NCONT2, NI, NJ, NL
      DOUBLE PRECISION  B1, P, Q, R, S, T
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DLARTG, DLASET, DROT,
     $                  DSCAL, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
      INFO = 0
C
C     Test the input scalar arguments.
C
      IF( NCONT.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.NCONT ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, NCONT ) ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return
C
         CALL XERBLA( 'SB01MD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( NCONT.EQ.0 .OR. N.EQ.0 )
     $   RETURN
C
C     Return if the system is not complete controllable.
C
      IF ( B(1).EQ.ZERO )
     $   RETURN
C
      IF ( NCONT.EQ.1 ) THEN
C
C        1-by-1 case.
C
         P = A(1,1) - WR(1)
         A(1,1) = WR(1)
         G(1) = P/B(1)
         Z(1,1) = ONE
         RETURN
      END IF
C
C     General case.  Save the contents of WI in DWORK.
C
      NCONT2 = 2*NCONT
      CALL DCOPY( NCONT, WI, 1, DWORK(NCONT2+1), 1 )
C
      B1 = B(1)
      B(1) = ONE
      L  = 0
      LL = 0
   20 CONTINUE
      L  = L + 1
      LL = LL + 1
      COMPL = DWORK(NCONT2+L).NE.ZERO
      IF ( L.NE.NCONT ) THEN
         LP1 = L + 1
         NL = NCONT - L
         IF ( LL.NE.2 ) THEN
            IF ( COMPL ) THEN
C
C              Compute complex eigenvector.
C
               DWORK(NCONT) = ONE
               DWORK(NCONT2) = ONE
               P = WR(L)
               T = DWORK(NCONT2+L)
               Q = T*DWORK(NCONT2+LP1)
               DWORK(NCONT2+L) = ONE
               DWORK(NCONT2+LP1) = Q
C
               DO 40 I = NCONT, LP1, -1
                  IM1 = I - 1
                  DWORK(IM1) = ( P*DWORK(I) + Q*DWORK(NCONT+I) -
     $               DDOT( NCONT-IM1, A(I,I), LDA, DWORK(I), 1 ) )
     $               /A(I,IM1)
                  DWORK(NCONT+IM1) = ( P*DWORK(NCONT+I) + DWORK(I) -
     $               DDOT( NCONT-IM1, A(I,I), LDA, DWORK(NCONT+I), 1 ) )
     $               /A(I,IM1)
   40          CONTINUE
C
            ELSE
C
C              Compute real eigenvector.
C
               DWORK(NCONT) = ONE
               P = WR(L)
C
               DO 60 I = NCONT, LP1, -1
                  IM1 = I - 1
                  DWORK(IM1) = ( P*DWORK(I) -
     $               DDOT( NCONT-IM1, A(I,I), LDA, DWORK(I), 1 ) )
     $               /A(I,IM1)
   60          CONTINUE
C
            END IF
         END IF
C
C        Transform eigenvector.
C
         DO 80 K = NCONT - 1, L, -1
            IF ( LL.NE.2 ) THEN
               R = DWORK(K)
               S = DWORK(K+1)
            ELSE
               R = DWORK(NCONT+K)
               S = DWORK(NCONT+K+1)
            END IF
            CALL DLARTG( R, S, P, Q, T )
            DWORK(K) = T
            IF ( LL.NE.2 ) THEN
               NJ = MAX( K-1, L )
            ELSE
               DWORK(NCONT+K) = T
               NJ = L - 1
            END IF
C
C           Transform  A.
C
            CALL DROT( NCONT-NJ+1, A(K,NJ), LDA, A(K+1,NJ), LDA, P, Q )
C
            IF ( COMPL .AND. LL.EQ.1 ) THEN
               NI = NCONT
            ELSE
               NI = MIN( K+2, NCONT )
            END IF
            CALL DROT( NI, A(1,K), 1, A(1,K+1), 1, P, Q )
C
            IF ( K.EQ.L ) THEN
C
C              Transform  B.
C
               T = B(K)
               B(K) = P*T
               B(K+1) = -Q*T
            END IF
C
C           Accumulate transformations.
C
            CALL DROT( NCONT, Z(1,K), 1, Z(1,K+1), 1, P, Q )
C
            IF ( COMPL .AND. LL.NE.2 ) THEN
               T = DWORK(NCONT+K)
               DWORK(NCONT+K) = P*T + Q*DWORK(NCONT+K+1)
               DWORK(NCONT+K+1) = P*DWORK(NCONT+K+1) - Q*T
            END IF
   80    CONTINUE
C
      END IF
C
      IF ( .NOT.COMPL ) THEN
C
C        Find one element of  G.
C
         K = L
         R = B(L)
         IF ( L.NE.NCONT ) THEN
            IF ( ABS( B(LP1) ).GT.ABS( B(L) ) ) THEN
               K = LP1
               R = B(LP1)
            END IF
         END IF
         P = A(K,L)
         IF ( K.EQ.L ) P = P - WR(L)
         P = P/R
C
         CALL DAXPY( LP1, -P, B, 1, A(1,L), 1 )
C
         G(L) = P/B1
         IF ( L.NE.NCONT ) THEN
            LL = 0
            GO TO 20
         END IF
      ELSE IF ( LL.EQ.1 ) THEN
         GO TO 20
      ELSE
C
C        Find two elements of  G.
C
         K = L
         R = B(L)
         IF ( L.NE.NCONT ) THEN
            IF ( ABS( B(LP1)).GT.ABS( B(L) ) ) THEN
               K = LP1
               R = B(LP1)
            END IF
         END IF
         P = A(K,L-1)
         Q = A(K,L)
         IF ( K.EQ.L ) THEN
            P = P - ( DWORK(NCONT+L)/DWORK(L-1) )*DWORK(NCONT2+L)
            Q = Q - WR(L) +
     $            ( DWORK(NCONT+L-1)/DWORK(L-1) )*DWORK(NCONT2+L)
         END IF
         P = P/R
         Q = Q/R
C
         CALL DAXPY( LP1, -P, B, 1, A(1,L-1), 1 )
         CALL DAXPY( LP1, -Q, B, 1, A(1,L), 1 )
C
         G(L-1) = P/B1
         G(L) = Q/B1
         IF ( L.NE.NCONT ) THEN
            LL = 0
            GO TO 20
         END IF
      END IF
C
C     Transform  G.
C
      CALL DGEMV( 'No transpose', NCONT, NCONT, ONE, Z, LDZ, G, 1,
     $            ZERO, DWORK, 1 )
      CALL DCOPY( NCONT, DWORK, 1, G, 1 )
      CALL DSCAL( NCONT, B1, B, 1 )
C
C     Annihilate A after the first subdiagonal.
C
      IF ( NCONT.GT.2 )
     $   CALL DLASET( 'Lower', NCONT-2, NCONT-2, ZERO, ZERO, A(3,1),
     $                LDA )
C
      RETURN
C *** Last line of SB01MD ***
      END
