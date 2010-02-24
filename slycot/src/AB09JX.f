      SUBROUTINE AB09JX( DICO, STDOM, EVTYPE, N, ALPHA, ER, EI, ED,
     $                   TOLINF, INFO )
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
C     To check stability/antistability of finite eigenvalues with
C     respect to a given stability domain.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the stability domain as follows:
C             = 'C':  for a continuous-time system;
C             = 'D':  for a discrete-time system.
C
C     STDOM   CHARACTER*1
C             Specifies whether the domain of interest is of stability
C             type (left part of complex plane or inside of a circle)
C             or of instability type (right part of complex plane or
C             outside of a circle) as follows:
C             = 'S':  stability type domain;
C             = 'U':  instability type domain.
C
C     EVTYPE  CHARACTER*1
C             Specifies whether the eigenvalues arise from a standard
C             or a generalized eigenvalue problem as follows:
C             = 'S':  standard eigenvalue problem;
C             = 'G':  generalized eigenvalue problem;
C             = 'R':  reciprocal generalized eigenvalue problem.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The dimension of vectors ER, EI and ED.  N >= 0.
C
C     ALPHA   (input) DOUBLE PRECISION
C             Specifies the boundary of the domain of interest for the
C             eigenvalues. For a continuous-time system
C             (DICO = 'C'), ALPHA is the boundary value for the real
C             parts of eigenvalues, while for a discrete-time system
C             (DICO = 'D'), ALPHA >= 0 represents the boundary value for
C             the moduli of eigenvalues.
C
C     ER, EI, (input) DOUBLE PRECISION arrays, dimension (N)
C     ED      If EVTYPE = 'S', ER(j) + EI(j)*i, j = 1,...,N, are
C             the eigenvalues of a real matrix.
C             ED is not referenced and is implicitly considered as
C             a vector having all elements equal to one.
C             If EVTYPE = 'G' or EVTYPE = 'R', (ER(j) + EI(j)*i)/ED(j),
C             j = 1,...,N, are the generalized eigenvalues of a pair of
C             real matrices. If ED(j) is zero, then the j-th generalized
C             eigenvalue is infinite.
C             Complex conjugate pairs of eigenvalues must appear
C             consecutively.
C
C     Tolerances
C
C     TOLINF  DOUBLE PRECISION
C             If EVTYPE = 'G' or 'R', TOLINF contains the tolerance for
C             detecting infinite generalized eigenvalues.
C             0 <= TOLINF < 1.
C
C     Error Indicator
C
C     INFO    INTEGER
C             =  0:  successful exit, i.e., all eigenvalues lie within
C                    the domain of interest defined by DICO, STDOM
C                    and ALPHA;
C             <  0:  if INFO = -i, the i-th argument had an illegal
C                    value;
C             =  1:  some eigenvalues lie outside the domain of interest
C                    defined by DICO, STDOM and ALPHA.
C     METHOD
C
C     The domain of interest for an eigenvalue lambda is defined by the
C     parameters ALPHA, DICO and STDOM as follows:
C        - for a continuous-time system (DICO = 'C'):
C               Real(lambda) < ALPHA if STDOM = 'S';
C               Real(lambda) > ALPHA if STDOM = 'U';
C        - for a discrete-time system (DICO = 'D'):
C               Abs(lambda) < ALPHA if STDOM = 'S';
C               Abs(lambda) > ALPHA if STDOM = 'U'.
C     If EVTYPE = 'R', the same conditions apply for 1/lambda.
C
C     CONTRIBUTORS
C
C     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001.
C
C     KEYWORDS
C
C     Stability.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      CHARACTER        DICO, EVTYPE, STDOM
      INTEGER          INFO, N
      DOUBLE PRECISION ALPHA, TOLINF
C     .. Array Arguments ..
      DOUBLE PRECISION ED(*), EI(*), ER(*)
C     .. Local Scalars
      LOGICAL          DISCR, RECEVP, STAB, STDEVP
      DOUBLE PRECISION ABSEV, RPEV, SCALE
      INTEGER          I
C     .. External Functions ..
      LOGICAL          LSAME
      DOUBLE PRECISION DLAPY2
      EXTERNAL         DLAPY2, LSAME
C     .. External Subroutines ..
      EXTERNAL         XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        ABS
C     .. Executable Statements ..
C
      INFO   = 0
      DISCR  = LSAME( DICO,   'D' )
      STAB   = LSAME( STDOM,  'S' )
      STDEVP = LSAME( EVTYPE, 'S' )
      RECEVP = LSAME( EVTYPE, 'R' )
C
C     Check the scalar input arguments.
C
      IF( .NOT. ( LSAME( DICO, 'C' ) .OR. DISCR ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( STAB .OR. LSAME( STDOM, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( STDEVP .OR. LSAME( EVTYPE, 'G' ) .OR.
     $                 RECEVP ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( DISCR .AND. ALPHA.LT.ZERO ) THEN
         INFO = -5
      ELSE IF( TOLINF.LT.ZERO .OR. TOLINF.GE.ONE ) THEN
         INFO = -9
      END IF
C
      IF( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'AB09JX', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      IF( STAB ) THEN
C
C        Check the stability of finite eigenvalues.
C
         SCALE = ONE
         IF( DISCR ) THEN
            DO 10 I = 1, N
               ABSEV = DLAPY2( ER(I), EI(I) )
               IF( RECEVP ) THEN
                  SCALE = ABSEV
                  ABSEV = ABS( ED(I) )
               ELSE IF( .NOT.STDEVP ) THEN
                  SCALE = ED(I)
               END IF
               IF( ABS( SCALE ).GT.TOLINF .AND.
     $            ABSEV.GE.ALPHA*SCALE ) THEN
                  INFO = 1
                  RETURN
               END IF
   10       CONTINUE
         ELSE
            DO 20 I = 1, N
               RPEV = ER(I)
               IF( RECEVP ) THEN
                  SCALE = RPEV
                  RPEV = ED(I)
               ELSE IF( .NOT.STDEVP ) THEN
                  SCALE = ED(I)
               END IF
               IF( ABS( SCALE ).GT.TOLINF .AND.
     $            RPEV.GE.ALPHA*SCALE ) THEN
                  INFO = 1
                  RETURN
               END IF
   20       CONTINUE
         END IF
      ELSE
C
C        Check the anti-stability of finite eigenvalues.
C
         IF( DISCR ) THEN
            DO 30 I = 1, N
               ABSEV = DLAPY2( ER(I), EI(I) )
               IF( RECEVP ) THEN
                  SCALE = ABSEV
                  ABSEV = ABS( ED(I) )
               ELSE IF( .NOT.STDEVP ) THEN
                  SCALE = ED(I)
               END IF
               IF( ABS( SCALE ).GT.TOLINF .AND.
     $            ABSEV.LE.ALPHA*SCALE ) THEN
                  INFO = 1
                  RETURN
               END IF
   30       CONTINUE
         ELSE
            DO 40 I = 1, N
               RPEV = ER(I)
               IF( RECEVP ) THEN
                  SCALE = RPEV
                  RPEV = ED(I)
               ELSE IF( .NOT.STDEVP ) THEN
                  SCALE = ED(I)
               END IF
               IF( ABS( SCALE ).GT.TOLINF .AND.
     $            RPEV.LE.ALPHA*SCALE ) THEN
                  INFO = 1
                  RETURN
               END IF
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C *** Last line of AB09JX ***
      END
