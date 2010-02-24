      SUBROUTINE MB02CU( TYPEG, K, P, Q, NB, A1, LDA1, A2, LDA2, B, LDB,
     $                   RNK, IPVT, CS, TOL, DWORK, LDWORK, INFO )
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
C     To bring the first blocks of a generator to proper form.
C     The positive part of the generator is contained in the arrays A1
C     and A2. The negative part of the generator is contained in B.
C     Transformation information will be stored and can be applied
C     via SLICOT Library routine MB02CV.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     TYPEG   CHARACTER*1
C             Specifies the type of the generator, as follows:
C             = 'D':  generator is column oriented and rank
C                     deficiencies are expected;
C             = 'C':  generator is column oriented and rank
C                     deficiencies are not expected;
C             = 'R':  generator is row oriented and rank
C                     deficiencies are not expected.
C
C     Input/Output Parameters
C
C     K       (input)  INTEGER
C             The number of rows in A1 to be processed.  K >= 0.
C
C     P       (input)  INTEGER
C             The number of columns of the positive generator.  P >= K.
C
C     Q       (input)  INTEGER
C             The number of columns in B containing the negative
C             generators.
C             If TYPEG = 'D',        Q >= K;
C             If TYPEG = 'C' or 'R', Q >= 0.
C
C     NB      (input)  INTEGER
C             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies
C             the block size to be used in the blocked parts of the
C             algorithm. If NB <= 0, an unblocked algorithm is used.
C
C     A1      (input/output)  DOUBLE PRECISION array, dimension
C             (LDA1, K)
C             On entry, the leading K-by-K part of this array must
C             contain the leading submatrix of the positive part of the
C             generator. If TYPEG = 'C', A1 is assumed to be lower
C             triangular and the strictly upper triangular part is not
C             referenced. If TYPEG = 'R', A1 is assumed to be upper
C             triangular and the strictly lower triangular part is not
C             referenced.
C             On exit, if TYPEG = 'D', the leading K-by-RNK part of this
C             array contains the lower trapezoidal part of the proper
C             generator and information for the Householder
C             transformations applied during the reduction process.
C             On exit, if TYPEG = 'C', the leading K-by-K part of this
C             array contains the leading lower triangular part of the
C             proper generator.
C             On exit, if TYPEG = 'R', the leading K-by-K part of this
C             array contains the leading upper triangular part of the
C             proper generator.
C
C     LDA1    INTEGER
C             The leading dimension of the array A1.  LDA1 >= MAX(1,K).
C
C     A2      (input/output)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K);
C             if TYPEG = 'R',                   dimension (LDA2, K).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-(P-K) part of this array must contain the (K+1)-st
C             to P-th columns of the positive part of the generator.
C             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of
C             this array must contain the (K+1)-st to P-th rows of the
C             positive part of the generator.
C             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-(P-K) part of this array contains information for
C             Householder transformations.
C             On exit, if TYPEG = 'R', the leading (P-K)-by-K part of
C             this array contains information for Householder
C             transformations.
C
C     LDA2    INTEGER
C             The leading dimension of the array A2.
C             If P = K,                   LDA2 >= 1;
C             If P > K and (TYPEG = 'D' or TYPEG = 'C'),
C                                         LDA2 >= MAX(1,K);
C             if P > K and TYPEG = 'R',   LDA2 >= P-K.
C
C     B       (input/output)  DOUBLE PRECISION array,
C             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q);
C             if TYPEG = 'R',                   dimension (LDB, K).
C             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-Q part of this array must contain the negative part
C             of the generator.
C             On entry, if TYPEG = 'R', the leading Q-by-K part of this
C             array must contain the negative part of the generator.
C             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading
C             K-by-Q part of this array contains information for
C             Householder transformations.
C             On exit, if TYPEG = 'R', the leading Q-by-K part of this
C             array contains information for Householder transformations.
C
C     LDB     INTEGER
C             The leading dimension of the array B.
C             If Q = 0,                  LDB >= 1;
C             if Q > 0 and (TYPEG = 'D' or TYPEG = 'C'),
C                                        LDB >= MAX(1,K);
C             if Q > 0 and TYPEG = 'R',  LDB >= Q.
C
C     RNK     (output)  INTEGER
C             If TYPEG = 'D', the number of columns in the reduced
C             generator which are found to be linearly independent.
C             If TYPEG = 'C' or TYPEG = 'R', then RNK is not set.
C
C     IPVT    (output)  INTEGER array, dimension (K)
C             If TYPEG = 'D', then if IPVT(i) = k, the k-th row of the
C             proper generator is the reduced i-th row of the input
C             generator.
C             If TYPEG = 'C' or TYPEG = 'R', this array is not
C             referenced.
C
C     CS      (output)  DOUBLE PRECISION array, dimension (x)
C             If TYPEG = 'D' and P = K,                   x = 3*K;
C             if TYPEG = 'D' and P > K,                   x = 5*K;
C             if (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K;
C             if (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K.
C             On exit, the first x elements of this array contain
C             necessary information for the SLICOT library routine
C             MB02CV (Givens and modified hyperbolic rotation
C             parameters, scalar factors of the Householder
C             transformations).
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION
C             If TYPEG = 'D', this number specifies the used tolerance
C             for handling deficiencies. If the hyperbolic norm
C             of two diagonal elements in the positive and negative
C             generators appears to be less than or equal to TOL, then
C             the corresponding columns are not reduced.
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = -17,  DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.
C             LDWORK >= MAX(1,4*K),         if TYPEG = 'D';
C             LDWORK >= MAX(1,MAX(NB,1)*K), if TYPEG = 'C' or 'R'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  if TYPEG = 'D', the generator represents a
C                   (numerically) indefinite matrix; and if TYPEG = 'C'
C                   or TYPEG = 'R', the generator represents a
C                   (numerically) semidefinite matrix.
C
C     METHOD
C
C     If TYPEG = 'C' or TYPEG = 'R', blocked Householder transformations
C     and modified hyperbolic rotations are used to downdate the
C     matrix [ A1  A2  sqrt(-1)*B ], cf. [1], [2].
C     If TYPEG = 'D', then an algorithm with row pivoting is used. In
C     the first stage it maximizes the hyperbolic norm of the active
C     row. As soon as the hyperbolic norm is below the threshold TOL,
C     the strategy is changed. Now, in the second stage, the algorithm
C     applies an LQ decomposition with row pivoting on B such that
C     the Euclidean norm of the active row is maximized.
C
C     REFERENCES
C
C     [1] Kailath, T. and Sayed, A.
C         Fast Reliable Algorithms for Matrices with Structure.
C         SIAM Publications, Philadelphia, 1999.
C
C     [2] Kressner, D. and Van Dooren, P.
C         Factorizations and linear system solvers for matrices with
C         Toeplitz structure.
C         SLICOT Working Note 2000-2, 2000.
C
C     NUMERICAL ASPECTS
C                               2
C     The algorithm requires 0(K *( P + Q )) floating point operations.
C
C     CONTRIBUTOR
C
C     D. Kressner, Technical Univ. Berlin, Germany, May 2001.
C
C     REVISIONS
C
C     V. Sima, Research Institute for Informatics, Bucharest, June 2001.
C     D. Kressner, Technical Univ. Berlin, Germany, July 2002.
C     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004.
C
C     KEYWORDS
C
C     Elementary matrix operations, Householder transformation, matrix
C     operations, Toeplitz matrix.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, P05
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, P05 = 0.05D0 )
C     .. Scalar Arguments ..
      CHARACTER          TYPEG
      INTEGER            INFO, K, LDA1, LDA2, LDB, LDWORK, NB, P, Q, RNK
      DOUBLE PRECISION   TOL
C     .. Array Arguments ..
      INTEGER            IPVT(*)
      DOUBLE PRECISION   A1(LDA1,*), A2(LDA2,*), B(LDB,*), CS(*),
     $                   DWORK(*)
C     .. Local Scalars ..
      LOGICAL            LCOL, LRDEF
      INTEGER            COL2, I, IB, IERR, IMAX, ITEMP, J, JJ, LEN,
     $                   NBL, PDW, PHV, POS, PST2, PVT, WRKMIN
      DOUBLE PRECISION   ALPHA, ALPHA2, BETA, C, DMAX, S, TAU1, TAU2,
     $                   TEMP, TEMP2
C     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAPY2, DNRM2
      EXTERNAL           IDAMAX, DLAPY2, DNRM2, LSAME
C     .. External Subroutines ..
      EXTERNAL           DAXPY, DGELQ2, DGEQR2, DLARF, DLARFB, DLARFG,
     $                   DLARFT, DLARTG, DROT, DSCAL, DSWAP, MA02FD,
     $                   XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
C
C     .. Executable Statements ..
C
C     Decode the scalar input parameters.
C
      INFO  = 0
      COL2  = P - K
      LRDEF = LSAME( TYPEG, 'D' )
      LCOL  = LSAME( TYPEG, 'C' )
      IF ( LRDEF ) THEN
         WRKMIN = MAX( 1, 4*K )
      ELSE
         WRKMIN = MAX( 1, NB*K, K )
      END IF
C
C     Check the scalar input parameters.
C
      IF ( .NOT.( LCOL .OR. LRDEF .OR. LSAME( TYPEG, 'R' ) ) ) THEN
         INFO = -1
      ELSE IF ( K.LT.0 ) THEN
         INFO = -2
      ELSE IF ( P.LT.K ) THEN
         INFO = -3
      ELSE IF ( Q.LT.0 .OR. ( LRDEF .AND. Q.LT.K ) ) THEN
         INFO = -4
      ELSE IF ( LDA1.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF ( ( ( P.EQ.K ) .AND. LDA2.LT.1 ) .OR.
     $          ( ( P.GT.K ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDA2.LT.MAX( 1, K ) ) ) .OR.
     $          ( ( P.GT.K ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDA2.LT.( P - K ) ) ) ) THEN
         INFO = -9
      ELSE IF ( ( ( Q.EQ.0 ) .AND. LDB.LT.1 ) .OR.
     $          ( ( Q.GT.0 ) .AND. ( LRDEF .OR. LCOL ) .AND.
     $            ( LDB.LT.MAX( 1, K ) ) ) .OR.
     $          ( ( Q.GT.0 ) .AND. .NOT.( LRDEF .OR. LCOL ) .AND.
     $            ( LDB.LT.Q ) ) ) THEN
         INFO = -11
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         DWORK(1) = DBLE( WRKMIN )
         INFO = -17
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'MB02CU', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( K.EQ.0 .OR. ( .NOT.LRDEF .AND. Q.EQ.0 .AND. P.EQ.K ) ) THEN
         IF ( LRDEF )
     $      RNK = 0
         RETURN
      END IF
C
      IF ( LRDEF ) THEN
C
C        Deficient generator.
C
         IF ( COL2.EQ.0 ) THEN
            PST2 = 2*K
         ELSE
            PST2 = 4*K
         END IF
C
C        Initialize partial hyperbolic row norms.
C
         RNK = 0
         PHV = 3*K
C
         DO 10 I = 1, K
            IPVT(I) = I
            DWORK(I) = DNRM2( K, A1(I,1), LDA1 )
   10    CONTINUE
C
         DO 20 I = 1, K
            DWORK(I) = DLAPY2( DWORK(I),
     $                         DNRM2( COL2, A2(I,1), LDA2 ) )
            DWORK(I+K) = DWORK(I)
   20    CONTINUE
C
         PDW = 2*K
C
         DO 30 I = 1, K
            PDW = PDW + 1
            DWORK(PDW) = DNRM2( Q, B(I,1), LDB )
   30    CONTINUE
C
C        Compute factorization.
C
         DO 90 I = 1, K
C
C           Determine pivot row and swap if necessary.
C
            PDW   = I
            ALPHA = ABS( DWORK(PDW) )
            BETA  = ABS( DWORK(PDW+2*K) )
            DMAX  = SIGN( SQRT( ABS( ALPHA - BETA ) )*
     $                    SQRT( ALPHA + BETA ), ALPHA - BETA )
            IMAX  = I
C
            DO 40 J = 1, K - I
               PDW   = PDW + 1
               ALPHA = ABS( DWORK(PDW) )
               BETA  = ABS ( DWORK(PDW+2*K) )
               TEMP  = SIGN( SQRT( ABS( ALPHA - BETA ) )*
     $                       SQRT( ALPHA + BETA ), ALPHA - BETA )
               IF ( TEMP.GT.DMAX ) THEN
                  IMAX = I + J
                  DMAX = TEMP
               END IF
   40       CONTINUE
C
C           Proceed with the reduction if the hyperbolic norm is
C           beyond the threshold.
C
            IF ( DMAX.GT.TOL ) THEN
C
               PVT = IMAX
               IF ( PVT.NE.I ) THEN
                  CALL DSWAP( K, A1(PVT,1), LDA1, A1(I,1), LDA1 )
                  CALL DSWAP( COL2, A2(PVT,1), LDA2, A2(I,1), LDA2 )
                  CALL DSWAP( Q, B(PVT,1), LDB, B(I,1), LDB )
                  ITEMP = IPVT(PVT)
                  IPVT(PVT) = IPVT(I)
                  IPVT(I) = ITEMP
                  DWORK(PVT) = DWORK(I)
                  DWORK(K+PVT) = DWORK(K+I)
                  DWORK(2*K+PVT) = DWORK(2*K+I)
               END IF
C
C              Generate and apply elementary reflectors.
C
               IF ( COL2.GT.1 ) THEN
                  CALL DLARFG( COL2, A2(I,1), A2(I,2), LDA2, TAU2 )
                  ALPHA2 = A2(I,1)
                  IF ( K.GT.I ) THEN
                     A2(I,1) = ONE
                     CALL DLARF( 'Right', K-I, COL2, A2(I,1), LDA2,
     $                           TAU2, A2(I+1,1), LDA2, DWORK(PHV+1) )
                  END IF
                  A2(I,1) = TAU2
               ELSE IF ( COL2.GT.0 ) THEN
                  ALPHA2  = A2(I,1)
                  A2(I,1) = ZERO
               END IF
C
               IF ( K.GT.I ) THEN
                  CALL DLARFG( K-I+1, A1(I,I), A1(I,I+1), LDA1, TAU1 )
                  ALPHA   = A1(I,I)
                  A1(I,I) = ONE
                  CALL DLARF( 'Right', K-I, K-I+1, A1(I,I), LDA1, TAU1,
     $                        A1(I+1,I), LDA1, DWORK(PHV+1) )
                  CS(PST2+I) = TAU1
               ELSE
                  ALPHA = A1(I,I)
               END IF
C
               IF ( COL2.GT.0 ) THEN
                  TEMP = ALPHA
                  CALL DLARTG( TEMP, ALPHA2, C, S, ALPHA )
                  IF ( K.GT.I )
     $               CALL DROT( K-I, A1(I+1,I), 1, A2(I+1,1), 1, C, S )
                  CS(2*K+I*2-1) = C
                  CS(2*K+I*2)   = S
               END IF
               A1(I,I) = ALPHA
C
               IF ( Q.GT.1 ) THEN
                  CALL DLARFG( Q, B(I,1), B(I,2), LDB, TAU2 )
                  BETA = B(I,1)
                  IF ( K.GT.I ) THEN
                     B(I,1) = ONE
                     CALL DLARF( 'Right', K-I, Q, B(I,1), LDB, TAU2,
     $                           B(I+1,1), LDB, DWORK(PHV+1) )
                  END IF
                  B(I,1) = TAU2
               ELSE IF ( Q.GT.0 ) THEN
                  BETA   = B(I,1)
                  B(I,1) = ZERO
               ELSE
                  BETA = ZERO
               END IF
C
C              Create hyperbolic Givens rotation.
C
               CALL MA02FD( A1(I,I), BETA, C, S, IERR )
               IF ( IERR.NE.0 ) THEN
C
C                 Error return:  This should not happen.
C
                  INFO = 1
                  RETURN
               END IF
C
C              Apply hyperbolic rotation.
C
               IF ( K.GT.I ) THEN
                  CALL DSCAL( K-I, ONE/C, A1(I+1,I), 1 )
                  CALL DAXPY( K-I,  -S/C, B(I+1,1), 1, A1(I+1,I), 1 )
                  CALL DSCAL( K-I,     C, B(I+1,1), 1 )
                  CALL DAXPY( K-I,    -S, A1(I+1,I), 1, B(I+1,1), 1 )
               END IF
               CS(I*2-1) = C
               CS(I*2)   = S
C
C              Downdate the norms in A1.
C
               DO 50  J = I + 1, K
                  TEMP  = ONE - ( ABS( A1(J,I) ) / DWORK(J) )**2
                  TEMP2 = ONE + P05*TEMP*
     $                          ( DWORK(J) / DWORK(K+J) )**2
                  IF ( TEMP2.EQ.ONE ) THEN
                     DWORK(J) = DLAPY2( DNRM2( K-I, A1(J,I+1), LDA1 ),
     $                                  DNRM2( COL2, A2(J,1), LDA2 ) )
                     DWORK(K+J)   = DWORK(J)
                     DWORK(2*K+J) = DNRM2( Q, B(J,1), LDB )
                  ELSE
                     IF ( TEMP.GE.ZERO ) THEN
                        DWORK(J) =  DWORK(J)*SQRT( TEMP )
                     ELSE
                        DWORK(J) = -DWORK(J)*SQRT( -TEMP )
                     END IF
                  END IF
   50          CONTINUE
C
               RNK = RNK + 1
            ELSE IF ( ABS( DMAX ).LT.TOL ) THEN
C
C              Displacement is positive semidefinite.
C              Do an LQ decomposition with pivoting of the leftover
C              negative part to find diagonal elements with almost zero
C              norm. These columns cannot be removed from the
C              generator.
C
C              Initialize norms.
C
               DO 60  J = I, K
                  DWORK(J)   = DNRM2( Q, B(J,1), LDB )
                  DWORK(J+K) = DWORK(J)
   60          CONTINUE
C
               LEN = Q
               POS = 1
C
               DO 80  J = I, K
C
C                 Generate and apply elementary reflectors.
C
                  PVT = ( J-1 ) + IDAMAX( K-J+1, DWORK(J), 1 )
C
C                 Swap rows if necessary.
C
                  IF ( PVT.NE.J ) THEN
                     CALL DSWAP( K, A1(PVT,1), LDA1, A1(J,1), LDA1 )
                     CALL DSWAP( COL2, A2(PVT,1), LDA2, A2(J,1), LDA2 )
                     CALL DSWAP( Q, B(PVT,1), LDB, B(J,1), LDB )
                     ITEMP = IPVT(PVT)
                     IPVT(PVT) = IPVT(J)
                     IPVT(J) = ITEMP
                     DWORK(PVT) = DWORK(J)
                     DWORK(K+PVT) = DWORK(K+J)
                  END IF
C
C                 Annihilate second part of the positive generators.
C
                  IF ( COL2.GT.1 ) THEN
                     CALL DLARFG( COL2, A2(J,1), A2(J,2), LDA2, TAU2 )
                     ALPHA2 = A2(J,1)
                     IF ( K.GT.J ) THEN
                        A2(J,1) = ONE
                        CALL DLARF( 'Right', K-J, COL2, A2(J,1), LDA2,
     $                              TAU2, A2(J+1,1), LDA2, DWORK(PHV+1))
                     END IF
                     A2(J,1) = TAU2
                  ELSE IF ( COL2.GT.0 ) THEN
                     ALPHA2  = A2(J,1)
                     A2(J,1) = ZERO
                  END IF
C
C                 Transform first part of the positive generators to
C                 lower triangular form.
C
                  IF ( K.GT.J ) THEN
                     CALL DLARFG( K-J+1, A1(J,J), A1(J,J+1), LDA1,
     $                            TAU1 )
                     ALPHA   = A1(J,J)
                     A1(J,J) = ONE
                     CALL DLARF( 'Right', K-J, K-J+1, A1(J,J), LDA1,
     $                           TAU1, A1(J+1,J), LDA1, DWORK(PHV+1) )
                     CS(PST2+J) = TAU1
                  ELSE
                     ALPHA = A1(J,J)
                  END IF
C
                  IF ( COL2.GT.0 ) THEN
                     TEMP = ALPHA
                     CALL DLARTG( TEMP, ALPHA2, C, S, ALPHA )
                     IF ( K.GT.J )
     $                  CALL DROT( K-J, A1(J+1,J), 1, A2(J+1,1), 1, C,
     $                             S )
                     CS(2*K+J*2-1) = C
                     CS(2*K+J*2)   = S
                  END IF
                  A1(J,J) = ALPHA
C
C                 Transform negative part to lower triangular form.
C
                  IF ( LEN.GT.1) THEN
                     CALL DLARFG( LEN, B(J,POS), B(J,POS+1), LDB, TAU2 )
                     BETA = B(J,POS)
                     IF ( K.GT.J ) THEN
                        B(J,POS) = ONE
                        CALL DLARF( 'Right', K-J, LEN, B(J,POS), LDB,
     $                              TAU2, B(J+1,POS), LDB, DWORK(PHV+1))
                     END IF
                     B(J,POS)  = BETA
                     CS(J*2-1) = TAU2
                  END IF
C
C                 Downdate the norms of the rows in the negative part.
C
                  DO 70  JJ = J + 1, K
                     IF ( DWORK(JJ).NE.ZERO ) THEN
                        TEMP  = ONE - ( ABS( B(JJ,POS) )
     $                                 / DWORK(JJ) )**2
                        TEMP  = MAX( TEMP, ZERO )
                        TEMP2 = ONE + P05*TEMP*
     $                               ( DWORK(JJ) / DWORK(K+JJ) )**2
                        IF ( TEMP2.EQ.ONE ) THEN
                           DWORK(JJ)   = DNRM2( LEN-1, B(JJ,POS+1), LDB)
                           DWORK(K+JJ) = DWORK(JJ)
                        ELSE
                           IF ( TEMP.GE.ZERO ) THEN
                              DWORK(JJ) =  DWORK(JJ)*SQRT( TEMP )
                           ELSE
                              DWORK(JJ) = -DWORK(JJ)*SQRT( -TEMP )
                           END IF
                        END IF
                     END IF
   70             CONTINUE
C
                  LEN = LEN - 1
                  POS = POS + 1
   80          CONTINUE
C
               RETURN
            ELSE
C
C              Error return:
C
C              Displacement is indefinite.
C              Due to roundoff error, positive semidefiniteness is
C              violated. This is a rather bad situation. There is no
C              meaningful way to continue the computations from this
C              point.
C
               INFO = 1
               RETURN
            END IF
   90    CONTINUE
C
      ELSE IF ( LCOL ) THEN
C
C        Column oriented and not deficient generator.
C
C        Apply an LQ like hyperbolic/orthogonal blocked decomposition.
C
         IF ( COL2.GT.0 ) THEN
            NBL = MIN( COL2, NB )
            IF ( NBL.GT.0 ) THEN
C
C              Blocked version.
C
               DO 110  I = 1, K - NBL + 1, NBL
                  IB = MIN( K-I+1, NBL )
                  CALL DGELQ2( IB, COL2, A2(I,1), LDA2, CS(4*K+I),
     $                         DWORK, IERR )
                  IF ( I+IB.LE.K ) THEN
                     CALL DLARFT( 'Forward', 'Rowwise', COL2, IB,
     $                            A2(I,1), LDA2, CS(4*K+I), DWORK, K )
                     CALL DLARFB( 'Right', 'No Transpose', 'Forward',
     $                            'Rowwise', K-I-IB+1, COL2, IB,
     $                            A2(I,1), LDA2, DWORK, K, A2(I+IB,1),
     $                            LDA2, DWORK(IB+1), K )
                  END IF
C
C                 Annihilate the remaining parts of A2.
C
                  DO 100  J = I, I + IB - 1
                     IF ( COL2.GT.1 ) THEN
                        LEN = MIN( COL2, J-I+1 )
                        CALL DLARFG( LEN, A2(J,1), A2(J,2), LDA2, TAU2 )
                        ALPHA2 = A2(J,1)
                        IF ( K.GT.J ) THEN
                           A2(J,1) = ONE
                           CALL DLARF( 'Right', K-J, LEN, A2(J,1), LDA2,
     $                                 TAU2, A2(J+1,1), LDA2, DWORK )
                        END IF
                        A2(J,1) = TAU2
                     ELSE
                        ALPHA2  = A2(J,1)
                        A2(J,1) = ZERO
                     END IF
                     ALPHA = A1(J,J)
                     CALL DLARTG( ALPHA, ALPHA2, C, S, A1(J,J) )
                     IF ( K.GT.J )
     $                  CALL DROT( K-J, A1(J+1,J), 1, A2(J+1,1), 1, C,
     $                             S )
                     CS(2*K+J*2-1) = C
                     CS(2*K+J*2)   = S
  100             CONTINUE
C
  110          CONTINUE
C
            ELSE
               I = 1
            END IF
C
C           Unblocked version for the last or only block.
C
            DO 120  J = I, K
               IF ( COL2.GT.1 ) THEN
                  CALL DLARFG( COL2, A2(J,1), A2(J,2), LDA2, TAU2 )
                  ALPHA2 = A2(J,1)
                  IF ( K.GT.J ) THEN
                     A2(J,1) = ONE
                     CALL DLARF( 'Right', K-J, COL2, A2(J,1), LDA2,
     $                           TAU2, A2(J+1,1), LDA2, DWORK )
                  END IF
                  A2(J,1) = TAU2
               ELSE
                  ALPHA2  = A2(J,1)
                  A2(J,1) = ZERO
               END IF
               ALPHA = A1(J,J)
               CALL DLARTG( ALPHA, ALPHA2, C, S, A1(J,J) )
               IF ( K.GT.J )
     $            CALL DROT( K-J, A1(J+1,J), 1, A2(J+1,1), 1, C, S )
               CS(2*K+J*2-1) = C
               CS(2*K+J*2)   = S
  120       CONTINUE
C
            PST2 = 5*K
         ELSE
            PST2 = 2*K
         END IF
C
C        Annihilate B with hyperbolic transformations.
C
         NBL = MIN( NB, Q )
         IF ( NBL.GT.0 ) THEN
C
C           Blocked version.
C
            DO 140  I = 1, K - NBL + 1, NBL
               IB = MIN( K-I+1, NBL )
               CALL DGELQ2( IB, Q, B(I,1), LDB, CS(PST2+I), DWORK,
     $                      IERR )
               IF ( I+IB.LE.K ) THEN
                  CALL DLARFT( 'Forward', 'Rowwise', Q, IB, B(I,1),
     $                         LDB, CS(PST2+I), DWORK, K )
                  CALL DLARFB( 'Right', 'No Transpose', 'Forward',
     $                         'Rowwise', K-I-IB+1, Q, IB, B(I,1),
     $                         LDB, DWORK, K, B(I+IB,1), LDB,
     $                         DWORK( IB+1 ), K )
               END IF
C
C              Annihilate the remaining parts of B.
C
               DO 130 J = I, I + IB - 1
                  IF ( Q.GT.1 ) THEN
                     CALL DLARFG( J-I+1, B(J,1), B(J,2), LDB, TAU2 )
                     ALPHA2 = B(J,1)
                     IF ( K.GT.J ) THEN
                        B(J,1) = ONE
                        CALL DLARF( 'Right', K-J, J-I+1, B(J,1), LDB,
     $                              TAU2, B(J+1,1), LDB, DWORK )
                     END IF
                     B(J,1) = TAU2
                  ELSE
                     ALPHA2 = B(J,1)
                     B(J,1) = ZERO
                  END IF
C
C                 Create hyperbolic rotation.
C
                  CALL MA02FD( A1(J,J), ALPHA2, C, S, IERR )
                  IF ( IERR.NE.0 ) THEN
C
C                    Error return:  The matrix is not positive definite.
C
                     INFO = 1
                     RETURN
                  END IF
C
C                 Apply hyperbolic rotation.
C
                  IF ( K.GT.J ) THEN
                     CALL DSCAL( K-J, ONE/C, A1(J+1,J), 1 )
                     CALL DAXPY( K-J,  -S/C, B(J+1,1), 1, A1(J+1,J), 1 )
                     CALL DSCAL( K-J,     C, B(J+1,1), 1 )
                     CALL DAXPY( K-J,    -S, A1(J+1,J), 1, B(J+1,1), 1 )
                  END IF
                  CS(J*2-1) = C
                  CS(J*2)   = S
  130          CONTINUE
C
  140       CONTINUE
C
         ELSE
            I = 1
         END IF
C
C        Unblocked version for the last or only block.
C
         DO 150  J = I, K
            IF ( Q.GT.1 ) THEN
               CALL DLARFG( Q, B(J,1), B(J,2), LDB, TAU2 )
               ALPHA2 = B(J,1)
               IF ( K.GT.J ) THEN
                  B(J,1) = ONE
                  CALL DLARF( 'Right', K-J, Q, B(J,1), LDB, TAU2,
     $                        B(J+1,1), LDB, DWORK )
               END IF
               B(J,1) = TAU2
            ELSE IF ( Q.GT.0 ) THEN
               ALPHA2 = B(J,1)
               B(J,1) = ZERO
            END IF
            IF ( Q.GT.0 ) THEN
C
C              Create hyperbolic rotation.
C
               CALL MA02FD( A1(J,J), ALPHA2, C, S, IERR )
               IF ( IERR.NE.0 ) THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
C              Apply hyperbolic rotation.
C
               IF ( K.GT.J ) THEN
                  CALL DSCAL( K-J, ONE/C, A1(J+1,J), 1 )
                  CALL DAXPY( K-J,  -S/C, B(J+1,1), 1, A1(J+1,J), 1 )
                  CALL DSCAL( K-J,     C, B(J+1,1), 1 )
                  CALL DAXPY( K-J,    -S, A1(J+1,J), 1, B(J+1,1), 1 )
               END IF
               CS(J*2-1) = C
               CS(J*2)   = S
            END IF
  150    CONTINUE
C
      ELSE
C
C        Row oriented and not deficient generator.
C
         IF ( COL2.GT.0 ) THEN
            NBL = MIN( NB, COL2 )
            IF ( NBL.GT.0 ) THEN
C
C              Blocked version.
C
               DO 170  I = 1, K - NBL + 1, NBL
                  IB = MIN( K-I+1, NBL )
                  CALL DGEQR2( COL2, IB, A2(1,I), LDA2, CS(4*K+I),
     $                         DWORK, IERR )
                  IF ( I+IB.LE.K ) THEN
                     CALL DLARFT( 'Forward', 'Columnwise', COL2, IB,
     $                            A2(1,I), LDA2, CS(4*K+I), DWORK, K )
                     CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                            'Columnwise', COL2, K-I-IB+1, IB,
     $                            A2(1,I), LDA2, DWORK, K, A2(1,I+IB),
     $                            LDA2, DWORK(IB+1), K )
                  END IF
C
C                 Annihilate the remaining parts of A2.
C
                  DO 160 J = I, I + IB - 1
                     IF ( COL2.GT.1 ) THEN
                        LEN = MIN( COL2, J-I+1 )
                        CALL DLARFG( LEN, A2(1,J), A2(2,J), 1, TAU2 )
                        ALPHA2 = A2(1,J)
                        IF ( K.GT.J ) THEN
                           A2(1,J) = ONE
                           CALL DLARF( 'Left', LEN, K-J, A2(1,J), 1,
     $                                 TAU2, A2(1,J+1), LDA2, DWORK )
                        END IF
                        A2(1,J) = TAU2
                     ELSE
                        ALPHA2  = A2(1,J)
                        A2(1,J) = ZERO
                     END IF
                     ALPHA = A1(J,J)
                     CALL DLARTG( ALPHA, ALPHA2, C, S, A1(J,J) )
                     IF ( K.GT.J )
     $                  CALL DROT( K-J, A1(J,J+1), LDA1, A2(1,J+1),
     $                             LDA2, C, S )
                     CS(2*K+J*2-1) = C
                     CS(2*K+J*2)   = S
  160             CONTINUE
C
  170          CONTINUE
C
            ELSE
               I = 1
            END IF
C
C           Unblocked version for the last or only block.
C
            DO 180  J = I, K
               IF ( COL2.GT.1 ) THEN
                  CALL DLARFG( COL2, A2(1,J), A2(2,J), 1, TAU2 )
                  ALPHA2 = A2(1,J)
                  IF ( K.GT.J ) THEN
                     A2(1,J) = ONE
                     CALL DLARF( 'Left', COL2, K-J, A2(1,J), 1, TAU2,
     $                           A2(1,J+1), LDA2, DWORK )
                  END IF
                  A2(1,J) = TAU2
               ELSE
                  ALPHA2  = A2(1,J)
                  A2(1,J) = ZERO
               END IF
               ALPHA = A1(J,J)
               CALL DLARTG( ALPHA, ALPHA2, C, S, A1(J,J) )
               IF ( K.GT.J )
     $            CALL DROT( K-J, A1(J,J+1), LDA1, A2(1,J+1), LDA2, C,
     $                       S )
               CS(2*K+J*2-1) = C
               CS(2*K+J*2)   = S
  180       CONTINUE
C
            PST2 = 5*K
         ELSE
            PST2 = 2*K
         END IF
C
C        Annihilate B with hyperbolic transformations.
C
         NBL = MIN( NB, Q )
         IF ( NBL.GT.0 ) THEN
C
C           Blocked version.
C
            DO 200  I = 1, K - NBL + 1, NBL
               IB = MIN( K-I+1, NBL )
               CALL DGEQR2( Q, IB, B(1,I), LDB, CS(PST2+I), DWORK,
     $                      IERR )
               IF ( I+IB.LE.K ) THEN
                  CALL DLARFT( 'Forward', 'Columnwise', Q, IB, B(1,I),
     $                         LDB, CS(PST2+I), DWORK, K )
                  CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                         'Columnwise', Q, K-I-IB+1, IB, B(1,I),
     $                         LDB, DWORK, K, B(1,I+IB), LDB,
     $                         DWORK( IB+1 ), K )
               END IF
C
C              Annihilate the remaining parts of B.
C
               DO 190  J = I, I + IB - 1
                  IF ( Q.GT.1 ) THEN
                     CALL DLARFG( J-I+1, B(1,J), B(2,J), 1, TAU2 )
                     ALPHA2 = B(1,J)
                     IF ( K.GT.J ) THEN
                        B(1,J) = ONE
                        CALL DLARF( 'Left', J-I+1, K-J, B(1,J), 1,
     $                               TAU2, B(1,J+1), LDB, DWORK )
                     END IF
                     B(1,J) = TAU2
                  ELSE
                     ALPHA2 = B(1,J)
                     B(1,J) = ZERO
                  END IF
C
C                 Create hyperbolic rotation.
C
                  CALL MA02FD( A1(J,J), ALPHA2, C, S, IERR )
                  IF ( IERR.NE.0 ) THEN
C
C                    Error return:  The matrix is not positive definite.
C
                     INFO = 1
                     RETURN
                  END IF
C
C                 Apply hyperbolic rotation.
C
                  IF ( K.GT.J ) THEN
                     CALL DSCAL( K-J, ONE/C, A1(J,J+1), LDA1 )
                     CALL DAXPY( K-J,  -S/C, B(1,J+1), LDB, A1(J,J+1),
     $                           LDA1 )
                     CALL DSCAL( K-J,     C, B(1,J+1), LDB )
                     CALL DAXPY( K-J,    -S, A1(J,J+1), LDA1, B(1,J+1),
     $                           LDB )
                  END IF
                  CS(J*2-1) = C
                  CS(J*2)   = S
  190          CONTINUE
C
  200       CONTINUE
C
         ELSE
            I = 1
         END IF
C
C        Unblocked version for the last or only block.
C
         DO 210  J = I, K
            IF ( Q.GT.1 ) THEN
               CALL DLARFG( Q, B(1,J), B(2,J), 1, TAU2 )
               ALPHA2 = B(1,J)
               IF ( K.GT.J ) THEN
                  B(1,J) = ONE
                  CALL DLARF( 'Left', Q, K-J, B(1,J), 1, TAU2,
     $                        B(1,J+1), LDB, DWORK )
               END IF
               B(1,J) = TAU2
            ELSE IF ( Q.GT.0 ) THEN
               ALPHA2 = B(1,J)
               B(1,J) = ZERO
            END IF
            IF ( Q.GT.0 ) THEN
C
C              Create hyperbolic rotation.
C
               CALL MA02FD( A1(J,J), ALPHA2, C, S, IERR )
               IF ( IERR.NE.0 ) THEN
C
C                 Error return:  The matrix is not positive definite.
C
                  INFO = 1
                  RETURN
               END IF
C
C              Apply hyperbolic rotation.
C
               IF ( K.GT.J ) THEN
                  CALL DSCAL( K-J, ONE/C, A1(J,J+1), LDA1 )
                  CALL DAXPY( K-J,  -S/C, B(1,J+1), LDB, A1(J,J+1), LDA1
     $                                                               )
                  CALL DSCAL( K-J,     C, B(1,J+1), LDB )
                  CALL DAXPY( K-J,    -S, A1(J,J+1), LDA1, B(1,J+1), LDB
     $                                                               )
               END IF
               CS(J*2-1) = C
               CS(J*2)   = S
            END IF
  210    CONTINUE
C
      END IF
C
C *** Last line of MB02CU ***
      END
