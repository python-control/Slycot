      SUBROUTINE MC03NY( NBLCKS, NRA, NCA, A, LDA, E, LDE, IMUK, INUK,
     $                   VEPS, LDVEPS, INFO )
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
C     To determine a minimal basis of the right nullspace of the
C     subpencil s*E(eps)-A(eps) using the method given in [1] (see
C     Eqs.(4.6.8), (4.6.9)).
C     This pencil only contains Kronecker column indices, and it must be
C     in staircase form as supplied by SLICOT Library Routine MB04VD.
C     The basis vectors are represented by matrix V(s) having the form
C
C                | V11(s) V12(s) V13(s)   . .   V1n(s) |
C                |        V22(s) V23(s)         V2n(s) |
C                |               V33(s)           .    |
C         V(s) = |                  .             .    |
C                |                      .         .    |
C                |                          .     .    |
C                |                              Vnn(s) |
C
C     where n is the number of full row rank blocks in matrix A(eps) and
C
C                                               k               j-i
C         Vij(s) = Vij,0 + Vij,1*s +...+ Vij,k*s +...+ Vij,j-i*s   . (1)
C
C     In other words, Vij,k is the coefficient corresponding to degree k
C     in the matrix polynomial Vij(s).
C     Vij,k has dimensions mu(i)-by-(mu(j)-nu(j)).
C     The coefficients Vij,k are stored in the matrix VEPS as follows
C     (for the case n = 3):
C
C         sizes      m1-n1    m2-n2   m2-n2    m3-n3   m3-n3   m3-n3
C
C             m1 { | V11,0 || V12,0 | V12,1 || V13,0 | V13,1 | V13,2 ||
C                  |       ||       |       ||       |       |       ||
C      VEPS = m2 { |       || V22,0 |       || V23,0 | V23,1 |       ||
C                  |       ||       |       ||       |       |       ||
C             m3 { |       ||       |       || V33,0 |       |       ||
C
C     where mi = mu(i), ni = nu(i).
C     Matrix VEPS has dimensions nrv-by-ncv where
C       nrv = Sum(i=1,...,n) mu(i)
C       ncv = Sum(i=1,...,n) i*(mu(i)-nu(i))
C
C     ==================================================================
C     REMARK: This routine is intended to be called only from the SLICOT
C             routine MC03ND.
C     ==================================================================
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     NBLCKS  (input) INTEGER
C             Number of full row rank blocks in subpencil
C             s*E(eps)-A(eps) that contains all Kronecker column indices
C             of s*E-A.  NBLCKS >= 0.
C
C     NRA     (input) INTEGER
C             Number of rows of the subpencil s*E(eps)-A(eps) in s*E-A.
C             NRA = nu(1) + nu(2) + ... + nu(NBLCKS).  NRA >= 0.
C
C     NCA     (input) INTEGER
C             Number of columns of the subpencil s*E(eps)-A(eps) in
C             s*E-A.
C             NCA = mu(1) + mu(2) + ... + mu(NBLCKS).  NCA >= 0.
C
C     A       (input/output) DOUBLE PRECISION array, dimension (LDA,NCA)
C     E       (input/output) DOUBLE PRECISION array, dimension (LDE,NCA)
C             On entry, the leading NRA-by-NCA part of these arrays must
C             contain the matrices A and E, where s*E-A is the
C             transformed pencil s*E0-A0 which is the pencil associated
C             with P(s) as described in [1] Section 4.6. The pencil
C             s*E-A is assumed to be in generalized Schur form.
C             On exit, these arrays contain no useful information.
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,NRA).
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,NRA).
C
C     IMUK    (input) INTEGER array, dimension (NBLCKS)
C             This array must contain the column dimensions mu(k) of the
C             full column rank blocks in the subpencil s*E(eps)-A(eps)
C             of s*E-A. The content of IMUK is modified by the routine
C             but restored on exit.
C
C     INUK    (input) INTEGER array, dimension (NBLCKS)
C             This array must contain the row dimensions nu(k) of the
C             full row rank blocks in the subpencil s*E(eps)-A(eps) of
C             s*E-A.
C
C     VEPS    (output) DOUBLE PRECISION array, dimension (LDVEPS,ncv)
C             Let nrv = Sum(i=1,...,NBLCKS) mu(i) = NCA,
C                 ncv = Sum(i=1,...,NBLCKS) i*(mu(i)-nu(i)).
C             The leading nrv-by-ncv part of this array contains the
C             column vectors of a minimal polynomial basis for the right
C             nullspace of the subpencil s*E(eps)-A(eps). (See [1]
C             Section 4.6.4.) An upper bound for ncv is (NRA+1)*NCA.
C
C     LDVEPS  INTEGER
C             The leading dimension of array VEPS.
C             LDVEPS >= MAX(1,NCA).
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             > 0:  if INFO = k, the k-th diagonal block of A had not a
C                   full row rank.
C
C     REFERENCES
C
C     [1] Th.G.J. Beelen, New Algorithms for Computing the Kronecker
C         structure of a Pencil with Applications to Systems and
C         Control Theory.
C         Ph.D.Thesis, Eindhoven University of Technology, 1987.
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTORS
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC03BY by Th.G.J. Beelen,
C     A.J. Geurts, and G.J.H.H. van den Hurk.
C
C     REVISIONS
C
C     Dec. 1997.
C
C     KEYWORDS
C
C     Elementary polynomial operations, Kronecker form, polynomial
C     matrix, polynomial operations, staircase form.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDE, LDVEPS, NBLCKS, NCA, NRA
C     .. Array Arguments ..
      INTEGER           IMUK(*), INUK(*)
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), VEPS(LDVEPS,*)
C     .. Local Scalars ..
      INTEGER           AC1, AC2, AR1, ARI, ARK, DIF, EC1, ER1, I, J, K,
     $                  MUI, NCV, NRV, NUI, SMUI, SMUI1, VC1, VC2, VR1,
     $                  VR2, WC1, WR1
C     .. Local Arrays ..
      DOUBLE PRECISION  DUMMY(1)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMM, DLASET, DSCAL, DTRTRS, XERBLA
C     .. Executable Statements ..
C
      INFO = 0
      IF( NBLCKS.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRA.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCA.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, NRA ) ) THEN
         INFO = -5
      ELSE IF( LDE.LT.MAX( 1, NRA ) ) THEN
         INFO = -7
      ELSE IF( LDVEPS.LT.MAX( 1, NCA ) ) THEN
         INFO = -11
      END IF
C
      IF ( INFO.NE.0 ) THEN
C
C        Error return.
C
         CALL XERBLA( 'MC03NY', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ( NBLCKS.EQ.0 .OR. NRA.EQ.0 .OR. NCA.EQ.0 )
     $   RETURN
C
C     Computation of the nonzero parts of W1 and W2:
C
C          | AH11 AH12 ... AH1n |       | EH11 EH12 ... EH1n |
C          |      AH22     AH2n |       |      EH22     EH2n |
C     W1 = |         .       .  |, W2 = |         .       .  |
C          |           .     .  |       |           .     .  |
C          |               AHnn |       |               EHnn |
C
C     with AHij = -pinv(Aii) * Aij, EHij = pinv(Aii) * Eij and EHii = 0,
C     AHij and EHij have dimensions mu(i)-by-mu(j), Aii = [ Oi | Ri ],
C     and
C       Ri is a regular nu(i)-by-nu(i) upper triangular matrix;
C       Oi is a not necessarily square null matrix.
C     Note that the first mu(i)-nu(i) rows in AHij and EHij are zero.
C     For memory savings, the nonzero parts of W1 and W2 are constructed
C     over A and E, respectively.
C
C     (AR1,AC1) denotes the position of the first element of the
C     submatrix Ri in matrix Aii.
C     EC1 is the index of the first column of Ai,i+1/Ei,i+1.
C
      EC1 = 1
      AR1 = 1
C
      DO 40 I = 1, NBLCKS - 1
         NUI = INUK(I)
         IF ( NUI.EQ.0 ) GO TO 60
         MUI = IMUK(I)
         EC1 = EC1 + MUI
         AC1 = EC1 - NUI
         CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', NUI,
     $                NCA-EC1+1, A(AR1,AC1), LDA, E(AR1,EC1), LDE,
     $                INFO )
         IF ( INFO.GT.0 ) THEN
            INFO = I
            RETURN
         END IF
C
         DO 20 J = 1, NUI
            CALL DSCAL( J, -ONE, A(AR1,AC1+J-1), 1 )
   20    CONTINUE
C
         CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', NUI,
     $                NCA-EC1+1, A(AR1,AC1), LDA, A(AR1,EC1), LDA,
     $                INFO )
         AR1 = AR1 + NUI
   40 CONTINUE
C
   60 CONTINUE
C
C     The contents of the array IMUK is changed for temporary use in
C     this routine as follows:
C
C        IMUK(i) = Sum(j=1,...,i) mu(j).
C
C     On return, the original contents of IMUK is restored.
C     In the same loop the actual number of columns of VEPS is computed.
C     The number of rows of VEPS is NCA.
C
C        NRV = Sum(i=1,...,NBLCKS) mu(i) = NCA,
C        NCV = Sum(i=1,...,NBLCKS) i*(mu(i)-nu(i)).
C
      SMUI = 0
      NCV = 0
C
      DO 80 I = 1, NBLCKS
         MUI = IMUK(I)
         SMUI = SMUI + MUI
         IMUK(I) = SMUI
         NCV = NCV + I*( MUI - INUK(I) )
   80 CONTINUE
C
      NRV = NCA
C
C     Computation of the matrix VEPS.
C
C     Initialisation of VEPS to zero.
C
      CALL DLASET( 'Full', NRV, NCV, ZERO, ZERO, VEPS, LDVEPS )
C                                                           | I |
C     Set Vii,0 = Kii in VEPS , i=1,...,NBLCKS, where Kii = |---|
C                                                           | O |
C     and I is an identity matrix of size mu(i)-nu(i),
C         O is a null matrix, dimensions nu(i)-by-(mu(i)-nu(i)).
C
C     WR1 := Sum(j=1,...,i-1) mu(j) + 1
C            is the index of the first row in Vii,0 in VEPS.
C     WC1 := Sum(j=1,...,i-1) j*(mu(j)-nu(j)) + 1
C            is the index of the first column in Vii,0 in VEPS.
C
      DUMMY(1) = ONE
      NUI = IMUK(1) - INUK(1)
      CALL DCOPY( NUI, DUMMY, 0, VEPS, LDVEPS+1 )
      WR1 = IMUK(1) + 1
      WC1 = NUI + 1
C
      DO 100 I = 2, NBLCKS
         NUI = IMUK(I) - IMUK(I-1) - INUK(I)
         CALL DCOPY( NUI, DUMMY, 0, VEPS(WR1,WC1), LDVEPS+1 )
         WR1 = IMUK(I) + 1
         WC1 = WC1 + I*NUI
  100 CONTINUE
C
C     Determination of the remaining nontrivial matrices in Vij,k
C     block column by block column with decreasing block row index.
C
C     The computation starts with the second block column since V11,0
C     has already been determined.
C     The coefficients Vij,k satisfy the recurrence relation:
C
C        Vij,k = Sum(r=i+1,...,j-k)   AHir*Vrj,k +
C              + Sum(r=i+1,...,j-k+1) EHir*Vrj,k-1,   i + k < j,
C
C              = EHi,i+1 * Vi+1,j,k-1                 i + k = j.
C
C     This recurrence relation can be derived from [1], (4.6.8)
C     and formula (1) in Section PURPOSE.
C
      VC1 = IMUK(1) - INUK(1) + 1
      ARI = 1
C
      DO 180 J = 2, NBLCKS
         DIF = IMUK(J) - IMUK(J-1) - INUK(J)
         ARI = ARI + INUK(J-1)
         ARK = ARI
C
C        Computation of the matrices Vij,k where i + k < j.
C        Each matrix Vij,k has dimension mu(i)-by-(mu(j) - nu(j)).
C
         DO 160 K = 0, J - 2
C
C           VC1, VC2 are the first and last column index of Vij,k.
C
            VC2 = VC1 + DIF - 1
            AC2 = IMUK(J-K)
            AR1 = ARK
            ARK = ARK - INUK(J-K-1)
C
            DO 120 I = J - K - 1, 1, -1
C
C              Compute the first part of Vij,k in decreasing order:
C              Vij,k := Vij,k + Sum(r=i+1,..,j-k) AHir*Vrj,k.
C              The non-zero parts of AHir are stored in
C              A(AR1:AR1+nu(i)-1,AC1:AC2) and Vrj,k are stored in
C              VEPS(AC1:AC2,VC1:VC2).
C              The non-zero part of the result is stored in
C              VEPS(VR1:VR2,VC1:VC2).
C
               VR2 = IMUK(I)
               AC1 = VR2 + 1
               VR1 = AC1 - INUK(I)
               AR1 = AR1 - INUK(I)
               CALL DGEMM( 'No transpose', 'No transpose', INUK(I),
     $                     DIF, AC2-VR2, ONE, A(AR1,AC1), LDA,
     $                     VEPS(AC1,VC1), LDVEPS, ONE, VEPS(VR1,VC1),
     $                     LDVEPS )
  120       CONTINUE
C
            ER1 = 1
C
            DO 140 I = 1, J - K - 1
C
C              Compute the second part of Vij,k+1 in normal order:
C              Vij,k+1 := Sum(r=i+1,..,j-k) EHir*Vrj,k.
C              The non-zero parts of EHir are stored in
C              E(ER1:ER1+nu(i)-1,EC1:AC2) and Vrj,k are stored in
C              VEPS(EC1:AC2,VC1:VC2).
C              The non-zero part of the result is stored in
C              VEPS(VR1:VR2,VC2+1:VC2+DIF), where
C              DIF = VC2 - VC1 + 1 = mu(j) - nu(j).
C              This code portion also computes Vij,k+1 for i + k = j.
C
               VR2 = IMUK(I)
               EC1 = VR2 + 1
               VR1 = EC1 - INUK(I)
               CALL DGEMM( 'No transpose', 'No transpose', INUK(I),
     $                     DIF, AC2-VR2, ONE, E(ER1,EC1), LDE,
     $                     VEPS(EC1,VC1), LDVEPS, ZERO, VEPS(VR1,VC2+1),
     $                     LDVEPS )
               ER1 = ER1 + INUK(I)
  140       CONTINUE
C
            VC1 = VC2 + 1
  160    CONTINUE
C
         VC1 = VC1 + DIF
  180 CONTINUE
C
C     Restore original contents of the array IMUK.
C
C     Since, at the moment:
C       IMUK(i) = Sum(j=1,...,i) mu(j),   (i=1,...,NBLCKS),
C     the original values are:
C       mu(i) = IMUK(i) - IMUK(i-1)  with IMUK(0 ) = 0.
C
      SMUI1 = 0
C
      DO 200 I = 1, NBLCKS
         SMUI = IMUK(I)
         IMUK(I) = SMUI - SMUI1
         SMUI1 = SMUI
  200 CONTINUE
C
      RETURN
C *** Last line of MC03NY ***
      END
