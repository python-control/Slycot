      SUBROUTINE MC03NX( MP, NP, DP, P, LDP1, LDP2, A, LDA, E, LDE )
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
C     Given an MP-by-NP polynomial matrix of degree dp
C                                    dp-1            dp
C     P(s) = P(0) + ... + P(dp-1) * s     + P(dp) * s            (1)
C
C     the routine composes the related pencil s*E-A where
C
C         | I              |           | O          -P(dp) |
C         |   .            |           | I .           .   |
C     A = |     .          |  and  E = |   . .         .   |.    (2)
C         |       .        |           |     . O       .   |
C         |         I      |           |       I  O -P(2)  |
C         |           P(0) |           |          I -P(1)  |
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
C     MP      (input) INTEGER
C             The number of rows of the polynomial matrix P(s).
C             MP >= 0.
C
C     NP      (input) INTEGER
C             The number of columns of the polynomial matrix P(s).
C             NP >= 0.
C
C     DP      (input) INTEGER
C             The degree of the polynomial matrix P(s).  DP >= 1.
C
C     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1)
C             The leading MP-by-NP-by-(DP+1) part of this array must
C             contain the coefficients of the polynomial matrix P(s)
C             in (1) in increasing powers of s.
C
C     LDP1    INTEGER
C             The leading dimension of array P.  LDP1 >= MAX(1,MP).
C
C     LDP2    INTEGER
C             The second dimension of array P.   LDP2 >= MAX(1,NP).
C
C     A       (output) DOUBLE PRECISION array, dimension
C             (LDA,(DP-1)*MP+NP)
C             The leading DP*MP-by-((DP-1)*MP+NP) part of this array
C             contains the matrix A as described in (2).
C
C     LDA     INTEGER
C             The leading dimension of array A.  LDA >= MAX(1,DP*MP).
C
C     E       (output) DOUBLE PRECISION array, dimension
C             (LDE,(DP-1)*MP+NP)
C             The leading DP*MP-by-((DP-1)*MP+NP) part of this array
C             contains the matrix E as described in (2).
C
C     LDE     INTEGER
C             The leading dimension of array E.  LDE >= MAX(1,DP*MP).
C
C     NUMERICAL ASPECTS
C
C     None.
C
C     CONTRIBUTOR
C
C     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997.
C     Supersedes Release 2.0 routine MC03BX by G.J.H.H. van den Hurk.
C
C     REVISIONS
C
C     -
C
C     KEYWORDS
C
C     Elementary polynomial operations, input output description,
C     polynomial matrix, polynomial operations.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      INTEGER           DP, LDA, LDE, LDP1, LDP2, MP, NP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), E(LDE,*), P(LDP1,LDP2,*)
C     .. Local Scalars ..
      INTEGER           H1, HB, HE, HI, J, K
C     .. External Subroutines ..
      EXTERNAL          DLACPY, DLASET, DSCAL
C     .. Executable Statements ..
C
      IF ( MP.LE.0 .OR. NP.LE.0 )
     $   RETURN
C
C     Initialisation of matrices A and E.
C
      H1 = DP*MP
      HB = H1 - MP
      HE = HB + NP
      CALL DLASET( 'Full', H1, HE, ZERO, ONE, A, LDA )
      CALL DLASET( 'Full', MP, HB, ZERO, ZERO, E, LDE )
      CALL DLACPY( 'Full', HB, HB, A, LDA, E(MP+1,1), LDE )
C
C     Insert the matrices P(0), P(1), ..., P(dp) at the right places
C     in the matrices A and E.
C
      HB = HB + 1
      CALL DLACPY( 'Full', MP, NP, P(1,1,1), LDP1, A(HB,HB), LDA )
      HI = 1
C
      DO 20 K = DP + 1, 2, -1
         CALL DLACPY( 'Full', MP, NP, P(1,1,K), LDP1, E(HI,HB), LDE )
         HI = HI + MP
   20 CONTINUE
C
      DO 40 J = HB, HE
         CALL DSCAL( H1, -ONE, E(1,J), 1 )
   40 CONTINUE
C
      RETURN
C *** Last line of MC03NX ***
      END
