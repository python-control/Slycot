      INTEGER FUNCTION UE01MD( ISPEC, NAME, OPTS, N1, N2, N3 )
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
C     To provide an extension of the LAPACK routine ILAENV to
C     machine-specific parameters for SLICOT routines.
C
C     The default values in this version aim to give good performance on
C     a wide range of computers. For optimal performance, however, the
C     user is advised to modify this routine. Note that an optimized
C     BLAS is a crucial prerequisite for any speed gains. For further
C     details, see ILAENV.
C
C     FUNCTION VALUE
C
C     UE01MD  INTEGER
C             The function value set according to ISPEC.
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     ISPEC   (input) INTEGER
C             Specifies the parameter to be returned as the value of
C             UE01MD, as follows:
C             = 1: the optimal blocksize; if the returned value is 1, an
C                  unblocked algorithm will give the best performance;
C             = 2: the minimum block size for which the block routine
C                  should be used; if the usable block size is less than
C                  this value, an unblocked routine should be used;
C             = 3: the crossover point (in a block routine, for N less
C                  than this value, an unblocked routine should be used)
C             = 4: the number of shifts, used in the product eigenvalue
C                  routine;
C             = 8: the crossover point for the multishift QR method for
C                  product eigenvalue problems.
C
C     NAME    (input) CHARACTER*(*)
C             The name of the calling subroutine, in either upper case
C             or lower case.
C
C     OPTS    (input) CHARACTER*(*)
C             The character options to the subroutine NAME, concatenated
C             into a single character string.
C
C     N1      (input) INTEGER
C     N2      (input) INTEGER
C     N3      (input) INTEGER
C             Problem dimensions for the subroutine NAME; these may not
C             all be required.
C
C     CONTRIBUTORS
C
C     D. Kressner, Technical Univ. Berlin, Germany, and
C     P. Benner, Technical Univ. Chemnitz, Germany, December 2003.
C
C     REVISIONS
C
C     V. Sima, June 2008 (SLICOT version of the HAPACK routine ILAHAP).
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3
C
C     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1, C3
      CHARACTER*2        C2
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, MAX
C
C     .. Executable Statements ..
C
      IF ( ISPEC.EQ.1 .OR. ISPEC.EQ.2 .OR. ISPEC.EQ.3 ) THEN
C
C        Convert NAME to upper case if the first character is lower
C        case.
C
         UE01MD = 1
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1:1 ) )
         IZ = ICHAR( 'Z' )
         IF ( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
C
C           ASCII character set.
C
            IF ( IC.GE.97 .AND. IC.LE.122 ) THEN
               SUBNAM( 1:1 ) = CHAR( IC-32 )
               DO 10 I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF( IC.GE.97 .AND. IC.LE.122 )
     $               SUBNAM( I:I ) = CHAR( IC-32 )
   10          CONTINUE
            END IF
C
         ELSE IF ( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
C
C           EBCDIC character set.
C
            IF ( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $           ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $           ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
               SUBNAM( 1:1 ) = CHAR( IC+64 )
               DO 20 I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF ( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $                 ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $                 ( IC.GE.162 .AND. IC.LE.169 ) )
     $               SUBNAM( I:I ) = CHAR( IC+64 )
   20          CONTINUE
            END IF
C
         ELSE IF ( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
C
C           Prime machines:  ASCII+128.
C
            IF ( IC.GE.225 .AND. IC.LE.250 ) THEN
               SUBNAM( 1:1 ) = CHAR( IC-32 )
               DO 30 I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF ( IC.GE.225 .AND. IC.LE.250 )
     $               SUBNAM( I:I ) = CHAR( IC-32 )
   30          CONTINUE
            END IF
         END IF
C
         C1 = SUBNAM( 1:1 )
         SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
         CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
         IF ( .NOT.( CNAME .OR. SNAME ) )
     $      RETURN
         C2 = SUBNAM( 4:5 )
         C3 = SUBNAM( 6:6 )
C
         IF ( ISPEC.EQ.1 ) THEN
C
C           Block size.
C
            NB = 1
            IF ( C2.EQ.'4S' .OR. C2.EQ.'4T' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NB = ILAENV( 1, 'DGEQRF', ' ', N1, N2, -1, -1 ) / 2
               ELSE IF ( C3.EQ.'T' ) THEN
                  NB = ILAENV( 1, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 4
               END IF
            ELSE IF ( C2.EQ.'4P' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NB = ILAENV( 1, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2
               END IF
            ELSE IF ( C2.EQ.'4W' .OR. C2.EQ.'4Q' ) THEN
               IF ( C3.EQ.'D' ) THEN
                  NB = ILAENV( 1, 'DORGQR', ' ', N1, N2, N3, -1 ) / 2
               ELSE IF ( C3.EQ.'B' ) THEN
                  NB = ILAENV( 1, 'DORMQR', ' ', N1, N2, N3, -1 ) / 2
               END IF
**          ELSE IF ( C2.EQ.'SH' ) THEN
**             IF ( C3.EQ.'PVB' ) THEN
**                NB = ILAENV( 1, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2
**             END IF
            END IF
            UE01MD = NB
         ELSE IF ( ISPEC.EQ.2 ) THEN
C
C           Minimum block size.
C
            NBMIN = 2
            IF ( C2.EQ.'4S' .OR. C2.EQ.'4T' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', N1, N2, -1,
     $                                    -1 ) / 2 )
               ELSE IF ( C3.EQ.'T' ) THEN
                  NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N1, N2, N1,
     $                                    -1 ) / 4 )
               END IF
            ELSE IF ( C2.EQ.'4P' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N1, N2, N1,
     $                                    -1 ) / 4 )
               END IF
            ELSE IF ( C2.EQ.'4W' .OR. C2.EQ.'4Q' ) THEN
               IF ( C3.EQ.'D' ) THEN
                  NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', N1, N2, N3,
     $                                    -1 ) / 2 )
               ELSE IF ( C3.EQ.'B' ) THEN
                  NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', ' ', N1, N2, N3,
     $                                    -1 ) / 2 )
               END IF
**          ELSE IF ( C2.EQ.'SH' ) THEN
**             IF ( C3.EQ.'PVB' ) THEN
**                NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N1, N2, N1,
**   $                                    -1 ) / 4 )
**             END IF
            END IF
            UE01MD = NBMIN
         ELSE IF ( ISPEC.EQ.3 ) THEN
C
C           Crossover point.
C
            NX = 0
            IF ( C2.EQ.'4S' .OR. C2.EQ.'4T' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NX = ILAENV( 3, 'DGEQRF', ' ', N1, N2, -1, -1 )
               ELSE IF ( C3.EQ.'T' ) THEN
                  NX = ILAENV( 3, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2
               END IF
            ELSE IF ( C2.EQ.'4P' ) THEN
               IF ( C3.EQ.'B' ) THEN
                  NX = ILAENV( 3, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2
               END IF
            ELSE IF ( C2.EQ.'4W' .OR. C2.EQ.'4Q' ) THEN
               IF ( C3.EQ.'D' ) THEN
                  NX = ILAENV( 3, 'DORGQR', ' ', N1, N2, N3, -1 )
               ELSE IF ( C3.EQ.'B' ) THEN
                  NX = ILAENV( 3, 'DORGQR', ' ', N1, N2, N3, -1 )
               END IF
**          ELSE IF ( C2.EQ.'SH' ) THEN
**             IF ( C3.EQ.'PVB' ) THEN
**                NX = ILAENV( 3, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2
**             END IF
            END IF
            UE01MD = NX
         END IF
      ELSE IF ( ISPEC.EQ.4 ) THEN
C
C        Number of shifts (used by MB03XP).
C
         UE01MD = ILAENV( 4, 'DHSEQR', OPTS, N1, N2, N3, -1 )
      ELSE IF ( ISPEC.EQ.8 ) THEN
C
C        Crossover point for multishift (used by MB03XP).
C
         UE01MD = ILAENV( 8, 'DHSEQR', OPTS, N1, N2, N3, -1 )
      ELSE
C
C        Invalid value for ISPEC.
C
         UE01MD = -1
      END IF
      RETURN
C *** Last line of UE01MD ***
      END
