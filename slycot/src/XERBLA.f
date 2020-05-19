      SUBROUTINE XERBLA(SRNAME, INFO)
C
C     SLYCOT
C     Override LAPACK XERBLA routine to noop instead
C     of PRINT and STOP
C
      CHARACTER*(*) SRNAME
      INTEGER INFO
C
      RETURN
C
      END
C
