      SUBROUTINE XERBLA(SRNAME, INFO)
C
C     SLYCOT
C     Override LAPACK XERBLA routine to raise Python Exception instead of
C     exiting the process
C
CF2PY INTENT(CALLBACK, HIDE) RAISE_XERBLA
C
      CHARACTER*(*) SRNAME
      INTEGER INFO
C
      EXTERNAL RAISE_XERBLA
C
      CALL RAISE_XERBLA(SRNAME, INFO)
C
      RETURN
C
      END
C