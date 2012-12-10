C     INPUTS: F TS-DENSITY; X-OBSERVATIONS; N SAMPLE SIZE
C
C     OUTPUTS: ICOMAX -maximum number of observations in interval
C              IRMAX- maximum number of consecutive intervals with no
C              observation
C

      SUBROUTINE DENLOCAL(CF,KNI,N,ICOMAX,IRMAX)
      INTEGER N,ICOMAX,IRMAX
      DOUBLE PRECISION CF(N)
      INTEGER KNI(N)
C
      INTEGER I,J,ICR,IC
C
C
      I=1
      J=1
      ICR=0
 20   CONTINUE
      IC=0

 30   CONTINUE
      IF(I.GT.N) GOTO 31
C
C     CALCULATE INTERVAL OCCUPANCIES IC
C
      IF(CF(I).LE.DBLE(J)/DBLE(N)) THEN
C
C     DETERMINE INTERVAL TO BE SQUEEZED FOR RUN LENGTH
C
         IF(ICR.GE.IRMAX) THEN
            KNI(I)=1
            KNI(I-1)=1
         ENDIF
         ICR=0
         IC=IC+1
         IF(I.LT.N) THEN
            I=I+1
            GOTO 30
         ENDIF
      ENDIF
 31   CONTINUE
C
C     CALCULATE MAXIMUM RUN LENGTH ICR
C
      IF(IC.EQ.0) ICR=ICR+1
      IF(IC.GE.ICOMAX) THEN
C
C     DETERMINE INTERVAL TO BE SQUEEZED FOR OCCUPANCY
C
         DO 35 LL=I-IC+1,I
            KNI(LL)=1
 35      CONTINUE
      ENDIF
      IF(I.EQ.N) RETURN
      J=J+1
      GOTO 20
C
      END
C

C
