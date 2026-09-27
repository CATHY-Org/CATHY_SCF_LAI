C
C**************************  ATMBAK ************************************
C
C  special handling of atmospheric boundary conditions for
C  the case where ATMTIM(1) < TIME <= ATMTIM(2) during back-stepping
C
C***********************************************************************
C
      SUBROUTINE ATMBAK(NNOD,TIME,IFATM,AREA,ATMPOT,ATMACT,
     1                  ATMTIM,ATMINP,IETO,DELTAT,SCF,VEG_TYPE)
C
C  SCF-VEG: VEG_TYPE added, SCF widened from scalar to array, indexed
C  by vegetation type, same as in ATMNXT/ATMONE/ETRAN.
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   I,J,IETO
      INTEGER   NNOD
      INTEGER   IFATM(NODMAX)
      INTEGER   VEG_TYPE(*)
      REAL*8    SLOPE,DSATM,ALPHA,SCF(*)
      REAL*8    TIME,GASDEV,DELTAT
      REAL*8    AREA(*),ATMPOT(*),ATMACT(*),ATMTIM(*),ATMINP(3,*)
C
C  we don't need to do anything if first time step or if atmospheric
C  inputs are homogeneous in time
C
      IF (ATMTIM(1) .GE. ATMTIM(2)) GO TO 800
      DO I=1,NNOD
         SLOPE=(ATMINP(2,I) - ATMINP(1,I))/(ATMTIM(2) - ATMTIM(1))
         IF (IETO.NE.0) SLOPE=0.0d0
         ATMPOT(I)=(ATMINP(1,I) + SLOPE*(TIME - ATMTIM(1)))*AREA(I)
         IF (IFATM(I) .EQ. 0) THEN
             IF (ATMPOT(I).GE.0.0d0) THEN
                 ATMACT(I)=ATMPOT(I)
             ELSE
                 ATMACT(I)=(1.0d0-SCF(VEG_TYPE(I)))*ATMPOT(I)
             END IF
         END IF
      END DO
C
  800 RETURN
      END
