C
C***********************************************************************
C             Subroutine ETRAN, computes root water uptake
C
C  SCF-CAP MODIFICATION: ATMACT (the soil-evaporation share of ATMPOT,
C  already computed by ATMNXT/ATMONE/ATMBAK) is now passed in so that
C  transpiration and evaporation are not withdrawn from the same
C  surface control volume without either mechanism being aware of the
C  other. When SCF is an intermediate value (0 < SCF < 1) and the
C  surface node has already reached the air-dry limit PMIN with an
C  active evaporation demand, the TOP root layer's transpiration
C  weight is withheld and its share is redistributed to deeper,
C  moister layers through the existing BTRAN/GX machinery, instead of
C  imposing a second, uncoordinated withdrawal on top of ATMACT.
C  (PMIN and SCF are already available here via SOILCHAR.H, exactly as
C  they are used, undeclared/unpassed, in ATMONE.F and elsewhere.)
C***********************************************************************
C
      SUBROUTINE ETRAN(N,NNOD,NSTR,ATMPOT,ATMACT,Z,PSI,PNODI,VEG_TYPE,
     1                  QTRANIE)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'SOILCHAR.H'

      INTEGER I,J,K,N,NNOD,NSTR
      INTEGER VEG_TYPE(*)
      REAL*8  SH2O,GX,ZERO,DZ,ZSURF,DEPTH,BETA
      REAL*8  S1,S2,GX1,GX2,TOPWT
      REAL*8  Z(*),PSI(*),PNODI(*)
      REAL*8  ATMPOT(*),ATMACT(*),QTRANIE(*)
      REAL*8  BTRAN(NNOD),BTRANI(N),ETP(NNOD),OMG(NNOD)
      DATA    ZERO/0.0d+00/
   
      CALL INIT0R(N,QTRANIE)
      CALL INIT0R(N,BTRANI)
      DO I = 1,NNOD
         BTRAN(I) = 0.0d0
         OMG(I)=0.0d0
         J = 1
         ZSURF = Z(I)
         DEPTH = 0.0d0
         IF (ATMPOT(I).LT.0.0d0) THEN
            ETP(I) = -1.0d0*SCF(VEG_TYPE(I))*ATMPOT(I)
         ELSE
            ETP(I) = 0.0d0
         END IF
         write(6,*) 'I=',I,' VEG_TYPE=',VEG_TYPE(I)
         write(6,*) 'ZROOT=',ZROOT(VEG_TYPE(I))
         DO WHILE (DEPTH.LE.ZROOT(VEG_TYPE(I)))
            K = (J-1)*NNOD+I
            S1 = PCANA(VEG_TYPE(I))
            S2 = PCANA(VEG_TYPE(I))+1.0D-03
            IF (J.EQ.1) THEN
               DZ = (ZSURF-Z(K+NNOD))/2.0d0
            ELSE
               DZ = (Z(K-NNOD)-Z(K+NNOD))/2.0d0
            END IF
            SH2O = PSI(K)
            GX1 = (SH2O-PCWLT(VEG_TYPE(I))) / 
     &            (PCREF(VEG_TYPE(I))-PCWLT(VEG_TYPE(I)))
            GX1 = MIN(1.0d0,MAX(0.0d0,GX1))
            GX2 = 1.0d0 - (SH2O-S1)/(S2-S1)
            GX2 = MIN(1.0d0,MAX(0.0d0,GX2))
            GX = MIN(GX1,GX2)
CM          IF (SH2O.GE.PCANA) GX = 0.0d0
CM          PZ = -7.044D-03*ZROOT(I)**4 + 1.261D-01*ZROOT(I)**3
CM   +           -8.101D-01*ZROOT(I)**2 + 6.849D+00*ZROOT(I) - 3.345D+00
            BETA = (1-DEPTH/ZROOT(VEG_TYPE(I)))*
     &             DEXP(-1.0d0*PZ(VEG_TYPE(I))*DEPTH/ZROOT(VEG_TYPE(I)))
            BTRANI(K) = MAX(ZERO,BETA*DZ*GX)
            BTRAN(I)  = BTRAN(I) + BETA*DZ
            OMG(I)    = OMG(I) + GX*BETA*DZ
            J = J + 1
            IF ((J-1)*NNOD+I.LE.NMAX) THEN
               DEPTH = ZSURF - Z((J-1)*NNOD+I)
            ELSE
               write(6,*)'DECREASE ZROOT, TOO CLOSE TO THE BOTTOM'
               stop
            END IF
cm          write(111,*)k,gx1,gx2,beta,btrani(k)
         END DO
C
C  SCF-CAP: node I's top root layer sits in the same surface control
C  volume as the ATMACT boundary node I. If ATMACT is already pulling
C  water there (evaporation demand present, i.e. SCF < 1) AND the node
C  has already reached the air-dry limit PMIN, don't also let
C  transpiration draw on that same, already-depleted top layer.
C  Withhold its weight here; BTRAN/OMG below then route the demand to
C  deeper layers instead of stacking a second withdrawal on PMIN-
C  limited storage.
C
         IF (ATMACT(I).LT.0.0D0 .AND. PSI(I).LE.PMIN) THEN
            TOPWT     = BTRANI(I)
            BTRAN(I)  = BTRAN(I) - TOPWT
            OMG(I)    = OMG(I)   - TOPWT
            BTRANI(I) = 0.0D0
         END IF
         write(6,*) 'I=',I,' BTRAN=',BTRAN(I),' OMG=',OMG(I),
     &               ' PSI=',PSI(I)
         BTRAN(I) = MAX(ZERO,BTRAN(I))
         IF (BTRAN(I).GT.0.0D0) THEN
            OMG(I) = OMG(I)/BTRAN(I)
         ELSE
C  SCF-CAP: no root-zone capacity left (e.g. a single-layer root zone
C  whose only layer was just withheld above); avoid a divide-by-zero
C  and let QTRANIE fall back to 0 for this node below.
            OMG(I) = 0.0D0
         END IF
      END DO
      DO I=1,NNOD
cm       write(666,*)i,etp(i),btran(i),zroot(i)
         IF (BTRAN(I).GT.0.0D0) THEN
            DO J=1,NSTR+1
               K = (J-1)*NNOD+I
               QTRANIE(K)=ETP(I)*BTRANI(K)/BTRAN(I)/
     &                    MAX(OMG(I),OMGC(VEG_TYPE(I)))
cm             write(111,*)k,qtranie(k),btrani(k),btran(i)
            END DO
         ELSE
            DO J=1,NSTR+1
               K = (J-1)*NNOD+I
               QTRANIE(K)=0.0D0
            END DO
         END IF
      END DO

      RETURN
      END
