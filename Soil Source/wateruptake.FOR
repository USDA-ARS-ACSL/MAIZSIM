c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||c
      Subroutine WaterUptake()
      Include 'public.ins'
      Include 'puplant.ins'
      Integer n,n1,SW10,SW20,SW30,SW40

                       
      t=sngl(Time)
      If (lInput.eq.0) goto 11
      EORSCS=0.0
      AWUPSS=0.0
      Total_Eor=0.0
      potential_T=0.0
      potential_T_EO=0.0
cS  DT adding a temporary output file for debugging the nitrogen uptake code 7/5/07
c     open(100, file='temp.txt')
   
c      write(100,'(1A)')' time        e   sink(e)     RootDenOld  
c     ! RootDenNew RootFluxOld RootFluxNew ErrorFlux'
 
cE
 11    continue
C
c  Routine calculations
C
c     
C
CYAP  No Plant = No Root Activity
CDt setting all water sinks to 0 after this point. if another water 
c   sink source is added later, this will have to be modified. now 
c   we assume all water uptake is only by the plants.
      If((NShoot.eq.0).OR.(LAI.lt.0.0001)) then
        Sink(:)=0.0
        Return 
      endif
CYAPEND
      

C     ------------------------------------------------------------------
C  CALCULATE WATER UPTAKE FROM EACH SOIL CELL FOR PAST PERIOD,
C  AND PREDICT WATER UPTAKE  FOR THE NEXT PERIOD OF THE DAY
C  
C  CALCULATE CURRENT WATER UPTAKE RATE FOR EACH CELL. ***
C
      PSIM=amax1(abs(PSISM),abs(PSIST))
C
      AWUPS = 0.0
      PCRTS = 0.0
      UPH2O = 0.0
      SW20 = 0
      SW30 = 0
c      EORSCF = EOR*SCF
      EORSCF = EOR
C NOTE 10.34
      DO 110, n=1,NumNP
        n1=iFavRoot(n)
        IF ((PSIS(n1).GT.PSILT).OR.(PSIS(n1).LT.PSILT)) THEN
          D20 = (PSIS(n1) - PSIL_)/(PSIS(n1) - PSILT)
        ELSE
          D20 = 0.0
        END IF
        IF (D20.LT.0.0) D20 = 0.0
        IF (SW20.EQ.1.OR.PDWR(n1).LE.0.0) THEN
          AWUP(n1) = (WUPM(n1) + WUPN(n1))*D20
cS 07/05/2007 DT, MK    
          FUP(n1,2) = WUPN(n1)
cE 07/05/2007 DT    

        ELSE

            PCRTS = PCRTS + (PDWR(n1)*CONVR)
          IF (PCRS.GT.PCRTS) THEN
            AWUP(n1) = (WUPM(n1) + WUPN(n1) + WUPT(n1))*D20
cS 07/05/2007 DT MK  
            FUP(n1,2)=WUPN(n1)+WUPT(n1)
cE 07/05/2007 DT
          ELSE
            SW20 = 1
            D21 = ((PDWR(n1)*CONVR)-PCRTS+PCRS)/CONVR   ! length per time

          AWUP(n1) = (WUPM(n1) + WUPN(n1)
     &       + (WUPT(n1)*D21/PDWR(n1)))*D20

cS 07/05/2007 DT MK

            FUP(n1,2)= WUPN(n1) + WUPT(n1)*D21/PDWR(n1)
cE 07/05/2007 DT

          END IF
        END IF
C
cS 07/05/2007 DT  MK  
        FUP(n1,1)=WUPM(n1)*D20
	  FUP(n1,2)=FUP(n1,2)*D20
cE 07/05/2007 DT
        IF (SW30.NE.0) THEN
          AWUP(n1) = 0.
cS 07/05/2007 DT  MK
		FUP(n1,1) = 0.0
		FUP(n1,2) = 0.0
cE 07/05/2007 DT
        ELSE
          AWUPS = AWUPS + AWUP(n1)
          IF (AWUPS.GT.EORSCF) THEN
            SW30 = 1
            AWUP(n1)=AWUP(n1)-AWUPS+EORSCF
cS 07/05/2007 DT   MK
            FUP(n1,1)=FUP(n1,1)-AWUPS+EORSCF
cE 07/05/2007 DT
            AWUPS = EORSCF
cS 07/05/2007 DT   MK
            IF(FUP(n1,1).lt.0) THEN
			FUP(n1,2) = FUP(n1,2)+FUP(n1,1)
			FUP(n1,1) = 0.0
            END IF 
cE 07/05/2007 DT
          END IF
        END IF
 110  CONTINUE


*
      Do n=1,NumNP
        Sink(n)=AWUP(n)/nodeArea(n)
      Enddo

      AWUPSS = AWUPSS + AWUPS*step
      EORSCS = EORSCS + EORSCF*step
      potential_T=potential_T + EOR*step
      potential_T_EO=potential_T_EO+EO*step ! this is calculated from penman and adjusted for plant cover
      Total_Eor=Total_Eor+EOR
c      Total_Pcrs=Total_Pcrs+PCRS/(poprow/100.)*step*2.0
C
cS MK, DT  07/05/2007 CALCULATE FLUX OF WATER UPTAKE BY EXISTING AND NEW ROOTS
        Do 191 n=1,NumNP
         VUP(n,1) = 0.0
	   VUP(n,2) = 0.0
	   FUP(n,1)=AMAX1(FUP(n,1),0.0)/nodeArea(n)
	   FUP(n,2)=AMAX1(FUP(n,2),0.0)/nodeArea(n)
	   IF (RDenM(n).GT.0) VUP(n,1) = FUP(n,1)/RDenM(n) ! water uptake per unit length of root for mature roots
	   IF (RDenY(n).GT.0) VUP(n,2) = FUP(n,2)/RDenY(n)
	   FUP(n,1)=amax1(RDenM(n),0.0)
	   FUP(n,2)=amax1(RDenY(n),0.0)
         IF (sink(n). GT.0.0) THEN	    
		balance=sink(n)-(VUP(n,1)*FUP(n,1)+VUP(n,2)*FUP(n,2))
          balance=sink(n)-(FUP(n,1)+FUP(n,2))
c		write(100,81) time, e, sink(e), 
c     !	FUP(e,1),FUP(e,2),VUP(e,1),VUP(e,2),balance     	
	    CONTINUE
	  ENDIF 	
cE MK 7/05/2007	
191	Continue
81    FORMAT(F14.5,2X,I6,10(2X,E10.3))

      Return
      End
