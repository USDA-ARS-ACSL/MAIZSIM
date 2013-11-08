*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Subroutine WaterUptake ()
      Include 'public.ins'
      Include 'puplant.ins'
      Integer e,e1,SW10,SW20,SW30,SW40

                       
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
      If((NShoot.eq.0).OR.(LAI<0.0001)) Return
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
      DO 110, e=1,NumEl
        e1=iFavRoot(e)
        IF ((PSIS(e1).GT.PSILT).OR.(PSIS(e1).LT.PSILT)) THEN
          D20 = (PSIS(e1) - PSIL_)/(PSIS(e1) - PSILT)
        ELSE
          D20 = 0.0
        END IF
        IF (D20.LT.0.0) D20 = 0.0
        IF (SW20.EQ.1.OR.PDWR(e1).LE.0.0) THEN
          AWUP(e1) = (WUPM(e1) + WUPN(e1))*D20
cS 07/05/2007 DT, MK    
          FUP(e1,2) = WUPN(e1)
cE 07/05/2007 DT    

        ELSE
          PCRTS = PCRTS + (PDWR(e1)*CONVR)
          IF (PCRS.GT.PCRTS) THEN
            AWUP(e1) = (WUPM(e1) + WUPN(e1) + WUPT(e1))*D20
cS 07/05/2007 DT MK  
            FUP(e1,2)=WUPN(e1)+WUPT(e1)
cE 07/05/2007 DT
          ELSE
            SW20 = 1
            D21 = ((PDWR(e1)*CONVR) - PCRTS + PCRS)/CONVR
            AWUP(e1) = (WUPM(e1) + WUPN(e1)
     &       + (WUPT(e1)*D21/PDWR(e1)))*D20
cS 07/05/2007 DT MK
                FUP(e1,2)= WUPN(e1) + WUPT(e1)*D21/PDWR(e1)
cE 07/05/2007 DT

          END IF
        END IF
C
cS 07/05/2007 DT  MK  
        FUP(e1,1)=WUPM(e1)*D20
	  FUP(e1,2)=FUP(e1,2)*D20
cE 07/05/2007 DT
        IF (SW30.NE.0) THEN
          AWUP(e1) = 0.
cS 07/05/2007 DT  MK
		FUP(e1,1) = 0.0
		FUP(e1,2) = 0.0
cE 07/05/2007 DT
        ELSE
          AWUPS = AWUPS + AWUP(e1)
          IF (AWUPS.GT.EORSCF) THEN
            SW30 = 1
            AWUP(e1) = AWUP(e1) - AWUPS + EORSCF
cS 07/05/2007 DT   MK
		   FUP(e1,1) = FUP(e1,1)- AWUPS + EORSCF
cE 07/05/2007 DT
            AWUPS = EORSCF
cS 07/05/2007 DT   MK
            IF(FUP(e1,1).lt.0) THEN
			FUP(e1,2) = FUP(e1,2)+FUP(e1,1)
			FUP(e1,1) = 0.0
            END IF 
cE 07/05/2007 DT
          END IF
        END IF
 110  CONTINUE


*
      Do e=1,NumEl
        Sink(e)=AWUP(e) / Area(e)
      Enddo

      AWUPSS = AWUPSS + AWUPS*step
      EORSCS = EORSCS + EORSCF*step
      potential_T=potential_T + EOR*step
      potential_T_EO=potential_T_EO+EO*step ! this is calculated from penman and adjusted for plant cover
      Total_Eor=Total_Eor+EOR
c      Total_Pcrs=Total_Pcrs+PCRS/(poprow/100.)*step*2.0
C
cS MK, DT  07/05/2007 CALCULATE FLUX OF WATER UPTAKE BY EXISTING AND NEW ROOTS
        Do 191 e=1,NumEl
         VUP(e,1) = 0.0
	   VUP(e,2) = 0.0
	   FUP(e,1)=AMAX1(FUP(e,1),0.0)/Area(e)
	   FUP(e,2)=AMAX1(FUP(e,2),0.0)
	   IF (RUTDEN(e).GT.0) VUP(e,1) = FUP(e,1)/RUTDEN(e)
	   IF (YRL(e).GT.0) VUP(e,2) = FUP(e,2)/YRL(e)
	   FUP(e,1)=amax1(RUTDEN(e),0.0)
	   FUP(e,2)=amax1(YRL(e)/Area(e),0.0)
         IF (sink(e). GT.0.0) THEN	    
		balance=sink(e)-(VUP(e,1)*FUP(e,1)+VUP(e,2)*FUP(e,2))
c		write(100,81) time, e, sink(e), 
c     !	FUP(e,1),FUP(e,2),VUP(e,1),VUP(e,2),balance     	
	    CONTINUE
	  ENDIF 	
cE MK 7/05/2007	
191	Continue
81    FORMAT(F14.5,2X,I6,10(2X,E10.3))

      Return
      End
