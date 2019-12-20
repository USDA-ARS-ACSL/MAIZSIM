c pot_ET is used to save the accumulated potential transpiration in the previous N steps
c p_ET is used to save the difference between the accumulated transpiration in the 
c previous N step and N+1 step; thus in fact it is the potential transpiration rate per time step
      subroutine Output ()
      Include 'public.ins'
      Include 'PuWeath.ins'
	  Include 'PuPlant.ins'

      Integer Plevel,e
      
      Dimension temp(1000)
      Real PotentialSoilEvaporation, WidthT,Period
      Real SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil
      Character*10 date, Stime
      Common /Outcom/NumPoint(2000),NodNum(2000,60),
     &        NumCell(1000),NumElPrint(1000,60),ModNum,
     &        Period, ETO, ActualTrans1,ActualSoilEvap,
     &        pot_ET,pot_ET_EO,flux, ActualInfil,timer,
     &        FluxAct,Runoff,PotentialSoilEvaporation,
     &        CumRain,
     &        SeasonRainfall, SeasonPotenTranspiration,
     &        SeasonActualTranspiration, SeasonPotenSoilEvap,
     &        SeasonActualSoilEvap,SeasonActualInfil

      t=Time
      
      If(lInput.eq.1) then
C
C  Initialize 
C
        If (DailyOutput.eq.1) Period=1.0
        If (HourlyOutput.eq.1) Period=1.0/(24.0)
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum) = time

C*** FROM HourlyOut.for BEGIN      
        timer=0
        flux=0.0
        fluxAct=0.0
        Runoff=0.0
        ActualSoilEvap=0
        ActualInfil=0
        PotentialSoilEvaporation=0.0
        CumRain=0.0
         SeasonRainfall =0.0
         SeasonActualInfil=0.0
         SeasonPotenTranspiration=0.0
        SeasonActualTranspiration=0.0
        SeasonPotenSoilEvap=0.0
        SeasonActualSoilEvap=0.0
	  ActualTrans1=0.0
	  pot_ET=0.0
	  pot_ET_EO=0.0
	  Open(81,file=NodeGraphics)
	  open(82,file=ElemGraphics)  
	  Open(83,file=SurfaceGraphics)
        Open(84,file=FluxGraphics)
         
      Write(81,5) "Date_time,", "Date,", "X,", "Y,", "hNew,","thNew,",
     &  "VZ,", "VX,","Q,", "ConcN,", "Temp,","GasConc"  
     
      Write(82,6) "Date_time,", "Date,", "X,", "Y,", "Node,", 
     &      "RMassM,",
     &      "RMassY,", "RDenM,", "RDenY,",
     &      "WaterSink,", "NitSink,", "GasSink,", "Area"
         
        Write(83,20) "Date_time,", "Date,", "PSoilEvap,", 
     &         "ASoilEVap,", "PET_PEN,",
     &         "PE_T_int,","transp,",
     &         "CumRain,","infil,",
     &         "FLuxAct,", "Runoff,",
     &         "cover,","PSIM,", "SeasPSoEv,",
     $         "SeasASoEv," ,"SeasPTran,",
     &         "SeasATran,","SeasRain,","SeasInfil" 
       Write(84,20) "Date_time,", "Date,","Flux,","SoilEA,",
     &    "WattSM,", "RNS,","RNC,","WATPOT,","WATACT,",
     &    "WATRAT,", "UpwardLW,","DifWat,","Cloud,",
     &     "Shadow,","DIfInt,","Rowinc,", "PARInt,","RADInt"
      Endif     !end initialization
C
C Routine calculations every time step
C

C Potential values:

*  add only if surface is not shaded
* note that VarBW(i,3) is the adjusted et or rainfall - i.e., sum of rain + et, not separated
* we may have to consider changing this in the future
* 6/2016 DT changed CumRain sum to use varBW(i,1) which is actual rain.
* Actual infiltration, after accounting for evap is calculated in ActualInfil
          Do i=1,NumBp
	          n=KXB(i)
             If (VarBW(i,3).GE.0) Then
      	        PotentialSoilEvaporation=PotentialSoilEvaporation+
     & VarBW(i,3)*width(i)*step
	         EndIf
	           CumRain=CumRain+varBW(i,1)*width(i)*step
	       End Do
	
C Actual values

CDT note that flow takes place in seepage faces when the boundary code is positive - this indicates a fixed
C  pressure head and flux is calculated CodeW >0 indicates head is set at 0 and flux is calculated.
         timer=timer+step
      	 do i=1,NumBp
	     n=KXB(i)
	     
	     
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))
     &             .or.(CodeW(n).eq.(1))) then
	        flux=flux+Q(n)*step
	        if (flux.gt.0) then 
	          iii=1
	          endif
	        else
	         if (CodeW(n).EQ.-2) then
	       	         iii=1
	       	  endif
           endif
c case for evaporation
           if ((abs(CodeW(n)).eq.4).and.(VarBw(i,3).GT.0)) then
              ActualSoilEvap=ActualSoilEvap-min(0.0,Q(n))*step
            endif
c case for rainfall
           if (abs(CodeW(n)).eq.4)
     !            then
              ActualInfil=ActualInfil+max(0.0,Q(n))*step
              Runoff=Runoff+RO(n)*step
              RO(n)=0.0   ! we will use this to store ponded water for next time period
              FluxAct=FluxAct+max(0.0,QAct(n))*step
            endif

	  EndDo
	  
	
C*** FROM HourlyOut.for END

11    If(abs(time-tNext(ModNum)).lt.0.001*Step.OR.lInput.eq.1) then
          iday=int(t)
        call caldat(iday,mm,id,iyyy) 
        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          Call ArrNodOut (Date,x,y,NumNP,NumEl,t,NumElD,hNew,
     !                   ThNew,vz,vx,Q,Conc(1,1),Tmpr,g(1,1))

        If(NShoot.ne.0) then
             Call ArrRootOut (Date,x,y,t, NumNP,RMassM,RMassY,
     !            RDenM, RDenY, Sink,
     !            cSink(1,1),gSink(1,1),NodeArea)
            !Call OutToDB(Sink(1))
        Endif
        tNext(ModNum) = time + Period
C      
* This code operates hourly to output data on the surface nodes
	  
	  WidthT=0.0
        Do i= 1, NumBP
           if (abs(CodeW(KXB(i))).eq.4) then
           	  WidthT=WidthT+width(i)
           	 endif
        Enddo


* Anything calculated from VarBW should not be divided by 24 
* output as mm per plant, essentially dividing by popslab

              PotentialSoilEvaporation=PotentialSoilEvaporation
     &        /EOMult/poprow*100
              ActualSoilEvap=ActualSoilEvap/EOMult/poprow*100
              ActualInfil=ActualInfil/EOMult/poprow*100
              FluxAct=FluxAct/EOMult/poprow*100
              Runoff=Runoff/EOMult/poprow*100
              CumRain=CumRain/EOMult/poprow*100
              Flux=Flux/EOMult/poprow*100
              


CDT distribute runoff back into Q This makes it uniform. implement this later

!c         do n=1, NumBP
!c           ROForInfil(n)=Runoff*Width(n)
!c          end do
          
* dt - value of awupss and pot_et are per half slab as calclated
* in the water uptake code. EOMULT accounts for whether the slab
* is a full or half (plant in the middle - full or edge -half)
* multiply by 1/poprow * 100  to multiply be the amount of space between 
* plants (cm) the scales the slab to the area corresponding to one
* plant - hence the output is per plant
* 

       
           CurrentTrans=(AWUPSS-ActualTrans1)/EOMult/poprow*100
           p_ET=(potential_T-pot_ET)/EOMult/poprow*100
c pot_ET_EO is the EO summation from the previous time step.           
           p_ET_EO=potential_T_EO-pot_ET_EO          ! already per plant, summation of EO from hourly weather and wateruptake
                                                     ! pot_ET_EO is Potential_T_EO from the previous time step
                                                     ! Potential_T_EO is the summation of EO over the simulation
                                                     ! EO is potential transpiration rate
c potential_T is calculated as a summation of EOR*Step so it is per slab and is cumulative over the simulation
C EOR is calculated from ET demand in the crop model
c  pot_ET is the value of potential_T from the previous hour
           
C EO is the potential transpiration rate from Penman adjusted for plant cover, calculated in the weather subroutine
c  potential_T_EO is the summation of EO which is soil surface evap from EO and the Penman equation.
c EOR_p is a rough estimate of potential ET. p_ET is a better estimate using time step to 
c  integrate rather than use an average
           EOR_p=EOR*period/EOMult/poprow*100 ! EOR comes from ET_demand from plant
 

             
C 
C    Accumulate seasonal amounts       
C 
          SeasonRainfall =SeasonRainfall+CumRain
          SeasonActualInfil=SeasonActualInfil+ActualInfil
          SeasonPotenTranspiration=SeasonPotenTranspiration+
     &       p_ET
          SeasonActualTranspiration=SeasonActualTranspiration+
     &       CurrentTrans
          SeasonPotenSoilEvap=SeasonPotenSoilEvap+
     &       PotentialSoilEvaporation
          SeasonActualSoilEvap=SeasonActualSoilEvap+ActualSoilEvap
         Write(83,10) t,Date,PotentialSoilEvaporation, 
     &               ActualSoilEvap,p_ET_EO,
     &              p_ET, CurrentTrans,CumRain,ActualInfil,FLuxAct,
     &              Runoff,cover, PSIM, SeasonPotenSoilEvap, 
     &              SeasonActualSoilEvap, SeasonPotenTranspiration,
     &               SeasonActualTranspiration,  SeasonRainfall,  
     &              SeasonActualInfil
        
        Write(84,13) t,date,flux,ActualSoilEvap,Wattsm(itime),
     &     RNS,RNC,WATPOT,WATACT, WATRAT,RNLU,
     &     Difwat(itime),Cloud,Shadow(itime),Difint(itime),
     &      RowInc(itime),Parint(itime), Radint(itime)
     

	  flux=0
	  ActualSoilEvap=0
	  ActualInfil=0
	  ActualTrans1=AWUPSS
	  pot_ET=potential_T
	  pot_ET_EO=Potential_T_EO
	  timer=0
	  FluxAct=0.0
	  Runoff=0.0
	  CumRain=0.0
	  PotentialSoilEvaporation=0.0
	  Endif  !end loop for this timestep
C*** FROM HourlyOut.for END
           
      Return
C10    Call errmes(im,il)
C format 5 is water, 6 is roots
  5   Format (1x,A12, T20,A7,T37,A2,T47,A2, T56,A5,T68,A6,T82,A3,T90,
     &     A3,T102,A2, T114, A6, T126, A5, T138, A6)   ! format for nodal data header
  6   Format (1x,A16,T20,A5,T36,A2,T46,A2,T56,A5,T65,A7,T78,A7, T90,A6,
     &        T102, A6, T112, A11, T124, A8, T136, A8, T146, A5 )
 10   Format (1x,F17.6,',', A15,',',16(F16.3, ','), F16.3)  !format for flux data
 13   Format (1x,F17.6,',',A15,',',15(F16.3, ','), F16.3)  ! format for second flux data
 20   Format (1x,19A16)              !format for flux header
      End        
C/***************************************************/
CDT !0/24/2018 modified to output nodal values for roots
      Subroutine ArrRootOut(Date,x,y,t,NumNP,RMassM,RMassY,
     !            RDenM, RDenY,
     !           SinkArr,cSinkArr,gSinkArr,Area)
      Character *10 Date
      Dimension RMassM(NumNP), RMassY(NumNP), RDenM(NumNP), 
     !           RDenY(NumNP),Area(NumNP),x(NumNP),y(NumNP),
     !         SinkArr(NumNP),cSinkArr(NumNP),gSinkArr(NumNP)

      do n=1,NumNP
        Write(82,140)t, Date, x(n),y(n), n,RMassM(n), RMassY(n),
     !      RDenM(n), RDenY(n),
     !      SinkArr(n),cSinkArr(n),gSinkArr(n),Area(n)
        Enddo
C      Endif     
130   format('     Time           x             y     Temp   RTWT_',
     !        '          Sink_       cSink       gSink'/)
 140   format(1x,f14.6,',', A12,',', 2(f11.2,','),i8,',',
     !   7(G11.3,','),G11.3)
      Return
      End      
C/***************************************************/
      subroutine ArrNodOut(Date,x,y,NumNP,NumEl,t,NumElD,Head,      
     !          ThetaArr, vzArr,vxArr,QArr,ConcArr,TemprArr,gArr)
      Character*10 Date
      Dimension Head(NumNP),x(NumNP),y(NumNP),
     !         ThetaArr(NumNP), vzArr(NumNP),vxArr(NumNP),QArr(NumNP),
     !         ConcArr(NumNP),TemprArr(NumNP),gArr(NumNP)

        !MyDate=ToCalendarDate(mm,id,iyyy)
        Do n=1,NumNP
          Write(81,140)t, date, x(n),y(n),Head(n),ThetaArr(n),
     !     vzArr(n),vxArr(n),QArr(n),ConcArr(n),TemprArr(n),
     !      gArr(n)
        Enddo
 140   format(1x, f14.6,',', a12,',', 2(f11.2,','),7(G11.3,','),G11.3)
      Return
      End
  
  
