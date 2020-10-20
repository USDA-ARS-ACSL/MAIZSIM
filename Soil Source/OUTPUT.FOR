c pot_ET is used to save the accumulated potential transpiration in the previous N steps
c p_ET is used to save the difference between the accumulated transpiration in the 
c previous N step and N+1 step; thus in fact it is the potential transpiration rate per time step
C summary output variables
!**                    cumRain is the sum of VarBw(,1) which is the rain that is input
!**                   Actual infil is   Q(i) which is rain - et from varBW(i,3)
!**                   Flux_act comes from Qact(i) which is the actual flux into or out of a node
!**                   FLux_drain is drainage from the bottom of the profile
!**                   Runoff is calculated as Actual_infil-FLux_act. For a time step
!**                       it may include any evaporation
      subroutine Output ()
      Include 'public.ins'
      Include 'puweath.ins'
      Include 'puplant.ins'
      include 'PuSurface.ins'
      Integer Plevel,e
      
      Dimension temp(1000)
      Real PotentialSoilEvaporation, WidthT,Period
      Real SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil
      Character*10 date, Stime
      double precision Rain_Flux,Surface_Flux,Bottom_Flux,
     &     Runoff_Flux,Runoff_Terminal_Flux,
     &     Runoff_Left_Sum,Runoff_Right_Sum,WaterStorage,
     &     WaterStorage_pre
      double precision
     &     ActualInfil,     ! cccz: original infiltration, Q(n) based, so it is a cumulative value.
     &     ActualInfil_02,  ! cccz: this is based on surface mass balance, so totally global.
     &     ActualInfil_03   ! cccz: this is based on Qact(n), so should be local but can be generated to global, same to the fluxact with a "relu filter"
      double precision
     &     Runoff,       ! cccz: original runoff base on RO: for each surface node, treat rainfall and surface runoff from upper stream as (potential) infiltraiton and this runoff represent the downstream water quantity, REPEATED CALCULATION INCLUDED   ----------  Good if you think soil surface is a "river"
     &     Runoff_02,    ! cccz: runoff leave the surface: for the whole surface, treat rainfall as input and count the water leave the surface, then distributed evenly along the surface, DO NOT INCLUDE REPEAT CALCULATION, BUT for furrow, CAN BE ZERO even surface water flow inside the domain ----------  Good if you think runoff effects out of your soil domain
     &     Runoff_03     ! cccz: runoff contribution @ each node: for each surface node, treat rainfall (exclude water income from upper-stream), the real runoff contribution for each node, then average over the soil surface ------- Good for the local runoff, but eventually averaged
      Common /Outcom/NumPoint(2000),NodNum(2000,60),
     &     NumCell(1000),NumElPrint(1000,60),ModNum,
     &     Period, ETO, ActualTrans1,ActualSoilEvap,
     &     pot_ET,pot_ET_EO,Flux_drain, timer,
     &     FluxAct,PotentialSoilEvaporation,
     &     CumRain,
     &     SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil,CumulativeFlux,
     &     Runoff,Runoff_02,Runoff_03,
     &     ActualInfil,ActualInfil_02,ActualInfil_03
      Common /SurfaceWaterBalance/ Rain_Flux,Surface_Flux,Bottom_Flux,
     &     Runoff_Flux,Runoff_Terminal_Flux,
     &     Runoff_Left_Sum,Runoff_Right_Sum,
     &     WaterStorage,WaterStorage_pre

        
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
        Flux_drain=0.0
        fluxAct=0.0
        
        ActualSoilEvap=0
        PotentialSoilEvaporation=0.0
        CumRain=0.0
        
cccz Initialize the infiltration block
        ActualInfil=0
        ActualInfil_02=0
        ActualInfil_03=0
cccz
        
cccz Initialize the runoff block
        Runoff=0.0
        Runoff_02=0.0
        Runoff_03=0.0
cccz
        
        SeasonRainfall =0.0
        SeasonActualInfil=0.0
        SeasonPotenTranspiration=0.0
        SeasonActualTranspiration=0.0
        SeasonPotenSoilEvap=0.0
        SeasonActualSoilEvap=0.0

	  ActualTrans1=0.0
	  pot_ET=0.0
	  pot_ET_EO=0.0
        
cccz initialize these only for runoff output file
        Rain_Flux=0.0
        Surface_Flux=0.0
        Bottom_Flux=0.0
        Runoff_Flux=0.0
        Runoff_Terminal_Flux=0.0
        Runoff_Left_Sum=0.0
        Runoff_Right_Sum=0.0
        
        
	  Open(81,file=NodeGraphics)
	  open(82,file=ElemGraphics)  
	  Open(83,file=SurfaceGraphics, recl=400)
        Open(84,file=FluxGraphics)
       open(99,file=MassBalanceRunoffFileOut,status='unknown',recl=250)
         
      Write(81,5) "Date_time,", "Date,", "X,", "Y,", "hNew,","thNew,",
     &  "VZ,", "VX,","Q,", "ConcN,", "Temp,","GasConc"  
     
      Write(82,6) "Date_time,", "Date,", "X,", "Y,", "Node,", 
     &      "RMassM,","RMassY,", "RDenM,", "RDenY,",
     &      "WaterSink,", "NitSink,", "GasSink,", "Area"
         
      Write(83,20) "Date_time,", "Date,", "PSoilEvap,", 
     &      "ASoilEVap,", "PET_PEN,","PE_T_int,","transp,",
     &      "CumRain,","infil,","FLuxAct,", "Drainage,","Runoff,",
     &      "cover,","PSIM,", "SeasPSoEv,","SeasASoEv," ,"SeasPTran,",
     &      "SeasATran,","SeasRain,",    
     &      "SeasInfil,", "Runoff_02,", "Runoff_03"

       
       Write(84,20) "Date_time,", "Date,","Flux,","SoilEA,",
     &   "RNS,", "WATTSM,","RNC,","WATPOT,","WATACT,",
     &    "WATRAT,", "UpwardLW,",
     &   "DifWat,","Cloud,","Shadow,","DIfInt,","Rowinc,", "PARInt,",
     &   "RADInt"
       
cccz initialize for surface runoff
       write(99,'(8A20)') 'time','Date','Rainfall','Inifiltration',
     &    'Runoff@end', 'WaterStorage', 
     &    'LeftDischarge','RightDischarge'
       
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

cccz We change the usage of VarBW and VarBW_Record
cccz VarbW is the one adjusted with the runoff module, while VarBW_Record is the one from weather calculation
        do i=1,NumBp
	    n=KXB(i)
cccz             If (VarBW(i,3).GE.0) Then
cccz      	      PotentialSoilEvaporation=PotentialSoilEvaporation+
cccz     & VarBW(i,3)*width(i)*step
cccz             EndIf
cccz one critical issue we need to face is varbW changed in runoff, recommand varBW_record             
cccz	           CumRain=CumRain+varBW(i,1)*width(i)*step
          if(VarBW_Record(i,3).ge.0.0D0) then
            PotentialSoilEvaporation=PotentialSoilEvaporation
     &        +VarBW_Record(i,3)*width(i)*step
          endif
cccz          if (codeW(n).eq.-4) Then
            CumRain=CumRain+VarBW_Record(i,1)*width(i)*step
cccz          endif
	  enddo
	
C Actual values

CDT note that flow takes place in seepage faces when the boundary code is positive - this indicates a fixed
C  pressure head and flux is calculated CodeW >0 indicates head is set at 0 and flux is calculated.
         timer=timer+step
      	 do i=1,NumBp
	     n=KXB(i)
	     
	     
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))
     &             .or.(CodeW(n).eq.(1))) then
	        Flux_drain=Flux_drain+Q(n)*step
	        
           endif
c case for evaporation

cccz changed varbw to varbw_record
          if ((abs(CodeW(n)).eq.4).and.(VarBw_Record(i,3).GT.0)) then
            ActualSoilEvap=ActualSoilEvap-min(0.0,Q(n))*step
          endif
c case for rainfall
          if (abs(CodeW(n)).eq.4) then
            ActualInfil=ActualInfil+max(0.0,Q(n))*step
            ActualInfil_03=ActualInfil_03+max(0.0,Qact(n))*step
cccz calculation of runoff
            Runoff=Runoff+RO(n)*step
            
cccz the calculation of runoff_03 is relatively complicated
cccz      first do the input flow, the amound of water from upper stream should be Q-varbw_record (infiltration-rainfall)
cccz      so eventually, it will be "Rain-Qact"
cccz      which means that point has positive contribution to the runoff (>0, increase the runoff quantity) or absorb runoff (<0)
cccz      this is meaningful, recall teh initiation of runoff is @ some time @ some node, we have "Rain>Qact"
           Runoff_03=Runoff_03
     &         +(RO(n)-(Q(n)-VarBW_Record(i,1)*width(i)))*step

c           RO(n)=0.0   ! we will use this to store ponded water for next time period
             FluxAct=FluxAct+max(0.0,QAct(n))*step   
c            FluxAct=FluxAct+QAct(n)*step
          endif
	  enddo

C*** FROM HourlyOut.for END

cccz
cccz Surface runoff calculation
          Surface_Flux=FluxAct
          Bottom_Flux=Flux_drain
          Runoff_Flux=Runoff
          Rain_Flux=CumRain
c we consider the runoff from the terminal edge of the soil profile

cccz calculation of terminal flux for mass balance in runoff module
          Runoff_Terminal_Flux=Runoff_Terminal_Flux
     &     +(abs(RunoffLeft)+abs(RunoffRight))*step                   ! Runoff from Saint-Venant Eq.
     &     +(abs(RunoffLeft_Efflux)+abs(RunoffRight_Efflux))*step     ! Runoff above the CriticalH_R (Force the Mass Balance)      
          Runoff_Left_Sum=Runoff_Left_Sum
     &     +abs(RunoffLeft+RunoffLeft_Efflux)*step
          Runoff_Right_Sum=Runoff_Right_Sum
     &     +abs(RunoffRight+RunoffRight_Efflux)*step
cccz calculation is done, just do an assignment         
          Runoff_02=Runoff_Terminal_Flux



11    If(abs(time-tNext(ModNum)).lt.0.001*Step.OR.lInput.eq.1) then
          iday=int(t)
        call caldat(iday,mm,id,iyyy) 
        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          Call ArrNodOut (Date,x,y,NumNP,NumEl,t,NumElD,hNew,
     !                   ThNew,vz,vx,Q,Conc(1,1),Tmpr,g(1,1))
        iitime=int(time)
        diff=abs(time-int(time))
        If((NShoot.ne.0).and.(abs((time-int(time))).gt.(0.99))) 
     !       then
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
        
cccz for runoff water storage calculation
cccz Calculate the total amount of water stay on soil surface
cccz because WATERSTORAGE is not a cumulative value, so we just do it once before the print.
         WaterStorage_pre=WaterStorage
         WaterStorage=0.0D0
         do i=1,SurNodeIndex-1
           WaterStorage=WaterStorage+0.5D0
     &      *(slopeCoord(i+1,1)-slopeCoord(i,1))*(h_Stay(i+1)+h_Stay(i))
         enddo 
cccz
cccz calculation is done, just do an assignment         
          Runoff_02=Runoff_Terminal_Flux
          ActualInfil_02=CumRain-
     &      (WaterStorage-WaterStorage_pre)-Runoff_Terminal_Flux



* Anything calculated from VarBW should not be divided by 24 
* output as mm per plant, essentially dividing by popslab
* VarBW is cm/cm2 per day. Any fluxes from the finite element code 
*  are cm of water in the width of the nodes per day. Thus these fluxes are not
* multiplied by width

              PotentialSoilEvaporation=PotentialSoilEvaporation
     &        /EOMult/poprow*100.0D0
              ActualSoilEvap=ActualSoilEvap/EOMult/poprow*100.0D0
              !ActualInfil=ActualInfil/EOMult/poprow*100
              !FluxAct=FluxAct/EOMult/poprow*100
              !Runoff=Runoff/EOMult/poprow*100
              !CumRain=CumRain/EOMult/poprow*100
              !bottomFlux=bottomFlux/EOMult/poprow*100
              !(/RowSP/2) gets us per cm2 now for 10 minute period
              !*6 gives us one hour
cccz              ActualInfil=ActualInfil/(RowSP/2)*10*24
cccz              FluxAct=FluxAct/(RowSP/2)*10*24
cccz              Runoff=Runoff/(RowSP/2)*10*24
cccz              CumRain=CumRain/(RowSP/2)*10*24
cccz              bottomFlux=bottomFlux/(RowSP/2)*10*24    !cm/cm2/day

      
              FluxAct=FluxAct/EOMult/poprow*100.0D0
              
cccz the averaged process for actual infiltraiton
              ActualInfil=ActualInfil/EOMult/poprow*100.0D0
              ActualInfil_02=ActualInfil_02/EOMult/poprow*100.0D0
              ActualInfil_03=ActualInfil_03/EOMult/poprow*100.0D0
              
cccz the averaged process for runoff
              Runoff=Runoff/EOMult/poprow*100.0D0
              Runoff_02=Runoff_02/EOMult/poprow*100.0D0
              Runoff_03=Runoff_03/EOMult/poprow*100.0D0  
cccz
              CumRain=CumRain/EOMult/poprow*100.0D0
              Flux_drain=Flux_drain/EOMult/poprow*100.0D0   !cccz unit [water:width(cm)*flux(cm/day)*timestep(day)]/EOMult[0.5 or 1]/[population/distance(cm)]
              CumulativeFlux=CumulativeFlux+Flux_drain


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
          SeasonActualInfil=SeasonActualInfil+ActualInfil_03
          SeasonPotenTranspiration=SeasonPotenTranspiration+
     &       p_ET
          SeasonActualTranspiration=SeasonActualTranspiration+
     &       CurrentTrans
          SeasonPotenSoilEvap=SeasonPotenSoilEvap+
     &       PotentialSoilEvaporation
          SeasonActualSoilEvap=SeasonActualSoilEvap+ActualSoilEvap
          
          
         Write(83,10) t,Date,PotentialSoilEvaporation, 
     &               ActualSoilEvap,p_ET_EO,
     &              p_ET, CurrentTrans,CumRain,ActualInfil_03,FLuxAct,
     &              Flux_drain, Runoff,                               ! cccz: Dennis, think about the runoff you like to print here.
     &              cover, PSIM, SeasonPotenSoilEvap,
     &              SeasonActualSoilEvap, SeasonPotenTranspiration,
     &               SeasonActualTranspiration,  SeasonRainfall,  
     &              SeasonActualInfil, Runoff_02, Runoff_03
       
cccz        Write(84,13) t,date,bottomflux,ActualSoilEvap,
cccz     &     CumulativeWattsm(itime),
cccz     &     RNS,RNC,WATPOT,WATACT, WATRAT,RNLU,
cccz     &     Difwat(itime),Cloud,Shadow(itime),Difint(itime),
cccz     &      RowInc(itime),Parint(itime), Radint(itime)
        
        Write(84,13) t,date,Flux_drain,ActualSoilEvap, Wattsm(itime),
     &     RNS,RNC,WATPOT,WATACT, WATRAT,RNLU,
     &     Difwat(itime),Cloud,Shadow(itime),Difint(itime),
     &      RowInc(itime),Parint(itime), Radint(itime)
        
        write(99,'(1F12.4,A12,F14.4,F14.4,F14.4,F14.4,F14.4,F14.4)')       
     &     time, date, Rain_Flux, Surface_Flux, Runoff_Terminal_Flux, 
     &     WaterStorage, Runoff_Left_Sum,Runoff_Right_Sum
        

          

	  Flux_drain=0
	  ActualSoilEvap=0
	  ActualTrans1=AWUPSS
	  pot_ET=potential_T
	  pot_ET_EO=Potential_T_EO
	  timer=0
	  FluxAct=0.0
	  CumRain=0.0
	  PotentialSoilEvaporation=0.0
cccz reset ActualInfil values        
	  ActualInfil=0.0
        ActualInfil_02=0.0
        ActualInfil_03=0.0
cccz
cccz reset runoff values
	  Runoff=0.0
        Runoff_02=0.0
        Runoff_03=0.0
cccz
        
cccz runoff water balance term
         Surface_Flux=0.0
         Bottom_Flux=0.0
         Runoff_Flux=0.0
         Runoff_Terminal_Flux=0.0
         Runoff_Left_Sum=0.0
         Runoff_Right_Sum=0.0
         Rain_Flux=0.0
cccz end

        
	  Endif  !end loop for this timestep
C*** FROM HourlyOut.for END
           
      Return
C10    Call errmes(im,il)
C format 5 is water, 6 is roots
  5   Format (1x,A12, T20,A7,T37,A2,T47,A2, T56,A5,T68,A6,T82,A3,T90,
     &     A3,T102,A2, T114, A6, T126, A5, T138, A6)   ! format for nodal data header
  6   Format (1x,A16,T20,A5,T36,A2,T46,A2,T56,A5,T65,A7,T78,A7, T90,A6,
     &        T102, A6, T112, A11, T124, A8, T136, A8, T146, A5 )
 10   Format (1x,F17.6,',', A15,',',19(F16.3, ','), F16.3)  !format for flux data
 13   Format (1x,F17.6,',',A15,',',15(F16.3, ','), F16.3)  ! format for second flux data
 20   Format (1x,22A16)              !format for flux header
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
  
  
