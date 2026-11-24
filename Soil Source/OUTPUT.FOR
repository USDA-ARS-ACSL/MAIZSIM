c pot_ET is used to save the accumulated potential transpiration in the previous N steps
c p_ET is used to save the difference between the accumulated transpiration in the 
c previous N step and N+1 step; thus in fact it is the potential transpiration rate per time step
c all output into the g05 file is now per cm2 - water is in mm
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
      Include 'nitvar.ins'
      include 'PuSurface.ins'
      Integer Plevel,e, hourCount
      
      Dimension temp(1000)
      Real PotentialSoilEvaporation, WidthT,Period
      Real SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil,
     &     SumArrayTh, MaxY, Cumulative_OM_N_sink
      Character*10 date, Stime
cccz output for runoff module
      double precision SurInput_Flux,Surface_Flux,Bottom_Flux,
     &     Runoff_Flux,Runoff_Terminal_Flux,Flux_drain,
     &     Runoff_Left_Sum,Runoff_Right_Sum,WaterStorage,
     &     WaterStorage_pre ,Surface_Gas_Flux,FluxGas,
     $     nitrogenFlux
      
      double precision
     &      ActualInfil,     ! cccz: original infiltration, Q(n) based, so it is a cumulative value.
     &      ActualInfil_02,  ! cccz: this is based on surface mass balance, so totally global.
     &      ActualInfil_03   ! cccz: this is based on Qact(n), so should be local but can be generated to global, same to the fluxact with a "relu filter"
      double precision
     &     Runoff,       ! cccz: original runoff base on RO: for each surface node, treat rainfall and surface runoff from upper stream as (potential) infiltraiton and this runoff represent the downstream water quantity, REPEATED CALCULATION INCLUDED   ----------  Good if you think soil surface is a "river"
     &     Runoff_02,    ! cccz: runoff leave the surface: for the whole surface, treat rainfall as input and count the water leave the surface, then distributed evenly along the surface, DO NOT INCLUDE REPEAT CALCULATION, BUT for furrow, CAN BE ZERO even surface water flow inside the domain ----------  Good if you think runoff effects out of your soil domain
     &     Runoff_03     ! cccz: runoff contribution @ each node: for each surface node, treat rainfall (exclude water income from upper-stream), the real runoff contribution for each node, then average over the soil surface ------- Good for the local runoff, but eventually averaged

      Common /Outcom/NumPoint(2000),NodNum(2000,60),
     &     NumCell(1000),NumElPrint(1000,60),ModNum,
     &     Period, ETO, ActualTrans1,ActualSoilEvap,
     &     pot_ET,pot_ET_EO,bottomFlux, timer,
     &     FluxAct,PotentialSoilEvaporation,
     &     CumRain,Cumulative_OM_N_sink,
     &     SeasonRainfall, SeasonPotenTranspiration,
     &     SeasonActualTranspiration, SeasonPotenSoilEvap,
     &     SeasonActualSoilEvap,SeasonActualInfil,CumulativeFlux,
     &     Runoff,Runoff_02,Runoff_03,
     &     ActualInfil,ActualInfil_02,ActualInfil_03,
     &     FluxGas(NumGD),Surface_Gas_Flux(NumGD),
     &     SumArrayTh(NumNPD), hourCount,
     &     nitrogenFlux,Flux_drain, 
     &      ThetaFull, MaxY
      Common /SurfaceWaterBalance/ 
     &     SurInput_Flux,Surface_Flux,Bottom_Flux,
     &     Runoff_Flux,Runoff_Terminal_Flux,
     &     Runoff_Left_Sum,Runoff_Right_Sum,
     &     WaterStorage,WaterStorage_pre

        
      t=Time
      
      If(lInput.eq.1) then
C
C  Initialize 
C
        If (DailyOutput.eq.1) Period=1.0
        If (HourlyOutput.eq.1) Period=1.0/24.0D0
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum) = time


C*** FROM HourlyOut.for BEGIN   
        SumArrayTh=ThNew
        hourCount=0
        timer=0
        Flux_drain=0.0
        fluxAct=0.0

cccz Initialize the runoff block
        Runoff=0.0
        Runoff_02=0.0
        Runoff_03=0.0
cccz
        
        ActualSoilEvap=0
        PotentialSoilEvaporation=0.0
        CumRain=0.0
        
cccz Initialize the infiltration block
        ActualInfil=0
        ActualInfil_02=0
        ActualInfil_03=0
cccz
        
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
        nitrogenFlux=0.0
        Cumulative_OM_N_sink=0.0
        
cccz initialize these only for runoff output file
        SurInput_Flux=0.0
        Surface_Flux=0.0
        Bottom_Flux=0.0
        Runoff_Flux=0.0
        Runoff_Terminal_Flux=0.0
        Runoff_Left_Sum=0.0
        Runoff_Right_Sum=0.0

C    For gas flux calculation
        FLuxGas=0.0
        Surface_Gas_Flux=0.0
        
C declare files for output
        
	  Open(81,file=NodeGraphics)                      !G03
	  open(82,file=ElemGraphics)                      !G04
	  Open(83,file=SurfaceGraphics, recl=450)         !G05
        Open(84,file=FluxGraphics)                      !G06
        Open(85,file=OrganicMatterGraphics,recl=300)    !G07
       open(99,file=MassBalanceRunoffFileOut, recl=250)
 
c NodeGrpahics       
      Write(81,5) "Date_time,", "Date,", "X,", "Y,", "hNew,","thNew,",
     &  "Vx,", "Vy,", "Q,", "NO3N,", "NH4N,","Temp,","Area,", 
     &  "CO2Conc," , "O2Conc,", "N2OConc"
c Element Graphics 
      Write(82,6) "Date_time,", "Date,", "X,", "Y,", "Node,", 
     &      "RMassM,","RMassY,", "RDenM,", "RDenY,",
     &      "WaterSink,", "NitSink,", "GasSink"
C G05 file  Surface Graphics       
      Write(83,20) "Date_time,", "Date,", "PSoilEvap,", 
     &      "ASoilEVap,", "PET_PEN,","PE_T_int,","transp,",
     &      "CumRain,","infil,","FLuxAct,", "Drainage,",
     &       "N_Leach,","Runoff,",
     &      "cover,","PSIM,", "SeasPSoEv,","SeasASoEv," ,"SeasPTran,",
     &      "SeasATran,","SeasRain,",    
     &      "SeasInfil,","ThetaAvail,","CO2FLux,","O2FLux,","N2OFLux"

c Flux Graphics       
       Write(84,20) "Date_time,", "Date,","Flux,","SoilEA,",
     &   "WATTSM,", "RNS,", "RNC,","WATPOT,","WATACT,",
     &    "WATRAT,", "UpwardLW,",
     &   "DifWat,","Cloud,","Shadow,","DIfInt,","Rowinc,", "PARInt,",
     &   "RADInt"

c Organic Matter Graphics

       Write(85,21) "Date_time,", "Date,", "X,", "Y,", "Humus_N,",
     &               "Humus_C,", "Litter_N,", "Litter_C,",
     &            "Manure_N,", "Manure_C,", "Root_N,", "Root_C"
c Runoff Mass Balance     
cccz initialize for surface runoff
       write(99,'(8A20)') 'time,','Date,','Rainfall,','Inifiltration,',
     &    'Runoff@end,', 'WaterStorage,', 
     &    'LeftDischarge,','RightDischarge'
       
cccz -----------------------------------------------------------------------
cccz initialize the mulch output file (analogy to the soil output g03)
       Open(2102,file=MassBalanceMulchFileOut)
       Write(2102,2103) "Date_time,","Date,","X,","Y,","hNew,",
     &  "thNew,","Temp"
       iday=int(time)
       call caldat(iday,mm,id,iyyy) 
       write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
       do n=1, numMulchNode
         Write(2102,2104) Time, Date, mulchNodeCoord(n,1),
     &     mulchNodeCoord(n,2),MulchNodehNew(n),MulchNodeThNew(n),
     &     MulchNodeTmpr(n)
       enddo
       


cccz -----------------------------------------------------------------------
       
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
        cRoot=0.0
        nRoot=0.0
cccz We change the usage of VarBW and Varbw_Air
cccz VarbW is the one adjusted with the runoff module, while Varbw_Air is the one from weather calculation
        do i=1,NumBp
	    n=KXB(i)
cccz one critical issue we need to face is varbW changed in runoff, recommand Varbw_Air             
cccz	           CumRain=CumRain+varBW(i,1)*width(i)*step
          if(Varbw_Air(i,3).ge.0.0D0) then
            PotentialSoilEvaporation=PotentialSoilEvaporation
     &        +Varbw_Air(i,3)*width(i)*step
          endif
cccz          if (codeW(n).eq.-4) Then
            CumRain=CumRain+Varbw_Air(i,1)*width(i)*step
            SurInput_Flux=SurInput_Flux+Varbw_Mulch(i,1)*width(i)*step
            Cumulative_OM_N_sink=Cumulative_OM_N_sink+cSink_OM(i,1)
     &         *width(i)* step*thNew(i)*NO3mass_2_Nmass
cccz         endif
	  enddo
	
C Actual values

CDT note that flow takes place in seepage faces when the boundary code is positive - this indicates a fixed
C  pressure head and flux is calculated CodeW >0 indicates head is set at 0 and flux is calculated.
         timer=timer+step
      	 do i=1,NumBp
	     n=KXB(i)
	     
C calculate N flux in drainage, only consider flow out of the domain (<0) units are g N/cm2 (convert from ug)	   
C make sure that the drainage is only calculated for the bottom boundary
        if (n.gt.nSurf) then
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))
     &             .or.(CodeW(n).eq.(1))) then
	        Flux_drain=Flux_drain+Q(n)*step
              nitrogenflux=nitrogenflux+min(0.0,Q(n))
     &             *conc(n,1)*step*NO3mass_2_Nmass/1.0e6              ! g N/cm2
      
           endif
        endif
c case for evaporation

cccz changed varbw to Varbw_Air
          if ((abs(CodeW(n)).eq.4).and.(Varbw_Air(i,3).GT.0)) then
            ActualSoilEvap=ActualSoilEvap-min(0.0,Q(n))*step
          endif
c case for rainfall
          if (abs(CodeW(n)).eq.4.or.(((pondingByFlux.eq.1).or.
     &         (pondingByHead.eq.1)))) then
            ActualInfil=ActualInfil+max(0.0,Qact(n))*step
            ActualInfil_03=ActualInfil_03+max(0.0,Qact(n))*step
cccz calculation of runoff
            Runoff=Runoff+max(Q(i)-Qact(i),0.0D0)*step
            
cccz the calculation of runoff_03 is relatively complicated
cccz      first do the input flow, the amound of water from upper stream should be Q-varbw_record (infiltration-rainfall)
cccz      so eventually, it will be "Rain-Qact"
cccz      which means that point has positive contribution to the runoff (>0, increase the runoff quantity) or absorb runoff (<0)
cccz      this is meaningful, recall teh initiation of runoff is @ some time @ some node, we have "Rain>Qact"
           Runoff_03=Runoff_03
     &         +(RO(n)-(Q(n)-Varbw_Mulch(i,1)*width(i)))*step    
            
c           RO(n)=0.0   ! we will use this to store ponded water for next time period
            FluxAct=FluxAct+max(0.0,QAct(n))*step         !QAct=cm3 at a node /time  , FluxAct =total cm3 of water in all the boundary nodes
c            FluxAct=FluxAct+QAct(n)*step
      endif
c Calculate the surface flux from gas mover   
      Do jjj=1,NumG
          if (CodeG(n).eq.-4. or. CodeG(n).eq.1) then
             FluxGas(jjj)=FluxGas(jjj)+QGas(n,jjj)*step                   !ugCo2 summed over all boundary nodes   QGas[ugco2/day]*[day], FLuxGas[ugco2 over the surface], this has both positive (soil to atm) and negative fluxes (in to the soil)
          end if
      end do
	  enddo

C*** FROM HourlyOut.for END

C Gas surface flux
          Surface_Gas_Flux=FluxGas

cccz
cccz Surface runoff calculation
          Surface_Flux=FluxAct
          Bottom_Flux=Flux_drain
          Runoff_Flux=Runoff
          SurInput_Flux=SurInput_Flux
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

          
          if(DailyOutput.eq.1) then
             SumArrayTh=SumArrayTh + ThNew
             hourCount=hourCount+1
           endif
11    If(abs(time-tNext(ModNum)).lt.0.001*Step.OR.lInput.eq.1) then
          iday=int(t)
          if (DailyOutput.eq.1) then
             SumArrayTh=SumArrayTh/Real(hourCount)
           endif
        call caldat(iday,mm,id,iyyy) 
        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          Call ArrNodOut (Date,x,y,NumNP,t,hNew,
     !      ThNew,vz,vx,Q,Conc(1,1),NH4,Tmpr,g(1,1), g(1,2),
     !      g(1,3),ugGasCm3air_to_ppm,NumG, NodeArea)
        iitime=int(time)
        diff=abs(time-int(time))
        If((NShoot.ne.0).and.(abs((time-int(time))).gt.(0.99))) 
     !       then
             Call ArrRootOut (Date,x,y,t, NumNP,RMassM,RMassY,
     !            RDenM, RDenY, Sink,
     !            cSink(1,1),gSink(1,1))
          ! add calcs of root mass in cell as ppm *(g/g)
            
               
        Endif 
          ! calculate N and C in roots (ug/cm3 soil) for output
          ! note that roots are g/cm3 of soil
          ! no need to use BD
         DO n=1, NumNP
            m=matnumN(n)
            cRoot(n)=RMassM(n)*cContentRootM
            cRoot(n)=cRoot(n)+RMassY(n)*cContentRootY
            nRoot(n)=RMassM(n)*nContentRootM
            nRoot(n)=nRoot(n)+RMassY(n)*nContentRootY
         end do
            cRoot=cRoot*1e6
            nRoot=nRoot*1e6
            
          Call ArrOrganicMatterOut(Date,x,y,NumNP,t,Ch,Nh,
     !       Cl,Nl, Cm, Nm, cRoot, nRoot)
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
     &      *(slopeCoord(i+1,1)-slopeCoord(i,1))*(h_Pond(i+1)+h_Pond(i))
         enddo 
cccz
cccz calculation is done, just do an assignment         
          Runoff_02=Runoff_Terminal_Flux
          ActualInfil_02=SurInput_Flux-
     &      (WaterStorage-WaterStorage_pre)-Runoff_Terminal_Flux


* Anything calculated from VarBW should not be divided by 24 
* output as mm per plant, essentially dividing by popslab
* VarBW is cm/cm2 per day. Any fluxes from the finite element code 
*  are cm of water in the width of the nodes per day. Thus these fluxes are not
* multiplied by width

cccz unit [water:width(cm)*flux(cm/day)*timestep(day)]/gridWidth * 10 mm/cm
              PotentialSoilEvaporation=PotentialSoilEvaporation
     &        /gridWidth*10.0
              ActualSoilEvap=ActualSoilEvap/gridWidth*10.0
              FluxAct=FluxAct/gridWidth*10.0
              
cccz the averaged process for actual infiltraiton
              ActualInfil=ActualInfil/gridWidth*10.0
              ActualInfil_02=ActualInfil_02/gridWidth*10.0
              ActualInfil_03=ActualInfil_03/gridWidth*10.0
              
cccz the averaged process for runoff
              Runoff = Runoff/gridWidth*10.0
              Runoff_02=Runoff_02/gridWidth*10.0
              Runoff_03=Runoff_03/gridWidth*10.0  
cccz
              CumRain=CumRain/gridWidth*10.0
              Flux_drain=Flux_drain/gridWidth*10.0   !cccz unit [water:width(cm)*flux(cm/day)*timestep(day)]/gridWidth *10 mm/cm
              CumulativeFlux=CumulativeFlux+Flux_drain

          
* dt - value of awupss and pot_et are per half slab as calclated
* in the water uptake code. 
* output is mm/cm2
* 

           Cumulative_OM_N_sink=Cumulative_OM_N_sink/gridWidth*10.0
           CurrentTrans=(AWUPSS-ActualTrans1)/gridWidth*10.0
           p_ET=(potential_T-pot_ET)/gridWidth*10.0
c pot_ET_EO is the EO summation from the previous time step.           
           p_ET_EO=potential_T_EO-pot_ET_EO          ! already per plant, summation of EO from hourly weather and wateruptake
                                                     ! pot_ET_EO is Potential_T_EO from the previous time step
                                                     ! Potential_T_EO is the summation of EO over the simulation
                                                     ! EO is potential transpiration rate
c potential_T is calculated as a summation of EOR*Step so it is per slab and is cumulative over the simulation
C EOR is calculated from ET demand in the crop model
c  pot_ET is the value of potential_T from the previous hour
           
C EO is the potential transpiration rate from Penman adjusted for plant cover, calculated in the weather subroutine
c  potential_T_EO is the summation of EO which is soil surface evap from EO and the Penman equation.vi
c EOR_p is a rough estimate of potential ET. p_ET is a better estimate using time step to 
c  integrate rather than use an average
c EOR is mm/cm2
           EOR_p=EOR*period/gridWidth*10.0 ! EOR comes from ET_demand from plant
 

             
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
C        Nitrogen flux conversion from g per slab to kg/ha   
         nitrogenFlux= nitrogenFlux/gridWidth*10**5        !10^5= 10^4 m2/ha x 10^4 cm2/m2 / 10^3 g/kg
C Calculate relative available water
              ThetaRelative=RelativeWaterContent(X,Y, ThAvail,ThFull,
     &            rtwt, NumNP, KX, NumElD,NumEl, NMat, MatNumN)
     
C     Gas unit conversion to ugCo2/cm2
          Surface_Gas_Flux=Surface_Gas_Flux/gridWidth    ! Unit ugCo2 per plant area
                
         Write(83,10) t,Date,PotentialSoilEvaporation, 
     &               ActualSoilEvap,p_ET_EO,
     &              p_ET, CurrentTrans,CumRain,ActualInfil,FLuxAct,
     &              Flux_drain, nitrogenFlux, Runoff_02,cover, PSIM, 
     &              SeasonPotenSoilEvap,
     &              SeasonActualSoilEvap, SeasonPotenTranspiration,
     &               SeasonActualTranspiration,  SeasonRainfall,  
     &              SeasonActualInfil, ThetaRelative,
     &              (Surface_Gas_Flux(jjj),jjj=1,NumG)
        
        
        Write(84,13) t,date,bottomflux,ActualSoilEvap, Wattsm(itime),
     &     RNS,RNC,WATPOT,WATACT, WATRAT,RNLU,
     &     Difwat(itime),Cloud,Shadow(itime),Difint(itime),
     &      RowInc(itime),Parint(itime), Radint(itime)
     &                  ! Unit ugCo2 per plant area              
        
        write(99,'(1F12.4,A12,F14.4,F14.4,F14.4,F14.4,F14.4,F14.4)')    
     &    time, date, SurInput_Flux, Surface_Flux, Runoff_Terminal_Flux,
     &    WaterStorage, Runoff_Left_Sum,Runoff_Right_Sum
        
        
cccz ---------------------- mulch ------------------------------------------
cccz write mulch output, when the mulch module is running "BoolMulchApply=1" (analogy to the soil output g03)
       if(BoolMulchApply.eq.1) then
        iday=int(time)
        call caldat(iday,mm,id,iyyy) 
        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
cccz output for node-base values
        do n=1,numMulchNode
         Write(2102,2104) Time,Date,mulchNodeCoord(n,1),
     &     mulchNodeCoord(n,2),MulchNodehNew(n),MulchNodeThNew(n),
     &     MulchNodeTmpr(n)
        enddo
      endif        
cccz ---------------------- finish mulch printing --------------------------
!        write(*,*) Cumulative_OM_N_sink   
          
        if (DailyOutput.eq.1) then
             SumArrayTh=0.0
             hourCount=0
           endif
	  Flux_drain=0
	  ActualSoilEvap=0
	  ActualInfil=0
	  ActualTrans1=AWUPSS
	  pot_ET=potential_T
	  pot_ET_EO=Potential_T_EO
	  timer=0
	  FluxAct=0.0
	  CumRain=0.0
        nitrogenFLux=0.0
	  PotentialSoilEvaporation=0.0

cccz reset runoff values
	  Runoff=0.0
        Runoff_02=0.0
        Runoff_03=0.0
cccz
cccz reset ActualInfil values        
	  ActualInfil=0.0
        ActualInfil_02=0.0
        ActualInfil_03=0.0
cccz
        
cccz runoff water balance term
         Surface_Flux=0.0
         Bottom_Flux=0.0
         Runoff_Flux=0.0
         Runoff_Terminal_Flux=0.0
         Runoff_Left_Sum=0.0
         Runoff_Right_Sum=0.0
         SurInput_Flux=0.0
cccz end

C Gas flux 
        FluxGas=0.0
        Surface_Gas_Flux=0.0 
        
	  Endif  !end loop for this timestep
C*** FROM HourlyOut.for END
           
      Return
C10    Call errmes(im,il)
C format 5 is water, 6 is roots Vy is the last on on this lin
  5   Format (1x,A12, T20,A7,T37,A2,T47,A2, T56,A5,T68,A6,T82,A3,T94,
     &     A3,T105,A2, T117, A6, T129, A6, T141, A5, 
     &     T152, A5,T165,A8,T176,A8,T187,A9 )   ! format for nodal data header
  6   Format (1x,A16,T20,A5,T36,A2,T46,A2,T56,A5,T65,A7,T78,A7, T90,A6,
     &        T102, A6, T112, A11, T124, A8, T136, A8)
 10   Format (1x,F17.6,',', A15,',',9(F16.3, ','), F14.5, ',',                  !relative water content  
     &        12(F16.3, ','),F16.3)  !format for flux data  [G03]
 13   Format (1x,F17.6,',',A15,',',15(F16.3, ','), F16.3)  ! format for second flux data [G06]
 20   Format (1x,25A16)              !format for flux header [G05]
 21   Format (1x,12A16)           ! format for Organic Matter header
2103  Format (1x,A12,T20,A12,T49,A2,T64,A2,T75,A10,T95,A10,T109,A10)
2104  Format (1x,F17.6,',',A15,',',4(F16.4,','),F16.4)     
      
      End        
C/***************************************************/
CDT !0/24/2018 modified to output nodal values for roots
      Subroutine ArrRootOut(Date,x,y,t,NumNP,RMassM,RMassY,
     !            RDenM, RDenY,
     !           SinkArr,cSinkArr,gSinkArr)
      Character *10 Date
      Dimension RMassM(NumNP), RMassY(NumNP), RDenM(NumNP), 
     !           RDenY(NumNP),x(NumNP),y(NumNP),
     !         SinkArr(NumNP),cSinkArr(NumNP),gSinkArr(NumNP)

      do n=1,NumNP
        Write(82,140)t, Date, x(n),y(n), n,RMassM(n), RMassY(n),
     !      RDenM(n), RDenY(n),
     !      SinkArr(n),cSinkArr(n),gSinkArr(n)
        Enddo
C      Endif     
130   format('     Time           x             y     Temp   RTWT_',
     !        '          Sink_       cSink       gSink'/)
 140   format(1x,f14.6,',', A12,',', 2(f11.2,','),i8,',',
     !   6(G11.3,','),G11.3)
      Return
      End      
C/***************************************************/
      subroutine ArrNodOut(Date,x,y,NumNP,t, Head,      
     !          ThetaArr, vzArr,vxArr,QArr,ConcArr,NH4Arr,TemprArr,
     !          gArr1,gArr2,gArr3,ugGasCm3air_to_ppm, NumG,areaArr)
      Character*10 Date
      Real x,y,t,Head,ThetaArr,vzArr,vxArr,QArr,ConcArr,NH4Arr,TemprArr,
     !      gArr1,gArr2,gArr3,ugGasCm3air_to_ppm, areaArr
      Integer NumNP,NumG, n

      Dimension Head(NumNP),x(NumNP),y(NumNP),
     !         ThetaArr(NumNP), vzArr(NumNP),vxArr(NumNP),QArr(NumNP),
     !         ConcArr(NumNP),TemprArr(NumNP),NH4Arr(NumNP),
     !         gArr1(NumNP), gArr2(NumNP), gArr3(NumNP),
     !         ugGasCm3air_to_ppm(NumG), areaArr(NumNP)
        !MyDate=ToCalendarDate(mm,id,iyyy)
        Do n=1,NumNP
          Write(81,140)t, date, x(n),y(n),Head(n),ThetaArr(n),
     !     vzArr(n),vxArr(n),QArr(n),ConcArr(n),NH4Arr(n), 
     !     TemprArr(n), areaArr(n),
     !     gArr1(n)*ugGasCm3air_to_ppm(1),
     !     gArr2(n)*ugGasCm3air_to_ppm(2), 
     !     gArr3(n)*ugGasCm3air_to_ppm(3)
          !convert the initial concentration of CO2 and O2 in [ug/cm3 air] to ppm
        Enddo
 140   format(1x, f14.6,',', a12,',', 2(f11.2,','),11(G11.3,','),G11.3)
      Return
      End
  
      subroutine ArrOrganicMatterOut(Date,x,y,NumNP,t,
     !   HumusCArr, HumusNArr, LitterCArr, LitterNArr,
     !    ManureCArr, ManureNArr, RootCArr, RootNArr)
      character*10 Date
      Integer NumNP
      Real   t
      Real HumusCArr, HumusNArr, LitterCArr, LitterNArr,
     !      ManureCArr, ManureNArr,RootCArr, RootNArr, 
     !      x, y 
      Dimension HumusCArr(NumNP), HumusNArr(NumNP), LitterCArr(NumNP), 
     !       LitterNArr(NumNP), ManureCArr(NumNP), ManureNArr(NumNP),
     !        RootCArr(NumNp), RootNArr(NumNP), 
     !       x(numNP),y(numNP)
c note output units are g/cm3 soil volume      
      
        Do n=1, NumNP
              Write(85,140) t, date, x(n),y(n), HumusNArr(n), 
     !         HumusCArr(n), 
     !        LitterNArr(n), LitterCArr(n), ManureNArr(n), 
     !        ManureCArr(n), RootNArr(n), RootCArr(n)
         EndDo 
  140   format (1x, F14.6,',', a12,',', 8(F11.2,','),F11.2, ',',F11.2)
       Return
       End
      
       
C calculate relative water content for profile to 100 cm    
C  returns the total cm of water in the domain - needs to be modified
C  by grid width to get cm of water per cm of soil
C this is  done after it is called

       
      function RelativeWaterContent(X,Y, ThetaAvail,
     !              ThetaFull,rtwt, NumNP, KX,  NumElD,
     !             NumEl, NMat, MatNumN)
      real X, Y, ThetaAvail, ThetaRel,ThetaFull,
     !    cii, bii, AE,ThetaAvailAll, ThetaFull1All, 
     !     rtwt
      integer NumNP, NMat, KX, NUS, K, i, j, l, n,
     !        NumElD,MatNumN
      
      
      Dimension Y(NumNP), X(NumNP), KX(NumElD,4), 
     !         ThetaAvail(NumNP), ThetaFull(NMat),
     !         rtwt(NumNP),MatNumN(NumNP),
     !         cii(3), bii(3)
            RelativeWaterContent=0.0
            ThetaFull1All=0.0
            ThetaAvailAll=0.0
      Do n=1,NumEl
		   NUS=4
		   if(KX(n,3).eq.KX(n,4)) NUS=3
c*         Loop on subelements
		   do k=1,NUS-2
			 i=KX(n,1)
			 j=KX(n,k+1)
			 l=KX(n,k+2)
			 Cii(1)=x(l)-x(j)
			 Cii(2)=x(i)-x(l)
			 Cii(3)=x(j)-x(i)
			 Bii(1)=y(j)-y(l)
			 Bii(2)=y(l)-y(i)
			 Bii(3)=y(i)-y(j)
			 AE=(Cii(3)*Bii(2)-Cii(2)*Bii(3))/2.
			
			 Thi=ThetaFull(MatNumN(i))
			 Thl=ThetaFull(MatNumN(l))
			 Thj=ThetaFull(MatNumN(j))
			 if (rtwt(i).le.1.0e-6) Thi=0.0
			 if (rtwt(j).le.1.0e-6) Thj=0.0
			 if (rtwt(l).le.1.0e-6) Thl=0.0
			 ThetaFullAll=ThetaFullAll+AE*(Thi+Thj+Thl)/3.
               
               Thi=ThetaAvail(MatNumN(i))
			 Thl=ThetaAvail(MatNumN(l))
			 Thj=ThetaAvail(MatNumN(j))
			 if (rtwt(i).le.1.0e-6) Thi=0.0
			 if (rtwt(j).le.1.0e-6) Thj=0.0
			 if (rtwt(l).le.1.0e-6) Thl=0.0
               ThetaAvailAll=ThetaAvailAll+AE*(Thi+Thj+Thl)/3.
     	   Enddo
       
      Enddo 
      RelativeWaterContent=
     !      min(1.0,ThetaAvailAll/(ThetaFullAll+.0001))    
        Return
        End