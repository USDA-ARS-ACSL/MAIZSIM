* This module read irrigation input details

      Subroutine Irrigation()
      Include 'public.ins'
      Include 'nitvar.ins'
      Include 'PuSurface.ins'
      Integer Max_Irrig_times,Irrig_times, IrrigHour
      Integer IrrigationApplied,HoursIrrigated
      Parameter (Max_Irrig_times=365)
      Character*10 Date(Max_Irrig_times)
      Character InString*132
      Real Amount_Irrig_Applied,AvgIrrRate,IrrigRate
      Parameter (PERIOD =1./24.)
      Common /Irrig/ tApplIrrig(Max_Irrig_times),
     !         StopIrrig(Max_Irrig_times),
     !         Amount_Irrig_Applied(Max_Irrig_times),
     !         Irrig_times,IrrigationApplied,ModNum,
     !         AvgIrrRate, IrrigRate(24),IrrigHour
      
      t=time
      IrrgStartHour=6.0                   ! irrigation start hour 6 AM   !if we want we can make this as an input by the user. 
      HoursIrrigated=0.0                  ! total number of hours of irrigation in a day
     
c
      If(lInput.eq.1) then
          Irrig_times=0                   ! Number of irrigations
          IrrigationApplied=0
          IrrigHour=1                     ! Array count for the irrigation rate
          ! Reading the irrigation details
          Open(40,file=IrrigationFile,err=20)
          Read(40,*,Err=20)
15        Read (40,'(A132)') InString
          if (InString(1:14).ne.'[Sprinkler]') goto 15
          Read(40,*,Err=20)
          Read(40,*,Err=20)
          Read(40,*,Err=20) AvgIrrRate    ! Average irrigation rate [cm/hour]
          Read(40,*,Err=20)
          Read(40,*,Err=20) Irrig_times   ! Number of days of irrigation
          Read(40,*,Err=20)
          Do i=1,Irrig_times
              Read(40,*,Err=20) Date(i), Amount_Irrig_Applied(i)  !irrigation applied = total irrigation per day in mm/day
          EndDo
          NumMod=NumMod+1
          ModNum=NumMod

          Do i=1,Irrig_times
              tApplIrrig(i)=julday(date(i))                   ! All the dates irrigation is applied stored in tApplIrrig
          EndDo
          
          Do i=1,24
              IrrigRate(i)=0.                                 ! To calculate the hourly irrigation rate every day
          Enddo

          if (Irrig_times.EQ.0) then
                tNext(ModNum)=1e22
           else 
                IrrigationApplied=IrrigationApplied+1 
           endif
                       ! this needs to be set to 1 to process the first irrigation event
          Close(40)
      Endif  ! end reading irrigation details


      If (Irrig_times.GT.0) then   !skip if no  irrigation (i.e., IrrigationApplied == 0)
        If (Abs(time-tApplIrrig(IrrigationApplied)).lt.0.001*Step) then !if day of irrigiation
         jj=IrrigationApplied
! dtIrrig is time in days to it takes for this amount to rain using the average hourly rainfall rate (IR)
! check to see if this amount will be reached in the current hour
! if it is too much for the hour, then the excess is carried over to the next hour
         
         
         DtIrrig=Amount_Irrig_Applied(jj)/10./AvgIrrRate  !gives number of hours in a day with irrigation [/10.0 to convert mm to cm]
         i=1
         precision_threshold=0.0001
        Do While (DtIrrig > precision_threshold)
           If((DtIrrig).gt.PERIOD*24.0) then            !period is 1/24
            IrrigRate(i)=AvgIrrRate             !cm/hour
            DtIrrig=DtIrrig-PERIOD*24             !how many more hours of irrigation is pending
            HoursIrrigated=HoursIrrigated+1
           Else                        !if dtIrrig=Period or < period
                                      ! WG
            IrrigRate(i)=DtIrrig*AvgIrrRate/(PERIOD*24.0)     !hour/daycm/hour /(1hour/24)= cm/hour
            DtIrrig=0.0                !no more irrigation
            HoursIrrigated=HoursIrrigated+1
            tNext(ModNum)=tApplIrrig(jj)+IrrgStartHour*period
            StopIrrig(jj)=tNext(ModNum)+HoursIrrigated*period
            IrrigHour=1
           Endif !! if DtIrrig>precision_threshold 
         i=i+1
         End Do
        Endif  !! if day of irrigiation
       End if  !! end of !skip if no  irrigation
   
      

655    If(Abs(time-tNext(ModNum)).lt.0.001*Step) then          !If its the day and time to irrigate
          jj=IrrigationApplied
          if (time.LT.StopIrrig(jj)) then                     !from start to the stop time
          
            !Here we pass the irrigation soil domain   (same format as in AutoIrrigate routine)
            do k=1, NumBp
              i=KXB(k)
              If(abs(CodeW(i)).eq.4) then             !checking if the surface boundary is either head or flux based
                  if (Q(i).lt.0) then                 !Flow out of the soil
                  CodeW(i)=4                          !Constant head
cccz since we treat irrigation as rainfall, we directly change the varbw terms
cccz so varbw_air will be rainfall+irrigation
              Varbw_Air(k,1)=Varbw_Air(k,1)+IrrigRate(IrrigHour)/period
          Varbw_Mulch(k,1)=Varbw_Mulch(k,1)+IrrigRate(IrrigHour)/period
              VarBW(k,1)=Varbw(k,1)+IrrigRate(IrrigHour)/period
cccz then re-calculation of Q, folloiwng VarBW and the starndard method show in weather
              VarBW_Air(k,3)=VarBW_Air(k,2)-VarBW_Air(k,1)
              VarBW_Mulch(k,3)=VarBW_Mulch(k,2)-VarBW_Mulch(k,1)
              VarBW(k,3)=VarBW(k,2)-VarBW(k,1)
              Q(i)=-Width(k)*VarBW(k,3)       
cccz now the irrigation is recorded by the "varbw flux temrs"
cccz                   if (Q(i).gt.0.0) CodeW(i)=-4   ! make sure bc changes if Qn goes > 0 (infiltration)
                  End if    !Q(i) <0
              End if   ! CodeW=4
            End Do
            !End of passing the flux      
            
              tNext(ModNum)=time + PERIOD                     !once irrigated in one hour, next irrigation is in next hour
              IrrigHour=IrrigHour+1                           !this count the number of irrigation in a day
          Else
              IrrigationApplied=IrrigationApplied+1           !when irrigation is stopped on a perticular day of irrigation
                  if (IrrigationApplied.gt.Irrig_times) then
                          tNext(ModNum)=1.E+32
                  else
                      tNext(ModNum)=tApplIrrig(IrrigationApplied)+  ! this set the tnext to the next day of irrigation start time
     !                IrrgStartHour*period
                  end if
          Endif
      End if
      
      Return
20    Stop 'Irrgation data error'

      End
        