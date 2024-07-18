* THis module will simulate a ponded head of water where the height changes over time due to infiltration
* The module has two sections when it executes. A ponding depth is applied when the ponding time comes
* otherwise the model has to check if ponding has been applied and in that case needs to 
* check if ponding has ended. If ponding has not ended there is no change to the head. if ponding has ened
* the program has to decrease the head according to infiltration. It checks if the head is zero.
* If the head is zero the program has to change the boundary conditions to no ponding

      
      Subroutine pondedIrrigationByHead()
      Include 'public.ins'
      Include 'nitvar.ins'
      Include 'pusurface.ins'
      
      Parameter (PERIOD= 1./24., pondTimeNum=20)
      Character (len=10) :: startDate(pondTimeNum), endDate(pondTimeNum)
      Character (len=132) :: InString
      Real :: pondedDepth(pondTimeNum)
      Real :: startTime(pondTimeNum), endTime(pondTimeNum)
      Real :: pondedTime, pondedRate  
      Real :: thisTime
      Integer:: pondStarted, pondEnd,jj, irrigationApplied,
     &         startHour(pondTimeNum), endHour(pondTimeNum),
     &         irrigationTimes, modNum
      
      
      Common /ponded_h/pondedDepth, pondedTime, 
     &       pondedRate, pondStarted, pondEnd,modNum, 
     &       startTime, endTime, irrigationApplied, irrigationTimes
      
        ! Intialize module   
      thisTime=time

      
       
      If(lInput.eq.1) then
          pondedTime=0.0
          irrigationApplied=0
          pondStarted=0
          pondEnd=0
          pondingByHead=0
          pondedTimes=0                   ! Number of irrigations
          pondingApplied=0               ! flag to indicate if ponding is applied
          IrrigHour=1                     ! Array count for the irrigation rate
          ! Reading the irrigation details
          Open(40,file=IrrigationFile,err=20)
          Read(40,*,Err=20)
  15      Read (40,'(A132)') InString
          if (InString(1:14).ne.'[Flood_H]') goto 15
          Read(40,*,Err=20)
          Read(40,*,Err=20)
          Read(40,*,Err=20) irrigationTimes   ! Number of days of irrigation
          Read(40,*,Err=20)
          Do i=1,irrigationTimes
              Read(40,*,Err=20) pondedDepth(i), startDate(i), 
     &             startHour(i), endDate(i), endHour(i)
          EndDo
          NumMod=NumMod+1
          ModNum=NumMod

          if (irrigationTimes.EQ.0) then
                tNext(ModNum)=1e22
            else
               
             Do i=1,irrigationTimes
                startTime(i)=julday(startDate(i)) + startHour(i)/24.        ! All the dates irrigation is applied stored in tApplIrrig
                endTime(i)=julday(endDate(i)) + endHour(i)/24.
             EndDo
             tNext(ModNum)=startTime(1)
             irrigationApplied=1
           End If  
          Close(40)
      Endif  ! end reading irrigation details 

!      
!      now check if the ponding is to be applied
         If(Abs(time-tNext(ModNum)).lt.0.001*Step) then
         if (pondStarted.eq.0) then
!Since we can be etiher starting or stopping ponding we need to check if we are starting or stopping
             pondingByHead=1
             pondStarted=1
             pondEnd=0
             jj=IrrigationApplied             
! change boundary conditions        
         Do n = 1, nsurf
           hNew(n) = pondedDepth(jj)
           CodeW(n) = 1
          end do
! now set tNext for the stop time
          tNext(ModNum)=endTime(jj)
! if there are more irrigations then increment counter          
          irrigationApplied=irrigationApplied+1
          else  ! case for ponding ending
          pondStarted=0
          pondEnd=1
          tnext(ModNum)=1e22
             
         End if ! pondStarted=0 
         
        End If  !if next time to apply pondingByHead or end ponding
          
! this is all that is done to initialize the start of pondingByHead           
         
! do within time step loop here while pondingByHead is applied
        if (pondingByHead.eq.1) then
            Qtotal=sum(Q(1:nsurf))/gridwidth*step
          if (pondEnd.eq.1) then
            hnew(1:nsurf)=max(0.0,hnew(1:nsurf)-Qtotal) ! decrease head by infiltration and allow to go to zero
            if (sum(hnew(1:nsurf)).le.0.0) then  !should consider doing this node by node and then changing ponding status when all nodes are below zero
              pondingByHead=0
              pondEnd=0
              do i=1,nsurf   ! need to move this to when head is zero
                CodeW(i)=-4  !soil surface will always be wet so rate will be flux dominated.
                end do
 ! set new times               
              if (irrigationApplied.gt.irrigationTimes) then
                 tnext(ModNum)=1e22
                else 
                 tNext(ModNum)=startTime(irrigationApplied)
              end if  !end else
            end if
          end if
          endif
          
          
      Return       
             
20    Stop 'Irrgation data error'     

      End
