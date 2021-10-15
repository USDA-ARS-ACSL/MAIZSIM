C this subroutine will carry out autmatic irrigation for the Agmip Simulations
C the trigger is if the profile has less than 50% available water down to 50 cm.
C irrigation amount is defined in the weather file and is based on rint = average infil rate.
C it is called once a day, at 5:00 am


       Subroutine SetBoundary()
         
       include 'public.ins'
       include 'puWeath.ins'
      
       Parameter (PERIOD =1./24.)
     
       
       character*20    Pond_StartTime, Pond_EndTime
       real    Pond_Head, Day, rPond_StartTime, rPond_EndTime                 
       Integer i, j, l, pos, myHour, GetHour
       character*10 cDate
       
       Common /Auto/ ModNum, Pond_StartTime, Pond_EndTime,
     &     Pond_Head,  rPond_StartTime, rPond_EndTime
       
     
     
       If(lInput.eq.1) then
C
C  Initialize 
C
        open(75,file="D:\2DSOIL_VS9.1\setBound.dat", status='old')
        read (75,*)
        read (75,*) Pond_StartTime, Pond_EndTime, Pond_Head
        NumMod=NumMod+1
        ModNum=NumMod
        
        ! pond_start and stop times are character day:hour
        ! these need to be converted to double based on the calendar
        Myhour=GetHour(Pond_StartTime)
        pos=index(Pond_StartTime,':')
        cDate=Pond_StartTime(1:pos-1)
        rPond_StartTime=JulDay(cDate)+Myhour/24.
        
        Myhour=GetHour(Pond_EndTime)
        pos=index(Pond_EndTime,':')
        cDate=Pond_EndTime(1:pos-1)
        rPond_EndTime= JulDay(cDate)+Myhour/24.
        
        tNext(ModNum) = rPond_StartTime
        
       End If

 ! ponding has ended - set back to original conditions      
       If (time.GT.rPond_EndTime.and.tNext(modNum).LT.1.0e6) then
        !set bc back to -4
         tNext(ModNum)=1.0e6
            Do i=1,NumBP
              n=KXB(i)
              If (CodeW(n).eq.1.and.(n.LT.20)) then
                  CodeW(n)=-4
              EndIf
            EndDo
            tNext(modNum)=1.0e6
       EndIf
       
11     If(abs(time-tNext(ModNum)).lt.0.001*Step) then
C 
c          add your code here
         Do n=1, NumNP
          If (abs(CodeW(n)).eq.4) then
              CodeW(n)=1
              hNew(n)=Pond_Head
          else
              hNew(n)=-150
          endif
         Enddo
c
c         Do i=1,NumBP
c          n=KXB(i)
c          If (CodeW(n).eq.-4) then
c                  CodeW(n)=1
c                  hNew(n)=Pond_Head
c          endif
c        Enddo



c  set tNext large so it does not run again.
           tNext(ModNum) = 1.0e5
           
           
           

       Endif  ! Tnext hourly calcs     
           
           Return
           end
