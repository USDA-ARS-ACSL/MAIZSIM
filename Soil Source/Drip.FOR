* This subroutine implements drip irrigation

* Note that the first day  for drip irrigation needs to be at least one day later than the 
* initial time in the file in the time folder. This allows the program to initialize the soil properties and fluxes

      Subroutine Drip()
      Include 'public.ins'
      Integer Max_times,Drip_times,Num_nodes,nAppl,applied,
     ! Max_nodes 
      Parameter (Max_times=75, Max_nodes=150)
      Parameter (PERIOD =1./24.)
      Real wAppl,tAppl_start,tAppl_stop,
     !       start_hour(Max_times),
     !       stop_hour(Max_times)
     
      Character*10 Start_Date(Max_times),Stop_Date(Max_times)

      Common /DripC/ tAppl_start(Max_times), tAppl_stop(Max_times), 
     !         wAppl(Max_times), Drip_times,
     !         Num_nodes(Max_nodes), nAppl(Max_times, Max_nodes),
     !         applied,ModNum
      t=time
      If(lInput.eq.1) then
        Open(40,file=DripFile,err=20)
        Read(40,*,Err=20) 
        Read(40,*,Err=20)
        Read(40,*,Err=20) Drip_times
        if (Drip_times.eq.0)  then
           NumMod=NumMod+1
           ModNum=NumMod
           tNext(ModNum)=1e22
        End if
        Read(40,*,Err=20)
        Do i=1,Drip_times
	    Read(40,*,Err=20) Start_Date(i), Start_hour(i), Stop_Date(i), 
     !            Stop_hour(i), wAppl(i), Num_nodes(i)
            Read(40,*,Err=20) 
            Read(40,*,Err=20) (nAppl(i,j),j=1,Num_nodes(i))
        EndDo
        NumMod=NumMod+1
        ModNum=NumMod
        Close(40)
        Do i=1,Drip_times
           tAppl_start(i)=julday(Start_Date(i)) + Start_hour(i)/24.0
           tAppl_stop(i)= julday(Stop_Date(i))  + Stop_hour(i)/24.0
          EndDo
        tNext(ModNum)=tAppl_start(1)
        if (Drip_Times.EQ.0) tNext(ModNum)=1e22
        applied=1
      Endif   ! end lInput.eq.1
      If(Abs(time-tNext(ModNum)).lt.0.001*Step) then

c case for start of dripping
         jj=Applied
         if (time.LT.tAppl_stop(jj)) then
            Do in=1,Num_nodes(jj)
               Do i=1,NumBP
                  n=KXB(i)
                  if (n.eq.nAppl(jj,in)) then
                    Q(nAppl(jj,in))=Q(nAppl(jj,in)) + 
     !                   wAppl(jj)*width(i)*24
                  endif
               EndDo
            EndDo
           tnext(ModNum)=time + PERIOD/8  ! run in 7.5 minute periods
        Else 
         tnext(ModNum)=1e+32
         applied=applied+1
         if (applied.gt.Drip_times) then
           tNext(ModNum)=1.E+32
         else
           tnext(ModNum)=tAppl_start(applied)
         endif
       Endif
      Endif
      Return
20    Stop 'Drip data error'
      End
        
