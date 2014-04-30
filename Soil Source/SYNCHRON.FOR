*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
       
      Subroutine Synchronizer()
      Include 'public.ins'

      character line*40
      Character*10 SDate,FDate
      double precision tFix,tRigid,dt,t,dtmax
*
      If (lInput.eq.1) then
        If (NumMod.eq.-1) then
          Open(42,file=TimeFile,status='old',err=10)
          im=1
          Read(42,*,ERR=10)
          im=im+1
          Read(42,*,ERR=10)
          im=im+1
          Read(42,*,ERR=10) SDate,Step,dtMin,dMul1,dMul2,Fdate
          Read(42,*,ERR=10)
          Read(42,*,ERR=10) DailyOutput,HourlyOutput
          Read(42,*,ERR=10)
          Read(42,*,ERR=10) DailyWeather,HourlyWeather
          Close(42)
          Time=Julday(Sdate)
          TFin=Julday(Fdate)
          dtOpt=Step
          Tinit=Time
          NumMod=0
          RunFlag=1
          Return
10        Call errmes(im,im-1)
        Else
          lInput=0
        Endif
      Else
        t=Time
        tRigid=1.E+30
        Do j=1,NumMod
          if(tNext(j).lt.tRigid) tRigid = tNext(j)
        Enddo
        if (time.gt.51.27) then
          hhh=1
          endif
        dtMax=dmin1(dtMx(1),dtMx(2),dtMx(3),dtMx(4))
        tFix=dmin1(tRigid,tTDB(1),tTDB(2),tTDB(3),tTDB(4),tFin)
        If(Iter.le.3.and.(tFix-t).ge.dMul1*dtOpt) 
     !  dtOpt=dmin1(dtMax,dMul1*dtOpt)
        If(Iter.ge.7)
     !  dtOpt=dmax1(dtMin,dMul2*dtOpt)
        dt=dmin1(dtOpt,tFix-t)
        If ((abs(dt/t).lt.1.E-20).or.(RunFlag.eq.0)) then
	    write(*,*) 'Finished at',time
          Stop
        Endif
        dt=dmin1((tFix-t)/anint((tFix-t)/dt),dtMax)
        If(tFix-t.ne.dt.and.dt.gt.(tFix-t)/2.) dt=(tFix-t)/2.
        Step = dt
      Endif
       WRITE(LINE,100) TIME
*       CALL COUP@(LINE,5,1,24)
*      write(*,*) Time
      Time = Time + Step
*
100    FORMAT(1x,G12.5)
      Return
      End
