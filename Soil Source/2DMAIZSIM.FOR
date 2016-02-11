*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     2DSOIL  - Numerical model of two-dimensional soil environment    *
*                           for field crops                            *
*                                                                      *
*                              Release 04                              *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Include 'public.ins'
      include 'puplant.ins'
	Include 'Puweath.ins'
	Include 'Nitvar.ins'
	
	character*132 RunDat
	integer len, status
	RunFile="Run.Dat"
	Call get_command_argument(1, RunDat, len, status)
	if (len.gt.0) then
	  RunFile=RunDat
	  endif
	  	
	
      Call Initialize()
      Call Get_Grid_and_Boundary()
1     Call Synchronizer()
C       
C SetSurfaceD --  reading daily weather
C SetSurfaceH --  reading hourly weather
      if ((HourlyWeather+DailyWeather.gt.1).or.
     !   ((HourlyWeather+DailyWeather).eq.0))then
             Write (*,*) "error in weather file type"
         
        Endif
      if (HourlyWeather.eq.1) Call SetSurfaceH()
      if (DailyWeather.eq.1) Call SetSurfaceD()                         
      Call SetTDB()
      Call AutoIrrigate()
      Call Mngm()

      Starter=1.0D0
      Call CROP(PCRL,MSW1, NumNP,NumSol,Sink(1),
     !      NumBP,TNext(1),NumMod,Starter)
      Call Carbon_Partitioning()   
      Call Rootgrow()
      Call WaterUptake()  
      Call SoluteUptake()
      Call GasUptake()
     
      Call WaterMover()
      Call SoluteMover()
      Call HeatMover()
      Call GasMover()
      Call SoilNitrogen()
      Call MacroChem()
c      Call MassBl() 
      if (outputSoilYes.gt.0) then 
         Call Output()
        endif

      
      Goto 1
      End         
