*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     2DSOIL  - Numerical model of two-dimensional soil environment    *
*                           for field crops                            *
*                                                                      *
*                              Release 3.0.0                              *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      Include 'public.ins'
      include 'puplant.ins'
	Include 'puweath.ins'
	Include 'nitvar.ins'
      Include 'PuSurface.ins'
      	
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
      Call Management()
      Call Irrigation()
      Call PondedIrrigationByHead()
      Call PondedIrrigationByFlux()
	Call TillageProcess
c	Call Drip()

      Starter=1.0D0
      Call CROP(PCRL,MSW1, NumNP,NumSol,Sink(1),
     !      NumBP,TNext(1),NumMod,Starter)
      Call Carbon_Partitioning()   
      Call Rootgrow()
      call Root_Mover()
      Call WaterUptake()
      Call SoluteUptake()
      Call GasUptake()
      
c      Call SetBoundary()
cccz here we start to add surface process
cccz surface mulch (physical perspective)
      Call Surface_Mulch_Process()
cccz surface water
      Call Surface_Water_Balance_Adjustment()
      Call SoilNitrogen()
      Call WaterMover()
      Call SoluteMover()
      Call HeatMover()
      Call GasMover()

C      Call MacroChem()
c      Call Water_Mass_Balance() 
      if (outputSoilYes.gt.0) then 
         Call Output()
        endif

      
      Goto 1
      End         
