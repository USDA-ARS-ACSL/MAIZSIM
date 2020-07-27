      Subroutine Initialize()
      Include 'public.ins'
      Include 'puweath.ins'
      Include 'puplant.ins'
	Character*10 Sowing, Ending, Date4,Date1  ! date1 is dummy for beginDay until we modify the itnerface
	Character*255 RootName,T1,T2
      Character*80 Indates,test
      integer iapos, Quote_COUNT, begindate
      
	Integer CurYear, JulDay
c    find root of file name in runfile TODO
c     for writing frequency to output
      Daily=0
      Hourly=0
      AutoIrrigateF=0
      open(9,file=RunFile,status='old', ERR=10)
	read(9,'(A132)')WeatherFile
	read(9,'(A132)')TimeFile
	read(9,'(A132)')BiologyFile
	read(9,'(A132)')ClimateFile
	read(9,'(A132)')NitrogenFile
	read(9,'(A132)')SoluteFile
	read(9,'(A132)')SoilFile
c      read(9,'(A132)')MulchFile
	read(9,'(A132)')ManagementFile
      read(9,'(A132)')DripFile
	read(9,'(A132)')WaterFile
	read(9,'(A132)')WaterBoundaryFile
	read(9,'(A132)')InitialsFile
	read(9,'(A132)')VarietyFile
	read(9,'(A132)')GeometryFile
      read(9,'(A132)')NodeGeomFile   
      read(9,'(A132)')MassBalanceFile
	read(9,'(A132)')PlantGraphics
      read(9,'(A132)')LeafGraphics 	
       read(9,'(A132)')NodeGraphics   
C15    Continue 
      read(9,'(A132)')ElemGraphics
C25    Continue
      read(9,'(A132)')SurfaceGraphics
C35    Continue
      read(9,'(A132)')FluxGraphics
C45    Continue
      read(9,'(A132)')MassBalanceFileOut
      read(9,'(A132)')MassBalanceRunoffFileOut
	close(9)
      Open(4,file='2DSOIL03.LOG')
c   end of temporary block
c  These 4 variables are for the iterative solver Orthomin
      ECNVRG=1.0d-6
	ACNVRG=1.0d-6
	RCNVRG=1.0d-6
	MaxItO=200
	AutoIrrigateF=0
c    Open and read initials file	
      Open(41, file=InitialsFile, status='old',err=9)
	  read(41,*,err=8)
	  read(41,*,err=8)
	  read(41,*,err=8) PopRow,RowSP, PopArea, rowAng, 
     &           xBStem, yBStem, CEC, EOMult
        read(41,*,err=8)
        read(41,*,err=8) LATUDE, Longitude, Altitude
        read(41,*,err=8)  
     
cdt 4/2015 fixed error here, variable was AutoIrrigate, added the 'F'        
        read(41,*,err=8) AutoIrrigateF
        read(41,*,err=8) 
        read(41,'(A80)',err=8) inDates
        beginDate=0
        date1='00/00/0000'
        write(test, '(A80)') inDates
        iapos = Quote_Count(inDates)
        if (iapos.eq.4) then
           read(test,*,err=8) Sowing, Ending, TimeStep        
         else if (iapos.eq.6) then
           read(inDates,*,err=8) Date1, Sowing, Ending, TimeStep 
         else 
            goto 8
          endif
  
        read(41,*,err=8)
        read(41,*,err=8)
        read(41,*,err=8) OutputSoilNo, OutPutSoilYes
       Close(41)
        if (OutPutSoilNo+OutPutSoilYes.gt.1) then
           Write(*,*) 'error in soil output flag'
           Goto 11
         endif
         if (date1.ne.'00/00/0000') beginDay=JulDay(date1)
         sowingDay=JulDay(Sowing)
         endDay=JulDay(Ending)
         Year=CurYear(Sowing)
         
         
      Do i=1,NumNPD
        hNew(i)=0.
c dt
        hNew_org(i)=0.
cdtend 8/28/98
        ThNew(i)=0.
        Vx(i)=0.
        Vz(i)=0.
        Q(i)=0.
        Tmpr(i)=0.
        Do j=1,NumSD
          Conc(i,j)=0.
        Enddo
        Do j=1,NumGD
          g(i,j)=0.
        Enddo
*
        Tmpr(i)=25.
        g(i,2)=0.21
*
        CodeW(i)=0
        CodeS(i)=0
        CodeT(i)=0
        CodeG(i)=0
      Enddo
*
      Do i=1,NumElD
        Sink(i)=0.
        RTWT(i)=0.
        Do j=1,NumSD
          cSink(i,j)=0.
        Enddo
        Do j=1,NumGD
          gSink(i,j)=0.
        Enddo
      Enddo
* 
      NumMod=-1
      NumBP =0
      NumSol=0
      NimG  =0
C AD NimG is not used in the entire solution      
      NSurf =0
      NVarBW=0
      NvarBS=0
      NvarBT=0
      NvarBG=0
      NShoot=0
      im=0
      Do i=1,NumBPD
        KXB(i)=0
        Width(i)=0.
      Enddo
*
      lInput=1
      Do i=1,NumModD
        tNext(i)=1.E+31
      Enddo
      Do i=1,4
        dtMx(i)=1.E+31
        tTDB(i)=1.E+31
        Movers(i)=0
      Enddo
      tatm=1.E+31
      
      Return
 8    Write(*,*) 'Error in initials file'
      goto 11     
 9    Write(*,*) 'initials file not found'
      goto 11     
10    Write(*,*)'Run.dat file not found'
11    continue
      End
      
      Function Quote_Count(instring)
      character*80 instring
      integer n,count, Quote_Count
      count=0
      Quote_count=0
      do n=1, len(instring)
       if (instring(n:n).eq."'") then
        count=count+1;
        endif
      end do
      Quote_count=count
      return
      end
      
      
       
        