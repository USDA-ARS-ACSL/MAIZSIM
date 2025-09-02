      Subroutine Initialize()
      Include 'public.ins'
      Include 'puweath.ins'
      Include 'puplant.ins'
	Character*10 Sowing, Ending, Date4,Date1  ! date1 is dummy for beginDay until we modify the itnerface
	Character*256 RootName,T1,T2
      Character*256 extract_path, path,logFile
      Character*80 Indates,test
      integer iapos, Quote_COUNT, begindate
      
	Integer CurYear, JulDay
c    find root of file name in runfile TODO
c     for writing frequency to output
      Daily=0
      Hourly=0
      AutoIrrigateF=0
      cContentRootM=0.40
      cContentRootY=0.40
      nContentRootM=0.012
      nContentRootY=0.19
      im=1
      logFile=""
      open(9,file=RunFile,status='old', ERR=10)
	read(9,5,err=6)WeatherFile
      im=im+1
	read(9,5,err=6)TimeFile
         im=im+1
	read(9,5,err=6)BiologyFile
         im=im+1
	read(9,5,err=6)ClimateFile
         im=im+1
	read(9,5,err=6)NitrogenFile
         im=im+1
	read(9,5,err=6)SoluteFile
         im=im+1
      read(9,5,err=6)ParamGasFile
         im=im+1
	read(9,5,err=6)SoilFile
         im=im+1
      read(9,5,err=6)MulchFile
         im=im+1
	read(9,5,err=6)ManagementFile
         im=im+1
      read(9,5,err=6)IrrigationFile
         im=im+1
      read(9,5,err=6)DripFile
         im=im+1
	read(9,5,err=6)WaterFile
         im=im+1
	read(9,5,err=6)WaterBoundaryFile
         im=im+1
	read(9,5,err=6)InitialsFile
         im=im+1
	read(9,5,err=6)VarietyFile
         im=im+1
	read(9,5,err=6)GeometryFile
         im=im+1
      read(9,5,err=6)NodeGeomFile   
         im=im+1
      read(9,5,err=6)MassBalanceFile
         im=im+1
	read(9,5,err=6)PlantGraphics
         im=im+1
      read(9,5,err=6)LeafGraphics 	
         im=im+1
       read(9,5,err=6)NodeGraphics   
          im=im+1
C15    Continue 
      read(9,5)ElemGraphics
         im=im+1
C25    Continue
      read(9,5)SurfaceGraphics
         im=im+1
C35    Continue
      read(9,5)FluxGraphics
         im=im+1
      read(9,5)OrganicMatterGraphics
         im=im+1
C45    Continue
      read(9,'(A132)',err=6)MassBalanceFileOut
         im=im+1
      read(9,'(A132)',err=6)MassBalanceRunoffFileOut
         im=im+1
      read(9,'(A132)',err=6)MassBalanceMulchFileOut
         im=im+1
	close(9)
      Path=extract_path(PlantGraphics)
      logFile=trim(Path)//'2DSOIL03.LOG'
      Open(4,file=logFile)
c   end of temporary block
c  These 4 variables are for the iterative solver Orthomin
      ECNVRG=1.0d-6
	ACNVRG=1.0d-6
	RCNVRG=1.0d-6
	MaxItO=200
	AutoIrrigateF=0
      im=1
c    Open and read initials file	
      Open(41, file=InitialsFile, status='old',err=9)
	  read(41,*,err=8)
        im=im+1
	  read(41,*,err=8)
         im=im+1
	  read(41,*,err=8) PopRow,RowSP, PopArea, rowAng, 
     &           xBStem, yBStem, CEC, EOMult
      im=im+1
        read(41,*,err=8)
         im=im+1
        read(41,*,err=8) LATUDE, Longitude, Altitude
         im=im+1
        read(41,*,err=8)  
         im=im+1
cdt 4/2015 fixed error here, variable was AutoIrrigate, added the 'F'        
cccz change here according to "GAS branch"
c        read(41,*,err=8) AutoIrrigateF
        read(41,*,err=8) AutoIrrAmt
         im=im+1
        if (AutoIrrAmt.GT.0)  AutoIrrigateF=1
        
        read(41,*,err=8) 
         im=im+1
        read(41,'(A80)',err=8) inDates
         im=im+1
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
         im=im+1
        read(41,*,err=8)
         im=im+1
        read(41,*,err=8) OutputSoilNo, OutPutSoilYes
         im=im+1
       Close(41)
        if (OutPutSoilNo+OutPutSoilYes.gt.1) then
           Write(*,*) 'error in soil output flag'
           Goto 11
         endif
         if (date1.ne.'00/00/0000') beginDay=JulDay(date1)
         sowingDay=JulDay(Sowing)
         endDay=JulDay(Ending)
         Year=CurYear(Sowing)

c dt
      hNew_org(:)=0.
      ThNew(:)=0.
      Vx(:)=0.
      Vz(:)=0.
      Q(:)=0.
      Tmpr(:)=25.
      Conc(:,:)=0.
      g(:,:)=0.
      QGas(:,:)=0.
*
      g(:,2)=0.
*
      CodeW(:)=0
      CodeS(:)=0
      CodeT(:)=0
      CodeG(:)=0
      Sink(:)=0.
      RTWT(:)=0.
      cSink(:,:)=0.
      gSink(:,:)=0.
	gSink_OM=0.
	gSink_rootY=0.
	gSink_rootM=0.
      gSink_N2O=0.
      SOMMassRatio(:)=1.0
      

* 
      NumMod=-1
      NumBP =0
      NumSol=0
      NSurf =0
      NVarBW=0
      NvarBS=0
      NvarBT=0
      NvarBG=0
      NShoot=0
      NumG=3                                    ! For now, only 3 gas is read (co2, o2, N2O)
    
        KXB(:)=0
        Width(:)=0.
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
 5    format(A256)     
      Return
 6    Write(*,*) 'Error in Run File on line:', im     
      goto 11
 8    Write(*,*) 'Error in initials file on line:',im
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
      
       function extract_path(filename)
       character *256 filename, path, extract_path
       integer :: i, len

       len = len_trim(filename)
       path = ""
    ! Find the last occurrence of the directory separator '\'
       
       do i = len, 1, -1
          if ((filename(i:i) == '\').OR.(filename(i:i) == '/')) then
              path = filename(1:i)
              if (filename(i:i) == '/') then   ! if windows
                  path=path // '/'
               else 
                 path=path // '\'               ! if linux
              end if
             exit
          end if
       end do

       extract_path = path
       end function extract_path
       
        