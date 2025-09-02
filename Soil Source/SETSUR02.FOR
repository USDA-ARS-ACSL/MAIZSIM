CDT 09/10/2014 added CO2 as an environmental variable
C     changed relative humidity to be input as a whole number and
C     convert to fraction internally to make
C     the routine consistent with the hourly one.

C**                                                                   **
C**  this routine taken from GLYCIM                                   **
C**  MSW1-  switch to indicate if daily wet bulb temperatures are     **
C**         available (=1 if YES).                                    **
C**                                                                   **
C**  MSW2-  switch to indicate if hourly wind is available (=1 if YES).**
C**  MSW3-  Daily rain intensities                                    **
C**  MSW4-  Chemical concentrations in rain water                     **
C**  MSW5-  1 for furrow irrigation                                   **
C**  MSW6-  1 if relative humidity is available                       **
C**  MSW7-  1 if Daily CO2 Measurements are available                 **
C**                                                                   **
C**  BSOLAR-factor for changing solar radiation units. Equals         **
C**         radiation in J m-2 divided by radiation in units used in  **
C**         weather data file.                                        **
C**                                                                   **
C**  BTEMP- factor for changing temperature units. Equals change in   **
C**         temperature, in units used in weather data file, that is  **
C**         equivalent to a one oC change.                            **
C**                                                                   **
C**  ATEMP- factor for changing temperature units. Equals temperature **
C**         used in weather data file that is equivalent to zero oC.  **
C**                                                                   **
C**  ERAIN- factor for changing rain units. Equals rainfall in mm     **
C**         day-1 divided by rainfall in units used in weather data   **
C**         file.                                                     **
C**                                                                   **
C**  BWIND- factor for changing wind units. Equals windspeed in       **
C**         km hr-1 divided by windspeed in units used in weather     **
C**         data file.                                                **
C**                                                                   **
C**                                                                   **
C** CO2   -  If daily values are not available read average           **
C**                                                                   **
C**  WIND-  if windspeed is unavailable read in average value for     **
C**         site. windspeed at 2 meters (km hr-1).                    **
C**                                                                   **
C**  (File 1, Line 3 to EOF, free format)                             **
C**      (RH WIND, TWET and TDRY CO2 are optional)                    **
C**                                                                   **
C**  JDAY-  day of year- Julian Date.  Only used for management       **
C**         of the weather file. Program converts calendar date       **
C**         to a day of year relative to 3/1/1900                     **
C**                                                                   **
C**  RI-    daily solar radiation integral (J m-2).                   **
C**                                                                   **
C**  TMAX-  maximum air temperature during the day (oC).              **
C**                                                                   **
C**  TMIN-  minimum air temperature during the day (oC).              **
C**                                                                   **
C**  RAIN-  rainfall (or irrigation) (mm day-1).                      **
C**                                                                   **
C**  WIND-  windspeed at 2 meters (km hr-1).                          **
C**                                                                   **
C**  RH-    Relative Humidity, (%, i.e., 90.2)                        **
C**                                                                   **
C**  IRAV-  avg Rainfall intensity for the day, (cm/hr)               **
C**                                                                   **
C**  TWET-  wet bulb temperature (oC).                                **
C**                                                                   **
C**  TDRY- dry bulb temperature corresponding to TWET (oC).           **
C**                                                                   **
C**  CO2  - CO2 concentration of the atmosphere ppm                   **
C**                                                                   **
C**  Daily values for radiation will be converted to watts m-2 as 
C**   soon as they are read in. 
C**                                                                   **
      Subroutine SetSurfaceD()
      Include 'public.ins'
      Include 'puplant.ins'
      Include 'puweath.ins'
      Include 'PuSurface.ins'
      
      Parameter (PERIOD =1./24.)
      integer jday, m, DayofYear,CurYear, ModNum,ThisYear
      Common /weatherD/ il,im,isol,ModNum, TMIN, TMAX,Gamma_psy,
     &   BEERS(24)
      
      double precision St,t
      real gamma_psy
      character*10 date
      Dimension CLIMAT(20),SDERP(9),SINALT(24),SINAZI(24),HRANG(24),
     &          SARANG(24),
     &           SOLALT(24),SOLAZI(24)
      Dimension xS(NumBPD),iS(NumBPD),kS(NumBPD)
      Data SDERP/0.3964E-0,0.3631E1,0.3838E-1,0.7659E-1,0.0000E0,
     & -0.2297E2,-0.3885E0,-0.1587E-0,-0.1021E-1/,PI /3.1415926/,
     &  TPI /6.2831852/,DEGRAD/0.017453293/,IPERD/24/
CYAP
      If(NSurf.eq.0) then
          RowSP=1.
          height=0.
         endif
CYAPEND
      t=Time
                
      St=0.1D0*Step
      If (lInput.eq.0) goto 11
      
cccz
cccz initialized surface water income  
      do i=1,NumBPD
         Varbw_Air(i,1)=0.0D0     ! cccz this is to record the weather data, because the varBW will change based on runoff
         Varbw_Air(i,2)=0.0D0
         Varbw_Air(i,3)=0.0D0  
         Varbw_Mulch(i,1)=0.0D0     ! cccz this is to record the weather data after mulch model, because the varBW will change based on runoff
         Varbw_Mulch(i,2)=0.0D0
         Varbw_Mulch(i,3)=0.0D0 
      enddo
            

csb: Gas transfer coefficient: surface gas flux change per unit gas content in the soil air at the soil surface (cm/day)
Csb: This value will be different for different gases
      GasTransfCoeff(1)=11920.
      GasTransfCoeff(2)=15400.      
      GasTransfCoeff(3)=12355.
      
      
      NumMod=NumMod+1
      ModNum=NumMod
      tNext(ModNum) = time
     
C
C  First and last days
C     
      JDFRST=idint(t+St)
      JDLAST=idint(tFin+St)
c
      im=160
      il=0

C   IS IT CLIMATE OR WEATHER FOLDER FILE
      Open (5,file=ClimateFile,status='old',ERR=10)
C
C  Read descriptors and conversion constants for weather data\
C
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10) LATUDE
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10) MSW1,MSW2,MSW3,MSW4,MSW5,MSW6,MSW7
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10) BSOLAR,BTEMP,ATEMP,ERAIN,BWIND,BIR
c      ISOL is definined in Grid_bnd
c       assume there is only one solute (NumSol=1)
      ISOL=NumSol*(MSW4)
 
c      the four basic variables are Daily wind, rain intensity daily concentration and co2     
c      we can have average values for the site for
      NCD=4-MSW2-MSW3-MSW7-MSW4*iSol
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read(5,*,ERR=10) (CLIMAT(i),i=1,NCD)
      IF(MSW2.eq.0) WINDA=CLIMAT(1)
      IF(MSW3.eq.0) IRAV =CLIMAT(2-MSW2)
c if MSW4=0 then we need an average value so have to
c  read at least one value from the line
      IF(MSW4.eq.0.AND.NumSol.gt.0) then
      Do i=1,NumSol
          CPREC(i)=CLIMAT(4-MSW2-MSW3-i)
        Enddo
      Endif
      IF(MSW7.eq.0) then
         CO2=CLIMAT(NCD) ! last one is always co2
         GAIR(1)=CO2   
      Endif
      
      close(5)      
C
C Total number of weather data
C     minimum of three columns for radiation, temperature and rainfall (should always be first three) for hourly
C    2*MSW1=TWET, TDRY; MSW2=wind; MSW3=rain rate; MSW4= number of solutes * MSW4, 
C     MSW6=RH, MSW7=CO2  
      NCD=4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6+MSW7
C
C     Nodal numbers of furrow nodes
C
      im1=im
      il1=il
      im=180
      il=0
      If(MSW5.gt.0) then
        Open(6,file='Furnod.dat',status='old',ERR=10)
        im=im+1
        il=il+1
        Read(6,*,ERR=10)
        im=im+1
        il=il+1
        Read(6,*,ERR=10)
        im=im+1
        il=il+1
        Read(6,*,ERR=10) NumFP
        im=im+1
        il=il+1
        Read(6,*,ERR=10)
        im=im+1
        il=il+1
        Read(6,*,ERR=10) (NumF(i),i=1,NumFP)
        im=im+1
        il=il+1
        Read(6,*,ERR=10)
        im=im+1
        il=il+1
        Read(6,*,ERR=10) (hFur(i),i=1,NumFP)
        IFUR=0
      Endif
      im=im1
      il=il1
C
C FIND CORRECT DAY IN WEATHER FILE
C
      Open (5,file=WeatherFile,status='old',ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
      im=im+1
      il=il+1
      Read (5,*,ERR=10)
1111  il=il+1
      Read (5,*,ERR=10,iostat=IOSTAT) iDum, date
      if (iostat<0) then
         write(*,*) 'Premature end of weather file'
         stop
      end if
      
      MDAY=julday(date)
cdt changed from LE to LT since jday was always gt. time
      If (MDAY.LT.JDFRST) GO TO 1111
C
C   CALCULATE A CLOUD COVER FACTOR FOR THIS LATITUDE 
C
      If (LATUDE.LE.25.0) THEN
          CLDFAC = 0.45 - (LATUDE*0.004)
      Else
          CLDFAC = 0.30 + (LATUDE*0.002)
      Endif
c
c..................... Routine calculations
c..................... Calculations made daily
   11 If(lInput.eq.1.OR.idint(time+St).eq.JDAY+1) then
c
c ZEROING DAILY ARRAYS
c
        Do i=1,24
          RINT(i)=0.
          RADINT(i)=0.
          WATTSM(i)=0.
        Enddo
c
c Read meteodata
c
        il=il-1
        BACKSPACE(5)
        il=il+1
        Read(5,*,ERR=10, iostat=IOSTAT) JDum,date,(CLIMAT(i),i=1,NCD)
        if (iostat<0) then
           write(*,*) 'Premature end of weather file'
           stop
        end if
           JDAY=Julday(date)
c    SInce some routines need the day of the calendar year we need a variable to hold this
C    since the julian day is referenced to a time longer in the past
        ThisYear=CurYear(date)
        write (date,'("01/01/",i4.4)') ThisYear  
        DayOfYear=JDAY-julday(date)
C
C   ADJUST UNITS FOR WEATHER DATA
C
c
        RI = CLIMAT(1)*BSOLAR
        TMAX = (CLIMAT(2) - ATEMP)/BTEMP
        TMIN = (CLIMAT(3) - ATEMP)/BTEMP
c   dt have to fix this later as model cannot handle freezing temperatures
        TMIN=max(TMIN,2.0)       
        TMAX=max(TMAX,TMIN)
        RAIN = CLIMAT(4)*ERAIN

        If (MSW2.GT.0) then 
          WIND = CLIMAT(5)*BWIND
        Else
          WIND=WINDA*BWIND
        Endif
        If (MSW1.GT.0) then
          TWET = (CLIMAT(5+MSW2) - ATEMP)/BTEMP
          TDRY = (CLIMAT(6+MSW2) - ATEMP)/BTEMP
        Endif
        If (MSW3.GT.0) then
          IR=CLIMAT(4+2*MSW1+MSW2+MSW3)*BIR
        Else
          IR=IRAV*BIR
        Endif
        If (MSW4.GT.0.AND.NumSol.ne.0) then
          Do i=1,NumSol
            CPREC(i)=CLIMAT(4+2*MSW1+MSW2+MSW3+i)
          Enddo
        Endif
        AutoIrrigAmt=IR   ! used in AutoIrrigation if carried out.
CDT

      If(MSW6.gt.0) then
c        Rel_Humid=Climat(4+2*MSW1+MSW2+MSW3+ISOL+MSW6)TMIN=max(TMIN,2.0)   
         Rel_Humid=min(Climat(4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6)
     &     /100.0,0.98)
          Rel_Humid=max(0.0,Rel_Humid)
c         Rel_Humid=CLIMAT(5)
       endif
      If(MSW7.gt.0) then
          CO2=Climat(4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6+MSW7)          !For gas mover, the surface co2 concetration is read from here [ppm]
          GAIR(1)=CO2                                                 !only for 1 gas now. we need to do that for all the gases
      EndIf
      
           GAIR(2)=209000    !this is the atmospheric O2 concentration in ppm  
csb O2 content of air=20.95% by volume, ie; 100 parts of air=20.95 parts of O2
csb 1 part of air=0.2095 parts of O2. in 10^6 parts of air=0.2095*10^6=209000 ppm
           GAIR(3)=0.332    !this is the atmospheric N2O concentration in ppm  [IPCC 2021,Technical summary]
csb N2O concentration in atm 332 ppb=0.33 ppm 
          

Cdt
c.........Radiation submodel
C
C  CALCULATE SOLAR DECLINATION
C
        XLAT = LATUDE*DEGRAD
        DEC = SDERP(1)
        Do I = 2,5
          N = I - 1
          J = I + 4
          D11 = N*0.01721*DayOfYear
          DEC = DEC + SDERP(I)*SIN(D11) + SDERP(J)*COS(D11)
        Enddo
        DEC = DEC*DEGRAD
C
C  CALCULATE SOME PARTS OF THE SOLAR RADIATION EQUATIONS THAT
C         OCCUR REPEATEDLY LATER
C
        D12 = SIN(XLAT)*SIN(DEC)
        D13 = COS(XLAT)*COS(DEC)
        D14 = D12 + D13
C
C  CALCULATE DAYLENGTH
C
        DAYLNG = ACOS((-0.014544 - D12)/D13)*7.6394
C       7.6394 = 180/3.1416/360*24*2
C
C  CALCULATE SOLAR RADIATION INCIDENT ON TOP OF THE EARTH'S
C         ATMOSPHERE AT SOLAR NOON
C
        RADVEC = 1 + (0.01674*SIN((DayOfYear - 93.5)*0.9863*DEGRAD))
        WATATM = 1325.4*D14/(RADVEC*RADVEC)
C
C  CALCULATE AN ATMOSPHERIC TRANSMISSION COEFFICIENT FOR THIS
C         LATITUDE AND TIME OF YEAR. ***
C
        If (DayOfYear.LT.145) then
          ATRANS = 0.68 + ((5.25E-5*LATUDE) - 0.1E-3)*(145 - DayOfYear)
        Else If (DayOfYear.LE.237) then
          ATRANS = 0.68
        Else
          If (LATUDE.LE.30.0) then
            D15 = (LATUDE*5.5E-5) - 0.1E-3
          Else
            D15 = 0.65E-3 + (LATUDE*3.0E-5)
          Endif
          ATRANS = 0.68 + D15*(DayOfYear - 237)
        Endif
C
C  CALCULATE POTENTIAL DIRECT + DIfFUSE RADIATION INCIDENT ON
C         CROP AT SOLAR NOON
C
        WATPOT = WATATM*0.5*(0.93 - (0.02/D14) + ATRANS**(1/D14))
C
C  CALCULATE ACTUAL RADIATION INCIDENT ON CROP AT SOLAR NOON
C  GIVEN DAILY INTEGRAL AND ASSUMING RADIATION FLUX DENSITY
C  VARIES AS A HALF SINE WAVE OVER THE PHOTOPERIOD
C
        WATACT = RI*4.363E-4/DAYLNG
C       4.363E-4 = 3.1416/3600/2
C
C  CALCULATE CLOUD COVER. ***
C
        WATRAT = WATACT/WATPOT
        If (WATRAT.GE.1.0) then
          CLOUD = 0.0
        Else
          CLOUD =(CLDFAC - SQRT((CLDFAC*CLDFAC) + (1.52*(1 - WATRAT))))
     &   /(-0.76)
C        1.52 = 4*0.38
C        0.76 = 2*0.38
          If (CLOUD.GE.1.0) CLOUD = 1.0
        Endif
C
C       DUSK AND DAWN
C
        DAWN = 12.0 - (DAYLNG/2.0)
        DUSK = 12.0 + (DAYLNG/2.0)
C
C  CALCULATE TOTAL RADIATION INCIDENT ON THE CROP AT EACH TIME
C
C  DIVIDE PHOTOPERIOD INTO IPERD EQUAL INCREMENTS AND CALCULATE
C  HOUR ANGLES BETWEEN THE MIDPOINTS AND SOLAR NOON
C
        IHPERD = IPERD/2
        Do 20, J = 1,IHPERD
          HRANG(J) = PI/12.0*(12 - J + 0.5)
 20     HRANG(IPERD - J + 1) = HRANG(J)
        DDIf = DAYLNG - (2.0*IfIX(DAYLNG/2))
        IUP = 13 - IfIX(DAYLNG/2)
        IDN = 12 + (13 - IUP)
        If (.NOT.((DDIf.GT.0.0).OR.(DDIf.LT.0.0)))then
          IDAWN = IUP
          IDUSK = IDN
        Else
          IDAWN = IUP - 1
          IDUSK = IDN + 1
        ENDIf
C
C  CALCULATE TOTAL RADIATION INCIDENT ON THE CROP AT EACH TIME
C
        Do I = 1,IPERD
          TIMH = I - 0.5
          WATTSM(I) = 0.0
          If(I.GE.IUP.AND.I.LE.IDN)
     &      WATTSM(I) = WATACT*SIN(PI*(TIMH - DAWN)/DAYLNG)
        Enddo
        If ((DDIf.GT.0.0).OR.(DDIf.LT.0.0)) then
          TEND = (1 - COS((PI/DAYLNG)*(DDIf/2)))
           WATTSM(IDAWN) = TEND*RI*1.389E-4
C          1.389E-4 = 1/2)/3600
           WATTSM(IDUSK) = WATTSM(IDAWN)
        ENDIf
c DT 10/12/2021 sometimes due to DST issues and weather station clocks,
c  radiation between dawn and dusk  can be 0 so need to adjust WATTSM for this
       Do I = 1, IPERD
        if (I.ge.IDAWN.and.I.le.IDUSK) then
          if (WATTSM(I).le.0.0) then
           WATTSM(I)=WATTSM(I)+0.1
           endif
        endif
       enddo
c..................... Temperature and vapour pressure submodel
C
C  READ NEXT DAY'S DATA FROM WEATHER FILE ***
C
        il=il+1
        Read (5,*,ERR=10,iostat=IOSTAT) MDAY, date,(CLIMAT(i),i=1,NCD)
        if (iostat<0) then
          write(*,*) 'Premature end of weather file'
          stop
        endif

c       Minimum tomorrow temperature
        MDAY=Julday(date)
        If (MDAY.LT.JDLAST) then
          TMINT = MAX((CLIMAT(3) - ATEMP)/BTEMP,2.0)
        Else
          TMINT = max(TMIN, 2.0)
        ENDIf
C
C  YESTERDAY DUSK TEMPERATURE
C
        If (lInput.GT.0) then
          TDUSKY = (TMAX + TMIN)/2.0
        Else
          TDUSKY = TDUSK
        Endif
C
C  CALCULATE THE TIME AFTER DAWN IN HOURS WHEN MAXIMUM 
C  TEMPERATURE IS REACHED
C
        D20 = 0.0945 - (WATACT*8.06E-5) + (TMAX*6.77E-4)
        D21 = TMAX/D20/WATACT
        If (D21.GE.1.0) D21 = 1.0
        TMAXHR = DAYLNG/PI*(PI - (ASIN(D21)))
C
C  CALCULATE SOME PARTS OF THE TEMPERATURE EQUATIONS THAT OCCUR
C  REPEATEDLY LATER
C
        D22 = ((TMAX - TMIN)/2.0)
        D23 = (PI/TMAXHR)
        D24 = 1.5*PI
C
C  CALCULATE AIR TEMPERATURE AT DUSK
C 
        TDUSK = (D22*(1.0 + SIN((D23*DAYLNG + D24)))) + TMIN
C
C  CALCULATE SOME ADDITIONAL PARTS OF THE TEMPERATURE EQUATIONS
C  THAT OCCUR REPEATEDLY LATER
C
        XTEMP = 5.0
        If (TMIN.LT.TDUSKY) then
          D25 = TDUSKY - TMIN + XTEMP
          D26 = ALOG(D25/XTEMP)/(2.0*DAWN)
        Else
          D27 = (TMIN - TDUSKY)/(2.0*DAWN)
        Endif
        If (TMINT.LT.TDUSK) then
          D28 = TDUSK - TMINT + XTEMP
          D29 = ALOG(D28/XTEMP)/(2.0*DAWN)
        Else
          D30 = (TMINT - TDUSK)/(2.0*DAWN)
        Endif
C
C  CALCULATE AIR TEMPERATURE AT EACH TIME
C
        Do 70, I = 1,IPERD
          TIMH = I - 0.5
          If (TIMH.GE.DAWN.AND.TIMH.LE.DUSK) then
            TAIR(I) = TMIN +
     &      (D22*(1.0 + SIN((D23*(TIMH - DAWN)) + D24)))
c             TAIR(I)=TMAX
          Else If(TIMH.LT.DAWN) then
            If (TMIN.LT.TDUSKY) then
              TAIR(I) = TMIN - XTEMP + (D25/EXP(D26*(DAWN + TIMH)))
            Else
              TAIR(I) = TDUSKY + (D27*(DAWN + TIMH))
            Endif
c            TAIR(I)=TMIN
c         else
          Else If(TIMH.GT.DUSK) then
            If (TMINT.LT.TDUSK) then
              TAIR(I) = TMINT - XTEMP + (D28/EXP(D29*(TIMH -DUSK)))
            Else
             TAIR(I) = TDUSK + (D30*(TIMH - DUSK))
             TAIR(I)=TMIN
            Endif
          Endif
 70     Continue

C
C  If DAILY WET BULB TEMPERATURES ARE NOT AVAILABLE, ASSUME
C  MINIMUM TEMPERATURE = DEW POINT TEMPERATURE AND CALCULATE
C  WATER VAPOR PRESSURE.  ASSUME A VALUE FOR THE PSYCHROMETRIC
C  "CONSTANT"
C
cdt added msw6 to if statement here

        If (MSW1.NE.1.and.MSW6.eq.0) then
          AVP = 0.61*EXP((17.27*TMIN)/(TMIN + 237.3))
          GAMMA_psy = 0.0645
        Else
cdt
           if(MSW1.eq.1) then      
C
C  SINCE DAILY WET BULB TEMPERATURES ARE AVAILABLE CALCULATE
C  SATURATION VAPOR PRESSURE AT THE WET BULB TEMPERATURES
C
          SVPW = 0.61*EXP((17.27*TWET)/(TWET + 237.3))
C
C  CALCULATE THE HUMIDITY RATIO, D31, LATENT HEAT OF EVAPORATION,
C  D32 AND THE PSYCHROMETRIC "CONSTANT"
C
          D31 = 0.622*(SVPW/(101.3 - SVPW))
          D32 = 2500.8 - (2.37*TWET)
          GAMMA_psy = 0.62*(1.006 + (1.846*D31))
     &   /((0.622 + D31)*(0.622 + D31)*D32)*101.3
C
C  CALCULATE ACTUAL WATER VAPOR PRESSURE
C
          AVP = SVPW - (GAMMA_psy*(TDRY - TWET))
C _todo   Calculate Relative humidity here if it is not input
        Endif !if MSW!.eq.1

CDT        
            if(MSW6.eq.1) then
C if relative humidity is available then calculate saturation vapor
C pressures at Tmin and Tmax and average the two values
C Also calculate TMean        
                TMean=(TMax+TMin)/2.0
                SVPW1= 0.61*EXP((17.27*TMin)/(TMin + 237.3))
                SVPW2= 0.61*EXP((17.27*TMax)/(TMax + 237.3))
                SVPW=(SVPW1+SVPW2)/2.0

                D31 = 0.622*(SVPW/(101.3 - SVPW))
                D32 = 2500.8 - (2.37*TMean)
                GAMMA_psy = 0.62*(1.006 + (1.846*D31))
     &          /((0.622 + D31)*(0.622 + D31)*D32)*101.3
                AVP= Rel_Humid*SVPW
          endif   ! end MSW6.eq.1
       
        Endif   !MSW1.NE.1.and.MSW6.eq.0
CDT
C
C  CALCULATE WATER STAURATION VAPOR PRESSURE AT AIR TEMPERATURE
C  AND THE SLOPE OF THE RELATIONSHIP AT THAT POINT
C
        Do 80, I = 1,IPERD
          SVPA = 0.61*EXP((17.27*TAIR(I))/(TAIR(I) + 237.3))
          DEL(I) = (0.61*EXP((17.27*(TAIR(I) + 1.))
     &   /(TAIR(I) + 1.0 + 237.3))) - SVPA
C
C  CALCULATE WATER VAPOR PRESSURE DEFICIT
C
          VPD(I) = SVPA - AVP
          If (VPD(I).LE.0.0) VPD(I) = 0.0
 80     Continue
c....................Light interception submodel
        If(NShoot.ne.0) then
C
C  CALCULATE SOLAR ALTITUDE AND SOLAR AZIMUTH FOR EACH TIME
C
          Do 30 I = 1,IPERD
            SINALT(I) = D12 + (D13*COS(HRANG(I)))
            SOLALT(I) = ASIN(SINALT(I))
            SINAZI(I) = max(-1.D0,(-COS(DEC)*
     &        SIN(HRANG(I))/COS(SOLALT(I))))
 30       Continue
          D18 = 0.0
          Do 60, I = 1,IHPERD
            SOLAZI(I) = PI + ASIN(SINAZI(I))
            if(.NOT.((D18.GT.0.0).OR.(D18.LT.0.0))) then
              D19 = HRANG(I) - 0.01
              TESAZI = PI + ((ASIN(-COS(DEC)*SIN(D19))
     &       /(COS(ASIN(D12 + (D13*COS(D19)))))))
              If (TESAZI.LT.SOLAZI(I)) then
                SOLAZI(I) = PI - SOLAZI(I)
              Else
                D18 = 1.0
              Endif
            Endif
            SOLAZI(IPERD - I + 1) = TPI - SOLAZI(I)
 60       Continue
C
C  CALCULATE PROPORTION OF TOTAL RADIATION THAT IS DIFFUSE AT
C  EACH TIME
C
          Do 40 I = 1,IPERD
            D16 = SINALT(I)
            If (D16.LT.0.01) then
              DIFWAT(I) = 1.0
            Else
              D17 = (0.93 - (0.02/D16) - ATRANS**(1/D16))
     &       /(0.93 - (0.02/D16) + ATRANS**(1/D16))
              If (CLOUD.GT.0.0) then
                DIFWAT(I) = 1.0 - ((1.0 - D17)*(1.0 - CLOUD)
     &         /(1.0 - ((CLDFAC+(0.38*CLOUD))*CLOUD)))
              Else
                DIFWAT(I) = D17
              Endif
              If (DIFWAT(I).LT.0.0) DIFWAT(I) = 0.0
              If (DIFWAT(I).GT.1.0) DIFWAT(I) = 1.0
            Endif
C
C     *** CALCULATE PHOTOSYNTHETICALLY ACTIVE RADIATION (PAR)
C         INCIDENT ON THE CROP ASSUMING THAT
C         PAR = 0.43 DIRECT + 0.57 DIFFUSE RADIATION. ***
C
           PAR(I) = WATTSM(I)
     &          *((DIFWAT(I)*0.57) + ((1.0 - DIFWAT(I))*0.43))
40         CONTINUE
C
C
C  If CLOUD COVER IS COMPLETE, OMIT CALCULATIONS PERTAINING
C  TO DIRECT RADIATION INTERCEPTION
C
          If (CLOUD.GE.1.0) then
            Do I = 1,IPERD
              DIRINT(I) = 0.0
            Enddo
          Endif
C
C  CALCULATE ANGLE BETWEEN ROW AND SOLAR AZIMUTH FOR EACH TIME
C
          Do I = 1,IPERD
            SARANG(I) = ABS(SOLAZI(I) - (ROWANG*DEGRAD))
            If (SARANG(I).GT.3.*PI/2.) SARANG(I)
     &       = ABS(SARANG(I) - TPI)
            If (SARANG(I).GT.PI/2.) SARANG(I)
     &       = ABS(SARANG(I) - PI)
C
C  CALCULATE PROPORTION OF DIRECT RADIATION INTERCEPTED BY ROWS
C  OF PLANTS ASSUMING THEY ARE OPAQUE CYLINDERS
C
            SHADOW(I) = HEIGHT/SIN(ATAN(TAN(SOLALT(I))
     &     /SIN(SARANG(I))))
            DIRINT(I) = SHADOW(I)/ROWSP
            If (DIRINT(I).GT.1.0) DIRINT(I) = 1.0
            If (DIRINT(I).LT.0.0) DIRINT(I) = 0.0
          Enddo
C
C  CALCULATE PROPORTION OF DIfFUSE RADIATION INTERCEPTED BY ROWS
C  OF PLANTS ASSUMING THEY ARE OPAQUE CYLINDERS
C  FIRST DIVIDE ROWSPACING INTO 10 EQUAL PARTS AND FIND THE
C  MIDPOINTS FROM ROW TO MID ROW
C
          DIfFIN = 0.0
          Do I = 1,10
            ROWINC(I) = ROWSP/20.0*(I - 0.5)
C
C  CALCULATE PROPORTION OF SKY OBSCURED BY "OPAQUE" ROWS AT EACH
C  POSITION FROM ROW TO MID-ROW
C
            DIfINT(I) = (ATAN(HEIGHT/2.0/(ROWSP - ROWINC(I)))
     &     +ATAN(HEIGHT/2.0/ROWINC(I)))/1.5708
            If (DIfINT(I).GT.1.0) DIfINT(I) = 1.0
C
C  INTEGRATE DIfFUSE RADIATION INTERCEPTION ACROSS THE ROW. ***
C
            DIfFIN = DIfFIN + DIfINT(I)
          Enddo
c
          DIfFIN = DIfFIN/10.0
C
C  CALCULATE LEAF AREA PER UNIT GROUND AREA COVERED BY CROP
C  CANOPY
C
C

C  CALCULATE AN EFFECTIVE LCAI ALLOWING FOR THE FACT THAT LIGHT
C  AT A LOW ANGLE TRAVERSES MORE LEAF LAYERS TO REACH THE SOIL
C


           LCAI=LAI
           
          Do 50, I = 1,IPERD
            If (SOLALT(I).GT.0.0) then
              ELCAI = LCAI/(SIN(ACOS(COS(SARANG(I))
     &       *COS(SOLALT(I)))))
              If (SHADOW(I).GT.ROWSP)
     &           ELCAI = ELCAI*SHADOW(I)/ROWSP
C
C  CALCULATE BEERS LAW CORRECTION FOR RADIATION INTERCEPTION
C
              If((ELCAI*CEC).LE.88.0) then
                BEERS(I) = (1 - EXP(-ELCAI*CEC))
              Else
                BEERS(I)=1.0
              Endif

C  CALCULATE PROPORTION OF PHOTOSYNTHETICALLY ACTIVE RADIATION
C   (PAR) INTERCEPTED BY THE CROP CANOPY ASSUMING THAT
C    PAR = 0.43 DIRECT + 0.57 DIFFUSE RADIATION. ***

      PARINT(I) = WATTSM(I)*((DIFFIN*0.57*DIFWAT(I))
     &     + (DIRINT(I)*0.43*(1.0 - DIFWAT(I))))
     &     /PAR(I)*BEERS(I)
C
C  CALCULATE PROPORTION OF TOTAL RADIATION INTERCEPTED BY THE
C  CROP CANOPY
C
              RADINT(I) = ((DIfFIN*DIfWAT(I)) + (DIRINT(I)
     &         *(1.0 - DIfWAT(I))))*BEERS(I)
            Else If (WATTSM(I).GT.0.0) then
              BEERS(I) = (1 - EXP(-LCAI*CEC))
              PARINT(I) = COVER*BEERS(I)
              RADINT(I) = COVER*BEERS(I)
            Else
               PARINT(I) = 0.0
               RADINT(I) = 0.0
            Endif
 50       Continue
        Endif
C
c..................End of radiation interception submodel
C
C
c..................Daily calculations for the precipitation
C
! dtRain is time in days to it takes for this amount to rain using the average hourly rainfall rate (IR)
! check to see if this amount will be reached in the current hour
! if it is too much for the hour, then the excess is carried over to the next hour
! note that RINT has to be a 24 hour rainfall rate - cm/day. it is input as cm over the hour.

        If(Rain.ge.0.) then
          DtRain=RAIN/IR 
          Do i=24,1,-1
            If(DtRain.ge.PERIOD*24.0) then
              RINT(i)=IR*24.0
              dtRain=dtRain-PERIOD*24.0
            Else
cdt changed this from dtrain/ir 11/14/96
              RINT(i)=DtRain*IR*24.0
              goto 55
            Endif
          Enddo
        Else
          QF=QF-RAIN
          Do i=1,NumFP
            k=NumF(i)
            CodeW(k)=1
            hNew(k)=hFur(i)
          Enddo
          IfUR=1
        Endif
c........................End of daily calculations
      Endif
C
C     Further we have hourly calculations
C
55    If(abs(time-tNext(ModNum)).lt.0.001*Step) then
          
c 55    If(abs(t-tAtm).gt.St) goto 1001
        ITIME=int(MOD(t+St,1.)/PERIOD+1)
C
C  FOR DAYLIGHT PERIODS:
C  CALCULATE NET UPWARD LONG-WAVE RADIATION ***
C
1000    Continue
        If (WATTSM(ITIME).GT.0.0) then
          RNLU = (1.6E-3*WATRAT*(100.0 - TAIR(ITIME)))*697.6
        Else
          RNLU=0.
        ENDIf
C
C  NO SHADE IF NO PLANT
C
        If(NShoot.eq.0) then
          SHADE = 0.
          COVER = 0.
c          RowSP = 1.0
cdt there was a problem here since it overwrites rowsp where it was
c initialized in crop - need to move the above code into a crop module
c even if a dummy one
       ENDIf
C
C  CALCULATE POTENTIAL WATER EVAPORATION RATE FROM THE SOIL
C  SURFACE
C
C  THE TOTAL RADIATION FALLING ON THE SOIL IS ASSUMED TO BE
C  CONCENTRATED ON THE EXPOSED SOIL CELLS AND AN EQUIVALENT
C  TOTAL RADIATION IS CALCULATED
C
         If (Cover.lt.0.97) then
C
             D11 = WATTSM(ITIME)

          Else
              D11=0.
          Endif
c
c  MARCH ALONG THE SURFACE TO FIND SOIL-ATMOSPHERE BOUNDARY NODES
c  AND THEIR LATERAL COORDINATES
c
        ic=0
        Do k=1,NumBP
          i=KXB(k)
          If(abs(CodeW(i)).eq.4) then
            ic=ic+1
            xS(ic)=x(i)
            iS(ic)=i
            kS(ic)=k
          Endif
        Enddo
        Call SORT03(ic,xS,iS,kS)
cccz Dennis: your basic assumption is horizontal distance from origin in the input is not sorted
cccz so if you use KS, then use the surface node index / node index based on KS
C
C initialize boundary values to 0
         Do i=1, NumBP
         do j=1, 4
           if (j.lt.4) VarBW(i,j)=0.0
           VarBT(i,j)=0.0
         Enddo
      Enddo

cccz -------------------------------------------------------------------------------------
cccz Dennis, please check here
c    The revision here is similar to the one in hourwea,
c    I put "If(PSh.gt.0.) then" after the calculation of "ESO"
c    the reason is the ELSE part also use ESO

cccz make some initialization
       ESO=0.0D0
       EPO=0.0D0      
        Do 90 i=1,ic
cDT 4-22-2011 assume 10% of et is minimum from a covered node
cccz this PSh expression is different from the hour.for one, and it has stronger lower boundary
          PSh=max(0.0,FSh(i,xBSTEM,SHADE,xS,ic))
C
C  CALCULATE THE ALBEDO OF THE EXPOSED SOIL CELLS
C
cccz parameter difference here, comparing to the hour.for
              LAMDAS = 0.3 - 0.5*ThNew(iS(i))
C
C  CALCULATE NET RADIATION ON THE EXPOSED SOIL CELLS
C
              RNS = ((1.0 - LAMDAS)*D11) - RNLU
              If(RNS.LE.0.0) RNS = 0.0
C
C  CALCULATE POTENTIAL EVAPORATION RATE (g m-2 hr-1)FROM EXPOSED SOIL CELLS
C
C note that Penman's original wind function was F(u)=0.26(1+0.54u)
C the 0.54 is a function of temperature and pressure
C  where windspeed was in m s-1 and et mm d-1 cm-2 and e in mb
C  109.375 = 0.26 * (10000*10)/24/10 (0.26 has probably been truncated)
C
C note the original Penman eqn wind function was 0.26*(1+u/160) where
C windspeed (u) was km/day here we use km/hour so you divide 1/160 by 24 
C gives you the .149

C See eqn 5.11 of http://mepas.pnnl.gov/mepas/formulations/source_term/5_0/5_13/5_13.html
              D12 = max(1.0,2.0-(2.0*COVER))
              If (D12.GE.1.0) D12 = 1.0
              ESO = ((DEL(ITIME)/GAMMA_psy*RNS*3600.0/(2500.8
     &     - (2.3668*TAIR(ITIME))))
     &     + (VPD(ITIME)*109.375*(1.0 + (0.149*WIND*D12))))
     &      /((DEL(ITIME)/GAMMA_psy) + 1.0)

c
c  IF THE NODE IS EXPOSED THEN
c
C
C   PROPORTIONAL TO THE AREA EXPOSED
          If(PSh.gt.0.) then
C
              VarBW(kS(i),2)=PSh*24.*ESO/10000.
c
c  ELSE IF THE NODE IS COVERED
cDT 4-22-2011 assume 10% of et is minimum from a covered node
c
          Else
            VarBW(kS(i),2)=0.05*24.*ESO/10000.
          Endif
cccz -----------------------------------------------------------------------------------
          
c                Endif for the surface nodes
 90     Continue
C
C........................ End evaporation 
C
C........................ Transpiration
C
        If (NShoot.ne.0) then
C
C  CALCULATE A CROP SURFACE ROUGHNESS PARAMETER
C
        If (ITIME.EQ.1) then
           ROUGH = 1.0/(ABS((HEIGHT/ROWSP) - 0.5) + 0.5)
           If (ROUGH.LT.1.0) ROUGH = 1.0
        Endif
C
C   ON STILL HOT DAYS, CONVECTION CURRENTS AT CROP LEVEL ASSURE
C   SOME AIR MOVEMENT
C
        WINDL = WIND
        If (TAIR(ITIME).GT.25.0.AND.WINDL.LE.0.36) WINDL = 0.36
C
C   CALCULATE POTENTIAL TRANSPIRATION RATE FOR THE CROP
C
C   ASSUME A VALUE FOR THE ALBEDO OF THE CROP
C
        LAMDAC = 0.23
C
C   THE TOTAL RADIATION INTERCEPTED BY THE CROP IS ASSUMED TO BE
C   SPREAD UNIfORMLY OVER THE AREA COVERED BY THE CROP AND AN
C   EQUIVALENT TOTAL RADIATION IS CALCULATED
C
C   CODE TO ACCOUNT FOR THE CASE WHERE THE PLANT
C   HAS NOT YET EMERGED AND COVER =0
        If (COVER.le.0.) then
          WATPL=0.0
        Else
          WATPL = WATTSM(ITIME)*RADINT(ITIME)/COVER
        Endif
C
C   CALCULATE NET RADIATION ON THE CROP
C
        RNC = ((1.0 - LAMDAC)*WATPL) - RNLU
        If (RNC.LE.0.0) RNC = 0.0
C
C   CALCULATE POTENTIAL TRANSPIRATION RATE FOR CROP ALLOWING FOR
C   INCOMPLETE GROUND COVER units are grams per m2 per hour
C
        EPO = ((DEL(ITIME)/GAMMA_psy*RNC*3600.0/(2500.8
     &   - (2.3668*TAIR(ITIME))))
     &   + (VPD(ITIME)*109.375*(ROUGH +(0.149*WINDL))))
     &   /((DEL(ITIME)/GAMMA_psy) + 1.0)
C
C   CODE ADDED TO PREVENT DIVISION BY o
C
        If (COVER.le.0.0.OR.EPO.LE.0.0) EPO=1.0E-4
C
C     *** CALCULATE POTENTIAL TRANSPIRATION RATE FOR AN AVERAGE PLANT **
C
cdt I took out etcorr, the ET correction factor
c  units are grams per plant per day adjusted for cover
        EO = 24.*EPO*cover/(poprow/rowsp*100)
        If(COVER.le.0.0) then
          TPot=0.
        Else
          TPot=EO/COVER
        Endif
      Endif
c................... End Transpiration .................................
c................... Precipitation and irrigation
      Do i=1,NumBP
      n=KXB(i)
      k=CodeW(n)
      If(K.eq.4.or.K.eq.-4) then
         VarBW(i,1)=RINT(ITIME)
       If(NumSol.ne.0) then
          do j=1,NumSol
cdt was varbs(n,j)
          VarBS(i,j)=CPREC(j)
          Enddo
       Endif
*         
         VarBW(i,3)=VarBW(i,2)-VarBW(i,1)
         Q(n)=-Width(i)*VarBW(i,3)        ! positive if rainfall-> going into the soil
         if (Q(n).gt.0.0)  then
           CodeW(n)=-4
          Endif
       Endif
       Enddo
c................... End of the precipitation and irrigation
C ............... Heat Balance is calculated here
C
C
CYAP
CDT - added FELWR based on changes from temperature paper
C FCSH is the constant for sensible heat transfer 

      FCSH=(4.0E-3+1.39E-3*Wind*max(1.0,2.0-(2.0*COVER)))
      FCSH=FCSH*4.1856*60.0*24.0
      
C FELWR is the coefficient for radiant heat lost when the 
C soil is warmer than the air.      
      FELWR=6.8E-03+9.0E-05*Tair(Itime)
      FELWR=FELWR*4.1856*60.0*24.0
* units are Joules cm-2 d-1 C-1
CYAPEND
C units are changed from cal cm-2 min-1 C-1 -->J cm-2 d-1 C-1
      cWat=4.1856
      TairN=Tair(Itime)
c  MARCH ALONG THE SURFACE TO FIND SOIL-ATMOSPHERE BOUNDARY NODES
c  AND THEIR LATERAL COORDINATES
c

      
c  MARCH ALONG THE SURFACE TO FIND SOIL-ATMOSPHERE BOUNDARY NODES
c  AND THEIR LATERAL COORDINATES
c
      Do 95 i=1,ic
        PSh=FSh(i,xBSTEM,SHADE,xS,ic)
        n=kS(i)
        k=codeT(n)
        Varbt(kS(i),4)=0.0
        VarBT(kS(i),3)=0.0
        VarBT(kS(i),2)=0.0

        If (k.eq.-4) then
C     if node is exposed         
          If((PSh.gt.0).or.(NShoot.eq.0)) then
* If evaporation
* calc heat flux due to latent heat of evaporation 
* given evaporation at it's potential rate
* units are Joules cm-2 d-1 C-1
* RNS*3600*24/10000 is total energy going to soil in J day-1 cm-2
C           FELWR=max(0.0,(tmpr(i)-TairN)*40.6)
c           If (VarBW(i,2).ge.0.0.and.RNS.gt.0.0) then
                VarBT(kS(i),4)=RNS*3600.0*24.0/10000.0

*  calculate convected heat  if soil is warmer than air
* FELWR should be zero when TAirN is warmer than soil
              if(tmpr(kS(i)).GT.TAirN) then
                VarBT(kS(i),2)=+FCSH+FELWR
                VarBT(kS(i),3)=VarBT(kS(i),2)*TAirN
              endif
C See if sensible heat transport works here in canopy              
              TFac=((0.058+1.7e-4*TAirN)+
     !            (0.052*EXP(.058*TAirN)))*0.004184*3600.0*24.0
               VarBT(kS(i),2)=VarBT(kS(i),2)+TFac
               VarBT(kS(i),3)=VarBT(kS(i),3)+TFac*TairN

         else
                    
*  the amount of heat added to the soil
*     is determined by air temperature.
*  units are millcal cm-2 sec-1 C-1
               VarBT(kS(i),2)=((0.058+1.7e-4*TAirN)+
     !                  (0.052*EXP(.058*TAirN)))
* the first term is thermal conductivity of the air. 
* the second term is the conductivity of water vapor
* Change units to Joules cm-2 d-1

               VarBT(kS(i),2)=VarBT(kS(i),2)*0.004184*3600.0*24.0
               VarBT(kS(i),3)=VarBT(kS(i),2)*TairN

Cdt and a part for radiative transfer               
               if(tmpr(kS(i)).GT.TAirN) then
                 VarBT(kS(i),2)=VarBT(kS(i),2)+(FCSH+FELWR)
                 VarBT(kS(i),3)=VarBT(kS(i),3)+(FCSH+FELWR)*TAirN
                  endif
c assume 5% of radiation reaches soil surface through canopy              
              fact=1.0-BEERS(ITIME)
              if (fact.LE.0) fact=0
              if (fact.GT.1.0) fact=1.0
              VarBT(kS(i),4)=fact*RNS*3600.0*24.0/10000.0
              
*  units are now J cm-2 d-1
* end if for crop canopy
          Endif
 
          VarBT(kS(i),1)=TAirN
        Endif
95    Continue
C ...............End of heat balance
c................... Gas movement
      Do i=1,NumBP
      n=KXB(i)
      k=CodeG(n)
      If(K.eq.4.or.K.eq.-4) then                      !-4: soil atmospheric boundary with gas flux instead of Gas content at the nodes
        do jjj=1,NumG
          VarBG(i,jjj,2)=GasTransfCoeff(jjj)                !GasTransfCoeff(jjj) : is the conductance of surface air layer to gas flow or rate constant of the gas exchange between the soil and the atmosphere [cm/day]
          VarBG(i,jjj,3)=GasTransfCoeff(jjj) 
     !       *GAIR(jjj)/ugGasCm3air_to_ppm(jjj)             ! convert the initial concentration of CO2 in [ppm]  to [ug co2/cm3 air] !GAIR: is the atmospheric CO2 concentration [ppm], see conversion details in grid_bnd
          VarBG(i,jjj,1)= GAIR(jjj)/ugGasCm3air_to_ppm(jjj) ! convert the initial concentration of CO2 in [ppm]  to [ug co2/cm3 air] !This is if BC=4, is gas content instead of flux, 
        Enddo
      Endif
      if (K.eq.1.or.K.eq.3.or.k.eq.6) then
          do jjj=1,NumG
              VarBG(i,jjj,2)=GasTransfCoeff(jjj)  
              VarBG(i,jjj,1)= GAIR(jjj)/ugGasCm3air_to_ppm(jjj)          ! convert the atm CO2 in [ppm]  to [ug co2/cm3 air]       
          end do
      end if
      Enddo
     
 
c................... End gas movement      
c................... This is the end of hourly calculations
      tNext(ModNum)=tNext(ModNum)+PERIOD
      
      do i=1,NumBPD
cccz extract useful variables for surface physical processes
cccz water input from air-mulch interface
         Varbw_Air(i,1)=VarBW(i,1)
         Varbw_Air(i,2)=VarBW(i,2)
         Varbw_Air(i,3)=VarBW(i,3)
cccz water input from mulch soil interface (wait to be adjusted by surface runoff)
         Varbw_Mulch(i,1)=VarBW(i,1)
         Varbw_Mulch(i,2)=VarBW(i,2)
         Varbw_Mulch(i,3)=VarBW(i,3) 
cccz renormalize heat fluxes at air-mulch interface
         Varbt_Air(i,1)=Varbt(i,1)
         Varbt_Air(i,2)=Varbt(i,2)
         Varbt_Air(i,3)=Varbt(i,3)
         Varbt_Air(i,4)=Varbt(i,4)
cccz renormalize heat fluxes at mulch-soil interface (wait to be adjusted by surface runoff)
         Varbt_Mulch(i,1)=Varbt(i,1)
         Varbt_Mulch(i,2)=Varbt(i,2)
         Varbt_Mulch(i,3)=Varbt(i,3)
         Varbt_Mulch(i,4)=Varbt(i,4)
         do jjj=1,NumG
cccz renormalize gas fluxes at air-mulch interface
          Varbg_Air(i,jjj,1)=Varbg(i,jjj,1)
          Varbg_Air(i,jjj,2)=Varbg(i,jjj,2)
          Varbg_Air(i,jjj,3)=Varbg(i,jjj,3)
cccz renormalize gas fluxes at mulch-soil interface (wait to be adjusted by surface runoff)
          Varbg_Mulch(i,jjj,1)=Varbg(i,jjj,1)
          Varbg_Mulch(i,jjj,2)=Varbg(i,jjj,2)
          Varbg_Mulch(i,jjj,3)=Varbg(i,jjj,3)
         enddo
      enddo      
      
      
      EndIf                   ! hourly calcs defined by tnext
 1001  Continue
c................... Furrow irrigation
c dt 4-22-2011 this needs to be put into its own block since we changed the tAtm variable
c  this won't work as written
      If(IfUR.eq.1) then
        Do i=1,NumFP
        k=NumF(i)
        QF=QF-Q(k)*Step
        Enddo
        If(QF.lt.0.0) then
          QF=0.
          IfUR=0
          Do i=1,NumFP
           k=NumF(i)
           CodeW(k)=-4
          Enddo
          tAtm=tAtm-PERIOD
          goto 1000  !hourly code would have been skipped over initially, so have to go back and do the hourly et calcs
        Endif
      Endif
c................... End of the furrow irrigation
      RETURN
10    call errmes(im,il)
       write(*,*) "Error in Weather File"
      END
C
      Subroutine SORT03(N,ARR1,iARR2,iARR3)
      Dimension ARR1(N),iARR2(N),iARR3(N)
      Do j=2,N
         A1=ARR1(j)
         iA2=iARR2(j)
         iA3=iARR3(j)
         Do i=j-1,1,-1
           if(ARR1(i).LE.A1) goto 10
           ARR1(i+1)=ARR1(i)
           iARR2(i+1)=iARR2(i)
           iARR3(i+1)=iARR3(i)
         Enddo
         i=0
10       ARR1(i+1)=A1
         iARR2(i+1)=iA2
         iARR3(i+1)=iA3
      Enddo
      Return
      End
C   Shaded-unshaded
      Function FSh(i, xStem, Shade, u, N)
      Dimension u(N)
      Fsh=1.
      If(Shade.le.0.) Return
      If (i.eq.1) then
        xm=u(i)
        xp=(u(1)+u(2))/2.
      Elseif(i.eq.N) then
        xm=(u(N-1)+u(N))/2.
        xp=u(N)
      Else
        xm=(u(i-1)+u(i))/2.
        xp=(u(i)+u(i+1))/2.
      Endif
      xL=xSTEM-Shade
      xR=xSTEM+Shade
      If(xm.le.xL) then
        If(xp.lt.xL) then
          FSh=1.
        Elseif(xp.le.xR) then
          FSh=(xL-xm)/(xp-xm)
        Else
          FSh=1.-(xR-xL)/(xp-xm)
        Endif
      Else
        If(xm.ge.xR) then
          FSh=1.
        Elseif(xp.le.xR) then
          FSh=0.
        Else
          FSh=(xp-xR)/(xp-xm)
        Endif
      Endif
      Return
      End
