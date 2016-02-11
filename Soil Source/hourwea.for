cdt this version had code added to calculate rain infiltration as 
c constant head infiltration for a period. I have to finish it to allow final 
c flux to just finish up the water.
CDT  9/10/2014 Added CO2 as a weather variable and took it out of the 
CDT   initials file. 

C**  MSW1-  switch to indicate if daily wet bulb temperatures are     **
C**         available (=1 if YES).                                    **
C**                                                                   **
C**  MSW2-  switch to indicate if hourly wind is available (=1 if YES).**
C**  MSW3-  Not used in hourly data (Daily rain intensities)          **
C** cdt don't think these two variables (above) are needed  for hourly data
C**  MSW4-  Chemical concentrations in rain water                     **
C**  MSW5-  1 for flood irrigation                                    **
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
C**                                                                   **
C**  WIND-  if windspeed is unavailable read in average value for     **
C**         site. windspeed at 2 meters (km hr-1).                    **
C**                                                                   **
C** CO2   -  If daily values are not available read average           **
C**                                                                   **
C**                                                                   **
C**  (File 1, Line 3 to EOF, free format)                             **
C**      (RH, WIND, TWET and TDRY CO2 are optional)                  **
C**                                                                   **
C**  JDAY-  day of year- Julian Date.  Only used for management       **
C**         of the weather file. Program converts calendar date       **
C**         to a day of year relative to 3/1/1900                     **
C**  DATE-  calendar date in quotes "01/01/2001"                      **
C**                                                                   **
C**  HOUR-  Hour of Day                                               **                                                                **
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
C**  TWET-  wet bulb temperature (oC).                                **
C**                                                                   **
C**  TDRY- dry bulb temperature corresponding to TWET (oC).           **
C**                                                                   **
C**  Hourly values for radiation will be converted to watts m-2 as 
C**   soon as they are read in. 
C**                                                                   **

c inputs hourly data
      Subroutine SetSurfaceH()
      Include 'Public.ins'
      Include 'Puplant.ins'
      Include 'Puweath.ins'
      Parameter (PERIOD =1./24.)
      integer jday,m,DayOfYear,CurYear,Modnum, ThisYear
      double precision St,t
      real Interval
      character*10 date
      integer HOUR
      Common /weather/ il,im,HRAIN(24),HSR(24),HTEMP(24),HTEMPY(24), 
     &     HWIND(24), Rel_Humid(24),AVP(24),isol,Date1,ModNum,
     &     Interval
      
      Dimension CLIMAT(20),SDERP(9),SINALT(24),SINAZI(24),HRANG(24),
     &           DIFWAT(24),SARANG(24),DIFINT(24),ROWINC(24),DIRINT(24),
     &           SHADOW(24),SOLALT(24),SOLAZI(24),
     &           HTM(24)
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
C
C  First and last days
C
      JDFRST=idint(t+St)
      JDLAST=idint(tFin+St)
c
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum)=time
        Interval=1/24.

c
      im=160
      il=0
C      Open (5,file='Weatherhrly.dat',status='old',ERR=10)
      Open (5,file=ClimateFile,status='old',ERR=10)
cdt
c       Do i=1,NumNP
c        PcodeW(i)=CodeW(i)
c         Enddo
C
C  Read descriptors and conversion constants for weather data
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
       ISOL=NumSol*Movers(2) ! Movers(2) is the solute mover
       Read (5,*,ERR=10) BSOLAR,BTEMP,ATEMP,ERAIN,BWIND,BIR
cdt 10/25/2007 changed from 1-MSW4
      ISOL=NumSol*(MSW4)*Movers(2)
cdt  9/10/2014 added CO2 as a weather variable. Made this line consistent
Cdt   with the same line from the daily weather file 
CDT    TODO this needs to be tested fully
CDT  1/2016 fix for dew point and co2
CDT  this line determines what is read from the last line in the climate file - this last line contains averages to be
CDT  used when daily/hourly values are not available. 
CDT at the minimum (for hourly) ncd will be three - wind, irav, CO2, because only these three have average values
CDT to be used in the program.  

CDT columns should be ordered as in the msw indices in the climate file
CDT  hourly bulb, hourly wind, hourl rain intensity, hourly concentration, hourly furrow, hourly
CDT  relative humidity, hourly CO2. For now we will only consider one solute in the 
CDT  rainfall, that is Nitrogen. We may have more than one column of concentrations if there is more than one solute
CDT  but this is not fully implemented yet 
CDT if MSW6 is 0 then there is no column for RH and it is calculated from minimum temperature - no average is used
CDT if MSW is one, then RH is not used even if present
CDT daily concentration has been modified to use only nitrogen in rainfall (no gasses).

      NCD=3-MSW2+ISOL-MSW3-MSW7
C      NCD=1-MSW2+ISOL+(NumG+1)*Movers(4)
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
      IF(MSW3.eq.0) IRAV= CLIMAT(2)
C     IF(MSW4.eq.0.AND.NumSol.gt.0) then
Cdt 02/03/2009 changed from 0-MSW4      
      IF(MSW4.eq.1.AND.NumSol.gt.0) then
        Do i=1,NumSol
          CPREC(i)=CLIMAT(2-MSW2-MSW3+i) ! should be 1 or 2 if no irav or winda and one or two solutes
        Enddo
      Endif
CDT 1/2016 no gasses for now
!      If(NumG.ne.0) then
!        PG=CLIMAT(3-MSW2+ISOL)
!        Do i=1,NumG
!          GAIR(i)=CLIMAT(3-MSW2-MSW3+ISOL+i)
!        Enddo
!      Endif
      IF(MSW7.eq.0) then
         CO2=CLIMAT(NCD)  ! last one is always co2
       Endif

      close(5)
C
C Total number of weather data
C     minimum of three columns for radiation, temperature and rainfall (should always be first three) for hourly
C
      NCD=3+2*MSW1+MSW2+ISOL*MSW4+MSW6+MSW7
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
      Read (5,*,ERR=10) iDum, date, iDum
      MDAY=julday(date)
C       MDAY = 38065
cdt changed from LE to LT since jday was always gt. time
      If (MDAY.LT.JDFRST) GO TO 1111
      backspace(5)

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


c..................... Input daily data
   11 If((linput.eq.1).OR.idint(t+St).eq.(JDAY+1)) then
      
      tAtm=t          !insure we call the weather code the next step
      
c       Do i=1,NumNP
c        CodeW(i)=PcodeW(i)
c         Enddo

c
c ZEROING DAILY ARRAYS
c
        Do i=1,24
          RINT(i)=0.
          RADINT(i)=0.
          WATTSM(i)=0.
          il=il-1
        Enddo
c
c Read meteodata
c
        Do m=1,24
          il=il+1
C          Read(5,*,ERR=10) JDAY,HOUR, DATE1, (CLIMAT(i),i=1,NCD)
          Read(5,*,ERR=10) JDAY, date, HOUR,  (CLIMAT(i),i=1,NCD)
           JDAY=julday(date)
c    SInce some routines need the day of the calendar year we need a variable to hold this
C    since the julian day is referenced to a time longer in the past
           ThisYear=CurYear(date)
           write (date,'("01/01/",i4.4)') ThisYear  
           DayOfYear=JDAY-julday(date)
        
           HSR(M)=max(0.0,climat(1))*BSOLAR/3600   ! convert to watts m-2 hard 
                                          !coded for hourly now

           if (lInput.ne.1) HTEMPY(M)=HTEMP(M)  ! save to yesterday's 
                                                !temperature
           HTEMP(M)=(climat(2)-ATEMP)/BTEMP     ! convert to celcius
           HRAIN(M) = CLIMAT(3)*erain         ! convert to cm rain in 
                                              ! this hour
           If (MSW2.gt.0) HWIND(M)=CLIMAT(4)
           If(MSW6.gt.0) then 
             Rel_Humid(m)=Min(Climat(3+2*MSW1+MSW2+ISOL+MSW6)
     #              /100.,0.98)
            endif
         enddo

       If(MSW7.gt.0) then
          CO2=Climat(4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6+MSW7)
       EndIf


C
C      ADJUST UNITS FOR HOURLY DATA
c      HSR(M) is converted to an hours worth of Watts m-2 (Joules m-2 h-1) 
C
202          RI = 0.0 
          SUM24 = 0.0
          TMAX = 0.0
          TMIN = 1000.0
          DO 1015 M = 1,24
c RI should be total J m-2 -- work in Watts * time = total energy
             RI = RI + amax1(0.0,HSR(M))*3600
             SUM24 = SUM24 + HRAIN(M)
             IF(TMAX .LT. HTEMP(M)) THEN
               TMAX = HTEMP(M)
             ENDIF
             IF(TMIN .GT. HTEMP(M)) THEN
               TMIN = HTEMP(M)
             ENDIF
             IF(MSW2 .GT. 0) THEN
               WIND = HWIND(M) * BWIND
             Else 
               WIND=WINDA*BWIND
             ENDIF
1015      CONTINUE
          RAIN = SUM24
        If (MSW1.GT.0) then
          TWET = (CLIMAT(5+MSW2) - ATEMP)/BTEMP
          TDRY = (CLIMAT(6+MSW2) - ATEMP)/BTEMP
        Endif
        If (MSW4.GT.0.AND.NumSol.ne.0) then
          Do i=1,NumSol
            CPREC(i)=CLIMAT(4+MSW2+2*MSW1+MSW3+i)
          Enddo
        Endif



c......................... Radiation submodel
C
C  CALCULATE SOLAR DECLINATION
C
        XLAT = LATUDE*DEGRAD
        DEC = SDERP(1)
        Do I = 2,5
          N = I - 1
          J = I + 4
          D11 = N*0.01721*JDAY
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
        WATACT = RI*4.363E-4/DAYLNG     ! Watts per m2
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
cdbg check the units here
     &      WATTSM(I) = HSR(I)
        Enddo
        If ((DDIf.GT.0.0).OR.(DDIf.LT.0.0)) then
           WATTSM(IDAWN) = HSR(IDAWN)
           WATTSM(IDUSK) = HSR(IDUSK)
        ENDIf

c..................... Temperature and vapour pressure submodel


C
C  ASSIGN AIR TEMPERATURE AT EACH TIME
C   
C FIND THE TIME AFTER DAWN IN HOURS WHEN MAXIMUM 
C  TEMPERATURE IS REACHED
C    
        tMax=-1e6
        Do 70, m = 1,IPERD
           TAIR(M) = (HTEMP(M)-ATEMP)/BTEMP
           if(m.eq.IDUSK) TDUSK=TAIR(m)
           if (TAIR(M).GT.TMAX) then
             TMAX=TAIR(m)
             TMAXHR=m
            End If
 70   CONTINUE
 

C
C  If DAILY WET BULB TEMPERATURES ARE NOT AVAILABLE, ASSUME
C  MINIMUM TEMPERATURE = DEW POINT TEMPERATURE AND CALCULATE
C  WATER VAPOR PRESSURE.  ASSUME A VALUE FOR THE PSYCHROMETRIC
C  "CONSTANT"
C
cdt added msw6 to if statement here

        If (MSW1.NE.1.and.MSW6.eq.0) then
          Do i=1,iperd
            AVP(i) = 0.61*EXP((17.27*TMIN)/(TMIN + 237.3))
           enddo
          GAMMA = 0.0645
        Else
cdt
           if(MSW1.eq.1) then      
C
C  SINCE DAILY WET BULB TEMPERATURES ARE AVAILABLE CALCULATE
C  SATURATION VAPOR PRESSURE AT THE WET BULB TEMPERATURES
C    http://www.pmel.org/HandBook/HBpage21.htm see this site for more info
C    on using Dewpoint to calculate RH
C
          SVPW = 0.61*EXP((17.27*TWET)/(TWET + 237.3))
C
C  CALCULATE THE HUMIDITY RATIO, D31, LATENT HEAT OF EVAPORATION,
C  D32 AND THE PSYCHROMETRIC "CONSTANT"
C
          D31 = 0.622*(SVPW/(101.3 - SVPW))
          D32 = 2500.8 - (2.37*TWET)
          GAMMA = 0.62*(1.006 + (1.846*D31))
     &   /((0.622 + D31)*(0.622 + D31)*D32)*101.3
C
C  CALCULATE ACTUAL WATER VAPOR PRESSURE
C
          Do i=1,iperd
             AVP(i) = SVPW - (GAMMA*(TDRY - TWET))
          enddo
        Endif

CDT        
            if(MSW6.eq.1) then
             Do i=1,iperd                   
                TMean=(TMax+TMin)/2.0
                SVPW= 0.61*EXP((17.27*HTEMP(i))/(HTEMP(i) + 237.3))
                D31 = 0.622*(SVPW/(101.3 - SVPW))
                D32 = 2500.8 - (2.37*HTEMP(i))
                GAMMA = 0.62*(1.006 + (1.846*D31))
     &          /((0.622 + D31)*(0.622 + D31)*D32)*101.3
                AVP(i)= Rel_Humid(i)*SVPW
              enddo
          endif
       
        Endif
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
          VPD(I) = SVPA - AVP(I)
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
            SINAZI(I) = -COS(DEC)*SIN(HRANG(I))/COS(SOLALT(I))
 30       Continue
          D18 = 0.0
          Do 60, I = 1,IHPERD
            SOLAZI(I) = PI + ASIN(SINAZI(I))
            If (.NOT.((D18.GT.0.0).OR.(D18.LT.0.0))) then
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
C  FIRST DIVIDE ROWSPACING INTO 20 EQUAL PARTS AND FIND THE
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
                BEERS = (1 - EXP(-ELCAI*CEC))
              Else
                BEERS=1.0
              Endif

C  CALCULATE PROPORTION OF PHOTOSYNTHETICALLY ACTIVE RADIATION
C   (PAR) INTERCEPTED BY THE CROP CANOPY ASSUMING THAT
C    PAR = 0.43 DIRECT + 0.57 DIFFUSE RADIATION. ***

      PARINT(I) = WATTSM(I)*((DIFFIN*0.57*DIFWAT(I))
     &     + (DIRINT(I)*0.43*(1.0 - DIFWAT(I))))
     &     /PAR(I)*BEERS
C
C  CALCULATE PROPORTION OF TOTAL RADIATION INTERCEPTED BY THE
C  CROP CANOPY
C
              RADINT(I) = ((DIfFIN*DIfWAT(I)) + (DIRINT(I)
     &         *(1.0 - DIfWAT(I))))*BEERS
            Else If (WATTSM(I).GT.0.0) then
              BEERS = (1 - EXP(-LCAI*CEC))
              PARINT(I) = COVER*BEERS
              RADINT(I) = COVER*BEERS
            Else
               PARINT(I) = 0.0
               RADINT(I) = 0.0
            Endif
 50       Continue
       

C
c..................End of radiation interception submodel


        Endif  !nshoot >0
C
c..................Daily calculations for the precipitation
C
          Do i=1,24
            RINT(i)=HRain(i)*24.0  ! scale hourly rainfall to a daily rate
          Enddo

        If(Rain.lt.0.) then
          QF=QF-Rain
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
55    If(Abs(Time-tNext(ModNum)).lt.0.01*Step.or.lInput.eq.1) then
       ITIME=Idint(dMOD(t+St,1.D0)/PERIOD+1)
       if(itime.eq.10) then
        iii=1;
        endif
        
       Do i=1, NumBP
         do j=1, 3
           VarBW(i,j)=0.0
         Enddo
        Enddo

c       Do i=1,NumNP
c        CodeW(i)=PCodeW(i)
c         Enddo

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
C
        Do 90 i=1,ic
          PSh=FSh(i,xBSTEM,SHADE,xS,ic)
c
c  IF THE NODE IS EXPOSED THEN
c
          If(PSh.gt.0.) then
C
C  CALCULATE THE ALBEDO OF THE EXPOSED SOIL CELLS
C
cdt changed 0.5 to 0.1
              LAMDAS = 0.25 - 0.05*ThNew(iS(i))
C
C  CALCULATE NET RADIATION ON THE EXPOSED SOIL CELLS
C
              RNS = ((1.0 - LAMDAS)*D11) - RNLU
              If(RNS.LE.0.0) RNS = 0.0
C
C  CALCULATE POTENTIAL EVAPORATION RATE FROM EXPOSED SOIL CELLS
C
              D12 = max(1.0,2.0 - (2.0*COVER))
              If (D12.GE.1.0) D12 = 1.0
              ESO = ((DEL(ITIME)/GAMMA*RNS*3600.0/(2500.8
     &     - (2.3668*TAIR(ITIME))))
     &     + (VPD(ITIME)*109.375*(1.0 + (0.149*WIND*D12))))
     &      /((DEL(ITIME)/GAMMA) + 1.0)
C
C   PROPORTIONAL TO THE AREA EXPOSED
C
            
              VarBW(kS(i),2)=PSh*24.*ESO/10000.
c
c  ELSE IF THE NODE IS COVERED
c
          Else
            VarBW(kS(i),2)=0.
          Endif
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
          WATPL = WATTSM(ITIME) !*RADINT(ITIME)/COVER /dt don't need when used in 2dsoil
        Endif
C
C   CALCULATE NET RADIATION ON THE CROP
C
        RNC = ((1.0 - LAMDAC)*WATPL) - RNLU
        If (RNC.LE.0.0) RNC = 0.0
C
C   CALCULATE POTENTIAL TRANSPIRATION RATE FOR CROP ALLOWING FOR
C   INCOMPLETE GROUND COVER
C
        EPO = ((DEL(ITIME)/GAMMA*RNC*3600.0/(2500.8
     &   - (2.3668*TAIR(ITIME))))
     &   + (VPD(ITIME)*109.375*(ROUGH +(0.149*WINDL))))
     &   /((DEL(ITIME)/GAMMA) + 1.0)
C
C   CODE ADDED TO PREVENT DIVISION BY o
C
        If (COVER.le.0.0.OR.EPO.LE.0.0) EPO=1.0E-4
C
C     *** CALCULATE POTENTIAL TRANSPIRATION RATE FOR AN AVERAGE PLANT **
C
cdt I took out etcorr, the ET correction factor - units are per plant
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
         Q(n)=-Width(i)*VarBW(i,3)
         If (Q(n).gt.0.0) then
          CodeW(n)=-4
          endif
       Endif
       Enddo
c................... End of the precipitation and irrigation
C ............... Heat Balance is calculated here
C
C
CYAP 
C FCSH is sensible heat flux when the ground is hotter than the air
       if(rowsp.le.0) then
         FCSH=4.0E-3+1.39E-3*Wind
        else
         FCSH=4.0E-3+1.39E-3*Wind*max(1.0,(2.-(2.*height/rowsp)))
       endif
       
c      FCSH=0.02
cdt 
* FELWR is reradiation when soil is warmer than the air.
      FCSH=FCSH*4.1856*60.0*24.0
      FELWR=6.8E-03+9.0E-05*Tair(Itime)
      FELWR=FELWR*4.1856*60.0*24.0
* units are Joules cm-2 d-1 C-1
CYAPEND
C units are changed from cal cm-2 min-1 C-1 -->J cm-2 d-1 C-1
      cWat=4.1856
      TairN=Tair(Itime)
      Do i=1,NumBP
        n=kxb(i)
        k=codeT(n)
        Varbt(i,4)=0.0
        VarBT(i,3)=0.0
        VarBT(i,2)=0.0
        If (k.eq.-4) then
          If(abs(x(n)-xBSTEM).gt.SHADE.or.NShoot.eq.0) then
* If evaporation
* calc heat flux due to latent heat of evaporation 
* given evaporation at it's potential rate
* units are J day-1
* RNS*3600*24/10000 is total energy going to soil in J day-1 cm-2
*  40.6 for mississippi 
c           FELWR=max(0.0,(tmpr(i)-TairN)*FELWR)
c           If (VarBW(i,2).ge.0.0.and.RNS.gt.0.0) then
                VarBT(i,4)=RNS*3600.0*24.0/10000.0

*  calculate convected heat  if soil is warmer than air
     
              VarBT(i,3)=(FCSH+FELWR)*TAirN
              VarBT(i,2)=FCSH+FELWR 
         else

*  the amount of heat added to the soil
*     is determined by air temperature.
*  units are millcal cm-2 min-1 C-1
               VarBT(i,2)=((0.058+1.7e-4*TAirN)+
     !                  (0.052*EXP(.058*TAirN)))
* the first term is thermal conductivity of the air. 
* the second term is the conductivity of water vapor
* Change units to Joules cm-2 d-1
cdt change this for a test
               VarBT(i,2)=VarBT(i,2)*4.1856*60.0*24.0/1000
cdt /1000.0
               VarBT(i,3)=VarBT(i,2)*TairN
*  units are now J cm-2 d-1
*  heat into soil with rain
* end if for crop canopy
        endif


* add heat input with rain
c       If(VarBW(i,3).lt.0.0) then
c         VarBT(i,2)=VarBT(i,2)+varBW(i,3)*cWat
c         VarBT(i,3)=VarBT(i,3)+VarBW(i,3)*cWat*TAirN
         
* should this be + varbw*width?
c               Endif

       VarBT(i,1)=TAirN
        Endif
       Enddo
C ...............End of heat balance
c................... Gas movement
      Do i=1,NumBP
      n=KXB(i)
      k=CodeW(n)
      If(K.eq.4.or.K.eq.-4) then
        do jjj=1,NumG
          VarBG(i,jjj,2)=PG*Width(i)
          VarBG(i,jjj,3)=PG*GAIR(jjj)*Width(i)
        Enddo
      Endif
      Enddo
c................... End gas movement   
c................... Furrow irrigation
c
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
        Endif
      Endif
c................... End of the furrow irrigation
   
c................... This is the end of hourly calculations

      tNext(ModNum)=Time+period
      Endif      ! hourly loop
c      

      
 

      RETURN
10    call errmes(im,il)
      END
C
     