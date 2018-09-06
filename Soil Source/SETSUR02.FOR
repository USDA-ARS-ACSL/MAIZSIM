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
C** cdt don't think these two variables (above) are needed  for hourly data
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
C**  IRAV-  avg Rainfall intensity for the day, (mm/day, hr?)              **
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
      Include 'Public.ins'
      Include 'Puplant.ins'
      Include 'Puweath.ins'
      Parameter (PERIOD =1./24.)
      integer jday, DayofYear,CurYear, ModNum,ThisYear
      Common /weather/ il,im,isol,ModNum
      
      double precision St,t
      character*10 date
      Dimension CLIMAT(20),SDERP(9),SINALT(24),SINAZI(24),HRANG(24),
     &           DIFWAT(24),SARANG(24),DIFINT(24),ROWINC(24),DIRINT(24),
     &           SHADOW(24),SOLALT(24),SOLAZI(24)
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
      ISOL=NumSol*Movers(2)
c      the four basic variables are Daily wind, rain intensity daily concentration and co2     
      NCD=4-MSW2-MSW3+ISOL+(NumG+1)*Movers(4)-MSW7
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
      IF(MSW3.eq.0) IRAV =CLIMAT(2-MSW2-MSW3)
      IF(MSW4.eq.1.AND.NumSol.gt.0) then
        Do i=1,NumSol
          CPREC(i)=CLIMAT(2-MSW2-MSW3+i)
        Enddo
      Endif
      If(NumG.ne.0) then
        PG=CLIMAT(3-MSW2-MSW3+ISOL)
        Do i=1,NumG
          GAIR(i)=CLIMAT(3-MSW2-MSW3+ISOL+i)
        Enddo
      Endif
      IF(MSW7.eq.0) then
         CO2=CLIMAT(NCD) ! last one is always co2
       Endif
      close(5)      
C
C Total number of weather data
C
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
      Read (5,*,ERR=10) iDum, date
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
        Read(5,*,ERR=10) JDum,date,(CLIMAT(i),i=1,NCD)
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
        TMIN=max(TMIN,3.0)        
        RAIN = CLIMAT(4)*ERAIN

        
c        If (rain.gt.0) then
c         do ii=1,NumBp
c           n=kxb(ii)
c           If (abs(codew(ii)).eq.4) then
c             hnew(n)=-1.0e3/(1+0.04*1.0e3)
c            EndIf
c          enddo
c         Endif
          




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
          IR=CLIMAT(4+MSW2+2*MSW1)*BIR
        Else
          IR=IRAV*BIR
        Endif
        If (MSW4.GT.0.AND.NumSol.ne.0) then
          Do i=1,NumSol
            CPREC(i)=CLIMAT(4+MSW2+2*MSW1+MSW3+i)
          Enddo
        Endif
        AutoIrrigAmt=IR   ! used in AutoIrrigation if carried out.
CDT

      If(MSW6.gt.0) then
c        Rel_Humid=Climat(4+2*MSW1+MSW2+MSW3+ISOL+MSW6)
         Rel_Humid=min(Climat(4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6)
     &     /100.0,0.98)
          Rel_Humid=max(0.0,Rel_Humid)
c         Rel_Humid=CLIMAT(5)
       endif
      If(MSW7.gt.0) then
          CO2=Climat(4+2*MSW1+MSW2+MSW3+ISOL*MSW4+MSW6+MSW7)
       EndIf
Cdt
c......................... Radiation submodel
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
c..................... Temperature and vapour pressure submodel
C
C  READ NEXT DAY'S DATA FROM WEATHER FILE ***
C
        il=il+1
        Read (5,*,ERR=10) MDAY, date,(CLIMAT(i),i=1,NCD)
c       Minimum tomorrow temperature
        MDAY=Julday(date)
        If (MDAY.LT.JDLAST) then
          TMINT = (CLIMAT(3) - ATEMP)/BTEMP
        Else
          TMINT = TMIN
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
          GAMMA = 0.0645
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
          GAMMA = 0.62*(1.006 + (1.846*D31))
     &   /((0.622 + D31)*(0.622 + D31)*D32)*101.3
C
C  CALCULATE ACTUAL WATER VAPOR PRESSURE
C
          AVP = SVPW - (GAMMA*(TDRY - TWET))
        Endif

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
                GAMMA = 0.62*(1.006 + (1.846*D31))
     &          /((0.622 + D31)*(0.622 + D31)*D32)*101.3
                AVP= Rel_Humid*SVPW
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
        Endif
C
c..................End of radiation interception submodel
C
C
c..................Daily calculations for the precipitation
C

        If(Rain.ge.0.) then
          DtRain=RAIN/IR
          Do i=24,1,-1
            If(DtRain.ge.PERIOD) then
              RINT(i)=IR
              dtRain=dtRain-PERIOD
            Else
cdt changed this from dtrain/ir 11/14/96
              RINT(i)=DtRain*IR/period
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
cdt        If (Cover.lt.0.97) then
C
          D11 = WATTSM(ITIME)
cdt      in glycim, radint is not used for soil evaportion
cdt     &   *(1.0 - RADINT(ITIME))*(1.-COVER)
          
cdt   original above was correct
cdt     &   *(1.0 - RADINT(ITIME))/(1.-COVER)
c
cdt        Else
cdt         D11=0.
cdt        Endif
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
C initialize boundary values to 0
         Do i=1, NumBP
         do j=1, 4
           if (j.lt.4) VarBW(i,j)=0.0
           VarBT(i,j)=0.0
         Enddo
        Enddo
        Do 90 i=1,ic
cDT 4-22-2011 assume 10% of et is minimum from a covered node

          PSh=max(0.0,FSh(i,xBSTEM,SHADE,xS,ic))
c
c  IF THE NODE IS EXPOSED THEN
c
          If(PSh.gt.0.) then
C
C  CALCULATE THE ALBEDO OF THE EXPOSED SOIL CELLS
C
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


C See eqn 5.11 of http://mepas.pnnl.gov/mepas/formulations/source_term/5_0/5_13/5_13.html
              D12 = max(1.0,2.0-(2.0*COVER))
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
cDT 4-22-2011 assume 10% of et is minimum from a covered node
c
c          Else
c            VarBW(kS(i),2)=0.
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
         if (Q(n).gt.0.0) CodeW(n)=-4
*
       Endif
       Enddo
c................... End of the precipitation and irrigation
C ............... Heat Balance is calculated here
C
C
CYAP
CDT - added FELWR based on changes from temperature paper
C FCSH is the constant for sensible heat transfer 
C FELWR is the coefficient for radiant heat lost when the 
C soil is warmer than the air.
      FCSH=(4.0E-3+1.39E-3*Wind*max(1.0,2.0-(2.0*COVER)))
      FCSH=FCSH*4.1856*60.0*24.0
      FELWR=6.8E-03+9.0E-05*Tair(Itime)
      FELWR=FELWR*4.1856*60.0*24.0
CYAPEND
C units are changed from cal cm-2 min-1 C-1 -->J cm-2 d-1 C-1
      cWat=4.1856
      TairN=Tair(Itime)
c    conversion factor 
c 8.64 = 3600 * 24 / 10,000 ; J m-2 s-1--> J cm-2 d-1
      Do i=1,NumBP
        n=kxb(i)
        k=codeT(n)
        VarBT(i,3)=0.0
        VarBT(i,2)=0.0
         VARBT(i,4)=0
        If (k.eq.-4) then
          If(abs(x(n)-xBSTEM).ge.SHADE.or.NShoot.eq.0) then
* If evaporation
* calc heat flux due to latent heat of evaporation 
* given evaporation at it's potential rate
* units are Joules cm-2 d-1 C-1
* RNS*3600*24/10000 is total energy going to soil in J day-1 cm-2
C           FELWR=max(0.0,(tmpr(i)-TairN)*40.6)
c           If (VarBW(i,2).ge.0.0.and.RNS.gt.0.0) then
                VarBT(i,4)=RNS*3600.0*24.0/10000.0

CDT            If (VarBW(i,3).ge.0.0) then
CDT               VarBT(i,3)=VarBW(i,3)*(2500.8-2.3668*TairN)
CDT            Endif
* calculate potential heat input due to precipitation
cdt heat input due to rain is done now for both shaded and bare together
* was else if
*             If (VarBw(i,3).lt.0.0) then
*                 VarBT(i,3)=VarBT(i,3)-VarBW(i,3)*CWat
*               Endif
*  calculate convected heat
cdt added FELWR from temperature paper results
cdt 11/2/2016 removed summation for varbt Now is consistent with hourly
                 VarBT(i,3)=(FCSH+FELWR)*TAirN  
                 VarBT(i,2)=FCSH+FELWR
                    
          else 
*  parameters beneath the canopy
*  the amount of heat added to the soil
*     is determined by air temperature.
*  units are millcal m d-1 C-1
              
               VarBT(i,2)=((0.058+1.7e-4*TAirN)+
     !                  (0.052*EXP(.058*TAirN)))

* Change units to Joules cm-2 d-1
cdt changed 3600 to 60 based on temperature paper
cdt 11/2/2016 check to see if we need the 100.
               VarBT(i,2)=VarBT(i,2)*4.1856/100.0*3600*24
* cdt for now use a constant 62.6 cal /sec /cm /C
* use 3600*4.184 as conversion
* units are J cm-1 d-1 oC (thermal cond of air)
*               VarBT(i,2)=2.05
               VarBT(i,3)=VarBT(i,2)*(TairN)
*  units are now J cm-2 d-1


* end if for crop canopy
          Endif
* Endif in or out of canopy
cdt modified statements after 'then' to add heat content of rain
*  heat into soil with rain 
             If(VarBW(i,3).lt.0.0) then          
                  VarBT(i,3)=VarBT(i,3)+VarBW(i,3)*cWat*TAirN
                  VarBT(i,2)=VarBT(i,2)+VarBW(i,3)*cWat
               Endif
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
c................... This is the end of hourly calculations
      tNext(ModNum)=tNext(ModNum)+PERIOD
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
