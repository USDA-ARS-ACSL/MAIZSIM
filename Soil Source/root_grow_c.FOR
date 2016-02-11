*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Subroutine Rootgrow ()
      Include 'public.ins'
      Include 'puplant.ins'
      Dimension iLoc(4)
      Integer e,e1,SW10,SW20,SW30,SW40
      Real*4  WeightNeighbor(4),LeftDist, RightDist, TopDist,BottomDist
      
                             
      t=sngl(Time)
      If (lInput.eq.0) goto 11
      HourlyCarboUsed=0          !hourly initializtion in crop model
C
C  Reading of the input files and initial calculations 
C
c      open(105, file='root.dat',recl=150)
c      write(105,'(10A12)')'time    ','elem   ' , 'AWR    ', 'PDWR     ', 
c     !  'pcrts    ' ,'pcrs    ','convr   ', 'rootwt    ', 'id',
c     !   'Matnum'
!        TotalRootWeight=0.0    !temporary
          
11    continue
C

C     ---------- Taken form R U T G R O   of  G L Y C I M
C
C  CALCULATE SOIL MECHANICAL RESISTANCE TO ROOT GROWTH IN EACH CELL
C
C  CALCULATE GROWTH BY NEW ROOTS IF ROOTS GREW AT POTENTIAL
C  RATE DURING PAST PERIOD. START WITH SOIL CELL WHERE CONDITIONS
C  ARE MOST FAVORABLE FOR GROWTH


C  No Plant = No Root Activity
      If(NShoot.eq.0.or.isGerminated.eq.0) Return
C

C calculate potential growth in root length (DT- May not need this loop - look into this)
      DO 30, e=1,NumEl
        e1=iFavRoot(e)
        PDRL = PDWR(e1)/RTWL ! potential rate of change in root length cm d-1
        DRL(e1) = (PDRL + PPDRL(e1))*Step/2.0  ! average here
        PPDRL(e1) = PDRL  ! value at previous calculation time

30    CONTINUE

C   Taken from ACTGRO of GLYCIM
C   THIS PART GROWS ROOTS IN THE MOST FAVORABLE CELLS UNTIL 
C   CARBON AVAILABLE IS USED
C
C   When the cells are small, the model does not always use all the carbon available
C   We can loop this until all carbon is used up or increase the growth rate. The latter
C   would require some iteration, it may be best to loop several times, the only problem is
C   that it would weight the first cells maybe too hightly
      PCRTS = 0.0
      SW40 = 0
      used=0.0
      DO 60 e = 1,NumEl
           e1=iFavRoot(e)
        IF (PDWR(e1).LE.0.0) THEN
          ADWR(e1) = 0.0
        ELSE
          PCRTS = PCRTS + (PDWR(e1)*CONVR)
        
          IF (PCRS.GE.PCRTS) THEN
            ADWR(e1) = PDWR(e1) 
            used=used+PDWR(e1)*CONVR
           ll=1
c           write(105,25) time, e1, ADWR(e1),PDWR(e1)*CONVR, PCRTS, 
c     !                        PCRS,Used,rtwt(e1), ll, matnumE(e1)
              ELSE IF ((SW40.GT.1).OR.(SW40.LT.1)) THEN
            SW40 = 1
            ADWR(e1) = ((PDWR(e1)*CONVR) - PCRTS + 
     !               PCRS)/CONVR
           ll=2
           used=used+ADWR(e1)*CONVR
c            write(105,25) time, e1, ADWR(e1),PDWR(e1)*CONVR, PCRTS, 
c     !                        PCRS,USED,rtwt(e1), ll, matnumE(e1)
          ELSE
            ADWR(e1) = 0.0
          END IF
        END IF
C
C  CALCULATE ROOT GROWTH RATE
        AWR(e1) = ADWR(e1)*Step
 60   CONTINUE
      fdif=pcrs-used
      if (fdif.gt.0) then
         iii=1;
        endif
C
CDT 03/13/2011 Added this line from GLYCIM. PCRS returns the actual amount of carbon used 
C If the roots don't use all the carbon allocated then pcrts will be less than pcrs
C if the roots use all the carbon then the code will overshoot pcrts slightly (the growth in the element 
C  is adjusted but not the summation of pcrts.
      PCRS=AMIN1(PCRTS,PCRS)
c    cdt HourlyCarboUsed is the summation of carbohydrate actually used to  grow roots
c       summed between calls to the plant model. It should be initialized to 0 hourly
c     in the plant model 
      HourlyCarboUsed=HourlyCarboUsed+PCRS*step
      
      
      Do e=1,NumEl
        ADRL(e) = 0.0
      Enddo
      
	DO 80	 e = 1,NumEl 
		 iCell = iFavRoot(e)
C   Note: GLYCIM calculates	RtMinWtPerUnitAra from RTWL,
C    RTMINW	= RTWL*SQRT(DEPTH*DEPTH	+ WIDTH*WIDTH)			
C Root weight must be more than	the	weight needed to cross the cell	diagonally for
C the	root to grow into	the	neiboring cell
	IF ((AWR(iCell).GT.0.0).OR.(AWR(iCell).LT.0.0)) THEN

C AD took out loop for triangular cells
C
C   MINIMUM	ROOT WEIGHT	TO CROSS A SOIL	CELL BOUNDARY for rectangular elements
C
cdt changed 2/17/2010 so that minimum length is defined by the distance
c   from the center of the cell to the boundary, (up, down, left, right)
c   In this manner, the weights are adjusted based on wether or not there 
c   is enough root lenght to go to the next cell
c   need four variables RootMinWtLeft, RootMinWtRight, RootMinWtTop
c   and RootMinWtBot
		
		RTMINW1	= RTWL*SQRT(HeightE(iCell)**2+ WidthE(iCell)**2)
		RTMINW = RtMinWtPerUnitArea	* Area(iCell)
		LeftDist=WidthE(iCell)/2.0
		RightDist=LeftDist
		TopDist=HeightE(iCell)/2.0
		BottomDist=TopDist
		
		RootMinWtRight =RTWL*RightDist
		RootMinWtLeft =RTWL*LeftDist
		RootMinWtTop  =RTWL*TopDist
		RootMinWtBot  =RTWL*BottomDist
		IF (RTWT(iCell).GT.RTMINW) THEN
		  iLeft	 = AdjCellNum(iCell,1)
		  iUpper = AdjCellNum(iCell,2)
		  iRight = AdjCellNum(iCell,3)
		  iLower = AdjCellNum(iCell,4)
		  If(RVR(iLeft).lt.RVR(iRight))	then
			D17	= Wl
			D18	= Wr
		  Else
			D17	= Wr
			D18	= Wl
		  Endif	
              
C AD took out code for triangular cells

C The	wider	the	source cell	is,	the	larger the weight to grow up or	down
C relative to growing left or right.	
C this (AWTot) is for testing how much root mass is transferred
          AWTot=0.0
		  WeightCell = RGCF(iCell) *(WidthE(iCell)+HeightE(iCell))/2

		  WeightLeft = RGCF(iLeft) * HeightE(iCell)	* D17
		  WeightRight= RGCF(iRight)* HeightE(iCell)	* D18
		  WeightUpper= Wa*RGCF(iUpper)*	WidthE(iCell)
		  WeightLower= Wb*RGCF(iLower)*	WidthE(iCell)
		  D19 =	   WeightCell 
     !			  +	   WeightLeft
     !			  +	   WeightRight
     !			  +	   WeightUpper
     !			  +	   WeightLower
		
		  IF (.NOT.((D19.GT.0.).OR.(D19.LT.0.))) GOTO 1212
		  D20 =	AWR(iCell)*	WeightCell /D19
		  AWTot=AWTot+D20
		  RTWT(iCell) =	RTWT(iCell)	+ D20
		  ADRL(iCell) =	ADRL(iCell)	+ D20/RTWL
		  
		  D21 =	AWR(iCell)*	WeightLeft /D19
		  AWTot=AWTot+D21
		  RTWT(iLeft) =	RTWT(iLeft)	+ D21
		  ADRL(iLeft) =	ADRL(iLeft)	+ D21/RTWL
		  
		  D22 =	AWR(iCell)*	WeightRight	/D19
		  RTWT(iRight) = RTWT(iRight) +	D22
		  AWTot=AWTot+D22
		  ADRL(iRight) = ADRL(iRight) +	D22/RTWL
		  
		  D23 =	AWR(iCell)*	WeightUpper	/D19
		  AWTot=AWTot+D23
		  RTWT(iUpper) = RTWT(iUpper) +	D23
		  ADRL(iUpper) = ADRL(iUpper) +	D23/RTWL
		  
		  D24 =	AWR(iCell)*	WeightLower	/D19
		  AWTot=AWTot+D24
		  RTWT(iLower) = RTWT(iLower) +	D24
		  ADRL(iLower) = ADRL(iLower) +	D24/RTWL
		  LL=1
		  
		  GO TO	80
		ENDIF   ! RTWT is not large enough to cross over cell
 1212		RTWT(iCell)	= RTWT(iCell) +	AWR(iCell)
		ADRL(iCell)	= ADRL(iCell) +	AWR(iCell)/RTWL
	END IF  !AWR >0
 80	CONTINUE
C
C   CALCULATE LENGTH OF YOUNG ROOTS AND ROOT DENSITY IN CELLS and total root mass***
C
      MaxRootDepth=1.0e6
      TotalRootWeight=0.0  
      DO 190, e = 1,NumEl ! added code to get
        i=KX(e,1)
        j=KX(e,2)
        k=KX(e,3)
        l=KX(e,4)
        YRL(e) = (YRL(e)*(1.0 - (ALPY*Step))) + ADRL(e)
        xx=1.0 - (ALPY*Step)
        xxy=xx*YRL(e)
        RUTDEN(e) = RUTDEN(e) + (YRL(e)*ALPY*Step/Area(e))
        TotalRootWeight=TotalRootWeight+RTWT(e)*CONVR   
        if (RTWT(e).GT.0) then
           MaxRootDepth=Amin1(y(i),y(j),y(k),y(l),MaxRootDepth)
         end if 
190   CONTINUE
      Total_Pcrs=Total_Pcrs+PCRS/(poprow/100.)*step*2.0   

      Return
10    Call errmes(im,il)
25    format (F12.6,1x,i6,2x,6G12.3,1x,2i5)
26    format (F12.6, 1x, 3i6,5F12.7, 2x,2i3)
27    format (F12.6,1x,F8.0,1x,F12.6,1x,3F9.3)
      Return
      End
