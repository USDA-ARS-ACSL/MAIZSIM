
CDT todo - need to take out weights and other input from the element geometry file.
C I am inputting dummy variables for the weights, we can 
      Subroutine Carbon_Partitioning()
      Include 'public.ins'
      Include 'puplant.ins'
      Dimension iLoc(4)
      Character InString*132
      Integer e,e1,SW10,SW20,SW30,SW40
      Real MECHR,RGCF1,RGCF2,RGCF3,RGCF4,RGCF5
c      Real PDWR_TOTAL   !for debugging
      Real xMean(NumNPD), yMean(NumNPD)  ! x and y coordinates of distance from seed
      integer InitialRootCell   ! array that holds element numbers of cells that are assigned an initial root density
      integer RootsInitiated
      double precision InitialRootWt
      Real RRRM,RRRY; ! root radial resistance for mature and young leaves (s-1 cm-1)     cccz bar day g-1
	Real RVRL;    !root vascular resistance s-1 cm -1                                   cccz bar day g-1
	Real RTWL;    !Root weight per unit length (g cm-1)                                 cccz g cm-1
	Real ALPM,ALPY; !relative growth rates for mature and young roots                   cccz day-1
	Real RtMinWtPerUnitArea; !minimum root weight per length (g cm-1)                   cccz g cm-2
	Real Wl,Wa,Wr,Wb; !Weighting parameters for root growth in grid used in Root_grow_c 
	integer ISink;        ! Determines mechanism of nitrogen uptake by root
      Real InitialRootCellArea
      integer t1,t2
c	0  'passive root mass uptake' no kinetics
c    1   Plant regulation diffusive and convective uptake
c    2  'convective diffusive uptake no plant regulation (Constant IMax)
c    3  'diffusive uptake only, no effect of water movement 
c    4  Michaelis-Menton uptake only
c    
      Real RootRadius;  !radius of root - avg in cm                                       cccz cm
	
      Common /CarbonPartition/ 
     !       VH2OC(NumNPD),                               ! cccz volumetric water content cm3 cm-3 
     !       CXT(NumNPD),                                 ! cccz oxygen concentration cm3 cm-3
     !       PSIRD(NumNPD),                               ! cccz soil water potential @ dawn bar
     !       TPRD(NumNPD),                                ! cccz root turger p @ dawn - threshold root growth bar
     !       RRRM,RRRY,                                   ! cccz concentration, g cm-3
     !       RVRL,                                        ! cccz root vascular resistance, bar day g-1 cm-1
     !       ALPM,                                        ! cccz potentional relative old/young root grow rate, day-1
     !       SC(NumNPD),                                  ! cccz axulliary for SR
     !       SR(NumNPD),                                  ! cccz soil resistance water flow bar day g-1
     !       OSMREG,                                      ! cccz osmoregulation
     !       PILOSM,                                      ! cccz leave osmotic potential @ dawn adjusted bar
     !       PSILZ,                                       ! cccz zero turgor leaf water potential bar
     !       PPSILT,                                      ! cccz previous leaf water potential shoot no growth bar
     !       PPSIL,                                       ! cccz previous leaf water potential bar
     !       SGTLT,                                       ! cccz shoot growing time lost proportion
     !       TPLT,                                        ! cccz leave turgor shoot no growth bar
     !       TPLD,                                        ! cccz leave turgor @ dawn bar
     !       BTPL,                                        ! cccz lowest leave turgor bar
     !       ET,                                           ! cccz not evapotranspiration
     !       InitialRootCell(NumNPD),
     !       RootsInitiated,
     !       iCellCount,
     !       InitialRootCellArea
    
      t=sngl(Time)
      If (lInput.eq.0) goto 11
c      open(108,file='rgcf.txt')
c      write(108,'(5A12)') 'Time','e','rgcf1','rgcf2','rgcf3'
C
C  Reading of the input files and initial calculations 
C
! initialize variables to deefault values
          RRRM=166.7D0;
		RRRY=31.3D0;
		RVRL=0.73D0;
		ALPM=0.350D0;
		ALPY=0.04D0;
		RTWL=1.059e-004; 
          RtMinWtPerUnitArea=2.0e-004;
          
      Do n=1,NumNP
        YRL(n)=0.0D0                          ! length of young root in soil cell cm (old version)
        RGCF(n)=0.0D0                         ! root growth reduction factor
        iFavRoot(n)=0.0D0
cccz start to move the nodeArea to GRID_BND file
cccz        nodeArea(n)=0.0D0
        RMassM(n)=0.0D0                       ! density of old root in soil cell g/cm2 (new version)
        RMassY(n)=0.0D0                       ! density of young root in soil cell g/cm2 (new version)
        RDenM(n)=0.0D0                        ! density of old root in soil cell cm/cm2 (new version)
        RDenY(n)=0.0D0                        ! density of young root in soil cell cm/cm2 (new version)
      Enddo
      
     	Open(40,file = VarietyFile, status = 'old',ERR = 10)

      im=220
      il=0
15       Read (40,'(A132)') InString
          if (InString(1:14).ne.'[SoilRoot]') goto 15
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) RRRM,RRRY,RVRL
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) ALPM,ALPY,RTWL,RtMinWtPerUnitArea
      close(40)

c DT 12/7/2018 removed element file reading     Only need to initialize root weight now
      Do n=1,NumNP
        RTWT(n)=0         ! dry weight of root for each soil node (area)
      Enddo
      Close(41)

cccz Assume that RTWT is node based root mass (g)
      TotalRootWeight=0.0D0      
      TotalArea=0.0D0
      do n=1,NumNP
          TotalRootWeight=TotalRootWeight+RTWT(n)*nodeArea(n)        ! RTWT root weigth in each cell g/cm^2
          TotalArea=TotalArea+nodeArea(n)
          RMassY(n)=RTWT(n)                        ! cccz initialized yong root density for each node g cm-2
          RDenY(n)=RMassY(n)/RTWL                  ! cccz initialized yong root density for each node cm cm-2
      enddo
      
      EORSCS=0.0                                  ! cumulation of EORSCF
      AWUPSS=0.0                                  ! cumulation of water uptake
      Total_Eor=0.0
      Total_Pcrs=0.0                              ! PCRS actual carbon rate  g day-1
      icount=0
      iCellCount=0
      RootsInitiated=0     
      xMin=1.E+32
      xMax=-1.E+32
      InitialRootCellArea=0.0
      PSISM=0.0D0                                 ! water potential for cells for water uptake bar
      
c      do n=1,NumNP
cccz RMassY has the unit g cm-2
cccz So TotalRootWt is not a weight but a sum of "mass-density"
c          TotalRootWt=TotalRootWt+RMassY(n)
c      enddo
      
      Do n = 1,NumNP
          PSIS(n)  = hNew(n)/1000.0D0  ! convert from cm to bars
          xDist = (x(n)-xBSTEM)
          yDist = (y(n)-yBSTEM)
          RVR(n)=RVRL*SQRT(xDist*xDist+yDist*yDist)   ! all the radius distance is use to calculate the radius resistance.
                                                      ! root resistance RVR: bar day g-1
        If ((abs(xDist).lt.2).and.(abs(yDist).lt.5)) then
           InitialRootCell(n)=1
           iCellCount=iCellCount+1
           InitialRootCellArea=InitialRootCellArea+nodeArea(n)
        Endif 
        If(RMassY(n).gt.0.0D0) then
            PSISM=PSISM+PSIS(n)
            icount=icount+1
        Endif

        MECHR = 80.03*abs(psis(n))**0.53
     &     *EXP(-8.08*(1.66 - BlkDn(MatNumN(n))))

        IF (MECHR.LE.0.0) MECHR = 0.0
        TPRD(n) = 2.0 + (0.1*PSIS(n) + (0.7*MECHR))
        PSIRD(n) = PSIS(n)
      Enddo

c   EOMult is defined in Crop since it is a Plant variable
c      EOMult=1.
c      If(xBSTEM.le.xMin.or.xBSTEM.ge.xMax) EOMult=0.5
cccz zhuangji EOMult is to mark if the plant is on the edge.
      
      PSISM=PSISM/max(1.0,float(icount))
      PPSIL = -0.5
      PPSILT = -0.5
      TPLT = 2.0
      PSILD = -0.5            ! leave water potential @ dawn bar
c since in the MAIZESIM model, the reference leaf water potential in the CGas_exchange::setleafpEffect
c function is set to be -1.9, so the initial value of leaf water potential should be assigned a higher
c value Yang 8/22/06
      IF (PSISM.LT.PSILD) PSILD = PSISM  ! leave pressure should be less than soil water potential
      PILOSM = -9.0D0                    ! PILOSM leave osmotic potential at dawn adjusted bar
      OSMREG = -0.25D0
      OSMFAC = 0.5
      PILD = ((PSILD + 2.0D0)*OSMFAC) + PILOSM  ! leaf osmotic pressure bar @ dawn
      TPLD = PSILD - PILD
      TPL = PPSIL - PILD
      BTPL = TPLD
      SGTLT = 0.0D0                   
11    continue
C
c  Routine calculations
C
c                     Cell variables
c     This piece of code averages nodal values over cells and
c     calculates the cell hydraulic conductivity.
c     It has no shoot-associated parameters
C
CYAP  No Plant = No Root Activity
      If((NShoot.eq.0).or.(isGerminated.eq.0)) Return
c       PDWR_TOTAL=0.0    ! for debugging purposes  AD and zhuangji  delete it
CYAPEND
      If ((isGerminated.eq.1).and.(RootsInitiated.eq.0)) then
         RootsInitiated=1
         InitialRootWt=InitialRootCarbo/InitialRootCellArea   ! calculates average root wt density
         ! Distribute this mass among cells with initial root distribution
         areaS=0
          Do n=1,NumNP
             if (InitialRootCell(n).gt.0) then
               RMassY(n)=RMassY(n)+InitialRootWt
               RDenY(n)=RMassY(n)/RTWL
               TotalRootWeight=TotalRootWeight+RMassY(n)*nodeArea(n)
               areaS=areaS+nodeArea(n)
              endif
           enddo
        Endif
CYAPEND
      Do n=1,NumNP
        PSIS(n)=hNew(n)/1000.0D0
        COND(n)=Con(n)
      Enddo

      
c     EO : potential transpiration rate
c     EO is in g/plant. Multiply g/plant x plant/m2 x 1 m2/10000 cm2 x rowsp
c     EOR : = EO*ROWSP if the plant is not at the border of the soil slab,
c           = 0.5 * EO*ROWSP, if the plant is at the left or right border.
c      EOR = EO*poprow/rowsp*100*rowsp*EOMult/10000
c        poprow*100/rowsp gives plants m-2
c      ET_demand is units of grams m-2 day-1
c     TPot is g plant-1 hr- poprow*100/Rowsp converts to g cm hr
c     EoMult*rowsp/10000/24 converts to g slab-1 (rowsp * 1cm strip for the finite element domain)
       EOR = ET_demand*poprow*(100/rowsp)*(EOMult*rowsp/10000) ! potential transportation rate per slab
c       EOR = ET_demand
c      PCRL = PCRL*12*poprow/rowsp*100*rowsp*EOMult/10000
c      PCRQ = PCRQ*12*poprow/rowsp*100*rowsp*EOMult/10000
      
c   PCRL in MAIZESIM is in the unit g/plant/hour, so mulitply 24 first to convert it to g/plant/day
c   then multiply it by plant/(m^2) to convert it to g/(m^2)/day
c   finally multiply it by rowsp*1(cm) to convert it into the soil slab we are dealing with Yang 9/24/06
             
c     
C     ---------- Taken form R U T G R O   of  G L Y C I M
C
C  CALCULATE SOIL MECHANICAL RESISTANCE TO ROOT GROWTH IN EACH CELL
C
      DO n=1,NumNP
        MECHR = 5.4*abs(psis(n))**0.25
     &     *EXP(-10.58*(1.7 - BlkDn(MatNumN(n))))
        IF (MECHR.LE.0.0D0) MECHR = 0.0D0
C
C  CALCULATE ROOT TURGOR PRESSURE IN EACH CELL AT DAWN MINUS THRESHOLD
C  TURGOR FOR GROWTH
C
        IF (ITIME.EQ.IDAWN) THEN
          TPRD(n) = 2.0D0 + (0.1D0*PSIS(n) + (0.7D0*MECHR))
          PSIRD(n) = PSIS(n)
        ENDIF
C
C  CALCULATE PROPORTIONIAL REDUCTION OF ROOT GROWTH CAUSED BY
C  MECHANICAL RESISTANCE AND SOIL WATER POTENTIAL ASSUMING ROOT
C  IS AT SAME POTENTIAL AS SOIL
C
         RGCF1 = ((TPRD(n)-MECHR)/2.0)-((PSIRD(n)-PSIS(n))/4.0)
C
C     CALCULATE PROPORTIONAL REDUCTION OF ROOT GROWTH CAUSED BY
C     SOIL TEMPERATURE (using a beta function)
C
        RGCF2=tempStress(Tmpr(n))
        

C
C   CALCULATE PROPORTIONIAL REDUCTION OF ROOT GROWTH CAUSED BY
C   PARTIAL PRESSURE OF SOIL OXYGEN
C
        RGCF3= oxyStress(g(n,2)*ugGasCm3air_to_ppm(2)/10000.)
        
        t1=(RGCF3.lt.1.0)
        
cccz the distance effects
        xDist = abs(x(n)-xBSTEM)
        yDist = abs(y(n)-yBSTEM)
        RGCF4=max(1.0D0-0.001*xDist,0.0D0)
     !     *max(1.0D0-0.0001*yDist,0.0D0)
cccz this effect is turned off in the gas one
        
C
C   CALCULATE PROPORTIONAL REDUCTION OF ROOT GROWTH FROM THE
C   PHYSICAL CAUSE THAT IS MOST LIMITING
C
cccz Zhuangji: AD 
c        RGCF5=(1.0-AMIN1(1.0,RMassY(n)/0.03))
         RGCF5=(1.0-AMIN1(1.0,(RMassY(n)+RMassM(n))/0.02))
         RGCF5=1.0D0
C this one for nitrogen. Use the ratio between actual uptake of N and IMax
    
         dif=amax1(0.0,Conc(n,1)-CMin0(1))
	   RGCF6 =  amax1(0.3,dif/ConstK(1))
	   RGCF6 = AMIN1(1.0, RGCF6)

c ccz zhuangji:AD removed the RGCF3 why? DT - O2 was not simulated 
C dt need to do more testing for RGCF6, add later
        RGCF(n) = AMIN1(RGCF1, RGCF2, RGCF3, RGCF4)
        
        RGCF(n) = AMAX1(RGCF(n),0.0)
        RGCF(n) = AMIN1(RGCF(n),1.0)

C
C   CALCULATE POTENTIAL RATE OF CHANGE IN ROOT DRY WEIGHT IN EACH
C   SOIL CELL

c
c  rtwt is gr in the cell
c  alpm a Percentage per day
c  RGCF has no units
c  pdwr = gr d-1 per cell
c  dt 03/18/11 testing the incorporation of area

       if ((isEmerged.lt.1).and.(isGerminated.gt.0)) RGCF(n)=1.0
 
cccz only new root or the total root? Dt - use total root for now
        PDWR(n)=(RMassY(n)+RMassM(n))*ALPM*RGCF(n)*nodeArea(n)
c        PDWR(n)=RMassY(n)*ALPM*RGCF(n)*nodeArea(n)
           
       ENDDO
C
C DT GLYCIM sorts on PDWR which is potential for root growth. RGCF is already
C    taken into account.
      Call SORT(NumNP,iFavRoot,RGCF,RVR)
C
C
C   Taken from TRADE  of GLYCIM
C NOTE 10.09
C   THIS ROUTINE MAINTAINS A FUNCTIONAL BALANCE BETWEEN ROOT AND
C   SHOOT BY ALLOCATING CARBON TO GROW ROOT AS NECESSARY TO MEET
C   TRANSPIRATION DEMAND. NEW ROOT HAS LESS RESISTANCE TO WATER UPTAKE
C   SINCE TRANSPIRATION RATE, PHOTOSYNTHESIS RATE AND HENCE POTENTIAL
C   ROOT GROWTH RATE ARE ROUGHLY PROPORTIONAL, LEAF WATER POTENTIAL IS
C   NOT AFFECTED MUCH BY TRANSPIRATION RATE. HOWEVER, LEAF WATER
C   POTENTIAL MUST FALL AS TRANSPIRATION INCREASES BECAUSE THE
C   THRESHOLD TURGOR PRESSURE FOR STOPPING SHOOT GROWTH AND DIVERTING
C   CARBOHYDRATES TO ROOTS DECREASES WITH TIME. HENCE THE ROUTINE ALSO
C   PREDICTS CHANGES IN TOTAL POTENTIAL, OSMOTIC POTENTIAL AND TURGOR
C   PRESSURE IN LEAVES DURING THE DAY.
C    
      DO n = 1,NumNP
        SR(n)   = 0.0
        WUPM(n) = 0.0
        WUPN(n) = 0.0
        WUPT(n) = 0.0
        DRL(n)  = 0.0
cS 07/05/2007 DT These variables are needed for nitrogen sink code    
        FUP(n,1) = 0.0 !WATER FLUX TO THE OLD ROOTS or old roots density 
	  FUP(n,2) = 0.0 !WATER FLUX TO THE NEW ROOTS or new roots density       
        VUP(n,1) = 0.0 !WATER FLOW TO THE OLD ROOTS
	  VUP(n,2) = 0.0 !WATER FLOW TO THE NEW ROOTS
cE 7/05/2007 DT 
      ENDDO
      SGTLI = 0.0
C
C    THE THRESHOLD LEAF WATER POTENTIAL AND TURGOR PRESSURE THAT
C    WILL JUST PREVENT ALL SHOOT EXPANSION ARE CALCULATED
C
C    CALCULATE DAWN LEAF WATER POTENTIALS ***
C
      IF (ITIME.EQ.IDAWN) THEN
c        PSILD = -2.0
c since in the MAIZESIM model, the reference leaf water potential in the CGas_exchange::setleafpEffect
c function is set to be -1.9, so the initial value of leaf water potential should be assigned a higher
c value Yang 8/22/06
        PSILD = - 0.5
        IF (PSISM.LT.PSILD) PSILD = PSISM
        PILOSM = PILOSM - OSMREG
        IF (PILOSM.LT.-12.0) PILOSM = -12.0
        IF (PILOSM.GT.-9.0) PILOSM = -9.0
        OSMREG = -0.25
        OSMFAC = 0.5
        PILD = ((PSILD + 2.0)*OSMFAC) + PILOSM
        TPLD = PSILD - PILD
        TPL = PPSIL - PILD
        BTPL = TPLD
        SGTLT = 0.0
       ENDIF
C
C   CALCULATE THE THRESHOLD LEAF WATER POTENTIAL AND TURGOR
C   PRESSURE THAT WILL JUST PREVENT ALL SHOOT GROWTH
C 

      D11 = MAX (TPL, TPLT)
      D12 = (D11 - 2.0D0)*(1.0D0 - EXP(-0.7D0*Step))
      IF (D12.LE.0.0D0) D12 = 0.0D0
      D13 = MAX (PPSIL, PPSILT)
      OSMFAC = 0.5D0
      PSILT = D13 - (D12/(1.0D0 - OSMFAC))
      PPSILT = PSILT
      TPLT = D11 - D12
      PSILZ = PSILT-(TPLT/(1.0D0 - OSMFAC))    ! leaf water potential @ zero turgor bar
C
C   WATER UPTAKE BY YOUNG AND OLD ROOTS UNDER VARIOUS
C   CIRCUMSTANCES IS CALCULATED
C
C   CALCULATE RADIUS OF SOIL CYLINDER THROUGH WHICH WATER MUST
C   TRAVEL TO REACH ROOTS. THIS IS APPROXIMATED BY A SIMPLE
C   FUNCTION OF SOIL WATER POTENTIAL TO AVOID ITERATION
C   SC IS THE HALF THE MEAN DISTANCE BETWEEN NEIGHBORING ROOTS
C   WE NEED TO SEE IF THIS CAN BE CALCULATED FROM ROOT DENSITY
C
      WUPMS = 0.
      WUPSI = 0.
      DO 50, n=1,NumNP
        
        D14 = PSIS(n) - PSILT               ! leaf pressure should be less to ensure water uptake
        IF (D14.GT.0.0D0) THEN
          SC(n) = 0.017D0 - (PSIS(n)*0.5D0)
C
C   CALCULATE RESISTANCE TO WATER FLOW IN SOIL CELL
C
        IF (COND(n).GT.0.0D0) THEN
            SR(n) = ALOG(SC(n)*SC(n)/3.0E-4)/(4.0*3.1416*COND(n)*1019.7)
C  3e-4 is 0.017*0.017; Equation is from Gardener, 1960
C
C   CALCULATE RATE OF WATER UPTAKE BY EXISTING ROOTS IF LEAF WATER
C   POTENTIAL IS AT THE THRESHOLD
C         

c If LAI=0, then all the resistance comes from soil;
c when LAI !=0, then the total resistance equals to resistance in the soil plus
c resistance of the xylum; here the xylum resistance is adopted as a constant: 6.67*1-E-3 molH2O m-2(leaf) s-1 MPa-1
c to convert this value into the resistance value for a slab, we have to divide it by LAI and the area of the slab and 
c convert MPa into bar, this is how we have 19.26//(LAI*rowsp)
c WUPM is rate of water uptake by mature roots g/day
c WUPN water uptake by new roots
            IF(LAI.EQ.0.0) Then
              WUPM(n) = (((PSIS(n) - PSILT)*nodeArea(n)
     &         *RDenM(n))/(SR(n) + RRRM + RVR(n)))
     &         *0.02275*(24.5 + Tmpr(n))
     
              WUPN(n) = (((PSIS(n) - PSILT)*nodeArea(n)
     &         *RDenY(n))/(SR(n) + RRRY + RVR(n)))
     &         *0.02275*(24.5 + Tmpr(n))
             ELSE
             
             WUPM(n)=0.0
             if (RDenM(n).GT.0.0) then
               WUPM(n)=(PSIS(n)-PSILT)/((SR(n)+RRRM+RVR(n))
     &           /(nodeArea(n)*RDenM(n)*0.02275*(24.5+Tmpr(n)))
     &           +19.26/(LAI*ROWSP))
             endif
             WUPN(n)=0.0
             if (RDenY(n).GT.0.0) then
             WUPN(n)=(PSIS(n)-PSILT)/((SR(n)+RRRY+RVR(n))
     &        /(nodeArea(n)*RDenY(n)*0.02275*(24.5+Tmpr(n)))
     &         +19.26/(LAI*ROWSP))
             endif
     
             END IF       

C  NOTE 10.17
c WUPMS total water uptake by mature roots
c WUPSI sum of water potential * rate of extraction in node
             WUPMS = WUPMS +(WUPM(n) + WUPN(n))
             WUPSI = WUPSI +(PSIS(n)*(WUPM(n) + WUPN(n)))
              
          ENDIF
        ENDIF
 50   CONTINUE
C
C  CALCULATE GROWTH BY NEW ROOTS IF ROOTS GREW AT POTENTIAL
C  RATE DURING PAST PERIOD. START WITH SOIL CELL WHERE CONDITIONS
C  ARE MOST FAVORABLE FOR GROWTH
C

      PCRTS = 0.       ! potential total root carbon use  g day-1 could include pcrq
      WUPTS = 0.       ! new root water uptake when all carbon goes to the root   g day-1
      WUPRS = 0.       !rate if water uptake when shoot growth potential has been satisfied.
      SW10 = 0
      DO 90, n=1,NumNP
        n1=iFavRoot(n)
        PDRL = PDWR(n1)/RTWL            ! cccz potential rate root length cm/day
        DRL(n1) = (PDRL + PPDRL(n1))*Step/2.0D0  ! unit cm
        PPDRL(n1) = PDRL
C
C   DETERMINE TOTAL WATER UPTAKE BY ALL NEW ROOTS AND WEIGHTED
C   AVERAGE SOIL WATER POTENTIAL WHEN:
C   I. ROOTS GET ALL THE CARBON TRANSLOCATED
C
        WUPT(n1)=0.0    ! rate of water uptake when water potential prevents shoot growth
        IF(LAI.EQ.0.0) Then
        WUPT(n1) = (((PSIS(n1) - PSILT)*DRL(n1))
     &   /(SR(n1) + RRRY + RVR(n1)))
     &   *0.02275*(24.5 + Tmpr(n1))
        Else
          if (DRL(n1).GT.0.0) then
            WUPT(n1)=(PSIS(n1)-PSILT)/((SR(n1)+RRRY+RVR(n1))
     &        /(DRL(n1)*0.02275*(24.5+Tmpr(n1)))
     &       +(19.26/(LAI*ROWSP)))
          endif
        EndIF

        IF (WUPT(n1).LE.0.0) WUPT(n1) = 0.0

        WUPTS = WUPTS + WUPT(n1)
        WUPSI = WUPSI + (PSIS(n1)*WUPT(n1))
        
        PCRTS = PCRTS + (PDWR(n1)*CONVR)    ! g day-1
C
C   II. ROOTS GET ONLY THE CARBON LEFT OVER FROM SHOOT GROWTH
C
c if carbon allocation changes, have to do it here too
        IF (PCRL.LE.PCRTS.AND.SW10.NE.1) THEN
          SW10 = 1
          D16 = 0.0D0
          IF (PDWR(n1).GT.0.0D0) D16 = WUPT(n1)*(1.0D0
     &     - ((PDWR(n1)*CONVR) - PCRTS + PCRL)/(PDWR(n1)*CONVR))
          WUPRS = WUPTS - D16
C NOTE 10.23
          IF ((WUPMS + WUPRS).GT.0.0D0) PSISM = (WUPSI - (PSIS(n1)*D16))
     &     /(WUPMS + WUPRS)                                               ! this is the average psi
        ENDIF
        IF (PCRQ.LE.PCRTS) GOTO 1111
 90   CONTINUE
C
      IF ((WUPMS + WUPTS).GT.0.0D0) PSIST = WUPSI/(WUPMS + WUPTS)         ! this is the average psi
      IF (SW10.EQ.0) THEN
        WUPRS = WUPTS
        PSISM = PSIST
      ENDIF
      GOTO 1112
C
 1111 D17 = 0.0
      IF (PDWR(n1).GT.0.0) D17 = WUPT(n1)*(1.0 - 
     &   ((PDWR(n1)*CONVR)- PCRTS + PCRQ)/(PDWR(n1)*CONVR))
      WUPTS = WUPTS - D17
C NOTE 10.24
      IF ((WUPMS + WUPTS).GT.0.) PSIST = (WUPSI - (PSIS(n1)*D17))
     & /(WUPMS + WUPTS)                                               ! this is the average psi
C
C   CALCULATE WATER UPTAKE RATE IF LEAF WATER POTENTIAL HAS NOT
C   RISEN ABOVE THE THRESHOLD WHICH JUST PREVENTED ALL SHOOT
C   EXPANSION IN THE LAST PERIOD
C
1112  WUPDS = (WUPMS + WUPRS)*(PSISM - PPSILT)/(PSISM - PSILT)
C
C   CALCULATE WATER UPTAKE RATE IF LEAF TURGOR PRESSURE = 2 BARS
C   OR = 0 BARS
C
      DPSI02 = 2.0/(1.0 - OSMFAC)
      WUP2S = (WUPMS + WUPTS)*(PSIST - PSILZ - DPSI02)/(PSIST - PSILT)
      WUP0S = (WUPMS + WUPTS)*(PSIST - PSILZ)/(PSIST - PSILT)
cccz??? the usage of WUPTS and WUPRS is confusing.        

C
C   CALCULATE POTENTIAL TRANSPIRATION RATE ALLOWING FOR NIGHT-TIME
C   STOMATAL CLOSURE.  THE STOMATA ARE ASSUMED TO CLOSE TIGHTER
C   WHEN THE PLANT IS EXPERIENCING WATER STRESS
C

c**************************************************************
c Stomatal openning and closing is determined in MAIZESIM, 
c so no need to for 2dsoil to determine closure/openning of
c stomatal Yang 8/24/06

c*****************************************************************
C
C   DETERMINE HOW POTENTIAL TRANSPIRATION RATE IS RELATED TO THE
C   VARIOUS POTENTIAL WATER UPTAKE RATES AND GO TO THE APPROPRIATE
C   SECTION OF THE CODE. CALCULATE LEAF WATER POTENTIAL,STOMATAL
C   CONDUCTANCE,THE PROPORTION OF TIME FOR WHICH THE SHOOT GROWS,
C   THE AMOUNT OF CARBON ACTUALLY NEEDED FOR ROOT GROWTH AND THE
C   AMOUNT OF OSMOREGULATION THAT WILL OCCUR
      D18 = WUPMS + WUPTS                                       ! all water uptake from mature roots and new roots
cccz??? I suppose this is WUPRS because that is the actural water uptake (in case carbon is not enough)
C
C   POTENTIAL TRANSPIRATION RATE IS SUCH THAT THE DEMAND FOR WATER
C   UPTAKE CAN BE MET WITHOUT SHOOT WATER POTENTIAL FALLING TO THE
C   POINT WHERE THE SHOOT LOSES TURGOR AND STOPS GROWING
C
      IF (EOR.LT.WUPDS) THEN                                      ! relatively small EOR
        PSIL_ = PSISM - ((PSISM - PSILT)*EOR/(WUPMS + WUPRS))
         SGT = 1.0
C
C   POTENTIAL TRANSPIRATION RATE IS SUCH THAT THE DEMAND FOR WATER
C   UPTAKE CANNOT BE MET WITHOUT SHOOT WATER POTENTIAL FALLING TO
C   THE POINT WHERE THE SHOOT LOSES TURGOR AND STOPS GROWING FOR
C   PART OF THE PERIOD.  WHEN THIS HAPPENS THE ROOTS GROW
C
      ELSE IF (EOR.LT.D18) THEN                                   ! relatively medium EOR
        D19 = (PPSIL - PSILT)/2.0
        PSIL_ = (PPSIL + PSILT)/2.0
c        SCF = 1.0
        IF (D19.LE.0.0) THEN
          SGT = 1.0
          GOTO 1113
        END IF
C
C   CALCULATE LEAF WATER POTENTIAL ITERATIVELY
C

      DO 100, J = 1,3
        WUPGS = (PSISM - PSIL_)/(PSISM - PSILT)
     &     *(WUPMS + WUPRS + ((WUPTS - WUPRS)
     &     *(PPSIL - PSIL_)/(D19*2)))
        IF (EOR.GT.WUPGS) PSIL_ = PSIL_ - (D19/2**J)
        IF (EOR.LT.WUPGS) PSIL_ = PSIL_ + (D19/2**J)
 100  CONTINUE
C
      SGT = (PSIL_ - PSILT)/(D19*2.0)
C
C   POTENTIAL TRANSPIRATION RATE IS BETWEEN THE POSSIBLE WATER
C   UPTAKE RATES WHEN LEAF TURGOR PRESSURE IS (1) SO LOW THAT ALL
C   SHOOT GROWTH IS PREVENTED, THE ROOTS GROW CONTINUOUSLY, AND
C   AND THE LEAVES OSMOREGULATE, AND (2) AT 2 BARS, WHERE STOMATAL
C   CLOSURE IS ASSUMED TO START
C
      ELSE IF (EOR.LT.WUP2S) THEN
        PSIL_ = PSIST - ((PSIST - PSILT)*EOR/(WUPMS + WUPTS))
        SCF = 1.0
        SGT = 0.0
        OSMREG = 0.25
C
C   SINCE POTENTIAL TRANSPIRATION RATE IS GREATER THAN THE
C   POSSIBLE WATER UPTAKE RATE WHEN LEAF TURGOR PRESSURE IS 2 BARS
C   STOMATAL CLOSURE WILL DECREASE TRANSPIRATION RATE.
C   CALCULATE LEAF WATER POTENTIAL AND STOMATAL CLOSURE FACTOR
C
CDT  04-25-2011 under severe stress, EOR, WUP2S, and WUP0S were all zero resulting
CDT  10-07-2011 under severe stress, EOR, WUP2S, and WUP0S were all zero resulting
C    in an undefined value for PSIL_. I added a case for this and set PSIL_ to -15.
      ELSE IF ((EOR+WUP2S+WUP0S).gt.0) THEN
        PSIL_ = PSILZ + DPSI02 - ((EOR - WUP2S)/((EOR/DPSI02)
     &   + ((WUP0S - WUP2S)/DPSI02)))
        SCF = 1.0 - ((PSILZ + DPSI02) - PSIL_)/DPSI02
        SGT = 0.0
        OSMREG = 0.25
      ELSE    ! case when all else fails
        PSIL_=-15
        SCF=0
        SGT=0
        OSMREG=0.1
      END IF
C
C   LEAF WATER POTENTIAL CANNOT RISE ABOVE ITS DAWN VALUE
C
 1113 IF (PSIL_.GT.PSILD) PSIL_ = PSILD
C
C  CALCULATE LEAF TURGOR PRESSURE
C
      PTPL = TPL
      TPL = ((PSIL_ - PPSIL)*(1.0 - OSMFAC)) + TPL
      IF (BTPL.GT.TPL) BTPL = TPL
      PPSIL = PSIL_
C
C  IF TURGOR PRESSURE IS LESS THAN 2 BARS THE SHOOT DOES NOT GROW
C
      IF (TPL.LT.2.0) THEN
        SGTLI = 1.0
        SGT = 0.0
        PCRS = PCRQ
C
C   IF TURGOR PRESSURE IS DECREASING AND IS LESS THAN 5 BARS
C   SHOOT GROWTH POTENTIAL IS FURTHER REDUCED
C
      ELSE IF (TPL.LT.PTPL) THEN
        IF (TPL.LE.5.0) SGT = SGT*(TPL - 2.0)/3.0
        SGTLT = SGTLT + 1.0 - SGT
        PCRS = PCRL
        IF (SGT.LT.1.0) PCRS = PCRL + ((PCRQ - PCRL)*(1.0 - SGT))
C
C     *** IF TURGOR PRESSURE IS INCREASING AND IS ABOVE 2.0 BARS,
C         SHOOT GROWTH POTENTIAL IS INCREASED. ***
C
      ELSE
        IF (BTPL.LT.2.0) BTPL = 2.0
        PCRS = PCRL
      END IF
     

      Return
10    Call errmes(im,il)
25    format(f9.2,1x,i6,1x,3F9.4)
      Return
      End
*********************************************************************
      SUBROUTINE SORT(NCR,LOCATE,VALUE,RVR)
C     --------------------------( SORT )--------------------------------
C     THIS SUBROUTINE SORTS THE CORRESPONDING ARRAYS "LOCATE" AND "VALUE"
C     SO THAT THE ELEMENTS OF "LOCATE" GIVE NUMBERS OF ELEMENTS OF "VALUE" 
C     IN ORDER FROM HIGHEST TO LOWEST VALUE. If neighbor values do not 
C     differ, preference is given to value with lower RVR.
C     IT USES THE QUICKSORT PROCEDURE PROPOSED BY C.A.R. HOARE AND
C     DESCRIBED IN MANY PROGRAMMING TEXTBOOKS.  Variables used only inside
C     this subroutine are not listed in the dictionary for GLYCIM.
C     ------------------------------------------------------------------
      DIMENSION LOCATE(NCR),VALUE(NCR),MARRAY(3500),NARRAY(3500),
     !          VALAux(3500),RVR(NCR)
C
      M = 1
      N = NCR
      KOUNT = 0
         DO I=1,NCR
         LOCATE(I)=I
         VALAUX(I)=VALUE(I)
         ENDDO
C
 1111 IF (M.LT.NCR) THEN
        IF (M.LT.N) THEN
C
C     *** SORT PART OF ARRAY BOUNDED BY M AND N ***
C
          I = M - 1
          J = N
          REF = VALAux (N)
 1112   CONTINUE
        I = I + 1
        IF (VALAux(I).GT.REF) GO TO 1112
 1113 CONTINUE
      J = J - 1
      IF (J.GT.1.AND.VALAux(J).LT.REF) GO TO 1113
      IF (I.LT.J) THEN
        D11 = VALAux(I)
        VALAux(I) = VALAux(J)
        VALAux(J) = D11
        M11 = LOCATE(I)
        LOCATE(I) = LOCATE(J)
        LOCATE(J) = M11
      ELSE
        GO TO 1114
      END IF
      GO TO 1112
 1114 D12 = VALAux(N)
      VALAux(N) = VALAux(I)
      VALAux(I) = D12
      M11 = LOCATE(N)
      LOCATE(N) = LOCATE(I)
      LOCATE(I) = M11
C
C     *** MANIPULATE STACK OF VALUES OF M AND N ***
C
      KOUNT = KOUNT + 1
      MARRAY(KOUNT) = I + 1
      NARRAY(KOUNT) = N
      N = I - 1
      ELSE
      M = MARRAY(KOUNT)
      N = NARRAY(KOUNT)
      KOUNT = KOUNT - 1
      END IF
      GO TO 1111
      END IF
C
C     *** FURTHER SORT THE CELLS SO THAT WITHIN GROUPS HAVING THE SAME
C       VALUE, THOSE NEARER THE PLANT STEM HAVE PRIORITY IN THE LIST ***
C
      I = 1
 1115 IF (I.LT.NCR) THEN
        IF (VALAux(I).GT.VALAux(I + 1)) THEN
          I = I + 1
        ELSE
c          IF (LOCATE(I).GT.LOCATE(I + 1)) THEN
          IF (RVR(LOCATE(I)).GT.RVR(LOCATE(I + 1))) THEN
            M12 = LOCATE(I)
            LOCATE(I) = LOCATE(I + 1)
            LOCATE(I + 1) = M12
            I = I - 1
            IF (I.LT.1) I = 1
          ELSE
            I = I + 1
          END IF
        END IF
        GO TO 1115
      END IF
C
      M=1
      RETURN
      END
C 
       Real Function oxyStress(O2)
C
CDT       function to simulate the effect of low oxygen content
C  in flooded soils on root growth of aradopsis, it is from
C van Dongen JT,et al.  Ann Bot. 2009 Jan;103(2):269-80.
C  O2 is input as a percentage , 0 is full stress, 1 is no stress
C  Using a MM equation
C returns a number from 0 to 1
       real O2,k, V, tmp
       k=4.95
       V=121.5
       tmp=(V * O2)/(k + O2)
       tmp=min(tmp,100.0)
       tmp=max(tmp,0.0)
       oxyStress=tmp/100.0
       
       RETURN
       END
       
       Real Function TempStress(t)
       
      
CDT    beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
       implicit none
	 real t, R_max, t_opt, t_ceil,f, g, alpha
       real t_base, beta
       
       t_base=5.0
       beta=1.0
       t_opt=29.1
       t_ceil=40.0
       R_max=1.0
       
       if((t.le.t_base).or.(t.ge.t_ceil)) then
         TempStress=0.0
         return
       endif
c   not sure this is necesssary, kind of a safety thing       
       if ((t_ceil.le.t_opt).or.(t_opt.le.t_base)) then
        TempStress=0.0
        return
       endif

	f = (t - t_base) / (t_opt - t_base)
	g = (t_ceil - t) / (t_ceil - t_opt)
	alpha = beta*(t_opt - t_base) / (t_ceil - t_opt)
      TempStress= R_max*(f**alpha)*(g**beta)
      TempStress=max(0.0,TempStress)
      TempStress=min(1.0,TempStress)
      return
      End