*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Subroutine Carbon_Partitioning()
      Include 'public.ins'
      Include 'puplant.ins'
      Dimension iLoc(4)
      Character InString*132
      Integer e,e1,SW10,SW20,SW30,SW40,iCellCount
      Real MECHR,RGCF1, RGCF2, RGCF3,RGCF4,RGCF5
      Real xMean(numElD), yMean(numElD)  ! x and y coordinates of distance to center of cell from seed
      integer InitialRootCell   ! array that holds element numbers of cells that are assigned an initial root density
      integer RootsInitiated
      Real InitialRootCellArea,InitialRootWt
      Real PDWR_TOTAL   !for debugging
      Real RRRM,RRRY; ! root radial resistance for mature and young leaves (s-1 cm-1)
	Real RVRL;    !root vascular resistance s-1 cm -1
	Real RTWL;    !Root weight per unit length (g cm-1)
	Real ALPM,ALPY; !relative growth rates for mature and young roots
	Real RtMinWtPerUnitArea; !minimum root weight per length (g cm-1)
	Real Wl,Wa,Wr,Wb; !Weighting parameters for root growth in grid used in Root_grow_c
	integer ISink;        ! Determines mechanism of nitrogen uptake by root
!	   /*      0  'passive root mass uptake' no kinetics
!             1   Plant regulation diffusive and convective uptake
!             2  'convective diffusive uptake no plant regulation (Constant IMax)
!             3  'diffusive uptake only, no effect of water movement 
!             4  Michaelis-Menton uptake only

!   */
      Real RootRadius;  !radius of root - avg in cm 
	
      Common /CarbonPartition/ 
     !       VH2OC(NumElD),
     !       CXT(NumElD),PSIRD(NumElD),TPRD(NumElD),
     !       RRRM,RRRY,RVRL,ALPM,Conc_e(NumElD),
     !       SR(NumElD),
     !       OSMREG,
     !       PILOSM,PSILZ,
     !       PPSILT,PPSIL,SGTLT,TPLT,TPLD,
     !       BTPL,ET,InitialRootCell(NumElD),
     !       RootsInitiated,iCellCount,InitialRootCellArea

      t=sngl(Time)
      If (lInput.eq.0) goto 11
c      open(108,file='rgcf.txt')
c      write(108,'(5A12)') 'Time','e','rgcf1','rgcf2','rgcf3'
C
C  Reading of the input files and initial calculations 
C
! initialize variables to deefault values
        RRRM=166.7;
		RRRY=31.3;
		RVRL=0.73;
		ALPM= 0.350;
		ALPY= 0.04;
		RTWL=1.059e-004; 
        RtMinWtPerUnitArea=2.0e-004;
		Wl=1.0;
		Wa=1.0;
		Wr=2.0;
	  Wb=10.0;
      Do i=1,NumElD
        YRL(i)=0.
        RUTDEN(i)=0.
        RGCF(i)=0.
        iFavRoot(i)=0.
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
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) Wl,Wa,Wr,Wb
       im=im+1
      il=il+1
      il=0
      Close (40)
      open(41,file=ElemGeomFile,status='old',ERR=10)
      im=im+1
      il=il+1
      Read(41,*,ERR=10)
      im=im+1
      il=il+1
      Read(41,*,ERR=10)
      im=im+1
      Do e=1,NumEl
        il=il+1
        Read(41,*,ERR=10) n,RTWT(n),(AdjCellNum(n,l),l=1,4)
  !      if ((n.eq.49).or.(n.eq.61)) then
  !      RTWT(n)=0.013
  !       else 
         RTWT(n)=0.
  !      end if
       
      Enddo
      Close(41)
*
      EORSCS=0.0
      AWUPSS=0.0
      Total_Eor=0.0
      Total_Pcrs=0.0
      TotalRootWt=0.0
      PSISM=0.
      icount=0
      iCellCount=0
      RootsInitiated=0
      xMin=1.E+32
      xMax=-1.E+32
      InitialRootCellArea=0.0
      
      Do e = 1,NumEl
        i=KX(e,1)
        j=KX(e,2)
        k=KX(e,3)
        l=KX(e,4)
        WidthE(e)=Amax1(abs(x(i)-x(j)),abs(x(i)-x(k)),abs(x(i)-x(l)))
        HeightE(e)=Amax1(abs(y(i)-y(j)),abs(y(i)-y(k)),abs(y(i)-y(l)))
        xMin=Amin1(x(i),x(j),x(k),x(l),xMin)
        xMax=Amax1(x(i),x(j),x(k),x(l),xMax)
!        YRL(e) = RTWT(e)/RTWL
        NumAng=3
        Do j=1,3
          iLoc(j)=KX(e,j)
        Enddo
        If (KX(e,4).ne.KX(e,3)) then
          iLoc(4)=KX(e,4)
          NumAng=4
        Endif
        VH2OC(e) = 0.
        PSIS(e)  = 0.
        TS(e)    = 0.
        CXT(e)   = 0.
        Conc_e(e) =0.
        xMean(e)=0.
        yMean(e)=0.
        InitialRootCell(e)=0
         Do i = 1,NumAng
          j = iLoc(i)
          VH2OC(e) = VH2OC(e) + ThNew(j)
          PSIS(e)  = PSIS(e)  + hNew(j)/1000. ! convert from cm to bars
          TS(e)    = TS(e)    + Tmpr(j)
          CXT(e)   = CXT(e)   + g(j,2)
          Conc_e(e) =Conc_e(e)+ Conc(j,1)
          xMean(e) = xMean(e) + (x(j)-xBSTEM)
          yMean(e) = yMean(e) + (y(j)-yBSTEM)
          
                 
        Enddo
	    VH2OC(e)  = VH2OC(e)/NumAng
	    PSIS(e)   = PSIS(e)/NumAng
	    TS(e)     = TS(e)/NumAng
	    CXT(e)    = CXT(e)/NumAng
	    Conc_e(e) =Conc_e(e)/NumAng
	    xMean(e) = xMean(e)/NumAng
	    yMean(e) = yMean(e)/NumAng
        RVR(e) = RVRL*SQRT(xMean(e)**2+yMean(e)**2)/4
        if ((abs(xMean(e)).lt.2).and.(yMean(e)>-15.0).and.
     !          (yMean(e)<0)) then
           InitialRootCell(e)=1
           iCellCount=iCellCount+1
           InitialRootCellArea=InitialRootCellArea+Area(e)
             endif 
        If(RTWT(e).gt.0.) then
          PSISM=PSISM+PSIS(e)
          icount=icount+1
        Endif
        MECHR = 80.03*abs(psis(e))**0.53
     &     *EXP(-8.08*(1.66 - BlkDn(MatNumE(e))))
        IF (MECHR.LE.0.0) MECHR = 0.0
        TPRD(e) = 2.0 + (0.1*PSIS(e) + (0.7*MECHR))
        PSIRD(e) = PSIS(e)
      Enddo

c   EOMult is defined in Crop since it is a Plant variable
c      EOMult=1.
c      If(xBSTEM.le.xMin.or.xBSTEM.ge.xMax) EOMult=0.5

c
      PSISM=PSISM/float(icount)
c      PPSIL = -2.0
c      PPSILT = -2.0
      PPSIL = -0.5
      PPSILT = -0.5
      TPLT = 2.0
c      PSILD = -2.0
      PSILD = -0.5
c since in the MAIZESIM model, the reference leaf water potential in the CGas_exchange::setleafpEffect
c function is set to be -1.9, so the initial value of leaf water potential should be assigned a higher
c value Yang 8/22/06
      IF (PSISM.LT.PSILD ) PSILD = PSISM
      PILOSM = -9.00
      OSMREG = -0.25
      OSMFAC = 0.5
      PILD = ((PSILD + 2.0)*OSMFAC) + PILOSM
      TPLD = PSILD - PILD
      TPL = PPSIL - PILD
      BTPL = TPLD
      SGTLT = 0.0                   
11    continue
C
c  Routine calculations
C
c                     Cell variables
*     This piece of code averages nodal values over cells and
*     calculates the cell hydraulic conductivity.
*     It has no shoot-associated parameters
C
CYAP  No Plant = No Root Activity
      If((NShoot.eq.0).or.(isGerminated.eq.0)) Return
       PDWR_TOTAL=0.0    ! for debugging purposes
CYAPEND
      If (isGerminated.and.RootsInitiated.eq.0) then
         RootsInitiated=1
         InitialRootWt=InitialRootCarbo/InitialRootCellArea   ! calculates average root wt density
         ! Distribute this mass among cells with initial root distribution
         AreaS=0
          Do E=1, NumEl
             if (InitialRootCell(e).gt.0) then
               RTWT(e)=InitialRootWt*area(e)
               YRL(e) = RTWT(e)/RTWL
                TotalRootWeight=TotalRootWeight+RTWT(e)
                AreaS=AreaS+area(e)
              endif
           enddo
        Endif
        
      Do e=1,NumEl
        NumAng=3
        Do j=1,3
          iLoc(j)=KX(e,j)
        Enddo
        If(KX(e,4).ne.KX(e,3)) then
          iLoc(4)=KX(e,4)
          NumAng=4
        Endif
        VH2OC(e) = 0.
        TS(e)    = 0.
        CXT(e)   = 0.
        COND(e)  = 0.
        Conc_e(e)= 0.
        Do i = 1,NumAng
          j = iLoc(i)
          VH2OC(e) = VH2OC(e) + ThNew(j)
          PSIS(e)  = PSIS(e)  + hNew(j)
          TS(e)    = TS(e)    + Tmpr(j)
          CXT(e)   = CXT(e)   + g(j,2)
          Conc_e(e) =Conc_e(e)+ Conc(j,1)
          COND(e)  = COND(e)  + Con(j)
        Enddo
        VH2OC(e) = VH2OC(e)/NumAng
        PSIS(e)  = PSIS(e)*0.001/NumAng
        TS(e)    = TS(e)   /NumAng
        CXT(e)   = CXT(e)  /NumAng
        Conc_e(e)= Conc_e(e)/NumAng
	  COND(e) = COND(e)/NumAng
    
      Enddo
*
c     EO : potential transpiration rate
c     EO is in g/plant. Multiply g/plant x plant/m2 x 1 m2/10000 cm2 x rowsp
c     EOR : = EO*ROWSP if the plant is not at the border of the soil slab,
c           = 0.5 * EO*ROWSP, if the plant is at the left or right border.
c      EOR = EO*poprow/rowsp*100*rowsp*EOMult/10000
c    If we are talking about g/m2, maybe we should not have the extra rowsp? Yang 
c        EOR = EO*poprow/rowsp*100*EOMult
c     note that rowSP/10000 is area per slab in m2
c        
c       1/popRow is distance between two plants along the row
c        rowsp is distance from row to row. thus the area of a plant is
c        rowsp/100 * 1/poprow
c       thus poprow*100/rowSP is area per plant in m2
c       1 over area per plant is poprow*100/rowsp
       EOR = ET_demand*poprow*(100/rowsp)*(EOMult*rowsp/10000)
c       EOR = ET_demand
c      PCRL = PCRL*12*poprow/rowsp*100*rowsp*EOMult/10000
c      PCRQ = PCRQ*12*poprow/rowsp*100*rowsp*EOMult/10000
      
c      EOR = ET_demand*7*(rowsp*EOMult/10000)
     
c      PCRL = PCRL*24*7*rowsp*EOMult/10000
c      PCRQ = PCRQ*24*7*rowsp*EOMult/10000
c   PCRL in MAIZESIM is in the unit g/plant/hour, so mulitply 24 first to convert it to g/plant/day
c   then multiply it by plant/(m^2) to convert it to g/(m^2)/day
c   finally multiply it by rowsp*1(cm) to convert it into the soil slab we are dealing with Yang 9/24/06
             
*     
C     ---------- Taken form R U T G R O   of  G L Y C I M
C
C  CALCULATE SOIL MECHANICAL RESISTANCE TO ROOT GROWTH IN EACH CELL
C
      DO e=1,NumEl
        MECHR = 10.03*abs(psis(e))**0.25
     &     *EXP(-10.58*(1.7 - BlkDn(MatNumE(e))))
        IF (MECHR.LE.0.0) MECHR = 0.0
C
C  CALCULATE ROOT TURGOR PRESSURE IN EACH CELL AT DAWN MINUS THRESHOLD
C  TURGOR FOR GROWTH
C
        IF (ITIME.EQ.IDAWN) THEN
          TPRD(e) = 2.0 + (0.1*PSIS(e) + (0.7*MECHR))
          PSIRD(e) = PSIS(e)
        ENDIF
C
C  CALCULATE PROPORTIONIAL REDUCTION OF ROOT GROWTH CAUSED BY
C  MECHANICAL RESISTANCE AND SOIL WATER POTENTIAL ASSUMING ROOT
C  IS AT SAME POTENTIAL AS SOIL
C
C        RGCF1=AMin1(1.,2.*1.66/BlkDn(MatNumE(e))*(-0.3)/psis(e))
          RGCF1 = ((TPRD(e) - MECHR)/2.0)
     &     - ((PSIRD(e) - PSIS(e))/4.0)
C
C     CALCULATE PROPORTIONAL REDUCTION OF ROOT GROWTH CAUSED BY
C     SOIL TEMPERATURE
C
        RGCF2=1.0
        IF (TS(e).GE.33.0) THEN
          RGCF2 = (TS(e)/33.)**(-1.66)
        ELSE IF (TS(e).LE.0.0) THEN
          RGCF2 = 0.0
        ELSE If (TS(e).lt.20.0) THEN
          RGCF2 = (TS(e)/20.) **(1.66)
        ENDIF


C
C   CALCULATE PROPORTIONIAL REDUCTION OF ROOT GROWTH CAUSED BY
C   PARTIAL PRESSURE OF SOIL OXYGEN
C
        RGCF3 = (CXT(e) - 0.02)*7.14
C
C   CALCULATE PROPORTIONAL REDUCTION OF ROOT GROWTH FROM THE
C   PHYSICAL CAUSE THAT IS MOST LIMITING
C
cdt test
c        RGCF2=1
c        RGCF4=Exp(psis(e)*0.6)
c        if (abs(psis(e)).lt.0.45) RGCF4=1.0
c        RGCF4=1.0
C limit here is 0.03 gr cm-3 dt changed to 0.01 for testing
        RGCF5=(1.0-AMIN1(1.0,RTWT(e)/Area(e)/0.03))
C this one for nitrogen. Use the ratio between actual uptake of N and IMax
    
         dif=amax1(0.0,Conc_e(e)-CMin0(1))
	   RGCF6 =  amax1(0.3,dif/ConstK(1))
	   RGCF6 = AMIN1(1.0, RGCF6)
        RGCF(e) = AMIN1(RGCF1, RGCF2, RGCF3,RGCF5)

c DT for debugging only
        !RGCF(e)=AMIN1(RGCF6,RGCF5)
        
        RGCF(e) = AMAX1(RGCF(e),0.0)
        RGCF(e) = AMIN1(RGCF(e),1.0)

c 
C if seed has just germinated and not emerged, then the RGCF is 1
C      
        if ((isEmerged.lt.1).and.(isGerminated.gt.0)) RGCF(e)=1.0
C
C
C   CALCULATE POTENTIAL RATE OF CHANGE IN ROOT DRY WEIGHT IN EACH
C   SOIL CELL

c
c rtwt is gr in the cell
c alpm a Percentage per day
c RGCF has no units
c pdwr = gr d-1 per cell
cdt 03/18/11 testing the incorporation of area
           PDWR(e) = RTWT(e)*ALPM*RGCF(e)*area(e) ! but the bigger the area, the smaller is pdwr - this may not work to divide....
           PDWR_TOTAL=PDWR_TOTAL+PDWR(e)
      ENDDO
c      do e=1, NumEL
c       write(108,25) time,e, rgcf1,rgcf2,rgcf3
c      enddo
C
CDT GLYCIM sorts on PDWR which is potential for root growth. RGCF is already
C   taken into account.
      Call SORT(NumEl,iFavRoot,RGCF,RVR)
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
      DO e = 1,NumEl
        SR(e)   = 0.0
        WUPM(e) = 0.0
        WUPN(e) = 0.0
        WUPT(e) = 0.0
        DRL(e)  = 0.0
cS 07/05/2007 DT These variables are needed for nitrogen sink code    
        FUP(E,1) = 0.0 !WATER FLUX TO THE OLD ROOTS or old roots density 
	  FUP(E,2) = 0.0 !WATER FLUX TO THE NEW ROOTS or new roots density       
        VUP(e,1) = 0.0 !WATER FLOW TO THE OLD ROOTS
	  VUP(e,2) = 0.0 !WATER FLOW TO THE NEW ROOTS
cE 7/05/2007 DT 
      ENDDO
      SGTLI = 0.0
C
C   THE THRESHOLD LEAF WATER POTENTIAL AND TURGOR PRESSURE THAT
C   WILL JUST PREVENT ALL SHOOT EXPANSION ARE CALCULATED
C
C   CALCULATE DAWN LEAF WATER POTENTIALS ***
C
      IF (ITIME.EQ.IDAWN) THEN
c        PSILD = -2.0
c since in the MAIZESIM model, the reference leaf water potential in the CGas_exchange::setleafpEffect
c function is set to be -1.9, so the initial value of leaf water potential should be assigned a higher
c value Yang 8/22/06
        PSILD = - 0.5
        IF (PSISM.LT.PSILD ) PSILD = PSISM
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
      D12 = (D11 - 2.0)*(1.0 - EXP(-0.7*Step))
      IF (D12.LE.0.0) D12 = 0.0
      D13 = MAX (PPSIL, PPSILT)
      OSMFAC = 0.5
      PSILT = D13 - (D12/(1.0 - OSMFAC))
      PPSILT = PSILT
      TPLT = D11 - D12
      PSILZ = PSILT - (TPLT/(1.0 - OSMFAC))
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
      DO 50, e=1,NumEl
        
        D14 = PSIS(e) - PSILT
        IF (D14.GT.0.0) THEN
          SC = 0.017 - (PSIS(e)*0.5)
C
C   CALCULATE RESISTANCE TO WATER FLOW IN SOIL CELL
C
          IF (COND(e).GT.0.0) THEN
                  SR(e) = ALOG(SC*SC/3.0E-4)/(4.0*3.1416*COND(e)*1019.7)
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
             
             IF(LAI.EQ.0.0) Then
              WUPM(e) = (((PSIS(e) - PSILT)*Area(e)
     &         *RUTDEN(e))/(SR(e) + RRRM + RVR(e)))
     &         *0.02275*(24.5 + TS(e))
              WUPN(e) = (((PSIS(e) - PSILT)*YRL(e))
     &         /(SR(e) + RRRY + RVR(e)))
     &         *0.02275*(24.5 + TS(e))
             ELSE
             WUPM(e)=(PSIS(e)-PSILT)/((SR(e)+RRRM+RVR(e))
     &       /(RUTDEN(e)*Area(e)*0.02275*(24.5+TS(e)))
     &       +(19.26/(LAI*ROWSP)))
    
             WUPN(e)=(PSIS(e)-PSILT)/((SR(e)+RRRY+RVR(e))/
     &      (YRL(e)*0.02275*(24.5+TS(e)))+19.26/(LAI*ROWSP))
            END IF
c 
c            WUPN(e)=(PSIS(e)-PSILT)/((SR(e)+RRRY+RVR(e))
c     &     /(YRL(e)*0.02275*(24.5+TS(e)))
c     &     +(9.63/(LAI*ROWSP)))
c            

C NOTE 10.17
              WUPMS = WUPMS + (WUPM(e) + WUPN(e))
              WUPSI = WUPSI + (PSIS(e)*(WUPM(e) + WUPN(e)))
          ENDIF
        ENDIF
 50   CONTINUE
C
C  CALCULATE GROWTH BY NEW ROOTS IF ROOTS GREW AT POTENTIAL
C  RATE DURING PAST PERIOD. START WITH SOIL CELL WHERE CONDITIONS
C  ARE MOST FAVORABLE FOR GROWTH
C
      PCRTS = 0.
      WUPTS = 0.
      SW10 = 0
      DO 90, e=1,NumEl
        e1=iFavRoot(e)
        PDRL = PDWR(e1)/RTWL
        DRL(e1) = (PDRL + PPDRL(e1))*Step/2.0
        PPDRL(e1) = PDRL
C
C   DETERMINE TOTAL WATER UPTAKE BY ALL NEW ROOTS AND WEIGHTED
C   AVERAGE SOIL WATER POTENTIAL WHEN:
C   I. ROOTS GET ALL THE CARBON TRANSLOCATED
C      
        IF(LAI.EQ.0.0) Then
        WUPT(e1) = (((PSIS(e1) - PSILT)*DRL(e1))
     &   /(SR(e1) + RRRY + RVR(e1)))
     &   *0.02275*(24.5 + TS(e1))
        Else
        WUPT(e1)=(PSIS(e1)-PSILT)/((SR(e1)+RRRY+RVR(e1))
     &     /(DRL(e1)*0.02275*(24.5+TS(e1)))
     &     +(19.26/(LAI*ROWSP)))
       End IF

        IF (WUPT(e1).LE.0.0) WUPT(e1) = 0.0
        WUPTS = WUPTS + WUPT(e1)
        WUPSI = WUPSI + (PSIS(e1)*WUPT(e1))
        PCRTS = PCRTS + (PDWR(e1)*CONVR)
C
C   II. ROOTS GET ONLY THE CARBON LEFT OVER FROM SHOOT GROWTH
C
c if carbon allocation changes, have to do it here too
        IF (PCRL.LE.PCRTS.AND.SW10.NE.1) THEN
          SW10 = 1
          D16 = 0.0
          IF (PDWR(e1).GT.0.) D16 = WUPT(e1)*(1.0
     &     - ((PDWR(e1)*CONVR) - PCRTS + PCRL)/(PDWR(e1)*CONVR))
          WUPRS = WUPTS - D16
C NOTE 10.23
          IF ((WUPMS + WUPRS).GT.0.) PSISM = (WUPSI - (PSIS(e1)*D16))
     &     /(WUPMS + WUPRS)
        ENDIF
        IF (PCRQ.LE.PCRTS) GOTO 1111
 90   CONTINUE
C
      IF ((WUPMS + WUPTS).GT.0.) PSIST = WUPSI/(WUPMS + WUPTS)
      IF (SW10.EQ.0) THEN
        WUPRS = WUPTS
        PSISM = PSIST
      ENDIF
      GOTO 1112
C
 1111 D17 = 0.0
      IF (PDWR(e1).GT.0.0) D17 = WUPT(e1)*(1.0 - ((PDWR(e1)*CONVR)
     & - PCRTS + PCRQ)/(PDWR(e1)*CONVR))
      WUPTS = WUPTS - D17
C NOTE 10.24
      IF ((WUPMS + WUPTS).GT.0.) PSIST = (WUPSI - (PSIS(e1)*D17))
     & /(WUPMS + WUPTS)
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
C
C   CALCULATE POTENTIAL TRANSPIRATION RATE ALLOWING FOR NIGHT-TIME
C   STOMATAL CLOSURE.  THE STOMATA ARE ASSUMED TO CLOSE TIGHTER
C   WHEN THE PLANT IS EXPERIENCING WATER STRESS
C

**************************************************************
c Stomatal openning and closing is determined in MAIZESIM, 
c so no need to for 2dsoil to determine closure/openning of
c stomatal Yang 8/24/06

c      IF (ITIME.LT.IDAWN.OR.ITIME.GT.IDUSK) THEN
c        SCF = 0.1
c        IF (PILOSM.LT.-9.0) SCF = 0.01
c        EOR = EOR*SCF
c      ENDIF
*****************************************************************
C
C   DETERMINE HOW POTENTIAL TRANSPIRATION RATE IS RELATED TO THE
C   VARIOUS POTENTIAL WATER UPTAKE RATES AND GO TO THE APPROPRIATE
C   SECTION OF THE CODE. CALCULATE LEAF WATER POTENTIAL,STOMATAL
C   CONDUCTANCE,THE PROPORTION OF TIME FOR WHICH THE SHOOT GROWS,
C   THE AMOUNT OF CARBON ACTUALLY NEEDED FOR ROOT GROWTH AND THE
C   AMOUNT OF OSMOREGULATION THAT WILL OCCUR
      D18 = WUPMS + WUPTS
C
C   POTENTIAL TRANSPIRATION RATE IS SUCH THAT THE DEMAND FOR WATER
C   UPTAKE CAN BE MET WITHOUT SHOOT WATER POTENTIAL FALLING TO THE
C   POINT WHERE THE SHOOT LOSES TURGOR AND STOPS GROWING
C
      IF (EOR.LT.WUPDS) THEN
        PSIL_ = PSISM - ((PSISM - PSILT)*EOR/(WUPMS + WUPRS))
c        SCF = 1.0
         SGT = 1.0
C
C   POTENTIAL TRANSPIRATION RATE IS SUCH THAT THE DEMAND FOR WATER
C   UPTAKE CANNOT BE MET WITHOUT SHOOT WATER POTENTIAL FALLING TO
C   THE POINT WHERE THE SHOOT LOSES TURGOR AND STOPS GROWING FOR
C   PART OF THE PERIOD.  WHEN THIS HAPPENS THE ROOTS GROW
C
      ELSE IF (EOR.LT.D18) THEN
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
 
c      IF (ITIME.LT.IDAWN.OR.ITIME.GT.IDUSK) THEN
c        SCF = 0.1
c        IF (PILOSM.LT.-9.0) SCF = 0.01
c        EOR = EOR/SCF
c      ENDIF
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
