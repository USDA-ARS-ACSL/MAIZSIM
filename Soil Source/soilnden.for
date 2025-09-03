      Subroutine SoilNitrogen()
C new comment
c second new comment
      Include 'public.ins'
      Include 'nitvar.ins'
      Include 'PuSurface.ins'
      common /nitrog/ModNum, ThOld(NumNPD)
      Real*4 BCh,BNh, BCl,BNl,BCm,BNm,BNH4,BNO3,BDENIT
      Real*4 BNO3_Added
      REAL*4 P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, 
     &   P12, P13, P14, P15
      REAL*4 Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11,
     &       Q12,Q13, Q14, Q15 
      REAL MolO2Available, MolCO2_OM

      REAL total  ! temp variable to hold total mass of organic components in a node 
      Real*8 wfPore,O2_mol  !water filled pore space, moles of O2

	 
      Logical enough
       t=time
       If(lInput.eq.1) then
        im=400
        il=0
        Open(40,File=NitrogenFile,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
         il=il+1
        Read(40,*,ERR=10)
        im=im+1
        Read(40,*) RowSp
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        Read(40,*,ERR=10)
        Do i=1,NMat
           il=il+1
           
           Read(40,*,ERR=10) m,kh0(m),kL0(m),km0(m),kn0(m),kd0(m),
     &      fe(m),fh(m),r0(m),rL(m),rm(m),fa(m),nq(m),cs(m)
     
	  Enddo    

        Close(40)
        NumSol=1
        !Call SetAbio(ew,et,ed,0,0.,0.)
        Call SetAbio_O2(ew,eO2,et,ed,0,0.,0.,0.D0,0.D0)    
        dtmx(4)=1./24
       Else  ! not the first time step
C
C               Soil water concentration is ug NO3 cm-3 of water  (multiply by water to get ug N per cm3 of soil) 
C               The units of N in this routine are mg N liter-1 of soil. 
C               The two are equivalent (both numerator and denominator differ by a factor of 1000)
C        
C           Conc is input as ppm, then converted to ug/cm3 in solmov initialization


        Do i=1,NumNP
        cSink_OM(i,1)=0.0
        m=MatNumN(i)
        BCh=Ch_Old(i)
        BCm=Cm_Old(i)
        BNm=Nm_Old(i)
        BCL=CL(i)
        BNL=NL(i)     
        BNh=Nh(i)      
        BNH4=NH4(i) 
        NH4_old(I)=NH4(i)
        BDENIT=Denit_old(i)
cccz        conc is ug NO3 per cm3 of water  so convert to ug/cm3 soil (mg L-1 of soil)
        NO3_OLD(i)=Conc(i,1)*thNew(i)*NO3mass_2_Nmass
        Nh_old(i)=Nh(i)    !Nh was modified in the mulch decomp model
        Ch_old(i)=Ch(i)
        NL_old(i)=NL(i)
        CL_old(i)=CL(i)
        Nm_old(i)=Nm(i)
        Cm_old(i)=Cm(i)
        
C

cccz note here the unit of NH4_Old is already changed
CEH:      Get current water-filled pore space in percent 
         !thSat: total porosity equal to the saturated water content       
        wfPore = AMin1((ThNew(i)+ThOld(i))/2./thSat(m), 1.0)  !water filled pore space 
        DO IT=1,3
CEH  Rate constants: Note 1)ew and eO2 are to account for dry and wet soil water condition respectively. 
          ! 2)current equation for ew is different from the old one in SetAbio
          ! 3)Use SOMMassRatio(i) to get O2 available for OM decomposition
          !   => g(n,2)*SOMMassRatio(i)*ugGasCm3air_to_ppm(2)/10000.  !10000 to convert ppm to percent
          O2_mol = g(i,2)/32.0*soilair(i)*SOMMassRatio(i)
          Call SetAbio_O2(ew,eO2,et,ed,    
     &       m,(ThNew(i)+ThOld(i))/2.,(Tmpr(i)+TmprOld(i))/2.
     &       ,O2_mol, wfPore*100) !*100 to convert to percent
          
!C  Rate constants

C  Dehumification rate (day-1)
          kh=kh0(m)*ew*et*eO2
C  Plant residue decomposition rate (day-1)
          kL=kL0(m)*ew*et*eO2
C  Organic fertilizer decomposition rate (day-1)
          km=km0(m)*ew*et*eO2
C  Nitrification rate (day-1)
          kn=kn0(m)*ew*et*eO2
C  Denitrification rate (day-1)
C calculate mean of current concentration and concentration from past time step
cccz note the unit is "ug N per g soil"
          Aux=(BNO3+NO3_Old(i))/2.
          kd=kd0(m)*ed*et*Aux/(Aux+cs(m))
Csb:kd is the potential denitrification rate
Csb:ed is the water content correction factor
Cs is the Michaelis-Menten constant. 

C  Carbon (P) and nitrogen (Q) fluxes
          P1 =    kh * (Ch_old(i)+BCh)/2.  ! rate is given as ug C per cm3 of area 
          Q1 =    kh * (Nh_old(i)+BNh)/2.  !N from humus pool to Nitrate pool (NH4)
          P2 =    kL * ((CL_old(i)+BCL)/2.) * fe(m) * fh(m)
          Q2 =    P2 / r0(m)               !N from litter to humus pool
          P12 =   km * ((Cm_old(i)+BCm)/2.) * fe(m) * fh(m)
          Q12 =   P12 / r0(m)              !N From organic fert to humus pool
          P3 =    kL * (CL_old(i)+BCL)/2.
          Q3 =    kL * (NL_old(i)+BNL)/2.   ! N from litter pool 
          P13 =   km * (Cm_old(i)+BCm)/2.
          Q13 =   km * (Nm_old(i)+BNm)/2.    ! N from Organic Fert pool to NH4	
C Added mineral nitrogen
          Added = (Q1+Q3+Q13)*Step          ! N going to mineral N via NH4
C Potentially immobilized and lost mineral nitrogen
          Q45pot =   kL * ((Cl_old(i)+BCL)/2.) * fe(m) / r0(m)
          Q1415pot =   km * ((Cm_old(i)+BCm)/2.) * fe(m) / r0(m)
          Q7 =  kd                           ! N lost through denitrification
          PotLost = (Q45pot + Q1415pot + Q7)*Step  !N lost from min N pool via immobilization and denitrification
C Present mineral nitrogen
CDT this uses N from previous time step
          
cccz the “present” makes sense now since those component are all based on "ug N per g soil"
          Present  = NO3_old(i)+NH4_old(i)
C 'Enough' is true when it is enough mineral nitrogen for immobilization
C Q45 is immobilization via NH4 and NO3 via the litter pool
C Q1415 is immobilization vie the organic fertilizer pool
          If(PotLost.Lt.Added+Present) then
             Enough=.true.
          else
             Enough=.false.
          endif
          If(Enough) then 
            Q45act  = Q45pot
            Q1415act= Q1415pot
          else
            Q45act  =fa(m)*Present*Q45pot  /(Q45pot+Q1415pot)
            Q1415act=fa(m)*Present*Q1415pot/(Q45pot+Q1415pot)
          endif

C
C Preferential immobilization of ammonium
C
C Nitrification
          Aux1=(NH4_Old(i)+BNH4)/2.  !8/28/2014 DT was BNH4 original typo
          Aux2=(BNO3+NO3_Old(i))/2.
          Q6 =  kn * AMAX1(0.,Aux1 - Aux2/nq(m)) ! nitrification
C Ammonium available for immobilization
          Avail = NH4_old(i) + (Q1+Q3+Q13-Q6)*Step  !all NH$ sources and sinks
C Enough is .true. if available ammonium cover all needs for immobilization
C dt 7-9-2007 Added code below to account for case when enough is false because there is almost no N in ths soil
C  then Q45act and Q1415act are zero. Have to avoid a divide by zero error.

          If((Q45act+Q1415act)*Step.LT.Avail) then
            Enough=.true.
          else
            Enough=.false.
          endif
          If(Enough) then
            Q4 = Q45act
            Q5 = 0.
            Q14= Q1415act
            Q15= 0.
          else
            Q4 = Avail    * Q45act/(Q45act+Q1415act)
             if (isnan(Q4)) Q4=0
            Q14= Avail    - Q4
            Q5=  Q45act   - Q4
            Q15= Q1415act - Q14
          endif
          P45   =  Q45act   * r0(m)
          P1415 =  Q1415act * r0(m)
C Incrementing nitrogen and carbon contents in compartments
          BCh = Ch_Old(i)+ Step * ( - P1 + P2 + P12)
          BNh = Nh_Old(i)+ Step * ( - Q1 + Q2 + Q12)
          BCL = CL_Old(i)+ Step * ( - P2 - P3 + P45)
          BNL = NL_Old(i)+ Step * ( - Q2 - Q3 + Q4 + Q5)
          BCm = Cm_Old(i)+ Step * ( - P12 - P13 + P1415)
          BNm = Nm_Old(i)+ Step * ( - Q12 - Q13 + Q14 + Q15)
          BNH4 = NH4_old(i) + Step *
     &             ( Q1 + Q3 + Q13 - Q4 - Q14 - Q6 )
          BNO3 = NO3_old(i) + Step * 
     &             ( Q6 - Q7 - Q5 - Q15)
          BNO3_Added=( Q6 - Q7 - Q5 - Q15)*step
          
          if (BNO3_Added.lt.-0) then
             jjj=1
             endif

c rate is given as ug C per cm3 of area (mg per L) then P rate transfer to C mol
C No negative values
		g(i,1) = AMAX1(g(i,1),0.)							   
          BCh = AMAX1(BCh,0.)
          BNh = AMAX1(BNh,0.)
          BCl = AMAX1(BCl,0.)
          BNl = AMAX1(BNl,0.)
          BCm = AMAX1(BCm,0.)
          BNm = AMAX1(BNm,0.)
          BNH4 = AMAX1(BNH4,0.)
          BNO3  = AMAX1(BNO3,0.)
          BDENIT=Denit_old(i)+Q7*Step
	 ENDDO


csun (P1+P3+P13):[ugC/cm3 volume]/day
csun Step:[day]
csun Soilair(i):        [cm3air/cm3 volume]
Csun  the release of C from organic matter (: ug C /cm3 soil day-1) to (ug CO2 cm-3 air)
cDT removed BlkDn, it is not needed. 
	   gSink_OM(i,1)=(P1+P3+P13-P45-P1415)*
     &         44.0/(12.0)                    ! calculate ug of CO2 produced
         gSink_OM(i,2)=-gSink_OM(i,1)/44.0*32.0      ! calculate umoles of O2 consumed - 1 mole CO2=1 mole O2
                                                 ! /44 converts from g to umoles of CO2 and *32.0 converts umoles of O2 to ug O2
                                          ! negative sign indicates O2 consumption
         !MolCO2_root=gSink_root(i,1)*step/44.0       ! umol of CO2 respired  in this node by roots
         MolCO2_OM=gsink_OM(i,1)*step/44.0       ! umol of CO2 respired  in this node by OM
         !Moles of O2 available for OM decomposition
         MolO2Available=g(i,2)/32.0*soilair(i)*SOMMassRatio(i)    ! convert ug/cm3 of O2 to total umol of O2 in the soil pore air this is the maximum amount of  O2 available for respiration in umol

         ! these next lines for debugging and testing, remove later
         if (MolO2Available.lt.MolCO2_OM) then
             iii=1
         endif

         gSink_OM(i,1)=gSink_OM(i,1)*amax1(1/soilair(i),0.0) ! convert to ug CO2 cm-3 air
         gSink_OM(i,2)=gSink_OM(i,2)*amax1(1/soilair(i),0.0) ! convert to ug O2 cm-3 air
         
         
Csb: Sink for N2O      
Csb: 1ugN=44/14 ug N2O	
Csb: unit of Q7 [ug N g soil day-1]  
Csb: unit of gsink_N2O [ug N2O cm=3 air]         
         
         gsink_N2O(i,1)=Q7*(44.0/14.0)*amax1(1.0/soilair(i),0.0)!ug N to ug N2O +convert to ug N2O cm-3 air	

  				 
C End of iterations
          Ch(i) = BCh
          Nh(i) = BNh
          CL(i) = BCL
          NL(i) = BNL
          Cm(i) = BCm
          Nm(i) = BNm
          Denit(i)=BDENIT
          NH4(i) = BNH4
          NO3_sol = BNO3

cccz the unit of NO3 and BNO3 are "ug N per cm3 soil", so need adjustment to calculate "Conc(i,1)"
CDT NO3_from_residue is in ug N per cm3 soil for the time step the sink is not adjusted for water content
CDT It has been multiplied by width in  the Mulch Decomp routine. It was also multiplied by step to get  a total
CDT amount for the time step. Here we have to divide by step to bring it back to per day units. We don't use thnew because this
CDT is not dissolved in water.
 
          cSink_OM(i,1)=-NO3_from_residue(i)/
     &        NO3mass_2_Nmass/step  
          cSink_OM(i,1) = cSink_OM(i,1) -
     &      (BNO3-NO3_old(i))/
     &       NO3mass_2_Nmass/step


          TotNitO=Nh_old(i)+Nl_old(i)+Nm_Old(i)+NO3_Old(i)+NH4_old(i)
          TotNit=Nh(i)+Nl(i)+Nm(i)+NO3_sol+NH4(i)  
          DTot=TotNit-TotNitO
          ThOld(i)=ThNew(i)
          TmprOld(i)=Tmpr(i)
        Enddo
      Endif


      call Nitrogen_Mass_Balance()
C      
C calculate proportion of roots and SOM in each node    
C
      Do i=1,NumNP
       total=( RMassY(i) + RMassM(i) ) * 1.0e6 + Ch(i) + CL(i) + Cm(i) ! roots are grams per cm3 of soil convert to ug per cm3
       RootRatioY(i)=amax1(RMassY(i)/total*1.0e6,0.0)
       RootRatioM(i)=amax1(RMassM(i)/total*1.0e6,0.0)
       SOMMassRatio(i)=1.0-RootRatioY(i)-RootRatioM(i)  ! can lump the SOM components together as they get a common ew for now
      
      EndDo
 
      

      Return
10    Call Errmes(im,il)
      End

C==================================================================
CEH: Rate correction factors for tempeature and soil water content AND O2 available
      Subroutine SetAbio_O2(ew,eO2,et,ed,m,theta,T,O2_mol,wfPore)
      Include 'public.ins'
      Common /Abio/ dThH,dThL,es,Th_m,tb,QT,dThD,Th_D
      Real*8 kMO2, kMW,O2_mol, wfPore

      !kMO2 = 33.157 ! Michaelis constant for oxygen  [%] from Sierra (2017) Biogeosciences 14.3 (2017): 703-710.
      kMO2 = 2.0   ! Michaelis constant for O2 concentration [O2 moles in soil pores]
      kMW = 76.026 ! Michaelis constant for water filled pore space  [%], theta/Th_sat*100 => mean values in Table 1
  
      If(lInput.eq.1) then
        im=420
        il=0
        Open(40,File=BiologyFile,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) dThH,dThL,es,Th_m  !moisture dependent paramters
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) tb,QT               !Temperature depenedence
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) dThD,Th_d           !Denitrification on water content
        Close(40)
        Return
      else
        eO2 = O2_mol/(kMO2 + O2_mol)    
        ew = wfPore/(kMW + wfPore)         
        eT=QT**((T-tb)/10.)               !eT and ed controls denitrification
        ThD=TUpperLimit(m)-dThD
        if(Theta.GT.ThD) then
          ed=((Theta-ThD)/(TUpperLimit(m)-ThD))**Th_d
        else
          ed=0. !denitrification rate correction factor
        endif
      endif
      return
10    Call Errmes(im,il)
      End
