cdt 7/9/07 added code to calculate loss of N in flux out of the domain 
Cdt N in 2DSOIL is nitrate (NO3) and the units are ug per cm3 or milligram/liter
       Subroutine Nitrogen_Mass_Balance()
        include 'public.ins'
        include 'nitvar.ins'
        include 'puplant.ins'
        Include 'PuSurface.ins'
        Dimension Bi(3),Ci(3)
        Character*10 Date
        Real AE,Profile_N,Manure_N,Litter_N,
     &    Min_N,Ammon_N,Org_N,Denitr,
     &    NO3_Mean,NhumusMean,NLitterMean,NManureMean,NAmmoniaMean,All_N
        Integer ModNum
C variables to hold mulch C and N totals
        Real t_CARB, t_CELL,t_LIGN,t_CARB_N,t_CELL_N, t_LIGN_N,
     &    t_ResidueMulch, t_ResidueMulch_N, t_ResidueMulch_CNR,
     &    t_Mulch_Thick,
     &    humusC, litterC, manureC, 
     &    All_C, CHumusMean,CLitterMean, CManureMean
     
        common /N_BAL/ModNum,CFlux,CFluxPrevious,C_RespirationOM,
     !    C_RespirationRoot

      If (lInput.eq.1) then
        open(91,file=MassBalanceFileOut,status='unknown',recl=520)
        write(91,8) 'Date_time,','Date,','Min_N,','Humus_N,',
     !     'Humus_C,','Manure_N,','Manure_C,','Litter_N,',
     !     'Litter_C,','OM_CO2_C,','Root_CO2_C,','Ammon_N,','All_N,',
     !     'All_C,','water,','Denitr,','CFlux,',
     !     'Mul_CARB_N,','Mul_CELL_N,','Mul_LIGN_N,',
     !     'Tot_res_N,',
     !     'Mul_Thick,', 'Mul_Mass,', 'Mul_CNR,'


8     Format (1x,A10,T20,A5,T40,A6,T60,A8,
     &  	  T80,A8,T104,A9,T124,A9,T144,A9,
     &        T164,A9,T186,A9,T209,A11,T229,A8,T249,A6, 
     &        T271,A6,T293,A7,T313,A7, 
     &        T333,A6,T353,A11,T376,A11, 
     &        T396,A11,
     &        T418,A10,T438,A10,T458,A9,T482,A8)	 
 
      NumMod=NumMod+1
	ModNum=NumMod
	tNext(ModNum)=time
	CFlux=0
	CFluxPrevious=0

      t_CARB=0.0D0
      t_CELL=0.0D0
      t_LIGN=0.0D0
      t_CARB_N=0.0D0
      t_CELL_N=0.0D0
      t_LIGN_N=0.0D0
      t_ResidueMulch=0.0D0
      t_ResidueMulch_N=0.0D0
      t_ResidueMulch_CNR=0.0D0
      t_Mulch_Thick=0.0D0
	C_RespirationOM=0.0D0
	C_RespirationRoot=0.0D0
       Endif
	     do i=1,NumBp
	     n=KXB(i)
C Cflux is loss of N in mg
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))) then
               
cccz	        Cflux=Cflux+Q(n)*conc(n,1)*step*14/62
cccz using a starndard conversion factor "NO3mass_2_Nmass"
cccz why "mg"??????????????????????????????????????????????????????
cDT only consider flow out.
              Cflux=Cflux+min(0.0,Q(n))*conc(n,1)*step*NO3mass_2_Nmass 
           endif
C for case of downward drainage and a constrant BC 
C (does not account for upward flow yet)
           if ((CodeW(n).eq.(1))) then
             if (Q(n).lt.0) then
cccz                 Cflux=Cflux+Q(n)*conc(n,1)*step*14/62 ! only consider chemical out- flux>0
cccz using a starndard conversion factor "NO3mass_2_Nmass"
cccz why "mg"??????????????????????????????????????????????????????
                 Cflux=Cflux+Q(n)*conc(n,1)*step*NO3mass_2_Nmass  ! only consider chemical out- flux>0
             endif
           endif
          EndDo
     
csun Calculate the co2 the units are ug
	   do i=1,NumNP
C gsink_OM and gSink_root unit is ug CO2 cm-3 air, need to calculate all C_RespirationOM (ugC/domain)		   
		   C_RespirationOM=C_RespirationOM+gsink_OM(i,1)*Step*12.0*nodeArea(i)
     !       *soilair(i)/(44.0)
             C_RespirationRoot=C_RespirationRoot+gSink_root(i,1)*
     !       Step*12.0*nodeArea(i)*soilair(i)/(44.0)
	   enddo
	  
        t=time
        if (Abs(time-tNext(ModNum)).lt.0.001*Step.or.lInput.ne.0) then
           tNext(ModNum)=tNext(ModNum)+1.0
           Profile_N=0.0
           Min_N=0.0
           Org_N=0.0
           Litter_N=0.0
           Manure_N=0.0
           Ammon_N=0.0
           W_Sum=0.0
           Denitr=0.0
           humusC=0.0
           litterC=0.0
           manureC=0.0

	   Sum=0.
           W_Sum=0.
	   Do n=1,NumEl
             NUS=4
             if(KX(n,3).eq.KX(n,4)) NUS=3
             Sum1=0.
             NDenitrifyMean=0.
             NO3_Mean=0.0
             NhumusMean=0.0
             NLitterMean=0.0
             NManureMean=0.0
             NAmmoniaMean=0.0
             CHumusMean=0.0
             CLitterMean=0.0
             CManureMean=0.0

*         Loop on subelements
             do k=1,NUS-2
               i=KX(n,1)
               j=KX(n,k+1)
               l=KX(n,k+2)
               Ci(1)=x(l)-x(j)
               Ci(2)=x(i)-x(l)
               Ci(3)=x(j)-x(i)
               Bi(1)=y(j)-y(l)
               Bi(2)=y(l)-y(i)
               Bi(3)=y(i)-y(j)
               AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
cccz Conc is NO3 as ug/cm3
               NO3_Mean=NO3_Mean
     &             +AE*(Conc(i,1)*ThNew(i)
     &                 +Conc(j,1)*ThNew(j)
     &                 +Conc(l,1)*ThNew(l))
     &                 /3.*NO3mass_2_Nmass  ! convert to N
               
               Sum1=Sum1+AE*(ThNew(i)+ThNew(j)+ThNew(l))/3.

cccz need to present all the N mass based using soil bulk density   
c  organic N is in ug/cm3 of volume, don't need BD order to sum over volume
               NhumusMean=NhumusMean+AE*(Nh(i)
     &             +Nh(j)
     &             +Nh(l))/3.
               NLitterMean=NLitterMean+AE*(NL(i)
     &             +NL(j)
     &             +NL(l))/3.
               NManureMean=NManureMean+AE*(Nm(i)
     &             +Nm(j)
     &             +Nm(l))/3.
     
               ChumusMean=ChumusMean+AE*(Ch(i)
     &             +Ch(j)
     &             +Ch(l))/3.
               CLitterMean=CLitterMean+AE*(CL(i)
     &             +CL(j)
     &             +CL(l))/3.
               CManureMean=CManureMean+AE*(Cm(i)
     &             +Cm(j)
     &             +Cm(l))/3.

               NAmmoniaMean=NAmmoniaMean+AE*(NH4(i)
     &             +NH4(j)
     &             +NH4(l))/3.
               NDenitrifyMean=NDenitrifyMean+AE*(Denit(i)
     &             +denit(j)
     &             +denit(l))/3.
       
             Enddo
             Min_N=Min_N+NO3_Mean
             Org_N=Org_N+NhumusMean
             Litter_N=Litter_N+NLitterMean
             Manure_N=Manure_N+NManureMean
             humusC=humusC + CHumusMean
             litterC=litterC + CLitterMean
             manureC=manureC + CManureMean
             Ammon_N=Ammon_N+NAmmoniaMean
             W_Sum=W_Sum+Sum1
             Denitr=Denitr+NDenitrifyMean
         Enddo
                    	   
Csun Carbon mass balance in soil means the soilC is consistent during time
	   All_C = humusC + litterC + manureC + C_RespirationOM
	   
C calculate total N in surface mulch residue. units are grams
C in the mulch model so multiply by 1e6
        t_CARB=0.0D0
        t_CELL=0.0D0
        t_LIGN=0.0D0
        t_CARB_N=0.0D0
        t_CELL_N=0.0D0
        t_LIGN_N=0.0D0
        t_ResidueMulch=0.0D0
        t_ResidueMulch_N=0.0D0
        t_ResidueMulch_CNR=0.0D0
        t_Mulch_Thick=0.0D0	  

        if(BoolMulchApply.eq.1) then
cccz make sure there is mulch applied
          do kk=1,SurNodeIndex_M-1
cccz should be SurNodeIndexM because there may horizontal merging
            do jj=1,mulchLayer
              t_CARB=t_CARB+CARB_mass(jj,kk)*1e6
              t_CELL=t_CELL+CELL_mass(jj,kk)*1e6
              t_LIGN=t_LIGN+LIGN_mass(jj,kk)*1e6
              t_CARB_N=t_CARB_N+CARB_N_mass(jj,kk)*1e6
              t_CELL_N=t_CELL_N+CELL_N_mass(jj,kk)*1e6
              t_LIGN_N=t_LIGN_N+LIGN_N_mass(jj,kk)*1e6
            enddo
          enddo
          t_ResidueMulch=t_CARB+t_CELL+t_LIGN
          t_ResidueMulch_N=t_LIGN_N+t_CARB_N+t_CELL_N
          t_ResidueMulch_CNR=t_ResidueMulch*0.41D0/t_ResidueMulch_N
          t_Mulch_Thick=mulchThick
        endif  

C calculate sum of N in all forms in soil
C factor is appropriate for mg/cm3 Total N is mg per slab (grid width x 1cm)
C Mineral ad organic N is ug/cm3. Do a summation over the simulation domain - total ug in a plant slab
C Fact now works while there is a plant 

          fact=1.0/(0.01*RowSp/100.*EOMult) !m2 of slab
          fact=fact*10000./1000./1000./1000.    !m2 of slab ->ha, ug-->Kg
          All_N=(Min_N+Org_N+Litter_N+Manure_N
     !           +Ammon_N)*fact
         iday=int(t)
	   call caldat(iday,mm,id,iyyy) 
          write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          write(91,10) 
     !       time,date,Min_N*fact,
     !       Org_N*fact,humusC*fact,Manure_N*fact,
     !       manureC*fact, Litter_N*fact, 
     !       litterC*fact, C_RespirationOM*fact,
     !       C_RespirationRoot*fact,Ammon_N*fact,
     !       All_N,All_C*fact,W_Sum,Denitr*fact,CFLux*fact,
     !       t_CARB_N*fact,t_CELL_N*fact,t_LIGN_N*fact,
     !        t_ResidueMulch_N*fact,
     !       t_Mulch_Thick,t_ResidueMulch*fact,t_ResidueMulch_CNR 
     	   
          CFluxPrevious=CFlux
csun reset a initial value of C_RespirationOM 
	    C_RespirationOM = 0.0D0
          C_RespirationRoot = 0.0D0
        endif
10    Format (1F12.4,',', A12,',', 21(F20.4, ','),F20.4)
      return
      end


