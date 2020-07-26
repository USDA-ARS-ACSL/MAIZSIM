cccz ******************************************
c     the mulch decomposition module
c     linked with mulch process module
c     10/29/2019
cccz ******************************************
      
cccc  mulch decomposition
       Subroutine MulchDecomp()
       include 'public.ins'
       include 'Pusurface.ins'
       double precision Frac_RDM, Frac_HCE, Frac_CEL, Frac_LIG
       double precision PC_mulch, PN_mulch, alpha_mulch,
     !   C_rate_RDM_0, C_rate_RDM_1, tb_RDM, QT_RDM,
     !   C_rate_HCE_0, C_rate_HCE_1, tb_HCE, QT_HCE,
     !   C_rate_CEL_0, C_rate_CEL_1, tb_CEL, QT_CEL,
     !   C_rate_LIG_0, C_rate_LIG_1, tb_LIG, QT_LIG,
     !   CO2mass_2_Cmass, ThetaV_2_ThetaM,
     !   C_Decomp_RDM,C_Decomp_HCE, C_Decomp_CEL, C_Decomp_LIG, 
     !   et_RDM, et_HCE, et_CEL, et_LIG, bbbb
       integer EleIndex,EleIndexUp,EleIndexDown
       double precision mulch_mass_temp
       double precision P_RDM_C,Q_RDM_N,P_HCE_C,Q_HCE_N,
     !   P_CEL_C,Q_CEL_N,P_LIG_C,Q_LIG_N,
     !   Mulch_Rate,Mulch_Avail,Mulch_AvailUp,Mulch_AvailDown
       double precision total_decomp_RDM,total_decomp_HCE,
     !   total_decomp_CEL,total_decomp_LIG,total_Length
       double precision Frac_Coef(11)
       double precision Mulch_Rate_RDM,Mulch_Rate_HCE,
     !   Mulch_Rate_CEL,Mulch_Rate_LIG,Mulch_Decompose_RDM_C,
     !   Mulch_Decompose_HCE_C,Mulch_Decompose_CEL_C,
     !   Mulch_Decompose_LIG_C,Mulch_Decompose_RDM_N,
     !   Mulch_Decompose_HCE_N,Mulch_Decompose_CEL_N,
     !   Mulch_Decompose_LIG_N
       Character InString*132
       common /MulchDecomposition/ 
     !   PC_mulch,PN_mulch,alpha_mulch,CO2mass_2_Cmass,ThetaV_2_ThetaM,
     !   C_rate_RDM_0, C_rate_RDM_1, tb_RDM, QT_RDM,
     !   C_rate_HCE_0, C_rate_HCE_1, tb_HCE, QT_HCE,
     !   C_rate_CEL_0, C_rate_CEL_1, tb_CEL, QT_CEL,
     !   C_rate_LIG_0, C_rate_LIG_1, tb_LIG, QT_LIG,
     !   Frac_RDM, Frac_HCE, Frac_CEL, Frac_LIG
       IF(lInput.eq.1) then
        Open(21,file=MulchFile,status='old',ERR=3116)
3115        Read (21,'(A132)') InString
            if (InString(1:23).ne.'[Mulch_Decomposition_1]') goto 3115
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110) PC_mulch, PN_mulch, alpha_mulch  ! g/g for C and N, alpha for feeding (upper layer to lower layer)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110) Frac_RDM, Frac_HCE, Frac_CEL, Frac_LIG
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110) C_rate_RDM_0,C_rate_RDM_1,
     &     C_rate_HCE_0,C_rate_HCE_1,C_rate_CEL_0,C_rate_CEL_1,
     &     C_rate_LIG_0,C_rate_LIG_1
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110) b_RDM,QT_RDM,b_HCE,
     &     QT_HCE,b_CEL,QT_CEL,b_LIG,QT_LIG
        close(21)
cccz need to calculate a TUpperL_mulch/TLowerL_mulch
        if(Van_Camp_mulch.eq.0) then
          bbbb=WQ_VANG_MULCH(-100.0D0,ParMulch,-1)
        elseif(Van_Camp_mulch.eq.1) then
          bbbb=WQ_CAMP_MULCH(-100.0D0,ParMulch,-1)
        endif
        CO2mass_2_Cmass=12.0107D0/44.01D0
        ThetaV_2_ThetaM=1.0D0/(1.0D0-f_mulch_pore)/rho_mulch
cccz initialize the mass of each portion (based on carbon)
        do jj=1,mulchLayer
         do kk=1,SurNodeIndex-1
          EleIndex=MulchEleMatrix(jj,kk)  
          mulch_mass_temp=MulchEleMarkerArea(EleIndex)
     &     *(1-f_mulch_pore)*rho_mulch
            RDM_mass(jj,kk)=mulch_mass_temp*PC_mulch*Frac_RDM
            HCE_mass(jj,kk)=mulch_mass_temp*PC_mulch*Frac_HCE
            CEL_mass(jj,kk)=mulch_mass_temp*PC_mulch*Frac_CEL
            LIG_mass(jj,kk)=mulch_mass_temp*PC_mulch*Frac_LIG
          enddo
        enddo
       ELSE
        total_decomp_RDM=0.0D0
        total_decomp_HCE=0.0D0
        total_decomp_CEL=0.0D0
        total_decomp_LIG=0.0D0
        total_Length=0.0D0
cccz calculate decomposition fraction based on
        bbbb=(1.0D0-alpha_mulch**mulchLayer)/(1.0D0-alpha_mulch)
        do n=1,mulchLayer
c         Frac_Coef(n)=alpha_mulch**(n-1)/bbbb
         Frac_Coef(n)=alpha_mulch**(n-1)
         Frac_Decomp(n)=0.0D0
        enddo

        do n=1,SurNodeIndex-1
         Theta_m=MulchEleThnew(1,n)*ThetaV_2_ThetaM
         Temp_m=MulchEleTmpr(1,n)
         EleIndex=MulchEleMatrix(1,n)
cccz decomposition dynamics ***************************************
         et_RDM=exp((Temp_m-tb_RDM)*QT_RDM)
         et_HCE=exp((Temp_m-tb_HCE)*QT_HCE)
         et_CEL=exp((Temp_m-tb_CEL)*QT_CEL)
         et_LIG=exp((Temp_m-tb_LIG)*QT_LIG)
         C_Decomp_RDM=(C_rate_RDM_0+C_rate_RDM_1*Theta_m)
     &     *CO2mass_2_Cmass*1.0D-4*24.0D0                     ! the decomposition of C g/cm^2/day
     &     *et_RDM
         C_Decomp_HCE=(C_rate_HCE_0+C_rate_HCE_1*Theta_m)
     &     *CO2mass_2_Cmass*1.0D-4*24.0D0 
     &     *et_HCE
         C_Decomp_CEL=(C_rate_CEL_0+C_rate_CEL_1*Theta_m)
     &     *CO2mass_2_Cmass*1.0D-4*24.0D0 
     &     *et_CEL
         C_Decomp_LIG=(C_rate_LIG_0+C_rate_LIG_1*Theta_m)
     &     *CO2mass_2_Cmass*1.0D-4*24.0D0 
     &     *et_LIG
cccz decomposition dynamics ***************************************
         P_RDM_C=C_Decomp_RDM
         P_HCE_C=C_Decomp_HCE
         P_CEL_C=C_Decomp_CEL
         P_LIG_C=C_Decomp_LIG
         Q_RDM_N=P_RDM_C/PC_mulch*PN_mulch
         Q_HCE_N=P_HCE_C/PC_mulch*PN_mulch
         Q_CEL_N=P_CEL_C/PC_mulch*PN_mulch
         Q_LIG_N=P_LIG_C/PC_mulch*PN_mulch          
cccz calculate the total loss amount (only the carbon fraction)
cccz only assume in the lowest layer, meybe the C decompose > C store
cccz but all the C decomposed should be from this layer, thus, I 
         Mulch_Rate_RDM
     &      =min(P_RDM_C*Step*widthPerMulchUnit(n),RDM_mass(1,n))   ! g
         Mulch_Rate_HCE
     &      =min(P_HCE_C*Step*widthPerMulchUnit(n),HCE_mass(1,n))
         Mulch_Rate_CEL
     &      =min(P_CEL_C*Step*widthPerMulchUnit(n),CEL_mass(1,n))
         Mulch_Rate_LIG
     &      =min(P_LIG_C*Step*widthPerMulchUnit(n),LIG_mass(1,n))        
         total_decomp_RDM=total_decomp_RDM
     &    +Mulch_Rate_RDM*widthPerMulchUnit(n)
         total_decomp_HCE=total_decomp_HCE
     &    +Mulch_Rate_HCE*widthPerMulchUnit(n)
         total_decomp_CEL=total_decomp_CEL
     &    +Mulch_Rate_CEL*widthPerMulchUnit(n)
         total_decomp_LIG=total_decomp_LIG
     &    +Mulch_Rate_HCE*widthPerMulchUnit(n)
         total_Length=total_Length+widthPerMulchUnit(n)
        enddo
cccz C decomposition rate g (average sense)
        Mulch_Decompose_RDM_C=total_decomp_RDM/total_Length
        Mulch_Decompose_HCE_C=total_decomp_HCE/total_Length
        Mulch_Decompose_CEL_C=total_decomp_CEL/total_Length
        Mulch_Decompose_LIG_C=total_decomp_LIG/total_Length
cccz mulch decomposition rate g (average sense, sum of C from different component)       
        Mulch_Decompose=(Mulch_Decompose_RDM_C+Mulch_Decompose_HCE_C
     &    +Mulch_Decompose_CEL_C+Mulch_Decompose_LIG_C)/PC_mulch
cccz N decomposition rate g (average sense)        
        Mulch_Decompose_RDM_N=Mulch_Decompose_RDM_C/PC_mulch*PN_mulch
        Mulch_Decompose_HCE_N=Mulch_Decompose_HCE_C/PC_mulch*PN_mulch
        Mulch_Decompose_CEL_N=Mulch_Decompose_CEL_C/PC_mulch*PN_mulch
        Mulch_Decompose_LIG_N=Mulch_Decompose_LIG_C/PC_mulch*PN_mulch
cccz based on the rectangular based design, we have to assume the mulch decomposed homogeneously
cccz so just calculate a fraction for each layer
        Mulch_Avail=rho_mulch*MulchEleMarkerArea(EleIndex)
     &      *(1.0D0-f_mulch_pore)
        Frac_Decomp(1)=min(Mulch_Decompose/Mulch_Avail,1.0D0)
        if(mulchLayer.gt.1) then
         do jj=2,mulchLayer
           EleIndexUp=MulchEleMatrix(jj,1) 
           EleIndexDown=MulchEleMatrix(jj-1,1) 
           Mulch_AvailUp=rho_mulch*MulchEleMarkerArea(EleIndexUp)
     &      *(1.0D0-f_mulch_pore)
           Mulch_AvailDown=rho_mulch*MulchEleMarkerArea(EleIndexDown)
     &      *(1.0D0-f_mulch_pore)
           Frac_Decomp(jj)=min(Mulch_Decompose
     &       *Frac_Coef(jj)/Mulch_AvailUp,1.0D0)
cccz update the decomposition due to the feeding from upper layer
           Frac_Decomp(jj-1)=(Frac_Decomp(jj-1)*Mulch_AvailDown
     &      -Mulch_Decompose*Frac_Coef(jj))/Mulch_AvailDown
cccz update the mass fraction of each carbon component
           do kk=1,SurNodeIndex-1
              RDM_mass(jj,kk)=RDM_mass(jj,kk)
     &          -Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_RDM
              HCE_mass(jj,kk)=HCE_mass(jj,kk)
     &          -Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_HCE
              CEL_mass(jj,kk)=CEL_mass(jj,kk)
     &          -Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_CEL
              LIG_mass(jj,kk)=LIG_mass(jj,kk)
     &          -Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_LIG
           enddo
           if(jj.eq.2) then
            do kk=1,SurNodeIndex-1
              RDM_mass(1,kk)=RDM_mass(1,kk)-Mulch_Decompose_RDM_C
     &          +Mulch_Decompose*Frac_Coef(2)*PC_mulch*Frac_RDM
              HCE_mass(1,kk)=HCE_mass(1,kk)-Mulch_Decompose_HCE_C
     &          +Mulch_Decompose*Frac_Coef(2)*PC_mulch*Frac_HCE
              CEL_mass(1,kk)=CEL_mass(1,kk)-Mulch_Decompose_CEL_C
     &          +Mulch_Decompose*Frac_Coef(2)*PC_mulch*Frac_CEL
              LIG_mass(1,kk)=LIG_mass(1,kk)-Mulch_Decompose_LIG_C
     &          +Mulch_Decompose*Frac_Coef(2)*PC_mulch*Frac_LIG
            enddo
           else
            do kk=1,SurNodeIndex-1
              RDM_mass(jj-1,kk)=RDM_mass(jj-1,kk)
     &          +Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_RDM
              HCE_mass(jj-1,kk)=HCE_mass(jj-1,kk)
     &          +Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_HCE
              CEL_mass(jj-1,kk)=CEL_mass(jj-1,kk)
     &          +Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_CEL
              LIG_mass(jj-1,kk)=LIG_mass(jj-1,kk)
     &          +Mulch_Decompose*Frac_Coef(jj)*PC_mulch*Frac_LIG
            enddo   
           endif
          enddo
         endif 
       ENDIF
       
       return
3116   Write(*,*) 'Soil file not found'      
3110   Call errmes(im,il)
       return
      END