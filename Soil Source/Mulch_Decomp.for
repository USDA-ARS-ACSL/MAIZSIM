cccz -----------------------------------------------------
c     the mulch decomposition module
c     linked with mulch process module
cccz -----------------------------------------------------
      
cccz -------------------------------------------
c     variable used in mulch decomposition module (some of the variables are defined in public ins)
c     RDM-Rapidly Decomposed Mulch, HCE-Hemi-celluloses Mulch, CEL-Cellulose Mulch, LIG-Lignin Mulch
c     INTEGER TYPE
c         
c         EleIndex            The Element index for the current element based on the grid established in "SurfaceMulchAdjustment.for"
c         InString_mulch*132  For temporary holding a string during input reading, the target should be "[mulch decomposition]"
c         il, im              Line number and cursor position
c     
c     FLOAT TYPE  [UNIT]
c
c         PC_mulch        [g/g]   Mass concentration of C in mulch, g of C/g of solid portion of mulch, should be 0-1, from literature, e.g. in (https://www.pnas.org/content/115/16/4033), PC=0.4368
c         PN_mulch        [g/g]   Mass concentration of N in mulch, g of N/g of solid portion of mulch, should be 0-1, from literature, e.g. in (https://www.pnas.org/content/115/16/4033), PN=0.0140
c         alpha_F         [1]     The mulch feeding parameter, if the decomposition on soil-mulch interface is "a", then the decompositon for one layer above should be "alpha_mulch*a"
c         Frac_CARB       [1]     Fraction of Carbohydrate over Total Organic Mulch, should be 0-1, change during the simulation
c         Frac_CELL       [1]     Fraction of Holo-Cellulose over Total Organic Mulch, should be 0-1, change during the simulation
c         Frac_LIGN       [1]     Fraction of Lignin over Total Organic Mulch, should be 0-1, change during the simulation
c         Frac_CARB_Init  [1]     Fraction of Initial Carbohydrate over Total Organic Mulch, should be 0-1, fixed before chemical decomposition
c         Frac_CELL_Init  [1]     Fraction of Initial Holo-Cellulose over Total Organic Mulch, should be 0-1, fixed before chemical decomposition
c         Frac_LIGN_Init  [1]     Fraction of Initial Lignin over Total Organic Mulch, should be 0-1, fixed before chemical decomposition 
c         FracN_CARB_Init [1]     Fraction of Initial N in Carbohydrate over Carbohydrate, should be 0-1, fixed before chemical decomposition
c         FracN_CELL_Init [1]     Fraction of Initial N in Holo-Cellulose over Holo-Cellulose, should be 0-1, fixed before chemical decomposition
c         FracN_LIGN_Init [1]     Fraction of Initial N in Lignin over Lignin, should be 0-1, fixed before chemical decomposition 
c         Humid_Factor    [1]     Humification factor, 0.125D0
c         K_CARB          [day-1] Coef of decompositon rate, reference factor for CARB, need adjusted
c         K_CELL          [day-1] Coef of decompositon rate, reference factor for CELL, need adjusted
c         K_LIGN          [day-1] Coef of decompositon rate, reference factor for LIGN, need adjusted
c         K_CARB_temp     [day-1] Coef of decompositon rate for CARB, after adjusted
c         K_CELL_temp     [day-1] Coef of decompositon rate for CELL, after adjusted
c         K_LIGN_temp     [day-1] Coef of decompositon rate for LIGN, after adjusted
c         mulch_mass_temp [g]     All the availiable soild mulch in the target mulch element
c         N_mass_temp     [g]     All the availiable N in the target mulch element
c         INKG            [g]     Microbe-available inorganic N in shallow soil layer
c         CNR             [1] the universal C/R ratio for an given element 
c         Frac_Decomp(Layer)  [1] the decomposition fraction of each layer, 0 means no decomposition, 1 means totally gone
c         CO2mass_2_Cmass     [1] convert CO2 mass to C mass, 12/44
c         NO3mass_2_Nmass     [1] convert NO3 mass to N mass, 14/62  
c         NH4mass_2_Nmass     [1] convert NH4 mass to N mass, 14/18
c         ThetaV_2_ThetaM     [1] convert volumetric water content to gravimetric water content, there ThetaV is based on the whole mulch not based on mulch solid, so should be devided by (1-f_mulch_pore) first, then use mulch bulk density
c         HNew_m          [MPa]   Mulch water potential
c         Theta_m         [g/g]   Mulch gravimetric water content, refer to the total mulch
c         Temp_m          [oC]    Mulch temperature
c
c ------------------------mass pools---------------------------------------------
c         CARB_mass       [g]     CARB mass of the mulch, CARB_mass+CELL_mass+LIGN_mass=mulch_mass_temp
c         CELL_mass       [g]     CELL mass of the mulch, CARB_mass+CELL_mass+LIGN_mass=mulch_mass_temp
c         LIGN_mass       [g]     LIGN mass of the mulch, CARB_mass+CELL_mass+LIGN_mass=mulch_mass_temp
c         CARB_N_mass     [g]     N within CARB, CARB_N_mass+CELL_N_mass+LIGN_N_mass=N_mass_temp
c         CELL_N_mass     [g]     N within CELL, CARB_N_mass+CELL_N_mass+LIGN_N_mass=N_mass_temp
c         LIGN_N_mass     [g]     N within LIGN, CARB_N_mass+CELL_N_mass+LIGN_N_mass=N_mass_temp
c
c ------------------------decomposition rate---------------------------------------------
c         CNRF            [0-1]   CNR effects on decomposition
c         MTRF            [0-1]   Moisture and Temperature effects on decomposition
c         N_im            [g/day] the immobilization of N before N reach the soil surface (N interception by mulch)
c         CARB_Decomp     [g/day] Decomposition speed for CARB, after add CNRF and MTRF
c         CELL_Decomp     [g/day] Decomposition speed for CELL, after add CNRF and MTRF
c         LIGN_Decomp     [g/day] Decomposition speed for LIGN, after add CNRF and MTRF
c         FOM_Decomp      [g/day] Decomposition speed for all the organic matter (mulch residue)
c         CARB_Decomp_N   [g/day] Gross Decomposition speed for N in CARB, after add CNRF and MTRF
c         CELL_Decomp_N   [g/day] Gross Decomposition speed for N in CELL, after add CNRF and MTRF
c         LIGN_Decomp_N   [g/day] Gross Decomposition speed for N in LIGN, after add CNRF and MTRF
c         FOMN_Decomp     [g/day] Gross/Actual Decomposition speed for N in organic matter (mulch residue N)
c         CARB_total_decomp_R   [g] Total residue decomposition of CARB for mulch-soil interface
c         CELL_total_decomp_R   [g] Total residue decomposition of CELL for mulch-soil interface
c         LIGN_total_decomp_R   [g] Total residue decomposition of LIGN for mulch-soil interface
c         CARB_total_decomp_N   [g] Total gross N mineralization of CARB for mulch-soil interface  
c         CELL_total_decomp_N   [g] Total gross N mineralization of CELL for mulch-soil interface 
c         LIGN_total_decomp_N   [g] Total gross N mineralization of LIGN for mulch-soil interface
c         Nim_total_decomp      [g] Total N mulch-interception
c         Nmine_total_decomp    [g] Total N mineralization, >0 mulch to soil; <0 soil to mulch
c         Nhumi_total_decomp    [g] Total N humification, >0 mulch to soil
c         Mulch_Decompose [g]     Total decomposition on mulch-soil interface based on MULCH (NOT C-BASED HERE)      
c         Mulch_Avail     [g]     Total availiable mulch for the current layer (horizontal level)
c         Mulch_AvailUp   [g]     Total availiable mulch for the current layer (horizontal level) when considering a two layered structure
c         Mulch_AvailDown [g]     Total availiable mulch for the lower layer (horizontal level) when considering a two layered structure
c         totalMulchWidth [cm]    the whole mulch width
c         CO2_to_Air      [nu g CO2/cm^3] the concentration of CO2 added to air based on mulch decomposition and mulch volume.
c     DUMMY VARIABLES
c         bbbb                    dummy variable, receive unwanted data.

      subroutine MulchDecomposition()
       include 'public.ins'
       include 'puplant.ins'
       include 'puweath.ins'
       include 'PuSurface.ins'
       include 'nitvar.ins'
       
       double precision 
     !   Frac_CARB_Init, Frac_CELL_Init, Frac_LIGN_Init,          !cccz initial fraction of each component of soil organic matter
     !   FracN_CARB_Init, FracN_CELL_Init, FracN_LIGN_Init,       !cccz INPUT: initial fraction of N within each residue pool
     !   Frac_CARB, Frac_CELL, Frac_LIGN,                         !cccz record the fraction of each category at soil surface (bottom mulch layer)
     !   Humid_Factor,
     !   Frac_Decomp_0(NumMulHLD),
     !   Frac_Decomp_Crit,Frac_Decomp_Crit_0,
     !   Crit_Level, Frac_Crit_Layer,
     !   CARB_mass_cont,CELL_mass_cont,LIGN_mass_cont,
     !   CARB_N_mass_cont,CELL_N_mass_cont,LIGN_N_mass_cont,
     !   CARB_mass_non_cont,CELL_mass_non_cont,LIGN_mass_non_cont,
     !   CARB_N_mass_non_cont,CELL_N_mass_non_cont,LIGN_N_mass_non_cont,
     !   CARB_mass_chg,CELL_mass_chg,LIGN_mass_chg,
     !   CARB_N_mass_chg,CELL_N_mass_chg,LIGN_N_mass_chg,
     !   thick_cont, thick_non_cont  
       double precision 
     !   PC_mulch, PN_mulch, alpha_F,
     !   K_CARB, K_CELL, K_LIGN,
     !   K_CARB_temp, K_CELL_temp, K_LIGN_temp,       
     !   ThetaV_2_ThetaM, 
     !   INKG, CNR,
     !   CNRF, MTRF, N_im
       double precision
     !   HNew_m, Theta_m, Temp_m
       integer 
     !   EleIndex,il,im, qLeft, qRight,
     !   Crit_Layer
       integer LocalFlag_MulchDecomp_Start,
     !   LocalFlag_MulchDecomp_FinalAssign  
       double precision 
     !   mulch_mass_temp, mulch_mass_N_temp, N_mass_temp
       double precision
     !   Mulch_Avail,Mulch_AvailUp,Mulch_AvailDown
       double precision 
     !   CARB_Decomp, CELL_Decomp, LIGN_Decomp, FOM_Decomp,  
     !   CARB_Decomp_N, CELL_Decomp_N, LIGN_Decomp_N, FOMN_Decomp,
     !   CARB_total_decomp_R, CELL_total_decomp_R,
     !   LIGN_total_decomp_R, CARB_total_decomp_N,
     !   CELL_total_decomp_N, LIGN_total_decomp_N,
     !   Nim_total_decomp,Nmine_total_decomp,Nhumi_total_decomp,
     !   Nim_cumu_decomp,Nmine_cumu_decomp,Nhumi_cumu_decomp,
     !   Mulch_Decompose,local_soil_mass,
     !   totalMulchC_final, totalMulchN_final  
       double precision N_Quan_Left,N_Quan_Right,
     !   Length_Left_Frac,Length_Right_Frac  
       Character InString_mulch*132
cccz just for paper writting
       double precision Tentative_time_output,
     !    Tentative_output,
     !    t_CNR, t_CARB, t_CELL, t_LIGN, 
     !    t_CARB_N, t_CELL_N, t_LIGN_N
cccz dummy variable, can be every thing, but explainable based on it neighbor lines       
       double precision aaaa,bbbb
       
cccztest just for test
c       double precision dummy_total_N_old,dummy_total_N_new
cccz ------------------------------------------------------------       
       
       common /MulchDecompositionBlock/ 
     !   PC_mulch,PN_mulch,alpha_F,
     !   ThetaV_2_ThetaM,
     !   Frac_CARB_Init, Frac_CELL_Init, Frac_LIGN_Init,
     !   FracN_CARB_Init,FracN_CELL_Init,FracN_LIGN_Init,
     !   Frac_CARB, Frac_CELL, Frac_LIGN,
     !   Humid_Factor,
     !   K_CARB, K_CELL,K_LIGN,
     !   INKG, Crit_Level,
     !   Nim_cumu_decomp,Nmine_cumu_decomp,Nhumi_cumu_decomp,
     !   Tentative_time_output,
     !   Tentative_output,
     !   LocalFlag_MulchDecomp_Start,
     !   LocalFlag_MulchDecomp_FinalAssign  
       
cccz ----------------------- start the initialization -------------------------

cccz initialization step 1: read the file and initize some parameter,
c    the residue mulch is not established yet.
c    the trager will be  "lInput.eq.1"
       

      IF(lInput.eq.1) then   
cccz make several initialization
        Frac_CARB_Init=0.20D0
        Frac_CELL_Init=0.70D0
        Frac_LIGN_Init=0.10D0
        FracN_CARB_Init=0.08D0            ! this is only an approximation
        FracN_CELL_Init=0.01D0
        FracN_LIGN_Init=0.01D0
        Humid_Factor=0.125D0
        K_CARB=0.36D0
        K_CELL=0.24D0
        K_LIGN=0.0228D0
        PC_mulch=0.41D0                   ! Concentrations (https://www.pnas.org/content/115/16/4033)
        PN_mulch=0.01D0                   ! Concentrations (https://www.pnas.org/content/115/16/4033)
        alpha_F=0.1D0
        INKG=0.5D0
        Crit_Level=0.45D0
        
        im=80
        il=0   
        Open(21,file=MulchFile,status='old',ERR=3116)
3115    Read (21,'(A132)') InString_mulch
cccz '[Mulch_Decomposition_2]' is really temporary
cccz some discussion of input is needed
        if (InString_mulch(1:23).ne.'[Mulch_Decomposition]') goto 3115
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
cccz g/g for C and N concentration, 
c      from literature, e.g. in (https://www.pnas.org/content/115/16/4033), PC=0.4368, PN=0.014,
cccz alpha_F for feeding (upper layer to lower layer)
        Read(21,*,ERR=3110) Crit_Level, alpha_F  
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
cccz fraction of three carbon forms easy decomposition -> very hard decomposition
        Read(21,*,ERR=3110) Frac_CARB_Init,Frac_CELL_Init,Frac_LIGN_Init
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
cccz fraction of three carbon forms easy decomposition -> very hard decomposition
        Read(21,*,ERR=3110) FracN_CARB_Init,
     &            FracN_CELL_Init, FracN_LIGN_Init
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
        Read(21,*,ERR=3110)
        im=im+1
        il=il+1
cccz coefficient for first order dynamics (decomposition equation)
        Read(21,*,ERR=3110) K_CARB, K_CELL, K_LIGN
        close(21)
cccz need to calculate a TUpperL_mulch/TLowerL_mulch
c    just need to run the function once 
          bbbb=WQ_CERES_MULCH(-100.0D0,0.1D0,0.1D0,-1)
          LocalFlag_MulchDecomp_Start=0 
          LocalFlag_MulchDecomp_FinalAssign=0
         if(Crit_Level.eq.0) then ! that means we do not want to have decomposition
             DecompCtrl=0
         endif    
      EndIf

c---------------------------- finish the input --------------------------------------------------------------------------------
cccz initialization step 2: REAL INITIALIZATION, occur once when the mulch starting time is hitted
c    the residue mulch is established at this time step
c    the trager will be  "MulchDecompIni.eq.1"

       IF(MulchDecompIni.eq.1
     &     .and.LocalFlag_MulchDecomp_Start.eq.0) then               ! double ensure this seciton will be only executed once
        LocalFlag_MulchDecomp_Start=1 
        LocalFlag_MulchDecomp_FinalAssign=0
cccz set up the critical level, the mulch will be divided into two layers, one is the contacting layer, one is the non cotacting layer
        Crit_Level=Crit_Level*mulchThick
cccz start the code based on CERES.for
c do a calculation on residue pool fraction
        Frac_CARB=Frac_CARB_Init
        Frac_CELL=Frac_CELL_Init
        Frac_LIGN=Frac_LIGN_Init
c some parameter (constant)        
c       CO2mass_2_Cmass=12.0107D0/44.01D0
c       NO3mass_2_Nmass=14.0067D0/62.0049D0
c       NH4mass_2_Nmass=14.0067D0/18.039D0
        ThetaV_2_ThetaM=1.0D0/(rho_mulch_b/rho_w)  !cccz we note here that density follows g m^-3, so we have to add rho_w here.
cccz initialize the mass/mass N for each category
        do kk=1,SurNodeIndex_M-1
         do jj=1,mulchLayer
          EleIndex=MulchEleMatrix(jj,kk)
cccz initialize the total mass, propotional to carbon
          mulch_mass_temp=MulchEleMarkerArea(EleIndex)        !cccz mulch area has unit cm^2*(slab thickness) for each element, and the width of the slab is assumed to be 1 cm, but the unit of mulch density is "rho_mulch" g m^-3
     &     *rho_mulch_b*1.0D-6                                !cccz there fore, we need /1.0D6 to convert m^3 to cm^3
          ! CARB
          aaaa=mulch_mass_temp*Frac_CARB_Init                 !cccz the mass pools of CARB the three categories of organic carbon
          CARB_mass(jj,kk)=aaaa
          CARB_N_mass(jj,kk)=aaaa*FracN_CARB_Init
          ! CELL
          aaaa=mulch_mass_temp*Frac_CELL_Init                 !cccz the mass pools of CELL the three categories of organic carbon
          CELL_mass(jj,kk)=aaaa
          CELL_N_mass(jj,kk)=aaaa*FracN_CELL_Init
          ! LIGN
          aaaa=mulch_mass_temp*Frac_LIGN_Init                 !cccz the mass pools of LIGN the three categories of organic carbon
          LIGN_mass(jj,kk)=aaaa
          LIGN_N_mass(jj,kk)=aaaa*FracN_LIGN_Init
         enddo
       enddo
       Nim_cumu_decomp=0.0D0
       Nmine_cumu_decomp=0.0D0
       Nhumi_cumu_decomp=0.0D0
cccz ------------------------ initial output for paper writting-----------------
cccz initialize the mulch output options
c       Tentative_time_output=1.0D0/24.0D0
c       Tentative_output=time+Tentative_time_output
c       Open(9998,file='D:\MaizsimRunoff\BGR_B\Mulch_decomp_paper.out')
c       Write(9998,9903) "Date_time", "Nim", "NDecomp", "Nhumi", "CNR", 
c     &    "CARB", "CELL", "LIGN", "CARB_N", "CELL_N", "LIGN_N"
c      Open(9999,file='D:\MaizsimRunoff\BGR_C_NB\Mulch_decomp_N.out')
c      Write(9999,9905) "Date_time", "RunningError1", "RunningError2"
cccz -----------------------------------------------------------------------
      return
      ENDIF
cccz ----------------------- finish the initialization ------------------------- 
      
cccz ----------------------- start the calculation -------------------------
      
cccz no decomposition occurs
      
      if(DecompCtrl.eq.0) goto 3117
      if(LocalFlag_MulchDecomp_Start.eq.0) goto 3117
      
cccz by the end of decomposition, assign mulch C and N to soil, just once
      if(BoolMulch_TotalDecomposed.eq.1
     &  .and.LocalFlag_MulchDecomp_FinalAssign.eq.0) then
          LocalFlag_MulchDecomp_FinalAssign=1
          goto 3118
      endif    
     
cccztest just for test
c      dummy_total_N_old=0.0D0 
c      dummy_total_N_new=0.0D0 
cccz ------------------------------------------------------------           
c      do jj=1,mulchLayer
c        do kk=1,SurNodeIndex_M-1 
c           dummy_total_N_old=dummy_total_N_old+CARB_N_mass(jj,kk)
c     &      +CELL_N_mass(jj,kk)+LIGN_N_mass(jj,kk)
c        enddo
c      enddo
cccz ------------------------------------------------------------ 
           
      
cccz if decomposition occurs
      
       CARB_total_decomp_R=0.0D0
       CELL_total_decomp_R=0.0D0
       LIGN_total_decomp_R=0.0D0
       CARB_total_decomp_N=0.0D0
       CELL_total_decomp_N=0.0D0
       LIGN_total_decomp_N=0.0D0
       Nim_total_decomp=0.0D0
       Nmine_total_decomp=0.0D0
       Nhumi_total_decomp=0.0D0
       Crit_Layer=0
       
cccz check the contact layer and non-contact layer
       bbbb=0.0D0
       Frac_Crit_Layer=-100.0D0
       do jj=1,mulchLayer
        bbbb=bbbb+thickPerLayer(jj)
        if(bbbb.le.Crit_Level) then
          Crit_Layer=jj
        else
         if(Frac_Crit_Layer.lt.-10.0D0) then
          Frac_Crit_Layer=1.0D0-(bbbb-Crit_Level)/thickPerLayer(jj)
          Frac_Crit_Layer=min(Frac_Crit_Layer,1.0D0)
          Frac_Crit_Layer=max(Frac_Crit_Layer,0.0D0)
         endif 
        endif 
       enddo  ! end mulchlayer loop 
       
cccz calculate decomposition fraction based on
       do jj=1,mulchLayer
         Frac_Decomp(jj)=0.0D0                         !cccz the actual decomposition factor for each layer, 0 means no decomposition, 1 means totally disappeared 
         Frac_Decomp_0(jj)=0.0D0
       enddo

       do kk=1,SurNodeIndex_M-1                          !cccz decomposition occured at mulch-soil interface, i.e., for the first layer of the mulch grid
        HNew_m=0.0D0
        Theta_m=0.0D0
        Temp_m=0.0D0
        mulch_mass_temp=0.0D0
        mulch_mass_N_temp=0.0D0
        CARB_mass_cont=0.0D0
        CELL_mass_cont=0.0D0
        LIGN_mass_cont=0.0D0
        CARB_N_mass_cont=0.0D0
        CELL_N_mass_cont=0.0D0
        LIGN_N_mass_cont=0.0D0
        thick_cont=0.0D0
        INKG=0.0D0
        
cccz --------------------- calculate the inorganic N from soil that availiable for microbe (mineralized) INKG ----------------------------
        do kkk=0,SurNodeIndex_M_Match(kk+1)-1
          kkkk=kkk+SurfMulchNodeSubIndex(kk)
          qLeft=SurfNodeSurfIndexH(kkkk)
          qRight=SurfNodeSurfIndexH(kkkk+1)
          if(kkkk.eq.1) then   ! left most point
            Length_Left_Frac=1.0D0
          else
            Length_Left_Frac=0.5D0
          endif
          if(kkkk.eq.SurNodeIndex) then   ! right most point
            Length_Right_Frac=1.0D0
          else
            Length_Right_Frac=0.5D0
          endif
          ! cccz we need to notice the unit of NO3 and NH4
          N_Quan_Left=(NO3_old(qLeft)
     &      +NH4(qLeft))
     &      *nodeArea(qLeft)*Length_Left_Frac
     &      *1.0D-6                                                ! convert "\mu g" to "g"
          N_Quan_Right=(NO3_old(qRight)
     &      +NH4(qRight))
     &      *nodeArea(qRight)*Length_Right_Frac       
     &      *1.0D-6                                                 ! convert "\mu g" to "g"
          INKG=INKG+N_Quan_Left+N_Quan_Right    
        enddo  ! end kkk=0,SurNodeIndex_M_Match(kk+1)-1
C ---------------------------------------------

        do jj=1,Crit_Layer
          thick_cont=thick_cont+thickPerLayer(jj)  
          HNew_m=HNew_m+log(abs(MulchElehnew(jj,kk)))
     &          *thickPerLayer(jj)        
          Theta_m=Theta_m+MulchEleThnew(jj,kk)*ThetaV_2_ThetaM
     &          *thickPerLayer(jj)      !cccz MulchEleThnew is mulch (not solid portion) based on volumetric water content for each element
          Temp_m=Temp_m+MulchEleTmpr(jj,kk)
     &          *thickPerLayer(jj)                         !cccz MulchEleTmpr is the mulch temperature
          EleIndex=MulchEleMatrix(jj,kk)
          mulch_mass_temp=mulch_mass_temp                         !cccz mulch area has unit cm^2*(slab thickness) for each element, and the width of the slab is assumed to be 1 cm, but the unit of mulch density is "rho_mulch" g m^-3
     &     +MulchEleMarkerArea(EleIndex)*rho_mulch_b*1.0D-6
          CARB_mass_cont=CARB_mass_cont+CARB_mass(jj,kk)
          CELL_mass_cont=CELL_mass_cont+CELL_mass(jj,kk)
          LIGN_mass_cont=LIGN_mass_cont+LIGN_mass(jj,kk)
          CARB_N_mass_cont=CARB_N_mass_cont+CARB_N_mass(jj,kk)
          CELL_N_mass_cont=CELL_N_mass_cont+CELL_N_mass(jj,kk)
          LIGN_N_mass_cont=LIGN_N_mass_cont+LIGN_N_mass(jj,kk)
        enddo !end jj loop 1 to Crit_Layer
C------------------------------------------        
        if(Crit_Layer.lt.mulchLayer.and.Frac_Crit_Layer.gt.0.0D0) then
         jj=Crit_Layer+1
         thick_cont=thick_cont+thickPerLayer(jj)*Frac_Crit_Layer
         HNew_m=HNew_m+log(abs(MulchElehnew(jj,kk)))
     &    *(thickPerLayer(jj)*Frac_Crit_Layer)        !cccz should convert “cm” water potential to “Mpa” later
         Theta_m=Theta_m+MulchEleThnew(jj,kk)*ThetaV_2_ThetaM
     &    *(thickPerLayer(jj)*Frac_Crit_Layer)     !cccz MulchEleThnew is mulch (not solid portion) based on volumetric water content for each element
         Temp_m=Temp_m+MulchEleTmpr(jj,kk)
     &    *(thickPerLayer(jj)*Frac_Crit_Layer)                      !cccz MulchEleTmpr is the mulch temperature
         EleIndex=MulchEleMatrix(jj,kk)
         mulch_mass_temp=mulch_mass_temp                         !cccz mulch area has unit cm^2*(slab thickness) for each element, and the width of the slab is assumed to be 1 cm, but the unit of mulch density is "rho_mulch" g m^-3
     &    +MulchEleMarkerArea(EleIndex)*rho_mulch_b*1.0D-6
     &    *Frac_Crit_Layer   
         CARB_mass_cont=CARB_mass_cont
     &    +CARB_mass(jj,kk)*Frac_Crit_Layer   
         CELL_mass_cont=CELL_mass_cont
     &    +CELL_mass(jj,kk)*Frac_Crit_Layer   
         LIGN_mass_cont=LIGN_mass_cont
     &    +LIGN_mass(jj,kk)*Frac_Crit_Layer   
         CARB_N_mass_cont=CARB_N_mass_cont
     &    +CARB_N_mass(jj,kk)*Frac_Crit_Layer   
         CELL_N_mass_cont=CELL_N_mass_cont
     &    +CELL_N_mass(jj,kk)*Frac_Crit_Layer   
         LIGN_N_mass_cont=LIGN_N_mass_cont
     &    +LIGN_N_mass(jj,kk)*Frac_Crit_Layer  
        endif ! end if Crit_Layer.lt.mulchLayer
              
        HNew_m=-exp(HNew_m/thick_cont)
        HNew_m=HNew_m*0.0000980665D0                  !cccz convert “cm” water potential to “Mpa”
        Theta_m=Theta_m/thick_cont
        Temp_m=Temp_m/thick_cont
        
        Frac_CARB=CARB_mass_cont/mulch_mass_temp
        Frac_CELL=CELL_mass_cont/mulch_mass_temp
        Frac_LIGN=LIGN_mass_cont/mulch_mass_temp
        
cccz ---------------------- decomposition dynamics for the total biomass (propotional to carbon) -------------------------------
        K_CARB_Temp=K_CARB*exp(-12.0D0*Frac_LIGN)
        K_CELL_Temp=K_CELL*exp(-12.0D0*Frac_LIGN)
        K_LIGN_Temp=K_LIGN
        
cccz determine the universal C/N ratio and CNRF factor
        CNR=PC_mulch*mulch_mass_temp
     &    /(CARB_N_mass_cont+CELL_N_mass_cont+LIGN_N_mass_cont
     &    +INKG)
        if(CNR.ge.13.0D0) then
          CNRF=min(exp(-0.693D0*(CNR-13.0D0)/13.0D0),1.0D0)
        else
          CNRF=1.0D0
        endif
cccz determine temperature/moisture effet MTRF
        if(Temp_m.ge.0.0D0) then
          MTRF=(0.384D0+0.018D0*Temp_m)
     &     *exp((0.142D0+0.628D0/Temp_m)*HNew_m)
          MTRF=min(max(MTRF,0.0D0),1.0D0)
        else
          MTRF=0.0D0
        endif
      
cccz Here we require "POOL_mass[g]" and "K_POOL_Temp[day^-1]". 
cccz We implicitly divide CARB_mass by area to use Resham's equation, then time the area to obtain mass         
        CARB_Decomp=K_CARB_Temp*CARB_mass_cont*MTRF*CNRF ! [g/day] 
        CELL_Decomp=K_CELL_Temp*CELL_mass_cont*MTRF*CNRF ! [g/day]
        LIGN_Decomp=K_LIGN_Temp*LIGN_mass_cont*MTRF*CNRF ! [g/day]
        FOM_Decomp=CARB_Decomp+CELL_Decomp+LIGN_Decomp
        
cccz ---------------------- decomposition dynamics for N (not propotional to carbon) -------------------------------
        CARB_Decomp_N=K_CARB_Temp*CARB_N_mass_cont*MTRF*CNRF        ! [g/day] 
        CELL_Decomp_N=K_CELL_Temp*CELL_N_mass_cont*MTRF*CNRF        ! [g/day]
        LIGN_Decomp_N=K_LIGN_Temp*LIGN_N_mass_cont*MTRF*CNRF        ! [g/day]
        FOMN_Decomp=CARB_Decomp_N+CELL_Decomp_N+LIGN_Decomp_N        ! the gross value from the above three
        
cccz convert the gross decomposition to net decomposition        
        N_im=min(FOM_Decomp*0.0213D0-FOMN_Decomp, INKG/Step)         ! [g/day]   
        N_im=max(N_im, 0.0D0)                                        ! [g/day] 
        FOMN_Decomp=FOMN_Decomp*(1.0D0-Humid_Factor)-N_im            ! [g/day] 
        FOMN_Humi=FOMN_Decomp*Humid_Factor                           ! [g/day] 

cccz ---------------------- sum over the horizontal scale, (equivalent to) take the average use the mass base sum --------------------------------
        CARB_total_decomp_R=CARB_total_decomp_R+
     &        min(CARB_Decomp*Step,CARB_mass_cont)
        CELL_total_decomp_R=CELL_total_decomp_R+
     &        min(CELL_Decomp*Step,CELL_mass_cont)
        LIGN_total_decomp_R=LIGN_total_decomp_R+
     &        min(LIGN_Decomp*Step,LIGN_mass_cont)
        CARB_total_decomp_N=CARB_total_decomp_N+
     &        min(CARB_Decomp_N*Step,CARB_N_mass_cont)
        CELL_total_decomp_N=CELL_total_decomp_N+
     &        min(CELL_Decomp_N*Step,CELL_N_mass_cont)
        LIGN_total_decomp_N=LIGN_total_decomp_N+
     &        min(LIGN_Decomp_N*Step,LIGN_N_mass_cont)
        Nim_total_decomp=Nim_total_decomp+N_im*Step
        Nmine_total_decomp=Nmine_total_decomp+FOMN_Decomp*Step
        Nhumi_total_decomp=Nhumi_total_decomp+FOMN_Humi*Step

       enddo    ! end kk loop 1 to SurNodeIndex_M-1
     
C ----------------------------------------       
      
       Nim_cumu_decomp=Nim_cumu_decomp+Nim_total_decomp               ! go back to CARB_N
       Nmine_cumu_decomp=Nmine_cumu_decomp+Nmine_total_decomp         ! N outlet
       Nhumi_cumu_decomp=Nhumi_cumu_decomp+Nhumi_total_decomp         ! N outlet

cccz mulch decomposition quantity for the whole layer (sum of all the residue pools), [g]
       Mulch_Decompose=CARB_total_decomp_R+CELL_total_decomp_R
     &   +LIGN_total_decomp_R
       
cccz based on the rectangular based design, we have to assume the mulch decomposed homogeneously
cccz so just calculate a fraction for the mulch-soil surface layer
       thick_cont=0.0D0
       do jj=1,Crit_Layer
         thick_cont=thick_cont+thickPerLayer(jj)
       enddo
       jj=Crit_Layer+1
       thick_cont=thick_cont+thickPerLayer(jj)*Frac_Crit_Layer
       Mulch_Avail=rho_mulch_b*1.0D-6
     &    *totalMulchWidth*thick_cont
       Frac_Decomp_Crit_0=min(Mulch_Decompose/Mulch_Avail,1.0D0)
       Frac_Decomp_Crit=Frac_Decomp_Crit_0
       do jj=1,Crit_Layer
        Frac_Decomp(jj)=Frac_Decomp_Crit_0
        Frac_Decomp_0(jj)=Frac_Decomp_Crit_0
       enddo
              
cccz C and N only show for the first layer
cccz for upper layers, the C and N are consider to transferred to lower layers as raw mulch material (inital C and N fraction)
       if(mulchLayer.gt.Crit_Layer) then
        thick_non_cont=max(mulchThick-thick_cont,0.0D0)
        Mulch_AvailUp=rho_mulch_b*1.0D-6
     &    *totalMulchWidth*thick_non_cont
        Mulch_AvailDown=rho_mulch_b*1.0D-6
     &    *totalMulchWidth*thick_cont
        bbbb=min(Mulch_Decompose*alpha_F/Mulch_AvailUp,1.0D0)
        do jj=Crit_Layer+1,mulchLayer
          Frac_Decomp(jj)=bbbb
          Frac_Decomp_0(jj)=bbbb
        enddo
        Frac_Decomp_Crit=max((Frac_Decomp_Crit_0*Mulch_AvailDown
     &      -bbbb*Mulch_AvailUp)/Mulch_AvailDown,0.0D0)
        do jj=1,Crit_Layer
cccz update the decomposition due to the feeding from upper layer
          Frac_Decomp(jj)=Frac_Decomp_Crit
        enddo
cccz update the mass fraction of each carbon component
cccz layer in contact level
        
        CARB_mass_non_cont=0.0D0
        CELL_mass_non_cont=0.0D0
        LIGN_mass_non_cont=0.0D0
        CARB_N_mass_non_cont=0.0D0
        CELL_N_mass_non_cont=0.0D0
        LIGN_N_mass_non_cont=0.0D0
        
        jj=Crit_Layer+1
        do kk=1,SurNodeIndex_M-1
          CARB_mass_non_cont=CARB_mass_non_cont+CARB_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
          CELL_mass_non_cont=CELL_mass_non_cont+CELL_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
          LIGN_mass_non_cont=LIGN_mass_non_cont+LIGN_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
          CARB_N_mass_non_cont=CARB_N_mass_non_cont+CARB_N_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
          CELL_N_mass_non_cont=CELL_N_mass_non_cont+CELL_N_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
          LIGN_N_mass_non_cont=LIGN_N_mass_non_cont+LIGN_N_mass(jj,kk)
     &     *(1.0D0-Frac_Crit_Layer)
        enddo
        if(mulchLayer.gt.Crit_Layer+1) then
         do jj=Crit_Layer+2,mulchLayer
          do kk=1,SurNodeIndex_M-1
           CARB_mass_non_cont=CARB_mass_non_cont+CARB_mass(jj,kk)
           CELL_mass_non_cont=CELL_mass_non_cont+CELL_mass(jj,kk)
           LIGN_mass_non_cont=LIGN_mass_non_cont+LIGN_mass(jj,kk)
           CARB_N_mass_non_cont=CARB_N_mass_non_cont+CARB_N_mass(jj,kk)
           CELL_N_mass_non_cont=CELL_N_mass_non_cont+CELL_N_mass(jj,kk)
           LIGN_N_mass_non_cont=LIGN_N_mass_non_cont+LIGN_N_mass(jj,kk)
          enddo
         enddo
        endif ! end if mulchLayer.gt.Crit_Layer+1
        bbbb=Frac_Decomp_0(Crit_Layer+1)
        CARB_mass_chg=CARB_total_decomp_R
     &      -bbbb*CARB_mass_non_cont
        CELL_mass_chg=CELL_total_decomp_R
     &      -bbbb*CELL_mass_non_cont
        LIGN_mass_chg=LIGN_total_decomp_R
     &      -bbbb*LIGN_mass_non_cont
        CARB_N_mass_chg=CARB_total_decomp_N
     &      -Nim_total_decomp                                     ! the immobilization part of N
     &      -bbbb*CARB_N_mass_non_cont
        CELL_N_mass_chg=CELL_total_decomp_N
     &      -bbbb*CELL_N_mass_non_cont
        LIGN_N_mass_chg=LIGN_total_decomp_N
     &      -bbbb*LIGN_N_mass_non_cont
        
        do jj=1,Crit_Layer
         do kk=1,SurNodeIndex_M-1
           bbbb=widthPerMulchUnit(kk)/totalMulchWidth
     &       *thickPerLayer(jj)/thick_cont
           CARB_mass(jj,kk)=CARB_mass(jj,kk)-CARB_mass_chg*bbbb
           CELL_mass(jj,kk)=CELL_mass(jj,kk)-CELL_mass_chg*bbbb
           LIGN_mass(jj,kk)=LIGN_mass(jj,kk)-LIGN_mass_chg*bbbb
c note N is not participate on the Frac_Decomp, hence not related to shrinking
           CARB_N_mass(jj,kk)=CARB_N_mass(jj,kk)-CARB_N_mass_chg*bbbb
           CELL_N_mass(jj,kk)=CELL_N_mass(jj,kk)-CELL_N_mass_chg*bbbb
           LIGN_N_mass(jj,kk)=LIGN_N_mass(jj,kk)-LIGN_N_mass_chg*bbbb
         enddo
        enddo
cccz right at the critical layer
        jj=Crit_Layer+1
        aaaa=Frac_Decomp(jj)*(1-Frac_Crit_Layer)
        do kk=1,SurNodeIndex_M-1
           bbbb=widthPerMulchUnit(kk)/totalMulchWidth
     &       *thickPerLayer(jj)*Frac_Crit_Layer/thick_cont 
           CARB_mass(jj,kk)=CARB_mass(jj,kk)
     &      -CARB_mass_chg*bbbb-CARB_mass(jj,kk)*aaaa
           CELL_mass(jj,kk)=CELL_mass(jj,kk)
     &      -CELL_mass_chg*bbbb-CELL_mass(jj,kk)*aaaa
           LIGN_mass(jj,kk)=LIGN_mass(jj,kk)
     &      -LIGN_mass_chg*bbbb-LIGN_mass(jj,kk)*aaaa
           CARB_N_mass(jj,kk)=CARB_N_mass(jj,kk)
     &      -CARB_N_mass_chg*bbbb-CARB_N_mass(jj,kk)*aaaa
           CELL_N_mass(jj,kk)=CELL_N_mass(jj,kk)
     &      -CELL_N_mass_chg*bbbb-CELL_N_mass(jj,kk)*aaaa
           LIGN_N_mass(jj,kk)=LIGN_N_mass(jj,kk)
     &      -LIGN_N_mass_chg*bbbb-LIGN_N_mass(jj,kk)*aaaa
         enddo
         
cccz layer in the non_contact layer
cccz decomposition is based on raw mulch (feeding lower layer)
cccz update the organic matter and N pools (this is the propotional part)
       if(mulchLayer.gt.Crit_Layer+1) then
        do jj=Crit_Layer+2,mulchLayer
         bbbb=1-Frac_Decomp(jj) 
         do kk=1,SurNodeIndex_M-1
           CARB_mass(jj,kk)=bbbb*CARB_mass(jj,kk)
           CELL_mass(jj,kk)=bbbb*CELL_mass(jj,kk)
           LIGN_mass(jj,kk)=bbbb*LIGN_mass(jj,kk)
           CARB_N_mass(jj,kk)=bbbb*CARB_N_mass(jj,kk) 
           CELL_N_mass(jj,kk)=bbbb*CELL_N_mass(jj,kk) 
           LIGN_N_mass(jj,kk)=bbbb*LIGN_N_mass(jj,kk)
         enddo
        enddo
        endif ! end if mulchLayer.gt.Crit_Layer+1
        jj=Crit_Layer+1
        Frac_Decomp(jj)=(1-Frac_Crit_Layer)*Frac_Decomp(jj)
     &    +Frac_Crit_Layer*Frac_Decomp_Crit

      else !else for nulchlayer > crit_layer
cccz update the mass fraction of each carbon component
cccz layers in contact level 
       do jj=1,mulchLayer       
        do kk=1,SurNodeIndex_M-1
          bbbb=widthPerMulchUnit(kk)/totalMulchWidth
     &       *thickPerLayer(jj)/thick_cont
          CARB_mass(jj,kk)=CARB_mass(jj,kk)-CARB_total_decomp_R*bbbb
          CELL_mass(jj,kk)=CELL_mass(jj,kk)-CELL_total_decomp_R*bbbb
          LIGN_mass(jj,kk)=LIGN_mass(jj,kk)-LIGN_total_decomp_R*bbbb
c note N is not participate on the Frac_Decomp, hence not related to shrinking
          CARB_N_mass(jj,kk)=CARB_N_mass(jj,kk)
     &      -(CARB_total_decomp_N-Nim_total_decomp)*bbbb
          CELL_N_mass(jj,kk)=CELL_N_mass(jj,kk)-CELL_total_decomp_N*bbbb
          LIGN_N_mass(jj,kk)=LIGN_N_mass(jj,kk)-LIGN_total_decomp_N*bbbb
         enddo          
        enddo
       endif  ! end mulchLayer.gt.Crit_Layer)

c ---------------------------- output of N (whole surface) --------------------------------------------------------
c the data exchange interface between mulch N to soil N
c should adjust the unit here before leaving this module, so soilN model can take them as source/sink
c        Nmine_total_decomp   [g]     total amount of mineralized N, can be >0, =0, <0 (soil to mulch N flux) -> add to soil mineral N pool
c        Nhumi_total_decomp   [g]     total amount of humified N -> add to soil humified N pool
cdt we need to use BD here to convert to ug/cm3 volume for both.
cccz we need to assign N to soil here, so use "SurNodeIndex" instead of "SurNodeIndex_M"
       local_soil_mass=0.0D0
       do i=1,NumNP
          NO3_from_residue(i)=0.0D0 
       enddo
      
       do kk=1,SurNodeIndex
          qLeft=SurfNodeSurfIndexH(kk)
          qRight=SurfNodeNodeIndexH(kk)
          local_soil_mass=nodeArea(qRight)*BlkDn(MatNumN(qRight))   ! can be "BlkDn(1)" is because we are always on soil surface
c    add the mineral N from mulch to an array, which will be added to BNO3 in soilnden.for
c soilNden needs ug/cm3 of total volume =  N * width/total width /soil mass *BD =ug cm3
          NO3_from_residue(qRight)=
     &       real(Nmine_total_decomp*Width(qLeft)/totalMulchWidth
     &        /nodeArea(qRight)
     &        *1.0D6)      
c    add the humified N from mulch to soil surface nodes 
          Nh(qRight)=Nh(qRight)
     &       +real(Nhumi_total_decomp*Width(qLeft)/totalMulchWidth
     &        /nodeArea(qRight)
     &        *1.0D6)    
       enddo
          
cccz ---------------------------- output of C (whole surface) -------------------------------------------------------- 
cccz the assumption is C will be release as CO2 and form a high CO2 layer soil surface
       CO2_to_Air=Mulch_Decompose*PC_mulch*(44.0095D0/12.0107D0)  ! decomposed mulch mass [g]* C mass frac [g/g]* convert C mass to CO2
     &    /(mulchThick*totalMulchWidth*1.0D0)                     ! the whole volume occupied by mulch [cm^3], 1.0D0 cm indicates the scale at the 3rd dimention
     &    *10D6                                                   ! convert unit [g/cm^3] to [nu g CO2/ cm^3]
       
       do n=1,SurNodeIndex
          PG=11920.0                                        !For test, we need to add this as an input, also check for the range of this value
     
          kSur=SurfNodeSurfIndexH(n)
          do jjj=1,NumG
              Varbg_Mulch(kSur,jjj,1)=Varbg_Air(kSur,jjj,1)+CO2_to_Air
              Varbg_Mulch(kSur,jjj,2)=PG
              Varbg_Mulch(kSur,jjj,3)=Varbg_Mulch(kSur,jjj,1)*PG
              Varbg(kSur,jjj,1)=Varbg_Mulch(kSur,jjj,1)
              Varbg(kSur,jjj,2)=Varbg_Mulch(kSur,jjj,2)
              Varbg(kSur,jjj,3)=Varbg_Mulch(kSur,jjj,3)
          enddo
      enddo
     
cccztest just for test
c      dummy_total_N_new=0.0D0 
cccz ------------------------------------------------------------           
c      do jj=1,mulchLayer
c        do kk=1,SurNodeIndex_M-1 
c           dummy_total_N_new=dummy_total_N_new+CARB_N_mass(jj,kk)
c     &      +CELL_N_mass(jj,kk)+LIGN_N_mass(jj,kk)
c        enddo
c      enddo
c      dummy_total_N_new=dummy_total_N_old
c     &   -dummy_total_N_new-Nmine_total_decomp-Nhumi_total_decomp
c      dummy_total_N_old=dummy_total_N_new/dummy_total_N_old
c      Write(9999,9906) time,dummy_total_N_new,dummy_total_N_old
cccz ------------------------------------------------------------    

       goto 3117
       
       
c ---------------------------- the final assignment of mulch C and mulch N to soil-----------------------------------------------       
c    by the end of mulch decomposition (this is the last step when the mulch model is kicked out).
c        totalMulchC_final   [g]     total amount of C in the mulch, based on resham's data should be PC_mulch~0.4-0.43
c                                    first record the whole mulch mass, then times PC 
c        totalMulchN_final   [g]     total amount of N in the mulch      
3118   bbbb=0.0D0
cccz determine the total mulch C and mulch N
       totalMulchC_final=0.0D0
       totalMulchN_final=0.0D0
       local_soil_mass=0.0D0
       do jj=1,mulchLayer
        do kk=1,SurNodeIndex_M-1  
          totalMulchC_final=totalMulchC_final+CARB_mass(jj,kk)
     &      +CELL_mass(jj,kk)+LIGN_mass(jj,kk)
          totalMulchN_final=totalMulchN_final+CARB_N_mass(jj,kk)
     &      +CELL_N_mass(jj,kk)+LIGN_N_mass(jj,kk)
          CARB_mass(jj,kk)=0.0D0
          CELL_mass(jj,kk)=0.0D0
          LIGN_mass(jj,kk)=0.0D0
          CARB_N_mass(jj,kk)=0.0D0
          CELL_N_mass(jj,kk)=0.0D0
          LIGN_N_mass(jj,kk)=0.0D0
        enddo
       enddo
       totalMulchC_final=totalMulchC_final*PC_mulch
cccz assign mulch C and N to soil pools - convert to ug/cm3 area 
       do kk=1,SurNodeIndex
          qLeft=SurfNodeSurfIndexH(kk)
          qRight=SurfNodeNodeIndexH(kk)
          local_soil_mass=nodeArea(qRight)*BlkDn(MatNumN(qRight))
c    add the C from mulch to soil surface nodes
          CL(qRight)=CL(qRight)
     &       +real(totalMulchC_final*Width(qLeft)/totalMulchWidth
     &       /nodeArea(qRight)
     &       *1.0D6)                                              ! convert to "ug per g of soil"
          CL_old(qRight)=CL(qRight)
c    add the N from mulch to soil surface nodes
          NL(qRight)=NL(qRight)
     &       +real(totalMulchN_final*Width(qLeft)/totalMulchWidth
     &        /nodeArea(qRight)
     &        *1.0D6)
          NL_old(qRight)=NL(qRight)
c    zero the N transfer variable when mulch is finished. 
          NO3_from_residue(qRight)=0.0D0
       enddo
       goto 3117
c -----------------------------------------------------------------------------------------------------------------------------
          
3117   bbbb=0.0D0
       return
3116   Write(*,*) 'Mulch Information is not Included'      
3110   Call errmes(im,il)
9903   Format (1x,A12,T20,A12,T49,A12,T64,A12,T76,A12,
     &  T95,A12,T109,A12,T129,A12,T149,A12,T169,A12,T189,A12)
9904   Format (1x,F17.6,',',9(F16.4,','),F16.4)   
9905   Format (1x,A12,T20,A12,T49,A12)
9906   Format (1x,F17.6, F17.6, F17.6)      
       return
      END