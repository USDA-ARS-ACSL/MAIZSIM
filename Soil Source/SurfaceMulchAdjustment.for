cccz -------------------------------------------      
c     Radiation Adjustment
c     Wind Speed Calculation
c     Vapor Flux
c     Convection Heat Flux, sensi_heat
c     Latent Heat Flux, latent_heat
cccz -------------------------------------------

cccz -------------------------------------------
c     variable used in mulch related module (some of the variables are defined in public ins)
c     scale: maximal layer 10; maximal horizontal nodes 50
c     note: GEO: geometrical property, PHY: physical property, CON: Control property
c
c     INTEGER TYPE
c         DiffusionRes,               [CON] =1 force to use diffusion for flux calculation
c         Van_Camp_mulch,             [CON] choose water retension curve equations
c         IterMulch,                  [CON] for counting the steps in Picard iteration
c         MaxIter,                    [CON] max number of iteration 
c         TimeShrink,                 [CON] indicator if "local time" was used, or "refine" was used
c         WeatherUpdate,              [CON] when weather module are called, update ambient conditions
c         mulchLayer,                 [GEO] number of mulch layer
c         numMulchNode,               [GEO] number of mulch node
c         numMulchEle,                [GEO] number of mulch element
c         MulchNodeMarker(550,2),     [GEO] mulch node array 1 horizontal position; 2 vertical position (1 is soil sruface, mulchLayer+1 is mulch surface)
c         MulchNodeMatrix(11,50),     [GEO] the first index is the vertical position; the second index is the horizontal position; (1,1) is the node in lower-left corner; report the node index
c         MulchEleMarker(550,6),      [GEO] mulch element array 1~4 mulch node index, 5 horizontal position, 6 vertical position (1 is soil surface, mulchLayer+1 is mulch surface)
c         MulchEleMatrix(11,50),      [GEO] the first index is the vertical position; the second index is the horizontal position; (1,1) is the node in lower-left corner; report the element index
c         SurNodeIndex,               [GEO] the number of nodes along soil surface
c         e,                          [GEO] temporery variable for element iteration
c         qLeft, qRight,              [GEO] indices to storge left/right surface nodes
c         SurfNodeNodeIndexH,         [GEO] index of surface node in NodeArray, from one side to another side
c         SurfNodeSurfIndexH,         [GEO] index of surface node in SurfaceNodeArray, from one side to another side
c         
c     FLOAT TYPE
c         maxhNewDiff,                [CON] used for establishing converging test in Picard iteration 
c         maxTmprDiff,                [CON] used for establishing converging test in Picard iteration
c         maxhNew,                    [CON] used for establishing converging test in Picard iteration
c         maxTmpr,                    [CON] used for establishing converging test in Picard iteration
c         errTol,                     [CON] error tolerance for Picard iteration
c         mulchVirtArea123,           [GEO] temporary vairable for ele-area calculation
c         hhhh,                       [GEO] temperary vairable for node vertical coordinate calculation
c         longEmission,               [PHY] temporary vairable for longwave radiation calculation
c         longEmissionAir,            [PHY] temporary vairable for longwave radiation calculation from atmosphere
c         longEmissionSoil,           [PHY] temporary vairable for longwave radiation calculation from soil surface
c         mulchThick,                 [GEO] mulch thickness (cm)
c         f_mulch_pore,               [GEO] mulch "air" space (%)
c         widthPerMulchUnit(50),      [GEO] horizontal width of each mulch section (cm)
c         thickPerLayer(10),          [GEO] thick per mulch layer (cm)
c         mulchNodeCoord(550,2),      [GEO] the coordinate of nodes along mulch-soil interface, 1 horizontal dir 2 vertical dir (cm)
c         MulchNodeMarkerArea(550),   [GEO] nodal area for mulch grid (cm^2) 
c         MulchEleMarkerArea(550),    [GEO] element area for mulch grid (cm^2)
c         slopeCoord(50,2),           [GEO] coordinates, 1 horizontal (x), 2 vertical (y) of nodes along soil surface-mulch interface
c         MulchNodeTmpr(550),         [PHY] nodal temperature
c         MulchNodehNew(550),         [PHY] nodal water potential
c         MulchNodeThNew(550),        [PHY] nodal water content
c         MulchEleTmpr(11,50),        [PHY] element-based temperature, 1 soil surface, mulchLayer for mulch-atmosphere interface
c         MulchElehNew(11,50),        [PHY] element-based water potential, 1 soil surface, mulchLayer for mulch-atmosphere interface
c         MulchEleThNew(11,50),       [PHY] element-based water content, 1 soil surface, mulchLayer for mulch-atmosphere interface
c         MulchElehNew_temp(11,50),   [PHY] element-based water potential for iteration 
c         MulchEleTmpr_temp(11,50),   [PHY] element-based temperature for iteration 
c         MulchEleThNew_temp(11,50),  [PHY] element-based water content for iteration 
c         MulchElehNew_temp0(11,50),  [PHY] element-based water potential for iteration   --- store local step data
c         MulchEleTmpr_temp0(11,50),  [PHY] element-based temperature for iteration       --- store local step data
c         MulchElehNew_temp2(11,50),  [PHY] element-based water potential for iteration   --- first to receive any updates
c         MulchEleTmpr_temp2(11,50),  [PHY] element-based temperature for iteration       --- first to receive any updates
c         MulchElehNew_temp3(11,50),  [PHY] element-based water potential for iteration   --- static var for multi-stage average
c         MulchEleTmpr_temp3(11,50),  [PHY] element-based temperature for iteration       --- static var for multi-stage average
c         shortRadDir(11),            [PHY] short wave downwards radiaiton at each nodal layer (mulchLayer+1 layers total) (W/m^2)
c         shortRadFirst(11,50),       [PHY] short wave upwards radiaiton (first order reflection) at each nodal layer (mulchLayer+1 layers total) (W/m^2)
c         longRadDown(11,50),         [PHY] long wave downwards radiaiton at each nodal layer (mulchLayer+1 layers total) (W/m^2)
c         longRadUp(11,50),           [PHY] long wave upwards radiaiton at each nodal layer (mulchLayer+1 layers total) (W/m^2)
c         netRad(11,50),              [PHY] net radiation at each nodal layer (mulchLayer+1 layers total) (W/m^2)
c         netRadEachLayer(11,50),     [PHY] net radiation of each mulch layer (mulchLayer layers total) [w/m^2 * (cm/100) = W/m]
c         RelaHumid_mulch(11,50),     [PHY] relative Humidity within each mulch element (%)
c         VaporSat_mulch(11,50),      [PHY] saturated vapor density for each mulch element (g/m^3)
c         VaporAct_mulch(11,50),      [PHY] actual vapor density for each mulch element (g/m^3)
c         VaporDiff_mulch_temp(11,50),[PHY] vapor diffusivity for each mulch element (g/m^2/day)
c         ParMulch(13,1),             [PHY] water characteristic curve for solid portion in mulching (plant tissue)
c         WaterCapaMulch_temp(11,50), [PHY] water capacity of each mulch element during iteration, i.e., the solid part values time (1-f_mulch_pore) (1/cm)
c         HeatCapa_Air_temp(11,50),   [PHY] heat capacity of mulch (air) of each mulch element during iteration (J/m^3/K)
c         HeatDiff_mulch_temp(11,50), [PHY] heat diffusion of mulch (air) of each mulch element during iteration (J/m^2/K/day)
c         RadSolarShort_Wea,          [PHY] incoming short wave radiation updated from weather module (W/m^2)
c         AirTemp_Wea,                [PHY] atmosphere tmperature (oC)
c         AirVaporP_Wea,              [PHY] atmosphere vapor pressure (kPa)
c         AirVPD_Wea,                 [PHY] atmosphere vapor deficit (kPa)
c         CloudCoverFactor_Wea,       [PHY] atmosphere cloud cover (%)
c         AirWind_Wea,                [PHY] atmosphere wind speed (km/hour)
c         RelaHumi_Wea,               [PHY] atmosphere relative humidity (%)
c         VaporSat_ambient,           [PHY] atmosphere saturated vapor density (g/m^3)
c         VaporAct_ambient,           [PHY] atmosphere actual vapor density (g/m^3)
c         HeatCapa_ambient,           [PHY] atmosphere (air) heat capacity (J/m^3/K)
c         Tmpr_Sur,                   [PHY] soil surface temperature
c         hNew_Sur,                   [PHY] soil surface water potential
c         VaporSat_Sur,               [PHY] soil surface saturated water vapor density
c         RelaHumid_Sur,              [PHY] soil surface saturated related humidity
c         VaporAct_Sur                [PHY] soil surface actural vapor density: VaporSat_Sur*RelaHumid_Sur
c         SurTemp_Compare(50),        [PHY] store the soil surface temperature after "SurfaceAdjPrior" for future comparison
c         SurhNew_Compare(50),        [PHY] store the soil surface water potential after "SurfaceAdjPrior" for future comparison
c         g_Heat_Sensi(30,50),        [PHY] sensible heat flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards flux (J/m^2/day)
c         g_Heat_Latent(30,50),       [PHY] latent heat flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards flux (J/m^2/day)
c         g_Vapor(30,50),             [PHY] water vapor flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards flux (g/m^2/day)
c         DeltaRshort,                [PHY] residue-area index of each mulch layer, Novak et al., (2000) (1)
c         Omega,                      [PHY] clumping index, which describes the way in which the residue elements are arranged, Novak et al., (2000) (1)
c         epsilon_m,                  [PHY] emissivity of mulch
c         alpha_m,                    [PHY] short wave reflectivity of mulch
c         epsilon_s,                  [PHY] emissivity of soil
c         alpha_s,                    [PHY] short wave reflectivity of soil
c         epsilon_a_0,                [PHY] ambient air emissivity
c         coef_epsilon_a,             [PHY] correction factors for ambient air emissivity
c         sigma,                      [PHY] Stefan-Boltzmann constant
c         co_disp_heig,               [PHY] coefficient of displacement height for wind speed (1)
c         co_rough_len,               [PHY] coefficient of roughness length for wind speed (1)
c         rho_mulch,                  [PHY] mulch density, assumed to be the density of wood (g/m^3)
c         rho_dryair,                 [PHY] density of dry air (g/m^3)
c         rho_w,                      [PHY] water density (g/m^3)
c         Tmpr_ref,                   [PHY] reference temperature for heat equation, especially the latent heat fluxes 
c         c_vap_s,                    [PHY] water vapor specific heat (J/g/K)
c         c_vap_v,                    [PHY] vaporization heat (J/g)
c         c_air_s,                    [PHY] specific heat of air, may need to change to enthalpy (J/g/K)
c         c_water_s                   [PHY] specific heat of water (J/g/K)
c         Dhm,                        [PHY] molecular diffusivity for sensible heat
c         nu_air,                     [PHY] air kinematic viscosity
c         Ra_Critical,                [PHY] critical number for Rayleigh number
c         GrRe_Critical               [PHY] critical number for Grashof/Reynold^2
c         LocalStep,                  [PHY] local time step for surface process modules
c         TotalTime,                  [PHY] total time run for the current module, i.e., the time step for all "other" modules
c         u_mulch(11),                [PHY] the wind speed (km/h) at mulch interior interface
c         u_above,                    [PHY] the wind speed (km/h) at mulch-ambient interface
c         u_soilsur,                  [PHY] the wind speed (km/h) at mulch-soil interface
c         RainFallInput(11,50),       [PHY] rainfall input for each mulch layer
c         inputPerLayer,              [PHY] iterative variable for rainfall
c         PriorStep,                  [GEO] save the previous time step
c         thresholdThick,             [GEO] threshold for merging two neighbor layers
c         thickPerLayer_temp(11),     [GEO] temporary array for thickness per layer
c         mergeIndex(11),             [GEO] input current layer index, output index after layer merging
c         TotalTime,                  [CON] control local time-step
c         LocalStep,                  [CON] control local time-step
c         Local_VarBW1,               [PHY] local variable for varBW, need average
c         Local_VarBW2,               [PHY] local variable for varBW, need average
c         Local_VarBW3,               [PHY] local variable for varBW, need average
c         Local_VarBT1,               [PHY] local variable for varBT, need average
c         Local_VarBT2,               [PHY] local variable for varBT, need average 
c         Local_VarBT3,               [PHY] local variable for varBT, need average
c         Local_VarBT4,               [PHY] local variable for varBT, need average
c         VarBW1_temp(11),            [PHY] local variable for varBW, after average
c         VarBW2_temp(11),            [PHY] local variable for varBW, after average
c         VarBW3_temp(11),            [PHY] local variable for varBW, after average
c         Q_temp(11),                 [PHY] local variable for Q, after average
c         VarBT1_temp(11),            [PHY] local variable for varBT, after average
c         VarBT2_temp(11),            [PHY] local variable for varBT, after average
c         VarBT3_temp(11),            [PHY] local variable for varBT, after average
c         VarBT4_temp(11),            [PHY] local variable for varBT, after average
c         OutputTime_Mulch,           [CON] output time (OUTPUT.for)
c         OutputTimeStep_Mulch        [CON] output time step (OUTPUT.for)
c         RDM_mass,                   [PHY] mass of Carbon in Rapidly decomposable material (RDM) for each element
c         HCE_mass,                   [PHY] mass of Carbon in Hemi-celluloses (HCE) for each element
c         CEL_mass,                   [PHY] mass of Carbon in Cellulose (CEL) for each element
c         LIG_mass,                   [PHY] mass of Carbon in Lignin (LIG) for each element
c         hhhh,                       [...] dummy variable to make some height
c         bbbb,                       [...] dummy variable
c         aaaa,                       [...] dummy variable to receive unimportant values/mark locations
        
       Subroutine Surface_Mulch_Process()
      
       include 'public.ins'
       include 'puplant.ins'
       include 'puweath.ins'
       include 'pusurface.ins'
      
       double precision aaaa
       double precision a11,a12,a21,a22,b11,b12
cccz parameter that connect surface water and mulch
       double precision h_Pond_max,DiffSubmerge,PerOccupation
       integer  SubmergeIndex
cccz parameter for control
       integer  DiffusionRes
cccz geometrical parameter for mulch nodes/elements
c       integer  numMulchNode,numMulchEle
       integer  e,IterMulch,qLeft,qRight
      ! temporary geometrical parameter (area and cumulated elevation calculation)
       double precision mulchVirtArea123,hhhh
cccz physical properties within mulch
       double precision                                           !cccz for mulch water release curve
     !   DeltaRshort,DeltaRlong,Omega,epsilon_m,alpha_m,alpha_s,
     !   epsilon_s,epsilon_a_0,coef_epsilon_a,                    !cccz radiation related parameters
     !   rho_dryair,rho_w,                                        !cccz density group (mulch density in head file)
     !   co_disp_heig,co_rough_len,u_above,u_soilsur,u_mulch(11), !cccz wind speed group
     !   c_vap_s,c_vap_v,c_air_s,c_water_s,Tmpr_ref,              !cccz heat group
     !   Dhm,nu_air,Ra_Critical,GrRe_Critical                     !cccz other variables
cccz physical process within mulch
       double precision
     !   netRad(11,50),netRadEachLayer(11,50),shortRadDir(11),
     !   shortRadFirst(11,50),longRadDown(11,50),longRadUp(11,50),     
     !   longEmission,longEmissionAir,longEmissionSoil,            !cccz temporary parameters for radiation calculation
     !   g_Heat_Sensi(30,50),g_Heat_Latent(30,50),g_Vapor(30,50)
cccz physical status of mulch
       double precision 
     !   MulchEleTmpr_temp(11,50),MulchElehNew_temp(11,50),
     !   MulchEleThNew_temp(11,50),                           !cccz for exterior Picard iteration
     !   MulchElehNew_temp2(11,50),MulchEleTmpr_temp2(11,50), !cccz for interior Picard iteration
     !   MulchElehNew_temp3(11,50),MulchEleTmpr_temp3(11,50), !cccz for cumulative average
     !   MulchElehNew_temp0(11,50),MulchEleTmpr_temp0(11,50)  !cccz for time domain subdivision
cccz parameters for converging tests
       integer  MaxIter,TimeShrink
       double precision maxhNewDiff,maxTmprDiff,maxhNew,maxTmpr,errTol
cccz parameters for weather/soil conditions
       double precision RainFallInput(11,50),inputPerLayer,            !ccczparameter for rainfall redistribution 
     !    VaporSat_ambient,VaporAct_ambient,HeatCapa_ambient,
     !    Tmpr_Sur,hNew_Sur,VaporSat_Sur,RelaHumid_Sur,VaporAct_Sur   !cccz derived surface values 
cccz parameters for shrinking the grid
       double precision PriorStep,          ! record the previous time step
     !    thresholdThick,thickPerLayer_temp(11)
       integer mergeIndex(11)
cccz some auxillary parameters
       double precision Phi_m,             ! the momentum correction factor (not included in this version yet, need complicated calculation later)
     !    karman,z_rough,d_displace,z_wind,hMulchInit,sigma
cccz some weather parameters
       double precision RadSolarShort_Wea,AirTemp_Wea,AirVaporP_Wea,
     !    AirVPD_Wea,CloudCoverFactor_Wea,AirWind_Wea,RelaHumi_Wea
cccz data exchange between mulch module and other modles
       double precision  Local_VarBW1, Local_VarBW2, Local_VarBW3,
     !    Local_VarBT1, Local_VarBT2, Local_VarBT3, Local_VarBT4,
     !    TotalTime, LocalStep
cccz data for output
cccz moved to output now
cccz      double precision OutputTime_Mulch, OutputTimeStep_Mulch
       double precision 
     !    RelaHumid_mulch_temp(11,50),VaporSat_mulch_temp(11,50),
     !    VaporDiff_mulch_temp(11,50),VaporAct_mulch_temp(11,50),
     !    HeatCapa_Air_temp(11,50),HeatDiff_mulch_temp(11,50),
     !    WaterCapaMulch_temp(11,50),Ra_mulch(11,50),GrRe_mulch(11,50),
     !    VarBW1_temp(11),VarBW2_temp(11),VarBW3_temp(11),Q_temp(11),
     !    VarBT1_temp(11),VarBT2_temp(11),VarBT3_temp(11),
     !    VarBT4_temp(11),RainFallInput_temp(11,50),
     !    WaterDifici_mulch_temp(11,50)
       integer mulch_np,mulch_nq, nNode,iday
cccz the organic matter decomposition parameters
       double precision DCarb,DCell,DLign,PerCarb,PerCell,PerLign,
     !    RecycleFactor,Fon,Fom,Fom_PctN,mulchBulkDensity,
     !    OrganC,InorganN
c      Character*10 Date
cccz the organic matter mass tracer
       double precision RDM_mass_temp(11,50),HCE_mass_temp(11,50),
     !    CEL_mass_temp(11,50),LIG_mass_temp(11,50)
cccz static variables
       Common /SurfaceAdj_Static/ DiffusionRes,
c     !    numMulchNode,numMulchEle,
     !    DeltaRshort,DeltaRlong,Omega,epsilon_m,alpha_m,alpha_s,
     !    epsilon_s,epsilon_a_0,coef_epsilon_a, 
     !    rho_dryair,rho_w,
     !    co_disp_heig,co_rough_len,u_above,u_soilsur,u_mulch,
     !    c_vap_s,c_vap_v,c_air_s,c_water_s,Tmpr_ref,
     !    Dhm,nu_air,Ra_Critical,GrRe_Critical,
     !    netRad,netRadEachLayer,shortRadDir,
     !    shortRadFirst,longRadDown,longRadUp,
     !    g_Heat_Sensi,g_Heat_Latent,g_Vapor,
     !    MulchEleTmpr_temp,MulchElehNew_temp,MulchEleThNew_temp,
     !    MulchElehNew_temp3,MulchEleTmpr_temp3,
     !    MaxIter,TimeShrink,errTol,
     !    RainFallInput,inputPerLayer,
     !    VaporSat_ambient,VaporAct_ambient,HeatCapa_ambient,
     !    PriorStep,thresholdThick,mergeIndex,
     !    RelaHumid_mulch_temp,VaporSat_mulch_temp,VaporDiff_mulch_temp,
     !    VaporAct_mulch_temp,HeatDiff_mulch_temp,WaterCapaMulch_temp,
     !    WaterDifici_mulch_temp,
     !    RadSolarShort_Wea,AirTemp_Wea,AirVaporP_Wea,
     !    CloudCoverFactor_Wea,AirWind_Wea,RelaHumi_Wea
cccz initial input values in case of reading failures
       Data mulchThick/10.0D0/,f_mulch_pore/1.0D0/,mulchLayer/5/,     
     &    DeltaRshort/0.28D0/,Omega/0.60D0/,epsilon_s/1.0D0/,
     &    epsilon_m/1.0D0/,epsilon_a_0/0.60D0/,coef_epsilon_a/5.95D-5/,
     &    alpha_m/0.46D0/,sigma/5.67D-8/,            
     &    co_disp_heig/0.87D0/,co_rough_len/0.079D0/,
     &    rho_mulch/800000.0D0/,rho_w/1000000.0D0/,rho_dryair/1225.0D0/,
     &    Tmpr_ref/20.0D0/,c_vap_s/1.862D0/,c_vap_v/2453.5D0/,
     &    Dhm/2.2D-5/,nu_air/1.5D-5/,c_air_s/0.718D0/,
     &    c_water_s/4.186D0/,
     &    Ra_Critical/1706.0D0/,GrRe_Critical/0.001D0/

       If (lInput.eq.1) then
                  
cccz first read the mulching information from file
       Open(20,file=MulchFile,status='old',ERR=2106)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100) MaxIter, errTol,DiffusionRes
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) DeltaRshort, DeltaRlong, 
     &      Omega, epsilon_m, alpha_m
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) epsilon_s, epsilon_a_0, coef_epsilon_a
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) mulchThick, mulchLayer, f_mulch_pore
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) rho_mulch, co_disp_heig, co_rough_len
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) rho_dryair
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) Tmpr_ref
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
CDT make parameter      
       Read(20,*,ERR=2100) c_vap_s, c_vap_v, c_air_s, c_water_s
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
       Read(20,*,ERR=2100) Van_Camp_mulch
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       i1=i1+1
cccz input selected moisture release curve      
       if(Van_Camp_mulch.eq.0) then                   ! van Genuchten Equation
        Read(20,*,ERR=2100) (ParMulch(i,1),i=1,13)      
       elseif(Van_Camp_mulch.eq.1) then               ! Campbell Equation
        Read(20,*,ERR=2100) (ParMulch(i,1),i=1,4)
       else
       endif
cccz       
       close(20)
                    
cccz  
cccz  transfer weather data for mulch
cccz  get short radiation, AirTemp at current time
      RadSolarShort_Wea=WATTSM(ITIME)
      AirTemp_Wea=TAIR(ITIME)
      AirVPD_Wea=VPD(ITIME)
      AirVaporP_Wea=
     &   0.61*EXP((17.27*AirTemp_Wea)/(AirTemp_Wea + 237.3))-AirVPD_Wea
      CloudCoverFactor_Wea=CLOUD
      AirWind_Wea=WIND             ! unit km/hour
      RelaHumi_Wea=AirVaporP_Wea/
     &   (0.61*EXP((17.27*AirTemp_Wea)/(AirTemp_Wea + 237.3)))
      
cccz start to establish the mulch geo configurations
       LayerHeight(1)=0.0D0
       do n=1,mulchLayer
        thickPerLayer(n)=mulchThick/mulchLayer     ! thickness per layer
        LayerHeight(n+1)=LayerHeight(n)+thickPerLayer(n)
       enddo
       thresholdThick=0.25D0*thickPerLayer(1)

cccz -----------------------------------------------------------------------------------------
c we keep this in both runoff and surface mulch in case

cccz here we start to sort the surface, 'SurNodeIndex' eventually = num of surface node,
cccz we keep it here because we want the surface computing system is relatively independent to other 2Dsoil module
cccz start with 'SurNodeIndex=0' 
       SurNodeIndex=0
cccz here we loop through the boundary node
       do i=1,NumBP
        n=KXB(i)
        k=CodeW(n)
cccz if surface condition (\pm 4) is found, then record
        If(K.eq.4.or.K.eq.-4) then
         SurNodeIndex=SurNodeIndex+1
         SurfNodeNodeIndexH(SurNodeIndex)=n
        endif
       enddo
          
cccz now we need to sort the surface node
cccz we have a basic assumption that the first node is always @ upper-left corner
cccz but we do not made additional assumption on the order of other node,
cccz or the shape of the surface.

cccz some aux variables, the meaning is straightforward from the following loop
      InElement=0       
      Dist=0.0D0
      do i=1,SurNodeIndex
       Dist=0.0D0
       do j=1,NumEL
        do kk=1,4
          if(SurfNodeNodeIndexH(i).eq.KX(j,kk)) InElement=1
        enddo
cccz the idea for the following loop
cccz choose the only one nearest point, next to the current step one
cccz 
        if(InElement.eq.1) then
         do kk=1,4
          if(abs(CodeW(KX(j,kk))).eq.4) then
           if(i.eq.1) then
            if(KX(j,kk).ne.SurfNodeNodeIndexH(i)) then
             SurfNodeNodeIndexH(i+1)=KX(j,kk)
            endif
           else
            if(KX(j,kk).ne.SurfNodeNodeIndexH(i).
     &        and.KX(j,kk).ne.SurfNodeNodeIndexH(i-1)) then
              if(Dist.le.0.0D0) then
                Dist=abs(X(SurfNodeNodeIndexH(i))-X(KX(j,kk)))
                SurfNodeNodeIndexH(i+1)=KX(j,kk) 
              else
              if(abs(X(SurfNodeNodeIndexH(i))-X(KX(j,kk))).lt.Dist) then
                 Dist=abs(X(SurfNodeNodeIndexH(i))-X(KX(j,kk)))                          
                 SurfNodeNodeIndexH(i+1)=KX(j,kk) 
              endif
              endif
            endif
           endif
          endif
         enddo
        endif
        InElement=0
       enddo
      enddo
cccz here we record the (x,y)-coordinate of surface node
cccz in the order given by SurfNodeNodeIndexH, i.e., from left to right
cccz we also find the peak elevation
      PEAK=0.0D0
      do i=1,SurNodeIndex
        slopeCoord(i,1)=X(SurfNodeNodeIndexH(i))
        slopeCoord(i,2)=Y(SurfNodeNodeIndexH(i))
        if(PEAK.lt.slopeCoord(i,2)) then
           PEAK=slopeCoord(i,2)
        endif
      enddo
cccz we also make an order based on surface node index
cccz SurfNodeNodeIndexH/SurfNodeSurfIndexH are important geometric properties
cccz which will be used for all the new developed surface modules
      do n=1,SurNodeIndex
       do i=1,NumBP
         if(SurfNodeNodeIndexH(n).eq.KXB(i)) then
           SurfNodeSurfIndexH(n)=i
         endif
       enddo
      enddo
      
cccz ---------------------------------------------------------------------------

      do n=1,SurNodeIndex-1
        widthPerMulchUnit(n)=abs(slopeCoord(n+1,1)-slopeCoord(n,1))   
      enddo

cccz Find appropriate initial values for mulch nodes    
CDT initial values of what - hydraulic properties?
cccz initialize the functions 
cccz (the parameters in those functions are initialized in the tests of linput)
      hMulchInit=4708.34749D0*(AirTemp_Wea+273.15D0)
     &   *log(100.*RelaHumi_Wea) ! 4708.34749 \approx 8.314D0/0.018D0/9.81D0*100.0D0
cccz we just run the following function one time
cccz to initialize the parameter within them
      if(Van_Camp_mulch.eq.0) then
        thMulchIni=WQ_VANG_MULCH(hMulchInit,ParMulch(:,1),lInput)
        thMulchIni=thMulchIni*(1.0D0-f_mulch_pore)                     ! the mulch theta value is for the whole mulch, not only the solid part. the converting parameter is 'f_mulch_pore'
        aaaa=WH_VANG_MULCH
     &       (thMulchIni/(1.0D0-f_mulch_pore),ParMulch(:,1),lInput)
        waterCapaIni=WC_VANG_MULCH(hMulchInit,ParMulch(:,1),lInput)
        waterCapaIni=waterCapaIni*(1.0D0-f_mulch_pore)                 ! the mulch water capa value is for the whole mulch, not only the solid part. the converting parameter is 'f_mulch_pore'
      elseif(Van_Camp_mulch.eq.1) then
        thMulchIni=WQ_CAMP_MULCH(hMulchInit,ParMulch(:,1),lInput)
        thMulchIni=thMulchIni*(1.0D0-f_mulch_pore)
        aaaa=WH_CAMP_MULCH
     &       (thMulchIni/(1.0D0-f_mulch_pore),ParMulch(:,1),lInput)
        waterCapaIni=WC_CAMP_MULCH(hMulchInit,ParMulch(:,1),lInput)
        waterCapaIni=waterCapaIni*(1.0D0-f_mulch_pore)
      else
      endif
cccz ----------------------------------------------------------------------------
cccz Establish the mulching grid (Node)     
      numMulchNode=0
      hhhh=0.0D0
      do i=1,mulchLayer+1
       if (i.ge.2) then
        hhhh=hhhh+thickPerLayer(i-1)
       endif
       do n=1,SurNodeIndex
        numMulchNode=numMulchNode+1
        if(i.eq.1) then
cccz Here we copy the soil surface node (establish the interface)
         mulchNodeCoord(numMulchNode,1)=slopeCoord(n,1)
         mulchNodeCoord(numMulchNode,2)=slopeCoord(n,2)
         MulchNodeMarker(numMulchNode,1)=n
         MulchNodeMarker(numMulchNode,2)=1
        else
cccz Stack upper nodes above soil surface, within mulch
cccz we always assmue the equal-distance layer at beginning
         mulchNodeCoord(numMulchNode,1)=slopeCoord(n,1)
         mulchNodeCoord(numMulchNode,2)=slopeCoord(n,2)+hhhh
         MulchNodeMarker(numMulchNode,1)=n
         MulchNodeMarker(numMulchNode,2)=i
        endif
        MulchNodeMatrix(i,n)=numMulchNode
        MulchNodeMarkerArea(numMulchNode)=0.0D0
        MulchNodeTmpr(numMulchNode)=AirTemp_Wea
        MulchNodehNew(numMulchNode)=hMulchInit        ! the initial water potential
        MulchNodeThNew(numMulchNode)=thMulchIni       ! the initial water content
       enddo
      enddo
cccz ---------------------------------------------------------------------          
cccz Establish the mulching grid (Element)
      numMulchEle=0
      do i=1,mulchLayer
       do n=1,SurNodeIndex-1
        numMulchEle=numMulchEle+1
        MulchEleMarker(numMulchEle,1)=MulchNodeMatrix(i,n)
        MulchEleMarker(numMulchEle,2)=MulchNodeMatrix(i,n+1)
        MulchEleMarker(numMulchEle,3)=MulchNodeMatrix(i+1,n+1)
        MulchEleMarker(numMulchEle,4)=MulchNodeMatrix(i+1,n)
        MulchEleMarker(numMulchEle,5)=n
        MulchEleMarker(numMulchEle,6)=i
        MulchEleMatrix(i,n)=numMulchEle
        MulchEleMarkerArea(numMulchEle)=0.0D0
       enddo
      enddo
cccz Calculate the area and geometry of each node and element    
      do e=1,numMulchEle
        n1=MulchEleMarker(e,1)
        n2=MulchEleMarker(e,2)
        n3=MulchEleMarker(e,3)
        n4=MulchEleMarker(e,4)
        x1=mulchNodeCoord(n1,1)
        x2=mulchNodeCoord(n2,1)
        x3=mulchNodeCoord(n3,1)
        x4=mulchNodeCoord(n4,1)
        y1=mulchNodeCoord(n1,2)
        y2=mulchNodeCoord(n2,2)
        y3=mulchNodeCoord(n3,2)
        y4=mulchNodeCoord(n4,2)                  
        if(x3.eq.x4.and.y3.eq.y4) then
cccz triangle shape element
cccz triangle element is hard to be implement here
cccz so far we decided not to implement triangular element
        else
cccz rectangular shape element
         CJ=x2-x1
         CK=x3-x1
         BJ=Y2-Y1
         BK=Y3-Y1   
         mulchVirtArea123=0.5D0*abs(CJ*BK-CK*BJ)                 
         CJ=x3-x1
         CK=x4-x1
         BJ=Y3-Y1
         BK=Y4-Y1
         mulchVirtArea123=mulchVirtArea123+0.5D0*abs(CJ*BK-CK*BJ)
         MulchEleMarkerArea(e)=mulchVirtArea123
         MulchNodeMarkerArea(n1)=MulchNodeMarkerArea(n1)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n2)=MulchNodeMarkerArea(n2)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n3)=MulchNodeMarkerArea(n3)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n4)=MulchNodeMarkerArea(n4)
     &    +0.25D0*MulchEleMarkerArea(e)
        endif                         
      enddo
cccz initial mulch tmpr, hNew and ThNew
      do e=1,numMulchEle
       n1=MulchEleMarker(e,1)
       n2=MulchEleMarker(e,2)
       n3=MulchEleMarker(e,3)
       n4=MulchEleMarker(e,4)
       kk=MulchEleMarker(e,6)
       jj=MulchEleMarker(e,5)
       MulchEleTmpr(kk,jj)=0.25D0     ! here we exchange the horizontal/vertical coordinate
     &  *(MulchNodeTmpr(n1)+MulchNodeTmpr(n2)
     &  +MulchNodeTmpr(n3)+MulchNodeTmpr(n4))
       MulchEleThNew(kk,jj)=0.25D0
     &  *(MulchNodeThNew(n1)+MulchNodeThNew(n2)
     &  +MulchNodeThNew(n3)+MulchNodeThNew(n4))
      enddo
      do kk=1,mulchLayer
       do jj=1,SurNodeIndex-1
       if(Van_Camp_mulch.eq.0) then
        MulchElehNew(kk,jj)=
     &   WH_VANG_MULCH(MulchEleThNew(kk,jj)/(1-f_mulch_pore),
     &   ParMulch(:,1),lInput)    
       elseif(Van_Camp_mulch.eq.1) then
        MulchElehNew(kk,jj)=
     &   WH_CAMP_MULCH(MulchEleThNew(kk,jj)/(1-f_mulch_pore),
     &   ParMulch(:,1),lInput) 
       else 
       endif
       enddo
      enddo
      
cccz -----------------------------------------------------------------------
cccz initialize the mulch output options
cccz we temporarily put it here, but will move it to the OUTPUT.for
c       OutputTimeStep_Mulch=1.0D0/24.0D0/60.0D0
c       OutputTime_Mulch=time+OutputTimeStep_Mulch
c       Open(2102,file=MassBalanceMulchFileOut)
c       Write(2102,2103) 'Date_time','Date','X','Y','hNew','thNew','Temp'
c       iday=int(time)
c       call caldat(iday,mm,id,iyyy) 
c       write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
c       do n=1, numMulchNode
c         Write(2102,2104) Time, Date, mulchNodeCoord(n,1),
c     &     mulchNodeCoord(n,2),MulchNodehNew(n),MulchNodeThNew(n),
c     &     MulchNodeTmpr(n)
c       enddo
cccz -----------------------------------------------------------------------
  
       VaporSat_ambient=exp(19.84D0-4975.9D0/(AirTemp_Wea+273.15D0))  ! Saturated Vapor density of ambient air (g/m^3)
       VaporAct_ambient=VaporSat_ambient*RelaHumi_Wea                 ! Actual Vapor density of ambient air (g/m^3)
       HeatCapa_ambient=1.85D0*VaporAct_ambient
     &   +1.006D0*(rho_dryair-VaporAct_ambient)                        ! Vapor Heat Capacity of ambient air (J/m^3/K)
          
cccz Calculate the shortwave downwards radiation in each interface among between adjacent mulch layers
cccz Based on the model in 'Simulating the radiation distribution within a barley-straw mulch' Novak 2000
cccz use random distribution Eq(3) for transmittivity

cccz need to consider the current surface ponded water heights 'h_Pond'
       h_Pond_max=0.0D0
       do n=1,SurNodeIndex
          h_Pond_max=max(h_Pond_max,h_Pond(n))
       enddo

cccz determine if the ponded water reach which 'element-based layer'
cccz in another word, node of layer 'SubmergeIndex' under the water surface
       SubmergeIndex=0
       DiffSubmerge=0
       if(h_Pond_max.le.0) then
         SubmergeIndex=0
       else
        do n=1,mulchLayer+1
         if(h_Pond_max.ge.LayerHeight(n)) then
           SubmergeIndex=n
           if(n.le.mulchLayer) then 
             DiffSubmerge=LayerHeight(n+1)-h_Pond_max
           endif
          endif
        enddo
       endif
       if(SubmergeIndex.eq.mulchLayer+1) then
         PerOccupation=3.0D0
       else
        if(DiffSubmerge.le.thresholdThick) then
          PerOccupation=2.0D0
        else
          PerOccupation=DiffSubmerge/thickPerLayer(SubmergeIndex)
        endif
       endif
        
       if(SubmergeIndex.le.0) then
cccz in the case there is no ponded water
         do kk=mulchLayer,1,-1
          if(kk.eq.mulchLayer) then
            shortRadDir(kk+1)=RadSolarShort_Wea
          elseif(kk.eq.(mulchLayer-1)) then
            shortRadDir(kk+1)=RadSolarShort_Wea*(1.0D0-DeltaRshort)
          else
          shortRadDir(kk+1)=shortRadDir(kk+2)*(1.0D0-Omega*DeltaRshort)
          endif
         enddo
         shortRadDir(1)=shortRadDir(2)*(1.0D0-Omega*DeltaRshort)
       else
cccz elimiated the short direct raditaiton to the node under water 
         do kk=mulchLayer,SubmergeIndex,-1
          if(kk.eq.mulchLayer) then
            shortRadDir(kk+1)=RadSolarShort_Wea
          elseif(kk.eq.(mulchLayer-1)) then
            shortRadDir(kk+1)=RadSolarShort_Wea*(1.0D0-DeltaRshort)
          else
          shortRadDir(kk+1)=shortRadDir(kk+2)*(1.0D0-Omega*DeltaRshort)
          endif
         enddo
         shortRadDir(SubmergeIndex)=
     &      shortRadDir(SubmergeIndex+1)*(1.0D0-Omega*DeltaRshort)
         do kk=1,SubmergeIndex-1
           shortRadDir(kk)=0.0D0
         enddo
       endif

       call MulchDecomp()
       
       return
      Endif
      
cccz **************************************************************************
cccz Non-initialization section ('lInput.ne.1') starts the real simulation
      
        
cccz  
cccz  transfer weather data for mulch
cccz  get short radiation, AirTemp at current time
      RadSolarShort_Wea=WATTSM(ITIME)
      AirTemp_Wea=TAIR(ITIME)
      AirVaporP_Wea=
     &   0.61*EXP((17.27*AirTemp_Wea)/(AirTemp_Wea + 237.3))-AirVPD_Wea
      CloudCoverFactor_Wea=CLOUD
      AirWind_Wea=WIND             ! unit km/hour
      RelaHumi_Wea=AirVaporP_Wea/
     &   (0.61*EXP((17.27*AirTemp_Wea)/(AirTemp_Wea + 237.3)))

      VaporSat_ambient=exp(19.84D0-4975.9D0
     &    /(AirTemp_Wea+273.15D0))                                    ! Saturated Vapor density of ambient air (g/m^3)
      VaporAct_ambient=VaporSat_ambient*RelaHumi_Wea                ! Actual Vapor density of ambient air (g/m^3)
      HeatCapa_ambient=1.85D0*VaporAct_ambient
     &    +1.006D0*(rho_dryair-VaporAct_ambient)                      ! Vapor Heat Capacity of ambient air (J/m^3/K)
cccz       Endif   
cccz ---------------------------------------------------------------------
cccz embed the code for changing grid
cccz adjust the mulch thickness
cccz if layer was shrunk, then there will be some water loss and energy loss (temperature loss)
cccz that propotional to the shrunk size
cccz Such assumption can be changed by actively recalculate the water content,
cccz but the temperature loss should be fine

cccz initialized the fraction factor for mulch shrinking
      do kk=1,mulchLayer
          Frac_Decomp(kk)=0.0D0
      enddo
      
      call MulchDecomp()
      mulchLayer_temp=mulchLayer
      do kk=1,mulchLayer
        thickPerLayer(kk)=thickPerLayer(kk)*(1.0D0-Frac_Decomp(kk))
        mergeIndex(kk)=kk
        thickPerLayer_temp(kk)=thickPerLayer(kk) 
      enddo
      
      kk=1
2201  if(mulchLayer_temp.gt.1.and.kk.lt.mulchLayer_temp) then
       if(thickPerLayer_temp(kk).lt.thresholdThick) then
         thickPerLayer_temp(kk)=thickPerLayer_temp(kk)
     &     +thickPerLayer_temp(kk+1)
         mulchLayer_temp=mulchLayer_temp-1
         mergeIndex(mulchLayer-mulchLayer_temp+1)=kk
         do jj=kk+1,mulchLayer_temp
          thickPerLayer_temp(jj)=thickPerLayer_temp(jj+1)
         enddo
         thickPerLayer_temp(mulchLayer_temp+1)=0.0D0
         do jj=mulchLayer-mulchLayer_temp+2,mulchLayer
          mergeIndex(jj)=mergeIndex(jj)-1
         enddo
         goto 2201
       else
        kk=kk+1
        goto 2201
       endif
      elseif(mulchLayer_temp.gt.1.and.kk.eq.mulchLayer_temp) then
       if(thickPerLayer_temp(kk).lt.thresholdThick) then  
        thickPerLayer_temp(kk-1)=thickPerLayer_temp(kk)
     &     +thickPerLayer_temp(kk-1)
        mulchLayer_temp=mulchLayer_temp-1
        mergeIndex(mulchLayer)=mergeIndex(mulchLayer-1)
        thickPerLayer_temp(kk)=0.0D0
       endif
      else          
      endif             
      PriorStep=Step

cccz Combine the physical parameters within the mulch
cccz initialize the temperatory physical coef: head, water content, temperature
cccz initialize the temperatory mulch coef: four-category of organic matter
      if(mulchLayer_temp.lt.mulchLayer) then
       do kk=1,mulchLayer
        do nn=1,SurNodeIndex-1
         MulchEleThNew_temp(kk,nn)=0.0D0
         MulchEleTmpr_temp(kk,nn)=0.0D0
         MulchElehNew_temp(kk,nn)=0.0D0
         RDM_mass_temp(kk,nn)=0.0D0
         HCE_mass_temp(kk,nn)=0.0D0
         CEL_mass_temp(kk,nn)=0.0D0
         LIG_mass_temp(kk,nn)=0.0D0
        enddo   
       enddo
cccz cumulate everything that belongs to 'one new' element
       do kk=1,mulchLayer
        do nn=1,SurNodeIndex-1
         MulchEleThNew_temp(mergeIndex(kk),nn)=
     &    MulchEleThNew_temp(mergeIndex(kk),nn)
     &    +MulchEleThNew(kk,nn)*thickPerLayer(kk)
         MulchEleTmpr_temp(mergeIndex(kk),nn)=
     &    MulchEleTmpr_temp(mergeIndex(kk),nn)
     &    +MulchEleTmpr(kk,nn)*thickPerLayer(kk)
         RDM_mass_temp(mergeIndex(kk),nn)=
     &    RDM_mass_temp(mergeIndex(kk),nn)+RDM_mass(kk,nn)
         HCE_mass_temp(mergeIndex(kk),nn)=
     &    HCE_mass_temp(mergeIndex(kk),nn)+HCE_mass(kk,nn)
         CEL_mass_temp(mergeIndex(kk),nn)=
     &    CEL_mass_temp(mergeIndex(kk),nn)+CEL_mass(kk,nn)
         LIG_mass_temp(mergeIndex(kk),nn)=
     &    LIG_mass_temp(mergeIndex(kk),nn)+LIG_mass(kk,nn)
        enddo   
       enddo
cccz finish normailizaiton & assignment
       do kk=1,mulchLayer_temp
        thickPerLayer(kk)=thickPerLayer_temp(kk)   
        do nn=1,SurNodeIndex-1
         MulchEleThNew_temp(kk,nn)=
     &    MulchEleThNew_temp(kk,nn)/thickPerLayer_temp(kk)
         MulchEleTmpr_temp(kk,nn)=
     &    MulchEleTmpr_temp(kk,nn)/thickPerLayer_temp(kk)
         MulchEleThNew(kk,nn)=MulchEleThNew_temp(kk,nn)
         MulchEleTmpr(kk,nn)=MulchEleTmpr_temp(kk,nn)
         
         RDM_mass(kk,nn)=RDM_mass_temp(kk,nn)
         HCE_mass(kk,nn)=HCE_mass_temp(kk,nn)
         CEL_mass(kk,nn)=CEL_mass_temp(kk,nn)
         LIG_mass(kk,nn)=LIG_mass_temp(kk,nn)
         
         if(Van_Camp_mulch.eq.0) then 
          MulchElehNew(kk,nn)=
     &     WH_VANG_MULCH(MulchEleThNew(kk,nn)
     &     /(1-f_mulch_pore),ParMulch(:,1),lInput)
         elseif(Van_Camp_mulch.eq.1) then 
          MulchElehNew(kk,nn)=
     &     WH_CAMP_MULCH(MulchEleThNew(kk,nn)
     &     /(1-f_mulch_pore),ParMulch(:,1),lInput)
         endif  
        enddo   
       enddo
cccz reset zeros for disappeared layers
       do kk=mulchLayer_temp+1,mulchLayer
        thickPerLayer(kk)=0.0D0   
        do nn=1,SurNodeIndex-1
          MulchEleThNew(kk,nn)=0.0D0
          MulchElehNew(kk,nn)=0.0D0
          MulchEleTmpr(kk,nn)=0.0D0
          RDM_mass(kk,nn)=0.0D0
          HCE_mass(kk,nn)=0.0D0
          CEL_mass(kk,nn)=0.0D0
          LIG_mass(kk,nn)=0.0D0
        enddo
       enddo
      endif 
      
cccz Establish the new geometery of the mulch    
      numMulchNode_temp=0
      hhhh=0.0D0
      do i=1,mulchLayer_temp+1
       if(i.ge.2) then
        hhhh=hhhh+thickPerLayer(i-1)
       endif
       do n=1,SurNodeIndex
        numMulchNode_temp=numMulchNode_temp+1
        if(i.eq.1) then
cccz Here we copy the soil surface node (establish the interface)
         mulchNodeCoord(numMulchNode_temp,1)=slopeCoord(n,1)
         mulchNodeCoord(numMulchNode_temp,2)=slopeCoord(n,2)
         MulchNodeMarker(numMulchNode_temp,1)=n
         MulchNodeMarker(numMulchNode_temp,2)=1
        else
cccz Stack upper nodes above soil surface, within mulch
cccz we always assmue the equal-distance layer at beginning
         mulchNodeCoord(numMulchNode_temp,1)=slopeCoord(n,1)
         mulchNodeCoord(numMulchNode_temp,2)=slopeCoord(n,2)+hhhh
         MulchNodeMarker(numMulchNode_temp,1)=n
         MulchNodeMarker(numMulchNode_temp,2)=i
        endif
         MulchNodeMatrix(i,n)=numMulchNode_temp
         MulchNodeMarkerArea(numMulchNode_temp)=0.0D0
       enddo
      enddo
      
      do i=numMulchNode_temp+1,numMulchNode
         mulchNodeCoord(i,1)=0.0D0
         mulchNodeCoord(i,2)=0.0D0
         MulchNodeMarker(i,1)=0
         MulchNodeMarker(i,2)=0
         MulchNodeMarkerArea(i)=0.0D0
      enddo
      
      do i=mulchLayer_temp+2,mulchLayer+1
       do n=1,SurNodeIndex
         MulchNodeMatrix(i,n)=0
       enddo
      enddo
          
cccz Establish the mulching grid (Element)
      numMulchEle_temp=0
      do i=1,mulchLayer_temp
       do n=1,SurNodeIndex-1
        numMulchEle_temp=numMulchEle_temp+1
        MulchEleMarker(numMulchEle_temp,1)=MulchNodeMatrix(i,n)
        MulchEleMarker(numMulchEle_temp,2)=MulchNodeMatrix(i,n+1)
        MulchEleMarker(numMulchEle_temp,3)=MulchNodeMatrix(i+1,n+1)
        MulchEleMarker(numMulchEle_temp,4)=MulchNodeMatrix(i+1,n)
        MulchEleMarker(numMulchEle_temp,5)=n
        MulchEleMarker(numMulchEle_temp,6)=i
        MulchEleMatrix(i,n)=numMulchEle_temp
        MulchEleMarkerArea(numMulchEle_temp)=0.0D0
       enddo
      enddo
      do i=numMulchEle_temp+1,numMulchEle
        MulchEleMarker(i,1)=0
        MulchEleMarker(i,2)=0
        MulchEleMarker(i,3)=0
        MulchEleMarker(i,4)=0
        MulchEleMarker(i,5)=0
        MulchEleMarker(i,6)=0
        MulchEleMarkerArea(i)=0.0D0
      enddo
      do i=mulchLayer_temp+1,mulchLayer
       do n=1,SurNodeIndex-1
         MulchEleMatrix(i,n)=0
       enddo
      enddo

cccz Calculate the area and geometry of each node and element    
      do e=1,numMulchEle_temp
       n1=MulchEleMarker(e,1)
       n2=MulchEleMarker(e,2)
       n3=MulchEleMarker(e,3)
       n4=MulchEleMarker(e,4)
       x1=mulchNodeCoord(n1,1)
       x2=mulchNodeCoord(n2,1)
       x3=mulchNodeCoord(n3,1)
       x4=mulchNodeCoord(n4,1)
       y1=mulchNodeCoord(n1,2)
       y2=mulchNodeCoord(n2,2)
       y3=mulchNodeCoord(n3,2)
       y4=mulchNodeCoord(n4,2)                  
       if(x3.eq.x4.and.y3.eq.y4) then
cccz triangle shape element
cccz triangle element is hard to be implement here
cccz so far we decided not to implement triangular element
       else
cccz rectangular shape element
        CJ=x2-x1
        CK=x3-x1
        BJ=Y2-Y1
        BK=Y3-Y1            
        mulchVirtArea123=0.5D0*abs(CJ*BK-CK*BJ)                 
        CJ=x3-x1
        CK=x4-x1
        BJ=Y3-Y1
        BK=Y4-Y1
        mulchVirtArea123=mulchVirtArea123+0.5D0*abs(CJ*BK-CK*BJ)
        MulchEleMarkerArea(e)=mulchVirtArea123
        MulchNodeMarkerArea(n1)=MulchNodeMarkerArea(n1)
     &   +0.25D0*MulchEleMarkerArea(e)
        MulchNodeMarkerArea(n2)=MulchNodeMarkerArea(n2)
     &   +0.25D0*MulchEleMarkerArea(e)
        MulchNodeMarkerArea(n3)=MulchNodeMarkerArea(n3)
     &   +0.25D0*MulchEleMarkerArea(e)
        MulchNodeMarkerArea(n4)=MulchNodeMarkerArea(n4)
     &   +0.25D0*MulchEleMarkerArea(e)
        endif                         
      enddo        
      mulchLayer=mulchLayer_temp
      numMulchEle=numMulchEle_temp
      numMulchNode=numMulchNode_temp
cccz finish update the geometrical parameters for the grid
cccz ------------------------------------------------------------------------------
      mulchThick=0.0D0
      do kk=1,mulchLayer
          mulchThick=mulchThick+thickPerLayer(kk)
      enddo
cccz one criterion for non-mulch, exit directly
      if(mulchThick.lt.0.25D0) then
         return
      endif

cccz ------------------------------------------------------------------------------
cccz determine the submerge depth            
      LayerHeight(1)=0.0D0
      do n=1,mulchLayer
        LayerHeight(n+1)=LayerHeight(n)+thickPerLayer(n)
      enddo
cccz need to consider the current surface ponded water heights 'h_Pond'
      h_Pond_max=0.0D0
      do n=1,SurNodeIndex
        h_Pond_max=max(h_Pond_max,h_Pond(n))
      enddo

cccz determine if the ponded water reach which 'element-based layer'
cccz in another word, node of layer 'SubmergeIndex' under the water surface
      SubmergeIndex=0
      if(h_Pond_max.le.0.0D0) then
        SubmergeIndex=0
      else
       do n=1,mulchLayer+1
        if(h_Pond_max.ge.LayerHeight(n)) then
          SubmergeIndex=n
          if(n.le.mulchLayer) then
            DiffSubmerge=LayerHeight(n+1)-h_Pond_max
          endif
        endif
       enddo
      endif
cccz PerOccupation.gt.1.0 means the layer is considered as a full layer of water.
      if(SubmergeIndex.eq.mulchLayer+1) then
          PerOccupation=3.0D0
      else
        if(DiffSubmerge.le.thresholdThick) then
           PerOccupation=2.0D0
        else
cccz this value should be smaller than 1, but larger than 'thresholdThick/thickPerLayer(SubmergeIndex)'
           PerOccupation=DiffSubmerge/thickPerLayer(SubmergeIndex)
        endif
      endif
cccz finished determine how much mulch is under water
cccz ------------------------------------------------------------------------------

cccz ----- short wave radiaiton ---------------------------------------------------
cccz Calculate the shortwave downwards radiation in each interface among between adjacent mulch layers
cccz Based on the model in 'Simulating the radiation distribution within a barley-straw mulch' Novak 2000
cccz use random distribution Eq(3) for transmittivity
cccz We repeat it here because it is outside of the 'initialization codes', i.e., 'lInput.eq.1'
      
      if(SubmergeIndex.le.0) then
cccz in the case there is no ponded water
       do kk=mulchLayer,1,-1
        if(kk.eq.mulchLayer) then
         shortRadDir(kk+1)=RadSolarShort_Wea
        elseif(kk.eq.(mulchLayer-1)) then
         shortRadDir(kk+1)=RadSolarShort_Wea*(1.0D0-DeltaRshort)
        else
         shortRadDir(kk+1)=shortRadDir(kk+2)*(1.0D0-Omega*DeltaRshort)
        endif
       enddo
       shortRadDir(1)=shortRadDir(2)*(1.0D0-Omega*DeltaRshort)
      elseif(SubmergeIndex.le.mulchLayer) then
cccz elimiated the short direct raditaiton to the node under water 
       do kk=mulchLayer,SubmergeIndex,-1
        if(kk.eq.mulchLayer) then
         shortRadDir(kk+1)=RadSolarShort_Wea
        elseif(kk.eq.(mulchLayer-1)) then
         shortRadDir(kk+1)=RadSolarShort_Wea*(1.0D0-DeltaRshort)
        else
         shortRadDir(kk+1)=shortRadDir(kk+2)*(1.0D0-Omega*DeltaRshort)
        endif
       enddo
       shortRadDir(SubmergeIndex)=
     &   shortRadDir(SubmergeIndex+1)*(1.0D0-Omega*DeltaRshort)
       do kk=1,SubmergeIndex-1
        shortRadDir(kk)=0.0D0
       enddo
      else
cccz the whole mulch is under water 
       shortRadDir(mulchLayer+1)=RadSolarShort_Wea
       do kk=1,mulchLayer
        shortRadDir(kk)=0.0D0
       enddo
      endif

        
cccz First order short wave reflective (upwards)
      if(SubmergeIndex.le.0) then      
cccz this is for original none-ponded water case
       do i=1,SurNodeIndex-1
        mulch_np=SurfNodeNodeIndexH(i)
        mulch_nq=SurfNodeNodeIndexH(i+1)
        alpha_s=0.25D0-0.025D0*(ThNew(mulch_np)+ThNew(mulch_nq))     ! surface albedo
        shortRadFirst(1,i)=shortRadDir(1)*alpha_s                    ! reflection @ soil surface
        shortRadFirst(2,i)=shortRadDir(1)*alpha_s*(1.0D0-DeltaRshort)! the first mulching interface
     &      +alpha_m*(shortRadDir(2)-shortRadDir(1))
        do kk=2,mulchLayer                                            ! loop from second mulching interface
         shortRadFirst(kk+1,i)=shortRadDir(1)*alpha_s                 ! add reflections from all the layers beneth the LOI
     &      *(1.0D0-DeltaRshort)*((1.0D0-Omega*DeltaRshort)**(kk-1))
         do jj=1,kk
          if(jj.eq.kk) then
            shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)
     &        +alpha_m*(shortRadDir(jj+1)-shortRadDir(jj))
          elseif(jj.eq.(kk-1)) then
            shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)+alpha_m
     &        *(shortRadDir(jj+1)-shortRadDir(jj))*(1.0D0-DeltaRshort)
          else
            shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)+alpha_m
     &        *(shortRadDir(jj+1)-shortRadDir(jj))*(1.0D0-DeltaRshort)
     &        *((1.0D0-Omega*DeltaRshort)**(kk-jj-1))
          endif
         enddo
        enddo
       enddo
      elseif(SubmergeIndex.le.mulchLayer) then
cccz this is revised when ponded water occurred
       do i=1,SurNodeIndex-1
        mulch_np=SurfNodeNodeIndexH(i)
        mulch_nq=SurfNodeNodeIndexH(i+1)
        if(h_Pond(i).eq.0.0D0) then
          alpha_np=0.25D0-0.05D0*ThNew(mulch_np)     ! surface albedo
        else
          alpha_np=0.06D0*f_mulch_pore+0.20D0*(1.0D0-f_mulch_pore) ! surface albedo; weighted average of water and mulch albedo
        endif
        if(h_Pond(i+1).eq.0.0D0) then
          alpha_nq=0.25D0-0.05D0*ThNew(mulch_nq)     ! surface albedo
        else
          alpha_nq=0.06D0*f_mulch_pore+0.20D0*(1.0D0-f_mulch_pore) ! surface albedo; weighted average of water and mulch albedo
        endif
        alpha_s=0.5D0*(alpha_np+alpha_nq)            ! this also include the albedo for water
        do kk=1,SubmergeIndex-1
          shortRadFirst(kk,i)=0.0D0
        enddo          
        shortRadFirst(SubmergeIndex,i)=                              ! reflection liquid surface
     &    shortRadDir(SubmergeIndex)*alpha_s                        
        shortRadFirst(SubmergeIndex+1,i)=                            ! the first mulching interface
     &    shortRadDir(SubmergeIndex)*alpha_s*(1.0D0-DeltaRshort)
     &    +alpha_m*(shortRadDir(SubmergeIndex+1)
     &    -shortRadDir(SubmergeIndex))
        do kk=SubmergeIndex+1,mulchLayer                              ! loop from second mulching interface
          shortRadFirst(kk+1,i)=shortRadDir(SubmergeIndex)*alpha_s     ! add reflections from all the layers beneth the LOI
     &      *(1.0D0-DeltaRshort)*((1.0D0-Omega*DeltaRshort)
     &      **(kk-SubmergeIndex))
          do jj=SubmergeIndex,kk
            if(jj.eq.kk) then
              shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)
     &          +alpha_m*(shortRadDir(jj+1)-shortRadDir(jj))
            elseif(jj.eq.(kk-1)) then
              shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)+alpha_m
     &          *(shortRadDir(jj+1)-shortRadDir(jj))*(1.0D0-DeltaRshort)
            else
              shortRadFirst(kk+1,i)=shortRadFirst(kk+1,i)+alpha_m
     &          *(shortRadDir(jj+1)-shortRadDir(jj))*(1.0D0-DeltaRshort)
     &          *((1.0D0-Omega*DeltaRshort)**(kk-jj-1))
            endif
          enddo
        enddo
       enddo
      else
cccz the mulch is totally submerged
        alpha_s=0.06D0*f_mulch_pore+0.20D0*(1.0D0-f_mulch_pore)        ! combine the albedo for water and mulch, based on volume fraction
        do i=1,SurNodeIndex-1
            shortRadFirst(mulchLayer+1,i)= RadSolarShort_Wea*alpha_s   ! reflection liquid surface
        enddo
        do kk=1,mulchLayer
         do i=1,SurNodeIndex-1   
          shortRadFirst(kk,i)=0.0D0
         enddo
        enddo  
      endif

cccz finish the short wave radiation
cccz -----------------------------------------------------------------

 
cccz The wind calculation along soil surface is only proceed once because the ambient condition within one iteration is stable
cccz 'An introduction to environmental biophysics' Campbell & Norman
cccz first determine the wind speed near soil surface
      
      Phi_m=0.0D0                                                     ! the momentum correction factor, need more information later
      karman=0.4D0
      z_rough=0.079D0*mulchThick                                      ! just think mulch is a 'thick' canopy
      d_displace=0.87D0*mulchThick                                    ! just think mulch is a 'thick' canopy
      z_wind=200.0D0                                                  ! two meters above the surface is where the wind was measured
      u_ast=karman*AirWind_Wea/log((z_wind-d_displace)/z_rough)       ! calcualate the friction velocity (km/hour)
      u_above=u_ast/karman*log((mulchThick-d_displace)/z_rough)       ! wind speed in mulch-air surface (km/hour)
      u_soilsur=2.1D0*u_ast*log(2.2D0*thickPerLayer(1)/z_rough/2.0D0) ! wind speed in mulch-soil surface (km/hour, Novak et al. literature seq)       
      u_dist=0.0D0
cccz we do not make difference between ponded/unponded;
cccz if ponded, we just do not use the wind speed in deep layers
      do n=1,mulchLayer
       u_dist=u_dist+thickPerLayer(n)
       u_mulch(n)=2.1D0*u_ast*log(2.2D0*u_dist/z_rough)               ! wind speed in mulch-soil surface (km/hour, Novak et al. literature seq)
      enddo
      u_above=(u_above+u_mulch(mulchLayer))/2.0D0
      
cccz ----------------------------------------------------------
cccz every thing above this point is not changed within the mulch iteration scale
cccz we assume that the surface water (on flat soil) has same temperature as surface soil
      do n=1,SurNodeIndex-1
       do k=1,mulchLayer
        MulchElehNew_temp(k,n)=MulchElehNew(k,n)
        MulchEleTmpr_temp(k,n)=MulchEleTmpr(k,n)
        MulchElehNew_temp0(k,n)=MulchElehNew(k,n)
        MulchEleTmpr_temp0(k,n)=MulchEleTmpr(k,n)
        MulchElehNew_temp3(k,n)=0.0D0
        MulchEleTmpr_temp3(k,n)=0.0D0
        MulchEleThNew_temp(k,n)=MulchEleThNew(k,n)
       enddo
      enddo
      
      do n=1,SurNodeIndex
       VarBW1_temp(n)=0.0D0
       VarBW2_temp(n)=0.0D0
       VarBW3_temp(n)=0.0D0
       VarBT1_temp(n)=0.0D0
       VarBT2_temp(n)=0.0D0
       VarBT3_temp(n)=0.0D0
       VarBT4_temp(n)=0.0D0
       Q_temp(n)=0.0D0
      enddo
      
      IterMulch=1
      TotalTime=0.0D0
      LocalStep=Step
      TimeShrink=0
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=0.0D0
       enddo
      enddo
      
2001  aaaa=1
cccz rainfall redistribution initialization
cccz the submerged part will also be initialized but kept 0 (water contributed to VARBW) during the calculation
      do k=1,SurNodeIndex-1
       kSur=SurfNodeSurfIndexH(k)
cccz we must use 'Varbw_Air' now since VarBW counts the infiltraiton due to runoff/ponded water
cccz also because the location is @ mulch-air interface
       RainFallInput_temp(mulchLayer+1,k)=Varbw_Air(kSur,1)*10000.0D0  ! recast rainfall unit to g/m^2/day
       do n=1,mulchLayer
        RainFallInput_temp(n,k)=0.0D0
       enddo
      enddo

cccz --------------- Long wave radiation ------------------------------
cccz Long wave radiation downwards (W/m^2*1cm or 1m)
cccz 10.0D0 change vapor pressure from kPa to mb
      epsilon_a=epsilon_a_0+coef_epsilon_a*10.0D0*AirVaporP_Wea
     &  *exp(1500.0D0/(AirTemp_Wea+273.15D0))
      epsilon_a=(1.0D0-0.84D0*CloudCoverFactor_Wea)*epsilon_a
     &  +0.84D0*CloudCoverFactor_Wea
      if(SubmergeIndex.le.0.0D0) then 
cccz no submerged part
       longEmissionAir=epsilon_a*sigma*
     &   ((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex-1
        longRadDown(mulchLayer+1,i)=longEmissionAir
        do kk=mulchLayer,1,-1
         longRadDown(kk,i)=longEmissionAir*(1.0D0-DeltaRlong)
     &    *((1.0D0-Omega*DeltaRlong)**(mulchLayer-kk))
        enddo
        do kk=mulchLayer,1,-1
         longEmission=epsilon_m*sigma
     &    *((MulchEleTmpr_temp(kk,i)+273.15D0)**4.0D0)
         longRadDown(kk,i)=longRadDown(kk,i)
     &    +longEmission*(1.0D0-(1.0D0-DeltaRlong))
         do jj=kk-1,1,-1
          longRadDown(jj,i)=longRadDown(jj,i)
     &     +longEmission*(1.0D0-DeltaRlong)*
     &     ((1.0D0-Omega*DeltaRlong)**(kk-1-jj)
     &    -(1.0D0-Omega*DeltaRlong)**(kk-jj))
         enddo
        enddo
       enddo
      elseif(SubmergeIndex.le.mulchLayer) then
cccz this is revised when ponded water occurred
       longEmissionAir=epsilon_a*sigma*
     &     ((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex-1
        longRadDown(mulchLayer+1,i)=longEmissionAir
        do kk=mulchLayer,SubmergeIndex,-1
         longRadDown(kk,i)=longEmissionAir*(1.0D0-DeltaRlong)
     &    *((1.0D0-Omega*DeltaRlong)**(mulchLayer-kk))
        enddo
        do kk=mulchLayer,SubmergeIndex,-1
         longEmission=epsilon_m*sigma
     &    *((MulchEleTmpr_temp(kk,i)+273.15D0)**4.0D0)
         longRadDown(kk,i)=longRadDown(kk,i)
     &    +longEmission*(1.0D0-(1.0D0-DeltaRlong))
         do jj=kk-1,SubmergeIndex,-1
          longRadDown(jj,i)=longRadDown(jj,i)
     &     +longEmission*(1.0D0-DeltaRlong)*
     &     ((1.0D0-Omega*DeltaRlong)**(kk-1-jj)
     &    -(1.0D0-Omega*DeltaRlong)**(kk-jj))
         enddo
        enddo
        do kk=1,SubmergeIndex-1
            longRadDown(kk,i)=0.0D0    ! we do not consider the long-rad emission 
        enddo
       enddo
      else
       longEmissionAir=epsilon_a*sigma*
     &   ((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex-1
        longRadDown(mulchLayer+1,i)=longEmissionAir
        do kk=1,mulchLayer
          longRadDown(kk,i)=0.0D0
        enddo
       enddo
      endif

cccz Long wave radiation upwards (W/m^2*1cm or 1m)
      if(SubmergeIndex.le.0.0D0) then 
cccz assume no submerged part
       do i=1,SurNodeIndex-1
        mulch_np=SurfNodeNodeIndexH(i)
        mulch_nq=SurfNodeNodeIndexH(i+1)
        longEmissionSoil=epsilon_s*sigma*
     &   ((0.5D0*(Tmpr(mulch_np)+Tmpr(mulch_nq))+273.15D0)**4.0D0)
        longRadUp(1,i)=longEmissionSoil
        do kk=2,mulchLayer+1
         longRadUp(kk,i)=longEmissionSoil*(1.0D0-DeltaRlong)
     &    *((1.0D0-Omega*DeltaRlong)**(kk-2))
        enddo
        do kk=2,mulchLayer+1
         longEmission=epsilon_m*sigma
     &    *((MulchEleTmpr_temp(kk-1,i)+273.15D0)**4.0D0)
         longRadUp(kk,i)=longRadUp(kk,i)
     &    +longEmission*(1.0D0-(1.0D0-DeltaRlong))
         do jj=kk+1,mulchLayer+1
          longRadUp(jj,i)=longRadUp(jj,i)
     &     +longEmission*(1.0D0-DeltaRlong)*
     &     ((1.0D0-Omega*DeltaRlong)**(jj-kk-1)
     &     -(1.0D0-Omega*DeltaRlong)**(jj-kk))
         enddo
        enddo
       enddo
      elseif(SubmergeIndex.le.mulchLayer) then
cccz this is revised when ponded water occurred
       do i=1,SurNodeIndex-1
        do kk=1,SubmergeIndex-1
            longRadUp(kk,i)=0.0D0    ! we do not consider the long-rad emission 
        enddo   
        mulch_np=SurfNodeNodeIndexH(i)
        mulch_nq=SurfNodeNodeIndexH(i+1)
cccz here we assume that the surface ponded water has the same temp as surface soil,
cccz ??? maybe not true, but let us fixed it later.
cccz also use 'epsilon_s=1.0D0' for water cause water is considered as a blackbody
        longEmissionSoil=epsilon_s*sigma*
     &   ((0.5D0*(Tmpr(mulch_np)+Tmpr(mulch_nq))+273.15D0)**4.0D0)
        longRadUp(SubmergeIndex,i)=longEmissionSoil
        do kk=SubmergeIndex+1,mulchLayer+1
         longRadUp(kk,i)=longEmissionSoil*(1.0D0-DeltaRlong)
     &    *((1.0D0-Omega*DeltaRlong)**(kk-(SubmergeIndex+1)))
        enddo
        do kk=SubmergeIndex+1,mulchLayer+1
         longEmission=epsilon_m*sigma
     &    *((MulchEleTmpr_temp(kk-1,i)+273.15D0)**4.0D0)
         longRadUp(kk,i)=longRadUp(kk,i)
     &    +longEmission*(1.0D0-(1.0D0-DeltaRlong))
         do jj=kk+1,mulchLayer+1
          longRadUp(jj,i)=longRadUp(jj,i)
     &     +longEmission*(1.0D0-DeltaRlong)*
     &     ((1.0D0-Omega*DeltaRlong)**(jj-kk-1)
     &     -(1.0D0-Omega*DeltaRlong)**(jj-kk))
         enddo
        enddo
       enddo
      else
       do i=1,SurNodeIndex-1
        mulch_np=SurfNodeNodeIndexH(i)
        mulch_nq=SurfNodeNodeIndexH(i+1)
cccz here we assume that the surface ponded water has the same temp as surface soil,
cccz ??? maybe not true, but let us fixed it later.
cccz also use 'epsilon_s=1.0D0' for water cause water is considered as a blackbody
        longEmissionSoil=epsilon_s*sigma*
     &    ((0.5D0*(Tmpr(mulch_np)+Tmpr(mulch_nq))+273.15D0)**4.0D0)
        longRadUp(mulchLayer+1,i)=longEmissionSoil
        do kk=1,mulchLayer
          longRadUp(kk,i)=0.0D0
        enddo
       enddo
      endif
cccz finish the calculation for long wave radiation
cccz -------------------------------------------------------------------------------
      
cccz --------- combine long/short wave radiation ----------------------------------
cccz Summarize the net radiation at each level
cccz (pay attention to the radiation unit, not based on per-day)
      do i=1,SurNodeIndex-1
       do j=1,mulchLayer+1
        netRad(j,i)=(shortRadDir(j)-shortRadFirst(j,i))
     &    +(longRadDown(j,i)-longRadUp(j,i))
       enddo
      enddo
cccz Summarize the net radiation received by each mulching layer (W/m^2)
      do i=1,SurNodeIndex-1
       do j=1,mulchLayer
        netRadEachLayer(j,i)=netRad(j+1,i)-netRad(j,i)
       enddo
      enddo
cccz finish combining the radiation
cccz now we have the radiation energy distribution for each layer
cccz -----------------------------------------------------------------------------------------
      
cccz ***************************
c now process the computaiton of fluxes and conversation laws
      
      if(SubmergeIndex.le.0.0D0) then 
cccz this is the unsubmerged case
          goto 2009
      elseif(SubmergeIndex.gt.0.0D0.and.PerOccupation.eq.3.0) then
cccz this is the totally submerged case
          goto 2400
      else
cccz this is the partially submerged case
          goto 2300
      endif

cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz NO PONDED WATER
2009  do n=1,SurNodeIndex-1
       do k=1,mulchLayer 
cccz calculate the 'vapor conditions'
        RelaHumid_mulch_temp(k,n)=exp(2.124D-4*MulchElehNew_temp(k,n) ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))
         RelaHumid_mulch_temp(k,n)=min(RelaHumid_mulch_temp(k,n),1.0D0)! make sure the relative humidity is of currect range      
         RelaHumid_mulch_temp(k,n)=max(RelaHumid_mulch_temp(k,n),0.0D0)
        VaporSat_mulch_temp(k,n)=exp(19.84D0-4975.9D0     ! Saturated Vapor density for each mulch element 
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))              ! (g/m^3)           
        VaporAct_mulch_temp(k,n)=VaporSat_mulch_temp(k,n) ! Actual Vapor density for each mulch element
     &   *RelaHumid_mulch_temp(k,n)                       ! (g/m^3)
cccz water difficit in each element
cccz including liquid water and water vapor to a certain degree of saturation
        WaterDifici_mulch_temp(k,n)=
     &    VaporSat_mulch_temp(k,n)-VaporAct_mulch_temp(k,n)
        WaterDifici_mulch_temp(k,n)=WaterDifici_mulch_temp(k,n)
     &    *(thickPerLayer(k)/100.0D0)/LocalStep           ! rate to saturate the vapor portion (g/m^2/day)
cccz water diffusivity calculation
cccz use diffusion type water fluxes
        VaporDiff_mulch_temp(k,n)=2.29D-5*
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**0.667D0)*f_mulch_pore*86400.0D0              ! Water Vapor diffusivity (m^2/day)
cccz water capacity (only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore)
cccz use soil water capacity, but recommend to use wood water retention curve in future.        
        if(Van_Camp_mulch.eq.0) then  
         WaterCapaMulch_temp(k,n)=
     &     WC_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        elseif(Van_Camp_mulch.eq.1) then
         WaterCapaMulch_temp(k,n)=
     &     WC_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        endif
cccz heat capacity and diffusivity calculation
cccz air density (rho_dryair) g m^-3, other parameters 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air_temp(k,n)=1.85D0*VaporAct_mulch_temp(k,n)
     &   +1.006D0*(rho_dryair-VaporAct_mulch_temp(k,n))                 ! heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=HeatCapa_Air_temp(k,n)*1.84896D0       ! heat diffusivity (J/m/K/day), 1.84896D0=2.14D-5(heat diff)*86400.0D0      
       enddo
      enddo

cccz calculate the water fluxes within the mulch grid
cccz for water, we assume only vapor flux occur, but allow some liquid water stored in the mulch solid materials.
cccz downward and leftward are positive flux (water income)
      do k=1,mulchLayer+1    
cccz for mulch-soil interface
       if(k.eq.1) then                                                
        do n=1,SurNodeIndex-1
cccz the soil surface vapor flux based on the 'potential evaporation' method
cccz maybe need to adjust the inwards/outwards flux to soil
         SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(1,n))
     &    /(MulchEleTmpr_temp(1,n)+237.3D0))
         DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(1,n)+1.0D0))
     &    /(MulchEleTmpr_temp(1,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch_temp(1,n))          ! calculate VPD for the first mulching layer
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, 
                                                                     ! should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(1,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur*max(netRad(1,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*MulchEleTmpr_temp(1,n))))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_soilsur))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 

cccz the soil surface vapor flux based on the 'vapor flux' method
cccz first calcualte the dimensionless parameters
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         hNew_Sur=0.5D0*(hNew(qLeft)+hNew(qRight))
         VaporSat_Sur=exp(19.84D0-4975.9D0/(Tmpr_Sur+273.15D0))
         RelaHumid_Sur=exp(2.124D-4*hNew_Sur/(Tmpr_Sur+273.15D0))          ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
          RelaHumid_Sur=min(RelaHumid_Sur,1.0D0)
          RelaHumid_Sur=max(RelaHumid_Sur,0.0D0)
         VaporAct_Sur=VaporSat_Sur*RelaHumid_Sur
         
         Ra_mulch(1,n)=9.81D0*abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur)
     &    *(0.5D0*thickPerLayer(1)/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur)/2.0D0)/Dhm
         GrRe_mulch(1,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur)
     &    /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur)/2.0D0)
     &    /((u_soilsur/3.6D0)**2.0D0)
     
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul=VaporDiff_mulch_temp(1,n)
     &      /(0.5D0*thickPerLayer(1)/100.0D0)                         ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
            VaporCond_Mul=0.055D0*((abs(Tmpr_Sur
     &      -MulchEleTmpr_temp(1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_soilsur/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_vapor_try2=VaporCond_Mul
     &    *(VaporAct_mulch_temp(1,n)-VaporAct_Sur)
cccz take the min between 'diffusive type' and 'potential' evporation     
cccz if the two values are of different directions, we assume zero flux fo numerical stable
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then
           g_vapor(1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(1,n)=0.0D0
         endif
        enddo      
        g_Vapor(2,1)=0.0D0                                            ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Vapor(2,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &     /(widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(1,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch_temp(1,n))
     &     *(VaporAct_mulch_temp(1,n)-VaporAct_mulch_temp(1,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
        enddo  
        
cccz for mulch-air interface        
       elseif(k.eq.(mulchLayer+1)) then                               
        do n=1,SurNodeIndex-1                                        ! use conduction/convection to calculate conductivity at atmosphere surface
         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    *((z_wind-mulchThick+0.5D0*thickPerLayer(k-1))/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)
     &    /((u_above/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul=sqrt(2.40D-5/((z_wind-mulchThick
     &      +0.5D0*thickPerLayer(k-1))/100.0D0)*86400.0D0
     &      *VaporDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0)) ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(AirTemp_Wea
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_above/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_Ambient-VaporAct_mulch_temp(k-1,n))
        enddo  

cccz for water flux within mulch
       else                                                           
        do n=1,SurNodeIndex-1
         Ra_mulch(k,n)=9.81D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    *((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
     &    /nu_air/((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0
     &    *sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))**2.0D0
     &    +(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/((u_mulch(k)/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch_temp(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch_temp(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(MulchEleTmpr_temp(k,n)
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_mulch(k)/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch_temp(k,n))
     &     *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
       endif
      enddo
           
cccz calculate the rainfall redistribution
      do n=1,SurNodeIndex-1
       k=mulchLayer+1
2008   if(k.gt.1) then
        inputPerLayer=
     &   min(RainFallInput_temp(k,n),WaterDifici_mulch_temp(k-1,n))
        if(inputPerLayer.lt.RainFallInput_temp(k,n)) then
         RainFallInput_temp(k-1,n)=RainFallInput_temp(k,n)-inputPerLayer
         RainFallInput_temp(k,n)=inputPerLayer
         k=k-1
         goto 2008
        else
         continue
        endif
       else
        continue
       endif
      enddo
     
cccz calculate the latent heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes based on the upwind direction of water
      do k=1,mulchLayer+1
cccz for soil surface layer
       if(k.eq.1) then    
        do n=1,SurNodeIndex-1
         if(g_Vapor(1,n).ge.0.0D0) then                               ! downwards vapor flow
           g_Heat_Latent(1,n)=c_vap_s*g_Vapor(1,n)
     &      *(MulchEleTmpr_temp(1,n)-Tmpr_ref)+c_vap_v*g_Vapor(1,n)
         else                                                         ! upwards vapor flow
           qLeft=SurfNodeNodeIndexH(n)
           qRight=SurfNodeNodeIndexH(n+1)
           Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
           g_Heat_Latent(1,n)=c_vap_s*g_Vapor(1,n)*(Tmpr_Sur-Tmpr_ref)
     &      +c_vap_v*g_Vapor(1,n)
         endif
        enddo
        g_Heat_Latent(2,1)=0.0D0
        g_Heat_Latent(2,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2,n).ge.0.0D0) then                               ! leftwards vapor flow
           g_Heat_Latent(2,n)=c_vap_s*g_Vapor(2,n)
     &      *(MulchEleTmpr_temp(1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2,n)
         else                                                         ! rightwards vapor flow
           g_Heat_Latent(2,n)=c_vap_s*g_Vapor(2,n)
     &      *(MulchEleTmpr_temp(1,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2,n)
         endif
        enddo
cccz for mulch-air interface layer
       elseif(k.eq.(mulchLayer+1)) then
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
cccz for other interier layers
       else
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                              ! leftwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n) 
         else                                                          ! rightwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
       endif   
      enddo
      
cccz calculate the sensible (diff) heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes along soil surface using the current soil surface temperature   
      do k=1,mulchLayer+1  
cccz for mulch-soil interface
       if(k.eq.1) then                                                  
        do n=1,SurNodeIndex-1
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul=HeatDiff_mulch_temp(1,n)
     &     /(thickPerLayer(1)/200.0D0)                                  ! the heat conductance (J/m/day/K)
         else
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul=0.050D0*((abs(Tmpr_Sur
     &     -MulchEleTmpr_temp(1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(1,n)
     &     *86400.0D0/41.4D0
          else
           HeatCond_Mul=0.135D0*sqrt((u_soilsur/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *HeatCapa_Air_temp(1,n)*86400.0D0/41.4D0
          endif
         endif
        g_Heat_Sensi(1,n)=HeatCond_Mul
     &    *(MulchEleTmpr_temp(1,n)-Tmpr_Sur)! the sensible heat conductance (J/m^2/day) 
        enddo
        g_Heat_Sensi(2,1)=0.0D0
        g_Heat_Sensi(2,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Heat_Sensi(2,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(1,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch_temp(1,n))
     &    *(MulchEleTmpr_temp(1,n)-MulchEleTmpr_temp(1,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
cccz for mulch-air interface 
       elseif(k.eq.(mulchLayer+1)) then                                 
        do n=1,SurNodeIndex-1
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul=sqrt(2.14D-5/((z_wind-mulchThick+0.05D0
     &      *thickPerLayer(k-1))/100.0D0)*HeatCapa_Air_temp(k-1,n)
     &      *86400.0D0                                                   ! the heat conductance (J/m/day/K)
     &      *HeatDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0))
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           HeatCond_Mul=0.050D0*((abs(AirTemp_Wea
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(k-1,n)
     &     *86400.0D0/41.4D0
          else
           HeatCond_Mul=0.135D0*sqrt((u_above/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *HeatCapa_Air_temp(k-1,n)*86400.0D0/41.4D0
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=HeatCond_Mul
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
cccz for mulch interior layers 
      else                                    
       do n=1,SurNodeIndex-1
        if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch_temp(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch_temp(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)      ! the heat conductance (J/m/day/K)
        else
         if(GrRe_mulch(k,n).ge.GrRe_Critical) then
          HeatCond_Mul=0.050D0*((abs(MulchEleTmpr_temp(k,n)
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*(HeatCapa_Air_temp(k-1,n)
     &     +HeatCapa_Air_temp(k,n))/2.0D0*86400.0D0/41.4D0
         else
          HeatCond_Mul=0.135D0*sqrt((u_mulch(k)/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *(HeatCapa_Air_temp(k-1,n)+HeatCapa_Air_temp(k,n))
     &     /2.0D0*86400.0D0/41.4D0
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=HeatCond_Mul*(MulchEleTmpr_temp(k,n)
     &   -MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex)=0.0D0
       do n=2,SurNodeIndex-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch_temp(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
       enddo 
       endif
      enddo
      
cccz finish the calculation of vapor and heat fluxes
cccz start to establish the equation systems which solve the temperature and water potential
      do n=1,SurNodeIndex-1
       do k=1,mulchLayer
        a11=WaterCapaMulch_temp(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0)              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    -WaterCapaMulch_temp(k,n)*VaporAct_mulch_temp(k,n))
     &    +c_water_s*rho_w*WaterCapaMulch_temp(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_air_s*rho_dryair+c_water_s*rho_w*MulchEleThNew_temp(k,n)
     &    +c_vap_s*VaporAct_mulch_temp(k,n)+(c_vap_v+c_vap_s
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref))*f_mulch_pore
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &   +RainFallInput_temp(k+1,n)/(thickPerLayer(k)/100.0D0)*LocalStep
        b12=(netRadEachLayer(k,n)*86400.0D0              
     &    +(g_Heat_Sensi(2*k+1,n)-g_Heat_Sensi(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Sensi(2*k,n+1)-g_Heat_Sensi(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0)
     &    +(g_Heat_Latent(2*k+1,n)-g_Heat_Latent(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Latent(2*k,n+1)-g_Heat_Latent(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +c_water_s*RainFallInput_temp(k+1,n)*(AirTemp_Wea-Tmpr_ref)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep

        SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)                    ! Cramer's rule for linear system
        SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
       enddo
      enddo
      
cccz finish the tempted calculation, start to check the convergence
      maxhNew=0.0D0
      maxTmpr=0.0D0
      maxhNewDiff=0.0D0
      maxTmprDiff=0.0D0
      if(IterMulch.le.mulchLayer) then                            ! the Picard iteration has to be processed for times >= mulch layer
       IterMulch=IterMulch+1
       do n=1,SurNodeIndex-1
        do k=1,mulchLayer
cccz protection
        if(MulchElehNew_temp2(k,n).gt.0.0D0) goto 2007
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) goto 2007  ! NAN case occurs
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) goto 2007  ! NAN case occurs
cccz protection
         MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
         MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
         if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         endif
        enddo
       enddo
       goto 2001
      else
       IF(IterMulch.le.MaxIter) then
        do n=1,SurNodeIndex-1
        do k=1,mulchLayer
cccz protection
        if(MulchElehNew_temp2(k,n).gt.0.0D0) goto 2007
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) goto 2007  ! NAN case occurs
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) goto 2007  ! NAN case occurs
cccz protection
        MulchElehNew_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchElehNew_temp3(k,n)
     &     +MulchElehNew_temp2(k,n))/dble(IterMulch-mulchLayer)
        MulchEleTmpr_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchEleTmpr_temp3(k,n)
     &     +MulchEleTmpr_temp2(k,n))/dble(IterMulch-mulchLayer)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp3(k,n)
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp3(k,n)
        maxhNewDiff=max(maxhNewDiff, 
     &    abs(MulchElehNew_temp(k,n)-MulchElehNew_temp2(k,n)))
        maxTmprDiff=max(maxTmprDiff,
     &    abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp2(k,n)))
        maxhNew=max(maxhNew,abs(MulchElehNew(k,n)))
        maxTmpr=max(maxTmpr,abs(MulchEleTmpr(k,n)))  
        MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
        MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
        if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
        elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
        endif
        enddo
        enddo
        IterMulch=IterMulch+1
        if(max(maxhNewDiff/maxhNew,maxTmprDiff/maxTmpr).gt.errTol) then
         goto 2001
        else
         TotalTime=TotalTime+LocalStep
         if(TotalTime-Step.gt.-0.01D0*dtMin) then
           if(TimeShrink.eq.0) then   ! if no 'local time step' was used, 'one step' was used for calculating water/energy exchange on soil surface
            goto 2002
           else                       ! if multiple 'local time step' was used, then 'cumulative averaged' water/energy exchange were calculated in 2003 and 2004
            goto 2003
           endif
         else
          do n=1,SurNodeIndex-1
           do k=1,mulchLayer
            MulchElehNew_temp0(k,n)=MulchElehNew_temp(k,n)
            MulchEleTmpr_temp0(k,n)=MulchEleTmpr_temp(k,n)
           enddo
          enddo
          goto 2003
2004      LocalStep=min(LocalStep*DMul1,Step-TotalTime)
cccz??? ask a question, why not reset the mulch iteration index
cccz let us reset it for now
          IterMulch=1
c          IterMulch=1
          goto 2001
         endif
        endif
       ELSE
2007     TimeShrink=1
         do n=1,SurNodeIndex-1
          do k=1,mulchLayer
          MulchElehNew_temp(k,n)=MulchElehNew_temp0(k,n)
          MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp0(k,n)
          enddo
         enddo
         LocalStep=max(LocalStep*DMul2,0.0001D0*dtMin)
         IterMulch=1
         goto 2001
       ENDIF
      endif
    
cccz finish the Picard iteration, but we have to update the water and energy fluxes across soil surface
cccz here we solve for the potential vapor flux at the soil surface
cccz if the actural one match the potential one, we do not need to solve the equation system later in the post process codes
cccz if not use the actural one and re-solve the equation system.
cccz first do the vapor flux calculation
2002  do n=1,SurNodeIndex-1
cccz the unit within this module will be g/day/m^2
cccz the unit taken in watermov module is g/day/cm^2
       if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,n)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(1,n)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3) 
       elseif(n.eq.(SurNodeIndex-1)) then
         kSurL=SurfNodeSurfIndexH(n)
         kSurR=SurfNodeSurfIndexH(n+1)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(1,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)
         
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,n)/10000.0D0
         Varbw_Mulch(kSurR,2)=-g_Vapor(1,n)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)
       else
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(1,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)  
       endif 
      enddo   
      
cccz second update the heat fluxes
      do n=1,SurNodeIndex-1
cccz the unit within this module will be J/day/m^2
cccz ********* ideas for the parameter assignment **************
c    VarBT(,1) is the temperature, as assign the temperature at mulch bottom surface (not important, I need fluxes more)
c    VarBT(,2) and VarBT(,3) are coupled, VarBT(,2) looks like a coefficient and VarBT(,3)=VarBT(,2)*temperature
c         The only temperature dependent thing is sensible heat in this code, so
c         VarBT(,2) sensible heat conductance,  VarBT(,3)= VarBT(,2)*(temperature at mulch bottom surface)
c    VarBT(,4) is the net radiation, we have short wave and long wave radiation here.
c    Latent heat is not presented in VarBT, it will depend on VarBW and soil temperature.
c    Make unit change to Cal/cm^2/min as shown in hourwea.for
cccz ********* end ****************
       if(n.eq.1) then
        kSurL=SurfNodeSurfIndexH(n)
        VarBT(kSurL,1)=MulchEleTmpr_temp(1,n)
        VarBT(kSurL,2)=HeatDiff_mulch_temp(1,n)*0.025104D0           ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*MulchEleTmpr_temp(1,n)
        VarBT(kSurL,4)=netRad(1,n)*8.64D0                           ! 8.64D0=86400D0/10000D0
       elseif(n.eq.(SurNodeIndex-1)) then
        kSurL=SurfNodeSurfIndexH(n)
        kSurR=SurfNodeSurfIndexH(n+1)
        VarBT(kSurL,1)=(MulchEleTmpr_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
       VarBT(kSurL,2)=(HeatDiff_mulch_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
        VarBT(kSurL,4)=(netRad(1,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
        VarBT(kSurR,1)=MulchEleTmpr_temp(1,n)
        VarBT(kSurR,2)=HeatDiff_mulch_temp(1,n)*0.025104D0            ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(1,n)
        VarBT(kSurR,4)=netRad(1,n)*8.64D0                            ! 8.64D0=86400D0/10000D0
       else
        kSurL=SurfNodeSurfIndexH(n)
        VarBT(kSurL,1)=(MulchEleTmpr_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
       VarBT(kSurL,2)=(HeatDiff_mulch_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
        VarBT(kSurL,4)=(netRad(1,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
       endif
      enddo
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput_temp(n,k)
       enddo
      enddo
      goto 2006
      
cccz need to estimate the soil-mulch interface water and energy flux for each time segment during the iteration
cccz use the same scheme showed following Index 2002
2003  do n=1,SurNodeIndex-1
       if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
         Local_VarBW2=-g_Vapor(1,n)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
       elseif(n.eq.(SurNodeIndex-1)) then
         kSurL=SurfNodeSurfIndexH(n)
         kSurR=SurfNodeSurfIndexH(n+1)
         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW2=-(g_Vapor(1,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
         Local_VarBW2=-g_Vapor(1,n)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n+1)=VarBW1_temp(n+1)+Local_VarBW1*LocalStep
         VarBW2_temp(n+1)=VarBW2_temp(n+1)+Local_VarBW2*LocalStep
         VarBW3_temp(n+1)=VarBW3_temp(n+1)+Local_VarBW3*LocalStep
         Q_temp(n+1)=Q_temp(n+1)+(-Width(kSurR)*Local_VarBW3)*LocalStep
       else
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW2=-(g_Vapor(1,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW3=Local_VarBW2-VarBW(kSurL,1)
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
       endif 
      enddo
      do n=1,SurNodeIndex-1
      if(n.eq.1) then
         Local_VarBT1=MulchEleTmpr_temp(1,n)
         Local_VarBT2=HeatDiff_mulch_temp(1,n)*0.025104D0           ! 0.025104D0=4.184D0*60.0D0/10000.0D0
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=netRad(1,n)*8.64D0                           ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
       elseif(n.eq.(SurNodeIndex-1)) then
        Local_VarBT1=(MulchEleTmpr_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        Local_VarBT2=(HeatDiff_mulch_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=(netRad(1,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
        VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
        VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
        VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        Local_VarBT1=MulchEleTmpr_temp(1,n)
        Local_VarBT2=HeatDiff_mulch_temp(1,n)*0.025104D0            ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=netRad(1,n)*8.64D0                            ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n+1)=VarBT1_temp(n+1)+Local_VarBT1*LocalStep
        VarBT2_temp(n+1)=VarBT2_temp(n+1)+Local_VarBT2*LocalStep
        VarBT3_temp(n+1)=VarBT3_temp(n+1)+Local_VarBT3*LocalStep
        VarBT4_temp(n+1)=VarBT4_temp(n+1)+Local_VarBT4*LocalStep
       else
        Local_VarBT1=(MulchEleTmpr_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        Local_VarBT2=(HeatDiff_mulch_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=(netRad(1,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
        VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
        VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
        VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
       endif
      enddo
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)
     &   +RainFallInput_temp(n,k)*LocalStep
       enddo
      enddo
      if(TotalTime-Step.gt.-0.01D0*dtMin) then
          goto 2005   ! finalize the iteration
      else
          goto 2004   ! process to the next local step
      endif

cccz use the 'cumulative averaged' water/energy exchange between 2003 and 2004
cccz determine the vapor and energy exchange on soil surface
2005  do n=1,SurNodeIndex
cccz water vapor part
       kSur=SurfNodeSurfIndexH(n)
       Varbw_Mulch(kSur,1)=VarBW1_temp(n)/TotalTime
       Varbw_Mulch(kSur,2)=VarBW2_temp(n)/TotalTime
       Varbw_Mulch(kSur,3)=VarBW3_temp(n)/TotalTime
       VarBW(kSur,1)=Varbw_Mulch(kSur,1)
       VarBW(kSur,2)=Varbw_Mulch(kSur,2)
       VarBW(kSur,3)=Varbw_Mulch(kSur,3)
       nNode=KXB(kSur)
       Q(nNode)=Q_temp(n)/TotalTime
cccz energy part
       VarBT(kSur,1)=VarBT1_temp(n)/TotalTime
       VarBT(kSur,2)=VarBT2_temp(n)/TotalTime
       VarBT(kSur,3)=VarBT3_temp(n)/TotalTime
       VarBT(kSur,4)=VarBT4_temp(n)/TotalTime
      enddo
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)/TotalTime
       enddo
      enddo
      goto 2006


2006  aaaa=2
cccz final assignment
      do n=1,SurNodeIndex-1
       do k=1,mulchLayer
         MulchElehNew(k,n)=MulchElehNew_temp(k,n)
         MulchEleTmpr(k,n)=MulchEleTmpr_temp(k,n)
         MulchEleThNew(k,n)=MulchEleThNew_temp(k,n)
       enddo
      enddo
 
      goto 2500
      
cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz PONDED WATER CASE, but the ponded depth does not exceed the mulch height
2300  aaaa=3
      do n=1,SurNodeIndex-1
cccz      ! under water part
       do k=1,SubmergeIndex-1
        RelaHumid_mulch_temp(k,n)=1.0D0
        VaporSat_mulch_temp(k,n)=1.0D6     ! density of water
        VaporAct_mulch_temp(k,n)=1.0D6     ! density of water
        WaterDifici_mulch_temp(k,n)=0.0D0
        VaporDiff_mulch_temp(k,n)=0.0D0
        WaterCapaMulch_temp(k,n)=0.0D0
        HeatCapa_Air_temp(k,n)=4.182D6     ! heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=0.0D0
       enddo
cccz      ! right @ water level 
       k=SubmergeIndex
       if(PerOccupation.eq.2.0D0) then
       ! the water level is high enough, this layer is now a water layer
        RelaHumid_mulch_temp(k,n)=1.0D0
        VaporSat_mulch_temp(k,n)=1.0D6     ! density of water
        VaporAct_mulch_temp(k,n)=1.0D6     ! density of water
        WaterDifici_mulch_temp(k,n)=0.0D0
        VaporDiff_mulch_temp(k,n)=0.0D0
        WaterCapaMulch_temp(k,n)=0.0D0
        HeatCapa_Air_temp(k,n)=4.182D6     ! heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=0.0D0
       else
       ! the water level is not high enough, calculate vapor for the rest of this layer 
        RelaHumid_mulch_temp(k,n)=exp(2.124D-4*MulchElehNew_temp(k,n) ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))
         RelaHumid_mulch_temp(k,n)=min(RelaHumid_mulch_temp(k,n),1.0D0)! make sure the relative humidity is of currect range      
         RelaHumid_mulch_temp(k,n)=max(RelaHumid_mulch_temp(k,n),0.0D0)
        VaporSat_mulch_temp(k,n)=exp(19.84D0-4975.9D0     ! Saturated Vapor density for each mulch element 
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))              ! (g/m^3)           
        VaporAct_mulch_temp(k,n)=VaporSat_mulch_temp(k,n) ! Actual Vapor density for each mulch element
     &   *RelaHumid_mulch_temp(k,n)                       ! (g/m^3)
          ! water difficit in each element, including liquid water and water vapor to a certain degree of saturation
        WaterDifici_mulch_temp(k,n)=
     &    VaporSat_mulch_temp(k,n)-VaporAct_mulch_temp(k,n)
        thickLayer_Sur=thickPerLayer(k)*PerOccupation    ! need to adjust the height
        WaterDifici_mulch_temp(k,n)=WaterDifici_mulch_temp(k,n)
     &     *(thickLayer_Sur/100.0D0)/LocalStep           ! rate to saturate the vapor portion (g/m^2/day)
          ! water diffusivity calculation, use diffusion type water fluxes
        VaporDiff_mulch_temp(k,n)=2.29D-5*
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**0.667D0)*f_mulch_pore*86400.0D0              ! Water Vapor diffusivity (m^2/day), (f_mulch_pore**0.667D0) is the turtorosity
          ! water capacity (only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore)
          ! use soil water capacity, but recommend to use wood water retention curve in future.        
        if(Van_Camp_mulch.eq.0) then  
         WaterCapaMulch_temp(k,n)=
     &     WC_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        elseif(Van_Camp_mulch.eq.1) then
         WaterCapaMulch_temp(k,n)=
     &     WC_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        endif
          ! heat capacity and diffusivity calculation
          ! air density (rho_dryair) g m^-3, other parameters 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air_temp(k,n)=1.85D0*VaporAct_mulch_temp(k,n)
     &   +1.006D0*(rho_dryair-VaporAct_mulch_temp(k,n))                ! heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=HeatCapa_Air_temp(k,n)*1.84896D0      ! heat diffusivity (J/m/K/day), 1.84896D0=2.14D-5(heat diff)*86400.0D0       
       endif      
      
cccz      ! above water 
       do k=SubmergeIndex+1,mulchLayer 
          ! calculate the 'vapor conditions'
        RelaHumid_mulch_temp(k,n)=exp(2.124D-4*MulchElehNew_temp(k,n) ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))
         RelaHumid_mulch_temp(k,n)=min(RelaHumid_mulch_temp(k,n),1.0D0)! make sure the relative humidity is of currect range      
         RelaHumid_mulch_temp(k,n)=max(RelaHumid_mulch_temp(k,n),0.0D0)
        VaporSat_mulch_temp(k,n)=exp(19.84D0-4975.9D0     ! Saturated Vapor density for each mulch element 
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))              ! (g/m^3)           
        VaporAct_mulch_temp(k,n)=VaporSat_mulch_temp(k,n) ! Actual Vapor density for each mulch element
     &   *RelaHumid_mulch_temp(k,n)                       ! (g/m^3)
          ! water difficit in each element, including liquid water and water vapor to a certain degree of saturation
        WaterDifici_mulch_temp(k,n)=
     &    VaporSat_mulch_temp(k,n)-VaporAct_mulch_temp(k,n)
        WaterDifici_mulch_temp(k,n)=WaterDifici_mulch_temp(k,n)
     &    *(thickPerLayer(k)/100.0D0)/LocalStep            ! rate to saturate the vapor portion (g/m^2/day)
          ! water diffusivity calculation, use diffusion type water fluxes
        VaporDiff_mulch_temp(k,n)=2.29D-5*
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**0.667D0)*f_mulch_pore*86400.0D0              ! Water Vapor diffusivity (m^2/day)
          ! water capacity (only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore)
          ! use soil water capacity, but recommend to use wood water retention curve in future.           
        if(Van_Camp_mulch.eq.0) then  
         WaterCapaMulch_temp(k,n)=
     &     WC_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        elseif(Van_Camp_mulch.eq.1) then
         WaterCapaMulch_temp(k,n)=
     &     WC_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput)
         WaterCapaMulch_temp(k,n)=
     &     WaterCapaMulch_temp(k,n)*(1.0D0-f_mulch_pore)
        endif
          ! heat capacity and diffusivity calculation
          ! air density (rho_dryair) g m^-3, other parameters 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air_temp(k,n)=1.85D0*VaporAct_mulch_temp(k,n)
     &   +1.006D0*(rho_dryair-VaporAct_mulch_temp(k,n))                ! heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=HeatCapa_Air_temp(k,n)*1.84896D0      ! heat diffusivity (J/m/K/day), 1.84896D0=2.14D-5(heat diff)*86400.0D0      
       enddo
      enddo

cccz calculate the water fluxes within the mulch grid
cccz for water, we assume only vapor flux occur, but allow some liquid water stored in the mulch solid materials.
cccz downward and leftward are positive flux (water income)
cccz also assume there is no vapor flux under ponded water surface
      if(PerOccupation.eq.2.0D0) then
       ! the water level is high enough, this layer is now a water layer
       ! zero all the under water layers
       do k=1,SubmergeIndex+1
         do n=1,SurNodeIndex
            g_vapor(2*k-1,n)=0.0D0
            g_vapor(2*k,n)=0.0D0
         enddo
       enddo
       
       do k=SubmergeIndex+1,mulchLayer+1    
       if(k.eq.(SubmergeIndex+1)) then                  
        do n=1,SurNodeIndex-1
cccz the soil surface vapor flux based on the 'potential evaporation' method
cccz maybe need to adjust the inwards/outwards flux to soil
         SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(k,n))
     &    /(MulchEleTmpr_temp(k,n)+237.3D0))
         DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(k,n)+1.0D0))
     &    /(MulchEleTmpr_temp(k,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch_temp(k,n))          ! calculate VPD for the first mulching layer
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, 
                                                                     ! should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(k,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur
     &    *max(netRad(k,n),0.0D0)*3600.0D0/(2500.8D0
     &    -(2.3668D0*MulchEleTmpr_temp(k,n))))+(VPD_Sur*109.375D0
     &    *(1.0D0+(0.149D0*u_mulch(k)))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 

cccz the soil surface vapor flux based on the 'vapor flux' method
cccz first calcualte the dimensionless parameters
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         hNew_Sur=0.0D0
         VaporSat_Sur=exp(19.84D0-4975.9D0/(Tmpr_Sur+273.15D0))
         RelaHumid_Sur=1.0D0                                   ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
         VaporAct_Sur=VaporSat_Sur*RelaHumid_Sur
         
         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur)
     &    *(0.5D0*thickPerLayer(k)/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur)/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur)
     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur)
     &    /2.0D0)/((u_mulch(k)/3.6D0)**2.0D0)
     
         if(Ra_mulch(k,n).le.Ra_Critical
     &     .or.DiffusionRes.eq.1) then
           VaporCond_Mul=VaporDiff_mulch_temp(k,n)
     &      /(0.5D0*thickPerLayer(k)/100.0D0)            ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
            VaporCond_Mul=0.055D0*((abs(Tmpr_Sur-MulchEleTmpr_temp(k,n))
     &      /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &      +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &      **0.25D0)*86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0*sqrt((u_mulch(k)/3.6D0)
     &      /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &      +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &      *86400.0D0/41.4D0
          endif
         endif
         g_vapor_try2=VaporCond_Mul
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_Sur)
cccz take the min between 'diffusive type' and 'potential' evporation     
cccz if the two values are of different directions, we assume zero flux fo numerical stable
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then
           g_vapor(2*k-1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(2*k-1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(2*k-1,n)=0.0D0
         endif
        enddo      
        g_Vapor(2*k,1)=0.0D0           ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(k,n-1)
     &    +widthPerMulchUnit(n)/VaporDiff_mulch_temp(k,n))
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
        enddo  
        
cccz for mulch-air interface        
       elseif(k.eq.(mulchLayer+1)) then                               
        do n=1,SurNodeIndex-1                                        ! use conduction/convection to calculate conductivity at atmosphere surface
         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    *((z_wind-mulchThick+0.5D0*thickPerLayer(k-1))/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)
     &    /((u_above/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul=sqrt(2.40D-5/((z_wind-mulchThick
     &      +0.5D0*thickPerLayer(k-1))/100.0D0)*86400.0D0
     &      *VaporDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0)) ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(AirTemp_Wea
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_above/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_Ambient-VaporAct_mulch_temp(k-1,n))
        enddo  

cccz for water flux within mulch
       else                                                           
        do n=1,SurNodeIndex-1
         Ra_mulch(k,n)=9.81D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    *((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
     &    /nu_air/((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0
     &    *sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))**2.0D0
     &    +(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/((u_mulch(k)/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch_temp(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch_temp(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(MulchEleTmpr_temp(k,n)
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_mulch(k)/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch_temp(k,n))
     &     *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
       endif
       enddo
      else   
       ! the water level is not high enough, the occupied layer is now a fraction of the original layers
       ! zero all the under water layers, could be more than needed but not a problem, because we calcualte it later
       do k=1,SubmergeIndex+1
         do n=1,SurNodeIndex
            g_vapor(2*k-1,n)=0.0D0
            g_vapor(2*k,n)=0.0D0
         enddo
       enddo    
       do k=SubmergeIndex,mulchLayer+1    
       if(k.eq.SubmergeIndex) then
        ! take the fraction of the thickness
        thickLayer_Sur=thickPerLayer(k)*PerOccupation
        do n=1,SurNodeIndex-1
cccz the soil surface vapor flux based on the 'potential evaporation' method
cccz maybe need to adjust the inwards/outwards flux to soil
         SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(k,n))
     &    /(MulchEleTmpr_temp(k,n)+237.3D0))
         DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(k,n)+1.0D0))
     &    /(MulchEleTmpr_temp(k,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch_temp(k,n))          ! calculate VPD for the first mulching layer
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, 
                                                                     ! should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(k,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur*max(netRad(k,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*MulchEleTmpr_temp(k,n))))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_mulch(k)))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 

cccz the soil surface vapor flux based on the 'vapor flux' method
cccz first calcualte the dimensionless parameters
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         hNew_Sur=0.0D0                                         ! because it is on the water surface
         VaporSat_Sur=exp(19.84D0-4975.9D0/(Tmpr_Sur+273.15D0))
         RelaHumid_Sur=1.0D0                                    ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
         VaporAct_Sur=VaporSat_Sur*RelaHumid_Sur
         
         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur)
     &    *(0.5D0*thickLayer_Sur/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur)/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur)
     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur)
     &    /2.0D0)/((u_mulch(k)/3.6D0)**2.0D0)
     
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul=VaporDiff_mulch_temp(k,n)
     &      /(0.5D0*thickLayer_Sur/100.0D0)            ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
            VaporCond_Mul=0.055D0*((abs(Tmpr_Sur
     &      -MulchEleTmpr_temp(k,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0*sqrt((u_mulch(k)/3.6D0)
     &      /sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_vapor_try2=VaporCond_Mul
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_Sur)
cccz take the min between 'diffusive type' and 'potential' evporation     
cccz if the two values are of different directions, we assume zero flux fo numerical stable
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then
           g_vapor(2*k-1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(2*k-1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(2*k-1,n)=0.0D0
         endif
        enddo      
        g_Vapor(2*k,1)=0.0D0           ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(k,n-1)
     &    +widthPerMulchUnit(n)/VaporDiff_mulch_temp(k,n))
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
        enddo  
        
cccz for mulch-air interface        
       elseif(k.eq.(mulchLayer+1)) then                               
        do n=1,SurNodeIndex-1                                        ! use conduction/convection to calculate conductivity at atmosphere surface
         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    *((z_wind-mulchThick+0.5D0*thickPerLayer(k-1))/100.0D0)/nu_air
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &    **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)
     &    /((u_above/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul=sqrt(2.40D-5/((z_wind-mulchThick
     &      +0.5D0*thickPerLayer(k-1))/100.0D0)*86400.0D0
     &      *VaporDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0)) ! the vapor conductance (m/day)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(AirTemp_Wea
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_above/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_Ambient-VaporAct_mulch_temp(k-1,n))
        enddo  

cccz for water flux within mulch
       else                                                           
        do n=1,SurNodeIndex-1
         Ra_mulch(k,n)=9.81D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    *((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
     &    /nu_air/((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/Dhm
         GrRe_mulch(k,n)=9.81D0
     &    *sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))**2.0D0
     &    +(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)
     &    +MulchEleTmpr_temp(k,n))/2.0D0)/((u_mulch(k)/3.6D0)**2.0D0)
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch_temp(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch_temp(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           VaporCond_Mul=0.055D0*((abs(MulchEleTmpr_temp(k,n)
     &      -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &      *86400.0D0/41.4D0
          else
           VaporCond_Mul=0.147D0
     &      *sqrt((u_mulch(k)/3.6D0)/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)
     &      -slopeCoord(n,2))**2.0D0)*100.0D0)*86400.0D0/41.4D0
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul
     &    *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch_temp(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch_temp(k,n))
     &     *(VaporAct_mulch_temp(k,n)-VaporAct_mulch_temp(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
       endif
       enddo    
      endif
       
           
cccz calculate the rainfall redistribution
cccz this part is universal because in both 'PerOccupation', the 'WaterDifici_mulch_temp' values is adjusted
cccz such that the under water part, 'WaterDifici_mulch_temp=0.0D0' and the rainfall directly pass to the soil surface
      do n=1,SurNodeIndex-1
       k=mulchLayer+1
2308   if(k.gt.1) then
        inputPerLayer=
     &   min(RainFallInput_temp(k,n),WaterDifici_mulch_temp(k-1,n))
        if(inputPerLayer.lt.RainFallInput_temp(k,n)) then
         RainFallInput_temp(k-1,n)=RainFallInput_temp(k,n)-inputPerLayer
         RainFallInput_temp(k,n)=inputPerLayer
         k=k-1
         goto 2308
        else
         continue
        endif
       else
        continue
       endif
      enddo
       
cccz calculate the latent heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes based on the upwind direction of water      
      if(PerOccupation.eq.2.0D0) then
       ! the water level is high enough, this layer is now a water layer
       ! zero all the under water layers
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex
          g_Heat_Latent(2*k-1,n)=0.0D0
          g_Heat_Latent(2*k,n)=0.0D0
        enddo
       enddo
       do k=SubmergeIndex+1,mulchLayer+1
cccz for water surface layer
       if(k.eq.(SubmergeIndex+1)) then    
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
           qLeft=SurfNodeNodeIndexH(n)
           qRight=SurfNodeNodeIndexH(n+1)
           Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(Tmpr_Sur-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                             ! leftwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &      *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         else                                                         ! rightwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &      *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
cccz for mulch-air interface layer
       elseif(k.eq.(mulchLayer+1)) then
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
cccz for other interier layers
       else
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                              ! leftwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n) 
         else                                                          ! rightwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
       endif   
       enddo
      else   
       ! the water level is not high enough, the occupied layer is now a fraction of the original layers
       ! zero all the under water layers, maybe more than needed
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex
          g_Heat_Latent(2*k-1,n)=0.0D0
          g_Heat_Latent(2*k,n)=0.0D0
        enddo
       enddo
       do k=SubmergeIndex,mulchLayer+1
cccz for soil surface layer
       if(k.eq.SubmergeIndex) then    
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
           qLeft=SurfNodeNodeIndexH(n)
           qRight=SurfNodeNodeIndexH(n+1)
           Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(Tmpr_Sur-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                             ! leftwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &      *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         else                                                         ! rightwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &      *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
cccz for mulch-air interface layer
       elseif(k.eq.(mulchLayer+1)) then
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
cccz for other interier layers
       else
        do n=1,SurNodeIndex-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                              ! leftwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n) 
         else                                                          ! rightwards vapor flow
           g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
       endif   
       enddo
      endif

      
cccz calculate the sensible (diff) heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes along soil surface using the current soil surface temperature   
      if(PerOccupation.eq.2.0D0) then
       ! the water level is high enough, this layer is now a water layer
       ! zero all the under water layers
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex
         g_Heat_Sensi(2*k-1,n)=0.0D0
         g_Heat_Sensi(2*k,n)=0.0D0
        enddo
       enddo 
       do k=SubmergeIndex+1,mulchLayer+1  
cccz for mulch-soil interface
       if(k.eq.(SubmergeIndex+1)) then                          
        do n=1,SurNodeIndex-1
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul=HeatDiff_mulch_temp(k,n)
     &     /(thickPerLayer(k)/200.0D0)                               ! the heat conductance (J/m/day/K)
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           HeatCond_Mul=0.050D0*((abs(Tmpr_Sur-MulchEleTmpr_temp(k,n))
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0+
     &     (slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)**0.25D0)
     &     *HeatCapa_Air_temp(k,n)*86400.0D0/41.4D0
          else
           HeatCond_Mul=0.135D0*sqrt((u_mulch(k)/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *HeatCapa_Air_temp(k,n)*86400.0D0/41.4D0
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=
     &    HeatCond_Mul*(MulchEleTmpr_temp(k,n)-Tmpr_Sur)              ! the sensible heat conductance (J/m^2/day) 
        enddo
        g_Heat_Sensi(2*k,1)=0.0D0
        g_Heat_Sensi(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Heat_Sensi(2*k,n)=
     &    (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(k,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch_temp(k,n))
     &    *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
cccz for mulch-air interface 
       elseif(k.eq.(mulchLayer+1)) then                                 
        do n=1,SurNodeIndex-1
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul=sqrt(2.14D-5/((z_wind-mulchThick+0.05D0
     &      *thickPerLayer(k-1))/100.0D0)*HeatCapa_Air_temp(k-1,n)
     &      *86400.0D0*HeatDiff_mulch_temp(k-1,n)                     ! the heat conductance (J/m/day/K)
     &      /(thickPerLayer(k-1)/200.0D0))
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           HeatCond_Mul=0.050D0*((abs(AirTemp_Wea
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(k-1,n)
     &     *86400.0D0/41.4D0
          else
           HeatCond_Mul=0.135D0*sqrt((u_above/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *HeatCapa_Air_temp(k-1,n)*86400.0D0/41.4D0
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=HeatCond_Mul
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
cccz for mulch interior layers 
      else                                    
       do n=1,SurNodeIndex-1
        if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch_temp(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch_temp(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)      ! the heat conductance (J/m/day/K)
        else
         if(GrRe_mulch(k,n).ge.GrRe_Critical) then
          HeatCond_Mul=0.050D0*((abs(MulchEleTmpr_temp(k,n)
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*(HeatCapa_Air_temp(k-1,n)
     &     +HeatCapa_Air_temp(k,n))/2.0D0*86400.0D0/41.4D0
         else
          HeatCond_Mul=0.135D0*sqrt((u_mulch(k)/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *(HeatCapa_Air_temp(k-1,n)+HeatCapa_Air_temp(k,n))
     &     /2.0D0*86400.0D0/41.4D0
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=HeatCond_Mul*(MulchEleTmpr_temp(k,n)
     &   -MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex)=0.0D0
       do n=2,SurNodeIndex-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch_temp(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo 
       endif
      enddo
      else
       ! the water level is not that high enough within that layer
       ! zero all the under water layers
       do k=1,SubmergeIndex+1
         do n=1,SurNodeIndex
            g_Heat_Sensi(2*k-1,n)=0.0D0
            g_Heat_Sensi(2*k,n)=0.0D0
         enddo
       enddo    
       do k=SubmergeIndex,mulchLayer+1  
cccz for mulch-soil interface
        if(k.eq.SubmergeIndex) then    
        thickLayer_Sur=thickPerLayer(k)*PerOccupation
         do n=1,SurNodeIndex-1
          qLeft=SurfNodeNodeIndexH(n)
          qRight=SurfNodeNodeIndexH(n+1)
          Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
          if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
            HeatCond_Mul=HeatDiff_mulch_temp(k,n)
     &        /(thickLayer_Sur/200.0D0)                       ! the heat conductance (J/m/day/K)
          else
           if(GrRe_mulch(k,n).ge.GrRe_Critical) then
            HeatCond_Mul=0.050D0*((abs(Tmpr_Sur
     &      -MulchEleTmpr_temp(k,n))/sqrt((slopeCoord(n+1,1)
     &      -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &      **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(k,n)
     &      *86400.0D0/41.4D0
           else
            HeatCond_Mul=0.135D0*sqrt((u_mulch(SubmergeIndex+1)/3.6D0)
     &      /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &      +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &      *HeatCapa_Air_temp(k,n)*86400.0D0/41.4D0
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=HeatCond_Mul*
     &      (MulchEleTmpr_temp(k,n)-Tmpr_Sur)                 ! the sensible heat conductance (J/m^2/day) 
        enddo
        g_Heat_Sensi(2*k,1)=0.0D0
        g_Heat_Sensi(2*k,SurNodeIndex)=0.0D0
        do n=2,SurNodeIndex-1
         g_Heat_Sensi(2*k,n)=
     &    (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(k,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch_temp(k,n))
     &    *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
cccz for mulch-air interface 
       elseif(k.eq.(mulchLayer+1)) then                                 
        do n=1,SurNodeIndex-1
         if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul=sqrt(2.14D-5/((z_wind-mulchThick+0.05D0
     &      *thickPerLayer(k-1))/100.0D0)*HeatCapa_Air_temp(k-1,n)
     &      *86400.0D0                                                 ! the heat conductance (J/m/day/K)
     &      *HeatDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0))
         else
          if(GrRe_mulch(k,n).ge.GrRe_Critical) then
           HeatCond_Mul=0.050D0*((abs(AirTemp_Wea
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(k-1,n)
     &     *86400.0D0/41.4D0
          else
           HeatCond_Mul=0.135D0*sqrt((u_above/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *HeatCapa_Air_temp(k-1,n)*86400.0D0/41.4D0
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=HeatCond_Mul
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
cccz for mulch interior layers 
      else                                    
       do n=1,SurNodeIndex-1
        if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch_temp(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch_temp(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)      ! the heat conductance (J/m/day/K)
        else
         if(GrRe_mulch(k,n).ge.GrRe_Critical) then
          HeatCond_Mul=0.050D0*((abs(MulchEleTmpr_temp(k,n)
     &     -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &     -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &     **2.0D0)*100.0D0)**0.25D0)*(HeatCapa_Air_temp(k-1,n)
     &     +HeatCapa_Air_temp(k,n))/2.0D0*86400.0D0/41.4D0
         else
          HeatCond_Mul=0.135D0*sqrt((u_mulch(k)/3.6D0)
     &     /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &     +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &     *(HeatCapa_Air_temp(k-1,n)+HeatCapa_Air_temp(k,n))
     &     /2.0D0*86400.0D0/41.4D0
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=HeatCond_Mul*(MulchEleTmpr_temp(k,n)
     &   -MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex)=0.0D0
       do n=2,SurNodeIndex-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch_temp(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch_temp(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
       enddo 
       endif
      enddo
      endif
             
cccz finish the calculation of vapor and heat fluxes
cccz start to establish the equation systems which solve the temperature and water potential
      if(PerOccupation.eq.2.0D0) then
       ! the water level is high enough, this layer is now a water layer
       ! fix all the values under water surface
       do k=1,SubmergeIndex
        do n=1,SurNodeIndex-1
         MulchElehNew_temp2(k,n)=0.0D0
         MulchElehNew_temp0(k,n)=0.0D0
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         MulchEleTmpr_temp2(k,n)=Tmpr_Sur
         MulchEleTmpr_temp0(k,n)=Tmpr_Sur
        enddo
       enddo
       do n=1,SurNodeIndex-1
        do k=SubmergeIndex+1,mulchLayer
         a11=WaterCapaMulch_temp(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
         a12=f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
         a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0)              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    -WaterCapaMulch_temp(k,n)*VaporAct_mulch_temp(k,n))
     &    +c_water_s*rho_w*WaterCapaMulch_temp(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
         a22=c_air_s*rho_dryair+c_water_s*rho_w*MulchEleThNew_temp(k,n)
     &    +c_vap_s*VaporAct_mulch_temp(k,n)+(c_vap_v+c_vap_s
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref))*f_mulch_pore
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
         b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +RainFallInput_temp(k+1,n)/(thickPerLayer(k)/100.0D0)
     &    *LocalStep
         b12=(netRadEachLayer(k,n)*86400.0D0              
     &    +(g_Heat_Sensi(2*k+1,n)-g_Heat_Sensi(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Sensi(2*k,n+1)-g_Heat_Sensi(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0)
     &    +(g_Heat_Latent(2*k+1,n)-g_Heat_Latent(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Latent(2*k,n+1)-g_Heat_Latent(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +c_water_s*RainFallInput_temp(k+1,n)*(AirTemp_Wea-Tmpr_ref)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep
      ! Cramer's rule for linear system
         SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)
         SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
         MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
         MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
        enddo
       enddo
      else
       do k=1,SubmergeIndex-1
        do n=1,SurNodeIndex-1
         MulchElehNew_temp2(k,n)=0.0D0
         MulchElehNew_temp0(k,n)=0.0D0
         qLeft=SurfNodeNodeIndexH(n)
         qRight=SurfNodeNodeIndexH(n+1)
         Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
         MulchEleTmpr_temp2(k,n)=Tmpr_Sur
         MulchEleTmpr_temp0(k,n)=Tmpr_Sur
        enddo
       enddo    
       k=SubmergeIndex
       thickLayer_Sur=thickPerLayer(k)*PerOccupation
       do n=1,SurNodeIndex-1
       ! for k=SubmergeIndex
        a11=WaterCapaMulch_temp(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0)              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    -WaterCapaMulch_temp(k,n)*VaporAct_mulch_temp(k,n))
     &    +c_water_s*rho_w*WaterCapaMulch_temp(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_air_s*rho_dryair+c_water_s*rho_w*MulchEleThNew_temp(k,n)
     &    +c_vap_s*VaporAct_mulch_temp(k,n)+(c_vap_v+c_vap_s
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref))*f_mulch_pore
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickLayer_Sur/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &   +RainFallInput_temp(k+1,n)/(thickLayer_Sur/100.0D0)*LocalStep
        b12=(netRadEachLayer(k,n)*86400.0D0              
     &    +(g_Heat_Sensi(2*k+1,n)-g_Heat_Sensi(2*k-1,n))
     &    /(thickLayer_Sur/100.0D0)
     &    +(g_Heat_Sensi(2*k,n+1)-g_Heat_Sensi(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0)
     &    +(g_Heat_Latent(2*k+1,n)-g_Heat_Latent(2*k-1,n))
     &    /(thickLayer_Sur/100.0D0)
     &    +(g_Heat_Latent(2*k,n+1)-g_Heat_Latent(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +c_water_s*RainFallInput_temp(k+1,n)*(AirTemp_Wea-Tmpr_ref)
     &    /(thickLayer_Sur/100.0D0)*LocalStep
      ! Cramer's rule for linear system
        SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)                    
        SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
       enddo 
       ! for other layer
       do k=SubmergeIndex+1,mulchLayer
        do n=1,SurNodeIndex-1
        a11=WaterCapaMulch_temp(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_temp(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0)              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    -WaterCapaMulch_temp(k,n)*VaporAct_mulch_temp(k,n))
     &    +c_water_s*rho_w*WaterCapaMulch_temp(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_air_s*rho_dryair+c_water_s*rho_w*MulchEleThNew_temp(k,n)
     &    +c_vap_s*VaporAct_mulch_temp(k,n)+(c_vap_v+c_vap_s
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref))*f_mulch_pore
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &   +RainFallInput_temp(k+1,n)/(thickPerLayer(k)/100.0D0)*LocalStep
        b12=(netRadEachLayer(k,n)*86400.0D0              
     &    +(g_Heat_Sensi(2*k+1,n)-g_Heat_Sensi(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Sensi(2*k,n+1)-g_Heat_Sensi(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0)
     &    +(g_Heat_Latent(2*k+1,n)-g_Heat_Latent(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)
     &    +(g_Heat_Latent(2*k,n+1)-g_Heat_Latent(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +c_water_s*RainFallInput_temp(k+1,n)*(AirTemp_Wea-Tmpr_ref)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep
      ! Cramer's rule for linear system
        SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)                    
        SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
       enddo
      enddo
      endif
           
cccz finish the tempted calculation, start to check the convergence
      maxhNew=0.0D0
      maxTmpr=0.0D0
      maxhNewDiff=0.0D0
      maxTmprDiff=0.0D0
      if(IterMulch.le.mulchLayer) then                            ! the Picard iteration has to be processed for times >= mulch layer
       IterMulch=IterMulch+1
       do n=1,SurNodeIndex-1
        do k=1,mulchLayer
cccz protection
        if(MulchElehNew_temp2(k,n).gt.0.0D0) goto 2307
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) goto 2307  ! NAN case occurs
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) goto 2307  ! NAN case occurs
cccz protection
         MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
         MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
        enddo
cccz water content calculation
          ! do water calculation
        do k=1,SubmergeIndex-1
          MulchEleThNew_temp(k,n)=1.0D0
        enddo
        k=SubmergeIndex
        if(PerOccupation.eq.2.0D0) then
cccz we avoid calculate the water content for water-filled layer
          MulchEleThNew_temp(k,n)=1.0D0
        else
         if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         endif   
        endif
        do k=SubmergeIndex+1,mulchLayer
         if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         endif 
        enddo
       enddo
       goto 2001
      else
       IF(IterMulch.le.MaxIter) then
        do n=1,SurNodeIndex-1
        do k=1,mulchLayer
cccz protection
        if(MulchElehNew_temp2(k,n).gt.0.0D0) goto 2307
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) goto 2307  ! NAN case occurs
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) goto 2307  ! NAN case occurs
cccz protection
        MulchElehNew_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchElehNew_temp3(k,n)
     &     +MulchElehNew_temp2(k,n))/dble(IterMulch-mulchLayer)
        MulchEleTmpr_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchEleTmpr_temp3(k,n)
     &     +MulchEleTmpr_temp2(k,n))/dble(IterMulch-mulchLayer)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp3(k,n)
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp3(k,n)
        maxhNewDiff=max(maxhNewDiff, 
     &    abs(MulchElehNew_temp(k,n)-MulchElehNew_temp2(k,n)))
        maxTmprDiff=max(maxTmprDiff,
     &    abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp2(k,n)))
        maxhNew=max(maxhNew,abs(MulchElehNew(k,n)))
        maxTmpr=max(maxTmpr,abs(MulchEleTmpr(k,n)))  
        MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
        MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
        enddo
cccz water content calculation
        do k=1,SubmergeIndex-1
          MulchEleThNew_temp(k,n)=1.0D0
        enddo
        k=SubmergeIndex
        if(PerOccupation.eq.2.0D0) then
cccz we avoid calculate the water content for water-filled layer
          MulchEleThNew_temp(k,n)=1.0D0
        else
         if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         endif   
        endif
        do k=SubmergeIndex+1,mulchLayer
         if(Van_Camp_mulch.eq.0) then
          MulchEleThNew_temp(k,n)=
     &      WQ_VANG_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         elseif(Van_Camp_mulch.eq.1) then
          MulchEleThNew_temp(k,n)=
     &      WQ_CAMP_MULCH(MulchElehNew_temp(k,n),ParMulch(:,1),lInput) 
          MulchEleThNew_temp(k,n)=MulchEleThNew_temp(k,n)
     &      *(1.0D0-f_mulch_pore)
         endif 
        enddo
        enddo
        IterMulch=IterMulch+1
        if(max(maxhNewDiff/maxhNew,maxTmprDiff/maxTmpr).gt.errTol) then
         goto 2001
        else
         TotalTime=TotalTime+LocalStep
         if(TotalTime-Step.gt.-0.01D0*dtMin) then
           if(TimeShrink.eq.0) then   ! if no 'local time step' was used, 'one step' was used for calculating water/energy exchange on soil surface
            goto 2302
           else                       ! if multiple 'local time step' was used, then 'cumulative averaged' water/energy exchange were calculated in 2003 and 2004
            goto 2303
           endif
         else
          do n=1,SurNodeIndex-1
           do k=1,mulchLayer
            MulchElehNew_temp0(k,n)=MulchElehNew_temp(k,n)
            MulchEleTmpr_temp0(k,n)=MulchEleTmpr_temp(k,n)
           enddo
          enddo
          goto 2303
2304      LocalStep=min(LocalStep*DMul1,Step-TotalTime)
cccz??? ask a question, why not reset the mulch iteration index
cccz let us reset it for now
          IterMulch=1
c          IterMulch=1
          goto 2001
         endif
        endif
       ELSE
2307    TimeShrink=1
cccz no update for water content
         do n=1,SurNodeIndex-1
          do k=1,mulchLayer
          MulchElehNew_temp(k,n)=MulchElehNew_temp0(k,n)
          MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp0(k,n)
          enddo
         enddo
         LocalStep=max(LocalStep*DMul2,0.0001D0*dtMin)
         IterMulch=1
         goto 2001
       ENDIF
      endif
   
cccz finish the Picard iteration, but we have to update the water and energy fluxes across soil surface
cccz here we solve for the potential vapor flux at the soil surface
cccz if the actural one match the potential one, we do not need to solve the equation system later in the post process codes
cccz if not use the actural one and re-solve the equation system.
cccz first do the vapor flux calculation
2302  if(PerOccupation.eq.2.0D0) then
          kk=2*SubmergeIndex+1
      else
          kk=2*SubmergeIndex-1
      endif
      do n=1,SurNodeIndex-1
cccz the unit within this module will be g/day/m^2
cccz the unit taken in watermov module is g/day/cm^2
       if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
cccz update 'Varbw_Mulch' for 'SurWater'
         ! all the water will pass to the soil surface
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,n)/10000.0D0
         ! the evaporation is calculated on certain height
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,n)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3) 
       elseif(n.eq.(SurNodeIndex-1)) then
         kSurL=SurfNodeSurfIndexH(n)
         kSurR=SurfNodeSurfIndexH(n+1)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &     -(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &     +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &     /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)
         
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,n)/10000.0D0         
         Varbw_Mulch(kSurR,2)=-g_Vapor(kk,n)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)
       else
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)  
       endif 
      enddo   
      
cccz second update the heat fluxes
cccz the unit within this module will be J/day/m^2
cccz ********* ideas for the parameter assignment **************
c    VarBT(,1) is the temperature, as assign the temperature at mulch bottom surface (not important, I need fluxes more)
c    VarBT(,2) and VarBT(,3) are coupled, VarBT(,2) looks like a coefficient and VarBT(,3)=VarBT(,2)*temperature
c         The only temperature dependent thing is sensible heat in this code, so
c         VarBT(,2) sensible heat conductance,  VarBT(,3)= VarBT(,2)*(temperature at mulch bottom surface)
c    VarBT(,4) is the net radiation, we have short wave and long wave radiation here.
c    Latent heat is not presented in VarBT, it will depend on VarBW and soil temperature.
c    Make unit change to Cal/cm^2/min as shown in hourwea.for
cccz ********* end ****************
      if(PerOccupation.eq.2.0D0) then
       kk=SubmergeIndex+1
      else
       kk=SubmergeIndex
      endif
      do n=1,SurNodeIndex-1
       if(n.eq.1) then
        kSurL=SurfNodeSurfIndexH(n)
        VarBT(kSurL,1)=MulchEleTmpr_temp(kk,n)
        VarBT(kSurL,2)=HeatDiff_mulch_temp(kk,n)*0.025104D0          ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*MulchEleTmpr_temp(kk,n)
        VarBT(kSurL,4)=netRad(kk,n)*8.64D0                           ! 8.64D0=86400D0/10000D0
       elseif(n.eq.(SurNodeIndex-1)) then
        kSurL=SurfNodeSurfIndexH(n)
        kSurR=SurfNodeSurfIndexH(n+1)
        VarBT(kSurL,1)=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        VarBT(kSurL,2)=(HeatDiff_mulch_temp(kk,n-1)
     &    *widthPerMulchUnit(n-1)+HeatDiff_mulch_temp(kk,n)
     &    *widthPerMulchUnit(n))/(widthPerMulchUnit(n-1)
     &    +widthPerMulchUnit(n))*0.025104D0                          ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
        VarBT(kSurL,4)=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0

        VarBT(kSurR,1)=MulchEleTmpr_temp(kk,n)
        VarBT(kSurR,2)=HeatDiff_mulch_temp(kk,n)*0.025104D0            ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(kk,n)
        VarBT(kSurR,4)=netRad(kk,n)*8.64D0                            ! 8.64D0=86400D0/10000D0
       else
        kSurL=SurfNodeSurfIndexH(n)
        VarBT(kSurL,1)=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        VarBT(kSurL,2)=(HeatDiff_mulch_temp(kk,n-1)
     &    *widthPerMulchUnit(n-1)+HeatDiff_mulch_temp(kk,n)
     &    *widthPerMulchUnit(n))/(widthPerMulchUnit(n-1)
     &    +widthPerMulchUnit(n))*0.025104D0                           ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
        VarBT(kSurL,4)=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
       endif
      enddo
      
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput_temp(n,k)
       enddo
      enddo
      goto 2306           ! to the final assignment
      
cccz need to estimate the soil-mulch interface water and energy flux for each time segment during the iteration
cccz use the same scheme showed following Index 2002
2303  if(PerOccupation.eq.2.0D0) then
       kk=SubmergeIndex+1
      else
       kk=SubmergeIndex
      endif
      do n=1,SurNodeIndex-1
       if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
         Local_VarBW2=-g_Vapor(kk,n)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
       elseif(n.eq.(SurNodeIndex-1)) then
         kSurL=SurfNodeSurfIndexH(n)
         kSurR=SurfNodeSurfIndexH(n+1)
         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW2=-(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
         
         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
         Local_VarBW2=-g_Vapor(kk,n)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n+1)=VarBW1_temp(n+1)+Local_VarBW1*LocalStep
         VarBW2_temp(n+1)=VarBW2_temp(n+1)+Local_VarBW2*LocalStep
         VarBW3_temp(n+1)=VarBW3_temp(n+1)+Local_VarBW3*LocalStep
         Q_temp(n+1)=Q_temp(n+1)+(-Width(kSurR)*Local_VarBW3)*LocalStep
       else
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW2=-(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Local_VarBW3=Local_VarBW2-VarBW(kSurL,1)
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
       endif 
      enddo
      do n=1,SurNodeIndex-1
      if(n.eq.1) then
         Local_VarBT1=MulchEleTmpr_temp(kk,n)
         Local_VarBT2=HeatDiff_mulch_temp(kk,n)*0.025104D0           ! 0.025104D0=4.184D0*60.0D0/10000.0D0
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=netRad(kk1,n)*8.64D0                           ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
       elseif(n.eq.(SurNodeIndex-1)) then
        Local_VarBT1=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        Local_VarBT2=(HeatDiff_mulch_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
        VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
        VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
        VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        Local_VarBT1=MulchEleTmpr_temp(kk,n)
        Local_VarBT2=HeatDiff_mulch_temp(kk,n)*0.025104D0            ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=netRad(kk,n)*8.64D0                            ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n+1)=VarBT1_temp(n+1)+Local_VarBT1*LocalStep
        VarBT2_temp(n+1)=VarBT2_temp(n+1)+Local_VarBT2*LocalStep
        VarBT3_temp(n+1)=VarBT3_temp(n+1)+Local_VarBT3*LocalStep
        VarBT4_temp(n+1)=VarBT4_temp(n+1)+Local_VarBT4*LocalStep
       else
        Local_VarBT1=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
        Local_VarBT2=(HeatDiff_mulch_temp(kk,n-1)*widthPerMulchUnit(n-1)
     &    +HeatDiff_mulch_temp(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*0.025104D0   ! 0.025104D0=4.184D0*60.0D0/10000.0D0
        Local_VarBT3=Local_VarBT2*Local_VarBT1
        Local_VarBT4=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
     &    +netRad(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
        VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
        VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
        VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
        VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
       endif
      enddo
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)
     &   +RainFallInput_temp(n,k)*LocalStep
       enddo
      enddo
      if(TotalTime-Step.gt.-0.01D0*dtMin) then
          goto 2305
      else
          goto 2304
      endif
      

cccz use the 'cumulative averaged' water/energy exchange between 2003 and 2004
cccz determine the vapor and energy exchange on soil surface
2305  do n=1,SurNodeIndex
cccz water vapor part
       kSur=SurfNodeSurfIndexH(n)
       Varbw_Mulch(kSur,1)=VarBW1_temp(n)/TotalTime
       Varbw_Mulch(kSur,2)=VarBW2_temp(n)/TotalTime
       Varbw_Mulch(kSur,3)=VarBW3_temp(n)/TotalTime
       VarBW(kSur,1)=Varbw_Mulch(kSur,1)
       VarBW(kSur,2)=Varbw_Mulch(kSur,2)
       VarBW(kSur,3)=Varbw_Mulch(kSur,3)
       nNode=KXB(kSur)
       Q(nNode)=Q_temp(n)/TotalTime
cccz energy part
       VarBT(kSur,1)=VarBT1_temp(n)/TotalTime
       VarBT(kSur,2)=VarBT2_temp(n)/TotalTime
       VarBT(kSur,3)=VarBT3_temp(n)/TotalTime
       VarBT(kSur,4)=VarBT4_temp(n)/TotalTime
      enddo
      do k=1,SurNodeIndex-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)/TotalTime
       enddo
      enddo
      goto 2306


2306  aaaa=2
cccz final assignment
      do n=1,SurNodeIndex-1
       do k=1,mulchLayer
         MulchElehNew(k,n)=MulchElehNew_temp(k,n)
         MulchEleTmpr(k,n)=MulchEleTmpr_temp(k,n)
         MulchEleThNew(k,n)=MulchEleThNew_temp(k,n)
       enddo
      enddo
      goto 2500

cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz ----------------------------------------------------------------------------------------------------------
cccz the mulch is totally submerged in surface water
cccz in this case, we assume 
c       (1) the soil surfae temp and ponded water temp are the same
c       (2) the mulch is totally saturated
c       (3) we only need to adjust the varBW2 (evaporation) for all varBW and varBT quantities
c       (4) no Picard iteration is needed.
cccz make assignment directly
2400  do n=1,SurNodeIndex-1
       qLeft=SurfNodeNodeIndexH(n)
       qRight=SurfNodeNodeIndexH(n+1)
       Tmpr_Sur=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))
       do k=1,mulchLayer 
         MulchElehNew(k,n)=0.0D0
         MulchElehNew_temp(k,n)=0.0D0
         MulchEleTmpr(k,n)=Tmpr_Sur
         MulchEleTmpr_temp(k,n)=Tmpr_Sur
         MulchEleThNew(k,n)=1.0D0
         MulchEleThNew_temp(k,n)=1.0D0
       enddo
      enddo
      
cccz then make some possible flux adjustments
      do n=1,SurNodeIndex-1
       do k=1,mulchLayer 
cccz calculate the 'vapor conditions'
        RelaHumid_mulch_temp(k,n)=1.0D0
        VaporSat_mulch_temp(k,n)=1.0D06                   ! water density (g/m^3)           
        VaporAct_mulch_temp(k,n)=1.0D06                       ! (g/m^3)
cccz water difficit in each element
        WaterDifici_mulch_temp(k,n)=0.0D0
cccz water diffusivity calculation
        VaporDiff_mulch_temp(k,n)=0.0D0           ! Water Vapor diffusivity (m^2/day)
cccz water capacity
        WaterCapaMulch_temp(k,n)=0.0D0
cccz heat capacity and diffusivity  
        HeatCapa_Air_temp(k,n)=0.0D0                ! there is no air for heat capacity (J/m^3/K)
        HeatDiff_mulch_temp(k,n)=0.0D0              ! there is no air for heat diffusivity (J/m/K/day)
       enddo
      enddo

cccz calculate the water fluxes within the mulch grid
cccz for water, we assume only vapor flux occur, but allow some liquid water stored in the mulch solid materials.
cccz downward and leftward are positive flux (water income)
      do k=1,mulchLayer
       do n=1,SurNodeIndex
        g_vapor(2*k-1,n)=0.0D0
        g_vapor(2*k,n)=0.0D0
       enddo
      enddo
cccz for water surface
      k=mulchLayer+1                                               
      do n=1,SurNodeIndex-1
cccz the soil surface vapor flux based on the 'potential evaporation' method
cccz maybe need to adjust the inwards/outwards flux to soil
        SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(k-1,n))
     &   /(MulchEleTmpr_temp(k-1,n)+237.3D0))
        DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(k-1,n)+1.0D0))
     &   /(MulchEleTmpr_temp(k-1,n)+1.0D0+237.3D0)))-SVPA_Sur
        VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch_temp(k-1,n))          ! calculate VPD for the first mulching layer
        D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, 
                                                                     ! should be the saturated vapor pressure of wet bulb temp in air
        D32=2500.8D0-2.37D0*MulchEleTmpr_temp(k-1,n)
        GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
        g_Vapor(2*k-1,n)=-((DEL_Sur/GAMMA_Sur*max(netRad(k,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*MulchEleTmpr_temp(k-1,n))))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_soilsur))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 
      enddo
     
cccz calculate the rainfall redistribution
cccz everything go to soil surface
      do n=1,SurNodeIndex-1
       RainFallInput_temp(1,n)=RainFallInput_temp(mulchLayer+1,n)
       do k=2,mulchLayer+1
         RainFallInput_temp(k,n)=0.0D0
       enddo
      enddo
           
cccz calculate the latent heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes based on the upwind direction of water
      do k=1,mulchLayer
       do n=1,SurNodeIndex
        g_Heat_Latent(2*k-1,n)=0.0D0
        g_Heat_Latent(2*k,n)=0.0D0
       enddo
      enddo
      k=mulchLayer+1
cccz for mulch-air interface layer
      do n=1,SurNodeIndex-1
       if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
         g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &    *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
       else                                                         ! upwards vapor flow
         g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &    *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
       endif
      enddo
      
cccz calculate the sensible (diff) heat flux at each boundary
cccz downward and leftward are positive flux (energy income)
cccz calculate the fluxes along soil surface using the current soil surface temperature  
      do k=1,mulchLayer
       do n=1,SurNodeIndex
        g_Heat_Sensi(2*k-1,n)=0.0D0
        g_Heat_Sensi(2*k,n)=0.0D0
       enddo
      enddo
      k=mulchLayer+1      
      do n=1,SurNodeIndex-1
       Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &   *((z_wind-mulchThick+0.5D0*thickPerLayer(k-1))/100.0D0)/nu_air
     &   /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
       GrRe_mulch(k,n)=9.81D0*sqrt((slopeCoord(n,1)-slopeCoord(n+1,1))
     &   **2.0D0+(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)/100.0D0
     &   *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &   /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)
     &   /((u_above/3.6D0)**2.0D0)
       if(Ra_mulch(k,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul=sqrt(2.14D-5/((z_wind-mulchThick+0.05D0
     &    *thickPerLayer(k-1))/100.0D0)*HeatCapa_Air_temp(k-1,n)
     &    *86400.0D0                                                   ! the heat conductance (J/m/day/K)
     &    *HeatDiff_mulch_temp(k-1,n)/(thickPerLayer(k-1)/200.0D0))
       else
        if(GrRe_mulch(k,n).ge.GrRe_Critical) then
          HeatCond_Mul=0.050D0*((abs(AirTemp_Wea
     &    -MulchEleTmpr_temp(k-1,n))/sqrt((slopeCoord(n+1,1)
     &    -slopeCoord(n,1))**2.0D0+(slopeCoord(n+1,2)-slopeCoord(n,2))
     &    **2.0D0)*100.0D0)**0.25D0)*HeatCapa_Air_temp(k-1,n)
     &    *86400.0D0/41.4D0
        else
          HeatCond_Mul=0.135D0*sqrt((u_above/3.6D0)
     &    /sqrt((slopeCoord(n+1,1)-slopeCoord(n,1))**2.0D0
     &    +(slopeCoord(n+1,2)-slopeCoord(n,2))**2.0D0)*100.0D0)
     &    *HeatCapa_Air_temp(k-1,n)*86400.0D0/41.4D0
        endif
       endif
       g_Heat_Sensi(2*k-1,n)=HeatCond_Mul
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
      enddo
    
cccz here we solve for the potential vapor flux at the soil surface
      kk=2*mulchLayer+1
      do n=1,SurNodeIndex-1
cccz the unit within this module will be g/day/m^2
cccz the unit taken in watermov module is g/day/cm^2
       if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,n)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,n)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3) 
       elseif(n.eq.(SurNodeIndex-1)) then
         kSurL=SurfNodeSurfIndexH(n)
         kSurR=SurfNodeSurfIndexH(n+1)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)
         
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,n)/10000.0D0
         Varbw_Mulch(kSurR,2)=-g_Vapor(kk,n)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)
       else
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)  
       endif 
      enddo   
      
cccz we do not assume varBT changes
      goto 2500
     
cccz convert element based values to node based values
2500     do n=1,numMulchNode
          MulchNodehNew(n)=0.0D0
          MulchNodeThNew(n)=0.0D0
          MulchNodeTmpr(n)=0.0D0
        enddo
        do e=1,numMulchEle
          n1=MulchEleMarker(e,1)
          n2=MulchEleMarker(e,2)
          n3=MulchEleMarker(e,3)
          n4=MulchEleMarker(e,4)
          kk=MulchEleMarker(e,6)
          jj=MulchEleMarker(e,5)
          MulchNodehNew(n1)=MulchNodehNew(n1)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchElehNew(kk,jj)
          MulchNodehNew(n2)=MulchNodehNew(n2)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchElehNew(kk,jj)
          MulchNodehNew(n3)=MulchNodehNew(n3)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchElehNew(kk,jj)
          MulchNodehNew(n4)=MulchNodehNew(n4)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchElehNew(kk,jj)
          MulchNodeThNew(n1)=MulchNodeThNew(n1)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleThNew(kk,jj)
          MulchNodeThNew(n2)=MulchNodeThNew(n2)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleThNew(kk,jj)
          MulchNodeThNew(n3)=MulchNodeThNew(n3)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleThNew(kk,jj)
          MulchNodeThNew(n4)=MulchNodeThNew(n4)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleThNew(kk,jj)
          MulchNodeTmpr(n1)=MulchNodeTmpr(n1)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleTmpr(kk,jj)
          MulchNodeTmpr(n2)=MulchNodeTmpr(n2)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleTmpr(kk,jj)
          MulchNodeTmpr(n3)=MulchNodeTmpr(n3)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleTmpr(kk,jj)
          MulchNodeTmpr(n4)=MulchNodeTmpr(n4)
     &     +0.25D0*MulchEleMarkerArea(e)*MulchEleTmpr(kk,jj)
        enddo
        do n=1,numMulchNode
         MulchNodehNew(n)=MulchNodehNew(n)/MulchNodeMarkerArea(n)
         MulchNodeThNew(n)=MulchNodeThNew(n)/MulchNodeMarkerArea(n)
         MulchNodeTmpr(n)=MulchNodeTmpr(n)/MulchNodeMarkerArea(n)
        enddo
cccz finishe the calculation          
      return
2106  Write(*,*) 'Mulch file not found'      
2100  Call errmes(im,il)
      Return  
      END 