cccz -------------------------------------------      
c     Radiation Adjustment
c     Wind Speed Calculation
c     Vapor Flux
c     Convection Heat Flux, sensi_heat
c     Latent Heat Flux, latent_heat
cccz -------------------------------------------

cccz -------------------------------------------
c     variable used in mulch related module (some of the variables are defined in public ins)
c     scale: maximal layer NumMulHLD; maximal horizontal nodes NumBPD
c     "..,.." or "..*.." means "NumMulHLD" "NumBPD"
c     note: GEO: geometrical property, PHY: physical property, BIO: biological property, CON: Control property; 
c     unit: boolean, index, or physical units
c
c     INTEGER TYPE
c         DiffusionRes,               [CON]   [0-1]       "DiffusionRes"=1 force to use diffusion for flux calculation; "DiffusionRes"=0 use diffusion-convection equation.
c         ForceForwards,              [CON]   [0-1]       When "layer_need_fix", "ForceForwards"=1 is a flag to jump out of Picard interation at the min time step, force a forward step.
c         IterMulch,                  [CON]   [1]         counting Picard iteration step and compare with max step "MaxIter"
c         MaxIter,                    [CON]   [1]         max number of Picard iteration 
c         TimeShrink,                 [CON]   [0-1]       indicator if "local time" was used, or "refined time step" was used. This variable invoke an alternative way to assign flux variables, "varbw/varbt"
c         LongWaveRadCtrl,            [CON]   [0-1-2]       indicator to force the "long wave radiation to be 0". Note the longwave radiation is always calculated, "LongWaveRadCtrl" invoke a process to assign them back to 0.
c         mulchLayer,                 [GEO]   [1]         number of mulch layer
c         numMulchNode,               [GEO]   [1]         number of mulch node
c         numMulchEle,                [GEO]   [1]         number of mulch element
c         MulchNodeMatrix(..,..),     [GEO]   [1]         Convert (vert,hori) location to numMulchNode: the first index is the vertical position; the second index is the horizontal position; (1,1) is the node in lower-left corner; report the node index
c         MulchEleMarker(..*..,6),    [GEO]   [1]         mulch element array 1~4 mulch node index, 5 horizontal position, 6 vertical position (1 is soil surface, mulchLayer+1 is mulch surface)
c         MulchEleMatrix(..,..),      [GEO]   [1]         Convert (vert,hori) location to numMulchEle: the first index is the vertical position; the second index is the horizontal position; (1,1) is the node in lower-left corner; report the element index
c         SurNodeIndex,               [GEO]   [number]    the number of nodes along soil surface
c         qLeft, qRight,              [GEO]   [index]     indices to left/right surface nodes, used in iteration
c         SurfNodeNodeIndexH(..),     [GEO]   [array]     index of surface node in NodeArray, from one side to another side
c         SurfNodeSurfIndexH(..),     [GEO]   [array]     index of surface node in SurfaceNodeArray, from one side to another side
c         SubmergeIndex,              [GEO]   [index]     0: ponded water not existed; >=1: the top layer (element based) that being filled or partially filled by ponded water 
c         layer_need_fix,             [GEO]   [index]     @ that layer, even at the smallest timestep, the diffusion eqution does not work, so we have to push a forward step at that layer, it is marker
c         DecompCtrl,                 [CON]   [0-1]       control if mulch decomposition is envoked, 0: by_pass, 1: envoked
c
c     FLOAT TYPE
c         maxhNewDiff,                [CON]   [cm]        max diff of water potential between to Picard iteration steps, used for establishing converging test in Picard iteration 
c         maxTmprDiff,                [CON]   [oC]        max diff of temperature between to Picard iteration steps, used for establishing converging test in Picard iteration
c         maxhNew,                    [CON]   [cm]        max (absolute) water potential in the current Picard iteration steps, used for establishing converging test in Picard iteration
c         maxTmpr,                    [CON]   [oC]        max temperature in the current Picard iteration steps, used for establishing converging test in Picard iteration
c         errTol,                     [CON]   [%]         error tolerance for Picard iteration
c         mulchThick,                 [GEO]   [cm]        the whole mulch thickness
c         totalMulchWidth,            [GEO]   [cm]        the whole mulch width
c         f_mulch_pore,               [GEO]   [0~1]       mulch "air" space, NOT the porespace within the mulch solid part
c         widthPerMulchUnit(..),      [GEO]   [cm]        HORIZONTAL width of each mulch section
c         thickPerLayer(..),          [GEO]   [cm]        thick per mulch layer
c         thresholdThick,             [GEO]   [cm]        threshold for smallest layer thick, empirical values, need to merging two neighbor layers or think that layer is fully filled by ponded water
c         thickPerLayer_temp(..),     [GEO]   [cm]        temporary array for mulch layer thickness, used for merging layers
c         LayerHeight(..),            [GEO]   [cm]        the height of each mulch layer counted from soil-mulch interface 
c         mulchNodeCoord(..*..,2),    [GEO]   [cm,cm]     the coordinate of nodes along mulch-soil interface, 1 horizontal dir 2 vertical dir (cm)
c         MulchNodeMarkerArea(..*..), [GEO]   [cm2]       nodal area for mulch grid
c         MulchEleMarkerArea(..*..),  [GEO]   [cm2]       element area for mulch grid
c         slopeCoord(..*..,2),        [GEO]   [cm,cm]     coordinates, 1 horizontal (x), 2 vertical (y) of nodes along soil surface-mulch interface
c         MulchNodeTmpr(..*..),       [PHY]   [oC]        nodal temperature; the output value from mulch physical simulation
c         MulchNodehNew(..*..),       [PHY]   [cm]        nodal water potential; the output value from mulch physical simulation
c         MulchNodeThNew(..*..),      [PHY]   [m3/m3]     nodal water content; the output value from mulch physical simulation
c         MulchEleTmpr(..,..),        [PHY]   [oC]        element-based temperature, 1 soil surface, mulchLayer for mulch-atmosphere interface: for computation
c         MulchElehNew(..,..),        [PHY]   [cm]        element-based water potential, 1 soil surface, mulchLayer for mulch-atmosphere interface: for computation
c         MulchEleThNew(..,..),       [PHY]   [m3/m3]     element-based water content, 1 soil surface, mulchLayer for mulch-atmosphere interface: for computation
c         MulchElehNew_temp(..,..),   [PHY]   [cm]        element-based water potential for iteration, used in the calculating temporary mulch physical properties
c         MulchEleTmpr_temp(..,..),   [PHY]   [oC]        element-based temperature for iteration, used in the calculating temporary mulch physical properties
c         MulchEleThNew_temp(..,..),  [PHY]   [m3/m3]     element-based water content for iteration, used in the calculating temporary mulch physical properties
c         MulchElehNew_temp0(..,..),  [PHY]   [cm]        element-based water potential for iteration   --- store for local step data
c         MulchEleTmpr_temp0(..,..),  [PHY]   [oC]        element-based temperature for iteration       --- store for local step data
c         MulchElehNew_temp2(..,..),  [PHY]   [cm]        element-based water potential for iteration   --- first to receive any updates, after solving the linear system of governing equations
c         MulchEleTmpr_temp2(..,..),  [PHY]   [oC]        element-based temperature for iteration       --- first to receive any updates, after solving the linear system of governing equations
c         MulchElehNew_temp3(..,..),  [PHY]   [cm]        element-based water potential for iteration   --- static variable for multi-iterantion moving average
c         MulchEleTmpr_temp3(..,..),  [PHY]   [oC]        element-based temperature for iteration       --- static variable for multi-iterantion moving average
c         meanMulchHnew,              [PHY]   [cm]        store a temporary horizontal averaged mulch hnew
c         shortRadDir(NumMulHLD),     [PHY]   [W/m2]      shortwave downwards radiaiton at each interface of element-based mulch layers, i.e., interface between to mulch layer (mulchLayer+1 layers total) [w/m^2 * (cm/100) = W/m]
c         shortRadFirst(..,..),       [PHY]   [W/m2]      short wave upwards radiaiton (first order reflection) at each interface of element-based mulch layers
c         longRadDown(..,..),         [PHY]   [W/m2]      long wave downwards radiaiton at each interface of element-based mulch layers
c         longRadUp(11,50),           [PHY]   [W/m2]      long wave upwards radiaiton at each interface of element-based mulch layers
c         netRad(NumMulHLD,NumBPD),   [PHY]   [W/m2]      net radiation at each interface of element-based mulch layers, i.e., interface between to mulch layer (mulchLayer+1 layers total) [w/m^2 * (cm/100) = W/m]
c         netRadEachLayer(..,..),     [PHY]   [W/m2]      net radiation of each element-based mulch layer (mulchLayer layers total) [w/m^2 * (cm/100) = W/m]
c         longEmission,               [PHY]   [W/m2]      longwave radiation at each mulch layer based on the temperature
c         longEmissionAir,            [PHY]   [W/m2]      longwave radiation calculation from atmosphere based on the temperature
c         longEmissionSoil,           [PHY]   [W/m2]      longwave radiation calculation from soil surface based on the temperature
c         RelaHumid_mulch(..,..),     [PHY]   [%]         relative Humidity within each mulch element
c         VaporSat_mulch(..,..),      [PHY]   [g/m3]      saturated vapor density for each mulch element
c         VaporAct_mulch_D(..,..),    [PHY]   [g/m3]      actual vapor density for each mulch element
c         VaporAct_mulch_P(..,..),    [PHY]   [Pa]        actual vapor density for each mulch element: NOTE the unit is not kpa but PA
c         VaporDiff_mulch(..,..),     [PHY]   [g/m2/day]  vapor diffusivity for each mulch element
c         ParMulch(13,1),             [PHY]   [complicated]   water characteristic curve for solid portion in mulching (plant tissue)
c         WaterCapa_mulch(..,..),     [PHY]   [1/cm]      water capacity of each mulch element during iteration based on mulch (NOT based on solid), i.e., the solid part values * (1-f_mulch_pore)
c         HeatCapa_Air(..,..),        [PHY]   [J/m3/K]    heat capacity of mulch (air) of each mulch element during iteration
c         HeatDiff_mulch(..,..),      [PHY]   [J/m/K/day]    heat diffusion of mulch (air) of each mulch element during iteration
c         HeatDiff_fabric,            [PHY]   [J/m/K/s]      heat diffusion of mulch (fabric solid) of each mulch element during iteration
c         RadSolarShort_Wea,          [PHY]   [W/m2]      incoming short wave radiation updated FROM weather module
c         AirTemp_Wea,                [PHY]   [oC]        atmosphere tmperatureupdated            FROM weather module
c         AirVaporS_Wea,              [PHY]   [kPa]       atmosphere saturated vapor pressure 
c         AirVaporP_Wea,              [PHY]   [kPa]       atmosphere vapor pressure               FROM weather module
c         AirVPD_Wea,                 [PHY]   [kPa]       atmosphere vapor deficit                FROM weather module
c         CloudCoverFactor_Wea,       [PHY]   [0~1]       atmosphere cloud cover                  FROM weather module
c         AirWind_Wea,                [PHY]   [km/h]      atmosphere wind speed                   FROM weather module
c         RelaHumi_Wea,               [PHY]   [0~1]         atmosphere relative humidity            FROM weather module
c         VaporSat_ambient,           [PHY]   [g/m3]      atmosphere saturated vapor density
c         VaporAct_ambient,           [PHY]   [g/m3]      atmosphere actual vapor density
c         HeatCapa_ambient,           [PHY]   [J/m3/K]    atmosphere (air) heat capacity
c         Phi_m                       [PHY]   [...]       momentum correction factor (maybe not necessary or need to install later, discuss and remove it later wzj???)
c         Tmpr_Sur,                   [PHY]   [oC]        soil surface temperature, averaged of nearby two nodes for one 
c         hNew_Sur,                   [PHY]   [cm]        soil surface water potential, averaged of nearby two nodes for one 
c         VaporSat_Sur,               [PHY]   [g/m3]      soil surface saturated water vapor density, averaged of nearby two nodes for one 
c         RelaHumid_Sur,              [PHY]   [1]         soil surface saturated related humidity
c         VaporAct_Sur,               [PHY]   [g/m3]      soil surface actural vapor density: VaporSat_Sur*RelaHumid_Sur
c         VaporAct_Sur_p,             [PHY]   [PA]        soil surface actural vapor pressure
c         g_Heat_Sensi(..,..),        [PHY]   [J/m2/day]  sensible heat flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards fluxes
c         g_Heat_Latent(..,..),       [PHY]   [J/m2/day]  latent heat flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards fluxes
c         g_Vapor(..,..),             [PHY]   [g/m2/day]  water vapor flux, odd rows for vertical, even rows for horizontal, positive for downwards and leftwards fluxes
c         g_Heat_CD_total(1,..)       [PHY]   [J/cm2/day/K]   the surface heat conductance to assign new values for "heatmow", used for varbt adjusting
c         meanG_vapor,                [PHY]   [g/m2/day]  horizontal average of water vapor flux for element-based layers, temporary
c         meanG_latent_heat,          [PHY]   [J/m2/day]  horizontal average of latent heat flux for element-based layers, temporary
c         VaporCond_Mul_D,            [PHY]   [m/day]     the diffusive conductance of the mulch for vapor, use: VaporCond_Mul_D*(vapor density gradient g/m3)=vapor flux, g/m2/day
c         VaporCond_Mul_C,            [PHY]   [g/m2/day/Pa]   the sum of convective conductance of the mulch for vapor, use: VaporCond_Mul_C*(vapor pressure gradient Pa)=vapor flux, g/m2/day, not use Pa not KPA
c         VaporCond_Mul_C_Fr,         [PHY]   [g/m2/day/Pa]   the free conductance of the mulch for vapor, 
c         VaporCond_Mul_C_Fo,         [PHY]   [g/m2/day/Pa]   the forced conductance of the mulch for vapor, 
c         HeatCond_Mul_D,             [PHY]   [J/m2/day/K]    the diffusive conductance of the mulch for heat, use: HeatCond_Mul_D*(tempearture gradient K)=heat flux, J/m2/day
c         HeatCond_Mul_C,             [PHY]   [J/m2/day/K]    the sum of convective conductance of the mulch for heat, use: HeatCond_Mul_C*(tempearture gradient K)=heat flux, J/m2/day   
c         HeatCond_Mul_C_Fr,          [PHY]   [J/m2/day/K]    the free conductance of the mulch for heat,
c         HeatCond_Mul_C_Fo,          [PHY]   [J/m2/day/K]    the forced conductance of the mulch for heat,
c         DeltaRshort,                [PHY]   [1]         residue-area index of each mulch layer, Novak et al. (2000); (1-DeltaRshort) means the radiaiton at this layer.
c         DeltaRlong,                 [PHY]   [1]         defined by this program, "DeltaRshort" is for shortwave; "DeltaRlong" is for longwave
c         Omega,                      [PHY]   [1]         clumping index, which describes the way in which the residue elements are arranged, Novak et al., (2000); (1-Omega*DeltaRshort) is the way how it is used
c         alpha_m,                    [PHY]   [1]         short wave reflectivity of mulch
c         alpha_s,                    [PHY]   [1]         short wave reflectivity of soil surface (sometimes submerged in water)
c         epsilon_m,                  [PHY]   [1]         emissivity of mulch
c         epsilon_s,                  [PHY]   [1]         emissivity of soil
c         epsilon_a_0,                [PHY]   [1]         ambient air emissivity
c         coef_epsilon_a,             [PHY]   [complicated]       correction factors for ambient air emissivity
c         sigma,                      [PHY]   [W/m2/T4]           Stefan-Boltzmann constant
c         co_disp_heig,               [PHY]   [1]         coefficient of displacement height (surface displacement) for wind speed
c         co_rough_len,               [PHY]   [1]         coefficient of roughness length for wind speed
c         z_rough                     [PHY]   [cm]        surface roughness in adjusting windspeed, just think mulch is a 'thick' canopy
c         d_displace                  [PHY]   [cm]        surface displacement in adjusting windspeed, just think mulch is a 'thick' canopy
c         z_wind,                     [PHY]   [cm]        wind measuring height
c         rho_mulch,                  [PHY]   [g/m3]      solid part mulch density, assumed to be the density of wood or grass   (defined as public variable, see PuSurface.ins)
c         rho_mulch_b,                [PHY]   [g/m3]      mulch bulk density, rho_mulch*(1-f_mulch_pore)                         (defined as public variable, see PuSurface.ins)
c         rho_dryair,                 [PHY]   [g/m3]      density of dry air                                                     (defined as public variable, see PuSurface.ins)
c         rho_w,                      [PHY]   [g/m3]      water density                                                          (defined as public variable, see PuSurface.ins)
c         Tmpr_ref,                   [PHY]   [oC]        reference temperature for heat equation, especially the latent heat fluxes, can be set as 0 or 20 or arbitrary, refer to the zero-potential surface for soil water
c         c_vap_s,                    [PHY]   [J/g/K]     water vapor specific heat
c         c_vap_v,                    [PHY]   [J/g]       vaporization heat, i.e., phase change from water liquid to vapor
c         c_air_s,                    [PHY]   [J/g/K]     specific heat of air, may need to change to enthalpy
c         c_water_s,                  [PHY]   [J/g/K]     specific heat of water liquid
c         c_mulch_s,                  [PHY]   [J/g/K]     specific heat of mulch solid portion, we can use wood or grass specific heat. Thus, in calculation, should use mulch solid properties, such as "rho_mulch"
c         Dhm,                        [PHY]   [m2/s]      molecular diffusivity for sensible heat, "Dhm"=2.2D-5 m2/s
c         nu_air,                     [PHY]   [m2/s]      air kinematic viscosity, "nu_air"=1.5D-5 m2/s
c         Ra_Critical,                [PHY]   [1]         critical number for Rayleigh number
c         GrRe_Critical,              [PHY]   [1]         critical number for Richardson number, i.e., Grashof/Reynold^2
c         a_free,                     [PHY]   [m/K^1/2/s] an empirical parameter for free convection, see Finderling et al. (2003), page 6-4
c         M_H2O,                      [PHY]   [kg/mol]    the molecular weight of water, i.e., "M_H2O"=0.018 kg/mol
c         R_gas,                      [PHY]   [J/K/mol]   gas constant, i.e., "R_gas"=8.314 J/K/mol
c         karman,                     [PHY]   [1]         von Karman constant, i.e., 0.4
c         u_mulch(NumMulHLD),         [PHY]   [m/s]       wind speed at mulch interior interface, i.e., interface between two element-based mulch layers
c         u_above,                    [PHY]   [m/s]       wind speed at mulch-air interface, read as (km/h) but convert to (m/s) in the simulation
c         u_soilsur,                  [PHY]   [m/s]       wind speed at mulch-soil surface, read as (km/h) but convert to (m/s) in the simulation (Novak et al. literature seq)
c         h_Pond_max,                 [PHY]   [cm]        connection with surface water flow, the maximum ponded depth
c         h_Pond_max_thre,            [PHY]   [cm]        connection with surface water flow, the ponded depth threshold, "h_Pond_max.ge.h_Pond_max_thre" then consider ponded water cases, so "h_Pond_max_thre" should be relatively small
c         RainFallInput(..,..),       [PHY]   [g/m2/day]  rainfall input for each mulch layer
c         RainFallInput_temp(..,..),  [PHY]   [g/m2/day]  temperatory rainfall input for each mulch layer
c         inputPerLayer,              [PHY]   [g/m2/day]  iterative variable for rainfall for each layer
c         PriorStep,                  [GEO]   [day]       save the previous time step (maybe not necessary, discuss and remove it later wzj???)
c         mergeIndex(..),             [GEO]   [index]     input current layer index, output index after layer merging
c         DiffSubmerge,               [GEO]   [cm]        if a layer is partially filled, the "DiffSubmerge"=the height of that layer above free water surface. If DiffSubmerge=0, means totally filled
c         PerOccupation,              [GEO]   [1 or index]    -1: no ponded water; 0-1: fraction of "DiffSubmerge" in the current layer; 2: that layer is fully submerged; 3: mulch is totally submerged
c         TotalTime,                  [CON]   [day]       when local time step is used, "TotalTime" means the full time step, i.e., "TotalTime"="step" in the global space
c         LocalStep,                  [CON]   [day]       the subdivided or “refined” local time-step
c         time_mulch_old,             [CON]   [day]       the time for the previous iteration (of the whole program), the propose is to prevent backward time-shift in watermov, heatmov or solmov
c         minLocalStep,               [CON]   [day]       the minimum local time-step
c         Local_VarBW1,               [PHY]   [g/cm2/day]     local variable for varBW, water flux across the mulch-soil surface due to rainfall and irrigation, need average based on local time step
c         Local_VarBW2,               [PHY]   [g/cm2/day]     local variable for varBW, water flux across the mulch-soil surface due to evaporation, need average based on local time step
c         Local_VarBW3,               [PHY]   [g/cm2/day]     local variable for varBW, sum of Local_VarBW1, Local_VarBW3,need average based on local time step
c         Local_VarBT1,               [PHY]   [K]             local variable for varBT, need average based on local time step, mulch temperature @ mulch-soil interface
c         Local_VarBT2,               [PHY]   [J/cm2/day/K]   local variable for varBT, need average based on local time step
c         Local_VarBT3,               [PHY]   [J/cm2/day]     local variable for varBT, need average based on local time step, =Local_VarBT2*Local_VarBT1
c         Local_VarBT4,               [PHY]   [J/cm2/day]     local variable for varBT, need average based on local time step, radiation
c         VarBW1_temp(..),            [PHY]   [g/cm2/day]     local variable for varBW, weighted sum based on local time step
c         VarBW2_temp(..),            [PHY]   [g/cm2/day]     local variable for varBW, weighted sum based on local time step
c         VarBW3_temp(..),            [PHY]   [g/cm2/day]     local variable for varBW, weighted sum based on local time step
c         Q_temp(..),                 [PHY]   [g/day]         local variable for Q, weighted sum based on local time step
c         VarBT1_temp(..),            [PHY]   [K]             local variable for varBT, weighted sum based on local time step
c         VarBT2_temp(..),            [PHY]   [J/cm2/day/K]   local variable for varBT, weighted sum based on local time step
c         VarBT3_temp(..),            [PHY]   [J/cm2/day]     local variable for varBT, weighted sum based on local time step
c         VarBT4_temp(..),            [PHY]   [J/cm2/day]     local variable for varBT, weighted sum based on local time step
c         thMulchSat,                 [PHY]   [m3/m3]     saturated mulch water content, mulch based calculation, NOT solid based calculation
c         hMulchInit,                 [PHY]   [cm]        the initial value of mulch water potential, based on weather condition (temp+relative humid), mulch-based not soild-based
c         thMulchIni,                 [PHY]   [m3/m3]     the initial value of mulch water content, based on "hMulchInit", mulch-based not soild-based
c         waterCapaIni,               [PHY]   [1/cm]      the initial value of mulch water capacity, based on "hMulchInit", mulch-based not soild-based
c         CARB_mass(..,..),           [BIO]   [g]         mass of organic matter in Carbohydrates (CARB) for each element,
c         CELL_mass(..,..),           [BIO]   [g]         mass of organic matter in Hemi-celluloses/Holo-cellulose (CELL) for each element,
c         LIGN_mass(..,..),           [BIO]   [g]         mass of organic matter in Lignin (LIGN) for each element,
c         CARB_mass_temp(..,..),      [BIO]   [g]         temporary carbonhydrates mass pool for layered shrinking and merging element
c         CELL_mass_temp(..,..),      [BIO]   [g]         temporary holo-cellulose mass pool for layered shrinking and merging element
c         LIGN_mass_temp(..,..),      [BIO]   [g]         temporary lignin mass pool for layered shrinking and merging element
c         CARB_mass_temp(..,..),      [BIO]   [g]         temporary carbonhydrates mass pool for layered shrinking and merging element
c         CELL_mass_temp(..,..),      [BIO]   [g]         temporary holo-cellulose mass pool for layered shrinking and merging element
c         LIGN_mass_temp(..,..),      [BIO]   [g]         temporary lignin mass pool for layered shrinking and merging element
c         LIGN_mass_frac(...),        [BIO]   [0~1]       layer-wise lignin mass fraction, used to control the water retention curve
c         Frac_Decomp(..(layer)),     [BIO]   [0~1]       the decomposition fraction of each layer, 0 means no decomposition, 1 means totally gone
c
c     DUMMY VARIABLE (varable not important or as coefficients, without clear physical meaning)
c         e,                          temporery variable for element iteration
c         hhhh,                       dummy variable to receive temperary vairable for node vertical coordinate calculation
c         aaaa,                       dummy variable to receive unimportant values/mark locations
c         a11,a12,a21,a22,b11,b12     dummy variable used in linear equation systems   
c         BackwardIndex,              dummy variable for the last step of layer merging
c -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


       Subroutine Surface_Mulch_Process()
      
       include 'public.ins'
       include 'puplant.ins'
       include 'puweath.ins'
       include 'PuSurface.ins'
       
cccz control variables
c    some dummy varibles that hold temperary values
       double precision aaaa,hhhh,a11,a12,a21,a22,b11,b12,
     !    ThNewtemp,hNewtemp,Temptemp 
      integer ModNum
c    parameters for shrinking the grid
       double precision PriorStep,          ! record the previous time step
     !    thickPerLayer_temp(NumMulHLD)
       double precision ThreHnew_HardFix, ThreHnew_UnderWater
       integer mergeIndex(NumMulHLD)
cccz parameter for control (iteration and time)
       integer  DiffusionRes,ForceForwards,                       !cccz Force to have only diffusive heat/vapor flow; Force the Picard iteration to go forwards
     !  MaxIter,IterMulch,TimeShrink,                             !cccz parameters for converging tests
     !  LongWaveRadCtrl,BackwardIndex
       double precision 
     !  maxhNewDiff,maxTmprDiff,maxhNew,maxTmpr,
     !  errTol,TotalTime, LocalStep, minLocalStep, time_mulch_old
cccz geometrical parameter for mulch nodes/elements
       integer  e,qLeft,qRight,nNode,InElement
cccz temporary geometrical parameter (area and cumulated elevation calculation)
       double precision Dist
cccz physical properties within mulch
       double precision     
c    radiation group (physical properties within mulch)
     !   DeltaRshort,DeltaRlong,Omega,epsilon_m,alpha_m,alpha_s,
     !   epsilon_s,epsilon_a_0,coef_epsilon_a,sigma,
c    wind speed group (physical properties within mulch)
     !   co_disp_heig,co_rough_len,                               
     !   u_above,u_soilsur,u_mulch(NumMulHLD),
     !   Phi_m, z_rough, d_displace, z_wind
c    thermal property group (physical properties within mulch)
     !   c_vap_s,c_vap_v,c_air_s,c_water_s,c_mulch_s,
     !   Tmpr_ref,HeatDiff_fabric,
c    mulch air stability group (physical properties within mulch) 
     !   Dhm,nu_air,Ra_Critical,GrRe_Critical,
     !   a_free,M_H2O,R_gas,karman
cccz physical process within mulch
       double precision
c    radiation group (physical process within mulch)     
     !   netRad(NumMulHLD,NumBPD),netRadEachLayer(NumMulHLD,NumBPD),
     !   shortRadDir(NumMulHLD),shortRadFirst(NumMulHLD,NumBPD),
     !   longRadDown(NumMulHLD,NumBPD),longRadUp(NumMulHLD,NumBPD),     
     !   longEmission,longEmissionAir,longEmissionSoil,
c     heat/water flux group (physical process within mulch)    
     !   g_Heat_Sensi(NumMulFluxD,NumBPD),                            !cccz sensible heat flux
     !   g_Heat_Latent(NumMulFluxD,NumBPD),                           !cccz latent heat flux
     !   g_Vapor(NumMulFluxD,NumBPD),                                 !cccz vapor flux
     !   g_Heat_CD_total(1,NumBPD),                                   !cccz parameter on soil surface, only used to reset varbt2
     !   meanG_vapor, meanG_latent_heat,                              !cccz some horizontal means
c     flux conductance group
     !   VaporCond_Mul_D, VaporCond_Mul_C, 
     !   VaporCond_Mul_C_Fr, VaporCond_Mul_C_Fo,
     !   HeatCond_Mul_D, HeatCond_Mul_C,
     !   HeatCond_Mul_C_Fr, HeatCond_Mul_C_Fo
cccz physical status of mulch
       double precision 
     !   MulchEleTmpr_temp(NumMulHLD,NumBPD),                         ! cccz for exterior Picard iteration
     !   MulchElehNew_temp(NumMulHLD,NumBPD),                         ! cccz for exterior Picard iteration
     !   MulchEleThNew_temp(NumMulHLD,NumBPD),                        ! cccz for exterior Picard iteration
     !   MulchElehNew_temp2(NumMulHLD,NumBPD),                        ! cccz for interior Picard iteration
     !   MulchEleTmpr_temp2(NumMulHLD,NumBPD),                        ! cccz for interior Picard iteration
     !   MulchElehNew_temp3(NumMulHLD,NumBPD),                        ! cccz for cumulative average
     !   MulchEleTmpr_temp3(NumMulHLD,NumBPD),                        ! cccz for cumulative average
     !   MulchElehNew_temp0(NumMulHLD,NumBPD),                        ! cccz for time domain subdivision
     !   MulchEleTmpr_temp0(NumMulHLD,NumBPD),                        ! cccz for time domain subdivision
     !   meanMulchHnew
cccz parameters for weather/soil conditions
       double precision
     !   RainFallInput(NumMulHLD,NumBPD),
     !   RainFallInput_temp(NumMulHLD,NumBPD),inputPerLayer,  !cccz parameter for rainfall redistribution 
     !   VaporSat_ambient,VaporAct_ambient,HeatCapa_ambient,
     !   Tmpr_Sur(NumMulHLD),hNew_Sur(NumMulHLD),                   !cccz derived surface values 
     !   VaporSat_Sur(NumMulHLD),RelaHumid_Sur(NumMulHLD),
     !   VaporAct_Sur(NumMulHLD),VaporAct_Sur_P(NumMulHLD),
     !   thMulchSat,hMulchInit,thMulchIni,
     !   RadSolarShort_Wea,AirTemp_Wea,AirVaporP_Wea,AirVaporS_Wea,
     !   AirVPD_Wea,CloudCoverFactor_Wea,AirWind_Wea,RelaHumi_Wea   
cccz parameter that connect surface water and mulch
       double precision 
     !   h_Pond_max,h_Pond_max_thre, 
     !   DiffSubmerge,PerOccupation
       integer  SubmergeIndex, layer_need_fix
cccz data exchange between mulch module and other models, via varbw and varbt
       double precision
     !    Local_VarBW1, Local_VarBW2, Local_VarBW3,
     !    Local_VarBT1, Local_VarBT2, Local_VarBT3, Local_VarBT4
cccz temperary physical properties
       double precision 
     !    RelaHumid_mulch(NumMulHLD,NumBPD), 
     !    VaporSat_mulch(NumMulHLD,NumBPD),
     !    VaporDiff_mulch(NumMulHLD,NumBPD),
     !    VaporAct_mulch_D(NumMulHLD,NumBPD),
     !    VaporAct_mulch_P(NumMulHLD,NumBPD),
     !    HeatCapa_Air(NumMulHLD,NumBPD),  
     !    HeatDiff_mulch(NumMulHLD,NumBPD),
     !    WaterCapa_mulch(NumMulHLD,NumBPD),
     !    WaterDifici_mulch(NumMulHLD,NumBPD),
     !    Ra_mulch(NumMulHLD,NumBPD),GrRe_mulch(NumMulHLD,NumBPD),
     !    VarBW1_temp(NumBPD),VarBW2_temp(NumBPD),
     !    VarBW3_temp(NumBPD),Q_temp(NumBPD),
     !    VarBT1_temp(NumBPD),VarBT2_temp(NumBPD),
     !    VarBT3_temp(NumBPD),VarBT4_temp(NumBPD)
cccz the organic matter mass tracer
       double precision 
     !    CARB_mass_temp(NumMulHLD,NumBPD),
     !    CELL_mass_temp(NumMulHLD,NumBPD),
     !    LIGN_mass_temp(NumMulHLD,NumBPD),
     !    CARB_N_mass_temp(NumMulHLD,NumBPD),
     !    CELL_N_mass_temp(NumMulHLD,NumBPD),
     !    LIGN_N_mass_temp(NumMulHLD,NumBPD),  
     !    LIGN_mass_frac(NumMulHLD)
cccz just for paper writting
c       double precision Tentative_time_output,
c     !    Tentative_output
cccz merge horizontal elements
       integer intSurK
cccz static variables
       Common /SurfaceAdj_Static/ DiffusionRes,LongWaveRadCtrl,
     !    DeltaRshort,DeltaRlong,Omega,epsilon_m,alpha_m,alpha_s,
     !    epsilon_s,epsilon_a_0,coef_epsilon_a, sigma,
     !    co_disp_heig,co_rough_len,u_above,u_soilsur,u_mulch,
     !    c_vap_s,c_vap_v,c_air_s,c_water_s,c_mulch_s,Tmpr_ref,
     !    Dhm,nu_air,Ra_Critical,GrRe_Critical,
     !    netRad,netRadEachLayer,shortRadDir,
     !    shortRadFirst,longRadDown,longRadUp,
     !    g_Heat_Sensi,g_Heat_Latent,g_Vapor,g_Heat_CD_total,
     !    MulchEleTmpr_temp,MulchElehNew_temp,MulchEleThNew_temp,
     !    MulchElehNew_temp3,MulchEleTmpr_temp3,thMulchSat,
     !    MaxIter,TimeShrink,errTol,
     !    RainFallInput,inputPerLayer,
     !    VaporSat_ambient,VaporAct_ambient,HeatCapa_ambient,
     !    PriorStep,mergeIndex,
     !    RelaHumid_mulch,VaporSat_mulch,VaporDiff_mulch,
     !    VaporAct_mulch_D,VaporAct_mulch_P,
     !    HeatDiff_mulch,HeatDiff_fabric,WaterCapa_mulch,
     !    WaterDifici_mulch,
     !    RadSolarShort_Wea,AirTemp_Wea,AirVaporP_Wea,AirVaporS_Wea,
     !    CloudCoverFactor_Wea,AirWind_Wea,RelaHumi_Wea,
     !    a_free,M_H2O,R_gas,karman,
     !    time_mulch_old, minLocalStep,
     !    h_Pond_max,h_Pond_max_thre,LIGN_mass_frac,
     !    ThreHnew_HardFix,ThreHnew_UnderWater,ModNum
       
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------
cccz -----------------------------Initialization---------------------------------------------
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------

cccz reset the varbt/varbw value before each
cccz reset everything based on the upper layer, i.e., air
	  do i=1,NumBP
	   Varbt(i,1)=Varbt_Air(i,1)
	   Varbt(i,2)=Varbt_Air(i,2)
	   Varbt(i,3)=Varbt_Air(i,3)
	   Varbt(i,4)=Varbt_Air(i,4)
	   Varbt_mulch(i,1)=Varbt_Air(i,1)
	   Varbt_mulch(i,2)=Varbt_Air(i,2)
	   Varbt_mulch(i,3)=Varbt_Air(i,3)
	   Varbt_mulch(i,4)=Varbt_Air(i,4)
	   Varbw(i,1)=Varbw_Air(i,1)
	   Varbw(i,2)=Varbw_Air(i,2)
	   Varbw(i,3)=Varbw_Air(i,3)
	   Varbw_mulch(i,1)=Varbw_Air(i,1)
	   Varbw_mulch(i,2)=Varbw_Air(i,2)
	   Varbw_mulch(i,3)=Varbw_Air(i,3)
cccz add gas code here to reset the flux, before future adjustment via mulch and runoff
       do jjj=1,NumG
        Varbg_Mulch(i,jjj,1)=Varbg_Air(i,jjj,1)
        Varbg_Mulch(i,jjj,2)=Varbg_Air(i,jjj,2)
        Varbg_Mulch(i,jjj,3)=Varbg_Air(i,jjj,3)
        Varbg(i,jjj,1)=Varbg_Air(i,jjj,1)
        Varbg(i,jjj,2)=Varbg_Air(i,jjj,2)
         Varbg(i,jjj,3)=Varbg_Air(i,jjj,3)
       enddo

	  enddo
        g_vapor(:,:)=0.0
       if (residueApplied.le.0) then
         return
       endif
cccz initialization step 1: at the very beginning, read files, create time variables
       If (lInput.eq.1) then
           
          ThreHnew_HardFix=-40000.0D0
          ThreHnew_UnderWater=-40000.0D0
          BoolMulch_TotalDecomposed=0

cccz pre-initialize some variables (or constant)
cccz       mulchThick=10.0D0                      !cccz, static var in "PuSurface.ins"               ! move this to management file
       f_mulch_pore=0.8D0                     !cccz, static var in "PuSurface.ins"
       rho_mulch=800000.0D0                   !cccz, static var in "PuSurface.ins"
       rho_w=1000000.0D0                      !cccz, static var in "PuSurface.ins"
       rho_dryair=1225.0D0                    !cccz, static var in "PuSurface.ins"
       nu_air=1.5D-5                          !cccz, static var in "SurfaceMulchAdjustment.for"
       Dhm=2.2D-5                             !cccz, static var in "SurfaceMulchAdjustment.for"
       c_air_s=0.718D0                        !cccz, static var in "SurfaceMulchAdjustment.for"
       c_vap_s=1.862D0                        !cccz, static var in "SurfaceMulchAdjustment.for"
       c_vap_v=2453.50D0                      !cccz, static var in "SurfaceMulchAdjustment.for"
       c_water_s=4.186D0                      !cccz, static var in "SurfaceMulchAdjustment.for"
       c_mulch_s=1.76D0                       !cccz, static var in "SurfaceMulchAdjustment.for"
       Ra_Critical=1706.0D0                   !cccz, static var in "SurfaceMulchAdjustment.for"
       GrRe_Critical=1.0D0                    !cccz, static var in "SurfaceMulchAdjustment.for"
       a_free=5.0D-2                          !cccz, static var in "SurfaceMulchAdjustment.for"
       M_H2O=0.018D0                          !cccz, static var in "SurfaceMulchAdjustment.for"
       R_gas=8.314D0                          !cccz, static var in "SurfaceMulchAdjustment.for"
       Tmpr_ref=20.0D0                        !cccz, static var in "SurfaceMulchAdjustment.for"
cccz       mulchLayer=5                           !cccz, static var in "PuSurface.ins"                ! move this to management file
       DeltaRshort=0.28D0                     !cccz, static var in "SurfaceMulchAdjustment.for"
       DeltaRlong=0.28D0                      !cccz, static var in "SurfaceMulchAdjustment.for"
       Omega=0.60D0                           !cccz, static var in "SurfaceMulchAdjustment.for"
       sigma=5.67D-8                          !cccz, static var in "SurfaceMulchAdjustment.for"
       epsilon_s=1.0D0                        !cccz, static var in "SurfaceMulchAdjustment.for"
       epsilon_m=1.0D0                        !cccz, static var in "SurfaceMulchAdjustment.for"
       alpha_m=1.0D0                          !cccz, static var in "SurfaceMulchAdjustment.for"
       epsilon_a_0=0.60D0                     !cccz, static var in "SurfaceMulchAdjustment.for"
       coef_epsilon_a=5.95D-5                 !cccz, static var in "SurfaceMulchAdjustment.for"
       co_disp_heig=0.87D0                    !cccz, static var in "SurfaceMulchAdjustment.for"
       co_rough_len=0.079D0                   !cccz, static var in "SurfaceMulchAdjustment.for"
       karman=0.4D0                           !cccz, static var in "SurfaceMulchAdjustment.for"
       HeatDiff_fabric=0.25D0                 !cccz, static var in "SurfaceMulchAdjustment.for"

cccz hold the current time           
       time_mulch_old=time
       minLocalStep=0.01D0*dtMin
       minLocalStep=dtMin
       
cccz ponded water initialization (very arbitrary)
       h_Pond_max=0.0D0
       h_Pond_max_thre=0.05D0
                  
cccz first read the mulching information from file
       im=90
       il=0
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
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100) Min_Hori_Ele_Size
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100)
       im=im+1
       il=il+1
       Read(20,*,ERR=2100) DiffusionRes, LongWaveRadCtrl, DecompCtrl
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
      Read(20,*,ERR=2100) DeltaRshort,DeltaRlong,Omega,epsilon_m,alpha_m
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
       Read(20,*,ERR=2100) MaxIter, errTol
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
       Read(20,*,ERR=2100) rho_mulch, f_mulch_pore, WaterStorageFrac
       close(20)  
        
cccz this two values must be zero when lInput=1, i.e., file-reading step
       MulchDecompIni=0
       BoolMulchApply=0
       
       rho_mulch_b=rho_mulch*(1-f_mulch_pore)
cccz need to compute the correct input
       if(Thick_Mass.eq.'t') then
          ! cccz no need to change the input
       elseif (Thick_Mass.eq.'m') then
          mulchThick=0.1D0*mulchThick                              ! convert unit from "t/ha" to "kg/m^2"
          mulchThick=mulchThick/rho_mulch_b*1000.0D0*100.0D0       ! convert "kg/m^2" to "cm"  (cccz unit checked)
       endif
       
cccz just run the following subroutines once to initialize them       
        aaaa=WQ_CERES_MULCH(-5000.0D0,rho_mulch_b,0.1D0,lInput)    ! the CERES model is based on bulk mulch instead of solid, so there is no need for a more conversion based on "f_mulch_pore"
        aaaa=WH_CERES_MULCH(0.01D0,rho_mulch_b,0.1D0,lInput)
        waterCapaIni=WC_CERES_MULCH(-5000.0D0,rho_mulch_b,0.1D0,lInput)
        thMulchSat=WQ_CERES_MULCH(-0.000001D0,rho_mulch_b,0.1D0,0)
        call MulchDecomposition()
       Endif   !lInput =1
       
cccz initialization step 2: Initialize the computation when mulch is applied
c "begin_mulch_date" is the date of mulch app if there is mulch. until the mulch app time
c we don't run the code
       If(Time.ge.(begin_mulch_date-0.001*Step)    ! add a buffer 0.001*Step here
     &   .and.boolMulchApply.eq.0) then 
           
       BoolMulchApply=1    
       BoolMulch_TotalDecomposed=0
cccz  transfer weather data from the weather model to mulch model
       RadSolarShort_Wea=WATTSM(ITIME)                     !cccz shortwave solar radiation W/m^2, based on "HSR" in radiation module
       AirTemp_Wea=TAIR(ITIME)                             !cccz air temperature, oC
       AirVPD_Wea=VPD(ITIME)                               !cccz vapor deficity, kPa
       AirVaporS_Wea=0.61D0*EXP((17.27D0*AirTemp_Wea)      !cccz saturated vapor deficity, kPa, use Tetens Eq (Monteith and Unsworth, 2008)
     &   /(AirTemp_Wea + 237.3D0))
       AirVaporP_Wea=AirVaporS_Wea-AirVPD_Wea              !cccz actual vapor deficity, "saturated one - deficit"
       CloudCoverFactor_Wea=CLOUD                          !cccz 0-1 cloud cover factor
       AirWind_Wea=WIND                                    !cccz air wind speed, km/hour, will change to m/s later in the computation
       RelaHumi_Wea=AirVaporP_Wea/AirVaporS_Wea            !cccz 0-1 relative humidity

       VaporSat_ambient=exp(19.84D0-4975.9D0/(AirTemp_Wea+273.15D0))  !cccz Saturated Vapor density of ambient air (g/m^3)  Kimball et al., 1976
       VaporAct_ambient=VaporSat_ambient*RelaHumi_Wea                 !cccz Actual Vapor density of ambient air (g/m^3)
       HeatCapa_ambient=1.855D0*VaporAct_ambient                      !cccz Air and Vapor Heat Capacity of ambient air (J/m^3/K)
     &   +1.0035D0*max(rho_dryair-VaporAct_ambient,0.0D0)   
      
cccz ---------------------------------mulch grid configurations---------------------------------
cccz 'SurNodeIndex' = num of surface node, we keep it here because we want the surface computing system is relatively independent to other 2Dsoil module
       SurNodeIndex=0
       do i=1,NumBP
        n=KXB(i)
        k=CodeW(n)
        If(K.eq.4.or.K.eq.-4) then
         SurNodeIndex=SurNodeIndex+1
         SurfNodeNodeIndexH(SurNodeIndex)=n
        endif
       enddo                                              !cccz after this step, we know how many surface point we have. The assumption is the index start at a corner, supported by triangle
          
cccz SORT the surface node
cccz BASIC ASSUMP: the first node is always @ a corner, NOT in the middle of a surface. No additional assumption on the order of other node, or the shape of the surface.

      InElement=0             !cccz aux variables: detect if a variable is in a element  
      Dist=0.0D0              !cccz aux variables: choose the closest one to be the next surface node.
      do i=1,SurNodeIndex
       Dist=0.0D0
       do j=1,NumEL
        do kk=1,4
          if(SurfNodeNodeIndexH(i).eq.KX(j,kk)) InElement=1
        enddo
cccz choose the only one nearest point, next to the current step one
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
cccz Record the (x,y)-coordinate of surface node, Ordered by SurfNodeNodeIndexH
cccz we also find the peak elevation (for runoff under complicated surface configuration)
      PEAK=0.0D0
      do i=1,SurNodeIndex
        slopeCoord(i,1)=X(SurfNodeNodeIndexH(i))
        slopeCoord(i,2)=Y(SurfNodeNodeIndexH(i))
        if(PEAK.lt.slopeCoord(i,2)) then
           PEAK=slopeCoord(i,2)
        endif
      enddo
cccz Make an order based on surface node index, SurfNodeNodeIndexH/SurfNodeSurfIndexH are important geometric properties
      do n=1,SurNodeIndex
       do i=1,NumBP
         if(SurfNodeNodeIndexH(n).eq.KXB(i)) then
           SurfNodeSurfIndexH(n)=i
         endif
       enddo
      enddo
      
cccz mulch layer thick and hori distance for mulch nodes
      LayerHeight(1)=0.0D0
      totalMulchWidth=0.0D0
      do n=1,mulchLayer
        thickPerLayer(n)=mulchThick/dble(mulchLayer)            ! thickness per layer
        LayerHeight(n+1)=LayerHeight(n)+thickPerLayer(n)
        LIGN_mass_frac(n)=0.1D0
      enddo
      thresholdThick=max(0.1D0*thickPerLayer(1),0.05D0)              ! when a layer is too small to be a layer
      bbbb=0.0D0
      intSurK=1
      SurfMulchNodeSubIndex(intSurK)=1
      intSurK=intSurK+1
      SurNodeIndex_M=SurNodeIndex
      do n=1,SurNodeIndex
        SurNodeIndex_M_Match(n)=1
      enddo  
      do n=1,SurNodeIndex-1  
        bbbb=bbbb+abs(slopeCoord(n+1,1)-slopeCoord(n,1))
        if (bbbb.lt.Min_Hori_Ele_Size) then
            SurfMulchNodeSubIndex(intSurK)=n+1
            SurNodeIndex_M=SurNodeIndex_M-1
           SurNodeIndex_M_Match(intSurK)=SurNodeIndex_M_Match(intSurK)+1
        else
            SurfMulchNodeSubIndex(intSurK)=n+1
            intSurK=intSurK+1
cccz maybe need this line to make the widthPermulchUnit increases monotonically            
c            Min_Hori_Ele_Size=max(Min_Hori_Ele_Size,bbbb)
            bbbb=0.0D0
        endif
      enddo
c      SurNodeIndex_M=SurNodeIndex_M+1
      
      do n=1,SurNodeIndex_M
        qleft=SurfMulchNodeSubIndex(n)
        slopeCoord_M(n,1)=slopeCoord(qleft,1)
        slopeCoord_M(n,2)=slopeCoord(qleft,2)  
      enddo    
            
      do n=1,SurNodeIndex_M-1  
        widthPerMulchUnit(n)
     &      =abs(slopeCoord_M(n+1,1)-slopeCoord_M(n,1))
        totalMulchWidth=totalMulchWidth+widthPerMulchUnit(n)
      enddo

      
cccz Establish the mulching grid (Node)   
      numMulchNode=0
      hhhh=0.0D0
      do i=1,mulchLayer+1
       if (i.ge.2) then
        hhhh=hhhh+thickPerLayer(i-1)
       endif
       do n=1,SurNodeIndex_M
        numMulchNode=numMulchNode+1
        if(i.eq.1) then                                           !cccz Here we copy the soil surface node (established the interface)
         mulchNodeCoord(numMulchNode,1)=slopeCoord_M(n,1)
         mulchNodeCoord(numMulchNode,2)=slopeCoord_M(n,2)
        else                                                      !cccz Stack upper nodes above soil surface but under mulch-air interface, in an equal-distance way.
         mulchNodeCoord(numMulchNode,1)=slopeCoord_M(n,1)
         mulchNodeCoord(numMulchNode,2)=slopeCoord_M(n,2)+hhhh
        endif
        MulchNodeMatrix(i,n)=numMulchNode
        MulchNodeMarkerArea(numMulchNode)=0.0D0
c        MulchNodeTmpr(numMulchNode)=AirTemp_Wea
c        MulchNodehNew(numMulchNode)=hMulchInit                    !cccz the initial water potential
c        MulchNodeThNew(numMulchNode)=thMulchIni                   !cccz the initial water content
       enddo
      enddo        
cccz Establish the mulching grid (Element)
      numMulchEle=0
      do i=1,mulchLayer
       do n=1,SurNodeIndex_M-1
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
cccz triangle element is hard to be implement here, so far we decided not to implement triangular element
CDT TODO I guess we need to address this at some point 
        else
cccz rectangular element
         MulchEleMarkerArea(e)=abs((x2-x1)*(y3-y2))
         MulchNodeMarkerArea(n1)=MulchNodeMarkerArea(n1)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n2)=MulchNodeMarkerArea(n2)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n3)=MulchNodeMarkerArea(n3)
     &    +0.25D0*MulchEleMarkerArea(e)
         MulchNodeMarkerArea(n4)=MulchNodeMarkerArea(n4)
     &    +0.25D0*MulchEleMarkerArea(e)
        endif         ! test for trianglar element                
      enddo

cccz after initialize some geometrical properties, before determining the physical issues, need to initialize the mulch bio properties   
      MulchDecompIni=1
      call MulchDecomposition()
      MulchDecompIni=0
      do kk=1,mulchLayer
        aaaa=0.0D0
        hhhh=0.0D0
        do nn=1,SurNodeIndex_M-1
          aaaa=aaaa+LIGN_mass(kk,nn)
          hhhh=hhhh+CARB_mass(kk,nn)+CELL_mass(kk,nn)+LIGN_mass(kk,nn)
        enddo
        LIGN_mass_Frac(kk)=aaaa/hhhh
       enddo   

cccz Initial values for Node based mulch properties
cccz mulch hydrological potentail (negative, cm), mulch temperature (C), all the hydro and thermo properties
cccz Run the water-retention curve or thermo functions once so all the common parameter within those functions can be initialized. 


cccz Note here for dry mulch (the wrr paper) the initial mulch moisture is dependent on air moisture
c      hMulchInit=4708.34749D0*(AirTemp_Wea+273.15D0)*log(RelaHumi_Wea) ! 4708.34749 \approx 8.314D0/0.018D0/9.81D0*100.0D0 (cm, negative value), 100.0D0 change unit froom "m" to "cm"
cccz trial initialization (If we start with fresh mulch)   
      hMulchInit=-40000.0D0
cccz
      
      thMulchIni=WQ_CERES_MULCH(hMulchInit,rho_mulch_b,0.1D0,lInput)            ! the CERES model is based on bulk mulch instead of solid, so there is no need for a more conversion based on "f_mulch_pore"
        aaaa=WH_CERES_MULCH(thMulchIni,rho_mulch_b,0.1D0,lInput)
      waterCapaIni=WC_CERES_MULCH(-5000.0D0,rho_mulch_b,0.1D0,lInput)
      waterCapaIni=WC_CERES_MULCH(hMulchInit,rho_mulch_b,0.1D0,lInput)
        thMulchSat=WQ_CERES_MULCH(-0.000001D0,rho_mulch_b,0.1D0,0)

cccz initial mulch tmpr, hNew and ThNew

cccz node-wise
      do kk=1,numMulchNode
        MulchNodeTmpr(kk)=AirTemp_Wea
        MulchNodehNew(kk)=hMulchInit                    !cccz the initial water potential
        MulchNodeThNew(kk)=thMulchIni                   !cccz the initial water content
      enddo
      
cccz element-wise
      do e=1,numMulchEle
       n1=MulchEleMarker(e,1)
       n2=MulchEleMarker(e,2)
       n3=MulchEleMarker(e,3)
       n4=MulchEleMarker(e,4)
       kk=MulchEleMarker(e,6)! here we exchange the horizontal/vertical coordinate, in MulchEleMarker hori coord is 5; verti coord is 6; for physical matrix, verti comes first
       jj=MulchEleMarker(e,5)! here we exchange the horizontal/vertical coordinate, in MulchEleMarker hori coord is 5; verti coord is 6; for physical matrix, verti comes first
       MulchEleTmpr(kk,jj)=0.25D0     
     &  *(MulchNodeTmpr(n1)+MulchNodeTmpr(n2)
     &  +MulchNodeTmpr(n3)+MulchNodeTmpr(n4))
       MulchEleThNew(kk,jj)=0.25D0
     &  *(MulchNodeThNew(n1)+MulchNodeThNew(n2)
     &  +MulchNodeThNew(n3)+MulchNodeThNew(n4))
      enddo
      do kk=1,mulchLayer
       do jj=1,SurNodeIndex_M-1
       
        MulchElehNew(kk,jj)=
     &   WH_CERES_MULCH(MulchEleThNew(kk,jj),rho_mulch_b,
     &   LIGN_mass_Frac(kk),lInput)    

       enddo
      enddo

cccz ---------------------------------   mulch_runoff relation   ---------------------------------
cccz need to consider the current surface ponded water heights 'h_Pond' to see if the soil surface was submerged
cccz here we need to use "SurNodeIndex" rather than "SurNodeIndex_M" because it is based on "non-shrunk" grid based on surface soil.
       h_Pond_max=0.0D0
       do n=1,SurNodeIndex
          h_Pond_max=max(h_Pond_max,h_Pond(n))
       enddo

cccz determine if the ponded water reach which 'element-based layer': and mark the node of layer using 'SubmergeIndex', right under the water surface
       SubmergeIndex=0
       DiffSubmerge=0
       PerOccupation=0
       if(h_Pond_max.lt.h_Pond_max_thre) then
         SubmergeIndex=0
         DiffSubmerge=0
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
       if(SubmergeIndex.eq.0) then                !cccz no ponding water
         PerOccupation=-1.0D0
       elseif(SubmergeIndex.eq.mulchLayer+1) then !cccz all the mulch is submerged
         PerOccupation=3.0D0                                      
       else
        if(DiffSubmerge.le.thresholdThick) then   !cccz the current layer is submerged so deep, such that this layer is numerically considered as fully submerged
         PerOccupation=2.0D0                                      
        else                                      !cccz the current layer is submerged, but not that deep, such that a fraction of this layer is for mulch_adjust computation (non-submerged fraction)
         PerOccupation=DiffSubmerge/thickPerLayer(SubmergeIndex)  
        endif
       endif
       
cccz ---------------------------------initialize the radiation---------------------------------
          
cccz Calculate the shortwave downwards radiation in each interface among between adjacent mulch layers
cccz Based on the model in 'Simulating the radiation distribution within a barley-straw mulch' Novak 2000, use random distribution Eq(3) for transmittivity
        
       if(SubmergeIndex.le.0) then                                !cccz in the case there is no ponded water, radiation exists throughout the layer
        do kk=mulchLayer,1,-1
         if(kk.eq.mulchLayer) then                                !cccz elimiated the short direct raditaiton to the node under water 
          shortRadDir(kk+1)=RadSolarShort_Wea
         elseif(kk.eq.(mulchLayer-1)) then
          shortRadDir(kk+1)=RadSolarShort_Wea*(1.0D0-DeltaRshort)
         else
          shortRadDir(kk+1)=shortRadDir(kk+2)*(1.0D0-Omega*DeltaRshort)
         endif
        enddo
        shortRadDir(1)=shortRadDir(2)*(1.0D0-Omega*DeltaRshort)
       elseif(SubmergeIndex.lt.mulchLayer) then                   !cccz in the case there is no ponded water, radiation exists throughout the layer
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
       else                                                       !cccz the ponded water is higher than the mulch, mulch is fully submerged
         shortRadDir(mulchLayer+1)=RadSolarShort_Wea
         do kk=1,mulchLayer
            shortRadDir(kk)=0.0D0
         enddo
       endif ! end sumbmergIndex <=0
       
cccz ------------------------ initial output for paper writting-----------------
cccz initialize the mulch output options
c       Tentative_time_output=1.0D0/24.0D0/60.0D0
c       Tentative_output=time+Tentative_time_output
c       Open(9999,file='D:\MaizsimRunoff\BGR_B\Mulch_minute_paper.out')
c       Write(9999,9901) "Date_time","Rainfall",
c     &       "RunoffLeft","RunoffRight","Infiltration",
c     &       "Ra","Ri"
c       iday=int(time)
c       call caldat(iday,mm,id,iyyy) 
c       write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
c       do n=1, numMulchNode
c         Write(9999,9902) Time, Date, mulchNodeCoord(n,1),
c     &     mulchNodeCoord(n,2),MulchNodehNew(n),MulchNodeThNew(n),
c     &     MulchNodeTmpr(n)
c       enddo
cccz -----------------------------------------------------------------------
       Endif ! end of second initialization routine
cccz
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------
cccz -----------------------------Initialization (Two Steps) Finished------------------------
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------
cccz
      
cccz
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------  
cccz -----------------------------Calculation Started----------------------------------------
c    ----------------------------------------------------------------------------------------
c    ----------------------------------------------------------------------------------------

c??? one criterion for non-mulch, exit directly, should think about this condition
cccz also in the mulch is totally decomposed, then we set "BoolMulchApply=0" such that the following code will not need to be run
      if(mulchThick.lt.thresholdThick) then
cccz if this thing happened, we assign the mulch C and N to surface soil
c    just need to run once
c    the trigger is "BoolMulch_TotalDecomposed==0" && "mulchThick<thresholdThick"          
        if(BoolMulch_TotalDecomposed.eq.0) then
          BoolMulch_TotalDecomposed=1
          call MulchDecomposition()
          return
        endif  
      endif        
cccz
cccz time is not invoked     
      if(BoolMulchApply.eq.0) return
      if(BoolMulch_TotalDecomposed.eq.1) return
cccz hold the current time
c      if(time.le.time_mulch_old) then
c       aaaa=1      
c      endif
      time_mulch_old=time
      
cccz_try set the min time step adaptively based on the main time step
      minLocalStep=0.00001D0*step
cccz_try end
      
cccz  
cccz  ------------------------------weather data------------------------------
      RadSolarShort_Wea=WATTSM(ITIME)                     !cccz shortwave solar radiation W/m^2, based on "HSR" in radiation module
      AirTemp_Wea=TAIR(ITIME)                             !cccz air temperature, oC
      AirVPD_Wea=VPD(ITIME)                               !cccz vapor deficity, kPa
      AirVaporS_Wea=0.61D0*EXP((17.27D0*AirTemp_Wea)      !cccz saturated vapor deficity, kPa, use Tetens Eq (Monteith and Unsworth, 2008)
     &  /(AirTemp_Wea + 237.3D0))
      AirVaporP_Wea=AirVaporS_Wea-AirVPD_Wea              !cccz actual vapor deficity, "saturated one - deficit"
      CloudCoverFactor_Wea=CLOUD                          !cccz 0-1 cloud cover factor
      AirWind_Wea=WIND                                    !cccz air wind speed, km/hour, will change to m/s later in the computation
      RelaHumi_Wea=AirVaporP_Wea/AirVaporS_Wea            !cccz 0-1 relative humidity

cccz conver the weather data (change a unit or define another quantity)
      VaporSat_ambient=exp(19.84D0-4975.9D0/(AirTemp_Wea+273.15D0))   ! Saturated Vapor density of ambient air (g/m^3)                           
      VaporAct_ambient=VaporSat_ambient*RelaHumi_Wea                  ! Actual Vapor density of ambient air (g/m^3)
      HeatCapa_ambient=1.855D0*VaporAct_ambient
     &   +1.0035D0*(rho_dryair-VaporAct_ambient)                      ! Air and Vapor Heat Capacity of ambient air (J/m^3/K)
      
cccz reset the varbt/varbw value before each
cccz reset everything based on the upper layer, i.e., air
      do i=1,NumBP
       Varbt(i,1)=Varbt_Air(i,1)
       Varbt(i,2)=Varbt_Air(i,2)
       Varbt(i,3)=Varbt_Air(i,3)
       Varbt(i,4)=Varbt_Air(i,4)
       Varbt_mulch(i,1)=Varbt_Air(i,1)
       Varbt_mulch(i,2)=Varbt_Air(i,2)
       Varbt_mulch(i,3)=Varbt_Air(i,3)
       Varbt_mulch(i,4)=Varbt_Air(i,4)
       Varbw(i,1)=Varbw_Air(i,1)
       Varbw(i,2)=Varbw_Air(i,2)
       Varbw(i,3)=Varbw_Air(i,3)
       Varbw_mulch(i,1)=Varbw_Air(i,1)
       Varbw_mulch(i,2)=Varbw_Air(i,2)
       Varbw_mulch(i,3)=Varbw_Air(i,3)
      enddo

cccz ------------------------------grid shrinking--------------------------
c    adjust the mulch thickness, based on mulch-mass (desnity) conservation, i.e., assume the mulch-based density is constant
c    merge layers if mulch layer < critical value

cccz initialized the fraction factor for mulch shrinking
      do kk=1,mulchLayer
         Frac_Decomp(kk)=0.0D0
      enddo
      
      if(DecompCtrl.eq.0) goto 2200
      call MulchDecomposition()
c      call MulchDecomp()
cccz_try test if involvement of decom cause any data type probem
c     temperoarily stop changing grid
c      do kk=1,mulchLayer
c          Frac_Decomp(kk)=0.02D0-0.001D0*dble(kk)
c      enddo
cccz_try

      mulchLayer_temp=mulchLayer
      do kk=1,mulchLayer
        thickPerLayer(kk)=thickPerLayer(kk)*(1.0D0-Frac_Decomp(kk))
        mergeIndex(kk)=kk
      enddo
      kk=1
      kkk=1
      thickPerLayer_temp(kk)=thickPerLayer(kkk) 
2201  if(mulchLayer.gt.1.and.kkk.lt.mulchLayer) then     !cccz there are more than 1 layer, and there is at least one layer above to be merged, i.e., we follow upwards merging trend
       if(thickPerLayer_temp(kk).lt.thresholdThick) then
         thickPerLayer_temp(kk)=thickPerLayer_temp(kk)
     &     +thickPerLayer(kkk+1)
         mulchLayer_temp=mulchLayer_temp-1
         mergeIndex(kkk+1)=kk
         do jj=kkk+2,mulchLayer
          mergeIndex(jj)=mergeIndex(jj)-1
         enddo
         kkk=kkk+1
         goto 2201
       else
        kk=kk+1
        kkk=kkk+1
        thickPerLayer_temp(kk)=thickPerLayer(kkk) 
        goto 2201
       endif
      elseif(mulchLayer.gt.1.and.kkk.eq.mulchLayer) then
       if(thickPerLayer_temp(kk).lt.thresholdThick) then  
        thickPerLayer_temp(kk-1)=thickPerLayer_temp(kk)
     &     +thickPerLayer_temp(kk-1)
        mulchLayer_temp=mulchLayer_temp-1
        BackwardIndex=0
        do jj=mulchLayer,1,-1
          if(mergeIndex(mulchLayer).ne.mergeIndex(jj)
     &         .and.BackwardIndex.eq.0) then
             BackwardIndex=1 
             do jjj=jj+1, mulchLayer
               mergeIndex(jjj)=mergeIndex(jj)
             enddo
          endif
        enddo
        thickPerLayer_temp(kk)=0.0D0
       endif
      else          
      endif             
      PriorStep=Step

cccz ------------------------------Combine the physical parameters for new grid-----------------------------
c    initialize the temperatory physical coef: water head, water content, temperature
c    initialize the temperatory mulch coef: four-category of organic matter

      if(mulchLayer_temp.lt.mulchLayer) then
       do kk=1,mulchLayer
        do nn=1,SurNodeIndex_M-1
         MulchEleThNew_temp(kk,nn)=0.0D0
         MulchEleTmpr_temp(kk,nn)=0.0D0
         MulchElehNew_temp(kk,nn)=0.0D0
         CARB_mass_temp(kk,nn)=0.0D0
         CELL_mass_temp(kk,nn)=0.0D0
         LIGN_mass_temp(kk,nn)=0.0D0
         CARB_N_mass_temp(kk,nn)=0.0D0
         CELL_N_mass_temp(kk,nn)=0.0D0
         LIGN_N_mass_temp(kk,nn)=0.0D0
        enddo   
       enddo
       do kk=1,mulchLayer                                     !cccz cumulate everything that belongs to 'one new' element
        do nn=1,SurNodeIndex_M-1
         MulchEleThNew_temp(mergeIndex(kk),nn)=
     &    MulchEleThNew_temp(mergeIndex(kk),nn)
     &    +MulchEleThNew(kk,nn)*thickPerLayer(kk)
         MulchEleTmpr_temp(mergeIndex(kk),nn)=
     &    MulchEleTmpr_temp(mergeIndex(kk),nn)
     &    +MulchEleTmpr(kk,nn)*thickPerLayer(kk)
         CARB_mass_temp(mergeIndex(kk),nn)=
     &    CARB_mass_temp(mergeIndex(kk),nn)+CARB_mass(kk,nn)
         CELL_mass_temp(mergeIndex(kk),nn)=
     &    CELL_mass_temp(mergeIndex(kk),nn)+CELL_mass(kk,nn)
         LIGN_mass_temp(mergeIndex(kk),nn)=
     &    LIGN_mass_temp(mergeIndex(kk),nn)+LIGN_mass(kk,nn)
         CARB_N_mass_temp(mergeIndex(kk),nn)=
     &    CARB_N_mass_temp(mergeIndex(kk),nn)+CARB_N_mass(kk,nn)
         CELL_N_mass_temp(mergeIndex(kk),nn)=
     &    CELL_N_mass_temp(mergeIndex(kk),nn)+CELL_N_mass(kk,nn)
         LIGN_N_mass_temp(mergeIndex(kk),nn)=
     &    LIGN_N_mass_temp(mergeIndex(kk),nn)+LIGN_N_mass(kk,nn)
        enddo   
       enddo
       do kk=1,mulchLayer_temp                                !cccz finish normailizaiton & assignment
        thickPerLayer(kk)=thickPerLayer_temp(kk)   
        do nn=1,SurNodeIndex_M-1
         MulchEleThNew_temp(kk,nn)=
     &    MulchEleThNew_temp(kk,nn)/thickPerLayer_temp(kk)
         MulchEleTmpr_temp(kk,nn)=
     &    MulchEleTmpr_temp(kk,nn)/thickPerLayer_temp(kk)         
         MulchEleThNew(kk,nn)=MulchEleThNew_temp(kk,nn)       !cccz "_temp" to no "_temp" assignment, say good bye to old grid
         MulchEleTmpr(kk,nn)=MulchEleTmpr_temp(kk,nn)
         CARB_mass(kk,nn)=CARB_mass_temp(kk,nn)
         CELL_mass(kk,nn)=CELL_mass_temp(kk,nn)
         LIGN_mass(kk,nn)=LIGN_mass_temp(kk,nn)
         CARB_N_mass(kk,nn)=CARB_N_mass_temp(kk,nn)
         CELL_N_mass(kk,nn)=CELL_N_mass_temp(kk,nn)
         LIGN_N_mass(kk,nn)=LIGN_N_mass_temp(kk,nn)
cccz leave the water potential uncalculated until the mass_fraction of lignin is updated  
        enddo   
       enddo
       do kk=mulchLayer_temp+1,mulchLayer                     !cccz reset zeros for disappeared layers
        thickPerLayer(kk)=0.0D0   
        do nn=1,SurNodeIndex_M-1
          MulchEleThNew(kk,nn)=0.0D0
          MulchElehNew(kk,nn)=0.0D0
          MulchEleTmpr(kk,nn)=0.0D0
          CARB_mass(kk,nn)=0.0D0
          CELL_mass(kk,nn)=0.0D0
          LIGN_mass(kk,nn)=0.0D0
          CARB_N_mass(kk,nn)=0.0D0
          CELL_N_mass(kk,nn)=0.0D0
          LIGN_N_mass(kk,nn)=0.0D0
        enddo
       enddo
      endif ! end mulch layer comparison
      
cccz ------------------------------new mulch geometry------------------------------
c    New node, new element arrays
c    Area of each element
      numMulchNode_temp=0                                         !cccz Establish the mulching grid (Node)
      hhhh=0.0D0
      do i=1,mulchLayer_temp+1
       if(i.ge.2) then
        hhhh=hhhh+thickPerLayer(i-1)
       endif
       do n=1,SurNodeIndex_M
        numMulchNode_temp=numMulchNode_temp+1
        if(i.eq.1) then                                           !cccz Here we copy the soil surface node (establish the interface)
         mulchNodeCoord(numMulchNode_temp,1)=slopeCoord_M(n,1)
         mulchNodeCoord(numMulchNode_temp,2)=slopeCoord_M(n,2)
        else                                                      !cccz Stack upper nodes above soil surface, within mulch, not necessary "equal-distance" as for the initialization
         mulchNodeCoord(numMulchNode_temp,1)=slopeCoord_M(n,1)
         mulchNodeCoord(numMulchNode_temp,2)=slopeCoord_M(n,2)+hhhh
        endif
         MulchNodeMatrix(i,n)=numMulchNode_temp
         MulchNodeMarkerArea(numMulchNode_temp)=0.0D0
       enddo
      enddo
      do i=numMulchNode_temp+1,numMulchNode
         mulchNodeCoord(i,1)=0.0D0
         mulchNodeCoord(i,2)=0.0D0
         MulchNodeMarkerArea(i)=0.0D0
      enddo
      do i=mulchLayer_temp+2,mulchLayer+1
       do n=1,SurNodeIndex_M
         MulchNodeMatrix(i,n)=0
       enddo
      enddo
      numMulchEle_temp=0
      do i=1,mulchLayer_temp                                      !cccz Establish the mulching grid (Element)
       do n=1,SurNodeIndex_M-1
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
       do n=1,SurNodeIndex_M-1
         MulchEleMatrix(i,n)=0
       enddo
      enddo   
      do e=1,numMulchEle_temp                 !cccz Calculate the area and geometry of each node and element 
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
cccz triangle element is hard to be implement here, so far we decided not to implement triangular element
       else
cccz rectangular element
        MulchEleMarkerArea(e)=abs((x2-x1)*(y3-y2))
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
      do n=1,mulchLayer
        LayerHeight(n)=0.0D0
      enddo        
      mulchLayer=mulchLayer_temp
      numMulchEle=numMulchEle_temp
      numMulchNode=numMulchNode_temp
      mulchThick=0.0D0                !cccz finally record the new mulch height
      do kk=1,mulchLayer
        mulchThick=mulchThick+thickPerLayer(kk)
        LayerHeight(kk+1)=mulchThick
      enddo
      
cccz -------------------------- layer based mulch-lignin fraction --------------------------------
       do kk=1,mulchLayer                     !cccz reset zeros for disappeared layers
        aaaa=0.0D0
        hhhh=0.0D0
        do nn=1,SurNodeIndex_M-1
          aaaa=aaaa+LIGN_mass(kk,nn)
          hhhh=hhhh+CARB_mass(kk,nn)+CELL_mass(kk,nn)+LIGN_mass(kk,nn)
        enddo
        LIGN_mass_Frac(kk)=aaaa/hhhh
       enddo
       do kk=1,mulchLayer
         do nn=1,SurNodeIndex_M-1
          
           MulchElehNew(kk,nn)=
     &      WH_CERES_MULCH(MulchEleThNew(kk,nn),
     &      rho_mulch_b,LIGN_mass_Frac(kk),lInput)   

         enddo
       enddo
c??? one criterion for non-mulch, exit directly, should think about this condition
2200  if(mulchThick.lt.thresholdThick) then
c         BoolMulch_TotalDecomposed=1
cccz we do do the computation one more time and go out of mulch in the next time step.         
cccz     return
      endif      

cccz -------------------------- runoff-mulch interactions ----------------------------------------    

      h_Pond_max=0.0D0                !cccz need to consider the current surface ponded water heights 'h_Pond'
cccz note this "SurNodeIndex" CANNOT be "SurNodeIndex_M" becaues it counts the runoff quantity, not mulch gird
      do n=1,SurNodeIndex
        h_Pond_max=max(h_Pond_max,h_Pond(n))
      enddo
      SubmergeIndex=0                 !cccz in another word, node of layer 'SubmergeIndex' under the water surface; no ponded water, then "SubmergeIndex=0"
      if(h_Pond_max.lt.h_Pond_max_thre) then
        SubmergeIndex=0
        DiffSubmerge=0.0D0
      else
       do n=1,mulchLayer+1
        if(h_Pond_max.ge.LayerHeight(n)) then
          SubmergeIndex=n
          if(n.le.mulchLayer) then
            DiffSubmerge=LayerHeight(n+1)-h_Pond_max
          else
            DiffSubmerge=0.0D0  
          endif
        endif
       enddo
      endif
      if(SubmergeIndex.eq.0) then                !cccz no ponding water
        PerOccupation=-1.0D0
      elseif(SubmergeIndex.eq.mulchLayer+1) then !cccz all the mulch is submerged
        PerOccupation=3.0D0                                      
      else
       if(DiffSubmerge.le.thresholdThick) then   !cccz the current layer is submerged so deep, such that this layer is numerically considered as fully submerged
         PerOccupation=2.0D0                                      
       else                                      !cccz the current layer is submerged, but not that deep, such that a fraction of this layer is for mulch_adjust computation (non-submerged fraction)
         PerOccupation=DiffSubmerge/thickPerLayer(SubmergeIndex)  
       endif
      endif

cccz -------------------Radiation(short wave radiaiton)-------------------------------------
c    Sshortwave downwards radiation in each interfaces of mulch layers, 'Simulating the radiation distribution within a barley-straw mulch' Novak 2000, Eq(3) for transmittivity
c    unit for all radiation terms "w m^-2=J m^-2 s^-1"

      if(SubmergeIndex.le.0) then                 !cccz in the case there is no ponded water
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
      elseif(SubmergeIndex.le.mulchLayer) then    !cccz partially filled with water, elimiated the short direct raditaiton to the node under water 
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
      else                                        !cccz the whole mulch is under water 
       shortRadDir(mulchLayer+1)=RadSolarShort_Wea
       do kk=1,mulchLayer
        shortRadDir(kk)=0.0D0
       enddo
      endif
cccz ---------------------------------------------------------------------------------------
c    First order short wave reflective (upwards)
c    unit for all radiation terms "w m^-2=J m^-2 s^-1"
      if(SubmergeIndex.le.0) then                                     !cccz this is for "none-ponded" water case
       do i=1,(SurNodeIndex_M-1)
        qLeft=SurfMulchNodeSubIndex(i)
        qRight=SurfMulchNodeSubIndex(i+1)
        ThNewtemp=0.0D0
        do jjj=qLeft,qRight
            ThNewtemp=ThNewtemp+ThNew(SurfNodeNodeIndexH(jjj))
        enddo
        ThNewtemp=ThNewtemp/(abs(qLeft-qRight)+1)
        alpha_s=0.25D0-0.05D0*ThNewtemp                               ! surface albedo
c        qLeft=SurfNodeNodeIndexH(i)
c        qRight=SurfNodeNodeIndexH(i+1)
c        alpha_s=0.25D0-0.025D0*(ThNew(qLeft)+ThNew(qRight))           ! surface albedo
        shortRadFirst(1,i)=shortRadDir(1)*alpha_s                     ! reflection @ soil surface
        shortRadFirst(2,i)=shortRadDir(1)*alpha_s                     ! the first mulching interface
     &    *(1.0D0-DeltaRshort)+alpha_m*(shortRadDir(2)-shortRadDir(1))
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
      elseif(SubmergeIndex.le.mulchLayer) then                        !cccz this is revised when ponded water occurred
       do i=1,SurNodeIndex_M-1
c        qLeft=SurfNodeNodeIndexH(i)
c        qRight=SurfNodeNodeIndexH(i+1)
c        if(h_Pond(i).eq.0.0D0) then
c          alpha_np=0.25D0-0.05D0*ThNew(qLeft)                      ! soil surface albedo
c        else
c          alpha_np=0.06D0*f_mulch_pore+alpha_m*(1.0D0-f_mulch_pore)   ! surface albedo; weighted average of water and mulch albedo
c        endif
c        if(h_Pond(i+1).eq.0.0D0) then
c          alpha_nq=0.25D0-0.05D0*ThNew(qRight)                      ! soil surface albedo
c        else
c          alpha_nq=0.06D0*f_mulch_pore+alpha_m*(1.0D0-f_mulch_pore)   ! surface albedo; weighted average of water and mulch albedo
c        endif
c        alpha_s=0.5D0*(alpha_np+alpha_nq)
           
        qLeft=SurfMulchNodeSubIndex(i)
        qRight=SurfMulchNodeSubIndex(i+1)
        alpha_s=0.0D0
        do jjj=qLeft,qRight
         if(h_Pond(jjj).eq.0.0D0) then 
           alpha_s=alpha_s
     &       +0.25D0-0.05D0*ThNew(SurfNodeNodeIndexH(jjj))            ! soil surface albedo
         else
           alpha_s=alpha_s
     &       +0.06D0*f_mulch_pore+alpha_m*(1.0D0-f_mulch_pore)        ! surface albedo; weighted average of water and mulch albedo 
         endif  
        enddo
        alpha_s=alpha_s/(abs(qLeft-qRight)+1)   
        do kk=1,SubmergeIndex-1                                       ! under water part
          shortRadFirst(kk,i)=0.0D0
        enddo          
        shortRadFirst(SubmergeIndex,i)=                               ! reflection liquid surface
     &    shortRadDir(SubmergeIndex)*alpha_s                        
        shortRadFirst(SubmergeIndex+1,i)=                             ! the first mulching interface
     &    shortRadDir(SubmergeIndex)*alpha_s*(1.0D0-DeltaRshort)
     &    +alpha_m*(shortRadDir(SubmergeIndex+1)
     &    -shortRadDir(SubmergeIndex))
        do kk=SubmergeIndex+1,mulchLayer                              ! loop from second mulching interface
          shortRadFirst(kk+1,i)=shortRadDir(SubmergeIndex)*alpha_s    ! add reflections from all the layers beneth the LOI
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
      else                                                            !cccz the mulch is totally submerged
       alpha_s=0.06D0*f_mulch_pore+alpha_m*(1.0D0-f_mulch_pore)       ! combine the albedo for water and mulch, based on volume fraction
       do i=1,SurNodeIndex_M-1
         shortRadFirst(mulchLayer+1,i)=RadSolarShort_Wea*alpha_s      ! reflection liquid surface
       enddo
       do kk=1,mulchLayer
        do i=1,SurNodeIndex_M-1   
         shortRadFirst(kk,i)=0.0D0
        enddo
       enddo  
      endif

cccz -----------------------------------Wind Speed Calculation----------------------------------------------------
cccz The wind calculation along soil surface is only proceed once because the ambient condition within one iteration is stable
cccz 'An introduction to environmental biophysics' Campbell & Norman, first determine the wind speed near soil surface
      Phi_m=0.0D0                                                     ! ??? the momentum correction factor, need more information later
      z_rough=co_rough_len*mulchThick                                 ! just think mulch is a 'thick' canopy
      d_displace=co_disp_heig*mulchThick                              ! just think mulch is a 'thick' canopy
      z_wind=200.0D0                                                  ! ??? mulch wind is always measured @ 10~60 cm height; two meters above the surface is where the wind was measured
      u_ast=karman*AirWind_Wea/log((z_wind-d_displace)/z_rough)       ! calcualate the friction velocity (km/hour)
      u_above=u_ast/karman*log((mulchThick-d_displace)/z_rough)       ! wind speed in mulch-air surface (km/hour)
      u_soilsur=0.21D0*u_ast                                          ! wind speed in mulch-soil surface (km/hour, Novak et al. literature seq)       
      u_dist=0.0D0
cccz we do not make difference between ponded/unponded;
cccz if ponded, we just do not use the wind speed in deep layers
      do n=1,mulchLayer
       u_dist=u_dist+thickPerLayer(n)
       u_mulch(n)=0.21D0*u_ast*exp(2.2D0*
     &     (u_dist-0.5D0*thickPerLayer(n))/mulchThick)                 ! mean wind speed within mulch layers (km/hour, Novak et al. literature seq)
      enddo
      u_mulch(mulchLayer+1)=0.21D0*u_ast*exp(2.2D0*(u_dist/mulchThick))! wind speed in mulch-air interface
      u_soilsur=u_soilsur*u_above/u_mulch(mulchLayer+1)/3.6D0          ! shift units to m/s
      do n=1,mulchLayer+1
       u_mulch(n)=u_mulch(n)*u_above/u_mulch(mulchLayer+1)/3.6D0
      enddo
      u_above=u_above/3.6D0
      
cccz ---------------------------Establish temporary variables-------------------------------
      do n=1,SurNodeIndex_M-1
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
cccz this part is for soil, so not using "SurNodeIndex_M" but use "SurNodeIndex" 
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
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=0.0D0
       enddo
      enddo
cccz ----------------------------------------- Soil surface calculation ----------------------------------------------------------------      
      if(SubmergeIndex.le.0) then                                 ! this is the unsubmerged case refer to "goto 2300"
       do n=1,SurNodeIndex_M-1
         qLeft=SurfMulchNodeSubIndex(n)
         qRight=SurfMulchNodeSubIndex(n+1)
         hNewtemp=0.0D0
         Temptemp=0.0D0
         do jjj=qLeft,qRight
             hNewtemp=hNewtemp+hNew(SurfNodeNodeIndexH(jjj))
             Temptemp=Temptemp+Tmpr(SurfNodeNodeIndexH(jjj))
         enddo
         hNewtemp=hNewtemp/(abs(qLeft-qRight)+1)
         Temptemp=Temptemp/(abs(qLeft-qRight)+1)
c         qLeft=SurfNodeNodeIndexH(n)
c         qRight=SurfNodeNodeIndexH(n+1)
c         Tmpr_Sur(n)=0.5D0*(Tmpr(qLeft)+Tmpr(qRight))  
c         hNew_Sur(n)=0.5D0*(hNew(qLeft)+hNew(qRight))
         Tmpr_Sur(n)=Temptemp
         hNew_Sur(n)=hNewtemp  
         VaporSat_Sur(n)=exp(19.84D0-4975.9D0/(Tmpr_Sur(n)+273.15D0))
         RelaHumid_Sur(n)=max(min(exp(2.124D-4*hNew_Sur(n)
     &     /(Tmpr_Sur(n)+273.15D0)),1.0D0),0.0D0)    ! rh=exp(gMh/RT), 2.124D-4=(9.81m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
         VaporAct_Sur(n)=VaporSat_Sur(n)*RelaHumid_Sur(n)
         VaporAct_Sur_P(n)=
     &     0.61D0*exp((17.27D0*Tmpr_Sur(n))/(Tmpr_Sur(n)+237.3D0))
     &     *RelaHumid_Sur(n)*1000.0D0                ! unit in Pa
       enddo   
      elseif(SubmergeIndex.gt.0.and.PerOccupation.eq.3.0D0) then    ! this is the totally submerged case refer to "goto 2500"
       do n=1,SurNodeIndex_M-1
         Tmpr_Sur(n)=AirTemp_Wea 
         hNew_Sur(n)=0.0D0
         VaporSat_Sur(n)=exp(19.84D0-4975.9D0/(Tmpr_Sur(n)+273.15D0))
cccz!! Note change relative humidity from 1 to 0.999999 can make ~0.05 changes in soil surface temperature.
         RelaHumid_Sur(n)=1.0D0          
         VaporAct_Sur(n)=VaporSat_Sur(n)*RelaHumid_Sur(n)
         VaporAct_Sur_P(n)=
     &     0.61D0*exp((17.27D0*Tmpr_Sur(n))/(Tmpr_Sur(n)+237.3D0))
     &     *RelaHumid_Sur(n)*1000.0D0                                    ! unit in Pa
       enddo 
      elseif(SubmergeIndex.eq.mulchlayer
     &   .and.PerOccupation.eq.2.0D0) then    ! this is the partially submerged case, but treated as "totally submerge" refer to "goto 2500"
       do n=1,SurNodeIndex_M-1
         Tmpr_Sur(n)=AirTemp_Wea  
         hNew_Sur(n)=0.0D0
         VaporSat_Sur(n)=exp(19.84D0-4975.9D0/(Tmpr_Sur(n)+273.15D0))
         RelaHumid_Sur(n)=1.0D0
         VaporAct_Sur(n)=VaporSat_Sur(n)*RelaHumid_Sur(n)
         VaporAct_Sur_P(n)=
     &     0.61D0*exp((17.27D0*Tmpr_Sur(n))/(Tmpr_Sur(n)+237.3D0))
     &     *RelaHumid_Sur(n)*1000.0D0                                    ! unit in Pa
       enddo 
      else                                                            ! this is the partially submerged case refer to "goto 2400"
       do n=1,SurNodeIndex_M-1
         Tmpr_Sur(n)=AirTemp_Wea
         hNew_Sur(n)=0.0D0
         VaporSat_Sur(n)=exp(19.84D0-4975.9D0/(Tmpr_Sur(n)+273.15D0))
         RelaHumid_Sur(n)=1.0D0
         VaporAct_Sur(n)=VaporSat_Sur(n)*RelaHumid_Sur(n)
         VaporAct_Sur_P(n)=
     &     0.61D0*exp((17.27D0*Tmpr_Sur(n))/(Tmpr_Sur(n)+237.3D0))
     &     *RelaHumid_Sur(n)*1000.0D0                                    ! unit in Pa
       enddo    
      endif
cccz ---------------------------------------------- End soil surface calculation -------------------------------------------------------------------

      
2001  aaaa=1 !cccz iteration starting point
cccz --------------------rainfall redistribution initialization------------------------------------------------------
c    the submerged part will also be initialized but kept 0 (water contributed to VARBW) during the calculation
      do k=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(k)
       qRight=SurfMulchNodeSubIndex(k+1)
       bbbb=0.0D0
       do jjj=qLeft,qRight
         bbbb=bbbb+Varbw_Air(SurfNodeSurfIndexH(jjj),1)
       enddo   
       RainFallInput_temp(mulchLayer+1,k)=                            ! 'Varbw_Air' for water income from weather model
     &      bbbb*10000.0D0/(abs(qLeft-qRight)+1)                      ! recast rainfall unit from g/cm^2/day to g/m^2/day  (this is the basic unit for water/vapor fluxes)
c       kSur=SurfNodeSurfIndexH(k)
c       RainFallInput_temp(mulchLayer+1,k)=Varbw_Air(kSur,1)           ! 'Varbw_Air' for water income from weather model
c     &      *10000.0D0                                                ! recast rainfall unit from g/cm^2/day to g/m^2/day  (this is the basic unit for water/vapor fluxes)
       do n=1,mulchLayer
        RainFallInput_temp(n,k)=0.0D0
       enddo
      enddo
cccz -------------------Radiation(long wave radiaiton)-------------------------------------
c    long wave emission is related to the mulch temperature, so need to be within picard iteration
c    unit W/m^2*1cm or 1m (slab width)
cccz give a choice to remove longwave radiation (LongWaveRadCtrl=0-1)      
      if(LongWaveRadCtrl.eq.0) then
cccz downwards radiation
      epsilon_a=epsilon_a_0+coef_epsilon_a*(10.0D0*AirVaporP_Wea)     ! the air pressure used here is mb=0.1 Kpa, 10.0D0 change vapor pressure from kPa to mb; "AirVaporP_Wea" is of kPa
     &  *exp(1500.0D0/(AirTemp_Wea+273.15D0))
      epsilon_a=(1.0D0-0.84D0*CloudCoverFactor_Wea)*epsilon_a
     &  +0.84D0*CloudCoverFactor_Wea
      if(SubmergeIndex.le.0.0D0) then                                 !cccz no submerged part
       longEmissionAir=epsilon_a*sigma*((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex_M-1
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
      elseif(SubmergeIndex.le.mulchLayer) then                        !cccz when ponded water occurred
       longEmissionAir=epsilon_a*sigma*((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex_M-1
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
        do kk=1,SubmergeIndex-1         ! we do not consider the long-rad emission under water surface 
            longRadDown(kk,i)=0.0D0    
        enddo
       enddo
      else                                                            !cccz the ponded water is so high that all the mulch is submerged
       longEmissionAir=epsilon_a*sigma*((AirTemp_Wea+273.15D0)**4.0D0)
       do i=1,SurNodeIndex_M-1
        longRadDown(mulchLayer+1,i)=longEmissionAir
        do kk=1,mulchLayer
          longRadDown(kk,i)=0.0D0
        enddo
       enddo
      endif
cccz upwards radiation (W/m^2*1cm or 1m)
      if(SubmergeIndex.le.0.0D0) then                                 !cccz no submerged part
       do i=1,SurNodeIndex_M-1
c        qLeft=SurfNodeNodeIndexH(i)
c        qRight=SurfNodeNodeIndexH(i+1)
        longEmissionSoil=epsilon_s*sigma*
     &   ((Tmpr_Sur(n)+273.15D0)**4.0D0)
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
      elseif(SubmergeIndex.le.mulchLayer) then                        !cccz when ponded water occurred
       do i=1,SurNodeIndex_M-1
        do kk=1,SubmergeIndex-1
           longRadUp(kk,i)=0.0D0                                      ! we do not consider the long-rad emission under (ponded) water surface 
        enddo   
c        qLeft=SurfNodeNodeIndexH(i)
c        qRight=SurfNodeNodeIndexH(i+1)
        longEmissionSoil=epsilon_s*sigma*                             !cccz also use 'epsilon_s=1.0D0' for water cause water is considered as a blackbody
     &   ((Tmpr_Sur(n)+273.15D0)**4.0D0)    !c??? assume that the surface ponded water has the same temp as surface soil,
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
      else                                                            !cccz the ponded water is so high that all the mulch is submerged
       do i=1,SurNodeIndex_M-1
c        qLeft=SurfNodeNodeIndexH(i)
c        qRight=SurfNodeNodeIndexH(i+1)
        longEmissionSoil=epsilon_s*sigma*                             !cccz also use 'epsilon_s=1.0D0' for water cause water is considered as a blackbody
     &    ((Tmpr_Sur(n)+273.15D0)**4.0D0)   !c??? assume that the surface ponded water has the same temp as surface soil,
        longRadUp(mulchLayer+1,i)=longEmissionSoil
        do kk=1,mulchLayer
          longRadUp(kk,i)=0.0D0
        enddo
       enddo
      endif
cccz give a choice to remove longwave radiation     
      elseif(LongWaveRadCtrl.eq.1) then
cccz_try take out the longwave radiation
c    this should be done eventually by a switch (boolean var)
      do i=1,SurNodeIndex_M-1
       do j=1,mulchLayer+1
        longRadDown(j,i)=0.0D0
        longRadUp(j,i)=0.0D0
       enddo
      enddo
cccz_try finish this trial
      else
      endif
      
cccz --------- combine long/short wave radiation ----------------------------------
c    For each mulch layer interface: net radiation pass accross that interface
c    For each mulch layer: net radiation stay within that layer (HEAT SOURCE TERM)
    
      do i=1,SurNodeIndex_M-1
       do j=1,mulchLayer+1
        netRad(j,i)=(shortRadDir(j)-shortRadFirst(j,i))
     &    +(longRadDown(j,i)-longRadUp(j,i))
       enddo
      enddo
      do i=1,SurNodeIndex_M-1
       do j=1,mulchLayer
        netRadEachLayer(j,i)=netRad(j+1,i)-netRad(j,i)
       enddo
      enddo
      
cccz -------------------------------vapor and heat flux computation--------------------------------------
c    process the computaiton of fluxes and conversation laws
c    separate the three casess: no ponded water; partially filled; mulch submerged under ponded water
      if(SubmergeIndex.le.0) then                                 ! this is the unsubmerged case
          goto 2300
      elseif(SubmergeIndex.gt.0.and.PerOccupation.eq.3.0D0) then    ! this is the totally submerged case
          goto 2500
      elseif(SubmergeIndex.eq.mulchlayer
     &   .and.PerOccupation.eq.2.0D0) then    ! this is the partially submerged case, but treated as "totally submerge"
          goto 2500
      else                                                            ! this is the partially submerged case
          goto 2400
      endif
cccz -------------------------------------------------------------------------------------------
cccz --------------------------------NO PONDED WATER (23XX Labels)-----------------------------------------------------------
cccz -------------------------------------------------------------------------------------------
2300  do n=1,SurNodeIndex_M-1
       do k=1,mulchLayer 
cccz 'vapor density'(g/m^3)/'vapor partial pressure'(Pa)
        RelaHumid_mulch(k,n)=max(min(exp(2.124D-4*MulchElehNew_temp(k,n)! rh=exp(gMh/RT), 2.124D-4=(9.81m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0)),1.0D0),0.0D0)
        VaporSat_mulch(k,n)=exp(19.84D0-4975.9D0                        ! Saturated Vapor density for each mulch element (g/m^3)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))   
        VaporAct_mulch_D(k,n)=VaporSat_mulch(k,n)                       ! Actual Vapor density for each mulch element (g/m^3)
     &   *RelaHumid_mulch(k,n)
        VaporAct_mulch_P(k,n)=
     &   (0.61D0*exp((17.27D0*MulchEleTmpr_temp(k,n))
     &   /(MulchEleTmpr_temp(k,n)+237.3D0)))                            ! present this in Pa, not kPa, 1000.0: change unit from kPa to Pa
     &   *RelaHumid_mulch(k,n)*1000.0D0

cccz water difficit in each element: including liquid water and water vapor to a certain degree of saturation
c    vapor part: reach to rh=1.0D0
c    solid/liquid part: no need to change, they will make balance with the vapor parts.
        WaterDifici_mulch(k,n)=
     &   (VaporSat_mulch(k,n)-VaporAct_mulch_D(k,n))
     &   *(thickPerLayer(k)/100.0D0)/LocalStep                   ! rate to saturate the vapor portion (g/m^2/day)
     &   *f_mulch_pore                                           ! cast the rate to gas phase of the mulch    
cccz water diffusivity calculation: intrinsic diffusivity * temperature effects (Kimball et al., 1976; Lai et al., 1976)
         VaporDiff_mulch(k,n)=2.35D-5*                        ! 2.35=(2.29_soil science+2.40_[campbell book])/2
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**1.667D0)*86400.0D0                       ! Water Vapor diffusivity (m^2/day): only temperature and mulch dependent
cccz water capacity: only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore
c    similar to the water capacity calculation for soil           
         WaterCapa_mulch(k,n)=
     &     WC_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   

cccz heat capacity and diffusivity calculation: air density (rho_dryair) is in g m^-3, heat capacity 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air(k,n)=1.85D0*VaporAct_mulch_D(k,n)
     &   +1.006D0*rho_dryair                                 ! heat capacity (J/m^3/K), weighted average of specific heat
        HeatDiff_mulch(k,n)=1.84896D0*HeatCapa_Air(k,n)! heat diffusivity (J/m/K/day), maybe more like a conductance (just need temp gradient), 1.84896D0=2.14D-5(heat diff)*86400.0D0 for the air part
     &   +HeatDiff_fabric*86400.0D0*(1-f_mulch_pore)             ! this record heat path through the solid/fabric tissue     
        
       enddo
cccz determine the overall Ra and Ri for the whole layer
         Ra_mulch(1,n)=9.81D0                                 ! Rayleigh num
     &    *abs(0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &     -0.5D0*(MulchEleTmpr_temp(1,n)+Tmpr_Sur(n)))
     &    *((mulchThick/100.0D0)**3.0D0)
     &    /nu_air
     &    /((546.30D0
     &    +0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &    +0.5D0*(MulchEleTmpr_temp(1,n)+Tmpr_Sur(n)))/2.0D0)
     &    /Dhm
         GrRe_mulch(1,n)=2.0D0*9.81D0
     &    *(mulchThick/100.0D0)                               ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
     &    *abs(0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &     -0.5D0*(MulchEleTmpr_temp(1,n)+Tmpr_Sur(n)))
     &    /(546.30D0
     &     +0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &     +0.5D0*(MulchEleTmpr_temp(1,n)+Tmpr_Sur(n)))
     &    /((u_mulch(mulchLayer)-u_mulch(1))**2.0D0)  
      enddo

cccz ------------------ water (vapor) fluxes within mulch----------------------------
c    only vapor flux occur, but allow some liquid water stored in the mulch solid materials.
c    downward and leftward fluxes are positive
c    unit of water vapor flux: g/m^2/day
      do k=1,mulchLayer+1    
       if(k.eq.1) then            !cccz for mulch-soil interface                                        
        do n=1,SurNodeIndex_M-1         !cccz the soil surface vapor flux based on the min('potential evaporation' method, 'convection/diffusion' method), i.e., how much soil can supply, how much mulch can transport.
         ! METHOD1: 'potential evaporation' method
         SVPA_Sur=0.61D0*exp((17.27D0*MulchEleTmpr_temp(1,n))
     &    /(MulchEleTmpr_temp(1,n)+237.3D0))
         DEL_Sur=(0.61D0*exp((17.27D0*(MulchEleTmpr_temp(1,n)+1.0D0))
     &    /(MulchEleTmpr_temp(1,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch(1,n))      ! calculate VPD for the first mulching layer
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))               ! we use actural vapor pressure for D31 here, should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(1,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0             ! since we are doing local analysis within th mulch, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur*max(netRad(1,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*MulchEleTmpr_temp(1,n))))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_soilsur))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                    ! vertical vapor flux near soil-mulch interface (g/m2/day), 24.0 changes unit from g/m^2/hour to g/m^2/day
         ! METHOD2: 'convection/diffusion' method
c         Ra_mulch(1,n)=9.81D0*abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n))  ! Rayleigh num
c     &    *((0.5D0*thickPerLayer(1)/100.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur(n))/2.0D0)/Dhm
c         GrRe_mulch(1,n)=2.0D0*9.81D0*(0.5D0*thickPerLayer(1)/100.0D0) ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n))
c     &    /(546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur(n))
c     &    /((u_mulch(1)-u_soilsur)**2.0D0)
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then    ! diffusive vapor flow always happens
           VaporCond_Mul_D=VaporDiff_mulch(1,n)
     &      /(thickPerLayer(1)/200.0D0) 
           VaporCond_Mul_C=0.0D0
         else
           VaporCond_Mul_D=VaporDiff_mulch(1,n)
     &      /(thickPerLayer(1)/200.0D0)   
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then                   ! the Ri is large, means Re^2 is small, thus turbulence is small, hence free convection is processed, forced convection is neglated
            VaporCond_Mul_C_Fr=a_free
     &        *sqrt(abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur(n))/2.0D0)
     &        *1000.0D0*86400D0! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
            VaporCond_Mul_C_Fo=0.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else                                                        ! both free and forced convection should be considered
           VaporCond_Mul_C_Fr=a_free
     &        *sqrt(abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur(n))/2.0D0)
     &        *1000.0D0*86400D0
           VaporCond_Mul_C_Fo=
     &        (karman**2.0D0)*u_soilsur
c     &        /(log((LayerHeight(1)+0.5D0*thickPerLayer(1))
c     &        /(0.079D0*(LayerHeight(1)+0.5D0*thickPerLayer(1))))
c     &        **2.0D0)
     &        *0.155D0      
     &        *M_H2O/R_gas
     &        /((546.30D0+MulchEleTmpr_temp(1,n)+Tmpr_Sur(n))/2.0D0)
     &        *1000.0D0*86400D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          endif
         endif
         g_vapor_try2=VaporCond_Mul_D*
     &    (VaporAct_mulch_D(1,n)-VaporAct_Sur(n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(1,n)-VaporAct_Sur_P(n))
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then      !cccz take the min between 'diffusive type' and 'potential' evporation, if the two values are of different directions, we assume zero flux fo numerical stable     
           g_vapor(1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(1,n)=0.0D0
         endif               
        enddo      
        g_Vapor(2,1)=0.0D0                                            ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Vapor(2,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &     /(widthPerMulchUnit(n-1)/VaporDiff_mulch(1,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch(1,n))
     &     *(VaporAct_mulch_D(1,n)-VaporAct_mulch_D(1,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
         if(g_Vapor(2,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2,n)=min(g_Vapor(2,n),WaterDifici_mulch(1,n-1))
         else
           g_Vapor(2,n)=max(g_Vapor(2,n),-WaterDifici_mulch(1,n))  
         endif
        enddo  
      elseif(k.eq.(mulchLayer+1)) then                               !cccz for mulch-air interface         
        do n=1,SurNodeIndex_M-1                                           ! use diffusion/convection to calculate conductivity at atmosphere surface    
c         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)   ! Rayleigh num
c     &    *(((thickPerLayer(k-1))/100.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
c        GrRe_mulch(k,n)=2.0D0*9.81D0*((z_wind-mulchThick)/100.0D0)       ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
c     &    /(546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)
c     &    /((AirWind_Wea-u_mulch(k))**2.0D0)                  
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0) 
           VaporCond_Mul_C=0.0D0
         else
           VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air 
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)  
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
            VaporCond_Mul_C_Fr=a_free                                 ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0            
            VaporCond_Mul_C_Fo=0.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else
           VaporCond_Mul_C_Fr=a_free
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0
           VaporCond_Mul_C_Fo=
     &        (karman**2.0D0)*u_mulch(k)
c     &        /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &        *0.155D0
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D*
     &    *(VaporAct_Ambient-VaporAct_mulch_D(k-1,n))
         AirVaporP_temp=AirVaporP_Wea*1000.0D0
         if(AirVaporP_temp.le.VaporAct_mulch_P(k-1,n)) then  ! upwards flow
          g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &     VaporCond_Mul_C*
     &     (AirVaporP_temp-VaporAct_mulch_P(k-1,n))
          if(localstep.le.minLocalStep) then
            g_Vapor(2*k-1,n)=max(g_Vapor(2*k-1,n),
     &       min(g_Vapor(2*k-3,n),0.0D0)-WaterDifici_mulch(k-1,n))  
          endif  
         else
          g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &     min(VaporCond_Mul_C*
     &     (AirVaporP_temp-VaporAct_mulch_P(k-1,n)),
     &     WaterDifici_mulch(k-1,n))    
         endif
        enddo  
       else                                                   !cccz for water flux within mulch                                        
        do n=1,SurNodeIndex_M-1
c         Ra_mulch(k,n)=9.81D0                                 ! Rayleigh num
c     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
c     &    *(((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)**3.0D0)
c     &    /nu_air/((546.30D0+MulchEleTmpr_temp(k-1,n)
c     &    +MulchEleTmpr_temp(k,n))/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0
c     &    *((thickperlayer(k)+thickperlayer(k-1))/200.0D0)    ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n))
c     &    /(546.30D0+MulchEleTmpr_temp(k,n)+MulchEleTmpr_temp(k-1,n))
c     &    /((u_mulch(k)-u_mulch(k-1))**2.0D0)  
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
           VaporCond_Mul_C=0.0D0
         else 
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)  
           VaporCond_Mul_C_Fo=0.0D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results) 
           VaporCond_Mul_C_Fo=
     &       (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)   
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &        *(f_mulch_pore**1.667D0)        ! cccz: something compared with the diffusivity, adjust the accessible for path, very arbitrary, like a guess ??????????????????
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D
     &    *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k-1,n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(k,n)-VaporAct_mulch_P(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch(k,n))
     &     *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
          if(g_Vapor(2*k,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2*k,n)=min(g_Vapor(2*k,n),WaterDifici_mulch(k,n-1))
          else
           g_Vapor(2*k,n)=max(g_Vapor(2*k,n),-WaterDifici_mulch(k,n))  
          endif
        enddo
       endif
      enddo
           
cccz ------------------------Rainfall Redistribution--------------------------------------
c    unit "g/m^2/day"
      do n=1,SurNodeIndex_M-1
       k=mulchLayer+1
2308   if(k.gt.1) then
        inputPerLayer=
     &   min(RainFallInput_temp(k,n),WaterDifici_mulch(k-1,n))
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
     
cccz -----------------------Latent heat flux----------------------------------------------
c    downward and leftward are positive, energy income
c    calculate the fluxes based on the upwind direction of water
c    unit "J/m^2/day"
      do k=1,mulchLayer+1                                             
       if(k.eq.1) then                                                !cccz for soil surface layer
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(1,n).ge.0.0D0) then                               ! downwards vapor flow
           g_Heat_Latent(1,n)=c_vap_s*g_Vapor(1,n)
     &      *(MulchEleTmpr_temp(1,n)-Tmpr_ref)+c_vap_v*g_Vapor(1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(1,n)=c_vap_s*g_Vapor(1,n)*(Tmpr_Sur(n)-Tmpr_ref)
     &      +c_vap_v*g_Vapor(1,n)
         endif
        enddo
        g_Heat_Latent(2,1)=0.0D0
        g_Heat_Latent(2,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         if(g_Vapor(2,n).ge.0.0D0) then                               ! leftwards vapor flow
           g_Heat_Latent(2,n)=c_vap_s*g_Vapor(2,n)
     &      *(MulchEleTmpr_temp(1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2,n)
         else                                                         ! rightwards vapor flow
           g_Heat_Latent(2,n)=c_vap_s*g_Vapor(2,n)
     &      *(MulchEleTmpr_temp(1,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2,n)
         endif
        enddo
       elseif(k.eq.(mulchLayer+1)) then                               !cccz for mulch-air interface layer
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
       else                                                           !cccz for interier layers
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
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
      
cccz ---------------Diffusion/Convective Heat Flux(Sensible Heat)----------------------
c    downward and leftward are positive, i.e., energy income
c    unit (J/m^2/day)
      do k=1,mulchLayer+1
       if(k.eq.1) then                                                 !cccz for mulch-soil interface                                                  
        do n=1,SurNodeIndex_M-1
         qLeft=SurfMulchNodeSubIndex(k)
         qRight=SurfMulchNodeSubIndex(k+1)
         bbbb=0.0D0
         do jjj=qLeft,qRight
           bbbb=bbbb+varbt_Air(SurfNodeSurfIndexH(jjj),2)
         enddo   
         varbt_Temp=bbbb/(abs(qLeft-qRight)+1)
         varbt_Temp=max(varbt_Temp,0.0D0)
c         qLeft=SurfNodeSurfIndexH(n)
c         qRight=SurfNodeSurfIndexH(n+1)
c         varbt_Temp=0.5D0*(varbt_Air(qLeft,2)+varbt_Air(qRight,2))
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=HeatDiff_mulch(1,n)
     &     /(thickPerLayer(1)/200.0D0)                                ! the heat conductance (J/m^2/day/K), we will assume the "3rd dimension" has scale 1 cm
          HeatCond_Mul_C=0.0D0
         else
          HeatCond_Mul_D=HeatDiff_mulch(1,n)
     &     /(thickPerLayer(1)/200.0D0)                                ! the heat conductance (J/m^2/day/K), we will assume the "3rd dimension" has scale 1 cm
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(1,n)
     &       *86400D0                                                 ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(1,n)
     &       *86400D0                                                 ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_soilsur
c     &       /(log((LayerHeight(1)+0.5D0*thickPerLayer(1))
c     &       /(0.079D0*(LayerHeight(1)+0.5D0*thickPerLayer(1))))
c     &       **2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(1,n)
     &       *86400D0                                                 ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
         endif
        g_Heat_Sensi(1,n)=
     &   100.0D0*sqrt(varbt_Temp*(HeatCond_Mul_D+HeatCond_Mul_C))
c     &   (2.0D0/(1.0D0/(HeatCond_Mul_D+HeatCond_Mul_C)+1.0D0/Soil_Cond)) ! do not need consider the 
     &   *(MulchEleTmpr_temp(1,n)-Tmpr_Sur(n))                     ! the sensible heat conductance (J/m^2/day) 
        g_Heat_CD_total(1,n)=sqrt(varbt_Temp                    !cccz the surface heat conductance to assign new values for "heatmow"
     &    *(HeatCond_Mul_C+HeatCond_Mul_D))                     ! @ this step, the unit is J/m^2/day/K
     &    /100.00                                               ! change unit to J/cm^2/day/K   (only need 100.0D0 because of sqrt, and varBT_Temp has the correct unit)
        enddo
        g_Heat_Sensi(2,1)=0.0D0
        g_Heat_Sensi(2,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Heat_Sensi(2,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch(1,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch(1,n))
     &    *(MulchEleTmpr_temp(1,n)-MulchEleTmpr_temp(1,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo 
       elseif(k.eq.(mulchLayer+1)) then                           ! cccz for mulch-air interface              
        do n=1,SurNodeIndex_M-1
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          HeatCond_Mul_C=0.0D0
         else
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          else
            HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(k-1,n)
     &       *86400D0                                             ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
      else                                                        !cccz for mulch interior layers                            
       do n=1,SurNodeIndex_M-1
        if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)
         HeatCond_Mul_C=0.0D0
        else
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)
         if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &       *86400D0   
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &       *86400D0        
           HeatCond_Mul_C_Fo=
     &      (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &      /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &      *0.155D0
     &      *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &      *86400D0 
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex_M)=0.0D0
       do n=2,SurNodeIndex_M-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
       enddo 
       endif
      enddo
      
      layer_need_fix=0
      
cccz -------------Establish Governing Equations (solve mulch temperature and moisture)------------------------------------- 
2310  do k=1,mulchLayer
       do n=1,SurNodeIndex_M-1
        a11=WaterCapa_mulch(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    +c_water_s*rho_w*WaterCapa_mulch(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_water_s*rho_w*MulchEleThNew_temp(k,n)+                ! partial T/partial t within mulch solid, i.e., liquid water
     &    c_mulch_s*rho_mulch_b+                                    ! the wood part, the f_mulch_pore should already be adjusted, i.e., (1-f_mulch_pore)
     &    (c_air_s*rho_dryair+c_vap_s*VaporAct_mulch_D(k,n))     ! partial T/partial t for the air space other than mulch solid
     &    *f_mulch_pore+
     &    (c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +RainFallInput_temp(k+1,n)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep
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
      
cccz -----------------Convergence Check (Picard Iteration)-----------------------------------------------
cccz_test
c      if(TotalTime.ge.3.0D-004) then
c         aaaa=1 
c      endif
cccz reset the cumulative smoother
      if(IterMulch.eq.1) then
       do k=1,mulchLayer
        do n=1,SurNodeIndex_M-1
         MulchElehNew_temp3(k,n)=0.0D0
         MulchEleTmpr_temp3(k,n)=0.0D0  
        enddo
       enddo
      endif
cccz reset/initialize the error measurements
      maxhNew=0.0D0                                
      maxTmpr=0.0D0
      maxhNewDiff=0.0D0
      maxTmprDiff=0.0D0

      if(ForceForwards.eq.1) then
       do kk=1,layer_need_fix  
        meanMulchHnew=0.0D0
        do nn=1,SurNodeIndex_M-1
         meanMulchHnew=meanMulchHnew
     &    +MulchElehNew_temp2(kk,nn)*widthPerMulchUnit(nn)
        enddo
       meanMulchHnew=min(meanMulchHnew/totalMulchWidth,ThreHnew_HardFix)
        do nn=1,SurNodeIndex_M-1
         MulchElehNew_temp2(kk,nn)=meanMulchHnew
         MulchElehNew_temp3(kk,nn)=meanMulchHnew
        enddo
       enddo 
      endif
      
      if(localstep.le.minLocalStep) then           ! when the timestep is so small, jump out of the Picard iteration
         IterMulch=mulchLayer+1  
      endif
      if(IterMulch.le.mulchLayer) then             ! the Picard iteration has to be processed for times >= mulch layer, all mulch layer has to be updated
       do k=1,mulchLayer                           ! cccz protection for 1. hnew>0; 2. hnew=NAN; 3. Tmpr=NAN
        do n=1,SurNodeIndex_M-1
        ForceForwards=0   
        layer_need_fix=0
        if(MulchElehNew_temp2(k,n).gt.0.0D0) then
         if(LocalStep.gt.minLocalStep) then
          goto 2307                              ! keep reducing the time step
         else 
          if(g_Vapor(2*k-1,n).le.0.0D0) then ! upwards flux cause positive head
            meanG_vapor=0.0D0
            meanG_latent_heat=0.0D0
            do nn=1,SurNodeIndex_M-1
             meanG_vapor=meanG_vapor
     &        +g_Vapor(1,nn)*widthPerMulchUnit(nn)
             meanG_latent_heat=meanG_latent_heat
     &        +g_Heat_Latent(1,nn)*widthPerMulchUnit(nn)
            enddo
            meanG_vapor=meanG_vapor/totalMulchWidth
            meanG_latent_heat=meanG_latent_heat/totalMulchWidth
            do kk=1,k                        ! set one layer of flux above
             do nn=1,SurNodeIndex_M-1 
              g_Vapor(2*kk+1,nn)=meanG_vapor
              g_Heat_Latent(2*kk+1,nn)=meanG_latent_heat
              g_Vapor(2*kk,nn)=0.0D0
              g_Heat_Latent(2*kk,nn)=0.0D0
             enddo
            enddo
          else                               ! downwards flux cause positive
            meanG_vapor=0.0D0
            meanG_latent_heat=0.0D0
            do nn=1,SurNodeIndex_M-1
             meanG_vapor=meanG_vapor
     &        +g_Vapor(2*k-1,nn)*widthPerMulchUnit(nn)
             meanG_latent_heat=meanG_latent_heat
     &        +g_Heat_Latent(2*k-1,nn)*widthPerMulchUnit(nn)
            enddo
            meanG_vapor=meanG_vapor/totalMulchWidth
            meanG_latent_heat=meanG_latent_heat/totalMulchWidth
            do nn=1,SurNodeIndex_M-1 
              g_Vapor(2*k+1,nn)=min(g_Vapor(2*k+1,nn),meanG_vapor)
              g_Heat_Latent(2*k+1,nn)=
     &         min(g_Heat_Latent(2*k+1,nn),meanG_vapor)
              g_Vapor(2*k,nn)=0.0D0
              g_Heat_Latent(2*k,nn)=0.0D0
            enddo  
          endif    
          ForceForwards=1
          layer_need_fix=k 
          goto 2310
         endif
        endif
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) then
            goto 2307 
        endif
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) then
            goto 2307 
        endif
        enddo
       enddo
       do n=1,SurNodeIndex_M-1
        do k=1,mulchLayer
c        MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
        MulchElehNew_temp(k,n)=
     &     min(MulchElehNew_temp2(k,n),ThreHnew_HardFix)    
        MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
        MulchEleThNew_temp(k,n)=
     &     WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   
        enddo
       enddo
       IterMulch=IterMulch+1
       goto 2001
      else   ! in this case, (IterMulch.gt.mulchLayer)
       IF(IterMulch.le.MaxIter) then
        do n=1,SurNodeIndex_M-1                                 ! cccz protection for 1. hnew>0; 2. hnew=NAN; 3. Tmpr=NAN
        do k=1,mulchLayer
        ForceForwards=0   
        layer_need_fix=0
        if(MulchElehNew_temp2(k,n).gt.0.0D0) then
         if(LocalStep.gt.minLocalStep) then
          goto 2307                              ! keep reducing the time step
         else 
          if(g_Vapor(2*k-1,n).le.0.0D0) then ! upwards flux cause positive head
            meanG_vapor=0.0D0
            meanG_latent_heat=0.0D0
            do nn=1,SurNodeIndex_M-1
             meanG_vapor=meanG_vapor
     &        +g_Vapor(1,nn)*widthPerMulchUnit(nn)
             meanG_latent_heat=meanG_latent_heat
     &        +g_Heat_Latent(1,nn)*widthPerMulchUnit(nn)
            enddo
            meanG_vapor=meanG_vapor/totalMulchWidth
            meanG_latent_heat=meanG_latent_heat/totalMulchWidth
            do kk=1,k                        ! set one layer of flux above
             do nn=1,SurNodeIndex_M-1 
              g_Vapor(2*kk+1,nn)=meanG_vapor
              g_Heat_Latent(2*kk+1,nn)=meanG_latent_heat
              g_Vapor(2*kk,nn)=0.0D0
              g_Heat_Latent(2*kk,nn)=0.0D0
             enddo
            enddo
          else                               ! downwards flux cause positive
            meanG_vapor=0.0D0
            meanG_latent_heat=0.0D0
            do nn=1,SurNodeIndex_M-1
             meanG_vapor=meanG_vapor
     &        +g_Vapor(2*k-1,nn)*widthPerMulchUnit(nn)
             meanG_latent_heat=meanG_latent_heat
     &        +g_Heat_Latent(2*k-1,nn)*widthPerMulchUnit(nn)
            enddo
            meanG_vapor=meanG_vapor/totalMulchWidth
            meanG_latent_heat=meanG_latent_heat/totalMulchWidth
            do nn=1,SurNodeIndex_M-1 
              g_Vapor(2*k+1,nn)=min(g_Vapor(2*k+1,nn),meanG_vapor)
              g_Heat_Latent(2*k+1,nn)=
     &         min(g_Heat_Latent(2*k+1,nn),meanG_vapor)
              g_Vapor(2*k,nn)=0.0D0
              g_Heat_Latent(2*k,nn)=0.0D0
            enddo  
          endif    
          ForceForwards=1
          layer_need_fix=k 
          goto 2310
         endif
        endif
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) then
            goto 2307
        endif
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) then
            goto 2307
        endif
        enddo
       enddo
       do n=1,SurNodeIndex_M-1
        do k=1,mulchLayer
        MulchElehNew_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchElehNew_temp3(k,n)
     &     +MulchElehNew_temp2(k,n))/dble(IterMulch-mulchLayer)
cccz
        MulchElehNew_temp3(k,n)=
     &     min(MulchElehNew_temp3(k,n),ThreHnew_HardFix) 
cccz        
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
        

        
        MulchEleThNew_temp(k,n)=
     &     WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   
        enddo
       enddo
       IterMulch=IterMulch+1
cccz_try force forwards here to get rid of Picard Iteration when step=minLocalStep
       if(localstep.le.minLocalStep) then
          maxhNewDiff=0.0D0
          maxTmprDiff=0.0D0
       endif
       if(max(maxhNewDiff/maxhNew,maxTmprDiff/maxTmpr).gt.errTol) then
        goto 2001
       else
        TotalTime=TotalTime+LocalStep
        if(TotalTime-Step.gt.-0.01D0*dtMin) then
         if(TimeShrink.eq.0) then   
           goto 2302      ! if no 'local time step' was used, 'one step' was used for calculating water/energy exchange on soil surface
         else                       
           goto 2303      ! if multiple 'local time step' was used, then 'cumulative averaged' water/energy exchange were calculated in 2003 and 2004
         endif
        else
         do n=1,SurNodeIndex_M-1
          do k=1,mulchLayer
            MulchElehNew_temp0(k,n)=MulchElehNew_temp(k,n)
            MulchEleTmpr_temp0(k,n)=MulchEleTmpr_temp(k,n)

          enddo
         enddo
         goto 2303
2304     LocalStep=min(LocalStep*DMul1,Step-TotalTime)
         IterMulch=1
         goto 2001
        endif
       endif
      ELSE
2307   TimeShrink=1
       do n=1,SurNodeIndex_M-1
        do k=1,mulchLayer
         MulchElehNew_temp(k,n)=MulchElehNew_temp0(k,n)
         MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp0(k,n)
         MulchEleThNew_temp(k,n)=
     &     WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   
        enddo
       enddo
       LocalStep=max(LocalStep*DMul2,minLocalStep)
       IterMulch=1
       goto 2001
      ENDIF
      endif  !if(IterMulch.le.mulchLayer)
    
cccz ----------------Update the water and energy fluxes across soil surface (Whole Step)----------------------------------
c    first do the vapor flux calculation
c    unit within this module will be g/day/m^2
c    unit taken in watermov module is g/day/cm^2
2302  do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0        !cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,2)=-g_Vapor(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)                           !cccz update 'VarBW' for 'WaterMov': May or maynot be changed later in surface_water
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)    
        elseif(n.eq.qLeft.and.nn.gt.1) then                           !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +RainFallInput_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Vapor(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)   
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurR,2)=-g_Vapor(1,nn)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)   
        else
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0        !cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,2)=-g_Vapor(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)                           !cccz update 'VarBW' for 'WaterMov': May or maynot be changed later in surface_water
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)       
        endif      
       enddo  
      enddo
             

       
c    second update the heat fluxes
cccz the unit within this module will be J/day/m^2
cccz ********* ideas for the parameter assignment **************
c    1. VarBT(,1) is the temperature, as assign the temperature at mulch bottom surface (not important, I need fluxes more)
c    2. VarBT(,2) and VarBT(,3) are coupled, VarBT(,2) looks like a coefficient and VarBT(,3)=VarBT(,2)*temperature
c         The only temperature dependent thing is sensible heat in this code, so
c         VarBT(,2) sensible heat conductance,  VarBT(,3)= VarBT(,2)*(temperature at mulch bottom surface)
c    3. VarBT(,4) is the net radiation, we have short wave and long wave radiation here.
c    Latent heat is not presented in VarBT, it will depend on VarBW and soil temperature.
c    Make unit change to Cal/cm^2/min as shown in hourwea.for
cccz ********* end ****************
       
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
          kSurL=SurfNodeSurfIndexH(n)
          VarBT(kSurL,1)=MulchEleTmpr_temp(1,nn)
          VarBT(kSurL,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurL,3)=VarBT(kSurL,2)*MulchEleTmpr_temp(1,nn)
          VarBT(kSurL,4)=netRad(1,nn)*8.64D0                        ! 8.64D0=86400D0/10000D0
          varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
          varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
          varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
          varBT_Mulch(kSurL,4)=VarBT(kSurL,4)  
        elseif(n.eq.qLeft.and.nn.gt.1) then                        !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
          kSurL=SurfNodeSurfIndexH(n)
       VarBT(kSurL,1)=(MulchEleTmpr_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +MulchEleTmpr_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
        VarBT(kSurL,2)=(g_Heat_CD_total(1,nn-1)*widthPerMulchUnit(nn-1)
     &     +g_Heat_CD_total(1,nn)*widthPerMulchUnit(nn))
     &     /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
        VarBT(kSurL,4)=(netRad(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +netRad(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))*8.64D0  ! 8.64D0=86400D0/10000D0
        varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
        varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
        varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
        varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
          kSurR=SurfNodeSurfIndexH(n)
          VarBT(kSurR,1)=MulchEleTmpr_temp(1,nn)  
          VarBT(kSurR,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(1,nn)
          VarBT(kSurR,4)=netRad(1,nn)*8.64D0                        ! 8.64D0=86400D0/10000D0
          varBT_Mulch(kSurR,1)=VarBT(kSurR,1)
          varBT_Mulch(kSurR,2)=VarBT(kSurR,2)
          varBT_Mulch(kSurR,3)=VarBT(kSurR,3)
          varBT_Mulch(kSurR,4)=VarBT(kSurR,4)  
        else
          kSurL=SurfNodeSurfIndexH(n)
          VarBT(kSurL,1)=MulchEleTmpr_temp(1,nn)
          VarBT(kSurL,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurL,3)=VarBT(kSurL,2)*MulchEleTmpr_temp(1,nn)
          VarBT(kSurL,4)=netRad(1,nn)*8.64D0                        ! 8.64D0=86400D0/10000D0
          varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
          varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
          varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
          varBT_Mulch(kSurL,4)=VarBT(kSurL,4)  
        endif      
       enddo  
      enddo 
      

      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput_temp(n,k)
       enddo
      enddo
cccz old code      
c      do k=1,SurNodeIndex-1
c       do n=1,mulchLayer+1
c        RainFallInput(n,k)=RainFallInput_temp(n,k)
c       enddo
c      enddo
      goto 2306
      
cccz ------------------------Update the water and energy fluxes across soil surface (For Each Time Segment)--------------------------      
c    for each time segment during the iteration, use the same scheme showed following 'Index 2302'
cccz first the water part
2303  do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then   
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(1,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
        elseif(n.eq.qLeft.and.nn.gt.1) then                           !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
        Local_VarBW1=(RainFallInput_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +RainFallInput_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Local_VarBW2=-(g_Vapor(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Vapor(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep  
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(1,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurR)*Local_VarBW3)*LocalStep  
        else
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(1,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
        endif      
       enddo  
      enddo

      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         Local_VarBT1=MulchEleTmpr_temp(1,nn)
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=netRad(1,nn)*8.64D0                           ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        elseif(n.eq.qLeft.and.nn.gt.1) then                        !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         Local_VarBT1=(MulchEleTmpr_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +MulchEleTmpr_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
         Local_VarBT2=(g_Heat_CD_total(1,nn-1)*widthPerMulchUnit(nn-1)
     &     +g_Heat_CD_total(1,nn)*widthPerMulchUnit(nn))
     &     /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=(netRad(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +netRad(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))*8.64D0       ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         Local_VarBT1=MulchEleTmpr_temp(1,nn)
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=netRad(1,nn)*8.64D0                            ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        else
         Local_VarBT1=MulchEleTmpr_temp(1,nn)
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=netRad(1,nn)*8.64D0                           ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        endif      
       enddo  
      enddo 
      
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)
     &   +RainFallInput_temp(n,k)*LocalStep
      enddo
      enddo

      if(TotalTime-Step.gt.-0.01D0*dtMin) then
          goto 2305   ! finalize the iteration
      else
          goto 2304   ! process to the next local step
      endif

cccz ---------'cumulative averaged' water/energy exchange between 2003 and 2004------------------------
c    because that is for each time segment, we have to put them together.
c    determine the vapor and energy exchange on soil surface
c    Note: We already released from the mulch hori element and we should use soil hori element "SurNodeIndex" here
2305  do n=1,SurNodeIndex
       kSur=SurfNodeSurfIndexH(n)                         !cccz: water vapor part
       Varbw_Mulch(kSur,1)=VarBW1_temp(n)/TotalTime
       Varbw_Mulch(kSur,2)=VarBW2_temp(n)/TotalTime
       Varbw_Mulch(kSur,3)=VarBW3_temp(n)/TotalTime
       VarBW(kSur,1)=Varbw_Mulch(kSur,1)
       VarBW(kSur,2)=Varbw_Mulch(kSur,2)
       VarBW(kSur,3)=Varbw_Mulch(kSur,3)
       nNode=KXB(kSur)
       Q(nNode)=Q_temp(n)/TotalTime
       VarBT(kSur,1)=VarBT1_temp(n)/TotalTime             !cccz: energy part
       VarBT(kSur,2)=VarBT2_temp(n)/TotalTime
       VarBT(kSur,3)=VarBT3_temp(n)/TotalTime
       VarBT(kSur,4)=VarBT4_temp(n)/TotalTime
       VarBT_Mulch(kSur,1)=VarBT(kSur,1)
       VarBT_Mulch(kSur,2)=VarBT(kSur,2)
       VarBT_Mulch(kSur,3)=VarBT(kSur,3)
       VarBT_Mulch(kSur,4)=VarBT(kSur,4)
      enddo
cccz always note "RainFallInput" is local mulch variable that should be of mulch "hori" dimension
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)/TotalTime
       enddo
      enddo
      goto 2306

cccz all the varibles here are for mulch "hori" scale
2306  aaaa=2                                              !cccz final assignment
      do n=1,SurNodeIndex_M-1
       do k=1,mulchLayer
         MulchElehNew(k,n)=MulchElehNew_temp(k,n)
         MulchEleTmpr(k,n)=MulchEleTmpr_temp(k,n)
         MulchEleThNew(k,n)=MulchEleThNew_temp(k,n)
       enddo
      enddo
      goto 2600                                           ! goto the node-based computation, then return 
      
cccz --------------------------------PARTIALLY PONDED WATER (24XX Labels)-----------------------------------------------------------
c    computation has three steps: under ponded water surface/ponded water surface/open air
2400  aaaa=3
      do n=1,SurNodeIndex_M-1               
       do k=1,SubmergeIndex-1             !cccz under water part
        RelaHumid_mulch(k,n)=1.0D0
        VaporSat_mulch(k,n)=1.0D6    ! density of liquid water, we do not use this value anyway
        VaporAct_mulch_D(k,n)=1.0D6
        VaporAct_mulch_P(k,n)=0.0D0
        WaterDifici_mulch(k,n)=0.0D0
        VaporDiff_mulch(k,n)=0.0D0
        WaterCapa_mulch(k,n)=0.0D0
        HeatCapa_Air(k,n)=4.182D6    ! heat capacity (J/m^3/K)
        HeatDiff_mulch(k,n)=0.0D0
       enddo 
       k=SubmergeIndex                    !cccz right @ water level
       if(PerOccupation.eq.2.0D0) then        ! the water level is high enough, this layer is now a water layer
        RelaHumid_mulch(k,n)=1.0D0
        VaporSat_mulch(k,n)=1.0D6        ! density of liquid water, we do not use this value anyway
        VaporAct_mulch_D(k,n)=1.0D6 
        VaporAct_mulch_P(k,n)=0.0D0
        WaterDifici_mulch(k,n)=0.0D0
        VaporDiff_mulch(k,n)=0.0D0
        WaterCapa_mulch(k,n)=0.0D0
        HeatCapa_Air(k,n)=4.182D6        ! heat capacity (J/m^3/K)
        HeatDiff_mulch(k,n)=0.0D0
       else                                   ! the water level is not high enough, calculate vapor for the rest of this layer 
        RelaHumid_mulch(k,n)=exp(2.124D-4*MulchElehNew_temp(k,n)     ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))
         RelaHumid_mulch(k,n)=min(RelaHumid_mulch(k,n),1.0D0)   ! make sure the relative humidity is of currect range      
         RelaHumid_mulch(k,n)=max(RelaHumid_mulch(k,n),0.0D0)
        VaporSat_mulch(k,n)=exp(19.84D0-4975.9D0                     ! Saturated Vapor density for each mulch element (g/m^3)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))    
        VaporAct_mulch_D(k,n)=VaporSat_mulch(k,n)                 ! Actual Vapor density for each mulch element (g/m^3)
     &   *RelaHumid_mulch(k,n) 
        VaporAct_mulch_P(k,n)=
     &   (0.61D0*exp((17.27D0*MulchEleTmpr_temp(k,n))
     &   /(MulchEleTmpr_temp(k,n)+237.3D0)))                            ! present this in Pa, not kPa, 1000.0: change unit from kPa to Pa
     &   *RelaHumid_mulch(k,n)*1000.0D0   
cccz water difficit in each element, including liquid water and water vapor to a certain degree of saturation
c    vapor part: reach to rh=1.0D0
c    solid/liquid part: no need to change, balance with the vapor part.
        WaterDifici_mulch(k,n)=
     &    (VaporSat_mulch(k,n)-VaporAct_mulch_D(k,n))             ! rate to saturate the vapor portion (g/m^2/day)
     &    *(thickPerLayer(k)*PerOccupation/100.0D0)/LocalStep             ! need to adjust the height "thickPerLayer(k)*PerOccupation"
     &    *f_mulch_pore         
cccz water vapor diffusivity calculation, use diffusion type water fluxes
        VaporDiff_mulch(k,n)=2.35D-5*
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**1.667D0)*86400.0D0                             ! Water Vapor diffusivity (m^2/day), (f_mulch_pore**0.667D0) is the turtorosity
cccz water capacity: only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore
c    similar to the water capacity calculation for soil        
         WaterCapa_mulch(k,n)=
     &     WC_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   
cccz heat capacity and diffusivity calculation, air density (rho_dryair) g m^-3, other parameters 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air(k,n)=1.85D0*VaporAct_mulch_D(k,n)
     &   +1.006D0*rho_dryair                                 ! heat capacity (J/m^3/K)
        HeatDiff_mulch(k,n)=HeatCapa_Air(k,n)*1.84896D0      ! heat diffusivity (J/m/K/day), more like a conductance 1.84896D0=2.14D-5(heat diff)*86400.0D0       
     &   +HeatDiff_fabric*86400.0D0*(1-f_mulch_pore)  
       endif      
       do k=SubmergeIndex+1,mulchLayer                    !cccz above water 
cccz 'vapor density (g/m^3)'/'vapor partial pressure (Pa)'
        RelaHumid_mulch(k,n)=max(min(exp(2.124D-4*MulchElehNew_temp(k,n)     ! rh=exp(gMh/RT), 2.124D-4=(8.91m/s^2)*(0.018kg/mol)/(100 cm/m)/(8.314)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0)),1.0D0),0.0D0)
        VaporSat_mulch(k,n)=exp(19.84D0-4975.9D0                     ! Saturated Vapor density for each mulch element (g/m^3)
     &   /(MulchEleTmpr_temp(k,n)+273.15D0))           
        VaporAct_mulch_D(k,n)=VaporSat_mulch(k,n)                 ! Actual Vapor density for each mulch element (g/m^3)
     &   *RelaHumid_mulch(k,n)
        VaporAct_mulch_P(k,n)=
     &   (0.61D0*exp((17.27D0*MulchEleTmpr_temp(k,n))
     &   /(MulchEleTmpr_temp(k,n)+237.3D0)))                            ! present this in Pa, not kPa, 1000.0: change unit from kPa to Pa
     &   *RelaHumid_mulch(k,n)*1000.0D0   
cccz water difficit in each element, including liquid water and water vapor to a certain degree of saturation
c    vapor part: reach to rh=1.0D0
c    solid/liquid part: no need to change, they will make balance with the vapor parts
        WaterDifici_mulch(k,n)=
     &    (VaporSat_mulch(k,n)-VaporAct_mulch_D(k,n))
     &    *(thickPerLayer(k)/100.0D0)/LocalStep               ! rate to saturate the vapor portion (g/m^2/day)
     &    *f_mulch_pore                                       ! cast the rates to gas phase of the mulch
cccz water diffusivity calculation, intrinsic diffusivity * temperature effects
        VaporDiff_mulch(k,n)=2.35D-5*                    ! 2.35=(2.29_soil science+2.40_[campbell book])/2
     &   ((1.0D0+MulchEleTmpr_temp(k,n)/273.15D0)**1.75D0)
     &   *(f_mulch_pore**1.667D0)*86400.0D0                   ! Water Vapor diffusivity (m^2/day)
cccz water capacity: only for the mulch 'solid material' part, but need to be casted into the full mulch elements use f_mulch_pore
c    use soil water capacity.
         WaterCapa_mulch(k,n)=
     &     WC_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &     LIGN_mass_Frac(k),lInput)   
cccz heat capacity and diffusivity calculation: air density (rho_dryair) g m^-3, other parameters 1.006D0, 1.85D0 are in 'kJ/kg/K' for air and water vapor  
        HeatCapa_Air(k,n)=1.85D0*VaporAct_mulch_D(k,n)       ! heat capacity (J/m^3/K)
     &   +1.006D0*rho_dryair
        HeatDiff_mulch(k,n)=HeatCapa_Air(k,n)*1.84896D0      ! heat diffusivity (J/m/K/day), 1.84896D0=2.14D-5(heat diff)*86400.0D0
     &   +HeatDiff_fabric*86400.0D0*(1-f_mulch_pore)
       enddo
       Ra_mulch(1,n)=9.81D0                                 ! Rayleigh num
     &  *abs(0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &  -0.5D0*(MulchEleTmpr_temp(SubmergeIndex,n)+Tmpr_Sur(n)))
     &  *((mulchThick/100.0D0)**3.0D0)
     &  /nu_air
     &  /((546.30D0
     &  +0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &  +0.5D0*(MulchEleTmpr_temp(SubmergeIndex,n)+Tmpr_Sur(n)))/2.0D0)
     &  /Dhm
       GrRe_mulch(1,n)=2.0D0*9.81D0
     &  *(mulchThick/100.0D0)                               ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
     &  *abs(0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &  -0.5D0*(MulchEleTmpr_temp(SubmergeIndex,n)+Tmpr_Sur(n)))
     &  /(546.30D0
     &  +0.5D0*(MulchEleTmpr_temp(mulchLayer,n)+AirTemp_Wea)
     &  +0.5D0*(MulchEleTmpr_temp(SubmergeIndex,n)+Tmpr_Sur(n)))
     &  /((u_mulch(mulchLayer)-u_mulch(SubmergeIndex))**2.0D0)  
      enddo

cccz ------------------ water (vapor) fluxes within mulch----------------------------
c    only vapor flux occur, but allow some liquid water stored in the mulch solid materials.
c    downward and leftward fluxes are positive
c    unit of water vapor flux: g/m^2/day
      if(PerOccupation.eq.2.0D0) then ! the water level is high enough, this layer is now a water layer, zero the vapor fluxes for all the under water layers
       do n=1,SurNodeIndex_M
        do k=1,SubmergeIndex+1
         g_vapor(2*k-1,n)=0.0D0
         g_vapor(2*k,n)=0.0D0
        enddo
       enddo
       do k=SubmergeIndex+1,mulchLayer+1    
       if(k.eq.(SubmergeIndex+1)) then    !cccz similar to mulch-soil interface, but now it is water-mulch interface             
        do n=1,SurNodeIndex_M-1             !cccz the soil surface vapor flux based on the min('potential evaporation' method, 'convection/diffusion' method), i.e., how much soil can supply, how much mulch can transport.
         ! METHOD1: 'potential evaporation' method
         SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(k,n))
     &    /(MulchEleTmpr_temp(k,n)+237.3D0))
         DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(k,n)+1.0D0))
     &    /(MulchEleTmpr_temp(k,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch(k,n))          ! calculate VPD for the mulching layer above soil surface
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(k,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur
     &    *max(netRad(k,n),0.0D0)*3600.0D0/(2500.8D0
     &    -(2.3668D0*MulchEleTmpr_temp(k,n))))+(VPD_Sur*109.375D0
     &    *(1.0D0+(0.149D0*u_mulch(k)))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near water-mulch interface (g/m2/day), 24.0 changes unit from g/m^2/hour to g/m^2/day 
         ! METHOD2: 'convection/diffusion' method
c         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))    ! Rayleigh num
c     &    *((0.5D0*thickPerLayer(k)/100.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0*(0.5D0*thickPerLayer(k)/100.0D0) ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))
c     &    /(546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))
c     &    /((u_mulch(k)-u_mulch(k-1))**2.0D0)
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then   ! diffusive vapor flow always happens
          VaporCond_Mul_D=VaporDiff_mulch(k,n)
     &      /(thickPerLayer(k)/200.0D0)
          VaporCond_Mul_C=0.0D0
         else
          VaporCond_Mul_D=VaporDiff_mulch(k,n)
     &      /(thickPerLayer(k)/200.0D0)
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then                   ! the Ri is large, means Re^2 is small, thus turbulence is small, hence free convection is processed, forced convection is neglected
            VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400D0                                        ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
            VaporCond_Mul_C_Fo=0.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)   
          else                                                        ! both free and forced convection should be considered
            VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400.0D0
            VaporCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log((LayerHeight(k)+0.5D0*thickPerLayer(k))
c     &       /(0.079D0*(LayerHeight(k)+0.5D0*thickPerLayer(k))))
c     &       **2.0D0)
     &       *0.155D0
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &       *(f_mulch_pore**1.667D0)            ! establish a path that similar to the diffusivity, using tortuosity and pore_spaces
          endif
         endif
         g_vapor_try2=VaporCond_Mul_D*
     &    (VaporAct_mulch_D(k,n)-VaporAct_Sur(n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(k,n)-VaporAct_Sur_P(n))
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then      !cccz take the min between 'diffusive type' and 'potential' evporation, if the two values are of different directions, we assume zero flux fo numerical stable     
           g_vapor(2*k-1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(2*k-1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(2*k-1,n)=0.0D0
         endif
        enddo      
        g_Vapor(2*k,1)=0.0D0           ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/VaporDiff_mulch(k,n-1)
     &    +widthPerMulchUnit(n)/VaporDiff_mulch(k,n))
     &    *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
         if(g_Vapor(2*k,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2*k,n)=
     &       min(g_Vapor(2*k,n),WaterDifici_mulch(k,n-1))
         else
           g_Vapor(2*k,n)=
     &       max(g_Vapor(2*k,n),-WaterDifici_mulch(k,n))  
         endif
        enddo   
       elseif(k.eq.(mulchLayer+1)) then                             !cccz for mulch-air interface                        
        do n=1,SurNodeIndex_M-1                                       ! use conduction/convection to calculate conductivity at atmosphere surface
c         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)   ! Rayleigh num
c     &    *(((thickPerLayer(k-1))/100.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0*((z_wind-mulchThick)/100.0D0)       ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
c     &    /(546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)
c     &    /((AirWind_Wea-u_mulch(k))**2.0D0)     
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
           VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0) 
           VaporCond_Mul_C=0.0D0
         else
           VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0) 
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
            VaporCond_Mul_C_Fr=a_free                                 ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0            
            VaporCond_Mul_C_Fo=0.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &        *(f_mulch_pore**1.667D0)        ! cccz: something compared with the diffusivity, adjust the accessible for path, very arbitrary, like a guess ??????????????????
          else
           VaporCond_Mul_C_Fr=a_free
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0
           VaporCond_Mul_C_Fo=
     &        (karman**2.0D0)*u_mulch(k)
c     &        /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &        *(f_mulch_pore**1.667D0)        ! cccz: something compared with the diffusivity, adjust the accessible for path, very arbitrary, like a guess ??????????????????
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D*
     &    *(VaporAct_Ambient-VaporAct_mulch_D(k-1,n))
         AirVaporP_temp=AirVaporP_Wea*1000.0D0
         if(AirVaporP_temp.le.VaporAct_mulch_P(k-1,n)) then  ! upwards flow
          g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &     VaporCond_Mul_C*
     &     (AirVaporP_temp-VaporAct_mulch_P(k-1,n))
          if(localstep.le.minLocalStep) then
            g_Vapor(2*k-1,n)=max(g_Vapor(2*k-1,n),
     &       min(g_Vapor(2*k-3,n),0.0D0)-WaterDifici_mulch(k-1,n))
          endif  
         else
          g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &     min(VaporCond_Mul_C*
     &     (AirVaporP_temp-VaporAct_mulch_P(k-1,n)),
     &     WaterDifici_mulch(k-1,n))    
         endif
        enddo  
       else              !cccz for water flux within mulch                                             
        do n=1,SurNodeIndex_M-1
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
           VaporCond_Mul_C=0.0D0
         else
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)  
           VaporCond_Mul_C_Fo=0.0D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results) 
           VaporCond_Mul_C_Fo=
     &       (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)   
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D
     &    *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k-1,n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(k,n)-VaporAct_mulch_P(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch(k,n))
     &     *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
          if(g_Vapor(2*k,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2*k,n)
     &     =min(g_Vapor(2*k,n),WaterDifici_mulch(k,n-1))
          else
           g_Vapor(2*k,n)
     &     =max(g_Vapor(2*k,n),-WaterDifici_mulch(k,n))  
          endif
        enddo
       endif
       enddo
      else   !cccz: the water level is not high enough, the occupied layer is now a fraction of the original layers
       do k=1,SubmergeIndex+1  ! zero all the under water layers, to "" could be more than needed but not a problem, because we calcualte it later
        do n=1,SurNodeIndex_M
         g_vapor(2*k-1,n)=0.0D0
         g_vapor(2*k,n)=0.0D0
        enddo
       enddo    
       do k=SubmergeIndex,mulchLayer+1    
       if(k.eq.SubmergeIndex) then    ! the "shrinked" water surface layer, but take the fraction of the thickness
        thickLayer_Sur=thickPerLayer(k)*PerOccupation
        if(k.eq.1) then
            u_SurWatermulch=u_soilsur
        else
            u_SurWatermulch=u_mulch(k-1)
        endif
        do n=1,SurNodeIndex_M-1
         ! METHOD1: 'potential evaporation' method
         SVPA_Sur=0.61D0*EXP((17.27D0*MulchEleTmpr_temp(k,n))
     &    /(MulchEleTmpr_temp(k,n)+237.3D0))
         DEL_Sur=(0.61D0*EXP((17.27D0*(MulchEleTmpr_temp(k,n)+1.0D0))
     &    /(MulchEleTmpr_temp(k,n)+1.0D0+237.3D0)))-SVPA_Sur
         VPD_Sur=SVPA_Sur*(1.0D0-RelaHumid_mulch(k,n))          ! calculate VPD for the first mulching layer
         D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))                   ! we use actural vapor pressure for D31 here, should be the saturated vapor pressure of wet bulb temp in air
         D32=2500.8D0-2.37D0*MulchEleTmpr_temp(k,n)
         GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
         g_Vapor_try1=-((DEL_Sur/GAMMA_Sur*max(netRad(k,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*MulchEleTmpr_temp(k,n))))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_mulch(k)))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 
         ! METHOD2: 'convection/diffusion' method
c         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))
c     &    *((thickLayer_Sur/200.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0*(thickLayer_Sur/200.0D0) ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))
c     &    /(546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))
c     &    /((u_mulch(k)-u_SurWatermulch)**2.0D0)
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then   ! diffusive vapor flow always happens
           VaporCond_Mul_D=VaporDiff_mulch(k,n)                  ! the vapor conductance (m/day)
     &      /(thickLayer_Sur/200.0D0) 
           VaporCond_Mul_C=0.0D0                                       
         else
           VaporCond_Mul_D=VaporDiff_mulch(k,n)                  ! the vapor conductance (m/day)
     &      /(thickLayer_Sur/200.0D0)    
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then                   ! the Ri is large, means Re^2 is small, thus turbulence is small, hence free convection is processed, forced convection is neglated
            VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400D0                                        ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
            VaporCond_Mul_C_Fo=0.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &       *(f_mulch_pore**1.667D0)            ! establish a path that similar to the diffusivity, using tortuosity and pore_spaces    
          else                                                        ! both free and forced convection should be considered
            VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400D0
            VaporCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log((LayerHeight(k)+0.5D0*thickPerLayer(k))
c     &       /(0.079D0*(LayerHeight(k)+0.5D0*thickPerLayer(k))))
c     &       **2.0D0)
     &       *0.155D0
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)+Tmpr_Sur(n))/2.0D0)
     &       *1000.0D0*86400D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
c     &       *(f_mulch_pore**1.667D0)            ! establish a path that similar to the diffusivity, using tortuosity and pore_spaces
          endif
         endif
         g_vapor_try2=VaporCond_Mul_D*
     &    (VaporAct_mulch_D(k,n)-VaporAct_Sur(n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(k,n)-VaporAct_Sur_P(n))
         if(g_Vapor_try1.lt.0.0D0.and.g_vapor_try2.lt.0.0D0) then    !cccz take the min between 'diffusive type' and 'potential' evporation, if the two values are of different directions, we assume zero flux fo numerical stable
           g_vapor(2*k-1,n)=max(g_Vapor_try1, g_Vapor_try2)
         elseif(g_Vapor_try1.gt.0.0D0.and.g_vapor_try2.gt.0.0D0) then
           g_vapor(2*k-1,n)=min(g_Vapor_try1, g_Vapor_try2)
         else
           g_vapor(2*k-1,n)=0.0D0
         endif
        enddo      
        g_Vapor(2*k,1)=0.0D0           ! horizontal vapor flux are in even rows, and has impermeable boundaries    
        g_Vapor(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/VaporDiff_mulch(k,n-1)
     &    +widthPerMulchUnit(n)/VaporDiff_mulch(k,n))
     &    *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)  ! horizontal vapor flux near soil-mulch interface (g/m2/day)
         if(g_Vapor(2*k,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2*k,n)=
     &       min(g_Vapor(2*k,n),WaterDifici_mulch(k,n-1))
         else
           g_Vapor(2*k,n)=
     &       max(g_Vapor(2*k,n),-WaterDifici_mulch(k,n))  
         endif
        enddo        
       elseif(k.eq.(mulchLayer+1)) then           !cccz for mulch-air interface                    
        do n=1,SurNodeIndex_M-1                                        ! use conduction/convection to calculate conductivity at atmosphere surface
c         Ra_mulch(k,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
c     &    *(((thickPerLayer(k-1))/100.0D0)**3.0D0)/nu_air
c     &    /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0*((z_wind-mulchThick)/100.0D0)  ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
c     &    /(546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)
c     &    /((AirWind_Wea-u_mulch(k))**2.0D0)     
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0) 
           VaporCond_Mul_C=0.0D0
         else
          VaporCond_Mul_D=sqrt(2.35D-5                               ! the vapor conductance (m/day); geometrical average of upper layer & open air
     &      *((1.0D0+AirTemp_Wea/273.15D0)**1.75D0)
     &      *86400.0D0*VaporDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0) 
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           VaporCond_Mul_C_Fr=a_free                                 ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400.0D0            
           VaporCond_Mul_C_Fo=0.0D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else
           VaporCond_Mul_C_Fr=a_free
     &        *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400.0D0
           VaporCond_Mul_C_Fo=
     &        (karman**2.0D0)*u_mulch(k)
c     &        /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &        *M_H2O/R_gas
     &        /((546.30D0+AirTemp_Wea+MulchEleTmpr_temp(k-1,n))/2.0D0)
     &        *1000.0D0*86400.0D0
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D*
     &    *(VaporAct_Ambient-VaporAct_mulch_D(k-1,n))
         AirVaporP_temp=AirVaporP_Wea*1000.0D0
         if(AirVaporP_temp.le.VaporAct_mulch_P(k-1,n)) then  ! upwards flow
           g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &      VaporCond_Mul_C*
     &      (AirVaporP_temp-VaporAct_mulch_P(k-1,n))
           if(localstep.le.minLocalStep) then
            g_Vapor(2*k-1,n)=max(g_Vapor(2*k-1,n),
     &       min(g_Vapor(2*k-3,n),0.0D0)-WaterDifici_mulch(k-1,n))  
           endif  
         else
          g_Vapor(2*k-1,n)=g_Vapor(2*k-1,n)+
     &     min(VaporCond_Mul_C*
     &     (AirVaporP_temp-VaporAct_mulch_P(k-1,n)),
     &     WaterDifici_mulch(k-1,n))    
         endif
        enddo  
       else                       !cccz for water flux within mulch                                                 
        do n=1,SurNodeIndex_M-1
c         Ra_mulch(k,n)=9.81D0
c     &    *abs(MulchEleTmpr_temp(k-1,n)-MulchEleTmpr_temp(k,n))
c     &    *(((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)**3.0D0)
c     &    /nu_air/((546.30D0+MulchEleTmpr_temp(k-1,n)
c     &    +MulchEleTmpr_temp(k,n))/2.0D0)/Dhm
c         GrRe_mulch(k,n)=2.0D0*9.81D0
c     &    *((thickperlayer(k)+thickperlayer(k-1))/200.0D0)    ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
c     &    *abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n))
c     &    /(546.30D0+MulchEleTmpr_temp(k,n)+MulchEleTmpr_temp(k-1,n))
c     &    /((u_mulch(k)-u_mulch(k-1))**2.0D0)  
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)
           VaporCond_Mul_C=0.0D0
         else
          VaporCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))/
     &     (thickPerLayer(k)/VaporDiff_mulch(k,n)+
     &     thickPerLayer(k-1)/VaporDiff_mulch(k-1,n))
     &     /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)   
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)  
           VaporCond_Mul_C_Fo=0.0D0
           VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          else
           VaporCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results) 
           VaporCond_Mul_C_Fo=
     &       (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *M_H2O/R_gas
     &       /((546.30D0+MulchEleTmpr_temp(k,n)
     &       +MulchEleTmpr_temp(k-1,n))/2.0D0)
     &       *1000.0D0*86400.0D0                                    ! change unit from "kg/m^3/s" to "g/m^3/day" (for the final flux results)   
            VaporCond_Mul_C=(VaporCond_Mul_C_Fr+VaporCond_Mul_C_Fo)
          endif
         endif
         g_Vapor(2*k-1,n)=VaporCond_Mul_D
     &    *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k-1,n))
     &    +VaporCond_Mul_C*
     &    (VaporAct_mulch_P(k,n)-VaporAct_mulch_P(k-1,n))
        enddo
        g_Vapor(2*k,1)=0.0D0
        g_Vapor(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1 
          g_Vapor(2*k,n)=(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &     (widthPerMulchUnit(n-1)/VaporDiff_mulch(k,n-1)
     &     +widthPerMulchUnit(n)/VaporDiff_mulch(k,n))
     &     *(VaporAct_mulch_D(k,n)-VaporAct_mulch_D(k,n-1))
     &     /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
          if(g_Vapor(2*k,n).ge.0.0D0) then  ! set a horizontal flux limter
           g_Vapor(2*k,n)
     &     =min(g_Vapor(2*k,n),WaterDifici_mulch(k,n-1))
          else
           g_Vapor(2*k,n)
     &     =max(g_Vapor(2*k,n),-WaterDifici_mulch(k,n))  
          endif
        enddo
       endif
       enddo    
      endif
               
cccz ------------------------Rainfall Redistribution--------------------------------------
c    the under ponded water, 'WaterDifici_mulch=0.0D0' and the rainfall directly pass to the soil surface
c    unit "g/m^2/day"
      do n=1,SurNodeIndex_M-1
       k=mulchLayer+1
2408   if(k.gt.1) then
        inputPerLayer=
     &   min(RainFallInput_temp(k,n),WaterDifici_mulch(k-1,n))
        if(inputPerLayer.lt.RainFallInput_temp(k,n)) then
         RainFallInput_temp(k-1,n)=RainFallInput_temp(k,n)-inputPerLayer
         RainFallInput_temp(k,n)=inputPerLayer
         k=k-1
         goto 2408
        else
         continue
        endif
       else
        continue
       endif
      enddo       
      
cccz -----------------------Latent heat flux----------------------------------------------
c    downward and leftward are positive, energy income
c    calculate the fluxes based on the upwind direction of water
c    unit "J/m^2/day"    
      if(PerOccupation.eq.2.0D0) then     ! the water level is high enough, this layer is now a water layer, zero all the under water layers
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex_M
         g_Heat_Latent(2*k-1,n)=0.0D0
         g_Heat_Latent(2*k,n)=0.0D0
        enddo
       enddo
       do k=SubmergeIndex+1,mulchLayer+1  
       if(k.eq.(SubmergeIndex+1)) then    !cccz for water surface layer
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
           g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(Tmpr_Sur(n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                             ! leftwards vapor flow
          g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         else                                                         ! rightwards vapor flow
          g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
       elseif(k.eq.(mulchLayer+1)) then     !cccz for mulch-air interface layer
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
       else                                 !cccz for other interier layers
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
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
      else   ! the water level is not high enough, the occupied layer is now a fraction of the original layers, zero all the under water layers, maybe more than needed
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex_M
         g_Heat_Latent(2*k-1,n)=0.0D0
         g_Heat_Latent(2*k,n)=0.0D0
        enddo
       enddo
       do k=SubmergeIndex,mulchLayer+1
       if(k.eq.SubmergeIndex) then    !cccz for soil surface layer 
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &      *(Tmpr_Sur(n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         if(g_Vapor(2*k,n).ge.0.0D0) then                             ! leftwards vapor flow
          g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         else                                                         ! rightwards vapor flow
          g_Heat_Latent(2*k,n)=c_vap_s*g_Vapor(2*k,n)
     &     *(MulchEleTmpr_temp(k,n-1)-Tmpr_ref)+c_vap_v*g_Vapor(2*k,n)
         endif
        enddo
       elseif(k.eq.(mulchLayer+1)) then   !cccz for mulch-air interface layer
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
       else                               !cccz for other interier layers
        do n=1,SurNodeIndex_M-1
         if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         else                                                         ! upwards vapor flow
          g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &     *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
         endif
        enddo
        g_Heat_Latent(2*k,1)=0.0D0
        g_Heat_Latent(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
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

cccz ---------------Diffusion/Convective Heat Flux(Sensible Heat)----------------------
c    downward and leftward are positive, i.e., energy income
c    unit (J/m^2/day) 
      if(PerOccupation.eq.2.0D0) then ! the water level is high enough, this layer is now a water layer, zero all the under water layers
       do k=1,SubmergeIndex+1
        do n=1,SurNodeIndex_M
         g_Heat_Sensi(2*k-1,n)=0.0D0
         g_Heat_Sensi(2*k,n)=0.0D0
        enddo
       enddo 
       do k=SubmergeIndex+1,mulchLayer+1  
       if(k.eq.(SubmergeIndex+1)) then     !cccz for mulch-soil interface                     
        do n=1,SurNodeIndex_M-1
         qLeft=SurfNodeSurfIndexH(n)
         qRight=SurfNodeSurfIndexH(n+1)
         varbt_Temp=0.5D0*(varbt_Air(qLeft,2)+varbt_Air(qRight,2))
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=HeatDiff_mulch(k,n)
     &     /(thickPerLayer(k)/200.0D0)              ! the heat conductance (J/m^2/day/K), we will assume the "3rd dimension" has scale 1 cm
          HeatCond_Mul_C=0.0D0
         else
          HeatCond_Mul_D=HeatDiff_mulch(k,n)
     &     /(thickPerLayer(k)/200.0D0)              ! the heat conductance (J/m^2/day/K), we will assume the "3rd dimension" has scale 1 cm
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log((LayerHeight(k)+0.5D0*thickPerLayer(k))
c     &       /(0.079D0*(LayerHeight(k)+0.5D0*thickPerLayer(k))))
c     &       **2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                              ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=
     &    100.0D0*sqrt(varbt_Temp*(HeatCond_Mul_D+HeatCond_Mul_C))
c     &   (2.0D0/(1.0D0/(HeatCond_Mul_D+HeatCond_Mul_C)+1.0D0/Soil_Cond)) ! do not need consider the 
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))                         ! the sensible heat conductance (J/m^2/day) 
         g_Heat_CD_total(1,n)=sqrt(varbt_Temp                    !cccz the surface heat conductance to assign new values for "heatmow"
     &    *(HeatCond_Mul_C+HeatCond_Mul_D))                     ! @ this step, the unit is J/m^2/day/K
     &    /100.00                                               ! change unit to J/cm^2/day/K   
        enddo
        g_Heat_Sensi(2*k,1)=0.0D0
        g_Heat_Sensi(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Heat_Sensi(2*k,n)=
     &    (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch(k,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch(k,n))
     &    *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
       elseif(k.eq.(mulchLayer+1)) then          ! cccz for mulch-air interface                       
        do n=1,SurNodeIndex_M-1
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          HeatCond_Mul_C=0.0D0
         else
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)   
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400.0D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(k-1,n)
     &       *86400.0D0                                             ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
      else          !cccz for mulch interior layers                           
       do n=1,SurNodeIndex_M-1
        if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)
         HeatCond_Mul_C=0.0D0
        else
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &   /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &   /HeatDiff_mulch(k-1,n))
     &   /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)
         if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &       *86400.0D0   
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &       *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &       *86400.0D0          
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &       *86400.0D0 
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex_M)=0.0D0
       do n=2,SurNodeIndex_M-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo 
       endif
      enddo
      else    ! the water level is not that high enough within that layer, zero all the under water layers
       do n=1,SurNodeIndex_M
        do k=1,SubmergeIndex+1
          g_Heat_Sensi(2*k-1,n)=0.0D0
          g_Heat_Sensi(2*k,n)=0.0D0
        enddo
       enddo    
       do k=SubmergeIndex,mulchLayer+1  
        if(k.eq.SubmergeIndex) then    !cccz for mulch-soil interface, should adjust the layer height
        thickLayer_Sur=thickPerLayer(k)*PerOccupation
         do n=1,SurNodeIndex_M-1
c          qLeft=SurfNodeNodeIndexH(n)
c          qRight=SurfNodeNodeIndexH(n+1)
c          Tmpr_Sur=AirTemp_Wea
          qLeft=SurfNodeSurfIndexH(n)
          qRight=SurfNodeSurfIndexH(n+1)
          varbt_Temp=0.5D0*(varbt_Air(qLeft,2)+varbt_Air(qRight,2))
          if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
            HeatCond_Mul_D=HeatDiff_mulch(k,n)
     &        /(thickLayer_Sur/200.0D0)                       ! the heat conductance (J/m/day/K)
            HeatCond_Mul_C=0.0D0
          else
            HeatCond_Mul_D=HeatDiff_mulch(k,n)
     &        /(thickLayer_Sur/200.0D0)                       ! the heat conductance (J/m/day/K)
           if(GrRe_mulch(1,n).ge.GrRe_Critical) then
            HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                       ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
            HeatCond_Mul_C_Fo=0.0D0
            HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
           else
            HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n)))
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                                 ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
            HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log((LayerHeight(k)+0.5D0*thickPerLayer(k))
c     &       /(0.079D0*(LayerHeight(k)+0.5D0*thickPerLayer(k))))
c     &       **2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(k,n)
     &       *86400.0D0                                                 ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
            HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
          endif
         g_Heat_Sensi(2*k-1,n)=
     &    100.0D0*sqrt(varbt_Temp*(HeatCond_Mul_D+HeatCond_Mul_C))
c     &   (2.0D0/(1.0D0/(HeatCond_Mul_D+HeatCond_Mul_C)+1.0D0/Soil_Cond)) ! do not need consider the 
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_Sur(n))                               ! the sensible heat conductance (J/m^2/day) 
         g_Heat_CD_total(1,n)=sqrt(varbt_Temp                    !cccz the surface heat conductance to assign new values for "heatmow"
     &    *(HeatCond_Mul_C+HeatCond_Mul_D))                     ! @ this step, the unit is J/m^2/day/K
     &    /100.00                                               ! change unit to J/cm^2/day/K
        enddo
        g_Heat_Sensi(2*k,1)=0.0D0
        g_Heat_Sensi(2*k,SurNodeIndex_M)=0.0D0
        do n=2,SurNodeIndex_M-1
         g_Heat_Sensi(2*k,n)=
     &    (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
     &    /(widthPerMulchUnit(n-1)/HeatDiff_mulch(k,n-1)
     &    +widthPerMulchUnit(n)/HeatDiff_mulch(k,n))
     &    *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &    /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
        enddo
       elseif(k.eq.(mulchLayer+1)) then   !cccz for mulch-air interface                               
        do n=1,SurNodeIndex_M-1
         if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          HeatCond_Mul_C=0.0D0
         else
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          if(GrRe_mulch(1,n).ge.GrRe_Critical) then
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400.0D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=0.0D0
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          else
           HeatCond_Mul_C_Fr=a_free
     &       *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &       *HeatCapa_Air(k-1,n)
     &       *86400.0D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C_Fo=
     &       (karman**2.0D0)*u_mulch(k)
c     &       /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &       *HeatCapa_Air(k-1,n)
     &       *86400.0D0                                             ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
           HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
          endif
         endif
         g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
        enddo
      else                                                        !cccz for mulch interior layers                             
       do n=1,SurNodeIndex_M-1
        if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &    /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &    /HeatDiff_mulch(k-1,n))
     &    /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)
         HeatCond_Mul_C=0.0D0
        else
         HeatCond_Mul_D=(thickPerLayer(k)+thickPerLayer(k-1))
     &    /(thickPerLayer(k)/HeatDiff_mulch(k,n)+thickPerLayer(k-1)
     &    /HeatDiff_mulch(k-1,n))
     &    /((thickPerLayer(k)+thickPerLayer(k-1))/200.0D0)         ! the heat conductance (J/m/day/K)   
         if(GrRe_mulch(1,n).ge.GrRe_Critical) then
          HeatCond_Mul_C_Fr=a_free
     &      *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &      *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &      *86400.0D0   
          HeatCond_Mul_C_Fo=0.0D0
          HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         else
          HeatCond_Mul_C_Fr=a_free
     &      *sqrt(abs(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n)))
     &      *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &      *86400.0D0          
          HeatCond_Mul_C_Fo=
     &      (karman**2.0D0)*(0.5D0*(u_mulch(k)+u_mulch(k-1)))
c     &      /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &       *0.155D0
     &      *(HeatCapa_Air(k-1,n)+HeatCapa_Air(k,n))/2.0D0
     &      *86400.0D0 
          HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
         endif
        endif
        g_Heat_Sensi(2*k-1,n)=(HeatCond_Mul_D+HeatCond_Mul_C)
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k-1,n))
       enddo
       g_Heat_Sensi(2*k,1)=0.0D0
       g_Heat_Sensi(2*k,SurNodeIndex_M)=0.0D0
       do n=2,SurNodeIndex_M-1
        g_Heat_Sensi(2*k,n)=
     &   (widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/
     &   (widthPerMulchUnit(n-1)/HeatDiff_mulch(k,n-1)
     &   +widthPerMulchUnit(n)/HeatDiff_mulch(k,n))
     &   *(MulchEleTmpr_temp(k,n)-MulchEleTmpr_temp(k,n-1))
     &   /((widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/200.0D0)
       enddo 
       endif
      enddo
      endif
      
      layer_need_fix=0
             
cccz -------------Establish Governing Equations (solve mulch temperature and moisture)------------------------------------- 
2410   if(PerOccupation.eq.2.0D0) then  ! the water level is high enough, this layer is now a water layer, fix all the values under water surface
       do n=1,SurNodeIndex_M-1
        do k=1,SubmergeIndex   
         MulchElehNew_temp2(k,n)=ThreHnew_UnderWater
         MulchElehNew_temp0(k,n)=ThreHnew_UnderWater
         MulchEleTmpr_temp2(k,n)=Tmpr_Sur(n)
         MulchEleTmpr_temp0(k,n)=Tmpr_Sur(n)
        enddo
       enddo
       do n=1,SurNodeIndex_M-1
        do k=SubmergeIndex+1,mulchLayer
         a11=WaterCapa_mulch(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
         a12=f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
         a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    +c_water_s*rho_w*WaterCapa_mulch(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
         a22=c_water_s*rho_w*MulchEleThNew_temp(k,n)                 ! partial T/partial t within mulch solid, i.e., liquid water
     &    +c_mulch_s*rho_mulch_b+                                    ! the wood part, the f_mulch_pore should already be adjusted, i.e., (1-f_mulch_pore)
     &    +(c_air_s*rho_dryair+c_vap_s*VaporAct_mulch_D(k,n))
     &    *f_mulch_pore
     &    +(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
         b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +RainFallInput_temp(k+1,n)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep
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
         SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)                     ! Cramer's rule for linear system
         SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
         MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
         MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
        enddo
       enddo
      else
       do k=1,SubmergeIndex-1
        do n=1,SurNodeIndex_M-1
         MulchElehNew_temp2(k,n)=ThreHnew_UnderWater
         MulchElehNew_temp0(k,n)=ThreHnew_UnderWater
         MulchEleTmpr_temp2(k,n)=Tmpr_Sur(n)
         MulchEleTmpr_temp0(k,n)=Tmpr_Sur(n)
        enddo
       enddo    
       k=SubmergeIndex
       thickLayer_Sur=thickPerLayer(k)*PerOccupation              ! adjust the pond surface element thickness
       do n=1,SurNodeIndex_M-1
        a11=WaterCapa_mulch(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    +c_water_s*rho_w*WaterCapa_mulch(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_water_s*rho_w*MulchEleThNew_temp(k,n)                  ! partial T/partial t within mulch solid, i.e., liquid water
     &    +c_mulch_s*rho_mulch_b+                                    ! the wood part, the f_mulch_pore should already be adjusted, i.e., (1-f_mulch_pore)
     &    +(c_air_s*rho_dryair+c_vap_s*VaporAct_mulch_D(k,n))
     &    *f_mulch_pore
     &    +(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickLayer_Sur/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +RainFallInput_temp(k+1,n)
     &    /(thickLayer_Sur/100.0D0)*LocalStep
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
        SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)           ! Cramer's rule for linear system               
        SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
       enddo 
       do k=SubmergeIndex+1,mulchLayer                     ! for other layer
        do n=1,SurNodeIndex_M-1
        a11=WaterCapa_mulch(k,n)*rho_w
     &    +f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))            ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
        a12=f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)     ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    -0.0002124D0*MulchElehNew_temp(k,n)
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        a21=(c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *(f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *0.02124D0/(MulchEleTmpr_temp(k,n)+273.15D0))              ! 0.02124D0\approx 9.81D0*0.018D0/8.314D0
     &    +c_water_s*rho_w*WaterCapa_mulch(k,n)
     &    *(MulchEleTmpr_temp(k,n)-Tmpr_ref)
        a22=c_water_s*rho_w*MulchEleThNew_temp(k,n)                ! partial T/partial t within mulch solid, i.e., liquid water
     &    +c_mulch_s*rho_mulch_b+                                  ! the wood part, the f_mulch_pore should already be adjusted, i.e., (1-f_mulch_pore)
     &    +(c_air_s*rho_dryair+c_vap_s*VaporAct_mulch_D(k,n))      ! partial T/partial t for the air space other than mulch solid
     &    *f_mulch_pore+
     &    (c_vap_v+c_vap_s*(MulchEleTmpr_temp(k,n)-Tmpr_ref))
     &    *f_mulch_pore*VaporAct_mulch_D(k,n)
     &    *(4975.9D0/((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0)
     &    -0.0002124D0*MulchElehNew_temp(k,n)                       ! 0.0002124D0\approx 9.81D0*0.018D0/8.314D0/100.0D0
     &    /((MulchEleTmpr_temp(k,n)+273.15D0)**2.0D0))
        b11=((g_Vapor(2*k+1,n)-g_Vapor(2*k-1,n))
     &    /(thickPerLayer(k)/100.0D0)+(g_Vapor(2*k,n+1)-g_Vapor(2*k,n))
     &    /(widthPerMulchUnit(n)/100.0D0))*LocalStep
     &    +RainFallInput_temp(k+1,n)
     &    /(thickPerLayer(k)/100.0D0)*LocalStep
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
        SolA=(a22*b11-a12*b12)/(a11*a22-a12*a21)           ! Cramer's rule for linear system               
        SolB=(a11*b12-a21*b11)/(a11*a22-a12*a21)
        MulchElehNew_temp2(k,n)=MulchElehNew_temp0(k,n)+SolA*100.0D0
        MulchEleTmpr_temp2(k,n)=MulchEleTmpr_temp0(k,n)+SolB
       enddo
      enddo
      endif
           
cccz -----------------Convergence Check (Picard Iteration)-----------------------------------------------
cccz_test
c      if(IterMulch.eq.1) then
c         aaaa=1 
c      endif
      
cccz reset the cumulative smoother
      if(IterMulch.eq.1) then
       do k=1,mulchLayer
        do n=1,SurNodeIndex_M-1
         MulchElehNew_temp3(k,n)=0.0D0
         MulchEleTmpr_temp3(k,n)=0.0D0  
        enddo
       enddo
      endif
cccz reset/initialize the error measurements
      maxhNew=0.0D0
      maxTmpr=0.0D0
      maxhNewDiff=0.0D0
      maxTmprDiff=0.0D0
      if(ForceForwards.eq.1) then
        if(PerOccupation.eq.2.0D0) then  ! the water level is high enough, this layer is now a water layer, fix all the values under water surface
          kklow=SubmergeIndex+1
        else
          kklow=SubmergeIndex
        endif
        meanMulchHnew=0.0D0
        do kk=kklow,layer_need_fix  
         meanMulchHnew=0.0D0
         do nn=1,SurNodeIndex_M-1
          meanMulchHnew=meanMulchHnew
     &      +MulchElehNew_temp2(kk,nn)*widthPerMulchUnit(nn)
         enddo
       meanMulchHnew=min(meanMulchHnew/totalMulchWidth,ThreHnew_HardFix)
         do nn=1,SurNodeIndex_M-1
          MulchElehNew_temp2(kk,nn)=meanMulchHnew
          MulchElehNew_temp3(kk,nn)=meanMulchHnew
         enddo
        enddo 
      endif
      
      if(localstep.le.minLocalStep) then           ! when the timestep is so small, jump out of the Picard iteration
         IterMulch=mulchLayer+1  
      endif

      if(IterMulch.le.mulchLayer) then             ! the Picard iteration has to be processed for times >= mulch layer, all mulch layer has to be updated
       do k=1,mulchLayer
        do n=1,SurNodeIndex_M-1                       ! cccz protection for 1. hnew>0; 2. hnew=NAN; 3. Tmpr=NAN
        ForceForwards=0
        layer_need_fix=0
        if(MulchElehNew_temp2(k,n).gt.0.0D0) then
         if(LocalStep.gt.minLocalStep) then
          goto 2407                              ! keep reducing the time step
         else
          if(PerOccupation.eq.2.0D0) then  ! the water level is high enough, this layer is now a water layer, fix all the values under water surface
           kklow=SubmergeIndex+1
          else
           kklow=SubmergeIndex
          endif  
          if(g_Vapor(2*k-1,n).le.0.0D0) then ! upwards flux cause positive head
           meanG_vapor=0.0D0
           meanG_latent_heat=0.0D0
           do nn=1,SurNodeIndex_M-1
            meanG_vapor=meanG_vapor
     &       +g_Vapor(2*kklow-1,nn)*widthPerMulchUnit(nn)
            meanG_latent_heat=meanG_latent_heat
     &       +g_Heat_Latent(2*kklow-1,nn)*widthPerMulchUnit(nn)
           enddo
           meanG_vapor=meanG_vapor/totalMulchWidth
           meanG_latent_heat=meanG_latent_heat/totalMulchWidth
           do kk=kklow,k                        ! set one layer of flux above
            do nn=1,SurNodeIndex_M-1 
             g_Vapor(2*kk+1,nn)=meanG_vapor
             g_Heat_Latent(2*kk+1,nn)=meanG_latent_heat
             g_Vapor(2*kk,nn)=0.0D0
             g_Heat_Latent(2*kk,nn)=0.0D0
            enddo
           enddo
          else                               ! downwards flux cause positive
           meanG_vapor=0.0D0
           meanG_latent_heat=0.0D0
           do nn=1,SurNodeIndex_M-1
            meanG_vapor=meanG_vapor
     &       +g_Vapor(2*k-1,nn)*widthPerMulchUnit(nn)
            meanG_latent_heat=meanG_latent_heat
     &       +g_Heat_Latent(2*k-1,nn)*widthPerMulchUnit(nn)
           enddo
           meanG_vapor=meanG_vapor/totalMulchWidth
           meanG_latent_heat=meanG_latent_heat/totalMulchWidth
           do nn=1,SurNodeIndex_M-1 
            g_Vapor(2*k+1,nn)=min(g_Vapor(2*k+1,nn),meanG_vapor)
            g_Heat_Latent(2*k+1,nn)=
     &         min(g_Heat_Latent(2*k+1,nn),meanG_vapor)
            g_Vapor(2*k,nn)=0.0D0
            g_Heat_Latent(2*k,nn)=0.0D0
           enddo  
          endif    
          ForceForwards=1
          layer_need_fix=k 
          goto 2410
         endif
        endif
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) then
            goto 2407
        endif
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) then
            goto 2407
        endif
        enddo
       enddo
       do n=1,SurNodeIndex_M-1
        do k=1,mulchLayer
c         MulchElehNew_temp(k,n)=MulchElehNew_temp2(k,n)
         MulchElehNew_temp(k,n)=
     &     min(MulchElehNew_temp2(k,n),ThreHnew_HardFix)     
         MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp2(k,n)
        enddo 
        do k=1,SubmergeIndex-1            !cccz water content calculation
          MulchEleThNew_temp(k,n)=f_mulch_pore     ! saturated part
        enddo
        k=SubmergeIndex
        if(PerOccupation.eq.2.0D0) then
          MulchEleThNew_temp(k,n)=f_mulch_pore   !cccz we avoid calculate the water content for water-filled layer
        else
          MulchEleThNew_temp(k,n)=
     &      WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &      LIGN_mass_Frac(k),lInput) 
        endif
        do k=SubmergeIndex+1,mulchLayer
          MulchEleThNew_temp(k,n)=
     &      WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &      LIGN_mass_Frac(k),lInput) 
        enddo
        enddo
        IterMulch=IterMulch+1
       goto 2001
      else
       IF(IterMulch.le.MaxIter) then
        do n=1,SurNodeIndex_M-1                       ! cccz protection for 1. hnew>0; 2. hnew=NAN; 3. Tmpr=NAN
        do k=1,mulchLayer
        ForceForwards=0
        layer_need_fix=0
        if(MulchElehNew_temp2(k,n).gt.0.0D0) then
         if(LocalStep.gt.minLocalStep) then
          goto 2407                              ! keep reducing the time step
         else
          if(PerOccupation.eq.2.0D0) then  ! the water level is high enough, this layer is now a water layer, fix all the values under water surface
           kklow=SubmergeIndex+1
          else
           kklow=SubmergeIndex
          endif  
          if(g_Vapor(2*k-1,n).le.0.0D0) then ! upwards flux cause positive head
           meanG_vapor=0.0D0
           meanG_latent_heat=0.0D0
           do nn=1,SurNodeIndex_M-1
            meanG_vapor=meanG_vapor
     &       +g_Vapor(2*kklow-1,nn)*widthPerMulchUnit(nn)
            meanG_latent_heat=meanG_latent_heat
     &       +g_Heat_Latent(2*kklow-1,nn)*widthPerMulchUnit(nn)
           enddo
           meanG_vapor=meanG_vapor/totalMulchWidth
           meanG_latent_heat=meanG_latent_heat/totalMulchWidth
           do kk=kklow,k                        ! set one layer of flux above
            do nn=1,SurNodeIndex_M-1 
             g_Vapor(2*kk+1,nn)=meanG_vapor
             g_Heat_Latent(2*kk+1,nn)=meanG_latent_heat
             g_Vapor(2*kk,nn)=0.0D0
             g_Heat_Latent(2*kk,nn)=0.0D0
            enddo
           enddo
          else                               ! downwards flux cause positive
           meanG_vapor=0.0D0
           meanG_latent_heat=0.0D0
           do nn=1,SurNodeIndex_M-1
            meanG_vapor=meanG_vapor
     &       +g_Vapor(2*k-1,nn)*widthPerMulchUnit(nn)
            meanG_latent_heat=meanG_latent_heat
     &       +g_Heat_Latent(2*k-1,nn)*widthPerMulchUnit(nn)
           enddo
           meanG_vapor=meanG_vapor/totalMulchWidth
           meanG_latent_heat=meanG_latent_heat/totalMulchWidth
           do nn=1,SurNodeIndex_M-1 
            g_Vapor(2*k+1,nn)=min(g_Vapor(2*k+1,nn),meanG_vapor)
            g_Heat_Latent(2*k+1,nn)=
     &         min(g_Heat_Latent(2*k+1,nn),meanG_vapor)
            g_Vapor(2*k,nn)=0.0D0
            g_Heat_Latent(2*k,nn)=0.0D0
           enddo  
          endif    
          ForceForwards=1
          layer_need_fix=k 
          goto 2410
         endif
        endif
        if(MulchElehNew_temp2(k,n).ne.MulchElehNew_temp2(k,n)) then
            goto 2407
        endif
        if(MulchEleTmpr_temp2(k,n).ne.MulchEleTmpr_temp2(k,n)) then
            goto 2407
        endif
        enddo
        enddo
        do n=1,SurNodeIndex_M-1 
        do k=1,mulchLayer
        MulchElehNew_temp3(k,n)=
     &    (dble(IterMulch-mulchLayer-1)*MulchElehNew_temp3(k,n)
     &     +MulchElehNew_temp2(k,n))/dble(IterMulch-mulchLayer)
cccz
        MulchElehNew_temp3(k,n)=
     &     min(MulchElehNew_temp3(k,n),ThreHnew_HardFix) 
cccz         
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
        do k=1,SubmergeIndex-1    !cccz water content calculation
          MulchEleThNew_temp(k,n)=f_mulch_pore
        enddo
        k=SubmergeIndex
        if(PerOccupation.eq.2.0D0) then
          MulchEleThNew_temp(k,n)=f_mulch_pore
        else
          MulchEleThNew_temp(k,n)=
     &      WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &      LIGN_mass_Frac(k),lInput)   
        endif
        do k=SubmergeIndex+1,mulchLayer
          MulchEleThNew_temp(k,n)=
     &      WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &      LIGN_mass_Frac(k),lInput) 
        enddo
        enddo
        IterMulch=IterMulch+1
cccz_try force forwards here to get rid of Picard Iteration when step=minLocalStep
        if(localstep.le.minLocalStep) then
          maxhNewDiff=0.0D0
          maxTmprDiff=0.0D0
        endif
        if(max(maxhNewDiff/maxhNew,maxTmprDiff/maxTmpr).gt.errTol) then
         goto 2001
        else
         TotalTime=TotalTime+LocalStep
         if(TotalTime-Step.gt.-0.01D0*dtMin) then
           if(TimeShrink.eq.0) then   ! if no 'local time step' was used, 'one step' was used for calculating water/energy exchange on soil surface
            goto 2402
           else                       ! if multiple 'local time step' was used, then 'cumulative averaged' water/energy exchange were calculated in 2003 and 2004
            goto 2403
           endif
         else
          do k=1,mulchLayer
           do n=1,SurNodeIndex_M-1
            MulchElehNew_temp0(k,n)=MulchElehNew_temp(k,n)
            MulchEleTmpr_temp0(k,n)=MulchEleTmpr_temp(k,n)
           enddo
          enddo
          goto 2403
2404      LocalStep=min(LocalStep*DMul1,Step-TotalTime)
          IterMulch=1
          goto 2001
         endif
        endif
       ELSE
2407    TimeShrink=1
         do n=1,SurNodeIndex_M-1
          do k=1,mulchLayer
          MulchElehNew_temp(k,n)=MulchElehNew_temp0(k,n)
          MulchEleTmpr_temp(k,n)=MulchEleTmpr_temp0(k,n)
          enddo
          do k=1,SubmergeIndex-1    !cccz water content calculation
            MulchEleThNew_temp(k,n)=f_mulch_pore
          enddo
          k=SubmergeIndex
          if(PerOccupation.eq.2.0D0) then
            MulchEleThNew_temp(k,n)=f_mulch_pore
          else
            MulchEleThNew_temp(k,n)=
     &       WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &       LIGN_mass_Frac(k),lInput)   
          endif
          do k=SubmergeIndex+1,mulchLayer
            MulchEleThNew_temp(k,n)=
     &        WQ_CERES_MULCH(MulchElehNew_temp(k,n),rho_mulch_b,
     &        LIGN_mass_Frac(k),lInput) 
          enddo
         enddo
         LocalStep=max(LocalStep*DMul2,0.0001D0*dtMin)
         IterMulch=1
         goto 2001
       ENDIF
      endif
      
cccz ----------------Update the water and energy fluxes across soil surface (Whole Step)----------------------------------
c    first do the vapor flux calculation
c    unit within this module will be g/day/m^2
c    unit taken in watermov module is g/day/cm^2
2402  if(PerOccupation.eq.2.0D0) then
          kk=2*SubmergeIndex+1
      else
          kk=2*SubmergeIndex-1
      endif
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)                                  !cccz update 'Varbw_Mulch' for 'SurWater', all the water will pass to the soil surface
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,nn)/10000.0D0                ! the evaporation is calculated on certain height
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)                          !cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)  
        elseif(n.eq.qLeft.and.nn.gt.1) then                           !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +RainFallInput_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(kk,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Vapor(kk,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)  
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,nn)/10000.0D0         
         Varbw_Mulch(kSurR,2)=-g_Vapor(kk,nn)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)
        else
         kSurL=SurfNodeSurfIndexH(n)                                  !cccz update 'Varbw_Mulch' for 'SurWater', all the water will pass to the soil surface
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,nn)/10000.0D0                ! the evaporation is calculated on certain height
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)                          !cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3) 
        endif      
       enddo  
      enddo
      

      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
          kSurL=SurfNodeSurfIndexH(n)
          VarBT(kSurL,1)=AirTemp_Wea
          VarBT(kSurL,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurL,3)=VarBT(kSurL,2)*AirTemp_Wea
          VarBT(kSurL,4)=0.0D0
          varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
          varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
          varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
          varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        elseif(n.eq.qLeft.and.nn.gt.1) then                        !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
          kSurL=SurfNodeSurfIndexH(n)
          VarBT(kSurL,1)=AirTemp_Wea
         VarBT(kSurL,2)=(g_Heat_CD_total(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Heat_CD_total(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
          VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
          VarBT(kSurL,4)=0.0D0       ! 8.64D0=86400D0/10000D0
          varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
          varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
          varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
          varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
          kSurR=SurfNodeSurfIndexH(n)
          VarBT(kSurR,1)=AirTemp_Wea
          VarBT(kSurR,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(kk,nn)
          VarBT(kSurR,4)=0.0D0                            ! 8.64D0=86400D0/10000D0
          varBT_Mulch(kSurR,1)=VarBT(kSurR,1)
          varBT_Mulch(kSurR,2)=VarBT(kSurR,2)
          varBT_Mulch(kSurR,3)=VarBT(kSurR,3)
          varBT_Mulch(kSurR,4)=VarBT(kSurR,4)
        else
          kSurL=SurfNodeSurfIndexH(n)
          VarBT(kSurL,1)=AirTemp_Wea
          VarBT(kSurL,2)=g_Heat_CD_total(1,nn)
          VarBT(kSurL,3)=VarBT(kSurL,2)*AirTemp_Wea
          VarBT(kSurL,4)=0.0D0
          varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
          varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
          varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
          varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        endif      
       enddo  
      enddo 
      
cccz old code      
c      if(PerOccupation.eq.2.0D0) then
c       kk=SubmergeIndex+1
c      else
c       kk=SubmergeIndex
c      endif
c      do n=1,SurNodeIndex-1
c       if(n.eq.1) then
c        kSurL=SurfNodeSurfIndexH(n)
c        VarBT(kSurL,1)=AirTemp_Wea
c        VarBT(kSurL,2)=g_Heat_CD_total(1,n)
c        VarBT(kSurL,3)=VarBT(kSurL,2)*AirTemp_Wea
c        VarBT(kSurL,4)=0.0D0
c        varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
c        varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
c        varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
c        varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
c       elseif(n.eq.(SurNodeIndex-1)) then
c        kSurL=SurfNodeSurfIndexH(n)
c        kSurR=SurfNodeSurfIndexH(n+1)
cc        VarBT(kSurL,1)=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
cc     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
cc        VarBT(kSurL,2)=(g_Heat_CD_total(1,n-1)*widthPerMulchUnit(n-1)
cc     &    +g_Heat_CD_total(1,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
cc        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
cc        VarBT(kSurL,4)=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
cc     &    +netRad(kk,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
c        VarBT(kSurL,1)=AirTemp_Wea
c        VarBT(kSurL,2)=(g_Heat_CD_total(1,n-1)*widthPerMulchUnit(n-1)
c     &    +g_Heat_CD_total(1,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
c        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
c        VarBT(kSurL,4)=0.0D0
cc        VarBT(kSurR,1)=MulchEleTmpr_temp(kk,n)
cc        VarBT(kSurR,2)=g_Heat_CD_total(1,n)
cc        VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(kk,n)
cc        VarBT(kSurR,4)=netRad(kk,n)*8.64D0                            ! 8.64D0=86400D0/10000D0
c        VarBT(kSurR,1)=AirTemp_Wea
c        VarBT(kSurR,2)=g_Heat_CD_total(1,n)
c        VarBT(kSurR,3)=VarBT(kSurR,2)*MulchEleTmpr_temp(kk,n)
c        VarBT(kSurR,4)=0.0D0                            ! 8.64D0=86400D0/10000D0
c        varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
c        varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
c        varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
c        varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
c        varBT_Mulch(kSurR,1)=VarBT(kSurR,1)
c        varBT_Mulch(kSurR,2)=VarBT(kSurR,2)
c        varBT_Mulch(kSurR,3)=VarBT(kSurR,3)
c        varBT_Mulch(kSurR,4)=VarBT(kSurR,4)
c       else
c        kSurL=SurfNodeSurfIndexH(n)
cc        VarBT(kSurL,1)=(MulchEleTmpr_temp(kk,n-1)*widthPerMulchUnit(n-1)
cc     &    +MulchEleTmpr_temp(kk,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
cc        VarBT(kSurL,2)=(g_Heat_CD_total(1,n-1)*widthPerMulchUnit(n-1)
cc     &    +g_Heat_CD_total(1,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
cc        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
cc        VarBT(kSurL,4)=(netRad(kk,n-1)*widthPerMulchUnit(n-1)
cc     &    +netRad(kk,n)*widthPerMulchUnit(n))
cc     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))*8.64D0       ! 8.64D0=86400D0/10000D0
c        VarBT(kSurL,1)=AirTemp_Wea
c        VarBT(kSurL,2)=(g_Heat_CD_total(1,n-1)*widthPerMulchUnit(n-1)
c     &    +g_Heat_CD_total(1,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))
c        VarBT(kSurL,3)=VarBT(kSurL,2)*VarBT(kSurL,1)
c        VarBT(kSurL,4)=0.0D0       ! 8.64D0=86400D0/10000D0
c        varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
c        varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
c        varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
c        varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
c       endif
c      enddo
      
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput_temp(n,k)
       enddo      
      enddo
cccz old code      
c      do k=1,SurNodeIndex-1
c       do n=1,mulchLayer+1
c        RainFallInput(n,k)=RainFallInput_temp(n,k)
c       enddo
c      enddo
      goto 2406           ! to the final assignment
      
cccz ------------------------Update the water and energy fluxes across soil surface (For Each Time Segment)--------------------------      
c    for each time segment during the iteration, use the same scheme showed following 'Index 2402'
2403  if(PerOccupation.eq.2.0D0) then
        kk=2*SubmergeIndex+1
      else
        kk=2*SubmergeIndex-1
      endif
      
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then   
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(kk,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
        elseif(n.eq.qLeft.and.nn.gt.1) then                           !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
        Local_VarBW1=(RainFallInput_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +RainFallInput_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Local_VarBW2=-(g_Vapor(kk,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Vapor(kk,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Local_VarBW3=Local_VarBW2-VarBW(kSurL,1)
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep  
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(kk,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurR)*Local_VarBW3)*LocalStep
        else
         kSurL=SurfNodeSurfIndexH(n)
         Local_VarBW1=RainFallInput_temp(1,nn)/10000.0D0
         Local_VarBW2=-g_Vapor(kk,nn)/10000.0D0
         Local_VarBW3=Local_VarBW2-Local_VarBW1
         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
        endif      
       enddo  
      enddo
      
cccz old code
c      do n=1,SurNodeIndex-1
c       if(n.eq.1) then
c         kSurL=SurfNodeSurfIndexH(n)
c         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
c         Local_VarBW2=-g_Vapor(kk,n)/10000.0D0
c         Local_VarBW3=Local_VarBW2-Local_VarBW1
c         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
c         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
c         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
c         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
c       elseif(n.eq.(SurNodeIndex-1)) then
c         kSurL=SurfNodeSurfIndexH(n)
c         kSurR=SurfNodeSurfIndexH(n+1)
c         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
c     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
c         Local_VarBW2=-(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
c     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
c         Local_VarBW3=Local_VarBW2-Local_VarBW1
c         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
c         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
c         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
c         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
c         Local_VarBW1=RainFallInput_temp(1,n)/10000.0D0
c         Local_VarBW2=-g_Vapor(kk,n)/10000.0D0
c         Local_VarBW3=Local_VarBW2-Local_VarBW1
c         VarBW1_temp(n+1)=VarBW1_temp(n+1)+Local_VarBW1*LocalStep
c         VarBW2_temp(n+1)=VarBW2_temp(n+1)+Local_VarBW2*LocalStep
c         VarBW3_temp(n+1)=VarBW3_temp(n+1)+Local_VarBW3*LocalStep
c         Q_temp(n+1)=Q_temp(n+1)+(-Width(kSurR)*Local_VarBW3)*LocalStep
c       else
c         kSurL=SurfNodeSurfIndexH(n)
c         Local_VarBW1=(RainFallInput_temp(1,n-1)*widthPerMulchUnit(n-1)
c     &    +RainFallInput_temp(1,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
c         Local_VarBW2=-(g_Vapor(kk,n-1)*widthPerMulchUnit(n-1)
c     &    +g_Vapor(kk,n)*widthPerMulchUnit(n))
c     &    /(widthPerMulchUnit(n-1)+widthPerMulchUnit(n))/10000.0D0
c         Local_VarBW3=Local_VarBW2-VarBW(kSurL,1)
c         VarBW1_temp(n)=VarBW1_temp(n)+Local_VarBW1*LocalStep
c         VarBW2_temp(n)=VarBW2_temp(n)+Local_VarBW2*LocalStep
c         VarBW3_temp(n)=VarBW3_temp(n)+Local_VarBW3*LocalStep
c         Q_temp(n)=Q_temp(n)+(-Width(kSurL)*Local_VarBW3)*LocalStep
c       endif 
c      enddo
      
    
      if(PerOccupation.eq.2.0D0) then
       kk=SubmergeIndex+1
      else
       kk=SubmergeIndex
      endif
      
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         Local_VarBT1=AirTemp_Wea
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=0.0D0  
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        elseif(n.eq.qLeft.and.nn.gt.1) then                        !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         Local_VarBT1=AirTemp_Wea
         Local_VarBT2=(g_Heat_CD_total(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Heat_CD_total(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=0.0D0                                            ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         Local_VarBT1=AirTemp_Wea
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=0.0D0                                            ! 8.64D0=86400D0/10000D0
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep    
        else
         Local_VarBT1=AirTemp_Wea
         Local_VarBT2=g_Heat_CD_total(1,nn)
         Local_VarBT3=Local_VarBT2*Local_VarBT1
         Local_VarBT4=0.0D0  
         VarBT1_temp(n)=VarBT1_temp(n)+Local_VarBT1*LocalStep
         VarBT2_temp(n)=VarBT2_temp(n)+Local_VarBT2*LocalStep
         VarBT3_temp(n)=VarBT3_temp(n)+Local_VarBT3*LocalStep
         VarBT4_temp(n)=VarBT4_temp(n)+Local_VarBT4*LocalStep
        endif      
       enddo  
      enddo 
      
      
  
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)
     &   +RainFallInput_temp(n,k)*LocalStep
       enddo
      enddo     
cccz old code      
c      do k=1,SurNodeIndex-1
c       do n=1,mulchLayer+1
c        RainFallInput(n,k)=RainFallInput(n,k)
c     &   +RainFallInput_temp(n,k)*LocalStep
c       enddo
c      enddo
      if(TotalTime-Step.gt.-0.01D0*dtMin) then
          goto 2405                               ! finalize the iteration
      else
          goto 2404                               ! process to the next local step
      endif
cccz ---------'cumulative averaged' water/energy exchange between 2403 and 2404------------------------
c    because that is for each time segment, we have to put them together.
c    determine the vapor and energy exchange on soil surface
2405  do n=1,SurNodeIndex
       kSur=SurfNodeSurfIndexH(n)
       Varbw_Mulch(kSur,1)=VarBW1_temp(n)/TotalTime   !cccz water vapor part
       Varbw_Mulch(kSur,2)=VarBW2_temp(n)/TotalTime
       Varbw_Mulch(kSur,3)=VarBW3_temp(n)/TotalTime
       VarBW(kSur,1)=Varbw_Mulch(kSur,1)
       VarBW(kSur,2)=Varbw_Mulch(kSur,2)
       VarBW(kSur,3)=Varbw_Mulch(kSur,3)
       nNode=KXB(kSur)
       Q(nNode)=Q_temp(n)/TotalTime
       VarBT(kSur,1)=VarBT1_temp(n)/TotalTime         !cccz energy part
       VarBT(kSur,2)=VarBT2_temp(n)/TotalTime
       VarBT(kSur,3)=VarBT3_temp(n)/TotalTime
       VarBT(kSur,4)=VarBT4_temp(n)/TotalTime
       varBT_Mulch(kSur,1)=VarBT(kSur,1)
       varBT_Mulch(kSur,2)=VarBT(kSur,2)
       varBT_Mulch(kSur,3)=VarBT(kSur,3)
       varBT_Mulch(kSur,4)=VarBT(kSur,4)
      enddo
      do k=1,SurNodeIndex_M-1
       do n=1,mulchLayer+1
        RainFallInput(n,k)=RainFallInput(n,k)/TotalTime
       enddo
      enddo
      goto 2406
      
2406  aaaa=2          !cccz final assignment
      do n=1,SurNodeIndex_M-1
       do k=1,mulchLayer
         MulchElehNew(k,n)=MulchElehNew_temp(k,n)
         MulchEleTmpr(k,n)=MulchEleTmpr_temp(k,n)
         MulchEleThNew(k,n)=MulchEleThNew_temp(k,n)
       enddo
      enddo
      goto 2600

cccz --------------------------------TOTALLY PONDED WATER (25XX Labels)-----------------------------------------------------------
c    the mulch is totally submerged in surface water (which should be unusual), in this case, we assume 
c       (1) the soil surfae temp and ponded water temp are the same
c       (2) the mulch is totally saturated
c       (3) we only need to adjust the varBW2 (evaporation) for all varBW and varBT quantities
c       (4) no Picard iteration is needed.
2500  do n=1,SurNodeIndex_M-1                   !cccz assign water content, temperature, water-head directly
c       qLeft=SurfNodeNodeIndexH(n)
c       qRight=SurfNodeNodeIndexH(n+1)
c       Tmpr_Sur=AirTemp_Wea
       do k=1,mulchLayer 
         MulchElehNew(k,n)=ThreHnew_UnderWater             ! this should be 0, however, we make it with a small value
         MulchElehNew_temp(k,n)=ThreHnew_UnderWater
         MulchEleTmpr(k,n)=Tmpr_Sur(n)
         MulchEleTmpr_temp(k,n)=Tmpr_Sur(n)
         MulchEleThNew(k,n)=f_mulch_pore
         MulchEleThNew_temp(k,n)=f_mulch_pore
       enddo
      enddo
      
      do n=1,SurNodeIndex_M-1                   !cccz make some flux adjustments
       do k=1,mulchLayer 
        RelaHumid_mulch(k,n)=1.0D0                   !'vapor conditions'
        VaporSat_mulch(k,n)=1.0D06                   ! water density (g/m^3)           
        VaporAct_mulch_D(k,n)=1.0D06
        WaterDifici_mulch(k,n)=0.0D0     !cccz water difficit in each element
        VaporDiff_mulch(k,n)=0.0D0       ! Water Vapor diffusivity (m^2/day)
        WaterCapa_mulch(k,n)=0.0D0
        HeatCapa_Air(k,n)=0.0D0          ! there is no air for heat capacity (J/m^3/K)
        HeatDiff_mulch(k,n)=0.0D0        ! there is no air for heat diffusivity (J/m/K/day)
       enddo
      enddo

      do k=1,mulchLayer                       !cccz assign 'dummy' water fluxes within the mulch grid
       do n=1,SurNodeIndex_M
        g_vapor(2*k-1,n)=0.0D0
        g_vapor(2*k,n)=0.0D0
       enddo
      enddo
      k=mulchLayer+1                          !cccz for water surface                           
      do n=1,SurNodeIndex_M-1
        !'potential evaporation' method
        SVPA_Sur=0.61D0*EXP((17.27D0*AirTemp_wea)/(AirTemp_wea+237.3D0))
        DEL_Sur=(0.61D0*EXP((17.27D0*(AirTemp_wea+1.0D0))
     &   /(AirTemp_wea+1.0D0+237.3D0)))-SVPA_Sur
        VPD_Sur=SVPA_Sur*(1.0D0-RelaHumi_Wea)          ! calculate VPD for the first mulching layer
        D31=0.622D0*(SVPA_Sur/(101.3D0-SVPA_Sur))      ! we use actural vapor pressure for D31 here, should be the saturated vapor pressure of wet bulb temp in air
        D32=2500.8D0-2.37D0*AirTemp_wea
        GAMMA_Sur=0.62D0*(1.006D0+(1.846D0*D31))
     &    /((0.622D0+D31)*(0.662D0+D31)*D32)*101.3D0                 ! since we are doing local analysis, suppose 'expose fraction PSh=1'
        g_Vapor(2*k-1,n)=-((DEL_Sur/GAMMA_Sur*max(netRad(k,n),0.0D0)
     &    *3600.0D0/(2500.8D0-(2.3668D0*AirTemp_wea)))
     &    +(VPD_Sur*109.375D0*(1.0D0+(0.149D0*u_above))))
     &    /((DEL_Sur/GAMMA_Sur)+1.0D0)*24.0D0                        ! vertical vapor flux near soil-mulch interface (g/m2/day) 
      enddo
     
      do n=1,SurNodeIndex_M-1                           ! cccz calculate the rainfall redistribution, everything go to soil
       RainFallInput_temp(1,n)=RainFallInput_temp(mulchLayer+1,n)
       do k=2,mulchLayer+1
         RainFallInput_temp(k,n)=0.0D0
       enddo
      enddo
           
      do k=1,mulchLayer                           ! cccz assign the latent heat flux
       do n=1,SurNodeIndex_M
        g_Heat_Latent(2*k-1,n)=0.0D0
        g_Heat_Latent(2*k,n)=0.0D0
       enddo
      enddo
      k=mulchLayer+1
      do n=1,SurNodeIndex_M-1                                         !cccz for mulch-air interface layer
       if(g_Vapor(2*k-1,n).ge.0.0D0) then                           ! downwards vapor flow
         g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &    *(AirTemp_Wea-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
       else                                                         ! upwards vapor flow
         g_Heat_Latent(2*k-1,n)=c_vap_s*g_Vapor(2*k-1,n)
     &    *(MulchEleTmpr_temp(k-1,n)-Tmpr_ref)+c_vap_v*g_Vapor(2*k-1,n)
       endif
      enddo

      do k=1,mulchLayer                           !cccz assign the sensible (diff) heat flux
       do n=1,SurNodeIndex_M
        g_Heat_Sensi(2*k-1,n)=0.0D0
        g_Heat_Sensi(2*k,n)=0.0D0
       enddo
      enddo
      k=mulchLayer+1      
      do n=1,SurNodeIndex_M-1
       Ra_mulch(1,n)=9.81D0*abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)   ! Rayleigh num
     &   *(((thickPerLayer(k-1))/100.0D0)**3.0D0)/nu_air
     &   /((546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)/2.0D0)/Dhm
       GrRe_mulch(1,n)=2.0D0*9.81D0*((z_wind-mulchThick)/100.0D0)       ! Richardson num = Grashof num / Reynolds num ^2, we use the absolute value here
     &    *abs(MulchEleTmpr_temp(k-1,n)-AirTemp_Wea)
     &    /(546.30D0+MulchEleTmpr_temp(k-1,n)+AirTemp_Wea)
     &    /((AirWind_Wea-u_above)**2.0D0)     
       qLeft=SurfNodeSurfIndexH(n)
       qRight=SurfNodeSurfIndexH(n+1)
       varbt_Temp=0.5D0*(varbt_Air(qLeft,2)+varbt_Air(qRight,2))
       if(Ra_mulch(1,n).le.Ra_Critical.or.DiffusionRes.eq.1) then
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
          HeatCond_Mul_C=0.0D0
       else
          HeatCond_Mul_D=sqrt(2.14D-5*HeatCapa_Air(k-1,n)*86400.0D0
     &      *HeatDiff_mulch(k-1,n))
     &      /(thickPerLayer(k-1)/100.0D0)                         ! the heat conductance (J/m/day/K)
        if(GrRe_mulch(1,n).ge.GrRe_Critical) then
          HeatCond_Mul_C_Fr=a_free
     &     *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &     *HeatCapa_Air(k-1,n)
     &     *86400D0                                               ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
          HeatCond_Mul_C_Fo=0.0D0
          HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
        else
          HeatCond_Mul_C_Fr=a_free
     &     *sqrt(abs(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)))
     &     *HeatCapa_Air(k-1,n)
     &     *86400D0                                              ! change unit from "W" to "J,day"  (i.e., J/m^3/day)
          HeatCond_Mul_C_Fo=
     &     (karman**2.0D0)*u_mulch(k)
c     &     /(log(LayerHeight(k)/(0.079D0*LayerHeight(k)))**2.0D0)
     &     *0.155D0
     &     *HeatCapa_Air(k-1,n)
     &     *86400D0                                              ! change unit from  "W" to "J,day"  (i.e., J/m^3/day)
          HeatCond_Mul_C=HeatCond_Mul_C_Fr+HeatCond_Mul_C_Fo
        endif
       endif
       g_Heat_Sensi(2*k-1,n)=
     &    100.0D0*sqrt(varbt_Temp*(HeatCond_Mul_D+HeatCond_Mul_C))
     &    *(AirTemp_Wea-MulchEleTmpr_temp(k-1,n)) 
       g_Heat_CD_total(1,n)=sqrt(varbt_Temp                    !cccz the surface heat conductance to assign new values for "heatmow"
     &    *(HeatCond_Mul_C+HeatCond_Mul_D))                    ! @ this step, the unit is J/m^2/day/K
     &    /100.00                                              ! change unit to J/cm^2/day/K   
      enddo
    
cccz ----------------Update the water and energy fluxes across soil surface (Whole Step Only for this case)----------------------------------
c    first do the vapor flux calculation
c    unit within this module will be g/day/m^2
c    unit taken in watermov module is g/day/cm^2
    
      kk=2*mulchLayer+1
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)          !cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,nn)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)  !cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3) 
        elseif(n.eq.qLeft.and.nn.gt.1) then                           !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurL,1)=
     &    (RainFallInput_temp(1,nn-1)*widthPerMulchUnit(nn-1)
     &    +RainFallInput_temp(1,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,2)=
     &    -(g_Vapor(kk,nn-1)*widthPerMulchUnit(nn-1)
     &    +g_Vapor(kk,nn)*widthPerMulchUnit(nn))
     &    /(widthPerMulchUnit(nn-1)+widthPerMulchUnit(nn))/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)          
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         Varbw_Mulch(kSurR,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurR,2)=-g_Vapor(kk,nn)/10000.0D0
         Varbw_Mulch(kSurR,3)=Varbw_Mulch(kSurR,2)-Varbw_Mulch(kSurR,1)
         VarBW(kSurR,1)=Varbw_Mulch(kSurR,1)
         VarBW(kSurR,2)=Varbw_Mulch(kSurR,2)
         VarBW(kSurR,3)=Varbw_Mulch(kSurR,3)
         nNode=KXB(kSurR)
         Q(nNode)=-Width(kSurR)*VarBW(kSurR,3)         
        else
         kSurL=SurfNodeSurfIndexH(n)          !cccz update 'Varbw_Mulch' for 'SurWater'
         Varbw_Mulch(kSurL,1)=RainFallInput_temp(1,nn)/10000.0D0
         Varbw_Mulch(kSurL,2)=-g_Vapor(kk,nn)/10000.0D0
         Varbw_Mulch(kSurL,3)=Varbw_Mulch(kSurL,2)-Varbw_Mulch(kSurL,1)
         VarBW(kSurL,1)=Varbw_Mulch(kSurL,1)  !cccz update 'VarBW' for 'WaterMov'
         VarBW(kSurL,2)=Varbw_Mulch(kSurL,2)
         VarBW(kSurL,3)=Varbw_Mulch(kSurL,3)
         nNode=KXB(kSurL)
         Q(nNode)=-Width(kSurL)*VarBW(kSurL,3)          
        endif      
       enddo  
      enddo
      
cccz we need to adjust the surface varbt variables, the varbt_air will keep the original values.
      
      kk=2*mulchLayer+1
      do nn=1,SurNodeIndex_M-1
       qLeft=SurfMulchNodeSubIndex(nn)
       qRight=SurfMulchNodeSubIndex(nn+1)
       do n=qLeft,qRight
        if(n.eq.1) then
         kSurL=SurfNodeSurfIndexH(n)
         VarBT(kSurL,1)=AirTemp_Wea
         VarBT(kSurL,2)=0.0D0
         VarBT(kSurL,3)=0.0D0
         VarBT(kSurL,4)=0.0D0
         varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
         varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
         varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
         varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        elseif(n.eq.qLeft.and.nn.gt.1) then                        !cccz "n.eq.qLeft.and.nn.gt.1" and "n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)" are mutually identical
         kSurL=SurfNodeSurfIndexH(n)
         VarBT(kSurL,1)=AirTemp_Wea
         VarBT(kSurL,2)=0.0D0
         VarBT(kSurL,3)=0.0D0
         VarBT(kSurL,4)=0.0D0       ! 8.64D0=86400D0/10000D0
         varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
         varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
         varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
         varBT_Mulch(kSurL,4)=VarBT(kSurL,4)
        elseif(n.eq.qRight.and.nn.lt.(SurNodeIndex_M-1)) then     
        elseif(n.eq.SurNodeIndex) then
         kSurR=SurfNodeSurfIndexH(n)
         VarBT(kSurR,1)=AirTemp_Wea
         VarBT(kSurR,2)=0.0D0
         VarBT(kSurR,3)=0.0D0
         VarBT(kSurR,4)=0.0D0
         varBT_Mulch(kSurR,1)=VarBT(kSurR,1)
         varBT_Mulch(kSurR,2)=VarBT(kSurR,2)
         varBT_Mulch(kSurR,3)=VarBT(kSurR,3)
         varBT_Mulch(kSurR,4)=VarBT(kSurR,4)
        else
         kSurL=SurfNodeSurfIndexH(n)
         VarBT(kSurL,1)=AirTemp_Wea
         VarBT(kSurL,2)=0.0D0
         VarBT(kSurL,3)=0.0D0
         VarBT(kSurL,4)=0.0D0
         varBT_Mulch(kSurL,1)=VarBT(kSurL,1)
         varBT_Mulch(kSurL,2)=VarBT(kSurL,2)
         varBT_Mulch(kSurL,3)=VarBT(kSurL,3)
         varBT_Mulch(kSurL,4)=VarBT(kSurL,4) 
        endif      
       enddo  
      enddo 
      goto 2600

cccz --------------------------------Final Assignment-----------------------------------------------------------
c    convert element based values to node based values
c    finish the calculation
2600    do n=1,numMulchNode
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
        
        
cccz ------------------------ output for paper writting-----------------
c      if(abs(time-Tentative_output).lt.1.0D-7
c     &      .or.time.gt.Tentative_output) then
c        Tentative_time_output=1.0D0/24.0D0/60.0D0
c        Tentative_output=time+Tentative_time_output
c        iday=int(time)
c        call caldat(iday,mm,id,iyyy) 
cc        write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy
c        TotalRa=0.0D0
c        TotalRi=0.0D0                                     
c        do n=1,SurNodeIndex-1
c          TotalRa=TotalRa+Ra_mulch(1,n)
c          TotalRi=TotalRi+GrRe_mulch(1,n)
c        enddo 
c        TotalRa=TotalRa/(SurNodeIndex-1)
c        TotalRi=TotalRi/(SurNodeIndex-1)
c        Runoff_mulch_record=abs(RunoffRight+RunoffRight_Efflux)
c     &    +abs(RunoffLeft+RunoffLeft_Efflux)
c        Runoff_mulch_record=Runoff_mulch_record/totalMulchWidth
c        Record_Infiltration_mulch=0.0D0
c        do k=1,munBP
c          Record_Infiltration_mulch=
c     &     Record_Infiltration_mulch+max(0.0,QAct(n))
c        enddo
c        Record_Infiltration_mulch=
c     &     Record_Infiltration_mulch/totalMulchWidth
c        Write(9999,9902) Time, Varbw_Air(1,1),Runoff_mulch_record,
c     &     Runoff_mulch_record,Record_Infiltration_mulch,
c     &     TotalRa,TotalRi
c      endif     
cccz -----------------------------------------------------------------------
      
      return
2106  Write(*,*) 'Mulch file not found'      
2100  Call errmes(im,il)
9901  Format (1x,A12,T20,A12,T49,A12,T64,A12,T76,A12,T95,A12,T109,A12)
9902  Format (1x,F17.6,',',5(F16.4,','),F16.4)    
      Return  
      END 