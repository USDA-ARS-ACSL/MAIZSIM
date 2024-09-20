  
      
cccz ************************************************************************************************************
c The development history
c There are three parts in the program, representing three trials we made for the runoff simulation
c Part I: Assume all of the runoff move to the next computation cell instantaneously. The Mass balance will work.
c Part II: Solve the fully coupled Saint-Venant, with some flat slope, the "h" will increase to a unreasonable value, so the diff
c         equation system may not always work
c Part III: Not only consider the coupled Saint-Venant, but also set a threshold for "h", the rest will go to the runoff.
c         The momentun equation will be used to estimate the flux speed.???
c ************************************************************************************************************

      Subroutine Surface_Water_Balance_Adjustment()
      include 'public.ins'
      include 'puweath.ins'
cccz      include 'puplant.ins'
      include 'PuSurface.ins'
      
      
      integer FlatSurface
      double precision totalLength,MaxPondFlatSur,
     !       AvePondFlatSur,MaxFlatSurHnew
      double precision  n_Stiff, g_accel, Add_Flux,
     !       RunoffLeft_Efflux_Old,RunoffRight_Efflux_Old,
     !       Slope_0_R, Slope_f, Slope_f_left, Slope_f_right,
     !       Slope_0_n,Slope_0_n_1,Fraction_n,Fraction_n_1,
     !       q_flux_square, h_flux_square,                        ! this are aux variables 
     !       q_flux_square1, h_flux_square1,
     !       q_flux_square2, h_flux_square2 
cccz SURNODE to index the node order from the whole node array      
cccz SURNODE_Sur to index the node order from the boundary node array
cccz Now use SurfNodeNodeIndexH/SurfNodeSurfIndexH to be consistent with the mulch module
      integer  FluxDir(15),iteration_Flux, iteration_Head, 
     !       FurtherCheck_h_Pond(NumBPD),Check_Exchange_V(NumBPD),
     !       InElement,EdgeReloc
      double precision iteration_Num,FluxLimit, HeadLimit
      double precision time_runoff_old
      double precision ll(NumBPD),H_L,H_M,H_R,L_L,L_M,L_R,Mean_Height_L,
     !       Mean_Height_R,Mean_Height,Ava_Height,Ava_Height_L,
     !       Ava_Height_R,Volume_L,Volume_R,Volume_T,Volume,
     !       Exchange_V(NumBPD), Exchange_V_temp(NumBPD)
      double precision RunoffInten(NumBPD),
     !    h_Pond_temp(NumBPD),q_Flux_temp(NumBPD),q_Flux_Node(NumBPD),
     !    h_Pond_Old(NumBPD),q_Flux_Old(NumBPD),q_Flux_Node_Old(NumBPD),
     !      Efflux_Old(NumBPD),
     !      RunoffLeft_old, RunoffRight_old,
     !      RunoffLeft_temp, RunoffRight_temp
cccz Zhuangji made the three variables to see if no runoff & surface is dry
cccz if so we can bypass this module to save some time.
      double precision TotalhStay_test,TotalqFlux_test,TotalR_test
      double precision CriticalH_R_Sur(NumNPD),
     !      BaseLeft,BaseRight,xBaseLeft,xBaseRight
      double precision varbw_Runoff_Loc_Record(NumBPD,3),   
     !     RunoffRight_Runoff_Loc_Record,
     !     RunoffLeft_Runoff_Loc_Record,
     !     Q_Runoff_Loc_Record(NumNPD) 
      integer TmprBCRecord(NumBPD),
     !     GasBCRecord(NumBPD)
      Common /SurWaterBalance/ FlatSurface,TmprBCRecord,GasBCRecord,
     !       ll,RunoffInten,
     !       h_Pond_temp,h_Pond_Old,
     !       q_Flux_temp,q_Flux_Node,q_Flux_Node_Old,q_Flux_Old,
     !       Efflux_Old, 
     !       RunoffLeft_old, RunoffRight_old,
     !       CriticalH_R_Sur,
     !       n_Stiff, g_accel,
     !       time_runoff_old,
     !       FluxLimit, HeadLimit,
     !       varbw_Runoff_Loc_Record,   
     !       RunoffRight_Runoff_Loc_Record,
     !       RunoffLeft_Runoff_Loc_Record,
     !       Q_Runoff_Loc_Record
      
      If (lInput.eq.1) then 
cccz Flat surface is 0, irregular surface is 1.
        FlatSurface=0
        time_runoff_old=time
        n_Stiff=0.5D0
        g_accel=0.0113425926D0
        g_accel=9.81D0*100.0D0/3600.0D0/3600.0D0/24.0D0/24.0D0
c -------------------- This part is an extension of weather module -----------------
c the basic usage of this portion is to provide an initialization
        
          
cccz CriticalH_R is the assumed surface height of ponded water before free flow occur
        CriticalH_R=0.1D0
cccz_this will be a new adjustment for runoff (mulch can hold more surface water)
      if (residueApplied.le.0) then 
         CriticalH_R=0.1D0  
      else
         CriticalH_R=mulchThick*WaterStorageFrac
      endif  
cccz_try
c        CriticalH_R=5.0D0
          
cccz CriticalH is for soil water module for 
        !CriticalH=5.101D0
          
cccz We treat runoff is a submodule of weather, but iterated with watermov
cccz but we do not need a module number because this runoff process is always consistent with water/heat_mov
cccz    ModNum=1
          
cccz Initially, we do not request the weather update for each hour based on "SurWeatherUpdate"
cccz But we actually do the weather update within the weather module, so remove such variable
cccz Because the weather module make the first calculation during the initialization

          
       do i=1,NumNPD
         RO(i)=0.0D0                 ! cccz: set the initial runoff from each soil node is zero
       enddo
          
       do i=1,NumBPD
         h_Pond(i)=0.0D0             ! cccz: initially, there is no ponded water (surface runoff water)
         q_Flux(i)=0.0D0             ! cccz: initially, there is surface water flux                
         h_Pond_temp(i)=0.0D0        ! cccz: zero ponded height during iteration
         q_Flux_temp(i)=0.0D0        ! cccz: zero surface flux during iteration
         q_Flux_Node(i)=0.0D0        ! cccz: initialize zero surface flux based on each node
         Efflux(i)=0.0D0             ! cccz: when h_Pond>CriticalH_R_Sur, then there will be free water flux
         RunoffInten(i)=0.0D0        ! cccz: the runoff intensity, RunoffInten=RO/width
         CriticalH_R_Sur(i)=CriticalH_R  ! cccz: Surface ponded height limit before free water flux
         ll(i)=0.0D0
       enddo
       RunoffLeft_Efflux=0.0D0         ! cccz: the runoff discharge from the left of the grid, Efflux portion
       RunoffRight_Efflux=0.0D0        ! cccz: the runoff discharge from the right of the grid, Efflux portion
          
c          We use a step-by-step method to read the grid on the surface.
c          As well as determine the horizontal distance, vertical distance between each node.
c ******************************************************************************************

cccz here we start to sort the surface, "SurNodeIndex" eventually = num of surface node,
cccz we keep it here because we want the surface computing system is relatively independent to other 2Dsoil module
cccz start with "SurNodeIndex=0" 
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
      
      do n=1,SurNodeIndex-1
        if(slopeCoord(n,2).ne.slopeCoord(n+1,2)) then
            FlatSurface=1
        endif
      enddo
      
cccz runoff discharge from the left/right edges of the grid
cccz current values; temp-used for iteration; old-record old values;          
      RunoffLeft=0.0D0
      RunoffRight=0.0D0 
      RunoffLeft_temp=0.0D0
      RunoffRight_temp=0.0D0
      RunoffLeft_old=0.0D0
      RunoffRight_old=0.0D0

cccz In the iteration we will not care about the error if the magnitude of the flux or head is below the limiter
      FluxLimit=0.00000000001D0
      HeadLimit=0.00000000001D0   
      
      do i=1,SurNodeIndex-1
        ll(i)=sqrt((abs(slopeCoord(i,1)-slopeCoord(i+1,1))**2.0D0)
     &         +(abs(slopeCoord(i,2)-slopeCoord(i+1,2))**2.0D0))
      enddo
      
cccz record the surface temperature boundary condition
      do i=1, SurNodeIndex
        n=SurfNodeNodeIndexH(i) 
        TmprBCRecord(i)=CodeT(n) 
        GasBCRecord(i)=CodeG(n) 
      enddo
      
      RunoffRight_Runoff_Loc_Record=RunoffRight
      RunoffLeft_Runoff_Loc_Record=RunoffLeft
      
               
      Endif
cccz initialization finished here
cccz ****************************************************************************************** 

C if we have ponded water for infiltration, runoff will not need to be calculated
        if ((PondingByFlux.eq.1).or.(PondingByHead.eq.1)) return
      
cccz sometime the program go backwards, then we should avoid multiple computation of runoff
cccz the simplist way is to track the time and do nothing if the time moves backwards

cccz illustrate some return conditions
cccz if the time step move backwards, return
      if(time.le.time_runoff_old) then
          do n=1,SurNodeIndex
            kk=SurfNodeSurfIndexH(n)
            nn=SurfNodeNodeIndexH(n)
            varbw(kk,1)=varbw_mulch(kk,1)
            varbw(kk,2)=varbw_mulch(kk,2)
            varbw(kk,3)=varbw_mulch(kk,3)
            Q(nn)=-Width(kk)*VarBW_mulch(kk,3)
          enddo    
          RunoffRight=RunoffRight_Runoff_Loc_Record
          RunoffLeft=RunoffLeft_Runoff_Loc_Record
          time_runoff_old=time
          return
      endif
      
      time_runoff_old=time
      
      do n=1,SurNodeIndex
          kk=SurfNodeSurfIndexH(n)
          nn=SurfNodeNodeIndexH(n)
          Q(nn)=-Width(kk)*VarBW_mulch(kk,3)
      enddo          
      
      CriticalH_R=0.1D0 
      if (residueApplied.le.0) then 
         CriticalH_R=0.1D0  
      else
        if(mulchThick.ge.thresholdThick) then
          CriticalH_R=mulchThick*WaterStorageFrac
        else
          CriticalH_R=0.1D0  
        endif 
      endif
      
      
cccz first save the data from previous step
cccz and make initialization
      TotalhStay_test=0.0D0
      TotalqFlux_test=0.0D0
      TotalR_test=0.0D0
      do n=1,SurNodeIndex
        i=SurfNodeSurfIndexH(n)           ! NumBP based indices
        k=SurfNodeNodeIndexH(n)           ! NumP based indices
        h_Pond_Old(n)=h_Pond(n)           ! Save the old ponded depth, h_Pond_Old is within this module
        q_Flux_Old(n)=q_Flux(n)           ! Save the old surface flux (not free one), q_Flux_Old is within this module
        Efflux_Old(n)=Efflux(n)           ! Save the old surface flux (free one), Efflux_Old is within this module
        q_Flux_Node_Old(n)=q_Flux_Node(n) ! Save the node flux (free one), q_Flux_Node_Old is within this module
        RunoffInten(n)=max(RO(k)/width(i),0.0D0)  ! calculate the runoff intensity, i.e. RO/width
        CriticalH_R_Sur(n)=h_Pond(n)      ! The surface critical height is (at least) the ponded depth              
        Exchange_V(n)=0.0D0               ! water exchange for smoothing process
        Check_Exchange_V(n)=0
        FurtherCheck_h_Pond(n)=1   
cccz made some addition to check
        TotalhStay_test=TotalhStay_test+h_Pond(n)
        TotalqFlux_test=TotalqFlux_test+q_Flux(n) 
        TotalR_test=TotalR_test+RunoffInten(n)
      enddo
      RunoffLeft_old=RunoffLeft
      RunoffRight_old=RunoffRight
      RunoffLeft_Efflux_Old=RunoffLeft_Efflux
      RunoffRight_Efflux_Old=RunoffRight_Efflux
     
      
cccz when entering this module, first check if the weather need to be updated.
cccz Varbw_Air records the "weather-based" input because "VarBW" should be changed when the runoff
cccz occurs, because watermov always take "VarBW" as the surface water flux.

          
cccz we have to adjust the rainfall based on the normal direction/vertical direction
cccz we always use the width, but for soil surface with slopes'
cccz the rainfall should be based on horizontal "cross-sectional" surface
cccz Dennis, think if this is necessary
cccz for very small slope, this can be neglected

cccz ------------------------------------------------------------
      if(FlatSurface.ne.0) then
       do n=1,SurNodeIndex
        i=SurfNodeSurfIndexH(n)
        if(n.eq.1) then
          Varbw_Mulch(i,1)=Varbw_Mulch(i,1)*0.5D0*
     &       (slopeCoord(2,1)-slopeCoord(1,1))/width(i)
        elseif(n.gt.1.and.n.lt.SurNodeIndex) then
          Varbw_Mulch(i,1)=Varbw_Mulch(i,1)*0.5D0*
     &       (slopeCoord(n+1,1)-slopeCoord(n-1,1))/width(i) 
        else
          Varbw_Mulch(i,1)=Varbw_Mulch(i,1)*0.5D0*
     &       (slopeCoord(n,1)-slopeCoord(n-1,1))/width(i)
        endif      
       enddo
      endif 
cccz ------------------------------------------------------------      
      
      
cccz the second by-pass condition   
cccz there is no surface ponding water/surface flux, nor runoff (efflux) from soil
      if((TotalhStay_test.le.0.0D0)
     &   .and.(TotalqFlux_test.le.0.0D0)
     &   .and.(TotalR_test.le.0.0D0)) then
cccz resume the original surface temperature boundary condition (in case it is recovered from surface runoff)
         do i=1, SurNodeIndex
           n=SurfNodeNodeIndexH(i) 
           CodeT(n)=TmprBCRecord(i)
           CodeG(n)=GasBCRecord(i)
         enddo  
         return
      endif
            
c ******************************************************************************************          
cccz Start the calculation

cccz if we have simple surface, then calculate it separately
      if(FlatSurface.eq.0) goto 1200
          
cccz need to evaluate criticalH_R for each node based on the hNew, horizontal/vertical distance, and h_Pond.
cccz the idea is look left-side, look right-side, see how large the "water pool" can be
      do n=1,SurNodeIndex
cccz two boundary points do not need this setting,
cccz and we should be prepared for the discharge from the ending point.
        if(n.eq.1) then
          CriticalH_R_Sur(1)=CriticalH_R
        elseif(n.eq.SurNodeIndex) then
          CriticalH_R_Sur(SurNodeIndex)=CriticalH_R
        else
          EdgeReloc=0
          BaseLeft=slopeCoord(n,2)
          BaseRight=slopeCoord(n,2)
          xBaseLeft=slopeCoord(n,1)
          xBaseRight=slopeCoord(n,1)
cccz we trace the point on the left-side
          do i=1,n-1
            k=n-i
cccz if "if condition" is true, k+1 is not the point we want,
cccz because it either lower than current point, or its left point.
            if(slopeCoord(k,2).ge.slopeCoord(k+1,2).or.
     &        BaseLeft.gt.slopeCoord(k+1,2)) then
            else
cccz if "if condition" is true, k+1 “maybe” a local bounday for pools
cccz but the pool now is so large that it can already pass the k+1 and 
cccz go left further
cccz so k+1 is not the left boundary of the pool we want
              if(h_Pond(k+1).gt.CriticalH_R) then 
              else
                BaseLeft=slopeCoord(k+1,2)   
                xBaseLeft=slopeCoord(k+1,1)
                EdgeReloc=1
              endif
             endif 
          enddo
          
cccz we make one more justification for the boundary node
cccz because the left node is not included in the previous loop, since always use k+1
cccz "EdgeReloc=1" means already found the left edge, so not in the if statement
cccz this "if" is designed for a big convex shape of the left side surface
          if(EdgeReloc.eq.0) then
           if(slopeCoord(1,2).ge.slopeCoord(2,2).and.
     &       BaseLeft.lt.slopeCoord(1,2)) then
                BaseLeft=slopeCoord(1,2) 
                xBaseLeft=slopeCoord(1,1)
           endif
          endif
cccz we trace the point on the right-side
          EdgeReloc=0
          do k=n+1,SurNodeIndex
cccz if this "if" is true, "k-1" is not the point we want
             if(slopeCoord(k,2).ge.slopeCoord(k-1,2).or.
     &            BaseRight.gt.slopeCoord(k-1,2)) then
             else
cccz if this "if" is true, "k-1" is already under water
cccz need to look further on the right side
               if(h_Pond(k-1).gt.CriticalH_R) then  
               else
                 BaseRight=slopeCoord(k-1,2)
                 xBaseRight=slopeCoord(k-1,1)
                 EdgeReloc=1
               endif
            endif 
          enddo 

cccz we make one more justification for the boundary node
cccz because the right node is not included in the previous loop, since always use k-1
cccz "EdgeReloc=1" means already found the right edge, so not in the if statement
cccz this "if" is designed for a big convex shape of the right side surface
          if(EdgeReloc.eq.0) then
            if(slopeCoord(SurNodeIndex,2).ge.
     &           slopeCoord(SurNodeIndex-1,2).and.
     &           BaseRight.lt.slopeCoord(SurNodeIndex,2)) then
                BaseRight=slopeCoord(SurNodeIndex,2) 
                xBaseRight=slopeCoord(SurNodeIndex,1)
            endif
          endif

cccz now we compare the left-searching and right-searching result                   
          if(xBaseLeft.eq.xBaseRight) then
            CriticalH_R_Sur(n)=min(BaseLeft,BaseRight)
     &        -slopeCoord(n,2)+CriticalH_R
          else
cccz in this case, left side is higher
            if(BaseLeft.gt.BaseRight) then
              CriticalH_R_Sur(n)=min(BaseLeft,BaseRight)
     &         -slopeCoord(n,2)+CriticalH_R
cccz here we make a linear interpolation, we allow a small angel from the higher end to the lower end
     &         +0.5*CriticalH_R*(xBaseRight-slopeCoord(n,1))
     &         /(xBaseRight-xBaseLeft)
cccz in this case, right side is higher
            elseif(BaseLeft.lt.BaseRight) then
              CriticalH_R_Sur(n)=min(BaseLeft,BaseRight)
     &         -slopeCoord(n,2)+CriticalH_R
cccz here we make a linear interpolation, we allow a small angel from the higher end to the lower end
     &         +0.5*CriticalH_R*(slopeCoord(n,1)-xBaseLeft)
     &         /(xBaseRight-xBaseLeft)
cccz both sides are equal-height
            else
              CriticalH_R_Sur(n)=min(BaseLeft,BaseRight)
     &         -slopeCoord(n,2)+CriticalH_R
            endif
          endif
         endif
        enddo
          
cccz define this value for numerical reasons
cccz the computer sometimes set "criticalH=-0.01" as "criticalH=-0.0099999"
cccz thus we need to buffer it
cccz          CriticalH_S=CriticalH-0.000001D0      

          
cccz start to solve the SV equations
cccz We need some hard judgement for the unit, we use "cm" and "day" 
        iteration_Num=0          
1101    iteration_Flux=0
        iteration_Head=0
        RunoffLeft_temp=0.0D0
        RunoffRight_temp=0.0D0
        do i=1,NumBPD
          h_Pond_temp(i)=0.0D0
          q_Flux_temp(i)=0.0D0
          q_Flux_Node(i)=0.0D0
          Efflux(i)=0.0D0
        enddo
        RunoffLeft_Efflux=0.0D0
        RunoffRight_Efflux=0.0D0
      
        do n=1,SurNodeIndex-1
cccz use iteration to partition the real runoff and amount of water stay 
cccz take the absolute value of the slope to ensure the following calculation
          Slope_0_R=abs((slopeCoord(n+1,2)-slopeCoord(n,2))
     &           /(slopeCoord(n,1)-slopeCoord(n+1,1)))
cccz start the explicit version for the q_flux
cccz we also check the right/left directions and made average for stability
          if(h_Pond(n+1).eq.0.0D0) then
             Slope_f_left=0.0D0
          else
             Slope_f_left=((n_Stiff*q_Flux(n))**2.0D0)
     &          /(h_Pond(n+1)**(3.333D0))
          endif
          if(h_Pond(n).eq.0.0D0) then
             Slope_f_right=0.0D0
          else
             Slope_f_right=((n_Stiff*q_Flux(n)**2.0D0)
     &         /(h_Pond(n)**(3.333D0)))
          endif
          Slope_f=0.50D0*(Slope_f_left+Slope_f_right)

cccz now we calcualte the fluxes
cccz start from the first (left) node
          if(n.eq.1) then
cccz first determine the width between two point,
cccz for the partial derivative along the edge.
c            ll=sqrt((abs(slopeCoord(1,1)-slopeCoord(2,1))**2.0D0)
c     &         +(abs(slopeCoord(1,2)-slopeCoord(2,2))**2.0D0))
cccz Water will flow from node 1 to node 2
cccz Always follow the upwind direction
           if((slopeCoord(1,2)+h_Pond(1)).ge.
     &        (slopeCoord(2,2)+h_Pond(2))) then
            if(h_Pond(1).eq.0.0D0) then
              q_flux_square=0.0D0
              h_flux_square=10.0D0       ! arbitrary none-zero number here
            else
              q_flux_square=q_Flux(1)
              h_flux_square=h_Pond(1)
            endif
cccz we can see for values labeled "_Old" is used to store old data
cccz which is not changed during the iteration
cccz this setting will be also important if the time-scale varies
cccz e.g., the mulch model

            q_Flux_temp(1)=q_Flux_Old(1)+
     &         (-(q_flux_square**2.0D0)/h_flux_square/ll(1)      ! there is no flux on the left of left edge
     &         +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1)
     &         +g_accel*h_Pond(1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(1)=max(q_Flux_temp(1),0.0D0)
            
cccz the direction implies runoff will not be discharged from the left side
            RunoffLeft_temp=0.0D0
           else
cccz Water will flow from Node 2 to Node 1, momentum source from q_flux(2), between Node 2 and Node 3
cccz this depends on whether Node 2 is a local maximum                   
            if ((slopeCoord(2,2)+h_Pond(2)).ge.
     &        (slopeCoord(3,2)+h_Pond(3))) then  
     
cccz now we calcualte water flow     
             if(h_Pond(2).eq.0.0D0) then
               q_flux_square=0.0D0
               h_flux_square=10.0D0                                   ! arbitrary number here
             else
               q_flux_square=q_Flux(1)
               h_flux_square=h_Pond(2)
             endif
       
             q_Flux_temp(1)=q_Flux_Old(1)+
     &          ((q_flux_square**2.0D0)/h_flux_square/ll(1)       ! there is no income momentum from q_flux(2)
     &          +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1)
     &          -g_accel*h_Pond(2)*max(Slope_0_R-Slope_f,0.0D0))*step
             q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
             
            else
cccz in this case, there is momentum source from q_flux(2)                 
             if(h_Pond(2).eq.0.0D0) then
               q_flux_square1=0.0D0
               h_flux_square1=10.0D0                                  ! arbitrary number here
             else
               q_flux_square1=q_Flux(1)
               h_flux_square1=h_Pond(2)
             endif
             if(h_Pond(3).eq.0.0D0) then
               q_flux_square2=0.0D0
               h_flux_square2=10.0D0                                  ! arbitrary number here
             else
               q_flux_square2=q_Flux(2)
               h_flux_square2=h_Pond(3)
             endif
                     
             q_Flux_temp(1)=q_Flux_Old(1)+
     &         (((q_flux_square1**2.0D0)/h_flux_square1
     &              -(q_flux_square2**2.0D0)/h_flux_square2)/ll(1)
     &         +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1) 
     &         -g_accel*h_Pond(2)*max(Slope_0_R-Slope_f,0.0D0))*step 
             q_Flux_temp(1)=min(q_Flux_temp(1),0.0D0)
           endif
cccz finish the calculation of q in the case Node 2 --> Node 1
cccz then calculate the runoff from the left edge
           if(h_Pond(1).eq.0.0D0) then
             q_flux_leftrunoff=0.0D0
             h_flux_leftrunoff=10.0D0                                 ! arbitrary number here
           else
             q_flux_leftrunoff=RunoffLeft
             h_flux_leftrunoff=h_Pond(1)
           endif
           if(h_Pond(2).eq.0.0D0) then
             q_flux_square=0.0D0
             h_flux_square=10.0D0                                     ! arbitrary number here
           else
             q_flux_square=q_Flux(1)
             h_flux_square=h_Pond(2)
           endif
           RunoffLeft_temp=RunoffLeft_old+
     &       (((q_flux_leftrunoff**2.0D0)/h_flux_leftrunoff
     &          -(q_flux_square**2.0D0)/h_flux_square)
     &          /width(SurfNodeSurfIndexH(n))
     &       +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1) 
     &       -g_accel*h_Pond(1)*max(Slope_0_R-Slope_f,0.0D0))*step 
           RunoffLeft_temp=min(RunoffLeft_temp,0.0D0)
          endif

cccz finished the calculation when n=1 (left-most) edge
cccz now calculate the flux in interior nodes
         elseif(n.gt.1.and.n.lt.SurNodeIndex-1) then

cccz first calculate the 
c          ll=sqrt(abs(slopeCoord(n,1)-slopeCoord(n+1,1))**2.0D0
c     &      +abs(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0)
     
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
cccz water will flow from Node n -> n+1, right-bound, momentun source may from n-1?
cccz always follow the upwind direction 
cccz in the first case, Node n is the local peak s.t. no momentun from the left side 
           if((slopeCoord(n,2)+h_Pond(n)).ge.
     &        (slopeCoord(n-1,2)+h_Pond(n-1))) then   
cccz the momentum is only from the current slope section, i.e., n
            if(h_Pond(n).eq.0.0D0) then
              q_flux_square=0.0D0
              h_flux_square=10.0D0                                    ! arbitrary none-zero number here
            else
              q_flux_square=q_Flux(n)
              h_flux_square=h_Pond(n)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        (-(q_flux_square**2.0D0)/h_flux_square/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
           
           else
cccz the momentum is from the left side, i.e., n-1                  
            if(h_Pond(n).eq.0.0D0) then
              q_flux_square1=0.0D0
              h_flux_square1=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square1=q_Flux(n)
              h_flux_square1=h_Pond(n)
            endif
            if(h_Pond(n-1).eq.0.0D0) then
              q_flux_square2=0.0D0
              h_flux_square2=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square2=q_Flux(n-1)
              h_flux_square2=h_Pond(n-1)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        ((-(q_flux_square1**2.0D0)/h_flux_square1
     &           +(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step  
            q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
           endif
         else
cccz Water will flow from Node n+1 -> n, momentun source of q_Flux(n+1)?
cccz always follow the upwind direction 
cccz now we assume n+1 is a local maximum, so no water flow Node n+2->n+1
cccz thus, there is no external momentum from q_Flux(n+1), which should be in an opposite direction
           if((slopeCoord(n+1,2)+h_Pond(n+1)).ge.
     &        (slopeCoord(n+2,2)+h_Pond(n+2))) then
             if(h_Pond(n+1).eq.0.0D0) then
               q_flux_square=0.0D0
               h_flux_square=10.0D0                                    ! arbitrary none-zero number here
             else
               q_flux_square=q_Flux(n)
               h_flux_square=h_Pond(n+1)
             endif
                                         
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        ((q_flux_square**2.0D0)/h_flux_square/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        -g_accel*h_Pond(n+1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
           else    
cccz now we assume n+2 is higher than n+1
cccz thus, there is external momentum from q_Flux(n+1), from Node n+2
            if(h_Pond(n+1).eq.0.0D0) then
              q_flux_square1=0.0D0
              h_flux_square1=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square1=q_Flux(n)
              h_flux_square1=h_Pond(n+1)
            endif
            if(h_Pond(n+2).eq.0.0D0) then
              q_flux_square2=0.0D0
              h_flux_square2=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square2=q_Flux(n+1)
              h_flux_square2=h_Pond(n+2)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        (((q_flux_square1**2.0D0)/h_flux_square1
     &           -(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &       -g_accel*h_Pond(n+1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                     
          endif
         endif

cccz now we discuss the right-most point, with flux q_flux(n=SurNodeIndex-1)
        else
       
c         ll=sqrt((abs(slopeCoord(n,1)-slopeCoord(n+1,1))**2.0D0)
c     &      +(abs(slopeCoord(n,2)-slopeCoord(n+1,2))**2.0D0))
             
cccz Finish the calculation for n=SurNodeIndex-1
cccz or in other words, we have "n=SurNodeIndex-1" here, so you will not see "n+2"
         if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
cccz water flow to the end point n=SurNodeIndex, need to calculate the momentun from n-2?
cccz n=SurNodeIndex-1 is the local maxima, no water flux from left side 
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &        (slopeCoord(n-1,2)+h_Pond(n-1))) then 
           if(h_Pond(n).eq.0.0D0) then
             q_flux_square=0.0D0
             h_flux_square=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square=q_Flux(n)
             h_flux_square=h_Pond(n)
           endif
                       
           q_Flux_temp(n)=q_Flux_Old(n)+
     &        (-(q_flux_square**2.0D0)/h_flux_square/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
           q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
          else
cccz n=SurNodeIndex-2 is higher, water flux from left side                   
           if(h_Pond(n).eq.0.0D0) then
             q_flux_square1=0.0D0
             h_flux_square1=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square1=q_Flux(n)
             h_flux_square1=h_Pond(n)
           endif    
           if(h_Pond(n-1).eq.0.0D0) then
             q_flux_square2=0.0D0
             h_flux_square2=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square2=q_Flux(n-1)
             h_flux_square2=h_Pond(n-1)
           endif
                  
           q_Flux_temp(n)=q_Flux_Old(n)+
     &       ((-(q_flux_square1**2.0D0)/h_flux_square1
     &           +(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &       +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
           q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
          endif

cccz in this case, runoff discharge on the right hand side occurs
cccz this calculation should be processed based on the SV equation
          if(h_Pond(n+1).eq.0.0D0) then
            q_flux_runoffright=0.0D0
            h_flux_runoffright=10.0D0       ! arbitrary number here
          else
            q_flux_runoffright=RunoffRight
            h_flux_runoffright=h_Pond(SurNodeIndex)
          endif    
          if(h_Pond(n).eq.0.0D0) then
            q_flux_square=0.0D0
            h_flux_square=10.0D0       ! arbitrary number here
          else
            q_flux_square=q_Flux(n)
            h_flux_square=h_Pond(n)
          endif                  

          RunoffRight_temp=RunoffRight_old+
     &      ((-(q_flux_runoffright**2.0D0)/h_flux_runoffright
     &         +(q_flux_square**2.0D0)/h_flux_square)
     &         /width(SurfNodeSurfIndexH(SurNodeIndex))
     &      +0.5D0*g_accel*(h_Pond(n)**2.0D0
     &         -h_Pond(SurNodeIndex)**2.0D0)/ll(n) 
     &      +g_accel*h_Pond(SurNodeIndex)*max(Slope_0_R-Slope_f,0.0D0))
     &         *step 
          RunoffRight_temp=max(RunoffRight_temp,0.0D0)
                             
         else

cccz water flow to from Node SurNodeIndex -> SurNodeIndex-1
cccz in this case, no right runoff discharge
cccz also, a good news is there is no right-hand node, thus, no external momentum source

          if(h_Pond(n+1).eq.0.0D0) then
            q_flux_square=0.0D0
            h_flux_square=10.0D0       ! arbitrary number here
          else
            q_flux_square=q_Flux(n)
            h_flux_square=h_Pond(SurNodeIndex)
          endif

          q_Flux_temp(n)=q_Flux_Old(n)+
     &       ((q_flux_square**2.0D0)/h_flux_square/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0
     &          -h_Pond(SurNodeIndex)**2.0D0)/ll(n) 
     &       -g_accel*h_Pond(SurNodeIndex)*max(Slope_0_R-Slope_f,0.0D0))
     &          *step
          q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
          RunoffRight_temp=0.0D0
        endif
       endif
      enddo
        
cccz finishe the flux estimation, now go the ponded water head estimation
cccz The following steps will make the h_temp goes infinitely high
cccz so after solving the SV equation, we need additional processes

      do n=1,SurNodeIndex
        i=SurfNodeSurfIndexH(n)
        k=SurfNodeNodeIndexH(n)
        if(n.eq.1) then
         h_Pond_temp(1)=h_Pond_Old(1)+step*
     &     ((-q_Flux_temp(1)+RunoffLeft_temp)/width(i)+RunoffInten(1))
         h_Pond_temp(1)=max(h_Pond_temp(1),0.0D0)

cccz calculate the amound of efflux -- which means "free flux"         
         if((slopeCoord(1,2)+h_Pond(1)).le.
     &      (slopeCoord(2,2)+h_Pond(2))) then
           if(h_Pond_temp(1).gt.CriticalH_R_Sur(1)) then
cccz for boundary point, it is easy
cccz the h_Pond should not exceed criticalH_R
             RunoffLeft_Efflux=RunoffLeft_Efflux
     &         -(h_Pond_temp(1)-CriticalH_R_Sur(1))*width(i)/step
             h_Pond_temp(1)=min(h_Pond_temp(1),CriticalH_R_Sur(1))
           endif
         else
           if(h_Pond_temp(1).gt.CriticalH_R_Sur(1)) then
             Efflux(1)=Efflux(1)
     &         +(h_Pond_temp(1)-CriticalH_R_Sur(1))*width(i)/step
             h_Pond_temp(1)=min(h_Pond_temp(1),CriticalH_R_Sur(1))
             RunoffLeft_temp=0.0D0
             RunoffLeft_Efflux=0.0D0
           endif 
         endif
         
cccz for the interior nodes, the efflux could be left/right/outwards/inwards
        elseif(n.gt.1.and.n.lt.SurNodeIndex) then
          h_Pond_temp(n)=h_Pond_Old(n)+step*
     &      ((-q_Flux_temp(n)+q_Flux_temp(n-1))/width(i)+RunoffInten(n))
          h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)     

          if((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
cccz Previuos time step will not take care of EFFlux(n-1)
           if((slopeCoord(n,2)+h_Pond(n)).le.
     &        (slopeCoord(n+1,2)+h_Pond(n+1))) then

cccz in this geometry, water will flow left in "n"             
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               Efflux(n-1)=Efflux(n-1)
     &           -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n)) 
             endif

cccz in this geometry, water will flow outwards like "/\" 
cccz but how the runoff partitioned to left/right sides
           else
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz first calculate all avail runoff water
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
cccz calculate left/right slope
              Slope_0_n_1=abs(((slopeCoord(n-1,2)+h_Pond(n-1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n-1,1)))
              Slope_0_n=abs(((slopeCoord(n+1,2)+h_Pond(n+1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n+1,1)))
cccz use the slope to determine the fraction
              Fraction_n_1=Slope_0_n_1/(Slope_0_n_1+Slope_0_n)
              Fraction_n=Slope_0_n/(Slope_0_n_1+Slope_0_n)
cccz calculate the Efflux based on the fraction
              Efflux(n-1)=Efflux(n-1)-Fraction_n_1*Add_Flux
              Efflux(n)=Efflux(n)+Fraction_n*Add_Flux
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
             endif
           endif
         elseif((slopeCoord(n,2)+h_Pond(n)).lt.
     &     (slopeCoord(n-1,2)+h_Pond(n-1))) then
cccz The q_Flux(n-1) are taken cared by previous steps, then just calculate new q_Flux(n)
cccz if and only if q_Flux(n)>0
cccz i.e., calculate rightwards efflux on the rightside
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Efflux(n)=Efflux(n)
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
            endif
cccz here we deal with a special case -- a valley
cccz if there is a valley, the efflux water will not go anywhere,
cccz but we should record it and be prepared to add it as water input for the next iteration
          else
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              q_Flux_Node(n)=(h_Pond_temp(n)-CriticalH_R_Sur(n))
     &           *width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
            endif
          endif

cccz now the elevation of the water level @ Node n and Node n-1
cccz was the same.

         else
           if((slopeCoord(n,2)+h_Pond(n)).lt.
     &       (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz the rightside is higher, momentum points leftwards
              Efflux(n-1)=Efflux(n-1)
     &           -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))   
            endif
           elseif((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz the rightside is lower, momentum points rightwards
              Efflux(n)=Efflux(n)
     &          +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step   
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n)) 
            endif
           else
cccz totally flat  
cccz first record the potential efflux quantity
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
cccz need to refer the momentum from previous step
              if(Efflux_old(n-1).gt.0.0D0) then
               if(Efflux_old(n).ge.0.0D0) then
                 Efflux(n)=Efflux(n)+Add_Flux
               else
cccz follow the side with larger momentum
                if(abs(Efflux_old(n)).lt.Efflux_old(n-1)) then
                  Efflux(n)=Efflux(n)+Add_Flux
                elseif(abs(Efflux_old(n)).gt.Efflux_old(n-1)) then
                  Efflux(n-1)=Efflux(n-1)-Add_Flux
                else
                  Efflux(n-1)=Efflux(n-1)-0.50D0*Add_Flux
                  Efflux(n)=Efflux(n)+0.50D0*Add_Flux
                endif
               endif
              elseif (Efflux_old(n-1).le.0.0D0) then
               if (Efflux_old(n).lt.0.0D0) then
                  Efflux(n-1)=Efflux(n-1)-Add_Flux
               else
cccz try average three point----OK
cccz                  Efflux(n-1)=Efflux(n-1)-0.50D0*Add_Flux
cccz                  Efflux(n)=Efflux(n)+0.50D0*Add_Flux
                  Efflux(n-1)=Efflux(n-1)-Add_Flux/3.0D0
                  Efflux(n)=Efflux(n)+Add_Flux/3.0D0
                  q_Flux_Node(n)=q_Flux_Node(n)+Add_Flux/3.0D0
               endif
              endif
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
             endif
           endif
          endif
                
cccz finally, we calcualte the h_Pond for the right-most point
         else
          h_Pond_temp(n)=h_Pond_Old(n)+step*
     &    ((-RunoffRight_temp+q_Flux_temp(n-1))/width(i)+RunoffInten(n))
          h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
cccz right runoff discharge case
          if((slopeCoord(n,2)+h_Pond(n)).le.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               RunoffRight_Efflux=RunoffRight_Efflux
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))    
            endif
          elseif((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Efflux(n-1)=Efflux(n-1)
     &          -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
              RunoffRight_temp=0.0D0
              RunoffRight_Efflux=0.0D0
            endif
          else              
         endif
        endif
       enddo

cccz finish update of head, need further work for setting the critical_H_Sur value for each node.
cccz Error esitimation for iteration, i.e.,
cccz the Picard's for SV equation
          
       iteration_Flux=0
       iteration_Head=0
       iteration_Num=iteration_Num+1
cccz   iteration_Num=1.0D0
cccz take average between the current iteration and previous iteration steps to seek for stability
       do n=1,SurNodeIndex
         q_Flux_temp(n)=q_Flux_temp(n)/iteration_Num+
     &     (iteration_Num-1)*q_Flux(n)/iteration_Num
         h_Pond_temp(n)=h_Pond_temp(n)/iteration_Num+
     &     (iteration_Num-1)*h_Pond(n)/iteration_Num
         h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
       enddo
       RunoffRight_temp=RunoffRight_temp/iteration_Num+
     &   (iteration_Num-1)*RunoffRight/iteration_Num
       RunoffLeft_temp=RunoffLeft_temp/iteration_Num+
     &   (iteration_Num-1)*RunoffLeft/iteration_Num
          
       do n=1,SurNodeIndex
        if(abs(q_Flux_temp(n)-q_Flux(n)).gt.0.01*abs(q_Flux(n))) then
          if(abs(q_Flux_temp(n)).gt.FluxLimit)  iteration_Flux=1
        endif
        if(abs(h_Pond_temp(n)-h_Pond(n)).gt.0.01*abs(h_Pond(n))) then
          if(abs(h_Pond_temp(n)).gt.HeadLimit)  iteration_Head=1
        endif
       enddo
          
       if(abs(RunoffLeft_temp-RunoffLeft).gt.0.01*abs(RunoffLeft)) then
         if(abs(RunoffLeft).gt.FluxLimit) iteration_Flux=1
       endif
       if(abs(RunoffRight_temp-RunoffRight).gt.0.01*abs(RunoffRight)) 
     &  then            
         if(abs(RunoffRight).gt.FluxLimit) iteration_Flux=1
       endif

       if (iteration_Flux.eq.1.or.iteration_Head.eq.1) then
        do n=1,SurNodeIndex
          q_Flux(n)=q_Flux_temp(n)
          h_Pond(n)=h_Pond_temp(n)
        enddo
          RunoffLeft=RunoffLeft_temp
          RunoffRight=RunoffRight_temp
        goto 1101
       else
        do n=1,SurNodeIndex
          q_Flux(n)=q_Flux_temp(n)
          h_Pond(n)=h_Pond_temp(n)
          FurtherCheck_h_Pond(n)=1
        enddo
        RunoffLeft=RunoffLeft_temp
        RunoffRight=RunoffRight_temp
       endif


cccz ask a question, why use "q_Flux_Node_Old" here?
       do n=1,SurNodeIndex
c         k=SurfNodeNodeIndexH(n)              ! for index in the whole node set
c         i=SurfNodeSurfIndexH(n)              ! for index in the boundary node set
c         VarBW(i,1)=Varbw_Mulch(i,1)+q_Flux_Node_Old(n)/width(i)
c         VarBW(i,3)=Varbw_Mulch(i,2)-VarBW(i,1)                 
c         Q(k)=-Width(i)*VarBW(i,3)                   

cccz there is more water can be used to redistribute/infiltration         
        if(h_Pond(n).ge.h_Pond_Old(n).and.h_Pond(n).ge.CriticalH_R) then
          FurtherCheck_h_Pond(n)=1
        else
          FurtherCheck_h_Pond(n)=0
        endif
       enddo

cccz need to further check h_Pond(n)
cccz fix a good h_Pond value based on the topology of the domain
       iteration_Num=0.0D0
1105   do n=1,NumBPD
        Exchange_V_temp(n)=0.0D0
        Check_Exchange_V(n)=0
       enddo
cccz the boundary points does not need to checked, because there is always runoff discharge
cccz the internal points need to be checked "pool-by-pool"
       do n=2,SurNodeIndex-1
         if(FurtherCheck_h_Pond(n).eq.1) then
cccz record the water level on the left/self/right          
           H_L=slopeCoord(n-1,2)+h_Pond(n-1)
           H_M=slopeCoord(n,2)+h_Pond(n)
           H_R=slopeCoord(n+1,2)+h_Pond(n+1)
cccz record the scale of width to the left/right node
           if(n.eq.2) then
            L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
           else
            L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-2,1))
           endif
           if(n.eq.(SurNodeIndex-1)) then
            L_R=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
           else
            L_R=0.5*(slopeCoord(n+2,1)-slopeCoord(n,1))
           endif
           L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
cccz valley case @ n              
           if(H_M.lt.H_L.and.H_M.lt.H_R) then
            Ava_Height=(H_R*L_R+H_L*L_L+H_M*L_M)/(L_R+L_L+L_M)    
cccz the amount of water need for the middle node
            Volume=max((Ava_Height-
     &         (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
cccz check avaliable water from both sides
cccz ask a question, here try to use "min((h_Pond_temp(n-1)-CriticalH_R),(H_L-Ava_Height))" in future
c            Ava_Volume_L=(h_Pond_temp(n-1)-CriticalH_R)*L_L
c            Ava_Volume_R=(h_Pond_temp(n+1)-CriticalH_R)*L_R
            Ava_Volume_L=min((h_Pond_temp(n-1)-CriticalH_R),
     &         (H_L-Ava_Height))*L_L
            Ava_Volume_R=min((h_Pond_temp(n+1)-CriticalH_R),
     &         (H_R-Ava_Height))*L_R
            Ava_Volume=Ava_Volume_L+Ava_Volume_R
cccz there is enough water 
cccz but we have to count the water contribution from both side
           if(Ava_Volume.ge.Volume) then
             if(Ava_Volume_L.le.(Volume/2.0D0)) then
               Volume_L=Ava_Volume_L
               Volume_R=Volume-Volume_L
             elseif(Ava_Volume_R.le.(Volume/2.0D0)) then
               Volume_R=Ava_Volume_R
               Volume_L=Volume-Volume_R
             else
               Volume_L=Volume/2.0D0
               Volume_R=Volume/2.0D0
             endif
cccz the water is not enough, bring all the water in
           else
             Volume_L=Ava_Volume_L
             Volume_R=Ava_Volume_R             
           endif
           if(Check_Exchange_V(n-1).eq.0) then
             Exchange_V_temp(n-1)=Volume_L
             Check_Exchange_V(n-1)=1
           else
cccz this process focused on the computing stability
            if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
            else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
            endif
           endif
           if(Check_Exchange_V(n).eq.0) then
              Exchange_V_temp(n)=-Volume_R
              Check_Exchange_V(n)=1
           else
            if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)
            else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)     
            endif
           endif
cccz this case indicate a leftwards slope     
          elseif(H_M.ge.H_L.and.H_M.le.H_R) then
            Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
            Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
            Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)
cccz convex shape or concave shape
            if(H_M.ge.Mean_Height) then
              Ava_Height_L=Mean_Height
              Ava_Height_R=Mean_Height_R
            else
              Ava_Height_L=Mean_Height_L
              Ava_Height_R=Mean_Height
            endif
cccz how much water I need
            Volume_L=max((Ava_Height_L-
     &        (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
cccz how much I have
            Volume_L=min(Volume_L,
     &        max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_L,
     &        h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
cccz similar idea to the left side     
            Volume_R=max((Ava_Height_R-
     &        (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
            Volume_R=min(Volume_R,
     &        max(min(slopeCoord(n+1,2)+h_Pond_temp(n+1)-Ava_Height_R,
     &        h_Pond_temp(n+1)-CriticalH_R)*L_R,0.0D0))
                       
            if(Check_Exchange_V(n-1).eq.0) then
              Exchange_V_temp(n-1)=-Volume_L
              Check_Exchange_V(n-1)=1
            else
             if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)
             else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L) 
             endif
            endif
                  
            if(Check_Exchange_V(n).eq.0) then
              Exchange_V_temp(n)=-Volume_R
              Check_Exchange_V(n)=1
            else
             if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)    
             else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)
             endif
            endif   
            
cccz rightwards slope
          elseif(H_M.le.H_L.and.H_M.ge.H_R) then
            Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
            Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
            Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)

            if(H_M.ge.Mean_Height) then
              Ava_Height_R=Mean_Height
              Ava_Height_L=Mean_Height_L
            else
              Ava_Height_R=Mean_Height_R
              Ava_Height_L=Mean_Height
            endif
          
            Volume_R=max((Ava_Height_R-
     &        (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
            Volume_R=min(Volume_R,
     &        max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_R,
     &        h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
            Volume_L=max((Ava_Height_L-
     &        (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
            Volume_L=min(Volume_L,
     &        max(min(slopeCoord(n-1,2)+h_Pond_temp(n-1)-Ava_Height_L,
     &        h_Pond_temp(n-1)-CriticalH_R)*L_L,0.0D0))
                       
            if(Check_Exchange_V(n-1).eq.0) then
              Exchange_V_temp(n-1)=Volume_L
              Check_Exchange_V(n-1)=1
            else
             if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
             else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
             endif
            endif
                  
            if(Check_Exchange_V(n).eq.0) then
              Exchange_V_temp(n)=Volume_R
              Check_Exchange_V(n)=1
            else
             if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
             else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
             endif
            endif   

cccz now we do the peak              
          else
            Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
            Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
            Mean_Height=(H_L*L_L+H_M*L_M+H_R*L_R)/(L_L+L_M+L_R)

            if(Mean_Height_L.lt.Mean_Height) then
              Ava_Height_L=Mean_Height
              Ava_Height_R=Mean_Height_R
            else
              Ava_Height_R=Mean_Height
              Ava_Height_L=Mean_Height_L
            endif
                  
            Volume_L=max((Ava_Height_L-
     &        (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
            Volume_R=max((Ava_Height_R-
     &        (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
cccz the total volume we need     
            Volume_T=Volume_L+Volume_R
cccz the total volume we have
            Volume=min(Volume_T,
     &        max(min(slopeCoord(n,2)+h_Pond_temp(n)-
     &                min(Mean_Height_L,Mean_Height_R),
     &        h_Pond_temp(n)-CriticalH_R)*L_L,0.0D0))
     
            if(Volume_T.eq.0.0D0) then
            else
              Volume_L=Volume*Volume_L/Volume_T
              Volume_R=Volume*Volume_R/Volume_T
            endif
                                    
            if(Check_Exchange_V(n-1).eq.0) then
              Exchange_V_temp(n-1)=-Volume_L
              Check_Exchange_V(n-1)=1
            else
             if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)
             else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L)
             endif
            endif
                  
            if(Check_Exchange_V(n).eq.0) then
              Exchange_V_temp(n)=Volume_R
              Check_Exchange_V(n)=1
            else
             if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
             else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
             endif
            endif   
cccz finish the calculation of flux exchange                           
          endif     
         endif
        enddo

cccz the water at the edge had to been drained.
cccz if Exchange_V_temp(1)>0, it will be wrong because that part of water should be in Efflux
cccz similar to the Exchange_V_temp(n), which cannot be <0
        Exchange_V_temp(1)=min(Exchange_V_temp(1),0.0D0)
        Exchange_V_temp(n)=max(Exchange_V_temp(n),0.0D0)
        error=0.0D0
              
        do n=1,SurNodeIndex
         if(n.eq.1) then
           L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
         elseif(n.eq.SurNodeIndex) then
           L_M=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
         else
           L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
         endif
          
         if(n.eq.1) then
           h_Pond_temp(n)=h_Pond_temp(n)-Exchange_V_temp(n)/L_M
           if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            RunoffLeft_Efflux=RunoffLeft_Efflux
     &       -(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
            h_Pond_temp(n)=CriticalH_R_Sur(n)
           endif
         elseif(n.eq.SurNodeIndex) then
           h_Pond_temp(n)=h_Pond_temp(n)+Exchange_V_temp(n-1)/L_M
           if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            RunoffRight_Efflux=RunoffRight_Efflux
     &        +(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
            h_Pond_temp(n)=CriticalH_R_Sur(n)
           endif
         else
           h_Pond_temp(n)=h_Pond_temp(n)+(Exchange_V_temp(n-1)
     &        -Exchange_V_temp(n))/L_M          
         endif             
         error=max(error,abs(h_Pond_temp(n)-h_Pond(n)))          
        enddo
cccz finish update the water level after post-adjustment (one iteration)     
        iteration_Num=iteration_Num+1
cccz determine if additional iteration is needed, max 500 iteration was allowed
cccz may need to take adjustment for the max iteration number
cccz because it cause some computing time.
       if((error.gt.0.01D0*CriticalH_R).and.(iteration_Num.lt.100)) then
        do n=1,SurNodeIndex
cccz ask a question,
cccz in order to maintian mass balance,
cccz should not use
cccz try the two way and think about it
cccz doable in the way without average, but for final calculation, have to use the one without average
cccz since mass balance has to be processed.
cccz      h_Pond(n)=(h_Pond_temp(n)+h_Pond(n)*(iteration_Num-1))/iteration_Num
         h_Pond(n)=h_Pond_temp(n)
         Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
        enddo
        h_Pond(1)=CriticalH_R_Sur(1)
        h_Pond(SurNodeIndex)=CriticalH_R_Sur(SurNodeIndex)
        goto 1105
       else
        do n=1,SurNodeIndex
          h_Pond(n)=h_Pond_temp(n)
          Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
        enddo
       endif
          
cccz after the h_Pond adjustment, we should have a better idea of CriticalH_R_Sur
cccz now set the new CriticalH_R_Sur
cccz and we should make h_Pond and flux calculation based on the new CriticalH_R_Sur
       do n=1,SurNodeIndex
         CriticalH_R_Sur(n)=max(h_Pond(n),CriticalH_R)
         Efflux(n)=0.0D0
         q_Flux_Node(n)=0.0D0
       enddo
          
       do n=2,SurNodeIndex-1
         k=SurfNodeNodeIndexH(n)
cccz         H_L=slopeCoord(n-1,2)+h_Pond_temp(n-1)
cccz         H_M=slopeCoord(n,2)+h_Pond_temp(n)
cccz         H_R=slopeCoord(n+1,2)+h_Pond_temp(n+1)
         H_L=slopeCoord(n-1,2)+h_Pond(n-1)
         H_M=slopeCoord(n,2)+h_Pond(n)
         H_R=slopeCoord(n+1,2)+h_Pond(n+1)
        if(H_M.lt.H_L.and.H_M.lt.H_R) then
         CriticalH_R_Sur(n)=(H_L*(slopeCoord(n+1,1)-slopeCoord(n,1))
     &    +H_R*(slopeCoord(n,1)-slopeCoord(n-1,1)))
     &    /(slopeCoord(n+1,1)-slopeCoord(n-1,1))-slopeCoord(n,2)
cccz the volume water needed at the center point to reach maximum volume 
         DeltaVolume=(max(slopeCoord(n-1,2)+h_Pond(n-1),
     &     slopeCoord(n+1,2)+h_Pond(n+1))-(slopeCoord(n,2)+h_Pond(n)))
     &     *0.5D0*(slopeCoord(n+1,1)-slopeCoord(n-1,1))               
cccz the runoff @ this point is enough to support is amount of water     
         if(DeltaVolume/step.lt.RO(k)) then
          CriticalH_R_Sur(n)=max(slopeCoord(n-1,2)+h_Pond(n-1),
     &     slopeCoord(n+1,2)+h_Pond(n+1))-slopeCoord(n,2)
         endif
        endif
cccz the target node is at the edge of the "local pool"
        if(H_M.le.H_L.and.H_M.lt.slopeCoord(n+1,2)) then
          CriticalH_R_Sur(n)=max(CriticalH_R_Sur(n-1)+slopeCoord(n-1,2)-
     &      slopeCoord(n,2)+CriticalH_R*0.1D0,CriticalH_R*1.1D0)
        elseif(H_M.le.H_R.and.H_M.lt.slopeCoord(n-1,2)) then
          CriticalH_R_Sur(n)=max(CriticalH_R_Sur(n+1)+slopeCoord(n+1,2)-
     &      slopeCoord(n,2)+CriticalH_R*0.1D0,CriticalH_R*1.1D0)
        endif
       enddo
       RunoffLeft_Efflux=0.0D0
       RunoffRight_Efflux=0.0D0 
       
cccz now set the new Efflux and the runoff
       do n=1,SurNodeIndex
        i=SurfNodeSurfIndexH(n)
        k=SurfNodeNodeIndexH(n)
        if(n.eq.1) then
         h_Pond_temp(n)=h_Pond_Old(n)+step*
     &    ((-q_Flux(n)+step*RunoffLeft)/width(i)+RunoffInten(n))
         h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
         if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz ask a question, why use h_Pond here. the answer is, it is on the edge
            h_Pond(n)=CriticalH_R_Sur(n)
         endif
        elseif(n.gt.1.and.n.lt.SurNodeIndex) then
         h_Pond_temp(n)=h_Pond_Old(n)+step*
     &     ((-q_Flux(n)+q_Flux(n-1))/width(i)+RunoffInten(n))
         h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
         if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz ask a question, why use h_Pond here
            h_Pond(n)=CriticalH_R_Sur(n)
         endif 
cccz ask a question, here for           
         CriticalH_R_Sur(n)=max(CriticalH_R_Sur(n),h_Pond_Old(n))
        else
         h_Pond_temp(n)=h_Pond_Old(n)+step*
     &    ((-RunoffRight+q_Flux(n-1))/width(i)+RunoffInten(n))
         h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
         if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            h_Pond(n)=CriticalH_R_Sur(n)
         endif 
        endif  
       enddo
         
cccz ask a question, seems we redo it here
cccz the answer is : because we revise the critical H
        do n=1,SurNodeIndex
         i=SurfNodeSurfIndexH(n)
         k=SurfNodeNodeIndexH(n)
         if(n.eq.1) then
          h_Pond_temp(n)=h_Pond_Old(n)+step*
     &      ((-q_Flux(n)+RunoffLeft)/width(i)+RunoffInten(n))
          h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
          if((slopeCoord(n,2)+h_Pond(n))
     &     .le.(slopeCoord(n+1,2)+h_Pond(n+1))) then
           if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            RunoffLeft_Efflux=RunoffLeft_Efflux
     &       -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
            h_Pond_temp(n)=CriticalH_R_Sur(n)
           endif
          else
           if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
             Efflux(n)=Efflux(n)+
     &         (h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
             h_Pond_temp(n)=CriticalH_R_Sur(n)
             RunoffLeft=0.0D0
             RunoffLeft_Effux=0.0D0
           endif 
          endif
         elseif(n.gt.1.and.n.lt.SurNodeIndex) then
           h_Pond_temp(n)=h_Pond_Old(n)+step*
     &       ((-q_Flux(n)+q_Flux(n-1))/width(i)+RunoffInten(n))
           h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
           if((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
c Previuos time step will not take care of q_Flux(n-1)
            if((slopeCoord(n,2)+h_Pond(n)).le.
     &        (slopeCoord(n+1,2)+h_Pond(n+1))) then
              if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               Efflux(n-1)=Efflux(n-1)
     &          -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=CriticalH_R_Sur(n)                   
              endif
            else
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
              Slope_0_n_1=abs(((slopeCoord(n-1,2)+h_Pond(n-1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n-1,1)))
              Slope_0_n=abs(((slopeCoord(n+1,2)+h_Pond(n+1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n+1,1)))
              Fraction_n_1=Slope_0_n_1/(Slope_0_n_1+Slope_0_n)
              Fraction_n=Slope_0_n/(Slope_0_n_1+Slope_0_n)

              Efflux(n-1)=Efflux(n-1)-Fraction_n_1*Add_Flux      
              Efflux(n)=Efflux(n)+Fraction_n*Add_Flux
              h_Pond_temp(n)=CriticalH_R_Sur(n)
             endif
            endif
                  
           elseif((slopeCoord(n,2)+h_Pond(n)).lt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
c The q_Flux(n-1) are taken cared by previous steps, then just calculate new q_Flux(n)
c if and only if q_Flux(n)>0
            if((slopeCoord(n,2)+h_Pond(n)).ge.
     &        (slopeCoord(n+1,2)+h_Pond(n+1))) then
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               Efflux(n)=Efflux(n)
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=CriticalH_R_Sur(n)
             endif
            else
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               q_Flux_Node(n)=
     &           (h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=CriticalH_R_Sur(n)
             endif
            endif

           else
            if((slopeCoord(n,2)+h_Pond(n)).lt.
     &         (slopeCoord(n+1,2)+h_Pond(n+1))) then
              if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
                Efflux(n-1)=Efflux(n-1)
     &            -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
                h_Pond_temp(n)=CriticalH_R_Sur(n)               
              endif
            elseif((slopeCoord(n,2)+h_Pond(n)).gt.
     &         (slopeCoord(n+1,2)+h_Pond(n+1))) then
              if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
                Efflux(n)=Efflux(n)
     &            +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
                h_Pond_temp(n)=CriticalH_R_Sur(n)
              endif
            else
              if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step   
              if (Efflux_old(n-1).gt.0.0D0) then
               if (Efflux_old(n).ge.0.0D0) then
                 Efflux(n)=Efflux(n)+Add_Flux
               else
                 if(abs(Efflux_old(n)).lt.Efflux_old(n-1)) then
                    Efflux(n)=Efflux(n)+Add_Flux
                 elseif(abs(Efflux_old(n)).gt.Efflux_old(n-1)) then
                    Efflux(n-1)=Efflux(n-1)-Add_Flux
                 else
                    Efflux(n-1)=Efflux(n-1)-0.50D0*Add_Flux
                    Efflux(n)=Efflux(n)+0.50D0*Add_Flux
                 endif
               endif
              elseif (Efflux_old(n-1).le.0.0D0) then
               if (Efflux_old(n).lt.0.0D0) then
                 Efflux(n-1)=Efflux(n-1)-Add_Flux
               else
                 Efflux(n-1)=Efflux(n-1)-0.50D0*Add_Flux
                 Efflux(n)=Efflux(n)+0.50D0*Add_Flux
               endif
              endif
              h_Pond_temp(n)=CriticalH_R_Sur(n)
              endif
            endif
           endif
          else
cccz the right-most point            
           h_Pond_temp(n)=h_Pond_Old(n)+step*
     &       ((-RunoffRight+q_Flux(n-1))/width(i)+RunoffInten(n))
           h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
           if((slopeCoord(n,2)+h_Pond(n)).le.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               RunoffRight_Efflux=RunoffRight_Efflux
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=CriticalH_R_Sur(n)
             endif
           else
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               Efflux(n-1)=Efflux(n-1)
     &           -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=CriticalH_R_Sur(n)
                        
               RunoffRight_temp=0.0D0
               RunoffRight_Efflux=0.0D0
             endif
           endif
          endif
         enddo
cccz update the ponded depth         
         do n=1,SurNodeIndex
           h_Pond(n)=max(h_Pond_temp(n),0.0D0)
         enddo
          
         do n=1,SurNodeIndex
          if (n.eq.1) then
            q_Flux_Node(n)=q_Flux_Node(n)
     &        -min(q_Flux(n)+Efflux(n),0.0D0)
     &        +max(RunoffLeft+RunoffLeft_Efflux,0.0D0)
          elseif (n.gt.1.and.n.lt.SurNodeIndex) then
            q_Flux_Node(n)=q_Flux_Node(n)
     &        +max(q_Flux(n-1)+Efflux(n-1),0.0D0)
     &        -min(q_Flux(n)+Efflux(n),0.0D0)
          else
            q_Flux_Node(n)=q_Flux_Node(n)
     &        +max(q_Flux(n-1)+Efflux(n-1),0.0D0)
     &        -min(RunoffRight+RunoffRight_Efflux,0.0D0)
          endif
         enddo
          
         do n=1,SurNodeIndex
          k=SurfNodeNodeIndexH(n)            ! for index in the whole node set
          i=SurfNodeSurfIndexH(n)            ! for index in the boundary node set
          VarBW(i,1)=Varbw_Mulch(i,1)+q_Flux_Node(n)/width(i) 
          VarBW(i,3)=Varbw_Mulch(i,2)-VarBW(i,1)                 ! no changes for varbw_mulch(k,2) or varbw(k,2), assumed equal
          Q(k)=-Width(i)*VarBW(i,3)                   
          If(Q(k).gt.0.0) then
           CodeW(k)=-4    
          endif
          if(Q(k).lt.Qact(k)) then
           Q_ave=max(Qact(k),0.0D0)
           h_Pond_temp(n)=max(h_Pond(n)-max(Q_ave-Q(k),0.0D0)
     &       *step/width(i),0.0D0)     
           VarBW(i,1)=VarBW(i,1)+(h_Pond(n)-h_Pond_temp(n))/step
           Q(k)=Q(k)+(h_Pond(n)-h_Pond_temp(n))/step*width(i)
           h_Pond(n)=max(h_Pond_temp(n),0.0D0)
           if(h_Pond(n).gt.CriticalH_R) then
              FurtherCheck_h_Pond(n)=1
           else
              FurtherCheck_h_Pond(n)=0
           endif
          else
           FurtherCheck_h_Pond(n)=0                 
          endif
cccz we need this furthercheck because we need to adjust 
cccz the ponded height based on the soil infiltrability
         enddo
          
c need to further check h_Pond(n)
c fix a good h_Pond value based on the topology of the domain
         iteration_Num=0.0D0
1106      do n=1,NumBPD
           Exchange_V_temp(n)=0.0D0
           Check_Exchange_V(n)=0
         enddo
          
         do n=2,SurNodeIndex-1
          if(FurtherCheck_h_Pond(n).eq.1) then
cccz         H_L=slopeCoord(n-1,2)+h_Pond_temp(n-1)
cccz         H_M=slopeCoord(n,2)+h_Pond_temp(n)
cccz         H_R=slopeCoord(n+1,2)+h_Pond_temp(n+1)
            H_L=slopeCoord(n-1,2)+h_Pond(n-1)
            H_M=slopeCoord(n,2)+h_Pond(n)
            H_R=slopeCoord(n+1,2)+h_Pond(n+1)
            if(n.eq.2) then
              L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
            else
              L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-2,1))
            endif
            if(n.eq.(SurNodeIndex-1)) then
              L_R=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
            else
              L_R=0.5*(slopeCoord(n+2,1)-slopeCoord(n,1))
            endif
            L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
cccz this is a valley shape
            if(H_M.lt.H_L.and.H_M.lt.H_R) then
              Ava_Height=(H_R*L_R+H_L*L_L+H_M*L_M)/(L_R+L_L+L_M)
              Volume=max((Ava_Height-
     &          (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
cccz check avaliable water from both sides
cccz ask a question, here try to use "min((h_Pond_temp(n-1)-CriticalH_R),(H_L-Ava_Height))" in future
c            Ava_Volume_L=(h_Pond_temp(n-1)-CriticalH_R)*L_L
c            Ava_Volume_R=(h_Pond_temp(n+1)-CriticalH_R)*L_R
            Ava_Volume_L=min((h_Pond_temp(n-1)-CriticalH_R),
     &         (H_L-Ava_Height))*L_L
            Ava_Volume_R=min((h_Pond_temp(n+1)-CriticalH_R),
     &         (H_R-Ava_Height))*L_R
              Ava_Volume=Ava_Volume_L+Ava_Volume_R
                  
              if(Ava_Volume.ge.Volume) then
                if(Ava_Volume_L.le.(Volume/2.0D0)) then
                  Volume_L=Ava_Volume_L
                  Volume_R=Volume-Volume_L
                elseif(Ava_Volume_R.le.(Volume/2.0D0)) then
                  Volume_R=Ava_Volume_R
                  Volume_L=Volume-Volume_R
                else
                  Volume_L=Volume/2.0D0
                  Volume_R=Volume/2.0D0
                endif
              else
               Volume_L=Ava_Volume_L
               Volume_R=Ava_Volume_R             
              endif
              if(Check_Exchange_V(n-1).eq.0) then
                Exchange_V_temp(n-1)=Volume_L
                Check_Exchange_V(n-1)=1
              else
               if(Exchange_V_temp(n-1).lt.0.0D0) then
                 Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
               else
                 Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
               endif
              endif
              if(Check_Exchange_V(n).eq.0) then
                Exchange_V_temp(n)=-Volume_R
                Check_Exchange_V(n)=1
              else
               if(Exchange_V_temp(n).lt.0.0D0) then
                 Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)
               else
                 Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)
               endif
              endif
cccz leftwards flow case     
             elseif(H_M.ge.H_L.and.H_M.le.H_R) then
               Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
               Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
               Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)
cccz always choose the higher value
               if(H_M.ge.Mean_Height) then
                 Ava_Height_L=Mean_Height
                 Ava_Height_R=Mean_Height_R
               else
                 Ava_Height_L=Mean_Height_L
                 Ava_Height_R=Mean_Height
               endif
     
               Volume_L=max((Ava_Height_L-
     &           (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
               Volume_L=min(Volume_L,
     &           max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_L,
     &           h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
     
               Volume_R=max((Ava_Height_R-
     &           (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
               Volume_R=min(Volume_R,
     &          max(min(slopeCoord(n+1,2)+h_Pond_temp(n+1)-Ava_Height_R,
     &          h_Pond_temp(n+1)-CriticalH_R)*L_R,0.0D0))
        
              if(Check_Exchange_V(n-1).eq.0) then
                Exchange_V_temp(n-1)=-Volume_L
                Check_Exchange_V(n-1)=1
              else
               if(Exchange_V_temp(n-1).lt.0.0D0) then
                Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)       
               else
                Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L)
               endif
              endif
              if(Check_Exchange_V(n).eq.0) then
                Exchange_V_temp(n)=-Volume_R
                Check_Exchange_V(n)=1
              else
               if(Exchange_V_temp(n).lt.0.0D0) then
                Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)
               else
                Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)
               endif
              endif
cccz rightwards flow case
            elseif(H_M.le.H_L.and.H_M.ge.H_R) then
              Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
              Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
              Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)
cccz always choose the higher value
              if(H_M.ge.Mean_Height) then
                Ava_Height_R=Mean_Height
                Ava_Height_L=Mean_Height_L
              else
                Ava_Height_R=Mean_Height_R
                Ava_Height_L=Mean_Height
              endif
              Volume_R=max((Ava_Height_R-
     &          (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
              Volume_R=min(Volume_R,
     &          max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_R,
     &          h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
                  
              Volume_L=max((Ava_Height_L-
     &          (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
              Volume_L=min(Volume_L,
     &          max(min(slopeCoord(n-1,2)+h_Pond_temp(n-1)-Ava_Height_L,
     &          h_Pond_temp(n-1)-CriticalH_R)*L_L,0.0D0))
                       
              if(Check_Exchange_V(n-1).eq.0) then
                Exchange_V_temp(n-1)=Volume_L
                Check_Exchange_V(n-1)=1
              else
               if(Exchange_V_temp(n-1).lt.0.0D0) then
                Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
               else
                Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
               endif
              endif

              if(Check_Exchange_V(n).eq.0) then
                Exchange_V_temp(n)=Volume_R
                Check_Exchange_V(n)=1
              else
                if(Exchange_V_temp(n).lt.0.0D0) then
                 Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
                else
                 Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
                endif
              endif       
cccz flat/peak case              
            else
              Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
              Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
              Mean_Height=(H_L*L_L+H_M*L_M+H_R*L_R)/(L_L+L_M+L_R)

              if(Mean_Height_L.lt.Mean_Height) then
                Ava_Height_L=Mean_Height
                Ava_Height_R=Mean_Height_R
              else
                Ava_Height_R=Mean_Height
                Ava_Height_L=Mean_Height_L
              endif
cccz water requirement                  
              Volume_L=max((Ava_Height_L-
     &          (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
              Volume_R=max((Ava_Height_R-
     &          (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
              Volume_T=Volume_L+Volume_R
cccz water supply          
              Volume=min(Volume_T,
     &          max(min(slopeCoord(n,2)+h_Pond_temp(n)-
     &          min(Mean_Height_L,Mean_Height_R),      
     &          h_Pond_temp(n)-CriticalH_R)*L_L,0.0D0))      
cccz water redistributed     
              if(Volume_T.eq.0.0D0) then
              else
                Volume_L=Volume*Volume_L/Volume_T
                Volume_R=Volume*Volume_R/Volume_T
              endif                      
              if(Check_Exchange_V(n-1).eq.0) then
                Exchange_V_temp(n-1)=-Volume_L
                Check_Exchange_V(n-1)=1
              else
               if(Exchange_V_temp(n-1).lt.0.0D0) then
                Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)
               else
                Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L)
               endif
              endif
              if(Check_Exchange_V(n).eq.0) then
                Exchange_V_temp(n)=Volume_R
                Check_Exchange_V(n)=1
              else
               if(Exchange_V_temp(n).lt.0.0D0) then
                Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
               else
                Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
               endif
              endif                
            endif   
          endif
         enddo

cccz the water right at the edge had to been drained.
cccz moreover, the edge points should not provide inwards water flux since it is already below CriticalHR
         Exchange_V_temp(1)=min(Exchange_V_temp(1),0.0D0)
         Exchange_V_temp(n)=max(Exchange_V_temp(n),0.0D0)
         error=0.0D0
              
         do n=1,SurNodeIndex
          if(n.eq.1) then
            L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
          elseif(n.eq.SurNodeIndex) then
            L_M=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
          else
            L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
          endif
          if(n.eq.1) then
            h_Pond_temp(n)=h_Pond_temp(n)-Exchange_V_temp(n)/L_M
cccz if "gt", the Exchange_V_temp(1)<0, and leftwards exchange occurs
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              RunoffLeft_Efflux=RunoffLeft_Efflux
     &          -(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
              h_Pond_temp(n)=CriticalH_R_Sur(n)
            endif
          elseif(n.eq.SurNodeIndex) then
            h_Pond_temp(n)=h_Pond_temp(n)+Exchange_V_temp(n-1)/L_M
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              RunoffRight_Efflux=RunoffRight_Efflux
     &          +(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
              h_Pond_temp(n)=CriticalH_R_Sur(n)
            endif
          else
            h_Pond_temp(n)=h_Pond_temp(n)+
     &        (Exchange_V_temp(n-1)-Exchange_V_temp(n))/L_M
          endif
          error=max(error,abs(h_Pond_temp(n)-h_Pond(n)))
         enddo
          
         iteration_Num=iteration_Num+1
          
         if((error.gt.0.01D0*CriticalH_R).and.(iteration_Num.lt.100)) 
     &       then
           do n=1,SurNodeIndex  
cccz              h_Pond(n)=(h_Pond_temp(n)+h_Pond(n)*(iteration_Num-1))
cccz     &               /iteration_Num
             h_Pond(n)=h_Pond_temp(n) 
             Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
           enddo
           h_Pond(1)=CriticalH_R_Sur(1)
           h_Pond(SurNodeIndex)=CriticalH_R_Sur(SurNodeIndex)
           goto 1106
         else
           do n=1,SurNodeIndex
             k=SurfNodeNodeIndexH(n)
             h_Pond(n)=h_Pond_temp(n)
             h_Pond(n)=max(h_Pond_temp(n),0.0D0)
             Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
             if(h_Pond(n).gt.CriticalH_R) then
cccz assign the boundary condition for infiltration
               hNew(k)=CriticalH+h_Pond(n)                 
             endif
           enddo
           endif
           
cccz now we should have a overall idea of the three types of water flow, all in unit (cm^2/day)
cccz "q_flux" is the most non-free flux which follow the SV equation
cccz "exchange_V" is the intermediately free for water level adjustment
cccz "efflux" is the very free flux
           
cccz finally, change the "exchange_V" (a volume) to a flux term
        do i=1,SurNodeIndex
           Exflux(n)=Exchange_V(n)/step
        enddo
        goto 1300


cccz ----------------------------------------------------------------------
c -------------------------------------------------------------------------
cccz start the flat surface one
cccz need to evaluate criticalH_R for each node based on the hNew, 
cccz for totally flat surface, the process can be very simple
 
cccz determine the criticalH_R while adding a small slope for the water surface
cccz first calculate the hori distance
1200  totalLength=abs(slopeCoord(SurNodeIndex,1)-slopeCoord(1,1))
      do n=1,SurNodeIndex
        CriticalH_R_Sur(n)=CriticalH_R
     &   +(1.0D0-abs(slopeCoord(n,1)-slopeCoord(1,1))/totalLength)
     &   *0.000001D0
      enddo
      
      iteration_Num=0          
1201  iteration_Flux=0
      iteration_Head=0
      RunoffLeft_temp=0.0D0
      RunoffRight_temp=0.0D0
      do i=1,NumBPD
         h_Pond_temp(i)=0.0D0
         q_Flux_temp(i)=0.0D0
         q_Flux_Node(i)=0.0D0
         Efflux(i)=0.0D0
      enddo
      RunoffLeft_Efflux=0.0D0
      RunoffRight_Efflux=0.0D0

          
      do n=1,SurNodeIndex-1
cccz use iteration to partition the real runoff and amount of water stay 
cccz take the absolute value of the slope to ensure the following calculation
        Slope_0_R=abs((slopeCoord(n+1,2)-slopeCoord(n,2))
     &      /(slopeCoord(n,1)-slopeCoord(n+1,1)))
cccz start the explicit version for the q_flux
cccz we also check the right/left directions and made average for stability
        if(h_Pond(n+1).eq.0) then
           Slope_f_left=0.0D0
        else
           Slope_f_left=((n_Stiff*q_Flux(n))**2.0D0)
     &       /(h_Pond(n+1)**(3.333D0))
        endif
        if(h_Pond(n).eq.0) then
           Slope_f_right=0.0D0
        else
           Slope_f_right=((n_Stiff*q_Flux(n)**2.0D0)
     &       /(h_Pond(n)**(3.333D0)))
        endif
        Slope_f=0.50D0*(Slope_f_left+Slope_f_right)

cccz now we calcualte the fluxes
cccz start from the first (left) node
        if(n.eq.1) then
cccz Water will flow from node 1 to node 2
cccz Always follow the upwind direction
         if((slopeCoord(1,2)+h_Pond(1)).ge.
     &    (slopeCoord(2,2)+h_Pond(2))) then
           if(h_Pond(1).eq.0.0D0) then
            q_flux_square=0.0D0
            h_flux_square=10.0D0       ! arbitrary none-zero number here
           else
            q_flux_square=q_Flux(1)
            h_flux_square=h_Pond(1)
           endif
cccz we can see for values labeled "_Old" is used to store old data
cccz which is not changed during the iteration
         
           q_Flux_temp(1)=q_Flux_Old(1)+
     &         (-(q_flux_square**2.0D0)/h_flux_square/ll(1)             ! there is no flux on the left of left edge
     &         +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1)   
     &         +g_accel*h_Pond(1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(1)=max(q_Flux_temp(1),0.0D0)
            
cccz the direction implies runoff will not be discharged from the left side
            RunoffLeft_temp=0.0D0
           else
cccz Water will flow from Node 2 to Node 1, momentum source from q_flux(2), between Node 2 and Node 3
cccz this depends on whether Node 2 is a local maximum                   
            if ((slopeCoord(2,2)+h_Pond(2)).ge.
     &        (slopeCoord(3,2)+h_Pond(3))) then  
cccz now we calcualte water flow     
             if(h_Pond(2).eq.0.0D0) then
               q_flux_square=0.0D0
               h_flux_square=10.0D0                                   ! arbitrary number here
             else
               q_flux_square=q_Flux(1)
               h_flux_square=h_Pond(2)
             endif       
             q_Flux_temp(1)=q_Flux_Old(1)+
     &          ((q_flux_square**2.0D0)/h_flux_square/ll(1)              ! there is no income momentum from q_flux(2)
     &          +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1)
     &          -g_accel*h_Pond(2)*max(Slope_0_R-Slope_f,0.0D0))*step
             q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
             
            else
cccz in this case, there is momentum source from q_flux(2)                 
             if(h_Pond(2).eq.0.0D0) then
               q_flux_square1=0.0D0
               h_flux_square1=10.0D0                                  ! arbitrary number here
             else
               q_flux_square1=q_Flux(1)
               h_flux_square1=h_Pond(2)
             endif
             if(h_Pond(3).eq.0.0D0) then
               q_flux_square2=0.0D0
               h_flux_square2=10.0D0                                  ! arbitrary number here
             else
               q_flux_square2=q_Flux(2)
               h_flux_square2=h_Pond(3)
             endif
                     
             q_Flux_temp(1)=q_Flux_Old(1)+
     &         (((q_flux_square1**2.0D0)/h_flux_square1
     &              -(q_flux_square2**2.0D0)/h_flux_square2)/ll(1)
     &         +0.5D0*g_accel*(h_Pond(1)**2.0D0-h_Pond(2)**2.0D0)/ll(1)           
     &         -g_accel*h_Pond(2)*max(Slope_0_R-Slope_f,0.0D0))*step         
             q_Flux_temp(1)=min(q_Flux_temp(1),0.0D0)
            endif

cccz finish the calculation of q in the case Node 2 --> Node 1
cccz then calculate the runoff from the left edge
cccz in this case we assume that the leftside is blocked for runoff, 
cccz i.e., no leftrunoff discharge      

            RunoffLeft_temp=0.0D0
          endif

cccz finished the calculation when n=1 (left-most) edge
cccz now calculate the flux in interior nodes
         elseif(n.gt.1.and.n.lt.SurNodeIndex-1) then
     
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
cccz water will flow from Node n -> n+1, right-bound, momentun source may from n-1?
cccz always follow the upwind direction 
cccz in the first case, Node n is the local peak s.t. no momentun from the left side 
           if((slopeCoord(n,2)+h_Pond(n)).ge.
     &        (slopeCoord(n-1,2)+h_Pond(n-1))) then   
cccz the momentum is only from the current slope section, i.e., n
            if(h_Pond(n).eq.0.0D0) then
              q_flux_square=0.0D0
              h_flux_square=10.0D0                                    ! arbitrary none-zero number here
            else
              q_flux_square=q_Flux(n)
              h_flux_square=h_Pond(n)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        (-(q_flux_square**2.0D0)/h_flux_square/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
           
           else
cccz the momentum is from the left side, i.e., n-1                  
            if(h_Pond(n).eq.0.0D0) then
              q_flux_square1=0.0D0
              h_flux_square1=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square1=q_Flux(n)
              h_flux_square1=h_Pond(n)
            endif
            if(h_Pond(n-1).eq.0.0D0) then
              q_flux_square2=0.0D0
              h_flux_square2=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square2=q_Flux(n-1)
              h_flux_square2=h_Pond(n-1)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        ((-(q_flux_square1**2.0D0)/h_flux_square1
     &           +(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step  
            q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
           endif
         else
cccz Water will flow from Node n+1 -> n, momentun source of q_Flux(n+1)?
cccz always follow the upwind direction 
cccz now we assume n+1 is a local maximum, so no water flow Node n+2->n+1
cccz thus, there is no external momentum from q_Flux(n+1), which should be in an opposite direction
           if((slopeCoord(n+1,2)+h_Pond(n+1)).ge.
     &        (slopeCoord(n+2,2)+h_Pond(n+2))) then
             if(h_Pond(n+1).eq.0.0D0) then
               q_flux_square=0.0D0
               h_flux_square=10.0D0                                    ! arbitrary none-zero number here
             else
               q_flux_square=q_Flux(n)
               h_flux_square=h_Pond(n+1)
             endif
                                         
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        ((q_flux_square**2.0D0)/h_flux_square/ll(n)
     &        +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &        -g_accel*h_Pond(n+1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
           else    
cccz now we assume n+2 is higher than n+1
cccz thus, there is external momentum from q_Flux(n+1), from Node n+2
            if(h_Pond(n+1).eq.0.0D0) then
              q_flux_square1=0.0D0
              h_flux_square1=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square1=q_Flux(n)
              h_flux_square1=h_Pond(n+1)
            endif
            if(h_Pond(n+2).eq.0.0D0) then
              q_flux_square2=0.0D0
              h_flux_square2=10.0D0                                   ! arbitrary none-zero number here
            else
              q_flux_square2=q_Flux(n+1)
              h_flux_square2=h_Pond(n+2)
            endif
                  
            q_Flux_temp(n)=q_Flux_Old(n)+
     &        (((q_flux_square1**2.0D0)/h_flux_square1
     &           -(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &       -g_accel*h_Pond(n+1)*max(Slope_0_R-Slope_f,0.0D0))*step
            q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
                     
          endif
         endif

cccz now we discuss the right-most point, with flux q_flux(n=SurNodeIndex-1)
        else
       
cccz Finish the calculation for n=SurNodeIndex-1
cccz or in other words, we have "n=SurNodeIndex-1" here, so you will not see "n+2"
         if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
cccz water flow to the end point n=SurNodeIndex, need to calculate the momentun from n-2?
cccz n=SurNodeIndex-1 is the local maxima, no water flux from left side 
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then 
           if(h_Pond(n).eq.0.0D0) then
             q_flux_square=0.0D0
             h_flux_square=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square=q_Flux(n)
             h_flux_square=h_Pond(n)
           endif
           q_Flux_temp(n)=q_Flux_Old(n)+
     &       (-(q_flux_square**2.0D0)/h_flux_square/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &       +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
           q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
          else
cccz n=SurNodeIndex-2 is higher, water flux from left side                   
           if(h_Pond(n).eq.0.0D0) then
             q_flux_square1=0.0D0
             h_flux_square1=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square1=q_Flux(n)
             h_flux_square1=h_Pond(n)
           endif    
           if(h_Pond(n-1).eq.0.0D0) then
             q_flux_square2=0.0D0
             h_flux_square2=10.0D0                                     ! arbitrary none-zero number here
           else
             q_flux_square2=q_Flux(n-1)
             h_flux_square2=h_Pond(n-1)
           endif
                  
           q_Flux_temp(n)=q_Flux_Old(n)+
     &       ((-(q_flux_square1**2.0D0)/h_flux_square1
     &           +(q_flux_square2**2.0D0)/h_flux_square2)/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0-h_Pond(n+1)**2.0D0)/ll(n) 
     &       +g_accel*h_Pond(n)*max(Slope_0_R-Slope_f,0.0D0))*step
           q_Flux_temp(n)=max(q_Flux_temp(n),0.0D0)
          endif

cccz in this case, runoff discharge on the right hand side occurs
cccz this calculation should be processed based on the SV equation
          if(h_Pond(n+1).eq.0.0D0) then
            q_flux_runoffright=0.0D0
            h_flux_runoffright=10.0D0       ! arbitrary number here
          else
            q_flux_runoffright=RunoffRight
            h_flux_runoffright=h_Pond(SurNodeIndex)
          endif    
          if(h_Pond(n).eq.0.0D0) then
            q_flux_square=0.0D0
            h_flux_square=10.0D0       ! arbitrary number here
          else
            q_flux_square=q_Flux(n)
            h_flux_square=h_Pond(n)
          endif                  

          RunoffRight_temp=RunoffRight_old+
     &      ((-(q_flux_runoffright**2.0D0)/h_flux_runoffright
     &         +(q_flux_square**2.0D0)/h_flux_square)
     &         /width(SurfNodeSurfIndexH(SurNodeIndex))
     &      +0.5D0*g_accel*(h_Pond(n)**2.0D0
     &         -h_Pond(SurNodeIndex)**2.0D0)/ll(n) 
     &      +g_accel*h_Pond(SurNodeIndex)*max(Slope_0_R-Slope_f,0.0D0))
     &         *step 
          RunoffRight_temp=max(RunoffRight_temp,0.0D0)
                             
         else

cccz water flow to from Node SurNodeIndex -> SurNodeIndex-1
cccz in this case, no right runoff discharge
cccz also, a good news is there is no right-hand node, thus, no external momentum source

          if(h_Pond(n+1).eq.0.0D0) then
            q_flux_square=0.0D0
            h_flux_square=10.0D0       ! arbitrary number here
          else
            q_flux_square=q_Flux(n)
            h_flux_square=h_Pond(SurNodeIndex)
          endif

          q_Flux_temp(n)=q_Flux_Old(n)+
     &       ((q_flux_square**2.0D0)/h_flux_square/ll(n)
     &       +0.5D0*g_accel*(h_Pond(n)**2.0D0
     &          -h_Pond(SurNodeIndex)**2.0D0)/ll(n) 
     &       -g_accel*h_Pond(SurNodeIndex)*max(Slope_0_R-Slope_f,0.0D0))
     &          *step
          q_Flux_temp(n)=min(q_Flux_temp(n),0.0D0)
          RunoffRight_temp=0.0D0
        endif
       endif
      enddo
        
cccz finishe the flux estimation, now go the ponded water head estimation
cccz The following steps will make the h_temp goes infinitely high
cccz so after solving the SV equation, we need additional processes

      do n=1,SurNodeIndex
        i=SurfNodeSurfIndexH(n)
        k=SurfNodeNodeIndexH(n)
        if(n.eq.1) then
         h_Pond_temp(1)=h_Pond_Old(1)+step*
     &     ((-q_Flux_temp(1)+RunoffLeft_temp)/width(i)+RunoffInten(1))
         h_Pond_temp(1)=max(h_Pond_temp(1),0.0D0)

cccz We do not allow the calculation for the left efflux        
        if(h_Pond_temp(1).gt.CriticalH_R_Sur(1)) then
          Efflux(1)=Efflux(1)
     &      +(h_Pond_temp(1)-CriticalH_R_Sur(1))*width(i)/step
          h_Pond_temp(1)=min(h_Pond_temp(1),CriticalH_R_Sur(1))
          RunoffLeft_temp=0.0D0
          RunoffLeft_Efflux=0.0D0
        endif 
         
cccz for the interior nodes, the efflux could be left/right/outwards/inwards
        elseif(n.gt.1.and.n.lt.SurNodeIndex) then
          h_Pond_temp(n)=h_Pond_Old(n)+step*
     &      ((-q_Flux_temp(n)+q_Flux_temp(n-1))/width(i)+RunoffInten(n))
          h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)     

          if((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
cccz Previuos time step will not take care of EFFlux(n-1)
           if((slopeCoord(n,2)+h_Pond(n)).le.
     &        (slopeCoord(n+1,2)+h_Pond(n+1))) then

cccz in this geometry, water will flow left in "n"             
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               Efflux(n-1)=Efflux(n-1)
     &           -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n)) 
             endif

cccz in this geometry, water will flow outwards like "/\" 
cccz but how the runoff partitioned to left/right sides
           else
             if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz first calculate all avail runoff water
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
cccz calculate left/right slope
              Slope_0_n_1=abs(((slopeCoord(n-1,2)+h_Pond(n-1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n-1,1)))
              Slope_0_n=abs(((slopeCoord(n+1,2)+h_Pond(n+1))
     &           -(slopeCoord(n,2)+h_Pond(n)))
     &           /(slopeCoord(n,1)-slopeCoord(n+1,1)))
cccz use the slope to determine the fraction
              Fraction_n_1=Slope_0_n_1/(Slope_0_n_1+Slope_0_n)
              Fraction_n=Slope_0_n/(Slope_0_n_1+Slope_0_n)
cccz calculate the Efflux based on the fraction
              Efflux(n-1)=Efflux(n-1)-Fraction_n_1*Add_Flux
              Efflux(n)=Efflux(n)+Fraction_n*Add_Flux
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
             endif
           endif
         elseif((slopeCoord(n,2)+h_Pond(n)).lt.
     &     (slopeCoord(n-1,2)+h_Pond(n-1))) then
cccz The q_Flux(n-1) are taken cared by previous steps, then just calculate new q_Flux(n)
cccz if and only if q_Flux(n)>0
cccz i.e., calculate rightwards efflux on the rightside
          if((slopeCoord(n,2)+h_Pond(n)).ge.
     &      (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Efflux(n)=Efflux(n)
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
            endif
cccz here we deal with a special case -- a valley
cccz if there is a valley, the efflux water will not go anywhere,
cccz but we should record it and be prepared to add it as water input for the next iteration
          else
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              q_Flux_Node(n)=(h_Pond_temp(n)-CriticalH_R_Sur(n))
     &           *width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
            endif
          endif

cccz now the elevation of the water level @ Node n and Node n-1
cccz was the same.

         else
           if((slopeCoord(n,2)+h_Pond(n)).lt.
     &       (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz the rightside is higher, momentum points leftwards
              Efflux(n-1)=Efflux(n-1)
     &           -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))   
            endif
           elseif((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n+1,2)+h_Pond(n+1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
cccz the rightside is lower, momentum points rightwards
              Efflux(n)=Efflux(n)
     &          +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step   
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n)) 
            endif
           else
cccz totally flat  
cccz first record the potential efflux quantity
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Add_Flux=(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
cccz need to refer the momentum from previous step
              if(Efflux_old(n-1).gt.0.0D0) then
               if(Efflux_old(n).ge.0.0D0) then
                 Efflux(n)=Efflux(n)+Add_Flux
               else
cccz follow the side with larger momentum
                if(abs(Efflux_old(n)).lt.Efflux_old(n-1)) then
                  Efflux(n)=Efflux(n)+Add_Flux
                elseif(abs(Efflux_old(n)).gt.Efflux_old(n-1)) then
                  Efflux(n-1)=Efflux(n-1)-Add_Flux
                else
                  Efflux(n-1)=Efflux(n-1)-0.50D0*Add_Flux
                  Efflux(n)=Efflux(n)+0.50D0*Add_Flux
                endif
               endif
              elseif (Efflux_old(n-1).le.0.0D0) then
               if (Efflux_old(n).lt.0.0D0) then
                  Efflux(n-1)=Efflux(n-1)-Add_Flux
               else

                  Efflux(n-1)=Efflux(n-1)-Add_Flux/3.0D0
                  Efflux(n)=Efflux(n)+Add_Flux/3.0D0
                  q_Flux_Node(n)=q_Flux_Node(n)+Add_Flux/3.0D0
               endif
              endif
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
             endif
           endif
          endif
                
cccz finally, we calcualte the h_Pond for the right-most point
         else
          h_Pond_temp(n)=h_Pond_Old(n)+step*
     &    ((-RunoffRight_temp+q_Flux_temp(n-1))/width(i)+RunoffInten(n))
          h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
cccz right runoff discharge case
          if((slopeCoord(n,2)+h_Pond(n)).le.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
               RunoffRight_Efflux=RunoffRight_Efflux
     &           +(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
               h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))    
            endif
          elseif((slopeCoord(n,2)+h_Pond(n)).gt.
     &       (slopeCoord(n-1,2)+h_Pond(n-1))) then
            if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
              Efflux(n-1)=Efflux(n-1)
     &          -(h_Pond_temp(n)-CriticalH_R_Sur(n))*width(i)/step
              h_Pond_temp(n)=min(h_Pond_temp(n),CriticalH_R_Sur(n))
              RunoffRight_temp=0.0D0
              RunoffRight_Efflux=0.0D0
            endif
          else              
         endif
        endif
       enddo

cccz finish update of head, need further work for setting the critical_H_Sur value for each node.
cccz Error esitimation for iteration, i.e.,
cccz the Picard's for SV equation
          
       iteration_Flux=0
       iteration_Head=0
       iteration_Num=iteration_Num+1
cccz   iteration_Num=1.0D0
cccz take average between the current iteration and previous iteration steps to seek for stability
       do n=1,SurNodeIndex
         q_Flux_temp(n)=q_Flux_temp(n)/iteration_Num+
     &     (iteration_Num-1)*q_Flux(n)/iteration_Num
         h_Pond_temp(n)=h_Pond_temp(n)/iteration_Num+
     &     (iteration_Num-1)*h_Pond(n)/iteration_Num
         h_Pond_temp(n)=max(h_Pond_temp(n),0.0D0)
       enddo
       RunoffRight_temp=RunoffRight_temp/iteration_Num+
     &   (iteration_Num-1)*RunoffRight/iteration_Num
       RunoffLeft_temp=RunoffLeft_temp/iteration_Num+
     &   (iteration_Num-1)*RunoffLeft/iteration_Num
          
       do n=1,SurNodeIndex
        if(abs(q_Flux_temp(n)-q_Flux(n)).gt.0.01*abs(q_Flux(n))) then
          if(abs(q_Flux_temp(n)).gt.FluxLimit)  iteration_Flux=1
        endif
        if(abs(h_Pond_temp(n)-h_Pond(n)).gt.0.01*abs(h_Pond(n))) then
          if(abs(h_Pond_temp(n)).gt.HeadLimit)  iteration_Head=1
        endif
       enddo
          
       if(abs(RunoffLeft_temp-RunoffLeft).gt.0.01*abs(RunoffLeft)) then
         if(abs(RunoffLeft).gt.FluxLimit) iteration_Flux=1
       endif
       if(abs(RunoffRight_temp-RunoffRight).gt.0.01*abs(RunoffRight)) 
     &  then            
         if(abs(RunoffRight).gt.FluxLimit) iteration_Flux=1
       endif

       if (iteration_Flux.eq.1.or.iteration_Head.eq.1) then
        do n=1,SurNodeIndex
          q_Flux(n)=q_Flux_temp(n)
          h_Pond(n)=h_Pond_temp(n)
        enddo
          RunoffLeft=RunoffLeft_temp
          RunoffRight=RunoffRight_temp
        goto 1201
       else
        do n=1,SurNodeIndex
          q_Flux(n)=q_Flux_temp(n)
          h_Pond(n)=h_Pond_temp(n)
          FurtherCheck_h_Pond(n)=1
        enddo
        RunoffLeft=RunoffLeft_temp
        RunoffRight=RunoffRight_temp
       endif

cccz Now we already solve the SV equation
cccz assign the flux as water input at each node
cccz for regular grid, the program can be ended here, after some arrangement
          
       do n=1,SurNodeIndex
        if (n.eq.1) then
          q_Flux_Node(n)=q_Flux_Node(n)
     &      -min(q_Flux(n)+Efflux(n),0.0D0)
     &      +max(RunoffLeft+RunoffLeft_Efflux,0.0D0)
        elseif (n.gt.1.and.n.lt.SurNodeIndex) then
          q_Flux_Node(n)=q_Flux_Node(n)
     &      +max(q_Flux(n-1)+Efflux(n-1),0.0D0)
     &      -min(q_Flux(n)+Efflux(n),0.0D0)
        else
          q_Flux_Node(n)=q_Flux_Node(n)
     &      +max(q_Flux(n-1)+Efflux(n-1),0.0D0)
     &      -min(RunoffRight+RunoffRight_Efflux,0.0D0)
        endif
       enddo
          
       do n=1,SurNodeIndex
        k=SurfNodeNodeIndexH(n)            ! for index in the whole node set
        i=SurfNodeSurfIndexH(n)            ! for index in the boundary node set
        VarBW(i,1)=Varbw_Mulch(i,1)+q_Flux_Node(n)/width(i) 
        VarBW(i,3)=Varbw_Mulch(i,2)-VarBW(i,1)                 
        Q(k)=-Width(i)*VarBW(i,3)                   
        If(Q(k).gt.0.0) then
         CodeW(k)=-4    
        endif
        if(Q(k).lt.Qact(k)) then
         Q_ave=max(Qact(k),0.0D0)
         h_Pond_temp(n)=max(h_Pond(n)-max(Q_ave-Q(k),0.0D0)
     &     *step/width(i),0.0D0)
         VarBW(i,1)=VarBW(i,1)+(h_Pond(n)-h_Pond_temp(n))/step
         Q(k)=Q(k)+(h_Pond(n)-h_Pond_temp(n))/step*width(i)
         h_Pond(n)=max(h_Pond_temp(n),0.0D0)
         if(h_Pond(n).gt.CriticalH_R) then
            FurtherCheck_h_Pond(n)=1
         else
            FurtherCheck_h_Pond(n)=0
         endif
        else
         FurtherCheck_h_Pond(n)=0                 
        endif
cccz we need this furthercheck because we need to adjust 
cccz the ponded height based on the soil infiltrability

       enddo
          
c need to further check h_Pond(n)
c fix a good h_Pond value based on the topology of the domain
       iteration_Num=0.0D0
1206   do n=1,NumBPD
        Exchange_V_temp(n)=0.0D0
        Check_Exchange_V(n)=0
       enddo
          
       do n=2,SurNodeIndex-1
        if(FurtherCheck_h_Pond(n).eq.1) then
          H_L=slopeCoord(n-1,2)+h_Pond(n-1)
          H_M=slopeCoord(n,2)+h_Pond(n)
          H_R=slopeCoord(n+1,2)+h_Pond(n+1)
          if(n.eq.2) then
            L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
          else
            L_L=0.5*(slopeCoord(n,1)-slopeCoord(n-2,1))
          endif
          if(n.eq.(SurNodeIndex-1)) then
            L_R=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
          else
            L_R=0.5*(slopeCoord(n+2,1)-slopeCoord(n,1))
          endif
          L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
cccz this is a valley shape
          if(H_M.lt.H_L.and.H_M.lt.H_R) then
           Ava_Height=(H_R*L_R+H_L*L_L+H_M*L_M)/(L_R+L_L+L_M)
           Volume=max((Ava_Height-
     &      (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
cccz check avaliable water from both sides
cccz ask a question, here try to use "min((h_Pond_temp(n-1)-CriticalH_R),(H_L-Ava_Height))" in future
c            Ava_Volume_L=(h_Pond_temp(n-1)-CriticalH_R)*L_L
c            Ava_Volume_R=(h_Pond_temp(n+1)-CriticalH_R)*L_R
           Ava_Volume_L=min((h_Pond_temp(n-1)-CriticalH_R),
     &      (H_L-Ava_Height))*L_L
           Ava_Volume_R=min((h_Pond_temp(n+1)-CriticalH_R),
     &      (H_R-Ava_Height))*L_R
           Ava_Volume=Ava_Volume_L+Ava_Volume_R
                  
           if(Ava_Volume.ge.Volume) then
            if(Ava_Volume_L.le.(Volume/2.0D0)) then
              Volume_L=Ava_Volume_L
              Volume_R=Volume-Volume_L
            elseif(Ava_Volume_R.le.(Volume/2.0D0)) then
              Volume_R=Ava_Volume_R
              Volume_L=Volume-Volume_R
            else
              Volume_L=Volume/2.0D0
              Volume_R=Volume/2.0D0
            endif
           else
            Volume_L=Ava_Volume_L
            Volume_R=Ava_Volume_R             
           endif
           if(Check_Exchange_V(n-1).eq.0) then
             Exchange_V_temp(n-1)=Volume_L
             Check_Exchange_V(n-1)=1
           else
            if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
            else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
            endif
           endif
           if(Check_Exchange_V(n).eq.0) then
             Exchange_V_temp(n)=-Volume_R
             Check_Exchange_V(n)=1
           else
            if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)
            else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)
            endif
           endif
cccz leftwards flow case     
          elseif(H_M.ge.H_L.and.H_M.le.H_R) then
            Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
            Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
            Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)
cccz always choose the higher value
           if(H_M.ge.Mean_Height) then
             Ava_Height_L=Mean_Height
             Ava_Height_R=Mean_Height_R
           else
             Ava_Height_L=Mean_Height_L
             Ava_Height_R=Mean_Height
           endif
     
           Volume_L=max((Ava_Height_L-
     &       (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
           Volume_L=min(Volume_L,
     &       max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_L,
     &       h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
     
           Volume_R=max((Ava_Height_R-
     &       (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
           Volume_R=min(Volume_R,
     &       max(min(slopeCoord(n+1,2)+h_Pond_temp(n+1)-Ava_Height_R,
     &       h_Pond_temp(n+1)-CriticalH_R)*L_R,0.0D0))
        
           if(Check_Exchange_V(n-1).eq.0) then
             Exchange_V_temp(n-1)=-Volume_L
             Check_Exchange_V(n-1)=1
           else
            if(Exchange_V_temp(n-1).lt.0.0D0) then
             Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)       
            else
             Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L)
            endif
           endif
           if(Check_Exchange_V(n).eq.0) then
             Exchange_V_temp(n)=-Volume_R
             Check_Exchange_V(n)=1
           else
            if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),-Volume_R)
            else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),-Volume_R)
            endif
           endif
cccz rightwards flow case
         elseif(H_M.le.H_L.and.H_M.ge.H_R) then
           Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
           Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
           Mean_Height=(H_L*L_L+H_R*L_R+H_M*L_M)/(L_R+L_M+L_L)
cccz always choose the higher value
           if(H_M.ge.Mean_Height) then
             Ava_Height_R=Mean_Height
             Ava_Height_L=Mean_Height_L
           else
             Ava_Height_R=Mean_Height_R
             Ava_Height_L=Mean_Height
           endif
           Volume_R=max((Ava_Height_R-
     &       (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
           Volume_R=min(Volume_R,
     &       max(min(slopeCoord(n,2)+h_Pond_temp(n)-Ava_Height_R,
     &       h_Pond_temp(n)-CriticalH_R)*L_M,0.0D0))
                  
           Volume_L=max((Ava_Height_L-
     &       (slopeCoord(n,2)+h_Pond_temp(n)))*L_M,0.0D0)
           Volume_L=min(Volume_L,
     &       max(min(slopeCoord(n-1,2)+h_Pond_temp(n-1)-Ava_Height_L,
     &       h_Pond_temp(n-1)-CriticalH_R)*L_L,0.0D0))
                       
           if(Check_Exchange_V(n-1).eq.0) then
             Exchange_V_temp(n-1)=Volume_L
             Check_Exchange_V(n-1)=1
           else
            if(Exchange_V_temp(n-1).lt.0.0D0) then
              Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),Volume_L)
            else
              Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),Volume_L)
            endif
           endif

           if(Check_Exchange_V(n).eq.0) then
             Exchange_V_temp(n)=Volume_R
             Check_Exchange_V(n)=1
           else
            if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
            else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
            endif
           endif       
cccz flat/peak case              
         else
           Mean_Height_L=(H_L*L_L+H_M*L_M)/(L_L+L_M)
           Mean_Height_R=(H_R*L_R+H_M*L_M)/(L_R+L_M)
           Mean_Height=(H_L*L_L+H_M*L_M+H_R*L_R)/(L_L+L_M+L_R)

          if(Mean_Height_L.lt.Mean_Height) then
            Ava_Height_L=Mean_Height
            Ava_Height_R=Mean_Height_R
          else
            Ava_Height_R=Mean_Height
            Ava_Height_L=Mean_Height_L
          endif
cccz water requirement                  
          Volume_L=max((Ava_Height_L-
     &      (slopeCoord(n-1,2)+h_Pond_temp(n-1)))*L_L,0.0D0)
          Volume_R=max((Ava_Height_R-
     &      (slopeCoord(n+1,2)+h_Pond_temp(n+1)))*L_R,0.0D0)
          Volume_T=Volume_L+Volume_R
cccz water supply          
          Volume=min(Volume_T,
     &      max(min(slopeCoord(n,2)+h_Pond_temp(n)-
     &      min(Mean_Height_L,Mean_Height_R),      
     &      h_Pond_temp(n)-CriticalH_R)*L_L,0.0D0))      
cccz water redistributed     
          if(Volume_T.eq.0.0D0) then
          else
            Volume_L=Volume*Volume_L/Volume_T
            Volume_R=Volume*Volume_R/Volume_T
          endif                      
          if(Check_Exchange_V(n-1).eq.0) then
            Exchange_V_temp(n-1)=-Volume_L
            Check_Exchange_V(n-1)=1
          else
           if(Exchange_V_temp(n-1).lt.0.0D0) then
             Exchange_V_temp(n-1)=max(Exchange_V_temp(n-1),-Volume_L)
           else
             Exchange_V_temp(n-1)=min(Exchange_V_temp(n-1),-Volume_L)
           endif
          endif
          if(Check_Exchange_V(n).eq.0) then
            Exchange_V_temp(n)=Volume_R
            Check_Exchange_V(n)=1
          else
            if(Exchange_V_temp(n).lt.0.0D0) then
              Exchange_V_temp(n)=max(Exchange_V_temp(n),Volume_R)
            else
              Exchange_V_temp(n)=min(Exchange_V_temp(n),Volume_R)
            endif
          endif                
         endif   
        endif
       enddo

cccz the water right at the edge had to been drained.
cccz moreover, the edge points should not provide inwards water flux since it is already below CriticalHR
       Exchange_V_temp(1)=min(Exchange_V_temp(1),0.0D0)
       Exchange_V_temp(n)=max(Exchange_V_temp(n),0.0D0)
       error=0.0D0
              
       do n=1,SurNodeIndex
        if(n.eq.1) then
          L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n,1))
        elseif(n.eq.SurNodeIndex) then
          L_M=0.5*(slopeCoord(n,1)-slopeCoord(n-1,1))
        else
          L_M=0.5*(slopeCoord(n+1,1)-slopeCoord(n-1,1))
        endif
        if(n.eq.1) then
          h_Pond_temp(n)=h_Pond_temp(n)-Exchange_V_temp(n)/L_M
cccz if "gt", the Exchange_V_temp(1)<0, and leftwards exchange occurs
          if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            RunoffLeft_Efflux=RunoffLeft_Efflux
     &        -(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
            h_Pond_temp(n)=CriticalH_R_Sur(n)
          endif
        elseif(n.eq.SurNodeIndex) then
          h_Pond_temp(n)=h_Pond_temp(n)+Exchange_V_temp(n-1)/L_M
          if(h_Pond_temp(n).gt.CriticalH_R_Sur(n)) then
            RunoffRight_Efflux=RunoffRight_Efflux
     &        +(h_Pond_temp(n)-CriticalH_R_Sur(n))*L_M/step
            h_Pond_temp(n)=CriticalH_R_Sur(n)
          endif
        else
          h_Pond_temp(n)=h_Pond_temp(n)+
     &      (Exchange_V_temp(n-1)-Exchange_V_temp(n))/L_M
        endif
        error=max(error,abs(h_Pond_temp(n)-h_Pond(n)))
       enddo
          
       iteration_Num=iteration_Num+1
          
       if((error.gt.0.01D0*CriticalH_R).and.(iteration_Num.lt.100)) then
         do n=1,SurNodeIndex  
           h_Pond(n)=h_Pond_temp(n) 
           Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
         enddo
         h_Pond(1)=CriticalH_R_Sur(1)
         h_Pond(SurNodeIndex)=CriticalH_R_Sur(SurNodeIndex)
         goto 1206
       else
         do n=1,SurNodeIndex
          k=SurfNodeNodeIndexH(n)
          h_Pond(n)=h_Pond_temp(n)
          h_Pond(n)=max(h_Pond_temp(n),0.0D0)
          Exchange_V(n)=Exchange_V(n)+Exchange_V_temp(n)
          if(h_Pond(n).gt.CriticalH_R) then
cccz assign the boundary condition for infiltration
            hNew(k)=CriticalH+h_Pond(n)                 
          endif
         enddo
         endif
         
cccz we add a hard smoothing factor when the max ponded depth exists but very low (smaller then CriticalH_R)
cccz this is just for flat surface now and we applied this in case when mulch exists (CriticalH_R can be relatively high)
       MaxPondFlatSur=0.0D0
       AvePondFlatSur=0.0D0
       MaxFlatSurHnew=-100.0D0
       do n=1,SurNodeIndex
          kk=SurfNodeNodeIndexH(n)
          MaxFlatSurHnew=max(hnew(kk),MaxFlatSurHnew)
          MaxPondFlatSur=max(h_Pond(n),MaxPondFlatSur)
       enddo
        
       if(MaxPondFlatSur.lt.CriticalH_R-0.0001D0) then
        do n=1,SurNodeIndex
         AvePondFlatSur=AvePondFlatSur
     &     +h_Pond(n)*width(SurfNodeSurfIndexH(n))
        enddo  
        AvePondFlatSur=AvePondFlatSur/totalLength
        do n=1,SurNodeIndex
         if(n.eq.1) then
           Exchange_V(n)=Exchange_V(n)
     &       +(h_Pond(n)-AvePondFlatSur)*width(SurfNodeSurfIndexH(n))  ! rightwards is +, leftwards is -
         elseif(n.eq.SurNodeIndex) then
         else
           Exchange_V(n)=Exchange_V(n)+(h_Pond(n)-AvePondFlatSur)
     &       *width(SurfNodeSurfIndexH(n))+Exchange_V(n-1)
         endif
         h_Pond(n)=AvePondFlatSur
        enddo             
       endif
       
           
cccz now we should have a overall idea of the three types of water flow, all in unit (cm^2/day)
cccz "q_flux" is the most non-free flux which follow the SV equation
cccz "exchange_V" is the intermediately free for water level adjustment
cccz "efflux" is the very free flux
           
cccz finally, change the "exchange_V" (a volume) to a flux term
       do i=1,SurNodeIndex
         Exflux(n)=Exchange_V(n)/step
       enddo

cccz reset the boundary temperature condition
       if(MaxFlatSurHnew.gt.CriticalH) then   ! then set the constant boundary condition
         do i=1, SurNodeIndex
           n=SurfNodeNodeIndexH(i)              
           CodeT(n)=4                         ! assumption: surface runoff has air temp in another module 
           CodeG(n)=0                         ! assumption: surface ponding induce zero-flux (impermeable) boundary condition 
         enddo 
       else  
         do i=1, SurNodeIndex
           n=SurfNodeNodeIndexH(i) 
           CodeT(n)=TmprBCRecord(i)
           CodeG(n)=GasBCRecord(i)
         enddo 
       endif

1300  RunoffRight_Runoff_Loc_Record=RunoffRight
      RunoffLeft_Runoff_Loc_Record=RunoffLeft 
       
      return
      
      end