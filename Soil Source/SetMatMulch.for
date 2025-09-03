cccz *********************************************************************
c     Record the mulch hydrological information
c     based on van Genucheton equation and Campbell format
c     Used in mulch surface water/energy module
cccz  10/29/2019
cccz *********************************************************************
      
      double precision function WQ_CERES_MULCH
     &        (hInput,rho_bulk,F_lign,lInput)
       double precision rho_bulk,F_lign,hInput
       integer lInput
       double precision F_lign_loc,a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,HH_m,rho_w_loc,WQM
       common /WQ_CERES_PRIOR/ a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,rho_w_loc
       If(lInput.eq.1) then
        Mpa_2_Watercm=10193.68D0
        rho_w_loc=1.0D6
        a_1=-20.086D0
        a_2=-0.249D0
        b_1=0.324D0
        b_2=0.124D0
        S_1=7.0993D0
        S_2=-0.079D0
        HMin_m=-200.0D0*Mpa_2_Watercm
        Hs_m=-10000.0D0
       Endif
       HH_m=max(dble(hInput),HMin_m)
       F_lign_loc=F_lign*100.0D0
       if(HH_m.lt.Hs_m) then
        HHm=HH_m/Mpa_2_Watercm
        WQM=((a_1*exp(a_2*F_lign_loc))/HHm)
     &       **(1.0D0/(b_1+b_2*F_lign_loc))
        WQM=WQM*rho_bulk/rho_w_loc
        WQ_CERES_Mulch=min(WQM,0.9D0)
       else
        WQM=S_1*exp(S_2*F_lign_loc)*rho_bulk/rho_w_loc
        WQ_CERES_MULCH=min(WQM,0.9D0)
       endif 
       return
      END
      
      double precision function WC_CERES_MULCH
     &        (hInput,rho_bulk,F_lign,lInput)
       double precision rho_bulk,F_lign,hInput
       integer lInput
       double precision F_lign_loc,a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,HH_m,rho_w_loc, WCM
       common /WC_CERES_PRIOR/ a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,rho_w_loc
       If(lInput.eq.1) then
        Mpa_2_Watercm=10193.68D0
        rho_w_loc=1.0D6
        a_1=-20.086D0
        a_2=-0.249D0
        b_1=0.324D0
        b_2=0.124D0
        S_1=7.0993D0
        S_2=-0.079D0
        HMin_m=-200.0D0*Mpa_2_Watercm
        Hs_m=-10000.0D0
       Endif
       HH_m=max(dble(hInput),HMin_m)
       F_lign_loc=F_lign*100.0D0
       if(HH_m.lt.Hs_m) then
        HHm=HH_m/Mpa_2_Watercm
        WCM=abs((1.0D0/(b_1+b_2*F_lign_loc))
     &     *(((a_1*exp(a_2*F_lign_loc))/HHm)
     &     **(1.0D0+1.0D0/(b_1+b_2*F_lign_loc)))
     &     /(a_1*exp(a_2*F_lign_loc))
     &     *rho_bulk/rho_w_loc)
        WC_CERES_MULCH=max(WCM,1.0D-37)
       else
        HH_m=Hs_m/Mpa_2_Watercm
        WCM=abs((1.0D0/(b_1+b_2*F_lign_loc))
     &     *(((a_1*exp(a_2*F_lign_loc))/HH_m)
     &     **(1.0D0+1.0D0/(b_1+b_2*F_lign_loc)))
     &     /(a_1*exp(a_2*F_lign_loc))
     &     *rho_bulk/rho_w_loc)
        WC_CERES_MULCH=max(WCM,1.0D-37)
       endif 
       return
      END
      
      double precision function WH_CERES_MULCH
     &      (theInput,rho_bulk,F_lign,lInput)
       double precision rho_bulk,F_lign,theInput
       integer lInput
       double precision F_lign_loc,a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,HH_m,rho_w_loc,WHM
       common /WH_CERES_PRIOR/ a_1,a_2,b_1,b_2,S_1,S_2,
     &     Mpa_2_Watercm,HMin_m,Hs_m,rho_w_loc
       If(lInput.eq.1) then
        Mpa_2_Watercm=10193.68D0
        rho_w_loc=1.0D6
        a_1=-20.086D0
        a_2=-0.249D0
        b_1=0.324D0
        b_2=0.124D0
        S_1=7.0993D0
        S_2=-0.079D0
        HMin_m=-200.0D0*Mpa_2_Watercm
        Hs_m=-10000.0D0
       Endif
       F_lign_loc=F_lign*100.0D0
       WHM=(a_1*exp(a_2*F_lign_loc))
     &    *((theInput*rho_w_loc/rho_bulk)**(-(b_1+b_2*F_lign_loc)))
     &    *Mpa_2_Watercm
       WHM=min(WHM,Hs_m)
       WH_CERES_MULCH=max(WHM,HMin_m)
       return
      END