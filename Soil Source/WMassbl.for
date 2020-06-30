cdt 7/9/07 added code to calculate loss of N in flux out of the domain 
Cdt N in 2DSOIL is nitrate (NO3) and the units are ug per cm3 or milligram/liter
      Subroutine Water_Mass_Balance()
      include 'Public.ins'
      include 'puplant.ins'
      include 'PuSurface.ins'
      
      Dimension Bi(3),Ci(3)
      Character*10 Date
      Real AE,
     &Mean1,Mean2,Mean3,Mean4,Mean5
      Integer ModNum
      common /SurW_BAL/ModNum,CFlux,CFluxPrevious,W_Sum

      If (lInput.eq.1) then
          open(91,file=MassBalanceFile,status='unknown',recl=150)
          write(91,'(3A12)') 'time','Date','water'

c CCCCCCCCCCCCCCCCCCCCCC Zhuangji's Runoff File CCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(99,file=MassBalanceRunoffFileOut,status='unknown',recl=1200)
      write(99,'(7A20)') 'time','Date','Rainfall','Inifiltration',
     &         'Runoff@end', 'WaterStorage', 'WaterDepth'
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          NumMod=NumMod+1
          ModNum=NumMod
          tNext(ModNum)=time
          
          W_Sum=0.0
          W_Sum_Old=0.0

c CCCCCCCCCCCCCCCCCCCCCC Zhuangji's Water Balance Term CCCCCCCCCCCCCCCCCCCC
        Surface_Flux=0.0
        Bottom_Flux=0.0
        Runoff_Flux=0.0
        Runoff_Terminal_Flux=0.0
        Soil_Water_Diff=0.0
        Bal_Z=0.0
        Bal_Rain_Runoff=0.0
        Rain_Flux=0.0
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Endif

      t=time
      if (Abs(time-tNext(ModNum)).lt.0.001*Step.or.lInput.ne.0) then
cccz comment this line out for Ex_1 but keep for all others 
cccz #ifdef EX_1       
cccz      X_Period=1.0/24.0D0/60.0D0/4.0D0
cccz #else      
      X_Period=1.0/24.0D0
cccz #endif
      tNext(ModNum)=tNext(ModNum)+X_period

c CCCCCCCCCCCCCCCCCCCCCCCCCCCC Zhuangji's Calcualate the soil water quantity CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      If (lInput.ne.1) then 
c           tNext(ModNum)=time+Step
           W_Sum_Old=W_Sum
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          
c process the water balance every 10 minutes
          W_Sum=0.0
          Sum=0.
          W_Sum=0.
          Do n=1,NumEl
              NUS=4
              if(KX(n,3).eq.KX(n,4)) NUS=3
              Sum1=0.
              

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
                  Sum1=Sum1+AE*(ThNew(i)+ThNew(j)+ThNew(l))/3.
              Enddo

              W_Sum=W_Sum+Sum1
          Enddo
          
c CCCCCCCCCCCCCCCCCCCCCCCCCCCC Zhuangji's Calcualate the soil water quantity CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Soil_Water_Diff=W_Sum-W_Sum_Old
      Bal_Z= Surface_Flux+Bottom_Flux-Soil_Water_Diff
      Bal_Rain_Runoff=Rain_Flux-Surface_Flux-Runoff_Flux
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C calculate sum of N in all forms in soil
C factor is appropriate for mg/cm3 Total N is mg per slab (grid width x 1cm)
C Mineral N is ug/cm3. Do a summation over the simulation domain - total ug in a plant slab
C Fact now works while there is a plant - need to make this more robust where we use the dimensions
C of the domain
        
          iday=int(t)
          call caldat(iday,mm,id,iyyy)
          write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy
          write(91,'(1F12.4,A12,F12.4,E14.4,F14.4,F14.4)') time,date,   
     &            W_Sum, Bal_Z, Bal_Rain_Runoff, Runoff_Flux
     
c CCCCCCCCCCCCCCCCCCCCCC Zhuangji's Runoff File CCCCCCCCCCCCCCCCCCCCCCCCCCC
           write(99,'(1F12.4,A12,F14.4,F14.4,F14.4,F14.4,
     &            F14.4,F14.4,F14.4,F14.4,F14.4,F14.4,F14.4, 
     &            F14.4,F14.4,F14.4,F14.4)') 
     &            time, date,
     &            Rain_Flux, Surface_Flux,
     &            Runoff_Terminal_Flux, WaterStorage, h_Stay(1),
     &            h_Stay(2), h_Stay(3), h_Stay(4), 
     &            h_Stay(5), h_Stay(6), h_Stay(7), 
     &            h_Stay(8), h_Stay(9), 
     &            Runoff_Left_Sum,Runoff_Right_Sum
     
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c CCCCCCCCCCCCCCCCCCCCCC Zhuangji's Water Balance Term CCCCCCCCCCCCCCCCCCCC
        Surface_Flux=0.0
        Bottom_Flux=0.0
        Runoff_Flux=0.0
        Runoff_Terminal_Flux=0.0
        Soil_Water_Diff=0.0
        Bal_Rain_Runoff=0.0
        Rain_Flux=0.0
        Bal_Z=0.0
        h_stay_water=0.0
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      endif

      return
      end


