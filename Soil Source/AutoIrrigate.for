C this subroutine will carry out autmatic irrigation for the Agmip Simulations
C the trigger is if the profile has less than 50% available water down to 50 cm.
C irrigation amount is defined in the weather file and is based on rint = average infil rate.
C it is called once a day, at 5:00 am


       Subroutine AutoIrrigate()
       include 'public.ins'
       include 'puweath.ins'
       include 'PuSurface.ins'
       Parameter (PERIOD =1./24.)

       Real ThetaAvail50, ThetaFull50,AvailWaterRatio
       Real Thi, Thj, Thl, ThFl_i,ThFl_j, ThFl_l
       Real Bii(3),Cii(3), AE, Sum1, Sum2
       Integer i, j, l
       Common /Auto/ ModNum
       
       
     
     
       If(lInput.eq.1) then
C
C  Initialize 
C
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum) = time + 1
cccz initialize to zero
        do i=1,NumNPD
          Qautoirrig(i)=0.0D0
        enddo
cccz
        
       End If

       if (AutoIrrigateF.eq.0.and.lInput.lt.1)  tNext(ModNum)=1.0e10
       
11     If(abs(time-tNext(ModNum)).lt.0.001*Step.and.
     & AutoIrrigateF.gt.0) then
cccz reset to zero every time to prevent "infinite irrigation"         
        do i=1,NumNPD
          Qautoirrig(i)=0.0D0
        enddo 
cccz
  
        ThetaAvail50=0.0
        ThetaFull50 =0.0
         Do n=1,NumEl
 		   NUS=4
		   if(KX(n,3).eq.KX(n,4)) NUS=3
		   Sum1=0.
		   Sum2=0.
*         Loop on subelements
		   do k=1,NUS-2
			 i=KX(n,1)
			 j=KX(n,k+1)
			 l=KX(n,k+2)
			 Cii(1)=x(l)-x(j)
			 Cii(2)=x(i)-x(l)
			 Cii(3)=x(j)-x(i)
			 Bii(1)=y(j)-y(l)
			 Bii(2)=y(l)-y(i)
			 Bii(3)=y(i)-y(j)
			 AE=(Cii(3)*Bii(2)-Cii(2)*Bii(3))/2.
			 Thi=0
			 Thj=0
			 Thl=0
			 ThFl_i=0
			 ThFl_j=0
			 ThFl_l=0
			 
			 If (y(1)-y(i).le.50) Then
			  ThFl_i=ThFull(MatNumN(i))
			  Thi=ThAvail(i)
			 End If
			 If (y(1)-y(j).le.50) Then
			 ThFl_j=ThFull(MatNumN(j))
			 Thj=ThAvail(j)
			 End If
			 If (y(1)-y(l).le.50) Then
			 ThFl_l=ThFull(MatNumN(l))
			 Thl=ThAvail(l)
			 End If
			 
			 ThetaFull50=ThetaFull50+AE*(ThFl_i+ThFl_j+ThFl_l)/3.
			 ThetaAvail50=ThetaAvail50+AE*(Thi+Thj+Thl)/3.
		   Enddo
          End Do
          if (thetaFull50.gt.0) AvailWaterRatio=ThetaAvail50/ThetaFull50
           
C      Auto irrigate if available water is less than 60%
C  irrigation goes into Q values
C TODO need to check code this was partially merged between the gas exchange and mulch versions
          If (AvailWaterRatio.lt.(0.60).and.lInput.lt.1) then
           do k=1, NumBp
            i=KXB(k)
            If(abs(CodeW(i)).eq.4) then
             if (Q(i).lt.0) then 
              CodeW(i)=4
cccz since we treat auto-irrigation as rainfall, we directly change the varbw terms
cccz so varbw_air will be rainfall+autoirrigation
              Varbw_Air(k,1)=Varbw_Air(k,1)+AutoIrrigAmt
              Varbw_Mulch(k,1)=Varbw_Mulch(k,1)+AutoIrrigAmt
              VarBW(k,1)=Varbw(k,1)+AutoIrrigAmt
cccz then re-calculation of Q, folloiwng VarBW and the starndard method show in weather
              VarBW_Air(k,3)=VarBW_Air(k,2)-VarBW_Air(k,1)
              VarBW_Mulch(k,3)=VarBW_Mulch(k,2)-VarBW_Mulch(k,1)
              VarBW(k,3)=VarBW(k,2)-VarBW(k,1)
              Q(i)=-Width(k)*VarBW(k,3)
cccz now the autoirrigation is recorded by the "varbw flux temrs"
cccz   Q(i)=Q(i)+AutoIrrigAmt*Width(k)   ! add to irrigation so we don't lose any rain if needed
cccz                   if (Q(i).gt.0.0) CodeW(i)=-4   ! make sure bc changes if Qn goes > 0 (infiltration)
             End if    !Q(i) <0
            End if   ! CodeW=4
          End Do
              
         End If  ! end auto irrigation
           tNext(ModNum) = time +PERIOD
           
         End if  ! Tnext hourly calcs     
           
           Return
           end