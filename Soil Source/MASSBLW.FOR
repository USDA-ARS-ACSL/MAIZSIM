      Subroutine MassBl()
*     File writes are to Unit 36- Conc_2d.out- opened in the
*     Output.For Subroutine
*     but if Movers(2)=0 file writes are to unit 31 -theta file
      Include 'public.ins'
      Integer out_unit,PLevel
      Real Mass_out,Init_mass
      character*20 Date(50), temp
      Common /MassBlCom/ Conc_old(50,12),Tprint(50),Mass_out(12),
     !                   PLevel,Total_Seepage(12),NumPrint,
     !                   Water_Out,Cum_Water_Out,Water_in,
     !                   Cum_water_in,BStorage,W_Sum,Init_mass,
     !                   out_unit,NumSol2,Cum_EO,Cum_Tpot,
     !                   P_Trans_out,P_evap_out,Evap_out,
     !                   Transpiration_sum,ModNum
      Dimension Bchk(NumElD),Bi(3),Ci(3),W_Bchk(NumElD)

      If(lInput.eq.1) then
        out_unit=86
        Open(40,file=MassBalanceFile,err=20)
        Read(40,'(A20)',Err=20) temp
        Read(40,'(A20)',Err=20) temp
        Read(40,*,Err=20) NumPrint
        Read(40,*,Err=20) 
        Read(40,*,Err=20) (Date(i),i=1,NumPrint)
        Close(40)
	  open(out_unit,file="MassBal.out")
	  Do i=1, Numprint
	     Tprint(i)=Julday(Date(i))
	  EndDo
	
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum)=time
        if (NumPrint.eq.0) tNext(ModNum)=1e32
        Cum_water_out=0.0
        Cum_water_in=0.0
        Water_in=0.0
        Water_Out=0.0
        Evap_out=0.0
        P_Evap_out=0.0
        P_Trans_out=0.0
        Cum_EO=0.0
        Cum_Tpot=0.0
        Transpiration_sum=0.0
        
        Plevel=1
       EndIf
       
	

* this ends the steps carried out on first call 
*to the subroutine 
	
* these steps are caried out every time step
	do i=1,NumBP
		 n=kxb(i)
		 k=CodeW(n)
		 if(k.eq.2.or.k.eq.-5) then
			Water_Out=Water_Out+Q(n)*step
		  endif
		 if(k.eq.1) then
			 Water_in=Water_in+Q(n)*step
		 endif
		 if(abs(k).eq.4) then
			 P_Evap_out=P_Evap_out+varBW(i,2)*step*Width(i)
			 if(Q(n).ge.0) then
				Water_in=Water_in+Q(n)*step
			 endif
			 if(Q(n).lt.0) then
			   Evap_out=Evap_out+Q(n)*step
			 endif  
		 endif
	enddo  

* Calculate total water uptake by crop for area under
* consideration if there is root water uptake
	Sum3=0.0
	if(NShoot.ne.0) then
		 do n=1,NumEL      
			Sum3=Sum3-Area(n)*Sink(n)*step
		  enddo
	endif
		Transpiration_sum=Transpiration_sum+sum3
		Cum_EO=Cum_EO+EO*step*24.0
		Cum_Tpot=Cum_TPot+Tpot*step*24.0
		P_Trans_out=P_Trans_out-eor*Step*24.0

* The following are carried out only when it is time to write

	If(Abs(Time-tNext(ModNum)).lt.0.001*Step.or.lInput.eq.1) then
        	BStorage=W_Sum
     		   Plevel=Plevel+1
		   if (PLevel.gt.NumPrint) then
		     tNext(ModNum)=1.E+32
		     else
		     tnext(ModNum)=Tprint(plevel)
		   endif
		
  
CYAP
	 Sum=0.
		 W_Sum=0.
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
			 Ci(1)=x(l)-x(j)
			 Ci(2)=x(i)-x(l)
			 Ci(3)=x(j)-x(i)
			 Bi(1)=y(j)-y(l)
			 Bi(2)=y(l)-y(i)
			 Bi(3)=y(i)-y(j)
			 AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
			 
			 Sum2=Sum2+AE*(ThNew(i)+ThNew(j)+ThNew(l))/3.
		   Enddo
		   W_Sum=W_Sum+Sum2
		   W_Bchk(n)=W_Sum
		   
		Enddo
		write(out_unit,*)
		write(out_unit,*)'       *************         '
		write(out_unit,*) ' Time=', date(pLevel),time
		Total_water_out=Water_Out+evap_out+
     !                Transpiration_sum
				 
c       Write(out_unit,*) (W_Bchk(n),n=1,NumEl)
		Write(out_unit,*) 'previous and current water in soil', 
     !                     BStorage,W_sum
		Write(out_unit,*)
		Write(out_unit,*) 'Water out during time step=',
     !                    Total_Water_Out
*        Write(out_unit,*) 'Cum water out=',Cum_Water_out
		Write(out_unit,*) 'Water in during time step=',Water_in
*        Write(out_unit,*) 'Cum water in=',Cum_water_in         
		Write(out_unit,*) 'Change in water=',W_sum-BStorage
		Write(out_unit,*) 'Compare next line to  change in water'
		Write(out_unit,*) 'Total mass balance for Water=',
     !                     Total_water_out+water_in
		Write(out_unit,*)
c        Write(out_unit,*) 'Water error=',-Cum_water_out-
c     !           (W_Sum-BStorage)+Cum_water_in
		write(out_unit,*) 'Drainage=',water_out
		write(out_unit,*) 'Transpiration',Transpiration_sum
		Write(out_unit,*) 'Potential Soil Evap for period=',P_Evap_out
		Write(out_unit,*) 'Actual Soil Evap for period=',Evap_out
		Write(out_unit,*) 'Poten Transpiration for period=',P_Trans_out
		Write(out_unit,*) 'Actual Trans=', Transpiration_sum

		write(out_unit,*)
		Water_in=0.0
		Water_Out=0.0
		P_Trans_out=0.0
		Cum_eo=0.0
     
		Evap_out=0.0
		P_Evap_out=0.0
		Transpiration_sum=0.0
	
		
	endif
	 Return
20     Stop 'MassBl data error'
	End
		
