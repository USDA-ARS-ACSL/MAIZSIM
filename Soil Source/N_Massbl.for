cdt 7/9/07 added code to calculate loss of N in flux out of the domain 
Cdt N in 2DSOIL is nitrate (NO3) and the units are ug per cm3 or milligram/liter
       Subroutine Nitrogen_Mass_Balance()
        include 'Public.ins'
        include 'nitvar.ins'
        include 'puplant.ins'
        Dimension Bi(3),Ci(3)
        Character*10 Date
        Real AE,Profile_N,Manure_N,Litter_N,
     &    Min_N,Ammon_N,Org_N,Denitr,
     &    Mean1,Mean2,Mean3,Mean4,Mean5,All_N
        Integer ModNum
        common /N_BAL/ModNum,CFlux,CFluxPrevious

      If (lInput.eq.1) then
        open(91,file=MassBalanceFileOut,status='unknown',recl=150)
        write(91,'(11A12)') 'time','Date','Min_N','Org_N','Manure_N',
     !    'Litter_N'
     !           ,'Ammon_N','All_N','water','Denitr', 'CFlux'
 
        NumMod=NumMod+1
	ModNum=NumMod
	tNext(ModNum)=time
	CFlux=0
	CFluxPrevious=0
 
       Endif
          do i=1,NumBp
	     n=KXB(i)
C Cflux is loss of N in mg
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))) then
	        Cflux=Cflux+Q(n)*conc(n,1)*step*14/62
           endif
C for case of downward drainage and a constrant BC 
C (does not account for upward flow yet)
           if ((CodeW(n).eq.(1))) then
             if (Q(n).lt.0) Cflux=Cflux+Q(n)*conc(n,1)*step*14/62 ! only consider chemical out- flux>0
           endif
	  EndDo
	  
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


	   Sum=0.
           W_Sum=0.
	   Do n=1,NumEl
             NUS=4
             if(KX(n,3).eq.KX(n,4)) NUS=3
             Sum1=0.
             Sum2=0.
             Mean1=0.0
             Mean2=0.0
             Mean3=0.0
             Mean4=0.0
             Mean5=0.0
  
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
               Mean1=Mean1+AE*(Conc(i,1)*ThNew(i)+Conc(j,1)*ThNew(j)
     &                 +Conc(l,1)*ThNew(l))/3.*14./62.
               Sum1=Sum1+AE*(ThNew(i)+ThNew(j)+ThNew(l))/3.
               Mean2=Mean2+AE*(Nh(i)+Nh(j)+Nh(l))/3.  ! organic N is mg cm3 area
               Mean3=Mean3+AE*(NL(i)+NL(j)+NL(l))/3.
               Mean4=Mean4+AE*(Nm(i)+Nm(j)+Nm(l))/3.
               Mean5=Mean5+AE*(NNH4(i)+NNH4(j)+NNH4(l))/3.
               Sum2 =Sum2 +AE*(Denit(i)+denit(j)+denit(l))/3.
             Enddo
             Min_N=Min_N+Mean1
             Org_N=Org_N+mean2
             Litter_N=Litter_N+Mean3
             Manure_N=Manure_N+Mean4
             Ammon_N=Ammon_N+Mean5
             W_Sum=W_Sum+Sum1
             Denitr=Denitr+sum2
         Enddo

C calculate sum of N in all forms in soil
C factor is appropriate for mg/cm3 Total N is mg per slab (grid width x 1cm)
C Mineral N is ug/cm3. Do a summation over the simulation domain - total ug in a plant slab
C Fact now works while there is a plant - need to make this more robust where we use the dimensions
C of the domain


          fact=1/(0.01*RowSp/100.*EOMult) !m2 of slab
          fact=fact*10000/1000/1000/1000    !m2 of slab ->ha, ug-->Kg
          All_N=(Min_N*fact+Org_N*fact+Litter_N*fact+Manure_N*fact
     !           +Ammon_N*fact)
         iday=int(t)
	   call caldat(iday,mm,id,iyyy) 
         write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          write(91,10) 
     !         time,date,Min_N*fact,
     !         Org_N*fact,Manure_N*fact,
     !         Litter_N*fact,Ammon_N*fact,All_N, W_Sum,Denitr*fact,
     !         CFLux*fact

          CFluxPrevious=CFlux
        endif
10    Format (1F12.4,',', A12,',', 8(F12.4, ','),F12.4)
      return
      end


