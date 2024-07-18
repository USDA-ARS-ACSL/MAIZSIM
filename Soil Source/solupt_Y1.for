c iSink values 
c         0  'passive root mass uptake' no kinetics
c         1   Plant regulation diffusive and convective uptake
c         2  'convective diffusive uptake no plant regulation (Constant IMax)
c         3  'diffusive uptake only, no effect of water movement 
c         4  Michaelis-Menton uptake only
c TODO take out CSink_M as the nodal value works since this is no longer element wise

      subroutine SoluteUptake()
      include 'public.ins'
      IncLude 'puplant.ins'
      Include 'puweath.ins'
      real MMUpN, bi(3),ci(3),PotNitrogen_t
      Character InString*132 ! to read input file
      Common  / SUP /  TotWSink,TotSSink,WincrSink,TotSSINK2
      Common  / Biom /  BM
      Common  / slpt /  Iflag
      Common / Mass_root_M / Isink,RootRadius,
     !     DlngR(NMatD), Disp(NumNPD,2),RootExchange(2),
     !	CSink_M(NumNPD),Cr_M(NumNPD,2), LastCr_M(NumNPD,2),
     !  TotalFUPM,TotalFUPY
   
      Common  / SoluteM /  ThOld(NumNPD), Dmol(NumSD),
     !                     DLng(NMatD,NumSD),Dtrn(NMatD,NumSD)
      Common  / PotNitr /  PotNitrogen_t
      Character * 40  UptakeN(6)
c      data UptakeN / 
c     !'Passive Mass Uptake (Mihaelson-Menton kinetics',
c     !'Passive Mass Uptake (Mihaelson-Menton kinetics)&convetcive flux',
c     !'Active  Mass Uptake  Cr = Croot',
c     !' ',
c     !' ',
c     !' ' / 
c       real qc(NumElD),DispM(2,NumElD)
C *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
	If(lInput.eq.1) then
c        Open(94,file = 'TotSink.res')
		TotWSink = 0.
		TotSSINK = 0.
		TotSSINK2 = 0.0
		Iflag = 0
		MMUpN = 0.0
		
		
   
C 
C  Reading of the input files and initial calculations 
C
		im = 430
		il = 0
		Open(40,file = VarietyFile, status = 'old',ERR = 10)
! look for section with soil nitrogen information
 15       Read (40,'(A132)') InString
          if (InString(1:14).ne.'[SoilNitrogen]') goto 15
                  
		
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10) Isink,RootRadius
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Do M = 1,2
		 Read(40, * ,ERR = 10) ConstI(m),constK(m), Cmin0(M)
		 il = il + 1
		enddo
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
	
		close(40)
          
cccz zhuangji updates the initial set for Cmin0, node based 
		do j=1, 2
		 do n=1, NumNP
		   LastCr_M(n,j)=Cmin0(j)
		  enddo
		Enddo 
          
c	  Open(101, file = 'tempC.txt', recl=300)
c		write(101, * )      UptakeN(Isink)
c		write(101,103)'          e   Csink        CsinkN       Croot      
c     $     Croot     CMean   alphaK    alphaK   alphaKK    AlphaKK
c     $      F_MM      SumSink' 
c		write(101,103)'             passive     Total    old     young    conc
c     $     young     old    young     old      sum       sumSink' 
103		format(a) 
      Endif
      
      If (NShoot.eq.0) then
          csink(:,1) =0.0
          return
       endif
      
      if(ITIME.EQ.1) iflag =  1
      SIncrSink = 0.0
      WIncrsink = 0.0
      SincrSink2 = 0.0
      PotNitrogen_t=NitroDemand*62.0/14.0  ! go from nitrogen to nitrate ug no3 day-1 slab-1
      TotalFUPY=0                            !Total root length
      TotalFUPM=0
      	Do M = 1,NMat
		 DlngR(M)=Dlng(M,1)
		enddo
       
       Do n=1,NumNP
          Csink_M(n)=Conc(n,1)   ! mean of concentration for element 
          CSink(n,1)=Conc(n,1)*Sink(n) 
          M=MatNumN(n)
          TotalFUPY=TotalFUPY+FUP(n,2)
          TotalFUPM=TotalFUPM+FUP(n,1)
          SIncrSink2=SIncrSink2+cSink(n,1)*step*14./62.*nodeArea(n)
          Disp(n,1)=Dmol(1)*ThNew(n)*Tau(hNew(n))+DlngR(M)*VUP(n,1)
		Disp(n,2)=Dmol(1)*ThNew(n)*Tau(hNew(n))+DlngR(M)*VUP(n,2)	
       Enddo
       
       
	 if (iSink.eq.1) then        ! plant regulation
         dC=0.1*NDemandError+0.05*cumulativeNDemandError 
	   ConstI(1)=amin1(75.0,ConstI(1)-dC)
	   ConstI(1)=amax1(ConstI(1),3.0)
	   ConstI(2)=ConstI(1)/2.0
	 ! change to 0 and wait for crop to be executed before
	 ! recalculating
	 !  NDemandError=0.0
	 !  cumulativeNDemandError=0.0
	   EndIf
C calculate nutrigen root uptake    
	 if (iSink.gt.0) Call massRootflux
c -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -     

cccz zhuangji for node based representation       
       do n = 1, NumNP  
         SIncrSink=SIncrSink+CSink(n,1)*Step*nodeArea(n)*14./62. ! go from Nitrate to Nitrogen 
         WIncrSink=WIncrSink+Sink(n)*nodeArea(n)*Step
       enddo

       TotWSink = TotWSink + WIncrSink   ! water sink
       TotSSINK = TotSSINK + SIncrSink   ! solute sink based on csink after dispersive calcs
       TotSSINK2 = TOTSSink2 + SIncrSink2 ! solute sink based on the active uptake before calling sub
c       If(ITIME.eq.24.and.iflag.eq.1) then
c         Write(94, * ) Jday,TotWSink,TotSSink,TotSSInk2
c         iflag = 0
c       Endif
       !if (iSink.EQ.0) SIncrSink=SIncrSink2   ! set passive uptake
      return
10    Write(*,*) 'Error in solute file'
      Call errmes(im,il)
      end
     
      subroutine massRootflux()
      include 'public.ins'
      IncLude 'puplant.ins'
      Include 'puweath.ins'
	real * 8 qsinkC(NumNPD,2),alphaK(NumNPD,2)
	Real LastCr_M,F_MM(NumNPD,2)
      integer e,Iroot(NumNPD,2)
	Real *8 determ,b,a,c,bet,beta2,beta2a,alf,alf2,part1,part2,cr,betR,
     !   gammaB, gammaB2,Pc1, Pc2, pi, TwoPi,Ro2,AlphaKK(2)
      Common / Mass_root_M / Isink,RootRadius,
     !   DlngR(NMatD), Disp(NumNPD,2),RootExchange(2),
     !	CSink_M(NumNPD),Cr_M(NumNPD,2), LastCr_M(NumNPD,2), 
     !    TotalFUPM,TotalFUPY
	Common  / PotNitr /  PotNitrogen_t
	  
	 pi = 3.14156
	 TwoPi = 6.28312
	 pi4 = 12.56
	 Ro2 = RootRadius * RootRadius                 ! root radius (Ro)
	 alpha2=0.05
C FUP is root density, VUP is water uptake per unit root length
c calculate and alternative alphaK depending on demand
       

       do j=1,2
           RootLenFact=0.0
           do n=1,NumNP
           RootLenFact=RootLenFact+FUP(n,j)*Cr_M(n,j)
           enddo
           if (RootLenFact.LE.0.0001) then  ! protect from zero values
            alphaKK(j)=0.0
            else
            alphaKK(j)=PotNitrogen_t/(TwoPi*RootRadius*RootLenFact)
          endif
       enddo         
	 do j = 1,2
	  do n = 1, NumNP 
	   Cr_M(n,j)=LastCr_M(n,j)
	   do k=1,3
         Iroot(n,j) = 0                 ! 1 if length >0
	   qsinkC(n,j) = 0.0 
	   if(FUP(n,j).gt.0.0001) THEN  ! don't want use the code if FUP is very small DT 8/31/2021
	    alphak(n,j)=(ConstI(j))/(ConstK(j)+Cr_M(n,j)) 
	    F_MM(n,j)=alphak(n,j)*Cr_M(n,j)
	    Rx = 1.0 / sqrt(pi * (FUP(n,1) + FUP(n,2)))
	    RootRatio=Rx/RootRadius
	    
	    SC = 0.017  -  (PSIS(n) * 0.5)  
	    betR = SC / RootRadius
	    RootRatio = dmin1(betR,RootRatio)
          if (RootRatio.le.1.0) then 
             RootRatio=1.1
           endif
	    RootRatio2=RootRatio*RootRatio
	     if (Disp(n,j).le.0.0) then
             GammaB=0
             PC=1.0 ! PC is not used but this is just to be safe
	      else
	      If ((VUP(n,j).le.1.0e-6).OR.(iSink.EQ.3)) Then
	       GammaB=alphaK(n,j)*RootRadius/(Disp(n,j))
	       PC=1.0-0.5*GammaB
	       PC2=Rx*Rx*GammaB*alog(RootRatio)/(Rx*Rx-RootRadius**2)
	       PC=1.0/(PC+PC2)
	       PC=amax1(0.00001,amin1(1.0,PC))  ! prevent PC from going to zero
	         Else
	       GammaB=RootRadius*VUP(n,j)/Disp(n,j)
	       GammaB2=2.0/(2.0-GammaB)
	       PC=(RootRatio**(2.0-GammaB)-1.0)/(RootRatio**2-1.0)
	       PC=AlphaK(n,j)+(VUP(n,j)-AlphaK(n,j))*GammaB2*PC
	       if (PC.le.0.0) then
	         PC=0.00001
	         Else
                PC=amax1(0.00001,amin1(1.0,VUP(n,j)/PC))	       
               EndIf
               if (isNAN(PC)) then
                iii=1
               endif
	       EndIF  !VUP <=0
	     EndIf   !Disp <= 0
	    qSinkC(n,j)=TwoPi*RootRadius
	    qSinkC(n,j)=qSinkC(n,j)*alphaK(n,j)*PC*Csink_M(n)	
	    Cr_M(n,j)=PC*CSink_M(n)         ! update value at root, may have to iterate      
	    Cr_M1=Cr_M(n,j)   
	    if (CSink_M(n).LT.amax1(0.01,Cr_M(n,j))) then
	       qSinkC(n,j)=0.0
	      endif
         End If   !FUP <0
      End Do   ! k
      End Do   ! j
      End Do   !e
      
	 do n = 1, NumNP
	  sumSinkN = 0.0
	  do j = 1,2
	   if (iSink.EQ.4) then
	   dif=amax1(0.0,CSink_m(n)-CMin0(j))
	     sumSinkN=sumSinkN + ConstI(j)*FUP(n,j)*dif
     !        /(ConstK(j)+dif)
          else
	   sumSinkN = sumSinkN + qsinkC(n,j)*FUP(n,j)
	   LastCr_M(n,j) = Cr_M(n,j)    
	   EndIF
	  enddo	  
	  If(CSink(n,1).gt.0)then
c	    write(101,102)time, e, CSink(e,1),qSinkC(e,1)+qSinkC(e,2),
c     !	 (Cr_M(e,j), j = 1,2) , Csink_m(e), (ConstI(j),j=1,2),
c     !     (AlphaK(e,j),j=1,2),F_MM(e,1)+F_MM(e,2),SumSinkN
        continue
	  endif
	    CSink(n,1) = amax1(sumSinkN,0.0)
	 enddo
c102    format(g14.8,1x,I5,12(2x,e10.3))
	return
	end
      
      Subroutine Matric(N,A,B,C,F)
	 real A( * ),B( * ),F( * ),C( * ),ALF(100)
	 N1 = N - 1
       R = 0.0
	 DO I = 1,N1
	  ALF(I) = B(I) / A(I)
	  R = R + ALF(I) * C(I)
	 ENDDO
	 CR = F(N) / R
	 DO I = 1,N1
        A(I) = F(I) / (ALF(I) + A(I))
	 ENDDO
	 RETURN
	END 