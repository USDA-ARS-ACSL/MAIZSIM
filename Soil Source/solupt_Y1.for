c iSink values 
c         0  'passive root mass uptake' no kinetics
c         1   Plant regulation diffusive and convective uptake
c         2  'convective diffusive uptake no plant regulation (Constant IMax)
c         3  'diffusive uptake only, no effect of water movement 
c         4  Michaelis-Menton uptake only

      subroutine SoluteUptake()
      include 'public.ins'
      IncLude 'Puplant.ins'
      Include 'Puweath.ins'
      real MMUpN, bi(3),ci(3),PotNitrogen_t
      Character InString*132 ! to read input file
      Common  / SUP /  TotWSink,TotSSink,WincrSink,TotSSINK2
      Common  / Biom /  BM
      Common  / slpt /  Iflag
      Common / Mass_root_M / Isink,RootRadius,
     !     DlngR(NMatD), Disp(NumElD,2),RootExchange(2),
     !	CSink_M(NumElD),Cr_M(NumElD,2), LastCr_M(NumElD,2),
     !  TotalFUPM,TotalFUPY
   
      Common  / SoluteM /  ThOld(NumNPD), Dmol(NumSD),
     !                     DLng(NMatD,NumSD),Dtrn(NMatD,NumSD)
      Common  / PotNitr /  PotNitrogen_t
      Character * 40  UptakeN(6)
c      data UptakeN / 
c     !'Passive Mass Uptake (Mihaelson-Menton kinetics'
c     !,'Passive Mass Uptake (Mihaelson-Menton kinetics)&convetcive flux'
c     !,'Active  Mass Uptake  Cr = Croot',
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
		do j=1, 2
		 do e=1, NumEl
		   LastCr_M(e,j)=Cmin0(j)
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
      
      If (NShoot.eq.0) return
      
      if(ITIME.EQ.1) iflag =  1
      SIncrSink =  0.
      WIncrsink = 0.
      SincrSink2 = 0.0
      PotNitrogen_t=NitroDemand*62.0/14.0  ! go from nitrogen to nitrate ug no3 day-1 slab-1
      TotalFUPY=0                            !Total root length
      TotalFUPM=0
      	Do M = 1,NMat
		 DlngR(M)=Dlng(M,1)
		enddo


     
c new mass uptake	
C
C calculate sum of N in all forms in soil using triangluar elements
C
C units of conc ug/liter 
	 Do n = 1,NumEl
           cSink(n,1) = 0.
           NUS = 4
           DmE=0
           if(KX(n,3).eq.KX(n,4)) NUS = 3
*          Loop on subelements
           AE = 0.0
           M = MatNumE(n)
           TotalFUPY=TotalFUPY+FUP(n,2)
           TotalFUPM=TotalFUPM+FUP(n,1)
		 do k = 1,NUS - 2
             i = KX(n,1)
             j = KX(n,k  +  1)
             l = KX(n,k  +  2)        
		     Ci(1) = x(l) - x(j)
             Ci(2) = x(i) - x(l)
             Ci(3) = x(j) - x(i)
             Bi(1) = y(j) - y(l)
             Bi(2) = y(l) - y(i)
             Bi(3) = y(i) - y(j)
             AE = AE  +  (Ci(3) * Bi(2) - Ci(2) * Bi(3)) / 2.
             CTot=Conc(i,1)+ Conc(j,1)+Conc(l,1)
             CSink(n,1) = CSink(n,1)  +  (Conc(i,1)  +  
     &             Conc(j,1)  +  Conc(l,1))/3.
    
            
		   DmE = DmE+  Dmol(1) * (ThNew(i) * Tau(hNew(i)) + 
     !	   ThNew(j) * Tau(hNew(j)) + ThNew(l) * Tau(hNew(l))) / 3.0     		       			
 		Enddo
      	Csink_M(n) = CSink(n,1)   ! mean of concentration for element              
		CSink(n,1) = CSink(n,1) * Sink(n) ! this is passive uptake
		                                  ! - if N uptake comes only with water
        SIncrSink2 = SIncrSink2 + cSink(n,1) * step * 14. / 62. * AE
		Disp(n,1) = DmE/(NUS-2.0) + DlngR(M) * VUP(n,1)
		Disp(n,2) = DmE/(NUS-2.0) + DlngR(M) * VUP(n,2)		
	 Enddo	 
	 if (iSink.eq.1) then        ! plant regulation
	 	 dC=0.1*NDemandError+0.05*cumulativeNDemandError 
	   ConstI(1)=amin1(75.0,ConstI(1)-dC)
	   ConstI(1)=amax1(ConstI(1),3.0)
	   ConstI(2)=ConstI(1)/2.0
	 ! change to 0 and wait for crop to be executed before
	 ! recalculating
	   NDemandError=0.0
	   cumulativeNDemandError=0.0
	   EndIf
C calculate nutrigen root uptake    
	 if (iSink.gt.0) Call massRootflux
c -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -     
       do n = 1,NumEl
         
         SIncrSink = SIncrSink + CSink(n,1) * Step* Area(n) *  14. / 62. ! go from Nitrate to Nitrogen
        
         WIncrSink = WIncrSink + Sink(n) * Area(n) * Step
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
10     Call errmes(im,il)
      end
     
      subroutine massRootflux()
      include 'public.ins'
      IncLude 'Puplant.ins'
      Include 'Puweath.ins'
	real * 8 qsinkC(NumElD,2),alphaK(NumElD,2)
	Real LastCr_M,F_MM(NumElD,2)
      integer e,Iroot(NumElD,2)
	Real *8 determ,b,a,c,bet,beta2,beta2a,alf,alf2,part1,part2,cr,betR,
     !   gammaB, gammaB2,Pc1, Pc2, pi, TwoPi,Ro2,AlphaKK(2)
      Common / Mass_root_M / Isink,RootRadius,
     !   DlngR(NMatD), Disp(NumElD,2),RootExchange(2),
     !	CSink_M(NumElD),Cr_M(NumElD,2), LastCr_M(NumElD,2) 
	Common  / PotNitr /  PotNitrogen_t
	  
	 pi = 3.14156
	 TwoPi = 6.28
	 pi4 = 12.56
	 Ro2 = RootRadius * RootRadius                 ! root radius (Ro)
	 alpha2=0.05
C FUP is root density, VUP is water uptake per unit root length
c calculate and alternative alphaK depending on demand
       

       do j=1,2
           RootLenFact=0.0
           do e=1,NumEL
           RootLenFact=RootLenFact+FUP(e,j)*Cr_M(e,j)
           enddo
           alphaKK(j)=PotNitrogen_t/(TwoPi*RootRadius*RootLenFact)
       enddo         
	 do j = 1,2
	  do e = 1, NumEl 
	   Cr_M(e,j)=LastCr_M(e,j)
	   do k=1,3
         Iroot(e,j) = 0                 ! 1 if length >0
	   qsinkC(e,j) = 0.0 
	   if(FUP(e,j).gt.0) THEN
	    alphak(e,j)=(ConstI(j))/(ConstK(j)+Cr_M(e,j)) 
	    F_MM(e,j)=alphak(e,j)*Cr_M(e,j)
	    Rx = 1.0 / sqrt(pi * (FUP(e,1) + FUP(e,2)))
	    RootRatio=Rx/RootRadius
	    
	    SC = 0.017  -  (PSIS(e) * 0.5)  
	    betR = SC / RootRadius
	    RootRatio = dmin1(betR,RootRatio)
	    RootRatio2=RootRatio*RootRatio
	     if (Disp(e,j).le.0.0) then
             GammaB=0
	      else
	      If ((VUP(e,j).le.0).OR.(iSink.EQ.3)) Then
	       GammaB=alphaK(e,j)*RootRadius/(Disp(e,j))
	       PC=1-0.5*GammaB
	       PC2=Rx*Rx*GammaB*alog(RootRatio)/(Rx*Rx-RootRadius**2)
	       PC=1.0/(PC+PC2)
	       PC=amax1(0.0,amin1(1.0,PC))
	         Else
	       GammaB=RootRadius*VUP(e,j)/Disp(e,j)
	       GammaB2=2.0/(2.0-GammaB)
	       PC=(RootRatio**(2.0-GammaB)-1.0)/(RootRatio**2-1.0)
	       PC=AlphaK(e,j)+(VUP(e,j)-AlphaK(e,j))*GammaB2*PC
	       if (PC.eq.0.0) then
	         PC=0.0
	         Else
                PC=amax1(0.0,amin1(1.0,VUP(e,j)/PC))	       
               EndIf
	       EndIF  !VUP <=0
	     EndIf   !Disp <= 0
	    if ((PC.lt.0).or.(PC.gt.1.0))then 
	    
	      iii=1
	      EndIf
	    qSinkC(e,j)=TwoPi*RootRadius
	    qSinkC(e,j)=qSinkC(e,j)*alphaK(e,j)*PC*Csink_M(e)	
	    Cr_M(e,j)=PC*CSink_M(e)         ! update value at root, may have to iterate      
	    Cr_M1=Cr_M(e,j)   
	    if (CSink_M(e).LT.amax1(0.01,Cr_M(e,j))) then
	       qSinkC(e,j)=0.0
	      endif
         End If   !FUP <0
      End Do   ! k
      End Do   ! j
      End Do   !e
      
	 do e = 1, NumEl
	  sumSinkN = 0.0
	  do j = 1,2
	   if (iSink.EQ.4) then
	   dif=amax1(0.0,CSink_m(e)-CMin0(j))
	     sumSinkN=sumSinkN + ConstI(j)*FUP(e,j)*dif
     !        /(ConstK(j)+dif)
          else
	   sumSinkN = sumSinkN + qsinkC(e,j)*FUP(e,j)
	   LastCr_M(e,j)=Cr_M(e,j)    
	   EndIF
	  enddo	  
	  If(CSink(e,1).gt.0 )then
c	    write(101,102)time, e, CSink(e,1),qSinkC(e,1)+qSinkC(e,2),
c     !	 (Cr_M(e,j), j = 1,2) , Csink_m(e), (ConstI(j),j=1,2),
c     !     (AlphaK(e,j),j=1,2),F_MM(e,1)+F_MM(e,2),SumSinkN
        continue
	  endif
	    CSink(e,1) = amax1(sumSinkN,0.0)
	 enddo
102    format(g14.8,1x,I5,12(2x,e10.3))
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
	


