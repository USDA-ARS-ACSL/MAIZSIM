******************************************************************
cdt sometimes, if head was near 1000, the program would create an array
c index for contab that was 101, instead of 100. I think this was 
c because of roundoff errors. I therefore subtracted 0.04 from
c the threshold value to force a number such as -999.9999 to act
c like it was more than 1000 and so no go into the table.

Cdt 11/29/2011 added ability to calculate available water
C added two new variables TUpperLimit, TLowerLimit and an array, TAVAIL,
C   and ThFUll = which is total available water
C dt 7/2014 added Bulk Density and soil texture to soil file. They are 
C  read here
C 12/ 2015 made TUpperLimit and TLowerLimit arguments so they could be used in the main program
C 
      subroutine SetMat(lInput,NumNP,hNew,hOld,NMat,
     &                   MatNum,Con,Cap,
     &                  BlkDn,hTemp,Explic,Theta,hTab1,hTabN,
     &                  hSat,ThSat,ThR, ThAvail,ThFull,
     &                  FracOM, FracSind, FracClay,
     &                  TupperLimit,TLowerLimit,soilair,SoilFile,
     &                  ThAMin,ThATr)
      Logical Explic
 !InitType allows us to input water content(w) or head (m) as an initial condition
      character SoilFile*132
      ! or character*1 InitType(5)
      Real BldDn,soilair,ThAMin,ThATr
      Parameter (NTabD=100,NPar=13,NMatD=15)
      character*1, dimension(NMatD)::InitType

      Common /HydPar/ Par(NPar,NMatD),
     !                hTab(NTabD),ConTab(NTabD,NMatD),
     !                CapTab(NTabD,NMatD),ConSat(NMatD),
     !                TheTab(NTabD,NMatD),alh1,dlh
   
      Dimension hNew(NumNP),hOld(NumNP),MatNum(NumNP),
     !          Con(NumNP),Cap(NumNP),hTemp(NumNP),Theta(NumNP),
     !          thR(NMatD),hSat(NMatD),thSat(NMatD),ThAvail(NumNP),
     !           ThFull(NMatD),BlkDn(NMatD),
     !          FracOM(NMatD), FracSind(NMatD), FracClay(NMatD),
     !          TUpperLimit(NMatD),TLowerLimit(NMatD),soilair(NumNP),
     !          ThAMin (NMatD),ThATr(NMatD)
      
c   parameters:
c         1  thr
c         2  ths
c         3  tha
c         4  th
c         5  Alpha
c         6  n
c         7  ks
c         8  kk
c         9  thk
c         10  BD
         
      If(lInput.eq.1) then
        im=60
        il=0
      Open(40,file=SoilFile,status='old',ERR=16)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read (40,*,ERR=10)
        im=im+1
        Do 11 M=1,NMat
          il=il+1
          Read (40,*,ERR=10) (Par(i,M),i=1,NPar),InitType(M) 
11      Continue
        Close(40)
        hTab1=-amin1(abs(hTab1),abs(hTabN))
        hTabN=-amax1(abs(hTab1),abs(hTabN))
c
        alh1=alog10(-hTab1)
        dlh=(alog10(-hTabN)-alog10(-hTab1))/(NTabD-1)
        Do 12 i=1,NTabD
          alh=alh1+(i-1)*dlh
          hTab(i)=-10**alh
12      Continue
        hTab(NTabD)=hTabN
c calculate MP from water content if wc is given as an initial
        do n=1, numNP
           if (InitType(MatNum(n)).eq."w") then  ! need to calculate hnew thnew is ok
             mm=MatNum(n)
             VG_m=(1-1/(par(6,mm)))*(-1)
             thetaRel=(hNew(n)-par(1,mm))/(par(2,mm)-par(1,mm))
             hNew(n)=(ThetaRel**(1/Vg_m))**(1/par(6,mm))-1
             hNew(n)=hNew(n)/par(5,mm) * (-1)
           endif
         enddo
         Do n=1,NumNP
            hTemp(n) = hNew(n)  ! need to assign here since we are creating new values of hnew
            hOld(n)=hTemp(n)
          Enddo
        Do 14 M=1,NMat
          Hr          =FH(0.0,Par(1,M))
          hSat(M)     =FH(1.0,Par(1,M))
          ConSat(M)   =FK(0.0,Par(1,M))
          thR(M)      =FQ(Hr ,Par(1,M))
          thSat(M)    =FQ(0.0,Par(1,M)) 
          
          
CDT FracSind is sand + silt. Since we don't input silt, we 
C     get it by difference
          BlkDn(M)=  Par(10,M)
          FracOM(M)   =Par(11,M)
c          FracSind(M) =1.0-(Par(12,M)+Par(13,M))
c          FracClay(M) =Par(13,M)
          FracSind(M) = Par(12,M)+Par(13,M) 
          FracClay(M) = 1-Par(12,M)-Par(13,M) 
          TUpperLimit(M)=FQ(-100.0,Par(1,M))
          TLowerLimit(M)=FQ(-150000.0,Par(1,M))
          ThFull(M)=TUpperLimit(M)-TLowerLimit(M)
          if (ThFull(M).lt.0) then
             iii=0;
            endif
          
          
          Do 13 i=1,NTabD
            ConTab(i,M)=FK(hTab(i),Par(1,M))
            CapTab(i,M)=FC(hTab(i),Par(1,M))
            TheTab(i,M)=FQ(hTab(i),Par(1,M))
           
13        Continue
          
          ThAMin(M)=0.01                           !Entrapped air, value can be estimated as a function of soil properties
          ThATr(M)=0.00                           !Threshold air filled porosity below which diffusion does not occur, value can be estimated as a function of soil properties
           
          
14      Continue
      Endif  !lInput=1
      Do 15 i=1,NumNP
        M=MatNum(i)
        hi1=amin1(hSat(M),hTemp(i))
        hi2=amin1(hSat(M), hNew(i))
        If(Explic) hi2=hi1
        hiM=0.1*hi1+0.9*hi2
        If(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Ci=ConSat(M)
        Else if(hiM.ge.(hTab(NTabD)+0.04).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          S1=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=ConTab(iT,M)+S1*(hiM-hTab(iT))
        Else
          Ci=FK(hiM,Par(1,M))
        Endif
        Con(i)=Ci
        hi1=hOld(i)
        hi2=hNew(i)
        If(Explic) hi2=hi1
        If(hi2.ge.hSat(M)) then
          Ci=0
          Ci=FC(hSat(M),PAR) !prevent capacitance from going to 0
          Ti=thSat(M)
        Else if(hi2.ge.(hTab(NTabD)+.04).and.hi2.le.hTab(1)) then
          iT=int((alog10(-hi2)-alh1)/dlh)+1
          S1=(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
          S2=(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=CapTab(iT,M)+S1*(hi2-hTab(iT))
          Ti=TheTab(iT,M)+S2*(hi2-hTab(iT))
        Else
          Ci=FC(hi2,Par(1,M))
          Ti=FQ(hi2,Par(1,M))
        Endif
        Cap(i)=Ci
        theta(i)=thR(M)+(Ti-thR(M))

csb soil air  equals total porosity (or saturated water content)- current water content
        soilair(i)= AMAX1(thSat(M)-theta(i),0.0001)
	
C calculate available water content
        If (theta(i).lt.TLowerLimit(M)) Then 
          ThAvail(i)=0;     
          else 
            ThAvail(i)=Amin1(theta(i),TUpperLimit(M))
     &                       -TLowerLimit(M)
            ThAvail(i)=Amax1(ThAvail(i),0.0)
         End If          
15    Continue
      Return
16    Write(*,*) 'Soil file not found'      
10    Call errmes(im,il)
      Return
      End

************************************************************************
      real function FK(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk
      real h,Par(12)
c      double precision PPar
      integer PPar
      
      Kr=1.0D0

      BPar=0.5D0
      PPar=2.0D0
      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      Ks=Par(7)
      Kk=Par(8)
      Qk=Par(9)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      Hk=-1./Alfa*(Qeek**(-1./m)-1.)**(1./n)
      if(h.lt.Hk) then
        Qee=(1.+(-Alfa*HH)**n)**(-m)
        Qe =(Qm-Qa)/(Qs-Qa)*Qee
        Qek=(Qm-Qa)/(Qs-Qa)*Qeek
        FFQ =1.-(1.-Qee**(1./m))**m
        FFQk=1.-(1.-Qeek**(1./m))**m
        if(FFQ.le.0.d0) FFQ=m*Qee**(1./m)
        Kr=((Qe/Qek)**Bpar)*((FFQ/FFQk)**PPar)*Kk/Ks
        FK=max(Ks*Kr,1.d-37)
        return
      endif
      if(h.ge.Hk.and.h.lt.Hs) then
        Kr=(1.-Kk/Ks)/(Hs-Hk)*(h-Hs)+1.
        FK=Ks*Kr
      endif
      if(h.ge.Hs) FK=Ks
      return
      end

************************************************************************

      real function FC(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(12)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      if(h.lt.Hs) then
        C1=(1.+(-Alfa*HH)**n)**(-m-1.)
        C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.)*C1
        FC=max(C2,1.d-37)
        return
      else
        FC=0.0001
       ! write(*,*) 'h > hs', h, Hs
      endif
      return
      end

************************************************************************

      real function FQ(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(12)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
      if(h.lt.Hs) then
        Qee=(1.+(-Alfa*HH)**n)**(-m)
        FQ=max(Qa+(Qm-Qa)*Qee,1.d-37)
        return
      else
        FQ=Qs
      endif
      return
      end

************************************************************************

      real function FH(Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real Qe,Par(12)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.-1./n
      HMin=-1.d300**(1./n)/max(Alfa,1.d0)
      QeeM=(1.+(-Alfa*HMin)**n)**(-m)
      Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
      FH=max(-1./Alfa*(Qee**(-1./m)-1.)**(1./n),-1.d37)
      return
      end

************************************************************************
