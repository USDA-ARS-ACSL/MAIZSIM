      subroutine GasMover()
      Include 'public.ins'
      Include 'puweath.ins'
      Double precision A,A_1,B,C,P,Sum,CGFLux
      Real bTort
      Integer gKod,newjjj
      Dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD),
     !          A_1(MBandD,NumNPD)
          
      Common /GasM/ ThOld(NumNPD),ThANew(NumNPD),ThAOld(NumNPD),
     !          ThTot(NMatD),bTort,
     !          Dair(NumGD),
     !          NLevel,epsi, S(3,3),Ci(3),Bi(3),List(3),
     !          Diff(NumNPD),Ac(NumNPD),Gc(NumNPD),Sc(NumNPD),
     !          g_1(NumNPD,NumGD)   
      
      CGFLux=0.0
      If (lInput.eq.0) goto 11
C  Reading of the input files and initial calculations
      im=130
      il=0	 
      Open(40,file=ParamGasFile, status='old',ERR=10)	 
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) NumG                      !Number of gases
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) epsi                  
      Nlevel=1
      If(epsi.lt.0.999) NLevel=2      
      im=im+1
      il=il+1
      Read(40,*,ERR=10)     
      im=im+1
      il=il+1
      Read(40,*,ERR=10) 
CDT use one value for the horizon. if we can determine a soil texture dependence we can use that
C   to tailor the value to a specific layer
      Read(40,*,ERR=10) bTort
      im=im+1
      il=il+1
      Read(40,*,ERR=10)     
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      Read(40,*,ERR=10) (Dair(jjj),jjj=1,NumG)                    !Diffusion coefficients of gases in free air [cm2/day]
      Close(40)
C

*	
      Do i=1,NumNP
        ThOld(i)=ThNew(i)
      Enddo
      DtMx(4)=1.E+30

! Representing saturated water content as the total porosity of the soil
      Do i=1,NMat
      ThTot(i)= thSat(i)                                    !Thtot: total porosity equal to the saturated water content
      End do 
           
11    continue 
      tOld = Time
      t=Time
      dt = Step
      Do 711 jjj=1,NumG	 	 
        call DiffCoeff(jjj,NumNP,NMat,NumGD,NMatD,Diff,Tmpr,
     !           ThAOld,ThOld,ThTot,ThATr,bTort,ThAMin,Dair,MatNumN)  !Returns Diff, the diffusion coefficient in soil air [cm2/day]
      xMul=1.0D0
      alf=1.0D0-epsi
      newjjj=MBand             
c	 
        Do 13 i=1,NumNP
          ThetaAir(i)=ThAOld(i)
          DS(i)=0.            
          B(i) =0.                            ! Coefficient vector
          Gc(i)=0.
          Sc(i)=0.
          g_1(i,jjj)=g(i,jjj) 
          If(lOrt) B1(i)=g(i,jjj)             ! Nodal gas concentration in ug/cm3 air
          If(epsi.lt.0.001D0) then
              if(lOrt) newjjj = IADD(i)
              A(newjjj,i)=0.0D0
          Else

              Do 12 j=1,MBandD
                  A(j,i)=0.0D0
12            Continue
          Endif
13	  Continue         

	  

C
C  Recasting sinks
C  If we are using the unit of [microgram co2/cm3 of air] as the unit of co2 gas in soil,
c  the unit of gsink should be [microgram co2/cm3/day]    
C  No need of the element level calculation, we can get that direcly
c  from the  co2 node and node area or directly with unit of gSink as [ug co2/cm3 air/day] 

cSB The sink should be determined for both CO2 and o2 gas
cSB Sink for O2 is developed by considering a stoichiometric relationship with co2 as the O2 is the source of the oxygen in co2
        Do i=1,NumNP
          if (jjj.eq.1) gSink(i,jjj)  = gSink_root(i,1) + gSink_OM(i,1)

cSB gSink (i,2): (i,2) refers to o2 gas
cSB gSink (i,2): sum the o2 sink from root and microbial respiration
cSB gSink_root(i,2), gsink_OM(i,2): these are estimated in the Rootgrow() and SoilNitrogen() [ug o2/cm3 air/day]
cSB gSink(i,jjj) [ug o2/cm3 air/day]

          if (jjj.eq.2) gSink(i,jjj)= -(gSink_root(i,jjj)+ 
     !     gsink_OM(i,jjj))
          
          if (jjj.eq.3) gSink(i,jjj)  = gSink_N2O(i,1) 

           Gc(i)= gSink(i,jjj)                                   ![ug co2/cm3 air/day]  
        Enddo
          
        
        
c
c  Assembling matrixes
c
        Do 21 Level=1,NLevel
          If(Level.eq.NLevel) then
            Call DiffCoeff(jjj,NumNP,NMat,NumGD,NMatD,Diff,Tmpr,
     !      ThANew,ThNew,ThTot,ThATr,bTort,ThAMin,Dair,MatNumN)
            Do 14 i=1,NumNP
              Ac(i)=-(ThAOld(i)*alf+ThANew(i)*epsi)           !volumetric air content [cm3 air/cm3 soil]
              DS(i)=0.
14          Continue
          Endif
          Do 15 i=1,NumNP
            F(i)=0.
15        Continue
C
C  Loop on elements
C
          Do 19 n=1,NumEl
            NUS=4
            If(KX(n,3).eq.KX(n,4)) NUS=3
C
C     Loop on subelements
C
            Do 18 k=1,NUS-2
              i=KX(n,1)
              j=KX(n,k+1)
              l=KX(n,k+2)
              List(1)=i                                               !boundary nodes of the element 
              List(2)=j            
              List(3)=l
              Ci(1)=x(l)-x(j)                                         !x and y are axis of the 2d domain
              Ci(2)=x(i)-x(l)
              Ci(3)=x(j)-x(i)
              Bi(1)=y(j)-y(l)
              Bi(2)=y(l)-y(i)
              Bi(3)=y(i)-y(j)
              AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.                         !Area of the element 
              If(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.         !Modifying factor to transform equations of planar to axisymmetric flow, cm
              GcE=(Gc(i)+Gc(j)+Gc(l))/3.                              !Average gsink=Gc= [ugCo2/cm3 air/day]
              Ec=(Diff(i)+Diff(j)+Diff(l))/3.                         !Average diff [cm2/day]
              FMul=xMul*AE/4.
              SMul1=-1./AE/4.*xMul
              Do 17 j1=1,3                                            !For each of the nodes
                i1=List(j1)                                           !Node number
                F(i1)=F(i1)+FMul*(GcE+Gc(i1))/3.                      !GcE: Average nodal extraction, GC=gsink= [ugCo2/cm3/day]*[cm3] F(i1)=[ugCo2/day] Eq: 9.11 2DSoil DOC. 
                If(Level.eq.NLevel) DS(i1)=DS(i1)+AE*Ac(i1)/3.        !AE: Element area [cm2], Ac=[cm3air/cm3Soil], DS[cm3air]     : this may be Q matrix    
                Do 16 j2=1,3
                  i2=List(j2)
                  S(j1,j2)=SMul1*Ec*(Bi(j1)*Bi(j2)+Ci(j1)*Ci(j2))     !S[cm2/day], EC is diff[cm2/day]
                  If(i1.eq.i2) 
     !              S(j1,j2)=S(j1,j2)-AE*(ThNew(i1)-ThOld(i1))/3./dt  !This may be G in the matrix eqn 
                  If(Level.ne.NLevel) then
                    B(i1)=B(i1)-alf*S(j1,j2)*g(i2,jjj)
                  Else
                     If (lOrt) Then
                    call Find(i1,i2,kk,NumNp,MBandD,IAD,IADN)
                    ib=kk
                  Else
                     iB=MBand+i2-i1
                  Endif
                     A(iB,i1)=A(iB,i1)+epsi*S(j1,j2)                  
                      
                  Endif
16              Continue
17            Continue
18          Continue
19        Continue
          
          A_1=A 
          
          Do 20 i=1,NumNP
            M=MatNumN(i)
            If(Level.ne.NLevel) then
              B(i)=B(i)-alf*F(i)
            Else
              if (lOrt) newjjj=IADD(i)
              A(newjjj,i)=A(newjjj,i)+DS(i)/dt                !DS(i)/dt=[cm3 air/day]   A=[cm3 air/day]                                
              B(i)=B(i)+DS(i)/dt*g(i,jjj)-epsi*F(i)           !DS(i)/dt*g(i,jjj) gas in the nodes=[cm3air/day]*ugCo2/cm3air    F(i)sink=[ugCO2/day]   B=[ugCO2/day]
            Endif
20        Continue
21      Continue
C
C     Boundary condition
        Do 114 i=1,NumNP
          If(CodeG(i).ne.0) then                      !If not impermeable
            Do 111 j=1,NumBP
              If(KXB(j).eq.i) then
                If(CodeG(i).gt.0) then                !If code is +ve            
                  gKod=1
                  If(CodeG(i).eq.1) then              !if constant concentration at the bottom, code=1, This value of concetration is from the initial nodal concentration   
                    gBnd=g(i,jjj)                     !g(i,jjj):this is already converted to [ugco2/cm3air] from [ppm] in Grid boundary        
                  Else
                     gBnd=VarBG(j,jjj,1)              !gBnd=CO2 from weather file with boundary code as 3  or 4  or 6 [ugco2/cm3air] 
                  Endif
                Else
                  gKod=3                              !Code is negative
                  If(CodeG(i).eq.-6) then             !-6: Time dependent gas flux within soil domain
                    Qg=VarBG(j,jjj,3)                 ![ugco2/cm3air  *  cm/day]              !VarBG(j,jjj,3)=PG*GAIR(jjj)/509.6=[ugco2/cm3 air * cm/day]   
                    Pg=VarBG(j,jjj,2)                 ![cm/day]                               !VarBG(j,jjj,2)=PG[cm/day]
                  Else                                !-4: Soil atmospheric boundary with gas flux prescribed[ugco2/cm3air] it can be codeG=-3 also
                    Qg=VarBG(j,jjj,3)*Width(j)        ![ugCo2/cm3 air  *  cm/day] * width [cm]        !VarBG(j,jjj,3)=PG*GAIR(jjj)/509.6=[ugco2/cm3 air * cm/day]        
                    Pg=VarBG(j,jjj,2)*Width(j)        ![cm/day] * width [cm]                          !VarBG(j,jjj,2)=PG[cm/day]         
                    
                    !QGas(j,jjj)=(Pg*g(j,jjj))-Qg      ![pg:cm/day] * [width:cm] *[ug co2/cm3air] = [ugCo2/day?], if positive flux goes out from soil to air    
                  Endif
                Endif
                Goto 112
              Endif
111         Continue
C
C     Dirichlet boundary condition
C
112      If(gKod.eq.1) then
            If(lOrt) then
	      A(IADD(i),i) = 1.d30
	      B(i) = 1.d30*gBnd                         !gBnd Constant concentration at the botttom [ugco2/cm3air]
          Else
              Do 113 j=1,2*MBand-1
                A(j,i)=0.
113           Continue
              A(MBand,i)=1.
              B(i)=gBnd
              end if
            Endif

*     Cauchy boundary condition
            If(gKod.eq.3) then                    !CodeG <0
              B(i)=B(i)-(Qg-alf*Pg*g(i,jjj))     
              if(lOrt) newjjj=IADD(i)
          A(newjjj,i)=A(newjjj,i)-EPSI*Pg           
            Endif
          Endif
114     Continue
        
C
c  Solve the global matrix equation for transport
     
*  Explicit formulation
        If(epsi.lt.0.001) then
          Do 22 i=1,NumNP
              if(lOrt) newjjj=IADD(i)
            B(i)=B(i)/A(newjjj,i)
22        Continue
           else if(lOrt) then
	       call ILU(A,NumNP,MBandD,IAD,IADN,IADD,A1)
             call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,
     !	            IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,4,
     !                MNorth,MaxItO,4)

        Else
c  Implicit formulation
          N1=NumNP-1
          Do 212 k=1,N1
            P=1./A(MBand,k)
            kk=k+1
            kc=MBand
            Do 211 i=kk,NumNP
              kc=kc-1
              If(kc.le.0) goto 212
              C=-P*A(kc,i)
              A(kc,i)=C
              ii=kc+1
              L=kc+MBand-1
              Do 211 j=ii,L
                jj=j+MBand-kc
                A(j,i)=A(j,i)+C*A(jj,k)
211         Continue
212       Continue
          Do 214 i=2,NumNP
            jj=MBand+1-i
            ii=1
            If(jj.le.0) then
              jj=1
              ii=i-MBand+1
            Endif
            Sum=0.
            Do 213 j=jj,MBand-1
              Sum=Sum+A(j,i)*B(ii)
              ii=ii+1
213         Continue
            B(i)=B(i)+Sum
214       Continue              
          B(NumNP)=B(NumNP)/A(MBand,NumNP)
          Do 216 k=1,N1
            i=NumNP-k
            jj=i
            m=min0(2*MBand-1,MBand+k)
            Sum=0.
            Do 215 j=MBand+1,m
              jj=jj+1
              Sum=Sum+A(j,i)*B(jj)
215         Continue
            B(i)=(B(i)-Sum)/A(MBand,i)          
216       Continue  
        Endif
*	
       Do i=1,NumNP
	  if(lOrt) B(i) = B1(i)
                  if (B(i).lt.0) then
                 B(i)=0.0
                 iii=1
            endif
        g(i,jjj)=sngl(B(i))     
       Enddo      
*       
       

!  Gas FLux calculation from the matrix
!  This should work for -4, 1, 3, 6 
      
      Do  n=1,NumNP
          If ((CodeG(n).eq.-4).or.(CodeG(n).eq.1)) then
              QG =F(n)+(DS(n)*(g(n,jjj)-g_1(n,jjj))/dt )
              do   j=1,IADN(n)
                  QG=QG+A_1(j,n)*g(IAD(j,n),jjj)       
              end do
              QGas(n,jjj)=QG
          End if
      
      end do
      
711   Continue
*
      Do i=1,NumNPD
        ThOld(i)=ThNew(i)
      Enddo
      
      
      Return
10    Call errmes(im,il)
      Return
      End
* 
      Subroutine DiffCoeff(jjj,NumNP,NMat,NumGD,NMatD,Diff,Tmpr,
     !           ThetaA,theta,ThTot,ThATr,bTort,ThAMin,Dair,MatNum)

      Dimension Dair(NumGD),ThTOt(NMat),ThATr(NMat),
     !          ThAMin(NMat),Diff(NumNP),MatNum(NumNP),
     !          Tmpr(NumNP),theta(NumNP),ThetaA(NumNP)
      
 

      Do 11 i=1,NumNP
        M=MatNum(i)
       ThetaA(i)=Amax1(ThTot(M)-theta(i)+ThAMin(M),ThAMin(M))        ! Get current air filled porosity, ThAmin is minimum value of airfilled pososity
        tau=Amax1(bTort*(ThetaA(i)-ThATr(M)),0.)                   ! bTort is toruosity change/unit of porosity change,  ThAtr is Threshold porosity, tau is tortousity factor
        Diff(i)=ThetaA(i)*Dair(jjj)*tau*(Tmpr(i)/273.15+1.)**2        ! Diff is nodal value of gas diffusion coefficient in soil air[cm2 soil/day], Dair is Diffusion coefficients of gases in free air [cm2/day]
11    Continue
      Return
      End