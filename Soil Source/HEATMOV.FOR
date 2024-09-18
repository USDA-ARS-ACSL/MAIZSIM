      subroutine HeatMover ()
      Include 'public.ins'
      Double precision A,B,C,P,Sum
      Real VSind, VClay, VOrgM, ThTot
      Integer tKod
      Dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD) 
      Common /HeatM/ ThOld(NumNPD),
     !                NLevel,dt,
     !                S(3,3),Ci(3),Bi(3),List(3),
     !                VxOld(NumNPD),VzOld(NumNPD),
     !                VxH(NumNPD),VzH(NumNPD),hOld(NumNPD),
     !                WeTab(3,2*NumElD),ThTot(NMatD),
     !                Ac(NumNPD),Fc(NumNPD),Qc(NumNPD),
     !                CWat,TCSzz(NumNpd),HCS(NumNPD),
     !                VSind(NMatD), VClay(NMATD),VOrgM(NMatD),
     !                RNS,Sc(NumNPD),Gc(NumNPD)
      If (lInput.eq.0) goto 11
C
C  Initialization here
C          
      Nlevel=1
      Do i=1,NMat
        VSind(i)=BlkDn(i)*FracSind(i)/2.65
        VClay(i)=BlkDn(i)*FracClay(i)/2.65
        VOrgM(i)=BlkDn(i)*FracOm(i)/1.0
        ThTot(i)=1-BlkDn(i)/2.65
      Enddo     
*       
      Do i=1,NumNP
        VxOld(i)=Vx(i)
        VzOld(i)=Vz(i)
        ThOld(i)=ThNew(i)
        hOld(i)=hNew(i)
      Enddo
* units are in Joules, not calories. This is the value of the 
* specific heat of water in J g-1 oC-1    
      cWat=1.0/.239
c volumetric heat capacity of water j/ cm3 same as specific heat

      Movers(3)=1
      Return
C
C  Routine calculations
C
11    Continue 
      tOld = Time
      t=Time
      dt = Step
      Call Thermal(HCS,MatNumN,NMatD,NumNP,VSind,VClay,VOrgM,ThNew,
     !               EPSI_Heat,ThTot,CWat,Tmpr,TcSxx,TCSzz,hCritA,hnew)
      xMul=1.
      alf=1.-EPSI_Heat
c
C      GR 1 line   
      newjjj = MBand
      Do 13 i=1,NumNP
C for the internal nodes A holds the heat capacity of the soil and the heat capacity of the water      
        Ac(i)= -(HCS(i)+cWat*(ThOld(i)*Alf+EPSI_Heat*Thnew(i)))
        Ds(i)=0.
        Fc(i)=0.
        Sc(i)=0.
        B(i) =0.
        Gc(i)=0.
C      GR
      if(lOrt) B1(i)= tmpr(i)

       If(EPSI_Heat.lt.0.001) then
          A(MBand,i)=0
        Else
          Do 12 j=1,MBandD
            A(j,i)=0.
12        Continue
        Endif
13    Continue
C
C   Calculate sinks 
C
cccz directly take the sink and nodearea 
C need to add cWat here
      Do n=1,NumNP
         Fc(n)=Sink(n)*cWat
         Sc(n)=NodeArea(n)
      Enddo
C
C Assembling matrixes
C
      Do 21 Level=1,NLevel

        If(Level.eq.NLevel) then
          Do i=1,NumNP
            VxH(i)=Vx(i)
            VzH(i)=Vz(i)
          Enddo
          Call Thermal(HCS,MatNumN,NMatD,NumNP,VSind,VClay,VOrgM,ThNew,
     &    EPSI_Heat,ThTot,CWat,Tmpr,TcSxx,TCSzz,hCritA,hnew)
          Do n=1,NumEl
          NUS=4
            If(KX(n,3).eq.KX(n,4)) NUS=3
            Do k=1,NUS-2
              i=KX(n,1)
              j=KX(n,k+1)
              l=KX(n,k+2)
              Ci(1)=x(l)-x(j)
              Ci(2)=x(i)-x(l)
              Ci(3)=x(j)-x(i)
              Bi(1)=y(j)-y(l)
              Bi(2)=y(l)-y(i)
              Bi(3)=y(i)-y(j)
              VxE=abs(Vx(i)+Vx(j)+Vx(l))/3.
              VzE=abs(Vz(i)+Vz(j)+Vz(l))/3.
            Enddo
          Enddo
c
        Else
          Do i=1,NumNP
            VxH(i)=VxOld(i)
            VzH(i)=VzOld(i)
          Enddo
        Endif
        Do 14 i=1,NumNP
          M=MatNumN(i)
          If(Level.ne.NLevel) then
            DPom = dt*cWat**2/6./(-Ac(i))       !CHECK UNITS
            TcSxx(i)=TcSxx(i)+VxH(i)*VxH(i)*DPom
            TcSzz(i)=TcSzz(i)+VzH(i)*VzH(i)*DPom
          Else
            DPom=dt*cWat**2/6./(-Ac(i))
            TcSxx(i)=TcSxx(i)-VxH(i)*VxH(i)*DPom
            TcSzz(i)=TcSzz(i)-VzH(i)*VzH(i)*DPom
            Fc(i)=0.
          Endif
14      Continue
        Do 15 i=1,NumNP
          F(i)=0.
          if(Level.eq.NLevel) DS(i)=0.
15      Continue
C
C  Loop on elements
C
        NumSEl=0
        Do 19 n=1,NumEl
          NUS=4
          If(KX(n,3).eq.KX(n,4)) NUS=3
*         Loop on subelements
          Do 18 k=1,NUS-2
            NumSEl=NumSEl+1
            i=KX(n,1)
            j=KX(n,k+1)
            l=KX(n,k+2)
            List(1)=i
            List(2)=j
            List(3)=l
            Ci(1)=x(l)-x(j)
            Ci(2)=x(i)-x(l)
            Ci(3)=x(j)-x(i)
            Bi(1)=y(j)-y(l)
            Bi(2)=y(l)-y(i)
            Bi(3)=y(i)-y(j)
            AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
            if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
            GcE= -tSink(n)
            Ec1=(TcSxx(i)+TcSxx(j)+TcSxx(l))/3.
            Ec3=(TcSzz(i)+TcSzz(j)+TcSzz(l))/3.
            VxE=(VxH(i)+VxH(j)+VxH(l))*cWat/3.
            VzE=(VzH(i)+VzH(j)+VzH(l))*cWat/3.
            AcE = (Ac(i)+Ac(j)+Ac(l))/3.
            FcE=(Fc(i)+Fc(j)+Fc(l))/3.
            GcE=(Gc(i)+Gc(j)+Gc(l))/3.
            FMul=xMul*AE/4.
            SMul1=-1./AE/4.*xMul
            SMul2=AE/20.*xMul
            Do 17 j1=1,3
              i1=List(j1)
              F(i1)=F(i1)+FMul*(GcE+Gc(i1))/3.
              If(Level.eq.NLevel) DS(i1)=DS(i1)+AE/4.0*(AcE+Ac(i1)/3.)
              Do 16 j2=1,3
                i2=List(j2)
                S(j1,j2)=SMul1*(Ec1*Bi(j1)*Bi(j2)+Ec3*Ci(j1)*Ci(j2))
                S(j1,j2)=S(j1,j2)-(Bi(j2)/8.*(VxE+(VxH(i1)*cWat)/3.)+
     !                    Ci(j2)/8.*(VzE+(VzH(i1)*cWat)/3.))*xMul
                ic=1
                if(i1.eq.i2) ic=2
                S(j1,j2)=S(j1,j2)+SMul2*ic*(FcE+(Fc(i1)+Fc(i2))/3.)
                If(Level.ne.NLevel) then
                  B(i1)=B(i1)-alf*S(j1,j2)*Tmpr(i2)
                Else
C     GR
	           if(lOrt) then
	              call Find(i1,i2,kk,NumNP, MBandD, IAD, IADN)
	               iB=kk
	           else
	                iB=MBand+i2-i1
	           endif
                  A(iB,i1)=A(iB,i1)+EPSI_Heat*S(j1,j2)
                Endif
16            Continue
17          Continue
18        Continue
19      Continue

        Do 20 i=1,NumNP
	    if (Level.ne.NLevel) then
               B(i)=B(i)-alf*F(i)	 
	else
	    if (lOrt) newjjj=IADD(i)
             A(newjjj,i)=A(newjjj,i)+DS(i)/dt           
             b(i)=b(i)+ds(i)/dt*Tmpr(i)-EPSI_Heat*F(i)
          Endif
		 
20      Continue
21    Continue
C
C     Boundary conditions
C

      Do 114 i=1,NumNP
      If(CodeT(i).ne.0) then
        Do 111 j=1,NumBP
          If(KXB(j).eq.i) then
            If(CodeT(i).gt.0) then
              tKod=1
              If(CodeT(i).eq.1) then 
                tBnd=tmpr(i)
              Else
                tBnd=VarBT(j,1)
              Endif
            Else
              tKod=3
              If(CodeT(i).eq.-6) then
                Qh=VarBT(j,3)
                Ph=VarBT(j,2)
              Else
                Qn=0.0
c account here for the increase in heat going into soil when
c evaporation proceeds at a rate less than potential and extra
c heat goes into the soil
c The difference is the heat of vaporization value of the water
c not evaporated.
                 If(Q(i).lt.0.0) then
                    Qn=VarBT(j,4)*width(j)-abs(Q(i))*
     &                  (2500.8-2.3668*VarBT(j,1))
                  else
                    If(Q(i).ge.0) then
c Z addition - account for the heat of vaporization of the part of the infiltrating water that evaporates
                    Qn=VarBT(j,4)-VarBW(j,2)*(2500.8-2.3668*VarBT(j,1))
                    Qn=Qn*width(j) !units are J day-1
                     endif ! end of if Q(i) ge 0
                endif ! end of if Q(i) le 0
c Ph and Qh are two components of convected heat, Qh has air temperature
c   incorporated into it. there is no convection in the case of ponded water
                if (hnew(i).le.0.000) then
                  Qh=VarBT(j,3)*Width(j)+Qn
                  Ph=VarBT(j,2)*Width(j)
                Endif
              Endif ! end of if CodeT(i) eq -6
            Endif ! end of if CodeT(i) gt 0
            goto 112
          Endif ! end of if KXB(j) eq i
111     Continue ! end of loop on boundary points
*     Dirichlet boundary condition
112     If(tKod.eq.1) then
          If(lOrt) then
	      A(IADD(i),i) = 1.d30
	      B(i) = 1.d30*tBnd
          Else
           Do 113 j=1,2*MBand-1
             A(j,i)=0.
113        Continue
           A(MBand,i)=1.
           B(i)=tBnd
	    Endif
        Endif
*     Cauchy boundary condition
        If(tKod.eq.3) then
          B(i)=B(i)-(Qh-alf*Ph*tmpr(i))
	    if(lOrt) newjjj=IADD(i)
          A(newjjj,i)=A(newjjj,i)-EPSI_Heat*Ph
        End if
      End if ! end of if codeT(i) ne 0
114   Continue ! end of loop on nodes
c   pond   on the soil-atmosphere surface
cMisha 18/9 2006
c dt 3/31/2011 I added width to the A and B matrices
C  SURFACE POND FOR TRANSPORT            
          Do 313 i=1,NumBP
          n=KXB(i)
          k=CodeW(n)
         If (iabs(k).eq.4) then		
C in the case of ponded infiltration, the amount of water that infiltrated (in cm) is the difference in head
C between the two time steps times the width of the element (total amount for the element)
C The units of A() is J/day/C and the units aof B() is J/day 
C then when you get the solution from B/A the units are Temperature (C)

			HNEWS=DMAX1(hNew(N),0.0D0)
		    HOLDS=DMAX1(hold(N),0.0D0)
		    newjjj=1	
		    if(lOrt) newjjj=IADD(n)   
              if (((HNEWS.GT.0.0).AND.(HOLDS.GT.0.0))
     &            .AND.(QAct(n).GT.0.0)) then
                   qtemp=(hnews-holds)/dt*width(i)
                   A(newjjj,n)=A(newjjj,n)-EPSI_Heat*
     &                QAct(n)*cwat*width(i)
                   B(n)=B(n)-(varBT(i,1)- tmpr(n)*alf)*
     &                QAct(n)*cWat*width(i)
              endif
          END IF
 313	  CONTINUE	 
C
c   pond   on the soil-atmosphere surface
cMisha 18/9 2006
C
C

C
C Solve the global matrix equation for transport
C
      If(EPSI_Heat.lt.0.001) then
	  
        Do 22 i=1,NumNP
	  if(lOrt) newjjj=IADD(i)
          B(i)=B(i)/A(newjjj,i)
22      Continue
C      Else
cSolve  
      else if(lOrt) then
	       call ILU(A,NumNP,MBandD,IAD,IADN,IADD,A1)
             call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,
     !	            IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,4,
     !                MNorth,MaxItO,3)
      Else
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
211           Continue
212       Continue
          Do 214 i=2,NumNP
            jj=MBand+1-i
            ii=1
            If(jj.le.0) then
              jj=1
              ii=i-MBand+1
            End if
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
      End if

      Do i=1,NumNP
	  if(lOrt) B(i) = B1(i)
        Tmpr(i)=sngl(B(i))
      Enddo      
*       
      
      
        VxOld(:)=Vx(:)
        VzOld(:)=Vz(:)
        ThOld(:)=ThNew(:)
        Hold(:)=Hnew(:)
      
      
      Return
      End
*
      Subroutine Thermal(HCS,MatNumN,NMatD,NumNP,VSind,
     !               VClay,VOrgM,ThNew,
     !               EPSI_Heat,ThTot,CWat,Tmpr,Tcsxx,TCSzz,hCritA,hNew)
C calculate thermal conducivity
      Dimension HCS(NumNP),MatNumN(NumNP),VSind(NMatD),VClay(NMatD),
     !          VOrgM(NMatD),ThNew(NumNP),ThTot(NMatD),Tmpr(NumNP),
     !          Tcsxx(NumNP),TCSzz(NumNP),hNew(NumNP)

      alf=1.0-EPSI_Heat
C     ---- THERMAL CONDUCTIVITY OF WATER SATURATED AIR ----
C  --- HCS is volumetric heat capacity of the mineral part of the soil soil J cm3-1 'C-1
C   using 0.46 cal/cm3/C for mineral and 0.6 for OM
C 4.186 converts calories to joules
      Do 10, i=1,NumNP
        m=MatNumN(i)
         HCS(i)=1.9225*(VSind(m)+VClay(m)) + 2.51*VOrgM(m)
        TCH2O=1.33+4.4E-03*tmpr(i) ! milical cm-1 s-1 Value checks OK
** calculate the thermal conduc of air at a temp of 30 degrees C
        TCAIR=0.058+1.7E-04*30.0  !milical cm-1 s-1 (value checks ok)
        TCVAP=0.052*EXP(0.058*tmpr(i))
        TCSAT=TCAIR + TCVAP*ThNew(i)/ThTot(m)
C
c  calc thermal conduc for air dry soil
        If ((hnew(i)).le.(hCritA+100)) then
C         ---- INTERMEDIATE PARAMETERS ----
          XSIND=((2.0/(1.0+(((20.4/TCair)-1.0)*0.125)))
     !      + (1.0/(1.0+(((20.4/TCair)-1.0)*0.75))))/3.
          XCLAY=((2.0/(1.0+(((7.0/TCair)-1.0)*0.125)))
     !      + (1.0/(1.0+(((7.0/TCair)-1.0)*0.75))))/3.
          XMUCK=((2.0/(1.0+(((0.6/TCair)-1.0)*0.5)))+1.0)/3.
          
C
C     ----- THERMAL CONDUCTIVITY (milliCAL/CM/SEC/'C) ----
          Tcsxx(i)=(((XSIND*VSind(m)*20.4)
     *      + (XCLAY*VCLAY(m)*7.0) + (XMUCK*VOrgM(m)*0.6)
     *      + (ThTot(m)*TCair))
     *      /((XSIND*VSind(m)) + (XCLAY*VCLAY(m))
     *      + (XMUCK*VOrgM(m)) + ThTot(m)))*1.25
        Else
C
C       ---- INTERMEDIATE PARAMETERS ----
C
          XSIND=((2.0/(1.0+(((20.4/TCH2O)-1.0)*0.125)))
     !        + (1.0/(1.0+(((20.4/TCH2O)-1.0)*0.75))))/3.
          XCLAY=((2.0/(1.0+(((7.0/TCH2O)-1.0)*0.125)))
     !        + (1.0/(1.0+(((7.0/TCH2O)-1.0)*0.75))))/3.
          XMUCK=((2.0/(1.0+(((0.6/TCH2O)-1.0)*0.5)))+1.0)/3.
C
          XGAIR=0.3333*ThNew(i)/ThTot(m)
          IF(XGAIR.GT.0.3333) XGAIR=0.3333
          XAIR=((2.0/(1.0+(((TCSAT/TCH2O)-1.0)*XGAIR))) +
     !      (1.0/(1.0+(((TCSAT/TCH2O)-1.0)*
     !      (1.0-2.0*XGAIR)))))/3.0
C
C     ----- THERMAL CONDUCTIVITY (10-3 CAL/CM/SEC/'C) ----
C
          Tcsxx(i)=((Thnew(i)*TCH2O) + (XSIND*VSind(m)*20.4)
     *    + (XCLAY*VCLAY(m)*7.0) + (XMUCK*VOrgM(m)*0.6)
     *    + (XAIR*(ThTot(m)-ThNew(i))*TCSAT))
     *    /(ThNew(i)+(XSIND*VSind(m)) + (XCLAY*VCLAY(m))
     *    + (XMUCK*VOrgM(m)) + (XAIR*(ThTot(m)-ThNew(i))))
        Endif
CDT debugging -check to see if 4.186 is needed, not sure it is.
C 0.004184 J per milical 
        Tcsxx(i)=TcSxx(i)*0.004184*24.*3600
        Tcszz(i)=Tcsxx(i)
10      Continue
C     -----------
      Return
      End
