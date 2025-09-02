      subroutine SoluteMover()
      Include 'public.ins'
      Include 'puplant.ins'
      Double precision A,B,C,P,Sum,dt1,dt2
      Integer cKod
      Dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD),
     &            Gc(NumNPD),Sc(NumNPD)
      Common /SoluteM/ ThOld(NumNPD),
     !       Dmol(NumSD),DLng(NMatD,NumSD),Dtrn(NMatD,NumSD),
     !                NLevel,dt,epsi,CourMax,lUpW,
     !                S(3,3),Wz(3),Wx(3),Ri(3),Ci(3),Bi(3),List(3),
     !                VxOld(NumNPD),VzOld(NumNPD),
     !       VxH(NumNPD),VzH(NumNPD),hOld(NumNPD),ConOld(NumNPD),
     !                Dispzz(NumNPD),Dispxx(NumNPD),Dispxz(NumNPD),
     !                WeTab(3,2*NumNPD)
      If (lInput.eq.0) goto 11
C 
C  Reading of the input files and initial calculations 
C
      im=70
      il=0
      Open(40,file=SoluteFile, status='old',ERR=10)   
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) NumSol
c      im=im+1
c      il=il+1
c      Read(40,*,ERR=10)
c      im=im+1
c      il=il+1
c      Read(40,*,ERR=10) NL
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) epsi,lUpW,CourMax!explicit if epsi =0
      Nlevel=1
      if(epsi.lt.0.999) NLevel=2      
      im=im+1
      il=il+1
      Read(40,*,ERR=10)     
      im=im+1
      il=il+1
      Read(40,*,ERR=10)     
      im=im+1
      il=il+1
      Do jjj =1,NumSol
         Read(40,*,err=10)SN, Dmol(jjj) !modelcular diffusion of each solute
         il=il+1
      enddo
      im=im+1
      il=il+1
      Read(40,*,ERR=10)     
      Do NS=1,NumSol
          Do M= 1,NMat
            il=il+1
            Read(40,*,ERR=10) SN,LN, Dlng(M,NS),Dtrn(M,NS) !transverse and longitudinal dispersion coeff (for each solute and soil layer)
          enddo
      enddo
c*       
      Do n=1,NumNP
        !convert concentration from ugNO3/g soil to ugNO3/cm2 water
        Conc(n,1)=Conc(n,1)/ThNew(n)*blkdn(MatNumN(n))
        ConOld(n)=Con(n) !store old state variable
        hOld(n)=hNew(n)
        VxOld(n)=Vx(n)
        VzOld(n)=Vz(n)
        ThOld(n)=ThNew(n)
      Enddo
      Movers(2)=1
      Return
C        
C  Routine calculations 
11    Continue 
      tOld = Time
      t=Time
      dt = Step
C   For each solute do 711      
      Do 711 jjj=1,NumSol                
        xMul=1.
        alf=1.-epsi
C  Matrix initialization
        newjjj = MBand          
	  Do 13 i=1,NumNP
          DS(i)=0.    !
          Gc(i)=0.    !solute conc sink term
          Sc(i)=0.    !node area term
          B(i) =0.    !right hand side of the vector
C        GR 
          if(lOrt) B1(i)= Conc(i,jjj) !orthogonal mesh? B1=conc
          If(epsi.lt.0.001) then      !if epsi<.001, means explicit scheme
	      if(lOrt) newjjj = IADD(i) !
	       A(newjjj,i) = 0.         !if explicit= A=0
          Else
C GR
            Do 12 j=1,MBandD          !For implicit scheme
              A(j,i)=0.               !A(MBand, node)
12          Continue
          Endif
13      Continue
C
C  Recasting sinks 
C
cccz
cccz zhuangji adopt node based expression
cccz zhuangji: the previous sink recast use a concentration-based weight for "sink average"
         Do n=1,NumNP
            Gc(n)=cSink(n,jjj)  + cSink_OM(n,jjj)      !solute sink at each node and gas
            Sc(n)=NodeArea(n)         !Node area at each node
         Enddo

C
C Assembling matrixes
C  calulating local element matrix (s(j1,j2= stiffness matrix (diffusion and advection) 
C  F(i1= source/sink term vector)
C  Assembling them in to global system matrix A(iB,i1) and right hand side B(i1)
C  Apply boundary condition- update A and B
C  AE- area of the element
C  GcE average sink term
C 
        Do 21 Level=1,NLevel
          If(lUpW.eq.1) then          !if upwind scheme
            If(Level.eq.NLevel) then  !use current velocoties if final level
              Do i=1,NumNP
                VxH(i)=Vx(i)
                VzH(i)=Vz(i)
              Enddo
            Else                      !Use old velocity for intermediate level
              Do i=1,NumNP
                VxH(i)=VxOld(i)
                VzH(i)=VzOld(i)
              Enddo
            Endif
            Call Disper(jjj,NumNP,NMat,Dispxx,Dispzz,Dispxz,VxH,VzH,
     !               ThNew,hNew,Dmol,Dlng,Dtrn,NumSD,NMatD,MatNumN)
            Call WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,VxH,VzH,Dispxx, ! weighing factor for finite element assembly
     !                    Dispzz,Dispxz)
          Endif


          Do 15 i=1,NumNP
            F(i)=0.                       !residual vector F
            If(Level.eq.NLevel) DS(i)=0.  !dispersion source term?
15        Continue
C
C  Numerical stability
C
c          CourMax=1.
          Courant=0.
          dtMx(2)=1.E+30  !max alow timestep
          dt1=1.E+30
          dt2=1.E+30
          NumSEl=0        !number of sub element
          Do 19 n=1,NumEl                 !n is each element
            NUS=4                         !quadrilateral elem
            If(KX(n,3).eq.KX(n,4)) NUS=3  !triangle elem
            Do 18 k=1,NUS-2
              NumSEl=NumSEl+1
              i=KX(n,1)
              j=KX(n,k+1)
              l=KX(n,k+2)
              List(1)=i
              List(2)=j
              List(3)=l
              Ci(1)=x(l)-x(j)                                 !element width, shape functions
              Ci(2)=x(i)-x(l)
              Ci(3)=x(j)-x(i)
              Bi(1)=y(j)-y(l)                                 !element height
              Bi(2)=y(l)-y(i)
              Bi(3)=y(i)-y(j)
              AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.                 !area of each element
              If(Level.eq.NLevel) then                        !compute max allwble time step size- based on velocity and element dim
                delX=amax1(abs(Ci(1)),abs(Ci(2)),abs(Ci(3)))
                delY=amax1(abs(Bi(1)),abs(Bi(2)),abs(Bi(3)))
                VxMax=amax1(abs(Vx(i))/ThNew(i),
     !                  abs(Vx(j))/ThNew(j),
     !                  abs(Vx(l))/ThNew(l))
                VzMax=amax1(abs(Vz(i))/ThNew(i),
     !                  abs(Vz(j))/ThNew(j),
     !                  abs(Vz(l))/ThNew(l))
                CourX=VxMax*dt/delX                               !compute courant number in x and y direction
                CourY=VzMax*dt/delY
                Courant=amax1(Courant,CourX,CourY)
                If(VxMax.gt.0.) dt1=CourMax*delX/VxMax
                If(VzMax.gt.0.) dt2=CourMax*delY/VzMax
                dtMx(2)=dmin1(dtMx(2),dt1,dt2)
              Endif
              If(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.     !if cylindrical ordinates
              GcE=(Gc(i)+Gc(j)+Gc(l))/3.                          !average of the solute concentration sink term                          
              If(Level.eq.NLevel) then                            !if final level- use current values 
                Yb=hNew(i)*Bi(1)+hNew(j)*Bi(2)+hNew(l)*Bi(3)
                Yc=hNew(i)*Ci(1)+hNew(j)*Ci(2)+hNew(l)*Ci(3)
                DmE=Dmol(jjj)*(ThNew(i)*Tau(hNew(i))+             !Tsu effect of the soil matrix on diffusion
     &                          ThNew(j)*Tau(hNew(j))+
     &                          ThNew(l)*Tau(hNew(l)))/3.         !average of (theta*molecular diffusion)0each node of the element 
                ConE=(Con(i)+Con(j)+Con(l))/3.                    !avg concentration
              Else
                Yb=hOld(i)*Bi(1)+hOld(j)*Bi(2)+hOld(l)*Bi(3)      !if not final level- use old values 
                Yc=hOld(i)*Ci(1)+hOld(j)*Ci(2)+hOld(l)*Ci(3)
                DmE=Dmol(jjj)*(ThOld(i)*Tau(hOld(i))+
     &                        ThOld(j)*Tau(hOld(j))+
     &                        ThOld(l)*Tau(hOld(l)))/3.
                ConE=(ConOld(i)+ConOld(j)+ConOld(l))/3.
              EndIf
              Yc=Yc+2.*AE
              Yb2=Yb*Yb
              Yc2=Yc*Yc
              FMul=xMul*AE/4.
              SMul1=-1./AE/4.*xMul
              M=MatNumE(n)
              DisL=Dlng(M,jjj)
              DisT=Dtrn(M,jjj)
              Vabs=Sqrt(Yb2+Yc2)
              If(lUpW.eq.1) then
                NS=NumSEl
                W1=WeTab(1,NS)
                W2=WeTab(2,NS)
                W3=WeTab(3,NS)
                Wx(1)=
     &          2.*VxH(i)*(W2-W3)+VxH(j)*(W2-2.*W3)+VxH(l)*(2.*W2-W3)
                Wx(2)=
     &          VxH(i)*(2.*W3-W1)+2.*VxH(j)*(W3-W1)+VxH(l)*(W3-2.*W1)
                Wx(3)=
     &          VxH(i)*(W1-2.*W2)+VxH(j)*(2.*W1-W2)+2.*VxH(l)*(W1-W2)
                Wz(1)=
     &          2.*VzH(i)*(W2-W3)+VzH(j)*(W2-2.*W3)+VzH(l)*(2.*W2-W3)
                Wz(2)=
     &          VzH(i)*(2.*W3-W1)+2.*VzH(j)*(W3-W1)+VzH(l)*(W3-2.*W1)
                Wz(3)=
     &          VzH(i)*(W1-2.*W2)+VzH(j)*(2.*W1-W2)+2.*VzH(l)*(W1-W2)
              Endif
              Do 17 j1=1,3
                i1=List(j1)
                F(i1)=F(i1) - FMul*(GcE+Gc(i1)/3.)                !Source/sink vector
                If(Level.eq.NLevel) DS(i1)=DS(i1)-AE/3.
                Do 16 j2=1,3
                  i2=List(j2)
                  S(j1,j2)=SMul1*(Bi(j1)*Bi(j2)+Ci(j1)*Ci(j2))*DmE  !sriffness matrix?- include diffusion and dispersion contributions
                  If(Vabs.gt.0.) then
                    S(j1,j2)=S(j1,j2)+Smul1*ConE*
     &              (Bi(j1)*Bi(j2)*(DisL*Yb2+DisT*Yc2)+
     &              (Bi(j1)*Ci(j2)+Bi(j2)*Ci(j1))*(DisL-DisT)*Yb*Yc+
     &              Ci(j1)*Ci(j2)*(DisL*Yc2+DisT*Yb2))/Vabs/2./AE
                    
                    S(j1,j2)=S(j1,j2) -
     &              (-Bi(j1)*Yb-Ci(j1)*Yc)
     &              *(ConE+Con(i2)/3.)/4.*SMul1*xMul
                    
                    If(lUpW.eq.1) S(j1,j2)=S(j1,j2)-xMul*
     !                (Bi(j2)/40.*Wx(j1)+Ci(j2)/40.*Wz(j1))
                  Endif
                  If(Level.ne.NLevel) then
                    B(i1)=B(i1)-alf*S(j1,j2)*Conc(i2,jjj)         !right hand side vector
                  Else
C GR
	              if (lOrt) then
	                 call Find(i1,i2,kk,NumNP,MBandD,IAD,IADN)
	                 iB=kk
	              else
	                 iB=MBand+i2-i1
				  endif
                    A(iB,i1)=A(iB,i1)+epsi*S(j1,j2)               !uses stiffness matrix
                  Endif
16              Continue
17            Continue
18          Continue  !for each node associated with an element
19        Continue    !for each element

          Do 20 i=1,NumNP
            If(Level.ne.NLevel) then
              B(i)=B(i)-alf*F(i)                                  !intermidiate level. B includes sink terms and previous time step conc
            Else
	        if (lOrt) newjjj=IADD(i)    

C              A(MBand,i)=A(MBand,i)+DS(i)*ThNew(i)/dt
              A(newjjj,i)=A(newjjj,i)+DS(i)*ThNew(i)/dt           !update A for the final level
              B(i)=B(i)+DS(i)*ThOld(i)/dt*Conc(i,jjj)-epsi*F(i)   !update right hand side
            Endif
20        Continue
21      Continue
C
C     Boundary conditions
C
        Do 114 i=1,NumNP
          Radd=0.
          cKod=0
          If(CodeW(i).ne.0) then
            Do 111 j=1,NumBP
              If(KXB(j).eq.i) then
                If(CodeS(i).eq.1) then
                  cKod=1
                  cBnd=Conc(i,jjj)
                  Goto 112
                Endif
                If(CodeS(i).eq.3.or.CodeS(i).eq.6) then
                  cKod=1
                  cBnd=VarBS(j,jjj)
                  Goto 112
                Endif
                If(CodeS(i).eq.-3.or.CodeS(i).eq.-6) then
                  cKod=3
                  cBnd=VarBS(j,jjj)
                  Goto 112
                Endif
                If(CodeS(i).eq.4.or.CodeS(i).eq.-4) then
                  cKod=3
                  If(Q(i).gt.0.) then
                    cBnd=VarBS(j,jjj)
                    if(cbnd.ne.0.0) write (*,*) cbnd
                  Else
                    cBnd=0.
                    If(CodeW(i).eq.-4.and.VarBW(j,1).gt.0.0)
     &                Radd=VarBW(j,1)*Width(j)*VarBS(j,jjj)
                      If (Radd.ne.0.0) write(*,*) Radd
                  Endif
                  Goto 112
                Endif
                If((CodeW(i)).eq.2.or.(CodeW(i).eq.-5).or.
     &                    (CodeW(i).eq.-7)) then   !CodeW=-7 is for free drainage - DT May 29, 2012
                  cKod=2
                  Goto 112
                Endif
                If((CodeW(i)).eq.1)
     &                    then   !bottom boundary and flux - DT Dec 14, 2017
                  if (Q(i).lt.0) cKod=2 ! only consider flux out for now
                  Goto 112
                Endif
              Endif
111         Continue
C
C  Dirichlet boundary condition
C
112         If(cKod.eq.1) then
              if(lOrt) then
	          A(IADD(i),i) = 1.d30
	          B(i)=1.d30*cBnd
	        else
                Do 113 j=1,2*MBand-1
                  A(j,i)=0.
113             Continue
                A(MBand,i)=1.d0
                B(i)=cBnd
	        endif
            Endif
C
C  Cauchy boundary condition
C
            If(cKod.eq.3) then
		    B(i)=B(i) - Q(i)*cBnd - Radd

	      endif
C
C  Neuman boundary condition
C
            If(cKod.eq.2) then !This might apply also to a unit gradient condition
	        if(lOrt) newjjj = IADD(i) 
              A(newjjj,i)=A(newjjj,i) + epsi*Q(i)  !
              B(i)=B(i) - alf*Q(i)*Conc(i,jjj)
            Endif
          Endif   
114     Continue

c   pond   on the soil-atmosphere surface
cMisha 18/9 2006
C  SURFACE POND FOR TRANSPORT            
          Do 313 i=1,NumBP
          n=KXB(i)
          k=CodeS(n)
         If (iabs(k).eq.4) then	                        ! surface ponding	
			HNEWS=DMAX1(hNew(N),0.0D0)                  ! ponding water height
		    HOLDS=DMAX1(hold(N),0.0D0)                  
		     newjjj=1	
		     if (lOrt) newjjj = IADD(n)             
              A(newjjj,n)=A(newjjj,n)+Width(i)*HNEWS/dt   !update for ponded water
              B(n)=B(n)+Width(i)*Conc(n,jjj)*HOLDS/dt     !update right hand side
          END IF
 313	  CONTINUE	 
C
c   pond   on the soil-atmosphere surface
cMisha 18/9 2006
C
C

C
C  Solve the global matrix equation for transport
C  AC=B
        If(epsi.lt.0.001) then
          Do 22 i=1,NumNP
C GR
	       if (lOrt) newjjj = IADD(i)
             B(i)=B(i)/A(newjjj,i)
C            B(i)=B(i)/A(MBand,i)
22        Continue
           else if(lOrt) then
	       call ILU(A,NumNP,MBandD,IAD,IADN,IADD,A1)
             call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,
     !	            IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,0,
     !                MNorth,MaxItO,2)
 
C  Solve the system of linear equations
C
        else
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
          Do  214 i=2,NumNP
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
C
        Do i=1,NumNP
C GR
          if(lOrt) B(i) = B1(i)
          if (B(i).lt.0) then
                 iii=1
            endif

          Conc(i,jjj)=sngl(B(i))
          Conc(i,jjj)=amax1(0.0,Conc(i,jjj))
        Enddo    
711   Continue
      Do i=1,NumNP
        ConOld(i)=Con(i)
        hOld(i)=hNew(i)
        VxOld(i)=Vx(i)
        VzOld(i)=Vz(i)
        ThOld(i)=ThNew(i)
      Enddo
*
      Return
10    Call errmes(im,il)
      Return
      End
* 
      Subroutine Disper(jjj,NumNP,NMat,Dispxx,Dispzz,Dispxz,VxH,VzH,
     !           theta,h,Dmol,Dlng,Dtrn,NumSD,NMatD,MatNum)
      Dimension VxH(NumNP),VzH(NumNP),
     !          Dmol(NumSD),Dlng(NMatD,NumSD),Dtrn(NMatd,NumSD),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),MatNum(NumNP),
     !          h(NumNP),theta(NumNP)
      Do 11 i=1,NumNP
        M=MatNum(i)
        Taun=Tau(h(i))
        Vabs=sqrt(VxH(i)*VxH(i)+VzH(i)*VzH(i))
        If(Vabs.gt.0.) then
          Dispxx(i)=Dlng(M,jjj)*VxH(i)*VxH(i)/Vabs+
     !              Dtrn(M,jjj)*VzH(i)*VzH(i)/Vabs+
     !              theta(i)*Dmol(jjj)*Taun
          Dispzz(i)=Dlng(M,jjj)*VzH(i)*VzH(i)/Vabs+
     !              Dtrn(M,jjj)*VxH(i)*VxH(i)/Vabs+
     !              theta(i)*Dmol(jjj)*Taun
          Dispxz(i)=(Dlng(M,jjj)-Dtrn(M,jjj))*VxH(i)*VzH(i)/Vabs
        Else
          Dispxx(i)=theta(i)*Dmol(jjj)*Taun
          Dispzz(i)=theta(i)*Dmol(jjj)*Taun
          Dispxz(i)=0.
        Endif
11    continue
      Return
      End
*
      subroutine WeFact(NumNP,NumEl,NumElD,x,y,KX,WeTab,VxH,VzH,Dispxx,
     !                  Dispzz,Dispxz)
      Dimension x(NumNP),y(NumNP),KX(NumElD,4),VxH(NumNP),VzH(NumNP),
     !          Dispxx(NumNP),Dispzz(NumNP),Dispxz(NumNP),
     !          WeTab(3,2*NumElD),Beta(3),List(3)
      Integer e

      TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))
      NumSEl=0
      Do 13 e=1,NumEl
        NCorn=4
        If(KX(e,3).eq.KX(e,4)) NCorn=3
        Do 12 n=1,NCorn-2
          NumSEl=NumSEl+1
          M1=KX(e,1)
          M2=KX(e,n+1)
          M3=KX(e,n+2)
          A=y(M2)-y(M1)
          B=x(M2)-x(M1)
          Beta(1)=atan2(A,B)
          A=y(M3)-y(M2)
          B=x(M3)-x(M2)
          Beta(2)=atan2(A,B)
          A=y(M1)-y(M3)
          B=x(M1)-x(M3)
          Beta(3)=atan2(A,B)
          List(1)=M1
          List(2)=M2
          List(3)=M3
          Do 11 j=1,3
            k=j-1
            If(k.eq.0) k=3
            WeTab(k,NumSEl)=0.
            M1=List(j)
            jp1=j+1
            If(j.eq.3) jp1=1
            M2=List(jp1)
            Vxx=(VxH(M1)+VxH(M2))/2.
            Vzz=(VzH(M1)+VzH(M2))/2.
            If(abs(Vxx).lt.1.e-30.and.abs(Vzz).lt.1.e-30) goto 11
            BetaV=atan2(Vzz,Vxx)
            Delta=abs(BetaV-Beta(j))
            If(Delta.gt.0.314.and.abs(Delta-3.1416).gt.0.314) goto 11
            ALeng=sqrt((x(M2)-x(M1))**2+(y(M2)-y(M1))**2)
            CBeta=cos(Beta(j))
            SBeta=sin(Beta(j))
            Val=Vxx*CBeta+Vzz*SBeta
            VV=sqrt(Vxx*Vxx+Vzz*Vzz)
            DLL=(Dispxx(M1)+Dispxx(M2))/2.
            DLT=(Dispxz(M1)+Dispxz(M2))/2.
            DTT=(Dispzz(M1)+Dispzz(M2))/2.
            DAL=abs(DLL*CBeta*CBeta+2.0*CBeta*SBeta*DLT+DTT*SBeta*SBeta)
            Vel=VAL*ALeng
            Disp=2.0*DAL
            aa=11.
            If(Disp.gt.0.) aa=abs(Vel/Disp)
            If(Disp.lt.1.e-30.or.abs(Vel).lt.0.001*VV.or.abs(aa).gt.10.)
     !            then
              If(abs(Vel).lt.0.001*VV) WeTab(k,NumSEl)=0.0
              If(Vel.gt.0.001*VV) WeTab(k,NumSEl)=1.0
              If(Vel.lt.-0.001*VV) WeTab(k,NumSEl)=-1.0
            Else
              WeTab(k,NumSEl)=1.0/TanH(Vel/Disp)-Disp/Vel
            Endif
11        Continue
12      Continue
13    Continue
      Return
      End
*
      Function Tau(h)
      If(h.gt.-0.2) then
        ttau=0.67
      Else
        ttau=0.59-0.12*alog10(abs(h))
      Endif
       If(ttau.lt.0.) ttau=0.
       Tau=ttau
       Return
      End