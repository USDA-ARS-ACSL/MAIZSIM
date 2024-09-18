*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
* codes that are positive are for prescribed pressure heads (flux is calculated). 
* Codes that are negative or 0 are for prescibed fluxes. 
* for constant or prescribed bc:
*       Dirichlet BC is prescribed or constant head
*       Neumann boundary condition is prescribed or constant flux
* CDT added a drainage boundary Nov 2007. This is like a seepage face but the nodes always have
* drainage and it is usually horizontal.
      subroutine WaterMover ()
      Include 'public.ins'
      Include 'puplant.ins'
      Include 'puweath.ins'
      include 'PuSurface.ins'
      
      Double precision A,B,C, B_1, A_1
      Double precision dt,dtOld,t,tOld,PI,DPI,F2
      real ATG,HSP
cccz move it to "PuSurface.ins" for public use 
cccz  Double precision CriticalH, CriticalH_R
      Logical Explic,ItCrit,FreeD
      Real  hOld_1(NumNPD)
      Real Dif(NumNPD)
      Integer trigger_Runoff, p_Runoff
      Dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD),
     !    Cap(NumNPD),ListE(NumElD),E(3,3),iLoc(3),Fc(NumNPD),
     !    Sc(NumNPD),B_1(NumNPD),ThOld_1(NumNPD),A_1(MBandD,NumNPD)
      Dimension Bii(3),Cii(3)
      Common /WaterM/ ThOld(NumNPD),hTemp(NumNPD),hOld(NumNPD),
     !                ConAxx(NumElD),ConAzz(NumElD),ConAxz(NumElD),
     !                MaxIt,TolTh,TolH,dt,dtOld,tOld,
     !                thR(NMatD),hSat(NMatD),
     !                isat(NumBPD),FreeD
      If (lInput.eq.0) goto 11  
        FreeD=.true.
        CriticalH=5.1D0
c       CriticalH_R=0.01D0
  
        hOld(:) = hNew(:)
        hTemp(:) = hOld(:)
        RO(:)=0.0
*
        ConAxz(:)=0.
        ConAxx(:)=1.
        ConAzz(:)=1.
      
      Explic=.false.
*
      im=50
      il=0
      Open(40,file=WaterFile, status='old',ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) 
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) MaxIt,TolTh,TolH,
     !                  hCritA,CriticalH,dtMx(1),hTab1,hTabN,EPSI_Heat,
     !                   EPSI_Solute
        close(40) 

        
      call IADMake(KX,NumNP,NumEl,NumElD,MBandD,IAD,IADN,IADD)
c  assign bulk density      
      call SetMat(lInput,NumNP,hNew,hOld,NMat,MatNumN,Con,Cap,
     !                  BlkDn,hTemp,Explic,ThNew,hTab1,hTabN,
     !                  hSat,ThSat,ThR, ThAvail,ThFull,
     !                  FracOM, FracSind, FracClay,
     !                  TupperLimit, TLowerLimit,soilair,
     !                  SoilFile,ThAMin,ThATr)
c  need to initialize arrays after assigning h values if water ws input
     
          hOld(:) = hNew(:)
          hTemp(:) = hOld(:)
      
             
      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                   ConAxx,ConAzz,ConAxz,Vx,Vz)
c   Calculate Total Available Water in Profile
     
      
      ThOld(:)=ThNew(:)

      dt=Step
      Movers(1)=1
      Return
C
C   Routine calculations
C       

11    continue 
      tOld = Time
      t=Time
      dt=Step
      
cccz set the auto irrgation part before the iteration
      do k=1, NumBp
        i=KXB(k)
        if(abs(CodeW(i)).eq.4) then
           Q(i)=Q(i)+Qautoirrig(i)
           if (Q(i).gt.0.0) CodeW(i)=-4  !cccz make sure bc changes if Qn goes > 0 (infiltration) after adding the autoirrigation
        endif
      enddo
      

     
c
c   Start of iteration loop
c     
1111  Iter=0
      Explic=.false.

        Fc(:)=0.
        Sc(:)=0.
        ThAvail(:)=0.0
cccz Let us initialize the RO (runoff array) 
        RO(:)=0.0
cccz

C
C  Recasting sinks 
C

cccz directly take the sink and nodearea 
         Fc(:)=Sink(:)
         Sc(:)=NodeArea(:)
    

C
C  Start of an iteration
C                
12    Continue

*
 

c             If (Iter.le.1) then
C
C  SetMat: hydraulic properties for every node based on
C  new values of pressure head
C
      call SetMat(lInput,NumNP,hNew,hOld,NMat,MatNumN,Con,Cap,
     !                  BlkDn, hTemp,Explic,ThNew,hTab1,hTabN,
     !                  hSat,ThSat,ThR, ThAvail,ThFull,
     !                  FracOM, FracSind, FracClay,
     !                  TupperLimit, TLowerLimit,soilair,
     !                  SoilFile,ThAMin,ThATr)

c
c  RESET: assembling of the matrixes
c
      xMul=1.

      !B(:)=0.
      !F(:)=0.
      Do 212 i=1,NumNP
 
        B(i)=0.
        F(i)=0.
C added 1 line
        if(lOrt) B1(i)=hNew(i) 
        If (Iter.eq.0) DS(i)=0.
        Do 211 j=1,MBandD
           A(j,i)=0.
211     Continue
212   Continue
C
C     Loop on elements
C
      Do 216 n=1,NumEL
        CondI=ConAxx(n)
        CondJ=ConAzz(n)
        CondK=ConAxz(n)
        NUS=4
        If (KX(n,3).eq.KX(n,4)) NUS=3
        Do 215 k=1,NUS-2
          i=KX(N,1)
          j=KX(N,k+1)
          l=KX(N,k+2)
          iLoc(1)=1
          iLoc(2)=k+1
          iLoc(3)=k+2
          Ci=x(l)-x(j)
          Cj=x(i)-x(l)
          Ck=x(j)-x(i)
          Bi=y(j)-y(l)
          Bj=y(l)-y(i)
          Bk=y(i)-y(j)
          AE=(Ck*Bj-Cj*Bk)/2.
          CapE=(CAP(i)+CAP(j)+CAP(l))/3.0
          ConE=(Con(i)+Con(j)+Con(l))/3.0
          If (KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          AMul=xMul*ConE/4./AE
          BMul=xMul*ConE/2.
          FMul=xMul*AE/12.
          If (Iter.eq.0) then
            SinkE=Fc(i)+Fc(j)+Fc(l)
            DS(i)=DS(i)+Fmul*(SinkE+Fc(i))
            DS(j)=DS(j)+Fmul*(SinkE+Fc(j))
            DS(l)=DS(l)+Fmul*(SinkE+Fc(l))
          Endif
          F(i)=F(i)+FMul*4.
          F(j)=F(j)+FMul*4.
          F(l)=F(l)+FMul*4.
          If (KAT.ge.1) then
            B(i)=B(i)+BMul*(CondK*Bi+CondJ*Ci)
            B(j)=B(j)+BMul*(CondK*Bj+CondJ*Cj)
            B(l)=B(l)+BMul*(CondK*Bk+CondJ*Ck)
          Endif
          E(1,1)=CondI*Bi*Bi+2.*CondK*Bi*Ci+CondJ*Ci*Ci
          E(1,2)=CondI*Bi*Bj+CondK*(Bi*Cj+Ci*Bj)+CondJ*Ci*Cj
          E(1,3)=CondI*Bi*Bk+CondK*(Bi*Ck+Ci*Bk)+CondJ*Ci*Ck
          E(2,1)=E(1,2)
          E(2,2)=CondI*Bj*Bj+2.*CondK*Bj*Cj+CondJ*Cj*Cj
          E(2,3)=CondI*Bj*Bk+CondK*(Bj*Ck+Bk*Cj)+CondJ*Cj*Ck
          E(3,1)=E(1,3)
          E(3,2)=E(2,3)
          E(3,3)=CondI*Bk*Bk+2.*CondK*Bk*Ck+CondJ*Ck*Ck
          Do 214 i=1,3
            iG=KX(n,iLoc(i))
            Do 213 j=1,3
              jG=KX(n,iLoc(j))
              if(lOrt) then
                call Find(iG,jG,kk,NumNP,MBandD,IAD,IADN)
                A(kk,iG)=A(kk,iG)+AMul*E(i,j)
              else
                iB=iG-jG+1
                if(iB.ge.1) A(iB,jG)=A(iB,jG)+AMul*E(i,j)
                end if
213         Continue
214       Continue
215     Continue
216   Continue
c
c     Determine boundary fluxes
c
        B_1=B       !dt save B matrix for flux calcs later
        A_1=A 
      Do 220 n=1,NumNP
                   
        If (CodeW(n).lt.1) goto 220
        QN=B(n)+DS(n)+F(n)*(ThNew(n)-ThOld(n))/dt
        if(lOrt) then
          do 1117 j=1,IADN(n)
            QN=QN+A(j,n)*hNew(IAD(j,n))
1117        continue
        else
        QN=QN+A(1,n)*hNew(n)
        Do 219 j=2,MBand
          k=n-j+1
          If (k.ge.1) then
217          QN=QN+A(j,k)*hNew(k)
          End if
          k=n+j-1
          If (k.le.NumNP) then
            QN=QN+A(j,n)*hNew(k)
          End if
219     Continue
        End if
        Q(n)=QN
220   Continue


c
c     Complete construction of RHS vector and form effective matrix
c
      Do 221 i=1,NumNP
        j=1
        if(lOrt) j=IADD(i)  
        A(j,i)=A(j,i)+F(i)*Cap(i)/dt
        B(i)=F(i)*Cap(i)*hNew(i)/dt-F(i)*(ThNew(i)-ThOld(i))/dt+
     !     Q(i)-B(i)-DS(i)

221   Continue
c
c     Modify conditions on seepage faces
c
      If (NSeep.ne.0) then
        Do 312 i=1,NSeep
          iCheck=0
          NS=NSP(i)
          Do 311 j=1,NS
            n=NP(i,j)
            If (CodeW(n).eq.-2) then
              If (hNew(n).lt.0.) then
                iCheck=1
              Else
                CodeW(n)=2
                hNew(n)=0.
              Endif
            Else
              If (iCheck.gt.0.or.Q(n).ge.0.) then
                CodeW(n)=-2
                Q(n)=0.
                iCheck=1
              End if
            Endif
311       Continue
312     Continue
      Endif
c
c     Modify conditions on Drainage boundaries
c
      If (NDrain.ne.0) then
        Do 3120 i=1,NDrain
          NDS=NDR(i)
c          iCheck=0
          Do 3110 j=1,NDS
            n=ND(i,j)
            If (CodeW(n).eq.-5) then
              if (hNew(n).ge.0) then
                 CodeW(n)=5
                 hNew(n)=0.
                endif
              Else
                If (Q(n).ge.0.) then
                  CodeW(n)=-5
                  Q(n)=0
              Endif
            Endif
3110       Continue
3120     Continue
      Endif

*     Free Drainage    
       if(FreeD) then      
        do i=1,NumBP
          n=KXB(i)
          k=CodeW(n)
          if(k.eq.-7) Q(n)=-Width(i)*Con(n)
         End do
      end if  
c
c     Modify potential surface flux boundaries
c
      If (NSurf.ne.0) then
        Do 313 i=1,NumBP
          n=KXB(i)
          k=CodeW(n)
          If (Explic.and.iabs(k).eq.4) then
             CodeW(n)=-iabs(k)
             Goto 313
          Endif


c
c   Critical surface pressure on the soil-atmosphere surface
c   valid for evaporation only
c
         If (K.eq.4) then
            If (abs(Q(n)).gt.abs(-VarBW(i,3)*Width(i))
     &                          .or.Q(n)*(-VarBW(i,3)).le.0) then
              CodeW(n)=-4
                 Q(n)=-VarBW(i,3)*Width(i)
            Endif

          Goto 3131
         Endif
c
c   Surface flux on on the soil-atmosphere surface
c
          If (K.eq.-4) then
            If (hNew(n).le.hCritA) then
              CodeW(n)=4
              hNew(n)=hCritA
               Goto 3131
            Endif
c            If (hNew(n).ge.hCritS) then
c              CodeW(n)=4
c              hNew(n)=hCritS
c           Endif
         Endif
3131     continue
c
c    pond  on  the soil-atmosphere surface
cMisha 18/9 2006
cMK----------------------------------------------------------------------------------     
 


		if ((CodeW(n).eq.-4).and.(q(n).gt.0)) then
c Ponded infiltration measurement is from Misha Kouznetzov
			HSP=0.009D0 !EMPIRICAL PARAMETER, HSP~=dz/3 - was 0.03
			PI=3.141592653589793238D0
			DPI=1.0d0/PI
c            ATG=0.0D0
			ATG=aTAN((hNew(N))/HSP)+PI/2.0d0 !Continuous Heaviside step function
c Delta is the derivative of the step function x hnew			
			Delta=HSP/(HSP*HSP+hNew(n)*hNew(n))*hNew(n) 
c F2 is a weighting function			
			F2=(ATG+delta)*DPI
			F2=Dmin1(F2,1.0D0)
			
		   IF(HNEW(N).LT.CriticalH) F2=1.0D-10
         		HNEWS=DMAX1(hNew(N),0.0D0)
  	      	HOLDS=DMAX1(hold(N),0.0D0)	
c update right and left sides of equation for flux due to the change in the ponded
c head (if any) 
	          j=1
			   if(lOrt) j=IADD(n)              
			    A(j,n)=A(j,n)+Width(i)*F2/dt
			    B(n)=B(n)+Width(i)*(F2*hNew(n)-(HNEWS-HOLDS))/dt          
			      
		 endif   !codew= -4

 
c   pond   on the soil-atmosphere surface
cMisha 18/9 2006
cMK-----------------     
CMK  	          
313       Continue
	Endif  !NSurf <> 0


c
c===== Dirich: constant head boundaries. Prscribed pressure heads
c===== are incorporated directly into matrices
c
      Do 412 n=1,NumNP
        If (CodeW(n).lt.1) goto 412
C added 4 lines
        if(lOrt) then
              A(IADD(n),n)=10.d30
                B(n)=10.d30*hNew(n)
         else
          Do 411 m=2,MBand
          k=n-m+1
          If (k.gt.0) then
            B(k)=B(k)-A(m,k)*hNew(n)
            A(m,k)=0.
          Endif
           l=n+m-1
          If (NumNP-l.ge.0) then
            B(l)=B(l)-A(m,n)*hNew(n)
          Endif
            A(m,n)=0.
411      Continue
          A(1,n)=1.
          B(n)=hNew(n)
C added 1 line
           end if
412      Continue

cdt this is where the old solver sits
c
c  Solving of the system of equations   
c

      if(lOrt) then
        call ILU (A,NumNP,MBandD,IAD,IADN,IADD,A1)
        call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,0,
     !                MNorth,MaxItO,1)
	endif
      if (.not.lOrt) then 
c*   Reduction
      Do 513 n=1,NumNP
        Do 512 m=2,MBand
        If (abs(A(m,n)).lt.1.e-30) goto 512
          C=A(m,n)/A(1,n)
          i=n+m-1
          If (i.gt.NumNP) goto 512
          j=0
          Do 511 k=m,MBand
            j=j+1
            A(j,i)=A(j,i)-C*A(k,n)
511       Continue
          A(m,n)=C
          B(i)=B(i)-A(m,n)*B(n)
512     Continue
        B(n)=B(n)/A(1,n)
513   Continue
c*   Back substitution
      n=NumNP
514   Do 515 k=2,MBand
        l=n+k-1
        If (l.gt.NumNP) goto 516
        B(n)=B(n)-A(k,n)*B(l)
515   Continue
516   n=n-1
      If (n.gt.0) goto 514
C
C   End of an iteration
C 
      EndIf 
cdt  conventional solver ends here


      Do 613 i=1,NumNP
          hTemp(i)=hNew(i)
            if(lOrt) B(i)=B1(i)
          hNew(i) =sngl(B(i))
 613   Continue
      Iter =Iter+1
      If (Explic) goto 619

      hMax=0.0
      do i=1,NumNP
         hAbs=Abs(hNew(i))
         if (hAbs.gt.hMax) hMax = hAbs
      enddo
C
C    Test for convergence
C
       ItCrit=.true.
       do 615 i=1,NumNp
            m=MatNumN(i)
            EpsTh=0.
            EpsH=0.
            hlev=hSat(m)-20.0
            if (hTemp(i).lt.hLev.and.hNew(i).lt.hLev) then
               Th=ThNew(i)+cap(i)*(hNew(i)-hTemp(i))
     !            /(ThSat(m)-ThR(m))
               EpsTh=abs(ThNew(i)-Th)
             else
               EpsH=abs(hNew(i)-hTemp(i))
c               Dif(i)=EpsH
             endif
             if (EpsTh.gt.TolTh.or.EpsH.gt.TolH)then
                   ItCrit=.false.
                   goto 616
             endif
615     continue
616   Continue
      If (.not.ItCrit) then
        If (Iter.lt.MaxIt.AND.hMax.lt.10.*abs(hCritA)) then
c adjust for runoff within an iteration
 

          Goto 12
        Elseif (dt.le.dtMin) then
          Explic=.true.
          Do 617 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
617       Continue
          Goto 12
        Else
C
C   Save new boundary conditions If any
C
          Do i=1,NumNP
            If (CodeW(i).gt.0) hOld(i)=hNew(i)
          Enddo

          Do 618 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
618       Continue
          dt=dmax1(dt/3.0D0,dtMin)
          dtOpt=dt
          t=tOld+dt
          goto 1111
        Endif
      Endif
c
c  end of iteration loops
c
619   Continue
      Do 20 i=1,NumNP
        If (CodeW(i).eq.99) then
          Q(i)=0.
          CodeW(i)=0
        Endif
        If (hNew(i).eq.hCritA) hNew(i)=0.999*hCritA
20    Continue
      Time = Time+dt-Step
      tOld=Time
      dtOld=dt
      Step=dt
      




c
c   Calculation of velocities
c
          
           
                         
      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                  ConAxx,ConAzz,ConAxz,Vx,Vz)
      call SetMat(lInput,NumNP,hNew,hOld,NMat,MatNumN,Con,Cap,
     !               BlkDn, hTemp,Explic,ThNew,hTab1,hTabN,
     !               hSat,ThSat,ThR, ThAvail,ThFull,
     !               FracOM, FracSind, FracClay,
     !               TupperLimit, TLowerLimit,soilair,
     !               SoilFile,ThAMin,ThATr)



c     
c     Final assignments
c


      Do i=1,NumNP
         hOld_1(i)=hOld(i)  !save h old and th old for runoff calculations
         ThOld_1(i)=ThOld(i)
         hOld(i) =hNew(i)  
         ThOld(i)=ThNew(i)
cdt it seems this if statement should be in the first part of this do block
c  hNew will always be the same as hOld?         
        If (CodeW(i).lt.1) then
          hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
          hNew(i) =hTemp(i)
        Else
          hTemp(i)=hNew(i)
        Endif
      Enddo
      
cdt - calculate actual boundary fluxes to see what we have
      Do 1299 n=1,NumNP
        If ((CodeW(n).eq.-4).or.(CodeW(n).eq.1)) then
        QN=B_1(n)+DS(n)+F(n)*(ThNew(n)-ThOld_1(n))/dt 
           do 1199 j=1,IADN(n)
              QN=QN+A_1(j,n)*hNew(IAD(j,n))
1199        continue
         QAct(n)=QN
        End if
  
1299   Continue

      
cdt - calculate available water content in root zone

      ThetaAvailRZ=0.0
      Do n=1,NumEl
		   NUS=4
		   if(KX(n,3).eq.KX(n,4)) NUS=3
		   Sum1=0.
		   Sum2=0.
c*         Loop on subelements
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
			 Thi=0.0
			 Thj=0.0
			 thl=0.0
			 if (rtwt(i).ge.1e-6) Thi=ThAvail(i)
			 if (rtwt(j).ge.1e-6) Thj=ThAvail(j)
			 if (rtwt(l).ge.1e-6) Thl=Thavail(l)
			 ThetaAvailRZ=ThetaAvailRZ+AE*(Thi+Thj+Thl)/3.
		   Enddo

		   
		Enddo
      



cccz start to calculate the runoff
cccz this is the water source part, i.e., the exfiltration from soil surface
c only calculate this when the surface nodes are atmospheric boundary nodes
      do k=1, NumBp
        i=KXB(k)
        if (((hnew(i).ge.CriticalH)).and.(abs(codeW(i)).eq.4)) then
          RO(i)=max(Q(i)-Qact(i),0.0D0)
          hNew(i)=CriticalH+h_Pond(k)         ! cccz could be CriticalH_R, but we force it to 
          hOld(i)=hNew(i)
        endif
       Enddo         
cccz turn this on for Ex_4 plastic mulching
cccz #ifdef EX_4P
cccz                 if (k.ge.2.and.k.le.7) then
cccz                    RO(k)=max(Q(k),0.0D0)
cccz                    hOld(k)=hNew(k)
cccz                endif
cccz #endif
cccz turn this on for Ex_4 plastic mulching

       
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Return
10    Call errmes(im,il)
      Return
      End
c*
      subroutine Veloc(NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,ConAxx,
     !                 ConAzz,ConAxz,Vx,Vz)
      Dimension hNew(NumNP),x(NumNP),y(NumNP),ListNE(NumNP),Con(NumNP),
     !          KX(NumElD,4),Vx(NumNP),Vz(NumNP),ConAxx(NumEl),
     !          ConAzz(NumEl),ConAxz(NumEl),List(3)
      Integer e
      Do 11 i=1,NumNP
        Vx(i)=0.
        Vz(i)=0.
11    Continue    
      Do 14 e=1,NumEl
        CAxx=ConAxx(e)
        CAzz=ConAzz(e)
        CAxz=ConAxz(e)
        NCorn=4
        If(KX(e,3).eq.KX(e,4)) NCorn=3
        Do 13 n=1,NCorn-2
          i=KX(e,1)
          j=KX(e,n+1)
          k=KX(e,n+2)
          List(1)=i
          List(2)=j
          List(3)=k
          vi=y(j)-y(k)
          vj=y(k)-y(i)
          vk=y(i)-y(j)
          wi=x(k)-x(j)
          wj=x(i)-x(k)
          wk=x(j)-x(i)
          Area=.5*(wk*vj-wj*vk)
          A=1./Area/2.
          Ai=CAxx*vi+CAxz*wi
          Aj=CAxx*vj+CAxz*wj
          Ak=CAxx*vk+CAxz*wk  
          Vxx=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k))+CAxz
          Ai=CAxz*vi+CAzz*wi
          Aj=CAxz*vj+CAzz*wj
          Ak=CAxz*vk+CAzz*wk
          Vzz=A*(Ai*hNew(i)+Aj*hNew(j)+Ak*hNew(k))+CAzz
          Do 12 m=1,3
            l=List(m)
            Vx(l)=Vx(l)-Con(l)*Vxx
            Vz(l)=Vz(l)-Con(l)*Vzz
12        Continue
13      Continue
14    Continue
      Do 15 i=1,NumNP
        Vx(i)=Vx(i)/ListNE(i)
        Vz(i)=Vz(i)/ListNE(i)
15    Continue
      Return
      End

