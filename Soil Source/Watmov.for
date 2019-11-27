*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
* codes that are positive are for prescribed pressure heads (flux is calculated). 
* Codes that are negative or 0 are for prescibed fluxes. 
* for constant or prescribed bc:
*       Dirichlet BC is prescribed or constant head
*       Neumann boundary condition is prescribed or constant flux
* CDT added a drainage boundary Nov 2007. This is like a seepage face but the nodes always have
* drainage and it is usually horizontal.
      subroutine WaterMover_New ()
      Include 'public.ins'
      Include 'puplant.ins'
      Double precision A,B,C
      Double precision dt,dtOld,t,tOld,PI,DPI,F2
      real ATG,HSP
      Logical Explic,ItCrit,FreeD
      Real  hOld_1(NumNPD)
      Dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD),
     !    Cap(NumNPD),ListE(NumElD),E(3,3),iLoc(3),Fc(NumNPD),
     !    Sc(NumNPD),B_1(NumNPD),ThOld_1(NumNPD),A_1(MBandD,NumNPD)
      Dimension Bii(3),Cii(3)
      Common /WaterM/ ThOld(NumNPD),hOld(NumNPD),hTemp(NumNPD),
     !                ConAxx(NumElD),ConAzz(NumElD),ConAxz(NumElD),
     !                MaxIt,TolTh,TolH,dt,dtOld,tOld,
     !                thR(NMatD),hSat(NMatD),thSat(NMatD),
     !                 isat(NumBPD),FreeD, CriticalH
      If (lInput.eq.0) goto 11  
       FreeD=.true.
       CriticalH=-0.01D0
       Do i=1,NumNP
          hOld(i) = hNew(i)
          hTemp(i) = hOld(i)
       Enddo
*
       Do i=1,NumEl
        ConAxz(i)=0.
        ConAxx(i)=1.
        ConAzz(i)=1.
      Enddo
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
     !                  hCritA,hCritS,dtMx(1),hTab1,hTabN,EPSI_Heat,
     !                   EPSI_Solute
        close(40) 
c      hCritS=1.0e9        !needs to be high for ponded infiltration to work
                        ! need to remove it from the parameter file
                         ! remove comment later when runoff is implemented

      call IADMake(KX,NumNP,NumEl,NumElD,MBandD,IAD,IADN,IADD)
      
      call SetMat(lInput,NumNP,hNew,hOld,NMat,MatNumN,Con,Cap,
     !                  BlkDn,hTemp,Explic,ThNew,hTab1,hTabN,
     !                  hSat,ThSat,ThR, ThAvail,ThFull,
     !                  FracOM, FracSind, FracClay,
     !                  TupperLimit, TLowerLimit, SoilFile)
c  assign bulk density
       
             
      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                   ConAxx,ConAzz,ConAxz,Vx,Vz)
c   Calculate Total Available Water in Profile
     
      
      Do i = 1,NumNP
        ThOld(i)=ThNew(i)
      Enddo
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
c
c   Start of iteration loop
c     
1111  Iter=0
      Explic=.false.
      Do i=1,NumNP
        Fc(i)=0.
        Sc(i)=0.
        ThAvail(i)=0.0
      Enddo
C
C  Recasting sinks 
C

cccz zhuangji, delete this part, for the sink is already casted
cccz based on node
c         Do n=1,NumEl
c           NUS=4
c           if(KX(n,3).eq.KX(n,4)) NUS=3
c*         Loop on subelements
c           do k=1,NUS-2
c             i=KX(n,1)
c             j=KX(n,k+1)
c             l=KX(n,k+2)
c             Ci=x(l)-x(j)
c             Cj=x(i)-x(l)
c             Cl=x(j)-x(i)
c             Bi=y(j)-y(l)
c             Bj=y(l)-y(i)
c             Bl=y(i)-y(j)
c             AE=(Cl*Bj-Cj*Bl)/2.
c             Fc(i)=Fc(i)+Sink(n)*AE/3.
c             Fc(j)=Fc(j)+Sink(n)*AE/3.
c             Fc(l)=Fc(l)+Sink(n)*AE/3.
c             Sc(i)=Sc(i)+AE/3.
c             Sc(j)=Sc(j)+AE/3.
c             Sc(l)=Sc(l)+AE/3.
c            enddo
c         enddo
c*      Do n=1,NumEl
c*        NUS=4
c*        If(KX(n,3).eq.KX(n,4)) NUS=3
c*        Do j=1,NUS
c*          Fc(KX(n,j))=Fc(KX(n,j))+Sink(n)*Area(n)/NUS
c*          Sc(KX(n,j))=Sc(KX(n,j))+Area(n)/NUS
c*        Enddo
c*      Enddo
c     Do i=1,NumNP
c        Fc(i)=Fc(i)/Sc(i)
c      Enddo

cccz directly take the sink and nodearea 
      Do n=1,NumNP
         Fc(n)=Sink(n)
         Sc(n)=NodeArea(n)
      Enddo
    
C calculate total available water in root zone.      
      ThetaFull=0.0
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
			 m=MatNumN(i)
			
			 thAWi=ThFull(MatNumN(i))
			 thAWl=ThFull(MatNumN(l))
			 thAWj=ThFull(MatNumN(j))
			 if (rtwt(i).le.1.0e-6) ThAWi=0.0
			 if (rtwt(j).le.1.0e-6) ThAWj=0.0
			 if (rtwt(l).le.1.0e-6) ThAWl=0.0
			 ThetaFull=ThetaFull+AE*(ThAWi+ThAWj+ThAWl)/3.
		   Enddo

		   
		Enddo

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
     !                  TupperLimit, TLowerLimit, SoilFile)

c
c  RESET: assembling of the matrixes
c
      xMul=1.
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
c                iCheck=1
c                else 
                 CodeW(n)=5
                 hNew(n)=0.
                endif
              Else
                If (Q(n).ge.0.) then
                  CodeW(n)=-5
                  Q(n)=0
c                  iCheck=1
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
c              if (VarBW(i,3).gt.0) then   ! evaporation case only remove comment later when runoff is implemented
                 Q(n)=-VarBW(i,3)*Width(i)
c               Endif
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
            If (hNew(n).ge.hCritS) then
              CodeW(n)=4
              hNew(n)=hCritS
            Endif
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
cccz
c         print*,time,',','Water Diff'
         call ILU (A,NumNP,MBandD,IAD,IADN,IADD,A1)
         call OrthoMin(A,B1,B,NumNP,MBandD,NumNPD,IAD,IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,0,
     !                MNorth,MaxItO)
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
c       if (time.gt.38064.055) then
c          viii=0
c         endif
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
      
c      CDT adjust for runoff here    
      Do i=1,NumBP
          n=KXB(i)
          If (CodeW(n).eq.-4) then
               if (hnew(n).ge.CriticalH) then
                 RO(n)=(hNew(n)-(CriticalH))/dt*Width(i)
                 hNew(n)=CriticalH
                 hOld(n)=hNew(n)
                endif
             endif
          Enddo
c
c   Calculation of velocities
c
cdt - moved this here              
      call Veloc(NumNP,NumEl,NumElD,hNew,x,y,KX,ListNE,Con,
     !                  ConAxx,ConAzz,ConAxz,Vx,Vz)
      call SetMat(lInput,NumNP,hNew,hOld,NMat,MatNumN,Con,Cap,
     !               BlkDn, hTemp,Explic,ThNew,hTab1,hTabN,
     !               hSat,ThSat,ThR, ThAvail,ThFull,
     !               FracOM, FracSind, FracClay,
     !               TupperLimit, TLowerLimit, SoilFile)



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
      
cdt - calculate available water content in root zone (cm3 per slab)

      ThetaAvail=0.0
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
			 ThetaAvail=ThetaAvail+AE*(Thi+Thj+Thl)/3.
		   Enddo
	   
		Enddo


cdt - calculate actual boundary fluxes to see what we have
      Do 1299 n=1,NumNP
        If (CodeW(n).eq.-4) then
        QN=B_1(n)+DS(n)+F(n)*(ThNew(n)-ThOld_1(n))/dt 
           do 1199 j=1,IADN(n)
              QN=QN+A_1(j,n)*hNew(IAD(j,n))
1199        continue
         QAct(n)=QN
        End if
  
1299   Continue

c calculate runoff to look at mass balances. RIght now we don't 
c   route ponded water to runoff
c disabled for now (4/12/2011) - need to finish the code       
c      Do i=1, NumBP
c        if (codew(i).eq.-4) then
c           n=KXB(i)
c           RO(n)=0.0
c            HSP=0.03 !EMPIRICAL PARAMETER, HSP~=dz/3
c			PI=3.141592653589793238D0
c			DPI=1.0d0/PI
c			ATG=aTAN(hNew(N)/HSP)+PI/2.d0 !Continuous Heaviside step function
cc Delta is the derivative of the step function x hnew			
c			Delta=HSP/(HSP*HSP+hNew(n)*hNew(n))*hNew(n) 
cc F2 is a weighting function			
c			F2=(ATG+delta)*DPI
c			F2=Dmin1(F2,1.0D0)
c			
c		   IF(HNEW(N).LT.-0.01) F2=1.0D-10
c          if((hNew(n).gt.(-0.01)).and.(F2.gt.0.98)) then
c            HNEWS=DMAX1(hNew(N),0.0D0)
c		      HOLDS=DMAX1(hold_1(N),0.0D0)
c		      RO(n)= Q(n)-QAct(n) !Width(i)*(-(HNEWS-HOLDS))/step
c              hNew(n)=-0.01
c              hOld(n)=hNew(n)
c        endif
c          endif
c       Enddo
       


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

