      Subroutine SoilNitrogen()
      Include 'Public.ins'
      Include 'NITVAR.ins'
      common /nitrog/ModNum
      
      Logical enough
       t=time
       If(lInput.eq.1) then
        im=400
        il=0
        Open(40,File=NitrogenFile,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
         il=il+1
        Read(40,*,ERR=10)
        im=im+1
        Read(40,*) RowSp
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        Read(40,*,ERR=10)
        Do i=1,NMat
           il=il+1
           Read(40,*,ERR=10) m,kh0(m),kL0(m),km0(m),kn0(m),kd0(m),
     &      fe(m),fh(m),r0(m),rL(m),rm(m),fa(m),nq(m),cs(m)
        Enddo         

C      initialize variables to zero here
        Do i=1, NumNP
          Ch_Old(i)=0
          Nh_Old(i)=0
          CL_Old(i)=0
          NL_Old(i)=0
          Cm_Old(i)=0
          Nm_Old(i)=0
          Denit_Old(i)=0
          NNH4_Old(i)=0
          NNO3_Old(i)=0
        EndDo
        
        Close(40)
        NumSol=1
        Call SetAbio(ew,et,ed,0,0.,0.)
        dtmx(4)=1./24
       Else

        Do i=1,NumNP
        m=MatNumN(i)
        BCh=Ch_Old(i)
        BNh=Nh_Old(i)
        BCL=CL_Old(i)
        BNL=NL_Old(i)
        BCm=Cm_Old(i)
        BNm=Nm_Old(i)
C                 ug cm3 * cm3 H2O cm-3 soil  * N / nitrate = ug N cm-3 volume
C                   multiply by 1000 to get mg N cm3 volume
        BNO3=Amax1(0.,Conc(i,1)*ThNew(i)*14./62.0/1000.0 )
        NNO3_OLD(i)=BNO3
        BNNH4=NNH4_Old(i) !8/28/2014 DT was BNH4 original typo
        BDENIT=Denit_old(i)
        DO IT=1,3
C  Rate constants
          Call SetAbio(ew,et,ed,
     &       m,(ThNew(i)+ThOld(i))/2.,(Tmpr(i)+TmprOld(i))/2.)
C  Dehumification rate (day-1)
          kh=kh0(m)*ew*et
C  Plant residue decomposition rate (day-1)
          kL=kL0(m)*ew*et
C  Organic fertilizer decomposition rate (day-1)
          km=km0(m)*ew*et 
C  Nitrification rate (day-1)
          kn=kn0(m)*ew*et
C  Denitrification rate (day-1)
C calculate mean of current concentration and concentration from past time step
          Aux=(BNO3+NNO3_Old(i))/2.
          kd=kd0(m)*ed*et*Aux/(Aux+cs(m))
C  Carbon and nitrogen fluxes
          P1 =    kh * (Ch_old(i)+BCh)/2. ! rate is given as mg C per cm3 of area (or mg per liter?)
          Q1 =    kh * (Nh_old(i)+BNh)/2.
          P2 =    kL * ((CL_old(i)+BCL)/2.) * fe(m) * fh(m)
          Q2 =    P2 / r0(m)
          P12 =   km * ((Cm_old(i)+BCm)/2.) * fe(m) * fh(m)
          Q12 =   P12 / r0(m)
          P3 =    kL * (CL_old(i)+BCL)/2.
          Q3 =    kL * (NL_old(i)+BNL)/2.
          P13 =   km * (Cm_old(i)+BCm)/2.
          Q13 =   km * (Nm_old(i)+BNm)/2.
C Added mineral nitrogen
          Added = (Q1+Q3+Q13)*Step
C Potentially immobilized and lost mineral nitrogen
          Q45pot =   kL * ((Cl_old(i)+BCL)/2.) * fe(m) / r0(m)
          Q1415pot =   km * ((Cm_old(i)+BCm)/2.) * fe(m) / r0(m)
          Q7 =  kd
          PotLost = (Q45pot + Q1415pot + Q7)*Step
C Present mineral nitrogen
CDT this uses N from previous time step
          Present  = NNO3_old(i)+NNH4_old(i)
C 'Enough' is true when it is enough mineral nitrogen for immobilization
          If(PotLost.Lt.Added+Present) then
             Enough=.true.
          else
             Enough=.false.
          endif
          If(Enough) then 
            Q45act  = Q45pot
            Q1415act= Q1415pot
          else
            Q45act  =fa(m)*Present*Q45pot  /(Q45pot+Q1415pot)
            Q1415act=fa(m)*Present*Q1415pot/(Q45pot+Q1415pot)
          endif
C
C Preferential immobilization of ammonium
C
C Nitrification
          Aux1=(NNH4_Old(i)+BNNH4)/2.  !8/28/2014 DT was BNH4 original typo
          Aux2=(BNO3+NNO3_Old(i))/2.
          Q6 =  kn * AMAX1(0.,Aux1 - Aux2/nq(m))
C Ammonium available for immobilization
          Avail = NNH4_old(i) + (Q1+Q3+Q13-Q6)*Step
C Enough is .true. if available ammonium cover all needs for immobilization
C dt 7-9-2007 Added code below to account for case when enough is false because there is almost no N in ths soil
C  then Q45act and Q1415act are zero. Have to avoid a divide by zero error.

          If((Q45act+Q1415act)*Step.LT.Avail) then
            Enough=.true.
          else
            Enough=.false.
          endif
          If(Enough) then
            Q4 = Q45act
            Q5 = 0.
            Q14= Q1415act
            Q15= 0.
          else
            Q4 = Avail    * Q45act/(Q45act+Q1415act)
             if (isnan(Q4)) Q4=0
            Q14= Avail    - Q4
            Q5=  Q45act   - Q4
            Q15= Q1415act - Q14
          endif
          P45   =  Q45act   * r0(m)
          P1415 =  Q1415act * r0(m)
C Incrementing nitrogen and carbon contents in compartments
          BCh = Ch_Old(i)+ Step * ( - P1 + P2 + P12)
          BNh = Nh_Old(i)+ Step * ( - Q1 + Q2 + Q12)
          BCL = CL_Old(i)+ Step * ( - P2 - P3 + P45)
          BNL = NL_Old(i)+ Step * ( - Q2 - Q3 + Q4 + Q5)
          BCm = Cm_Old(i)+ Step * ( - P12 - P13 + P1415)
          BNm = Nm_Old(i)+ Step * ( - Q12 - Q13 + Q14 + Q15)
          BNNH4 = NNH4_old(i) + Step *
     &             ( Q1 + Q3 + Q13 - Q4 - Q14 - Q6 )
          BNO3 = NNO3_old(i) + Step * 
     &             ( Q6 - Q7 - Q5 - Q15)
C No negative values
          BCh = AMAX1(BCh,0.)
          BNh = AMAX1(BNh,0.)
          BCl = AMAX1(BCl,0.)
          BNl = AMAX1(BNl,0.)
          BCm = AMAX1(BCm,0.)
          BNm = AMAX1(BNm,0.)
          BNNH4 = AMAX1(BNNH4,0.)
          BNO3  = AMAX1(BNO3,0.)
          BDENIT=Denit_old(i)+Q7*Step
          ENDDO
C End of iterations
          Ch(i) = BCh
          Nh(i) = BNh
          CL(i) = BCL
          NL(i) = BNL
          Cm(i) = BCm
          Nm(i) = BNm
          Denit(i)=BDENIT
          NNH4(i) = BNNH4
          NNO3_sol = BNO3
          Conc(i,1)=NNO3_sol*62./ThNew(i)/14.0*1000.0 
          TotNitO=Nh_old(i)+Nl_old(i)+Nm_Old(i)+NNO3_Old(i)+NNH4_old(i)
          TotNit=Nh(i)+Nl(i)+Nm(i)+NNO3_sol+NNH4(i)  
          DTot=TotNit-TotNitO
          Nh_old(i)=Nh(i)
          Ch_old(i)=Ch(i)
          NL_old(i)=NL(i)
          CL_old(i)=CL(i)
          Nm_old(i)=Nm(i)
          Cm_old(i)=Cm(i)
          Denit_old(i)=Denit(i)
          NNH4_old(i)=NNH4(i)
          NNO3_old(i)=amax1(0.,Conc(i,1)*ThNew(i)*14./62./1000.)
CDT the following is OK if we sum BNO3 above          
cdt          NNO3_old(i)=NNO3_sol
          ThOld(i)=ThNew(i)
          TmprOld(i)=Tmpr(i)
        Enddo
      Endif


      call Nitrogen_Mass_Balance()


      Return
10    Call Errmes(im,il)
      End
C
      Subroutine SetAbio(ew,et,ed,m,theta,T)
      Include 'public.ins'
      Common /Abio/ ThSat(NMatD),ThW(NMatD),
     &              dThH,dThL,es,Th_m,tb,QT,dThD,Th_D
      If(lInput.eq.1) then
        im=420
        il=0
        Open(40,File=BiologyFile,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
C
        Read(40,*,err=10)NM
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        Do i=1,NMat
           il=il+1
           Read(40,*,ERR=10) l,ThSat(l),ThW(l)
        Enddo
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) dThH,dThL,es,Th_m
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) tb,QT
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) dThD,Th_d
        Close(40)
        Return
      else
        ThH=ThSat(m)-dThH
        ThL=ThW(m)  +dThL
        if(theta.gt.ThH) then
          ew=es+(1.-es)*((ThSat(m)-Theta)/(ThSat(m)-ThH))**Th_m
        elseif(theta.lt.ThL) then
          if(Theta.le.ThW(m)) then
            ew=0.
          else
            ew=((Theta-ThW(m))/(ThL-ThW(m)))**Th_m
          endif
        else
          ew=1.
        endif
        eT=QT**((T-tb)/10.)
        ThD=ThSat(m)-dThD
        if(Theta.GT.ThD) then
          ed=((Theta-ThD)/(ThSat(m)-ThD))**Th_d
        else
          ed=0.
        endif
      endif
      return
10    Call Errmes(im,il)
      End

