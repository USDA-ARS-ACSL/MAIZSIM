*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
CDT Nov, 2007 added a drainage boundary. It is input similar to seepage. The difference is
C*     that drainage always takes place from a drainage boundary - it is usually horizontal
C      seepage faces are usually used in ditches.
      Subroutine Get_Grid_and_Boundary()
      Include 'public.ins'
      Include 'nitvar.ins'
      Include 'PuSurface.ins'
*
      Integer e
      Logical lConst
      Dimension idummy(4)
      
*
cccz some parameter (constant) in the Pubsurface.ins      
      CO2mass_2_Cmass=12.0107D0/44.01D0
      NO3mass_2_Nmass=14.0067D0/62.0049D0
      NH4mass_2_Nmass=14.0067D0/18.039D0
      
csb: following are for converting the ppm to ugco2/cm3 air. This value needs to be added for all the gases      
      ugGasCm3air_to_ppm(1)=509.6                            !See below text: convert the atm CO2 in [ppm]  to [ug co2/cm3 air] 
      ugGasCm3air_to_ppm(2)=700.6                            !See below text: convert the atm O2 in [ppm]  to [ug O2/cm3 air]      
      ugGasCm3air_to_ppm(3)=508.6                            !See below text: convert the atm N2O in [ppm]  to [ug N2O/cm3 air]

!Unit conversion of [ppm] to [micro gram co2/ cm3 of air]
!1 ppm [parts per million on a volume basis]= 1 umol CO2/mol of air
!At STP gas occupies:  1mol/22.4 Liter = 1 mol/22400 cm3
!1 ppm = 1 umol CO2/22400 cm3 =  4.46*10-5 umol CO2/cm3 air ? 4.56*10-5 = 1/22400  
!1umol CO2= 44.009 ug co2
!1ppm=4.46*10-5*44.009  ug CO2/cm3 air
!1ppm= 1.962*10-3 ug CO2/cm3 air
!1ppm=(1/509.6) ug CO2/cm3 air    
      
!For oxygen 
!1umol O2= 32 ug o2
!1ppm=4.46*10-5*32  ug O2/cm3 air
!1ppm=(1/700.6) ug O2/cm3 air 
      
!For N2O 
!1umol N2O= 44.08 ug N2O
!1ppm=4.46*10-5*44.08  ug N2O/cm3 air
!1ppm=(1/508.6) ug N2O/cm3 air 
csb: Gas transfer coefficient: surface gas flux change per unit gas content in the soil air at the soil surface (cm/day)

      
      im=20
      il=0
*
      Open(40,file=GeometryFile,status='old',ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) KAT,NumNP,NumEl,NumBP,IJ, NMat
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
C
C Reading nodal information
C
      Do i=1,NumNP
        il=il+1
        Read(40,*,ERR=10) n,x(n),y(n),MatNumN(n)
      Enddo
*  Reading element information
*
      Num=0
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      Do e=1,NumEL
        il=il+1
        Read(40,*,ERR=10) Num,(KX(Num,i),i=1,4),MatNuME(Num)
cdt 10/21/2002 added this line to allow input of 3 node elements
        if(KX(Num,4).eq.0) KX(Num,4)=KX(Num,3)
      Enddo
C
C  Reading boundary information
C
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      Do i=1,NumBP
        il=il+1
        Read(40,*,ERR=10) n,(idummy(j),j=1,4),Width(i)
	KXB(i)=n
        CodeW(n)=idummy(1)
        CodeS(n)=idummy(2)
        CodeT(n)=idummy(3)      !CodeT=1 means constant temperature
        CodeG(n)=idummy(4)      !CodeG=1 means constant gas concentration
        Do j=1,4
          If(idummy(j).ne.0) Movers(j)=1
        Enddo
      Enddo
C
C  Reading seepage face information
C
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) NSeep
      If (NSeep.ne.0) then
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) (NSP(i),i=1,NSeep)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        Do i=1,NSeep
          il=il+1
          Read(40,*,ERR=10) (NP(i,j),j=1,NSP(i))
        Enddo
      EndIf
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10)
      im=im+1
      il=il+1
      Read(40,*,ERR=10) NDrain
      If (NDrain.ne.0) then
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        il=il+1
        Read(40,*,ERR=10) (NDR(i),i=1,NDrain)
        im=im+1
        il=il+1
        Read(40,*,ERR=10)
        im=im+1
        Do i=1,NDrain
          il=il+1
          Read(40,*,ERR=10) (ND(i,j),j=1,NDR(i))
        Enddo
      EndIf
      Close(40)
C
C
C       
      Do n=1,NumNP
        ListN(n)=1
        ListNE(n)=0
        nodeArea(n)=0.0D0
      Enddo
C
C  Area of elements and bandwidth MBand
C
      lOrt=.true.
      lConst=.true.
      MBand=1
      xMul=1.
      Do e=1,NumEl
        Area(e)=0.
        NUS=4
        If(KX(e,3).eq.KX(e,4)) NUS=3
        Do k=1,NUS-2
          MB=1
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          CJ=x(i)-x(l)
          CK=x(j)-x(i)
          BJ=y(l)-y(i)
          BK=y(i)-y(j)
          AE=(CK*BJ-CJ*BK)/2.
          If(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          Area(e)=Area(e)+xMul*AE
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(l)=ListNE(l)+1
          If(abs(i-j).gt.MB) MB=abs(i-j)
          If(abs(i-l).gt.MB) MB=abs(i-l)
          If(abs(j-l).gt.MB) MB=abs(j-l)
          If(MB.gt.MBand) MBand=MB
          If(e.eq.1.and.k.eq.1) then
            MB1=MB
          Else
            If(MB1.ne.MB) lConst=.false.
          Endif
        Enddo
      Enddo
cccz do the computaiton of nodeArea here, moved from the carbon_partitioning      
      Do e=1,NumEl
          If(KX(e,3).eq.KX(e,4)) then
              i=KX(e,1)
              j=KX(e,k+1)
              l=KX(e,k+2)
              nodeArea(i)=nodeArea(i)+Area(e)/3.0D0
              nodeArea(j)=nodeArea(j)+Area(e)/3.0D0
              nodeArea(l)=nodeArea(l)+Area(e)/3.0D0
          Else
              i=KX(e,1)
              j=KX(e,2)
              l=KX(e,3)
              k=KX(e,4)
              nodeArea(i)=nodeArea(i)+Area(e)/4.0D0
              nodeArea(j)=nodeArea(j)+Area(e)/4.0D0
              nodeArea(l)=nodeArea(l)+Area(e)/4.0D0        
              nodeArea(k)=nodeArea(k)+Area(e)/4.0D0 
          Endif
      Enddo

      
      MBand=MBand+1
      	
      If(.not.lConst) IJ=NumNP
C
C  Numbers of boundary nodes where boundary conditions are not constant
C
      Do j=1,NumBP
        i=KXB(j)
        k=abs(CodeW(i))
        If(k.eq.3.or.k.eq.6) NVarBW=NVarBW+1
        If(k.eq.4) NSurf = NSurf +1
        k=abs(CodeS(i))
        If(k.eq.3.or.k.eq.6) NVarBS=NVarBS+1
        If(abs(CodeT(i)).eq.3) NVarBT=NvarBT+1
        If(abs(CodeG(i)).eq.3) NVarBG=NvarBG+1    !CodeG=3 when time dependent gas flux
        NVarBG=0                                  ! For now since we take the boundary CO2 from the weather file itself
      Enddo

C calculate total width of surface and bottom boundary nodes
C for scaling fluxes 
C at this point we assume the that bottom and top boundaries have the same width
      gridWidth=0.0
      do i=1, nsurf
         gridWidth=gridWidth+width(i)
      end do 
      
      
        NumSol=1

        im=450
        il=0
      Open(41,File=NodeGeomFile,ERR=10)
        im=im+1
        il=il+1
        Read(41,*)
        im=im+1
        Read(41,*)
        im=im+1
        
cccz ------------------ The Carbon/Nitrogen Component ---------------------------------
c
c     Nh(NumNP)   [ug/cm3 soil]     ELEMENTAL N in soil humus
c     Nl(NumNP)   [ug/cm3 soil]     ELEMENTAL N in soil litter
c     Nm(NumNP)   [ug/cm3 soil]     ELEMENTAL N in soil organic fertilizer
c     NH4(NumNP) [ug/cm3 soil]     Mineral N pool in the form of ammonium
c     NO3(NumNP) [ug/cm3 soil]     Mineral N pool in the form of nitrate
c
c     Cl(NumNP)   [ug/cm3]     ELEMENTAL C in the litter component of soil organic matter    
c     Cm(NumNP)   [ug/cm3]     ELEMENTAL C in manure component of soil organic matter  
c     Ch(NumNP)   [ug/cm3]     ELEMENTAL C in humus component of soil organic matter  
cccz ----------------------------------------------------------------------------------

        Do j=1,NumNP
          il=il+1
          Read(41,*,ERR=10)i,Nh(i),Ch(i),NL(i),CL(i),Nm(i),Cm(i),NH4(i)
     &      ,(Conc(i,jjj),jjj=1,NumSol),Tmpr(i),hNew(i),
     &      (g(i,jjj),jjj=1,NumG), RTWT(i)                         ! g(i,jjj):Read the initial CO2, O2 and N2O concentration in the soil

C     right now conc is input as ug NO3 per g soil. further conversions will be done in the solute mover initialization routine when we have water contents.

C   ammonia is input as ppm or ug /g soil 

cccz --------------------------------------------------------------------------------------------------------------------
cccz the unit of Conc is initially ug/g soil of NO3. Units of NO3_old(i) is "ug N as NO3 per cm3 volume", 
c     multiply by BD to get ug/ cm3 volume. then by molar mass to get N
          NO3_old(i)=Conc(i,1)*NO3mass_2_Nmass*Blkdn(MatNumN(i))
cdt       NH4 is input as ug/g soil (ppm). it is converted to ug /cm3 of soil for calculations          
cccz the unit of NH4(i) and NH4_OLD(i) should be "ug N per cm3 volume", conversion needed, i.e.,
          NH4(i)=NH4(i)*NH4mass_2_Nmass*Blkdn(MatNumN(i))
          NH4_old(i)=NH4(i)
cccz --------------------------------------------------------------------------------------------------------------------
          TmprOld(i)=Tmpr(i)
          Denit(i)=0.0D0
          Denit_old(i)=0.0D0
          NO3_from_residue(i)=0.0D0
        Enddo
      close(41)
      
!Unit conversion of initial gas concentrations in the nodes  
      
      Do j=1,NumNP
          Do jjj=1,NumG
                  g(j,jjj)= g(j,jjj)/ugGasCm3air_to_ppm(jjj) !jjj=1,2,3: CO2,O2, N2O:  gas convert the initial concentration of CO2 or O2 or N2O in [ppm]  to [ug co2 or O2 /cm3 air]               
          end do
      end do

  
      
      
! READING FROM Nodal_N.dat alternative file.  BY GR EDN      
      Write(*,*) 'Reading Grid Done'
      Return
10    Call errmes(im,il)
      End

