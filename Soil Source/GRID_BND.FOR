*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
CDT Nov, 2007 added a drainage boundary. It is input similar to seepage. The difference is
C*     that drainage always takes place from a drainage boundary - it is usually horizontal
C      seepage faces are usually used in ditches.
      Subroutine Get_Grid_and_Boundary()
      Include 'public.ins'
      Include 'nitvar.ins'
*
      Integer e
      Logical lConst
      Dimension idummy(4)
*
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
        CodeT(n)=idummy(3)
        CodeG(n)=idummy(4)
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
      Enddo
C
C  Area of elements and bandwidth MBand
C
      lOrt=.false.
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
      MBand=MBand+1
      if(MBand.gt.MBandD.or.(2*MBand-1.gt.MBandD)) lOrt=.true.
	
      If(.not.lConst) IJ=NumNP
CDT note SWMS_2D used 200   (lChem.and. -> this went before 2 and after .or)
      if(MBand.gt.10.or.NumNP.gt.200) lOrt=.true.
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
        If(abs(CodeG(i)).eq.3) NVarBG=NvarBG+1
      Enddo

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

        Do j=1,NumNP
          il=il+1
          Read(41,*,ERR=10)i,Nh(i),Ch(i),NL(i),CL(i),Nm(i),Cm(i),NNH4(i)
     &      ,(Conc(i,jjj),jjj=1,NumSol),Tmpr(i),hNew(i)
     
          Ch_old(i)=Ch(i)
          Nh_old(i)=Nh(i)
          CL_old(i)=CL(i)
          NL_old(i)=NL(i)
          Cm_old(i)=Cm(i)
          Nm_old(i)=Nm(i)
          NNO3_old(i)=Amax1(0.,Conc(i,1)*ThNew(i)*14./62.)
          NNH4_old(i) =  NNH4(i)
          ThOld(i)=ThNew(i)
          TmprOld(i)=Tmpr(i)
          Denit(i)=0.
          Denit_old(i)=0.
        Enddo
      close(41)
! READING FROM Nodal_N.dat alternative file.  BY GR EDN      
      Write(*,*) 'Reading Grid Done'
      Return
10    Call errmes(im,il)
      End

