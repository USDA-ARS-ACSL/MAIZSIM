      Subroutine SetTDB()
      Include 'public.ins'
      Dimension Dummy(NumSurfDatD,3)
      Common /std/im1,im2,im3,im4,il1,il2,il3,il4
      If(NVarBW.ne.0) then
        If(lInput.eq.1) then
          im1=300
          il1=0
          Open(6,file=WaterBoundaryFile,status='old',ERR=200)        
          im1=im1+1
          il1=il1+1
          Read(6,*,ERR=100)
          im1=im1+1
          il1=il1+1
          Read(6,*,ERR=100)
        Endif
        If(abs(Time-tTDB(1)).lt.0.001*Step.or.lInput.eq.1) then
          im1=im1+1
          il1=il1+1
          Read(6,*,ERR=100) tTDB(1)
          im1=im1+1
          Do 11 j=1,NVarBW
            il1=il1+1
            Read(6,*,ERR=100) n,VarB
            Do 12 i=1,NumBP
              k=KXB(i)
              If(k.eq.n) then
                If(CodeW(n).eq.3.or.CodeW(n).eq.6) hNew(n)=VarB
                If(CodeW(n).eq.-3) Q(n)=Width(i)*VarB
                If(CodeW(n).eq.-6) Q(n)=VarB
                Goto 11
              Endif
12          Continue
11        Continue
        Endif
      Endif
      If(NVarBS.ne.0) then
        If(lInput.eq.1) then
          im2=325
          il2=0
          Open(7,file='VarBS.dat',status='old',ERR=200)
          im2=im2+1
          il2=il2+1
          Read(7,*,ERR=200)
          im2=im2+1
          il2=il2+1
          Read(7,*,ERR=200)
        Endif
        If(abs(Time-tTDB(2)).lt.0.001*Step.or.lInput.eq.1) then
          im2=im2+1
          il2=il2+1
          Read(7,*,ERR=200) tTDB(2)
          im2=im2+1
          Do 21 j=1,NVarBS
            il2=il2+1
            Read(7,*,ERR=200) n,(Dummy(jjj,1),jjj=1,NumSol)
            Do 22 i=1,NumBP
              k=KXB(i)
              If(k.eq.n) then
                Do jjj=1,NumSol
                  VarBS(i,jjj)=Dummy(jjj,1)
                Enddo
                If (CodeS(n).eq.3) then
                  Do jjj=1,NumSol
                    Conc(n,jjj)=Dummy(jjj,1)
                  Enddo
                Endif
                Goto 21
              Endif
22          Continue
21        Continue
        Endif
      Endif
      If(NVarBT.ne.0) then
        If(lInput.eq.1) then
          im3=350
          il3=0
          Open(8,file='VarBT.dat',status='old',ERR=300)
          im3=im3+1
          il3=il3+1
          Read(8,*,ERR=300)
          im3=im3+1
          il3=il3+1
          Read(8,*,ERR=300)
        Endif
        If(abs(Time-tTDB(3)).lt.0.001*Step.or.lInput.eq.1) then
          im3=im3+1
          il3=il3+1
          Read(8,*,ERR=300) tTDB(3)
          im3=im3+1
          Do 31 j=1,NVarBT
            il3=il3+1
            Read(8,*,ERR=300) n,VarB1,VarB2,VarB3
            Do 32 i=1,NumBP
              k=KXB(i)
              If(k.eq.n) then
                VarBT(i,1)=VarB1
                VarBT(i,2)=VarB2
                VarBT(i,3)=VarB3
                If(CodeT(n).eq.3.or.CodeT(n).eq.6) Tmpr(n)=VarB1
                Goto 31
              Endif
32          Continue
31        Continue
        Endif
      Endif
      If(NVarBG.ne.0) then
        If(lInput.eq.1) then
          im4=375
          il4=0
          Open(9,file='VarBG.dat',status='old',ERR=400)
          im4=im4+1
          il4=il4+1
          Read(9,*,ERR=400)
          im4=im4+1
          il4=il4+1
          Read(9,*,ERR=400)
        Endif
        If(abs(Time-tTDB(4)).lt.0.001*Step.or.lInput.eq.1) then
          im4=im4+1
          il4=il4+1
          Read(9,*,ERR=400) tTDB(4)
          im4=im4+1
          Do 41 j=1,NVarBG
            il4=il4+1
            Read(9,*,ERR=400) n,((Dummy(jjj,m),m=1,3),jjj=1,NumG)
            Do 42 i=1,NumBP
              k=KXB(i)
              If(k.eq.n) then
                Do jjj=1,NumG
                  Do m=1,3
                    VarBG(i,jjj,m)=Dummy(jjj,m)
                  Enddo
                Enddo
                If (CodeG(n).eq.3) then
                  Do jjj=1,NumG
                    g(n,jjj)=Dummy(jjj,1)
                  Enddo
                Endif
                Goto 41
              Endif
42          Continue
41        Continue
        Endif
      Endif
      Return
100   Call errmes(im1,il1)
200   Call errmes(im2,il2)
300   Call errmes(im3,il3)
400   Call errmes(im4,il4)
      Return
      End

