c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      subroutine GasMover ()
      include 'public.ins'

      double precision A,B,C,P,Sum
      integer gKod

      dimension A(MBandD,NumNPD),B(NumNPD),F(NumNPD),DS(NumNPD),
     !                S(3,3),Ci(3),Bi(3),List(3),
     !                Diff(NumNPD),Ac(NumNPD),Gc(NumNPD),Sc(NumNPD)
	 
      Common /GasM/ ThOld(NumNPD),ThANew(NumNPD),ThAOld(NumNPD),
     !                ThTot(NMatD),ThATr(NMatD),bTort(NMatD),
     !                ThAMin(NMatD),Dair(NumGD),
     !                NLevel,epsi
      
      
      C     -------------------------( SOLOXY )-------------------------------
C     THIS SUBROUTINE CALCULATES THE OXYGEN CONCENTRATION IN THE
C     ROOT ZONE.
C     ------------------------------------------------------------------
      INCLUDE 'COMMON96.FOR'
      DIMENSION PORS(20,10),DPR(0:21,0:11),
     &           OXLOSS(20,10),H(20,0:11),V(0:21,10)
     &           ,alp(0:20),bet(0:20)
       DTSEC=86400./24.
       OWID=1./WIDTH/WIDTH
       ODEP=1./DEPTH/DEPTH
       OWID2=OWID/2.
       ODEP2=ODEP/2.
       ODTSEC=1./DTSEC
       EITH = 8./3.
       DO L=1,NL
         DO K=1,NKH
           PORS(L,K)=AMAX1(TOPORS(L,K)-VH2OC(L,K),0.08*TOPORS(L,K))
           D11=0.66*(PORS(L,K)-0.1)*0.178*(TS(L,K)/273.+1.)**2
           DPR(L,K)=AMAX1(D11,0.)
           OXLOSS(L,K)=RESPS(L,K)*10.*8.31*(273.+TS(L,K))*ODTSEC
     &               /44./DEPTH/WIDTH
         ENDDO
       ENDDO
       DO L=1,NL
         H(L,0)=DPR(L,1)*OWID
         H(L,NKH)=DPR(L,NKH)*OWID
         DO K=1,NKH-1
           H(L,K)=(DPR(L,K)+DPR(L,K+1))*OWID2
         ENDDO
C Forth
         alp(1)=1.
         bet(1)=0.
         Do K=1,NKH
           A=H(L,K-1)
           B=H(L,K)
                 E=PORS(L,K)*ODTSEC
           C=A+B+E
           F=CXT(L,K)*E
           D=1./(C-Alp(K)*A)
           alp(K+1)=B*D
           bet(K+1)=(F+A*Bet(K))*D
         Enddo
         CXT(L,NKH)=Bet(NKH+1)/(1.-Alp(NKH+1))
C Back
         Do K=NKH,2,-1
            CXT(L,K-1)=Bet(K)+Alp(K)*CXT(L,K)
         Enddo
      ENDDO
      DO K=1,NKH
C 8/3 = 1./(3/4)/(1/2)
        V(0,K)=DPR(1,K)*ODEP*EITH
        V(NL,K)=DPR(NL,K)*ODEP
        DO L=1,NL-1
          V(L,K)=(DPR(L,K)+DPR(L+1,K))*ODEP2
        ENDDO
C Forth
        alp(1)=0.
        bet(1)=0.21
        Do L=1,NL
          A=V(L-1,K)
          B=V(L,K)
                E=PORS(L,K)*ODTSEC
          C=A+B+E
          F=CXT(L,K)*E-OXLOSS(L,K)
          D=1./(C-Alp(L)*A)
          Alp(L+1)=B*D
          Bet(L+1)=(F+A*Bet(L))*D
        Enddo
        CXT(NL,K)=Bet(NL+1)/(1.-Alp(NL+1))
C Back
        Do L=NL,2,-1
          CXT(L-1,K)=Bet(L)+Alp(L)*CXT(L,K)
        Enddo
      ENDDO
       DO L=1,NL
       DO K=1,NKH
         RESPS(L,K)=0.
         IF(CXT(L,K).LT.0.) CXT(L,K)=0.
       ENDDO
       ENDDO
      RETURN
      END