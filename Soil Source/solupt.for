      subroutine SoluteUptake()
      include 'public.ins'
      IncLude 'puplant.ins'
      Include 'puweath.ins'
      real MMUpN, bi(3),ci(3)
      Common /SUP/ TotWSink,TotSSink,WincrSink,TotSSINK2
      Common /Biom/ BM
      Common /slpt/ Iflag
      If(lInput.eq.1) then
        Open(94,file='TotSink.res')
        TotWSink=0.
        TotSSINK=0.
        TotSSINK2=0.0
        Iflag=0
        MMUpN=0.0
      Endif
      if(ITIME.EQ.1) iflag=1
      SIncrSink=0.
      WIncrsink=0.
      SincrSink2=0.0
C calculate sum of N in all forms in soil using triangluar elements
	 Do n=1,NumEl
           cSink(n,1)=0.
           NUS=4
           if(KX(n,3).eq.KX(n,4)) NUS=3
*         Loop on subelements
           AE=0.0
           do k=1,NUS-2
             i=KX(n,1)
             j=KX(n,k+1)
             l=KX(n,k+2)
             Ci(1)=x(l)-x(j)
             Ci(2)=x(i)-x(l)
             Ci(3)=x(j)-x(i)
             Bi(1)=y(j)-y(l)
             Bi(2)=y(l)-y(i)
             Bi(3)=y(i)-y(j)
             AE=AE+(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
             CSink(n,1)=CSink(n,1)+(Conc(i,1)+
     &             Conc(j,1)+Conc(l,1))/3.
           Enddo
          CSink(n,1)=CSink(n,1)*Sink(n)
          SIncrSink=SIncrSink+cSink(n,1)*step*14./62.*AE

       Enddo


      do n=1,NumEl
	 NUS=4
         sum1=0.0
	 if(KX(n,3).eq.KX(n,4)) NUS=3
	 do j=1,NUS
           Sum1=Sum1+Conc(KX(n,j),1)*Sink(n)/NUS
c           IF(Conc(KX(n,j),1).GT.0.000001) THEN
c              MMUpN = 9.00/(1.+0.000000248/(1000.*Conc(KX(n,j),1)))
c              MMUpN=MMUpN*RTWT(n)
c             else
c               MMUpN=0.0
c            ENDIF
c            cSink(n,1)=cSink(n,1)+MMUpN/1000./NUS
	 enddo
         SIncrSink2=SIncrSink2+sum1*Area(n)*Step*14./62.
         WIncrSink=WIncrSink+Sink(n)*Area(n)*Step

      enddo
      TotWSink=TotWSink+WIncrSink
      TotSSINK=TotSSINK+SIncrSink
      TotSSINK2=TOTSSink2+SIncrSink2
      If(ITIME.eq.24.and.iflag.eq.1) then
        Write(94,*) Jday,TotWSink,TotSSink,TotSSInk2
        iflag=0
      Endif
      return
      end
