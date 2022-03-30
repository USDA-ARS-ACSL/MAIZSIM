*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Subroutine Rootgrow()
      Include 'public.ins'
      Include 'puplant.ins'
      Dimension iLoc(4)
      Integer n,n1,e,SW10,SW20,SW30,SW40
      Real*4  WeightNeighbor(4),LeftDist, RightDist, TopDist,BottomDist
      
                             
      t=sngl(Time)
      If (lInput.eq.0) goto 11
      HourlyCarboUsed=0          !hourly initializtion in crop model
C
C  Reading of the input files and initial calculations 
C
c      open(105, file='root.dat',recl=150)
c      write(105,'(10A12)')'time    ','elem   ' , 'AWR    ', 'PDWR     ', 
c     !  'pcrts    ' ,'pcrs    ','convr   ', 'rootwt    ', 'id',
c     !   'Matnum'
!        TotalRootWeight=0.0    !temporary
          
11    continue
C

C     ---------- Taken form R U T G R O   of  G L Y C I M
C
C  CALCULATE SOIL MECHANICAL RESISTANCE TO ROOT GROWTH IN EACH CELL
C
C  CALCULATE GROWTH BY NEW ROOTS IF ROOTS GREW AT POTENTIAL
C  RATE DURING PAST PERIOD. START WITH SOIL CELL WHERE CONDITIONS
C  ARE MOST FAVORABLE FOR GROWTH


C  No Plant = No Root Activity
      If(NShoot.eq.0.or.isGerminated.eq.0) Return
C



C   Taken from ACTGRO of GLYCIM
C   THIS PART GROWS ROOTS IN THE MOST FAVORABLE CELLS UNTIL 
C   CARBON AVAILABLE IS USED
C
C   When the cells are small, the model does not always use all the carbon available
C   We can loop this until all carbon is used up or increase the growth rate. The latter
C   would require some iteration, it may be best to loop several times, the only problem is
C   that it would weight the first cells maybe too hightly
      
      xMul=1.0D0
  
      PCRTS=0.0D0
      SW40=0
      used=0.0D0
      
      DO 60 n = 1,NumNP
           n1=iFavRoot(n)
        IF (PDWR(n1).LE.0.0D0) THEN
          ADWR(n1) = 0.0
        ELSE
          PCRTS = PCRTS + (PDWR(n1)*CONVR)
        
          IF (PCRS.GE.PCRTS) THEN
            ADWR(n1) = PDWR(n1)
            used=PCRTS
            ll=1
c           write(105,25) time, e1, ADWR(e1),PDWR(e1)*CONVR, PCRTS, 
c     !                        PCRS,Used,rtwt(e1), ll, matnumE(e1)
          ELSE IF ((SW40.GT.1).OR.(SW40.LT.1)) THEN
            SW40 = 1
            ADWR(n1)=((PDWR(n1)*CONVR)-PCRTS+PCRS)/CONVR
            ll=2

            used=used+ADWR(n1)*CONVR
c            write(105,25) time, e1, ADWR(e1),PDWR(e1)*CONVR, PCRTS, 
c     !                        PCRS,USED,rtwt(e1), ll, matnumE(e1)
          ELSE
            ADWR(n1) = 0.0
          END IF
        END IF

cccz  CALCULATE ROOT GROWTH RATE
        AWR(n1) = ADWR(n1)*Step
        
 60   CONTINUE
      fdif=pcrs-used
c      if (fdif.gt.0) then
c         iii=1;
c        endif
C
CDT 03/13/2011 Added this line from GLYCIM. PCRS returns the actual amount of carbon used 
C If the roots don't use all the carbon allocated then pcrts will be less than pcrs
C if the roots use all the carbon then the code will overshoot pcrts slightly (the growth in the element 
C  is adjusted but not the summation of pcrts.
      PCRS=AMIN1(PCRTS,PCRS)
c    cdt HourlyCarboUsed is the summation of carbohydrate actually used to  grow roots
c       summed between calls to the plant model. It should be initialized to 0 hourly
c     in the plant model 
      HourlyCarboUsed=HourlyCarboUsed+PCRS*step
      
      
cccz zhuangji: since we use the nodebased representation,
cccz   we can directly calculate the root weight based on the
cccz   node RMassY and RMassM results
      TotalRootWeight=0.0D0
      do n=1,NumNP
          TotalRootWeight=TotalRootWeight
     &      +nodeArea(n)*(RMassM(n)+RMassY(n))
      enddo

c   Actual growth

25    format (F12.6,1x,i6,2x,6G12.3,1x,2i5)
      Return
      End
*
      Function f1(PSI)
c
c  Root growth correction by soil matric potential
c
      bet=0.0028
      psi_thr=-500.0D0
      psi_ths=-150.0D0
      If (PSI.GT.psi_ths) then
          aux=1.0D0
      elseif (PSI.LT.psi_thr) then
          aux=0.0D0
      else
          pi=4.D0*DATAN(1.D0)
          aux=0.5D0*sin(pi/(psi_ths-psi_thr)
     !      *(PSI-0.5D0*(psi_ths+psi_thr)))+0.5D0
      endif
      f1=min(max(aux,0.0D0),1.0D0)
      Return
      End
*
CDT took out functions for f1 and f2 that adjusted  root diffusion coefficients based
C  on soil properties. This is all done here via carbon allocation
C AD took out subroutine SORT, and kept it in carbon_partitioning
cccz add f1 f2 and f3 back for root diffusion coef correction
      
      Function f2(TMP)
c
c  Root growth correction by soil temperature
c
      q=1.0D0
      p=10000.0D0
      u=18000.0D0
      t0=295.0D0
      t=TMP+273.0D0
      tlim=290.0D0
      if (t.gt.t0) then
          aux=max((1+exp(q-u/t0))*exp(p/t-p/t0)/(1+exp(q-u/t)),1.0D0)
      elseif (t.lt.tlim) then
          t=t0+tlim-t
          aux=max((1+exp(q-u/t0))*exp(p/t-p/t0)/(1+exp(q-u/t)),1.0D0)
      else
          aux=1.0D0
      endif
      f2=max(min(aux,1.0D0),0.0D0)
      Return
      End
*
      Function f3(RTD)
c
c  Root growth correction by root density
c
      Ymax=200.0D0
      aux=AMAX1(1.-RTD/Ymax,0.0D0)
      f3=aux
      Return
      End
