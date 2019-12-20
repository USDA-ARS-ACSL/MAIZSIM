
* I adapted this to add a nitrogen fertilizer. In this case
* Nitrogen is added to the concentration for node i
* This module allows multiple fertilizer and water applications
* cAppl is total mg N  (1 mg/m2 = 1kg/ha) applied to surface area of grid of nodes speficied
* Solute (N) in 2dsoil is ug/cm3 or mg/l
* to get applied N in Kg/ha 1/(.01*rowsp*0.5)*applied/1000/1000*10000
* where N is applied to all the nodes in the surface and the plant
* is on the left of the domain.

* Note that the first management day needs to be at least one day later than the 
* initial time in the file in the time folder. Otherwise the applied fertilizer
* will be divided by 0 value of thNew. 
* dt - 9/10/2018 added organic fertilizer as residue (Litter)

      Subroutine Mngm()
      Include 'public.ins'
      Include 'NitVar.ins'
      Integer Max_fert_times,Fert_times,Num_nodes,nAppl,applied
      Parameter (Max_fert_times=75, Max_fert_nodes=300)
      Dimension Bi(3),Ci(3)
      Character*10 Date(Max_fert_times)
      Real AppliedAmount, AppliedAmount_NL, AppliedAmount_CL

      Common /Mngmnt/ tAppl(Max_fert_times), cAppl(Max_fert_times),
     !         cAppl_CL(Max_fert_times), cAppl_NL(Max_fert_times),
     !         Num_nodes(Max_fert_nodes), nAppl(Max_fert_times,
     !         Max_fert_nodes),Fert_times,applied,ModNum
      t=time
      If(lInput.eq.1) then
        Open(40,file=ManagementFile,err=20)
        Read(40,*,Err=20) 
        Read(40,*,Err=20)
        Read(40,*,Err=20) Fert_times
        Read(40,*,Err=20)
        Do i=1,Fert_times
	    Read(40,*,Err=20) Date(i), cAppl(i), cAppl_CL(i),  cAppl_NL(i),
     !                      Num_nodes(i)
            Read(40,*,Err=20) 
            Read(40,*,Err=20) (nAppl(i,j),j=1,Num_nodes(i))
        EndDo
        NumMod=NumMod+1
        ModNum=NumMod
        Close(40)
        Do i=1,Fert_times
           tAppl(i)=julday(date(i))
          EndDo
        tNext(ModNum)=tAppl(1)
        if (Fert_Times.EQ.0) tNext(ModNum)=1e22
      applied=1
      Endif
      If(Abs(time-tNext(ModNum)).lt.0.001*Step) then
        jj=applied
        AETotal=0
        Do in=1,Num_nodes(jj)
           Do e=1,NumEl
             NUS=4
             If(KX(e,3).eq.KX(e,4)) NUS=3
             Do k=1,NUS-2
c find node numbers of vertices of triangle
               i=KX(e,1)
               j=KX(e,k+1)
               l=KX(e,k+2)
               Ci(1)=x(l)-x(j)
               Ci(2)=x(i)-x(l)
               Ci(3)=x(j)-x(i)
               Bi(1)=y(j)-y(l)
               Bi(2)=y(l)-y(i)
               Bi(3)=y(i)-y(j)
               AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
               if (i.eq.nAppl(jj,in)) then
                  AETotal=AETotal+AE/3
                endif
                if (j.eq.nAppl(jj,in)) then
                  AETotal=AETotal+AE/3
                endif
                if (l.eq.nAppl(jj,in)) then
                  AETotal=AETotal+AE/3
                endif           
             enddo
            enddo              
          Enddo  ! Num_nodes loop
 
C convert mg N to total ug NO3 in the soil domain. 
            AppliedAmount=cAppl(jj)/AETotal
            AppliedAmount_NL=cAppl_NL(jj)/AETotal
            AppliedAmount_CL=cAppl_CL(jj)/AETotal
       Do in=1,Num_nodes(jj)
         Conc(nAppl(jj,in),1)=AppliedAmount/ ThNew(nAppl(jj,in))           
     &            *1000*62./14. +conc(nAppl(jj,in),1) 
         Nl(nAppl(jj,in))=AppliedAmount_NL*1000
         Cl(nAppl(jj,in))=AppliedAmount_CL*1000
       Enddo
c             
             
c  application at surface is mg/m2 need to find area of elements that the fertilizer was 
c  applied to and then distribute over element. Some nodes may be part of two elements. 
c  need to account for this. 
        applied=applied+1
        if (applied.gt.Fert_times) then
           tNext(ModNum)=1.E+32
         else
           tnext(ModNum)=tAppl(applied)
        endif
      Endif
      Return
20    Stop 'Mngm data error'
      End
        