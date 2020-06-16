
* I adapted this to add a nitrogen fertilizer. In this case
* Nitrogen is added to the concentration for node i
* This module allows multiple fertilizer and water applications
* Amount_Applied is total mg N  (1 mg/m2 = 1kg/ha) applied to surface area of grid of nodes speficied
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
      Include 'nitvar.ins'
      Integer Max_fert_times,Fert_times,Num_nodes,node_Applied,applied
      Parameter (Max_fert_times=75, Max_fert_nodes=300)
      Dimension Bi(3),Ci(3)
      Character*10 Date(Max_fert_times)
      Real AppliedAmount, AppliedAmount_NL, AppliedAmount_CL, MaxY,
     !    y_limit, fert_depth

      Common /Mngmnt/ tAppl(Max_fert_times),
     !         Amount_Applied(Max_fert_times),
     !         Amount_Applied_CL(Max_fert_times), 
     !        Amount_Applied_NL(Max_fert_times),
     !         Num_nodes(Max_fert_nodes), node_Applied(Max_fert_times,
     !         Max_fert_nodes),Fert_times,applied,ModNum, 
     !         fert_depth(Max_fert_times)
      t=time
      If(lInput.eq.1) then
        Open(40,file=ManagementFile,err=20)
        Read(40,*,Err=20) 
        Read(40,*,Err=20)
        Read(40,*,Err=20) Fert_times
        Read(40,*,Err=20)
        Do i=1,Fert_times
            Read(40,*,Err=20) Date(i), Amount_Applied(i), fert_depth(i),
     &     Amount_Applied_CL(i), Amount_Applied_NL(i)
                         
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
        
cdt modified code to find the nodes in an applicatino depth rather than input the nodes
     
c fert depth range must be calclated from Y at surface and given depth       
        ! find max y - surface node
        maxY=0
        do i=1, NumNp
          if (y(i).gt.maxY) maxY=y(i)  
        end do
        y_limit=maxY-fert_depth(jj)
c now assign node numbers to the list of applied nodes where the depth
c  is equal to or less than the fertilizer applied depth

          ii=0 ! counter for Num_nodes
          do i=1, NumNP
          if(y(i).ge.y_limit) then
             ii=ii+1
             node_applied(jj,ii)=i
           End If
          End Do
          
          Num_nodes(jj)=ii
        
        
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
               if (i.eq.node_Applied(jj,in)) then
                  AETotal=AETotal+AE/3
                endif
                if (j.eq.node_Applied(jj,in)) then
                  AETotal=AETotal+AE/3
                endif
                if (l.eq.node_Applied(jj,in)) then
                  AETotal=AETotal+AE/3
                endif           
             enddo
            enddo              
          Enddo  ! Num_nodes loop
 
C convert mg N to total ug NO3 in the soil domain. 
            AppliedAmount=Amount_Applied(jj)/AETotal
            AppliedAmount_NL=Amount_Applied_NL(jj)/AETotal
            AppliedAmount_CL=Amount_Applied_CL(jj)/AETotal
       Do in=1,Num_nodes(jj)
         Conc(node_Applied(jj,in),1)=
     &            AppliedAmount/ ThNew(node_Applied(jj,in))           
     &            *1000*62./14. +conc(node_Applied(jj,in),1) 
         Nl(node_Applied(jj,in))=AppliedAmount_NL*1000
         Cl(node_Applied(jj,in))=AppliedAmount_CL*1000
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
        