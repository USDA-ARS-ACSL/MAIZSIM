* This is for tillage processes
* The input for the tillage processes will be read and processed on the day of tillage
* Switch for tillage application: tillageApplied [0: no till, 1: tillage]

		  Subroutine TillageProcess()
		  Include 'public.ins'
		  include 'nitvar.ins'

		  Double precision tillage_date
		  Real till_Depth, y_tillageLimit,tillArea
		  Real vwcBefore,vwcAfter,ThAvg,ThTill
		  Real concBefore,concAfter,concAvg
		  Real NhBefore,NhAfter,NhAvg
		  Real NLBefore,NLAfter,NLAvg
		  Real NmBefore,NmAfter,NmAvg
		  Real ChBefore,ChAfter,ChAvg
		  Real CLBefore,CLAfter,CLAvg
		  Real CmBefore,CmAfter,CmAvg

		  Real NH4Before,NH4After, NH4Avg
		  Real DenitBefore,DenitAfter,DenitAvg
            
		  Real TmprBefore,TmprAfter, TmprAvg
		  Real HeadBefore,HeadAfter, HeadAvg
            

		  Integer tillageApplied,node_tillApplied
		  Integer Num_tillNodes
		  character*10 till_Date
		  Character*132 InString
		  Parameter(Max_till_nodes=300)

		  Common/Till/tillage_date, till_Depth,tillageApplied,
     !        ModNum,node_tillApplied(Max_till_nodes),
     !        ThTill(Max_till_nodes)

		  !---Start of tillage input reading
		  If (lInput.eq.1) then
			  Open(40,file=ManagementFile,err=20)
			  !For reading the tillage input
16              Read (40,'(A132)') InString
			  if (InString(1:14).ne.'[Tillage]') goto 16

			  read(40,*,Err=20)                                           ! Read tillage information
			  read(40,*,Err=20) tillageApplied                            ! 1: tillage applied, 0: No-till
			  NumMod=NumMod+1
			  ModNum=NumMod

			  if (tillageApplied.gt.0) then                               ! If tillage is applied
				  read(40,*,Err=20)
				  read(40,*,Err=20) till_Date, till_Depth                 ! Only one date for tillage
				  tillage_date=julday(till_Date)                          ! Tillage_date,till_Depth,tillageApplied [Common/Till]
				  tNext(ModNum)=tillage_date
				  close(40)
			  else
				  tNext(ModNum)=1e22
				  close(40)
			  end if
		  end if
		  !---End of tillage input reading


		  !---The following is to be carried out on the day of tillage
            If(Abs(time-tNext(ModNum)).lt.0.001*Step) then

            !Get the maximum y depth of the domain
                maxY=0
                do i=1, NumNp
                    if (y(i).gt.maxY) maxY=y(i)
                end do

                !Based on the soil domain, get the max Y from surface to apply tillage
                y_tillageLimit=maxY-till_Depth

                !Get the nodes which should be considered for the tillage
                ii=0                                                  !Counter for Num_nodes
                do i=1, NumNP
                    if(y(i).ge.y_tillageLimit) then
                        ii=ii+1
                        node_tillApplied(ii)=i                        !These are the nodes tillage effects, i is the node number, ii is the count
                    End If
                End Do

                Num_tillNodes=ii                                      !Total number of nodes that tillage effects
                !selecting nodes-end

                !calculate total areas of these nodes
                tillArea=0
                do i=1, Num_tillNodes
                    tillArea=tillArea+nodeArea(node_tillApplied(i))   !Total tillage area [cm2]
                end do

            !*** MIXING BOWL PROCESS 
            ! Currently framed by assuming that the soil has the same material within the tillage depth
            ! Area weighted average redistribution of state variables over the tillage depth
            ! Area associated with each nodes are already estimated in Grid_Bnd (nodeArea(i) [cm2])
            ! State variables: (a) Water content(ThNew[-]), (b) pressure head(hNew[cm]),(c)temperature (Tmpr(i))
            ! State variables: (C(L,H,M), N(L,H,M), conc(n,1),NH4(n) 
            ! Perform the area weighted average over the unit of C and N in ug/cm3 of soil])
            ! P and K can be added later

            !---WC Redistribution----------------
            !Water content (vol of water/volume of soil); volume of soil=2Darea*unit width [cm3]
            !Total volume of water before mixing
                vwcBefore=0
                do i=1, Num_tillNodes
                    vwcBefore=vwcBefore+nodeArea(node_tillApplied(i))
     !          *ThNew(i)
                    ThTill(i)=ThNew(i)                                !Saved the initial water content before mixing to be accounted in conc redistribution
                end do
              !Average mositure content
		      ThAvg=vwcBefore/tillArea
		      !Reassign this water content back to the tilled nodes
                do i=1, Num_tillNodes
                    ThNew(i)=ThAvg
                end do

              !Find the difference in the total volume of water before and after [for Debug]
		      vwcAfter=0
			    do i=1, Num_tillNodes
				    vwcAfter=vwcAfter+nodeArea(node_tillApplied(i))
     !          *ThNew(i)
			    end do
		  !---End of Water content redistribution-------------------


		    !---Conc Redistribution------------
		      !Total conc in the tillage region before mixing, Conc(i,1) refers to N in N03 [ug/cm3 water]
		      concBefore=0
                do i=1, Num_tillNodes
                    concBefore=concBefore+nodeArea(node_tillApplied(i))
     !          *Conc(i,1)*ThTill(i)                                  !Conc[ug/cm3 water]*water content[cm3water/cm3soil]=concBefore[ug/cm3soil]over the tillage area
                    !Here the water content before mixing is used, it will be ressigned based on the new water content after mixing
                end do
		      !Average concentration
		      concAvg=concBefore/tillArea                                   !concBefore[ug/slab]/tillArea[cm2]= concAvg[ug/cm3soil]
		      !Reassign this concentration back to the tilled nodes
                do i=1, Num_tillNodes
                    Conc(i,1)=ConcAvg/ThNew(i)                        !concAvg[ug/cm3soil]/ThNew[cm3water/cm3soil]=Conc(i,1)[ug/cm3water]
                end do
		      !Find the difference in the total conc before and after [For Debug]
		      concAfter=0
                do i=1, Num_tillNodes
                    ConcAfter=concAfter+nodeArea(node_tillApplied(i))
     !          *Conc(i,1)*ThNew(i)
                end do
		    !---End of Conc redistribution----------------------


      
		      !----Nh,NL,Nm,Ch,CL,Cm,NH4,Denit Redistribution [All units are ug/g of soil, nitrogen is in N form)------------
		      NhBefore=0
		      NLBefore=0
		      NmBefore=0
		      ChBefore=0
		      CLBefore=0
		      CmBefore=0
		      NH4Before=0
		      DenitBefore=0
                
cdt no need to use BD as units remain as g per unit volume
                do i=1, Num_tillNodes                                 
                    NhBefore=NhBefore+nodeArea(node_tillApplied(i))!Nh[ug N/Unit volume]=Nh[ug/cm3soil]over the tillage area
     !              *Nh(i)                          
                    NLBefore=NLBefore+nodeArea(node_tillApplied(i))
     !              *NL(i)
                    NmBefore=NmBefore+nodeArea(node_tillApplied(i))
     !              *Nm(i)
                    ChBefore=ChBefore+nodeArea(node_tillApplied(i))
     !              *Ch(i)
                    CLBefore=CLBefore+nodeArea(node_tillApplied(i))
     !              *CL(i)
                    CmBefore=CmBefore+nodeArea(node_tillApplied(i))
     !              *Cm(i)
                    NH4Before=NH4Before+nodeArea(node_tillApplied(i))
     !              *NH4(i)
                    DenitBefore=DenitBefore
     !              +nodeArea(node_tillApplied(i))*denit(i)
                end do
		  !Average concentration
		      NhAvg=NhBefore/tillArea   !NhBefore[ug N/slab]/tillArea[cm2]= NhAvg[ug/unit volume]
		      NLAvg=NLBefore/tillArea
		      NmAvg=NmBefore/tillArea

		      ChAvg=ChBefore/tillArea
		      CLAvg=CLBefore/tillArea
		      CmAvg=CmBefore/tillArea

		      NH4Avg=NH4Before/tillArea
		      DenitAvg=DenitBefore/tillArea

		  !Reassign this concentration back to the tilled nodes
                do i=1, Num_tillNodes
                    Nh(i)=NhAvg                !NhAvg[ug/unitVol]
                    NL(i)=NLAvg
                    Nm(i)=NmAvg

                    Ch(i)=ChAvg
                    CL(i)=CLAvg
                    Cm(i)=CmAvg

                    NH4(i)=NH4Avg
                    Denit(i)=DenitAvg
                end do
		  !Find the difference in the total  N before and after [For Debug]
		      NhAfter=0
		      NLAfter=0
		      NmAfter=0
		      ChAfter=0
		      CLAfter=0
		      CmAfter=0
		      NH4After=0
		      DenitAfter=0
                do i=1, Num_tillNodes
                    NhAfter=NhAfter+nodeArea(node_tillApplied(i))
     !          *Nh(i)                          
                    NLAfter=NLAfter+nodeArea(node_tillApplied(i))
     !          *NL(i)
                    NmAfter=NmAfter+nodeArea(node_tillApplied(i))
     !          *Nm(i)
                    ChAfter=ChAfter+nodeArea(node_tillApplied(i))
     !          *Ch(i)
                    CLAfter=CLAfter+nodeArea(node_tillApplied(i))
     !          *CL(i)
                    CmAfter=CmAfter+nodeArea(node_tillApplied(i))
     !          *Cm(i)
                    NH4After=NH4After+nodeArea(node_tillApplied(i))
     !          *NH4(i)
                    DenitAfter=DenitAfter+nodeArea(node_tillApplied(i))
     !          *denit(i)
                end do
		      !---End of Nh,NL,Nm,Ch,CL,Cm,NH4,Denit Redistribution-------
		  
      
      
                !---Temperature Redistribution------------
                TmprBefore=0
                do i=1, Num_tillNodes
c                    write(*,*) 'tb',i,Tmpr(i)
                    TmprBefore=TmprBefore+nodeArea(node_tillApplied(i))
     !              *Tmpr(i)                                           
                end do

                !Average temperature
                TmprAvg=TmprBefore/tillArea                            

                !Reassign this temperature to the tilled nodes
                do i=1, Num_tillNodes
                    Tmpr(i)=TmprAvg                                    
                end do

                TmprAfter=0                               ![For Debug]
                do i=1, Num_tillNodes
                    TmprAfter=TmprAfter+nodeArea(node_tillApplied(i))
     !              *Tmpr(i)
                end do
                !---End of temperature redistribution-----

      
      
                !---Head Redistribution----------------
                HeadBefore=0
                do i=1, Num_tillNodes
                    HeadBefore=HeadBefore+nodeArea(node_tillApplied(i))
     !              *LOG10(abs(hNew(i)))                          !SB: incase of ponding or saturated contn, h will be 0 or +ve?
                end do
                !Average head
                HeadAvg=HeadBefore/tillArea

                !Reassign head to the tilled nodes
                do i=1, Num_tillNodes
                    hNew(i)=-1.0*10**HeadAvg
                end do

                HeadAfter=0
                do i=1, Num_tillNodes
                    HeadAfter=HeadAfter+nodeArea(node_tillApplied(i))
     !              *LOG10(abs(hNew(i)))
                end do
                !---End of pressure head redistribution------
     
       
          tNext(ModNum)=1.E+32                                    !Since tillage is applied only once
      end if  !Tillage mixing end
		  Return
20    Stop 'Tillage data error'
      end


