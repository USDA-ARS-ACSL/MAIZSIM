* Source file ORTHOFEM.FOR |||||||||||||||||||||||||||||||||||||||||||||
*
*                            ORTHOFEM
*
*                           VERSION 1.02
*
*                      FORTRAN SUBROUTINES FOR
*                   ORTHOMIN OR CONJUGATE GRADIENT 
*               MATRIX SOLUTION ON FINITE-ELEMENT GRIDS
*
*     ***********************************************************
*
*                           CARL A. MENDOZA
*
*                      WITH CONTRIBUTIONS FROM:
*                            RENE THERRIEN
*
*                     BASED ON AN ORIGINAL CODE BY:
*                         FRANK W. LETNIOWSKI
*
*              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
*                       UNIVERSITY OF WATERLOO
*                         WATERLOO, ONTARIO
*                          CANADA, N2L 3G1
*
*                    LATEST UPDATE: JANUARY 1991
*
*
*     ***********************************************************
*
*                   COPYRIGHT (c) 1989, 1990, 1991
*                            E.A. SUDICKY
*                            C.A. MENDOZA
*              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
*
*          DUPLICATION OF THIS PROGRAM, OR ANY PART THEREOF,
*          WITHOUT THE EXPRESS PERMISSION OF THE COPYRIGHT
*          HOLDERS IS STRICTLY FORBIDDEN
*
*     ***********************************************************
*
*               Modified for SWMS_2D code by Jirka Simunek
*                            april 1993
*
*     ***********************************************************
*
*                             DISCLAIMER
*
*     ALTHOUGH GREAT CARE HAS BEEN TAKEN IN PREPARING THIS CODE 
*     AND THE ACCOMPANYING DOCUMENTATION, THE AUTHORS CANNOT BE
*     HELD RESPONSIBLE FOR ANY ERRORS OR OMISSIONS.  THE USER IS 
*     EXPECTED TO BE FAMILIAR WITH THE FINITE-ELEMENT METHOD,
*     PRECONDITIONED ITERATIVE TECHNIQUES AND FORTRAN PROGRAMING.
*
*     *********************************************************** 
*
*             A USER'S GUIDE IS AVAILABLE -> CONSULT IT!
*
*     THESE SUBROUTINES SOLVE A BANDED (OR SPARSE) MATRIX USING:
*       - PRECONDITIONED ORTHOMIN FOR ASYMMETRIC MATRICES, OR
*       - PRECONDITIONED CONJUGATE GRADIENT FOR SYMMETRIC MATRICES
*         (FULL MATRIX STORAGE REQUIRED)
*
*     PRECONDITIONING IS BY INCOMPLETE LOWER-UPPER DECOMPOSITION
*       - ONLY ONE FACTORIZATION (GAUSSIAN ELIMINATION) IS PERFORMED
*       - EQUIVALENT TO DKR FACTORIZATION
*
*     THE SUBROUTINES ARE DESIGNED FOR FINITE-ELEMENT GRIDS
*       - ARBITRARY ELEMENT SHAPES AND NUMBERING MAY BE USED
*         - NUMBERING MAY, HOWEVER, AFFECT EFFICIENCY
*           - TRY TO MINIMIZE THE BANDWIDTH AS MUCH AS POSSIBLE
*         - ALL ELEMENTS MUST HAVE THE SAME NUMBER OF LOCAL NODES
*
*
*     THE FOLLOWING ROUTINES ARE CALLED FROM THE SOURCE PROGRAM:
*       IADMAKE (IN,NN,NE,NLN,MNLN,MAXNB,IAD,IADN,IADD)
*         -> ASSEMBLE ADJACENCY MATRIX
*       FIND (I,J,K,NN,MAXNB,IAD,IADN)
*         -> LOCATE MATRIX POSITION FOR A NODAL PAIR (ASSEMBLY)
*       ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
*         -> DECOMPOSE GLOBAL MATRIX
*       ORTHOMIN (R,C,GT,NNR,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,
*                 RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
*                 NORTH,MNORTH,MAXIT)
*         -> SOLVE DECOMPOSED MATRIX
*
*     THESE ROUTINES CALL OTHER ROUTINES (LOCATED DIRECTLY BELOW THE 
*     APPROPRIATE PRIMARY ROUTINE IN THE CODE)
* 
*     THE FOLLOWING ARRAYS MUST BE DEFINED IN THE SOURCE PROGRAM
*     (THESE ARRAYS ARE PASSED TO THE SOLVER SUBROUTINES):
*
*     IN(MNLN,MAXNE) - INCIDENCE MATRIX (ELEMENTAL NODE DEFINITION)
*
*     GT(MAXNN)      - RIGHT-HAND-SIDE VECTOR
*     C(MAXNN)       - SOLUTION VECTOR
*     R(MAXNB,MAXNN) - GLOBAL MATRIX TO BE SOLVED
*
*     ARRAY DIMENSIONING PARAMETERS
*
*     MAXNN  - MAXIMUM NUMBER OF NODES
*     MAXNE  - MAXIMUM NUMBER OF ELEMENTS
*     MNLN   - MAXIMUM NUMBER OF LOCAL NODES (IN AN ELEMENT)
*     MAXNB  - MAXIMUM NUMBER OF NODES ADJACENT TO A PARTICULAR NODE
*              (INCLUDING ITSELF).
*            - IE. THE MAXIMUM NUMBER OF INDEPENDENT NODES THAT A 
*              PARTICULAR NODE SHARES AN ELEMENT WITH.
*            - THIS WILL BE IDENTICALLY EQUIVALENT TO THE MAXIMUM 
*              NUMBER OF NONZERO ENTRIES IN A ROW OF THE FULL MATRIX.
*     MNORTH - MAXIMUM NUMBER OF ORTHOGONALIZATIONS PERFORMED
*              (AT LEAST MNORTH = 1 REQUIRED FOR CONJUGATE GRADIENT)
*
*
*     ORTHOMIN ARRAY SPACE/VARIABLES
*
*     NORTH  - NUMBER OF ORTHOGONALIZATIONS TO PERFORM
*            - SET NORTH=0 FOR SYMMETRIC MATRICES (CONJUGATE GRADIENT)
*     ECNVRG - RESIDUAL CONVERGENCE TOLERANCE
*     ACNVRG - ABSOLUTE CONVERGENCE TOLERANCE
*     RCNVRG - RELATIVE CONVERGENCE TOLERANCE
*     MAXIT  - MAXIMUM NUMBER OF ITERATIONS TO PERFORM
*     ITERP  - NUMBER OF ITERATIONS PERFORMED
*
*     B(MAXNB,MAXNN) - ILU DECOMPOSED MATRIX
*     Q(MAXNN)   - SEARCH DIRECTION Q
*     RQ(MAXNN)  - PRODUCT OF R AND Q
*     VRV(MAXNN) - EITHER V OR PRODUCT OF R AND V
*     RES(MAXNN) - RESIDUAL
*
*     QI(MAXNN,MNORTH)  - STORAGE OF Q'S
*     RQI(MAXNN,MNORTH) - STORAGE OF PRODUCTS OF R AND Q
*     RQIDOT(MNORTH)    - STORAGE OF DOT PRODUCTS OF RQ AND RQ
*
*     RESV   - PREVIOUS VALUE OF RES V DOT PRODUCT (CONJUGATE GRADIENT)
*
*     IAD(MAXNB,MAXNN) - ADJACENCY MATRIX (NODAL CONNECTIONS)
*
*     IADN(MAXNN) - NUMBER OF ADJACENT NODES IN IAD (SELF-INCLUSIVE)
*
*     IADD(MAXNN) - POSITION OF DIAGONAL IN ADJACENCY MATRIX
*
*
*     OTHER PARAMETERS PASSED FROM SOURCE PROGRAM
*
*     NN  - NUMBER OF NODES
*     NE  - NUMBER OF ELEMENTS
*     NLN - NUMBER OF LOCAL NODES IN AN ELEMENT
*
*
*     APPROXIMATE REAL STORAGE SPACE FOR ORTHOMIN AND MATRIX EQUATION
*
*       ((6 + 2*MAXNB + 2*MNORTH)*MAXNN)*(8 BYTES)
*
************************************************************************

      subroutine IADMake(KX,NumNP,NumEl,NumElD,MaxNB,IAD,IADN,IADD)

*     Generate the adjacency matrix for nodes from the element 
*     indidence matrix

*     Requires subroutine Insert

      dimension KX(NumElD,4),IAD(MaxNB,NumNP),IADN(NumNP),IADD(NumNP)
      integer e

*     Determine independent adjacency within each element
*     version for SWMS_2D

      IADN = 0
      IADD = 0
      IAD  = 0
!      do 12 i=1,NumNP
!        IADN(i)=0
!        IADD(i)=0
!        do 11 j=1,MaxNB
!          IAD(j,i)=0
!11      continue
!12    continue

      do 14 e=1,NumEl
        NCorn=4
        if(KX(e,3).eq.KX(e,4)) NCorn=3
        do 13 n=1,NCorn-2
          i=KX(e,1)
          j=KX(e,n+1)
          k=KX(e,n+2)

          call Insert(i,j,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(j,i,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(i,k,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(k,i,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(j,k,kk,NumNP,MaxNB,IAD,IADN)
          call Insert(k,j,kk,NumNP,MaxNB,IAD,IADN)
13      continue
14    continue

*     Determine self-adjacency terms

      do 15 i=1,NumNP
        call Insert(i,i,kk,NumNP,MaxNB,IAD,IADN)
*     Store self-adjacency position
        IADD(i)=kk
15    continue
      return
      end

************************************************************************

      SUBROUTINE INSERT (I,J,K,NN,MAXNB,IAD,IADN)

*     ADD J TO THE ADJACENCY LIST FOR I

*     RETURNS THE POSITION K WHERE IT HAS BEEN ADDED, OR WHERE IT 
*     WAS ALREADY IN THE LIST.

      DIMENSION IAD(MAXNB,NN),IADN(NN)

      LOGICAL FOUND
      FOUND = .FALSE.

*     DETERMINE NUMBER OF NODES ALREADY IN ADJACENCY LIST

      N = IADN(I)
      K = N + 1

*     DETERMINE WHETHER ALREADY IN LIST

      DO 10 L=1,N
        INODE = IAD(L,I)
        IF (INODE.GE.J) THEN
          K = L
          IF (INODE.EQ.J) FOUND = .TRUE.
          GO TO 15
        ENDIF
   10 CONTINUE

   15 CONTINUE

*     PLACE IN LIST (NUMERICAL ORDER)

      IF (FOUND) THEN
       CONTINUE
      ELSE
        IF ((N+1).GT.MAXNB) THEN
          WRITE (* ,601) I,MAXNB
          WRITE (50,601) I,MAXNB
  601     FORMAT (//5X,'ERROR IN IADMAKE: NODE ',I5,' HAS > '
     *            ,I5,' ADJACENCIES')
          STOP
        ENDIF

        IADN(I) = N + 1
        DO 20 L=(N+1),(K+1),(-1)
          IAD(L,I) = IAD(L-1,I)
   20   CONTINUE
       IAD(K,I) = J
      ENDIF

      RETURN
      END

************************************************************************

      SUBROUTINE FIND (I,J,K,NN,MAXNB,IAD,IADN)

*     FOR NODE I, DETERMINE THE 'BAND' (K) RELATED TO ITS ADJACENCY TO
*     NODE J.

*     IF NODE NOT ADJACENT, RETURN 0 AS THE 'BAND'

C      DIMENSION IAD(MAXNB,NN),IADN(NN)
C GR
      DIMENSION IAD(MAXNB,NN),IADN(I)

      K = 0
      N = IADN(I)

      DO 10 L=1,N
        INODE = IAD(L,I)

*     EXIT THE LOOP IF AT OR PAST THE REQUIRED POSITION

        IF (INODE.GE.J) THEN
          IF (INODE.EQ.J) K = L
          GO TO 20
        ENDIF
   10 CONTINUE

   20 CONTINUE

      RETURN
      END

************************************************************************

      SUBROUTINE ILU (R,NN,MAXNB,IAD,IADN,IADD,B)

*     INCOMPLETE LOWER-UPPER DECOMPOSITION OF MATRIX R INTO B
*     ONE STEP OF GAUSSIAN ELIMINATION PERFORMED
*     DIAGONAL DOMINANCE IS ASSUMED - NO PIVOTING PERFORMED
*     REQUIRES FUNCTION DU

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)

*     INITIALIZE B

      DO 10 I=1,NN
        DO 10 J=1,MAXNB
          B(J,I) = 0.0D0
   10 CONTINUE

*     LOOP OVER NODES

      DO 20 I=1,NN

*     DETERMINE NUMBER OF BANDS/POSITION OF DIAGONAL IN THIS ROW

        N = IADN(I)
        K = IADD(I)

*     LOWER TRIANGULAR MATRIX

        DO 30 J=1,(K-1)
          SUM = R(J,I)
          ICUR = IAD(J,I)
          DO 40 L=1,(J-1)
            INODE = IAD(L,I)
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
   40     CONTINUE
          B(J,I) = SUM
   30   CONTINUE

*     DIAGONAL

        SUM = R(K,I)
        DO 50 L=1,(K-1)
          INODE = IAD(L,I)
          SUM = SUM - B(L,I)*DU(INODE,I,NN,MAXNB,IAD,IADN,IADD,B)
   50   CONTINUE
        D = 1.0D0/SUM
        B(K,I) = D

*     UPPER TRIANGULAR MATRIX
*       - ACTUALLY D*U TO OBTAIN UNIT DIAGONAL

        DO 60 J=(K+1),N
          SUM = R(J,I)
          ICUR = IAD(J,I)
          DO 70 L=1,(K-1)
            INODE = IAD(L,I)
            SUM = SUM - B(L,I)*DU(INODE,ICUR,NN,MAXNB,IAD,IADN,IADD,B)
   70     CONTINUE
          B(J,I) = D*SUM
   60   CONTINUE

   20 CONTINUE

      RETURN
      END

************************************************************************

      FUNCTION DU (I,INODE,NN,MAXNB,IAD,IADN,IADD,B)

*     SEARCHES THE I'TH ROW OF THE UPPER DIAGONAL MATRIX
*     FOR AN ADJACENCY TO THE NODE 'INODE'

*     RETURNS CORRESPONDING VALUE OF B (OR ZERO)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN)

      TEMP = 0.0D0
      N = IADN(I)
      K = IADD(I)

      IF (I.EQ.INODE) THEN
        TEMP = 1.0D0
        GO TO 20
      ENDIF

      DO 10 J=(K+1),N
        IF (INODE.EQ.IAD(J,I)) THEN
          TEMP = B(J,I)
          GO TO 20
        ENDIF
   10 CONTINUE

   20 CONTINUE
      DU = TEMP

      RETURN
      END

************************************************************************

      SUBROUTINE ORTHOMIN (R,C,GT,NN,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,
     !                     RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
     !                     NORTH,MNORTH,MAXIT, modNum)

*     ORTHOMIN OR CONJUGATE GRADIENT ACCELERATION/SOLUTION

*     CONJUGATE GRADIENT (SYMMETRIC MATRIX) IF NORTH=0
*     (HOWEVER, NOTE THAT MNORTH MUST BE AT LEAST 1)

*     REQUIRES FUNCTIONS SDOT,SDOTK,SNRM
*     REQUIRES SUBROUTINES LUSOLV,MATM2,SAXPYK,SCOPY,SCOPYK 

      IMPLICIT REAL*8 (A-H,O-Z)
      integer modnum
      DIMENSION R(MAXNB,NN),C(NN),GT(NN),IAD(MAXNB,NN),IADN(NN),
     !          IADD(NN),B(MAXNB,NN),VRV(NN),RES(NN),RQI(MAXNN,MNORTH),
     !          RQ(NN),Q(NN),RQIDOT(MNORTH),QI(MAXNN,MNORTH)

*     INITIALIZE RESIDUAL VECTOR

      CALL MATM2 (RES,R,C,NN,IAD,IADN,MAXNB)

*     Solution for homogeneous system of equations - Modified by Simunek
      iHomog=0
      DO 10 I=1,NN
        RES(I) = GT(I) - RES(I)
        if(abs(GT(i)).gt.1.d-300) iHomog=1
   10 CONTINUE
      if(iHomog.eq.0) then
        do 11 i=1,NN
	  C(i)=GT(i)
11      continue
        return
      end if     

*     LOOP OVER A MAXIMUM OF MAXIT ITERATIONS

      NORCUR = 0

      DO 100 ITER=1,MAXIT

*     INVERT LOWER/UPPER MATRICES

        CALL SCOPY (NN,RES,VRV)

        CALL LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)

*     COPY V INTO Q

        CALL SCOPY (NN,VRV,Q)

*     CALCULATE PRODUCT OF R AND V

        CALL MATM2 (VRV,R,Q,NN,IAD,IADN,MAXNB)

*     COPY RV INTO RQ

        CALL SCOPY (NN,VRV,RQ)

*     RES V DOT PRODUCT (CONJUGATE GRADIENT)

        IF (NORTH.EQ.0) THEN

          DOT = SDOT(NN,RES,Q)
          IF (NORCUR.EQ.0) RESV = DOT

        ENDIF

*     LOOP OVER PREVIOUS ORTHOGONALIZATIONS

        K = 1

   20   IF (K.GT.NORCUR) GO TO 30

*     DETERMINE WEIGHTING FACTOR (CONJUGATE GRADIENT)

          IF (NORTH.EQ.0) THEN

            ALPHA = DOT/RESV
            RESV = DOT

*     DETERMINE WEIGHTING FACTOR (ORTHOMIN)

          ELSE

            DOT = SDOTK(NN,K,RQI,VRV,MAXNN,MNORTH)
            ALPHA = -DOT/RQIDOT(K)

          ENDIF

*     SUM TO OBTAIN NEW Q AND RQ

          CALL SAXPYK (NN,ALPHA,K,QI,Q,MAXNN,MNORTH)
          CALL SAXPYK (NN,ALPHA,K,RQI,RQ,MAXNN,MNORTH)

          K = K + 1
          GO TO 20
   30   CONTINUE

*     CALCULATE WEIGHTING FACTOR (CONJUGATE GRADIENT)

        IF (NORTH.EQ.0) THEN

          DOT = SDOT(NN,Q,RQ)
          OMEGA = RESV/DOT

*     CALCULATE WEIGHTING FACTOR (ORTHOMIN)

        ELSE

          DOT = SDOT(NN,RES,RQ)
          RQNORM = SDOT(NN,RQ,RQ)
          OMEGA = DOT/RQNORM

        ENDIF

*     SAVE VALUES FOR FUTURE ORTHOGONALIZATIONS

        NORCUR = NORCUR + 1
        IF (NORCUR.GT.NORTH) NORCUR = 1

        CALL SCOPYK (NN,NORCUR,Q,QI,MAXNN,MNORTH)
        CALL SCOPYK (NN,NORCUR,RQ,RQI,MAXNN,MNORTH)
        RQIDOT(NORCUR) = RQNORM
 
*     UPDATE SOLUTION/RESIDUAL VECTORS

        CALL SAXPYK (NN,OMEGA,NORCUR,QI,C,MAXNN,MNORTH)
        CALL SAXPYK (NN,-OMEGA,NORCUR,RQI,RES,MAXNN,MNORTH)

*     DETERMINE CONVERGENCE PARAMETERS

        RESMAX = SNRM(NN,RES)
        DXNORM = DABS(OMEGA)*SNRM(NN,Q)
        CNORM = SNRM(NN,C)
        XRATIO = DXNORM/CNORM

*     ITERATION (DEBUG) OUTPUT

*     STOP ITERATING IF CONVERGED

        IF (RESMAX.LT.ECNVRG) GO TO 200
        IF (XRATIO.LT.RCNVRG) GO TO 200
        IF (DXNORM.LT.ACNVRG) GO TO 200

  100 CONTINUE

*     TERMINATE IF TOO MANY ITERATIONS

      WRITE (* ,602) MAXIT,RESMAX,XRATIO,DXNORM
      WRITE (70,602) MAXIT,RESMAX,XRATIO,DXNORM
      write (*,*) 'modNum=',modNum

 602  FORMAT (///5X,'ORTHOMIN TERMINATES -- TOO MANY ITERATIONS',
     *          /8X,'MAXIT  = ',I5,
     *          /8X,'RESMAX = ',E12.4,
     *          /8X,'XRATIO = ',E12.4,
     *          /8X,'DXNORM = ',E12.4)

      STOP

  200 CONTINUE

*     RETURN NUMBER OF ITERATIONS REQUIRED

      ITERP = ITER

      RETURN
      END

************************************************************************

      SUBROUTINE LUSOLV (NN,MAXNB,IAD,IADN,IADD,B,VRV)

*     LOWER DIAGONAL MATRIX INVERSION BY FORWARD SUBSTITUTION
*     UPPER DIAGONAL MATRIX INVERSION BY BACKWARD SUBSTITUTION
*     LOWER/UPPER MATRICES ARE IN B
*     RIGHT-HAND-SIDE VECTOR IS IN VRV AT START
*     'SOLUTION' IS RETURNED IN VRV UPON EXIT

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IAD(MAXNB,NN),IADN(NN),IADD(NN),B(MAXNB,NN),VRV(NN)

*     LOWER INVERSION
    
      DO 20 I=1,NN
        SUM = VRV(I)
        K = IADD(I)
        DO 30 J=1,(K-1)
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   30   CONTINUE
   
        VRV(I) = B(K,I)*SUM
   
   20 CONTINUE

*     UPPER INVERSION

      DO 40 I=NN,1,-1
        SUM = VRV(I)
        N = IADN(I)
        K = IADD(I)
        DO 50 J=(K+1),N
          INODE = IAD(J,I)
          SUM = SUM - B(J,I)*VRV(INODE)
   50   CONTINUE

        VRV(I) = SUM

   40 CONTINUE

      RETURN
      END

************************************************************************

      SUBROUTINE MATM2 (S1,R,P,NN,IAD,IADN,MAXNB)

*     MULTIPLY MATRIX R BY VECTOR P TO OBTAIN S1

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S1(NN),P(NN),R(MAXNB,NN),IAD(MAXNB,NN),IADN(NN)

      DO 30 I=1,NN
        SUM = 0.0D0
        N = IADN(I)
        DO 40 J=1,N
          INODE = IAD(J,I)
          SUM = SUM + R(J,I)*P(INODE)
   40   CONTINUE

        S1(I) = SUM
   30 CONTINUE

      RETURN
      END

************************************************************************

      FUNCTION SDOT (NN,R,B)

*     OBTAIN DOT PRODUCT OF R AND B

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(NN),B(NN)

      SDOT= SUM(R * B)
!      SDOT = 0.0D0
!      DO 100 L=1,NN
!        SDOT = SDOT + R(L)*B(L)
!100   CONTINUE

      RETURN
      END

************************************************************************

!      FUNCTION SDOTK2 (NN,K,R,B,MAXNN,MNORTH)
!
!*     OBTAIN DOT PRODUCT OF R AND B
!
!      IMPLICIT REAL*8 (A-H,O-Z)
!      INTEGER MaxNN, MNORTH    
!      DIMENSION R(MAXNN,MNORTH),B(NN)
!
!      SDOTK2 = 0.0D0
!      DO 100 L=1,NN
!        SDOTK2 = SDOTK2 + R(L,K)*B(L)
!100   CONTINUE
!
!      RETURN
!      END
      
      FUNCTION SDOTK(NN, K, R, B,MAXNN, MNORTH) 

      IMPLICIT NONE

      INTEGER :: NN, K, MAXNN, MNORTH       
      REAL*8 :: SDOTK, R(NN, *), B(NN)

      SDOTK = SUM(R(:, K) * B)

      END FUNCTION SDOTK

************************************************************************

      FUNCTION SNRM (NN,R)

*     COMPUTE MAXIMUM NORM OF R

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(NN)

      SNRM = 0.0D0
      SNRM = MAXVAL(DABS(R))
 !     DO 100 L=1,NN
 !       TEMP = DABS(R(L))
 !     IF (TEMP.GT.SNRM) SNRM = TEMP
 !100   CONTINUE

      RETURN
      END

************************************************************************

      SUBROUTINE SAXPYK (NN,SA,K,FX,FY,MAXNN,MNORTH)

*     MULTIPLY VECTOR FX BY SCALAR SA AND ADD TO VECTOR FY

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FX(MAXNN,MNORTH),FY(NN)

       IF (NN > 0) THEN
          FY(1:NN) = SA * FX(1:NN, K) + FY(1:NN)
        ENDIF

!      IF (NN.GT.0) THEN
!        DO 100 I=1,NN
!           FY(I) = SA*FX(I,K) + FY(I)
!100     CONTINUE
!      ENDIF

      RETURN
      END

************************************************************************
      SUBROUTINE SCOPY (NN,FX,FY)

*     COPY VECTOR FX INTO VECTOR FY

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FX(NN),FY(NN)
      
       IF (NN.GT.0) THEN
         FY=FX
       ENDIF
      
!
!      IF (NN.GT.0) THEN
!        DO 100 I=1,NN
!           FY(I) = FX(I)
!100     CONTINUE
!      ENDIF

      RETURN
      END

************************************************************************

      SUBROUTINE SCOPYK (N,K,FX,FY,MAXNN,MNORTH)

*     COPY VECTOR FX INTO VECTOR FY
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FX(N),FY(MAXNN,MNORTH)

       IF (N > 0) THEN
         FY(1:N, K) = FX(1:N)
       ENDIF
!      IF (N.GT.0) THEN
!            DO 100 I=1,N
!           FY(I,K) = FX(I)
!100     CONTINUE
!      ENDIF

      RETURN
      END

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||