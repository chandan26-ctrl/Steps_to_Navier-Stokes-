!****************************************************************************                                                                         *
! *   FILE         = 2D_POISSON.F90                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************  

     PROGRAM TWOD_POISSON
     IMPLICIT NONE
     REAL, ALLOCATABLE :: P(:,:),PN(:,:),B(:,:)
     INTEGER :: NY,I,J,K
     REAL :: DX,DY
     INTEGER, PARAMETER :: NX = 80, NT = 800!PSEUDO TIME
    
     NY=NX
    
     ALLOCATE(P(NX+1,NY+1),PN(NX+1,NY+1),B(NX+1,NY+1))
     DX=2.0/NX
     DY=2.0/NY
     
     ! INITIAL CONDITION
     DO I=1,NX
       DO J=1,NY
         P(I,J)=0.0
         B(I,J)=0.0
       END DO
     END DO

   
     ! BOUNDARY CONDITION
     DO I=1,NX
        P(I,1)=0.0
        P(I,NY)=0.0
     END DO
     DO J=1,NY
        P(1,J)=0.0
        P(NX,J)=0.0
     END DO
     B(NX/4,NY/4)=100.0
     B(3*NX/4,3*NY/4)=-100.0 !! INITAL SPIKE 


     OPEN(1,FILE='INITIAL.dat')
      WRITE(1,*) " TITLE =  ""FLOW-FIELD"" "
      WRITE(1,*) "VARIABLES = X, Y, P"
      WRITE(1,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NY, ",F=POINT"
     DO I=1,NX
        DO J=1,NY
        WRITE(1,*)I,J,P(I,J)
        END DO
     END DO
     CLOSE(1)

     
     DO K=1,NT
      DO I=2,NX-1
        DO J=2,NY-1
           PN(I,J)=((P(I+1,J)+P(I-1,J))*DY*DY)+((P(I,J+1)+P(I,J-1))*DX*DX)-(B(I,J)*DX*DX*DY*DY)/ &
                     (2*(DX**2+DY**2))
        END DO
      END DO
  
      DO I=2,NX-1
        DO J=2,NY-1
         P(I,J)=PN(I,J)
        END DO
      END DO

     END DO


      OPEN(2,FILE='FINAL.dat')
      WRITE(2,*) " TITLE = ""FLOW-FIELD"" "
      WRITE(2,*)  "VARIABLES = X, Y, P"
      WRITE(2,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NY, ",F=POINT"
     DO I=1,NX
       DO J=1,NY
        WRITE(2,*)I,J,P(I,J)
       END DO
     END DO
     PRINT*, " FILES FOR INITIAL AND FINAL DATA ARE GENERATED."
     PRINT*, " USE PARAVIEW OR TECPLOT OR ANY OTHER TOOL TO VISUALIZE."
    
     CLOSE(2)


     END PROGRAM TWOD_POISSON
     
