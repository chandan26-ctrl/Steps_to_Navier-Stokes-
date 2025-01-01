!****************************************************************************                       
! *   FILE         = 2D_DIFFUSION.F90                                       *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************    

     PROGRAM TWOD_DIFFUSION
     IMPLICIT NONE
     REAL,ALLOCATABLE :: U(:,:), UN(:,:)
     REAL,PARAMETER :: NU=0.05,SIGMA=0.25
     INTEGER :: NY,I,J,K,NNODES
     INTEGER, PARAMETER :: NT = 50, NX = 80
   
     REAL :: DX,DY,DT
     
      NY = NX
     ALLOCATE(U(NX+1,NY+1),UN(NX+1,NY+1))
     DX=2/(NX-1.0) !!0.05
     DY=2/(NY-1.0)
     DT=SIGMA*DX*DY/NU !!0.0125
     DO I=1,NX-1
       DO J=1,NY-1
        U(I,J)=1.0
       END DO
     END DO
     
     DO I=(NX-1)/4,(NX-1)/5+50
      DO J=(NY-1)/4,(NY-1)/5+50
        U(I,J)=2.0
      END DO
     END DO

     OPEN(1,FILE='INITIAL.dat')
      WRITE(1,*) " TITLE =  ""FLOW-FIELD"" "
      WRITE(1,*) "VARIABLES = X, Y, U"
      WRITE(1,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
        DO J=1,(NY-1)
        WRITE(1,*)I,J,U(I,J)
        END DO
     END DO
     CLOSE(1)

    
     DO K=1,NT
       DO I=2,NX-2
         DO J=2,NY-2
          UN(I,J)= U(I,J)+(NU*DT/DX**2)*(U(I+1,J)-2*U(I,J)+U(I-1,J))+ &
                    (NU*DT/DY**2)*(U(I,J+1)-2*U(I,J)+U(I,J-1))
         END DO
       END DO
       DO I=2,NX-2
         DO J=2,NY-2
           U(I,J)=UN(I,J)
         END DO
       END DO
     END DO



      OPEN(2,FILE='FINAL.dat')
      WRITE(2,*) " TITLE = ""FLOW-FIELD"" "
      WRITE(2,*)  "VARIABLES = X, Y, U"
      WRITE(2,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
       DO J=1,(NY-1)
        WRITE(2,*)I,J,U(I,J)
       END DO
     END DO
     PRINT*, " FILES FOR INITIAL AND FINAL DATA ARE GENERATED."
     PRINT*, " USE PARAVIEW OR TECPLOT OR ANY OTHER TOOL TO VISUALIZE."

     CLOSE(2)
    END PROGRAM TWOD_DIFFUSION
         
