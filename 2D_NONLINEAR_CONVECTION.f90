!****************************************************************************                                                                         *
! *   FILE         = 2D_NONLINAER_CONVECTION.F90                            *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************
     PROGRAM TWOD_NONLINEAR_CONVECTION 
     IMPLICIT NONE 
     INTEGER, PARAMETER :: NNODES = 100
     REAL :: U(NNODES,NNODES),UN(NNODES,NNODES),V(NNODES,NNODES),VN(NNODES,NNODES)
     REAL :: C,SIGMA,DX,DY,DT
     INTEGER :: I,J,K,NT,NX,NY
   
     
     !NX- GRID IN X DIRECTION , NY- GRID IN Y DIRECTION
     NX = 80
     NY = 80

     !TOTAL TIME  = NT*DT
     NT =  50
     C = 1.0
     SIGMA = 0.2
     
     DX=2.0/(NX-1.0)
     DY=2.0/(NY-1.0)
     DT=SIGMA*DX
     DO I=1,(NX-1)
        DO J=1,(NY-1)
         U(I,J)=1.0
         V(I,J)=1.0
        END DO
     END DO
     DO I=(NX/5),(NX/5+20)
       DO J=(NY/5),(NY/5+20)
         U(I,J)=2.0
         V(I,J)=2.0
       END DO
     END DO
     OPEN(1,FILE='INITIAL.dat')
      WRITE(1,*) " TITLE =  ""FLOW-FIELD"" "
      WRITE(1,*) "VARIABLES = X, Y, U, V"
      WRITE(1,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
        DO J=1,(NY-1)
        WRITE(1,*)I,J,U(I,J),V(I,J)
        END DO
     END DO
     CLOSE(1)
     
     DO K=1,NT
     DO I=2,(NX-1)
       DO J=2,(NY-1)
         UN(I,J)=U(I,J)-(U(I,J)*DT/DX)*(U(I,J)-U(I-1,J))-(V(I,J)*DT/DY)*(U(I,J)-U(I,J-1))
         VN(I,J)=V(I,J)-(U(I,J)*DT/DX)*(V(I,J)-V(I-1,J))-(V(I,J)*DT/DY)*(V(I,J)-V(I,J-1))
       END DO
     END DO
     DO I=2,(NX-1)
       DO J=2,(NY-1)
         U(I,J)=UN(I,J)
         V(I,J)=VN(I,J)
        END DO
     END DO
     END DO

      OPEN(2,FILE='FINAL.dat')
      WRITE(2,*) " TITLE = ""FLOW-FIELD"" "
      WRITE(2,*)  "VARIABLES = X, Y, U, V"
      WRITE(2,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
       DO J=1,(NY-1)
        WRITE(2,*)I,J,U(I,J),V(I,J)
       END DO
     END DO
     PRINT*, " FILES FOR INITIAL AND FINAL DATA ARE GENERATED."
     PRINT*, " USE PARAVIEW OR TECPLOT OR ANY OTHER TOOL TO VISUALIZE."


     CLOSE(2)
    END PROGRAM TWOD_NONLINEAR_CONVECTION

