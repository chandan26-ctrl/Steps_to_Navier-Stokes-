!****************************************************************************                                                                         *
! *   FILE         = 2D_LINAER_CONVECTION.F90                               *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************

     PROGRAM TWOD_LINEAR_CONVECTION     
     IMPLICIT NONE 

     REAL,DIMENSION(82,82) :: U,UN
     REAL, PARAMETER:: C=1.0,SIGMA=0.2
     INTEGER, PARAMETER :: NT = 100, NX = 81, NY = 81
     REAL :: DX,DY,DT
     INTEGER :: I,J,K

     DX=2.0/(NX-1.0)
     DY=2.0/(NY-1.0)
     DT=SIGMA*DX


     ! A SQUARE INITIAL FUNCTION

     DO I=1,80
        DO J=1,80
         U(I,J)=1.0
        END DO
     END DO

     DO I=20,40
       DO J=20,40
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
     DO I=2,80
       DO J=2,80
         UN(I,J)=U(I,J)-(C*DT/DX)*(U(I,J)-U(I-1,J))-(C*DT/DY)*(U(I,J)-U(I,J-1))
         !!! IF CHANGE IN DX OR DT IS REQUIRED, IT SHOULD SATISFY THE CFL CONDITION. 
       END DO
     END DO

     DO I=2,80
       DO J=2,80
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
    END PROGRAM TWOD_LINEAR_CONVECTION

