!****************************************************************************                                                                         *
! *   FILE         = 1D_NONLINAER_CONVECTION.F90                            *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************

     PROGRAM ONED_NONLINEAR_CONVECTION
     IMPLICIT NONE
 
     REAL,DIMENSION(41) :: U,UN
     REAL, PARAMETER :: DT=0.025,NX=41.0,C=1.0
     INTEGER, PARAMETER :: NT = 20
     REAL :: DX
     INTEGER :: I,J

     DX=2.0/(NX-1.0)

     
     !INITIAL CONDITION AS A STEP FUNCTION
     DO I=0,9
         U(I)=1.0          !         -----------
     END DO                !         |         |
                           !         |         |
     DO I=10,20            !         |         |
         U(I)=2.0          ! ---------         ----------
     END DO

     DO I=21,40
        U(I)=1.0
     END DO

     OPEN(2,FILE="DATA.DAT")
     DO J= 0,NT
       PRINT*,"AT TIME", DT*J,"SECONDS, VARIABLES = ""X"",""Y"" "
       WRITE(2,*),"AT TIME", DT*J,"SECONDS, VARIABLES = ""X"",""Y"" "
       DO I=1,40
          UN(I)=U(I)-(U(I)*DT/DX)*(U(I)-U(I-1)) !!! IF CHANGE IN DX OR DT IS REQUIRED, IT SHOULD SATISFY THE CFL CONDITION. 
          PRINT*,I,UN(I)
       END DO

       DO I=1,40
        U(I)=UN(I)
       END DO
     END DO
    CLOSE(2)

    END PROGRAM ONED_NONLINEAR_CONVECTION
