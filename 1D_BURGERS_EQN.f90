!****************************************************************************                                                                         *
! *   FILE         = 1D_BURGERS_EQN.F90                                     *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!****************************************************************************
     PROGRAM ONED_BURGERS_EQN
     IMPLICIT NONE
     REAL,DIMENSION(102) :: U,UN,X
     REAL, PARAMETER :: SIGMA=0.2,NU =0.07
     INTEGER, PARAMETER :: NT = 100, NX = 101    
     DOUBLE PRECISION :: PI
     REAL :: DX,DT
     INTEGER :: I,J

     PI=4.D0*DATAN(1.D0)
     DX=2.0*PI/(NX-1.0)
     DT = DX*NU

      X(1)=0.0
     DO I=2,NX-1
          X(I)=X(I-1)+DX
     END DO

       X(NX)=2*PI


     !INITIAL CONDITION 
     DO I=1,NX
         U(I)=((-2*NU)*((-X(I)*EXP(-X(I)**2/(4*NU))/(2*NU))-2*(X(I)-2*PI)*EXP(-(X(I) &
         -2*PI)**2/(4*NU))/(4*NU))/(EXP(-X(I)**2/(4*NU))+EXP(-(X(I)-2*PI)**2/(4*NU))))+4
     END DO

     PRINT*, "INITIAL" 
     DO I=1,NX
         PRINT*, X(I),U(I)
     END DO

     DO J= 0,NT
      DO I=2,NX-1
       UN(I)=U(I)-(U(I)*DT/DX)*(U(I)-U(I-1))+(NU*DT/DX**2)*(U(I+1)-2*U(I)+U(I-1))
      END DO
      DO I=2,NX-1
        U(I)=UN(I)
      END DO
     END DO

     PRINT*, "FINAL"
     DO I=1,NX      
         PRINT*, X(I),U(I)
     END DO

     OPEN(3,FILE='DATA.DAT')
     DO I=1,NX    
         WRITE(3,*) X(I),U(I)
     END DO
     CLOSE(3)

     END PROGRAM ONED_BURGERS_EQN
