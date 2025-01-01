!****************************************************************************                  
! *   FILE         = 2D_LAPLACE.F90                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!**************************************************************************** 

     PROGRAM TWOD_LAPLACE
     IMPLICIT NONE
     REAL, ALLOCATABLE :: P(:,:),PN(:,:),ERROR(:,:)
     INTEGER :: NY,ITA,I,J,K
     REAL :: DX,DY,XER,RMS
      INTEGER, PARAMETER :: NX = 80
     NY=NX
     DX=2.0/(NX-1.0)
     DY=2.0/(NY-1.0)
     ALLOCATE(P(NX+1,NY+1),PN(NX+1,NY+1),ERROR(NX+1,NY+1))
    
    !INITIAL CONDITION      
     DO I=1,NX-1
       DO J=1,NY-1
          P(I,J)=0.0
       END DO
     END DO
   
    !BOUNDARY CONDITIONS
     DO J=1,NY-1
        P(1,J)=0.0
     END DO
     DO J=1,NY-1
        P(NX-1,J)= DY*J
     END DO
     DO I=1,NX-1
        P(I,1)=P(I,2)
     END DO
     DO I=1,NX-1
        P(I,NY-1)=P(I,NY-2)
     END DO


     OPEN(1,FILE='INITIAL.DAT')
      WRITE(1,*) " TITLE =  ""FLOW-FIELD"" "
      WRITE(1,*) "VARIABLES = X, Y, P"
      WRITE(1,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
        DO J=1,(NY-1)
        WRITE(1,*)I,J,P(I,J)
        END DO
     END DO
     CLOSE(1)

     CLOSE(1)
 
     ITA=0
71   ITA=ITA+1     
   
     DO I=2,NX-2
        DO J=2,NY-2
          PN(I,J)=((DY**2)*(P(I+1,J)+P(I-1,J))+(DX**2)*(P(I,J+1)+P(I,J-1)))/(2*(DX**2+DY**2))
        END DO
     END DO

     DO I=2,NX-2
        PN(I,1)=PN(I,2)
     END DO
     DO I=2,NX-2
        PN(I,NY-1)=PN(J,NY-2)
     END DO
      
     XER=0.0
     DO I=2,NX-2
        DO J=2,NY-2
          ERROR(I,J)=(PN(I,J)-P(I,J))*(PN(I,J)-P(I,J))
          XER=XER+ERROR(I,J)
        END DO
     END DO
     PRINT*,XER

     RMS=SQRT(XER/9801.0)
    
     IF(RMS.GT.1E-6) THEN
         DO I=2,NX-2
           DO J=1,NY-1
              P(I,J)=PN(I,J)
           END DO
         END DO
         GOTO 71
      END IF
      
      OPEN(2,FILE='FINAL.DAT')
      WRITE(2,*) " TITLE = ""FLOW-FIELD"" "
      WRITE(2,*)  "VARIABLES = X, Y, P"
      WRITE(2,*) "ZONE T=  ""N= 0"", I= ", NX-1, ", J= ", NY-1, ",F=POINT"
     DO I=1,(NX-1)
       DO J=1,(NY-1)
        WRITE(2,*)I,J,P(I,J)
       END DO
     END DO
     PRINT*, " FILES FOR INITIAL AND FINAL DATA ARE GENERATED."
     PRINT*, " USE PARAVIEW OR TECPLOT OR ANY OTHER TOOL TO VISUALIZE."
             
     CLOSE(2)
     END PROGRAM TWOD_LAPLACE



    
          




      

