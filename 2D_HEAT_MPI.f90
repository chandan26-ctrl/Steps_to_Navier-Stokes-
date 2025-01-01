!****************************************************************************
! *   FILE         = 2D_HEAT_MPI.F90                                        *
! *   DESCRIPTION  = 2D DOMAIN DECOMPOSITION (RUNS ONLY ON 4 PROCESSORS)    *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
!*                                                                          *
!*          T=100.                                                          *
!*         _________                    ________|_______      |Y            *
!*         |       |                   |        |        |    |             *
!*         |       |                   |    2   |    3   |    |             *
!* T=400.  |T=100.0| T =600.         __|________|________|__  |______X      *
!*         |       |                   |        |        |                  *
!*         |_______|                   |    0   |     1  |                  *
!*          T=1000.0                   |________|________|                  *
!*                                              |                           *
!****************************************************************************
     PROGRAM TWOD_HEAT_MPI
     IMPLICIT NONE
     INCLUDE 'mpif.h'
     INTEGER, PARAMETER :: NX=500,NY=500
     REAL, ALLOCATABLE :: T(:,:),TN(:,:), X(:), Y(:)
     REAL, ALLOCATABLE :: TSENDL(:),TSENDR(:),TRECVR(:),TRECVL(:)
     REAL, ALLOCATABLE :: TSENDT(:),TSENDB(:),TRECVT(:),TRECVB(:)
     REAL :: DX,DY,S,RMS,G_RMS, X0,Y0, XMAX, YMAX
     INTEGER :: I,J,K,M
     INTEGER :: NUMTASKS,MYID,IERR,NUMWORKER
     INTEGER :: NMAX, NSTART,NEND,NREM,ICT,ITA,NYMAX
     INTEGER :: NXREM,NYREM,DEST,MM
     DOUBLE PRECISION :: T1, T2
     INTEGER :: STATUS(MPI_STATUS_SIZE)

     CHARACTER(LEN=20) :: FILE_ID
     CHARACTER(LEN=50) :: FILE_NAME

     CALL MPI_INIT(IERR)
     CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMTASKS, IERR)
     NUMWORKER=NUMTASKS/2
     
      DX=5.0/(NX-1-1)
      DY=5.0/(NY-1-1)
     
      NMAX=NX/NUMWORKER
      NYMAX=NY/NUMWORKER

      NXREM=NX-NMAX*NUMWORKER
      IF(MYID.LT.NXREM) THEN
       NMAX=NMAX+1
      END IF
      NYREM=NX-NYMAX*NUMWORKER
      IF(MYID.LT.NYREM) THEN
         NYMAX=NYMAX+1
      END IF
      IF(MYID.EQ.0) THEN 
       DO DEST=1,NUMWORKER-1
            CALL MPI_SEND(NYMAX,1,MPI_INTEGER,DEST,10, MPI_COMM_WORLD,IERR)
       END DO
      END IF

      IF (MYID.LT.NUMWORKER.AND.MYID.GT.0) THEN
           CALL MPI_RECV(NYMAX,1,MPI_INTEGER,0,10,MPI_COMM_WORLD,STATUS,IERR)
      ENDIF
    

     

      IF(MYID.EQ.0) THEN 
       DO DEST=NUMWORKER,NUMTASKS-NUMWORKER,NUMWORKER
           CALL MPI_SEND(NMAX,1,MPI_INTEGER,DEST,20,MPI_COMM_WORLD,IERR)
       END DO
      END IF

      IF(MOD(MYID,NUMWORKER).EQ.0 .AND. MYID.NE.0) THEN
          CALL MPI_RECV(NMAX,1, MPI_INTEGER,0,20,MPI_COMM_WORLD,STATUS,IERR)
      END IF


      IF (MYID.GT.0 .AND. MYID.LT.NUMWORKER) THEN
           MM=NUMWORKER+MYID
         DO DEST=1,NUMWORKER-1
          CALL MPI_SEND(NMAX,1,MPI_INTEGER,MM,30,MPI_COMM_WORLD,IERR)
           MM=MM+NUMWORKER
         END DO
      END IF

      IF(MYID.GT.NUMWORKER.AND. MOD(MYID,NUMWORKER).NE.0) THEN
          IF(MOD((MYID-1),NUMWORKER).EQ.0 .AND.(MYID-1).NE.0) THEN
            CALL MPI_RECV(NMAX,1,MPI_INTEGER,1,30,MPI_COMM_WORLD,STATUS,IERR)
          END IF
      END IF
          
          
          
         NMAX=NMAX+1
         NYMAX=NYMAX+1  
        
       PRINT*,'PROCESSOR ID:',MYID,NMAX,NYMAX
   
        ALLOCATE(T(NMAX,NYMAX))
        ALLOCATE(TN(NMAX,NYMAX))
        ALLOCATE(TSENDL(NYMAX),TRECVL(NYMAX),TSENDR(NYMAX),TRECVR(NYMAX))
        ALLOCATE(TSENDT(NMAX), TRECVT(NMAX), TSENDB(NMAX), TRECVB(NMAX))
        ALLOCATE(X(NMAX), Y(NY))

!  GENERATE GRID
         X(1) = 0.0
         Y(1) = 0.0
       DO I =2, NMAX 
           X(I)  = X(I-1)+DX
       END DO
       DO I = 2, NYMAX
          Y(I) = Y(I-1) +DY
       END DO
      

       IF (MYID .EQ. 0) THEN
        DO I =1,3
          XMAX = X(NMAX)
          YMAX = Y(NYMAX)
          CALL MPI_SEND ( XMAX, 1, MPI_REAL, I,10, MPI_COMM_WORLD, IERR)
          CALL MPI_SEND ( YMAX, 1, MPI_REAL, I,10, MPI_COMM_WORLD, IERR)
        END DO
       ELSE
         CALL MPI_RECV(X0, 1, MPI_REAL, 0,10, MPI_COMM_WORLD, STATUS, IERR)
         CALL MPI_RECV(Y0, 1, MPI_REAL, 0,10, MPI_COMM_WORLD, STATUS, IERR)
         IF (MYID .EQ. NUMTASKS-1) THEN
            X(1) = X0- 2.0*DX
            Y(1) = Y0 - 2.0*DY
          DO I =2, NMAX
             X(I) =X(I-1) +DX
          END DO
          DO I =2, NYMAX
             Y(I) =Y(I-1) +DY
          END DO
         ELSE IF (MYID .EQ. 1) THEN
            X(1) = X0- 2.0*DX
            DO I =2, NMAX
             X(I) =X(I-1) +DX
          END DO
         ELSE
              Y(1) = Y0 - 2.0*DY
          DO I =2, NYMAX
             Y(I) =Y(I-1) +DY
          END DO
        END IF
      END IF


 !INITIAL CONDITION
      DO I=1,NMAX
       DO J=1,NYMAX
          T(I,J)=100.0
       END DO
      END DO

      IF (MYID.LT.NUMWORKER) THEN
        DO I=1,NMAX
          T(I,1)=1000.0
        END DO
        IF(MYID.EQ.0) THEN
          DO J=1,NYMAX
           T(1,J)=400.0
          END DO
        ELSE
         DO J=1,NYMAX
          T(NMAX,J)=600.0
         END DO
        END IF
      END IF

      IF (MYID .GE.NUMWORKER) THEN
         DO I=1,NMAX
            T(I,NYMAX)=100.0
         END DO
        IF (MYID.EQ.NUMWORKER) THEN
           DO J=1,NYMAX
             T(1,J)=400.0
           END DO
        ELSE
          DO J=1,NYMAX
             T(NMAX,J)=600.0
          END DO
        END IF
     END IF

      T1 = MPI_Wtime()
      ITA=0
        
71    ITA=ITA+1
      DO J=2,NYMAX-1
      DO I=2,NMAX-1
         TN(I,J) = (DY**2*(T(I+1,J)+T(I-1,J)) + DX**2*(T(I,J+1)+T(I,J-1)))/(2.0*DY**2+2.0*DX**2)
      END DO
      END DO
    
      IF(MOD(MYID,NUMWORKER).EQ.0) THEN
         ICT=1
        DO J=2,NYMAX-1
           TSENDR(ICT)=TN(NMAX-1,J)
           ICT=ICT+1
        END DO
      ELSE
        ICT=1
        DO J=2,NYMAX-1
          TSENDL(ICT)=TN(2,J)
          ICT=ICT+1
        END DO
      END IF

       

      IF(MYID.LT.NUMWORKER) THEN
        ICT=1
       DO I=2,NMAX-1
         TSENDT(ICT)=TN(I,NYMAX-1)
         ICT=ICT+1
       END DO
      ELSE
       ICT=1
       DO I=2,NMAX-1
         TSENDB(ICT)=TN(I,2)
         ICT=ICT+1
       END DO
      END IF

      IF(MOD(MYID,NUMWORKER).EQ.0) THEN
        CALL MPI_SEND(TSENDR,NYMAX-2,MPI_REAL,MYID+1,10,MPI_COMM_WORLD,IERR)
      ELSE
        CALL MPI_SEND(TSENDL,NYMAX-2,MPI_REAL,MYID-1,20,MPI_COMM_WORLD,IERR)
      END IF
      

      IF(MOD(MYID,NUMWORKER).EQ.0) THEN
          CALL MPI_RECV(TRECVL,NYMAX-2,MPI_REAL,MYID+1,20,MPI_COMM_WORLD,STATUS,IERR)
      ELSE
          CALL MPI_RECV(TRECVR,NYMAX-2,MPI_REAL,MYID-1,10,MPI_COMM_WORLD,STATUS,IERR)
      END IF
        
      IF(MYID.LT.NUMWORKER) THEN
        CALL MPI_SEND(TSENDT,NMAX-2,MPI_REAL,MYID+NUMWORKER,30,MPI_COMM_WORLD,IERR)
      ELSE
        CALL MPI_SEND(TSENDB,NMAX-2,MPI_REAL,MYID-NUMWORKER,40,MPI_COMM_WORLD,IERR)
      END IF

      IF(MYID.LT.NUMWORKER) THEN
         CALL MPI_RECV(TRECVB,NMAX-2,MPI_REAL,MYID+NUMWORKER,40,MPI_COMM_WORLD,STATUS,IERR)
      ELSE
         CALL MPI_RECV(TRECVT,NMAX-2,MPI_REAL,MYID-NUMWORKER,30,MPI_COMM_WORLD,STATUS,IERR)
      END IF

      IF(MOD(MYID,NUMWORKER).EQ.0) THEN
         ICT=1
        DO J=2,NYMAX-1
         T(NMAX,J)=TRECVL(ICT)
        ICT=ICT+1
        END DO
      ELSE
         ICT=1
        DO J=2,NYMAX-1
           T(1,J)=TRECVR(ICT)
         ICT=ICT+1
        END DO
      END IF
   
      IF(MYID.LT.NUMWORKER) THEN
         ICT=1
       DO I=2,NMAX-1
          T(I,NYMAX)=TRECVB(ICT)
          ICT=ICT+1
       END DO
      ELSE
        ICT=1
       DO I=2,NMAX-1
          T(I,1)=TRECVT(ICT)
          ICT=ICT+1
       END DO
      END IF
  
   S=0.0
   DO I=2,NMAX-1
      DO J=2,NYMAX-1
      S=S+(TN(I,J)-T(I,J))**2
      END DO
   END DO
    S=S/((NMAX-2)*(NY-2))
    RMS=SQRT(S)

    DO I=2,NMAX-1
     DO J=2,NYMAX-1
       T(I,J)=TN(I,J)
     END DO
    END DO

   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(RMS,G_RMS,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,IERR)

   IF(MYID.EQ.0) THEN
    PRINT*, 'ERROR=', G_RMS,'ITERATION',ITA
   END IF

   IF(G_RMS.GT.1.0E-5) GOTO 71
   IF(MYID.EQ.0) PRINT*,'CONVERGED'
    T2 = MPI_Wtime()
   IF (MYID .EQ. 0) THEN
    PRINT*, 'ELAPSED TIME:', T2-T1     
   END IF


! WRITE FILES FOR TECPLOT FROM EACH PROCESSOR
   WRITE(FILE_ID, '(I6)') MYID+1
   FILE_NAME = 'PRINT' // TRIM(ADJUSTL(FILE_ID)) // '.dat'

   OPEN(FILE=TRIM(FILE_NAME), UNIT = 4)
   

   WRITE(4,*) " TITLE = ""FLOW-FIELD 2D"" "
   WRITE(4,*)  "VARIABLES =  X, Y, T"
   
   IF(MYID.GT.0 .AND. MYID.LT.NUMTASKS-1)THEN
   WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NYMAX, ",F=POINT"
        DO J=1,NYMAX 
        DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
       END DO
   ELSE 
     IF(MYID.EQ.0) THEN
      WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NYMAX, ",F=POINT"
         DO J=1,NYMAX 
          DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
         END DO
        END DO
      ENDIF
      IF(MYID.EQ.NUMTASKS-1) THEN
       WRITE(4,*) "ZONE T=  ""N= 0"", I= ", NMAX, ", J= ", NYMAX, ",F=POINT"
         DO J=1,NYMAX
          DO I=1,NMAX
           WRITE(4,*)X(I),Y(J),T(I,J)
          END DO
         END DO
      END IF
    END IF
    CLOSE (4)
 
   CALL MPI_FINALIZE(IERR)

   
   END PROGRAM  

