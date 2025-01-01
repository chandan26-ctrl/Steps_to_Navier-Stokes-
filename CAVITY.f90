!****************************************************************************                                                                         *
! *   FILE         = CAVITY.F90                                             *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *                                                                        *
!**************************************************************************** 
PROGRAM CAVITY
    IMPLICIT NONE
    INTEGER, PARAMETER :: NX = 101, NY = 101, NT = 700, NIT = 50
    REAL(8), PARAMETER :: XMIN = 0.D0, YMIN = 0.D0, XMAX = 1.D0, YMAX = 1.D0
    REAL(8)            :: DX, DY
    REAL(8), DIMENSION(NX,NY) :: U, UN, V, VN
    REAL(8), DIMENSION(NX,NY) :: P, PN, B
    REAL(8), DIMENSION(NX) :: X
    REAL(8), DIMENSION(NY) :: Y
    REAL(8), PARAMETER :: RHO = 1.D0, VIS = 0.1D0, DT = 0.0001D0
    INTEGER :: I, J, K

    DX = (XMAX - XMIN) / (NX-1)
    DY = (YMAX - YMIN) / (NY-1)

    CALL  GRIDGEN(X, Y, DX, DY, NX, NY)
   

    U = 0.D0
    V = 0.D0
    P = 0.D0
    B = 0.D0

    DO K = 1, NT   ! LOOP FOR TIME
        UN = U
        VN = V

        CALL BUILDUPB(B, RHO, DT, U, V, DX, DY, NX, NY)
        CALL PRESPOISSON(P, DX, DY, B, NX, NY, NIT, RHO)

        CALL NAVIERSTOKES(U, UN, V, VN, P, NX, NY, DX, DY, DT, RHO, VIS)
        CALL BOUNDARYCONDITIONS(U, V, NX, NY)

    ENDDO

    CALL WRITEFILE(U, V, P, X, Y, NX, NY)


    CONTAINS


    SUBROUTINE BUILDUPB(B, RHO, DT, U, V, DX, DY, NX, NY)
        INTEGER, INTENT(IN) :: NX, NY
        REAL(8), DIMENSION(NX,NY), INTENT(INOUT) :: B
        REAL(8), DIMENSION(NX,NY), INTENT(IN) :: U, V
        REAL(8), INTENT(IN) :: RHO, DT, DX, DY
        INTEGER :: I, J
        
        DO J = 2, NY-1
            DO I = 2, NX-1
               B(I,J) = ( (U(I+1,J)-U(I-1,J))/2.D0/DX + &
               & (V(I,J+1)-V(I,J-1))/2.D0/DY ) / DT - &
               & ( (U(I+1,J)-U(I-1,J))/2.D0/DX )**2 - &
               & 2.D0*(U(I,J+1)-U(I,J-1))/2.D0/DY * (V(I+1,J)-V(I,J-1))/2.D0/DX - &
               & ( (V(I,J+1)-V(I,J-1))/2.D0/DY)**2
            ENDDO
        ENDDO


    END SUBROUTINE BUILDUPB


    SUBROUTINE PRESPOISSON(P, DX, DY, B, NX, NY, NIT, RHO)
        INTEGER, INTENT(IN) :: NX, NY, NIT
        REAL(8), DIMENSION(NX, NY), INTENT(INOUT) :: P
        REAL(8), DIMENSION(NX, NY) :: PN
        REAL(8), DIMENSION(NX,NY), INTENT(IN) :: B
        REAL(8), INTENT(IN) :: DX, DY, RHO
        INTEGER :: I, J, K

        PN = P
!        DO K = 1, NIT
!            PN = P
            DO J = 2, NY-1
                DO I = 2, NX-1 
                    P(I,J) = ((PN(I+1,J) + PN(I-1,J))*DY**2 + &
                           & (PN(I,J+1) + PN(I,J-1))*DX**2 - &
                           & B(I,J)*RHO* DX**2 * DY**2) / (DX**2 + DY**2) / 2.D0
                ENDDO
             ENDDO
!        ENDDO

        DO J = 1, NY
            P(NX,J) = P(NX-1,J)   ! DP/DY = 0 AT X = 2
            P(1,J) = P(2,J)       ! DP/DY = 0 AT X = 0
        ENDDO

        DO I = 1, NX
            P(I,1) = P(I,2)       ! DP/DX = 0 AT Y = 0
            P(I,NY) = 0.D0        ! P = 0 AT Y = 2
        ENDDO

    END SUBROUTINE PRESPOISSON


    SUBROUTINE GRIDGEN(X, Y, DX, DY, NX, NY)
        INTEGER, INTENT(IN) :: NX, NY
        REAL(8), DIMENSION(NX), INTENT(INOUT) :: X
        REAL(8), DIMENSION(NY), INTENT(INOUT) :: Y
        REAL(8), INTENT(IN) :: DX, DY
 
        X(1) = 0.D0
        DO I = 2, NX
            X(I) = X(I-1) + DX
        END DO

        Y(1) = 0.D0
        DO I = 2, NY
            Y(I) = Y(I-1) + DY
        END DO
    END SUBROUTINE GRIDGEN


    SUBROUTINE NAVIERSTOKES(U, UN, V, VN, P, NX, NY, DX, DY, DT, RHO, VIS)
        INTEGER, INTENT(IN) :: NX, NY
        REAL(8), DIMENSION(NX,NY), INTENT(INOUT) :: U, V
        REAL(8), DIMENSION(NX,NY), INTENT(IN) :: UN, VN, P
        REAL(8), INTENT(IN) :: DX, DY, DT
        REAL(8), INTENT(IN) :: RHO, VIS

        INTEGER :: I, J

        DO J = 2, NY-1
            DO I = 2, NX-1
                U(I,J) = UN(I,J) - UN(I,J)*DT/DX*( UN(I,J)-UN(I-1,J) ) - &
                       & VN(I,J)*DT/DY*( UN(I,J)-UN(I,J-1) ) - &
                       & 1.D0/RHO*( P(I+1,J)-P(I-1,J) )*DT/2.D0/DX + &
                       & VIS*DT/DX**2*( UN(I+1,J)-2.D0*UN(I,J)+UN(I-1,J) ) + &
                       & VIS*DT/DY**2*( UN(I,J+1)-2.D0*UN(I,J)+UN(I,J-1) )

                V(I,J) = VN(I,J) - UN(I,J)*DT/DX*( VN(I,J)-VN(I-1,J) ) - &
                       & VN(I,J)*DT/DY*( VN(I,J)-VN(I,J-1) ) - &
                       & 1.D0/RHO*( P(I,J+1)-P(I,J-1) )*DT/2.D0/DY + &
                       & VIS*DT/DX**2*( VN(I+1,J)-2.D0*VN(I,J)+VN(I-1,J) ) + &
                       & VIS*DT/DY**2*( VN(I,J+1)-2.D0*VN(I,J)+VN(I,J-1) )
            ENDDO
        ENDDO 
    END SUBROUTINE NAVIERSTOKES


    SUBROUTINE BOUNDARYCONDITIONS(U, V, NX, NY)
        INTEGER, INTENT(IN) :: NX, NY
        REAL(8), DIMENSION(NX,NY), INTENT(INOUT) :: U, V
        INTEGER :: I, J

        DO I = 1, NX
            U(I,1) = 0.D0
            U(I,NY) = 1.D0

            V(I,1) = 0.D0
            V(I,NY) = 0.D0
        ENDDO

        DO J = 1, NY
           U(1,J) = 0.D0
           U(NX,J) = 0.D0

           V(1,J) = 0.D0
           V(NX,J) = 0.D0
        ENDDO

    END SUBROUTINE BOUNDARYCONDITIONS


    SUBROUTINE WRITEFILE(U, V, P, X, Y, NX, NY)
        INTEGER, INTENT(IN) :: NX, NY
        REAL(8), DIMENSION(NX,NY), INTENT(IN) :: U, V, P
        REAL(8), DIMENSION(NX), INTENT(IN) :: X
        REAL(8), DIMENSION(NY), INTENT(IN) :: Y



      OPEN(7, FILE='UVP.dat')
      WRITE(7,*) " TITLE = ""FLOW-FIELD_UVP"" "
      WRITE(7,*)  "VARIABLES = X, Y, U, V, P"
      WRITE(7,*) "ZONE T=  ""N= 0"", I= ", NX, ", J= ", NY, ",F=POINT"
        
        DO J = 1, NY
            DO I = 1, NX
                WRITE(7, *), X(I), Y(J), U(I,J), V(I,J), P(I,J)
            END DO
        END DO
        CLOSE(7)


    END SUBROUTINE WRITEFILE

END PROGRAM CAVITY 

