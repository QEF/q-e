MODULE wyckoff
USE kinds,  ONLY : DP
USE space_group, ONLY : sym_brav, find_equivalent_tau
IMPLICIT NONE

INTEGER :: nattot
REAL(DP), ALLOCATABLE :: tautot(:,:), extfortot(:,:)
INTEGER, ALLOCATABLE :: ityptot(:), if_postot(:,:)

SAVE
PRIVATE

PUBLIC sup_spacegroup, clean_spacegroup, nattot, tautot, ityptot, extfortot, &
       if_postot

CONTAINS

   SUBROUTINE sup_spacegroup(tau,ityp,extfor,if_pos,space_group_number,not_eq,&
              uniqueb,rhombohedral,choice,ibrav)
      INTEGER, INTENT(IN) :: space_group_number, choice
      LOGICAL, INTENT (IN) :: uniqueb, rhombohedral
      INTEGER, INTENT (INOUT) ::  not_eq
      INTEGER, INTENT(OUT) :: ibrav
      REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: tau, extfor
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ityp
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: if_pos
      INTEGER :: i,k,l,n,sym_n
      INTEGER,DIMENSION(:),allocatable :: msym_n 
      character(LEN=1) :: unique
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: inco
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: outco   


      ALLOCATE(inco(3,not_eq))
      ALLOCATE(msym_n(not_eq))
     
     inco=tau

     !conv from uniqueb,rhombohedral,choice to unique
         unique='1'
        IF ((uniqueb).or.(.not.rhombohedral).or.(choice==2)) then
            unique='2'
        END IF
     !select ibrav and number of symmetries
     CALL sym_brav(space_group_number,sym_n,ibrav)
     
     do i=1,not_eq
        msym_n(i)=sym_n
     end do
  
     IF (((ibrav==12).or.(ibrav==13)).and.(unique=='2')) ibrav=-ibrav
     
     ALLOCATE(outco(3,sym_n,not_eq))
     
     !make symmetries, convert coordinates, esclusion
     DO i=1,not_eq
        CALL find_equivalent_tau(space_group_number, inco, outco,i,unique)
     END DO

     call ccord(outco,sym_n,not_eq,ibrav,unique)
     do i=1,not_eq
       call zerone(outco,sym_n,not_eq,i)
     end do
     call esclusion(outco,sym_n,msym_n,not_eq)

     nattot=SUM(msym_n)
     
     ALLOCATE(tautot(3,nattot))
     ALLOCATE(ityptot(nattot))
     ALLOCATE(extfortot(3,nattot))
     ALLOCATE(if_postot(3,nattot))

     !conversione tra outco e tau
     l=0
     DO i=1,not_eq
        IF (i/=1) THEN
           l=l+msym_n(i-1)
        END IF
        !
        DO k=1,msym_n(i)
           tautot(1,k+l)=outco(1,k,i)
           tautot(2,k+l)=outco(2,k,i)
           tautot(3,k+l)=outco(3,k,i)
           ityptot(k+l) = ityp(i)
           if_postot(:,k+l) = if_pos(:,i)
           extfortot(:,k+l) = extfor(:,i)
        END DO
     END DO

     DEALLOCATE(inco)
     DEALLOCATE(outco)
     DEALLOCATE(msym_n)

     RETURN
END SUBROUTINE sup_spacegroup

SUBROUTINE clean_spacegroup

DEALLOCATE(tautot)
DEALLOCATE(ityptot)
DEALLOCATE(extfortot)
DEALLOCATE(if_postot)

RETURN
END SUBROUTINE clean_spacegroup

   SUBROUTINE ccord(outco,sym_n,not_eq,ibrav,unique)
   IMPLICIT NONE
      REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: outco
      INTEGER, INTENT(in) :: ibrav,sym_n,not_eq
      CHARACTER, INTENT(in) :: unique
      INTEGER :: i,k
      REAL(DP) :: tmpx, tmpy, tmpz

      Cambio: SELECT CASE (ibrav)
      CASE (2) !fcc
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)
            
            outco(1,i,k)=-tmpx-tmpy+tmpz
            outco(2,i,k)=tmpx+tmpy+tmpz
            outco(3,i,k)=-tmpx-tmpz+tmpy
         END DO
      END DO

      CASE (3) !bcc
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)
            
            outco(1,i,k)=tmpx+tmpz
            outco(2,i,k)=tmpy-tmpx
            outco(3,i,k)=tmpz-tmpy
         END DO
      END DO

      CASE (5) !Only for trigonal
         IF (unique=='2') THEN
            DO k=1,not_eq
               DO i=1,sym_n
                  tmpx=outco(1,i,k)
                  tmpy=outco(2,i,k)
                  tmpz=outco(3,i,k)
         
                  outco(1,i,k)=tmpx-tmpy+tmpz
                  outco(2,i,k)=tmpy+tmpz
                  outco(3,i,k)=tmpz-tmpx
               END DO
            END DO
         END IF   

      CASE (7) !Body Centred Tetragonal
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)
            
            outco(1,i,k)=tmpx-tmpy
            outco(2,i,k)=tmpy+tmpz
            outco(3,i,k)=tmpz-tmpx
         END DO
      END DO

      CASE (9) !Base Centrata ORTHORHOMBIC C
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)
            
            outco(1,i,k)=tmpx+tmpy
            outco(2,i,k)=tmpy-tmpx
            outco(3,i,k)=tmpz
         END DO
      END DO

      CASE (91) !Base Centrata ORTHORHOMBIC A
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)

            outco(1,i,k)=tmpx
            outco(2,i,k)=tmpy+tmpz
            outco(3,i,k)=tmpy-tmpz
         END DO
      END DO

      CASE (10) !Tutte le faccie centrate ORTHORHOMBIC
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)

            outco(1,i,k)=tmpx-tmpy+tmpz
            outco(2,i,k)=tmpx+tmpy-tmpz
            outco(3,i,k)=-tmpx+tmpy+tmpz
         END DO
      END DO
      
      CASE (11) !Corpo Centrato ORTHORHOMBIC
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
            tmpz=outco(3,i,k)
            
            outco(1,i,k)=tmpx+tmpz
            outco(2,i,k)=tmpy-tmpx
            outco(3,i,k)=tmpz-tmpy
         END DO
      END DO

      CASE (13) !Centrato C unique MONOCLINO
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpz=outco(3,i,k)

            outco(1,i,k)=tmpx-tmpz
            outco(3,i,k)=tmpz+tmpx
         END DO
      END DO

      CASE (-13) !Centrato B unique MONOCLINO
      DO k=1,not_eq
         DO i=1,sym_n
            tmpx=outco(1,i,k)
            tmpy=outco(2,i,k)
      
            outco(1,i,k)=tmpx-tmpy
            outco(2,i,k)=tmpy+tmpx
         END DO
      END DO
      
      END SELECT cambio
   END SUBROUTINE ccord

!Translation in order to have 0<=x,y,z<=1
   SUBROUTINE zerone(outco,sym_n, not_eq,k)
   IMPLICIT NONE
      REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: outco
      INTEGER, INTENT(in) :: sym_n,not_eq, k
      INTEGER :: i
         DO i=1, sym_n
            DO
            IF (outco(1,i,k)>=1.0_DP) THEN
               outco(1,i,k)=outco(1,i,k)-1.0_DP
            else
               EXIT
            END IF
            END DO

            DO
            IF (outco(2,i,k)>=1.0_DP) THEN
               outco(2,i,k)=outco(2,i,k)-1.0_DP
            else
               EXIT
            END IF
            END DO

            DO
            IF (outco(3,i,k)>=1.0_DP) THEN
               outco(3,i,k)=outco(3,i,k)-1.0_DP
            else
               EXIT
            END IF
            END DO

            DO
            IF (outco(1,i,k)<0.0_DP) THEN
               outco(1,i,k)=outco(1,i,k)+1.0_DP
            else
               EXIT
            END IF
            END DO

            DO
            IF (outco(2,i,k)<0.0_DP) THEN
               outco(2,i,k)=outco(2,i,k)+1.0_DP
            else
               EXIT
            END IF
            END DO

            DO
            IF (outco(3,i,k)<0.0_DP) THEN
               outco(3,i,k)=outco(3,i,k)+1.0_DP
            else
               EXIT
            END IF
            END DO
         END DO
         RETURN
   END SUBROUTINE zerone

   SUBROUTINE esclusion(outco,sym_n,msym_n,not_eq)
      IMPLICIT NONE
      REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: outco
      INTEGER, DIMENSION (:), INTENT(INOUT) :: msym_n
      INTEGER, INTENT(in) :: not_eq, sym_n
      INTEGER :: i,l,k,j
      REAL(DP), DIMENSION(:,:,:),allocatable :: temp
      LOGICAL :: bol
      REAL :: eps

      eps=1.D-6

      ALLOCATE(temp(3,sym_n,not_eq))

      DO k=1,not_eq
         l=0
         DO j=1,sym_n
            bol=.false.
            i=j+1
            DO while (i<=sym_n)
               IF ((abs(outco(1,j,k)-outco(1,i,k))<eps).and.&
                   (abs(outco(2,j,k)-outco(2,i,k))<eps).and.&
                   (abs(outco(3,j,k)-outco(3,i,k))<eps)) THEN
                  bol=.true.
               END IF
            i=i+1
            END DO
         
            IF (.not.bol) THEN
               l=l+1
               temp(1,l,k)=outco(1,j,k)
               temp(2,l,k)=outco(2,j,k)
               temp(3,l,k)=outco(3,j,k)
            END IF
         END DO
         msym_n(k)=l
      END DO

      outco=temp

      DEALLOCATE(temp)

      RETURN
END SUBROUTINE esclusion

END MODULE wyckoff
