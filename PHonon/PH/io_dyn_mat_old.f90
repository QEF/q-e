!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
Module dynamicalq
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  COMPLEX(DP), ALLOCATABLE :: phiq(:,:,:,:,:)
  REAL(DP), ALLOCATABLE ::  tau(:,:), zeu(:,:,:)
  INTEGER, ALLOCATABLE ::  ityp(:)
  !
end module
!-----------------------------------------------------------------------
subroutine write_dyn_on_file (xq, phi, nat, iudyn)
  !-----------------------------------------------------------------------
  USE kinds, only : DP
  implicit none
  ! input variables
  integer :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell
  complex(DP) :: phi (3, 3, nat, nat)
  !  the dynamical matrix
  real(DP) :: xq (3)
  ! the q vector
  ! local variables

  integer :: na, nb, icar, jcar
  ! counters on atoms
  ! cartesian coordinate counters
  write (iudyn, 9000) (xq (icar), icar = 1, 3)
  do na = 1, nat
     do nb = 1, nat
        write (iudyn, '(2i5)') na, nb
        do icar = 1, 3
!           write (iudyn, '(3e24.12)') (phi(icar,jcar,na,nb), jcar=1,3)
           write (iudyn, '(3(2f12.8,2x))') (phi(icar,jcar,na,nb), jcar=1,3)
        enddo
     enddo
  enddo

  return
9000 format(/,5x,'Dynamical  Matrix in cartesian axes', &
       &       //,5x,'q = ( ',3f14.9,' ) ',/)
end subroutine write_dyn_on_file


  SUBROUTINE write_old_dyn_mat_head(iudyn)
!
!  This routine is here for compatibility with the old code.
!  It will be removed when the xml file format of the dynamical matrix
!  will be tested.
!
  USE constants, ONLY: amu_ry
  USE ions_base, ONLY : ntyp => nsp, nat, ityp, tau, atm, amass
  USE cell_base, ONLY : ibrav, celldm, at
  USE run_info, ONLY : title

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iudyn
  INTEGER :: nt, na, i, j

  WRITE (iudyn, '("Dynamical matrix file")')
  WRITE (iudyn, '(a)') title
  WRITE (iudyn, '(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  IF (ibrav==0) THEN
     WRITE (iudyn,'("Basis vectors")')
     WRITE (iudyn,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  END IF
  DO nt = 1, ntyp
     WRITE (iudyn, * ) nt, ' ''', atm (nt) , ' '' ', amu_ry*amass(nt)
  ENDDO
  DO na = 1, nat
     WRITE (iudyn, '(2i5,3f18.10)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  ENDDO

  RETURN
  END SUBROUTINE write_old_dyn_mat_head
!

!----------------------------------------------------------------------------
SUBROUTINE read_dyn_from_file( nqs, xq, epsil, lrigid, &
                               ntyp, nat, ibrav, celldm, at, atm, amass )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE dynamicalq, ONLY: phiq, tau, ityp, zeu
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps8=1.D-8
  ! I/O variables
  LOGICAL :: lrigid
  INTEGER :: nqs, ntyp, nat, ibrav
  REAL(DP) :: epsil(3,3)
  REAL(DP) :: xq(3,48), celldm(6), at(3,3), amass(ntyp)
  CHARACTER(LEN=3) atm(ntyp)
  ! local variables
  INTEGER :: ntyp1,nat1,ibrav1,ityp1
  INTEGER :: i, j, na, nb, nt, ios
  REAL(DP) :: tau1(3), amass1, at1(3,3), celldm1(6), q2
  REAL(DP) :: phir(3),phii(3)
  CHARACTER(LEN=75) :: line
  CHARACTER(LEN=3)  :: atm1
  LOGICAL, SAVE :: first =.TRUE.
  !
  IF (ionode) THEN
     READ(1,*)
     READ(1,*)
  ENDIF
  IF (first) THEN
     !
     ! read cell information from file
     !
     IF (ionode) THEN
        READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
        if (ibrav==0) then
           read (1,'(a)') atm1  ! for compatibility
           read (1,*) ((at(i,j),i=1,3),j=1,3)
        end if
     END IF
     CALL mp_bcast(ntyp, ionode_id, intra_image_comm)
     CALL mp_bcast(nat, ionode_id, intra_image_comm)
     CALL mp_bcast(ibrav, ionode_id, intra_image_comm)
     CALL mp_bcast(celldm, ionode_id, intra_image_comm)
     IF (ibrav==0) THEN
        CALL mp_bcast(at, ionode_id, intra_image_comm)
     ENDIF

     IF (ntyp.GT.nat) CALL errore('read_dyn_from_file','ntyp.gt.nat!!',ntyp)
     DO nt = 1,ntyp
        IF (ionode) READ(1,*) i,atm(nt),amass(nt)
        CALL mp_bcast(i, ionode_id, intra_image_comm)
        IF (i.NE.nt) CALL errore('read_dyn_from_file','wrong data read',nt)
     END DO
     CALL mp_bcast(atm, ionode_id, intra_image_comm)
     CALL mp_bcast(amass, ionode_id, intra_image_comm)
     ALLOCATE ( ityp(nat), tau(3,nat) )
     DO na=1,nat
        IF (ionode) READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
        CALL mp_bcast(i, ionode_id, intra_image_comm)
        IF (i.NE.na) CALL errore('read_dyn_from_file','wrong data read',na)
     END DO
     CALL mp_bcast(ityp, ionode_id, intra_image_comm)
     CALL mp_bcast(tau, ionode_id, intra_image_comm)
     !
     ALLOCATE ( phiq (3,3,nat,nat,48), zeu (3,3,nat) )
     !
     first=.FALSE.
     lrigid=.FALSE.
     !
  ELSE
     !
     ! check cell information with previous one
     !
     IF (ionode) READ(1,*) ntyp1,nat1,ibrav1,(celldm1(i),i=1,6)
     CALL mp_bcast(ntyp1, ionode_id, intra_image_comm)
     CALL mp_bcast(nat1, ionode_id, intra_image_comm)
     CALL mp_bcast(ibrav1, ionode_id, intra_image_comm)
     CALL mp_bcast(celldm1, ionode_id, intra_image_comm)
     IF (ntyp1.NE.ntyp) CALL errore('read_dyn_from_file','wrong ntyp',1)
     IF (nat1.NE.nat) CALL errore('read_dyn_from_file','wrong nat',1)
     IF (ibrav1.NE.ibrav) CALL errore('read_dyn_from_file','wrong ibrav',1)
     DO i=1,6
        IF( abs (celldm1(i)-celldm(i)) > eps8 ) &
             CALL errore('read_dyn_from_file','wrong celldm',i)
     END DO
     if (ibrav==0) then
         IF (ionode) read (1,'(a)') atm1 ! for compatibility
         IF (ionode) read (1,*) ((at1(i,j),i=1,3),j=1,3)
         CALL mp_bcast(at1, ionode_id, intra_image_comm)
         do i=1,3
            do j=1,3
               if( abs (at1(i,j)-at(i,j)) > eps8) &
                 CALL errore('read_dyn_from_file','wrong at(i,j)',i+3*(j-1))
            end do
         end do
     end if
     DO nt = 1,ntyp
        IF (ionode) READ(1,*) i,atm1,amass1
        CALL mp_bcast(i, ionode_id, intra_image_comm)
        CALL mp_bcast(atm1, ionode_id, intra_image_comm)
        CALL mp_bcast(amass1, ionode_id, intra_image_comm)
        IF (i.NE.nt) CALL errore('read_dyn_from_file','wrong data read',nt)
        IF (atm1.NE.atm(nt)) CALL errore('read_dyn_from_file','wrong atm',nt)
        IF (abs(amass1-amass(nt)) > eps8 ) &
             CALL errore('read_dyn_from_file','wrong amass',nt)
     END DO
     DO na=1,nat
        IF (ionode) READ(1,*) i,ityp1,(tau1(j),j=1,3)
        CALL mp_bcast(i, ionode_id, intra_image_comm)
        CALL mp_bcast(ityp1, ionode_id, intra_image_comm)
        CALL mp_bcast(tau1, ionode_id, intra_image_comm)
        IF (i.NE.na) CALL errore('read_dyn_from_file','wrong data read',na)
        IF (ityp1.NE.ityp(na)) CALL errore('read_dyn_from_file','wrong ityp',na)
        IF ( abs (tau1(1)-tau(1,na)) > eps8 .OR. &
             abs (tau1(2)-tau(2,na)) > eps8 .OR. &
             abs (tau1(3)-tau(3,na)) > eps8 ) &
             CALL errore('read_dyn_from_file','wrong tau',na)
     END DO
  END IF
  !
  !
  nqs = 0
100 CONTINUE
  IF (ionode) THEN
     READ(1,*,iostat=ios)
     IF(ios==0) READ(1,'(a)',iostat=ios) line
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF(ios==0) CALL mp_bcast(line, ionode_id, intra_image_comm)
  !
  IF (ios/=0 .or. line(6:14).NE.'Dynamical') THEN
     IF (nqs.EQ.0) CALL errore('read_dyn_from_file',' stop with nqs=0 !!',1)
     q2 = xq(1,nqs)**2 + xq(2,nqs)**2 + xq(3,nqs)**2
     IF (q2.NE.0.d0) RETURN
     DO WHILE (line(6:15).NE.'Dielectric')
        IF (ionode) READ(1,'(a)',iostat=ios) line
        CALL mp_bcast(ios, ionode_id, intra_image_comm)
        IF (ios /=0) GOTO 200
        CALL mp_bcast(line,ionode_id, intra_image_comm)
     END DO
     lrigid=.TRUE.
     IF (ionode) THEN
        READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
        READ(1,*)
        READ(1,*)
        READ(1,*)
     ENDIF
     CALL mp_bcast(epsil,ionode_id, intra_image_comm)
     WRITE (stdout,*) 'macroscopic fields =',lrigid
     WRITE (stdout,'(3f10.5)') ((epsil(i,j),j=1,3),i=1,3)
     IF (ionode) THEN
        DO na=1,nat
           READ(1,*)
           READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
           WRITE (stdout,*) ' na= ', na
           WRITE (stdout,'(3f10.5)') ((zeu(i,j,na),j=1,3),i=1,3)
        END DO
     END IF
     CALL mp_bcast(zeu,ionode_id, intra_image_comm)
     RETURN
200  WRITE (stdout,*) ' Dielectric Tensor not found'
     lrigid=.FALSE.
     RETURN
  END IF
  !
  nqs = nqs + 1
  IF (ionode) THEN
     READ(1,*)
     READ(1,'(a)') line
     READ(line(11:75),*) (xq(i,nqs),i=1,3)
     READ(1,*)
  ENDIF
  CALL mp_bcast(xq(:,nqs), ionode_id, intra_image_comm)
  !
  DO na=1,nat
     DO nb=1,nat
        IF (ionode) READ(1,*) i,j
        CALL mp_bcast(i, ionode_id, intra_image_comm)
        CALL mp_bcast(j, ionode_id, intra_image_comm)
        IF (i.NE.na) CALL errore('read_dyn_from_file','wrong na read',na)
        IF (j.NE.nb) CALL errore('read_dyn_from_file','wrong nb read',nb)
        DO i=1,3
           IF (ionode) READ (1,*) (phir(j),phii(j),j=1,3)
           CALL mp_bcast(phir, ionode_id, intra_image_comm)
           CALL mp_bcast(phii, ionode_id, intra_image_comm)
           DO j = 1,3
              phiq (i,j,na,nb,nqs) = CMPLX(phir(j),phii(j),kind=DP)
           END DO
        END DO
     END DO
  END DO
  !
  go to 100
  !
END SUBROUTINE read_dyn_from_file
!
