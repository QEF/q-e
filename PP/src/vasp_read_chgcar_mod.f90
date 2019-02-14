!--------------------------------------------------------------------
!
! 
! Program written by Yang Jiao, Sep 2017, GPL, No warranties.
!
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE vasp_read_chgcar
  !----------------------------------------------------------------------
  ! 
  ! ... this module contains methods to read data produced by VASP 
  ! ... tested on VASP.5.3.3
  !
  !
  USE kinds, ONLY : DP
  !
!  USE io_files,      ONLY : iunpun
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_images,     ONLY : intra_image_comm
  USE mp,            ONLY : mp_bcast, mp_sum, mp_max
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  ! internal data to be set
  ! 
  INTEGER          :: iunit !, ounit
  INTEGER          :: iunchg  = 41
  !
  ! end of declarations
  !
  PUBLIC :: vaspread_rho
  !
  !
  CONTAINS
  !
!----------------------------------------------
! ... read CHARCAR
!----------------------------------------------
    !--------------------------------------------------------------------
    SUBROUTINE vaspread_rho(ierr)
    !-----------------------------------------------------------------------
      ! This subroutine will read information from VASP output
      !    rho from CHARCAR
      !
      USE cell_base,     ONLY : omega
      USE fft_base,      ONLY : dfftp
      USE fft_interfaces, ONLY : fwfft
      USE scatter_mod,   ONLY : scatter_grid
      USE scf,           ONLY : scf_type, create_scf_type
      USE scf,           ONLY : rho, create_scf_type
      USE ions_base,     ONLY : nat, atm
      USE io_files,      ONLY : tmp_dir
      USE lsda_mod,      ONLY : nspin
      USE wavefunctions, ONLY : psic
      USE gvect,         ONLY : ngm
      !
      IMPLICIT NONE
      INTEGER,                 INTENT(out)  :: ierr
      !
      INTEGER                  :: ngxf, ngyf, ngzf, nalloc
      INTEGER                  :: ispin, iat, iz, ixy, nread
      REAL(DP), ALLOCATABLE    :: rho_r_(:,:), atomom(:)
      REAL(DP), ALLOCATABLE    :: rho_r_up(:), rho_r_dn(:)
      CHARACTER(LEN=80)        :: errmsg
      !
      ierr = 0
!      CALL create_scf_type( rho )
      ALLOCATE( rho_r_(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nspin) )
      rho_r_=0._DP
      IF ( ionode ) THEN
         WRITE(stdout,'(5X,a)') "Read charge density from CHGCAR"
         OPEN(UNIT=iunchg,FILE=TRIM(tmp_dir)//'CHGCAR',STATUS='old',ACTION='read',IOSTAT=ierr)
         CALL vaspread_rhoheader(ierr)
         DO ispin = 1, nspin
            IF(ispin==2) THEN
               ALLOCATE(atomom(nat))
               READ(iunchg,*,IOSTAT=ierr) (atomom(iat),iat=1,nat)
               DEALLOCATE(atomom)
            END IF
            READ( iunchg, '(3I5)') ngxf, ngyf, ngzf   
            IF((ngxf.NE.dfftp%nr1).OR.(ngyf.NE.dfftp%nr2).OR.(ngzf.NE.dfftp%nr3)) THEN
               errmsg = 'Dimension in CHGCAR not compatible'
            END IF
            nread=0
            nalloc=ngxf*ngyf
            DO iz=1, ngzf
               DO ixy=1, nalloc
                  nread=nread+1
                  IF(MOD(nread,5)==0) THEN
                     READ( iunchg, '(1X,E17.11)', ADVANCE='YES' ) rho_r_(nread,ispin)
                  ELSE
                     READ( iunchg, '(1X,E17.11)', ADVANCE='NO' ) rho_r_(nread,ispin)
                  END IF
               END DO
            END DO
            IF(MOD(nread,5)/=0) READ( iunchg, '(1X)' )
            CALL vaspread_rhoaugocc(ierr)
         END DO

         CLOSE(iunchg)
         IF(nspin==2) THEN 
            ALLOCATE(rho_r_up(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
            ALLOCATE(rho_r_dn(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
            rho_r_up=0.5_DP*(rho_r_(:,1)+rho_r_(:,2))
            rho_r_dn=0.5_DP*(rho_r_(:,1)-rho_r_(:,2))
            rho_r_(:,1)=rho_r_up
            rho_r_(:,2)=rho_r_dn
            DEALLOCATE(rho_r_up, rho_r_dn)
         END IF
      END IF 
!      CALL mp_bcast( atm,             ionode_id, intra_image_comm )
      DO ispin = 1, nspin
         CALL scatter_grid(dfftp, rho_r_(:,ispin), rho%of_r(:,ispin))
      END DO
      DEALLOCATE( rho_r_ )
      rho%of_r=rho%of_r/omega
!
      ! ... bring rho to G-space
      Do ispin = 1, nspin
         !
         psic(:) = rho%of_r(:,ispin)
         CALL fwfft ('Rho', psic, dfftp)
         rho%of_g(:,ispin) = psic(dfftp%nl(:))
         !
      END DO
      !
      CONTAINS
      !
      !------------------------------------------------------------------------
      SUBROUTINE vaspread_rhoheader(ierr)
        !------------------------------------------------------------------------
        USE ions_base,     ONLY : nsp
        IMPLICIT NONE
        INTEGER,                 INTENT(out)  :: ierr
        INTEGER              :: i 
        INTEGER              :: nityp(10)
        CHARACTER(LEN=3)     :: atm_(10)
        !
        ierr = 0
        READ( iunchg, * )
        DO i = 1, 4
           READ( iunchg, * )
        END DO
        READ( iunchg, '(20A5)' ) atm_(1:nsp)
        READ( iunchg, '(20I6)' ) nityp(1:nsp)
        READ( iunchg, * ) 
        DO i = 1, nat
           READ( iunchg, * )
        END DO
        READ( iunchg, * )
      END SUBROUTINE vaspread_rhoheader
      !------------------------------------------------------------------------
      SUBROUTINE vaspread_rhoaugocc(ierr)
        !------------------------------------------------------------------------
        ! read augmentation occupancies
        !
        USE ions_base,     ONLY : nat
        IMPLICIT NONE
        INTEGER,                 INTENT(out)  :: ierr
        INTEGER        :: i, iat
        INTEGER        :: ni_read, nelements_read, lmmax 
        REAL(DP), ALLOCATABLE    :: buffer_f(:)
        !
        ierr = 0
        lmmax=20   ! ???
        ALLOCATE(buffer_f(lmmax*lmmax))
        DO iat = 1, nat
           READ( iunchg, '(24X,2I4)', IOSTAT=ierr ) ni_read, nelements_read
           IF (ni_read /= iat .OR. ierr/=0 ) THEN
               ierr=1
           ENDIF
           IF(nelements_read > lmmax*lmmax) THEN
              errmsg='vaspread_rhoaugocc running out of buffer'
              STOP
           END IF
           IF(ierr==0) READ(iunchg,*,IOSTAT=ierr) (buffer_f(i),i=1,nelements_read)
        END DO

      END SUBROUTINE vaspread_rhoaugocc
    END SUBROUTINE vaspread_rho
       
END MODULE vasp_read_chgcar

    !
