!
! Copyright (C) 2014-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE YAMBO
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY :  DP
  !----------------------------------------------------------------------------
  !
  LOGICAL :: elph_yambo    = .FALSE.
  LOGICAL :: dvscf_yambo   = .FALSE.
  !
  CHARACTER(300) :: yambo_elph_file_name = " "
  !
END MODULE YAMBO
!
!-----------------------------------------------------------------------
SUBROUTINE debye_waller(y_grad_at_gamma,y_pol_vec,ibnd,jbnd,ik,mu)
  !-----------------------------------------------------------------------
  !
  ! debye_waller term
  !
  USE kinds,      ONLY:DP
  USE wvfct,      ONLY:nbnd
  USE ions_base,  ONLY:nat
  USE control_lr, ONLY:lgamma
  USE el_phon,    ONLY:el_ph_mat
  USE phcom,      ONLY:u,dyn
  !
  IMPLICIT NONE
  INTEGER     :: ibnd,jbnd,ik,mu
  COMPLEX(DP) :: y_grad_at_gamma(nat,3),y_pol_vec(nat,3)
  !
  ! Work Space
  !
  INTEGER     :: i,ia,icart,idisp
  !
  if (lgamma) then
    !
    y_grad_at_gamma(:,:)=(0._DP,0._DP)
    !
    ! At gamma The <\grad> is stored in y_grad_at_gamma
    ! after transforming from displacment basis in cartesian
    !
    do i=1,3*nat
      ia =   (i - 1) / 3 + 1
      icart = i - 3 * (ia - 1)
      do idisp=1,3*nat
         y_grad_at_gamma(ia,icart)=y_grad_at_gamma(ia,icart)+&
&                                  el_ph_mat(ibnd,jbnd,ik,idisp)*conjg(u(i,idisp))
      enddo
    enddo
  endif
  !
  ! Polarization vector (/sqrt(atom_mass))
  !
  do i=1,3*nat
    ia =   (i - 1) / 3 + 1
    icart = i - 3 * (ia - 1)
    y_pol_vec(ia,icart) = dyn(i,mu)
  enddo
  !
end subroutine debye_waller
!
!-----------------------------------------------------------------------
SUBROUTINE elph_yambo_eval_and_IO( )
  !-----------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE ions_base,   ONLY : nat
  USE wvfct,       ONLY : nbnd, et
  USE el_phon,     ONLY : el_ph_mat
  USE klist,       ONLY : xk
  USE qpoint,      ONLY : xq, nksq, ikks
  USE modes,       ONLY : u,nmodes
  USE dynmat,      ONLY : w2,dyn
  USE io_global,   ONLY : ionode, ionode_id
  USE xml_io_base, ONLY : create_directory
  USE mp_images,   ONLY : intra_image_comm
  USE mp,          ONLY : mp_bcast
  USE control_lr,  ONLY : lgamma
  USE control_ph,  ONLY : current_iq
  USE cell_base,   ONLY : alat
  USE YAMBO,       ONLY : yambo_elph_file_name
  !
  IMPLICIT NONE
  !
  REAL(dp)       :: yambo_kpts(3,nksq)
  CHARACTER(300) :: y_file_name
  COMPLEX(DP)    :: gkkp_disk(nbnd,nbnd,nmodes),y_grad_at_gamma(nbnd,nbnd,nat,3),&
&                   y_pol_vec(3*nat,nat,3)
  !
  INTEGER  :: ik, ikk, ikq,  ibnd, jbnd,  mu, i,j
  LOGICAL  :: exst
  !
  CHARACTER(LEN=256) :: elph_dir
  !
  elph_dir='elph_dir/'
  IF (ionode) INQUIRE(file=TRIM(elph_dir), EXIST=exst)
  CALL mp_bcast(exst, ionode_id, intra_image_comm)
  IF (.NOT.exst) CALL create_directory( elph_dir )
  !
  WRITE (6, '(5x,"electron-phonon interaction (to be used in YAMBO)  ..."/)')
  !
  IF ( ionode ) open(unit=99,file=trim(yambo_elph_file_name),form='unformatted')
  IF ( ionode ) write (99) nmodes,nksq,nbnd
  DO ik=1,nksq
     ikk = ikks(ik)
     yambo_kpts(:,ik)=xk(:,ikk)
  END DO
  IF ( ionode ) write (99) alat,xq,yambo_kpts
  IF ( ionode ) write (99) w2(:)
  !
  DO ik=1,nksq
    ikk = ik
    ikq = ik
    IF (.not.lgamma) THEN
      ikk = 2 * ik - 1
      ikq = ikk + 1
    ENDIF
    !
    gkkp_disk=(0.d0, 0.d0)
    !
    DO ibnd=1,nbnd
      DO jbnd=1,nbnd
        DO mu=1,3*nat
          !
          call debye_waller(y_grad_at_gamma(ibnd,jbnd,:,:),y_pol_vec(mu,:,:),&
&                           ibnd,jbnd,ik,mu)
          !
          ! Here we have transition xk(:,ikq) (k+q) -> xk(:,ikk) (k), with ME
          ! el_ph_mat(i,j,k,I)= <\psi(k+q) n_i|dV_{SCF}/du^q_{i a}|\psi(k) n_j>
          !
          ! Note that as wu are small I use
          !
          ! gkk(k,i,j,mu) = \sum_{IJ} el_ph_mat(i,j,k,I) u^*(I,J) dyn(J,mu)
          !
          DO j=1,3*nat
            DO i=1,3*nat
              gkkp_disk(ibnd,jbnd,mu)=gkkp_disk(ibnd,jbnd,mu)+&
&                                     el_ph_mat(ibnd,jbnd,ik,i)*conjg(u(j,i))*dyn(j,mu)
            ENDDO
          ENDDO
        ENDDO !mu
      ENDDO !jbnd
    ENDDO !ibnd
    !
    IF ( ionode ) write(99) gkkp_disk(:,:,:)
    IF ( ionode ) write(99) y_pol_vec(:,:,:)
    IF ( ionode .and. lgamma ) write(99) y_grad_at_gamma(:,:,:,:)
    IF ( ionode ) write (99) et(:nbnd,ikk)
    IF ( ionode ) write (99) et(:nbnd,ikq)
    !
  ENDDO !ik
  IF ( ionode ) close(99)
  !
  RETURN
END SUBROUTINE elph_yambo_eval_and_IO
