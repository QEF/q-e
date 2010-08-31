! FOR GWW
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Created by Joe Stenuit - 05/06/2009
!
!----------------------------------------------------------------------------
SUBROUTINE energies_u_gamma( e_u )
  !----------------------------------------------------------------------------
  !
  ! #ifdef __GWW
  ! This routine computes the expectation value of the Hubbard potential
  ! The result is e_u
  ! ... output:
  !       e_u   expection values due to the Hubbard U
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : rytoev
  USE uspp,                 ONLY : okvan, vkb, nkb
  USE wvfct,                ONLY : igk, npwx, npw, nbnd
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE becmod,               ONLY : bec_type, becp, allocate_bec_type, calbec
  !
#ifdef EXX
  USE exx,      ONLY : vexx !Suriano
  USE funct,    ONLY : exx_is_active
#endif

  !
  IMPLICIT NONE
  !
  ! ... input/output arguments
  !
  REAL(kind=DP), intent(out) :: e_u(nbnd)
  !
  INTEGER  :: ibnd, ig
  COMPLEX(kind=DP), ALLOCATABLE :: hupsi(:,:)
  !
  IF ( .NOT. gamma_only ) THEN
     !
     write(stdout,*) 'NOT YET IMPLEMENTED !'
     stop
     !
  END IF
  !
  !
  !if required calculates also the KS energies
  write(stdout,*) 'Start subroutine energies_u_gamma'
  write(stdout,*) 'nkb=', nkb, 'nbnd=', nbnd
  call flush_unit(stdout)
  !
  call allocate_bec_type (nkb, nbnd,becp)
  !
  IF ( nkb > 0 )  CALL init_us_2( npw, igk, xk(1,1), vkb )
  call calbec(npw, vkb, evc, becp, nbnd) !! N.B. nbnd is optional !
  !
  allocate(hupsi(npwx,nbnd))
  !
  ! since in vhpsi, the result is added to hupsi =>
  ! we have to initialize all the hupsi to zero
  hupsi(:,:)=0.d0
  !
  write(stdout,*) "ubound hupsi=", ubound(hupsi,1), ' X ', ubound(hupsi,2)
  write(stdout,*) "ubound evc=", ubound(evc,1), ' X ', ubound(evc,2)
  call flush_unit(stdout)
  !
  CALL vhpsi( npwx, npw, nbnd, evc, hupsi )
  !
  write(stdout,*) "CALL vhpsi done"
  write(stdout,*) "ubound hupsi=", ubound(hupsi,1), ' X ', ubound(hupsi,2)
  write(stdout,*) "ubound evc=", ubound(evc,1), ' X ', ubound(evc,2)
  call flush_unit(stdout)
  !
  do ibnd=1,nbnd
      write(stdout,*) 'hupsi(1', ibnd ,')= ', hupsi(1,ibnd)
  enddo
  call flush_unit(stdout)
  !
  e_u(:)=0.d0
  !
  !!! for info:
  write(stdout,*) 'npwx=', npwx, ' and evc(1,1)=', evc(1,1)
  write(stdout,*) 'gstart=', gstart, ' and nbnd=', nbnd, ' and npw=', npw
  call flush_unit(stdout)
  do ibnd=1,nbnd
      !
      do ig=1,npw
          e_u(ibnd)=e_u(ibnd)+2.d0*dble(conjg(evc(ig,ibnd))*hupsi(ig,ibnd))
      enddo
      if(gstart==2) then
          e_u(ibnd)=e_u(ibnd)-dble(conjg(evc(1,ibnd))*hupsi(1,ibnd))
      endif
  enddo
  !
  call mp_sum(e_u(:))
  !
  do ibnd=1,nbnd
      write(stdout,*) '<ener_u>(', ibnd,')= ', e_u(ibnd)*rytoev
  enddo
  call flush_unit(stdout)
  deallocate(hupsi,becp%r)
  !
  return
  !
END SUBROUTINE energies_u_gamma
  !
  !
SUBROUTINE write_energies_u(e_u_in)
  !
  ! #ifdef __GWW
  !
  USE kinds,      ONLY : DP
  USE wannier_gw, ONLY : num_nbnds, nbnd_normal
  USE io_files,   ONLY : find_free_unit, prefix
  USE io_global,  ONLY : ionode
  USE wvfct,      ONLY : nbnd
  !
  implicit none
  !
  REAL(kind=DP), INTENT(in) :: e_u_in(nbnd)!exchange and correlation energies
  !
  INTEGER :: iunu, iw
  !
  !
  if(ionode) then
     iunu = find_free_unit()
     !
     open(unit=iunu,file=trim(prefix)//'.hubbard_u',status='unknown',form='unformatted')
     !open(unit=iunu,file=trim(prefix)//'.hubbard_u',status='unknown',form='unformatted')
     !
     write(iunu) nbnd_normal
     !
     do iw=1,nbnd_normal
        write(iunu) e_u_in(iw)
     enddo
     !
     close(iunu)
  endif
  !
  ! #endif __GWW
  !
  return
  !
END SUBROUTINE write_energies_u
