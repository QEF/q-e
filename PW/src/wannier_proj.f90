! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

subroutine wannier_proj(ik, wan_func)
  ! This routine computes <phi_i|S|psi_j> for all eigenvectors
  ! for current k-point

  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout 
  USE io_files
  USE wannier_new,      ONLY : wan_in, nwan, use_energy_int
  USE ions_base,        ONLY : nat, ityp
  USE wvfct,            ONLY : nbnd, npwx, et
  USE klist,            ONLY : ngk
  USE lsda_mod,         ONLY : lsda, isk
  USE constants,        ONLY : rytoev
  USE basis,            ONLY : swfcatom
  USE control_flags,    ONLY : gamma_only
  USE uspp_param,       ONLY : upf 
  USE wavefunctions_module, ONLY : evc
  USE gvect,                ONLY : gstart
  USE noncollin_module, ONLY : npol
  USE buffers,          ONLY : get_buffer, save_buffer

  
  implicit none
  ! input-output
  INTEGER, intent(in) :: ik
  COMPLEX(DP), intent(out) :: wan_func(npwx,nwan)
  !
  COMPLEX(DP), ALLOCATABLE :: pp(:,:)
  COMPLEX(DP), ALLOCATABLE :: trialwf(:,:)
  
  INTEGER :: current_spin, npw, i,j,k, ierr, ibnd, iwan
  
  REAL(DP), EXTERNAL :: ddot
  COMPLEX(DP) :: zdotc

  ALLOCATE(trialwf(npwx,nwan))
  ALLOCATE(pp(nwan, nbnd))
  
  current_spin = 1
  IF (lsda) current_spin  = isk(ik)
  npw = ngk(ik)
  
  ! Read current wavefunctions
  !
  evc = ZERO
  ! See comment in PP/src/openfil.f90 why davcio and not get_buffer
  ! call get_buffer ( evc, nwordwfc, iunwfc, ik )  
  call davcio ( evc, 2*nwordwfc, iunwfc, ik, -1 )  
  ! Reads ortho-atomic wfc
  ! You should prepare data using orthoatwfc.f90
  swfcatom = ZERO
  CALL get_buffer (swfcatom, nwordatwfc, iunsat, ik)
  
  ! generates trial wavefunctions as a sum of ingredients
  trialwf = ZERO
  do iwan=1, nwan
     do j=1,wan_in(iwan,current_spin)%ning
        do k=1,npw
           trialwf(k,iwan) = trialwf(k,iwan) + &
                CMPLX(wan_in(iwan,current_spin)%ing(j)%c,0.d0,KIND=DP) * &
                swfcatom(k,wan_in(iwan,current_spin)%ing(j)%iatomwfc)
        end do
     end do
  end do
  
  ! computes <\Psi|\hat S|\phi> for all \Psi and \phi
  ! later one should select only few columns 
  pp = ZERO
  DO ibnd = 1, nbnd
     DO iwan = 1, nwan
        pp (iwan, ibnd) = zdotc (npw, trialwf (1, iwan), 1, evc (1, ibnd), 1)
     ENDDO
  ENDDO

  ! And now we should nullify few elements
  do iwan=1, nwan
     do ibnd=1, nbnd
        if(use_energy_int) then
           if( et(ibnd,ik) < wan_in(iwan,current_spin)%bands_from ) pp(iwan,ibnd) = ZERO
           if( et(ibnd,ik) > wan_in(iwan,current_spin)%bands_to ) pp(iwan,ibnd) = ZERO
        else
           if( (ibnd < INT(wan_in(iwan,current_spin)%bands_from)) &
                .OR. ( ibnd > INT(wan_in(iwan,current_spin)%bands_to) )) then
              pp(iwan,ibnd) = ZERO
              ! write(stdout,'(5x,"nullify component for band",i3," of wannier",i3)') ibnd,iwan
           end if
        end if
     end do
  end do

  ! Orthogonalize pp
  CALL ortho_wfc(nwan,nbnd,pp,ierr)
  IF (ierr .NE. 0) &
    write(stdout,'(/5x,a43,i4,a12,i2/)') 'Wannier orthogonalization failed on k-point', ik, 'with status', ierr

  !And write ortho-pp to file
  call save_buffer( pp, nwordwpp, iunwpp, ik)

  wan_func = ZERO
  
  call ZGEMM('N', 'C', npw, nwan, nbnd, ONE, evc, &
       npwx, pp, nwan, ZERO, wan_func, npwx)

  !And dump wannier to file
  call save_buffer( wan_func, nwordwf, iunwf, ik)

  DEALLOCATE(trialwf)
  DEALLOCATE(pp)
  
  RETURN
  !
END SUBROUTINE wannier_proj
