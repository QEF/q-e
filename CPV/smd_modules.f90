!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Y. Kanai
!
! mod. smd_variables  ==> see Modules/path_variables.f90
! mod. smd_rep
! mod. smd_ene
!
!--------------------------------------------------------------------------
!
!
!
MODULE smd_rep
 !
  use parameters, only: natx,nacx
 !
  implicit none
 !
 !
 !--- For ionic variables ---------
 !
 TYPE ION
  real(8) :: tau0(3,natx)
  real(8) :: taup(3,natx)
  real(8) :: taum(3,natx)
  real(8) :: taus(3,natx)
  real(8) :: tausp(3,natx)
  real(8) :: tausm(3,natx)
  real(8) :: fion(3,natx)
  real(8) :: fionm(3,natx)
  real(8) :: vels(3,natx)
  real(8) :: velsm(3,natx)
  real(8) :: tan(3,natx)
  real(8) :: cdm0(3)
  real(8) :: cdmvel(3)
  real(8) :: acc(nacx)
 END TYPE ION
 !
 !
 !--- For electronic varibles --------
 !
 TYPE ELE
  complex(8), pointer :: c0(:,:)
  complex(8), pointer :: cm(:,:)
  complex(8), pointer :: phi(:,:)
  real(8), pointer :: lambda(:,:)
  real(8), pointer :: lambdam(:,:)
  real(8), pointer :: lambdap(:,:)
  real(8), pointer :: bec(:,:)
  real(8), pointer :: rhovan(:,:,:)
 END TYPE ELE
 !
 !
 !--- VARIABLES --------------- 
 !
 TYPE(ION), allocatable, target :: rep(:)
 TYPE(ELE), allocatable :: rep_el(:)
 !
 CONTAINS
    subroutine deallocate_smd_rep()
      IF( ALLOCATED( rep    ) ) DEALLOCATE( rep    )
      IF( ALLOCATED( rep_el ) ) DEALLOCATE( rep_el )
  end subroutine deallocate_smd_rep
 !
 !
END MODULE smd_rep

MODULE smd_ene
 !
  implicit none
 !
  real(8), allocatable :: etot_ar(:)
  real(8), allocatable :: ekin_ar(:)
  real(8), allocatable :: eht_ar(:)
  real(8), allocatable :: epseu_ar(:)
  real(8), allocatable :: exc_ar(:)
  real(8), allocatable :: esr_ar(:)
 !
 CONTAINS
  subroutine deallocate_smd_ene()
      IF( ALLOCATED( etot_ar  ) ) DEALLOCATE( etot_ar  )
      IF( ALLOCATED( ekin_ar  ) ) DEALLOCATE( ekin_ar  )
      IF( ALLOCATED( eht_ar   ) ) DEALLOCATE( eht_ar   )
      IF( ALLOCATED( epseu_ar ) ) DEALLOCATE( epseu_ar )
      IF( ALLOCATED( exc_ar   ) ) DEALLOCATE( exc_ar   )
      IF( ALLOCATED( esr_ar   ) ) DEALLOCATE( esr_ar   )
  end subroutine deallocate_smd_ene
 !
END MODULE smd_ene


