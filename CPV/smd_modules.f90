!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Y. Kanai
!
! mod. smd_variables
! mod. smd_rep
! mod. smd_ene
!
!--------------------------------------------------------------------------
!
!
MODULE smd_variables
 !
  implicit none
 !
 !
 ! smx = max integer of sm_p
 ! smwout = index for replica files
 ! mi  = a parameter for 
 !       polynomial interpolation 
 !	 # of replicas used for 
 !	 interpolation
 ! 
 integer, parameter :: smx = 20
 integer, parameter :: smwout = 20
 integer, parameter :: mi = 4
 !
 !
 ! sm_p = 0 .. SM_P replica 
 ! smcp = regular CP calculation
 ! smlm = String method w/ Lagrange Mult.
 ! smopt= CP for 2 replicas, initial & final
 ! linr = linear interpolation
 ! polm = polynomial interpolation  
 ! kwnp = # of points used in polm
 ! tol  = tolrance on const in terms of 
 !	  [alpha(k) - alpha(k-1)] - 1/sm_P   
 ! smfreq = frequency of calculating Lag. Mul 
 ! max_ite = # of such iteration allowed
 !
 integer :: sm_p
 logical :: smcp,smlm,smopt
 logical :: linr,polm,stcd
 integer :: kwnp
 integer :: codfreq, forfreq, smwfreq
 real(kind=8) :: tol
 integer :: lmfreq
 integer :: maxlm
 !
 !
 !--- TYPE pointer ------------
 !
 TYPE ptr
  real(kind=8), pointer :: d3(:,:,:)
 END TYPE ptr
 !
 !
 !
END MODULE smd_variables
!
!
MODULE smd_rep
 !
  use parameters, only: nsx,natx,nacx
 !
  implicit none
 !
 !
 !--- For ionic variables ---------
 !
 TYPE ION
  real(kind=8) :: tau0(3,natx,nsx)
  real(kind=8) :: taup(3,natx,nsx)
  real(kind=8) :: taum(3,natx,nsx)
  real(kind=8) :: taus(3,natx,nsx)
  real(kind=8) :: tausp(3,natx,nsx)
  real(kind=8) :: tausm(3,natx,nsx)
  real(kind=8) :: fion(3,natx,nsx)
  real(kind=8) :: fionm(3,natx,nsx)
  real(kind=8) :: vels(3,natx,nsx)
  real(kind=8) :: velsm(3,natx,nsx)
  real(kind=8) :: tan(3,natx,nsx)
  real(kind=8) :: cdm0(3)
  real(kind=8) :: cdmvel(3)
  real(kind=8) :: acc(nacx)
 END TYPE ION
 !
 !
 !--- For electronic varibles --------
 !
 TYPE ELE
  complex(kind=8), pointer :: c0(:,:)
  complex(kind=8), pointer :: cm(:,:)
  complex(kind=8), pointer :: phi(:,:)
  real(kind=8), pointer :: lambda(:,:)
  real(kind=8), pointer :: lambdam(:,:)
  real(kind=8), pointer :: lambdap(:,:)
  real(kind=8), pointer :: bec(:,:)
  real(kind=8), pointer :: rhovan(:,:,:)
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
  end subroutine
 !
 !
END MODULE smd_rep

MODULE smd_ene
 !
  implicit none
 !
  real(kind=8), allocatable :: etot_ar(:)
  real(kind=8), allocatable :: ekin_ar(:)
  real(kind=8), allocatable :: eht_ar(:)
  real(kind=8), allocatable :: epseu_ar(:)
  real(kind=8), allocatable :: exc_ar(:)
  real(kind=8), allocatable :: esr_ar(:)
 !
 CONTAINS
  subroutine deallocate_smd_ene()
      IF( ALLOCATED( etot_ar  ) ) DEALLOCATE( etot_ar  )
      IF( ALLOCATED( ekin_ar  ) ) DEALLOCATE( ekin_ar  )
      IF( ALLOCATED( eht_ar   ) ) DEALLOCATE( eht_ar   )
      IF( ALLOCATED( epseu_ar ) ) DEALLOCATE( epseu_ar )
      IF( ALLOCATED( exc_ar   ) ) DEALLOCATE( exc_ar   )
      IF( ALLOCATED( esr_ar   ) ) DEALLOCATE( esr_ar   )
  end subroutine
 !
END MODULE smd_ene


