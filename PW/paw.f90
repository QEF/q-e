!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE paw

  USE kinds, ONLY: DP
  USE parameters, ONLY: nbrx, npsx, ndmx
  !
  ! ... These parameters are needed for the paw variables
  !
  SAVE
  !
  REAL(DP) :: &
       paw_betar(ndmx,nbrx,npsx)  ! radial beta_{mu} functions
  INTEGER :: &
       paw_nh(npsx),             &! number of beta functions per atomic type
       paw_nbeta(npsx),          &! number of beta functions
       paw_kkbeta(npsx),         &! point where the beta are zero
       paw_lll(nbrx,npsx)         ! angular momentum of the beta function
  INTEGER :: &
       paw_nhm,              &! max number of different beta functions per atom
       paw_nkb,              &! total number of beta functions, with st.fact.
       paw_nqxq,             &! size of interpolation table
       paw_lmaxkb,           &! max angular momentum
       paw_lmaxq,              &! max angular momentum + 1 for Q functions
       paw_nqx                ! number of interpolation points
  INTEGER, ALLOCATABLE ::&
       paw_indv(:,:),        &! correspondence of betas atomic <-> soli
       paw_nhtol(:,:),       &! correspondence n <-> angular momentum
       paw_nhtom(:,:),       &! correspondence n <-> magnetic angular m
       paw_nl(:,:),          &! number of projectors for each l
       paw_iltonh(:,:,:)        ! corresp l, num <--> n for each type
  complex(DP), ALLOCATABLE, TARGET :: &
       paw_vkb(:,:),         &   ! all beta functions in reciprocal space
       paw_becp(:,:)             !  products of wavefunctions and proj
  REAL(DP), ALLOCATABLE :: &
       paw_tab(:,:,:)              ! interpolation table for PPs
#ifdef USE_SPLINES
  !<ceres>
  REAL(DP), ALLOCATABLE :: &
       paw_tab_d2y(:,:,:)          ! for cubic splines
  !</ceres>
#endif
  !
  type wfc_label
     integer  :: na , &   ! Atom number
          nt ,        &   ! Type
          n  ,        &   ! Chi index
          l  ,        &   ! l
          m  ,        &   ! m
          nrc             ! indice of core radius in mesh
     real(DP) :: rc  ! paw core radius
  end type wfc_label

  type at_wfc
     type(wfc_label)          :: label
     integer                  :: kkpsi
!     real(DP)            :: rmt   = 0.0_DP ! Like FLAPW or LMTO Muffin Tinradius
     real(DP)  , pointer :: psi(:)
  end type at_wfc

  type(at_wfc),pointer :: aephi(:,:), psphi(:,:) ! Atom

CONTAINS

  subroutine paw_wfc_init(phi)
    ! 
    ! Initialize default values for labe end kkpsi
    !

    type(at_wfc) :: phi(:,:)

    phi%label%na = 0
    phi%label%nt = 0
    phi%label%n  = 0
    phi%label%l  = -99
    phi%label%m  = -99
    phi%label%nrc = 0
    phi%kkpsi    = 0

    return
  end subroutine paw_wfc_init

  subroutine read_recon(filerec)

  !
  ! Read all-electron and pseudo atomic wavefunctions 
  !  needed for PAW reconstruction
  !

  use read_upf_module, only: scan_begin, scan_end
  USE ions_base,          ONLY : ntyp => nsp
  use atom, only: mesh 
  use kinds, only: DP
  use parameters, only : ntypx
  USE io_global,  ONLY : stdout
  implicit none

  character (len=256) :: filerec(ntypx)
  integer :: j,i,jtyp,kkphi,nbetam

  do jtyp=1,ntyp
     open(14,file=filerec(jtyp))
     call scan_begin(14,'PAW',.true.)
     read(14,*) paw_nbeta(jtyp)
     call scan_end(14,'PAW')
     close(14)
  enddo
  nbetam=maxval(paw_nbeta)
  allocate( psphi(ntyp,nbetam) )
  allocate( aephi(ntyp,nbetam) )

  call paw_wfc_init(psphi)
  call paw_wfc_init(aephi)


  recphi_read: do jtyp=1,ntyp
     open(14,file=filerec(jtyp))
     write (stdout,*) "N_AEwfc atom",jtyp,":",paw_nbeta(jtyp)
     recphi_loop: do i=1,paw_nbeta(jtyp)
        allocate(aephi(jtyp,i)%psi(maxval(mesh(1:ntyp))))
        aephi(jtyp,i)%label%nt=jtyp
        aephi(jtyp,i)%label%n=i
        call scan_begin(14,'REC',.false.)
        call scan_begin(14,'kkbeta',.false.)
        read(14,*)  kkphi
        call scan_end(14,'kkbeta')
        aephi(jtyp,i)%kkpsi=kkphi
        call scan_begin(14,'L',.false.)        
        read(14,*)  aephi(jtyp,i)%label%l
        call scan_end(14,'L')
        call scan_begin(14,'REC_AE',.false.)
        read(14,*) (aephi(jtyp,i)%psi(j),j=1,kkphi)
        call scan_end(14,'REC_AE')
        allocate (psphi(jtyp,i)%psi(maxval(mesh(1:ntyp))))
        psphi(jtyp,i)%label%nt=jtyp
        psphi(jtyp,i)%label%n=i
        psphi(jtyp,i)%label%l=aephi(jtyp,i)%label%l
        psphi(jtyp,i)%kkpsi=kkphi
        call scan_begin(14,'REC_PS',.false.)
        read(14,*) (psphi(jtyp,i)%psi(j),j=1,kkphi)
        call scan_end(14,'REC_PS')
        call scan_end(14,'REC')
     end do recphi_loop
     close(14)
  end do recphi_read

end subroutine read_recon


END MODULE paw
