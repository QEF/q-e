!
! Copyright (C) 2004 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module is USEd, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

  IMPLICIT NONE
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

!
!---------------------------------------------------------------------
subroutine set_pseudo_upf (is, upf)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !
  ! PWSCF modules
  !
  USE parameters, ONLY: ndmx
  USE atom,  ONLY: zmesh, mesh, dx, r, rab, &
       chi, oc, nchi, lchi, jchi, rho_at, rho_atc, nlcc
  USE pseud, ONLY: lloc, lmax, zp
  USE uspp_param, ONLY: vloc_at, dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, &
       rinner, nbeta, kkbeta, lll, jjj, psd, tvanp
  USE funct, ONLY: set_dft_from_name, dft_is_meta, dft_is_hybrid
  !
  USE ions_base, ONLY: zv, ntyp => nsp
  USE spin_orb, ONLY: lspinorb
  USE pseudo_types
  !
!  USE paw, ONLY : paw_nbeta, aephi, psphi, gipaw_ae_vloc, gipaw_ps_vloc, &
!                  vloc_present, gipaw_data_in_upf_file, &
!                  gipaw_ncore_orbital, gipaw_core_orbital
  USE paw, ONLY : paw_recon
  !
  implicit none
  !
  integer :: is
  !
  !     Local variables
  !
  integer :: nb, mb, ijv, ir
  TYPE (pseudo_upf) :: upf
  !
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  call set_dft_from_name( upf%dft )
  !
#if defined (EXX)
#else
  IF ( dft_is_hybrid() ) &
    CALL errore( 'upf_to_internals ', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  mesh(is) = upf%mesh
  IF ( mesh(is) > ndmx ) &
     CALL errore('upf_to_internals', 'too many grid points', 1)
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)
  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  lmax(is) = upf%lmax
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  do nb = 1, upf%nbeta
     do mb = nb, upf%nbeta
        ijv = mb * (mb-1) / 2 + nb
        qfunc (1:upf%mesh, ijv, is) = upf%qfunc(1:upf%mesh, nb, mb)
     end do
  end do
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )
  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)

  if (lspinorb.and..not.upf%has_so) &
     call infomsg ('upf_to_internal','At least one non s.o. pseudo', -1)
   
  if (upf%has_so) then
     jchi(1:upf%nwfc, is) = upf%jchi(1:upf%nwfc)
     jjj(1:upf%nbeta, is) = upf%jjj(1:upf%nbeta)
  else
     jchi(1:upf%nwfc, is) = 0.d0
     jjj(1:upf%nbeta, is) = 0.d0
  endif
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  !!! TEMP
  lloc(is) = 0
  !!!
  vloc_at(1:upf%mesh,is) = upf%vloc(1:upf%mesh)

  zv(is) = zp(is)  !!! maybe not needed: it is done in setup
  
  !<apsi>
  IF ( upf%has_gipaw ) THEN
     IF ( .NOT. ALLOCATED ( paw_recon ) ) THEN
        ALLOCATE ( paw_recon(ntyp) )
     END IF
     
     paw_recon(is)%paw_nbeta = upf%gipaw_wfs_nchannels
     paw_recon(is)%vloc_present = .TRUE.
     paw_recon(is)%gipaw_data_in_upf_file = .TRUE.
     
     paw_recon(is)%gipaw_ncore_orbital = upf%gipaw_ncore_orbitals
     ALLOCATE ( paw_recon(is)%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
     paw_recon(is)%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals) &
          = upf%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals)
     
     ALLOCATE ( paw_recon(is)%gipaw_ae_vloc(upf%mesh) )
     ALLOCATE ( paw_recon(is)%gipaw_ps_vloc(upf%mesh) )
     paw_recon(is)%gipaw_ae_vloc(:upf%mesh) = upf%gipaw_vlocal_ae(:upf%mesh)
     paw_recon(is)%gipaw_ps_vloc(:upf%mesh) = upf%gipaw_vlocal_ps(:upf%mesh)
     
     ALLOCATE ( paw_recon(is)%aephi(upf%gipaw_wfs_nchannels) )
     ALLOCATE ( paw_recon(is)%psphi(upf%gipaw_wfs_nchannels) )
     
     DO nb = 1, upf%gipaw_wfs_nchannels
        ALLOCATE ( paw_recon(is)%aephi(nb)%psi(mesh(is)) )
        paw_recon(is)%aephi(nb)%label%nt = is
        paw_recon(is)%aephi(nb)%label%n = nb
        paw_recon(is)%aephi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%aephi(nb)%label%m = 
        paw_recon(is)%aephi(nb)%label%nrc = upf%mesh
        paw_recon(is)%aephi(nb)%kkpsi = upf%mesh
        paw_recon(is)%aephi(nb)%label%rc = -1.0_dp
        paw_recon(is)%aephi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ae(:upf%mesh,nb)
        
        ALLOCATE ( paw_recon(is)%psphi(nb)%psi(mesh(is)) )
        paw_recon(is)%psphi(nb)%label%nt = is
        paw_recon(is)%psphi(nb)%label%n = nb
        paw_recon(is)%psphi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%psphi(nb)%label%m = 
        paw_recon(is)%psphi(nb)%label%nrc = upf%mesh
        paw_recon(is)%psphi(nb)%kkpsi = upf%mesh
        paw_recon(is)%psphi(nb)%label%rc = upf%gipaw_wfs_rcutus(nb)
        paw_recon(is)%psphi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ps(:upf%mesh,nb)
     END DO
  ELSE
     paw_recon(is)%gipaw_data_in_upf_file = .FALSE.
  END IF
  !</apsi>
  
end subroutine set_pseudo_upf

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
