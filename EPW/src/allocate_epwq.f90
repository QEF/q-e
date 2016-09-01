  !                                                                        
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino     
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/allocate_phq - Quantum-ESPRESSO group
  ! 09/2009 There is a lot of excess in this file.  
  !
  !----------------------------------------------------------------------- 
  subroutine allocate_epwq
  !-----------------------------------------------------------------------
  !!
  !! Dynamical allocation of arrays: quantities needed for the linear
  !! response problem
  !!
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  USE ions_base,    ONLY : nat, ntyp => nsp
  USE pwcom,        ONLY : npwx, nbnd, ngm, nspin, nks
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions_module,  ONLY: evc
  USE spin_orb,     ONLY : lspinorb
  USE control_lr,   ONLY : lgamma, nbnd_occ
  USE phcom,        ONLY : evq, dvpsi, dpsi, vlocq,&
                           dmuxc, npertx 
  USE phus,         ONLY : int1, int1_nc, int2, int2_so, &
                           int4, int4_nc, int5, int5_so, becsum_nc, &
                           alphasum, alphasum_nc, alphap
  USE lr_symm_base, ONLY : rtau               
  USE qpoint,       ONLY : eigqts
  USE lrus,         ONLY : becp1, int3, int3_nc
  USE elph2,        ONLY : elph, el_ph_mat
  USE becmod,       ONLY : becp, allocate_bec_type
  USE uspp_param,   ONLY : nhm
  USE uspp,         ONLY : okvan, nkb
! SP creation of ffts
  USE units_ph,     ONLY : this_dvkb3_is_on_file, this_pcxpsi_is_on_file
  USE modes,        ONLY : u, npert, name_rap_mode, num_rap_mode
  USE fft_base,     ONLY : dtgs
  USE klist,        ONLY : nks

  implicit none
  ! SP: Had to add these allocations becaue they are now private in QE 5.0.
  !     See private 'becp_nc' variable from espresso/Modules/becmod.f90. 
  INTEGER :: ik, ipol
!  COMPLEX(DP), ALLOCATABLE ::  &
!       becp_nc(:,:,:)   !  <beta|psi> for real (at Gamma) wavefunctions for spinors

  !
  !  allocate space for the quantities needed in EPW
  !
  IF (lgamma) THEN
     !
     !  q=0  : evq is a pointers to evc
     !
     evq  => evc
  ELSE
     !
     !  q!=0 : evq is ALLOCATEd and calculated at point k+q
     !
     ALLOCATE (evq ( npwx*npol, nbnd))    
  ENDIF
  !
  ALLOCATE (dvpsi ( npwx*npol, nbnd))    
  ALLOCATE ( dpsi ( npwx*npol, nbnd))    
  !
  ALLOCATE (vlocq ( ngm, ntyp))    
! SP: nrxx is not used in QE 5 ==> tg_nnr is the maximum among nnr
!     This SHOULD have the same dim as nrxx had.
!  ALLOCATE (dmuxc ( nrxx, nspin, nspin))  
! SP: Again a new change in QE (03/08/2016)  
!  ALLOCATE (dmuxc ( dffts%tg_nnr, nspin, nspin))    
  ALLOCATE (dmuxc ( dtgs%tg_nnr, nspin, nspin))    
  !
  ALLOCATE (eigqts ( nat))    
  ALLOCATE (rtau ( 3, 48, nat))    
  ALLOCATE (u ( 3 * nat, 3 * nat))    
!  ALLOCATE (t (npertx, npertx, 48,3 * nat))    
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat ))
  ALLOCATE (npert ( 3 * nat))    
  IF (okvan) THEN
     ALLOCATE (int1 ( nhm, nhm, 3, nat, nspin))    
     ALLOCATE (int2 ( nhm, nhm, 3, nat, nat))    
     ALLOCATE (int3 ( nhm, nhm, npertx, nat, nspin))    
     ALLOCATE (int4 ( nhm * (nhm + 1)/2, 3, 3, nat, nspin))    
     ALLOCATE (int5 ( nhm * (nhm + 1)/2, 3, 3, nat , nat))    
     IF (noncolin) THEN
        ALLOCATE(int1_nc( nhm, nhm, 3, nat, nspin))
        ALLOCATE(int3_nc( nhm, nhm, npertx, nat, nspin))
        ALLOCATE(int4_nc( nhm, nhm, 3, 3, nat, nspin))
        ALLOCATE(becsum_nc( nhm*(nhm+1)/2, nat, npol, npol))
        ALLOCATE(alphasum_nc( nhm*(nhm+1)/2, 3, nat, npol, npol))
        IF (lspinorb) THEN
           ALLOCATE(int2_so( nhm, nhm, 3, nat, nat, nspin))
           ALLOCATE(int5_so( nhm, nhm, 3, 3, nat, nat, nspin))
        END IF
     END IF
     ALLOCATE (alphasum ( nhm * (nhm + 1)/2, 3, nat , nspin))    
     ALLOCATE (this_dvkb3_is_on_file(nks))    
     this_dvkb3_is_on_file(:)=.false.
  ENDIF
  ALLOCATE (this_pcxpsi_is_on_file(nks,3))
  this_pcxpsi_is_on_file(:,:)=.false.
! SP : from new PHonon/PH/allocate_phq.f90 
  ALLOCATE (becp1(nks))
  ALLOCATE (alphap(3,nks))
  ALLOCATE(nbnd_occ(nks))
  DO ik=1,nks
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbnd, alphap(ipol,ik) )
     ENDDO
  END DO
  CALL allocate_bec_type ( nkb, nbnd, becp )

  IF (elph) ALLOCATE (el_ph_mat( nbnd, nbnd, nks, 3*nat))    
END SUBROUTINE allocate_epwq
