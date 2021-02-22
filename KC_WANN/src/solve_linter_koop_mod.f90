! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"

MODULE solve_linter_koop_mod

CONTAINS
!
!-----------------------------------------------------------------------
subroutine solve_linter_koop ( spin_ref, i_ref, delta_vr, drhog_scf, delta_vg, drhor_scf)
  !-----------------------------------------------------------------------
  !
  !!    Driver routine for the solution of the linear system which defines the 
  !!    change of the wavefunction due to an external perturbation. It work 
  !!    exactly as solve_linter in PH but for a genenral perturbation dv delta_vr.
  !!    In genereal, it performs the following tasks:
  !!     a) computes the bare potential term Delta V | psi >, add the scf contribution
  !!     b) applies P_c^+ (orthogonalization to valence states)
  !!     c) calls cgstabsolve_all to solve the linear system
  !!     d) computes Delta rho
  !
  ! ## FIXME delta_vg passed only for debug reason (to test symmetries) Needs to be removed once the
  !          problem will be solved
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout
  USE check_stop,            ONLY : check_stop_now
  USE wavefunctions,         ONLY : evc, psic
  USE klist,                 ONLY : lgauss, wk, xk, igk_k, ngk
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : fwfft, invfft, fft_interpolate
  USE gvect,                 ONLY : gstart
  USE gvecs,                 ONLY : doublegrid, ngms
  USE becmod,                ONLY : calbec
  USE wvfct,                 ONLY : npw, npwx, nbnd, et
  USE uspp_param,            ONLY : nhm
  USE control_ph,            ONLY : niter_ph, nmix_ph, tr2_ph, &
                                    convt, &
                                    alpha_mix, flmixdpot
  USE control_lr,            ONLY : lgamma, nbnd_occ
  USE units_ph,              ONLY : iudwf, lrdwf, iubar, lrbar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE buffers,               ONLY : save_buffer, get_buffer
  USE eqv,                   ONLY : dvpsi, dpsi, evq
  USE qpoint,                ONLY : npwq, nksq, ikks, ikqs
  USE uspp,                  ONLY : vkb, okvan  
  USE mp,                    ONLY : mp_sum
  USE mp_global,             ONLY : intra_pool_comm, inter_pool_comm 
  USE noncollin_module,      ONLY : npol
  USE control_kc_wann
  USE dv_of_drho_lr
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm
  !
  !USE cell_base,            ONLY : omega
  !
  implicit none
  !
  integer, parameter :: npe = 1
  !
  complex(dp), intent(in) ::delta_vr (dffts%nnr, nspin)
  complex(dp), intent(out) ::  drhog_scf (ngms, nspin)
  complex(dp), intent(out), optional :: drhor_scf(dffts%nnr,nspin)
  !
  ! input: the imaginary frequency
  real(DP) :: aux_avg (2), averlt
  !
  !complex(DP) :: drhoscf (nrxx, nspin, 1)
  complex(DP), allocatable :: drhoscf (:,:,:)
  real(DP) , allocatable   :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  real(DP) :: thresh, anorm, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight
  ! dos_ef: DOS at Fermi energy (in calculation for a metal)
  ! weight: weight of k-point     
  complex(DP), allocatable, target :: dvscfin(:,:,:)
  ! change of the scf potential 
  complex(DP), pointer :: dvscfins (:,:,:)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:,:), dvscfout (:,:,:)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:),&
       dbecsum (:,:,:,:), aux(:), aux1 (:,:), aux2(:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  complex(DP), allocatable :: drhoc(:)
  ! drhoc: response core charge density
  REAL(DP), allocatable :: becsum1(:,:,:)
  ! for 
  !
  logical :: conv_root,  & ! true if linear system is converged
             lmetq0        ! true if xq=(0,0,0) in a metal

  integer :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ltaver,     & ! average counter
             ig,         & ! counter on G vectors
             ir,         & ! counter on mesh points
             is,         & ! counter on spin polarizations
             ipert,      & ! counter on perturbations
             spin_ref,   & ! the spin of the reference orbital
             i_ref,      & ! the orbital we want to keep fix
             nrec          ! the record number for dvpsi and dpsi
           


  real(DP) :: tcpu, get_clock ! timing variables

  external cch_psi_all, ccg_psi
  external ch_psi_all,  cg_psi

  !!## DEBUG
  complex(DP) :: drhok(dffts%nnr)
  complex(DP) :: delta_vg(ngms,nspin)
  !!## DEBUG 
  !
  call start_clock ('solve_linter')
  !
  allocate( h_diag(npwx, nbnd) )

  allocate (dvscfin (dfftp%nnr, nspin , npe))    
  if (doublegrid) then
     allocate (dvscfins (dffts%nnr , nspin , npe))    
  else
     dvscfins => dvscfin
  endif
  !allocate (drhoscf  (dfftp%nnr, nspin, 1) )
  allocate (drhoscf  (dffts%nnr, nspin, 1) ) !! NsC
  allocate (drhoscfh (dfftp%nnr, nspin , npe))    
  allocate (dvscfout (dfftp%nnr, nspin , npe))    
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin , npe))    
  allocate (aux ( dffts%nnr ))    
  allocate (aux1 ( dffts%nnr, npol ))    
  allocate (aux2(npwx*npol, nbnd))
  allocate (drhoc(dfftp%nnr))
  !
  IF (kc_at_ks .AND. fix_orb) WRITE (stdout, '("FREEZING ORBITAL #", i4, 3x , "spin", i4)') i_ref, spin_ref
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  lmetq0 = lgauss.and.lgamma
  if (lmetq0) then
     allocate ( ldos ( dfftp%nnr, nspin) )    
     allocate ( ldoss( dffts%nnr, nspin) )    
     allocate (becsum1 ( (nhm * (nhm + 1))/2 , nat , nspin))
     call localdos ( ldos , ldoss , becsum1, dos_ef )
  endif
  !
  !   The outside loop is over the iterations
  !
  dr2=0.d0
  iter0 = 0 ! We do not have a restart procedure yet: NsC
  do kter = 1, niter_ph
     iter = kter + iter0
     !
     ltaver = 0
     !
     lintercall = 0
     drhoscf(:,:,:) = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     !
     do ik = 1, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        if (lsda) current_spin = isk (ikk)
        !
        ! read unperturbed wavefunctions psi(k) and psi(k+q)
        !
        if (nksq.gt.1) then
           if (lgamma) then
              call get_buffer (evc, lrwfc, iuwfc, ikk)
           else
              call get_buffer (evc, lrwfc, iuwfc, ikk)
              call get_buffer (evq, lrwfc, iuwfc, ikq)
           endif
        endif
        !
        ! compute beta functions and kinetic energy for k-point ikq
        ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
        !
        CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
        CALL g2_kin (ikq) 
        !
        ! compute preconditioning matrix h_diag used by cgsolve_all
        !
        CALL h_prec (ik, evq, h_diag)
        !
        ipert = 1
        nrec = (ipert - 1) * nksq + ik
        !
        ! compute the right hand side of the linear system due to
        ! the perturbation, dvscfin used as work space
        !
        IF (iter == 1 ) THEN 
           ! 
           ! ... IF first iteration compute dv_bare*psi and set dvscf to zero
           !
           dvpsi(:,:) = (0.D0, 0.D0)
           do ibnd = 1, nbnd_occ (ik)
              aux(:) = (0.d0, 0.d0)
              do ig = 1, npw
                 aux (dffts%nl(igk_k(ig,ikk)))=evc(ig,ibnd)
              enddo
              CALL invfft ('Wave', aux, dffts)
              do ir = 1, dffts%nnr
                  aux(ir)=aux(ir)*delta_vr(ir,current_spin) 
              enddo
              !
              CALL fwfft ('Wave', aux, dffts)
              do ig = 1, npwq
                 dvpsi(ig,ibnd)=aux(dffts%nl(igk_k(ig,ikq)))
              enddo
              !
           enddo
           if (okvan) then
              call errore('solve_linter_koop', 'USPP not implemented yet', 1)
           endif
           !
           call save_buffer (dvpsi, lrbar, iubar, nrec)        
           dvscfin(:,:,1) = (0.D0, 0.D0)
           !
           thresh = 1.d-2
           thresh = 1.d-6
           !
           !
        ELSE 
           ! ... Else read dv_bare*psi from file 
           ! ... The scf potential from previous iteration (stored in dvscfin) is applied to the wfc
           ! ... add it to dv_bare*psi
           !
           !
           call get_buffer (dvpsi, lrbar, iubar, nrec)
           !
           aux1(:,:) = (0.D0,0.D0)
           aux2(:,:) = (0.D0,0.D0)
           !
           DO ibnd = 1, nbnd_occ (ikk)
              !
              ! 
              call cft_wave (ik, evc (1, ibnd), aux1, +1)
              call apply_dpot(dffts%nnr,aux1, dvscfins(1,1,ipert), current_spin)
              call cft_wave (ik, aux2 (1, ibnd), aux1, -1)
              !
           ENDDO
           !
           dvpsi=dvpsi+aux2
           do ibnd = 1, nbnd_occ (ik)
              !
              !WRITE(*,*) 'ibnd=', ibnd
              !WRITE(*,*) dvpsi(1:3,ibnd)
           enddo
           !
           thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
           thresh = min (1.d-2 * sqrt (dr2), 1.d-6)
           !
        ENDIF
        !
        !
        ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
        ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
        !
        ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+. 
        ! And finally |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
        !
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .false.)
        !
        weight = wk (ikk); conv_root = .true.; dpsi(:,:)=(0.d0,0.d0)
        call cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi, &
                          h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                          anorm, nbnd_occ(ikk), npol)
        !
        if (nksq.gt.1) call davcio ( dpsi, lrdwf, iudwf, ik, +1)
        !
        ltaver = ltaver + lter
        lintercall = lintercall + 1
        if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4,  &
             &              " solve_linter_iu: root not converged ",e10.3)') &
             &              ik , ibnd, anorm
          
        !
        ! ...            Calculates drho, sum over k             ...
        ! ... Eventually freeze the orbital i_ref given in input ...
        !
        IF (kc_at_KS .AND. fix_orb) THEN
          IF (spin_ref == current_spin) THEN 
             dpsi(:,i_ref) = (0.D0, 0.D0)  !! disregard the variation of the orbital i_ref 
          ENDIF
        ENDIF
        !
        call incdrhoscf (drhoscf(1,current_spin,1), weight, ik, &
                         dbecsum(1,1,current_spin,1), dpsi)
       !




!!## DEBUG 
!       !! This is to check which k points are effectively equivalent (DEBUG/UNDERSTAND)
!       !! only sum over the band
!       drhok = CMPLX(0.D0,0.D0,kind=DP)
!       call incdrhoscf ( drhok, weight, ik, &
!                         dbecsum(1,1,current_spin,1), dpsi)
!       aux(:) = drhok(:)
!       CALL fwfft ('Rho', aux, dffts)
!       !! drhog_scf used as workspace
!       drhog_scf(:,1) = aux(dffts%nl(:))
!       !!WRITE(*,'("NICOLA", 6f22.18)') drhok(1:3) 
!       !!WRITE(*,'("NICOLA", 6f22.18)') dpsi(1:3,1) 
!       !!WRITE(*,'("NICOLA", 6f22.18)')  drhog_scf(1:3,1)
!       !!WRITE(*,'("NICOLA", 6f22.18)')  delta_vg(1:3,1)
!       WRITE(*,'("NICOLA", i5, f12.4, 3x, 2f12.8)') ik, weight, sum (CONJG(drhog_scf(:,1)) * delta_vg(:,1))*omega
!!!## DEBUG




       !
     enddo ! on k-points
     !
#ifdef __MPI
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     call mp_sum (dbecsum, intra_pool_comm)
     ! 
#endif
     !
     if (doublegrid) then
        do is = 1, nspin
           !call cinterpolate (drhoscfh(1,is,1), drhoscf(1,is,1), 1)
           call fft_interpolate (dffts, drhoscf(:,is,1), dfftp, drhoscfh(:,is,1))
        enddo
     else
        call zcopy (npe*nspin*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     endif
     !
     ! if q=0, make sure that charge conservation is guaranteed
     !
     if ( lgamma ) then
        psic(:) = drhoscfh(:, nspin, npe)
        CALL fwfft ('Rho', psic, dfftp)
        !CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
        if ( gstart==2) psic(dfftp%nl(1)) = (0.d0, 0.d0)
        CALL invfft ('Rho', psic, dfftp)
        !CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
        drhoscfh(:, nspin, npe) = psic(:)
     endif
     !
     !    Now we compute for all perturbations the total charge and potential
     !
     !call addusddens (drhoscfh, dbecsum, irr, imode0, npe, 0)
     !
#ifdef __MPI
     !
     !   Reduce the delta rho across pools
     !
     call mp_sum (drhoscf, inter_pool_comm)
     call mp_sum (drhoscfh, inter_pool_comm)
     !
#endif
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !
     do ipert = 1, npe
        !
        call zcopy (dfftp%nnr*nspin,drhoscfh(1,1,ipert),1,dvscfout(1,1,ipert),1)
        ! NB: always call with imode=0 to avoid call to addcore in dv_of_drho for 
        !     nlcc pseudo. The call is not needed since we are not moving atoms!!
        !
        call dv_of_drho (dvscfout(1,1,ipert), .false.)
     enddo
     !
     !
     ! ... On output in dvscfin we have the mixed potential
     !
     call mix_potential (2*npe*dfftp%nnr*nspin, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
     !WRITE(mpime+1000, '(1i5,es10.3,1l1,1i5)') my_pool_id, dr2, convt, iter
     !
     ! check that convergent have been reached on ALL processors in this image
     CALL check_all_convt(convt)

     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin
              !call cinterpolate (dvscfin(1,is,ipert), dvscfins(1,is,ipert), -1)
              call fft_interpolate (dffts, drhoscf(:,is,ipert), dfftp, drhoscfh(:,is,ipert))
           enddo
        enddo
     endif
     !
     !
#ifdef __MPI
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     call mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
     tcpu = get_clock ('KC_WANN')

     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npe
     WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     ! 
     FLUSH( stdout )
     !
     if (check_stop_now()) call stop_smoothly_ph (.false.)
     if (convt) goto 155
     !
  enddo  ! loop over iteration
  !
155 iter0=0
  !
  ipert =1  
  ! The density variation in G-space
  !
  DO is = 1, nspin
     !aux(:) = drhoscf(:,is,ipert)
     aux(:) = drhoscfh(:,is,ipert)
     CALL fwfft ('Rho', aux, dffts)
     drhog_scf(:,is) = aux(dffts%nl(:))
  ENDDO
  !
  IF(present(drhor_scf)) drhor_scf = drhoscf(:,:,ipert)
  !
  ! The induced density in G space
  !
  if (lmetq0) deallocate (ldoss)
  if (lmetq0) deallocate (ldos)
  deallocate (h_diag)
  deallocate (aux)
  deallocate (aux1)
  deallocate (aux2)
  deallocate (drhoc)
  deallocate (dbecsum)
  deallocate (drhoscf )
  deallocate (dvscfout)
  deallocate (drhoscfh)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)

  call stop_clock ('solve_linter')
  return
end subroutine solve_linter_koop

END MODULE
