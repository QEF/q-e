!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine project (ipol)
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is READ from file if this_pcxpsi_is_on_file(ik,ipol)=.true.
  ! otherwise dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set)
  !
  USE io_files,        ONLY :nwordwfc
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba2
  USE io_global,       ONLY : stdout,ionode
  USE klist,           ONLY : xk
  USE gvect,           ONLY : g,gstart
  USE wvfct,           ONLY : npw, npwx, nbnd, g2kin, et
  USE klist,           ONLY : igk_k
  USE wavefunctions_module, ONLY: evc

  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : bec_type, becp, calbec, &
                              allocate_bec_type, deallocate_bec_type
  USE uspp,            ONLY : okvan, nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
  USE ramanm,          ONLY : eth_rps
  USE eqv,             ONLY : dpsi, dvpsi
  USE lrus,            ONLY : becp1
  USE qpoint,          ONLY : nksq
  USE units_ph,        ONLY : this_pcxpsi_is_on_file, lrcom, iucom, &
                              lrebar, iuebar
  USE control_lr,      ONLY : nbnd_occ,alpha_pv
  use mp, ONLY: mp_sum,mp_min,mp_max
  USE mp_global,     ONLY : inter_pool_comm,intra_pool_comm 
  USE eqv,                  ONLY : evq
  use hartree_mod, only :init_linear
  
  implicit none
  !
  real(DP):: emin,emax
  integer, intent(IN) :: ipol
  integer :: ik,iun
  integer, external   :: find_free_unit
  real(DP),allocatable :: eprec(:,:) 
  !
  ! Local variables
  !
  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0,  &
             nrec, is, js, ijs, kbnd,ipw
  ! counters
  
  type (bec_type) ::becp0
  real(DP), allocatable  :: h_diag (:,:)
  ! the diagonal part of h_scf
  type(bec_type) :: becp2 ! the scalar products
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter
  logical :: conv_root
  ! true if convergence has been achieved
  COMPLEX(DP), EXTERNAL :: zdotc
  real(DP), EXTERNAL ::ddot
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:),work(:,:)
  real(DP) ::emme(nbnd,nbnd)  
  REAL(DP), ALLOCATABLE  :: eprec_real(:)
  logical ::l_test,exst
  character(len=256) ::prefix_restart_cg,c_ipol
 
  external ch_psi_all, cg_psi
  !debug
  allocate(evq(npwx,nbnd))
  evq=evc
  !

  call start_clock ('project')
  prefix_restart_cg='cg_restart'
  ik=1
  write(c_ipol,'(I0)') ipol
  prefix_restart_cg=trim(prefix_restart_cg)//'.'//trim(c_ipol)
  ALLOCATE( aux1( npwx, nbnd ) )
  allocate( work(npwx,nbnd))
  dpsi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  !
  call allocate_bec_type ( nkb, nbnd, becp2) 

  ! calculate the commutator [H,x_ipol]  psi > and store it in dpsi
  ! dvpsi used as workspace
!aggiunto
  call allocate_bec_type ( nkb, nbnd, becp0)
  CALL calbec (npw, vkb, evc, becp0) 

  call commutator_Hx_psi (ik, nbnd, becp0, becp2, ipol, dpsi, dvpsi )
  call deallocate_bec_type(becp0)
  !
  !    orthogonalize dpsi to the valence subspace: ps = <evc|dpsi>
  !    Apply -P^+_c
  !    NB it uses dvpsi as workspace
  !
  l_test=.true.
  !ortogonalizzazione manuale
  if (l_test) then  
     emme=0.d0
     call dgemm ('T','N',nbnd,nbnd,2*npw,2.d0,evc,2*npwx,dpsi,2*npwx,0.d0,emme,nbnd)
     if (gstart==2) then
         do ibnd=1,nbnd
            do jbnd=1,nbnd
               emme(ibnd,jbnd)=emme(ibnd,jbnd)-dble(conjg(evc(1,ibnd))*dpsi(1,jbnd))
            end do
         end do
     end if
     call mp_sum(emme, intra_pool_comm)
!     do ibnd=1,nbnd
!        do jbnd=1,nbnd
!            dpsi(1:npw,ibnd)=dpsi(1:npw,ibnd)-evc(1:npw,jbnd)*emme(jbnd,ibnd)
!        end do
!     end do
      call dgemm ('N','N',2*npw,nbnd,nbnd,-1.d0,evc,2*npwx,emme,nbnd,1.d0,dpsi,2*npwx) 
  end if
  !dpsi=-dpsi
  !
  !   dpsi contains P^+_c [H-eS,x] psi_v for the three crystal polarizations
  !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
  !
  
  do ig = 1, npw
!     g2kin (ig) = SUM((xk(1:3,ik) +g (1:3, igk (ig)) ) **2) *tpiba2
       g2kin (ig) = SUM((g (1:3, igk_k (ig,1)) ) **2) *tpiba2
  enddo

  allocate (h_diag( npwx*npol, nbnd))
  h_diag=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!eprec stuff (al momento eprec=1, meglio non riesco a farlo andare)

  allocate(eprec(nbnd,1))
  allocate(eprec_real(nbnd))
!  aux1=(0.d0,0.d0)
!  DO ig = 1, npw
!     aux1 (ig,1:nbnd) = g2kin (ig) * evc (ig, 1:nbnd)
!  END DO
!  DO ibnd=1,nbnd
!       work =0.d0
!       DO ig=1,npw
!          work(ig,ibnd)=g2kin(ig)*evc(ig,ibnd)
!       ENDDO

!        eprec_real(ibnd)=2.0d0*DDOT(2*npw,evc(1,ibnd),1,work(1,ibnd),1)
          !
!          IF(gstart==2) THEN
!             eprec_real(ibnd)=eprec_real(ibnd)-DBLE(evc(1,ibnd))*DBLE(work(1,ibnd))
!          ENDIF
          !
!          eprec_real(ibnd)=1.35d0*eprec_real(ibnd)

!      eprec (ibnd,1) = 1.35d0 * zdotc(npwx*npol,evc(1,ibnd),1,aux1(1,ibnd),1)
!      eprec (ibnd,1)=eprec (ibnd,ik) * 2.d0
!      if (gstart==2) then
!         eprec(ibnd,ik)=eprec(ibnd,ik)-1.35d0*conjg(evc(1,ibnd))*aux1(1,ibnd)
!      end if

!#ifdef __MPI
!        CALL mp_sum ( eprec_real, intra_pool_comm )
!#endif
!  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! preconditioning
  do ibnd = 1, nbnd
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) )
     enddo
  enddo
  deallocate(eprec_real)
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calcolo di alpha_pv

  emin = et (1, 1)
 
  do ibnd = 1, nbnd
     emin = min (emin, et (ibnd, 1) )
  enddo

#ifdef __MPI
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif
  emax = et (1, 1)
  do ibnd = 1, nbnd
      emax = max (emax, et (ibnd, 1) )
  enddo
#ifdef __MPI
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
  alpha_pv = 2.d0 * (emax - emin)

!
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!fine conto di alpha_pv
!
!!!!!!!!!!!!!!!!!!!!!!!! inizializzazione partenza per cg

  if (trim(init_linear)=='restart') then
     iun=find_free_unit()
     call diropn_due (prefix_restart_cg,iun, 'save', 2*nwordwfc, exst )
     print*,'PROJECT,DIROPN',exst,prefix_restart_cg,trim(c_ipol)
     call davcio (dvpsi, 2*nwordwfc, iun, 1, - 1)
     close(iun)
  else 
      dvpsi(:,:)=(0.d0,0.d0)         
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  do ibnd = 1, nbnd
!     do ig = 1, npw
!        h_diag (ig, ibnd) = 1.d0 / max(1.d0,(g2kin(ig)-et(ibnd,1)+alpha_pv )/eprec(ibnd,ik))
!     enddo
!  enddo

  !
  eth_rps      = 1.D-10
  thresh = eth_rps
  !
  call cgsolve_all (ch_psi_all, cg_psi, et (1, 1), dpsi, dvpsi, &
       h_diag, npwx, npw, thresh, 1, lter, conv_root, anorm, &
       nbnd, npol)
  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm
  !CALL flush_unit( stdout )

!!!!!!!!!! scriviamo la soluzione su disco per un restart successivo
  if ((trim(init_linear)=='restart').or.(trim(init_linear)=='scratch')) then
      iun=find_free_unit()
      call diropn_due (prefix_restart_cg,iun, 'save', 2*nwordwfc, exst )
      print*,'PROJECT,DIROPN',exst,prefix_restart_cg,trim(c_ipol)
      call davcio (dvpsi, 2*nwordwfc, iun, 1, + 1)
      close(iun)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  deallocate (h_diag)
  deallocate (eprec)
  deallocate(aux1)
  deallocate(work)
  !
  ! we have now obtained P_c x |psi>.
  ! In the case of USPP this quantity is needed for the Born
  ! effective charges, so we save it to disc
  !
  ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
  ! therefore we apply S again, and then subtract the additional term
  ! furthermore we add the term due to dipole of the augmentation charges.
  !
  IF (nkb > 0) call deallocate_bec_type (becp2)

  deallocate(evq)

  call stop_clock ('project')
  return
end subroutine project
