!!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!----------------------------------------------------------------------
!
! OBM :
!  090709: Complete overhaul, trying to mimick what Brent Walker has done to 
!          obsolote dv_psi_e of PH code. 
!
subroutine lr_dvpsi_e(ik,ipol,dvpsi)
!----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !
  ! On output: dvpsi contains P_c^+ x | psi_ik > in crystal axis 
  !            (projected on at(*,ipol) )
  !
  ! dvpsi is COMPUTED and WRITTEN on file (vkb,evc,igk must be set) !OBM: This is now handled elesewhere
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba, at
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk
  USE gvect,           ONLY : g, gstart
  USE wvfct,           ONLY : npw, npwx, nbnd, igk, g2kin, et
  USE wavefunctions_module, ONLY: evc
  USE lsda_mod,        ONLY : current_spin
  USE spin_orb,        ONLY : lspinorb
  USE noncollin_module,ONLY : noncolin, npol
  USE uspp,            ONLY : okvan, nkb, vkb, qq, qq_so, deeq, deeq_nc
  USE uspp_param,      ONLY : nh
  use control_flags,   only : gamma_only
  !USE ramanm,          ONLY : eth_rps
  !USE eqv,             ONLY : d0psi, eprec
  !USE phus,            ONLY : becp1, becp1_nc !becp1 is calculated here for every k point
  !USE qpoint,          ONLY : npwq, nksq
  !USE units_ph,        ONLY : this_pcxpsi_is_on_file, lrcom, iucom, &
  !                            lrebar, iuebar !related to removed reading part
  USE control_ph,      ONLY : nbnd_occ

  USE mp_global,       ONLY: intra_pool_comm
  USE mp,              ONLY: mp_sum

  USE realus,                ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma,npw_k, real_space_debug
 
   USE lr_variables,   ONLY : lr_verbosity, evc0
   USE io_global,      ONLY : stdout
!DEBUG
  USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,check_vector_gamma,check_vector_f
  !
  !
  implicit none
  !
  integer, intent(IN) :: ipol, ik 
  !
  complex(kind=dp),intent(out) :: dvpsi(npwx,nbnd) 
  real(kind=dp) :: atnorm

  !

  !
  ! Local variables
  !  !
  complex(kind=dp),allocatable :: d0psi(:,:) 
  !

  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0, nrec
  ! counters

  real(DP), allocatable  :: gk (:,:), h_diag (:,:), eprec(:)

  ! the derivative of |k+G|
  real(DP) ::   anorm, thresh
  ! preconditioning cut-off
  ! the desired convergence of linter

  logical :: conv_root
  ! true if convergence has been achieved

  complex(DP), allocatable :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),  &
       work (:,:), spsi(:,:), psc(:,:,:,:), aux(:)
  real(kind=dp), external :: ddot
  complex(DP), external :: ZDOTC
  ! the scalar products
  external lr_ch_psi_all, cg_psi
  !
  !obm debug
  real(DP) ::obm_debug 
  !
  call start_clock ('lr_dvpsi_e') 

  If (lr_verbosity > 5) WRITE(stdout,'("<lr_dvpsi_e>")')

  allocate(d0psi(npwx*npol,nbnd))
  d0psi=(0.d0, 0.d0)
  dvpsi=(0.d0, 0.d0)
  !if (this_pcxpsi_is_on_file(ik,ipol)) then
  !   nrec = (ipol - 1)*nksq + ik
  !   call davcio(dvpsi, lrebar, iuebar, nrec, -1)
  !   call stop_clock ('dvpsi_e')
  !   return
  !end if
  !
  allocate (work ( npwx, MAX(nkb,1)))
  allocate (aux ( npwx*npol ))
  allocate (gk ( 3, npwx))    
  allocate (h_diag( npwx*npol, nbnd))
  !OBM!!!! eprec is also calculated on the fly for each k point     
  allocate(eprec(nbnd))
  !OBM!!!!
  evc(:,:)=evc0(:,:,ik) 
  if (nkb > 0) then
     allocate (dvkb (npwx, nkb), dvkb1(npwx, nkb))
     dvkb (:,:) = (0.d0, 0.d0)
     dvkb1(:,:) = (0.d0, 0.d0)
  end if
  do ig = 1, npw_k(ik)
     gk (1, ig) = (xk (1, ik) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, ik) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, ik) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  enddo 
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e g2kin:",F15.8)') SUM(g2kin(:))
  endif


  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw_k(ik)
        d0psi (ig, ibnd) = (at(1, ipol) * gk(1, ig) + &
             at(2, ipol) * gk(2, ig) + &
             at(3, ipol) * gk(3, ig) ) &
             *(0.d0,-2.d0)*evc (ig, ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, npw_k(ik)
           d0psi (ig+npwx, ibnd) = (at(1, ipol) * gk(1, ig) + &
                at(2, ipol) * gk(2, ig) + &
                at(3, ipol) * gk(3, ig) ) &
                 *(0.d0,-2.d0)*evc (ig+npwx, ibnd)
        end do
     END IF
  enddo
  !!OBM debug
  !      obm_debug=0
  !     do ibnd=1,nbnd
  !        !
  !        obm_debug=obm_debug+ZDOTC(npwx*npol,d0psi(:,ibnd),1,d0psi(:,ibnd),1)
  !        !
  !     enddo
  !     print *, "lr_dvpsi_e d0psi kinetic contribution", obm_debug
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e d0psi kinetic contribution:")')
       do ibnd=1,nbnd
          call check_vector_gamma(d0psi(:,ibnd))
       enddo
  endif
  !!obm_debug



!
! Uncomment this goto and the continue below to calculate 
! the matrix elements of p without the commutator with the
! nonlocal potential.
!
!  goto 111
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  call gen_us_dj (ik, dvkb)
  call gen_us_dy (ik, at (1, ipol), dvkb1) 
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e dvkb:")') 
       jkb=0
      do nt = 1, ntyp
       do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
                call check_vector_f(dvkb(:,jkb))
           enddo
        endif
       enddo
      enddo
  endif
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e dvkb1:")') 
       jkb=0
      do nt = 1, ntyp
       do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
                call check_vector_f(dvkb1(:,jkb))
           enddo
        endif
       enddo
      enddo
  endif



  do ig = 1, npw_k(ik)
     if (g2kin (ig) < 1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     else
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     endif
  enddo

  jkb = 0
  work=(0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
              do ig = 1, npw_k(ik)
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, ipol) * gk (1, ig) + &
                       at (2, ipol) * gk (2, ig) + &
                       at (3, ipol) * gk (3, ig) )
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate (gk)
  !OBM!!!be careful, from bwalker, why?!!!!! 
  work(:,:)=(0.0d0,1.0d0)*work(:,:)
  !OBM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e non-local contribution:")')
       jkb=0
      do nt = 1, ntyp
       do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
                call check_vector_gamma(work(:,jkb))
           enddo
        endif
       enddo
      enddo
  endif


  if(gamma_only) then
     call lr_dvpsi_e_gamma()
  else if (noncolin) then
     call lr_dvpsi_e_noncolin()
  else
     call lr_dvpsi_e_k()
  endif
  !!OBM debug
  !      obm_debug=0
  !     do ibnd=1,nbnd
  !        !
  !        obm_debug=obm_debug+ZDOTC(npwx*npol,dvpsi(:,ibnd),1,dvpsi(:,ibnd),1)
  !        !
  !     enddo
  !     print *, "lr_dvpsi_e norm of dvpsi being returned=", obm_debug
  ! 
  if (lr_verbosity > 10 ) then
       write(stdout,'("lr_dvpsi_e d0psi after lr_dvpsi_e case specific calc:")')
       do ibnd=1,nbnd
          call check_vector_gamma(d0psi(:,ibnd))
       enddo
       write(stdout,'("lr_dvpsi_e dvpsii after lr_dvpsi_e case specific calc:")')
       do ibnd=1,nbnd
          call check_vector_gamma(dvpsi(:,ibnd))
       enddo
  endif

  !!obm_debug

  IF (nkb > 0) THEN
     deallocate (dvkb1, dvkb)
  END IF
  !OBM!!!! End of Brent's seperation
  
  deallocate (h_diag)
  deallocate (work)
  deallocate (aux)
  deallocate (eprec)
  deallocate (d0psi)
  
  !OBM!!!! Addendum to PH dvpsi
  if (okvan) then
    allocate (spsi ( npwx*npol, nbnd))    
    call sm1_psi(.true.,ik,npwx,npw_k(ik),nbnd,dvpsi,spsi)
    dvpsi(:,:) = spsi(:,:)
    deallocate(spsi)  
    !!OBM debug
    If (lr_verbosity > 7) then 
        obm_debug=0
       do ibnd=1,nbnd
          !
          obm_debug=obm_debug+ZDOTC(npwx*npol,dvpsi(:,ibnd),1,dvpsi(:,ibnd),1)
          !
       enddo
#ifdef __PARA 
       call mp_sum(obm_debug, intra_pool_comm)
#endif 
       WRITE(stdout,'("lr_dvpsi_e: dvpsi after sm1_psi:",E15.5)') obm_debug
    endif 
  endif
  
  ! 
  ! For some ibrav the crystal axes are not normalized
  ! Here we include the correct normalization
  ! for Lanczos initial wfcs
  atnorm=dsqrt(at(1,ipol)**2+at(2,ipol)**2+at(3,ipol)**2)
  !
  dvpsi(:,:)=dvpsi(:,:)/atnorm
  !
  !nrec = (ipol - 1)*nksq + ik
  !call davcio(dvpsi, lrebar, iuebar, nrec, 1)
  !this_pcxpsi_is_on_file(ik,ipol) = .true.
  call stop_clock ('lr_dvpsi_e')
  !
  return
  !
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!k-point modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_dvpsi_e_k()
  !OBM!! The k point part is essentially the same as the phonon dvpsi_e, with noncolin part removed
  USE becmod,          ONLY : bec_type, becp, calbec
  
  implicit none
  
  complex(kind=dp), allocatable :: becp1(:,:),becp2(:,:)  
  
  IF (nkb > 0) THEN
        allocate (becp1 (nkb, nbnd))
        allocate (becp2 (nkb, nbnd))
  END IF
     becp1=0.0
     becp2=0.0
     becp%k=0.0
     
     ! OBM : Change calbec to use calbec_type as soon as possible
     !
     call calbec (npw_k(ik), vkb, evc, becp1)
     
     if(nkb>0) then !If there are no beta functions, below loop returns an error
       call calbec (npw_k(ik), work, evc, becp2)
    
     ijkb0 = 0
     allocate (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
     do nt = 1, ntyp
        do na = 1, nat
           if (nt == ityp (na)) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd_occ (ik)
                         !OBM!!! Notice that in the original version ikb is 1, I 
                         !       think this was a bug of the original code
                          ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                              (deeq(ih,jh,na,current_spin) &
                              -et(ibnd,ik)*qq(ih,jh,nt))
                          ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd) * &
                              (-1.d0,0.d0)*(deeq(ih,jh,na,current_spin)&
                              -et(ibnd,ik)*qq(ih,jh,nt)) 
                         !original from Phonon
                         !ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         ! ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         !OBM!!! But why?
                    enddo
                 enddo
              enddo
              ijkb0=ijkb0+nh(nt)
           end if
        end do
     end do
     if (ikb /= nkb .OR. jkb /= nkb) call errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the 
     !       above loop as ZAXPY
     ! 
     
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik), nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik), nkb, &
             (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        deallocate (ps2)
     endif
     !print *, "new version k point, before ortho", d0psi(1:3,1)
     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand
     
     IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
     CALL s_psi (npwx, npw_k(ik), nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi,.false.)
     !d0psi=-d0psi
     
    !print *, "new version k point, after ortho", d0psi(1:3,1)
    
     
     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps 
     thresh = 1.d-5
     h_diag=0.d0
     
     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        deallocate (work)
        allocate (work(npwx,nbnd))
        do ibnd=1,nbnd 
          ! 
          work = 0.d0
          ! 
          conv_root = .true. 
          ! 
          do ig=1,npw_k(ik) 
             work(ig,1)=g2kin(ig)*evc(ig,ibnd) 
          enddo
          ! 
          eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1) 
          ! 
       enddo
       ! 
#ifdef __PARA 
       call mp_sum(eprec, intra_pool_comm)
#endif 
       !OBM!!!
       !print *, "eprec", eprec     
       do ibnd = 1, nbnd_occ (ik)
          do ig = 1, npw_k(ik)
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          enddo
          IF (noncolin) THEN
             do ig = 1, npw_k(ik)
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             enddo
          END IF
       enddo
       !
       !print *, h_diag(1:10,1)
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw_k(ik),& 
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       call lr_cgsolve_all (lr_ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw_k(ik), thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !print *, "after cg_solve_all"
       !CALL lr_normalise( dvpsi(:,:), anorm)
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born 
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       if (okvan) then
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          allocate (spsi ( npwx*npol, nbnd))    
          CALL calbec (npw_k(ik), vkb, dvpsi, becp )
          CALL s_psi(npwx,npw_k(ik),nbnd,dvpsi,spsi)
          call DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          deallocate (spsi)
          call lr_adddvepsi_us_k(becp1,becp2,ipol,ik,dvpsi) 
          !call adddvepsi_us(becp2,ipol,ik)
       endif
      
       !OBM!!!! Addendum to PH/dvpsi_e
       ! 
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace 
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd)
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.) !very strange, check: wothout the last false, normconserving pps fail
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization
       IF (nkb > 0) then
             deallocate(becp1)
             deallocate(becp2)
       endif
       !
       !OBM!!!!
  end subroutine lr_dvpsi_e_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!noncolin modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lr_dvpsi_e_noncolin()
  !OBM!! essentially the same as the phonon dvpsi_e, only with non-collinear part
  USE becmod,          ONLY : bec_type, becp, calbec
  
  implicit none
  
  complex(kind=dp), allocatable :: becp1_nc(:,:,:), becp2_nc(:,:,:)
  
  call errore ('lr_dvpsi_e', 'non collinear not implemented', 1) !well, everything is here, but have to understand the formulation, do not have the time
  IF (nkb > 0) THEN
        allocate (becp1_nc (nkb, npol, nbnd))
        allocate (becp2_nc (nkb, npol, nbnd))
  END IF
 
     !OBM!!! This part is seperated as gamma_only and k_point in Brent's version
        !OBM!!becp1 is also calculated
        call calbec (npw_k(ik), vkb, evc, becp1_nc)
        if(nkb>0) call calbec (npw_k(ik), work, evc, becp2_nc)
    
     ijkb0 = 0
        allocate (psc ( nkb, npol, nbnd, 2))
        psc=(0.d0,0.d0)
     do nt = 1, ntyp
        do na = 1, nat
           if (nt == ityp (na)) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd_occ (ik)
                          IF (lspinorb) THEN 
                             psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+(0.d0,-1.d0)* &
                                (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+       &
                                 becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)- &
                                          et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                             psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+(0.d0,-1.d0)*  &
                                (becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+       &
                                 becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)- &
                                          et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                             psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+(0.d0,-1.d0)* &
                                (becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,1,nt) )+      &
                                becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,2)-  &
                                          et(ibnd,ik)* qq_so(ih,jh,2,nt) ) )
                             psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+(0.d0,-1.d0)*  &
                                (becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,3)  &
                                    -et(ibnd,ik)*qq_so(ih,jh,3,nt) )+      &
                                becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4)-  &
                                          et(ibnd,ik)* qq_so(ih,jh,4,nt) ) )
                          ELSE
                             psc(ikb,1,ibnd,1)=psc(ikb,1,ibnd,1)+ (0.d0,-1.d0)* &
                                 ( becp2_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                                -et(ibnd,ik)*qq(ih,jh,nt)) + &
                                   becp2_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                             psc(ikb,2,ibnd,1)=psc(ikb,2,ibnd,1)+ (0.d0,-1.d0)* &
                                 ( becp2_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+  &
                                   becp2_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                             psc(ikb,1,ibnd,2)=psc(ikb,1,ibnd,2)+ (0.d0,-1.d0)* &
                                 ( becp1_nc(jkb,1,ibnd)*(deeq_nc(ih,jh,na,1) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                   becp1_nc(jkb,2,ibnd)*deeq_nc(ih,jh,na,2) )
                             psc(ikb,2,ibnd,2)=psc(ikb,2,ibnd,2)+ (0.d0,-1.d0)* &
                                 ( becp1_nc(jkb,2,ibnd)*(deeq_nc(ih,jh,na,4) &
                                                -et(ibnd,ik)*qq(ih,jh,nt))+ &
                                   becp1_nc(jkb,1,ibnd)*deeq_nc(ih,jh,na,3) )
                          END IF
                    enddo
                 enddo
              enddo
              ijkb0=ijkb0+nh(nt)
           end if
        end do
     end do
     if (ikb /= nkb .OR. jkb /= nkb) call errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the 
     !       above loop as ZAXPY
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik)*npol, nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, psc(1,1,1,1), nkb, (1.d0,0.d0), &
             d0psi, npwx )
        CALL ZGEMM( 'N', 'N', npw_k(ik), nbnd_occ(ik)*npol, nkb, &
             (1.d0,0.d0),work(1,1), npwx, psc(1,1,1,2), nkb, (1.d0,0.d0), &
             d0psi, npwx )
        deallocate (psc)
    
     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand
     IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd )
     CALL s_psi (npwx, npw_k(ik), nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi, .false.)
     !d0psi=-d0psi
     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps 
     thresh = 1.d-5
     h_diag=0.d0
     
     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        deallocate (work)
        allocate (work(npwx,nbnd))
        do ibnd=1,nbnd 
          ! 
          work = 0.d0
          ! 
          conv_root = .true. 
          ! 
          do ig=1,npw_k(ik) 
             work(ig,1)=g2kin(ig)*evc(ig,ibnd) 
          enddo
          ! 
          eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1) 
          ! 
       enddo
       ! 
#ifdef __PARA 
       call mp_sum(eprec, intra_pool_comm)
#endif 
       !OBM!!!
       !print *, "eprec", eprec
       do ibnd = 1, nbnd_occ (ik)
          do ig = 1, npw_k(ik)
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          enddo
             do ig = 1, npw_k(ik)
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             enddo
       enddo
       !
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw_k(ik),& 
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       call lr_cgsolve_all (lr_ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw_k(ik), thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born 
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       if (okvan) then
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          allocate (spsi ( npwx*npol, nbnd))    
          CALL calbec (npw_k(ik), vkb, dvpsi, becp )
          CALL s_psi(npwx,npw_k(ik),nbnd,dvpsi,spsi)
          call DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          deallocate (spsi)
          !OBM!!! non collinear in lr_addd_vespi is still not working
          !call adddvepsi_us(becp2_nc,ipol,ik)
       endif     
       !OBM!!!! Addendum to PH/dv_psi_e
       ! 
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace 
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw_k(ik), vkb, evc, becp, nbnd )
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.)
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization
       IF (nkb > 0) then
             deallocate(becp1_nc)
             deallocate(becp2_nc)
       endif
       !
       !OBM!!!!
  end subroutine lr_dvpsi_e_noncolin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!GAMMA point modifications
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_dvpsi_e_gamma
  !OBM!! modified for gamma point algorithms
  USE becmod,          ONLY : bec_type,becp, calbec
  
  implicit none
  
  real(kind=dp), allocatable :: becp1(:,:),becp2(:,:)  
  
  IF (nkb > 0) THEN
        allocate (becp1 (nkb, nbnd))
        allocate (becp2 (nkb, nbnd))
  END IF
     becp1=0.0
     becp2=0.0
     becp%r=0.0
 
     call calbec (npw, vkb, evc, becp1)
 
     if(nkb>0) then
  
      call calbec (npw, work, evc, becp2)
    
     ijkb0 = 0
     allocate (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
     do nt = 1, ntyp
        do na = 1, nat
           if (nt == ityp (na)) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd_occ (ik)
                         !OBM!!! Notice that in the original version ikb is 1
                          ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                              (deeq(ih,jh,na,current_spin) &
                              -et(ibnd,ik)*qq(ih,jh,nt))
                          ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd) * &
                              (-1.d0,0.d0)*(deeq(ih,jh,na,current_spin)&
                              -et(ibnd,ik)*qq(ih,jh,nt)) 
                         !original from Phonon
                         !ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin) &
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         ! ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                         !     (0.d0,-1.d0)*(deeq(ih,jh,na,current_spin)&
                         !     -et(ibnd,ik)*qq(ih,jh,nt))
                         !OBM!!! But why?
                    enddo
                 enddo
              enddo
              ijkb0=ijkb0+nh(nt)
           end if
        end do
     end do
     if (ikb /= nkb .OR. jkb /= nkb) call errore ('lr_dvpsi_e', 'unexpected error',1)
     !OBM!!! In the archaic version Brent used, this part was embedded in the 
     !       above loop as ZAXPY
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
             (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
             (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
             d0psi(1,1), npwx )
        deallocate (ps2)
     endif
     !
     !    orthogonalize d0psi to the valence subspace: ps = <evc|d0psi>
     !    Apply -P^+_c
     !OBM!! lr_ortho no longer calculates sevc, do it beforhand
     !print *, "new version,gamma before ortho", d0psi(1:3,1)
     IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
     CALL s_psi (npwx, npw, nbnd, evc, dvpsi)
     CALL lr_ortho(d0psi, evc, ik, ik, dvpsi,.false.) !Strange bug, without last .false., norm conserving pps fail
     !d0psi=-d0psi
     !print *, "new version,gamma after ortho", d0psi(1:3,1)
     
     !
     !   d0psi contains P^+_c [H-eS,x] psi_v for the polarization direction ipol
     !   Now solve the linear systems (H-e_vS)*P_c(x*psi_v)=P_c^+ [H-e_vS,x]*psi_v
     !
     !OBM!!! thresh = eth_rps 
     thresh = 1.d-5
     h_diag=0.d0
     
     !OBM!!! on the fly calculation of eprec
     ! I am using the work as a working memory, it seems not to be used anymore
        deallocate (work)
        allocate (work(npwx,nbnd))
        do ibnd=1,nbnd 
          ! 
          work = 0.d0
          ! 
          conv_root = .true. 
          ! 
          do ig=1,npw
             work(ig,1)=g2kin(ig)*evc(ig,ibnd) 
          enddo
          ! 
          eprec(ibnd)=2.0d0*DDOT(2*npw,evc(1,ibnd),1,work,1)                  
          !                                                                    
          if(gstart==2) eprec(ibnd)=eprec(ibnd)-dble(evc(1,ibnd))*dble(work(1,ibnd))
          !                                                                    
          eprec(ibnd)=1.35d0*eprec(ibnd)                                       
          !OBM!!! was : eprec(ibnd)=1.35d0*ZDOTC(npw_k(ik),evc(1,ibnd),1,work,1) 
          ! 
       enddo
       ! 
#ifdef __PARA 
       call mp_sum(eprec, intra_pool_comm)
#endif 
       !print *, eprec
       !OBM!!!
       !print *, "eprec", eprec
       do ibnd = 1, nbnd_occ (ik)
          do ig = 1, npw
             h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
          enddo
          IF (noncolin) THEN
             do ig = 1, npw
                h_diag (ig+npwx, ibnd) = 1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
             enddo
          END IF
       enddo
       !print *, h_diag(1:10,1)
       !
       !OBM!!! upto here, the dvpsi was used as a scratch
       !
       dvpsi(:,:) = (0.d0, 0.d0)
       !
       !OBM!!!! Original was
       !call lr_cgsolve_all(et(1,ik),d0psi,dvpsi,h_diag,npwx,npw,& 
       !            thresh,ik,lter,conv_root,anorm,nbnd_occ (ik),lr_alpha_pv)
       call lr_cgsolve_all (lr_ch_psi_all, cg_psi, et (1, ik), d0psi, dvpsi, &
            h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
            nbnd_occ(ik), 1)
       !
       if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
            & " linter: root not converged ",e10.3)') &
            ik, ibnd, anorm
       !
       CALL flush_unit( stdout )
       !print *, "after cg_solve_all"
       !CALL lr_normalise( dvpsi(:,:), anorm)
       !
       !
       ! we have now obtained P_c x |psi>.
       ! In the case of USPP this quantity is needed for the Born 
       ! effective charges, so we save it to disc
       !
       ! In the US case we obtain P_c x |psi>, but we need P_c^+ x | psi>,
       ! therefore we apply S again, and then subtract the additional term
       ! furthermore we add the term due to dipole of the augmentation charges.
       !
       if (okvan) then
          !
          ! for effective charges
          !
          ! nrec = (ipol - 1) * nksq + ik
          ! call davcio (dvpsi, lrcom, iucom, nrec, 1)
          !
          allocate (spsi ( npwx*npol, nbnd))    
          CALL calbec (npw, vkb, dvpsi, becp )
          CALL s_psi(npwx,npw,nbnd,dvpsi,spsi)
          call DCOPY(2*npwx*npol*nbnd,spsi,1,dvpsi,1)
          deallocate (spsi)
          call lr_adddvepsi_us_gamma(becp1,becp2,ipol,ik,dvpsi) 
          !call adddvepsi_us(becp2,ipol,ik)
       endif
      
       !OBM!!!! Addendum to PH dvpsi_e
       ! 
       work = 0.d0 !Reset working space in any case
       ! orthogonalize dvpsi to the valence subspace 
       !
       !OBM!!!! due to modifications, lr_ortho requires sevc as input, putting it into work
       IF (okvan) CALL calbec ( npw, vkb, evc, becp, nbnd)
       CALL s_psi (npwx, npw_k(ik), nbnd, evc, work)
       CALL lr_ortho(dvpsi, evc, ik, ik, work,.false.)
       !dvpsi=-dvpsi
       !OBM!!!! end of orthogonolization
       
       IF (nkb > 0) then
             deallocate(becp1)
             deallocate(becp2)
       endif
       !
       !OBM!!!!
  end subroutine lr_dvpsi_e_gamma

end subroutine lr_dvpsi_e


