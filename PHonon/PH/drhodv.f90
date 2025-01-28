!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine drhodv (nu_i0, nper, drhos)
  !-----------------------------------------------------------------------
  !! This subroutine computes the electronic term
  !! \( 2\langle d\psi|dv - e ds|\psi\rangle \)
  !! of the dynamical matrix - Eq. (B35) of PRB 64, 235118 (2001).
  !! The contribution of the nonlocal potential is calculated in
  !! \(\texttt{drhodvnl}\), the contribution of the local potential
  !! in \(\texttt{drhodvloc}\).  
  !! Note that \(\text{drhos}\) contain only the smooth part of the
  !! induced charge density, calculated in solve linter.  
  !! February 2020: the routine has been generalized to 
  !! address also the case noncolin.and.domag. For the 
  !! theoretical background please refer to: 
  !! Phys. Rev. B 100, 045115 (2019).
  !
  !civn: at the time I am doing this OpenACC hardly manages offloading arrays of data structures.
  !      probably in future we can remove bectmp and directly use dbecq(mu) and dalpq(ipol,mu)
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE klist,     ONLY : xk, ngk, igk_k
  USE gvect,     ONLY : g
  USE cell_base, ONLY : tpiba
  USE lsda_mod,  ONLY : current_spin, lsda, isk, nspin
  USE wvfct,     ONLY : npwx, nbnd
  USE uspp,      ONLY : nkb, vkb, deeq_nc, okvan
  USE becmod,    ONLY : calbec, bec_type, becscal, becupdate, &
                        allocate_bec_type, deallocate_bec_type, &
                        allocate_bec_type_acc, deallocate_bec_type_acc
  USE fft_base,  ONLY : dffts
  USE io_global, ONLY : stdout
  USE buffers,   ONLY : get_buffer
  USE noncollin_module, ONLY : noncolin, domag, npol, nspin_mag

  USE dynmat,   ONLY : dyn, dyn_rec
  USE modes,    ONLY : u
  USE units_lr, ONLY : lrdwf, iudwf

  USE eqv,      ONLY : dpsi
  USE qpoint,   ONLY : nksq, ikks, ikqs
  USE lrus,     ONLY : becp1
  USE phus,     ONLY : alphap, int1_nc
  USE qpoint_aux, ONLY : becpt, alphapt
  USE lr_nc_mag,  ONLY : int1_nc_save, deeq_nc_save

  USE mp_pools,         ONLY : inter_pool_comm
  USE mp,               ONLY : mp_sum
  USE uspp_init,        ONLY : init_us_2
  USE control_flags,    ONLY : offload_type
  USE control_lr,       ONLY : lmultipole

  implicit none

  integer :: nper
  !! input: number of perturbations of this represent
  integer :: nu_i0
  !! input: the initial position of the mode
  complex(DP) :: drhos(dffts%nnr,nspin_mag,nper)
  !! the change of density due to perturbations
  !
  ! ... local variables
  !
  integer :: mu, ik, ikq, ig, nu_i, nu_j, na_jcart, ibnd, nrec, &
       ipol, ikk, ipert, npw, npwq, isolv, nsolv
  ! counters
  ! ikk: record position for wfc at k

  complex(DP) :: fact, ps, dynwrk (3 * nat, 3 * nat), &
       wdyn (3 * nat, 3 * nat)
  complex(DP), allocatable ::  aux (:,:)
  ! work space

#if defined(__CUDA)
  TYPE (bec_type) :: bectmp
#endif
  TYPE (bec_type), ALLOCATABLE :: dbecq(:), dalpq(:,:)
  !
  !   Initialize the auxiliary matrix wdyn
  !
  call start_clock ('drhodv')
#if defined(__CUDA)
  Call allocate_bec_type_acc( nkb, nbnd, bectmp ) 
#endif
  ALLOCATE (dbecq(nper))
  ALLOCATE (dalpq(3,nper))
  DO ipert=1,nper
     call allocate_bec_type ( nkb, nbnd, dbecq(ipert) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbnd, dalpq(ipol,ipert) )
     ENDDO
  END DO
  allocate (aux   ( npwx*npol , nbnd))
  !
  dynwrk(:,:) = (0.d0, 0.d0)
  wdyn  (:,:) = (0.d0, 0.d0)
  !
  !   We need a sum over all k points ...
  !
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2
  !$acc data copyin(xk, dpsi) create(aux( 1:npwx*npol , 1:nbnd)) 
  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     if (lsda) current_spin = isk (ikk)
     call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb, .true.)
     DO isolv=1, nsolv
        do mu = 1, nper
           !
           IF (.NOT. lmultipole) THEN
             nrec = (mu - 1) * nksq + ik + (isolv-1) * nksq * nper
             if (nksq > 1 .or. nper > 1 .OR. nsolv==2) then
                               call get_buffer(dpsi, lrdwf, iudwf, nrec)
                               !$acc update device(dpsi)
             endif
           ELSE
             nrec = ik
             call get_buffer(dpsi, lrdwf, iudwf, nrec)
             !$acc update device(dpsi)
           ENDIF
           !
#if defined(__CUDA)
           call calbec( offload_type, npwq, vkb, dpsi, bectmp )
           call becupdate( offload_type, dbecq, mu, nper, bectmp )
#else
           call calbec( offload_type, npwq, vkb, dpsi, dbecq(mu) )
#endif
           do ipol = 1, 3
#if defined(__CUDA)
              !$acc parallel loop collapse(2)
              do ibnd = 1, nbnd
                 do ig = 1, npwq
                    aux (ig,ibnd) = (0.d0,0.d0)
                 end do
              end do
#else
              aux=(0.d0,0.d0)
#endif
              !$acc parallel loop collapse(2) 
              do ibnd = 1, nbnd
                 do ig = 1, npwq
                    aux (ig, ibnd) = dpsi (ig, ibnd) * &
                         (xk (ipol, ikq) + g (ipol, igk_k(ig,ikq) ) )
                 enddo
              end do
              if (noncolin) then
                 !$acc parallel loop collapse(2)
                 do ibnd = 1, nbnd
                    do ig = 1, npwq
                       aux (ig+npwx, ibnd) = dpsi (ig+npwx, ibnd) * &
                            (xk (ipol, ikq) + g (ipol, igk_k(ig,ikq) ) )
                    enddo
                 enddo
              endif
#if defined(__CUDA)
              call calbec( offload_type, npwq, vkb, aux, bectmp )
              call becupdate( offload_type, dalpq, ipol, 3, mu, nper, bectmp )
#else
              call calbec( offload_type, npwq, vkb, aux, dalpq(ipol,mu) )
#endif
           enddo

        enddo
        fact = CMPLX(0.d0, tpiba,kind=DP)
        DO ipert=1,nper
           DO ipol=1,3
              CALL becscal(fact,dalpq(ipol,ipert),nkb,nbnd)
           ENDDO
        ENDDO
        IF (isolv==1) THEN
           call drhodvnl (ik, ikk, nper, nu_i0, dynwrk, becp1, alphap, &
                                                              dbecq, dalpq)
        ELSE
           IF (okvan) THEN
              deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,2)
              !$acc update device(deeq_nc)
              int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,2)
           ENDIF
           call drhodvnl (ik, ikk, nper, nu_i0, dynwrk, becpt, alphapt, &
                                                         dbecq, dalpq)
           IF (okvan) THEN
              deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
              !$acc update device(deeq_nc)
              int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,1)
           ENDIF
        ENDIF
     ENDDO
  enddo
  !$acc end data
  !
  !   put in the basis of the modes
  !
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        ps = (0.0d0, 0.0d0)
        do na_jcart = 1, 3 * nat
           ps = ps + dynwrk (nu_i, na_jcart) * u (na_jcart, nu_j)
        enddo
        wdyn (nu_i, nu_j) = wdyn (nu_i, nu_j) + ps
     enddo

  enddo
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call mp_sum ( wdyn, inter_pool_comm )
  IF (nsolv==2) wdyn=wdyn/2.0_DP
  !
  ! add the contribution of the local part of the perturbation
  !
   call drhodvloc (nu_i0, nper, drhos, wdyn)
  !
  ! add to the rest of the dynamical matrix
  !
  !      WRITE( stdout,*) 'drhodv dyn, wdyn'
  !      call tra_write_matrix('drhodv dyn',dyn,u,nat)
  !      call tra_write_matrix('drhodv wdyn',wdyn,u,nat)

  dyn (:,:) = dyn (:,:) + wdyn (:,:)
  dyn_rec(:,:) = dyn_rec(:,:) + wdyn(:,:)

  deallocate (aux)
#if defined(__CUDA)
  call deallocate_bec_type_acc( bectmp )
#endif
  do ipert=1,nper
     do ipol=1,3
        call deallocate_bec_type ( dalpq(ipol,ipert) )
     enddo
  end do
  deallocate (dalpq)
  do ipert=1,nper
     call deallocate_bec_type ( dbecq(ipert) )
  end do
  deallocate(dbecq)

  call stop_clock ('drhodv')
  return
end subroutine drhodv
