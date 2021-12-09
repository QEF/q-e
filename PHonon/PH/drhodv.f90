!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine drhodv (nu_i0, nper, drhoscf)
  !-----------------------------------------------------------------------
  !! This subroutine computes the electronic term
  !! \( 2\langle d\psi|dv - e ds|\psi\rangle \)
  !! of the dynamical matrix - Eq. (B35) of PRB 64, 235118 (2001).
  !! The contribution of the nonlocal potential is calculated in
  !! \(\texttt{drhodvnl}\), the contribution of the local potential
  !! in \(\texttt{drhodvloc}\).  
  !! Note that \(\text{drhoscf}\) contain only the smooth part of the
  !! induced charge density, calculated in solve linter.  
  !! February 2020: the routine has been generalized to 
  !! address also the case noncolin.and.domag. For the 
  !! theoretical background please refer to: 
  !! Phys. Rev. B 100, 045115 (2019).
  !
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE klist,     ONLY : xk, ngk, igk_k
  USE gvect,     ONLY : g
  USE cell_base, ONLY : tpiba
  USE lsda_mod,  ONLY : current_spin, lsda, isk, nspin
  USE wvfct,     ONLY : npwx, nbnd
  USE uspp,      ONLY : nkb, vkb, deeq_nc, okvan
  USE becmod,    ONLY : calbec, bec_type, becscal, allocate_bec_type, &
                        deallocate_bec_type
  USE fft_base,  ONLY : dfftp
  USE io_global, ONLY : stdout
  USE buffers,   ONLY : get_buffer
  USE noncollin_module, ONLY : noncolin, domag, npol, nspin_mag

  USE dynmat,   ONLY : dyn, dyn_rec
  USE modes,    ONLY : u
  USE units_lr, ONLY : lrdwf, iudwf

  USE eqv,      ONLY : dpsi
  USE qpoint,   ONLY : nksq, ikks, ikqs
  USE control_lr, ONLY : lgamma
  USE lrus,     ONLY : becp1
  USE phus,     ONLY : alphap, int1_nc
  USE qpoint_aux, ONLY : becpt, alphapt
  USE nc_mag_aux, ONLY : int1_nc_save, deeq_nc_save

  USE mp_pools,         ONLY : inter_pool_comm
  USE mp,               ONLY : mp_sum
  USE uspp_init,        ONLY : init_us_2

  implicit none

  integer :: nper
  !! input: number of perturbations of this represent
  integer :: nu_i0
  !! input: the initial position of the mode
  complex(DP) :: drhoscf(dfftp%nnr,nspin_mag,nper)
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

  TYPE (bec_type), POINTER :: dbecq(:), dalpq(:,:)
  !
  !   Initialize the auxiliary matrix wdyn
  !
  call start_clock ('drhodv')

  ALLOCATE (dbecq(nper))
  ALLOCATE (dalpq(3,nper))
  DO ipert=1,nper
     call allocate_bec_type ( nkb, nbnd, dbecq(ipert) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbnd, dalpq(ipol,ipert) )
     ENDDO
  END DO
  allocate (aux   ( npwx*npol , nbnd))
  dynwrk(:,:) = (0.d0, 0.d0)
  wdyn  (:,:) = (0.d0, 0.d0)
  !
  !   We need a sum over all k points ...
  !
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2
  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     if (lsda) current_spin = isk (ikk)
     call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     DO isolv=1, nsolv
        do mu = 1, nper
           nrec = (mu - 1) * nksq + ik + (isolv-1) * nksq * nper
           if (nksq > 1 .or. nper > 1 .OR. nsolv==2) &
                             call get_buffer(dpsi, lrdwf, iudwf, nrec)
           call calbec (npwq, vkb, dpsi, dbecq(mu) )
           do ipol = 1, 3
              aux=(0.d0,0.d0)
              do ibnd = 1, nbnd
                 do ig = 1, npwq
                    aux (ig, ibnd) = dpsi (ig, ibnd) * &
                         (xk (ipol, ikq) + g (ipol, igk_k(ig,ikq) ) )
                 enddo
                 if (noncolin) then
                    do ig = 1, npwq
                       aux (ig+npwx, ibnd) = dpsi (ig+npwx, ibnd) * &
                            (xk (ipol, ikq) + g (ipol, igk_k(ig,ikq) ) )
                    enddo
                 endif
              enddo
              call calbec (npwq, vkb, aux, dalpq(ipol,mu) )
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
              int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,2)
           ENDIF
           call drhodvnl (ik, ikk, nper, nu_i0, dynwrk, becpt, alphapt, &
                                                         dbecq, dalpq)
           IF (okvan) THEN
              deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
              int1_nc(:,:,:,:,:)=int1_nc_save(:,:,:,:,:,1)
           ENDIF
        ENDIF
     ENDDO
  enddo
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
   call drhodvloc (nu_i0, nper, drhoscf, wdyn)
  !
  ! add to the rest of the dynamical matrix
  !
  !      WRITE( stdout,*) 'drhodv dyn, wdyn'
  !      call tra_write_matrix('drhodv dyn',dyn,u,nat)
  !      call tra_write_matrix('drhodv wdyn',wdyn,u,nat)

  dyn (:,:) = dyn (:,:) + wdyn (:,:)
  dyn_rec(:,:) = dyn_rec(:,:) + wdyn(:,:)

  deallocate (aux)

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
