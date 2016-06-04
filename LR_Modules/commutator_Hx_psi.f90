!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine commutator_Hx_psi (ik, nbnd_occ, becp1, becp2, ipol, dpsi)
  !----------------------------------------------------------------------
  !
  ! On output: dpsi contains [H,x_ipol] | psi_ik > in crystal axis 
  !            (projected on at(*,ipol) )
  !
  ! vkb and evc must be properly set for the appropriate k-point
  ! in addition becp1 must be set equal to becp1 = <vkb|evc>
  ! as it is done in PH/phq_init.f90 for the k-point ik
  ! NB: here the last index of becp1 is missing, hence it refers 
  !     to a single k-point
  !
  !    CALL calbec (npw, vkb, evc, becp1(:,:) )
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba, at
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE io_global,       ONLY : stdout
  USE klist,           ONLY : xk, igk_k, ngk
  USE gvect,           ONLY : g
  USE wvfct,           ONLY : npwx, nbnd, et
  USE wavefunctions_module, ONLY: evc
  USE lsda_mod,        ONLY : nspin
  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : becp, bec_type, calbec
  USE uspp,            ONLY : nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
  USE control_flags,   ONLY : gamma_only

  implicit none
  COMPLEX(DP), INTENT(OUT)    :: dpsi(npwx*npol,nbnd)
  TYPE(bec_type), INTENT(IN)  :: becp1 ! dimensions ( nkb, nbnd )
  TYPE(bec_type), INTENT(INOUT) :: becp2 ! dimensions ( nkb, nbnd )
  !
  INTEGER, INTENT(IN) :: ik, nbnd_occ, ipol
  !
  ! Local variables
  !
  integer :: npw, ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0,  &
             nrec, is, js, ijs
  ! counters

  real(DP), allocatable  :: gk (:,:), g2k(:)
  ! the derivative of |k+G|
  complex(DP), allocatable :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),  &
       work (:,:), psc(:,:,:,:), deff_nc(:,:,:,:)
  REAL(DP), allocatable :: deff(:,:,:)
  !

  CALL start_clock ('commutator_Hx_psi')
  dpsi=(0.d0, 0.d0)
  !
  npw = ngk(ik)
  allocate (gk(3, npw), g2k(npw) )    
  do ig = 1, npw
     gk (1:3, ig) = (xk (1:3, ik) + g (1:3, igk_k(ig,ik) ) ) * tpiba
     g2k (ig) = SUM(gk (1:3, ig) **2 )
  enddo
  !
  ! this is  the kinetic contribution to [H,x]:  -2i (k+G)_ipol * psi
  !
  do ibnd = 1, nbnd_occ
     do ig = 1, npw
        dpsi(ig,ibnd) = SUM(at(1:3,ipol)*gk(1:3,ig))*(0.d0,-2.d0)*evc (ig,ibnd)
     enddo
     IF (noncolin) THEN
        do ig = 1, npw
           dpsi (ig+npwx, ibnd) = (at(1, ipol) * gk(1, ig) + &
                at(2, ipol) * gk(2, ig) + &
                at(3, ipol) * gk(3, ig) ) &
                 *(0.d0,-2.d0)*evc (ig+npwx, ibnd)
        end do
     END IF
  enddo
!
! Uncomment this goto and the continue below to calculate 
! the matrix elements of p without the commutator with the
! nonlocal potential.
!
!  goto 111
  !
  ! and this is the contribution from nonlocal pseudopotentials
  !
  if (nkb == 0) go to 111
  !
  allocate (work ( npwx, nkb) )
  IF (noncolin) THEN
     allocate (deff_nc (nhm, nhm, nat, nspin))
  ELSE
     allocate (deff (nhm, nhm, nat ))
  END IF
  allocate (dvkb (npwx, nkb), dvkb1(npwx, nkb))
  dvkb (:,:) = (0.d0, 0.d0)
  dvkb1(:,:) = (0.d0, 0.d0)
 
  call gen_us_dj (ik, dvkb)
  call gen_us_dy (ik, at (1, ipol), dvkb1)
  do ig = 1, npw
     if (g2k (ig) < 1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     else
        gk (1, ig) = gk (1, ig) / sqrt (g2k (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2k (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2k (ig) )
     endif
  enddo

  jkb = 0
  work=(0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
              do ig = 1, npw
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, ipol) * gk (1, ig) + &
                       at (2, ipol) * gk (2, ig) + &
                       at (3, ipol) * gk (3, ig) )
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate (g2k, gk)

  ! In the case of gamma point systems becp2 is real
  ! so we have to include a factor of i before calling
  ! calbec otherwise we would be stuck with the wrong component
  ! of becp2 later on.
  IF (gamma_only) work=(0.0_DP,1.0_DP)*work
  CALL calbec (npw, work, evc, becp2)

  IF (noncolin) THEN
     allocate (psc ( nkb, npol, nbnd, 2))
     psc=(0.d0,0.d0)
  ELSE
     allocate (ps2 ( nkb, nbnd, 2))
     ps2=(0.d0,0.d0)
  END IF
  DO ibnd = 1, nbnd_occ
     IF (noncolin) THEN
        CALL compute_deff_nc(deff_nc,et(ibnd,ik))
     ELSE
        CALL compute_deff(deff,et(ibnd,ik))
     ENDIF
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (nt == ityp (na)) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    IF (noncolin) THEN
                       ijs=0
                       DO is=1, npol
                          DO js = 1, npol
                             ijs=ijs+1
                             psc(ikb,is,ibnd,1)=psc(ikb,is,ibnd,1)+  &
                                       (0.d0,-1.d0)*    &
                                  becp2%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs) 
                             psc(ikb,is,ibnd,2)=psc(ikb,is,ibnd,2)+ &
                                     (0.d0,-1.d0)* &
                                 becp1%nc(jkb,js,ibnd)*deff_nc(ih,jh,na,ijs) 
                          END DO
                       END DO
                    ELSEIF (gamma_only) THEN
                       ! Note the different prefactors due to the factor 
                       ! of i introduced to work(:,:), as becp[1,2] are
                       ! real.
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1) + becp2%r(jkb,ibnd) * &
                            (1.0d0, 0.0d0)*deff(ih,jh,na) 
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) + becp1%r(jkb,ibnd)* &
                            (-1.0d0, 0.0d0)*deff(ih,jh,na)
                    ELSE
                       ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1) + becp2%k(jkb,ibnd) * &
                            (0.0d0,-1.0d0)*deff(ih,jh,na) 
                       ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) + becp1%k(jkb,ibnd)* &
                            (0.0d0,-1.0d0)*deff(ih,jh,na)
                    END IF
                 enddo
              enddo
              ijkb0=ijkb0+nh(nt)
           end if
        enddo  ! na
     end do  ! nt
  end do ! nbnd
  if (ikb /= nkb .OR. jkb /= nkb) call errore ('commutator_Hx_psi', 'unexpected error',1)
  IF (noncolin) THEN
     CALL zgemm( 'N', 'N', npw, nbnd_occ*npol, nkb, &
          (1.d0,0.d0), vkb(1,1), npwx, psc(1,1,1,1), nkb, (1.d0,0.d0), &
          dpsi, npwx )
     CALL zgemm( 'N', 'N', npw, nbnd_occ*npol, nkb, &
          (1.d0,0.d0),work(1,1), npwx, psc(1,1,1,2), nkb, (1.d0,0.d0), &
          dpsi, npwx )
  ELSE
     CALL zgemm( 'N', 'N', npw, nbnd_occ, nkb, &
          (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
          dpsi(1,1), npwx )
     CALL zgemm( 'N', 'N', npw, nbnd_occ, nkb, &
          (1.d0,0.d0),work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
          dpsi(1,1), npwx )
  ENDIF

  IF (noncolin) THEN
     deallocate (psc)
     deallocate (deff_nc)
  ELSE
     deallocate (ps2)
     deallocate (deff)
  END IF
  deallocate (work)

  IF (nkb > 0) THEN
     deallocate (dvkb1, dvkb)
  END IF

  111 continue

  call stop_clock ('commutator_Hx_psi')
  return
end subroutine commutator_Hx_psi
