! This subroutine calculates the commutator [r,V_nonloc]


subroutine gen_beta_simple (qk, npw_max, dvkb)
  !----------------------------------------------------------------------
  !
  !   Calculates the beta function pseudopotentials with
  !   the derivative of the Bessel functions
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE klist,      ONLY : ngk
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  USE uspp,       ONLY : nkb, indv, nhtol, nhtolm
  USE us,         ONLY : nqx, tab, tab_d2y, dq, spline_ps
  USE m_gth,      ONLY : mk_dffnl_gth
  USE splinelib
  USE uspp_param, ONLY : upf, lmaxkb, nbetam, nh
  USE io_global, ONLY : stdout
  !
  implicit none
  !
  real(kind=DP), dimension(3) :: qk
  integer :: npw_max
  complex(DP), intent(out) :: dvkb (npw_max, nkb)
  !
  ! local variables
  !
  integer :: ikb, nb, ih, ig, i0, i1, i2, i3 , nt
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on G vectors
  ! index of the first nonzero point in the r
  ! counter on atomic type

  real(DP) :: arg, px, ux, vx, wx
  ! argument of the atomic phase factor

  complex(DP) :: phase, pref
  ! atomic phase factor
  ! prefactor

  integer :: na, l, iig, lm, iq
  real(DP), allocatable :: djl (:,:,:), ylm (:,:), q (:), gk (:,:)
  real(DP) ::  qt
  complex(DP), allocatable :: sk (:)
  real(DP), allocatable :: xdata(:)

  call start_clock('gen_beta1')

  if (nkb.eq.0) return

  call start_clock('stres_us31')

  allocate (djl( npw_max , nbetam , ntyp))
  allocate (ylm( npw_max ,(lmaxkb + 1) **2))
  allocate (gk( 3, npw_max))
  allocate (q( npw_max))
  do ig = 1, npw_max
     iig = ig
     gk (1,ig) = qk(1) + g(1, iig)
     gk (2,ig) = qk(2) + g(2, iig)
     gk (3,ig) = qk(3) + g(3, iig)
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  call stop_clock('stres_us31')
  call start_clock('stres_us32')
  call ylmr2 ((lmaxkb+1)**2, npw_max, gk, q, ylm)
  call stop_clock('stres_us32')
  call start_clock('stres_us33')

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  do nt = 1, ntyp
     do nb = 1, upf(nt)%nbeta
        if ( upf(nt)%is_gth ) then
           call mk_dffnl_gth( nt, nb, npw_max, q, djl(1,nb,nt) )
           cycle
        endif
        do ig = 1, npw_max
           qt = sqrt(q (ig)) * tpiba
           if (spline_ps) then
             djl(ig,nb,nt) = splint_deriv(xdata, tab(:,nb,nt), &
                                                 tab_d2y(:,nb,nt), qt)
           else
             px = qt / dq - int (qt / dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = qt / dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             if (i3 <= nqx) then ! Approximation
                djl(ig,nb,nt) = ( tab (i0, nb, nt) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                               tab (i1, nb, nt) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                               tab (i2, nb, nt) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                               tab (i3, nb, nt) * (+ux*vx-px*vx-px*ux)/6.d0 )/dq
             else
                djl(ig,nb,nt) = 0.d0  ! Approximation
             endif
           endif
        enddo
     enddo
  enddo
  call stop_clock('stres_us33')
  call start_clock('stres_us34')

  deallocate (q)
  deallocate (gk)

  allocate (sk( npw_max))
  ikb = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           arg = (qk (1) * tau(1,na) + &
                  qk (2) * tau(2,na) + &
                  qk (3) * tau(3,na) ) * tpi
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           do ig = 1, npw_max
              iig = ig
              sk (ig) = eigts1 (mill (1,iig), na) * &
                        eigts2 (mill (2,iig), na) * &
                        eigts3 (mill (3,iig), na) * phase
           enddo
           do ih = 1, nh (nt)
              nb = indv (ih, nt)
              l = nhtol (ih, nt)
              lm= nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0) **l
              !
              do ig = 1, npw_max
                 dvkb (ig, ikb) = djl (ig, nb, nt) * sk (ig) * ylm (ig, lm) &
                      * pref
              enddo
           enddo
        endif
     enddo

  enddo
  call stop_clock('stres_us34')
  call stop_clock('gen_beta1')

  if (ikb.ne.nkb) call errore ('gen_us_dj', 'unexpected error', 1)
  deallocate (sk)
  deallocate (ylm)
  deallocate (djl)
  if (spline_ps) deallocate(xdata)
  return
end subroutine gen_beta_simple



subroutine gen_beta_simple_2 (qk, npw_max, u, dvkb)
  !----------------------------------------------------------------------
  !
  !  Calculates the kleinman-bylander pseudopotentials with the
  !  derivative of the spherical harmonics projected on vector u
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : tpi
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE klist,      ONLY : ngk, igk_k
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  USE uspp,       ONLY : nkb, indv, nhtol, nhtolm
  USE us,         ONLY : nqx, tab, tab_d2y, dq, spline_ps
  USE splinelib
  USE uspp_param, ONLY : upf, lmaxkb, nbetam, nh
  !
  implicit none
  !
  real(kind=DP), dimension(3) :: qk
  integer :: npw_max
  real(DP) :: u (3)

  complex(DP) :: dvkb (npw_max, nkb)
  integer :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
       i3, ig
  real(DP), allocatable :: gk(:,:), q (:)
  real(DP) :: px, ux, vx, wx, arg

  real(DP), allocatable :: vkb0 (:,:,:), dylm (:,:), dylm_u (:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u

  complex(DP), allocatable :: sk (:)
  complex(DP) :: phase, pref

  integer :: iq
  real(DP), allocatable :: xdata(:)
  !
  !
  call start_clock('gen_beta2')

  dvkb(:,:) = (0.d0, 0.d0)
  if (lmaxkb.le.0) return

  allocate ( vkb0(npw_max,nbetam,ntyp), dylm_u(npw_max,(lmaxkb+1)**2), gk(3,npw_max) )
  allocate ( q(npw_max) )

  do ig = 1, npw_max
     iig = ig
     gk (1, ig) = qk (1) + g (1,iig)
     gk (2, ig) = qk (2) + g (2,iig)
     gk (3, ig) = qk (3) + g (3,iig)
     q (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo

  allocate ( dylm(npw_max,(lmaxkb+1)**2) )
  dylm_u(:,:) = 0.d0
  do ipol = 1, 3
     call dylmr2  ((lmaxkb+1)**2, npw_max, gk, q, dylm, ipol)
     call daxpy (npw_max * (lmaxkb + 1) **2, u (ipol), dylm, 1, dylm_u, 1)
  enddo
  deallocate (dylm)

  do ig = 1, npw_max
     q (ig) = sqrt ( q(ig) ) * tpiba
  end do

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

  do nt = 1, ntyp
     ! calculate beta in G-space using an interpolation table
     do nb = 1, upf(nt)%nbeta
        do ig = 1, npw_max
           if (spline_ps) then
             vkb0(ig,nb,nt) = splint(xdata, tab(:,nb,nt), &
                                     tab_d2y(:,nb,nt), q(ig))
           else
             px = q (ig) / dq - int (q (ig) / dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = q (ig) / dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             if (i3<=nqx) then ! DEBUG
                vkb0 (ig, nb, nt) = tab (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                                    tab (i1, nb, nt) * px * vx * wx / 2.d0 - &
                                   tab (i2, nb, nt) * px * ux * wx / 2.d0 + &
                                    tab (i3, nb, nt) * px * ux * vx / 6.d0
             else
                vkb0 (ig, nb, nt) = 0.d0 ! DEBUG
             endif
           endif
        enddo
     enddo
  enddo

  deallocate (q)
  allocate ( sk(npw_max) )

  ikb = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           arg = (qk (1) * tau (1, na) + qk (2) * tau (2, na) &
                + qk (3) * tau (3, na) ) * tpi
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           do ig = 1, npw_max
              iig = ig
              sk (ig) = eigts1 (mill (1,iig), na) * &
                        eigts2 (mill (2,iig), na) * &
                        eigts3 (mill (3,iig), na) * phase
           enddo
           do ih = 1, nh (nt)
              nb = indv (ih, nt)
              l = nhtol (ih, nt)
              lm = nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0) **l
              !
              do ig = 1, npw_max
                 dvkb (ig, ikb) = vkb0(ig, nb, nt) * sk(ig) * dylm_u(ig, lm) &
                      * pref / tpiba
              enddo
           enddo
        endif
     enddo
  enddo

  call stop_clock('gen_beta2')

  if (ikb.ne.nkb) then
     WRITE( stdout, * ) ikb, nkb
     call errore ('gen_us_dy', 'unexpected error', 1)
  endif

  deallocate ( sk )
  deallocate ( vkb0, dylm_u, gk )
  if (spline_ps) deallocate(xdata)

  return
end subroutine gen_beta_simple_2


subroutine commutator_Hx_psi_simple (qk, npw_max, nbnd_occ, becp1, becp2, ipol, dpsi, wfc_e, vkb_max)
  !----------------------------------------------------------------------
  !
  ! On output: dpsi contains [H,x_ipol] | psi_ik > in crystal axis
  !            (projected on at(*,ipol) )
  !
  ! vkb,evc,igk must be properly set for the appropriate k-point
  ! in addition becp1 must be set equal to becp1 = <vkb|evc>
  ! as it is done in PH/phq_init.f90 for the k-point ik
  ! NB: here the last index of becp1 is missing, hence it refers
  !     to a single k-point
  !
  !    CALL calbec (npw, vkb, evc, becp1(:,:) )
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : tpiba
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp
  USE io_global,       ONLY : stdout
  USE gvect,           ONLY : g
  USE wvfct,           ONLY : nbnd, et
  USE wavefunctions, ONLY: evc
  USE lsda_mod,        ONLY : nspin
  USE noncollin_module,ONLY : noncolin, npol
  USE becmod,          ONLY : becp, bec_type, calbec
  USE uspp,            ONLY : nkb
  USE uspp_param,      ONLY : nh, nhm
  USE control_flags,   ONLY : gamma_only
  USE uspp,          ONLY:  deeq, deeq_nc
  USE wvfct,           ONLY : npwx
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp, ONLY : mp_sum , mp_bcast

  implicit none
  integer :: npw_max
  INTEGER, INTENT(IN) :: nbnd_occ, ipol
  COMPLEX(DP), INTENT(OUT)    :: dpsi(npw_max*npol,nbnd_occ)
  TYPE(bec_type), INTENT(IN)  :: becp1 ! dimensions ( nkb, nbnd )
  TYPE(bec_type), INTENT(INOUT) :: becp2 ! dimensions ( nkb, nbnd )
  COMPLEX(kind=DP), INTENT(IN)  :: vkb_max(npw_max, nbnd_occ) 
  !
  COMPLEX(kind=DP), DIMENSION(npol*npw_max,nbnd_occ) :: wfc_e
  REAL(kind=DP), DIMENSION(3) :: qk
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: g2kin
  REAL(kind=DP), dimension(3,3) :: att
  INTEGER :: nbegin, nend , nbnd_loc
  !
  ! Local variables
  !
  integer :: ig, na, ibnd, jbnd, ikb, jkb, nt, lter, ih, jh, ijkb0,  &
             nrec, is, js, ijs , ii
  ! counters

  real(DP), allocatable  :: gk (:,:)
  ! the derivative of |k+G|
  complex(DP), allocatable :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:),  &
       work (:,:), psc(:,:,:,:), deff_nc(:,:,:,:)
  REAL(DP), allocatable :: deff(:,:,:)
  !

  CALL start_clock ('commut_Hx_psi')
  dpsi=(0.d0, 0.d0)
  !
  allocate(g2kin(npw_max))
  allocate (gk ( 3, npw_max))
  do ig = 1, npw_max
     gk (1:3, ig) = (qk (1:3) + g (1:3, ig ) ) * tpiba
     g2kin (ig) = SUM(gk (1:3, ig) **2 )
  enddo
  !
  !
  ! the contribution from nonlocal pseudopotentials
  !
  if (nkb == 0) go to 111
  !
  allocate (work ( npw_max, nkb) )
  IF (noncolin) THEN
     allocate (deff_nc (nhm, nhm, nat, nspin))
  ELSE
     allocate (deff (nhm, nhm, nat ))
  END IF
  allocate (dvkb (npw_max, nkb), dvkb1(npw_max, nkb))
  dvkb (:,:) = (0.d0, 0.d0)
  dvkb1(:,:) = (0.d0, 0.d0)

  att = 0.0d0
  do ii = 1,3
    att(ii,ii) = 1.0d0
  enddo

  call gen_beta_simple (qk, npw_max, dvkb)
  call gen_beta_simple_2 (qk, npw_max, att(1,ipol), dvkb1)

  do ig = 1, npw_max
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
              do ig = 1, npw_max
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (att (1, ipol) * gk (1, ig) + &
                       att (2, ipol) * gk (2, ig) + &
                       att (3, ipol) * gk (3, ig) )
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate (gk)
  

  ! In the case of gamma point systems becp2 is real
  ! so we have to include a factor of i before calling
  ! calbec otherwise we would be stuck with the wrong component
  ! of becp2 later on.
  IF (gamma_only) work=(0.0_DP,1.0_DP)*work
  work(npwx+1:npw_max,1:nkb) = 0.d0          ! DEBUG

  CALL calbec (npw_max, work, wfc_e, becp2)

  IF (noncolin) THEN
     allocate (psc ( nkb, npol, nbnd_occ, 2))
     psc=(0.d0,0.d0)
  ELSE
     allocate (ps2 ( nkb, nbnd_occ, 2))
     ps2=(0.d0,0.d0)
  END IF

  call start_clock('commut_nbnd')

!  if(nbnd_loc*nproc < (nbnd_occ)) nbnd_loc = nbnd_loc + 1  ! DEBUG
!  nbegin = mpime*nbnd_loc + 1
!  nend = nbegin+nbnd_loc - 1
!  if(nend > nbnd_occ) nend=nbnd_occ  ! DEBUG

  ! Parallelization over nbnd_occ
  nbnd_loc= (nbnd_occ)/nproc  ! nbnd_loc = num of bands per processor
  nbegin = mpime*nbnd_loc + 1
  nend = nbegin+nbnd_loc - 1
  if(mpime == (nproc-1)) nend=nbnd_occ  ! last processor takes all the remaining bands
  !
  DO ibnd = nbegin , nend
     IF (noncolin) THEN
        deff_nc=deeq_nc
     ELSE
        deff(:,:,:) = deeq(:,:,:,1)
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
        !write(*,*) nt , na , ikb , jkb ! DEBUG
     end do  ! nt
  end do ! nbnd
  if (ikb /= nkb .OR. jkb /= nkb) call errore ('commutator_Hx_psi_simple', 'unexpected error',1)
  !
  IF (noncolin) THEN
     call mp_sum(psc,world_comm)
  ELSE
     call mp_sum(ps2,world_comm)
  ENDIF
  !
  call stop_clock('commut_nbnd')
  !
  call start_clock('commut_zgemm')
  IF (noncolin) THEN
     CALL zgemm( 'N', 'N', npw_max, nbnd_occ*npol, nkb, &
          (1.d0,0.d0), vkb_max(1,1), npw_max, psc(1,1,1,1), nkb, (1.d0,0.d0), &
          dpsi, npw_max )
     CALL zgemm( 'N', 'N', npw_max, nbnd_occ*npol, nkb, &
          (1.d0,0.d0),work(1,1), npw_max, psc(1,1,1,2), nkb, (1.d0,0.d0), &
          dpsi, npw_max )
  ELSE
     CALL zgemm( 'N', 'N', npw_max, nbnd_occ, nkb, &
          (1.d0,0.d0), vkb_max(1,1), npw_max, ps2(1,1,1), nkb, (1.d0,0.d0), &
          dpsi(1,1), npw_max )
     CALL zgemm( 'N', 'N', npw_max, nbnd_occ, nkb, &
          (1.d0,0.d0),work(1,1), npw_max, ps2(1,1,2), nkb, (1.d0,0.d0), &
          dpsi(1,1), npw_max )
  ENDIF
  call stop_clock('commut_zgemm')

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

  deallocate(g2kin)
  call stop_clock ('commut_Hx_psi')
  return
end subroutine commutator_Hx_psi_simple
