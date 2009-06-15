MODULE paw_postproc

    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info
    IMPLICIT NONE

    PUBLIC :: PAW_make_ae_charge
    PRIVATE

    CONTAINS

SUBROUTINE PAW_make_ae_charge(rho)
   USE paw_onecenter,     ONLY : paw_rho_lm
   USE atom,              ONLY : g => rgrid
   USE ions_base,         ONLY : nat, ityp, tau
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf
   USE scf,               ONLY : scf_type
   USE fft_base,          ONLY : dfftp
   USE mp_global,         ONLY : me_pool
   USE splinelib,         ONLY : spline, splint
   USE gvect,             ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
   USE cell_base,         ONLY : at, bg, alat

   TYPE(scf_type), INTENT(INOUT) :: rho
   TYPE(paw_info)          :: i                     ! minimal info on atoms
   INTEGER                 :: ipol                  ! counter on x,y,z
   INTEGER                 :: ir                    ! counter on grid point
   INTEGER                 :: is                    ! spin index
   INTEGER                 :: lm                    ! counters on angmom and radial grid
   INTEGER                 :: j,k,l, idx, idx0
   INTEGER                 :: ia
   REAL(DP),ALLOCATABLE    :: wsp_lm(:,:,:), ylm_posi(:,:), d1y(:), d2y(:)
   REAL(DP),ALLOCATABLE    :: rho_lm(:,:,:), rho_lm_ae(:,:,:), rho_lm_ps(:,:,:)
   REAL(DP)                :: posi(3), first, second
   REAL(DP)                :: inv_nr1, inv_nr2, inv_nr3, distsq
   INTEGER, EXTERNAL :: ldim_block, gind_block

   CALL start_clock('PAW_ae')
   ! Some initialization
   !
   inv_nr1 = 1.D0 / DBLE( nr1 )
   inv_nr2 = 1.D0 / DBLE( nr2 )
   inv_nr3 = 1.D0 / DBLE( nr3 )
   !
   ! I cannot parallelize on atoms, because it is already parallelized
   ! on charge slabs
   !
   atoms: DO ia = 1, nat
      !
      i%a = ia                      ! atom's index
      i%t = ityp(ia)                ! type of atom ia
      i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
      i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
      i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
      !
      ifpaw: IF (upf(i%t)%tpawp) THEN
         !
         ! Arrays are allocated inside the cycle to allow reduced
         ! memory usage as differnt atoms have different meshes
         ALLOCATE(rho_lm_ae(i%m,i%l**2,nspin), &
                  rho_lm_ps(i%m,i%l**2,nspin) )
         ALLOCATE(rho_lm(i%m,i%l**2,nspin), &
                  ylm_posi(1,i%l**2),       &
                  wsp_lm(i%m, i%l**2,nspin)  )
         !
         ! Compute rho spherical harmonics expansion from becsum and pfunc
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, upf(i%t)%qfuncl)
         !
         ! Add the components, with some soft step function
         DO is=1,nspin
         DO lm = 1,i%l*2
         DO ir = 1, i%m
!  soften the charge with a step function beyond kkbeta
! (stupid, because there AE-PS charge is exactly zero)
!                 rho_lm(ir,lm,is) = soft_step( g(i%t)%r(ir),  &
!                                              0.9_dp * g(i%t)%r(upf(i%t)%kkbeta), &
!                                              0.1_dp *  g(i%t)%r(upf(i%t)%kkbeta) ) &
!                                * ( rho_lm_ae(ir,lm,is)-rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir)

!  do not soften 
               rho_lm(ir,lm,is) = ( rho_lm_ae(ir,lm,is)-rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir)
!                write(1000*ia+lm, '(2f12.6)') g(i%t)%r(ir), rho_lm(ir,lm,is)
!                if( rho_lm(ir,lm,is) > 1._dp) rho_lm(ir,lm,is) = 1._dp
         ENDDO
         ENDDO
         ENDDO
         ! deallocate asap
         DEALLOCATE(rho_lm_ae, rho_lm_ps)
         !
         ! now let's rock!!!
         ALLOCATE( d1y(upf(i%t)%kkbeta), d2y(upf(i%t )%kkbeta) )
         DO is = 1,nspin
         DO lm = 1, i%l**2
             CALL radial_gradient(rho_lm(1:upf(i%t)%kkbeta,lm,is), d1y, &
                                  g(i%t)%r, upf(i%t)%kkbeta, 1)
             CALL radial_gradient(d1y, d2y, g(i%t)%r, upf(i%t)%kkbeta, 1)
             !
             first  = d1y(1) ! first derivative in first point
             second = d2y(1) ! second derivative in first point
             ! prepare interpolation
             CALL spline( g(i%t)%r(:), rho_lm(:,lm,is), first, second, wsp_lm(:,lm,is) )
         ENDDO
         ENDDO
         DEALLOCATE(d1y, d2y)
         !
#if defined (__PARA)
         idx0 = nrx1*nrx2 * SUM ( dfftp%npp(1:me_pool) )
#else
         idx0 = 0
#endif
         rsp_point : &
         DO ir = 1, nrxx
            !
            ! three dimensional indices (i,j,k)
            idx   = idx0 + ir - 1
            k     = idx / (nrx1*nrx2)
            idx   = idx - (nrx1*nrx2)*k
            j     = idx / nrx1
            idx   = idx - nrx1*j
            l     = idx
            !
            ! ... do not include points outside the physical range!
            IF ( l >= nr1 .or. j >= nr2 .or. k >= nr3 ) CYCLE rsp_point
            !
            DO ipol = 1, 3
               posi(ipol) = DBLE( l )*inv_nr1*at(ipol,1) + &
                            DBLE( j )*inv_nr2*at(ipol,2) + &
                            DBLE( k )*inv_nr3*at(ipol,3)
            END DO
            !
            ! find the distance of real-space grid's point ir w.r.t
            ! closer periodic image of atom ia
            posi(:) = posi(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi, bg, -1 )
            posi(:) = posi(:) - ANINT( posi(:) )
            CALL cryst_to_cart( 1, posi, at, 1 )
            !
            posi(:) = posi(:) * alat
            distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
            ! don't consider points too far from the atom:
            IF ( distsq > g(i%t)%r2(upf(i%t)%kkbeta) ) &
               CYCLE rsp_point
            !
            ! I have to generate the atomic charge on point posi(:), which means
            ! sum over l and m components rho_lm_ae-rho_lm_ps
            ! I also have to interpolate the radial function at distance |posi(:)|
            !
            ! prepare spherical harmonics
            CALL ylmr2( i%l**2, 1, posi, distsq, ylm_posi )
            DO is = 1,nspin
            DO lm = 1, i%l**2
                ! do interpolation
                rho%of_r(ir,is)= rho%of_r(ir,is) + ylm_posi(1,lm) &
                                * splint(g(i%t)%r(:) , rho_lm(:,lm,is), wsp_lm(:,lm,is), SQRT(distsq) )
            ENDDO
            ENDDO
         END DO rsp_point

         DEALLOCATE(rho_lm, ylm_posi, wsp_lm)
         !
      ENDIF ifpaw
   ENDDO atoms

    ! put energy back in the output variable
   CALL stop_clock('PAW_ae')

END SUBROUTINE PAW_make_ae_charge

#ifdef __COMPILE_THIS_UNUSED_FUNCTION
FUNCTION soft_step(x,x0,k)
    !
    REAL(DP) :: soft_step
    REAL(DP),INTENT(IN) :: x
    REAL(DP),OPTIONAL,INTENT(IN) :: x0,k
    REAL(DP) :: x0_,k_
    !
    k_ = 1._dp
    if(present(k) ) k_ = k
    !
    x0_ = 0._dp
    if(present(x0) ) x0_ = x0
    !
    soft_step = 1._dp - 1._dp / (1._dp+exp(-2._dp*(x-x0_)/k_))
    !write(*,'(4f12.6)') soft_step, x,x0_,k_
    !
END FUNCTION
#endif

END MODULE paw_postproc






