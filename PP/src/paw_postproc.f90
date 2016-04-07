MODULE paw_postproc

    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info
    IMPLICIT NONE

    PUBLIC :: PAW_make_ae_charge
    PRIVATE

    CONTAINS

! Reconstruct the valence (withcore=.false.) or all-electron
! (withcore=.true.) density using the PAW transformation and the
! information in the UPF file.
SUBROUTINE PAW_make_ae_charge(rho,withcore)
   USE constants,         ONLY : sqrtpi
   USE paw_onecenter,     ONLY : paw_rho_lm
   USE atom,              ONLY : g => rgrid
   USE ions_base,         ONLY : nat, ityp, tau
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf
   USE scf,               ONLY : scf_type
   USE fft_base,          ONLY : dfftp
   USE mp_global,         ONLY : me_pool
   USE splinelib,         ONLY : spline, splint
   USE cell_base,         ONLY : at, bg, alat

   TYPE(scf_type), INTENT(inout) :: rho
   LOGICAL, INTENT(IN) :: withcore
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

   ! Some initialization
   !
   inv_nr1 = 1.D0 / dble(  dfftp%nr1 )
   inv_nr2 = 1.D0 / dble(  dfftp%nr2 )
   inv_nr3 = 1.D0 / dble(  dfftp%nr3 )
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
         ! memory usage as different atoms have different meshes
         ALLOCATE(rho_lm_ae(i%m,i%l**2,nspin), &
                  rho_lm_ps(i%m,i%l**2,nspin) )
         ALLOCATE(rho_lm(i%m,i%l**2,nspin), &
                  ylm_posi(1,i%l**2),       &
                  wsp_lm(i%m, i%l**2,nspin)  )
         !
         ! Compute rho spherical harmonics expansion from becsum and pfunc
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, &
              upf(i%t)%qfuncl)
         !
         DO is=1,nspin
            DO lm = 1,i%l**2
               DO ir = 1, i%m
                  rho_lm(ir,lm,is) = ( rho_lm_ae(ir,lm,is) - &
                       rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir)
               ENDDO
            ENDDO
            !
            ! add core charge
            !
            IF (withcore) THEN
               DO ir = 1, i%m
                  rho_lm(ir,1,is) = rho_lm(ir,1,is) + &
                     upf(i%t)%paw%ae_rho_atc(ir) / nspin * (2._DP * sqrtpi)
               ENDDO
            ENDIF
         ENDDO

         ! deallocate asap
         DEALLOCATE(rho_lm_ae, rho_lm_ps)
         !
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
               CALL spline( g(i%t)%r(:), rho_lm(:,lm,is), first, second, &
                    wsp_lm(:,lm,is) )
            ENDDO
         ENDDO
         DEALLOCATE(d1y, d2y)
         !
         ! idx0 = starting index of real-space FFT arrays for this processor
         idx0 =  dfftp%nr1x* dfftp%nr2x * dfftp%ipp(me_pool+1)
         !
         rsp_point : DO ir = 1, dfftp%nr1x*dfftp%nr2x * dfftp%npl
            !
            ! three dimensional indices (i,j,k)
            idx   = idx0 + ir - 1
            k     = idx / ( dfftp%nr1x* dfftp%nr2x)
            idx   = idx - ( dfftp%nr1x* dfftp%nr2x)*k
            j     = idx /  dfftp%nr1x
            idx   = idx -  dfftp%nr1x*j
            l     = idx
            !
            ! ... do not include points outside the physical range!
            IF ( l >=  dfftp%nr1 .or. j >=  dfftp%nr2 .or. k >=  dfftp%nr3 ) CYCLE rsp_point
            !
            DO ipol = 1, 3
               posi(ipol) = dble( l )*inv_nr1*at(ipol,1) + &
                            dble( j )*inv_nr2*at(ipol,2) + &
                            dble( k )*inv_nr3*at(ipol,3)
            ENDDO
            !
            ! find the distance of real-space grid's point ir w.r.t
            ! closer periodic image of atom ia
            !
            posi(:) = posi(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi, bg, -1 )
            posi(:) = posi(:) - anint( posi(:) )
            CALL cryst_to_cart( 1, posi, at, 1 )
            !
            posi(:) = posi(:) * alat
            distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
            ! don't consider points too far from the atom:
            IF ( distsq > g(i%t)%r2(upf(i%t)%kkbeta) ) &
                 CYCLE rsp_point
            !
            ! generate the atomic charge on point posi(:), which means
            ! sum over l and m components rho_lm_ae-rho_lm_ps
            ! interpolate the radial function at distance |posi(:)|
            !
            ! prepare spherical harmonics
            CALL ylmr2( i%l**2, 1, posi, distsq, ylm_posi )
            DO is = 1,nspin
               DO lm = 1, i%l**2
                  ! do interpolation
                  rho%of_r(ir,is)= rho%of_r(ir,is) + ylm_posi(1,lm) &
                       * splint(g(i%t)%r(:) , rho_lm(:,lm,is), &
                       wsp_lm(:,lm,is), sqrt(distsq) )
               ENDDO
            ENDDO
         ENDDO rsp_point
         !
         DEALLOCATE(rho_lm, ylm_posi, wsp_lm)
         !
      ENDIF ifpaw
   ENDDO atoms

 END SUBROUTINE PAW_make_ae_charge

END MODULE paw_postproc

