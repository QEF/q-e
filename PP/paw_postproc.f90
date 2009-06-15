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
   USE cell_base,         ONLY : at, bg

   TYPE(scf_type), INTENT(INOUT) :: rho
   TYPE(paw_info)          :: i                     ! minimal info on atoms
   INTEGER                 :: ipol                  ! counter on x,y,z
   INTEGER                 :: ir                    ! counter on grid point
   INTEGER                 :: is                    ! spin index
   INTEGER                 :: lm                    ! counters on angmom and radial grid
   INTEGER                 :: j,k,l, idx, idx0
   INTEGER                 :: na_loc, ia_e, ia_s, ia
   REAL(DP),ALLOCATABLE    :: wsp(:), ylm_posi(:,:)
   REAL(DP),ALLOCATABLE    :: rho_lm(:,:,:), rho_lm_ae(:,:,:), rho_lm_ps(:,:,:)
   REAL(DP)                :: posi(3)
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
                  rho_lm_ps(i%m,i%l**2,nspin), &
                  rho_lm(i%m,i%l**2,nspin),    &
                  ylm_posi(1,i%l**2),       &
                  wsp(i%m)  )
         !
         ! Compute rho spherical harmonics expansion from becsum and pfunc
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, upf(i%t)%qfuncl)
         !
         DO is=1,nspin
         DO lm = 1,i%l*2
            rho_lm(:,lm,is) = ( rho_lm_ae(:,lm,is)-rho_lm_ps(:,lm,is) ) * g(i%t)%rm2(:)
         ENDDO
         ENDDO
         ! deallocate asap
         DEALLOCATE(rho_lm_ae, rho_lm_ps)
         !
         ! now let's rock!!!
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
                ! prepare interpolation
                CALL spline( g(i%t)%r(:), rho_lm(:,lm,is), 0._dp, 0._dp, wsp )
                ! do interpolation
                rho%of_r(ir,is)= rho%of_r(ir,is) &
                               + ylm_posi(1,lm) &
                                * splint(g(i%t)%r(:) , rho_lm(:,lm,is), wsp, SQRT(distsq) )
            ENDDO
            ENDDO
         END DO rsp_point

         DEALLOCATE(rho_lm, ylm_posi, wsp)
         !
      ENDIF ifpaw
   ENDDO atoms

    ! put energy back in the output variable
   CALL stop_clock('PAW_ae')

END SUBROUTINE PAW_make_ae_charge


END MODULE paw_postproc