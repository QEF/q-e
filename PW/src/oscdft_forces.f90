MODULE oscdft_forces
#if defined (__OSCDFT)
   USE kinds,               ONLY : DP
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type

   PRIVATE
   PUBLIC oscdft_forces_type,&
          oscdft_derivatives_gamma_type,&
          oscdft_derivatives_k_type,&
          oscdft_derivatives_gamma_alloc,&
          oscdft_derivatives_gamma_dealloc,&
          oscdft_derivatives_k_alloc,&
          oscdft_derivatives_k_dealloc

   TYPE oscdft_derivatives_gamma_type
      ! nwfcO is oscdft_ctx%forces%wfcO%n
      REAL(DP), ALLOCATABLE :: dwfatpsi(:,:),& ! size(nwfcO,nbnd)
                               wfatdpsi(:,:),& ! size(nwfcO,nbnd)
                               wfatbeta(:,:),& ! size(nwfcO,nkb)
                               wfatdbeta(:,:),& ! size(nwfcO,nkb)
                               betapsi(:,:),& ! size(nkb,nbnd)
                               dbetapsi(:,:) ! size(nkb,nbnd)
   END TYPE oscdft_derivatives_gamma_type

   TYPE oscdft_derivatives_k_type
      ! nwfcO is oscdft_ctx%forces%wfcO%n
      COMPLEX(DP), ALLOCATABLE :: dwfatpsi(:,:),& ! size(nwfcO,nbnd)
                                  wfatdpsi(:,:),& ! size(nwfcO,nbnd)
                                  wfatbeta(:,:),& ! size(nwfcO,nkb)
                                  wfatdbeta(:,:),& ! size(nwfcO,nkb)
                                  betapsi(:,:),& ! size(nkb,nbnd)
                                  dbetapsi(:,:) ! size(nkb,nbnd)
   END TYPE oscdft_derivatives_k_type

   TYPE oscdft_forces_type
      TYPE(oscdft_wavefunction_type)      :: wfcO
      TYPE(oscdft_derivatives_gamma_type) :: deriv_gamma
      TYPE(oscdft_derivatives_k_type)     :: deriv_k
   END TYPE oscdft_forces_type

   CONTAINS
      SUBROUTINE oscdft_derivatives_gamma_alloc(derivs, nwfcO, nbnd, nkb)
         IMPLICIT NONE

         TYPE(oscdft_derivatives_gamma_type), INTENT(INOUT) :: derivs
         INTEGER,                             INTENT(IN)    :: nwfcO, nbnd, nkb


         ALLOCATE(derivs%dwfatpsi(nwfcO,nbnd),&
                  derivs%wfatdpsi(nwfcO,nbnd),&
                  derivs%wfatbeta(nwfcO,nkb),&
                  derivs%wfatdbeta(nwfcO,nkb),&
                  derivs%betapsi(nkb,nbnd),&
                  derivs%dbetapsi(nkb,nbnd))

         !$acc enter data create(derivs%dwfatpsi,&
         !$acc&                  derivs%wfatdpsi,&
         !$acc&                  derivs%wfatbeta,&
         !$acc&                  derivs%wfatdbeta,&
         !$acc&                  derivs%betapsi,&
         !$acc&                  derivs%dbetapsi)

         ! not needed
#if 0
         !$acc kernels
         derivs%dwfatpsi  = 0.D0
         derivs%wfatdpsi  = 0.D0
         derivs%wfatbeta  = 0.D0
         derivs%wfatdbeta = 0.D0
         derivs%betapsi   = 0.D0
         derivs%dbetapsi  = 0.D0
         !$acc end kernels
#endif
      END SUBROUTINE oscdft_derivatives_gamma_alloc

      SUBROUTINE oscdft_derivatives_gamma_dealloc(derivs)
         IMPLICIT NONE

         TYPE(oscdft_derivatives_gamma_type), INTENT(INOUT) :: derivs

         !$acc exit data delete(derivs%dwfatpsi,&
         !$acc&                 derivs%wfatdpsi,&
         !$acc&                 derivs%wfatbeta,&
         !$acc&                 derivs%wfatdbeta,&
         !$acc&                 derivs%betapsi,&
         !$acc&                 derivs%dbetapsi)
         DEALLOCATE(derivs%dwfatpsi,&
                    derivs%wfatdpsi,&
                    derivs%wfatbeta,&
                    derivs%wfatdbeta,&
                    derivs%betapsi,&
                    derivs%dbetapsi)
      END SUBROUTINE oscdft_derivatives_gamma_dealloc

      SUBROUTINE oscdft_derivatives_k_alloc(derivs, nwfcO, nbnd, nkb)
         IMPLICIT NONE

         TYPE(oscdft_derivatives_k_type), INTENT(INOUT) :: derivs
         INTEGER,                         INTENT(IN)    :: nwfcO, nbnd, nkb

         ALLOCATE(derivs%dwfatpsi(nwfcO,nbnd),&
                  derivs%wfatdpsi(nwfcO,nbnd),&
                  derivs%wfatbeta(nwfcO,nkb),&
                  derivs%wfatdbeta(nwfcO,nkb),&
                  derivs%betapsi(nkb,nbnd),&
                  derivs%dbetapsi(nkb,nbnd))

         !$acc enter data create(derivs%dwfatpsi,&
         !$acc&                  derivs%wfatdpsi,&
         !$acc&                  derivs%wfatbeta,&
         !$acc&                  derivs%wfatdbeta,&
         !$acc&                  derivs%betapsi,&
         !$acc&                  derivs%dbetapsi)

         ! not needed
#if 0
         !$acc kernels
         derivs%dwfatpsi  = (0.D0, 0.D0)
         derivs%wfatdpsi  = (0.D0, 0.D0)
         derivs%wfatbeta  = (0.D0, 0.D0)
         derivs%wfatdbeta = (0.D0, 0.D0)
         derivs%betapsi   = (0.D0, 0.D0)
         derivs%dbetapsi  = (0.D0, 0.D0)
         !$acc end kernels
#endif
      END SUBROUTINE oscdft_derivatives_k_alloc

      SUBROUTINE oscdft_derivatives_k_dealloc(derivs)
         IMPLICIT NONE

         TYPE(oscdft_derivatives_k_type), INTENT(INOUT) :: derivs

         !$acc exit data delete(derivs%dwfatpsi,&
         !$acc&                 derivs%wfatdpsi,&
         !$acc&                 derivs%wfatbeta,&
         !$acc&                 derivs%wfatdbeta,&
         !$acc&                 derivs%betapsi,&
         !$acc&                 derivs%dbetapsi)
         DEALLOCATE(derivs%dwfatpsi,&
                    derivs%wfatdpsi,&
                    derivs%wfatbeta,&
                    derivs%wfatdbeta,&
                    derivs%betapsi,&
                    derivs%dbetapsi)
      END SUBROUTINE oscdft_derivatives_k_dealloc
#endif
END MODULE oscdft_forces
