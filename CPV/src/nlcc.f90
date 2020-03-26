!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   subroutine core_charge_ftr( tpre )
!=----------------------------------------------------------------------------=!
     !
     !  Compute the fourier transform of the core charge, from the radial
     !  mesh to the reciprocal space
     !
     use kinds,              ONLY : DP
     use ions_base,          ONLY : nsp
     use atom,               ONLY : rgrid
     use uspp,               ONLY : nlcc_any
     use uspp_param,         ONLY : upf
     use smallbox_gvec,      ONLY : ngb, gb
     use small_box,          ONLY : omegab, tpibab
     use pseudo_base,        ONLY : compute_rhocg
     use pseudopotential,    ONLY : tpstab, rhoc1_sp, rhocp_sp
     use cell_base,          ONLY : omega, tpiba2, tpiba
     USE splines,            ONLY : spline
     use gvect,              ONLY : gg, gstart
     USE core,               ONLY : rhocb, rhocg, drhocg
     USE fft_base,           ONLY: dfftp
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: tpre
     !
     INTEGER :: is, ig
     REAL(DP) :: xg, cost1
     !
     !
     IF( .NOT. nlcc_any ) RETURN
     !
     IF( .NOT. ALLOCATED( rgrid ) ) &
        CALL errore( ' core_charge_ftr ', ' rgrid not allocated ', 1 )
     IF( .NOT. ALLOCATED( upf ) ) &
        CALL errore( ' core_charge_ftr ', ' upf not allocated ', 1 )
     !
     do is = 1, nsp
        !
        if( upf(is)%nlcc ) then
           !
              CALL compute_rhocg( rhocb(:,is), rhocb(:,is), rgrid(is)%r, &
                  rgrid(is)%rab, upf(is)%rho_atc(:), gb, omegab, tpibab**2, &
                  rgrid(is)%mesh, ngb, 0 )
              !
           IF( tpre ) THEN
              !
              IF( tpstab ) THEN
                 !
                 cost1 = 1.0d0/omega
                 !
                 IF( gstart == 2 ) THEN
                    rhocg (1,is) = rhoc1_sp(is)%y( 1 ) * cost1
                    drhocg(1,is) = 0.0d0
                 END IF
                 DO ig = gstart, SIZE( rhocg, 1 )
                    xg = SQRT( gg(ig) ) * tpiba
                    rhocg (ig,is) = spline( rhoc1_sp(is), xg ) * cost1
                    drhocg(ig,is) = spline( rhocp_sp(is), xg ) * cost1
                 END DO
                 !
              ELSE

                 CALL compute_rhocg( rhocg(:,is), drhocg(:,is), rgrid(is)%r, &
                                     rgrid(is)%rab, upf(is)%rho_atc(:), gg, &
                                     omega, tpiba2, rgrid(is)%mesh, dfftp%ngm, 1 )

              END IF
              !
           END IF
           !
        endif
        !
     end do

     return
   end subroutine core_charge_ftr



!-----------------------------------------------------------------------
      subroutine add_cc( rhoc, rhog, rhor )
!-----------------------------------------------------------------------
!
! add core correction to the charge density for exch-corr calculation
!
      USE kinds,              ONLY: DP
      use electrons_base,     only: nspin
      use control_flags,      only: iverbosity
      use io_global,          only: stdout
      use mp_global,          only: intra_bgrp_comm
      use cell_base,          only: omega
      USE mp,                 ONLY: mp_sum

      ! this isn't really needed, but if I remove it, ifc 7.1
      ! gives an "internal compiler error"
      use gvect, only: gstart
      USE fft_interfaces,     ONLY: fwfft
      USE fft_base,           ONLY: dfftp
      USE fft_helper_subroutines, ONLY: fftx_add_threed2oned_gamma
!
      implicit none
      !
      REAL(DP),    INTENT(IN)   :: rhoc( dfftp%nnr )
      REAL(DP),    INTENT(INOUT):: rhor( dfftp%nnr, nspin )
      COMPLEX(DP), INTENT(INOUT):: rhog( dfftp%ngm,  nspin )
      !
      COMPLEX(DP), ALLOCATABLE :: wrk1( : )
!
      integer :: ig, ir, iss, isup, isdw
      REAL(DP) :: rsum
      !
      IF( iverbosity > 1 ) THEN
         rsum = SUM( rhoc ) * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         CALL mp_sum( rsum, intra_bgrp_comm )
         WRITE( stdout, 10 ) rsum 
10       FORMAT( 3X, 'Core Charge = ', D14.6 )
      END IF
      !
      ! In r-space:
      !
      if ( nspin .eq. 1 ) then
         iss=1
         call daxpy(dfftp%nnr,1.d0,rhoc,1,rhor(1,iss),1)
      else
         isup=1
         isdw=2
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isup),1)
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isdw),1)
      end if 
      !
      ! rhoc(r) -> rhoc(g)  (wrk1 is used as work space)
      !
      allocate( wrk1( dfftp%nnr ) )

      wrk1(:) = rhoc(:)

      call fwfft('Rho',wrk1, dfftp )
      !
      ! In g-space:
      !
      if (nspin.eq.1) then
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,iss) )
      else
         wrk1 = wrk1 * 0.5d0
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,isup) )
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,isdw) )
      end if

      deallocate( wrk1 )
!
      return
      end subroutine add_cc


!
!-----------------------------------------------------------------------
      subroutine force_cc(irb,eigrb,vxc,fion1)
!-----------------------------------------------------------------------
!
!     core correction force: f = \int V_xc(r) (d rhoc(r)/d R_i) dr
!     same logic as in newd - uses box grid. For parallel execution:
!     the sum over node contributions is done in the calling routine
!
      USE kinds,             ONLY: DP
      use electrons_base,    only: nspin
      use smallbox_gvec,     only: gxb, ngb
      use smallbox_subs,     only: fft_oned2box, boxdotgrid
      use cell_base,         only: omega
      use ions_base,         only: nsp, na, nat, ityp
      use small_box,         only: tpibab
      use uspp_param,        only: upf
      use core,              only: rhocb
      use fft_interfaces,    only: invfft
      use fft_base,          only: dfftb, dfftp
      use gvect,             only: gstart

      implicit none

! input
      integer, intent(in)        :: irb(3,nat)
      complex(dp), intent(in):: eigrb(ngb,nat)
      real(dp), intent(in)   :: vxc(dfftp%nnr,nspin)
! output
      real(dp), intent(inout):: fion1(3,nat)
! local
      integer :: iss, ix, ig, is, ia, nfft, isa
      real(dp) :: fac, res
      complex(dp) ci, facg
      complex(dp), allocatable :: qv(:), fg1(:), fg2(:)
      real(dp), allocatable :: fcc(:,:)

#if defined(_OPENMP)
      INTEGER :: itid, mytid, ntids, omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
!
      call start_clock( 'forcecc' )
      ci = (0.d0,1.d0)

      fac = omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3*nspin)

!$omp parallel default(none) &      
!$omp          shared(nsp, na, ngb, eigrb, dfftb, irb, ci, rhocb, &
!$omp                 gxb, nat, fac, upf, vxc, nspin, tpibab, fion1, ityp ) &
!$omp          private(mytid, ntids, is, ia, nfft, ig, qv, fg1, fg2, itid, res, ix, fcc, facg, iss )


      allocate( fcc( 3, nat ) )
      allocate( qv( dfftb%nnr ) )
      allocate( fg1( ngb ) )
      allocate( fg2( ngb ) )

      fcc(:,:) = 0.d0

#if defined(_OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
      itid  = 0
#endif

      do ia = 1, nat

         is = ityp(ia)

         if( .not. upf(is)%nlcc ) then
            cycle
         end if 

         nfft = 1

#if defined(__MPI)
         if ( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) ) &
            CYCLE
#endif

#if defined(_OPENMP)
         IF ( mytid /= itid ) THEN
            itid = MOD( itid + 1, ntids )
            CYCLE
         ELSE
            itid = MOD( itid + 1, ntids )
         END IF
#endif

         do ix=1,3
            if (nfft.eq.2) then
               do ig=1,ngb
                  facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                  fg1(ig) = eigrb(ig,ia  )*facg
                  fg2(ig) = eigrb(ig,ia+1)*facg
               end do
               CALL fft_oned2box( qv, fg1, fg2 )
            else
               do ig=1,ngb
                  facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                  fg1(ig) = eigrb(ig,ia)*facg
               end do
               CALL fft_oned2box( qv, fg1 )
            end if
!
            call invfft( qv, dfftb, ia )
            !
            ! note that a factor 1/2 is hidden in fac if nspin=2
            !
            do iss=1,nspin
               res = boxdotgrid(irb(:,ia),1,qv,vxc(:,iss))
               fcc(ix,ia) = fcc(ix,ia) + fac * res
               if (nfft.eq.2) then
                  res = boxdotgrid(irb(:,ia+1),2,qv,vxc(:,iss))
                  fcc(ix,ia+1) = fcc(ix,ia+1) + fac * res 
               end if
            end do
         end do
      end do

!
!$omp critical
      do ia = 1, nat
        fion1(:,ia) = fion1(:,ia) + fcc(:,ia)
      end do
!$omp end critical

      deallocate( qv )
      deallocate( fg1 )
      deallocate( fg2 )
      deallocate( fcc )

!$omp end parallel
!
      call stop_clock( 'forcecc' )

      return
      end subroutine force_cc


!
!-----------------------------------------------------------------------
      subroutine set_cc( rhoc )
!-----------------------------------------------------------------------
!
!     Calculate core charge contribution in real space, rhoc(r)
!     Same logic as for rhov: use box grid for core charges
! 
      use kinds, only: dp
      use ions_base,         only: nsp, na, nat, ityp
      use uspp_param,        only: upf
      use smallbox_gvec,     only: ngb
      use smallbox_subs,     only: fft_oned2box, box2grid
      use control_flags,     only: iprint
      use core,              only: rhocb
      use fft_interfaces,    only: invfft
      use fft_base,          only: dfftb, dfftp
      USE cp_main_variables, ONLY: irb, eigrb
      USE mp_global,         ONLY: nproc_bgrp, me_bgrp, inter_bgrp_comm, my_bgrp_id, nbgrp
      USE mp,                ONLY: mp_sum

      implicit none
! output
      real(dp), intent(out)  :: rhoc(dfftp%nnr)
! local
      integer :: ig, is, ia, isa, iia
      INTEGER :: nabox, iabox( nat )
      complex(dp), PARAMETER :: ci = (0.d0,1.d0)
      complex(dp), allocatable :: wrk1(:)
      complex(dp), allocatable :: qv(:), fg1(:), fg2(:)

      INTEGER :: mytid, ntids
#if defined(_OPENMP)
      INTEGER :: omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
!
      call start_clock( 'set_cc' )

      allocate( wrk1 ( dfftp%nnr ) )

      nabox = 0
      DO ia = 1, nat
         IF( .NOT. upf(ityp(ia))%nlcc ) CYCLE
         IF( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2 ( ia ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) CYCLE 
         nabox = nabox + 1
         iabox( nabox ) = ia
      END DO
!
!$omp parallel default(none) &      
!$omp          shared( nsp, na, ngb, eigrb, dfftb, irb, rhocb, &
!$omp                 nat, wrk1, ityp, nabox, iabox ) &
!$omp          private( mytid, ntids, is, ia, iia, ig, qv, fg1 )

      allocate( qv ( dfftb%nnr ) )
      allocate( fg1 ( ngb ) )

!$omp workshare
      wrk1 = (0.d0, 0.d0)
!$omp end workshare
!
#if defined(_OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
#else
      mytid = 0
      ntids = 1
#endif

      do iia = 1, nabox
         IF( MOD( iia - 1, ntids ) == mytid ) THEN
            ia = iabox(iia)
            is = ityp(ia)
            fg1 = eigrb(1:ngb,ia  )*rhocb(1:ngb,is)
            CALL fft_oned2box( qv, fg1 )
            call invfft( qv, dfftb, ia )
            call box2grid(irb(:,ia),1,qv,wrk1)
         END IF
      end do
!
      deallocate( fg1  )
      deallocate( qv  )

!$omp end parallel

      CALL mp_sum( wrk1, inter_bgrp_comm ) 

      call dcopy( dfftp%nnr, wrk1, 2, rhoc, 1 )

      deallocate( wrk1 )
!
      call stop_clock( 'set_cc' )
!
      return
   end subroutine set_cc

