!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define DEVICEATTR ,DEVICE
#else
#define DEVICEATTR
#endif

!-----------------------------------------------------------------------
   SUBROUTINE rhoofr_cp &
      ( nfi, c_bgrp, c_d, bec_bgrp, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, &
        enl, denl, ekin, dekin, tstress, ndwwf )
!-----------------------------------------------------------------------
!
!  this routine computes:
!  rhor  = normalized electron density in real space
!  ekin  = kinetic energy
!  dekin = kinetic energy term of QM stress
!
!    rhor(r) = (sum over ib) fi(ib) |psi(r,ib)|^2
!
!    Using quantities in scaled space
!    rhor(r) = rhor(s) / Omega
!    rhor(s) = (sum over ib) fi(ib) |psi(s,ib)|^2 
!
!    fi(ib) = occupation numbers
!    psi(r,ib) = psi(s,ib) / SQRT( Omega ) 
!    psi(s,ib) = INV_FFT (  c0(ig,ib)  )
!
!    ib = index of band
!    ig = index of G vector
!  ----------------------------------------------
!     the normalized electron density rhor in real space
!     the kinetic energy ekin
!     subroutine uses complex fft so it computes two ft's
!     simultaneously
!
!     rho_i,ij = sum_n < beta_i,i | psi_n >< psi_n | beta_i,j >
!     < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                   2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
!
!     e_v = sum_i,ij rho_i,ij d^ion_is,ji
!
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: iprint, iverbosity, thdyn, tpre, trhor, ndr
      USE ions_base,          ONLY: nat
      USE gvect,              ONLY: gstart, ig_l2g
      USE smallbox_gvec,      ONLY: ngb
      USE uspp_param,         ONLY: nh, nhm
      USE cell_base,          ONLY: omega
      USE electrons_base,     ONLY: nspin, nbsp_bgrp, ispin_bgrp, f_bgrp, f_d
      USE constants,          ONLY: pi, fpi
      USE mp,                 ONLY: mp_sum
      USE io_global,          ONLY: stdout, ionode
      USE mp_global,          ONLY: intra_bgrp_comm, nbgrp, inter_bgrp_comm, &
                                    me_bgrp, nproc_bgrp, root_bgrp
      USE xc_lib,             ONLY: xclib_dft_is
      USE cg_module,          ONLY: tcg
      USE cp_interfaces,      ONLY: stress_kin, enkin
      USE fft_interfaces,     ONLY: fwfft, invfft
      USE fft_base,           ONLY: dffts, dfftp
      USE cp_interfaces,      ONLY: checkrho, ennl, calrhovan, dennl
      USE cp_main_variables,  ONLY: iprint_stdout, idesc, irb, eigrb
      USE wannier_base,       ONLY: iwf
      USE exx_module,         ONLY: rhopr 
      USE input_parameters,   ONLY: tcpbo ! BS
      USE io_base,            ONLY: read_rhog
      USE io_files,           ONLY: restart_dir
      USE fft_rho
      USE fft_helper_subroutines, ONLY: c2psi_gamma
      USE mp,                 ONLY: mp_barrier
      USE mp_world,           ONLY: mpime, world_comm
      !
      IMPLICIT NONE
      INTEGER nfi
      REAL(DP) bec_bgrp(:,:)
      REAL(DP) dbec(:,:,:,:)
      REAL(DP) rhovan(:, :, : )
      REAL(DP) rhor(:,:)
      REAL(DP) drhor(:,:,:,:)
      REAL(DP) rhos(:,:)
      REAL(DP) enl, ekin
      REAL(DP) denl(3,3), dekin(6)
      COMPLEX(DP) rhog( :, : )
      COMPLEX(DP) drhog( :, :, :, : )
      COMPLEX(DP) c_bgrp( :, : )
      COMPLEX(DP) DEVICEATTR :: c_d( :, : )
      LOGICAL, OPTIONAL, INTENT(IN) :: tstress
      INTEGER, OPTIONAL, INTENT(IN) :: ndwwf

      ! local variables

      INTEGER  :: iss, isup, isdw, iss1, iss2, i, ir, ig, k
      REAL(DP) :: rsumr(2), rsumg(2), sa1, sa2, detmp(6), mtmp(3,3)
      REAL(DP) :: rnegsum, rmin, rmax, rsum
      COMPLEX(DP) :: ci,fp,fm
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: psis, drhovan
#endif
#endif
      COMPLEX(DP), ALLOCATABLE :: psis(:)
      REAL(DP), ALLOCATABLE :: drhovan(:,:,:,:,:)
      CHARACTER(LEN=320) :: filename

      LOGICAL, SAVE :: first = .TRUE.
      LOGICAL :: ttstress
      
      !
      CALL start_clock( 'rhoofr' )

      ttstress = tpre
      IF( PRESENT( tstress ) ) ttstress = tstress

      ci = ( 0.0d0, 1.0d0 )

      rhor = 0.d0
      rhos = 0.d0
      rhog = (0.d0, 0.d0)
      !
      !  calculation of kinetic energy ekin
      !
#if defined (__CUDA)
      ekin = enkin( c_d, f_d, nbsp_bgrp )
#else
      ekin = enkin( c_bgrp, f_bgrp, nbsp_bgrp )
#endif
      !
      IF( nbgrp > 1 ) &
         CALL mp_sum( ekin, inter_bgrp_comm )
      !
      IF( ttstress ) THEN
         !
         ! ... compute kinetic energy contribution
         !
         CALL stress_kin( dekin, c_bgrp, f_bgrp )
         !
         IF( nbgrp > 1 ) &
            CALL mp_sum( dekin, inter_bgrp_comm )
         !
      END IF

      IF( PRESENT( ndwwf ) ) THEN
         !
         !     called from WF, compute only of rhovan
         !
         CALL calrhovan( rhovan, bec_bgrp, iwf )
         !
      ELSE
         !
         !     calculation of non-local energy
         !
         CALL ennl( enl, rhovan, bec_bgrp )
         !
         IF( nbgrp > 1 ) THEN
            CALL mp_sum( enl, inter_bgrp_comm )
            CALL mp_sum( rhovan, inter_bgrp_comm )
         END IF
         !
      END IF
      !
      IF( ttstress ) THEN
         !
         ALLOCATE( drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 ) )
         !
         CALL dennl( bec_bgrp, dbec, drhovan, denl, idesc ) 
         !
         IF( nbgrp > 1 ) THEN
            CALL mp_sum( denl, inter_bgrp_comm )
            CALL mp_sum( drhovan, inter_bgrp_comm )
         END IF
         !
      END IF
      !    
      !    warning! trhor and thdyn are not compatible yet!   
      !
      COMPUTE_CHARGE: IF( trhor .AND. ( .NOT. thdyn ) ) THEN
         !
         !     ==================================================================
         !     non self-consistent charge: charge density is read from unit 47
         !     ==================================================================
         !
         ! Lingzhu Kong
         !
         ! FIXME: in non-scf calculations, the charge density must be read at
         ! FIXME: the beginning, the potential computed and no longer updated.
         !
         IF( first ) THEN
            CALL errore('rhoofr','option trhor unverified, please report',1)
            filename = TRIM( restart_dir(ndr) ) // "charge-density"
            CALL read_rhog ( filename, root_bgrp, intra_bgrp_comm, &
                 ig_l2g, nspin, rhog )
            !
            !^^ ... TEMPORARY FIX  (newlsda) ...
            IF ( nspin==2 ) THEN
               rhog(:,1) = ( rhog(:,1) + rhog(:,2) )*0.5d0
               rhog(:,2) = rhog(:,1) - rhog(:,2)
            ENDIF
            !^^.......................
            !
            CALL rho_g2r ( dfftp, rhog, rhor )
            rhopr = rhor
            first = .FALSE.
         ELSE
            rhor = rhopr
         END IF

         CALL rho_r2g( dfftp, rhor, rhog )

      ELSE
         !
         !     ==================================================================
         !     self-consistent charge
         !     ==================================================================
         !
         !     important: if n is odd then nx must be .ge.n+1 and c(*,n+1)=0.
         ! 

         IF ( MOD( nbsp_bgrp, 2 ) /= 0 ) THEN
            !
            IF( SIZE( c_bgrp, 2 ) < nbsp_bgrp + 1 ) &
               CALL errore( ' rhoofr ', ' c second dimension too small ', SIZE( c_bgrp, 2 ) )
            !
            c_bgrp( :, nbsp_bgrp + 1 ) = ( 0.d0, 0.d0 )
#if defined (__CUDA)
            c_d( :, nbsp_bgrp + 1 ) = ( 0.d0, 0.d0 )
#endif
            !
         ENDIF
         !
         IF( PRESENT( ndwwf ) ) THEN
            !
            ! Wannier function, charge density from state iwf
            !
            ALLOCATE( psis( dffts%nnr ) ) 
            !
            CALL c2psi_gamma( dffts, psis, c_bgrp(:,iwf) )
            !
            CALL invfft('Wave',psis, dffts )
            !
            rhos(1:dffts%nnr,1) = rhos(1:dffts%nnr,1) + f_bgrp(iwf) * ( DBLE(psis(:)))**2 / omega
            !
            DEALLOCATE( psis )
            !
         ELSE 
            !
#if defined (__CUDA)
            CALL loop_over_states_gpu()
#else
            CALL loop_over_states()
#endif
            !
         END IF
         !
         !     smooth charge in g-space is put into rhog(ig)
         !
         CALL rho_r2g( dffts, rhos, rhog )
         !
         rhog(dffts%ngm+1:,:) = 0.0d0
         !
         CALL rho_g2r( dfftp, rhog, rhor )
         !
         IF ( xclib_dft_is('meta') ) THEN
            CALL kedtauofr_meta( c_bgrp ) ! METAGGA
         END IF
         !
         !     add vanderbilt contribution to the charge density
         !     drhov called before rhov because input rho must be the smooth part
         !
         IF ( ttstress ) THEN
            CALL drhov( irb, eigrb, rhovan, drhovan, rhog, rhor, drhog, drhor )
            DEALLOCATE( drhovan )
         END IF
         !
         CALL rhov( rhovan, rhog, rhor )

      ENDIF COMPUTE_CHARGE
!
      IF( PRESENT( ndwwf ) ) THEN
         !
         CALL errore('cp_rhoofr','old_write_rho no longer implemented',1)
         !
      END IF
!
!     here to check the integral of the charge density
!
! BS: I have turned off computing and printing integrated electronic density at
! every iprint_stdout steps during CP-BO calculations ... 
!     IF( ( iverbosity > 1 ) .OR. ( nfi == 0 ) .OR. &
!         ( MOD(nfi, iprint_stdout) == 0 ) .AND. ( .NOT. tcg ) ) THEN
      IF( ( iverbosity > 1 ) .OR. ( nfi == 0 ) .OR. &
          ( MOD(nfi, iprint_stdout) == 0 ) .AND. ( .NOT. tcg ) .AND. (.NOT.tcpbo )) THEN

         IF( iverbosity > 1 ) THEN
            CALL checkrho( dfftp%nnr, nspin, rhor, rmin, rmax, rsum, rnegsum )
            rnegsum = rnegsum * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
            rsum    = rsum    * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
            WRITE( stdout,'(a,4(1x,f12.6))')                                     &
     &     ' rhoofr: rmin rmax rnegsum rsum  ',rmin,rmax,rnegsum,rsum
         END IF

         CALL sum_charge( rsumg, rsumr )

         IF ( nspin == 1 ) THEN
           WRITE( stdout, 10) rsumg(1), rsumr(1)
         ELSE
           WRITE( stdout, 20) rsumg(1), rsumr(1), rsumg(2), rsumr(2)
         ENDIF

      ENDIF

10    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 )
20    FORMAT( /, 3X, 'from rhoofr: total integrated electronic density', &
            & /, 3X, 'spin up', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 , &
            & /, 3X, 'spin down', &
            & /, 3X, 'in g-space = ', f13.6, 3x, 'in r-space =', f13.6 )
!
      CALL stop_clock( 'rhoofr' )

!
      RETURN


   CONTAINS   
      !
      !
      SUBROUTINE sum_charge( rsumg, rsumr )

         !
         REAL(DP), INTENT(OUT) :: rsumg( : )
         REAL(DP), INTENT(OUT) :: rsumr( : )
         INTEGER :: iss
         !
         DO iss=1,nspin
            rsumg(iss)=omega*DBLE(rhog(1,iss))
            rsumr(iss)=SUM(rhor(:,iss),1)*omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         END DO

         IF ( gstart .NE. 2 ) THEN
            ! in the parallel case, only one processor has G=0 !
            rsumg( 1:nspin ) = 0.0d0
         END IF

         CALL mp_sum( rsumg( 1:nspin ), intra_bgrp_comm )
         CALL mp_sum( rsumr( 1:nspin ), intra_bgrp_comm )

         RETURN
      END SUBROUTINE

      !

      SUBROUTINE loop_over_states
         !
         USE parallel_include
         USE fft_helper_subroutines
         !
         !        MAIN LOOP OVER THE EIGENSTATES
         !           - This loop is also parallelized within the task-groups framework
         !           - Each group works on a number of eigenstates in parallel
         !
         IMPLICIT NONE
         !
         INTEGER :: from, i, eig_index, eig_offset, ii, tg_nr3
         !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: tmp_rhos
#endif
#endif
         REAL(DP), ALLOCATABLE :: tmp_rhos(:,:)

         IF( fftx_ntgrp(dffts) > 1 ) THEN
            CALL errore('rhoofr','Task group not supported',1)
         END IF

         ALLOCATE( psis( dffts%nnr_tg ) ) 
         !
         CALL tg_get_group_nr3( dffts, tg_nr3 )
         !
         ALLOCATE( tmp_rhos ( dffts%nr1x * dffts%nr2x * tg_nr3, nspin ) )
         !
         tmp_rhos = 0_DP

         do i = 1, nbsp_bgrp, 2 * fftx_ntgrp(dffts)

            CALL c2psi_gamma( dffts, psis, c_bgrp(:,i), c_bgrp(:,i+1) )

            CALL invfft('Wave', psis, dffts )
            !
            ! Now the first proc of the group holds the first two bands
            ! of the 2*nogrp bands that we are processing at the same time,
            ! the second proc. holds the third and fourth band
            ! and so on
            !
            ! Compute the proper factor for each band
            !
            ii=dffts%mype2+1
            !
            ! Remember two bands are packed in a single array :
            ! proc 0 has bands ibnd   and ibnd+1
            ! proc 1 has bands ibnd+2 and ibnd+3
            ! ....
            !
            ii = 2 * ii - 1

            IF( ii + i - 1 < nbsp_bgrp ) THEN
               iss1=ispin_bgrp( ii + i - 1 )
               sa1 =f_bgrp( ii + i - 1 )/omega
               iss2=ispin_bgrp( ii + i )
               sa2 =f_bgrp( ii + i )/omega
            ELSE IF( ii + i - 1 == nbsp_bgrp ) THEN
               iss1=ispin_bgrp( ii + i - 1 )
               sa1 =f_bgrp( ii + i - 1 )/omega
               iss2=iss1
               sa2=0.0d0
            ELSE
               iss1=ispin_bgrp( nbsp_bgrp )
               sa1 = 0.0d0
               iss2=iss1
               sa2 =0.0d0
            END IF
            !
            !Compute local charge density
            !
            !This is the density within each orbital group...so it
            !coresponds to 1 eignestate for each group and there are
            !NOGRP such groups. Thus, during the loop across all
            !occupied eigenstates, the total charge density must me
            !accumulated across all different orbital groups.
            !

            !This loop goes through all components of charge density that is local
            !to each processor. In the original code this is nnr. In the task-groups
            !code this should be equal to the total number of planes
            !
            do ir = 1, dffts%nr1x*dffts%nr2x*tg_nr3
               tmp_rhos(ir,iss1) = tmp_rhos(ir,iss1) + sa1*( real(psis(ir)))**2
               tmp_rhos(ir,iss2) = tmp_rhos(ir,iss2) + sa2*(aimag(psis(ir)))**2
            end do
            !
         END DO

         IF( nbgrp > 1 ) THEN
            CALL mp_sum( tmp_rhos, inter_bgrp_comm )
         END IF

         CALL tg_reduce_rho( rhos, tmp_rhos, dffts )

         DEALLOCATE( tmp_rhos )
         DEALLOCATE( psis ) 

         RETURN
      END SUBROUTINE loop_over_states

#if defined (__CUDA)

      SUBROUTINE loop_over_states_gpu
         !
         USE parallel_include
         USE fft_helper_subroutines
         USE control_flags, ONLY : many_fft
         USE cudafor
         !
         !        MAIN LOOP OVER THE EIGENSTATES
         !           - This loop is also parallelized within the task-groups framework
         !           - Each group works on a number of eigenstates in parallel
         !
         IMPLICIT NONE
         !
         INTEGER :: from, i, ig, eig_index, eig_offset, ii, tg_nr3, ioff
         !
         REAL(DP), DEVICE, ALLOCATABLE    :: rhos_d(:,:)
         COMPLEX(DP), DEVICE, ALLOCATABLE :: psis(:)
         COMPLEX(DP), DEVICE, ALLOCATABLE :: ptmp(:,:)
         INTEGER,     DEVICE, POINTER     :: nl_d(:), nlm_d(:)

         ALLOCATE( psis( dffts%nnr * many_fft ) )  ! dffts%nnr * many_fft
         ALLOCATE( rhos_d ( SIZE(rhos,1), SIZE(rhos,2) ) )
         !
         rhos_d = 0_DP
         nl_d => dffts%nl_d
         nlm_d => dffts%nlm_d

         do i = 1, nbsp_bgrp, 2 * many_fft

            psis = 0.0d0

            ioff = 0
            DO ii = i, i + 2 * many_fft - 1, 2
              IF( ii < nbsp_bgrp ) THEN
!$cuf kernel do(1)
                 do ig = 1, dffts%ngw
                    psis( nlm_d( ig ) + ioff) = CONJG( c_d( ig, ii ) ) + ci * conjg( c_d( ig, ii+1 ))
                    psis( nl_d( ig )  + ioff) = c_d( ig, ii ) + ci * c_d( ig, ii+1 )
                 end do
              ELSE IF( ii == nbsp_bgrp ) THEN
!$cuf kernel do(1)
                 do ig = 1, dffts%ngw
                    psis( nlm_d( ig ) + ioff) = CONJG( c_d( ig, ii ) )
                    psis( nl_d( ig )  + ioff) = c_d( ig, ii )
                 end do
              END IF
              ! CALL c2psi_gamma( dffts, psis, c_bgrp(:,ii), c_bgrp(:,ii+1) )
              ioff = ioff + dffts%nnr
              
            END DO

            CALL invfft('Wave', psis, dffts, many_fft )

            ioff = 0
            DO ii = i, i + 2 * many_fft - 1, 2
              IF( ii < nbsp_bgrp ) THEN
                iss1=ispin_bgrp( ii )
                sa1 =f_bgrp( ii )/omega
                iss2=ispin_bgrp( ii + 1 )
                sa2 =f_bgrp( ii + 1 )/omega
!$cuf kernel do(1)
                do ir = 1, dffts%nnr
                   rhos_d(ir,iss1) = rhos_d(ir,iss1) + sa1*( real(psis(ir + ioff)))**2
                   rhos_d(ir,iss2) = rhos_d(ir,iss2) + sa2*(aimag(psis(ir + ioff)))**2
                end do
              ELSE IF( ii == nbsp_bgrp ) THEN
                iss1=ispin_bgrp( ii )
                sa1 =f_bgrp( ii )/omega
                iss2=iss1
                sa2=0.0d0
!$cuf kernel do(1)
                do ir = 1, dffts%nnr
                   rhos_d(ir,iss1) = rhos_d(ir,iss1) + sa1*( real(psis(ir + ioff)))**2
                END DO
              END IF
              ioff = ioff + dffts%nnr
            END DO
            !
         END DO

         rhos = rhos_d

         IF( nbgrp > 1 ) THEN
            CALL mp_sum( rhos, inter_bgrp_comm )
         END IF

         DEALLOCATE( rhos_d )
         DEALLOCATE( psis ) 
         DEALLOCATE( ptmp ) 

         RETURN
      END SUBROUTINE loop_over_states_gpu
#endif


!-----------------------------------------------------------------------
   END SUBROUTINE rhoofr_cp
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------
   SUBROUTINE checkrho_x(nnr,nspin,rhor,rmin,rmax,rsum,rnegsum)
!----------------------------------------------------------------------
!
!     check \int rho(r)dr and the negative part of rho
!
      USE kinds,     ONLY: DP
      USE mp,        ONLY: mp_sum
      USE mp_global, ONLY: intra_bgrp_comm

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nnr, nspin
      REAL(DP) rhor(nnr,nspin), rmin, rmax, rsum, rnegsum
      !
      REAL(DP) roe
      INTEGER ir, iss
!
      rsum   =0.0d0
      rnegsum=0.0d0
      rmin   =100.d0
      rmax   =0.0d0 
      DO iss = 1, nspin
         DO ir = 1, nnr
            roe  = rhor(ir,iss)
            rsum = rsum + roe
            IF ( roe < 0.0d0 ) rnegsum = rnegsum + roe
            rmax = MAX( rmax, roe )
            rmin = MIN( rmin, roe )
         END DO
      END DO
      CALL mp_sum( rsum, intra_bgrp_comm )
      CALL mp_sum( rnegsum, intra_bgrp_comm )
      RETURN
   END SUBROUTINE checkrho_x


!-----------------------------------------------------------------------
SUBROUTINE drhov(irb,eigrb,rhovan,drhovan,rhog,rhor,drhog,drhor)
!-----------------------------------------------------------------------
!     this routine calculates arrays drhog drhor, derivatives wrt h of:
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     Same logic as in routine rhov.
!     On input rhor and rhog must contain the smooth part only !!!
!     Output in (drhor, drhog)
!
      USE kinds,                    ONLY: DP
      USE control_flags,            ONLY: iprint
      USE ions_base,                ONLY: na, nsp, nat, ityp
      USE uspp_param,               ONLY: nhm, nh, upf
      USE uspp,                     ONLY: nkbus
      USE electrons_base,           ONLY: nspin
      USE smallbox_gvec,            ONLY: ngb
      USE smallbox_subs,            ONLY: fft_oned2box, box2grid
      USE cell_base,                ONLY: ainv
      USE qgb_mod,                  ONLY: qgb, dqgb
      USE fft_interfaces,           ONLY: fwfft, invfft
      USE fft_base,                 ONLY: dfftb, dfftp
      USE mp_global,                ONLY: my_bgrp_id, nbgrp, inter_bgrp_comm
      USE mp,                       ONLY: mp_sum
      USE fft_helper_subroutines,   ONLY: fftx_add_threed2oned_gamma

      IMPLICIT NONE
! input
      INTEGER,     INTENT(IN) ::  irb(3,nat)
      REAL(DP),    INTENT(IN) ::  rhor(dfftp%nnr,nspin)
      REAL(DP),    INTENT(IN) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      REAL(DP),    INTENT(IN) ::  drhovan(nhm*(nhm+1)/2,nat,nspin,3,3)
      COMPLEX(DP), INTENT(IN) ::  eigrb(ngb,nat), rhog(dfftp%ngm,nspin)
! output
      REAL(DP),    INTENT(OUT) :: drhor(dfftp%nnr,nspin,3,3)
      COMPLEX(DP), INTENT(OUT) :: drhog(dfftp%ngm,nspin,3,3)
! local
      INTEGER i, j, isup, isdw, iv, jv, ig, ijv, is, iss, ia, ir, ijs, itid
      REAL(DP) :: asumt, dsumt
      COMPLEX(DP) fp, fm, ci
      COMPLEX(DP), ALLOCATABLE :: v(:)
      COMPLEX(DP), ALLOCATABLE:: dqgbt(:,:)
      COMPLEX(DP), ALLOCATABLE :: qv(:)
      COMPLEX(DP), ALLOCATABLE :: fg1(:), fg2(:)
#if defined(_OPENMP)
      INTEGER  :: omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
!
!$omp parallel do collapse(3) default(none), private(i,j,iss,ir,ig), shared(nspin,dfftp,drhor,drhog,rhor,rhog,ainv) 
      DO j=1,3
         DO i=1,3
            DO iss=1,nspin
               DO ir=1,dfftp%nnr
                  drhor(ir,iss,i,j)=-rhor(ir,iss)*ainv(j,i)
               END DO
               DO ig=1,dfftp%ngm
                  drhog(ig,iss,i,j)=-rhog(ig,iss)*ainv(j,i)
              END DO
            END DO
         END DO
      END DO
!$omp end parallel do

      IF ( nkbus <= 0 ) THEN
         GO TO 1000
      END IF

      ALLOCATE( v( dfftp%nnr ) )

      ci =( 0.0d0, 1.0d0 )

      IF( nspin == 1 ) THEN
         !  
         !  nspin=1 : two fft at a time, one per atom, if possible
         ! 
         DO i=1,3
            DO j=1,3

               v(:) = (0.d0, 0.d0)

!$omp parallel default(none) &
!$omp          shared(nat, ityp, ngb, nh, eigrb, dfftb, irb, v, &
!$omp                 i, j, dqgb, qgb, nhm, rhovan, drhovan, upf ) &
!$omp          private( is, ia, iv, jv, ijv, ig, iss, &
!$omp                  qv, fg1, fg2, itid, dqgbt, dsumt, asumt )

               ALLOCATE( qv( dfftb%nnr ) )
               ALLOCATE( dqgbt( ngb, 2 ) )
               ALLOCATE( fg1( ngb ) )
               ALLOCATE( fg2( ngb ) )

               itid  = 0
               iss=1

               DO ia=1,nat
                  is=ityp(ia)
                  IF( upf(is)%tvanp ) THEN

#if defined(__MPI)
                     IF ( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) ) THEN
                        CYCLE
                     END IF
#endif

#if defined(_OPENMP)
                     IF ( omp_get_thread_num() /= itid ) THEN
                        itid = MOD( itid + 1, omp_get_num_threads() )
                        CYCLE
                     ELSE
                        itid = MOD( itid + 1, omp_get_num_threads() )
                     END IF
#endif

                     dqgbt(:,:) = (0.d0, 0.d0) 
                     DO iv=1,nh(is)
                        DO jv=iv,nh(is)
                           ijv = (jv-1)*jv/2 + iv
                           IF(iv.NE.jv) THEN
                              asumt = 2.0d0 *  rhovan( ijv, ia, iss )
                              dsumt = 2.0d0 * drhovan( ijv, ia, iss, i, j )
                           ELSE
                              asumt =  rhovan( ijv, ia, iss )
                              dsumt = drhovan( ijv, ia, iss, i, j )
                           ENDIF
                           DO ig=1,ngb
                              dqgbt(ig,1)=dqgbt(ig,1) + asumt*dqgb(ig,ijv,is,i,j)
                              dqgbt(ig,1)=dqgbt(ig,1) + dsumt*qgb(ig,ijv,is)
                           END DO
                        END DO
                     END DO
                     !     
                     ! add structure factor
                     !
                     fg1 = eigrb(1:ngb,ia   )*dqgbt(1:ngb,1)
                     CALL fft_oned2box( qv, fg1 )
                     !
                     CALL invfft( qv, dfftb, ia )
                     !
                     !  qv = US contribution in real space on box grid
                     !       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
                     !
                     !  add qv(r) to v(r), in real space on the dense grid
                     !
                     CALL box2grid( irb(:,ia), 1, qv, v )
                  END IF
               END DO

               DEALLOCATE( fg1 )
               DEALLOCATE( fg2 )
               DEALLOCATE( dqgbt )
               DEALLOCATE( qv )
!
!$omp end parallel


               iss = 1

               DO ir=1,dfftp%nnr
                  drhor(ir,iss,i,j) = drhor(ir,iss,i,j) + DBLE(v(ir))
               END DO
!
               CALL fwfft( 'Rho', v, dfftp )
               CALL fftx_add_threed2oned_gamma( dfftp, v, drhog(:,iss,i,j) )
!
            ENDDO
         ENDDO
!
      ELSE
         !
         !     nspin=2: two fft at a time, one for spin up and one for spin down
         ! 
         isup=1
         isdw=2
         DO i=1,3
            DO j=1,3
               v(:) = (0.d0, 0.d0)
               ALLOCATE( qv( dfftb%nnr ) )
               ALLOCATE( dqgbt( ngb, 2 ) )
               ALLOCATE( fg1( ngb ) )
               ALLOCATE( fg2( ngb ) )
               DO ia=1,nat
                  is=ityp(ia)
                  IF( upf(is)%tvanp ) THEN
#if defined(__MPI)
                     IF ( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) ) CYCLE
#endif
                     DO iss=1,2
                        dqgbt(:,iss) = (0.d0, 0.d0)
                        DO iv= 1,nh(is)
                           DO jv=iv,nh(is)
                              ijv = (jv-1)*jv/2 + iv
                              asumt=rhovan(ijv,ia,iss)
                              dsumt =drhovan(ijv,ia,iss,i,j)
                              IF(iv.NE.jv) THEN
                                 asumt =2.d0*asumt
                                 dsumt=2.d0*dsumt
                              ENDIF
                              DO ig=1,ngb
                                 dqgbt(ig,iss)=dqgbt(ig,iss)  +         &
     &                               (asumt*dqgb(ig,ijv,is,i,j) +         &
     &                               dsumt*qgb(ig,ijv,is))
                              END DO
                           END DO
                        END DO
                     END DO
                     !     
                     ! add structure factor
                     !
                     fg1 = eigrb(1:ngb,ia)*dqgbt(1:ngb,1)
                     fg2 = eigrb(1:ngb,ia)*dqgbt(1:ngb,2)
                     CALL fft_oned2box( qv, fg1, fg2 )

                     CALL invfft(qv, dfftb, ia )
                     !
                     !  qv is the now the US augmentation charge for atomic species is
                     !  and atom ia: real(qv)=spin up, imag(qv)=spin down
                     !
                     !  add qv(r) to v(r), in real space on the dense grid
                     !
                     CALL box2grid(irb(:,ia),qv,v)
                     !
                  END IF
               END DO

               DEALLOCATE( dqgbt )
               DEALLOCATE( qv )
               DEALLOCATE( fg1 )
               DEALLOCATE( fg2 )
!
               DO ir=1,dfftp%nnr
                  drhor(ir,isup,i,j) = drhor(ir,isup,i,j) + DBLE(v(ir))
                  drhor(ir,isdw,i,j) = drhor(ir,isdw,i,j) +AIMAG(v(ir))
               ENDDO
!
               CALL fwfft('Rho', v, dfftp )
               CALL fftx_add_threed2oned_gamma( dfftp, v, drhog(:,isup,i,j), drhog(:,isdw,i,j) )

            END DO
         END DO
      ENDIF

      DEALLOCATE( v )
!
1000  CONTINUE
      RETURN
END SUBROUTINE drhov

!
!-----------------------------------------------------------------------
SUBROUTINE rhov(rhovan,rhog,rhor)
!-----------------------------------------------------------------------
!     Add Vanderbilt contribution to rho(r) and rho(g)
!
!        n_v(g) = sum_i,ij rho_i,ij q_i,ji(g) e^-ig.r_i
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      USE kinds,                    ONLY: dp
      USE ions_base,                ONLY: nat, na, nsp, ityp
      USE io_global,                ONLY: stdout
      USE mp_global,                ONLY: intra_bgrp_comm, inter_bgrp_comm
      USE mp,                       ONLY: mp_sum
      USE uspp_param,               ONLY: nh, nhm, upf
      USE uspp,                     ONLY: deeq, nkbus
      USE electrons_base,           ONLY: nspin
      USE smallbox_gvec,            ONLY: ngb
      USE smallbox_subs,            ONLY: fft_oned2box, box2grid
      USE cell_base,                ONLY: omega
      USE small_box,                ONLY: omegab
      USE control_flags,            ONLY: iprint, iverbosity, tpre
      USE qgb_mod,                  ONLY: qgb
      USE fft_interfaces,           ONLY: fwfft, invfft
      USE fft_base,                 ONLY: dfftb, dfftp, dfftb
      USE fft_helper_subroutines,   ONLY: fftx_add_threed2oned_gamma
      USE cp_main_variables,        ONLY: irb, eigrb, iabox, nabox
!
      IMPLICIT NONE
      !
      REAL(DP),    INTENT(IN) ::  rhovan(nhm*(nhm+1)/2,nat,nspin)
      ! 
      REAL(DP),     INTENT(inout):: rhor(dfftp%nnr,nspin)
      COMPLEX(DP),  INTENT(inout):: rhog(dfftp%ngm,nspin)
!
      INTEGER, PARAMETER :: isup = 1
      INTEGER, PARAMETER :: isdw = 2
      INTEGER     :: iv, jv, ig, ijv, is, iss, ia, ir, i, j, iia
      REAL(DP)    :: sumrho
      COMPLEX(DP) :: fp, fm, ca
      COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: qgbt, v, qv
#endif
#endif
      COMPLEX(DP), ALLOCATABLE :: qgbt(:,:)
      COMPLEX(DP), ALLOCATABLE :: v(:)
      COMPLEX(DP), ALLOCATABLE :: qv(:)
      COMPLEX(DP), ALLOCATABLE :: fg1(:), fg2(:)

      INTEGER  :: mytid, ntids
#if defined(_OPENMP)
      INTEGER  :: omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif

      !  Quick return if this sub is not needed
      !
      IF ( nkbus <= 0 ) THEN
         GO TO 1000
      END IF

      CALL start_clock( 'rhov' )
!
!
      ALLOCATE( v( dfftp%nnr ) )

      ! private variable need to be initialized, otherwise
      ! outside the parallel region they have an undetermined value
      !
      mytid = 0
      ntids = 1
!
      IF(nspin.EQ.1) THEN
         ! 
         !     nspin=1 : two fft at a time, one per atom, if possible
         !

!$omp parallel default(none) &
!$omp          shared(na, ngb, nh, rhovan, qgb, eigrb, dfftb, iverbosity, omegab, irb, v, &
!$omp                 stdout, rhor, dfftp, upf, nsp, ityp, nat, nspin, iabox, nabox, inter_bgrp_comm ) &
!$omp          private(mytid, ntids, is, ia, iia, iv, jv, ijv, sumrho, qgbt, ig, ca, qv, ir )

!$omp workshare
         v (:) = (0.d0, 0.d0)
!$omp end workshare

#if defined(_OPENMP)
         mytid = omp_get_thread_num()  ! take the thread ID
         ntids = omp_get_num_threads() ! take the number of threads
#endif

         ALLOCATE( qgbt( ngb, nspin ) )
         ALLOCATE( qv( dfftb%nnr ) )

         DO iia = 1, nabox
            IF( MOD( iia - 1, ntids ) == mytid ) THEN
               ia = iabox(iia)
               is = ityp(ia)
               qgbt(:,1) = (0.d0, 0.d0)
               DO iv= 1,nh(is)
                  DO jv=iv,nh(is)
                     ijv = (jv-1)*jv/2 + iv
                     sumrho=rhovan(ijv,ia,1)
                     IF(iv.NE.jv) sumrho=2.d0*sumrho
                     DO ig=1,ngb
                        qgbt(ig,1)=qgbt(ig,1) + sumrho*qgb(ig,ijv,is)
                     END DO
                  END DO
               END DO
               !
               ! add structure factor
               !
               qgbt(1:ngb,1) = eigrb(1:ngb,ia)*qgbt(1:ngb,1)
               CALL fft_oned2box( qv, qgbt(:,1) )

               CALL invfft( qv, dfftb, ia )
               !
               !  qv = US augmentation charge in real space on box grid
               !       for atomic species is, real(qv)=atom ia, imag(qv)=atom ia+1
               !
               !  add qv(r) to v(r), in real space on the dense grid
               !
               CALL  box2grid(irb(:,ia),1,qv,v)
            END IF
         END DO

         DEALLOCATE(qv)
         DEALLOCATE(qgbt)
         !
         !  rhor(r) = total (smooth + US) charge density in real space
         !
!$omp barrier
!$omp master
         !
         CALL mp_sum( v, inter_bgrp_comm )
         !
!$omp end master
!$omp barrier
!$omp do 
         DO ir=1,dfftp%nnr
            rhor(ir,1)=rhor(ir,1)+DBLE(v(ir))        
         END DO
!$omp end do
!$omp end parallel

         CALL fwfft('Rho',v, dfftp )
         !
         !  rhog(g) = total (smooth + US) charge density in G-space
         !
         CALL fftx_add_threed2oned_gamma( dfftp, v, rhog(:,1) )

      ELSE
         !
         !     nspin=2: two fft at a time, one for spin up and one for spin down
         !
!$omp parallel default(none) &
!$omp          shared(na, ngb, nh, rhovan, qgb, eigrb, dfftb, iverbosity, omegab, irb, v, &
!$omp                 stdout, rhor, dfftp, upf, nsp, ityp, nat, nspin, iabox, nabox, inter_bgrp_comm ) &
!$omp          private(mytid, ntids, is, ia, iia, iv, jv, ijv, sumrho, qgbt, ig, ca, qv, fg1, fg2, ir )


!$omp workshare
         v (:) = (0.d0, 0.d0)
!$omp end workshare

#if defined(_OPENMP)
         mytid = omp_get_thread_num()  ! take the thread ID
         ntids = omp_get_num_threads() ! take the number of threads
#endif

         ALLOCATE( qgbt( ngb, 2 ) )
         ALLOCATE( qv( dfftb%nnr ) )
         ALLOCATE( fg1( ngb ) )
         ALLOCATE( fg2( ngb ) )

         DO iia = 1, nabox
            IF( MOD( iia - 1, ntids ) == mytid ) THEN
               ia = iabox(iia)
               is = ityp(ia)
               DO iss=1,2
                  qgbt(:,iss) = (0.d0, 0.d0)
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        sumrho=rhovan(ijv,ia,iss)
                        IF(iv.NE.jv) sumrho=2.d0*sumrho
                        DO ig=1,ngb
                           qgbt(ig,iss)=qgbt(ig,iss)+sumrho*qgb(ig,ijv,is)
                        END DO
                     END DO
                  END DO
               END DO
               !     
               ! add structure factor
               !
               fg1 = eigrb(1:ngb,ia)*qgbt(1:ngb,1)
               fg2 = eigrb(1:ngb,ia)*qgbt(1:ngb,2)
               CALL fft_oned2box( qv, fg1, fg2 )
               !
               CALL invfft( qv,dfftb,ia)
               !
               !  qv is the now the US augmentation charge for atomic species is
               !  and atom ia: real(qv)=spin up, imag(qv)=spin down
               !
               !  add qv(r) to v(r), in real space on the dense grid
               !
               CALL box2grid(irb(:,ia),qv,v)

            END IF
         END DO
         DEALLOCATE(qgbt)
         DEALLOCATE( qv )
         DEALLOCATE( fg1 )
         DEALLOCATE( fg2 )
!$omp barrier
!$omp master
         CALL mp_sum( v, inter_bgrp_comm )
!$omp end master
!$omp barrier
!$omp do 
         DO ir=1,dfftp%nnr
            rhor(ir,isup)=rhor(ir,isup)+DBLE(v(ir)) 
            rhor(ir,isdw)=rhor(ir,isdw)+AIMAG(v(ir)) 
         END DO
!$omp end do 
!$omp end parallel
!
         IF( iverbosity > 1 ) THEN
            ca = SUM(v)
            CALL mp_sum( ca, intra_bgrp_comm )
            WRITE( stdout,'(a,2f12.8)') 'rhov:in n_v  ',omega*ca/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         ENDIF
!
         CALL fwfft('Rho',v, dfftp )
!
         IF( iverbosity > 1 ) THEN
            WRITE( stdout,*) 'rhov: smooth up',omega*rhog(1,isup)
            WRITE( stdout,*) 'rhov: smooth dw',omega*rhog(1,isdw)
            WRITE( stdout,*) 'rhov: vander up',omega*DBLE(v(1))
            WRITE( stdout,*) 'rhov: vander dw',omega*AIMAG(v(1))
            WRITE( stdout,*) 'rhov: all up',                                  &
     &           omega*(rhog(1,isup)+DBLE(v(1)))
            WRITE( stdout,*) 'rhov: all dw',                                  &
     &           omega*(rhog(1,isdw)+AIMAG(v(1)))
         ENDIF
!
         CALL fftx_add_threed2oned_gamma( dfftp, v, rhog(:,isup), rhog(:,isdw) )
!
         IF( iverbosity > 1 ) THEN
            WRITE( stdout,'(a,f12.8,/,a,f12.8)')                   &
     &        ' rhov: n_v(g=0) up   = ',omega*DBLE (rhog(1,isup)), &
     &        ' rhov: n_v(g=0) down = ',omega*DBLE(rhog(1,isdw))
         END IF
!
      ENDIF

      DEALLOCATE( v )

      CALL stop_clock( 'rhov' )

1000  CONTINUE
!
      RETURN

CONTAINS

      SUBROUTINE print_rhov()
           IF( nspin == 1 ) THEN
               IF( iverbosity > 1 ) THEN
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom g-sp = ',         &
     &                 omegab*DBLE(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: 1-atom r-sp = ',         &
     &                 omegab*DBLE(ca)/(dfftb%nr1*dfftb%nr2*dfftb%nr3)
               ENDIF
           ELSE
               IF( iverbosity > 1 ) THEN
                  ca = SUM(qv)
                  WRITE( stdout,'(a,f12.8)') ' rhov: up   g-space = ',        &
     &                 omegab*DBLE(qgbt(1,1))
                  WRITE( stdout,'(a,f12.8)') ' rhov: up r-sp = ',             &
     &                 omegab*DBLE(ca)/(dfftb%nr1*dfftb%nr2*dfftb%nr3)
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw g-space = ',          &
     &                 omegab*DBLE(qgbt(1,2))
                  WRITE( stdout,'(a,f12.8)') ' rhov: dw r-sp = ',             &
     &                 omegab*AIMAG(ca)/(dfftb%nr1*dfftb%nr2*dfftb%nr3)
               ENDIF
           ENDIF
      END SUBROUTINE
END SUBROUTINE rhov

SUBROUTINE rhoofr_host &
      ( nfi, c_bgrp, irb, eigrb, bec_bgrp, dbec, rhovan, rhor, drhor, rhog, drhog, rhos, &
        enl, denl, ekin, dekin, tstress, ndwwf )
         USE kinds,         ONLY: DP
#if defined (__CUDA)
         USE cudafor
#endif
         USE cp_interfaces
         IMPLICIT NONE
         INTEGER nfi
         COMPLEX(DP) c_bgrp( :, : )
         INTEGER irb( :, : )
         COMPLEX(DP) eigrb( :, : )
         REAL(DP) bec_bgrp(:,:)
         REAL(DP) dbec(:,:,:,:)
         REAL(DP) rhovan(:, :, : )
         REAL(DP) rhor(:,:)
         REAL(DP) drhor(:,:,:,:)
         COMPLEX(DP) rhog( :, : )
         COMPLEX(DP) drhog( :, :, :, : )
         REAL(DP) rhos(:,:)
         REAL(DP) enl, ekin
         REAL(DP) denl(3,3), dekin(6)
         LOGICAL, OPTIONAL, INTENT(IN) :: tstress
         INTEGER, OPTIONAL, INTENT(IN) :: ndwwf
         COMPLEX(DP), ALLOCATABLE DEVICEATTR :: c(:,:) 
         ALLOCATE( c, SOURCE=c_bgrp )
         CALL rhoofr(nfi, c_bgrp, c, bec_bgrp, dbec, rhovan, rhor, &
              drhor, rhog, drhog, rhos, enl, denl, ekin, dekin, tstress, ndwwf )
         DEALLOCATE( c ) 
END SUBROUTINE rhoofr_host
