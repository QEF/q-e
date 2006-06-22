!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Wed Apr  5 23:04:18 MDT 2000
!  ----------------------------------------------

#include "f_defs.h"


  MODULE stress

       USE kinds
       USE control_flags, ONLY: timing

       IMPLICIT NONE
       PRIVATE
       SAVE

       PUBLIC :: pstress, print_stress_time
       PUBLIC :: pseudo_stress, compute_gagb, stress_har
       PUBLIC :: stress_hartree, add_drhoph, stress_local

       INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
       INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)

       REAL(DP),  DIMENSION(3,3), PARAMETER :: delta = reshape &
         ( (/ 1.0_DP, 0.0_DP, 0.0_DP, &
              0.0_DP, 1.0_DP, 0.0_DP, &
              0.0_DP, 0.0_DP, 1.0_DP  &
            /), (/ 3, 3 /) )

! ...  dalbe(:) = delta(alpha(:),beta(:))
       REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)

       REAL(DP) :: timeek = 0.0d0
       REAL(DP) :: timeex = 0.0d0
       REAL(DP) :: timeesr = 0.0d0
       REAL(DP) :: timeeh = 0.0d0
       REAL(DP) :: timeel = 0.0d0
       REAL(DP) :: timeenl = 0.0d0
       REAL(DP) :: timetot = 0.0d0
       INTEGER   :: timcnt = 0

       REAL(DP), EXTERNAL  :: cclock

     CONTAINS


!------------------------------------------------------------------------------!

   SUBROUTINE pstress( strvxc, rhoeg, vxc, pail, desr, bec, c0, cdesc, occ, &
                       eigr, sfac, grho, v2xc, box, edft )

      !  this routine computes stress tensor from dft total energy

      USE cell_module,          ONLY: boxdimensions
      USE energies,             ONLY: dft_energy_type
      USE ions_base,            ONLY: nsp
      USE mp_global,            ONLY: intra_image_comm
      USE mp,                   ONLY: mp_sum
      USE wave_types,           ONLY: wave_descriptor
      USE cell_base,            ONLY: tpiba2
      USE io_global,            ONLY: ionode
      USE exchange_correlation, ONLY: stress_xc
      USE control_flags,        ONLY: iprsta
      USE reciprocal_vectors,   ONLY: gx, gstart
      USE gvecp,                ONLY: ngm
      USE gvecs,                ONLY: ngs
      USE local_pseudo,         ONLY: vps, dvps, rhops
      USE atom,                 ONLY: nlcc
      USE core,                 ONLY: drhocg
      USE electrons_base,       ONLY: nspin

      IMPLICIT NONE

      ! ... declare subroutine arguments

      REAL(DP) :: pail(:,:), desr(:), strvxc
      REAL(DP) :: grho(:,:,:), v2xc(:,:,:)
      REAL(DP) :: bec(:,:)
      COMPLEX(DP) :: rhoeg(:,:), vxc(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      REAL(DP), INTENT(IN) :: occ(:,:)
      COMPLEX(DP), INTENT(IN) :: c0(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (boxdimensions), INTENT(IN) :: box
      COMPLEX(DP), INTENT(IN) :: eigr(:,:)
      TYPE (dft_energy_type) :: edft

      ! ... declare other variables

      COMPLEX(DP), ALLOCATABLE :: rhog( : )
      COMPLEX(DP), ALLOCATABLE :: drhog( :, : )
      REAL(DP),    ALLOCATABLE :: gagb(:,:)
      REAL(DP), DIMENSION (6) :: dekin, deht, deps, denl, dexc, dvdw
      REAL(DP), DIMENSION (3,3) :: paiu
      REAL(DP) :: s1, s2, s3, s4, s5, s6, s7, s8, s0
      REAL(DP) :: omega, ehr
      INTEGER  :: k, ig, is

      ! ... end of declarations
      !

      IF( .NOT. cdesc%gamma ) &
        CALL errore( ' pstress ', ' k-point stress not yet implemented ', 1 )

      omega = box%deth
      ehr   = edft%eht - edft%esr + edft%eself

      
      IF(timing) s0 = cclock()

      ALLOCATE( gagb( 6, ngm ) )
      !
      CALL compute_gagb( gagb, gx, ngm, tpiba2 )

      IF(timing) s1 = cclock()

      ! ... compute kinetic energy contribution

      CALL stress_kin( dekin, c0, cdesc, occ, gagb )

      IF(timing) s2 = cclock()

      ! ... compute exchange & correlation energy contribution

      CALL stress_xc(dexc, strvxc, sfac, vxc, grho, v2xc, gagb, nlcc, drhocg, box)

      IF(timing) s3 = cclock()

      ALLOCATE( drhog( ngm, 6 ), rhog( ngm ) )
      
      ! sum up spin components
      !  
      rhog( 1:ngm ) = rhoeg( 1:ngm, 1 )
      IF( nspin > 1 ) rhog( 1:ngm ) = rhog( 1:ngm ) + rhoeg( 1:ngm, 2 )
      
      ! add drho_e / dh       
      !
      DO k = 1, 6
         drhog( 1:ngm, k ) = - rhog( 1:ngm ) * dalbe( k )
      END DO

      CALL stress_local( deps, edft%epseu, gagb, sfac, rhog, drhog, omega )

      !    CALL pseudo_stress(deps, edft%epseu, gagb, sfac, dvps, rhoeg, box%deth)

      IF(timing) s4 = cclock()

      !    IF(tvdw) THEN
      !       CALL vdw_stress(c6, iesr, stau0, dvdw, na, nax, nsp)
      !    END IF

      IF(timing) s5 = cclock()

      ! ... compute hartree energy contribution

      ! add Ionic pseudo charges  rho_I
      !
      DO is = 1, nsp
         DO ig = gstart, ngs
            rhog( ig ) = rhog( ig ) + sfac( ig, is ) * rhops( ig, is )
         END DO
      END DO
      
      ! add drho_I / dh
      !
      CALL add_drhoph( drhog, sfac, gagb )

      CALL stress_hartree(deht, ehr, sfac, rhog, drhog, gagb, omega )

      !    CALL stress_har(deht, ehr, sfac, rhoeg, gagb, box%deth)

      DEALLOCATE( drhog, rhog )

      IF(timing) s6 = cclock()

      ! ... compute enl (non-local) contribution

      CALL stress_nl(denl, gagb, c0, cdesc, occ, eigr, bec, edft%enl, box%a, box%m1 )

      IF(timing) s7 = cclock()

      IF( iprsta >= 2 ) THEN
         CALL stress_debug(dekin, deht, dexc, desr, deps, denl, box%m1 )
      END IF

      ! ... total stress (pai-lowercase)

      DO k = 1, 6
        paiu(alpha(k),beta(k)) = -( dekin(k) + deht(k) + dexc(k) + desr (k) + deps(k) + denl(k) )
        paiu(beta(k),alpha(k)) = paiu(alpha(k),beta(k))
      END DO

      pail(:,:) = matmul( paiu(:,:), box%m1(:,:) )
    
      CALL mp_sum( pail, intra_image_comm )
  
      DEALLOCATE( gagb )

      IF( timing ) THEN
        s8 = cclock()
        timeek  = (s2 - s1) + timeek
        timeeh  = (s3 - s2) + timeeh
        timeex  = (s4 - s3) + timeex
        timeesr = (s5 - s4) + timeesr 
        timeel  = (s6 - s5) + timeel
        timeenl = (s7 - s6) + timeenl
        timetot = (s8 - s0) + timetot
        timcnt = timcnt + 1
      END IF

 50   FORMAT(6X,3(F20.12))
 60   FORMAT(6X,6(F20.12))
100   FORMAT(6X,A3,10X,F8.4)

      RETURN
   END SUBROUTINE pstress


!------------------------------------------------------------------------------!


      SUBROUTINE print_stress_time( iunit )
        USE io_global, ONLY: ionode
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iunit
        IF( timing .AND. timcnt > 0 ) THEN
          timeek  = timeek/timcnt
          timeeh  = timeeh/timcnt
          timeex  = timeex/timcnt
          timeesr = timeesr/timcnt
          timeel  = timeel/timcnt
          timeenl = timeenl/timcnt
          timetot = timetot/timcnt
          IF(ionode) THEN
            WRITE( iunit, 999 ) timeek, timeex, timeesr, timeeh, timeel, timeenl, timetot
          END IF
        END IF
        timeek = 0.0d0 
        timeex = 0.0d0
        timeesr = 0.0d0
        timeeh = 0.0d0 
        timeel = 0.0d0
        timeenl = 0.0d0
        timetot = 0.0d0
        timcnt = 0
999     FORMAT(1X,7(1X,F9.3))

        RETURN
      END SUBROUTINE print_stress_time


!------------------------------------------------------------------------------!


   SUBROUTINE stress_nl(denl, gagb, c0, cdesc, occ, eigr, bec, enl, ht, htm1)

      !  this routine computes nl part of the stress tensor from dft total energy

      USE constants, ONLY: pi, au_gpa
      USE pseudopotential, ONLY: nlin_stress, nlin, nspnl, nsanl
      USE ions_base, ONLY: nsp, na
      USE spherical_harmonics, ONLY: set_dmqm, set_fmrm, set_pmtm
      USE mp_global, ONLY: intra_image_comm, root_image, me_image
      USE io_global, ONLY: stdout
      USE wave_types, ONLY: wave_descriptor
      USE cell_base, ONLY: tpiba2, omega
      USE control_flags, ONLY: force_pairing
      USE reciprocal_vectors, ONLY: gstart, gzero, g, gx
      USE uspp_param, only: nh, lmaxkb, nbeta, nhm
      USE uspp, only: nhtol, nhtolm, indv, nkb
      USE electrons_base, only: nupdwn, iupdwn, nspin
      USE cdvan, only: dbec
      USE cvan, only: ish
      USE mp, only: mp_sum

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP), INTENT(IN) :: occ(:,:)
      COMPLEX(DP), INTENT(IN) :: c0(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      REAL(DP), INTENT(OUT) :: denl(:)
      REAL(DP), INTENT(IN) :: gagb(:,:)
      COMPLEX(DP), INTENT(IN) :: eigr(:,:)
      REAL(DP) :: bec(:,:)
      REAL(DP), INTENT(IN) :: enl
      REAL(DP), INTENT(IN) :: ht(3,3)
      REAL(DP), INTENT(IN) :: htm1(3,3)

! ... declare functions
      REAL(DP)  DDOT

! ... declare other variables
      INTEGER :: is, l, ll, al, be, s, k
      INTEGER :: ir, kk, m, mm, isa, ig, iy, iv, iyy, ih, ihh
      INTEGER :: ia, in, i, iss, nx, ispin, ngw, j, inl
      INTEGER :: iss_wfc, ispin_wfc, mi(16), igh(0:3)
      REAL(DP)  xg,xrg,arg,wnd,wnd1,wnd2,temp,tt1,fac,tt2
      REAL(DP)  temp2, fg, gmod, anm, sgn
      REAL(DP)  pm(3,3), pmtm(6,3,3)
      REAL(DP)  dm(6,5), dmqm(6,5,5)
      REAL(DP)  fm(3,3,3,7), fmrm(6,7,7)
      REAL(DP)  facty(16)
      REAL(DP)  denl_new(3,3), detmp(3,3), denl_gpa(3,3)
      LOGICAL :: new_stress = .false.

      COMPLEX(DP), ALLOCATABLE :: auxc(:,:)
      REAL(DP), ALLOCATABLE :: wnla(:,:,:)
      REAL(DP), ALLOCATABLE :: wnl(:,:,:,:)
      REAL(DP), ALLOCATABLE :: wsg(:,:)
      REAL(DP), ALLOCATABLE :: btmp(:,:,:,:)
      REAL(DP), ALLOCATABLE :: fnls(:,:)
      REAL(DP), ALLOCATABLE :: fnlb(:,:,:,:)
      REAL(DP), ALLOCATABLE :: gspha(:,:)
      REAL(DP), ALLOCATABLE :: gwtmp(:)
      REAL(DP), PARAMETER :: twothird = 2.0d0/3.0d0
      COMPLEX(DP), PARAMETER :: uimag = (0.0d0,1.0d0)
! ... (i)^l
      COMPLEX(DP), PARAMETER :: csign(0:3) = (/ (1.0d0, 0.0d0), &
          (0.0d0,1.0d0), (-1.0d0,0.0d0), (0.0d0, -1.0d0) /)


!  end of declarations
!  ----------------------------------------------


      IF( new_stress ) then

         DO iss = 1, nspin
            !
            iss_wfc = iss
            IF( force_pairing ) iss_wfc = 1
            !
            ALLOCATE( btmp( nkb, nupdwn( iss ), 3, 3 ) )
            ! 
            CALL caldbec( cdesc%ngwl, nkb, nupdwn( iss ), 1, nspnl, eigr(1,1), &
                          c0( 1, 1, iss_wfc ), btmp( 1, 1, 1, 1 ), .false. )
            !
            DO j = 1, 3
               DO i = 1, 3
                  DO in = iupdwn( iss ), iupdwn( iss ) + nupdwn( iss ) - 1
                     dbec( :, in , i, j ) = btmp( :, in - iupdwn( iss ) + 1, i, j )
                  END DO
               END DO
            END DO
            !
            DEALLOCATE( btmp )
            !
         END DO
 
         CALL dennl( bec, denl_new )

      end if


      IF(gzero) THEN
        denl = - enl * dalbe
      ELSE
        denl = 0.0_DP
      END IF

      ngw = cdesc%ngwl
      
      ! ... initialize array wnla
 
      ALLOCATE( wsg ( nhm, nsp ) )
      ALLOCATE( wnl ( ngw, MAXVAL( nbeta( 1:nsp ) ), nsp, 1 ) )
      ALLOCATE( wnla( ngw, MAXVAL( nbeta( 1:nsp ) ), nsp ) )
      ALLOCATE( fnlb( nsanl, MAXVAL( nh( 1:nsp ) ), MAXVAL( nupdwn ), nspin ) )
      !
      fac = sqrt( omega ) / ( 2.0d0 * 4.0d0 * pi )
      !
      DO iss = 1, nspin
         DO in = 1, nupdwn( iss )
            isa = 0
            DO is = 1, nspnl
               DO ia = 1, na( is )
                  DO ih = 1, nh( is )
                     l   = nhtol ( ih, is )
                     inl = ish(is) + (ih-1) * na(is) + ia
                     sgn = 1.0d0
                     IF( MOD( l, 2 ) /= 0 ) sgn = -1.0d0   !  ( -1)^l
                     fnlb( isa + ia, ih, in, iss ) = sgn * fac * bec( inl, in + iupdwn( iss ) - 1 )
                     ! WRITE(6,*) ' i ', iss, in, is, ia, ih, l
                     ! WRITE(6,*) ' v ', fnlb( isa + ia, ih, in, iss ) , fnl(1,iss)%r( isa+ia, ih, in )
                     ! WRITE(6,*) ' r ', fnlb( isa + ia, ih, in, iss ) / fnl(1,iss)%r( isa+ia, ih, in )
                  END DO
               END DO
               isa = isa + na( is )
            END DO
         END DO
      END DO
      !
      CALL nlin( wsg, wnl )
      CALL nlin_stress( wnla )

      ALLOCATE( gwtmp( ngw ) )
      ALLOCATE( gspha( ngw, (lmaxkb+1)**2 ) )

      CALL ylmr2( (lmaxkb+1)**2, ngw, gx, g, gspha )

      DO iy = 1, (lmaxkb+1)**2
        DO ig = gstart, ngw
          gspha(ig,iy) = gspha(ig,iy) / (g(ig)*tpiba2)
        END DO
      END DO

      CALL set_pmtm( pm, pmtm )
      CALL set_dmqm( dm, dmqm )
      CALL set_fmrm( fm, fmrm )

      mi( 1 ) = 1

      mi( 2 ) = 2    !   im( 1 ) = 3
      mi( 3 ) = 3    !   im( 2 ) = 1
      mi( 4 ) = 1    !   im( 3 ) = 2

      mi( 5 ) = 3    !   im( 1 ) = 5
      mi( 6 ) = 4    !   im( 2 ) = 3
      mi( 7 ) = 2    !   im( 3 ) = 1
      mi( 8 ) = 5    !   im( 4 ) = 2
      mi( 9 ) = 1    !   im( 5 ) = 4

      mi( 10 ) = 4    !   im( 1 ) = 7
      mi( 11 ) = 5    !   im( 2 ) = 5
      mi( 12 ) = 3    !   im( 3 ) = 3
      mi( 13 ) = 6    !   im( 4 ) = 1
      mi( 14 ) = 2    !   im( 5 ) = 2
      mi( 15 ) = 7    !   im( 6 ) = 4
      mi( 16 ) = 1    !   im( 7 ) = 6


      SPIN_LOOP: DO ispin = 1, nspin

        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1

        nx = cdesc%nbl( ispin )

        IF( nx < 1 ) CYCLE SPIN_LOOP

        iss = 1

        SPECIES: DO is = 1, nspnl

          ALLOCATE(fnls(na(is),nx))
          ALLOCATE(auxc(ngw,na(is)))

          DO kk = 1, 6
            !
            igh(0:3) = -1

            DO ih = 1, nh( is )
        
              iy  = nhtolm( ih, is )
              iv  = indv  ( ih, is )
              l   = nhtol ( ih, is )
              anm = 2*l + 1


              ! WRITE(6,*) 'DEBUG ih, iy, iv, l = ', ih, iy, iv, l

              gwtmp(1) = 0.0d0
              DO ig = gstart, ngw
                gwtmp( ig ) = gagb(kk,ig) * gspha(ig,iy) * ( anm * wnl(ig,iv,is,1) - wnla(ig,iv,is) )
              END DO

              IF( igh(l) < 0 ) igh(l) = ih
              IF ( l == 1 ) THEN
              ELSE IF ( l == 2 ) THEN
                DO ig = gstart, ngw
                  gwtmp(ig) = gwtmp(ig) - 2.0d0/3.0d0 * dm( kk, mi( iy ) ) * wnl(ig,iv,is,1)
                END DO
              ELSE IF ( l == 3 ) THEN
                al = alpha(kk)
                be = beta(kk)
                DO ig = gstart, ngw
                  fg = 0.0d0
                  gmod = SQRT( g(ig) )
                  DO s = 1, 3
                    fg = fg + 3.0d0/5.0d0 * fm(be,s,s,mi(iy)) * gx(al,ig) / gmod
                  END DO
                  DO s = 1, 3
                    fg = fg + 6.0d0/5.0d0 * fm(be,s,al,mi(iy)) * gx(s,ig) / gmod
                  END DO
                  gwtmp(ig) = gwtmp(ig) - fg * wnl(ig,iv,is,1)
                END DO
              END IF

              DO ihh = igh(l), igh(l) + 2*l
                iyy = nhtolm( ihh, is )
                IF ( l == 0 ) THEN
                  facty( ihh ) = 0.0d0
                ELSE IF( l == 1 ) THEN
                  facty( ihh ) =  pmtm(kk, mi( iy ), mi( iyy ) )
                ELSE IF( l == 2 ) THEN
                  facty( ihh ) =  dmqm(kk, mi( iy ), mi( iyy ) )
                ELSE IF( l == 3 ) THEN
                  facty( ihh ) =  fmrm(kk, mi( iy ), mi( iyy ) )
                END IF
              END DO
              !
              DO ia = 1, na(is)
                auxc(1,ia) = CMPLX(0.0d0,0.0d0)
                DO ig = gstart, ngw
                  auxc(ig,ia) = csign(l) * gwtmp(ig) * eigr(ig,ia+iss-1)
                END DO
              END DO
    
              CALL DGEMM( 'T', 'N', na(is), nx, 2*ngw, 1.0d0, auxc(1,1), &
                2*ngw, c0(1,1,ispin_wfc), 2 * cdesc%ldg, 0.0d0, fnls(1,1), na(is) )

              DO in = 1, nx
                !
                fac = 2.0d0 * occ( in, ispin ) * wsg( ih, is)
                !
                DO ia = 1, na(is)
                  isa = iss + ia - 1
                  temp2 = 0.d0
                  IF( me_image == root_image ) THEN
                    DO ihh = igh(l), igh(l) + 2*l
                      temp2 = temp2 + facty( ihh ) * fnlb( isa, ihh, in, ispin )
                    END DO
                  END IF
                  tt1 = fnlb( isa, ih, in, ispin )
                  tt2 = - l * temp2 + 2.d0 * fnls( ia, in )
                  denl(kk) = denl(kk) + fac * tt1 * tt2
                END DO
              END DO

            END DO

          END DO
          !
          DEALLOCATE(auxc)
          DEALLOCATE(fnls)

          iss = iss + na(is)
          !
        END DO SPECIES


      END DO SPIN_LOOP

      DEALLOCATE( gwtmp )
      DEALLOCATE( gspha )
      DEALLOCATE( fnlb )
      DEALLOCATE( wnla )
      DEALLOCATE( wnl )
      DEALLOCATE( wsg )

      IF( new_stress ) THEN
        DO k=1,6
          detmp(alpha(k),beta(k)) = denl(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        denl_gpa = detmp
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "FROM stress_nl derivative of e(nl)"
        CALL mp_sum( detmp, intra_image_comm )
        CALL mp_sum( denl_new, intra_image_comm )
        CALL mp_sum( denl_gpa, intra_image_comm )
        WRITE( stdout,*) "denl FPMD"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
        WRITE( stdout,*) "denl CP"
        WRITE( stdout,5555) ((denl_new(i,j),j=1,3),i=1,3)
        WRITE( stdout,*) "Stress nl (GPa)"
        WRITE( stdout,5555) ((-denl_gpa(i,j)*au_gpa/omega,j=1,3),i=1,3)
      END IF

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)


      RETURN
      END SUBROUTINE stress_nl



!------------------------------------------------------------------------------!



   SUBROUTINE pseudo_stress( deps, epseu, gagb, sfac, dvps, rhoeg, omega )
      !
      USE ions_base,          ONLY: nsp
      USE reciprocal_vectors, ONLY: gstart
      USE gvecs,              ONLY: ngs
      USE electrons_base,     ONLY: nspin

      REAL(DP),     INTENT(IN)  :: omega
      REAL(DP),     INTENT(OUT) :: deps(:)
      REAL(DP),     INTENT(IN)  :: gagb(:,:)
      COMPLEX(DP),  INTENT(IN)  :: rhoeg(:,:)
      COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
      REAL(DP),     INTENT(IN)  :: dvps(:,:)
      REAL(DP),     INTENT(IN)  :: epseu

      INTEGER     :: ig,k,is, ispin
      COMPLEX(DP) :: rhets, depst(6)
      COMPLEX(DP), ALLOCATABLE :: rhoe( : )
      COMPLEX(DP), ALLOCATABLE :: drhoe( :, : )
      !
      ALLOCATE( drhoe( ngs, 6 ), rhoe( ngs ) )

      rhoe( 1:ngs ) = rhoeg( 1:ngs, 1 )
      IF( nspin > 1 ) rhoe( 1:ngs ) = rhoe( 1:ngs ) + rhoeg( 1:ngs, 2 )

      DO k = 1, 6
         drhoe( 1:ngs, k ) = - rhoe( 1:ngs ) * dalbe( k )
      END DO

      CALL stress_local( deps, epseu, gagb, sfac, rhoe, drhoe, omega )

      DEALLOCATE( drhoe, rhoe )

      RETURN
   END SUBROUTINE pseudo_stress

!  ----------------------------------------------

   SUBROUTINE stress_local( deps, epseu, gagb, sfac, rhoe, drhoe, omega )
      !
      USE ions_base,          ONLY: nsp
      USE reciprocal_vectors, ONLY: gstart
      USE gvecs,              ONLY: ngs
      USE electrons_base,     ONLY: nspin
      USE local_pseudo,       ONLY: vps, dvps

      REAL(DP),     INTENT(IN)  :: omega
      REAL(DP),     INTENT(OUT) :: deps(:)
      REAL(DP),     INTENT(IN)  :: gagb(:,:)
      COMPLEX(DP),  INTENT(IN)  :: rhoe(:)
      COMPLEX(DP),  INTENT(IN)  :: drhoe(:,:)
      COMPLEX(DP),  INTENT(IN)  :: sfac(:,:)
      REAL(DP),     INTENT(IN)  :: epseu

      INTEGER     :: ig,k,is, ispin
      COMPLEX(DP) :: dsvp, svp, depst(6)
      REAL(DP)    :: wz
      !
      depst = (0.d0,0.d0)

      wz = 2.0d0

      DO ig = gstart, ngs
         svp = 0.0d0
         DO is = 1, nsp
            svp = svp + sfac( ig, is ) * vps( ig, is )
         END DO
         depst = depst + wz * CONJG( drhoe( ig, : ) ) * svp
      END DO
      IF( gstart == 2 ) THEN
         svp = 0.0d0
         DO is = 1, nsp
            svp = svp + sfac( 1, is ) * vps( 1, is )
         END DO
         depst = depst + CONJG( drhoe( 1, : ) ) * svp
      END IF

      DO ig = gstart, ngs
         dsvp = 0.0d0
         DO is = 1, nsp
            dsvp = dsvp + sfac( ig, is ) * dvps( ig, is )
         END DO
         DO k = 1, 6
            depst( k ) = depst( k ) - wz * 2.0d0 * CONJG( rhoe( ig ) ) * dsvp * gagb( k, ig )
         END DO
      END DO

      deps = omega * DBLE( depst )

      RETURN
   END SUBROUTINE stress_local


!  ----------------------------------------------


      SUBROUTINE stress_kin(dekin, c0, cdesc, occ, gagb) 

!  this routine computes the kinetic energy contribution to the stress 
!  tensor
!
!  dekin(:) = - 2 (sum over i) f(i) * 
!    ( (sum over g) gagb(:,g) CONJG( c0(g,i) ) c0(g,i) )
!                       
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE gvecw,              ONLY: ecsig, ecfix, ecutz, ngw
      USE wave_types,         ONLY: wave_descriptor
      USE constants,          ONLY: pi
      USE control_flags,      ONLY: force_pairing
      USE reciprocal_vectors, ONLY: gstart, g
      USE cell_base,          ONLY: tpiba2
      USE electrons_base,     ONLY: nspin

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP), INTENT(OUT) :: dekin(:)
      COMPLEX(DP), INTENT(IN) :: c0(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      REAL(DP), INTENT(IN) :: occ(:,:)
      REAL(DP) gagb(:,:)

! ... declare other variables
      REAL(DP)  :: sk(6), scg, efac
      REAL(DP), ALLOCATABLE :: arg(:)
      INTEGER    :: ib, ig, ispin, ispin_wfc

! ... end of declarations
!  ----------------------------------------------

      dekin = 0.0_DP
      ALLOCATE( arg( ngw ) ) 

      efac = 2.0d0 * ecutz / ecsig / SQRT(pi)
      IF( efac > 0.0d0 ) THEN
        DO ig = gstart, ngw
          arg(ig) = 1.0d0 + efac * exp( -( ( tpiba2 * g(ig) - ecfix ) / ecsig )**2 )
        END DO
      ELSE
        arg = 1.0d0
      END IF

      ! ... compute kinetic energy contribution

      DO ispin = 1, nspin
        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1
        DO ib = 1, cdesc%nbl( ispin )
          sk = 0.0_DP
          DO ig = gstart, ngw
            scg = arg(ig) * CONJG( c0(ig,ib,ispin_wfc) ) * c0(ig,ib,ispin_wfc)
            sk(1)  = sk(1) + scg * gagb(1,ig)
            sk(2)  = sk(2) + scg * gagb(2,ig)
            sk(3)  = sk(3) + scg * gagb(3,ig)
            sk(4)  = sk(4) + scg * gagb(4,ig)
            sk(5)  = sk(5) + scg * gagb(5,ig)
            sk(6)  = sk(6) + scg * gagb(6,ig)
          END DO
          dekin = dekin  + occ(ib,ispin) * sk
        END DO
      END DO
      dekin = - 2.0_DP * dekin
      DEALLOCATE(arg) 
      RETURN
      END SUBROUTINE stress_kin


!  ----------------------------------------------

   SUBROUTINE add_drhoph( drhot, sfac, gagb )
      !
      USE kinds,        ONLY: DP
      USE gvecs,        ONLY: ngs
      USE ions_base,    ONLY: nsp, rcmax
      USE local_pseudo, ONLY: rhops
      !
      COMPLEX(DP), INTENT(INOUT) :: drhot( :, : )
      COMPLEX(DP), INTENT(IN) :: sfac( :, : )
      REAL(DP),    INTENT(IN) :: gagb( :, : )
      !
      INTEGER     :: ij, is, ig
      COMPLEX(DP) :: drhop
      !
      DO ij = 1, 6
         IF( dalbe( ij ) > 0.0d0 ) THEN
            DO is = 1, nsp
               DO ig = 1, ngs
                  drhot(ig,ij) = drhot(ig,ij) - sfac(ig,is)*rhops(ig,is)
               ENDDO
            END DO
         END IF
      END DO
      DO ig = 1, ngs
         drhop = 0.0d0
         DO is = 1, nsp
           drhop = drhop - sfac( ig, is ) * rhops(ig,is) * rcmax(is)**2 * 0.5D0
         END DO
         DO ij = 1, 6
             drhot(ig,ij) = drhot(ig,ij) - drhop * gagb( ij, ig )
         END DO
      END DO
      RETURN
   END SUBROUTINE add_drhoph

!  ----------------------------------------------


   SUBROUTINE stress_har(deht, ehr, sfac, rhoeg, gagb, omega ) 

      use ions_base,          only: nsp, rcmax
      USE cell_module,        only: boxdimensions
      use mp_global,          ONLY: me_image, root_image
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecs,              ONLY: ngs
      USE gvecp,              ONLY: ngm
      USE local_pseudo,       ONLY: rhops
      USE electrons_base,     ONLY: nspin

      IMPLICIT NONE

      REAL(DP)    :: omega, DEHT(:), EHR, gagb(:,:)
      COMPLEX(DP) :: RHOEG(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)

      COMPLEX(DP)    DEHC(6)
      COMPLEX(DP)    RHOP,DRHOP
      COMPLEX(DP)    RHET,RHOG,RHETS,RHOGS
      COMPLEX(DP)    CFACT
      COMPLEX(DP), ALLOCATABLE :: rhot(:), drhot(:,:)
      REAL(DP)       hgm1

      INTEGER       ig, is, k, ispin


      ALLOCATE( rhot( ngm ) )
      ALLOCATE( drhot( ngm, 6 ) )

      ! sum up spin components
      !
      DO ig = gstart, ngm
         rhot( ig ) = rhoeg( ig, 1 )
         IF( nspin > 1 ) rhot( ig ) = rhot( ig ) + rhoeg( ig, 2 )
      END DO
      !
      ! add Ionic pseudo charges  rho_I
      !
      DO is = 1, nsp
         DO ig = gstart, ngs
            rhot( ig ) = rhot( ig ) + sfac( ig, is ) * rhops( ig, is )
         END DO
      END DO

      ! add drho_e / dh
      !
      DO k = 1, 6
         IF( dalbe( k ) > 0.0d0 ) THEN
            drhot( :, k ) = - rhoeg( :, 1 )
            IF( nspin > 1 ) drhot( :, k ) = drhot( :, k ) + rhoeg( :, 2 )
         ELSE
            drhot( :, k ) = 0.0d0
         END IF
      END DO

      ! add drho_I / dh
      !
      CALL add_drhoph( drhot, sfac, gagb )

      CALL stress_hartree(deht, ehr, sfac, rhot, drhot, gagb, omega ) 

      DEALLOCATE( rhot, drhot )

      RETURN
   END SUBROUTINE stress_har

!  ----------------------------------------------

   SUBROUTINE stress_hartree(deht, ehr, sfac, rhot, drhot, gagb, omega ) 

      ! This subroutine computes: d E_hartree / dh  =
      !   E_hartree * h^t + 
      !   4pi omega rho_t * CONJG( rho_t ) / G^2 / G^2 * G_alpha * G_beta +
      !   4pi omega Re{ CONJG( rho_t ) * drho_t / G^2 }
      ! where:
      !   rho_t  = rho_e + rho_I
      !   drho_t = d rho_t / dh = -rho_e + d rho_hard / dh  + d rho_I / dh

      use ions_base,          only: nsp, rcmax
      USE cell_module,        only: boxdimensions
      use mp_global,          ONLY: me_image, root_image
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecs,              ONLY: ngs
      USE gvecp,              ONLY: ngm
      USE local_pseudo,       ONLY: rhops
      USE electrons_base,     ONLY: nspin

      IMPLICIT NONE

      REAL(DP)    :: omega, DEHT(:), EHR, gagb(:,:)
      COMPLEX(DP) :: rhot(:)  ! total charge: Sum_spin ( rho_e + rho_I )
      COMPLEX(DP) :: drhot(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)

      COMPLEX(DP)    DEHC(6)
      COMPLEX(DP)    CFACT
      REAL(DP), ALLOCATABLE :: hgm1( : )
      REAL(DP)    :: wz

      INTEGER       ig, is, k, iss

      DEHC  = (0.D0,0.D0)
      DEHT  = 0.D0

      wz = 2.0d0

      ALLOCATE( hgm1( ngm ) )

      hgm1( 1 ) = 0.0d0
      DO ig = gstart, ngm
         hgm1( ig ) = 1.D0 / g(ig) / tpiba2
      END DO

      ! Add term  rho_t * CONJG( rho_t ) / G^2 * G_alpha * G_beta / G^2

      DO ig = gstart, ngm
         cfact = rhot( ig ) * CONJG( rhot( ig ) ) * hgm1( ig ) ** 2 
         dehc = dehc + cfact * gagb(:,ig)
      END DO

      ! Add term  2 * Re{ CONJG( rho_t ) * drho_t / G^2 }

      DO ig = gstart, ngm
         DO k = 1, 6
            dehc( k ) = dehc( k ) +  rhot( ig ) * CONJG( drhot( ig, k ) ) * hgm1( ig )
         END DO
      END DO

      ! term:  E_h * h^t

      if ( me_image == root_image ) then
        deht = wz * fpi * omega * DBLE(dehc) + ehr * dalbe
      else
        deht = wz * fpi * omega * DBLE(dehc)
      end if

      DEALLOCATE( hgm1 )

      RETURN
      END SUBROUTINE stress_hartree

!  ----------------------------------------------
!  ----------------------------------------------
!  ----------------------------------------------


      SUBROUTINE stress_debug(dekin, deht, dexc, desr, deps, denl, htm1)
        USE io_global, ONLY: stdout
        USE mp_global, ONLY: intra_image_comm
        USE mp,        ONLY: mp_sum
        REAL(DP) :: dekin(:), deht(:), dexc(:), desr(:), deps(:), denl(:)
        REAL(DP) :: detot( 6 ), htm1(3,3)
        REAL(DP) :: detmp(3,3)

        INTEGER :: k, i, j

        detot = dekin + deht + dexc + desr + deps + denl

        DO k=1,6
          detmp(alpha(k),beta(k)) = detot(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(tot)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = dekin(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(kin)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deht(k) + desr(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(electrostatic)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deht(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(h)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = desr(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(sr)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = deps(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(ps)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = denl(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(nl)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

        DO k=1,6
          detmp(alpha(k),beta(k)) = dexc(k)
          detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        END DO
        CALL mp_sum( detmp, intra_image_comm )
        detmp = MATMUL( detmp(:,:), htm1(:,:) )
        WRITE( stdout,*) "derivative of e(xc)"
        WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

      RETURN
      END SUBROUTINE stress_debug

!=----------------------------------------------------------------------------=!


      SUBROUTINE compute_gagb( gagb, gx, ngm, tpiba2 )

         ! ... compute G_alpha * G_beta  

         USE kinds, ONLY: DP
         INTEGER,  INTENT(IN)  :: ngm
         REAL(DP), INTENT(IN)  :: gx(:,:)
         REAL(DP), INTENT(OUT) :: gagb(:,:)
         REAL(DP), INTENT(IN)  :: tpiba2

         INTEGER :: k, ig
      
         DO k = 1, 6
            DO ig = 1, ngm          
               gagb( k, ig ) = gx( alpha( k ), ig ) * gx( beta( k ), ig ) * tpiba2
            END DO
         END DO

      END SUBROUTINE compute_gagb


!=----------------------------------------------------------------------------=!

  END MODULE stress

!=----------------------------------------------------------------------------=!
