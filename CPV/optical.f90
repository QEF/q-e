!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
      MODULE optical_properties

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        INTEGER :: nfreq
        REAL(DP) :: maxdie ! Hartree units
        REAL(DP) :: ddie  ! Hartree units
        REAL(DP) :: temperature  ! Kelvin
        REAL(DP), PARAMETER :: small_freq = 1.d-6      ! Hartree units

        COMPLEX(DP), ALLOCATABLE :: dielec_total(:)
        REAL(DP), ALLOCATABLE :: sigma_total(:)
        INTEGER, ALLOCATABLE :: n_total(:)

        PUBLIC :: optical_setup, optical_closeup, opticalp, write_dielec
 
      CONTAINS

        SUBROUTINE optical_setup(woptical, noptical, boptical)
          USE constants, ONLY: au
          REAL(DP), INTENT(IN) :: woptical
          REAL(DP), INTENT(IN) :: boptical
          INTEGER, INTENT(IN) :: noptical
          IF( noptical < 1 ) THEN
            CALL errore(' optical_properties: optical_setup ',' noptical out of range ',noptical)
          END IF
          IF( woptical < small_freq ) THEN
            CALL errore(' optical_properties: optical_setup ',' woptical out of range ',INT(woptical))
          END IF
          nfreq = noptical
          maxdie = woptical / au
          ddie = maxdie / DBLE(nfreq)
          temperature = boptical
          ALLOCATE( dielec_total(nfreq), sigma_total(nfreq), n_total(nfreq) )
          dielec_total = 0.0d0
          sigma_total = 0.0d0
          n_total = 0
          RETURN
        END SUBROUTINE optical_setup

        SUBROUTINE optical_closeup
          IF( ALLOCATED( dielec_total ) ) DEALLOCATE( dielec_total )
          IF( ALLOCATED( sigma_total ) ) DEALLOCATE( sigma_total )
          IF( ALLOCATED( n_total ) ) DEALLOCATE( n_total )
          RETURN
        END SUBROUTINE optical_closeup

        SUBROUTINE opticalp(nfi, box, atoms, c0, wfill, occ, ce, wempt, vpot, eigr, bec )

          USE cell_module, ONLY: boxdimensions
          USE wave_types, ONLY: wave_descriptor
          USE pseudo_projector, ONLY: projector, allocate_projector, deallocate_projector
          USE pseudopotential, ONLY: nsanl
          USE nl, ONLY: nlsm1_s
          USE forces, ONLY: dforce_all
          USE brillouin, ONLY: kpoints, kp
          USE electrons_module, ONLY: ei, ei_emp
          USE kohn_sham_states, ONLY: kohn_sham
          USE constants, ONLY: au, pi, k_boltzman_au, au_to_ohmcmm1
          USE cell_base, ONLY: tpiba2
          USE mp, ONLY: mp_sum
          USE mp_global, ONLY: group
          USE io_global, ONLY: ionode
          USE atoms_type_module, ONLY: atoms_type
          USE io_files, ONLY: dielecunit, dielecfile
          USE control_flags, ONLY: force_pairing
          USE reciprocal_vectors, ONLY: gx, g
          USE reciprocal_space_mesh, ONLY: gkx_l
          USE pseudopotential, ONLY: nspnl
          USE uspp_param, ONLY: nhm
          USE uspp, ONLY: nkb

          INTEGER, INTENT(IN) :: nfi
          TYPE(boxdimensions), INTENT(IN) :: box
          TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
          COMPLEX(DP), INTENT(IN) :: c0(:,:,:,:)
          COMPLEX(DP), INTENT(INOUT) :: ce(:,:,:,:)
          TYPE(wave_descriptor), INTENT(IN) :: wempt, wfill
          REAL(DP), INTENT(IN) :: occ(:,:,:)
          REAL (DP), INTENT(in) ::   vpot(:,:)
          REAL (DP) ::   bec(:,:)
          COMPLEX(DP) :: eigr(:,:)

          INTEGER :: nspin, ispin, ngw, nb_l, nk, ngw_g
          INTEGER :: ie, if, nf, ne, ik, idie, ig, ierr
          COMPLEX(DP), ALLOCATABLE :: eforce(:,:,:)
          REAL (DP), ALLOCATABLE ::   bece(:,:)
          REAL(DP) :: curr(3), currt, wef
          COMPLEX(DP) :: ccurr(3), ccurrt
          COMPLEX(DP), ALLOCATABLE :: diet(:), cf(:,:,:,:)
          INTEGER :: cif, cie
          REAL(DP), ALLOCATABLE :: sigma(:), fi(:), eig(:), ff(:,:)
          INTEGER, ALLOCATABLE :: ndiet(:)
          REAL(DP) :: FACT, ef, beta
          LOGICAL :: gamma_symmetry, gzero

! ...     SUBROUTINE BODY

          CALL errore( ' opticalp ', ' not working ', 1 )

          IF( force_pairing ) &
            CALL errore( ' opticalp ', ' force_pairing not implemented ', 1 )

          ALLOCATE( cf( SIZE(c0,1), SIZE(c0,2), SIZE(c0,3), SIZE(c0,4) ) )
          cf = c0

          nk    = wfill%nkl
          nspin = wfill%nspin

          beta = 1.0d0 / ( k_boltzman_au * temperature )

          ALLOCATE( diet(nfreq), ndiet(nfreq), sigma(nfreq) )
          diet  = 0.0d0
          sigma = 0.0d0
          dielec_total = 0.0d0
          sigma_total = 0.0d0
          ndiet = 0

          DO ispin = 1, nspin

            ngw    = wfill%ngwl
            nb_l   = wfill%nbl( ispin )

            ALLOCATE( eforce( ngw,  nb_l, nk ) )
            !
            CALL nlsm1 ( nb_l, 1, nspnl, eigr, cf(1,1,1,ispin), bec )

            CALL dforce_all( ispin, cf(:,:,1,ispin), wfill, occ(:,1,ispin), eforce(:,:,1), &
                             vpot(:,ispin), eigr, bec )

            CALL kohn_sham( ispin, cf(:,:,:,ispin), wfill, eforce )

            ngw  = wempt%ngwl
            nb_l = wempt%nbl( ispin )

            ALLOCATE( ff( nb_l, nk ) )
            DO ik = 1, nk
              ff( 1:nb_l, ik ) = 2.0d0 / nspin
            END DO

            ALLOCATE( bece( nkb, nb_l ) ) 
            !
            CALL nlsm1 ( nb_l, 1, nspnl, eigr, ce(1,1,1,ispin), bece )
            !
            CALL dforce_all( ispin, ce(:,:,1,ispin), wempt, ff( :, 1), eforce(:,:,1), &
                             vpot(:,ispin), eigr, bece )
            !
            CALL kohn_sham( ispin, ce(:,:,:,ispin), wempt, eforce )

            DEALLOCATE( bece )
            DEALLOCATE( eforce )
            DEALLOCATE( ff )

          END DO

          fact = 2.0d0 * pi / ( 3.0d0 * box%deth )
          gamma_symmetry = wfill%gamma

          ierr = 0
          IF( ionode ) THEN
             OPEN(UNIT=dielecunit, FILE=dielecfile, STATUS='unknown', POSITION='append', IOSTAT=ierr)
          END IF
          CALL mp_sum( ierr )
          IF( ierr /= 0 ) &
            CALL errore(' opticalp ', ' opening file '//TRIM(dielecfile), 1 )


#if defined __CONDUCTIVITY

          WRITE( dielecunit, * ) 'STEP: ',nfi

          DO ispin = 1, nspin

            nf   = wfill%nbl( ispin )
            ne   = wempt%nbl( ispin )
            ngw  = wfill%ngwl
            nk   = kp%nkpt

            ALLOCATE( fi( nf + ne ), eig( nf + ne ) )

            DO ik = 1, nk

              ef = ( ei_emp(1,ik,ispin) + ei(nf,ik,ispin) ) / 2.0
              DO if = 1, nf
                fi( if ) = 1.0 / nspin / ( exp( beta * (ei(if,ik,ispin) - ef) ) + 1 )
                eig( if ) = ei( if, ik, ispin )
                IF( ionode ) WRITE( dielecunit, fmt = '(I4,2F12.6)' ) if, fi(if), eig(if)
              END DO
              DO ie = nf+1, ne+nf
                fi( ie ) = 1.0 / nspin / ( exp( beta * (ei_emp(ie-nf,ik,ispin) - ef) ) + 1 )
                eig( ie ) = ei_emp( ie-nf, ik, ispin )
                IF( ionode ) WRITE( dielecunit, fmt = '(I4,2F12.6)' ) ie, fi(ie), eig(ie)
              END DO

              DO if = 1, (nf + ne - 1)
                DO ie = if + 1, (nf + ne)

                  ! frequencies in atomic units
                  wef  = eig(ie) - eig(if)
                  ! discretize the frequency
                  idie = wef / ddie + 1

                  IF( (wef > small_freq) .AND. (idie <= nfreq) ) THEN

                    IF( ie <= nf ) THEN
                      cie = ie
                    ELSE
                      cie = ie-nf
                    END IF

                    IF( if <= nf ) THEN
                      cif = if
                    ELSE
                      cif = if-nf
                    END IF

                    IF( gamma_symmetry ) THEN
                      curr = 0.0d0
                      DO ig = 1, ngw
                        curr(1) = curr(1) + gx(1,ig) * &
                          AIMAG( ce( ig, cie, ik, ispin ) * CONJG( cf( ig, cif, ik, ispin ) )
                        curr(2) = curr(2) + gx(2,ig) * &
                          AIMAG( ce( ig, cie, ik, ispin ) * CONJG( cf( ig, cif, ik, ispin ) )
                        curr(3) = curr(3) + gx(3,ig) * &
                          AIMAG( ce( ig, cie, ik, ispin ) * CONJG( cf( ig, cif, ik, ispin ) )
                      END DO
                      ! parallel sum of curr
                      CALL mp_sum( curr, group ) 
                      ! the factor 4.0d0 accounts for gamma symmetry
                      currt = 4.0d0 * (fi(if)-fi(ie)) * ( curr(1)**2 + curr(2)**2 + curr(3)**2 )
                      currt = currt * tpiba2  / wef 
                      ! update dielectric tensor
                      diet(idie)  = diet(idie)  + CMPLX(0.0d0, currt) / wef
                      sigma(idie) = sigma(idie) + currt
                      ndiet(idie) = ndiet(idie) + 1
                    END IF 
                  END IF
                END DO
              END DO

!              DO if = 1, nf
!                DO ie = 1, ne
!                  ! frequencies in atomic units
!                  wef  = ei_emp(ie,ik,ispin) - ei(if,ik,ispin) 
!                  ! discretize the frequency
!                  idie = wef / ddie + 1
!                  IF( (wef .GT. small_freq) .AND. (idie .LT. nfreq) ) THEN
!                    IF( gamma_symmetry ) THEN
!                      curr = 0.0d0
!                      DO ig = 1, ngw
!                        curr(1) = curr(1) + gx(1,ig) * &
!                          AIMAG( ce(ik,ispin)%w(ig, ie) * CONJG( cf(ik,ispin)%w(ig, if) ) )
!                        curr(2) = curr(2) + gx(2,ig) * &
!                          AIMAG( ce(ik,ispin)%w(ig, ie) * CONJG( cf(ik,ispin)%w(ig, if) ) )
!                        curr(3) = curr(3) + gx(3,ig) * &
!                          AIMAG( ce(ik,ispin)%w(ig, ie) * CONJG( cf(ik,ispin)%w(ig, if) ) )
!                      END DO
!                      ! parallel sum of curr
!                      CALL mp_sum( curr, group ) 
!                      ! the factor 4.0d0 accounts for gamma symmetry
!                      currt = 4.0d0 * ( curr(1)**2 + curr(2)**2 + curr(3)**2 )
!                      currt = currt * tpiba2 ! / wef 
!                      ! update dielectric tensor
!                      diet(idie)  = diet(idie)  + CMPLX(0.0d0, currt) ! / wef
!                      sigma(idie) = sigma(idie) + currt
!                      ndiet(idie) = ndiet(idie) + 1
!                    END IF 
!                  END IF
!                END DO
!              END DO


            END DO
            DEALLOCATE( fi, eig )
          END DO

#else

          WRITE( dielecunit, * ) 'STEP: ',nfi
          DO ispin = 1, nspin
            nf   = wfill%nbl( ispin )
            ne   = wempt%nbl( ispin )
            ngw  = wfill%ngwl
            nk   = kp%nkpt
            DO ik = 1, nk
              DO ie = 1, ne
                DO if = 1, nf
                  wef  = ei_emp(ie,ik,ispin) - ei(if,ik,ispin)
                  IF( gamma_symmetry ) THEN
                    curr = 0.0d0
                    DO ig = 1, ngw
                      curr(1) = curr(1) + gx(1,ig) * &
                          AIMAG( ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) ) )
                      curr(2) = curr(2) + gx(2,ig) * &
                          AIMAG( ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) ) )
                      curr(3) = curr(3) + gx(3,ig) * &
                          AIMAG( ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) ) )
                    END DO
                    CALL mp_sum( curr, group )
                    currt = 2.0d0 * ( curr(1)**2 + curr(2)**2 + curr(3)**2 )
                  ELSE
                    ccurr = 0.0d0
                    DO ig = 1, ngw
                      ccurr(1) = ccurr(1) + gkx_l(1, ig, ik) * &
                          ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) )
                      ccurr(2) = ccurr(2) + gkx_l(2, ig, ik) * &
                          ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) )
                      ccurr(3) = ccurr(3) + gkx_l(3, ig, ik) * &
                          ce( ig, ie, ik, ispin ) * CONJG( cf( ig, if, ik, ispin ) )
                    END DO
                    CALL mp_sum( ccurr, group )
                    ccurrt = ccurr(1)*CONJG(ccurr(1)) + ccurr(2)*CONJG(ccurr(2)) + ccurr(3)*CONJG(ccurr(3))
                    WRITE( dielecunit ,100 ) ispin, ik, ie, if, wef, ccurrt
  100               FORMAT(4I5,1D14.6,3X,2D14.6)
                  END IF
                END DO
              END DO
            END DO
          END DO

#endif

          ierr = 0
          IF( ionode ) THEN
             CLOSE(UNIT=dielecunit, IOSTAT=ierr)
          END IF
          CALL mp_sum( ierr )
          IF( ierr /= 0 ) &
            CALL errore(' opticalp ', ' opening file '//TRIM(dielecfile), 1 )


          ! accumulate statistical values
          WHERE( ndiet > 0 ) 
            dielec_total = fact * diet 
            sigma_total = fact * sigma 
            n_total = ndiet
          END WHERE

          DEALLOCATE( diet, ndiet, sigma )
          DEALLOCATE( cf )

          RETURN
        END SUBROUTINE opticalp


        SUBROUTINE write_dielec (nfi, tm)
          USE constants, ONLY: au, au_to_ohmcmm1
          USE io_files, ONLY: dielecunit, dielecfile
          USE io_global, ONLY: ionode
          USE mp, ONLY: mp_sum

          INTEGER, INTENT(IN) :: nfi
          REAL(DP), INTENT(IN) :: tm
          REAL(DP) :: w
          INTEGER :: i, ierr

#if defined __CONDUCTIVITY

          ierr = 0
          IF( ionode ) THEN
             OPEN(UNIT=dielecunit, FILE=dielecfile, STATUS='unknown', POSITION='append', IOSTAT=ierr)
          END IF
          CALL mp_sum( ierr )
          IF( ierr /= 0 ) &
            CALL errore(' write_dielec ', ' opening file '//TRIM(dielecfile), 1 )

          WRITE( dielecunit, 30 ) nfi, tm
          DO I = 1, SIZE(dielec_total)
            w = (DBLE(i)-0.5d0) * ddie
            ! WRITE(dielecunit,100) &
            !   w * au, dielec_total(i) / w / w, sigma_total(i) * au_to_ohmcmm1 / w, n_total(i)
            WRITE(dielecunit,100) &
              w * au, dielec_total(i), sigma_total(i) * au_to_ohmcmm1, n_total(i)
          END DO

          ierr = 0
          IF( ionode ) THEN
             CLOSE(UNIT=dielecunit, IOSTAT=ierr)
          END IF
          CALL mp_sum( ierr )
          IF( ierr /= 0 ) &
            CALL errore(' write_dielec ', ' closing file '//TRIM(dielecfile), 1 )

! ...       write to stdout
          WRITE( stdout,40) 
          WRITE( stdout,50) 
          DO I = 1, SIZE(dielec_total)
            w = (DBLE(i)-0.5d0) * ddie
            ! WRITE( stdout,110) w * au, sigma_total(i) * au_to_ohmcmm1 / w, n_total(i)
            WRITE( stdout,110) w * au, sigma_total(i) * au_to_ohmcmm1 / ddie, n_total(i)
          END DO

#endif

 30       FORMAT(2X,'STEP:',I7,1X,F10.2)
 40       FORMAT(/,3X,'Frequency dependent electronic conductivity',/)
 50       FORMAT(3X,'frequency (eV)  conductivity ( Ohm^-1 cm^-1 )')
 90       FORMAT(3X,'Dielectric function')
 100      FORMAT(3X,F12.6,3D16.8,I5)
 110      FORMAT(3X,F12.6,D16.8,I5)
          RETURN
        END SUBROUTINE write_dielec
        

      END MODULE optical_properties
