!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"



MODULE optical_properties


        USE kinds, ONLY: DP

        IMPLICIT NONE
        SAVE

        PRIVATE

        INTEGER  :: nfreq
        REAL(DP) :: maxdie = 10.d0                  ! Hartree units
        REAL(DP) :: ddie   = 0.01d0                 ! Hartree units
        REAL(DP) :: temperature = 300.0d0           ! Kelvin


        PUBLIC ::  opticalp
 

CONTAINS


   SUBROUTINE opticalp(nfi, omega, c0, wfill, ce, wempt, vpot, eigr, vkb, bec )

          USE wave_types,         ONLY: wave_descriptor
          USE pseudopotential,    ONLY: nsanl
          USE cp_interfaces,      ONLY: dforce, kohn_sham
          USE electrons_module,   ONLY: ei, ei_emp, n_emp, nupdwn_emp, iupdwn_emp, nb_l, n_emp_l
          USE electrons_base,     ONLY: nupdwn, iupdwn, nspin
          USE constants,          ONLY: autoev, pi, k_boltzmann_au, au_to_ohmcmm1, eps8
          USE cell_base,          ONLY: tpiba2
          USE mp,                 ONLY: mp_sum
          USE mp_global,          ONLY: intra_image_comm
          USE io_global,          ONLY: ionode, stdout
          USE control_flags,      ONLY: force_pairing
          USE reciprocal_vectors, ONLY: gx, g
          USE uspp_param,         ONLY: nhm
          USE uspp,               ONLY: nkb
          USE read_pseudo_module_fpmd, ONLY: nspnl

          INTEGER,  INTENT(IN) :: nfi
          REAL(DP), INTENT(IN) :: omega
          COMPLEX(DP), INTENT(IN) :: c0(:,:)
          COMPLEX(DP), INTENT(INOUT) :: ce(:,:)
          TYPE(wave_descriptor), INTENT(IN) :: wempt, wfill
          REAL (DP), INTENT(in) ::   vpot(:,:)
          REAL (DP) ::   bec(:,:)
          COMPLEX(DP) :: eigr(:,:)
          COMPLEX(DP) :: vkb(:,:)

          INTEGER :: iss, ngw
          INTEGER :: ie, if, nf, ne, idie, ig, ierr, i, nbt, nsst
          COMPLEX(DP), ALLOCATABLE :: eforce(:,:)
          REAL (DP), ALLOCATABLE ::   bece(:,:)
          REAL(DP) :: curr(3), currt, wef, w, wg2, p2
          COMPLEX(DP) :: ccurr(3), ccurrt
          COMPLEX(DP), ALLOCATABLE :: diet(:), cf(:,:)
          INTEGER :: cif, cie
          REAL(DP), ALLOCATABLE :: sigma(:), fi(:), eig(:), ff(:), dipole(:,:,:,:)
          REAL(DP), ALLOCATABLE :: epsilon2(:), dos(:)
          INTEGER, ALLOCATABLE :: ndiet(:)
          REAL(DP) :: FACT, ef, beta, sumrule(3)


          ! ...     SUBROUTINE BODY

          IF( force_pairing ) &
            CALL errore( ' opticalp ', ' force_pairing not implemented ', 1 )

          nfreq = INT( maxdie / ddie + 1.0d0 )
          !
          beta = 1.0d0 / ( k_boltzmann_au * temperature )
          !
          fact = 2.0d0 * pi / ( 3.0d0 * omega )

          ngw  = wfill%ngwl

          nbt = wfill%nbl( 1 ) + wempt%nbl( 1 )
          IF( nspin > 1 ) nbt = MAX( nbt, wfill%nbl( 2 ) + wempt%nbl( 2 ) )

          ALLOCATE( diet(nfreq), ndiet(nfreq), sigma(nfreq) )
          ALLOCATE( cf( SIZE(c0,1), SIZE(c0,2) ) )
          ALLOCATE( dipole( 3, nbt, nbt, nspin ) )
          ALLOCATE( epsilon2( nfreq ) )
          ALLOCATE( dos( nfreq ) )

          diet  = 0.0d0
          sigma = 0.0d0
          ndiet = 0
          dipole = 0.0d0
          epsilon2 = 0.0d0
          dos = 0.0d0

          cf = c0

          ALLOCATE( eforce( ngw,  MAX( SIZE( c0, 2 ), SIZE( ce, 2 ) ) ) )

          nsst = nupdwn( 1 )
          IF( nspin == 2 ) nsst = nsst + nupdwn( 2 )
          CALL nlsm1 ( nsst, 1, nspnl, eigr, cf, bec )
          !
          nsst = nupdwn_emp( 1 )
          IF( nspin == 2 ) nsst = nsst + nupdwn_emp( 2 )
          ALLOCATE( bece( nkb, nsst ) ) 
          CALL nlsm1 ( nsst, 1, nspnl, eigr, ce, bece )

          DO iss = 1, nspin

             ALLOCATE( ff( MAX( nupdwn( iss ), nupdwn_emp( iss ) ) ) )

             ff = 2.0d0 / nspin
             !
             CALL dforce( cf, ff, eforce, vpot(:,iss), vkb, bec, nupdwn(iss), iupdwn(iss) )

             CALL kohn_sham( cf, ngw, eforce, nupdwn(iss), nb_l(iss), iupdwn(iss) )
             !
             CALL dforce( ce, ff, eforce, vpot(:,iss), vkb, bece, nupdwn_emp(iss), iupdwn_emp(iss) )
             !
             CALL kohn_sham( ce, ngw, eforce, nupdwn_emp(iss), n_emp_l(iss), iupdwn_emp(iss) )

             DEALLOCATE( ff )

          END DO

          DEALLOCATE( bece )
          DEALLOCATE( eforce )


          WRITE( stdout, * ) 'epsilon: '

          sumrule = 0.0d0

          DO iss = 1, nspin

            nf   = wfill%nbl( iss )
            ne   = wempt%nbl( iss )

            ALLOCATE( fi( nf + ne ), eig( nf + ne ) )

            ef = ( ei_emp(1,iss) + ei(nf,iss) ) / 2.0

            DO if = 1, nf
               fi( if ) = 2.0 / nspin
               eig( if ) = ei( if, iss )
               ! IF( ionode ) WRITE( stdout, fmt = '(I4,2F12.6)' ) if, fi(if), eig(if)
            END DO

            DO ie = nf+1, ne+nf
               fi( ie ) = 0.0d0
               eig( ie ) = ei_emp( ie-nf, iss )
               ! IF( ionode ) WRITE( stdout, fmt = '(I4,2F12.6)' ) ie, fi(ie), eig(ie)
            END DO

            DO if = 1, nf
               !
               DO ie = nf + 1, (nf + ne)

                  ! frequencies in atomic units
                  !
                  wef  = eig(ie) - eig(if)
                  !
                  ! discretize the frequency
                  !
                  idie = wef / ddie + 1

                  IF( wef > eps8 ) THEN

                      cie = ie-nf
                      cif = if

                      ccurr = 0.0d0

                      DO ig = 1, ngw
                        ccurr(1) = ccurr(1) + gx(1,ig) * &
                           ce( ig, cie + iupdwn_emp(iss) - 1 ) * CONJG( cf( ig, cif + iupdwn(iss) - 1 ) )
                        ccurr(2) = ccurr(2) + gx(2,ig) * &
                           ce( ig, cie + iupdwn_emp(iss) - 1 ) * CONJG( cf( ig, cif + iupdwn(iss) - 1 ) )
                        ccurr(3) = ccurr(3) + gx(3,ig) * &
                           ce( ig, cie + iupdwn_emp(iss) - 1 ) * CONJG( cf( ig, cif + iupdwn(iss) - 1 ) )
                      END DO
                      !
                      CALL mp_sum( ccurr, intra_image_comm ) 
                      !
                      ! the factor 4.0d0 accounts for gamma symmetry
                      !
                      wg2 = 4.0d0
                      !
                      curr = AIMAG( ccurr )
                      !
                      !dipole( :, if, ie, iss ) = wg2 * tpiba2 * DBLE( ccurr(:) * CONJG( ccurr(:) ) )
                      !dipole( :, ie, if, iss ) = wg2 * tpiba2 * DBLE( ccurr(:) * CONJG( ccurr(:) ) )
                      dipole( :, if, ie, iss ) = wg2 * tpiba2 * curr(:)**2
                      dipole( :, ie, if, iss ) = wg2 * tpiba2 * curr(:)**2 
                      !
                      p2 = DBLE( dipole( 1, if, ie, iss ) + dipole( 2, if, ie, iss ) + dipole( 3, if, ie, iss ) )
                      !
                      !
                      currt = wg2 * (fi(if)-fi(ie)) * ( curr(1)**2 + curr(2)**2 + curr(3)**2 )
                      currt = currt * tpiba2  / wef 
                      !
                      ! update dielectric tensor
                      !
                      IF( idie <= nfreq ) THEN
                          
                         diet(idie)  = diet(idie)  + CMPLX(0.0d0, currt) / wef
                         sigma(idie) = sigma(idie) + currt
                         ndiet(idie) = ndiet(idie) + 1
                         epsilon2(idie) = epsilon2(idie) + 4.0d0 * pi * pi * p2 * fi(if) / wef**2 / 3.0d0
                      END IF

                      sumrule = sumrule + fi(if) * 2.0d0 * dipole( :, if, ie, iss ) / wef

                  END IF
                  !
               END DO
               !
            END DO

            DEALLOCATE( fi, eig )

          END DO


          diet  = fact * diet 
          sigma = fact * sigma 

          ! ...       write to stdout

          IF( ionode ) THEN

             WRITE( stdout,*) '  sumrule = ', sumrule(1)
             WRITE( stdout,*) '  sumrule = ', sumrule(2)
             WRITE( stdout,*) '  sumrule = ', sumrule(3)
             WRITE( stdout,*) '  2p/omega  = ', 2.0d0 * pi / omega
             WRITE( stdout,*) '  sumrule = ', sumrule(1) * 2.0d0 * pi / omega
             WRITE( stdout,*) '  sumrule = ', sumrule(2) * 2.0d0 * pi / omega
             WRITE( stdout,*) '  sumrule = ', sumrule(3) * 2.0d0 * pi / omega
             DO i = 1, nfreq
                IF( ndiet(i) > 0 ) THEN
                   w = (DBLE(i)-0.5d0) * ddie
                   WRITE( 100, * ) w * autoev, epsilon2(i) / omega
                END IF
             END DO
             !WRITE( stdout,40) 
             !WRITE( stdout,50) 
             !DO I = 1, SIZE(diet)
             !  w = (DBLE(i)-0.5d0) * ddie
             !  WRITE( stdout,110) w * autoev, sigma(i) * au_to_ohmcmm1 / ddie, ndiet(i)
             !  ! WRITE(stdout,100) &
             !  !   w * autoev, diet(i) / w / w, sigma(i) * au_to_ohmcmm1 / w, ndiet(i)
             !END DO

          END IF



          DEALLOCATE( diet, ndiet, sigma )
          DEALLOCATE( cf )
          DEALLOCATE( dipole )
          DEALLOCATE( dos )

 30       FORMAT(2X,'STEP:',I7,1X,F10.2)
 40       FORMAT(/,3X,'Frequency dependent electronic conductivity',/)
 50       FORMAT(3X,'frequency (eV)  conductivity ( Ohm^-1 cm^-1 )')
 90       FORMAT(3X,'Dielectric function')
 100      FORMAT(3X,F12.6,3D16.8,I5)
 110      FORMAT(3X,F12.6,D16.8,I5)

          RETURN

   END SUBROUTINE opticalp



END MODULE optical_properties
