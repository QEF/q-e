!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE energies
        
        USE io_global,  ONLY : stdout
        USE kinds
        USE funct,      ONLY : dft_is_hybrid, get_exx_fraction
        IMPLICIT NONE
        SAVE

        PRIVATE

        TYPE dft_energy_type
          REAL(DP)  :: ETOT
          REAL(DP)  :: SKIN
          REAL(DP)  :: EMKIN
          REAL(DP)  :: EHT
          REAL(DP)  :: EH
          REAL(DP)  :: SELF_EHTE
          REAL(DP)  :: EHTE
          REAL(DP)  :: EHTI
          REAL(DP)  :: EPSEU
          REAL(DP)  :: ENL
          REAL(DP)  :: ENT
          REAL(DP)  :: exx              
          REAL(DP)  :: VXC
          REAL(DP)  :: EXC
          REAL(DP)  :: SELF_VXC
          REAL(DP)  :: SELF_EXC
          REAL(DP)  :: ESELF
          REAL(DP)  :: ESR
          REAL(DP)  :: EVDW
          REAL(DP)  :: EBAND
          REAL(DP)  :: EKIN
          REAL(DP)  :: ATOT     ! Ensemble DFT
          REAL(DP)  :: ENTROPY  ! Ensemble DFT
          REAL(DP)  :: EGRAND   ! Ensemble DFT
          REAL(DP)  :: VAVE     ! Ensemble DFT
          REAL(DP)  :: EEXTFOR  ! Energy of the external forces
        END TYPE

        REAL(DP)  :: EHTE = 0.0_DP
        REAL(DP)  :: SELF_EHTE = 0.0_DP
        REAL(DP)  :: EHTI = 0.0_DP
        REAL(DP)  :: EH = 0.0_DP
        REAL(DP)  :: EHT = 0.0_DP
        REAL(DP)  :: SELF_EXC = 0.0_DP
        REAL(DP)  :: SELF_VXC = 0.0_DP
        REAL(DP)  :: EKIN = 0.0_DP
        REAL(DP)  :: ESELF = 0.0_DP
        REAL(DP)  :: EVDW = 0.0_DP
        REAL(DP)  :: EPSEU = 0.0_DP
        REAL(DP)  :: ENT = 0.0_DP
        REAL(DP)  :: ETOT = 0.0_DP
        REAL(DP)  :: ENL = 0.0_DP
        REAL(DP)  :: ESR = 0.0_DP
        REAL(DP)  :: EXC = 0.0_DP
        REAL(DP)  :: VXC = 0.0_DP
        REAL(DP)  :: exx = 0.0_DP                  
        REAL(DP)  :: EBAND = 0.0_DP
        REAL(DP)  :: ATOT = 0.0_DP
        REAL(DP)  :: ENTROPY = 0.0_DP
        REAL(DP)  :: EGRAND = 0.0_DP
        REAL(DP)  :: VAVE = 0.0_DP    ! average potential
        REAL(DP)  :: EEXTFOR = 0.0_DP  ! Energy of the external forces
        
        REAL(DP) :: enthal = 0.0_DP, ekincm

        PUBLIC :: dft_energy_type, total_energy, eig_total_energy, &
                  print_energies, debug_energies

        PUBLIC :: etot, eself, enl, ekin, epseu, esr, eht, exc, ekincm, exx  
        PUBLIC :: self_exc, self_ehte

        PUBLIC :: atot, entropy, egrand, enthal, vave

        PUBLIC :: eextfor

      CONTAINS

! ---------------------------------------------------------------------------- !

        SUBROUTINE total_energy( edft )

          TYPE (dft_energy_type) :: edft

          eself      = edft%eself
          epseu      = edft%epseu
          ent        = edft%ent
          enl        = edft%enl
          evdw       = edft%evdw
          esr        = edft%esr
          ekin       = edft%ekin
          vxc        = edft%vxc
          ehti       = edft%ehti
          ehte       = edft%ehte
          self_ehte  = edft%self_ehte
          self_exc   = edft%self_exc
          self_vxc   = edft%self_vxc
          exc        = edft%exc
          eht        = edft%eht

          etot  = ekin + eht + epseu + enl + exc + evdw - ent
          !
          edft%etot = etot

          RETURN
        END SUBROUTINE total_energy

! ---------------------------------------------------------------------------- !

        SUBROUTINE eig_total_energy(ei)
          IMPLICIT NONE
          REAL(DP), INTENT(IN) :: ei(:)
          INTEGER :: i
          REAL(DP) etot_band, EII
          eband = 0.0_DP
          do i = 1, SIZE(ei)
            eband = eband + ei(i) * 2.0_DP
          end do
          EII = ehti + ESR - ESELF
          etot_band = eband - ehte + (exc-vxc) + eii
          WRITE( stdout,200) etot_band, eband, ehte, (exc-vxc), eii
 200      FORMAT(' *** TOTAL ENERGY : ',F14.8,/ &
                ,'     eband        : ',F14.8,/ &
                ,'     eh           : ',F14.8,/ &
                ,'     xc           : ',F14.8,/ &
                ,'     eii          : ',F14.8)
          RETURN
        END SUBROUTINE eig_total_energy

! ---------------------------------------------------------------------------- !

        SUBROUTINE print_energies( tsic, iprsta, edft, sic_alpha, sic_epsilon, textfor )
          LOGICAL, INTENT(IN) :: tsic
          TYPE (dft_energy_type), OPTIONAL, INTENT(IN) :: edft
          INTEGER, OPTIONAL, INTENT(IN) :: iprsta
          REAL(DP), OPTIONAL, INTENT(IN) :: sic_alpha, sic_epsilon
          LOGICAL, OPTIONAL, INTENT(IN) :: textfor

          IF( PRESENT ( edft ) ) THEN
              WRITE( stdout,  * )
              WRITE( stdout,  * )
              WRITE( stdout,  1 ) edft%etot
              WRITE( stdout,  2 ) edft%ekin
              WRITE( stdout,  3 ) edft%eht
              WRITE( stdout,  4 ) edft%eself   !   self interaction of the pseudocharges NOT SIC!
              WRITE( stdout,  5 ) edft%esr
              WRITE( stdout,  9 ) edft%epseu
              WRITE( stdout, 10 ) edft%enl
              WRITE( stdout, 11 ) edft%exc
              IF( PRESENT( iprsta ) ) THEN
                 IF( iprsta > 1 ) THEN
                    WRITE( stdout,  * )
                    WRITE( stdout,  6 ) edft%eh
                    WRITE( stdout,  7 ) edft%ehte
                    WRITE( stdout,  8 ) edft%ehti
                    WRITE( stdout, 12 ) edft%evdw
                    WRITE( stdout, 13 ) edft%emkin
                 END IF
              END IF
          ELSE
             !
             WRITE( stdout,100) etot, ekin, eht, esr, eself, epseu, enl, exc, vave

!====================================================================================
!exx_wf related
             if(dft_is_hybrid()) then
                WRITE( stdout,101) -exx*get_exx_fraction(), etot
             end if
!====================================================================================

          END IF
          !
          IF( tsic ) THEN
             !
             IF( .NOT. PRESENT( sic_alpha ) .OR. .NOT. PRESENT( sic_epsilon ) ) &
                CALL errore( ' print_energies ', ' sic without parameters? ', 1 )

             WRITE( stdout, fmt = "('Sic contributes in Mauri&al. approach:')" )
             WRITE( stdout, fmt = "('--------------------------------------')" )
             !
             !  qui e' da aggiungere i due parametetri alpha_si e si_epsilon che determinano "quanto"
             !  correggo lo exc e hartree
             !           
             WRITE( stdout, 14 ) self_ehte, sic_epsilon
             WRITE( stdout, 15 ) self_exc, sic_alpha
          END IF
          !
          IF( PRESENT( textfor ) ) THEN
             IF( textfor ) WRITE( stdout, 16 ) eextfor
          END IF
          !
          CALL plugin_print_energies()
          !
1         FORMAT(6X,'                total energy = ',F18.10,' Hartree a.u.')
2         FORMAT(6X,'              kinetic energy = ',F18.10,' Hartree a.u.')
3         FORMAT(6X,'        electrostatic energy = ',F18.10,' Hartree a.u.')
4         FORMAT(6X,'                       eself = ',F18.10,' Hartree a.u.')
5         FORMAT(6X,'                         esr = ',F18.10,' Hartree a.u.')
6         FORMAT(6X,'              hartree energy = ',F18.10,' Hartree a.u.')
7         FORMAT(6X,'                hartree ehte = ',F18.10,' Hartree a.u.')
8         FORMAT(6X,'                hartree ehti = ',F18.10,' Hartree a.u.')
9         FORMAT(6X,'      pseudopotential energy = ',F18.10,' Hartree a.u.')
10        FORMAT(6X,'  n-l pseudopotential energy = ',F18.10,' Hartree a.u.')
11        FORMAT(6X,' exchange-correlation energy = ',F18.10,' Hartree a.u.')
12        FORMAT(6X,'        van der waals energy = ',F18.10,' Hartree a.u.')
13        FORMAT(6X,'        emass kinetic energy = ',F18.10,' Hartree a.u.')
14        FORMAT(6X,'            hartree sic_ehte = ',F18.10,' Hartree a.u.', 1X, 'corr. factor = ',F6.3)
15        FORMAT(6X,' sic exchange-correla energy = ',F18.10,' Hartree a.u.', 1X, 'corr. factor = ',F6.3)
16        FORMAT(6X,'       external force energy = ',F18.10,' Hartree a.u.')

  100 format(//'                total energy = ',f20.11,' Hartree a.u.'/ &
     &         '              kinetic energy = ',f14.5,' Hartree a.u.'/ &
     &         '        electrostatic energy = ',f14.5,' Hartree a.u.'/ &
     &         '                         esr = ',f14.5,' Hartree a.u.'/ &
     &         '                       eself = ',f14.5,' Hartree a.u.'/ &
     &         '      pseudopotential energy = ',f14.5,' Hartree a.u.'/ &
     &         '  n-l pseudopotential energy = ',f14.5,' Hartree a.u.'/ &
     &         ' exchange-correlation energy = ',f14.5,' Hartree a.u.'/ &
     &         '           average potential = ',f14.5,' Hartree a.u.'//)
  101 format(//'                  exx energy = ',F14.5,' Hartree a.u.'/ &
     &         '       total energy with exx = ',F14.5,' Hartree a.u.' / )
          RETURN
        END SUBROUTINE print_energies


! ---------------------------------------------------------------------------- !

        SUBROUTINE debug_energies( edft )
          TYPE (dft_energy_type), OPTIONAL, INTENT(IN) :: edft
          IF( PRESENT ( edft ) ) THEN
            WRITE( stdout,2) edft%ETOT, edft%EKIN, edft%EHT, &
              edft%ESELF, edft%ESR, edft%EH, &
              edft%EPSEU, edft%ENL, edft%EXC, edft%VXC, edft%EVDW, edft%EHTE, &
              edft%EHTI, edft%ENT, edft%EBAND, (edft%EXC-edft%VXC), &
              (edft%EHTI+edft%ESR-edft%ESELF), &
              edft%EBAND-edft%EHTE+(edft%EXC-edft%VXC)+(edft%EHTI+edft%ESR-edft%ESELF)
          ELSE
            WRITE( stdout,2) ETOT, EKIN, EHT, ESELF, ESR, EH, EPSEU, ENL, EXC, VXC, &
              EVDW, EHTE, EHTI, ENT, EBAND, (EXC-VXC), (EHTI+ESR-ESELF), &
              EBAND-EHTE+(EXC-VXC)+(EHTI+ESR-ESELF)
          END IF
2         FORMAT(/,/ &
            ,6X,' ETOT .... = ',F18.10,/ &
            ,6X,' EKIN .... = ',F18.10,/ &
            ,6X,' EHT ..... = ',F18.10,/ &
            ,6X,' ESELF ... = ',F18.10,/ &
            ,6X,' ESR ..... = ',F18.10,/ &
            ,6X,' EH ...... = ',F18.10,/ &
            ,6X,' EPSEU ... = ',F18.10,/ &
            ,6X,' ENL ..... = ',F18.10,/ &
            ,6X,' EXC ..... = ',F18.10,/ &
            ,6X,' VXC ..... = ',F18.10,/ &
            ,6X,' EVDW .... = ',F18.10,/ &
            ,6X,' EHTE .... = ',F18.10,/ &
            ,6X,' EHTI .... = ',F18.10,/ &
            ,6X,' ENT ..... = ',F18.10,/ &
            ,6X,' EBAND ... = ',F18.10,/ &
            ,6X,' EXC-VXC ............................. = ',F18.10,/ &
            ,6X,' EHTI+ESR-ESELF ...................... = ',F18.10,/ &
            ,6X,' EBAND-EHTE+(EXC-VXC)+(EHTI+ESR-ESELF) = ',F18.10)
          RETURN
        END SUBROUTINE debug_energies


      END MODULE energies
