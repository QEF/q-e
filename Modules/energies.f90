!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE energies
        
        USE io_global,  ONLY : stdout
        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

        TYPE dft_energy_type
          REAL(dbl)  :: ETOT
          REAL(dbl)  :: SKIN
          REAL(dbl)  :: EMKIN
          REAL(dbl)  :: SELF_SXC
          REAL(dbl)  :: SXC
          REAL(dbl)  :: EHT
          REAL(dbl)  :: EH
          REAL(dbl)  :: SELF_EHTE
          REAL(dbl)  :: EHTE
          REAL(dbl)  :: EHTI
          REAL(dbl)  :: EPSEU
          REAL(dbl)  :: ENL
          REAL(dbl)  :: ENT
          REAL(dbl)  :: VXC
          REAL(dbl)  :: EXC
          REAL(dbl)  :: ESELF
          REAL(dbl)  :: ESR
          REAL(dbl)  :: EVDW
          REAL(dbl)  :: EBAND
          REAL(dbl)  :: EKIN
          REAL(dbl)  :: ATOT     ! Ensamble DFT
          REAL(dbl)  :: ENTROPY  ! Ensamble DFT
          REAL(dbl)  :: EGRAND   ! Ensamble DFT
          REAL(dbl)  :: VAVE   ! Ensamble DFT
        END TYPE

        REAL(dbl)  :: EHTE = 0.0_dbl
        REAL(dbl)  :: SELF_EHTE = 0.0_dbl
        REAL(dbl)  :: EHTI = 0.0_dbl
        REAL(dbl)  :: EH = 0.0_dbl
        REAL(dbl)  :: EHT = 0.0_dbl
        REAL(dbl)  :: SXC = 0.0_dbl
        REAL(dbl)  :: SELF_SXC = 0.0_dbl
        REAL(dbl)  :: EKIN = 0.0_dbl
        REAL(dbl)  :: ESELF = 0.0_dbl
        REAL(dbl)  :: EVDW = 0.0_dbl
        REAL(dbl)  :: EPSEU = 0.0_dbl
        REAL(dbl)  :: ENT = 0.0_dbl
        REAL(dbl)  :: ETOT = 0.0_dbl
        REAL(dbl)  :: ENL = 0.0_dbl
        REAL(dbl)  :: ESR = 0.0_dbl
        REAL(dbl)  :: EXC = 0.0_dbl
        REAL(dbl)  :: VXC = 0.0_dbl
        REAL(dbl)  :: SELF_VXC = 0.0_dbl
        REAL(dbl)  :: EBAND = 0.0_dbl
        REAL(dbl)  :: ATOT = 0.0_dbl
        REAL(dbl)  :: ENTROPY = 0.0_dbl
        REAL(dbl)  :: EGRAND = 0.0_dbl
        REAL(dbl)  :: VAVE = 0.0_dbl    ! average potential
        
        REAL(KIND=dbl) :: enthal, ekincm

        PUBLIC :: dft_energy_type, total_energy, eig_total_energy, &
                  print_energies, debug_energies

        PUBLIC :: etot, eself, enl, ekin, epseu, esr, eht, exc, ekincm

        PUBLIC :: atot, entropy, egrand, enthal, vave

      CONTAINS

! ---------------------------------------------------------------------------- !

!        SUBROUTINE total_energy( edft, omega, eexc, vvxc, eh, eps, &
!          self_ehte_in, self_sxc_in, self_vxc_in, nnr)

        SUBROUTINE total_energy( edft, omega, vvxc, eps, self_vxc_in, nnr)

          TYPE (dft_energy_type) :: edft
          REAL(dbl), INTENT(IN) :: OMEGA, VVXC
          REAL(dbl) :: VXC
          REAL(dbl) :: self_vxc_in
          COMPLEX(dbl), INTENT(IN) :: EPS
          INTEGER, INTENT(IN) :: nnr 

          eself = edft%eself
          ent   = edft%ent
          enl   = edft%enl
          evdw  = edft%evdw
          esr   = edft%esr
          ekin  = edft%ekin
          sxc   = edft%sxc
          ehti  = edft%ehte
          ehte  = edft%ehte
          self_ehte  = edft%ehte
          self_sxc   = edft%self_sxc

          self_vxc = self_vxc_in

          EXC   = edft%sxc * omega / DBLE(NNR) !EEXC * omega / DBLE(NNR)
          VXC   = VVXC * omega / DBLE(NNR)

          edft%exc  = exc
          edft%vxc  = vxc

          !EHT   = REAL( eh ) + esr - eself
          edft%eht =  edft%eh + esr - eself ! = eht
          EHT = edft%eht

          EPSEU = DBLE(eps)
          edft%epseu = epseu

          ETOT  = EKIN + EHT + EPSEU + ENL + EXC + EVDW - ENT
          edft%etot = etot

          RETURN
        END SUBROUTINE total_energy

! ---------------------------------------------------------------------------- !

        SUBROUTINE eig_total_energy(ei)
          IMPLICIT NONE
          REAL(dbl), INTENT(IN) :: ei(:)
          INTEGER :: i
          REAL(dbl) etot_band, EII
          eband = 0.0d0
          do i = 1, SIZE(ei)
            eband = eband + ei(i) * 2.0d0 
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

        SUBROUTINE print_energies( tsic, iprsta, edft )
          LOGICAL, INTENT(IN) :: tsic
          TYPE (dft_energy_type), OPTIONAL, INTENT(IN) :: edft
          INTEGER, OPTIONAL, INTENT(IN) :: iprsta
          REAL ( dbl ) :: EHT

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
              IF( tsic ) THEN
                WRITE( stdout, fmt = "('Sic contributes:')" )
                WRITE( stdout, fmt = "('----------------')" )
                WRITE( stdout, 14 ) edft%self_ehte
                WRITE( stdout, 15 ) edft%self_sxc
                WRITE( stdout, 16 ) vxc
                WRITE( stdout, 17 ) self_vxc
              END IF
            ! ETOT  = EKIN + EHT + EPSEU + ENL + EXC + EVDW - ENT
          ELSE
             !
 999         WRITE( stdout,100) etot,ekin,eht,esr,eself,epseu,enl,exc,vave
             !
          END IF
1         FORMAT(    6X,'                TOTAL ENERGY = ',F18.10,' A.U.')
2         FORMAT(    6X,'              KINETIC ENERGY = ',F18.10,' A.U.')
3         FORMAT(    6X,'        ELECTROSTATIC ENERGY = ',F18.10,' A.U.')
4         FORMAT(    6X,'                       ESELF = ',F18.10,' A.U.')
5         FORMAT(    6X,'                         ESR = ',F18.10,' A.U.')
6         FORMAT(    6X,'              HARTREE ENERGY = ',F18.10,' A.U.')
7         FORMAT(    6X,'                HARTREE EHTE = ',F18.10,' A.U.')
8         FORMAT(    6X,'                HARTREE EHTI = ',F18.10,' A.U.')
9         FORMAT(    6X,'      PSEUDOPOTENTIAL ENERGY = ',F18.10,' A.U.')
10        FORMAT(    6X,'  N-L PSEUDOPOTENTIAL ENERGY = ',F18.10,' A.U.')
11        FORMAT(    6X,' EXCHANGE-CORRELATION ENERGY = ',F18.10,' A.U.')
12        FORMAT(    6X,'        VAN DER WAALS ENERGY = ',F18.10,' A.U.')
13        FORMAT(    6X,'        EMASS KINETIC ENERGY = ',F18.10,' A.U.')
14        FORMAT(    6X,'            HARTREE SIC_EHTE = ',F18.10,' A.U.')
15        FORMAT(    6X,' SIC EXCHANGE-CORRELA ENERGY = ',F18.10,' A.U.')
16        FORMAT(    6X,'     EXCHANGE-CORRELA POTENT = ',F18.10,' A.U.')
17        FORMAT(    6X,' SIC EXCHANGE-CORRELA POTENT = ',F18.10,' A.U.')

  100 format(//'                total energy = ',f14.5,' a.u.'/         &
     &         '              kinetic energy = ',f14.5,' a.u.'/         &
     &         '        electrostatic energy = ',f14.5,' a.u.'/         &
     &         '                         esr = ',f14.5,' a.u.'/         &
     &         '                       eself = ',f14.5,' a.u.'/         &
     &         '      pseudopotential energy = ',f14.5,' a.u.'/         &
     &         '  n-l pseudopotential energy = ',f14.5,' a.u.'/         &
     &         ' exchange-correlation energy = ',f14.5,' a.u.'/         &
     &         '           average potential = ',f14.5,' a.u.'//)

          RETURN
        END SUBROUTINE print_energies


! ---------------------------------------------------------------------------- !

        SUBROUTINE debug_energies( edft )
          TYPE (dft_energy_type), OPTIONAL, INTENT(IN) :: edft
          IF( PRESENT ( edft ) ) THEN
            WRITE( stdout,2) edft%ETOT, edft%EKIN, edft%EHT, &
              edft%ESELF, edft%ESR, edft%EH, &
              edft%EPSEU, edft%ENL, edft%SXC, edft%VXC, edft%EVDW, edft%EHTE, &
              edft%EHTI, edft%ENT, edft%EBAND, (edft%SXC-edft%VXC), &
              (edft%EHTI+edft%ESR-edft%ESELF), &
              edft%EBAND-edft%EHTE+(edft%SXC-edft%VXC)+(edft%EHTI+edft%ESR-edft%ESELF)
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
