!
! Copyright (C) 2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------
MODULE qe_dft_list
  !----------------------------------------------------------------------------------
  !! Contains the list of functionals available in QE, both single terms (family+type)
  !! and combinations.
  !! Extra parameters of QE functionals are set in 'xclib_set_auxiliary_flag' and 
  !! subsequent routines 
  !
  ! NOTE: when a dft term slot is filled with 'xxxx' it means that the term
  ! is included in one of the full dft combinations available, but it cannot
  ! be used by itself.
  !
  USE dft_setting_params,  ONLY: notset
  
  SAVE
  !
  ! -- single DFT terms (family-type)
  INTEGER, PARAMETER :: nxc=10, ncc=14, ngcx=46, ngcc=14, nmeta=6
  CHARACTER(LEN=4)   :: dft_LDAx_name(0:nxc),  dft_LDAc_name(0:ncc),  &
                        dft_GGAx_name(0:ngcx), dft_GGAc_name(0:ngcc), &
                        dft_MGGA_name(0:nmeta)
  !
  ! -- total DFTs
  INTEGER, PARAMETER :: n_dft=41
  CHARACTER(LEN=10)  :: dft_name(n_dft)
  CHARACTER(LEN=10)  :: dft_name2(n_dft)
  INTEGER, DIMENSION(n_dft,6) :: dft_IDs
  !
  !
  ! - LDA exchange
  DATA dft_LDAx_name / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', &
                       'B3LP', 'KZK', 'xxxx', 'xxxx' /
  ! - LDA correlation
  DATA dft_LDAc_name / 'NOC', 'PZ', 'VWN', 'LYP', 'PW',   'WIG', 'HL',  &
                       'OBZ', 'OBW', 'GL', 'KZK', 'xxxx', 'B3LP','xxxx',&
                       'xxxx' /
  ! - GGA exchange
  DATA dft_GGAx_name / 'NOGX', 'B88',  'GGX',  'PBX',  'REVX', 'HCTH',  &
                       'OPTX', 'xxxx', 'PB0X', 'B3LP', 'PSX',  'WCX',   &
                       'HSE',  'RW86', 'PBE',  'xxxx', 'C09X', 'SOX',   &
                       'xxxx', 'Q2DX', 'GAUP', 'PW86', 'B86B', 'OBK8',  &
                       'OB86', 'EVX',  'B86R', 'CX13', 'X3LP', 'CX0',   &
                       'R860', 'CX0P', 'AHCX', 'AHF2', 'AHPB', 'AHPS',  &
                       'CX14', 'CX15', 'BR0',  'CX16', 'C090', 'B86X',  &
                       'B88X', 'BEEX', 'HHNX', 'W31X', 'W32X' /
  ! - GGA correlation
  DATA dft_GGAc_name / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH',     &
                       'NONE', 'B3LP','PSC', 'PBE',  'xxxx','xxxx',     &
                       'Q2DC', 'xxxx','BEEC' /
  ! - MGGA exchange
  DATA dft_MGGA_name / 'NONE', 'TPSS', 'M06L', 'TB09', 'META', 'SCAN',  &
                       'SCA0' /
  !
  !
  ! ---- Full DFTs ----
  !
  DATA dft_name(1)     / 'PZ' /
  DATA dft_name2(1)    / 'LDA' /
  DATA dft_IDs(1,1:6)  / 1,1,0,0,0,0 /    ! sla+pz
  !
  DATA dft_name(2)     / 'PW' /
  DATA dft_name2(2)    / 'none' /
  DATA dft_IDs(2,1:6)  / 1,4,0,0,0,0 /    ! sla+pw
  !
  DATA dft_name(3)     / 'VWN-RPA' /
  DATA dft_name2(3)    / 'none' /
  DATA dft_IDs(3,1:6)  / 1,11,0,0,0,0 /   ! sla+xxxx[vwn1_rpa]
  !
  DATA dft_name(4)     / 'OEP' /
  DATA dft_name2(4)    / 'none' /
  DATA dft_IDs(4,1:6)  / 4,0,0,0,0,0 /    ! oep
  !
  DATA dft_name(5)     / 'KLI' /
  DATA dft_name2(5)    / 'none' /
  DATA dft_IDs(5,1:6)  / 10,0,0,0,0,0 /   ! kli
  !
  DATA dft_name(6)     / 'HF' /
  DATA dft_name2(6)    / 'none' /
  DATA dft_IDs(6,1:6)  / 5,0,0,0,0,0 /    ! hf
  !
  DATA dft_name(7)     / 'PBE' /
  DATA dft_name2(7)    / 'none' /
  DATA dft_IDs(7,1:6)  / 1,4,3,4,0,0 /    ! sla+pw+pbx+pbc
  !
  DATA dft_name(8)     / 'B88' /
  DATA dft_name2(8)    / 'none' /
  DATA dft_IDs(8,1:6)  / 1,1,1,0,0,0 /    ! sla+pz+b88
  !
  DATA dft_name(9)     / 'BP' /
  DATA dft_name2(9)    / 'none' /
  DATA dft_IDs(9,1:6)  / 1,1,1,1,0,0 /    ! sla+pz+b88+p86
  !----
  DATA dft_name(10)    / 'PW91' /
  DATA dft_name2(10)   / 'none' /
  DATA dft_IDs(10,1:6) / 1,4,2,2,0,0 /    ! sla+pw+ggx+ggc
  !
  DATA dft_name(11)    / 'REVPBE' /
  DATA dft_name2(11)   / 'none' /
  DATA dft_IDs(11,1:6) / 1,4,4,4,0,0 /    ! sla+pw+revx+pbc
  !
  DATA dft_name(12)    / 'PBESOL' /
  DATA dft_name2(12)   / 'none' /
  DATA dft_IDs(12,1:6) / 1,4,10,8,0,0 /   ! sla+pw+psx+psc
  !
  DATA dft_name(13)    / 'BLYP' /
  DATA dft_name2(13)   / 'none' /
  DATA dft_IDs(13,1:6) / 1,3,1,3,0,0 /    ! sla+lyp+b88+blyp
  !
  DATA dft_name(14)    / 'OPTBK88' /
  DATA dft_name2(14)   / 'none' /
  DATA dft_IDs(14,1:6) / 1,4,23,1,0,0 /   ! sla+pw+obk8+p86
  !
  DATA dft_name(15)    / 'OPTB86B' /
  DATA dft_name2(15)   / 'none' /
  DATA dft_IDs(15,1:6) / 1,4,24,1,0,0 /   ! sla+pw+ob86+p86
  !
  DATA dft_name(16)    / 'PBC' /
  DATA dft_name2(16)   / 'none' /
  DATA dft_IDs(16,1:6) / 1,4,0,4,0,0 /    ! sla+pw+pbc
  !
  DATA dft_name(17)    / 'HCTH' /
  DATA dft_name2(17)   / 'none' /
  DATA dft_IDs(17,1:6) / 0,0,5,5,0,0 /    ! nox+noc+hcth+hcth
  !
  DATA dft_name(18)    / 'OLYP' /
  DATA dft_name2(18)   / 'none' /
  DATA dft_IDs(18,1:6) / 0,3,6,3,0,0 /    ! nox+lyp+optx+blyp
  !
  DATA dft_name(19)    / 'WC' /
  DATA dft_name2(19)   / 'none' /
  DATA dft_IDs(19,1:6) / 1,4,11,4,0,0 /   ! sla+pw+wcx+pbc
  !
  DATA dft_name(20)    / 'PW86PBE' /
  DATA dft_name2(20)   / 'none' /
  DATA dft_IDs(20,1:6) / 1,4,21,4,0,0 /   ! sla+pw+pw86+pbc
  !
  DATA dft_name(21)    / 'B86BPBE' /
  DATA dft_name2(21)   / 'none' /
  DATA dft_IDs(21,1:6) / 1,4,22,4,0,0 /   ! sla+pw+b86b+pbc
  !
  DATA dft_name(22)    / 'PBEQ2D' /
  DATA dft_name2(22)   / 'Q2D' /
  DATA dft_IDs(22,1:6) / 1,4,19,12,0,0 /  ! sla+pw+q2dx+q2dc
  !
  DATA dft_name(23)    / 'SOGGA' /
  DATA dft_name2(23)   / 'none' /
  DATA dft_IDs(23,1:6) / 1,4,17,4,0,0 /   ! sla+pw+sox+pbec
  !
  DATA dft_name(24)    / 'EV93' /
  DATA dft_name2(24)   / 'none' /
  DATA dft_IDs(24,1:6) / 1,4,25,0,0,0 /   ! sla+pw+evx+nogc
  !
  DATA dft_name(25)    / 'RPBE' /
  DATA dft_name2(25)   / 'none' /
  DATA dft_IDs(25,1:6) / 1,4,44,4,0,0 /   ! sla+pw+hhnx+pbc
  !
  DATA dft_name(26)    / 'PBE0' /
  DATA dft_name2(26)   / 'none' /
  DATA dft_IDs(26,1:6) / 6,4,8,4,0,0 /    ! pb0x+pw+pb0x+pbc
  !
  DATA dft_name(27)    / 'B86BPBEX' /
  DATA dft_name2(27)   / 'none' /
  DATA dft_IDs(27,1:6) / 6,4,41,4,0,0 /   ! sla+pw+b86x+pbc
  !
  DATA dft_name(28)    / 'BHAHLYP' /
  DATA dft_name2(28)   / 'BHANDHLYP' /
  DATA dft_IDs(28,1:6) / 6,4,42,3,0,0 /   ! pb0x+pw+b88x+blyp
  !
  DATA dft_name(29)    / 'HSE' /
  DATA dft_name2(29)   / 'none' /
  DATA dft_IDs(29,1:6) / 1,4,12,4,0,0 /   ! sla+pw+hse+pbc
  ! NOTE ABOUT HSE: there are two slight deviations with respect to the HSE06
  ! functional as it is in Gaussian code (that is considered as the reference
  ! in the chemistry community):
  ! - The range separation in Gaussian is precisely 0.11 bohr^-1,
  !   instead of 0.106 bohr^-1 in this implementation
  ! - The gradient scaling relation is a bit more complicated
  !   [ see: TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
  !          J. Chem. Phys. 131, 044108 (2009) ]
  ! These two modifications accounts only for a 1e-5 Ha difference for a
  ! single He atom. Info by Fabien Bruneval.
  !----
  DATA dft_name(30)    / 'GAUP' /
  DATA dft_name2(30)   / 'GAUPBE' /
  DATA dft_IDs(30,1:6) / 1,4,20,4,0,0 /   ! sla+pw+gaup+pbc
  !
  DATA dft_name(31)    / 'B3LYP' /
  DATA dft_name2(31)   / 'none' /
  DATA dft_IDs(31,1:6) / 7,12,9,7,0,0 /   ! b3lp+b3lp+b3lp+b3lp
  !
  DATA dft_name(32)    / 'B3LYP-V1R' /
  DATA dft_name2(32)   / 'none' /
  DATA dft_IDs(32,1:6) / 7,13,9,7,0,0 /   ! b3lp+xxxx[b3lyp_v1r]+b3lp+b3lp
  !
  DATA dft_name(33)    / 'X3LYP' /
  DATA dft_name2(33)   / 'none' /
  DATA dft_IDs(33,1:6) / 9,14,28,13,0,0 / ! xxxx[x3lyp_ldax]+xxxx[x3lyp_ldac]+x3lp+x3lc
  !
  DATA dft_name(34)    / 'TPSS' /
  DATA dft_name2(34)   / 'none' /
  DATA dft_IDs(34,1:6) / 1,4,7,6,1,0 /
  !
  DATA dft_name(35)    / 'TPSS-only' /
  DATA dft_name2(35)   / 'none' /
  DATA dft_IDs(35,1:6) / 0,0,0,0,1,0 /
  !
  DATA dft_name(36)    / 'M06L' /
  DATA dft_name2(36)   / 'none' /
  DATA dft_IDs(36,1:6) / 0,0,0,0,2,0 /
  !
  DATA dft_name(37)    / 'TB09' /
  DATA dft_name2(37)   / 'none' /
  DATA dft_IDs(37,1:6) / 0,0,0,0,3,0 /
  !
  DATA dft_name(38)    / 'SCAN' /
  DATA dft_name2(38)   / 'none' /
  DATA dft_IDs(38,1:6) / 0,0,0,0,5,0 /    ! scan[calls Libxc SCAN]
  !
  DATA dft_name(39)    / 'SCAN0' /
  DATA dft_name2(39)   / 'none' /
  DATA dft_IDs(39,1:6) / 0,0,0,0,6,0 /
  !
  DATA dft_name(40)    / 'PZ+META' /
  DATA dft_name2(40)   / 'LDA+META' /
  DATA dft_IDs(40,1:6) / 1,1,0,0,4,0 /
  !
  ! +meta activates MGGA even without MGGA-XC
  DATA dft_name(41)    / 'PBE+META' /
  DATA dft_name2(41)   / 'none' /
  DATA dft_IDs(41,1:6) / 1,4,3,4,4,0 /
  !
  !
CONTAINS
  !
  !------------------------------------------------------------------
  SUBROUTINE get_IDs_from_shortname( name, IDs )
    !---------------------------------------------------------------
    !! Get ID numbers of each family-kind term from the DFT shortname.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(INOUT) :: IDs(6)
    INTEGER :: i
    !
    IDs = notset
    !
    DO i = 1, n_dft
      IF (name==dft_name(i) .OR. name==dft_name2(i)) THEN
        IDs(:) = dft_IDs(i,:)
        EXIT
      ENDIF  
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------
  SUBROUTINE get_shortname_from_IDs( IDs, name )
    !---------------------------------------------------------------
    !! Get the DFT shortname from ID numbers of each family-kind term.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: IDs(6)
    CHARACTER(LEN=*), INTENT(INOUT) :: name
    INTEGER :: i
    !
    DO i = 1, n_dft
      IF (ALL(IDs(:)==dft_IDs(i,:))) THEN
        name = dft_name(i)
        EXIT
      ENDIF  
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE
  !
  !
END MODULE

