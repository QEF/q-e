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
  ! NOTE: when a dft term slot is filled with 'xxxx' it usually means that the term
  ! is included in one of the full dft combinations available, but it cannot be
  ! called by itself.
  !
  ! NOTE: The references and descriptions of each dft are stored in a different
  ! module in file 'qe_dft_refs.f90', which is used by the program 'xc_infos'.
  ! If you modify something here, please check consistency there too.
  !
  USE dft_setting_params,  ONLY: notset
  !
  SAVE
  !
  INTEGER, PARAMETER :: nxc=10, ncc=14, ngcx=50, ngcc=14, nmeta=8
  CHARACTER(LEN=4)   :: dft_LDAx_name(0:nxc),  dft_LDAc_name(0:ncc),  &
                        dft_GGAx_name(0:ngcx), dft_GGAc_name(0:ngcc), &
                        dft_MGGA_name(0:nmeta)
  !
  TYPE dft_label
     CHARACTER(LEN=10) :: name
     CHARACTER(LEN=10) :: name2
     INTEGER :: IDs(6)
  END TYPE dft_label
  !
  INTEGER, PARAMETER :: n_dft=42
  TYPE(dft_label) :: dft_full(n_dft)
  !
  !
  !
  ! LDA exchange terms
  DATA dft_LDAx_name / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', & ! 0 to  6
                       'B3LP', 'KZK', 'xxxx', 'xxxx' /                    ! 7 "  10
  ! LDA correlation terms
  DATA dft_LDAc_name / 'NOC', 'PZ', 'VWN', 'LYP', 'PW',   'WIG', 'HL',  & ! 0 to  6
                       'OBZ', 'OBW', 'GL', 'KZK', 'xxxx', 'B3LP','xxxx',& ! 7 "  13
                       'xxxx' /                                           !14
  ! GGA exchange terms
  DATA dft_GGAx_name / 'NOGX', 'B88',  'GGX',  'PBX',  'REVX', 'HCTH',  & ! 0 to  5
                       'OPTX', 'xxxx', 'PB0X', 'B3LP', 'PSX',  'WCX',   & ! 6 "  11
                       'HSE',  'RW86', 'PBE',  'xxxx', 'C09X', 'SOX',   & !12 "  17 
                       'xxxx', 'Q2DX', 'GAUP', 'PW86', 'B86B', 'OBK8',  & !18 "  23
                       'OB86', 'EVX',  'B86R', 'CX13', 'X3LP', 'CX0',   & !24 "  29
                       'R860', 'CX0P', 'AHCX', 'AHF2', 'AHPB', 'AHPS',  & !30 "  35
                       'CX14', 'CX15', 'BR0',  'CX16', 'C090', 'B86X',  & !36 "  41
                       'B88X', 'BEEX', 'HHNX', 'W31X', 'W32X', 'AHBR',  & !42 "  47 
                       'EHPB', 'HJPB', 'HJPS' /                           !48 "  50 
  ! GGA correlation terms
  DATA dft_GGAc_name / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH',     & ! 0 to  5
                       'NONE', 'B3LP','PSC', 'PBE' , 'xxxx','xxxx',     & ! 6 "  11
                       'Q2DC', 'xxxx','BEEC' /                            !12 "  14
  ! MGGA exchange+correlation terms
  DATA dft_MGGA_name / 'NONE', 'TPSS', 'M06L', 'TB09', 'NONE', 'SCAN',  & ! 0 to  5
                       'SCA0', 'xxxx', 'xxxx' /                           ! 6 "   8
  !
  !
  ! ---- Full DFTs (except vdW-DFs, those are still under Modules) ----
  !
  DATA dft_full(1)%name      / 'PZ'  /
  DATA dft_full(1)%name2     / 'LDA' /
  DATA dft_full(1)%IDs(1:6)  / 1,1,0,0,0,0 /    ! sla+pz
  !
  DATA dft_full(2)%name      / 'PW'   /
  DATA dft_full(2)%name2     / 'none' /
  DATA dft_full(2)%IDs(1:6)  / 1,4,0,0,0,0 /    ! sla+pw
  !
  DATA dft_full(3)%name      / 'VWN-RPA' /
  DATA dft_full(3)%name2     / 'none'    /
  DATA dft_full(3)%IDs(1:6)  / 1,11,0,0,0,0 /   ! sla+xxxx[vwn1_rpa]
  !
  DATA dft_full(4)%name      / 'OEP'  /
  DATA dft_full(4)%name2     / 'none' /
  DATA dft_full(4)%IDs(1:6)  / 4,0,0,0,0,0 /    ! oep
  !
  DATA dft_full(5)%name      / 'KLI'  /
  DATA dft_full(5)%name2     / 'none' /
  DATA dft_full(5)%IDs(1:6)  / 10,0,0,0,0,0 /   ! kli
  !
  DATA dft_full(6)%name      / 'HF'   /
  DATA dft_full(6)%name2     / 'none' /
  DATA dft_full(6)%IDs(1:6)  / 5,0,0,0,0,0 /    ! hf
  !
  DATA dft_full(7)%name      / 'PBE'  /
  DATA dft_full(7)%name2     / 'none' /
  DATA dft_full(7)%IDs(1:6)  / 1,4,3,4,0,0 /    ! sla+pw+pbx+pbc
  !
  DATA dft_full(8)%name      / 'B88'  /
  DATA dft_full(8)%name2     / 'none' /
  DATA dft_full(8)%IDs(1:6)  / 1,1,1,0,0,0 /    ! sla+pz+b88
  !
  DATA dft_full(9)%name      / 'BP'   /
  DATA dft_full(9)%name2     / 'none' /
  DATA dft_full(9)%IDs(1:6)  / 1,1,1,1,0,0 /    ! sla+pz+b88+p86
  !----
  DATA dft_full(10)%name     / 'PW91' /
  DATA dft_full(10)%name2    / 'none' /
  DATA dft_full(10)%IDs(1:6) / 1,4,2,2,0,0 /    ! sla+pw+ggx+ggc
  !
  DATA dft_full(11)%name     / 'REVPBE' /
  DATA dft_full(11)%name2    / 'none'   /
  DATA dft_full(11)%IDs(1:6) / 1,4,4,4,0,0 /    ! sla+pw+revx+pbc
  !
  DATA dft_full(12)%name     / 'PBESOL' /
  DATA dft_full(12)%name2    / 'none'   /
  DATA dft_full(12)%IDs(1:6) / 1,4,10,8,0,0 /   ! sla+pw+psx+psc
  !
  DATA dft_full(13)%name     / 'BLYP' /
  DATA dft_full(13)%name2    / 'none' /
  DATA dft_full(13)%IDs(1:6) / 1,3,1,3,0,0 /    ! sla+lyp+b88+blyp
  !
  DATA dft_full(14)%name     / 'OPTBK88' /
  DATA dft_full(14)%name2    / 'none'    /
  DATA dft_full(14)%IDs(1:6) / 1,4,23,1,0,0 /   ! sla+pw+obk8+p86
  !
  DATA dft_full(15)%name     / 'OPTB86B' /
  DATA dft_full(15)%name2    / 'none'    /
  DATA dft_full(15)%IDs(1:6) / 1,4,24,1,0,0 /   ! sla+pw+ob86+p86
  !
  DATA dft_full(16)%name     / 'PBC'  /
  DATA dft_full(16)%name2    / 'none' /
  DATA dft_full(16)%IDs(1:6) / 1,4,0,4,0,0 /    ! sla+pw+pbc
  !
  DATA dft_full(17)%name     / 'HCTH' /
  DATA dft_full(17)%name2    / 'none' /
  DATA dft_full(17)%IDs(1:6) / 0,0,5,5,0,0 /    ! nox+noc+hcth+hcth
  !
  DATA dft_full(18)%name     / 'OLYP' /
  DATA dft_full(18)%name2    / 'none' /
  DATA dft_full(18)%IDs(1:6) / 0,3,6,3,0,0 /    ! nox+lyp+optx+blyp
  !
  DATA dft_full(19)%name     / 'WC'   /
  DATA dft_full(19)%name2    / 'none' /
  DATA dft_full(19)%IDs(1:6) / 1,4,11,4,0,0 /   ! sla+pw+wcx+pbc
  !
  DATA dft_full(20)%name     / 'PW86PBE' /
  DATA dft_full(20)%name2    / 'none'    /
  DATA dft_full(20)%IDs(1:6) / 1,4,21,4,0,0 /   ! sla+pw+pw86+pbc
  !
  DATA dft_full(21)%name     / 'B86BPBE' /
  DATA dft_full(21)%name2    / 'none'    /
  DATA dft_full(21)%IDs(1:6) / 1,4,22,4,0,0 /   ! sla+pw+b86b+pbc
  !
  DATA dft_full(22)%name     / 'PBEQ2D' /
  DATA dft_full(22)%name2    / 'Q2D'    /
  DATA dft_full(22)%IDs(1:6) / 1,4,19,12,0,0 /  ! sla+pw+q2dx+q2dc
  !
  DATA dft_full(23)%name     / 'SOGGA' /
  DATA dft_full(23)%name2    / 'none'  /
  DATA dft_full(23)%IDs(1:6) / 1,4,17,4,0,0 /   ! sla+pw+sox+pbec
  !
  DATA dft_full(24)%name     / 'EV93' /
  DATA dft_full(24)%name2    / 'none' /
  DATA dft_full(24)%IDs(1:6) / 1,4,25,0,0,0 /   ! sla+pw+evx+nogc
  !
  DATA dft_full(25)%name     / 'RPBE' /
  DATA dft_full(25)%name2    / 'none' /
  DATA dft_full(25)%IDs(1:6) / 1,4,44,4,0,0 /   ! sla+pw+hhnx+pbc
  !
  DATA dft_full(26)%name     / 'PBE0' /
  DATA dft_full(26)%name2    / 'none' /
  DATA dft_full(26)%IDs(1:6) / 6,4,8,4,0,0 /    ! pb0x+pw+pb0x+pbc
  !
  DATA dft_full(27)%name     / 'B86BPBEX' /
  DATA dft_full(27)%name2    / 'none'     /
  DATA dft_full(27)%IDs(1:6) / 6,4,41,4,0,0 /   ! sla+pw+b86x+pbc
  !
  DATA dft_full(28)%name     / 'BHAHLYP'   /
  DATA dft_full(28)%name2    / 'BHANDHLYP' /
  DATA dft_full(28)%IDs(1:6) / 6,4,42,3,0,0 /   ! pb0x+pw+b88x+blyp
  !
  DATA dft_full(29)%name     / 'HSE'  /
  DATA dft_full(29)%name2    / 'none' /
  DATA dft_full(29)%IDs(1:6) / 1,4,12,4,0,0 /   ! sla+pw+hse+pbc
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
  DATA dft_full(30)%name     / 'GAUP'   /
  DATA dft_full(30)%name2    / 'GAUPBE' /
  DATA dft_full(30)%IDs(1:6) / 1,4,20,4,0,0 /   ! sla+pw+gaup+pbc
  !
  DATA dft_full(31)%name     / 'B3LYP' /
  DATA dft_full(31)%name2    / 'none'  /
  DATA dft_full(31)%IDs(1:6) / 7,12,9,7,0,0 /   ! b3lp+b3lp+b3lp+b3lp
  !
  DATA dft_full(32)%name     / 'B3LYP-V1R' /
  DATA dft_full(32)%name2    / 'none'      /
  DATA dft_full(32)%IDs(1:6) / 7,13,9,7,0,0 /   ! b3lp+xxxx[b3lyp_v1r]+b3lp+b3lp
  !
  DATA dft_full(33)%name     / 'X3LYP' /
  DATA dft_full(33)%name2    / 'none'  /
  DATA dft_full(33)%IDs(1:6) / 9,14,28,13,0,0 / ! xxxx[x3lyp_ldax]+xxxx[x3lyp_ldac]+x3lp+x3lc
  !
  DATA dft_full(34)%name     / 'TPSS' /
  DATA dft_full(34)%name2    / 'none'      /
  DATA dft_full(34)%IDs(1:6) / 0,0,0,0,1,0 /
  !
  DATA dft_full(35)%name     / 'M06L' /
  DATA dft_full(35)%name2    / 'none' /
  DATA dft_full(35)%IDs(1:6) / 0,0,0,0,2,0 /
  !
  DATA dft_full(36)%name     / 'TB09' /
  DATA dft_full(36)%name2    / 'none' /
  DATA dft_full(36)%IDs(1:6) / 0,0,0,0,3,0 /
  !
  DATA dft_full(37)%name     / 'SCAN' /
  DATA dft_full(37)%name2    / 'none' /
  DATA dft_full(37)%IDs(1:6) / 0,0,0,0,5,0 /    ! scan[wrapper to Libxc SCAN]
  !
  DATA dft_full(38)%name     / 'SCAN0' /
  DATA dft_full(38)%name2    / 'none'  /
  DATA dft_full(38)%IDs(1:6) / 0,0,0,0,6,0 /
  !
  DATA dft_full(39)%name     / 'R2SCAN' /
  DATA dft_full(39)%name2    / 'none' /
  DATA dft_full(39)%IDs(1:6) / 0,0,0,0,7,0 /
  ! 'AH series for vdW-DFs and analytical-hole PBE(sol)-AH: JPCM 34, 025902 (2022)
  DATA dft_full(40)%name     / 'PBE-AH' /
  DATA dft_full(40)%name2    / 'none' /
  DATA dft_full(40)%IDs(1:6) / 1,4,34,4,0,0 /   ! P.H.: sla+pw+ahpb+pbc; Differs from EHPB-based, incl HSE
  !
  DATA dft_full(41)%name     / 'PBESOL-AH' /
  DATA dft_full(41)%name2    / 'none' /
  DATA dft_full(41)%IDs(1:6) / 1,4,35,8,0,0 /   ! P.H.: sla+pw+ahps+psc
  !
  DATA dft_full(42)%name     / 'RSCAN' /
  DATA dft_full(42)%name2    / 'none' /
  DATA dft_full(42)%IDs(1:6) / 0,0,0,0,8,0 /
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
    DO i = 1, n_dft
      IF (name==dft_full(i)%name .OR. name==dft_full(i)%name2) THEN
        IDs(:) = dft_full(i)%IDs(:)
        EXIT
      ENDIF  
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------
  SUBROUTINE get_shortname_from_IDs( IDs, name, id_full )
    !---------------------------------------------------------------
    !! Get the DFT shortname from ID numbers of each family-kind term.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: IDs(6)
    CHARACTER(LEN=*), INTENT(INOUT) :: name
    INTEGER, INTENT(OUT), OPTIONAL :: id_full
    INTEGER :: i
    !
    DO i = 1, n_dft
      IF (ALL( IDs(:)==dft_full(i)%IDs(:) )) THEN
        name = dft_full(i)%name
        IF (PRESENT(id_full)) id_full = i
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

