!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


   MODULE update

     IMPLICIT NONE
     SAVE

     PRIVATE

     PUBLIC :: updatevariables

   CONTAINS


      subroutine updatevariables(t_diis, TNOSEE, atoms_m, atoms_0, atoms_p, &
         cm, c0, cp, cdesc, htm2, htm, ht0, htp, tfor, thdyn, TNOSEP, edft)

      USE kinds
      USE energies, ONLY: dft_energy_type
      USE cell_module, ONLY: updatecell 
      USE nose_ions, ONLY: update_nose_ions 
      USE nose_electrons, ONLY: update_nose_electrons 
      USE wave_functions, ONLY: update_wave_functions 
      USE ions_module, ONLY: update_ions 
      USE cell_module, ONLY: boxdimensions
      USE atoms_type_module, ONLY: atoms_type
      USE wave_types, ONLY: wave_descriptor

      IMPLICIT NONE 
!
! ... arguments                                                         
      TYPE (atoms_type) :: atoms_m, atoms_0, atoms_p
      COMPLEX(dbl), INTENT(OUT)   :: cm(:,:,:,:)
      COMPLEX(dbl), INTENT(INOUT) :: c0(:,:,:,:)
      COMPLEX(dbl), INTENT(IN)    :: cp(:,:,:,:) 
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (boxdimensions) :: htm2, htm, ht0, htp 
      TYPE (dft_energy_type) :: edft
                                                                        
      LOGICAL, INTENT(IN) ::  t_diis 
      LOGICAL, INTENT(IN) ::  tnosee, tnosep 
      LOGICAL, INTENT(IN) ::  tfor, thdyn 

      INTEGER :: is
                                                                        
! ... subroutine body                                                   
      if(.not. t_diis) then 
        call update_wave_functions(cm, c0, cp, cdesc) 
        IF(TNOSEE) THEN 
          call update_nose_electrons 
        END IF 
      else
        IF(.not.tfor) THEN
          cm = c0
        END IF
      end if 
      IF(tfor) THEN 
        call update_ions(atoms_m, atoms_0, atoms_p) 
      END IF 
      IF(thdyn) THEN 
        CALL updatecell(htm2, htm, ht0, htp) 
      END IF 
      IF(tfor) THEN 
        IF(TNOSEP) THEN 
          call update_nose_ions 
        END IF 
      END IF 

      return 
      END  SUBROUTINE                                   

   END MODULE
