!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE cell_module
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl
      USE cell_base, ONLY: boxdimensions, alat, celldm, a1, a2, a3
      USE cell_base, ONLY: s_to_r, cell_init, pbcs, gethinv
      USE cell_base, ONLY: r_to_s, pbc, get_cell_param, dgcell, updatecell
      USE cell_base, ONLY: wc => wmass, press
      USE cell_base, ONLY: at, ainv
      USE cell_base, ONLY: ibrav, tcell_base_init
!
      IMPLICIT NONE
      SAVE
!
      PRIVATE

!
      PUBLIC :: gethinv, boxdimensions, pbc, get_cell_param, &
         updatecell, dgcell, movecell, r_to_s, s_to_r,  &
         get_lattice_vectors, pbcs, get_celldm, &
         cell_init, alat, press, at


!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!


        SUBROUTINE get_lattice_vectors(a1_out,a2_out,a3_out)
          REAL(dbl), intent(out) :: a1_out(3), a2_out(3), a3_out(3)
            a1_out   = a1
            a2_out   = a2
            a3_out   = a3
          RETURN
        END SUBROUTINE get_lattice_vectors

!------------------------------------------------------------------------------!

        SUBROUTINE get_celldm( ibrav_out, celldm_out)
          REAL(dbl), intent(out) :: celldm_out(6)
          INTEGER, intent(out) :: ibrav_out
            ibrav_out  = ibrav
            celldm_out = celldm
          RETURN
        END SUBROUTINE get_celldm

!------------------------------------------------------------------------------!

        SUBROUTINE movecell( tsdc, box_tm1, box_t0, box_tp1, velh )

          USE time_step, ONLY: delt
          USE cell_base, ONLY: cell_verlet, cell_steepest, iforceh, cell_move, &
                               frich
          USE cell_nose, ONLY: vnhh
          USE control_flags, ONLY: tnoseh

          LOGICAL :: tsdc
          TYPE (boxdimensions) :: box_tm1, box_t0, box_tp1
          REAL(dbl) :: velh(3,3)
          
          REAL(dbl) :: fcell(3,3)

          IF( wc == 0.0d0 ) &
            CALL errore( ' movecell ',' cell mass is 0 ! ', 1 )

          ! Force on the cell
          !
          fcell = box_t0%pail(:,:)
          fcell = fcell - box_t0%omega * box_t0%m1 * press
          fcell = fcell / wc

          CALL cell_move( box_tp1%hmat, box_t0%hmat, box_tm1%hmat, delt, &
            iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )

          box_tp1%a = TRANSPOSE( box_tp1%hmat(:,:) )
          CALL gethinv( box_tp1 )
          box_tp1%g    = MATMUL( box_tp1%a(:,:), box_tp1%hmat(:,:) )
          box_tp1%gvel = ( box_tp1%g(:,:) - box_tm1%g(:,:) ) / ( 2.0d0 * delt )
        
          RETURN
        END SUBROUTINE MOVECELL


!
!------------------------------------------------------------------------------!
   END MODULE cell_module
!------------------------------------------------------------------------------!
