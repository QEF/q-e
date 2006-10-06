!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

        SUBROUTINE movecell_x( tsdc, box_tm1, box_t0, box_tp1, velh )

          USE kinds,         ONLY: DP
          USE time_step,     ONLY: delt
          USE cell_base,     ONLY: iforceh, cell_move, frich, gethinv
          USE cell_base,     ONLY: wc => wmass, press, boxdimensions
          USE cell_nose,     ONLY: vnhh
          USE control_flags, ONLY: tnoseh

          IMPLICIT NONE

          LOGICAL :: tsdc
          TYPE (boxdimensions) :: box_tm1, box_t0, box_tp1
          REAL(DP) :: velh(3,3)
          
          REAL(DP) :: fcell(3,3)

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
        END SUBROUTINE movecell_x
