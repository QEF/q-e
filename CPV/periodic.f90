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
         metric_print_info, updatecell, dgcell, &
         movecell, r_to_s, s_to_r,  &
         get_lattice_vectors, pbcs, get_celldm, &
         cell_init, alat, press, at


!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

        subroutine metric_print_info( unit )

          USE constants, ONLY: au_gpa
          USE control_flags, ONLY: thdyn, tsdc, tzeroc, tbeg
          USE control_flags, ONLY: tnoseh
          USE input_parameters, ONLY: rd_ht
          USE io_global, ONLY: stdout
          USE cell_base, ONLY: press, frich, greash
          USE cell_nose, ONLY: temph, qnh



          INTEGER, INTENT(IN) :: unit
          INTEGER :: i,j
          REAL(dbl) :: htm1(3,3)
          REAL(dbl) :: h(3,3)
          REAL(dbl) :: ht(3,3)
          REAL(dbl) :: omega

          IF ( tbeg ) THEN
            ht = rd_ht
          ELSE
            DO i = 1, 3
              ht(1,i) = a1(i)
              ht(2,i) = a2(i)
              ht(3,i) = a3(i)
            END DO
          END IF

          h = TRANSPOSE( ht )
          CALL invmat (3, ht, htm1, omega )

          WRITE(unit,545) 
          IF ( tbeg ) THEN
            WRITE(unit,546) 
          ELSE
            WRITE(unit,547) 
          END IF
          WRITE(unit,550) ibrav, omega
          WRITE(unit,551)

! ...     The matrix "ht" in FPMD correspond to the transpose of matrix "at" in PW
!
          DO I = 1, 3
            WRITE(unit,555) i, ( ht(i,j), j = 1, 3 )
          END DO
          WRITE(unit,552)

! ...     The matrix "htm1" in FPMD correspond to the matrix "bg" in PW
!
          DO I = 1, 3
            WRITE(unit,556) i, ( htm1(j,i), j = 1, 3 )   
          END DO

          IF( .NOT. thdyn ) THEN
            WRITE( stdout,525)
          ELSE
            IF( tsdc ) THEN
              WRITE( stdout,526)
            ELSE
              IF( frich /= 0.0d0 ) THEN
                WRITE( stdout,528) frich
              ELSE
                WRITE( stdout,527)
              END IF
            END IF
            WRITE( stdout,530) press * au_gpa, wc
            IF(tzeroc) THEN
              WRITE( stdout,563)
            ENDIF
            WRITE( stdout,565)
          END IF

        if( thdyn ) then
          WRITE( stdout,600)
          if( thdyn ) then
            ! if( thdiag ) WRITE( stdout,608)
            if( tnoseh ) then
               frich = 0.0d0  ! if nose is in use friction must be 0
               WRITE( stdout,604) temph, qnh, press * au_gpa
            else
               WRITE( stdout,602) frich, greash, press * au_gpa
            endif
          else
            WRITE( stdout,606)
          endif
        endif

 600  format( 3X, 'internal stress tensor calculated')
 602  format( 3X, 'cell parameters dynamics with frich = ',f7.4,            &
     &        3X, 'and greash = ',f7.4,/                                    &
     &        3X, 'external pressure = ',f11.7,'(gpa)'//)
 604  format( 3X, 'cell parameters dynamics with nose` temp. control:'/     &
     &        3X, 'cell temperature required = ',f10.5,'(kelvin)',          &
     &        3X, 'nose` mass = ',f10.3,/                                   &
     &        3X, 'external pressure = ',f11.7,'(gpa)'//)
 606  format( 3X, 'cell parameters are not allowed to move'//)
 608  format( 3X, 'frozen off-diagonal cell parameters'//)


  545     FORMAT(//,3X,'Cell Dynamics Parameters (from STDIN)',/ &
                   ,3X,'-------------------------------------')
  550     FORMAT(   3X,'Lattice index (ibrav) = ',I3,' Volume = ',F12.4)
  546     FORMAT(   3X,'Simulation cell read from STDIN')
  547     FORMAT(   3X,'Starting cell generated from CELLDM')
  551     FORMAT(   3X,'Ht0  : ')
  552     FORMAT(   3X,'Htm1 : ')
  555     FORMAT(   3X,'A',I1,' = ',4X,3(1X,F8.4))
  556     FORMAT(   3X,'B',I1,' = ',4X,3(1X,F8.4))
  525     FORMAT(   3X,'Constant VOLUME Molecular dynamics')
  526     FORMAT(   3X,'Volume dynamics with steepest descent')
  527     FORMAT(   3X,'Volume dynamics with newton equations')
  528     FORMAT(   3X,'Volume dynamics with verlet and friction ',F8.4) 
  530     FORMAT(   3X,'Constant PRESSURE Molecular dynamics:',/ &
                   ,3X,'External pressure (GPa) = ',F11.2,/ &
                   ,3X,'Volume mass             = ',F11.2)
  563     FORMAT(   3X,'Zero initial momentum for cell variables')
  565     FORMAT(   3X,'Volume dynamics: ', &
                       'the temperature is not controlled')

          return
        end subroutine metric_print_info

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
          box_tp1%g = MATMUL( box_tp1%a(:,:), box_tp1%hmat(:,:) )
        
          RETURN
        END SUBROUTINE MOVECELL


!
!------------------------------------------------------------------------------!
   END MODULE cell_module
!------------------------------------------------------------------------------!
