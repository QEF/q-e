!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/f_defs.h"

!
!=======================================================================
!
                        program main
!
!=======================================================================
!***  Molecular Dynamics using Density-Functional Theory   ****
!***  this is a Car-Parrinello program using Vanderbilt pseudopotentials
!***********************************************************************
!***  based on version 11 of cpv code including ggapw 07/04/99
!***  copyright Alfredo Pasquarello 10/04/1996
!***  parallelized and converted to f90 by Paolo Giannozzi (2000),
!***  using parallel FFT written for PWSCF by Stefano de Gironcoli
!***  PBE added by Michele Lazzeri (2000)
!***  variable-cell dynamics by Andrea Trave (1998-2000)
!***********************************************************************
!***  appropriate citation for use of this code:
!***  Car-Parrinello method    R. Car and M. Parrinello, PRL 55, 2471 (1985) 
!***  current implementation   A. Pasquarello, K. Laasonen, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRL 69, 1982 (1992);
!***                           K. Laasonen, A. Pasquarello, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRB 47, 10142 (1993).
!***  implementation gga       A. Dal Corso, A. Pasquarello, A. Baldereschi,
!***                           and R. Car, PRB 53, 1180 (1996).
!***********************************************************************
!***  
!***  f90 version, with dynamical allocation of memory
!***  Variables that do not change during the dynamics are in modules
!***  (with some exceptions) All other variables are passed as arguments
!***********************************************************************
!***
!*** fft : uses machine's own complex fft routines, two real fft at the time
!*** ggen: g's only for positive halfspace (g>)
!*** all routines : keep equal c(g) and [c(-g)]*
!***
!***********************************************************************
!    general variables:
!     delt           = delta t
!     emass          = electron mass (fictitious)
!     dt2bye         = 2*delt/emass
!***********************************************************************
!
      use input_cp, only: read_input_file, iosys_pseudo
      use io_global, ONLY: io_global_start
      use mp_global, ONLY: mp_global_start
      use mp, ONLY: mp_end
      use para_mod, ONLY: me, mygroup, nproc
      use io_files, only: psfile, pseudo_dir
      use ions_base, only: nsp
      use control_flags, only: lneb
!
      implicit none
!
!     program starts here
!
!     Initialize processors IDs
      call startup()
      call io_global_start( (me-1), 0 )
      call mp_global_start(0, (me-1), mygroup, nproc )

      CALL init_clocks( .TRUE. )
      CALL start_clock( 'CP' )

      !
      !  readin the input file
      !

      call read_input_file( lneb )

      !
      !  copy pseudopotential input parameter into internal variables
      !  and read in pseudopotentials and wavefunctions files
      !

      call iosys_pseudo( psfile, pseudo_dir, nsp )

      IF( lneb ) THEN
        call neb_loop( 1 )
      ELSE
        call cpr_loop( 1 )
      END IF

      CALL stop_clock( 'CP' )
      CALL print_clock( 'CP' )

      call mp_end()

      stop
      end program


    SUBROUTINE cpr_loop( nloop )
      USE input_parameters, ONLY: nat, tprnfor
      USE input_parameters, ONLY: ion_positions, rd_pos
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nloop
      INTEGER :: iloop
      REAL(kind=8), ALLOCATABLE :: tau( :, : )
      REAL(kind=8), ALLOCATABLE :: fion( :, : )
      REAL(kind=8) :: etot

      IF( nat > 0 ) THEN
        ALLOCATE( tau( 3, nat ), fion( 3, nat ) )
      ELSE
        CALL errore( ' cpr_loop ', ' nat less or equal 0 ', 1 )
      END IF

      ! ... set tprnfor = .true. to get atomic forces
      ! ... even if the atoms do not move

      ! ... set ion_positions = 'from_input'
      ! ... and rd_pos = +your_positions+
      ! ... to force cprmain to compute forces for  
      ! ... +your_position+ configuration

      DO iloop = 1, nloop
        call cprmain( tau(1,1), fion(1,1), etot)
        call memstat( 1 )
      END DO

      DEALLOCATE( tau, fion )

      RETURN
    END SUBROUTINE


   SUBROUTINE neb_loop( iloop )

     USE kinds
     USE io_global,        ONLY: ionode, stdout
     USE neb_variables,    ONLY: conv_neb
     USE neb_variables,    ONLY: neb_deallocation
     USE neb_routines,     ONLY: initialize_neb, search_mep, iosys_neb
     USE io_routines,      ONLY: write_output
     USE ions_base,        ONLY: deallocate_ions_base

     IMPLICIT NONE

     INTEGER :: iloop

     ! ... stdout is connected to a file ( specific for each image )
     ! ... via unit 17

     IF( ionode ) THEN
       !
       stdout = 17
       !
     END IF

     CALL iosys_neb()

     CALL initialize_neb( 'FP' )
     !
     ! ... this routine does all the NEB job
     !
     CALL search_mep()
     !
     ! ... output is written
     !
     CALL write_output()

     CALL deallocate_ions_base()
     !
     ! ... stdout is reconnected to standard output
     !
     stdout = 6
     !

     RETURN
   END SUBROUTINE neb_loop
