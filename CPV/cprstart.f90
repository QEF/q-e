!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

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
      use input, only: read_input_file, iosys_pseudo
      use io_global, ONLY: io_global_start, io_global_getionode
      use mp_global, ONLY: mp_global_start
      use mp, ONLY: mp_end, mp_start, mp_env
      use io_files, only: psfile, pseudo_dir
      use ions_base, only: nsp
      use control_flags, only: lneb, lsmd, lwf, program_name
      use environment, ONLY: environment_start, environment_end
!
      implicit none

      INTEGER :: mpime, nproc, gid, root, ionode_id
      LOGICAL :: ionode

!
!     program starts here
!

! ... initialize MPI (parallel processing handling)

      root = 0
      CALL mp_start()
      CALL mp_env( nproc, mpime, gid )
      CALL mp_global_start( root, mpime, gid, nproc )

! ... mpime = processor number, starting from 0
! ... nproc = number of processors
! ... gid   = group index
! ... root  = index of the root processor

      program_name = 'CPVC'

! ... initialize input output

      CALL io_global_start( mpime, root )
      CALL io_global_getionode( ionode, ionode_id )


      CALL environment_start( )

      !
      !  readin the input file
      !

      call read_input_file( lneb, lsmd, lwf )

      !
      !  copy pseudopotential input parameter into internal variables
      !  and read in pseudopotentials and wavefunctions files
      !

      call iosys_pseudo( psfile, pseudo_dir, nsp )

      IF( lneb ) THEN
        call neb_loop( 1 )
      ELSE IF ( lsmd ) THEN
        call smd_loop( 1 )
      ELSE IF ( lwf ) THEN
        call wf_loop( 1 )
      ELSE
        call cpr_loop( 1 )
      END IF

      CALL environment_end( )

      call mp_end()

      stop
      end program
