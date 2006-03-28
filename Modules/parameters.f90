!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------------
MODULE parameters
  !---------------------------------------------------------------------------
  !
  IMPLICIT NONE
  SAVE
  !
  INTEGER, PARAMETER :: &
       ntypx  = 10,     &! max number of different types of atom
       npsx   = ntypx,  &! max number of different PPs (obsolete)
       npk    = 40000,  &! max number of k-points               
       lmaxx  = 3,      &! max non local angular momentum       
       nchix  = 6,      &! max number of atomic wavefunctions per atom
       ndmx   = 2000     ! max number of points in the atomic radial mesh
  !
  INTEGER, PARAMETER :: &
       nbrx = 14,          &! max number of beta functions
       lqmax= 2*lmaxx+1,   &! max number of angular momenta of Q
       nqfx = 8             ! max number of coefficients in Q smoothing
  !
  INTEGER, PARAMETER :: cp_lmax = lmaxx + 1  ! max number of channels
                                             ! (s,p,d,f)
  INTEGER, PARAMETER :: nacx    = 10         ! max number of averaged 
                                             ! quantities saved to the restart
  INTEGER, PARAMETER :: nsx     = ntypx      ! max number of species
  INTEGER, PARAMETER :: natx    = 5000       ! max number of atoms
  INTEGER, PARAMETER :: npkx    = npk        ! max number of K points
  INTEGER, PARAMETER :: ncnsx   = 101        ! max number of constraints
  INTEGER, PARAMETER :: nspinx  = 2          ! max number of spinors
  !
  INTEGER, PARAMETER :: nhclm   = 4  ! max number NH chain length, nhclm can be
                                     ! easily increased since the restart file 
                                     ! should be able to handle it, perhaps
                                     ! better to align nhclm by 4
  !
  INTEGER, PARAMETER :: max_num_of_images = 50 ! max number of images in "path"
                                               ! calculations ( NEB or SMD )
  !
  INTEGER, PARAMETER :: max_nconstr = 100
  !
  INTEGER, PARAMETER  ::  MAXCPU = 2**17  ! Maximum number of CPU
  INTEGER, PARAMETER  ::  MAXGRP = 128    ! Maximum number of task-groups
  !
END MODULE parameters
