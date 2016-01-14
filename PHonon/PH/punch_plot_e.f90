!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE punch_plot_e()
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output the change of the charge density,
  !     due to an electric field in a real space mesh which can be read
  !     by chdens.f which cuts a bidimensional plane to plot contour level
  !     or selects a line for a usual line plot. The routine produces
  !     3 files with the change of charge density due to perturbations in
  !     three cartesian directions. The names of the files are
  !     in the variable fildrho given in input.
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, atm, zv, tau
  USE io_global,  ONLY : stdout, ionode
  USE run_info, ONLY : title
  USE fft_base,   ONLY : dfftp
  USE scatter_mod,   ONLY : gather_grid
  USE gvect,      ONLY : gcutm
  USE gvecs,    ONLY : dual
  USE cell_base,  ONLY : bg, ibrav, celldm
  USE lsda_mod,   ONLY : lsda
  USE noncollin_module, ONLY : nspin_mag
  USE output,     ONLY : fildrho
  USE units_ph,   ONLY : iudrho, lrdrho
  USE gvecw,      ONLY : ecutwfc
  IMPLICIT NONE

  INTEGER :: plot_num, iunplot, ios, ipol, jpol, na, ir, nt
  ! type of plot (not used)
  ! unit of the plot file
  ! integer variable for I/O contr
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on mesh points

  CHARACTER :: caux * 1, filin * 80, which*2
  ! used to compose the name
  ! complete name of the file

  REAL(DP), ALLOCATABLE  :: raux (:)
  ! auxiliary vector

  COMPLEX(DP), ALLOCATABLE :: aux (:,:,:), aux1 (:,:,:)
  ! auxiliary space to rotate the
  ! induced charge

#if defined (__MPI)
  ! auxiliary vector
  REAL(DP), ALLOCATABLE :: raux1 (:)
#endif

  IF (fildrho.EQ.' ') RETURN
  WRITE( stdout, '(/5x,"Calling punch_plot_e" )')
  WRITE( stdout, '(5x,"Writing on file  ",a)') fildrho
  !
  !    reads drho from the file
  !
  ALLOCATE (aux  (dfftp%nnr,nspin_mag,3))
  ALLOCATE (aux1 (dfftp%nnr,nspin_mag,3))
  ALLOCATE (raux (dfftp%nnr))
  !
  !     reads the delta_rho on the aux variable
  !
  DO ipol = 1, 3
     CALL davcio_drho (aux (1,1,ipol), lrdrho, iudrho, ipol, - 1)
  ENDDO
  !
  !     rotate the charge and transform to cartesian coordinates
  !
  aux1(:,:,:) = (0.0d0, 0.0d0)
  DO ipol = 1, 3
     DO jpol = 1, 3
        CALL daxpy (2 *dfftp%nnr, bg (ipol, jpol), aux (1,1,jpol), 1, &
             aux1 (1,1,ipol), 1)
     ENDDO
  ENDDO
  !
  !     write on output the change of the charge
  !
  iunplot = 4
  which='_e'
  DO ipol = 1, 3
     WRITE (caux, '(i1)') ipol
     filin = TRIM(fildrho) //which//caux
     !
     IF ( ionode ) THEN
        !
        OPEN (unit = iunplot, file = filin, status = 'unknown', err = &
             100, iostat = ios)

100     CALL errore ('plotout', 'opening file'//filin, ABS (ios) )
        REWIND (iunplot)
        !
        !    Here we write some needed quantities
        !
        ! not used
        plot_num = - 1
        WRITE (iunplot, '(a)') title
        WRITE (iunplot, '(8i8)') dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, &
             ntyp
        WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
        WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, plot_num
        WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
                                (nt, atm (nt), zv (nt), nt=1, ntyp)
        WRITE (iunplot, '(i4,3x,3f14.10,3x,i2)') (na, &
          (tau (jpol, na), jpol = 1, 3), ityp (na), na = 1, nat)
        !
     ENDIF
     !
     !      plot of the charge density
     !
     raux (:) =  DBLE (aux1 (:,1, ipol) )
     IF (lsda) CALL daxpy (dfftp%nnr, 1.d0, aux1 (1,2, ipol), 2, raux, 1)
     !
#if defined (__MPI)
     ALLOCATE (raux1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
     CALL gather_grid (dfftp, raux, raux1)
     IF ( ionode ) WRITE (iunplot, '(5(1pe17.9))') &
          (raux1 (ir) , ir = 1, dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
     DEALLOCATE (raux1)
#else
     WRITE (iunplot, '( 5( 1pe17.9 ) )') (raux (ir) , ir = 1, dfftp%nnr)
#endif
     IF (ionode) CLOSE (unit = iunplot)
  ENDDO
  DEALLOCATE (raux)
  DEALLOCATE (aux1)
  DEALLOCATE (aux)
  RETURN

END SUBROUTINE punch_plot_e
