!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE punch_plot_ph()
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output the change of the charge density,
  !     due to a perturbation ubare defined in the routine phq_setup. It
  !     can be read by chdens.f which cuts a bidimensional plane
  !     to plot contour levels,
  !     or selects a line for a usual line plot. The routine produces
  !     a file with the name in the variable fildrho# given in input.
  !
  USE kinds,      ONLY : DP
  USE fft_base, ONLY : dfftp
  USE gvect,      ONLY : gcutm
  USE gvecs,    ONLY : dual
  USE cell_base,  ONLY : ibrav, celldm
  USE lsda_mod,   ONLY : nspin, lsda
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, atm, zv, tau
  USE io_global,  ONLY : stdout, ionode
  USE fft_base,   ONLY : grid_gather
  USE output,     ONLY : fildrho
  USE units_ph,   ONLY : iudrho, lrdrho
  USE modes,      ONLY : ubar, u, npert, npertx, nirr
  USE partial,    ONLY : comp_irr
  USE wvfct,      ONLY : ecutwfc
  USE printout_base, ONLY : title

  IMPLICIT NONE

  INTEGER :: iunplot, ios, ipert, irr, na, ir, nt, imode0, plot_num, jpol
  ! unit of the plot file
  ! integer variable for I/O contr
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on mesh points
  ! counter on atom types
  ! the starting mode
  ! compatibility variable
  ! counter on polarization

  CHARACTER(len=80) :: filin
  ! complete name of the file

  REAL(DP), ALLOCATABLE :: raux (:)
  ! auxiliary vector

  COMPLEX(DP) :: ps
  COMPLEX(DP), EXTERNAL :: zdotc
  COMPLEX(DP), ALLOCATABLE :: aux (:,:,:), aux1 (:,:)
  ! the scalar product
  ! scalar product function
  ! auxiliary space to rotate the
  ! induced charge
#if defined (__PARA)
  ! auxiliary vector
  REAL(DP), ALLOCATABLE :: raux1 (:)
#endif

  IF (fildrho.EQ.' ') RETURN
  WRITE( stdout, '(/5x,"Calling punch_plot_ph" )')
  WRITE( stdout, '(5x,"Writing on file  ",a)') fildrho
  !
  !    reads drho from the file
  !
  ALLOCATE (aux  (  dfftp%nnr,nspin,npertx))
  ALLOCATE (aux1 (  dfftp%nnr,nspin))
  ALLOCATE (raux (  dfftp%nnr))
  !
  !
  !     reads the delta_rho on the aux variable
  !
  aux1(:,:) = (0.d0, 0.d0)
  imode0 = 0
  DO irr = 1, nirr
     IF (comp_irr (irr) .EQ.1) THEN
        DO ipert = 1, npert (irr)
           CALL davcio_drho (aux (1, 1, ipert), lrdrho, iudrho, imode0 + &
                ipert, - 1)
        ENDDO
#if defined (__PARA)
        CALL psymdvscf (npert (irr), irr, aux)
#else
        CALL symdvscf (npert (irr), irr, aux)
#endif
        DO ipert = 1, npert (irr)
           ps = zdotc (3 * nat, ubar, 1, u (1, imode0 + ipert), 1)
           CALL zaxpy (dfftp%nnr * nspin, ps, aux (1, 1, ipert), 1, aux1, 1)
        ENDDO
     ENDIF
     imode0 = imode0 + npert (irr)
  ENDDO
  !
  !     write on output the change of the charge
  !
  iunplot = 4
  filin = TRIM(fildrho)
  !
  IF ( ionode ) THEN
     !
     OPEN (unit = iunplot, file = filin, status = 'unknown', err = &
          100, iostat = ios)

100  CALL errore ('plotout', 'opening file'//filin, ABS (ios) )
     REWIND (iunplot)
     !
     !       Here we write some information quantity which are always necessa
     !
     plot_num = 0
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
  CALL dcopy (dfftp%nnr, aux1 (1, 1), 2, raux, 1)

  IF (lsda) CALL daxpy (dfftp%nnr, 1.d0, aux1 (1, 2), 2, raux, 1)

#if defined (__PARA)
  ALLOCATE (raux1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
  CALL grid_gather (raux, raux1)
  IF ( ionode ) WRITE (iunplot, * ) (raux1 (ir), ir = 1, dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
  DEALLOCATE (raux1)
#else
  WRITE (iunplot, * ) (raux (ir), ir = 1, dfftp%nnr)
#endif

  IF (ionode) CLOSE (unit = iunplot)
  DEALLOCATE (raux)
  DEALLOCATE (aux1)
  DEALLOCATE (aux)
  RETURN
END SUBROUTINE punch_plot_ph
