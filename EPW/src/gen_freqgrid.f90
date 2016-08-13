  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !
  !-----------------------------------------------------------------------
  SUBROUTINE gen_freqgrid_raxis
  !-----------------------------------------------------------------------
  !!
  !! Automatic generation of the frequency-grid for real-axis calculations.
  !!
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : nswfc, nswc, pwc, wsfc, wscut, lunif
  USE eliashbergcom, ONLY : nsw, ws, dws
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw
  !
  ! define a grid ws in 2 step sizes
  ! 1. a fine grid of nswfc points between (0,wsfc)
  ! 2. a rough grid of nswc points between (wsfc,wscut).
  ! above wsfc the gap function varies slowly
  !
  ! nswfc = nr. of grid points between (0,wsfc)
  ! nswc  = nr. of grid points between (wsfc,wscut)
  !
  WRITE(stdout,'(a)') '    '
  WRITE(stdout,'(5x,a,i6,a)') 'Total number of nsw = ', nsw, ' grid-points are divided in:'
  WRITE(stdout,'(5x,a,i6,a,f12.6,a,f12.6)') 'nswfc = ', nswfc, '  from ', 0.0, ' to ', wsfc  
  WRITE(stdout,'(5x,a,i6,a,f12.6,a,f12.6)') 'nswc  = ', nswc,  '  from ', wsfc, ' to ', wscut
  WRITE(stdout,'(a)') '    '
  !
  IF ( .not. ALLOCATED(ws) )  ALLOCATE( ws(nsw) )
  IF ( .not. ALLOCATED(dws) ) ALLOCATE( dws(nsw) )
  ws(:) = 0.d0
  dws(:) = 0.d0
  !
  DO iw = 1, nswfc
     dws(iw) = wsfc / dble(nswfc)
     ws(iw) = dble(iw) * dws(iw)
  ENDDO
  DO iw = nswfc + 1, nsw
     dws(iw) = ( wscut - wsfc ) / dble(nswc)
     IF ( lunif ) THEN 
        ws(iw) = wsfc + dble(iw) * dws(iw)
     ELSE 
        ! RM this needs to be checked
        ws(iw) = wsfc + dble( iw/nswc )**pwc * (wscut - wsfc)
     ENDIF
  ENDDO
  !
  IF ( .not. lunif ) THEN 
     DO iw = nswfc+1, nsw-1
        dws(iw) = ws(iw+1) - ws(iw)
     ENDDO
     dws(nsw) = dws(nsw-1)
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE gen_freqgrid_raxis
  !
  !-----------------------------------------------------------------------
  SUBROUTINE gen_freqgrid_iaxis( itemp )
  !-----------------------------------------------------------------------
  !
  ! Automatic generation of the frequency-grid for imaginary-axis calculations.
  !
  !
  ! input
  !
  ! itemp  - temperature point
  !
  USE constants,     ONLY : pi
  USE epwcom,        ONLY : nqstep, lpade, lacon 
  USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, wsph, dwsph, estemp, wsphmax
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, n, imelt
  !
  ! frequency-grid for imaginary-axis
  ! nsiw(itemp) = nr. of grid points between (0,wscut) 
  !
  ! memory allocated for wsi and ws
  imelt = nsiw(itemp) + nsw 
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(wsi) )  ALLOCATE( wsi(nsiw(itemp)) )
  wsi(:) = 0.d0
  DO iw = 1, nsiw(itemp)
     n = iw - 1
     wsi(iw) = dble(2*n+1) * pi * estemp(itemp) 
     !WRITE(*,*) iw, wsi(iw)
  ENDDO
  !
  ! frequency-grid for real-axis ( Pade approximants and analytic continuation)
  !
  IF ( lpade .OR. lacon ) THEN
     IF ( .not. ALLOCATED(ws) )  ALLOCATE( ws(nsw) )
     ws(:) = 0.d0
     DO iw = 1, nsw
        IF ( iw .le. nqstep ) THEN 
           ws(iw) = wsph(iw)
        ELSE
           ws(iw) = wsphmax + dble(iw-nqstep)*dwsph
        ENDIF
        !WRITE(*,*) iw, ws(iw), wsph(iw)
     ENDDO
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE gen_freqgrid_iaxis
