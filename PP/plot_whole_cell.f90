subroutine plot_whole_cell (alat, at, nat, tau, atm, ityp, &
     nr1, nr2, nr3, nrx1, nrx2, nrx3, rho, output_format, ounit)
  USE kinds, only : DP
  implicit none
  integer          :: nat, ityp (nat), output_format, ounit
  integer          :: nrx1, nrx2, nrx3, nr1, nr2, nr3
  character(len=3) :: atm(*)
  real(kind=DP)    :: alat, tau (3, nat), at (3, 3), rho(2, nrx1,nrx2,nrx3)

  if ( output_format .eq. 3 ) then
     !
     ! XCRYSDEN FORMAT
     !
     call xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
     call xsf_fast_datagrid_3d &
          (rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, at, alat, ounit)

  elseif ( output_format .eq. 4 ) then
     !
     ! gOpenMol format
     !

     ! not yet implemented
     ! add code here ...
  else
     call errore('plot_whole_cell', 'wrong output_format', 1)
  endif
end subroutine plot_whole_cell
