!
! Copyright (C) 2004-2007 Quantm-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
program ld1
  !---------------------------------------------------------------
  !
  !     atomic self-consistent local-density program
  !     atomic rydberg units are used : e^2=2, m=1/2, hbar=1
  !     psi(r) = rR(r), where R(r) is the radial part of the wfct
  !     rho(r) = psi(r)^2 => rho(r) = (true charge density)*(4\pi r^2)
  !                       The same applies to the core charge
  !---------------------------------------------------------------
  !
  USE global_version,    ONLY : version_number
  USE io_files,          ONLY : nd_nmbr
  USE mp,                ONLY : mp_barrier, mp_end
  USE ld1inc,            ONLY : iswitch, write_coulomb
  !
  implicit none
  character :: day*9, hour*9
  CHARACTER (LEN=9) :: code = 'LD1'
  !
  !   write initialization information
  !
  call startup( nd_nmbr, code, version_number )
  !
  !    read input, possible pseudopotential and set the main variables
  !
  call ld1_readin ( )
  call ld1_setup ( )
  !
  !   three possible working mode:
  !
  if (iswitch.eq.1) then
     !
     !   all-electron calculation
     !
     call all_electron(.true.,1)
     if ( write_coulomb ) call write_fake_pseudo ( )
     !
  elseif (iswitch.eq.2) then
     !
     !   pseudopotential test
     !
     call run_test ( )
     call ld1_writeout ( )
     !
  elseif (iswitch.eq.3) then
     !
     !  pseudopotential generation and test
     !
     call all_electron(.false.,1)
     call gener_pseudo ( )
     call run_test ( )
     call ld1_writeout ( )
     !
  else
     call errore('ld1','iswitch not implemented',1)
  endif
  call mp_barrier()
  call mp_end()

end program ld1
!
!---------------------------------------------------------------------
subroutine write_fake_pseudo
  !---------------------------------------------------------------------
  !
  use io_global, only : ionode, ionode_id
  use mp,        only : mp_bcast
  use ld1inc, only : file_pseudopw, rel, grid, &
                     etot, zed, nwf,  el,  ll,  oc,  psi,  rho, &
                     etots,zval,nwfts,elts,llts,octs,phits,rhos
  implicit none

  integer :: ios,   &  ! I/O control
             iunps     ! the unit with the pseudopotential
  character (len=2) :: atom
  character (len=2), external :: atom_name
  
  atom = atom_name(nint(zed))
  IF ( atom(1:1) == ' ' ) THEN
     file_pseudopw = atom(2:2) // '.UPF'
  ELSE
     file_pseudopw = TRIM(atom) // '.UPF'
  END IF
  iunps=28
  if (ionode) &
     open(unit=iunps, file=trim(file_pseudopw), status='unknown',  &
          form='formatted', err=50, iostat=ios)
50  call mp_bcast(ios, ionode_id)
  call errore('write_fake_pseudo','opening file_pseudopw',abs(ios))

  if ( rel==2 ) call errore('write_fake_pseudo','you cannot be serious!!!',rel)
  if (ionode) then
    !
    call write_pseudo_comment(iunps)  
    zval = zed
    etots= etot
    call write_pseudo_header(iunps)  
    call write_pseudo_mesh(iunps)
    nwfts = nwf
    elts(1:nwfts) = el(1:nwf)
    llts(1:nwfts) = ll(1:nwf)
    octs(1:nwfts) = oc(1:nwf)
    phits(1:grid%mesh,1:nwfts) = psi(1:grid%mesh,1,1:nwf )
    call write_pseudo_pswfc(iunps)
    rhos(1:grid%mesh,1)= rho(1:grid%mesh,1) 
    call write_pseudo_rhoatom(iunps)  
    !
    close(iunps)
    !
  endif
  !
  return
end subroutine write_fake_pseudo
