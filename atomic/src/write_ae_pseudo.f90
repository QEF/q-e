!
!---------------------------------------------------------------------
subroutine write_ae_pseudo
  !---------------------------------------------------------------------
  !
  ! This routine generates a UPF file with a local Coulomb potential
  ! and the all-electron starting wave-functions. It allows to make
  ! an all-electron calculation with plane-waves.
  !
  use kinds,     only : dp
  use io_global, only : ionode, ionode_id
  use mp,        only : mp_bcast
  use ld1inc, only : file_pseudopw, rel, grid, &
                     etot, zed, nwf,  el,  ll,  oc,  psi,  rho, &
                     etots,zval,nwfts,elts,llts,octs,phits,rhos, &
                     nwfs, lloc, rcloc, iswitch, nlcc, lpaw, ecutrho, &
                     ecutwfc, nbeta, lmax
  implicit none

  integer :: ios,   &  ! I/O control
             iunps     ! the unit with the pseudopotential
  character (len=2) :: atom
  character (len=2), external :: atom_name

  IF (iswitch /= 1 ) call errore('write_ae_pseudo','wrong iswitch',1)
  
  atom = atom_name(nint(zed))
  IF ( atom(1:1) == ' ' ) THEN
     file_pseudopw = atom(2:2) // '.UPF'
  ELSE
     file_pseudopw = TRIM(atom) // '.UPF'
  END IF
!  iunps=28
!  if (ionode) &
!     open(unit=iunps, file=trim(file_pseudopw), status='unknown',  &
!          form='formatted', err=50, iostat=ios)
!50  call mp_bcast(ios, ionode_id)
!  call errore('write_ae_pseudo','opening file_pseudopw',abs(ios))

  if ( rel==2 ) call errore('write_ae_pseudo','you cannot be serious!!!',rel)
  if (ionode) then
    !
    lloc = 0
    rcloc = 0.0_DP
    nwfs = 0
!    call write_pseudo_comment(iunps)  
    zval = zed
    etots= etot
    nwfts = nwf
    nbeta = 0
    nlcc = .false.
    ecutwfc=0.0_DP
    ecutrho=0.0_DP
    lpaw=.false.
    lmax=0
    elts(1:nwfts) = el(1:nwf)
    llts(1:nwfts) = ll(1:nwf)
    octs(1:nwfts) = oc(1:nwf)
    phits(1:grid%mesh,1:nwfts) = psi(1:grid%mesh,1,1:nwf )
!    call write_pseudo_header(iunps)  
!    call write_pseudo_mesh(iunps)
!    call write_pseudo_pswfc(iunps)
    rhos(1:grid%mesh,1)= rho(1:grid%mesh,1) 
!    call write_pseudo_rhoatom(iunps)  
    call ld1_writeout()
    !
!    close(iunps)
    !
  endif
  !
  return
end subroutine write_ae_pseudo
