!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!     
!---------------------------------------------------------------------
subroutine ld1_writeout
  !---------------------------------------------------------------------
  !
  !     This routine writes on output the quantities which defines 
  !     a multiprojector pseudopotential. It can be in the
  !     Vanderbilt form or in the norm-conserving form
  !
  use radial_grids, only: ndmx
  use io_global, only : qestdin, stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use mp_world,  only : world_comm
  use ld1inc, only : file_pseudopw, zed, grid, &
                     nconf , lpaw, rel, pawsetup, pseudotype, &
                     rhoc, vnl, phits, vpsloc, & 
                     elts, llts, octs, rcut, etots, nwfts, &
                     lmax, lloc, zval, nlc, nnl, alps, alpc, alc, cc, nlcc, &
                     iswitch
  use funct, only : get_dft_name
  use paw_type, only : deallocate_pseudo_paw
  use open_close_input_file, only: close_input_file

  implicit none

  integer :: &
       ios,   &  ! I/O control
       iunps     ! the unit with the pseudopotential

  logical, external :: matches
  logical :: oldformat
  character (len=20) :: dft_name
  
  if (file_pseudopw == ' ') return

  if (iswitch/=4.and.nconf > 1) &
       call errore('ld1_writeout','more than one test configuration',1)

  if ( (( rel == 2) .or. lpaw) &
       .and. .not. matches('.UPF',file_pseudopw) &
       .and. .not. matches('.upf',file_pseudopw) ) then
     file_pseudopw=trim(file_pseudopw)//'.UPF'
  end if

  oldformat = .not. matches('.UPF',file_pseudopw) .and. &
              .not. matches('.upf',file_pseudopw)

  iunps=28
  if (ionode) &
     open(unit=iunps, file=trim(file_pseudopw), status='unknown',  &
          form='formatted', err=50, iostat=ios)
50  call mp_bcast(ios, ionode_id, world_comm)
  call errore('ld1_writeout','opening file_pseudopw',abs(ios))

  if (ionode) then
     if (oldformat) then
        !
        if (pseudotype == 1) then
          dft_name = get_dft_name()
          !
          ! write in CPMD format 
          if ( matches('.psp',file_pseudopw) ) then
             call write_cpmd &
                  (iunps,zed,grid%xmin,grid%dx,grid%mesh,ndmx,grid%r,grid%r2,  &
                  dft_name,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
                  rhoc,vnl,phits,vpsloc,elts,llts,octs,rcut,etots,nwfts)
          else
          !
          ! write old "NC" format (semilocal)
          !
             call write_pseudo &
                  (iunps,zed,grid%xmin,grid%dx,grid%mesh,ndmx,grid%r,grid%r2,  &
                  dft_name,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
                  rhoc,vnl,phits,vpsloc,elts,llts,octs,etots,nwfts)
          end if
        else
          !
          ! write old "RRKJ" format (nonlocal)
          !
           call write_rrkj ( iunps )
        end if
        !
     else
        !
        call export_upf(iunps, qestdin)
        !
        if(lpaw) call deallocate_pseudo_paw( pawsetup )
        !
     endif
     !
     close(iunps)
     ! close input data unit if not done previously
     ios = close_input_file ( )
  endif
  !
  return
end subroutine ld1_writeout

!---------------------------------------------------------------------
subroutine write_rrkj (iunps)
  !---------------------------------------------------------------------
  !
  use ld1inc, only : title, pseudotype, rel, nlcc, zval, etots, lmax, &
                     els, nns, lls, rcut, rcutus, betas, phis, grid, &
                     nwfs, nbeta, bmat, qq, qvan, ikk, rhoc, rhos, &
                     vpsloc, ocs, rcloc
  use funct, only: get_iexch, get_icorr, get_igcx, get_igcc, dft_is_nonlocc
  implicit none
  !
  integer, intent(in):: iunps ! I/O unit
  !
  integer :: nb, mb, & ! counters on beta functions
             ios,    & ! I/O control
             ir        ! counter on mesh points
  integer :: iexch, icorr, igcx, igcc, inlc
  logical :: nonlocc
  !
  !
  nonlocc = dft_is_nonlocc()
  if (nonlocc) &
     CALL errore('write_rrkj','non-local functional not implemented yet', 1)

  write( iunps, '(a75)', err=100, iostat=ios ) title
  !
  write( iunps, '(i5)',err=100, iostat=ios ) pseudotype
  if (rel > 0) then
     write( iunps, '(2l5)',err=100, iostat=ios ) .true., nlcc
  else
     write( iunps, '(2l5)',err=100, iostat=ios ) .false., nlcc
  endif
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()
  inlc  = 0
  write( iunps, '(4i5)',err=100, iostat=ios ) iexch, icorr, igcx, igcc

  write( iunps, '(2e17.11,i5)') zval, etots, lmax
  write( iunps, '(4e17.11,i5)',err=100, iostat=ios ) &
       grid%xmin,grid%rmax,grid%zmesh,grid%dx,grid%mesh

  write( iunps, '(2i5)', err=100, iostat=ios ) nwfs, nbeta
  write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
       ( rcut(nb), nb=1,nwfs )
  write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
       ( rcutus(nb), nb=1,nwfs )
  do nb=1,nwfs
     write(iunps,'(a2,2i3,f6.2)',err=100,iostat=ios) &
          els(nb), nns(nb), lls(nb), ocs(nb)
  enddo
  do nb=1,nbeta
     write ( iunps, '(i6)',err=100, iostat=ios ) ikk(nb)
     write ( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
          ( betas(ir,nb), ir=1,ikk(nb))
     do mb=1,nb
        write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
             bmat(nb,mb)
        if (pseudotype == 3) then
           write(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                qq(nb,mb)
           write(iunps,'(1p4e19.11)',err=100,iostat=ios) & 
                (qvan(ir,nb,mb),ir=1,grid%mesh)
        endif
     enddo
  enddo
  !
  !   writes the local potential 
  !
  write( iunps, '(1p4e19.11)',err=100, iostat=ios ) rcloc, &
       ( vpsloc(ir), ir=1,grid%mesh )
  !
  !   writes the atomic charge
  !
  write( iunps, '(1p4e19.11)',err=100, iostat=ios )  &
       ( rhos(ir,1), ir=1,grid%mesh )
  !
  !   If present writes the core charge
  !
  if ( nlcc ) then 
     write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
          ( rhoc(ir), ir=1,grid%mesh )
  endif
  !
  !    Writes the wavefunctions of the atom
  !      
  write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
       ((phis(ir,nb),ir=1,grid%mesh),nb=1,nwfs)
100 call errore('write_rrkj','Writing pseudopw file',abs(ios))
  !
end subroutine write_rrkj

