!     
!---------------------------------------------------------------------
subroutine ld1_writeout
  !---------------------------------------------------------------------
  !
  !     This routine writes on output the quantities which defines 
  !     a multiprojector pseudopotential. It can be in the
  !     Vanderbilt form or in the norm-conserving form
  !
  use ld1inc
  use funct
  use atomic_paw, only : paw_io
  implicit none

  integer :: &
       ios,   &  ! I/O control
       iunps     ! the unit with the pseudopotential

  logical, external :: matches
  logical :: oldformat
  
  if (file_pseudopw == ' ') return

  if (nconf > 1) &
       call errore('ld1_writeout','more than one test configuration',1)

  if (.not. lpaw) then

     if (rel == 2 .and. .not. matches('UPF',file_pseudopw) &
                  .and. .not. matches('upf',file_pseudopw) ) then
        file_pseudopw=trim(file_pseudopw)//'.UPF'
     end if

     oldformat = .not. matches('UPF',file_pseudopw) .and. &
                 .not. matches('upf',file_pseudopw)

  else
     if ( .not. matches('PAW',file_pseudopw) .and. &
          .not. matches('paw',file_pseudopw) ) then
        file_pseudopw=trim(file_pseudopw)//'.PAW'
     end if
  end if

  iunps=28
  open(unit=iunps, file=trim(file_pseudopw), status='unknown',  &
       form='formatted', err=50, iostat=ios)
50 call errore('ld1_writeout','opening file_pseudopw',abs(ios))

  if (lpaw) then
     !
     ! write experimental PAW setup
     !
     call paw_io(pawsetup,iunps,"OUT")
     !
  else if (oldformat) then
     !
     if (pseudotype == 1) then
       !
       ! write old "NC" format (semilocal)
       !
        call write_pseudo &
             (iunps,zed,xmin,dx,mesh,ndm,r,r2,  &
             dft,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
             rhoc,vnl,phis,vpsloc,elts,llts,octs,etots,nwfts)
     else
       !
       ! write old "RRKJ" format (nonlocal)
       !
        call write_rrkj ( iunps )
     end if
     !
  else
     !
     if (pseudotype == 1) then
        !
        !  prepare for writing UPF file 
        !
        if (rel == 2 ) then
           call copy_ncpp_so ()
        else
           call copy_ncpp ()
        end if
        !
     endif
     !
     call write_upf(iunps)
     !
  endif
  !
  close(iunps)
  !
  return
end subroutine ld1_writeout

!---------------------------------------------------------------------
subroutine write_rrkj (iunps)
  !---------------------------------------------------------------------
  !
  use ld1inc
  use funct
  implicit none
  !
  integer, intent(in):: iunps ! I/O unit
  !
  integer :: nb, mb, & ! counters on beta functions
             ios,    & ! I/O control
             ir        ! counter on mesh points
  !
  !
  write( iunps, '(a75)', err=100, iostat=ios ) title
  !
  write( iunps, '(i5)',err=100, iostat=ios ) pseudotype
  if (rel > 0) then
     write( iunps, '(2l5)',err=100, iostat=ios ) .true., nlcc
  else
     write( iunps, '(2l5)',err=100, iostat=ios ) .false., nlcc
  endif
  write( iunps, '(4i5)',err=100, iostat=ios ) iexch, icorr, &
       igcx, igcc

  write( iunps, '(2e17.11,i5)') zval, etots, lmax
  write( iunps, '(4e17.11,i5)',err=100, iostat=ios ) &
       xmin,rmax,zmesh,dx,mesh

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
                (qvan(ir,nb,mb),ir=1,mesh)
        endif
     enddo
  enddo
  !
  !   writes the local potential 
  !
  write( iunps, '(1p4e19.11)',err=100, iostat=ios ) rcloc, &
       ( vpsloc(ir), ir=1,mesh )
  !
  !   writes the atomic charge
  !
  write( iunps, '(1p4e19.11)',err=100, iostat=ios )  &
       ( rhos(ir,1), ir=1,mesh )
  !
  !   If present writes the core charge
  !
  if ( nlcc ) then 
     write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
          ( rhoc(ir), ir=1,mesh )
  endif
  !
  !    Writes the wavefunctions of the atom
  !      
  write( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
       ((phis(ir,nb),ir=1,mesh),nb=1,nwfs)
100 call errore('ld1_writeout','Writing pseudopw file',abs(ios))
  !
end subroutine write_rrkj

!---------------------------------------------------------------------
subroutine copy_ncpp_so ()
  !---------------------------------------------------------------------
  !
  use ld1inc
  implicit none
  !
  integer :: n,     & ! counter on wavefunctions
             nch,   & ! counter on chi functions
             l,     & ! counter on angular momentum
             ir       ! counter on mesh points
  real(kind=dp), external :: int_0_inf_dr
  real(kind=dp) :: aux(ndm)
  !
  if ( lloc == 0) then
     do ir=1,mesh
        vpsloc(ir)=vnlo(ir,lloc,1)
     enddo
  else if ( lloc > 0 ) then
     do ir=1,mesh
        vpsloc(ir) = ((lloc+1.0_dp)*vnlo(ir,lloc,2) + &
             lloc*vnlo(ir,lloc,1)) / (2.0_dp*lloc + 1.0_dp)
     enddo
  endif
  nbeta=0
  do l=0,lmax
     if (l /= lloc) then
        nbeta=nbeta+1
        nch=0
        do n=1,nwfts
           if (llts(n) == l .and. abs(jjts(n)-l+0.5_dp) < 1e-3_dp) nch=n
        enddo
        if (l==0) nch=1
        if (nch == 0) call errore('copy_ncpp_so','jj not found',1)
        do ir=1,mesh
           betas(ir,nbeta) = (vnlo(ir,l,1)-vpsloc(ir)) * phis(ir,nch) 
        enddo
        lls(nbeta)=llts(nch)
        jjs(nbeta)=jjts(nch)
        ikk(nbeta)=mesh
        do ir=mesh-1,1,-1
           if (abs(betas(ir,nbeta)) < 1.e-11_dp)then
              ikk(nbeta)=ir
           else
              goto 203
           endif
        enddo
203     continue
        do ir = 1, mesh
           aux (ir) = phis(ir, nch) * betas (ir, nbeta)
        enddo
        bmat(nbeta,nbeta)=1.0_dp/int_0_inf_dr(aux,r,r2,dx,mesh,2*(l+1))
        if (l /= 0) then
           nbeta=nbeta+1
           nch=0
           do n=1,nwfts
              if (llts(n) == l.and.abs(jjts(n)-l-0.5_dp) < 1e-3_dp) nch=n
           enddo
           if (nch == 0) call errore('convert','jj not found',1)
           do ir=1,mesh
              betas(ir,nbeta) = (vnlo(ir,l,2)-vpsloc(ir)) * phis(ir,nch) 
           enddo
           lls(nbeta)=llts(nch)
           jjs(nbeta)=jjts(nch)
           ikk(nbeta)=mesh
           do ir=mesh-1,1,-1
              if (abs(betas(ir,nbeta)) < 1.e-11_dp) then
                 ikk(nbeta)=ir
              else
                 goto 204
              endif
           enddo
204        continue
           do ir = 1, mesh
              aux(ir) = phis(ir,nch)*betas(ir,nbeta)
           enddo
           bmat(nbeta,nbeta)=1.0_dp/int_0_inf_dr(aux,r,r2,dx,mesh,2*(l+1))
        endif
     endif
  end do

  return
end subroutine copy_ncpp_so

!---------------------------------------------------------------------
subroutine copy_ncpp ()
  !---------------------------------------------------------------------
  !
  use ld1inc
  implicit none
  !
  integer :: n,     & ! counter on wavefunctions
             nch,   & ! counter on chi functions
             l,     & ! counter on angular momentum
             ir       ! counter on mesh points
  real(kind=dp), external :: int_0_inf_dr
  real(kind=dp) :: aux(ndm)
  !
  if ( lloc > -1) then
     do ir=1,mesh
        vpsloc(ir)=vnl (ir,lloc)
     enddo
  endif
  nbeta=0
  do l=0,lmax
     if (l /= lloc) then
        nch=0
        do n=1,nwfts
           if (llts(n) == l) nch=n
        enddo
        if (nch == 0) call errore('copy_ncpp','l not found',1)
        !
        nbeta=nbeta+1
        do ir=1,mesh
           betas(ir,nbeta) = (vnl(ir,l)-vpsloc(ir)) * phis(ir, nch)
        enddo
        lls(nbeta)=l
        ikk(nbeta)=mesh
        do ir=mesh-1,1,-1
           if (abs(betas(ir,nbeta)) < 1.e-11_dp)then
              ikk(nbeta)=ir
           else
              goto 203
           endif
        enddo
203     continue
        do ir = 1, mesh
           aux (ir) = phis(ir, nch) * betas (ir, nbeta)
        enddo
        bmat(nbeta,nbeta)=1.0_dp/int_0_inf_dr(aux,r,r2,dx,mesh,2*(l+1))
     endif
  end do

  return
end subroutine copy_ncpp
