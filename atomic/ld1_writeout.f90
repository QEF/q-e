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
  implicit none

  integer :: &
       nb,mb, &  ! counters on beta functions
       ir,    &  ! counters on mesh points
       ios,   &  ! I/O control
       n,     &  ! counter on mesh points
       nch,   &  ! counter on chi functions
       l,     &  ! counter on angular momenta
       iunps     ! the unit with the pseudopotential

  real(kind=dp) :: int_0_inf_dr, aux(ndm)

  logical :: matches

  
  if (file_pseudopw.eq.' ') return

  if (nconf.gt.1) &
       call errore('ld1_writeout','more than one test configuration',1)

  if (rel.eq.2.and..not.matches('UPF',file_pseudopw)) &
       file_pseudopw=trim(file_pseudopw)//'.UPF'

  if (matches('UPF',file_pseudopw).or.matches('upf',file_pseudopw)) goto 201

  iunps=28
  open(unit=iunps, file=file_pseudopw, status='unknown',  &
       form='formatted', err=50, iostat=ios)
50 call errore('ld1_writeout','opening file_pseudopw',abs(ios))

  if (pseudotype.eq.1) then
     call write_pseudo &
          (iunps,zed,xmin,dx,mesh,ndm,r,r2,  &
          dft,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
          rhoc,vnl,phis,vpsloc,llts,octs,etots,nwfts)

     close(iunps)
     return
  endif

  write( iunps, '(a75)', err=100, iostat=ios ) title

  write( iunps, '(i5)',err=100, iostat=ios ) pseudotype
  if (rel.gt.0) then
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
        if (pseudotype.eq.3) then
           write(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                qq(nb,mb)
           write(iunps,'(1p4e19.11)',err=100,iostat=ios) & 
                (qvan(n,nb,mb),n=1,mesh)
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
  close(iunps)

201 continue


  if (matches('UPF',file_pseudopw).or.matches('upf',file_pseudopw)) then
     if (rel.eq.2.and.pseudotype.eq.1) then
        !
        !  prepare for writing UPF file 
        !
        if (iswitch.eq.2) then 
           if (lloc.ne.-1) then
              if (lloc==0) then
                 do ir=1,mesh
                    vpsloc(ir)=vnlo(ir,lloc,1)
                 enddo
              else
                 do ir=1,mesh
                    vpsloc(ir)=((lloc+1.0_dp)*vnlo(ir,lloc,2)+lloc*vnlo(ir,lloc,1))/&
                         (2.0_dp*lloc + 1.0_dp)
                 enddo
              endif
           endif
           nbeta=0
           do l=0,lmax
              if (l.ne.lloc) then
                 nbeta=nbeta+1
                 nch=0
                 do n=1,nwfts
                    if (llts(n).eq.l.and.abs(jjts(n)-l+0.5_dp).lt.1e-3_dp) nch=n
                 enddo
                 if (l==0) nch=1
                 if (nch.eq.0) call errore('convert','jj not found',1)
                 do ir=1,mesh
                    betas(ir,nbeta) = (vnlo(ir,l,1)-vpsloc(ir)) * phis(ir,nch) 
                 enddo
                 lls(nbeta)=llts(nch)
                 jjs(nbeta)=jjts(nch)
                 ikk(nbeta)=mesh
                 do ir=mesh-1,1,-1
                    if (abs(betas(ir,nbeta)).lt.1.e-11_dp)then
                       ikk(nbeta)=ir
                    else
                       goto 203
                    endif
                 enddo
203              continue
                 do ir = 1, mesh
                    aux (ir) = phis(ir, nch) * betas (ir, nbeta)
                 enddo
                 bmat(nbeta,nbeta)=1.0_dp/int_0_inf_dr(aux,r,r2,dx,mesh,2*(l+1))
                 if (l.ne.0) then
                    nbeta=nbeta+1
                    nch=0
                    do n=1,nwfts
                       if (llts(n).eq.l.and.abs(jjts(n)-l-0.5_dp).lt.1e-3_dp) nch=n
                    enddo
                    if (nch.eq.0) call errore('convert','jj not found',1)
                    do ir=1,mesh
                       betas(ir,nbeta) = (vnlo(ir,l,2)-vpsloc(ir)) * phis(ir,nch) 
                    enddo
                    lls(nbeta)=llts(nch)
                    jjs(nbeta)=jjts(nch)
                    ikk(nbeta)=mesh
                    do ir=mesh-1,1,-1
                       if (abs(betas(ir,nbeta)).lt.1.e-11_dp) then
                          ikk(nbeta)=ir
                       else
                          goto 204
                       endif
                    enddo
204                 continue
                    do ir = 1, mesh
                       aux(ir) = phis(ir,nch)*betas(ir,nbeta)
                    enddo
                    bmat(nbeta,nbeta)=1.0_dp/int_0_inf_dr(aux,r,r2,dx,mesh,2*(l+1))
                 endif
              endif
           enddo
        endif
     endif

     iunps=28
     open(unit=iunps, file=trim(file_pseudopw), status='unknown',  &
          form='formatted', err=51, iostat=ios)
51   call errore('ld1_writeout','opening file_pseudopw',abs(ios))

     call write_upf(iunps)
     close(iunps)
  endif
  return
end subroutine ld1_writeout
