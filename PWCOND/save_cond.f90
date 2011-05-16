!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine save_cond (lwrite, lsr, ef, nrz, nocros, noins,  &
                      norb, r, rab, betar)
!
!  This subroutine writes/reads variables needed for PWCOND
!  so that the punch file from PW calculations is not needed.
!
  use kinds, only : DP
  USE parameters, only : npsx
  use radial_grids, only: ndmx
  USE cell_base, ONLY : alat, tpiba, tpiba2, at, bg
  use lsda_mod, only: nspin
  USE noncollin_module, ONLY : noncolin, npol
  use spin_orb, only : lspinorb
  use cond, only : sarea, nrx, nry, norbf, tblml, crosl, taunewl, &
     zpseul, zpseul_nc, zl, vppotl, tblms, cross, taunews, zpseus,&
     zpseus_nc, zs, vppots, tblmr, crosr, taunewr, zpseur,        &
     zpseur_nc, zr, vppotr, iofspin, nbrx, save_file


  implicit none

  integer :: lsr, nrz, nocros, noins, norb, i, j, k, l, m
  logical :: lwrite
  REAL(DP) :: ef, r(ndmx,npsx), rab(ndmx,npsx),    &
                   betar(ndmx,nbrx,npsx)
  integer, ALLOCATABLE :: ind(:,:), tblm(:,:), cros(:,:)
  REAL(DP), ALLOCATABLE :: z(:), zpseu(:,:,:), re(:,:,:,:), &
                                im(:,:,:,:), c(:), taunew(:,:)
  COMPLEX(DP), ALLOCATABLE :: vppot(:,:,:,:), zpseu_nc(:,:,:,:)


  character(len=2) :: ext

  call start_clock('save_cond')

  if(lsr.eq.1) then
    ext='.l'
  elseif(lsr.eq.2) then
    ext='.s'
  elseif(lsr.eq.3) then
    ext='.r'
  endif

  if (lwrite) then
    allocate( vppot(nrz, nrx * nry, npol, npol) )
    allocate( z(nrz+1) )
    allocate( taunew(4,norb) )
    allocate( tblm(4,norb) )
    allocate( cros(norb, nrz) )
    if (noncolin) then
      allocate(zpseu_nc(2, norb, norb, nspin))
    else
      allocate( zpseu(2, norb, norb) )
    endif
    if(lsr.eq.1) then
      vppot = vppotl
      z = zl
      taunew = taunewl
      tblm = tblml
      cros = crosl
      if (noncolin) then
         zpseu_nc = zpseul_nc
      else
         zpseu = zpseul
      endif
    elseif(lsr.eq.2) then
      vppot = vppots
      z = zs
      taunew = taunews
      tblm = tblms
      cros = cross
      if (noncolin) then
         zpseu_nc = zpseus_nc
      else
         zpseu = zpseus
      endif
    elseif(lsr.eq.3) then
      vppot = vppotr
      z = zr
      taunew = taunewr
      tblm = tblmr
      cros = crosr
      if (noncolin) then
         zpseu_nc = zpseur_nc
      else
         zpseu = zpseur
      endif
    endif
    open (3,file=trim(save_file)//ext,form='formatted', &
                      status='unknown')
    write(3,*) nspin, npol, noncolin, lspinorb
    if(nspin.eq.2) write(3,*) iofspin
    write(3,*) alat, tpiba, tpiba2
    write(3,'(6f20.14)') ((at(i,j),i=1,3),j=1,3)
    write(3,'(6f20.14)') ((bg(i,j),i=1,3),j=1,3)
    write(3,*) sarea
    write(3,*) ef
    write(3,*) nrx, nry, nrz
    write(3,*) nocros, noins, norb, norbf
    write(3,'(40i3)') ((tblm(i,j),i=1,4),j=1,norb)
    write(3,'(120i1)') ((cros(j,i),i=1,nrz),j=1,norb)
    write(3,'(6f20.14)') ((taunew(i,j),i=1,4),j=1,norb)
!   write zpseu
    if(noncolin) then
      write(3,'(6f20.14)') (((( DBLE(zpseu_nc(i,j,k,l)),i=1,2),     &
            j=1,norb),k=1,norb),l=1,nspin)
      write(3,'(6f20.14)') ((((AIMAG(zpseu_nc(i,j,k,l)),i=1,2),     &
            j=1,norb),k=1,norb),l=1,nspin)
    else
      allocate( ind(3,2*norb*norb) )
      allocate( c(2*norb*norb) )
      m=0
      do i=1, 2
        do j=1, norb
          do k=1, norb
            if(abs(zpseu(i,j,k)).gt.1.d-12) then
              m = m+1
              ind(1,m) = i
              ind(2,m) = j
              ind(3,m) = k
              c(m) = zpseu(i,j,k)
            endif
          enddo
        enddo
      enddo
      write(3,*) m
      write(3,'(25i5)') ((ind(i,j),i=1,3),j=1,m)
      write(3,'(6f20.14)') (c(i),i=1,m)
      deallocate(ind)
      deallocate(c)
    endif
    write(3,'(6f20.14)') (z(i), i=1, nrz+1)
    write(3,'(6f20.14)') (((( DBLE(vppot(i,j,k,l)),i=1,nrz),        &
                             j=1,nrx*nry),k=1,npol),l=1,npol)
    write(3,'(6f20.14)') ((((AIMAG(vppot(i,j,k,l)),i=1,nrz),        &
                             j=1,nrx*nry),k=1,npol),l=1,npol)
!   write r
    allocate( ind(2,npsx*ndmx) )
    allocate( c(npsx*ndmx) )
    m=0
    do i=1, ndmx
      do j=1, npsx
          if(abs(r(i,j)).gt.1.d-12) then
            m = m+1
            ind(1,m) = i
            ind(2,m) = j
            c(m) = r(i,j)
          endif
      enddo
    enddo
    write(3,*) m
    write(3,'(25i5)') ((ind(i,j),i=1,2),j=1,m)
    write(3,'(6f20.14)') (c(i),i=1,m)
    deallocate(ind)
    deallocate(c)
!   write rab
    allocate( ind(2,npsx*ndmx) )
    allocate( c(npsx*ndmx) )
    m=0
    do i=1, ndmx
      do j=1, npsx
          if(abs(rab(i,j)).gt.1.d-12) then
            m = m+1
            ind(1,m) = i
            ind(2,m) = j
            c(m) = rab(i,j)
          endif
      enddo
    enddo
    write(3,*) m
    write(3,'(25i5)') ((ind(i,j),i=1,2),j=1,m)
    write(3,'(6f20.14)') (c(i),i=1,m)
    deallocate(ind)
    deallocate(c)
!   write betar
    allocate( ind(3,npsx*nbrx*ndmx) )
    allocate( c(npsx*nbrx*ndmx) )
    m=0
    do i=1, ndmx
      do j=1, nbrx
        do k=1, npsx
          if(abs(betar(i,j,k)).gt.1.d-12) then
            m = m+1
            ind(1,m) = i
            ind(2,m) = j
            ind(3,m) = k
            c(m) = betar(i,j,k)
          endif
        enddo
      enddo
    enddo
    write(3,*) m
    write(3,'(25i5)') ((ind(i,j),i=1,3),j=1,m)
    write(3,'(6f20.14)') (c(i),i=1,m)
    deallocate(ind)
    deallocate(c)
    close(unit=3)
  else
    open (3,file=trim(save_file)//ext,form='formatted', &
                      status='unknown')
    read(3,*) nspin, npol, noncolin, lspinorb
    if(nspin.eq.2) read(3,*) iofspin
    read(3,*) alat, tpiba, tpiba2
    read(3,'(6f20.14)') ((at(i,j),i=1,3),j=1,3)
    read(3,'(6f20.14)') ((bg(i,j),i=1,3),j=1,3)
    read(3,*) sarea
    read(3,*) ef
    read(3,*) nrx, nry, nrz
    read(3,*) nocros, noins, norb, norbf
    allocate( vppot(nrz, nrx * nry, npol, npol) )
    allocate( z(nrz+1) )
    allocate( taunew(4,norb) )
    allocate( tblm(4,norb) )
    allocate( cros(norb, nrz) )
    if (noncolin) then
       allocate(zpseu_nc(2, norb, norb, nspin))
    else
       allocate( zpseu(2, norb, norb) )
    endif
    read(3,'(40i3)') ((tblm(i,j),i=1,4),j=1,norb)
    read(3,'(120i1)') ((cros(j,i),i=1,nrz),j=1,norb)
    read(3,'(6f20.14)') ((taunew(i,j),i=1,4),j=1,norb)
!   read zpseu
    if(noncolin) then
      allocate ( re(2,norb,norb,nspin) )
      allocate ( im(2,norb,norb,nspin) )
      read(3,'(6f20.14)') ((((re(i,j,k,l),i=1,2),     &
            j=1,norb),k=1,norb),l=1,nspin)
      read(3,'(6f20.14)') ((((im(i,j,k,l),i=1,2),     &
            j=1,norb),k=1,norb),l=1,nspin)
      zpseu_nc = CMPLX(re,im,kind=DP)
      deallocate(re)
      deallocate(im)
    else
      read(3,*) m
      allocate( ind(3,m) )
      allocate( c(m) )
      read(3,'(25i5)') ((ind(i,j),i=1,3),j=1,m)
      read(3,'(6f20.14)') (c(i),i=1,m)
      zpseu = 0.d0
      do i=1, m
        zpseu(ind(1,i),ind(2,i),ind(3,i)) = c(i)
      enddo
      deallocate(ind)
      deallocate(c)
    endif
!-------------
    read(3,'(6f20.14)') (z(i), i=1, nrz+1)
    allocate ( re(nrz,nrx*nry,npol,npol) )
    allocate ( im(nrz,nrx*nry,npol,npol) )
    read(3,'(6f20.14)') ((((re(i,j,k,l),i=1,nrz),     &
          j=1,nrx*nry),k=1,npol),l=1,npol)
    read(3,'(6f20.14)') ((((im(i,j,k,l),i=1,nrz),     &
          j=1,nrx*nry),k=1,npol),l=1,npol)
    vppot = CMPLX(re,im,kind=DP)
    deallocate(re)
    deallocate(im)
!   read r
    read(3,*) m
    allocate( ind(2,m) )
    allocate( c(m) )
    read(3,'(25i5)') ((ind(i,j),i=1,2),j=1,m)
    read(3,'(6f20.14)') (c(i),i=1,m)
    r = 0.d0
    do i=1,m
     r(ind(1,i),ind(2,i))=c(i)
    enddo
    deallocate(ind)
    deallocate(c)
!   read rab
    read(3,*) m
    allocate( ind(2,m) )
    allocate( c(m) )
    read(3,'(25i5)') ((ind(i,j),i=1,2),j=1,m)
    read(3,'(6f20.14)') (c(i),i=1,m)
    rab = 0.d0
    do i=1,m
     rab(ind(1,i),ind(2,i))=c(i)
    enddo
    deallocate(ind)
    deallocate(c)
!   read betar
    read(3,*) m
    allocate( ind(3,m) )
    allocate( c(m) )
    read(3,'(25i5)') ((ind(i,j),i=1,3),j=1,m)
    read(3,'(6f20.14)') (c(i),i=1,m)
    betar = 0.d0
    do i=1,m
     betar(ind(1,i),ind(2,i),ind(3,i))=c(i)
    enddo
    deallocate(ind)
    deallocate(c)
    close(unit=3)
    if(lsr.eq.1) then
      allocate( vppotl(nrz, nrx * nry, npol, npol) )
      allocate( zl(nrz+1) )
      allocate( taunewl(4,norb) )
      allocate( tblml(4,norb) )
      allocate( crosl(norb, nrz) )
      if (noncolin) then
        allocate(zpseul_nc(2, norb, norb, nspin))
      else
        allocate( zpseul(2, norb, norb) )
      endif
      vppotl = vppot
      zl = z
      taunewl = taunew
      tblml = tblm
      crosl = cros
      if (noncolin) then
        zpseul_nc = zpseu_nc
      else
        zpseul = zpseu
      endif
    elseif(lsr.eq.2) then
      allocate( vppots(nrz, nrx * nry, npol, npol) )
      allocate( zs(nrz+1) )
      allocate( taunews(4,norb) )
      allocate( tblms(4,norb) )
      allocate( cross(norb, nrz) )
      if (noncolin) then
        allocate(zpseus_nc(2, norb, norb, nspin))
      else
        allocate( zpseus(2, norb, norb) )
      endif
      vppots = vppot
      zs = z
      taunews = taunew
      tblms = tblm
      cross = cros
      if (noncolin) then
        zpseus_nc = zpseu_nc
      else
        zpseus = zpseu
      endif
    elseif(lsr.eq.3) then
      allocate( vppotr(nrz, nrx * nry, npol, npol) )
      allocate( zr(nrz+1) )
      allocate( taunewr(4,norb) )
      allocate( tblmr(4,norb) )
      allocate( crosr(norb, nrz) )
      if (noncolin) then
        allocate(zpseur_nc(2, norb, norb, nspin))
      else
        allocate( zpseur(2, norb, norb) )
      endif
      vppotr = vppot
      zr = z
      taunewr = taunew
      tblmr = tblm
      crosr = cros
      if (noncolin) then
        zpseur_nc = zpseu_nc
      else
        zpseur = zpseu
      endif
    endif
  endif

    deallocate( vppot )
    deallocate( z )
    deallocate( taunew )
    deallocate( tblm )
    deallocate( cros )
    if (noncolin) then
      deallocate(zpseu_nc)
    else
      deallocate( zpseu )
    endif

    call stop_clock('save_cond')

end subroutine save_cond



