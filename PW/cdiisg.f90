!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cdiisg (ndim, ndmx, nvec, nvecx, buflen, btype, psi, &
     ethr, e, notcnv, iter, keep_flag)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * psi = 0
  !
  !     where h is a complex hermitean matrix, e is a real scalar,
  !     and S is a complex hermitean matrix and psi is a complex vector
  !
  !     Gabriele Cipriani 3/00
  !
#include "machine.h"
  use parameters, only : DP
  implicit none
  logical :: keep_flag
  ! if to keep the old wfc in the subspac

  integer :: ndim, ndmx, nvec, nvecx, buflen
  ! dimension of the matrix to be
  ! diagonalized
  ! leading dimension of matrix psi,
  ! as declared in the calling pgm unit
  !      number of sought eigeipairs
  ! max  ""
  ! wfc buffer lenght
  complex(kind=DP) :: psi (ndmx, nvecx)
  ! wfc
  real(kind=DP) :: ethr
  ! energy threshold for convergence.
  ! root improvement is stopped, when two
  ! consecutive estimates of the root dif
  ! by less than ethr.
  integer :: btype (nvecx)
  ! on OUTPUT
  !     complex*16 psi             ! the first nvec columns contain the
  !                                ! refined estimates of the eigenvectors
  real(kind=DP) :: e (nvec)
  real(kind=DP), allocatable ::  et (:)
  ! energies
  ! idem

  integer :: iter, notcnv
  ! number of iterations performed.
  ! number of unconverged (valence) roots
  ! LOCAL variables
  !
  ! parameters
  !
  integer :: maxter, diis_steps, cbnd_steps
  ! maximum number of iterations
  ! maximum number  diis refinements
  ! maximum number  diis refinements for co
  real(kind=DP) :: pi, min_tstep, max_tstep, small, big, shift, cfact
  !
  ! minimum diis step lenght
  ! maximum diis step lenght
  ! diagonal shift for rmat
  ! convergence factor for the norm of the

  parameter (maxter = 1, diis_steps = 4, cbnd_steps = 2, pi = &
       3.14159265359d0, min_tstep = 0.1d0, max_tstep = 1.d0, small = &
       1.d-10, big = 1.d10, shift = 1.d-12, cfact = 0.01d0)
  integer :: kter
  integer, allocatable  :: conv (:)
  integer ::nc, n, nb, ns, wfoff, &
       bufoff, nopt, bufwfc, i, j, ib, nvectot, nbuf, ibuf, jbuf, ipos
  integer, allocatable ::  buf2wfc (:), wfc2buf (:,:)
  ! counter on iterations
  ! conv flag
  ! bnd type
  ! counter on wfc buffers
  ! counter on wfc
  ! counter on wfc in the buffer
  ! counter on diis steps
  ! first wfc of a given buffer
  ! buffer offset
  ! number of wfc still to optimize
  ! number of wfc in the buffer
  ! dimension of the iterative subspace
  ! indexing table for the buffe
  ! idem
  complex(kind=DP), allocatable :: vc (:,:), vcc (:,:),&
       rmat (:,:,:), smat (:,:,:), psip (:,:), psis (:,:), &
       pres (:,:), dwfc (:,:), dres (:,:), eau(:), cstep(:)
  complex(kind=DP) ::  ZDOTC, hh (2, 2), ss (2, 2), hv (2, 2)
  !
  !
  ! residues matrix
  ! overlap  matrix
  ! the product of H and
  ! the product of S and
  ! preconditioned residu
  ! iterative space
  ! residues
  ! scalar product routin
  ! auxiliary complex var
  !  diis trial step leng

  real(kind=DP), allocatable :: nrm (:), ew (:), ep (:), np (:), ec (:)
  real(kind=DP) :: h (2), tstep

  external ZDOTC
  !
  ! allocate the work arrays
  !

  call start_clock ('cdiisg')
  allocate (vc(  diis_steps , diis_steps))    
  allocate (vcc( diis_steps , buflen))    
  allocate (psip(  ndmx , buflen))    
  allocate (psis(  ndmx , buflen))    
  allocate (pres(  ndmx , buflen))    
  allocate (dwfc(  ndmx , diis_steps * buflen))    
  allocate (dres(  ndmx , diis_steps * buflen))    
  allocate (ec( diis_steps))    
  allocate (ep( nvec))    
  allocate (ew( buflen))    
  allocate (et( 2 * nvec))    
  allocate (np( nvec))    
  allocate (eau( buflen))    
  allocate (nrm( buflen))    
  allocate (rmat(  diis_steps,  diis_steps,  buflen))    
  allocate (smat(  diis_steps,  diis_steps,  buflen))    
  allocate (conv( nvec))    
  allocate (cstep( buflen))    
  allocate (buf2wfc(  buflen * diis_steps))    
  allocate (wfc2buf(  buflen , diis_steps))    

  if (diis_steps.lt.cbnd_steps) call errore ('cdiisg', 'wrong diis_steps', 1)
  !
  do n = 1, nvec
     conv (n) = 0
     ep (n) = 0.d0
     np (n) = 0.d0
  enddo
  do n = 1, buflen
     ew (n) = 0.d0
  enddo
  do n = 1, 2 * nvec
     et (n) = 0.d0
  enddo
  do j = 1, diis_steps
     do i = 1, buflen
        buf2wfc (i * j) = 0
        wfc2buf (i, j) = 0
     enddo
  enddo
  ! number of wfc buffers
  nbuf = nvec / buflen

  if (mod (nvec, buflen) .ne.0) nbuf = nbuf + 1
  ! DIIS

  do kter = 1, maxter
     iter = kter

     if (kter.eq.1) nvectot = nvec

     call rotate_wfc (ndim, ndmx, nvectot, nvectot, et, psi)
     !        write(*,'(4f12.6)')(et(n),n=1,nvec)
     ! save bands
     if (keep_flag) then
        call ZCOPY (ndmx * nvec, psi (1, 1), 1, psi (1, nvec + 1), &
             1)
        nvectot = 2 * nvec
     else
        nvectot = nvec
     endif
     ! buffers' loop

     do nc = 1, nbuf
        wfoff = (nc - 1) * buflen + 1
        bufwfc = min (nvec - wfoff + 1, buflen)
        ! starting indexing
        nopt = bufwfc
        bufoff = 1

        ipos = bufwfc
        do ib = 1, bufwfc
           buf2wfc (ib) = ib
           wfc2buf (ib, 1) = ib
        enddo
        ! workspaces:
        call setv (2 * diis_steps * diis_steps * buflen, 0.d0, rmat, 1)
        call setv (2 * diis_steps * diis_steps * buflen, 0.d0, smat, 1)
        call setv (2 * ndmx * diis_steps * buflen, 0.d0, dwfc, 1)
        call setv (2 * ndmx * diis_steps * buflen, 0.d0, dres, 1)
        call setv (2 * ndmx * buflen, 0.d0, psip, 1)
        call setv (2 * ndmx * buflen, 0.d0, psis, 1)

        call ZCOPY (ndmx * nopt, psi (1, wfoff), 1, dwfc (1, 1), 1)
        ! diis optimization:

        do ns = 1, diis_steps
           ! residues
           call h_psi (ndmx, ndim, nopt, dwfc (1, bufoff), psip, psis)
           !
           do ib = 1, nopt
              ibuf = bufoff + ib - 1
              nrm (ib) = 0.d0
              nrm (ib) = DREAL (ZDOTC (ndim, dwfc (1,ibuf),1,psis (1,ib),1) )
           enddo
#ifdef __PARA
           call reduce (nopt, nrm)
#endif
           ! rescale
           do ib = 1, nopt
              ibuf = bufoff + ib - 1
              nrm (ib) = 1.d0 / dsqrt (nrm (ib) )
              call DSCAL (2 * ndmx, nrm (ib), dwfc (1, ibuf), 1)
              call DSCAL (2 * ndmx, nrm (ib), psip (1, ib), 1)
              call DSCAL (2 * ndmx, nrm (ib), psis (1, ib), 1)
              eau (ib) = - ZDOTC (ndim, dwfc (1, ibuf), 1, psip (1, ib),1)
           enddo
#ifdef __PARA
           call reduce (2 * nopt, eau)
#endif
           ! energies
           do ib = 1, nopt
              ew (ib) = - DREAL (eau (ib) )
           enddo
           ! residues
           call ZCOPY (ndmx * nopt, psip, 1, dres (1, bufoff), 1)
           do ib = 1, nopt
              ibuf = bufoff + ib - 1
              call ZAXPY (ndmx, eau (ib), psis (1,ib),1,dres (1,ibuf),1)
           enddo
           ! update rmat/smat
           do ib = 1, nopt
              ibuf = bufoff + ib - 1
              nb = buf2wfc (ibuf)
              do j = 1, ns
                 jbuf = wfc2buf (nb, j)
                 rmat (j, ns, nb) = ZDOTC (ndim,dres(1,jbuf),1,dres(1,ibuf),1)
#ifdef __PARA
                 call reduce (2, rmat (j, ns, nb) )
#endif
                 rmat (ns, j, nb) = conjg (rmat (j, ns, nb) )
                 if (ns.eq.1) then
                    smat (ns, ns, nb) = (1.d0, 0.d0)
                 else
                    smat (ns, j, nb) = ZDOTC(ndim,psis(1,ib),1,dwfc(1,jbuf),1)
#ifdef __PARA
                    call reduce (2, smat (ns, j, nb) )
#endif
                    smat (j, ns, nb) = conjg (smat (ns, j, nb) )
                 endif
              enddo
           enddo
           ! diagonal shift
           do ib = 1, nopt
              ibuf = bufoff + ib - 1
              nb = buf2wfc (ibuf)
              rmat (ns, ns, nb) = rmat (ns, ns, nb) + shift
              smat (ns, ns, nb) = smat (ns, ns, nb) + shift
           enddo
           ! dump
           !          write(6,*)
           !          do ib=1,nopt
           !             ibuf = bufoff+ib-1
           !             nb = buf2wfc(ibuf)
           !             n = wfoff+nb-1
           !             write(6,"('n ns rmat ew de ',2i5,3e14.6)")n,ns,
           !     +            DREAL(rmat(ns,ns,nb)),ew(ib),
           !     +            dabs(ew(ib)-ep(n))
           !          enddo
           !          write(6,*)

           if (ns.eq.1) then
              ! first iteration
              ! ep np
              do ib = 1, nopt
                 ibuf = bufoff + ib - 1
                 nb = buf2wfc (ibuf)
                 n = wfoff + nb - 1
                 ep (n) = ew (ib)
                 np (n) = rmat (ns, ns, nb)
              enddo
              ! preconditioned residue
              call ZCOPY (ndmx * nopt, dres (1, 1), 1, pres, 1)
              call g_psi (ndmx, ndim, nopt, pres, ew)
              ! minimize rayleigh quotient
              call h_psi (ndmx, ndim, nopt, pres, psip, psis)
              if (.true.) then
                 do ib = 1, nopt
                    ibuf = bufoff + ib - 1
                    nb = buf2wfc (ibuf)

                    n = wfoff + nb - 1
                    call setv (8, 0.d0, hh, 1)
                    call setv (8, 0.d0, ss, 1)
                    hh (1, 1) = DCMPLX (ew (ib), 0.d0)
                    hh (1, 2) = ZDOTC (ndim, dwfc (1, ib), 1, psip (1, ib),1)
                    hh (2, 1) = conjg (hh (1, 2) )
                    hh (2, 2) = ZDOTC (ndim, pres (1, ib), 1, psip (1, ib),1)
                    ss (1, 1) = smat (1, 1, ib)
                    ss (1, 2) = ZDOTC (ndim, dwfc (1, ib), 1, psis (1, ib),1)
                    ss (2, 1) = conjg (ss (1, 2) )
                    ss (2, 2) = ZDOTC (ndim, pres (1, ib), 1, psis (1, ib),1)
#ifdef __PARA
                    call reduce (8, hh)
                    call reduce (8, ss)
                    hh (1, 1) = DCMPLX (ew (ib), 0.d0)
                    ss (1, 1) = smat (1, 1, ib)
#endif

                    call cdiaghg (2, 2, hh, ss, 2, h, hv)
                    if (abs (hv (1, 1) ) .lt.small) then
                       cstep (ib) = DCMPLX(min_tstep, 0.d0)
                    else
                       cstep (ib) = hv (2, 1) / hv (1, 1)
                       !                          write(6,*)'cstep0 ',n,cstep(ib)
                       tstep = abs (cstep (ib) )
                       cstep (ib) = cstep (ib) / tstep
                       tstep = max (tstep, min_tstep)
                       tstep = min (tstep, max_tstep)
                       cstep (ib) = cstep (ib) * tstep
                    endif
                    !                       write(6,*)'cstep1 ',n,cstep(ib)
                 enddo
              else
                 do ib = 1, nopt
                    cstep (ib) = ( - 1.d0, 0.d0)
                 enddo
              endif
              ! next dwfc
              if (diis_steps.eq.1) then
                 do ib = 1, nopt
                    ibuf = bufoff + ib - 1
                    nb = buf2wfc (ibuf)
                    n = wfoff + nb - 1
                    call ZCOPY (ndmx, dwfc (1, ib), 1, psi (1, n), 1)
                    call ZAXPY (ndmx, cstep (ib), pres(1,ib),1,psi (1, n),1)
                 enddo
              else
                 do ib = 1, nopt
                    ibuf = bufoff + ib - 1
                    nb = buf2wfc (ibuf)
                    ipos = ipos + 1
                    call ZCOPY (ndmx, dwfc (1, ib), 1, dwfc (1, ipos), 1)
                    call ZAXPY (ndmx,cstep(ib),pres(1,ib),1,dwfc(1,ipos),1)
                    wfc2buf (nb, ns + 1) = ipos
                    buf2wfc (ipos) = nb
                 enddo
              endif
              !
              bufoff = bufoff + nopt

              nopt = ipos - bufoff + 1
              ! ns > 1
           else
              ! check convergence
              do ib = 1, nopt
                 ibuf = bufoff + ib - 1
                 nb = buf2wfc (ibuf)
                 n = wfoff + nb - 1
                 if (dabs (ew (ib) - ep (n) ) .lt.ethr .or. &
                     DREAL (rmat (ns, ns, nb) ) .lt.cfact * np (n) .or. &
                     ((btype (n).eq.1) .and. (ns.eq.cbnd_steps)) ) conv(n)=1
                 ep (n) = ew (ib)
              enddo
              ! solve diis eigenproblem
              do ib = 1, nopt
                 ibuf = bufoff + ib - 1
                 nb = buf2wfc (ibuf)
                 call cdiaghg (ns, ns, rmat (1, 1, nb), smat (1, 1, nb), &
                      diis_steps, ec, vc)
                 if (ec (1) .lt.0.d0.and. (dabs (ec (1) ) .gt.small) ) then
                    write ( * , * ) ' *********** ec *********** '
                    write ( *, * ) (ec (i), i = 1, diis_steps)
                    call errore ('cdiisg', 'ec(1) < 0', 1)
                 endif
                 do j = 1, ns
                    vcc (j, nb) = vc (j, 1)
                 enddo
              enddo
              ! new residue
              call setv (2 * ndmx * nopt, 0.d0, pres, 1)
              do ib = 1, nopt
                 ibuf = bufoff + ib - 1
                 nb = buf2wfc (ibuf)
                 do j = 1, ns
                    jbuf = wfc2buf (nb, j)
                    call ZAXPY (ndmx,vcc(j,nb),dres(1,jbuf),1,pres(1,ib),1)
                 enddo
              enddo
              !
              call g_psi (ndmx, ndim, nopt, pres, ew)
              ! new wfc
              do ib = 1, nopt
                 ibuf = bufoff + ib - 1
                 nb = buf2wfc (ibuf)
                 n = wfoff + nb - 1
                 if (conv (n) .eq.1.or.ns.eq.diis_steps) then
                    call setv (2 * ndmx, 0.d0, psi (1, n), 1)
                    do j = 1, ns
                       jbuf = wfc2buf (nb, j)
                       call ZAXPY (ndmx,vcc(j,nb),dwfc(1,jbuf),1,psi(1,n),1)
                    enddo
                    call ZAXPY (ndmx,cstep(ib),pres(1,ib),1,psi(1,n),1)
                 else
                    ipos = ipos + 1
                    do j = 1, ns
                       jbuf = wfc2buf (nb, j)
                       call ZAXPY (ndmx, vcc(j,nb), dwfc(1,jbuf), 1, &
                            dwfc(1,ipos),1)
                    enddo
                    call ZAXPY (ndmx, cstep (ib), pres (1, ib), 1, &
                         dwfc (1,ipos), 1)
                    wfc2buf (nb, ns + 1) = ipos
                    buf2wfc (ipos) = nb
                 endif
              enddo
              !
              bufoff = bufoff + nopt
              !
              nopt = ipos - bufoff + 1
              if (nopt.eq.0) goto 10
           endif
           ! ns loop
        enddo
        !
10      continue
        ! nc loop
     enddo
     ! update convergence flags
     notcnv = 0
     do n = 1, nvec
        if (btype (n) .eq.0.and.conv (n) .ne.1) notcnv = notcnv + 1
     enddo
     ! exit?
     if (notcnv.eq.0) goto 100
     ! kter loop
  enddo

100 continue
  !
  call rotate_wfc (ndim, ndmx, nvectot, nvectot, et, psi)
  !
  do n = 1, nvec
     e (n) = et (n)
  enddo

  deallocate (wfc2buf)
  deallocate (buf2wfc)
  deallocate (cstep)
  deallocate (conv)
  deallocate (smat)
  deallocate (rmat)
  deallocate (nrm)
  deallocate (eau)
  deallocate (np)
  deallocate (et)
  deallocate (ew)
  deallocate (ep)
  deallocate (ec)
  deallocate (dres)
  deallocate (dwfc)
  deallocate (pres)
  deallocate (psis)
  deallocate (psip)
  deallocate (vcc)
  deallocate (vc)

  call stop_clock ('cdiisg')
  return
end subroutine cdiisg

