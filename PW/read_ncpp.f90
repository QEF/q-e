!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine read_ncpp (np, iunps)
  !-----------------------------------------------------------------------
  !
  USE kinds, only: dp
  USE parameters, ONLY: nchix, lmaxx, ndm
  use atom,  only: zmesh, mesh, xmin, dx, r, rab, vnl, chi, oc, nchi, &
       lchi, rho_at, rho_atc, numeric
  use char, only: psd
  use nl_c_c,only: nlcc, a_nlcc, b_nlcc, alpha_nlcc
  use pseud, only: cc, alpc, zp, aps, alps, nlc, nnl, lmax, lloc, bhstype
  use funct, only: dft, which_dft
  implicit none

  integer :: iunps, np
  real(kind=DP) :: x
  integer :: nb, ios, i, l, ir
  !
  read (iunps, '(a)', end=300, err=300, iostat=ios) dft
  if (dft (1:2) .eq.'**') dft = 'PZ'
  call which_dft (dft)
  !
  read (iunps, *, err=300, iostat=ios) psd(np), zp(np), lmax(np), nlc(np), &
                                       nnl(np), nlcc(np), lloc(np), bhstype(np)
  if (nlc(np).gt.2 .or. nnl(np).gt.3) &
                            call errore ('read_ncpp', 'Wrong nlc or nnl', np)
  if (nlc(np)*nnl(np).lt.0) call errore ('read_ncpp', 'nlc*nnl < 0 ? ', np)
  if (zp(np).le.0d0) call errore ('read_ncpp', 'Wrong zp ', np)
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric (np) = nlc (np) .le.0.and.nnl (np) .le.0

  if (lloc (np).eq.-1000) lloc (np) = lmax (np)
  if (lloc (np).lt.0 .or. lmax(np).lt.0 .or. &
      .not.numeric(np) .and. (lloc(np).gt.min(lmax(np)+1,lmaxx+1) .or. &
      lmax(np).gt.max(lmaxx,lloc(np))) .or. &
      numeric(np) .and. (lloc(np).gt.lmax(np) .or. lmax(np).gt.lmaxx) ) &
                      call errore ('read_ncpp', 'wrong lmax and/or lloc', np)
  if (.not.numeric (np) ) then
     !
     !   read here pseudopotentials in analytic form
     !
     read (iunps, *, err=300, iostat=ios) (alpc(i,np), i=1,2), (cc(i,np), i=1,2)
     if (abs (cc(1,np)+cc(2,np)-1.d0) .gt. 1.0d-6) &
          call errore ('read_ncpp', 'wrong pseudopotential coefficients', 1)
     do l = 0, lmax (np)
        read (iunps, *, err=300, iostat=ios) (alps(i,l,np), i=1,3), &
                                             (aps(i,l,np),  i=1,6)
     enddo
     if (nlcc (np) ) then
        read (iunps, *, err=300, iostat=ios) a_nlcc(np), b_nlcc(np), &
                                             alpha_nlcc(np)
        if (alpha_nlcc(np).le.0.d0) call errore('read_ncpp','nlcc but alpha=0',np)
     endif
  endif
  read (iunps, *, err=300, iostat=ios) zmesh(np), xmin(np), dx(np), &
                                       mesh(np), nchi(np)
  if (mesh(np).gt.ndm .or. mesh(np).le.0) &
                         call errore ('read_ncpp', 'mesh too big', np)
  if ( nchi(np).gt.nchix .or. &
      (nchi(np).lt.lmax(np)   .and. lloc(np).eq.lmax(np)) .or. & 
      (nchi(np).lt.lmax(np)+1 .and. lloc(np).ne.lmax(np)) ) &
                         call errore ('read_ncpp', 'wrong no. of wfcts', np)
  !
  !  Here pseudopotentials in numeric form are read
  !
  if (numeric (np) ) then
     do l = 0, lmax (np)
        read (iunps, '(a)', err=300, iostat=ios)
        read (iunps, *, err=300, iostat=ios) (vnl(ir,l,np), ir=1,mesh(np) )
     enddo
     !
     ! and the local part is subtracted
     !
     do l = 0, lmax (np)
        if (l.ne.lloc (np) ) then
           do ir = 1, mesh (np)
              vnl(ir,l,np) = vnl(ir,l,np) - vnl(ir,lloc(np),np)
           enddo
        endif
     enddo
     if (nlcc (np) ) then
        !**            read( iunps, '(a)', err=300, iostat=ios )
        read (iunps, *, err=300, iostat=ios) (rho_atc(ir,np), ir=1,mesh(np))
     endif
  endif
  !
  !  Here pseudowavefunctions in numeric form are read
  !
  do nb = 1, nchi (np)
     read (iunps, '(a)', err=300, iostat=ios)
     read (iunps, *, err=300, iostat=ios) lchi(nb,np), oc(nb,np)
     !
     !     Test lchi and occupation numbers
     !
     if (nb.le.lmax(np) .and. lchi(nb,np)+1.ne.nb) &
                      call errore ('read_ncpp', 'order of wavefunctions', 1)
     if (lchi(nb,np).gt.lmaxx .or. lchi(nb,np).lt.0) &
                      call errore ('read_ncpp', 'wrong lchi', np)
     if (oc(nb,np).lt.0.d0 .or. oc(nb,np).gt.2.d0*(2*lchi(nb,np)+1)) &
                      call errore ('read_ncpp', 'wrong oc', np)
     read (iunps, *, err=300, iostat=ios) ( chi(ir,nb,np), ir=1,mesh(np) )
  enddo
  !
  !    compute the radial mesh
  !
  r (0, np) = 0.d0
  rab (0, np) = 0.d0
  do ir = 1, mesh (np)
     x = xmin (np) + dble (ir - 1) * dx (np)
     r (ir, np) = exp (x) / zmesh (np)
     rab (ir, np) = dx (np) * r (ir, np)
  enddo
  !
  !    compute the atomic charges
  !
  rho_at(:,np) = 0.d0
  do nb = 1, nchi (np)
     if (oc (nb, np) .ne.0.d0) then
        do ir = 1, mesh (np)
           rho_at(ir,np) = rho_at(ir,np) + oc(nb,np) * chi(ir,nb,np)**2
        enddo
     endif
  enddo

  return

300 call errore ('read_ncpp', 'pseudo file is empty or wrong', abs (np) )
end subroutine read_ncpp

