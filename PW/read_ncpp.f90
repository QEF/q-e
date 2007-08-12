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
  USE parameters, ONLY: nchix, lmaxx
  use radial_grids, only: ndmx
  use atom,  only: rgrid, chi, oc, &
       nchi, lchi, rho_at, rho_atc, numeric, nlcc
  use pseud, only: cc, alpc, aps, alps, nlc, nnl, lmax, lloc, &
       a_nlcc, b_nlcc, alpha_nlcc
  use uspp_param, only: zp, vloc_at, betar, kkbeta, nbeta, lll, dion, psd
  use funct, only: set_dft_from_name, dft_is_meta, dft_is_hybrid
  implicit none
  !
  integer :: iunps, np
  !
  real(DP) :: x, vll
  real(DP), allocatable:: vnl(:,:)
  real(DP), parameter :: rcut = 10.d0, e2 = 2.d0
  real(DP), external :: erf
  integer :: nb, ios, i, l, ir
  logical :: bhstype
  character (len=20) :: dft_name
  !
  !====================================================================
  ! read norm-conserving PPs
  !
  read (iunps, '(a)', end=300, err=300, iostat=ios) dft_name
  if (dft_name (1:2) .eq.'**') dft_name = 'PZ'
  read (iunps, *, err=300, iostat=ios) psd(np), zp(np), lmax(np), nlc(np), &
                                       nnl(np), nlcc(np), lloc(np), bhstype
  if (nlc(np) > 2 .or. nnl(np) > 3) &
       call errore ('read_ncpp', 'Wrong nlc or nnl', np)
  if (nlc(np)*nnl(np) < 0) call errore ('read_ncpp', 'nlc*nnl < 0 ? ', np)
  if (zp(np) <= 0d0 .or. zp(np) > 100 ) &
       call errore ('read_ncpp', 'Wrong zp ', np)
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric (np) = (nlc (np) <= 0) .and. (nnl (np) <= 0)
  !
  if (lloc (np) == -1000) lloc (np) = lmax (np)
  if (lloc (np) < 0 .or. lmax(np) < 0 .or. &
       .not.numeric(np) .and. (lloc(np) > min(lmax(np)+1,lmaxx+1) .or. &
       lmax(np) > max(lmaxx,lloc(np))) .or. &
       numeric(np) .and. (lloc(np) > lmax(np) .or. lmax(np) > lmaxx) ) &
       call errore ('read_ncpp', 'wrong lmax and/or lloc', np)
  if (.not.numeric (np) ) then
     !
     !   read here pseudopotentials in analytic form
     !
     read (iunps, *, err=300, iostat=ios) &
          (alpc(i,np), i=1,2), (cc(i,np), i=1,2)
     if (abs (cc(1,np)+cc(2,np)-1.d0) > 1.0d-6) &
          call errore ('read_ncpp', 'wrong pseudopotential coefficients', 1)
     do l = 0, lmax (np)
        read (iunps, *, err=300, iostat=ios) (alps(i,l,np), i=1,3), &
                                             (aps(i,l,np),  i=1,6)
     enddo
     if (nlcc (np) ) then
        read (iunps, *, err=300, iostat=ios) &
             a_nlcc(np), b_nlcc(np), alpha_nlcc(np)
        if (alpha_nlcc(np) <= 0.d0) call errore('read_ncpp','alpha_nlcc=0',np)
     endif
  endif
  read (iunps, *, err=300, iostat=ios) rgrid(np)%zmesh, rgrid(np)%xmin, rgrid(np)%dx, &
                                       rgrid(np)%mesh, nchi(np)
  if (rgrid(np)%mesh > ndmx .or. rgrid(np)%mesh <= 0) &
       call errore ('read_ncpp', 'mesh too big', np)
  if ( nchi(np) > nchix .or. &
       (nchi(np) < lmax(np)   .and. lloc(np) == lmax(np)) .or. & 
       (nchi(np) < lmax(np)+1 .and. lloc(np) /= lmax(np)) ) &
       call errore ('read_ncpp', 'wrong no. of wfcts', np)
  !
  !  Here pseudopotentials in numeric form are read
  !
  allocate (vnl(rgrid(np)%mesh, 0:lmax(np)))
  if (numeric (np) ) then
     do l = 0, lmax (np)
        read (iunps, '(a)', err=300, iostat=ios)
        read (iunps, *, err=300, iostat=ios) (vnl(ir,l), ir=1,rgrid(np)%mesh )
     enddo
     if (nlcc (np) ) then
        read (iunps, *, err=300, iostat=ios) (rho_atc(ir,np), ir=1,rgrid(np)%mesh)
     endif
  endif
  !
  !  Here pseudowavefunctions (in numeric form) are read
  !
  do nb = 1, nchi (np)
     read (iunps, '(a)', err=300, iostat=ios)
     read (iunps, *, err=300, iostat=ios) lchi(nb,np), oc(nb,np)
     !
     !     Test lchi and occupation numbers
     !
     if (nb <= lmax(np) .and. lchi(nb,np)+1 /= nb) &
          call errore ('read_ncpp', 'order of wavefunctions', 1)
     if (lchi(nb,np) > lmaxx .or. lchi(nb,np) < 0) &
                      call errore ('read_ncpp', 'wrong lchi', np)
     if (oc(nb,np) < 0.d0 .or. oc(nb,np) > 2.d0*(2*lchi(nb,np)+1)) &
          call errore ('read_ncpp', 'wrong oc', np)
     read (iunps, *, err=300, iostat=ios) ( chi(ir,nb,np), ir=1,rgrid(np)%mesh )
  enddo
  !
  !====================================================================
  ! PP read: now setup 
  !
  call set_dft_from_name( dft_name )
  !
#if defined (EXX)
#else
  IF ( dft_is_hybrid() ) &
    CALL errore( 'read_ncpp ', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  !    compute the radial mesh
  !
  do ir = 1, rgrid(np)%mesh
     x = rgrid(np)%xmin + DBLE (ir - 1) * rgrid(np)%dx
     rgrid(np)%r(ir) = exp (x) / rgrid(np)%zmesh 
     rgrid(np)%rab(ir) = rgrid(np)%dx * rgrid(np)%r(ir)
  enddo
  do ir = 1, rgrid(np)%mesh
     if ( rgrid(np)%r(ir) > rcut) then
        kkbeta(np) = ir
        go to 5
     end if
  end do
  kkbeta(np) = rgrid(np)%mesh
  !
  ! ... force kkbeta to be odd for simpson integration (obsolete?)
  !
5 kkbeta(np) = 2 * ( ( kkbeta(np) + 1 ) / 2) - 1
  !
  vloc_at (:, np) = 0.d0
  if (.not. numeric(np)) then
     !
     ! bring analytic potentials into numerical form
     !
     IF ( nlc(np) == 2 .AND. nnl(np) == 3 .AND. bhstype ) &
          CALL bachel( alps(1,0,np), aps(1,0,np), 1, lmax(np) )
     !
     do i = 1, nlc (np)
        do ir = 1, kkbeta(np)
           vloc_at (ir, np) = vloc_at (ir, np) - zp(np) * e2 * &
                 cc (i, np) * erf ( sqrt (alpc (i, np)) * rgrid(np)%r(ir) ) &
                 / rgrid(np)%r(ir)
        end do
     end do
     do l = 0, lmax (np)
        vnl (:, l) = vloc_at (1:rgrid(np)%mesh,np)
        do i = 1, nnl (np)
           vnl (:, l) = vnl (:, l) + e2 * (aps (i, l, np) + &
                   aps (i + 3, l, np) * rgrid(np)%r (:) **2) * &
                   exp ( - rgrid(np)%r(:) **2 * alps (i, l, np) )
        enddo
     enddo
     ! core corrections are still analytic!
     !!! numeric(np) =.true.
  end if
  !
  ! assume l=lloc as local part and subtract from the other channels
  !
  if (lloc (np) <= lmax (np) ) vloc_at (1:rgrid(np)%mesh, np) = vnl (:, lloc (np))
  ! lloc > lmax is allowed for PP in analytical form only
  ! it means that only the erf part is taken as local part 
  do l = 0, lmax (np)
     if (l /= lloc(np)) vnl (:, l) = vnl(:, l) - vloc_at(1:rgrid(np)%mesh, np)
  enddo
  !
  !    compute the atomic charges
  !
  rho_at(:,np) = 0.d0
  do nb = 1, nchi (np)
     if (oc (nb, np) > 0.d0) then
        do ir = 1, rgrid(np)%mesh
           rho_at(ir,np) = rho_at(ir,np) + oc(nb,np) * chi(ir,nb,np)**2
        enddo
     endif
  enddo
  !====================================================================
  ! convert to separable (KB) form
  !
  dion (:,:,np) = 0.d0
  nb = 0
  do l = 0, lmax (np)
     if (l /= lloc (np) ) then
        nb = nb + 1
        ! betar is used here as work space
        do ir = 1, kkbeta(np)
           betar (ir, nb, np) = chi(ir, l+1, np) **2 * vnl(ir, l)
        end do
        call simpson (kkbeta (np), betar (1, nb, np), rgrid(np)%rab, vll )
        dion (nb, nb, np) = 1.d0 / vll
        ! betar stores projectors  |beta(r)> = |V_nl(r)phi(r)>
        do ir = 1, kkbeta (np)
           betar (ir, nb, np) = vnl (ir, l) * chi (ir, l + 1, np)
        enddo
        lll (nb, np) = l
     endif
  enddo
  nbeta (np) = nb
  deallocate (vnl)

  return

300 call errore ('read_ncpp', 'pseudo file is empty or wrong', abs (np) )
end subroutine read_ncpp

