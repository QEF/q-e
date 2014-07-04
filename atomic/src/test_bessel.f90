!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine test_bessel ( )
  !---------------------------------------------------------------
  !
  !     diagonalization of the pseudo-atomic hamiltonian 
  !     in a basis of spherical Bessel functions: j_l (qr)
  !     with zero boundary conditions at the border R: 
  !     j_l (q_l R) = 0
  !     The basis sets contains all q_l's up to a kinetic
  !     energy cutoff Ecut, such that q^2 <= Ecut (in Ry a.u.)
  !     rm is the radius R of the box
  !
  use io_global, only : stdout
  use kinds, only : dp
  use constants, only: pi
  use ld1inc, only: lmax, lmx, grid
  use ld1inc, only: ecutmin, ecutmax, decut, rm
  !
  implicit none
  !
  real(kind=dp), parameter :: eps = 1.0d-4
  real(kind=dp) ::  ecut   ! kinetic-energy cutoff (equivalent to PW)
  !
  real(kind=dp), allocatable :: q(:,:) ! quantized momenta in the spherical box
  !
  integer :: nswx, &            ! max number of quantized momenta q
             nsw (0:lmx), &     ! number of q such that q^2 < Ecut
             mesh_, &           ! mesh_ is such that r (mesh_) <= R
             ncut, nc
  !
  if ( ecutmin < eps .or. ecutmax < eps .or. ecutmax < ecutmin + eps .or. &
       decut < eps .or. rm < 5.0_dp ) return
  !
  write (stdout, "(/,5x,14('-'), ' Test with a basis set of Bessel functions ',&
       & 10('-'),/)")
  !
  write (stdout, "(5x,'Box size (a.u.) : ',f6.1)" ) rm
  ncut = nint ( ( ecutmax-ecutmin ) / decut ) + 1
  !
  !  we redo everything for each cutoff: not really a smart implementation
  !
  do nc = 1, ncut
     !
     ecut = ecutmin + (nc-1) * decut
     !
     !   estimated max number of q needed for all values of l
     !
     nswx = int ( sqrt ( ecut*rm**2/pi**2 ) ) + lmax + 1
     !
     !   find grid point mesh_ such that r(mesh_) < R
     !   note that mesh_ must be odd to perform simpson integration
     !
     do mesh_ = 1, grid%mesh
        if ( grid%r(mesh_) >= rm ) go to 20
     end do
     call errore ('test_bessel','r(mesh) < Rmax', mesh_)
20   continue
     mesh_ = 2 * ( (mesh_-1) / 2 ) + 1 
     !
     !   find quantized momenta q
     !
     allocate ( q ( nswx, 0:lmax ) )
     call q_fill ( nswx, lmax, rm, ecut, nsw, q )
     !
     !   fill and diagonalize Kohn-Sham pseudo-hamiltonian
     !
     write (stdout, "(/5x,'Cutoff (Ry) : ',f6.1)" ) ecut
     call h_diag ( mesh_, nswx, nsw, lmax, q )
     !
     deallocate ( q )
     !
  end do
  !
  write (stdout, "(/,5x,14('-'), ' End of Bessel function test ',24('-'),/)")
  !
end subroutine test_bessel
!
!-----------------------------------------------------------------------
subroutine q_fill ( nswx, lmax, rm, ecut, nsw, q )
  !-----------------------------------------------------------------------
  !
  !   find quantized momenta in a box of radius R = rm
  !
  use kinds, only : dp
  use constants, only : pi
  implicit none
  integer, intent(in) :: nswx, lmax
  real(kind=dp), intent(in)  :: rm, ecut
  integer, intent(out) :: nsw(0:lmax)
  real(kind=dp), intent(out) :: q(nswx,0:lmax)
  !
  integer :: i, l, iret
  real(kind=dp), parameter :: eps=1.0d-10
  real(kind=dp), external :: find_root
  !
  !   l = 0 : j_0 (qR)=0  =>  q_i = i * (2pi/R)
  !
  do i = 1, nswx
     q(i,0) = i*pi
  end do
  !
  !   l > 0 : zeros for j_l (x) are found between two consecutive
  !           zeros for j_{l-1} (x)
  !
  do l = 1, lmax
     do i = 1, nswx-l
        q(i,l) = find_root ( l, q(i,l-1), q(i+1,l-1), eps, iret ) 
        if (iret /= 0)  call errore('q_fill','root not found',l)
     end do
  end do
  !
  do l = 0, lmax
     do i = 1, nswx-l
        q (i,l) = q(i,l) / rm
        if ( q(i,l)**2 > ecut ) then
           nsw(l) = i-1
           goto 20
        endif
     enddo
     call errore('q_fill','nswx is too small',nswx)
20   continue
  end do
  !
  return
end subroutine q_fill
!
!-----------------------------------------------------------------------
function find_root   ( l, xt1, xt2, eps, iret )
  ! ----------------------------------------------------------------------
  !
  use kinds, only : dp
  implicit none
  !
  integer, intent (in) :: l
  real(kind=dp), intent (in) :: xt1, xt2, eps
  !
  integer, intent (out) :: iret
  real(kind=dp) :: find_root
  !
  real(kind=dp) :: x1, x2, x0, f1, f2, f0
  !
  x1 = MIN (xt1, xt2) 
  x2 = MAX (xt1, xt2)
  !
  call sph_bes ( 1, 1.0_dp, x1, l, f1 )
  call sph_bes ( 1, 1.0_dp, x2, l, f2 )
  !
  iret = 0
  !
  if ( sign(f1,f2) == f1 ) then
     find_root = 0.0_dp
     iret = 1
     return
  end if
!
  do while ( abs(x2-x1) > eps )
     x0 = 0.5*(x1+x2)
     call sph_bes ( 1, 1.0_dp, x0, l, f0 )
     if ( sign(f0,f1) == f0 ) then
        x1 = x0
        f1 = f0
     else
        x2 = x0
        f2 = f0
     end if
  end do
  !
  find_root = x0
  !
  return
end function find_root
!
!-----------------------------------------------------------------------
subroutine h_diag  ( mesh_, nswx, nsw, lmax, q )
  !-----------------------------------------------------------------------
  !
  ! diagonalize the radial Kohn-Sham hamiltonian 
  ! in the basis of spherical bessel functions
  ! Requires the self-consistent potential from a previous calculation!
  !
  use io_global, only : stdout
  use kinds, only : dp
  use ld1inc, only: lmx, grid
  use ld1inc, only: nbeta, betas, qq, ddd, vpstot, vnl, lls, jjs, &
       nspin, rel, pseudotype
  implicit none
  !
  ! input
  integer, intent (in) :: mesh_, lmax, nswx, nsw(0:lmx)
  real(kind=dp), intent (in) :: q(nswx,0:lmax)
  ! local
  real(kind=dp), allocatable :: &
       h(:,:), s(:,:),     & ! hamiltonian and overlap matrix
       chi (:,:), enl (:), & ! eigenvectors and eigenvalues
       s0(:),              & ! normalization factors for jl(qr)
       betajl(:,:),        & ! matrix elements for beta
       betajl_(:,:),       & ! work space:  \sum_m D_lm beta_m
       vaux (:), jlq (:,:), work (:) ! more work space
  real(kind=dp), external :: int_0_inf_dr
  real(kind=dp) :: j
  character(len=2), dimension (2) :: spin = (/ 'up', 'dw' /)
  integer :: n_states = 3, l, n, m, nb, mb, ind, is, nj
  !
  !
  allocate ( h (nswx, nswx), s0(nswx), chi (nswx, n_states), enl(nswx) )
  allocate ( jlq ( mesh_, nswx ), work (mesh_), vaux (mesh_) )
  if ( pseudotype > 2 ) allocate ( s(nswx, nswx) )
  !
  write(stdout,"( 20x,3(7x,'N = ',i1) )" ) (n, n=1,n_states)
  !
  do l=0,lmax
     !
     !  nj: number of j-components in fully relativistic case
     !
     if ( rel == 2 .and. l > 0  ) then
        nj=2
     else
        nj=1
     endif
     !
     do n = 1, nsw(l)
        !
        call sph_bes ( mesh_, grid%r, q(n,l), l, jlq(1,n) )
        !
        !  s0 = < j_l(qr) | j_l (qr) >
        !
        work (:) = ( jlq(:,n) * grid%r(1:mesh_) ) ** 2
        s0(n) = sqrt ( int_0_inf_dr ( work, grid, mesh_, 2*l+2 ) )
        !
     end do
     !
     !    vaux is the local + scf potential ( + V_l for semilocal PPs)
     !    note that nspin = 2 only for lsda; nspin = 1 in all other cases
     !
     do is  = 1, nspin
        !
        do ind = 1, nj
           !
           if ( rel == 2 .and. l > 0 ) then
              j = l + (ind-1.5_dp) ! this is J=L+S
           else if ( rel == 2 .and. l == 0 ) then
              j = 0.5_dp
           else
              j = 0.0_dp
           end if
           !
           if (pseudotype == 1) then
              vaux(:) = vpstot (1:mesh_, is) + vnl (1:mesh_, l, ind)
           else
              vaux(:) = vpstot (1:mesh_, is)
              allocate ( betajl ( nswx, nbeta ), betajl_ ( nswx, nbeta ) )
           endif
           !
           h (:,:) = 0.0_dp
           !
           do n = 1, nsw(l)
              !
              !  matrix elements for vaux
              !
              do m = 1, n
                 work (:) = jlq(:,n) * jlq(:,m) * vaux(1:mesh_) * grid%r2(1:mesh_)
                 h(m,n) = int_0_inf_dr ( work, grid, mesh_, 2*l+2 ) &
                      / s0(n) / s0(m)
              end do
              !
              !    kinetic energy
              !
              h(n,n) =  h(n,n) + q(n,l)**2
              !
              !    betajl(q,n) = < beta_n | j_l(qr) >
              !
              if ( pseudotype > 1 ) then
                 do nb = 1, nbeta
                    if ( lls (nb) == l .and.  abs(jjs (nb) - j) < 0.001_8 ) then
                       work (:) = jlq(:,n) * betas( 1:mesh_, nb ) * grid%r(1:mesh_)
                       betajl (n, nb) = 1.0_dp / s0(n) * &
                            int_0_inf_dr ( work, grid, mesh_, 2*l+2 )
                    end if
                 end do
              end if
           end do
           !
           !  separable PP
           !
           if ( pseudotype > 1 ) then
              !
              !    betajl_(q,m) = \sum_n D_mn * < beta_n | j_l(qr) >
              !
              betajl_ (:,:) = 0.0_dp
              do mb = 1, nbeta
                 if ( lls (mb) == l .and. abs(jjs (mb) - j) < 0.001_8 ) then
                    do nb = 1, nbeta
                       if ( lls (nb) == l .and. abs(jjs (nb) - j) < 0.001_8 ) then
                          betajl_ (:, mb) = betajl_ (:, mb) + &
                               ddd (mb,nb,is) * betajl (:, nb)
                       end if
                    end do
                 end if
              end do
              !
              !    matrix elements for nonlocal (separable) potential
              !
              do n = 1, nsw(l)
                 do m = 1, n
                    do nb = 1, nbeta
                       if ( lls (nb) == l .and. abs(jjs (nb) - j) < 0.001_8 ) then
                          h(m,n) = h(m,n) + betajl_ (m, nb) * betajl (n, nb)
                       end if
                    end do
                 end do
              end do
              !
              !  US PP
              !
              if ( pseudotype > 2 ) then
                 !
                 !    betajl_(q,m) = \sum_n D_mn * < beta_n | j_l(qr) >
                 !
                 betajl_ (:,:) = 0.0_dp
                 do mb = 1, nbeta
                    if ( lls (mb) == l .and. abs(jjs (mb) - j) < 0.001_8 ) then
                       do nb = 1, nbeta
                          if ( lls (nb) == l .and. &
                               abs(jjs (nb) - j) < 0.001_8 ) then
                             betajl_ (:, mb) = betajl_ (:, mb) + &
                                  qq (mb,nb) * betajl (:, nb)
                          end if
                       end do
                    end if
                 end do
                 !
                 !    overlap matrix S
                 !
                 s (:,:) = 0.0_dp
                 do n = 1, nsw(l)
                    do m = 1, n
                       do nb = 1, nbeta
                          if ( lls (nb) == l .and. &
                               abs(jjs (nb) - j) < 0.001_8 ) then
                             s(m,n) = s(m,n) + betajl_ (m, nb) * betajl (n, nb)
                          end if
                       end do
                    end do
                    s(n,n) = s(n,n) + 1.0_dp
                 end do
              end if
              !
              deallocate ( betajl_, betajl )
              !
           end if
           !
           if ( pseudotype > 2 ) then
              call rdiags ( nsw(l), h, s, nswx, n_states, enl, chi, nswx)
           else
              call rdiagd ( nsw(l), h, nswx, n_states, enl, chi, nswx)
           end if
           !
           if ( nspin == 2 ) then
              write(stdout, &
                  "( 5x,'E(L=',i1,',spin ',a2,') =',4(f10.4,' Ry') )" ) &
                   l, spin(is), (enl(n), n=1,n_states)
           else if ( rel == 2 ) then
              write(stdout, &
                  "( 5x,'E(L=',i1,',J=',f3.1,') =',4(f10.4,' Ry') )" ) &
                   l, j, (enl(n), n=1,n_states)
           else 
              write(stdout, &
                  "( 5x,'E(L=',i1,') =',5x,4(f10.4,' Ry') )" ) &
                   l, (enl(n), n=1,n_states)
           end if
        end do
        !
     end do
  end do
  !
  if ( pseudotype > 2 ) deallocate ( s )
  deallocate ( vaux, work, jlq )
  deallocate ( enl, chi, s0, h )
  !
  return
end subroutine h_diag
!
!------------------------------------------------------------------
subroutine rdiagd ( n, h, ldh, m, e, v, ldv )
  !------------------------------------------------------------------
  !
  !     LAPACK driver for eigenvalue problem H*x = ex
  !
  use kinds, only : dp
  implicit none
  !
  integer, intent (in) :: &
       n,      & ! dimension of the matrix to be diagonalized
       ldh,    & ! leading dimension of h, as declared in the calling pgm
       m,      & ! number of roots to be searched
       ldv       ! leading dimension of the v matrix
  real(kind=dp), intent (inout) :: &
       h(ldh,n) ! matrix to be diagonalized, UPPER triangle
                ! DESTROYED ON OUTPUT
  real(kind=dp), intent (out) :: e(m), v(ldv,m) ! eigenvalues and eigenvectors 
  !
  integer :: mo, lwork, info
  integer, allocatable :: iwork(:), ifail(:)
  real(kind=dp)  :: vl, vu
  real(kind=dp), allocatable  :: work(:)
  !
  lwork = 8*n
  allocate ( work (lwork), iwork(5*n), ifail(n) )
  v (:,:) = 0.0_dp
  !
  call DSYEVX ( 'V', 'I', 'U', n, h, ldh, vl, vu, 1, m, 0.0_dp, mo, e,&
               v, ldv, work, lwork, iwork, ifail, info )
  !
  if ( info > 0) then
     call errore('rdiagd','failed to converge',info)
  else if(info < 0) then
     call errore('rdiagd','illegal arguments',-info)
  end if
  deallocate ( ifail, iwork, work )
  !
  return
end subroutine rdiagd
!
!----------------------------------------------------------------------------
SUBROUTINE rdiags( n, h, s, ldh, m, e, v, ldv )
  !----------------------------------------------------------------------------
  !
  !     LAPACK driver for generalized eigenvalue problem H*x = e S*x
  !
  use kinds, only : dp
  IMPLICIT NONE
  !
  integer, intent (in) :: &
       n,      & ! dimension of the matrix to be diagonalized
       ldh,    & ! leading dimension of h, as declared in the calling pgm
       m,      & ! number of roots to be searched
       ldv       ! leading dimension of the v matrix
  real(kind=dp), intent (inout) :: &
       h(ldh,n), & ! matrix to be diagonalized, UPPER triangle
       s(ldh,n)    ! overlap matrix - both destroyed on output
  !
  real(kind=dp), intent (out) :: e(m), v(ldv,m) ! eigenvalues and eigenvectors 
  !
  ! ... LOCAL variables
  !
  integer :: mo, lwork, info
  integer, allocatable :: iwork(:), ifail(:)
  real(kind=dp), allocatable  :: work(:)

  lwork = 8 * n
  allocate ( work (lwork), iwork(5*n), ifail(n) )
  v (:,:) = 0.0_dp
  !
  CALL DSYGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
       0.0_dp, 0.0_dp, 1, m, 0.0_dp, mo, e, v, ldh, work, lwork, &
       iwork, ifail, info )
  !             
  if ( info > n) then
     call errore('rdiags','failed to converge (factorization)',info-n)
  else if(info > 0) then
     call errore('rdiags','failed to converge: ',info)
  else if(info < 0) then
     call errore('rdiags','illegal arguments',-info)
  end if
  deallocate ( ifail, iwork, work )  
  !
  return
  !
END SUBROUTINE rdiags
