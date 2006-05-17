!
!---------------------------------------------------------------
subroutine test_bessel ( )
  !---------------------------------------------------------------
  !
  !     diagonalization of the pseudo-atomic hamiltonian 
  !     in a basis of spherical Bessel functions: j_l (qr)
  !     with zero boundary conditions at the border R: j_l (qR) = 0
  !
  !     units: atomic rydberg units
  !     ecut (Ry)   = energy cutoff (Ry)
  !     rm          = radius of the box         (default  30 a.u.)
  !
  !
  use ld1inc, only: lmax, lmx, mesh, r, r2, dx
  !
  implicit none
  !
  real*8, parameter :: pi=3.141592653589793d0
  real*8 ::  ecut =40.0,      & ! kinetic-energy cutoff (equivalent to PW)
             rm =  30.d0         ! radius R of the box
  !
  real*8, allocatable :: q(:,:)  ! quantized momenta in the spherical box
  !
  integer :: nswx, &            ! max number of quantized momenta q
             nsw (0:lmx), &     ! number of q such that q^2 < Ecut
             mesh_              ! mesh_ is such that r (mesh_) <= R
  !
  !   estimated max number of q needed for all values of l
  !
  nswx = int ( sqrt ( ecut*rm**2/pi**2 ) ) + lmax + 1
  !
  !   find grid point mesh_ such that r(mesh_) < R
  !   note that mesh_ must be odd to perform simpson integration
  !
  do mesh_ = 1, mesh
     if ( r(mesh_) >= rm ) go to 20
  end do
  call errore ('test_bessel','r(mesh) < Rmax', mesh_)
20 continue
  mesh_ = 2 * ( (mesh_-1) / 2 ) + 1 
  !
  !   find quantized momenta q
  !
  allocate ( q ( nswx, 0:lmax ) )
  call q_fill ( nswx, lmax, rm, ecut, nsw, q )
  !
  !   fill and diagonalize Kohn-Sham pseudo-hamiltonian
  !
  write (6, "(/,5x,14('-'), ' Test with a basis set of Bessel functions ',&
       & 10('-'),/)")
  write (6, "(5x,'Cutoff (Ry) : ',f6.1,'   Box size (a.u.) : ',f6.1)" ) &
       ecut, rm
  call h_diag ( mesh_, r, r2, dx, nswx, nsw, lmax, q )
  write (6, "(/,5x,14('-'), ' End of Bessel function test ',24('-'),/)")
  !
end subroutine test_bessel
!
!-----------------------------------------------------------------------
subroutine q_fill ( nswx, lmax, rm, ecut, nsw, q )
  !-----------------------------------------------------------------------
  !
  !   find quantized momenta in a box of radius R = rm
  !
  implicit none
  integer, intent(in) :: nswx, lmax
  real*8, intent(in)  :: rm, ecut
  integer, intent(out) :: nsw(0:lmax)
  real*8, intent(out) :: q(nswx,0:lmax)
  !
  integer :: i, l, iret
  real*8, parameter :: pi=3.141592653589793d0, eps=1.0e-10
  real*8, external :: find_root
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
  implicit none
  !
  integer, intent (in) :: l
  real*8, intent (in) :: xt1, xt2, eps
  !
  integer, intent (out) :: iret
  real*8 :: find_root
  !
  real*8 :: x1, x2, x0, f1, f2, f0
  !
  x1 = MIN (xt1, xt2) 
  x2 = MAX (xt1, xt2)
  !
  call sph_bes ( 1, 1.d0, x1, l, f1 )
  call sph_bes ( 1, 1.d0, x2, l, f2 )
  !
  iret = 0
  !
  if ( sign(f1,f2) == f1 ) then
     iret = 1
     return
  end if
!
  do while ( abs(x2-x1) > eps )
     x0 = 0.5*(x1+x2)
     call sph_bes ( 1, 1.d0, x0, l, f0 )
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
subroutine h_diag  ( mesh_, r, r2, dx, nswx, nsw, lmax, q )
  !-----------------------------------------------------------------------
  !
  ! diagonalize the radial Kohn-Sham hamiltonian 
  ! in the basis of spherical bessel functions
  ! Requires the self-consistent potential from a previous calculation!
  !
  use ld1inc, only: lmx, nbeta, betas, qq, ddd, vpstot, vnl, lls, jjs, &
       nspin, rel, pseudotype
  implicit none
  !
  ! input
  integer, intent (in) :: mesh_, lmax, nswx, nsw(0:lmx)
  real*8, intent (in) :: r(mesh_), r2(mesh_), dx
  real*8, intent (in) :: q(nswx,0:lmax)
  ! local
  real*8, allocatable :: h(:,:), s0(:), & ! hamiltonian and overlap matrix
                         chi (:,:), enl (:), & ! eigenvectors, eigenvalues
                         betajl(:,:), betajld(:,:), vaux (:), & ! workspace
                         jlq (:,:), work (:) ! workspace
  real*8, external :: int_0_inf_dr
  real*8 :: j
  character(len=2), dimension (2) :: spin = [ 'up', 'dw' ]
  integer :: n_states = 3, l, n, m, nb, mb, ind, is, nj
  !
  !
  allocate ( h (nswx, nswx),  s0(nswx), chi (nswx, n_states), enl(nswx) )
  allocate ( jlq ( mesh_, nswx ), work (mesh_), vaux (mesh_) )
  !
  write(6,"( 20x,3(7x,'N = ',i1) )" ) (n, n=1,n_states)
  !
  do l=0,lmax
     !
     !  nj: number of j-components in fully relativistic case
     !
     if ( rel < 2 .or. l == 0 ) then
        nj=1
     else
        nj=2
     endif
     !
     do n = 1, nsw(l)
        !
        call sph_bes ( mesh_, r, q(n,l), l, jlq(1,n) )
        !
        !  s0 is the orthonormalization factor for j_l
        !
        work (:) = ( jlq(:,n) * r(1:mesh_) ) ** 2
        s0(n) = sqrt ( int_0_inf_dr ( work, r, r2, dx, mesh_, 2*l+2 ) )
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
           if ( nj  == 2 ) then
              j = l + (ind-1.5d0) ! this is J=L+S
           else
              j = 0.d0
           end if
           !
           if (pseudotype == 1) then
              vaux(:) = vpstot (1:mesh_, is) + vnl (1:mesh_, l, ind)
           else
              vaux(:) = vpstot (1:mesh_, is)
              allocate ( betajl ( nswx, nbeta ), betajld ( nswx, nbeta ) )
           endif
           !
           h (:,:) = 0.d0
           !
           do n = 1, nsw(l)
              !
              !  matrix elements for vaux
              !
              do m = 1, n
                 work (:) = jlq(:,n) * jlq(:,m) * vaux(1:mesh_) * r2(1:mesh_)
                 h(m,n) = int_0_inf_dr ( work, r, r2, dx, mesh_, 2*l+2 ) &
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
                       work (:) = jlq(:,n) * betas( 1:mesh_, nb ) * r(1:mesh_)
                       betajl (n, nb) = 1.d0 / s0(n) * &
                            int_0_inf_dr ( work, r, r2, dx, mesh_, 2*l+2 )
                    end if
                 end do
              end if
           end do
           !
           !    betajld(q,m) = \sum_n D_mn * < beta_n | j_l(qr) >
           !
           if ( pseudotype > 1 ) then
              betajld (:,:) = 0.d0
              do mb = 1, nbeta
                 if ( lls (mb) == l .and. abs(jjs (mb) - j) < 0.001_8 ) then
                    do nb = 1, nbeta
                       if ( lls (nb) == l .and. abs(jjs (nb) - j) < 0.001_8 ) then
                          betajld (:, mb) = betajld (:, mb) + &
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
                          h(m,n) = h(m,n) + betajld (m, nb) * betajl (n, nb)
                       end if
                    end do
                 end do
              end do
              deallocate ( betajld, betajl )
           end if
           !
           call rdiagd ( nsw(l), h, nswx, n_states, enl, chi, nswx)
           !
           if ( nspin == 2 ) then
              write(6,"( 5x,'E(L=',i1,',spin 'a2,') =',4(f10.4,' Ry') )" ) &
                   l, spin(is), (enl(n), n=1,n_states)
           else if ( rel == 2 ) then
              write(6,"( 5x,'E(L=',i1,',J=',f3.1,') =',4(f10.4,' Ry') )" ) &
                   l, j, (enl(n), n=1,n_states)
           else 
              write(6,"( 5x,'E(L=',i1,') =',5x,4(f10.4,' Ry') )" ) &
                   l, (enl(n), n=1,n_states)
           end if
        end do
        !
     end do
  end do
  !
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
  !     LAPACK driver for matrix diagonalizer
  !
  implicit none
  !
  integer, intent (in) :: &
       n,      & ! dimension of the matrix to be diagonalized
       ldh,    & ! leading dimension of h, as declared in the calling pgm
       m,      & ! number of roots to be searched
       ldv       ! leading dimension of the v matrix
  real*8, intent (inout) :: &
       h(ldh,n) ! matrix to be diagonalized, UPPER triangle
  !
  real*8, intent (out) :: e(m), v(ldv,m) ! eigenvalues and eigenvectors 
  !
  integer :: i, j, mo, lwork, info, iwork(5*n), ifail(n)
  real*8  :: vl, vu, work(8*n)
  !
  lwork = 8*n
  v (:,:) = 0.d0

  call DSYEVX ( 'V', 'I', 'U', n, h, ldh, vl, vu, 1, m, 0.0d0, mo, e,&
               v, ldv, work, lwork, iwork, ifail, info )
  
  if ( info > 0) then
     call errore('rdiagd','failed to converge',info)
  else if(info < 0) then
     call errore('rdiagd','illegal arguments',-info)
  end if

  return
end subroutine rdiagd
