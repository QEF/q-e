!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine rgd_blk (nr1,nr2,nr3,nat,dyn,q,tau,epsil,zeu,bg,omega,sign)
  !-----------------------------------------------------------------------
  ! compute the rigid-ion (long-range) term for q
  ! The long-range term used here, to be added to or subtracted from the
  ! dynamical matrices, is exactly the same of the formula introduced in:
  ! X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is
  ! implemented: the Ewald parameter alpha must be large enough to
  ! have negligible r-space contribution
  !
  use kinds, only: dp
  use constants, only: pi, fpi, e2
  implicit none
  integer ::  nr1, nr2, nr3    !  FFT grid
  integer ::  nat              ! number of atoms
  complex(DP) :: dyn(3,3,nat,nat) ! dynamical matrix
  real(DP) &
       q(3),           &! q-vector
       tau(3,nat),     &! atomic positions
       epsil(3,3),     &! dielectric constant tensor
       zeu(3,3,nat),   &! effective charges tensor
       at(3,3),        &! direct     lattice basis vectors
       bg(3,3),        &! reciprocal lattice basis vectors
       omega,          &! unit cell volume
       sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term
  !
  ! local variables
  !
  real(DP):: geg                    !  <q+G| epsil | q+G>
  integer :: na,nb, i,j, m1, m2, m3
  integer :: nr1x, nr2x, nr3x
  real(DP) :: alph, fac,g1,g2,g3, facgd, arg, gmax
  real(DP) :: zag(3),zbg(3),zcg(3), fnat(3)
  complex(dp) :: facg
  !
  ! alph is the Ewald parameter, geg is an estimate of G^2
  ! such that the G-space sum is convergent for that alph
  ! very rough estimate: geg/4/alph > gmax = 14
  ! (exp (-14) = 10^-6)
  !
  gmax= 14.d0
  alph= 1.0d0
  geg = gmax*alph*4.0d0

  ! Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
  ! Only for dimensions where periodicity is present, e.g. if nr1=1
  ! and nr2=1, then the G-vectors run along nr3 only.
  ! (useful if system is in vacuum, e.g. 1D or 2D)
  !
  if (nr1 == 1) then
     nr1x=0
  else
     nr1x = int ( sqrt (geg) / &
                  sqrt (bg (1, 1) **2 + bg (2, 1) **2 + bg (3, 1) **2) ) + 1
  endif
  if (nr2 == 1) then
     nr2x=0
  else
     nr2x = int ( sqrt (geg) / &
                  sqrt (bg (1, 2) **2 + bg (2, 2) **2 + bg (3, 2) **2) ) + 1
  endif
  if (nr3 == 1) then
     nr3x=0
  else
     nr3x = int ( sqrt (geg) / &
                  sqrt (bg (1, 3) **2 + bg (2, 3) **2 + bg (3, 3) **2) ) + 1
  endif
  !
  if (abs(sign) /= 1.0_DP) &
       call errore ('rgd_blk',' wrong value for sign ',1)
  !
  fac = sign*e2*fpi/omega
  do m1 = -nr1x,nr1x
  do m2 = -nr2x,nr2x
  do m3 = -nr3x,nr3x
     !
     g1 = m1*bg(1,1) + m2*bg(1,2) + m3*bg(1,3)
     g2 = m1*bg(2,1) + m2*bg(2,2) + m3*bg(2,3)
     g3 = m1*bg(3,1) + m2*bg(3,2) + m3*bg(3,3)
     !
     geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
            g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
            g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
     !
     if (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) then
        !
        facgd = fac*exp(-geg/alph/4.0d0)/geg
        !
        do na = 1,nat
           zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
           fnat(:) = 0.d0
           do nb = 1,nat
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              zcg(:) = g1*zeu(1,:,nb) + g2*zeu(2,:,nb) + g3*zeu(3,:,nb)
              fnat(:) = fnat(:) + zcg(:)*cos(arg)
           end do
           do j=1,3
              do i=1,3
                 dyn(i,j,na,na) = dyn(i,j,na,na) - facgd * &
                                  zag(i) * fnat(j)
              end do
           end do
        end do
     end if
     !
     g1 = g1 + q(1)
     g2 = g2 + q(2)
     g3 = g3 + q(3)
     !
     geg = (g1*(epsil(1,1)*g1+epsil(1,2)*g2+epsil(1,3)*g3)+      &
            g2*(epsil(2,1)*g1+epsil(2,2)*g2+epsil(2,3)*g3)+      &
            g3*(epsil(3,1)*g1+epsil(3,2)*g2+epsil(3,3)*g3))
     !
     if (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) then
        !
        facgd = fac*exp(-geg/alph/4.0d0)/geg
        !
        do nb = 1,nat
           zbg(:)=g1*zeu(1,:,nb)+g2*zeu(2,:,nb)+g3*zeu(3,:,nb)
           do na = 1,nat
              zag(:)=g1*zeu(1,:,na)+g2*zeu(2,:,na)+g3*zeu(3,:,na)
              arg = 2.d0*pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                              g2 * (tau(2,na)-tau(2,nb))+             &
                              g3 * (tau(3,na)-tau(3,nb)))
              !
              facg = facgd * CMPLX(cos(arg),sin(arg),kind=DP)
              do j=1,3
                 do i=1,3
                    dyn(i,j,na,nb) = dyn(i,j,na,nb) + facg *      &
                                     zag(i) * zbg(j)
                 end do
              end do
           end do
        end do
     end if
  end do
  end do
  end do
  !
  return
  !
end subroutine rgd_blk
!
!-----------------------------------------------------------------------
subroutine nonanal(nat, nat_blk, itau_blk, epsil, q, zeu, omega, dyn )
  !-----------------------------------------------------------------------
  !     add the nonanalytical term with macroscopic electric fields
  !
  use kinds, only: dp
  use constants, only: pi, fpi, e2
 implicit none
 integer, intent(in) :: nat, nat_blk, itau_blk(nat)
 !  nat: number of atoms in the cell (in the supercell in the case
 !       of a dyn.mat. constructed in the mass approximation)
 !  nat_blk: number of atoms in the original cell (the same as nat if
 !       we are not using the mass approximation to build a supercell)
 !  itau_blk(na): atom in the original cell corresponding to
 !                atom na in the supercell
 !
 complex(DP), intent(inout) :: dyn(3,3,nat,nat) ! dynamical matrix
 real(DP), intent(in) :: q(3),  &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nat_blk),   &! effective charges tensor
      &       omega            ! unit cell volume
 !
 ! local variables
 !
 real(DP) zag(3),zbg(3),  &! eff. charges  times g-vector
      &       qeq              !  <q| epsil | q>
 integer na,nb,              &! counters on atoms
      &  na_blk,nb_blk,      &! as above for the original cell
      &  i,j                  ! counters on cartesian coordinates
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
!print*, q(1), q(2), q(3)
 if (qeq < 1.d-8) then
    write(6,'(5x,"A direction for q was not specified:", &
      &          "TO-LO splitting will be absent")')
    return
 end if
 !
 do na = 1,nat
    na_blk = itau_blk(na)
    do nb = 1,nat
       nb_blk = itau_blk(nb)
       !
       do i=1,3
          !
          zag(i) = q(1)*zeu(1,i,na_blk) +  q(2)*zeu(2,i,na_blk) + &
                   q(3)*zeu(3,i,na_blk)
          zbg(i) = q(1)*zeu(1,i,nb_blk) +  q(2)*zeu(2,i,nb_blk) + &
                   q(3)*zeu(3,i,nb_blk)
       end do
       !
       do i = 1,3
          do j = 1,3
             dyn(i,j,na,nb) = dyn(i,j,na,nb)+ fpi*e2*zag(i)*zbg(j)/qeq/omega
!             print*, zag(i),zbg(j),qeq, fpi*e2*zag(i)*zbg(j)/qeq/omega
          end do
       end do
    end do
 end do
 !
 return
end subroutine nonanal

!-----------------------------------------------------------------------
subroutine nonanal_ifc(nat, nat_blk, itau_blk, epsil, q, zeu, omega, dyn, nr1,nr2,nr3,f_of_q )
  !-----------------------------------------------------------------------
  !     add the nonanalytical term with macroscopic electric fields
  !
  use kinds, only: dp
  use constants, only: pi, fpi, e2
 implicit none
 integer, intent(in) :: nat, nat_blk, itau_blk(nat), nr1,nr2,nr3
 !  nat: number of atoms in the cell (in the supercell in the case
 !       of a dyn.mat. constructed in the mass approximation)
 !  nat_blk: number of atoms in the original cell (the same as nat if
 !       we are not using the mass approximation to build a supercell)
 !  itau_blk(na): atom in the original cell corresponding to
 !                atom na in the supercell
 !
 complex(DP), intent(inout) :: dyn(3,3,nat,nat),f_of_q(3,3,nat,nat) ! dynamical matrix
 real(DP), intent(in) :: q(3),  &! polarization vector
      &       epsil(3,3),     &! dielectric constant tensor
      &       zeu(3,3,nat_blk),   &! effective charges tensor
      &       omega            ! unit cell volume
 !
 ! local variables
 !
 real(DP) zag(3),zbg(3),  &! eff. charges  times g-vector
      &       qeq              !  <q| epsil | q>
 integer na,nb,              &! counters on atoms
      &  na_blk,nb_blk,      &! as above for the original cell
      &  i,j                  ! counters on cartesian coordinates
 !
IF ( q(1)==0.d0 .AND. &
     q(2)==0.d0 .AND. &
     q(3)==0.d0 ) return
 !
 qeq = (q(1)*(epsil(1,1)*q(1)+epsil(1,2)*q(2)+epsil(1,3)*q(3))+    &
        q(2)*(epsil(2,1)*q(1)+epsil(2,2)*q(2)+epsil(2,3)*q(3))+    &
        q(3)*(epsil(3,1)*q(1)+epsil(3,2)*q(2)+epsil(3,3)*q(3)))
 !
!print*, q(1), q(2), q(3)
 if (qeq < 1.d-8) then
    write(6,'(5x,"A direction for q was not specified:", &
      &          "TO-LO splitting will be absent")')
    return
 end if

 do na = 1,nat
    na_blk = itau_blk(na)
    do nb = 1,nat
       nb_blk = itau_blk(nb)
       !
       do i=1,3
          !
          zag(i) = q(1)*zeu(1,i,na_blk) +  q(2)*zeu(2,i,na_blk) + &
                   q(3)*zeu(3,i,na_blk)
          zbg(i) = q(1)*zeu(1,i,nb_blk) +  q(2)*zeu(2,i,nb_blk) + &
                   q(3)*zeu(3,i,nb_blk)
       end do
       !
       do i = 1,3
          do j = 1,3
!             dyn(i,j,na,nb) = dyn(i,j,na,nb)+ fpi*e2*zag(i)*f_of_q*zbg(j)/qeq/omega/(nr1*nr2*nr3)
              f_of_q(i,j,na,nb)=fpi*e2*zag(i)*zbg(j)/qeq/omega/(nr1*nr2*nr3)
!             print*, i,j,na,nb, dyn(i,j,na,nb)
          end do
       end do
    end do
 end do
 !
 return
end subroutine nonanal_ifc
!
!-----------------------------------------------------------------------
subroutine dyndiag (nat,ntyp,amass,ityp,dyn,w2,z)
  !-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On input:  amass = masses, in amu
  !   On output: w2 = energies, z = displacements
  !
  use kinds, only: dp
  use constants, only: amu_ry
  implicit none
  ! input
  integer nat, ntyp, ityp(nat)
  complex(DP) dyn(3,3,nat,nat)
  real(DP) amass(ntyp)
  ! output
  real(DP) w2(3*nat)
  complex(DP) z(3*nat,3*nat)
  ! local
  real(DP) diff, dif1, difrel
  integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
  complex(DP), allocatable :: dyn2(:,:)
  !
  !  fill the two-indices dynamical matrix
  !
  nat3 = 3*nat
  allocate(dyn2 (nat3, nat3))
  !
  do na = 1,nat
     do nb = 1,nat
        do ipol = 1,3
           do jpol = 1,3
              dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = dyn(ipol,jpol,na,nb)
           end do
        end do
     end do
  end do
  !
  !  impose hermiticity
  !
  diff = 0.d0
  difrel=0.d0
  do i = 1,nat3
     dyn2(i,i) = CMPLX( DBLE(dyn2(i,i)),0.d0,kind=DP)
     do j = 1,i - 1
        dif1 = abs(dyn2(i,j)-CONJG(dyn2(j,i)))
        if ( dif1 > diff .and. &
             max ( abs(dyn2(i,j)), abs(dyn2(j,i))) > 1.0d-6) then
           diff = dif1
           difrel=diff / min ( abs(dyn2(i,j)), abs(dyn2(j,i)))
        end if
        dyn2(i,j) = 0.5d0* (dyn2(i,j)+CONJG(dyn2(j,i)))
        dyn2(j,i) = CONJG(dyn2(i,j))
     end do
  end do
  if ( diff > 1.d-6 ) write (6,'(5x,"Max |d(i,j)-d*(j,i)| = ",f9.6,/,5x, &
       & "Max |d(i,j)-d*(j,i)|/|d(i,j)|: ",f8.4,"%")') diff, difrel*100
  !
  !  divide by the square root of masses
  !
  do na = 1,nat
     nta = ityp(na)
     do nb = 1,nat
        ntb = ityp(nb)
        do ipol = 1,3
           do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = &
                  dyn2((na-1)*3+ipol, (nb-1)*3+jpol) / &
                  (amu_ry*sqrt(amass(nta)*amass(ntb)))
          end do
       end do
    end do
 end do
 !
 !  diagonalisation
 !
 call cdiagh2(nat3,dyn2,nat3,w2,z)
 !
 deallocate(dyn2)
 !
 !  displacements are eigenvectors divided by sqrt(amass)
 !
 do i = 1,nat3
    do na = 1,nat
       nta = ityp(na)
       do ipol = 1,3
          z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)/ sqrt(amu_ry*amass(nta))
       end do
    end do
 end do
 !
 return
end subroutine dyndiag
!
!-----------------------------------------------------------------------
subroutine cdiagh2 (n,h,ldh,e,v)
  !-----------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
  use kinds, only: dp
  implicit none
  !
  ! on INPUT
  integer          n,       &! dimension of the matrix to be diagonalized
       &           ldh       ! leading dimension of h, as declared
  ! in the calling pgm unit
  complex(DP)  h(ldh,n)  ! matrix to be diagonalized
  !
  ! on OUTPUT
  real(DP)     e(n)      ! eigenvalues
  complex(DP)  v(ldh,n)  ! eigenvectors (column-wise)
  !
  ! LOCAL variables (LAPACK version)
  !
  integer          lwork,   &! aux. var.
       &           ILAENV,  &! function which gives block size
       &           nb,      &! block size
       &           info      ! flag saying if the exec. of libr. routines was ok
  !
  real(DP), allocatable::   rwork(:)
  complex(DP), allocatable:: work(:)
  !
  !     check for the block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  if (nb.lt.1) nb=max(1,n)
  if (nb.eq.1.or.nb.ge.n) then
     lwork=2*n-1
  else
     lwork = (nb+1)*n
  endif
  !
  ! allocate workspace
  !
  call zcopy(n*ldh,h,1,v,1)
  allocate(work (lwork))
  allocate(rwork (3*n-2))
  call ZHEEV('V','U',n,v,ldh,e,work,lwork,rwork,info)
  call errore ('cdiagh2','info =/= 0',abs(info))
  ! deallocate workspace
  deallocate(rwork)
  deallocate(work)
  !
  return
end subroutine cdiagh2
