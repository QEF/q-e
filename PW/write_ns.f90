!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_ns
  !-----------------------------------------------------------------------

  use pwcom
  implicit none
  integer :: is, na, nt, m1, m2, ldim
  ! counter on spin component
  ! counter on atoms and their type
  ! counters on d components
  integer, parameter :: ldmx = 7
  complex(kind=DP) :: f (ldmx, ldmx), vet (ldmx, ldmx)
  real(kind=DP) :: lambda (ldmx), nsum

  write (*,*) 'enter write_ns'

  if ( 2 * Hubbard_lmax + 1 .gt. ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('U(',nt,') =', Hubbard_U(nt) * rytoev, nt=1,ntyp)
  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('alpha(',nt,') =', Hubbard_alpha(nt) * rytoev, nt=1,ntyp)

  nsum = 0.d0
  do is = 1, nspin
     do na = 1, nat
        nt = ityp (na)
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then
           ldim = 2 * Hubbard_l(nt) + 1
           do m1 = 1, ldim
              nsum = nsum + nsnew (na, is, m1, m1)
              do m2 = 1, ldim
                 f (m1, m2) = nsnew (na, is, m1, m2)
              enddo
           enddo
           call cdiagh(ldim, f, ldmx, lambda, vet)
           write(6,'(a,x,i2,2x,a,x,i2)') 'atom', na, 'spin', is
           write(6,'(a,7f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,ldim)
           write(6,*) 'eigenvectors'
           do m2 = 1, ldim
              write(6,'(i2,2x,7(f10.7,x))') m2,(dreal(vet(m1,m2)),m1=1,ldim)
           end do
           write(6,*) 'occupations'
           do m1 = 1, ldim
              write (6,'(7(f6.3,x))') (nsnew(na,is,m1,m2),m2=1,ldim)
           end do
        endif
     enddo
  enddo

  write (6, '(a,x,f11.7)') 'nsum =', nsum
  write (*,*) 'exit write_ns'
  return
end subroutine write_ns
