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
  integer :: is, na, nt, m1, m2
  ! counter on spin component
  ! counter on atoms and their type
  ! counters on d components
  complex(kind=DP) :: f (5, 5), vet (5, 5)


  real(kind=DP) :: lambda (5), nsum
  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('U(',nt,') =', Hubbard_U(nt) * rytoev, nt=1,ntyp)
  write (6,'(6(a,i2,a,f8.4,6x))') &
        ('alpha(',nt,') =', Hubbard_alpha(nt) * rytoev, nt=1,ntyp)

  nsum = 0.d0
  do is = 1, nspin
     do na = 1, nat
        nt = ityp (na)
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then
           do m1 = 1, 5
              nsum = nsum + nsnew (na, is, m1, m1)
              do m2 = 1, 5
                 f (m1, m2) = nsnew (na, is, m1, m2)
              enddo
           enddo
           call cdiagh(5,f,5,lambda,vet)
           write(6,'(a,x,i2,2x,a,x,i2)') 'atom', na, 'spin', is
           write(6,'(a,5f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,5)
           write(6,*) 'eigenvectors'
           do m2 = 1,5
              write(6,'(i2,2x,5(f10.7,x))') &
                        m2,(dreal(vet(m1,m2)),m1=1,5)
           end do
           write(6,*) 'occupations'
           write (6,'(5(f6.3,x))') &
                 ((nsnew(na,is,m1,m2),m2=1,5),m1=1,5)
        endif
     enddo
  enddo

  write (6, '(a,x,f11.7)') 'nsum =', nsum
  return
end subroutine write_ns
