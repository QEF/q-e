!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_ns
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha
  !
  implicit none
  !
  integer :: is, na, nt, m1, m2, ldim
  ! counter on spin component
  ! counters on atoms and their type
  ! counters on d components
  integer, parameter :: ldmx = 7
  complex(DP) :: f (ldmx, ldmx), vet (ldmx, ldmx)
  real(DP) :: lambda (ldmx), nsum, nsuma

  WRITE( stdout,*) 'enter write_ns'

  if ( 2 * Hubbard_lmax + 1 > ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

  WRITE( stdout,'(6(a,i2,a,f8.4,6x))') &
        ('U(',nt,') =', Hubbard_U(nt) * rytoev, nt=1,ntyp)
  WRITE( stdout,'(6(a,i2,a,f8.4,6x))') &
        ('alpha(',nt,') =', Hubbard_alpha(nt) * rytoev, nt=1,ntyp)

  nsum = 0.d0
  do na = 1, nat
     nt = ityp (na)
     if (Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0) then
        ldim = 2 * Hubbard_l(nt) + 1
        nsuma = 0.d0
        do is = 1, nspin
           do m1 = 1, ldim
              nsuma = nsuma + rho%ns (m1, m1, is, na)
           end do
        end do
        if (nspin.eq.1) nsuma = 2.d0 * nsuma
        WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)]= ",f11.7)') na, nsuma
        nsum = nsum + nsuma
        do is = 1, nspin
           do m1 = 1, ldim
              do m2 = 1, ldim
                 f (m1, m2) = rho%ns (m1, m2, is, na)
              enddo
           enddo
           call cdiagh(ldim, f, ldmx, lambda, vet)
           WRITE( stdout,'("atom ",i4,3x,"spin ",i2)') na, is
           WRITE( stdout,'(a,7f10.7)') 'eigenvalues: ',(lambda(m1),m1=1,ldim)
           WRITE( stdout,*) 'eigenvectors'
           do m2 = 1, ldim
              WRITE( stdout,'(i2,2x,7(f10.7,1x))') m2,( DBLE(vet(m1,m2)),m1=1,ldim)
           end do
           WRITE( stdout,*) 'occupations'
           do m1 = 1, ldim
              WRITE( stdout,'(7(f6.3,1x))') (rho%ns(m1,m2,is,na),m2=1,ldim)
           end do
        enddo
     endif
  enddo

  WRITE( stdout, '(a,1x,f11.7)') 'nsum =', nsum
  WRITE( stdout,*) 'exit write_ns'
  return
end subroutine write_ns
