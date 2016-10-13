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
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                         Hubbard_alpha, lda_plus_u_kind, Hubbard_J0, &
                         Hubbard_beta
  !
  implicit none
  !
  integer :: is, na, nt, m1, m2, ldim
  ! counter on spin component
  ! counters on atoms and their type
  ! counters on d components
  integer, parameter :: ldmx = 7
  complex(DP) :: f (ldmx, ldmx), vet (ldmx, ldmx)
  real(DP) :: lambda (ldmx), nsum, nsuma(2)

  WRITE (stdout,*) '--- enter write_ns ---'

  if ( 2 * Hubbard_lmax + 1 > ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

!--
! output of +U parameters
!
  write (stdout,*) 'LDA+U parameters:'
  if (lda_plus_u_kind.eq.0) then
    do nt = 1, ntyp
      if (Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0) then
        if (Hubbard_J0(nt) /= 0.d0 .or. Hubbard_beta(nt) /=0.d0) then
          write (stdout,'(a,i2,a,f12.8)') 'U(',nt,')     =', Hubbard_U(nt)*rytoev 
          write (stdout,'(a,i2,a,f12.8)') 'J0(',nt,')     =', Hubbard_J0(nt)*rytoev 
          write (stdout,'(a,i2,a,f12.8)') 'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
          write (stdout,'(a,i2,a,f12.8)') 'beta(',nt,') =', Hubbard_beta(nt)*rytoev
        else
          write (stdout,'(a,i2,a,f12.8)') 'U(',nt,')     =', Hubbard_U(nt)*rytoev
          write (stdout,'(a,i2,a,f12.8)') 'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
        end if
      endif
    enddo
  else
    do nt = 1, ntyp
      if (Hubbard_U(nt) /= 0.d0) then
        if (Hubbard_l(nt).eq.0) then
          write (stdout,'(a,i2,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
        elseif (Hubbard_l(nt).eq.1) then
         write (stdout,'(2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                              'J(',nt,') =', Hubbard_J(1,nt)*rytoev
        elseif (Hubbard_l(nt).eq.2) then
         write (stdout,'(3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                              'J(',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                              'B(',nt,') =', Hubbard_J(2,nt)*rytoev
        elseif (Hubbard_l(nt).eq.3) then
         write (stdout,'(4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                              'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                              'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                              'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
        endif
      endif
    enddo
  endif
!-- 

  nsum = 0.d0
  do na = 1, nat
     nt = ityp (na)
     if (Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0) then
        ldim = 2 * Hubbard_l(nt) + 1
        nsuma = 0.d0
        do is = 1, nspin
           do m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%ns (m1, m1, is, na)
           enddo
           nsum = nsum + nsuma(is)
        enddo

        if (nspin.eq.1) then
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] = ",f9.5)') &
                              na, 2.d0*nsuma(1)
        else
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] (up, down, total) = ",3f9.5)') &
                              na, nsuma(1), nsuma(2), nsuma(1) + nsuma(2)
        endif 

        do is = 1, nspin
           do m1 = 1, ldim
              do m2 = 1, ldim
                 f (m1, m2) = rho%ns (m1, m2, is, na)
              enddo
           enddo
           call cdiagh(ldim, f, ldmx, lambda, vet)

           if (nspin.ne.1) write( stdout,'("   spin ",i2)') is
           WRITE( stdout,*) '   eigenvalues: '
           WRITE( stdout,'(7f7.3)') (lambda(m1), m1=1, ldim)
  
           WRITE( stdout,*) '   eigenvectors:'
           do m1 = 1, ldim
             WRITE( stdout,'(7f7.3)') ( REAL(vet(m1,m2))**2 + &
                                      AIMAG(vet(m1,m2))**2, m2=1, ldim )
           enddo

           WRITE( stdout,*) '   occupations:'
           do m1 = 1, ldim
             WRITE( stdout,'(7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim )
           enddo
        enddo

        if (nspin.ne.1) write(stdout,'("atomic mag. moment = ",f12.6)') &
                      nsuma(1) - nsuma(2) 
     endif
  enddo

  if (nspin.eq.1) nsum = 2.d0 * nsum 
  WRITE( stdout, '(a,1x,f11.6)') 'N of occupied +U levels =', nsum
  WRITE( stdout,*) '--- exit write_ns ---'
  return
end subroutine write_ns


subroutine write_ns_nc
  !---------------------------------
  ! Noncollinear version (A. Smogunov). 
  !

  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY : npol
  USE io_global,  ONLY : stdout
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_alpha, &
                         Hubbard_U, Hubbard_J
  !
  implicit none
  !
  integer :: is, js, i, na, nt, m1, m2, ldim
  integer, parameter :: ldmx = 7
  complex(DP) :: f (2*ldmx, 2*ldmx), vet (2*ldmx, 2*ldmx)
  real(DP) :: lambda (2*ldmx), nsum, nsuma(2), ns, mx, my, mz

  WRITE (stdout,*) '--- enter write_ns ---'

  if ( 2 * Hubbard_lmax + 1 > ldmx ) &
       call errore ('write_ns', 'ldmx is too small', 1)

!--
! output of +U parameters
!
  write (stdout,*) 'LDA+U parameters:'
  do nt = 1, ntyp
    if (Hubbard_U(nt) /= 0.d0) then
      if (Hubbard_l(nt).eq.0) then
        write (stdout,'(a,i2,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
      elseif (Hubbard_l(nt).eq.1) then
       write (stdout,'(2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                            'J(',nt,') =', Hubbard_J(1,nt)*rytoev
      elseif (Hubbard_l(nt).eq.2) then
       write (stdout,'(3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                            'J(',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                            'B(',nt,') =', Hubbard_J(2,nt)*rytoev
      elseif (Hubbard_l(nt).eq.3) then
       write (stdout,'(4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                            'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                            'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                            'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
      endif
    endif
  enddo
!-- 

  nsum = 0.d0
  do na = 1, nat
     nt = ityp (na)
     if (Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0) then
        ldim = 2 * Hubbard_l(nt) + 1
        nsuma = 0.d0
        do is = 1, npol
           i = is**2
           do m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%ns_nc(m1, m1, i, na)
           end do
        end do
        nsum = nsum + nsuma(1) + nsuma(2) 

        WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] (up, down, total) = ",3f9.5)') &
                              na, nsuma(1), nsuma(2), nsuma(1) + nsuma(2)  

        do m1 = 1, ldim
          do m2 = 1, ldim
            f(m1, m2)           = rho%ns_nc(m1, m2, 1, na)
            f(m1, ldim+m2)      = rho%ns_nc(m1, m2, 2, na)
            f(ldim+m1, m2)      = rho%ns_nc(m1, m2, 3, na)
            f(ldim+m1, ldim+m2) = rho%ns_nc(m1, m2, 4, na)
          enddo
        enddo 

        call cdiagh(2*ldim, f, 2*ldmx, lambda, vet)
        WRITE( stdout,*) 'eigenvalues: '
        WRITE( stdout,'(14f7.3)') (lambda(m1), m1=1, 2*ldim)

        WRITE( stdout,*) 'eigenvectors:'
        do m1 = 1, 2*ldim
          WRITE( stdout,'(14f7.3)') ( REAL(vet(m1,m2))**2 + &
                                      AIMAG(vet(m1,m2))**2, m2=1, 2*ldim )
        enddo

        WRITE( stdout,*) 'occupations, | n_(i1, i2)^(sigma1, sigma2) |:'
        do m1 = 1, 2*ldim
          WRITE( stdout,'(14f7.3)') ( sqrt(REAL(f(m1,m2))**2+ &
                                           AIMAG(f(m1,m2))**2), m2=1, 2*ldim)
        enddo

!-- calculate the spin moment on +U atom 
!
        mx = 0.d0
        my = 0.d0
        mz = 0.d0
        do m1 = 1, 2 * Hubbard_l(nt) + 1
          mx = mx + DBLE( rho%ns_nc(m1, m1, 2, na) + rho%ns_nc(m1, m1, 3, na) )
          my = my + 2.d0 * AIMAG( rho%ns_nc(m1, m1, 2, na) )
          mz = mz + DBLE( rho%ns_nc(m1, m1, 1, na) - rho%ns_nc(m1, m1, 4, na) )
        enddo
        write(stdout,'("atomic mx, my, mz = ",3f12.6)') mx, my, mz
!--

     endif
  enddo

  WRITE( stdout, '(a,1x,f11.6)') 'N of occupied +U levels =', nsum
  WRITE( stdout,*) '--- exit write_ns ---'

  return
end subroutine write_ns_nc

