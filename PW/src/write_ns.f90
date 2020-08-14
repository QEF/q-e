!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_ns
  !-----------------------------------------------------------------------
  !! Prints on output ns infos.
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                         Hubbard_alpha, lda_plus_u_kind, Hubbard_J0, is_hubbard, &
                         Hubbard_beta, Hubbard_U_back, Hubbard_alpha_back, &
                         is_hubbard_back, ldim_back, reserv, reserv_back
  !
  IMPLICIT NONE
  !
  INTEGER :: is, na, nt, m1, m2, ldim
  ! counter on spin component
  ! counters on atoms and their type
  ! counters on d components
  INTEGER, PARAMETER :: ldmx = 10
  COMPLEX(DP) :: f(ldmx,ldmx), vet(ldmx,ldmx)
  REAL(DP) :: lambda(ldmx), nsum, nsuma(2), rsrv
  !
  WRITE (stdout,*) '--- enter write_ns ---'
  !
  IF ( 2*Hubbard_lmax+1 > ldmx ) &
           CALL errore( 'write_ns', 'ldmx is too small', 1 )
  !
  ! ... output of +U parameters
  !
  WRITE (stdout,*) 'LDA+U parameters:'
  !
  IF (lda_plus_u_kind == 0) THEN
     !
     DO nt = 1, ntyp
        IF (is_hubbard(nt)) THEN
          IF (Hubbard_U(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)')          & 
                       'U(',nt,')     =', Hubbard_U(nt)*rytoev
          IF (Hubbard_J0(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)')         &
                       'J0(',nt,')     =', Hubbard_J0(nt)*rytoev
          IF (Hubbard_alpha(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)')      &
                       'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
          IF (Hubbard_beta(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)')       &
                       'beta(',nt,') =', Hubbard_beta(nt)*rytoev
        ENDIF
        IF (is_hubbard_back(nt)) THEN
          IF (Hubbard_U_back(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)')     &
                       'U_back(',nt,')     =', Hubbard_U_back(nt)*rytoev
          IF (Hubbard_alpha_back(nt) /= 0.d0) WRITE(stdout,'(a,i2,a,f12.8)') &
                       'alpha_back(',nt,') =', Hubbard_alpha_back(nt)*rytoev
      ENDIF
     ENDDO
     !
  ELSEIF (lda_plus_u_kind == 1) THEN
     !
     DO nt = 1, ntyp
        IF (Hubbard_U(nt) /= 0.d0) THEN
           IF (Hubbard_l(nt) == 0) THEN
              WRITE(stdout,'(a,i2,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
           ELSEIF (Hubbard_l(nt) == 1) THEN
              WRITE(stdout,'(2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,    &
                                                  'J(',nt,') =', Hubbard_J(1,nt)*rytoev
           ELSEIF (Hubbard_l(nt) == 2) THEN
              WRITE(stdout,'(3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,    &
                                                  'J(',nt,') =', Hubbard_J(1,nt)*rytoev,  &
                                                  'B(',nt,') =', Hubbard_J(2,nt)*rytoev
           ELSEIF (Hubbard_l(nt) == 3) THEN
              WRITE(stdout,'(4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                                  'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                                  'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                                  'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
           ENDIF
        ENDIF
     ENDDO
     !
  ENDIF
  ! 
  nsum = 0.d0
  rsrv = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (is_hubbard(nt)) THEN
        !
        ldim = 2 * Hubbard_l(nt) + 1
        nsuma = 0.d0
        !
        DO is = 1, nspin
           DO m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%ns(m1,m1,is,na)
              IF (reserv(nt)) rsrv = rsrv + rho%ns (m1, m1, is, na)
           ENDDO
           nsum = nsum + nsuma(is)
        ENDDO
        !
        ! Write to the output file
        !
        IF (nspin==1) THEN
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] = ",f9.5)') &
                              na, 2.d0*nsuma(1)
        ELSE
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] (up, down, total) = ",3f9.5)') &
                              na, nsuma(1), nsuma(2), nsuma(1)+nsuma(2)
        ENDIF 
        !
        ! Diagonalize the occupation matrix
        !
        DO is = 1, nspin
           DO m1 = 1, ldim
              DO m2 = 1, ldim
                 f(m1,m2) = rho%ns(m1,m2,is,na)
              ENDDO
           ENDDO
           !
           CALL cdiagh( ldim, f, ldmx, lambda, vet )
           !
           IF (nspin /= 1) write( stdout,'("   spin ",i2)') is
           WRITE( stdout,*) '   eigenvalues: '
           WRITE( stdout,'(7f7.3)') (lambda(m1), m1=1, ldim)
           !
           WRITE( stdout,*) '   eigenvectors:'
           DO m1 = 1, ldim
              WRITE( stdout,'(7f7.3)') ( REAL(vet(m1,m2))**2 + &
                                         AIMAG(vet(m1,m2))**2, m2=1, ldim )
           ENDDO
           !
           WRITE( stdout,*) '   occupations:'
           DO m1 = 1, ldim
              WRITE( stdout,'(7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim )
           ENDDO
        ENDDO
        !
        IF (nspin /= 1) WRITE( stdout,'("atomic mag. moment = ",f9.5)') &
                               nsuma(1) - nsuma(2) 
     ENDIF
     !
     ! Background part
     !
     IF (is_hubbard_back(nt)) THEN
        !
        ldim = ldim_back(nt) !2 * Hubbard_l_back(nt) + 1
        nsuma = 0.d0
        !
        DO is = 1, nspin
           DO m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%nsb (m1, m1, is, na)
              IF (reserv_back(nt)) rsrv = rsrv + rho%nsb (m1, m1, is, na)
           ENDDO
           nsum = nsum + nsuma(is)
        ENDDO
        !
        ! Write to the output file
        !
        IF (nspin.EQ.1) THEN
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] = ",f9.5)') &
                              na, 2.d0*nsuma(1)
        ELSE
           WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] (up, down, total) = ",3f9.5)') &
                              na, nsuma(1), nsuma(2), nsuma(1) + nsuma(2)
        ENDIF
        !
        ! Diagonalize
        !
        DO is = 1, nspin
           !
           DO m1 = 1, ldim
              DO m2 = 1, ldim
                 f (m1, m2) = rho%nsb (m1, m2, is, na)
              ENDDO
           ENDDO
           !
           CALL cdiagh(ldim, f, ldmx, lambda, vet)
           !
           IF (nspin.NE.1) WRITE( stdout,'("   spin ",i2)') is
           WRITE( stdout,*) '   eigenvalues: '
           WRITE( stdout,'(7f7.3)') (lambda(m1), m1=1, ldim)
           !
           WRITE( stdout,*) '   eigenvectors:'
           DO m1 = 1, ldim
             WRITE( stdout,'(7f7.3)') ( DBLE(vet(m1,m2))**2 + &
                                      AIMAG(vet(m1,m2))**2, m2=1, ldim )
           ENDDO
           !
           WRITE( stdout,*) '   occupations:'
           DO m1 = 1, ldim
             WRITE( stdout,'(7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim )
           ENDDO
           !
        ENDDO
        !
        IF (nspin.NE.1) WRITE(stdout,'(''atomic mag. moment = '',f9.5)') &
                        nsuma(1) - nsuma(2)
        !
     ENDIF
     !
  ENDDO ! na
  !
  !
  IF (nspin==1) nsum = 2.d0 * nsum 
  !
  WRITE( stdout, '(a,1x,f11.6)') 'N of occupied +U levels =', nsum
  !
  IF (rsrv.GT.0.d0) &
     WRITE(stdout,'("Total occupation of reservoir states = ",x,f11.6)') rsrv
  !
  WRITE( stdout,*) '--- exit write_ns ---'
  !
  RETURN
  !
END SUBROUTINE write_ns
!
!-----------------------------------------------------------------------
SUBROUTINE write_ns_nc
  !-----------------------------------------------------------------------
  !! Prints on output ns infos. Noncollinear version (A. Smogunov).
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : rytoev
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp
  USE noncollin_module,  ONLY : npol
  USE io_global,         ONLY : stdout
  USE scf,               ONLY : rho
  USE ldaU,              ONLY : Hubbard_lmax, Hubbard_l, Hubbard_alpha, &
                                Hubbard_U, Hubbard_J
  !
  IMPLICIT NONE
  !
  INTEGER :: is, js, i, na, nt, m1, m2, ldim
  INTEGER, PARAMETER :: ldmx = 7
  COMPLEX(DP) :: f(2*ldmx,2*ldmx), vet(2*ldmx,2*ldmx)
  REAL(DP) :: lambda(2*ldmx), nsum,nsuma(2), ns, mx, my, mz
  !
  WRITE (stdout,*) '--- enter write_ns ---'
  !
  IF ( 2 * Hubbard_lmax + 1 > ldmx ) &
       CALL errore( 'write_ns', 'ldmx is too small', 1 )
  !
  ! ... output of +U parameters
  !
  WRITE (stdout,*) 'LDA+U parameters:'
  !
  DO nt = 1, ntyp
     IF (Hubbard_U(nt) /= 0.d0) THEN
        IF (Hubbard_l(nt)==0) THEN
           WRITE(stdout,'(a,i2,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
        ELSEIF (Hubbard_l(nt)==1) THEN
           WRITE(stdout,'(2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev, &
                                               'J(',nt,') =', Hubbard_J(1,nt)*rytoev
        ELSEIF (Hubbard_l(nt)==2) THEN
           WRITE(stdout,'(3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                               'J(',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                               'B(',nt,') =', Hubbard_J(2,nt)*rytoev
        ELSEIF (Hubbard_l(nt)==3) THEN
           WRITE(stdout,'(4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                               'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                               'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                               'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
        ENDIF
     ENDIF
  ENDDO
  !
  !
  nsum = 0.d0
  DO na = 1, nat
     nt = ityp (na)
     IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
        ldim = 2 * Hubbard_l(nt) + 1
        nsuma = 0.d0
        DO is = 1, npol
           i = is**2
           DO m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%ns_nc(m1,m1,i,na)
           ENDDO
        ENDDO
        nsum = nsum + nsuma(1) + nsuma(2) 
        !
        WRITE( stdout,'("atom ",i4,3x,"Tr[ns(na)] (up, down, total) = ",3f9.5)') &
                              na, nsuma(1), nsuma(2), nsuma(1) + nsuma(2)  
        !
        DO m1 = 1, ldim
           DO m2 = 1, ldim
              f(m1, m2)           = rho%ns_nc(m1,m2,1,na)
              f(m1, ldim+m2)      = rho%ns_nc(m1,m2,2,na)
              f(ldim+m1, m2)      = rho%ns_nc(m1,m2,3,na)
              f(ldim+m1, ldim+m2) = rho%ns_nc(m1,m2,4,na)
           ENDDO
        ENDDO 
        !
        CALL cdiagh( 2*ldim, f, 2*ldmx, lambda, vet )
        !
        WRITE( stdout,*) 'eigenvalues: '
        WRITE( stdout,'(14f7.3)') (lambda(m1), m1=1, 2*ldim)
        !
        WRITE( stdout,*) 'eigenvectors:'
        !
        DO m1 = 1, 2*ldim
          WRITE( stdout,'(14f7.3)') ( REAL(vet(m1,m2))**2 + &
                                      AIMAG(vet(m1,m2))**2, m2=1, 2*ldim )
        ENDDO
        !
        WRITE( stdout,*) 'occupations, | n_(i1, i2)^(sigma1, sigma2) |:'
        !
        DO m1 = 1, 2*ldim
           WRITE( stdout,'(14f7.3)') ( SQRT(REAL(f(m1,m2))**2 + &
                                       AIMAG(f(m1,m2))**2), m2=1, 2*ldim)
        ENDDO
        !
        ! ... calculate the spin moment on +U atom 
        !
        mx = 0.d0
        my = 0.d0
        mz = 0.d0
        DO m1 = 1, 2 * Hubbard_l(nt) + 1
          mx = mx + DBLE( rho%ns_nc(m1,m1,2,na) + rho%ns_nc(m1,m1,3,na) )
          my = my + 2.d0 * AIMAG( rho%ns_nc(m1,m1,2,na) )
          mz = mz + DBLE( rho%ns_nc(m1,m1,1,na) - rho%ns_nc(m1,m1,4,na) )
        ENDDO
        WRITE(stdout,'("atomic mx, my, mz = ",3f12.6)') mx, my, mz
        !
     ENDIF
  ENDDO
  !
  WRITE( stdout, '(a,1x,f11.6)') 'N of occupied +U levels =', nsum
  WRITE( stdout,*) '--- exit write_ns ---'
  !
  RETURN
  !
END SUBROUTINE write_ns_nc

!-----------------------------------------------------------------------
SUBROUTINE write_nsg 
  !-----------------------------------------------------------------------
  !
  ! Generalized ns (i.e. nsg) for DFT+U+V
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : rytoev
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE io_global,  ONLY : stdout
  USE scf,        ONLY : rho
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                         Hubbard_alpha, lda_plus_u_kind, Hubbard_J0, &
                         Hubbard_beta, iso_sys, is_hubbard, is_hubbard_back, &
                         ldim_u, reserv, reserv_back, Hubbard_V, &
                         at_sc, neighood, Hubbard_U_back,Hubbard_alpha_back, &
                         nsgnew, backall, ldim_back
  !
  IMPLICIT NONE
  INTEGER :: is, na, nt, m1, m2, ldim
  ! counter on spin component
  ! counters on atoms and their type
  ! counters on d components
  INTEGER :: ldmx, ldmax
  INTEGER :: na1, na2, ldm1, ldm2, ldim1,ldim2, viz, nt1, nr, nt2
  !counter on atom, 
  !keeping track of hubbard space dimension
  !number of the neighbor atom
  !type
  !number of reservoir states
  COMPLEX(DP), ALLOCATABLE :: f(:,:) , vet(:,:)
  REAL(DP) :: nsum, nsuma(2)
  REAL(DP) :: totsum, sumr, sumr_b, norm
  !total occupations, reservoir and back res.
  COMPLEX(DP) :: aux !temp
  REAL(DP), ALLOCATABLE :: lambda(:)
  ! 
  WRITE (stdout,*) '--- enter write_nsg ---'
  !
  ldmx  = 0
  ldmax = 0
  DO na = 1, nat
     nt = ityp(na)
     ldmx = ldmx + ldim_u(nt)
     ldmax = MAX(ldmax,ldim_u(nt))
  ENDDO
  !
  IF ( 2*Hubbard_lmax+1 > ldmx .OR. &
       2*Hubbard_lmax+1 > ldmax ) &
       CALL errore ('write_nsg', 'ldmx is too small', 1)
  !
  IF (iso_sys) THEN 
     ALLOCATE (f(ldmx,ldmx),vet(ldmx,ldmx),lambda(ldmx))
  ELSE
     ALLOCATE (f(ldmax,ldmax),vet(ldmax,ldmax),lambda(ldmax))
  ENDIF
  !
  DO nt = 1, ntyp
     IF (Hubbard_alpha(nt) /= 0.d0 .OR. Hubbard_alpha_back(nt) /= 0.d0 .OR. &
         Hubbard_J0(nt) /= 0.d0) THEN
        WRITE (stdout,*) 'DFT+U+V parameters:'
        WRITE (stdout,'(a,i2,a,f12.8)') 'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
        WRITE (stdout,'(a,i2,a,f12.8)') 'alpha_back(',nt,') =', Hubbard_alpha_back(nt)*rytoev
        WRITE (stdout,'(a,i2,a,f12.8)') 'J0(',nt,') =', Hubbard_J0(nt)*rytoev
     ENDIF
  ENDDO
  !
  ! Calculation of occupations
  !
  DO is = 1, nspin
     !
     ! Construct the occupation matrix to be diagonalized
     !
     IF (iso_sys) THEN
        !
        ! This is the case for isolated systems
        !
        ldm1 = 0
        !
        DO na1 = 1, nat
           !
           nt1 = ityp(na1)
           !
           IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
              !
              DO m1 = 1, ldim_u(nt1)
                 !
                 ldm2 = 0
                 !
                 DO na2 = 1, nat
                    ! 
                    nt2 = ityp(na2)
                    !
                    IF ( is_hubbard(nt2) .OR. is_hubbard_back(nt2) ) THEN
                       !
                       DO m2 = 1, ldim_u(nt2)
                          f(ldm1+m1,ldm2+m2) = (0.d0, 0.d0)
                       ENDDO
                       !
                       DO viz = 1, neighood(na1)%num_neigh
                          IF (neighood(na1)%neigh(viz).EQ.na2) THEN
                             DO m2 = 1, ldim2
                                f(ldm1+m1,ldim2+m2) = nsgnew(m2,m1,viz,na1,is)
                             ENDDO
                             GO TO 3
                          ENDIF
                       ENDDO
                       !
3                      CONTINUE
                       !
                       ldm2 = ldm2 + ldim_u(nt2)
                       !
                    ENDIF
                    !
                 ENDDO !na2
                 ! 
              ENDDO !m1
              !
              ldm1 = ldm1 + ldim_u(nt1)
              ! 
           ENDIF
           !
        ENDDO !na1 
        !
        ! Diagonalize the occupation matrix, we need the eigenvalues
        !
        CALL cdiagh(ldmx, f, ldmx, lambda, vet)
        !
        ! Report the eigenvalues
        !
        nsum = 0.d0
        DO m1 = 1, ldmx   
           !
           WRITE(stdout,'(a6,x,i2,x,a8,x,i3,x,a13,x,f8.4)')    &
                'Spin: ', is, ' State: ', m1, ' Occupation: ', lambda(m1) 
           !
           nsum = nsum + lambda(m1)
           ldm1 = 0 
           !
           DO na1 = 1, nat
              !
              nt1 = ityp(na1)
              !
              IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
                 !
                 ldim1 = ldim_u(nt1)
                 norm = 0.d0
                 !
                 DO m2 = 1, ldim1
                    norm = norm + DBLE(vet(ldm1+m2,m1))**2 + AIMAG(vet(ldm1+m2,m1))**2
                 ENDDO
                 !
                 IF (norm.GE.1.d-5) THEN
                    WRITE(stdout,'(a7,x,i3,x,a10,x,i2)') ' Atom: ', na1, ' n States: ', ldim1
                    WRITE(stdout,*) 'norm of eigenvector = ', norm
                    WRITE(stdout,'(a7,9(x,f10.6))') &
                         'real: ', (DBLE(vet(ldm1+m2,m1)),m2=1,ldim1)
                    WRITE(stdout,*)
                 ENDIF
                 !
                 ldm1 = ldm1 + ldim1
                 !
              ENDIF
              !
           ENDDO !na1
           !
        ENDDO !m1
        WRITE(stdout,*) 'Sum = ', nsum
        !
     ELSE
        !
        ! This is the case for solids
        !
        DO na1 = 1, nat
           !
           f(:,:) = (0.d0, 0.d0)
           !
           nt1 = ityp(na1)
           !
           IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
              !
              ldim1 = ldim_u(nt1)
              !
              DO viz = 1, neighood(na1)%num_neigh
                 na2 = neighood(na1)%neigh(viz)
                 IF (na2.EQ.na1) THEN
                    DO m1 = 1, ldim1
                       DO m2 = 1, ldim1
                          f(m1,m2) = nsgnew(m2,m1,viz,na1,is)
                       ENDDO
                    ENDDO
                    GO TO 4
                 ENDIF
              ENDDO
              !
4             CONTINUE
              !
              ! This will only print the on-site blocks of the matrix. 
              ! The diagonalization will not give components on states of other atoms. 
              ! To be improved for periodic systems.
              !
              CALL cdiagh(ldim1, f, ldmax, lambda, vet)
              !
              WRITE(stdout,'(a7,x,i3,x,a8,x,i1)') ' Atom: ', na1, ' Spin: ', is
              !
              WRITE(stdout,*) '   eigenvalues and eigenvectors of the occupation matrix:'
              DO m1 = 1, ldim1
                WRITE( stdout,'(f7.3)') lambda(m1)
                WRITE( stdout,'(9f7.3)') ( DBLE(vet(m2,m1)), m2=1, ldim1 )
              ENDDO
              !
              WRITE(stdout,*) '   occupation matrix before diagonalization:'
              DO m1 = 1, ldim1
                 WRITE( stdout,'(9f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim1 )
              ENDDO
              !
           ENDIF
           !
        ENDDO !na1
        !
     ENDIF !iso_sys
     !
  ENDDO !is 
  !
  ! THE CASE OF THE RESERVOIR  
  !
  nr = 0
  sumr = 0.d0
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        !
        DO viz = 1, neighood(na1)%num_neigh
           !
           na2 = neighood(na1)%neigh(viz)
           !
           IF (na2.EQ.na1) THEN
              !
              totsum = 0.d0
              DO is = 1, nspin
                 DO m1 = 1, 2*Hubbard_l(nt1)+1 
                    totsum = totsum + DBLE(nsgnew(m1,m1,viz,na1,is))
                 ENDDO
              ENDDO
              !
              IF (nspin.EQ.1) totsum = 2.d0 * totsum
              !
              IF (reserv(nt1)) THEN
                 nr = nr + 1               
                 sumr = sumr + totsum
              ENDIF
              !
              WRITE( stdout,'(a,1x,i4,2x,a,f11.7)') 'atom', na1, ' Tr[ns(na)]= ', totsum
              !
              IF (nspin.EQ.2) THEN
                 totsum = 0.d0
                 DO m1 = 1, 2*Hubbard_l(nt1)+1
                    totsum = totsum + DBLE(nsgnew(m1,m1,viz,na1,1)-nsgnew(m1,m1,viz,na1,2))
                 ENDDO
                 WRITE( stdout,'(a,1x,i4,2x,a,f11.7)') 'atom', na1, ' Mag[ns(na)]= ', totsum
              ENDIF
              !
           ENDIF !na1 = na2
           !
        ENDDO !viz
        !
     ENDIF
     !
  ENDDO !na1
  !
  ! THE CASE OF THE BACKGROUND STATES WITH V: RESERV_BACK
  !
  nr = 0
  sumr_b = 0.d0
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard_back(nt1) ) THEN 
        !
        DO viz = 1, neighood(na1)%num_neigh
           !
           na2 = neighood(na1)%neigh(viz)
           !
           IF ( na2.EQ.na1 ) THEN
              !
              totsum = 0.d0
              DO is = 1, nspin
                 DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                    totsum = totsum + DBLE(nsgnew(m1,m1,viz,na1,is))
                 ENDDO
              ENDDO
              !
              IF (nspin.EQ.1) totsum = 2.d0 * totsum
              !
              IF (reserv_back(nt1)) THEN
                 nr = nr + 1
                 sumr_b = sumr_b + totsum
              ENDIF
              !
              WRITE( stdout,'(a,1x,i4,2x,a,f11.7)') 'atom', na1, ' Tr[ns(na)]= ', totsum
              !  
              IF (nspin.EQ.2) THEN
                 totsum = 0.d0
                 DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                    totsum = totsum + DBLE(nsgnew(m1,m1,viz,na1,1)-nsgnew(m1,m1,viz,na1,2))
                 ENDDO
                 WRITE( stdout,'(a,1x,i4,2x,a,f11.7)') 'atom', na1, ' Mag[ns(na)]= ', totsum
              ENDIF
              !
           ENDIF !na1=na2
           !
        ENDDO !viz
        !
     ENDIF
     !
  ENDDO !na1
  !
  IF (nr.GT.0) WRITE(stdout, '(a,1x,f11.7)') 'Occ Reservoir =', sumr + sumr_b    
  !
  DEALLOCATE(f,vet,lambda) 
  WRITE( stdout,*) '--- exit write_nsg ---'
  !
  RETURN
  !
END SUBROUTINE write_nsg 

