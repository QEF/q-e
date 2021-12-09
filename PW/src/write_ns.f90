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
  USE ldaU,       ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, ldim_u,  &
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
  COMPLEX(DP), ALLOCATABLE :: f(:,:) , vet(:,:)
  REAL(DP), ALLOCATABLE :: lambda(:)
  REAL(DP) :: nsum, nsuma(2), rsrv
  !
  ! Print Hubbard parameters
  !
  WRITE (stdout,'(/5x,"Hubbard parameters (eV):")')
  !
  IF (lda_plus_u_kind == 0) THEN
     !
     DO nt = 1, ntyp
        IF (is_hubbard(nt)) THEN
          IF (Hubbard_U(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)')          & 
                       'U(',nt,') =', Hubbard_U(nt)*rytoev
          IF (Hubbard_J0(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)')         &
                       'J0(',nt,') =', Hubbard_J0(nt)*rytoev 
          IF (Hubbard_alpha(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)')      &
                       'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
          IF (Hubbard_beta(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)')       &
                       'beta(',nt,') =', Hubbard_beta(nt)*rytoev
        ENDIF
        IF (is_hubbard_back(nt)) THEN
          IF (Hubbard_U_back(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)')     &
                       'U_back(',nt,') =', Hubbard_U_back(nt)*rytoev
          IF (Hubbard_alpha_back(nt) /= 0.d0) WRITE(stdout,'(5x,a,i3,a,f8.4)') &
                       'alpha_back(',nt,') =', Hubbard_alpha_back(nt)*rytoev
      ENDIF
     ENDDO
     !
  ELSEIF (lda_plus_u_kind == 1) THEN
     !
     DO nt = 1, ntyp
        IF (Hubbard_U(nt) /= 0.d0) THEN
           IF (Hubbard_l(nt) == 0) THEN
              WRITE(stdout,'(5x,a,i3,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
           ELSEIF (Hubbard_l(nt) == 1) THEN
              WRITE(stdout,'(5x,2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,    &
                                                  'J(',nt,') =', Hubbard_J(1,nt)*rytoev
           ELSEIF (Hubbard_l(nt) == 2) THEN
              WRITE(stdout,'(5x,3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,    &
                                                     'J(',nt,') =', Hubbard_J(1,nt)*rytoev,  &
                                                     'B(',nt,') =', Hubbard_J(2,nt)*rytoev
           ELSEIF (Hubbard_l(nt) == 3) THEN
              WRITE(stdout,'(5x,4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                                     'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                                     'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                                     'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
           ENDIF
        ENDIF
     ENDDO
     !
  ENDIF
  ! 
  WRITE (stdout,'(/5x,17("="), " HUBBARD OCCUPATIONS ",16("="))')
  !
  nsum = 0.d0
  rsrv = 0.d0
  !
  ! Calculation of occupations
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF (is_hubbard(nt)) THEN
        !
        ldim = 2 * Hubbard_l(nt) + 1
        !
        ! Compute the trace of the occupation matrix and the magnetization
        !
        WRITE( stdout,'(5x,21("-")," ATOM ",i4,1x,22("-"))') na
        !
        nsuma = 0.d0
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
           WRITE( stdout,'(5x,"Tr[ns] = ",f9.5)') 2.d0*nsuma(1)
        ELSE
           WRITE( stdout,'(5x,"Tr[ns] (up, down, total) = ",3f9.5)') &
                              nsuma(1), nsuma(2), nsuma(1)+nsuma(2)
           WRITE( stdout,'(5x,"Atomic magnetic moment = ",f9.5)') nsuma(1) - nsuma(2)
        ENDIF 
        !
        ALLOCATE (f(ldim,ldim), vet(ldim,ldim), lambda(ldim))
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
           CALL cdiagh( ldim, f, ldim, lambda, vet )
           !
           IF (nspin /= 1) WRITE( stdout,'(5x,"SPIN ",i2)') is
           WRITE( stdout,'(5x,"eigenvalues:")')
           WRITE( stdout,'(5x,7f7.3)') (lambda(m1), m1=1, ldim)
           !
           WRITE( stdout,'(5x,"eigenvectors (columns):")')
           DO m1 = 1, ldim
              WRITE( stdout,'(5x,7f7.3)') ( DBLE(vet(m1,m2)), m2=1, ldim )
           ENDDO
           !
           WRITE( stdout,'(5x,"occupation matrix ns (before diag.):")')
           DO m1 = 1, ldim
              WRITE( stdout,'(5x,7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim )
           ENDDO
        ENDDO
        !
        DEALLOCATE (f, vet, lambda)
        !
     ENDIF
     !
     ! Background part
     !
     IF (is_hubbard_back(nt)) THEN
        !
        ldim = ldim_back(nt) !2 * Hubbard_l_back(nt) + 1
        !
        WRITE( stdout,'(5x," Background part ")')
        !
        nsuma = 0.d0
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
           WRITE( stdout,'(5x,"Tr[ns] = ",f9.5)') 2.d0*nsuma(1)
        ELSE
           WRITE( stdout,'(5x,"Tr[ns] (up, down, total) = ",3f9.5)') &
                              nsuma(1), nsuma(2), nsuma(1) + nsuma(2)
           WRITE(stdout,'(5x,"Atomic magnetic moment = ",f9.5)') nsuma(1) - nsuma(2)
        ENDIF
        !
        ALLOCATE (f(ldim,ldim), vet(ldim,ldim), lambda(ldim))
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
           CALL cdiagh(ldim, f, ldim, lambda, vet)
           !
           IF (nspin /= 1) WRITE( stdout,'(5x,"SPIN ",i2)') is
           WRITE( stdout,'(5x,"eigenvalues:")')
           WRITE( stdout,'(5x,7f7.3)') (lambda(m1), m1=1, ldim)
           !
           WRITE( stdout,'(5x,"eigenvectors (columns):")')
           DO m1 = 1, ldim
             WRITE( stdout,'(5x,7f7.3)') ( DBLE(vet(m1,m2)), m2=1, ldim )
           ENDDO
           !
           WRITE( stdout,'(5x,"occupation matrix ns (before diag.):")')
           DO m1 = 1, ldim
             WRITE( stdout,'(5x,7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim )
           ENDDO
           !
        ENDDO
        !
        DEALLOCATE (f, vet, lambda)
        !
     ENDIF
     !
  ENDDO ! na
  !
  !
  IF (nspin==1) nsum = 2.d0 * nsum 
  !
  WRITE( stdout, '(/5x,a,1x,f9.4)') 'Number of occupied Hubbard levels =', nsum
  !
  IF (rsrv.GT.0.d0) &
     WRITE(stdout,'(5x,"Total occupation of reservoir states = ",x,f11.6)') rsrv
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
  COMPLEX(DP), ALLOCATABLE :: f(:,:) , vet(:,:)
  REAL(DP), ALLOCATABLE :: lambda(:)
  REAL(DP) :: nsum,nsuma(2), ns, mx, my, mz
  !
  ! Print Hubbard parameters
  !
  WRITE (stdout,'(/5x,"Hubbard parameters (eV):")')
  !
  DO nt = 1, ntyp
     IF (Hubbard_U(nt) /= 0.d0) THEN
        IF (Hubbard_l(nt)==0) THEN
           WRITE(stdout,'(5x,a,i3,a,f12.8)') 'U(',nt,') =', Hubbard_U(nt) * rytoev
        ELSEIF (Hubbard_l(nt)==1) THEN
           WRITE(stdout,'(5x,2(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev, &
                                                  'J(',nt,') =', Hubbard_J(1,nt)*rytoev
        ELSEIF (Hubbard_l(nt)==2) THEN
           WRITE(stdout,'(5x,3(a,i3,a,f9.4,3x))') 'U(',nt,') =', Hubbard_U(nt)*rytoev,   &
                                                  'J(',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                                  'B(',nt,') =', Hubbard_J(2,nt)*rytoev
        ELSEIF (Hubbard_l(nt)==3) THEN
           WRITE(stdout,'(5x,4(a,i3,a,f9.4,3x))') 'U (',nt,') =', Hubbard_U(nt)*rytoev,   &
                                                  'J (',nt,') =', Hubbard_J(1,nt)*rytoev, &
                                                  'E2(',nt,') =', Hubbard_J(2,nt)*rytoev, &
                                                  'E3(',nt,') =', Hubbard_J(3,nt)*rytoev
        ENDIF
     ENDIF
  ENDDO
  !
  WRITE (stdout,'(/5x,17("="), " HUBBARD OCCUPATIONS ",16("="))')
  !
  ! Calculation of occupations
  !
  nsum = 0.d0
  DO na = 1, nat
     ! 
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
        !     
        ldim = 2 * Hubbard_l(nt) + 1
        !
        ! Compute the trace of the occupation matrix and the magnetization
        !
        WRITE( stdout,'(5x,21("-")," ATOM ",i4,1x,22("-"))') na
        !
        nsuma = 0.d0
        DO is = 1, npol
           i = is**2
           DO m1 = 1, ldim
              nsuma(is) = nsuma(is) + rho%ns_nc(m1,m1,i,na)
           ENDDO
        ENDDO
        nsum = nsum + nsuma(1) + nsuma(2) 
        !
        WRITE( stdout,'(5x,"Tr[ns] (up, down, total) = ",3f9.5)') &
                              nsuma(1), nsuma(2), nsuma(1) + nsuma(2)  
        !
        ALLOCATE (f(2*ldim,2*ldim), vet(2*ldim,2*ldim), lambda(2*ldim))
        !
        ! Diagonalize the occupation matrix
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
        CALL cdiagh( 2*ldim, f, 2*ldim, lambda, vet )
        !
        WRITE( stdout,'(5x,"eigenvalues:")')
        WRITE( stdout,'(5x,14f7.3)') (lambda(m1), m1=1, 2*ldim)
        !
        WRITE( stdout,'(5x,"eigenvectors (columns):")')
        DO m1 = 1, 2*ldim
          WRITE( stdout,'(5x,14f7.3)') ( DBLE(vet(m1,m2)), m2=1, 2*ldim )
        ENDDO
        !
        WRITE( stdout,'(5x,"occupations, | n_(i1, i2)^(sigma1, sigma2) |:")')
        ! IT: Should be enough to print only DBLE(f(m1,m2))
        DO m1 = 1, 2*ldim
           WRITE( stdout,'(5x,14f7.3)') ( DSQRT(DBLE(f(m1,m2))**2 + &
                                          AIMAG(f(m1,m2))**2), m2=1, 2*ldim)
        ENDDO
        !
        DEALLOCATE (f, vet, lambda)
        !
        ! ... calculate the spin moment on +U atom 
        !
        mx = 0.d0
        my = 0.d0
        mz = 0.d0
        DO m1 = 1, ldim
          mx = mx + DBLE( rho%ns_nc(m1,m1,2,na) + rho%ns_nc(m1,m1,3,na) )
          my = my + 2.d0 * AIMAG( rho%ns_nc(m1,m1,2,na) )
          mz = mz + DBLE( rho%ns_nc(m1,m1,1,na) - rho%ns_nc(m1,m1,4,na) )
        ENDDO
        WRITE(stdout,'(5x,"Atomic magnetic moment mx, my, mz = ",3f12.6)') mx, my, mz
        !
     ENDIF
  ENDDO
  !
  WRITE( stdout, '(/5x,a,1x,f9.4)') 'Number of occupied Hubbard levels =', nsum
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
  INTEGER :: na1, na2, ldm1, ldm2, ldim1,ldim2, viz, nt1, nr, nt2, ldmx
  !counter on atom, 
  !keeping track of hubbard space dimension
  !number of the neighbor atom
  !type
  !number of reservoir states
  COMPLEX(DP), ALLOCATABLE :: f(:,:) , vet(:,:)
  REAL(DP), ALLOCATABLE :: lambda(:)
  REAL(DP) :: nsum, nsum_iso, nsuma(2), rsrv, norm
  ! 
  !
  IF (iso_sys) THEN 
     ldmx  = 0
     DO na = 1, nat
        nt = ityp(na)
        ldmx = ldmx + ldim_u(nt)
     ENDDO
     IF ( 2*Hubbard_lmax+1 > ldmx ) &
       CALL errore ('write_nsg', 'ldmx is too small', 1)
     ALLOCATE (f(ldmx,ldmx), vet(ldmx,ldmx), lambda(ldmx))
  ENDIF
  !
  DO nt = 1, ntyp
     IF (Hubbard_alpha(nt) /= 0.d0 .OR. Hubbard_alpha_back(nt) /= 0.d0 .OR. &
         Hubbard_J0(nt) /= 0.d0) THEN
        WRITE (stdout,'(/5x,"Hubbard parameters (eV):")')
        WRITE (stdout,'(5x,a,i3,a,f12.8)') 'alpha(',nt,') =', Hubbard_alpha(nt)*rytoev
        WRITE (stdout,'(5x,a,i3,a,f12.8)') 'alpha_back(',nt,') =', Hubbard_alpha_back(nt)*rytoev
        WRITE (stdout,'(5x,a,i3,a,f12.8)') 'J0(',nt,') =', Hubbard_J0(nt)*rytoev
     ENDIF
  ENDDO
  !
  WRITE (stdout,'(/5x,17("="), " HUBBARD OCCUPATIONS ",16("="))')
  !
  ! Construct the occupation matrix to be diagonalized
  !
  IF (iso_sys) THEN
     !
     ! This is the case for isolated systems
     !
     DO is = 1, nspin
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
                             DO m2 = 1, ldim_u(nt2) ! TODO: Check 
                                f(ldm1+m1,ldm2+m2) = nsgnew(m2,m1,viz,na1,is)
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
        nsum_iso = 0.d0
        DO m1 = 1, ldmx   
           !
           WRITE(stdout,'(5x,a6,x,i2,x,a8,x,i3,x,a13,x,f8.4)')    &
                'Spin: ', is, ' State: ', m1, ' Occupation: ', lambda(m1) 
           !
           nsum_iso = nsum_iso + lambda(m1)
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
                    WRITE(stdout,'(5x,a7,x,i3,x,a10,x,i2)') ' Atom: ', na1, ' n States: ', ldim1
                    WRITE(stdout,'(5x,"norm of eigenvector =",1x,f10.6)') norm
                    WRITE(stdout,'(5x,a7,9(x,f10.6))') &
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
        WRITE(stdout,'(5x,"Sum of eigenvalues =",1x,f10.6)') nsum_iso
        !
     ENDDO ! is
     !
  ELSE
     !
     ! This is the case for solids
     ! Note: this case is adapted to print the same information and in the same format
     ! as in the DFT+U case (see write_ns)
     !
     nsum = 0.d0
     rsrv = 0.d0
     !
     DO na1 = 1, nat
        !
        nt1 = ityp(na1)
        !
        IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
           !
           ldim1 = 2*Hubbard_l(nt1)+1
           !
           ! Compute the trace of the occupation matrix and the magnetization
           !
           WRITE( stdout,'(5x,21("-")," ATOM ",i4,1x,22("-"))') na1
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              !
              IF (na2.EQ.na1) THEN
                 !
                 nsuma = 0.d0
                 DO is = 1, nspin
                    DO m1 = 1, ldim1
                       nsuma(is) = nsuma(is) + DBLE(nsgnew(m1,m1,viz,na1,is))
                       IF (reserv(nt1)) rsrv = rsrv + DBLE(nsgnew(m1,m1,viz,na1,is))
                    ENDDO
                    nsum = nsum + nsuma(is)
                 ENDDO
                 !
                 IF (nspin.EQ.1) THEN
                    WRITE( stdout,'(5x,"Tr[ns] = ",f9.5)') 2.d0*nsuma(1)
                 ELSE
                    WRITE( stdout,'(5x,"Tr[ns] (up, down, total) = ",3f9.5)') &
                                  nsuma(1), nsuma(2), nsuma(1)+nsuma(2)
                    WRITE( stdout,'(5x,"Atomic magnetic moment = ",f9.5)') nsuma(1) - nsuma(2)
                 ENDIF
                 !
              ENDIF ! na1 = na2
              !
           ENDDO ! viz
           !
           ALLOCATE (f(ldim1,ldim1), vet(ldim1,ldim1), lambda(ldim1))
           !
           DO is = 1, nspin   
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
              CALL cdiagh(ldim1, f, ldim1, lambda, vet)
              !
              IF (nspin /= 1) WRITE( stdout,'(5x,"SPIN ",i2)') is
              WRITE( stdout,'(5x,"eigenvalues:")')
              WRITE( stdout,'(5x,7f7.3)') (lambda(m1), m1=1, ldim1)
              !
              WRITE( stdout,'(5x,"eigenvectors (columns):")')
              DO m1 = 1, ldim1
                 WRITE( stdout,'(5x,7f7.3)') ( DBLE(vet(m1,m2)), m2=1, ldim1 )
              ENDDO
              !
              WRITE( stdout,'(5x,"occupation matrix ns (before diag.):")')
              DO m1 = 1, ldim1
                 WRITE( stdout,'(5x,7f7.3)') ( DBLE(f(m1,m2)), m2=1, ldim1 )
              ENDDO
              !
           ENDDO ! is
           !
           DEALLOCATE (f, vet, lambda)
           !
        ENDIF
        !
     ENDDO !na1
     !
  ENDIF !iso_sys
  !
  IF (iso_sys) THEN
     !
     ! Compute the trace of the occupation matrix and the magnetization
     !     
     nsum = 0.d0
     rsrv = 0.d0
     !
     DO na1 = 1, nat
        !
        nt1 = ityp(na1)
        !
        IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
           !
           ldim1 = 2*Hubbard_l(nt1)+1
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              !
              IF (na2.EQ.na1) THEN
                 !
                 nsuma = 0.d0
                 DO is = 1, nspin
                    DO m1 = 1, ldim1
                       nsuma(is) = nsuma(is) + DBLE(nsgnew(m1,m1,viz,na1,is))
                       IF (reserv(nt1)) rsrv = rsrv + DBLE(nsgnew(m1,m1,viz,na1,is))
                    ENDDO
                    nsum = nsum + nsuma(is)
                 ENDDO
                 !
                 IF (nspin.EQ.1) THEN
                    WRITE( stdout,'(5x,"Tr[ns] = ",f9.5)') 2.d0*nsuma(1)
                 ELSE
                    WRITE( stdout,'(5x,"Tr[ns] (up, down, total) = ",3f9.5)') &
                                   nsuma(1), nsuma(2), nsuma(1)+nsuma(2)
                    WRITE( stdout,'(5x,"Atomic magnetic moment = ",f9.5)') nsuma(1) - nsuma(2)
                 ENDIF
                 !
              ENDIF ! na1 = na2
              !
           ENDDO ! viz
           !
        ENDIF
        !
     ENDDO !na1
     !
  ENDIF
  !
  ! THE CASE OF THE BACKGROUND STATES WITH V: RESERV_BACK
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard_back(nt1) ) THEN 
        !
        ! Compute the trace of the occupation matrix and the magnetization
        !
        DO viz = 1, neighood(na1)%num_neigh
           !
           na2 = neighood(na1)%neigh(viz)
           !
           IF ( na2.EQ.na1 ) THEN
              !
              nsuma = 0.d0
              DO is = 1, nspin
                 DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                    nsuma(is) = nsuma(is) + DBLE(nsgnew(m1,m1,viz,na1,is))
                    IF (reserv_back(nt1)) rsrv = rsrv + DBLE(nsgnew(m1,m1,viz,na1,is))
                 ENDDO
                 nsum = nsum + nsuma(is)
              ENDDO
              !
              IF (nspin.EQ.1) THEN
                 WRITE( stdout,'(5x,"Tr[ns] back = ",f9.5)') 2.d0*nsuma(1)
              ELSE
                 WRITE( stdout,'(5x,"Tr[ns] back (up, down, total) = ",3f9.5)') &
                                nsuma(1), nsuma(2), nsuma(1)+nsuma(2)
                 WRITE( stdout,'(5x,"Atomic magnetic moment (back) = ",f9.5)') nsuma(1) - nsuma(2)
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
  IF (nspin==1) nsum = 2.d0 * nsum
  !
  WRITE( stdout, '(/5x,a,1x,f9.4)') 'Number of occupied Hubbard levels =', nsum
  !
  IF (rsrv.GT.0.d0) &
     WRITE(stdout,'(5x,"Total occupation of reservoir states = ",x,f11.6)') rsrv
  !
  IF (iso_sys) DEALLOCATE (f, vet, lambda) 
  !
  RETURN
  !
END SUBROUTINE write_nsg 

!-----------------------------------------------------------------------
SUBROUTINE read_ns()
  !---------------------------------------------------------------------
  !
  ! This routine was written for the final SCF after vc-relax (M. Cococcioni).
  ! The occupations ns/nsg also need to be read in order to reproduce
  ! the right electronic ground state. When using Hubbard corrections
  ! this might be different (i.e. have different ordering of states) 
  ! from that simply obtained from the superposition of free ions.
  ! In other words the KS Hamiltonian (and its ground state) is also 
  ! functional of the Hubbard interaction parameters.
  !
  USE kinds,              ONLY : DP
  USE mp,                 ONLY : mp_bcast
  USE mp_images,          ONLY : intra_image_comm
  USE io_global,          ONLY : ionode, ionode_id
  USE scf,                ONLY : rho, v
  USE ldaU,               ONLY : lda_plus_u_kind, nsg, v_nsg, hub_back
  USE noncollin_module,   ONLY : noncolin
  USE io_files,           ONLY : restart_dir
  !
  IMPLICIT NONE
  INTEGER :: iunocc, iunocc1, ierr
  CHARACTER (LEN=256) :: dirname
  REAL(DP) :: eth, eth1
  !
  dirname = restart_dir()
  !
  IF ( ionode ) THEN
     !
     OPEN ( NEWUNIT=iunocc, FILE = TRIM(dirname) // 'occup.txt', &
            FORM='formatted', STATUS='old', IOSTAT=ierr )
     IF (lda_plus_u_kind.EQ.0) THEN
        READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
        IF (hub_back) THEN
           READ( UNIT = iunocc, FMT = * , iostat = ierr) rho%nsb
        ENDIF
     ELSEIF (lda_plus_u_kind.EQ.1) THEN
        IF (noncolin) THEN
           READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns_nc
        ELSE
           READ( UNIT = iunocc, FMT = *, iostat = ierr ) rho%ns
        ENDIF
     ELSEIF (lda_plus_u_kind.EQ.2) THEN
        READ( UNIT = iunocc, FMT = * , iostat = ierr) nsg
     ENDIF
     CLOSE(UNIT=iunocc,STATUS='keep')
     !
  ELSE
     !
     IF (lda_plus_u_kind.EQ.0) THEN
        rho%ns(:,:,:,:) = 0.D0
        IF (hub_back) rho%nsb(:,:,:,:) = 0.D0
     ELSEIF (lda_plus_u_kind.EQ.1) THEN
        IF (noncolin) THEN
           rho%ns_nc(:,:,:,:) = 0.D0
        ELSE
           rho%ns(:,:,:,:) = 0.D0
        ENDIF
     ELSEIF (lda_plus_u_kind.EQ.2) THEN
        nsg(:,:,:,:,:) = (0.d0, 0.d0)
     ENDIF
     !
  ENDIF
  !
  CALL mp_bcast( ierr, ionode_id, intra_image_comm )
  !
  IF (lda_plus_u_kind.EQ.0) THEN
     CALL mp_bcast(rho%ns, ionode_id, intra_image_comm)
     CALL v_hubbard (rho%ns, v%ns, eth)
     IF (hub_back) THEN
        CALL mp_bcast(rho%nsb, ionode_id, intra_image_comm)
        CALL v_hubbard_b (rho%nsb, v%nsb, eth1)
        eth = eth + eth1
     ENDIF
  ELSEIF (lda_plus_u_kind.EQ.1) THEN
     IF (noncolin) THEN
        CALL mp_bcast(rho%ns_nc, ionode_id, intra_image_comm)
        CALL v_hubbard_full_nc (rho%ns_nc, v%ns_nc, eth)
     ELSE
        CALL mp_bcast(rho%ns, ionode_id, intra_image_comm)
        CALL v_hubbard_full (rho%ns, v%ns, eth)
     ENDIF
  ELSEIF (lda_plus_u_kind.EQ.2) THEN
     CALL mp_bcast(nsg, ionode_id, intra_image_comm)
     CALL v_hubbard_extended (nsg, v_nsg, eth)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE read_ns
