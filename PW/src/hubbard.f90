!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE hub_summary()
    !--------------------------------------------------------------------------
    !!
    !! This routine writes a summary about Hubbard parameters
    !!
    USE ldaU,               ONLY : lda_plus_u_kind, Hubbard_projectors, is_hubbard, &
                                   Hubbard_U, Hubbard_J0, Hubbard_U2, Hubbard_alpha, &
                                   Hubbard_beta, Hubbard_alpha_back, is_hubbard_back,  &
                                   Hubbard_J, Hubbard_l, Hubbard_n, Hubbard_Um, &
                                   lda_plus_u, Hubbard_alpha_m, hub_pot_fix, apply_U, &
                                   orbital_resolved, Hubbard_Um_nc, Hubbard_alpha_m_nc
    USE lsda_mod,           ONLY : nspin 
    USE noncollin_module,   ONLY : lspinorb
    USE ions_base,          ONLY : ntyp => nsp
    USE io_global,          ONLY : stdout
    USE constants,          ONLY : rytoev, eps16
    USE ions_base,          ONLY : atm
    USE upf_utils,          ONLY : l_to_spdf
    !
    IMPLICIT NONE
    INTEGER :: nt, is, ldim, m1
    !
    WRITE(stdout,'(5x,a)') 'Hubbard projectors: ' // TRIM(Hubbard_projectors)
    IF (lda_plus_u_kind == 0 .AND. (.NOT. orbital_resolved)) THEN
       !Printing of orbital-resolved DFT+U is defined at the end of the routine
       WRITE( stdout, '(5x,"Hubbard parameters of DFT+U (Dudarev formulation) in eV:")')
       DO nt = 1, ntyp
          IF (is_hubbard(nt)) THEN
            CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
            CALL write_hub_param (nt, Hubbard_J0(nt), 'J0', 1)
            CALL write_hub_param (nt, Hubbard_alpha(nt), 'alpha', 1)
            CALL write_hub_param (nt, Hubbard_beta(nt), 'beta', 1)
          ENDIF
          IF (is_hubbard_back(nt)) THEN
            CALL write_hub_param (nt, Hubbard_U2(nt), 'U', 2)
            CALL write_hub_param (nt, Hubbard_alpha_back(nt), 'alpha', 2)
        ENDIF
       ENDDO
    ELSEIF (lda_plus_u_kind == 1) THEN
       WRITE( stdout, '(5x,"Hubbard parameters of DFT+U (Liechtenstein formulation) in eV:")')
       DO nt = 1, ntyp
          IF (Hubbard_U(nt) /= 0.d0) THEN
             IF (Hubbard_l(nt) == 0) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
             ELSEIF (Hubbard_l(nt) == 1) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
             ELSEIF (Hubbard_l(nt) == 2) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
                CALL write_hub_param (nt, Hubbard_J(2,nt), 'B', 1)
             ELSEIF (Hubbard_l(nt) == 3) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
                CALL write_hub_param (nt, Hubbard_J(2,nt), 'E2', 1)
                CALL write_hub_param (nt, Hubbard_J(3,nt), 'E3', 1)
             ENDIF
          ENDIF
       ENDDO
       IF (lspinorb) THEN
           WRITE(stdout, '(5x,"DFT+U on averaged j=l+1/2,l-1/2 radial WFs")')
       ENDIF
    ELSEIF (lda_plus_u_kind == 2) THEN
       ! Info about the Hubbard V is printed by the routine alloc_neighborhood
       IF (ANY(Hubbard_J0(:)/=0.d0) .OR. ANY(Hubbard_alpha(:)/=0.d0) .OR. &
           ANY(Hubbard_beta(:)/=0.d0) .OR. ANY(Hubbard_alpha_back(:)/=0.d0)) &
           WRITE( stdout, '(5x,"Hubbard parameters of DFT+U+V (Dudarev formulation) in eV:")')
       DO nt = 1, ntyp
          IF (is_hubbard(nt)) THEN
            CALL write_hub_param (nt, Hubbard_J0(nt), 'J0', 1)
            CALL write_hub_param (nt, Hubbard_alpha(nt), 'alpha', 1)
            CALL write_hub_param (nt, Hubbard_beta(nt), 'beta', 1)
          ENDIF
          IF (is_hubbard_back(nt)) THEN
            CALL write_hub_param (nt, Hubbard_alpha_back(nt), 'alpha', 2)
        ENDIF
       ENDDO
    ENDIF
    !
    IF (orbital_resolved) THEN
    ! ... currently not using "write_hub_param" because the latter expects 
    ! ... hub_parameter to be a scalar, but orbital-resolved Hubbard parameters
    ! ... are arrays and should be printed as such
    !
      WRITE( stdout, '(5x,"Orbital-resolved Hubbard parameters in eV:")')
      DO nt = 1, ntyp
         IF (is_hubbard(nt)) THEN
               IF (ANY(ABS(Hubbard_Um(:,:,nt)) .GT. eps16)) THEN
                  WRITE(stdout,'(5x,a,i1,a)')  &
                     'U' // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
                     l_to_spdf(Hubbard_l(nt),.FALSE.) // ')'
                  !
                  ldim = 2*Hubbard_l(nt)+1
                  DO is = 1, nspin
                     WRITE( stdout,'(5x,"spin-channel ",i2,": ", 8f8.4)') &
                              is,(Hubbard_Um(m1,is,nt)*rytoev, m1=1, ldim)
                  ENDDO
               ENDIF
               !
               IF (ANY(ABS(Hubbard_alpha_m(:,:,nt)) .GT. eps16)) THEN
                  WRITE(stdout,'(5x,a,i1,a)')  &
                     'ALPHA' // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
                     l_to_spdf(Hubbard_l(nt),.FALSE.) // ')'
                  !
                  DO is = 1, nspin
                     WRITE( stdout,'(5x,"spin-channel ",i2,": ", 8f8.4)') &
                              is,(Hubbard_alpha_m(m1,is,nt)*rytoev, m1=1, ldim)
                  ENDDO
               ENDIF
               !
               IF (ANY(ABS(Hubbard_Um_nc(:,nt)) .GT. eps16)) THEN
                  WRITE(stdout,'(5x,a,i1,a)')  &
                     'U' // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
                     l_to_spdf(Hubbard_l(nt),.FALSE.) // ')'
                  !
                  ldim = 2*Hubbard_l(nt)+1
                  WRITE( stdout,'(5x,15f8.4)') &
                           (Hubbard_Um_nc(m1,nt)*rytoev, m1=1, 2*ldim)
               ENDIF
               !
               IF (ANY(ABS(Hubbard_alpha_m_nc(:,nt)) .GT. eps16)) THEN
                  WRITE(stdout,'(5x,a,i1,a)')  &
                     'ALPHA' // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
                     l_to_spdf(Hubbard_l(nt),.FALSE.) // ')'
                  !
                  ldim = 2*Hubbard_l(nt)+1
                  WRITE( stdout,'(5x,15f8.4)') &
                           (Hubbard_alpha_m_nc(m1,nt)*rytoev, m1=1, 2*ldim)
               ENDIF
         ENDIF
      ENDDO
    ENDIF
    !
    WRITE(stdout, '(/5x,"Internal variables: lda_plus_u =",l, ", lda_plus_u_kind = ",i1)') &
            lda_plus_u, lda_plus_u_kind
    WRITE(stdout, '(/5x,"Internal variables: orbital_resolved =",l, ", hub_pot_fix =",l, ", apply_U = ",l)') &
            orbital_resolved, hub_pot_fix, apply_U
    !
    RETURN
    !
END SUBROUTINE hub_summary

!------------------------------------------------------------------------------
SUBROUTINE write_hub_param (nt, hub_parameter, hub_name, flag)
    !--------------------------------------------------------------------------
    !
    USE kinds,          ONLY : DP
    USE ions_base,      ONLY : atm
    USE constants,      ONLY : rytoev
    USE io_global,      ONLY : stdout
    USE ldaU,           ONLY : Hubbard_n, Hubbard_l, Hubbard_n2, Hubbard_l2, &
                               Hubbard_n3, Hubbard_l3, backall
    USE upf_utils,      ONLY : l_to_spdf
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nt
    CHARACTER(len=*), INTENT(IN) :: hub_name
    INTEGER, INTENT(IN) :: flag ! 1: first Hubbard channel
                                ! 2: second Hubbard channel
    REAL(DP) :: hub_parameter
    !
    IF (hub_parameter /= 0.d0) THEN
       IF (flag==1) THEN
          WRITE(stdout,'(5x,a,i1,a,f8.4)')  &
              hub_name // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
              l_to_spdf(Hubbard_l(nt),.FALSE.) // ') =', hub_parameter*rytoev
       ELSEIF (flag==2) THEN
          IF (.NOT.backall(nt)) THEN
             ! In this case there is one Hubbard channel for background states.
             WRITE(stdout,'(5x,a,i1,a,f8.4)')  &
                hub_name // '(' // TRIM(atm(nt)) // '-', Hubbard_n2(nt), &
                l_to_spdf(Hubbard_l2(nt),.FALSE.) // ') =', hub_parameter*rytoev
          ELSE
             ! In this case there are two Hubbard channels for background states.
             WRITE(stdout,'(5x,a,i1,a,i1,a,f8.4)')  &
                hub_name // '(' // TRIM(atm(nt))  // '-', Hubbard_n2(nt), &
                l_to_spdf(Hubbard_l2(nt),.FALSE.) // '-', Hubbard_n3(nt), &
                l_to_spdf(Hubbard_l3(nt),.FALSE.) // ') =', hub_parameter*rytoev
          ENDIF
       ENDIF
    ENDIF
    !
    RETURN
    !
END SUBROUTINE write_hub_param

!------------------------------------------------------------------------------
SUBROUTINE determine_hubbard_occ ( nt, lflag )
    !-----------------------------------------------------------------------
    !!
    !! This routine sets the starting Hubbard occupations based on the data
    !! read from the pseudopotentials.
    !! lflag=1 - first Hubbard manifold (main)
    !! lflag=2 - second Hubbard manifold
    !! lflag=3 - third Hubbard manifold
    !!
    USE kinds,        ONLY : DP
    USE uspp_param,   ONLY : upf
    USE ldaU,         ONLY : is_hubbard, is_hubbard_back, Hubbard_n, Hubbard_l, &
                             Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3,    &
                             hubbard_occ
    USE io_global,    ONLY : stdout
    USE upf_utils,    ONLY : l_to_spdf, capital
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nt          ! atomic type
    INTEGER, INTENT(IN) :: lflag       ! Hubbard channel
    !
    CHARACTER(LEN=2), ALLOCATABLE :: label(:)
    CHARACTER(LEN=2) :: label_aux
    CHARACTER(LEN=2) :: label_hub
    INTEGER :: i, & ! runs over all pseudo-atomic orbitals for the atomic type nt
               ldim
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    LOGICAL :: first
    !
    IF ( upf(nt)%nwfc < 1 ) THEN
       CALL errore('determine_hubbard_occ', 'no atomic wavefunctions in &
               &pseudopotential file for species #' // upf(nt)%psd // new_line('a') // &
               &'use a pseudopotential file with atomic wavefunctions!', 1)
    ENDIF
    !
    IF (lflag==1) THEN
       label_hub = TRIM(int_to_char(Hubbard_n(nt))) // l_to_spdf(Hubbard_l(nt))
    ELSEIF (lflag==2) THEN
       label_hub = TRIM(int_to_char(Hubbard_n2(nt))) //l_to_spdf(Hubbard_l2(nt))
    ELSEIF (lflag==3) THEN
       label_hub = TRIM(int_to_char(Hubbard_n3(nt))) //l_to_spdf(Hubbard_l3(nt))
    ELSE
       CALL errore('determine_hubbard_occ','Not allowed value of lflag',lflag)
    ENDIF
    !
    ldim = upf(nt)%nwfc
    ALLOCATE(label(ldim))
    !
    first=.true.
    DO i = 1, ldim
       ! Label of the i-th atomic orbital for the atomic type nt
       ! (if lowercase, then capitalize)
       label_aux = upf(nt)%els(i)
       label(i) = label_aux(1:1) // capital(label_aux(2:2))
       IF (label(i)==label_hub) THEN
          IF (first) THEN
             hubbard_occ(nt,lflag) = upf(nt)%oc(i)     
             first=.false.
          ELSE
             ! In FR-PPs the same Hubbard manifold appears twice,
             ! so we need to sum up all the electrons
             hubbard_occ(nt,lflag) = hubbard_occ(nt,lflag) &
                                    + upf(nt)%oc(i)
          ENDIF
       ENDIF
    ENDDO
    !
    IF (hubbard_occ(nt,lflag) < 0.d0) THEN
       WRITE(stdout,'(5x,a,8(1x,a))') TRIM(upf(nt)%psd) // &
               " pseudopotential contains the orbitals: ", (label(i),i=1,ldim)
       WRITE(stdout,'(5x,2a)') "Requested Hubbard manifold from the input: ", label_hub
       CALL errore('determine_hubbard_occ',&
               &'Mismatch between the requested and available manifolds',1)
    ENDIF
    !
    DEALLOCATE (label)
    !
    RETURN
    !
END SUBROUTINE determine_hubbard_occ
!
!----------------------------------------------------------------------------
SUBROUTINE order_eigenvecs(order, eigvec1, eigvec2, ldim)
!----------------------------------------------------------------------------
! This subroutine compares two arrays containing column-wise
! eigenvectors by computing the dot products between these vectors.
! Since the dot product of two almost identical vectors is close to 1
! while that of orthogonal vectors is 0, this information is used
! to keep track of the magnetic quantum orbitals
! during orbital-resolved DFT+U calculations.
! The idea was first proposed in JCTC 15,4871 (2019)
!
USE kinds,     ONLY : DP
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN)  :: ldim
COMPLEX(DP), INTENT(IN) :: eigvec1(ldim, ldim), eigvec2(ldim, ldim)
! column-wise eigenvectors to be compared
INTEGER, INTENT(OUT) :: order(ldim)
! vector describing the order of the first array of eigenvectors
! with respect to the second one. For example, if the 
! same eigenvector is stored in column 4 of eigvec1 and 
! in column 2 of eigvec2, then order(4)=2.
!
! local variables
INTEGER :: m1, m2, m3
REAL(DP) :: tmp, d
!
order = 0
!
IF (ALL(eigvec1(:,:) == eigvec2(:,:))) THEN
   DO m1 = 1, ldim
      order(m1) = m1
   ENDDO
ELSE
   DO m1 = 1, ldim
      tmp = 0.0
      DO m2 = 1, ldim
         d = 0.0
         DO m3 = 1, ldim
            ! compute the dot product
            d = d + REAL(eigvec1(m3, m1)) * REAL(eigvec2(m3, m2))
         ENDDO
         !
         IF (ABS(d) > tmp .AND. ALL(order .NE. m2)) THEN
            ! establish the order number: based on the rank of
            ! the dot product 
            tmp = ABS(d)
            order(m1) = m2
         ENDIF
#if defined(__DEBUG)
         WRITE(stdout, '("m1 :", i3, " ,dot product = ", f7.3, " ,order = ", i3)') m1, tmp, order(m1)
#endif
      ENDDO
      ! if no reference eigenvectors were stored yet, use default order
      IF ( tmp == 0.0 ) order(m1) = m1
   ENDDO
ENDIF
END SUBROUTINE order_eigenvecs
!
!----------------------------------------------------------------------------
SUBROUTINE diag_ns( ldim, ns, lambda, eigenvecs_current )
!---------------------------------------------------------------------
!
!! Diagonalizes the ns matrix and returns the
!! resulting eigenvalues & eigenvectors.
!
USE kinds,                ONLY : DP
USE lsda_mod,             ONLY : nspin
USE io_global,            ONLY : stdout
IMPLICIT NONE
!
INTEGER, INTENT(IN)        :: ldim
!! Number of magnetic quantum orbitals (2l+1)
REAL(DP), INTENT(IN)       :: ns(ldim,ldim,nspin)
!! Occupation matrix (undiagonalized)
REAL(DP), INTENT(OUT)      :: lambda(ldim,nspin)
!! Occupation eigenvalues
COMPLEX(DP), INTENT(OUT)   :: eigenvecs_current(ldim,ldim,nspin)
!! Eigenvectors of current iteration
!
!  ... local variables
!
REAL(DP)                   :: psum
COMPLEX(DP)                :: f(ldim,ldim)
INTEGER                    :: is,  m1, m2
!
DO is = 1, nspin
   !
   f(:,:) = CMPLX(0.d0,0.d0, kind=dp)
   eigenvecs_current(:,:,is) = CMPLX(0.d0,0.d0, kind=dp)
   lambda(:,is) = 0.d0
   !
   DO m1 = 1, ldim
      DO m2 = 1, ldim
         f(m1,m2) = ns(m1,m2,is)
      ENDDO
   ENDDO
   ! Make the occupation matrix hermitian
   DO m1 = 1, ldim
      DO m2 = m1, ldim
         psum = ABS( f(m1,m2) - f(m2,m1) )
         IF (psum > 1.d-10) THEN
            WRITE( stdout, * ) m1, m2, is
            WRITE( stdout, * ) f(m1,m2)
            WRITE( stdout, * ) f(m2,m1)
            CALL errore( 'diag_ns', 'non-hermitean occupation matrix', 1 )
         ELSE
            f(m1,m2) = 0.5d0 * (f(m1,m2) + f(m2,m1) )
            f(m2,m1) = f(m1,m2)
         ENDIF
      ENDDO
   ENDDO
   !
   ! Diagonalize the occupation matrix and save the eigenvectors in an array
   !
   CALL cdiagh( ldim, f, ldim, lambda(:,is), eigenvecs_current(:,:,is) )
   !
ENDDO
!
RETURN
END SUBROUTINE diag_ns
!
!----------------------------------------------------------------------------
SUBROUTINE diag_ns_nc( ldim, ns_nc, lambda, eigenvecs_current )
   !---------------------------------------------------------------------
   !
   !! Diagonalizes the noncollinear (complex) ns matrix and
   !! returns the resulting eigenvalues & eigenvectors.
   !
   USE kinds,                ONLY : DP
   USE lsda_mod,             ONLY : nspin
   USE io_global,            ONLY : stdout
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)        :: ldim
   !! Number of magnetic quantum orbitals (2l+1)
   COMPLEX(DP), INTENT(IN)    :: ns_nc(ldim,ldim,4)
   !! Occupation matrix (undiagonalized)
   REAL(DP), INTENT(OUT)      :: lambda(2*ldim)
   !! Occupation eigenvalues
   COMPLEX(DP), INTENT(OUT)   :: eigenvecs_current(2*ldim,2*ldim)
   !! Eigenvectors of current iteration
   !
   !  ... local variables
   !
   COMPLEX(DP)                :: f(2*ldim,2*ldim)
   INTEGER                    :: m1, m2
   !
   !
   f(:,:) = CMPLX(0.d0,0.d0, kind=dp)
   eigenvecs_current(:,:) = CMPLX(0.d0,0.d0, kind=dp)
   lambda(:) = 0.d0
   !
   DO m1 = 1, ldim
      DO m2 = 1, ldim
         f(m1, m2)           = ns_nc(m1,m2,1)
         f(m1, ldim+m2)      = ns_nc(m1,m2,2)
         f(ldim+m1, m2)      = ns_nc(m1,m2,3)
         f(ldim+m1, ldim+m2) = ns_nc(m1,m2,4)
      ENDDO
   ENDDO 
   !
   ! Diagonalize the occupation matrix and save the eigenvectors in an array
   !
   CALL cdiagh( 2*ldim, f, 2*ldim, lambda(:), eigenvecs_current(:,:) )
   !
   !
RETURN
END SUBROUTINE diag_ns_nc
!
!----------------------------------------------------------------------------
SUBROUTINE alpha_m_trace( ns )
!---------------------------------------------------------------------
!
!! Sums up and prints orbital occupation eigenvalues of states
!! perturbed by nonzero Hubbard_alpha_m. Also prints occupation eigenvalues
!! of other states having the same Hubbard U manifold as a perturbed species. 
!! Useful for the self-consistent computation of Hubbard parameters.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_n, &
                                 Hubbard_alpha_m, lambda_ns, &
                                 eigenvecs_ref, Hubbard_Um
USE lsda_mod,             ONLY : nspin
USE io_global,            ONLY : stdout
USE constants,            ONLY : rytoev, eps16
USE upf_utils,            ONLY : l_to_spdf
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN)     :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! Occupation matrix
!
!  ... local variables
!
COMPLEX(DP), ALLOCATABLE :: eigenvecs_current(:,:,:)
!
INTEGER, ALLOCATABLE     :: eigval_list(:), unpert_ats(:)
INTEGER                  :: i, is, na, nt, m, m1, m2, ldim, eigval_add
INTEGER                  :: order(2*Hubbard_lmax+1,nspin), m_order
!
REAL(DP)                 :: tralpha, u
!
LOGICAL                  :: has_second_manifold
!
CHARACTER(len=6)         :: manifold
!
!
ALLOCATE( unpert_ats(0) )
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
WRITE( stdout,'(/5x,"Occupations of perturbed manifolds (active HUBBARD ALPHA):")')
DO na = 1, nat
   !
   nt = ityp (na)
   !
   WRITE(manifold, "(i1,a1)") Hubbard_n(nt), l_to_spdf(Hubbard_l(nt),.FALSE.)
   ldim = 2*Hubbard_l(nt)+1
   tralpha = 0.0
   ALLOCATE( eigval_list(0), eigenvecs_current(ldim,ldim,nspin) )
   !
   IF ( ANY(Hubbard_alpha_m(:,:,nt) /= 0.d0) ) THEN
      ! ... case a) the species is directly affected by Hubbard_alpha_m:
      !             print the sum over the perturbed manifold's occupation eigenvalues
      !
      ! diagonalize the occupation matrix
      CALL diag_ns( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:ldim,1:ldim,na), eigenvecs_current )
      !
      DO is = 1, nspin
         !
         IF ( ANY(eigenvecs_ref(:,:,is,na) /= 0.0d0) ) THEN
            ! order the eigenstates
            order(:,is) = 0
            CALL order_eigenvecs( order(1:ldim,is), eigenvecs_current(1:ldim,1:ldim,is), &
                  eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         ELSE
            ! if this routine is called, there should always be reference 
            ! eigenvectors because if Hubbard_alpha /=0 the Hubbard 
            ! corrections are applied starting from the first iteration
            CALL errore( 'alpha_m_trace', 'missing reference eigenvectors', 1 )
            !
         ENDIF
         !
         DO m1 = 1, ldim    
            !
            IF (ABS(Hubbard_alpha_m(order(m1,is),is,nt)) >= eps16) THEN
               ! sum up the occupation of the perturbed states
               tralpha = tralpha + lambda_ns(order(m1,is),is,na)
               u = Hubbard_Um(order(m1,is),is,nt)*rytoev
               !
               ! Replicate the input sturcture of the HUBBARD card
               eigval_add = 0
               IF ( is .EQ.2 ) eigval_add = ldim
               !
               eigval_list = [ eigval_list,order(m1,is)+eigval_add ]
            ENDIF
            !
         ENDDO
      ENDDO
      !
      IF (nspin ==1) tralpha = tralpha*2
      !
      WRITE( stdout,'(/5x,"@ ATOM: ",i3," | MANIFOLD: ",a2," | U: ", f4.2, &
            & " | OCCUPATION: ", f10.8," | EIGVALS:", 21i3)') na,manifold,u,tralpha,(eigval_list(m),m=1,SIZE(eigval_list))

      !
   ELSEIF ( ANY(Hubbard_Um(:,:,nt) .NE. 0.d0) .AND. ALL(Hubbard_alpha_m(:,:,nt) .EQ. 0.d0) ) THEN
      ! ...case b) the species is NOT directly affected by Hubbard_alpha_m
      !            but (eventually, checked below) has a Hubbard_U manifold 
      !
      has_second_manifold = .TRUE.
      !
      IF ( .NOT. ANY(na == unpert_ats) ) THEN
         ! store the atom index of the unperturbed atom and the type
         ! of its perturbed counterpart in arrays (unless already there)
         unpert_ats = [ unpert_ats,na ]
      ENDIF
      !
   ENDIF
   !
   DEALLOCATE( eigval_list )
   DEALLOCATE( eigenvecs_current)
   !
ENDDO
!
IF (has_second_manifold) THEN
   ! ...case b) continuation: now print the occupations
   !            of the states with Hubbard_U manifolds
   !            equal to those of perturbed species
   WRITE( stdout,'(/5x,"Occupations of unperturbed manifolds (active HUBBARD U):")')
   !
   DO i = 1, SIZE(unpert_ats)
      !
      na = unpert_ats(i)
      nt = ityp (na)
      !
      WRITE(manifold, "(i1,a1)") Hubbard_n(nt), l_to_spdf(Hubbard_l(nt),.FALSE.)
      ldim = 2*Hubbard_l(nt)+1
      tralpha = 0.0
      ALLOCATE( eigval_list(0) )
      ALLOCATE( eigenvecs_current(ldim,ldim,nspin))
      !
      ! diagonalize the occupation matrix
      CALL diag_ns( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:ldim,1:ldim,na), eigenvecs_current )
      !
      DO is = 1, nspin
         !
         IF ( ANY(eigenvecs_ref(1:ldim,1:ldim,is,na) .NE. 0.0d0) ) THEN
            order(1:ldim,is) = 0
            CALL order_eigenvecs( order(1:ldim,is), eigenvecs_current(1:ldim,1:ldim,is), &
                  eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         ELSE
            CALL errore( 'alpha_m_trace', 'missing reference eigenvectors', 1 )
         ENDIF
         !
         DO m1 = 1, ldim
            !
            ! check if Um of the resepective eigenstate /= 0.0
            ! if yes, sum this state's occupation eigenvalue
            !
            IF ( Hubbard_Um(order(m1,is),is,nt) .NE. 0.d0 ) THEN
               tralpha = tralpha + lambda_ns(order(m1,is),is,na)
               u = Hubbard_Um(order(m1,is),is,nt)*rytoev
               !
               eigval_add = 0
               IF ( is .EQ.2 ) eigval_add = ldim
               !
               eigval_list = [ eigval_list,order(m1,is) + eigval_add ]
            ENDIF
            !
         ENDDO
      ENDDO
      !
      IF (nspin ==1) tralpha = tralpha*2
      WRITE( stdout,'(/5x,"@ ATOM: ",i3," | MANIFOLD: ",a2," | U: ", f4.2, &
      & " | OCCUPATION: ", f10.8," | EIGVALS:", 21i3)') na,manifold,u,tralpha,(eigval_list(m),m=1,SIZE(eigval_list))
      !
   DEALLOCATE( eigval_list )
   DEALLOCATE( eigenvecs_current)
   ENDDO
ENDIF
!
DEALLOCATE( unpert_ats )
!
RETURN
!
END SUBROUTINE alpha_m_trace
!-----------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE alpha_m_nc_trace( ns_nc )
!---------------------------------------------------------------------
!
!! Sums up and prints orbital occupation eigenvalues of states
!! perturbed by nonzero Hubbard_alpha_m. Also prints occupation eigenvalues
!! of other states having the same Hubbard U manifold as a perturbed species. 
!! Useful for the self-consistent computation of Hubbard parameters.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_n, &
                                 Hubbard_alpha_m_nc, lambda_ns, &
                                 eigenvecs_ref, Hubbard_Um_nc
USE lsda_mod,             ONLY : nspin
USE io_global,            ONLY : stdout
USE constants,            ONLY : rytoev, eps16
USE upf_utils,            ONLY : l_to_spdf
!
IMPLICIT NONE
!
COMPLEX(DP), INTENT(IN)     :: ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,4,nat)
!! Occupation matrix
!
!  ... local variables
!
COMPLEX(DP), ALLOCATABLE :: eigenvecs_current(:,:)
!
INTEGER, ALLOCATABLE     :: eigval_list(:), unpert_ats(:)
INTEGER                  :: i, is, na, nt, m, m1, m2, ldim
INTEGER                  :: order(4*Hubbard_lmax+2), m_order
!
REAL(DP)                 :: tralpha, u
!
LOGICAL                  :: has_second_manifold
!
CHARACTER(len=6)         :: manifold
!
!
ALLOCATE( unpert_ats(0) )
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
WRITE( stdout,'(/5x,"Occupations of perturbed manifolds (active HUBBARD ALPHA):")')
DO na = 1, nat
   !
   nt = ityp (na)
   !
   WRITE(manifold, "(i1,a1)")Hubbard_n(nt), l_to_spdf(Hubbard_l(nt),.FALSE.)
   ldim = 2*Hubbard_l(nt)+1
   tralpha = 0.0
   ALLOCATE( eigval_list(0), eigenvecs_current(2*ldim,2*ldim) )
   !
   IF ( ANY( Hubbard_alpha_m_nc(:,nt) /= 0.d0) ) THEN
      ! ... case a) the species is directly affected by Hubbard_alpha_m:
      !             print the sum over the perturbed manifold's occupation eigenvalues
      !
      ! diagonalize the occupation matrix
      CALL diag_ns_nc( ldim, ns_nc(1:ldim,1:ldim,:,na), lambda_ns(1:2*ldim,1,na), &
                  eigenvecs_current(1:2*ldim,1:2*ldim) )
      !
      IF ( ANY(eigenvecs_ref(:,:,1,na) .NE. 0.0d0) ) THEN
         ! order the eigenstates
         order(:) = 0
         CALL order_eigenvecs( order(:), eigenvecs_current(1:2*ldim,1:2*ldim), &
               eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
      ELSE
         ! if this routine is called, there should always be reference 
         ! eigenvectors because if Hubbard_alpha /=0 the Hubbard 
         ! corrections are applied starting from the first iteration
         CALL errore( 'alpha_m_nc_trace', 'missing reference eigenvectors', 1 )
         !
      ENDIF
      !
      DO m1 = 1, 2*ldim    
         !
         IF (ABS(Hubbard_alpha_m_nc(order(m1),nt)) .GE. eps16) THEN
            ! sum up the occupation of the perturbed states
            tralpha = tralpha + lambda_ns(order(m1),1,na)
            u = Hubbard_Um_nc(order(m1),nt)*rytoev
            !
            ! Replicate the input sturcture of the HUBBARD card
            eigval_list = [ eigval_list,order(m1) ]
         ENDIF
         !
      ENDDO
      !
      WRITE( stdout,'(/5x,"@ ATOM: ",i3," | MANIFOLD: ",a2," | U: ", f4.2, &
            & " | OCCUPATION: ", f10.8," | EIGVALS:", 21i3)') na,manifold,u,tralpha,(eigval_list(m),m=1,SIZE(eigval_list))
      !
   ELSEIF ( ANY(Hubbard_Um_nc(:,nt) .NE. 0.d0) .AND. &
            ALL(Hubbard_alpha_m_nc(:,nt) .EQ. 0.d0) ) THEN
      ! ...case b) the species is NOT directly affected by Hubbard_alpha_m
      !            but (eventually, checked below) has a Hubbard_U manifold 
      !
      has_second_manifold = .TRUE.
      !
      IF ( .NOT. ANY(na == unpert_ats) ) THEN
         ! store the atom index of the unperturbed atom and the type
         ! of its perturbed counterpart in arrays (unless already there)
         unpert_ats = [ unpert_ats,na ]
      ENDIF
      !
   ENDIF
   !
   DEALLOCATE( eigval_list )
   DEALLOCATE( eigenvecs_current)
   !
ENDDO
!
IF (has_second_manifold) THEN
   ! ...case b) continuation: now print the occupations
   !            of the states with Hubbard_U manifolds
   !            equal to those of perturbed species
   WRITE( stdout,'(/5x,"Occupations of unperturbed manifolds (active HUBBARD U):")')
   !
   DO i = 1, SIZE(unpert_ats)
      !
      na = unpert_ats(i)
      nt = ityp (na)
      !
      WRITE(manifold, "(i1,a1)")Hubbard_n(nt), l_to_spdf(Hubbard_l(nt),.FALSE.)
      ldim = 2*Hubbard_l(nt)+1
      tralpha = 0.0
      ALLOCATE( eigval_list(0) )
      ALLOCATE( eigenvecs_current(2*ldim,2*ldim) )
      !
      ! diagonalize the occupation matrix
      CALL diag_ns_nc( ldim, ns_nc(1:ldim,1:ldim,:,na), lambda_ns(1:2*ldim,1,na), &
                        eigenvecs_current(1:2*ldim,1:2*ldim) )
      !
      IF ( ANY(eigenvecs_ref(:,:,1,na) .NE. 0.0d0) ) THEN
         order(:) = 0
         CALL order_eigenvecs( order(:), eigenvecs_current(1:2*ldim,1:2*ldim), &
               eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
      ELSE
         CALL errore( 'alpha_m_trace', 'missing reference eigenvectors', 1 )
      ENDIF
      !
      DO m1 = 1, 2*ldim
         !
         ! check if Um of the resepective eigenstate /= 0.0
         ! if yes, sum this state's occupation eigenvalue
         !
         IF ( Hubbard_Um_nc(order(m1),nt) .NE. 0.d0 ) THEN
            tralpha = tralpha + lambda_ns(order(m1),1,na)
            u = Hubbard_Um_nc(order(m1),nt)*rytoev
            !
            eigval_list = [ eigval_list,order(m1) ]
         ENDIF
         !
      ENDDO
      !
      WRITE( stdout,'(/5x,"@ ATOM: ",i3," | MANIFOLD: ",a2," | U: ", f4.2, &
      & " | OCCUPATION: ", f10.8," | EIGVALS:", 21i3)') na,manifold,u,tralpha,(eigval_list(m),m=1,SIZE(eigval_list))
      !
   DEALLOCATE( eigval_list )
   DEALLOCATE( eigenvecs_current)
   ENDDO
ENDIF
!
DEALLOCATE( unpert_ats )
!
RETURN
   !
END SUBROUTINE alpha_m_nc_trace
   !-----------------------------------------------------------------------
