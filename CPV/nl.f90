!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 14 08:09:24 MET 1999
!  ----------------------------------------------

#include "f_defs.h"

      MODULE nl

        USE kinds
        USE spherical_harmonics
        USE cp_types
        USE nl_base

        IMPLICIT NONE
        SAVE

        PRIVATE

        INTERFACE nlsm1
          MODULE PROCEDURE nlsm1_s, nlsm1_v
        END INTERFACE
        INTERFACE nlsm2
          MODULE PROCEDURE nlsm2_s, nlsm2_v
        END INTERFACE

        PUBLIC :: nlsm1, nlrh_m

      CONTAINS


!  ----------------------------------------------
!  BEGIN manual

      REAL(dbl) FUNCTION nlrh_m(c0, cdesc, tforce, atoms, occ, kp, fnl, wsg, wnl, eigr)

!  this routine computes:
!  fnl: Kleinman-Bylander pseudopotential terms (see nlsm1,nlsm1_kp)
!  enl: nonlocal potential contribution to total energy (see ene_nl,ene_nl_kp)
!  nonlocal potential contribution to forces on ions
!
!    fion(n,ia) = -2 (sum over ib,igh,ik) kp%weight(ik) occ%s(ib,ik)
!                 Re { conjugate(dfnl(ia,igh,ib,ik)) fnl(ia,igh,ib,ik) }
!                 ps%wsg(igh,is)
!
!    kp%weight(ik) = weight of k point
!    occ%s(ib,ik) = occupation number
!    dfnl(ia,igh,ib,ik) = derivative of fnl with respect to R(n,ia)  n = x,y,z
!    fnl(ia,igh,ib,ik) = Kleinman-Bylander factor (see nlsm1,nlsm1_kp)
!    ps%wsg(igh,is) = inverse denominator in Kleinman-Bylander's formula (see nlset)
!  ----------------------------------------------
!  END manual

! ... include modules
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector, allocate_projector, deallocate_projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(dbl)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE (kpoints),        INTENT(IN)     :: kp            ! G and k vectors
      COMPLEX(dbl),           INTENT(INOUT)  :: c0(:,:,:,:)  ! wave functions
      TYPE (wave_descriptor), INTENT(IN)    :: cdesc         ! wave functions descriptor
      REAL(dbl),           INTENT(IN)     :: occ(:,:,:)      ! occupations
      LOGICAL,               INTENT(IN)     :: tforce        ! if .TRUE. compute forces on ions
      TYPE(atoms_type),      INTENT(INOUT)  :: atoms         ! ions structure
      REAL(dbl),             INTENT(IN)     :: wsg(:,:)      ! KB inverse denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector),      INTENT(OUT)    :: fnl(:,:)      ! KB  factors
      REAL(dbl),             INTENT(INOUT)   :: wnl(:,:,:,:) ! KB products <Y phi V | exp(i(k+G) dot r)>

! ... declare other variables
      INTEGER :: ispin, ispin_wfc

!  end of declarations
!  ----------------------------------------------

      CALL nl_projector_m(c0, cdesc, atoms, kp, fnl, wsg, wnl, eigr)
      nlrh_m = nl_energy_m(fnl, atoms, occ, kp, wsg)
      ! WRITE( stdout,*) 'DEBUG nlrh ', SUM( fnl(1,1)%r ), SUM( wsg ), SUM( wnl )
      IF( tforce ) THEN
        DO ispin = 1, cdesc%nspin
          ispin_wfc = ispin
          IF( force_pairing ) ispin_wfc = 1
          CALL nl_ionic_forces_v( ispin, c0( :, :, :, ispin_wfc ), cdesc, atoms, occ(:,:,ispin), &
            kp, fnl(:,ispin), wsg, wnl, eigr)
        END DO
      END IF
      !  .. WRITE( stdout,*) 'DEBUG NLRH:',atoms%for(:,1:atoms%nat)
      RETURN
      END FUNCTION nlrh_m


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE nl_projector_m(c0, cdesc, atoms, kp, fnl, wsg, wnl, eigr)

!  this routine computes:
!  fnl: Kleinman-Bylander pseudopotential terms (see nlsm1,nlsm1_kp)
!
!  ----------------------------------------------
!  END manual

! ... include modules
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing, gamma_only
      USE reciprocal_space_mesh, ONLY: gkx_l, gk_l
      USE reciprocal_vectors, ONLY: gx, g

      IMPLICIT NONE

! ...   declare subroutine arguments
      COMPLEX(dbl)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE(atoms_type),      INTENT(INOUT)  :: atoms         ! ions structure
      TYPE (kpoints), INTENT(IN)  :: kp       ! G and k vectors
      COMPLEX(dbl),    INTENT(INOUT)  :: c0(:,:,:,:)  ! wave functions
      TYPE (wave_descriptor),  INTENT(IN)  :: cdesc  ! wave functions
      REAL(dbl),      INTENT(IN)  :: wsg(:,:) ! Kleinman-Bylander inverse
                                              ! denominators
                                              ! <Y phi | V | phi Y>**(-1)
! ...   Kleinman-Bylander factors, with and without k points
      TYPE (projector), INTENT(OUT) :: fnl(:,:)
! ...   Kleinman-Bylander products <Y phi V | exp(i(k+G) dot r)>
      REAL(dbl), INTENT(IN) :: wnl(:,:,:,:)
! ...   declare other variables
      INTEGER      :: ik, ispin, ispin_wfc
!  end of declarations
!  ----------------------------------------------

      DO ispin = 1, cdesc%nspin
        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1
        DO ik = 1, cdesc%nkl
          IF( gamma_only ) THEN
            CALL nlsm1_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:, :, ik, ispin_wfc), cdesc, g, gx, fnl(ik, ispin))
          ELSE
            CALL nlsm1_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:, :, ik, ispin_wfc), cdesc, &
                          gk_l(:,ik), gkx_l(:,:,ik), fnl(ik, ispin))
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE nl_projector_m


!  ----------------------------------------------
!  BEGIN manual

      REAL(dbl) FUNCTION nl_energy_m(fnl, atoms, occ, kp, wsg)

!  this routine computes:
!  enl: nonlocal potential contribution to total energy (see ene_nl,ene_nl_kp)
!  ----------------------------------------------
!  END manual

! ... include modules
      USE brillouin, ONLY: kpoints
      USE pseudopotential, ONLY: nspnl
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (kpoints), INTENT(IN)    :: kp       ! G and k vectors
      REAL(dbl),    INTENT(IN)    :: occ(:,:,:) ! occupations
      TYPE(atoms_type),  INTENT(IN) :: atoms    ! ions structure
      REAL(dbl),      INTENT(IN)    :: wsg(:,:) ! Kleinman-Bylander inverse
                                                ! denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector) :: fnl(:,:) ! ... Kleinman-Bylander factors

! ... declare other variables
      REAL(dbl)    :: enl
      INTEGER      :: ik, ispin

!  end of declarations
!  ----------------------------------------------

! ... compute nonlocal contribution to total energy
      enl = 0.0d0
      DO ispin = 1, SIZE( fnl, 2 )
        DO ik = 1, SIZE( fnl, 1 )
          enl = enl + kp%weight(ik) * &
             ene_nl(fnl(ik,ispin), wsg, occ(:,ik,ispin), nspnl, atoms%na(:))
        END DO
      END DO
      nl_energy_m = enl
      RETURN
      END FUNCTION nl_energy_m



!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE nl_ionic_forces_m(c0, cdesc, atoms, occ, kp, fnl, wsg, wnl, eigr)

!  this routine computes:
!  fnl/fnlk: Kleinman-Bylander pseudopotential terms (see nlsm1,nlsm1_kp)
!  enl: nonlocal potential contribution to total energy (see ene_nl,ene_nl_kp)
!  nonlocal potential contribution to forces on ions
!
!    fion(n,ia) = -2 (sum over ib,igh,ik) kp%weight(ik) occ%s(ib,ik)
!                 Re { conjugate(dfnlk(ia,igh,ib,ik)) fnlk(ia,igh,ib,ik) }
!                 ps%wsg(igh,is)
!
!    kp%weight(ik) = weight of k point
!    occ%s(ib,ik) = occupation number
!    dfnlk(ia,igh,ib,ik) = derivative of fnlk with respect to R(n,ia)
!                          n = x,y,z
!    fnlk(ia,igh,ib,ik) = Kleinman-Bylander factor (see nlsm1,nlsm1_kp)
!    ps%wsg(igh,is) = inverse denominator in Kleinman-Bylander's formula
!                     (see nlset)
!  ----------------------------------------------
!  END manual

! ... include modules
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(dbl)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE (kpoints), INTENT(IN)  :: kp       ! G and k vectors
      COMPLEX(dbl),    INTENT(INOUT)  :: c0(:,:,:,:)    ! wave functions
      TYPE (wave_descriptor),    INTENT(IN)  :: cdesc    ! wave functions desc.
      REAL(dbl),    INTENT(IN)  :: occ(:,:,:)      ! occupations
      TYPE(atoms_type),  INTENT(INOUT) :: atoms    ! ions structure
      REAL(dbl),      INTENT(IN)  :: wsg(:,:)  ! KB inverse denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector), INTENT(INOUT) :: fnl(:,:) ! KB factors
      REAL(dbl), INTENT(INOUT) :: wnl(:,:,:,:)    ! KB products <Y phi V | exp(i(k+G) dot r)>
      INTEGER :: ispin, ispin_wfc
      !
      DO ispin = 1, cdesc%nspin
        !
        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1
        !
        CALL nl_ionic_forces_v( ispin, c0( :, :, :, ispin_wfc), cdesc, atoms, occ(:,:,ispin), &
          kp, fnl(:,ispin), wsg, wnl, eigr)
        !
      END DO
      !
      RETURN
      END SUBROUTINE nl_ionic_forces_m

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE nl_ionic_forces_v( ispin, c0, cdesc, atoms, occ, kp, fnl, wsg, wnl, eigr)

!  this routine computes:
!  nonlocal potential contribution to forces on ions
!
!    fion(n,ia) = -2 (sum over ib,igh,ik) kp%weight(ik) occ%s(ib,ik)
!                 Re { conjugate(dfnl(ia,igh,ib,ik)) fnl(ia,igh,ib,ik) } ps%wsg(igh,is)
!
!    kp%weight(ik) = weight of k point
!    occ%s(ib,ik) = occupation number
!    fnl(ia,igh,ib,ik) = Kleinman-Bylander factor 
!    dfnl(ia,igh,ib,ik) = derivative of fnl with respect to R(n,ia) n = x,y,z
!    ps%wsg(igh,is) = inverse denominator in Kleinman-Bylander's formula  (see nlset)
!  ----------------------------------------------
!  END manual

! ... include modules
      USE pseudopotential, ONLY: nspnl, nsanl, ngh
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector, allocate_projector, deallocate_projector
      USE atoms_type_module, ONLY: atoms_type
      USE reciprocal_space_mesh, ONLY: gkx_l, gk_l
      USE control_flags, ONLY: gamma_only
      USE reciprocal_vectors, ONLY: gx, g

      IMPLICIT NONE

! ... declare subroutine arguments
      COMPLEX(dbl)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE(atoms_type),       INTENT(INOUT) :: atoms ! ions structure
      TYPE (kpoints),         INTENT(IN)    :: kp        ! K points
      COMPLEX(dbl),            INTENT(INOUT) :: c0(:,:,:)     ! wave functions
      TYPE (wave_descriptor), INTENT(IN)    :: cdesc     ! wave functions desc.
      REAL(dbl),            INTENT(IN)    :: occ(:,:)       ! occupations
      REAL(dbl),              INTENT(IN)    :: wsg(:,:)  ! KB inverse denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector),       INTENT(INOUT) :: fnl(:)   ! KB factors
      REAL(dbl),              INTENT(INOUT) :: wnl(:,:,:,:)    ! KB products <Y phi V | exp(i(k+G) dot r)>
      INTEGER, INTENT(IN) :: ispin

! ... declare other variables
      TYPE (projector) :: dfnl
      COMPLEX(dbl) :: ctmp
      REAL(dbl)    :: temp, fac, tt, fsum, enl
      INTEGER      :: me, ib, ia, k, isa, igh, is, ik

!  end of declarations
!  ----------------------------------------------


! ... compute nonlocal contribution to forces on ions
      KAPPA: DO ik = 1, cdesc%nkl
        CARTE: DO k = 1, 3  ! x,y,z directions
          CALL allocate_projector( dfnl, nsanl, cdesc%nbl( ispin ), ngh, cdesc%gamma )
          IF( gamma_only ) THEN
            CALL nlsm2_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:,:,ik), cdesc, g, gx, dfnl, k)
          ELSE
            CALL nlsm2_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:,:,ik), cdesc, &
                          gk_l(:,ik), gkx_l(:,:,ik), dfnl, k)
          END IF
          CHANN: DO igh = 1, ngh
            BANDE: DO ib = 1, cdesc%nbl( ispin )
              isa=0
              SPECS: DO is = 1, nspnl
                temp = 2.d0 * wsg(igh,is) * occ(ib,ik) * kp%weight(ik)
                IF( fnl(ik)%gamma_symmetry ) THEN
                  DO ia = 1, atoms%na(is)
                    isa = isa + 1
                    tt = dfnl%r(isa,igh,ib) * fnl(ik)%r(isa,igh,ib)
                    atoms%for(k,isa) = atoms%for(k,isa) - tt * temp
                  END DO
                ELSE
                  DO ia = 1, atoms%na(is)
                    isa = isa+1
                    tt = REAL( CONJG( dfnl%c(isa,igh,ib) ) * fnl(ik)%c(isa,igh,ib) )
                    atoms%for(k,isa) = atoms%for(k,isa) - tt * temp
                  END DO
                END IF
              END DO SPECS
            END DO BANDE
          END DO CHANN
          CALL deallocate_projector( dfnl )
        END DO CARTE
      END DO KAPPA
      RETURN
      END SUBROUTINE nl_ionic_forces_v


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE nlsm1_s( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, fnl)

!  this routine computes the Kleinman-Bylander factors of the nonlocal
!  part of the pseudopotentials
!
!    fnl(ia,ib,igh) = (sum over ig) c0%w(ib,ig) exp(i G dot R(ia))
!                     < vnl(is,igh) phi(igh) Y(igh) | exp(i G dot r) >
!
!    wnl(is,igh) = nonlocal part of pseudopotentials
!    phi(igh) = reference (isolated-atom) state for pseudopotentials
!    Y(igh) = s, p_x, p_y, p_z ... orbitals
!    R(ia) = positions of ions
!
!    ia = index of ion
!    ib = index of band
!    ig = index of G vector
!    igh = index of orbital
!  ----------------------------------------------
!  END manual

! ... declare modules

      USE pseudopotential, ONLY: l2ind, nspnl, tl, lm1x
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE mp
      USE mp_global, ONLY: nproc, mpime, group
      USE atoms_type_module, ONLY: atoms_type

      IMPLICIT NONE 

! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: ispin
      COMPLEX(dbl),      INTENT(INOUT) :: c(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL(dbl)                       :: wnl(:,:,:)
      REAL(dbl)                       :: g2(:), gx(:,:)
      TYPE (projector), INTENT(OUT)   :: fnl
      COMPLEX(dbl) :: eigr(:,:)

! ... declare other variables
      INTEGER :: is, igh, isa, iss, ia, ig, nb
      INTEGER :: l, ll, m, ngw, lda, ldw, ldf
      REAL(dbl), ALLOCATABLE :: gwork(:)
      COMPLEX(dbl), ALLOCATABLE :: gxtmp(:)
      COMPLEX(dbl), ALLOCATABLE :: auxc(:,:)
      COMPLEX(dbl), PARAMETER :: ONE  = (1.0d0,0.0d0)
      COMPLEX(dbl), PARAMETER :: ZERO = (0.0d0,0.0d0)
! ... i^l
      COMPLEX(dbl), PARAMETER :: csign(0:3) = (/ (1.0d0, 0.0d0), &
          (0.0d0,1.0d0), (-1.0d0,0.0d0), (0.0d0,0.0d0) /)

!  end of declarations
!  ----------------------------------------------

      ngw = cdesc%ngwl
      nb  = cdesc%nbl( ispin )
      IF( cdesc%gamma ) THEN
        lda = 2*ngw
        ldw = 2*SIZE(c, 1)
        ldf = SIZE(fnl%r, 1) * SIZE(fnl%r, 2)
        fnl%r = 0.0d0
      ELSE
        lda = ngw
        ldw = SIZE(c, 1)
        ldf = SIZE(fnl%c, 1) * SIZE(fnl%c, 2)
        fnl%c = 0.0d0
      END IF

      ! WRITE( stdout,fmt="('DEBUG nlsm1 ',10I6)" ) ngw, nb, lda, ldw, ldf

      ALLOCATE( gwork(ngw), gxtmp(ngw) )

! ... angular momentum l = 0
! ... orbital: s
! ... angular momentum l = 1
! ... orbitals: p_x, p_y, p_z
! ... angular momentum l = 2
! ... orbitals: d_{z^2}, d_{x^2-y^2}, d_{xy}, d_{yz}, d_{zx}

      igh = 0
      DO l = 0, lm1x
        IF(tl(l)) THEN
          DO m = -l, l
            igh = igh + 1
            CALL spharm( gwork, gx, g2, ngw, l, m )
            iss = 1
            DO is = 1, nspnl
              ll  = l2ind( l + 1, is )
              IF(ll.gt.0) THEN
                ALLOCATE( auxc( ngw, atoms%na(is) ) )
                gxtmp(1:ngw) = csign(l) * wnl(1:ngw,ll,is) * gwork(1:ngw)
                IF( cdesc%gamma .AND. cdesc%gzero ) gxtmp(1) = gxtmp(1) * 0.5d0
                DO ia = 1, atoms%na(is)
                  auxc(1:ngw,ia) = gxtmp(1:ngw) * eigr(1:ngw,iss+ia-1)
                END DO
                ! DO ia = 1, atoms%na(is)
                ! DO ig = 1, ngw
                ! WRITE( stdout7,fmt="(3I5,2D18.8)") is, ia, ig, auxc(ig,ia)
                ! END DO
                ! END DO

                ! WRITE( stdout,*) 'DEBUG nlsm_s', ll, ' ', m, ' ', is, SUM(auxc), atoms%na(is), SUM(c)
                IF ( cdesc%gamma ) THEN
                  ! CALL DGEMUL( auxc(1,1), lda, 'T', c(1,1), ldw, 'N', &
                  !   fnl%r(iss,igh,1), ldf, atoms%na(is), 2*ngw, nb )
                  CALL DGEMM( 'T', 'N', atoms%na(is), nb, 2*ngw, 1.0d0, &
                    auxc(1,1), lda, c(1,1), ldw, 0.0d0, fnl%r(iss,igh,1), ldf)
                ELSE
                  CALL ZGEMM( 'C', 'N', atoms%na(is), nb, ngw, one, auxc(1,1), lda, &
                    c(1,1), ldw, zero, fnl%c(iss,igh,1), ldf )
                END IF
                DEALLOCATE(auxc)
                ! WRITE( stdout,*) 'DEBUG nlsm_s', tl(l), l, SUM( fnl%r )
              END IF
              iss = iss + atoms%na(is)
            END DO
          END DO
        END IF
      END DO
      DEALLOCATE(gwork, gxtmp)

! ... since G vectors only span half space, multiply results by two
      IF ( cdesc%gamma ) THEN
        CALL DSCAL( size( fnl%r ), 2.0d0, fnl%r(1,1,1), 1 )
      !  WRITE( stdout,*) ' DEBUG nlsm1: ', SIZE(fnl%r,1), SIZE(fnl%r,2), SIZE(fnl%r,3) ! DEBUG
        CALL mp_sum( fnl%r, group )
      ELSE
        CALL mp_sum( fnl%c, group )
      END IF
      
      RETURN
      END SUBROUTINE nlsm1_s

      SUBROUTINE nlsm1_v( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, fnl)
!! ...   declare modules
        USE wave_types, ONLY: wave_descriptor
        USE pseudo_projector, ONLY: projector
        USE atoms_type_module, ONLY: atoms_type
        IMPLICIT NONE
!! ...   declare subroutine arguments
        COMPLEX(dbl),    INTENT(INOUT) :: c(:,:,:)
        TYPE (wave_descriptor),  INTENT(IN) :: cdesc
        INTEGER, INTENT(IN) :: ispin
        TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
        REAL(dbl) :: wnl(:,:,:,:)
        REAL(dbl) :: g2(:,:), gx(:,:,:)
        TYPE (projector), INTENT(OUT) :: fnl(:)
        COMPLEX(dbl) :: eigr(:,:)
        INTEGER :: i
          DO i = 1, SIZE(fnl)
            CALL nlsm1_s( ispin, wnl(:,:,:,i), atoms, eigr, c(:,:,i), cdesc, g2(:,i), &
              gx(:,:,i), fnl(i) )
          END DO
        RETURN
      END SUBROUTINE

      SUBROUTINE nlsm2_v( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, dfnl, kk)
!! ...   declare modules
        USE wave_types, ONLY: wave_descriptor
        USE pseudo_projector, ONLY: projector
        USE atoms_type_module, ONLY: atoms_type
        IMPLICIT NONE
!! ...   declare subroutine arguments
        INTEGER, INTENT(IN) :: ispin
        COMPLEX(dbl),    INTENT(IN) :: c(:,:,:)
        TYPE (wave_descriptor),  INTENT(IN) :: cdesc
        TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
        REAL(dbl) :: wnl(:,:,:,:)
        REAL(dbl) :: g2(:,:), gx(:,:,:)
        TYPE (projector), INTENT(OUT) :: dfnl(:)
        COMPLEX(dbl) :: eigr(:,:)
        INTEGER, INTENT(IN) :: kk
        INTEGER :: i
          DO i = 1, SIZE(dfnl)
            CALL nlsm2_s( ispin, wnl(:,:,:,i), atoms, eigr, c(:,:,i), cdesc, g2(:,i), &
              gx(:,:,i), dfnl(i), kk)
          END DO
        RETURN
      END SUBROUTINE


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE nlsm2_s( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, dfnl, kk)

!  this routine computes the derivatives of the Kleinman-Bylander
!  factors fnl, to be used for Hellmann-Feynman forces evaluation
!
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE pseudopotential, ONLY: l2ind, lm1x, nspnl, tl
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE cell_base, ONLY: tpiba

      IMPLICIT   NONE

! ... declare subroutine arguments
      COMPLEX(dbl), INTENT(IN) :: c(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (projector), INTENT(OUT) :: dfnl
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL(dbl), INTENT(IN) :: wnl(:,:,:), g2(:), gx(:,:)
      COMPLEX(dbl), INTENT(IN) :: eigr(:,:)
      INTEGER, INTENT(IN) :: kk
      INTEGER, INTENT(IN) :: ispin

! ... declare other variables
      REAL(dbl), ALLOCATABLE :: gwork(:)
      INTEGER :: is, ia, igh, isa, ig, iss, ll, l, m, ngw, nb, lda, ldw, ldf
      COMPLEX(dbl), ALLOCATABLE :: auxc(:,:), gxtmp(:)
      COMPLEX(dbl), PARAMETER :: ONE  = (1.0d0,0.0d0), ZERO = (0.0d0,0.0d0)
! ... (-i) * i^l
      COMPLEX(dbl), PARAMETER :: csign(0:3) = (/ (0.0d0,-1.0d0), &
        (1.0d0,0.0d0), (0.0d0,1.0d0), (0.0d0,0.0d0) /)

!  end of declarations
!  ----------------------------------------------

      ngw = cdesc%ngwl
      nb  = cdesc%nbl( ispin )
      IF( cdesc%gamma ) THEN
        lda = 2*SIZE(c, 1)
        ldw = 2*SIZE(c, 1)
        ldf = SIZE(dfnl%r, 1) * SIZE(dfnl%r, 2)
        dfnl%r = 0.0d0
      ELSE
        lda = SIZE(c, 1)
        ldw = SIZE(c, 1)
        ldf = SIZE(dfnl%c, 1) * SIZE(dfnl%c, 2)
        dfnl%c = 0.0d0
      END IF

      ALLOCATE(gwork(ngw), gxtmp(ngw))
      igh = 0
      ANGMO: DO l = 0, lm1x
        IF(tl(l)) THEN
          MAGNE: DO m = -l, l
            igh = igh + 1
            iss = 1
            CALL spharm(gwork, gx, g2, NGW, L, M)
            gwork(1:ngw) = tpiba * gx(kk,1:ngw) * gwork(1:ngw)
            SPECS: DO is = 1, nspnl
              ll  = l2ind(l + 1,is)
              IF(ll.gt.0) THEN
                ALLOCATE(auxc(ngw,atoms%na(is)))
                gxtmp(1:ngw) = csign(l) * wnl(1:ngw,ll,is) * gwork(1:ngw)
                DO ia = 1, atoms%na(is)
                  auxc(1:ngw,ia) = gxtmp(1:ngw) * eigr(1:ngw,iss + ia - 1)
                END DO
                IF( cdesc%gamma ) THEN
                  ! CALL DGEMUL(auxc(1,1), lda, 'T', c%w(1,1), ldw, 'N', &
                  !   dfnl%r(iss,igh,1), ldf, atoms%na(is), 2*ngw, nb)
                  CALL DGEMM('T', 'N', atoms%na(is), nb, 2*ngw, 1.0d0, auxc(1,1), lda, &
                    c(1,1), ldw, 0.0d0, dfnl%r(iss,igh,1), ldf)
                ELSE
                  CALL ZGEMM('C', 'N', atoms%na(is), nb, ngw, one, auxc(1,1), lda, &
                    c(1,1), ldw, zero, dfnl%c(iss,igh,1), ldf)
                END IF
                DEALLOCATE(auxc)
              END IF
              iss = iss + atoms%na(is)
            END DO SPECS
          END DO MAGNE
        END IF
      END DO ANGMO
      IF( cdesc%gamma ) CALL DSCAL(size(dfnl%r),2.0d0,dfnl%r(1,1,1),1)
      DEALLOCATE(gwork, gxtmp)
      RETURN
      END SUBROUTINE nlsm2_s


      END MODULE nl
