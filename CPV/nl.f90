!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!----------------------------------------------
 MODULE nl
!----------------------------------------------

        USE kinds
        USE spherical_harmonics
        USE nl_base

        IMPLICIT NONE
        SAVE

        PRIVATE


        PUBLIC :: nlrh_m, nlsm1_s

!----------------------------------------------
 CONTAINS
!----------------------------------------------



   REAL(DP) FUNCTION nlrh_m( c0, cdesc, tforce, atoms, occ, bec, becdr, eigr )

      !  this routine computes:
      !  Kleinman-Bylander pseudopotential terms (see nlsm1,nlsm1_kp)
      !  enl: nonlocal potential contribution to total energy (see ene_nl,ene_nl_kp)
      !  nonlocal potential contribution to forces on ions, see nlsm2
      !
      ! ... include modules

      USE brillouin,         ONLY: kpoints, kp
      USE wave_types,        ONLY: wave_descriptor
      USE pseudo_projector,  ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags,     ONLY: force_pairing
      USE pseudopotential,   ONLY: nspnl
      USE electrons_base,    ONLY: iupdwn, nupdwn
      USE uspp,              ONLY: becsum, nkb

      IMPLICIT NONE

      ! ... declare subroutine arguments

      COMPLEX(DP)                          :: eigr(:,:)     ! exp(i G dot r)
      COMPLEX(DP),           INTENT(INOUT) :: c0(:,:,:,:)   ! wave functions
      TYPE (wave_descriptor), INTENT(IN)    :: cdesc         ! wave functions descriptor
      REAL(DP),           INTENT(IN)       :: occ(:,:,:)    ! occupations
      LOGICAL,               INTENT(IN)     :: tforce        ! if .TRUE. compute forces on ions
      TYPE(atoms_type),      INTENT(INOUT)  :: atoms         ! ions structure
      REAL(DP)                             :: bec(:,:)
      REAL(DP)                             :: becdr(:,:,:)

      REAL(DP)    :: ennl
      EXTERNAL     :: ennl

      ! ... declare other variables
      !
      INTEGER      :: iss, iss_wfc, i, j
      REAL(DP)    :: etmp
      REAL(DP), ALLOCATABLE :: btmp( :, :, : )
      REAL(DP), ALLOCATABLE :: fion( :, : )

!  end of declarations
!  ----------------------------------------------

      DO iss = 1, cdesc%nspin
         !
         iss_wfc = iss
         IF( force_pairing ) iss_wfc = 1
         !
         CALL nlsm1 ( cdesc%nbl( iss ), 1, nspnl, eigr(1,1),    &
                      c0( 1, 1, 1, iss_wfc ), bec(1, iupdwn( iss ) ) )
         !
         IF( tforce ) THEN
            !
            ALLOCATE( btmp( nkb, nupdwn( iss ), 3 ) ) 
            !
            CALL nlsm2( cdesc%ngwl, nkb, nupdwn( iss ), eigr(1,1), &
                        c0( 1, 1, 1, iss_wfc ), btmp( 1, 1, 1 ), .false. )
            !
            DO i = 1, 3
               DO j = iupdwn( iss ), iupdwn( iss ) + nupdwn( iss ) - 1
                  becdr( :, j , i ) = btmp( :, j - iupdwn( iss ) + 1, i ) 
               END DO
            END DO
            !
            DEALLOCATE( btmp )
            !
         END IF
         !
      END DO
      
      nlrh_m = ennl( becsum, bec )

      IF( tforce ) THEN
         !
         CALL force_nl( atoms%for, bec, becdr )
         !
      END IF

      ! CALL nl_projector_m(c0, cdesc, atoms, fnl, wsg, wnl, eigr)
      !
      ! nlrh_m = nl_energy_m(fnl, atoms, occ, wsg)
      !
      ! WRITE(6,*) ' DEBUG nlrh (enl) = ', etmp - nlrh_m, tforce
      !
      ! IF( tforce ) THEN
        !
        ! ALLOCATE( fion( 3, atoms%nat ) )
        !
        ! fion = atoms%for
        !
        ! DO iss = 1, cdesc%nspin
        !   iss_wfc = iss
        !   IF( force_pairing ) iss_wfc = 1
        !   CALL nl_ionic_forces_v( iss, c0( :, :, :, iss_wfc ), cdesc, atoms, occ(:,:,iss), &
        !     fnl(:,iss), wsg, wnl, eigr)
        ! END DO
        !
        ! DO i = 1, atoms%nat
        !    WRITE(6,*) 'X = ', fion( 1, i ) - atoms%for( 1, i )
        !    WRITE(6,*) 'Y = ', fion( 2, i ) - atoms%for( 2, i )
        !    WRITE(6,*) 'Z = ', fion( 3, i ) - atoms%for( 3, i )
        ! END DO
        !
        ! DEALLOCATE( fion )
        !
      ! END IF
      !
      RETURN
   END FUNCTION nlrh_m

!  ----------------------------------------------


   SUBROUTINE nl_projector_m(c0, cdesc, atoms, fnl, wsg, wnl, eigr)

      !  this routine computes:
      !  fnl: Kleinman-Bylander pseudopotential terms (see nlsm1,nlsm1_kp)
      !

      USE brillouin, ONLY: kpoints, kp
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing, gamma_only
      USE reciprocal_vectors, ONLY: gx, g

      IMPLICIT NONE

      ! ...   declare subroutine arguments
      COMPLEX(DP)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE(atoms_type),      INTENT(INOUT)  :: atoms         ! ions structure
      COMPLEX(DP),    INTENT(INOUT)  :: c0(:,:,:,:)  ! wave functions
      TYPE (wave_descriptor),  INTENT(IN)  :: cdesc  ! wave functions
      REAL(DP),      INTENT(IN)  :: wsg(:,:) ! Kleinman-Bylander inverse
                                              ! denominators
                                              ! <Y phi | V | phi Y>**(-1)
      ! ...   Kleinman-Bylander factors, with and without k points
      TYPE (projector), INTENT(OUT) :: fnl(:,:)
      ! ...   Kleinman-Bylander products <Y phi V | exp(i(k+G) dot r)>
      REAL(DP), INTENT(IN) :: wnl(:,:,:,:)
      ! ...   declare other variables
      INTEGER      :: ik, ispin, ispin_wfc

      DO ispin = 1, cdesc%nspin
        ispin_wfc = ispin
        IF( force_pairing ) ispin_wfc = 1
        DO ik = 1, cdesc%nkl
           CALL nlsm1_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:, :, ik, ispin_wfc), cdesc, g, gx, fnl(ik, ispin))
        END DO
      END DO
      RETURN
   END SUBROUTINE nl_projector_m


!  ----------------------------------------------

   REAL(DP) FUNCTION nl_energy_m(fnl, atoms, occ, wsg)

      !  this routine computes:
      !  enl: nonlocal potential contribution to total energy (see ene_nl,ene_nl_kp)

      ! ... include modules
      USE brillouin, ONLY: kpoints, kp
      USE pseudopotential, ONLY: nspnl
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type

      IMPLICIT NONE

      ! ... declare subroutine arguments
      REAL(DP),    INTENT(IN)    :: occ(:,:,:) ! occupations
      TYPE(atoms_type),  INTENT(IN) :: atoms    ! ions structure
      REAL(DP),      INTENT(IN)    :: wsg(:,:) ! Kleinman-Bylander inverse
                                                ! denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector) :: fnl(:,:) ! ... Kleinman-Bylander factors

      ! ... declare other variables
      REAL(DP)    :: enl
      INTEGER      :: ik, ispin


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



   SUBROUTINE nl_ionic_forces_v( ispin, c0, cdesc, atoms, occ, fnl, wsg, wnl, eigr)

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

      ! ... include modules
      USE pseudopotential, ONLY: nspnl, nsanl
      USE brillouin, ONLY: kpoints, kp
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector, allocate_projector, deallocate_projector
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: gamma_only
      USE reciprocal_vectors, ONLY: gx, g
      USE uspp_param, ONLY: nhm

      IMPLICIT NONE

      COMPLEX(DP)  :: eigr(:,:)          ! exp(i G dot r)
      TYPE(atoms_type),       INTENT(INOUT) :: atoms ! ions structure
      COMPLEX(DP),            INTENT(INOUT) :: c0(:,:,:)     ! wave functions
      TYPE (wave_descriptor), INTENT(IN)    :: cdesc     ! wave functions desc.
      REAL(DP),            INTENT(IN)    :: occ(:,:)       ! occupations
      REAL(DP),              INTENT(IN)    :: wsg(:,:)  ! KB inverse denominators <Y phi | V | phi Y>**(-1)
      TYPE (projector),       INTENT(INOUT) :: fnl(:)   ! KB factors
      REAL(DP),              INTENT(IN) :: wnl(:,:,:,:)    ! KB products <Y phi V | exp(i(k+G) dot r)>
      INTEGER, INTENT(IN) :: ispin

      ! ... declare other variables
      TYPE (projector) :: dfnl
      COMPLEX(DP) :: ctmp
      REAL(DP)    :: temp, fac, tt, fsum, enl
      INTEGER      :: me, ib, ia, k, isa, igh, is, ik


      ! ... compute nonlocal contribution to forces on ions
      KAPPA: DO ik = 1, cdesc%nkl
        CARTE: DO k = 1, 3  ! x,y,z directions
          CALL allocate_projector( dfnl, nsanl, cdesc%nbl( ispin ), nhm, cdesc%gamma )
          CALL nlsm2_s( ispin, wnl(:,:,:,ik), atoms, eigr, c0(:,:,ik), cdesc, g, gx, dfnl, k)
          CHANN: DO igh = 1, nhm
            BANDE: DO ib = 1, cdesc%nbl( ispin )
              isa=0
              SPECS: DO is = 1, nspnl
                temp = 2.d0 * wsg(igh,is) * occ(ib,ik) * kp%weight(ik)
                IF( fnl(ik)%gamma_symmetry ) THEN
                  DO ia = 1, atoms%na(is)
                    isa = isa + 1
                    tt = dfnl%r(isa,igh,ib) * fnl(ik)%r(isa,igh,ib)
                    atoms%for(k,isa) = atoms%for(k,isa) - tt * temp
                    ! WRITE( 6, * ) 'DD tt, temp = ', igh, ib, is, ia, tt, temp
                  END DO
                ELSE
                  DO ia = 1, atoms%na(is)
                    isa = isa+1
                    tt = DBLE( CONJG( dfnl%c(isa,igh,ib) ) * fnl(ik)%c(isa,igh,ib) )
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

      ! ... declare modules

      USE pseudopotential, ONLY: nspnl
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: intra_image_comm
      USE atoms_type_module, ONLY: atoms_type
      USE uspp_param, only: nh, lmaxkb
      USE uspp, only: nhtol, nhtolm, indv
      USE uspp, only: beta
      USE cell_base, only: omega
      USE constants, only: pi

      IMPLICIT NONE 

      ! ... declare subroutine arguments
      INTEGER, INTENT(IN) :: ispin
      COMPLEX(DP),      INTENT(INOUT) :: c(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL(DP), INTENT(IN)           :: wnl(:,:,:)
      REAL(DP)                       :: g2(:), gx(:,:)
      TYPE (projector), INTENT(OUT)   :: fnl
      COMPLEX(DP) :: eigr(:,:)

      ! ... declare other variables
      INTEGER :: is, igh, isa, iss, ia, ig, nb, iy
      INTEGER :: l, ll, m, ngw, lda, ldw, ldf, ih, iv
      REAL(DP), ALLOCATABLE :: gwork(:,:)
      REAL(DP) :: ftmp
      COMPLEX(DP), ALLOCATABLE :: gxtmp(:)
      COMPLEX(DP), ALLOCATABLE :: auxc(:,:)
      COMPLEX(DP), PARAMETER :: ONE  = (1.0d0,0.0d0)
      COMPLEX(DP), PARAMETER :: ZERO = (0.0d0,0.0d0)
      ! ... i^l
      COMPLEX(DP), PARAMETER :: csign(0:3) = (/ (1.0d0, 0.0d0), &
          (0.0d0,1.0d0), (-1.0d0,0.0d0), (0.0d0,-1.0d0) /)

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

      !      write(6,*) 'debug = ', omega
      ftmp = sqrt( 4.0d0 * ( 4.0d0 * pi )**2 / omega )
      !      write(6,*) 'debug = ', ftmp
      !      beta(ig,iv,is) = ftmp * wnl( ig, iv, is) * ylm( ig, iy )
      !      wsg = ftmp**2 * dvan
      !      bec = ftmp * fnl

      ALLOCATE( gwork( ngw, (lmaxkb+1)**2 ), gxtmp(ngw) )

      ! ... angular momentum l = 0
      ! ... orbital: s
      ! ... angular momentum l = 1
      ! ... orbitals: p_x, p_y, p_z
      ! ... angular momentum l = 2
      ! ... orbitals: d_{z^2}, d_{x^2-y^2}, d_{xy}, d_{yz}, d_{zx}

      CALL ylmr2( (lmaxkb+1)**2, ngw, gx, g2, gwork )

      iss = 1
      DO is = 1, nspnl
        ALLOCATE( auxc( ngw, atoms%na(is) ) )
        DO ih = 1, nh( is )
          iy  = nhtolm( ih, is )
          iv  = indv  ( ih, is )
          ll  = nhtol ( ih, is ) + 1
          l   = ll - 1
          igh = ih
          gxtmp(1:ngw) = csign(l) * wnl(1:ngw, iv, is) * gwork(1:ngw, iy )
          IF( cdesc%gamma .AND. cdesc%gzero ) gxtmp(1) = gxtmp(1) * 0.5d0
          DO ia = 1, atoms%na(is)
            auxc(1:ngw,ia) = gxtmp(1:ngw) * eigr(1:ngw,iss+ia-1)
          END DO
          IF ( cdesc%gamma ) THEN
            CALL DGEMM( 'T', 'N', atoms%na(is), nb, 2*ngw, 1.0d0, &
                 auxc(1,1), lda, c(1,1), ldw, 0.0d0, fnl%r(iss,igh,1), ldf)
          ELSE
            CALL ZGEMM( 'C', 'N', atoms%na(is), nb, ngw, one, auxc(1,1), lda, &
                 c(1,1), ldw, zero, fnl%c(iss,igh,1), ldf )
          END IF
        END DO
        iss = iss + atoms%na(is)
        DEALLOCATE(auxc)
      END DO
      DEALLOCATE(gwork, gxtmp)

      ! ... since G vectors only span half space, multiply results by two
      IF ( cdesc%gamma ) THEN
        CALL DSCAL( size( fnl%r ), 2.0d0, fnl%r(1,1,1), 1 )
        CALL mp_sum( fnl%r, intra_image_comm )
      ELSE
        CALL mp_sum( fnl%c, intra_image_comm )
      END IF

      RETURN
   END SUBROUTINE nlsm1_s


!  ----------------------------------------------


   SUBROUTINE nlsm2_s( ispin, wnl, atoms, eigr, c, cdesc, g2, gx, dfnl, kk)

      !  this routine computes the derivatives of the Kleinman-Bylander
      !  factors fnl, to be used for Hellmann-Feynman forces evaluation
      !

      ! ... declare modules
      USE pseudopotential, ONLY: nspnl
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE cell_base, ONLY: tpiba
      USE uspp_param, only: nh, lmaxkb
      USE uspp, only: nhtol, nhtolm, indv

      IMPLICIT   NONE

      ! ... declare subroutine arguments
      COMPLEX(DP), INTENT(IN) :: c(:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (projector), INTENT(OUT) :: dfnl
      TYPE(atoms_type), INTENT(INOUT) :: atoms ! ions structure
      REAL(DP), INTENT(IN) :: wnl(:,:,:), g2(:), gx(:,:)
      COMPLEX(DP), INTENT(IN) :: eigr(:,:)
      INTEGER, INTENT(IN) :: kk
      INTEGER, INTENT(IN) :: ispin

      ! ... declare other variables
      REAL(DP), ALLOCATABLE :: gwork(:,:)
      INTEGER :: is, ia, igh, isa, ig, iss, ll, l, m, ngw, nb, lda, ldw, ldf
      INTEGER :: iy, ih, iv
      COMPLEX(DP), ALLOCATABLE :: auxc(:,:), gxtmp(:)
      COMPLEX(DP), PARAMETER :: ONE  = (1.0d0,0.0d0), ZERO = (0.0d0,0.0d0)
      ! ... (-i) * i^l
      COMPLEX(DP), PARAMETER :: csign(0:3) = (/ (0.0d0,-1.0d0), &
        (1.0d0,0.0d0), (0.0d0,1.0d0), (-1.0d0,0.0d0) /)


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

      ALLOCATE(gwork(ngw, (lmaxkb+1)**2 ), gxtmp(ngw))

      CALL ylmr2( (lmaxkb+1)**2, ngw, gx, g2, gwork )
      !
      DO iy = 1, (lmaxkb+1)**2 
        gwork(1:ngw,iy) = tpiba * gx(kk,1:ngw) * gwork(1:ngw,iy)
      END DO

      iss = 1
      SPECS: DO is = 1, nspnl
        ALLOCATE(auxc(ngw,atoms%na(is)))
        LM: DO ih = 1, nh( is )
          iv  = indv  ( ih, is )
          iy  = nhtolm( ih, is )
          ll  = nhtol ( ih, is ) + 1
          l   = ll - 1
          igh = ih
          ! write( 6, * ) 'DEBUG = ', SUM( wnl( :, iv, is ) ), SUM( gwork( :, iy ) )
          gxtmp(1:ngw) = csign(l) * wnl(1:ngw,iv,is) * gwork(1:ngw,iy)
          DO ia = 1, atoms%na(is)
            auxc(1:ngw,ia) = gxtmp(1:ngw) * eigr(1:ngw,iss + ia - 1)
          END DO
          IF( cdesc%gamma ) THEN
             CALL DGEMM('T', 'N', atoms%na(is), nb, 2*ngw, 1.0d0, auxc(1,1), 2*ngw, &
                c(1,1), ldw, 0.0d0, dfnl%r(iss,igh,1), ldf)
           ELSE
             CALL ZGEMM('C', 'N', atoms%na(is), nb, ngw, one, auxc(1,1), ngw, &
                c(1,1), ldw, zero, dfnl%c(iss,igh,1), ldf)
           END IF
        END DO LM
        DEALLOCATE(auxc)
        iss = iss + atoms%na(is)
      END DO SPECS
      IF( cdesc%gamma ) CALL DSCAL(size(dfnl%r),2.0d0,dfnl%r(1,1,1),1)
      !write( 6, * ) 'DEBUG ==== ', SUM( dfnl%r )
      DEALLOCATE(gwork, gxtmp)
      RETURN
   END SUBROUTINE nlsm2_s


!----------------------------------------------
 END MODULE nl
!----------------------------------------------
