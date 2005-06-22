!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
!  Last modified: Tue Nov  2 08:03:11 MET 1999
!  ----------------------------------------------

#include "f_defs.h"

   MODULE pseudopotential


! ...   declare modules
        USE kinds
        USE parameters, ONLY: cp_lmax
        USE cp_types, ONLY: pseudo, allocate_pseudo
        USE splines, ONLY: spline_data
        USE read_pseudo_module_fpmd, ONLY: nspnl

        IMPLICIT NONE
        SAVE

!  declare module-scope variables
        INTEGER :: nsanl   ! number of atoms of the non local species
        INTEGER :: nbetax  ! number of atoms of the non local species

        TYPE (spline_data), ALLOCATABLE ::  vps_sp(:)
        TYPE (spline_data), ALLOCATABLE :: dvps_sp(:)
        TYPE (spline_data), ALLOCATABLE ::  wnl_sp(:,:)
        TYPE (spline_data), ALLOCATABLE :: wnla_sp(:,:)
        TYPE (spline_data), ALLOCATABLE :: rhoc1_sp(:)
        TYPE (spline_data), ALLOCATABLE :: rhocp_sp(:)
        REAL(dbl), ALLOCATABLE :: xgtab(:)
        REAL(dbl)              :: xgtabmax 
        INTEGER :: pstab_size

        LOGICAL :: tpstab, tpstab_first

        PRIVATE

        PUBLIC :: pseudopotential_setup, formf, nlin, nlin_stress
        PUBLIC :: pseudo_wave_info, pseudopotential_init
        PUBLIC :: deallocate_pseudopotential
        PUBLIC :: nspnl, nsanl
        PUBLIC :: pseudopotential_initval, pseudopotential_indexes
        PUBLIC :: compute_dvan, compute_betagx, compute_qradx

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  ----------------------------------------------


   SUBROUTINE pseudopotential_setup( nsp, tpstab_inp, pstab_size_inp )

        USE ions_base, ONLY: zv
        USE pseudo_types, ONLY: pseudo_ncpp
        USE read_pseudo_module_fpmd, ONLY: ap

        INTEGER, INTENT(IN) :: nsp, pstab_size_inp
        LOGICAL, INTENT(IN) :: tpstab_inp

        INTEGER :: i

        !  initialize the global array zv containing the 
        !  value of the ionic charges

        DO i = 1, nsp
          zv( i ) = ap( i )%zv
        END DO

        !  set the sizes for the spline tables

        tpstab       = tpstab_inp
        pstab_size   = pstab_size_inp

     RETURN
   END SUBROUTINE pseudopotential_setup



!  ----------------------------------------------



   SUBROUTINE pseudopotential_initval()

        USE splines, ONLY: nullify_spline
        USE ions_base, ONLY: zv, nsp
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE read_pseudo_module_fpmd, ONLY: ap
        USE splines, ONLY: kill_spline
        USE constants, ONLY: pi
        use uspp_param, only: nbeta

        INTEGER :: i, is, il, j

        ! ...   set flags selecting angular momentum

        nbetax = MAXVAL( nbeta( 1:nsp ) )

        IF( tpstab ) THEN
          !
          IF( ALLOCATED( vps_sp ) ) THEN
            DO i = 1, size(vps_sp)
              CALL kill_spline( vps_sp(i), 'a' )
            END DO
            DEALLOCATE( vps_sp )
          END IF
          ALLOCATE( vps_sp(nsp))
          !
          IF(ALLOCATED(dvps_sp)) THEN
            DO i = 1, size(dvps_sp)
              CALL kill_spline(dvps_sp(i),'a')
            END DO
            DEALLOCATE(dvps_sp)
          END IF
          ALLOCATE( dvps_sp(nsp))
          !
          IF(ALLOCATED(rhoc1_sp)) THEN
            DO i = 1, size(rhoc1_sp)
              CALL kill_spline(rhoc1_sp(i),'a')
            END DO
            DEALLOCATE(rhoc1_sp)
          END IF
          ALLOCATE( rhoc1_sp(nsp))
          !
          IF(ALLOCATED(rhocp_sp)) THEN
            DO i = 1, size(rhocp_sp)
              CALL kill_spline(rhocp_sp(i),'a')
            END DO
            DEALLOCATE(rhocp_sp)
          END IF
          ALLOCATE( rhocp_sp(nsp))
          !
          IF(ALLOCATED(wnl_sp)) THEN
            DO i = 1, size(wnl_sp,2)
              DO j = 1, size(wnl_sp,1)
                CALL kill_spline(wnl_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnl_sp)
          END IF
          ALLOCATE( wnl_sp(nbetax,nsp))
          !
          IF(ALLOCATED(wnla_sp)) THEN
            DO i = 1, size(wnla_sp,2)
              DO j = 1, size(wnla_sp,1)
                CALL kill_spline(wnla_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnla_sp)
          END IF
          ALLOCATE( wnla_sp(nbetax,nsp))
          !
          DO is = 1, nsp
            CALL nullify_spline( vps_sp( is ) )
            CALL nullify_spline( dvps_sp( is ) )
            CALL nullify_spline( rhoc1_sp( is ) )
            CALL nullify_spline( rhocp_sp( is ) )
            DO il = 1, nbetax
              CALL nullify_spline( wnl_sp( il, is ) )
              CALL nullify_spline( wnla_sp( il, is ) )
            END DO
          END DO
          !
          tpstab_first = .TRUE.
          !
        END IF

        call compute_dvan()

     RETURN
   END SUBROUTINE pseudopotential_initval




!  ----------------------------------------------


!     
!     calculate array  dvan(iv,jv,is)
!    
!  rw**2 * vrps   = [ ( Vpsnl(r) - Vpsloc(r) )* Rps(r) * r^2 ]
!                 = [ DVpsnl(r) * Rps(r) * r^2 ]
!  dion           = (2l+1) / < Rps(r) | DVpsnl(r) | Rps(r) >


   SUBROUTINE compute_dvan()
     use uspp,       only: dvan, nhtolm, indv
     use uspp_param, only: nhm, nh, dion
     use ions_base,  only: nsp
     use atom,       only: numeric
     implicit none
     integer :: is, iv, jv
     real(dbl) :: fac
     !
     if( allocated( dvan ) ) deallocate( dvan )
     allocate( dvan( nhm, nhm, nsp ) )
     dvan(:,:,:) =0.d0
     !
     do is = 1, nsp
       if ( .not. numeric( is ) ) then
         fac = 1.0d0
       else
         !     fac converts ry to hartree
         fac = 0.5d0
       end if
       do iv=1,nh(is)
         do jv=1,nh(is)
           if ( nhtolm(iv,is) == nhtolm(jv,is) ) then
             dvan( iv, jv, is ) = fac * dion( indv(iv,is), indv(jv,is), is )
           endif
         end do
       end do
     end do
     RETURN
   END SUBROUTINE compute_dvan



!  ----------------------------------------------



   SUBROUTINE pseudopotential_indexes()

      use parameters, only: lmaxx    !
      use ions_base,  only: nsp, &   !  number of specie
                            na       !  number of atoms for each specie
      use cvan,       only: ish      !
      use uspp,       only: nkb, &   !
                            nkbus    !
      use core,       only: nlcc_any !
      use uspp_param, only: nbeta,  &!
                            lmaxkb, &!
                            lll,    &!
                            nhm,    &!
                            nh,     &!
                            tvanp,  &!
                            nqlc,   &!
                            lmaxq    !
      use uspp,       only: nhtol,  &!
                            nhtolm, &!
                            indv     !
      use atom,       only: nlcc     !


      IMPLICIT NONE
     
      INTEGER :: is, iv, ind, il, lm
      !     ------------------------------------------------------------------
      !     find  number of beta functions per species, max dimensions,
      !     total number of beta functions (all and Vanderbilt only)
      !     ------------------------------------------------------------------
      lmaxkb   = -1
      nhm      =  0
      nkb      =  0
      nkbus    =  0
      nlcc_any = .false.
      !
      do is = 1, nsp
         ind = 0
         do iv = 1, nbeta(is)
            lmaxkb = max( lmaxkb, lll( iv, is ) )
            ind = ind + 2 * lll( iv, is ) + 1
         end do
         nh(is) = ind
         nhm    = max( nhm, nh(is) )
         ish(is)=nkb
         nkb = nkb + na(is) * nh(is)
         if( tvanp(is) ) nkbus = nkbus + na(is) * nh(is)
         nlcc_any = nlcc_any .OR. nlcc(is)
      end do
      if (lmaxkb > lmaxx) call errore('nlinit ',' l > lmax ',lmaxkb)
      lmaxq = 2*lmaxkb + 1
      !
      ! the following prevents an out-of-bound error: nqlc(is)=2*lmax+1
      ! but in some versions of the PP files lmax is not set to the maximum
      ! l of the beta functions but includes the l of the local potential
      !
      do is=1,nsp
         nqlc(is) = MIN ( nqlc(is), lmaxq )
      end do
      if (nkb <= 0) call errore('nlinit ','not implemented ?',nkb)

      if( allocated( nhtol ) ) deallocate( nhtol )
      if( allocated( indv  ) ) deallocate( indv )
      if( allocated( nhtolm  ) ) deallocate( nhtolm )
      !
      allocate(nhtol(nhm,nsp))
      allocate(indv (nhm,nsp))
      allocate(nhtolm(nhm,nsp))

      !     ------------------------------------------------------------------
      !     definition of indices nhtol, indv, nhtolm
      !     ------------------------------------------------------------------
      !
      do is = 1, nsp
         ind = 0
         do iv = 1, nbeta(is)
            lm = lll(iv,is)**2
            do il = 1, 2*lll( iv, is ) + 1
               lm = lm + 1
               ind = ind + 1
               nhtolm( ind, is ) = lm
               nhtol( ind, is ) = lll( iv, is )
               indv( ind, is ) = iv
            end do
         end do
      end do


      RETURN
   END SUBROUTINE


!  ----------------------------------------------

      SUBROUTINE pseudopotential_init( ps, na, nsp, kp )

! ...   declare modules
        USE brillouin, ONLY: kpoints
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE local_pseudo, ONLY: allocate_local_pseudo
        USE read_pseudo_module_fpmd, ONLY: ap
        USE gvecp, ONLY: ngm
        USE gvecw, ONLY: ngw
        USE core,  ONLY: nlcc_any, rhoc, drhoc
        use uspp_param, only: lmaxkb, nhm
        
        IMPLICIT NONE

        TYPE (pseudo) :: ps
        TYPE (kpoints), INTENT(IN) :: kp
        INTEGER, INTENT(IN) :: na(:), nsp

        INTEGER :: i

! ...     Calculate the number of atoms with non local pseudopotentials
          nsanl = SUM( na(1:nspnl) )

          ! WRITE( *, * ) 'DEBUG: nlcc_any = ', nlcc_any
          IF( nlcc_any ) THEN
            ! WRITE( *, * ) 'DEBUG: allocating rhoc = ', ngm, nsp
            ALLOCATE( rhoc( ngm, nsp ) )
            ALLOCATE( drhoc( ngm, nsp ) )
          END IF
          CALL allocate_local_pseudo( ngm, nsp )
          CALL allocate_pseudo(ps, nsp, ngw, nbetax, nhm, kp%nkpt)

          ps%ap => ap

        RETURN
      END SUBROUTINE pseudopotential_init

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE deallocate_pseudopotential

          USE splines, ONLY: kill_spline
          USE uspp,    ONLY: dvan

          INTEGER :: i, j

          IF( ALLOCATED( dvan ) ) DEALLOCATE( dvan )
          IF( ALLOCATED( vps_sp  ) ) THEN
            DO i = 1, size(vps_sp)
              CALL kill_spline(vps_sp(i),'a')
            END DO
            DEALLOCATE(vps_sp)
          END IF
          IF( ALLOCATED(dvps_sp) ) THEN
            DO i = 1, size(dvps_sp)
              CALL kill_spline(dvps_sp(i),'a')
            END DO
            DEALLOCATE(dvps_sp)
          END IF
          IF( ALLOCATED(rhoc1_sp) ) THEN
            DO i = 1, size(rhoc1_sp)
              CALL kill_spline(rhoc1_sp(i),'a')
            END DO
            DEALLOCATE(rhoc1_sp)
          END IF
          IF( ALLOCATED(rhocp_sp) ) THEN
            DO i = 1, size(rhocp_sp)
              CALL kill_spline(rhocp_sp(i),'a')
            END DO
            DEALLOCATE(rhocp_sp)
          END IF
          IF( ALLOCATED(wnl_sp) ) THEN
            DO i = 1, size(wnl_sp,2)
              DO j = 1, size(wnl_sp,1)
                CALL kill_spline(wnl_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnl_sp)
          END IF
          IF( ALLOCATED(wnla_sp) ) THEN
            DO i = 1, size(wnla_sp,2)
              DO j = 1, size(wnla_sp,1)
                CALL kill_spline(wnla_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnla_sp)
          END IF
          IF( ALLOCATED( xgtab ) )  DEALLOCATE( xgtab )

          RETURN
        END SUBROUTINE deallocate_pseudopotential


!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE formf(ht, kp, ps)

!  this routine computes:
!    the form factors of:
!      pseudopotential (vps)
!      ionic pseudocharge (rhops)
!      core corrections to the pseudopotential (rhoc)
!    the derivatives with respect to cell degrees of freedom of:
!      pseudopotential form factors (dvps)
!      core correction to the pseudopotential (drhoc)
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: na, nsp, rcmax
      USE cell_module, ONLY: boxdimensions
      USE cell_base, ONLY: tpiba2
      USE brillouin, ONLY: kpoints
      USE pseudotab_base, ONLY: chkpstab
      USE pseudotab_base, ONLY: formftab_base
      USE pseudotab_base, ONLY: corecortab_base
      USE pseudo_base, ONLY: formfn
      USE pseudo_base, ONLY: compute_rhops, compute_rhocg
      USE reciprocal_vectors, ONLY: g, ngm
      USE local_pseudo, ONLY: vps, dvps, rhops, drhops
      use uspp_param, only: lmaxkb, nhm, nh
      use uspp, only: dvan
      USE constants, ONLY: pi
      USE core, ONLY: rhoc, drhoc

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (pseudo)  ps
      TYPE (kpoints), INTENT(IN) :: kp
      TYPE (boxdimensions), INTENT(IN)    :: ht

! ... declare other variables
      INTEGER is, isc, ist, ig, igh, i, l, iv
      REAL(dbl)  omega, s1, s2, s3, s4
      REAL(dbl), ALLOCATABLE :: vloc(:)

!  end of declarations
!  ----------------------------------------------

      omega = ht%deth

      ps%wsg = 0.0d0
      do is = 1, size( ps%wsg, 2 )
        do igh = 1, nh( is )
          ps%wsg( igh, is) = 4.0d0 * ( 4.0d0 * pi ) ** 2 * dvan( igh, igh, is ) / omega
        end do
      end do

      IF(tpstab .AND. (.NOT.tpstab_first)) THEN
! ...   check the consistency of the tab with respect the g vectors
        tpstab_first = chkpstab(g, xgtabmax)
      END IF

      IF(tpstab .AND. tpstab_first) THEN
! ...   build the pseudopotential tabs
        CALL build_pstab() 
        tpstab_first = .FALSE.
      END IF


! ... handle local part
      DO is = 1, nsp

        CALL compute_rhops( rhops(:,is), drhops(:,is), ps%ap(is)%zv, &
             rcmax(is), g, omega, tpiba2, ngm, .true. )

        IF( ps%ap(is)%tnlcc ) THEN
          IF(tpstab) THEN
            CALL corecortab_base(g, rhoc(:,is), drhoc(:,is), &
                   rhoc1_sp(is), rhocp_sp(is), xgtabmax, omega) 
          ELSE
            CALL compute_rhocg( rhoc(:,is), drhoc(:,is), ps%ap(is)%rw, &
                   ps%ap(is)%rab, ps%ap(is)%rhoc, g, omega, tpiba2, ps%ap(is)%mesh, ngm, 1 )

          END IF
        END IF

! ...   numeric pseudopotential
        IF(tpstab) THEN
          CALL formftab_base(g, vps(:,is), dvps(:,is), &
               vps_sp(is), dvps_sp(is), xgtabmax, omega )
        ELSE
          !
          ALLOCATE( vloc( ps%ap(is)%mesh ) )
          !
          vloc = ps%ap(is)%vloc * 2.0d0
          CALL formfn( vps(:,is), dvps(:,is), ps%ap(is)%rw, ps%ap(is)%rab, &
                      vloc, ps%ap(is)%zv, rcmax(is), g, omega, &
                      tpiba2, 0.0d0, ps%ap(is)%mesh, ngm, .false., .true. )
          dvps(:,is) = -dvps(:,is)
          !
          ! WRITE(6,*) ' DEBUG vloc = ', SUM( vloc( : ) )
          ! WRITE(6,*) ' DEBUG mesh = ', ps%ap(is)%mesh
          ! WRITE(6,*) ' DEBUG rcma = ', rcmax(is)
          ! WRITE(6,*) ' DEBUG zv   = ', ps%ap(is)%zv
          ! WRITE(6,*) ' DEBUG rw   = ', SUM( ps%ap(is)%rw( : ) )
          ! WRITE(6,*) ' DEBUG rab  = ', SUM( ps%ap(is)%rab( : ) )
          ! WRITE(6,*) ' DEBUG dvps = ', SUM( dvps( :, is ) )
          !
          DEALLOCATE( vloc )

        END IF

      END DO


! ... handle nonlocal part
      CALL nlin(kp, ps%wnl)


      RETURN
 1000 FORMAT(3E16.8)
      END SUBROUTINE formf

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE build_pstab()

        USE ions_base, ONLY: nsp, rcmax
        USE constants, ONLY: pi, fpi
        USE cell_base, ONLY: tpiba, tpiba2
        USE splines, ONLY: init_spline, allocate_spline
        USE mp, ONLY: mp_max
        USE mp_global, ONLY: mpime, group, nproc
        USE pseudo_base, ONLY: formfn
        USE pseudo_base, ONLY: nlin_base
        USE pseudo_base, ONLY: nlin_stress_base
        USE pseudo_base, ONLY: compute_rhocg
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE read_pseudo_module_fpmd, ONLY: ap
        USE reciprocal_vectors, ONLY: g

        IMPLICIT NONE

        REAL(dbl), ALLOCATABLE :: fintl(:,:)
        REAL(dbl), ALLOCATABLE :: vloc(:)
        INTEGER :: ig, is, mmax, lloc, nval, l, ll
        REAL(dbl)  :: xg, xgmax, xgmin, dxg, res
        LOGICAL :: tnum, tnlcc

        tnlcc = .false.
        DO is = 1, nsp
          tnlcc = tnlcc .or. ap(is)%tnlcc
        END DO

        IF(.NOT.ALLOCATED(xgtab))     ALLOCATE(xgtab(pstab_size))
        nval = pstab_size

        xgmin = 0.0d0 
        xgmax = tpiba * SQRT( MAXVAL( g ) )  
        CALL mp_max(xgmax, group)
        xgmax = xgmax + (xgmax-xgmin)
        dxg   = (xgmax - xgmin) / REAL(nval-1)
        DO ig = 1, SIZE( xgtab )
          xgtab(ig) = xgmin + REAL(ig-1) * dxg
        END DO
        xgtabmax = xgtab( SIZE( xgtab ) )
        xgtab = xgtab**2 / tpiba**2

        DO is = 1, nsp

          tnlcc  = ap(is)%tnlcc 

          CALL allocate_spline( vps_sp(is), pstab_size, xgmin, xgmax )
          CALL allocate_spline( dvps_sp(is), pstab_size, xgmin, xgmax )

          ALLOCATE( vloc( ap(is)%mesh ) )
          vloc = ap(is)%vloc * 2.0d0
          CALL formfn( vps_sp(is)%y, dvps_sp(is)%y, ap(is)%rw, ap(is)%rab, &
                      vloc, ap(is)%zv, rcmax(is), xgtab, 1.0d0, &
                      tpiba2, 0.0d0, ap(is)%mesh, pstab_size, .false., .true. )
          dvps_sp(is)%y = -dvps_sp(is)%y
          DEALLOCATE( vloc )

          CALL init_spline( vps_sp(is) )
          CALL init_spline( dvps_sp(is) )

          IF(tnlcc) THEN
            CALL allocate_spline( rhoc1_sp(is), pstab_size, xgmin, xgmax )
            CALL allocate_spline( rhocp_sp(is), pstab_size, xgmin, xgmax )
            ! CALL corecor_base(ap(is), xgtab, rhoc1_sp(is)%y, rhocp_sp(is)%y, 1.0d0)
            CALL compute_rhocg( rhoc1_sp(is)%y, rhocp_sp(is)%y, ap(is)%rw, &
                   ap(is)%rab, ap(is)%rhoc, xgtab, 1.0d0, tpiba2, ap(is)%mesh, pstab_size, 1 )
            CALL init_spline( rhoc1_sp(is) )
            CALL init_spline( rhocp_sp(is) )
          END IF

! ...     Initialize Tables for array WNL
! ...     This section entered only if the specie is non local

          NONLOCAL: IF( is <= nspnl ) THEN

            DO l = 1, ap(is)%nbeta
              CALL allocate_spline(  wnl_sp(l,is), pstab_size, xgmin, xgmax )
              CALL allocate_spline( wnla_sp(l,is), pstab_size, xgmin, xgmax )
            END DO

            ALLOCATE( fintl( SIZE( xgtab ), SIZE( wnl_sp, 1) ) )
            CALL nlin_base(ap(is), xgtab(:), fintl)
            DO l = 1, ap(is)%nbeta
              wnl_sp(l,is)%y = fintl(:,l)
            END DO
            CALL nlin_stress_base(ap(is), xgtab, fintl)
            DO l = 1, ap(is)%nbeta
              wnla_sp(l,is)%y = fintl(:,l)
            END DO
            DEALLOCATE(fintl)

            DO l = 1, ap(is)%nbeta
              CALL init_spline( wnl_sp(l,is) )
              CALL init_spline( wnla_sp(l,is) )
            END DO

          END IF NONLOCAL

        END DO

        RETURN
      END SUBROUTINE  build_pstab
!  ----------------------------------------------
!  ----------------------------------------------

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE nlin(kp, wnl)

!  this routine computes the temporary arrays twnl
!  to be used by nlrh and dforce
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: nsp
      USE brillouin, ONLY: kpoints
      USE pseudotab_base, ONLY: nlintab_base
      USE pseudo_base, ONLY: nlin_base
      USE read_pseudo_module_fpmd, ONLY: ap
      USE reciprocal_space_mesh, ONLY: gk_l
      USE reciprocal_vectors, ONLY: g
      USE control_flags, ONLY: gamma_only

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (kpoints), INTENT(IN) :: kp
      REAL(dbl) :: wnl(:,:,:,:)

! ... declare other variables
      INTEGER  is, ik

!  end of declarations
!  ----------------------------------------------

      wnl = 0.0d0
      DO is = 1, nspnl
        DO ik = 1, kp%nkpt
          IF( gamma_only ) THEN
            IF( tpstab ) THEN
              CALL nlintab_base(g, wnl(:,:,is,ik), ap(is)%nbeta, wnl_sp(:,is), xgtabmax)
            ELSE
              CALL nlin_base(ap(is), g, wnl(:,:,is,ik))
            END IF
          ELSE
            IF( tpstab ) THEN
              CALL nlintab_base(gk_l(:,ik), wnl(:,:,is,ik), ap(is)%nbeta, wnl_sp(:,is), xgtabmax)
            ELSE
              CALL nlin_base(ap(is), gk_l(:,ik), wnl(:,:,is,ik))
            END IF
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE nlin

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE nlin_stress(wnla)

!  this routine computes the temporary arrays twnl
!  to be used by nlrh and dforce.
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: nsp
      USE pseudotab_base, ONLY: nlintab_base
      USE pseudo_base, ONLY: nlin_stress_base
      USE read_pseudo_module_fpmd, ONLY: ap
      USE reciprocal_vectors, ONLY: g

      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(dbl) :: wnla(:,:,:)
      INTEGER :: is

!  end of declarations
!  ----------------------------------------------
      wnla = 0.0d0
      DO IS = 1, NSPNL
        IF (tpstab) THEN
          CALL nlintab_base(g, wnla(:,:,is), ap(is)%nbeta, wnla_sp(:,is), xgtabmax)
        ELSE
          CALL nlin_stress_base(ap(is), g, wnla(:,:,is))
        END IF
      END DO
      RETURN
      END SUBROUTINE nlin_stress

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE pseudo_wave_info(oc_out,nchan_out,mesh_out,dx_out, &
                                  rw_out,rps_out)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE read_pseudo_module_fpmd, ONLY: ap

! ...   declare subroutine arguments
        REAL(dbl), INTENT(OUT) :: oc_out(:,:)
        INTEGER, INTENT(OUT) :: nchan_out(:)
        INTEGER, INTENT(OUT) :: mesh_out(:)
        REAL(dbl), INTENT(OUT) :: dx_out(:)
        REAL(dbl), INTENT(OUT) :: rw_out(:,:)
        REAL(dbl), INTENT(OUT) :: rps_out(:,:,:)
        INTEGER i,j,k

!  end of declarations
!  ----------------------------------------------

        DO i = 1, size(ap,1)
          nchan_out(i) = ap(i)%nrps
          mesh_out(i) = ap(i)%mesh
          dx_out(i) = ap(i)%dx
          DO j=1,size(ap(i)%oc)
            oc_out(j,i) = ap(i)%oc(j)
          END DO
          DO j=1,size(ap(i)%rw)
            rw_out(j,i) = ap(i)%rw(j)
          END DO
          DO k=1,size(ap(i)%rps,2)
            DO j=1,size(ap(i)%rps,1)
              rps_out(j,i,k) = ap(i)%rps(j,k)
            END DO
          END DO
        END DO

        RETURN
      END SUBROUTINE pseudo_wave_info


! ----------------------------------------------

! calculation of array  betagx(ig,iv,is)

! ----------------------------------------------


    SUBROUTINE compute_betagx( tpre )
      !
      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: nh, kkbeta, betar, nhm, nbeta
      USE atom,       ONLY: r, numeric, rab
      USE uspp,       ONLY: nhtol, indv
      USE betax,      only: refg, betagx, mmx, dbetagx
      USE cvan,       only: oldvan
      USE qrl_mod,    only: qrl, cmesh
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ltmp, i0, ir
      REAL(dbl), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), jltmp(:)
      REAL(dbl) :: xg, xrg
      !
      IF( ALLOCATED( betagx  ) ) DEALLOCATE( betagx )
      IF( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
      ALLOCATE( betagx ( mmx, nhm, nsp ) )
      ALLOCATE( dbetagx( mmx, nhm, nsp ) )

      !
      do is = 1, nsp
         !
         if ( tpre ) then
            allocate( dfint( kkbeta( is ) ) )
            allocate( djl  ( kkbeta( is ) ) )
            allocate( jltmp( kkbeta( is ) ) )
         end if
         allocate( fint ( kkbeta( is ) ) )
         allocate( jl   ( kkbeta( is ) ) )
         !
         do iv = 1, nh(is)
            !
            l = nhtol(iv,is) + 1
            !
            do il = 1, mmx
               !
               xg = sqrt( refg * (il-1) )
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl )
!
               if( tpre )then
                  !
                  ltmp=l-1
                  !
                  ! r(i0) is the first point such that r(i0) >0
                  !
                  i0 = 1
                  if ( r(1,is) < 1.0d-8 ) i0 = 2
                  ! special case q=0
                  if ( xg < 1.0d-8 ) then
                     if (l == 1) then
                        ! Note that dj_1/dx (x=0) = 1/3
                        jltmp(:) = 1.0d0/3.d0
                     else
                        jltmp(:) = 0.0d0
                     end if
                  else
                     call sph_bes (kkbeta(is)+1-i0, r(i0,is), xg, ltmp-1, jltmp )
                  end if
                  do ir = i0, kkbeta(is)
                     xrg = r(ir,is) * xg
                     djl(ir) = jltmp(ir) * xrg - l * jl(ir)
                  end do
                  if ( i0 == 2 ) djl(1) = djl(2)
                  !
               endif
               !
               !     beta(ir)=r*beta(r)
               !
               do ir = 1, kkbeta(is)
                  fint(ir) = r(ir,is) * betar( ir, indv(iv,is), is ) * jl(ir)
               end do
               if (oldvan(is)) then
                  call herman_skillman_int(kkbeta(is),cmesh(is),fint,betagx(il,iv,is))
               else
                  call simpson_cp90(kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               endif
               ! 
               if(tpre) then
                  do ir = 1, kkbeta(is)
                     dfint(ir) = r(ir,is) * betar( ir, indv(iv,is), is ) * djl(ir)
                  end do
                  if (oldvan(is)) then
                     call herman_skillman_int(kkbeta(is),cmesh(is),dfint,dbetagx(il,iv,is))
                  else
                     call simpson_cp90(kkbeta(is),dfint,rab(1,is),dbetagx(il,iv,is))
                  end if
               endif
               !
            end do
         end do
!
         deallocate(jl)
         deallocate(fint)
         if (tpre) then
            deallocate(jltmp)
            deallocate(djl)
            deallocate(dfint)
         end if
         !
      end do

      RETURN
    END SUBROUTINE compute_betagx


!     ---------------------------------------------------------------
!     calculation of array qradx(igb,iv,jv,is)
!
!       qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
!
!     ---------------------------------------------------------------


    SUBROUTINE compute_qradx( tpre )
      !
      use io_global,  only: stdout
      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: nh, kkbeta, betar, nhm, nqlc, qqq, nbrx, lmaxq, nbeta
      USE atom,       ONLY: r, numeric, rab
      USE uspp,       ONLY: nhtol, indv
      USE betax,      only: refg, qradx, mmx, dqradx
      USE cvan,       only: oldvan, ish, nvb
      USE qrl_mod,    only: qrl, cmesh
      use gvecb,      only: ngb
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ltmp, i0, ir, jv
      REAL(dbl), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), jltmp(:)
      REAL(dbl) :: xg, xrg

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
      ALLOCATE(  qradx( mmx, nbrx, nbrx, lmaxq, nsp ) )
      ALLOCATE( dqradx( mmx, nbrx, nbrx, lmaxq, nsp ) )

      DO is = 1, nvb
         !
         IF ( tpre ) THEN
            ALLOCATE( dfint( kkbeta(is) ) )
            ALLOCATE( djl  ( kkbeta(is) ) )
            ALLOCATE( jltmp( kkbeta(is) ) )
         END IF
         allocate( fint( kkbeta(is) ) )
         allocate( jl  ( kkbeta(is) ) )
         !
         !     qqq and beta are now indexed and taken in the same order
         !     as vanderbilts ppot-code prints them out
         !
         WRITE( stdout,*) ' nlinit  nh(is), ngb, is, kkbeta, lmaxq = ', &
     &        nh(is), ngb, is, kkbeta(is), nqlc(is)
         !
         do l = 1, nqlc( is )
            !
            do il = 1, mmx
               !
               xg = sqrt( refg * (il-1) )
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl)
               !
               if(tpre) then
                  !
                  ltmp = l - 1
                  !
                  ! r(i0) is the first point such that r(i0) >0
                  !
                  i0 = 1
                  if ( r(1,is) < 1.0d-8 ) i0 = 2
                  ! special case q=0
                  if ( xg < 1.0d-8 ) then
                     if (l == 1) then
                        ! Note that dj_1/dx (x=0) = 1/3
                        jltmp(:) = 1.0d0/3.d0
                     else
                        jltmp(:) = 0.0d0
                     end if
                  else
                     call sph_bes (kkbeta(is)+1-i0, r(i0,is), xg, ltmp-1, jltmp )
                  end if
                  do ir = i0, kkbeta(is)
                     xrg = r(ir,is) * xg
                     djl(ir) = jltmp(ir) * xrg - l * jl(ir)
                  end do
                  if (i0.eq.2) djl(1) = djl(2)
               endif
               !
               do iv = 1, nbeta(is)
                  do jv = iv, nbeta(is)
                     !
                     !      note qrl(r)=r^2*q(r)
                     !
                     do ir=1,kkbeta(is)
                        fint(ir)=qrl(ir,iv,jv,l,is)*jl(ir)
                     end do
                     if (oldvan(is)) then
                        call herman_skillman_int(kkbeta(is),cmesh(is),fint,qradx(il,iv,jv,l,is))
                     else
                        call simpson_cp90(kkbeta(is),fint,rab(1,is),qradx(il,iv,jv,l,is))
                     end if
                     !
                     qradx(il,jv,iv,l,is)=qradx(il,iv,jv,l,is)
                     !
                     if( tpre ) then
                        do ir = 1, kkbeta(is)
                           dfint(ir) = qrl(ir,iv,jv,l,is) * djl(ir)
                        end do
                        if ( oldvan(is) ) then
                           call herman_skillman_int(kkbeta(is),cmesh(is),dfint,dqradx(il,iv,jv,l,is))
                        else
                           call simpson_cp90(kkbeta(is),dfint,rab(1,is),dqradx(il,iv,jv,l,is))
                        end if
                     end if
                     !
                  end do
               end do
            end do
         end do
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    qqq '
         do iv=1,nbeta(is)
            WRITE( stdout,'(8f9.4)') (qqq(iv,jv,is),jv=1,nbeta(is))
         end do
         WRITE( stdout,*)
         !
         deallocate(jl)
         deallocate(fint)
         if (tpre) then
            deallocate(jltmp)
            deallocate(djl)
            deallocate(dfint)
         end if
         !
      end do

      RETURN
    END SUBROUTINE compute_qradx


!  ----------------------------------------------
   END MODULE pseudopotential
!  ----------------------------------------------

