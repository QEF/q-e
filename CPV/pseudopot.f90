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

!  BEGIN manual

   MODULE pseudopotential

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE readpot(na,nsp,tcc)
!  SUBROUTINE defngh(nsp)
!  INTEGER FUNCTION l2ind( lw, is )
!  SUBROUTINE pseudopotential_setup(nsp,psfile)
!  SUBROUTINE deallocate_pseudopotential
!  SUBROUTINE formf(ht,ps)
!  SUBROUTINE nlin(kp,wnl)
!  SUBROUTINE nlin_stress(wnla)
!  SUBROUTINE pseudo_wave_info(oc_out,nchan_out,mesh_out,dx_out, &
!                              rw_out,rps_out)
!  REAL(dbl) FUNCTION zvpseudo(is)
!  ----------------------------------------------
!  END manual

! ...   declare modules
        USE kinds
        USE parameters, ONLY: cp_lmax
        USE cp_types, ONLY: pseudo, allocate_pseudo
        USE splines, ONLY: spline_data
        USE read_pseudo_module_fpmd, ONLY: nspnl, l2ind

        IMPLICIT NONE
        SAVE

!  declare module-scope variables
        INTEGER :: nsanl   ! number of atoms of the non local species
        INTEGER :: lnlx    ! maximum number of different non local components
        INTEGER :: lmax    ! maximum value of the angular momentum
        INTEGER :: lm1x    ! maximum value of the angular momentum, of non local components
        INTEGER :: ngh     ! actual number of spherical harmonics 

        REAL(dbl), ALLOCATABLE :: wsginit(:,:)

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

        LOGICAL :: ts, tp, td, tf, tl(0:(cp_lmax-1))

        PRIVATE

        PUBLIC :: pseudopotential_setup, formf, nlin, nlin_stress
        PUBLIC :: pseudo_wave_info, pseudopotential_init
        PUBLIC :: deallocate_pseudopotential
        PUBLIC :: l2ind
        PUBLIC :: nspnl, tl, lm1x, ngh, nsanl, ts, tp, td, tf, lnlx, lmax

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------


      SUBROUTINE pseudopotential_setup(nsp, tpstab_inp, pstab_size_inp, raggio_inp)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

        USE splines, ONLY: nullify_spline
        USE pseudo_base, ONLY: nlset_base
        USE ions_base, ONLY: zv
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE read_pseudo_module_fpmd, ONLY: ap
        USE splines, ONLY: kill_spline

        INTEGER, INTENT(IN) :: nsp, pstab_size_inp
        LOGICAL, INTENT(IN) :: tpstab_inp
        REAL(dbl), INTENT(IN) :: raggio_inp(:)

        INTEGER :: i, is, il, l, j

!  end of declarations
!  ----------------------------------------------

        DO i = 1, nsp
          ap(i)%raggio = raggio_inp(i)
          IF( ap(i)%raggio <= 0.0d0 ) THEN
            CALL errore(' pseudopotential_setup ',' ion_radius less than 0 ',-1)
          END IF
        END DO

        !  initialize the global array zv containing the 
        !  value of the ionic charges

        DO i = 1, nsp
          zv( i ) = ap( i )%zv
        END DO

! ...   set flags selecting angular momentum
        lnlx = 0
        ts=.FALSE.
        tp=.FALSE.
        td=.FALSE.
        tf=.FALSE.
        tl=.FALSE.
        DO is = 1, nspnl
          lnlx = MAX( lnlx, ap(is)%lnl )
          DO l = 1, ap(is)%lnl
            IF (ap(is)%indl(l).EQ.1) ts = .TRUE.  ! state = s
            IF (ap(is)%indl(l).EQ.2) tp = .TRUE.  ! state = p
            IF (ap(is)%indl(l).EQ.3) td = .TRUE.  ! state = d
            IF (ap(is)%indl(l).EQ.4) tf = .TRUE.  ! state = f
          END DO
        END DO
        tl(0) = ts
        tl(1) = tp
        tl(2) = td
        tl(3) = tf

! ...   count orbitals
        ngh = 0
        IF (ts) ngh = ngh + 1
        IF (tp) ngh = ngh + 3
        IF (td) ngh = ngh + 5
        IF (tf) ngh = ngh + 7

        ! maximum value of l (lmax = 1(s),2(p),3(d),4(f) )
        ! lmax - 1, number of non-local components (Projectors)
        ! one of the component is added to the local part

        lm1x = 0
        IF (ts) lm1x = 0
        IF (tp) lm1x = 1
        IF (td) lm1x = 2
        IF (tf) lm1x = 3
        lmax  = lm1x
        DO is = 1, nspnl
          lmax = MAX( lmax, ( ap(is)%lloc - 1 ) )
        END DO

        tpstab = tpstab_inp
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
          ALLOCATE( wnl_sp(lnlx,nsp))
          !
          IF(ALLOCATED(wnla_sp)) THEN
            DO i = 1, size(wnla_sp,2)
              DO j = 1, size(wnla_sp,1)
                CALL kill_spline(wnla_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnla_sp)
          END IF
          ALLOCATE( wnla_sp(lnlx,nsp))
          !
          DO is = 1, nsp
            CALL nullify_spline( vps_sp( is ) )
            CALL nullify_spline( dvps_sp( is ) )
            CALL nullify_spline( rhoc1_sp( is ) )
            CALL nullify_spline( rhocp_sp( is ) )
            DO il = 1, lnlx
              CALL nullify_spline( wnl_sp( il, is ) )
              CALL nullify_spline( wnla_sp( il, is ) )
            END DO
          END DO
          !
          tpstab_first = .TRUE.
          pstab_size   = pstab_size_inp
          !
        END IF

        IF( ALLOCATED( wsginit ) ) DEALLOCATE( wsginit )
        ALLOCATE(wsginit(ngh,nsp))
        wsginit = 0.0d0
        DO is = 1, nspnl
          CALL nlset_base(ap(is), wsginit(:,is))
        END DO

        RETURN
      END SUBROUTINE pseudopotential_setup

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE pseudopotential_init( ps, na, nsp, kp )

! ...   declare modules
        USE brillouin, ONLY: kpoints
        USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
        USE read_pseudo_module_fpmd, ONLY: ap
        USE gvecp, ONLY: ngm
        USE gvecw, ONLY: ngw
        
        IMPLICIT NONE

        TYPE (pseudo) :: ps
        TYPE (kpoints), INTENT(IN) :: kp
        INTEGER, INTENT(IN) :: na(:), nsp

        LOGICAL :: tcc(nsp)
        INTEGER :: i

! ...     Calculate the number of atoms with non local pseudopotentials
          nsanl = SUM( na(1:nspnl) )

          tcc = .FALSE.
          WHERE ( ap(1:nsp)%tnlcc ) tcc = .TRUE.

          CALL allocate_pseudo(ps, nsp, ngm, ngw, kp%nkpt, lnlx, ngh, tcc)

          ps%ap => ap

        RETURN
      END SUBROUTINE pseudopotential_init

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE deallocate_pseudopotential

          USE splines, ONLY: kill_spline

          INTEGER :: i, j

          IF( ALLOCATED( wsginit ) ) DEALLOCATE( wsginit )
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
!      core corrections to the pseudopotential (cc(:)%rhoc1)
!    the derivatives with respect to cell degrees of freedom of:
!      pseudopotential form factors (dvps)
!      core correction to the pseudopotential (cc(:)%rhocp)
!  ----------------------------------------------

! ... declare modules
      USE ions_base, ONLY: na, nsp
      USE cell_module, ONLY: boxdimensions
      USE cell_base, ONLY: tpiba2
      USE brillouin, ONLY: kpoints
      USE pseudotab_base, ONLY: chkpstab
      USE pseudotab_base, ONLY: formftab_base
      USE pseudotab_base, ONLY: corecortab_base
      USE pseudo_base, ONLY: formfn
      USE pseudo_base, ONLY: compute_rhops, compute_rhocg
      USE reciprocal_vectors, ONLY: g, ngm

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (pseudo)  ps
      TYPE (kpoints), INTENT(IN) :: kp
      TYPE (boxdimensions), INTENT(IN)    :: ht

! ... declare other variables
      INTEGER is, isc, ist, ig, igh, i
      REAL(dbl)  omega, s1, s2, s3, s4
      REAL(dbl), ALLOCATABLE :: vps(:), dvps(:), vloc(:)

!  end of declarations
!  ----------------------------------------------

      omega = ht%deth

      do is = 1, size(ps%wsg,2)
        do igh = 1, size(ps%wsg,1)
          ps%wsg(igh, is) = wsginit(igh, is) / omega
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

        CALL compute_rhops( ps%rhops(:,is), ps%drhops(:,is), ps%ap(is)%zv, &
             ps%ap(is)%raggio, g, omega, tpiba2, ngm, .true. )

        IF( ps%ap(is)%tnlcc ) THEN
          IF(tpstab) THEN
            CALL corecortab_base(g, ps%rhoc1(:,is), ps%rhocp(:,is), &
                   rhoc1_sp(is), rhocp_sp(is), xgtabmax, omega) 
          ELSE
            ! CALL corecor_base(ps%ap(is), g, ps%rhoc1(:,is), ps%rhocp(:,is), omega)
            CALL compute_rhocg( ps%rhoc1(:,is), ps%rhocp(:,is), ps%ap(is)%rw, &
                   ps%ap(is)%rab, ps%ap(is)%rhoc, g, omega, tpiba2, ps%ap(is)%mesh, ngm, 1 )

          END IF
        END IF

! ...   numeric pseudopotential
        IF(tpstab) THEN
          CALL formftab_base(g, ps%vps(:,is), ps%dvps(:,is), &
               vps_sp(is), dvps_sp(is), xgtabmax, omega )
        ELSE
          !
          ALLOCATE( vloc( ps%ap(is)%mesh ) )
          !
          vloc = ps%ap(is)%vloc * 2.0d0
          CALL formfn( ps%vps(:,is), ps%dvps(:,is), ps%ap(is)%rw, ps%ap(is)%rab, &
                      vloc, ps%ap(is)%zv, ps%ap(is)%raggio, g, omega, &
                      tpiba2, 0.0d0, ps%ap(is)%mesh, ngm, .false., .true. )
          ps%dvps(:,is) = -ps%dvps(:,is)
          !
          DEALLOCATE( vloc )

        END IF

        ! DEBUG
        ! IF( is == 2 ) THEN
        !   DO i = 1, SIZE(ps%vps,1)
        !   ps%dvps(1,is) = 0.0
        !   WRITE( stdout,fmt="(I5,2F14.6)" ) i,ps%vps(i,is),ps%dvps(i,is)
        !   END DO
        ! END IF

        ! WRITE( stdout,fmt="(/,'* FORMF TIMING: ',F18.8)" ) (s2-s1)
        ! WRITE( stdout,fmt="(/,'* FORMF TIMING: ',F18.8)" ) (s3-s2)
        ! WRITE( stdout,fmt="(/,'* FORMF TIMING: ',F18.8)" ) (s4-s3)

      END DO


! ... handle nonlocal part
      CALL nlin(kp, ps%wnl)


      RETURN
 1000 FORMAT(3E16.8)
      END SUBROUTINE formf

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE build_pstab()

        USE ions_base, ONLY: nsp
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
                      vloc, ap(is)%zv, ap(is)%raggio, xgtab, 1.0d0, &
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

            DO l = 1, ap(is)%lnl
              CALL allocate_spline(  wnl_sp(l,is), pstab_size, xgmin, xgmax )
              CALL allocate_spline( wnla_sp(l,is), pstab_size, xgmin, xgmax )
            END DO

            ALLOCATE( fintl( SIZE( xgtab ), SIZE( wnl_sp, 1) ) )
            CALL nlin_base(ap(is), xgtab(:), fintl)
            DO l = 1, ap(is)%lnl
              wnl_sp(l,is)%y = fintl(:,l)
            END DO
            CALL nlin_stress_base(ap(is), xgtab, fintl)
            DO l = 1, ap(is)%lnl
              wnla_sp(l,is)%y = fintl(:,l)
            END DO
            DEALLOCATE(fintl)

            DO l = 1, ap(is)%lnl
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
              CALL nlintab_base(g, wnl(:,:,is,ik), ap(is)%lnl, wnl_sp(:,is), xgtabmax)
            ELSE
              CALL nlin_base(ap(is), g, wnl(:,:,is,ik))
            END IF
          ELSE
            IF( tpstab ) THEN
              CALL nlintab_base(gk_l(:,ik), wnl(:,:,is,ik), ap(is)%lnl, wnl_sp(:,is), xgtabmax)
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
          CALL nlintab_base(g, wnla(:,:,is), ap(is)%lnl, wnla_sp(:,is), xgtabmax)
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

!  ----------------------------------------------
   END MODULE pseudopotential
!  ----------------------------------------------

