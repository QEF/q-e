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

  USE kinds,   ONLY: DP
  USE splines, ONLY: spline_data
  USE betax,   ONLY: mmx
  USE read_pseudo_module_fpmd, ONLY: nspnl

  IMPLICIT NONE
  SAVE

  !  declare module-scope variables

  INTEGER :: nsanl   ! number of atoms of the non local species

  TYPE (spline_data), ALLOCATABLE ::  vps_sp(:)
  TYPE (spline_data), ALLOCATABLE :: dvps_sp(:)
  !
  TYPE (spline_data), ALLOCATABLE ::  wnl_sp(:,:)
  TYPE (spline_data), ALLOCATABLE :: wnla_sp(:,:)
  !
  TYPE (spline_data), ALLOCATABLE :: rhoc1_sp(:)
  TYPE (spline_data), ALLOCATABLE :: rhocp_sp(:)
  !
  REAL(DP), ALLOCATABLE :: xgtab(:)

  LOGICAL               :: tpstab = .TRUE.

  PRIVATE

  PUBLIC :: nlin, nlin_stress
  PUBLIC :: deallocate_pseudopotential
  PUBLIC :: nspnl, nsanl
  PUBLIC :: pseudopotential_indexes
  PUBLIC :: compute_dvan, compute_betagx, compute_qradx
  PUBLIC :: interpolate_beta, interpolate_qradb, exact_beta
  PUBLIC :: rhoc1_sp, rhocp_sp, build_cctab, tpstab, chkpstab
  PUBLIC :: build_pstab, vps_sp, dvps_sp
  PUBLIC :: check_tables

  !  ----------------------------------------------

CONTAINS

  !  ----------------------------------------------



   SUBROUTINE compute_dvan()
     !     
     !     calculate array  dvan(iv,jv,is)
     !    
     !  rw**2 * vrps   = [ ( Vpsnl(r) - Vpsloc(r) )* Rps(r) * r^2 ]
     !                 = [ DVpsnl(r) * Rps(r) * r^2 ]
     !  dion           = (2l+1) / < Rps(r) | DVpsnl(r) | Rps(r) >

     use uspp,       only: dvan, nhtolm, indv
     use uspp_param, only: nhm, nh, dion
     use ions_base,  only: nsp
     use atom,       only: numeric
     implicit none
     integer :: is, iv, jv
     real(DP) :: fac
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



   SUBROUTINE pseudopotential_indexes( nlcc_any )

      use parameters, only: lmaxx    !
      use ions_base,  only: nsp, &   !  number of specie
                            na       !  number of atoms for each specie
      use cvan,       only: ish      !
      use uspp,       only: nkb, &   !
                            nkbus    !
      use uspp_param, only: nbeta,  &!
                            lmaxkb, &!
                            lll,    &!
                            nhm,    &!
                            nbetam, &!
                            nh,     &!
                            tvanp,  &!
                            nqlc,   &!
                            lmaxq    !
      use uspp,       only: nhtol,  &!
                            nhtolm, &!
                            indv     !
      use atom,       only: nlcc     !


      IMPLICIT NONE
     
      LOGICAL, INTENT(OUT) :: nlcc_any
      !
      INTEGER :: is, iv, ind, il, lm
      !     ------------------------------------------------------------------
      !     find  number of beta functions per species, max dimensions,
      !     total number of beta functions (all and Vanderbilt only)
      !     ------------------------------------------------------------------
      lmaxkb   = -1
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
         ish(is)=nkb
         nkb = nkb + na(is) * nh(is)
         if( tvanp(is) ) nkbus = nkbus + na(is) * nh(is)
         nlcc_any = nlcc_any .OR. nlcc(is)
      end do
      nhm    = MAXVAL( nh(1:nsp) )
      nbetam = MAXVAL(nbeta(1:nsp))
      if (lmaxkb > lmaxx) call errore(' pseudopotential_indexes ',' l > lmax ',lmaxkb)
      lmaxq = 2*lmaxkb + 1
      !
      ! the following prevents an out-of-bound error: nqlc(is)=2*lmax+1
      ! but in some versions of the PP files lmax is not set to the maximum
      ! l of the beta functions but includes the l of the local potential
      !
      do is=1,nsp
         nqlc(is) = MIN ( nqlc(is), lmaxq )
      end do
      if (nkb <= 0) call errore(' pseudopotential_indexes ',' not implemented ?',nkb)

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

      ! ...     Calculate the number of atoms with non local pseudopotentials
      !
      nsanl = SUM( na(1:nspnl) )

      RETURN
   END SUBROUTINE


!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE deallocate_pseudopotential

          USE splines,      ONLY: kill_spline
          USE local_pseudo, ONLY: deallocate_local_pseudo
          USE uspp,         ONLY: dvan

          INTEGER :: i, j

          CALL deallocate_local_pseudo()
          !
          IF( ALLOCATED( dvan ) ) DEALLOCATE( dvan )
          IF( ALLOCATED( xgtab ) )  DEALLOCATE( xgtab )
          !
          IF( ALLOCATED( vps_sp  ) ) THEN
            DO i = 1, size(vps_sp)
              CALL kill_spline(vps_sp(i),'a')
            END DO
            DEALLOCATE(vps_sp)
          END IF
          !
          IF( ALLOCATED(dvps_sp) ) THEN
            DO i = 1, size(dvps_sp)
              CALL kill_spline(dvps_sp(i),'a')
            END DO
            DEALLOCATE(dvps_sp)
          END IF
          !
          IF( ALLOCATED(rhoc1_sp) ) THEN
            DO i = 1, size(rhoc1_sp)
              CALL kill_spline(rhoc1_sp(i),'a')
            END DO
            DEALLOCATE(rhoc1_sp)
          END IF
          !
          IF( ALLOCATED(rhocp_sp) ) THEN
            DO i = 1, size(rhocp_sp)
              CALL kill_spline(rhocp_sp(i),'a')
            END DO
            DEALLOCATE(rhocp_sp)
          END IF
          !
          IF( ALLOCATED(wnl_sp) ) THEN
            DO i = 1, size(wnl_sp,2)
              DO j = 1, size(wnl_sp,1)
                CALL kill_spline(wnl_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnl_sp)
          END IF
          !
          IF( ALLOCATED(wnla_sp) ) THEN
            DO i = 1, size(wnla_sp,2)
              DO j = 1, size(wnla_sp,1)
                CALL kill_spline(wnla_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnla_sp)
          END IF

          RETURN
        END SUBROUTINE deallocate_pseudopotential

!  ----------------------------------------------

   LOGICAL FUNCTION chkpstab(hg, xgtabmax)
      !
      USE mp,            ONLY: mp_max
      USE io_global,     ONLY: stdout
      USE cell_base,     ONLY: tpiba
      USE control_flags, ONLY: iprsta
      !
      IMPLICIT none
      !
      REAL(DP), INTENT(IN) :: hg(:)
      REAL(DP), INTENT(IN) :: xgtabmax
      REAL(DP) :: xgmax

      chkpstab = .FALSE.
      !
      xgmax = tpiba * SQRT( MAXVAL( hg ) )
      CALL mp_max(xgmax)
      !
      IF( xgmax > xgtabmax ) THEN
         chkpstab = .TRUE.
         IF( iprsta > 2 ) &
            WRITE( stdout, fmt='(  "CHKPSTAB: recalculate pseudopotential table" )' )
      END IF
      !
      RETURN
   END FUNCTION chkpstab


!  ----------------------------------------------

   SUBROUTINE compute_xgtab( xgmin, xgmax, xgtabmax )
      !
      USE cell_base, ONLY: tpiba, tpiba2
      USE mp, ONLY: mp_max
      USE mp_global, ONLY: mpime, group, nproc
      USE reciprocal_vectors, ONLY: g
      !
      REAL(DP), INTENT(OUT)  :: xgmax, xgmin, xgtabmax
      !
      INTEGER    :: ig, nval
      REAL(DP)  :: xg, dxg, res
      !
      IF( .NOT. ALLOCATED( xgtab ) )     ALLOCATE( xgtab( mmx ) )
      nval = mmx
      !
      xgmin = 0.0d0
      xgmax = tpiba * SQRT( MAXVAL( g ) )
      CALL mp_max(xgmax, group)
      xgmax = xgmax + (xgmax-xgmin)
      dxg   = (xgmax - xgmin) / DBLE(nval-1)
      !
      DO ig = 1, SIZE( xgtab )
         xgtab(ig) = xgmin + DBLE(ig-1) * dxg
      END DO
      !
      xgtabmax = xgtab( SIZE( xgtab ) )
      xgtab = xgtab**2 / tpiba**2
      !
      RETURN
   END SUBROUTINE compute_xgtab

!  ----------------------------------------------

   SUBROUTINE build_pstab( )

      USE atom,          ONLY : mesh, r, rab, numeric
      USE ions_base,     ONLY : nsp, rcmax, zv
      USE cell_base,     ONLY : tpiba, tpiba2
      use bhs,           ONLY : rc1, rc2, wrc2, wrc1, rcl, al, bl, lloc
      USE splines,       ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,   ONLY : formfn, formfa
      USE uspp_param,    only : vloc_at, oldvan
      USE control_flags, only : tpre
      use reciprocal_vectors, ONLY : g, gstart

      IMPLICIT NONE

      INTEGER    :: is, ig
      REAL(DP)  :: xgmax, xgmin
      LOGICAL    :: compute_tab
      REAL(DP)  :: xgtabmax = 0.0d0
      !
      compute_tab = chkpstab( g, xgtabmax ) 
      !
      IF( ALLOCATED( vps_sp ) ) THEN
         !
         IF( .NOT. compute_tab ) return
         !
         DO is = 1, nsp
            CALL kill_spline( vps_sp(is), 'a' )
            CALL kill_spline(dvps_sp(is),'a')
         END DO
         DEALLOCATE( vps_sp )
         DEALLOCATE(dvps_sp)
         !
      END IF
      !
      CALL compute_xgtab( xgmin, xgmax, xgtabmax )
      !
      ALLOCATE( vps_sp(nsp))
      ALLOCATE( dvps_sp(nsp))
      !
      DO is = 1, nsp

         CALL nullify_spline( vps_sp( is ) )
         CALL nullify_spline( dvps_sp( is ) )

         CALL allocate_spline( vps_sp(is), mmx, xgmin, xgmax )
         CALL allocate_spline( dvps_sp(is), mmx, xgmin, xgmax )

         if ( numeric(is) ) then

            call formfn( vps_sp(is)%y, dvps_sp(is)%y, r(:,is), rab(:,is), vloc_at(:,is), &
                         zv(is), rcmax(is), xgtab, 1.0d0, tpiba2, mesh(is), &
                         mmx, oldvan(is), tpre )

         else

            !     bhs pseudopotentials
            !
            call formfa( vps_sp(is)%y, dvps_sp(is)%y, rc1(is), rc2(is), wrc1(is), wrc2(is), &
                         rcl(:,is,lloc(is)), al(:,is,lloc(is)), bl(:,is,lloc(is)),    &
                         zv(is), rcmax(is), xgtab, 1.0d0, tpiba2, mmx, 2 , tpre )

         end if

         ! WRITE( 13, "(3D16.8)" ) ( xgtab(ig), vps_sp(is)%y(ig), dvps_sp(is)%y(ig), ig = 1, mmx )

         CALL init_spline( vps_sp(is) )
         CALL init_spline( dvps_sp(is) )

      END DO

      RETURN
   END SUBROUTINE  build_pstab

!  ----------------------------------------------

   SUBROUTINE build_cctab( )

      USE atom,        ONLY : mesh, r, rab, nlcc, rho_atc
      USE ions_base,   ONLY : nsp, rcmax
      USE cell_base,   ONLY : tpiba, tpiba2
      USE splines,     ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base, ONLY : compute_rhocg
      USE uspp_param,  ONLY : kkbeta
      use reciprocal_vectors, ONLY : g, gstart

      IMPLICIT NONE

      INTEGER   :: is
      REAL(DP) :: xgmax, xgmin
      LOGICAL    :: compute_tab
      REAL(DP)  :: xgtabmax = 0.0d0
      !
      compute_tab = chkpstab( g, xgtabmax )
      !
      IF( ALLOCATED( rhoc1_sp ) ) THEN
         !
         IF( .NOT. compute_tab ) return
         !
         DO is = 1, nsp
            CALL kill_spline(rhoc1_sp(is),'a')
            CALL kill_spline(rhocp_sp(is),'a')
         END DO
         DEALLOCATE(rhoc1_sp)
         DEALLOCATE(rhocp_sp)
         !
      END IF
      !
      CALL compute_xgtab( xgmin, xgmax, xgtabmax )
      !
      ALLOCATE( rhoc1_sp(nsp))
      ALLOCATE( rhocp_sp(nsp))
      !
      DO is = 1, nsp

         CALL nullify_spline( rhoc1_sp( is ) )
         CALL nullify_spline( rhocp_sp( is ) )

         IF( nlcc( is ) ) THEN
            !
            CALL allocate_spline( rhoc1_sp(is), mmx, xgmin, xgmax )
            CALL allocate_spline( rhocp_sp(is), mmx, xgmin, xgmax )
            !
            CALL compute_rhocg( rhoc1_sp(is)%y, rhocp_sp(is)%y, r(:,is), &
                 rab(:,is), rho_atc(:,is), xgtab, 1.0d0, tpiba2, kkbeta(is), mmx, 1 )
            !
            CALL init_spline( rhoc1_sp(is) )
            CALL init_spline( rhocp_sp(is) )
            !
         END IF

      END DO

      RETURN
   END SUBROUTINE  build_cctab


!  ----------------------------------------------


   SUBROUTINE build_nltab( )

        ! ...     Initialize Tables for array WNL

        USE ions_base,          ONLY: nsp, rcmax
        USE cell_base,          ONLY: tpiba, tpiba2
        USE splines,            ONLY: init_spline, allocate_spline, kill_spline, nullify_spline
        USE pseudo_base,        ONLY: nlin_base
        USE pseudo_base,        ONLY: nlin_stress_base
        USE pseudo_types,       ONLY: pseudo_ncpp, pseudo_upf
        USE reciprocal_vectors, ONLY: g, gstart
        USE uspp_param,         ONLY: nbeta, nbetam
        USE read_pseudo_module_fpmd, ONLY: ap

        IMPLICIT NONE

        REAL(DP), ALLOCATABLE :: fintl(:,:)
        INTEGER    :: is, l, i
        REAL(DP)  :: xgmax, xgmin
        LOGICAL    :: compute_tab
        REAL(DP)  :: xgtabmax = 0.0d0

        compute_tab = chkpstab( g, xgtabmax )

        IF( ALLOCATED( wnl_sp ) ) THEN

           IF( .NOT. compute_tab ) return

           DO l = 1, nbetam
              DO i = 1, nspnl
                 CALL kill_spline( wnl_sp( l, i ), 'a' )
                 CALL kill_spline( wnla_sp( l, i ), 'a' )
              END DO
           END DO
           !
           DEALLOCATE( wnl_sp )
           DEALLOCATE( wnla_sp )
           !
        END IF

        CALL compute_xgtab( xgmin, xgmax, xgtabmax )

        ALLOCATE( wnl_sp( nbetam, nspnl ) )
        ALLOCATE( wnla_sp( nbetam, nspnl ) )

        DO is = 1, nspnl

           DO l = 1, nbetam
              CALL nullify_spline(  wnl_sp( l, is ) )
              CALL nullify_spline( wnla_sp( l, is ) )
           END DO

           DO l = 1, nbeta( is )
              CALL allocate_spline(  wnl_sp(l,is), mmx, xgmin, xgmax )
              CALL allocate_spline( wnla_sp(l,is), mmx, xgmin, xgmax )
           END DO

           ALLOCATE( fintl( SIZE( xgtab ), SIZE( wnl_sp, 1) ) )
           !
           CALL nlin_base(ap(is), xgtab(:), fintl)
           !
           DO l = 1, nbeta( is )
              wnl_sp( l, is )%y = fintl(:,l)
           END DO
           !
           CALL nlin_stress_base( ap(is), xgtab, fintl )
           
           DO l = 1, nbeta( is )
               wnla_sp( l, is )%y = fintl(:,l)
           END DO
           !
           DO l = 1, nbeta( is )
              CALL init_spline(  wnl_sp( l, is ) )
              CALL init_spline( wnla_sp( l, is ) )
           END DO

           DEALLOCATE(fintl)

        END DO

        RETURN
   END SUBROUTINE  build_nltab


!  ----------------------------------------------

!  

!  ----------------------------------------------


   SUBROUTINE nlin( wsg, wnl )

      !  this routine computes the temporary arrays twnl
      !  to be used by nlrh and dforce
      !

      USE ions_base,      ONLY: nsp
      USE cell_base,      ONLY: omega, tpiba
      USE brillouin,      ONLY: kpoints, kp
      USE pseudo_base,    ONLY: nlin_base
      USE control_flags,  ONLY: gamma_only
      use uspp,           only: dvan
      use uspp_param,     only: nh, nbeta
      use constants,      only: pi
      USE splines,        ONLY: spline
      USE read_pseudo_module_fpmd, ONLY: ap
      USE reciprocal_vectors,      ONLY: g, gstart

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: wsg( :, : )
      REAL(DP), INTENT(OUT) :: wnl( :, :, :, : )

      ! ... declare other variables
      !
      REAL(DP) :: xg
      INTEGER   :: iv, is, ik, ig, l

      !  end of declarations

      wsg = 0.0d0
      do is = 1, size( wsg, 2 )
        do iv = 1, nh( is )
          wsg( iv, is) = 4.0d0 * ( 4.0d0 * pi ) ** 2 * dvan( iv, iv, is ) / omega
        end do
      end do

      IF( tpstab ) THEN
         !
         CALL build_nltab( )
         !
      END IF
     
      wnl = 0.0d0

      DO is = 1, nspnl
         !
         DO ik = 1, kp%nkpt
            !
            IF( tpstab ) THEN
               !
               DO l = 1, nbeta( is )
                  !
                  IF( gstart == 2 ) THEN
                     wnl(1,l,is,ik) = wnl_sp( l, is )%y(1)
                  END IF
                  !
                  DO ig = gstart, SIZE( wnl, 1 )
                     xg = SQRT( g(ig) ) * tpiba
                     wnl(ig,l,is,ik) = spline( wnl_sp( l, is ), xg )
                  END DO
                  !
               END DO
               !
            ELSE
               !
               CALL nlin_base( ap(is), g, wnl(:,:,is,ik) )
               !
            END IF
            !
         END DO
         !
      END DO

      RETURN
   END SUBROUTINE nlin

!  ----------------------------------------------

   SUBROUTINE nlin_stress( wnla )

      !  this routine computes the temporary arrays wnla
      !  to be used by stress subroutine.
      !
      !  Note: subroutine nlin should be called first
      !  
      USE ions_base,               ONLY: nsp
      USE cell_base,               ONLY: tpiba
      USE pseudo_base,             ONLY: nlin_stress_base
      USE splines,                 ONLY: spline
      USE uspp_param,              ONLY: nbeta
      USE read_pseudo_module_fpmd, ONLY: ap
      USE reciprocal_vectors,      ONLY: g, gstart

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: wnla(:,:,:)
      !
      INTEGER   :: is, l, ig
      REAL(DP) :: xg
      !
      !  end of declarations
      !  
      wnla = 0.0d0
      !
      DO is = 1, nspnl
         !
         IF ( tpstab ) THEN
            !
            DO l = 1, nbeta( is )
               !
               IF( gstart == 2 ) THEN
                  wnla(1,l,is) = wnla_sp( l, is )%y(1)
               END IF
               !
               DO ig = gstart, SIZE( wnla, 1 )
                  xg = SQRT( g(ig) ) * tpiba
                  wnla(ig,l,is) = spline( wnla_sp( l, is ), xg )
               END DO
               !
            END DO
            !
         ELSE
            !
            CALL nlin_stress_base(ap(is), g, wnla(:,:,is))
            !
         END IF
         !
      END DO
      !
      RETURN
   END SUBROUTINE nlin_stress


! ----------------------------------------------

! calculation of array  betagx(ig,iv,is)

! ----------------------------------------------


   SUBROUTINE compute_betagx( tpre )
      !

      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: nh, kkbeta, betar, nhm, nbeta, oldvan
      USE atom,       ONLY: r, numeric, rab
      USE uspp,       ONLY: nhtol, indv
      USE betax,      only: refg, betagx, mmx, dbetagx
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ltmp, i0, ir
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), jltmp(:)
      REAL(DP) :: xg, xrg
      !
      IF( ALLOCATED( betagx  ) ) DEALLOCATE( betagx )
      IF( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
      ALLOCATE( betagx ( mmx, nhm, nsp ) )
      IF (tpre) ALLOCATE( dbetagx( mmx, nhm, nsp ) )

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
                  call herman_skillman_int(kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               else
                  call simpson_cp90(kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               endif
               ! 
               if(tpre) then
                  do ir = 1, kkbeta(is)
                     dfint(ir) = r(ir,is) * betar( ir, indv(iv,is), is ) * djl(ir)
                  end do
                  if (oldvan(is)) then
                     call herman_skillman_int(kkbeta(is),dfint,rab(1,is),dbetagx(il,iv,is))
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
!     calculation of array qradx(igb,iv,jv,is) for interpolation table
!     (symmetric wrt exchange of iv and jv: a single index ijv is used)
!
!       qradx(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
!     
!     ---------------------------------------------------------------


    SUBROUTINE compute_qradx( tpre )
      !
      use io_global,  only: stdout
      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: nh, kkbeta, betar, nhm, nbetam, nqlc, qqq, &
           lmaxq, nbeta, oldvan
      USE atom,       ONLY: r, numeric, rab
      USE uspp,       ONLY: nhtol, indv
      USE betax,      only: refg, qradx, mmx, dqradx
      USE cvan,       only: ish, nvb
      use gvecb,      only: ngb
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ltmp, i0, ir, jv, ijv
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), jltmp(:), &
           qrl(:,:,:)
      REAL(DP) :: xg, xrg

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
      ALLOCATE(  qradx( mmx, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      IF (tpre) ALLOCATE( dqradx( mmx, nbetam*(nbetam+1)/2, lmaxq, nsp ) )

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
         ALLOCATE ( qrl(kkbeta(is), nbeta(is)*(nbeta(is)+1)/2, nqlc(is)) )
         call fill_qrl (is,  SIZE(qrl, 1), SIZE(qrl, 2), SIZE(qrl, 3), qrl)
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
               ijv = 0
               do iv = 1, nbeta(is)
                  do jv = iv, nbeta(is)
                     ijv = ijv + 1
                     !
                     !      note qrl(r)=r^2*q(r)
                     !
                     do ir=1,kkbeta(is)
                        fint(ir) = qrl(ir,ijv,l)*jl(ir)
                     end do
                     if (oldvan(is)) then
                        call herman_skillman_int &
                             (kkbeta(is),fint,rab(1,is),qradx(il,ijv,l,is))
                     else
                        call simpson_cp90 &
                             (kkbeta(is),fint,rab(1,is),qradx(il,ijv,l,is))
                     end if
                     !
                     if( tpre ) then
                        do ir = 1, kkbeta(is)
                           dfint(ir) = qrl(ir,ijv,l) * djl(ir)
                        end do
                        if ( oldvan(is) ) then
                           call herman_skillman_int &
                                (kkbeta(is),dfint,rab(1,is),dqradx(il,ijv,l,is))
                        else
                           call simpson_cp90 &
                                (kkbeta(is),dfint,rab(1,is),dqradx(il,ijv,l,is))
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
         DEALLOCATE ( qrl )
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

! ----------------------------------------------

! check table size against cell variations

! ----------------------------------------------

    LOGICAL FUNCTION check_tables( )
      !
      USE kinds, ONLY : DP
      USE betax, ONLY : refg
      USE mp,    ONLY : mp_max
      USE gvecw, ONLY: ngw
      USE cell_base, ONLY: tpiba2
      USE small_box, ONLY: tpibab
      USE gvecb,     ONLY: gb, ngb
      USE reciprocal_vectors, ONLY: g
      !
      IMPLICIT NONE
      !
      REAL(DP) :: gg, ggb, gmax
      !
      gg  = MAXVAL( g( 1:ngw ) )
      gg  = gg * tpiba2 / refg
      !
      ggb = MAXVAL( gb( 1:ngb ) )
      ggb = ggb * tpibab * tpibab / refg
      !
      gmax = MAX( gg, ggb )
      !
      CALL mp_max( gmax )
      !
      check_tables = .FALSE.
      IF( ( INT( gmax ) + 2 ) > mmx ) check_tables = .TRUE.
      !
      RETURN
    END FUNCTION
    

! ----------------------------------------------

! interpolate array beta(ig,iv,is)

! ----------------------------------------------


    SUBROUTINE interpolate_beta( tpre )

      USE kinds, ONLY : DP
      USE control_flags, only: iprsta
      USE constants, only: pi, fpi
      USE io_global, only: stdout
      USE gvecw, only: ngw
      USE ions_base, only: nsp
      USE reciprocal_vectors, only: g, gx, gstart
      USE uspp_param, only: lmaxq, nqlc, lmaxkb, kkbeta, nbeta, nh
      USE uspp, only: qq, nhtolm, beta
      USE cell_base, only: ainv, omega, tpiba2, tpiba
      USE betax, ONLY : refg, betagx, dbetagx
      USE cdvan, ONLY : dbeta

      LOGICAL, INTENT(IN) :: tpre
 
      REAL(DP), ALLOCATABLE ::  ylm(:,:), dylm(:,:,:,:)
      REAL(DP) :: c, gg, betagl, dbetagl
      INTEGER   :: is, iv, lp, ig, jj, i, j

      ALLOCATE( ylm( ngw, (lmaxkb+1)**2 ) )
      CALL ylmr2 ( (lmaxkb+1)**2, ngw, gx, g, ylm)
      !
      !
      do is = 1, nsp
         !   
         !   calculation of array  beta(ig,iv,is)
         !  
         if( iprsta .ge. 4 ) WRITE( stdout,*)  '  beta  '
         c = fpi / sqrt(omega)
         do iv = 1, nh(is)
            lp = nhtolm( iv, is )
            do ig = 1, ngw
               gg = g( ig ) * tpiba * tpiba / refg
               jj = int( gg ) + 1
               betagl = betagx( jj+1, iv, is ) * ( gg - DBLE(jj-1) ) + betagx( jj, iv, is ) * ( DBLE(jj) - gg )
               beta( ig, iv, is ) = c * ylm( ig, lp ) * betagl
            end do
         end do
      end do

      if (tpre) then
         !
         !     calculation of array dbeta required for stress, variable-cell
         !
         allocate( dylm( ngw, (lmaxkb+1)**2, 3, 3 ) )
         !
         call dylmr2_( (lmaxkb+1)**2, ngw, gx, g, ainv, dylm )
         !
         do is = 1, nsp
            if( iprsta .ge. 4 ) WRITE( stdout,*)  '  dbeta  '
            c = fpi / sqrt(omega)
            do iv = 1, nh(is)
               lp = nhtolm(iv,is)
               betagl = betagx(1,iv,is)
               do i=1,3
                  do j=1,3
                     dbeta(1,iv,is,i,j)=-0.5*beta(1,iv,is)*ainv(j,i)    &
     &                                 +c*dylm(1,lp,i,j)*betagl
                  enddo
               enddo
               do ig=gstart,ngw
                  gg=g(ig)*tpiba*tpiba/refg
                  jj=int(gg)+1
                  betagl = betagx(jj+1,iv,is)*(gg-DBLE(jj-1)) +         &
     &                     betagx(jj,iv,is)*(DBLE(jj)-gg)
                  dbetagl= dbetagx(jj+1,iv,is)*(gg-DBLE(jj-1)) +        &
     &                     dbetagx(jj,iv,is)*(DBLE(jj)-gg)
                  do i=1,3
                     do j=1,3
                        dbeta(ig,iv,is,i,j)=                            &
     &                    -0.5*beta(ig,iv,is)*ainv(j,i)                 &
     &                    +c*dylm(ig,lp,i,j)*betagl                     &
     &                    -c*ylm (ig,lp)*dbetagl*gx(i,ig)/g(ig)         &
     &                    *(gx(1,ig)*ainv(j,1)+                         &
     &                      gx(2,ig)*ainv(j,2)+                         &
     &                      gx(3,ig)*ainv(j,3))
                     end do
                  end do
               end do
            end do
         end do
         !
         deallocate(dylm)
         !
      end if
      !
      deallocate(ylm)

      RETURN
    END SUBROUTINE interpolate_beta


! ----------------------------------------------

! interpolate array qradb(ig,iv,is)

! ----------------------------------------------


    SUBROUTINE interpolate_qradb( tpre )
      !
      use control_flags, only: iprint, iprsta
      use io_global, only: stdout
      use gvecw, only: ngw
      use cell_base, only: ainv
      use cvan, only: nvb
      use uspp, only: qq, nhtolm, beta
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use uspp_param, only: lmaxq, nqlc, lmaxkb, kkbeta, nbeta, nbetam, nh
      use qradb_mod, only: qradb
      use qgb_mod, only: qgb
      use gvecb, only: gb, gxb, ngb
      use small_box,  only: omegab, tpibab
      use dqrad_mod, only: dqrad
      use dqgb_mod, only: dqgb
      USE betax, ONLY: qradx, dqradx, refg, mmx
!
      implicit none

      LOGICAL, INTENT(IN) :: tpre

      integer  is, l, lp, ig, ir, iv, jv, ijv, i,j, jj, ierr
      real(8), allocatable:: fint(:), jl(:), dqradb(:,:,:,:)
      real(8), allocatable:: ylmb(:,:), dylmb(:,:,:,:)
      complex(8), allocatable:: dqgbs(:,:,:)
      real(8) xg, c, betagl, dbetagl, gg
!
!
      allocate( ylmb( ngb, lmaxq*lmaxq ), STAT=ierr )
      IF( ierr  /= 0 ) &
        CALL errore(' interpolate_qradb ', ' cannot allocate ylmb ', 1 )
!
      qradb(:,:,:,:,:) = 0.d0
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

      do is = 1, nvb
         !
         !     calculation of array qradb(igb,iv,jv,is)
         !
         if( iprsta .ge. 4 ) WRITE( stdout,*)  '  qradb  '
         !
         c = fpi / omegab
         !
         ijv=0
         do iv= 1,nbeta(is)
            do jv=iv,nbeta(is)
               ijv=ijv+1
               do ig=1,ngb
                  gg=gb(ig)*tpibab*tpibab/refg
                  jj=int(gg)+1
                  do l=1,nqlc(is)
                     if(jj.ge.mmx) then
                        qradb(ig,iv,jv,l,is)=0.
                     else
                        qradb(ig,iv,jv,l,is)=                           &
     &                       c*qradx(jj+1,ijv,l,is)*(gg-DBLE(jj-1))+  &
     &                       c*qradx(jj,ijv,l,is)*(DBLE(jj)-gg)
                     endif
                     qradb(ig,jv,iv,l,is)=qradb(ig,iv,jv,l,is)
                  enddo
               enddo
            enddo
         enddo
!
!     ---------------------------------------------------------------
!     stocking of qgb(igb,ijv,is) and of qq(iv,jv,is)
!     ---------------------------------------------------------------
         ijv=0
         do iv= 1,nh(is)
            do jv=iv,nh(is)
!
!       compact indices because qgb is symmetric:
!       ivjv:  11 12 13 ... 22 23...
!       ijv :   1  2  3 ...
!
               ijv=ijv+1
               call qvan2b(ngb,iv,jv,is,ylmb,qgb(1,ijv,is) )
!
               qq(iv,jv,is)=omegab*DBLE(qgb(1,ijv,is))
               qq(jv,iv,is)=qq(iv,jv,is)
!
            end do
         end do

      end do
!
      if (tpre) then
!     ---------------------------------------------------------------
!     arrays required for stress calculation, variable-cell dynamics
!     ---------------------------------------------------------------
         allocate(dqradb(ngb,nbetam*(nbetam+1)/2,lmaxq,nsp))
         allocate(dylmb(ngb,lmaxq*lmaxq,3,3))
         allocate(dqgbs(ngb,3,3))
         dqrad(:,:,:,:,:,:,:) = 0.d0
         !
         call dylmr2_(lmaxq*lmaxq, ngb, gxb, gb, ainv, dylmb)
         !
         do is=1,nvb
            !
            ijv=0
            do iv= 1,nbeta(is)
               do jv=iv,nbeta(is)
                  ijv=ijv+1
                  do l=1,nqlc(is)
                     do ig=1,ngb
                        gg=gb(ig)*tpibab*tpibab/refg
                        jj=int(gg)+1
                        if(jj.ge.mmx) then
                           dqradb(ig,ijv,l,is) = 0.
                        else
                           dqradb(ig,ijv,l,is) =  &
                                dqradx(jj+1,ijv,l,is)*(gg-DBLE(jj-1)) +  &
                                dqradx(jj,ijv,l,is)*(DBLE(jj)-gg)
                        endif
                     enddo
                     do i=1,3
                        do j=1,3
                           dqrad(1,iv,jv,l,is,i,j) = &
                                -qradb(1,iv,jv,l,is) * ainv(j,i)
                           dqrad(1,jv,iv,l,is,i,j) = &
                                dqrad(1,iv,jv,l,is,i,j)
                           do ig=2,ngb
                              dqrad(ig,iv,jv,l,is,i,j) =                &
     &                          -qradb(ig,iv,jv,l,is)*ainv(j,i)         &
     &                          -c*dqradb(ig,ijv,l,is)*               &
     &                          gxb(i,ig)/gb(ig)*                       &
     &                          (gxb(1,ig)*ainv(j,1)+                   &
     &                           gxb(2,ig)*ainv(j,2)+                   &
     &                           gxb(3,ig)*ainv(j,3))
                              dqrad(ig,jv,iv,l,is,i,j) =                &
     &                          dqrad(ig,iv,jv,l,is,i,j)
                           enddo
                        enddo
                     enddo
                  end do
               enddo
            enddo
            !
            ijv=0
            !
            do iv= 1,nh(is)
               do jv=iv,nh(is)
                  !
                  !       compact indices because qgb is symmetric:
                  !       ivjv:  11 12 13 ... 22 23...
                  !       ijv :   1  2  3 ...
                  !
                  ijv=ijv+1
                  call dqvan2b(ngb,iv,jv,is,ylmb,dylmb,dqgbs )
                  do i=1,3
                     do j=1,3
                        do ig=1,ngb
                           dqgb(ig,ijv,is,i,j)=dqgbs(ig,i,j)
                        enddo
                     enddo
                  enddo
               end do
            end do
         end do
         deallocate(dqgbs)
         deallocate(dylmb)
         deallocate(dqradb)
      end if
      deallocate(ylmb)

      RETURN
    END SUBROUTINE interpolate_qradb


! ----------------------------------------------

! compute array beta without interpolation

! ----------------------------------------------


    SUBROUTINE exact_beta( tpre )

      USE control_flags, only : iprsta
      USE kinds,         ONLY : DP
      USE constants,     only : pi, fpi
      USE io_global,     only : stdout
      USE gvecw,         only : ngw
      USE ions_base,     only : nsp
      USE uspp_param,    only : lmaxq, nqlc, lmaxkb, kkbeta, nbeta, nh, &
           betar, nhm, oldvan
      USE uspp,          only : qq, nhtolm, beta, nhtol, indv
      USE cell_base,     only : ainv, omega, tpiba2, tpiba
      USE cdvan,         ONLY : dbeta
      USE atom,          ONLY : r, numeric, rab
      USE reciprocal_vectors, only : g, gx, gstart

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: tpre
 
      REAL(DP), ALLOCATABLE ::  ylm(:,:), dylm(:,:,:,:)
      REAL(DP) :: c, gg, betagl, dbetagl
      INTEGER :: is, iv, lp, ig, jj, i, j
      INTEGER :: l, il, ltmp, i0, ir
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), jltmp(:)
      REAL(DP), ALLOCATABLE :: betagx ( :, :, : ), dbetagx( :, :, : )
      REAL(DP) :: xg, xrg

      ALLOCATE( ylm( ngw, (lmaxkb+1)**2 ) )
      ALLOCATE( betagx ( ngw, nhm, nsp ) )
      IF (tpre) ALLOCATE( dbetagx( ngw, nhm, nsp ) )

      CALL ylmr2 ( (lmaxkb+1)**2, ngw, gx, g, ylm)

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
            do il = 1, ngw
               !
               xg = sqrt( g( il ) * tpiba * tpiba )
               call sph_bes (kkbeta(is), r(1,is), xg, l-1, jl )
               !
               if( tpre )then
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
                  call herman_skillman_int(kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               else
                  call simpson_cp90(kkbeta(is),fint,rab(1,is),betagx(il,iv,is))
               endif
               ! 
               if(tpre) then
                  do ir = 1, kkbeta(is)
                     dfint(ir) = r(ir,is) * betar( ir, indv(iv,is), is ) * djl(ir)
                  end do
                  if (oldvan(is)) then
                     call herman_skillman_int(kkbeta(is),dfint,rab(1,is),dbetagx(il,iv,is))
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
      !
      do is = 1, nsp
         !   
         !   calculation of array  beta(ig,iv,is)
         !  
         if( iprsta .ge. 4 ) WRITE( stdout,*)  '  beta  '
         c = fpi / sqrt(omega)
         do iv = 1, nh(is)
            lp = nhtolm( iv, is )
            do ig = 1, ngw
               betagl = betagx( ig, iv, is ) 
               beta( ig, iv, is ) = c * ylm( ig, lp ) * betagl
            end do
         end do
      end do

      if (tpre) then
         !
         !     calculation of array dbeta required for stress, variable-cell
         !
         allocate( dylm( ngw, (lmaxkb+1)**2, 3, 3 ) )
         !
         call dylmr2_( (lmaxkb+1)**2, ngw, gx, g, ainv, dylm )
         !
         do is = 1, nsp
            if( iprsta .ge. 4 ) WRITE( stdout,*)  '  dbeta  '
            c = fpi / sqrt(omega)
            do iv = 1, nh(is)
               lp = nhtolm(iv,is)
               betagl = betagx(1,iv,is)
               do i=1,3
                  do j=1,3
                     dbeta(1,iv,is,i,j)=-0.5*beta(1,iv,is)*ainv(j,i)    &
     &                                 +c*dylm(1,lp,i,j)*betagl
                  enddo
               enddo
               do ig=gstart,ngw
                  betagl = betagx(ig,iv,is)
                  dbetagl= dbetagx(ig,iv,is)
                  do i=1,3
                     do j=1,3
                        dbeta(ig,iv,is,i,j)=                            &
     &                    -0.5*beta(ig,iv,is)*ainv(j,i)                 &
     &                    +c*dylm(ig,lp,i,j)*betagl                     &
     &                    -c*ylm (ig,lp)*dbetagl*gx(i,ig)/g(ig)         &
     &                    *(gx(1,ig)*ainv(j,1)+                         &
     &                      gx(2,ig)*ainv(j,2)+                         &
     &                      gx(3,ig)*ainv(j,3))
                     end do
                  end do
               end do
            end do
         end do
         !
         deallocate(dylm)
         !
      end if
      !
      deallocate(ylm)
      IF( ALLOCATED( betagx  ) ) DEALLOCATE( betagx )
      IF( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )

      RETURN
    END SUBROUTINE exact_beta

!  ----------------------------------------------
   END MODULE pseudopotential
!  ----------------------------------------------
