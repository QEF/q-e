!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"


   SUBROUTINE compute_dvan_x()
     !     
     !     calculate array  dvan(iv,jv,is)
     !    
     !  rw**2 * vrps   = [ ( Vpsnl(r) - Vpsloc(r) )* Rps(r) * r^2 ]
     !                 = [ DVpsnl(r) * Rps(r) * r^2 ]
     !  dion           = (2l+1) / < Rps(r) | DVpsnl(r) | Rps(r) >

     USE kinds,      ONLY: DP
     use uspp,       only: dvan, nhtolm, indv
     use uspp_param, only: upf, nhm, nh
     use ions_base,  only: nsp
     use atom,       only: numeric
     !
     implicit none
     !
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
             dvan( iv, jv, is ) = fac * upf(is)%dion( indv(iv,is), indv(jv,is) )
           endif
         end do
       end do
     end do
     RETURN
   END SUBROUTINE compute_dvan_x



!------------------------------------------------------------------------------!



   SUBROUTINE pseudopotential_indexes_x( nlcc_any )

      use parameters, only: lmaxx    !
      use ions_base,  only: nsp, &   !  number of specie
                            na       !  number of atoms for each specie
      use cvan,       only: ish      !
      use uspp,       only: nkb, &   !
                            nkbus    !
      use uspp_param, only: upf,  &!
                            lmaxkb, &!
                            nhm,    &!
                            nbetam, &!
                            nh,     &!
                            lmaxq    !
      use uspp,       only: nhtol,  &!
                            nhtolm, &!
                            indv     !
      use atom,       only: nlcc     !

      use pseudopotential,         ONLY: nsanl
      USE read_pseudo_module_fpmd, ONLY: nspnl

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
         do iv = 1, upf(is)%nbeta
            lmaxkb = max( lmaxkb, upf(is)%lll( iv ) )
            ind = ind + 2 * upf(is)%lll( iv ) + 1
         end do
         nh(is) = ind
         ish(is)=nkb
         nkb = nkb + na(is) * nh(is)
         if(  upf(is)%tvanp ) nkbus = nkbus + na(is) * nh(is)
         nlcc_any = nlcc_any .OR. nlcc(is)
      end do
      nhm    = MAXVAL( nh(1:nsp) )
      nbetam = MAXVAL( upf(1:nsp)%nbeta )
      if (lmaxkb > lmaxx) call errore(' pseudopotential_indexes ',' l > lmax ',lmaxkb)
      lmaxq = 2*lmaxkb + 1
      !
      ! the following prevents an out-of-bound error: nqlc(is)=2*lmax+1
      ! but in some versions of the PP files lmax is not set to the maximum
      ! l of the beta functions but includes the l of the local potential
      !
      do is=1,nsp
          upf(is)%nqlc = MIN (  upf(is)%nqlc, lmaxq )
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
         do iv = 1,  upf(is)%nbeta
            lm =  upf(is)%lll(iv)**2
            do il = 1, 2* upf(is)%lll( iv ) + 1
               lm = lm + 1
               ind = ind + 1
               nhtolm( ind, is ) = lm
               nhtol( ind, is ) =  upf(is)%lll( iv )
               indv( ind, is ) = iv
            end do
         end do
      end do

      ! ...     Calculate the number of atoms with non local pseudopotentials
      !
      nsanl = SUM( na(1:nspnl) )

      RETURN
   END SUBROUTINE pseudopotential_indexes_x


!------------------------------------------------------------------------------!


   LOGICAL FUNCTION chkpstab_x(hg, xgtabmax)
      !
      USE kinds,         ONLY: DP
      USE mp,            ONLY: mp_max
      USE io_global,     ONLY: stdout
      USE mp_global,     ONLY: intra_image_comm
      USE cell_base,     ONLY: tpiba
      USE control_flags, ONLY: iprsta
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: hg(:)
      REAL(DP), INTENT(IN) :: xgtabmax
      REAL(DP) :: xgmax

      chkpstab_x = .FALSE.
      !
      xgmax = tpiba * SQRT( MAXVAL( hg ) )
      CALL mp_max( xgmax, intra_image_comm )
      !
      IF( xgmax > xgtabmax ) THEN
         chkpstab_x = .TRUE.
         IF( iprsta > 2 ) &
            WRITE( stdout, fmt='(  "CHKPSTAB: recalculate pseudopotential table" )' )
      END IF
      !
      RETURN
   END FUNCTION chkpstab_x


!------------------------------------------------------------------------------!


   SUBROUTINE compute_xgtab_x( xgmin, xgmax, xgtabmax )
      !
      USE kinds,              ONLY : DP
      USE cell_base,          ONLY : tpiba, tpiba2
      USE mp,                 ONLY : mp_max
      USE mp_global,          ONLY : intra_image_comm
      USE reciprocal_vectors, ONLY : g
      USE pseudopotential,    ONLY : xgtab
      USE betax,              ONLY : mmx
      !
      IMPLICIT NONE
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
      CALL mp_max(xgmax, intra_image_comm)
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
   END SUBROUTINE compute_xgtab_x


!------------------------------------------------------------------------------!


   SUBROUTINE build_pstab_x( )

      USE kinds,              ONLY : DP
      USE atom,               ONLY : rgrid, numeric
      USE ions_base,          ONLY : nsp, rcmax, zv
      USE cell_base,          ONLY : tpiba, tpiba2
      use bhs,                ONLY : rc1, rc2, wrc2, wrc1, rcl, al, bl, lloc
      USE splines,            ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,        ONLY : formfn, formfa
      USE uspp_param,         only : upf, oldvan
      USE control_flags,      only : tpre
      use reciprocal_vectors, ONLY : g, gstart
      USE cp_interfaces,      ONLY : compute_xgtab, chkpstab
      USE pseudopotential,    ONLY : vps_sp, dvps_sp, xgtab
      USE betax,              ONLY : mmx

      IMPLICIT NONE

      INTEGER   :: is, ig
      REAL(DP)  :: xgmax, xgmin
      LOGICAL   :: compute_tab
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

            call formfn( vps_sp(is)%y, dvps_sp(is)%y, rgrid(is)%r, &
                         rgrid(is)%rab, upf(is)%vloc(1:rgrid(is)%mesh), zv(is),   &
                         rcmax(is), xgtab, 1.0d0, tpiba2, rgrid(is)%mesh, &
                         mmx, oldvan(is), tpre )

         else

            !     bhs pseudopotentials
            !
            call formfa( vps_sp(is)%y, dvps_sp(is)%y, rc1(is), rc2(is), &
                         wrc1(is), wrc2(is), rcl(:,is,lloc(is)), &
                         al(:,is,lloc(is)), bl(:,is,lloc(is)), zv(is), &
                         rcmax(is), xgtab, 1.0d0, tpiba2, mmx, 2 , tpre )

         end if

         ! WRITE( 13, "(3D16.8)" ) ( xgtab(ig), vps_sp(is)%y(ig), dvps_sp(is)%y(ig), ig = 1, mmx )

         CALL init_spline( vps_sp(is) )
         CALL init_spline( dvps_sp(is) )

      END DO

      RETURN
   END SUBROUTINE  build_pstab_x


!------------------------------------------------------------------------------!


   SUBROUTINE build_cctab_x( )

      USE kinds,              ONLY : DP
      USE atom,               ONLY : rgrid, nlcc, rho_atc
      USE ions_base,          ONLY : nsp, rcmax
      USE cell_base,          ONLY : tpiba, tpiba2
      USE splines,            ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,        ONLY : compute_rhocg
      USE reciprocal_vectors, ONLY : g, gstart
      USE cp_interfaces,      ONLY : compute_xgtab, chkpstab
      USE pseudopotential,    ONLY : rhoc1_sp, rhocp_sp, xgtab
      USE betax,              ONLY : mmx

      IMPLICIT NONE

      INTEGER  :: is
      REAL(DP) :: xgmax, xgmin
      LOGICAL  :: compute_tab
      REAL(DP) :: xgtabmax = 0.0d0
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
            CALL compute_rhocg( rhoc1_sp(is)%y, rhocp_sp(is)%y, rgrid(is)%r, &
                 rgrid(is)%rab, rho_atc(:,is), xgtab, 1.0d0, tpiba2, rgrid(is)%mesh, mmx, 1 )
            !
            CALL init_spline( rhoc1_sp(is) )
            CALL init_spline( rhocp_sp(is) )
            !
         END IF

      END DO

      RETURN
   END SUBROUTINE  build_cctab_x


!------------------------------------------------------------------------------!


   SUBROUTINE build_nltab_x( )

      ! ...     Initialize Tables for array WNL

      USE kinds,                   ONLY : DP
      USE ions_base,               ONLY : nsp, rcmax
      USE cell_base,               ONLY : tpiba, tpiba2
      USE splines,                 ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,             ONLY : nlin_base
      USE pseudo_base,             ONLY : nlin_stress_base
      USE reciprocal_vectors,      ONLY : g, gstart
      USE uspp_param,              ONLY : upf, nbetam
      USE read_pseudo_module_fpmd, ONLY : nspnl
      USE cp_interfaces,           ONLY : compute_xgtab, chkpstab
      USE pseudopotential,         ONLY : wnl_sp, wnla_sp, xgtab
      USE betax,                   ONLY : mmx
      

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

           DO l = 1,  upf(is)%nbeta
              CALL allocate_spline(  wnl_sp(l,is), mmx, xgmin, xgmax )
              CALL allocate_spline( wnla_sp(l,is), mmx, xgmin, xgmax )
           END DO

           ALLOCATE( fintl( SIZE( xgtab ), SIZE( wnl_sp, 1) ) )
           !
           CALL nlin_base( upf(is), xgtab(:), fintl)
           !
           DO l = 1,  upf(is)%nbeta
              wnl_sp( l, is )%y = fintl(:,l)
           END DO
           !
           CALL nlin_stress_base( upf(is), xgtab, fintl )
           
           DO l = 1,  upf(is)%nbeta
               wnla_sp( l, is )%y = fintl(:,l)
           END DO
           !
           DO l = 1,  upf(is)%nbeta
              CALL init_spline(  wnl_sp( l, is ) )
              CALL init_spline( wnla_sp( l, is ) )
           END DO

           DEALLOCATE(fintl)

        END DO

        RETURN
   END SUBROUTINE  build_nltab_x


!------------------------------------------------------------------------------!


   SUBROUTINE nlin_x( wsg, wnl )

      !  this routine computes the temporary arrays twnl
      !  to be used by nlrh and dforce
      !

      USE kinds,                   ONLY : DP
      USE ions_base,               ONLY: nsp
      USE cell_base,               ONLY: omega, tpiba
      USE pseudo_base,             ONLY: nlin_base
      USE control_flags,           ONLY: gamma_only
      USE uspp,                    ONLY: dvan
      USE uspp_param,              ONLY: upf, nh
      USE constants,               ONLY: pi
      USE splines,                 ONLY: spline
      USE read_pseudo_module_fpmd, ONLY: nspnl
      USE reciprocal_vectors,      ONLY: g, gstart
      USE cp_interfaces,           ONLY: build_nltab
      USE pseudopotential,         ONLY: wnl_sp, tpstab

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
      ik  = 1

      DO is = 1, nspnl
            !
            IF( tpstab ) THEN
               !
               DO l = 1, upf(is)%nbeta
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
               CALL nlin_base( upf(is), g, wnl(:,:,is,ik) )
               !
            END IF
            !
      END DO

      RETURN
   END SUBROUTINE nlin_x


!------------------------------------------------------------------------------!


   SUBROUTINE nlin_stress_x( wnla )

      !  this routine computes the temporary arrays wnla
      !  to be used by stress subroutine.
      !
      !  Note: subroutine nlin should be called first
      !  
      USE kinds,                   ONLY : DP
      USE ions_base,               ONLY : nsp
      USE cell_base,               ONLY : tpiba
      USE pseudo_base,             ONLY : nlin_stress_base
      USE splines,                 ONLY : spline
      USE uspp_param,              ONLY : upf
      USE read_pseudo_module_fpmd, ONLY : nspnl
      USE reciprocal_vectors,      ONLY : g, gstart
      USE pseudopotential,         ONLY : wnla_sp, tpstab

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
            DO l = 1, upf(is)%nbeta
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
            CALL nlin_stress_base( upf(is), g, wnla(:,:,is))
            !
         END IF
         !
      END DO
      !
      RETURN
   END SUBROUTINE nlin_stress_x


!------------------------------------------------------------------------------!


   SUBROUTINE compute_betagx_x( tpre )
      !
      ! calculation of array  betagx(ig,iv,is)
      !
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nsp
      USE uspp_param, ONLY : upf, nh, nhm, oldvan
      USE atom,       ONLY : rgrid, numeric
      USE uspp,       ONLY : nhtol, indv
      USE betax,      only : refg, betagx, mmx, dbetagx
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ir, nr
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:)
      REAL(DP) :: xg
      !
      IF( ALLOCATED( betagx  ) ) DEALLOCATE( betagx )
      IF( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
      !
      ALLOCATE( betagx ( mmx, nhm, nsp ) )
      IF ( tpre ) ALLOCATE( dbetagx( mmx, nhm, nsp ) )
      !
      do is = 1, nsp
         !
         nr = upf(is)%kkbeta
         !
         if ( tpre ) then
            allocate( dfint( nr ) )
            allocate( djl  ( nr ) )
         end if
         !
         allocate( fint ( nr ) )
         allocate( jl   ( nr ) )
         !
         do iv = 1, nh(is)
            !
            l = nhtol(iv,is)
            !
            do il = 1, mmx
               !
               xg = sqrt( refg * (il-1) )
               call sph_bes ( nr, rgrid(is)%r, xg, l, jl )
!
               if( tpre )then
                  !
                  call sph_dbes1 ( nr, rgrid(is)%r, xg, l, jl, djl)
                  !
               endif
               !
               !     beta(ir)=r*beta(r)
               !
               do ir = 1, nr
                  fint(ir) = rgrid(is)%r(ir) * jl(ir) * &
                             upf(is)%beta( ir, indv(iv,is) ) 
               end do
               if (oldvan(is)) then
                  call herman_skillman_int(nr,fint,rgrid(is)%rab,betagx(il,iv,is))
               else
                  call simpson_cp90(nr,fint,rgrid(is)%rab,betagx(il,iv,is))
               endif
               ! 
               if(tpre) then
                  do ir = 1, nr
                     dfint(ir) = rgrid(is)%r(ir) * djl(ir) * &
                                 upf(is)%beta( ir, indv(iv,is) )
                  end do
                  if (oldvan(is)) then
                     call herman_skillman_int(nr,dfint,rgrid(is)%rab,dbetagx(il,iv,is))
                  else
                     call simpson_cp90(nr,dfint,rgrid(is)%rab,dbetagx(il,iv,is))
                  end if
               endif
               !
            end do
         end do
!
         deallocate(jl)
         deallocate(fint)
         !
         if (tpre) then
            deallocate(djl)
            deallocate(dfint)
         end if
         !
      end do
      RETURN
   END SUBROUTINE compute_betagx_x


!------------------------------------------------------------------------------!


   SUBROUTINE compute_qradx_x( tpre )
      !
      !     calculation of array qradx(igb,iv,jv,is) for interpolation table
      !     (symmetric wrt exchange of iv and jv: a single index ijv is used)
      !
      !       qradx(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
      !     
      !
      !
      USE kinds,         ONLY : DP
      use io_global,     only : stdout
      USE ions_base,     ONLY : nsp
      USE uspp_param,    ONLY : upf, nh, nhm, nbetam, lmaxq, oldvan
      USE atom,          ONLY : rgrid, numeric
      USE uspp,          ONLY : indv
      USE betax,         only : refg, qradx, mmx, dqradx
      USE cvan,          only : ish, nvb
      use gvecb,         only : ngb
      USE cp_interfaces, ONLY : fill_qrl
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ir, jv, ijv, ierr
      INTEGER :: nr
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), qrl(:,:,:)
      REAL(DP) :: xg

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
      !
      ALLOCATE(  qradx( mmx, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      !
      IF ( tpre ) ALLOCATE( dqradx( mmx, nbetam*(nbetam+1)/2, lmaxq, nsp ) )

      DO is = 1, nvb
         !
         !     qqq and beta are now indexed and taken in the same order
         !     as vanderbilts ppot-code prints them out
         !
         WRITE( stdout,*) ' nlinit  nh(is), ngb, is, kkbeta, lmaxq = ', &
     &        nh(is), ngb, is, upf(is)%kkbeta, upf(is)%nqlc
         !
         nr = upf(is)%kkbeta
         !
         IF ( tpre ) THEN
            ALLOCATE( djl  ( nr ) )
            ALLOCATE( dfint( nr ) )
         END IF
         !
         ALLOCATE( fint( nr ) )
         ALLOCATE( jl  ( nr ) )
         ALLOCATE( qrl( nr, upf(is)%nbeta*(upf(is)%nbeta+1)/2, upf(is)%nqlc) )
         !
         call fill_qrl ( is, qrl )
         !
         do l = 1, upf(is)%nqlc
            !
            do il = 1, mmx
               !
               xg = sqrt( refg * DBLE(il-1) )
               !
               call sph_bes ( nr, rgrid(is)%r, xg, l-1, jl(1) )
               !
               if( tpre ) then
                  !
                  call sph_dbes1 ( nr, rgrid(is)%r, xg, l-1, jl, djl)
                  !
               endif
               !
               ! 
               do iv = 1, upf(is)%nbeta
                  do jv = iv, upf(is)%nbeta
                     ijv = jv * ( jv - 1 ) / 2 + iv
                     !
                     !      note qrl(r)=r^2*q(r)
                     !
                     do ir = 1, nr
                        fint( ir ) = qrl( ir, ijv, l ) * jl( ir )
                     end do
                     if (oldvan(is)) then
                        call herman_skillman_int &
                             (nr,fint(1),rgrid(is)%rab,qradx(il,ijv,l,is))
                     else
                        call simpson_cp90 &
                             (nr,fint(1),rgrid(is)%rab,qradx(il,ijv,l,is))
                     end if
                     !
                     if( tpre ) then
                        do ir = 1, nr
                           dfint(ir) = qrl(ir,ijv,l) * djl(ir)
                        end do
                        if ( oldvan(is) ) then
                           call herman_skillman_int &
                                (nr,dfint(1),rgrid(is)%rab,dqradx(il,ijv,l,is))
                        else
                           call simpson_cp90 &
                                (nr,dfint(1),rgrid(is)%rab,dqradx(il,ijv,l,is))
                        end if
                     end if
                     !
                  end do
               end do
               !
               !
            end do
         end do
         !
         DEALLOCATE (  jl )
         DEALLOCATE ( qrl )
         DEALLOCATE ( fint  )
         !
         if ( tpre ) then
            DEALLOCATE(djl)
            DEALLOCATE ( dfint )
         end if
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    qqq '
         !
         do iv=1,upf(is)%nbeta
            WRITE( stdout,'(8f9.4)') (upf(is)%qqq(iv,jv),jv=1,upf(is)%nbeta)
         end do
         WRITE( stdout,*)
         !
      end do

      RETURN
    END SUBROUTINE compute_qradx_x

!------------------------------------------------------------------------------!

    SUBROUTINE exact_qradb_x( tpre )
      !
      USE kinds,         ONLY : DP
      use io_global,  only: stdout
      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: upf, nh, nhm, nbetam, lmaxq, oldvan
      use uspp_param, only: lmaxkb
      USE atom,       ONLY: rgrid, numeric
      USE uspp,       ONLY: indv
      use uspp,       only: qq, beta
      USE betax,      only: refg, qradx, mmx, dqradx
      USE cvan,       only: ish, nvb
      use gvecb,      only: ngb
      use control_flags, only: iprint, iprsta
      use cell_base,  only: ainv
      use constants,  only: pi, fpi
      use qradb_mod,  only: qradb
      use qgb_mod,    only: qgb
      use gvecb,      only: gb, gxb
      use small_box,  only: omegab, tpibab
      use dqrad_mod,  only: dqrad
      use dqgb_mod,   only: dqgb
      USE cp_interfaces, ONLY: fill_qrl
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ir, jv, ijv, ierr
      INTEGER :: ig, i,j, jj, nr
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), qrl(:,:,:)
      REAL(DP) :: xg, c, betagl, dbetagl, gg
      REAL(DP), ALLOCATABLE :: dqradb(:,:,:,:)
      REAL(DP), ALLOCATABLE :: ylmb(:,:), dylmb(:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: dqgbs(:,:,:)

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
      !
      ALLOCATE(  qradx( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      !
      IF ( tpre ) ALLOCATE( dqradx( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )

      DO is = 1, nvb
         !
         !     qqq and beta are now indexed and taken in the same order
         !     as vanderbilts ppot-code prints them out
         !
         WRITE( stdout,*) ' nlinit  nh(is), ngb, is, kkbeta, lmaxq = ', &
     &        nh(is), ngb, is, upf(is)%kkbeta, upf(is)%nqlc
         !
         nr = upf(is)%kkbeta
         !
         IF ( tpre ) THEN
            ALLOCATE( djl  ( nr ) )
            ALLOCATE( dfint( nr ) )
         END IF
         !
         ALLOCATE( fint( nr ) )
         ALLOCATE( jl  ( nr ) )
         ALLOCATE( qrl( nr, upf(is)%nbeta*(upf(is)%nbeta+1)/2, upf(is)%nqlc) )
         !
         call fill_qrl ( is, qrl )
         ! qrl = 0.0d0
         !
         do l = 1, upf(is)%nqlc
            !
            do il = 1, ngb
               !
               xg = sqrt( gb( il ) * tpibab * tpibab )
               !
               call sph_bes ( nr, rgrid(is)%r, xg, l-1, jl(1) )
               !
               if( tpre ) then
                  !
                  call sph_dbes1 ( nr, rgrid(is)%r, xg, l-1, jl, djl)
                  !
               endif
               !
               ! 
               do iv = 1, upf(is)%nbeta
                  do jv = iv, upf(is)%nbeta
                     ijv = jv * ( jv - 1 ) / 2 + iv
                     !
                     !      note qrl(r)=r^2*q(r)
                     !
                     do ir = 1, nr
                        fint( ir ) = qrl( ir, ijv, l ) * jl( ir )
                     end do
                     if (oldvan(is)) then
                        call herman_skillman_int &
                             (nr,fint(1),rgrid(is)%rab,qradx(il,ijv,l,is))
                     else
                        call simpson_cp90 &
                             (nr,fint(1),rgrid(is)%rab,qradx(il,ijv,l,is))
                     end if
                     !
                     if( tpre ) then
                        do ir = 1, nr
                           dfint(ir) = qrl(ir,ijv,l) * djl(ir)
                        end do
                        if ( oldvan(is) ) then
                           call herman_skillman_int &
                                (nr,dfint(1),rgrid(is)%rab,dqradx(il,ijv,l,is))
                        else
                           call simpson_cp90 &
                                (nr,dfint(1),rgrid(is)%rab,dqradx(il,ijv,l,is))
                        end if
                     end if
                     !
                  end do
               end do
               !
               !
            end do
         end do
         !
         DEALLOCATE (  jl )
         DEALLOCATE ( qrl )
         DEALLOCATE ( fint  )
         !
         if ( tpre ) then
            DEALLOCATE(djl)
            DEALLOCATE ( dfint )
         end if
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    qqq '
         !
         do iv=1, upf(is)%nbeta
            WRITE( stdout,'(8f9.4)') (upf(is)%qqq(iv,jv),jv=1, upf(is)%nbeta)
         end do
         WRITE( stdout,*)
         !
      end do

      allocate( ylmb( ngb, lmaxq*lmaxq ), STAT=ierr )
      IF( ierr  /= 0 ) &
        CALL errore(' exact_qradb ', ' cannot allocate ylmb ', 1 )
!
      qradb(:,:,:,:) = 0.d0
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

      do is = 1, nvb
         !
         !     calculation of array qradb(igb,iv,jv,is)
         !
         if( iprsta .ge. 4 ) WRITE( stdout,*)  '  qradb  '
         !
         c = fpi / omegab
         !
         do iv= 1, upf(is)%nbeta
            do jv = iv, upf(is)%nbeta
               ijv = jv*(jv-1)/2 + iv
               do ig=1,ngb
                  do l=1,upf(is)%nqlc
                     qradb(ig,ijv,l,is)= c*qradx(ig,ijv,l,is)
                  enddo
               enddo
            enddo
         enddo
         !
         !     ---------------------------------------------------------------
         !     stocking of qgb(igb,ijv,is) and of qq(iv,jv,is)
         !     ---------------------------------------------------------------
         !
         do iv= 1,nh(is)
            do jv=iv,nh(is)
               !
               !       compact indices because qgb is symmetric
               !
               ijv = jv*(jv-1)/2 + iv
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
         dqrad(:,:,:,:,:,:) = 0.d0
         !
         call dylmr2_(lmaxq*lmaxq, ngb, gxb, gb, ainv, dylmb)
         !
         do is=1,nvb
            !
            do iv= 1, upf(is)%nbeta
               do jv=iv, upf(is)%nbeta
                  ijv = jv*(jv-1)/2 + iv
                  do l=1,upf(is)%nqlc
                     do ig=1,ngb
                        dqradb(ig,ijv,l,is) =  dqradx(ig,ijv,l,is)
                     enddo
                     do i=1,3
                        do j=1,3
                           dqrad(1,ijv,l,is,i,j) = &
                                -qradb(1,ijv,l,is) * ainv(j,i)
                           do ig=2,ngb
                              dqrad(ig,ijv,l,is,i,j) =                  &
     &                          -qradb(ig,ijv,l,is)*ainv(j,i)           &
     &                          -c*dqradb(ig,ijv,l,is)*                 &
     &                          gxb(i,ig)/gb(ig)*                       &
     &                          (gxb(1,ig)*ainv(j,1)+                   &
     &                           gxb(2,ig)*ainv(j,2)+                   &
     &                           gxb(3,ig)*ainv(j,3))
                           enddo
                        enddo
                     enddo
                  end do
               enddo
            enddo
            !
            do iv= 1,nh(is)
               do jv=iv,nh(is)
                  !
                  !       compact indices because qgb is symmetric
                  !
                  ijv = jv*(jv-1)/2 + iv
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

      deallocate( ylmb )

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )

      RETURN
    END SUBROUTINE exact_qradb_x


!------------------------------------------------------------------------------!


    LOGICAL FUNCTION check_tables_x( )
      !
      ! check table size against cell variations
      !
      !
      USE kinds,              ONLY : DP
      USE betax,              ONLY : refg, mmx
      USE mp,                 ONLY : mp_max
      USE mp_global,          ONLY : intra_image_comm
      USE gvecw,              ONLY : ngw
      USE cell_base,          ONLY : tpiba2
      USE small_box,          ONLY : tpibab
      USE gvecb,              ONLY : gb, ngb
      USE reciprocal_vectors, ONLY : g
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
      CALL mp_max( gmax, intra_image_comm )
      !
      check_tables_x = .FALSE.
      IF( ( INT( gmax ) + 2 ) > mmx ) check_tables_x = .TRUE.
      !
      RETURN
    END FUNCTION check_tables_x
    

!------------------------------------------------------------------------------!


    SUBROUTINE interpolate_beta_x( tpre )
      !
      ! interpolate array beta(ig,iv,is)
      !
      !
      USE kinds, ONLY : DP
      USE control_flags, only: iprsta
      USE constants, only: pi, fpi
      USE io_global, only: stdout
      USE gvecw, only: ngw
      USE ions_base, only: nsp
      USE reciprocal_vectors, only: g, gx, gstart
      USE uspp_param, only: upf, lmaxq, lmaxkb, nh
      USE uspp, only: qq, nhtolm, beta
      USE cell_base, only: ainv, omega, tpiba2, tpiba
      USE betax, ONLY : refg, betagx, dbetagx
      USE cdvan, ONLY : dbeta

      IMPLICIT NONE

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
            do ig = gstart, ngw
               gg = g( ig ) * tpiba * tpiba / refg
               jj = int( gg ) + 1
               betagl = betagx( jj+1, iv, is ) * ( gg - DBLE(jj-1) ) + betagx( jj, iv, is ) * ( DBLE(jj) - gg )
               beta( ig, iv, is ) = c * ylm( ig, lp ) * betagl
            end do
            if( gstart == 2 ) then
               beta( 1, iv, is ) = c * ylm( 1, lp ) * betagx( 1, iv, is )
            end if
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
                     dbeta( 1, iv, is, i, j ) = -0.5d0 * beta( 1, iv, is ) * ainv( j, i )    &
     &                                          - c * dylm( 1, lp, i, j ) * betagl         ! SEGNO
                  enddo
               enddo
               do ig = gstart, ngw
                  gg = g(ig) * tpiba * tpiba / refg
                  jj=int(gg)+1
                  betagl = betagx(  jj+1, iv, is ) * ( gg - DBLE(jj-1) ) +         &
     &                     betagx(  jj  , iv, is ) * ( DBLE(jj) - gg )
                  dbetagl= dbetagx( jj+1, iv, is ) * ( gg - DBLE(jj-1) ) +        &
     &                     dbetagx( jj  , iv, is ) * ( DBLE(jj) - gg )
                  do i=1,3
                     do j=1,3
                        dbeta( ig, iv, is, i, j ) =                            &
     &                    - 0.5d0 * beta( ig, iv, is ) * ainv( j, i )          &
     &                    - c * dylm( ig, lp, i, j ) * betagl                  &  ! SEGNO
     &                    - c * ylm ( ig, lp )       * dbetagl * gx( i, ig ) / g( ig )         &
     &                    * ( gx( 1, ig ) * ainv( j, 1 ) + gx( 2, ig ) * ainv( j, 2 ) + gx( 3, ig ) * ainv( j, 3 ) )
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
    END SUBROUTINE interpolate_beta_x


!------------------------------------------------------------------------------!


   SUBROUTINE interpolate_qradb_x( tpre )
      !
      ! interpolate array qradb(ig,iv,is)
      !
      !
      USE kinds,         ONLY : DP
      use control_flags, only: iprint, iprsta
      use io_global, only: stdout
      use gvecw, only: ngw
      use cell_base, only: ainv
      use cvan, only: nvb
      use uspp, only: qq, nhtolm, beta
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use uspp_param, only: upf, lmaxq, lmaxkb, nbetam, nh
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

      integer  is, l, ig, ir, iv, jv, ijv, i,j, jj, ierr
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
      qradb(:,:,:,:) = 0.d0
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

      do is = 1, nvb
         !
         !     calculation of array qradb(igb,iv,jv,is)
         !
         if( iprsta .ge. 4 ) WRITE( stdout,*)  '  qradb  '
         !
         c = fpi / omegab
         !
         do iv= 1, upf(is)%nbeta
            do jv = iv, upf(is)%nbeta
               ijv = jv*(jv-1)/2 + iv
               do l=1, upf(is)%nqlc
                  qradb(1,ijv,l,is) = c * qradx(1,ijv,l,is)
               end do
               do ig=2,ngb
                  gg=gb(ig)*tpibab*tpibab/refg
                  jj=int(gg)+1
                  do l=1,upf(is)%nqlc
                     if(jj.ge.mmx) then
                        qradb(ig,ijv,l,is)=0.d0
                     else
                        qradb(ig,ijv,l,is)=                           &
     &                       c*qradx(jj+1,ijv,l,is)*(gg-DBLE(jj-1))+  &
     &                       c*qradx(jj,ijv,l,is)*(DBLE(jj)-gg)
                     endif
                  enddo
               enddo
            enddo
         enddo
!
!     ---------------------------------------------------------------
!     stocking of qgb(igb,ijv,is) and of qq(iv,jv,is)
!     ---------------------------------------------------------------
         do iv= 1,nh(is)
            do jv=iv,nh(is)
!
!       compact indices because qgb is symmetric
!
               ijv = jv*(jv-1)/2 + iv
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
         dqrad(:,:,:,:,:,:) = 0.d0
         !
         call dylmr2_( lmaxq*lmaxq, ngb, gxb, gb, ainv, dylmb )
         !
         do is=1,nvb
            !
            do iv= 1, upf(is)%nbeta
               do jv=iv, upf(is)%nbeta
                  ijv = jv*(jv-1)/2 + iv
                  do l=1,upf(is)%nqlc
                     dqradb(1,ijv,l,is) =  dqradx(1,ijv,l,is)
                     do ig=2,ngb
                        gg=gb(ig)*tpibab*tpibab/refg
                        jj=int(gg)+1
                        if(jj.ge.mmx) then
                           dqradb(ig,ijv,l,is) = 0.d0
                        else
                           dqradb(ig,ijv,l,is) =  &
                                dqradx(jj+1,ijv,l,is)*(gg-DBLE(jj-1)) +  &
                                dqradx(jj,ijv,l,is)*(DBLE(jj)-gg)
                        endif
                     enddo
                     do i=1,3
                        do j=1,3
                           dqrad(1,ijv,l,is,i,j) = - qradb(1,ijv,l,is) * ainv(j,i)
                           do ig=2,ngb
                              dqrad(ig,ijv,l,is,i,j) =                  &
     &                          - qradb(ig,ijv,l,is)*ainv(j,i)          &  
     &                          - c * dqradb(ig,ijv,l,is)*              &
     &                          gxb(i,ig)/gb(ig)*                       &
     &                          (gxb(1,ig)*ainv(j,1)+                   &
     &                           gxb(2,ig)*ainv(j,2)+                   &
     &                           gxb(3,ig)*ainv(j,3))
                           enddo
                        enddo
                     enddo
                  end do
               enddo
            enddo
            !
            do iv= 1,nh(is)
               do jv=iv,nh(is)
                  !
                  !       compact indices because qgb is symmetric
                  !
                  ijv = jv*(jv-1)/2 + iv
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
    END SUBROUTINE interpolate_qradb_x


!------------------------------------------------------------------------------!


    SUBROUTINE exact_beta_x( tpre )
      !
      ! compute array beta without interpolation
      !
      !
      USE control_flags, only : iprsta
      USE kinds,         ONLY : DP
      USE constants,     only : pi, fpi
      USE io_global,     only : stdout
      USE gvecw,         only : ngw
      USE ions_base,     only : nsp
      USE uspp_param,    only : upf, lmaxq, lmaxkb, nh, nhm, oldvan
      USE uspp,          only : qq, nhtolm, beta, nhtol, indv
      USE cell_base,     only : ainv, omega, tpiba2, tpiba
      USE cdvan,         ONLY : dbeta
      USE atom,          ONLY : rgrid, numeric
      USE reciprocal_vectors, only : g, gx, gstart

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: tpre
 
      REAL(DP), ALLOCATABLE ::  ylm(:,:), dylm(:,:,:,:)
      REAL(DP) :: c, gg, betagl, dbetagl
      INTEGER :: is, iv, lp, ig, jj, i, j, nr
      INTEGER :: l, il, ir
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:)
      REAL(DP), ALLOCATABLE :: betagx ( :, :, : ), dbetagx( :, :, : )
      REAL(DP) :: xg

      ALLOCATE( ylm( ngw, (lmaxkb+1)**2 ) )
      ALLOCATE( betagx ( ngw, nhm, nsp ) )
      IF (tpre) ALLOCATE( dbetagx( ngw, nhm, nsp ) )

      CALL ylmr2 ( (lmaxkb+1)**2, ngw, gx, g, ylm)

      !
      do is = 1, nsp
         !
         nr = upf(is)%kkbeta
         !
         if ( tpre ) then
            allocate( dfint( nr ) )
            allocate( djl  ( nr ) )
         end if
         !
         allocate( fint ( nr ) )
         allocate( jl   ( nr ) )
         !
         do iv = 1, nh(is)
            !
            l = nhtol(iv,is)
            !
            do il = 1, ngw
               !
               xg = sqrt( g( il ) * tpiba * tpiba )
               call sph_bes (nr, rgrid(is)%r, xg, l, jl )
               !
               if( tpre )then
                  !
                  call sph_dbes1 ( nr, rgrid(is)%r, xg, l, jl, djl)
                  !
               endif
               !
               !     beta(ir)=r*beta(r)
               !
               do ir = 1, nr
                  fint(ir) = rgrid(is)%r(ir) * jl(ir) * &
                             upf(is)%beta( ir, indv(iv,is) )
               end do
               if (oldvan(is)) then
                  call herman_skillman_int(nr,fint,rgrid(is)%rab,betagx(il,iv,is))
               else
                  call simpson_cp90(nr,fint,rgrid(is)%rab,betagx(il,iv,is))
               endif
               ! 
               if(tpre) then
                  do ir = 1, nr
                     dfint(ir) = rgrid(is)%r(ir) * djl(ir) * &
                                 upf(is)%beta( ir, indv(iv,is) )
                  end do
                  if (oldvan(is)) then
                     call herman_skillman_int(nr,dfint,rgrid(ir)%rab,dbetagx(il,iv,is))
                  else
                     call simpson_cp90(nr,dfint,rgrid(is)%rab,dbetagx(il,iv,is))
                  end if
               endif
               !
            end do
         end do
!
         deallocate(jl)
         deallocate(fint)
         !
         if (tpre) then
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
                     dbeta(1,iv,is,i,j)=-0.5d0*beta(1,iv,is)*ainv(j,i)    &
     &                                 -c*dylm(1,lp,i,j)*betagl  ! SEGNO
                  enddo
               enddo
               do ig=gstart,ngw
                  betagl = betagx(ig,iv,is)
                  dbetagl= dbetagx(ig,iv,is)
                  do i=1,3
                     do j=1,3
                        dbeta(ig,iv,is,i,j)=                            &
     &                    -0.5d0*beta(ig,iv,is)*ainv(j,i)                 &
     &                    -c*dylm(ig,lp,i,j)*betagl                     &  ! SEGNO
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
    END SUBROUTINE exact_beta_x
!
!
!------------------------------------------------------------------------------!
!
!
   SUBROUTINE fill_qrl_x( is, qrl )
      !
      ! fill l-components of Q(r) as in Vanderbilt's approach
      !
      USE uspp_param, ONLY: upf
      USE atom,       ONLY: rgrid
      USE kinds,      ONLY: DP
      USE io_global,  ONLY: stdout
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN)  :: is
      REAL(DP), INTENT(OUT) :: qrl( :, :, : )
      !
      INTEGER :: iv, jv, ijv, lmin, lmax, l, ir, i
      INTEGER :: dim1, dim2, dim3
      !
      dim1 = SIZE( qrl, 1 )
      dim2 = SIZE( qrl, 2 )
      dim3 = SIZE( qrl, 3 )
      !
      IF ( upf(is)%kkbeta > dim1 ) &
           CALL errore ('fill_qrl', 'bad 1st dimension for array qrl', 1)
      !
      qrl = 0.0d0
      !
      do iv = 1,  upf(is)%nbeta
         !
         do jv = iv,  upf(is)%nbeta
            !
            ijv = (jv-1)*jv/2 + iv
            !
            IF ( ijv > dim2) &
                 CALL errore ('fill_qrl', 'bad 2nd dimension for array qrl', 2)

            ! notice that L runs from 1 to Lmax+1

            lmin = ABS (upf(is)%lll(jv) - upf(is)%lll(iv)) + 1
            lmax = upf(is)%lll(jv) + upf(is)%lll(iv) + 1

            ! WRITE( stdout, * ) 'QRL is, jv, iv = ', is, jv, iv
            ! WRITE( stdout, * ) 'QRL lll jv, iv = ', upf(is)%lll(jv), upf(is)%lll(iv)
            ! WRITE( stdout, * ) 'QRL lmin, lmax = ', lmin, lmax
            ! WRITE( stdout, * ) '---------------- '

            IF ( lmin < 1 .OR. lmax > dim3) THEN
                 WRITE( stdout, * ) ' lmin, lmax = ', lmin, lmax
                 CALL errore ('fill_qrl', 'bad 3rd dimension for array qrl', 3)
            END IF


             do l = lmin, lmax
               do ir = 1, upf(is)%kkbeta
                  if ( rgrid(is)%r(ir) >= upf(is)%rinner(l) ) then
                     qrl(ir,ijv,l)=upf(is)%qfunc(ir,ijv)
                  else
                     qrl(ir,ijv,l)=upf(is)%qfcoef(1,l,iv,jv)
                     do i = 2, upf(is)%nqf
                        qrl(ir,ijv,l)=qrl(ir,ijv,l) +      &
                             upf(is)%qfcoef(i,l,iv,jv)*rgrid(is)%r(ir)**(2*i-2)
                     end do
                     qrl(ir,ijv,l) = qrl(ir,ijv,l) * rgrid(is)%r(ir)**(l+1)
                  end if
               end do
            end do
         end do
      end do
      RETURN
   END SUBROUTINE fill_qrl_x
