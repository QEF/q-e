!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!



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
       !     fac converts ry to hartree
       fac = 0.5d0
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



   SUBROUTINE pseudopotential_indexes_x( )

      use parameters, only: lmaxx    !
      use ions_base,  only: nsp, &   !  number of specie
                            na       !  number of atoms for each specie
      use uspp,       only: nkb, &   !
                            nkbus    !
      use uspp_param, only: ish,    &!
                            upf,    &!
                            lmaxkb, &!
                            nhm,    &!
                            nbetam, &!
                            nh,     &!
                            lmaxq    !
      use uspp,       only: nhtol,  &!
                            nhtolm, &!
                            indv     !

      IMPLICIT NONE
     
      !
      INTEGER :: is, iv, ind, il, lm
      !     ------------------------------------------------------------------
      !     find  number of beta functions per species, max dimensions,
      !     total number of beta functions (all and Vanderbilt only)
      !     ------------------------------------------------------------------
      lmaxkb   = -1
      nkb      =  0
      nkbus    =  0
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

      RETURN
   END SUBROUTINE pseudopotential_indexes_x


!------------------------------------------------------------------------------!


   SUBROUTINE compute_xgtab_x( xgmin, xgmax )
      !
      USE kinds,              ONLY : DP
      USE cell_base,          ONLY : tpiba, tpiba2
      USE mp,                 ONLY : mp_max
      USE mp_global,          ONLY : intra_bgrp_comm
      USE gvect,              ONLY : gg
      USE pseudopotential,    ONLY : xgtab
      USE betax,              ONLY : mmx, refg
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(OUT)  :: xgmax, xgmin
      !
      INTEGER   :: ig
      REAL(DP)  :: xg, dxg, res
      !
      IF( ALLOCATED( xgtab ) )  & 
         DEALLOCATE( xgtab )

      ALLOCATE( xgtab( mmx ) )
      !
      xgmin = 0.0d0
      xgmax = SQRT( refg * mmx )
      dxg   = (xgmax - xgmin) / DBLE( mmx - 1 )
      !
      DO ig = 1, SIZE( xgtab )
         xgtab(ig) = xgmin + DBLE(ig-1) * dxg
      END DO
      !
      xgtab = xgtab**2 / tpiba**2
      !
      RETURN
   END SUBROUTINE compute_xgtab_x


!------------------------------------------------------------------------------!


   SUBROUTINE build_pstab_x( )

      USE kinds,              ONLY : DP
      USE atom,               ONLY : rgrid
      USE ions_base,          ONLY : nsp, rcmax, zv
      USE cell_base,          ONLY : tpiba, tpiba2
      USE splines,            ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,        ONLY : formfn, formfa
      USE uspp_param,         only : upf, oldvan
      USE control_flags,      only : tpre
      use gvect, ONLY : gg, gstart
      USE cp_interfaces,      ONLY : compute_xgtab
      USE pseudopotential,    ONLY : vps_sp, dvps_sp, xgtab
      USE local_pseudo,       ONLY : vps0
      USE betax,              ONLY : mmx

      IMPLICIT NONE

      INTEGER   :: is, ig
      REAL(DP)  :: xgmax, xgmin
      LOGICAL   :: compute_tab
      !
      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' build_pstab_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' build_pstab_x ', ' upf not allocated ', 1 )
      !
      IF( ALLOCATED( vps_sp ) .AND. ALLOCATED( dvps_sp ) ) THEN
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
      IF(  ALLOCATED( vps_sp ) .OR. ALLOCATED( dvps_sp ) ) THEN
         CALL errore( ' build_pstab_x ', ' inconsistent allocation ', 1 )
      END IF
      !
      CALL compute_xgtab( xgmin, xgmax )
      !
      ALLOCATE( vps_sp(nsp))
      ALLOCATE( dvps_sp(nsp))
      !
      DO is = 1, nsp

         CALL nullify_spline( vps_sp( is ) )
         CALL nullify_spline( dvps_sp( is ) )

         CALL allocate_spline( vps_sp(is), mmx, xgmin, xgmax )
         CALL allocate_spline( dvps_sp(is), mmx, xgmin, xgmax )

         call formfn( rgrid(is)%r, rgrid(is)%rab, &
                      upf(is)%vloc(1:rgrid(is)%mesh), zv(is), rcmax(is), &
                      xgtab, 1.0d0, tpiba2, rgrid(is)%mesh, mmx, oldvan(is),&
                      tpre, vps_sp(is)%y, vps0(is), dvps_sp(is)%y )
         ! obsolete BHS format
         !call formfa( vps_sp(is)%y, dvps_sp(is)%y, rc1(is), rc2(is), &
         !             wrc1(is), wrc2(is), rcl(:,is,lloc(is)), &
         !             al(:,is,lloc(is)), bl(:,is,lloc(is)), zv(is), &
         !             rcmax(is), xgtab, 1.0d0, tpiba2, mmx, 2 , tpre )

         ! WRITE( 13, "(3D16.8)" ) ( xgtab(ig), vps_sp(is)%y(ig), dvps_sp(is)%y(ig), ig = 1, mmx )

         CALL init_spline( vps_sp(is) )
         CALL init_spline( dvps_sp(is) )

      END DO

      RETURN
   END SUBROUTINE  build_pstab_x


!------------------------------------------------------------------------------!


   SUBROUTINE build_cctab_x( )

      USE kinds,              ONLY : DP
      USE atom,               ONLY : rgrid
      USE uspp_param,         ONLY : upf
      USE ions_base,          ONLY : nsp, rcmax
      USE cell_base,          ONLY : tpiba, tpiba2
      USE splines,            ONLY : init_spline, allocate_spline, kill_spline, nullify_spline
      USE pseudo_base,        ONLY : compute_rhocg
      USE gvect, ONLY : gg, gstart
      USE cp_interfaces,      ONLY : compute_xgtab
      USE pseudopotential,    ONLY : rhoc1_sp, rhocp_sp, xgtab
      USE betax,              ONLY : mmx

      IMPLICIT NONE

      INTEGER  :: is
      REAL(DP) :: xgmax, xgmin
      LOGICAL  :: compute_tab
      !
      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' build_cctab_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' build_cctab_x ', ' upf not allocated ', 1 )
      !
      IF( ALLOCATED( rhoc1_sp ) .AND. ALLOCATED( rhocp_sp ) ) THEN
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
      IF(  ALLOCATED( rhoc1_sp ) .OR. ALLOCATED( rhocp_sp ) ) THEN
         CALL errore( ' build_cctab_x ', ' inconsistent allocation ', 1 )
      END IF
      !
      CALL compute_xgtab( xgmin, xgmax )
      !
      ALLOCATE( rhoc1_sp(nsp))
      ALLOCATE( rhocp_sp(nsp))
      !
      DO is = 1, nsp

         CALL nullify_spline( rhoc1_sp( is ) )
         CALL nullify_spline( rhocp_sp( is ) )

         IF( upf(is)%nlcc ) THEN
            !
            CALL allocate_spline( rhoc1_sp(is), mmx, xgmin, xgmax )
            CALL allocate_spline( rhocp_sp(is), mmx, xgmin, xgmax )
            !
            CALL compute_rhocg( rhoc1_sp(is)%y, rhocp_sp(is)%y, rgrid(is)%r, &
                 rgrid(is)%rab, upf(is)%rho_atc(:), xgtab, 1.0d0, tpiba2, &
                 rgrid(is)%mesh, mmx, 1 )
            !
            CALL init_spline( rhoc1_sp(is) )
            CALL init_spline( rhocp_sp(is) )
            !
         END IF

      END DO

      RETURN
   END SUBROUTINE  build_cctab_x


!------------------------------------------------------------------------------!


   SUBROUTINE compute_betagx_x( tpre )
      !
      ! calculation of array  betagx(ig,iv,is)
      !
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nsp
      USE uspp_param, ONLY : upf, nh, nhm, oldvan
      USE atom,       ONLY : rgrid
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
      CALL start_clock('betagx')
      !
      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' compute_betagx_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' compute_betagx_x ', ' upf not allocated ', 1 )
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
         do iv = 1, nh(is)
            !
            l = nhtol(iv,is)
            !
!$omp parallel default(none), private( dfint, djl, fint, jl, il, xg, ir ), &
!$omp shared( tpre, nr, mmx, refg, l, is, rgrid, upf, indv, iv, betagx, dbetagx, oldvan )
            if ( tpre ) then
               allocate( dfint( nr ) )
               allocate( djl  ( nr ) )
            end if
            !
            allocate( fint ( nr ) )
            allocate( jl   ( nr ) )
            !
!$omp do
            interp_tab : do il = 1, mmx
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
            end do interp_tab
!$omp end do
            !
            deallocate(jl)
            deallocate(fint)
            !
            if (tpre) then
               deallocate(djl)
               deallocate(dfint)
            end if
            !
!$omp end parallel
         !
         end do
         !
      end do
      CALL stop_clock('betagx')
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
      USE uspp_param,    ONLY : upf, nh, nhm, nbetam, lmaxq, oldvan, ish, nvb
      USE atom,          ONLY : rgrid
      USE uspp,          ONLY : indv
      USE betax,         only : refg, qradx, mmx, dqradx
      use smallbox_gvec,         only : ngb
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

      CALL start_clock('qradx')

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' compute_qradx_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' compute_qradx_x ', ' upf not allocated ', 1 )

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
         ALLOCATE( qrl( nr, upf(is)%nbeta*(upf(is)%nbeta+1)/2, upf(is)%nqlc) )
         !
         call fill_qrl ( is, qrl )
         !
         do l = 1, upf(is)%nqlc
            !
!$omp parallel default(none), private( djl, dfint, fint, jl, il, iv, jv, ijv, xg, ir ), &
!$omp shared( tpre, nr, mmx, refg, rgrid, l, upf, qrl, oldvan, qradx, dqradx, is )
            IF ( tpre ) THEN
               ALLOCATE( djl  ( nr ) )
               ALLOCATE( dfint( nr ) )
            END IF
            !
            ALLOCATE( fint( nr ) )
            ALLOCATE( jl  ( nr ) )
!$omp do
            interp_tab : do il = 1, mmx
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
            end do interp_tab
!$omp end do
            !
            DEALLOCATE (  jl )
            DEALLOCATE ( fint  )
            !
            if ( tpre ) then
               DEALLOCATE(djl)
               DEALLOCATE ( dfint )
            end if
            !
!$omp end parallel
         !
         end do
         !
         DEALLOCATE ( qrl )
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    qqq '
         !
         do iv=1,upf(is)%nbeta
            WRITE( stdout,'(8f9.4)') (upf(is)%qqq(iv,jv),jv=1,upf(is)%nbeta)
         end do
         WRITE( stdout,*)
         !
      end do

      CALL stop_clock('qradx')

      RETURN
    END SUBROUTINE compute_qradx_x

!------------------------------------------------------------------------------!

    SUBROUTINE exact_qradb_x( tpre )
      !
      USE kinds,         ONLY : DP
      use io_global,  only: stdout
      USE ions_base,  ONLY: nsp
      USE uspp_param, ONLY: upf, nh, nhm, nbetam, lmaxq, oldvan
      use uspp_param, only: lmaxkb, ish, nvb
      USE atom,       ONLY: rgrid
      USE uspp,       ONLY: indv
      use uspp,       only: qq, beta
      USE betax,      only: refg, qradx, mmx, dqradx
      use smallbox_gvec,      only: ngb
      use control_flags, only: iprint, iverbosity
      use cell_base,  only: ainv
      use constants,  only: pi, fpi
      use qgb_mod,    only: qgb, dqgb
      use smallbox_gvec,      only: gb, gxb
      use small_box,  only: omegab, tpibab
      USE cp_interfaces, ONLY: fill_qrl
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: tpre
      !
      INTEGER :: is, iv, l, il, ir, jv, ijv, ierr
      INTEGER :: ig, i,j, jj, nr
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:), qrl(:,:,:)
      REAL(DP) :: xg, c, betagl, dbetagl
      REAL(DP), ALLOCATABLE :: dqradb(:,:,:,:)
      REAL(DP), ALLOCATABLE :: dqrad( :, :, :, :, :, : )
      REAL(DP), ALLOCATABLE :: qradb(:,:,:,:)
      REAL(DP), ALLOCATABLE :: ylmb(:,:), dylmb(:,:,:,:)
      COMPLEX(DP), ALLOCATABLE :: dqgbs(:,:,:)

      IF( nvb < 1 ) &
         return

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' exact_qradb_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' exact_qradb_x ', ' upf not allocated ', 1 )

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )
      !
      ALLOCATE(  qradx( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      !
      IF ( tpre ) ALLOCATE( dqradx( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )

      ALLOCATE( qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      qradb(:,:,:,:) = 0.d0


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
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

      do is = 1, nvb
         !
         !     calculation of array qradb(igb,iv,jv,is)
         !
         !     qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
         !
         if( iverbosity > 2 ) WRITE( stdout,*)  '  qradb  '
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
               call qvan2b(ngb,iv,jv,is,ylmb,qgb(1,ijv,is),qradb )
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
         allocate( dqrad( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp, 3, 3 ) )
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
                  call dqvan2b(ngb,iv,jv,is,ylmb,dylmb,dqgbs,dqrad,qradb )
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
         deallocate(dqrad)
         deallocate(dqgbs)
         deallocate(dylmb)
         deallocate(dqradb)
      end if

      deallocate( ylmb )
      deallocate( qradb )

      IF( ALLOCATED(  qradx ) ) DEALLOCATE(  qradx )
      IF( ALLOCATED( dqradx ) ) DEALLOCATE( dqradx )

      RETURN
    END SUBROUTINE exact_qradb_x


!------------------------------------------------------------------------------!


    LOGICAL FUNCTION check_tables_x( gmax )
      !
      ! check table size against cell variations
      !
      !
      USE kinds,              ONLY : DP
      USE betax,              ONLY : refg, mmx
      USE mp,                 ONLY : mp_max
      USE mp_global,          ONLY : intra_bgrp_comm
      USE gvecw,              ONLY : ngw
      USE cell_base,          ONLY : tpiba2
      USE small_box,          ONLY : tpibab
      USE smallbox_gvec,      ONLY : gb, ngb
      USE gvect, ONLY : gg, ngm
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(OUT) :: gmax
      REAL(DP) :: g2, g2b
      !
      g2  = MAXVAL( gg( 1:ngm ) )
      !
      g2  = g2 * tpiba2 / refg
      !
      IF( ALLOCATED( gb ) ) THEN
         !
         g2b = MAXVAL( gb( 1:ngb ) )
         g2b = g2b * tpibab * tpibab / refg
         gmax = MAX( g2, g2b )
         !
      ELSE
         !
         gmax = g2
         !
      END IF
      !
      CALL mp_max( gmax, intra_bgrp_comm )
      !
      check_tables_x = .FALSE.
      IF( INT( gmax ) + 2 >= mmx ) check_tables_x = .TRUE.
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
      USE control_flags, only: iverbosity
      USE constants, only: pi, fpi
      USE io_global, only: stdout
      USE gvecw, only: ngw
      USE ions_base, only: nsp
      USE gvect, only: gg, g, gstart
      USE uspp_param, only: upf, lmaxq, lmaxkb, nh
      USE uspp, only: qq, nhtolm, beta, dbeta
      USE cell_base, only: ainv, omega, tpiba2, tpiba
      USE betax, ONLY : refg, betagx, dbetagx

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: tpre
 
      REAL(DP), ALLOCATABLE ::  ylm(:,:), dylm(:,:,:,:)
      REAL(DP) :: c, g2, betagl, dbetagl
      INTEGER   :: is, iv, lp, ig, jj, i, j

      ALLOCATE( ylm( ngw, (lmaxkb+1)**2 ) )
      CALL ylmr2 ( (lmaxkb+1)**2, ngw, g, gg, ylm)
      !
      !
      do is = 1, nsp
         !   
         !   calculation of array  beta(ig,iv,is)
         !  
         if( iverbosity > 2 ) WRITE( stdout,*)  '  beta  '
         c = fpi / sqrt(omega)
         do iv = 1, nh(is)
            lp = nhtolm( iv, is )
            do ig = gstart, ngw
               g2 = gg( ig ) * tpiba * tpiba / refg
               jj = int( g2 ) + 1
               betagl = betagx( jj+1, iv, is ) * ( g2 - DBLE(jj-1) ) + betagx( jj, iv, is ) * ( DBLE(jj) - g2 )
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
         call dylmr2_( (lmaxkb+1)**2, ngw, g, gg, ainv, dylm )
         !
         do is = 1, nsp
            if( iverbosity > 2 ) WRITE( stdout,*)  '  dbeta  '
            c = fpi / sqrt(omega)
            do iv = 1, nh(is)
               lp = nhtolm(iv,is)
               if( ngw > 0 ) then
                  betagl = betagx(1,iv,is)
                  do i=1,3
                     do j=1,3
                        dbeta( 1, iv, is, i, j ) = -0.5d0 * beta( 1, iv, is ) * ainv( j, i )    &
     &                                             - c * dylm( 1, lp, i, j ) * betagl         ! SEGNO
                     enddo
                  enddo
               end if
               do ig = gstart, ngw
                  g2 = gg(ig) * tpiba * tpiba / refg
                  jj=int(g2)+1
                  betagl = betagx(  jj+1, iv, is ) * ( g2 - DBLE(jj-1) ) +         &
     &                     betagx(  jj  , iv, is ) * ( DBLE(jj) - g2 )
                  dbetagl= dbetagx( jj+1, iv, is ) * ( g2 - DBLE(jj-1) ) +        &
     &                     dbetagx( jj  , iv, is ) * ( DBLE(jj) - g2 )
                  do i=1,3
                     do j=1,3
                        dbeta( ig, iv, is, i, j ) =                            &
     &                    - 0.5d0 * beta( ig, iv, is ) * ainv( j, i )          &
     &                    - c * dylm( ig, lp, i, j ) * betagl                  &  ! SEGNO
     &                    - c * ylm ( ig, lp )       *dbetagl * g(i,ig)/gg(ig)&
     &                    * ( g( 1, ig ) * ainv( j, 1 ) + g( 2, ig ) * ainv( j, 2 ) + g( 3, ig ) * ainv( j, 3 ) )
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
      use control_flags, only: iprint, iverbosity
      use io_global, only: stdout
      use gvecw, only: ngw
      use cell_base, only: ainv
      use uspp, only: qq, nhtolm, beta
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use uspp_param, only: upf, lmaxq, lmaxkb, nbetam, nh, nvb
      use qgb_mod, only: qgb, dqgb
      use smallbox_gvec, only: gb, gxb, ngb
      use small_box,  only: omegab, tpibab
      USE betax, ONLY: qradx, dqradx, refg, mmx
!
      implicit none

      LOGICAL, INTENT(IN) :: tpre

      integer  is, l, ig, ir, iv, jv, ijv, i,j, jj, ierr
      real(dp), allocatable:: fint(:), jl(:), dqradb(:,:,:,:)
      real(dp), allocatable:: ylmb(:,:), dylmb(:,:,:,:)
      REAL(DP), ALLOCATABLE :: dqrad( :, :, :, :, :, : )
      REAL(DP), ALLOCATABLE :: qradb( :, :, :, : )
      complex(dp), allocatable:: dqgbs(:,:,:)
      real(dp) xg, c, betagl, dbetagl, g2
!
      !
      if( nvb < 1 ) &
         return

      allocate( qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ), STAT=ierr )
      IF( ierr  /= 0 ) &
        CALL errore(' interpolate_qradb ', ' cannot allocate qradb ', 1 )
      !
      qradb(:,:,:,:) = 0.d0
!
      allocate( ylmb( ngb, lmaxq*lmaxq ), STAT=ierr )
      IF( ierr  /= 0 ) &
        CALL errore(' interpolate_qradb ', ' cannot allocate ylmb ', 1 )
!
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)

      do is = 1, nvb
         !
         !     calculation of array qradb(igb,iv,jv,is)
         !
         !     qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
         !
         if( iverbosity > 2 ) WRITE( stdout,*)  '  qradb  '
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
                  g2=gb(ig)*tpibab*tpibab/refg
                  jj=int(g2)+1
                  do l=1,upf(is)%nqlc
                     if(jj.ge.mmx) then
                        qradb(ig,ijv,l,is)=0.d0
                     else
                        qradb(ig,ijv,l,is)=                           &
     &                       c*qradx(jj+1,ijv,l,is)*(g2-DBLE(jj-1))+  &
     &                       c*qradx(jj,ijv,l,is)*(DBLE(jj)-g2)
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
               call qvan2b(ngb,iv,jv,is,ylmb,qgb(1,ijv,is),qradb )
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
         allocate( dqrad( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp, 3, 3 ) )
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
                        g2=gb(ig)*tpibab*tpibab/refg
                        jj=int(g2)+1
                        if(jj.ge.mmx) then
                           dqradb(ig,ijv,l,is) = 0.d0
                        else
                           dqradb(ig,ijv,l,is) =  &
                                dqradx(jj+1,ijv,l,is)*(g2-DBLE(jj-1)) +  &
                                dqradx(jj,ijv,l,is)*(DBLE(jj)-g2)
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
                  call dqvan2b(ngb,iv,jv,is,ylmb,dylmb,dqgbs,dqrad,qradb )
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
         deallocate(dqrad)
         deallocate(dqgbs)
         deallocate(dylmb)
         deallocate(dqradb)
      end if
      deallocate(ylmb)
      deallocate(qradb)

      RETURN
    END SUBROUTINE interpolate_qradb_x


!------------------------------------------------------------------------------!


    SUBROUTINE exact_beta_x( tpre )
      !
      ! compute array beta without interpolation
      !
      !
      USE control_flags, only : iverbosity
      USE kinds,         ONLY : DP
      USE constants,     only : pi, fpi
      USE io_global,     only : stdout
      USE gvecw,         only : ngw
      USE ions_base,     only : nsp
      USE uspp_param,    only : upf, lmaxq, lmaxkb, nh, nhm, oldvan
      USE uspp,          only : qq, nhtolm, beta, nhtol, indv, dbeta
      USE cell_base,     only : ainv, omega, tpiba2, tpiba
      USE atom,          ONLY : rgrid
      USE gvect, only : gg, g, gstart

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: tpre
 
      REAL(DP), ALLOCATABLE ::  ylm(:,:), dylm(:,:,:,:)
      REAL(DP) :: c, g2, betagl, dbetagl
      INTEGER :: is, iv, lp, ig, jj, i, j, nr
      INTEGER :: l, il, ir
      REAL(DP), ALLOCATABLE :: dfint(:), djl(:), fint(:), jl(:)
      REAL(DP), ALLOCATABLE :: betagx ( :, :, : ), dbetagx( :, :, : )
      REAL(DP) :: xg

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' exact_beta_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' exact_beta_x ', ' upf not allocated ', 1 )

      ALLOCATE( ylm( ngw, (lmaxkb+1)**2 ) )
      ALLOCATE( betagx ( ngw, nhm, nsp ) )
      IF (tpre) ALLOCATE( dbetagx( ngw, nhm, nsp ) )

      CALL ylmr2 ( (lmaxkb+1)**2, ngw, g, gg, ylm)

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
               xg = sqrt( gg( il ) * tpiba * tpiba )
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
         if( iverbosity > 2 ) WRITE( stdout,*)  '  beta  '
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
         call dylmr2_( (lmaxkb+1)**2, ngw, g, gg, ainv, dylm )
         !
         do is = 1, nsp
            if( iverbosity > 2 ) WRITE( stdout,*)  '  dbeta  '
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
     &                    -0.5d0*beta(ig,iv,is)*ainv(j,i)               &
     &                    -c*dylm(ig,lp,i,j)*betagl                     &  ! SEGNO
     &                    -c*ylm (ig,lp)*dbetagl*g(i,ig)/gg(ig)        &
     &                    *(g(1,ig)*ainv(j,1)+                         &
     &                      g(2,ig)*ainv(j,2)+                         &
     &                      g(3,ig)*ainv(j,3))
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
      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' fill_qrl_x ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' fill_qrl_x ', ' upf not allocated ', 1 )

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
                  IF( upf(is)%tvanp ) THEN
                    ! BEWARE: index l in upf%qfuncl(l) runs from 0 to lmax, not from 1 to lmax+1
                    qrl(ir,ijv,l)=upf(is)%qfuncl(ir,ijv,l-1)
                  ENDIF
               end do
            end do
         end do
      end do
      RETURN
   END SUBROUTINE fill_qrl_x
