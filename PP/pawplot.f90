!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE paw_postproc_

    USE kinds,          ONLY : DP
    USE paw_variables,  ONLY : paw_info
    IMPLICIT NONE

    PUBLIC :: PAW_make_ae_charge_
    PRIVATE

    CONTAINS

SUBROUTINE PAW_make_ae_charge_ ( rho, flag, nx, r, rhopaw )
   USE paw_onecenter,     ONLY : paw_rho_lm
   USE atom,              ONLY : g => rgrid
   USE constants,         ONLY : fpi
   USE ions_base,         ONLY : nat, ityp, tau
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nh, nhm, upf
   USE scf,               ONLY : scf_type
   USE fft_base,          ONLY : dfftp
   USE mp_global,         ONLY : me_pool
   USE splinelib,         ONLY : spline, splint
   USE gvect,             ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
   USE cell_base,         ONLY : at, bg, alat
   !
   TYPE(scf_type), INTENT(IN) :: rho           ! only rho%bec is actually needed
   INTEGER, INTENT (IN)       :: flag          ! -1=core, 0 =valence, 1=both
   INTEGER, INTENT (IN)       :: nx            ! number of points in space
   REAL (dp), INTENT(IN)      :: r(3,nx)       ! points in space (alat units)
   REAL (dp), INTENT(OUT)  :: rhopaw(nx,nspin) ! PAW charge
   !
   TYPE(paw_info)          :: i                ! minimal info on atoms
   INTEGER                 :: ip               ! counter on x,y,z
   INTEGER                 :: ir               ! counter on grid point
   INTEGER                 :: is               ! spin index
   INTEGER                 :: lm               ! counters on angmom and radial grid
   INTEGER                 :: j,k,l, idx, idx0
   INTEGER                 :: ia
   REAL(DP),ALLOCATABLE    :: wsp_lm(:,:,:), ylm_posi(:,:), d1y(:), d2y(:)
   REAL(DP),ALLOCATABLE    :: rho_lm(:,:,:), rho_lm_ae(:,:,:), rho_lm_ps(:,:,:)
   REAL(DP)                :: posi(3), first, second, distsq

   !
   rhopaw (:,:) = 0.0_dp
   !
   ! Currently assuming parallelization on input data points
   ! (no parallelization on atoms)
   !
   atoms: DO ia = 1, nat
      !
      i%a = ia                      ! atom's index
      i%t = ityp(ia)                ! type of atom ia
      i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
      i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
      i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia
      !
      ifpaw: IF (upf(i%t)%tpawp) THEN
         !
         ! Arrays are allocated inside the cycle to allow reduced
         ! memory usage as different atoms have different meshes
         ALLOCATE(rho_lm_ae(i%m,i%l**2,nspin), &
                  rho_lm_ps(i%m,i%l**2,nspin) )
         ALLOCATE(rho_lm(i%m,i%l**2,nspin), &
                  ylm_posi(1,i%l**2),       &
                  wsp_lm(i%m, i%l**2,nspin)  )
         !
         ! Compute rho spherical harmonics expansion from becsum and pfunc
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
         CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, &
              upf(i%t)%qfuncl)
         !
         DO is=1,nspin
            IF ( flag >= 0 ) THEN
               DO lm = 1,i%l**2
                  DO ir = 1, i%m
                     rho_lm(ir,lm,is) = ( rho_lm_ae(ir,lm,is) - &
                          rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir)
                  ENDDO
               ENDDO
            ELSE
               rho_lm(:,:,is) = 0.0_dp
            END IF
            !
            ! add core charge (divide by Y_00=1/sqrt(4pi) to get l=0 component)
            ! 
            IF ( abs(flag) == 1 ) THEN
               DO ir = 1, i%m
                  rho_lm(ir,1,is) = rho_lm(ir,1,is) + &
                       sqrt( fpi ) * upf(i%t)%paw%ae_rho_atc(ir) / nspin
               ENDDO
            END IF
         ENDDO

         ! deallocate asap
         DEALLOCATE(rho_lm_ae, rho_lm_ps)
         !
         ALLOCATE( d1y(upf(i%t)%kkbeta), d2y(upf(i%t )%kkbeta) )
         DO is = 1,nspin
            DO lm = 1, i%l**2
               CALL radial_gradient(rho_lm(1:upf(i%t)%kkbeta,lm,is), d1y, &
                                    g(i%t)%r, upf(i%t)%kkbeta, 1)
               CALL radial_gradient(d1y, d2y, g(i%t)%r, upf(i%t)%kkbeta, 1)
               !
               first  = d1y(1) ! first derivative in first point
               second = d2y(1) ! second derivative in first point
               ! prepare interpolation
               CALL spline( g(i%t)%r(:), rho_lm(:,lm,is), first, second, &
                    wsp_lm(:,lm,is) )
            ENDDO
         ENDDO
         DEALLOCATE(d1y, d2y)
         !
         rsp_point : DO ir = 1, nx
            !
            posi(:) = r(:,ir)
            !
            ! find the distance of real-space grid's point ir w.r.t
            ! closer periodic image of atom ia
            !
            posi(:) = posi(:) - tau(:,ia)
            CALL cryst_to_cart( 1, posi, bg, -1 )
            posi(:) = posi(:) - ANINT( posi(:) )
            CALL cryst_to_cart( 1, posi, at, 1 )
            !
            posi(:) = posi(:) * alat
            distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
            !
            ! don't consider points too far from the atom
            ! (criterion not valid in principle if core charge is present)
            !
            IF ( abs(flag) == 1 .AND. &
                 distsq > g(i%t)%r2(upf(i%t)%kkbeta) ) CYCLE rsp_point
            !
            ! generate the atomic charge on point posi(:), which means
            ! sum over l and m components rho_lm_ae-rho_lm_ps
            ! interpolate the radial function at distance |posi(:)|
            !
            ! prepare spherical harmonics
            CALL ylmr2( i%l**2, 1, posi, distsq, ylm_posi )
            DO is = 1,nspin
               DO lm = 1, i%l**2
                  ! do interpolation
                  rhopaw(ir,is)= rhopaw(ir,is) + ylm_posi(1,lm) &
                       * splint(g(i%t)%r(:) , rho_lm(:,lm,is), &
                       wsp_lm(:,lm,is), SQRT(distsq) )
               ENDDO
            ENDDO
         END DO rsp_point
         !
         DEALLOCATE(rho_lm, ylm_posi, wsp_lm)
         !
      ENDIF ifpaw
   ENDDO atoms
   
 END SUBROUTINE PAW_make_ae_charge_
 
END MODULE paw_postproc_
!
!-----------------------------------------------------------------------
PROGRAM PAWplot
  !-----------------------------------------------------------------------
  !
  !    Plot PAW charge density
  !
  USE kinds,      ONLY : dp
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE mp_global,  ONLY : mp_startup
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start
  USE lsda_mod,   ONLY : nspin, current_spin
  USE cell_base,  ONLY : bg
  USE gvect,      ONLY : ngm, nr1,nr2,nr3, nrx1,nrx2,nrx3, nrxx, nl
  USE scf,        ONLY : rho
  USE io_files,   ONLY : tmp_dir, prefix, trimcheck
  USE noncollin_module, ONLY : noncolin
  USE paw_variables,    ONLY : okpaw
  USE paw_postproc_,    ONLY : PAW_make_ae_charge_
  !
  IMPLICIT NONE
  CHARACTER(len=256) :: outdir, filplot
  CHARACTER(len=16)  :: plot
  INTEGER :: spin_component, nx,ny,nz, flag, ios, is
  REAL(dp) :: e1(3), e2(3), e3(3), x0(3)
  REAL(dp), ALLOCATABLE :: rhoplot(:), rhopaw(:,:), r(:,:)
  COMPLEX(dp), ALLOCATABLE :: rhog(:)
  LOGICAL, external :: matches
  LOGICAL :: onedim, twodim, tredim
  !
  NAMELIST / inputpp / outdir, prefix, spin_component, &
       filplot, e1, e2, e3, x0, nx, ny, nz, plot
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'POST-PROC' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  filplot = 'pawcharge.dat' 
  plot    = 'valence'
  spin_component = 0
  e1(:)         = 0.d0
  e2(:)         = 0.d0
  e3(:)         = 0.d0
  x0(:)         = 0.d0
  nx            = 0
  ny            = 0
  nz            = 0
  !
  ios = 0
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  END IF
  call mp_bcast (ios, ionode_id)
  IF ( ios /= 0) &
       CALL errore ('postproc', 'reading inputpp namelist', ABS(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( e1, ionode_id )
  CALL mp_bcast( e2, ionode_id )
  CALL mp_bcast( e3, ionode_id )
  CALL mp_bcast( x0, ionode_id )
  CALL mp_bcast( nx, ionode_id )
  CALL mp_bcast( ny, ionode_id )
  CALL mp_bcast( nz, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( plot, ionode_id )
  CALL mp_bcast( spin_component, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  ALLOCATE ( rhog(ngm) )
  !
  !      plot of the charge density - select rho(G)
  !
  IF (noncolin) THEN
     rhog (:) = rho%of_g(:,1)
  ELSE
     IF (spin_component == 0) THEN
        rhog (:) = rho%of_g(:,1)
        DO is = 2, nspin
           rhog(:) = rhog (:) + rho%of_g(:,is)
        ENDDO
     ELSE
        IF (nspin == 2) current_spin = spin_component
        rhog (:) = rho%of_g(:,current_spin)
     ENDIF
  ENDIF
  !  
  tredim = ( e3(1)**2 + e3(2)**2 + e3(3)**2 > 1d-6 )
  twodim = ( e2(1)**2 + e2(2)**2 + e2(3)**2 > 1d-6 ) .AND. .NOT. tredim
  onedim = ( e1(1)**2 + e1(2)**2 + e1(3)**2 > 1d-6 ) .AND. .NOT. twodim
  !
  if ( onedim ) then
     !
     !     One-dimensional plot
     !
     if (nx <= 0 )   call errore ('chdens', 'wrong nx', 1)
     ALLOCATE ( rhoplot(nx) )
     IF ( okpaw ) THEN
        WRITE (stdout, '(5x,"Reconstructing all-electron charge (PAW)")')
        ALLOCATE ( rhopaw(nx,nspin), r(3,nx) )
        DO is = 1, nx
           r(:, is) = x0 (:) + (is-1) * e1(:) / (nx-1)
        END DO
        !
        IF ( matches ('core',plot) .AND. matches ('valence',plot) ) THEN
            flag = 1
        ELSE IF ( matches ('core',plot) ) THEN
            flag =-1
        ELSE
            flag = 0 
        END IF 
        CALL PAW_make_ae_charge_ (rho, flag, nx, r, rhopaw )
        !
        IF (spin_component == 0 .and. nspin ==2 ) THEN
           rhoplot(:) = rhopaw(:,1)+ rhopaw(:,2)
        ELSE
           IF (nspin == 2) current_spin = spin_component
           rhoplot(:) = rhopaw(:,current_spin)
        END IF
        DEALLOCATE ( r, rhopaw )
     ELSE
        rhoplot(:) = 0.0_dp
     END IF
     !
     call plot_1d_ (nx, x0, e1, rhog, rhoplot, flag, filplot )
     !
     DEALLOCATE ( rhoplot )
     !
  else if ( twodim ) then
     if ( abs(e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)) > 1d-6) &
          call errore ('pawplot', 'e1 and e2 are not orthogonal', 1)
     if ( nx <= 0 .OR. ny <= 0 )   call errore ('chdens', 'wrong nx or ny', 1)
  else if (tredim) then
     if ( nx <= 0 .OR. ny <= 0 .OR. nz <=0 ) &
          call errore ('chdens', 'wrong nx or ny or nz', 1)
  END IF
  !
  DEALLOCATE (rhog)
  
end PROGRAM PAWPLOT
!
!-----------------------------------------------------------------------
subroutine plot_1d_ (nx, x0, e, rhog, rhoplot, flag, filplot )
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  use constants, only:  pi
  USE io_global, only : stdout, ionode
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE gvect,      ONLY : g, gstart, ngm
  USE control_flags, ONLY : gamma_only

  implicit none
  integer, intent(in) :: nx, flag
  ! number of points along the line
  ! flag=-1: do not add smooth term
  real(DP), intent(in) :: e (3), x0 (3)
  ! vector defining the line
  ! origin of the line
  complex(DP), intent(in) :: rhog (ngm)
  ! rho in G space
  character(LEN=*), intent(in) :: filplot
  real(DP), intent(inout) :: rhoplot(nx)
  !
  integer :: i, ig, ounit
  real(DP) :: rhosum(nx), rhomin, rhomax, x(3), deltax, arg
  !
  do i = 1, nx
     x(:) = x0 (:) + (i-1) * e (:) / (nx-1)
     !
     !     for each point we compute the charge from the Fourier components
     !
     rhosum(i) = 0.0_dp
     do ig = gstart, ngm
        !
        !     NB: G are in 2pi/alat units, r are in alat units
        !
        arg = 2.0_dp*pi* ( x(1)*g(1,ig) + x(2)*g(2,ig) + x(3)*g(3,ig) )
        rhosum(i) = rhosum(i) + DBLE ( rhog (ig) ) * cos (arg) - &
                               AIMAG ( rhog (ig) ) * sin (arg)
     enddo
     if ( gamma_only )  rhosum(i) = 2.0_dp * rhosum(i)
     if ( gstart == 2 ) rhosum(i) = rhosum(i) + DBLE( rhog (1) )
     !
  enddo
  call mp_sum( rhosum, intra_pool_comm )
  !
  IF ( flag /= -1) rhoplot (:) = rhoplot(:) + rhosum(:)
  !
  !    Here we check the value of the resulting charge
  !
  rhomin =  MINVAL ( rhoplot(:) )
  rhomax =  MAXVAL ( rhoplot(:) )
  
  write(stdout, '(5x,"Min, Max charge: ",2f12.6)') rhomin, rhomax
  !
  !       we print the charge on output
  !
  if (ionode) then
     if (filplot /= ' ') then
        ounit = 1
        open (unit=ounit, file=filplot, form='formatted', status='unknown')
        WRITE( stdout, '(/5x,"Writing data to be plotted to file ",a)') &
             TRIM(filplot)
     else
        ounit = 6
     endif
     !
     deltax = sqrt(e(1)**2+e(2)**2+e(3)**2) / (nx - 1)
     do i = 1, nx
        write (ounit, '(2f20.10)') deltax*DBLE(i-1), rhoplot(i)
     enddo
     if (ounit == 1) CLOSE (unit = ounit, status='keep')
  end if
  
  return

end subroutine plot_1d_

