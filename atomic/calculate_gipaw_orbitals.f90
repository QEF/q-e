!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
SUBROUTINE calculate_gipaw_orbitals
  !------------------------------------------------------------
  
  use ld1_parameters, ONLY : nwfsx
  use ld1inc,         only : nld, file_logder, zed, grid, &
       wfc_ae_recon, wfc_ps_recon, nspin, nwftsc, lltsc, enltsc, &
       nwfs, eltsc, el, nwf, enl, rel, vpot, nbeta, pseudotype, &
       vpstot, nstoaets, enlts, vnl, lls, betas, ddd, qq, jjs, els, &
       rcutus, ikk, nwfts, verbosity
  USE io_global,     ONLY : stdout
  use kinds,         only : dp
  use radial_grids,  only : ndmx, series
  implicit none
  
  !<apsi> from lderiv.f90
  integer ::        &
       lam,    &   ! the angular momentum
       nc,     &   ! counter on logarithmic derivatives
       idum,   &   ! integer variable for lschps
       is,     &   ! counter on spin
       ierr,   &   ! used for allocation control
       ios,    &   ! used for I/O control
       n,ie        ! generic counter 
  
  real(DP) ::            &
       aux(ndmx),         & ! the square of the wavefunction
       aux_dir(ndmx,2),   & ! the square of the wavefunction
       ze2,              & ! the nuclear charge in Ry units
       e,                & ! the eigenvalue
       j                   ! total angular momentum for log_der
  
  integer  ::       &
       nbf           ! number of b functions

  real(DP) ::  &
       jam,     &    ! the total angular momentum
       lamsq,            & ! combined angular momentum
       b(0:3),c(4),      & ! used for starting guess of the solution 
       b0e, rr1,rr2,     & ! auxiliary
       xl1, x4l6, ddx12, &
       x6l12, x8l20

  real(DP) :: &
       vaux(ndmx),     &  ! auxiliary: the potential 
       al(ndmx)           ! the known part of the differential equation

  real(DP), external ::           &
       int_0_inf_dr

  integer :: &
       ib,jb,iib,jjb, &  ! counters on beta functions
       nst,nstop,     &  ! auxiliary for integrals
       ind               ! counters on index
  
  REAL ( dp ) :: r1, r2
  
  INTEGER :: ns, outer_r_radial_index, n_idx, ik, i, r_write_max_index
  INTEGER :: ik_stop
  REAL ( dp ) :: factor, rcut_match, max_val
  REAL ( dp ) :: wfc_ae_recon2(ndmx,nwfsx)
  REAL ( dp ) :: wfc_ps_recon2(ndmx,nwfsx)
  REAL ( dp ) :: de
  REAL ( dp ), ALLOCATABLE :: f_ae(:), f_ps(:)
  
  IF ( verbosity == 'high' ) THEN
     WRITE ( stdout, '( 3A )' ) &
          "     ------------------------", &
          " (GI)PAW reconstruction ", &
          "-------------------------"
  END IF
  if ( nld > nwfsx ) &
       call errore ( 'calculate_gipaw_orbitals', 'nld is too large', 1 )
  
  ze2=-zed*2.0_dp
  
  ! Choose a radius for the normalisation of all-electron wave functions
  ! In principle the value for radius is arbitrary [apsi]
  ik_stop = 10
  DO n = 1, grid%mesh
     IF ( grid%r(n) < 1.5 ) ik_stop = n
  END DO
  
  do is=1,nspin
     
     do ns = 1, nwftsc(1)
        lam = lltsc(ns,1)
        j = 0.0_dp
        
        DO n = 1, grid%mesh
           IF ( grid%r(n) < 15.0 ) THEN
              outer_r_radial_index = n
           END IF
        END DO
        
        n_idx = -1
        IF ( ABS ( enltsc(ns,1) ) > 1e-8 ) THEN
           e = enltsc(ns,1)
           DO n = 1, nwf
              IF ( eltsc(ns,1) == el(n) ) THEN
                 n_idx = n
              END IF
           END DO
        ELSE
           DO n = 1, nwf
              IF ( eltsc(ns,1) == el(n) ) THEN
                 e = enl(n)
                 n_idx = n
              END IF
           END DO
        END IF
        
        IF ( verbosity == 'high' ) THEN
           WRITE ( stdout, '( /, 5X, 3A )') "========================= ", &
                TRIM(el(n_idx)), &
                " ========================="
           WRITE ( stdout, * ) "     AE: e(ref), l, n_idx ", e, lam, n_idx
        END IF
        
           !
           !    integrate outward up to ik_stop
           !
           if (rel == 1) then
              call lschps(3,zed,grid,idum,outer_r_radial_index,1,lam, &
                   e,aux,vpot(1,is),nstop)
           else if (rel == 2) then
              call dir_outward(ndmx,outer_r_radial_index,lam,j,e,grid%dx,&
                   aux_dir,grid%r,grid%rab,vpot(1,is))
              aux(:)=aux_dir(:,1)
           else
              call intref(lam,e,outer_r_radial_index,grid,vpot(1,is),ze2,aux)
           endif
           
           wfc_ae_recon(:grid%mesh,n_idx) = aux(:grid%mesh)
           
           ! Set the maximum to be +-1
           wfc_ae_recon(:grid%mesh,n_idx) = wfc_ae_recon(:grid%mesh,n_idx) &
                / MAXVAL ( ABS ( wfc_ae_recon(:ik_stop,n_idx) ) )
           
     end do
     
  enddo
  !</apsi> from lderiv.f90
  
  !<apsi> from lderivps.f90
  
  do is = 1, nspin
     
     IF ( .NOT. rel < 2 ) THEN
        CALL errore ( 'calculate_gipaw_orbitals', &
             'not implemented for rel >= 2', rel )
     END IF
     
     do ns = 1, nwftsc(1)
        lam = lltsc(ns,1)
        jam=0.0_dp
        
        xl1=lam+1
        x4l6=4*lam+6
        x6l12=6*lam+12
        x8l20=8*lam+20
        ddx12=grid%dx*grid%dx/12.0_dp
        nst=(lam+1)**2  
        nbf=nbeta
        if (pseudotype == 1) then
           if (rel < 2 .or. lam == 0 .or. abs(jam-lam+0.5_dp) < 0.001_dp) then
              ind=1
           else if (rel==2 .and. lam>0 .and. abs(jam-lam-0.5_dp)<0.001_dp) then
              ind=2
           endif
           do n=1,grid%mesh
              vaux(n) = vpstot(n,is) + vnl(n,lam,ind)
           enddo
           nbf=0.0
        else
           do n=1,grid%mesh
              vaux(n) = vpstot(n,is)
           enddo
        endif
        
        do n=1,4
           al(n)=vaux(n)-ze2/grid%r(n)
        enddo
        call series(al,grid%r,grid%r2,b)
        
        IF ( ABS ( enltsc(ns,1) ) > 1e-8 ) THEN
           e = enltsc(ns,1)
        ELSE
           e = enlts(ns)
        END IF
        
        DO n = 1, nwftsc(1)
           IF ( eltsc(n,1) == el(nstoaets(ns)) ) THEN
              n_idx = n
           END IF
        END DO
        IF ( verbosity == 'high' ) THEN
           WRITE ( stdout, '( /, 5X, 3A )' ) "========================= ", &
                TRIM(eltsc(ns,1)), &
                " ========================="
           WRITE ( stdout, * ) "     PS: e(ref), e(eig), n_idx ", &
                e, enlts(ns), n_idx
        END IF
        
        rcut_match = -1.0_dp
        DO n = 1, nwfs
           ! If this one has already a core radius...
           IF ( els(n) == el(nstoaets(ns)) ) THEN
              rcut_match = rcutus(n)
           END IF
        END DO
        IF ( rcut_match < 0.0_dp ) THEN
           DO n = 1, nwfs
              ! If there is one with the same l...
              IF ( lltsc(ns,1) == lls(n) ) THEN
                 rcut_match = rcutus(n)
              END IF
           END DO
        END IF
        IF ( rcut_match < 0.0_dp ) THEN
           max_val = -1.0_dp
           DO n = 1, grid%mesh
              IF ( grid%r(n) < 5.0_dp &
                   .AND. ABS ( wfc_ae_recon(n,n_idx) ) > max_val ) THEN
                 rcut_match = grid%r(n)
                 max_val = ABS ( wfc_ae_recon(n,n_idx) )
              END IF
           END DO
           ! In case over divergence
           rcut_match = MIN ( rcut_match, 5.0_dp )
        END IF
        ik = 10
        DO n = 1, grid%mesh
           IF ( grid%r(n) < rcut_match ) ik = n
        END DO
        IF ( verbosity == 'high' ) THEN
           WRITE ( stdout, * ) "     r(cut): ", rcut_match, ik, outer_r_radial_index
        END IF
        
        DO n = 1, grid%mesh
           IF ( grid%r(n) < 15.0 ) THEN
              outer_r_radial_index = n
           END IF
        END DO
        
           lamsq=(lam+0.5_dp)**2
           !
           !     b) find the value of solution s in the first two points
           !
           b0e=b(0)-e
           c(1)=0.5_dp*ze2/xl1
           c(2)=(c(1)*ze2+b0e)/x4l6
           c(3)=(c(2)*ze2+c(1)*b0e+b(1))/x6l12
           c(4)=(c(3)*ze2+c(2)*b0e+c(1)*b(1)+b(2))/x8l20
           r1 = grid%r(1)
           r2 = grid%r(2)
           rr1=(1.0_dp+r1*(c(1)+r1*(c(2)+r1*(c(3)+r1*c(4)))))*r1**(lam+1)
           rr2=(1.0_dp+r2*(c(1)+r2*(c(2)+r2*(c(3)+r2*c(4)))))*r2**(lam+1)
           aux(1)=rr1/grid%sqr(1)
           aux(2)=rr2/grid%sqr(2)
           
           do n=1,grid%mesh
              al(n)=( (vaux(n)-e)*grid%r2(n) + lamsq )*ddx12
              al(n)=1.0_dp-al(n)
           enddo
           
           call integrate_outward (lam,jam,e,grid%mesh,ndmx,grid,al,b,aux, &
                betas,ddd,qq,nbf,nwfsx,lls,jjs,ikk,outer_r_radial_index)
           wfc_ps_recon(:grid%mesh,n_idx) = aux(:grid%mesh) &
                * SQRT ( grid%r(:grid%mesh) )
           
           DO n = 1, grid%mesh
              IF ( grid%r(n) > 5.0 ) EXIT
           END DO
           
           IF ( ABS ( wfc_ps_recon(ik,n_idx) ) < 1e-5 ) THEN
              WRITE ( stdout, * ) "     Warning: ", wfc_ps_recon(ik,n_idx), ns
              CALL errore ( "calculate_gipaw_orbitals", &
                   "safer to stop here...", ik )
           END IF
           
           factor = wfc_ae_recon(ik,nstoaets(n_idx)) / wfc_ps_recon(ik,n_idx)
           wfc_ps_recon(:grid%mesh,n_idx) = wfc_ps_recon(:grid%mesh,n_idx) &
                * factor
           
           IF ( verbosity == 'high' ) THEN
              WRITE ( stdout, * ) "     SCALE: ", &
                   wfc_ae_recon(ik,nstoaets(n_idx)) &
                   / wfc_ps_recon(ik,n_idx), grid%r(ik), factor
              WRITE ( stdout, * ) "     SCALE: ", &
                   wfc_ae_recon(ik+5,nstoaets(n_idx)) &
                   / wfc_ps_recon(ik+5,n_idx),grid%r(ik+5)
              
              ALLOCATE ( f_ae(grid%mesh), f_ps(grid%mesh) )
              
              f_ae = wfc_ae_recon(:grid%mesh,nstoaets(n_idx)) ** 2
              f_ps = wfc_ps_recon(:grid%mesh,n_idx) ** 2
              
              ! Test the norm
              nst = ( lam + 1 ) ** 2
              IF ( verbosity == 'high' ) THEN
                 WRITE ( stdout, '(A,3F12.8)' ) "     NORM: ", &
                      int_0_inf_dr ( f_ae, grid, ik, nst ), &
                      int_0_inf_dr ( f_ps, grid, ik, nst )
              END IF
              
              DEALLOCATE ( f_ae, f_ps )
           END IF 
           
        end do
     
  end do
  !</apsi> from lderivps.f90
  
  IF ( verbosity == 'high' ) THEN
     WRITE ( stdout, '( 3A )' ) &
          "     ---------------------", &
          " End of (GI)PAW reconstruction ", &
          "---------------------"
  END IF
  
END SUBROUTINE calculate_gipaw_orbitals
