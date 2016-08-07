  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
  ! Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  SUBROUTINE analytic_cont_aniso_iaxis_to_raxis( itemp, iter, conv ) 
  !-----------------------------------------------------------------------
  !
  ! This routine does the analytic continuation of the anisotropic Eliashberg equations 
  ! from the imaginary-axis to the real axis
  ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
  !
  ! input
  !
  ! itemp  - temperature point
  ! iter   - iteration number
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE phcom,         ONLY : nmodes
  USE elph2,         ONLY : wqf, wf
  USE epwcom,        ONLY : nqstep, degaussq, nsiter, conv_thr_racon, fsthick, & 
                            lpade, eps_acustic
  USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, wsph, gap, Agap, Gp, Gm, ADsumi, AZsumi, &                           
                            Delta, Znorm, ADelta, ADeltap, AZnorm, AZnormp, g2, lacon_fly, & 
                            a2fij, wkfs, dosef, ixkqf, ixqfs, nqfs, w0g, nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : pi, ci
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER :: i, iw, iwp, iwph, itemp, iter, ik, iq, iq0, ibnd, jbnd, imode, & 
             lower_bnd, upper_bnd, imelt
  REAL(kind=DP) :: rgammap, rgammam, absdelta, reldelta, errdelta, weight, a2f_
  REAL(kind=DP), EXTERNAL :: w0gauss
  COMPLEX(kind=DP) :: esqrt, root
  COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  LOGICAL :: conv
  CHARACTER (len=256) :: cname
  !
  ! SP: Need initialization
  a2f_ = 0.0_DP
  !
  IF ( iter .eq. 1 ) THEN
     !
     ! get the size of required allocated memory for 
     ! Delta, Znorm, Deltaold, ADelta, ADeltap, AZnorm, AZnormp, Gp, Gm 
     IF ( lpade ) THEN                 
        imelt = 2 * ( 1 + 2 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
     ELSE
        imelt = 2 * ( 3 + 4 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
     ENDIF
     CALL mem_size_eliashberg( imelt )
     !
     IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
     IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
     IF ( .not. ALLOCATED(ADelta) )   ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
     IF ( .not. ALLOCATED(ADeltap) )  ALLOCATE( ADeltap(nbndfs,nkfs,nsw) )
     IF ( .not. ALLOCATED(AZnorm) )   ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
     IF ( .not. ALLOCATED(AZnormp) )  ALLOCATE( AZnormp(nbndfs,nkfs,nsw) )
     IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
     ADeltap(:,:,:) = (0.d0, 0.d0)
     AZnormp(:,:,:) = (1.d0, 0.d0)
     Deltaold(:) = (0.d0, 0.d0)
     IF ( lpade ) THEN 
        ADeltap(:,:,:) = ADelta(:,:,:)
        Deltaold(:) = Delta(:)
     ELSE
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 ADeltap(ibnd,ik,:) = Agap(ibnd,ik,itemp) 
              ENDIF
           ENDDO ! ibnd
        ENDDO ! ik
        Deltaold(:) = gap(itemp)
     ENDIF
     !
     IF ( .not. ALLOCATED(Gp) ) ALLOCATE( Gp(nsw,nqstep) )
     IF ( .not. ALLOCATED(Gm) ) ALLOCATE( Gm(nsw,nqstep) )
     DO iw = 1, nsw ! loop over omega
        DO iwp = 1, nqstep ! loop over omega_prime
           CALL gamma_acont( ws(iw), ws(iwp), estemp(itemp), rgammap, rgammam )
           Gp(iw,iwp) = rgammap
           Gm(iw,iwp) = rgammam
        ENDDO
     ENDDO
     CALL kernel_aniso_iaxis_analytic_cont( itemp )
     CALL eliashberg_memlt_aniso_acon
     IF ( .not. lacon_fly ) CALL evaluate_a2fij
  ENDIF
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  !
  DO ik = lower_bnd, upper_bnd
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iq = 1, nqfs(ik)
              ! iq0 - index of q-point on the full q-mesh
              iq0 = ixqfs(ik,iq)
              DO jbnd = 1, nbndfs
                 IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                    !
                    IF ( lacon_fly ) THEN ! evaluate a2fij on the fly
                       DO imode = 1, nmodes
                          IF ( wf(imode,iq0) .gt. eps_acustic ) THEN
                             DO iwph = 1, nqstep
                                weight  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / degaussq, 0 ) / degaussq
                                a2f_ = weight * dosef * g2(ik,iq,ibnd,jbnd,imode)
                             ENDDO ! iwph
                          ENDIF ! wf
                       ENDDO ! imode
                    ENDIF
                    !
                    weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                    DO iw = 1, nsw ! loop over omega
                       DO iwp = 1, nqstep ! loop over omega_prime
                          !
                          i = iw + iwp - 1
                          IF ( i .le. nsw ) THEN
                             root = sqrt(   AZnormp(jbnd,ixkqf(ik,iq0),i)**2.d0 & 
                                          * ( ws(i)**2.d0 - ADeltap(jbnd,ixkqf(ik,iq0),i)**2.d0 ) )
                             IF ( aimag(root) .lt. 0.d0 ) THEN 
                                esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / conjg(root)
                             ELSE  
                                esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / root
                             ENDIF
                             IF ( lacon_fly ) THEN 
                                esqrt = esqrt * weight * Gp(iw,iwp) * a2f_
                             ELSE
                                esqrt = esqrt * weight * Gp(iw,iwp) * a2fij(ik,iq,ibnd,jbnd,iwp) 
                             ENDIF
                             AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) - ws(i) * esqrt 
                             ADelta(ibnd,ik,iw) = ADelta(ibnd,ik,iw) - ADeltap(jbnd,ixkqf(ik,iq0),i) * esqrt
                          ENDIF
                          ! 
                          i = abs(iw - iwp) + 1
                          root = sqrt(   AZnormp(jbnd,ixkqf(ik,iq0),i)**2.d0 & 
                                       * ( ws(i)**2.d0 - ADeltap(jbnd,ixkqf(ik,iq0),i)**2.d0 ) )
                          IF ( aimag(root) .lt. 0.d0 ) THEN 
                             esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / conjg(root)
                          ELSE  
                             esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / root
                          ENDIF
                          esqrt = esqrt * weight * Gm(iw,iwp) * a2fij(ik,iq,ibnd,jbnd,iwp)
                          IF ( iw .lt. iwp ) THEN 
                             AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) - ws(i) * esqrt 
                          ELSE
                             AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) + ws(i) * esqrt 
                          ENDIF
                          ADelta(ibnd,ik,iw) = ADelta(ibnd,ik,iw) + ADeltap(jbnd,ixkqf(ik,iq0),i) * esqrt
                       ENDDO ! iwp
                    ENDDO ! iw
                 ENDIF
              ENDDO ! jbnd
           ENDDO ! iq
           DO iw = 1, nsw ! loop over omega
              AZnorm(ibnd,ik,iw) = - estemp(itemp) * AZsumi(ibnd,ik,iw) + ci * AZnorm(ibnd,ik,iw) * dwsph
              ADelta(ibnd,ik,iw) =   estemp(itemp) * ADsumi(ibnd,ik,iw) + ci * ADelta(ibnd,ik,iw) * dwsph
           ENDDO ! iw
        ENDIF
     ENDDO ! ibnd
  ENDDO ! ik
  !
  ! collect contributions from all pools 
  CALL mp_sum( AZnorm, inter_pool_comm )
  CALL mp_sum( ADelta, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF (mpime .eq. ionode_id) THEN
    absdelta = 0.d0
    reldelta = 0.d0
    DO iw = 1, nsw ! loop over omega
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                Znorm(iw) = Znorm(iw) + weight * AZnorm(ibnd,ik,iw)
                Delta(iw) = Delta(iw) + weight * ADelta(ibnd,ik,iw)
                AZnorm(ibnd,ik,iw) = 1.d0 + pi * AZnorm(ibnd,ik,iw) / ws(iw)
                ADelta(ibnd,ik,iw) = pi * ADelta(ibnd,ik,iw) / AZnorm(ibnd,ik,iw)
             ENDIF
          ENDDO ! ibnd                   
       ENDDO ! ik
       Znorm(iw) = 1.0d0 + pi * Znorm(iw) / ws(iw)
       Delta(iw) = pi * Delta(iw) / Znorm(iw)
       reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
       absdelta = absdelta + abs( Delta(iw) ) 
    ENDDO ! iw
    errdelta = reldelta / absdelta
    Deltaold(:) = Delta(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
                 '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
                 '   Re[Delta(1)] = ', real(Delta(1))
    !
    IF ( errdelta .lt. conv_thr_racon ) conv = .true.
    IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
       cname = 'acon'
       CALL eliashberg_write_cont_raxis( itemp, cname )
    ENDIF
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
       CALL errore('analytic_cont_aniso_iaxis_to_raxis','increase nsiter or reduce conv_thr_racon',1)
    ENDIF
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_bcast( conv, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF ( conv .OR. iter .eq. nsiter ) THEN
     !
     IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
     IF( ALLOCATED(Gp) )       DEALLOCATE(Gp)
     IF( ALLOCATED(Gm) )       DEALLOCATE(Gm)
     IF( ALLOCATED(ADsumi) )   DEALLOCATE(ADsumi)
     IF( ALLOCATED(AZsumi) )   DEALLOCATE(AZsumi)
     !
     ! remove memory allocated for Deltaold, Gp, Gm, ADsumi, AZsumi
     imelt = 2 * nsw + 2 * nqstep * nsw + 2 * ( upper_bnd - lower_bnd + 1 ) * nbndfs * nsw
     CALL mem_size_eliashberg( -imelt )
     !
     IF ( .not. lacon_fly ) THEN
        !
        IF ( ALLOCATED(a2fij) ) DEALLOCATE(a2fij)
        !
        ! remove memory allocated for a2fij
        imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
        CALL mem_size_eliashberg( -imelt )
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE analytic_cont_aniso_iaxis_to_raxis
  !
  !-----------------------------------------------------------------------
  SUBROUTINE pade_cont_aniso_iaxis_to_raxis( itemp, N, conv )
  !-----------------------------------------------------------------------
  !
  ! This routine uses pade approximants to continue the anisotropic Eliashberg equations 
  ! from the imaginary-axis to the real-axis
  !
  ! input
  !
  ! itemp  - temperature point
  ! N      - number Pade approximants
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : fsthick
  USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Agap, Delta, Znorm, & 
                            ADelta, AZnorm, ADeltai, AZnormi, &              
                            wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : cone, ci
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, N, ik, ibnd, lower_bnd, upper_bnd, imelt
  REAL(DP) :: absdelta, reldelta, errdelta, weight
  COMPLEX(DP) :: omega, padapp
  COMPLEX(DP), ALLOCATABLE :: a(:), b(:), z(:), u(:), v(:), Deltaold(:)
  LOGICAL :: conv
  CHARACTER (len=256) :: cname
  !
  ! get the size of required allocated memory for 
  ! a, b, z, u, v, Delta, Znorm, Deltaold, ADelta, AZnorm
  imelt = 2 * 5 * N + 2 * ( 3 + 2 * nbndfs * nkfs ) * nsw
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
  IF ( .not. ALLOCATED(ADelta) )   ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
  IF ( .not. ALLOCATED(AZnorm) )   ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
  IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
  IF ( .not. ALLOCATED(a) )        ALLOCATE( a(N) )
  IF ( .not. ALLOCATED(b) )        ALLOCATE( b(N) )
  IF ( .not. ALLOCATED(z) )        ALLOCATE( z(N) )
  IF ( .not. ALLOCATED(u) )        ALLOCATE( u(N) )
  IF ( .not. ALLOCATED(v) )        ALLOCATE( v(N) )
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
  Deltaold(:) = gap(itemp)
  absdelta = 0.d0
  reldelta = 0.d0
  a(:) = (0.d0, 0.d0)
  b(:) = (0.d0, 0.d0)
  z(:) = (0.d0, 0.d0)
  u(:) = (0.d0, 0.d0)
  v(:) = (0.d0, 0.d0)
  !
  CALL fkbounds( nkfs, lower_bnd, upper_bnd )
  !
  DO ik = lower_bnd, upper_bnd
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iw = 1, N
              z(iw) = ci * wsi(iw)
              u(iw) = cone * ADeltai(ibnd,ik,iw) 
              v(iw) = cone * AZnormi(ibnd,ik,iw)
           ENDDO
           CALL pade_coeff( N, z, u, a )
           CALL pade_coeff( N, z, v, b )
           DO iw = 1, nsw
              omega = cone * ws(iw)
              CALL pade_eval( N, z, a, omega, padapp )
              ADelta(ibnd,ik,iw) = padapp
              CALL pade_eval( N, z, b, omega, padapp )
              AZnorm(ibnd,ik,iw) = padapp
           ENDDO
        ENDIF
     ENDDO ! ibnd
  ENDDO ! ik
  !
  ! collect contributions from all pools 
  CALL mp_sum( AZnorm, inter_pool_comm )
  CALL mp_sum( ADelta, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF (mpime .eq. ionode_id) THEN
    DO iw = 1, nsw ! loop over omega
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                Znorm(iw) = Znorm(iw) + weight * AZnorm(ibnd,ik,iw)
                Delta(iw) = Delta(iw) + weight * ADelta(ibnd,ik,iw)
             ENDIF
          ENDDO ! ibnd                   
       ENDDO ! ik
       reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) )
       absdelta = absdelta + abs( Delta(iw) )
    ENDDO ! iw
    errdelta = reldelta / absdelta
    !
    IF ( errdelta .gt. 0.d0 ) THEN
       conv = .true.
       WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'pade = ', N, & 
              '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
              '   Re[Delta(1)] = ', real(Delta(1))
       cname = 'pade'
       CALL eliashberg_write_cont_raxis( itemp, cname )
    ENDIF
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_bcast( conv, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
  IF( ALLOCATED(a) )        DEALLOCATE(a)
  IF( ALLOCATED(b) )        DEALLOCATE(b)
  IF( ALLOCATED(z) )        DEALLOCATE(z)
  IF( ALLOCATED(u) )        DEALLOCATE(u)
  IF( ALLOCATED(v) )        DEALLOCATE(v)
  !
  ! remove memory allocated for Deltaold, a, b, z, u, v
  imelt = 2 * ( nsw + 5 * N )
  CALL mem_size_eliashberg( -imelt )
  !
  RETURN
  !
  END SUBROUTINE pade_cont_aniso_iaxis_to_raxis
  !
  !-----------------------------------------------------------------------
