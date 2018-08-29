  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE superconductivity
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutine linked with superconductivity using  
  !! the isotropic or anisotropic Eliashberg formalism. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !                                                                            
    !----------------------------------------------------------------------- 
    SUBROUTINE eliashberg_init
    !-----------------------------------------------------------------------
    !
    ! This routine initializes the control variables needed to solve the eliashberg 
    ! equations
    !
    USE kinds,           ONLY : DP
    USE io_global,       ONLY : stdout
    USE epwcom,          ONLY : eliashberg, nkf1, nkf2, nkf3, nsiter, nqstep, &
                                nqf1, nqf2, nqf3, ntempxx, nswi, nswfc, nswc, &
                                nstemp, muc, lreal, lpade, liso, limag, laniso, &
                                lacon, kerwrite, kerread, imag_read, fila2f, temps, &
                                wsfc, wscut, tempsmin, tempsmax, rand_q, rand_k
    USE constants_epw,   ONLY : pi, kelvin2eV, eps6
    USE eliashbergcom,   ONLY : estemp, nsw, nsiw, dwsph, wsphmax, wsph
    !
    IMPLICIT NONE
    !
    INTEGER :: itemp, iwph, imelt
    REAL(DP) :: dtemp
    !
    IF ( eliashberg .AND. liso .AND. laniso ) CALL errore('eliashberg_init', &
         'liso or laniso needs to be true',1)
    IF ( .not.eliashberg .AND. liso ) CALL errore('eliashberg_init', &
         'liso requires eliashberg true',1)
    IF ( .not.eliashberg .AND. laniso ) CALL errore('eliashberg_init', &
         'laniso requires eliashberg true',1)
    IF ( laniso .and. (fila2f .ne. ' ') ) &
         CALL errore('eliashberg_init', 'anisotropic case can not use fila2f',1)
    IF ( eliashberg .AND. lreal .AND. laniso ) CALL errore('eliashberg_init', &
         'lreal is implemented only for the isotriopic case',1)
    IF ( eliashberg .AND. lreal .AND. limag )  & 
         CALL errore('eliashberg_init', 'lreal or limag needs to be true',1)
    IF ( eliashberg .AND. lreal .AND. lacon )  &
         CALL errore('eliashberg_init', 'lreal or lacon needs to be true',1)
    IF ( eliashberg .AND. lreal .AND.  lpade ) &
         CALL errore('eliashberg_init', 'lreal or lpade needs to be true',1)
    IF ( eliashberg .AND. imag_read .AND. .not.limag .AND. .not.laniso ) &
         CALL errore('eliashberg_init', 'imag_read requires limag true and laniso true',1)
    IF ( eliashberg .AND. lpade .AND. .not.limag ) &                
         CALL errore('eliashberg_init', 'lpade requires limag true',1)
    IF ( eliashberg .AND. lacon .AND. (.not.limag .OR. .not.lpade ) ) & 
         CALL errore('eliashberg_init', 'lacon requires both limag and lpade true',1)
    IF ( eliashberg .AND. lreal .AND. (kerread .AND. kerwrite) ) & 
         CALL errore('eliashberg_init', 'kerread cannot be used with kerwrite',1)
    IF ( eliashberg .AND. lreal .AND. (.not.kerread .AND. .not.kerwrite) ) &
         CALL errore('eliashberg_init', 'kerread or kerwrite must be true',1)
    IF ( eliashberg .AND. lreal .AND. wsfc .gt. wscut ) CALL errore('eliashberg_init', &
         'wsfc should be .lt. wscut',1)
    IF ( eliashberg .AND. lreal .AND. wsfc .lt. 0.d0 ) CALL errore('eliashberg_init', &
         'wsfc should be .gt. 0.d0',1)
    IF ( eliashberg .AND. nswi .gt. 0 .AND. .not.limag ) &
         CALL errore('eliashberg_init', 'nswi requires limag true',1)
    IF ( eliashberg .AND. nswi .lt. 0 ) CALL errore('eliashberg_init', &
         'nswi should be .gt. 0',1)
    IF ( eliashberg .AND. wscut .lt. 0.d0 ) &
         CALL errore('eliashberg_init', 'wscut should be .gt. 0.d0',1)
    IF ( eliashberg .AND. nstemp .lt. 1 ) CALL errore('eliashberg_init', &
         'wrong number of nstemp',1)
    IF ( eliashberg .AND. maxval(temps(:)) .gt. 0.d0 .AND. & 
         tempsmin .gt. 0.d0 .AND. tempsmax .gt. 0.d0 ) &
         CALL errore('eliashberg_init', & 
         'define either (tempsmin and tempsmax) or temp(:)',1)
    IF ( eliashberg .AND. tempsmax .lt. tempsmin ) &
         CALL errore('eliashberg_init', & 
         'tempsmax should be greater than tempsmin',1)
    IF ( eliashberg .AND. nsiter .lt. 1 ) CALL errore('eliashberg_init', &
         'wrong number of nsiter',1)
    IF ( eliashberg .AND. muc .lt. 0.d0 ) CALL errore('eliashberg_init', &
         'muc should be .ge. 0.d0',1) 
    IF ( eliashberg .and. (rand_k .OR. rand_q ) .and. (fila2f .eq. ' ') ) &
         CALL errore('eliashberg_init', 'eliashberg requires a uniform grid when fila2f is not used',1)
    IF ( eliashberg .and. (mod(nkf1,nqf1) .ne. 0 .OR. mod(nkf2,nqf2) &
         .ne. 0 .OR. mod(nkf3,nqf3) .ne. 0 ) .and. (fila2f .eq. ' ') ) &
         CALL errore('eliashberg_init', &
         'eliashberg requires nkf1,nkf2,nkf3 to be multiple of nqf1,nqf2,nqf3 when fila2f is not used',1)
    !
    DO itemp = 1, ntempxx
      IF (temps(itemp) .gt. 0.d0) THEN
        nstemp = itemp
      ENDIF
    ENDDO
    !
    IF ( .not. ALLOCATED(estemp) ) ALLOCATE( estemp(nstemp) )
    estemp(:) = 0.d0
    !
    ! go from K to eV
    IF ( maxval(temps(:)) .gt. 0.d0 ) THEN
      DO itemp= 1, nstemp 
        estemp(itemp) = temps(itemp) * kelvin2eV
      ENDDO
    ELSE
      IF ( nstemp .eq. 1 ) THEN
        estemp(1) = tempsmin * kelvin2eV
      ELSE
        dtemp = ( tempsmax - tempsmin ) * kelvin2eV / dble(nstemp-1)
        DO itemp = 1, nstemp
          estemp(itemp) = tempsmin * kelvin2eV + dble(itemp-1) * dtemp
        ENDDO
      ENDIF
    ENDIF
    !
    IF ( lreal ) THEN
      !
      IF ( ABS(wsfc) < eps6 .OR. ABS(wscut) < eps6 .OR. nswfc .eq. 0 .OR. nswc .eq. 0 ) THEN 
        wsfc  =  5.d0 * wsphmax
        wscut = 15.d0 * wsphmax
        nswfc = 5 * nqstep
        nswc  = 2 * nqstep
      ENDIF
      nsw = nswfc + nswc  
      WRITE(stdout,'(5x,a7,f12.6,a11,f12.6)') 'wsfc = ', wsfc, '   wscut = ', wscut
      WRITE(stdout,'(5x,a8,i8,a10,i8,a9,i8)') 'nswfc = ', nswfc, '   nswc = ', nswc, & 
                                               '   nsw = ', nsw 
      IF ( nsw .eq. 0 ) CALL errore('eliashberg_setup','wrong number of nsw',1)
      !
    ELSEIF ( limag ) THEN
      !
      IF ( .not. ALLOCATED(nsiw) ) ALLOCATE( nsiw(nstemp) )
      nsiw(:) = 0
      !
      IF ( nswi .gt. 0 ) THEN
        nsiw(:) = nswi
      ELSEIF ( wscut .gt. 0.d0 ) THEN
        DO itemp = 1, nstemp
           nsiw(itemp) = int(0.5d0 * ( wscut / pi / estemp(itemp) - 1.d0 )) + 1
        ENDDO
      ELSEIF ( nswi .gt. 0 .AND. wscut .gt. 0.d0 ) THEN
        nsiw(:) = nswi
        WRITE(stdout,'(5x,a)') 'when nswi .gt. 0, wscut is not used for limag=.true.'
      ENDIF
      !
      IF ( ABS(wscut) < eps6 ) THEN 
        wscut = 10.d0 * wsphmax
      ENDIF
      ! 
      IF ( lpade .OR. lacon ) THEN 
        nsw = nqstep * nint(wscut/wsphmax)
        IF ( nsw .eq. 0 ) CALL errore('eliashberg_setup','wrong number of nsw',1)
      ENDIF
      !
    ENDIF
    !
    ! create phonon grid 
    !
    !dwsph = wsphmax / dble(nqstep-1)
    dwsph = wsphmax / dble(nqstep)
    IF ( .not. ALLOCATED(wsph) ) ALLOCATE( wsph(nqstep) )
    wsph(:) = 0.d0
    DO iwph = 1, nqstep
      !wsph(iwph) = dble(iwph-1) * dwsph
      wsph(iwph) = dble(iwph) * dwsph
    ENDDO
    !
    ! memory allocated for wsph, estemp
    imelt = nqstep + nstemp
    CALL mem_size_eliashberg( imelt )
    !
    ! memory allocated for nsiw
    imelt = nstemp
    CALL mem_integer_size_eliashberg( imelt )
    !
    RETURN
    !
    END SUBROUTINE eliashberg_init
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2f_lambda
    !-----------------------------------------------------------------------
    !
    ! computes the isotropic spectral function a2F(w), total lambda, and 
    ! distribution of lambda
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_epw,        ONLY : iua2ffil, iudosfil, iufillambda, iufillambdaFS
    USE io_files,      ONLY : prefix
    USE phcom,         ONLY : nmodes
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE elph2,         ONLY : nqtotf, wqf, wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq, delta_qsmear, nqsmear, & 
                              degaussw, nkf1, nkf2, nkf3
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, ixkqf, ixqfs, nqfs, w0g, ekfs, ef0, dosef, wsph, &
                              wkfs, dwsph, a2f_iso, ixkff
    USE constants_epw, ONLY : ryd2ev
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE superconductivity_aniso, ONLY : lambdar_aniso_ver1
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ik, iq, iq0, iwph, ibnd, jbnd, imode, lower_bnd, upper_bnd, &
         ismear, ibin, nbin, nbink, i, j, k
    REAL(DP) :: weight, weightq, l_sum, lambda_eph, lambda_max(npool), & 
         sigma, dbin, dbink, x1, x2, x3
    REAL(DP), ALLOCATABLE :: a2f(:,:), phdos(:,:), l_a2f(:), lambda_k(:,:), &
         lambda_k_bin(:), lambda_pairs(:), a2f_modeproj(:,:), phdos_modeproj(:,:)
    REAL(DP), EXTERNAL :: w0gauss
    CHARACTER (len=256) :: name1
    ! 
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    npool = 1
    my_pool_id = 0
#endif
    !
    ! degaussq is read from the input file in meV and converted to Ryd in epw_readin.f90
    ! go from Ryd to eV
    degaussq = degaussq * ryd2ev
    delta_qsmear = delta_qsmear * ryd2ev
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    IF ( .not. ALLOCATED(a2f_iso) ) ALLOCATE( a2f_iso(nqstep) )
    IF ( .not. ALLOCATED(a2f) )     ALLOCATE( a2f(nqstep,nqsmear) )
    IF ( .not. ALLOCATED(a2f_modeproj) ) ALLOCATE( a2f_modeproj(nmodes,nqstep) )
    a2f_iso(:) = 0.d0
    a2f(:,:) = 0.d0
    a2f_modeproj(:,:) = 0.d0
    !
    ! RM - the 0 index in k is required when printing out values of lambda_k 
    ! When the k-point is outside the Fermi shell, ixkff(ik)=0
    IF ( .not. ALLOCATED(lambda_k) ) ALLOCATE(lambda_k(0:nkfs,nbndfs))
    lambda_k(:,:) = 0.d0
    !
    l_sum = 0.d0
    lambda_max(:) = 0.d0
    DO ismear = 1, nqsmear
       sigma = degaussq + (ismear-1) * delta_qsmear
       DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN 
                DO iq = 1, nqfs(ik)
                   ! iq0 - index of q-point on the full q-mesh
                   iq0 = ixqfs(ik,iq)
                   DO jbnd = 1, nbndfs
                      IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                         weight = wkfs(ik) * wqf(iq) * w0g(ibnd,ik) * w0g(jbnd,ixkqf(ik,iq0))
                         lambda_eph = 0.d0
                         DO imode = 1, nmodes
                            IF ( wf(imode,iq0) .gt. eps_acustic ) THEN
                               IF ( ismear .eq. 1 ) THEN 
                                  lambda_eph = lambda_eph + g2(ik,iq,ibnd,jbnd,imode) / wf(imode,iq0)
                               ENDIF
                               DO iwph = 1, nqstep
                                  weightq  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / sigma, 0 ) / sigma
                                  a2f(iwph,ismear) = a2f(iwph,ismear) + weight * weightq * g2(ik,iq,ibnd,jbnd,imode)
                                  IF ( ismear .eq. 1 ) THEN
                                     a2f_modeproj(imode,iwph) = a2f_modeproj(imode,iwph) +&
                                         weight * weightq * g2(ik,iq,ibnd,jbnd,imode)
                                  ENDIF
                               ENDDO ! iwph
                            ENDIF ! wf
                         ENDDO ! imode
                         IF ( ismear .eq. 1 .AND. lambda_eph .gt. 0.d0 ) THEN
                            l_sum = l_sum + weight * lambda_eph
                            weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) 
                            lambda_k(ik,ibnd) = lambda_k(ik,ibnd) + weight * lambda_eph
                            IF ( lambda_eph .gt. lambda_max(my_pool_id+1) ) THEN
                               lambda_max(my_pool_id+1) = lambda_eph
                            ENDIF
                         ENDIF
                      ENDIF ! ekq
                   ENDDO ! jbnd
                ENDDO ! iq
             ENDIF ! ekk
          ENDDO ! ibnd
       ENDDO ! ik
    ENDDO ! ismear
    !
    a2f(:,:) = 0.5d0 * a2f(:,:) / dosef
    a2f_modeproj(:,:) = 0.5d0 * a2f_modeproj(:,:) / dosef
    l_sum = l_sum / dosef
    lambda_k(:,:) = 2.d0 * lambda_k(:,:)
    lambda_max(:) = 2.d0 * dosef * lambda_max(:)
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( l_sum, inter_pool_comm )
    CALL mp_sum( a2f, inter_pool_comm )
    CALL mp_sum( a2f_modeproj, inter_pool_comm )
    CALL mp_sum( lambda_max, inter_pool_comm )
    CALL mp_sum( lambda_k, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      OPEN( unit = iua2ffil, file = TRIM(prefix)//".a2f", form = 'formatted')
      OPEN( unit = iudosfil, file = TRIM(prefix)//".phdos", form = 'formatted')
      !
      IF ( .not. ALLOCATED(phdos) )     ALLOCATE( phdos(nqstep,nqsmear) )
      IF ( .not. ALLOCATED(phdos_modeproj) ) ALLOCATE( phdos_modeproj(nmodes,nqstep) )
      phdos(:,:) = 0.d0
      phdos_modeproj(:,:) = 0.d0
      !
      DO ismear = 1, nqsmear
         sigma = degaussq + (ismear-1) * delta_qsmear
         DO iq = 1, nqtotf
            DO imode = 1, nmodes
               IF ( wf(imode,iq) .gt. eps_acustic ) THEN
                  DO iwph = 1, nqstep
                     weightq  = w0gauss( ( wsph(iwph) - wf(imode,iq)) / sigma, 0 ) / sigma
                     phdos(iwph,ismear) = phdos(iwph,ismear) + wqf(iq) * weightq
                     IF ( ismear .eq. 1 ) THEN
                        phdos_modeproj(imode,iwph) = phdos_modeproj(imode,iwph) + wqf(iq) * weightq
                     ENDIF
                  ENDDO ! iwph
               ENDIF ! wf
            ENDDO ! imode
         ENDDO ! iq
      ENDDO ! ismear
      !
      IF ( .not. ALLOCATED(l_a2f) )   ALLOCATE( l_a2f(nqsmear) )
      l_a2f(:) = 0.d0
      !
      DO ismear = 1, nqsmear
         DO iwph = 1, nqstep
            l_a2f(ismear) = l_a2f(ismear) + a2f(iwph,ismear) / wsph(iwph)
            ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
            IF (ismear .eq. nqsmear) WRITE (iua2ffil,'(f12.7,15f12.7)') wsph(iwph)*1000.d0, a2f(iwph,:)
            IF (ismear .eq. nqsmear) WRITE (iudosfil,'(f12.7,15f15.7)') wsph(iwph)*1000.d0, phdos(iwph,:)/1000.d0
         ENDDO
         l_a2f(ismear) = 2.d0 * l_a2f(ismear) * dwsph
      ENDDO
      !
      WRITE(iua2ffil,*) "Integrated el-ph coupling"
      WRITE(iua2ffil,'("  #         ", 15f12.7)') l_a2f(:)
      WRITE(iua2ffil,*) "Phonon smearing (meV)" 
      WRITE(iua2ffil,'("  #         ", 15f12.7)') ( (degaussq+(ismear-1)*delta_qsmear)*1000.d0,ismear=1,nqsmear )
      WRITE(iua2ffil,'(" Electron smearing (eV)", f12.7)') degaussw
      WRITE(iua2ffil,'(" Fermi window (eV)", f12.7)') fsthick
      WRITE(iua2ffil,'(" Summed el-ph coupling ", f12.7)') l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      a2f_iso(:) = a2f(:,1)
      OPEN( unit = iua2ffil, file = TRIM(prefix)//".a2f_iso", form = 'formatted')
      OPEN( unit = iudosfil, file = TRIM(prefix)//".phdos_proj", form = 'formatted')
      DO iwph = 1, nqstep
         ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
         WRITE(iua2ffil,'(f12.7,100f12.7)') wsph(iwph)*1000.d0, a2f_iso(iwph), a2f_modeproj(:,iwph)
         WRITE(iudosfil,'(f12.7,100f15.7)') wsph(iwph)*1000.d0, phdos(iwph,1)/1000.d0, phdos_modeproj(:,iwph)/1000.d0
      ENDDO
      WRITE(iua2ffil,'(a,f18.7,a,f18.7)') 'lambda_int = ', l_a2f(1), '   lambda_sum = ',l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      IF ( ALLOCATED(phdos) )          DEALLOCATE( phdos )
      IF ( ALLOCATED(phdos_modeproj) ) DEALLOCATE( phdos_modeproj )
      IF ( ALLOCATED(l_a2f) )          DEALLOCATE( l_a2f )
      !
    ENDIF
    !
    CALL mp_bcast( a2f_iso, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(a2f) )            DEALLOCATE( a2f )
    IF ( ALLOCATED(a2f_modeproj) )   DEALLOCATE( a2f_modeproj )
    !
    nbink = int( 1.25d0 * maxval(lambda_k(:,:)) / 0.005d0 )
    dbink = 1.25d0 * maxval(lambda_k(:,:)) / dble(nbink)
    IF ( .not. ALLOCATED(lambda_k_bin) ) ALLOCATE ( lambda_k_bin(nbink) )
    lambda_k_bin(:) = 0.d0
    !
    !SP : Should be initialized
    nbin = 0
    dbin = 0.0_DP
    !
    IF ( iverbosity == 2 ) THEN
       nbin = int( 1.25d0 * maxval(lambda_max(:)) / 0.005d0 )
       dbin = 1.25d0 * maxval(lambda_max(:)) / dble(nbin)
       IF ( .not. ALLOCATED(lambda_pairs) ) ALLOCATE ( lambda_pairs(nbin) )
       lambda_pairs(:) = 0.d0
    ENDIF
    ! 
    WRITE(stdout,'(5x,a13,f21.7,a18,f21.7)') 'lambda_max = ', maxval(lambda_max(:)), '   lambda_k_max = ', maxval(lambda_k(:,:))
    WRITE(stdout,'(a)') ' '
    !
    lambda_k(:,:) = 0.d0
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, 0.d0, lambda_eph )
                      lambda_k(ik,ibnd) = lambda_k(ik,ibnd) +  weight * lambda_eph
                      IF ( iverbosity == 2 ) THEN
                         DO ibin = 1, nbin
                            sigma = 1.d0 * dbin
                            weight = w0gauss( ( lambda_eph - dble(ibin) * dbin ) / sigma, 0 ) / sigma
                            lambda_pairs(ibin) = lambda_pairs(ibin) + weight
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
             DO ibin = 1, nbink
                sigma = 1.d0 * dbink
                weight = w0gauss( ( lambda_k(ik,ibnd) - dble(ibin) * dbink ) / sigma, 0 ) / sigma
                lambda_k_bin(ibin) = lambda_k_bin(ibin) + weight
             ENDDO
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools 
    CALL mp_sum( lambda_k, inter_pool_comm )
    IF ( iverbosity .eq. 2 ) THEN  
      CALL mp_sum( lambda_pairs, inter_pool_comm )
    ENDIF
    CALL mp_sum( lambda_k_bin, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      ! SP: Produced if user really wants it 
      IF ( iverbosity .eq. 2 ) THEN
        OPEN(unit = iufillambda, file = TRIM(prefix)//".lambda_aniso", form = 'formatted')
        WRITE(iufillambda,'(2a12,2a7)') '# enk-e0[eV]','  lambda_nk','# kpt','# band'
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 WRITE(iufillambda,'(2f12.7,2i7)') ekfs(ibnd,ik) - ef0, lambda_k(ik,ibnd), ik, ibnd
              ENDIF
           ENDDO
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      OPEN(unit = iufillambda, file = TRIM(prefix)//".lambda_k_pairs", form = 'formatted')
      WRITE(iufillambda,'(a12,a30)') '# lambda_nk','  \rho(lambda_nk) scaled to 1'
      DO ibin = 1, nbink
        WRITE(iufillambda,'(2f21.7)') dbink*dble(ibin), lambda_k_bin(ibin)/maxval(lambda_k_bin(:))
      ENDDO
      CLOSE(iufillambda)
      !
      ! SP: Produced if user really wants it 
      IF ( iverbosity == 2 ) THEN  
        OPEN( unit = iufillambda, file = TRIM(prefix)//".lambda_pairs", form = 'formatted')
      WRITE(iufillambda,'(a12,a30)') "# lambda_nk,n'k'", "  \rho(lambda_nk,n'k') scaled to 1"
        DO ibin = 1, nbin
          WRITE(iufillambda,'(2f21.7)') dbin*dble(ibin), lambda_pairs(ibin)/maxval(lambda_pairs(:))
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
      !
      ! RM - If the k-point is outside the Fermi shell,
      ! ixkff(ik)=0 and lambda_k(0,ibnd) = 0.0
      !
      IF ( iverbosity .eq. 2 ) THEN
        !
        DO ibnd = 1, nbndfs
          !
          IF ( ibnd < 10 ) THEN
            WRITE(name1,'(a,a8,i1,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSEIF ( ibnd < 100 ) THEN
            WRITE(name1,'(a,a8,i2,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSEIF( ibnd < 1000 ) THEN
            WRITE(name1,'(a,a8,i3,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSE 
            CALL errore( 'eliashberg_setup', 'Too many bands ',1)  
          ENDIF  
          !  
          OPEN(iufillambdaFS, file=name1, form='formatted')
          WRITE(iufillambdaFS,*) 'Cubfile created from EPW calculation'
          WRITE(iufillambdaFS,*) 'lambda'
          WRITE(iufillambdaFS,'(i5,3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf1, (bg(i,1)/DBLE(nkf1),i=1,3)
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf2, (bg(i,2)/DBLE(nkf2),i=1,3)
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf3, (bg(i,3)/DBLE(nkf3),i=1,3)
          WRITE(iufillambdaFS,'(i5,4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS,'(6f12.6)') ( lambda_k(ixkff(ik),ibnd), ik=1,nkf1*nkf2*nkf3 )
          CLOSE(iufillambdaFS)
        ENDDO
        !
      ENDIF
      !
      ! SP & RM : Write on file the lambda close to the Fermi surface along with 
      ! Cartesian coordinate, band index, energy distance from Fermi level
      ! and lambda value.
      !
      OPEN(unit = iufillambdaFS, file = TRIM(prefix)//".lambda_FS", form='formatted')
      WRITE(iufillambdaFS,'(a75)') '#               k-point                  Band Enk-Ef [eV]            lambda'
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = k + (j-1)*nkf3 + (i-1)*nkf2*nkf3
               IF ( ixkff(ik) .gt. 0 ) THEN
                  DO ibnd = 1, nbndfs
                     ! SP: Here take a 0.2 eV interval around the FS.
                     IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. fsthick ) THEN
                     !IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. 0.2 ) THEN
                        x1 = bg(1,1)*(i-1)/nkf1+bg(1,2)*(j-1)/nkf2+bg(1,3)*(k-1)/nkf3
                        x2 = bg(2,1)*(i-1)/nkf1+bg(2,2)*(j-1)/nkf2+bg(2,3)*(k-1)/nkf3
                        x3 = bg(3,1)*(i-1)/nkf1+bg(3,2)*(j-1)/nkf2+bg(3,3)*(k-1)/nkf3
                        WRITE(iufillambdaFS,'(3f12.6,i8,f12.6,f24.15)') x1, x2, x3, ibnd, &
                                         ekfs(ibnd,ixkff(ik))-ef0, lambda_k(ixkff(ik),ibnd)
                     ENDIF
                  ENDDO ! ibnd
               ENDIF
            ENDDO  ! k
         ENDDO ! j
      ENDDO ! i
      CLOSE(iufillambdaFS)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(lambda_k) )     DEALLOCATE(lambda_k)
    IF ( ALLOCATED(lambda_pairs) ) DEALLOCATE(lambda_pairs)
    IF ( ALLOCATED(lambda_k_bin) ) DEALLOCATE(lambda_k_bin)
    !
    RETURN
    !
    END SUBROUTINE evaluate_a2f_lambda
    ! 
    !----------------------------------------------------------------------- 
    SUBROUTINE estimate_tc_gap
    !-----------------------------------------------------------------------
    !
    ! this subroutine estimates the Tc using Allen-Dynes formula and 
    ! the BCS superconducting gap as the initial guess for Delta 
    !  
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, muc, tempsmin, tempsmax, temps
    USE eliashbergcom, ONLY : wsph, dwsph, a2f_iso, gap0
    USE constants_epw, ONLY : kelvin2eV
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !  
    IMPLICIT NONE
    !  
    INTEGER :: iwph
    REAL(DP):: l_a2f, logavg, tc
    !
    IF ( mpime .eq. ionode_id ) THEN
      l_a2f  = 0.0d0 
      logavg = 0.0d0
      DO iwph = 1, nqstep
         l_a2f  = l_a2f  + a2f_iso(iwph) / wsph(iwph)
         logavg = logavg + a2f_iso(iwph) * log(wsph(iwph)) / wsph(iwph)
      ENDDO
      l_a2f  = l_a2f  * 2.d0 * dwsph
      logavg = logavg * 2.d0 * dwsph
      logavg = exp( logavg / l_a2f )
      WRITE(stdout,'(5x,a,f12.7)') 'Electron-phonon coupling strength = ', l_a2f
      WRITE(stdout,'(a)') ' '
      !
      ! Allen-Dynes estimate of Tc
      !
      tc = logavg / 1.2d0 * exp( - 1.04d0 * ( 1.d0 + l_a2f ) &
                                 / ( l_a2f - muc * ( 1.d0 + 0.62d0 * l_a2f ) ) )
      !
      ! initial guess for the gap edge using BCS superconducting ratio 3.52
      !
      gap0 = 3.52d0 * tc / 2.d0
      IF ( gap0 .le. 0.d0 ) CALL errore('estimate_tc_gap', &
         'initial guess for gap edge should be .gt. 0.d0',1)
      !
      ! tc in K
      !
      tc = tc / kelvin2eV
      WRITE(stdout,'(5x,a,f15.7,a,f10.5)') 'Estimated Allen-Dynes Tc = ', tc, ' K for muc = ', muc
      WRITE(stdout,'(a)') '  '
      WRITE(stdout,'(5x,a,f15.7,a)') 'Estimated BCS superconducting gap = ', gap0, ' eV'
      !
      IF ( tempsmin .gt. 1.3d0*tc .OR. minval(temps(:)) .gt. 1.3d0*tc ) THEN
         CALL errore('eliashberg_init','tempsmin or minval(temps) .gt. estimated Allen-Dynes 1.3*Tc',-1)
      ELSEIF ( tempsmax .gt. tc .OR. maxval(temps(:)) .gt. tc ) THEN
         WRITE(stdout,'(a)') '  '
         WRITE(stdout,'(5x,a)') 'WARNING WARNING WARNING '
         WRITE(stdout,'(a)') '  '
         WRITE(stdout,'(5x,a)') 'The code will crash for tempsmax much larger than Allen-Dynes Tc'
      ELSEIF ( tempsmax .gt. 1.5d0*tc .OR. maxval(temps(:)) .gt. 1.5d0*tc ) THEN
         CALL errore('eliashberg_init','tempsmax or maxval(temps) .gt. estimated Allen-Dynes 1.5*Tc',-1)
      ENDIF
      !
    ENDIF
    CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    RETURN
    !
    END SUBROUTINE estimate_tc_gap
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mem_size_eliashberg( imelt )
    !-----------------------------------------------------------------------
    !
    !  subroutine estimates the amount of memory taken up or 
    !  released by different arrays 
    !  if imelt > 0 memory is added
    !  if imelt < 0 memory is subtracted  
    !
    USE io_global,     ONLY : stdout
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : max_memlt
    USE eliashbergcom, ONLY : memlt_pool
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: imelt
    REAL(DP) :: rmelt
    ! 
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    my_pool_id = 0
#endif
    !
    rmelt = 0.d0
    rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
    rmelt = rmelt + memlt_pool(my_pool_id+1)
    !
    memlt_pool(:) = 0.d0
    memlt_pool(my_pool_id+1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum( memlt_pool, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
       CALL errore('mem_size_eliashberg', 'Size of required memory exceeds max_memlt',1)
    ELSEIF( maxval(memlt_pool(:)) .gt. 0.5d0*max_memlt ) THEN 
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE mem_size_eliashberg
    !                            
    !-----------------------------------------------------------------------
    SUBROUTINE mem_integer_size_eliashberg( imelt )
    !-----------------------------------------------------------------------
    !
    !  subroutine estimates the amount of memory taken up or 
    !  released by different arrays 
    !  if imelt > 0 memory is added
    !  if imelt < 0 memory is subtracted  
    !
    USE io_global,     ONLY : stdout
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : max_memlt
    USE eliashbergcom, ONLY : memlt_pool
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: imelt
    REAL(DP) :: rmelt
    !
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    my_pool_id = 0
#endif  
    !
    rmelt = 0.d0
    rmelt = dble(imelt) * 4.d0 / 1073741824.d0 ! 4 bytes per number, value in Gb
    rmelt = rmelt + memlt_pool(my_pool_id+1)
    !
    memlt_pool(:) = 0.d0
    memlt_pool(my_pool_id+1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum( memlt_pool, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
       CALL errore('mem_integer_size_eliashberg', 'Size of required memory exceeds max_memlt',1)
    ELSEIF( maxval(memlt_pool(:)) .gt. 0.5d0*max_memlt ) THEN
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE mem_integer_size_eliashberg
    !                            
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic 
    !! Eliashberg equations on the real-axis.  
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_epw,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim
    USE eliashbergcom, ONLY : nsw, Delta, Deltap, gap, estemp
    USE constants_epw, ONLY : kelvin2eV, ci
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp, iter
    REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv 
    CHARACTER (len=256) :: filgap
    !
    CALL start_clock( 'iso_raxis' ) 
    !
    WRITE(stdout,'(5x,a)') 'Solve isotropic Eliashberg equations on real-axis'
    !
    CALL gen_freqgrid_raxis
    !
    DO itemp = 1, nstemp ! loop over temperature
       !
       WRITE(stdout,'(a)') '    '
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
       WRITE(stdout,'(a)') '    '
       iter = 1
       conv = .false.
       DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL integrate_eliashberg_iso_raxis( itemp, iter, conv )
          rdeltain(:) = real(Deltap(:))
          cdeltain(:) = aimag(Deltap(:))
          rdeltaout(:) = real(Delta(:))
          cdeltaout(:) = aimag(Delta(:))
          CALL mix_broyden ( nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
          CALL mix_broyden2( nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
          Deltap(:) = rdeltain(:) + ci * cdeltain(:)
          iter = iter + 1
       ENDDO ! iter
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a,f10.6,a,a,f10.6,a)') &
                    'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K ', &
                    '  gap_edge(', itemp, ') = ', gap(itemp), ' eV ', &
                    '  Re[Delta(1)] = ', real(Delta(1)), ' eV '
       WRITE(stdout,'(a)') '    '
       tcpu = get_clock( 'iso_raxis' )
       WRITE( stdout,'(5x,a,i3,a,f8.1,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
       !
       IF ( conv ) THEN
          WRITE(stdout,'(a)') '    '
          CALL print_clock( 'iso_raxis' )
          WRITE(stdout,'(a)') '    '
       ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'iso_raxis' )
          CALL print_clock( 'iso_raxis' )
          CALL errore('integrate_eliashberg_iso_raxis','converged was not reached',1)
          RETURN
       ENDIF
       !
    ENDDO ! itemp
    filgap = TRIM(prefix) // '.gap'
    OPEN(iufilgap, file=filgap, status='unknown')
    DO itemp = 1, nstemp ! loop over temperature
       WRITE(iufilgap,'(2f12.6)') estemp(itemp)/kelvin2eV, gap(itemp)
    ENDDO
    CLOSE(iufilgap)
    !    
    CALL stop_clock( 'iso_raxis' )
    ! 
    RETURN
    !
    END SUBROUTINE eliashberg_iso_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE integrate_eliashberg_iso_raxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations on the real-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilker, iufilgap
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : nswfc, nqstep, nsiter, muc, conv_thr_raxis, &
                              kerwrite, kerread, nstemp
    USE eliashbergcom, ONLY : nsw, estemp, ws, dws, gap0, gap, fdwp, Kp, Km, & 
                              Delta, Deltap, Znorm
    USE constants_epw, ONLY : kelvin2eV, ci
    USE superconductivity_iso, ONLY : kernel_raxis
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged  
    ! 
    ! Local variables
    INTEGER :: iw, iwp
    REAL(DP) :: dstep, a, b, c, d, absdelta, reldelta, errdelta, temp
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    REAL(DP) :: eps=1.0d-6
    COMPLEX(DP) :: kernelp, kernelm, esqrt
    COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    REAL(DP), EXTERNAL :: wgauss
    LOGICAL :: lgap
    CHARACTER(len=256) :: name1, name2
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsw) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsw) )
    !
    IF ( iter .eq. 1 ) THEN 
       IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
       IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
       IF ( .not. ALLOCATED(Deltap) ) ALLOCATE( Deltap(nsw) )
       IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
       gap(itemp) = 0.d0
       Deltap(:)  = (0.d0, 0.d0)
       Deltap(:)  = gap0 
       IF ( .not. ALLOCATED(fdwp) ) ALLOCATE( fdwp(nsw) )
       IF ( .not. ALLOCATED(Kp) )   ALLOCATE( Kp(nsw,nsw) )
       IF ( .not. ALLOCATED(Km) )   ALLOCATE( Km(nsw,nsw) )
    ENDIF
    Delta(:) = (0.d0, 0.d0)
    Znorm(:) = (0.d0, 0.d0)
    !
    temp = estemp(itemp) / kelvin2eV
    IF ( temp .lt. 10.d0 ) THEN  
       WRITE(name2,'(a,a7,f4.2)') TRIM(prefix),'.ker_00', temp
    ELSEIF ( temp .ge. 10.d0 ) THEN 
       WRITE(name2,'(a,a6,f5.2)') TRIM(prefix),'.ker_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN 
       WRITE(name2,'(a,a5,f6.2)') TRIM(prefix),'.ker_', temp
    ENDIF
    OPEN(iufilker, file=name2, form='unformatted')
    !
    IF ( iter .eq. 1 ) THEN
       IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
       Deltaold(:) = gap0
    ENDIF          
    absdelta = 0.d0
    reldelta = 0.d0
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nsw ! loop over omega_prime
        IF ( iter .eq. 1 ) THEN
          IF ( iw .eq. 1 ) THEN
            IF ( ABS(estemp(itemp)) <  eps ) THEN
               fdwp(iwp) = 0.d0
            ELSE
               fdwp(iwp) = wgauss( -ws(iwp) / estemp(itemp), -99 )
            ENDIF
          ENDIF
          !
          ! read the kernels from file if they were calculated before otherwise calculate them
          IF ( kerread ) THEN 
            READ(iufilker) a, b, c, d
            Kp(iw,iwp) = a + ci*b
            Km(iw,iwp) = c + ci*d
          ENDIF
          IF ( kerwrite ) THEN 
            CALL kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
            Kp(iw,iwp) = kernelp
            Km(iw,iwp) = kernelm
            WRITE(iufilker) real(Kp(iw,iwp)), aimag(Kp(iw,iwp)), &
                            real(Km(iw,iwp)), aimag(Km(iw,iwp))
          ENDIF
        ENDIF
        !
        ! this step is performed at each iter step only for iw=1 since it is independ of w(iw)
        IF ( iw .eq. 1 ) THEN
           esqrt = 1.d0 / sqrt( ws(iwp)**2.d0 - Deltap(iwp)**2.d0 )
           wesqrt(iwp) =  real( ws(iwp) * esqrt )
           desqrt(iwp) =  real( Deltap(iwp) * esqrt )
        ENDIF
        !
        ! end points contribute only half ( trapezoidal integration rule )
        IF ( (iwp .eq. 1) .OR. (iwp .eq. nsw) ) THEN
           dstep = 0.5d0 * dws(iwp) 
        ! boundary points contribute half from left and half from right side
        ELSEIF ( iwp .eq. nswfc ) THEN
           dstep = 0.5d0 * ( dws(iwp) + dws(iwp+1) )
        ELSE
           dstep = dws(iwp)
        ENDIF 
        Znorm(iw) = Znorm(iw) + dstep * wesqrt(iwp) * Km(iw,iwp)
        Delta(iw) = Delta(iw) + dstep * desqrt(iwp) &
                  * ( Kp(iw,iwp) - muc*( 1.d0 - 2.d0*fdwp(iwp) ) )
      ENDDO ! iwp
      Znorm(iw) = 1.d0 - Znorm(iw) / ws(iw)
      Delta(iw) = Delta(iw) / Znorm(iw)
      reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) * dws(iw)
      absdelta = absdelta + abs( Delta(iw) ) * dws(iw)
    ENDDO ! iw 
    CLOSE(iufilker)
    errdelta = reldelta / absdelta
    Deltaold(:) = Delta(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, '   error = ', errdelta, & 
                                  '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1)) 
    !
    IF ( errdelta .lt. conv_thr_raxis) conv = .true.
    IF ( errdelta .lt. conv_thr_raxis .OR. iter .eq. nsiter ) THEN
      IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a8,f4.2)') TRIM(prefix),'.gapr_00', temp
      ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a7,f5.2)') TRIM(prefix),'.gapr_0', temp
      ELSEIF ( temp .ge. 100.d0 ) THEN
        WRITE(name1,'(a,a6,f6.2)') TRIM(prefix),'.gapr_', temp
      ENDIF
      OPEN(iufilgap, file=name1, form='formatted')
      !
      WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
      lgap = .true.
      ! DO iw = 1, nsw
      DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
        IF ( lgap .AND. iw .lt. (nqstep) .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
             ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN
           gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                      / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )
           lgap = .false.
        ENDIF
        WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                    real(Delta(iw)), aimag(Delta(iw))
      ENDDO
      CLOSE(iufilgap)
      IF ( lgap ) & 
         gap(itemp) =  real(Delta(1))
      gap0 = gap(itemp)
    ENDIF
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
       IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
       CALL errore('integrate_eliashberg_iso_raxis','increase nsiter or reduce conv_thr_raxis',-1)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE integrate_eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE dos_quasiparticle( itemp )
    !-----------------------------------------------------------------------
    !!
    !! Computes the quasiparticle density of states in the superconducting state
    !!
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iuqdos
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : lreal, limag, liso, laniso, fsthick
    USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, dws, Delta, ADelta, & 
                              wkfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV, ci
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    INTEGER :: iw, ik, ibnd
    REAL(kind=DP) :: degaussw0, temp, weight
    REAL(kind=DP), ALLOCATABLE :: dos_qp(:)
    COMPLEX(kind=DP) :: omega
    CHARACTER (len=256) :: fildos
    !
    ! SP: This needs to be initialized
    degaussw0 = 0.0_DP
    !
    IF ( lreal ) THEN 
       degaussw0 = 1.d0 * dws(1)
    ELSEIF ( limag ) THEN 
       degaussw0 = 1.d0 * dwsph
    ENDIF
    !
    temp = estemp(itemp) / kelvin2eV
    IF ( temp .lt. 10.d0 ) THEN
       WRITE(fildos,'(a,a8,f4.2)') TRIM(prefix),'.qdos_00', temp
    ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0 ) THEN
       WRITE(fildos,'(a,a7,f5.2)') TRIM(prefix),'.qdos_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN
       WRITE(fildos,'(a,a6,f6.2)') TRIM(prefix),'.qdos_', temp
    ENDIF
    OPEN(iuqdos, file=fildos, form='formatted')
    !
    IF ( .not. ALLOCATED(dos_qp) ) ALLOCATE( dos_qp(nsw) )
    dos_qp(:) = 0.d0          
    !
    IF ( laniso ) THEN
       WRITE(iuqdos,'(5a20)') 'w [eV]', 'N_S/N_F'
       DO iw = 1, nsw 
          omega = ws(iw) + ci*degaussw0
          DO ik = 1, nkfs
             DO ibnd = 1, nbndfs
                IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                   weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik)
                   dos_qp(iw) = dos_qp(iw) + weight & 
                              * real( omega / sqrt( omega*omega - ADelta(ibnd,ik,iw)*ADelta(ibnd,ik,iw) ) ) 
                ENDIF
             ENDDO
          ENDDO
          WRITE(iuqdos,'(2ES20.10)') ws(iw), dos_qp(iw)
       ENDDO
    ELSEIF ( liso ) THEN 
       WRITE(iuqdos,'(5a20)') 'w [eV]', 'N_S/N_F'
       DO iw = 1, nsw
          omega = ws(iw) + ci*degaussw0
          dos_qp(iw) = dos_qp(iw) + real( omega / sqrt( omega*omega - Delta(iw)*Delta(iw) ) ) 
          WRITE(iuqdos,'(2ES20.10)') ws(iw), dos_qp(iw)
       ENDDO
    ENDIF
    CLOSE(iuqdos)
    !
    IF ( ALLOCATED(dos_qp) ) DEALLOCATE(dos_qp)
    !
    RETURN
    !
    END SUBROUTINE dos_quasiparticle
    !
    !-----------------------------------------------------------------------
    SUBROUTINE free_energy( itemp )
    !-----------------------------------------------------------------------
    !!
    !! Computes the free energy difference between the superconducting and normal
    !! states
    !!
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufe
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : liso, laniso, fsthick
    USE eliashbergcom, ONLY : estemp, wsi, nsiw, ADeltai, AZnormi, NAZnormi, &
                              Deltai, Znormi, NZnormi, &
                              wkfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : pi, kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    ! 
    ! Local variables
    INTEGER :: iw, ik, ibnd 
    REAL(DP) :: dFE, omega, temp, weight
    CHARACTER (len=256) :: filfe
    !
    temp = estemp(itemp) / kelvin2eV
    IF ( temp .lt. 10.d0 ) THEN
       WRITE(filfe,'(a,a6,f4.2)') TRIM(prefix),'.fe_00', temp
    ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0 ) THEN
       WRITE(filfe,'(a,a5,f5.2)') TRIM(prefix),'.fe_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN
       WRITE(filfe,'(a,a4,f6.2)') TRIM(prefix),'.fe_', temp
    ENDIF
    OPEN(iufe, file=filfe, form='formatted')
    !
    dFE = 0.d0
    IF ( laniso ) THEN
       DO iw = 1, nsiw(itemp)
          DO ik = 1, nkfs
             DO ibnd = 1, nbndfs
                IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                   weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik)
                   omega = sqrt( wsi(iw)*wsi(iw) + ADeltai(ibnd,ik,iw)*ADeltai(ibnd,ik,iw) )
                   dFE = dFE - weight * ( omega - wsi(iw) ) & 
                       * ( AZnormi(ibnd,ik,iw) - NAZnormi(ibnd,ik,iw) * wsi(iw) / omega )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ELSEIF ( liso ) THEN
       DO iw = 1, nsiw(itemp) 
          omega = sqrt( wsi(iw)*wsi(iw) + Deltai(iw)*Deltai(iw) )
          dFE = dFE - ( omega - wsi(iw) ) &
              * ( Znormi(iw) - NZnormi(iw) * wsi(iw) / omega )
       ENDDO
    ENDIF
    dFE = dFE * pi * estemp(itemp)
    WRITE(iufe,'(2ES20.10)') temp, dFE
    CLOSE(iufe)
    !
    RETURN
    !
    END SUBROUTINE free_energy
    !
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_cont_raxis( itemp, cname )
    !-----------------------------------------------------------------------
    !
    !
    ! This routine writes to files results from the solutions of the Eliashberg
    ! equations on the real-axis 
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso
    USE eliashbergcom, ONLY : nsw, estemp, ws, gap, Agap, &
                              Delta, Znorm, ADelta, AZnorm, &
                              nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    USE io_eliashberg, ONLY : gap_distribution_FS
    !
    IMPLICIT NONE
    !
    INTEGER :: iw, itemp, ik, ibnd
    REAL(DP) :: temp
    LOGICAL :: lgap
    CHARACTER(len=256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    IF ( laniso ) THEN 
       IF ( iverbosity .eq. 2 ) THEN
          IF ( temp .lt. 10.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a9,f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
          ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0  ) THEN
             WRITE(name1,'(a,a1,a4,a8,f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
          ELSEIF ( temp .ge. 100.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a7,f6.2)') TRIM(prefix), '.', cname, '_aniso_', temp
          ENDIF
          OPEN(iufilgap, file=name1, form='formatted')
          WRITE(iufilgap,'(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]',&
                                                            'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       ENDIF
       !
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                lgap = .true.
                ! DO iw = 1, nsw
                DO iw = 1, nsw-1   ! FG: this change is to prevent segfault in ws(iw+1) and ADelta(*,*,iw+1)
                   IF ( lgap .AND. iw .lt. nqstep .AND. real(ADelta(ibnd,ik,iw)) .gt. 0.d0 &
                        .AND. real(ADelta(ibnd,ik,iw+1)) .gt. 0.d0 &
                        .AND. ( ws(iw) - real(ADelta(ibnd,ik,iw)) )*( ws(iw+1) - real(ADelta(ibnd,ik,iw+1)) ) .lt. 0.d0 ) THEN
                      Agap(ibnd,ik,itemp) = (   ( real(ADelta(ibnd,ik,iw))   - ws(iw)   ) * ws(iw+1) &
                                              - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) * ws(iw) ) &
                                          / ( ( real(ADelta(ibnd,ik,iw)) - ws(iw) ) - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) )
                      lgap = .false.
                   ENDIF
                   IF ( iverbosity .eq. 2 ) THEN
                      WRITE(iufilgap,'(6ES20.10)') ws(iw), ekfs(ibnd,ik)-ef0, &
                                     real(AZnorm(ibnd,ik,iw)), aimag(AZnorm(ibnd,ik,iw)), &
                                     real(ADelta(ibnd,ik,iw)), aimag(ADelta(ibnd,ik,iw))
                   ENDIF
                ENDDO ! iw
                IF ( lgap ) & 
                   Agap(ibnd,ik,itemp) = real(ADelta(ibnd,ik,1))
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik
       IF ( iverbosity .eq. 2 ) CLOSE(iufilgap)
       !
       CALL gap_distribution_FS ( itemp, cname )
       !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF ( ( laniso .AND. iverbosity .eq. 2 ) .OR. liso ) THEN
       IF ( temp .lt. 10.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a7,f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
       ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0  ) THEN
          WRITE(name1,'(a,a1,a4,a6,f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
       ELSEIF ( temp .ge. 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a5,f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
       ENDIF
       OPEN(iufilgap, file=name1, form='formatted')
       WRITE(iufilgap,'(5a20)') 'w [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       lgap = .true.
       ! DO iw = 1, nsw
       DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
          IF ( lgap .AND. iw .lt. nqstep .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
               ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN
             gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                        / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )
             lgap = .false.
          ENDIF
          WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                       real(Delta(iw)), aimag(Delta(iw))
       ENDDO ! iw
       CLOSE(iufilgap)
       IF ( lgap ) & 
          gap(itemp) = real(Delta(1))
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE eliashberg_write_cont_raxis
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_memlt_aniso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! Estimate the memory requirements for anisotropic Eliashberg equations 
    !! on imaginary axis
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : max_memlt
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, nqfs, limag_fly, memlt_pool
    USE mp_global, ONLY : inter_pool_comm, my_pool_id
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: itemp, lower_bnd, upper_bnd, imelt
    REAL(DP) :: rmelt
    !
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    my_pool_id = 0
#endif  
    !
    limag_fly = .false.
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of the AKeri kernels that need to stored in each pool
    ! imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
    ! rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
    ! RM - avoid problems when imelt is greater than (2^31)-1 (singed long integer) 
    !
    imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:))
    rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
    imelt = nbndfs**2 * ( 2 * nsiw(itemp) )
    rmelt = dble(imelt) * rmelt 
    rmelt = rmelt + memlt_pool(my_pool_id+1)
    !
    memlt_pool(:) = 0.d0
    memlt_pool(my_pool_id+1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum( memlt_pool, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
       limag_fly = .true.
       !
       ! remove memory required for AKeri 
       ! imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
       ! RM - avoid problems when imelt is greater than (2^31)-1 (singed long integer)
       !
       imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:))
       rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
       imelt = nbndfs**2 * ( 2 * nsiw(itemp) )
       rmelt = dble(imelt) * rmelt
       rmelt = - rmelt + memlt_pool(my_pool_id+1)
       !
       memlt_pool(:) = 0.d0
       memlt_pool(my_pool_id+1) = rmelt
       !
    ! collect contributions from all pools
    CALL mp_sum( memlt_pool, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
  
    ENDIF
    !
    IF ( limag_fly ) THEN
       WRITE(stdout,'(/,5x,a/)') "AKeri is calculated on the fly since its size exceedes max_memlt"
    ELSE
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE eliashberg_memlt_aniso_iaxis
    !  
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_memlt_aniso_acon
    !-----------------------------------------------------------------------
    !!  
    !! Estimate the memory requirements for the anisotropic Eliashberg funtion
    !! used for analytic continuation from imaginary to real axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, max_memlt
    USE eliashbergcom, ONLY : nkfs, nbndfs, nqfs, lacon_fly, memlt_pool
    USE mp_global, ONLY : inter_pool_comm, my_pool_id
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: lower_bnd, upper_bnd, imelt
    REAL(DP) :: rmelt
    !
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    my_pool_id = 0
#endif  
    !
    lacon_fly = .false.
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of a2fij that need to stored in each pool
    imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
    rmelt = dble(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
    rmelt = rmelt + memlt_pool(my_pool_id+1)
    !
    memlt_pool(:) = 0.d0
    memlt_pool(my_pool_id+1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum( memlt_pool, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( maxval(memlt_pool(:)) .gt. max_memlt ) THEN
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of required memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
       lacon_fly = .true.
       !
       ! remove memory required for a2fij
       imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
       CALL mem_size_eliashberg( -imelt )
    ENDIF
    !
    IF ( lacon_fly ) THEN
       WRITE(stdout,'(/,5x,a/)') "a2fij is calculated on the fly since its size exceedes max_memlt"
    ELSE
       WRITE(stdout,'(/,5x,a,a,f9.4,a)') "Size of allocated memory per pool :", &
            " ~= ", maxval(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE eliashberg_memlt_aniso_acon
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_raxis
    !-----------------------------------------------------------------------
    !!
    !! Automatic generation of the frequency-grid for real-axis calculations.
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nswfc, nswc, pwc, wsfc, wscut, lunif
    USE eliashbergcom, ONLY : nsw, ws, dws
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iw
    !
    ! define a grid ws in 2 step sizes
    ! 1. a fine grid of nswfc points between (0,wsfc)
    ! 2. a rough grid of nswc points between (wsfc,wscut).
    ! above wsfc the gap function varies slowly
    !
    ! nswfc = nr. of grid points between (0,wsfc)
    ! nswc  = nr. of grid points between (wsfc,wscut)
    !
    WRITE(stdout,'(a)') '    '
    WRITE(stdout,'(5x,a,i6,a)') 'Total number of nsw = ', nsw, ' grid-points are divided in:'
    WRITE(stdout,'(5x,a,i6,a,f12.6,a,f12.6)') 'nswfc = ', nswfc, '  from ', 0.0, ' to ', wsfc  
    WRITE(stdout,'(5x,a,i6,a,f12.6,a,f12.6)') 'nswc  = ', nswc,  '  from ', wsfc, ' to ', wscut
    WRITE(stdout,'(a)') '    '
    !
    IF ( .not. ALLOCATED(ws) )  ALLOCATE( ws(nsw) )
    IF ( .not. ALLOCATED(dws) ) ALLOCATE( dws(nsw) )
    ws(:) = 0.d0
    dws(:) = 0.d0
    !
    DO iw = 1, nswfc
       dws(iw) = wsfc / dble(nswfc)
       ws(iw) = dble(iw) * dws(iw)
    ENDDO
    DO iw = nswfc + 1, nsw
       dws(iw) = ( wscut - wsfc ) / dble(nswc)
       IF ( lunif ) THEN 
          ws(iw) = wsfc + dble(iw) * dws(iw)
       ELSE 
          ! RM this needs to be checked
          ws(iw) = wsfc + dble( iw/nswc )**pwc * (wscut - wsfc)
       ENDIF
    ENDDO
    !
    IF ( .not. lunif ) THEN 
       DO iw = nswfc+1, nsw-1
          dws(iw) = ws(iw+1) - ws(iw)
       ENDDO
       dws(nsw) = dws(nsw-1)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE gen_freqgrid_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gen_freqgrid_iaxis( itemp )
    !-----------------------------------------------------------------------
    !
    ! Automatic generation of the frequency-grid for imaginary-axis calculations.
    !
    !
    ! input
    !
    ! itemp  - temperature point
    !
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : nqstep, lpade, lacon 
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, wsph, dwsph, estemp, wsphmax
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iw, itemp, n, imelt
    !
    ! frequency-grid for imaginary-axis
    ! nsiw(itemp) = nr. of grid points between (0,wscut) 
    !
    ! memory allocated for wsi and ws
    imelt = nsiw(itemp) + nsw 
    CALL mem_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(wsi) )  ALLOCATE( wsi(nsiw(itemp)) )
    wsi(:) = 0.d0
    DO iw = 1, nsiw(itemp)
       n = iw - 1
       wsi(iw) = dble(2*n+1) * pi * estemp(itemp) 
       !WRITE(*,*) iw, wsi(iw)
    ENDDO
    !
    ! frequency-grid for real-axis ( Pade approximants and analytic continuation)
    !
    IF ( lpade .OR. lacon ) THEN
       IF ( .not. ALLOCATED(ws) )  ALLOCATE( ws(nsw) )
       ws(:) = 0.d0
       DO iw = 1, nsw
          IF ( iw .le. nqstep ) THEN 
             ws(iw) = wsph(iw)
          ELSE
             ws(iw) = wsphmax + dble(iw-nqstep)*dwsph
          ENDIF
          !WRITE(*,*) iw, ws(iw), wsph(iw)
       ENDDO
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE gen_freqgrid_iaxis
    ! 
    ! General note on PADE:
    ! Lebegue, Arnaud, Alouani, and Blochel [PRB 67, 155208 (2003)]
    ! state that when they use Pade of order N = 12 (resulting in
    ! numerator or order (N-2)/2 = 5 and denominator N/2 = 6),
    ! they obtain extremely stable fits and the quasiparticle energies
    ! are essentially identical to those obtained using the contour
    ! deformation method.
    !
    ! using this sub:
    !
    ! integer :: N
    ! complex(DP) :: z(N), u(N), a(N), w, padapp
    !
    ! call pade_coeff ( N, z, u, a)
    ! call pade_eval ( N, z, a, w, padapp)
    !
    !-----------------------------------------------------------
    subroutine pade_coeff ( N, z, u, a)
    !-----------------------------------------------------------
    ! N-point Pade' approximant - find the Pade' coefficients
    !
    ! This subroutine uses the recursive algorithm described in
    ! HJ Vidberg and JW Serene, "Solving the Eliashberg equations
    ! by means of N-point Pade' approximants", J Low Temp Phys
    ! 29, 179 (1977). The notation adopted here is the same as
    ! in the above manuscript.
    !
    ! input
    !
    ! N      - order of the Pade' approximant
    ! z(1:N) - points at which the original function is known
    ! u(1:N) - values of the function at the z points
    !
    ! output
    !
    ! a(1:N) - coefficients of the continued fraction
    !-----------------------------------------------------------
    !
    USE kinds, ONLY : DP
#if ! defined(__GFORTRAN__) || (__GNUC__ > 4 )
    USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
    !
    implicit none
    INTEGER :: N
    INTEGER :: i 
    INTEGER :: p
    REAL(kind=DP) :: ar
    !! Real part
    REAL(kind=DP) :: ai
    !! Complex part
    COMPLEX(kind=DP) :: z(N), u(N)
    COMPLEX(kind=DP) :: g(N,N), a(N)
    ! g(p,i) = g_p (z_i) in the notation of Vidberg and Serene
    COMPLEX(kind=DP) :: tmp1, tmp2
    !
    do p = 1, N
      if (p.eq.1) then
        do i = 1, N
           g (p,i) = u(i)
        enddo
      else
        do i = p, N
        !  g (p,i) = ( g(p-1,p-1) - g(p-1,i) ) / &
         !           ( ( z(i) - z(p-1) ) * g (p-1,i) )
           !
           ! this seems necessary to avoid nasty NaN when
           ! still don't quite understand why the procedure
           ! becomes unstable - certainly it happens only
           ! when u(:) is very small
           !
  !if(abs(g(p-1,i)) .eq. 0) then
  !       write(6,'(4x, "fitting parameter too small. g(p-1,i)= ",2f9.5)')g(p-1,i)
  !       stop
  !end if
  !
           tmp1 = g(p-1,p-1)/g(p-1,i)
           tmp2 = g(p-1,i)/g(p-1,i)
           g (p,i) = ( tmp1 - tmp2 ) / ( z(i) - z(p-1) )
           !
        enddo
      endif
      a(p) = g (p,p)
      !
      ! check whether a(p) is not NaN
      !
      ar = real(a(p))
      ai = aimag(a(p))
#if ! defined(__GFORTRAN__) || (__GNUC__ > 4 )
      IF ( IEEE_IS_NAN(ar) .or. IEEE_IS_NAN(ai) ) THEN
  !     write(6,*) (z(i),i=1,N)
  !     write(6,*) (u(i),i=1,N)
  !     write(6,*) (a(i),i=1,N)
        write(6,*) 'one or more coefficients are NaN'
  !     call errore('pade_coeff','one or more coefficients are NaN',1)
      ENDIF
#endif
      !
    enddo
    !
    end subroutine pade_coeff
    !
    !-----------------------------------------------------------
    subroutine pade_eval ( N, z, a, w, padapp)
    !-----------------------------------------------------------
    ! N-point Pade' approximant - evaluate the Pade' approximant
    !
    ! This subroutine uses the recursive algorithm described in
    ! HJ Vidberg and JW Serene, "Solving the Eliashberg equations
    ! by means of N-point Pade' approximants", J Low Temp Phys
    ! 29, 179 (1977). The notation adopted here is the same as
    ! in the above manuscript.
    !
    ! input
    !
    ! N      - order of the Pade' approximant
    ! z(1:N) - points at which the original function is known
    ! a(1:N) - coefficients of the continued fraction
    ! w      - point at which we need the approximant
    !
    ! output
    !
    ! padapp - value of the approximant at the point w
    !-----------------------------------------------------------
    !
    USE kinds,          ONLY : DP
    implicit none
    integer :: N
    complex(DP) :: a(N), z(N), acap(0:N), bcap(0:N)
    complex(DP) :: w, padapp
    integer :: i
    !
    acap(0) = 0.d0
    acap(1) = a(1)
    bcap(0) = 1.d0
    bcap(1) = 1.d0
    !
    do i = 2, N
      acap(i) = acap(i-1) + (w-z(i-1)) * a(i) * acap(i-2)
      bcap(i) = bcap(i-1) + (w-z(i-1)) * a(i) * bcap(i-2)
    enddo
    padapp = acap(N)/bcap(N)
    !
    end subroutine pade_eval
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gamma_acont( omega, omegap, temp, rgammap, rgammam )
    !-----------------------------------------------------------------------
    !!
    !! computes gammam(w,wp)  (notes RM)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : eps6
    ! 
    IMPLICIT NONE
    !
    REAL(kind=DP), INTENT (in) :: omega
    !! frequency w at point iw on the real-axis
    REAL(kind=DP), INTENT (in) :: omegap
    !! frequency w' at point iwp on the real-axis
    REAL(kind=DP), INTENT (in) :: temp
    !! temperature in eV
    REAL(kind=DP), INTENT (out) :: rgammap
    !! -bose_einstein( w' ) - fermi_dirac(  w + w' )
    REAL(kind=DP), INTENT (out) :: rgammam
    !! bose_einstein( w' ) + fermi_dirac( -w + w' )
    ! 
    rgammap = 0.d0
    rgammam = 0.d0
    IF ( ABS(temp) < eps6 ) THEN
       rgammap = 0.d0
       rgammam = 1.d0
    ELSEIF ( omegap .gt. 0.d0 ) THEN 
       rgammap = 0.5d0 * (   tanh( 0.5d0 * ( omega + omegap ) / temp ) &
                           - 1.d0 / tanh( 0.5d0 * omegap / temp ) )
       rgammam = 0.5d0 * (   tanh( 0.5d0 * ( omega - omegap ) / temp ) &
                           + 1.d0 / tanh( 0.5d0 * omegap / temp ) )
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE gamma_acont
    !                         
    !-----------------------------------------------------------------------
    !                                                                            
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_eliashberg
    !!
    USE epwcom,        ONLY : liso, laniso, limag
    USE eliashbergcom, ONLY : wsph, estemp, gap, wsph, agap
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(wsph) )  DEALLOCATE(wsph)
    IF( ALLOCATED(estemp)) DEALLOCATE(estemp)
    !
    IF ( liso ) THEN 
       IF ( limag ) THEN
          CALL deallocate_eliashberg_iso_iaxis
       ENDIF
       CALL deallocate_eliashberg_iso_raxis
    ENDIF
    !
    IF ( laniso ) THEN
       IF ( limag ) THEN
          CALL deallocate_eliashberg_aniso_iaxis
       ENDIF
       CALL deallocate_eliashberg_aniso_raxis
       IF( ALLOCATED(gap))  DEALLOCATE(gap)
       IF( ALLOCATED(Agap)) DEALLOCATE(Agap)
    ENDIF
    !
    CALL deallocate_elphon
    !
    RETURN
    !
    END SUBROUTINE deallocate_eliashberg
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_iso_iaxis
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_eliashberg_iso_iaxis
    !!
    !----------------------------------------------------------------------
    !
    USE eliashbergcom
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(wsi) )     DEALLOCATE(wsi)
    IF( ALLOCATED(Deltai) )  DEALLOCATE(Deltai)
    IF( ALLOCATED(Deltaip) ) DEALLOCATE(Deltaip)
    IF( ALLOCATED(Znormi) )  DEALLOCATE(Znormi)
    IF( ALLOCATED(NZnormi) ) DEALLOCATE(NZnormi)
    IF( ALLOCATED(Keri) )    DEALLOCATE(Keri)
    !
    RETURN
    !
    END SUBROUTINE deallocate_eliashberg_iso_iaxis
    !                         
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_iso_raxis
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_eliashberg_iso_raxis
    !!
    USE epwcom, ONLY : lreal, limag, lacon
    USE eliashbergcom
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(ws) )    DEALLOCATE(ws)
    !
    IF( ALLOCATED(Delta))  DEALLOCATE(Delta)
    IF( ALLOCATED(Deltap)) DEALLOCATE(Deltap)
    IF( ALLOCATED(Znorm))  DEALLOCATE(Znorm)
    IF( ALLOCATED(Znormp)) DEALLOCATE(Znormp)
    IF( ALLOCATED(gap))    DEALLOCATE(gap)
    !
    IF ( lreal ) THEN
       IF( ALLOCATED(dws) )   DEALLOCATE(dws)
       IF( ALLOCATED(fdwp) )  DEALLOCATE(fdwp)
       IF( ALLOCATED(bewph) ) DEALLOCATE(bewph)
       IF( ALLOCATED(Kp))     DEALLOCATE(Kp)    
       IF( ALLOCATED(Km))     DEALLOCATE(Km) 
    ENDIF
    !
    IF ( limag .AND. lacon ) THEN
       IF( ALLOCATED(Gp))     DEALLOCATE(Gp)
       IF( ALLOCATED(Gm))     DEALLOCATE(Gm)
       IF( ALLOCATED(Dsumi) ) DEALLOCATE(Dsumi)
       IF( ALLOCATED(Zsumi) ) DEALLOCATE(Zsumi)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE deallocate_eliashberg_iso_raxis
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_aniso_iaxis
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_eliashberg_aniso_iaxis
    !!
    USE eliashbergcom
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(wsi) )      DEALLOCATE(wsi)
    !
    IF( ALLOCATED(Deltai) )   DEALLOCATE(Deltai)
    IF( ALLOCATED(Znormi) )   DEALLOCATE(Znormi)
    ! 
    IF( ALLOCATED(ADeltai) )  DEALLOCATE(ADeltai)
    IF( ALLOCATED(ADeltaip) ) DEALLOCATE(ADeltaip)
    IF( ALLOCATED(AZnormi) )  DEALLOCATE(AZnormi)
    IF( ALLOCATED(NAZnormi) ) DEALLOCATE(NAZnormi)
    !
    RETURN
    !
    END SUBROUTINE deallocate_eliashberg_aniso_iaxis
    !                                
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_eliashberg_aniso_raxis
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by allocate_eliashberg_aniso_raxis
    !!
    USE eliashbergcom
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(ws))       DEALLOCATE(ws) 
    !
    IF( ALLOCATED(Delta))    DEALLOCATE(Delta)
    IF( ALLOCATED(Znorm))    DEALLOCATE(Znorm)
    !
    IF( ALLOCATED(ADelta) )  DEALLOCATE(ADelta)
    IF( ALLOCATED(ADeltap) ) DEALLOCATE(ADeltap)
    IF( ALLOCATED(AZnorm) )  DEALLOCATE(AZnorm)
    IF( ALLOCATED(AZnormp) ) DEALLOCATE(AZnormp)
    !
    RETURN
    !
    END SUBROUTINE deallocate_eliashberg_aniso_raxis
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_elphon
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by electron-phonon
    !!
    USE elph2,         ONLY : wf, wqf
    USE eliashbergcom, ONLY : ekfs, xkfs, wkfs, xkff, g2, a2f_iso, w0g, ixkff, ixkqf, ixqfs, nqfs
    !
    IMPLICIT NONE
    !
    IF( ALLOCATED(wf) )      DEALLOCATE(wf)
    IF( ALLOCATED(wqf) )     DEALLOCATE(wqf)
    IF( ALLOCATED(ekfs) )    DEALLOCATE(ekfs)
    IF( ALLOCATED(xkfs) )    DEALLOCATE(xkfs)
    IF( ALLOCATED(wkfs) )    DEALLOCATE(wkfs)
    IF( ALLOCATED(xkff) )    DEALLOCATE(xkff)
    IF( ALLOCATED(g2) )      DEALLOCATE(g2)
    IF( ALLOCATED(a2f_iso) ) DEALLOCATE(a2f_iso)
    IF( ALLOCATED(w0g) )     DEALLOCATE(w0g)
    IF( ALLOCATED(ixkff) )   DEALLOCATE(ixkff)
    IF( ALLOCATED(ixkqf) )   DEALLOCATE(ixkqf)
    IF( ALLOCATED(ixqfs) )   DEALLOCATE(ixqfs)
    IF( ALLOCATED(nqfs) )    DEALLOCATE(nqfs)
    !
    RETURN
    !
    END SUBROUTINE deallocate_elphon
    !                                                                            
    !----------------------------------------------------------------------
    ! 
  END MODULE superconductivity
